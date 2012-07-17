#!/usr/bin/env python
"""
General script to convert LOFAR single station correlation ACC files to Measurement Sets
"""

#based on mkant.py by Stephen Bourke, JIVE 
#and aipy by Aaron Parsons

#todo:
# HBA specific setup

import pyrap.tables as pt
from pyrap.tables import tablecreatescalarcoldesc as cldsc
from pyrap.tables import tablecreatearraycoldesc as clarrdsc

import numpy
import ephem

import sys,os
import struct
import time

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-o', '--out', dest='ms_out', default=None,
        help='Output MS filename prefix, if mode is \'single\' default: <array type>_, if mode is \'merge\' default: <array_type>_m<subband0>')
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help='AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help='AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help='Station RCU Mode, usually 3,5,6,7, default: 3(LBA High)')
    o.add_option('-n','--name', dest='station_name', default='LOFAR',
        help='Station name, default: LOFAR')
    o.add_option('-i','--int', dest='int_time', default=1.0, type='float',
        help='Integration length per subband, in seconds, default: 1.0')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help='Select which subband to select out for the MS. Can use a list 200,400,411 or a range 200_210, default:0')
    o.add_option('-m', '--mode', dest='mode', default='single',
        help='Output data format: \'single\' outputs a MS for each subband, \'merge\' outputs a MS with the subbands merged to a single integration time, default: single')
    o.add_option('-t', '--utc', dest='utcOffset', default=0., type='int',
        help='Offset in local clock from UTC where observation was taken as the file timestamp is in local units, e.g. +1 during BST. default:0')
    opts, args = o.parse_args(sys.argv[1:])
acc_files=args

rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000., 'sfreq': 0.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000., 'sfreq': 0.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000., 'sfreq': 0.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000., 'sfreq': 0.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000., 'sfreq': 0.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000., 'sfreq': 100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':80000000., 'sfreq': 1560000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000., 'sfreq': 200000000.}] #7

class tStep:
    def __init__(self):
        self.curtime=time.time()
        self.lasttime=self.curtime
        self.dt=self.curtime-self.lasttime
    def update(self):
        self.lasttime=self.curtime
        self.curtime=time.time()
        self.dt=self.curtime-self.lasttime
        print self.dt

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    init = [map(int, option.split('_')) for option in arg]
    rv = []
    for i in init:
        if len(i) == 1:
            rv.append(i[0])
        elif len(i) == 2:
            rv.extend(range(i[0],i[1]+1))
    return rv

def aeff(freq):
    """Effective collecting area of a LOFAR dipole based on the frequency, in Hz"""
    wl=ephem.c/freq
    d=5.5   #for the 'sparse' international stations this is the average distance between LBA dipoles in meters, core=3.33, dutch=7.2
    if freq<100000000.0:
        return min((ephem.pi*(d**2)/4),(wl**2)/3)
    else: return min((1.5625,(wl**2)/3))

def phsCenterSrc(obs, t):
    """return an ephem FixedBody source based on the time offset from the obs"""
    src=ephem.FixedBody()
    t0=obs.date
    obs.date=t
    src._ra=obs.sidereal_time()
    src._dec=obs.lat
    obs.date=t0
    return src
def eq2top_m(ha, dec):
    """Return the 3x3 matrix converting equatorial coordinates to topocentric
    at the given hour angle (ha) and declination (dec)."""
    sin_H, cos_H = numpy.sin(ha), numpy.cos(ha)
    sin_d, cos_d = numpy.sin(dec), numpy.cos(dec)
    zero = numpy.zeros_like(ha)
    map =  numpy.array([[    sin_H    ,       cos_H  ,       zero  ],
                        [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                        [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
    if len(map.shape) == 3: map = map.transpose([2, 0, 1])
    return map
def get_baseline(i, j, src, obs):
    """Return the baseline corresponding to i,j""" 
    bl = j - i
    try:
        if src.alt < 0:
            raise PointingError('Phase center below horizon')
        m=src.map
    except(AttributeError):
        ra,dec = src._ra,src._dec
        #opposite HA since the we want to move the source at zenith away to phase to the original zenith source
        m = eq2top_m(ra-obs.sidereal_time(), dec)
        #normal HA
        #m = eq2top_m(obs.sidereal_time() - ra, dec)
    return numpy.dot(m, bl).transpose()
def gen_uvw(i, j, src, obs, f):
    """Compute uvw coordinates of baseline relative to provided FixedBody"""
    x,y,z = get_baseline(i,j,src,obs)
    afreqs = numpy.reshape(f, (1,f.size))
    afreqs = afreqs/ephem.c #1/wavelength
    if len(x.shape) == 0: return numpy.array([x*afreqs, y*afreqs, z*afreqs])
    x.shape += (1,); y.shape += (1,); z.shape += (1,)
    return numpy.array([numpy.dot(x,afreqs), numpy.dot(y,afreqs), numpy.dot(z,afreqs)])
def gen_phs(i, j, src, obs, f):
    """Return phasing that is multiplied to data to point to src."""
    u,v,w = gen_uvw(i,j,src,obs,f)
    #dw for refraction, ignored
    dw = 0
    #assume no sources are resolved
    res = 1
    #assume no frequency depenedent phase offset
    o = 0
    phs = res * numpy.exp(-1j*2*numpy.pi*(w + dw + o))
    return phs.squeeze()
def phs2src(data, i, j, src, obs, f):
    """Apply phasing to zenith-phased data to point to src."""
    return data * gen_phs(i, j, src, obs, f)

sbs=convert_arg_range(opts.subband)

rcumode=opts.rcumode
array_type=rcuInfo[rcumode]['array_type']
if opts.ms_out==None:
    if opts.mode=='merge':
        ms_out_pre=array_type.lower()+'_m'
    else:
        ms_out_pre=array_type.lower()+'_'
else: ms_out_pre=opts.ms_out

timeTrack=tStep()

print 'Generating Antenna Tables...'
#Parse AntennaArray file
fh=open(opts.ant_array)
readLatLon=False
for line in fh:
    if line.lower().startswith(array_type.lower()):
        readLatLon=True
        continue
    if readLatLon:
        lon=line.split(' ')[2]
        lat=line.split(' ')[3]
        elev=line.split(' ')[4]
        readLatLon=False
        continue
fh.close()
#Parse AntennaField file
fh=open(opts.ant_field)
readXYZ=False
readArrayHDR=False
readArray=False
linecnt=0
nelem=0
ants=[]
for line in fh:
    if line.lower().startswith(array_type.lower()):
        readXYZ=True
        continue
    if readXYZ:
        arr_x=line.split(' ')[2]
        arr_y=line.split(' ')[3]
        arr_z=line.split(' ')[4]
        readXYZ=False
        readArrayHDR=True
        continue
    if readArrayHDR:
        nelem=int(line.split(' ')[0])
        readArrayHDR=False
        readArray=True
        continue
    if readArray and linecnt < nelem:
        cl=' '.join(line.split())
        ants.append(map(float,cl.split(' ')))
        linecnt+=1
fh.close()

ant_list=[]
mount='ALT-AZ'
name_pre=array_type
subbands=512
bw=rcuInfo[rcumode]['bw']
df=(bw/subbands) #bandwidth / 512 subbands
start_freq=rcuInfo[rcumode]['sfreq']
avg_freq=numpy.average(numpy.array(sbs))*df
ant_size=aeff(avg_freq+start_freq)
nants=len(ants)
for a,ant in enumerate(ants):
    ant_parms = {'POSITION': [ant[0]+float(arr_x),ant[1]+float(arr_y),ant[2]+float(arr_z)], 'MOUNT': mount, 'STATION': opts.station_name, 'NAME': name_pre+str(a)}
    ant_parms['DISH_DIAMETER'] = ant_size
    ant_list.append(ant_parms)

# Define columns
offset_desc = clarrdsc('OFFSET', value=float(), ndim=1, shape=[3])
position_desc = clarrdsc('POSITION', value=float(), ndim=1, shape=[3])
type_desc = cldsc('TYPE', value=str())
dish_desc = cldsc('DISH_DIAMETER', value=float())
flag_desc = cldsc('FLAG_ROW', value=bool())
mount_desc = cldsc('MOUNT', value=str())
name_desc = cldsc('NAME', value=str())
station_desc = cldsc('STATION', value=str())

desc = pt.tablecreatedesc([offset_desc, position_desc, type_desc,
        dish_desc, flag_desc, mount_desc, name_desc, station_desc])

# Create and populate our table
ant_table='acc2msAntTable.ms'
table = pt.table(ant_table, desc, nrow=len(ant_list), readonly=False)
table.putcolkeywords('OFFSET', {'QuantumUnits': ['m', 'm', 'm'],
    'MEASINFO': {'Ref': 'ITRF', 'type': 'position'}})
table.putcolkeywords('POSITION', {'QuantumUnits': ['m', 'm', 'm'],
    'MEASINFO': {'Ref': 'ITRF', 'type': 'position'}})
table.putcolkeywords('DISH_DIAMETER', {'QuantumUnits': ['m']})
table[:] = {'DISH_DIAMETER': 25.0, 'OFFSET': [0.0,0.0,0.0], 'TYPE': 'GROUND-BASED'}
for i in range(len(ant_list)):
    table[i] = ant_list[i]
table.close()
print 'Done\n'

#Select subband from each ACC file
print 'Selecting Out Subband(s) %s from ACC Files...'%(str(sbs))
nantpol=nants*2
fullAcc=numpy.zeros((nantpol,nantpol,len(acc_files)*len(sbs)),dtype=complex)
for fid,fn in enumerate(acc_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(acc_files))
    fnArr=fn.split('/')[-1]
    #fnArr=fn.split('_')
    fnArr=fnArr.split('_')
    accSize=fnArr[3].split('x')
    nsb=int(accSize[0])
    nantpol=int(accSize[1])
    accComp=numpy.fromfile(fn,dtype='complex')
    accComp=accComp.reshape(nsb,nantpol,nantpol)
    for sid,sb in enumerate(sbs):
        accSB=accComp[sb]
        fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(int(fnArr[1][0:2])-opts.utcOffset),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
        fTime=time.mktime(fTimeStruct)
        fTime-=(subbands-sb)*opts.int_time #fTimeStruct is the time for the last subband (511)
        if fid==0 and sid==0:
            fullAcc[:,:,0]=accSB
            initTimeStruct=time.localtime(fTime)
            cfgStartTime=[time.strftime('%Y/%m/%d/%H:%M:%S',initTimeStruct)]
            ephemTime=[time.strftime('%Y/%m/%d %H:%M:%S',initTimeStruct)]
            fullTimes=[fTime]
        else:
            fullAcc[:,:,fid*len(sbs)+sid]=accSB
            initTimeStruct=time.localtime(fTime)
            cfgStartTime.append(time.strftime('%Y/%m/%d/%H:%M:%S',initTimeStruct))
            ephemTime.append(time.strftime('%Y/%m/%d %H:%M:%S',initTimeStruct))
            fullTimes.append(fTime)
print 'Done\n'
fullTimes=numpy.array(fullTimes)
fullTimes=fullTimes.reshape(len(acc_files),fullTimes.shape[0]/len(acc_files))
fullAcc=fullAcc.reshape(fullAcc.shape[0],fullAcc.shape[1],len(acc_files),fullAcc.shape[2]/len(acc_files))
if len(acc_files)==1: cfgStepTime=opts.int_time     #set the correct integration time if only porcessing a single ACC file
else: cfgStepTime=fullTimes[1,0]-fullTimes[0,0]     #Use the time offset between the first two files as the time interval for makems, fix integration time at the end

#phase center
obs=ephem.Observer()
obs.long=lon
obs.lat=lat
obs.elevation=float(elev)
obs.epoch=2000.0
if opts.mode.startswith('single'):
    for sid,sb in enumerate(sbs):
        print 'Subband %i'%(sb)
        obs.date=ephemTime[sid]
        #strSrcDec=str(obs.lat)
        #strSrcRA=str(obs.sidereal_time())
        strSrcDec='%frad'%obs.lat
        strSrcRA='%frad'%obs.sidereal_time()
        print '\tPhase Center:'
        print '\tRA: ', strSrcRA
        print '\tDec: ', strSrcDec
        src=ephem.FixedBody()
        src._ra=obs.sidereal_time()
        src._dec=obs.lat

        ##Cyngus A
        #src=ephem.FixedBody()
        #src._ra='19:59:28.30'
        #src._dec='40:44:02.0'
        #print '\tPhase Center:'
        #print '\tRA: %frad'%src._ra
        #print '\tDec: %frad'%src._dec
        #strSrcDec='%frad'%src._dec
        #strSrcRA='%frad'%src._ra
        
        #generate makems config file
        print '\tGenerating configuration file for makems...'
        cfgStr='NParts=1\n'
        cfgStr+='NBands=1\n'
        cfgStr+='NFrequencies=1\n'
        #Set the start frequency based on the subband number and array type
        cfgStartFreq= sb * df + start_freq
        if cfgStartFreq==0.: cfgStartFreq=df
        cfgStr+='StartFreq=%f\n'%cfgStartFreq
        cfgStr+='StepFreq=%f\n'%df
        cfgStr+='StartTime=%s\n'%cfgStartTime[sid]
        cfgStr+='StepTime=%s\n'%cfgStepTime
        cfgStr+='NTimes=%i\n'%len(acc_files)
        cfgStr+='RightAscension=%s\n'%strSrcRA
        cfgStr+='Declination=%s\n'%strSrcDec
        cfgStr+='WriteAutoCorr=T\n'
        cfgStr+='WriteImagerColumns=T\n'
        ms_out='%s%i.ms'%(ms_out_pre,sb)
        cfgStr+='MSName=%s\n'%ms_out
        cfgStr+='AntennaTableName=%s\n'%ant_table
        cfgFN='tempMakems.cfg'
        cfgFH=open(cfgFN,'w')
        cfgFH.write(cfgStr)
        cfgFH.close()
        print '\tDone\n'
        
        #run makems
        print 'Running makems...'
        os.system('which makems')
        os.system('makems %s > makems.log 2>&1'%cfgFN)
        while not os.path.isdir('%s_p0'%ms_out): time.sleep(1)
        pwd=os.getcwd()
        if os.path.exists(pwd+'/'+ms_out): os.system('rm -r %s/%s'%(pwd,ms_out))
        os.system('mv %s_p0 %s'%(ms_out,ms_out))
        print 'Done\n'
        
        #phase to source
        nint=len(acc_files)
        print 'Applying Phase corrections...'
        tSB=ephemTime[sid::len(sbs)]
        for tid,t in enumerate(tSB):
            obs.date=ephemTime[tid]
            src.compute(obs)
            for i in range(nants):
                for j in range(nants):
                    if i==j: continue
                    xyz_i, xyz_j = numpy.array(ant_list[i]['POSITION']), numpy.array(ant_list[j]['POSITION'])
                    fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid]=phs2src(fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid], xyz_i, xyz_j, src, obs, numpy.array(cfgStartFreq))
        print 'done\n'
        
        ##apply phase correction for multiple files
        #nint=len(acc_files)
        #if nint>1:
        #    print 'Applying Phase corrections...'
        #    tSB=ephemTime[sid::len(sbs)]
        #    for tid,t in enumerate(tSB):
        #        if tid==0: continue
        #        src=phsCenterSrc(obs,t)
        #        src.compute(obs)
        #        for i in range(nants):
        #            for j in range(nants):
        #                if i==j: continue
        #                xyz_i, xyz_j = numpy.array(ant_list[i]['POSITION']), numpy.array(ant_list[j]['POSITION'])
        #                fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid]=phs2src(fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid], xyz_i, xyz_j, src, obs, numpy.array(cfgStartFreq))
        #    print 'done\n'

        #fill ms with acc data
        print 'Filling DATA column with Subband data...'
        ms=pt.table(ms_out,readonly=False)
        dc=ms.col('DATA')
        ant1c=ms.col('ANTENNA1')
        ant2c=ms.col('ANTENNA2')
        dc_data=dc.getcol()
        ant1c_data=ant1c.getcol()
        ant2c_data=ant2c.getcol()
        ncorr=len(ant1c_data)/len(acc_files)
        for j in range(nint):
            for i in range(ncorr):
                #ant1=ant1c_data[i]
                #ant2=ant2c_data[i]
                #use top half of correlation matrix
                ant1=ant2c_data[j*ncorr+i]
                ant2=ant1c_data[j*ncorr+i]
                #CORR_TYPE: [9,10,11,12] is [XX,XY,YX,YY]
                dc_data[j*ncorr+i,0]=[  fullAcc[2*ant1,2*ant2,j,sid],
                                        fullAcc[2*ant1,2*ant2+1,j,sid],
                                        fullAcc[2*ant1+1,2*ant2,j,sid],
                                        fullAcc[2*ant1+1,2*ant2+1,j,sid]]
        #correct EXPOSURE and INTERVAL for multiple integration sets
        if nint>1:
            exposure_col=ms.col('EXPOSURE')
            interval_col=ms.col('INTERVAL')
            exp_data=exposure_col.getcol()
            int_data=interval_col.getcol()
            exp_data=numpy.ones_like(exp_data)*opts.int_time
            int_data=numpy.ones_like(int_data)*opts.int_time
            ms.putcol('EXPOSURE', exp_data)
            ms.putcol('INTERVAL', int_data)
        ms.putcol('DATA', dc_data)
        ms.close()
        print 'Done\n'

elif opts.mode.startswith('merge'):
    merge_sid=0
    print 'Using Subband %i as merge reference'%(sbs[merge_sid])
    obs.date=ephemTime[merge_sid]
    strSrcDec=str(obs.lat)
    strSrcRA=str(obs.sidereal_time())
    print '\tPhase Center:'
    print '\tRA: ', strSrcRA
    print '\tDec: ', strSrcDec

    #generate makems config file
    print '\tGenerating configuration file for makems...'
    cfgStr='NParts=1\n'
    cfgStr+='NBands=1\n'
    cfgStr+='NFrequencies=%i\n'%(len(sbs))
    #Set the start frequency based on the subband number and array type
    df=(bw/subbands) #100 MHz / 512 subbands
    cfgStartFreq= sbs[merge_sid] * df + start_freq
    cfgStr+='StartFreq=%f\n'%cfgStartFreq
    cfgStr+='StepFreq=%f\n'%df
    cfgStr+='StartTime=%s\n'%cfgStartTime[merge_sid]
    cfgStr+='StepTime=%s\n'%cfgStepTime
    cfgStr+='NTimes=%i\n'%len(acc_files)
    cfgStr+='RightAscension=%s\n'%strSrcRA
    cfgStr+='Declination=%s\n'%strSrcDec
    cfgStr+='WriteAutoCorr=T\n'
    cfgStr+='WriteImagerColumns=T\n'
    ms_out='%s%i.ms'%(ms_out_pre,sbs[merge_sid])
    cfgStr+='MSName=%s\n'%ms_out
    cfgStr+='AntennaTableName=%s\n'%ant_table
    cfgFN='tempMakems.cfg'
    cfgFH=open(cfgFN,'w')
    cfgFH.write(cfgStr)
    cfgFH.close()
    print '\tDone\n'
    
    #run makems
    print 'Running makems...'
    os.system('which makems')
    #os.system('makems %s'%cfgFN)
    os.system('makems %s 2&> test.log'%cfgFN)
    pwd=os.getcwd()
    if os.path.exists(pwd+'/'+ms_out): os.system('rm -r %s/%s'%(pwd,ms_out))
    os.system('mv %s_p0 %s'%(ms_out,ms_out))
    print 'Done\n'
    
    #apply phase correction for multiple files
    nint=len(acc_files)
    if nint>1:
        print 'Applying Phase corrections...'
        for sid,sb in enumerate(sbs):
            tSB=ephemTime[sid::len(sbs)]
            for tid,t in enumerate(tSB):
                if tid==0: continue
                src=phsCenterSrc(obs,t)
                src.compute(obs)
                for i in range(nants):
                    for j in range(nants):
                        if i==j: continue
                        xyz_i, xyz_j = numpy.array(ant_list[i]['POSITION']), numpy.array(ant_list[j]['POSITION'])
                        fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid]=phs2src(fullAcc[2*i:2*i+1,2*j:2*j+1,tid,sid], xyz_i, xyz_j, src, obs, numpy.array(cfgStartFreq))
        print 'done\n'
    
    #fill ms with acc data
    print 'Filling DATA column with Subband data...'
    ms=pt.table(ms_out,readonly=False)
    dc=ms.col('DATA')
    ant1c=ms.col('ANTENNA1')
    ant2c=ms.col('ANTENNA2')
    dc_data=dc.getcol()
    ant1c_data=ant1c.getcol()
    ant2c_data=ant2c.getcol()
    ncorr=len(ant1c_data)/len(acc_files)
    nint=len(acc_files)
    for j in range(nint):
        for i in range(ncorr):
            #ant1=ant1c_data[i]
            #ant2=ant2c_data[i]
            #use top half of correlation matrix
            ant1=ant2c_data[j*ncorr+i]
            ant2=ant1c_data[j*ncorr+i]
            #CORR_TYPE: [9,10,11,12] is [XX,XY,YX,YY]
            for sid,sb in enumerate(sbs):
                dc_data[j*ncorr+i,sid]=[fullAcc[2*ant1,2*ant2,j,sid],
                                        fullAcc[2*ant1,2*ant2+1,j,sid],
                                        fullAcc[2*ant1+1,2*ant2,j,sid],
                                        fullAcc[2*ant1+1,2*ant2+1,j,sid]]
    #correct EXPOSURE and INTERVAL for multiple integration sets
    if nint>1:
        exposure_col=ms.col('EXPOSURE')
        interval_col=ms.col('INTERVAL')
        exp_data=exposure_col.getcol()
        int_data=interval_col.getcol()
        exp_data=numpy.ones_like(exp_data)*opts.int_time
        int_data=numpy.ones_like(int_data)*opts.int_time
        ms.putcol('EXPOSURE', exp_data)
        ms.putcol('INTERVAL', int_data)
    ms.putcol('DATA', dc_data)
    ms.close()
    print 'Done\n'

#meta data clean-up
print 'Metadata cleanup...'
ms=pt.table(ms_out+'/ANTENNA',readonly=False,ack=True)
dd=ms.col('DISH_DIAMETER')
dd_data=dd.getcol()
an=ms.col('NAME')
an_data=an.getcol()
for i in range(len(ant_list)):
    an_data[i]=ant_list[i]['NAME']
    dd_data[i]=ant_list[i]['DISH_DIAMETER']
ms.putcol('NAME',an_data)
ms.putcol('DISH_DIAMETER',dd_data)
ms.close()
print 'Done\n'

