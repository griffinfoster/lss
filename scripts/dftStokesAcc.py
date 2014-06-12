#!/usr/bin/env python
"""
Perform a DFT on ACC data to form a Stokes dirty image
"""

import numpy
import ephem
import pylab

import sys,os
import struct
import time

import lss

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help='AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help='AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help='Station RCU Mode, usually 3,5,6,7, default: 3(LBA High)')
    o.add_option('-s', '--subband', dest='subband', default='0',
        help='Select which subband to image, default:0')
    o.add_option('-p', '--pixels', dest='pixels', default=128, type='int',
        help='Width of image in pixels, default: 128')
    o.add_option('-C','--cal',dest='calfile',default=None,
        help='Apply a calibration soultion file to the data.')
    o.add_option('-S','--save',dest='savefig',action='store_true',
        help='Save the figure as a png')
    opts, args = o.parse_args(sys.argv[1:])
acc_files=args

rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7

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

#Parse AntennaArray file
fh=open(opts.ant_array)
rcumode=opts.rcumode
readLatLon=False
for line in fh:
    if line.lower().startswith(rcuInfo[rcumode]['array_type'].lower()):
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
    if line.lower().startswith(rcuInfo[rcumode]['array_type'].lower()):
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
sbs=convert_arg_range(opts.subband)
subbands=512
bw=rcuInfo[rcumode]['bw']
df=(bw/subbands) #100 MHz / 512 subbands
avg_freq=numpy.average(numpy.array(sbs))*df
freqs=numpy.array(sbs)*df
if rcumode==5: freqs+=100000000.
elif rcumode==6: freqs+=160000000.
elif rcumode==7: freqs+=200000000.
nants=len(ants)
int_time=1.0

#Select subband from each ACC file
print 'Selecting Out Subband(s) %s from ACC Files...'%(str(sbs))
npols=2
nantpol=nants*npols
fullAcc=numpy.zeros((nantpol,nantpol,len(acc_files)*len(sbs)),dtype=complex)
for fid,fn in enumerate(acc_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(acc_files))
    fnBase=fn.split('/')[-1]
    fnArr=fnBase.split('_')
    #accSize=fnArr[3].split('x')
    #nsb=int(accSize[0])
    accComp=numpy.fromfile(fn,dtype='complex')
    nsb=accComp.shape[0]/(192*192)
    accComp=accComp.reshape(nsb,nantpol,nantpol)
    for sid,sb in enumerate(sbs):
        #accSB=accComp[sb]
        if nsb==1: accSB=accComp[0]
        else: accSB=accComp[sb]
        fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(fnArr[1][0:2]),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
        fTime=time.mktime(fTimeStruct)
        #fTime-=(subbands-sb)*int_time #fTimeStruct is the time for the last subband (511)
        if fid==0 and sid==0:
            fullAcc[:,:,0]=accSB
            initTimeStruct=time.localtime(fTime)
            ephemTime=[time.strftime('%Y/%m/%d %H:%M:%S',initTimeStruct)]
            fullTimes=[fTime]
        else:
            fullAcc[:,:,fid*len(sbs)+sid]=accSB
            initTimeStruct=time.localtime(fTime)
            ephemTime.append(time.strftime('%Y/%m/%d %H:%M:%S',initTimeStruct))
            fullTimes.append(fTime)
print 'Done\n'
fullTimes=numpy.array(fullTimes)
fullTimes=fullTimes.reshape(len(acc_files),fullTimes.shape[0]/len(acc_files))
fullAcc=fullAcc.reshape(fullAcc.shape[0],fullAcc.shape[1],len(acc_files),fullAcc.shape[2]/len(acc_files))
ephemTime=numpy.array(ephemTime)
ephemTime=ephemTime.reshape(len(acc_files),len(sbs))

obs=ephem.Observer()
obs.long=lon
obs.lat=lat
obs.elevation=float(elev)
obs.epoch=2000.0
xyz=[]
for a in ants: xyz.append([a[0]+float(arr_x),a[1]+float(arr_y),a[2]+float(arr_z)])
xyz=numpy.array(xyz)

pixels=opts.pixels
px=[pixels,pixels]
fov=numpy.pi    #Field of View in radians
res=fov/px[0]   #pixel resolution
print 'res:',res
pol=0
xxAcc=fullAcc[0::2,0::2,:]
xyAcc=fullAcc[0::2,1::2,:]
yxAcc=fullAcc[1::2,0::2,:]
yyAcc=fullAcc[1::2,1::2,:]

#read cal file if included
if not (opts.calfile is None):
    cal=numpy.fromfile(opts.calfile, dtype='complex')
    cal=numpy.reshape(cal,(subbands,nants,npols))
    calx=cal[:,:,0]
    caly=cal[:,:,1]
    for sid,sb in enumerate(sbs):
        calxSB=calx[sb]
        calxSB=numpy.reshape(calxSB,(96,1))
        gains=numpy.conj(calxSB) * numpy.transpose(calxSB)
        polAcc=numpy.multiply(polAcc,gains)

earth_radius=(xyz[0,0]**2+xyz[0,1]**2+xyz[0,2]**2)**.5
for fid,fn in enumerate(acc_files):
    for sid,sb in enumerate(sbs):
        obs.date=ephemTime[fid,sid]
        src=ephem.FixedBody()
        src._ra=obs.sidereal_time()
        src._dec=obs.lat
        src.compute(obs)
        uvw=lss.xyz2uvw(xyz, src, obs, freqs[sid])
        #uvw plot
        #pylab.plot(uvw[:,:,0],uvw[:,:,1],'.')
        stokesAcc=numpy.array([xxAcc[:,:,fid,sid],xyAcc[:,:,fid,sid],yxAcc[:,:,fid,sid],yyAcc[:,:,fid,sid]])
        stokesIm=lss.dftImage(stokesAcc,uvw,px,res,mask=False,weighting=False,stokes=True)

        pylab.subplot(2,2,1)
        pylab.imshow(stokesIm[:,:,0])
        pylab.xlabel('Pixels (E-W)')
        pylab.ylabel('Pixels (N-S)')
        pylab.title('I')
        pylab.colorbar()
        pylab.subplot(2,2,2)
        pylab.imshow(stokesIm[:,:,1])
        pylab.xlabel('Pixels (E-W)')
        pylab.ylabel('Pixels (N-S)')
        pylab.title('Q')
        pylab.colorbar()
        pylab.subplot(2,2,3)
        pylab.imshow(stokesIm[:,:,2])
        pylab.xlabel('Pixels (E-W)')
        pylab.ylabel('Pixels (N-S)')
        pylab.title('U')
        pylab.colorbar()
        pylab.subplot(2,2,4)
        pylab.imshow(stokesIm[:,:,3])
        pylab.xlabel('Pixels (E-W)')
        pylab.ylabel('Pixels (N-S)')
        pylab.title('V')
        pylab.colorbar()

if opts.savefig: pylab.savefig('dft_s%i.png'%sbs[0])
pylab.show()

