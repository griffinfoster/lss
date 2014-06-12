#!/usr/bin/env python
"""
Perform a DFT on ACC data to form a dirty image
"""

#toDo
#multiple subbands
#multiple times/ACCs
#include w component as phase
#polarizations
#plot: psf, unclib, calib, uv

import numpy
import ephem
import pylab

import sys,os
import struct
import time

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
def xyz2uvw(xyz, src, obs, f):
    """Return an array of UVW values"""
    uvw=numpy.zeros((xyz.shape[0],xyz.shape[0],3))
    for i in range(xyz.shape[0]):
        for j in range(xyz.shape[0]):
            if i==j: continue
            uvw[i,j]=gen_uvw(xyz[i], xyz[j], src, obs, f)[:,0,0]
    return uvw
def dft2(d,k,l,u,v):
    """compute the 2d DFT for position (k,l) based on (d,uvw)"""
    return numpy.sum(d*numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
    #psf:
    #return numpy.sum(numpy.exp(-2.*numpy.pi*1j*((u*k) + (v*l))))
def dftImage(d,uvw,px,res,mask=False):
    """return a DFT image"""
    nants=uvw.shape[0]
    im=numpy.zeros((px[0],px[1]),dtype=complex)
    mid_k=int(px[0]/2.)
    mid_l=int(px[1]/2.)
    u=uvw[:,:,0]
    v=uvw[:,:,1]
    w=uvw[:,:,2]
    u/=mid_k
    v/=mid_l
    start_time=time.time()
    for k in range(px[0]):
        for l in range(px[1]):
            im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
            if mask:        #mask out region beyond field of view
                rad=(((k-mid_k)*res)**2 + ((l-mid_l)*res)**2)**.5
                if rad > mid_k*res: im[k,l]=0
                #else: im[k,l]=dft2(d,(k-mid_k),(l-mid_l),u,v)
    print time.time()-start_time
    return im

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
    accSize=fnArr[3].split('x')
    nsb=int(accSize[0])
    accComp=numpy.fromfile(fn,dtype='complex')
    accComp=accComp.reshape(nsb,nantpol,nantpol)
    for sid,sb in enumerate(sbs):
        accSB=accComp[sb]
        fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(fnArr[1][0:2]),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
        fTime=time.mktime(fTimeStruct)
        fTime-=(subbands-sb)*int_time #fTimeStruct is the time for the last subband (511)
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
polAcc=fullAcc[::2,::2,:]

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
        uvw=xyz2uvw(xyz, src, obs, freqs[sid])
        #uvw plot
        #pylab.plot(uvw[:,:,0],uvw[:,:,1],'.')
        im=dftImage(polAcc[:,:,fid,sid],uvw,px,res,mask=True)
        im=im.real
        #normalize to peak value
        #im/=numpy.max(im)
        im=numpy.log(im)
        pylab.imshow(im)
        pylab.xlabel('Pixels (E-W)')
        pylab.ylabel('Pixels (N-S)')
        pylab.title('DFT Image')
        #pylab.contour(im)
if opts.savefig: pylab.savefig('dft_s%i.png'%sbs[0])
pylab.show()

