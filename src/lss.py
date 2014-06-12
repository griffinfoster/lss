"""LOFAR Single Station Functions
"""

#TODO: fft convolution functions
#TODO: fft weighting
#TODO: fft scaling
#TODO: save with hermitian

import numpy
import ephem
import pylab

import sys,os
import struct
import time

"""
FT Functions
"""

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

def dft2(d,m,l,u,v,psf=False):
    """compute the 2d DFT for position (m,l) based on (d,uvw)"""
    if psf: return numpy.sum(numpy.exp(-2.*numpy.pi*1j*((u*m) + (v*l))))/u.size
    else: return numpy.sum(d*numpy.exp(-2.*numpy.pi*1j*((u*m) + (v*l))))/u.size

def dftImage(d,uvw,px,res,mask=False,weighting=False,stokes=False):
    """return a DFT image"""
    nants=uvw.shape[0]
    if stokes: im=numpy.zeros((px[0],px[1],4),dtype=complex)
    else: im=numpy.zeros((px[0],px[1]),dtype=complex)
    maskIm=numpy.zeros((px[0],px[1]),dtype=bool)
    mid_m=int(px[0]/2.)
    mid_l=int(px[1]/2.)
    u=numpy.array(uvw[:,:,0])
    v=numpy.array(uvw[:,:,1])
    w=numpy.array(uvw[:,:,2])
    u/=mid_m
    v/=mid_l
    start_time=time.time()
    for m in range(px[0]):
        for l in range(px[1]):
            #weigthing to account for missing numpy.sqrt(1-l^2-m^2) in flat-field approximation
            if weighting: scale=numpy.sqrt(1.-((float(l-mid_l)/mid_l)**2.)-((float(m-mid_m)/mid_m)**2.))
            else: scale=1.
            if stokes:
                im[m,l,0]=dft2(d[0],(m-mid_m),(l-mid_l),u,v)*scale
                im[m,l,1]=dft2(d[1],(m-mid_m),(l-mid_l),u,v)*scale
                im[m,l,2]=dft2(d[2],(m-mid_m),(l-mid_l),u,v)*scale
                im[m,l,3]=dft2(d[3],(m-mid_m),(l-mid_l),u,v)*scale
            else: im[m,l]=dft2(d,(m-mid_m),(l-mid_l),u,v)*scale
            if mask:        #mask out region beyond field of view
                rad=(((m-mid_m)*res)**2 + ((l-mid_l)*res)**2)**.5
                if rad > mid_m*res: maskIm[m,l]=True
    if stokes:
        stokesIm=numpy.zeros_like(im,dtype=float)
        stokesIm[:,:,0]=(im[:,:,0]+im[:,:,3]).real
        stokesIm[:,:,1]=(im[:,:,0]-im[:,:,3]).real
        stokesIm[:,:,2]=(im[:,:,1]+im[:,:,2]).real
        stokesIm[:,:,3]=(im[:,:,1]-im[:,:,2]).imag
        im=stokesIm
    im=numpy.flipud(numpy.rot90(im))
    print time.time()-start_time
    if mask: return im,maskIm
    else: return im

def fftImage(d,uvw,px,res,mask=False,conv='fast'):
    """Grid visibilities and perform an FFT to return an image"""
    start_time=time.time()
    
    mid_m=int(px[0]/2.)
    mid_l=int(px[1]/2.)
    u=numpy.array(uvw[:,:,0])
    v=numpy.array(uvw[:,:,1])
    w=numpy.array(uvw[:,:,2])
    
    gridVis=numpy.zeros((px[0],px[1]),dtype=complex)
    uflat=u.flatten()
    vflat=v.flatten()
    dflat=d.flatten()

    #umax=px[0]/4.3
    #vmax=px[1]/4.3
    umax=px[0]/5.
    vmax=px[1]/5.
    #print umax
    udelta=umax*2./px[0]
    vdelta=vmax*2./px[1]
    gridUV=numpy.mgrid[-1.*umax:umax:(2*umax/px[0]),-1.*vmax:vmax:(2*vmax/px[1])]
    gridVis=numpy.zeros((px[0],px[1]),dtype=complex)
    if conv.startswith('fast'):
        for did,dd in enumerate(dflat):
            #simple, rectangular convolution function
            uu=int(uflat[did]/udelta)
            vv=int(vflat[did]/vdelta)
            gridVis[uu,vv]+=dd
    else:
        convFunc=convRect(udelta,vdelta)
        #stopped here

    im=numpy.fft.fftshift(numpy.fft.fft2(gridVis))
    im=numpy.rot90(numpy.fliplr(im))
    
    print time.time()-start_time
    if mask: return im,maskIm
    else: return im

def convRect(udelta,vdelta):
    """Return a rectangular convolution function"""
    return lambda uu,vv: (1./udelta)*(1./vdelta) if (uu/udelta <= .5 and vv/vdelta <= .5) else 0.

def convGauss(res,alpha=0.75):
    """Return a Gaussian convolution function"""
    return lambda uu,vv: ((1./(alpha*res*numpy.sqrt(numpy.pi)))**2.)*numpy.exp(((-1.*uu)/(alpha*res))**2.)*numpy.exp(((-1.*vv)/(alpha*res))**2.)
