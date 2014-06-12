#!/usr/bin/env python
"""
Plot an All Sky image of a LOFAR LBA fits file
"""

import sys,os
import pyfits as pf
import numpy as n

import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as p
from mpl_toolkits.basemap import Basemap

def readFITS(fn,hdr=False):
    """Read a FITS image file and returns a numpy array
    """
    hdulist=pf.open(fn)
    im=hdulist[0].data
    hdulist.close()
    if hdr:
        return im[0,0],getFITSInfo(fn)
    else: return im[0,0]

def getFITSInfo(fn):
    """Parse the FITS header for pointing and pixel size information
    return [RA,DEC], pixel resolution, pixel of [RA,DEC]
    """
    hdulist=pf.open(fn)
    hdr=hdulist[0].header
    #CTYPE1: RA---[PROJ], projection SIN/TAN/ARC
    #CRVAL1: reference RA position in degrees
    #CRPIX1: location of reference pixel
    #CDELT1: delta RA/pixel size in degrees
    #CTYPE2: DEC--[PROJ], projection SIN/TAN/ARC
    #CRVAL2: reference DEC position in degrees
    #CRPIX2: location of reference pixel
    #CDELT2: delta DEC/pixel size in degrees
    ra=hdr['CRVAL1']
    dra=hdr['CDELT1']
    raPix=hdr['CRPIX1']
    dec=hdr['CRVAL2']
    ddec=hdr['CDELT2']
    decPix=hdr['CRPIX2']
    hdulist.close()
    return {'ra':ra,'dec':dec,'dra':dra,'ddec':ddec,'raPix':raPix,'decPix':decPix}

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-c','--contour',dest='contour',default=False,action='store_true',
        help='Plot a contour map')
    o.add_option('-m','--mask',dest='mask',action='store_true', default=False,
        help='Mask out the data beyond the horizon limits')
    o.add_option('-s','--savefig',dest='savefig',default=None,
        help='Save the figure to a file')
    o.add_option('-n', '--ncontour', dest='ncontour', default=20, type='int',
        help='Number of contour lines to use in plot, default: 20')
    o.add_option('-t', '--title', dest='title', default='All-Sky LBA Image',
        help='Plot title')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]
    print 'Plotting:',fn

    im,hdr=readFITS(fn,hdr=True)
    print hdr

    fig=p.figure(figsize=(6,6))
    #fig=p.figure()
    #ax=fig.add_axes([0,0,1,1])
    m=Basemap(projection='ortho',lon_0=hdr['ra'],lat_0=hdr['dec'],resolution='l')
    m.drawmapboundary(linewidth=.5)
    parallels=n.arange(-90.,120.,15.)
    m.drawparallels(parallels)
    meridians=n.arange(0.,360.,30.)
    m.drawmeridians(meridians)
   
    x,y=m.makegrid(im.shape[0],im.shape[1])
    pos_thresh=1e10
    im=n.fliplr(im)
    if opts.mask: im=n.ma.masked_where(x+y>pos_thresh,im)
    if opts.contour:
        x0,y0=m(x,y)
        cs=m.contour(x0,y0,im,opts.ncontour)
        print cs.levels
    else:
        m.imshow(im)

    #class A sources
    src_ras=[299.86791,350.84583,83.63333,187.705833]
    src_decs=[40.733888,58.810833,22.01444,12.39111]
    src_names=['CYG','CAS','TAU','VIR']
    ##flip RAs
    #src_ras=hdr['ra']-(n.array(src_ras)-hdr['ra'])
    #src_ras=n.array(src_ras)+7.
    src_ras=n.array(src_ras)
    src_decs=n.array(src_decs)
    #print src_ras,src_decs
    src_x,src_y=m(src_ras,src_decs)
    m.plot(src_x,src_y,'k+')
    #print src_x,src_y
    xoffset=100000
    yoffset=100000
    for sid,sn in enumerate(src_names):
        if src_x[sid] > pos_thresh or src_y[sid] > pos_thresh: continue
        p.text(src_x[sid]+xoffset,src_y[sid]+yoffset,src_names[sid])

    #labels
    p.xlabel('S')
    p.ylabel('E')
    p.title(opts.title)
    #print parallels,meridians
    #RA labels
    raX,raY=m(meridians,n.zeros_like(meridians))
    #print raX,raY
    for lid,label in enumerate(raX):
        if raX[lid] > pos_thresh: continue
        raStr='%ih'%int(meridians[lid]/15.)
        p.text(raX[lid],raY[lid],raStr,size='small')
    #DEC labels
    decX,decY=m(n.ones_like(parallels)*hdr['ra'],parallels)
    for lid,label in enumerate(decY):
        if decY[lid] > pos_thresh: continue
        decStr='%i$^{\circ}$'%int(parallels[lid])
        p.text(decX[lid]+xoffset,decY[lid]+yoffset,decStr,size='small')

    if opts.savefig is not None:
        p.savefig(opts.savefig)
    p.show()

