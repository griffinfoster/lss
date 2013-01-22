#! /usr/bin/env python
"""
This is a general-purpose script for making faceted, spherical maps (stored in
Healpix FITS format) from individual "flat" maps stored in FITS files.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse, ephem
import pylab as p

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.fits')
o.set_description(__doc__)
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
o.add_option('-n', '--nobeam', dest='nobeam', action='store_true',
    help='Do not apply a Gaussian downweighting')
o.add_option('--min_alt', dest='min_alt', default=20.0, type='float',
    help='Mask out pixels below a minimum altitude (degrees), default=20.0')
o.add_option('--nside', dest='nside', type='int', default=256,
    help='NSIDE parameter for map, if creating a new file.')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
opts, args = o.parse_args(sys.argv[1:])

# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=opts.nside)
skymap.set_interpol(opts.interpolate)

prev_dra = None
for i, filename in enumerate(args):
    img, kwds = a.img.from_fits(filename)
    img = img.squeeze()
    #transpose and flip image to correct direction
    img=n.fliplr(n.flipud(img.T))
    ##shift the image so the minimum is zero
    print 'MIN',n.min(img),n.max(img)
    #img+=20000.
    #img=img-n.min(img)
    # Read ra/dec of image center, which are stored in J2000
    assert(kwds['epoch'] == 2000)
    s = ephem.Equatorial(kwds['ra']*a.img.deg2rad, kwds['dec']*a.img.deg2rad, 
        epoch=ephem.J2000)
    # To precess the entire image to J2000, we actually need to precess the
    # coords of the center back to the epoch of the observation (obs_date),
    # get the pixel coordinates in that epoch, and then precess the coords of
    # each pixel.  This is because extrapolating pixel coords from the J2000 
    # center assumes the J2000 ra/dec axes, which may be tilted relative to
    # the ra/dec axes of the epoch of the image.
    print kwds['obs_date']
    obs_date=kwds['obs_date']
    dd,tt=obs_date.split('T')
    obs_date=dd+' '+tt
    #s = ephem.Equatorial(s, epoch=kwds['obs_date'])
    s = ephem.Equatorial(s, epoch=obs_date)
    ra, dec = s.get()
    print '-----------------------------------------------------------'
    print 'Reading file %s (%d / %d)' % (filename, i + 1, len(args))
    print kwds
    print 'Pointing (ra, dec):', ra, dec
    print 'Image Power:', n.abs(img).sum()
    if prev_dra != kwds['d_ra']:
        prev_dra = kwds['d_ra']
        DIM = img.shape[0]
        RES = 1. / (kwds['d_ra'] * a.img.deg2rad * DIM)
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        # Define a weighting for gridding data into the skymap
        map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
        if opts.nobeam: map_wgts=n.ones_like(map_wgts)
        #add a mask to pixels below minimum altitude
        min_alt=opts.min_alt    #degrees
        mask_ind=n.where(tz < n.sin(min_alt*a.img.deg2rad))
        map_wgts.mask[mask_ind]=True
        #p.imshow(tz)
        #p.colorbar()
        #p.subplot(121)
        #p.imshow(map_wgts)
        #p.colorbar()
        #p.subplot(122)
        #p.imshow(img)
        #p.colorbar()
        #p.imshow(img*map_wgts)
        #p.colorbar()
        map_wgts.shape = (map_wgts.size,)
        valid = n.logical_not(map_wgts.mask)
        map_wgts = map_wgts.compress(valid)
    # Get coordinates of image pixels in the epoch of the observation
    ex,ey,ez = im.get_eq(ra, dec, center=(DIM/2,DIM/2))
    ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
    img = img.flatten(); img = img.compress(valid)
    #p.show()
    # Precess the pixel coordinates to the (J2000) epoch of the map
    #m = a.coord.convert_m('eq','eq', 
    #    iepoch=kwds['obs_date'], oepoch=ephem.J2000)
    m = a.coord.convert_m('eq','eq', 
        iepoch=obs_date, oepoch=ephem.J2000)
    ex,ey,ez = n.dot(m, n.array([ex,ey,ez])) 
    # Put the data into the skymap
    skymap.add((ex,ey,ez), map_wgts, img)
skymap.to_fits(opts.map, clobber=True)

