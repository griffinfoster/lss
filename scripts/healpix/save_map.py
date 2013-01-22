#!/usr/bin/env python
"""
Save a HEALPIX FITS map
"""

import sys,os
import numpy as n

import healpy as h

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_MAPS')
    o.set_description(__doc__)
    o.add_option('-s','--save',dest='savemap',default=None,
        help='Save map to file')
    opts, args = o.parse_args(sys.argv[1:])

    m=None
    w=None
    for fn in args:
        print 'Opening:',fn
        if m is None: m,w,hdr=h.read_map(fn,field=(0,1),h=True)
        else:
            m0,w0,hdr=h.read_map(fn,field=(0,1),h=True)
            m+=m0
            w+=w0

    if not(opts.savemap is None):
        ofn=opts.savemap
        h.write_map(ofn,[m,w,w],dtype=n.float64,coord='E')

