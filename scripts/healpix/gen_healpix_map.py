#!/usr/bin/env python
"""
Generate healpix maps using AIPY's mk_map.py
"""

import sys,os
import numpy as n

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_FILES')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])
fits_files=args

mk_map_path='/home/griffin/node/projects/lss/healpix'
mk_map_path+='/mk_map.py'
nside=128
fwidth=400.

for fid,fn in enumerate(fits_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(fits_files))
    fnBase=fn.split('/')[-1]
    fnArr=fnBase.split('.')
    healpix_fn="%s.heal.fits"%(fnArr[0])
    execStr="--nside=%i --fwidth=%f -n -m %s %s"%(nside,fwidth,healpix_fn,fn)
    print execStr
    os.system('%s %s'%(mk_map_path,execStr))

