#!/usr/bin/env python
"""
Use meqtrees-pipeliner to apply a self calibration for LOFAR single subband snapshot data
"""

import sys,os

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] DIRECTORY/MS')
    o.set_description(__doc__)
    o.add_option('-s', '--subband', dest='subband', default='all',
        help='Select which subbands for image generation. Can use a list 200,400,411 or a range 200_210, default:all')
    o.add_option('-p', '--profile', dest='profile',
        help='MeqTrees TDL profile to use')
    o.add_option('-S', '--section', dest='section',
        help='Section name of the TDL profile to use')
    o.add_option('-P', '--prefix', dest='prefix',
        help='Prefix of files to process')
    opts, args = o.parse_args(sys.argv[1:])
rootdir=args[0]

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

mt=8    #number of core for meqtrees multithreading
pwd=os.getcwd()
fullpath='%s/%s'%(pwd,rootdir)

plpath='/home/griffin/meqtrees/trunk/Timba/PyApps/src'
plpath+='/meqtree-pipeliner.py'
calpath='/home/griffin/meqtrees/trunk/Cattery/Calico'
calpath+='/calico-wsrt-tens-lofar.py'

sbStr=opts.subband
subbands=512
if opts.subband.startswith('all'): sbStr='0_%i'%(subbands-1)
sbs=convert_arg_range(sbStr)

subdirs=os.listdir(fullpath)
for s in subdirs:
    snum=int(s[1:])
    if snum in sbs:
        fns=os.listdir('%s/%s'%(fullpath,s))
        for fn in fns:
            if fn.startswith(opts.prefix):
                filepath='%s/%s/%s'%(fullpath,s,fn)
                execStr='%s -c %s --mt %i @%s ms_sel.msname=%s %s =cal_G_diag.G_diag_clear_meptables =cal_G_diag =make_dirty_image'%(plpath, opts.profile, mt, opts.section,filepath,calpath)
                print execStr
                os.system(execStr)

