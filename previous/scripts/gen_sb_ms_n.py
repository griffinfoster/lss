#!/usr/bin/env python
"""
Generate single subband measurement set files from ACC files, store in directory structure
"""

import sys,os
import numpy as n

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-F', '--ant_field', dest='ant_field', default=None,
        help='AntennaField.conf file for the LOFAR station of the ACC files, default: None')
    o.add_option('-A', '--ant_array', dest='ant_array', default=None,
        help='AntennaArray.conf file for the LOFAR station geographical coordinates, default: None')
    o.add_option('-r', '--rcumode', dest='rcumode', default='3',
        help='RCU Mode: default: 3 (LBA High)')
    o.add_option('-n','--name', dest='station_name', default='LOFAR',
        help='Station name, default: LOFAR')
    o.add_option('-s', '--subband', dest='subband', default='all',
        help='Select which subband to select out for the MS. Can use a list 200,400,411 or a range 200_210, default:all')
    o.add_option('-t', '--utc', dest='utcOffset', default=0., type='int',
        help='Offset in local clock from UTC where observation was taken as the file timestamp is in local units, e.g. +1 during BST. default:0')
    o.add_option('-D', '--dir', dest='rootdir',
        help='Directory name to store subband MS files, default: <array_type>_<rcumode> in the current directory')
    opts, args = o.parse_args(sys.argv[1:])
acc_files=args

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

acc2ms_path='/home/griffin/node/projects/lss/acc'
acc2ms_path+='/acc2ms.py'
rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
            {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
            {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
            {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
            {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
            {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
            {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
            {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7

#Setup directory structure
subbands=512
sbStr=opts.subband
rcumode=int(opts.rcumode)
if opts.subband.startswith('all'): sbStr='0_%i'%(subbands-1)
if opts.rootdir is None: rootdir=rcuInfo[rcumode]['mode']
else: rootdir=opts.rootdir
pwd=os.getcwd()
fullpath=pwd + '/' + rootdir
print 'Using %s as the root directory'%fullpath
if not os.path.exists(fullpath):
    os.makedirs(fullpath)
sbs=convert_arg_range(sbStr)
for sb in sbs:
    if sb < 10: sbpath='%s/s00%i'%(fullpath,sb)
    elif sb < 100: sbpath='%s/s0%i'%(fullpath,sb)
    else: sbpath='%s/s%i'%(fullpath,sb)
    if not os.path.exists(sbpath): os.makedirs(sbpath)

for fid,fn in enumerate(acc_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(acc_files))
    fnArr=fn.rsplit('/',1)
    shortFn=fnArr[-1]
    shortFnArr=shortFn.split('_')
    ms_prefix='%s/%s_%s_r%i_s'%(fnArr[0],shortFnArr[0],shortFnArr[1],rcumode)
    execStr="-F %s -A %s -r %s -s %s -t %i -o %s %s"%(opts.ant_field,opts.ant_array,rcumode,sbStr,opts.utcOffset,ms_prefix,fn)
    os.system('%s %s'%(acc2ms_path,execStr))
    print "Moving MS files to subband directories"
    for sb in sbs:
        sbMS='%s%i.ms'%(ms_prefix,sb)
        if sb < 10: sbID='00%i'%sb
        elif sb < 100: sbID='0%i'%sb
        else: sbID=str(sb)
        os.system("mv %s %s/s%s/"%(sbMS,fullpath,sbID))

