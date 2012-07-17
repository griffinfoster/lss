#!/usr/bin/env python

import numpy
import pylab
import time,sys,struct,math

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-c', dest='chan_index', default='all',
        help='Select which channels to plot. Options are <ch_i>,...,<ch_j>, or a range <ch_i>_<ch_j>. Default=all')
    o.add_option('-d', '--decimate', dest='decimate', default=1,
        help='Decimate in time by N samples to speed up plotting, Default=None')
    o.add_option('-m', '--mode', dest='mode', default='log',
        help='Plotting mode: linear, log, default:log')
    o.add_option('-t', '--time', dest='time', default='all',
        help='Select which time samples to plot, <t_i> or <t_i>,<t_j>,... or if in waterfall mode <t_0>_<t_k>. Default: all times')
    o.add_option('-w','--water', dest='water', default=False, action='store_true',
        help='Produce a waterfall plot of a time range using -t <t_0>_<t_n>')
    o.add_option('--chan', dest='chan_time', action='store_true',
        help='Plot individual channels as a function of time')
    o.add_option('--share', dest='share', action='store_true',
        help='Share plots in a single frame.')
    o.add_option('--nolegend', dest='nolegend', default=False, action='store_true',
        help='Turn off the legend')
    opts, args = o.parse_args(sys.argv[1:])
sst_files=args

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

decimate=int(opts.decimate)
subbands=512
m2 = int(math.sqrt(len(sst_files)))
m1 = int(math.ceil(float(len(sst_files)) / m2))
for fid,fn in enumerate(sst_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(sst_files))
    fn_parse=fn.split('/')[-1]
    fnArr=fn_parse.split('_')
    fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(fnArr[1][0:2]),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
    rcuID=fnArr[3].split('.')
    rcuID=rcuID[0]
    fh=open(fn,'rb')
    sst_str=fh.read()
    nflts=len(sst_str)/8 #number of 8 byte floats
    sst=struct.unpack('<%id'%nflts, sst_str)
    sst=numpy.array(sst)
    if opts.mode.startswith('log'):
        sst=numpy.log10(sst)
    fh.close()
    autos=len(sst)/subbands
    
    if opts.chan_index.startswith('all'): chans = numpy.arange(subbands)
    else: chans=numpy.array(convert_arg_range(opts.chan_index))
    if opts.time.startswith('all'): ts = numpy.arange(autos)
    else: ts=numpy.array(convert_arg_range(opts.time))
    if decimate>1:ts=ts[::decimate]
    
    sst=numpy.reshape(sst,(autos,subbands))
    sst=sst[:,chans]
    sst=sst[ts]
    
    if not opts.share:
        pylab.subplot(m2, m1, fid+1)
        pylab.title(rcuID)
    if opts.water:
        pylab.imshow(sst, aspect='auto')
        pylab.colorbar()
    else:
        if opts.chan_time:
            for cid,c in enumerate(chans):
                pylab.plot(sst[:,cid],label='%s:%s'%(rcuID,c))
        else:
            for tid,t in enumerate(ts):
                pylab.plot(sst[tid],label='%s:%s'%(rcuID,t))
if not opts.water and not opts.nolegend:
    pylab.legend()
pylab.show()
