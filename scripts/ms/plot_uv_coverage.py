#!/usr/bin/python
"""Extract UVW positions and plot them"""

import numpy as n
import pylab as p
import sys
from pyrap.tables import *

if __name__ == '__main__':                                                                                                                          
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS_FILE')
    o.set_description(__doc__)
    o.add_option('-r', '--rcumode', dest='rcumode', default=3, type='int',
        help='Station RCU Mode, usually 3,5,6,7, default: 3(LBA High)')
    o.add_option('-d', '--uvdist', dest='uvdist', action='store_true',
        help='Plot a histor gram of the number of samples based on their UV distance')
    o.add_option('-S', '--savefig', dest='savefig', action='store_true',
        help='Save figure as a PNG')
    opts, args = o.parse_args(sys.argv[1:])

    rcuInfo = [ {'mode':'OFF', 'rcuID':0, 'array_type':'LBA', 'bw':100000000.},            #0
                {'mode':'LBL_HPF10MHZ', 'rcuID':1, 'array_type':'LBA', 'bw':100000000.},   #1
                {'mode':'LBL_HPF30MHZ', 'rcuID':2, 'array_type':'LBA', 'bw':100000000.},   #2
                {'mode':'LBH_HPF10MHZ', 'rcuID':3, 'array_type':'LBA', 'bw':100000000.},   #3
                {'mode':'LBH_HPF30MHZ', 'rcuID':4, 'array_type':'LBA', 'bw':100000000.},   #4
                {'mode':'HBA_110_190MHZ', 'rcuID':5, 'array_type':'HBA', 'bw':100000000.}, #5
                {'mode':'HBA_170_230MHZ', 'rcuID':6, 'array_type':'HBA', 'bw':100000000.}, #6
                {'mode':'HBA_210_290MHZ', 'rcuID':7, 'array_type':'HBA', 'bw':100000000.}] #7
    opts, args = o.parse_args(sys.argv[1:])

    for i in args:
        ms = table(i,readonly=False)
        uvw=ms.getcol('UVW')
        ms.close()

    if opts.uvdist:
        uvdist=n.sqrt((uvw[:,0]**2)+(uvw[:,1]**2))
        p.hist(uvdist,bins=75,alpha=.5)
        p.xlabel('UV Distance (m)')
        p.ylabel('# Baselines')
        p.title('UV Distribution')
        if opts.rcumode<5:
            figname='lba_uv_distribution.pdf'
        else:
            figname='hba_uv_distribution.pdf'
    else:
        fig=p.figure(figsize=(8.5,8))
        axes=p.axes()
        if opts.rcumode<5:
            p.plot(uvw[:,0],uvw[:,1],'k,')
            p.title('LBA Snapshot UV Coverage')
            p.xlim(-70,70)
            p.ylim(-70,70)
            figname='lba_uv_coverage.pdf'
        else:
            p.plot(uvw[:,0],uvw[:,1],'k,')
            p.plot(-1*uvw[:,0],-1*uvw[:,1],'k,')
            p.title('HBA Snapshot UV Coverage')
            figname='hba_uv_coverage.pdf'
        
        p.xlabel("U (m)")
        p.ylabel("V (m)")
        ax=p.gca()
        p.grid(True)

    if opts.savefig: p.savefig(figname)
    p.show()

