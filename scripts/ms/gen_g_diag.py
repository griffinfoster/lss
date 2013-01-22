#!/usr/bin/env python
"""
Generate a set of antenna gains terms based on an average from MeqTrees calibrated LOFAR single station measurement sets, store to a pickle file
"""

import pylab as p
import numpy as n
import Owlcat.ParmTables
import sys,os
import cPickle as pickle

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    nants=96    #HARDCODE:LOFAR station elements

    g_diag_xx={}
    g_diag_yy={}
    for a in range(nants):
        g_diag_xx['ant%i'%a]=[]
        g_diag_yy['ant%i'%a]=[]
    for ms in args:
        print 'Accessing:', ms
        fmep=ms+'/'+'G_diag.fmep'
        ptab=Owlcat.ParmTables.ParmTab(fmep)
        for a in range(nants):
            rr=ptab.funkset('G:%i:xx:r'%a).array()
            ii=ptab.funkset('G:%i:xx:i'%a).array()
            g_diag_xx['ant%i'%a].append(rr.flatten()[0]+ii.flatten()[0]*1j)
            rr=ptab.funkset('G:%i:yy:r'%a).array()
            ii=ptab.funkset('G:%i:yy:i'%a).array()
            g_diag_yy['ant%i'%a].append(rr.flatten()[0]+ii.flatten()[0]*1j)

    gain_xx_mean=[]
    gain_yy_mean=[]
    gain_xx_median=[]
    gain_yy_median=[]
    
    for a in range(nants):
        gains_xx=n.array(g_diag_xx['ant%i'%a])
        gains_yy=n.array(g_diag_yy['ant%i'%a])
        gain_xx_mean.append(n.mean(gains_xx))
        gain_yy_mean.append(n.mean(gains_yy))
        gain_xx_median.append(n.median(gains_xx))
        gain_yy_median.append(n.median(gains_yy))

    gain_xx_mean=n.array(gain_xx_mean)
    gain_yy_mean=n.array(gain_yy_mean)
    gain_xx_median=n.array(gain_xx_median)
    gain_yy_median=n.array(gain_yy_median)

    ms_base=args[0].split('/')[-1]
    bs_base=ms_base.split('.')[0]
    date=ms_base.split('_')[0]
    sb=ms_base.split('_')[-1]
    pfn='gains_'+sb+'_'+date+'.pkl'
    gains={ 'mean_xx':gain_xx_mean,
            'mean_yy':gain_yy_mean,
            'median_xx':gain_xx_median,
            'median_yy':gain_yy_median,
            'args':args}
    fh=open(pfn,'wb')
    pickle.dump(gains,fh)
    fh.close()
    print 'Gains written to',pfn

