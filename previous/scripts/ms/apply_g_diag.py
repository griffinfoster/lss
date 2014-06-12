#!/usr/bin/env python
"""
Apply a set of antenna gains terms to LOFAR single station measurement sets
"""

import pylab as p
import numpy as n
import sys,os
import cPickle as pickle
from pyrap.tables import *

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS')
    o.set_description(__doc__)
    o.add_option('-g', '--gains', dest='gain_fn', default=None,
        help='Gain pickle file to be applied to DATA and written to CORRECTED_DATA, default: None')
    o.add_option('-m', '--mode', dest='mode', default='mean',
        help='Gain values mode computed which should be applied: mean, median, default: mean')
    opts, args = o.parse_args(sys.argv[1:])

    #load gains pickle file
    fh=open(opts.gain_fn,'rb')
    gains=pickle.load(fh)
    fh.close()
    gxx=gains[opts.mode+'_xx']
    gyy=gains[opts.mode+'_yy']
    print gains['args']
    
    newcolname='CORRECTED_DATA'
    for i in args:
        ms = table(i,readonly=False)
        data=ms.getcol('DATA')
        ant1=ms.getcol('ANTENNA1')
        ant2=ms.getcol('ANTENNA2')

        if not(newcolname in ms.colnames()):
            print '%s column does not exist, creating one...'%newcolname,
            coldmi=ms.getdminfo('DATA')
            coldmi["NAME"]=newcolname
            ms.addcols( maketabdesc( makearrcoldesc(newcolname,data,shape=[1,4],valuetype='complex') ), coldmi)
            print 'done'

        #generate corrected data
        cdata=n.zeros_like(data)
        for idx,d in enumerate(data):
            cdata[idx]=n.array([d[0,0] / (gxx[ant1[idx]]*n.conjugate(gxx[ant2[idx]])),
                                d[0,1] / (gxx[ant1[idx]]*n.conjugate(gyy[ant2[idx]])),
                                d[0,2] / (gyy[ant1[idx]]*n.conjugate(gxx[ant2[idx]])),
                                d[0,3] / (gyy[ant1[idx]]*n.conjugate(gyy[ant2[idx]])) ])
        ms.putcol(newcolname,cdata)
        ms.close()


