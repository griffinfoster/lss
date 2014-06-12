#!/usr/bin/python
"""
Plot the correlation matrix
"""

import numpy as n
import pylab as p
import sys
from pyrap.tables import *

if __name__ == '__main__':                                                                                                                          
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS_FILE')
    o.set_description(__doc__)
    o.add_option('-m', '--mode', dest='mode', default='lin',
        help='Plotting mode: linear,log,real,imag,phs, default:linear')

    opts, args = o.parse_args(sys.argv[1:])

    colname='CORRECTED_DATA'

    ms = table(args[0],readonly=False)
    ant1=ms.getcol('ANTENNA1')
    ant2=ms.getcol('ANTENNA2')
    data=ms.getcol(colname)

    nantpol=96*2    #HARDCODE for LOFAR international station
    corr=n.zeros((nantpol,nantpol),dtype=complex)
    for i,ant in enumerate(ant1):
        corr[2*ant1[i],2*ant2[i]]=data[i,0,0]
        corr[2*ant1[i],2*ant2[i]+1]=data[i,0,1]
        corr[2*ant1[i]+1,2*ant2[i]]=data[i,0,2]
        corr[2*ant1[i]+1,2*ant2[i]+1]=data[i,0,3]
    ms.close()

    #flag autos
    for i in range(nantpol): corr[i,i]=0

    print n.max(corr)

    if opts.mode.startswith('real'): pwrSB=corr.real
    elif opts.mode.startswith('imag'): pwrSB=corr.imag
    elif opts.mode.startswith('lin'): pwrSB=n.abs(corr)
    elif opts.mode.startswith('log'): pwrSB=n.log(n.abs(corr))
    elif opts.mode.startswith('phs'): pwrSB=n.angle(corr)
    
    p.imshow(pwrSB)
    p.colorbar()
    p.show()
    
