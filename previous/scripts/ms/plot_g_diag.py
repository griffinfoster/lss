#!/usr/bin/env python
"""
Plot complex gain terms from MeqTrees calibrated LOFAR single station measurement sets
"""

import pylab as p
import numpy as n
import Owlcat.ParmTables
import sys,os

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    nants=96    #HARDCODE:LOFAR station elements

    g_diag={}
    for a in range(nants):
        g_diag['ant%i'%a]=[]
    for mid,ms in enumerate(args):
        print mid, 'Accessing:', ms
        fmep=ms+'/'+'G_diag.fmep'
        ptab=Owlcat.ParmTables.ParmTab(fmep)
        for a in range(nants):
            rr=ptab.funkset('G:%i:yy:r'%a).array()
            ii=ptab.funkset('G:%i:yy:i'%a).array()
            g_diag['ant%i'%a].append(rr.flatten()[0]+ii.flatten()[0]*1j)
    for a in range(nants):
        gains=n.array(g_diag['ant%i'%a])
        p.subplot(211)
        p.plot(n.abs(gains),label='ant%i'%a)
        #p.plot(gains.real,label='ant%i'%a)
    for a in range(nants):
        gains=n.array(g_diag['ant%i'%a])
        p.subplot(212)
        p.plot(n.angle(gains),label='ant%i'%a)
        #p.plot(gains.imag,label='ant%i'%a)
    p.show()

