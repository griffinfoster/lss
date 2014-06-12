#!/usr/bin/env python

import numpy
import pylab
import time,sys,struct

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-s', '--subband', dest='subband', default=0, type='int',
        help='Subband (0-512), default:0')
    o.add_option('-m', '--mode', dest='mode', default='lin',
        help='Plotting mode: linear,log,real,imag,phs, default:linear')
    opts, args = o.parse_args(sys.argv[1:])
acc_files=args

for fid,fn in enumerate(acc_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(acc_files))
    fnBase=fn.split('/')[-1]
    fnArr=fnBase.split('_')
    accSize=fnArr[3].split('x')
    accComp=numpy.fromfile(fn,dtype='complex')
    nantpol=192
    nsb=accComp.shape[0]/(nantpol*nantpol)
    accComp=accComp.reshape(nsb,nantpol,nantpol)
    if nsb==1: accSB=accComp[0]
    else: accSB=accComp[opts.subband]
    
    fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(fnArr[1][0:2]),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
    fTime=time.mktime(fTimeStruct)
    
    #zero autos
    for i in range(nantpol): accSB[i,i]=0
    
    if opts.mode.startswith('real'): pwrSB=accSB.real
    elif opts.mode.startswith('imag'): pwrSB=accSB.imag
    elif opts.mode.startswith('lin'): pwrSB=numpy.abs(accSB)
    elif opts.mode.startswith('log'): pwrSB=numpy.log(numpy.abs(accSB))
    elif opts.mode.startswith('phs'): pwrSB=numpy.angle(accSB)
    
    pylab.imshow(pwrSB)
    pylab.colorbar()
    print numpy.min(pwrSB), numpy.max(pwrSB)
    pylab.show()

