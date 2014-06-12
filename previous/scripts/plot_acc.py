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
        help='Select which subband to select out for the MS, default:0')
    o.add_option('-m', '--mode', dest='mode', default='lin',
        help='Plotting mode: linear,log,real,imag,phs, default:linear')
    opts, args = o.parse_args(sys.argv[1:])
acc_files=args

for fid,fn in enumerate(acc_files):
    print '\t' + fn + '(%i of %i)'%(fid+1,len(acc_files))
    fnBase=fn.split('/')[-1]
    fnArr=fnBase.split('_')
    accSize=fnArr[3].split('x')
    nsb=int(accSize[0])
    nantpol=int(accSize[1])
    accComp=numpy.fromfile(fn,dtype='complex')
    print accComp.shape
    accComp=accComp.reshape(nsb,nantpol,nantpol)
    accSB=accComp[opts.subband,:,:]
    accSB=accSB.reshape(nantpol,nantpol,1)
    print accSB.shape,accComp.shape
    
    fTimeStruct=time.struct_time((int(fnArr[0][0:4]),int(fnArr[0][4:6]),int(fnArr[0][6:8]),int(fnArr[1][0:2]),int(fnArr[1][2:4]),int(fnArr[1][4:6]),0,0,0))
    fTime=time.mktime(fTimeStruct)
    
    #zero autos
    for i in range(nantpol): accSB[i,i,:]=0
    
    if opts.mode.startswith('real'): pwrSB=accSB[:,:,0].real
    elif opts.mode.startswith('imag'): pwrSB=accSB[:,:,0].imag
    elif opts.mode.startswith('lin'): pwrSB=numpy.abs(accSB[:,:,0])
    elif opts.mode.startswith('log'): pwrSB=numpy.log(numpy.abs(accSB[:,:,0]))
    elif opts.mode.startswith('phs'): pwrSB=numpy.angle(accSB[:,:,0])
    
    pylab.imshow(pwrSB)
    pylab.colorbar()
    print numpy.min(pwrSB), numpy.max(pwrSB)
    pylab.show()

