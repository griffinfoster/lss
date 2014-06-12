#!/usr/bin/env python
"""
Generate an OSKAR 2 telescope model from a list of LOFAR stations
"""

#TODO: Missing RS210-AntennaArrays.conf
#TODO: Missing RS305-AntennaArrays.conf
#TODO: Missing RS310-AntennaArrays.conf
#TODO: Missing RS407-AntennaArrays.conf
#TODO: Missing RS409-AntennaArrays.conf
#TODO: Missing FI609-AntennaArrays.conf
#TODO: Generate antpos without AntennaArrays files?

#Current available LOFAR Stations
terps=['CS002','CS003','CS004','CS005','CS006','CS007'] #superterp stations
cs=['CS001','CS011','CS013','CS017','CS021','CS024','CS026','CS028','CS030','CS031','CS032','CS101','CS103','CS201','CS301','CS302','CS401','CS501'] #core stations
rs=['RS106','RS205','RS208','RS306','RS307','RS406','RS503','RS508','RS509'] #remote stations
ins=['DE601','DE602','DE603','DE604','DE605','FR606','SE607','UK608'] #international stations

#Hardcoded data directory
lofarData='/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/'
hbaTile='/home/griffin/AABeamGen/data/LOFAR/HBAtile.txt'

import os,sys,shutil
import time
import numpy as n
import lofarConfig

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] ACC_FILES')
    o.set_description(__doc__)
    o.add_option('-s','--stations',dest='stations',default='core',
        help='List of stations to use in telescope model, options: comma separated list, superterp, core (superterp + core stations), remote (superterp + core stations + remote stations), all . default: ccore')
    o.add_option('-m','--mode',dest='mode',default='lba',
        help='LOFAR station type: LBA or HBA , default: LBA')
    o.add_option('--single',dest='single',action='store_true',
        help='Create an OSKAR directory structure where all elements are stations, useful for single station or superterp simulation')
    o.add_option('-n','--name',dest='name',default=None,
        help='Telescope name, default:telescope+timestamp')
    o.add_option('--csmode',dest='csmode',default=None,
        help='Add an extra layout option for core stations: LBA_INNER, LBA_OUTER, LBA_X, LBA_Y, LBA_SPARSE0, LBA_SPARSE1, HBA_0, HBA_1, HBA0, HBA1')
    o.add_option('--rsmode',dest='rsmode',default=None,
        help='Add an extra layout option for remote stations: LBA_INNER, LBA_OUTER, LBA_X, LBA_Y, LBA_SPARSE0, LBA_SPARSE1')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.stations.lower().startswith('super'): stations=terps
    elif opts.stations.lower().startswith('core'): stations=terps+cs
    elif opts.stations.lower().startswith('remote'): stations=terps+cs+rs
    elif opts.stations.lower().startswith('all'): stations=terps+cs+rs+ins
    else: stations=opts.stations.upper().split(',')

    mode=opts.mode.upper()
    if opts.csmode is None: csmode=mode
    else: csmode=opts.csmode.upper()
    if opts.rsmode is None: rsmode=mode
    else: rsmode=opts.rsmode.upper()

    if not os.path.exists(lofarData):
        print 'ERROR: hardcoded data directory path does not exist, edit this script'
        exit(1)

    if opts.name is None:
        rootDir='telescope%i'%int(time.time())
    else: rootDir=opts.name

    logStr=''
    logStr+='Generating an OSKAR telescope\n'
    logStr+='Stations: %s\n'%stations
    logStr+='Mode: '+mode+'\n'
    logStr+='Core Station Mode: '+csmode+'\n'
    logStr+='Remote Stayion Mode: '+rsmode+'\n'
    logStr+='Directory: '+rootDir+'\n'
    if opts.single: logStr+='Single Station Mode: True\n'
    print logStr

    if not os.path.exists(rootDir): os.makedirs(rootDir)
    stationLocs=[]
    if opts.single:
        #use the first station in the list as the reference station
        refStation=lofarConfig.lofarStation(stations[0],lofarData+stations[0]+'-AntennaField.conf',lofarData+stations[0]+'-AntennaArrays.conf')
        for sid,s in enumerate(stations):
            ls=lofarConfig.lofarStation(s,lofarData+s+'-AntennaField.conf',lofarData+s+'-AntennaArrays.conf')
            offset=lofarConfig.relativeStationOffset(refStation,ls,mode=mode)
            if ls.stype=='international':
                antpos=ls.antArrays.antpos[mode]
                modeStr=mode
            elif ls.stype=='remote':
                antpos=ls.antArrays.antpos[rsmode]
                modeStr=rsmode
            elif ls.stype=='core':
                antpos=ls.antArrays.antpos[csmode]
                modeStr=csmode
            stationLocs.append(antpos[:,0,:]+offset)
        stationLocs=n.array(stationLocs)
        stationLocs=stationLocs.reshape((stationLocs.shape[0]*stationLocs.shape[1],stationLocs.shape[2]))
        n.savetxt(rootDir+'/layout.txt', stationLocs, fmt='%f', header='single station mode, %s'%stations)

        stationDir=rootDir+'/station'
        os.makedirs(stationDir)
        n.savetxt(stationDir+'/layout.txt', n.array([[0.,0.,0.]]), fmt='%f')
        if mode=='HBA':
            #generate tile structure
            tileDir=stationDir+'/tile'
            os.makedirs(tileDir)
            shutil.copy(hbaTile, tileDir+'/layout.txt')
            
    else:
        for sid,s in enumerate(stations):
            stationDir=rootDir+'/station%03d'%sid
            os.makedirs(stationDir)
            ls=lofarConfig.lofarStation(s,lofarData+s+'-AntennaField.conf',lofarData+s+'-AntennaArrays.conf')
            stationLocs.append(ls.antField.location[mode])
            if ls.stype=='international':
                antpos=ls.antArrays.antpos[mode]
                modeStr=mode
            elif ls.stype=='remote':
                antpos=ls.antArrays.antpos[rsmode]
                modeStr=rsmode
            elif ls.stype=='core':
                antpos=ls.antArrays.antpos[csmode]
                modeStr=csmode
            n.savetxt(stationDir+'/layout.txt', antpos[:,0,:], fmt='%f', header='%s %s'%(ls.name,modeStr))

            if mode=='HBA':
                #generate tile structure
                tileDir=stationDir+'/tile'
                os.makedirs(tileDir)
                shutil.copy(hbaTile, tileDir+'/layout.txt')
            
        stationLocs=n.array(stationLocs)
        #write layout_etrf.txt
        n.savetxt(rootDir+'/layout_ecef.txt', stationLocs, fmt='%f', header='%s'%stations)
    
    #write lofarGenTelescope.log
    logStr+=time.asctime()
    logFile=open(rootDir+'/lofarGenTelescope.log','w')
    logFile.write(logStr)
    logFile.close()

