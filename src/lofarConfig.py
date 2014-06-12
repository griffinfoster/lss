"""
Functions and classes to read and parse LOFAR station configuration files
"""

#TODO: use rotation matrix in antennaField?
#TODO: drop antennaArray

import numpy as n

class lofarStation():
    def __init__(self,name,affn,aafn):
        self.name=name
        self.antField=antennaField(name,affn)
        self.antArrays=antennaArrays(name,aafn)
        if name.lower().startswith('cs'): self.stype='core'
        elif name.lower().startswith('rs'): self.stype='remote'
        else: self.stype='international'

class antennaField():
    def __init__(self, name, fn):
        self.name=name
        self.rotMatrix={}
        self.normVec={}
        self.antpos={}
        self.location={}
        fh=open(fn)
        lines=[]
        lines=fh.read().split('\n')
        fh.close()
        cleanStr=''
        for l in lines:
            if l=='' or l.startswith('#'): continue
            cleanStr+=" ".join(l.split())+" "
        cleanStr=" ".join(cleanStr.split())
        lastMode=None
        for l in cleanStr.split(']'):
            if l=='':continue
            infoStr,dataStr=l.split('[')
            infoStr=infoStr.strip()
            if infoStr.lower().startswith('normal'):
                name,mode,length=infoStr.split(' ')
                self.normVec[mode]=n.array(map(float,dataStr.strip().split(' ')))
            elif infoStr.lower().startswith('rotation'):
                name,mode,l0,fill,l1=infoStr.split(' ')
                self.rotMatrix[mode]=n.array(map(float,dataStr.strip().split(' '))).reshape((int(l0),int(l1)))
            elif infoStr.lower().startswith('lba') or infoStr.lower().startswith('hba'):
                mode,length=infoStr.split(' ')
                self.location[mode]=n.array(map(float,dataStr.strip().split(' ')))
                lastMode=mode
            else: #antenna positions
                l0,f0,l1,f1,l2=infoStr.split(' ')
                self.antpos[lastMode]=n.array(map(float,dataStr.strip().split(' '))).reshape((int(l0),int(l1),int(l2)))
        #print self.rotMatrix.keys(), self.normVec.keys(),self.antpos.keys(),self.location.keys()

class antennaArrays():
    def __init__(self, name, fn):
        self.name=name
        self.antpos={}
        self.location={}
        fh=open(fn)
        lines=[]
        lines=fh.read().split('\n')
        fh.close()
        cleanStr=''
        for l in lines:
            if l=='' or l.startswith('#'): continue
            cleanStr+=" ".join(l.split())+" "
        cleanStr=" ".join(cleanStr.split())
        lastMode=None
        for l in cleanStr.split(']'):
            if l=='':continue
            infoStr,dataStr=l.split('[')
            infoStr=infoStr.strip()
            if infoStr.lower().startswith('lba') or infoStr.lower().startswith('hba'):
                mode,length=infoStr.split(' ')
                self.location[mode]=n.array(map(float,dataStr.strip().split(' ')))
                lastMode=mode
            else: #antenna positions
                l0,f0,l1,f1,l2=infoStr.split(' ')
                self.antpos[lastMode]=n.array(map(float,dataStr.strip().split(' '))).reshape((int(l0),int(l1),int(l2)))
        #print self.antpos.keys(),self.location.keys()

def rotationMatrix(alpha,beta,gamma):
    """Generic rotation matrix to apply to an XYZ co-ordinate"""
    return n.matrix([[n.cos(beta)*n.cos(gamma), n.cos(gamma)*n.sin(alpha)*n.sin(beta)-n.cos(alpha)*n.sin(gamma), n.cos(alpha)*n.cos(gamma)*n.sin(beta)+n.sin(alpha)*n.sin(gamma)],
                     [n.cos(beta)*n.sin(gamma), n.cos(alpha)*n.cos(gamma)+n.sin(alpha)*n.sin(beta)*n.sin(gamma), -1.*n.cos(gamma)*n.sin(alpha)+n.cos(alpha)*n.sin(beta)*n.sin(gamma)],
                     [-1.*n.sin(beta),          n.cos(beta)*n.sin(alpha),                                        n.cos(alpha)*n.cos(beta)]])

def rotMatrixfromXYZ(station,mode='LBA'):
    """Return a rotation matrix which will rotate a station to (0,0,1)"""
    loc=station.antField.location[mode]
    longRotMat=rotationMatrix(0.,0.,-1.*n.arctan(loc[1]/loc[0]))
    loc0=n.dot(longRotMat,loc)
    latRotMat=rotationMatrix(0.,n.arctan(loc0[0,2]/loc0[0,0]),0.)
    return n.dot(latRotMat,longRotMat)

def applyRotMatrix(station,rm,mode='LBA'):
    """Apply a rotation matrix to a station location, returns new location after apply rotation matrix"""
    loc=station.antField.location[mode]
    return n.dot(rm,loc)

def relativeStationOffset(s0,s1,mode='LBA'):
    """Return the relative offset of station s1 from station s0"""
    rotMat=rotMatrixfromXYZ(s0,'LBA')
    s0loc=applyRotMatrix(s0,rotMat,mode)
    s1loc=applyRotMatrix(s1,rotMat,mode)
    return n.array(s1loc-s0loc)[0][::-1]

if __name__ == '__main__':
    print 'Running test cases'

    #antfield=antennaField('CS013','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS013-AntennaField.conf')
    #antfield=antennaField('RS208','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/RS208-AntennaField.conf')
    #antfield=antennaField('UK608','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/UK608-AntennaField.conf')

    #antArrys=antennaArrays('CS013','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS013-AntennaArrays.conf')
    #antArrys=antennaArrays('RS208','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/RS208-AntennaArrays.conf')
    #antArrys=antennaArrays('UK608','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/UK608-AntennaArrays.conf')
    CS013=lofarStation('CS013','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS013-AntennaField.conf','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS013-AntennaArrays.conf')
    print CS013.name
    RS208=lofarStation('RS208','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/RS208-AntennaField.conf','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/RS208-AntennaArrays.conf')
    print RS208.name
    UK608=lofarStation('UK608','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/UK608-AntennaField.conf','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/UK608-AntennaArrays.conf')
    print UK608.name

    CS002=lofarStation('CS002','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS002-AntennaField.conf','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS002-AntennaArrays.conf')
    CS003=lofarStation('CS003','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS003-AntennaField.conf','/home/griffin/AABeamGen/data/LOFAR/StaticMetaData/CS003-AntennaArrays.conf')
    rotMat=rotMatrixfromXYZ(CS002,'LBA')
    #print applyRotMatrix(CS002,rotMat,'LBA')
    #print applyRotMatrix(CS003,rotMat,'LBA')
    #print applyRotMatrix(CS013,rotMat,'LBA')
    #print applyRotMatrix(RS208,rotMat,'LBA')
    print relativeStationOffset(CS002,CS003)
    print relativeStationOffset(CS002,CS013)
    print relativeStationOffset(CS002,RS208)
    print relativeStationOffset(CS002,UK608)

    #print '\n'
    #aa=CS002.antField.location['LBA']
    #print aa
    #print n.sqrt(aa[0]**2.+aa[1]**2.+aa[2]**2.)
    #print CS002.antField.rotMatrix['LBA']
    #ab=n.dot(n.linalg.inv(CS002.antField.rotMatrix['LBA']),CS002.antField.location['LBA'])
    #print ab
    #print n.sqrt(ab[0]**2.+ab[1]**2.+ab[2]**2.)

    print 'Made it through without any errors.'
