#!/usr/bin/env python
"""
Batch image measurement sets to produce dirty images
"""

import sys,os

imager_script='/home/griffin/meqtrees/trunk/Cattery/Meow/make_dirty_image.py'

data='CORRECTED_DATA'
mode='channel'
weight='natural'
stokes='I'
npix='512'
prefervelocity='False'
cellsize='808.59375arcsec'
spwid='0'
field='0'
padding='1.000000'
image_viewer='none'
chanmode='channel'
nchan='1'
chanstart='0'
chanstep='1'
img_nchan='1'
img_chanstart='0'
img_chanstep='1'

args0=' data=CORRECTED_DATA ms=/media/mcnuggets/lss/data/lba_20101105/apply_gains/LBH_HPF10MHZ/s400/20101105_144800_r3_s400.ms mode=channel weight=natural stokes=I npix=512 prefervelocity=False cellsize=808.59375arcsec spwid=0 field=0 padding=1.000000 image_viewer=none chanmode=channel nchan=1 chanstart=0 chanstep=1 img_nchan=1 img_chanstart=0 img_chanstep=1 fits=20101105_144800_r3_s400.ms.CORRECTED_DATA.channel.1ch.fits image=20101105_144800_r3_s400.ms.CORRECTED_DATA.channel.1ch.img remove_image'

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] MS')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    #exclude baselines using these antennas
    exclude=[37,58]
    nants=96
    select="\'(("
    firstSelect=True
    for ant1 in range(nants-1):
        for ant2 in range(nants-ant1-1):
            if (ant1 in exclude) or (ant2+ant1+1 in exclude): continue
            else:
                if firstSelect:
                    select+='(ANTENNA1==%i&&ANTENNA2==%i)'%(ant1,ant2+ant1+1)
                    firstSelect=False
                else:
                    select+='||(ANTENNA1==%i&&ANTENNA2==%i)'%(ant1,ant2+ant1+1)
    select+="))\'"

    exec_str=''
    for i in args:
        exec_str=''
        abs_i=os.path.abspath(i)
        fn_base=abs_i.split('/')[-1]
        exec_str+=' data=%s'%data
        exec_str+=' ms=%s'%abs_i
        exec_str+=' mode=%s'%mode
        exec_str+=' weight=%s'%weight
        exec_str+=' stokes=%s'%stokes
        exec_str+=' npix=%s'%npix
        exec_str+=' prefervelocity=%s'%prefervelocity
        exec_str+=' cellsize=%s'%cellsize
        exec_str+=' spwid=%s'%spwid
        exec_str+=' field=%s'%field
        exec_str+=' padding=%s'%padding
        exec_str+=' image_viewer=%s'%image_viewer
        exec_str+=' chanmode=%s'%chanmode
        exec_str+=' nchan=%s'%nchan
        exec_str+=' chanstart=%s'%chanstart
        exec_str+=' chanstep=%s'%chanstep
        exec_str+=' img_nchan=%s'%img_nchan
        exec_str+=' img_chanstart=%s'%img_chanstart
        exec_str+=' img_chanstep=%s'%img_chanstep
        exec_str+=' select=%s'%select
        exec_str+=' fits=%s'%fn_base+'.CORRECTED_DATA.channel.1ch.fits'
        exec_str+=' image=%s'%fn_base+'.CORRECTED_DATA.channel.1ch.img'
        exec_str+=' remove_image'
        #print exec_str

        os.system('python %s %s'%(imager_script,exec_str))

