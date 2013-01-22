lss
===

LOFAR Single Station correlation scripts  
Last Updated: 22.01.13  
Contact: griffin.foster@gmail.com  

This is a collection of python scripts useful for plotting LOFAR station correlation data products and converting them to standard Measurement Sets.

Use the '-h' flag for script help

Requirements
===

* Numpy  
* PyEphem (rhodesmill.org/pyephem)  
* pyrap (code.google.com/p/pyrap)  
* Matplotlib/pylab w/ Basemap  
* healpy (https://github.com/healpy/healpy)  
* MeqTrees (1.2+) 
* makems (Part of the MeqTrees repo)  

Test Data
===

Example LBA and HBA correlation files is available at:  

* LBA (RCU 3): http://www-astro.physics.ox.ac.uk/~FosterG/example/20120513_052251_acc_512x192x192.dat 
* HBA (RCU 5): http://www-astro.physics.ox.ac.uk/~FosterG/example/20111205_162258_acc_512x192x192.dat 

Script Examples
===

Extract Subbandss from ACC files and generate Measurement Sets:

`gen_sb_ms_n.py -F AntennaField.conf -A AntennaArrays.conf -r 3 -s 103,155,200,245,300,350,400,420 -t 1 20120513_052251_acc_512x192x192.dat`

MeqTrees batch calibration:

`mt_calico-wsrt-tens-lofar.py -s 400 -p tdl_lba_s400.conf -S calico-wsrt-tens-lofar -P 20120612_04 LBH_HPF10MHZ`

Plot *G* Gains:

`plot_g_diag.py LBH_HPF10MHZ/s400/*.ms`

Generate *G* gains file:

`gen_g_diag.py LBH_HPF10MHZ/s400/*.ms`

Apply *G* gains:

`apply_g_diag.py -g gains_s400.ms_20120612.pkl LBH_HPF10MHZ/s400/*.ms`

Generate dirty images:

`batch_dirty_image.py LBH_HPF10MHZ/s400/*.ms`

Generate healpix maps:

`gen_healpix_map.py *.fits`

