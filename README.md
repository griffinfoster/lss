lss
===

LOFAR Single Station Software 

Last Updated: 12.06.14 

Contact: griffin.foster@gmail.com  

This is a collection of python scripts useful for manipulating LOFAR station correlation data products and converting them to standard Measurement Sets.

Use the '-h' flag for script help

Requirements
===

* Numpy  
* PyEphem (rhodesmill.org/pyephem)  
* pyrap (code.google.com/p/pyrap)  
* Matplotlib/pylab w/ Basemap  
* makems (Part of the MeqTrees repo)  

Installation
===

Currently the package doesn't have a installation module, to use the software add this to your .bashrc file:

`export PYTHONPATH=$PYTHONPATH:/home/[user]/sw/lss/src`

`export LSSDATA=/home/[user]/sw/lss/data`

Test Data
===

Example LBA and HBA correlation files is available at:  

* LBA (RCU 3): http://www-astro.physics.ox.ac.uk/~FosterG/example/20120513_052251_acc_512x192x192.dat 
* HBA (RCU 5): http://www-astro.physics.ox.ac.uk/~FosterG/example/20111205_162258_acc_512x192x192.dat 

Script Examples
===

Extract Subbandss from ACC files and generate Measurement Sets:

`gen_sb_ms_n.py -F AntennaField.conf -A AntennaArrays.conf -r 3 -s 103,155,200,245,300,350,400,420 -t 1 20120513_052251_acc_512x192x192.dat`

To Do
===

Add the following scripts:

* derive polarization calibration file
* derive intensity calibration file
* plot calibration file
* apply calibration files
* plot XYZ,  UVW
* plot PSF
* convert XST/ACC -> MS
* plot SST file
* convert DFT image -> FITS
* convert FFT image -> FITS
* plot all sky FITS file
* convert XST/ACC files -> HDF5 container
* convert SST files -> HDF5 container

FFT Plotter:

* add convultion functions
* correct u,v -> l,m scaling for hemisphere imaging
* add sqrt(1-l^2-m^2) weighting

AABeamGen:

* Generate antpos without AntennaArrays files
* test with OSKAR

Installation and Setup:

* remove hard-coded paths
* make a proper python module

