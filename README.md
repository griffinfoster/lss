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

Currently the package doesn't have a installation module, to use the software add the src directory to the $PYTHONPATH, e.g. add export PYTHONPATH=$PYTHONPATH:/home/[user]/sw/lss/src to your .bashrc file

Test Data
===

Example LBA and HBA correlation files is available at:  

* LBA (RCU 3): http://www-astro.physics.ox.ac.uk/~FosterG/example/20120513_052251_acc_512x192x192.dat 
* HBA (RCU 5): http://www-astro.physics.ox.ac.uk/~FosterG/example/20111205_162258_acc_512x192x192.dat 

Script Examples
===

Extract Subbandss from ACC files and generate Measurement Sets:

`gen_sb_ms_n.py -F AntennaField.conf -A AntennaArrays.conf -r 3 -s 103,155,200,245,300,350,400,420 -t 1 20120513_052251_acc_512x192x192.dat`

