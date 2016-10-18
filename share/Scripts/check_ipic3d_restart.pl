#!/usr/bin/perl 

foreach $file (glob "restart*.hdf"){
    $_ = `h5dump $file | head -1`;
    print "Error opening $file:$_\n" unless /^HDF5/;
};
