#!/usr/bin/perl -pi -s

# Usage:
#
#    fixmpi.pl [-mpi=ModMpiOrig] *.f90 *.F90
#
#  The default mpi module is ModMpi
#
# Replace 
# 
#   implicit none
#   ! at most 2 lines here
#   include '../Common/mpif90.h'
#
# with
#
#   use ModMpi
#   implicit none
#   ! atmost 2 lines here
#
#  and the remaining 
#
#   include '../Common/mpif90.h' 
#
# strings with
#
#   use ModMpi

BEGIN{$mpi = "ModMpi" unless $mpi}

if(/^ *implicit +none/i){
    $implnone=$_;
    for $i (1..3){
	$_=<>;
	if(/^( *)include +['"].*mpif90.h["']/i){
	    $_="$1use $mpi\n$implnone";
	    last;
	}else{
	    $implnone.=$_;
	    $_=$implnone;
	}
    }
}

# Replace the remaining strings
s/include +['"].*mpif90.h["']/use $mpi/gi;
