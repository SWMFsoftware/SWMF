#!/usr/bin/perl -pi~
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Purpose:
# Replace old #NONCONSERVATIVE command with new one.
#
# For example
#
# #NONCONSERVATIVE
# -1.0
#
# is replaced with
#
# #NONCONSERVATIVE
# F
#
# and
#
# #NONCONSERVATIVE
# 6.0
#
# is replaced with
#
# #NONCONSERVATIVE
# T
#
# #CONSERVATIVECRITERIA
# 1
# r
# 6.0
#
# Usage:
#
#    FixNonconservative.pl PARAM1.in PARAM2.in ...
#
# The original files will be saved into PARAM1.in~ PARAM2.in~ ...

if($read){
    if(not /^\s*(T|F|.true.|.false.)/i){
	if(/^\s*\-/){
	    $_="F\n";
	}else{
	    print "T\t\t\tUseNonConervative

#CONSERVATIVECRITERIA
1\t\t\tnConservCrit
r\t\t\tTypeConservCrit
";
	}
    }
    $read=0;
}

$read=1 if /^#NONCONSERVATIVE/;
