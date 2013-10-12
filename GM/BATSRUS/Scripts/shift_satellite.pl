#!/usr/bin/perl -p -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
# Shift the satellite trajectory by a constant vector (dx, dy, dz)
# Usage: shift_satellite.pl -dx=1.3 -dy=-1.2 -dz=0.1 ORIGINAL.sat > SHIFTED.sat
s/([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)$/sprintf("%10.5f %10.5f %10.5f",$1+$dx,$2+$dy,$3+$dz)/e;


