#!/usr/bin/perl
#^CFG COPYRIGHT UM

# Change the byte order of every 4 bytes

# Usage:  Fix4Endian.pl InFile > Outfile
#         Fix4Endian.pl < InFile > Outfile

# No end of line record
undef $/;

# Read the whole file into $_
$_=<>;

# Reverse every 4 bytes
for($i=0; $i<length(); $i+=4){
    substr($_,$i,4)=reverse substr($_,$i,4);
}

print;
