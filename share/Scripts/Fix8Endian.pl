#!/usr/bin/perl
#^CFG COPYRIGHT UM

# Change the endianness (byte order) in a Fortran file,
# which contains 8 byte reals and integers only

# Usage:  Fix8Endian.pl InFile > Outfile
#         Fix8Endian.pl < InFile > Outfile

# No end of line record
undef $/;

# Read the whole file into $_
$_=<>;

# initialize pointer for the string
$i=0;
while ( $i < length() ){
    # Get length of record
    $len = unpack('L',substr($_,$i,4));

    # Check if length is reasonable
    die "At position $i record length $len is too large?!\n" 
	if $i+$len > length();

    # Reverse leading 4 byte length marker
    $lenfixed = reverse(substr($_,$i,4));
    substr($_,$i,4)=$lenfixed;

    # Reverse 8 byte reals/integers
    for($j=$i+4; $j<$i+$len; $j+=8){
	substr($_,$j,8)=reverse(substr($_,$j,8));
    }

    # Reverse trailing 4 byte length marker
    substr($_,$j,4)=$lenfixed;

    # Step to the next record
    $i = $j + 4;
}

print;
