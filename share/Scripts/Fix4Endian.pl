#!/usr/bin/perl -s
#^CFG COPYRIGHT UM
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Binary Data Manipulation}
#
#!ROUTINE: Fix4Endian.pl - change the byte order of every 4 bytes
#!DESCRIPTION:
# Change the byte order in a binary file which normally consists of 
# single precision reals and 4 byte integers and logicals only. 
# The -only=A:B flag allows to limit the replacement to the bytes between
# the A-th and B-th positions (positions are starting from 1).
# The -except=M:N flag allows to exclude the bytes between 
# the B-th and C-th positions (usually a character string).
#\begin{verbatim}
# Usage:  
#         Fix4Endian.pl [-only=A:B] [-except=C:D] < InFile > Outfile
#         perl -si Fix4Endian.pl [-only=A:B] [-except=C:D] File1 File2 ...
#\end{verbatim}
#!REVISION HISTORY:
# 07/13/2001 G. Toth - initial version
# 08/18/2004 G. Toth - added -except=M:N option for RCM restart files.
#EOP
#BOC
# No end of line record
undef $/;

# Extract range
if($only){
    $_ = $only;
    ($start,$finish) = /(\d+).(\d+)/;
    die "Incorrect format for -only=$only\n"
        unless 0 < $start and $start < $finish;
}else{
    $start  =  1; 
    $finish = -1;
}

# Extract excpetion range if present
if($except){
    $_ = $except;
    ($min,$max) = /(\d+).(\d+)/;
    die "Incorrect format for -except=$except\n" 
	unless 0 < $min and $min <= $max;
}else{
    $min = -1;
    $max = -1;
}

# Read the whole file into $_
$_=<>;
$length = length();

if($start > $length){
    die "Starting byte number in -only=$only exceeds file length\n";
}

if($finish < 0){
    $finish = $length;
}elsif($finish > $length){
    $finish = $length;
    warn "Final byte number in -only=$only exceeds file length\n";
}

if($min > $length){
    warn "Minimum byte number in -except=$except exceeds file length\n";
}

if($max > $length){
    $max = $length;
    warn "Maximum byte number in -except=$except exceeds file length\n";
}

# Reverse every 4 bytes
for($i = $start-1; $i < $finish; $i += 4){
    if($i+1 >= $min and $i < $max){
	$i = $max;
	last if $i >= $finish;
    }
    substr($_,$i,4)=reverse substr($_,$i,4);
}

print;
#EOC
