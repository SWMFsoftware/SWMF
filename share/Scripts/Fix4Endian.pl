#!/usr/bin/perl -s
#^CFG COPYRIGHT UM
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Binary Data Manipulation}
#
#!ROUTINE: Fix4Endian.pl - change the byte order of every 4 bytes
#!DESCRIPTION:
# Change the byte order in a binary file which consists of single precision
# reals and 4 byte integers only. The -except=M:N flag allows to exclude
# the bytes between the M-th to N-th positions (usually a character string).
#\begin{verbatim}
# Usage:  
#         Fix4Endian.pl [-except=M:N] < InFile > Outfile
#         perl -is Fix4Endian.pl [-except=M:N] File1 File2 ...
#\end{verbatim}
#!REVISION HISTORY:
# 07/13/2001 G. Toth - initial version
# 08/18/2004 G. Toth - added -except=M:N option for RCM restart files.
#EOP
#BOC
# No end of line record
undef $/;

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

if($max > length()){
    $max = length();
    warn "Maximum byte number in -except=$except exceeds file length\n";
}

# Reverse every 4 bytes
for($i=0; $i<length(); $i+=4){
    $i = $max if $i+1 >= $min and $i+1 <= $max;
    warn "i=$i\n";
    substr($_,$i,4)=reverse substr($_,$i,4);
}

print;
#EOC
