#!/usr/bin/perl
#^CFG COPYRIGHT UM
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Binary Data Manipulation}
#
#!ROUTINE: Fix4Endian.pl - change the byte order of every 4 bytes
#!DESCRIPTION:
# Change the byte order in a binary file which consists of single precision
# reals and 4 byte integers only.
#\begin{verbatim}
# Usage:  Fix4Endian.pl InFile > Outfile
#         Fix4Endian.pl < InFile > Outfile
#\end{verbatim}
#!REVISION HISTORY:
# 07/13/2001 G. Toth - initial version
#EOP
#BOC
# No end of line record
undef $/;

# Read the whole file into $_
$_=<>;

# Reverse every 4 bytes
for($i=0; $i<length(); $i+=4){
    substr($_,$i,4)=reverse substr($_,$i,4);
}

print;
#EOC
