#!/usr/bin/perl -i
use strict;
use warnings;

#read the argument line
foreach (@ARGV) { 
   if ((/^-h$/)||(/^-help/)) { #print help message
     print "-extraflag= to set extra flagc to do compiling in this dir\n";
     exit;
   }
  elsif (/^-extraflag=(.*)$/) {
    my $options=$1; 

    #replace FLAGC variable in makefile    
    `mv makefile makefile.bak`;
    
    open(fMakefileIn,"<makefile.bak") || die "Cannot open makefile.bak";
    open(fMakefileOut,">makefile") || die "Cannot open makefile";
    
    my $line;
    while ($line=<fMakefileIn>) {
      if ($line=~m/FLAGC/) {
        print fMakefileOut "FLAGC+ = $options\n";
      }
      else {
        print fMakefileOut "$line";
      }
    }
    
    close(fMakefileIn);
    close(fMakefileOut);
    next;
  } 
 
  else {
    print "option $_ is unknown\n";
    exit;
  }
}




