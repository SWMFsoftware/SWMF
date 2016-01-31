#!/usr/bin/perl -i
use strict;
use warnings;

my %mode=(val =>"GAS",'flag' => 0);
my %OutputFileNumber=('val' => 0,'flag' => 0);
my %nDustSpec=('val' => 1,'flag' => 0);
my %nDustGroup=('val' => 1,'flag' => 0);
my %nGasSpec=('val' => 2,'flag' => 0);

#read the argument line
foreach (@ARGV) { 
   if ((/^-h$/)||(/^-help/)) { #print help message
     print "The script sets the post processing code of the CG model\n";
     print "-mode=[gas,dust]\tsets the gas or dust\n";
     print "-out=[number]\t\tthe number of the output file that will be read\n"; 
     print "-nndust=[number]\tthe number of the dust species\n"; 
     print "-ngas=[number]\t\tthe number of the gas species\n";
     print "-ngroup=[number]\tthe number of the dust size groups\n";
     
     exit;
   }
  elsif (/^-mode=(.*)$/i)  {$mode{val}=uc($1); $mode{flag}=1; next;}
  elsif (/^-out=(.*)$/i)   {$OutputFileNumber{val}=$1; $OutputFileNumber{'flag'}=1; next;}
  elsif (/^-ndust=(.*)$/i) {$nDustSpec{val}=$1; $nDustSpec{flag}=1; next;}
  elsif (/^-ngas=(.*)$/i) {$nGasSpec{val}=$1; $nGasSpec{flag}=1; next;}
  elsif (/^-ngroup=(.*)$/i) {$nDustGroup{val}=$1; $nDustGroup{flag}=1; next;}
  else {
    print "option $_ is unknown\n";
    exit;
  }
}

#process the configuration settings
if ($mode{flag} == 1) {
  `echo "\n#undef _MODE_\n#define _MODE_ _"$mode{val}"_MODE_" >> config.h`;
}

if ($nDustSpec{flag}==1) {
  `echo "#undef _DUST_CASE_\n#define _DUST_CASE_ _DUST_CASE__"$nDustSpec{val}"SPEC_"$nDustGroup{val}"GROUP_" >> config.h`;
  `echo "#undef _DUST_SPEC_NUMBER_\n#define _DUST_SPEC_NUMBER_\ "$nDustSpec{val}"\n" >> config.h`;
}

if ($nDustGroup{flag}==1) {
  `echo "#undef _DUST_CASE_\n#define _DUST_CASE_ _DUST_CASE__"$nDustSpec{val}"SPEC_"$nDustGroup{val}"GROUP_" >> config.h`;
  `echo "#undef _DUST_GROUP_NUMBER_\n#define _DUST_GROUP_NUMBER_\ "$nDustGroup{val}"\n" >> config.h`;
}

if ($nGasSpec{flag}==1) {
  `echo "#undef _GAS_SPEC_NUMBER_\n#define _GAS_SPEC_NUMBER_\ "$nGasSpec{val}"\n" >> config.h`;
}

if ($OutputFileNumber{flag}==1) {
  `echo "#undef _OUTPUT_FILE_NUMBER_\n#define _OUTPUT_FILE_NUMBER_ "$OutputFileNumber{val}  >> config.h`;
}

if (! -e "config.h") { 
  `echo "\n" >> config.h`;
}
