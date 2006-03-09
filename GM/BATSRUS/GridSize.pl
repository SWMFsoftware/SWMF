#!/usr/bin/perl -s -i

# Put switches into properly named variables
my $Help         = $h or $H or $help or $Help;
my $nCells       = $c;
my $MaxBlock     = $b;
my $GridSize     = $g;
my $ShowHuman    = $s;
my $ShowCompact  = not($c or $b or $i or $g or $s);

use strict;

# Print help if required
&print_help if $Help;

# Set error string for error messages
my $ERROR = 'GridSize_ERROR:';
my $WARNING = 'GridSize_WARNING:';

# Set the name of the file that contains the grid size
my $NameFile="src/ModSize.f90";

# Check parameters
die "$ERROR Do not use -g together with -b or -c\n"
    if $GridSize and ($MaxBlock or $nCells);

# Set grid size variables based on the arguments
my ($nI, $nJ, $nK);

my $MaxImplBlock;

if($GridSize=~/^\d+,\d+,\d+,\d+
   ,\d+  #^CFG IF IMPLICIT
   $/x){
    ($nI,$nJ,$nK,$MaxBlock,$MaxImplBlock)= split(',', $GridSize);
}elsif($GridSize){
    die "$ERROR -g=$GridSize should be ".
	"5".  #^CFG IF IMPLICIT
	#"4". #^CFG IF NOT IMPLICIT
	" integers separated with commas\n";
}

if($nCells=~/^\d+$/){
    $nI=$nCells; $nJ=$nCells; $nK=$nCells;
}elsif($nCells=~/^\d+,\d+,\d+/){
    ($nI,$nJ,$nK)= split(',',$nCells);
}elsif($nCells){
    die "$ERROR -c=$nCells should be a single integer ".
	"or 3 integers separated with commas\n";
}

# Read previous grid size and set IsNewGrid to true if there is change
my $IsNewGrid;
&read_grid_size;

# Check the grid size (to be set)
die "$ERROR nI=$nI must be 4 or more\n" if $nI < 2;
die "$ERROR nJ=$nJ must be 4 or more\n" if $nJ < 2;
die "$ERROR nK=$nK must be 4 or more\n" if $nK < 2;

die "$ERROR nI=$nI must be an even integer\n" if $nI!=int($nI) or $nI%2!=0; 
die "$ERROR nJ=$nJ must be an even integer\n" if $nJ!=int($nJ) or $nJ%2!=0; 
die "$ERROR nK=$nK must be an even integer\n" if $nK!=int($nK) or $nK%2!=0; 

warn "$WARNING nI=$nI nJ=$nJ nK=$nK does not allow AMR\n" if $nI == 2 or $nJ == 2 or $nK==2;

die "$ERROR MaxBlock=$MaxBlock must be a positive integer\n" 
    if $MaxBlock<1 or $MaxBlock != int($MaxBlock);

#^CFG IF IMPLICIT BEGIN
die "$ERROR MaxImplBlock=$MaxImplBlock must be a positive integer\n" 
    if $MaxImplBlock<1 or $MaxImplBlock != int($MaxImplBlock);

die "$ERROR MaxImplBlock=$MaxImplBlock cannot exceed MaxBlock=$MaxBlock\n"
    if $MaxImplBlock > $MaxBlock;
#^CFG END IMPLICIT

# Write the new grid into $NameFile if necessary
if($IsNewGrid){
    &set_grid_size;
    # Check if it worked be rereading info. Now IsNewGrid should be false.
    &read_grid_size;
    if($IsNewGrid){
	&show_grid_size;
	die "$ERROR incorrect grid information was saved";
    }
}

&show_grid_size if $ShowHuman;

print "GridSize.pl -g=$nI,$nJ,$nK,$MaxBlock",
    ",$MaxImplBlock"                           #^CFG IF IMPLICIT
    ,"\n" if $ShowCompact;

exit 0;

#############################################################################

sub read_grid_size{

    # Local variables to be read
    my ($nI_, $nJ_, $nK_, $MaxBlock_, $MaxImplBlock_);

    # Read size of the grid from $NameFile
    open(MODSIZE,$NameFile) or die "$ERROR could not open $NameFile\n";
    while(<MODSIZE>){
	next if /^\s*!/; # skip commented out lines
        $nI_=$1           if /\bnI\s*=\s*(\d+)/i;
	$nJ_=$1           if /\bnJ\s*=\s*(\d+)/i;
	$nK_=$1           if /\bnK\s*=\s*(\d+)/i;
        $MaxBlock_=$1     if /\bnBLK\s*=\s*(\d+)/i;
	$MaxImplBlock_=$1 if /\bMaxImplBLK\s*=[^0-9]*(\d+)/i;
    }
    close MODSIZE;

    die "$ERROR could not read nI from $NameFile\n" unless length($nI_);
    die "$ERROR could not read nJ from $NameFile\n" unless length($nJ_);
    die "$ERROR could not read nK from $NameFile\n" unless length($nK_);
    die "$ERROR could not read MaxBlock from $NameFile\n" 
	unless length($MaxBlock_);
    die "$ERROR could not read MaxImplBlock from $NameFile\n" #^CFG IF IMPLICIT
	unless length($MaxImplBlock_);                        #^CFG IF IMPLICIT

    $nI           = $nI_           unless $nI;
    $nJ           = $nJ_           unless $nJ;
    $nK           = $nK_           unless $nK;
    $MaxBlock     = $MaxBlock_     unless $MaxBlock;
    $MaxImplBlock = $MaxImplBlock_ unless $MaxImplBlock;

    $IsNewGrid = "$nI,$nJ,$nK,$MaxBlock,$MaxImplBlock" ne 
	"$nI_,$nJ_,$nK_,$MaxBlock_,$MaxImplBlock_";
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size into $NameFile...\n";

    @ARGV = ($NameFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	s/\b(nBLK\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	s/\b(MaxImplBLK\s*=[^0-9]*)(\d+)/$1$MaxImplBlock/i;
	print;
    }

}

##############################################################################

sub show_grid_size{

    print "Number of cells in a block        : nI=$nI, nJ=$nJ, nK=$nK\n";
    print "Max. number of blocks/PE          : MaxBlock    =$MaxBlock\n";
    #^CFG IF IMPLICIT BEGIN
    print "Max. number of implicit blocks/PE : MaxImplBlock=$MaxImplBlock\n";
    #^CFG END IMPLICIT

}

#############################################################################

sub print_help{

    print "
Purpose:

    Get and set grid size information from and into src/ModSize.f90
    and check if settings are correct.
    If called without any switch the current settings are shown in
    a compact format which is suitable for processing by other codes.

Usage:

    GridSize.pl [-h] [-c=n|nI,nJ,nK] [-b=MaxBlock] 
                [-g=nI,nJ,nK,MaxBlock",
    ",MaxImplBlock", #^CFG IF IMPLICIT
    "] [-s]

  -h              print help message and stop

  -c=n            set the block size to be n*n*n

  -c=nI,nJ,nK     set the block size to be nI*nJ*nK

  -b=MaxBlock     set the maximum number of blocks per processor to MaxBlock

  -g=nI,nJ,nK...  set the block size and the number of blocks together

  -s              show current settings in human readable form

Examples:

  Show current settings in a compact form:

GridSize.pl

  Set block size to 6x6x6 and the number of blocks to 500 
  and show current settings in a human readable form:

GridSize.pl -c=6 -b=100 -s

  Set blocks size to 8x8x8, number of blocks to 400",
" and implicit blocks to 100", #^CFG IF IMPLICIT
":

GridSize.pl -g=8,8,8,400",
",100", #^CFG IF IMPLICIT
"\n\n";
    exit 0;
}
