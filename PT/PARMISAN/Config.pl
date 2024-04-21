#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing                                                        
$^I = "";

# Add local directory to search                                                 
push @INC, ".";

use strict;

our $Component = "PT";
our $Code = "PARMISAN";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments= @ARGV;


my $config     = "share/Scripts/Config.pl";
# get util and share
my $GITCLONE = "git clone"; my $GITDIR = "git\@github.com:SWMFsoftware/";

if (-f $config or -f "../../$config"){
}else{
    `$GITCLONE $GITDIR/share.git; $GITCLONE $GITDIR/util.git`;
}



if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# Variables inherited from share/Scripts/Config.pl
our %Remaining; # Unprocessed arguments
our $ERROR;
our $WARNING;
our $Help;
our $Verbose;
our $Show;
our $ShowGridSize;
our $NewGridSize;
our $NewGhostCell;


&print_help if $Help;

my $Src = 'src';

# Grid size variables
my $NameSizeFile = "$Src/ModSize.f90";
my $GridSize;
my $nX;
my $nP;

# Read previous grid size
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};
    if(/^-g$/)                {$ShowGridSize=1;                next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


# Set new grid size
&set_grid_size if ($NewGridSize and $NewGridSize ne $GridSize);

# Show current settings
my $Settings = &current_settings; print $Settings if $Show;

# Show grid
print "$GridSize\n" if $ShowGridSize;

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/;
	$nX   =$1 if /\bnVertexMax\s*=[^0-9]*(\d+)/i;
	$nP   =$1 if /\bnParticle\s*=[^0-9]*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read nVertexMax from $NameSizeFile\n" 
	unless length($nP);                         

    $GridSize = "$nP";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

      if($GridSize =~ /^[1-9]\d*,[1-9]\d*$/){
	  ($nX,$nP) = split(',', $GridSize);
    }elsif($GridSize =~ /^[1-9]\d*$/){
	$nX = $GridSize;
    }elsif($GridSize){
	die "$ERROR -g=$GridSize must be one to three integers\n";
    }
    # Check the grid size (to be set)
    die "$ERROR nParticleMax=$nX must be positive\n" if $nX<=0;

    print "Writing new grid size $GridSize into ".
	"$NameSizeFile ...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nVertexMax\s*=[^0-9]*)(\d+)/$1$nX/i;
	s/\b(nParticle\s*=[^0-9]*)(\d+)/$1$nP/i if length($nP)>0;
	print;
    }
}

#############################################################################

sub print_help{

    print "
Additional options for PARMISAN/Config.pl:

-g=nX,nP
                Set grid size. 
                nX is maximum number of particles per field line,
\n";
    exit 0;
}

#############################################################################

sub current_settings{
    $Settings  = "Number of nodes per line   : nNode=$nX\n";
    $Settings .= "Number of particles per line =$nP\n";
	
}

