#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

use strict;

our $Component = "GM";
our $Code = "BATSRUS";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments = @ARGV;

my $config     = "share/Scripts/Config.pl";
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

# Force settings without checking
my $Force;

# Equation and user module variables
my $Src         = 'src';
my $SrcUser     = 'srcUser';
my $UserMod     = "$Src/ModUser.f90";
my $UserModSafe = "$Src/ModUser.f90.safe";
my $SrcEquation = 'srcEquation';
my $EquationMod = "$Src/ModEquation.f90";
my $EquationModSafe = "$Src/ModEquation.f90.safe";
my $Equation;
my $UserModule;

# Grid size variables
my $NameSizeFile = "$Src/ModSize.f90";
my $NameBatlFile = "srcBATL/BATL_size.f90";
my $GridSize;
my $GhostCell;
my ($nI, $nJ, $nK, $MaxBlock, $MaxImplBlock, $nG);

# additional variable information
my $nWave;
my $nWaveNew;
my $nMaterial;
my $nMaterialNew;

# For SC/BATSRUS and IH/BATSRUS src/ is created during configuration of SWMF
if(not -d $Src){exit 0};

# Read previous grid size, equation and user module
&get_settings;

foreach (@Arguments){
    if(/^-f$/)                {$Force=1;                       next};
    if(/^-e$/)                {$Equation=1;                    next};
    if(/^-u$/)                {$UserModule=1;                  next};
    if(/^-e=(.*)$/)           {$Equation=$1;                   next};
    if(/^-u=(.*)$/)           {$UserModule=$1;                 next};
    if(/^-s$/)                {$Show=1;                        next};
    if(/^-nWave=(.*)$/i)      {
	# Check the number of wave bins (to be set)
	die "$ERROR nWave=$1 must be 1 or more\n" if $1 < 1;
	$nWaveNew=$1;
	next};
    if(/^-nMaterial=(.*)$/i)  {
	# Check the number of material level indices (to be set)
	die "$ERROR nMaterial=$1 must be 1 or more\n" if $1 < 1;
	$nMaterialNew=$1;
	next};
    if(/^-ng$/i)              {print "ng=$GhostCell\n"; next};
    if(/^-ng=(.*)$/i)         {$NewGhostCell=$1; next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&set_grid_size if ($NewGridSize  and $NewGridSize ne $GridSize)
    or            ($NewGhostCell and $NewGhostCell ne $GhostCell);

# Show grid size in a compact form if requested
print "Config.pl -g=$nI,$nJ,$nK,$MaxBlock,$MaxImplBlock\n",
    if $ShowGridSize and not $Show;

# Set or list the equations
&set_equation if $Equation;

# Set additional variable information
open(FILE, $EquationMod);
while(<FILE>){
    next if /^\s*!/; # skip commented out lines
    $nWave=$1        if /\bnWave\s*=\s*(\d+)/i;
    $nMaterial=$1    if /\bnMaterial\s*=\s*(\d+)/i;
}
close FILE;
die "$ERROR nWave was not found in equation module\n" if $nWaveNew and not $nWave;
&set_nwave     if $nWaveNew and $nWaveNew ne $nWave;
die "$ERROR nMaterial was not found in equation module\n" if $nMaterialNew and not $nMaterial;
&set_nmaterial if $nMaterialNew and $nMaterialNew ne $nMaterial;

# Set or list the user modules
&set_user_module if $UserModule;

my $Settings = &current_settings; print $Settings if $Show;

# (Re)Create Makefile.RULES file(s) based on current settings
&create_makefile_rules($Settings);

exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $MaxBlock=$1     if /\bMaxBlock\s*=\s*(\d+)/i;
	$MaxImplBlock=$1 if /\bMaxImplBLK\s*=[^0-9]*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read MaxBlock from $NameSizeFile\n" 
	unless length($MaxBlock);

    die "$ERROR could not read MaxImplBlock from $NameSizeFile\n" 
	unless length($MaxImplBlock);                         

    # Make sure that BATL_size.f90 is up-to-date
    `make $NameBatlFile`;
    open(FILE, $NameBatlFile) or die "$ERROR could not open $NameBatlFile\n";
    while(<FILE>){
        next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
	$nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
	$nK=$1           if /\bnK\s*=\s*(\d+)/i;
	$GhostCell=$1    if /\bnG\s*=\s*(\d)/;
    }
    close FILE;

    die "$ERROR could not read nI from $NameBatlFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameBatlFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameBatlFile\n" unless length($nK);
    die "$ERROR could not read nG from $NameBatlFile\n" 
	unless length($GhostCell);

    $GridSize = "$nI,$nJ,$nK,$MaxBlock,$MaxImplBlock";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

    if($GridSize =~ /^[1-9]\d*,[1-9]\d*,[1-9]\d*,[1-9]\d*,[1-9]\d*$/){
	($nI,$nJ,$nK,$MaxBlock,$MaxImplBlock) = split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be 5".
	    " positive integers separated with commas\n";
    }

    $GhostCell = $NewGhostCell if $NewGhostCell;
    die "$ERROR -ng=$GhostCell must be 2,3,4 or 5\n" if $GhostCell!~/^[2345]/;

    # Check the grid size (to be set)
    die "$ERROR nK=$nK must be 1 if nJ is 1\n"         if $nJ==1 and $nK>1;
    die "$ERROR nI=$nI must be an even integer\n"      if           $nI%2!=0;
    die "$ERROR nJ=$nJ must be 1 or an even integer\n" if $nJ>1 and $nJ%2!=0;
    die "$ERROR nK=$nK must be 1 or an even integer\n" if $nK>1 and $nK%2!=0;

    die "$ERROR MaxImplBlock=$MaxImplBlock cannot exceed MaxBlock=$MaxBlock\n"
	if $MaxImplBlock > $MaxBlock;

    if(not $Force){
	die "$ERROR nI=$nI nJ=$nJ nK=$nK does not allow AMR\n" 
	    if $nI == 2 or $nJ == 2 or $nK==2;
	die "$ERROR -ng=$GhostCell should not exceed nI/2=$nI/2\n"
	    if $GhostCell > $nI/2;
	die "$ERROR -ng=$GhostCell should not exceed nJ/2=$nJ/2\n"
	    if $GhostCell > $nJ/2 and $nJ>1;
	die "$ERROR -ng=$GhostCell should not exceed nK/2=$nK/2\n"
	    if $GhostCell > $nK/2 and $nK>1;
    }


    print "Writing new grid size $GridSize and $GhostCell ghost cells into ".
	"$NameSizeFile and $NameBatlFile...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(MaxBlock\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	s/\b(MaxImplBLK\s*=[^0-9]*)(\d+)/$1$MaxImplBlock/i;
	print;
    }

    @ARGV = ($NameBatlFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	s/\b(nG\s*=[^0-9]*)\d/$1$GhostCell/i;
	print;
    }

}

##############################################################################

sub set_equation{

    if($Equation eq '1'){
	my @EquationModules;
	chdir $SrcEquation;
	@EquationModules = sort(glob("ModEquation?*.f90"));
	for (@EquationModules){s/^ModEquation//; s/\.f90$//;}
	print "Available Equations:\n   ",join("\n   ",@EquationModules),"\n";
	chdir "..";
	return;
    }

    my $File = "$SrcEquation/ModEquation$Equation.f90";
    die "$ERROR File $File does not exist!\n" unless -f $File;

    # Check if there is any change in the equation module
    if(-f $EquationMod){
	open(FILE1,$File); open(FILE2,$EquationMod);
	my $IsSame = 1; my $line1; my $line2;
	while($line1=<FILE1> and $line2=<FILE2>){
	    # Ignore the nWave and nMaterial definitions
	    next if $line1=~/nMaterial|nWave/ and $line2=~/nMaterial|nWave/;
	    if($line1 ne $line2){
		$IsSame = 0;
		last;
	    }
	}
	close(FILE1); close(FILE2);

	# Do not overwrite the equation module if there are no differences
	return if $IsSame;
    }

    `cp $EquationMod $EquationModSafe` if -f $EquationMod; # save previous eq.
    print "cp $File $EquationMod\n" if $Verbose;
    `cp $File $EquationMod`;
}

##############################################################################

sub set_nwave{

    $nWave = $nWaveNew;

    print "Writing new nWave = $nWaveNew into $EquationMod...\n";

    my $nWaveTwo = sprintf("%02d", $nWaveNew);

    @ARGV = ($EquationMod);

    my $prev;
    while(<>){
        if(/^\s*!/){print; next} # Skip commented out lines
        if(m/\&\s*\n/){         # Concatenate continuation lines
            $prev .= $_;
            next;
        }
        $_ = $prev . $_;
	$prev = "";
        s/\b(nWave\s*=[^0-9]*)(\d+)/$1$nWaveNew/i;
        s/I\([^\)]+\)/I($nWaveTwo)/m if /NamePrimitiveVar\s*\=/;
        print;
    }
}

#############################################################################

sub set_nmaterial{

    $nMaterial = $nMaterialNew;

    print "Writing new nMaterial = $nMaterialNew into $EquationMod...\n";

    @ARGV = ($EquationMod);

    my $prev;
    while(<>){
        if(/^\s*!/){print; next} # Skip commented out lines
        if(m/\&\s*\n/){         # Concatenate continuation lines
            $prev .= $_;
            next;
        }
        $_ = $prev . $_;
	$prev = "";
        s/\b(nMaterial\s*=[^0-9]*)(\d+)/$1$nMaterialNew/i;
        s/M\([^\)]+\)/M($nMaterialNew)/m if /NamePrimitiveVar\s*\=/;
        print;
    }
}

#############################################################################

sub set_user_module{

    if($UserModule eq '1'){
	my @UserModules;
	chdir $SrcUser;
	@UserModules = ("Default", sort(glob("*.f90")));
	for (@UserModules){s/^ModUser//; s/\.f90$//;}
	print "Available User Modules:\n   ",join("\n   ",@UserModules),"\n";
	chdir "..";
	return;
    }

    my $File;
    if($UserModule eq "Default"){
	$File = "$Src/ModUserDefault.f90";
    }else{
	$File = "$SrcUser/ModUser$UserModule.f90";
    }
    die "$ERROR File $File does not exist!\n" unless -f $File;
    return if -f $UserMod and not `diff $File $UserMod`;
    `cp $UserMod $UserModSafe` if -f $UserMod; # save previous user module
    print "cp $File $UserMod\n" if $Verbose;
    `cp $File $UserMod`;
}

#############################################################################

sub current_settings{

    $Settings = 
	"Number of cells in a block        : nI=$nI, nJ=$nJ, nK=$nK\n";
    $Settings .= 
	"Max. number of blocks/PE          : MaxBlock=$MaxBlock\n";
    $Settings .= 
	"Max. number of implicit blocks/PE : MaxImplBlock=$MaxImplBlock\n";
    $Settings .=
	"Number of ghost cell layers       : nG=$GhostCell\n";
    $Settings .=
	"Number of wave bins               : nWave=$nWave\n" if $nWave;

    $Settings .=
	"Number of materials               : nMaterial=$nMaterial\n" 
	if $nMaterial;

    open(FILE, $UserMod) or die "$ERROR Could not open $UserMod\n";
    my $Module='???';
    my $Version='???';
    while(<FILE>){
	if(/NameUserModule/){
	    $Module = $';
	    $Module = <FILE> if $Module =~ /\&\s*$/; # read continuation line
	    $Module =~ s/\s*=\s*//;   # remove equal sign
	    $Module =~ s/^\s*[\'\"]//;   # remove leading quotation mark
	    $Module =~ s/[\'\"].*\n//;  # remove trailing quotation marks
	}
	$Version = $1 if /VersionUserModule\s*=\s*([\d\.]+)/;
    }
    close(FILE);
    $Settings .= "UserModule = $Module, ver $Version\n";

    open(FILE, $EquationMod) or die "$ERROR Could not open $EquationMod\n";
    my $Equation='???';
    my $prev;
    while(<FILE>){
	if(s/\&\s*\n//){
	    $prev .= $_;
	    next;
	}
	$_ = $prev . $_;
	$prev = "";
	next unless /NameEquation\s*=\s*[\'\"]([^\'\"]*)/;
	$Equation = $1; last;
    }

    $Settings .= "Equation   = $Equation\n";

}

#############################################################################

sub print_help{

    print "
Additional options for BATSRUS/Config.pl:

-g=NI,NJ,NK,MAXBLK,MAXIMPLBLK     
                Set grid size. NI, NJ and NK are the number of cells 
                in the I, J and K directions, respectively. 
                If NK = 1, the 3rd dimension is ignored: 2D grid.
                If NJ=NK=1, the 2nd and 3rd dimensions are ignored: 1D grid.
                In non-ignored dimensions NI, NJ, NK have to be even integers.
                To allow AMR, the number of cells has to be 4 or more in all 
                non-ignored directions. 
                MAXBLK is the maximum number of blocks per processor.
                MAXIMPLBLK is the maximum number of implicitly integrated 
                blocks per processor. Cannot be larger than MAXBLK.

-ng             Print current setting for number of ghost cell layers.

-ng=GHOSTCELL   Set number of ghost cell layers around grid blocks. 
                The value GHOSTCELL has to be 2, 3, 4 or 5, but  
                not more than NI/2, NJ/2 (if NJ>1), and NK/2 (if NK>1).
                Default value is GHOSTCELL=2.

-f              Force the -g and -ng settings even if GHOSTCELL 
                exceeds NI/2, NJ/2 (if NJ>1), NK/2 (if NK>1) limits.

-e              List all available equation modules.

-e=EQUATION     Select equation EQUATION. 

-u              List all the available user modules.

-u=USERMODULE   Select the user module USERMODULE. 

-nWave=NWAVE
                Set the number of wave bins used for radiation or wave
                turbulence to NWAVE for the selected EQUATION module.

-nMaterial=NMATERIAL
                Set the number of material levels to NMATERIAL
                for the selected EQUATION module.

Examples for BATSRUS/Config.pl:

List available options for equations and user modules:

    Config.pl -e -u

Select the MHD equations, the Default user module:

    Config.pl -e=MHD -u=Default

Select the CRASH equation and user modules and 
set the number of materials to 5 and number of radiation groups to 30:

    Config.pl -e=Crash -u=Crash -nMaterial=5 -nWave=30

Set block size to 8x8x8, number of blocks to 400 and implicit blocks to 100
and number of ghost cells to 2:

    Config.pl -g=8,8,8,400,100 -ng=2

Show settings for BATSRUS:

    Config.pl -s
\n";
    exit 0;


}

