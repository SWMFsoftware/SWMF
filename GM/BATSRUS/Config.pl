#!/usr/bin/perl -i
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

&print_help if $Help;

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
my $Setvar;

# Grid size variables
my $NameSizeFile = "$Src/ModSize.f90";
my $NameBatlFile = "srcBATL/BATL_size.f90";
my $GridSize;
my ($nI, $nJ, $nK, $MaxBlock);
my $MaxImplBlock;

# additional variable information
my $nWave;

# For SC/BATSRUS and IH/BATSRUS src/ is created during configuration of SWMF
if(not -d $Src){exit 0};

# Read previous grid size, equation and user module
&get_settings;

foreach (@Arguments){
    if(/^-e$/)                {$Equation=1;                    next};
    if(/^-u$/)                {$UserModule=1;                  next};
    if(/^-e=(.*)$/)           {$Equation=$1;                   next};
    if(/^-u=(.*)$/)           {$UserModule=$1;                 next};
    if(/^-s$/)                {$Show=1;                        next};
    if(/^-dynamic$/)          {`cd $Src; make DYNAMIC`;        next};
    if(/^-static$/)           {`cd $Src; make STATIC`;         next};
    if(/^-setvar=(.*)$/)      {$Setvar=$1;                     next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;
# Show grid size in a compact form if requested
print "Config.pl -g=$nI,$nJ,$nK,$MaxBlock",
    ",$MaxImplBlock"                           #^CFG IF IMPLICIT
    ,"\n" if $ShowGridSize and not $Show;


# Set or list the equations
&set_equation if $Equation;

# Set additional variable information
&set_var if $Setvar;

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

    #^CFG IF IMPLICIT BEGIN
    die "$ERROR could not read MaxImplBlock from $NameSizeFile\n" 
	unless length($MaxImplBlock);                         
    #^CFG END IMPLICIT

    open(FILE, $NameBatlFile) or die "$ERROR could not open $NameBatlFile\n";
    while(<FILE>){
        next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
	$nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
	$nK=$1           if /\bnK\s*=\s*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read nI from $NameBatlFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameBatlFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameBatlFile\n" unless length($nK);

    $GridSize = "$nI,$nJ,$nK,$MaxBlock";
    $GridSize .= ",$MaxImplBlock";                            #^CFG IF IMPLICIT

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^[1-9]\d*,[1-9]\d*,[1-9]\d*,[1-9]\d*
       ,[1-9]\d*  #^CFG IF IMPLICIT
       $/x){
	($nI,$nJ,$nK,$MaxBlock,$MaxImplBlock)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be ".
	    "5".  #^CFG IF IMPLICIT
	    #"4". #^CFG IF NOT IMPLICIT
	    " positive integers separated with commas\n";
    }

    # Check the grid size (to be set)
    die "$ERROR nK=$nK must be 1 if nJ is 1\n"         if $nJ==1 and $nK>1;
    die "$ERROR nI=$nI must be an even integer\n"      if           $nI%2!=0;
    die "$ERROR nJ=$nJ must be 1 or an even integer\n" if $nJ>1 and $nJ%2!=0;
    die "$ERROR nK=$nK must be 1 or an even integer\n" if $nK>1 and $nK%2!=0;

    warn "$WARNING nI=$nI nJ=$nJ nK=$nK does not allow AMR\n" 
	if $nI == 2 or $nJ == 2 or $nK==2;

    #^CFG IF IMPLICIT BEGIN
    die "$ERROR MaxImplBlock=$MaxImplBlock cannot exceed MaxBlock=$MaxBlock\n"
	if $MaxImplBlock > $MaxBlock;
    #^CFG END IMPLICIT

    print "Writing new grid size $GridSize into ".
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
    return if -f $EquationMod and not `diff $File $EquationMod`;
    `cp $EquationMod $EquationModSafe` if -f $EquationMod; # save previous eq.
    print "cp $File $EquationMod\n" if $Verbose;
    `cp $File $EquationMod`;
}

##############################################################################

sub set_var{

    # Create a %Setvars hash from the Setvar string
    my %Setvars;
    my @Setvars = split( /,/, $Setvar);

    my $i;
    for ($i=0; $i<=$#Setvars; ++$i){
	my %OneSetvar;
	my @OneSetvar = split( /=/, $Setvars[$i]);

        my $Var = $OneSetvar[0];
        my $Value = $OneSetvar[1];

	if($Var eq "nWave"){
	    $nWave = $Value;
            # Check the number of wave bins (to be set)
	    die "$ERROR nWave=$nWave must be 1 or more\n" if $nWave < 1;
	}else{
	    die "$ERROR unkown variable $Var in setvar\n";
	}

	print "Writing new $Var = $Value into $EquationMod...\n";

	@ARGV = ($EquationMod);

	while(<>){
	    if(/^\s*!/){print; next} # Skip commented out lines
	    s/\b($Var\s*=[^0-9]*)(\d+)/$1$Value/i;
	    print;
	}
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
    #^CFG IF IMPLICIT BEGIN
    $Settings .= 
	"Max. number of implicit blocks/PE : MaxImplBlock=$MaxImplBlock\n";
    #^CFG END IMPLICIT

    # Check if the large arrays are static or dynamic
    my $IsDynamic = `grep IsDynamic $Src/ModAdvance.f90 | grep true` ?
	"yes" : "no";
    $Settings .= 
	"Allocation of large arrays        : IsDynamic=$IsDynamic\n";

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

    $Settings     .= "Equation   = $Equation\n";

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
                1D and 2D grids are still experimental and limited.
                In non-ignored dimensions NI, NJ, NK have to be even integers.
                To allow AMR, the number of cells has to be 4 or more in all 
                non-ignored directions. 
                MAXBLK is the maximum number of blocks per processor.
                MAXIMPLBLK is the maximum number of implicitly integrated 
                blocks per processor. Cannot be larger than MAXBLK.

-e              List all available equation modules.

-e=EQUATION     Select equation EQUATION. 

-u              List all the available user modules.

-u=USERMODULE   Select the user module USERMODULE. 

-dynamic        Use dynamic allocation for large arrays.

-static         Use static allocation for large arrays.

-setvar=nWave=NWAVE
                Set additional variables for the selected EQUATION module.
                NWAVE is the number of wave bins used for radiation or wave
                turbulence. Resetting the EQUATION module will remove the
                stored NWAVE.

Examples for BATSRUS/Config.pl:

List available options for equations and user modules:

    Config.pl -e -u

Select the MHD equations, the Default user module:

    Config.pl -e=MHD -u=Default

Set block size to 8x8x8, number of blocks to 400",
" and implicit blocks to 100", #^CFG IF IMPLICIT
":

    Config.pl -g=8,8,8,400",
",100", #^CFG IF IMPLICIT
"

Show settings for BATSRUS:

    Config.pl -s
\n";
    exit 0;


}

