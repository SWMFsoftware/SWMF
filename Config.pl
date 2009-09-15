#!/usr/bin/perl -i
use strict;

# Run the shared Config.pl script first
our $Component       = '';
our $Code            = 'SWMF';
our $MakefileDefOrig = 'CON/Makefile.def';
our @Arguments = @ARGV;

my $config     = "share/Scripts/Config.pl";
require $config or die "Could not find $config!\n";

# These were set by the shared Config.pl script
our %Remaining;         # Remaining arguments after processed by $Config.pl
our $DryRun;
our $MakefileDef;
our $MakefileConf;
our $DIR;
our $Installed;
our $Show;
our $Help;
our $WARNING;
our $ERROR;

&print_help if $Help;

my $ValidComp    = '[A-Z][A-Z]';      # Pattern to match component ID-s
my $ConfigScript = 'Config.pl';

my @Version = ('unknown'); # Component versions selected in $MakefileDef
my %Version;               # Hash table to get version for a component

# Obtain current settings
&get_settings;

my $ListVersions;
my $ShowShort;
my @NewVersion;
my $GridSize;
my $Options;

# Set actions based on the switches
foreach (@Arguments){
    if(/^-l(ist)?$/)          {$ListVersions=1;                 next};
    if(/^-s$/)                {$ShowShort=1;                    next};
    if(/^-v=(.*)/)            {push(@NewVersion,split(/,/,$1)); next};
    if(/^-g=(.*)/)            {$GridSize.="$1,";                next};
    if(/^-o=(.*)/)            {$Options.=",$1";                 next};

    print "$WARNING Unknown switch: $_\n" if $Remaining{$_};
}

# Make sure that all the $Version/$MakefileDef files are correct
&set_version_makefile_comp if $Installed;

die "$ERROR SWMF is not installed.\n" unless $Installed;

if($ListVersions){
    &list_versions;
    exit 0;
}

&set_versions if @NewVersion;

if($Installed){
    # Create IH/BATSRUS if needed
    &shell_command("make IHBATSRUS")
	if $Version{"IH"} eq "BATSRUS" and not -f "IH/BATSRUS/src/Makefile";

    # Create OH/BATSRUS if needed
    &shell_command("make OHBATSRUS")
	if $Version{"OH"} eq "BATSRUS" and not -f "OH/BATSRUS/src/Makefile";

    # Create SC/BATSRUS if needed
    &shell_command("make SCBATSRUS")
	if $Version{"SC"} eq "BATSRUS" and not -f "SC/BATSRUS/src/Makefile";

    # Create LC/BATSRUS if needed
    &shell_command("make LCBATSRUS")
	if $Version{"LC"} eq "BATSRUS" and not -f "LC/BATSRUS/src/Makefile";

}

&set_grid_size if $GridSize;

&set_options if $Options;

if($Show or $ShowShort){
    &get_settings;
    &show_settings;
}

exit 0;

##############################################################################

sub list_versions{

    # Find all components from the directory structure
    my @Components;
    @Components = grep {-d and /$ValidComp/} glob("[A-Z][A-Z]");

    # Find all component versions from the directory structure
    my %Versions;
    my $Comp;
    foreach $Comp (@Components){
	$Versions{$Comp}=join(' ',grep {-d and not /\/CVS$/} glob("$Comp/*"));
    }

    # Read selected versions from $MakefileDef
    open(MAKEFILE, $MakefileDef)
	or die "$ERROR could not open $MakefileDef\n";

    while(<MAKEFILE>){
	if(/^\s*([A-Z][A-Z])_VERSION\s*=\s*(\w+)/){

	    my $Comp = $1;
	    my $Version = $2;

	    # Put the selected version to the front
	    $Versions{$Comp} =~ s/\b$Comp\/$Version\b// or
		warn "$WARNING there is no directory correspoding ".
		"to the selected version $Version\n";
	    $Versions{$Comp} = "$Comp/$Version $Versions{$Comp}";
	}
    }
    close(MAKEFILE);

    my $Comp;
    print "\nSelected version    Other versions\n","-" x 79,"\n";
    foreach $Comp (sort keys %Versions){
	my $Version;
	foreach $Version (split(' ',$Versions{$Comp})){
	    printf "%-20s","$Version";
	}
	print "\n";
    }
    print "\n";

}

##############################################################################
sub get_settings{

    return if not $Installed;

    # Set defaults/initial values
    @Version = ();

    # Read information from $MakefileDef
    open(MAKEFILE, $MakefileDef)
	or die "$ERROR could not open $MakefileDef\n";

    while(<MAKEFILE>){
	if(/^\s*([A-Z][A-Z])_VERSION\s*=\s*(\w+)/){
	    # Store version in both the array and the hash table
	    push(@Version, "$1/$2");
	    $Version{$1}=$2;
	}
    }
    close(MAKEFILE);
}

##############################################################################

sub show_settings{

    die "SWMF is not installed\n" unless $Installed;

    print "The selected component versions and settings are:\n";
    my $Comp;
    foreach $Comp (sort keys %Version){
	if($Show){
	    print "\n$Comp/$Version{$Comp}";
	}else{
	    printf "%-30s","\t$Comp/$Version{$Comp}";
	}
	if(-x "$Comp/$Version{$Comp}/$ConfigScript"){
	    if($Show){
		print ":\n";
		print `cd $Comp/$Version{$Comp}; ./$ConfigScript -s`;
	    }else{
		my $Grid = `cd $Comp/$Version{$Comp}; ./$ConfigScript -g`;
		chop $Grid; $Grid =~ s/.*=//;
		print "grid: $Grid" if $Grid;
	        print "\n";
	    }
	}else{
	    print "\n";
	}
    }
    print "\n";
}

##############################################################################
sub set_versions{

    # Check if there is any change and whether the version directory exists
    my $compversion; # Stores information in COMPONENT/VERSION format
    my $comp;        # COMPONENT
    my $version;     # VERSION
    my $change;      # true if some version has changed

    if($NewVersion[0] =~ /^empty$/i){
	# Remove the first element of the array
	shift(@NewVersion);

	# select the Empty version for all components that are not set
	foreach $comp (keys %Version){
	    push(@NewVersion, "$comp/Empty") 
		unless grep(/^$comp\//, @NewVersion);
	}
    }
    foreach $compversion (@NewVersion){

	die("$ERROR directory $compversion does not exist.")
	    unless -d $compversion;

	($comp,$version) = split(/\//,$compversion);
	if($Version{$comp} ne $version){
	    $change = 1;
	    $Version{$comp}=$version;
	}
    }

    die "$ERROR IH/BATSRUS_share requires GM/BATSRUS(_conf) ".
	"and not GM/$Version{GM}\n"
        if $Version{IH} eq 'BATSRUS_share' and $Version{GM} !~ /^BATSRUS/;

    die "$ERROR non Empty UA version requires non Empty IE version\n"
        if $Version{UA} ne 'Empty' and $Version{IE} eq 'Empty';

    return unless $change;

    print "Modifying versions in $MakefileDef\n";
    return if $DryRun;

    # Set component versions in $MakefileDef
    @ARGV = ($MakefileDef);
    my %Found;
    while(<>){
	$Found{$2} = 1 if
	    s/^\s*(($ValidComp)_VERSION)\s*=\s*\w+/$1 = $Version{$2}/;
	print;
    }

    # Check if all components have been found
    foreach $comp (keys %Version){
	die "$ERROR could not find ${comp}_VERSION line in $MakefileDef\n"
	    unless $Found{$comp}
    }

    # Set equation and user routines for IH/BATSRUS_share if needed
    &shell_command("cd GM/BATSRUS; ./Config.pl -u=Ih -e=Mhd")
	if $Version{"IH"} eq "BATSRUS_share";

    @Version = @NewVersion;
}

##############################################################################

sub set_grid_size{

    # Create a %GridSize{CompID} hash from the GridSize string
    my %GridSize;
    while( $GridSize =~ s/^($ValidComp)[:,]([\d,]+)// ){
	my $Comp=$1;
	my $Grid=$2;
	$Grid =~ s/,$//;            # remove trailing comma
	$GridSize{$Comp} = $Grid;
    }
    
    die "$ERROR -g=...$GridSize has incorrect format\n" if $GridSize;

    # Set the size of the grid for all the components listed in %GridSize
    my $Comp;
    foreach $Comp (sort keys %GridSize){
	die "$ERROR -g=..$Comp:.. version is not known for this component\n"
	    unless $Version{$Comp};

	&shell_command("cd $Comp/$Version{$Comp}; ".
		       "./$ConfigScript -g=$GridSize{$Comp}");
    }
}

##############################################################################

sub set_options{

    # Create a %Options{CompID} hash from the Options string
    my %Options;
    my @Options = split( /,($ValidComp):/, $Options);

    die "$ERROR -o=...$Options has incorrect format\n" unless @Options;

    my $i;
    for ($i=1; $i<=$#Options; $i+=2){
	my $Comp = $Options[$i];
	my $Option = $Options[$i+1];
	$Option =~ s/,/ \-/g;   # Replace comma with ' -' for multiple options
	$Options{$Comp}=$Option;
    }

    # Set the options for all the components listed in %Options
    my $Comp;
    foreach $Comp (sort keys %Options){
	die "$ERROR -o=..$Comp:.. version is not known for this component\n"
	    unless $Version{$Comp};

	&shell_command("cd $Comp/$Version{$Comp}; ".
		       "./$ConfigScript -$Options{$Comp}");
    }
}

##############################################################################

sub set_version_makefile_comp{

    # Set $MakefileDef (if it exists) in all component versions

    # Collect all $MakefileDef files in all component versions
    my @File;
    @File = glob("[A-Z][A-Z]/*/$MakefileDef");

    my $File;
    foreach $File (@File){
	# Set the version based on the file name
	$File =~ m|^([A-Z][A-Z]/[^/]+)|;
	my $Version = $1;
	# Put the proper include command into $MakefileDef and $MakefileConf
	&shell_command("echo include $DIR/$MakefileDef > $File");
	&shell_command("echo include $DIR/$MakefileConf > ".
		       "$Version/$MakefileConf");
    }
}

##############################################################################
#!QUOTE: \chapter{Makefiles and Scripts}
#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{The SMWF Scripts}
#!QUOTE: \subsection{Installation and Configuration with Config.pl and Configure.pl}
#!ROUTINE: Config.pl - (un)installation and configuration of SWMF
#!DESCRIPTION:
# The Config.pl provides a single uniform interface towards 
# installation, configuration and uninstallation.
#
#!REVISION HISTORY:
# 10/29/2003 G. Toth - initial version
#                      several extensions and modifications
# 01/20/2007           renamed from SetSWMF.pl to Config.pl and the
#                      core of the script is moved into share/Scripts/Config.pl
#06/01/2008  R. Oran   Modifications to include new OH component
#EOP
sub print_help{

    print 
#BOC
"Additional options for SWMF/config.pl:

-g=ID:GRIDSIZE set the size of the grid to GRIDSIZE for the component 
               identified with ID. This flag can occur multiple times and/or 
               multiple grid sizes can be given in a comma separated list.

-o=ID:OPTIONS  set options OPTIONS for the component identified with ID.
               This flag can occur multiple times and/or options for
               multiple components can be given in a comma separated list.

-l  -list      List available component versions

-v=VERSION     select component version VERSION. This flag can occur 
               multiple times and/or multiple versions can be given
               in a comma separated list.
               If IH/BATSRUS or SC/BATSRUS is selected for the first time, 
               make IHBATSRUS or make SCBATSRUS is done, respectively.
               If the first VERSION is set to 'Empty' without specifying the
               component, the Empty version is selected for all the components
               that are not listed explicitly.

Examples for the SWMF Config.pl:

Select the empty version for all components except GM/BATSRUS:

    Config.pl -v=Empty,GM/BATSRUS

Select IH/BATSRUS, IM/RCM2 and UA/GITM2 component versions:

    Config.pl -v=IH/BATSRUS -v=IM/RCM2,UA/GITM2

Set the grid size for GM and IH:

    Config.pl -g=GM:8,8,8,400,100 -g=IH:6,6,6,800,1

Show the currently set options for the GM and SC component:

    Config.pl -o=GM:s,SC:s

Show the available user modules and equations for GM and IH

    Config.pl -o=GM:e,u -o=IH:e,u

Set equation and user module for GM:

    Config.pl -o=GM:e=Mhd,u=Default"
#EOC
    ,"\n\n";
    exit 0;
}

##############################################################################
