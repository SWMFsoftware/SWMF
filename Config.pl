#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

use strict;

# Include current dir perl modules
# Needed for perl5 change of not including current dir
use Cwd qw( abs_path );
use File::Basename qw( dirname );
use lib dirname(abs_path($0));

# Run the shared Config.pl script first
our $Component       = '';
our $Code            = 'SWMF';
our $MakefileDef     = 'Makefile.def';
our $MakefileConf    = 'Makefile.conf';
our $MakefileDefOrig = 'CON/Makefile.def';
our @Arguments = @ARGV;
our $CloneOnly;

# Figure out remote git server
my $remote = `git config remote.origin.url`; $remote =~ s/\/SWMF(\.git)?\n//;
my $SWMFsoftware = ($remote =~ /SWMFsoftware/);

my $config = "share/Scripts/Config.pl";
my $gitclone = `pwd`; chop $gitclone;
if($SWMFsoftware){
    $gitclone .= "/share/Scripts/gitclone -s";
}else{
    $gitclone .= "/share/Scripts/githubclone";
}

# Git clone missing directories as needed. Start with share/ to get $gitclone.
if (not -f $config){
    `git clone $remote/share; git clone $remote/util`;
}

#print "remote=$remote SWMFsoftware=$SWMFsoftware\n";

# Fix component Makefile.def and .conf files if they exist
&set_version_makefile_comp;

# Hash of model names with corresponding component name
my %component;

if($SWMFsoftware){
    %component = (
	"BATSRUS"       => "GM", 
	"Ridley_serial" => "IE", 
	"CIMI"          => "IM", 
	"HEIDI"         => "IM", 
	"RCM2"          => "IM",
	"ALTOR"         => "PC",
	"FLEKS"         => "PC", 
	"FLEKS_PT"      => "PT", 
	"DGCPM"         => "PS", 
	"AMPS_PT"       => "PT", 
	"PARMISAN"      => "PT",
	"PWOM"          => "PW", 
	"RBE"           => "RB",
	"MFLAMPA"       => "SP",
	"MGITM"         => "UA",
	"GITM2"         => "UA");
}else{
    %component = (
	"BATSRUS"       => "GM", 
	"Ridley_serial" => "IE", 
	"RBE"           => "RB",
	"RCM2"          => "IM");
}	

my $History;
my @models;
my $Sleep = $ENV{GITHUBSLEEP};

my $AddAmrex;
foreach (@Arguments){
    if( /^-(install|clone)/){
	$CloneOnly = 1 if /^-clone/;
	if( /^-(install|clone)=(.*)/){
	    @models = split(/,/, $2);
	    my $model;
	    foreach $model (@models){
		die "ERROR: Unknown model=$model listed in -install=...\n"
		    unless $component{$model};
	    }
	}else{
	    @models = sort(keys(%component));
	}
	next;
    }
    if( /^-history$/ or /^-date=/){$History = 1; next;}
    if( /^-sleep=(\d+)/){$Sleep = $1; next;}
    $AddAmrex = 1 if /^-v=.*FLEKS/;
    $AddAmrex = 0 if /^-amrex/;
}
# Add -amrex if FLEKS is switched on
push(@Arguments, '-amrex') if $AddAmrex;

$gitclone .= " -history" if $History;

my $repo;
foreach $repo ("util", @models){
    my $component = ($component{$repo} or ".");
    my $repo1 = $repo;
    $repo1 =~ s/AMPS_P[TC]/AMPS/; # remove _PC, _PT
    $repo1 =~ s/FLEKS_PT/FLEKS/; # remove _PT
    my $model = "$component/$repo1";
    `cd $component; $gitclone $repo1` unless -d $model;
    if($model eq "GM/BATSRUS"){
	`cd GM/BATSRUS; $gitclone srcBATL`
			    if not -d "$model/srcBATL";
	`cd GM/BATSRUS; $gitclone srcUserExtra`
			    if not -d "$model/srcUserExtra" and $SWMFsoftware;
    }
}

my $config = "share/Scripts/Config.pl";
require $config or die "Could not find $config!\n";

if($CloneOnly){
    # The share/Scripts/Config.pl is needed before this to possibly set date
    print "Done cloning git repositories\n";
    exit 0;
}

# These were set by the shared Config.pl script
our %Remaining;         # Remaining arguments after processed by $Config.pl
our $DryRun;
our $DIR;
our $Installed;
our $Show;
our $Help;
our $WARNING;
our $ERROR;
our $NewHdf5;

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
my $Options;

# Set actions based on the switches
foreach (@Arguments){
    if(/^-l(ist)?$/)          {$ListVersions=1;                 next;}
    if(/^-s$/)                {$ShowShort=1;                    next;}
    if(/^-v=(.*)/)            {push(@NewVersion,split(/,/,$1)); next;}
    if(/^-o=(.*)/)            {$Options.=",$1";                 next;}
    if(/^-sleep=\d+/)         {next;}
    
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

# Set HDF5 macro for IPIC3D2 if necessary
if($Version{PC} eq 'IPIC3D2' and ($NewHdf5 or @NewVersion)){
    &shell_command("cd PC/IPIC3D2; ./Config.pl -sethdf5");
}

if($Installed){
    # Create EE/BATSRUS if needed
    &shell_command("make EEBATSRUS")
	if $Version{"EE"} eq "BATSRUS" and not -f "EE/BATSRUS/src/Makefile";

    # Create IH/BATSRUS if needed
    &shell_command("make IHBATSRUS")
	if $Version{"IH"} eq "BATSRUS" and not -f "IH/BATSRUS/src/Makefile";

    # Create OH/BATSRUS if needed
    &shell_command("make OHBATSRUS")
	if $Version{"OH"} eq "BATSRUS" and not -f "OH/BATSRUS/src/Makefile";

    # Create SC/BATSRUS if needed
    &shell_command("make SCBATSRUS")
	if $Version{"SC"} eq "BATSRUS" and not -f "SC/BATSRUS/src/Makefile";

}

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
	    unless $compversion =~ /\/\w/ and -d $compversion;

	($comp,$version) = split(/\//,$compversion);
	if($Version{$comp} ne $version){
	    $change = 1;
	    $Version{$comp}=$version;

	    # The wrapper of the new component version should be recompiled
	    foreach my $src ("src", "srcInterface"){
		my $wrapper = "$compversion/$src/$comp"."_wrapper.f90";
		&shell_command("touch $wrapper") if -e $wrapper;
	    }
	}
    }
    return unless $change;

    # CON/Interface needs to be recompiled (uses modules from component wrappers)
    &shell_command("touch CON/Interface/src/*.f90");

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

    @Version = @NewVersion;
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
	$Option =~ s/,([a-z])/ \-$1/ig;  # Replace ',x' with ' -x' for multiple options
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

    # Fix $MakefileDef and .conf (if they exist) in all component versions

    # Collect $MakefileDef files in all component versions
    my @File;
    @File = glob("[A-Z][A-Z]/*/$MakefileDef");

    my $pwd = `pwd`; chop $pwd;
    my $File;
    foreach $File (@File){
	# Set the version based on the file name
	$File =~ m|^(([A-Z][A-Z])/[^/]+)|;
	my $Version = $1;
	my $Component = $2;
	# Write correct information into $MakefileDef
	`echo include $pwd/$MakefileDef > $File`;
	`echo MYDIR=$pwd/$Version >> $File`;
	`echo COMPONENT=$Component >> $File`;
	# If there is $MakefileDef, there should be a $MakefileConf too
	`echo include $pwd/$MakefileConf > $Version/$MakefileConf`;
    }

    # Fix DIR and OS settings
    if(-f $MakefileDef){
	@ARGV = ($MakefileDef);
	while(<>){
	    s/^((MY)?DIR\s*=).*/$1 $pwd/;
	    s/^(OS\s*=).*\n/"$1 ".`uname`/e;
	    print;
	}
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
"Additional options for SWMF/Config.pl:

-clone         Use git clone to get all missing models and/or share and util.
-clone=MOD1,MOD2 Use git clone to get the models listed.
	       The model names should be separated with commas, but no space.
	       Models that are already present will not be removed.

-install=MOD1,MOD2 Perform the same as -clone and also (re)install.

-history       Get the models from Git with full development history. 
               Default is no history to reduce size and download time.

-sleep=VALUE   Sleep VALUE number of seconds after each git clone, so
	       the server does not reject the ssh connection. Only useful
	       in combination with the -install and -clone options. 
               Default is set by the \$GITHUBSLEEP environment variable.

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

Clone share, util and all models with full version history, wait 5s in between:

    Config.pl -clone -history -sleep=5

Clone share, util and the listed models and set checkout date:

    Config.pl -clone=BATSRUS,Ridley_serial,RCM2,RBE -date='2019-01-01 19:00'

(Re)install BATSRUS and AMPS_PT, which is the AMPS model as a PT component,
with full version history:

    Config.pl -install=BATSRUS,AMPS_PT -history

Select the empty version for all components except GM/BATSRUS:

    Config.pl -v=Empty,GM/BATSRUS

Select IH/BATSRUS, IM/RCM2 and UA/GITM2 component versions:

    Config.pl -v=IH/BATSRUS -v=IM/RCM2,UA/GITM2

Set the grid size for GM and IH:

    Config.pl -o=GM:g=8,8,8,IH:g=6,6,6

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
