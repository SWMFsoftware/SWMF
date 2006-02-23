#!/usr/bin/perl -i
use strict;

# Pattern to match component ID-s
my $ValidComp = 'IH|GM|IE|IM|RB|SC|SP|UA';

&print_help if not @ARGV;

# Directories and files used
my $MakefileConf     = 'Makefile.conf';
my $MakefileConfOrig = 'share/build/Makefile.';
my $MakefileDef      = 'Makefile.def';
my $MakefileDefOrig  = 'CON/Makefile.def';
my $MakefileDefCvs   = 'Makefile.def.orig';
my $GridSizeScript   = 'GridSize.pl';
my $OptionScript     = 'Options.pl';

# Default precision for installation
my $DefaultPrecision = 'double';

# Warning and error messages
my $WARNING = "!!! SetSWMF_WARNING:";
my $ERROR   = "!!! SetSWMF_ERROR:";

# Global variables for the settings
my $Installed;             # true if SWMF is installed ($MakefileConf exists)
my $OS='unknown';          # operating system in $MakefileDef
my $DIR='unknown';         # main directory for SWMF in $MakefileDef
my $Compiler='unknown';    # Non default F90 compiler in $MakefileConf
my $MpiVersion='';         # Non default MPI version for mpif90.h
my $Precision='unknown';   # Precision set in $MakefileConf
my @Version;               # Component versions selected in $MakefileDef
@Version = ('unknown');
my %Version;               # Hash table to get version for a component

# Obtain current settings
&get_settings;

# Default values for the various actions
my $Install;
my $ReInstall;
my $Show;
my $ShowLong;
my $ListVersions;
my $NewPrecision;
my $DryRun;
my @NewVersion;
my $GridSize;
my $Options;
my $Uninstall;
my $IsCompilerSet;
my $Quiet;
my $Debug;

# Set actions based on the switches
my @switch = @ARGV;
foreach (@switch){
    if(/^-c=(.*)/)            {$Compiler=$1; $IsCompilerSet=1;  next};
    if(/^-m=(.*)/)            {$MpiVersion=$1;                  next};
    if(/^-d(ry)?(run)?$/)     {$DryRun=1;                       next};
    if(/^-h(elp)?$/i)         {&print_help};
    if(/^-p=(single|double)$/){$NewPrecision=$1;                next};
    if(/^-i(nstall)?$/)       {$Install=1;$ReInstall=/-install/;next};
    if(/^-l(ist)?$/)          {$ListVersions=1;                 next};
    if(/^-s(how)?$/)          {$Show=1; $ShowLong=/show/;       next};
    if(/^-q(uiet)?$/)         {$Quiet=1;                        next};
    if(/^-D(ebug)?$/)         {$Debug=1;                        next};
    if(/^-v=(.*)/)            {push(@NewVersion,split(/,/,$1)); next};
    if(/^-uninstall$/)        {$Uninstall=1;                    next};
    if(/^-g=(.*)/)            {$GridSize.="$1,";                next};
    if(/^-o=(.*)/)            {$Options.=",$1";                 next};

    print "$WARNING Unknown switch: $_\n";
}

# Make sure that all the $Version/$MakefileDef files are correct
&set_version_makefile_comp if $Installed;

if($Uninstall){
    if(not $Installed){
	die "$ERROR SWMF is not installed.\n";
    }else{
	&shell_command("make distclean");
	exit 0;
    }
}

# Execute the actions in the appropriate order
&install_swmf if $Install;

die "$ERROR SWMF is not installed.\n" unless $Installed;

if($ListVersions){
    &list_versions;
    exit 0;
}

&set_precision;

&set_versions if @NewVersion;

&set_grid_size if $GridSize;

&set_options if $Options;

if($Show){
    &get_settings;
    &show_settings;
}

exit 0;
##############################################################################

sub list_versions{

    # Read information from $MakefileDef
    open(MAKEFILE, $MakefileDef)
	or die "$ERROR could not open $MakefileDef\n";

    my %Versions;
    while(<MAKEFILE>){
	if(/^(\#)?\s*([A-Z][A-Z])_VERSION\s*=\s*(\w+)/){
	    # Store version in both the array and the hash table
	    if(not $1 and $Versions{$2}){
		# Put the selected version to the front
		$Versions{$2} ="$3,$Versions{$2}";
	    }else{
		# Append the version
		$Versions{$2}.=$3.",";
	    }
	}
    }
    close(MAKEFILE);

    my $Comp;
    print "\nSelected version    Other versions\n","-" x 79,"\n";
    foreach $Comp (sort keys %Versions){
	my $Version;
	foreach $Version (split(',',$Versions{$Comp})){
	    printf "%-20s","$Comp/$Version";
	}
	print "\n";
    }
    print "\n";

}

##############################################################################
sub get_settings{

    $Installed = -e $MakefileConf;

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = 'single';
    @Version = ();

    # Read information from $MakefileDef
    open(MAKEFILE, $MakefileDef)
	or die "$ERROR could not open $MakefileDef\n";

    while(<MAKEFILE>){
	$OS  = $1 if /^\s*OS\s*=\s*(\w+)/;
	$DIR = $1 if /^\s*SWMF_ROOT\s*=\s*(\S+)/;
	if(/^\s*([A-Z][A-Z])_VERSION\s*=\s*(\w+)/){
	    # Store version in both the array and the hash table
	    push(@Version, "$1/$2");
	    $Version{$1}=$2;
	}
    }
    close(MAKEFILE);

    # Read information from $MakefileConf
    open(MAKEFILE, $MakefileConf)
	or die "$ERROR could not open $MakefileConf\n";

    while(<MAKEFILE>){
	$Compiler = $1 if /^\s*COMPILE.f90\s*=\s*(\S+)/;
	$Precision = 'double' if /^\s*PRECISION\s*=\s*(\-r8|\-real_size\s*64|\-\-dbl)/;
    }

    close(MAKEFILE);
}

##############################################################################

sub show_settings{

    die "SWMF is not installed\n" unless $Installed;

    print "\n";
    print "The SWMF is installed in directory $DIR.\n";
    print "The installation is for the $OS operating system.\n";
    print "The selected F90 compiler is $Compiler.\n";
    print "The default precision for reals is $Precision precision.\n";
    print "The selected component versions and grid sizes are:\n";
    my $Comp;
    foreach $Comp (sort keys %Version){
	if($ShowLong){
	    print "\n$Comp/$Version{$Comp}";
	}else{
	    printf "%-30s","\t$Comp/$Version{$Comp}";
	}
	if(-x "$Comp/$Version{$Comp}/$GridSizeScript"){
	    if($ShowLong){
		my $Grid = `cd $Comp/$Version{$Comp}; ./$GridSizeScript -s`;
		print ":\n$Grid" if $Grid;
		print `cd $Comp/$Version{$Comp}; ./$OptionScript -s`
		    if -x "$Comp/$Version{$Comp}/$OptionScript";
	    }else{
		my $Grid = `cd $Comp/$Version{$Comp}; ./$GridSizeScript`;
		$Grid =~ s/$GridSizeScript(\s*-g=)?//;
		print "grid: $Grid" if $Grid;
	    }
	}else{
	    print "\n";
	}
    }
    print "\n";


}

##############################################################################

sub install_swmf{

    if($Installed){
	print "SWMF is already installed.\n";
	return unless $ReInstall; 
	print "Reinstalling SWMF...\n";
    }

    # Obtain $OS and $DIR
    $OS  = `uname`    or die "$ERROR Could not obtain OS\n";
    chomp $OS;
    $DIR = `/bin/pwd` or die "$ERROR Could not obtain DIR\n";
    chomp $DIR;

    print "Creating $MakefileDef with OS=$OS, SWMF_ROOT=$DIR\n";
    if(not $DryRun){
	open(MAKEFILE,">$MakefileDef")
	    or die "$ERROR Could not open $MakefileDef for writing\n";

	print MAKEFILE "OS = $OS\nSWMF_ROOT = $DIR\n";
	close(MAKEFILE);
    }

    &shell_command("cat $MakefileDefOrig >> $MakefileDef");

    # Create $MakefileConf
    my $makefile = $MakefileConfOrig.$OS;
    $makefile .= $Compiler if $IsCompilerSet;
    die("$ERROR $makefile does not exist.") unless -e $makefile;
    &shell_command("cp $makefile $MakefileConf");

    # Read version and other info from main Makefile.def
    &get_settings;

    # Make sure that all the $Version/$MakefileDef files are correct
    &set_version_makefile_comp;

    # Initialize CON and the components
    my $command = "make install";
    $command .= " COMPILER='$Compiler'" if $IsCompilerSet;
    $command .= " MPIVERSION='$MpiVersion'" if $MpiVersion;
    &shell_command($command);

    # Set initial precision for reals
    $NewPrecision = $DefaultPrecision unless $NewPrecision;
    &set_precision("init");

    # Now SWMF is installed
    $Installed = 1 unless $DryRun;
}

##############################################################################

sub set_precision{

    # Set the precision for reals
    # If called with a true argument, initial setting is done, otherwise
    # an existing setting is changed

    my $init = $_[0];

    # Return if there is nothing to do
    return unless $init or ($NewPrecision and $NewPrecision ne $Precision);

    # Precision will be NewPrecision after changes
    $Precision = $NewPrecision;

    # clean the distribution unless initial installation is done
    &shell_command("make clean") unless $init;

    print "Setting PRECISION variable to $Precision precision in ".
	"$MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    if($OS eq "OSF1"){
		# Comment out line if in conflict with new precision
                $_="#$_" if /^\s*PRECISION\s*=\s*\-real_size\s+64/ and $Precision eq 'single'
                    or      /^\s*PRECISION\s*=\s*\-real_size\s+32/ and $Precision eq 'double';

                # Uncomment line if agrees with new precision
                s/#// if /^\s*#\s*PRECISION\s*=\s*\-real_size\s+64/ and $Precision eq 'double'
                    or   /^\s*#\s*PRECISION\s*=\s*\-real_size\s+32/ and $Precision eq 'single';

	    }else{
		# Comment out line if in conflict with new precision
		$_="#$_" if /^\s*PRECISION\s*=\s*\-/ and $Precision eq 'single'
		    or      /^\s*PRECISION\s*=\s*$/  and $Precision eq 'double';
		# Uncomment line if agrees with new precision
		s/#// if /^\s*#\s*PRECISION\s*=\s*\-/ and $Precision eq 'double'
		    or   /^\s*#\s*PRECISION\s*=\s*$/  and $Precision eq 'single';
	    }
	    print;
	}
    }
}

##############################################################################

sub set_versions{

    # Check if there is any change and whether the version directory exists
    my $compversion; # Stores information in COMPONENT/VERSION format
    my $comp;        # COMPONENT
    my $version;     # VERSION
    my $change;      # true if some version has changed
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
	# Skip uninteresting lines
	if(not /_VERSION/){print; next};

	# Check line against all the new versions
	foreach $compversion (@NewVersion){

	    ($comp,$version) = split(/\//,$compversion);
	    $comp.='_VERSION';

	    # Comment out all versions for this component
	    $_="#".$_ if /^\s*$comp/;
	    # Uncomment the good version
	    if(/^\s*#\s*$comp\s*=\s$version\b/){
	       s/#//;
	       $Found{$compversion}=1;
	   }
	}
	print;
    }

    # Check if all new versions have been found
    foreach $compversion (@NewVersion){
	die "$ERROR could not find version $compversion in $MakefileDef\n"
	    unless $Found{$compversion}
    }

    # Create IH/BATSRUS if needed
    &shell_command("make IHBATSRUS")
	if $Version{"IH"} eq "BATSRUS" and not -f "IH/BATSRUS/src/Makefile";

    # Create SC/BATSRUS if needed
    &shell_command("make SCBATSRUS")
	if $Version{"SC"} eq "BATSRUS" and not -f "SC/BATSRUS/src/Makefile";

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

	my $Script = "$Comp/$Version{$Comp}/$GridSizeScript";
	die "$ERROR -g=..$Comp:.. no executable $Script was found\n"
	    unless -x $Script;

	&shell_command("cd $Comp/$Version{$Comp}; ".
		       "./$GridSizeScript -g=$GridSize{$Comp}");
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

	my $Script = "$Comp/$Version{$Comp}/$OptionScript";
	die "$ERROR -o=..$Comp:.. no executable $Script was found\n"
	    unless -x $Script;

	&shell_command("cd $Comp/$Version{$Comp}; ".
		       "./$OptionScript -$Options{$Comp}");
    }
}

##############################################################################

sub set_version_makefile_comp{

    # Set $MakefileDef (if it exists) in all component versions

    # Only echo shell commands in debug mode
    my $QuietOrig = $Quiet; $Quiet = not $Debug;

    # Collect all $MakefileDef files in all component versions
    my @File;
    @File = glob("[A-Z][A-Z]/*/$MakefileDef");

    my $File;
    foreach $File (@File){
	# Set the version based on the file name
	$File =~ m|^([A-Z][A-Z]/[^/]+)|;
	my $Version = $1;
	# Save original file if not yet saved
	&shell_command("mv $File $Version/$MakefileDefCvs")
	    unless -f "$Version/$MakefileDefCvs";
	# Put the proper include command into $MakefileDef and $MakefileConf
	&shell_command("echo include $DIR/$MakefileDef > $File");
	&shell_command("echo include $DIR/$MakefileConf > ".
		       "$Version/$MakefileConf");
    }
    $Quiet = $QuietOrig;
}

##############################################################################

sub shell_command{

    my $command = join(' ',@_);
    print "$command\n" unless $Quiet;

    return if $DryRun;

    system($command)
	and die "$ERROR Could not execute command=$command\n";
}

##############################################################################
#!QUOTE: \chapter{Makefiles and Scripts}
#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{The SMWF Scripts}
#!QUOTE: \subsection{Installation and Configuration with SetSWMF.pl and Configure.pl}
#!ROUTINE: SetSWMF.pl - (un)installation and configuration of SWMF
#!DESCRIPTION:
# The SetSWMF.pl provides a single uniform interface towards 
# installation, configuration and uninstallation.
#
#!REVISION HISTORY:
# 10/29/2003 G. Toth - initial version
#                      several extensions and modifications
#EOP
sub print_help{

    print 
#BOC
"SetSWMF.pl can be used for installing and setting various options for SWMF.
This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.
This script will also be used by the GUI to interact with SWMF.

Usage: SetSWMF.pl [-h] [-q] [-D] [-d] [-l] [-s]
                  [-i[nstall] [-c=COMPILER] [-m=MPIVERSION]]
                  [-p=PRECISION] 
                  [-v=VERSION[,VERSION2,...] 
                  [-g=ID:GRIDSIZE[,ID:GRIDSIZE,...]
                  [-o=ID:OPTIONS[,ID:OPTIONS,...]
                  [-uninstall]

Options:

-c=COMPILER    copy Makefile.conf for a non-default F90 compiler COMPILER
               during installation

-d  -dry       dry run (do not modify anything, just show actions)

-D  -Debug     debug information shown

-g=ID:GRIDSIZE set the size of the grid to GRIDSIZE for the component 
               identified with ID. This flag can occur multiple times and/or 
               multiple grid sizes can be given in a comma separated list.

-o=ID:OPTIONS  set options OPTIONS for the component identified with ID.
               This flag can occur multiple times and/or options for
               multiple components can be given in a comma separated list.

-h  -help      show this help message

-i  -install   install SWMF (create Makefile.conf, Makefile.def, make install)
               The short form returns with a warning if Makefile.conf 
               already exists.
               The long form -install redoes the installation in any case.

-l  -list      List available component versions

-m=MPIVERSION  copy mpif90_OSMPIVERSION into mpif90.h during installation

-p=PRECISION   set precision (in Makefile.conf, copy mpif90.h and make clean)
               Possible values are 'single' and 'double'.

-q  -quiet     quiet execution

-s  -show      show current settings. The long form -show gives more details.

-uninstall     uninstall SWMF (make distclean)

-v=VERSION     select component version VERSION. This flag can occur 
               multiple times and/or multiple versions can be given
               in a comma separated list.
               If IH/BATSRUS or SC/BATSRUS is selected for the first time, 
               make IHBATSRUS or make SCBATSRUS is done, respectively.

Examples of use:

Show current settings: 

    SetSWMF.pl -s

Show current settings with more details: 

    SetSWMF.pl -show

Install SWMF with the ifort compiler and Altix MPI and select single precision:

    SetSWMF.pl -i -c=ifort -m=Altix -p=single

Reinstall SWMF with new compiler:

    SetSWMF.pl -install -c=pgf90

Uninstall SWMF (if this fails, run SetSWMF.pl -install first):

    SetSWMF.pl -uninstall

Select IH/BATSRUS, IM/Empty and UA/Empty component versions:

    SetSWMF.pl -v=IH/BATSRUS -v=IM/Empty,UA/Empty

Set the grid size for GM and IH:

    SetSWMF.pl -g=GM:8,8,8,400,100 -g=IH:6,6,6,800,1

Show the currently set for the GM and SC component:

    SetSWMF.pl -o=GM:show,SC:show

Show the available user modules and equations for GM and IH

    SetSWMF.pl -o=GM:equation,usermodule -o=IH:equation,usermodule

Set equation and user module for GM:

    SetSWMF.pl -o=GM:e=Mhd,u=Default"
#EOC
    ,"\n\n";
    exit 0;
}

##############################################################################
