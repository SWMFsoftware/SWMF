#!/usr/bin/perl -i
my @Switch = @Arguments;
my $Code   = $CodeName;
my $MakefileDefOrig  = $OurMakefileDef;

use strict;

&print_help if not @Switch;

# Directories and files used
my $MakefileConf     = 'Makefile.conf';
my $MakefileConfOrig = 'share/build/Makefile.';
my $MakefileDef      = 'Makefile.def';
my $GridSizeScript   = 'GridSize.pl';
my $OptionScript     = 'Options.pl';

# Default precision for installation
my $DefaultPrecision = 'double';

# Warning and error messages
my $WARNING = "!!! $Code:config.pl WARNING:";
my $ERROR   = "!!! $Code:config.pl ERROR:";

# Global variables for the settings
my $Installed;             # true if code is installed ($MakefileConf exists)
my $OS='unknown';          # operating system in $MakefileDef
my $DIR='unknown';         # main directory for code in $MakefileDef
my $Compiler='unknown';    # Non default F90 compiler in $MakefileConf
my $MpiVersion='';         # Non default MPI version for mpif90.h
my $Precision='unknown';   # Precision set in $MakefileConf

# Obtain current settings
&get_settings;

# Default values for the various actions
my $Install;
my $ReInstall;
my $Show;
my $ShowLong;
my $NewPrecision;
my $DryRun;
my $GridSize;
my $Uninstall;
my $IsCompilerSet;
my $Quiet;
my $Debug;

# Set actions based on the switches
foreach (@Switch){
    if(/^-c=(.*)/)            {$Compiler=$1; $IsCompilerSet=1;  next};
    if(/^-m=(.*)/)            {$MpiVersion=$1;                  next};
    if(/^-d(ry)?(run)?$/)     {$DryRun=1;                       next};
    if(/^-h(elp)?$/i)         {&print_help};
    if(/^-p=(single|double)$/){$NewPrecision=$1;                next};
    if(/^-i(nstall)?$/)       {$Install=1;$ReInstall=/-install/;next};
    if(/^-s(how)?$/)          {$Show=1; $ShowLong=/show/;       next};
    if(/^-q(uiet)?$/)         {$Quiet=1;                        next};
    if(/^-D(ebug)?$/)         {$Debug=1;                        next};
    if(/^-uninstall$/)        {$Uninstall=1;                    next};
    if(/^-g=(.*)/)            {$GridSize.="$1,";                next};

    print "$WARNING Unknown switch: $_\n";
}

if($Uninstall){
    if(not $Installed){
	die "$ERROR: $Code is not installed.\n";
    }else{
	&shell_command("make distclean");
	exit 0;
    }
}

# Execute the actions in the appropriate order
&install_code if $Install;

&set_precision;

if($Show){
    &get_settings;
    &show_settings;
}

##############################################################################
sub get_settings{

    $Installed = -e $MakefileConf;

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = 'single';

    # Read information from $MakefileDef
    open(MAKEFILE, $MakefileDef)
	or die "$ERROR could not open $MakefileDef\n";

    while(<MAKEFILE>){
	$OS  = $1 if /^\s*OS\s*=\s*(\w+)/;
	$DIR = $1 if /^\s*SWMF_ROOT\s*=\s*(\S+)/;
    }
    close(MAKEFILE);

    # Read information from $MakefileConf
    open(MAKEFILE, $MakefileConf)
	or die "$ERROR could not open $MakefileConf\n";

    while(<MAKEFILE>){
	$Compiler = $1 if /^\s*COMPILE.f90\s*=\s*(\S+)/;
	$Precision = 'double' if 
	    /^\s*PRECISION\s*=\s*(\-r8|\-real_size\s*64|\-\-dbl|\-qautodbl=dbl4)/;
    }

    close(MAKEFILE);
}

##############################################################################

sub show_settings{

    die "SWMF is not installed\n" unless $Installed;

    print "\n";
    print "$Code is installed in directory $DIR.\n";
    print "The installation is for the $OS operating system.\n";
    print "The selected F90 compiler is $Compiler.\n";
    print "The default precision for reals is $Precision precision.\n";

    print "\n";

}

##############################################################################
sub install_code{

    if($Installed){
	print "$Code is already installed.\n";
	return unless $ReInstall; 
	print "Reinstalling $Code...\n";
    }

    # Obtain $OS and $DIR
    $OS  = `uname`    or die "$ERROR Could not obtain OS\n";
    chomp $OS;
    $DIR = `/bin/pwd` or die "$ERROR Could not obtain DIR\n";
    chomp $DIR;

    print "Creating $MakefileDef with OS=$OS, DIR=$DIR\n";
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

    # Read info from main Makefile.def
    &get_settings;

    # Install the code
    my $command = "make install";
    $command .= " COMPILER='$Compiler'" if $IsCompilerSet;
    $command .= " MPIVERSION='$MpiVersion'" if $MpiVersion;
    &shell_command($command);

    # Set initial precision for reals
    $NewPrecision = $DefaultPrecision unless $NewPrecision;
    &set_precision("init");

    # Now code is installed
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

exit 0;
##############################################################################

sub shell_command{

    my $command = join(' ',@_);
    print "$command\n" unless $Quiet;

    return if $DryRun;

    system($command)
	and die "$ERROR Could not execute command=$command\n";
}

##############################################################################
#BOP
#!QUOTE: \subsection{Installation and Configuration with config.pl}
#!ROUTINE: config.pl - (un)installation and configuration of SWMF/components
#!DESCRIPTION:
# The config.pl provides a single uniform interface towards 
# installation, configuration and uninstallation for the SWMF and its
# components.
#
#!REVISION HISTORY:
# 12/16/2006 G. Toth - initial version based on SetSWMF.pl
#EOP
sub print_help{

    print 
#BOC
"config.pl can be used for installing and setting various options for SWMF
or its components. The core of the script is in share/Scripts/config.pl,
and this is used by the config.pl scripts in the main SWMF and component 
directories. This help describes the options/features of the core script.
Additional features (if any) will be shown below.

This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.

Usage: config.pl [-h] [-q] [-D] [-d] [-s]
                  [-i[nstall] [-c=COMPILER] [-m=MPIVERSION]]
                  [-p=PRECISION] 
                  [-uninstall]

Options:

-c=COMPILER    copy Makefile.conf for a non-default F90 compiler COMPILER
               during installation

-d  -dry       dry run (do not modify anything, just show actions)

-D  -Debug     debug information shown

-h  -help      show this help message

-i  -install   install code (create Makefile.conf, Makefile.def, make install)
               The short form returns with a warning if Makefile.conf 
               already exists.
               The long form -install redoes the installation in any case.

-m=MPIVERSION  copy mpif90_OSMPIVERSION into mpif90.h during installation

-p=PRECISION   set precision (in Makefile.conf, copy mpif90.h and make clean)
               Possible values are 'single' and 'double'.

-q  -quiet     quiet execution

-s  -show      show current settings. The long form -show gives more details.

-uninstall     uninstall code (make distclean)

Examples of use:

Show current settings: 

    config.pl -s

Show current settings with more detail: 

    config.pl -show

Install code with the ifort compiler and Altix MPI and select single precision:

    config.pl -i -c=ifort -m=Altix -p=single

Reinstall code with new compiler:

    config.pl -install -c=pgf90

Uninstall code (if this fails, run config.pl -install first):

    config.pl -uninstall"
#EOC
    ,"\n\n";
    exit 0;
}

##############################################################################

1;
