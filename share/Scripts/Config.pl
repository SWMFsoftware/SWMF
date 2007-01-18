#!/usr/bin/perl -i
use strict;

# Default compiler per machine or OS
my %Compiler = ("Linux"   => "f95",
	     "Darwin"  => "f95",
	     "OSF1"    => "f90",
	     "IRIX64"  => "f90",
	     "palm"    => "ifort",
	     "cfe1"    => "ifort",
	     "cfe2"    => "ifort",
	     "cfe3"    => "ifort"
	     );

# Default MPI library per machine or OS
my %MpiVersion = ("Linux"   => "mpich",
	"Darwin"  => "mpich",
	"OSF1"    => "mpich",
	"IRIX64"  => "SGI",
	"palm"    => "ifort",
	"cfe1"    => "ifort",
	"cfe2"    => "ifort",
	"cfe3"    => "ifort"
	);

my $WARNING_='share/Scripts/config.pl WARNING:';
my $ERROR_  ='share/Scripts/config.pl ERROR:';

# Obtain $OS, $DIR, and the machine name and provide it to caller script
our $OS  = `uname`    or die "$ERROR_ could not obtain OS\n"; chop $OS;
our $DIR = `/bin/pwd` or die "$ERROR_ could not obtain DIR\n"; chop $DIR;
our $Machine = `hostname -s`; chop $Machine;

# These are either obtained from the calling script or set here
our $Component;             # The SWMF component the code is representing
our $Code;                  # The name of the code
($Component, $Code) = ($DIR =~ /([A-Z][A-Z])\/([^\/]+)$/)
    unless $Code;

# Strings for the error and warning messages for the caller script
our $ERROR   = "$Code/config.pl ERROR:";
our $WARNING = "$Code/config.pl WARNING:";

# Obtain the default compiler for this machine / OS
our $Compiler;
$Compiler = $Compiler{$Machine} or $Compiler = $Compiler{$OS} or
    die "$ERROR_ default compiler is not known for OS=$OS\n";

# Obtain the default MPI version for mpif90.h
our $MpiVersion;
$MpiVersion = $MpiVersion{$Machine} or $MpiVersion = $MpiVersion{$OS} or
    die "$ERROR_ default MPI version is not known for OS=$OS\n";

# These are always obtained from the calling script
our $MakefileDefOrig;       # Original Makefile.def 
our @Arguments;             # Arguments obtained from the caller script

# The arguments not handled by this script are provided to the caller
our %Remaining;

# These file names are provided to the calling script
our $MakefileDef      = 'Makefile.def';
our $MakefileConf     = 'Makefile.conf';
our $MakefileConfOrig = 'share/build/Makefile.';
our $MpiHeader        = 'share/Library/src/mpif90.h';
our $MpiHeaderOrig    = 'share/include/mpif90_';

# These options are set here and provided to the calling script
our $Help;                  # Print help message
our $Verbose;               # Verbose information is printed if true
our $Show;                  # Show information
our $DryRun;                # True if no change is actually made
our $Precision='unknown';   # Precision set in $MakefileConf
our $Installed;             # true if code is installed ($MakefileConf exists)
our $Install;               # True if code is (re)installed
our $Uninstall;             # True if code is uninstalled
our $ShowGridSize;          # Show grid size for caller code
our $NewGridSize;           # New grid size to be set in caller code

# Default precision for installation
my $DefaultPrecision = 'double';

# Global variables for the settings
my $IsComponent=0;         # True if code is installed as a component of SWMF
my $NewPrecision;
my $NewOptimize;
my $NewDebug;
my $IsCompilerSet;
my $Debug;
my $Optimize;

# Obtain current settings
&get_settings_;

# Show current settings if no arguments are given.
$Show = 1 if not @Arguments;

# Set actions based on the switches
foreach (@Arguments){
    if(/^-dryrun$/)           {$DryRun=1;                       next};
    if(/^-verbose$/i)         {$Verbose=1;                      next};
    if(/^-h(elp)?$/i)         {$Help=1;                         next};
    if(/^-show$/i)            {$Show=1;                         next};
    if(/^-(single|double)$/i) {$NewPrecision=lc($1);            next};
    if(/^-install(=.*)?$/)    {my $value=$1;
                               $IsComponent=1 if $value =~ /^=c/i;
			       $IsComponent=0 if $value =~ /^=s/i;
			       $Install=1;                      next};
    if(/^-uninstall$/i)       {$Uninstall=1;                    next};
    if(/^-compiler=(.*)$/i)   {$Compiler=$1; $IsCompilerSet=1;  next};
    if(/^-mpi=(.*)$/i)        {$MpiVersion=$1;                  next};
    if(/^-standalone$/i)      {$IsComponent=0;                  next};
    if(/^-component$/i)       {$IsComponent=1;                  next};
    if(/^-debug$/i)           {$NewDebug="yes";                 next};
    if(/^-nodebug$/i)         {$NewDebug="no";                  next};
    if(/^-O[0-4]$/i)          {$NewOptimize=$_;                 next};  
    if(/^-g(rid)?$/)          {$ShowGridSize=1;                 next};
    if(/^-g(rid)?=([\d,]+)$/) {$NewGridSize=$+;                 next};

    $Remaining{$_}=1;
}

&print_help_ if $Help;

if($Uninstall){
    if(not $Installed){
	warn "$ERROR_ $Code is not installed.\n";
	exit 1;
    }else{
	&shell_command("cd share; make distclean")
	    if -d "share" and not $IsComponent;
	&shell_command("make distclean");
	exit 0;
    }
}

# Execute the actions in the appropriate order
&install_code_ if $Install;

# Change precision of reals if required
if($NewPrecision and $NewPrecision ne $Precision){
    &shell_command("make clean");
    &set_precision_;
}

# Change debugging flags if required
&set_debug_ if $NewDebug and $NewDebug ne $Debug;

# Change optimization level if required
&set_optimization_ if $NewOptimize and $NewOptimize ne $Optimize;

if($Show){
    &get_settings_;
    &show_settings_;
}

##############################################################################
sub get_settings_{

    $Installed = (-e $MakefileConf and -e $MakefileDef);

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = "single";

  TRY:{
      # Read information from $MakefileDef
      open(MAKEFILE, $MakefileDef)
	  or die "$ERROR_ could not open $MakefileDef\n";

      while(<MAKEFILE>){
	  if(/^\s*include\s+(.*$MakefileDef)\s*$/){
	      $MakefileDef = $1;
	      $IsComponent = 1;
	      close MAKEFILE;
	      redo TRY;
	  }
	  $OS         = $1 if /^\s*OS\s*=\s*(\w+)/;
      }
      close(MAKEFILE);
  }

    $Debug = "no";
  TRY:{
      # Read information from $MakefileConf
      open(MAKEFILE, $MakefileConf)
	  or die "$ERROR_ could not open $MakefileConf\n";

      while(<MAKEFILE>){
	  if(/^\s*include\s+(.*$MakefileConf)\s*$/){
	      $MakefileConf = $1;
	      close MAKEFILE;
	      redo TRY;
	  }
	  $Compiler = $+ if /^\s*COMPILE.f90\s*=\s*(\$\{CUSTOMPATH_F\})?(\S+)/;
	  $Precision = lc($1) if /^\s*PRECISION\s*=\s*(SINGLE|DOUBLE)PREC/;
          $Debug = "yes" if /^\s*DEBUG\s*=\s*\$\{DEBUGFLAG\}/;
          $Optimize = $1 if /^\s*OPT[0-4]\s*=\s*(-O[0-4])/;
      }
  }
    close(MAKEFILE);
}

##############################################################################

sub show_settings_{

    if(not $Installed){
	print "$Code is not installed\n";
	exit 0;
    }

    print "\n";
    if($IsComponent){
	print "$Code is installed in directory $DIR the $Component component.\n";
    }else{
	print "$Code is installed in directory $DIR.\n";
    }
    print "The installation is for the $OS operating system.\n";
    print "The selected F90 compiler is $Compiler.\n";
    print "The default precision for reals is $Precision precision.\n";
    print "The maximum optimization level is $Optimize\n";
    print "Debugging flags: $Debug\n";

    print "\n";

}

##############################################################################
sub install_code_{

    my $Text = $Installed ? "Reinstalling $Code" : "Installing $Code";
    $Text .= " as an SWMF component" if $IsComponent;  
    print "$Text\n";

    if($IsComponent){
	my $dir = $DIR; $dir =~ s|/[^/]*/[^/]*$||;  # go two directories up
	my $makefile = "$dir/$MakefileDef";          # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > $MakefileDef");

	$makefile = "$dir/$MakefileConf"; # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > $MakefileConf");
    }else{
	die "$ERROR_ original $MakefileDef is not given\n" unless
	    $MakefileDefOrig;
	die "$ERROR_ $MakefileDefOrig is missing\n" unless
	    -f $MakefileDefOrig;
	&shell_command("echo OS=$OS > $MakefileDef");
	&shell_command("echo SWMF_ROOT=$DIR >> $MakefileDef");
	&shell_command("echo ${Component}DIR=$DIR >> $MakefileDef");
	&shell_command("cat $MakefileDefOrig >> $MakefileDef");
	&shell_command("cat $MakefileConfOrig$OS.$Compiler > $MakefileConf");
	&shell_command("cat $MpiHeaderOrig${OS}_$MpiVersion.h > $MpiHeader");
    }

    # Read info from main Makefile.def
    &get_settings_;

    # Set initial precision for reals
    $NewPrecision = $DefaultPrecision unless $NewPrecision;
    &set_precision_;

    # Install the code
    &shell_command("cd share; make install") 
	if -d "share" and not $IsComponent;
    &shell_command("make install");

    # Now code is installed
    $Installed = 1 unless $DryRun;
}

##############################################################################

sub set_precision_{

    # Set the precision for reals in $MakefileConf

    # Precision will be NewPrecision after changes
    $Precision = $NewPrecision;

    my $PREC = uc($Precision)."PREC";
    print "Setting PRECISION variable to $PREC in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    s/^(\s*PRECISION\s*=\s*)(SINGLE|DOUBLE)PREC/$1$PREC/;
	    print;
	}
    }
}

##############################################################################

sub set_debug_{

    # Set the debug compilation flags in $MakefileConf

    # Debug will be NewDebug after changes
    $Debug = $NewDebug;

    my $DEBUG; $DEBUG = '${DEBUGFLAG}' if $Debug eq "yes";
    print "Setting debugging flags to '$Debug' in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    s/^(\s*DEBUG\s*=).*/$1 $DEBUG/;
	    print;
	}
    }
}

##############################################################################

sub set_optimization_{

    # Set the optimization flags in $MakefileConf
    $Optimize = $NewOptimize;

    my $Level=$Optimize; $Level =~ s/-O//;
    print "Setting maximum optimization flag to $Optimize in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    if (/^\s*OPT([0-4])\s*=\s*/){
		if($1 > $Level){
		    $_ = "OPT$1 = -O$Level\n";
		}else{
		    $_ = "OPT$1 = -O$1\n";
		}
	    }
	    print;
	}
    }
}

##############################################################################

sub shell_command{

    my $command = join(' ',@_);
    print "$command\n" if $Verbose;

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
sub print_help_{

    print 
#BOC
"config.pl can be used for installing and setting various options for SWMF
or its components. The core of the script is in share/Scripts/config.pl,
and this is used by the config.pl scripts in the main SWMF and component 
directories. This help describes the options/features of the core script.
Additional features (if any) will be shown below.

This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.

Usage: config.pl [-help] [-verbose] [-show] [-dryrun] 
                 [-install[=s|=c] [-compiler=COMP] [-mpi=VERSION]] [-uninstall]
                 [-single|-double] [-debug|-nodebug] [-O0|-O1|-O2|-O3|-O4]

If called without arguments, the current settings are shown.

Information:

-h  -help      show help message.
-dryrun        dry run (do not modify anything, just show actions).
-show          show current settings.
-verbose       show verbose information.

(Un/Re)installation:

-uninstall     uninstall code (make distclean)

-install=c     (re)install code as an SWMF component (c)
-install=s     (re)install code as a stand-alone (s) code
-install       install code as a stand-alone if it is not yet installed,
               or reinstall the same way as it was installed originally:
               (re)creates Makefile.conf, Makefile.def, make install

-compiler=COMP copy Makefile.conf for a non-default F90 compiler COMP
-mpi=VERSION   copy share/include/mpif90_OS_VERSION 
               into share/Library/src/mpif90.h

Compilation:

-single        set precision to single in Makefile.conf and make clean
-double        set precision to double in Makefile.conf and make clean

-debug         select debug options for the compiler in Makefile.conf
-nodebug       do not use debug options for the compiler in Makefile.conf
-O0            set all optimization levels to -O0
-O1            set optimization levels to at most -O1
-O2            set optimization levels to at most -O2
-O3            set optimization levels to at most -O3
-O4            set maximum optimization level

Examples of use:

Show current settings: 

    config.pl

Show current settings with more detail: 

    config.pl -show

Install code with the ifort compiler and Altix MPI and select single precision:

    config.pl -install -compiler=ifort -mpi=Altix -single

Set optimization level to -O0 and switch on debugging flags:

    config.pl -debug -O0

Set optimization level to -03 and switch off debugging flags:

    config.pl -nodebug -O3

Uninstall code (if this fails, run config.pl -install first):

    config.pl -uninstall"
#EOC
    ,"\n\n";
}

##############################################################################

1;
