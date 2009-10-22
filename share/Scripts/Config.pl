#!/usr/bin/perl -i
use strict;

# Default compiler per machine or OS
my %Compiler = (
		"Linux"     => "f95",
		"Darwin"    => "f95",
		"OSF1"      => "f90",
		"IRIX64"    => "f90",
		"AIX"       => "xlf90",
		"palm"      => "ifort",
		"cfe"       => "ifort",
		"pfe"       => "ifort",
		"sysx"      => "xlf90",
		"nyx-login" => "pgf90",
		"hera"      => "mpiifort",
		);

# Default MPI library per machine or OS
my %MpiVersion = (
		  "Linux"     => "mpich",
		  "Darwin"    => "mpich",
		  "OSF1"      => "mpich",
		  "AIX"       => "IBM",
		  "IRIX64"    => "SGI",
		  "palm"      => "Altix",
		  "cfe"       => "Altix",
		  "pfe"       => "Altix",
		  "grendel"   => "openmpi",
		  "nyx-login" => "openmpi",
		  "hera"      => "mvapich",
		  );

my $WARNING_='share/Scripts/Config.pl WARNING:';
my $ERROR_  ='share/Scripts/Config.pl ERROR:';

# Obtain $OS, $DIR, and the machine name and provide it to caller script
our $OS  = `uname`    or die "$ERROR_ could not obtain OS\n"; chop $OS;
our $DIR = `/bin/pwd` or die "$ERROR_ could not obtain DIR\n"; chop $DIR;
our $Machine = `hostname -s`; chop $Machine; 

# remove numbers from the machine name
$Machine =~ s/\d+$//; 

# These are either obtained from the calling script or set here
our $Component;             # The SWMF component the code is representing
our $Code;                  # The name of the code
($Component, $Code) = ($DIR =~ /([A-Z][A-Z])\/([^\/]+)$/)
    unless $Code;

# Strings for the error and warning messages for the caller script
our $ERROR   = "$Code/Config.pl ERROR:";
our $WARNING = "$Code/Config.pl WARNING:";

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
our $MakefileConfOrig = 'share/build/Makefile';
our $MpiHeader        = 'share/Library/src/mpif90.h';
our $MpiHeaderOrig    = 'share/include/mpif90_';
our $MakefileRules    = 'Makefile.RULES';
our $MakefileDepend   = 'Makefile.DEPEND';

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
my $NewMpi;
my $IsCompilerSet;
my $Debug;
my $Mpi;
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
    if(/^-mpi$/i)             {$NewMpi="yes";                   next};
    if(/^-nompi$/i)           {$NewMpi="no";                    next};
    if(/^-O[0-4]$/i)          {$NewOptimize=$_;                 next};  
    if(/^-g(rid)?$/)          {$ShowGridSize=1;                 next};
    if(/^-g(rid)?=([\d,]+)$/) {$NewGridSize=$+;                 next};

    if(/^-g(rid)?=(.*)$/ and $IsComponent){
	die "$ERROR: incorrect grid size -g=$+\n"};

    $Remaining{$_}=1;
}

if(not $MakefileDefOrig and not $IsComponent){
    warn "$WARNING $Code cannot be used in stand alone mode!\n ".
	"Switching to component mode...\n";
    $IsComponent = 1;
}

&print_help_ if $Help;

if($Uninstall){
    if(not $Installed){
	warn "$ERROR_ $Code is not installed.\n";
	exit 1;
    }else{
	&shell_command("cd share; make distclean")
	    if -d "share" and not $IsComponent;
	&shell_command("cd util; make distclean")
	    if -d "util" and not $IsComponent;
	&shell_command("make allclean");
	&shell_command("rm -f Makefile.def Makefile.conf data ".
		       "src*/$MakefileDepend src*/$MakefileRules");
	exit 0;
    }
}

# Execute the actions in the appropriate order
&install_code_ if $Install;

if(-f $MakefileDef and not $IsComponent){
    my @Stat = stat $MakefileDef;
    my $Time = $Stat[9];
    my @Stat = stat $MakefileDefOrig;
    my $TimeOrig = $Stat[9];
    die "$ERROR $MakefileDefOrig is newer than $MakefileDef !\n".
	"   Reinstall or merge changes into $MakefileDef !\n"
	if $Time < $TimeOrig;
}

# Change precision of reals if required
if($NewPrecision and $NewPrecision ne $Precision){
    &shell_command("make clean");
    &set_precision_;
}

# Change debugging flags if required
&set_debug_ if $NewDebug and $NewDebug ne $Debug;

# Link with MPI vs. NOMPI library if required
&set_mpi_ if $NewMpi and $NewMpi ne $Mpi;

# Change optimization level if required
&set_optimization_ if $NewOptimize and $NewOptimize ne $Optimize;

if($Show){
    &get_settings_;
    &show_settings_;
}

# Recreate Makefile.RULES with the current settings
&create_makefile_rules;

##############################################################################
sub get_settings_{

    $Installed = (-e $MakefileConf and -e $MakefileDef);

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = "unknown";

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
    $Mpi   = "yes";
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

	  $Precision = lc($1) if /^\s*PRECISION\s*=.*(SINGLE|DOUBLE)PREC/;
          $Debug = "yes" if /^\s*DEBUG\s*=\s*\$\{DEBUGFLAG\}/;
	  $Mpi   = "no"  if /^\s*MPILIB\s*=.*\-lNOMPI/;
          $Optimize = $1 if /^\s*OPT[0-4]\s*=\s*(-O[0-4])/;
      }
  }
    close(MAKEFILE);

    open(MPIHEADER, $MpiHeader) or open(MPIHEADER, "../../$MpiHeader") or
	return;

    my $IsFound = 0;
    while(<MPIHEADER>){
	next unless /MPI_HEADER_FILE\s*=.*_(\w+)\.h/;
	$IsFound = 1;
	$MpiVersion = $1;
	last;
    }
    close(MPIHEADER);
    warn "$WARNING_ could not find MPI_HEADER_FILE string in $MpiHeader\n"
	unless $IsFound;
}

##############################################################################

sub show_settings_{

    if(not $Installed){
	print "$Code is not installed\n";
	exit 0;
    }

    print "\n";
    print "$Code is installed in directory $DIR\n";

    if($IsComponent){
	print "    as the $Component component.\n";
    }else{
	print "    as a stand-alone code.\n";
    }
    print "The installation is for the $OS operating system.\n";
    print "The selected F90 compiler is $Compiler.\n";
    print "The selected MPI library is  $MpiVersion.\n";
    print "The default precision for reals is $Precision precision.\n";
    print "The maximum optimization level is $Optimize\n";
    print "Debugging flags: $Debug\n";
    print "Linked with MPI: $Mpi\n";

    print "\n";

}

##############################################################################
sub install_code_{

    my $Text = $Installed ? "Reinstalling $Code" : "Installing $Code";
    $Text .= " as a $Component component" if $IsComponent;  
    print "$Text\n";

    if($IsComponent){
	my $dir = $DIR; $dir =~ s|/[^/]*/[^/]*$||;  # go two directories up
	my $makefile = "$dir/Makefile.def";          # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > Makefile.def");

	$makefile = "$dir/Makefile.conf"; # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > Makefile.conf");
    }else{
	die "$ERROR_ original $MakefileDef is not given\n" unless
	    $MakefileDefOrig;
	die "$ERROR_ $MakefileDefOrig is missing\n" unless
	    -f $MakefileDefOrig;
	&shell_command("echo OS=$OS > $MakefileDef");
	&shell_command("echo ${Component}DIR=$DIR >> $MakefileDef");
	&shell_command("echo COMPILER=$Compiler >> $MakefileDef");
	&shell_command("cat $MakefileDefOrig >> $MakefileDef");

	my $Makefile = "$MakefileConfOrig.$OS.$Compiler";
	if(-f $Makefile){
	    &shell_command("cat $Makefile > $MakefileConf");
	}else{
	    # Try to use generic Makefile with provided compiler
	    warn "$WARNING_: $Makefile was not found,".
		" using generic $MakefileConfOrig.conf\n";
	    $Makefile = "$MakefileConfOrig.conf";
	    open(IN, $Makefile) or die "$ERROR_ $Makefile is missing\n";
	    open(OUT, ">$MakefileConf") 
		or die "$ERROR_ could not open $MakefileConf\n";
	    while(<IN>){
		s/_COMPILER_/$Compiler/;
		s/_OS_/$OS/;
		print OUT $_;
	    }
	    close IN; close OUT;
	}
	
	my $Header = "$MpiHeaderOrig${OS}_$MpiVersion.h";
	if(-f $Header){
	    &shell_command("cat $Header > $MpiHeader");
	}else{
	    warn "$WARNING_: $Header was not found,".
		" using generic ${MpiHeaderOrig}mpich.h\n";
	    &shell_command("cat ${MpiHeaderOrig}mpich.h > $MpiHeader");
	}
    }

    # Read info from main Makefile.def
    &get_settings_;

    # Set initial precision for reals
    $NewPrecision = $DefaultPrecision unless $NewPrecision;
    &set_precision_ if $NewPrecision ne $Precision;

    # Create Makefile.RULES as needed
    &create_makefile_rules;

    # Link data directory if possible
    &link_smwf_data;

    # Install the code
    &shell_command("cd share; make install") 
	if -d "share" and not $IsComponent;
    &shell_command("cd util; make install") 
	if -d "util" and not $IsComponent;
    &shell_command("make install");

    # Now code is installed
    $Installed = 1 unless $DryRun;
}

##############################################################################

sub set_precision_{

    # Set the precision for reals in $MakefileConf

    # Precision will be NewPrecision after changes
    $Precision = $NewPrecision;

    my $PREC = '${'.uc($Precision).'PREC}';
    print "Setting PRECISION variable to $PREC in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    s/^(\s*PRECISION\s*=\s*).*/$1$PREC/;
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

sub set_mpi_{

    # Select the MPI or NOMPI library in $MakefileConf

    # $Mpi will be $NewMpi after changes
    $Mpi = $NewMpi;

    print "Selecting MPI library in $MakefileConf\n" if $Mpi eq "yes";
    print "Selecting NOMPI library in $MakefileConf\n" if $Mpi eq "no";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    # Comment/uncomment MPILIB definitions
	    if(/MPILIB\s*=/){
		s/^\s*M/\#M/ if /lNOMPI/ eq ($Mpi eq "yes");
		s/^\#\s*M/M/ if /lNOMPI/ eq ($Mpi eq "no");
	    }
	    # Modify LINK.f90 definition
	    if(/^\s*LINK.f90\s*=.*mpif90/){
		s/\{CUSTOMPATH_MPI\}/\{COMPILE.f90\}\# \t/ if $Mpi eq "no";
		s/\{COMPILE.f90\}\#\s*/\{CUSTOMPATH_MPI\}/ if $Mpi eq "yes";
	    }
	    print;
	}
    }
    &shell_command("make NOMPI") if $Mpi eq "no";
    print "Remove executable and make it to link with the (NO)MPI library!\n";
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

sub create_makefile_rules{

    my @InFile = glob("src*/$MakefileRules.all */*/src*/$MakefileRules.all");

    return unless @InFile;

    # Hash for general configuration settings
    my %Settings = (OS         => $OS, 
		    Compiler   => $Compiler, 
		    MpiVersion => $MpiVersion, 
		    Debug      => $Debug,
		    Machine    => $Machine,
		    Precision  => $Precision);

    # Add settings from the caller Config.pl script
    my $Settings = shift;
    $Settings{$1}=$2 while $Settings =~ s/(\w+)\s*=\s*([^,\n]+)//;

    # Create Makefile.RULES from Makefile.RULES.all in all src* directories
    my $InFile;
    foreach $InFile (@InFile){

	$InFile =~ /([\w,\/]*)\//;
	my $SrcDir = $1;

	# Open Makefile.RULES.all for reading
	open(INFILE, $InFile) or die "$ERROR_: could not open $InFile\n";

	# Open Makefile.RULES for writing
	my $OutFile = "$SrcDir/$MakefileRules";
	open(OUTFILE, ">$OutFile") or die "$ERROR_: could not open $OutFile\n";

	# Evaluate conditional rules in Makefile.RULES.all
	my $Condition;
	while(<INFILE>){
	    next if /^\#/ or /^\s/; 
	    $Condition = $_;

	    # Replace $xxx variables with their actually set values
	    my $key;
	    foreach $key (keys %Settings){
		$Condition =~ s/\$$key/"$Settings{$key}"/g;
	    }

	    # Skip rules unless condition is true
	    next unless eval($Condition);

	    # Create compilation rules
	    my $Rule;
	    while($Rule = <INFILE>){
		last unless $Rule =~ s/^\t//;

		$Rule =~ /([\w\.]+)\s*$/;

		my $SrcFile = $1;
		my $ObjectFile = $SrcFile; $ObjectFile =~ s/\.\w+$/.o/;

		print OUTFILE "$ObjectFile: $SrcFile\n\t$Rule\n";
	    }
	}
	close INFILE;
	close OUTFILE;
    }
}

##############################################################################

sub link_smwf_data{

    return if -d "data";
    my $SwmfData = "";
    foreach ("$ENV{HOME}/SWMF_data", "/csem1/SWMF_data"){
	next unless -d $_;
	$SwmfData = $_;
	last;
    }
    my $DataDir;
    if($IsComponent){
	$DataDir = "$SwmfData/$Component/$Code/data";
    }else{
	$DataDir = "$SwmfData/data";
    }
    &shell_command("ln -s $DataDir data") if -d $DataDir;
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
#!QUOTE: \subsection{Installation and Configuration with Config.pl}
#!ROUTINE: Config.pl - (un)installation and configuration of SWMF/components
#!DESCRIPTION:
# The Config.pl provides a single uniform interface towards 
# installation, configuration and uninstallation for the SWMF and its
# components.
#
#!REVISION HISTORY:
# 12/16/2006 G. Toth - initial version based on SetSWMF.pl
#EOP
sub print_help_{

    print 
#BOC
"Config.pl can be used for installing and setting various options for SWMF
or its components. The core of the script is in share/Scripts/Config.pl,
and this is used by the Config.pl scripts in the main SWMF and component 
directories. This help message starts with the options/features/examples 
of the core script, and continues with the additional features (if any)
of the SWMF/component script (starting with the text 'Additional ...').

This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.

Usage: Config.pl [-help] [-verbose] [-show] [-dryrun] 
                 [-install[=s|=c] [-compiler=COMP] [-mpi=VERSION]] [-uninstall]
                 [-single|-double] [-debug|-nodebug] [-mpi|-nompi]
                 [-O0|-O1|-O2|-O3|-O4]

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
               Only works together with -install flag

-mpi=VERSION   copy share/include/mpif90_OS_VERSION 
               into share/Library/src/mpif90.h
               Only works together with -install flag

Compilation:

-single        set precision to single in Makefile.conf and make clean
-double        set precision to double in Makefile.conf and make clean

-debug         select debug options for the compiler in Makefile.conf
-nodebug       do not use debug options for the compiler in Makefile.conf
-nompi         compile and link with the NOMPI library for serial execution
-mpi           compile and link with the MPI library for parallel execution
-O0            set all optimization levels to -O0
-O1            set optimization levels to at most -O1
-O2            set optimization levels to at most -O2
-O3            set optimization levels to at most -O3
-O4            set maximum optimization level

Examples of use:

Show current settings: 

    Config.pl

Show current settings with more detail: 

    Config.pl -show

Install code with the g95 compiler and Intel MPI and select single precision:

    Config.pl -install -compiler=g95 -mpi=Intel -single

Set optimization level to -O0, switch on debugging flags and link with NOMPI:

    Config.pl -debug -O0 -nompi

Set optimization level to -03, switch off debugging flags and link with MPI:

    Config.pl -nodebug -O3 -mpi

Uninstall code (if this fails, run Config.pl -install first):

    Config.pl -uninstall"
#EOC
    ,"\n\n";
}

##############################################################################

1;
