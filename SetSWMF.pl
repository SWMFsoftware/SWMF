#!/usr/bin/perl -i -w
use strict;

# Pattern to match component ID-s
my $ValidComp = 'IH|GM|IE|IM|UA';

&print_help if not @ARGV;

# Directories and files used
my $MakefileConf     = 'Makefile.conf';
my $MakefileConfOrig = 'share/build/Makefile.';
my $MakefileComp     = 'Makefile.def';
my $MakefileCompOrig = 'CON/Makefile.def';
my $MakefileCompCvs  = 'Makefile.def.orig';
my $GridSizeScript   = 'GridSize.pl';

# Default precision for installation
my $DefaultPrecision = 'double';

# Warning and error messages
my $WARNING = "!!! SetSWMF_WARNING:";
my $ERROR   = "!!! SetSWMF_ERROR:";

# Global variables for the settings
my $Installed;             # true if SWMF is installed ($MakefileConf exists)
my $OS='unknown';          # operating system in $MakefileComp
my $DIR='unknown';         # main directory for SWMF in $MakefileComp
my $Compiler='unknown';    # F90 compiler in $MakefileConf
my $Precision='unknown';   # Precision set in $MakefileConf
my @Version;               # Component versions selected in $MakefileComp
@Version = ('unknown');
my %Version;               # Hash table to get version for a component

# Obtain current settings
&get_settings;

# Default values for the various actions
my $Install;
my $Show;
my $NewPrecision;
my $DryRun;
my @NewVersion;
my $GridSize;
my $Uninstall;
my $IsCompilerSet;
my $Quiet;
my $Debug;

# Set actions based on the switches
my @switch = @ARGV;
foreach (@switch){
    if(/^-c=(.*)/)            {$Compiler=$1; $IsCompilerSet=1;  next};
    if(/^-d(ry)?(run)?$/)     {$DryRun=1;                       next};
    if(/^-h(elp)?$/i)         {&print_help};
    if(/^-p=(single|double)$/){$NewPrecision=$1;                next};
    if(/^-i(nstall)?$/)       {$Install=1;                      next};
    if(/^-s(how)?$/)          {$Show=1;                         next};
    if(/^-q(uiet)?$/)         {$Quiet=1;                        next};
    if(/^-D(ebug)?$/)         {$Debug=1;                        next};
    if(/^-v=(.*)/)            {push(@NewVersion,split(/,/,$1)); next};
    if(/^-uninstall$/)        {$Uninstall=1;                    next};
    if(/^-g=(.*)/)            {$GridSize.="$1,";                next};

    print "$WARNING Unknown switch: $_\n";
}

if($Uninstall){
    if(not $Installed){
	print "$WARNING SWMF is not installed.\n";
    }else{
	&shell_command("make distclean");
    }
    exit 0;
}

# Execute the actions in the appropriate order
&install_swmf if $Install;

if(not $Installed){
    print "$WARNING SWMF is not installed.\n";
    exit 0;
}

&set_precision;

&set_versions if @NewVersion;

&set_grid_size if $GridSize;

if($Show){
    &get_settings;
    &show_settings;
}

exit 0;
##############################################################################

sub get_settings{

    $Installed = -e $MakefileConf;

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = 'single';
    @Version = ();

    # Read information from $MakefileComp
    open(MAKEFILE, $MakefileComp)
	or die "$ERROR could not open $MakefileComp\n";

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
	$Precision = 'double' if /^\s*PRECISION\s*=\s*(\-r8|\-real_size\s*64)/;
    }

    close(MAKEFILE);
}

##############################################################################

sub show_settings{

    my $Show = 'SetSWMF.pl ';

    if(not $Installed){
	print "SWMF is not installed\n";
	return;
    }

    print "\n";
    print "The SWMF is installed in directory $DIR.\n";
    print "The installation is for the $OS operating system.\n";
    print "The selected F90 compiler is $Compiler.\n";
    print "The default precision for reals is $Precision precision.\n";
    print "The selected component versions and grid sizes are:\n";
    my $Comp;
    foreach $Comp (sort keys %Version){
	print "\t$Comp/$Version{$Comp}\t\t";
	if(-x "$Comp/$Version{$Comp}/$GridSizeScript"){
	    my $Grid = `cd $Comp/$Version{$Comp}; ./$GridSizeScript`;
	    $Grid =~ s/$GridSizeScript(\s*-g=)?//;
	    print "grid: $Grid" if $Grid;
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
	return;
    }

    # Obtain $OS and $DIR
    $OS  = `uname`    or die "$ERROR Could not obtain OS\n";
    chomp $OS;
    $DIR = `/bin/pwd` or die "$ERROR Could not obtain DIR\n";
    chomp $DIR;

    print "Creating $MakefileComp with OS=$OS, SWMF_ROOT=$DIR\n";
    if(not $DryRun){
	open(MAKEFILE,">$MakefileComp")
	    or die "$ERROR Could not open $MakefileComp for writing\n";

	print MAKEFILE "OS = $OS\nSWMF_ROOT = $DIR\n";
	close(MAKEFILE);
    }

    &shell_command("cat $MakefileCompOrig >> $MakefileComp");

    # Create $MakefileConf
    my $makefile = $MakefileConfOrig.$OS;
    $makefile .= $Compiler if $IsCompilerSet;
    die("$ERROR $makefile does not exist.") unless -e $makefile;
    &shell_command("cp $makefile $MakefileConf");

    # Read version and other info from main Makefile.COMP
    &get_settings;

    # set Makefile.COMP in versions to point to SWMF
    &set_version_makefile_comp;

    # Initialize CON and the components
    if($IsCompilerSet){
	&shell_command("make install COMPILER='$Compiler'");
    }else{
	&shell_command("make install");
    }

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

    if(not $change){
	&set_version_makefile_comp; return;
    }

    print "Modifying versions in $MakefileComp\n";
    return if $DryRun;

    # Set component versions in $MakefileComp
    @ARGV = ($MakefileComp);
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
	die "$ERROR could not find version $compversion in $MakefileComp\n"
	    unless $Found{$compversion}
    }

    # Create IH/BATSRUS if needed
    &shell_command("make IHBATSRUS")
	if $Version{"IH"} eq "BATSRUS" and not -f "IH/BATSRUS/src/Makefile";

    @Version = @NewVersion;

    &set_version_makefile_comp;

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

sub set_version_makefile_comp{

    # Set Makefile.COMP (if it exists) in the share/ directories
    # of the selected versions to include the global SWMF Makefile.COMP. 

    my $Version;
    my $QuietOrig = $Quiet; $Quiet = not $Debug;
    foreach $Version (sort @Version){
	my $file = "$Version/$MakefileComp";
	next unless -f $file;
	&shell_command("mv $Version/$MakefileComp $Version/$MakefileCompCvs")
	    unless -f "$Version/$MakefileCompCvs";
	&shell_command("echo include $DIR/$MakefileComp > $file");
	$file    = "$Version/$MakefileConf";
	&shell_command("echo include $DIR/$MakefileConf > $file");
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

sub print_help{

    print "
SetSWMF.pl can be used for installing and setting various options for SWMF.
This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.
This script will also be used by the GUI to interact with SWMF.

Usage: SetSWMF.pl [-h] [-q] [-D] [-d] [-i [-c=COMPILER]] [-p=PRECISION] [-s]
                  [-v=VERSION[,VERSION2,...] [-g=ID:GRIDSIZE[,ID:GRIDSIZE,...]
                  [-uninstall]

Options:

-c=COMPILER    create Makefile.conf with a non-default F90 compiler COMPILER

-d  -dry       dry run (do not modify anything, just show actions)

-D  -Debug     debug information shown

-g=ID:GRIDSIZE set the size of the grid to GRIDSIZE for the component 
               identified with ID. This flag can occur multiple times and/or 
               multiple grid sizes can be given in a comma separated list.

-h  -help      show this help message

-i  -install   install SWMF (create Makefile.conf, Makefile.COMP, make install)

-p=PRECISION   set precision (in Makefile.conf, copy mpif90.h and make clean)
               Possible values are 'single' and 'double'.

-q  -quiet     quiet execution

-s  -show      show current settings

-uninstall     uninstall SWMF (make distclean)

-v=VERSION     select component verion VERSION. This flag can occur 
               multiple times and/or multiple versions can be given
               in a comma separated list.
               If IH/BATSRUS is selected for the first time, 
               make IHBATSRUS is done.

Examples of use:

Show current settings: 

    SetSWMF.pl -s

Install SWMF with the pgf90 compiler and select single precision:

    SetSWMF.pl -i -c=pgf90 -p=single

Select IH/BATSRUS, IM/Empty and UA/Empty component versions:

    SetSWMF.pl -v=IH/BATSRUS -v=IM/Empty,UA/Empty

Set the grid size for GM and IH:

    SetSWMF.pl -g=GM:8,8,8,400,100 -g=IH:6,6,6,800,1

";
    exit 0;
}

##############################################################################
