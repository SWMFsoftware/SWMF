#!/usr/bin/perl -s

my $Help        = ($h or $H or $help);
my $Mode        = ($m or $mode or "auto");
my $OutputOnly  = ($o or $output);
my $InputOnly   = ($i or $input);
my $CheckOnly   = ($c or $check);
my $Verbose     = ($v or $verbose);
my $TimeUnit    = ($t or $u or $timeunit or $unit);
my $Repeat      = ($r or $repeat);
my $Wait        = ($w or $wait or 120);
my $RestartTree = $ARGV[0];        # Name of the restart tree directory
$RestartTree =~ s/\/+$//;          # Remove trailing /

use strict;

&print_help if $Help;

my $WARNING ="WARNING in Restart.pl:";           # Warning message string
my $ERROR   ="ERROR in Restart.pl:";             # Error message string
my $HELP    ="\nType Restart.pl -h for help.\n"; # Help message string


# Figure out if we are in a framework mode or not
my $Framework = 0;

if($Mode =~ /^a/i){
    $Framework = 1 if glob('SWMF*.exe') or -d 'STDOUT' or -f 'RESTART.in';
}elsif($Mode =~ /^f/){
    $Framework = 1;
}elsif($Mode !~ /^s/){
    die "$ERROR invalid mode (should be 'a', 'f' or 's'): -m=$Mode\n";
}

# Check for illegal combination of switches
die "$ERROR at most one argument can be specified!$HELP" if $#ARGV > 0;
die "$ERROR cannot use -i and -o together!$HELP" if $InputOnly and $OutputOnly;
die "$ERROR cannot use -i and -r together!$HELP" if $InputOnly and $Repeat;
die "$ERROR cannot use -c and -r together!$HELP" if $CheckOnly and $Repeat;
die "$ERROR restart tree must be specified with -i option!$HELP" 
    if $InputOnly and not $RestartTree;
die "$ERROR restart tree cannot be specified with -r option!$HELP" 
    if $Repeat and $RestartTree;

# Declare global variables
my $RestartOutFile = "RESTART.out";# Name of SWMF output restart file
my $RestartInFile  = "RESTART.in"; # Name of SWMF input restart file
my $SimulationTime = -1;           # Simulation time
my $nStep          = -1;           # Number of steps

# List of input and output restart directory name(s) for each component
# Alternative names should be separated by commas without space.
my %RestartOutDir = (
		     GM => "GM/restartOUT",
                     LC => "LC/restartOUT", 
                     SC => "SC/restartOUT",
		     IH => "IH/restartOUT",
                     OH => "OH/restartOUT",
		     IM => "IM/restartOUT",
		     PW => "PW/restartOUT",
		     RB => "RB/restartOUT",
		     UA => "UA/restartOUT,UA/RestartOUT" );

my %RestartInDir =  (
		     GM => "GM/restartIN",
                     LC => "LC/restartIN",
		     SC => "SC/restartIN",
		     IH => "IH/restartIN",
                     OH => "OH/restartIN",
		     IM => "IM/restartIN",
		     PW => "PW/restartIN",
		     RB => "RB/restartIN",
		     UA => "UA/restartIN,UA/RestartIN" );

# Hashes for the actually found directories
my %RestartOutDirFound;
my %RestartInDirFound;

# The name of the restart header file (if any) for each component.
my %HeaderFile   =  (
		     GM => "restart.H",
                     LC => "restart.H",
		     SC => "restart.H",
                     OH => "restart.H",
		     IH => "restart.H" );

# List possible time units and corresponding number of seconds
my %UnitSecond = ("ns" => 1e-9,      # nano second
                  "us" => 1e-6,      # micro second
		  "ms" => 0.001,     # millisecond
		  "s" => 1,          # second
		  "m" => 60,         # minute
		  "h" => 3600,       # hour
		  "d" => 86400,      # day
		  "y" => 31536000    # year
		  );

# Check the time unit parameter if given
die "$ERROR time unit $TimeUnit is unknown!\n" 
    if $TimeUnit and not $UnitSecond{$TimeUnit};

# Figure out the standalone component if not in framework mode
if($Framework){
    print "Restart.pl running in framework mode\n" if $Verbose;
}else{
    my $Comp;
    foreach (keys %HeaderFile){
	$Comp = $_;
	last if -d $Comp;
    }
    die "$ERROR could not find a usable component directory\n" unless -d $Comp;

    # Use this file name to check component output restart files
    $RestartOutFile = "$RestartOutDir{$Comp}/$HeaderFile{$Comp}";

    # Use this file name to check files in the input restart tree
    $RestartInFile  = "$Comp/$HeaderFile{$Comp}";
    print "Restart.pl running in standalone mode for component $Comp\n" 
	if $Verbose;
}

LOOP:{
    if($Repeat){
	if(not -f $RestartOutFile){
	    # If there is no new file wait $Repeat seconds
	    print "sleep $Repeat\n" if $Verbose;
	    sleep $Repeat;
	    redo LOOP;
	}
	# Check if the output restart files are old enogh to be moved
	my @stat = stat($RestartOutFile);
	my $age = time - $stat[9];
	if($age < $Wait){
	    my $wait = $Wait - $age;
	    print "sleep $wait\n" if $Verbose;
	    sleep $Wait - $age;
	}
    }

    # Initialize simulation time and number of steps to impossible values
    $SimulationTime = -1;
    $nStep          = -1;

    # Create restart tree if required
    if(not $InputOnly){
	&create_tree_check;
	&create_tree unless $CheckOnly;
    }

    # Link restart tree if required
    if(not $OutputOnly){
	&link_tree_check;
	&link_tree unless $CheckOnly;
    }

    # Loop if required
    redo LOOP if $Repeat;
}

exit 0;

##############################################################################
sub get_time_step{
    my $File = shift;

    my $Time = -1;
    my $Step = -1;
    open(FILE, $File) or die "$ERROR could not open file $File\n";
    while(<FILE>){
	if(/\#TIMESIMULATION/){
	    # Read in simulation time
	    $Time = <FILE>; chop($Time);
	    $Time =~ s/^\s+//; # Remove leading spaces
	    $Time =~ s/\s.*//; # Remove anything after a space
	    $Time += 0;        # Convert to a number
	}
	if(/\#NSTEP/){
	    # Read in number of steps
	    $Step = <FILE>; chop($Step);
	    $Step =~ s/^\s+//; # Remove leading spaces
	    $Step =~ s/\s.*//; # Remove anything after a space
	    $Step += 0;        # Convert to a number
	}
    }
    die "$ERROR could not find simulation time in $File!\n" if $Time < 0;
    die "$ERROR could not find time step in $File!\n" if $Step < 0;

    print "# Restart.pl read Time=$Time Step=$Step from $File\n" if $Verbose;

    # Save time and step if not yet specified
    $SimulationTime = $Time if $SimulationTime < 0;
    $nStep          = $Step if $nStep < 0;

    # Check if times are consistent in a time accurate run
    die "$ERROR in file $File time $Time differs from ".
	"simulation time $SimulationTime!\n" 
	if $SimulationTime and abs($Time - $SimulationTime) > 0.01;
}

##############################################################################
sub create_tree_check{

    # Check the SWMF restart file
    die "$ERROR could not find restart file $RestartOutFile!\n" 
	unless -f $RestartOutFile;

    # Obtain time/step from the restart file
    &get_time_step($RestartOutFile);

    # Set the name of restart tree if not specified in the command line
    if(not $ARGV[0]){
	# Check if it is a time accurate run
	if($SimulationTime){
	    # If the time unit is not set try to guess it from simulation time
	    if(not $TimeUnit){
		my $Unit;
		$TimeUnit = "ns"; 
		foreach $Unit (sort {$UnitSecond{$a} <=> $UnitSecond{$b}} 
			       keys %UnitSecond){
		    $TimeUnit = $Unit if $SimulationTime >= $UnitSecond{$Unit};
		}
	    }
	    # Use the simulation time for time accurate runs
	    $RestartTree = sprintf("RESTART_t%6.2f%s", 
				   $SimulationTime/$UnitSecond{$TimeUnit},
				   $TimeUnit);
	}else{
	    # Use the time step number for steady state runs
	    $RestartTree = sprintf "RESTART_n%6d", $nStep;
	}

	# Replace spaces with zeros
	$RestartTree =~ s/ /0/g;

	print "# Restart.pl set restart tree name to $RestartTree/.\n"
	    if $Verbose;
    }

    # Check the restart tree directory
    die "$ERROR restart tree $RestartTree is in the way!\n" if -d $RestartTree;

    # Check output restart directories for all components
    my $Comp;
  COMPONENT:
    foreach $Comp (sort keys %RestartOutDir){
	next COMPONENT unless -d $Comp;

	my $Dirs = $RestartOutDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir};
	if(not -d $Dir){
	    print "$WARNING could not find directory $Dirs!\n";
	    next COMPONENT;
	}

	opendir(DIR,$Dir) or die "$ERROR could not open directory $Dir!\n";
	my @Content = readdir(DIR);
	closedir(DIR);
	if($#Content < 2){
	    print "$WARNING directory $Dir is empty!\n";
	    next COMPONENT;
        }

	# Store the output directory for this component
	$RestartOutDirFound{$Comp} = $Dir;

	if($Framework){
	    # Check if header file exists and check the simulation time
	    my $HeaderFile = $HeaderFile{$Comp};
	    &get_time_step("$Dir/$HeaderFile") if $HeaderFile;
	}

	print "# Restart.pl has checked $Dir\n" if $Verbose;
    }

    print "# Restart.pl has checked output restart file and directories.\n";
}
##############################################################################
sub create_tree{

    # Create restart directory
    print "mkdir $RestartTree\n" if $Verbose;
    mkdir $RestartTree,0777 
	or die "$ERROR restart tree $RestartTree could not be created!\n";

    if($Framework){
	# Move the SWMF restart file
	my $File = "$RestartTree/$RestartOutFile";
	print "mv $RestartOutFile $File\n" if $Verbose;
	rename $RestartOutFile, $File or 
	    die "$ERROR could not move $RestartOutFile into $File!";
    }

    # Move the output restart directories of the components into the tree
    # and create empty output restart directories
    my $Comp;
    foreach $Comp (sort keys %RestartOutDirFound){
	my $Dir=$RestartOutDirFound{$Comp};

	print "mv $Dir $RestartTree/$Comp\n" if $Verbose;
	rename $Dir, "$RestartTree/$Comp" or 
	    die "$ERROR could not move $Dir into $RestartTree/$Comp!\n";

	print "mkdir $Dir\n" if $Verbose;
	mkdir $Dir, 0777 or die "$ERROR could not create directory $Dir!\n";
    }

    print "# Restart.pl has created restart tree $RestartTree/.\n";
}
##############################################################################
sub link_tree_check{

    # If the create phase was checked only the the tree is not created
    my $NoTreeCheck = ($CheckOnly and not $InputOnly);

    # Check the tree
    die "$ERROR restart tree $RestartTree is missing!\n" 
	unless (-d $RestartTree or $NoTreeCheck);

    if($Framework){
	# Check for an existing restart file
	die "$ERROR file $RestartInFile is in the way!\n" if 
	    (-f $RestartInFile and not -l $RestartInFile);
    }

    # Check the SWMF/component restart header file in the restart tree
    my $File;
    if($Framework){
	$File = "$RestartTree/$RestartOutFile";
	die "$ERROR could not find SWMF restart file $File!\n" 
	    unless (-f $File or $NoTreeCheck);
    }else{
	$File = "$RestartTree/$RestartInFile";
	die "$ERROR could not find restart header file $File!\n" 
	    unless (-f $File or $NoTreeCheck);
    }

    # Set the step and the simulation time
    &get_time_step($File) unless $NoTreeCheck;

    my $Comp;
  COMPONENT:
    foreach $Comp (sort keys %RestartInDir){
	next COMPONENT unless -d $Comp;
	if(not -d "$RestartTree/$Comp"){
	    print "$WARNING could not find restart directory ".
		"$RestartTree/$Comp!\n";
	    next COMPONENT;
	}
		
	my $Dirs = $RestartInDir{$Comp};
	my $Dir;
	if($Dirs =~ /,/){
	    # Figure out which restart in directory should be used
	    foreach (split /,/,$Dirs){$Dir=$_ if -d $_ or -l $_};
	    die "$ERROR could not find input restart directory/link $Dirs!\n" 
		unless $Dir;
	}else{
	    # There is only one possible directory name, use that one
	    $Dir = $Dirs;
	}

	# Store the input directory for this component
	$RestartInDirFound{$Comp} = $Dir;

	if($Framework){
	    # Check if the header file exists and check the simulation time
	    my $HeaderFile = $HeaderFile{$Comp};
	    &get_time_step("$RestartTree/$Comp/$HeaderFile") if $HeaderFile;
	}

	print "# Restart.pl has checked $Dir\n" if $Verbose;
    }

    print "# Restart.pl has checked  input restart file and directories.\n";
}
##############################################################################
sub link_tree{

    if($Framework){
	# Remove existing input restart link
	if(-l $RestartInFile){
	    print "rm -f $RestartInFile\n" if $Verbose;
	    unlink $RestartInFile or 
		die "$ERROR could not remove link $RestartInFile!\n";
	}

	# Link in the SWMF restart file in the restart tree
	my $File = "$RestartTree/$RestartOutFile";
	print "ln -s $File $RestartInFile\n" if $Verbose;
	symlink $File, $RestartInFile or 
	    die "$ERROR could not link $File to $RestartInFile!\n";
    }

    my $Comp;
    foreach $Comp (sort keys %RestartInDirFound){
	my $Dir=$RestartInDirFound{$Comp};

	# Remove input restart link or directory
	if(-l $Dir){
	    print "rm -f $Dir\n" if $Verbose;
	    unlink $Dir or die "$ERROR could not remove link $Dir!\n";
	}elsif(-d $Dir){
	    print "rmdir $Dir\n" if $Verbose;
	    rmdir $Dir or die "$ERROR could not remove directory $Dir!\n";
	}

	# Link input restart directory in the restart tree
	print "ln -s ../$RestartTree/$Comp $Dir\n" if $Verbose;
	symlink "../$RestartTree/$Comp", $Dir or 
	    die "$ERROR could not link $RestartTree/$Comp to $Dir!\n";
    }
    print "# Restart.pl has linked  restart tree $RestartTree/.\n";
}
##############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Restarting Runs with Restart.pl}
#!ROUTINE: Restart.pl - create and link restart file trees in SWMF
#!DESCRIPTION:
# This script is copied into the run directory and it should be executed there.
# The script can collect the output restart files from individual components 
# into one directory tree: the restart tree. The script can also link the 
# input restart directories of the components to the restart tree. 
# Finally the script can run continuously (in the background) and save
# the output restart files into restart trees as they are created during 
# the run. This way multiple restart trees can be saved from a single run.
# The script checks the consistency of the restart files.
#
#!REVISION HISTORY:
# 12/21/2004 G. Toth - initial version
# 02/10/2005 G. Toth - allow restart info for a subset of components: 
#                      warns about missing restart dirs but does not die.
#EOP

sub print_help{
    print
#BOC
'Purpose:

    Collect current output restart directories into one tree (and link to it).
    Link input restart files to an existing directory tree. 
    Create multiple restart trees while the SWMF is running.

Usage:

    Restart.pl -h

    Restart.pl [-o] [-t=s|m|h|d|y] [-m=a|f|s] [-c] [-v] [DIR]

    Restart.pl -i [-m=a|f|s] [-c] [-v] DIR

    Restart.pl -r=REPEAT [-w=WAIT] [-o] [-t=s|m|h|d|y] [-v] &

    -h -help    Print help message and exit.

    -o -output  Create restart tree from output directories but do not link.
                Cannot be used together with the -i switch.
                Default is to link input directories to the tree as well.

    -i -input   Link input restart directories to an existing restart tree.
                The name of the restart tree must be specified.
                Cannot be used together with the -o switch.
                Default is to create the restart tree first and then link.

    -c -check   Check but do not actually create or link.
                Default is to create and link as specified by -i and -o.

    -r=REPEAT   Repeat creating (and linking unless -o is used) of the 
    -repeat=... restart tree every REPEAT seconds. This can be used to
                store multiple copies of the restart tree.
                Cannot be used together with the -i or -c switches.
                Cannot specify the name of the directory tree.
                Default is to create the restart tree only once.

    -w=WAIT     Wait WAIT seconds before moving the output restart files.
    -wait=...   This switch can only be used with the -r switch.
                Default is to wait 2 minutes.

    -t=UNIT     Time unit to form the name of the restart tree from the
    -time=...   simulation time (only matters for time accurate run).
    -u=UNIT     The UNIT can be given as one of the following characters:
    -unit=...   s, m, h, d, y corresponding to seconds, minute, hour, day and
                year respectively. The -t option has no effect if the 
                name of the restart tree is specified by the parameter DIR.
                The default time unit is the largest unit which does not 
                exceed the simulation time.

    -m=MODE     Mode of operation can be "auto", "framework" or "standalone".
    -mode=...   Only the first character is significant. The default mode is
                "auto", when Restart.pl tries to decide if the run directory 
                was created for the framework or a standalone component.

    -v -verbose Print verbose information.

    DIR         Name of the restart directory tree. This argument
                must be specified if the -i switch is used. Otherwise
                the default name is RESTART_n012345 for steady state runs
                and RESTART_t012.34u for time accurate runs, where the
                numbers should be replaced with the actual time step and
                simulation time, and the "u" with the actual time unit.

Examples:

    Check the output and input restart files and directories:

Restart.pl -c

    Create restart tree from current results and link input to it:

Restart.pl

    Create restart trees every 10 minutes, use hours as the 
    time unit for the simulation time in the restart tree names,
    and save output and error messages (if any) into Restart.log:

Restart.pl -o -r=600 -t=h >& Restart.log &

    Check linking to the existing RESTART_t002.00h tree:

Restart.pl -i -c RESTART_t002.00h

    Link to the existing RESTART_t002.00h tree and print verbose info:

Restart.pl -v -i RESTART_t002.00h

    Use Restart.pl in standalone mode and wait 60 seconds before moving files

Restart.pl -m=standalone -w=60'
#EOC
    ,"\n\n";
    exit;
}
##############################################################################
