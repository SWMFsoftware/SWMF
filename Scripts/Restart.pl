#!/usr/bin/perl -s

my $Help       = $h or $H or $help;
my $OutputOnly = $o;
my $InputOnly  = $i;
my $CheckOnly  = $c;
my $Verbose    = $v;
my $TimeUnit   = $t;

use strict;

&print_help if $Help;

my $ERROR ="ERROR in Restart.pl:"; # Error message string

if($InputOnly and $OutputOnly){
    print "$ERROR cannot use -i and -o together!\n".
	"Type Restart.pl -h for help.\n";
    exit 1;
}

my $RestartTree = $ARGV[0];        # Name of the restart tree directory
$RestartTree =~ s/\/+$//;          # Remove trailing /

if($InputOnly and not $RestartTree){
    print "$ERROR restart tree must be specified with -i option!\n".
	"Type Restart.pl -h for help.\n";
    exit 1;
}

my $RestartOutFile = "RESTART.out";# Name of SWMF output restart file
my $RestartInFile  = "RESTART.in"; # Name of SWMF input restart file
my $SimulationTime = -1;           # Simulation time
my $nStep          = -1;           # Number of steps

# List of input and output restart directory name(s) for each component
# Multiple names should be separated by commas without space.
my %RestartOutDir = (
		     GM => "GM/restartOUT",
		     SC => "SC/restartOUT",
		     IH => "IH/restartOUT",
		     IM => "IE/restartOUT",
		     UA => "UA/RestartOUT,UA/restartOUT" );

my %RestartInDir =  (
		     GM => "GM/restartIN",
		     SC => "SC/restartIN",
		     IH => "IH/restartIN",
		     IM => "IE/restartIN",
		     UA => "UA/RestartIN,UA/restartIN" );

my $HeaderFile = "restart.H"; # the name of the restart header files

# List possible time units and corresponding number of seconds
my %UnitSecond = (
		  "s" => 1,          # second
		  "m" => 60,         # minute
		  "h" => 3600,       # hour
		  "d" => 86400,      # day
		  "y" => 31536000    # year
		  );

# Check the time unit parameter if given
die "$ERROR time unit $TimeUnit is unknown!\n" 
    if $TimeUnit and not $UnitSecond{$TimeUnit};

if(not $InputOnly){
    &create_tree_check;
    &create_tree unless $CheckOnly;
}
if(not $OutputOnly){
    &link_tree_check;
    &link_tree unless $CheckOnly;
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
	    $Time = <FILE>;
	    $Time =~ s/^\s+//; # Remove leading spaces
	    $Time =~ s/\s.*//; # Remove anything after a space
	}
	if(/\#NSTEP/){
	    # Read in number of steps
	    $Step = <FILE>;
	    $Step =~ s/^\s+//; # Remove leading spaces
	    $Step =~ s/\s.*//; # Remove anything after a space
	}
    }
    die "$ERROR could not find simulation time in file $File!\n" if $Time < 0;
    die "$ERROR could not find time step in file $File!\n" if $Step < 0;

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

    # Form the name of restart tree if not specified already
    if(not $RestartTree){
	# Check if it is a time accurate run
	if($SimulationTime){
	    # If the time unit is not set try to guess it from simulation time
	    if(not $TimeUnit){
		my $Unit;
		foreach $Unit (keys %UnitSecond){
		    $TimeUnit = $Unit if $SimulationTime >= $UnitSecond{$Unit};
		}
	    }
	    $RestartTree = sprintf("RESTART_t%6.2f%s", 
				   $SimulationTime/$UnitSecond{$TimeUnit},
				   $TimeUnit);
	}else{
	    $RestartTree = sprintf "RESTART_n%6d", $nStep;
	}
	$RestartTree =~ s/ /0/g;
	print "# Restart.pl set restart tree name to $RestartTree/.\n"
	    if $Verbose;
    }

    # Check the restart tree directory
    die "$ERROR restart tree $RestartTree is in the way!\n" if -d $RestartTree;

    # Check output restart directories for alll components
    my $Comp;
    foreach $Comp (sort keys %RestartOutDir){
	next unless -d $Comp;

	my $Dirs = $RestartOutDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir};
	die "$ERROR could not find directory $Dirs!\n" unless -d $Dir;

	opendir(DIR,$Dir) or die "$ERROR could not open directory $Dir!\n";
	my @Content = readdir(DIR);
	closedir(DIR);
	die "$ERROR directory $Dir is empty!\n" unless $#Content > 1;

	# Check consistency of the simulation time
	my $File = "$Dir/$HeaderFile";
	&get_time_step($File) if -f $File;

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

    # Move the SWMF restart file
    my $File = "$RestartTree/$RestartOutFile";
    print "mv $RestartOutFile $File\n" if $Verbose;
    rename $RestartOutFile, $File or 
	die "$ERROR could not move $RestartOutFile into $File!";

    # Move the output restart directories of the components into the tree
    # and create empty output restart directories
    my $Comp;
    foreach $Comp (sort keys %RestartOutDir){
	next unless -d $Comp;
	my $Dirs = $RestartOutDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir};

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

    # Check for an existing restart file
    die "$ERROR file $RestartInFile is in the way!\n" if 
	(-f $RestartInFile and not -l $RestartInFile);

    # Check the SWMF restart file in the restart tree
    my $File = "$RestartTree/$RestartOutFile";
    die "$ERROR could not find restart file $File!\n" 
	unless (-f $File or $NoTreeCheck);

    # Set the step and the simulation time
    &get_time_step($File) unless $NoTreeCheck;

    my $Comp;
    foreach $Comp (sort keys %RestartInDir){
	next unless -d $Comp;
	my $Dirs = $RestartInDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir or -l $Dir};

	die "$ERROR could not find input restart directory/link $Dirs!\n" 
	    unless -d $Dir or -l $Dir;

	die "$ERROR could not find restart directory $RestartTree/$Comp!\n" 
	    unless (-d "$RestartTree/$Comp" or $NoTreeCheck);

	# Check the consistency of the simulation time
	my $File = "$RestartTree/$Comp/$HeaderFile";
	&get_time_step($File) if -f $File;

	print "# Restart.pl has checked $Dir\n" if $Verbose;
    }

    print "# Restart.pl has checked  input restart file and directories.\n";
}
##############################################################################
sub link_tree{

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

    my $Comp;
    foreach $Comp (sort keys %RestartInDir){
	next unless -d $Comp;
	my $Dirs = $RestartInDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir or -l $Dir};

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
sub print_help{
    print "
Purpose:

    Collect current output restart files into one directory tree and/or 
    link input restart files to a directory tree.

Usage:

    Restart.pl [-h] [-v] [-c] [-o] [-t=s|m|h|d|y] [[-i] DIR]

    -h         Print help message and exit.

    -v         Print verbose information.

    -c         Check only, do not actually create or link.
               Default is to create and link as specified by -i and -o.

    -o         Create restart tree from output directories only.
               Cannot be used together with the -i switch.
               Default is to link input directories to the tree as well.

    -i         Link input restart directories to the restart tree only.
               Cannot be used together with the -o switch.
               The name of the restart tree must be specified.
               Default is to create the restart tree first.

    -t=UNIT    Time unit to form the name of the restart tree.
               The UNIT can be given as one of the following characters:
               s, m, h, d, y corresponding to seconds, minute, hour, day and
               year respectively. The -t option has no effect if the 
               name of the restart tree is specified by the parameter DIR.
               The default time unit is the largest unit which does not 
               exceed the simulation time.

    DIR        Name of the restart directory tree. This argument
               must be specified if the -i switch is used. Otherwise
               the default name is RESTART_n012345 for steady state runs
               and RESTART_t012.34u for time accurate runs, where the
               numbers should be replaced with the actual time step and
               simulation time, and the 'u' with the actual time unit.

Examples:

    Check the output and input restart file and directories:

Restart.pl -c

    Create restart tree from current results and link input to it:

Restart.pl

    Create restart tree from current results and use hours as time unit:

Restart.pl -o -u=h

    Check linking to the existing RESTART_t002.00h tree

Restart.pl -i -c RESTART_t002.00h

    Link to the existing RESTART_t002.00h tree and print verbose info:

Restart.pl -v -i RESTART_t002.00h

";
    exit;
}
##############################################################################
