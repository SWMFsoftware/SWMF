#!/usr/bin/perl -s

my $Help       = $h or $H or $help;
my $OutputOnly = $o;
my $InputOnly  = $i;
my $CheckOnly  = $c;
my $Verbose    = $v;

use strict;

my $ERROR ="ERROR in Restart.pl:"; # Error message string
my $RestartTree    = "RESTART";    # Default directory for restart tree
my $RestartOutFile = "RESTART.out";# Default name for SWMF restart file
my $RestartInFile  = "RESTART.in"; # Default name for SWMF restart file

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

&print_help if $Help;

if($InputOnly and $OutputOnly){
    print "$ERROR do not use -i and -o together! ".
	"Type Restart.pl -h for help.";
    exit 1;
}

if($ARGV[0]){
    $RestartTree = $ARGV[0];
    $RestartTree =~ s/\/+$//; # Remove trailing /
}

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
sub create_tree_check{

    # Check the SWMF restart file
    die "$ERROR could not find restart file $RestartOutFile!\n" 
	unless -f $RestartOutFile;

    # Check the restart directory
    die "$ERROR restart tree $RestartTree is in the way!\n" if -d $RestartTree;

    # Check output restart directories for alll components
    my $Comp;
    foreach $Comp (sort keys %RestartOutDir){
	next unless -d $Comp;

	my $Dirs = $RestartOutDir{$Comp};
	my $Dir;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir};
	print "# Restart.pl is checking $Dirs\n" if $Verbose;
	die "$ERROR could not find directory $Dirs!\n" unless -d $Dir;

	opendir(DIR,$Dir) or die "$ERROR could not open directory $Dir!\n";
	my @Content = readdir(DIR);
	closedir(DIR);
	die "$ERROR directory $Dir is empty!\n" unless $#Content > 1;
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

    my $Comp;
    foreach $Comp (sort keys %RestartInDir){
	next unless -d $Comp;
	my $Dirs = $RestartInDir{$Comp};
	my $Dir;
	print "# Restart.pl is checking $Dirs\n" if $Verbose;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir or -l $Dir};

	die "$ERROR could not find input restart directory/link $Dirs!\n" 
	    unless -d $Dir or -l $Dir;

	die "$ERROR could not find restart directory $RestartTree/$Comp!\n" 
	    unless (-d "$RestartTree/$Comp" or $NoTreeCheck);
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
    print "# Restart.pl has linked to restart tree $RestartTree/.\n";
}
##############################################################################
sub print_help{
    print "
Purpose:

    Collect current output restart files into one directory tree and/or 
    link input restart files to a directory tree.

Usage:

    Restart.pl [-h] [-v] [-c] [-o|-i] [DIR]

    -h         Print help message and exit.
    -v         Print verbose information.
    -c         Check only, do not actually create or link.
               Default is to create and link as specified by -i and -o.
    -o         Create restart tree from output directories only.
               Default is to link input directories to the tree as well.
    -i         Link input restart directories to the restart tree only.
               Default is to create the restart tree first.
    DIR        Name of the restart directory tree. 
               Default is RESTART.

Examples:

    Create RESTART tree from current results and link input to it:

Restart.pl

    Check the output and input restart file and directories:

Restart.pl -c

    Create RESTART_2hr tree from current output restart file and directories
    and link input restart file and directories to it. Print verbose info:

Restart.pl -v RESTART_2hr

    Create the default RESTART tree but do not link to it:

Restart.pl -o

    Check linking to the existing RESTART_4hr tree

Restart.pl -i -c RESTART_4hr

    Link to the existing RESTART_4hr tree

Restart.pl -i RESTART_4hr

";
    exit;
}
##############################################################################
