#!/usr/bin/perl -s

my $Help       = $h or $H or $help;
my $OutputOnly = $o;
my $InputOnly  = $i;
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

$RestartTree = $ARGV[0] if $ARGV[0];

&create_tree if not $InputOnly;
&link_tree   if not $OutputOnly;

exit 0;
##############################################################################
sub create_tree{
    die "$ERROR restart tree $RestartTree is in the way!\n" if -d $RestartTree;
    mkdir $RestartTree,0777 
	or die "$ERROR restart tree $RestartTree could not be created!\n";

    # Take care of the SWMF restart file
    die "$ERROR could not find restart file $RestartOutFile\n" 
	unless -f $RestartOutFile;
    my $File = "$RestartTree/$RestartOutFile";
    rename $RestartOutFile, $File or 
	die "$ERROR could not move $RestartOutFile into $File";

    my $Comp;
    foreach $Comp (sort keys %RestartOutDir){
	next unless -d $Comp;
	my $Dirs = $RestartOutDir{$Comp};
	my $Dir;
	print "Comp=$Comp checking Dirs=$Dirs\n" if $Verbose;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir};
	die "$ERROR could not find directory $Dir!\n" unless -d $Dir;
	opendir(DIR,$Dir) or die "$ERROR could not open directory $Dir\n";
	my @Content = readdir(DIR);
	closedir(DIR);
	die "$ERROR directory $Dir is empty\n" unless  $#Content > 0;
	rename $Dir, "$RestartTree/$Comp" or 
	    die "$ERROR could not move $Dir into $RestartTree/$Comp\n";
	mkdir $Dir, 0777 or die "$ERROR could not create directory $Dir\n";
    }
}
##############################################################################
sub link_tree{

    die "$ERROR restart tree $RestartTree is missing!\n" 
	unless -d $RestartTree;

    unlink $RestartInFile if -l $RestartInFile;
    die "$ERROR file $RestartInFile is in the way\n" if -f $RestartInFile;

    # Link in the SWMF restart file in the restart tree
    my $File = "$RestartTree/$RestartOutFile";
    die "$ERROR could not find restart file $File\n" unless -f $File;
    symlink $File, $RestartInFile or 
	die "$ERROR could not link $File to $RestartInFile";

    my $Comp;
    foreach $Comp (sort keys %RestartInDir){
	next unless -d $Comp;
	my $Dirs = $RestartInDir{$Comp};
	my $Dir;
	print "Comp=$Comp checking Dirs=$Dirs\n" if $Verbose;
	foreach (split /,/,$Dirs){$Dir=$_; last if -d $Dir or -l $Dir};
	die "$ERROR could not find directory/link $Dir!\n" unless 
	    -d $Dir or -l $Dir;
	unlink $Dir if -l $Dir;
	if(-d $Dir){
	    rmdir $Dir or die "$ERROR could not remove directory $Dir\n";
	}
	symlink "../$RestartTree/$Comp", $Dir or 
	    die "$ERROR could not link $RestartTree/$Comp to $Dir\n";
    }
}
##############################################################################
sub print_help{
    print "
Purpose:

    Collect current output restart files into one directory tree and/or 
    link input restart files to a directory tree.

Usage:

    Restart.pl [-h|-i|-o] [-v] [DIR]

    -h         Print help message and exit.
    -v         Print verbose information.
    -o         Create restart tree from output directories only.
               Default is to link input directories to the tree as well.
    -i         Link input restart directories to the restart tree only.
               Default is to create the restart tree first.
    DIR        Name of the restart directory tree. 
               Default is RESTART.

Examples:

    Create RESTART tree from current results and link input to it:

Restart.pl

    Create RESTART_2hr tree from current results and link input to it:

Restart.pl RESTART_2hr

    Only create the default RESTART tree

Restart.pl -o

    Link to the existing RESTART_4hr tree

Restart.pl -i RESTART_4hr

";
    exit;
}
##############################################################################
