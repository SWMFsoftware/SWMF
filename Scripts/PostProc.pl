#!/usr/bin/perl -s

my $Help    = ($h or $H or $help);
my $Verbose = ($v or $verbose);
my $Gzip    = ($g or $gzip);

use strict;

&print_help if $Help;

my $Pwd = `pwd`; chop $Pwd;

my $Dir;
foreach $Dir ("GM", "IH", "SC"){
    next unless -d "$Dir/IO2";
    print "cd $Dir\n" if $Verbose;
    chdir $Dir or die "Could not change directory to $Dir\n";
    &shell("./pIDL");
    if($Gzip){
	&shell("./pTEC g");
    }else{
	&shell("./pTEC p r");
    }
    chdir $Pwd;
}    

$Dir = "IE";
if(-d "$Dir/ionosphere"){
    print "cd $Dir\n" if $Verbose;
    chdir $Dir or die "Could not change directory to $Dir\n";
    if($Gzip){
	&shell("./pION -g");
    }else{
	&shell("./pION");
    }
    chdir $Pwd;
}


exit 0;

#############################################################

sub shell{
    my $command = join(" ",@_);
    print "$command\n" if $Verbose;
    my $result = `$command`;
    print $result if $Verbose;
}

##############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Post-process plot files with Scripts/PostProc.pl}
#!ROUTINE: PostProc.pl - post-process plot files of the components
#!DESCRIPTION:
# This script is copied into the run directory and it should be executed there.
# The script post processes the plot files created by the components.
#
#!REVISION HISTORY:
# 02/12/2005 G. Toth - initial version
#EOP

sub print_help{
    print
#BOC
'Purpose:

   Post-process the plot files.

Usage:

   PostProc.pl [-h] [-v] [-g]

   -h -help    Print help message and exit.

   -v -verbose Print verbose information.

   -g -gzip    Gzip the ASCII files.

Examples:

   Post-process the plot files:

PostProc.pl

   Post-process the plot files, compress them and print verbose info:

PostProc.pl -g -v'
#EOC
    ,"\n\n";
    exit;
}
##############################################################################
