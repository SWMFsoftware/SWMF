#!/usr/bin/perl -s

$Verbose = $v;

die "Usage: RemoveFiles.pl [-v] FILENAMESTART\n" unless $ARGV[0];
opendir(DIR,'.');
@files = grep {/^$ARGV[0]/ and not -d} readdir(DIR);
closedir(DIR);
die "No files match $ARGV[0] !\n" unless @files;
print "Matching files: @files\n" if $Verbose;
print "Remove all ",$#files+1," files starting with $ARGV[0] [y/n] ?\n";
$reply = <STDIN>;
die "Cancelled removal !\n" if $reply !~ /^y|yes$/;
unlink @files;
