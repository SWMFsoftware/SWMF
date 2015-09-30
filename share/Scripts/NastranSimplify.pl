#!/usr/bin/perl
use strict;

if($#ARGV < 1 or $ARGV[0] =~ /-h/i){
    print "
Convert raw NASTRAN file into a better formatted NASTRAN file
that starts with comments that provide the number of triangles, 
edges and nodes for the code reading the file. The code also
checks if the shape satisfies the expected topological constraint
that the number_of_faces - number_of_edges + number_of_vertices=2.

Usage: NastranSimplify.pl INPUT.bdf OUTPUT.bdf
";
    exit;
}

my $infile  = $ARGV[0];
my $outfile = $ARGV[1];

open(FILE, $infile) or die "Could not open input file $infile\n";

my @lines = <FILE>;

close(FILE);

print "Number of lines in $infile: ",$#lines+1,"\n";

my $npoint;
my $ntriangle;

my %edge;

foreach (@lines){
    $npoint++ if /GRID/;
    $ntriangle++ if /CTRIA/;

    if(/CTRIA3\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)/){
	$edge{"$1;$2"}++ if $1 < $2;
	$edge{"$2;$1"}++ if $1 > $2;
	$edge{"$1;$3"}++ if $1 < $3;
	$edge{"$3;$1"}++ if $1 > $3;
	$edge{"$2;$3"}++ if $2 < $3;
	$edge{"$3;$2"}++ if $2 > $3;
	die "Bad indexes in line $_" if $1 == $2 or $1 == $3 or $2 == $3;
    }
}

print "There are $npoint points and $ntriangle triangles in $infile\n";

my $nedge =  scalar(keys %edge);
print "There are $nedge edges\n";
print "Face-edge+vertex=", $ntriangle-$nedge+$npoint," should be 2\n";

foreach (keys %edge){
    print "Edge $_ is repeated $edge{$_} times\n" if $edge{$_} != 2;
}

open(FILE, ">$outfile")  or die "Could not open output file $outfile\n";;

print FILE "\$ Created by NastranSimplify.pl\n";
print FILE "\$ $npoint    POINTS\n";
print FILE "\$ $ntriangle TRIANGLES\n";

foreach (@lines){
    s/\s*\*N\d+\s*\n//;
    s/\*N\d+\s*/ /;
    s/ +/ /g;
    print FILE $_;
}

close(FILE);
