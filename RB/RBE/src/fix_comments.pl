#!/usr/bin/perl
undef $/;

$file = $ARGV[0] or die "Usage: fix_comments.pl FILE.f"
open(FILE, $file);
$_ = <FILE>;
s/\n     \*/\&\n      /gm;
s/^c/\!/g;
close(FILE);

$file =~ s/f$/f90/;

open(FILE, "$file");
print FILE $_;
close(FILE);
