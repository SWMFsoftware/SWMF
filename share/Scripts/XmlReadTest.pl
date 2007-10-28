#!/usr/bin/perl -s

my $compact = $c;

use strict;

require Data::Dumper or
    die "ERROR in XmlReadTest: Data::Dumper package is not installed\n";
require "XmlRead.pl" or
    die "ERROR in XmlReadTest: XmlRead.pl is not found\n";

my $Result;
$Result = &XmlRead( join('',<>) );

$Data::Dumper::Indent=0 if $compact;

print Data::Dumper->Dump([$Result],["Result"]);

exit 0;
