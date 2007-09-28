#!/usr/bin/perl -s

# Read an XML enhanced parameter file and recover the pure text format
# Type ParamXmlToText.pl -h for help.

# Read command line options
my $Help        = $h; undef $h;

use strict;

# Error string
my $ERROR   = 'ParamXmlToTex_ERROR:';
my $WARNING = 'ParamXmltToTex_WARNING:';

# Set default values
my $XmlFile  = ($ARGV[0] or 'run/PARAM.in.xml');

# Print help message and exit if -h switch was used
&print_help if $Help;

open(INFILE,$XmlFile) or
    die "$ERROR Could not open XML parameter file $XmlFile!\n";

my $TextFile = $XmlFile.".txt";
$TextFile =~ s/\.xml//i;

open(TEXTFILE,">$TextFile") or
    die "$ERROR Could not open output text file $TextFile!\n";

# Initialize some variables that need to be remembered
my $iSession;
my $NameSession;
my $NameSection;
my $IsCommand;
my $Output;

while($_=<INFILE>){
    
    if(/<SESSION .*NAME=\"([^\"]*)/){
	$iSession++;
	$NameSession=($1 or $iSession);
	$Output .= "Begin session: $NameSession\n\n";
    }elsif(/<\/SESSION>/){
	$Output .= "\nEnd session: $NameSession\n#END ".("\#" x 60)."\n";
    }elsif(/<SECTION .*NAME=\"(\w*)/){
	$NameSection = $1;
	$Output .= "\#BEGIN_COMP $NameSection ".("-" x 40)."\n\n"
	    if $NameSection;
    }elsif(/<\/SECTION>/){
	$Output .= "\#END_COMP $NameSection   ".("-" x 40)."\n\n"
	    if $NameSection;
    }elsif(/<ITEM TYPE=\"(COMMAND|USERINPUT)/){
	$IsCommand=1;
    }elsif(/<\/ITEM>/){
	$Output .= "\n" if $IsCommand;
	$IsCommand=0;
    }elsif(not /^\s*</){
	$Output .= $_;
    }
}

# Replace #END with #RUN when a session follows
$Output =~ s/\#END (\#+\nBegin session:)/\#RUN $1/;

# Replace double empty lines with a single one
1 while ($Output =~ s/\n\n\n/\n\n/g);

print TEXTFILE $Output;
close TEXTFILE;

exit 0;

##############################################################################
#BOP
#!QUOTE: \subsection{Convert XML extended parameter File to plain text}
#!ROUTINE: ParamXmlToText.pl - convert a XML parameter file to text
#!DESCRIPTION:
# This script is used internally by the parameter editor.
# 
#!REVISION HISTORY:
# 9/28/2007 G.Toth - initial version
#
#EOP
sub print_help{

    print 
#BOC
"Purpose:

     Read an XML enhanced parameter file and output a plain text file.
     The output filename ends with a .txt extension.

Usage:

  ParamTextToXml.pl [-h] [XMLFILE]

  -h            print help message and stop

  XMLFILE       The file containing the XML enhanced input parameters. 
                The default file name is run/PARAM.in.xml
Examples:

    Create run/PARAM.in.txt from run/PARAM.in.xml:

ParamXmlToText.pl

    Create run_new/PARAM.in.txt from run_new/PARAM.in.xml:

ParamXmlToText.pl run_new/PARAM.in.xml"
    ,"\n\n";
    exit 0;
}
