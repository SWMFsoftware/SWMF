#!/usr/bin/perl -s
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Automated Documentation}
#!ROUTINE: XmlToTex.pl - generate Latex documentation from XML definitions of input parameters
#!DESCRIPTION:
# Create Latex documentation based on the XML definitions of the 
# input parameters typically found in the PARAM.XML files.
# This script allows to store the parameter descriptions in a single XML file
# which is suitable for automated parameter checking and GUI generation,
# yet provide the same information in a well formatted and printable manual
# as well. The specific format of the PARAM.XML files is described by the
# share/Scripts/CheckParam.pl script and the manual.
#
#!REVISION HISTORY:
# 03/22/2004 G.Toth - initial version
# 08/13/2004          preserve _{.. for subscripts in Latex.
#EOP
if($h|$help|$H|$Help){
    print 
#BOC
"Purpose:

   Convert XML description of input commands into LaTex description.

Usage:

   XmlToTex [-h] [infile] [> outfile]

-h       This help message

infile   Input file. Default is reading from STDIN.

outfile  Output file. Default is writing to STDOUT.

Example: 

  share/Scripts/XmlToTex.pl Param/PARAM.XML > Param/PARAM.xmltex"
#EOC
    ,"\n\n";
    exit 0;
}
#BOC
use strict;
my ($verbatim, $comment, $rule, $start);

while(<>){

    # Skip lines with lots of exclamation marks used for human reading
    next if /!!!!!!!!/;

    # commandList --> section
    if(/^\s*<commandList\s+name=[\'\"]([^\'\"]+)/){
	$_="\\clearpage\n\\section\{Input Commands for the $1\}\n\n";
        $start = 1;
    }

    # commandgroup --> subsection
    if(/^\s*<commandgroup\s+name=[\'\"]([^\'\"]+)/){
	$_="\\subsection\{\u\L$1\}\n\n";
    }

    # skip things before the commandList
    next unless $start; 

    # command --> subsubsection
    if(/^\s*<command/){
	/name=[\'\"]([^\'\"]+)/;
	my $command   = "\#$1";
	my @path = split('/',$ARGV);
	my $comp = 'CON'; 
	$comp="$path[$#path-2]/$path[$#path-1]"  # UA/GITM2
	    if $path[$#path-1] ne 'Param';       # Param/PARAM.XML is for CON

	$comp =~ s/GM\/BATSRUS/GM,SC,IH\/BATSRUS/; # GM --> GM,SC,IH

	my $index = "$comp\!$command";# Form the index term
	$_="\\subsubsection\{$command command\}\\index\{$index\}\n\n";
    }

    # #COMMAND or #COMMAND ID --> verbatim
    if(/^\#[A-Z0-9_]+( [A-Z][A-Z])?\s*$/){
	$_='\begin'.'{verbatim}'."\n$_";
	$verbatim = 1;
    }

    # verbatim part ends with an empty line
    if($verbatim and /^\s*$/){
	$_='\end'.'{verbatim}'."\n";
	$verbatim = 0;
    }

    # One line XML comments are replaced with Tex comments
    s/<!--(.*?)-->/\%$1/;

    # Beginning of multiline XML comment
    $comment = 1 if s/<!--.*//;

    # End of XML comment
    $comment = 0 if s/.*-->//;

    # Beginning of rule
    $rule=1 if /\s*<rule/;

    # End of rule
    $rule=0 if s/.*<\/rule>//;

    # Skip XML lines, XML comments and rules
    next if /[<>]/ or $comment or $rule;

    # Replace special XML characters with Tex
    s/\&lt;/\$<\$/g;
    s/\&gt;/\$>\$/g;
    s/\&amp;/\\&/g;

    # Replace TAB characters with spaces
    1 while s/\t+/' ' x (length($&) * 8 - length($`) % 8)/e;

    if(not $verbatim){
	# Put \ before special Tex character #
	s/\#/\\\#/g;
	# Put \ before special Tex character _ but not before _{
	# This allows the use of _{1} subscript in math mode
	s/(\_[^\{])/\\$1/g;
	# Do not put two \ before # and _
	s/\\\\([\#_])/\\$1/g;
    }

    # Remove ! from the beginning of the lines
    s/^! ?//;

    # Print out the line
    print;
}
#EOC
