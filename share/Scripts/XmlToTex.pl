#!/usr/bin/perl -s

if($h|$help|$H|$Help){
    print "
Purpose:

   Convert XML description of input commands into LaTex description.

Usage:

   XmlToTex [-h] [infile] [> outfile]

-h       This help message

infile   Input file. Default is reading from STDIN.

outfile  Output file. Default is writing to STDOUT

Example: 

  XmlToTex ../../Param/PARAM.XML > SWMF_commands.tex
";
    exit 0;
}

use strict;
my ($verbatim, $comment, $rule, $start);

while(<>){

    # Skip lines with lots of exclamation marks used for human reading
    next if /!!!!!!!!/;

    # commandList --> section
    if(/^\s*<commandList\s+name=[\'\"]([^\'\"]+)/){
	$_="\\section\{Input Commands for the $1\}\n\n";
        $start = 1;
    }

    # commandgroup --> subsection
    if(/^\s*<commandgroup\s+name=[\'\"]([^\'\"]+)/){
	$_="\\subsection\{\u\L$1\}\n\n";
    }

    next unless $start; # skip things before the first command group

    # command --> subsubsection
    if(/^\s*<command/){
	/name=[\'\"]([^\'\"]+)/;
	$_="\\subsubsection\{\#$1 command\}\n\n";
    }

    # #COMMAND or #COMMAND ID --> verbatim
    if(/^\#[A-Z0-9_]+( [A-Z][A-Z])?\s*$/){
	$_='\begin{verbatim}'."\n$_";
	$verbatim = 1;
    }

    # verbatim part ends with an empty line
    if($verbatim and /^\s*$/){
	$_='\end{verbatim}'."\n";
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

    # Put \ before special Tex characters, but not 2 backslashes
    if(not $verbatim){
	s/([\#_])/\\$1/g unless $verbatim;
	s/\\\\([\#_])/\\$1/g;
    }

    # Remove ! from the beginning of the lines
    s/^! ?//;

    # Print out the line
    print;
}
