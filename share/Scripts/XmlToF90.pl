#!/usr/bin/perl -s


my $Help     = ($h or $help); undef $h; 
my $Commands = $c;            undef $c;

use strict;

&print_help if $Help or not @ARGV;

require 'share/Scripts/XmlRead.pl';

my $ERROR = "ERROR in XmlToF90.pl:";

my $InputFile  = $ARGV[0];
my $OutputFile = $ARGV[1];

open(IN, $InputFile) or die "$ERROR could not open input file $InputFile\n";

my $Input = join('',<IN>);
close(IN);

my $Indent  = "    ";
my $Indent1 = "    ";
my $F90;
$F90 = $Indent."select case(NameCommand)\n"; 
my $Xml = &XmlRead($Input);

&process_elements($Xml);

$F90 .= $Indent."end select\n";

print $F90;

exit 0;

##############################################################################
sub process_elements{

    # recursive subroutine that processes the XML file

    my $content = shift;

    foreach my $element (@$content){

        next if $element->{type} eq 't'; # Ignore elements of type text

	my $name = lc( $element->{"name"} );

	if($name eq 'commandlist'){
	    &process_elements($element->{content});
	}elsif($name eq 'commandgroup'){
	    $F90 .= $Indent."! $element->{attrib}->{name}\n";
	    &process_elements($element->{content});
	}elsif($name eq 'command'){
	    my $Attrib = $element->{attrib};
            my $Name   = "\"\#$Attrib->{name}\"";
	    my $Alias  = $Attrib->{alias};
	    foreach my $Alias (split ",", $Attrib->{alias}){
		$Name .= ", \"\#$Alias\"";
	    }
	    $F90 .= $Indent."case($Name)\n";
	    $Indent .= $Indent1;
	    my $If = $Attrib->{if};
	    if($If =~ /IsFirstSession/){
		$F90 .= $Indent."if(.not.is_first_session())CYCLE\n";
	    }
	    &process_elements($element->{content});
	    $Indent =~ s/$Indent1//;
	}elsif($name eq 'parameter'){
	    my $Attrib = $element->{attrib};
	    my $Name   = $Attrib->{name};
	    my $If     = perl_to_f90($Attrib->{if});
	    $F90 .= $Indent."if($If) &\n".$Indent1 if $If;
	    $F90 .= $Indent."call read_var('$Name', $Name)\n";
	}elsif($name eq 'if'){
	    my $Expr = perl_to_f90($element->{attrib}->{expr});
	    $F90 .= $Indent."if($Expr)then\n";
	    $Indent .= $Indent1;
	    &process_elements($element->{content});
	    $Indent =~ s/$Indent1//;
	    $F90 .= $Indent."end if\n";
	}elsif($name eq 'for'){
	    my $Attrib = $element->{attrib};
	    my $From = perl_to_f90($Attrib->{from});
	    my $To   = perl_to_f90($Attrib->{to});
	    my $Index = (perl_to_f90($Attrib->{name}) or "i");
	    $F90 .= $Indent."do $Index = $From, $To\n";
	    $Indent .= $Indent1;
	    &process_elements($element->{content});
	    $Indent =~ s/$Indent1//;
	    $F90 .= $Indent."end do\n";
	}
    }
}
###############################################################################
sub perl_to_f90{

    $_ = shift;

    # replace special variables
    s/\$_command/NameCommand/g;

    # remove dollar signs
    s/\$//g;

    # convert logical relations
    s/ eq / == /g;
    s/ and / .and. /g;
    s/ or / .or. /g;
    s/ != / \/= /g;
    s/\bnot /.not. /g;

    # replace string matching (this is not exact!)
    s/(\w+)\s*=~\s*\/([^\/]+)\/i?/index($1, "$2") > 0/g;

    return $_;
}

###############################################################################
#BOP
#!ROUTINE: XmlToF90.pl - generate F90 source from XML definitions of input parameters
#!DESCRIPTION:
# Generate F90 source code based on the XML definitions of the
# input parameters typically found in the PARAM.XML files.
# This script allows to store the parameter descriptions in a single XML file
# which is suitable for automated parameter checking, manual and GUI 
# generation, as well as generating the F90 code that reads in the parameters.
# The specific format of the PARAM.XML files is described by the
# share/Scripts/CheckParam.pl script and the manual.
#
#!REVISION HISTORY:
# 05/27/2008 G.Toth - initial version
#EOP
sub print_help{

    print 
#BOC
"Purpose:

   Convert XML description of input commands into F90 source code.

Usage:

   XmlToF90 [-h] [-c=COMMANDS] infile [outfile]

-h          This help message

-c=COMMANDS COMMANDS is a comma separated list of commands to be transformed
            from XML to F90. Default is transforming all the commands.

infile      Input XML file.

outfile     Output F90 file. Default is writing to STDOUT.
            If output file already exists, only update it.
            This allows hand made modifications to be preserved.

Examples:

            Generate F90 source code for a few commands and write to STDOUT:

share/Scripts/XmlToF90.pl -c=MAGNETICAXIS,ROTATIONAXIS Param/PARAM.XML

            Update all commands in an existing F90 source file

share/Scripts/XmlToF90.pl PARAM.XML src/set_parameters.f90"
#EOC
    ,"\n\n";
    exit 0;
}
