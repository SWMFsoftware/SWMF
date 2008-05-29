#!/usr/bin/perl -s

my $Help     = ($h or $help); undef $h; 
my $Commands = $c;            undef $c;
my $Force    = $f;            undef $f;

use strict;

&print_help if $Help or not @ARGV;

my $ERROR = "ERROR in XmlToF90.pl:";

die "$ERROR cannot use -f and -c switches together!\n" if $Force and $Commands;
die "$ERROR input file argument is missing!\n" unless @ARGV;

require 'share/Scripts/XmlRead.pl';

my $InputFile  = $ARGV[0];
my $OutputFile = $ARGV[1];

open(IN, $InputFile) or die "$ERROR could not open input file $InputFile\n";
my $Input = join('',<IN>);
close(IN);

my $Update = ($OutputFile and -f $OutputFile and not $Force);

my $NewFile = ($OutputFile and not ($Update or $Commands));

my $Indent     = (' ' x 5);  # initial indentation for case statements
my $Indent1    = (' ' x 3);  # incremental indentation
my $IndentCont = (' ' x 5);  # indentation for continuation lines

my $F90;                # F90 source code written by process_xml
my $iPart;              # part index for multi-part parameters
my $ForeachName;        # name of the index variable in a foreach loop
my $ForeachValue;       # value of the index variable in a foreach loop
my %VariableType;       # Hash for variable types
my %DefaultValue;       # Hash for default values

# put in UseStrict
$VariableType{"UseStrict"} = "logical";
$DefaultValue{"UseStrict"} = "F";

my $Xml = &XmlRead($Input); # XML tree data structure created from XML text

&process_xml($Xml);

if($NewFile){
    my $Subroutine = $OutputFile;
    $Subroutine = $1 if $Subroutine =~ /\/([^\/]+)$/;  # remove path
    $Subroutine =~ s/\..*$//;                          # remove extension

    # replace variables with their default values if possible
    foreach my $Variable (keys %DefaultValue){
	$_ = $DefaultValue{$Variable};
	while(/\$(\w+)/ and $DefaultValue{$1}){
	    s/\$(\w+)/$DefaultValue{$1}/;
	    $DefaultValue{$Variable} = $_;
	}
    }

    my $Declarations;
    foreach my $Variable (sort keys %VariableType){
	my $Type  = $VariableType{$Variable};
	my $Value = $DefaultValue{$Variable};

	$Type =~ s/strings?/character(len=100)/;

	# Fix value
	$Value =~ s/T/.true./ if $Type eq "logical";
	$Value =~ s/F/.false./ if $Type eq "logical";
	$Value =  "'$Value'" if $Type =~ /character/ and length($Value);
	$Value .= ".0" if $Type eq "real" and $Value =~ /^\d+$/;

	# Add declaration
	if(length($Value)){
	    $Declarations .= "  $Type :: $Variable = $Value\n";
	}else{
	    $Declarations .= "  $Type :: $Variable\n";
	}
    }
    $F90 = 
"subroutine $Subroutine

  use ModReadParam, ONLY: i_session_read, read_line, read_command, read_var
  use ModUtilities, ONLY: split_string

  implicit none

  integer :: iSession
  character (len=100) :: NameCommand
  character(len=*), parameter:: NameSub = '$Subroutine'

$Declarations

  !---------------------------------------------------------------------------
  iSession = i_session_read()
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
$F90
     case default
        !if(iProc==0) then
        write(*,*) NameSub // ' WARNING: unknown command ' // &
             trim(NameCommand),' !'
        if(UseStrict)call CON_stop('Correct PARAM.in!')
        !end if
     end select
  end do

contains

  logical function is_first_session()
    is_first_session = iSession == 1
    if(iSession > 1)then
       ! if(iProc==0) then
       write(*,*) NameSub // ' WARNING: command ',trim(NameCommand), &
            ' can be used in first session only!'
       if(UseStrict)call CON_stop('Correct PARAM.in!')
       ! end if
    end if
  end function is_first_session

end subroutine $Subroutine
";
}
if($OutputFile){
    open(OUT, ">$OutputFile") or 
	die "$ERROR could not open output file $OutputFile\n";
    print OUT $F90;
    close(OUT);
}else{
    print $F90;
}

exit 0;

##############################################################################
sub process_xml{

    # recursive subroutine that processes the XML file

    my $content = shift;

    foreach my $element (@$content){

        next if $element->{type} eq 't'; # Ignore elements of type text

	my $name = lc( $element->{"name"} );

	if($name eq 'commandlist'){
	    &process_xml($element->{content});
	}elsif($name eq 'commandgroup'){
	    $F90 .= "\n!!! $element->{attrib}->{name}\n";
	    &process_xml($element->{content});
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
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	}elsif($name eq 'parameter'){
	    my $Attrib = $element->{attrib};
	    my $Name   = perl_to_f90($Attrib->{name});
	    my $Type   = $Attrib->{type};
	    my $Case   = $Attrib->{case};
	    my $If     = perl_to_f90($Attrib->{if});

	    # Store variable name and type
	    $VariableType{$Name} = $Type;
	    $DefaultValue{$Name} = $Attrib->{default};;

	    # Create line
	    $F90 .= $Indent."if($If) &\n".$IndentCont if $If;
	    $F90 .= $Indent."call read_var('$Name', $Name)\n";

	    $F90 =~ s/\)\n$/, IsUpperCase=.true.)\n/ if $Case eq "upper";
	    $F90 =~ s/\)\n$/, IsLowerCase=.true.)\n/ if $Case eq "lower";

	    if($Type eq "strings"){
		my $MaxPart = $Attrib -> {max};
		$VariableType{"nStringPart"} = "integer";
		$VariableType{"StringPart_I(100)"} = "string";
		$F90 .= $Indent . "call split_string" . 
		    "($Name, $MaxPart, StringPart_I, nStringPart)\n";
		$iPart = 0;
		&process_xml($element->{content});
	    }
	}elsif($name eq 'part'){
	    my $Name = $element->{attrib}->{name};
	    $VariableType{$Name} = "string";
	    $iPart++;
	    $F90 .= $Indent."$Name = StringPart_I($iPart)\n";
	}elsif($name eq 'if'){
	    my $Expr = perl_to_f90($element->{attrib}->{expr});
	    $F90 .= $Indent."if($Expr)then\n";
	    $Indent .= $Indent1;
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	    $F90 .= $Indent."end if\n";
	}elsif($name eq 'for'){
	    my $Attrib = $element->{attrib};
	    my $From  =  perl_to_f90($Attrib->{from});
	    my $To    =  perl_to_f90($Attrib->{to});
	    my $Index = (perl_to_f90($Attrib->{name}) or "i");
	    $VariableType{$Index} = "integer";
	    $F90 .= $Indent."do $Index = $From, $To\n";
	    $Indent .= $Indent1;
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	    $F90 .= $Indent."end do\n";
	}elsif($name eq 'foreach'){
	    my $Attrib = $element->{attrib};
	    $ForeachName = $Attrib->{name};
	    foreach (split(/,/, $Attrib->{values})){
		$ForeachValue = $_;
		&process_xml($element->{content});
	    }
	    $ForeachName = ''; $ForeachValue = '';
	}
    }
}

###############################################################################

sub perl_to_f90{

    $_ = shift;

    # replace special variables provided by the CheckParam.pl script
    s/\$_command/NameCommand/ig;
    s/\$_namecomp/NameComp/ig;

    # replace foreach variable with actual value
    s/\$$ForeachName/$ForeachValue/g if $ForeachName;

    # remove all dollar signs from variable names
    s/\$//g;

    # convert relation operator
    s/ eq / == /g;
    s/ ne / \/= /g;
    s/ and / .and. /g;
    s/ or / .or. /g;
    s/ != / \/= /g;
    s/\bnot /.not. /g;

    # replace string matching (this is not quite right!)
    s/(\w+)\s*=~\s*\/([^\/]+)\/i?/index($1, "$2") > 0/g;

    s/(\w+)\s*\!~\s*\/([^\/]+)\/i?/index($1, "$2") < 1/g;

    # Remove \b from patterns (this is not quite right!)
    s/\\b//g;

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

   XmlToF90 [-h] [-f | -c=COMMANDS] infile [outfile]

-h          This help message

-c=COMMANDS COMMANDS is a comma separated list of commands to be transformed
            from XML to F90. Default is transforming all the commands.

-f          Force creating a new output file even if it already exists. 
            This flag cannot be combined with the -c=COMMAND switch.
            Default is to update the (selected) commands only in an existing
            file.

infile      Input XML file.

outfile     Output F90 file. Default is writing to STDOUT.

Examples:

            Generate F90 source code for a few commands and write to STDOUT:

share/Scripts/XmlToF90.pl -c=MAGNETICAXIS,ROTATIONAXIS Param/PARAM.XML

            Update all commands in an existing F90 source file

share/Scripts/XmlToF90.pl PARAM.XML src/set_parameters.f90"
#EOC
    ,"\n\n";
    exit 0;
}
