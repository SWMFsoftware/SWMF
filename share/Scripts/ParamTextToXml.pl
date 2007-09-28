#!/usr/bin/perl -s

# Read a parameter file and add some XML info useful for the parameter editor.
# Type ParamTextToXml.pl -h for help.

# Read command line options
my $Help        = $h; undef $h;

use strict;

# Set framework mode to true if CON/Control directory is present
my $Framework = -d "CON/Control";

# Default list of possible components (for SWMF)
my $ValidComp = 'SC|IH|SP|GM|IM|PS|PW|RB|IE|UA';

# Name of file containing component versions (SWMF)
my $MakefileDef = "Makefile.def";

# Error string
my $ERROR   = 'ParamTextToXml_ERROR:';
my $WARNING = 'ParamTextToXml_WARNING:';

# Set default values (needed in the help message too)
my $ParamFile  = ($ARGV[0] or 'run/PARAM.in');

my $ExpandParamScript= 'share/Scripts/ExpandParam.pl';

# Print help message and exit if -h switch was used
&print_help if $Help;

die "$ERROR could not find $ParamFile\n" unless -f $ParamFile;

my $ParamExpand = "${ParamFile}.expand";

# Flatten the parameter file
if(not -x $ExpandParamScript){
    warn "$WARNING $ExpandParamScript is not found: ".
	"cp $ParamFile $ParamExpand\n";
    my $error=`cp $ParamFile $ParamExpand`;
    die "$ERROR in cp $ParamFile $ParamExpand: $error\n"
	if $error;
}else{
    my $error=`$ExpandParamScript $ParamFile > $ParamExpand`;
    die "$ERROR in $ExpandParamScript $ParamFile > $ParamExpand: $error\n"
	if $error;
}

# Set list of component versions
my $Components;
if($Framework){
    open(INFILE,$MakefileDef) or
	print "$ERROR could not open file $MakefileDef\n";

    while(<INFILE>){
	next unless /^(\w\w)_VERSION\s*=\s*(\w+)/;
	$Components .= "$1/$2," unless $2 eq "Empty";
    }
    $Components =~ s/,$//;
    close(INFILE);
}

open(INFILE,$ParamExpand) or
    die "$ERROR Could not open expanded parameter file $ParamExpand!\n";

my $XmlFile = "${ParamFile}.xml";

open(XMLFILE,">$XmlFile") or
    die "$ERROR Could not open output XML file $XmlFile!\n";

# Initialize for checking input parameter file
my $nLine        = 0;   # Line number in the expanded file
my $nSession=1;         # current session number
my @NameSession;        # names of the sessions
my $Section="";         # name of component section
my $IsCommand=0;        # distingish commands from comments and user input
my $UserInput="";       # true inside #BEGINUSERINPUT and #ENDUSERINPUT
my $iItem=0;            # command index in the current section

my @Output;             # Array of hashes of arrays (session/section/item)
my @Clipboard;          # The lines after the #END command

my @Line;               # Collect lines for outputfile

while($_=<INFILE>){

    if(/^\#END\b/){
	# Read the rest of the file into the 'clip board'
	@Clipboard = <INFILE>;
	last;
    }
    
    $nLine++;

    if($UserInput and /^\#(BEGIN_COMP|END_COMP|RUN|USERINPUTBEGIN)\b/){
	&print_error( " for command $_".
		      "\tthis command cannot occur after ".
		      "#USERINPUTBEGIN at line $UserInput");
	$UserInput = 0;
    }

    # Check for BEGIN_COMP and END_COMP commands
    if(/^\#BEGIN_COMP\b/){
	if(not $Framework){
	    &print_error(" for command $_".
			 "\tshould not be used in stand-alone mode");
	    next;
	}
	if($Section){
	    &print_error(" for command $_".
			 "\talready had BEGIN_COMP $Section");
	    next;
	}
	# Figure out which component is beginning here
	($Section) = /BEGIN_COMP ([A-Z][A-Z])/ or
	    &print_error(" for command $_".
			 "\tincorrectly formatted BEGIN_COMP command");
	    
	# Set the item index to the last one of this session
	$iItem = $#{ $Output[$nSession]{$Section} } + 1;

    }elsif(/^\#END_COMP/){
	if(not $Framework){
	    &print_error(" for command $_".
			 "\tshould not be used in stand-alone mode");
	    next;
	}

	if(not $Section){
	    &print_error(" for command $_".
                         "\tmissing \#BEGIN_COMP command");
            next;
        }

	# Extract name of the component from #END_COMP ID
	my $Comp;
	($Comp) = /END_COMP ([A-Z][A-Z])/ or
	    &print_error(" for command $_".
			 "\tincorrectly formatted END_COMP command");

	# Check if the component name matches
	if($Comp ne $Section){
	    &print_error(" for command $_\tcomponent does not match".
			 " BEGIN_COMP $Section");
	}

	# In any case return to CON
	$Section = '';

	# Set item index for CON
	$iItem = $#{ $Output[$nSession]{$Section} } + 1;
 
    }elsif(/^\#RUN/){		# Check if a new session has started
	# Check if the required commands are defined and
	# if the parameters are correct for the session

	if($Section and $Framework){
	    print "Error in session $nSession ending at line $nLine ".
                "in file $ParamExpand:\n".
                "\tBEGIN_COMP $Section does not have matching END_COMP\n";
	    $Section = '';
	}
	$nSession++;
	$iItem=0;

    }elsif( /^(Begin|End)\s+session:\s*(.*)/i ){
	$NameSession[$nSession] = $2;
    }else{
	# Analyze line for various items
	if(/^\#USERINPUTBEGIN/){
	    $iItem++;
	    $UserInput = $nLine+1;
	    $IsCommand=1;
	}elsif(/\#USERINPUTEND/){
	    if(not $UserInput){
		&print_error(" for command $_\tthere is no matching".
			     " USERINPUTBEGIN command");
	    }else{
		$UserInput = 0;
	    }
	}elsif(/^\#/ and not $UserInput){
	    $iItem++;
	    $IsCommand=1;
	}elsif( /^\s*\n$/ and $IsCommand and not $UserInput){
	    $iItem++;
	    $IsCommand=0;
	}
	# Store line
	$Output[$nSession]{$Section}[$iItem].=$_;
    }
}

# Write the XMLFILE
if($Framework){ 
    print XMLFILE "\t\t\t<MODE FRAMEWORK=\"1\" COMPONENTS=\"$Components\"/>\n";
}else{
    print XMLFILE "\t\t\t<MODE FRAMEWORK=\"0\"/>\n";
}

my $iSession;
for $iSession (1..$nSession){
    print XMLFILE
	"\t\t\t<SESSION NAME=\"$NameSession[$iSession]\" VIEW=\"MAX\">\n";
    my $SessionRef = $Output[$iSession];

    foreach $Section (sort keys %$SessionRef){
	print XMLFILE "\t\t\t\t<SECTION NAME=\"$Section\" VIEW=\"MAX\">\n";
	my $SectionRef = $Output[$iSession]{$Section};

	for $iItem (0..$#{ $SectionRef }){
	    my $Item = $Output[$iSession]{$Section}[$iItem];

	    $Item =~ s/^\s*//m;   # Remove starting empty lines
	    next unless $Item;    # Ignore empty items

	    my $Type;
	    if($Item =~ /^\#USERINPUTBEGIN/){
		$Type="USERINPUT";
	    }elsif($Item =~ /^\#/){
		$Type="COMMAND";
	    }else{
		$Type="COMMENT";
	    }
	    print XMLFILE "\t\t\t\t\t<ITEM TYPE=\"$Type\" VIEW=\"MAX\">\n";
	    print XMLFILE $Item;
	    print XMLFILE "\t\t\t\t\t</ITEM>\n";
	}
	print XMLFILE "\t\t\t\t</SECTION>\n";
    }
    print XMLFILE "\t\t\t</SESSION>\n";
}

my $Clipboard;
$Clipboard = join("",@Clipboard);
$Clipboard =~ s/^\n*//m;

print XMLFILE "
\t\t\t<EDITOR SELECT=\"ALL SESSIONS\" INSERT=\"NEW SESSION\" ".
    "FILE=\"\" ABC=\"0\"/>
\t\t\t<CLIPBOARD SESSION=\"NONE\" SECTION=\"NONE\" TYPE=\"NONE\">
$Clipboard\t\t\t</CLIPBOARD>\n";

close XMLFILE;

exit 0;

##############################################################################
sub print_error{
    print "Error at line $nLine in file $ParamExpand",@_,"\n";
}

##############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Convert Parameter File to XML for Editor}
#!ROUTINE: ParamTextToXml.pl - convert a parameter file to XML for editor
#!DESCRIPTION:
# This script is used internally by the parameter editor.
# 
#!REVISION HISTORY:
# 9/27/2007 G.Toth - initial version
#
#EOP
sub print_help{

    print 
#BOC
"Purpose:

     Read a parameter file (and all included files) and output a single
     parameter file with extra XML tags that is used by the parameter editor.
     The output is written into a file named as the parameter file with 
     an extra .xml extension. If the CON/Control directory is present,
     the code assumes that it is run for the SWMF, otherwise it is running
     in stand-alone mode. In framework mode the code reads the component
     versions from Makefile.def. In stand alone mode the #BEGIN_COMP 
     and #END_COMP commands are ignored.
     This script uses the ExpandParam.pl script.

Usage:

  ParamTextToXml.pl [-h] [PARAMFILE]

  -h            print help message and stop

  PARAMFILE     The plain text file containing the input parameters 
                and include commands for other files (if any).
                The default file name is run/PARAM.in.

Examples:

    Create run/PARAM.in.xml from run/PARAM.in:

ParamTextToXml.pl

    Create run_new/PARAM.in.xml from run_new/PARAM.in:

ParamTextToXml.pl run_new/PARAM.in"
    ,"\n\n";
    exit 0;
}
