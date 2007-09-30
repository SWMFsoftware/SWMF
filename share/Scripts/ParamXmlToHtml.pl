#!/usr/bin/perl -s

# Read an XML enhanced parameter file and write out an HTML page for the editor
# Type ParamXmlToHtml.pl -h for help.

# This subroutine is outside the scope of strict and the "my" variables
# This is a safety feature, so the global variables can not be changed
sub eval_comp{eval("package COMP; $_[0]")}

# Read command line options
my $Debug     = $D; undef $D;
my $Help      = $h; undef $h;

use strict;

# Read optional argument
my $XmlFile   = ($ARGV[0] or 'run/PARAM.in.xml');

# Print help message and exit if -h switch was used
&print_help if $Help;

# Error string
my $ERROR = 'ParamXmlToHtml_ERROR';

# Output files have fixed names
my $ParamHtmlFile  = "param.html";
my $EditorHtmlFile = "editor.html";

# Name of the original parameter file
my $ParamFile = $XmlFile;
$ParamFile =~ s/\.xml//;

# Read Perl dump of XML description: no need for the XML reader and faster
my $TreeFile  = 'PARAM.pl';

open(XMLFILE, $XmlFile) 
    or die "$ERROR: could not open input file $XmlFile\n";
open(PARAM, ">$ParamHtmlFile") 
    or die "$ERROR: could not open output file $ParamHtmlFile\n";
open(EDITOR, ">$EditorHtmlFile") 
    or die "$ERROR: could not open output file $EditorHtmlFile\n";

# Global variables
my $Framework;                                   # True in framework mode
my $ValidComp = "SC,IH,SP,GM,IM,PW,RB,IE,UA,PS"; # List of components
my $nSession = 0;                                # Number of sessions

my @SessionRef;    # Data structure read from the XML param file
my %Editor;        # Hash for the editor parameters
my %Clipboard;     # Hash for the clipboard

&read_xml_file;

&write_editor_html;

&write_param_html;

exit 0;

##############################################################################

sub read_xml_file{

    my $iItem;
    my $SessionView;
    my $Section;
    my $SectionView;

    while($_ = <XMLFILE>){
	if(/<MODE FRAMEWORK=\"(\d)\" (COMPONENTS=\"([^\"]*)\")?/){
	    $Framework = $1;
	    $ValidComp = $3 if $Framework and $3;
	    print "Framework=$Framework ValidComp=$ValidComp\n" if $Debug;
	}elsif(/<SESSION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $nSession++;
	    my $Name=($1 or "Session $nSession");
	    $SessionView=$2;
	    $SessionRef[$nSession]{NAME} = $Name;
	    $SessionRef[$nSession]{VIEW} = $SessionView;
	    print "nSession=$nSession Name=$Name View=$SessionView\n" 
		if $Debug;
	}elsif(/<SECTION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    my $Name=$1;
	    $SectionView = $2;
	    $Section = ($Name or "CON");
	    $SessionRef[$nSession]{SECTION}{$Section}{VIEW} = $SectionView;
	    $iItem=0;

	    print "Section=$Section View=$SectionView\n" if $Debug;

	    # Do not read details if the section is minimized
	    if($SectionView eq "MIN" or $SessionView eq "MIN"){
		$_=<XMLFILE> while not /<\/SECTION>/;
	    }		
	}elsif(/<ITEM TYPE=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $iItem++;
	    my $Type=$1;
	    my $View=$2;

	    print "iItem=$iItem Type=$Type View=$View\n" if $Debug;

	    my $ItemRef=
		$SessionRef[$nSession]{SECTION}{$Section}{ITEM}[$iItem];
	    $ItemRef->{TYPE} = $Type;
	    $ItemRef->{VIEW} = $View;
	    while($_=<XMLFILE>){
		last if /<\/ITEM>/;
		$ItemRef->{BODY} .= $_ unless $View eq "MIN";
	    }
	}elsif(/<EDITOR 
	       \ SELECT=\"([^\"]*)\"
	       \ INSERT=\"([^\"]*)\"
	       \ FILE  =\"([^\"]*)\"
	       \ ABC   =\"([^\"]*)\"/x){
	    $Editor{SELECT} = $1;
	    $Editor{INSERT} = $2;
	    $Editor{FILE}   = $3;
	    $Editor{ABC}    = $4;
	}elsif(/<CLIPBOARD
	       \ SESSION=\"([^\"]+)\"
	       \ SECTION=\"([^\"]+)\"
	       \ TYPE   =\"([^\"]+)\"/x){
	    $Clipboard{SESSION} = $1;
	    $Clipboard{SECTION} = $2;
	    $Clipboard{TYPE}    = $3;
	    while($_=<XMLFILE>){
		last if /<\/CLIPBOARD>/;
		$Clipboard{BODY} .= $_;
	    }
	}
    }
    close XMLFILE;
}

##############################################################################
sub command_list{

    my $ParamXml = "PARAM.XML";
    my $Session = $Editor{SELECT};

    if($Framework){
	my $Section;
	($Session, $Section) = split(/\//, $Session, 2);
	if($Section eq "CON"){
	    $ParamXml = "Param/PARAM.XML";
	}else{
	    $ParamXml = "$Section/PARAM.XML";
	}
    }
    my $IsFirstSession = $Session eq $SessionRef[1]{NAME};

    open(FILE, $ParamXml) or die "$ERROR: could not open $ParamXml\n";

    my $CommandList;
    if($Editor{ABC}){
	# Read commands and sort them alphabetically
	my @Command;
	while($_=<FILE>){
	    push(@Command,$1) if /^<command\s+name=\"([^\"]+)\"/;
	}
	$CommandList = 
	    "    <OPTION>\#".join("\n    <OPTION>\#",sort @Command)."\n";
    }else{
	# Read commands and command groups and form OPTIONs and OPTGROUPs
	while($_=<FILE>){
	    $CommandList .= "    <OPTION>\#$1\n" 
		if /^<command\s+name=\"([^\"]+)\"/;
	    $CommandList .= 
		"    </OPTGROUP>\n".
		"    <OPTGROUP LABEL=\"$1\">>\n" 
		if /^<commandgroup\s+name=\"([^\"]+)\"/;
	}
	# Remove first </OPTGROUP> and add a last one
	$CommandList =~ s/^    <\/OPTGROUP>\n//;
	$CommandList .=   "    <\/OPTGROUP>\n";
    }
    close FILE;

    return $CommandList;

    #my $tree=do($TreeFile);
    #die "Could not read command definitions from $TreeFile\n" unless $tree;
    #print "Description tree was read from file $TreeFile\n" if $Debug;
    #my $CommandList = $tree->[0];
    #if($CommandList->{type} ne 'e' or $CommandList->{name} ne 'commandList'){
    #    die "$ERROR first node should be an element named commandList\n";
    #}
    
}
##############################################################################

sub write_editor_html{

    # Read the parameter definitions
    # $CommandList = &read_tree($TreeFile) 

    print "Starting write_editor_html\n" if $Debug;

    my $SessionSection="    <OPTION>ALL SESSIONS\n";

    my $iSession;
    for $iSession (1..$nSession){
	print "iSession=$iSession\n" if $Debug;

	my $SessionRef = $SessionRef[$iSession];
	my $SessionName = $SessionRef->{NAME};
        $SessionSection .= "    <OPTION>$SessionName\n";
	next unless $Framework;
	my $Section;
	for $Section (sort keys %{$SessionRef->{SECTION}}){
	    $SessionSection .= 
		"      <OPTION>".('&nbsp;' x 3)."$SessionName/$Section\n";
	}
    }

    # Add SELECTED
    my $Selected = $Editor{SELECT};
    $SessionSection =~ s/<OPTION>(.*$Selected)/<OPTION SELECTED>$1/;

    print "SessionSection=$SessionSection\n" if $Debug;

    my $InsertList;

    if($Selected eq "ALL SESSIONS"){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SESSION\n" if $Clipboard{TYPE} eq 
	    "SESSION";
        $InsertList .= "    <OPTION>NEW SESSION\n";

    }elsif($Selected !~ /\// and $Framework){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SECTION\n" if $Clipboard{TYPE} eq 
	    "SECTION";
        $InsertList .= "    <OPTION>Section CON\n";
	my $Comp;
	for $Comp (split ',', $ValidComp){
	    $InsertList .= "     <OPTION>Section $Comp\n";
	}
    }else{
        $InsertList  =     "    <OPTION>FILE\n";
	$InsertList .=     "    <OPTION>PASTE COMMAND/COMMENT\n" if 
	    $Clipboard{TYPE} eq "ITEM";
	$InsertList .=     "    <OPTION>NEW COMMENT\n";
    	if($Editor{ABC}){
	    $InsertList .= "    <OPTION>COMMANDS ALPHABETICALLY\n";
	}else{
	    $InsertList .= "    <OPTION>COMMANDS BY GROUP\n";
	}
	$InsertList     .= &command_list;
    }

    # Add SELECTED
    my $Insert = $Editor{INSERT};
    if($Insert){
	$InsertList =~ s/<OPTION>(.*$Insert)/<OPTION SELECTED>$1/;
    }

    my $InsertItem;

    if($Insert eq "FILE"){
	my @Files;
	@Files = glob("run/PARAM*.in* Param/PARAM*.in.*");

	my $Files;
	$Files = "    <OPTION>".join("\n    <OPTION>",@Files) if @Files;
	$InsertItem=
"  <SELECT NAME=FILE>
    <OPTION>SELECT FILE
$Files
  </SELECT>
";
	# Add SELECTED
	my $File = $Editor{FILE};
	$InsertItem =~ s/<OPTION>(.*$File\n)/<OPTION SELECTED>$1/ if $File;
    }elsif($Insert =~ /^PASTE/){
	$InsertItem = "
  <TEXTAREA ROWS=2 COLS=40 READONLY>
$Clipboard{BODY}
  </TEXTAREA>
";
    }
    # Add checkbox if COMMAND list
    if($InsertList =~ /COMMAND/){
	$InsertItem .= 
	    " <INPUT TYPE=CHECKBOX NAME=ALPHABETIC VALUE=1>Abc\n";
	$InsertItem =~ s/>Abc$/ CHECKED>Abc/ if $Editor{ABC};
    }

    chop $SessionSection;
    chop $InsertList;
    chop $InsertItem;
    my $Editor = "
<body bgcolor=WHITE>

<FORM ALIGN=RIGHT>
  <TABLE WIDTH=100%>
    <TR>
      <TD ALIGN=LEFT>
         <FONT COLOR=red>
$ParamFile
         </FONT>
      </TD>
      <TD ALIGN=CENTER>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=CHECK>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=SAVE>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=SAVE AS>
      </TD>
      <TD ALIGN=RIGHT>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=EXIT>
      </TD>
    </TR>
  </TABLE>
  </FORM>
<hr>
  <FORM>
  <TABLE BORDER=0>
    <TR>
      <TD ALIGN=CENTER>Session/section</TD>
      <TD ALIGN=CENTER COLSPAN=2>Insert item
      </TD>
    </TR>
    <TR>
      <TD>
  <SELECT NAME=SESSION>
$SessionSection
  </SELECT>
      </TD>
      <TD>
  <SELECT NAME=INSERT>
$InsertList
  </SELECT>
      </TD>
      <TD>
$InsertItem
      </TD>
    </TR>
  </TABLE>
  </FORM>
  <hr>
</body>
";
    print EDITOR $Editor;
    close EDITOR;
}

##############################################################################

sub write_param_html{
    close PARAM;
}

##############################################################################
#BOP
#!ROUTINE: ParamXmlToHtml.pl - convert XML enhanved parameter to HTML
#!DESCRIPTION:
# This script is called internally from the parameter editor.
# 
#!REVISION HISTORY:
# 09/29/2007 G.Toth - initial version
#EOP
sub print_help{

    print 
#BOC
"Purpose:

     Read an XML enhanced parameter file and write HTML files 
     param.html and editor.html used by the parameter editor. 
     The PARAM.pl files are also read for the list and 
     description of commands. 

Usage:

  ParamXmlToHtml.pl [-h] [PARAMXMLFILE]

  -h            print help message and stop

  PARAMXMLFILE  Name of the input XML enhanced parameter file. 
                The default name is run/PARAM.in.xml

Examples:

    Convert run/PARAM.in.xml into editor.html and param.html

ParamXmlToHtml.pl

    Convert run_new/PARAM.in.xml into editor.html and param.html

ParamXmlToHtml.pl run_new/PARAM.in.xml"
#EOC
    ,"\n\n";
    exit 0;
}
