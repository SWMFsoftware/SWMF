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
my $IndexHtmlFile  = "index.html";
my $ParamHtmlFile  = "param.html";
my $EditorHtmlFile = "editor.html";
my $ManualHtmlFile = "manual.html";
my $ImageDir       = "share/Scripts";

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
open(MANUAL, ">$ManualHtmlFile") 
    or die "$ERROR: could not open output file $ManualHtmlFile\n";

# Global variables
my $Framework;                                   # True in framework mode
my $ValidComp = "SC,IH,SP,GM,IM,PW,RB,IE,UA,PS"; # List of components
my %CompVersion;                                 # Hash for component versions
my $nSession = 0;                                # Number of sessions

my @SessionRef;    # Data structure read from the XML param file
my %Editor;        # Hash for the editor parameters
my %Clipboard;     # Hash for the clipboard
my $CommandXml;    # XML description of the command
my $CommandExample;# XML description of the command
my $CommandText;   # Normal text description of the command

my %TableColor = ("command" => "\#CCCCCC",
		  "comment" => "\#DDDDDD",
		  "userinput"=>"\#CCCCCC");

&read_xml_file;

&write_index_html unless -f $IndexHtmlFile;

&write_editor_html;

&write_manual_html;

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
	    if($Framework and $3){
		$ValidComp = $3;
		%CompVersion = split( /[,\/]/ , $3 );
		print "ValidComp=$ValidComp\n" if $Debug;
		print "CompVersion=",join(', ',%CompVersion),"\n" if $Debug;
	    }
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

	    # Create new array element for this item
	    $SessionRef[$nSession]{SECTION}{$Section}{ITEM}[$iItem]{TYPE} = 
		$Type;

	    my $ItemRef=
		$SessionRef[$nSession]{SECTION}{$Section}{ITEM}[$iItem];
	    $ItemRef->{VIEW} = $View;
	    $ItemRef->{HEAD} = <XMLFILE>; # First line of body is always read

	    while($_=<XMLFILE>){
		last if /<\/ITEM>/;
		next if $ItemRef->{TAIL} and $View eq "MIN";
		$ItemRef->{TAIL} .= $_;
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
	    $ParamXml = "$Section/$CompVersion{$Section}/PARAM.XML";
	}
    }
    my $IsFirstSession = $Session eq $SessionRef[1]{NAME};

    open(FILE, $ParamXml) or die "$ERROR: could not open $ParamXml\n";

    my $CommandList;
    my $Command = $Editor{INSERT};

    if($Editor{ABC}){
	# Read commands and sort them alphabetically
	my @Command;
	while($_=<FILE>){
	    if(/^<command\s+name=\"([^\"]+)\"/){
		push(@Command,$1);
		&read_command_info if "\#$1" eq $Command;
	    }
	}
	$CommandList = 
	    "    <OPTION>\#".join("\n    <OPTION>\#",sort @Command)."\n";
    }else{
	# Read commands and command groups and form OPTIONs and OPTGROUPs
	while($_=<FILE>){
	    if(/^<command\s+name=\"([^\"]+)\"/){
		$CommandList .= "    <OPTION>\#$1\n";
		&read_command_info if "\#$1" eq $Command;
	    }
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

sub read_command_info{

    my $CommandInfo;
    while($_=<FILE>){
	last if /<\/command>/;
	$CommandInfo .= $_;
    }

    ($CommandXml, $CommandInfo)    = split(/\n\#/, $CommandInfo, 2);

    ($CommandExample, $CommandText) = split(/\n\n/, $CommandInfo, 2);

    $CommandExample = "\#$CommandExample\n";

    print "
CommandXml=
$CommandXml

CommandExample=
$CommandExample

CommandText=
$CommandText" if $Debug;

}

##############################################################################

sub write_index_html{

    open(FILE, ">$IndexHtmlFile") 
	or die "$ERROR: could not open $IndexHtmlFile\n";
    print FILE
"<FRAMESET ROWS=120,1*>
  <FRAME SRC=\"./$EditorHtmlFile\" NAME=EDITOR FRAMEBORDER=1>
  <FRAMESET COLS=75%,25%>
    <FRAME SRC=\"./$ParamHtmlFile#SELECTED\"  NAME=PARAMFILE FRAMEBORDER=1>
    <FRAME SRC=\"./$ManualHtmlFile\"  NAME=MANUAL FRAMEBORDER=1>
  </FRAMESET>
</FRAMESET>
";
    close(FILE);
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
      <TD ALIGN=LEFT WIDTH=20%>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=CHECK>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=SAVE>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=\"SAVE AS\">
      </TD>
      <TD ALIGN=CENTER WIDTH=60%>
         <FONT COLOR=red>
$ParamFile
         </FONT>
      </TD>
      <TD ALIGN=RIGHT WIDTH=20%>
<INPUT TYPE=SUBMIT NAME=SUBMIT VALUE=EXIT>
      </TD>
    </TR>
  </TABLE>
<hr>
  <CENTER>
  View: <SELECT NAME=SESSION>
$SessionSection
  </SELECT>
&nbsp;Insert: 
  <SELECT NAME=INSERT>
$InsertList
  </SELECT>
$InsertItem
  </CENTER>
  </FORM>
</body>
";
    print EDITOR $Editor;
    close EDITOR;
}

##############################################################################

sub write_manual_html{

    my $Manual;
    if($CommandExample){
	$Manual =  $CommandText;
	$Manual =~ s/\n\n/\n<p>\n/g;
	$Manual =  "<H1>Manual</H1>\n<PRE>\n$CommandExample\n</PRE>\n$Manual";
	$Manual =~ s/\n$//;
    }elsif($Editor{INSERT} =~ /^PASTE/){
	$Manual = "<H1>Clipboard</H1>\n<PRE>\n$Clipboard{BODY}\n</PRE>";
    }
    $Manual = $CommandText if 

    print MANUAL
"<BODY BGCOLOR=WHITE>
$Manual
</BODY>
";
    close MANUAL;
}

##############################################################################

sub write_param_html{

    my $Param = '  <head>
  <style type="text/css">
  a {text-decoration: none;}
  </style>
  </head>

  <BODY BGCOLOR=#DDDDD TEXT=BLACK LINK=BLUE VLINK=BLUE>
';

    my $InsertSessionButton;

    $InsertSessionButton = 
"  <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert session\"></A>
"                 if $Editor{SELECT} eq "ALL SESSIONS";

    ########################## SESSION #################################
    my $iSession;
    for $iSession (1..$nSession){
	my $SessionView =  $SessionRef[$iSession]{VIEW};
	my $SessionName =  $SessionRef[$iSession]{NAME};
	my $SessionGivenName;
	$SessionGivenName=": $SessionName" unless 
	    $SessionName =~ /^Session \d/;

	my $InsertSectionButton;

	$InsertSectionButton =
"    <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert section\"></A>
" 	    if $Editor{SELECT} !~ /(ALL SESSIONS|\/)/;

	# Place anchor to selected session
	$Param .= "<A NAME=SELECTED></A>\n" 
	    if $Editor{SELECT} eq $SessionName;

	$Param .=
"<hr COLOR=BLACK>
  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#CCCCCC>
    <TR>
      <TD WIDTH=20 ALIGN=LEFT>
";
	if($SessionView eq "MIN"){
	    $Param .=
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize session\"></A>
";
	}else{
	    $Param .=
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize session\"></A>
";
	}
	$Param .=
"      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
Session $iSession$SessionGivenName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSessionButton
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy session\"></A>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove session\"></A>
      </TD>
    </TR>
  </TABLE>
";
	next if $SessionView eq "MIN";

	######################## SECTION ###############################

	my $Section;
	for $Section (sort keys %{ $SessionRef[$iSession]{SECTION} }){
	    my $SectionRef  = $SessionRef[$iSession]{SECTION}{$Section};
	    my $SectionView = $SectionRef->{VIEW};

	    my $InsertItemButton;
	    $InsertItemButton = 
"     <A HREF=\"$IndexHtmlFile\"
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert\"></A>
"                   if $Editor{SELECT} =~ /\/$Section$/;


	    # Place anchor to selected session
	    $Param .= "<A NAME=SELECTED></A>\n" 
		if $Editor{SELECT} eq "$SessionName/$Section";

	    $Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#CCCCCC>
    <TR>
      <TD WIDTH=20>
      </TD>
      <TD COLSPAN=2>
<hr COLOR=GREY>
      </TD>
    </TR>
    <TR>
      <TD ALIGN=CENTER>
";
	    if($SectionView eq "MIN"){
		$Param .=
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize section\"></A>
";
	    }else{
		$Param .=
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize section\"></A>
";
	    }
	    $Param .=
"      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
Section: $Section
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy section\"></A>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove section\"></A>
      </TD>
    </TR>
  </TABLE>
";
	    next if $SectionView eq "MIN";
	    
###################### ITEM LOOP ############################################

	    my $iItem;
	    for $iItem (1..$#{ $SectionRef->{ITEM} }){

		my $ItemRef  = $SectionRef->{ITEM}[$iItem];
		my $ItemView = $ItemRef->{VIEW};
		my $ItemType = lc($ItemRef->{TYPE});
		my $ItemHead = $ItemRef->{HEAD};
		my $ItemTail = $ItemRef->{TAIL};

		my $TableColor = $TableColor{$ItemType};

		if($ItemType eq "userinput"){
		    # Fix the content of the user input
		    $ItemHead = "#USERINPUT";
		    $ItemTail =~ s/\n*\#USERINPUTEND.*//m;
		}
		my $nLine = ($ItemTail =~ s/\n/\n/g);

		if($ItemView eq "EDIT"){
		    $nLine += 2;
		    if($ItemType eq "comment"){
			$ItemTail = $ItemHead.$ItemTail;
			$ItemHead = "";
			$nLine++;
		    }
		    $Param .= "
  <FORM TARGET=_parent>
  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#BBEEFF>
    <TR>
      <TD WIDTH=20>
      </TD>
      <TD WIDTH=380><FONT COLOR=BLUE>
$ItemHead
      </FONT></TD>
      <TD ALIGN=RIGHT>
<IMG SRC=$ImageDir/button_save.gif TITLE=Save>
         &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<IMG SRC=$ImageDir/button_remove.gif TITLE=\"Cancel\">
      </TD>
    </TR>
    <TR>
      <TD>
      </TD>
       <TD COLSPAN=2>
          <TEXTAREA COLS=60 ROWS=$nLine>
$ItemTail
          </TEXTAREA>
       </TD>
    </TR>
  </TABLE>
  </FORM>
";
		    next;
		}

		my $CopyButton =
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy $ItemType\"></A>
";
		my $RemoveButton =
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove $ItemType\"></A>
";

		# Show first line with usual buttons
		$ItemHead = "<PRE>$ItemHead</PRE>" if $ItemType eq "comment";

		$Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=20 ALIGN=RIGHT>
";
		if($ItemTail){
		    if($ItemView eq "MIN"){
			$Param .= 
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
      ><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize $ItemType\"></A>
";
		    }else{
			$Param .= 
"      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
      ><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize $ItemType\"></A>
";
		    }
		}
		$Param .= 
"      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
$ItemHead
      </A></TD>
      <TD ALIGN=RIGHT ALIGN=TOP>
$InsertItemButton$CopyButton$RemoveButton
      </TD>
    </TR>
  </TABLE>
";

		if($ItemView eq "MIN"){
                    $Param .= "<p>\n" if $ItemType ne "comment";
		    next;
                }

		if($ItemType eq "comment" or $ItemType eq "userinput"){

		    $Param .= 
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=20>
      </TD>
      <TD WIDTH=380 COLSPAN=2><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
<PRE>$ItemTail</PRE>
      </A></TD>
    </TR>
";
		    if($ItemType eq "userinput"){
			$Param .= 
"    <TR>
      <TD WIDTH=20 ALIGN=RIGHT>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
      ><IMG SRC=$ImageDir/button_minimize_up.gif TITLE=\"Minimize $ItemType\"></A>
      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
#USERINPUT
      </A></TD>
      <TD ALIGN=RIGHT>
$CopyButton$RemoveButton
      </TD>
    </TR>
";
			$Param .= 
"  </TABLE>
";
		    }
		}else{  #command type item
                    $Param .= 
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#CCCCCC>
";
		    my $iLine;
		    for $iLine (0..$nLine-1){
			$ItemTail =~ s/(.*)\n//;
			my $Value;
			my $Comment;
			($Value,$Comment) = split(/\t|\s\s\s+/, $1, 2);

			$Param .= 
"    <TR>
      <TD WIDTH=20>
      </TD>
      <TD WIDTH=380>
$Value
      </TD>
      <TD>
$Comment
      </TD>
    </TR>
";

		    } # end body line loop

		    $Param .= 
"  </TABLE>
";
		} #endif item type

		$Param .= "<p>\n" if $ItemType ne "comment";

	    } # end item loop

	    ###### End section #########

	    $Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#CCCCCC>
    <TR>
      <TD WIDTH=20 ALIGN=CENTER>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_minimize_up.gif TITLE=\"Minimize section\"></A>
      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
Section: $Section
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertItemButton
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy section\"></A>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove section\"></A>
      </TD>
    </TR>
  </TABLE>
";
	} # Section loop

	###### End session #########

	$Param .= 
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=#CCCCCC>
    <TR>
      <TD>
      </TD>
      <TD COLSPAN=2>
<hr COLOR=GREY>
      </TD>
    </TR>
    <TR>
      <TD WIDTH=20 ALIGN=CENTER>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_minimize_up.gif TITLE=\"Minimize session\"></A>
      </TD>
      <TD WIDTH=380><A HREF=\"$IndexHtmlFile\" TARGET=_parent>
Session $iSession$SessionGivenName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy session\"></A>
      <A HREF=\"$IndexHtmlFile\" TARGET=_parent
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove session\"></A>
      </TD>
    </TR>
  </TABLE>
";

    } # session loop

    $Param .= 
"<hr COLOR=BLACK>\n
  <TABLE WIDTH=100% BGCOLOR=#CCCCCC><TR><TD ALIGN=RIGHT>
$InsertSessionButton
  </TD></TR></TABLE>
</BODY>\n";

    print PARAM $Param;
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
     param.html, editor.html and manual.html used by the GUI
     parameter editor. The PARAM.XML files are also read for 
     the list and description of commands. 

Usage:

  ParamXmlToHtml.pl [-h] [PARAMXMLFILE]

  -h            print help message and stop

  PARAMXMLFILE  Name of the input XML enhanced parameter file. 
                The default name is run/PARAM.in.xml

Examples:

    Convert run/PARAM.in.xml

ParamXmlToHtml.pl

    Convert run_new/PARAM.in.xml

ParamXmlToHtml.pl run_new/PARAM.in.xml"
#EOC
    ,"\n\n";
    exit 0;
}
