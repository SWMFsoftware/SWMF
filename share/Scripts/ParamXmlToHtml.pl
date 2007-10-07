#!/usr/bin/perl -s

# Read an XML enhanced parameter file and write out an HTML page for the editor
# Type ParamXmlToHtml.pl -h for help.

# This subroutine is outside the scope of strict and the "my" variables
# This is a safety feature, so the global variables can not be changed
sub eval_comp{eval("package COMP; $_[0]")}

# Read command line options
my $Debug     = $D; undef $D;
my $Help      = $h; undef $h;
my $Submit    = $submit; undef $submit;

use strict;

# Read optional argument
my $XmlFile   = ($ARGV[0] or 'run/PARAM.in.xml');

# Print help message and exit if -h switch was used
&print_help if $Help;

# Error string
my $ERROR = 'ParamXmlToHtml_ERROR';

# Output files have fixed names
my $IndexHtmlFile  = "index.php";
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
my $Framework = 0;                               # True in framework mode
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

&modify_xml_data if $Submit;

&write_xml_file if $Submit;

&write_index_html unless -f $IndexHtmlFile;

&write_editor_html;

&write_manual_html;

&write_param_html;

exit 0;

##############################################################################

sub read_xml_file{

    my $iItem;
    my $iSection;
    my $iSession;
    
    my $SessionView;

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
	    $iSection=0;
	    my $Name=($1 or "Session $nSession");
	    $SessionView=$2;
	    $SessionRef[$nSession]{NAME} = $Name;
	    $SessionRef[$nSession]{VIEW} = $SessionView;
	    print "nSession=$nSession Name=$Name View=$SessionView\n" 
		if $Debug;
	}elsif(/<SECTION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $iSection++;
	    $iItem=0;
	    my $Name=$1;
	    my $SectionView = $2;
	    $SessionRef[$nSession]{SECTION}[$iSection]{VIEW} = $SectionView;
	    $SessionRef[$nSession]{SECTION}[$iSection]{NAME} = $Name;

	    print "Section=$iSection name=$Name View=$SectionView\n" if $Debug;

	    # Do not read details if the section is minimized
	    if(not $Submit 
	       and ($SectionView eq "MIN" or $SessionView eq "MIN")){
		$_=<XMLFILE> while not /<\/SECTION>/;
	    }		
	}elsif(/<ITEM TYPE=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $iItem++;
	    my $Type=$1;
	    my $View=$2;

	    # Create new array element for this item
	    $SessionRef[$nSession]{SECTION}[$iSection]{ITEM}[$iItem]{TYPE} = 
		$Type;

	    my $ItemRef=
		$SessionRef[$nSession]{SECTION}[$iSection]{ITEM}[$iItem];
	    $ItemRef->{VIEW} = $View;
	    $ItemRef->{HEAD} = <XMLFILE>; # First line of body is always read

	    while($_=<XMLFILE>){
		last if /<\/ITEM>/;
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

sub modify_xml_data{


    warn "Submit = $Submit\n";

    $_ = $Submit;

    if( /^CHECK\b/ ){
	my $Error = `./TestParam.pl`;
	warn $Error if $Error; # to be written
    }elsif( /^SAVE\b/ ){
	warn "share/Scripts/ParamXmlToText.pl $XmlFile\n";
	`share/Scripts/ParamXmlToText.pl $XmlFile`;
    }elsif( /^SAVE AS\b/ ){
	# to be written
    }elsif( /^EXIT\b/ ){
	# kill the job
	kill(-9, getpgrp);
    }elsif( /^ABC_ON/ ){
	$Editor{ABC}=1;
    }elsif( /^ABC_OFF/ ){
	$Editor{ABC}=0;
    }else{
	my %Form;
	%Form = split (/[:;]/, $Submit);
	warn join(', ',%Form),"\n";

	$_ = $Form{action};
	my $id= $Form{id};

	if( /select_session/ ){
	    $Editor{INSERT} = "new" unless 
		($Editor{SELECT} =~ s/(\d+)/$1/g) == ($id =~ s/(\d+)/$1/g);
	    $Editor{SELECT} = $id;
	    return;
	}elsif( /select_insert/ ){
	    $Editor{INSERT} = $id;
	    return;
	}elsif( /select_file/ ){
	    $Editor{FILE} = $id;
	    if(open(MYFILE, $id)){
		$Clipboard{BODY} = join('',<MYFILE>);
		if($Editor{SELECT} eq "all"){
		    $Clipboard{TYPE}="SESSION";
		}elsif($Editor{SELECT} =~ /^\d+$/){
		    $Clipboard{TYPE}="SECTION";
		}else{
		    $Clipboard{TYPE}="COMMAND";
		}
	    }
	    return;
	}

	my $iSession = $1;
	my $iSection = $2;
	my $iItem    = $3;

	($iSession,$iSection,$iItem) = split(/,/,$id);

	warn "iSession=$iSession iSection=$iSection iItem=$iItem\n";
	
	my $SessionRef = $SessionRef[$iSession];
	my $NameSession= $SessionRef->{NAME};
	my $SectionRef = $SessionRef->{SECTION}[$iSection];
	my $NameSection= ($SectionRef->{NAME} or "CON");
	my $ItemRef    = $SectionRef->{ITEM}[$iItem];

	# View the stuff represented by id if id consist of numbers
	$Editor{SELECT}=$id if $id =~ /^[\d,]+$/;

	if( /minimize_session/ ){
	    $SessionRef->{VIEW}="MIN";
	}elsif( /minimize_section/ ){
	    $SectionRef->{VIEW}="MIN";
	}elsif( /minimize_item/ ){
	    $ItemRef->{VIEW}="MIN";
	}elsif( /maximize_session/ ){
	    $SessionRef->{VIEW}="MAX";
	}elsif( /maximize_section/ ){
	    $SectionRef->{VIEW}="MAX";
	}elsif( /maximize_item/ ){
	    $ItemRef->{VIEW}="MAX";
	}elsif( /edit_session/ ){
	    $SessionRef->{VIEW}="EDIT";
	}elsif( /edit_item/ ){
	    $ItemRef->{VIEW}="EDIT";
	}elsif( /remove_session/ or /copy_session/ ){
	    $Clipboard{SESSION} = $iSession;
	    $Clipboard{SECTION} = "NONE";
	    $Clipboard{TYPE}    = "SESSION";
	    $Clipboard{BODY}    = "Session $id $NameSession: should be here\n";
	    $Editor{SELECT}     = "all";
	    $Editor{INSERT}     = "PASTE SESSION";
	    if(/remove_session/){
		splice(@SessionRef,$iSession,1);
		$nSession--;
	    }
	}elsif( /remove_section/ or /copy_section/ ){
	    $Clipboard{SESSION} = $iSession;
	    $Clipboard{SECTION} = $NameSection;
	    $Clipboard{TYPE}    = "SECTION";
	    $Clipboard{BODY}    = "Section $id $NameSection: should be here\n";
	    $Editor{SELECT}     = $iSession;
	    $Editor{INSERT}     = "PASTE SECTION";
	    if( /remove_section/ ){
		splice (@{$SessionRef->{SECTION}}, $iSection, 1);
	    }
	}elsif( /remove_item/ or /copy_item/ ){
	    $Clipboard{SESSION} = $iSession;
	    $Clipboard{SECTION} = $NameSection;
	    $Clipboard{TYPE}    = $ItemRef->{TYPE};
	    $Clipboard{BODY}    = $ItemRef->{HEAD}.$ItemRef->{TAIL};
	    $Editor{SELECT}     = "$iSession,$iSection";
	    $Editor{INSERT}     = "PASTE COMMAND/COMMENT";
	    if( /remove_item/ ){
		splice (@{$SectionRef->{ITEM}}, $iItem, 1);
	    }
	}elsif( /insert_session/ ){
	    my $NewSessionRef;
	    $NewSessionRef->{VIEW} = "MAX";
	    $NewSessionRef->{NAME} = "Session $iSession";
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{HEAD} = $Clipboard{BODY};
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{TYPE} = "COMMENT";
	    $Editor{SELECT} = "all";
	    splice (@SessionRef, $iSession, 0, $NewSessionRef);
	    $nSession++;
	}elsif( /insert_section/ ){
	    my $NewSectionRef;
	    $NewSectionRef->{VIEW} = "MAX";
	    $NewSectionRef->{NAME} = $Clipboard{SECTION};
	    $NewSectionRef->{ITEM}[1]{HEAD} = $Clipboard{BODY};
	    $NewSectionRef->{ITEM}[1]{TYPE} = "COMMENT";

	    splice (@{$SessionRef->{SECTION}}, $iSection, 0, $NewSectionRef);
	}elsif( /insert_item/ ){
	    my $NewItemRef;
	    $NewItemRef->{VIEW} = "MAX";
	    $NewItemRef->{TYPE} = $Clipboard{TYPE};
            ($NewItemRef->{HEAD}, $NewItemRef->{TAIL}) 
		= split(/\n/, $Clipboard{BODY}, 2);
	    $NewItemRef->{HEAD}.="\n";

	    splice (@{$SectionRef->{ITEM}}, $iItem, 0, $NewItemRef);
	}elsif( /cancel_item_edit/ ){
	    $ItemRef->{VIEW}="MAX";
	}
    }
}

##############################################################################

sub write_xml_file{

    open(XMLFILE, ">$XmlFile") 
	or die "$ERROR: could not open $XmlFile for output!\n";

    # Write the XMLFILE
    print XMLFILE "\t\t\t<MODE FRAMEWORK=\"$Framework\"";
    print XMLFILE " COMPONENTS=\"$ValidComp\"" if $Framework;
    print XMLFILE "/>\n";

    my $iSession;
    for $iSession (1..$#SessionRef){
	print XMLFILE
	    "\t\t\t<SESSION NAME=\"$SessionRef[$iSession]{NAME}\" ".
	    "VIEW=\"$SessionRef[$iSession]{VIEW}\">\n";

	my $iSection;
	for $iSection (1..$#{ $SessionRef[$iSession]{SECTION} }){
	    my $SectionRef = $SessionRef[$iSession]{SECTION}[$iSection];

	    # Skip empty session ($iItem starts with 1)
	    next unless $#{ $SectionRef->{ITEM} } > 0;

	    print XMLFILE "\t\t\t\t<SECTION NAME=\"$SectionRef->{NAME}\"".
		" VIEW=\"$SectionRef->{VIEW}\">\n";
	    
	    my $iItem;
	    for $iItem (1..$#{ $SectionRef->{ITEM} }){

		my $Item = $SectionRef->{ITEM}[$iItem];

		print XMLFILE "\t\t\t\t\t<ITEM TYPE=\"$Item->{TYPE}\"".
		    " VIEW=\"$Item->{VIEW}\">\n";
		print XMLFILE $Item->{HEAD}, $Item->{TAIL};
		print XMLFILE "\t\t\t\t\t</ITEM>\n";
	    }
	    print XMLFILE "\t\t\t\t</SECTION>\n";
	}
	print XMLFILE "\t\t\t</SESSION>\n";
    }

    print XMLFILE "
<EDITOR SELECT=\"$Editor{SELECT}\" INSERT=\"$Editor{INSERT}\" FILE=\"$Editor{FILE}\" ABC=\"$Editor{ABC}\"/>
<CLIPBOARD SESSION=\"$Clipboard{SESSION}\" SECTION=\"$Clipboard{SECTION}\" TYPE=\"$Clipboard{TYPE}\">
$Clipboard{BODY}</CLIPBOARD>
";
    close(XMLFILE)
}

##############################################################################
sub command_list{

    my $ParamXml = "PARAM.XML";
    my $iSession;
    my $iSection;
    ($iSession, $iSection) = split( /,/, $Editor{SELECT});

    if($Framework){
	my $Section = $SessionRef[$iSession]{SECTION}[$iSection]{NAME};
	if($Section){
	    $ParamXml = "$Section/$CompVersion{$Section}/PARAM.XML";
	}else{
	    $ParamXml = "Param/PARAM.XML";
	}
    }
    my $IsFirstSession = ($iSession == 1);

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
$CommandText" if $Debug =~ /read_command_info/;

}

##############################################################################

sub write_index_html{

    open(FILE, ">$IndexHtmlFile") 
	or die "$ERROR: could not open $IndexHtmlFile\n";
    print FILE
"
<?php Exec('share/Scripts/ParamXmlToHtml.pl') ?>
<#@ `share/Scripts/ParamXmlToHtml.pl -submit='\$FORM{submit}'`; #>
<FRAMESET ROWS=120,1*>
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

    my $DoDebug = ($Debug =~ /write_editor_html/);

    print "Starting write_editor_html\n" if $DoDebug;

    my $SessionSection="    <OPTION VALUE=all>ALL SESSIONS\n";

    my $iSession;
    for $iSession (1..$nSession){
	print "iSession=$iSession\n" if $DoDebug;

	my $SessionRef = $SessionRef[$iSession];
	my $SessionName = $SessionRef->{NAME};
        $SessionSection .= "    <OPTION VALUE=$iSession>$SessionName\n";
	next unless $Framework;
	my $iSection;
	for $iSection (1..$#{$SessionRef->{SECTION}}){
	    print "iSession,iSection=$iSession,$iSection\n" if $DoDebug;
	    my $SectionName = 
		($SessionRef->{SECTION}[$iSection]{NAME} or "CON");
	    $SessionSection .= 
		"      <OPTION VALUE=$iSession,$iSection>"
		.('&nbsp;' x 3)."$SessionName/$SectionName\n";
	}
    }

    # Add SELECTED
    my $Selected = $Editor{SELECT};
    $Selected =~ s/^(\d+,\d+),\d+$/$1/; # Chop off item index if present
    $SessionSection =~ s/(VALUE=$Selected)/$1 SELECTED/;

    print "SessionSection=$SessionSection\n" if $DoDebug;

    my $InsertList;

    if($Selected eq "all"){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SESSION\n"
	    if $Clipboard{TYPE} eq "SESSION";
        $InsertList .= "    <OPTION>NEW SESSION\n";

    }elsif($Selected =~ /^\d+$/ and $Framework){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SECTION\n"
	    if $Clipboard{TYPE} eq "SECTION";
        $InsertList .= "    <OPTION>Section CON\n";
	my $Comp;
	for $Comp (split ',', $ValidComp){
	    $InsertList .= "     <OPTION>Section $Comp\n";
	}
    }else{
        $InsertList  =     "    <OPTION>FILE\n";
	$InsertList .=     "    <OPTION>PASTE COMMAND/COMMENT\n"
	    if $Clipboard{TYPE} =~ /COMMAND|COMMENT|USERINPUT/;
	$InsertList .=     "    <OPTION>COMMENT\n";
    	if($Editor{ABC}){
	    $InsertList .= "    <OPTION>COMMANDS ALPHABETICALLY\n";
	}else{
	    $InsertList .= "    <OPTION>COMMANDS BY GROUP\n";
	}
	$InsertList     .= &command_list;
    }

    # Add SELECTED
    my $Insert = $Editor{INSERT};

    if( not ($InsertList =~s/OPTION>$Insert/OPTION SELECTED>$Insert/)){
	$InsertList =~ s/OPTION>(COMMANDS|NEW|Section CON)/OPTION SELECTED>$1/;
	$Insert = "new";
	$Editor{INSERT} = "new";
    }

    my $InsertItem;

    if($Insert eq "FILE"){
	my @Files;
	@Files = glob("run/PARAM*.in* Param/PARAM*.in.*");

	my $Files;
	$Files = "    <OPTION>".join("\n    <OPTION>",@Files) if @Files;
	$InsertItem=
"  <SELECT NAME=file onChange=\"dynamic_select('editor','file')\">
    <OPTION>SELECT FILE
$Files
  </SELECT>
";
	# Add SELECTED
	my $File = $Editor{FILE};
	$InsertItem =~ s/<OPTION(.*$File\n)/<OPTION SELECTED$1/ if $File;
    }
    # Add checkbox if COMMAND list
    if($InsertList =~ /COMMAND/){
	$InsertItem .= 
"  <INPUT TYPE=CHECKBOX NAME=abc VALUE=1
      onChange=\"parent.location.href='$IndexHtmlFile?submit=ABC_ON'\"
   >abc
";
	if($Editor{ABC}){
	    $InsertItem =~ s/>abc$/ CHECKED>abc/;
	    $InsertItem =~ s/ABC_ON'/ABC_OFF'/;
	}
    }

    chop $SessionSection;
    chop $InsertList;
    chop $InsertItem;

    my $Editor = "
<html>
<head>
  <script language=javascript type=text/javascript>
  <!--
  function dynamic_select(NameForm, NameElement){
    elem = document.forms[NameForm][NameElement];
    parent.location.href = '$IndexHtmlFile?submit=action:select_'
        + NameElement + ';id:' 
        + escape(elem.options[elem.selectedIndex].value);
  }
  // -->
  </script>
</head>
<body bgcolor=WHITE>
  <FORM NAME=editor ACTION=$IndexHtmlFile TARGET=_parent>
  <TABLE WIDTH=100%>
    <TR>
      <TD ALIGN=LEFT WIDTH=20%>
<INPUT TYPE=SUBMIT NAME=submit VALUE=CHECK>
<INPUT TYPE=SUBMIT NAME=submit VALUE=SAVE>
<INPUT TYPE=SUBMIT NAME=submit VALUE=\"SAVE AS\">
      </TD>
      <TD ALIGN=CENTER WIDTH=60%>
         <FONT COLOR=red>
$ParamFile
         </FONT>
      </TD>
      <TD ALIGN=RIGHT WIDTH=20%>
<INPUT TYPE=SUBMIT NAME=submit VALUE=EXIT>
      </TD>
    </TR>
  </TABLE>
<hr>
  <CENTER>
  View: 
  <SELECT NAME=session onChange=\"dynamic_select('editor','session')\">
$SessionSection
  </SELECT>
&nbsp;Insert: 
  <SELECT NAME=insert onChange=\"dynamic_select('editor','insert')\">
$InsertList
  </SELECT>
$InsertItem
  </CENTER>
  </FORM>
</body>
</html>
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
    }elsif($Editor{INSERT} eq "FILE" and -f $Editor{FILE}){
	$Manual = "<H1>$Editor{FILE}</H1>\n<PRE>$Clipboard{BODY}\n</PRE>";
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

    my $Param = "  <head>
  <style type=\"text/css\">
  a {text-decoration: none;}
  </style>
  </head>

  <BODY BGCOLOR=\#DDDDDD TEXT=BLACK LINK=BLUE VLINK=BLUE>
";

    my $MinMaxSessionButton;
    my $MinMaxSectionButton;
    my $MinMaxItemButton;
    my $InsertSessionButton;
    my $InsertSectionButton;
    my $InsertItemButton;
    my $CopySessionButton;
    my $CopySectionButton;
    my $CopyItemButton;
    my $RemoveSessionButton;
    my $RemoveSectionButton;
    my $RemoveItemButton;


    ########################## SESSION #################################
    my $iSession;
    for $iSession (1..$nSession){
	my $SessionView =  $SessionRef[$iSession]{VIEW};
	my $SessionName =  $SessionRef[$iSession]{NAME};
	my $SessionGivenName;
	$SessionGivenName=": $SessionName" unless 
	    $SessionName =~ /^Session \d/;

	my $Action = "A TARGET=_parent ".
	    "HREF=$IndexHtmlFile?submit=id:$iSession;action";

	if($SessionView eq "MIN"){
	    $MinMaxSessionButton = "      <$Action:maximize_session
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize session\"></A>
";
	}else{
	    $MinMaxSessionButton = "      <$Action:minimize_session
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize session\"></A>
";
	}

	$InsertSessionButton = "    <$Action:insert_session
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert session\"></A>
"                 if $Editor{SELECT} eq "all";

	$CopySessionButton = "      <$Action:copy_session
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy session\"></A>
";
	$RemoveSessionButton = "      <$Action:remove_session
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove session\"></A>
";

	# Place anchor to selected session
	$Param .= "<A NAME=SELECTED>anchor</A>\n" if $Editor{SELECT} eq $iSession;

	$Param .=
"<hr COLOR=BLACK>
  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#CCCCCC>
    <TR>
      <TD WIDTH=20 ALIGN=LEFT>
$MinMaxSessionButton
      </TD>
      <TD WIDTH=380><$Action:select_session>
Session $iSession$SessionGivenName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSessionButton$CopySessionButton$RemoveSessionButton
      </TD>
    </TR>
  </TABLE>
";
	next if $SessionView eq "MIN";

	######################## SECTION ###############################

	my $iSection;
	my $nSection = $#{ $SessionRef[$iSession]{SECTION} };
	for $iSection (1..$nSection){
	    my $SectionRef  = $SessionRef[$iSession]{SECTION}[$iSection];
	    my $SectionView = $SectionRef->{VIEW};
	    my $SectionName = ($SectionRef->{NAME} or "CON");

	    my $Action = "A TARGET=_parent ".
		"HREF=$IndexHtmlFile?submit=id:$iSession,$iSection;action";

	    if($SectionView eq "MIN"){
		$MinMaxSectionButton = "      <$Action:maximize_section
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize section\"></A>
";
	    }else{
		$MinMaxSectionButton = "      <$Action:minimize_section
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize section\"></A>
";
	    }

	    $InsertSectionButton = "    <$Action:insert_section
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert section\"></A>
" 	    if $Editor{SELECT} =~ /^\d+$/;

	    $CopySectionButton = "      <$Action:copy_section
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy section\"></A>
";
	    $RemoveSectionButton = "      <$Action:remove_section
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove section\"></A>
";

	    my $InsertItemButton;

	    # Place anchor to selected section
	    $Param .= "<A NAME=SELECTED>anchor</A>\n" 
		if $Editor{SELECT} eq "$iSession,$iSection";

	    $Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#CCCCCC>
    <TR>
      <TD WIDTH=20>
      </TD>
      <TD COLSPAN=2>
<hr COLOR=GREY>
      </TD>
    </TR>
    <TR>
      <TD ALIGN=CENTER>
$MinMaxSectionButton
      </TD>
      <TD WIDTH=380><$Action:select_section>
Section: $SectionName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton$CopySectionButton$RemoveSectionButton
    </TR>
  </TABLE>
";
	    next if $SectionView eq "MIN";
	    
###################### ITEM LOOP ############################################

	    my $Action = "A TARGET=_parent HREF=$IndexHtmlFile?submit=".
		"id:$iSession,$iSection,0;action";

	    if($SectionName eq $Clipboard{SECTION} and 
	       $Clipboard{TYPE} =~ /COMMAND|COMMENT|USERINPUT/ ){
		$InsertItemButton = "    <$Action:insert_item
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert item\"></A>
";
	    }else{
		$InsertItemButton="";
	    }

	    my $iItem;
	    my $nItem = $#{ $SectionRef->{ITEM} };
	    for $iItem (1..$nItem){

		my $ItemRef  = $SectionRef->{ITEM}[$iItem];
		my $ItemView = $ItemRef->{VIEW};
		my $ItemType = lc($ItemRef->{TYPE});
		my $ItemHead = $ItemRef->{HEAD};
		my $ItemTail = $ItemRef->{TAIL};

		my $TableColor = $TableColor{$ItemType};

		$Action =~ s/\d+;action$/$iItem;action/;
		$InsertItemButton =~ s/\d+;action:/$iItem;action:/;

		if($ItemType eq "userinput"){
		    # Fix the content of the user input
		    $ItemHead = "\#USERINPUT";
		    $ItemTail =~ s/\n*\#USERINPUTEND.*//m;
		}
		my $nLine = ($ItemTail =~ s/\n/\n/g);

		# Place anchor to selected item
		$Param .= "<A NAME=SELECTED>anchor</A>\n" 
		    if $Editor{SELECT} eq "$iSession,$iSection,$iItem";

		if($ItemView eq "EDIT"){
		    $nLine += 2;
		    if($ItemType eq "comment"){
			$ItemTail = $ItemHead.$ItemTail;
			$ItemHead = "";
			$nLine++;
		    }
		    $Param .= "
  <FORM TARGET=_parent>
  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#BBEEFF>
    <TR>
      <TD WIDTH=20>
      </TD>
      <TD WIDTH=380><FONT COLOR=BLUE>
$ItemHead
      </FONT></TD>
      <TD ALIGN=RIGHT>
        <$Action:save_item_edit
><IMG SRC=$ImageDir/button_save.gif TITLE=\"Save item\"></A>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <$Action:cancel_item_edit
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Cancel\"></A>
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
		    next; # done with editor
		}

		# Create buttons
		if($ItemTail){
		    if($ItemView eq "MIN"){
			$MinMaxItemButton = "      <$Action:maximize_item
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize $ItemType\"></A>
";
		    }else{
			$MinMaxItemButton = "      <$Action:minimize_item
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize $ItemType\"></A>
";
		    }
		}else{
		    $MinMaxItemButton = "";
		}

		$CopyItemButton = "      <$Action:copy_item
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy $ItemType\"></A>
";
		$RemoveItemButton = "      <$Action:remove_item
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove $ItemType\"></A>
";

		# Show first line with usual buttons
		$ItemHead = "<PRE>$ItemHead</PRE>" if $ItemType eq "comment";

		$Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=20 ALIGN=RIGHT>
$MinMaxItemButton
      </TD>
      <TD WIDTH=380><$Action:edit_item>
$ItemHead
      </A></TD>
      <TD ALIGN=RIGHT ALIGN=TOP>
$InsertItemButton$CopyItemButton$RemoveItemButton
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
      <TD WIDTH=380 COLSPAN=2><$Action:edit_item>
<PRE>$ItemTail</PRE>
      </A></TD>
    </TR>
";
		    if($ItemType eq "userinput"){
			$MinMaxItemButton =~ s/minimize\.gif/minimize_up.gif/;
			$Param .= 
"    <TR>
      <TD WIDTH=20 ALIGN=RIGHT>
$MinMaxItemButton
      </TD>
      <TD WIDTH=380><$Action:edit_item>
\#USERINPUT
      </A></TD>
      <TD ALIGN=RIGHT>
$CopyItemButton$RemoveItemButton
      </TD>
    </TR>
";
		    }
		    $Param .= 
"  </TABLE>
";
		}else{  #command type item
                    $Param .= 
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#CCCCCC>
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

		    } # end item body line loop

		    $Param .= 
"  </TABLE>
";
		} #endif item type

		$Param .= "<p>\n" if $ItemType ne "comment";

	    } # end item loop

	    ###### End section #########

	    $iItem=$nItem+1; $iItem=1 if $iItem==0;
	    $InsertItemButton =~ s/id:[\d,]+/id:$iSession,$iSection,$iItem/;
	    $MinMaxSectionButton =~ s/minimize\.gif/minimize_up.gif/;

	    $Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#CCCCCC>
    <TR>
      <TD WIDTH=20 ALIGN=CENTER>
$MinMaxSectionButton
      </TD>
      <TD WIDTH=380><$Action:select_section>
Section: $SectionName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertItemButton$CopySectionButton$RemoveSectionButton
      </TD>
    </TR>
  </TABLE>
";
	} # Section loop

	###### End session #########

	$InsertSectionButton =~ s/id:[\d,]+/id:$iSession,$iSection/;
	$MinMaxSessionButton =~ s/minimize\.gif/minimize_up.gif/;
	$Param .=
"  <TABLE BORDER=0 WIDTH=100% BGCOLOR=\#CCCCCC>
    <TR>
      <TD>
      </TD>
      <TD COLSPAN=2>
<hr COLOR=GREY>
      </TD>
    </TR>
    <TR>
      <TD WIDTH=20 ALIGN=LEFT>
$MinMaxSessionButton
      </TD>
      <TD WIDTH=380><$Action:select_session>
Session $iSession$SessionGivenName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton$CopySessionButton$RemoveSessionButton
      </TD>
    </TR>
  </TABLE>
";
    } # session loop

    $iSession=$nSession+1;
    $InsertSessionButton =~ s/id:[\d,]+/id:$iSession/;

    $Param .= 
"<hr COLOR=BLACK>\n
  <TABLE WIDTH=100% BGCOLOR=\#CCCCCC><TR><TD ALIGN=RIGHT>
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
