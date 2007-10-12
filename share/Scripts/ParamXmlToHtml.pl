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

# Print help message and exit if -h switch was used
&print_help if $Help;

# Error string
my $ERROR = 'ParamXmlToHtml_ERROR';

# Fixed file names
my $IndexPhpFile   = "index.php";
my $ParamHtmlFile  = "param.html";
my $EditorHtmlFile = "editor.html";
my $ManualHtmlFile = "manual.html";
my $JumpHtmlFile   = "jump.html";
my $ConfigFile     = "$ENV{HOME}/ParamEditor.conf";
my $ImageDir       = "share/Scripts";

# Variables that can be modified in the included config file
our $RedoFrames       =0;           # rewrite index.php (jump.html)
our $DoSafariJumpFix  =0;           # work around the Safari bug
our $FrameHeights     ='15%,85%';   # heights of top frame and lower frames
our $FrameWidths      ='60%,40%';   # widths of the left and right frames
our $TopBgColor       ='#DDDDDD';   # background color for the top frame
our $TopFileNameFont  ='COLOR=RED'; # font used for file name in top frame
our $FileNameEditorWidth = 40;      # width (chars) of filename input box
our $TopTableWidth    ='100%'   ;  # width of stuff above the line in top frame
our $TopLine          ='<HR>'   ;  # separator line in top frame

our $RightBgColor     ='WHITE';    # background color for the right side frame

our $LeftTableWidth   ='100%'   ;  # width of table in left frame
our $LeftColumn1Width =  '30'   ;  # width of minimize/maximize button column
our $LeftColumn2Width = '380'   ;  # width of column with commands/comments
our $LeftBgColor      ='#DDDDDD';  # background color for the left side frame
our $SessionBgColor   ='#CCCCCC';  # background color for session markers
our $SessionLine ='<HR COLOR=BLACK NOSHADE SIZE=4>'; # session separator line
our $SessionEditorSize=30;         # size (characters) for session name input
our $SectionBgColor   ='#CCCCCC';  # background color for section markers
our $SectionLine ='<HR COLOR=GREY NOSHADE SIZE=2>'; # section separator line
our $SectionColumn1Width=8      ;  # The space before the section marker
our $ItemEditorBgColor='#BBEEFF';  # background color command/comment editor
our $ItemEditorWidth  =60;         # width (chars) for command/comment editor
our $CommandBgColor   ='#CCCCCC';  # background color for commands
our $ParameterBgColor ='#CCCCCC';  # background color for command parameters
our $CommentBgColor   ='#DDDDDD';  # background color for comments
our $UserInputBgColor ='#CCCCCC';  # background color for user input commands

# Allow user to modify defaults
do $ConfigFile;

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

my $CheckResult;   # Result from TestParam.pl script

my %TableColor = ("command"   => $CommandBgColor,
		  "comment"   => $CommentBgColor,
		  "userinput" => $UserInputBgColor);


# Set name of parameter file
my $ParamFile = "run/PARAM.in"; 
if($ARGV[0]){
    $ParamFile = $ARGV[0];
}elsif(open(FILE, $ParamHtmlFile)){
    while($_ = <FILE>){
	next unless /<TITLE>([^\<]+)/;
	$ParamFile = $1;
	last;
    }
    close(FILE);
}
my $XmlFile   = "$ParamFile.xml"; # Name of the XML enhanced parameter file
my $TextFile  = "$ParamFile.txt"; # Name of the temporary text parameter file

if($ParamFile){
    `share/Scripts/ParamTextToXml.pl $ParamFile` unless -f $XmlFile;
    &read_xml_file;
}

open(PARAM, ">$ParamHtmlFile") 
    or die "$ERROR: could not open output file $ParamHtmlFile\n";
open(EDITOR, ">$EditorHtmlFile") 
    or die "$ERROR: could not open output file $EditorHtmlFile\n";
open(MANUAL, ">$ManualHtmlFile") 
    or die "$ERROR: could not open output file $ManualHtmlFile\n";

&modify_xml_data if $Submit;

&write_index_php if $RedoFrames or not -f $IndexPhpFile;

&write_jump_html if $DoSafariJumpFix and ($RedoFrames or not -f $JumpHtmlFile);

&write_editor_html;

&write_manual_html;

&write_param_html;

&write_xml_file if $Submit;

exit 0;

##############################################################################

sub read_xml_file{

    open(XMLFILE, $XmlFile) 
	or die "$ERROR: could not open input file $XmlFile\n";

    $nSession   = 0;
    @SessionRef = ();
    %Editor     = ();
    %Clipboard  = ();
    my $iSection;
    my $iItem;

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
	    my $Name=$1;
	    my $View=$2;
	    $SessionRef[$nSession]{NAME} = $Name;
	    $SessionRef[$nSession]{VIEW} = $View;
	    print "nSession=$nSession Name=$Name View=$View\n" if $Debug;
	}elsif(/<SECTION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $iSection++;
	    $iItem=0;
	    my $Name=$1;
	    my $View = $2;
	    $SessionRef[$nSession]{SECTION}[$iSection]{VIEW} = $View;
	    $SessionRef[$nSession]{SECTION}[$iSection]{NAME} = $Name;

	    print "Section=$iSection name=$Name View=$View\n" if $Debug;

	    # Do not read details if the section is minimized
	    #if(not $Submit 
	    #   and ($SectionView eq "MIN" or $SessionView eq "MIN")){
	    #   $_=<XMLFILE> while not /<\/SECTION>/;
	    #}		
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

    my %Form;
    %Form = split (/[:;]/, $Submit);
    #warn "Submit=$Submit\n";
    #warn join(', ',%Form),"\n";

    if($_ = $Form{submit}){
	if( /^CHECK\b/ ){
	    `share/Scripts/ParamXmlToText.pl $XmlFile`;
	    my $TestScript = ($Framework ? "Scripts" : ".")."/TestParam.pl";
	    $CheckResult = `$TestScript $TextFile 2>&1`;
	    $CheckResult = "No errors found" unless $CheckResult;
	}elsif( /^SAVE$/ ){
	    `share/Scripts/ParamXmlToText.pl $XmlFile`;
	    `cp $TextFile $ParamFile`;
	}elsif( /^SAVE AS$/ ){
	    $Form{FILENAME} =~ s/\.(xml|txt)$//;
	    if(open(FILE,">$Form{FILENAME}")){
		close(FILE);
		$ParamFile = $Form{FILENAME};
		rename($XmlFile, "$ParamFile.xml");
		$XmlFile   = "$ParamFile.xml";
		$TextFile  = "$ParamFile.txt";
		`share/Scripts/ParamXmlToText.pl $XmlFile`;
		`cp $TextFile $ParamFile`;
	    }else{
		$Editor{READFILENAME}="SAVE AS";
		$Editor{NEWFILENAME}="Could not open $Form{FILENAME}"
		    if $Form{FILENAME};
	    }
	}elsif( /^OPEN$/ ){
	    $Form{FILENAME} =~ s/\.(xml|txt)$//;
	    if(open(FILE,"$Form{FILENAME}")){
		close(FILE);
		# Make a safety save
		`share/Scripts/ParamXmlToText.pl $XmlFile` if -f $XmlFile;
		$ParamFile = $Form{FILENAME};
		$XmlFile   = "$ParamFile.xml";
                $TextFile  = "$ParamFile.txt";
		`share/Scripts/ParamTextToXml.pl $ParamFile`;
		&read_xml_file;
		warn "command1=$SessionRef[1]{SECTION}[1]{ITEM}[1]{HEAD}\n";
	    }else{
		$Editor{READFILENAME}="OPEN";
		$Editor{NEWFILENAME}="Could not open $Form{FILENAME}"
		    if $Form{FILENAME};
	    }
	}elsif( /^CANCEL$/ ){
	    # Do nothing
	}elsif( /^SAVE AND EXIT$/ ){
	    # save the file then kill the job
	    `share/Scripts/ParamXmlToText.pl $XmlFile` if -f $XmlFile;
	    `cp $TextFile $ParamFile`;
	    kill(-9, getpgrp);
	}elsif( /^EXIT$/ ){
	    # make a safety save then kill the job
	    `share/Scripts/ParamXmlToText.pl $XmlFile`;
	    kill(-9, getpgrp);
	}elsif( /^ABC_ON$/ ){
	    $Editor{ABC}=1;
	}elsif( /^ABC_OFF$/ ){
	    $Editor{ABC}=0;
	}elsif( /^SAVE SESSION NAME$/ ){
	    my $iSession = $Form{id};
	    $SessionRef[$iSession]{VIEW}="MAX";
	    $SessionRef[$iSession]{NAME}=$Form{name};
	    $Editor{SELECT}=$iSession;
	}
    }elsif($_ = $Form{action}){
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
	
	# The following actions
	my $iSession;
	my $iSection;
	my $iItem;
	($iSession,$iSection,$iItem) = split(/,/,$id);
	
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
	    $NewSessionRef->{VIEW} = "EDIT";
	    $NewSessionRef->{NAME} = "";
	    $NewSessionRef->{SECTION}[1]{VIEW}="MAX";
	    if($Editor{INSERT} eq "PASTE SESSION"){
		$NewSessionRef->{SECTION}[1]{ITEM}[1]{HEAD} = $Clipboard{BODY};
	    }elsif($Editor{INSERT} eq "NEW SESSION"){
		$NewSessionRef->{SECTION}[1]{ITEM}[1]{HEAD} = "New session\n";
	    }
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{TYPE} = "COMMENT";
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{VIEW} = "MAX";
	    splice (@SessionRef, $iSession, 0, $NewSessionRef);
	    $nSession++;
	}elsif( /insert_section/ ){
	    my $NewSectionRef;
	    $NewSectionRef->{VIEW} = "MAX";
	    if($Editor{INSERT} eq "PASTE SECTION"){
		$NewSectionRef->{NAME} = $Clipboard{SECTION};
		$NewSectionRef->{ITEM}[1]{HEAD} = $Clipboard{BODY};
	    }elsif($Editor{INSERT} =~ /Section (\w+)/){
		$NewSectionRef->{NAME} = $1;
		$NewSectionRef->{ITEM}[1]{VIEW} = "MAX";
		$NewSectionRef->{ITEM}[1]{TYPE} = "COMMENT";
		$NewSectionRef->{ITEM}[1]{HEAD} = "New $1 section\n";
	    }
	    splice (@{$SessionRef->{SECTION}}, $iSection, 0, $NewSectionRef);
	}elsif( /insert_item/ ){
	    my $NewItemRef;
	    $NewItemRef->{TYPE} = $Clipboard{TYPE};
	    if($Editor{INSERT} =~ /^PASTE/){
		$NewItemRef->{TYPE} = $Clipboard{TYPE};
		$NewItemRef->{VIEW} = "MAX";
		($NewItemRef->{HEAD}, $NewItemRef->{TAIL}) 
		    = split(/\n/, $Clipboard{BODY}, 2);
		$NewItemRef->{HEAD}.="\n";
	    }elsif($Editor{INSERT} =~ /^COMMENT/){
		$NewItemRef->{TYPE} = "COMMENT";
		$NewItemRef->{VIEW} = "EDIT";
	    }elsif($Editor{INSERT} =~ /^\#USERINPUTBEGIN/){
		$NewItemRef->{TYPE} = "USERINPUT";
		$NewItemRef->{VIEW} = "EDIT";
		$NewItemRef->{HEAD} = "\#USERINPUT";
	    }elsif($Editor{INSERT} =~ /^\#/){
		$NewItemRef->{TYPE} = "COMMAND";
		$NewItemRef->{VIEW} = "EDIT";
		$NewItemRef->{HEAD} = "$Editor{INSERT}\n";
		$NewItemRef->{TAIL} = "CommandExample\n";
	    }else{
		# invalid choice like "COMMANDS BY GROUP";
		return;
	    }
	    splice (@{$SectionRef->{ITEM}}, $iItem, 0, $NewItemRef);
	}elsif( /Save command|Save userinput|Save comment/ ){
	    $ItemRef->{VIEW}="MAX";
	    if(/Save comment/){
		($ItemRef->{HEAD}, $ItemRef->{TAIL})=split(/\n/,$Form{text},2);
		$ItemRef->{HEAD} .= "\n";
	    }else{
		$ItemRef->{TAIL}=$Form{text};
	    }
	    $ItemRef->{TAIL} =~ s/\s*\n/\n/g;    # clean up line endings
	    $ItemRef->{TAIL} =~ s/\n\n+/\n\n/g;  # remove multiple empty lines
	    $ItemRef->{TAIL} =~ s/\n*$/\n/;      # close with a single \n
	}elsif( /CANCEL/ ){
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

sub write_index_php{

    open(FILE, ">$IndexPhpFile") 
	or die "$ERROR: could not open $IndexPhpFile\n";

    my $ParamLink = ($DoSafariJumpFix ? $JumpHtmlFile : "$ParamHtmlFile#HERE");

    print FILE
"
<?php Exec('share/Scripts/ParamXmlToHtml.pl') ?>
<#@ \$form=join(';',%FORM);                             #>
<#@ `share/Scripts/ParamXmlToHtml.pl -submit='\$form'`; #>
<FRAMESET ROWS=$FrameHeights>
  <FRAME SRC=\"./$EditorHtmlFile\" NAME=EDITOR FRAMEBORDER=1 SCROLLING=no 
                                                                   NORESIZE>
  <FRAMESET COLS=$FrameWidths>
    <FRAME SRC=\"./$ParamLink\"  NAME=PARAMFILE FRAMEBORDER=1 NORESIZE>
    <FRAME SRC=\"./$ManualHtmlFile\"  NAME=MANUAL FRAMEBORDER=1 NORESIZE>
  </FRAMESET>
</FRAMESET>
";
    close(FILE);
}

##############################################################################

sub write_jump_html{

    open(FILE, ">$JumpHtmlFile")
	or die "$ERROR: could not open $JumpHtmlFile\n";

    print FILE
"<HEAD> 
<SCRIPT language=JavaScript> 
   window.location='$ParamHtmlFile#HERE';
</SCRIPT>
</HEAD>
<BODY BGCOLOR=$LeftBgColor>
</BODY>
";
    close(FILE);
}

##############################################################################

sub write_editor_html{

    my $DoDebug = ($Debug =~ /write_editor_html/);

    print "Starting write_editor_html\n" if $DoDebug;

    my $EditButtons;

    if(not -f $ParamFile){
	print EDITOR "
<BODY BGCOLOR=$TopBgColor>
  <FORM NAME=editor ACTION=$IndexPhpFile TARGET=_parent>
    <CENTER>
        <FONT COLOR=RED>Please provide name of input parameter file:</RED><BR>
        <INPUT TYPE=TEXT SIZE=$FileNameEditorWidth NAME=FILENAME>
        <INPUT TYPE=SUBMIT NAME=submit VALUE=OPEN>&nbsp;&nbsp;&nbsp;&nbsp;
        <INPUT TYPE=SUBMIT NAME=submit VALUE=EXIT>
    </CENTER>
  </FORM>
</BODY>
";
	close(EDITOR);
	exit 0;
    }

    if($Editor{READFILENAME}){
	$EditButtons = 
"      <TD COLSPAN=2 ALIGN=CENTER>
        <FONT $TopFileNameFont>New parameter file name:</FONT>
        <INPUT TYPE=TEXT SIZE=$FileNameEditorWidth
	      NAME=FILENAME VALUE=\"$Editor{NEWFILENAME}\">
        <INPUT TYPE=SUBMIT NAME=submit VALUE=\"$Editor{READFILENAME}\">
        &nbsp
        <INPUT TYPE=SUBMIT NAME=submit VALUE=CANCEL>
      </TD>
      <TD>
        <INPUT TYPE=SUBMIT NAME=submit VALUE=EXIT>
      </TD>
";
    }else{
	$EditButtons = 
"      <TD ALIGN=LEFT>
<INPUT TYPE=SUBMIT NAME=submit VALUE=CHECK>
<INPUT TYPE=SUBMIT NAME=submit VALUE=SAVE>
<INPUT TYPE=SUBMIT NAME=submit VALUE=\"SAVE AS\">
<INPUT TYPE=SUBMIT NAME=submit VALUE=OPEN>
      </TD>
      <TD ALIGN=CENTER>
<FONT $TopFileNameFont>$ParamFile</FONT>
      </TD>
      <TD ALIGN=RIGHT>
<INPUT TYPE=SUBMIT NAME=submit VALUE=\"SAVE AND EXIT\">&nbsp;&nbsp;&nbsp;
<INPUT TYPE=SUBMIT NAME=submit VALUE=EXIT>
      </TD>
";
    }
    my $SessionSection="    <OPTION VALUE=all>ALL SESSIONS\n";

    my $iSession;
    for $iSession (1..$nSession){
	print "iSession=$iSession\n" if $DoDebug;

	my $SessionRef = $SessionRef[$iSession];
	my $SessionName = ($SessionRef->{NAME} or "Session $iSession");
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
      onChange=\"parent.location.href='$IndexPhpFile?submit=ABC_ON'\"
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
<HTML>
<HEAD>
  <SCRIPT LANGUAGE=javascript TYPE=text/javascript>
  <!--
  function dynamic_select(NameForm, NameElement){
    elem = document.forms[NameForm][NameElement];
    parent.location.href = '$IndexPhpFile?action=select_'
        + NameElement + '&id=' 
        + escape(elem.options[elem.selectedIndex].value);
  }
  // -->
  </SCRIPT>
  <BASE TARGET=_parent>
</HEAD>
<BODY BGCOLOR=$TopBgColor>
  <FORM NAME=editor ACTION=$IndexPhpFile>
  <TABLE WIDTH=$TopTableWidth>
    <TR>
$EditButtons
    </TR>
    <TR>
      <TD COLSPAN=3>
$TopLine
      </TD>
    </TR>
    <TR>
      <TD COLSPAN=3 ALIGN=CENTER>
  View: 
  <SELECT NAME=session onChange=\"dynamic_select('editor','session')\">
$SessionSection
  </SELECT>
&nbsp;Insert: 
  <SELECT NAME=insert onChange=\"dynamic_select('editor','insert')\">
$InsertList
  </SELECT>
$InsertItem
      </TD>
    </TR>
  </TABLE>
  </FORM>
</BODY>
</HTML>
";
    print EDITOR $Editor;
    close EDITOR;
}

##############################################################################

sub write_manual_html{

    my $Manual;
    if($CheckResult){
	$Manual = "<H1>Checking for Errors</H1>\n".
	    "<FONT COLOR=RED><PRE>$CheckResult</PRE></FONT>\n";
	if($CheckResult !~ /no errors/i and open(FILE, $TextFile)){
	    my @TextFile;
	    @TextFile = <FILE>;
	    close(FILE);
	    my $iError = 1;
	    while($CheckResult =~ 
		  s[Error at line (\d+)(.*)\n((\t.*\n)+)]
		  [<A HREF=\#ERROR$iError>ERROR at line $1$2</A>\n$3]){
		my $iLine = $1;
		my $Error = $3; $Error =~ s/\t//g; chop $Error;
		$iLine-- while $iLine > 0 and $TextFile[$iLine] !~ /^\#/;
		$TextFile[$iLine] = 
		    "<A NAME=ERROR$iError><FONT COLOR=RED>$Error</FONT></A>".
		    "\n$TextFile[$iLine]";
		$iError++;
	    }
	    $Manual .= "<PRE>\n". join('', @TextFile) . "</PRE>\n";
	}
    }elsif($CommandExample){
	$Manual =  $CommandText;
	$Manual =~ s/\n\n/\n<p>\n/g;
	$Manual =  "<H1>Manual</H1>\n<PRE>\n$CommandExample\n</PRE>\n$Manual";
	$Manual =~ s/\n$//;
    }elsif($Editor{INSERT} =~ /^PASTE/){
	$Manual = "<H1>Clipboard</H1>\n<PRE>\n$Clipboard{BODY}\n</PRE>";
    }elsif($Editor{INSERT} eq "FILE" and -f $Editor{FILE}){
	$Manual = "<H1>$Editor{FILE}</H1>\n<PRE>$Clipboard{BODY}\n</PRE>";
    }

    print MANUAL
"<BODY BGCOLOR=$RightBgColor>
$Manual
</BODY>
";
    close MANUAL;
}

##############################################################################

sub write_param_html{

    my $SelectedSectionName;
    if($Editor{SELECT} =~ /^(\d+),(\d+)/){
	$SelectedSectionName = ($SessionRef[$1]{SECTION}[$2]{NAME} or "CON");
    }

    my $Param = 
"  <HEAD>
    <TITLE>$ParamFile</TITLE>
    <STYLE TYPE=\"text/css\">
      A {text-decoration: none;}
      A IMG {border: none;}
    </STYLE>
    <BASE TARGET=_parent>
  </HEAD>
  <BODY BGCOLOR=$LeftBgColor TEXT=BLACK LINK=BLUE VLINK=BLUE>
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

	my $Action = "A HREF=$IndexPhpFile?id=$iSession\&action";

	my $SessionTagTop;
	my $SessionTagBot;

	if($SessionView eq "EDIT"){
	    $SessionTagTop = "
<FORM NAME=action ACTION=$IndexPhpFile>
Session $iSession: 
<INPUT NAME=name TYPE=TEXT SIZE=$SessionEditorSize VALUE=\"$SessionName\">
<INPUT NAME=id TYPE=HIDDEN VALUE=$iSession>
<INPUT TYPE=SUBMIT NAME=submit VALUE=\"SAVE SESSION NAME\"></FORM>";
	    $SessionTagBot = "Session $iSession";
	}else{
	    if($SessionName){
		$SessionTagTop = "<$Action=select_session 
                    TITLE=\"Select session\">Session $iSession</A>:\&nbsp;".
		    "<$Action=edit_session 
                    TITLE=\"Edit session name\">$SessionName</A>";
	    }else{
		$SessionTagTop = "<$Action=edit_session TITLE=\"edit session\">Session $iSession</A>";
	    }
	    $SessionTagBot = $SessionTagTop;
	}

	if($SessionView eq "MIN"){
	    $MinMaxSessionButton = "      <$Action=maximize_session
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize session\"></A>
";
	}else{
	    $MinMaxSessionButton = "      <$Action=minimize_session
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize session\"></A>
";
	}

	$InsertSessionButton = "    <$Action=insert_session
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert session\"></A>
"                 if $Editor{SELECT} eq "all";

	$CopySessionButton = "      <$Action=copy_session
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy session\"></A>
";
	$RemoveSessionButton = "      <$Action=remove_session
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove session\"></A>
";

	# Place anchor to selected session
	$Param .= "<DIV><A NAME=HERE></A></DIV>\n" 
	    if $Editor{SELECT} eq $iSession;

	$Param .=
"$SessionLine
  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$SessionBgColor>
    <TR>
      <TD ALIGN=LEFT>
$MinMaxSessionButton
$SessionTagTop
      </TD>
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

	    my $Action = "A HREF=$IndexPhpFile?id=$iSession,$iSection\&action";

	    if($SectionView eq "MIN"){
		$MinMaxSectionButton = "      <$Action=maximize_section
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize section\"></A>
";
	    }else{
		$MinMaxSectionButton = "      <$Action=minimize_section
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize section\"></A>
";
	    }

	    $InsertSectionButton = "    <$Action=insert_section
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert section\"></A>
" 	    if $Editor{SELECT} =~ /^\d+$/;

	    $CopySectionButton = "      <$Action=copy_section
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy section\"></A>
";
	    $RemoveSectionButton = "      <$Action=remove_section
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove section\"></A>
";

	    # Place anchor to selected section
	    $Param .= "<DIV><A NAME=HERE></A></DIV>\n" 
		if $Editor{SELECT} eq "$iSession,$iSection";

	    $Param .=
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$SectionBgColor>
    <TR>
      <TD WIDTH=$SectionColumn1Width>
      </TD>
      <TD COLSPAN=2>
$SectionLine
      </TD>
    </TR>
    <TR>
      <TD WIDTH=$SectionColumn1Width>
      </TD>
      <TD ALIGN=LEFT>
$MinMaxSectionButton
      <$Action=select_section>
Section $SectionName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton$CopySectionButton$RemoveSectionButton
    </TR>
  </TABLE>
";
	    next if $SectionView eq "MIN";
	    
###################### ITEM LOOP ############################################

	    my $Action = 
		"A HREF=$IndexPhpFile?id=$iSession,$iSection,0\&action";

	    if($SelectedSectionName eq $SectionName){
#	    if($SectionName eq $Clipboard{SECTION} and 
#	       $Clipboard{TYPE} =~ /COMMAND|COMMENT|USERINPUT/ ){
		$InsertItemButton = "    <$Action=insert_item
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

		$Action =~ s/\d+\&action$/$iItem\&action/;
		$InsertItemButton =~ s/\d+\&action=/$iItem\&action=/;

		if($ItemType eq "userinput"){
		    # Fix the content of the user input
		    $ItemHead = "\#USERINPUT";
		    $ItemTail =~ s/\n*\#USERINPUTEND.*//m;
		}
		my $nLine = ($ItemTail =~ s/\n/\n/g);

		# Place anchor to selected item
		$Param .= "<DIV><A NAME=HERE></A></DIV>\n"  
		    if $Editor{SELECT} eq "$iSession,$iSection,$iItem";

		if($ItemView eq "EDIT"){
		    $nLine += 2;
		    if($ItemType eq "comment"){
			$ItemTail = $ItemHead.$ItemTail;
			$ItemHead = "";
		    }elsif($ItemTail =~ /CommandExample/){
			($ItemHead,$ItemTail) = split(/\n/,$CommandExample,2);
		    }
		    $nLine = ($ItemTail =~ s/\n/\n/g) + 2;

		    $Param .= "
  <FORM NAME=item_editor ACTION=$IndexPhpFile>
  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$ItemEditorBgColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width>
      </TD>
      <TD WIDTH=$LeftColumn2Width><FONT COLOR=BLUE>
$ItemHead
      </FONT></TD>
      <TD ALIGN=RIGHT>
<INPUT TYPE=SUBMIT name=action value=\"Save $ItemType\">
<INPUT TYPE=SUBMIT name=action value=CANCEL>
<INPUT TYPE=HIDDEN name=id value=$iSession,$iSection,$iItem>
      </TD>
    </TR>
    <TR>
      <TD>
      </TD>
       <TD COLSPAN=2>
<TEXTAREA NAME=text COLS=$ItemEditorWidth ROWS=$nLine>
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
			$MinMaxItemButton = "      <$Action=maximize_item
><IMG SRC=$ImageDir/button_maximize.gif TITLE=\"Maximize $ItemType\"></A>
";
		    }else{
			$MinMaxItemButton = "      <$Action=minimize_item
><IMG SRC=$ImageDir/button_minimize.gif TITLE=\"Minimize $ItemType\"></A>
";
		    }
		}else{
		    $MinMaxItemButton = "";
		}

		$CopyItemButton = "      <$Action=copy_item
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy $ItemType\"></A>
";
		$RemoveItemButton = "      <$Action=remove_item
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove $ItemType\"></A>
";

		# Show first line with usual buttons
		$ItemHead = "<PRE>$ItemHead</PRE>" if $ItemType eq "comment";

		$Param .=
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width ALIGN=RIGHT>
$MinMaxItemButton
      </TD>
      <TD WIDTH=$LeftColumn2Width><$Action=edit_item>
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
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width>
      </TD>
      <TD WIDTH=$LeftColumn2Width COLSPAN=2><$Action=edit_item>
<PRE>$ItemTail</PRE>
      </A></TD>
    </TR>
";
		    if($ItemType eq "userinput"){
			$MinMaxItemButton =~ s/minimize\.gif/minimize_up.gif/;
			$Param .= 
"    <TR>
      <TD WIDTH=$LeftColumn1Width ALIGN=RIGHT>
$MinMaxItemButton
      </TD>
      <TD WIDTH=$LeftColumn2Width><$Action=edit_item>
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
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$ParameterBgColor>
";
		    my $iLine;
		    for $iLine (0..$nLine-1){
			$ItemTail =~ s/(.*)\n//;
			my $Value;
			my $Comment;
			($Value,$Comment) = split(/\t|\s\s\s+/, $1, 2);

			$Param .= 
"    <TR>
      <TD WIDTH=$LeftColumn1Width>
      </TD>
      <TD WIDTH=$LeftColumn2Width>
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
	    $InsertItemButton =~ s/id=[\d,]+/id=$iSession,$iSection,$iItem/;
	    $MinMaxSectionButton =~ s/minimize\.gif/minimize_up.gif/;

	    $Param .=
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$SectionBgColor>
    <TR>
      <TD WIDTH=$SectionColumn1Width>
      </TD>
      <TD ALIGN=LEFT>
$MinMaxSectionButton
      <$Action=select_section>
Section $SectionName
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertItemButton$CopySectionButton$RemoveSectionButton
      </TD>
    </TR>
  </TABLE>
";
	} # Section loop

	###### End session #########

	$iSection = $nSection+1;
	$InsertSectionButton =~ s/id=[\d,]+/id=$iSession,$iSection/;
	$MinMaxSessionButton =~ s/minimize\.gif/minimize_up.gif/;
	$Param .=
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$SectionBgColor>
    <TR>
      <TD WIDTH=$SectionColumn1Width>
      </TD>
      <TD>
$SectionLine
      </TD>
    </TR>
  </TABLE>
  <TABLE WIDTH=$LeftTableWidth BGCOLOR=$SessionBgColor>
    <TR>
      <TD ALIGN=LEFT>
$MinMaxSessionButton
      <$Action=select_session>
$SessionTagBot
      </A></TD>
      <TD ALIGN=RIGHT>
$InsertSectionButton$CopySessionButton$RemoveSessionButton
      </TD>
    </TR>
  </TABLE>
";
    } # session loop

    $iSession=$nSession+1;
    $InsertSessionButton =~ s/id=[\d,]+/id=$iSession/;

    $Param .= "$SessionLine\n";
    $Param .=
"  <TABLE WIDTH=$LeftTableWidth BGCOLOR=$LeftBgColor><TR><TD ALIGN=RIGHT>
$InsertSessionButton
  </TD></TR></TABLE>
"       if $InsertSessionButton;

    $Param .= "</BODY>\n";

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

     Many parameters of the HTML page can be customized by creating a 
     ParamEditor.conf file in the home directory. Type 

grep '^our' share/Scripts/ParamXmlToHtml.pl

     to see the variables that can be modified using the same syntax, 
     but different values. Note the use of semicolons at the end of lines.

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
