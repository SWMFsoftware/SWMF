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
&print_help if $Help or $#ARGV > 1;

# Error string
my $ERROR = 'ParamConvert_ERROR';

# If there are two arguments do conversion between two files
&convert_type if $#ARGV == 1;

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
our $ErrorBgColor     ='RED';      # background color for errors
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

my $JumpHere;      # location ID for the param.html#HERE anchor
my $Here = "<A NAME=HERE><DIV> </DIV></A>\n"; # The HERE anchor

my %TableColor = ("command"   => $CommandBgColor,
		  "comment"   => $CommentBgColor,
		  "error"     => $ErrorBgColor,
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
my $XmlFile  = "$ParamFile.xml"; # Name of the XML enhanced parameter file
my $BakFile  = "$ParamFile.bak"; # Name of the backup text parameter file

if($ARGV[0] or not -f $XmlFile){
    # Read parameter file if argument is present or the XML file is missing
    &set_framework_components;
    &read_text(&expand_param($ParamFile), "ReadClipboard");
    &write_xml_file($XmlFile, "DefaultEditorAndClipboard");
    `cp $ParamFile $ParamFile~` if -f $ParamFile;
}else{
    # Read XML file
    &read_xml_file($XmlFile);
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

&write_xml_file($XmlFile) if $Submit;

exit 0;

##############################################################################

sub read_xml_file{

    my $File = shift;

    my $DoDebug = ($Debug =~ /read_xml_file/);
    warn "start read_xml_file from XmlFile=$File\n" if $DoDebug;

    open(XMLFILE, $File) 
	or die "$ERROR: could not open input file $File\n";

    $nSession   = 0;
    @SessionRef = ();
    %Editor     = ();
    %Clipboard  = ();
    my $iSection;
    my $iItem;

    while($_ = <XMLFILE>){

	warn "Line $. = $_" if $DoDebug;

	if(/<MODE FRAMEWORK=\"(\d)\" (COMPONENTS=\"([^\"]*)\")?/){
	    $Framework = $1;
	    if($Framework and $3){
		$ValidComp = $3;
		%CompVersion = split( /[,\/]/ , $3 );
		print "ValidComp=$ValidComp\n" if $DoDebug;
		print "CompVersion=",join(', ',%CompVersion),"\n" if $DoDebug;
	    }
	}elsif(/<SESSION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $nSession++;
	    $iSection=0;
	    my $Name=$1;
	    my $View=$2;
	    $SessionRef[$nSession]{NAME} = $Name;
	    $SessionRef[$nSession]{VIEW} = $View;
	    print "nSession=$nSession Name=$Name View=$View\n" if $DoDebug;
	}elsif(/<SECTION NAME=\"([^\"]*)\" VIEW=\"([\w]+)\"/){
	    $iSection++;
	    $iItem=0;
	    my $Name=$1;
	    my $View = $2;
	    $SessionRef[$nSession]{SECTION}[$iSection]{VIEW} = $View;
	    $SessionRef[$nSession]{SECTION}[$iSection]{NAME} = $Name;

	    print "Section=$iSection name=$Name View=$View\n" if $DoDebug;

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
	    while($_ = <XMLFILE>){
		last if /<\/ITEM>/;
		$ItemRef->{BODY} .= $_;
	    }

	    warn "iItem=$iItem Type=$Type View=$View Body=$ItemRef->{BODY}" 
		if $DoDebug;

	}elsif(/<EDITOR 
	       \ SELECT=\"([^\"]*)\"
	       \ INSERT=\"([^\"]*)\"
	       \ FILE  =\"([^\"]*)\"
	       \ ABC   =\"([^\"]*)\"/x){
	    $Editor{SELECT} = $1;
	    $Editor{INSERT} = $2;
	    $Editor{FILE}   = $3;
	    $Editor{ABC}    = $4;
	}elsif(/<CLIPBOARD SECTION=\"([^\"]*)\" TYPE=\"([^\"]*)\"/){
	    $Clipboard{SECTION} = $1;
	    $Clipboard{TYPE}    = $2;
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

    my $DoDebug = ($Debug =~ /modify_xml_data/);

    my %Form;
    %Form = split (/[:;]/, $Submit);
    warn "Submit=$Submit\n" if $DoDebug;
    warn join(', ',%Form),"\n" if $DoDebug;

    if($_ = $Form{submit}){
	if( /^CHECK\b/ ){
	    open(FILE, ">$BakFile") 
		or die "$ERROR: could not open $BakFile\n";
	    my $Text = &write_text;
	    my $nErrorOld = ($Text =~ s/^(ERROR|WARNING).*\n+//mg);
	    print FILE $Text;
	    close FILE;

	    my $TestScript = ($Framework ? "Scripts" : ".")."/TestParam.pl";
	    $CheckResult = `$TestScript $BakFile 2>&1`;
	    if(not $CheckResult){
		$CheckResult = "No errors found";
		&read_text($Text) if $nErrorOld;
		return;
	    }
	    my @Text; @Text = split("\n", $Text);
	    my $iLine;
	    my $ErrorWarning;
	    my $Section;
	    for (split(/\n/, $CheckResult)){
		if(/TestParam_ERROR: parameter errors for (\w+)/){
		    $Section = " for $1" if $Framework;
		}elsif(/(Error|Warning) at line (\d+)/){
		    $ErrorWarning = uc($1);
		    $iLine = $2;

		    if(/for session (\d+):/){
			$iLine-- while $Text[$iLine] !~ /^Begin session/;

			# Move error to #BEGIN_COMP XX if possible
			if($Section =~ /for ([A-Z][A-Z])$/){
			    my $Comp = $1;
			    my $jLine = $iLine;

			    $jLine++ while $Text[$jLine] !~ 
				/^(End session|\#BEGIN_COMP $Comp)/;

			    $iLine = $jLine if $1 ne "End session";
			}
		    }else{
			# Move error above the last command
			$iLine-- while $iLine > 0 and $Text[$iLine] !~ /^\#/;
			$iLine--;
		    }
		}elsif(/IsFirstSession is false/){
		    $Text[$iLine] .= " it is not in the first session";
		}elsif(/^\t+(.*)/){
		    $Text[$iLine] .= ($ErrorWarning ? 
				      "\n$ErrorWarning$Section: $1" : 
				      "<BR>         $1");
		    $ErrorWarning = "";
		}
	    }
	    $Text = join("\n",@Text)."\n";
	    &read_text($Text);
	}elsif( /^SAVE$/ ){
	    &write_text_file($ParamFile);
	}elsif( /^SAVE AS$/ ){
	    $Form{FILENAME} =~ s/\.(xml|txt)$//;
	    my $File = $Form{FILENAME}; $File =~ s/^\s+//; $File =~ s/\s+$//;
	    if(open(FILE,">$File")){
		close(FILE);
		$ParamFile = $File;
		rename($XmlFile, "$ParamFile.xml");
		$XmlFile   = "$ParamFile.xml";
		$BakFile  = "$ParamFile.bak";
		&write_text_file($ParamFile);
	    }else{
		$Editor{READFILENAME}="SAVE AS";
		$Editor{NEWFILENAME}="Could not open $File" if $File;
	    }
	}elsif( /^OPEN$/ ){
	    $Form{FILENAME} =~ s/\.(xml|txt)$//;
	    my $File = $Form{FILENAME}; $File =~ s/^\s+//; $File =~ s/\s+$//;
	    if($File){
		# Make a safety save
		&write_text_file($BakFile);
		unlink($XmlFile);
		unlink($BakFile) unless `diff $BakFile $ParamFile 2>&1`;

		# Use new filenames
		$ParamFile = $File;
		$XmlFile   = "$ParamFile.xml";
		$BakFile  = "$ParamFile.bak";
		`cp $ParamFile $ParamFile~` if -f $ParamFile;

		# Read new file
		&read_text(&expand_param($ParamFile), "ReadClipboard");
	    }else{
		$Editor{READFILENAME}="OPEN";
		$Editor{NEWFILENAME}="";
	    }
	}elsif( /^REOPEN$/ ){
	    # Make a safety save and reread file
	    &write_text_file($BakFile);
	    unlink($BakFile) unless `diff $BakFile $ParamFile 2>&1`;
	    &read_text(&expand_param($ParamFile), "ReadClipboard");
	}elsif( /^READ FILE FOR INSERT$/ ){
	    my $File = $Form{FILENAME}; $File =~ s/^\s+//; $File =~ s/\s+$//;
	    if(-f $File and open(MYFILE, $File)){
		$Clipboard{BODY}    = join('',<MYFILE>);
		$Clipboard{TYPE}    = "FILE";
		$Clipboard{SECTION} = "";
		close MYFILE;
		$Editor{FILE} = $File;
	    }else{
		$Editor{READFILENAME}= "READ FILE FOR INSERT";
		$Editor{NEWFILENAME} = "could not open $File" if $File;
	    }
	}elsif( /^CANCEL$/ ){
	    # Do nothing
	    $Editor{FILE} = "" if $Editor{FILE} eq "ENTER FILENAME";
	}elsif( /^SAVE AND EXIT$/ ){
	    # save the file then kill the job
	    &write_text_file($ParamFile);
	    unlink($XmlFile, $BakFile, $IndexPhpFile, $EditorHtmlFile, 
		   $ParamHtmlFile, $ManualHtmlFile, $JumpHtmlFile);
	    kill(-9, getpgrp);
	}elsif( /^EXIT$/ ){
	    # make a safety save then kill the job
	    &write_text_file($BakFile);
	    unlink($XmlFile, $IndexPhpFile, $EditorHtmlFile, 
		   $ParamHtmlFile, $ManualHtmlFile, $JumpHtmlFile);
	    unlink($BakFile) unless `diff $BakFile $ParamFile 2>&1`;
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
	    $JumpHere = $iSession;
	}
    }elsif($_ = $Form{action}){

	my $id= $Form{id};

	if( /select_session/ ){
	    # Change insert item unless it is the same level as before
	    $Editor{INSERT} = "NEW" unless 
		($Editor{SELECT} =~ s/(\d+)/$1/g) == ($id =~ s/(\d+)/$1/g);
	    $Editor{SELECT} = $id;
	    $JumpHere       = $id;

	    if($id =~ /^all(.*)/){
		&set_view($1);
	    }elsif($id =~ /^(\d+)$/){
		&set_view("session$1");
	    }elsif($id =~ /^(\d+,\d+)/){
		&set_view("section$1");
	    }
	    return;
	}elsif( /select_insert/ ){
	    $Editor{INSERT} = $id;
	    return;
	}elsif( /select_file/ ){
	    $Editor{FILE} = $id;
	    if(open(MYFILE, $id)){
		$Clipboard{BODY}    = join('',<MYFILE>);
		$Clipboard{TYPE}    = "FILE";
		$Clipboard{SECTION} = "";
		close MYFILE;
	    }else{
		$Clipboard{TYPE} = '';
		$Clipboard{BODY} = '';
		$Clipboard{SECTION} = '';
		$Editor{READFILENAME} = "READ FILE FOR INSERT" if
		    $id = "ENTER FILENAME";
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
	$Editor{SELECT} = $id if ($id =~ /^[\d,]+$/ and /edit_|insert_/);
	$JumpHere       = $id if $id =~ /^[\d,]+$/;

	my $Reread = 0; 

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
	    # Set the insert item to the edited command so the manual shows up
	    ($Editor{INSERT}) = split("\n", $ItemRef->{BODY}) if 
		$ItemRef->{TYPE} eq "COMMAND";
	}elsif( /remove_session/ or /copy_session/ ){
	    $Clipboard{SECTION} = "";
	    $Clipboard{TYPE}    = "SESSION";
	    $Clipboard{BODY}    = &write_text($iSession);
	    $Editor{SELECT}     = "allsessions";
	    $Editor{INSERT}     = "PASTE SESSION";
	    if(/remove_session/){
		splice(@SessionRef,$iSession,1);
		$nSession--;
	    }
	}elsif( /remove_section/ or /copy_section/ ){
	    $Clipboard{SECTION} = $NameSection;
	    $Clipboard{TYPE}    = "SECTION";
	    $Clipboard{BODY}    = &write_text($iSession, $iSection);
	    $Editor{SELECT}     = $iSession;
	    $Editor{INSERT}     = "PASTE SECTION";
	    if( /remove_section/ ){
		splice(@{$SessionRef->{SECTION}}, $iSection, 1);
	    }
	}elsif( /remove_item/ or /copy_item/ ){
	    $Clipboard{TYPE}    = $ItemRef->{TYPE};
	    $Clipboard{SECTION} = $NameSection;
	    $Clipboard{SECTION} = "" if $Clipboard{TYPE} eq "COMMENT";
	    $Clipboard{BODY}    = $ItemRef->{BODY};
	    $Editor{SELECT}     = "$iSession,$iSection";
	    $Editor{INSERT}     = "PASTE COMMAND/COMMENT";
	    if( /remove_item/ ){
		splice(@{$SectionRef->{ITEM}}, $iItem, 1);
	    }
	}elsif( /insert_session/ ){
	    my $NewSessionRef;
	    $NewSessionRef->{VIEW} = "EDIT";
	    $NewSessionRef->{NAME} = "";
	    $NewSessionRef->{SECTION}[1]{VIEW} = "MAX";
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{TYPE} = "COMMENT";
	    $NewSessionRef->{SECTION}[1]{ITEM}[1]{VIEW} = "MAX";

	    if($Editor{INSERT} =~ /PASTE|FILE/ ){
		$NewSessionRef->{SECTION}[1]{ITEM}[1]{BODY} = $Clipboard{BODY};
		$Reread = 1;
	    }elsif($Editor{INSERT} =~ /NEW/){
		$NewSessionRef->{SECTION}[1]{ITEM}[1]{BODY} = "New session\n";
	    }
	    if($nSession > 0){
		splice(@SessionRef, $iSession, 0, $NewSessionRef);
	    }else{
		$SessionRef[1] = $NewSessionRef;
	    }
	    $nSession++;
	}elsif( /insert_section/ ){
	    my $NewSectionRef;
	    $NewSectionRef->{VIEW} = "MAX";
	    $NewSectionRef->{ITEM}[1]{VIEW} = "MAX";
	    $NewSectionRef->{ITEM}[1]{TYPE} = "COMMENT";
	    if($Editor{INSERT} =~ /PASTE|FILE/ ){
		$NewSectionRef->{NAME} = $Clipboard{SECTION} if 
		    $Clipboard{SECTION} ne "CON";
		$NewSectionRef->{ITEM}[1]{BODY} = $Clipboard{BODY};
		$Reread = 1;
	    }elsif($Editor{INSERT} =~ /Section (\w+)/){
		$NewSectionRef->{NAME} = $1 if $1 ne "CON";
		$NewSectionRef->{ITEM}[1]{BODY} = "New $1 section\n";
	    }
	    if( @{$SessionRef->{SECTION}} ){
		splice(@{$SessionRef->{SECTION}}, $iSection,0, $NewSectionRef);
	    }else{
		$SessionRef->{SECTION}[1] = $NewSectionRef;
	    }
	}elsif( /insert_item/ ){
	    my $NewItemRef;
	    $NewItemRef->{TYPE} = $Clipboard{TYPE};
	    if($Editor{INSERT} =~ /^PASTE|FILE/){
		$NewItemRef->{TYPE} = $Clipboard{TYPE};
		$NewItemRef->{TYPE} = "COMMENT" if $Clipboard{TYPE} eq "FILE";
		$NewItemRef->{VIEW} = "MAX";
		$NewItemRef->{BODY} = $Clipboard{BODY};
		$Reread = 1 if $NewItemRef->{TYPE} eq "COMMENT" and 
		    $NewItemRef->{BODY} =~ /^\#/m;
	    }elsif($Editor{INSERT} =~ /^COMMENT/){
		$NewItemRef->{TYPE} = "COMMENT";
		$NewItemRef->{VIEW} = "EDIT";
	    }elsif($Editor{INSERT} =~ /^\#USERINPUT/){
		$NewItemRef->{TYPE} = "USERINPUT";
		$NewItemRef->{VIEW} = "EDIT";
	    }elsif($Editor{INSERT} =~ /^\#/){
		$NewItemRef->{TYPE} = "COMMAND";
		$NewItemRef->{VIEW} = "EDIT";
		$NewItemRef->{BODY} = "CommandExample\n";
	    }else{
		# invalid choice like "COMMANDS BY GROUP";
		return;
	    }
	    if(@{$SectionRef->{ITEM}}){
		splice(@{$SectionRef->{ITEM}}, $iItem, 0, $NewItemRef);
	    }else{
		$SectionRef->{ITEM}[1] = $NewItemRef;
	    }
	}elsif( /Save command|Save userinput|Save comment/ ){
	    $ItemRef->{VIEW}="MAX";
	    $Form{text} =~ s/^\s*//;        # remove leading space
	    $Form{text} =~ s/[ \t]+\n/\n/g; # clean up line endings

	    if(/Save command/){
		# Keep first line (name of command), and set the rest
		$ItemRef->{BODY} =~ /\n/;
		$ItemRef->{BODY} = $`."\n".$Form{text};

	    }elsif(/Save comment/){
		$ItemRef->{BODY} = $Form{text};
		$ItemRef->{BODY} = "no comment\n" 
		    if $ItemRef->{BODY} =~ /^\s*$/;
		$Reread = 1 if $ItemRef->{BODY} =~ /^\#/m;
	    }else{
		$ItemRef->{BODY} = $Form{text};
	    }
	    $ItemRef->{BODY} =~ s/\n*$/\n/;      # exactly one \n at end
	}elsif( /CANCEL/ ){
	    $ItemRef->{VIEW}="MAX";
	    $ItemRef->{BODY} = "no comment\n" if 
		$ItemRef->{TYPE} eq "COMMENT" and $ItemRef->{BODY} =~ /^\s*$/;
	}

	# Reread text if we inserted something unusual like a file or a comment
	# containing \# at the beginning of the line
	if($Reread){
	    &read_text(&write_text);
	    $Editor{SELECT} = "all";
	    $JumpHere = "";
	}
    }
}

##############################################################################

sub write_xml_file{

    my $File       = shift;
    my $UseDefault = shift;

    open(XMLFILE, ">$File") 
	or die "$ERROR: could not open XML file $File for output!\n";

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
		print XMLFILE $Item->{BODY};
		print XMLFILE "\t\t\t\t\t</ITEM>\n";
	    }
	    print XMLFILE "\t\t\t\t</SECTION>\n";
	}
	print XMLFILE "\t\t\t</SESSION>\n";
    }

    if($UseDefault){
	$Editor{SELECT} = "all";
	$Editor{INSERT} = "NEW SESSION";
	$Editor{FILE}    = "";
	$Editor{ABC}    = "0";
	$Clipboard{SECTION} = "";
	$Clipboard{TYPE}    = "FILE";
    }

    print XMLFILE "
<EDITOR SELECT=\"$Editor{SELECT}\" INSERT=\"$Editor{INSERT}\" FILE=\"$Editor{FILE}\" ABC=\"$Editor{ABC}\"/>
<CLIPBOARD SECTION=\"$Clipboard{SECTION}\" TYPE=\"$Clipboard{TYPE}\">
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

    open(FILE, $ParamXml) or return ": no $ParamXml\n";

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
<?php Exec('share/Scripts/ParamConvert.pl') ?>
<#@ \$form=join(';',%FORM);                             #>
<#@ `share/Scripts/ParamConvert.pl -submit='\$form'`; #>
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

    if($Editor{READFILENAME}){
	$EditButtons = 
"      <TD COLSPAN=2 ALIGN=CENTER>
        <FONT $TopFileNameFont>Enter file name:</FONT>
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
<INPUT TYPE=SUBMIT NAME=submit VALUE=REOPEN>
<INPUT TYPE=SUBMIT NAME=submit VALUE=OPEN>
      </TD>
      <TD ALIGN=CENTER>
<FONT $TopFileNameFont>$ParamFile</FONT>
      </TD>
      <TD ALIGN=RIGHT>
<INPUT TYPE=SUBMIT NAME=submit VALUE=\"SAVE AND EXIT\">&nbsp;&nbsp;&nbsp;
<INPUT TYPE=SUBMIT NAME=submit VALUE=EXIT>&nbsp;&nbsp;&nbsp;
<A HREF=share/Scripts/ParamEditorHelp.html TARGET=HELP>Help</A>
      </TD>
";
    }
    my $SessionSection="
    <OPTION VALUE=all        >ALL
    <OPTION VALUE=allitems   >ALL ITEMS
    <OPTION VALUE=allsections>ALL SECTIONS
    <OPTION VALUE=allsessions>ALL SESSIONS
";

    my $iSession;
    for $iSession (1..$nSession){
	print "iSession=$iSession\n" if $DoDebug;

	my $SessionRef = $SessionRef[$iSession];
	my $SessionName = "Session $iSession";
	$SessionName .= ": $SessionRef->{NAME}" if $SessionRef->{NAME};
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

    if($Selected =~ /^all(\w*)$/ ){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SESSION\n"
	    if $Clipboard{TYPE} =~ /SESSION|FILE/;
        $InsertList .= "    <OPTION>NEW SESSION\n";

    }elsif($Selected =~ /^\d+$/ and $Framework){

	$InsertList  = "    <OPTION>FILE\n";
	$InsertList .= "    <OPTION>PASTE SECTION\n"
	    if $Clipboard{TYPE} =~ /SECTION|FILE/;
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

    if( not ($InsertList =~ s/OPTION>$Insert$/OPTION SELECTED>$Insert/m)){
	$InsertList =~ 
	    s/OPTION>(COMMANDS|NEW SESSION|Section CON)/OPTION SELECTED>$1/;
	$Insert = $1;
	$Editor{INSERT} = $1;
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
    <OPTION>ENTER FILENAME
$Files
  </SELECT>
";
	# Add SELECTED if file has been loaded already
	my $File = $Editor{FILE};
	$InsertItem =~ s/<OPTION(.*$File\n)/<OPTION SELECTED$1/ 
	    if $File and ($File eq "ENTER FILENAME" or 
			  $Clipboard{TYPE} eq "FILE" and $Clipboard{BODY});
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
    }elsif($CommandExample){
	$Manual =  $CommandText;
	$Manual =~ s/\n\n/\n<p>\n/g;
	$Manual =  "<H1>Manual</H1>\n<PRE>\n$CommandExample\n</PRE>\n$Manual";
	$Manual =~ s/\n$//;
    }elsif( $Editor{INSERT} =~ /^PASTE/ ){
	$Manual = "<H1>Clipboard: $Clipboard{SECTION} $Clipboard{TYPE}</H1>\n".
	    "<PRE>\n$Clipboard{BODY}\n</PRE>";
    }elsif( $Editor{INSERT} eq "FILE" and $Clipboard{TYPE} eq "FILE" and
	    $Editor{FILE} ){
	$Manual = "<H1>FILE: $Editor{FILE}</H1>\n".
	    "<PRE>\n$Clipboard{BODY}\n</PRE>";
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
    $JumpHere = $Editor{SELECT} unless $JumpHere;

    # Jump to previous item if item index is larger than 2
    $JumpHere =~ s/(\d+,\d+,)([2-9]|\d\d+)$/"$1".($2-1)/e;

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
"                 if $Editor{SELECT} =~ /^all/;

	$CopySessionButton = "      <$Action=copy_session
><IMG SRC=$ImageDir/button_copy.gif TITLE=\"Copy session\"></A>
";
	$RemoveSessionButton = "      <$Action=remove_session
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove session\"></A>
";

	# Place anchor to selected session
	$Param .= $Here if $JumpHere eq $iSession;

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
	    $Param .= $Here if $JumpHere eq "$iSession,$iSection";

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
      <$Action=select_session TITLE=\"Select section\">
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

		$ItemRef->{BODY} =~ s/CommandExample/$CommandExample/;
		my $ItemHead;
		my $ItemTail;
		if($ItemType eq "userinput"){
		    $ItemHead = '#USERINPUT';
		    $ItemTail = $ItemRef->{BODY};
		}else{
		    ($ItemHead, $ItemTail) = split("\n",$ItemRef->{BODY},2);
		}
		$ItemType = "error" if $ItemHead =~ /^(ERROR|WARNING)/;

		my $TableColor = $TableColor{$ItemType};

		$Action =~ s/\d+\&action$/$iItem\&action/;
		$InsertItemButton =~ s/\d+\&action=/$iItem\&action=/;

		my $nLine = ($ItemTail =~ s/\n/\n/g);

		# Place anchor to selected item
		$Param .= $Here if $JumpHere eq "$iSession,$iSection,$iItem";

		if($ItemView eq "EDIT"){

		    if($ItemType eq "comment"){
			$ItemHead = "Comment:";
			$ItemTail = $ItemRef->{BODY};
		    }elsif($ItemType eq "userinput"){
			$ItemTail .= ("\n" x 10);
		    }
		    $nLine = ($ItemTail =~ s/\n/\n/g) + 2;

		    $Param .= "
  <FORM NAME=item_editor ACTION=$IndexPhpFile>
  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$ItemEditorBgColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width>
      </TD>
      <TD ALIGN=LEFT><FONT COLOR=BLUE>
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

		my $ActionEditItem = "$Action=edit_item ".
		    "TITLE=\"Edit $ItemType in session $iSession/$SectionName\"";

		$ActionEditItem = "A NAME=ERROR" if $ItemType eq "error";

		# Create buttons
		if($ItemTail and $ItemType ne "error"){
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

		$CopyItemButton = "" if $ItemType eq "error";

		$RemoveItemButton = "      <$Action=remove_item
><IMG SRC=$ImageDir/button_remove.gif TITLE=\"Remove $ItemType\"></A>
";

		# Show first line with usual buttons

		$Param .=
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width ALIGN=RIGHT>
$MinMaxItemButton
      </TD>
      <TD ALIGN=LEFT><$ActionEditItem>
$ItemHead
      </A></TD>
      <TD ALIGN=RIGHT ALIGN=TOP>
$InsertItemButton$CopyItemButton$RemoveItemButton
      </TD>
    </TR>
  </TABLE>
";

		if( ($ItemView eq "MIN" and $ItemType ne "error") 
		    or not $ItemTail){
                    $Param .= "<p>\n" if $ItemType !~ /comment|error/;
		    next;
                }

		if($ItemType ne "command"){

		    $ItemTail =~ s/\n/<BR>\n/g;

		    $Param .= 
"  <TABLE BORDER=0 WIDTH=$LeftTableWidth BGCOLOR=$TableColor>
    <TR>
      <TD WIDTH=$LeftColumn1Width>
      </TD>
      <TD COLSPAN=2><$ActionEditItem>
$ItemTail
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
      <TD ALIGN=LEFT><$ActionEditItem>
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

		$Param .= "<p>\n" if $ItemType !~ /comment|error/;

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
      <$Action=select_session TITLE=\"Select section\">
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
	
	$iSection = ($nSection+1 or 1);
	$InsertSectionButton = 
"    <A HREF=$IndexPhpFile?id=$iSession,$iSection\&action=insert_section
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert section\"></A>
" 	    if $Editor{SELECT} =~ /^\d+$/;

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

    $Param .= "$SessionLine\n";

    $Param .=
"  <TABLE WIDTH=$LeftTableWidth BGCOLOR=$LeftBgColor><TR><TD ALIGN=RIGHT>
    <A HREF=$IndexPhpFile?id=".($nSession+1)."\&action=insert_session
><IMG SRC=$ImageDir/button_insert.gif TITLE=\"Insert session\"></A>
  </TD></TR></TABLE>
"                 if $Editor{SELECT} =~ /^all/;

    $Param .= "</BODY>\n";

    print PARAM $Param;
    close PARAM;
}

##############################################################################
sub set_view{
    $_ = @_[0];

    my $SessionView = ( /session|\d/ ? "MIN" : "MAX");
    my $SectionView = ( /section|\d/ ? "MIN" : "MAX");
    my $ItemView    = ( /items/      ? "MIN" : "MAX");

    my $iSession;
    my $iSection;
    my $iItem;

    for $iSession (1..$nSession){
	if($SessionRef[$iSession]{VIEW} ne "EDIT"){
	    $SessionRef[$iSession]{VIEW} = $SessionView;

	    $SessionRef[$iSession]{VIEW} = "MAX" if /ion$iSession/;
	}
	next if /session/;

	for $iSection (1..$#{ $SessionRef[$iSession]{SECTION} }){
	    my $SectionRef = $SessionRef[$iSession]{SECTION}[$iSection];

	    $SectionRef->{VIEW} = $SectionView;

	    $SectionRef->{VIEW} = "MAX" if /ion$iSession,$iSection/;

	    next if /section/;
	    
	    for $iItem (1...$#{ $SectionRef->{ITEM} }){
		my $ItemRef = $SectionRef->{ITEM}[$iItem];

		$ItemRef->{VIEW} = $ItemView;
	    }
	}
    }

}

##############################################################################
sub convert_type{

    my $InputFile  = $ARGV[0];
    my $OutputFile = $ARGV[1];
    die "$ERROR: input and output filenames are the same: $InputFile\n"
	if $InputFile eq $OutputFile;

    my $InputType  = ( ($InputFile  =~ /\.(expand|xml)$/) ? $1 : "txt");
    my $OutputType = ( ($OutputFile =~ /\.xml$/) ? "xml" : "expand");

    die "$ERROR: input and output file types are the same: $InputType\n"
	if $InputType eq $OutputType;

    die "$ERROR: input file $InputFile does not exist\n"
	unless -f $InputFile;

    &set_framework_components;

    # Read input file
    if($InputType eq "xml"){
	&read_xml_file($InputFile);
    }else{
	# Expand file (this is needed for files with #INCLUDE)
	&read_text(&expand_param($InputFile), "ReadClipboard");
    }

    # Write output file
    if($OutputType eq "xml"){
	&write_xml_file($OutputFile, "UseDefaultEditorAndClipboard");
    }else{
	&write_text_file($OutputFile);
    }
    exit 0;
}

##############################################################################
sub set_framework_components{

    # Set framework mode to true if CON/Control directory is present
    $Framework = -d "CON/Control";

    return unless $Framework;

    # Set list of components
    my $MakefileDef = "Makefile.def";
    open(INFILE,$MakefileDef) or
	print "$ERROR could not open file $MakefileDef\n";

    $ValidComp = '';
    while(<INFILE>){
	next unless /^(\w\w)_VERSION\s*=\s*(\w+)/;
	$ValidComp .= "$1/$2," unless $2 eq "Empty";
	last if /TIMING_VERSION/;
    }
    $ValidComp =~ s/,$//;

    close(INFILE);
}
##############################################################################
sub read_text{

    my $ReadClipboard;
    my @Text;
    if($_[0] eq "ReadClipboard"){
	$ReadClipboard = 1;
    }else{
	@Text = split(/\n/, shift);
	$ReadClipboard = shift;
    }

    # Delete previous values if any
    @SessionRef = ();

    # Initialize for checking input parameter file
    $nSession=1;            # number of sessions
    my $nLine=0;            # Line number in the text
    my $iSection=1;         # index of section in session
    my $iItem=0;            # index of item in section
    my $Section="";         # name of component section
    my $IsCommand=0;        # true while reading a command
    my $IsComment=0;        # true while reading a comment
    my $UserInput=0;        # true between #USERINPUTBEGIN and #USERINPUTEND

    my $SectionRef;         # Pointer to current section in $SessionRef

    # We assume that there is at least one session with one section inside it
    $SessionRef[1]{VIEW} = "MAX";               # default view of session 1
    $SessionRef[1]{NAME} = "";                  # default name of session 1
    $SessionRef[1]{SECTION}[1]{VIEW} = "MAX";   # default view of section 1
    $SessionRef[1]{SECTION}[1]{NAME} = "";      # default name of section 1
    $SectionRef = $SessionRef[1]{SECTION}[1];

    if(not @Text){
	# Add a comment if there is no text passed
	$SectionRef->{ITEM}[1]{VIEW} = "MAX";
	$SectionRef->{ITEM}[1]{TYPE} = "COMMENT";
	$SectionRef->{ITEM}[1]{BODY} = "New parameter file\n";

	return;
    }

    while(@Text){

	$_ = shift(@Text);

	if(/^\#END\b/){

	    # Remove the last section if it has no items
	    pop( @{ $SessionRef[$nSession]{SECTION} } )
		if $#{ $SectionRef->{ITEM} } < 0;
	    
	    # Remove the last session if it has no sections
	    if( $#{ $SessionRef[$nSession]{SECTION} } < 1){
		pop( @SessionRef );
		$nSession--;
	    }
	    
	    return unless $ReadClipboard;

	    # Read the rest of the text into the 'clipboard'
	    $Clipboard{BODY} = join("\n", @Text) . "\n" if @Text;
	    $Clipboard{BODY} =~ s/^\n+$/\n/;
	    return;
	}

	$nLine++;

	if($UserInput and /^\#(BEGIN_COMP|END_COMP|RUN|USERINPUTBEGIN)\b/){
	    &print_error( $nLine, " for command $_".
			  "\tthis command cannot occur after ".
			  "#USERINPUTBEGIN at line $UserInput");
	    $UserInput = 0;
	}

	# Check for BEGIN_COMP and END_COMP commands
	if(/^\#BEGIN_COMP\b/){
	    if(not $Framework){
		&print_error( $nLine, " for command $_".
			     "\tshould not be used in stand-alone mode");
		next;
	    }
	    if($Section){
		&print_error( $nLine, " for command $_".
			     "\talready had BEGIN_COMP $Section");
		next;
	    }
	    # Figure out which component is beginning here
	    ($Section) = /BEGIN_COMP ([A-Z][A-Z])/ or
		&print_error( $nLine, " for command $_".
			     "\tincorrectly formatted BEGIN_COMP command");

	    $iSection++ if $#{ $SectionRef->{ITEM} } >= 0;
	    $SessionRef[$nSession]{SECTION}[$iSection]{NAME} = $Section;
	    $SessionRef[$nSession]{SECTION}[$iSection]{VIEW} = "MAX";
	    $SectionRef = $SessionRef[$nSession]{SECTION}[$iSection];
	    $iItem = 0;
	    $IsComment = 0; $IsCommand = 0; $UserInput = 0;

	}elsif(/^\#END_COMP/){
	    if(not $Framework){
		&print_error( $nLine, " for command $_".
			     "\tshould not be used in stand-alone mode");
		next;
	    }

	    if(not $Section){
		&print_error( $nLine, 
			      " for command $_".
			      "\tmissing \#BEGIN_COMP command");
		next;
	    }

	    # Extract name of the component from #END_COMP ID
	    my $Comp;
	    ($Comp) = /END_COMP ([A-Z][A-Z])/ or
		&print_error( $nLine, 
			      " for command $_".
			      "\tincorrectly formatted END_COMP command");

	    # Check if the component name matches
	    if($Comp ne $Section){
		&print_error( $nLine, 
			      " for command $_\tcomponent does not match".
			      " BEGIN_COMP $Section");
	    }

	    # In any case return to next CON section
	    $iSection++ if $#{ $SectionRef->{ITEM} } >= 0;
	    $SessionRef[$nSession]{SECTION}[$iSection]{NAME} = "";
	    $SessionRef[$nSession]{SECTION}[$iSection]{VIEW} = "MAX";
	    $SectionRef = $SessionRef[$nSession]{SECTION}[$iSection];
	    $Section = '';
	    $iItem   = 0;
	    $IsComment = 0; $IsCommand = 0; $UserInput = 0;

	}elsif(/^\#RUN/){	# Check if a new session has started
	    # Check if the required commands are defined and
	    # if the parameters are correct for the session

	    if($Section and $Framework){
		print "Error in session $nSession ending at line $nLine ".
		    "in expanded file $ParamFile:\n".
		    "\tBEGIN_COMP $Section does not have matching END_COMP\n";
		$Section = '';
	    }

	    # Remove the last section if it has no items
	    pop( @{ $SessionRef[$nSession]{SECTION} } )
		if $#{ $SectionRef->{ITEM} } < 0;

	    # Note: sections are indexed from 1
	    $nSession++ if $#{ $SessionRef[$nSession]{SECTION} } > 0;

	    $SessionRef[$nSession]{VIEW} = "MAX";
	    $SessionRef[$nSession]{SECTION}[1]{VIEW} = "MAX";
	    $SectionRef = $SessionRef[$nSession]{SECTION}[1];
	    $iSection=1;
	    $iItem=0;
	    $IsComment = 0; $IsCommand = 0; $UserInput = 0;

	}elsif( /^(Begin|End)\s+session:\s*(.*)/i ){
	    # Session names are stored as comments
	    # Begin session: or End session:
	    my $Name = $2;
	    $SessionRef[$nSession]{NAME}= $Name unless $Name =~ /^\d+$/;
	    $IsComment = 0; $IsCommand = 0; $UserInput = 0;

	}elsif(/^\#USERINPUTBEGIN/){
	    $iItem++;
	    $SectionRef->{ITEM}[$iItem]{VIEW}="MAX";
	    $SectionRef->{ITEM}[$iItem]{TYPE}="USERINPUT";
	    $UserInput = $nLine+1;
	    $IsCommand = 0;
	    $IsComment = 0;
	}elsif(/\#USERINPUTEND/){
	    if(not $UserInput){
		&print_error( $nLine, 
			      " for command $_\tthere is no matching".
			      " USERINPUTBEGIN command");
	    }
	    $UserInput = 0;
	}else{

	    # Analyze line for various items
	    if(/^\#/ and not $UserInput){
		$iItem++;
		$SectionRef->{ITEM}[$iItem]{VIEW}="MAX";
		$SectionRef->{ITEM}[$iItem]{TYPE}="COMMAND";
		$IsCommand=1;
		$IsComment=0;
	    }elsif( /^\s*$/ ){
		# Empty line closes commands and comments
		$IsCommand = 0;
		$IsComment = 0;
		next unless $UserInput;
	    }elsif( /^(ERROR|WARNING)/ or 
		    (not $IsCommand and not $UserInput and not $IsComment)){
		# non-empty line starts new comment
		$iItem++;
		$IsComment = 1;
		$SectionRef->{ITEM}[$iItem]{VIEW}="MAX";
		$SectionRef->{ITEM}[$iItem]{TYPE}="COMMENT";
	    }

	    if($iItem == 0){
		print "Error iItem=0 for line: $_";
	    }

	    # Store line
	    $SectionRef->{ITEM}[$iItem]{BODY} .= "$_\n";
	}
    }
}
##############################################################################
sub print_error{
    my $nLine = shift;
    print "$ERROR at line $nLine in expanded file",@_,"\n";
}
##############################################################################
sub write_text{

    my $iSession0 = shift;
    my $iSection0 = shift;

    my $Output;

    my @Sessions;
    if($iSession0){
	@Sessions = ($iSession0);
    }else{ 
	@Sessions = (1..$#SessionRef);
    }
    my $iSession;
    for $iSession (@Sessions){
	my $NameSession = ($SessionRef[$iSession]{NAME} or $iSession);
	$Output .= "Begin session: $NameSession\n\n" unless $iSession0;

	my @Sections;
	if($iSection0){
	    @Sections = ($iSection0);
	}else{
	    @Sections = (1..$#{ $SessionRef[$iSession]{SECTION} });
	}
	my $iSection;
	for $iSection (@Sections){
	    my $SectionRef = $SessionRef[$iSession]{SECTION}[$iSection];

	    # Skip empty section ($iItem starts with 1)
	    next unless $#{ $SectionRef->{ITEM} } > 0;

	    my $NameSection = $SectionRef->{NAME};

	    $Output .= "\#BEGIN_COMP $NameSection ".("-" x 40)."\n\n"
		if $NameSection and not $iSection0;

	    my $iItem;
	    for $iItem (1..$#{ $SectionRef->{ITEM} }){

		my $ItemRef = $SectionRef->{ITEM}[$iItem];
		my $ItemType= $ItemRef->{TYPE};

		if($ItemType eq "USERINPUT"){
		    $Output .= "#USERINPUTBEGIN ----------------\n".
			$ItemRef->{BODY}.
			"#USERINPUTEND   ----------------\n\n";
		}else{
		    $Output .= $ItemRef->{BODY};
		    $Output .= "\n" if $ItemType eq "COMMAND";
		}
	    }
	    $Output .= "\#END_COMP $NameSection   ".("-" x 40)."\n\n"
		if $NameSection and not $iSection0;
	}
	$Output .= "End session: $NameSession\n\#END ".("\#" x 60)."\n" 
	    unless $iSession0;
    }

    # Replace #END with #RUN when a session follows
    $Output =~ s/\#END (\#+\nBegin session:)/\#RUN $1/g;

    return $Output;
}
##############################################################################
sub write_text_file{
    my $FileName = shift;
    my $NoClipboard = shift;

    open(FILE, ">$FileName") 
	or die "$ERROR: write_text_file could not open $FileName\n";
    print FILE &write_text; 
    print FILE $Clipboard{BODY} unless $NoClipboard or not $Clipboard{BODY};
    close FILE;
}
##############################################################################
sub expand_param{

    my $basefile = @_[0];

    return unless -f $basefile;

    # Check if basefile is in local directory or not
    if($basefile =~ /\/([^\/]+)$/){

	# Change to directory so that the include files can be read
	my $dir = $`;
	$basefile = $1;
	my $pwd = `pwd`; chop $pwd;
	chdir $dir or die "$ERROR: could not cd $dir\n";

	# Expand the file recursively
	my $result = &process_file($basefile, 'fh00');

	# Go back to original directory
	chdir $pwd;

	return $result;
    }else{
	# Expand the file recursively
	return &process_file($basefile, 'fh00');
    }
}
##############################################################################

sub process_file {
    no strict;
    local($filename, $input) = @_;
    local($output);

    $input++;
    open($input, $filename) or die"$ERROR: cannot open $filename: $!\n";
    while (<$input>) {
	# Stop reading if #END command is read
        last if /^#END\b/;

	# Check for #INCLUDE
        if (/^#INCLUDE\b/) {
	    # Read file name following #INCLUDE
            $includefile=<$input>;
	    # process include file recursively
	    $output .= &process_file($includefile,$input);
        }else{
	    # Print line as it is otherwise
	    $output .= $_;
	}
    }
    if($input eq "fh01"){
	$output .= '#END ' . ('#' x 60) ."\n". join('',<$input>);
    }
    close $input;

    return $output;
}

##############################################################################
#BOP
#!ROUTINE: ParamConvert.pl - convert between parameter file formats
#!DESCRIPTION:
# This script is usually called internally from the parameter editor.
# 
#!REVISION HISTORY:
# 10/19/2007 G.Toth - initial version integrated from ParamXmlToHtml.pl, 
#                     ParamTextToXml.pl, ParamXmlToText.pl and ExpandParam.pl
#                      
#EOP
sub print_help{

    print
#BOC
"Purpose:

     Convert between various formats of the input parameter file.
     The most important use of this script is to convert the input
     parameter file into several HTML and PHP files that serve as
     the parameter editor graphical user interface (GUI).

     The GUI can be customized by creating a ParamEditor.conf file in 
     the user's home directory. Type 

grep '^our' share/Scripts/ParamConvert.pl

     to see the variables that can be modified using the same syntax, 
     but different values. Note the use of semicolons at the end of lines.

     Depending on the various applications, this script is reading
     several files, including Makefile.def to get the list of components
     for the SWMF, the PARAM.XML files of the SWMF and the components to
     get the list and description of commands. The script also executes
     the TestParam.pl script to check the correctness of the parameter file.
     For this reason it is normally executed from the main directory of the
     SWMF or the stand-alone physics model.

Usage:

  ParamConvert.pl -h

    -h            print help message and stop.

  ParamConvert.pl [-submit=FORM] [PARAMFILE]

    -submit=FORM  FORM is a semi-colon separated list of form variables and 
                  their values. This parameter is normally passed by the 
                  parameter editor GUI.

    PARAMFILE     Name of the plain text input parameter file to be edited. 
                  By default the input file name is obtained from the HTML 
                  <TITLE> of the param.html file. If there is no param.html
                  file present, the default PARAMFILE is 'run/PARAM.in'.

  ParamConvert.pl INPUTFILE OUTPUTFILE

    INPUTFILE     Name of the input parameter file. The extensions  
                  .expand (expanded file with no \#INCLUDE files) and 
                  .xml    (XML enhanced input parameter file)
                  are recognized. Everything else is taken as a 
                  plain parameter file with possible \#INCLUDE files.

    OUTPUTFILE    Name of the output parameter file. The extension
                  .xml (XML enhanced input parameter file)
                  is recognized, everything else is taken as an
                  expanded parameter file with no \#INCLUDE files.

Examples:

    Create GUI files (index.php, editor.html, param.html...) from run/UAM.in:

share/Scripts/ParamConvert.pl run/UAM.in

    Execute some action of the GUI for the file in the TITLE of param.html:

share/Scripts/ParamConvert.pl -submit='submit;SAVE AND EXIT'

    Convert the plain text file with included files into a single file:

share/Scripts/ParamConvert.pl run/PARAM.in run/PARAM.expand

    Convert the expanded text file into an XML enhanced file:

share/Scripts/ParamConvert.pl run/PARAM.expand run/PARAM.xml

    Convert the XML enhanced text file into a plain (expanded) text file:

share/Scripts/ParamConvert.pl run/PARAM.in.xml run/PARAM.in"
#EOC
    ."\n\n";
    exit 0;
}

