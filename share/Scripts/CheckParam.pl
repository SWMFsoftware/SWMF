#!/usr/bin/perl -s

# Read a parameter file and verify its correctness based on an XML description.
# Type CheckParam.pl -h for help.

# This subroutine is outside the scope of strict and the "my" variables
# This is a safety feature, so the global variables can not be changed
sub eval_comp{eval("package COMP; $_[0]")}

# Read command line options
my $Debug       = $D; undef $D;
my $Help        = $h; undef $h;
my $Interactive = $i; undef $i;
my $Verbose     = $v; undef $v; 
my $ReadTree    = $r; undef $r;
my $SaveTree    = $s; undef $s;
my $XmlFile     = $x; undef $x;
my $NameComp    = $c; undef $c;
my $Components  = $C; undef $C;
my $Precision   = $p; undef $p;
my $GridSize    = $g; undef $g;
my $nProc       = $n; undef $n;
my $StandAlone  = $S; undef $S;

use strict;

# Pattern to match component ID-s
my $ValidComp = 'SC|IH|GM|IE|IM|UA';

# Error string
my $ERROR = 'CheckParam_ERROR:';

# Set default values (needed in the help message too)
my $XmlFileDefault   = 'Param/PARAM.XML';
my $InputFileDefault = 'run/PARAM.in';

# Print help message and exit if -h switch was used
&print_help if $Help;

# Reset filenames to defaults if needed
$XmlFile  = $XmlFileDefault unless $XmlFile;

# Default tree file has the same name as the XML file with .pl extension
if($SaveTree eq '1'){
    $SaveTree = $XmlFile;
    $SaveTree =~ s/xml$/pl/i or 
	die "$ERROR could not replace .XML extension with .pl for XML file\n".
	    "$XmlFile\n";
}


# Check the correctness of the command line arguments
&check_arguments;

if($SaveTree){
    # Convert XML description into Perl tree data structure and exit
    &parse_xml($XmlFile) or 
	die "$ERROR XML::Parser package is not installed\n";
    exit 1;
}

if(not $ReadTree){
    # Default tree file has the same name as the XML file with .pl extension
    $ReadTree = $XmlFile;
    $ReadTree =~ s/xml$/pl/i or $ReadTree = '';
}
    
# Initialize for checking input parameter file
my $InputFile;
$InputFile       = $ARGV[0] or 
    $InputFile   = $InputFileDefault;
my $nLine        = 0; # Line number in the current file
my $IncludeLevel = 0; # Each include file increases the include level by 1
my @nLine;            # Array storing line numbers for all include levels
my @InputFile;        # Array storing file names for all include levels
my $FileHandle;       # File handle for the current file

# these global variables are needed for checking parts of strings without
# complaints and for showing possible options if all string parts fail
my $DontShowParamError = 0; 
my $optionvalues;

# Read the parameter definitions
my $commandList;
$commandList = &parse_xml($XmlFile) or
    $commandList = &read_tree($ReadTree) or
    die "$ERROR neither the XML::Parser::EasyTree package is installed\n".
    "\tnor the tree file $ReadTree was found!\n";

# Store values in the COMP package to be used and changed by the XML rules
&init_comp;

# Set the default values defined at the top level of the XML file
foreach my $node (@{$commandList->{content}}){
    &set_value($node) if $node->{type} eq 'e' and $node->{name} eq 'set';
}

# Find commands
my $commandName;        # Current command name (can be an alias)
my $realName;           # The real name of the command (not the alias)
my %realName;           # Hash for real names  $realName{$alias}=$realName
my %command;            # Hash for tree nodes  $command{$realName}=$node
my %commandText;        # Hash for (help) text $commandText{$realName}=$text
my @required;           # Array of required commands @required=($realName,...)

&find_commands($commandList);

# Read and check the parameter file from STDIN or $ARGV
my %defined;            # stores the line number for each command read
my %definedSessionLast; # stores the commands read in last session
my $nSession=1;         # current session number
my $paramName;          # name of the current parameter
my $paramType;          # type of the current parameter
my $paramValue;         # value of the current parameter
my $InsideComp;         # the component for which parameters are read
my $UserInput;          # set to true between #USERINPUTBEGIN and #USERINPUTEND

# Set the file handle for the top level file (or STDIN for interactive mode)
if($Interactive){
    $FileHandle= \*STDIN;
    $InputFile = 'STDIN';
}else{
    $FileHandle="F$IncludeLevel";
    # Change into local directory if necessary
    $InputFile =~ s/(.*)\///;
    if($1){
	chdir $1 or
	    die "$ERROR Could not change directory to $1\n";
	print "chdir $1\n" if $Debug;
    }

    no strict;
    open($FileHandle,$InputFile) or
	die "$ERROR Could not open parameter input file $InputFile!\n";
}

while($_=&read_line){

    # Read command of form #COMMANDNAME
    if(/^\#(\w+)/){
	$commandName=$1;
    }else{
	$commandName='';
	next;
    }

    # Check and store the command
    next unless &check_command($commandName);

    # Get the real name for aliases (e.g. #BODY instead of #MAGNETOSPHERE)
    $realName = $realName{$commandName};

    print "Read command name=$commandName\n" if $Debug;

    # Read the parameters for the command
    &process_elements($command{$realName});

}
# Check if the final session has the required commands defined and
# if the parameters are correct
&check_session;

exit 0;

##############################################################################
sub check_arguments{

    # Check command line arguments

    if($XmlFile and not -f $XmlFile){
	die "$ERROR Could not find XML file $XmlFile\n" 
	    if $XmlFile and not -f $XmlFile;
    }

    if($ReadTree and $SaveTree){
	die "$ERROR do not read (-r) and save (-s) tree files "
	    ."at the same time!\n"
	    ."For help type CheckParam.pl -h\n";
    }

    if($XmlFile eq $SaveTree){
	die "$ERROR do not overwrite XML file $XmlFile "
	    ."by saving tree into the same file!\n"
	    ."Choose a different TREEFILE in -s=TREEFILE";
    }

    if($NameComp and $NameComp !~ /^$ValidComp$/){
	die "$ERROR -c=$NameComp is not among valid component ID-s"
	    ." $ValidComp\n";
    }

    if($Components){
	foreach (split ',',$Components){
	    if(not /^$ValidComp$/){
		die "$ERROR -C=$_ is not among valid component ID-s"
		    ." $ValidComp\n";
	    }
	}
    }

    if($Precision and not $Precision =~ /^single|double$/){
	die "$ERROR -p=$Precision is not 'single' or 'double'\n";
    }

    if($GridSize and not $GridSize =~ /^\d+(,\d+)*$/){
	die "$ERROR -g=$GridSize is not"
	    ." a comma separated list of integers\n";
    }

    if(length($nProc) and not $nProc =~ /^[1-9]\d*$/){
	die "$ERROR -n=$nProc is not a positive integer\n";
    }

    if($Interactive and $ARGV[0]){
	die "$ERROR -i should not be combined"
	    ." with the file argument $ARGV[0]\n";
    }
}

##############################################################################
sub init_comp{

    # Set variables in the COMP package for use in the XML files

    # Set COMP::_IsStandAlone according to the -S switch
    $COMP::_IsStandAlone = $StandAlone;

    # Set COMP::_IsFirstSession to true
    $COMP::_IsFirstSession = 1;

    # Set the size of the grid 
    if($GridSize){
	@COMP::_GridSize = split(',',$GridSize);
    }

    # Set the number of processors
    if($nProc){
	$COMP::_nProc = $nProc;
    }

    # Set COMP::_nByteReal for easy check on precision
    $COMP::_nByteReal = $Precision eq 'single' ? 4 : 8;

    if($Components){
	# Set COMP::_Components for the XML rules to check components
	$COMP::_Components = $Components;

	# Create a COMP::_Registered hash for registered components
	foreach (split ',',$Components){
	    $COMP::_Registered{$_} = 1;
	}

	# Initialize COMP::_UsedComp hash with the registered components
	%COMP::_UsedComp = %COMP::_Registered;
    }
}
##############################################################################
sub parse_xml{
    # Parse the XML file and return a pointer to the parsed tree
    my $XmlFile=$_[0];

    # Try loading the modules
    eval('
    use XML::Parser;
    use XML::Parser::EasyTree;
    $XML::Parser::EasyTree::Noempty=1;
    open(FILE,$XmlFile) or die "$ERROR could not open XML file $XmlFile!\n";
    close(FILE);
    ');
    if($@){
	# If it did not succeed return false value and try something else
	print "XML description could not be parsed.\n" if $Debug;
	return 0 if $@;
    }

    my $parser = new XML::Parser(Style => 'EasyTree');
    my $tree   = $parser->parsefile($XmlFile);
    my $commandList = $tree->[0];

    if($commandList->{type} ne 'e' or $commandList->{name} ne 'commandList'){
	die "$ERROR first node should be an element named commandList\n";
    }

    if($SaveTree){
	# Save $tree into file $SaveTree

	# Try Saving the tree with Data::Dumper module
	require Data::Dumper or
	    die "$ERROR Data::Dumper package is not installed\n";
	$Data::Dumper::Indent=0;
	$Data::Dumper::Variable="tree";

	open TREE, ">$SaveTree" or 
	    die "$ERROR could not open tree file $SaveTree for writing\n";
	print TREE "#^"."CFG FILE _FALSE_\n";
	print TREE Data::Dumper->Dump([$tree],["tree"]);
	close TREE or
	    die "$ERROR could not close saved tree file $SaveTree\n";
	print "XML tree read from $XmlFile has been saved into $SaveTree\n";
	exit 0;
    }

    return $commandList;
}

##############################################################################
sub read_tree{
    # read a tree file named $_[0] containing the parsed XML tree
    my $ReadTree=$_[0];

    my $tree=do($ReadTree);
    return 0 unless length($tree)>0;
    print "Description tree was read from file $ReadTree\n" if $Debug;
    my $commandList = $tree->[0];
    if($commandList->{type} ne 'e' or $commandList->{name} ne 'commandList'){
	die "$ERROR first node should be an element named commandList\n";
    }
    return $commandList;
}

##############################################################################
sub set_value{
    # set values defined by <set name=".." type=".." value="..."> XML tag
    my $node  = $_[0];
    my $name  = &extract($node->{attrib}->{name}, "string");
    my $type  = &extract($node->{attrib}->{type}, "string");
    my $value = &extract($node->{attrib}->{value}, $type);

    &store($name,$value);
}

##############################################################################
sub find_commands{
    # a recursive traverse of the XML tree finds all the commands
    # and stores them in the global %command hash and  @required array
    my $node = $_[0];
    if($node->{type} eq 'e'){
	if($node->{name} eq 'command'){
	    $realName = $node->{attrib}->{name};
	    $realName{$realName}=$realName;
	    $command{$realName} =$node;

	    print "found command $realName\n" if $Debug;
	    
	    if($Verbose){
		# Store the text of the command
		foreach my $element (@{$node->{content}}){
		    $commandText{$realName}=$element->{content} if
			$element->{type} eq 't';
		}
	    }

	    # Store the aliases into %realName
	    foreach $commandName (split(',',$node->{attrib}->{alias})){
		$realName{$commandName}=$realName;
	    }

	    if($node->{attrib}->{required}){
		my $required=&extract($node->{attrib}->{required},"logical");
		push(@required,$realName) if $required;
	    }

	}
	foreach my $element (@{$node->{content}}){
	    &find_commands($element);
	}
    }
}

##############################################################################
sub read_line{

    # Read and return next line from the parameter file(s)

    # Loop until we run out of EOF (end of file), #END and #INCLUDE commands

    # Pre 5.6 versions of Perl do not allow <$FileHandle> with strict
    no strict;
    while((not $_=<$FileHandle>) or /^\#(END|INCLUDE)\b/){

	if($_){
	    # We found an #END or an #INCLUDE command
	    $nLine++;
	    print "At line $nLine in file $InputFile script command $_"
		if $Debug;
	}else{
	    print "EOF at line $nLine in file $InputFile\n" if $Debug;
	}

	if(/INCLUDE/){
	    # Read the parameter of the INCLUDE command
	    $commandName=$1;
	    $realName=$realName{$commandName};
	    no strict;
	    if($paramValue = <$FileHandle>){
		use strict;
		$nLine++; 
		chop $paramValue;
		$paramValue =~ s/^ +//;  # Remove leading spaces
		# Remove everything after 3 spaces or a TAB
		$paramValue =~ s/   .*//;
		$paramValue =~ s/\t.*//;

		print "read script parameter = $paramValue\n" if $Debug;
	    }else{
		use strict;
		print "Error at line $nLine in file $InputFile:\n".
		    "\tend of file after script command \#$commandName\n";
		return 0 unless &previous_file;
	    }
	    use strict;
	    my $file;
	    # The parameter of the INCLUDE command is the file name itself
	    $file = $paramValue;
	    if(not $file){
		$paramName = "NameIncludeFile";
		$paramType = "string";
		&param_error("is missing!");
		return 1;
	    }
	    # Try to open include file with new file handle
	    my $FileHandleOld = $FileHandle;
	    $FileHandle       = "F".($IncludeLevel+1);
	    no strict;
	    if(-f $file and open($FileHandle, $file)){
		use strict;

		print "Opened include file '$file'\n" if $Debug;

		# Save current file name and line number
		$InputFile[$IncludeLevel] = $InputFile;
		$nLine[$IncludeLevel]     = $nLine;
		
		# Move to next include level
		$IncludeLevel++;
		$InputFile                = $file;
		$nLine                    = 0;
	    }else{
		use strict;
		print "Error at line $nLine in file $InputFile for ".
		    "command $_".
			"\tCould not open include file $file\n";
		$FileHandle=$FileHandleOld;
		return 1;
	    }
	}else{
	    # File either ended or an #END command was found
	    # Close file and return to previous file or return 0
	    return 0 unless &previous_file;
	}
    }
    use strict;

    if($UserInput and /^\#(BEGIN_COMP|END_COMP|RUN|USERINPUTBEGIN)\b/){
	print "Error at line $nLine in file $InputFile for ".
	    "command $_".
	    "\tthis command cannot occur after #USERINPUTBEGIN at line $UserInput\n";
	    $UserInput = 0;
    }

    # Check for BEGIN_COMP and END_COMP commands
    if(/^\#BEGIN_COMP\b/ and not $StandAlone){
	if($InsideComp){
	    print "Error at line $nLine in file $InputFile for ".
		"command $_".
		"\tAlready had BEGIN_COMP $InsideComp\n";
	}else{
	    # Figure out which component is beginning here
	    ($InsideComp) = /BEGIN_COMP ([A-Z][A-Z])/ or
		print "Error at line $nLine in file $InputFile for ".
		"command $_".
		"\tIncorrectly formatted BEGIN_COMP command\n";

	    # Check if the component is registered and ON
	    if($Components){
		if(not $COMP::_Registered{$InsideComp}){
		    print "Error at line $nLine in file $InputFile for ".
			"command $_".
			"\tComponent $InsideComp is not registered\n";
		}
		if(not $COMP::_UsedComp{$InsideComp}){
		    print "Error at line $nLine in file $InputFile for ".
			"command $_".
			"\tRegistered component $InsideComp is OFF.\n";
		}
	    }

	    # Return an empty line
	    $_="\n";
	}
    }elsif(/^\#END_COMP/ and not $StandAlone){
	# Extract name of the component from #END_COMP ID
	my $Comp;
	($Comp) = /END_COMP ([A-Z][A-Z])/ or
	    print "Error at line $nLine in file $InputFile for ".
	    "command $_".
	    "\tIncorrectly formatted END_COMP command\n";

	# Check if the component name matches
	if($Comp ne $InsideComp){
	    print "Error at line $nLine in file $InputFile for ".
		"command $_".
		"\tComponent does not match BEGIN_COMP $InsideComp\n";
	}else{
	    # END_COMP matches BEGIN_COMP, so we are back to CON params
	    $InsideComp = '';
	    # Return an empty line
	    $_="\n";
	}
    }elsif(/^\#RUN/){ # Check if a new session has started
	# Check if the required commands are defined and
	# if the parameters are correct for the session

	print "Session $nSession is complete\n" if $Debug;

	&check_session;
	undef %definedSessionLast;
	$nSession++;
	$COMP::_IsFirstSession=0;
    }elsif(/^\#USERINPUTBEGIN/){
	$UserInput = $nLine+1;
    }elsif(/^\#USERINPUTEND/){
	if(not $UserInput){
	    print "Error at line $nLine in file $InputFile for ".
		"command $_".
		"\tthere is no matching #USERINPUTBEGIN command\n";
	}
	$UserInput = 0;
    }
    
    $nLine++;

    # Return the line only for the selected component outside user input
    if($InsideComp eq $NameComp and not $UserInput){
	return $_;
    }else{
	return "\n";
    }
}
##############################################################################
sub previous_file{
    # Select previous file if it exists

    if($IncludeLevel > 0){
	no strict;
	close $FileHandle;
	use strict;

	# Restore previous include level
	$IncludeLevel--;
	$InputFile = $InputFile[$IncludeLevel];
	if($InputFile eq 'STDIN'){
	    $FileHandle = \*STDIN;
	}else{
	    $FileHandle="F$IncludeLevel";
	}
	$nLine    = $nLine[$IncludeLevel];
	return 1;
    }else{
	return 0;
    }
}

##############################################################################
sub check_if{
    my $node = $_[0];
    my $if;

    if($node->{name} =~ /^(if|rule)$/ )
    {
	# The condition is in the compulsary attribute "expr" 
	# for tags <if...> and <rule ...>

	if(not $if=$node->{attrib}->{expr})
	{
	    print "ERROR: attribute expr= is missing from tag $node->{name}".
		"in the XML description of $commandName\n";
	    return 1;
	}
    }
    else
    {
	# The condition is in the optional attribute "if" for other tags
	$if=$node->{attrib}->{if};
    }

    # Evaluate the conditional expression
    if($if){
	my $check=&extract($if,"logical");

	# Rules need to be processed if condition is false
	if($node->{name} eq 'rule'){$check = 1-$check};

	print "if=$if check=$check\n" if $Debug;

	return $check;
    }else{
	return 1;
    }
}

##############################################################################
sub process_elements{
    my $node = $_[0];
    foreach my $element (@{$node->{content}}){
	next if $element->{type} eq 't';
	next unless &check_if($element);
	my $name = $element->{"name"};
	if($name eq 'parameter')
	{
	    if(&read_parameter($node,$element)){
		&check_value_format and
		    &check_param_value($element,$paramValue) and
			&store($paramName,$paramValue);
	    }else{
		&param_error("could not be read from file");
	    }
	}
	elsif($name eq 'set')
	{
	    &set_value($element);
	}
	elsif($name eq 'defined')
	{
	    &check_command($element->{attrib}->{name}); # check and store
	}
	elsif($name eq 'if')
	{
	    &process_elements($element);
	}
	elsif($name eq 'rule')
	{
	    print "Error at line $nLine in file $InputFile:\n".
		"\t$element->{attrib}->{expr}\n",
		"\tis false! Rule description:\n",
		"$element->{content}->[0]->{content}\n";
	    print "Command description\n$commandText{$realName}\n"
		if $Verbose;
	}
	elsif($name eq 'for')
	{
	    my $index=&extract($element->{attrib}->{name},"string");
	    my $from =&extract($element->{attrib}->{from},"integer");
	    my $to   =&extract($element->{attrib}->{to},  "integer");
	    print "FOR index $index= $from to $to\n" if $Debug;
	    foreach my $value ($from..$to){
		# Set the index name in package
                &store($index,$value) if $index;
		&process_elements($element);
	    }
	}
	elsif($name eq 'foreach')
	{
	    my $index =&extract($element->{attrib}->{name},"string");
	    my $values=&extract($element->{attrib}->{values},"string");
	
	    print "FOREACH index=$index values=$values\n" if $Debug;
	    foreach my $value (split(',',$values)){
		# Set the index name in package
		&store($index,$value);
		&process_elements($element);
	    }
	}
    }
}

##############################################################################
sub read_parameter{
    my $command=$_[0];
    my $param=$_[1];

    my $attrib=$param->{attrib};

    $paramName=&extract($attrib->{name},"string");
    $paramType=lc($attrib->{type});

    # read the value of the parameter
    print "$paramName=" if $Debug;
    $paramValue = &read_line or return 0;
    chop $paramValue;

    # Remove leading spaces
    $paramValue =~ s/^ +//;

    # Remove everything after 3 spaces or a TAB
    $paramValue =~ s/   .*//;
    $paramValue =~ s/\t.*//;

    print "reading $paramName = $paramValue\n" if $Debug;

    return 1;
}

##############################################################################
sub check_value_format{
    my $bad=0;

    if($paramType =~ /logical/)
    {
	# T F .true. .false.
	$bad = &bad_value unless 
	    $paramValue =~ s/^\s*(f|\.false\.)$/0/i or 
	    $paramValue =~ s/^\s*(t|\.true\.)$/1/i;
    }
    elsif($paramType =~ /integer/)
    {
	# 3 +30 -124
	$bad = &bad_value unless $paramValue =~ /^\s*[\+\-]?\d+$/;
    }
    elsif($paramType =~ /real/)
    {

	# replace exponent D with E
	$paramValue =~ s/d/e/i;

	# 3 -3 -3. .21 -3.21 -3.21e2 -3.21D+21
	$bad = &bad_value
	    unless $paramValue =~ 
		/^\s*[\+\-]?(\d+(\.\d*)?|\.\d+)(e[\+\-]?\d+)?/i;
    }
    elsif($paramType =~ /string/)
    {
	# strings cannot have wrong format
    }
    else
    {
	print "For command \#$commandName parameter $paramName ".
	    "has unknown type $paramType !?";
	$bad=1;
    }

    return(not $bad);
}

##############################################################################
sub check_param_value{
    # Check value
    my $node =$_[0];
    my $value=$_[1];

    my $content = $node->{content};
    my $attrib  = $node->{attrib};
    my $type    = &extract($attrib->{type},"string");

    # Set the values for the parts
    my $good=1;
    if($type eq 'strings')
    {
	# check if the number of parts is correct
	my @parts  =split(' ',$paramValue);
	my $nPart  =scalar @parts;
	my $min    =&extract($attrib->{min});
	if($nPart < $min){
	    &param_error("has $nPart parts ",
			 "which is less than minimum $min");
	    return 0;
	}
	my $max=&extract($attrib->{max});
	if($max and $max < $nPart){
	    &param_error("has $nPart parts ",
			 "which is greater than maximum $max");
	    return 0;
	}
	my $duplicate=&extract($attrib->{duplicate},"logical");
	if($duplicate and $nPart < $max){	    
	    foreach my $iPart ($nPart..$max-1){
		$parts[$iPart]=$parts[$nPart-1];
	    }
	    $nPart=$max;
	    print "nPart=$nPart parts=",@parts,"\n";
	}
	my $ordered=&extract($attrib->{ordered},"logical");
	my $iPart=0;
	foreach my $element (@$content)
	{
	    next unless $element->{type} eq 'e';
	    if($element->{name} eq 'part'){
		my $name     = &extract($element->{attrib}->{name},"string");
		my $required = &extract($element->{attrib}->{required},
					"logical");
		my $multiple = &extract($element->{attrib}->{multiple},
					"multiple");
		if($ordered){
		    if(not &check_param_value($element,$parts[$iPart])){
			print "\tin part $iPart named $name\n";
			$good=0;
		    }else{
			&store($name,$value);
		    }
		}else{
		    # Unordered name parts. Check one by one.
		    my $any=0;
		    $DontShowParamError=1;
		    my $store;
		    foreach my $jPart (0..$nPart-1){
			my $part=$parts[$jPart];
			if(length($part)>0 and 
			   &check_param_value($element,$part)){
			    $store.="$part ";
			    $parts[$jPart] = "";
			    $any = 1;
			    last unless $multiple;
			}
		    }
		    chop $store;
		    &store($name,$store) if $any;
		    $DontShowParamError=0;
		    if($required and not $any){
			&param_error("has value '$paramValue'",
				     " with no part matching\n",
				     "\t$optionvalues for $name");
			$good=0;
		    }
		}
		$iPart++;
	    }	    
	}
	if(join('',@parts) and not $ordered){
	    &param_error("has value '$paramValue' containing >> ",
			 join(' ',@parts),
			 " <<\n\twhich does not match anything");
	}
	return $good;
    }

    my $good=0;
    if($attrib->{input} eq 'select')
    {
	# Read options
	$optionvalues="";
	my $optionvalue;
	my $shown=0;

      OPTION: foreach my $element (@$content){
	    next OPTION unless $element->{type} eq 'e';
	    next unless &check_if($element);
	    if($element->{name} eq 'option')
	    {
		if(length($element->{attrib}->{value})>0)
		{
		    $optionvalue = 
			&extract($element->{attrib}->{value},$type);
		}
		else
		{
		    $optionvalue = 
			&extract($element->{attrib}->{name},$type);
		}
		# Collect the possible option values
		$optionvalues .= "$optionvalue,";

		if( length($value)>0 and 
		    (($type eq 'string' and $optionvalue =~/\b$value\b/) or
		     ($type eq 'integer' and $optionvalue == $value) or
		     ($type eq 'real'    and abs($optionvalue-$value)<1e-7)
		     ))
		{
		    $good=1;
		    last OPTION;
		}
	    }
	    elsif($element->{name} eq 'optioninput')
	    {
		$good=&check_range($element,$value,$type);
		if(not $good){
		    chop $optionvalues;
		    print "\tand is not found among options $optionvalues\n";
		    $shown = 1;
		}
	    }
	}
	if(not $good and not $shown){
	    chop $optionvalues;
	    if(not length($value)){
		 &param_error("is missing! Possible options are ",
			      $optionvalues);
	     }else{
		 &param_error("has value $value not listed among options ",
			      $optionvalues);
	     }
	}
    }
    else
    {
	$good=&check_range($node,$value,$type);
    }
    return $good;
}

##############################################################################
sub check_range{
    my $node =$_[0];
    my $value=$_[1];
    my $type =$_[2];

    my $min=&extract($node->{attrib}->{min},$type);
    if(length($min)>0){
	if($value < $min){
	    &param_error("has value $value which is less than minimum $min");
	    return 0;
	}
    }
    my $max=&extract($node->{attrib}->{max},$type);
    if(length($max)>0){
	if($value > $max){
	    &param_error("has value $value which is greater ",
			 "than maximum $max");
	    return 0;
	}
    }
    my $length=&extract($node->{attrib}->{length},"integer");
    if(length($length)>0){
	if(length($value) > $length){
	    &param_error("has length ",length($value),
		" which is longer than maximum length $length");
	    return 0;
	}
    }
    return 1;
}

##############################################################################
sub bad_value{

    &param_error("has an incorrectly formatted value '$paramValue'");

    return 1;
}

##############################################################################
sub check_command{
    my $commandName=$_[0];

    my $realName = $realName{$commandName};
    my $command  = $command{$realName};

    # Check if the command is descibed in the XML file or not
    if( not $command ){
	print "Error at line $nLine in file $InputFile:\n".
	    "\tcommand \#$commandName is unknown\n";
	return 0;
    }

    # Check if the command is currently available
    if( not &check_if($command) ){
	print "Error at line $nLine in file $InputFile:\n".
	    "\tcommand \#$commandName is not available\n",
	    "\tbecause condition $command->{attrib}->{if} is false\n";
	return 0;
    }


    # Check if command was defined before in this session
    if($definedSessionLast{$realName} and not $command->{attrib}->{multiple}){
	print "Warning at line $nLine in file $InputFile:\n".
	    "\tcommand \#$commandName has already been defined ",
	    "in this session:\n",
	    "\t$definedSessionLast{$realName}\n";
    }

    # Store the command
    $defined{$realName}="\#$commandName at line $nLine in file $InputFile";
    $definedSessionLast{$realName}=$defined{$realName};

    print "Command \#$commandName defined at line $nLine in file $InputFile\n"
	if $Debug;

    return 1;
}

##############################################################################
sub check_session{

    print "InsideComp: $InsideComp\n" if $Debug;

    if($InsideComp){
	print "Error in session $nSession ending at line $nLine ".
                "in file $InputFile:\n".
                "\tBEGIN_COMP $InsideComp does not have matching END_COMP\n";
    }


    print "required commands: ",join(',',@required),"\n" if $Debug;

    foreach $realName (@required){
	if(not $defined{$realName}){
	    print "Error in session $nSession ending at line $nLine ".
		"in file $InputFile:\n".
		"\trequired command \#$realName was not defined\n";
	}
    }

    # Check the rules
    foreach my $node (@{$commandList->{content}}){
	next unless $node->{type} eq 'e';
	if($node->{name} eq 'rule'){

	    print "checking global rule $node->{attrib}->{expr}\n"
		if $Debug;

	    next unless check_if($node);

	    print "Error at line $nLine in file $InputFile in ".
		"session $nSession:\n",
		"\t$node->{attrib}->{expr}\n",
		"\tis false! Rule description:\n",
		"$node->{content}->[0]->{content}\n";
	}
    }
}
##############################################################################
sub store{
    # Usage:
    #      store($name,$value);
    # Store value $value into variable $name in package COMP

    my $name=$_[0];
    my $value=$_[1];
    &eval_comp("\$$name='$value'");
    print "eval_comp(\$$name='$value')\n" if $Debug;
    print "ERROR at line $nLine in file $InputFile:\n".
	"\tfor command $commandName\n",
	"\tcould not evaluate \$$name='$value':\n\t$@\n"
	    if $@;
}
##############################################################################
sub extract{
    # Usage:
    #        $result = &extract($expression, $type);
    #
    # $result will get the value extracted from $expression.
    # The variables containing a $ sign are substituted from package
    # COMP and the expressions are evaluated
    # The type given by the second argument can be 
    # "logical", "integer", "real" or "string"
    #
    # Example:
    #       $result = &extract('$nI+2','integer')
    #
    # the function will return the value of
    #
    #     $COMP::nI + 2
    #
    # into $result.

    my $value  = $_[0]; return if length($value) == 0;
    my $type   = $_[1];

    my $value_ = $value;

    if($value_ =~ /\$/){
	if($type eq "string"){
	    $value='"'.$value.'"';
	}else{
	    $value='('.$value.')';
	}
	&eval_comp("\$_value_=$value;");
	$value_ = $COMP::_value_;

        print "ERROR at line $nLine in file $InputFile:\n".
		"\tcommand \#$commandName could not evaluate\n",
		"\tsetting \$value_=$value:\n\t$@\n"
		    if $@;
    }

    if($type eq "logical"){
	# replace T or .true. with 1, F or .false. with 0
	$value_ =~ s/^(f|.false.)$/0/i;
	$value_ =~ s/^(t|.true.)$/1/i;
	# Convert possible empty string (=false) into a more readable 0
	$value_ = $value_+0;
	if($value_ !~ /^[01]$/){
	    print "ERROR at line $nLine in file $InputFile:\n".
		"\tcommand \#$commandName contains incorrect logical\n",
		"\t$value = $value_\n";
	    $value_=0;
	}
    }

    return $value_;

}
##############################################################################
sub param_error{
    return if $DontShowParamError;
    print "Error at line $nLine in file $InputFile:\n".
	"\tParameter $paramName of type $paramType in ",
	"command \#$commandName\n",
	"\t@_\n";
    print "Command description:\n$commandText{$realName}" if $Verbose;
}
##############################################################################
sub print_help{

    print "
Purpose:

     Read a parameter file and verify its correctness based on an 
     XML description. In most applications this scipt is called by 
     TestParam.pl.

Usage:

  CheckParam.pl [-h] [-v] [-D] [-x=XMLFILE] [-r=TREEFILE] [-s[=TREEFILE]] 
                [-S] [-c=ID] [-C=IDLIST] [-p=PRECISION] [-i] [PARAMFILE]

  -h            print help message and stop

  -D            print debug information

  -v            verbosity adds command description to every error message.

  -x=XMLFILE    Use XMLFILE for the XML description. 
                If the -x switch is omitted the $XmlFileDefault file is used.

  -r=TREEFILE   Read TREEFILE containing the parameter description in a Perl 
                tree data structure if no XML description is found or the
                XML::Parser::Easytree module is not installed.
                The default tree file has the same name as the XML file 
                but with a .pl extension.
              
  -s[=TREEFILE] Save the parsed XML description into TREEFILE in form of a Perl
                tree data structure and exit. 
                If the -s switch is given without a TREEFILE value, 
                the tree is saved into the default tree file, which has
                the same name as the XML file but with a .pl extension.

  -S            Check parameters assuming a stand alone mode. In stand alone 
                mode the #BEGIN_COMP and #END_COMP commands are ignored.
                Also sets \$_IsStandAlone to true to be used in the XML rules.

  -c=ID         Check the parameters for the component defined by ID,
                which is the two-character name of the component.
                The default is to check the CON parameters.

  -C=IDLIST     List of registered component ID-s separated with commas.
                The default is not checking the registration.

  -g=GRIDSIZE   List of grid dimensions separated with commas.
                The default is not checking the grid size.

  -p=PRECISION  The default precision for real numbers. Possible values are
                'single' and 'double'. Default value is 'double'.

  -i            Interactive mode. The parameters are read from STDIN.

  PARAMFILE     The file containing the parameters to test. 
                The default file name is run/PARAM.in.
                In interactive mode (-i) or when the XML tree is saved (-s)
                the PARAMFILE is ignored.

Examples:

    Check CON parameters in run/PARAM.in for correctness with verbose info:

CheckParam.pl -C='GM,IH,IE' -v

    Check GM parameters in run/PARAM_new.in for correctness:

CheckParam.pl -x=GM/BATSRUS/PARAM.XML -c=GM run/PARAM_new.in

    Convert GM/BATSRUS/PARAM.XML into the tree file GM/BATSRUS/PARAM.pl:

CheckParam.pl -x=GM/BATSRUS/PARAM.XML -s

    Set the XML file name to an empty string and read a Perl tree file:

CheckParam.pl -x= -r=mytree.pl

    Check lines typed through standard input with debug info:

CheckParam.pl -D -i
#DESCRIPTION
Just a test.
...
Ctrl-D
";
    exit 0;
}
