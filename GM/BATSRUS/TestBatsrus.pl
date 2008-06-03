#!/usr/bin/perl
#^CFG COPYRIGHT UM
#^CFG FILE TESTING

use strict;

my $Self        ="TestBatsrus.pl";
my $TESTSUITE   ="Param/TESTSUITE";
my $TESTSUITErun="GM/$TESTSUITE";
my $testtmp     ="test.tmp";
my $exe_default ="BATSRUS.exe"; 
my $run_default ="mpirun -np 2";
my $rundir      ="run";

# Figure out directory from which the call is made
my $Dir = $0; $Dir =~ s/[^\/]*$//; $Dir = './' unless $Dir;

# Correct path to the test suite
$TESTSUITE=$Dir.$TESTSUITE;

# Find all the swithes in the test suite directory
opendir(DIR, $TESTSUITE)
    or die "ERROR: Could not open directory $TESTSUITE\n";
my @switch = grep {-f "$TESTSUITE/$_/DEFAULT"} readdir(DIR);

closedir(DIR);
@switch = sort @switch;

my %choice;     # hash of arrays holding possible values for switches

foreach my $switch (@switch){
    opendir(DIR, "$TESTSUITE/$switch")
	or die "ERROR: cannot open directory $TESTSUITE/$switch\n";	
    my @dir=readdir(DIR);
    close(DIR);
    my @option = grep {not /^(\.|\.\.|CFG|CVS|DEFAULT|.*\~)$/} @dir;
    my $defaultfile="$TESTSUITE/$switch/DEFAULT";
    open(DEFAULT,$defaultfile) or die "ERROR: cannot open $defaultfile\n";
    my $default = <DEFAULT>; $default=~s/\n//;
    close(DEFAULT);
    if($option[0] ne $default){
	# Make the default value the first option
	# Find the index of the default value
	my $i;
	for($i=0; $i<=$#option; $i++){
	    last if $option[$i] eq $default;
	}
	if($i>$#option){
	    print "ERROR: $defaultfile contains default value $default\n";
	    print "       which is not among choices=",
	    join(',',@option),"\n";
	    die   "       Correct $defaultfile or add choice!\n";
	}
	# Remove default from @option and put it into the first element
	splice @option, $i, 1;
	unshift @option, $default;
    }
    $choice{$switch}=[@option];
}

# Set default values
my %select;  # hash of arrays holding selected values for switches
my $switch;
foreach $switch (@switch){
    $select{$switch}=[ $choice{$switch}[0] ];
}

# Find the unique values which can be used as a shortcut
my %switch;  # hash containing the switch(es) corresponding to a value
my @choice;  # array of choices for a switch
my $value;   # one value
foreach $switch (@switch){
    foreach $value (@{$choice{$switch}}){
	if($switch{$value}){
	    $switch{$value}.=",$switch";
	}else{
	    $switch{$value}=$switch;
	}
    }
}

# Read arguments and modify the %select array

my $Verbose;          # verbosity
my $run=$run_default; # run command
my $exe=$exe_default; # executable to run
my $flg;              # flags for the executable
my $Sleep;            # number of seconds to sleep after each run
my $norun;            # true if no runs are done
my @value;            # array of values for a switch

foreach my $arg (@ARGV){

    if($arg !~ /^-/){
	print "ERROR: incorrect argument $arg\n";
	&print_usage;
    }

    # Process normal switches first
    &print_help if $arg =~ /^-(h|help)$/;

    if($arg=~/^-(l|list)$/){
	&print_hash(%choice); exit;
    }

    if($arg eq "-v"           ){$Verbose= 1; next}
    if($arg =~ /^-v=(\d)$/    ){$Verbose=$1; next}

    if($arg =~ /^-norun$/     ){$norun= 1  ; next}
    if($arg =~ /^-norun=(\d)$/){$norun=$1  ; next} 

    if($arg =~ /^-r=(.*)/     ){$run=$1    ; next}

    if($arg =~ /^-x=(.*)/     ){$exe=$1    ; next}

    if($arg =~ /^-f=(.*)/     ){$flg=$1    ; next}

    if($arg =~ /^-s=(.*)/     ){$Sleep=$1  ; next}

    if($arg =~ /^-d=(.*)/     ){$rundir=$1 ; next}

    if( ($switch,$value) = ($arg=~/-([A-Z]\w*)=(.*)$/) ){
	# -Switch=value
	if(not $choice{$switch} ){
	    print "ERROR: unknown switch -$switch ! Valid switches:\n   -";
	    print join(" -",sort keys %choice),"\n";
	    &print_usage;
	}
	@value  = split(',',$value);
	@choice = @{$choice{$switch}};
	foreach $value (@value){
	    if(not grep {$value eq $_} @choice){
		print "ERROR: value $value is unknown for switch -$switch ! ";
		print "Valid values:\n    ",join(",",@choice),"\n";
		&print_usage;
	    }
	}
	# Set selection
	$select{$switch}=[@value];
    }else{
	# -value
	($value) = ($arg=~/^-(.*)$/);
	if(not $switch{$value}){
	    print "ERROR: unknown value -$value !\n";
	    &print_usage;
	}
	if($switch{$value}=~/,/){
	    print "ERROR: -$value is not a unique value.\n",
  	          "Possible switches are $switch{$value}\n",
	          " Use -Switch=$value format!\n";
	    &print_usage;
	}
	# set selection
	$select{$switch{$value}} = [$value];
    }
}

# &print_hash(%select) if $Verbose;

# Check if we are in the main directory
if(-f "Makefile"){
    # Make run directory if necessary
    if(not -d $rundir){
	print "Making $rundir directory...\n";
	&execute("make rundir RUNDIR=$rundir STANDALONE=YES");
	die "ERROR: $rundir could not be made with make rundir\n" 
	    unless -d $rundir;
    }
    if(not $norun){
	# Check if executable exists and compile if necessary
	if($exe ne $exe_default and $exe ne "SWMF.exe"){
	    die "ERROR: No executable $rundir/$exe was found!\n" 
		unless -x "$rundir/$exe";
	}elsif(not -f not -f "src/$exe" and not "bin/$exe" and not -f "../../bin/$exe"){
	    print "Compiling $exe...\n";
	    &execute("make");
	    die "ERROR: $exe could not be compiled with make\n"
		    unless -x "bin/$exe" or -x "src/$exe" 
		    or -x "../../bin/$exe";
	}
	# Compile postprocessing executables if necessary
	my @plottype=@{$select{"Plottype"}};
	foreach my $PLT ("IDL","SPH"){
	    if(not -f "src/Post$PLT.exe" and 
	       not -f "bin/Post$PLT.exe" and 
	       not "../../bin/Post$PLT.exe"){
		print "Compiling Post$PLT.exe...\n";
		&execute("make P$PLT");
		warn "WARNING: Post$PLT.exe could not be compiled".
			" with make P$PLT\n"
		    unless -x "src/Post$PLT.exe" or -x "bin/Post$PLT.exe";
	    }
	}
    }
    chdir $rundir;
}else{
    # There is no Makefile here. Check if executable is here if needed
    if(not $norun){
	# Check if executable is there
	die "ERROR: No executable $exe was found!\n" unless -x $exe;
    }
}

# Find the number to start numbering the test.??? directories
my $n = 0;
opendir(DIR,'./');
while (my $filename = readdir(DIR)) {
    if ($filename =~ /test\.(\d\d\d)/) {
	if ($1+1 > $n) {$n = $1+1; }
    }
}
close (DIR);

# Set $SWMF to true if SWMF is run
my $SWMF = 0; $SWMF = 1 if $exe =~ /SWMF/;

# Touch RESTART.out if the executable is not SWMF
&execute("touch RESTART.out") if not $SWMF;

my $num  = "0" x (3-length($n)) . $n; # string identifying the test run
my %value;       # hash containing one set of values
&do_loop(0);     # loop over the values for the 1st switch first
                 # and continue recursively
exit;

##########################################################################
sub do_loop{
    # recursive subroutine looping through all value combinations

    my $iswitch = $_[0];
    my $switch = $switch[$iswitch];

    print "do_loop: switch=$switch\n" if $Verbose;

    my @option= @{$select{$switch}};

    print "do_loop: values=",join(',',@option),"\n" if $Verbose;

    foreach $value (@option) {
	$value{$switch}=$value;
	if($iswitch < $#switch){
	    &do_loop($iswitch+1);
	}else{
	    &do_test;
	    $num++;
	}
    }
}
##########################################################################
sub do_test{
    
    my $switch;
    my $value;
    my $nondefault;
    foreach my $switch (sort keys %value){
	$value=$value{$switch};
	if($value ne $choice{$switch}[0]){
	    $nondefault .= "-$switch=$value ";
	}
    }
    if($nondefault){
	chop($nondefault);
	print "test.$num: $nondefault\n";
    }else{
	print "test.$num: DEFAULT\n";
    }
    &execute("rm -rf $testtmp") if -d $testtmp;
    &execute("mkdir $testtmp");

    # Save switch settings into test.*/SWITCHES
    my $filename="$testtmp/SWITCHES";
    open(OUT,">$filename") or die "ERROR: cannot open $filename\n";
    print OUT "$Self $nondefault\n\n";
    foreach $switch (@switch){
	print OUT "-$switch=$value{$switch}\n";
    }
    close(OUT);
    
    &execute("cp $TESTSUITErun/PARAM.in $testtmp/PARAM.in");
    &execute("cp $TESTSUITErun/LAYOUT.in $testtmp/LAYOUT.full") if $SWMF;
    foreach $switch (@switch){
	&execute("cp $TESTSUITErun/$switch/$value{$switch} $testtmp/$switch");
    }

    # expand the parameter file
    chdir $testtmp or die "ERROR: cannot chdir $testtmp\n";
    open OUT, ">PARAM.expand" or die "ERROR: cannot open PARAM.expand\n";
    &expand_param("PARAM.in", "fh00");
    close OUT;

    if($SWMF){
	# shrink the LAYOUT file
	open IN, "LAYOUT.full" 
	    or die "ERROR: cannot open LAYOUT.full";
	open OUT,">LAYOUT.in"  
	    or die "ERROR: cannot open LAYOUT.in for writing";
	while(<IN>){
	    print OUT unless /^IE/ 
		and not $value{"Inner"} eq "full_ionosphere";
	}
	close IN;
	close OUT;
    }
    chdir "..";

    my $result;
    if($norun){
	$result=0;
    }else{
	# Run the code
	&execute("cp $testtmp/PARAM.expand PARAM.in");
	&execute("cp $testtmp/LAYOUT.in LAYOUT.in") if $SWMF;

	$result = &execute("$run ./$exe $flg 1> log.$num 2> error.$num");

	&execute("sleep $Sleep") if $Sleep > 0;

	chdir "GM" or die "ERROR: Could not enter directory GM\n";
	opendir(DIR,"IO2") or die "ERROR: Could not open directory GM/IO2\n";
	my @dir=readdir(DIR);
	closedir(DIR);
	if(grep /\.h$/, @dir){
	    print "    Post Processing IDL files.\n";
	    execute("./pIDL >> ../log.$num");
	}
	if(grep /\.T$/, @dir){
	    print "    Post Processing TEC files.\n";
	    &execute("./pTEC p r T >> ../log.$num");
	}
	if(grep /\.S$/, @dir){
	    print "    Post Processing TEC SPH files.\n";
	    &execute("./pTEC p r S >> ../log.$num");
	}
	chdir "..";

	if($SWMF){
	    chdir "IE" or die "ERROR: Could not change to directory IE\n";
	    opendir(DIR,"ionosphere") or 
		die "ERROR: Could not open directory IE/ionosphere\n";
	    my @dir=readdir(DIR);
	    closedir(DIR);
	    if(grep /_b1\./, @dir){
		print "    Post Processing IONO files.\n";
		execute("./pION -r >> ../log.$num");
	    }
	    chdir "..";
	}

	# Move output files into the temporary test directory
	my $indir;
	my $outdir;
	foreach $indir ( 'GM/IO2'
			,'IE/ionosphere'    #^CFG IF IONOSPHERE
			){
	    next unless -d $indir;
	    $outdir = $indir; $outdir =~ s/^(GM|IE)\///; 
	    $outdir = "$testtmp/$outdir";
	    &execute("mkdir $outdir");
	    opendir(DIR,$indir);
	    while (my $filename = readdir(DIR)) {
		rename "$indir/$filename", "$outdir/$filename"
		    unless $filename =~ /^\./;
	    }
	    close (DIR);
	}
	# Move error and log files into the temporary test directory
	&execute("mv error.$num log.$num $testtmp");
    }
    # Rename the temporary test directory into test.$num
    &execute("/bin/rm -rf test.$num");
    &execute("mv $testtmp test.$num");

    if ($result == 0) {
	print "test.$num: finished. Check or compare results!\n";
		    
    } else {
	print "test.$num: FAILED FAILED FAILED !!!!\n";
	if ($result == 256 or $result == 65280) {
	    print "     This was caused by a valid STOP in the code, which\n";
	    print "     may mean that there was a problem with the test, or\n";
	    print "     may indicate that there was a negative pressure.\n";
	} else {
	    print "     This was caused by an error which occured, such as\n";
	    print "     an Arithmetic exception.\n";
	}
    }
    
}
##############################################################################

sub execute{

    my $command=join(' ',@_);

    print "    $command\n" if $Verbose;
    my $result = system $command;
    
    return $result;
}

##############################################################################

sub expand_param {
    no strict 'refs';

    my $filename=@_[0];
    my $input=@_[1];
    $input++;
    open($input, $filename) || die"Can't open $filename: $!\n";
    while (<$input>) {
        last if /^#END\b/;

        if (/^#INCLUDE\b/) {
            $filename=<$input>;
	    chop $filename;
	    # Remove trailing comments following white space
	    $filename =~ s/(\S)\s.*/$1/;
	    &expand_param($filename,$input);
        }else{
	    print OUT $_;
	}
    }
    close $input;
}

##############################################################################
sub print_usage{

    print "
Usage:

    $Self -switch1 -switch2=value -switch3=value1,value2...

For a list of switches and values type

    $Self -l

For detailed information on usage type
  
    $Self -h

";
    exit 1;
}
##############################################################################

#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{GM/BATSRUS: Scripts Used by BATSRUS and SWMF}
#!QUOTE: \subsection{Testing the GM and IE Components}
#!ROUTINE: TestBatsrus.pl - test Batsrus with a set of parameters
#!DESCRIPTION:
# Run BATSRUS.exe or SWMF.exe with various parameters in the run/ directory
# and put results into test directories run/test.000, run/test.001, etc.
# Start with no or empty run/ directory to avoid loss of data.
#
#!REVISION HISTORY:
# 07/10/03 G.Toth TestBatsrus.pl modified to TestSWMF.pl
# 03/11/04 G.Toth Renamed and modified to TestBatsrus.pl so it works
#                 both in stand alone BATSRUS and in the framework.
#EOP

sub print_help{

    print "
Purpose:

    Run $exe_default with various parameters in the run/ directory
    and put results into test directories run/test.000, run/test.001, etc.
    Start with no or empty run/ directory to avoid loss of data.

",
#BOC
"Usage:

    TestBatsrus.pl [-h] [-l] [-v] [-norun] [-d=DIR]
           [-r=RUNCOMMAND] [-x=EXECUTABLE] [-f=FLAGS]
           [-value1] [-Switch2=value2] [-Switch3=value3,value4]...

  -h             help message
  -l             list switches and possible values and stop
  -v             verbose
  -norun         produce test.*/ directories with parameter files but no run.
  -d=DIR         set the name of the run directory. DEFAULT is -d='run'
  -r=RUNCOMMAND  run command. DEFAULT is -r='$run_default'
  -x=EXECUTABLE  executable.  DEFAULT is -x=$exe_default
  -f=FLAGS       flags for the executables. DEFAULT is empty string.
  -s=SLEEP       sleep SLEEP seconds between tests. DEFAULT is 0 seconds.

  -Switch=value1,value2...
                 run the code with Switch=value1, then Switch=value2, etc.
                 all switch names start with capital letter.
                 e.g.: -Solver=rusanov,roe

  -value         run the code with Switch=value. This short notation
                 is allowed only if there is no ambiguity, ie. 
                 the value name must be unique, e.g. -rusanov
                 This short notation can only be used to set a single value
                 for the switch.

The code will run with all possible combinations of the switch values.
Undefined switches have a default value. 

The switches, the possible values and defaults are defined by
the $TESTSUITE directory tree. The switches correspond to the capitalized 
directories in $TESTSUITE which contain a DEFAULT file.
Each switch directory contains files named according to the possible
values for the switch. These files contain the actual parameter commands. 
The file DEFAULT must be present in the switch directory and it contains 
the name of the default value for the switch.

The default values are shown as the first choice in the list produced by

   TestBatsrus.pl -l"
#EOC
   ,"\n\n";
    exit;
}

##############################################################################
sub print_hash{

    my %hash = @_;

    print "
switch              values
--------------------------------------------------------------------------
";
    foreach my $switch (sort keys %hash){
	my @option= @{$hash{$switch}};
	printf "%-20s", "-$switch";
	print join(',',@option),"\n";
    }
}
