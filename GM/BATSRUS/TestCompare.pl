#!/usr/bin/perl -s
#^CFG COPYRIGHT UM
#^CFG FILE TESTING

$Help=$h; &print_help if $Help;

$Tar=$tar;

if($Tar){
    # tar up a directory
    &print_help if $#ARGV != 0;
    $dir = $ARGV[0];
    $list = 'test.*/log.* test.*/error.* test.*/SWITCHES test.*/PARAM.expand '.
	' test.*/IO2/*.log test.*/IO2/*.sat test.*/ionosphere';
    print "cd $dir;\ntar -cf ../TEST.tar $list\n";
    print `cd $dir;  tar -cf ../TEST.tar $list`;
    exit 0;
}

&compare_speed if $speed;

# Normal use with two directories
&print_help if $#ARGV != 1;
$Strict=$s;
$Quiet =$q;
$Omit  =$o;
$Verbose=$v;
$tol=0.001 unless $tol;
$dir1=$ARGV[0];
$dir2=$ARGV[1];

opendir(DIR1, $dir1) or die "ERROR: cannot open directory $dir1 !\n";
opendir(DIR2, $dir2) or die "ERROR: cannot open directory $dir2 !\n";

@test1 = grep /test\.\d\d\d/, readdir DIR1;
@test2 = grep /test\.\d\d\d/, readdir DIR2;

die "ERROR no test.xxx directories in $dir1 !\n" unless @test1;
die "ERROR no test.xxx directories in $dir2 !\n" unless @test2;
warn "WARNING different number of test directories in $dir1 and $dir2 !\n"
    unless $#test1 == $#test2;

if($#test1 >= $#test2){
    $dir =$dir1;
    @test=@test1;
}else{
    $dir =$dir2;
    @test=@test2;
}
print "Number of tests to compare is ",$#test+1,"\n";

@parfile=("SWITCHES", "PARAM.expand");

$pwd=`pwd`; chop $pwd; # store current directory

TEST: foreach $test (sort @test){

    $anyerror=0;
    foreach $testdir ("$dir1/$test","$dir2/$test"){
	if(not -d $testdir){
	    print "$test: $testdir is missing!\n";
	    $anyerror=1;
	}
    }
    next TEST if $anyerror;

    # Compare parameter files

  PARFILE: foreach $parfile (@parfile){
	$parfile1="$dir1/$test/$parfile";
	$parfile2="$dir2/$test/$parfile";

	close(PAR1);
	close(PAR2);
	$ok1=open (PAR1, $parfile1);
	$ok2=open (PAR2, $parfile2);

	print "$test: $parfile1 could not be opened !\n" unless $ok1;
	print "$test: $parfile2 could not be opened !\n" unless $ok2;
	next TEST unless $ok1 and $ok2;

	if($parfile eq "SWITCHES"){
	    $switch1=<PAR1>; $switch1 =~ s/^test.*pl //i;
	    $switch2=<PAR2>; $switch2 =~ s/^test.*pl //i;

	    if($switch1 ne $switch2){
		&report("$test: Different options in SWITCHES:\n",
		      $switch1,$switch2) unless $Quiet;
		next TEST if $Strict;
	    }
	    next PARFILE;
	}

	undef @line1;
	while($line1=<PAR1>){
	    push(@line1,$line1) unless $line1 =~ /(^\s*$|^[\!\^])/;
	} 

      LINE: while($line2=<PAR2>){
	    next LINE if $line2 =~ /(^\s*$|^[\!\^])/;
	    if(not $line1 = shift(@line1)){
		&report("$test: $parfile1 has fewer lines than $parfile2\n")
		    unless $Quiet;
		next TEST if $Strict;
		last PARFILE;
	    }
	    if($line1 ne $line2){
	        &report("$test: $parfile files differ at line $.:\n",
	          "    $line1    $line2\n")  unless $Quiet;
		next TEST if $Strict;
		last PARFILE;
	    }
	}
	if(@line1){
	    &report("$test: $parfile1 has more lines than $parfile2\n")
		unless $Quiet;
	    next TEST if $Strict;
	    last PARFILE;
	}
    }
    close(PAR1);
    close(PAR2);

  ERRFILE: foreach $errorstar ("$dir1/$test/error.*","$dir2/$test/error.*"){
	@errorfile=glob($errorstar);
	if($#errorfile != 0){
	    &report("$test: no $errorstar found!\n");
	    next ERRFILE;
	}
	$errorfile=$errorfile[0];
	if(not open(ERR, $errorfile)){
	    &report("$test: $errorfile could not be opened !\n");
	    next ERRFILE;
	}
	$error = <ERR>;
	close(ERR);
	if(length($error)>0){
	    &report("$test: $errorfile=",$error);
	    $anyerror=1;
	}
    }
    next TEST if $anyerror and $Strict;

    if(not $Quiet){
	# Compare number and size of files in the IO2 directories
	chdir "$dir1/$test/IO2";
	$wc1 = `wc -c *`;
	chdir "$pwd/$dir2/$test/IO2";
	$wc2 = `wc -c *`;
	chdir $pwd;
	if($wc1 ne $wc2){
	    @wc1=split(/\n/,$wc1);
	    @wc2=split(/\n/,$wc2);
	    if($#wc1 != $#wc2){
		&report("$test !!! number of files in IO2-s differ !!!\n");
		&report("$dir1/$test/IO2:\n$wc1$dir2/$test/IO2:\n$wc2");
	    }else{
		# Report size differences except for compressed TEC files 
		# (.dat.gz) which may differ due to small differences 
		# in the results
		for($i=0; $i<$#wc1-1; $i++){
		    if($wc1[$i] ne $wc2[$i] and $wc1[$i]!~/dat\.gz/){
			&report("$test: IO2 files differ:\n",
				$wc1[$i],"\n",$wc2[$i],"\n");
		    }
		}
	    }
	}
    }

    # Compare logfiles
    chdir "$dir/$test";
    @logfile=glob("IO2/*.log IO2/*.sat ionosphere/*.idl ionosphere/*.tec");
    chdir $pwd;


  LOGFILE: foreach $logfile (@logfile){

	$logfile1="$dir1/$test/$logfile";
	$logfile2="$dir2/$test/$logfile";

	close(LOG1);
	close(LOG2);
	$ok1=open(LOG1, $logfile1);
	$ok2=open(LOG2, $logfile2);

        # This logfile is missing from both
	next LOGFILE if not $ok1 and not $ok2; 

	&report("$test: $logfile1 could not be opened !\n") unless $ok1;
	&report("$test: $logfile2 could not be opened !\n") unless $ok2;
	next LOGFILE unless $ok1 and $ok2;

	$n1=`wc -l $logfile1`;
	$n2=`wc -l $logfile1`;

	next LOGFILE if $n1 == 0 and $n2 == 0; # Empty logfiles

	if($n1 != $n2){
	    &report("$test: !!! $logfile number of lines $n1 /= $n2 !!!\n");
	    next LOGFILE;
	}

	$diffmax=0;

	if ($logfile =~ /iono/){
	    $nheadline=100;
	    @logvar=("Theta","Psi","SigmaH","SigmaP","Jr","Phi");
	}else{
	    $nheadline=2;
	}


        LINE: while($line1=<LOG1> and $line2=<LOG2>){
	    if($.<=$nheadline){
		$line1 =~ s|GM/Param/TESTSUITE|Param/TESTSUITE|;
		$line2 =~ s|GM/Param/TESTSUITE|Param/TESTSUITE|;
		$line1 =~ s/\-(0\.0+)\b/ $1/;
		$line2 =~ s/\-(0\.0+)\b/ $1/;
		if($line1 ne $line2){
		    &report("$test $logfile, different head lines:\n",
			    "    $line1    $line2\n");
		    next LOGFILE if $Strict;
		}
		
		@logvar=split(' ',$line1) if $logfile !~ /iono/ and $.==2;
		$nheadline = $. if $line1 =~ /^BEGIN|ZONE/;
		next LINE;
	    }
	    @item1=split(' ',$line1);
	    @item2=split(' ',$line2);
	    if($#item1 != $#item2){
		&report("$test1 at line $. number of columns differ:\n",
			"    $line1    $line2\n");
		next LOGFILE;
	    }
	    for($i=0; $i<=$#item1; $i++){
		$item1=$item1[$i];
		$item2=$item2[$i];
		$diff = abs($item1-$item2);
		if ($diff > $tol){
		    &report("$test $logfile: diff=$diff for ",$logvar[$i],
			    " at line $. !!!\n");
		    next LOGFILE;
		}
		if($diff > $diffmax){
		    $diffmax=$diff;
		    $linenum=$.;
		    $itemnum=$i;
		}
	    }
	}
	close(LOG1);
	close(LOG2);

	if($diffmax >= $Omit){
	    if($diffmax>0){
		&report("$test $logfile: diffmax=$diffmax for ",
			$logvar[$itemnum]," at line $linenum\n");
	    }else{
		&report("$test $logfile: identical\n");
	    }
	}
    }
}

exit;

##############################################################################
sub report{

    if($Verbose and $switch1){
	print "\n$switch1";
	undef $switch1;
    }
    print @_;
}

##############################################################################

sub compare_speed{

    &print_help if $#ARGV == -1;

    my @dir=@ARGV; # list of directories

    my $What = $speed; $What = 'BATSRUS\|SWMF' if $What eq "1"; #speed of what?

    my @speedsum; # sum of speeds for tests completed in all directories
    my $speedsum; # number of tests completed in all directories

    print "Timings for $What\n";
    print "tst";
    foreach (@dir){printf "%8s",substr($_,0,7)};
    print "     diff" if $#dir==1;
    print "   switches\n";
    print "=" x 79,"\n";
    foreach $number ("000".."050"){
	my $test = "test.$number";

	my @speed;
	my $switch;
	my $missing;

	foreach $dir (@dir){
	    if(-d "$dir/$test"){
		my $myswitch;
		$myswitch = `head -1 $dir/$test/SWITCHES`; chop($myswitch);
		$myswitch =~ s/Test(Batsrus|SWMF).pl\s*//;
		if(not $switch){
		    $switch = $myswitch;
		}else{
		    warn "WARNING: switch=$switch ne myswitch=$myswitch ".
			"in $dir/$test/SWITCHES\n"
			if $switch ne $myswitch;
		}
		$_ = `grep "$What" $dir/$test/log* | tail -1`;
		($speed) = /(\d+\.\d\d)/;
		$speed   = sprintf("%8.2f",$speed);
	    }elsif(-f "$dir/log.$number"){
		$_ = `grep "$What" $dir/log.$number | tail -1`;
		($speed) = /(\d+\.\d\d)/;
		$speed   = sprintf("%8.2f",$speed);
	    }else{
		$speed = "    --- "; $missing++;
	    }
	    push @speed, $speed;
	}
	next if $missing > $#dir; # nothing to report

	print "$number",join("",@speed);
	if($#dir==1){
	    # Calculate speed difference for 2 directories
	    if($speed[0] > 0 and $speed[1] > 0){
		printf "%9.2f",$speed[1]-$speed[0];
	    }else{
		print "     --- ";
	    }
	}
	if(length($switch)>70){
	    print "   ",substr($switch,0,70),"...\n";
	}else{
	    print "   $switch\n";
	}

	next if $missing; # do not add up if any result is missing

	for $i (0..$#dir){$speedsum[$i] += $speed[$i]} # add up speeds
	$speedsum++;
    }
    if(@speedsum){
	print "-" x 79,"\n";
	print "Sum";
	for $i (0..$#dir){
	    printf "%8.2f",$speedsum[$i];
	}
	if($#dir==1){
	    # Calculate speed difference for 2 directories
	    if($speedsum[1]>0 and $speedsum[1]>0){
		printf "%9.2f",$speedsum[1]-$speedsum[0];
	    }else{
		print "     --- ";
	    }
	}
	print "   for $speedsum compeleted tests\n";
    }
    exit 0;
}

##############################################################################
#BOP
#!ROUTINE: TestCompare.pl - a tool for quantitative comparison of test results
#!DESCRIPTION:
# This script compares two sets of test results obtained with different
# versions or on different platforms. 

#!REVISION HISTORY:
# 01/27/02 G.Toth - initial version developed for BATSRUS
# 07/10/03 G.Toth - adapted to SWMF
# 04/01/04 G.Toth - added -speed option to compare speeds
#EOP

sub print_help{

    print "
Purpose: 
Compare log files between complete test suites.
Compare execution speeds between test suites.
Tar up important log files from a test suite.

",
#BOC
"Usage:

TestCompare.pl [-h] [-q] [-s] [-v] [-o=L] [-tol=X] Dir1 Dir2

or

TestCompare.pl -speed[=WHAT] Dir1 [Dir2] [Dir3] ...

or

TestCompare.pl -tar Dir

-h      Help message.
-q      Quiet: suppress messages about less important differences.
-s      Strict mode: insist on identical parameter files.
-v      Verbose: show the non-default options for each test.
-o=L    Omit tests with max difference below L, default is 0.
-tol=X  Print difference above X as an error, default is 0.001

Dir1    Directory 1 containing test runs
Dir2    Directory 2 containing test runs

-speed[=WHAT]
        Extract and list execution speeds for all tests of the test suite.
        The timings will be compared for 'WHAT' which is the (abbreviated) 
        name found in the timing report. Default value is 'BATSRUS|SWMF'.
        If called with more than one directory the speeds for a given
        test are listed next to each other which makes comparison easy.
        When 2 directories are compared, the difference is also printed.

Dir1,Dir2,Dir2...    
        The directories containing test.0xx/log.0xx or log.0xx files,
        where xx goes from 00 to 35.

-tar    Tar up files relevant for comparison into TEST.tar
Dir     Directory containing test runs to be tarred up

Examples:

      TestCompare.pl run1 run2

      TestCompare.pl -q -v -o=1e-12 run1 run2

      TestCompare.pl -speed run1 run2

      TestCompare.pl -speed=exch_msgs run1 run2 run3

      TestCompare.pl -tar run"
#EOC
,"\n\n";
    exit;
}
