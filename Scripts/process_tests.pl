#!/usr/bin/perl -s

# Number of days to show on web page
my $nDay = ($n or $nday or 7);

use strict;

my $SiteDir = `pwd`; chop($SiteDir);

# Conditions for merging into stable branch
my $MinScore = 600;  # minimum for total score: most tests ran
my $MinRate  = 0.95; # minimum 95% success rate
my $merge_stable =   # command to merge master into stable branch
    'cd SWMF && share/Scripts/gitall "checkout stable && ' . 
    'git pull --depth=100 && git merge master && git push"';

# weights for each platform to calculate skill scores
my %WeightMachine = (
    "pleiades"     => "1.0",  # ifort pleiades
    "pgi"          => "0.0",  # pgif90 pleiades
    "gfortran"     => "1.0",  # gfortran optimized
    "grid"         => "1.0",  # nagfor debug
    "mesh"         => "1.0",  # nagfor optimized
    );

# Describe machine in the Html table
my %HtmlMachine = (
    "pleiades"     => "ifort<br>pleiades",
    "pgi"          => "pgf90<br>pleiades",
    "gfortran"     => "gfortran<br>optimized",
    "grid"         => "nagfor<br>debug",
    "mesh"         => "nagfor<br>optimized",
    );

# List of platforms
my @machines = sort keys %WeightMachine;

my $ERROR = "ERROR in process_tests.pl";

my $templatefile = "process_tests.html";

my $resfile = "test_swmf.res";
my $htmlfile = "test_swmf.html";
my $logfile = "test_swmf.log";
my $indexfile = "index.html";
my $codefile = "code.diff";
my $manerrorfile = "manual.err";
my $paramerrorfile = "param.err";
my $lastpassfile = "lastpass.txt";

my %tests;
my %resfile;
my %result;

my $day;
my $machine;

# Take last seven days.
my $days = `ls -r -d SWMF_TEST_RESULTS/*/*/* | head -$nDay`;

my @days = split "\n", $days;

# Extract results for all days and convert logfile into an HTML file
foreach $day (@days){

    #print "day = $day\n";
    next unless -d $day;
    chdir $day;
    #print "pwd = ",`pwd`,"\n";

    foreach $machine (@machines){
	my $file = "$machine/$resfile";

	next unless -f $file;

	open RESULTS, $file or die "$ERROR: could not open $file\n";

	my $file2 = "$machine/$htmlfile";

	open HTML, ">$file2";
	print HTML 
	    "<h1>List of ${day}'s difference files for $machine</h1>\n".
	    "<pre>\n";
	while(<RESULTS>){
	    last if /===========/;
            s/Domain Users//;
	    my @item = split(' ',$_);
	    my $test = $item[-1]; $test =~ s/\.diff$//; 
	    $test =~ s/test(\d_)/test0$1/; # rename testN to test0N
	    my $size = $item[4];

	    $tests{$test} = 1;
	    if($size){
		$result{$day}{$test}{$machine}=
		    "<font color=red>failed</font>";
	    }else{
		$result{$day}{$test}{$machine}=
		    "<font color=green>passed</font>";
	    }
	    print HTML $_;
	}
	print HTML 
	    "</pre>\n".
	    "<H1>Head of ${day}'s difference files for $machine</H1><pre>\n";

	# Count number of failed test for the GM/BATSRUS functionality tests
	my $testfunc="GM/BATSRUS/test_func";
	my $test;
	my $stage;
	my $nfail;
	my $isfunc;
	while(<RESULTS>){
	    if( s{==>\ (\S+)\.diff\ <==}
		{</pre><H2><A NAME=$1>$1 ($day: $machine)</H2><pre>}x){
		my $newtest = $1;

		# specify failure for previous test
		# Indicate number of failures for functionality test
		if($isfunc){
		    $stage = "$nfail test";
		    $stage .= "s" if $nfail > 1;
		}
		$result{$day}{$test}{$machine} =~ s/failed/$stage/ if $stage;

		$test = $newtest;
		$stage = "result";
		$isfunc = ($test eq $testfunc);
	    }

	    # read last stage
	    $stage = $1 if /(compile|rundir|run)\.\.\.$/;
	    $stage = "run" if /could not open/;

	    # Count number of failed tests for BATSRUS functionality test suite
	    $nfail++ if $isfunc and /^\-/;

	    print HTML $_;
	}
	# Fix the stage for the last test
	$result{$day}{$test}{$machine} =~ s/failed/$stage/ if $stage;

	close RESULTS;
	close HTML;
    }
    chdir $SiteDir;
}

# update last pass information
our %lastpass;
my $machine;
my $key;
my $day = @days[0];
my $date = $day; $date =~ s/SWMF_TEST_RESULTS\///;

# read in last pass information
require $lastpassfile if -f $lastpassfile;

# update information with last days' results
foreach $machine (@machines){
    my $test;
    foreach $test (sort keys %tests){
	foreach $day (@days){
	    if($result{$day}{$test}{$machine} =~ /passed/i){
		$key = sprintf("%-50s : %-10s", $test, $machine);
		my $date = $day; $date =~ s/SWMF_TEST_RESULTS\///;
		$lastpass{$key} = $date;
		last;
	    }
	}
    }
}

# save latest information
open LASTPASS, ">$lastpassfile" or die "$ERROR: could not open $lastpassfile\n";
print LASTPASS "\%lastpass = (\n";
foreach $key (sort keys %lastpass){
    print LASTPASS "'$key' => '$lastpass{$key}',\n";
}
print LASTPASS ");\n";
close LASTPASS;

# Avoid creating a gigantic index.html file
exit 0 if $nDay > 10;

# Create result table
my $Table = "<hr>\n<center>\n";

my %change;
foreach $day (@days){

    next unless -d $day;

    my $dayname = $day;
    $dayname =~ 
	s/SWMF_TEST_RESULTS\/\d\d\d\d\//Test results for 7pm /;
    $dayname .= '. Score: ';

    my $MaxScores;
    my $Scores;

    # Start table with first row containing the machine names
    $Table .=	
	"<h3>$dayname<FONT COLOR=GREEN>_SCORE_</FONT></h3>\n".
	"<p>\n".
	"<table border=3>\n".
	"  <tr>".
	"   <td><b>compiler/options/platform</b><br><br>test</td>";
    my $machine;
    foreach $machine (@machines){
	my $file = "$day/$machine/$htmlfile";
	$Table .= "    <td><b>$HtmlMachine{$machine}</b></br>\n";
	if(-s $file){
	    $Table .= "    <A HREF=$file TARGET=swmf_test_results>results".
		"</A></br>\n";
	}else{
	    $Table .= "    <font color=red>no results</font></br>\n";
	};
	# Make a row for the log files
    	$file = "$day/$machine/$logfile";
    	if(-s $file){
    	    $Table .= "    <A HREF=$file TARGET=swmf_test_log>log file</A>"
		. "</td>\n";
    	}else{
    	    $Table .= "    <font color=red>no log</font></td>\n";
    	};
    }
    $Table .= "  </tr>\n";
    
    # Print a row for each test
    my $test;
    foreach $test (sort keys %tests){
	$test =~ s/notest/<font color=red>notest<\/font>/;
	$Table .= "<tr>\n  <td>$test</td>\n";
	if($test !~ /notest|see_gm_batsrus/){
	    foreach $machine (@machines){
		my $result = $result{$day}{$test}{$machine};

		my $key = sprintf("%-50s : %-10s", $test, $machine);
		my $lastpass = $lastpass{$key};
		$lastpass = "?/?/?" unless $lastpass;

		if($result){

		    # Assign numeric value to the result
		    my $score; 
		    my $MaxScore;
		    if($result =~ /passed/i){
			$score = 1.0;
		    }elsif($result =~ /result/i){
			$score = 0.5;
		    }elsif($result =~ /run/i){
			$score = 0.1;
		    }elsif($result =~ /(\d+) tests/i){
			$score = (36-$1)/36.0;
		    }else{
			$score = 0;
		    }

		    # Add up results multiplied by various weights
		    my $WeightMachine = $WeightMachine{$machine};
		    $MaxScores += $WeightMachine;
		    $Scores    += $WeightMachine*$score;
		}

		# Indicate changes with upper case
		if($day eq $days[0] or $day eq $days[1]){
		    my $resultnew = $result{$days[0]}{$test}{$machine};
		    my $resultold = $result{$days[1]}{$test}{$machine};

		    $result = uc($result) if $resultnew ne $resultold;
		}

		# Add HTML link for failed tests
		$result = 
		    "<A HREF=\"$day/$machine/$htmlfile\#$test\" ".
		    "target=swmf_test_results title=$lastpass>".
		    $result.
		    "</HREF>" if $result !~ /passed/i;

		$Table .= "<td>$result</td>\n";
	    }
	}
	$Table .= "  </tr>\n";
    }

    # Calculate success rate score and save it into file and table
    my $score = sprintf("%.1f", 100*$Scores/($MaxScores+1e-30)). '%';

    print "day=$day score=$score MaxScores=$MaxScores\n";
    if($MaxScores > $MinScore and $Scores > $MinRate*$MaxScores){
	$Table =~ s/_SCORE_/$score/; # stable score is green
	# Merge last day into stable branch if it has not been done yet
	if($day eq $days[0] and not -f "$day/stable.txt"){
	    print "$merge_stable\n";
	    print `$merge_stable 2>&1`;
	}
	unlink "$day/unstable.txt";
	open SCORE, ">$day/stable.txt";          # indicates stable score
    }else{
	$Table =~ s/GREEN\>_SCORE_/RED\>$score/; # unstable score is red
	unlink "$day/stable.txt";
	open SCORE, ">$day/unstable.txt";        # indicate unstable score
    }
    print SCORE "ALL: $score\n";                 # set/update score
    close SCORE;

    $Table .= "  </tr>\n</table>\n";
}
$Table .= "</center>\n";

open FILE, ">$indexfile" or die "$ERROR: could not open $indexfile\n";

open TEMPLATEFILE, $templatefile;
while(<TEMPLATEFILE>){
    last  if /_TABLES_/;
    print FILE $_;
}

# Check for errors in PARAM.in and PARAM.XML files
foreach my $machine (@machines){
    my $file = $days[0]."/$machine/$logfile";
    next unless -s $file;
    my $Error = `grep -C2 '^Error at line' $file`;
    $Error .= `grep -C2 '^TestParam_ERROR' $file`;
    $Error .= `grep -C2 '^XML ERROR' $file`;
    $Error .= `grep -C2 '/PARAM.in_orig_' $file`;
    next unless $Error;
    open ERR, ">$paramerrorfile";
    print ERR $Error;
    close ERR;    
    print FILE "
<h3><A HREF=$paramerrorfile TARGET=swmf_param_error>
<font color=red>Errors in PARAM.in and/or PARAM.XML files</font></A> See the 
<A HREF=$file TARGET=swmf_test_log>logfile</A> for more info.
</h3>
";
    last;
}

if(-s $manerrorfile){
    print FILE "
<h3><A HREF=$manerrorfile TARGET=swmf_man_error>
<font color=red>Errors in creating the manual.</font></A> See the 
<A HREF=manual.log TARGET=swmf_man_error>manual creation logfile</A> 
for more info.
</h3>
";
}

if(-s $codefile){
    print FILE "
<h3><A HREF=$codefile TARGET=swmf_code_changes>
Source code changed: diff -r SWMF SWMF_yesterday
</A></h3>
";
}

print FILE $Table;

while(<TEMPLATEFILE>){
    print FILE $_;
}
close TEMPLATEFILE;

print FILE "<h3>Weighting scheme</h3>\n";

print FILE "
<b><pre>
Weight - machine
------------------------------------------
";

foreach my $machine (sort keys %WeightMachine){
    print FILE "   $WeightMachine{$machine} - $machine\n";
}


print FILE "</pre></b>\n";
close FILE;

exit 0;
