#!/usr/bin/perl
use strict;

my %WeightMachine = (
    "pleiades"     => "1.0",
    "columbia"     => "1.0",
    "grendel"      => "0.5",
    "grid"         => "1.0",
    "mesh"         => "1.0",
    "nyx"          => "0.2",
    "nyx_pgf90"    => "1.0",
    "xena"         => "0.1",
    );

my @ScoreTypes = ("ALL", "CCHM", "CWMM", "CRASH");
my %WeightTest = (

    "CCHM:test2_sc"                       => "1.0",
    "CCHM:test2_ih"                       => "1.0",
    "CCHM:test_lc"                        => "1.0",
    "CCHM:GM/BATSRUS/test_corona"         => "1.0",

    "CWMM:test1_gm"                       => "1.0",
    "CWMM:test1_ie"                       => "1.0",
    "CWMM:test1_pw"                       => "0.5",
    "CWMM:test3_gm"                       => "1.0",
    "CWMM:test3_ie"                       => "1.0",
    "CWMM:test3_im"                       => "1.0",
    "CWMM:test4_gm"                       => "1.0",
    "CWMM:test4_ie"                       => "1.0",
    "CWMM:test4_im"                       => "1.0",
    "CWMM:test_pw"                        => "0.5",
    "CWMM:test_rb"                        => "0.5",
    "CWMM:GM/BATSRUS/test_magnetometer"   => "0.1",
    "CWMM:PW/PWOM/test_Earth"             => "0.5",
    "CWMM:RB/RBE/test"                    => "0.5",

    "CRASH:GM/BATSRUS/test_eosgodunov"    => "1.0",
    "CRASH:GM/BATSRUS/test_graydiffusion" => "1.0", 
    "CRASH:GM/BATSRUS/test_hyades2d"      => "1.0",
    "CRASH:GM/BATSRUS/test_levelset"      => "1.0",

    );

my $ERROR = "ERROR in process_tests.pl";
my @machines;

my $templatefile = "process_tests.html";

my $resfile = "test_swmf.res";
my $htmlfile = "test_swmf.html";
my $logfile = "test_swmf.log";
my $indexfile = "index.html";
my $changefile = "test.diff";
my $codefile = "code.diff";
my $manerrorfile = "manual.err";

my %tests;
my %resfile;
my %result;

my $day;
my $machine;

# Take last seven days.
my $days = `ls -r -d SWMF_TEST_RESULTS/*/*/* | head -7`;
my @days = split "\n", $days;

# Extract results for all days and convert logfile into an HTML file
foreach $day (@days){

    next unless -d $day;
    chdir $day;

    # The results for each machine are in different directories
    @machines = grep {-d} glob('*') unless @machines;

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
	    my @item = split(' ',$_);
	    my $test = $item[-1]; $test =~ s/\.diff$//;
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
    chdir "../../../..";
}

my $Table = "<hr>\n<center>\n";

my %change;
foreach $day (@days){

    next unless -d $day;

    my $dayname = $day;
    $dayname =~ 
	s/SWMF_TEST_RESULTS\/\d\d\d\d\//Test results and scores for 7pm /;

    my %MaxScores;
    my %Scores;

    # Start table with first row containing the machine names
    $Table .=	
	"<h3>$dayname<br><FONT COLOR=GREEN>_SCORE_</FONT></h3>\n".
	"<p>\n".
	"<table border=3>\n".
	"  <tr>".
	"   <td>test / machine</td>";
    my $machine;
    foreach $machine (@machines){
	my $file = "$day/$machine/$htmlfile";
	$Table .= "    <td><b>$machine</b></br>\n";
	if(-s $file){
	    $Table .= "    <A HREF=$file TARGET=swmf_test_results>results".
		"</A></br>\n";
	}else{
	    $Table .= "    <font color=red>no results</font></br>\n";
	};
	# Make a row for the log files
    	$file = "$day/$machine/$logfile";
    	if(-s $file){
    	    $Table .= "    <A HREF=$file TARGET=swmf_test_log>log file".
    		"</A></td>\n";
    	}else{
    	    $Table .= "    <font color=red>no log</font></td>\n";
    	};
    }
    $Table .= "  </tr>\n";

    #$Table .= "  <tr>\n<td>log files</td>\n";
    #foreach $machine (@machines){
    #}

    
    # Print a row for each test
    my $test;
    foreach $test (sort keys %tests){
	$test =~ s/notest/<font color=red>notest<\/font>/;
	$Table .= "<tr>\n  <td>$test</td>\n";
	if($test !~ /notest|see_gm_batsrus/){
	    foreach $machine (@machines){
		my $result = $result{$day}{$test}{$machine};

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
		    $MaxScores{"ALL"} += $WeightMachine;
		    $Scores{"ALL"}    += $WeightMachine*$score;
		    foreach my $type (@ScoreTypes){
			$MaxScore = $WeightMachine*$WeightTest{"$type:$test"};
			$MaxScores{$type} += $MaxScore;
			$Scores{$type}    += $MaxScore*$score;
		    }
		}

		if($day eq $days[0] or $day eq $days[1]){
		    my $resultnew = $result{$days[0]}{$test}{$machine};
		    my $resultold = $result{$days[1]}{$test}{$machine};

		    # remove HTML (font colors)
		    $resultnew =~ s/<[^>]*>//g; 
		    $resultold =~ s/<[^>]*>//g;

		    $resultnew .= " failed" 
			unless $resultnew =~ /passed|failed/;
		    $resultold .= " failed" 
			unless $resultold =~ /passed|failed/;

		    if($resultnew ne $resultold){
			$change{"$test on $machine     today: $resultnew\n"
				    ."$test on $machine yesterday: $resultold\n\n"}++;
			$result = uc($result);
		    }
		}

		# Add HTML link for failed tests
		$result = 
		    "<A HREF=\"$day/$machine/$htmlfile\#$test\" ".
		    "target=swmf_test_results>".
		    $result.
		    "</HREF>" if $result !~ /passed/i;

		$Table .= "<td>$result</td>\n";
	    }
	}
	$Table .= "  </tr>\n";
    }

    # Calculate score per centages and save them into file and table
    my $score; 
    foreach my $type (@ScoreTypes){
	$score .= "$type: " 
	    . sprintf("%.1f", 100*$Scores{$type}/($MaxScores{$type}+1e-30))
	    . '%, ';
    }
    chop($score); chop($score);
    $Table =~ s/_SCORE_/$score/;

    open SCORE, ">$day/all_scores.txt";
    $score =~ s/, /\n/g;
    print SCORE "$score\n";
    close SCORE;

    $Table .= "  </tr>\n</table>\n";
}
$Table .= "</center>\n";

open FILE, ">$changefile" or die "$ERROR: could not open $changefile\n";
print FILE sort keys %change;
close FILE;

open FILE, ">$indexfile" or die "$ERROR: could not open $indexfile\n";

open TEMPLATEFILE, $templatefile;
while(<TEMPLATEFILE>){
    last  if /_TABLES_/;
    print FILE $_;
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

if(-s $changefile){
    print FILE "
<h3><A HREF=test.diff TARGET=swmf_test_summary>
Summary of test differences between $days[0] and $days[1]
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


my $oldtype;
foreach my $typetest (sort keys %WeightTest){
    my $test = $typetest;
    $test =~ s/(\w+)://;
    my $type = $1;
    if($type ne $oldtype){
	print FILE "
PROJECT $type:
   weight - test
   ---------------------------------------
";
	$oldtype = $type;
    }
    print FILE "      $WeightTest{$typetest} - $test\n";
}
print FILE "</pre></b>\n";
close FILE;

exit 0;
