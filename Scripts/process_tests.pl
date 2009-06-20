#!/usr/bin/perl
use strict;

my $ERROR = "ERROR in process_tests.pl";
my @machines;


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

my @days = 
    ('Today','Yesterday','2days','3days','4days','5days','6days','7days');

# Extract results for both days and convert logfile into an HTML file
foreach $day (@days){

    next unless -d $day;
    chdir $day;

    # The results for each machine are in different directories
    @machines = grep {-d} glob('*');

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
    chdir "..";
}

`cat process_tests.html > $indexfile`;

my $Table = "<hr>\n<center>\n";

my %change;
foreach $day (@days){

    next unless -d $day;

    my $dayname = $day;
    $dayname =~ s/(\d)days/$1 days ago/;

    # Start table with first row containing the machine names
    $Table .=	
	"<h3>$dayname _SCORE_</h3>\n".
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

    my $MaxScore = 0;
    my $score = 0;
    
    # Print a row for each test
    my $test;
    foreach $test (sort keys %tests){
	$test =~ s/notest/<font color=red>notest<\/font>/;
	$Table .= "<tr>\n  <td>$test</td>\n";
	if($test !~ /notest|see_gm_batsrus/){
	    foreach $machine (@machines){
		my $result = $result{$day}{$test}{$machine};

		$MaxScore += 1;
		if($result =~ /passed/){
		    $score += 1
		}

		if($day =~ /Today|Yesterday/){
		    my $resultnew = $result{"Today"}{$test}{$machine};
		    my $resultold = $result{"Yesterday"}{$test}{$machine};

		    # remove HTML (font colors)
		    $resultnew =~ s/<[^>]*>//g; 
		    $resultold =~ s/<[^>]*>//g;

		    $resultnew .= " failed" unless $resultnew =~ /passed|failed/;
		    $resultold .= " failed" unless $resultold =~ /passed|failed/;

		    if($resultnew ne $resultold){
			$change{"$test on $machine     today: $resultnew\n".
				    "$test on $machine yesterday: $resultold\n\n"}++;
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

    my $skill = int(100*$score/$MaxScore);
    $Table =~ s/_SCORE_/passed\/total=$score\/$MaxScore=$skill\%/;

    $Table .= "  </tr>\n</table>\n";
}
$Table .= "</center>\n";

open FILE, ">$changefile" or die "$ERROR: could not open $changefile\n";
print FILE sort keys %change;
close FILE;

open FILE, ">>$indexfile" or die "$ERROR: could not open $indexfile\n";

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
Summary of test differences between Today and Yesterday
</A></h3>
";
}

print FILE $Table;
close FILE;

