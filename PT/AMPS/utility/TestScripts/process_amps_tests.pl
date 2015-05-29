#!/usr/bin/perl
use strict;


my $ERROR = "ERROR in process_tests.pl";
my @machines;

my $templatefile = "template.html";

my $resfile = "test_amps.res";
my $htmlfile = "test_amps.html";
my $logfile = "test_amps.log";
my $indexfile = "index.html";
#my $changefile = "test.diff";
#my $codefile = "code.diff";
my $manerrorfile = "manual.err";

my %tests;
my %resfile;
my %result;

my $day;
my $machine;

# Take last seven days.
my $days = `ls -r -d AMPS_TEST_RESULTS/*/*/* | head -7`;
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
        my $test;
        my $stage;
        while(<RESULTS>){
            if( s{==>\ (\S+)\.diff\ <==}
                {</pre><H2><A NAME=$1>$1 ($day: $machine)</H2><pre>}x){
                my $newtest = $1;
                # specify failure for previous test
                $result{$day}{$test}{$machine} =~ s/failed/$stage/ if $stage;

                $test = $newtest;
                $stage = "result";
            }

            # read last stage
            $stage = $1 if /(compile|rundir|run)\.\.\.$/;
            $stage = "run" if /could not open/;

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
	s/AMPS_TEST_RESULTS\/\d\d\d\d\//Test results and scores for midnight /;

    my %MaxScores;
    my %Scores;

    # Start table with first row containing the machine names
    $Table .=	
	"<h3>$dayname</h3>\n".
	"<p>\n".
	"<table border=3>\n".
	"  <tr>".
	"   <td>test / machine</td>";
    my $machine;
    foreach $machine (@machines){
	my $file = "$day/$machine/$htmlfile";
	$Table .= "    <td><b>$machine</b></br>\n";
	if(-s $file){
	    $Table .= "    <A HREF=$file TARGET=amps_test_results>results".
		"</A></br>\n";
	}else{
	    $Table .= "    <font color=red>no results</font></br>\n";
	};
	# Make a row for the log files
    	$file = "$day/$machine/$logfile";
    	if(-s $file){
    	    $Table .= "    <A HREF=$file TARGET=amps_test_log>log file".
    		"</A></td>\n";
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
		    "target=amps_test_results>".
		    $result.
		    "</HREF>" if $result !~ /passed/i;

		$Table .= "<td>$result</td>\n";
	    }
	}
	$Table .= "  </tr>\n";
    }

    $Table .= "  </tr>\n</table>\n";
}
$Table .= "</center>\n";


open FILE, ">$indexfile" or die "$ERROR: could not open $indexfile\n";

open TEMPLATEFILE, $templatefile;
while(<TEMPLATEFILE>){
    last  if /_TABLES_/;
    print FILE $_;
}

if(-s $manerrorfile){
    print FILE "
<h3><A HREF=$manerrorfile TARGET=amps_man_error>
<font color=red>Errors in creating the manual.</font></A> See the 
<A HREF=manual.log TARGET=amps_man_error>manual creation logfile</A> 
for more info.
</h3>
";
}

print FILE $Table;

while(<TEMPLATEFILE>){
    print FILE $_;
}
close TEMPLATEFILE;

print FILE "</pre></b>\n";
close FILE;

exit 0;
