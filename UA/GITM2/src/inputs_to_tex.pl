#!/usr/bin/perl -s

$help = $h;

open(INFILE,"<set_inputs.f90");

open(OUTFILE,">set_inputs.tex");

while (<INFILE>) {

    if (/Incorrect format/) {

	$IsDone=0;
	while (!$IsDone) {
	    $line = <INFILE>;
	    if (!($line =~ /\#/)) {
		if ($line =~ /(write\(\*,.*\)\s)[\'\"](.+)[\'\"]/) {
		    print OUTFILE "$2\n";
		}
	    } else {
		$IsDone = 1;
	    }
#	    print OUTFILE "$2\n" if (!$IsDone);
	}

	print OUTFILE "\\begin{verbatim}\n";
	print OUTFILE "$2\n" if ($line =~ /(write\(\*,.*\)\s)[\'\"](.+)[\'\"]/);
	$IsDone=0;
	while (!$IsDone) {

	    $line = <INFILE>;
	    $IsDone=1 if (!($line =~ /(write\(\*,.*\)\s)[\'\"](.+)[\'\"]/));
	    print OUTFILE "$2\n" if (!$IsDone);

	}

	print OUTFILE "\\end{verbatim}\n";
	print OUTFILE "\n";

    }

}

close(INFILE);
close(OUTFILE);

exit(1);
