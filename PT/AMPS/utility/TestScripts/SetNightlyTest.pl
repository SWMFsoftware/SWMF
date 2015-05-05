#!/usr/bin/perl
# The scrpit sets up a nightly tests of AMPS on the current machine
use strict;


# the server name
my $server = 'dborovik@herot.engin.umich.edu';

# the full hostname of the current machine
my $hostname = `hostname -f`;

# home directory
my $home = $ENV{"HOME"};

# the machine is from a list of machines with special scripts
# (Pleiades, Yelowstone, etc.)
my $IsSpecial='';
# the name for the test (if chosen)
our $TestName;

# pfe number on Pleiades to set tests on
my $PFE = "26";

# name of the script to be installed
my $script="";

# read CRONTAB file to check whether the tests have already been installed
my $cron = `crontab -l 2> /dev/null`;

if($cron =~ m/run_test_amps/) {
    print 
	"WARNING: nightly tests may already be set up in your CRONTAB file!\n"
	.
	"Please check/edit it with command\n\tcrontab -e\n"
	.
	"and look for line with execution of script\n\trun_test_maps.*.sh\n"
	.
	"Remove tests\' execution manually if you want to reinstall it.\n"
	.
	"Tests\' installation has been ABORTED!\n";
    return 1;
};

# set MAILTO to empty string unless it has already been set to some value
$cron = "MAILTO=\"\"\n" . $cron unless ($cron =~ m/^MAILTO/);

# add test execution to CRONTAB based on the current machine
if($hostname =~ m/^pfe(.*).nas.nasa.gov/){
    # Pleiades
    $IsSpecial = '1';
    if($1 ne $PFE){
	print 
	    "WARNING: nightly tests on Pleiades should be set on pfe$PFE,\n"
	    .
	    "but currently you are on $hostname."
	    .
	    "SSH to pfe$PFE and try again."
	    .
	    "Tests\' installation has been ABORTED!\n";
	return 1;	
    }
    else {
	# edit the CRONTAB file & choose the script to install
	$cron = $cron . "30 0 * * * \$HOME/bin/run_test_amps.pleiades.sh\n";
	$script = "run_test_amps.pleiades.sh"
    }
}
else{ 
    # any other machine
    # edit the CRONTAB file & choose the script to install
    $cron = $cron . "30 0 * * * \$HOME/bin/run_test_amps.sh\n";
    $script = "run_test_amps.sh"
}

system('mkdir -p $HOME/bin');
system("cp ./utility/TestScripts/$script $home/bin");

# check if a name for a test has been chosen
$TestName = `hostname -s` unless($TestName);
chomp($TestName);
unless($IsSpecial) {
    # read the script's content
    open(my $fin,'<',"$home/bin/$script");
    local $/ = undef;
    my $fcontent = <$fin>;
    close $fin;
    # put a test's name into the content of script
    $fcontent =~ s/Sites\/Current\/\`hostname -s\`/Sites\/Current\/$TestName/g;
    # write the changed content into the script
    open(my $fout,'>',"$home/bin/$script");
    print $fout $fcontent;
    close $fout;
    system("ssh $server mkdir -p Sites/Current/$TestName");
}


# write the new CRONTAB content to a temporary file
open(my $fh,'>','./tmp_crontab');
print $fh $cron;
close $fh;

# set the CRONTAB
`crontab ./tmp_crontab`;

# remove the temporary file
`rm -rf tmp_crontab`;

return 1;
