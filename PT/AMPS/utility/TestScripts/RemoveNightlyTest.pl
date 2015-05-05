#!/usr/bin/perl
# The scrpit sets up a nightly tests of AMPS on the current machine
use strict;


# the server name
my $server = 'dborovik@herot.engin.umich.edu';

# the full hostname of the current machine
my $hostname_s = `hostname -s`;
chomp($hostname_s);

# home directory
my $home = $ENV{"HOME"};

# the name of the test as installed on the current machine
my $TestName;

# name of the script
my $script="";

# remove the test from the CRONTAB file
my $cron = `crontab -l 2> /dev/null`;
$cron =~ s/\n(.*)run_test_amps(.*)\n/\n/g or return 1;

# script name on this machine
$script = "run_test_amps$2";

{
# read the script's content
    open(my $fin,'<',"$home/bin/$script") or die "ERROR: nightly test script is not found at $home/bin!";
    local $/ = undef;
    my $fcontent = <$fin>;
    close $fin;
# extract the test's name from the content of script
    $fcontent =~ m/Sites\/Current\/(.*)\//;
    $TestName = $1;
    if($TestName eq '`hostname -s`'){
	$TestName = $hostname_s;
    }
    unless($TestName) {die "ERROR: Test name cannot be recovered from $home/bin/$script!";}
    system("ssh $server rm -rf Sites/Current/$TestName");
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
