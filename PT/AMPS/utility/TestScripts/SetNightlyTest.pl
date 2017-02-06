#!/usr/bin/perl
# The scrpit sets up a nightly tests of AMPS on the current machine
use strict;

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

#------------------------------------------------------------------------------
# Proceed with tests' installation

# the server name
my $server = 'vtenishe@herot.engin.umich.edu';

# the full hostname of the current machine
my $hostname = `hostname -f`;

# home directory
my $home = $ENV{"HOME"};

# pfe number on Pleiades to set tests on
my $PFE = "26";

# login number on Stampede to set tests on
my $SLOGIN = "1";

# date on the machine
my $date = `date`;

# hash for installing a test time based on time zone
my %TZ2CRON_run = (
    "EDT" => "30 23", "EST" => "30 23", 
    "CDT" => "30 22", "CST" => "30 22",
    "MDT" => "30 21", "MST" => "30 21",
    "PDT" => "30 20", "PST" => "30 20");
my %TZ2CRON_scp = (
    "EDT" => "30 10", "EST" => "30 10", 
    "CDT" => "30 9", "CST" => "30 9",
    "MDT" => "30 8", "MST" => "30 8",
    "PDT" => "30 7", "PST" => "30 7");

# name of the script to be installed
my $run_script="run_test_amps.sh";
my $scp_script="scp_test_amps.sh";

# the name for the test (if chosen)
our $TestName;
# check if a name for the test has been chosen
$TestName = `hostname -s` unless($TestName);
chomp($TestName);


# list of compilers to use
our @Compilers;

# list of supported compilers' 
#naming is case sensitive and the same as in %Blocks hash table
my @CompilersSupported = ("GNU", "Intel","PGI");

# HASH TABLE with blocks in the $run_script
# the body of the named script is divided into blocks of the format
#  #>BlockName ##############
#  #command_1               #
#  # ...                    #
#  #command_last           <#
# certain blocks should be uncommented at the test installation,
# e.g. if test is GNU compiled, then block #>GNU <# should be uncommented,
#      if test is set on Pleiades, then #>Pleiades <# should be uncommented

# the table below is initialized to 0 for each block
# further, this scripts determines which blocks should be activated
my %Blocks = (
    "GNUAll"       => "0",
    "GNUOther"     => "0",
    "IntelAll"     => "0",
    "IntelOther"   => "0",
    "PGIAll"       => "0",
    "PGIOther"     => "0",
    "Pleiades"     => "0",
    "Yellowstone"  => "0",
    "Stampede"     => "0",
    "Valeriy"      => "0"
    );

# fill the hash table %Blocks
if($hostname =~ m/^pfe(.*).nas.nasa.gov/){
    # Pleiades
    $Blocks{'Pleiades'} = "1";
    $TestName = "pleiades";
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
}
elsif($hostname =~ m/^yslogin(.*)/){ 
    # Yellowstone
    $Blocks{"Yellowstone"} = "1";
    $TestName = "yellowstone";
}
elsif($hostname =~ m/login(.).stampede.tacc.utexas.edu/){
    # Stampede
    $Blocks{"Stampede"} = "1";
    $TestName = "stampede";
    if($1 ne $SLOGIN){
	print 
	    "WARNING: nightly tests on Stampede should be set on login$SLOGIN,\n"
	    .
	    "but currently you are on $hostname."
	    .
	    "SSH to login$SLOGIN and try again."
	    .
	    "Tests\' installation has been ABORTED!\n";
	return 1;	
    }
}
elsif($hostname =~ m/srbwks2014-0079.engin.umich.edu/){
    # Valeriy
    $Blocks{"Valeriy"} = "1";
    $TestName = "valeriy"
}

# set hashes for compilers (already converted to lower case)
if($Compilers[0] eq "all"){
    foreach my $CompilerSupported(@CompilersSupported){
	$Blocks{$CompilerSupported."All"  } = "1";
	unless($Blocks{"Pleiades"} or 
	       $Blocks{"Yellowstone"} or
	       $Blocks{"Stampede"}){
	    $Blocks{$CompilerSupported."Other"  } = "1";
	}
    }
}
else{
  LoopCompiler:
    foreach my $Compiler (@Compilers){
	foreach my $CompilerSupported(@CompilersSupported){
	    if ($Compiler eq lc($CompilerSupported)){
		$Blocks{$CompilerSupported."All"  } = "1";
		unless($Blocks{"Pleiades"} or 
		       $Blocks{"Yellowstone"} or
		       $Blocks{"Stampede"}){
		    $Blocks{$CompilerSupported."Other"  } = "1";
		}
		next LoopCompiler;
	    }
	}
	die "ERROR: Compiler ".$Compiler." is not recognized!";
    }
}

# copy $run_script to $HOME/bin
system('mkdir -p $HOME/bin');
system("cp ./utility/TestScripts/$run_script $home/bin");
system("cp ./utility/TestScripts/$scp_script $home/bin");

&prepare_script("$home/bin/$run_script");
&prepare_script("$home/bin/$scp_script");

#------------------------------------------------------------------------------
# Set the CRONTAB script

# set MAILTO to empty string unless it has already been set to some value
$cron = "MAILTO=\"\"\n" . $cron unless ($cron =~ m/^MAILTO/);

# set scripts execution
my $TZ = `date +%Z`; chomp($TZ);
$cron .= $TZ2CRON_run{$TZ}." * * * \$HOME/bin/$run_script\n";
$cron .= $TZ2CRON_scp{$TZ}." * * * \$HOME/bin/$scp_script\n";


# write the new CRONTAB content to a temporary file
open(my $fh,'>','./tmp_crontab');
print $fh $cron;
close $fh;

# set the CRONTAB
`crontab ./tmp_crontab`;

# remove the temporary file
#`rm -rf tmp_crontab`;



#------------------------------------------------------------------------------
# create test's folder on the server
if($Compilers[0] eq "all"){
    foreach my $CompilerSupported (@CompilersSupported){
	system("ssh $server mkdir -p Sites/Current/$TestName\_".
	       lc($CompilerSupported));
    }
}
else{
    foreach my $Compiler (@Compilers){
	system("ssh $server mkdir -p Sites/Current/$TestName\_$Compiler");
    }
}

return 1;


sub prepare_script {
# introduce changes in $run_script based on hash table $Blocks                 
    foreach my $Block(keys %Blocks){
        next unless ($Blocks{$Block});
        open(my $fh,'<',$_[0]);
        my @lines = <$fh>;
        close($fh);
        my $IsInBlock ='0';
        open(my $fh,'>',$_[0]);
        foreach my $line (@lines){
            if($line =~ m/>$Block/){
                $IsInBlock = '1';
            }
            elsif($IsInBlock){
		if ($line =~ m/<#/){$IsInBlock ='0'; $line =~ s/<#//g;}
                $line =~ s/^#//g;
            }
	    # put a test's name into the content of script 
	    $line =~ s/\/Current\/\`hostname -s\`/\/Current\/$TestName/g;

            print $fh $line;
        }
        close($fh);
    }
    # make script file executable
    system("chmod u+x $_[0]");
    
}
