#!/usr/bin/perl
# The scrpit builds Makefile.test specific for the current machine
use strict;

# the hostname of the current machine
my $hostname = `hostname -s`;
chomp($hostname);

#path to the Makefile.test source
my $path="MakefileTest";

#table of the nightly tests
my $fTable="$path/Table";
#content of the table
my @Table=read_content($fTable);

#base for the Makefile.test
my $fBase="$path/makefile.test.base";
my @Base=read_content($fBase);

#individual tests for the Makefile.test
my $fApp="$path/makefile.test.app";
my @App=read_content($fApp);

my $fFinal="Makefile.test";
my @Final;
my @FinalApps;


#parameters of tests
my $Name;
my $Keys;
my $Outs;

#process table with test description
while(@Table){
    my $ref;
    ($ref,$Name,$Keys,$Outs) = get_next_test(@Table);
    @Table = @$ref;
    $Outs="test_$Name" unless($Outs);
    next unless($Name);

    for (my $i = 0; $i<=@Base-1; $i++){
	# general part with all tests
	if($Base[$i]=~m/<APP>/){
	    $Final[$i]=$Final[$i].$Base[$i];
	    $Final[$i]=~s/<APP>/$Name/g;
	    $Final[$i]=~s/<APPKEYS>/$Keys/g;
	    $Final[$i]=~s/<APPOUTS>/$Outs/g;
	}
	else{
	    $Final[$i]=$Base[$i];
	}
    }
    # application specific blocks
    my @lines = @App;
    for (my $i = 0; $i<=@lines-1; $i++){
	$lines[$i]=~s/<APP>/$Name/g;
	$lines[$i]=~s/<APPKEYS>/$Keys/g;
	$lines[$i]=~s/<APPOUTS>/$Outs/g;
    }
    push(@FinalApps,@lines);
}

#build the final Makefile
push(@Final,"\n\n");
push(@Final,@FinalApps);
#write it
&write_content($fFinal,@Final);

return 1;


sub read_content {
    # read the content of a file
    open(my $fh,'<',$_[0]);
    my @lines = <$fh>;
    close($fh);
    @lines;
}

sub write_content {
    # write the content of a file
    open(my $fh,'>',shift(@_));
    print $fh @_;
    close($fh);
}


sub get_next_test{
    #returns test description extracted from content of fTable
    my $IsTestBlock='';
    #number of lines to remove
    my $nRemove=0;
    #test description
    my $Name=''; my $Keys=''; my $Outs='';
    #error flag
    my $ErrorRead='';
    foreach my $line (@_){
	#remove spaces at the end of $line
	$line =~ s/\s*$//g;
	if($line =~ m/<#$/){
	    $IsTestBlock='0';
	    last;
	}
	if($IsTestBlock){
	    #extract test description: Name, Key, Outs
	    if($line =~ m/Name=(.*)/){
		unless($Name){$Name=$1;}
		else{
		    #there can be only ONE definition for NAME
		    $ErrorRead='1';
		}
	    }
	    elsif($line =~ m/Keys=(.*)/){$Keys="$Keys$1," if($1);}
	    elsif($line =~ m/Outs=(.*)/){$Outs="$Outs$1," if($1);}
	    else{$ErrorRead='1';}
	}
	if($line =~ m/^#>/){
	    $IsTestBlock='1';
	}
	$nRemove++;
    }
    #remove processed lines
    splice(@_,0,$nRemove+1);
    
    unless($ErrorRead){
	my @Keys = split(/,/,$Keys) if($Keys);$Keys='';
	my @Outs = split(/,/,$Outs) if($Outs);$Outs='';
	unless($ErrorRead){
	    #process parameters for this machine
	    $Name = process_option($Name);
	    #process Keys
	    foreach my $Key (@Keys){
		$Key = process_option($Key);
		$Keys="$Keys$Key " if($Key);
	    }
	    #process Outs
	    foreach my $Out (@Outs){
		$Out = process_option($Out);
		$Outs="$Outs$Out " if($Out);
	    }
    }
	(\@_,$Name,$Keys,$Outs);
    }
    else{
	(\@_,'','','');
    }
}

sub process_option{
    my $Machines;
    my $OptionName=$_[0];
    if($OptionName =~ m/(.*)\((.*)\)/){
	# check if current machine is in the list of
	$OptionName = $1;
	$Machines   = $2;
	# first, count number of occurences
	my $nOccur = () = $Machines =~ /$hostname/g;
	if($nOccur > 1) {return('');}
	if($nOccur==0){
	    #check if all machines are removed
	    if($Machines =~ m/^0/){return('');}
	}
	if($nOccur==1){
	    if($Machines =~ m/\+$hostname/){
		#do nothing
	    }
	    elsif($Machines =~ m/-$hostname/){return('');}
	    else{return('');}
	}
    }
    return($OptionName);
}
