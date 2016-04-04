#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'PT';
our $Code            = 'AMPS';
our $MakefileDefOrig = 'Makefile.def.amps';
our @Arguments       = @ARGV;

my $Application;
my $MpiLocation;
my $Mpirun="mpirun";
my $TestRunProcessorNumber=4;

# name for the nightly tests
our $TestName;
our @Compilers;

#create Makefile.local
if (! -e "Makefile.local") {
  `touch Makefile.local`;
}

# build AMPS' Makefile.test
foreach (@Arguments) { 
    if(/^-install/) {
      require "utility/TestScripts/BuildTest.pl";
     
      #create .general.conf if not exists
      `touch .general.conf` unless (-e ".general.conf");
    }
}


# Run the shared Config.pl script
my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

our %Remaining;   # Arguments not handled by share/Scripts/Config.pl

#sort the argument line to push the 'application' argument first
my $cntArgument=0;

foreach (@Arguments) {
  if (/^-application=(.*)/i) {
    if ($cntArgument!=0) {
      my $t;
      
      $t=$Arguments[0];
      $Arguments[0]=$Arguments[$cntArgument];
      $Arguments[$cntArgument]=$t;
      
      last;
    }
  }
  
  $cntArgument++;
}

#process the 'help' argument of the script
foreach (@Arguments) { 
   if ((/^-h$/)||(/^-help/)) { #print help message
     print "Config.pl can be used for installing and setting of AMPS\n";
     print "AMPS' settings can be defined in the SWMF's Config.pl: Config.pl -o=PT:spice-path=path,spice-kernels=path,ices-path=path,application=casename\n\n";
     print "Usage: Config.pl [-help] [-show] [-spice-path] [-spice-kernels] [-ices-path] [-cplr-data-path] [-application]\n";
     print "\nInformation:\n";
     print "-h -help\t\t\tshow help message.\n";
     print "-s -show\t\t\tshow current settings.\n";
     print "-append\t\t\t\tadd new settings variables into .amsp.conf instead of re-defining their values.\n";
     print "-spice-path=PATH\t\tpath to the location of the SPICE installation.\n";
     print "-spice-kernels=PATH\t\tpath to the location of the SPICE kernels.\n";
     print "-ices-path=PATH\t\t\tpath to the location of ICES\n";
     print "-application=case\t\tthe name of the model case to use.\n";
     print "-boost-path=PATH\t\tthe path to the boost library.\n";
     print "-kameleon-path=PATH\t\tthe path for the Kameleon library.\n";
     print "-cplr-data-path=PATH\t\tthe path to the data files used to import other model results by PIC::CPLR::DATAFILE.\n";
     print "-model-data-path=PATH\t\tthe path to the user-defined input data files (the files are read directly by the user).\n";
     print "-output-path=PATH\t\tthe directory where AMPS' output files will be saved.\n";
     print "-batl-path=PATH\t\t\tthe path to the BATL directory\n";
     print "-swmf-path=PATH\t\t\tthe path to the SWMF directory\n";
     print "-set-test(=NAME)\/comp\t\tinstall nightly tests (e.g. comp=gnu,intel|pgi|all)\n";
     print "-rm-test\/comp\t\t\tremove nightly tests\n";
     print "-amps-test=[on,off]\t\ttells the code that a nightly test is executed\n";
     print "-openmp=[on,off]\t\twhen \"on\" use OpenMP and MPI libraries for compiling AMPS\n";
     print "-link-option=-opt1,-opt2\tadd options \"-opt1 -opt2\" to linker\n";
     exit;
   }
   
   
   if ((/^-s$/i)||(/^-show/)) { #print AMPS settings
     print "\nAMPS' Setting options:\n";
      
     if (-e ".amps.conf") {
       my @f;
       
       open (FSETTINGS,"<.amps.conf");
       @f=<FSETTINGS>;
       close(FSETTINGS);
       
       print @f;
       print "\n";
     }
     else {
       print "No custom settings\n";
     }
     
     exit;
   } 
}

#process the configuration settings
foreach (@Arguments) {
  if (/^-application=(.*)/i) {
    my $application = lc($1);
  
    #default values
    add_line_amps_conf("TESTMODE=off");
    add_line_amps_conf("BATL=nobatl");
    add_line_amps_conf("KAMELEON=nokameleon");
    
    `cp -f input/$application.* input/species.input .`;

     #remove path from the name of the input file is such exists
     my @list;

     $application=~s/\// /g;
     @list=split(' ',$application);
     $application=$list[$#list];


    add_line_amps_conf("InputFileAMPS=$application.input");   
    add_line_amps_conf("APPLICATION=$application");

    `echo "InputFileAMPS=$application.input" >> Makefile.local`;
    next
  };
    
  if (/^-mpi=(.*)$/i)        {$MpiLocation=$1;                next}; 
  if (/^-np=(.*)$/i)         {$TestRunProcessorNumber=$1;     next};
  if (/^-spice-path=(.*)$/i)      {
      add_line_amps_conf("SPICE=$1");
      `echo SPICE=$1 >> Makefile.local`;
      next}; 
  if (/^-spice-kernels=(.*)$/i)    {
      add_line_amps_conf("SPICEKERNELS=$1");
      next};
  if (/^-ices-path=(.*)$/i)       {
      add_line_amps_conf("ICESLOCATION=$1");
      next};
  
  if (/^-boost-path=(.*)$/i)       {
      add_line_amps_conf("BOOST=$1");
      next};
  if (/^-kameleon-path=(.*)$/i)       {
      add_line_amps_conf("KAMELEON=$1");
      next};
      
      
      
      
  if (/^-openmp=(.*)$/i) {
    my $t;
    $t=lc($1);
    add_line_amps_conf("OPENMP=$1");
       
    if ($t eq "on") {
      add_line_general_conf("#undef _COMPILATION_MODE_ \n#define _COMPILATION_MODE_ _COMPILATION_MODE__HYBRID_\n");
      `echo OPENMP=on >> Makefile.local`;
      next;
    }
    elsif ($t eq "off") {
      add_line_general_conf("#undef _COMPILATION_MODE_ \n#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_\n");
      `echo OPENMP=off >> Makefile.local`;
      next;
    }
    else {
      die "Option is unrecognized: -openmpfgf=($1)";
    }
  }   
      
  
  if (/^-batl-path=(.*)$/i)       {
      add_line_amps_conf("BATL=$1");
      next};
  if (/^-swmf-path=(.*)$/i)       {
      add_line_amps_conf("SWMF=$1");
      next};  
  
  if (/^-cplr-data-path=(.*)$/i)       {
      add_line_amps_conf("CPLRDATA=$1");
      next};  #path to the data files used in the PIC::CPLR::DATAFILE file readers
  if (/^-model-data-path=(.*)$/i)  {
      add_line_amps_conf("MODELINPUTDATA=$1");
      next}; #the path to the data file used by the user direactly (not through the PIC::CPLR::DATAFILE readers)

  if (/^-output-path=(.*)$/i)       {
      add_line_amps_conf("OUTPUT=$1");
      next};  #the directory where AMPS' output files will be located 

  #compile for the nightly test
  if (/^-amps-test=(.*)$/i)  {
    my $t;

    $t=lc($1);
    add_line_amps_conf("TESTMODE=$t");
    `echo TESTMODE=$t >> Makefile.local`;
    next;
  }; 

  #set the APPEND flog to append .amps.conf instead of re-defining of AMPS' settings
  if (/^-append$/i) { 
    `echo APPEND >> .amps.conf`;
    next;
  }

  # set nightly test:
  #  -set-name=NAME/comp - tests' name NAME cannot be empty
  #  -set-name/comp      - tests' name will be `hostname -s`
  if (/^-set-test=(.*)$/i)      {
      die "ERROR: test name is empty!" unless ($1);
      my $CompilersRaw;
      ($TestName, $CompilersRaw) = split("\/",$1,2);
      die "ERROR: no compiler is indicated for tests!" unless ($CompilersRaw);
      @Compilers = split (',',lc $CompilersRaw);
      require "./utility/TestScripts/SetNightlyTest.pl";
      next}; 
  if (/^-set-test\/(.*)$/i)      {
      die "ERROR: no compiler is indicated for tests!" unless ($1);
      $TestName=''; 
      @Compilers = split (',',lc $1);
      require "./utility/TestScripts/SetNightlyTest.pl";     
      next}; 

  if (/^-rm-test$/i)      {require "./utility/TestScripts/RemoveNightlyTest.pl";     next}; 
  if(/^-link-option=(.*)$/){
      my $options=$1; $options =~ s/,/ /g;
      open(my $fh,'<',"Makefile"); my @lines = <$fh>; close($fh);
      open(my $fh,'>',"Makefile");
      foreach my $line (@lines){
	  if($line =~ m/^EXTRALINKEROPTIONS=/){
	      $line = "EXTRALINKEROPTIONS=$options\n";
	  }
	  print $fh $line;
      }
      close($fh);
      next;
  };
  
  warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}


exit 0;

#=============================== Add a line to .general.conf
sub add_line_general_conf {
  #create .general.conf if not exists
  `touch .general.conf` unless (-e ".general.conf");
  
  open (SETTINGS,">>",".general.conf") || die "Cannot open .general.conf\n";
  print SETTINGS "$_[0]\n";
  close SETTINGS;
}

#=============================== Add a line to .amps.conf
# USAGE:
# add_line_amps_conf($Newline)
#  $Newline has a form 'PARAMETER=value'
# if PARAMETER has been defined before, it is OVERWRITTEN with the new value
sub add_line_amps_conf{
    # create Makefile.local if it doesn't exist
    `touch .amps.conf` unless (-e ".amps.conf");
    
    # check that there are exactly 2 arguments
    my $nArg = @_;
    die "ERROR: add_line_amps_conf takes exactly 1 argument" 
	unless ($nArg==1);

    # get the name of parameter being defined
    my $Parameter = '';
    if($_[0] =~ m/(.*?)=(.*)/){
	$Parameter = $1;
    }
    else{die "Trying to add invalid line to .amps.conf: $_[0]\n";}

    # trim the parameter name as well as newline symbol at the end
    $Parameter =~ s/^\s+//; $Parameter =~ s/\s+\n*$//;

    # read the current content
    open (SETTINGS,"<",".amps.conf") || 
	die "Cannot open .amps.conf\n";
    my @settingslines=<SETTINGS>;
    close(SETTINGS);

    #check if the APPEND keywork is present in the settings
    my $AppendFlag=0;

    foreach (@settingslines) {
      if ($_ =~ /APPEND/){
        $AppendFlag=1;
        last;
      }
    }

    # check if this parameter has already been defined
    my $IsPresent = 0;

    if ($AppendFlag == 0) {
      foreach (@settingslines) {
        if ($_ =~ /$Parameter\s*=.*/){
          $_ = $_[0]; chomp($_); $_ = "$_\n";
          $IsPresent = 1;
          last;
        }
      }
    }
    
    # if the parameter hasn't been defined before, add its definition
    unless ($IsPresent) {
	# temporary variable for storing a line to be printed
	my $tmp; $tmp = $_[0]; chomp($tmp);# $tmp = "$tmp\n";
	push(@settingslines,"$tmp\n") 
    }

    
    # write changed content of .amps.conf
    open (SETTINGS,">",".amps.conf");   
    print SETTINGS  @settingslines;
    close (SETTINGS);
}
