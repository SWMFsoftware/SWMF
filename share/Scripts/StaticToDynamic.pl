#!/usr/bin/perl -s

my $Help  = ($h or $H or $help);
my $Debug = $D;
my $types = ($t or 'real|logical|integer');     # types of the variable
my $large = ($l or 'nBLK|MaxBlock|MaxImplBLK'); # pattern for large variables
my $mindim= ($d or 3);                          # minimum dimension 

use strict;

&print_help if $Help;

my $cont;            # true if previous line is continued
my $line;            # current line with continuation lines
my %dimension;       # hash contains dimension for variables
my $type;            # the type of the variables in the current declaration
my @allocated;       # list of allocated variables
my $begin;           # true if a configuration is started
my $config;          # current configuration directive
my $configend;       # the config to be written at the end of declaration
my %config;          # hash for configuration of allocate and deallocate
my $allocate;        # true if allocation was written
my $deallocate;      # true if deallocation was written
my $contains;        # true if 'contains' statement has been found

my $WARNING = "StaticToDynamic.pl WARNING:";
my $ERROR = "StaticToDynamic.pl ERROR:";

while(<>){

    # Set $contains variable
    $contains = 1 if /^\s*contains\b/;

    if($contains){
	# Add allocate statements into subroutine allocate_*
	if(/^\s*end\s+subroutine\s+allocate_/i){
	    $allocate=1;
	    my $variable;
	    foreach $variable (@allocated){
		print "    allocate($variable".
		    "$dimension{$variable},stat=iError) $config{$variable}\n";
	    }
	    print "\n";
	}

	# Add deallocate statements into subroutine deallocate_*
	if($contains and /^\s*end\s+subroutine\s+deallocate_/i){
	    $deallocate=1;
	    my $variable;
	    foreach $variable (@allocated){
		print "    deallocate($variable) $config{$variable}\n";
	    }
	    print "\n";
	}
	print; next;
    }

    # process lines which are before the contains statement

    # Initialize configuration for this line
    $configend = '';
    $config = '' unless $begin;

    # store configuration
    if(/(\^CFG\s+IF\s+(NOT\s+)?\w+)(\s+BEGIN)?/i){
	$config = "!$1";
	$begin  = $3;
	$configend = $config unless $begin;
	print "config=$config begin=$begin configend=$configend\n" if $Debug;
    }

    # process end of configuration
    if(/\^CFG\s+END\s/){
	warn "$WARNING \^CFG END without \^CFG ... BEGIN!?\n" unless $begin;
	$config= '';
	$begin = '';
	print "config=$config begin=$begin configend=$configend\n" if $Debug;
    }

    # Check line for variable declaration
    if(not $cont){
	# Ignore lines that do not look like a declaration
	if(not /^(\s*($types|common))\b/i){
	    print; next;
	}
	# Store the type of the variables
	$type = $1;
    }

    # Collect continuation lines together
    if($cont){
	# Append continuation line
	$line .= $_;
    }else{
	# Just set $line
	$line = $_;
    }

    # Check if the current line is a continuation line
    if($line =~ /\&\s*(\!.*)?$/){
	$cont = 1; next;
    }else{
	$cont = 0;
    }

    $_=$line;

    # Remove all common statements
    next if $type =~ /common/i;

    # Do nothing if no large arrays are declared
    if(not /$large/i){
	print; next;
    }

    my @variables;

    # Check for "dimension(size1,size2,size3) :: "
    if(/\bdimension\s*(\([^\(]+\))\s*(::)?\s*/i){
	my $dimension = $1;
	my $variables = $';

	# Check if number of dimensions is larger than $mindim
	my $ndim = ($dimension =~ s/,/,/g) + 1;
	if($ndim < $mindim){
	    print; next;
	}

	# Remove comments and junk from variable list
	$variables =~ s/\!.*//ig;
	$variables =~ s/[\s\&\n]//g;
	@variables = split(/,/,$variables);

	# Set the dimensionality of the variables
	my $variable;
	foreach $variable (@variables){
	    $dimension{$variable} = $dimension;
	}
    }else{
	# Match the beginning of the variable declaration
	/^\s*$type\s*(::)?\s*/;

	# Add a comma to the end of the variable list for easier parsing
	my $variables = $'.',';

	# Remove comments and junk from variable list
	$variables =~ s/\!.*//ig;
	$variables =~ s/[\s\&\n]//g;

	# Extract variables one by one from the list
	while($variables =~ s/^\s*(\w+)\s*(\([^\)]+\))?\s*,//){
	    my $variable = $1;
	    my $dimension= $2;
	    push(@variables,$variable);
	    $dimension{$variable} = $dimension;
	}
    }

    if($Debug){
	print "VARIABLES=";
	my $variable;
	foreach $variable (@variables){
	    print $variable,'[',$dimension{$variable},']';
	}
	print "\n";
    }

    # Write out allocatable or normal declaration
    my $variable;
    foreach $variable (@variables){
	my $dimension = $dimension{$variable};
	if($dimension =~ /$large/){
	    my $ndim = ($dimension =~ s/[^\(\),]+/:/g);
	    if($ndim >= $mindim){
		print "$type, allocatable :: $variable$dimension $configend\n";
		# store variable and its configuration
		push(@allocated,$variable);
		$config{$variable}=$config;
		next;
	    }
	}
	print "$type :: $variable","$dimension{$variable} $configend\n";
    }
}

warn "$ERROR subroutine allocate_* was not found in $ARGV !!!\n" 
    unless $allocate;
warn "$WARNING subroutine deallocate_* was not found in $ARGV !!!\n" 
    if $Debug and not $deallocate;

exit 0;

##############################################################################

sub print_help{

    print "

Purpose: Transform a module with static variable declarations into a 
         module with dynamic (allocatable) declarations.
         The original module must contain a subroutine named 'allocate_*'
         and it may contain a subroutine named 'deallocate_*'.
         These subroutines can be empty in the static version of the module.

Usage:   StaticToDynamic.pl [-h] [-D] [-t=TYPES] [-d=MINDIM] [-l=LARGE] 
               [INFILE] [> OUTFILE]

   -h         - this help message

   -D         - print debugging information

   -t=TYPES   - TYPES contains a | separated list of the types of variables 
                to be transformed. Default value is -t='real|logical|integer'

   -d=MINDIM  - MINDIM is the minimum dimensionality of variables to be 
                transformed. Default values is -d=3.

   -l=LARGE   - LARGE is a pattern matching the large variable.
                Default values is -l='nBLK|MaxBlock|MaxImplBLK'.

   INFILE     - the file to be transformed. Default is reading from STDIN.

   OUTFILE    - the output file. By default output goes to STDOUT.

Examples:

  Convert ModAdvance_static.f90 to a dynamic ModAdvance.f90:

StaticToDynamic.pl ModAdvance_static.f90 > ModAdvance.f90

  Convert only the real and double precision arrays of at least 4 dimensions 
  matching 'nBLK' or 'MaxBlock' in their declaration:

StaticToDynamic.pl -t='real|double precision' -d=4 -l='nBLK|MaxBlock'

";
    exit;
}
