#!/usr/bin/perl
#^CFG COPYRIGHT UM
use strict;

# Default values
my $Output = "Makefile.DEPEND"; # Default output
my $Help;                       # No help unless needed
my @search;                     # Array of search

# Error and warning messages
my $ERROR   = "ERROR in depend.pl:";
my $WARNING = "WARNING in depend.pl:";

# Translation between directory names and Makefile definitions
my %defdir = ('/share/Library/src' => 'SHAREDIR',
	      '/CON/Library/src'   => 'LIBRARYDIR',
	      '/CON/Coupler/src'   => 'COUPLERDIR',
	      '/CON/Interface/src' => 'INTERFACEDIR');

# Read flags
while($ARGV[0] =~ /-/){
    my $flag = shift(@ARGV);
    if($flag =~ /^-o=/){$Output=$'};  # -o=Makefile.test
    if($flag =~ /^-h/i){$Help=1};     # -h -help -H -Help
    if($flag =~ /^(-s|-p|-I)=?/){        # -s=path1,path2 -Ipath
	my $path = $';  $path = shift(@ARGV) unless $path;

	if($path !~ /[:,]/){
	    # Add 'variable:' for directories defined in %defdir hash
	    my $dir;
	    foreach $dir (keys %defdir){
		$path = "$defdir{$dir}:$path" if $path =~ /$dir/;
	    }
	}

	# Store the path
	push(@search,split(/,/,$path));
    }
}

$Help   = 1 if $#ARGV<0;                    # No source files

#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Source Code Manipulation}
#
#!ROUTINE: depend.pl - automatic generation of Fortran source dependencies
#
#!DESCRIPTION:
# Create a makefile with dependencies based on the 
# \begin{verbatim}
#   include 'file'
#   include "file"
#   use SomeModule
# \end{verbatim}
# statements. The script takes care of differencies in capitalization,
# it associates files with the modules found in them, it also finds
# modules in the search path. All in all this script figures out dependencies
# in an automated fashion for Fortran codes which can be used in Makefile-s
# or to get a dependency tree for its own sake.
#
#!REVISION HISTORY:
# 07/10/01 G.Toth - initial version written for BATSRUS.
# 08/20/03 G.Toth - added search path options with intelligent file association
# 03/20/04 G.Toth - added multiple -Ipath option so the compiler flags 
#                   can be used without any change
# 07/30/04 I.Sokolov - added search for include files in the search path
# 01/20/05 G.Toth - improved module -> object file association scheme
#EOP
if($Help){
    print 
#BOC
'Usage: depend.pl [-h] [-o=filename] [-s=path] [-Ipath] file1 file2

Options:

-h             this help message

-Ipath -s=path look for modules in the comma separated list of directories.
               The directory name can be preceeded with an environment
               variable name and a colon. This flag can be given multiple
               times. The -s= format is kept for backwards compatibility.

-o=filename    write dependencies into filename (default is Makefile.DEPEND)

Examples of use:

In a makefile:  depend.pl -s="SEARCHDIR:${SEARCHDIR},../src" ${OBJECTS}

                SEARCH_EXTRA = -I${LIBRARYDIR} -I../src
                depend.pl ${SEARCH} ${SEARCH_EXTRA} ${OBJECTS}

Interactively:  depend.pl -o=Makefile.test main.o ModMain.o'
#EOC
    ,"\n\n";
    exit 1;
}

# Open output file now so the code dies fast if there is an error
open(OUTPUT,">$Output") or 
    die "$ERROR could not open dependency file $Output !!!\n";

# Collect the modules from the search path
my %env;        # Directory name --> '${SHELLVAR}'
my %modulefile; # Module object  --> File name containing the module
my $dir;
foreach $dir (@search){

    if($dir =~ /:/){
	# Split environment name from dir name if a colon is present
	my $env;
	($env,$dir)=split(':',$dir);

	# Store environment variable for this directory in suitable form
        $env{$dir}='${'.$env.'}';
    }

    -d $dir or die "$ERROR $dir is not a directory\n";
    opendir(DIR,$dir) or die "$ERROR: could not open directory $dir\n";

    my @source; # List of fortran 90 files
    @source = grep /\.(f|f90|h)$/i, readdir DIR;
    closedir DIR;

    my $file; # Actual F90 file
    foreach $file (@source){
	open FILE,"$dir/$file" or die "$ERROR: could not open $dir/$file\n";

	# Form object name from source file name
	my $objectfile = $file; $objectfile =~ s/\.f90$/.o/i;

	while(<FILE>){

	    if(/^\s*module\s+(\w+)/i){
		my $module = uc($1); # capitalize module name (ignore case)
		my $object = $module.'.O'; # capitalized object file name

		# The object file must exist already.
		# If there are multiple source+object files for the same module
		# use the source file with the name matching the module name.

		# Remove IH_ or SC_ from the name of the module for matching
		# with the file name (the file name is never renamed)
		$module =~ s/^(IH|SC)_//;

		if(-e "$dir/$objectfile" and 
		   $modulefile{$object} !~ /\/$module\.o$/i){
		    # If not, store the filename into %modulefile
		    if($env{$dir}){
			# Store the name using the environment variable
			$modulefile{$object}="$env{$dir}/$objectfile";
		    }else{
			# Store the full path
			$modulefile{$object}="$dir/$objectfile";
		    }
		}
	    }
	}
	close FILE;
    }
}

my @base;    # List of base names (without extension)
my %use;     # Base name --> space separated list of used module objects
my %include; # Base name --> space separated list of include files

my $object;  # Name of object file
OBJECT: 
    while($object=shift(@ARGV)){

    my $base=$object;
    # Skip files in other directories
    next OBJECT if $base=~/^\.\./;

    # Skip files which do not have the .o extension
    $base=~s/\.o$// or next OBJECT;

    my $file; # Name of the source file corresponding to the object file

    # Try different extensions
    my $ext;
  SOURCE: 
    foreach $ext ('.f90','.F90','.f','.F'){
	if(-e "$base$ext"){
	    $file = "$base$ext";
	    last SOURCE;
	}
    }
    if(not $file){
	print "$WARNING source file not found, skipping $object !!!\n";
	next OBJECT;
    }

    if(not open(FILE, $file)){
	print "$WARNING error opening file $file !!!\n";
	next OBJECT;
    }

    # Store list of names without extension
    push(@base, $base);

    # Build dependency list $depend for file $file
    my $depend;
    while($_ = <FILE>){
	# Collect module file names corresponding to the modules
	if(/^\s*module\s+(\w+)/i){
	    my $module=uc("$1\.o"); # capitalize module name (ignore case)
	    $modulefile{$module}=$object;
	}
	# Check for 'use module'
	if(/^\s*use\s+(\w+)/i){
	    my $module="$1.o";

	    # Append module object to the %use hash if it is not yet listed
	    $use{$base}.=" $module" 
		unless $use{$base}=~/ $module\b/i;
	}
	# Check for 'include "filename"'
	if(/^\s*include\s+[\"\']([^\'\"]+)/i){
	    my $include=$1;
	    my $includeorig=$include;
	    # If include file is not found check the search path
	    if(not -e $include){
		my $dir;
		foreach $dir (@search){
		    if( -e "$dir/$include"){
			if($env{$dir}){
			    $include="$env{$dir}/$include";
			}else{
			    $include="$dir/$include";
			}
			last;
		    }
		}
	    }
	    # Append include file to the %include hash if it is not yet listed
	    $include{$base} .= " $include" 
		unless $include{$base} =~ / $include\b/;
	}
    }
}

my $base; # Name of base file
foreach $base (@base){
    my $use; # Space separeted list of used module objects
    $use = $use{$base};
    if($use){
	# Correct module names to file names
	my @use = split(' ',$use);
	my $mfile;
	map {if($mfile=$modulefile{uc($_)}){$_=$mfile}} (@use);

        # Exclude dependency on itself and compiler provided modules
        map { $_='' if $_ eq "$base.o" or /F90_UNIX_IO/i or /ESMF_Mod.o/i
	      or /netcdf.o/i or /ezspline/i or /ezcdf.o/i}
	(@use);
	    
        # Make string out of array
	$use=' '.join(' ',@use);
    }
    my $depend; # Space separated list of include files and used module objects
    $depend = $include{$base}.$use;

    # Get rid of leading and trailing spaces
    $depend =~ s/^ *//; $depend =~ s/ *$//;

    # Replace space separator with continuation lines and tabs
    $depend =~ s/ +/ \\\n\t/g;

    # Write dependency rule into output file
    print OUTPUT "$base\.o:\\\n\t$depend\n\n" if $depend;

}
close OUTPUT;

exit;
