#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
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
	push(@search, split(/,/,$path));
    }
}

$Help   = 1 if $#ARGV<0;                    # No source files

#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Source Code Manipulation}
#
#!ROUTINE: depend.pl - automatic generation of Fortran/C source dependencies
#
#!DESCRIPTION:
# Create a makefile with dependencies based on the 
# \begin{verbatim}
#   include 'file'
#   include "file"
#   #include "file"
#   #include <file>
#   use SomeModule
# \end{verbatim}
# statements. The script takes care of differencies in capitalization for Fortran,
# it associates files with the modules found in them, it also finds modules and
# header files in the search path. All in all this script figures out dependencies
# in an automated fashion for Fortran and C codes which can be used in Makefile-s
# or to get a dependency tree for its own sake.
#
#!REVISION HISTORY:
# 07/10/01 G.Toth - initial version written for BATSRUS.
# 08/20/03 G.Toth - added search path options with intelligent file association
# 03/20/04 G.Toth - added multiple -Ipath option so the compiler flags 
#                   can be used without any change
# 07/30/04 I.Sokolov - added search for include files in the search path
# 01/20/05 G.Toth - improved module -> object file association scheme
# 03/10/14 L.Daldorff and G.Toth - added C/C++ code processing
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

# Collect the modules and header files from the search path
my %env;        # Directory name --> '${SHELLVAR}'
my %modulefile; # Module object  --> File name containing the module
my $dir;
my %headeruse;  # filename --> header --> 1 if filename depends on header
my %headerdir;  # headerfilename --> directory it is in

# Loop over the search directories to find F90 modules
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

    my @source; # List of Fortran files
    @source = grep /\.(f|f90)$/i, readdir DIR;
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

&process_header_files;

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
    foreach $ext ('.f90','.F90','.f','.F','.cpp','.cxx','.cc','.c'){
	if(-e "$base$ext"){
	    $file = "$base$ext";
	    last SOURCE;
	}
    }
    if(not $file){
	print "$WARNING source file not found, skipping $object !!!\n"
	    unless $object =~ /^(main\.o|advect_main\.o|game_of_life\.o)$/;
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
	if(/^\s*include\s+[\"\']([^\'\"]+)/i and not /\bmpif.h\b/){
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
	if(/^\#include\s+[<\"]([^>\"]+)/i){
	    # Found C type header file. Find the directory it is from.
	    my $headerfile = $1;
	    my $dir = $headerdir{$1};

	    # Nothing to do if it is not a know directory
	    next unless $dir;

	    # Nothing to do if header file is already listed
	    next if $include{$base} =~ / $dir\/$headerfile\b/;

	    # Add the header file and header files that depend on this header file
	    my $headerfile2;
	    foreach $headerfile2 (sort ($headerfile, keys %{$headeruse{$headerfile}})){
		my $hfile = "$headerdir{$headerfile2}/$headerfile2";
		$include{$base} .= " $hfile" unless $include{$base} =~ / $hfile\b/;
	    }
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
  	      or /netcdf.o/i or /ezspline/i or /ezcdf.o/i or /^hdf5.o/i}
  	(@use);

	# Make string out of array
  	$use = ' '.join(' ',@use);
    }
    my $depend;
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

################################################################################
sub process_header_files{

    # Find header files
    my $dir;
    my $file;
    my $headerfile;

    # Loop over the object files and collect directories
    my $object;
    my %objectdir;	 # hash of directories containing object files
    foreach $object (@ARGV){
	$dir = ".";		       # assume local directory
	$dir = $` if $object =~ m#/[^/]+#; # check for / in the name
	$objectdir{$dir} = 1;
    }

    # Loop over object directories, find header files and store their dependencies
    foreach $dir (@search, keys %objectdir){
	next if $dir =~ /:/;
	opendir(DIR,$dir) or die "$ERROR: could not open object directory $dir\n";
	my @headerfiles = grep /\.h(xx)?$/i, readdir DIR;
	closedir DIR;
	foreach $file (@headerfiles){
	    # store headerfile into headerdir (assume that names are unique!)
	    $headerdir{$file} = $dir;

	    # See if the header file depends on other files
	    open(FILE, "$dir/$file") 
		or die "$ERROR: could not open header file $dir/$file\n";

	    # Find #include statements and store dependency
	    while (<FILE>){
		$headeruse{$file}{$1} = 1 if /^#include\s+[<\"]([^\">]+)/i;
	    }
	    close FILE;
	}
    }

    # Remove headers that are not in the search or object directories
    foreach $file (keys %headeruse){
	foreach $headerfile (keys %{$headeruse{$file}}){
	    delete $headeruse{$file}{$headerfile} unless $headerdir{$headerfile};
	}
    }

    # Add nested dependencies (loop as long as needed)
    my $headerfile2;
    LEVEL:{
	my $change = 0;
	foreach $file (keys %headeruse){
	    foreach $headerfile (keys %{$headeruse{$file}}){
		foreach $headerfile2 (keys %{$headeruse{$headerfile}}){
		    next if $headeruse{$file}{$headerfile2};
		    $headeruse{$file}{$headerfile2} = 1;
		    $change = 1;
		}
	    }
	}
	redo if $change;
    }

    #print "HEADER FILES:\n";
    #foreach $file (sort keys %headeruse){
    #    print "$file uses: ";
    #    foreach $headerfile (sort keys %{ $headeruse{$file} }){
    #	print " $headerdir{$headerfile}/$headerfile";
    #    }
    #    print "\n";
    #}

}
