#!/usr/bin/perl
#^CFG COPYRIGHT UM

# Default values
my $Output = "Makefile.DEPEND"; # Default output
my $Help;                       # No help unless needed
my @search;                     # Array of search

# Translation between directory names and Makefile definitions
my %defdir = ('/share/Library/src' => SHAREDIR,
	      '/CON/Library/src'   => LIBRARYDIR,
	      '/CON/Coupler/src'   => COUPLERDIR,
	      '/CON/Interface/src' => INTERFACEDIR);

# Read flags
while($ARGV[0] =~ /-/){
    my $flag = shift(@ARGV);
    if($flag =~ /^-o=/){$Output=$'}; # -o=Makefile.test
    if($flag =~ /^-h/i){$Help=1};    # -h -help -H -Help
    if($flag =~ /^(-s|-I)=?/){        # -s=path1,path2 -Ipath
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

#!QUOTE: \section{share/Scripts: for SWMF and Physics Modules}
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

# Collect the modules from the search path
foreach $dir (@search){
   
    if($dir =~ /:/){
	# Split environment name from dir name if a colon is present
	($env,$dir)=split(':',$dir);

	# Store environment variable for this directory in suitable form
        $env{$dir}='${'.$env.'}';
    }

    -d $dir or die "depend.pl ERROR: $dir is not a directory\n";
    opendir(DIR,$dir) or die "depend.pl ERROR: ".
	"could not open directory $dir\n";

    @source = grep /\.f90$/i, readdir DIR;
    closedir DIR;

    foreach $file (@source){
	open FILE,"$dir/$file" or die "depend.pl ERROR: ".
	    "could not open $dir/$file\n";
	while(<FILE>){

	    if(/^\s*module\s+(\w+)/i){
		$module = uc($1); # capitalize module name (ignore case)
		$object = $module.'.O'; 

		
                # Check if the file with the same name as the module
		# has been found already

		# Remove IH_ or SC_ from the name of the module for matching
		$modmatch = $module; $modmatch =~ s/^(IH|SC)_//;

		if($modulefile{$object} !~ /\/$modmatch\.o$/i){
		    # If not, store the filename into %modulefile
		    $file =~ s/\.f90$/.o/i;
		    if($env{$dir}){
			# Store the name using the environment variable
			$modulefile{$object}="$env{$dir}/$file";
		    }else{
			# Store the full path
			$modulefile{$object}="$dir/$file";
		    }
		}
	    }
	}
	close FILE;
    }
}

open(OUTPUT,">$Output") or die "depend.pl ERROR: ".
    "could not open dependency file $Output !!!\n";

OBJECT: 
    while($object=shift(@ARGV)){

    $base=$object;
    # Skip files in other directories
    next OBJECT if $base=~/^\.\./;

    # Skip files which do not have the .o extension
    $base=~s/\.o$// or next OBJECT;

    $file = '';

  SOURCE: 
    foreach $ext ('.f90','.F90','.f','.F'){
	if(-e "$base$ext"){
	    $file = "$base$ext";
	    last SOURCE;
	}
    }
    if(not $file){
	print "ERROR: source file not found, skipping $object !!!\n";
	next OBJECT;
    }

    if(not open(FILE,$file)){
	print "Error opening file $file !!!\n";
	next OBJECT;
    }

    # Store list of names without extension
    push(@base,$base);

    # Build dependency list $depend for file $file
    $depend="";
    while($_=<FILE>){
	# Collect module file names corresponding to the modules
	if(/^\s*module\s+(\w+)/i){
	    $module=uc("$1\.o"); # capitalize module name (ignore case)
	    $modulefile{$module}=$object;
	}
	# Check for 'use module'
	if(/^\s*use\s+(\w+)/i){
	    $module="$1.o";
	    $use{$base}.=" $module" 
		unless $use{$base}=~/ $module\b/i;
	    ### print "base=$base uses $use{$base}\n";
	}
	# Check for 'include "filename"'
	if(/^\s*include\s+[\"\']?([^\'\"\s]+)/i){
	    $include=$1;
	    $includeorig=$include;
	    if(not -e $include){
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
	    $include{$base}.=" $include" 
		unless $include{$base}=~/ $include\b/;
	}
    }
}

foreach $base (@base){
    $use=$use{$base};
    if($use){
	# Correct module names to file names
	@use=split(' ',$use);
	map {if($mfile=$modulefile{uc($_)}){$_=$mfile}} (@use);

        # Exclude dependency on itself and compiler provided modules
        map { $_='' if $_ eq "$base.o" or /F90_UNIX_IO/i}(@use);
	    
        # Make string out of array
	$use=' '.join(' ',@use);
    }
    $depend=$include{$base}.$use;

    # Get rid of leading and trailing spaces
    $depend =~ s/^ *//; $depend =~ s/ *$//;

    # Replace space separator with continuation lines and tabs
    $depend =~ s/ +/ \\\n\t/g;

    # Write dependency rule into output file
    print OUTPUT "$base\.o:\\\n\t$depend\n\n" if $depend;

}
close OUTPUT;

exit;
