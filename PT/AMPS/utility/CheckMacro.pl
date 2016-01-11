#!/usr/bin/perl
# -----------------------------------------------------------------------------
# This script goes through the given folder, finds all .cpp and .h files 
# and collects all used C++ preprocessor macros
# -----------------------------------------------------------------------------

use strict;
use warnings;

# names of folders to be skipped
my @FolderSkip = ('.', '..', 'CVS');

# extensions of files to be considered
my @FileType = ('h', 'cpp', 'c', 'dfn');

if( scalar(@ARGV) == 0){
    die 'There has to be at least one command line argument (see help)!';
} 
    

# check if user wants to add the check in-place
my $DoInPlace = 0;
if(defined $ARGV[1]){
    $DoInPlace = 1 if ($ARGV[1] =~ m/-in-place/i);
}

# Print help
# -----------------------------------------------------------------------------
if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help'){
    print "This script adds a check for all used " .
	"C++ preprocessor macros.\n" .
	"USAGE:\n".
	" ./CheckMacro.pl FOLDER\t\t\tTest is stored in a file MacroCheck.h\n".
	" ./CheckMacro.pl FOLDER -in-place\tTest is added in place\n";
    exit 0;
}

# Execute the script
# -----------------------------------------------------------------------------
# remember the current location
my $path = `pwd`;
chomp($path);

my $DirIn = $ARGV[0];
$DirIn =~ s/\/$//;

# container for macros
my @Macros;

# container for a file being processed
my @SourceLines;

# last found line with directive #if
my $iLastIf = -1;

# collect macros
&process_directory($DirIn);

if(!$DoInPlace){
    # put a check for each macro in a header file
    my $HeaderOut = 'MacroCheck.h';
    my @HeaderLines;
    push(@HeaderLines,"#ifndef _PIC_MACRO_CHECK_\n".
	 "#define _PIC_MACRO_CHECK_\n#endif\n");
    foreach my $Macro (@Macros){
	push(@HeaderLines, 
	     "#ifndef $Macro\n".
	     "#error ERROR: $Macro is used but not defined\n".
	     "#endif\n");
    }

    open(FILE,'>', $HeaderOut);
    print FILE @HeaderLines;
    close(FILE);
}
# =============================================================================

sub process_directory{
    # check correctness
    die("Trying to process file $path/$_[0] as a folder!") 
	unless (-d "$path/$_[0]");

    # move to the given directory
    opendir(my $DIR, "$path/$_[0]") or die $!;
    $path = "$path/$_[0]";

    # collect all files
    while(my $file = readdir($DIR)){
	# check if this is a folder
	if(-d "$path/$file"){
	    # check if this folder is in the list of exceptions
	    next if (grep {$_ eq $file} @FolderSkip);
	    # process this directory as well
	    &process_directory($file);
	}
	else {
	    # ignore temporary files
	    next if($file =~ m/^.*~/ || $file =~ m/^#.*#$/);
	    # ignore any file that is not a C++ code
	    next unless (grep {$file =~ m/.*\.$_$/} @FileType);
	    &process_file("$path\/$file");
	}
    }
    # move out of the given directory
    closedir($DIR);
    $path =~ s/(.*)\/(.*)?/$1/;
}

sub process_file{
    # check correctness
    die("Trying to process folder $_[0] as a file!") if (-d "$_[0]");

    # process this file

    # read content
    open(my $FILE, '<', $_[0]);
    my @lines = <$FILE>;
    @SourceLines = @lines;
    close($FILE);

    # find all macros in this file
    # -----------------------------
    # flag to mark a comment section of form /* */
    my $IsComment   = 0;
    # go through the file
    foreach (my $i=0; $i<=@lines-1; $i++){
	# remove a comment part of the line
	$lines[$i] =~ s/\/\/.*//;
	$lines[$i] =~ s/\/\*.*?\*\///g;
	if($lines[$i] =~ m/.*\*\// && $IsComment){
	    $lines[$i] =~ s/.*\*\///; $IsComment = 0;}
	next if($IsComment);
	if($lines[$i] =~ m/\/\*.*/){
	    $lines[$i] =~ s/\/\*.*//; $IsComment = 1;}

	# find if there is a macro being used in this line
	if($lines[$i] =~ m/\s*#(el)?if (.*)?==(.*)/){
	    # check if it is #if
	    $iLastIf = $i if (! defined $1);
	    # extract macros from matched blocks
	    &extract_macro($2);
	    &extract_macro($3);
	}
    }
    # add a check within the file itself
    if($DoInPlace){
	open(my $FILE, '>', $_[0]);
	print $FILE @SourceLines;
	close($FILE);
    }

}


sub extract_macro{
    # extract EXACTLY ONE macro from the input string
    my $Macro = $_[0];
    # check of string consists of whitespace and digits only => skip
    if($Macro =~ m/^[\s0-9]*$/){return;}
    # check if this is a macro-function
    if($Macro =~ m/(.*)\((.*)\)/){
	&extract_macro($1);
	return;
    }
    # trim: remove leading and trailing whitespaces
    $Macro =~ s/^\s+|\s+$//g;
    
    if($DoInPlace){
	# add the check within the file itself in front of last found #if
	$SourceLines[$iLastIf] = 
	    "#ifndef $Macro\n".
	    "#error ERROR: $Macro is used but not defined\n".
	    "#endif\n".
	    $SourceLines[$iLastIf];
    }
    else{
	# add to list of macros (unless is present already)
	push(@Macros,$Macro) unless(grep {$Macro eq $_}@Macros);
    }
    
}


