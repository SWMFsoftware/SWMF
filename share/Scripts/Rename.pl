#!/usr/bin/perl -s
#^CFG COPYRIGHT UM
#BOP
#!ROUTINE: Rename.pl - rename many variables/strings in many source files
#!INTERFACE
#!DESCRIPTION:
# This script makes variable renaming relatively safe and easy in 
# a large number of source files. It can rename many variables at
# the same time. Case insensitivity is taken care of, sub strings
# are not replaced, and conflicting renaming rules are checked for.
# The renaming rules are stored in an associative array, which
# is an easy to edit Perl code. 
#
#!REVISION HISTORY:
# 07/13/2001 G.Toth gtoth@umich.edu - initial version
# 07/23/2004 G.Toth allow protecting some lines
#EOP

$Help=$h; 
$Input=$i;
$List=$l; 
$Check=$c;  
$Rename=$r; 
$Undo=$u; 
$Debug=$d;
$Quiet=$q;
$Warning=$w;

if($Check and $Warning){
    print "Either check (-c) or no warnings (-w) \n\n";
    $Help=1;
}

if($List + $Check + $Rename + $Undo != 1){
    print "Specify exactly one of the option out of  -l, -c, -r, -u !\n\n";
    $Help=1;
}

$Input="RenameList.pl" unless $Input;

if($Help or not ($List or $Check or $Rename or $Undo) ){
    print 
#BOC
'Usage: Rename.pl [options...] file1 file2 ...

Options (specify them separately as  -a -b  and not as  -ab !):

-h             help (this message)

-i=inputfile   input file (default is RenameList.pl)

-l             list names to be changed
-c             check input file and source files, do not rename variables
-r             replace old names with new
-u             undo replacements

-d             debug info is printed
-q             quiet run, errors and warnings printed only
-w             warning messages are suppressed

You have to specify exactly one of -l, -c, -r, or -u. 

For replace (-r) and undo (-u) the original file is put into filename~
if any replacements were done. Individual lines can be protected against
renaming by adding a "!do not rename" trailing comment. 

Typical usage:

Rename.pl -c              #  check the replacement rules
Rename.pl -c *.f90        #  check the source files
Rename.pl -r *.f90        #  do replacements'
#EOC
,"\n";

    exit;
}

require $Input or die "Could not read input file $Input !!!\n";

@oldname=sort keys %newname;

# Define upper case version of old names as the keys of %NEWNAME
foreach $oldname (@oldname){
    $NEWNAME{uc($oldname)}=$newname{$oldname};
}

# Check for multiple occurances of new names and produce inverse lookup
$error=0;
foreach $oldname (@oldname){
    $newname=$newname{$oldname};
    $NEWNAME=uc($newname);
    if($otheroldname=$OLDNAME{$NEWNAME}){
	print "ERROR: both $oldname \& $otheroldname --> $newname !!!\n";
	$error=1;
    }else{
	$oldname{$newname}=$oldname;
        $OLDNAME{$NEWNAME}=$oldname;
    }
}

exit if $error;

# Check for new names that would be renamed next time
foreach $oldname (@oldname){
    $newname=$newname{$oldname};
    $NEWNAME=uc($newname);

    if($newnewname=$NEWNAME{$NEWNAME}){
        $NEWNEWNAME=uc($newnewname);
        # Do not worry if only the case is different
	next if $NEWNEWNAME eq $NEWNAME;

        # Get the old name that should really be renamed into $newnewname 
        $oldnewnewname=$OLDNAME{$NEWNEWNAME};

	if($oldnewnewname eq $newname){
            # Even the cases are the same
	    print "WARNING: $oldname  --> $newname -->  $newnewname !!!\n"
		unless $Warning;
	}else{
	    print "Warning: $oldname --> $newname, ".
		"$oldnewnewname --> $newnewname !!!\n" unless $Warning;
	    $DANGER{$NEWNAME}=1;
	}
    }
}

if($List){
    # list old and new names in alphabetical order 
    foreach $oldname (@oldname){
	print "$oldname => $newname{$oldname}\n";
    }
}

if($Rename or $Check){
    &rename(%newname);
}elsif($Undo){
    &rename(%oldname);
};

exit;

##########################################################################
sub rename{

    %rename=@_;

    @old=sort keys   %rename;
    @new=sort values %rename;

    while($file=shift(@ARGV)){
	if(not open(FILE,$file)){
	    print "Error opening file $file !!!\n";
	    next;
	}
	read(FILE,$text,-s $file);
	close(FILE);
	print "old text=\n$text\n" if $Debug;

	# Protect lines with containing the '!DO NOT RENAME' string and 
	# case(' or case(" statements.
	$icase=0; @case=();
	while($text =~ s/^(.*\!\s*do\ not\ rename.*|
			   \s*case\s*\(\s*['"].*)/_\[CASE$icase\]_/imx){
	    print " replacing case $icase\n" if $Debug;
	    $case[$icase++]=$1;
	}

	foreach $oldname (@old){
	    $newname=$rename{$oldname};
	    next if lc($newname) eq lc($oldname); # Only capitalization changes
	    if($text=~/\b$newname\b/i){
		print "  warning: variable $newname occurs in file $file !\n"
		    unless $Warning;
	    }
	}

	# Take next file if check only
	next if $Check;
	
	print "Renaming variables...\n" if $Debug;

	# Replace old names with tokens of the form _[NUMBER]_
        $count=0;
	for($i=0;$i<=$#old;$i++){
            $oldname=$old[$i];
            if($Undo or $DANGER{uc($oldname)}){
                # Replace only with strict case agreement
		$count += ($text=~s/\b$oldname\b/_\[$i\]_/g);
	    }else{
		# Replace with relaxed case checking
		$count += ($text=~s/\b$oldname\b/_\[$i\]_/ig);
	    }
	}
	if($count==0){
	    print "No variables to rename in file $file\n" unless $Quiet;
	    next;
	}
	print "tok text=\n$text\n" if $Debug;
        # Replace tokens with new names
	for($i=0;$i<=$#old;$i++){
	    $text=~s/_\[$i\]_/$rename{$old[$i]}/g;
	}
	print "New text=\n$text\n" if $Debug;

	# Put back the case(' and case(" lines
	for($icase=0; $icase<=$#case; $icase++){
	    $text =~ s/_\[CASE$icase\]_/$case[$icase]/i;
	}

        # Replace the file with the modified text
	rename $file,"$file~";
	open(FILE,">$file");
	print FILE $text;
	close(FILE);

	print "Finished $count replacements in file $file\n" unless $Quiet;
	$countall += $count;
	$nfile += 1;
    }

    print "Finished  $countall replacement(s) in $nfile file(s)\n"
	if $countall and not $Quiet;
}
