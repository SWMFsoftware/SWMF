#!/usr/bin/perl
#^CFG COPYRIGHT UM

#BOP
#!ROUTINE: ExpandParam.pl
#!DESCRIPTION:
# Expands a parameter file by including all the nested include files.
# If an include file is not found, an error is produced.
# Usage:
#\begin{verbatim}
# ExpandParam.pl > PARAM.in
# ExpandParam.pl PARAM.in > PARAM.expand
#\end{verbatim}
#!REVISION HISTORY:
# 01/20/2002 G. Toth - initial version
#EOP

$ERROR = "ExpandParam.pl_ERROR";

# Default filename
$basefile='PARAM.in';

# Overwrite with argument if it exists
$basefile=$ARGV[0] if $ARGV[0];

# Extract directoryname if any

if($basefile =~ /\/([^\/]+)$/){
    $dir = $`;
    $basefile = $1;
    chdir $dir or die "$ERROR: could not cd $dir\n";
}

# Expand the file recursively
&process_file($basefile, 'fh00');

exit;
##############################################################################

sub process_file {
    local($filename, $input) = @_;
    $input++;
    open($input, $filename) or die"$ERROR: can't open $filename: $!\n";
    while (<$input>) {
	# Stop reading if #END command is read
        last if /^#END\b/;

	# Check for #INCLUDE
        if (/^#INCLUDE\b/) {
	    # Read file name following #INCLUDE
            $includefile=<$input>;
	    # process include file recursively
	    &process_file($includefile,$input);
        }else{
	    # Print line as it is otherwise
	    print;
	}
    }
    if($basefile eq $filename){
	print "\#END\n",<$input>;
    }
    close $input;
}
