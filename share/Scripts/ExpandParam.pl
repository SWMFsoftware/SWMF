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

# Default filename
$input='PARAM.in';

# Overwrite with argument if it exists
$input=$ARGV[0] if $ARGV[0];

# Expand the file recursively
&process_file($input, 'fh00');

exit;
##############################################################################

sub process_file {
    local($filename, $input) = @_;
    $input++;
    open($input, $filename) || die"Can't open $filename: $!\n";
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
    close $input;
}
