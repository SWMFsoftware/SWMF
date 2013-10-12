#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
# Convert/test the endianness of VAC/BATSRUS binary data/plot file(s).
# Convert from or to Cray format (8-byte) integers.

# Store input flags in long variable names
my $Help    = ($h or $help);
my $Test    = ($t or $test);
my $Convert = ($c or $convert);
my $Inplace = ($i or $inplace);
my $Endian  = ($e or $endian);
my $tocray  = ($C or $cray);
my $Quiet   = ($q or $quiet);

### use strict;

my $nArg = @ARGV;

my $zero=pack('L',0);

undef $/; # ignore end of line

my $HELP="\nType 'FixEndian.pl -h' for help!\n";
my $ERROR="ERROR in FixEndian.pl";

&print_help if $Help;

die "$ERROR: use exactly one of the -t(est), -i(nplace) or -c(onvert) flags!".
    $HELP unless $Convert + $Inplace + $Test == 1;

die "$ERROR: the -i(nplace) flag requires at least one file name!$HELP"
    if $Inplace and $nArg == 0;

die "$ERROR: the -c(onvert) flag requires exactly two filenames!$HELP"
    if $Convert and $nArg != 2;

# Global variables
my  $B79 = unpack 'N', pack('L',79);   # Big endian value
my  $L79 = unpack 'V', pack('L',79);   # Little endian value

my $B500 = unpack 'N', pack('L',500);   # Big endian value
my $L500 = unpack 'V', pack('L',500);   # Little endian value


my $machine = &test_machine;        # endianness of this computer

if($Convert){
    &converttofile($ARGV[1]) if &readhead($ARGV[0]);
}else{
    my $file;
    while($file=shift(@ARGV)){
	if(&readhead($file) and $Inplace){
	    if(&converttofile("$file~")){
		rename("$file~",$file) or
		    die "Could not move $file~ to $file: $!\n";
		print "\nMoved $file~ to $file\n\n" unless $Quiet;
	    }
	}
    }
}

print "\nScript FixEndian.pl finished running\n\n" 
   if $nArg > 1 and not $Quiet and not $Convert;

exit;

##############################################################################
sub print_help{

print "
Purpose:

   Show or convert the endianness of binary VAC/BATSRUS type data/plot files 
   and/or convert between 4 and 8-byte integer formats.

Usage: 

   FixEndian.pl -t [-q] [filename1 filename2...]
   FixEndian.pl -i [-q] [-e=ENDIAN] [-C] filename1 [filename2...]
   FixEndian.pl -c [-q] [-e=ENDIAN] [-C] filefrom fileto

-t -test     Test and show endianness of the machine and listed file(s) if any.
             Also show header information for the file(s) unless in quiet mode.

-i -inplace  Convert file(s) in-place (uses FILENAME~ as a temporary file)

-c -convert  Convert filefrom to fileto.

-q -quiet    Suppress printing file header and other verbose information.

-e=ENDIAN    Convert to the endianness defined by ENDIAN. Possible values are:
             'big', 'little', 'machine', 'not-machine', 'keep', 'swap'
             Only the first character is significant. The value 'machine'
             means that the file is converted to the endianness of the 
             machine the script is running on, while 'not-machine' converts
             to the opposite endianness of the machine. The value 'keep'
             keeps the endianness of the file (but may convert the
             integers from 4 to 8 bytes or vice-versa), while the value 'swap'
             changes the endianness of the file.
             Default is to 'swap' the endianness.

-C -cray     Convert integers to Cray format (8-byte).

Examples:

Test the endianness (and Cray format) of multiple files:

   FixEndian.pl -t *.out

Convert a file to another file with opposite endianness:

   FixEndian.pl -c oldfilename newfilename

In-place conversion of the endiannes to the machine's endianness 
for multiple files in quiet mode:

   FixEndian.pl -i -q -e=m *.out

Keep the endiannes and convert between Cray (8-byte) and normal (4-byte)
integer formats: 

   FixEndian.pl -e=k -cray -c oldfilename crayfilename
   FixEndian.pl -e=k -c       crayfilename newfilename
";
exit;
}

##############################################################################
sub test_machine{

    if(79 == $B79){
	$machine='big';
    }elsif(79 == $L79){
	$machine='little';
    }else{
	die "Machine is neither little nor big endian type!\n";
    }
    print "\nThis is a $machine endian machine.\n\n" unless $Quiet;
    return $machine;

}

##############################################################################
sub readhead{

    my $file = shift;

    open(FROM,"$file") or die "Could not open $file: $!\n";

    if(-T $file){
	print "$file seems to be an ASCII file!\n"; return(0);
    }

    if(200 > -s $file){
	print "$file is too small to be a binary VAC/BATSRUS data/plot file!\n";
	return(0);
    }

###    my $line1;

    read FROM, $line1, 4;  # read len1

###    my $l1;
###    my $headline;
###    my $l2;

    ($l1) = unpack 'L', $line1;

    if(    $l1 == $B79){
	$from='big';    $strlen =  79
    }elsif($l1 == $L79){
	$from='little'; $strlen =  79
    }elsif($l1 == $B500){
	$from='big'; $strlen = 500
    }elsif($l1 == $L500){
	$from='little'; $strlen = 500
    }else{
	print "l1=$l1 B79=$B79 L79=$L79 B500=$B500 L500=$L500\n";
	print "$file is not a VAC/BATSRUS binary data/plot file in ".
	    "UNIX F77 format!\n";
	return(0);
    }

    if($Endian =~ /^b/){
	$to = 'big';
    }elsif($Endian =~ /^l/){
	$to = 'little'
    }elsif($Endian =~ /^m/){
	$to = $machine;
    }elsif($Endian =~ /^n/){
	$to = ($machine eq 'big') ? 'little' : 'big';
    }elsif($Endian =~ /^k/){
	$to = $from;
    }else{
	$to = ($from eq 'big') ? 'little' : 'big';
    }

    # Keep endianness?
    $Keep = ($from eq $to);

    $convertfrom = $from ne $machine;
    $convertto   = $to   ne $machine;

    # read rest of the first line and the second record length
    read FROM, $line1, $strlen+8;  # read headline,len1,len2
    ($headline,$l1,$l2)=unpack 'a'.$strlen.'LL', $line1;

    $len1=pack("L",$strlen);
    $len1=&convstr($len1) if $convertto;        # length of line1 for output

    $ll2= $convertfrom ? &convint($l2) : $l2;   # length of line2 from input

    # Decide file type based on the length of the 2nd input line
    if($ll2 == 24){
	$precision="double"; $r='d'; $rbyte=8; $ibyte=4; $double=1; $cray=0;
    }elsif($ll2 == 20){
	$precision="single"; $r='f'; $rbyte=4; $ibyte=4; $double=0; $cray=0;
    }elsif($ll2 == 40){
	$precision="double"; $r='d'; $rbyte=8; $ibyte=8; $double=1; $cray=1;
	$type="Cray ";
    }elsif($ll2 == 36){
	$precision="single"; $r='f'; $rbyte=4; $ibyte=8; $double=0; $cray=1;
	$type="Cray ";
    }else{
	print "ll2=$ll2\n";
	print "$file is not a VAC/BATSRUS binary data/plot file ".
	    "in UNIX F77 format!\n";
	return(0);
    }

    $ibyteout= $tocray ? 8 : 4;
    $typeout=  $tocray ? "Cray " : "";

###    my $line2;
    read FROM, $line2, $ll2;           # read line2 
    read FROM, $l2, 4;                 # read closing length

    $len2=pack("L",4*$ibyteout+$rbyte); # calculate length of line 2 for output
    $len2=&convstr($len2) if $convertto;

    &line2fromcray if $cray;        # Remove extra bytes for Cray input file
    $line2conv=&convstr($line2);    # Convert endiannes
    &swaptbytes if $double;         # Swap bytes of time variable "T" if needed
    $line2_= $Keep ? $line2 : $line2conv; # Line 2 for output
    &line2tocray if $tocray;        # Insert extra bytes for Cray output file

    # Read values from line 2

    ($it,$t,$ndimini,$neqpar,$nw)=unpack "l${r}lll",
    $convertfrom ? $line2conv : $line2; 

    $ndim=abs($ndimini);

    $ll3=$ibyte*$ndim;
    read FROM, $l3, 4;                 # len3
    read FROM, $line3, $ll3;           # line3: n1, n2..
    read FROM, $l3, 4;                 # len3

    $len3=pack("L",$ibyteout*$ndim);   # calculate length of line 3 for output
    $len3=&convstr($len3) if $convertto;

    &line3fromcray if $cray;              # remove extra Cray bytes
    $line3conv=&convstr($line3);          # convert endiannes
    $line3_= $Keep ? $line3 : $line3conv; # line 3 for output
    &line3tocray if $tocray;              # insert extra Cray bytes

    @nx= unpack "l*", $convertfrom ? $line3conv : $line3; # Read nx from line 3
    $nxs=1; foreach $nx (@nx){$nxs*=$nx};

    $ll4=$rbyte*$neqpar;               # length of line4 from input
    read FROM, $l4, 4;
    read FROM, $line4, $ll4;           # line4: eqpar
    read FROM, $l4, 4;

    $len4=pack("L",$ll4);    # calculate length of line 4 for output
    $len4=&convstr($len4) if $convertto;

    $line4conv =
	$double ? &convdble($line4) : &convstr($line4); # convert endianness
    $line4_= $Keep ? $line4 : $line4conv;               # line 4 for output
    @eqpar= unpack "$r*", $convertfrom ? $line4conv : $line4;  # read eqpar

    read FROM, $l5, 4;
    read FROM, $varnames, $strlen;              # line5: varnames
    read FROM, $l5, 4;

    $len5=$len1;                                # length of line 5 

    $llx=$rbyte*$nxs*$ndim;                      # length of the X array
    $lenx=pack("L",$llx);                        # length of X for output
    $lenx=&convstr($lenx) if $convertto;

    $llw=$rbyte*$nxs;                            # length of W(*,*..,iw)
    $lenw=pack("L",$llw);                        # length of W for output
    $lenw=&convstr($lenw) if $convertto;

    print "$file is a $from endian $type$precision precision file.\n"
	if $Test or not $Quiet;

    unless($Quiet){
	print "\nheadline=",substr($headline,0,70),"\n";
	print "it=$it t=$t ndim=$ndimini neqpar=$neqpar nw=$nw\n";
	print "nx=",join(',',@nx),"\n";
	print "eqpar=",join(',',@eqpar),"\n";
	print "names=",substr($varnames,0,70),"\n\n";

	$pictsize=(6+$nw)*8+2*$strlen+(4+$ndim)*$ibyte+
	    (1+$neqpar+$nxs*($ndim+$nw))*$rbyte;
	$filesize= -s $file;
	print "file=$filesize bytes / snapshot=$pictsize bytes = ",
	$filesize/$pictsize," snapshot(s)\n\n";
    }

    # Check if there is anything to do
    return(0) if $Keep and $cray == $tocray;

    return(1);
}

##############################################################################
sub converttofile{

    local($file)=$_[0];
    open(TO,">$file") or print STDERR "Could not open $file: $!\n", return(0);

    print "$file will be a $to endian $typeout$precision precision file...\n\n"
	    unless $Quiet;

    &printhead;
    &dobody;
    &dorest;

    print "\nFinished conversion to $file\n" unless $Quiet;

    return(1);
}

##############################################################################
sub printhead{

    print TO   $len1,$headline,$len1;
    print TO   $len2,$line2_,$len2;
    print TO   $len3,$line3_,$len3;
    print TO   $len4,$line4_,$len4;
    print TO   $len5,$varnames,$len5;

}

##############################################################################
sub dobody{

    read FROM, $lx, 4;
    $llx = unpack("L",$lx);
    $llx = &convint($llx) if $convertfrom;

    read FROM, $xx, $llx;             # read x(*,*..)
    read FROM, $lx, 4;
    $lenx=pack("L",$llx);
    $lenx = &convstr($lenx) if $convertto;
    print TO $lenx;
    print TO $Keep ? $xx : ($double ? &convdble($xx) : &convstr($xx));
    print TO $lenx;

    for($iw=0; $iw<$nw; $iw++){
	read FROM, $lw, 4;
	$llw = unpack("L",$lw);
	$llw = &convint($llw) if $convertfrom;

	read FROM, $ww, $llw;         # read len,w(*,*..,iw),len
	read FROM, $lw, 4;
	$lenw = pack("L",$llw);
	$lenw = &convstr($lenw) if $convertto;
	print TO $lenw;
	print TO $Keep ? $ww : ($double ? &convdble($ww): &convstr($ww));
	print TO $lenw;
    }
}

##############################################################################
sub dorest{

    # Read rest of snapshots using the header info

    $ipict=0;
    print "Converted snapshot ",++$ipict,"\n" unless $Quiet;

    while(read FROM,$l1,4){                  # read len1

	read FROM, $headline, $strlen;       # line1: headline
	read FROM, $l1l2, 8;                 # read len1, len2
	read FROM, $line2, $ll2;             # line2: it,t,ndim,neqpar,nw 
	read FROM, $l2l3, 8;                 # read len2, len3
	read FROM, $line3, $ll3;             # line3: nx1, nx2..
	read FROM, $l3l4, 8;                 # read len3, len4
	read FROM, $line4, $ll4;             # line4: eqpar1,eqpar2..
	read FROM, $l4l5, 8;                 # len4,len5
	read FROM, $varnames, $strlen;       # line5: varnames
	read FROM, $l5, 4;                   # len5

	# Convert line 2

	&line2fromcray if $cray;
        if($Keep){
	    $line2_=$line2;
	}else{
	    $line2conv=&convstr($line2);
	    &swaptbytes if $double;
	    $line2_= $line2conv;
	}

	# Convert line 3

	&line3fromcray if $cray;                    # remove extra Cray bytes
	$line3_= $Keep ? $line3 : &convstr($line3); # convert endiannes
	&line3tocray if $tocray;                    # insert extra Cray bytes

	# Convert line 4

        if($Keep){
	    $line4_=$line4;
	}else{
	    $line4_= $double ? &convdble($line4) : &convstr($line4);
	}

	&printhead;
	
	&dobody;
	
	print "Converted snapshot ",++$ipict,"\n" unless $Quiet;
    }
}

##############################################################################
sub convint{
    unpack('L',reverse pack('L',$_[0]));
}

##############################################################################
sub convstr{
    local($i);
    $_=$_[0];
    for($i=0; $i<length(); $i+=4){
        substr($_,$i,4)=reverse substr($_,$i,4);
    }
    $_;
}

##############################################################################
sub convdble{
    local($i);
    $_=$_[0];

    for($i=0; $i<length(); $i+=8){
        substr($_,$i,8)=reverse(substr($_,$i,8));
    }

    $_;

}
##############################################################################
sub line2fromcray{
    # get rid of extra bytes of Cray format in line 2: it,t,ndim,neqpar,nw
    $line2=substr($line2,4,4).substr($line2,8,$rbyte).
	substr($line2,-20,4).substr($line2,-12,4).substr($line2,-4,4);
}
##############################################################################
sub line2tocray{
    # insert extra bytes of Cray format in line 2: it,t,ndim,neqpar,nw
    $line2_=$zero.substr($line2_,0,4).substr($line2_,4,$rbyte).
            $zero.substr($line2_,-12,4).$zero.substr($line2_,-8,4).
            $zero.substr($line2_,-4,4);
}
##############################################################################
sub line3fromcray{
    # get rid of extra bytes of Cray format in line 3 containing the nx array
    if($ndim==1){
	$line3=substr($line3,4,4);
    }elsif($ndim==2){
	$line3=substr($line3,4,4).substr($line3,12,4);
    }else{
	$line3=substr($line3,4,4).substr($line3,12,4).substr($line3,20,4);
    }
}
##############################################################################
sub line3tocray{
    # insert extra bytes of Cray format in line 3 containing the nx array
    if($ndim==1){
	$line3_=$zero.$line3_;
    }elsif($ndim==2){
	$line3_=$zero.substr($line3_,0,4).$zero.substr($line3_,4,4);
    }else{
	$line3_=$zero.substr($line3_,0,4).$zero.substr($line3_,4,4).
                $zero.substr($line3_,8,4);
    }
}
##############################################################################
sub swaptbytes{
    # swap 2 bytes of double precision time "T" in the converted line 2

    substr($line2conv,4,8)=substr($line2conv,8,4).substr($line2conv,4,4);
}
