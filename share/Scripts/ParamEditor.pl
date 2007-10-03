#!/usr/bin/perl -s
#
# *******************************************************************
# *                                                                 *
# *  WEBLINK  Copyright (C) 1996 Jozsef Hollosi                     *
# *                                                                 *
# *  WEBLINK Version 1.04                                           *
# *                                                                 *
# *  All Rights Reserved. Unpublished rights reserved under         *
# *  the copyright laws of the United States.                       *
# *                                                                 *
# *  The software contained on this media is proprietary to and     *
# *  embodies the confidential technology developed by J. Hollosi.  *
# *  Possession, use, duplication or dissemination of the software  *
# *  and media is authorized only pursuant to a valid written       *
# *  license agreement from Jozsef Hollosi.                         *
# *                                                                 *
# *  Specifically, WEBLINK can be freely modified and distributed   *
# *  in original or modified form with the VAC simulation package,  *
# *  as long as this copyright notice is unchanged.                 *
# *                                                                 *
# *******************************************************************
#
# FOR WEBLINK USERS:
#
# To run weblink as an interface type:
#
#   weblink.pl [-browser=browser] [-stdout=xterminalname]
#              [-srcfile=sourcefile] [-default=defaultfile] [-basedir=dirname]
#
# FOR WEBLINK PROGRAMMERS:
#
# Weblink is a standalone webserver/runtime environment. It works
# as a standard webserver, with an extension server-side language,
# called HyPerlText. It enables you to insert Perl statements directly
# into your HTML text.
#
# Syntax:
#
# eval:       <#@ expr#>
# eval+print: <#= expr#>
# if:         <#? expr#> ... <#/?#>
# if not:     <#! expr#> ... <#/!#>
# loop:       <#loop#> ... <#break expr#> ... <#/loop#>
#
# limitations: ifs and ifnots can be embedded, but loops can't
#
# variables:
# - the system keeps infinite state of all variables during the session
#   EXCEPT:
#           %FORM and @FORM has the latest submitted form variables
#
# - all variables and function names starting with _ are reserved
#          $_agent is either netscape or mosaic
#
# - all other variables are controlled by the user
#
# - subroutine for running external programs in the background:
#   &_background($command,$nohup);
#   where $nohup=0 means that the job is killed when you exit weblink,
#   while $nohup=1 lets the job run further.
#
# --------------------------------------------------------------------------

# Put switches which come from command line into _VARIABLES

$_BROWSER = $browser;
$_STDOUT_XTERM = $stdout;
$_BASEDIR = ($basedir or ".");
$_SRCFILE = ($srcfile or "index.php");

# Execute the $default file which sets some variables, among others 
# the $stdout variable which contains the name of the xterminal that 
# should show the STDOUT. Empty string implies stdout into the original window.

$default = "Weblink/default" unless $default;
if(-e $default){
    require $default;
}
else {
    # print STDERR "Warning: Customization file $default was not found!\n";
}
unless(-e $_SRCFILE){
    die "Source file $srcfile was not found!\n";
}
# Put switches which come from $default into _VARIABLES unless already set

$_BROWSER = $browser unless $_BROWSER;
$_STDOUT_XTERM = $stdout unless $_STDOUT_XTERM;

$_STDOUT_XTERM = "original" unless $_STDOUT_XTERM; # To be safe.

# Start work

&_system_setup;

&_start_browser;


while (1) {
    accept(_SOCK,"websocket");
    $_SOCK = 1;
    select(_SOCK); $| = 1;
    select(STDOUT); $| = 1;

    $_file = &_read_input;

    #&_debug("file = $_file\n");

    &_process($_file);

    shutdown(_SOCK, 2);
    close(_SOCK);
}

exit(0);

# -------------------------------------------------------------------

sub _system_setup
{
    $|=1;
    $SIG{'PIPE'} = '_sigpipe';
    $SIG{'CHLD'} = 'IGNORE';
    &_redirect_stdout($_STDOUT_XTERM);
}

sub _start_server
{
    $_websocket = &_open_socket("websocket", 3131);
    if ( fork() ) {

	print "Location: http://". $ENV{"SERVER_NAME"}. ":$_websocket/\n\n";
    	exit;
    }
### foreach (3..16) { syscall(6, $_); }
    close(STDIN); close(STDOUT); close(STDERR);
    alarm(300);
}

sub _start_browser
{
    if($ENV{"OSTYPE"} eq "darwin"){
	# On the MAC use the generic open statement
	$_browser = 'open';
    }else{
	# On other platforms find the browser
	@_bro=("netscape","mozilla"); unshift(@_bro,$_BROWSER) if $_BROWSER;
      TRY: 
	foreach $_bro (@_bro){
	    foreach $_dir (split(/:/, $ENV{"PATH"}.
				 ":/usr/bin/X11".
				 ":/usr/local/bin".
				 ":/usr/local/gnu/bin".
				 ":/opt/bin".
				 ":/opt/gnu/bin"
				 )) {
		if ( -x $_dir."/$_bro" ) {
		    $_browser = $_dir."/$_bro";
		    last TRY;
		}
	    }
	}
    }

    if ( defined($_browser) ) {
	print STDERR "Starting browser $_browser...\n";
	$_websocket = &_open_socket("websocket", 3131);
	if ( fork() ) {
	    return;
	}
	else {
	    if ( fork() ) {
		sleep(10);
		exit;
	    }
	    else {
		exec("$_browser http://127.0.0.1:$_websocket");
	    }
	}
    }
    else {
	print STDERR "ERROR: browser $_BROWSER not found.\n";
	exit(1);
    }
}

sub _redirect_stdout
{
    local($_XTERM)=@_;

    # find the new $_XTERM program

    if ( $_XTERM ne "original") {
	foreach $_dir (split(/:/, $ENV{"PATH"}.
			     ":/usr/bin/X11:/usr/X11R6/bin".
			     ":/usr/X11R5/bin:/usr/X386/bin")) {
	    if ( -x $_dir."/$_XTERM" ) {
		$_xterm = $_dir."/$_XTERM";
		last;
	    }
	}

	if ( defined($_xterm) ) {
	    if ( ! $_xtermsocket ) {
		$_xtermsocket = &_open_socket("xterm", 3114);
	    }
	    if ( fork() ) {
		# nothing
	    }
	    else {
		if ( fork() ) {
		    sleep(10);
		    exit;
		}
		else {
		    exec($_xterm, "-e","telnet", "localhost", $_xtermsocket);
		}
	    }
	    accept(STDOUT, "xterm");
	}
	else {
	    print "WARNING: $_XTERM is not found in your PATH\n";
	}
    }
    open(STDERR, ">&STDOUT");
    select(STDERR); $|=1; select(STDOUT); $|=1;
    $_STDOUT = 1;
    &_stdout("You will see STDOUT and STDERR from WEBLINK here.\n");
}

sub _getline
{
    local($_str, $_ch);
    while (1) {
	if ( ! sysread(_SOCK, $_ch, 1) ) { return($_str); }
	$_str .= $_ch;
	if ( $_ch eq "\n" ) { return($_str); }
    }
}

sub _sigpipe
{
    #&_debug("SIGPIPE called\n");
    $_SOCK = 0;
}

sub _sigchld
{
    #&_debug("SIGCHLD called\n");
}

sub _read_input
{
    #&_debug("read_input\n");
    ($_type, $_file) = split(" ", &_getline);
    #&_debug("type = $_type\n");
    #&_debug("file = $_file\n");

    $_CONTENT_LENGTH = 0;
    while ( $_line = &_getline ) {
	#&_debug("line read = $_line\n");
	$_line =~ s/[\n\r]+$//;
	if ( $_line =~ /^$/ ) { last; }
	if ( $_line =~ /^User-Agent:\s+(.+)/i ) {
	    $_agent    = $1;
	    if ( $_agent =~ /^Mozilla.(\d)/ ) {
		$_agentversion = $1;
		if ( $_agentversion >= 2 ) {
		    $_agent = "netscape2";
		}
		else {
		    $_agent = "netscape";
		}
	    }
	    else {
		$_agent = "mosaic";
	    }
	}
	if ( $_line =~ /^Content-Length:\s+(\d+)/i ) {
	    $_CONTENT_LENGTH = $1;
	}
    }

    #&_debug("finished header\n");

    $_getform  = '';
    $_postform = '';
    if ( $_file =~ /^(.*)\?(.*)$/ ) {
	$_file = $1;
	$_getform = $2;
    }
    if ( $_CONTENT_LENGTH > 0 ) {
	#&_debug("before post content read $_CONTENT_LENGTH\n");
	read(_SOCK, $_postform, $_CONTENT_LENGTH);
	#&_debug("read post content len=".length($_postform)."\n");
    }
    undef(@FORM);
    undef(%FORM);
    &_formvalues( $_getform . "&" . $_postform );
    return($_file);
}

sub _formvalues
{
    if ( $_[0] eq "&" ) { return; }
    local($_str) = $_[0];
    local($_var, $_val);
    $_str =~ s/&/&&/g;
    $_str = "&".$_str."&";
    while ( $_str =~ /&([\w\d]+)=([^&]*)&/g ) {
        $_var = $1;
        $_val = $2;
        push(@FORM, $_var);
        $FORM{$_var} = $_val;
    }
    foreach $_var ( @FORM ) {
        $FORM{$_var} = &_unhtml($FORM{$_var});
    }
}

sub _unhtml
{
    local($_str);
    $_str = '';
    $_[0] =~ s/\+/ /g;
    foreach $_piece ( split(/(%..)/, $_[0]) ) {
        if ( $_piece =~ /%(..)/ ) {
            $_num = hex($1);
            $_str .= sprintf("%c", $_num);
        }
        else {
            $_str .= $_piece;
        }
    }
    #TOTH: replace web newline "\x0D\x0A" by Perl newline "\n"
    $_str=~s/\x0D\x0A/\n/g;

    return($_str);
}

sub _print_result_ok
{
    &_printout("HTTP/1.0 200 Document follows\n");
    &_printout("MIME-Version: 1.0\n");
    &_printout("Server: WEBLINK/1.0\n");
    &_printout("Content-Type: $_[0]\n");
    if ( $_[1] > 0 ) {
	&_printout("Content-Length: $_[1]\n");
    }
    &_printout("\n");
}

sub _print_result_error
{
    &_printout("HTTP/1.0 404 Not found\n");
    &_printout("MIME-Version: 1.0\n");
    &_printout("Server: WEBLINK/1.0\n");
    &_printout("Content-Type: text/html\n\n");
    &_printout("
<HTML>
<HEAD>
<TITLE>Error</TITLE>
</HEAD>
<BODY>
<H1>Error 404</H1>

Not found - file does not exist or is read protected

</BODY>
</HTML>
");
}

sub _open_socket
{
    $_LSOCK = $_[0];
    local($_tryport) = $_[1];
    $_AF_INET     = 2;
    $_PF_INET     = 2;

    $_SOL_SOCKET     = 0xffff;
    $_SO_REUSEADDR   = 4;
    $_sockaddr = 'S n a4 x8';
    $_pat = 'S n C4 x8';
    ($_name, $aliases, $proto) = getprotobyname('tcp');
    while (1) {
	$_thisport = pack($_pat, $_AF_INET, $_tryport, 127,0,0,1);
	# Try to open socket streams 1, 2 or die
	socket($_LSOCK, $_PF_INET, 1, $proto) ||
	socket($_LSOCK, $_PF_INET, 2, $proto) ||
	die "Error: could not open socket.\n";

	setsockopt($_LSOCK, $_SOL_SOCKET, $_SO_REUSEADDR, "\001\0\0\0");
	if ( bind($_LSOCK,$_thisport) ) { last; }
	$_tryport++;
    }
    listen($_LSOCK,5);
    return($_tryport);
}

sub _process
{
    #&_debug("process = $_[0]\n");

    $_SOURCE = "$_BASEDIR/$_[0]";
    if ( -d $_SOURCE ) {
	$_SOURCE .= "$_SRCFILE";
	if ( ! -f $_SOURCE ) {
	    $_SOURCE =~ s/\/+[^\/]+$//;
	    opendir(DIR, $_SOURCE);
	    foreach (sort grep(!/^\./, readdir(DIR))) {
		&_printout("<a href=$_");
		if ( -d "$_SOURCE/$_" ) {
		    &_printout("/");
		}
		&_printout(">$_</a><br>\n");
	    }
	    closedir(DIR);
	    return;
	}
    }

    if ( ! open(_SOURCE, $_SOURCE) ) {
	&_print_result_error;
	return;
    }

    if ( $_SOURCE =~ /(gif|jpeg)$/ ) {
	&_print_result_ok("image/$1", 0);
	while ( read(_SOURCE, $_readbuffer, 10000) ) {
	    &_printout($_readbuffer);
	}
	return;
    }
    elsif ( $_SOURCE =~ /\.(php|html?)$/ ) {
	&_print_result_ok("text/html", 0);
    }
    else {
	&_print_result_ok("text/plain", 0);
    }

    $_part = '';
    while ( $_line = <_SOURCE> ) {
	&_processline($_line);
    }
    close(_SOURCE);
}

sub _processline
{
    local($_line) = $_[0];

    foreach $_token ( split(/(<#|#>)/, $_line) ) {
	if ( $_token eq "<#" ) {
	    $_part = $_token;
	    $_inpart = 1;
	    next;
	}
	elsif ( $_token eq "#>" ) {
	    $_part .= $_token;
	    $_inpart = 0;
	}
	else {
	    if ( $_inpart ) {
		$_part .= $_token;
		next;
	    }
	    $_part = $_token;
	}

	if ( $_part =~ /^<#break(\s+\S.*)?#>$/i ) {
	    $_args = $1;
	    if ( $_readloop ) { $_loopcontent .= $_part; next; }
	    if ( $_falsemode ) { next; }
	    if ( ! $_inloop ) { next; }
	    if ( ( length($_args) > 0 ) && ( ! &_val($_args) ) ) { next; }
	    $_inloop = 0;
	    #&_debug("breaking out\n");
	    return;
	}
	if ( $_part =~ /^<#\/loop#>$/i ) {
	    $_readloop = 0;
	    $_inloop = 1;
	    if ( $_falsemode ) { next; }
	    while ($_inloop) {
		&_processline($_loopcontent);
	    }
	    next;
	}
	if ( $_readloop ) {
	    $_loopcontent .= $_part;
	    next;
	}
	if ( $_part =~ /^<#loop#>$/i ) {
	    $_readloop = 1;
	    $_loopcontent = '';
	    next;
	}
	if ( $_part =~ /^<#\/.?#>$/ ) {
	    if ( $_level > 0 ) { $_level--; }
	    $_falsemode = $_falsemode[$_level];
	    next;
	}
	#&_debug("PART ==$_part==\n");
	if ( $_part =~ /^<#(.)\s+((.|\n)+)#>$/ ) {
	    $_action = $1;
	    $_args   = $2;
	    #&_debug("action=$_action args=$_args\n");
	    if ( ( $_action eq "?" ) || ( $_action eq "!" ) ) {
		$_level++;
		if ( $_falsemode ) {
		    ;
		}
		elsif ( ( $_action eq "?" ) && ( ! &_val($_args) ) ) {
		    $_falsemode = 1;
		}
		elsif ( ( $_action eq "!" ) && ( &_val($_args) ) ) {
		    $_falsemode = 1;
		}
		$_falsemode[$_level] = $_falsemode;
		next;
	    }
	    if ( $_action eq "=" ) {
		if ( ! $_falsemode ) { &_printout(&_val($_args)); }
		next;
	    }
	    if ( $_action eq "@" ) {
		if ( ! $_falsemode ) { &_val($_args); }
		next;
	    }
	}
	if ( $_part =~ /^<#/ ) { next; }
	if ( ! $_falsemode ) { &_printout($_part); }
    }
}

sub _val
{
    local($_str) = $_[0];
    
    $_str =~ s/^\s+//;
    $_str =~ s/\s+$//;
    
    #&_debug("val = $_str\n");
    return(eval $_str);
}

sub _printout
{
    if ( $_SOCK ) {
	syswrite(_SOCK, $_[0], length($_[0]));
	#&_debug("printout: >>> $_[0] <<<\n");
    }
}

sub _stdout
{
    local($_numwritten);
    if ( $_STDOUT ) {
	$_numwritten = syswrite(STDOUT, $_[0], length($_[0]));
	if ( $_numwritten != length($_[0]) ) {
	    close(STDOUT);
	    $_STDOUT = 0;
	    #&_debug("closed xterm\n");
	}
    }
}

sub _debug
{
    print STDERR "debug: ", $_[0] if $debug;
}

# ---------- special user subroutine ----------

sub _background
{
    local($command,$nohup)=@_;

    #&_debug("background: $command\n");
    if ( fork() ) {
	return;
    }
    else {
	if ( fork() ) {
	    sleep(10);
	    exit;
	}
    }

    setpgrp(0,$$) if $nohup;
    exec($command);
    exit;
}

############## END ############
