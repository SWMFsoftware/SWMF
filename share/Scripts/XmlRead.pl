#!/usr/bin/perl -s

my $Debug = $D;

use strict;

my @Result;
@Result = &XmlRead( join('',<>) );

exit 0;
##############################################################################
sub XmlRead{

    local($_) = shift;  # XML text

    no strict;
    local(@Xml); # Array of XML elements

    # Remove all comments
    while( /<!--/ ){
	my $before = $`;
	/-->/ or die "missing -->\n";
	$_ = $before . $';
    }

    while( s/^([^<]*)<//){

	# check for text before next <
	my $text = $1;
	if($text){
	    my $ElementRef;
	    $ElementRef->{type} = "t";
	    $ElementRef->{body} = $text;
	    warn "text = $text\n" if $Debug;
	    push(@Xml, $ElementRef);
	}

	# Now read the XML comment or tag
	my $ElementRef;
	if(s/^(\w+)\s*//){
	    my $tag = $1; $_ = $';
	    $ElementRef->{type} = "e";
	    $ElementRef->{name} = $tag;

	    s/^([^>]*)>// or die "<$tag without closing >\n";

	    my $attributes = $1; $_ = $';

	    # Check if the last character is a /
	    my $endtag = ($attributes =~ s/\/$//);

	    print "tag=$tag, attributes=$attributes, endtag=$endtag\n"
		if $Debug;

	    # Read attributes
	    while($attributes=~s/^\s*(\w+)\s*=\s*(\'[^\']*\'|\"[^\"]*\")\s*//){
		my $attribute = $1;
		my $value = $2; $value =~ s/^[\'\"]//; $value =~ s/[\'\"]$//;
		$ElementRef->{attrib}{$attribute} = $value;

		print "attribute=$attribute, value=$value\n" if $Debug;
	    }

	    die "could not read tag $tag\'s attributes=$attributes\n" 
		if $attributes;

	    if(not $endtag){
		# Check for nested <tag...> of the same type and </tag>
		my $content;
		my $level = 1; # there was an opening tag so far
		while(/(<$tag\b|<\/$tag>)/){
		    $content .= $`; $_ = $';
		    my $tag2 = $1;
		    if($tag2 =~ /<\//){
			$level--;
			last unless $level; # exit if we got back to level 0
		    }else{
			$level++;
		    }
		    $content .= $tag2;
		}
		die "missing </$tag> near".substr($content,0,100)."\n"
		    if $level;

		$ElementRef->{content} = &XmlRead($content) if $content;
	    }
	}else{
	    die "incorrect <\n";
	}
	push(@Xml, $ElementRef);
    }

    # Remove trailing space
    s/\s*$//;
    
    # Check for closing text
    if($_){
	my $ElementRef;
	$ElementRef->{type} = "t";
	$ElementRef->{body} = $_;
	warn "final text = $_\n" if $Debug;
	push(@Xml, $ElementRef);
    }

    return @Xml;
}
