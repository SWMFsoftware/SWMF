#!/usr/bin/perl
no strict; # because of local

$Debug = 0;

# XML special characters
my %unxml = ( 'gt' => '>',
	      'lt' => '<',
	      'amp'=> '&',
	      'quot' => '"',
	      'apos' => "'",
	      );

sub XmlRead{

    local($_) = shift;  # XML text
    local(@Xml); # Array of XML elements

    # Remove all comments
    while( /<!--/ ){
	my $before = $`;
	/-->/ or die "missing -->\n";
	$_ = $before . $';
    }

    # Find XML elements <...>
    while( s/^([^<]*)<//){

	# check for text before next <
	my $text = $1;
	if($text){
	    $text =~ s/\&(\w+);/$unxml{$1}/g;
	    push(@Xml, { 'type' => 't', 'content' => $text } );
	    warn "text = $text\n" if $Debug;
	}

	# Now read the XML comment or tag
	my $ElementRef;
	if(s/^(\w+)\s*//){
	    my $tag = $1; $_ = $';
	    $ElementRef->{type} = "e";
	    $ElementRef->{name} = $tag;

	    # Read attributes
	    while( s/^\s*(\w+)\s*=\s*(\'[^\']*\'|\"[^\"]*\")\s*// ){
		my $attribute = $1;
		my $value = $2; $value =~ s/^[\'\"]//; $value =~ s/[\'\"]$//;
		$value =~ s/\&(\w+);/$unxml{$1}/g;

		$ElementRef->{attrib}{$attribute} = $value;

		warn "attribute=$attribute, value=$value\n" if $Debug;
	    }

	    # Check if the last character is a /
	    my $endtag = (s|^/||);

	    # Check if the closing > is there
	    s/^>// or die "$ERROR: missing > for <$tag ...\n";

	    warn "tag=$tag, endtag=$endtag\n" if $Debug;

	    if(not $endtag){
		# Check for nested <tag...> of the same type and </tag>
		my $content;
		my $level = 1; # there was an opening tag so far
		while(/(<$tag\b[^>]*>|<\/$tag>)/){
		    $content .= $`; $_ = $';
		    my $tag2 = $1;
		    if($tag2 eq "</$tag>"){
			$level--;
			last unless $level; # exit if we got back to level 0
		    }elsif($tag2 !~ /\/>$/){
			$level++;
		    }
		    $content .= $tag2; # still inside, keep searching for end
		}
		die "missing </$tag> near".substr($content,0,100)."\n"
		    if $level;

		$ElementRef->{content} = &XmlRead($content) 
		    if length($content);
	    }
	}else{
	    die "incorrect <\n";
	}
	push(@Xml, $ElementRef);
    }

    # Check for closing text
    if( $_ ne "\n" ){
	s/\&(\w+);/$unxml{$1}/g;
	push(@Xml, { 'type' => 't', 'content' => $_ } );
	warn "final text = $_\n" if $Debug;
    }

    # Return a reference to the @Xml array
    return [ @Xml ]; 
}

1;
