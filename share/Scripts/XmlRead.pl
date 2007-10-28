#!/usr/bin/perl
no strict; # because of local

$Debug = 0;

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
	    push(@Xml, { 'type' => 't', 'content' => $text } );
	    warn "text = $text\n" if $Debug;
	}

	# Now read the XML comment or tag
	my $ElementRef;
	if(s/^(\w+)\s*//){
	    my $tag = $1; $_ = $';
	    $ElementRef->{type} = "e";
	    $ElementRef->{name} = $tag;

	    # Find closing > but ignore quoted text like if="$x>1"
	    my $attributes;
	    while( /([\'\">])/ ){
		my $q = $1; $attributes .= $`; $_ = $';

		last if $q eq ">"; # found matching >

		$attributes .= $q; # start of some quoted text

		# Check for matching quotation mark
		if( /$q/ ){
		    $attributes .= $`.$q; $_ = $';
		}else{
		    die "$ERROR: no matching $q for <$tag $attributes in $_\n";
		}
	    }

	    # Check if the last character is a /
	    my $endtag = ($attributes =~ s/\/$//);

	    warn "tag=$tag, attributes=$attributes, endtag=$endtag\n"
		if $Debug;

	    # Read attributes
	    while($attributes=~s/^\s*(\w+)\s*=\s*(\'[^\']*\'|\"[^\"]*\")\s*//){
		my $attribute = $1;
		my $value = $2; $value =~ s/^[\'\"]//; $value =~ s/[\'\"]$//;
		$ElementRef->{attrib}{$attribute} = $value;

		warn "attribute=$attribute, value=$value\n" if $Debug;
	    }

	    die "could not read tag $tag\'s attributes=$attributes\n" 
		if $attributes;

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
	push(@Xml, { 'type' => 't', 'content' => $_ } );
	warn "final text = $_\n" if $Debug;
    }

    # Return a reference to the @Xml array
    return [ @Xml ]; 
}

1;
