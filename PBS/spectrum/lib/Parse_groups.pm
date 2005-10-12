package Parse_groups;
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(extract_q, cmu_to_chroma,make_prop_name); # Subroutines
our $VERSION = 1.00;		# Version number

sub extract_q{
#
#  First we have to remove the leading spaces.  Note that the up arrow
#  here denotes that this is at the beginning of the line
    $_[0] =~ s/^\s+//;
    ($number_prop, $snk[0], $snk[1], $snk[2], 
     $snk_len_offset, $src[0], $src[1], $src[2], $src_len_offset) =
	split(/\s+/,$_[0]);
#
# The src_len and snk_len are given in Colin's code as offsets
# We must turn them into "real" lengths
# Note that the one here comes from offset past the first argument.
#

#
#  Loop over the flavours of the quarks, recognising that isospin

    for($ctr = 0; $ctr < 3; $ctr++){

	$flav_in = $_[$ctr+1];

	if( $flav_in eq "u" || $flav_in eq "d"){
	    $flav[$ctr] = "ud";}
	else{
	    $flav[$ctr] = $flav_in;
	}
	print "Flav is $flav[$ctr]\n";
   }


    $src_len = $_[$src_len_offset+4];
    $snk_len = $_[$snk_len_offset+4];

#
#  Here we assume that the first two propagators are isospin invariant
#  Third quark propagator is s quark, and therefore might be different
#  in mass

#
# One subtlety is that, if sink or source are zero, we should put displacement
# lengths to be zero to avoid double-counting

    $ctr=0;

    while ( $ctr < 3){
	if ( $snk[$ctr] != 0){
	    $snk_out = $snk_len;
	}
	else{
	    $snk_out = 0;
	}

	if ( $src[$ctr] != 0){
	    $src_out = $src_len;
	}
	else{
	    $src_out = 0;
	}

	$prop[$ctr] = "$flav[$ctr] $snk[$ctr] $src[$ctr] $snk_out $src_out";
	$ctr = $ctr + 1;
    }

    return @prop;
}

sub cmu_to_chroma{
#
# Colin's output files and the chroma code expect the displacement
# lengths to be of the form +/- 1, 2, 3 and len, whilst the chroma
# input file expects them to be of the form 0, 1, 2 +/len
# Thus we need to translate between the two
#
    $dir_in=$_[0];
    $dir_out=$dir_in;
    $dir_out =~ s/-//;		# Remove any minus signs
    if($dir_out != 0){
	$dir_out--;			# decrement by one
    }

    $len_in=$_[1];
    $len_out=$len_in;

    if ($dir_in < 0){		# Negative directions go to negative length
	$len_out=-$len_out;
    }
    if ($dir_in == 0){		# No derivative is zero length
	$len_out=0;
    }
    
    return ($dir_out, $len_out);
}

sub make_prop_name{
#
#  This subroutine appends the source and sink information to make
#  the propagator name
#  Note that the information is in CHROMA FORMAT
    my $prop_root = $_[0];
    my $flav = $_[1];
    my $snk_dir = $_[2];
    my $snk_len = $_[3];
    my $src_dir = $_[4];
    my $src_len = $_[5];

#
#  We encode the name in the form root_snk_src
#  where snk and src might be xp3, xm3, N for derivative in
#  plus and minus, or now derivative term
    
    if($snk_len == 0){
	$snk_label="nnn";
    }
    else{
	if($snk_dir == 0){
	    $snk_label = "x";
	}
	elsif($snk_dir == 1){
	    $snk_label = "y";
	}
	else{
	    $snk_label = "z";
	}
	
	if($snk_len < 0){
	    $pm = "m";
	    $snk_len_mod = $snk_len;
	    $snk_len_mod =~ s/-//;
	}
	else{
	    $pm = "p";
	    $snk_len_mod = $snk_len;
	}

	$snk_label = $snk_label . $pm; # append p or m
	$snk_label = $snk_label . $snk_len_mod;
    }

    if($src_len == 0){
	$src_label="nnn";
    }
    else{
	if($src_dir == 0){
	    $src_label = "x";
	}
	elsif($src_dir == 1){
	    $src_label = "y";
	}
	else{
	    $src_label = "z";
	}
	
	if($src_len < 0){
	    $pm = "m";
	    $src_len_mod = $src_len;
	    $src_len_mod =~ s/-//;
	}
	else{
	    $pm = "p";
	    $src_len_mod = $src_len;
	}

	$src_label = $src_label . $pm; # append p or m
	$src_label = $src_label . $src_len_mod;
    }



    return $prop_root."__".$flav."_".$snk_label."_".$src_label;
}
	    

