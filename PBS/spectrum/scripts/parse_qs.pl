#!/usr/bin/perl
# This program creates the data for the qqq input files
#

use lib "/home/adam/research/lhpc/chroma/PBS/spectrum/lib/";	# This needs to be changed
use Parse_groups;

#
#
#  First check the command line arguments are correct



die "Usage parse_qs.pl <qqq prop file> <displacement table> <src_template> <prop_template> <ud_mass> <s_mass>
<sink_template> <prop_root> <config_number> <cfg_template>" unless $#ARGV eq 9;

#
#  Now check the file exists

die "File $ARGV[0] does not exist\n" unless -f $ARGV[0];

$qqq_prop = $ARGV[0];
$disp_file = $ARGV[1];
$src_template = $ARGV[2];
$prop_template = $ARGV[3];
$ud_mass = $ARGV[4];
$s_mass = $ARGV[5];
$snk_template = $ARGV[6];
$prop_root = $ARGV[7];
$cfg = $ARGV[8];
$cfg_template = $ARGV[9];

#
#  The first step is to open the smearing length file, and
#  translate the smearing offsets into real smearing lengths
#

open(DISP, "< $disp_file");

$disp_number = 0;

while (<DISP>){

#
#  Extract each line, and map to the various smearing leng

    chomp;

    $disp_len[$disp_number] = $_;
    $disp_number = $disp_number + 1;
}

print "disp_len is @disp_len\n";

#
#  We now read the first line of the propagator, and compute the propagators
#  accordingly

open(QQQPROP, "< $qqq_prop");

# Read the first line to give the number of qqq propagators

$_ = <QQQPROP>; chomp;
($number_qqq, $flav[0], $flav[1], $flav[2]) = split;

print "number_qqq is $number_qqq";

#
#  We now have to parse the lines of the file, and convert the various
#  numbers into quark propagator descriptions
#

$ctr = 0;

while (<QQQPROP>){
#
#  We now need to extract the propagator information from the file
#

#  First truncate the line
    chomp;
# We now pull out the source and sink of the propagators
    $qqq_info[$ctr] = $_;
    $ctr = $ctr + 1;
}

#
#  So the qqq info is now stored in a file, and we have to extract the
#  Propagators we need

die "$qqq_prop has wrong number of lines" unless $ctr eq $number_qqq;

#
#  Loop over the various qqq propagators

open(PROP_LIST, "> prop_tmp");	# List of propagators

for($ctr = 0; $ctr < $number_qqq; $ctr++){

    ($prop[0], $prop[1], $prop[2]) = 
	Parse_groups::extract_q($qqq_info[$ctr], @flav, @disp_len);

#
#  We now write out the various propagator specs, and find the unique
#  list of propagators we need - we assume the root uniquely specifies
#  the propagator properties, including the corresponding template file

    for($prop_ctr = 0; $prop_ctr < 3; $prop_ctr++){
	print PROP_LIST 
	    "$prop[$prop_ctr]\n";
    }

#
#  Return the number of qqqs computed

#    return ($number_qqq);
}

#
#  We now close the propagator files
close(PROP_LIST);

#
#  We now have to pull out the unique lines in the file using "uniq" 
#  and "sort", nifty little unix commands
#

system("sort prop_tmp | uniq > prop_list");
#system("rm prop_tmp");

#
#  We now have a unique list of propagators that we need to generate.
#  The next job is to pull them out in terms of sources and sinks
#

#  The way to do is this first to pull out all the sources
#  Begin with the u propagators

open(PROP_LIST, "< prop_list");
open(SRC_LIST, "> src_list_tmp");

#
#  We now read in the file, pulling out the source information
#  and writing the input file for the sink

$ctr=0;
$sink_ctr=0;

while(<PROP_LIST>){
    
    chomp;			# truncate the line
    $_ =~ s/^\s+//;		# remove leading spaces
    ($flav, $snk, $src, $snk_len, $src_len) 
	= split(/\s+/,$_);	# Split the line

    ($snk_out, $snk_len_out) = Parse_groups::cmu_to_chroma($snk, $snk_len);
    ($src_out, $src_len_out) = Parse_groups::cmu_to_chroma($src, $src_len);

    # Note here that setting the sink length to zero is equivalent to 
    # not smearing

    $propname_in = Parse_groups::make_prop_name($prop_root."_PS", $flav,
				  $snk_out, 0, $src_out, $src_len_out);
    $propname_in = $propname_in . ".cfg$cfg";

    $propname_out = Parse_groups::make_prop_name($prop_root."_SS",
						 $flav, $snk_out,
						 $snk_len_out,
						 $src_out,
						 $src_len_out);
						 $propname_out =
						 $propname_out
						 . ".cfg$cfg";

    # We now see if the sink-smeared propagator exists, and if so
    # we generate a sink-smearing elem for it
    #

    if(! -e $propname_out){

	open(IN, "< $snk_template");

#
#  If this is the first call to the file, we open the file, and write
#  header information
#
	if($sink_ctr == 0){
	    open(OUT, "> sink_ini.$sink_ctr"); # Open the file
#
#  Write the header of the input file
	    print OUT '<?xml version="1.0"?>';
	    print OUT "\n";
	    print OUT "<sink_smearing>\n";
	    print OUT "<annotation>\n";
	    print OUT "</annotation>\n";
	    print OUT "<Param>\n";
	}



	while(<IN>){
	    s/_PROP_NAME_IN/$propname_in/;
	    s/_PROP_NAME_OUT/$propname_out/;
	    s/_DISP_LENGTH/$snk_len_out/; # Replace the source length
	    s/_DISP_DIR/$snk_out/;	# Replace source direction
	    print OUT $_;
	}
	close(IN);
	$sink_ctr++;		# Update the sink counter
    }
    

#  Now write out the sources

    print SRC_LIST "$flav $src $src_len\n";
    $ctr++;
}

#
#  We conclude by writing out the tail information, assuming we've written
#  an ini file
#

if( $sink_ctr != 0){

    open(IN, " < $cfg_template");

    print OUT "</Param>\n";

    while (<IN>){
	print OUT $_;
    }

    close(IN);

    print OUT "</sink_smearing>\n";
    close(OUT);
}


#
#  Close the files, and remove duplicate entries
close(SRC_LIST);
close(PROP_LIST);

system("sort src_list_tmp | uniq > src_list; rm src_list_tmp");

#
#  We now generate the source files on the basis of this file

open(SRC_LIST,"< src_list");


$ctr=0;				# Set the counter to zero
$src_ctr=0;

while(<SRC_LIST>){

    chomp;			# Truncate the line
    $_ =~ s/^\s+//;		# remove leading spaces

    ($flav, $src, $src_len) = split(/\s+/,$_);

#  There are different templates for ud and for s quarks, since the
#  masses are different

    if($flav eq "ud"){
	$mass = $ud_mass;
    }
    else{
	$mass = $s_mass;
    }

    ($src_out, $src_len_out) = Parse_groups::cmu_to_chroma($src, $src_len);

    $propname = Parse_groups::make_prop_name($prop_root."_PS", $flav, 0, 0, $src_out, $src_len_out);
    $propname = $propname . ".cfg$cfg";

    print "DEBUG propname is $propname\n";

    # We now see if this exists, and if not create the prop and src input files

    if(! -e $propname){

	open(IN,"< $prop_template");	# Open the template file for read
	open(OUT,"> prop_ini.$src_ctr");

	while(<IN>){
	    s/_MASS/$mass/;
	    s/_SOURCE_NAME/source.$src_ctr/;
	    s/_PROP_NAME/$propname/;
	    print OUT $_;
	}
	close(IN);
	close(OUT);

#
#  We no read in the template file, and write it out the appropriate
#  source file
#
	open(IN,"< $src_template");	# Open the template file for read
	open(OUT,"> src_ini.$src_ctr");
    
	while(<IN>){
	    s/_DISP_LENGTH/$src_len_out/; # Replace the source length
	    s/_DISP_DIR/$src_out/;	# Replace source direction
	    s/_SOURCE_NAME/source.$src_ctr/;
	    print OUT $_;
	}
	close(IN);
	close(OUT);
	$src_ctr++;		# Finally update the source sounter
    }

    $ctr++;
} 

$number_ud_source=$ctr;


	
