#!/usr/bin/perl
# This program creates the data for the chroma xml input file
#

use lib "/home/dgr/lib/";	# This needs to be changed
use Parse_groups;

#
#
#  First check the command line arguments are correct

die "Usage parse_baryons.pl <qqq prop file> <displacement table> <src_template> <prop_template> <ud_mass> <s_mass> <sink_template>
<prop_root> <config_number> <qqq_template> <cfg_template>" unless $#ARGV eq 10;

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
$qqq_template = $ARGV[9];
$cfg_template = $ARGV[10];

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

open(PROP_LIST, "> prop_tmp");	# List of propagators
open(QQQ_LIST, "> qqq_tmp");	# List of qqq files
open(QQQPROP, "< $qqq_prop");

# Read the first line to give the number of qqq propagators

while(<QQQPROP>){

    chomp;			# Remove extraneous newlines etc

    $len = split;		# Number of arguments in line

    if($len == 4){
#
#  This is the line that gives the flavour, and number of correlators

	($number_qqq, $flav[0], $flav[1], $flav[2]) = split;
	print "Starting new isospin: ";
	print "Flavours are $flav[0], $flav[1], $flav[2]\n";
	print "Number of qqqs is $number_qqq\n\n";
	$ctr = 0;
    }
    else{
#
#  We now have to parse the lines of the file, and convert the various
#  numbers into quark propagator descriptions
#
	$qqq_info = $_;
	$ctr = $ctr + 1;

	($prop[0], $prop[1], $prop[2]) = 
	    Parse_groups::extract_q($qqq_info, @flav, @disp_len);

	print QQQ_LIST "$prop[0]:$prop[1]:$prop[2]\n";	# Write the QQQ list out

#
#  We now write out the various propagator specs, and find the unique
#  list of propagators we need - we assume the root uniquely specifies
#  the propagator properties, including the corresponding template file
	
	for($prop_ctr = 0; $prop_ctr < 3; $prop_ctr++){
	    print PROP_LIST 
		"$prop[$prop_ctr]\n";
	}
    }
}
#
#  We now close the propagator and QQQ files
close(PROP_LIST);
close(QQQ_LIST);
#
#  We now have to pull out the unique lines in the file using "uniq" 
#  and "sort", nifty little unix commands
#

system("sort prop_tmp | uniq > prop_list");
system("sort qqq_tmp | uniq > qqq_list");

#
# Next we pull out the source information
open(PROP_LIST, "< prop_list");
open(SRC_LIST, "> src_list_tmp");

while(<PROP_LIST>){
    chomp;			# truncate the line
    $_ =~ s/^\s+//;		# remove leading spaces
    ($flav, $snk, $src, $snk_len, $src_len) 
	= split(/\s+/,$_);	# Split the line
    print SRC_LIST "$flav $src $src_len\n";
}

#
#  Close the files, and remove duplicate entries
close(SRC_LIST);
close(PROP_LIST);

system("sort src_list_tmp | uniq > src_list; rm src_list_tmp");

#
#  prop_list contains the unique list of smeared-sink propagators
#  src_list contans the unique list of smear-source propagators
#

#
#  The next job is to compute them within chroma, by producing the
#  appropriate chroma input files

#
# Begin by writing the chroma header

$chroma_ini = "DATA";
open(OUT,">$chroma_ini");
print OUT '<?xml version="1.0"?>';
print OUT "\n";
print OUT "<chroma>\n";
print OUT "<annotation>\n";
print OUT "QQQ inputs\n";
print OUT "</annotation>\n";
print OUT "<Param>\n";
print OUT "  <InlineMeasurements>\n";

#
#  We now read the file, pulling out the source and propagator information
#  and appending to the input file
#

open(SRC_LIST,"< src_list");

$src_ctr=0;			# Number of sources we need to compute
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

    $propname = Parse_groups::make_prop_name($prop_root."_PS", $flav, 
					     0, 0, $src_out, $src_len_out);
    $propname = $propname . ".cfg$cfg";

#
#  First write the src information to the template file
    open(IN,"< $src_template");	# Open the template file for read
    while(<IN>){		# Loop through the file, setting source
	s/_DISP_LENGTH/$src_len_out/; # Replace the source length
	s/_DISP_DIR/$src_out/;	# Replace source direction
	s/_SOURCE_NAME/source.$src_ctr/;
	print OUT $_;
    }
    close(IN);

#  Now write input to the chroma file
#

    open(IN,"< $prop_template"); # Open the template file for read
    while(<IN>){
	s/_MASS/$mass/;
	s/_SOURCE_NAME/source.$src_ctr/;
	s/_PROP_NAME/$propname/;
	print OUT $_;
    }
    close(IN);

    $src_ctr++;			# Update the source counter
}

#
#  Now we go through and SMEAR the propagators

$sink_ctr=0;
open(PROP_LIST, "< prop_list");
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

    open(IN, "< $snk_template");

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

#
#  Finally, we need to write out the qqq files
#  In the clear light of day, this can be rationalized somewhat
#  Basically, this repeats the start of the code, but without the uniq construction
#

open(QQQ_LIST, "< qqq_list");

while(<QQQ_LIST>){
    chomp;
    @prop = split(/:+/);

#
#  We now form the propagator from the various information above
#

    $qqq_name = "qqq_".$prop_root."_SS";

    for($prop_ctr = 0; $prop_ctr < 3; $prop_ctr++){
	print "prop_ctr is $prop_ctr: prop is $prop[$prop_ctr]\n";
#
#  Get the source and sink information
	($flavour,$snk, $src, $snk_len, $src_len) = 
	    split(/\s+/,$prop[$prop_ctr]);
#
#  Convert to Chroma format
	($snk_out, $snk_len_out) = Parse_groups::cmu_to_chroma($snk, $snk_len);
	($src_out, $src_len_out) = Parse_groups::cmu_to_chroma($src, $src_len);

#
#  Form the propagator names, with the configuration extension
	$propname[$prop_ctr] =
	    Parse_groups::make_prop_name($prop_root ."_SS", $flavour,
					 $snk_out, $snk_len_out,
					 $src_out, $src_len_out);
	$propname[$prop_ctr] = $propname[$prop_ctr] . ".cfg$cfg";
#
#  Append the source and sink information to the qqq name
	$qqq_name = Parse_groups::make_prop_name($qqq_name, $flavour,
						 $snk_out, $snk_len_out,
						 $src_out, $src_len_out);

    }

    $qqq_name =	$qqq_name . ".cfg$cfg";	# Append configuration number


    open(IN,"<$qqq_template");	# Open the template file for read
    while (<IN>){		# Scan through template file
	s/_PROPAGATOR1/$propname[0]/;
	s/_PROPAGATOR2/$propname[1]/;
	s/_PROPAGATOR3/$propname[2]/;
	s/_QQQNAME/$qqq_name/;

	print OUT $_;
    }
    close(IN);
}
close(QQQ_LIST);


#
#  Finally, the config file
open(IN, "< $cfg_template");
while (<IN>){
    print OUT $_;
}

print OUT "</chroma>\n";
close(OUT);
