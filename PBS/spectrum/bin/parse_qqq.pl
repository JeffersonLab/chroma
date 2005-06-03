#!/usr/bin/perl
# This program creates the data for the qqq input files
# 
# Note that this version now writes out a SINGLE input file qqq_ini.0
# 

use lib "/home/dgr/lib/";	# This needs to be changed
use Parse_groups;

#
#
#  First check the command line arguments are correct

die "Usage parse_qqq.pl <qqq prop file> <displacement table> <qqq_template> <prop_root> <config_number> <cfg_template>" unless $#ARGV eq 5;

#
#  Now check the file exists

die "File $ARGV[0] does not exist\n" unless -f $ARGV[0];

$qqq_prop = $ARGV[0];
$disp_file = $ARGV[1];
$qqq_template = $ARGV[2];
$prop_root = $ARGV[3];
$cfg = $ARGV[4];
$cfg_template = $ARGV[5];

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


#print "number_qqq is $number_qqq";

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
#  We now create the output file, and write the header
$out_qqq = "qqq_ini.".0;	# output qqq filename
open(OUT,">$out_qqq");

#
#  Write the header of the input file
print OUT '<?xml version="1.0"?>';
print OUT "\n";
print OUT "<qqq_w>\n";
print OUT "<annotation>\n";
print OUT "</annotation>\n";
print OUT "<Param>\n";

#
#  Loop over the various qqq propagators


for($ctr = 0; $ctr < $number_qqq; $ctr++){

    ($prop[0], $prop[1], $prop[2]) = 
	Parse_groups::extract_q($qqq_info[$ctr], @flav, @disp_len);

#
#  We now form the propagator from the various information above
#

    $qqq_name = "qqq_".$prop_root."_SS";

    for($prop_ctr = 0; $prop_ctr < 3; $prop_ctr++){
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

#
#  We now open the cfg template file, and write to the output file

open(IN, " < $cfg_template");

print OUT "</Param>\n";

while (<IN>){
    print OUT $_;
}

close(IN);

print OUT "</qqq_w>\n";
close(OUT);

print $number_qqq;
