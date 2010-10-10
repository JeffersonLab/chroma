#!/usr/bin/perl
#
#  This little script parses Colin's group-theory input files and prints out
#  the independent elements
#
use IO::Handle;

#
# First verify the command-line arguments
#

die "Usage: parse_diquark.pl <qqq_props> <first quark> <second quark>" unless $#ARGV eq 2;

#
# Check the file exists

for($ctr = 0; $ctr <= $#ARGV - 2; $ctr++){
die "File $ARGV[$ctr] does not exist\n" unless -f $ARGV[0];
$qqq_file[$ctr] = $ARGV[$ctr];
}

$first=$ARGV[$#ARGV-1];
$second=$ARGV[$#ARGV];


open(TMPQUARK, " > tmpquark");

for($file_ctr = 0; $file_ctr <= $#ARGV - 2; $file_ctr++){
    open(QQQPROP, "< $qqq_file[$file_ctr]");
    $ctr = 0;

    $_=<QQQPROP>;

    chomp;

    ($number_qqq, $flav[0], $flav[1], $flav[2]) = split;

    print "Number of qqqs is $number_qqq\n";

#
# Now split up the file

    while(<QQQPROP>){
	
	chomp;
	($offset, $sink[0], $sink[1], $sink[2], 
	 $len_snk, $src[0], $src[1], $src[2], $len_src) = split;

#
# We now want to write out the information that specifies the diquark

	print TMPQUARK "$flav[$first] $flav[$second] $sink[$first] $sink[$second] $len_snk $src[$first] $src[$second] $len_src\n";
    }
}
close(TMPQUARK);

system("sort tmpquark | uniq > diquark");

