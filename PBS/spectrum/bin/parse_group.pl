#!/usr/bin/perl
#
#  This little script parses Colin's group-theory input files and prints out
#  the independent elements
#
use IO::Handle;

#
# First verify the command-line arguments
#

die "Usage: parse_groups.pl <group file1> <group file 2>..." unless $#ARGV ge 0;

#
# Check the file exists

for($ctr = 0; $ctr <= $#ARGV; $ctr++){
die "File $ARGV[$ctr] does not exist\n" unless -f $ARGV[0];
$group_file[$ctr] = $ARGV[$ctr];
}

#
#  The first step is t
$max_qqq=-1;			# the maximum qqq

for($file_ctr = 0; $file_ctr <= $#ARGV; $file_ctr++){
    open(GROUP, "< $group_file[$file_ctr]");


#
#  Reset the counters

$line_ctr = 0;
$row = 0;
$col = 0;

$number_terms=0;
$ctr = 0;

$_=<GROUP>;

chomp;

($number_ops) = split;

print "Number of operators is $number_ops\n";

#
# We now loop over the row and column

while($col < $number_ops){
    while($row < $number_ops){
	$ctr = 0;		# set the term counter to zero
	$_ = <GROUP>;
	chomp;
	($i, $j, $number_terms) = split;
	
#
#  Now loop over all the terms
	while($ctr < $number_terms){
	    $_= <GROUP>;
	    chomp;
	    ($prop, $dirac, $waste) = split;
	    if($prop > $max_qqq){
		$max_qqq = $prop;
		$fh[$prop] = IO::Handle->new();
		open($fh[$prop], "> dirac.$prop");
	    }

	    print { $fh[$prop] }  "$prop:$dirac\n";
		 
#	    print "prop is $prop, dirac is $dirac\n";
	    $ctr = $ctr + 1;
	}
	$ctr = 0;
	$row = $row + 1;
    }
    $col = $col + 1;
    $row = 0;
}
close(GROUP);
}

#
#  Now remove the duplicated entries from the file, and sort
print "About to uniqify\n";
for($prop = 0; $prop <= $max_qqq; $prop++){
    close($fh[$prop]);
    print "About to sort\n";
    system("sort dirac.$prop | uniq > uniqdirac.$prop");
}
