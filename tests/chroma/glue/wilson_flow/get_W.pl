#!/usr/bin/perl -w

use File::Basename;
my $basedir = dirname($0);

require "${basedir}/flow.pm" ;
require "${basedir}/statistics.pm" ;
my $tt = shift ;

my $c=0;
while($file = shift){
    ($rt,$rG) = get_flow($file);
    @t = @{$rt};
    @G = @{$rG};
    @W=compute_W(\@t,\@G);
    for $i (1..($#t-1)){
	if($t[$i]==$tt){
	    print "$W[$i-1]\n";
	}
    }
}
