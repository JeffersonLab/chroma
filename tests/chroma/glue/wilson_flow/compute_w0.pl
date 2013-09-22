#!/usr/bin/perl -w

use File::Basename;
$basedir = dirname($0);

require "${basedir}/flow.pm" ;
require "${basedir}/statistics.pm" ;

my $c=0;
while($file = shift){
    ($rt,$rG) = get_flow($file);
    @t = @{$rt};
    @G = @{$rG};
    @W=compute_W(\@t,\@G);
    for $i (1..($#t-1)){
	$ensW[$i-1][$c] = $W[$i-1];
    }
    $c=$c+1;
}
#$Nconf=$c ;

($avW0,$erW0) = compute_w0(0.3,\@t,\@ensW);
print "W0: $avW0, $erW0\n";



