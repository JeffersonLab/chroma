#!/usr/bin/perl -w

while(<>){
    if(/^WFLOW/){
	if(!/time/){
	    chop ;
	    ($tag,$t,$ss,$st)=split ;
	    push @time,$t ;
	    push @act,($ss+$st)/2.0;
	    $tag="";
	}
    }
}
$eps = $time[1]-$time[0] ;
 
foreach $i (1..$#time-1){
    #$w = $time[$i]*0.5*($time[$i+1]**2*$act[$i+1] - $time[$i-1]**2*$act[$i-1] )/$eps;
    $w = $time[$i]*(2.0*$time[$i]*$act[$i] + $time[$i]**2*0.5*($act[$i+1]-$act[$i-1])/$eps);
    #push @w0,$w;
    print "$time[$i] $w\n";
}

