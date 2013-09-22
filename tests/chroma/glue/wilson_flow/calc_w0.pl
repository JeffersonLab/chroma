#!/usr/bin/perl 

my (@time, @act);

while(<>){
    if(/^WFLOW/){
	if(!/time/){
	    chop ;
	    my ($tag,$t,$ss,$st)=split ;
	    push @time,$t ;
	    push @act,($ss+$st)/2.0;
	    $tag="";
	}
    }
}
my $eps = $time[1]-$time[0] ;
 
foreach $i (1..$#time-1){
    my $w = $time[$i]*(2.0*$time[$i]*$act[$i] + $time[$i]**2*0.5*($act[$i+1]-$act[$i-1])/$eps);
    print "$time[$i] $w\n";
}

