# This package contains ini files for Wilson flow analysis
#
# September 2013
#
# Kostas Orginos
#
#
require Exporter ;
@ISA = qw(Exporter);
@EXPORT = qw();

$VERSION=1.0 ;

sub get_flow{
    my $file = shift ;
    open(F,$file) or die "Can't open file $file\n";
    my (@gact_4i, @gact_ij) ;
    my ( @t,@G ) ;
    while(<F>){
	chop;
	if(/wflow_step/){
	    s/<wflow_step>//;
	    s/<.wflow_step>//;
	    @t = split ;
	    next;
	}
	if(/wflow_gact4i/){
	    s/<wflow_gact4i>//;
	    s/<.wflow_gact4i>//;
	    @gact_4i = split ;
	    next;
	}
	if(/wflow_gactij/){
	    s/<wflow_gactij>//;
	    s/<.wflow_gactij>//;
	    @gact_ij = split ;
	    next ;
	}
    }
    close(F);
    my $i ;
    for $i (0..$#t){
	push @G, ($gact_ij[$i]+$gact_4i[$i])/2 ;
#	print "$t[$i] $G[$i]\n";
    }

    return \@t,\@G;
}

sub compute_W{
    
    my ($rt,$rG) ;
    $rt = shift ;
    @t = @{$rt} ;
    $rG = shift ;
    @G = @{$rG} ;
    my $ww ;
    my @W ;
    my $i ;
    for $i (1..($#t-1)){
#	print "$i $t[$i] $G[$i]\n";
	$ww=2*$t[$i]**2*$G[$i] + $t[$i]**2*($G[$i+1]-$G[$i-1])/($t[$i+1]-$t[$i-1]);
	push @W,$ww ;
    }

    return @W ;
}


sub compute_w0{
    my $r ;
    my @t ;
    my @ensW ;
    my $cW0 = shift ;
    $r = shift ;
    @t = @{$r} ;
    $r = shift ;
    @ensW = @{$r} ;

    my $Ncnf = $#{$ensW[0]}+1;
#    print "$Ncnf\n";

    for $i (1..($#t-1)){
	#jackknife list
	for $c (0..$Ncnf-1){
	    @jW = @{$ensW[$i-1]};
       	    splice @jW,$c,1;
#	    print "$c $#jW\n";
	    $jensW[$i-1][$c] = statistics::mean(@jW) ;
	}
    }
    for $c (0..$Ncnf-1){
	#search 
	for $i (1..($#t-2)){
	    if($jensW[$i][$c]>$cW0){
		my $s=($jensW[$i][$c] - $jensW[$i-1][$c])/($t[$i] - $t[$i-1]);
		$w0[$c] = ($cW0 - $jensW[$i-1][$c])/$s + $t[$i-1] ; 
		last ;
	    } 
	}
    }
    
    # print the ensemble file
    if (1){
      my $ens_file = "w0.dat";
      open(FILE, "> $ens_file");
      printf FILE "%d 1 0 1 1\n", $Ncnf;
      for $c (0..$Ncnf-1){
        printf FILE "0 %s\n", $w0[$c];
      }
      close(FILE);
    }

    # averages
    my $avW0 = statistics::mean(@w0) ; 
    my $erW0 = statistics::std(@w0)*($Ncnf-1)/sqrt($Ncnf) ; 

    return $avW0,$erW0 ;
}

