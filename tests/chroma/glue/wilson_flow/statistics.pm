# This package contains simple statistics subroutines 
#
# October 1998
#
# Kostas Orginos
#
#
package statistics ;
#require Exporter ;
#@ISA = qw(Exporter);
#@EXPORT = qw( mean );

$VERSION=1.0 ;

sub mean{ # returns the mean 
    my ($mean, $count, $v)  ;

    $mean = 0 ;
    $count = 0 ;
    foreach $v ( @_ ){
	$mean += $v ;
	$count++ ;
    }
    $mean /= $count ; 
    return ($mean) ;
}

sub var{ # returns the variance  
    my ($mean, $var,$vm,$vv,$count) ;

    $var = 0 ;
    $vm = 0 ;
    $count = 0 ;
    $mean = mean(@_) ;
    foreach $v ( @_ ){
	$vv = $v - $mean ;
	$vm += $vv ;
	$var += $vv*$vv ;
	$count++ ;
    }
#see page 607 numerical recipes
    $var = ($var - $vm*$vm/$count)/($count-1) ; 
    return ($var) ;
}

sub std{ # returns the standard deviation 
    return (sqrt(var(@_))) ;
}

sub err{ # returns the standard deviation of the mean (error)  
    return (sqrt(var(@_)/@_)) ;
}

sub kurt{ # returns the kurtosis  
    my ($mean, $std, $kurt, $vv, $v, $count) ;
    
    $count = 0 ;
    $kurt = 0 ;
    $mean = mean(@_) ;
    $std = std(@_) ;
    foreach $v ( @_ ){
	$vv = ($v - $mean)/$std ;
	$kurt += $vv**4 ;
	$count++ ;
    } 
    $kurt /= $count ;  
    return ($kurt-3) ;
}

sub jackknife{ # does the jackkife and returns the  jackknife mean and error
    my( $mean,$err, $i, @list, @jack );
    @jack = () ;
    foreach $i ( 0..$#_ ){
	@list = @_ ;
	splice @list, $i, 1 ;
	push @jack, mean(@list) ;
    }
    ($mean, $err) = jackknife_list(@jack) ;
    return ( $mean, $err ) ;
}

sub jackknife_list{
# only calculates  the  jackknife mean and error given a jackknife list
    my( $mean,$err, $N );
    $N = $#_ ;
    $mean = mean(@_) ;
    $err  = $N*std(@_)/sqrt($N+1)  ;
    return ( $mean, $err ) ;
}

sub jack_kurt{ # returns the  jackknife kurtosis and error
    my( $mean,$err, $N, $i, @list, @jack );
    @jack = () ;
    $N = $#_ ;
    foreach $i ( 0..$N ){
	@list = @_ ;
	splice @list, $i, 1 ;
	push @jack, kurt(@list) ;
    }
    $mean = mean(@jack) ;
    $err  = $N*std(@jack)/sqrt(@jack)  ;
    return ( $mean, $err ) ;
}

sub autocorr{ # returns the autocorrelation of a list of data points
    my @auto ;
    my @sh = @_ ;
    my $N = $#sh ;
    my @c ;
    
    $auto[0] = 1.0 ;
    for $i ( 1..$N/2 ){
	pop @_ ;
	shift @sh ;
	@c = () ;
	for $j ( 0..$#sh ){
	    $c[$j] = $_[$j]*$sh[$j] ;
	} 
	$auto[$i] = (mean(@c) - mean(@_) * mean(@sh))/(std(@_) * std(@sh));
    }
    return @auto ;
}

sub corr{ # returns the correlation of two lists of data points
    my $l1 = shift ;
    my $l2 = shift ;
    my $N1 =  $#{$l1} ;
    my $N2 =  $#{$l2} ;
    
    die "ERROR! In statistics::corr: Lists do not have the same length\n"
	if ($N1!=$N2) ;
    for $i (0..$N1){
	$pr[$i] = $$l1[$i]*$$l2[$i] ;
    }
    $c = (mean(@pr) - mean(@{$l1}) * mean(@{$l2}))/(std(@{$l1}) * std(@{$l2}));
    return $c*($N1+1)/$N1 ;
}

sub c_mean{ # returns the complex mean 
    my @re = () ;
    my @im = () ;
    foreach $v ( @_ ){
	push @re, $$v[0] ;
	push @im, $$v[1] ;

    }
    return ( mean(@re), mean(@im)) ;
}

sub c_var{ # returns the variance  
    my @re = () ;
    my @im = () ;
    foreach $v ( @_ ){
	push @re, $$v[0] ;
	push @im, $$v[1] ;

    }
    return ( var(@re), var(@im)) ;
}

sub c_std{ # returns the standard deviation 
    my @re = () ;
    my @im = () ;
    foreach $v ( @_ ){
	push @re, $$v[0] ;
	push @im, $$v[1] ;

    }
    return ( std(@re), std(@im)) ;
}

sub c_err{ # returns the standard deviation of the mean (error)  
    my @re = () ;
    my @im = () ;
    foreach $v ( @_ ){
	push @re, $$v[0] ;
	push @im, $$v[1] ;

    }
    return ( err(@re), err(@im)) ;
}
