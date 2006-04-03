#!/usr/bin/perl
#
# $Id: strip_spectrum_w.pl,v 3.0 2006-04-03 04:59:22 edwards Exp $
#
# SZIN-style stripper for namelist output of spectrum_w.  Namelist output
# is presented to the stripper on the STDIN.
#
# Usage
#   strip_spectrum_w.pl <tempdir>
#

eval "exec /usr/bin/perl -S $0 $*" if $running_under_some_shell ;

### $[ = 1;             # set array base to 1 (DEPRECATED)
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

die "Usage: $0 <tmpdir>\n" unless $#ARGV == 0 ;

$dir = $ARGV[0] ;

die "Can't find the directory $dir\n" unless -d $dir ;
chdir $dir || die "Can't cd to $dir: $!\n" ;

$number = 0 ;
$valid_grp_name = 0 ;

# Scan through the intro stuff to determine the number of kappas and the
# lattice size

$mom2_max = 0 ;
$FermTypeP = 0 ;
$run = "" ;
$max_retry = 100 ;

$found = 0 ;
while (($found == 0) && defined($_ = <STDIN>)) {

  chop ;          # strip record separator
  @Fld = split ;  # split(/\s+/, $_)

  if ($Fld[0] =~ /&param/) {
    $end = 0 ;
    while (($end == 0) && defined($_ = <STDIN>)) {
      chop ;
      @Fld = split ;  # split(/\s+/, $_);

      if ($Fld[0] eq 'FermTypeP')      { $FermTypeP          = $Fld[2] ; }
      if ($Fld[0] eq 'numKappa')       { $numKappa           = $Fld[2] ; }
      if ($Fld[0] =~ /Kappa/)          { $Kappa{$Fld[1]}     = $Fld[4] ; }
      if ($Fld[0] eq 'j_decay')        { $j_decay            = $Fld[2] ; }
      if ($Fld[0] eq 'Pt_src')         { $Pt_src             = $Fld[2] ; }
      if ($Fld[0] eq 'Sl_src')         { $Sl_src             = $Fld[2] ; }
      if ($Fld[0] eq 'Pt_snk')         { $Pt_snk             = $Fld[2] ; }
      if ($Fld[0] eq 'Sl_snk')         { $Sl_snk             = $Fld[2] ; }
      if ($Fld[0] eq 'Wvf_kind')       { $Wvf_kind           = $Fld[2] ; }
      if ($Fld[0] =~ /wvf_param/)      { $wvf_param{$Fld[1]} = $Fld[4] ; }
      if ($Fld[0] =~ /WvfIntPar/)      { $WvfIntPar{$Fld[1]} = $Fld[4] ; }
      if ($Fld[0] eq 'mom2_max')       { $mom2_max           = $Fld[2] ; }
      if ($Fld[0] eq 'avg_equiv_mom')  { $avg_equiv_mom      = $Fld[2] ; }

      if ($Fld[0] eq '&END' ) { $end = 1 ; }

    }  # end while (($end == 0) && defined($_ = <STDIN>))
  }  # end if ($Fld[0] =~ /&param/)

  if ($Fld[0] eq '&lattis') {
    $end = 0 ;
    while (($end == 0) && defined($_ = <STDIN>)) {
      chop ;
      @Fld = split ;  # split(/\s+/, $_);

      if ($Fld[0] =~ /nrow/) { $nrow{$Fld[1]} = $Fld[4] ; }

      if ($Fld[0] eq '&END' ) { $end = 1 ; }
    }  # end while (($end == 0) && defined($_ = <STDIN>))
    $found = 1 ;
  }  # end if ($Fld[0] eq '&lattis')
}  # end while (($found == 0) && defined($_ = <STDIN>))


$L_t = $nrow{$j_decay} ;
$L_s = $nrow{0} ;

printf "L=%d  numKappa=%d  dir=%s\n", $L_t, $numKappa, $dir ;

# Clean up some of the names
$FermTypeP =~ s/\"//g;
$Wvf_kind  =~ s/\"//g;



# HACK - in some NML the string was written out as an int...
if ($FermTypeP eq "WILSON" || $FermTypeP eq "1") {
  printf "Using Wilson Fermions\n";
  $meson_particle{'Wilson_Mesons', 'mesprop[', 0, ']'} = 'a0'    ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 1, ']'} = 'rho_x' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 2, ']'} = 'rho_y' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 3, ']'} = 'b1_z'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 4, ']'} = 'rho_z' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 5, ']'} = 'b1_y'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 6, ']'} = 'b1_x'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 7, ']'} = 'pion'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 8, ']'} = 'a0'    ;
  $meson_particle{'Wilson_Mesons', 'mesprop[', 9, ']'} = 'rho_x' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',10, ']'} = 'rho_y' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',11, ']'} = 'a1_z'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',12, ']'} = 'rho_z' ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',13, ']'} = 'a1_y'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',14, ']'} = 'a1_x'  ;
  $meson_particle{'Wilson_Mesons', 'mesprop[',15, ']'} = 'pion'  ;
} else {
  die "Unknown or unsupported FermTypeP = $FermTypeP\n" ;
}

$meson_smear_state{0}  = '_1' ;
$meson_smear_state{1}  = '_1' ;
$meson_smear_state{2}  = '_1' ;
$meson_smear_state{3}  = '_1' ;
$meson_smear_state{4}  = '_1' ;
$meson_smear_state{5}  = '_1' ;
$meson_smear_state{6}  = '_1' ;
$meson_smear_state{7}  = '_2' ;
$meson_smear_state{8}  = '_2' ;
$meson_smear_state{9}  = '_2' ;
$meson_smear_state{10} = '_2' ;
$meson_smear_state{11} = '_1' ;
$meson_smear_state{12} = '_2' ;
$meson_smear_state{13} = '_1' ;
$meson_smear_state{14} = '_1' ;
$meson_smear_state{15} = '_1' ;

$state1_col{'mesprop['} = 3 ;

$state2_col{'mesprop['} = 4 ;

$value_col{'mesprop['}  = 6 ;

#
# stringize sink_mom
#
$mommax_int = int(sqrt($mom2_max)+0.5) ;
foreach $qx ( -$mommax_int .. $mommax_int ) {
  next if ($avg_equiv_mom && ($qx < 0)) ;
  foreach $qy ( -$mommax_int .. $mommax_int ) {
    next if ($avg_equiv_mom && (($qy < 0) || ($qy > $qx))) ;
    foreach $qz ( -$mommax_int .. $mommax_int ) {
      next if ($avg_equiv_mom && (($qz < 0) || ($qz > $qy))) ;

      $mom2 = $qx*$qx + $qy*$qy + $qz*$qz ;

      next if ($mom2 > $mom2_max) ;

      if ($mom2 == 0)
      {
	$mom_name{$qx, $qy, $qz} = "" ;
      }
      else
      {
	$mom_name{$qx, $qy, $qz} = "_px" . $qx . "_py" . $qy . "_pz" . $qz ;
      }
    }
  }
}

#
# Stringize names
#
for($i = 0; $i < $numKappa; ++$i) {
  # Kappa values as an integer
  $Kappa_i = int(10000*$Kappa{$i} + 0.5) ;

  # Width names. Remove trailing zeros
  $wvf_param_i =  $wvf_param{$i} ;
  $wvf_param_i =~ s/\./p/        ;
  $wvf_param_i =~ s/[0]*$//      ;
  $wvf_param_i =~ s/p$//         ;

# HACK - in some NML the string was written as an int
  if ($Wvf_kind eq "GAUGE_INV_GAUSSIAN" || $Wvf_kind eq "1") {
    $wvf_param_i = 'G' . $wvf_param_i ;
  } else {
    die "Unknown or unsupported Wvf_kind = $Wvf_kind\n" ;
  }

  $Kappa_orig{$i} = $Kappa{$i}   ;
  $Kappa{$i}      = $Kappa_i     ;
  $wvf_param{$i}  = $wvf_param_i ;

  print $i, $Kappa{$i}, $wvf_param{$i} ;
}  # end for ($i)

#
# Create Meson source-sink smearing names
#
printf "Creating meson names\n";

for($i1=0; $i1 < $numKappa; ++$i1) {

  for($i2=0; $i2 < $numKappa; ++$i2) {

    if ( $Kappa_orig{$i1} > $Kappa_orig{$i2} ) {next} ;

    # Deal with possible Kappa combinations
    if ( $Kappa{$i1} == $Kappa{$i2} ) {
      $diag = 1 ;
      $Kappa_12 = sprintf(".D%d", $Kappa{$i1}) ;

#     die "Same smearing wvf_param in same kappa\n"
#         if ($i1 != $i2 && $wvf_param{$i1} == $wvf_param{$i2}) ;
    } else {
      $diag = 0 ;
      $Kappa_12 = sprintf(".H%d_%d", $Kappa{$i1}, $Kappa{$i2}) ;
    }

##  printf "Kappa_12 = %s\n", $Kappa_12 ;

#    # Make a shorter wvf_param name if possible
#    if ( $wvf_param{$i1} eq $wvf_param{$i2} ) {
#      $wvf_param12 = 'D' . $wvf_param{$i1} ;
#    } else {
#      $wvf_param12 = $wvf_param{$i1} . '_' . $wvf_param{$i2} ;
#    }

    $wvf_param12 = $wvf_param{$i1} . '_' . $wvf_param{$i2} ;

    foreach $m (keys %meson_smear_state) {

      $ms = $meson_smear_state{$m} ;

      if ( $Pt_src == 1 && $Pt_snk == 1 ) {
        $meson_source_sink{$m,$i1,$i2,$i1,$i2,'Point','Point'} =
          $Kappa_12 . '.P' . $ms . '.P' . $ms . '.PP' ;
      }

      if ( $diag == 1 ) {

        if ( $Pt_src == 1 && $Sl_snk == 1 ) {
          $meson_source_sink{$m,$i1,$i1,$i1,$i2,'Point','Shell'} =
            $Kappa_12 . '.P' . $ms . '.D' . $wvf_param{$i2} . $ms . '.PS' ;
        }

        if ( $Sl_src == 1 && $Pt_snk == 1) {
          $meson_source_sink{$m,$i1,$i1,$i1,$i1,'Shell','Point'} =
            $Kappa_12 . '.D' . $wvf_param{$i1} . $ms . '.P' . $ms . '.SP' ;
        }
        if ( $Sl_src == 1 && $Sl_snk == 1) {
          $meson_source_sink{$m,$i1,$i1,$i2,$i2,'Shell','Shell'} =
            $Kappa_12 . '.D' . $wvf_param{$i1} . $ms . '.D' . $wvf_param{$i2}
            . $ms . '.SS' ;
        }

      } else {  # if ($diag != 1)

        if ( $Pt_src == 1 && $Sl_snk == 1 ) {
          $meson_source_sink{$m,$i1,$i2,$i1,$i2,'Point','Shell'} =
            $Kappa_12 . '.P' . $ms . '.' . $wvf_param12 . $ms . '.PS' ;
        }

        if ( $Sl_src == 1 && $Pt_snk == 1 ) {
          $meson_source_sink{$m,$i1,$i2,$i1,$i2,'Shell','Point'} =
            $Kappa_12 . '.' . $wvf_param12 . $ms . '.P' . $ms . '.SP' ;
        }
        if ( $Sl_src == 1 && $Sl_snk == 1 ) {
          $meson_source_sink{$m,$i1,$i2,$i1,$i2,'Shell','Shell'} =
            $Kappa_12 . '.' . $wvf_param12 . $ms . '.' . $wvf_param12 . $ms
            . '.SS';
        }

      }  # end if ($diag == 1)
    }  # end foreach $m
  }  # end for $i2
}  # end for $i1

#
# Create the names of the files
#
printf "Creating meson file names\n" ;

foreach $combined (keys %meson_source_sink) {
  @cs = split($;, $combined) ;

  foreach $X (keys %meson_particle) {
    @xs = split($;, $X) ;

    if ( $xs[2] != $cs[0] ) { next ; }

    foreach $q (keys %mom_name) {

      @qs = split($;, $q) ;

      $name = ( $meson_particle{$X} . $mom_name{$q}
                . $meson_source_sink{$combined} ) ;

      $meson_particle_file{$xs[0], $xs[1], $xs[2], $xs[3], $qs[0], $qs[1],
        $qs[2], $cs[1], $cs[2], $cs[3], $cs[4], $cs[5], $cs[6]} = $name ;

      if ( -f $name ) { `rm $name` ; }
    }
  }
}

#
# Now, scan for propagators
#
printf "Searching for propagators....\n";

while (<STDIN>) {
  chop ;          # strip record separator
  @Fld = split ;  # split(/\s+/, $_) ;

  if ( $Fld[0] =~ /&End_(Staggered_h|Wilson_h|H)adron_measurements/ ) {
    $end = 0 ;
    while (($end == 0) && defined($_ = <STDIN>)) {
      chop ;
      @Fld = split ;  # split(/\s+/, $_) ;

      if ($Fld[0] eq '&END') { $end = 1 ; }
    }

    # Keep track of the number of propagators
    $number++ ;
    printf "Printed propagators for configuration $number\n" ;
  }

  if ( $Fld[0] =~ /&(Staggered_h|Wilson_h|H)adron_measurements/ ) {
    $end = 0 ;
    while (($end == 0) && defined($_ = <STDIN>)) {
      chop ;
      @Fld = split ;  # split(/\s+/, $_) ;

      if ($Fld[0] eq 'loop') { $m_loop_1 = $Fld[2] ; }
      if ($Fld[0] eq '&END') { $end = 1 ; }
    }
    $m_loop_2  = $m_loop_1  ;
    $m_smear_1 = $m_smear_2 = $m_loop_1  ;
    print $m_loop_1, 'Diag' ;
  }

  if ($Fld[0] =~ /&((Point)|(Shell))[-_]((Point)|(Shell))_/) {
    # This is the source and sink combination of the particle
    $flag =  $Fld[0] ;
    $flag =~ s/-/_/g ;
    $flag =~ s/\&//g ;

    @stuff  = split(/_/, $flag) ;
    $source = $stuff[0] ;
    $sink   = $stuff[1] ;

    $flag =~ s/^[^_]*_// ;
    $flag =~ s/^[^_]*_// ;
    $type =  $flag ;

    $end = 0 ;
    while (($end == 0) && defined($_ = <STDIN>)) {
      chop ;
      @Fld = split ;  # split(/\s+/, $_) ;

      if ( $Fld[0] eq '&END' ) {
        $end = 1 ;
      } elsif ( $Fld[0] eq 'sink_mom[' ) {
        $sink_mom[$Fld[1]] = $Fld[4] ;
      } else {
        $name = $Fld[0];

        if ($state1_col{$name} eq '') {
          $state1 = 'null' ;
          $state2 = 'null' ;
          $t      = 'null' ;
          $valuet = ''     ;
        } elsif ($state2_col{$name} eq '') {
          $state1 = 'null' ;
          $state2 = 'null' ;
          $t      = 'null' ;
          $valuet = ''     ;
        } else {
          $state1 = $Fld[$state1_col{$name}] ;
          $state2 = $Fld[$state2_col{$name}] ;
          $col    = $value_col{$name}        ;

          # Skip if not a known file
          if ($name =~ /mesprop/) {
            # Mesons
            $file = $meson_particle_file{$type, $name, $state1, $state2,
                      $sink_mom[0], $sink_mom[1], $sink_mom[2],
                      $m_loop_1, $m_loop_2, $m_smear_1, $m_smear_2,
                      $source, $sink} ;
            if ( $file eq '' ) { next ; }
            ++$meson_particle_tag{$type, $name, $state1, $state2,
                $sink_mom[0], $sink_mom[1], $sink_mom[2],
                $m_loop_1, $m_loop_2, $m_smear_1, $m_smear_2,
                $source, $sink} ;
          }

          # Open the file
          if ( ! -f $file ) {

            $i = 1 ;
            $err = open(FILE, ">" . $file) ;
            while ($err == 0 && $i < $max_retry) {
              printf "Retrying ($i) opening file $file\n" ;
              sleep 1 ;
              ++$i ;
              $err = open(FILE, ">" . $file) ;
            }
            die "After $i tries, error opening file $file\n" if $err == 0 ;

            # Write header
            printf FILE "XXXXXX %d 0 %d 1\n", $L_t, $L_s ;

          } else {  # if (-f $file)

            $i = 1 ;
            $err = open(FILE, ">>" . $file) ;
            while ($err == 0 && $i < $max_retry) {
              printf "Retrying ($i) opening file $file\n" ;
              sleep 1 ;
              ++$i ;
              $err = open(FILE, ">>" . $file) ;
            }
            die "After $i tries, error appending file $file\n"
              if $err == 0 ;

          }  # end if ( ! -f $file )

          # Write the single propagator
          $t = $Fld[1] ;
          die "t = $t\n" if $t != 0 ;
          printf FILE "0 %s\n", $Fld[$col] ;
          for ($t=1; ($t < $L_t) && defined($_ = <STDIN>); ++$t) {
            chop ;          # strip record separator
            @Fld = split ;  # split(/\s+/, $_) ;

            printf FILE "%d %s\n", $t, $Fld[$col] ;
          }

          close(FILE) ;
        }  # end if ($state1_col)
      }  # end if ($Fld[0])
    }
  }
}

printf "Found  %d  propagators\n", $number ;

# Prepend a line containing the number of propagators found
printf "Prepending number of meson propagators found.\n" ;
foreach $X (keys %meson_particle_tag) {

  $name = $meson_particle_file{$X} ;

  if ( $name ne "" ) {
    printf "Header  %s\n", $name ;
    $i = 1 ;
    $err = open(FILE, "+<" . $name) ;
    while ($err == 0 && $i < $max_retry) {
      printf "Retrying opening file $name\n" ;
      sleep 1 ;
      ++$i ;
      $err = open(FILE, "+<" . $name) ;
    }
    die "After $i tries, error opening file $name\n" if $err == 0 ;
    printf FILE "%6d", $number ;
    close(FILE) ;
  }
# else {
#   @xs = split($;, $X) ;
#   printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s: %s\n\n",
#     $xs[0], $xs[1], $xs[2], $xs[3],
#     $xs[4], $xs[5], $xs[6], $xs[7], $xs[8], $xs[9], $valuet ;
# }

}  # end foreach $X


