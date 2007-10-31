#!/usr/bin/perl
#
#
# Usage
#   wallformfac_delta.pl
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: $0\n" unless scalar($@ARGV) == 1;

die "config.pl does not exist\n" unless -f "config.pl";

do './config.pl';

#### Example input
#$nam = 'nonlocal';
#$cur = 'n';
#$nam = 'local';
#$cur = 'l';

#$t_src = 7;
#$t_snk = 25;
#$spext = 'D600.DG2p5_1.P_1.SP';
#$ssext = 'D600.DG2p5_1.DG2p5_1.SS';
#$norm = '1';

#$t_src = 5;
#$t_snk = 20;
#$spext = 'D1480.DG3_1.P_1.SP';
#$ssext = 'D1480.DG3_1.DG3_1.SS';
#$norm = '(2*0.1480)';
#### End of example

die "Put the lattice spacing \'a\' in the config.pl file\n" unless defined(a);

$pi = 3.14159265359;
$fmtoGeV = 0.200;

$p_f_sq = 0.0;
foreach $i (0 .. 2)
{
  $p_f[$i] = $sink_mom[$i];

  $p_f_sq += $p_f[$i]**2;
}

@cp_f = &canonical_momenta(*p_f);

$mommax_int = int(sqrt($mom2_max)+0.5) ;

# Vector form factors
$numGamma = 4;
$Vector[0] = 1;
$Vector[1] = 2;
$Vector[2] = 4;
$Vector[3] = 8;

# Initialize an anti-symmetric matrix
foreach $i (0 .. 2)
{
  foreach $j (0 .. 2)
  {
    foreach $k (0 .. 2)
    {
      $eps{$i,$j,$k} = 0;
    }
  }
}

$eps{0,1,2} = $eps{1,2,0} = $eps{2,0,1} = 1;
$eps{1,0,2} = $eps{0,2,1} = $eps{2,1,0} = -1;

#------------------------------------------------------------------------------
#
# stringize sink_mom
#
$avg_equiv_mom = 1;

$mommax_int = int(sqrt($mom2_max)+0.5) ;
foreach $x (-$mommax_int .. $mommax_int)
{
  foreach $y (-$mommax_int .. $mommax_int)
  {
    foreach $z (-$mommax_int .. $mommax_int)
    {
      $p2 = $x*$x + $y*$y + $z*$z;
      next if ($p2 > $mom2_max);

      $mom[0] = $x;
      $mom[1] = $y;
      $mom[2] = $z;

      @p = &canonical_momenta(*mom);

      $mom2 = $x*$x + $y*$y + $z*$z ;

      next if ($mom2 > $mom2_max) ;

      if ($mom2 == 0)
      {
	$delta_x_sp{$x, $y, $z} = "delta_x.$spext" ;
	$delta_x_sw{$x, $y, $z} = "delta_x.$swext" ;

	$delta_y_sp{$x, $y, $z} = "delta_y.$spext" ;
	$delta_y_sw{$x, $y, $z} = "delta_y.$swext" ;

	$delta_z_sp{$x, $y, $z} = "delta_z.$spext" ;
	$delta_z_sw{$x, $y, $z} = "delta_z.$swext" ;

#	$delta_wp{$x, $y, $z} = "delta.$wpext" ;
#	$delta_ws{$x, $y, $z} = "delta.$wsext" ;
#	$delta_ww{$x, $y, $z} = "delta.$wwext" ;
      }
      else
      {
	$mom_x_name = "delta_x" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$mom_y_name = "delta_y" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$mom_z_name = "delta_z" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$delta_x_sp{$x, $y, $z} = $mom_x_name . ".$spext" ;
	$delta_z_sw{$x, $y, $z} = $mom_z_name . ".$swext" ;

#	$delta_wp{$x, $y, $z} = $mom_name . ".$wpext" ;
#	$delta_ws{$x, $y, $z} = $mom_name . ".$wsext" ;
#	$delta_ww{$x, $y, $z} = $mom_name . ".$wwext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Normalizations
print "Delta form-factor";

# Assume zero momenta delta exist
if (-f delta_x.$ssext) {exit(1);}
if (-f delta_x.$swext) {exit(1);}

# Use this as the insertion point - it is midway
$t_ins = int(($t_snk - $t_src) / 2);
print "t_ins = $t_ins";

#
# Extract the energy of each mom. state. Use a crude exp eff. mass
#
foreach $qx ( -$mommax_int .. $mommax_int ) {
  next if ($avg_equiv_mom && ($qx < 0)) ;
  foreach $qy ( -$mommax_int .. $mommax_int ) {
    next if ($avg_equiv_mom && (($qy < 0) || ($qy > $qx))) ;
    foreach $qz ( -$mommax_int .. $mommax_int ) {
      next if ($avg_equiv_mom && (($qz < 0) || ($qz > $qy))) ;

      $mom2 = $qx*$qx + $qy*$qy + $qz*$qz ;

      next if ($mom2 > $mom2_max) ;

      if (-f $delta_x_sp{$qx,$qy,$qz})
      {
	print "found for ", $delta_sp{$qx,$qy,$qz};

	@q = ($qx, $qy, $qz);

	$delta_energy{$qx, $qy, $qz} = "energy." . $delta_x_sp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("$delta_energy{0,0,0}","$delta_x_sp{0,0,0}",$t_ins);
	  ($mass_g, $mass_g_err) = &calc("$delta_energy{0,0,0}");
	  $mass = &compute_2pt_ener($mass_g, *q);
	  $mass_err = $mass_g_err;
	}
	else
	{
	  &meff("$delta_energy{$qx, $qy, $qz}","$delta_x_sp{$qx,$qy,$qz}",$t_ins);
	  ($mass, $mass_err) = &calc("$delta_energy{$qx, $qy, $qz}");
	}

	$delta_mass{$qx, $qy, $qz} = $mass;
	$delta_mass_err{$qx, $qy, $qz} = $mass_err;
	print "mass = ",$delta_mass{$qx, $qy, $qz};
	print "mass_name = ",$delta_energy{$qx, $qy, $qz};
      }
    }
  }
}

#------------------------------------------------------------------------------

# Normalizations
# Assume zero momenta delta exist
if (-f delta_x.$ssext) {exit(1);}
if (-f delta_x.$swext) {exit(1);}

&ensbc("delta_norm = extract($delta_x_sw{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");


# Quark electric charges within certain systems
# The flip-flopping allows for the "u" to look like a "d" in the neutron.
$u_charge{"D"} = '(2/3)';
$d_charge{"D"} = '(-1/3)';

#$u_charge{"N"} = '(-1/3)';
#$d_charge{"N"} = '(2/3)';

#------------------------------------------------------------------------------
#
$mes = "delta";
$u = "u";
$sigma = 1;
$tau = 1;

# Construct necessary real parts
# Average over all momenta
foreach $qz (-$mommax_int .. $mommax_int)
{
  foreach $qy (-$mommax_int .. $mommax_int)
  {
    foreach $qx (-$mommax_int .. $mommax_int)
    {
      @q = ($qx, $qy, $qz);

      $qsq = &compute_psq(*q);

      if ($qsq > $mom2_max) {next;}

      # Construct p_i using mom. conservation
      foreach $i (0 .. 2)
      {
	$p_i[$i] = -$q[$i] + $p_f[$i];     # note sign convention on q
      }
      $p_i_sq = &compute_psq(*p_i);

      @cp_i = &canonical_momenta(*p_i);

      print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_i_sq = $p_i_sq, p_f=[$p_f[0],$p_f[1],$p_f[2]]";

      $fc_u = "${nam}_cur3ptfn_DELTA_f0_u_p0_snk${sigma}_g8_src${tau}_qx$q[0]_qy$q[1]_qz$q[2]";
      printf "Looking for file %s\n",$fc_u;
      if (! -f $fc_u) {next;}

      $fc_d = "${nam}_cur3ptfn_DELTA_f0_d_p0_snk${sigma}_g8_src${tau}_qx$q[0]_qy$q[1]_qz$q[2]";
      printf "Looking for file %s\n",$fc_d;
      if (! -f $fc_d) {next;}

      printf "Looking for file %s\n", "$delta_x_sp{$cp_f[0],$cp_f[1],$cp_f[2]}";
      if (! -f "$delta_x_sp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

      printf "Looking for file %s\n", "$delta_x_sp{$cp_i[0],$cp_i[1],$cp_i[2]}";
      if (! -f "$delta_x_sp{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}

      $delta_disp = -(($fmtoGeV/$a)**2)*&compute_disp_pipf_sq($delta_mass{0,0,0},*p_i,*p_f);
      printf "delta mass = %g +- %g,  qsq (via vector disp) = %g\n", 
      $delta_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $delta_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, $delta_disp;

      $ff_u = "re_DELTA_f0_u_p0_snk${sigma}_g8_src${tau}_qx$q[0]_qy$q[1]_qz$q[2]";
      $ff_d = "re_DELTA_f0_d_p0_snk${sigma}_g8_src${tau}_qx$q[0]_qy$q[1]_qz$q[2]";

      &realpart($fc_u,$ff_u);
      &realpart($fc_d,$ff_d);

      $var = "$norm*($delta_x_sp{$cp_f[0],$cp_f[1],$cp_f[2]}) / ($delta_x_sp{$cp_i[0],$cp_i[1],$cp_i[2]} * delta_norm)";

# (2 * $delta_x_energy{$cp_f[0],$cp_f[1],$cp_f[2]} / ($delta_x_energy{$cp_i[0],$cp_i[1],$cp_i[2]} + $delta_x_energy{$cp_f[0],$cp_f[1],$cp_f[2]}))

      # Use some number of significant digits to uniquely identity the floating point qsq
      $qsq_int = int(10000*$delta_disp);

      if ($delta_cnt{$qsq_int} == 0)
      {
	$delta_cnt{$qsq_int} = 1;

	&ensbc("${mes}_r_mu3_q${qsq_int} = ((2/3)*$ff_u - (1/3)*$ff_d)*$var");
      }
      else
      {
	++$delta_cnt{$qsq_int};

	&ensbc("${mes}_r_mu3_q${qsq_int} = ${mes}_r_mu3_q${qsq_int} + ((2/3)*$ff_u - (1/3)*$ff_d)*$var");
      }
    }
  }
}

# Normalize
print "normalize";
foreach $qsq_int (keys %delta_cnt)
{
  if ($delta_cnt{$qsq_int} > 0)
  {
    &ensbc("${mes}_r_mu3_q${qsq_int} = ${mes}_r_mu3_q${qsq_int} / $delta_cnt{$qsq_int}");
  }
}



#
# Print Meson Electric form factors
#
print "Printing delta electric form-factors";
foreach $mes ("delta")
{
  $t_ext = $t_snk - $t_src + 1;
  $t_ext_m1 = $t_ext - 1;

  foreach $qsq_int (keys %delta_cnt)
  {
    $qsq = $qsq_int / 10000;

    open(FOO,"> ${mes}_r_mu3_q${qsq}.ax");
    print FOO '#e c \cr';
    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
    close(FOO);
    
    open(FOO,"> ${mes}_r_mu3_q${qsq}_norm.ax");
    print FOO '#e c \cr';
    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
    close(FOO);
    
#    system("calc ${mes}_r_mu3_q${qsq_int} | head -$t_ext >> ${mes}_r_mu3_q${qsq}.ax");
#    system("calcbc \"${mes}_r_mu3_q${qsq_int} / delta_r_mu3_q0\" | head -$t_ext_m1 > ${mes}_r_mu3_q${qsq}_norm.ax");
    system("calc ${mes}_r_mu3_q${qsq_int} >> ${mes}_r_mu3_q${qsq}.ax");
    system("calcbc \"${mes}_r_mu3_q${qsq_int} / delta_r_mu3_q0\" > ${mes}_r_mu3_q${qsq}_norm.ax");
  }
}


exit(0);



sub ensbc
{
  local($line) = @_;
  
  print "ensbc: ${line}";

  open(ENSBC, "| ensbc");
  print ENSBC $line;
  close(ENSBC);
}

sub calc
{
  local($line) = @_;
  
  $ret = `echo "calc($line)" | ensbc`;
  chop $ret;
  ($junk, $val, $err, $junk) = split(' ', $ret);

  return ($val, $err);
}

sub calc_prop
{
  local($line) = @_;
  
  $ret = `echo "calc($line)" | ensbc`;
  chop $ret;
  ($junk, $val, $err, $junk) = split(' ', $ret);

  printf "Norm of $line\n%g %g\n", $val, $err;
}

sub realpart
{
  local($f1,$f2) = @_;
  
  system("/home/edwards/bin/realpart.pl < $f1 > $f2");
}

sub imagpart
{
  local($f1,$f2) = @_;
  
  system("/home/edwards/bin/imagpart.pl < $f1 > $f2");
}

sub abs
{
  local($x) = @_;

  return ($x < 0) ? -$x : $x;
}


sub meff
{
  local($outfile,$infile,$timeslice) = @_;
  
  open(ENSBC, "| ensbc");
  print "meff: $outfile = log(extract($infile,$timeslice)/extract($infile,$timeslice+1))";
  print ENSBC "$outfile = log(extract($infile,$timeslice)/extract($infile,$timeslice+1))";
  close(ENSBC);
}


# Compute norm of  |p|^2
sub compute_psq
{
  local(*p) = @_;
  local($psq);
 
  $psq = 0;
  foreach $i (0 .. 2)
  {
    $psq += $p[$i]**2;
  }

  return $psq;
}


# Compute 2-pt energy via dispersion relation
sub compute_2pt_ener
{
  local($m,*p) = @_;

  local($p_sq) = &compute_psq(*p);
  local($E) = sqrt($m**2 + $p_sq*(2*$pi/$L_s)**2);
#  printf "\t%g %g  %g\n",$E_i*200/$a, $E_f*200/$a, $m*200/$a;
  return $E;
}
  
  
# This should be correct - uses vectors for  q = p_i - p_f
sub compute_disp_pipf_sq
{
  local($m,*p_i,*p_f) = @_;
  local($Qsq);

  local($q_sq) = 0;
  local($i);

  local($pi_sq) = &compute_psq(*p_i);
  local($pf_sq) = &compute_psq(*p_f);
  foreach $i (0 .. 2)
  {
    $q_sq  += ($p_i[$i] - $p_f[$i])**2;
  }

  local($E_i) = sqrt($m**2 + $pi_sq*(2*$pi/$L_s)**2);
  local($E_f) = sqrt($m**2 + $pf_sq*(2*$pi/$L_s)**2);
  
  $Qsq = ($E_i-$E_f)**2 - $q_sq*(2*$pi/$L_s)**2;
#  printf "\t%g %g  %g\n",$E_i*200/$a, $E_f*200/$a, $m*200/$a;
  return $Qsq;
}
  

sub reverse_sort
{
  if ($a < $b)
  {
    return 1;
  }
  elsif ($a == $b)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}
 

#
# Canonicalize momenta
#
sub canonical_momenta
{
  local(*qi) = @_;

  return sort reverse_sort (&abs($qi[0]), &abs($qi[1]), &abs($qi[2]));
}

  

  
