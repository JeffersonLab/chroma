#!/usr/bin/perl
#
# $Id: wallformfac_pion.pl,v 3.1 2007-10-31 14:14:38 edwards Exp $
#
# Usage
#   formfact.pl
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: $0\n" unless (scalar($@ARGV) == 1);

die "config.pl does not exist\n" unless -f "config.pl";

do './config.pl';

# Examples
# example of smearing names
# Mass
# $Mass = "D-7104";
# 3-pt is  AB
# wave channel
# $L = 1;
# $A = 'W';  $Aext = 'W';
# $B = 'S';  $Bext = "DG1p2";
# $C = 'S';  $Cext = "DG1p2";
# $X = 'W';  $Xext = $Aext  # used for energy only
#### End of example

die "Put the lattice spacing \'a\' in the config.pl file\n" unless defined(a);

$pi = 3.14159265359;
$fmtoGeV = 0.200;

#printf "source_mom=[%d,%d,%d]\n",$source_mom[0], $source_mom[1], $source_mom[2];
#printf "sink_mom=[%d,%d,%d]\n",$sink_mom[0], $sink_mom[1], $sink_mom[2];

print "source_mom=[$source_mom[0],$source_mom[1],$source_mom[2]]";
print "sink_mom=[$sink_mom[0],$sink_mom[1],$sink_mom[2]]";

# Check if source or sink momenta set
#die "Cannot fix source and sink momenta" if (defined(source_mom) && defined(sink_mom)); 
#die "Neither source nor sink momenta set" if ((! defined(source_mom)) && (! defined(sink_mom))); 

#$which_mom = 1 if (defined(source_mom));
#$which_mom = 2 if (defined(sink_mom));

die "momenta not set" if ($which_mom == 0);
print "which_mom = $which_mom";

if ($which_mom == 1)
{
  print "Source mom is fixed\n";
  foreach $i (0 .. 2)
  {
    $p_i[$i] = $source_mom[$i];
  }

  @cp_i = &canonical_momenta(*p_i);
}
else
{
  print "Sink mom is fixed\n";
  foreach $i (0 .. 2)
  {
    $p_f[$i] = $sink_mom[$i];
  }

  @cp_f = &canonical_momenta(*p_f);
}

# Extensions
# $abext = "${Mass}.${Aext}_${L}.${Bext}_${L}.${A}${B}";   # not used
$apext = "${Mass}.${Aext}_${L}.P_${L}.${A}P";
$cpext = "${Mass}.${Cext}_${L}.P_${L}.${C}P";
$cbext = "${Mass}.${Cext}_${L}.${Bext}_${L}.${C}${B}";
$xpext = "${Mass}.${Xext}_${L}.P_${L}.${X}P";


# Vector form factors
$Nd = 4;
$numGamma = 4;
$Vector[0] = 1;
$Vector[1] = 2;
$Vector[2] = 4;
$Vector[3] = 8;

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
	$pion_ap{$x, $y, $z} = "pion.$apext" ;
	$pion_cp{$x, $y, $z} = "pion.$cpext" ;
	$pion_cb{$x, $y, $z} = "pion.$cbext" ;
	$pion_xp{$x, $y, $z} = "pion.$xpext" ;
      }
      else
      {
	$mom_name = "pion" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$pion_ap{$x, $y, $z} = $mom_name . ".$apext" ;
	$pion_cp{$x, $y, $z} = $mom_name . ".$cpext" ;
	$pion_cb{$x, $y, $z} = $mom_name . ".$cbext" ;
	$pion_xp{$x, $y, $z} = $mom_name . ".$xpext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Normalizations
print "Pion Electric form-factor";

# Assume zero momenta pion exist
if (-f pion.$cbext) {exit(1);}

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

      print "q=[$qx,$qy,$qz], qsq = $mom2   $pion_xp{$qx,$qy,$qz}";

      if (-f $pion_xp{$qx,$qy,$qz})
      {
	print "found for ", $pion_xp{$qx,$qy,$qz};

	@q = ($qx, $qy, $qz);

	$pion_energy{$qx, $qy, $qz} = "energy." . $pion_xp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("foo","$pion_xp{0,0,0}",$t_ins);
	  &compute_2pt_ener_disp("$pion_energy{$qx,$qy,$qz}", "foo", *q);
	}
	else
	{
	  &meff("$pion_energy{$qx, $qy, $qz}","$pion_xp{$qx,$qy,$qz}",$t_ins);
	}

	($mass, $mass_err) = &calc("$pion_energy{$qx, $qy, $qz}");

	$pion_mass{$qx, $qy, $qz} = $mass;
	$pion_mass_err{$qx, $qy, $qz} = $mass_err;
	print "pion mass = ",$pion_mass{$qx, $qy, $qz};
	print "pion mass_name = ",$pion_energy{$qx, $qy, $qz};

	$mom_tag{$qz, $qy, $qz} = 1;
      }
    }
  }
}


# Terms needed for electric form factors
#
# NOTE:  
#   s0  is  u  quark contribution to E nucleon form-fac
#   s1  is  d  quark contribution to E nucleon form-fac
#   s10 is  (u or d)  quark contribution to pion form-fac
#
print "Electric";
$s = "PION_f0_d_p0";
$mes = "PION";
$k = 3;
$g = $Vector[$k];

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

      if ($which_mom == 1)
      {
        print "Reconstruct p_f";
	# Construct p_f using mom. conservation
	foreach $i (0 .. 2)
	{
	  $p_f[$i] =  $q[$i] + $p_i[$i];     # note sign convention on q
	}
	@cp_f = &canonical_momenta(*p_f);
      }
      else
      {
        print "Reconstruct p_i";
	# Construct p_i using mom. conservation
	foreach $i (0 .. 2)
	{
	  $p_i[$i] = -$q[$i] + $p_f[$i];     # note sign convention on q
	}
	@cp_i = &canonical_momenta(*p_i);
      }

      print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_f=[$p_f[0],$p_f[1],$p_f[2]]";
      print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  cp_i=[$cp_i[0],$cp_i[1],$cp_i[2]], cp_f=[$cp_f[0],$cp_f[1],$cp_f[2]]";

      printf "Looking for file %s\n","${nam}_cur3ptfn_${s}_snk15_g8_src_15_qx$q[0]_qy$q[1]_qz$q[2]";
      if (! -f "${nam}_cur3ptfn_${s}_snk15_g8_src15_qx$q[0]_qy$q[1]_qz$q[2]") {next;}

      printf "Looking for file %s\n", "$pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}";
      if (! -f "$pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

      printf "Looking for file %s\n", "$pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]}";
      if (! -f "$pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}

      printf "Found file %s\n","${nam}_cur3ptfn_${s}_snk15_g8_src15_qx$q[0]_qy$q[1]_qz$q[2]";

      $pion_disp = -(($fmtoGeV/$a)**2)*&compute_disp_pipf_sq($pion_mass{0,0,0},*p_i,*p_f);
      printf "pion_mass_i = %g +- %g,  pion_mass_f = %g +- %g,  qsq (via vector disp) = %g\n", 
      $pion_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $pion_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, 
      $pion_mass{$cp_f[0],$cp_f[1],$cp_f[2]}, $pion_mass_err{$cp_f[0],$cp_f[1],$cp_f[2]}, 
      $pion_disp;

      &realpart("${nam}_cur3ptfn_${s}_snk15_g8_src15_qx$q[0]_qy$q[1]_qz$q[2]","${cur}_${s}_mu3_$q[0]$q[1]$q[2]");
      &ensbc("pion_norm=extract($pion_cb{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

      $var = "$norm*(${cur}_${s}_mu3_$q[0]$q[1]$q[2] * $pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) * (2 * $pion_energy{$cp_f[0],$cp_f[1],$cp_f[2]} / ($pion_energy{$cp_i[0],$cp_i[1],$cp_i[2]} + $pion_energy{$cp_f[0],$cp_f[1],$cp_f[2]}))/ ($pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]} * pion_norm)";

      # hack - for p_f=0 this will average over all q momenta of the current
      $var_sum1 = "${cur}_${s}_mu3_$q[0]$q[1]$q[2]";
      $var_sum2 = "$norm*(${cur}_${s}_mu3_$q[0]$q[1]$q[2] * $pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]})";
      $var_sum3 = "$norm*(${cur}_${s}_mu3_$q[0]$q[1]$q[2] * $pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) * (2 * $pion_energy{$cp_f[0],$cp_f[1],$cp_f[2]} / ($pion_energy{$cp_i[0],$cp_i[1],$cp_i[2]} + $pion_energy{$cp_f[0],$cp_f[1],$cp_f[2]}))";


      # Use some number of significant digits to uniquely identity the floating point qsq
      $qsq_int = int(10000*$pion_disp);

      print "qsq_int = ", $qsq_int;

      if ($pion_cnt{$qsq_int} == 0)
      {
	$pion_cnt{$qsq_int} = 1;

	&ensbc("${mes}_r_mu3_q${qsq_int} = $var");
	&ensbc("sum1_q${qsq_int} = $var_sum1");
	&ensbc("sum2_q${qsq_int} = $var_sum2");
	&ensbc("sum3_q${qsq_int} = $var_sum3");
      }
      else
      {
	++$pion_cnt{$qsq_int};

	&ensbc("${mes}_r_mu3_q${qsq_int} = ${mes}_r_mu3_q${qsq_int} + $var");
	&ensbc("sum1_q${qsq_int} = sum1_q${qsq_int} + $var_sum1");
	&ensbc("sum2_q${qsq_int} = sum2_q${qsq_int} + $var_sum2");
	&ensbc("sum3_q${qsq_int} = sum3_q${qsq_int} + $var_sum3");
      }
    }
  }
}

# Normalize
print "normalize";
foreach $qsq_int (keys %pion_cnt)
{
  if ($pion_cnt{$qsq_int} > 0)
  {
    &ensbc("${mes}_r_mu3_q${qsq_int} = ${mes}_r_mu3_q${qsq_int} / $pion_cnt{$qsq_int}");
    &ensbc("sum1_q${qsq_int} = sum1_q${qsq_int} / $pion_cnt{$qsq_int}");
    &ensbc("sum2_q${qsq_int} = sum2_q${qsq_int} / $pion_cnt{$qsq_int}");
    &ensbc("sum3_q${qsq_int} = sum3_q${qsq_int} / $pion_cnt{$qsq_int}");
  }
}

#
# Extract Z_V for analysis
#
print "Extract Z_V";
&ensbc("Z_V = 1 / extract(${mes}_r_mu3_q0,$t_ins)");


#
# Print Meson Electric form factors
#
print "Printing pion electric form-factors";
foreach $mes ("PION")
{
  $t_ext = $t_snk - $t_src + 1;
  $t_ext_m1 = $t_ext - 1;

  foreach $qsq_int (keys %pion_cnt)
  {
    $qsq = $qsq_int / 10000;

    print "qsq_int=", $qsq_int;

    open(FOO,"> ${mes}_r_mu3_q${qsq_int}.ax");
    print FOO '#e c \cr';
    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
    close(FOO);
    
    open(FOO,"> ${mes}_r_mu3_q${qsq_int}_norm.ax");
    print FOO '#e c \cr';
    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
    close(FOO);
    
#    system("calc ${mes}_r_mu3_q${qsq_int} | head -$t_ext >> ${mes}_r_mu3_q${qsq_int}.ax");
#    system("calcbc \"${mes}_r_mu3_q${qsq_int} / pion_r_mu3_q0\" | head -$t_ext_m1 > ${mes}_r_mu3_q${qsq_int}_norm.ax");
    system("calc ${mes}_r_mu3_q${qsq_int} >> ${mes}_r_mu3_q${qsq_int}.ax");
    system("calcbc \"${mes}_r_mu3_q${qsq_int} / ${mes}_r_mu3_q0\" > ${mes}_r_mu3_q${qsq_int}_norm.ax");

    # Define the FF at the midpoint insertion
    ($ff, $ff_err) = &calc("extract(${mes}_r_mu3_q${qsq_int} / ${mes}_r_mu3_q0,$t_ins)");
    open(FOO,"> ${mes}_r_mu3_q${qsq_int}_ff.ax");
    printf FOO "%g  %g %g\n", $qsq, $ff, $ff_err;
    close(FOO);
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
sub compute_2pt_ener_disp
{
  local($f,$m,*p) = @_;

  local($p_sq) = &compute_psq(*p);
  &ensbc("$f = sqrt($m*$m + $p_sq*(2*$pi/$L_s)*(2*$pi/$L_s))");
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

  
