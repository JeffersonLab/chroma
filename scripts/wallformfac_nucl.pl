#!/usr/bin/perl
#
#
# Usage
#   formfact.pl
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: formfact.pl\n" unless scalar($@ARGV) == 1;

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
	$proton_sp{$x, $y, $z} = "proton.$spext" ;
	$proton_ss{$x, $y, $z} = "proton.$ssext" ;
	$proton_sw{$x, $y, $z} = "proton.$swext" ;  printf "proton_sw = %s\n", $proton_sw{0,0,0};
	$proton_wp{$x, $y, $z} = "proton.$wpext" ;
	$proton_ws{$x, $y, $z} = "proton.$wsext" ;
	$proton_ww{$x, $y, $z} = "proton.$wwext" ;
      }
      else
      {
	$mom_name = "proton" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$proton_sp{$x, $y, $z} = $mom_name . ".$spext" ;
	$proton_ss{$x, $y, $z} = $mom_name . ".$ssext" ;
	$proton_sw{$x, $y, $z} = $mom_name . ".$swext" ;
	$proton_wp{$x, $y, $z} = $mom_name . ".$wpext" ;
	$proton_ws{$x, $y, $z} = $mom_name . ".$wsext" ;
	$proton_ww{$x, $y, $z} = $mom_name . ".$wwext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Normalizations
print "Nucleon form-factor";

# Assume zero momenta pion exist
if (-f proton.$ssext) {exit(1);}
if (-f proton.$swext) {exit(1);}

&ensbc("proton_norm=extract($proton_sw{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

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

      if (-f $proton_sp{$qx,$qy,$qz})
      {
	print "found for ", $proton_sp{$qx,$qy,$qz};

	$proton_energy{$qx, $qy, $qz} = "energy." . $proton_sp{$qx, $qy, $qz};
	&meff("$proton_energy{$qx, $qy, $qz}","$proton_sp{$qx,$qy,$qz}",$t_ins);
	($mass, $mass_err) = &calc("$proton_energy{$qx, $qy, $qz}");

	$proton_mass{$qx, $qy, $qz} = $mass;
	$proton_mass_err{$qx, $qy, $qz} = $mass_err;
	print "mass = ",$proton_mass{$qx, $qy, $qz};
	print "mass_name = ",$proton_energy{$qx, $qy, $qz};
      }
    }
  }
}

#------------------------------------------------------------------------------

# Terms needed for electric form factors
#
# NOTE:  
#   f0  is  u  quark contribution to E nucleon form-fac
#   f1  is  d  quark contribution to E nucleon form-fac
#
print "Electric";
foreach $s ('f0', 'f1')
{
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

	$qsq = $qx*$qx + $qy*$qy + $qz*$qz;

	next if ($qsq > $mom2_max);

	print "qsq = $qsq";

	if (! -f "${nam}_cur3ptfn_h1_${s}_g8_qx${qx}_qy${qy}_qz${qz}") {next;}

	printf "Found file %s\n","${nam}_cur3ptfn_h1_${s}_g8_qx${qx}_qy${qy}_qz${qz}";

	&realpart("${nam}_cur3ptfn_h1_${s}_g8_qx${qx}_qy${qy}_qz${qz}","${cur}_${s}_mu3_${qx}${qy}${qz}");

	if ($cur_cnt{$s,3,$qsq} == 0)
	{
	  &ensbc("${cur}_${s}_mu3_q${qsq} = ${cur}_${s}_mu3_${qx}${qy}${qz}");
	  $cur_cnt{$s,3,$qsq} = 1;
	  ++$tot_cnt{$s};
	}
	else
	{
	  &ensbc("${cur}_${s}_mu3_q${qsq} = ${cur}_${s}_mu3_q${qsq} + ${cur}_${s}_mu3_${qx}${qy}${qz}");
	  ++$cur_cnt{$s,3,$qsq};
	  ++$tot_cnt{$s};
	}
      }
    }
  }

  # Normalize
  print "normalize";
  foreach $i (0 .. $mom2_max)
  {
    if ($cur_cnt{$s,3,$i} > 0)
    {
      &ensbc("${cur}_${s}_mu3_q${i} = ${cur}_${s}_mu3_q${i} / $cur_cnt{$s,3,$i}");
    }
  }
}

# Terms needed for magnetic form factors
print "Magnetic";
foreach $s ('f2', 'f3')
{
  # Use  \Gamma_k = \sigma_3*(1+\gamma_4)  which corresponds to k = 2 in szin
  $k = 2;

  # Loop over spatial directions
  foreach $j (0 .. 2)
  {
    $g = $Vector[$j];

    # Another loop over spatial directions
    foreach $l (0 .. 2)
    {
      if ($eps{$j,$k,$l} == 0) {next;}

      # Construct necessary real parts
      # Average over all momenta
      foreach $qz (-$mommax_int .. $mommax_int)
      {
	foreach $qy (-$mommax_int .. $mommax_int)
	{
	  foreach $qx (-$mommax_int .. $mommax_int)
	  {
	    @q = ($qx, $qy, $qz);
	    
	    if ($q[$l] == 0) {next;}

	    $qsq = $qx*$qx + $qy*$qy + $qz*$qz;

	    if ($qsq > $mom2_max) {next;}

	    print "qsq = $qsq";

	    if (! -f "${nam}_cur3ptfn_h1_${s}_g${g}_qx${qx}_qy${qy}_qz${qz}") {next;}

	    &realpart("${nam}_cur3ptfn_h1_${s}_g${g}_qx${qx}_qy${qy}_qz${qz}","${cur}_${s}_mu${j}_${qx}${qy}${qz}");

	    $e = $eps{$j,$k,$l} / $q[$l];

	    print "qx=$qx, qy=$qy, qz=$qz, e=$e";

	    if ($cur_cnt{$s,$j,$qsq} == 0)
	    {
	      &ensbc("${cur}_${s}_mu${j}_q${qsq} = $e * ${cur}_${s}_mu${j}_${qx}${qy}${qz}");
	      $cur_cnt{$s,$j,$qsq} = 1;
	      ++$tot_cnt{$s};
	    }
	    else
	    {
	      &ensbc("${cur}_${s}_mu${j}_q${qsq} = ${cur}_${s}_mu${j}_q${qsq} + $e * ${cur}_${s}_mu${j}_${qx}${qy}${qz}");
	      ++$cur_cnt{$s,$j,$qsq};
	      ++$tot_cnt{$s};
	    }
	  }
	}
      }

      # Normalize
      print "normalize";
      foreach $i (1 .. $mom2_max)
      {
	if ($cur_cnt{$s,$j,$i} > 0)
	{
	  &ensbc("${cur}_${s}_mu${j}_q${i} = ${cur}_${s}_mu${j}_q${i} / $cur_cnt{$s,$j,$i}");
	}
      }
    }

  }
}


# Normalizations
# Assume zero momenta proton exist
if (-f proton.$ssext) {exit(1);}
if (-f proton.$swext) {exit(1);}

&ensbc("had_f1_q0_sw=2*$proton_sw{0,0,0}");
&ensbc("had_f1_q0_sp=2*$proton_sp{0,0,0}");
&ensbc("had_f1_q1_sp=2*$proton_sp{1,0,0}");
&ensbc("had_f1_q2_sp=2*$proton_sp{1,1,0}");
&ensbc("had_f1_q3_sp=2*$proton_sp{1,1,1}");

#&ensbc("proton_norm=extract(had_f1_q0_ss, $t_snk - $t_src)");
&ensbc("proton_norm=extract($proton_sw{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");


# Quark electric charges within certain systems
# The flip-flopping allows for the "u" to look like a "d" in the neutron.
$u_charge{"P"} = '(2/3)';
$d_charge{"P"} = '(-1/3)';

$u_charge{"N"} = '(-1/3)';
$d_charge{"N"} = '(2/3)';

#------------------------------------------------------------------------------
#
# Only if f0  was found can we have electric form-factors
#
printf "Found %d files of f0\n", $tot_cnt{"f0"} ;
if ($tot_cnt{"f0"} > 0)
{
#
# Baryon Electric form factors
#
  print "Computing nucleon electric form-factors";
  foreach $nuc ("P", "N")
  {
    &ensbc("${nuc}_j_mu3_000 = $norm*($u_charge{$nuc}*${cur}_f0_mu3_000 + $d_charge{$nuc}*${cur}_f1_mu3_000)");
    &ensbc("${nuc}_j_mu3_100 = $norm*($u_charge{$nuc}*${cur}_f0_mu3_100 + $d_charge{$nuc}*${cur}_f1_mu3_100)");
    &ensbc("${nuc}_j_mu3_110 = $norm*($u_charge{$nuc}*${cur}_f0_mu3_110 + $d_charge{$nuc}*${cur}_f1_mu3_110)");
    &ensbc("${nuc}_j_mu3_111 = $norm*($u_charge{$nuc}*${cur}_f0_mu3_111 + $d_charge{$nuc}*${cur}_f1_mu3_111)");

    foreach $i (0 .. $mom2_max)
    {
      &ensbc("${nuc}_j_mu3_q$i = $norm*($u_charge{$nuc}*${cur}_f0_mu3_q$i + $d_charge{$nuc}*${cur}_f1_mu3_q$i)");
    }
    
    &ensbc("${nuc}_r_mu3_000 = (${nuc}_j_mu3_000 * had_f1_q0_sp) / (proton_norm * had_f1_q0_sp)");
    &ensbc("${nuc}_r_mu3_100 = (${nuc}_j_mu3_100 * had_f1_q0_sp) / (proton_norm * had_f1_q1_sp)");
    &ensbc("${nuc}_r_mu3_110 = (${nuc}_j_mu3_110 * had_f1_q0_sp) / (proton_norm * had_f1_q2_sp)");
    &ensbc("${nuc}_r_mu3_111 = (${nuc}_j_mu3_111 * had_f1_q0_sp) / (proton_norm * had_f1_q3_sp)");
    
    foreach $i (0 .. $mom2_max)
    {
      &ensbc("${nuc}_r_mu3_q${i} = (${nuc}_j_mu3_q${i} * had_f1_q0_sp) / (proton_norm * had_f1_q${i}_sp)");
    }

    $t_ext = $t_snk - $t_src + 1;

#   system("calc ${nuc}_r_mu3_000 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu3_000.png");
#   system("calc ${nuc}_r_mu3_100 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu3_100.png");
#   system("calc ${nuc}_r_mu3_110 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu3_110.png");
#   system("calc ${nuc}_r_mu3_111 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu3_111.png");
    system("(echo '#e c \cr'; calc ${nuc}_r_mu3_000 | head -$t_ext) > ${nuc}_r_mu3_000.ax");
    system("(echo '#e c \cr'; calc ${nuc}_r_mu3_100 | head -$t_ext) > ${nuc}_r_mu3_100.ax");
    system("(echo '#e c \cr'; calc ${nuc}_r_mu3_110 | head -$t_ext) > ${nuc}_r_mu3_110.ax");
    system("(echo '#e c \cr'; calc ${nuc}_r_mu3_111 | head -$t_ext) > ${nuc}_r_mu3_111.ax");

    $t_ext_m1 = $t_ext - 1;

    foreach $i (0 .. $mom2_max)
    {
      system("(echo '#e c \cr'; calc ${nuc}_r_mu3_q${i} | head -$t_ext) > ${nuc}_r_mu3_q${i}.ax");
      system("(echo '#e c \cr'; calcbc \"${nuc}_r_mu3_q${i} / P_r_mu3_q0\" | head -$t_ext_m1) > ${nuc}_r_mu3_q${i}_norm.ax");
    }

    if (0)
    {
      # Extract Z_V and other ratios by simple average
      $ti = 2;
      $tf = 13;

      $ttot = $tf - $ti + 1;

      foreach $i (0 .. $mom2_max)
      {
	$t = $ti;
	
	&ensbc("${nuc}_r_mu3_q${i}_avg = extract(${nuc}_r_mu3_q${i},$t)");
	
	while ($t < $tf)
	{
	  ++$t;
	  &ensbc("${nuc}_r_mu3_q${i}_avg = ${nuc}_r_mu3_q${i}_avg + extract(${nuc}_r_mu3_q${i},$t)");
	}

	# The first average - Z_V - is the norm. for the non-zero mom. guys
	&ensbc("${nuc}_r_mu3_q${i}_avg = ${nuc}_r_mu3_q${i}_avg / $ttot");
      }
    }
  }
}


#------------------------------------------------------------------------------
#
# Only if f2  was found can we have magnetic form-factors
#
printf "Found %d files of f2\n", $tot_cnt{"f2"} ;
if ($tot_cnt{"f2"} > 0)
{
  #
  # Magnetic form factors
  #
  print "Computing nucleon magnetic form-factors";
  foreach $nuc ("P", "N")
  {
    # Use  \Gamma_k = \sigma_3*(1+\gamma_4)  which corresponds to k = 2 in szin
    $k = 2;
    
    # Loop over spatial directions
    foreach $j (0 .. 2)
    {
      $g = $Vector[$j];
      
      if ($j == $k) {next;}

      &ensbc("${nuc}_j_mu${j}_100 = $norm*($u_charge{$nuc}*${cur}_f2_mu${j}_100 + $d_charge{$nuc}*${cur}_f3_mu${j}_100)");
      &ensbc("${nuc}_j_mu${j}_110 = $norm*($u_charge{$nuc}*${cur}_f2_mu${j}_110 + $d_charge{$nuc}*${cur}_f3_mu${j}_110)");
      &ensbc("${nuc}_j_mu${j}_111 = $norm*($u_charge{$nuc}*${cur}_f2_mu${j}_111 + $d_charge{$nuc}*${cur}_f3_mu${j}_111)");
    
      foreach $i (1 .. $mom2_max)
      {
	&ensbc("${nuc}_j_mu${j}_q$i = $norm*($u_charge{$nuc}*${cur}_f2_mu${j}_q$i + $d_charge{$nuc}*${cur}_f3_mu${j}_q$i)");
      }

      &ensbc("${nuc}_r_mu${j}_100 = (${nuc}_j_mu${j}_100 * had_f1_q0_sp) / (proton_norm * had_f1_q1_sp)");
      &ensbc("${nuc}_r_mu${j}_110 = (${nuc}_j_mu${j}_110 * had_f1_q0_sp) / (proton_norm * had_f1_q2_sp)");
      &ensbc("${nuc}_r_mu${j}_111 = (${nuc}_j_mu${j}_111 * had_f1_q0_sp) / (proton_norm * had_f1_q3_sp)");

      foreach $i (1 .. $mom2_max)
      {
	&ensbc("${nuc}_r_mu${j}_q${i} = (${nuc}_j_mu${j}_q${i} * had_f1_q0_sp) / (proton_norm * had_f1_q${i}_sp)");
      }

      $t_ext = $t_snk - $t_src + 1;

#     system("calc ${nuc}_r_mu${j}_100 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu${j}_100.png");
#     system("calc ${nuc}_r_mu${j}_110 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu${j}_110.png");
#     system("calc ${nuc}_r_mu${j}_111 | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu${j}_111.png");
      system("(echo '#e c \cr'; calc ${nuc}_r_mu${j}_100 | head -$t_ext) > ${nuc}_r_mu${j}_100.ax");
      system("(echo '#e c \cr'; calc ${nuc}_r_mu${j}_110 | head -$t_ext) > ${nuc}_r_mu${j}_110.ax");
      system("(echo '#e c \cr'; calc ${nuc}_r_mu${j}_111 | head -$t_ext) > ${nuc}_r_mu${j}_111.ax");

      $t_ext_m1 = $t_ext - 1;

      foreach $i (1 .. $mom2_max)
      {
#       system("calc ${nuc}_r_mu${j}_q${i} | head -$t_ext | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu${j}_q${i}.png");
#       system("calcbc \"${nuc}_r_mu${j}_q${i} / P_r_mu3_q0\" | head -$t_ext_m1 | axis -e -c '\cr' |plot -Tpng > ${nuc}_r_mu${j}_q${i}_norm.png");
	system("(echo '#e c \cr'; calc ${nuc}_r_mu${j}_q${i} | head -$t_ext) > ${nuc}_r_mu${j}_q${i}.ax");
	system("(echo '#e c \cr'; calcbc \"${nuc}_r_mu${j}_q${i} / P_r_mu3_q0\" | head -$t_ext_m1) > ${nuc}_r_mu${j}_q${i}_norm.ax");
      }
    }
  }
}


exit(0);


sub ensbc
{
  local($line) = @_;
  
  print "ensbc: ", $line;

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
  

# This is not really correct, should use vectors of p's instead of norms
# So, only really works for say p_f=0  hence  q = p_i
sub compute_ener_qsq
{
  local($E_i,$E_f,$qsq) = @_;
  local($Qsq);

  $Qsq = ($E_f-$E_i)**2 - $qsq*(2*$pi/$L_s)**2;   # note only q=p_i
#  printf "\t%g %g\n", $E_i*200/$a, $E_f*200/$a;
  return $Qsq;
}
  
  
# This is not really correct, should use vectors of p's instead of norms
# So, only really works for say p_f=0  hence  q = p_i
sub compute_disp_qsq
{
  local($m,$qsq,$psq_i,$psq_f) = @_;
  local($Qsq);

  local($E_i) = sqrt($m**2 + $psq_i*(2*$pi/$L_s)**2);
  local($E_f) = sqrt($m**2 + $psq_f*(2*$pi/$L_s)**2);
  
  $Qsq = ($E_f-$E_i)**2 - $qsq*(2*$pi/$L_s)**2;   # note only q=p_i
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

  
