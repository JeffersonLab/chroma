#!/usr/bin/perl
#
# $Id: wallformfac_rhopi.pl,v 1.2 2004-06-24 07:21:53 edwards Exp $
#
# Usage
#   wallformfac_rhopi.pl
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: formfact.pl\n" unless scalar($@ARGV) == 1;

die "config.pl does not exist\n" unless -f "config.pl";

do './config.pl';

#### Example input
#$mom2_max = 5
#@sink_mom = (0,0,0);
#$a = 0.1;
#$L_s = 8;

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

# Vector form factors
$Nd = 4;
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
	$pion_sp{$x, $y, $z} = "pion.$spext" ;
	$pion_sw{$x, $y, $z} = "pion.$swext" ;  printf "pion_sw = %s\n", $pion_sw{0,0,0};
	$pion_wp{$x, $y, $z} = "pion.$ssext" ;
	$pion_ws{$x, $y, $z} = "pion.$ssext" ;

	$rho_sp{$x, $y, $z} = "rho.$spext" ;
	$rho_sw{$x, $y, $z} = "rho.$swext" ;

	$rho_x_sp{$x, $y, $z} = "rho_x.$spext" ;
	$rho_x_sw{$x, $y, $z} = "rho_x.$swext" ;
	$rho_y_sp{$x, $y, $z} = "rho_y.$spext" ;
	$rho_y_sw{$x, $y, $z} = "rho_y.$swext" ;
	$rho_z_sp{$x, $y, $z} = "rho_z.$spext" ;
	$rho_z_sw{$x, $y, $z} = "rho_z.$swext" ;
      }
      else
      {
	$mom_name = "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];
	$pion_sp{$x, $y, $z} = $mom_name . ".$spext" ;
	$pion_sw{$x, $y, $z} = $mom_name . ".$swext" ;
	$pion_wp{$x, $y, $z} = $mom_name . ".$ssext" ;
	$pion_ws{$x, $y, $z} = $mom_name . ".$ssext" ;

	$rho_sp{$x, $y, $z} = "rho" . $mom_name . ".$spext" ;
	$rho_sw{$x, $y, $z} = "rho" . $mom_name . ".$swext" ;

	$rho_x_sp{$x, $y, $z} = "rho_x" . $mom_name . ".$spext" ;
	$rho_x_sw{$x, $y, $z} = "rho_x" . $mom_name . ".$swext" ;
	$rho_y_sp{$x, $y, $z} = "rho_y" . $mom_name . ".$spext" ;
	$rho_y_sw{$x, $y, $z} = "rho_y" . $mom_name . ".$swext" ;
	$rho_z_sp{$x, $y, $z} = "rho_z" . $mom_name . ".$spext" ;
	$rho_z_sw{$x, $y, $z} = "rho_z" . $mom_name . ".$swext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Normalizations
print "Rho Electric form-factor";

# Assume zero momenta pion exist
if (-f pion.$ssext) {exit(1);}
if (-f pion.$swext) {exit(1);}

&ensbc("pion_norm=extract($pion_sw{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");
&ensbc("rho_norm=extract($rho_sw{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

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

      if (-f $rho_x_sp{$qx,$qy,$qz})
      {
	# Average the rho
	&ensbc("$rho_sp{$qx,$qy,$qz} = ($rho_x_sp{$qx,$qy,$qz} + $rho_y_sp{$qx,$qy,$qz} + $rho_z_sp{$qx,$qy,$qz})/3");
      }
	
      if (-f $pion_sp{$qx,$qy,$qz})
      {
	print "found for ", $pion_sp{$qx,$qy,$qz};

	@q = ($qx, $qy, $qz);

	$pion_energy{$qx, $qy, $qz} = "energy." . $pion_sp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("foo","$pion_sp{0,0,0}",$t_ins);
	  &compute_2pt_ener_disp("$pion_energy{$qx,$qy,$qz}", "foo", *q);
	}
	else
	{
	  &meff("$pion_energy{$qx, $qy, $qz}","$pion_sp{$qx,$qy,$qz}",$t_ins);
	}

	($mass, $mass_err) = &calc("$pion_energy{$qx, $qy, $qz}");

	$pion_mass{$qx, $qy, $qz} = $mass;
	$pion_mass_err{$qx, $qy, $qz} = $mass_err;
	print "pion mass = ",$pion_mass{$qx, $qy, $qz};
	print "pion mass_name = ",$pion_energy{$qx, $qy, $qz};

	$rho_energy{$qx, $qy, $qz} = "energy." . $rho_sp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("foo","$rho_sp{0,0,0}",$t_ins);
	  &compute_2pt_ener_disp("$rho_energy{$qx,$qy,$qz}", "foo", *q);
	}
	else
	{
	  &meff("$rho_energy{$qx, $qy, $qz}","$rho_sp{$qx,$qy,$qz}",$t_ins);
	}

	($mass, $mass_err) = &calc("$rho_energy{$qx, $qy, $qz}");

	$rho_mass{$qx, $qy, $qz} = $mass;
	$rho_mass_err{$qx, $qy, $qz} = $mass_err;
	print "rho mass = ",$rho_mass{$qx, $qy, $qz};
	print "rho mass_name = ",$rho_energy{$qx, $qy, $qz};
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
$mes = "RHO_PI";
$k = 3;
$g = $Vector[$k];

print "Electric";
foreach $h ('RHO_PI')
{
  foreach $s ('d')
  {
    # Loop over projection directions
    foreach $k (0 .. 2)
    {
      $gk = $Vector[$k];

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

		$qsq = &compute_psq(*q);

		if ($qsq > $mom2_max) {next;}

		# Construct p_i using mom. conservation
		foreach $i (0 .. 2)
		{
		  $p_i[$i] = -$q[$i] + $p_f[$i];     # note sign convention on q
		}
		$p_i_sq = &compute_psq(*p_i);

		@cp_i = &canonical_momenta(*p_i);

                printf "\nTEST: s=$s j=$j k=$k l=$l, proj=$proj g=$g, q=[$qx,$qy,$qz] qsq_int=$qsq_int\n";


                if ($p_i[$l] == 0) {next;}

                printf "\nNEW: s=$s j=$j k=$k l=$l, proj=$proj g=$g\n";

                print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_i_sq = $p_i_sq, p_f=[$p_f[0],$p_f[1],$p_f[2]]";

		printf "Looking for file %s\n","${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]";
		if (! -f "${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]") {next;}

		printf "Looking for file %s\n", "$pion_sp{$cp_f[0],$cp_f[1],$cp_f[2]}";
		if (! -f "$pion_sp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

		printf "Looking for file %s\n", "$pion_sp{$cp_i[0],$cp_i[1],$cp_i[2]}";
		if (! -f "$pion_sp{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}

		printf "Found file %s\n","${nam}_cur3ptfn_${s}_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]";

		&realpart("${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]","${cur}_${h}_f0_${s}_s{$k}_mu${j}_$q[0]$q[1]$q[2]");
		&realpart("${nam}_cur3ptfn_${h}_f1_${s}_p0_snk${gk}_g${g}_src15_qx$q[0]_qy$q[1]_qz$q[2]","${cur}_${h}_f1_${s}_s{$k}_mu${j}_$q[0]$q[1]$q[2]");

		$qsq_dim = -&compute_disp_pipf_sq($rho_mass{0,0,0},*p_i,$pion_mass{0,0,0},*p_f);
		$pion_disp = -(($fmtoGeV/$a)**2)*$qsq_dim;
		printf "pion mass = %g +- %g,  rho mass = %g +- %g,   qsq (via vector disp) = %g, qsq (GeV^2) = %g\n", 
		$pion_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $pion_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, 
		$rho_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $rho_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, 
		$qsq_dim, $pion_disp;

		# Use some number of significant digits to uniquely identity the floating point qsq
		$qsq_int = int(10000*$pion_disp);

		$e = ($eps{$j,$k,$l} / $p_i[$l]);

		print "qsq_int = $qsq_int,  qx=$qx, qy=$qy, qz=$qz, e=$e";

		&ensbc("foo0 = ($e) * ($norm)*(${cur}_${h}_f0_${s}_s${k}_mu${j}_$q[0]$q[1]$q[2] * $pion_sp{$cp_f[0],$cp_f[1],$cp_f[2]}) * ($pion_energy{0,0,0} + $rho_energy{0,0,0}) / ($rho_sp{$cp_i[0],$cp_i[1],$cp_i[2]} * pion_norm)");
		
		&ensbc("foo1 = ($e) * ($norm)*(${cur}_${h}_f1_${s}_s${k}_mu${j}_$q[0]$q[1]$q[2] * $rho_sp{$cp_f[0],$cp_f[1],$cp_f[2]}) * ($pion_energy{0,0,0} + $rho_energy{0,0,0}) / ($pion_sp{$cp_i[0],$cp_i[1],$cp_i[2]} * rho_norm)");
		
		&ensbc("foo_0_${s}s${k}m${j}q${qx}${qy}${qz} = foo0");
		&ensbc("foo_1_${s}s${k}m${j}q${qx}${qy}${qz} = foo1");
		
		if (! defined($rhopi_cnt{$h}{$k}{$j}{$qsq_int}))
		{
		  &ensbc("RHO_PI_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = foo0 * foo1");
		  $rhopi_cnt{$h}{$k}{$j}{$qsq_int} = 1;
		}
		else
		{
		  &ensbc("RHO_PI_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = RHO_PI_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} + foo0 * foo1");
		  ++$rhopi_cnt{$h}{$k}{$j}{$qsq_int};
		}

		print "now rhopi_cnt($h,$k,$j,$qsq_int) = ", $rhopi_cnt{$h}{$k}{$j}{$qsq_int};
	      }
	    }
	  }
	}
      }
    }
  }
}


# Normalize
print "Normalize";
foreach $h (keys %rhopi_cnt)
{
  foreach $k (keys %{$rhopi_cnt{$h}})
  {
    foreach $j (keys %{$rhopi_cnt{$h}{$k}})
    {
      foreach $qsq_int (keys %{$rhopi_cnt{$h}{$k}{$j}})
      {
	printf "norm(h=%s,proj=%s,j=%s,qsq=%s)=%s\n",$h,$k,$j,$qsq_int,$rhopi_cnt{$h}{$k}{$j}{$qsq_int};
	if ($rhopi_cnt{$h}{$k}{$j}{$qsq_int} > 0)
	{
	  # Correct for double counting by multiplying by 2
	  &ensbc("P_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = 2* P_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} / $rhopi_cnt{$h}{$k}{$j}{$qsq_int}");
	  &ensbc("N_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = 2* N_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} / $rhopi_cnt{$h}{$k}{$j}{$qsq_int}");
	}
      }
    }
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
foreach $h (keys %rhopi_cnt)
{
  print "Printing baryon electric and magnetic form-factors";
  foreach $nuc ("P", "N")
  {
    $t_ext = $t_snk - $t_src + 1;
    $t_ext_m1 = $t_ext - 1;

    foreach $k (keys %{$rhopi_cnt{$h}})
    {
      foreach $j (keys %{$rhopi_cnt{$h}{$k}})
      {
	foreach $qsq_int (keys %{$rhopi_cnt{$h}{$k}{$j}})
	{
	  print "keys: qsq=",$qsq_int, " cnt=",$rhopi_cnt{$h}{$k}{$j}{$qsq_int};
	  if ($rhopi_cnt{$h}{$k}{$j}{$qsq_int} > 0)
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
      }
    }
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
  local($m_i,*p_i,$m_f,*p_f) = @_;
  local($Qsq);

  local($q_sq) = 0;
  local($i);

  local($pi_sq) = &compute_psq(*p_i);
  local($pf_sq) = &compute_psq(*p_f);
  foreach $i (0 .. 2)
  {
    $q_sq  += ($p_i[$i] - $p_f[$i])**2;
  }

  local($E_i) = sqrt($m_i**2 + $pi_sq*(2*$pi/$L_s)**2);
  local($E_f) = sqrt($m_f**2 + $pf_sq*(2*$pi/$L_s)**2);
  
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

  
