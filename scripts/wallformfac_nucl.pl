#!/usr/bin/perl
#
#
# Usage
#   formfact.pl
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: $0\n" unless scalar($@ARGV) == 1;

die "config.pl does not exist\n" unless -f "config.pl";

do './config.pl';

#### Example input
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

print "source_mom=[$source_mom[0],$source_mom[1],$source_mom[2]]";
print "sink_mom=[$sink_mom[0],$sink_mom[1],$sink_mom[2]]";

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
$numGamma = 4;
$Vector[0] = 1;
$Vector[1] = 2;
$Vector[2] = 4;
$Vector[3] = 8;

$DR_sign[0] = 1;
$DR_sign[1] = -1;
$DR_sign[2] = 1;
$DR_sign[3] = 1;

$proj_sign[0] = 1;
$proj_sign[1] = -1;
$proj_sign[2] = 1;
$proj_sign[3] = -1;


# Quark electric charges within certain systems
# The flip-flopping allows for the "u" to look like a "d" in the neutron.
$P_charge{"u"} = '(2/3)';
$P_charge{"d"} = '(-1/3)';

$N_charge{"u"} = '(-1/3)';
$N_charge{"d"} = '(2/3)';

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
	$proton_ap{$x, $y, $z} = "proton.$apext" ;
	$proton_cp{$x, $y, $z} = "proton.$cpext" ;
	$proton_cb{$x, $y, $z} = "proton.$cbext" ;
	$proton_xp{$x, $y, $z} = "proton.$xpext" ;
      }
      else
      {
	$mom_name = "proton" . "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];

	$proton_ap{$x, $y, $z} = $mom_name . ".$apext" ;
	$proton_cp{$x, $y, $z} = $mom_name . ".$cpext" ;
	$proton_cb{$x, $y, $z} = $mom_name . ".$cbext" ;
	$proton_xp{$x, $y, $z} = $mom_name . ".$xpext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Normalizations
print "Nucleon form-factor";

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

      if (-f $proton_xp{$qx,$qy,$qz})
      {
	print "found for ", $proton_xp{$qx,$qy,$qz};

	@q = ($qx, $qy, $qz);

	$proton_energy{$qx, $qy, $qz} = "energy." . $proton_xp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("foo","$proton_xp{0,0,0}",$t_ins);
	  &compute_2pt_ener_disp("$proton_energy{$qx,$qy,$qz}", "foo", *q);
	}
	else
	{
	  &meff("$proton_energy{$qx, $qy, $qz}","$proton_xp{$qx,$qy,$qz}",$t_ins);
	}

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
#
# Terms needed for electric form factors
#
$proj = 0;

print "Electric";
#foreach $h ('fred')
foreach $h ('NUCL')
{
  foreach $s ('u', 'd')
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

	  print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_i_sq = $p_i_sq, p_f=[$p_f[0],$p_f[1],$p_f[2]]";

	  printf "Looking for file %s\n","${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g8_src_15_qx$q[0]_qy$q[1]_qz$q[2]";
	  if (! -f "${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g8_src15_qx$q[0]_qy$q[1]_qz$q[2]") {next;}

	  printf "Looking for file %s\n", "$proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}";
	  if (! -f "$proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

	  printf "Looking for file %s\n", "$proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]}";
	  if (! -f "$proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}

	  &realpart("${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g8_src15_qx$q[0]_qy$q[1]_qz$q[2]",
		    "${cur}_${h}_f0_${s}_p0_mu3_$q[0]$q[1]$q[2]");

	  $qsq_dim = -&compute_disp_pipf_sq($proton_mass{0,0,0},*p_i,*p_f);
	  $proton_disp = (($fmtoGeV/$a)**2)*$qsq_dim;
	  $tau = $qsq_dim / (4*$proton_mass{$cp_i[0],$cp_i[1],$cp_i[2]}
			     *$proton_mass{$cp_i[0],$cp_i[1],$cp_i[2]});
	  printf "proton mass = %g +- %g,  qsq (via vector disp) = %g, qsq (GeV^2) = %g\n", 
	  $proton_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $proton_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, $qsq_dim, 
	  $proton_disp;

	  # Normalizations
	  &ensbc("proton_norm=extract($proton_cb{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

	  # Use some number of significant digits to uniquely identity the floating point qsq
	  $qsq_int = int(10000*$proton_disp);

	  open(FOO, "> tau_${qsq_int}");
	  print FOO $tau;
	  close(FOO);

	  print "qsq_int = ", $qsq_int, "    tau = ", $tau;
	  print "proton_cnt = ", $proton_cnt{$h}{0}{3}{$qsq_int};

	  &ensbc("foo = $norm*(${cur}_${h}_f0_${s}_p0_mu3_$q[0]$q[1]$q[2] * $proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) / ($proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]} * proton_norm)");
	  
	  if (! defined($proton_cnt{$h}{0}{3}{$qsq_int}))
	  {
	    &ensbc("P_${cur}_r_${h}_p0_mu3_q${qsq_int} = $P_charge{$s}*foo");
	    &ensbc("N_${cur}_r_${h}_p0_mu3_q${qsq_int} = $N_charge{$s}*foo");

	    $proton_cnt{$h}{0}{3}{$qsq_int} = 1;
	  }
	  else
	  {
	    &ensbc("P_${cur}_r_${h}_p0_mu3_q${qsq_int} = P_${cur}_r_${h}_p0_mu3_q${qsq_int} + $P_charge{$s}*foo");
	    &ensbc("N_${cur}_r_${h}_p0_mu3_q${qsq_int} = N_${cur}_r_${h}_p0_mu3_q${qsq_int} + $N_charge{$s}*foo");

	    ++$proton_cnt{$h}{0}{3}{$qsq_int};
	  }

	  print "now proton_cnt = ", $proton_cnt{$h}{0}{3}{$qsq_int};
	}
      }
    }
  }
}


#
# Terms needed for magnetic form factors
#
print "Magnetic";
#foreach $h ('fred')
foreach $h ('NUCL')
{
  foreach $s ('u', 'd')
  {
    # Loop over projection directions
    foreach $k (0 .. 2)
    {
      # Use  \Gamma_k = \sigma_3*(1+\gamma_4)  which corresponds to k = 2 in szin
      $proj = $k + 1;

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

		$proton_disp = -(($fmtoGeV/$a)**2)*&compute_disp_pipf_sq($proton_mass{0,0,0},*p_i,*p_f);
		$qsq_int = int(10000*$proton_disp);
		printf "\nTEST: s=$s j=$j k=$k l=$l, proj=$proj g=$g, q=[$qx,$qy,$qz] qsq_int=$qsq_int\n";


		if ($p_i[$l] == 0) {next;}

		printf "\nNEW: s=$s j=$j k=$k l=$l, proj=$proj g=$g\n";

		print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_i_sq = $p_i_sq, p_f=[$p_f[0],$p_f[1],$p_f[2]]";

		printf "Looking for file %s\n","${nam}_cur3ptfn_${h}_f0_${s}_p${proj}_snk15_g${g}_src_15_qx$q[0]_qy$q[1]_qz$q[2]";
		if (! -f "${nam}_cur3ptfn_${h}_f0_${s}_p${proj}_snk15_g${g}_src15_qx$q[0]_qy$q[1]_qz$q[2]") {next;}

		printf "Looking for file %s\n", "$proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}";
		if (! -f "$proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

		printf "Looking for file %s\n", "$proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]}";
		if (! -f "$proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}


		&realpart("${nam}_cur3ptfn_${h}_f0_${s}_p${proj}_snk15_g${g}_src15_qx$q[0]_qy$q[1]_qz$q[2]",
			  "${cur}_${h}_f0_${s}_p${proj}_mu${j}_$q[0]$q[1]$q[2]");

		$proton_disp = -(($fmtoGeV/$a)**2)*&compute_disp_pipf_sq($proton_mass{0,0,0},*p_i,*p_f);
		printf "proton mass = %g +- %g,  qsq (via vector disp) = %g\n", 
		$proton_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $proton_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, $proton_disp;

		$e = (-$eps{$j,$k,$l} / $p_i[$l]);
#		$e = $proj_sign[$proj] * $DR_sign[$j] * (-$eps{$j,$k,$l} / $p_i[$l]);

		# Normalizations
		&ensbc("proton_norm=extract($proton_cb{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

		# Use some number of significant digits to uniquely identity the floating point qsq
		$qsq_int = int(10000*$proton_disp);

		print "qsq_int = $qsq_int,  qx=$qx, qy=$qy, qz=$qz, e=$e";

		&ensbc("foo = ($e) * ($norm)*(${cur}_${h}_f0_${s}_p${proj}_mu${j}_$q[0]$q[1]$q[2] * $proton_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) * ($proton_energy{$cp_i[0],$cp_i[1],$cp_i[2]} + $proton_energy{0,0,0}) / ($proton_ap{$cp_i[0],$cp_i[1],$cp_i[2]} * proton_norm)");
		
		#------------- HACK -------------
		&ensbc("foo = 2 * foo");
		#------------- END OF HACK -------------
		
		&ensbc("foo_${s}p${proj}m${j}q${qx}${qy}${qz} = foo");
		
		if (! defined($proton_cnt{$h}{$proj}{$j}{$qsq_int}))
		{
		  &ensbc("P_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = $P_charge{$s}* foo");
		  &ensbc("N_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = $N_charge{$s}* foo");

		  $proton_cnt{$h}{$proj}{$j}{$qsq_int} = 1;
		}
		else
		{
		  &ensbc("P_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = P_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} + $P_charge{$s}* foo");
		  &ensbc("N_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = N_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} + $N_charge{$s}* foo");

		  ++$proton_cnt{$h}{$proj}{$j}{$qsq_int};
		}

		print "now proton_cnt($h,$proj,$j,$qsq_int) = ", $proton_cnt{$h}{$proj}{$j}{$qsq_int};
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
foreach $h (keys %proton_cnt)
{
  foreach $proj (keys %{$proton_cnt{$h}})
  {
    foreach $j (keys %{$proton_cnt{$h}{$proj}})
    {
      foreach $qsq_int (keys %{$proton_cnt{$h}{$proj}{$j}})
      {
	printf "norm(h=%s,proj=%s,j=%s,qsq=%s)=%s\n",$h,$proj,$j,$qsq_int,$proton_cnt{$h}{$proj}{$j}{$qsq_int};
	if ($proton_cnt{$h}{$proj}{$j}{$qsq_int} > 0)
	{
	  # Correct for double counting by multiplying by 2
	  &ensbc("P_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = 2* P_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} / $proton_cnt{$h}{$proj}{$j}{$qsq_int}");
	  &ensbc("N_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} = 2* N_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} / $proton_cnt{$h}{$proj}{$j}{$qsq_int}");
	}
      }
    }
  }
}


#
# Extract Z_V for later analysis
#
print "Extract Z_V";
&ensbc("Z_V = 1 / extract(P_${cur}_r_NUCL_p0_mu3_q0 ,$t_ins)");


#------------------------------------------------------------------------------
#
#
# Print Baryon electric and magnetic form factors
#
foreach $h (keys %proton_cnt)
{
  print "Printing baryon electric and magnetic form-factors";
  foreach $nuc ("P", "N")
  {
    $t_ext = $t_snk - $t_src + 1;
    $t_ext_m1 = $t_ext - 1;

    foreach $proj (keys %{$proton_cnt{$h}})
    {
      foreach $j (keys %{$proton_cnt{$h}{$proj}})
      {
	foreach $qsq_int (keys %{$proton_cnt{$h}{$proj}{$j}})
	{
	  print "keys: qsq=",$qsq_int, " cnt=",$proton_cnt{$h}{$proj}{$j}{$qsq_int};
	  if ($proton_cnt{$h}{$proj}{$j}{$qsq_int} > 0)
	  {
	    $qsq = $qsq_int / 10000;

	    print "qsq_int=", $qsq_int;

	    open(FOO,"> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}.ax");
	    print FOO '#e c \cr';
	    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
	    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
	    close(FOO);
	    
	    open(FOO,"> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}_norm.ax");
	    print FOO '#e c \cr';
	    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
	    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
	    close(FOO);
	    
	    system("calc ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} >> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}.ax");
	    system("calcbc \"${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} * Z_V\" >> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}_norm.ax");

#	    system("calcbc \"${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} / ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q0\" >> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}_norm.ax");

	    # Define the FF at the midpoint insertion
	    ($ff, $ff_err) = &calc("extract(${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} * Z_V,$t_ins)");
	    open(FOO,"> ${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}_ff.ax");
	    printf FOO "%g  %g %g\n", $qsq, $ff, $ff_err;
	    close(FOO);

	    $f = "${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int}_ff";
	    &ensbc("$f = extract(${nuc}_${cur}_r_${h}_p${proj}_mu${j}_q${qsq_int} * Z_V,$t_ins - 1, $t_ins + 1)");
	    system("~/bin/i386-linux/polyfit -t 0 -T 2 -p 0 -E eig.jknf -i $f -o foo.jknf");
	    &ensbc("${f}_fit = extract(foo.jknf, 0)");
	    ($ff, $ff_err) = &calc("${f}_fit");
	    open(FOO,"> ${f}_fit.ax");
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

sub realpart_rev
{
  local($f1,$f2) = @_;
  
  system("/home/edwards/bin/realpart.pl < $f1 | /home/edwards/bin/reverse_prop.pl > $f2");
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

  
