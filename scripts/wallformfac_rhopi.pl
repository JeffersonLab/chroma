#!/usr/bin/perl
#
# $Id: wallformfac_rhopi.pl,v 3.1 2007-10-31 14:14:38 edwards Exp $
#
# Usage
#   wallformfac_rhopi.pl
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

die "Put the lattice spacing \'a\' in the config.pl file\n" unless defined($a);

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
	$pion_ap{$x, $y, $z} = "pion.$apext" ;
	$pion_cp{$x, $y, $z} = "pion.$cpext" ;
	$pion_cb{$x, $y, $z} = "pion.$cbext" ;
	$pion_xp{$x, $y, $z} = "pion.$xpext" ;

	$rho_ap{$x, $y, $z} = "rho.$apext" ;
	$rho_cp{$x, $y, $z} = "rho.$cpext" ;
	$rho_cb{$x, $y, $z} = "rho.$cbext" ;
	$rho_xp{$x, $y, $z} = "rho.$xpext" ;

	$rho_x_ap{$x, $y, $z} = "rho_x.$apext" ;
	$rho_x_cp{$x, $y, $z} = "rho_x.$cpext" ;
	$rho_x_cb{$x, $y, $z} = "rho_x.$cbext" ;
	$rho_x_xp{$x, $y, $z} = "rho_x.$xpext" ;

	$rho_y_ap{$x, $y, $z} = "rho_y.$apext" ;
	$rho_y_cp{$x, $y, $z} = "rho_y.$cpext" ;
	$rho_y_cb{$x, $y, $z} = "rho_y.$cbext" ;
	$rho_y_xp{$x, $y, $z} = "rho_y.$xpext" ;

	$rho_z_ap{$x, $y, $z} = "rho_z.$apext" ;
	$rho_z_cp{$x, $y, $z} = "rho_z.$cpext" ;
	$rho_z_cb{$x, $y, $z} = "rho_z.$cbext" ;
	$rho_z_xp{$x, $y, $z} = "rho_z.$xpext" ;
      }
      else
      {
	$mom_name = "_px" . $p[0] . "_py" . $p[1] . "_pz" . $p[2];

	$pion_ap{$x, $y, $z} = "pion" . $mom_name . ".$apext" ;
	$pion_cp{$x, $y, $z} = "pion" . $mom_name . ".$cpext" ;
	$pion_cb{$x, $y, $z} = "pion" . $mom_name . ".$cbext" ;
	$pion_xp{$x, $y, $z} = "pion" . $mom_name . ".$xpext" ;

	$rho_ap{$x, $y, $z} = "rho" . $mom_name . ".$apext" ;
	$rho_cp{$x, $y, $z} = "rho" . $mom_name . ".$cpext" ;
	$rho_cb{$x, $y, $z} = "rho" . $mom_name . ".$cbext" ;
	$rho_xp{$x, $y, $z} = "rho" . $mom_name . ".$xpext" ;

	$rho_x_ap{$x, $y, $z} = "rho_x" . $mom_name . ".$apext" ;
	$rho_x_cp{$x, $y, $z} = "rho_x" . $mom_name . ".$cpext" ;
	$rho_x_cb{$x, $y, $z} = "rho_x" . $mom_name . ".$cbext" ;
	$rho_x_xp{$x, $y, $z} = "rho_x" . $mom_name . ".$xpext" ;

	$rho_y_ap{$x, $y, $z} = "rho_y" . $mom_name . ".$apext" ;
	$rho_y_cp{$x, $y, $z} = "rho_y" . $mom_name . ".$cpext" ;
	$rho_y_cb{$x, $y, $z} = "rho_y" . $mom_name . ".$cbext" ;
	$rho_y_xp{$x, $y, $z} = "rho_y" . $mom_name . ".$xpext" ;

	$rho_z_ap{$x, $y, $z} = "rho_z" . $mom_name . ".$apext" ;
	$rho_z_cp{$x, $y, $z} = "rho_z" . $mom_name . ".$cpext" ;
	$rho_z_cb{$x, $y, $z} = "rho_z" . $mom_name . ".$cbext" ;
	$rho_z_xp{$x, $y, $z} = "rho_z" . $mom_name . ".$xpext" ;
      }
    }
  }
}

#------------------------------------------------------------------------------
#
# Extract the energy of each mom. state. Use a crude exp eff. mass
#
print "Rho Electric form-factor";

# Use this as the insertion point - it is midway
$t_ins = int(($t_snk - $t_src) / 2);
print "t_ins = $t_ins";

foreach $qx ( -$mommax_int .. $mommax_int ) {
  next if ($avg_equiv_mom && ($qx < 0)) ;
  foreach $qy ( -$mommax_int .. $mommax_int ) {
    next if ($avg_equiv_mom && (($qy < 0) || ($qy > $qx))) ;
    foreach $qz ( -$mommax_int .. $mommax_int ) {
      next if ($avg_equiv_mom && (($qz < 0) || ($qz > $qy))) ;

      $mom2 = $qx*$qx + $qy*$qy + $qz*$qz ;

      next if ($mom2 > $mom2_max) ;

      # Average the rho
      if (-f $rho_x_xp{$qx,$qy,$qz})
      {
	&ensbc("$rho_xp{$qx,$qy,$qz} = ($rho_x_xp{$qx,$qy,$qz} + $rho_y_xp{$qx,$qy,$qz} + $rho_z_xp{$qx,$qy,$qz})/3");
      }
	
      if (-f $rho_x_cb{$qx,$qy,$qz})
      {
	&ensbc("$rho_cb{$qx,$qy,$qz} = ($rho_x_cb{$qx,$qy,$qz} + $rho_y_cb{$qx,$qy,$qz} + $rho_z_cb{$qx,$qy,$qz})/3");
      }
	
      # meson masses
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

	# rho mass
	$rho_energy{$qx, $qy, $qz} = "energy." . $rho_xp{$qx, $qy, $qz};
	if ($mom2 > 0)
	{
	  # Use dispersion relation
	  &meff("foo","$rho_xp{0,0,0}",$t_ins);
	  &compute_2pt_ener_disp("$rho_energy{$qx,$qy,$qz}", "foo", *q);
	}
	else
	{
	  &meff("$rho_energy{$qx, $qy, $qz}","$rho_xp{$qx,$qy,$qz}",$t_ins);
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

                printf "\nTEST: s=$s j=$j k=$k l=$l, g=$g, q=[$qx,$qy,$qz] qsq=$qsq\n";

                if ($p_i[$l] == 0) {next;}   # Check this, may not be correct with p_i fixed!!

                printf "\nNEW: s=$s j=$j k=$k l=$l, g=$g\n";

                print "q=[$q[0],$q[1],$q[2]], qsq = $qsq,  p_i=[$p_i[0],$p_i[1],$p_i[2]], p_f=[$p_f[0],$p_f[1],$p_f[2]]";

		printf "Looking for file %s\n","${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]";
		if (! -f "${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]") {next;}

		printf "Looking for file %s\n", "$pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}";
		if (! -f "$pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}") {next;}

		printf "Looking for file %s\n", "$pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]}";
		if (! -f "$pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]}") {next;}

		printf "Found file %s\n","${nam}_cur3ptfn_${s}_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]";

		$f0 = "${cur}_${h}_f0_${s}_s${k}_mu${j}_$q[0]$q[1]$q[2]";
		&imagpart("${nam}_cur3ptfn_${h}_f0_${s}_p0_snk15_g${g}_src${gk}_qx$q[0]_qy$q[1]_qz$q[2]",$f0);
		$f1 = "${cur}_${h}_f1_${s}_s${k}_mu${j}_$q[0]$q[1]$q[2]";
		&imagpart("${nam}_cur3ptfn_${h}_f1_${s}_p0_snk${gk}_g${g}_src15_qx$q[0]_qy$q[1]_qz$q[2]",$f1);

		$qsq_dim = &compute_pipf_sq($rho_mass{$cp_i[0],$cp_i[1],$cp_i[2]},*p_i,
					    $pion_mass{$cp_f[0],$cp_f[1],$cp_f[2]},*p_f);
 		$rhopi_disp = -(($fmtoGeV/$a)**2)*$qsq_dim;
#		$rhopi_disp =  (($fmtoGeV/$a)**2)*$qsq_dim;   # do not flip sign!!!!
		printf "pion mass = %g +- %g,  rho mass = %g +- %g,   qsq (via vector disp) = %g, qsq (GeV^2) = %g\n", 
		$pion_mass{$cp_f[0],$cp_f[1],$cp_f[2]}, $pion_mass_err{$cp_f[0],$cp_f[1],$cp_f[2]}, 
		$rho_mass{$cp_i[0],$cp_i[1],$cp_i[2]}, $rho_mass_err{$cp_i[0],$cp_i[1],$cp_i[2]}, 
		$qsq_dim, $rhopi_disp;


		# Normalizations
		&ensbc("pion_norm=extract($pion_cb{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");
		&ensbc("rho_norm=extract($rho_cb{$p_f[0],$p_f[1],$p_f[2]}, $t_snk - $t_src)");

		# Use some number of significant digits to uniquely identity the floating point qsq
		$qsq_int = int(10000*$rhopi_disp);

		$e = ($eps{$j,$k,$l} / $p_i[$l]);

		print "qsq_int = $qsq_int,  qx=$qx, qy=$qy, qz=$qz, e=$e";

		&ensbc("foo0 = ($e) * ($norm)*($f0 * $pion_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) * ($pion_energy{0,0,0} + $rho_energy{0,0,0}) / ($rho_ap{$cp_i[0],$cp_i[1],$cp_i[2]} * pion_norm)");
		
		&ensbc("foo1 = ($e) * ($norm)*($f1 * $rho_cp{$cp_f[0],$cp_f[1],$cp_f[2]}) * ($pion_energy{0,0,0} + $rho_energy{0,0,0}) / ($pion_ap{$cp_i[0],$cp_i[1],$cp_i[2]} * rho_norm)");
		
		&ensbc("foo_0_${s}s${k}m${j}q${qx}${qy}${qz} = foo0");
		&ensbc("foo_1_${s}s${k}m${j}q${qx}${qy}${qz} = foo1");
		
		&ensbc("foo4 = - foo0 * foo1");
		&ensbc("foo3 = sqrt(foo4) * foo4/abs(foo4)"); # make sure to sqrt a pos number - fold sign into result

		if (! defined($rhopi_cnt{$h}{$k}{$j}{$qsq_int}))
		{
		  &ensbc("${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = foo3");
		  $rhopi_cnt{$h}{$k}{$j}{$qsq_int} = 1;
		}
		else
		{
		  &ensbc("${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} + foo3");
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
	printf "norm(h=%s,k=%s,j=%s,qsq=%s)=%s\n",$h,$k,$j,$qsq_int,$rhopi_cnt{$h}{$k}{$j}{$qsq_int};
	if ($rhopi_cnt{$h}{$k}{$j}{$qsq_int} > 0)
	{
	  &ensbc("${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} = ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} / $rhopi_cnt{$h}{$k}{$j}{$qsq_int}");
	}
      }
    }
  }
}



#
# Extract Z_V for analysis
#
print "Will expect Z_V";
die "Ooops, file Z_V not found" if (! -f "Z_V");

#print "Extract Z_V";
#&ensbc("Z_V = 1 / extract(${mes}_${cur}_r_${h}_s${k}_mu${j}_q0,$t_ins)");


#
# Print Meson Electric form factors
#
print "Printing pion electric form-factors";
foreach $h (keys %rhopi_cnt)
{
  print "Printing rho->pi form-factors";
  foreach $mes ("RHO_PI")
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

	    open(FOO,"> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}.ax");
	    print FOO '#e c \cr';
	    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
	    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
	    close(FOO);
	    
	    open(FOO,"> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}_norm.ax");
	    print FOO '#e c \cr';
	    printf FOO "! a = %s fm = %g GeV^{-1}\n", $a, $fmtoGeV/$a;
	    printf FOO "! Qsq = %g GeV^{2}\n", $qsq;
	    close(FOO);
	    
	    system("calc ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} >> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}.ax");
	    system("calcbc \"${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} * Z_V\" >> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}_norm.ax");

#	    system("calcbc \"${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} / ${mes}_${cur}_r_${h}_s${k}_mu${j}_q0\" >> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}_norm.ax");

	    # Define the FF at the midpoint insertion
	    ($ff, $ff_err) = &calc("extract(${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} * Z_V,$t_ins)");
	    open(FOO,"> ${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}_ff.ax");
	    printf FOO "%g  %g %g\n", $qsq, $ff, $ff_err;
	    close(FOO);

  	    $f = "${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int}_ff";
  	    &ensbc("$f = extract(${mes}_${cur}_r_${h}_s${k}_mu${j}_q${qsq_int} * Z_V,$t_ins - 1, $t_ins + 1)");
# 	    system("~/bin/i386-linux/polyfit -t 0 -T 2 -p 0 -E eig.jknf -i $f -o foo.jknf");
# 	    &ensbc("${f}_fit = extract(foo.jknf, 0)");
# 	    ($ff, $ff_err) = &calc("${f}_fit");
# 	    open(FOO,"> ${f}_fit.ax");
# 	    printf FOO "%g  %g %g\n", $qsq, $ff, $ff_err;
# 	    close(FOO);
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
sub compute_pipf_sq
{
  local($E_i,*p_i,$E_f,*p_f) = @_;
  local($Qsq);

  local($q_sq) = 0;
  local($i);

  foreach $i (0 .. 2)
  {
    $q_sq  += ($p_i[$i] - $p_f[$i])**2;
  }

  $Qsq = ($E_i-$E_f)**2 - $q_sq*(2*$pi/$L_s)**2;
#  printf "DISP: \t%g %g  %g %g  %g\n",$E_i, $E_f, ($E_i-$E_f)**2, $q_sq*(2*$pi/$L_s)**2, $Qsq;
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

  
