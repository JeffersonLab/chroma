#!/usr/local/bin/perl
#
# $Id: strip_bar3ptfn.pl,v 1.1 2003-05-09 23:22:11 edwards Exp $
#
# Strip out all data of bar3ptfn.m output
#
# Usage
#   strip_bar3ptfn.pl  tempdir
#

$[ = 0;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

#die "Usage: strip_bar3ptfn.pl  temp_directory\n" unless scalar($@ARGV) == 1;

$dir = $ARGV[0];

die "Can't find the directory $dir\n" unless -d $dir;
chdir $dir || die "Can't cd to $dir: $!\n";

$number = 0;
$valid_grp_name = 0;

# Scan through the intro stuff to determine the number of kappas and the
# lattice size

# Default V2 format. V3 now defines these things
$out_version = 2;
$numSeq_src = 4;
$Seq_src[0] = 0;
$Seq_src[1] = 1;
$Seq_src[2] = 2;
$Seq_src[3] = 3;
$numGamma = 4;
$Gamma_list[0] = 1;
$Gamma_list[1] = 2;
$Gamma_list[2] = 4;
$Gamma_list[3] = 8;


# other stuff
$FermTypeP = 0;
$Z3_src = 0;

$found = 0;
while ($found == 0)
{
  $_ = <STDIN>;
  chop;	# strip record separator
  @Fld = split(' ', $_);
  
  if ( $Fld[0] eq '&Output_version' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'out_version' ) {$out_version = $Fld[2];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param1' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'Nd' ) {$Nd = $Fld[2];}
      if ( $Fld[0] eq 'FermTypeP' ) {$FermTypeP = $Fld[2];}
      if ( $Fld[0] eq 'numKappa' ) {$numKappa = $Fld[2];}
      if ( $Fld[0] eq 'Kappa(' ) {$Kappa[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param2' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'FermAct' ) {$FermAct = $Fld[2];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param3' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'j_decay' ) {$j_decay = $Fld[2];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param5' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      # New format
      if ( $Fld[0] eq 'Z3_src' ) {$Z3_src = $Fld[2];}
      if ( $Fld[0] eq 'Pt_src' ) {$Pt_src = $Fld[2];}
      if ( $Fld[0] eq 'Sl_src' ) {$Sl_src = $Fld[2];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param6' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      # New format
      if ( $Fld[0] eq 'Pt_snk' ) {$Pt_snk = $Fld[2];}
      if ( $Fld[0] eq 'Sl_snk' ) {$Sl_snk = $Fld[2];}
      if ( $Fld[0] eq 't_sink' ) {$t_snk = $Fld[2];}
      if ( $Fld[0] eq 'sink_mom(' ) {$sink_mom[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }


  if ( $Fld[0] eq '&param7' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'Wvf_kind' ) {$Wvf_kind = $Fld[2];}
      if ( $Fld[0] eq 'wvf_param(' ) {$wvf_param[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&param9' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      # V3 format
      if ( $Fld[0] eq 'numSeq_src' ) {$numSeq_src = $Fld[2];}
      if ( $Fld[0] eq 'Seq_src(' ) {$Seq_src[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq 'numGamma' ) {$numGamma = $Fld[2];}
      if ( $Fld[0] eq 'Gamma_list(' ) {$Gamma_list[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '&lattis' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'nrow(' ) {$nrow[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq 't_srce(' ) {$t_srce[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
    $found = 1;
  }
}

printf "nrow = [%d,%d,%d,%d]\n", $nrow[0], $nrow[1], $nrow[2], $nrow[3];

$L_t = $nrow[$j_decay];
$L_s = $nrow[0];
printf "L=%d  numKappa=%d  dir=%s\n", $L_t, $numKappa, $dir;
printf "numSeq_src=%d  numGamma=%d\n", $numSeq_src, $numGamma;
foreach $x (0 .. $numSeq_src-1)
{ 
  $i = $Seq_src[$x];
  print "Seq_src[$x] = $i";
}

#
# Setup file names
#
$hand_inc = 0;

foreach $x (0 .. $numSeq_src-1)
{ 
  $i = $Seq_src[$x];

  # Get a unique handle
  $s = 'seq_hadron';
  $file = $s . "_s" . $i;
	    
  $bar_tag{$s, $i} = 0;
  $bar_file{$s, $i} = $file;

  if ( -f $file ) {unlink($file);}
}


foreach $x (0 .. $numSeq_src-1)
{ 
  $i = $Seq_src[$x];

  foreach $gg (0 .. $numGamma-1)
  { 
    $g = $Gamma_list[$gg];

    foreach $x (-2 .. 2)
    {
      foreach $y (-2 .. 2)
      {
	foreach $z (-2 .. 2)
	{
	  $p2 = $x*$x + $y*$y + $z*$z;
	  next if ($p2 > 5);

	  foreach $s ('local_cur3ptfn', 'nonlocal_cur3ptfn')
	  {
	    # Get a unique handle
	    $file = $s . "_s" . $i . "_g" . $g . "_qx" . $x . "_qy" . $y . "_qz" . $z;
	    
	    $bar_tag{$s, $i, $g, $x, $y, $z} = 0;
	    $bar_file{$s, $i, $g, $x, $y, $z} = $file;

	    if ( -f $file ) {unlink($file);}
	  }
	}
      }
    }
  }
}


#
# Stringize smearing names
#
for($i=0; $i < $numKappa; ++$i)
{
  $Kappa_orig[$i] = $Kappa[$i];

  # Kappa values as an integer
  $Kappa_i = int(10000*$Kappa[$i] + 0.5);

  # Width names. Remove trailing zeros
  $wvf_param_i = $wvf_param[$i];
  $wvf_param_i =~ s/\./p/;
  $wvf_param_i =~ s/[0]*$//;
  $wvf_param_i =~ s/p$//;

  if ( $Wvf_kind == 1 ) 
  {
    $wvf_param_i = 'g' . $wvf_param_i;
  }
  elsif ( $Wvf_kind == 2 )
  {
    $wvf_param_i = 'e' . $wvf_param_i;
  }
  elsif ( $Wvf_kind == 3 )
  {
    $wvf_param_i = 'G' . $wvf_param_i;
  }
  elsif ( $Wvf_kind == 4 )
  {
    $wvf_param_i = 'w' . $wvf_param_i;
  }
    
  $Kappa_orig[$i] = $Kappa[$i];
  $Kappa[$i] = $Kappa_i;
  $wvf_param[$i] = $wvf_param_i;

  print $i,$Kappa[$i],$wvf_param[$i];
}


#
# Write a info (configuration) file  to use with  formfact.pl
#
open(CONF, "> config.pl") || die "Error opening configuration file for writing\n";

die "Currently only support 1 Kappa in config.pl\n" unless $numKappa == 1;

$t_src = $t_srce[$j_decay];

print CONF "\$numKappa = 1;";
$Kappa_12 = sprintf("D%d", $Kappa[0]);
$ms = "_1";

if ( $Pt_src == 1 )
{
  $spext = $Kappa_12 . '.P' . $ms . '.P' . $ms . '.PP';
  $ssext = $spext;
}
if ( $Sl_src == 1 )
{
  $spext = $Kappa_12 . '.D' . $wvf_param[0] . $ms . '.P' . $ms . '.SP';
  $ssext = $Kappa_12 . '.D' . $wvf_param[0] . $ms . '.D' . $wvf_param[0] . $ms . '.SS';
}

if ($FermAct == 0) 
{
  $norm = "2*$Kappa_orig[0]";
  $nam = 'nonlocal';
  $cur = 'n';
} 
elsif ($FermAct == 23) 
{
  $norm = "1";
  $nam = 'local';
  $cur = 'l';
}
else
{
  die "Unknown fermion action FermAct = $FermAct\n";
}


printf CONF "\$nam = \'$nam\';\n";
printf CONF "\$cur = \'$cur\';\n";

printf CONF "\$L_s  = $L_s;\n";
printf CONF "\$L_t  = $L_t;\n";
print CONF "\$t_src = $t_src;";
print CONF "\$t_snk = $t_snk;";
printf CONF "\$spext = \'$spext\';\n";
printf CONF "\$ssext = \'$ssext\';\n";
printf CONF "\$norm  = \'$norm\';\n";

printf CONF "\@sink_mom  = ($sink_mom[0]";
for($i=1; $i < $Nd-1; ++$i)
{
  printf CONF ",$sink_mom[$i]";
}
printf CONF ")\n";

close(CONF);


#
# Now, scan for correlation functions
#
printf "Searching for correlation functions....\n";

while (<STDIN>)
{
  chop;	# strip record separator
  @Fld = split(' ', $_);

  if ( $Fld[0] eq '&End_Wilson_3Pt_fn_measurements' ||  
       $Fld[0] eq '&End_Wilson_Baryon_3Pt_fn_measurements' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }

    # Keep track of the number of propagators
    $number++;
    printf "printed correlators for configuration  $number\n";
  }


  if ( $Fld[0] eq '&Wilson_Baryon_3Pt_fn_measurements' || $Fld[0] eq '&Wilson_Meson_3Pt_fn_measurements' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'seq_src' ) {$seq_src = $Fld[2];$seq_src_value = $seq_src;}
      if ( $Fld[0] eq 'seq_src_value' ) {$seq_src_value = $Fld[2];$seq_src = $seq_src_value;}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }

  if ( $Fld[0] eq '' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'seq_src' ) {$seq_src = $Fld[2];$seq_src_value = $seq_src;}
      if ( $Fld[0] eq 'seq_src_value' ) {$seq_src_value = $Fld[2];$seq_src = $seq_src_value;}
      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }


  if ( $Fld[0] eq '&Wilson_hadron_2Pt_fn' )
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq 'seq_hadron_0' ) 
      {
	$s = 'seq_hadron';
	$file = $bar_file{$s, $seq_src};
	$tag = ++$bar_tag{$s, $seq_src};

	if ($tag == 1)
	{
#	  print "Open seq_hadron $file";

	  $bar_hand{$s, $seq_src} = $hand = ("foo" . ++$hand_inc);
	  open($hand, ">" . $file) || die "Unable to open seq_hadron file $file\n";
	      
	  # Write header
	  printf $hand "XXXXXX 1 1 %d 1\n", $L_t, $L_s;
	}
	else
	{
#	  print "Reopen seq_hadron $file";

	  $hand = $bar_hand{$s, $seq_src};
	  open($hand, ">>" . $file) || die "Unable to reopen seq_hadron file $file\n";
	}

	printf $hand "0 %s %s\n", $Fld[3], $Fld[5];
	close($hand);
      }

      if ( $Fld[0] eq '&END' ) {$end = 1;}
    }
  }


  if ( $Fld[0] eq '&Wilson_Current_3Pt_fn')   # old format
  {
    $_ = <STDIN>; chop; @Fld = split(' ', $_);
    die "no mu\n" if ( $Fld[0] ne 'mu' );
    $mu = $Fld[2];
    $_ = <STDIN>; chop; @Fld = split(' ', $_);
    die "no j_decay\n" if ( $Fld[0] ne 'j_decay' );
    foreach $i (1 .. $Nd-1)
    {
      $_ = <STDIN>; chop; @Fld = split(' ', $_);
      die "no inser_mom\n" if ( $Fld[0] ne 'inser_mom(' );
      $inser_mom[$Fld[1]] = $Fld[4];
    }

    $gamma_value = 1 << $mu;

    $s = 'local_cur3ptfn';
    $file = $bar_file{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
    $tag = ++$bar_tag{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

    if ($tag == 1)
    {
#      print "Open local_cur3ptfn $file";

      $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]} = 
	$hand = ("foo" . ++$hand_inc);
      open($hand, ">" . $file) || die "Unable to open local_cur3ptfn file $file\n";
	      
      # Write header
      printf $hand "XXXXXX %d 1 %d 1\n", $L_t, $L_s;
    }
    else
    {
#      print "Reopen local_cur3ptfn $file";

      $hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
      open($hand, ">>" . $file) || die "Unable to reopen local_cur3ptfn file $file\n";
    }

    foreach $t (0 .. $L_t-1)
    {
      $_ = <STDIN>; chop; @Fld = split(' ', $_);
      $tt = $Fld[1];
      printf $hand "%d %s %s\n", $tt, $Fld[5], $Fld[7];
    }

    close($hand);

    $s = 'nonlocal_cur3ptfn';
    $file = $bar_file{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
    $tag = ++$bar_tag{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

    if ($tag == 1)
    {
#      print "Open nonlocal_cur3ptfn $file";

      $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]} = 
	$hand = ("foo" . ++$hand_inc);
      open($hand, ">" . $file) || die "Unable to open nonlocal_cur3ptfn file $file\n";
	      
      # Write header
      printf $hand "XXXXXX %d 1 %d 1\n", $L_t, $L_s;
    }
    else
    {
#      print "Reopen nonlocal_cur3ptfn $file";

      $hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
      open($hand, ">>" . $file) || die "Unable to reopen nonlocal_cur3ptfn file $file\n";
    }

    $hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

    foreach $t (0 .. $L_t-1)
    {
      $_ = <STDIN>; chop; @Fld = split(' ', $_);
      $tt = $Fld[1];
      printf $hand "%d %s %s\n", $tt, $Fld[5], $Fld[7];
    }

    close($hand);

    $_ = <STDIN>; chop; @Fld = split(' ', $_);
    die "no END\n" if ( $Fld[0] ne '&END' );
  }


  if ( $Fld[0] eq '&Wilson_Local_Current_3Pt_fn')
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq '&END' ) {$end = 1;}
      if ( $Fld[0] eq 'gamma_value' ) {$gamma_value = $Fld[2];}
      if ( $Fld[0] eq 'j_decay' ) {$j_decay = $Fld[2];}
      if ( $Fld[0] eq 'inser_mom(' ) {$inser_mom[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq 'local_cur3ptfn(' )
      {
	$s = 'local_cur3ptfn';
	$file = $bar_file{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
	$tag = ++$bar_tag{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

	if ($tag == 1)
	{
#	  print "Open local_cur3ptfn $file";

	  $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]} = 
	    $hand = ("foo" . ++$hand_inc);
	  open($hand, ">" . $file) || die "Unable to open local_cur3ptfn file $file\n";
	      
	  # Write header
	  printf $hand "XXXXXX %d 1 %d 1\n", $L_t, $L_s;
	}
	else
	{
#         print "Reopen local_cur3ptfn $file";

	  $hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
	  open($hand, ">>" . $file) || die "Unable to reopen local_cur3ptfn file $file\n";
	}

	printf $hand "%d %s %s\n", $Fld[1], $Fld[5], $Fld[7];
	foreach $t (1 .. $L_t-1)
	{
	  $_ = <STDIN>; chop; @Fld = split(' ', $_);
	  printf $hand "%d %s %s\n", $Fld[1], $Fld[5], $Fld[7];
	}

	close($hand);
      }
    }
  }

  if ( $Fld[0] eq '&Wilson_NonLocal_Current_3Pt_fn')
  {
    $end = 0;
    while ($end == 0)
    {
      $_ = <STDIN>;
      chop;
      @Fld = split(' ', $_);

      if ( $Fld[0] eq '&END' ) {$end = 1;}
      if ( $Fld[0] eq 'gamma_value' ) {$gamma_value = $Fld[2];}
      if ( $Fld[0] eq 'j_decay' ) {$j_decay = $Fld[2];}
      if ( $Fld[0] eq 'inser_mom(' ) {$inser_mom[$Fld[1]] = $Fld[4];}
      if ( $Fld[0] eq 'nonlocal_cur3ptfn(' )
      {
	$s = 'nonlocal_cur3ptfn';
	$file = $bar_file{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
	$tag = ++$bar_tag{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

	if ($tag == 1)
	{
#         print "Open nonlocal_cur3ptfn $file";

	  $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]} = 
	    $hand = ("foo" . ++$hand_inc);
	  open($hand, ">" . $file) || die "Unable to open nonlocal_cur3ptfn file $file\n";
	      
	  # Write header
	  printf $hand "XXXXXX %d 1 %d 1\n", $L_t, $L_s;
	}
	else
	{
#         print "Reopen nonlocal_cur3ptfn $file";

	  $hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};
	  open($hand, ">>" . $file) || die "Unable to reopen nonlocal_cur3ptfn file $file\n";
	}
	
	$hand = $bar_hand{$s, $seq_src, $gamma_value, $inser_mom[0], $inser_mom[1], $inser_mom[2]};

	printf $hand "%d %s %s\n", $Fld[1], $Fld[5], $Fld[7];
	foreach $t (1 .. $L_t-1)
	{
	  $_ = <STDIN>; chop; @Fld = split(' ', $_);
	  printf $hand "%d %s %s\n", $Fld[1], $Fld[5], $Fld[7];
	}

	close($hand);
      }
    }
  }
}


printf "Found  %d  propagators\n", $number;

# Prepend a line containing the number of propagators found
printf "prepending number of measurements found\n";
foreach $X (keys %bar_tag)
{
  $tag  = $bar_tag{$X};

  if ( $tag > 0 )
  {
    $name = $bar_file{$X};

#    printf "Header  %s\n", $name;
    
    open(FILE, "+<" . $name) || die "could not reopen $name for rewrite\n";
    printf FILE "%6d", $number;
    close(FILE);
  }
#
#  else
#  {
#   @xs = split($;, $X);
#   printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n\n",
#      $xs[1], $xs[2], $xs[2], $xs[4], 
#      $xs[4], $xs[6], $xs[7], $xs[8], $xs[9], $xs[10];
#  }  

}

