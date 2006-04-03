# $Id: config.pl,v 3.0 2006-04-03 04:59:22 edwards Exp $
#
# This is a sample config.pl to use in the various formfac
# analysis scripts.
#
$mom2_max = 1;
#@sink_mom = (0,0,0);
@source_mom = (0,0,0);
$a = 0.1;
$L_s = 4;

$nam = 'nonlocal';
$cur = 'n';

$t_src = 1;
$t_snk = 6;

$Mass = "D-7104";

# wave channel
$L = 1;
$A = 'W';  $Aext = 'W';
$B = 'S';  $Bext = "DG1p2";
$C = 'S';  $Cext = "DG1p2";
$X = 'W';  $Xext = $Aext;    # used for energy only

$norm = '1';

