# $Id: config.pl,v 1.3 2004-09-11 21:58:46 edwards Exp $
#
# This is a sample config.pl to use in the various formfac
# analysis scripts.
#
$mom2_max = 1;
@sink_mom = (0,0,0);
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
$X = 'W';  $Xext = $Aext  # used for energy only



$spext = "DG3_1.P_1.SP';
$ssext = 'D-7104.DG3_1.DG3_1.SS';
$swext = 'D-7104.DG3_1.W_1.SW';
$wpext = 'D-7104.W_1.P_1.WP';
$wsext = 'D-7104.W_1.DG3_1.WS';
$wwext = 'D-7104.W_1.W_1.WW';
$norm = '1';

