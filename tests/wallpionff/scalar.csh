#!/bin/tcsh

set dir = $HOME/qcd/chroma/scalar/mainprogs/main

set cfgtransf_EXEC = $dir/cfgtransf
set make_source_EXEC = $dir/make_source
set prop_EXEC = $dir/propagator
set qproptransf_EXEC = $dir/qproptransf
set spectro_EXEC = $dir/spectrum_w
set seqprop_EXEC = $dir/seqprop
set wallformfac_EXEC = $dir/wallformfac

set run = ""

/bin/rm -f scalar.output
touch scalar.output

/bin/cp make_source_sh.ini DATA
$run $make_source_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT make_source_sh.xml

/bin/cp make_source_wl.ini DATA
$run $make_source_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT make_source_wl.xml

/bin/cp prop_sh.ini DATA
$run $prop_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT prop_sh.xml

/bin/cp prop_wl.ini DATA
$run $prop_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT prop_wl.xml

/bin/cp spectrum_sh.ini DATA
$run $spectro_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT spectro_sh.xml

/bin/cp spectrum_wl.ini DATA
$run $spectro_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT spectro_wl.xml

/bin/cp wallformfac_v1.xml.ini DATA
$run $wallformfac_EXEC < /dev/null >>& play.output
/bin/mv XMLDAT wallformfac.xml

