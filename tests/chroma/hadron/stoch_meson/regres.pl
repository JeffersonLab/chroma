#
#  $Id: regres.pl,v 1.1 2006-05-09 17:47:18 edwards Exp $
#
#  This is the portion of a script this is included recursively
#

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#
@regres_list = 
    (
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/stoch_meson/stoch_meson.ini.xml" , 
	 output      => "stoch_meson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_meson/stoch_meson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_meson/stoch_meson.out.xml" ,
     }
     );
