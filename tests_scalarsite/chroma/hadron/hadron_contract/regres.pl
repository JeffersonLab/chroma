#
#  $Id: regres.pl,v 1.2 2007-08-11 22:48:12 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/hadron_contract/meson_2pt/meson_2pt.ini.xml" , 
	 output      => "meson_2pt.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/hadron_contract/meson_2pt/meson_2pt.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/hadron_contract/meson_2pt/meson_2pt.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/hadron_contract/stoch_condensates/stoch_condensates.ini.xml" , 
	 output      => "stoch_condensates.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/hadron_contract/stoch_condensates/stoch_condensates.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/hadron_contract/stoch_condensates/stoch_condensates.out.xml" ,
     }
     );
