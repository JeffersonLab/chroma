#
#  $Id: regres.pl,v 1.3 2008-11-12 00:25:34 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/stoch_group_meson/stoch_group_meson.ini.xml" , 
	 aux_files   => ["$test_dir/chroma/hadron/stoch_group_meson/two_displace"] , 
	 output_dir  => "stoch_group_meson",
	 output      => "stoch_group_meson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_group_meson/stoch_group_meson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_group_meson/stoch_group_meson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "make_meson_ops" , 
	 input       => "$test_dir/chroma/hadron/stoch_group_meson/make_meson_ops.ini.xml" , 
	 aux_files   => ["$test_dir/chroma/hadron/stoch_group_meson/test_meson"] , 
	 output_dir  => "stoch_group_meson",
	 output      => "make_meson_ops.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_group_meson/make_meson_ops.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_group_meson/make_meson_ops.out.xml" ,
     }
     );
