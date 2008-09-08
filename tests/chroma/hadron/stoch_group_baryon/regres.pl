#
#  $Id: regres.pl,v 1.2 2008-09-08 21:38:02 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/stoch_group_baryon/stoch_group_baryon.ini.xml" , 
	 aux_files   => ["$test_dir/chroma/hadron/stoch_group_baryon/Nucleon_elem"] , 
	 output      => "stoch_group_baryon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_group_baryon/stoch_group_baryon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_group_baryon/stoch_group_baryon.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "make_baryon_ops" , 
	 input       => "$test_dir/chroma/hadron/stoch_group_baryon/make_baryon_ops.ini.xml" , 
	 aux_files   => ["$test_dir/chroma/hadron/stoch_group_baryon/Nucleon_G1g_1"] , 
	 output      => "make_baryon_ops.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_group_baryon/make_baryon_ops.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_group_baryon/make_baryon_ops.out.xml" ,
     }
     );
