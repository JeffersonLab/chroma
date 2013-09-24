#
#  $Id: regres.pl,v 1.1 2006-11-07 23:11:13 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/static_light/static_light.ini.xml" , 
	 output      => "static_light.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/static_light/static_light.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/static_light/static_light.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/static_light/static_light_pot.ini.xml" , 
	 output      => "static_light_pot.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/static_light/static_light_pot.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/static_light/static_light_pot.out.xml" ,
     }
     );
