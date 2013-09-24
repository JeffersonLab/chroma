#
#  $Id: regres.pl,v 3.1 2006-06-11 06:30:36 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/noisy_building_blocks/nbb.ini.xml" , 
	 output      => "nbb.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/noisy_building_blocks/nbb.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/noisy_building_blocks/nbb.out.xml" ,
     }
     );
