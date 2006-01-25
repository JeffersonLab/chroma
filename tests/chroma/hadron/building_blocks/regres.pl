#
#  $Id: regres.pl,v 2.1 2006-01-25 05:52:34 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/building_blocks/bb-v2.ini.xml" , 
	 output      => "bb-v2.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/building_blocks/bb-v2.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/building_blocks/bb-v2.out.xml" ,
     }
     );
