#
#  $Id: regres.pl,v 1.1 2007-06-13 14:14:08 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/hadron_contract/hadron-meson.ini.xml" , 
	 output      => "hadron-meson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/hadron_contract/hadron-meson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/hadron_contract/hadron-meson.out.xml" ,
     }
     );
