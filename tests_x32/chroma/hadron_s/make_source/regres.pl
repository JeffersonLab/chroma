#
#  $Id: regres.pl,v 1.1 2006-11-17 02:10:16 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron_s/make_source/make_source_s.ini.xml" , 
	 output      => "make_source_s.candidate.xml",
	 metric      => "$test_dir/chroma/hadron_s/make_source/make_source_s.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron_s/make_source/make_source_s.out.xml" ,
     }
     );
