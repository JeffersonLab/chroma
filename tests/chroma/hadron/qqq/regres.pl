#
#  $Id: regres.pl,v 3.2 2007-02-25 22:49:21 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/qqq/qqq.ini.xml" , 
	 output      => "qqq.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/qqq/qqq.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/qqq/qqq.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/qqq/qqdiquark.ini.xml" , 
	 output      => "qqdiquark.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/qqq/qqdiquark.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/qqq/qqdiquark.out.xml" ,
     }
     );
