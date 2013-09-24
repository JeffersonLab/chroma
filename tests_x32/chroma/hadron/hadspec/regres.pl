#
#  $Id: regres.pl,v 1.1 2006-07-10 19:03:30 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/hadspec/hadspec.ini.xml" , 
	 output      => "hadspec.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/hadspec/hadspec.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/hadspec/hadspec.out.xml" ,
     }
     );
