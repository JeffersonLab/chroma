#
#  $Id: regres.pl,v 2.1 2005-12-25 18:29:39 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/mres/mres-v1.ini.xml" , 
	 output      => "mres-v1.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/mres/mres-v1.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/mres/mres-v1.out.xml" ,
     }
     );
