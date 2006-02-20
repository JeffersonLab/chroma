#
#  $Id: regres.pl,v 1.1 2006-02-20 17:37:58 edwards Exp $
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
	 input       => "$test_dir/chroma/glue/wilslp/wilslp.ini.xml" , 
	 output      => "wilslp.candidate.xml",
	 metric      => "$test_dir/chroma/glue/wilslp/wilslp.metric.xml" ,
	 controlfile => "$test_dir/chroma/glue/wilslp/wilslp.out.xml" ,
     }
     );
