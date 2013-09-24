#
#  $Id: regres.pl,v 1.1 2008-04-24 14:04:43 edwards Exp $
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
	 input       => "$test_dir/chroma/glue/qtop_naive/qtop_naive.ini.xml" , 
	 output      => "qtop_naive.candidate.xml",
	 metric      => "$test_dir/chroma/glue/qtop_naive/qtop_naive.metric.xml" ,
	 controlfile => "$test_dir/chroma/glue/qtop_naive/qtop_naive.out.xml" ,
     }
     );
