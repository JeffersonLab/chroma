#
#  $Id: regres.pl,v 1.1 2006-04-19 02:59:34 edwards Exp $
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
	 input       => "$test_dir/chroma/schrfun/schrfun.ini.xml" , 
	 output      => "schrfun.candidate.xml",
	 metric      => "$test_dir/chroma/schrfun/schrfun.metric.xml" ,
	 controlfile => "$test_dir/chroma/schrfun/schrfun.out.xml" ,
     }
     );
