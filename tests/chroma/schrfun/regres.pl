#
#  $Id: regres.pl,v 1.2 2007-08-27 22:47:00 edwards Exp $
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
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/schrfun/schrfun-hadspec.ini.xml" , 
	 output      => "schrfun-hadspec.candidate.xml",
	 metric      => "$test_dir/chroma/schrfun/schrfun-hadspec.metric.xml" ,
	 controlfile => "$test_dir/chroma/schrfun/schrfun-hadspec.out.xml" ,
     }
     );
