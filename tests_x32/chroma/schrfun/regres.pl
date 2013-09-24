#
#  $Id: regres.pl,v 1.4 2007-08-28 17:24:59 uid3790 Exp $
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
	 input       => "$test_dir/chroma/schrfun/schrfun-pcac.ini.xml" , 
	 output      => "schrfun-pcac.candidate.xml",
	 metric      => "$test_dir/chroma/schrfun/schrfun-pcac.metric.xml" ,
	 controlfile => "$test_dir/chroma/schrfun/schrfun-pcac.out.xml" ,
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
