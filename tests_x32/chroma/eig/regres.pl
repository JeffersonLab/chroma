#
#  $Id: regres.pl,v 3.1 2006-09-21 19:45:55 edwards Exp $
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
	 input       => "$test_dir/chroma/eig/eigbnds.ini.xml" , 
	 output      => "eigbnds.candidate.xml",
	 metric      => "$test_dir/chroma/eig/eigbnds.metric.xml" ,
	 controlfile => "$test_dir/chroma/eig/eigbnds.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/eig/ritz_Hw-v1.ini.xml" , 
	 output      => "ritz_Hw-v1.candidate.xml",
	 metric      => "$test_dir/chroma/eig/ritz_Hw-v1.metric.xml" ,
	 controlfile => "$test_dir/chroma/eig/ritz_Hw-v1.out.xml" ,
     }
     );
