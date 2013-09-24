#
#  $Id: regres.pl,v 1.1 2006-09-21 20:18:15 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/fermstate/fermstate.ini.xml" , 
	 output      => "fermstate.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/fermstate/fermstate.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/fermstate/fermstate.out.xml" ,
     }
     );
