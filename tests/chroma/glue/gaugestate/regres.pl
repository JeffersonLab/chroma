#
#  $Id: regres.pl,v 1.2 2006-09-26 13:47:27 edwards Exp $
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
	 input       => "$test_dir/chroma/glue/gaugestate/gaugestate.ini.xml" , 
	 output      => "gaugestate.candidate.xml",
	 metric      => "$test_dir/chroma/glue/gaugestate/gaugestate.metric.xml" ,
	 controlfile => "$test_dir/chroma/glue/gaugestate/gaugestate.out.xml" ,
     }
     );
