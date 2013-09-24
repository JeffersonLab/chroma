#
#  $Id: regres.pl,v 1.2 2006-05-17 20:19:49 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/stoch_baryon/stoch_baryon.ini.xml" , 
	 output      => "stoch_baryon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_baryon/stoch_baryon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_baryon/stoch_baryon.out.xml" ,
     }
     );
