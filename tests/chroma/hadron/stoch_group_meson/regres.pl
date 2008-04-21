#
#  $Id: regres.pl,v 1.1 2008-04-21 03:18:52 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/stoch_baryon/stoch_group_baryon.ini.xml" , 
	 output      => "stoch_group_baryon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/stoch_baryon/stoch_group_baryon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/stoch_baryon/stoch_group_baryon.out.xml" ,
     }
     );
