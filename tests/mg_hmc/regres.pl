#
#  $Id: regres.pl,v 3.4 2007-04-16 15:49:49 bjoo Exp $
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
	 execute     => "mg_hmc" , 
	 input       => "$test_dir/hmc/hmc-hfa.ini.xml" , 
	 log         => "hmc-hfa.candidate.xml",
	 metric      => "$test_dir/mg_hmc/hmc-hfa.metric.xml" ,
	 controlfile => "$test_dir/mg_hmc/hmc-hfa.log.xml" ,
     }
     );
