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
	 execute     => "hmc" , 
	 input       => "$test_dir/hmc/hmc.prec_wilson.ini.xml" , 
	 log         => "hmc.prec_wilson.candidate.xml",
	 metric      => "$test_dir/hmc/hmc.prec_wilson.metric.xml" ,
	 controlfile => "$test_dir/hmc/hmc.prec_wilson.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "hmc" , 
	 input       => "$test_dir/hmc/hmc.prec_wilson_repro.ini.xml" , 
	 log         => "hmc.prec_wilson_repro.candidate.xml",
	 metric      => "$test_dir/hmc/hmc.prec_wilson_repro.metric.xml" ,
	 controlfile => "$test_dir/hmc/hmc.prec_wilson_repro.log.xml" ,
     }
     );
