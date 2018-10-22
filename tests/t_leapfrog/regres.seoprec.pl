#
#  $Id: regres.tprec.pl,v 1.1 2008-01-09 19:03:06 bjoo Exp $
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
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.seoprec_logdet_diag.ini.xml" , 
	 log         => "t_leapfrog.seoprec_logdet_diag.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.seoprec_logdet_diag.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.seoprec_logdet_diag.log.xml" ,
     },
  


     );
