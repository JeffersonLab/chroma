#
#  $Id: regres.pl,v 1.1 2005-12-18 23:51:30 edwards Exp $
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
	 input       => "$test_dir/t_leapfrog/t_leapfrog.prec_wilson.ini.xml" , 
	 output      => "t_leapfrog.prec_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.prec_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.prec_wilson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.prec_parwilson.ini.xml" , 
	 output      => "t_leapfrog.prec_parwilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.prec_parwilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.prec_parwilson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.prec_dwf.ini.xml" , 
	 output      => "t_leapfrog.prec_dwf.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.prec_dwf.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.prec_dwf.out.xml" ,
     }
     );
