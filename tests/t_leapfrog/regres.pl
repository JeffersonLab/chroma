#
#  $Id: regres.pl,v 1.3 2006-01-15 21:09:46 bjoo Exp $
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
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_wilson.ini.xml" , 
	 output      => "t_leapfrog.unprec_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_wilson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_hasenbusch_wilson.ini.xml" , 
	 output      => "t_leapfrog.unprec_hasenbusch_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_hasenbusch_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_hasenbusch_wilson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.prec_hasenbusch_wilson.ini.xml" , 
	 output      => "t_leapfrog.prec_hasenbusch_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.prec_hasenbusch_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.prec_hasenbusch_wilson.out.xml" ,
     }
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.prec_1flav_wilson.ini.xml" , 
	 output      => "t_leapfrog.prec_1flav_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.prec_1flav_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.prec_1flav_wilson.out.xml" ,
     }
     );
