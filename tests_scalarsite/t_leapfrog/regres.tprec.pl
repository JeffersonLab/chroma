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
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_aniso_apbc.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_aniso_apbc.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_aniso_apbc.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_schr1link.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_schr1link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT_schr1link.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson.log.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_schr1link.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson_schr1link.log.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_schr1link.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_schr2link.ini.xml" , 
	 log         => "t_leapfrog.unprec_s_cprec_t_wilson_schr2link.log.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_schr2link.log.xml" ,
     },  
    {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_aniso_apbc.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_aniso_apbc.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_aniso_apbc.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_schr1link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_schr1link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_logdetTT_schr1link.log.xml" ,
     },  
    {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_schr1link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson_schr1link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_schr1link.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_schr2link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_wilson_schr2link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_wilson_schr2link.log.xml" ,
     },  

    {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover_logdetTT.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_aniso_apbc.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_aniso_apbc.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_aniso_apbc.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_schr1link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_schr1link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson_logdetTT.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_logdetTT_schr1link.log.xml" ,
     },  
    {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover.log.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_schr1link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover_schr1link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_schr1link.log.xml" ,
     },  
     {
	 exec_path   => "$top_builddir/mainprogs/tests" , 
	 execute     => "t_leapfrog" , 
	 input       => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_schr2link.ini.xml" , 
	 log         => "t_leapfrog.iluprec_s_cprec_t_clover_schr2link.candidate.xml",
	 metric      => "$test_dir/t_leapfrog/t_leapfrog.unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/t_leapfrog/t_leapfrog.iluprec_s_cprec_t_clover_schr2link.log.xml" ,
     },  



     );
