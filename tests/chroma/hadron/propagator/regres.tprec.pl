#
#  $Id: regres.tprec.pl,v 3.2 2008-01-08 17:52:08 bjoo Exp $
#
#  This is the portion of a script this is included recursively
#

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#
@regres_list = 
    (

# Commented these out so they don't break nightlies. They work but need
# time to be local...  Will fix later by mucking with the run script

    {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson.ini.xml" , 
	 output      => "unprec_s_cprec_t_wilson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson.out.xml" ,
     },   
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr1link.ini.xml" , 
	 output      => "unprec_s_cprec_t_wilson_schr1link.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr1link.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr1link.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link.ini.xml" , 
	 output      => "unprec_s_cprec_t_wilson_schr2link.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link_space.ini.xml" , 
	 output      => "unprec_s_cprec_t_wilson_schr2link_space.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link_space.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_s_cprec_t_wilson_schr2link_space.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson.ini.xml" ,
        output      => "iluprec_s_cprec_t_wilson.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr1link.ini.xml" ,
        output      => "iluprec_s_cprec_t_wilson_schr1link.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr1link.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr1link.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link.ini.xml" ,
        output      => "iluprec_s_cprec_t_wilson_schr2link.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link_space.ini.xml" ,
        output      => "iluprec_s_cprec_t_wilson_schr2link_space.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link_space.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_wilson_schr2link_space.out.xml" ,
     },
    {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover.ini.xml" , 
	 output      => "iluprec_s_cprec_t_clover.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover.out.xml" ,
     },
    {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr1link.ini.xml" , 
	 output      => "iluprec_s_cprec_t_clover_schr1link.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr1link.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr1link.out.xml" ,
     },
    {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link.ini.xml" , 
	 output      => "iluprec_s_cprec_t_clover_schr2link.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link.out.xml" ,
     },
    {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link_space.ini.xml" , 
	 output      => "iluprec_s_cprec_t_clover_schr2link_space.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link_space.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/iluprec_s_cprec_t_clover_schr2link_space.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_wilson.ini.xml" ,
        output      => "eo3dprec_s_cprec_t_wilson.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_wilson.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_wilson.out.xml" ,
     },
     {
        exec_path   => "$top_builddir/mainprogs/main" ,
        execute     => "chroma" ,
        input       => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_clover.ini.xml" ,
        output      => "eo3dprec_s_cprec_t_clover.candidate.xml",
        metric      => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_clover.metric.xml" ,
        controlfile => "$test_dir/chroma/hadron/propagator/eo3dprec_s_cprec_t_clover.out.xml" ,
     }



     );
