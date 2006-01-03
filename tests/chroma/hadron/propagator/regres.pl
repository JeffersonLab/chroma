#
#  $Id: regres.pl,v 1.3 2006-01-03 05:03:46 bjoo Exp $
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
	 input       => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.ini.xml" , 
	 output      => "prec_wilson-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.out.xml" ,
     },
#     {
#	 exec_path   => "$top_builddir/mainprogs/main" , 
#	 execute     => "chroma" , 
#	 input       => "$test_dir/chroma/hadron/propagator/prec_clover-v9.ini.xml" , 
#	 output      => "prec_clover-v9.candidate.xml",
#	 metric      => "$test_dir/chroma/hadron/propagator/prec_clover-v9.metric.xml" ,
#	 controlfile => "$test_dir/chroma/hadron/propagator/prec_clover-v9.out.xml" ,
#     },
#     {
#	 exec_path   => "$top_builddir/mainprogs/main" , 
#	 execute     => "chroma" , 
#	 input       => "$test_dir/chroma/hadron/propagator/unprec_clover-v9.ini.xml" , 
#	 output      => "unprec_clover-v9.candidate.xml",
#	 metric      => "$test_dir/chroma/hadron/propagator/unprec_clover-v9.metric.xml" ,
#	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_clover-v9.out.xml" ,
#     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.ini.xml" , 
	 output      => "prec_wilson-twisted-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.ini.xml" , 
	 output      => "unprec-ovlap-zolo4d-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.ini.xml" , 
	 output      => "unprec_hamberwu-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.ini.xml" , 
	 output      => "prec_dwf-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.ini.xml" , 
	 output      => "prec_parwilson-v9.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.out.xml" ,
     }
     );
