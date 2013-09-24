#
#  $Id: regres.pl,v 1.1 2006-11-17 02:11:05 edwards Exp $
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
	 execute     => "t_meas_wilson_flow" , 
	 input       => "$test_dir/chroma/tests/t_meas_wilson_flow/t_meas_wilson_flow.in.xml" , 
	 output      => "t_meas_wilson_flow.candidate.xml",
	 metric      => "$test_dir/chroma/tests/t_meas_wilson_flow/sink_smearing_s.metric.xml" ,
	 controlfile => "$test_dir/chroma/tests/t_meas_wilson_flow/t_meas_wilson_flow.ref.xml" ,
     }
     );
