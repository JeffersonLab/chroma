#
#  $Id: regres.pl,v 3.0 2006-04-03 04:59:33 edwards Exp $
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
	 execute     => "spectrum_s" , 
	 input       => "$test_dir/spectrum_s/spectrum_s_reg.xml" , 
	 output      => "spectrum_s.candidate.xml",
	 metric      => "$test_dir/spectrum_s/spectrum_s.metric.xml" ,
	 controlfile => "$test_dir/spectrum_s/spectrum_s_reg_test_output.xml" ,
     }
     );
