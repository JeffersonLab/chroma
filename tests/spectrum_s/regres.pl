#
#  $Id: regres.pl,v 3.1 2006-08-24 21:07:11 edwards Exp $
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
	 input       => "$test_dir/spectrum_s/spectrum_s.ini.xml" , 
	 output      => "spectrum_s.candidate.xml",
	 metric      => "$test_dir/spectrum_s/spectrum_s.metric.xml" ,
	 controlfile => "$test_dir/spectrum_s/spectrum_s.out.xml" ,
     }
     );
