#
#  $Id: regres.pl,v 3.1 2007-03-14 18:21:47 bjoo Exp $
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
	 input       => "$test_dir/chroma/hadron_s/spectrum_s/spectrum_s.ini.xml" , 
	 output      => "inline_spectrum_s.candidate.xml",
	 metric      => "$test_dir/chroma/hadron_s/spectrum_s/spectrum_s.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron_s/spectrum_s/spectrum_s.out.xml" ,
     }
     );
