#
#  $Id: regres.pl,v 1.2 2005-12-25 06:51:32 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/spectrumOct/spectrumOct-v11.ini.xml" , 
	 output      => "spectrumOct-v11.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/spectrumOct/spectrumOct-v11.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/spectrumOct/spectrumOct-v11.out.xml" ,
     }
     );
