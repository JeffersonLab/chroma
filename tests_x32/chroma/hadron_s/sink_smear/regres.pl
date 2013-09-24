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
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron_s/sink_smear/sink_smearing_s.ini.xml" , 
	 output      => "sink_smearing_s.candidate.xml",
	 metric      => "$test_dir/chroma/hadron_s/sink_smear/sink_smearing_s.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron_s/sink_smear/sink_smearing_s.out.xml" ,
     }
     );
