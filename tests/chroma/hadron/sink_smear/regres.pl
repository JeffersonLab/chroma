#
#  $Id: regres.pl,v 1.1 2006-02-22 03:23:48 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/sink_smear/sink_smearing-v5.ini.xml" , 
	 output      => "sink_smearing-v5.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/sink_smear/sink_smearing-v5.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/sink_smear/sink_smearing-v5.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/sink_smear/sink-displace.ini.xml" , 
	 output      => "sink-displace.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/sink_smear/sink-displace.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/sink_smear/sink-displace.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/sink_smear/sink-deriv.ini.xml" , 
	 output      => "sink-deriv.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/sink_smear/sink-deriv.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/sink_smear/sink-deriv.out.xml" ,
     }
     );
