#
#  $Id: regres.pl,v 1.1 2008-05-02 21:01:22 bjoo Exp $
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
	 input       => "$test_dir/chroma/io/usqcd_ddpairs_prop/usqcd_prop_readwrite_internal.ini.xml" , 
	 output      => "usqcd_prop_readwrite_internal.candidate.xml",
	 metric      => "$test_dir/chroma/io/usqcd_ddpairs_prop/usqcd_prop_readwrite_internal.metric.xml" ,
	 controlfile => "$test_dir/chroma/io/usqcd_ddpairs_prop/usqcd_prop_readwrite_internal.out.xml" ,
     }
     );
