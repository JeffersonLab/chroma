#
#  $Id: regres.pl,v 1.1 2005-12-25 18:41:54 edwards Exp $
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
	 input       => "$test_dir/chroma/smear/link_smear/link_smear.ini.xml" , 
	 output      => "link_smear.candidate.xml",
	 metric      => "$test_dir/chroma/smear/link_smear/link_smear.metric.xml" ,
	 controlfile => "$test_dir/chroma/smear/link_smear/link_smear.out.xml" ,
     }
     );
