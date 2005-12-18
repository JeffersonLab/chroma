#
#  $Id: regres.pl,v 1.1 2005-12-18 23:51:30 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-v2.ini.xml" , 
	 output      => "seqsource-v2.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-v2.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-v2.out.xml" ,
     }
     );
