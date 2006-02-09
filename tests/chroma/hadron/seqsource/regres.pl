#
#  $Id: regres.pl,v 1.3 2006-02-09 02:19:57 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-meson-v1.ini.xml" , 
	 output      => "seqsource-meson-v1.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-meson-v1.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-meson-v1.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-baryon-v1.ini.xml" , 
	 output      => "seqsource-baryon-v1.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-baryon-v1.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-baryon-v1.out.xml" ,
     }
     );
