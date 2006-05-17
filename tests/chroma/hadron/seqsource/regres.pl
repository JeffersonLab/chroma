#
#  $Id: regres.pl,v 3.1 2006-05-17 20:18:19 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-meson.ini.xml" , 
	 output      => "seqsource-meson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-meson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-meson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-baryon.ini.xml" , 
	 output      => "seqsource-baryon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-baryon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-baryon.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/seqsource/seqsource-photon.ini.xml" , 
	 output      => "seqsource-photon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/seqsource/seqsource-photon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/seqsource/seqsource-photon.out.xml" ,
     },
     );
