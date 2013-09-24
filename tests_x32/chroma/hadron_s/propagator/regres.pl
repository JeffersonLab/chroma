#
#  $Id: regres.pl,v 1.2 2006-12-07 18:36:53 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron_s/propagator/asqtad.ini.xml" , 
	 output      => "asqtad.candidate.xml",
	 metric      => "$test_dir/chroma/hadron_s/propagator/asqtad.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron_s/propagator/asqtad.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron_s/propagator/klein_gordon.ini.xml" , 
	 output      => "klein_gordon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron_s/propagator/klein_gordon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron_s/propagator/klein_gordon.out.xml" ,
     }
     );
