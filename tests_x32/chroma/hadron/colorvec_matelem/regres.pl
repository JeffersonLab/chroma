#
#  $Id: regres.pl,v 1.2 2008-08-06 15:27:21 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/colorvec_matelem/colorvec_matelem.ini.xml" , 
	 output      => "colorvec_matelem.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/colorvec_matelem/colorvec_matelem.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/colorvec_matelem/colorvec_matelem.out.xml" ,
     }

     );
