#
#  $Id: regres.pl,v 3.1 2006-12-10 01:53:24 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/qqqNucNuc/qqqNucNuc.ini.xml" , 
	 output      => "qqqNucNuc.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/qqqNucNuc/qqqNucNuc.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/qqqNucNuc/qqqNucNuc.out.xml" ,
     }
     );
