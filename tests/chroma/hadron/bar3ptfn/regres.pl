#
#  $Id: regres.pl,v 3.1 2006-04-14 18:43:14 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/bar3ptfn/bar3ptfn.ini.xml" , 
	 output      => "bar3ptfn.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/bar3ptfn/bar3ptfn.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/bar3ptfn/bar3ptfn.out.xml" ,
     }
     );
