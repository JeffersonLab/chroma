#
#  $Id: regres.pl,v 1.1 2005-12-20 20:10:34 edwards Exp $
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
	 execute     => "purgaug" , 
	 input       => "$test_dir/purgaug/purgaug.ini.xml" , 
	 output      => "purgaug.candidate.xml",
	 metric      => "$test_dir/purgaug/purgaug.metric.xml" ,
	 controlfile => "$test_dir/purgaug/purgaug.out.xml" ,
     }
     );
