#
#  $Id: regres.pl,v 3.1 2006-04-19 02:36:15 edwards Exp $
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
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "purgaug" , 
	 input       => "$test_dir/purgaug/purgaug.sfnonpt.ini.xml" , 
	 output      => "purgaug.sfnonpt.candidate.xml",
	 metric      => "$test_dir/purgaug/purgaug.sfnonpt.metric.xml" ,
	 controlfile => "$test_dir/purgaug/purgaug.sfnonpt.out.xml" ,
     }
     );
