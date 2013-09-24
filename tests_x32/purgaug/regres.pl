#
#  $Id: regres.pl,v 3.3 2008-01-23 15:44:06 edwards Exp $
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
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "purgaug" , 
	 input       => "$test_dir/purgaug/purgaug.2+2.ini.xml" , 
	 output      => "purgaug.2+2.candidate.xml",
	 metric      => "$test_dir/purgaug/purgaug.2+2.metric.xml" ,
	 controlfile => "$test_dir/purgaug/purgaug.2+2.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "purgaug" , 
	 input       => "$test_dir/purgaug/purgaug.2+2.1loop.ini.xml" , 
	 output      => "purgaug.2+2.1loop.candidate.xml",
	 metric      => "$test_dir/purgaug/purgaug.2+2.1loop.metric.xml" ,
	 controlfile => "$test_dir/purgaug/purgaug.2+2.1loop.out.xml" ,
     }
     );
