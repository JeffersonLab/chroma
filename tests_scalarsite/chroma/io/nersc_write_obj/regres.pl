#
#  $Id: regres.pl,v 1.1 2006-04-27 01:53:41 edwards Exp $
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
	 input       => "$test_dir/chroma/io/nersc_write_obj/nersc_write_obj.ini.xml" , 
	 output      => "nersc_write_obj.candidate.xml",
	 metric      => "$test_dir/chroma/io/nersc_write_obj/nersc_write_obj.metric.xml" ,
	 controlfile => "$test_dir/chroma/io/nersc_write_obj/nersc_write_obj.out.xml" ,
     }
     );
