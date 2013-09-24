#
#  $Id: regres.pl,v 3.1 2006-04-13 19:25:28 edwards Exp $
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
	 input       => "$test_dir/chroma/io/szin_read_obj/szin_read_obj.ini.xml" , 
	 output      => "szin_read_obj.candidate.xml",
	 metric      => "$test_dir/chroma/io/szin_read_obj/szin_read_obj.metric.xml" ,
	 controlfile => "$test_dir/chroma/io/szin_read_obj/szin_read_obj.out.xml" ,
     }
     );
