#
#  $Id: regres.pl,v 3.2 2006-12-02 18:08:26 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/building_blocks/bb-baryon.ini.xml" , 
	 output      => "bb-baryon.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/building_blocks/bb-baryon.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/building_blocks/bb-baryon.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/building_blocks/bb-meson.ini.xml" , 
	 output      => "bb-meson.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/building_blocks/bb-meson.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/building_blocks/bb-meson.out.xml" ,
     },
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/hadron/building_blocks/bb-deriv.ini.xml" , 
	 output      => "bb-deriv.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/building_blocks/bb-deriv.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/building_blocks/bb-deriv.out.xml" ,
     }
     );
