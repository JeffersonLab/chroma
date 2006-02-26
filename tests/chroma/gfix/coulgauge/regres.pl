#
#  $Id: regres.pl,v 1.1 2006-02-26 04:32:39 edwards Exp $
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
	 input       => "$test_dir/chroma/gfix/coulgauge/coulgauge.ini.xml" , 
	 output      => "coulgauge.candidate.xml",
	 metric      => "$test_dir/chroma/gfix/coulgauge/coulgauge.metric.xml" ,
	 controlfile => "$test_dir/chroma/gfix/coulgauge/coulgauge.out.xml" ,
     }
     );
