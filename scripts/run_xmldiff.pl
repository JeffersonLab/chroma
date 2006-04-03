#
#  $Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/scripts/run_xmldiff.pl,v 3.0 2006-04-03 04:59:22 edwards Exp $
#
#  This is wrapper script to run the xmldiff application from
#  a makefile
#
#  Homepage for xmldiff application
#  https://forge.nesc.ac.uk/projects/xmldiff/
#
#
# More work:
#   At the moment this script assumes that the 
#   xmldiff application is in the users path.
#   Perhaps, this should be incorporated in the autoconf
#   tool chain.

# location of xmldiff
$xmldiff = "xmldiff" ;


# check whether xmldiff is available

$xml_check = `xmldiff` ; 
chop ($xml_check) ; 

if( $xml_check =~ /Usage/ )
{
#    print "xmldiff found\n" ; 
}
else
{
    print "Error:".$0." needs the xmldiff utility in your path\n" ; 
    print "Download it from http://forge.nesc.ac.uk/projects/xmldiff/ \n" ; 
    exit (0) ;
}



# at some stage this should be factored into another perl script

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#

 %HoH = (
	  t_propagator_fuzz_s => {
                       input       => "../../tests/t_asqtad_prop/INPUT_t_propagator_fuzz_s" , 
                       output      => "t_propagator_fuzz_s.xml",
                       metric       => "../../tests/t_asqtad_prop/t_propagator_fuzz_s_METRIC.xml" ,
                       controlfile  => "../../tests/t_asqtad_prop/t_propagator_fuzz_s.xml" ,
		   },

	  t_propagator_s => {
                       input  => "../../tests/t_asqtad_prop/DATA_v2" , 
                       output      => "t_propagator_s.xml",
                       metric       => "../../tests/t_asqtad_prop/t_propagator_s_METRIC.xml" ,
                       controlfile  => "../../tests/t_asqtad_prop/t_propagator_s.xml" ,
		   },


	  t_disc_loop_s => {
                       input       => "../../tests/t_asqtad_prop/DISC_DATA_v2" , 
                       output      =>  "t_disc_loop_s.xml",
                       metric       => "../../tests/t_asqtad_prop/t_disc_loop_s_METRIC.xml" ,
                       controlfile  => "../../tests/t_asqtad_prop/t_disc_loop_s.xml" ,
		   },


	  t_lower_tests  => {
                       input      => "NOTHING",
                       output      => "t_lower_tests.xml",
                       metric       => "../../tests/t_unit_tests/t_lower_tests_METRIC.xml" ,
                       controlfile  => "../../tests/t_unit_tests/t_lower_tests.xml" ,
		   },


#
#    free field test of the wilson inverter
#

	  t_propagator_w  => {
                       input      => "../../tests/propagator_w/INPUT_W.xml",
                       output      => "t_propagator_w.xml",
                       metric       => "../../tests/propagator_w/t_propagator_w_METRIC.xml" ,
                       controlfile  => "../../tests/propagator_w/t_propagator_w.xml" ,
		   },



#	  t_mesons_w    => {
#                       output       => "./t_mesons_w.xml",
#                       metric       => "control/metric_t_mesons_w.xml" ,
#                       controlfile  => "control/t_mesons_w.xml" ,
#		   },
#	 t_io    => {
#                       output       =>  "./t_io.xml",
#                       metric       =>  "control/metric_t_io.xml" ,
#                       controlfile  =>  "control/t_io.xml" ,
#		   },
	 ) ; 


$role = "output" ; 



print "Test   conclusion\n" ; 
# run the tests
foreach $execute ( keys %HoH ) 
{
    print "$execute \t\t";
    $candidate =  $HoH{$execute}{"output"} ; 
    $metric    =  $HoH{$execute}{"metric"} ; 
    $control   =  $HoH{$execute}{"controlfile"} ; 
    $input   =  $HoH{$execute}{"input"} ; 

##    if( ! -f  $execute )
##    {
	$log_make  = "log_make_".$execute  ; 
	system("make $execute >& $log_make ") ; 
##    }


    if(  -f  $execute )
    {
	if( $input ne "NOTHING" )
        {
	    $in_arg = "-i ".$input ; 
	}
	else
        {
	    $in_arg = "" ; 
        }
	

	$log = "log_".$execute  ; 
	$status_tmp = system("$execute ".$in_arg." -o $candidate >& $log") / 256 ; 
	if( $status_tmp != 0  ) 
        {
		print "   RUN_FAIL\n"  ; 
	}
	else
	{
	    $log_xml = "log_xml_".$execute  ; 

	    $status_tmp = system("$xmldiff $control $candidate $metric $log_xml ") ; 

	    $status = $status_tmp   ;   ## some perl feature

	    if( $status == 0 ) 
	    {
		print "   PASS\n"  ; 
	    }
	    else
	    {
		print "   FAIL\n"  ; 
	    }
	}
    }
    else
    {
	printf("   FAIL (compile)\n"); 

    }



}
