#!/usr/bin/perl
#
#  $Id: run_chroma_xmldiff.pl,v 1.2 2005-12-18 21:06:52 edwards Exp $
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

die "Usage: run_chroma_xmldiff.pl  top_srcdir  top_builddir\n" unless scalar(@ARGV) == 2;

$top_srcdir = $ARGV[0];
$top_builddir = $ARGV[1];
die "$top_srcdir does not exist" unless -d $top_srcdir;
die "$top_builddir does not exist" unless -d $top_builddir;

printf "Source directory = %s\n", $top_srcdir;
printf "Build directory = %s\n", $top_builddir;


# location of xmldiff
$xmldiff = "/usr/local/bin/xmldiff" ;
#$xmldiff = "$top_srcdir/scripts/xmldiff" ;

if( ! -x $xmldiff)
{
    print "Error:".$0." needs the xmldiff utility in your path\n" ; 
    print "Download it from http://forge.nesc.ac.uk/projects/xmldiff/ \n" ; 
    exit(1) ;
}

$test_dir = "$top_srcdir/tests";
$regres_dir = "$top_builddir/regres";

printf "Regression test directory = %s\n", $regres_dir;

# at some stage this should be factored into another perl script

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#

@HoH = (
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.ini.xml" , 
	    output      => "prec_wilson-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.out.xml" ,
	},
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.ini.xml" , 
	    output      => "prec_wilson-twisted-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/prec_wilson-twisted-v9.out.xml" ,
	},
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.ini.xml" , 
	    output      => "unprec-ovlap-zolo4d-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/unprec-ovlap-zolo4d-v9.out.xml" ,
	},
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.ini.xml" , 
	    output      => "unprec_hamberwu-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/unprec_hamberwu-v9.out.xml" ,
	},
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.ini.xml" , 
	    output      => "prec_dwf-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/prec_dwf-v9.out.xml" ,
	},
	{
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.ini.xml" , 
	    output      => "prec_parwilson-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/prec_parwilson-v9.out.xml" ,
	},
	); 


# Clear out the regression test dir
use File::Path;
if ( -d $regres_dir )
{
    rmtree([$regres_dir]);
}
#printf "regres=$regres_dir\n";
mkpath([$regres_dir], 0, 0755);
    

# run the tests
printf "\n%-15s %-40s         %s\n", "Program", "Candidate","Conclusion";
for $h (@HoH) 
{
    my($exec_path) =  $h->{"exec_path"} ; 
    my($execute)   =  $h->{"execute"} ; 
    my($candidate) =  $h->{"output"} ; 
    my($outdir)    =  $h->{"output"} ; 
    my($metric)    =  $h->{"metric"} ; 
    my($control)   =  $h->{"controlfile"} ; 
    my($input)     =  $h->{"input"} ; 

    my $canddir = "$regres_dir/$outdir";

#    printf "exec_path=$exec_path\n";
#    printf "canddir=$canddir\n";

    if (-d $canddir)
    {
	rmtree([$canddir]);
    }
    mkpath([$canddir], 0, 0755);
    chdir($canddir) || die "error cd $canddir : $!\n";

    printf "%-15s %-40s      ", $execute, $candidate;
    if (-x "$exec_path/$execute")
    {
	if( $input ne "NOTHING" )
        {
	    $in_arg = "-i ".$input ; 
	}
	else
        {
	    $in_arg = "" ; 
        }
	

	my($log) = "$canddir/$execute"  ; 
	my($exe) = "$exec_path/$execute ".$in_arg." -o $candidate 2>${log}.err > ${log}.out"; 
#	print $exe;
	my($status_tmp) = system("$exe") / 256 ; 
	if( $status_tmp != 0  ) 
        {
	    print "   RUN_FAIL\n"  ; 
	}
	else
	{
	    my($log_xml) = "$canddir/xmldiff.log";

	    my($xml_exe) = "$xmldiff $control $canddir/$candidate $metric $log_xml"; 
#	    print $xml_exe;
	    my($status_tmp) = system("$xml_exe") ; 
	    my($status) = $status_tmp   ;   ## some perl feature

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
	printf("FAIL (compile)\n"); 

    }


}
