#!/usr/bin/perl
#
#  $Id: run_chroma_xmldiff.pl,v 1.1 2005-12-17 18:51:23 edwards Exp $
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

printf "top_srcdir = %s\n", $top_srcdir;
printf "top_builddir = %s\n", $top_builddir;


# location of xmldiff
$xmldiff = "/usr/local/bin/xmldiff" ;
#$xmldiff = "$top_srcdir/scripts/xmldiff" ;

if( ! -x $xmldiff)
{
    print "Error:".$0." needs the xmldiff utility in your path\n" ; 
    print "Download it from http://forge.nesc.ac.uk/projects/xmldiff/ \n" ; 
    exit(1) ;
}

# at some stage this should be factored into another perl script

$test_dir = "$top_srcdir/tests";
$regres_dir = "$top_builddir/regres";

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#

%HoH = (
	1 => {
	    exec_path   => "$top_builddir/mainprogs/main" , 
	    execute     => "chroma" , 
	    input       => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.ini.xml" , 
	    output      => "prec_wilson-v9.candidate.xml",
	    metric      => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.metric.xml" ,
	    controlfile => "$test_dir/chroma/hadron/propagator/prec_wilson-v9.out.xml" ,
	},
	) ; 


# Clear out the regression test dir
use File::Path;
if ( -d $regres_dir )
{
    rmtree([$regres_dir]);
}
printf "regres=$regres_dir\n";
mkpath([$regres_dir], 1, 0755);
    

# run the tests
print "Test   conclusion\n" ; 
foreach $a ( keys %HoH ) 
{
    &xmldiff($regres_dir,*HoH,$a);
}


exit(0);

sub xmldiff
{
    local($rundir,*obj,$k) = @_;

    my($exec_path) =  $obj{$k}{"exec_path"} ; 
    my($execute)   =  $obj{$k}{"execute"} ; 
    my($candidate) =  $obj{$k}{"output"} ; 
    my($outdir)    =  $obj{$k}{"output"} ; 
    my($metric)    =  $obj{$k}{"metric"} ; 
    my($control)   =  $obj{$k}{"controlfile"} ; 
    my($input)     =  $obj{$k}{"input"} ; 

    use File::Path;

    my $canddir = "$rundir/$outdir";
    if (-d $canddir)
    {
	rmtree([$canddir]);
    }
    mkpath([$canddir], 1, 0755);
    chdir($canddir) || die "error cd $canddir : $!\n";

    print "$execute \t\t";
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
	$exe = "$exec_path/$execute ".$in_arg." -o $candidate 2>${log}.err > ${log}.out"; 
#	print $exe;
	my($status_tmp) = system("$exe") / 256 ; 
	if( $status_tmp != 0  ) 
        {
	    print "   RUN_FAIL\n"  ; 
	}
	else
	{
	    my($log_xml) = "$canddir/xmldiff.log";

	    my($status_tmp) = system("$xmldiff $control $candidate $metric $log_xml") ; 
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
	printf("   FAIL (compile)\n"); 

    }


}
