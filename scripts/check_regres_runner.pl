#!/usr/bin/perl
#
#  $Id: check_regres_runner.pl,v 1.4 2008-11-12 15:31:59 edwards Exp $
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
#   xmldiff is in a fixed location.
#   Perhaps, this should be incorporated in the autoconf
#   tool chain.

die "Usage: $0   top_srcdir  top_builddir xmldiff\n" unless scalar(@ARGV) == 3;

use File::Basename;

die "$top_srcdir does not exist" unless -d $ARGV[0];
die "$top_builddir does not exist" unless -d $ARGV[1];


$top_srcdir = &abs_path($ARGV[0]);
$top_builddir = &abs_path($ARGV[1]);
$xmldiff = $ARGV[2];

$test_dir = "$top_srcdir/tests";
$regres_dir = "$top_builddir/regres";

#
# Clear out the regression test dir
#
use File::Path;

# The list of regression dirs comes from here
require "$test_dir/regres.pl";

open(XMLOUT, ">$top_builddir/regres_report.xml");
print XMLOUT "<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n";
print XMLOUT "<testResults>\n";

$num_tests = 0;
$num_fails = 0;
$num_successes = 0;
for $file (&regresDirs())
{
    unless ($return = do "$file")
    {
	die "could not parse file $top_regres: $@" if $@;
	die "could not do $top_regres: $!"  unless defined $return;
	die "could not run $top_regres"     unless $return;
    }


    #
    # Create Load Leveler Script for the tests
    #
    for $h (@regres_list) 
    {
	my($exec_path) =  $h->{"exec_path"} ; 
	my($execute)   =  $h->{"execute"} ; 
	my($candidate);
	my($outdir);
	my($out_file);
	my($log_file);
        my($outdir);

	if ( $h->{"output"} ne "" && $h->{"log"} ne "" )
        {
	    die "Found both an output and a log entry - need only 1\n";
        }

        if ( $h->{"output"} ne "" )
        {
	    $candidate =  $h->{"output"} ; 
	    $out_file  =  $candidate;
	    $log_file  =  "XMLLOG";
        }
        elsif ( $h->{"log"} ne "" )
        {
	    $candidate =  $h->{"log"} ; 
	    $out_file  =  "XMLDAT";
	    $log_file  =  $candidate;
        }
        else
        {
	    die "Did not find either an  output or a log entry\n";
        }

        if ( $h->{"output_dir"} ne "" )
        {
          $outdir    =  $h->{"output_dir"};
        }
        else
        {
          $outdir    =  $candidate;
        }

	my($metric)    =  $h->{"metric"} ; 
	my($control)   =  $h->{"controlfile"} ; 
	my($input)     =  $h->{"input"} ; 

	my $canddir = "$regres_dir/$outdir";
	my($error_string) = "";

	printf "Checking: $candidate";

	chdir($canddir) || die "error cd $canddir : $!\n";


	my($log_xml) = "$canddir/xmldiff.log";
	my($xml_exe) = "$xmldiff $control $canddir/$candidate $metric $log_xml"; 
#	    print $xml_exe;
	my($status_tmp) = system("$xml_exe") ; 
	my($status) = $status_tmp   ;   ## some perl feature
	
	if( $status == 0 ) 
	{
	    $error_string = "PASS"; 
	    ++$num_successes;
	}
	else
	{
	    $error_string = "FAIL";
	    ++$num_fails;
	}

	++$num_tests;

	print "\t".$error_string."\n";
	print XMLOUT "<test program=\'$execute\' candidate=\'$candidate\'>";
	print XMLOUT $error_string;
	print XMLOUT "</test>\n";
    }

}
print XMLOUT "<summary tests=\'".$num_tests."\' successes=\'".$num_successes."\' failures=\'".$num_fails."\' />\n";
print XMLOUT "</testResults>\n";
close XMLOUT;


sub abs_path
{
    my($dir) = @_;
    die "directory  $dir  does not exit\n" unless -d $dir;
    my($abs) = `cd $dir ; pwd`; 
    chomp $abs;
    return $abs;
}
