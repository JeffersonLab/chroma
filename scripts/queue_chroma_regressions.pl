#!/usr/bin/perl
#
#  $Id: queue_chroma_regressions.pl,v 1.3 2007-10-31 14:14:38 edwards Exp $
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

die "Usage: $0  top_srcdir  top_builddir  queue_command\n" unless scalar(@ARGV) == 3;

use File::Basename;

die "$top_srcdir does not exist" unless -d $ARGV[0];
die "$top_builddir does not exist" unless -d $ARGV[1];
die "$run does not exist" unless -x $ARGV[2];

$top_srcdir = &abs_path($ARGV[0]);
$top_builddir = &abs_path($ARGV[1]);
$queue = &abs_path(dirname($ARGV[2])) . "/" . basename($ARGV[2]);
$start_time = time();
printf "\nStart date = %s\n\n", scalar(localtime($start_time));
printf "Source directory = %s\n", $top_srcdir;
printf "Build directory  = %s\n", $top_builddir;
printf "Run command      = %s\n", $queue;

$test_dir = "$top_srcdir/tests";
$regres_dir = "$top_builddir/regres";
printf "Regression test directory = %s\n", $regres_dir;

#
# Clear out the regression test dir
#
use File::Path;
if ( -d $regres_dir )
{
    rmtree([$regres_dir]);
}
#printf "regres=$regres_dir\n";
mkpath([$regres_dir], 0, 0755);

# The list of regression dirs comes from here
require "$test_dir/regres.pl";

printf "\nRun through regression test list\n";
$num_errors = 0;

for $file (&regresDirs())
{
    printf "\nSource $file\n";
    unless ($return = do "$file")
    {
	die "could not parse file $top_regres: $@" if $@;
	die "could not do $top_regres: $!"  unless defined $return;
	die "could not run $top_regres"     unless $return;
    }

    #
    # Run the tests
    #
    printf "\n%-15s  %-10s  %s\n", "Program", "Conclusion", "Candidate";
    for $h (@regres_list) 
    {
	my($exec_path) =  $h->{"exec_path"} ; 
	my($execute)   =  $h->{"execute"} ; 
	my($candidate);
	my($outdir);
	my($out_file);
	my($log_file);

        if ( $h->{"output"} ne "" && $h->{"log"} ne "" )
        {
          die "Found both an output and a log entry - need only 1\n";
        }

        if ( $h->{"output"} ne "" )
        {
	  $candidate =  $h->{"output"} ; 
	  $outdir    =  $h->{"output"} ; 
	  $out_file  =  $candidate;
	  $log_file  =  "XMLLOG";
        }
        elsif ( $h->{"log"} ne "" )
        {
	  $candidate =  $h->{"log"} ; 
	  $outdir    =  $h->{"log"} ; 
	  $out_file  =  "XMLDAT";
	  $log_file  =  $candidate;
        }
        else
        {
          die "Did not find either an  output or a log entry\n";
        }

	my($metric)    =  $h->{"metric"} ; 
	my($control)   =  $h->{"controlfile"} ; 
	my($input)     =  $h->{"input"} ; 

	my $canddir = "$regres_dir/$outdir";
	my($error_string) = "";

#       printf "exec_path=$exec_path\n";
#       printf "canddir=$canddir\n";

	if (-d $canddir)
	{
	    rmtree([$canddir]);
	}
	mkpath([$canddir], 0, 0755);
	chdir($canddir) || die "error cd $canddir : $!\n";

	printf "%-15s", $execute;
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
	    my($exe) = "$queue $exec_path/$execute ".$in_arg." -o $out_file -l $log_file < /dev/null 2>${log}.err > ${log}.out"; 
	    my($status_tmp) = system("$exe") / 256 ; 
	    if( $status_tmp != 0  ) 
	    {
		$error_string = "RUN_FAIL"; 
		++$num_errors;
	    }
	    else
	    {
		$error_string = "QUEUED";
	    }
	}
	else
	{
	    $error_string = "FAIL (compile)";
	    ++$num_errors;
	}

	printf "  %-10s  %s\n", $error_string, $candidate;
    }
}

$end_time = time();
printf "\nEnd date = %s\n", scalar(localtime($end_time));
$time_diff = $end_time - $start_time;
if ($time_diff > 1) 
{
  printf "Regression queue time = %d secs = %.1f min\n", $time_diff, $time_diff/60.0;
}

printf("\nTotal of $num_errors failures\n"); 
exit($num_errors);


sub abs_path
{
    my($dir) = @_;
    die "directory  $dir  does not exit\n" unless -d $dir;
    my($abs) = `cd $dir ; pwd`; 
    chomp $abs;
    return $abs;
}
