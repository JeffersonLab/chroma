#!/usr/bin/perl
#
#  $Id: run_chroma_xmldiff.pl,v 1.10 2006-01-03 04:02:44 bjoo Exp $
#
#  This is wrapper script to run the xmldiff application from
#  a makefile
#
#  Homepage for xmldiff application
#  https://forge.nesc.ac.uk/projects/xmldiff/

#
# More work:
#   At the moment this script assumes that the 
#   xmldiff is in a fixed location.
#   Perhaps, this should be incorporated in the autoconf
#   tool chain.

die "Usage: run_chroma_xmldiff.pl  top_srcdir  top_builddir\n" unless scalar(@ARGV) == 2;

$top_srcdir = &abs_path($ARGV[0]);
$top_builddir = &abs_path($ARGV[1]);
die "$top_srcdir does not exist" unless -d $top_srcdir;
die "$top_builddir does not exist" unless -d $top_builddir;

printf "Source directory = %s\n", $top_srcdir;
printf "Build directory = %s\n", $top_builddir;


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
    printf "\n%-15s %-40s         %s\n", "Program", "Candidate","Conclusion";
    for $h (@regres_list) 
    {
	my($exec_path) =  $h->{"exec_path"} ; 
	my($execute)   =  $h->{"execute"} ; 
	my($candidate) =  $h->{"output"} ; 
	my($outdir)    =  $h->{"output"} ; 
	my($metric)    =  $h->{"metric"} ; 
	my($control)   =  $h->{"controlfile"} ; 
	my($input)     =  $h->{"input"} ; 

	my $canddir = "$regres_dir/$outdir";

#       printf "exec_path=$exec_path\n";
#       printf "canddir=$canddir\n";

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
		++$num_errors;
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
		    ++$num_errors;
		}
	    }
	}
	else
	{
	    printf("FAIL (compile)\n"); 
	    ++$num_errors;
	}

    }
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
