#!/usr/bin/perl
#
#  $Id: create_regres_runner.pl,v 1.2 2007-10-31 14:14:38 edwards Exp $
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

die "Usage: $0  top_srcdir  top_builddir runner\n" unless scalar(@ARGV) == 3;

use File::Basename;

die "$top_srcdir does not exist" unless -d $ARGV[0];
die "$top_builddir does not exist" unless -d $ARGV[1];


$top_srcdir = &abs_path($ARGV[0]);
$top_builddir = &abs_path($ARGV[1]);
$run = $ARGV[2];

$test_dir = "$top_srcdir/tests";
$regres_dir = "$top_builddir/regres";

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

print "#! /bin/bash\n";
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
	my($candidate) =  $h->{"output"} ; 
	my($outdir)    =  $h->{"output"} ; 
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

	if( $input ne "NOTHING" )
	{
	    $in_arg = "-i ".$input ; 
	}
	else
	{
	    $in_arg = "" ; 
	}
	
	
	my($log) = "$canddir/$execute"  ; 
	my($exe) = "pushd $canddir ; echo Doing $candidate; $run $exec_path/$execute ".$in_arg." -o $candidate 2>${log}.err > ${log}.out; popd \n"; 
	print $exe;
    }
}

sub abs_path
{
    my($dir) = @_;
    die "directory  $dir  does not exit\n" unless -d $dir;
    my($abs) = `cd $dir ; pwd`; 
    chomp $abs;
    return $abs;
}
