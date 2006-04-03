# -*- perl -*-
# $Id: link.pl,v 3.0 2006-04-03 04:59:22 edwards Exp $

# This is supposed to do the same as a 'make check' in the build directory
# except the code is compiled with the installed header files and linked  
# with the installed libraries.
# To find the code to build, it looks in directories in the source tree 
# called 'examples' for a 'Makefile.am', from which it gets the targets from
# the 'check_PROGRAMS' variable. 

use File::Find;

$#ARGV==-1 and die "Usage: $0 module";

$module = $ARGV[0];

$srcdir = "../../$module" ; 

# Find directories that contain a Makefile.am with check_PROGRAMS 
# or bin_PROGRAMS defined

File::Find::find({wanted => \&wanted}, $srcdir);

sub wanted {

    # The Chroma other_libs programs tend to be specialised, so let's
    # ignore them..

    $module eq "chroma" && ("$File::Find::dir" =~ /other_libs/) && return;
    

    /Makefile\.am$/ || return;
    (open MFAM, "<$ENV{PWD}/$File::Find::name") || warn "Cannot open $File::Find::name $!\n";

    while(<MFAM>){
	/^\s*(check|bin)_PROGRAMS\s*=/ || next;
 	push @srcexdirs, "$File::Find::dir";
	last;
    }
    close MFAM;
}


$srcdir =~ s|\+|\\+|g;           # '+' needs to be escaped for the regexp

# Loop over all suitable named directories in the module source tree

$status = 0;
for $srcexdir (@srcexdirs){

    ($linkdir = $srcexdir) =~ s|$srcdir/||; 



# Create a simple makefile in a build directory. This depends on the module
# so we give this job to a local script

    $msg = `bash ../../makemake.sh $srcexdir $linkdir`;
    if($? != 0) {
	warn "Failed to make Makefile in $linkdir. $msg\n";
	next;
    }
    
    
    # Get the build targets from the original Makefile

    $makefile = "$srcexdir/Makefile.am";
    open MAKEFILE, "<$makefile" 
	or warn "Cannot read $makefile. $!\n" and next;

    $in_targets = 0;
    while (<MAKEFILE>){
	if(/^\s*(check|bin)_PROGRAMS\s*=/){
	    $in_targets = 1;
	    @targets = split;
	    shift @targets; shift @targets;
	    $targets[$#targets] eq "\\" && pop @targets or last;
	    next
	}
	$in_targets or next;
	push @targets, split;
	/\\\s*$/ and pop @targets or $in_targets = 0;
    }

    close MAKEFILE;

    # Build the targets

    for $target (@targets){
	@args = ("gmake", "-C", "$linkdir", "BIN=$target");
	$status += system(@args) ; 
    }


}

exit $status % 255;
