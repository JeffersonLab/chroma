#!/usr/bin/perl
#
# $Id: build_libchroma_sources.pl,v 3.0 2006-04-03 04:59:22 edwards Exp $
#
# Build the  libchroma_a_SOURCES  lines in  chroma/lib/Makefile.am
#
# Usage
#   cd chroma/lib
#   build_libchroma_sources.pl
#

print "libchroma_a_SOURCES =";
# All subdirs
@subdirs = ("actions", "io", "meas", "util");

# Go through each dir and find fermion independent sources
foreach $i (@subdirs)
{
  system("find $i -name \"*.cc\" -print | grep -v '_[ws].cc' | fmt -w 65 | awk '{printf \" \\\\\\n\\t%s\", \$0}'");
}

print "\n\n";

# Now print only the Wilson specific sources
print "# Wilson specific sources\n";
print "libchroma_a_SOURCES +=";

# Go through each dir and find Wilson specific headers
foreach $i (@subdirs)
{
  system("find $i -name \"*_w.cc\" -print | fmt -w 65 | awk '{printf \" \\\\\\n\\t%s\", \$0}'");
}

print "\n\n";

# Now print only the Staggered specific headers
print "# Staggered specific sources\n";
print "libchroma_a_SOURCES +=";

# Go through each dir and find Staggered specific sources
foreach $i (@subdirs)
{
  system("find $i -name \"*_s.cc\" -print | fmt -w 65 | awk '{printf \" \\\\\\n\\t%s\", \$0}'");
}

print "\n\n";


exit(0);
