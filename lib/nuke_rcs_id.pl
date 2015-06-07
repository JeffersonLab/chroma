#!/usr/bin/perl

open(FF, 'find . -name "*.{h,cc}" -print|') || die "error opening files\n";

while($_ = <FF>)
{
  chomp; my $f = $_;

  system("cat $f | grep -v '^\/\/ \$Id: ' > ff ; mv ff $f");
}
