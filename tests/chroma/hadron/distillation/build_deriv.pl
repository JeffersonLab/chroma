#!/usr/bin/perl

# always print empty list
my @foo = ();
printf "          <elem></elem>\n";

# 1-deriv
foreach $n (1 .. 3) {
  push(@foo, "$n");
}

&printlist(@foo);

# 2 up to 6 deriv
foreach my $N (2 .. 6) {
  &printlist(&recurs($N, @foo));
}
 
exit 0;

sub printlist
{
  my (@list) = @_;
  foreach my $f (@list) {
    printf "          <elem>${f}</elem>\n", $f;
  }
}
  
sub recurs
{
  my (@list) = @_;
  my @new_list = ();
  my $n = shift(@list);

  if ($n <= 1) {
    return @list;
  }

  foreach $ll (@list) {
    foreach $n (1 .. 3) {
      push(@new_list, "$ll $n");
    }
  }

  return &recurs($n - 1, @new_list);
}


