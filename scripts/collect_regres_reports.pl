#!/usr/bin/perl
#
# Collect regres reports from jlab-standard-chroma-build
#

$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

die "Usage: $0  <top-level nightly or tagged dir>  <target directory>\n" unless scalar(@ARGV) == 2;

$top_dir    = $ARGV[0];
$target_dir = $ARGV[1];

die "Top-level build directory does not exist\n" unless -d $top_dir;
die "Target directory\n" unless -d $target_dir;

chdir($top_dir);

use File::Basename;

open(BUILD, "find $top_dir -name regres_report.html -print |");
open(INDEX, "> $target_dir/index.html");

print INDEX <<EOF;
<HTML>
<HEAD><TITLE>Regression Reports</TITLE></HEAD>

<H1 ALIGN="center">Regression Reports</H1>

EOF

while(<BUILD>)
{
  $input_file = $_;
  chomp $input_file;

  $output_file = basename($input_file);
  $basedir = dirname($input_file);

  printf "file=$input_file\n";

  $final_dir = "$target_dir/$basedir";

  system("mkdir -p $final_dir; /bin/cp -f $file $final_dir");

  # Highly annoying - I need to fix the stupid <meta> line since 
  # XSL spits out something that is not valid XML!
  open(IN, "< $input_file");
  open(OUT, "> $basedir/$output_file");
  while(<IN>)
  {
    chomp;
    $_ =~ s/><\/head>/\/><\/head>/;

    print OUT, $_;
  }
  
  close(IN);
  close(OUT);

  print INDEX <<EOF;
<A HREF="http://lqcd.jlab.org/~edwards/$target_dir/$basedir/$output_file">
${file}</a>
EOF

}

print INDEX <<EOF;
<P>
<FONT SIZE="-1">Last modified: `date`</FONT>

</BODY>
</HTML>
EOF

close(INDEX);
close(BUILD);

exit(0);
