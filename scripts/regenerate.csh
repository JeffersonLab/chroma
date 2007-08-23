#!/bin/tcsh
# 
# Regenerate regression output within a directory tree.

if ($#argv != 2) then
  echo "Usage: $0  <executable>  <regression directory>"
  exit 1
endif

set exe = $1
set dir = $2

foreach f (`find $dir -name "*.metric.xml" -print`)
set ini=`echo $f | sed 's/\.metric\.xml/.ini.xml/'`
set out=`echo $f | sed 's/\.metric\.xml/.out.xml/'`
echo $exe -i $ini -o $out 
$exe  -i $ini -o $out
end
