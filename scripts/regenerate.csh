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

set tmp_status=$status
if ($tmp_status != 0) then
  echo "Error in $0, found error returned by $exe in regression input $input"  
  exit 1
endif
  
end
