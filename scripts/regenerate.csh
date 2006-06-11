#!/bin/tcsh

set exe = $1
set dir = $2

foreach f (`find $dir -name "*.metric.xml" -print`)
set ini=`echo $f | sed 's/\.metric\.xml/.ini.xml/'`
set out=`echo $f | sed 's/\.metric\.xml/.out.xml/'`
echo $exe -i $ini -o $out 
$exe  -i $ini -o $out > output
end
