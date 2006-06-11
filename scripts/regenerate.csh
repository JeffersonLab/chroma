#!/bin/tcsh

foreach f (`find . -name "*.metric.xml" -print`)
set ini=`echo $f | sed 's/\.metric\.xml/.ini.xml/'`
set out=`echo $f | sed 's/\.metric\.xml/.out.xml/'`
echo -i $ini -o $out 
/home/edwards/qcd/chroma/scalar/mainprogs/main/chroma.exe -i $ini -o $out > output
end
