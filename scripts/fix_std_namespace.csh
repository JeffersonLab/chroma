#!/bin/tcsh

set f=$1

echo $f

cat $f |\
sed 's/ cout/ std::cout/g' |\
sed 's/ cerr/ std::cerr/g' |\
sed 's/string/std::string/g' |\
sed 's/flush/std::flush/g' |\
sed 's/ostd::string/ostring/g' |\
sed 's/istd::string/istring/g' |\
sed 's/ostringstream/std::ostringstream/g' |\
sed 's/istringstream/std::istringstream/g' |\
sed 's/endl/std::endl/g' |\
sed 's/map/std::map/g' |\
sed 's/vector/std::vector/g' |\
sed 's/setw/std::setw/g' |\
sed 's/setfill/std::setfill/g' |\
sed 's/std::std/std/g' |\
sed 's/std::string_in/string_in/g' |\
sed 's/std::string_out/string_out/g' |\
sed 's/include <std::string>/include <string>/' |\
sed 's/include <std::string.h>/include <string.h>/' |\
sed 's/include <std::map>/include <map>/' |\
sed 's/include <std::vector>/include <vector>/' |\
sed 's/error_std::string/error_string/g' |\
sed 's/_std::map/_map/g' |\
sed 's/_std::map/_map/g' |\
sed 's/_objstd::map/_objmap/g' |\
sed 's/std::vector_/vector_/g' |\
sed 's/_std::vector/_vector/g' |\
sed 's/_std::vector/_vector/g' |\
grep -v 'using namespace std;' > ${f}.tmp

/bin/mv ${f}.tmp $f
