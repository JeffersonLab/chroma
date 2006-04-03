# $Id: makemake_qdp.sh,v 3.0 2006-04-03 04:59:22 edwards Exp $

# Write a simple makefile using the module's config script.

if [ $# -ne 2 ]
then 
    echo "Usage: $0 src-directory build-directory"
    exit 1
fi

module=qdp++
srcdir=$1
builddir=$2

config=../install/bin/$module-config
if [ ! -f $config ] 
then
    echo "Cannot find $config"
    exit 1
fi

if [ ! -d $srcdir ]
then 
    echo "Cannot find source directory $srcdir"
    exit 1;
fi

if [ -d $builddir ] 
then
    rm -r $builddir
fi

mkdir -p $builddir


# Some of the test sources are not stand-alone, in the sense that they need
# to be built with some other auxiliary files. Rather than try and figure
# out what needs what, I'll just compile all of them - apart from junk.cc

for f in `find $srcdir -name "*.cc" -o -name "*.c"`
do
    [ `grep -c "main *(" $f` -ne 0 ] && continue       # Ignore 'main' files
    f=`basename $f`
    [ "$f" == "junk.cc" ] && continue                  # Ignore this file!
    [ "$f" == "linalg2.c" ] && continue                # and this one!
    auxsrc="$auxsrc `basename $f`"
done

# In the QDP libraries we have C and C++ code. The C will be compiled
# with the implicit $(CC), but, since we must explicitly define the linking
# rule (because we're linking with libraries), the linking is done with the
# C++ compiler. Since this is just going to hand over the objects to the 
# loader, we should be all right.

echo CXX = `$config --cxx`          >  $builddir/Makefile
echo LDFLAGS = `$config --ldflags`>> $builddir/Makefile
echo CXXFLAGS = `$config --cxxflags`>> $builddir/Makefile
echo CFLAGS = `$config --cxxflags`>> $builddir/Makefile
echo LIBS = `$config --libs`      >> $builddir/Makefile

echo "VPATH = $PWD/$srcdir" 	  >> $builddir/Makefile
echo "AUX := \$(basename $auxsrc)" >> $builddir/Makefile
echo 'AUX := $(addsuffix .o, $(AUX))' >> $builddir/Makefile
echo 'OBJ = $(BIN).o'  	  >> $builddir/Makefile
echo '${BIN}: ${OBJ} $(AUX)'      	  >> $builddir/Makefile
printf "\t\${CXX} -o \$@ \${OBJ} \$(AUX) \${LDFLAGS} \${LIBS}\n" >> $builddir/Makefile



