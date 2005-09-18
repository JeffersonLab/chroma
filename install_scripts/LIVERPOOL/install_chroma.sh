#!/usr/local/bin/bash
#
# Create an installed and tagged chroma/qdp++ installation
# on a Liverpool machine.
#
# This is a script to automate pulling
# the chroma/qdp++ package from JLAB
# and to compile it.
#
# See:
# http://www.jlab.org/~edwards/lattice_cvs.html
#
# The "cvs login" command only needs
# to be executed once, because the password
# is stored in: ~/.cvspass
#
#
#
export wot_tag="LIVTAG_1.00"


export CVSROOT=:pserver:anonymous@cvs.jlab.org:/group/lattice/cvsroot
export ARCH="scalar"

#
#
echo "Creating CHROMA distribution: "${wot_tag}

export up=`cd .. ; pwd`
echo $up
qdp_install=$up"/chroma_install/"${ARCH}_${wot_tag}


#
#  compile in a separate directory
#
export src_dir="src_"${wot_tag}

if test -d $src_dir
then
 echo "Error $src_dir already exists"
 exit 1
else
 mkdir -p $src_dir
fi


cd ${src_dir}

here=`pwd`
export log=${here}"/log_compilation"


if test -d $qdp_install 
then
 echo "Error $qdp_install already exists"
 exit 1
else
 mkdir -p $qdp_install
fi


#
# get & compile qdp++
#
cvs co qdp++  >> ${log} 2>&1
cd qdp++
./configure --enable-parallel-arch=${ARCH}  --enable-Ns=4 --enable-Nc=3 --prefix=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000" >> ${log} 2>&1
make >> ${log} 2>&1
make install  >> ${log} 2>&1
cd ..


#
# get chroma
#

cvs co chroma

cd chroma
aclocal

aclocal
automake
./configure  --prefix=${qdp_install}  --with-qdp=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000" 

make

make install

cd ..

#
#  write some simple test information
#

cd ${qdp_install}

info="compilation_information"

echo "Compilation information"   >> ${info}
echo "-----------------------"   >> ${info}
date >> ${info}
machine=`uname -a`
echo "Machine: $machine"  >> ${info}
compiler=`whoami`
echo "Compiler: $compiler"  >> ${info}
echo "Notes"  >> ${info}

exit 0 


