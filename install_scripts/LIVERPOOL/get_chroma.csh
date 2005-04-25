#!/usr/local/bin/bash
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
# This is the script I give people to get both
# qdp++/chroma. This is the version we use on the 
# linux boxes. Don't use it on a SUN. Only old people
# use suns, not smart young things.
#
# To set up the correct version of autoconf/automake
# and pickup the correct version of xmldiff 
# you should source /RAID/lattice_qcdlib/setup_qcd.sh
#
# where RAID is the main raid disk for that machine.
# Fo example, on rusty RAID=raidr and on iasc RAID=raidz
#
#
#

export CVSROOT=:pserver:anonymous@cvs.jlab.org:/group/lattice/cvsroot

##qdp_install="/raidz/mcneile/projects/large_Nc/qdp++_install"
here=`pwd`
qdp_install=${here}"/qdp++_install"

# only needs to be done once
##cvs login

#
# get & compile qdp++
#
cvs co qdp++
cd qdp++
./configure --enable-parallel-arch=scalar  --enable-Ns=4 --enable-Nc=3 --prefix=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000"
make
make install 
cd ..


#
# get chroma
#

cvs co chroma

cd chroma

##./configure --with-qdp=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000" --enable-fermion-type=staggered

##./configure --with-qdp=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000" --enable-fermion-type=staggered

./configure --with-qdp=${qdp_install} CXXFLAGS=" -O2 -finline-limit=50000" 

make


cd ..


exit 0 


