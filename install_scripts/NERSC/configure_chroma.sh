#!/usr/bin/bash

ROOT=/u1/bjoo/Devel/SciDAC/install
CHROMA_HOME=$ROOT/chroma_wilson_single
QMP_HOME=$ROOT/qmp_generic
QDP_HOME=$ROOT/qdp_parscalar_wilson_1.9.1
BAGEL_HOME=$ROOT/bagel_single

./configure CXXFLAGS="-qrtti=all -I$BAGEL_HOME -I$QMP_HOME/include" \
            CFLAGS="-I$BAGEL_HOME -I$QMP_HOME/include" \
            LDFLAGS="-L$BAGEL_HOME -L$QMP_HOME/lib"  \
            LIBS="-lwfmpowerIIIs -lqmp" \
	    --with-qdp=$QDP_HOME \
	    --enable-pab-wilson-dslash=noarch \
	    --host=powerpc-ibm-aix --build=none \
	    --prefix=$CHROMA_HOME
