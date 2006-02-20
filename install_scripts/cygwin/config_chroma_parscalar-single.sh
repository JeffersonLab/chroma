#!/bin/csh

../configure --prefix=/usr/local/share/chroma/parscalar-single --with-qdp=/usr/local/share/qdp++/parscalar-single --enable-sse-wilson-dslash CXXFLAGS="" LIBS="-lgmp"  --enable-sse-dwf-cg  --with-gmp=/usr
