#!/bin/tcsh

../configure --prefix=/usr/local/share/chroma/scalar --with-qdp=/usr/local/share/qdp++/scalar CXXFLAGS="" LIBS="-lgmp" --enable-sse-wilson-dslash  --enable-gmp
