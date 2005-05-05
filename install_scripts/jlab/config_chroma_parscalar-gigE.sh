#!/bin/sh

#../configure --prefix=/usr/local/chroma/parscalar-gigE -with-qdp=/usr/local/qdp++/parscalar-gigE CXXFLAGS=-static --enable-sse-wilson-dslash  --enable-sse-dwf-cg  --enable-gmp  LIBS="-lgmp"

../configure --prefix=/usr/local/chroma/parscalar-gigE -with-qdp=/usr/local/qdp++/parscalar-gigE CXXFLAGS=-static --enable-sse-wilson-dslash   --enable-gmp  LIBS="-lgmp"
