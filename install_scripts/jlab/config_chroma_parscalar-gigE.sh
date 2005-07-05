#!/bin/sh

#../configure --prefix=/usr/local/chroma/parscalar-gigE -with-qdp=/usr/local/qdp++/parscalar-gigE CXXFLAGS=-static --enable-sse-wilson-dslash  --enable-sse-dwf-cg  --with-gmp=/usr

../configure --prefix=/usr/local/chroma/parscalar-gigE -with-qdp=/usr/local/qdp++/parscalar-gigE CXXFLAGS=-static --enable-sse-wilson-dslash   --enable-opt-cfz-linop  --with-gmp=/usr
