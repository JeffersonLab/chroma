#!/bin/sh

../configure --prefix=/usr/local/chroma/parscalar-gm --with-qdp=/usr/local/qdp++/parscalar-gm CXXFLAGS="-static" --enable-sse-wilson-dslash --enable-sse-dwf-cg

