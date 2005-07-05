#!/bin/sh

../configure --prefix=/usr/local/chroma/parscalar-single --with-qdp=/usr/local/qdp++/parscalar-single CXXFLAGS="-static" --enable-sse-wilson-dslash   --enable-opt-cfz-linop  --enable-sse-dwf-cg  --with-gmp=/usr

