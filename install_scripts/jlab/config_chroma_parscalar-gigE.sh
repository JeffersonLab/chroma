#!/bin/sh

../configure --prefix=/usr/local/chroma/parscalar-gigE -with-qdp=/usr/local/qdp++/parscalar-gigE CXXFLAGS=-static --enable-sse-wilson-dslash 
