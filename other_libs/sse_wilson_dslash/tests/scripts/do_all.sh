make distclean 
../sse_wilson_dslash/configure  \
	--with-qdp=/home/bjoo/Devel/QCD/jlab-standard-chroma-build/install/qdp++/qdp1-25-1/parscalar-ib-ofed1-ompi \
	CXXFLAGS="-g" CFLAGS="-g -O3"
make clean; make ; make check
./test_all.sh
