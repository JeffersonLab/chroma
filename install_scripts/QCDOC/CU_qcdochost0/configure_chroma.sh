../../chroma/configure CXXFLAGS="" \
	        CFLAGS="" \
                LDFLAGS="" \
                LIBS="" \
                --with-qdp=/home/bj/install/QOS_QMP_v2/qdp_parscalar_single_ddr \
		--enable-pab-wilson-dslash=noarch \
	        --host=powerpc-gnu-elf --build=sparc-sun-solaris2.9 \
		--prefix=/home/bj/install/QOS_QMP_v2/chroma_wilson_single_ddr
