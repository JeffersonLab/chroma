nobase_include_HEADERS=sse_align.h \
	sse_config.h \
	sse_dslash.h \
	sse_dslash_qdp_packer.h \
	sse_dslash_qdp_packer_3d.h \
	sse_dslash_3d.h

if SCALAR_SINGLE
nobase_include_HEADERS += types32.h shift_tables_scalar.h dispatch_scalar.h site_dslash_32bit_scalar.h
endif

if SCALAR_DOUBLE
nobase_include_HEADERS += types64.h shift_tables_scalar.h dispatch_scalar.h site_dslash_64bit_scalar.h
endif

if PARSCALAR_SINGLE
nobase_include_HEADERS += types32.h shift_tables_parscalar.h dispatch_parscalar.h mvv_recons_32bit.h decomp.h decomp_hvv.h recons.h 
endif

if PARSCALAR_DOUBLE
nobase_include_HEADERS += types64.h shift_tables_parscalar.h  dispatch_parscalar.h mvv_recons_64bit.h decomp.h decomp_hvv.h recons.h
endif


nodist_include_HEADERS=sse_config_internal.h
