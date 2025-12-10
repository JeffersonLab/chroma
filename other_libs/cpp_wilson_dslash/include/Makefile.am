nobase_include_HEADERS= allocate.h \
	cpp_dslash.h \
	cpp_clover.h \
	cache.h \
	dslash_config.h \
	sse_sign_32bit.h \
	sse_sign_64bit.h \
	cpp_dslash_types.h \
	cpp_clover_types.h


if SCALAR
nobase_include_HEADERS += \
	cpp_dslash_scalar.h \
	dispatch_scalar.h \
	shift_table_scalar.h \
	shift_table_3d_scalar.h\
	cpp_dslash_scalar_32bit.h \
	cpp_dslash_scalar_64bit.h \
	cpp_dslash_scalar_32bit_sse.h \
	cpp_dslash_scalar_32bit_c.h \
	cpp_dslash_scalar_64bit_sse.h \
	cpp_dslash_scalar_64bit_c.h \
	cpp_clover_scalar.h \
	cpp_clover_site_apply_32bit.h \
	cpp_clover_site_apply_64bit.h
endif

if PARSCALAR
nobase_include_HEADERS += cpp_dslash_parscalar.h \
	dispatch_parscalar.h \
	shift_table_parscalar_types.h \
	shift_table_parscalar.h \
	shift_table_3d_parscalar.h \
	tables_parscalar.h \
	cpp_dslash_parscalar_utils_32bit.h \
	cpp_dslash_parscalar_utils_64bit.h \
	cpp_dslash_matvec32bit_c.h \
	cpp_dslash_matvec64bit_c.h \
	cpp_dslash_parscalar_decomp_32bit_c.h \
	cpp_dslash_parscalar_decomp_hvv_32bit_c.h \
	cpp_dslash_parscalar_mvv_recons_32bit_c.h \
	cpp_dslash_parscalar_recons_32bit_c.h
	cpp_dslash_parscalar_decomp_64bit_c.h \
	cpp_dslash_parscalar_decomp_hvv_64bit_c.h \
	cpp_dslash_parscalar_mvv_recons_64bit_c.h \
	cpp_dslash_parscalar_recons_64bit_c.h \
	cpp_dslash_parscalar_decomp_32bit_sse2.h \
	cpp_dslash_parscalar_decomp_hvv_32bit_sse2.h \
	cpp_dslash_parscalar_mvv_recons_32bit_sse2.h \
	cpp_dslash_parscalar_recons_32bit_sse2.h
	cpp_dslash_parscalar_decomp_64bit_sse2.h \
	cpp_dslash_parscalar_decomp_hvv_64bit_sse2.h \
	cpp_dslash_parscalar_mvv_recons_64bit_sse2.h \
	cpp_dslash_parscalar_recons_64bit_sse2.h \
	cpp_clover_parscalar.h \
	cpp_clover_site_apply_32bit.h \
	cpp_clover_site_apply_32bit_c.h \
	cpp_clover_site_apply_32bit_sse.h \
	cpp_clover_site_apply_64bit.h \
	cpp_clover_site_apply_64bit_c.h \
	cpp_clover_site_apply_64bit_sse.h \
	cpp_clover_odd_odd_parscalar_32bit.h \
	cpp_clover_odd_odd_parscalar_64bit.h

endif

if BUILD_QDP_PACKERS
nobase_include_HEADERS += cpp_dslash_qdp_packer.h
endif
nodist_include_HEADERS=dslash_config_internal.h
