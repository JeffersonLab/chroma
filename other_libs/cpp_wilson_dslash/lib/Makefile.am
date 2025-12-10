# $Id: Makefile.am,v 1.4 2009-07-17 19:14:47 bjoo Exp $
#
# Path for all files used in building the library
# This file was automatically generated

TOPSRCDIR=@top_srcdir@
TOPBUILDDIR=@top_builddir@
#
# Include Path
#
INCFLAGS = -I$(TOPSRCDIR)/include -I$(TOPBUILDDIR)/include
#
# Any other compiler flags
# 
AM_CFLAGS = $(INCFLAGS) @CFLAGS@ @DEFS@ @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @DEFS@ @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_CC  = $(CC)
AM_CXX = $(CXX)

.SUFFIXES=.c .h .cc .s. o

.c.s:
	$(CC) $(AM_CFLAGS) -S $< 

lib_LIBRARIES = libdslash.a


libdslash_a_SOURCES = 

#
# An AM_CONDITIONAL if statement to decide which fermion type
#

if SCALAR
libdslash_a_SOURCES += shift_table_scalar.cc cpp_dslash_scalar_32bit.cc \
	cpp_dslash_3d_scalar_32bit.cc \
	cpp_dslash_scalar_64bit.cc \
	cpp_dslash_3d_scalar_64bit.cc \
	cpp_clover_scalar_32bit.cc \
	cpp_clover_scalar_64bit.cc

if BUILD_NOTHREADS
libdslash_a_SOURCES += dispatch_scalar.cc
endif 

if BUILD_QMT
libdslash_a_SOURCES += dispatch_scalar_qmt.cc
endif
if BUILD_OMP
libdslash_a_SOURCES += dispatch_scalar_openmp.cc
endif

if BUILD_QDP_PACKERS
libdslash_a_SOURCES += qdp_packer_nopad.cc qdp_packer_3d.cc
endif
endif


if PARSCALAR
libdslash_a_SOURCES += cpp_dslash_parscalar_32bit.cc  cpp_dslash_parscalar_64bit.cc  cpp_dslash_parscalar_3d_32bit.cc cpp_dslash_parscalar_3d_64bit.cc cpp_dslash_parscalar_utils_32bit.cc cpp_dslash_parscalar_utils_64bit.cc

# Removed: cpp_clover_parscalar_32bit.cc cpp_clover_parscalar_64bit.cc

if BUILD_NOTHREADS
libdslash_a_SOURCES += dispatch_parscalar.cc
endif 

if BUILD_QMT
libdslash_a_SOURCES += dispatch_parscalar_qmt.cc
endif
if BUILD_OMP
libdslash_a_SOURCES += dispatch_parscalar_openmp.cc
endif

if BUILD_QDP_PACKERS
libdslash_a_SOURCES += qdp_packer_nopad.cc \
	qdp_packer_3d.cc
endif
endif
