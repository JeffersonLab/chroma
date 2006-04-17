#!/usr/bin/bash

# Get build functions
. ./make_functions.sh

QDP_VERSION="1.18.1"
CHROMA_VERSION="2.5.5"
HERE=`/usr/bin/pwd`
#SRCDIR=/qcdoc/sfw/chroma/v2/src
SRCDIR=${HERE}
INSTALL_ROOT=${HERE}/../qcdoc/v2.6.0-CJ-recompiled/aix5.2f
QOS=/qcdoc/sfw/qos/devel/v2.6.0-CJ-recompiled/aix5.2f
PRECISION=double
HOST_SYS=powerpc-gnu-elf
BUILD_SYS=none
BAGEL_VERSION=1.3.3
BAGEL_QDP_VERSION=1.1.4
BAGEL_DSLASH_VERSION=1.3.3
LIBXML_VERSION=2.6.6
BAGEL_QDP=${HERE}/bagel_qdp
BAGEL_INSTALL_DIR=${INSTALL_ROOT}/bagel
BAGEL_CPU=ppc440
BAGEL_ALLOC=qalloc
BAGEL_COMM=qmp
BAGEL_WILSON_DIR=${INSTALL_ROOT}/bagel_wilson_dslash_${PRECISION}_${BAGEL_CPU}_${BAGEL_ALLOC}_${BAGEL_COMM}
BAGELQDP_INSTALLDIR=${INSTALL_ROOT}/bagel_qdp_${PRECISION}_${BAGEL_CPU}

LIBXML_SRCDIR=${SRCDIR}/libxml2-${LIBXML_VERSION}
LIBXML=${INSTALL_ROOT}/libxml2-${LIBXML_VERSION}
QDP_PARALLEL_ARCH=parscalar
QDP_DO_BLAS=yes
CHROMA_DO_PAB_DSLASH=yes
CHROMA_DO_GMP=yes
CHROMA_GMPDIR=/qcdoc/sfw/packages/gmp/qos-2.5.9a
# Munge directory names
QDP_INSTALLDIR=${INSTALL_ROOT}/qdp-${PRECISION}-${QDP_VERSION}

if test "X${QDP_DO_BAGEL}X" == "XyesX";
then 
	QDP_INSTALLDIR=${QDP_INSTALLDIR}-bagel_qdp;
fi

CHROMA_INSTALLDIR=${INSTALL_ROOT}/chroma-${PRECISION}-${CHROMA_VERSION}
if test "X${CHROMA_DO_PAB_DSLASH}X" == "XyesX";
then
        CHROMA_INSTALLDIR=${CHROMA_INSTALLDIR}-bagel_dslash;
fi

##
## The actual building is done here
##
##

## Build BAGEL
build_bagel ${SRCDIR}/bagel-${BAGEL_VERSION} ${BAGEL_INSTALL_DIR}

source ${QOS}/scripts/setup.sh
# Build Wilson Dslash

build_bagel_wilson_dslash ${SRCDIR}/bagel_wilson_dslash-${BAGEL_DSLASH_VERSION} ${BAGEL_WILSON_DIR} ${BAGEL_INSTALL_DIR} ${PRECISION} ${BAGEL_COMM} ${BAGEL_ALLOC}  ${BAGEL_CPU} ${HOST_SYS} ${BUILD_SYS} ${QOS}

build_bagel_qdp ${SRCDIR}/bagel_qdp-${BAGEL_QDP_VERSION} ${BAGEL_INSTALL_DIR} ${BAGELQDP_INSTALLDIR} ${PRECISION} ${BAGEL_CPU} ${HOST_SYS} ${BUILD_SYS}

build_libxml ${LIBXML_SRCDIR} ${LIBXML} ${HOST_SYS} ${BUILD_SYS}

## Build QDP++
build_qdp  ${SRCDIR}/qdp-${QDP_VERSION} ${QDP_INSTALLDIR} ${QOS} ${LIBXML} ${PRECISION} ${QDP_DO_BLAS} ${HOST_SYS} ${BUILD_SYS} ${BAGELQDP_INSTALLDIR}

## Build Chroma
build_chroma ${SRCDIR}/chroma-${CHROMA_VERSION} ${CHROMA_INSTALLDIR} ${QDP_INSTALLDIR} ${HOST_SYS} ${BUILD_SYS} ${CHROMA_DO_PAB_DSLASH} ${BAGEL_WILSON_DIR} ${CHROMA_DO_GMP} ${CHROMA_GMPDIR}

# find . -name "*" -type d -exec chmod ugo+rx {} \; -print
