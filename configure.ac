AC_INIT(chroma,1.22,edwards@jlab.org)
AC_PREREQ([2.57])
AC_CONFIG_AUX_DIR(./config)
AC_CANONICAL_HOST

AM_INIT_AUTOMAKE

dnl stupid templates to shut up stupid autoheader
dnl AC_DEFINE([BUILD_PAB_WILSON_DSLASH],[],[Build Pab Wilson Dslash])
dnl AC_DEFINE([BUILD_SSE_WILSON_DSLASH],[],[Build SSE wilson Dslash])
dnl AC_DEFINE([CHROMA_BUILD_STAGGERED],[], [Build Staggered version of Chroma])
dnl AC_DEFINE([CHROMA_BUILD_WILSON],[],[Build Wilson version of Chroma])

AC_CONFIG_HEADERS([lib/chroma_config_internal.h])

# AC_PRO_MAKE_SET        # (RGE) Causes error: do not find in autoconf docs
AC_PROG_RANLIB
##################################
# Check for the correct AR
##################################

AC_CHECK_TOOL(AR, ar, [ar])

dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl George T. Fleming, 12 February 2003
dnl Complete rewrite of configure switches based on new config scripts.

dnl --with-qdp=DIR
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
    [--with-qdp=DIR],
    [Build Chroma on top of QDP++ where QDP++ is installed in DIR]
  ),
  [QDPXX_HOME="$with_qdp"]
)

AC_ARG_ENABLE(sse_wilson_dslash,
  AC_HELP_STRING( 
    [--enable-sse-wilson-dslash],
    [Build and Use the SSE2 Wilson Dslash Library])
)

AC_ARG_ENABLE(sse_dwf_cg,
  AC_HELP_STRING( 
    [--enable-sse-dwf-cg],
    [Build and Use the SSE DWF CG Library])
)

AC_ARG_ENABLE(pab_wilson_dslash,
  AC_HELP_STRING(
    [--enable-pab-wilson-dslash],
    [Use PABs Wilson Dslash Library])
)

AC_ARG_ENABLE(asqtad_cps_inverter,
   AC_HELP_STRING(
    [--enable-asqtad-cps-inverter],
    [Wrap the CPS asqtad inverter])
)

AC_ARG_ENABLE(opt_cfz_linop,
  AC_HELP_STRING( 
    [--enable-opt-cfz-linop],
    [Build and use the generic optimized Continued Frac. Zolo. (CFZ) Linear Oper.])
)

AC_ARG_ENABLE(gmp,
   AC_HELP_STRING(
    [--enable-gmp],
    [Enable use of GMP - GNU Multi-Precision library])
)

if test "X${QDPXX_HOME}X" = "XX" ; then
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
else
  AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
fi

if test "X${QDPXX_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
fi

AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
AC_SUBST(CXX, "`${QDPXX_CONFIG} --cxx`")
AC_MSG_NOTICE([QDP++ compiler: ${CXX}])
AC_SUBST(QDPXX_CXXFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
AC_MSG_NOTICE([QDP++ compile flags: ${QDPXX_CXXFLAGS}])
AC_SUBST(QDPXX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
AC_MSG_NOTICE([QDP++ linking flags: ${QDPXX_LDFLAGS}])
AC_SUBST(QDPXX_LIBS,     "`${QDPXX_CONFIG} --libs`")
AC_MSG_NOTICE([QDP++ libraries flags: ${QDPXX_LIBS}])
AC_SUBST(QDPXX_ND,       "`${QDPXX_CONFIG} --Nd`")
if test "X${QDPXX_ND}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine QDP++ spacetime dimension])
else
  AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
fi
AC_SUBST(QDPXX_NC,       "`${QDPXX_CONFIG} --Nc`")
AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])
AC_SUBST(QDPXX_NS,       "`${QDPXX_CONFIG} --Ns`")
if test "X${QDPXX_NS}X" = "XX" ; then
  AC_MSG_ERROR([Cannot determine number of spin components in QDP++])
else
  AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])
fi

AC_SUBST(QDPXX_PARALLEL_ARCH, "`${QDPXX_CONFIG} --parallel-arch`")
AC_MSG_NOTICE([QDP++ parallel arch: ${QDPXX_PARALLEL_ARCH}])

AC_SUBST(QDPXX_PRECISION, "`${QDPXX_CONFIG} --precision`")
AC_MSG_NOTICE([QDP++ precision: ${QDPXX_PRECISION}])

dnl SSE Wilson Dslash stuff
case "$enable_sse_wilson_dslash" in
  yes) 
        AC_MSG_NOTICE( [Configuring SSE Wilson Dslash ] )
	AC_CONFIG_SUBDIRS(other_libs/intel_sse_wilson_dslash)
	AC_SUBST(SSE_DSLASH_DIR, [intel_sse_wilson_dslash])
	AC_DEFINE([BUILD_SSE_WILSON_DSLASH],[],[ Build SSE Dslash Op ])
	;;
  *)
	AC_SUBST(SSE_DSLASH_DIR)
	AC_MSG_NOTICE( [Not building SSE Dslash ] )
        ;;
esac 
AC_SUBST(SSE_DSLASH_ENABLED,[$enable_sse_wilson_dslash]	)	

AM_CONDITIONAL(BUILD_SSE_WILSON_DSLASH,
  [test "x${enable_sse_wilson_dslash}x" = "xyesx" ])


dnl SSE DWF CG stuff
case "$enable_sse_dwf_cg" in
  yes) 
        AC_MSG_NOTICE([Configuring SSE DWF CG])
	AC_CONFIG_SUBDIRS(other_libs/intel_sse_dwf_cg)
	AC_SUBST(SSE_DWF_CG_DIR, [intel_sse_dwf_cg])
	AC_DEFINE([BUILD_SSE_DWF_CG],[],[Build SSE DWF CG])
	;;
  *)
	AC_SUBST(SSE_DWF_CG_DIR)
	AC_MSG_NOTICE([Not building SSE DWF CG])
        ;;
esac 
AC_SUBST(SSE_DWF_CG_ENABLED,[$enable_sse_dwf_cg])

AM_CONDITIONAL(BUILD_SSE_DWF_CG,
  [test "x${enable_sse_dwf_cg}x" = "xyesx" ])


dnl QCDOC Wilson Dslash stuff
case "$enable_pab_wilson_dslash" in
  noarch)
	ac_build_pab_wilson_dslash="yes" 
        AC_MSG_NOTICE( [Configuring PAB Wilson Dslash ] )
	AC_DEFINE([BUILD_PAB_WILSON_DSLASH],[],[ Build Peters Dslash Op])
	;;
  *)
	ac_build_pab_wilson_dslash="no"
	AC_MSG_NOTICE( [Not building PAB Dslash ] )
        ;;
esac 

AM_CONDITIONAL(BUILD_PAB_WILSON_DSLASH,
  [test "x${ac_build_pab_wilson_dslash}x" = "xyesx" ])

dnl QCDOC Asqtad Inverter stuff
case "$enable_asqtad_cps_inverter" in
 yes)
        AC_MSG_NOTICE([Configuring Asqtad CPS Inverter Wrapper])
	AC_DEFINE([BUILD_CPS_ASQTAD_INVERTER],[],[Build CPS ASQTAD INVERTER])
	;;
 *) 
        AC_MSG_NOTICE([Not Configuring Asqtad CPS Inverter Wrapper])
        ;;	
esac

AM_CONDITIONAL(BUILD_ASQTAD_CPS_INVERTER_WRAPPER,
  [test "x${enable_asqtad_cps_inverter}x" = "xyesx" ])


dnl Generically optimized CFZ Linop
case "$enable_opt_cfz_linop" in
  yes)
	ac_build_opt_cfz_linop="yes" 
        AC_MSG_NOTICE( [Configuring optimized Continued Frac. (CFZ) LinOp] )
	AC_DEFINE([BUILD_CFZ_CFZ_LINOP],[],[ Build optimized Continued Frac. (CFZ) LinOp])
	;;
  *)
	ac_build_opt_cfz_linop="no"
	AC_MSG_NOTICE( [Not building optimized CFZ LinOp] )
        ;;
esac 

AM_CONDITIONAL(BUILD_OPT_CFZ_LINOP,
  [test "x${enable_opt_cfz_linop}x" = "xyesx" ])


dnl Enable GMP
case "$enable_gmp" in
  yes) 
        AC_MSG_NOTICE( [Enabling GNU Multi-Precision Library] )
	AC_DEFINE([BUILD_GMP_REMEZ],[],[ Build GMP Remez code ])
	;;
  *)
	AC_MSG_NOTICE( [Not building GMP Remez code ] )
        ;;
esac 

AM_CONDITIONAL(BUILD_GMP_REMEZ,
  [test "x${enable_gmp}x" = "xyesx" ])


dnl Try to compile a QDP++ program to check the --with options
AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])

dnl This function is defined in acinclude.m4
dnl it tries to try to compile the program 
PAC_QDPXX_LINK_CXX_FUNC(
  ${QDPXX_CXXFLAGS},
  ${QDPXX_LDFLAGS},
  ${QDPXX_LIBS},
  ,
  ,
  [qdpxx_link_ok=yes],
  [qdpxx_link_ok=no]
)
if test "X${qdpxx_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with QDP++.
    Use --with-qdp++=<dir> to select a working version.])
fi 


##################################
# Check for programs to build docs
##################################
AC_CHECK_PROG(LATEX, latex, latex, true)
AC_CHECK_PROG(DVIPS, dvips, dvips, true)
AC_CHECK_PROG(DVIPDF, dvipdf, dvipdf, true)
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen, true)

#######################
# Produce output
#######################

#
# Set up the doxygen file
# Substitute the right directories
#
# The docs stuff needs more work. Turn it all off for now.
# AC_CONFIG_FILES(docs/chromadoc)

dnl Make the Makefiles
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(chroma-config)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(mainprogs/Makefile)
AC_CONFIG_FILES(mainprogs/tests/Makefile)
AC_CONFIG_FILES(mainprogs/main/Makefile)
# AC_CONFIG_FILES(docs/Makefile)
AC_CONFIG_FILES(other_libs/Makefile)
AC_OUTPUT
