AC_INIT(chroma,0.2,edwards@jlab.org)
AC_PREREQ([2.57])
AC_CONFIG_AUX_DIR(./config)
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

AM_INIT_AUTOMAKE
# AC_PRO_MAKE_SET        # (RGE) Causes error: do not find in autoconf docs
AC_PROG_RANLIB

dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl George T. Fleming, 12 February 2003
dnl Complete rewrite of configure switches based on new config scripts.


dnl --enable-fermion-type
AC_ARG_ENABLE(
  [fermion-type],
  AC_HELP_STRING(
    [--enable-fermion-type=<fermion type>],
    [Set the fermion type (wilson or staggered) default is wilson]
  ),
  [
    case "$enable_fermion_type" in
      STAG*) ac_fermion_type="staggered" ;;
      Stag*) ac_fermion_type="staggered" ;;
      stag*) ac_fermion_type="staggered" ;;

      WIL*)  ac_fermion_type="wilson" ;;
      Wil*)  ac_fermion_type="wilson" ;;
      wil*)  ac_fermion_type="wilson" ;;

      *)     ac_fermion_type="$enable_fermion_type" ;;
    esac
  ],
  [ac_fermion_type="wilson"]
) 


dnl --with-qdp=DIR
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
    [--with-qdp=DIR],
    [Build QDP++ Apps on top of QDP++ where QDP++ is installed in DIR]
  ),
  [QDPXX_HOME="$with_qdp"]
)

AC_ARG_ENABLE(sse_wilson_dslash,
  AC_HELP_STRING( 
    [--enable-sse-wilson-dslash],
    [Build and Use the SSE2 Wilson Dslash Library ]),
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

dnl
dnl
case $ac_fermion_type in
  staggered)
    if test $QDPXX_NS -ne 1 ; then
      AC_MSG_ERROR([STAGGERED-like fermions should have Ns=1. Check QDP++])
    else 
      AC_MSG_NOTICE([Configuring QDP++ Apps with STAGGERED-like fermions])
    fi
    ;;
  wilson)
    AC_MSG_NOTICE([Configuring QDP++ Apps with WILSON-like fermions])
    ;;
  *)
    AC_MSG_ERROR([${ac_fermion_type}-like fermions not supported yet.])
    ;;
esac

case "$enable_sse_wilson_dslash" in
  yes) 
        if test "x$ac_fermion_type" != "xwilson";
	then	
	    AC_MSG_ERROR([ Attempting to build Assembler Wilson Dslash with NON Wilson fermion type ])
	fi
        AC_MSG_NOTICE( [Configuring SSE Wilson Dslash ] )
	AC_CONFIG_SUBDIRS(other_libs/intel_sse_wilson_dslash)
	AC_SUBST(SSE_DSLASH_DIR, [intel_sse_wilson_dslash])
	;;
  no)
	AC_SUBST(SSE_DSLASH_DIR)
	AC_MSG_NOTICE( [Not building SSE Dslash ] )
        ;;
  *)
	AC_MSG_ERROR( [Unknown value for --enable-sse-wilson-dslash ] )
  ;;
esac 


AM_CONDITIONAL(BUILD_SSE_WILSON_DSLASH,
  [test "x${enable_sse_wilson_dslash}x" = "xyesx" ])

AM_CONDITIONAL(BUILD_STAGGERED_LIBS,
  [test "X${ac_fermion_type}X" = "XstaggeredX"])
AM_CONDITIONAL(BUILD_WILSON_LIBS,
  [test "X${ac_fermion_type}X" = "XwilsonX"])

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
AC_CONFIG_FILES(docs/chromadoc)

dnl Make the Makefiles
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(mainprogs/Makefile)
AC_CONFIG_FILES(mainprogs/tests/Makefile)
AC_CONFIG_FILES(mainprogs/main/Makefile)
AC_CONFIG_FILES(docs/Makefile)
AC_CONFIG_FILES(other_libs/Makefile)
AC_OUTPUT
