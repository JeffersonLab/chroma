AC_INIT(sse_wilson_dslash, 1.0.0, bjoo@jlab.org)
AC_CONFIG_AUX_DIR(config)


AC_CANONICAL_HOST

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

dnl initialise automake
AM_INIT_AUTOMAKE

dnl put defines in include sse/config.h
AC_CONFIG_HEADERS([include/dslash_config_internal.h])
dnl dnl **************************************************************
dnl dnl ***                 Get User Options                       ***
dnl dnl **************************************************************

dnl learn the parallel architecture
AC_ARG_ENABLE(parallel-arch,
  AC_HELP_STRING([--enable-parallel-arch=<arch>],
    [Build Assembler for parallel architecture <arch>]),
  [ PARALLEL_ARCH="${enableval}"]
)

AC_ARG_ENABLE(testing,
   AC_HELP_STRING([--enable-testing],
   [Enable make check to automatically run tests, rather than just build them]),   [ do_tests="${enableval}" ]
)

AC_ARG_ENABLE(site-tests,
  AC_HELP_STRING([--enable-site-tests],
   [Enable building of site-ops tests. These still use inline assembler so may not be suitable for all compilers ]),
   [ site_tests="${enableval}"]
)


dnl check for QMP -- needed for parallel configures
AC_ARG_WITH(qmp,
  AC_HELP_STRING(
    [--with-qmp=DIR],
    [Build Assembler on top of QMP, where QMP is installed in DIR.]
  ),
  [QMP_HOME="${with_qmp}"]
)


dnl Use Jie Chens threads library
AC_ARG_WITH(qmt,
   AC_HELP_STRING(
    [--with-qmt=<DIR>],
    [Enable use of QMT Threads library installed in DIR]
   ),
   [ qmt_enabled="yes"
     QMT_HOME=${with_qmt} ],
   [ qmt_enabled="no" ]
)

dnl Use OpenMP threading
AC_ARG_ENABLE(openmp,
   AC_HELP_STRING(
    [--enable-openmp],
    [Enable building of OpenMP dispatcher. Pass flags necessary for OpenMP on CFLAGS]
   ),
   [ omp_enabled="yes" ],
   [ omp_enabled="no" ]
)

AC_ARG_ENABLE(nocompute,
	AC_HELP_STRING(
	  [--enable-nocompute],
	  [Disable the compute parts of dslash. This allows the timing of the routine with only communications.]
        ),
	[ nocompute_enabled="yes" ],
	[ nocompute_enabled="no" ]
)

AC_ARG_ENABLE(nocomms,
	AC_HELP_STRING(
	  [--enable-nocomms],
	  [Disable the comms parts of dslash. This allows one to gauge the effect of communications.]
        ),
	[ nocomms_enabled="yes" ],
	[ nocomms_enabled="no" ]
)

AC_ARG_ENABLE(sse2,
   AC_HELP_STRING(
     [--enable-sse],
     [Enables the SSE optimized code. This needs you to have a compiler which can cope with intel style xmmintrin intrinsics, and this option may not be suitable for non SSE Compatible platforms]
   ),
   [ sse_enabled=${enableval} ],
   [ sse_enabled="no"]
)

AC_ARG_ENABLE(sse3,
   AC_HELP_STRING(
     [--enable-sse],
     [Enables the SSE optimized code. This needs you to have a compiler which can cope with intel style xmmintrin intrinsics, and this option may not be suitable for non SSE Compatible platforms]
   ),
   [ sse3_enabled=${enableval} ],
   [ sse3_enabled="no"]
)

dnl with QDP -- QDP will override everything else
dnl since a parallel QDP will know its architecture
dnl precition and have all the flags for QMP set already
dnl in its CXXFLAGS etc etc
dnl 
dnl This flag and all the other flags are complimentary
dnl to each other therefore. Also --with-qdp is also 
dnl passed down the command line when compiling Chroma,
dnl so a sub-configure will also work
dnl we can also put the QDP++ packers in here too...
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
     [--with-qdp=DIR],
     [ Build Assembler with QDP++, where QDP++ is installed in DIR. ]
  ),
  [ QDPXX_HOME="${with_qdp}"
    QDPXX_GIVEN="yes"
  ],
  [ QDPXX_GIVEN="no" ]
)



dnl dnl ******************** DONE WITH USER INPUT *********************
dnl dnl ***        Now Configuration                               ****
dnl dnl ***************************************************************
if test "X${site_tests}X" = "XyesX" ; 
then 
 AC_MSG_NOTICE([ Building Site Tests ] )
else 
 AC_MSG_NOTICE([ Not building site tests ] )
fi
AM_CONDITIONAL(BUILD_SITE_TESTS, [ test "X${site_tests}X" == "XyesX" ] )

if test "X${do_tests}X" == "XyesX" ;
then
 AC_MSG_NOTICE([ Enabling automatic running of tests for make check ])
else
 AC_MSG_NOTICE([ make check will make but not run the tests ])
fi

AM_CONDITIONAL(RUN_TESTS, [ test "X${do_tests}X" == "XyesX" ] )
dnl 
dnl     ***************************************************************
dnl     ***       There are two distinct branches depending         ***
dnl     ***      on whether we are compiling with QDP++ or not      ***
dnl     ***************************************************************

dnl 
AM_CONDITIONAL(BUILD_QDP_TESTS, [ test "X${QDPXX_GIVEN}X" == "XyesX" ] )





if test "X${QDPXX_GIVEN}X" = "XyesX"; then 
     dnl Find the QDP++ config program 
     if test "X${QDPXX_HOME}X" = "XX" ; 
     then
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
     else
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
     fi

     if test "X${QDPXX_CONFIG}X" = "XX" ; then
	AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
     fi

     AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
     AC_MSG_CHECKING([if QDP++ uses checkerboarded layout])
     layout=`${QDPXX_CONFIG} --layout`
     if test "X${layout}X" == "Xcb2X" -o "X${layout}X" == "Xcb3dX" ; then 
	AC_MSG_RESULT([yes])
     else
        AC_MSG_RESULT([no])
	AC_MSG_ERROR([This Dslash Operator only works for checkerboarded QDP++. Sorry])
     fi


     AC_SUBST(SSEXX_AUX_CFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
     AC_MSG_NOTICE([QDP++ compile flags: ${SSEXX_AUX_CFLAGS}])

     AC_SUBST(SSEXX_AUX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
     AC_MSG_NOTICE([QDP++ linking flags: ${SSEXX_AUX_LDFLAGS}])

     AC_SUBST(SSEXX_AUX_LIBS,     "`${QDPXX_CONFIG} --libs`")
     AC_MSG_NOTICE([QDP++ libraries flags: ${SSEXX_AUX_LIBS}])

     QDPXX_ND="`${QDPXX_CONFIG} --Nd`"
     AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
     
     QDPXX_NC="`${QDPXX_CONFIG} --Nc`"
     AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])

     QDPXX_NS="`${QDPXX_CONFIG} --Ns`"
     AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])

     PARALLEL_ARCH="`${QDPXX_CONFIG} --parallel-arch`"

     AC_MSG_NOTICE([QDP++ has parallel arch: ${PARALLEL_ARCH} ] )

     dnl This function is defined in acinclude.m4
     dnl it tries to try to compile the program 
     PAC_QDPXX_LINK_CXX_FUNC(
	${SSEXX_AUX_CFLAGS},
	${SSEXX_AUX_LDFLAGS},
	${SSEXX_AUX_LIBS},
	,
	,
	[qdpxx_link_ok=yes],
	[qdpxx_link_ok=no]
     )

     AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
     if test "X${qdpxx_link_ok}X" = "XyesX" ; then 
	AC_MSG_RESULT(yes)
     else 
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([Cannot compile/link a program with QDP++.
    Use --with-qdp++=<dir> to select a working version.])
    fi 

    AC_DEFINE(SSE_USE_QDPXX, [1], [Dslash uses QDPXX])

    dnl *********************************************************
    dnl ***              end of section with QDP++            ***
    dnl *********************************************************
else 
     dnl ***************************************************************
     dnl *** Section without QDP++ -- Need to check for QMP possibly ***
     dnl ***************************************************************
     case "${PARALLEL_ARCH}" in
	  scalar)
		AC_MSG_NOTICE([ Scalar Build! Not checking for QMP ])
		AC_SUBST(SSEXX_AUX_CFLAGS)
		AC_SUBST(SSEXX_AUX_LDFLAGS)
		AC_SUBST(SSEXX_AUX_LIBS)
		;;

	  parscalar)
		AC_MSG_NOTICE([ Parscalar build! Checking for QMP ])

		dnl look for qmp-config program
		if test "X${QMP_HOME}X" = "XX" ; then
		   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [])
		else
		   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMP_HOME}/bin:${PATH}])
		fi

		dnl look for flags needed for QMP
		if test "X${QMP_CONFIG}X" != "XX" ; then
		   AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
		   AC_SUBST(SSEXX_AUX_CFLAGS,  "`${QMP_CONFIG} --cflags`")
		   AC_MSG_NOTICE([QMP compile flags: ${SSEXX_AUX_CFLAGS}])

		   AC_SUBST(SSEXX_AUX_LDFLAGS, "`${QMP_CONFIG} --ldflags`")
		   AC_MSG_NOTICE([QMP linking flags: ${SSEXX_AUX_LDFLAGS}])

		   AC_SUBST(SSEXX_AUX_LIBS,    "`${QMP_CONFIG} --libs`")
		   AC_MSG_NOTICE([QMP libraries flags: ${SSEXX_AUX_LIBS}])
		else 
		   AC_MSG_WARN([QMP configuration program qmp-config not found.])
		   AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS before configure])
		fi

		dnl verify we can compile and link against QMP, if needed
		AC_PROG_CXX
		PAC_QMP_LINK_CC_FUNC(
			[${SSEXX_AUX_CFLAGS}],
			[${SSEXX_AUX_LDFLAGS}],
			[${SSEXX_AUX_LIBS}],
			,
			,
			[qmp_link_ok=yes],
			[qmp_link_ok=no]
		)
		
		AC_MSG_CHECKING([if we can compile/link of a simple QMP C program])
		if test "X${qmp_link_ok}X" = "XyesX" ; then
		   AC_MSG_RESULT(yes)
		else
		   AC_MSG_RESULT(no)
		   AC_MSG_ERROR([Cannot compile/link a basic QMP C program! Check QMP_CFLAGS, QMP_LDFLAGS, QMP_LIBS.])
		fi
		;;
	*)
		AC_MSG_ERROR([ Unknown value for --enable-parallel-arch ])
		;;
	esac
fi

case ${PARALLEL_ARCH} in
  scalar|parscalar)
	;;
  *)
	AC_MSG_ERROR([Unknown Parallel Architecture: ${PARALLEL_ARCH} ] ) 
	;; 
esac


dnl ************************************************************************
dnl **** Compile Defines to select right stuff. Always SSE now             *
dnl ************************************************************************

case ${PARALLEL_ARCH} in
scalar)
	AC_DEFINE_UNQUOTED(CPP_DSLASH_SCALAR, [1], [Scalar Arch])
	;;
parscalar)
	AC_DEFINE_UNQUOTED(CPP_DSLASH_PARSCALAR, [1], [Parscalar Arch])
	;;
*)
	AC_MSG_ERROR([Unsupported parallel arch])
	;;
esac
	

AM_CONDITIONAL(SCALAR, [ test "X${PARALLEL_ARCH}X" == "XscalarX" ] )
AM_CONDITIONAL(PARSCALAR, [ test "X${PARALLEL_ARCH}X" == "XparscalarX" ] )
AM_CONDITIONAL(BUILD_QDP_PACKERS, [ test "X${QDPXX_GIVEN}X" == "XyesX" ])				 

dnl NOCOMPUTE DIagnostic
if test "X${nocompute_enabled}X" == "XyesX";
then 
     AC_MSG_NOTICE([Disabling The Compute Parts of the Dslash])
     AC_DEFINE([SSEDSLASH_4D_NOCOMPUTE],[1],[Disable computation in 4D kernels])
fi

if test "X${nocomms_enabled}X" == "XyesX";
then 
     AC_MSG_NOTICE([Disabling the Comms in the Dslash])
     AC_DEFINE([SSEDSLASH_4D_NOCOMMS],[1],[Disable computation in 4D kernels])
fi


dnl QMT Threads stuff
if test "X${qmt_enabled}X" == "XyesX"; 
then
	AC_MSG_NOTICE([Configuring QMT Threading])
	if test "X${omp_enabled}X" == "XyesX";
	then 
	  AC_MSG_ERROR([Cannot have OpenMP and QMT threading defined simultaneously])
	fi

	AC_DEFINE([DSLASH_USE_QMT_THREADS], [1], [ Use QMT Threads library ])
	AC_SUBST(QMT_CFLAGS, "-I${QMT_HOME}/include")
	AC_SUBST(QMT_LDFLAGS, "-L${QMT_HOME}/lib")
	AC_SUBST(QMT_LIBS, "-lqmt -lpthread -lm")
else 
	AC_SUBST(QMT_CFLAGS, "")
	AC_SUBST(QMT_LDFLAGS, "")
	AC_SUBST(QMT_LIBS, "")

fi

if test "X${omp_enabled}X" == "XyesX";
then
	AC_MSG_NOTICE([Configuring OpenMP Threading])
	if test "X${qmt_enabled}X" == "XyesX";
        then 
          AC_MSG_ERROR([Cannot have OpenMP and QMT threading defined simultaneously])
	fi
	  
	AC_DEFINE([DSLASH_USE_OMP_THREADS], [1], [ Use OpenMP Threads ])
fi

if test "X${sse_enabled}X" == "XyesX";
then 
     AC_MSG_NOTICE([Configuring SSE2 Code]);
     AC_DEFINE([DSLASH_USE_SSE2],[1], [ Use SSE2 Code in kernels])
else
     AC_MSG_NOTICE([Not Configuring SSE2 Code]);
fi

if test "X${sse3_enabled}X" == "XyesX";
then 
     AC_MSG_NOTICE([Configuring SSE2 and SSE3 Code]);
     AC_DEFINE([DSLASH_USE_SSE2],[1], [ Use SSE2 Code in kernels])
     AC_DEFINE([DSLASH_USE_SSE2],[1], [ Use SSE3 Code in kernels])
else
     AC_MSG_NOTICE([Not Configuring SSE3 Code]);
fi

AM_CONDITIONAL(BUILD_QMT, [test "x${qmt_enabled}x" = "xyesx" ])
AM_CONDITIONAL(BUILD_OMP, [test "x${omp_enabled}x" = "xyesx" ])
AM_CONDITIONAL(BUILD_NOTHREADS, [ test "x${omp_enabled}x" = "xnox"  -a  "x${qmt_enabled}x" = "xnox" ])
dnl dnl produce output
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_OUTPUT
