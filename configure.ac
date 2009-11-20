AC_INIT(chroma,3.36.1,edwards@jlab.org)
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
AC_PROG_CXX(g++)
AC_PROG_CC(gcc)
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

dnl CPP_WILSON_DSLASH_OPTIONS
AC_ARG_ENABLE(cpp_wilson_dslash,
   AC_HELP_STRING(
    [--enable-cpp-wilson-dslash],
    [Build and Use the CPP Wilson Dslash Library. Can also specify --enable-sse2 or --enable-sse3 and requires a QMP location in parscalar more with --with-qmp]   )
)

AC_ARG_ENABLE(sse2,
   AC_HELP_STRING(
    [--enable-sse2],
    [Enable SSE2 support in the CPP Wilson Dslash package]
   )
)

AC_ARG_ENABLE(sse3,
   AC_HELP_STRING(
    [--enable-sse3],
    [Enable SSE3 support in the CPP Wilson Dslash package]
   )
)

AC_ARG_ENABLE(cg-dwf,
	AC_HELP_STRING([--enable-cg-dwf=<cpuarch> <cpuarch>=sse,bluelight,altivec]),
	[cg_dwf_cpu_arch=${enableval}
	 build_cg_dwf="yes"],
	[cg_dwf_cpu_arch="none"
	 build_cg_dwf="no"]
)

AC_ARG_ENABLE(cg-dwf-lowmem,
	AC_HELP_STRING([--enable-cg-dwf-lowmem=<yes/no> Enable Low/High Memory mode of CG-DWF inverter]), 
	[ cg_dwf_lowmem="yes"],
	[ cg_dwf_lowmem="no" ]
)

AC_ARG_ENABLE(cg-solver-restart,
	AC_HELP_STRING([ Enable a single restart in the CG linop syssolvers ]),
	[cg_do_one_restart="yes"],
	[cg_do_one_restart="no"]
)

AC_ARG_WITH(mdwf,
	AC_HELP_STRING([--with-mdwf=<install location>]),
	[mdwf_path=${withval}
	 build_mdwf="yes"],
	[mdwf_path=""
	 build_mdwf="no"]
)
AC_ARG_ENABLE(asqtad_level3_inverter,
   AC_HELP_STRING(
    [--enable-asqtad-level3-inverter],
    [Wrap the level3 asqtad inverter])
)


AC_ARG_WITH(bagel_wilson_dslash,
  AC_HELP_STRING(
    [--with-bagel-wilson-dslash=DIR],
    [Use Peter Boyles BAGEL Wilson Dslash Library installed in DIR]
  ),
    [bagel_wilson_dslash_enabled="yes"
     BAGEL_WILSON_DSLASH_HOME="${with_bagel_wilson_dslash}"]
)

AC_ARG_ENABLE(bagel_wilson_dslash_sloppy,
	AC_HELP_STRING(
        [--enable-bagel-wilson-dslash-sloppy], 
	[Enable Sloppy Precision in BAGEL Dslash -- different internal and external precisions]),
	[bagel_wilson_dslash_sloppy="yes"],
	[bagel_wilson_dslash_sloppy="no"]
)


AC_ARG_WITH(bagel_clover,
  AC_HELP_STRING(
    [--with-bagel-clover=DIR],
    [Use the BAGEL Clover Apply Library installed in DIR]
  ),
    [bagel_clover_enabled="yes"
     BAGEL_CLOVER_HOME="${with_bagel_clover}"]
)

AC_ARG_ENABLE(ssed_clover,
	AC_HELP_STRING(
        [--enable-ssed-clover], 
	[Enable SSE Double Prec clover operator]),
	[ssed_clover_enabled="${enableval}"],
	[ssed_clover_enabled="no"]
)

AC_ARG_ENABLE(lapack,
  AC_HELP_STRING([--enable-lapack=<type>],
    [Use the type of lapack bindings specified by <type>. Currently, this can be lapack]),
  [LAPACK_BINDING="${enableval}"]
)

AC_ARG_ENABLE(opt-eigcg,
  AC_HELP_STRING([--enable-opt-eigcg],
    [Build optimized eigcg]),
  [opt_eigcg_enabled="yes"]
)

AC_ARG_ENABLE(opt_cfz_linop,
  AC_HELP_STRING( 
    [--enable-opt-cfz-linop],
    [Build and use the generic optimized Continued Frac. Zolo. (CFZ) Linear Oper.])
)

AC_ARG_ENABLE(enable_generic_scalarsite_bicgstab_kernels,
	AC_HELP_STRING(
	[--enable-generic-scalarsite-bicgstab-kernels],
	[Switch on kernels to reduce Reliable BiCGStab BLAS memory bandwidth on threaded machines])
)

AC_ARG_ENABLE(enable_sse_scalarsite_bicgstab_kernels,
	AC_HELP_STRING(
	[--enable-sse-scalarsite-bicgstab-kernels],
	[Switch on SSE kernels to reduce Reliable BiCGStab BLAS memory bandwidth on threaded machines])
)

AC_ARG_ENABLE(testcase-runner,
  AC_HELP_STRING([--enable-testcase-runner=script],
    [Use <script> to run testcases: trivial|cobalt|6n_mpirun_rsh|7n_mpirun_rsh]),
    [chroma_testcase_runner=${enable_testcase_runner}],
    [chroma_testcase_runner="trivial_runner.sh"]
)

AC_ARG_WITH(gmp,
   AC_HELP_STRING(
    [--with-gmp=<DIR>],
    [Enable use of GMP - GNU Multi-Precision library, installed in DIR]
   ),
    [gmp_enabled="yes"
     GMP_HOME=${with_gmp}],
    [gmp_enabled="no"]
)

AC_ARG_WITH(qmt,
   AC_HELP_STRING(
    [--with-qmt=<DIR>],
    [Enable use of QMT Threads library installed in DIR]
   ),
   [ qmt_enabled="yes"
     QMT_HOME=${with_qmt} ],
   [ qmt_enabled="no" ]
)

AC_ARG_WITH(quda,
   AC_HELP_STRING(
      [--with-quda=<DIR>],
      [Use QUDA Library installed in DIR]
   ),
   [ quda_enabled="yes"
     QUDA_HOME=${with_quda} ],
   [ quda_enabled="no" ]
)

AC_ARG_WITH(cuda,
   AC_HELP_STRING(
      [--with-cuda=<DIR>],
      [Use CUDA installation in DIR]
   ),
   [ CUDA_HOME=${with_cuda} ],
   [ CUDA_HOME="/usr/local/cuda" ]
)

dnl level 3 libraries for asqtad inverter
AC_ARG_WITH(asqtadlevel3,
  AC_HELP_STRING(
    [--with-asqtadlevel3=DIR],
    [Build Chroma on top of the level3 asqtad libraries, where the libraries installed in DIR.]
  ),
  [QDPCVER_HOME="$with_asqtadlevel3"]
)

if test "X${QDPCVER_HOME}X" = "XX" ; then
  AC_MSG_NOTICE([QDPCVER_HOME variable not set.])
else
  AC_SUBST(QDPC_DIR, "${QDPCVER_HOME}")
fi


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

dnl Try to compile a QDP++ program to check the --with options
dnl This function is defined in acinclude.m4
AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
PAC_QDPXX_LINK_CXX_FUNC(
  ${QDPXX_CXXFLAGS},
  ${QDPXX_LDFLAGS},
  ${QDPXX_LIBS},
  ,
  ,
  [qdpxx_link_ok=yes],
  [qdpxx_link_ok=no]
)
if test "X${qdpxx_link_ok}X" = "XyesX"; 
then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with QDP++. Use --with-qdp++=<dir> to select a working version.])
fi 


dnl ************************************************************************
dnl **** SSE Wilson dslash                                              ****
dnl ************************************************************************
case "$enable_sse_wilson_dslash" in
  yes) 
       if test "X${enable_cpp_wilson_dslash}X" = "XyesX";
       then 
          AC_MSG_ERROR([Cannot compile and link both CPP Wilson Dslash and SSE Wilson Dslash. Disable one of them])
       fi
        AC_MSG_NOTICE([Configuring SSE Wilson Dslash])
	AC_CONFIG_SUBDIRS(other_libs/sse_wilson_dslash)
        AC_MSG_NOTICE([Finished configuring SSE Wilson Dslash])
	AC_SUBST(SSE_DSLASH_DIR, [sse_wilson_dslash])
	AC_DEFINE([BUILD_SSE_WILSON_DSLASH],[],[ Build SSE Dslash Op ])
        AC_SUBST(SSE_DSLASH_ENABLED,[yes]	)	  
	;;
  *)
	AC_SUBST(SSE_DSLASH_DIR)
	AC_MSG_NOTICE( [Not building SSE Dslash ] )
        AC_SUBST(SSE_DSLASH_ENABLED,[no]	)
        ;;
esac 


AM_CONDITIONAL(BUILD_SSE_WILSON_DSLASH,
  [test "x${enable_sse_wilson_dslash}x" = "xyesx" ])


dnl ************************************************************************
dnl **** CPP DSlash Stuff                                              *****
dnl ************************************************************************
case "$enable_cpp_wilson_dslash" in
 yes)
       if test "X${enable_sse_wilson_dslash}X" = "XyesX";
       then 
          AC_MSG_ERROR([Cannot compile and link both CPP Wilson Dslash and SSE Wilson Dslash. Disable one of them])
       fi
        AC_MSG_NOTICE([Configuring CPP Wilson Dslash])
	AC_CONFIG_SUBDIRS(other_libs/cpp_wilson_dslash)
        AC_MSG_NOTICE([Finished configuring CPP Wilson Dslash])
	AC_SUBST(CPP_DSLASH_DIR, [cpp_wilson_dslash])
	AC_DEFINE([BUILD_CPP_WILSON_DSLASH],[],[ Build CPP Dslash Ops ])
        AC_SUBST(CPP_DSLASH_ENABLED,[yes]	)	  
	;;
  *)
	AC_SUBST(CPP_DSLASH_DIR)
	AC_MSG_NOTICE( [Not building CPP Dslash ] )
        AC_SUBST(CPP_DSLASH_ENABLED,[no]	)
        ;;
esac 
AM_CONDITIONAL(BUILD_CPP_WILSON_DSLASH,
  [test "x${enable_cpp_wilson_dslash}x" = "xyesx" ])

dnl ************************************************************************
dnl **** Generic Scalarsite BiCGStab Stuff
dnl ************************************************************************
case "$enable_generic_scalarsite_bicgstab_kernels" in
 yes)
        AC_MSG_NOTICE([Enabling Generic Scalarsite BiCGStab Kernels])
	AC_DEFINE([BUILD_SCALARSITE_BICGSTAB],[],[Build Scalarsite BICGStab Kernels])
	AC_DEFINE([BUILD_GENERIC_SCALARSITE_BICGSTAB],[],[Use Generic Kernels ])
	;;
  *)
        ;;
esac 

dnl ************************************************************************
dnl **** Generic Scalarsite BiCGStab Stuff
dnl ************************************************************************
case "$enable_sse_scalarsite_bicgstab_kernels" in
 yes)
        AC_MSG_NOTICE([Enabling SSE Scalarsite BiCGStab Kernels])
	AC_DEFINE([BUILD_SCALARSITE_BICGSTAB],[],[Build Scalarsite BICGStab Kernels])
	AC_DEFINE([BUILD_SSE_SCALARSITE_BICGSTAB],[],[Use Generic Kernels ])
	;;
  *)
        ;;
esac 
AM_CONDITIONAL(BUILD_SCALARSITE_BICGSTAB,
  [test "x${enable_sse_scalarsite_bicgstab_kernels}x" = "xyesx" -o "x${enable_generic_scalarsite_bicgstab_kernels}x"  ])




dnl ************************************************************************
dnl **** CG DWF stuff                                                   ****
dnl ************************************************************************
if [ test "x${build_cg_dwf}x" = "xyesx" ];
then
        AC_MSG_NOTICE([Configuring CG DWF])
	if [ test "x${build_mdwf}x" = "xyesx" ];
	then 
	  AC_MSG_ERROR([Cannot have both CG-DWF    and MDWF configured at the same time])
        fi

	case ${cg_dwf_cpu_arch} in
	sse)
		AC_MSG_NOTICE([Configuring CG-DWF for SSE architecture])
		;;
	bluelight)
		AC_MSG_NOTICE([Configuring CG-DWF for BlueLight architecture])
		;;
	altivec)
		AC_MSG_NOTICE([Configuring CG-DWF for AltiVec architecture])
		;;
	*)
		AC_MSG_ERROR([Unknown cpu architecture. Use --enable-cg-dwf=<arch> to select a known one])
		;;
	esac

        AC_MSG_NOTICE([Configuring cg-dwf])
	AC_CONFIG_SUBDIRS(other_libs/cg-dwf)
        AC_MSG_NOTICE([Finished configuring cg-dwf])
	AC_SUBST(CG_DWF_DIR, [cg-dwf])
	AC_DEFINE([BUILD_CG_DWF],[1],[Build CG DWF])
	AC_SUBST(CG_DWF_ENABLED,[yes])
else
	AC_SUBST(CG_DWF_DIR)
	AC_MSG_NOTICE([Not building CG DWF])
fi


 
AM_CONDITIONAL(BUILD_CG_DWF,
	[test "x${build_cg_dwf}x" = "xyesx" ])

AM_CONDITIONAL(BUILD_CG_DWF_SSE,
	[test "x${cg_dwf_cpu_arch}x" = "xssex" \
         -a "x${build_cg_dwf}x"="xyesx"])

AM_CONDITIONAL(BUILD_CG_DWF_BLUELIGHT,
	[test "x${cg_dwf_cpu_arch}x" = "xbluelightx" \
	 -a "x${build_cg_dwf}x"="xyesx" ])

AM_CONDITIONAL(BUILD_CG_DWF_ALTIVEC,
	[test "x${cg_dwf_cpu_arch}x" = "xaltivecx" \
	 -a "x${build_cg_dwf}x"="xyesx" ])

if [ test "x${cg_dwf_lowmem}x" = "xyesx" ];
then
	AC_MSG_NOTICE([Configuring CG-DWF in Lower Memory Mode])
	AC_DEFINE([CHROMA_USE_CG_DWF_LOWMEM],[1], [Switch in Lower Memory File])
else
	AC_MSG_NOTICE([Configuring CG-DWF for Higher Memory Mode])
fi

if [ test "x${cg_do_one_restart}x" = "xyesx" ];
then
	 AC_MSG_NOTICE([Enabling restart in CG Syssolvers ])	
	 AC_DEFINE([CHROMA_DO_ONE_CG_RESTART],[1],[Enable Restart in linop syssolvers])
fi
 

dnl ************************************************************************
dnl **** Bagel wilson dslash                                            ****
dnl ************************************************************************
if test "X${bagel_wilson_dslash_enabled}X" = "XyesX";
then
   AC_MSG_NOTICE( [Configuring with BAGEL Wilson Dslash] )   
   AC_DEFINE([BUILD_PAB_WILSON_DSLASH],[], [ Use Peter Boyles BAGEL Wilson Dslash library])
   if test "X${BAGEL_WILSON_DSLASH_HOME}X" = "XyesX";
   then 
     AC_PATH_PROG([BAGEL_WILSON_DSLASH_CONFIG], [wfm-config], [])
   else
     AC_PATH_PROG([BAGEL_WILSON_DSLASH_CONFIG], [wfm-config], [], [${BAGEL_WILSON_DSLASH_HOME}/bin:${PATH}])
   fi     

   if test "X${BAGEL_WILSON_DSLASH_CONFIG}X" = "XX"; 
   then
     AC_MSG_WARN([wfm-config script for bagel_wilson_dslash not found])
     AC_SUBST(BAGEL_WILSON_DSLASH_CXXFLAGS, ${BAGEL_WILSON_DSLASH_CXXFLAGS})
     AC_SUBST(BAGEL_WILSON_DSLASH_LDFLAGS, ${BAGEL_WILSON_DSLASH_LDFLAGS})
     AC_SUBST(BAGEL_WILSON_DSLASH_LIBS, ${BAGEL_WILSON_DSLASH_LIBS})
   else 
     AC_MSG_NOTICE([Found bagel_wilson_dslash configuration program ${BAGEL_WILSON_DSLASH_CONFIG}])
     AC_SUBST(BAGEL_WILSON_DSLASH_CXXFLAGS, "`${BAGEL_WILSON_DSLASH_CONFIG} --cxxflags`")
     AC_MSG_NOTICE([bagel_wilson_dslash compile flags: ${BAGEL_WILSON_DSLASH_CXXFLAGS}])
     AC_SUBST(BAGEL_WILSON_DSLASH_LDFLAGS, "`${BAGEL_WILSON_DSLASH_CONFIG} --ldflags`")
     AC_MSG_NOTICE([bagel_wilson_dslash link flags : ${BAGEL_WILSON_DSLASH_LDFLAGS}])
     AC_SUBST(BAGEL_WILSON_DSLASH_LIBS,     "`${BAGEL_WILSON_DSLASH_CONFIG} --libs`")
     AC_MSG_NOTICE([bagel_wilson_dslash libraries flags: ${BAGEL_WILSON_DSLASH_LIBS}])
   fi
else
   AC_SUBST(BAGEL_WILSON_DSLASH_CXXFLAGS, "")
   AC_SUBST(BAGEL_WILSON_DSLASH_LDFLAGS, "")
   AC_SUBST(BAGEL_WILSON_DSLASH_LIBS, "")
   AC_MSG_NOTICE( [Not using BAGEL Wilson Dslash] )  
fi

if test "X${bagel_wilson_dslash_enabled}X" = "XyesX";
then
  AC_MSG_CHECKING([if we can compile/link a simple BAGEL Wilson Dslash program])
  PAC_BAGEL_WFM_LINK_CXX_FUNC(
    ${BAGEL_WILSON_DSLASH_CXXFLAGS},
    ${BAGEL_WILSON_DSLASH_LDFLAGS},
    ${BAGEL_WILSON_DSLASH_LIBS},
    [ WilsonArg wil;],
    [ wfm_vec_init(&wil); wfm_vec_end(&wil);],
    [bagel_wilson_dslash_link_ok=yes],
    [bagel_wilson_dslash_link_ok=no]
  )
  if test "X${bagel_wilson_dslash_link_ok}X" = "XyesX"; 
  then 
    AC_MSG_RESULT(yes)
  else 
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot compile/link a program with Bagel Wilson Dslash. Use --with-bagel-wilson-dslash=<dir> to select a working version.])
  fi 

  if test "X${bagel_wilson_dslash_sloppy}X" = "XyesX";
  then 
    AC_MSG_RESULT([Enabling Sloppy Precision in Bagel Wilson Dslash])
    AC_DEFINE([CHROMA_USE_SLOPPY_BAGEL_DSLASH], 
              [1], 
	      [Enable Sloppy Precision in Bagel Wilson Dslash])
  fi
else
  if test "X${bagel_wilson_dslash_sloppy}X" = "XyesX";
  then 
    AC_MSG_WARN([Enabling Sloppy Precision has no effect if not using --with-bagel-wilson-dslash option])
  fi
fi


AM_CONDITIONAL(BUILD_PAB_WILSON_DSLASH, [test "x${bagel_wilson_dslash_enabled}x" = "xyesx" ])


dnl ************************************************************************
dnl **** Bagel clover                                                   ****
dnl ************************************************************************
if test "X${bagel_clover_enabled}X" = "XyesX";
then
   AC_MSG_NOTICE( [Configuring with BAGEL Clover Apply] )   
   AC_DEFINE([BUILD_BAGEL_CLOVER_TERM],[], [ Use the BAGEL Clover Term Apply library])
   if test "X${BAGEL_CLOVER_HOME}X" = "XyesX";
   then 
     AC_PATH_PROG([BAGEL_CLOVER_CONFIG], [bagel-clover-config], [])
   else
     AC_PATH_PROG([BAGEL_CLOVER_CONFIG], [bagel-clover-config], [], [${BAGEL_CLOVER_HOME}/bin:${PATH}])
   fi     

   if test "X${BAGEL_CLOVER_CONFIG}X" = "XX"; 
   then
     AC_MSG_WARN([bagel-clover-config script for bagel_wilson_dslash not found])
     AC_SUBST(BAGEL_CLOVER_DSLASH_CXXFLAGS, ${BAGEL_CLOVER_CXXFLAGS})
     AC_SUBST(BAGEL_CLOVER_DSLASH_LDFLAGS, ${BAGEL_CLOVER_LDFLAGS})
     AC_SUBST(BAGEL_CLOVER_DSLASH_LIBS, ${BAGEL_CLOVER_LIBS})
   else 
     AC_MSG_NOTICE([Found bagel-clover-config configuration program ${BAGEL_CLOVER_CONFIG}])
     AC_SUBST(BAGEL_CLOVER_CXXFLAGS, "`${BAGEL_CLOVER_CONFIG} --cxxflags`")
     AC_MSG_NOTICE([Bagel_Clover compile flags: ${BAGEL_CLOVER_CXXFLAGS}])
     AC_SUBST(BAGEL_CLOVER_LDFLAGS, "`${BAGEL_CLOVER_CONFIG} --ldflags`")
     AC_MSG_NOTICE([Bagel_Clover link flags : ${BAGEL_CLOVER_LDFLAGS}])
     AC_SUBST(BAGEL_CLOVER_LIBS,     "`${BAGEL_CLOVER_CONFIG} --libs`")
     AC_MSG_NOTICE([Bagel Clover libraries flags: ${BAGEL_CLOVER_LIBS}])
   fi
else
   AC_SUBST(BAGEL_CLOVER_CXXFLAGS, "")
   AC_SUBST(BAGEL_CLOVER_LDFLAGS, "")
   AC_SUBST(BAGEL_CLOVER_LIBS, "")
   AC_MSG_NOTICE( [Not using BAGEL Clover] )  
fi


if test "X${bagel_clover_enabled}X" = "XyesX";
then
  AC_MSG_CHECKING([if we can compile/link a simple BAGEL Clover program])
  PAC_BAGEL_CLOVER_LINK_CXX_FUNC(
    ${BAGEL_CLOVER_CXXFLAGS},
    ${BAGEL_CLOVER_LDFLAGS},
    ${BAGEL_CLOVER_LIBS},
    [],
    [bagel_clover(0, 0, 0, 0, 1);],
    [bagel_clover_link_ok=yes],
    [bagel_clover_link_ok=no]
  )
  if test "X${bagel_clover_link_ok}X" = "XyesX"; 
  then 
    AC_MSG_RESULT(yes)
  else 
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot compile/link a program with Bagel Clover. Use --with-bagel-clover=<dir> to select a working version.])
  fi 

fi

AM_CONDITIONAL(BUILD_BAGEL_CLOVER_APPLY, [test "x${bagel_clover_enabled}x" = "xyesx" ])

dnl ************************************************************************
dnl **** SSE Double Prec clover                                          ****
dnl ************************************************************************
if test "X${ssed_clover_enabled}X" = "XyesX";
then
   AC_MSG_NOTICE( [Building with SSED Clover Apply] )   
   AC_DEFINE([BUILD_SSED_CLOVER_TERM],[], [ Use SSED Clover Term Apply ])
else
   AC_MSG_NOTICE( [Not using SSE DP Clover] )  
fi

AM_CONDITIONAL(BUILD_SSED_CLOVER_APPLY, [test "x${ssed_clover_enabled}x" = "xyesx" ])


dnl ************************************************************************
dnl **** Use Lapack                                                     ****
dnl ************************************************************************
case "${LAPACK_BINDING}" in
lapack)
	AC_MSG_NOTICE([Configuring to use Lapack.])
	AC_DEFINE(BUILD_LAPACK,[],[Using Lapack])
        AC_SUBST(LAPACK_ENABLED,[yes])
	;;
*)
	AC_MSG_NOTICE([Will not use Lapack])
        AC_SUBST(LAPACK_ENABLED,[no])
	;;
esac

AM_CONDITIONAL(BUILD_LAPACK, [ test "X${LAPACK_ENABLED}X" = "XyesX" ])

dnl ************************************************************************
dnl **** Use optimized eigcg                                            ****  
dnl ************************************************************************   
case "${opt_eigcg_enabled}" in
yes)
        AC_MSG_NOTICE([Configuring to use optimized eigcg])
        AC_DEFINE(BUILD_OPT_EIGCG,[],[Configuring to use optimized eigcg])
        AC_SUBST(OPT_EIGCG_ENABLED,[yes])
        ;;
*)
        AC_MSG_NOTICE([Will not build optimized eigCG.])
        AC_SUBST(OPT_EIGCG_ENABLED,[no])
        ;;
esac

AM_CONDITIONAL(BUILD_OPT_EIGCG, [ test "X${OPT_EIGCG_ENABLED}X" = "XyesX" ])


dnl ************************************************************************
dnl **** Level3 Asqtad inverter stuff                                   ****
dnl ************************************************************************
case "$enable_asqtad_level3_inverter" in
 yes)
        AC_MSG_NOTICE([Configuring Asqtad level3 Inverter Wrapper])
	AC_DEFINE([BUILD_CPS_ASQTAD_INVERTER],[],[Build CPS ASQTAD INVERTER])
	;;
 *) 
        AC_MSG_NOTICE([Not Configuring Asqtad level3 Inverter Wrapper])
        ;;	
esac

AM_CONDITIONAL(BUILD_ASQTAD_LEVEL3_INVERTER_WRAPPER,
  [test "x${enable_asqtad_level3_inverter}x" = "xyesx" ])


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
if test "X${gmp_enabled}X" = "XyesX";
then 
        AC_MSG_NOTICE( [Enabling GNU Multi-Precision Library] )
	AC_DEFINE([BUILD_GMP_REMEZ],[],[ Build GMP Remez code ])
	if test "X${GMP_HOME}X" != "XyesX";
	then
           dnl if we were given an install path set up the flags
	   AC_SUBST(GMP_CXXFLAGS, "-I${GMP_HOME}/include")
	   AC_SUBST(GMP_LDFLAGS, "-L${GMP_HOME}/lib")
	   AC_SUBST(GMP_LIBS, "-lgmp")
	else
	   dnl if we were not given an install path use existing env variables/standard include
           dnl link paths
           AC_SUBST(GMP_CXXFLAGS, ${GMP_CXXFLAGS})
	   AC_SUBST(GMP_LDFLAGS, ${GMP_LDFLAGS})

	   dnl if stuff is on the standard include and lib paths we still need the library
	   dnl invocation
	   if test "X${GMP_LIBS}X" = "XX";
	   then 
	     AC_SUBST(GMP_LIBS, "-lgmp")
	   else 
	     AC_SUBST(GMP_LIBS, ${GMP_LIBS})
           fi
	fi
	AC_MSG_NOTICE([GMP CXXFLAGS: ${GMP_CXXFLAGS}])
	AC_MSG_NOTICE([GMP LDFLAGS: ${GMP_LDFLAGS}])
	AC_MSG_NOTICE([GMP LIBS: ${GMP_LIBS}])
        AC_MSG_CHECKING([if we can compile/link a simple GMP program])	
	dnl check we can link against gmp
	PAC_GMP_LINK_CXX_FUNC(
	 [ ${GMP_CXXFLAGS} ],
	 [ ${GMP_LDFLAGS} ],
	 [ ${GMP_LIBS} ],
	 [ mpf_t x; ],
	 [ mpf_init(x); ],
         [gmp_working="yes"],
	 [gmp_working="no"]
        )
	
	if test "X${gmp_working}X" = "XyesX";
	then 
	  AC_MSG_RESULT(yes)
	else
	  AC_MSG_RESULT(no)
	  AC_MSG_ERROR([Couldn't link a simple GMP program. Use --with-gmp=<dir> to set a working one])
	fi
else
	AC_MSG_NOTICE( [Not building GMP Remez code ] )
	AC_SUBST(GMP_CXXFLAGS, "")
	AC_SUBST(GMP_LDFLAGS,"")
	AC_SUBST(GMP_LIBS, "")
fi

AM_CONDITIONAL(BUILD_GMP_REMEZ,
  [test "x${gmp_enabled}x" = "xyesx" ])

dnl QMT Threads stuff
if test "X${qmt_enabled}X" == "XyesX"; 
then
	AC_MSG_NOTICE([Configuring QMT Threading])
	AC_DEFINE([CHROMA_USE_QMT_THREADS], [1], [ Use QMT Threads library ])
	AC_SUBST(QMT_CXXFLAGS, "-I${QMT_HOME}/include")
	AC_SUBST(QMT_LDFLAGS, "-L${QMT_HOME}/lib")
	AC_SUBST(QMT_LIBS, "-lqmt -lpthread -lm")
else 
	AC_SUBST(QMT_CXXFLAGS, "")
	AC_SUBST(QMT_LDFLAGS, "")
	AC_SUBST(QMT_LIBS, "")

fi
AM_CONDITIONAL(BUILD_QMT, [test "x${qmt_enabled}x" = "xyesx" ])

case "${chroma_testcase_runner}" in
trivial)
	AC_MSG_NOTICE([Using Trivial Testcase Runner: trivial_runner.sh])
	chroma_testcase_runner_script="trivial_runner.sh"
	;;
6n_mpirun_rsh)
	AC_MSG_NOTICE([Using 6n_mpirun_rsh Testcase Runner: 6n_mpirun_rsh_runner.sh])
	chroma_testcase_runner_script="6n_mpirun_rsh_runner.sh"
	;;
7n_mpirun_rsh)
	AC_MSG_NOTICE([Using 7n_mpirun_rsh Testcase Runner: 7n_mpirun_rsh_runner.sh])
	chroma_testcase_runner_script="7n_mpirun_rsh_runner.sh"
	;;
cobalt)
	AC_MSG_NOTICE([Using the BG/L Cobalt Testcase Runner: cobalt_runner.sh])
	chroma_testcase_runner_script="cobalt_runner.sh"
	;;
*)
	AC_MSG_WARN([Unknown value for chroma_testcase_runner. Will use trivial])
	chroma_testcase_runner_script="trivial_runner.sh"
	;;
esac



dnl
dnl Always build qdp-lapack. If no --enable-lapack option is given,
dnl it will build a dummy lib with stubs. If the enable option is
dnl given, it will build something that can be used by chroma
dnl and not fry when call
dnl
AC_MSG_NOTICE([Configuring qdp-lapack])
AC_CONFIG_SUBDIRS(other_libs/qdp-lapack)
AC_MSG_NOTICE([Finished configuring qdp-lapack])
AC_SUBST(QDP_LAPACK_DIR, [qdp-lapack])

dnl
dnl ************************************************************
dnl MDWF stuff
dnl ************************************************************
dnl
if [ test "x${build_mdwf}x" = "xyesx" ];
then
  
  AC_MSG_NOTICE([Configuring to use MDWF DWF Solver package])

  dnl sanity check
  if [ test "x${build_cg_dwf}x" = "xyesx" ];
  then
    AC_MSG_ERROR([Cannot have both MDWF and CG-DWF enabled together])
  fi

  dnl set a preprocessor macro
  AC_DEFINE([BUILD_MDWF],[], [ Build with the MDWF (Moebious DFW) library ])

  dnl find the config package
  if test "X${mdwf_path}X" = "XX" ; then
     AC_PATH_PROG(MDWF_CONFIG, [mdwf-config], [])
  else
     AC_PATH_PROG(MDWF_CONFIG, [mdwf-config], [], [${mdwf_path}/bin:${PATH}])
  fi

  if test "X${MDWF_CONFIG}X" = "XX" ; then
    AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
  fi

  dnl set some flags for the makefiles
  AC_MSG_NOTICE([Found MDWF configuration program ${MDWF_CONFIG}])
  AC_SUBST(MDWF_CXXFLAGS, "`${MDWF_CONFIG} --cflags`")
  AC_MSG_NOTICE([MDWF compile flags: ${MDWF_CXXFLAGS}])
  AC_SUBST(MDWF_LDFLAGS,  "`${MDWF_CONFIG} --ldflags`")
  AC_MSG_NOTICE([MDWF linking flags: ${MDWF_LDFLAGS}])
  AC_SUBST(MDWF_LIBS,     "`${MDWF_CONFIG} --libs`")
  AC_MSG_NOTICE([MDWF libraries flags: ${MDWF_LIBS}])
  
  dnl for the chroma-config file
  AC_SUBST(MDWF_ENABLED, [yes])
fi

dnl set a conditional in automake for MDWF Specific files etc.
AM_CONDITIONAL(BUILD_MDWF,
	[test "x${build_mdwf}x" = "xyesx" ])

##################################
# Check for programs to build docs
##################################
AC_CHECK_PROG(LATEX, latex, latex, true)
AC_CHECK_PROG(DVIPS, dvips, dvips, true)
AC_CHECK_PROG(DVIPDF, dvipdf, dvipdf, true)
AC_CHECK_PROG(DOXYGEN, doxygen, doxygen, true)
AC_PATH_PROG(XMLDIFF, [xmldiff],[xmldiff],[$PATH:/usr/local/bin])

if test "X${XMLDIFF}X" = "XxmldiffX";
then 
     AC_MSG_WARN([xmldiff not found in your PATH. make xcheck may fail])
fi

###################################
# 
# QUDA Linkage etc
# 
###################################
if test "X${quda_enabled}X" = "XyesX" ;
then
   AC_MSG_NOTICE([ Using QUDA Library in ${QUDA_HOME} ])

   case "$host" in 
   x86_64*)
        CUDA_LIBDIR="${CUDA_HOME}/lib64"	
	;;
   *)	
   CUDA_LIBDIR="${CUDA_HOME}/lib"
   ;;	
   esac
	
   AC_MSG_NOTICE([ Linking to CUDA Runtime in ${CUDA_LIBDIR}])

   AC_DEFINE([BUILD_QUDA],[], [Build with the BU QUDA Library])  
   AC_SUBST(QUDA_CXXFLAGS, "-I${QUDA_HOME} -I${QUDA_HOME}/include -I${CUDA_HOME}/include" )
   AC_SUBST(QUDA_LDFLAGS,  "-L${CUDA_LIBDIR} -L${QUDA_HOME} -L${QUDA_HOME}/lib" )
   AC_SUBST(QUDA_LIBS,     "-lquda -lcudart" )
fi

AM_CONDITIONAL(BUILD_QUDA,
	[test "x${quda_enabled}x" = "xyesx" ])

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
AC_CONFIG_FILES(docs/Makefile)
AC_CONFIG_FILES(other_libs/Makefile)

dnl run all the scripts through AC CONFIG files to put
dnl them in the build directory. We can even potentially install them
dnl if we want later. Also this allows us to configure complicated runners 
dnl later with extre --enable options etc
AC_CONFIG_FILES(scripts/test_runners/trivial_runner.sh)
AC_CONFIG_FILES(scripts/test_runners/6n_mpirun_rsh_runner.sh)
AC_CONFIG_FILES(scripts/test_runners/7n_mpirun_rsh_runner.sh)
AC_CONFIG_FILES(scripts/test_runners/cobalt_runner.sh)

dnl Copy - dont symlink - the chosen testcase runner to the top builddir for
dnl for the fixture runner script
AC_CONFIG_COMMANDS([RUN],
		   [ cp scripts/test_runners/${the_testcase_runner} ./RUN ; chmod u+x ./RUN ],[the_testcase_runner=${chroma_testcase_runner_script}])

AC_OUTPUT
