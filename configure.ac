AC_INIT(szin,1.2,edwards@jlab.org)
AC_CONFIG_AUX_DIR(./config)
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_CXX(g++ CC)

dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl Balint Joo 13/12/2002 
dnl Options for finding out where MPI/GM/QMP is
dnl folowing George Flemings QMP mods

dnl --enable-fermion-type
AC_ARG_ENABLE(fermion-type,
	      AC_HELP_STRING([--enable-fermion-type=<fermion type>],
	        [ Set the fermion type (wilson or staggered) default is wilson ]),
	      [ ac_fermion_type="$enable_fermion_type" ],
	      [ ac_fermion_type="wilson" ] ) 


dnl --with-qdp-cxxflags
AC_ARG_WITH( qdp-cxxflags,
  AC_HELP_STRING([--with-qdp-cxxflags=QDP_CXXFLAGS],
    [To pass include path to QDP header files via compiler flag
      -I]),
  [QDP_CXXFLAGS="$with_qdp_cxxflags"])
AC_SUBST(QDP_CXXFLAGS)

dnl --with-qdp-ldflags
AC_ARG_WITH( qdp-ldflags,
  AC_HELP_STRING([--with-qdp-ldflags=QDP_LDFLAGS],
    [To pass QDP library search path via compiler flag -L]),
  [QDP_LDFLAGS="$with_qdp_ldflags"])
AC_SUBST(QDP_LDFLAGS)

dnl --with-qdp-libs 
AC_ARG_WITH( qdp-libs,
  AC_HELP_STRING([--with-qdp-libs=QDP_LIBS],
    [To pass optional list of QDP libraries to linker via compiler flag -l]),
  [QDP_LIBS="$with_qdp_libs"])
AC_SUBST(QDP_LIBS)

dnl now the qmp-flags
dnl --with-qmp-cxxflags
AC_ARG_WITH( qmp-cxxflags,
  AC_HELP_STRING([--with-qmp-cxxflags=QMP_CXXFLAGS],
    [To pass optional include path to QMP header files via compiler flag
      -I]),
  [QMP_CXXFLAGS="$with_qmp_cxxflags"] )
AC_SUBST(QMP_CXXFLAGS)

dnl --with-qmp-ldflags
AC_ARG_WITH( qmp-ldflags, 
  AC_HELP_STRING([--with-qmp-ldflags=QMP_LDFLAGS],
    [To pass optional link path to QMP libraries via compiler flag -L]),
    [QMP_LDFLAGS="$with_qmp_ldflags"])
AC_SUBST(QMP_LDFLAGS)

dnl --with-qmp-libs
AC_ARG_WITH( qmp-libs,
  AC_HELP_STRING([-with-qmp-libs=QMP_LIBS],
    [To pass optional libraries to QMP libraries via compiler flag -l]),
    [QMP_LIBS="$with_qmp_libs"])
AC_SUBST(QMP_LIBS)

dnl Now the QMP Comms Flags

dnl --with-qmp-comms-type
AC_ARG_WITH(qmp-comms-type,
  AC_HELP_STRING([--with-qmp-comms-type=TYPE],
    [Build QMP on top of MPI or GM.  If MPI over GM is desired, choose MPI.
      Default is MPI.]),
  [case "${with_qmp_comms_type}" in
    MPI) ;;
    mpi) with_qmp_comms_type=MPI ;;
    GM)  ;;
    gm) with_qmp_comms_type=GM ;;
    *) AC_MSG_ERROR([bad value "${with_qmp_comms_type}"
      for --with-qmp-comms-type]) ;;
    esac],
  [with_qmp_comms_type=MPI])

dnl --with-qmp-comms-cxxflags
AC_ARG_WITH( qmp-comms-cxxflags,
  AC_HELP_STRING([--with-qmp-comms-cxxflags=QMP_COMMS_CXXFLAGS],
    [To pass optional include path to comms header files via compiler flag
      -I]),
  [QMP_COMMS_CXXFLAGS="$with_qmp_comms_cxxflags"] )
AC_SUBST(QMP_COMMS_CXXFLAGS)

dnl --with-qmp-comms-ldflags
AC_ARG_WITH( qmp-comms-ldflags,
  AC_HELP_STRING([--with-qmp-comms-ldflags=QMP_COMMS_LDFLAGS],
    [To pass optional library search path to comms libraries via compiler
      flag -L]),
  [QMP_COMMS_LDFLAGS="$with_qmp_comms_ldflags"])
AC_SUBST(QMP_COMMS_LDFLAGS)

dnl --with-qmp-comms-libs
AC_ARG_WITH( qmp-comms-libs,
  AC_HELP_STRING([--with-qmp-comms-libs=QMP_COMMS_LIBS],
    [To pass optional list of comms libraries to linker via compiler flag -l]),
  [QMP_COMMS_LIBS="$with_qmp_comms_libs"])
AC_SUBST(QMP_COMMS_LIBS)


dnl
dnl Check the fermion type
dnl
case $ac_fermion_type in
  wilson) 
    AC_MSG_NOTICE([ Configuring SZIN++ with WILSON-like fermions ])
    ;;
  staggered)
    AC_MSG_NOTICE([ Configuring SZIN++ with STAGGERED-like fermions ])
    ;;
  *)
    AC_MSG_ERROR([ must specify either wilson or staggered with --enable-fermion-type])
    ;;
esac

dnl Try to compile a QDP++ program to check the --with options
AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])

dnl This function is defined in acinclude.m4
dnl it tries to try to compile the program 
PAC_QDP_LINK_CXX_FUNC( ${QDP_CXXFLAGS}, ${QDP_LDFLAGS}, ${QDP_LIBS}, dnl
                       ${QMP_CXXFLAGS}, ${QMP_LDFLAGS}, ${QMP_LIBS}, dnl
		       ${QMP_COMMS_CXXFLAGS}, ${QMP_COMMS_LDFLAGS},  dnl
                       ${QMP_COMMS_LIBS}, , , [qdp_link_ok=yes], [qdp_link_ok=no] )
if test "X${qdp_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR( [ Cannot compile/link a program with QDP please check/set 
   the values of QDP_CXXFLAGS, QDP_LDFLAGS, QDP_LIBS, QMP_CXXFLAGS, QMP_LDFLAGS
   QMP_LIBS, QMP_COMMS_CXXFLAGS, QMP_COMMS_LDFLAGS, QMP_COMMS_LIBS ] )
fi 

dnl Set a conditional to construct the list of sources for libszin.a
AM_CONDITIONAL(BUILD_WILSON_LIBS, test "X${ac_fermion_type}X" = "XwilsonX" )

dnl Make the Makefiles
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(mainprogs/Makefile)
AC_CONFIG_FILES(mainprogs/tests/Makefile)
AC_CONFIG_FILES(mainprogs/main/Makefile)
AC_OUTPUT
