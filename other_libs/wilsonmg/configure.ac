AC_INIT(wilsonmg,1.0.0, sdcohen@uw.edu)
AC_PREREQ([2.57])
AC_CONFIG_AUX_DIR(./config)
AC_CONFIG_HEADERS([include/wilsonmg_config_internal.h])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE

AC_PROG_CC(gcc)
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, [ar])

AC_ARG_WITH(qdpc,
  AC_HELP_STRING(
    [--with-qdpc=DIR], 
    [Use QDPC in directory DOR]
  ), 
  [QDPC_DIR="${with_qdpc}"],
  [QDPC_DIR="not_given"]
)

AC_ARG_WITH(qmp,
  AC_HELP_STRING(
    [--with-qmp=DIR], 
    [Use QMP in directory DIR]
  ), 
  [QMP_DIR="${with_qmp}"],
  [QMP_DIR="not_given"]
)

AC_ARG_WITH(qla,
  AC_HELP_STRING(
    [--with-qla=DIR], 
    [Use QLA in directory DOR]
  ), 
  [QLA_DIR="${with_qla}"],
  [QLA_DIR="not_given"]
)

AC_ARG_WITH(qopqdp,
  AC_HELP_STRING(
    [--with-qopqdp=DIR], 
    [Use QOPQDP in directory DOR]
  ), 
  [QOPQDP_DIR="${with_qopqdp}"],
  [QOPQDP_DIR="not_given"]
)

AC_ARG_WITH(qio,
  AC_HELP_STRING(
    [--with-qio=DIR], 
    [Use QIO for header includes in directory DOR]
  ), 
  [QIO_DIR="${with_qio}"],
  [QIO_DIR="not_given"]
)

dnl if we got here then we have the requisite directories	
AC_SUBST(QOP_MG_CFLAGS, "-I${QMP_DIR}/include -I${QLA_DIR}/include -I${QIO_DIR}/include -I${QDPC_DIR}/include -I${QOPQDP_DIR}/include")

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_OUTPUT
