dnl George T. Fleming, 13 February 2003
dnl
dnl Descended originally from mpich-1.2.4/mpe
dnl
dnl PAC_QDPXX_LINK_CXX_FUNC(
dnl   QDPXX_CXXFLAGS,
dnl   QDPXX_LDFLAGS,
dnl   QDPXX_LIBS,
dnl   QDPXX_VARS,
dnl   QDPXX_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  QDPXX_CXXFLAGS for the necessary includes paths (-I)
dnl  QDPXX_LDFLAGS  for the necessary library search paths (-L)
dnl  QDPXX_LIBS     for the libraries (-l<lib> etc)
dnl  QDPXX_VARS     for the declaration of variables needed
dnl                 to call QDPXX_FUNC code fragment
dnl  QDPXX_FUNC     for the body of a QDP++ function call or even general code
dnl                 fragment on which to run a compile/link test.
dnl                 If QDPXX_VARS and QDPXX_FUNC are empty, a basic test
dnl                 of compiling and linking a QDP++ program is run.
dnl
AC_DEFUN(
  PAC_QDPXX_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl - so input variables can be CXXFLAGS, LDFLAGS or LIBS
    pac_QDPXX_CXXFLAGS="$1"
    pac_QDPXX_LDFLAGS="$2"
    pac_QDPXX_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_QDPXX_CXXFLAGS"
    LDFLAGS="$LDFLAGS $pac_QDPXX_LDFLAGS"
    LIBS="$LIBS $pac_QDPXX_LIBS"
    AC_TRY_LINK(
      [
        #include <qdp.h>
        using namespace QDP;
      ], [
        int argc ; char **argv ;
        // Turn on the machine
        QDP_initialize(&argc, &argv) ;
        // Create the layout
        const int foo[] = {2,2,2,2} ;
        multi1d<int> nrow(Nd) ;
        nrow = foo ; // Use only Nd elements
        Layout::setLattSize(nrow) ;
        Layout::create() ;
        $4 ;
        $5 ;
        QDP_finalize() ;
      ],
      [pac_qdpxx_working=yes],
      [pac_qdpxx_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_qdpxx_working}X"="XyesX";
    then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)

AC_DEFUN(
  PAC_BAGEL_WFM_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_BAGEL_WFM_CFLAGS="$1"
    pac_BAGEL_WFM_LDFLAGS="$2"
    pac_BAGEL_WFM_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_BAGEL_WFM_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_BAGEL_WFM_LDFLAGS"
    LIBS="$LIBS $pac_BAGEL_WFM_LIBS"
    AC_TRY_LINK(
      [
        #include "wfm.h"
      ],
      [
        int argc ; char **argv ;
	$4
	$5
      ],
      [pac_bagel_wfm_working=yes],
      [pac_bagel_wfm_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_bagel_wfm_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)

AC_DEFUN(
  PAC_GMP_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_GMP_CFLAGS="$1"
    pac_GMP_LDFLAGS="$2"
    pac_GMP_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_GMP_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_GMP_LDFLAGS"
    LIBS="$LIBS $pac_GMP_LIBS"
    AC_TRY_LINK(
      [
        #include <gmp.h>
      ],
      [
        int argc ; char **argv ;
	$4
	$5
      ],
      [pac_gmp_working=yes],
      [pac_gmp_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_gmp_working}X" = "XyesX" ; 
    then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)

AC_DEFUN(
  PAC_BAGEL_CLOVER_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_BAGEL_CLOVER_CFLAGS="$1"
    pac_BAGEL_CLOVER_LDFLAGS="$2"
    pac_BAGEL_CLOVER_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_BAGEL_CLOVER_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_BAGEL_CLOVER_LDFLAGS"
    LIBS="$LIBS $pac_BAGEL_CLOVER_LIBS"
    AC_TRY_LINK(
      [
        #include <bagel_clover.h>
      ],
      [
        int argc ; char **argv ;
	$4
	$5
      ],
      [pac_bagel_clover_working=yes],
      [pac_bagel_clover_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_bagel_clover_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)

