dnl Balint Joo, 13/12/2002
dnl
dnl Stole this from mpich-1.2.4/mpe
dnl
dnl PAC_QDP_LINK_CXX_FUNC(QDP_CXXFLAGS, QDP_LDFLAGS, QDP_LIBS,
dnl                       QMP_CXXFLAGS, QMP_LDFLAGS, QMP_LIBS,
dnl                       QMP_COMMS_CXXFLAGS, QMP_COMMS_LFLAGS, QMP_COMMS_LIBS,
dnl                       QDP_VARS, QDP_FUNC,
dnl                       [action if working], [action if not working] )
dnl
dnl  QDP_CXXFLAGS       is the include option (-I) for QDP includes
dnl  QDP_LDFLAGS        is the link path (-L) option for QDP libraries
dnl  QDP_LIBS           is the library (-l) option for QDP libraries
dnl  QMP_CXXFLAGS       is the include option (-I) for QMP includes
dnl  QMP_LDFLAGS        is the link path (-L) option for QMP libraries
dnl  QMP_LIBS           is the library (-l) option for QMP libaries
dnl  QMP_COMMS_CXXFLAGS is the include option (-I) for QMP_COMMS includes
dnl                     such as MPI 
dnl  CMP_COMMS_LDFLAGS  is the link path (-L) option for the QMP_COMMS library
dnl                     such as MPI
dnl  QMP_COMMS_LIBS     is the library (-l) option for the QMP_COMMS libraries
dnl                     such as MPI
dnl - QMP_VARS    is the the declaration of variables needed to call QMP_FUNC
dnl - QMP_FUNC    is the body of QMP function call to be checked for existence
dnl               e.g.  QMP_VARS="QMP_u32_t foo;"
dnl                     QMP_FUNC="foo = QMP_get_SMP_count();"
dnl               if QMP_FUNC is empty, assume linking with basic MPI program.
dnl               i.e. check if QMP definitions are valid
dnl
AC_DEFUN(PAC_QDP_LINK_CXX_FUNC,[
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_QDP_CXXFLAGS="$1 $4 $7"
    pac_QDP_LDFLAGS="$2 $5 $8"
    pac_QDP_LIBS="$3 $6 $9"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_QDP_CXXFLAGS"
    LDFLAGS="$LDFLAGS $pac_QDP_LDFLAGS"
    LIBS="$LIBS $pac_QDP_LIBS"
    AC_TRY_LINK( [#include <qdp.h>
		 using namespace QDP; ], [
    int argc; char **argv;
    // Turn on the machine
    QDP_initialize(&argc, &argv);
    // Create the layout
    const int foo[] = {2,2,2,2};
    multi1d<int> nrow(Nd);
    nrow = foo;  // Use only Nd elements
    Layout::setLattSize(nrow);
    Layout::create();
    $10;
    $11;
    QDP_finalize(); ], pac_qdp_working=yes, pac_qdp_working=no )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "$pac_qdp_working" = "yes" ; then
       ifelse([$12],,:,[$12])
    else
       ifelse([$13],,:,[$13])
    fi
])
