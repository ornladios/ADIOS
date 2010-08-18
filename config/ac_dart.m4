#
#
# AC_PROG_DART
#
# Test for DART installation
# and set $DART to the correct value.
#
#
dnl @synopsis AC_PROG_DART
dnl
dnl This macro test if DART is installed. If DART
dnl is installed, it set $DART to the right value
dnl
dnl @version 1.0
dnl @author Norbert Podhorszki, pnorbert@ornl.gov
dnl
AC_DEFUN([AC_PROG_DART],[

AM_CONDITIONAL(HAVE_DART,true)

AC_ARG_WITH(dart,
        [AS_HELP_STRING([--with-dart=DIR],
           [Build the DART transport method. Point to the DART installation.])],
        [DART_LDFLAGS="-L$withval/lib";
         DART_CPPFLAGS="-I$withval/include";],
        [with_dart=check])

dnl If --without-dart was given set HAVE_DART to false and do nothing more
if test "x$with_dart" == "xno"; then

   AM_CONDITIONAL(HAVE_DART,false)

else

    dnl allow args --with-dart incdir and --with-dart-libdir
    AC_ARG_WITH(dart-incdir,
                [  --with-dart-incdir=<location of dart includes>],
                [DART_INCDIR=$withval
                 with_dart=detailed])
    
    AC_ARG_WITH(dart-libdir,
                [  --with-dart-libdir=<location of dart library>],
                [DART_LIBDIR=$withval
                 with_dart=detailed])
    
    
    dnl If we know DART_DIR, then we can know DART_INCDIR.
    dnl We don't overwrite DART_INCDIR.
    if test -n "${DART_DIR}" -a -z "${DART_INCDIR}"; then
            DART_INCDIR="${DART_DIR}/include";
    dnl We may have DART denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DART}" -a -d "${DART}"; then
            DART_INCDIR="${DART}/include"
    fi

    dnl If we know DART_DIR, then we can know DART_LIBDIR.
    dnl We don't overwrite DART_LIBDIR.
    if test -n "${DART_DIR}" -a -z "${DART_LIBDIR}"; then
            DART_LIBDIR="${DART_DIR}/lib";
    dnl We may have DART denoting the dir (e.g. on ewok BUT on franklin it contains all flags)
    elif test -n "${DART}" -a -d "${DART}"; then
            DART_LIBDIR="${DART}/lib"
    fi

    dnl Add "-I" to DART_INCDIR.
    if test -n "${DART_INCDIR}"; then
            DART_CPPFLAGS="-I${DART_INCDIR}"
    else
            ac_dart_ok=no
    fi

    dnl Add "-L" to DART_LIBDIR.
    if test -n "${DART_LIBDIR}"; then
            DART_LDFLAGS="-L${DART_LIBDIR}"
    else
            ac_dart_ok=no
    fi

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -ldart2 -lspaces"
    LDFLAGS="$LDFLAGS $DART_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $DART_CPPFLAGS"
    
    if test -z "${HAVE_DART_TRUE}"; then
            AC_CHECK_HEADERS(dart.h,
                    ,
                    [AM_CONDITIONAL(HAVE_DART,false)])
    fi
    
    # Check for the Mini-XML library and headers
    AC_TRY_COMPILE([#include "dart.h"],
            [int err; err = dart_init(1,1);],
            [DART_LIBS="-ldart2 -lspaces"],
            [AM_CONDITIONAL(HAVE_DART,false)])
    
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(DART_LIBS)
    AC_SUBST(DART_LDFLAGS)
    AC_SUBST(DART_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_DART_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_DART,1,[Define if you have the DART.])],[$1])
            :
    else
            $2
            :
    fi

fi
])dnl AC_DART
