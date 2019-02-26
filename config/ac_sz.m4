#
#
# AC_SZ
#
#
#
dnl @synopsis AC_SZ
dnl
dnl This macro test if sz (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef SZ
dnl     #include "sz.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_SZ],[

AC_MSG_NOTICE([=== checking for SZ ===])

AM_CONDITIONAL(HAVE_SZ,true)

AC_ARG_WITH(sz,
        [  --with-sz=DIR      Location of SZ library],
        [SZ_LDFLAGS="-L$withval";
         SZ_LIBS="-L$withval/lib -lSZ -lzlib -lzstd";
         SZ_CPPFLAGS="-I$withval/include";],
        [with_sz=no])

if test "x$with_sz" == "xno"; then

   AM_CONDITIONAL(HAVE_SZ,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lSZ"
    LDFLAGS="$LDFLAGS $SZ_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $SZ_CPPFLAGS"
    
    # Check for the SZ library and headers
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(SZ_LIBS)
    AC_SUBST(SZ_LDFLAGS)
    AC_SUBST(SZ_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_SZ_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_SZ,1,[Define if you have SZ.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_SZ
