#
#
# AC_LUSTRE
#
#
#
dnl @synopsis AC_LUSTRE
dnl
dnl This macro test if dmalloc is to be used. 
dnl Use in C code:
dnl     #ifdef DMALLOC
dnl     #include "dmalloc.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, UT
dnl
AC_DEFUN([AC_LUSTRE],[

AC_MSG_NOTICE([=== checking for Lustre ===])

AM_CONDITIONAL(HAVE_LUSTRE,true)

AC_ARG_WITH(lustre,
        [  --with-lustre=DIR      Location of lustre library],
        [LUSTRE_LDFLAGS="-L$withval/lib";
         LUSTRE_LIBS="-llustreapi";
         LUSTRE_CPPFLAGS="-I$withval/include";],
        [with_netcdf=check])

if test "x$with_lustre" == "xno"; then

   AM_CONDITIONAL(HAVE_LUSTRE,false)

else

save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
save_LDFLAGS="$LDFLAGS"
LIBS="$LIBS -llustreapi"
LDFLAGS="$LDFLAGS $LUSTRE_LDFLAGS"
CPPFLAGS="$CPPFLAGS $LUSTRE_CPPFLAGS"

dnl if test -z "${HAVE_DMALLOC_TRUE}"; then
dnl        AC_CHECK_HEADERS(dmalloc.h,
dnl                ,
dnl                [AM_CONDITIONAL(HAVE_DMALLOC,false)])
dnl fi

# Check for the lustre library and headers
dnl AC_TRY_COMPILE([#include "dmalloc.h"],
dnl        [char * s; s=malloc(sizeof(char)*10); free(s);],
dnl        [DMALLOC_LIBS="-llustreapi"],
dnl        [AM_CONDITIONAL(HAVE_LUSTRE,false)])

LIBS="$save_LIBS"
LDFLAGS="$save_LDFLAGS"
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(LUSTRE_LIBS)
AC_SUBST(LUSTRE_LDFLAGS)
AC_SUBST(LUSTRE_CPPFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test -z "${HAVE_LUSTRE_TRUE}"; then
        ifelse([$1],,[AC_DEFINE(HAVE_LUSTRE,1,[Define if you have LUSTRE.])],[$1])
        :
else
        $2
        :
fi
fi
])dnl AC_LUSTRE
