#
#
# AC_SIRIUS
#
#
#
dnl @synopsis AC_SIRIUS
dnl
dnl This macro test if sirius (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef SIRIUS
dnl     #include "sirius.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_SIRIUS],[

AC_MSG_NOTICE([=== checking for SIRIUS ===])

AM_CONDITIONAL(HAVE_SIRIUS,true)

AC_ARG_WITH(sirius,
        [  --with-sirius=DIR      Location of SIRIUS library],
        [SIRIUS_LDFLAGS="-L$withval/lib";
         SIRIUS_LIBS="-lsplitdoubles /opt/gcc/4.9.0/snos/lib64/libstdc++.a";
         SIRIUS_CPPFLAGS="-I$withval/include";],
        [with_sirius=no])

if test "x$with_sirius" == "xno"; then

   AM_CONDITIONAL(HAVE_SIRIUS,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lsplitdoubles /opt/gcc/4.9.0/snos/lib64/libstdc++.a"
    LDFLAGS="$LDFLAGS $SIRIUS_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $SIRIUS_CPPFLAGS"
    
    # Check for the SIRIUS library and headers
    dnl AC_TRY_COMPILE([struct obd_uuid {char uuid[40];};int fd, num_ost;struct obd_uuid uuids[1024];],
    dnl        [llapi_lov_get_uuids(fd, uuids, &num_ost);],
    dnl        [LUSTRE_LIBS="-llustreapi"],
    dnl        [AM_CONDITIONAL(HAVE_LUSTRE,false)])
    
    dnl LIBS="$save_LIBS"
    dnl LDFLAGS="$save_LDFLAGS"
    dnl CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(SIRIUS_LIBS)
    AC_SUBST(SIRIUS_LDFLAGS)
    AC_SUBST(SIRIUS_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_SIRIUS_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_SIRIUS,1,[Define if you have SIRIUS.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_SIRIUS
