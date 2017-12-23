#
#
# AC_JSON
#
#
#
dnl @synopsis AC_JSON
dnl
dnl This macro test if json (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef JSON
dnl     #include "json.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_JSON],[

AC_MSG_NOTICE([=== checking for JSON ===])

AM_CONDITIONAL(HAVE_JSON,true)

AC_ARG_WITH(json,
        [  --with-json=DIR      Location of JSON library],
        [JSON_LDFLAGS="-L$withval/lib";
         JSON_LIBS="-ljson-c";
         JSON_CPPFLAGS="-I$withval/include";],
        [with_json=no])

if test "x$with_json" == "xno"; then

   AM_CONDITIONAL(HAVE_JSON,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -ljson $EXTRA_LIBS"
    LDFLAGS="$LDFLAGS $JSON_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $JSON_CPPFLAGS"
    
    # Check for the JSON library and headers
    LIBS="$save_LIBS $EXTRA_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(JSON_LIBS)
    AC_SUBST(JSON_LDFLAGS)
    AC_SUBST(JSON_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_JSON_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_JSON,1,[Define if you have JSON.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_JSON
