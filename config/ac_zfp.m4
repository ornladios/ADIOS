#
#
# AC_ZFP
#
#
#
dnl @synopsis AC_ZFP
dnl
dnl This macro tests if ZFP is to be used.
dnl Use in C code:
dnl     #ifdef ZFP
dnl     #include "zfp.h"
dnl     #endif
dnl
dnl @version 0.1
dnl @author Eric Suchyta
dnl

AC_DEFUN([AC_ZFP],[

AC_MSG_NOTICE([=== checking for ZFP ===])
AM_CONDITIONAL(HAVE_ZFP,true)
AM_CONDITIONAL(BUILD_ZFP,false)


AC_ARG_WITH(zfp,
        [  --with-zfp=DIR      Location of ZFP library],
        [],
        [with_zfp=builtin])


if test "x$with_zfp" == "xno"; then
   AM_CONDITIONAL(HAVE_ZFP,false)

elif test "x$with_zfp" == "xbuiltin"; then
   AM_CONDITIONAL(BUILD_ZFP,true)

else
    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"


	if test -z "${ZFP_LIBS}"; then
		ZFP_LIBS="-lzfp"
	fi

	if test -z "${ZFP_LDFLAGS}"; then
		ZFP_LDFLAGS="-L$with_zfp/lib"
	fi

	if test -z "${ZFP_CPPFLAGS}"; then
		ZFP_CPPFLAGS="-I$with_zfp/inc"
	fi


    LIBS="$LIBS $ZFP_LIBS"
    LDFLAGS="$LDFLAGS $ZFP_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ZFP_CPPFLAGS"

    if test -z "${HAVE_ZFP_TRUE}"; then
           AC_CHECK_HEADERS(zfp.h, , [AM_CONDITIONAL(HAVE_ZFP,false)])
    fi

    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    AC_SUBST(ZFP_LIBS)
    AC_SUBST(ZFP_LDFLAGS)
    AC_SUBST(ZFP_CPPFLAGS)

    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ZFP_TRUE}"; then
            ifelse([$1], , [AC_DEFINE(HAVE_ZFP,1,[Define if you have ZFP.])], [$1])
            :
    else
            $2
            :
    fi
fi
])
dnl AC_ZFP
