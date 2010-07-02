#
#
# AC_PROG_MXML
#
# Test for Mini-XML
# and set $MXML to the correct value.
#
#
dnl @synopsis AC_PROG_MXML
dnl
dnl This macro test if mini-XML is installed. If mini-XML
dnl is installed, it set $MXML to the right value
dnl
dnl @version 1.0
dnl @author Jay Lofstead lofstead@cc.gatech.edu
dnl
AC_DEFUN([AC_PROG_MXML],[

AM_CONDITIONAL(HAVE_MXML,true)

AC_ARG_WITH(mxml,
        [  --with-mxml=DIR      Location of Mini-XML library],
        [MXML_LDFLAGS="-L$withval/lib";
         MXML_CPPFLAGS="-I$withval/include";])

save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
save_LDFLAGS="$LDFLAGS"
LIBS="$LIBS -lmxml"
if test -n "$MXML_LDFLAGS"; then
    LDFLAGS="$LDFLAGS $MXML_LDFLAGS"
elif test -n "$MXML_LIB"; then
    LDFLAGS="$LDFLAGS $MXML_LIB"
    MXML_LDFLAGS="$MXML_LIB"
fi
if test -n "$MXML_CPPFLAGS"; then
    CPPFLAGS="$CPPFLAGS $MXML_CPPFLAGS"
elif test -n "$MXML_INC"; then
    CPPFLAGS="$CPPFLAGS $MXML_INC"
    MXML_CPPFLAGS="$MXML_INC"
fi 

if test -z "${HAVE_MXML_TRUE}"; then
        AC_CHECK_HEADERS(mxml.h,
                ,
                [AM_CONDITIONAL(HAVE_MXML,false)])
fi

# Check for the Mini-XML library and headers
AC_TRY_COMPILE([#include "mxml.h"],
        [mxml_node_t * n; mxmlWalkNext (n, n, MXML_DESCEND);],
        [MXML_LIBS="-lmxml"],
        [AM_CONDITIONAL(HAVE_MXML,false)])

LIBS="$save_LIBS"
LDFLAGS="$save_LDFLAGS"
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(MXML_LIBS)
AC_SUBST(MXML_LDFLAGS)
AC_SUBST(MXML_CPPFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test -z "${HAVE_MXML_TRUE}"; then
        ifelse([$1],,[AC_DEFINE(HAVE_MXML,1,[Define if you have the MXML.])],[$1])
        :
else
        $2
        :
fi
])dnl AC_MXML
