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

if test -z "${HAVE_MXML_TRUE}"; then
    # Check for the Mini-XML library and headers
    AC_MSG_CHECKING([if mxml code can be linked])
    AC_TRY_LINK([#include "mxml.h"],
        [mxml_node_t * n; 
         char *buffer;
         char *value;
         n = mxmlLoadString (0, buffer, MXML_TEXT_CALLBACK);
         mxmlWalkNext (n, n, MXML_DESCEND);
         value = mxmlElementGetAttr (n, "value");
         mxmlRelease (n);
        ],
        [MXML_LIBS="-lmxml"
         AC_MSG_RESULT(yes)
        ],
        [AM_CONDITIONAL(HAVE_MXML,false)
         AC_MSG_RESULT(no)
        ])

    dnl If Linking above failed, one reason might be that mxml uses pthreads and
    dnl the compiler does not use it by default. Try getting phtreads
    if test -z "${HAVE_MXML_FALSE}"; then
        # Check for the Mini-XML library and headers
        AC_REQUIRE([ACX_PTHREAD])
        LDFLAGS="$LDFLAGS $PTHREAD_LDFLAGS $PTHREAD_LIBS"
        AC_MSG_CHECKING([if mxml code can be linked using pthreads])
        AC_TRY_LINK([#include "mxml.h"],
            [mxml_node_t * n; 
             char *buffer;
             char *value;
             n = mxmlLoadString (0, buffer, MXML_TEXT_CALLBACK);
             mxmlWalkNext (n, n, MXML_DESCEND);
             value = mxmlElementGetAttr (n, "value");
             mxmlRelease (n);
            ],
            [MXML_LDFLAGS="$MXML_LDFLAGS $PTHREAD_LDFLAGS"
             MXML_LIBS="-lmxml $PTHREAD_LIBS"
             AM_CONDITIONAL(HAVE_MXML,true)
             AC_MSG_RESULT(yes)
            ],
            [AM_CONDITIONAL(HAVE_MXML,false)
             AC_MSG_RESULT(no)
            ])
    fi
fi


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
