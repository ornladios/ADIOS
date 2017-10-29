#
#
# AC_ZMQ
#
#
#
dnl @synopsis AC_ZMQ
dnl
dnl This macro test if zmq (for Titan only) is to be used. 
dnl Use in C code:
dnl     #ifdef ZMQ
dnl     #include "zmq.h"
dnl     #endif
dnl
dnl @version 1.0
dnl @author Qing Liu, ORNL
dnl
AC_DEFUN([AC_ZMQ],[

AC_MSG_NOTICE([=== checking for ZMQ ===])

AM_CONDITIONAL(HAVE_ZMQ,true)

AC_ARG_WITH(zmq,
        [  --with-zmq=DIR      Location of ZMQ library],
        [ZMQ_LDFLAGS="-L$withval/lib";
         ZMQ_LIBS="-lzmq";
         ZMQ_CPPFLAGS="-I$withval/include";],
        [with_zmq=no])

if test "x$with_zmq" == "xno"; then

   AM_CONDITIONAL(HAVE_ZMQ,false)

else

    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lzmq"
    LDFLAGS="$LDFLAGS $ZMQ_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $ZMQ_CPPFLAGS"
    
    # Check for the ZMQ library and headers
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    AC_SUBST(ZMQ_LIBS)
    AC_SUBST(ZMQ_LDFLAGS)
    AC_SUBST(ZMQ_CPPFLAGS)
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_ZMQ_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_ZMQ,1,[Define if you have ZMQ.])],[$1])
            :
    else
            $2
            :
    fi
fi
])dnl AC_ZMQ
