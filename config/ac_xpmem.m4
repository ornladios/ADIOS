dnl @synopsis AC_XPMEM([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro tries to find out how to compile programs that
dnl use the Xpmem API.
dnl
dnl On success, it defines HAVE_XPMEM and sets XPMEM_LIBS
dnl to any libraries that are needed for linking
dnl Xpmem 
dnl
dnl If you want to compile everything with Xpmem, you should set:
dnl
dnl     LIBS="$XPMEM_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a Portals
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_PORTALS.
dnl
dnl @version $Id: acx_mpi.m4 676 2006-05-16 20:44:08Z raoldfi $
dnl @author Ron A. Oldfield <raoldfi@sandia.gov>

AC_DEFUN([AC_XPMEM], [

AM_CONDITIONAL(HAVE_XPMEM,true)

AC_ARG_WITH(xpmem,
	[  --with-xpmem=DIR      Location of xpmem library],
	[XPMEM_LDFLAGS="-L$withval/lib";
	 XPMEM_CPPFLAGS="-I$withval/include";])

dnl AC_LANG_PUSH([C++])

save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
save_LDFLAGS="$LDFLAGS"
LIBS="$LIBS -lxpmem"
LDFLAGS="$LDFLAGS $XPMEM_LDFLAGS"
CPPFLAGS="$CPPFLAGS $XPMEM_CPPFLAGS"

if test -z "${HAVE_XPMEM_TRUE}"; then
        AC_CHECK_HEADERS(xpmem.h,
		,
		[AM_CONDITIONAL(HAVE_XPMEM,false)])
fi

# Check for the Xpmem library and headers for SDDF
AC_TRY_COMPILE([#include <Attributes.h>],
	[Attributes attr;attr.clearEntries();],
	[XPMEM_LIBS="-lXpmem"],
	[AM_CONDITIONAL(HAVE_XPMEM,false)])

LIBS="$save_LIBS"
LDFLAGS="$save_LDFLAGS"
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(XPMEM_LIBS)
AC_SUBST(XPMEM_LDFLAGS)
AC_SUBST(XPMEM_CPPFLAGS)

dnl AC_LANG_POP([C++])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test -z "${HAVE_XPMEM_TRUE}"; then
        ifelse([$1],,[AC_DEFINE(HAVE_XPMEM,1,[Define if you have the Xpmem.])],[$1])
        :
else
        $2
        :
fi
])dnl AC_XPMEM
