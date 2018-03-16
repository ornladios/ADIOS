
AC_DEFUN([AC_FLEXPATH], [


AC_MSG_NOTICE([=== checking for FLEXPATH ===])

FP_SRCDIR=""
FP_CPPFLAGS=""
FP_LDFLAGS=""
FP_LIBS=""

ac_flexpath_ok=yes

temptest=enable

flexpath_dir=""


AC_ARG_WITH(flexpath, 
	[  --with-flexpath=DIR 	Location of FlexPath], 
	[ ac_with_flexpath=$withval], [with_flexpath=no])

if test "x$with_flexpath" = "xno"; then
	ac_flexpath_ok=no
	temptest=disable

elif test x"$with_flexpath" != xyes -a x"$with_flexpath" != xcheck; then
		
AC_MSG_CHECKING(got with flexpath argument $with_flexpath)
#	with_evpath=$with_flexpath

fi

if test "x$ac_flexpath_ok" != "xno"; then

    flexpath_dir=$withval
    datatap_dir=$withval

    CERCS_REQUIRE_PACKAGE(evpath, evpath.h, libevpath.la)

    if test -z "$cercs_cv_evpath_link_dir" ;then
        unset cercs_cv_evpath_link_dir
        CERCS_REQUIRE_PACKAGE(evpath, evpath.h, libevpath.a)
    fi

    if test -n "$cercs_cv_evpath_link_dir" -a -n "$cercs_cv_evpath_include_arg";then
	if (test -x "$withval/bin/evpath_config") ; then
	  FP_CONFIG="$withval/bin/evpath_config"
	else
	  FP_CONFIG="$cercs_cv_evpath_link_dir/../bin/evpath_config"
	fi
	if (test -x "$FP_CONFIG") ; then
	  FP_LIBS=$(${FP_CONFIG} -s)
	  FP_CFLAGS=$(${FP_CONFIG} -c)
	  FP_CPPFLAGS=$(${FP_CONFIG} -c)
	else
	  FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_evpath_link_dir"
	  FP_LIBS="$FP_LIBS -levpath"
	  FP_CFLAGS="$FP_CFLAGS $cercs_cv_evpath_include_arg"
	  FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_evpath_include_arg"
	fi
    else
	echo "FLEXPATH couldn't find evpath -  Not building flexpath"
	ac_flexpath_ok=no
    fi
fi

AC_SUBST(FP_LIBS)
AC_SUBST(FP_CPPFLAGS)
AC_SUBST(FP_LDFLAGS)


if test "x$ac_flexpath_ok" = "xyes"; then
   HAVE_FLEXPATH=1
   NO_FLEXPATH=0
else
   HAVE_FLEXPATH=0
   NO_FLEXPATH=1
fi
AM_CONDITIONAL(HAVE_FLEXPATH, test "x$ac_flexpath_ok" = "xyes")
AC_DEFINE_UNQUOTED(NO_FLEXPATH, $NO_FLEXPATH, [Flexpath is disabled])
AC_DEFINE_UNQUOTED(HAVE_FLEXPATH, $HAVE_FLEXPATH, [Flexpath is enabled])

]) dnl AC_FLEXPATH
