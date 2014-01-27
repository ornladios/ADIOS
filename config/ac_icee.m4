
AC_DEFUN([AC_ICEE], [


AC_MSG_NOTICE([=== checking for ICEE ===])

FP_SRCDIR=""
FP_CPPFLAGS=""
FP_LDFLAGS=""
FP_LIBS=""

ac_icee_ok=yes

temptest=enable

icee_dir=""


AC_ARG_WITH(icee, 
	[  --with-icee=DIR 	Location of icee], 
	[ ac_with_icee=$withval], [with_icee=no])

if test "x$with_icee" = "xno"; then
	ac_icee_ok=no
	temptest=disable

elif test x"$with_icee" != xyes -a x"$with_icee" != xcheck; then
		
AC_MSG_CHECKING(got with icee argument $with_icee)
#	with_evpath=$with_icee

fi

if test "x$ac_icee_ok" != "xno"; then

    icee_dir=$withval
    datatap_dir=$withval

    CERCS_REQUIRE_PACKAGE(evpath, evpath.h, libevpath.a)
    CERCS_REQUIRE_PACKAGE(ffs, ffs.h,libffs.a)
    CERCS_REQUIRE_PACKAGE(atl, atl.h,libatl.a)
    CERCS_REQUIRE_PACKAGE(dill, dill.h, libdill.a)
    CERCS_REQUIRE_PACKAGE(cercs_env, cercs_env.h, libcercs_env.a)


    if test -n "$cercs_cv_evpath_link_dir" -a -n "$cercs_cv_evpath_include_arg";then
	FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_evpath_link_dir"
	FP_LIBS="$FP_LIBS -levpath"
	FP_CFLAGS="$FP_CFLAGS $cercs_cv_evpath_include_arg"
	FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_evpath_include_arg"
    else
	echo "ICEE couldn't find evpath -  Not building icee"
	ac_icee_ok=no
    fi
    if test -n "$cercs_cv_ffs_link_dir" -a -n "$cercs_cv_ffs_include_arg"; then
	FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_ffs_link_dir"
	FP_LIBS="$FP_LIBS -lffs"
	FP_CFLAGS="$FP_CFLAGS $cercs_cv_ffs_include_arg"
	FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_ffs_include_arg"
    else 
	echo "ICEE couldn't find ffs -  Not building icee"
	ac_icee_ok=no
    fi
    if test -n "$cercs_cv_atl_link_dir" -a -n "$cercs_cv_atl_include_arg"; then
	FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_atl_link_dir"
	FP_LIBS="$FP_LIBS -latl"
	FP_CFLAGS="$FP_CFLAGS $cercs_cv_atl_include_arg"
	FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_atl_include_arg"
    else 
	ac_icee_ok=no
	echo "ICEE couldn't find atl -  Not building icee"
    fi
    if test -n "$cercs_cv_dill_link_dir" -a -n "$cercs_cv_dill_include_arg"; then
	FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_dill_link_dir"
	FP_LIBS="$FP_LIBS -ldill"
	FP_CFLAGS="$FP_CFLAGS $cercs_cv_dill_include_arg"
	FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_dill_include_arg"
    else 
	ac_icee_ok=no
	echo "ICEE couldn't find dill -  Not building icee"
    fi
    if test -n "$cercs_cv_cercs_env_link_dir" -a -n "$cercs_cv_cercs_env_include_arg"; then
	FP_LDFLAGS="$FP_LDFLAGS -L$cercs_cv_cercs_env_link_dir"
	FP_LIBS="$FP_LIBS -lcercs_env"
	FP_CFLAGS="$FP_CFLAGS $cercs_cv_cercs_env_include_arg"
	FP_CPPFLAGS="$FP_CPPFLAGS $cercs_cv_cercs_env_include_arg"
    else 
	ac_icee_ok=no
	echo "ICEE couldn't find cercs_env -  Not building icee"
    fi

fi

AC_SUBST(FP_LIBS)
AC_SUBST(FP_CPPFLAGS)
AC_SUBST(FP_LDFLAGS)


if test "x$ac_icee_ok" = "xyes"; then
   HAVE_ICEE=1
   NO_ICEE=0
else
   HAVE_ICEE=0
   NO_ICEE=1
fi
AM_CONDITIONAL(HAVE_ICEE, test "x$ac_icee_ok" = "xyes")
AC_DEFINE_UNQUOTED(NO_ICEE, $NO_ICEE, [icee is disabled])
AC_DEFINE_UNQUOTED(HAVE_ICEE, $HAVE_ICEE, [icee is enabled])

]) dnl AC_ICEE
