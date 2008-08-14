dnl ######################################################################
dnl
dnl Finds netCDF
dnl
dnl ######################################################################

AC_DEFUN([AC_NETCDF],
[
AM_CONDITIONAL(HAVE_NETCDF,true)

AC_ARG_WITH(netcdf,
            [  --with-netcdf=<location of NetCDF installation>],
            [NETCDF_DIR=$withval])

AC_ARG_WITH(netcdf-incdir,
            [  --with-netcdf-incdir=<location of NetCDF includes>],
            [NETCDF_INCDIR=$withval])

AC_ARG_WITH(netcdf-libdir,
            [  --with-netcdf-libdir=<location of NetCDF library>],
            [NETCDF_LIBDIR=$withval])

dnl If we know NETCDF_DIR, then we can know NETCDF_INCDIR.
dnl We don't overwrite NETCDF_INCDIR.
if test -n "${NETCDF_DIR}" -a -z "${NETCDF_INCDIR}"; then
        NETCDF_INCDIR="${NETCDF_DIR}/include";
else
	ac_netcdf_ok=no
fi

dnl If we know NETCDF_DIR, then we can know NETCDF_LIBDIR.
dnl We don't overwrite NETCDF_LIBDIR.
if test -n "${NETCDF_DIR}" -a -z "${NETCDF_LIBDIR}"; then
        NETCDF_LIBDIR="${NETCDF_DIR}/lib";
else
	ac_netcdf_ok=no
fi

dnl Add "-I" to NETCDF_INCDIR.
if test -n "${NETCDF_INCDIR}"; then
	NETCDF_CPPFLAGS="-I${NETCDF_INCDIR}"
else
	ac_netcdf_ok=no
fi

dnl Add "-L" to NETCDF_LIBDIR.
if test -n "${NETCDF_LIBDIR}"; then
	NETCDF_LDFLAGS="-L${NETCDF_LIBDIR}"
else
	ac_netcdf_ok=no
fi

save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
save_LDFLAGS="$LDFLAGS"
LIBS="$LIBS -lnetcdf"
LDFLAGS="$LDFLAGS $NETCDF_LDFLAGS"
CPPFLAGS="$CPPFLAGS $NETCDF_CPPFLAGS"

if test -z "${HAVE_NETCDF_TRUE}"; then
        AC_CHECK_HEADERS(netcdf.h,
                ,
                [AM_CONDITIONAL(HAVE_NETCDF,false)])
fi

AC_TRY_COMPILE([#include "netcdf.h"],
        [int ncid;
         nc_create("a.nc", NC_CLOBBER, &ncid);
         nc_close(ncid);
        ],
        [NETCDF_LIBS="-lnetcdf"],
        [AM_CONDITIONAL(HAVE_NETCDF,false)])

LIBS="$save_LIBS"
LDFLAGS="$save_LDFLAGS"
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(NETCDF_LIBS)
AC_SUBST(NETCDF_LDFLAGS)
AC_SUBST(NETCDF_CPPFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test -z "${HAVE_NETCDF_TRUE}"; then
        ifelse([$1],,[AC_DEFINE(HAVE_NETCDF,1,[Define if you have NETCDF.])],[$1])
        :
else
        $2
        :
fi

])
