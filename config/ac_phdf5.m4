dnl ######################################################################
dnl
dnl Finds PHDF5
dnl
dnl ######################################################################

AC_DEFUN([AC_PHDF5],
[
AM_CONDITIONAL(HAVE_PHDF5,true)

AC_ARG_WITH(phdf5,
            [  --with-phdf5=<location of PHDF5 installation>],
            [PHDF5_DIR=$withval])

AC_ARG_WITH(phdf5-incdir,
            [  --with-phdf5-incdir=<location of PHDF5 includes>],
            [PHDF5_INCDIR=$withval])

AC_ARG_WITH(phdf5-libdir,
            [  --with-phdf5-libdir=<location of PHDF5 library>],
            [PHDF5_LIBDIR=$withval])

dnl If we know PHDF5_DIR, then we can know PHDF5_INCDIR.
dnl We don't overwrite PHDF5_INCDIR.
if test -n "${PHDF5_DIR}" -a -z "${PHDF5_INCDIR}"; then
        PHDF5_INCDIR="${PHDF5_DIR}/include";
else
	ac_phdf5_ok=no
fi

dnl If we know PHDF5_DIR, then we can know PHDF5_LIBDIR.
dnl We don't overwrite PHDF5_LIBDIR.
if test -n "${PHDF5_DIR}" -a -z "${PHDF5_LIBDIR}"; then
        PHDF5_LIBDIR="${PHDF5_DIR}/lib";
else
	ac_phdf5_ok=no
fi

dnl Add "-I" to PHDF5_INCDIR.
if test -n "${PHDF5_INCDIR}"; then
	PHDF5_CPPFLAGS="-I${PHDF5_INCDIR}"
else
	ac_phdf5_ok=no
fi

dnl Add "-L" to PHDF5_LIBDIR.
if test -n "${PHDF5_LIBDIR}"; then
	PHDF5_LDFLAGS="-L${PHDF5_LIBDIR}"
else
	ac_phdf5_ok=no
fi

save_CPPFLAGS="$CPPFLAGS"
save_LIBS="$LIBS"
save_LDFLAGS="$LDFLAGS"
LIBS="$LIBS -lhdf5"
LDFLAGS="$LDFLAGS $PHDF5_LDFLAGS"
CPPFLAGS="$CPPFLAGS $PHDF5_CPPFLAGS"

if test -z "${HAVE_PHDF5_TRUE}"; then
        AC_CHECK_HEADERS(hdf5.h,
                ,
                [AM_CONDITIONAL(HAVE_PHDF5,false)])
fi

AC_TRY_COMPILE([#include "hdf5.h"],
        [hid_t file_id;
         herr_t status;
         file_id = H5Fcreate("a.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
         status = H5Fclose(file_id);
        ],
        [PHDF5_LIBS="-lhdf5"],
        [AM_CONDITIONAL(HAVE_PHDF5,false)])

LIBS="$save_LIBS"
LDFLAGS="$save_LDFLAGS"
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(PHDF5_LIBS)
AC_SUBST(PHDF5_LDFLAGS)
AC_SUBST(PHDF5_CPPFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test -z "${HAVE_PHDF5_TRUE}"; then
        ifelse([$1],,[AC_DEFINE(HAVE_PHDF5,1,[Define if you have PHDF5.])],[$1])
        :
else
        $2
        :
fi

])
