dnl ######################################################################
dnl
dnl Finds HDF5
dnl
dnl ######################################################################

AC_DEFUN([AC_HDF5],
[
AC_MSG_NOTICE([=== checking for HDF5 ===])

AM_CONDITIONAL(HAVE_HDF5,true)

AC_ARG_WITH(hdf5,
            [  --with-hdf5=<location of HDF5 installation>],
            [HDF5_DIR=$withval], [with_hdf5=check])

dnl If --without-hdf5 was given set HAVE_HDF5 to false and do nothing more
if test "x$with_hdf5" == "xno"; then

   AM_CONDITIONAL(HAVE_HDF5,false)

else

   dnl allow args --with-hdf5 incdir and --with-hdf5-libdir

   AC_ARG_WITH(hdf5-incdir,
                [  --with-hdf5-incdir=<location of HDF5 includes>],
                [HDF5_INCDIR=$withval
                 with_hdf5=detailed])
    
   AC_ARG_WITH(hdf5-libdir,
                [  --with-hdf5-libdir=<location of HDF5 library>],
                [HDF5_LIBDIR=$withval
                 with_hdf5=detailed])

   AC_ARG_WITH(hdf5-libs,
                [  --with-hdf5-libs=<linker flags besides -L<hdf5_libdir>, e.g. -lhdf5 -lhdf5_hl -lz>],
                [HDF5_LIBS=$withval
                 with_hdf5=detailed])
    
    
    ac_use_cray_hdf5=no  dnl will set to yes if we will use CRAY_HDF5_DIR below

    dnl If we know HDF5_DIR, then we can know HDF5_INCDIR.
    dnl If we know CRAY_HDF5_DIR, then we leave HDF5_INCDIR empty.
    dnl We don't overwrite HDF5_INCDIR.
    if test -z "${HDF5_INCDIR}"; then
        if test -n "${HDF5_DIR}"; then
            HDF5_INCDIR="${HDF5_DIR}/include";
        elif test -n "${CRAY_HDF5_DIR}"; then
            HDF5_INCDIR="";
            ac_use_cray_hdf5=yes
        fi
    fi
    
    dnl If we know HDF5_DIR, then we can know HDF5_LIBDIR.
    dnl If we know CRAY_HDF5_DIR, then we leave HDF5_LIBDIR empty.
    dnl We don't overwrite HDF5_LIBDIR.
    if test -z "${HDF5_LIBDIR}"; then
        if test -n "${HDF5_DIR}"; then
            HDF5_LIBDIR="${HDF5_DIR}/lib";
        elif test -n "${CRAY_HDF5_DIR}"; then
            HDF5_LIBDIR="";
            ac_use_cray_hdf5=yes
        fi
    fi
    
    if test -n "${HDF5_CLIB}"; then
        HDF5_CPPFLAGS="${HDF5_CLIB}"
        dnl echo " --- HDF5_CLIB was defined. HDF5_CPPFLAGS=${HDF5_CPPFLAGS}"
    elif test -n "${HDF5_INCDIR}"; then
        dnl Add "-I" to HDF5_INCDIR.
        HDF5_CPPFLAGS="-I${HDF5_INCDIR}"
    else
        ac_hdf5_ok=no
    fi

    if test -n "${HDF5_CLIB}"; then
        HDF5_LDFLAGS="${HDF5_CLIB}"
        dnl echo " --- HDF5_CLIB was defined. HDF5_LDFLAGS=${HDF5_CPPFLAGS}"
    elif test -n "${HDF5_LIBDIR}"; then
        dnl Add "-L" to HDF5_LIBDIR.
        HDF5_LDFLAGS="-L${HDF5_LIBDIR} ${HDF5_LIBS}"
    else
        ac_hdf5_ok=no
    fi
    
    dnl if hdf5 libs are not defined (and not Cray hdf5 lib), then guess and define it
    if test -z "${HDF5_LIBS}"; then
        if test "${ac_use_cray_hdf5}" != "yes"; then
            dnl default HDF5 lib is usually just -lhdf5
            HDF5_LIBS="-lhdf5"
        fi
    fi
    
    save_CC="$CC"
    save_CPPFLAGS="$CPPFLAGS"
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS $HDF5_LIBS"
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
    CC="$MPICC"
    
    if test -z "${HAVE_HDF5_TRUE}"; then
        AC_CHECK_HEADERS(hdf5.h,
            ,
            [if test "x$with_hdf5" != xcheck; then
                AC_MSG_FAILURE( [--with-hdf5 was given, but test for hdf5.h failed])
             fi
             AM_CONDITIONAL(HAVE_HDF5,false)])
    fi
    
    if test -z "${HAVE_HDF5_TRUE}"; then
        AC_MSG_CHECKING([if hdf5 code can be compiled])
        AC_TRY_COMPILE([#include "hdf5.h"],
            [hid_t file_id;
             herr_t status;
             file_id = H5Fcreate("a.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
             status = H5Fclose(file_id);
            ],
            [AC_MSG_RESULT(yes)],
            [AC_MSG_RESULT(no)
             if test "x$with_hdf5" != xcheck; then
                AC_MSG_FAILURE( [--with-hdf5 was given, but compile test failed])
             fi
             AM_CONDITIONAL(HAVE_HDF5,false)
            ])
    
        AC_SUBST(HDF5_LIBS)
        AC_SUBST(HDF5_LDFLAGS)
        AC_SUBST(HDF5_CPPFLAGS)
    fi
    
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    CC="$save_CC"
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_HDF5_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_HDF5,1,[Define if you have HDF5.])],[$1])
            :
    else
            $2
            :
    fi

fi

])
