dnl ######################################################################
dnl
dnl Finds PHDF5
dnl
dnl ######################################################################

AC_DEFUN([AC_PHDF5],
[
AM_CONDITIONAL(HAVE_PHDF5,true)

AC_ARG_WITH([phdf5],
            [  --with-phdf5=<location of PHDF5 installation>],
            [PHDF5_DIR=$withval], [with_phdf5=check])

dnl If --without-phdf5 was given set HAVE_PHDF5 to false and do nothing more
if test "x$with_phdf5" == "xno"; then

   AM_CONDITIONAL(HAVE_PHDF5,false)

else

    dnl allow args --with-phdf5 incdir and --with-phdf5-libdir

    AC_ARG_WITH(phdf5-incdir,
                [  --with-phdf5-incdir=<location of PHDF5 includes>],
                [PHDF5_INCDIR=$withval
                 with_phdf5=detailed])
    
    AC_ARG_WITH(phdf5-libdir,
                [  --with-phdf5-libdir=<location of PHDF5 library>],
                [PHDF5_LIBDIR=$withval
                 with_phdf5=detailed])
    
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
    elif test -n "${HDF5_CLIB}"; then
            PHDF5_CPPFLAGS="${HDF5_CLIB}"
            dnl    echo " --- HDF5_CLIB was defined. PHDF5_CPPFLAGS=${HDF5_CPPFLAGS}"
    else
            ac_phdf5_ok=no
    fi
    
    dnl Add "-L" to PHDF5_LIBDIR.
    if test -n "${PHDF5_LIBDIR}"; then
            PHDF5_LDFLAGS="-L${PHDF5_LIBDIR}"
    elif test -n "${HDF5_CLIB}"; then
            PHDF5_LDFLAGS="${HDF5_CLIB}"
            dnl    echo " --- HDF5_CLIB was defined. PHDF5_LDFLAGS=${HDF5_CPPFLAGS}"
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
        [if test "x$with_phdf5" != xcheck; then
           AC_MSG_FAILURE( [--with-phdf5 was given, but test for hdf5.h failed])
         fi
         AM_CONDITIONAL(HAVE_PHDF5,false)
        ])
    fi
    
    if test -z "${HAVE_PHDF5_TRUE}"; then
        AC_MSG_CHECKING([if phdf5 code can be compiled])
        AC_TRY_COMPILE([#include "hdf5.h"],
            [hid_t file_id;
             herr_t status;
#ifdef H5_HAVE_PARALLEL
             file_id = H5Fcreate("a.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
             status = H5Fclose(file_id);
#else
             /* This must deliberately fail */
             file_id = THE_HDF5_INSTALLATION_FOUND_IS_NOT_PARALLEL_HDF5 
#endif
            ],
            [AC_MSG_RESULT(yes)
	     PHDF5_LIBS="-lhdf5"],
            [AC_MSG_RESULT(no)
	     if test "x$with_phdf5" != xcheck; then
               AC_MSG_FAILURE( [--with-phdf5 was given, but compile test failed])
             fi
             AM_CONDITIONAL(HAVE_PHDF5,false)
            ])
        
        AC_SUBST(PHDF5_LIBS)
        AC_SUBST(PHDF5_LDFLAGS)
        AC_SUBST(PHDF5_CPPFLAGS)
    fi
    
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CPPFLAGS="$save_CPPFLAGS"
    
    # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
    if test -z "${HAVE_PHDF5_TRUE}"; then
            ifelse([$1],,[AC_DEFINE(HAVE_PHDF5,1,[Define if you have PHDF5.])],[$1])
            :
    else
            $2
            :
    fi
    
fi

])
