dnl @synopsis AC_RDMACM([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro tries to find out how to compile programs that
dnl use the rdmacm rmda_cma.h API (together with ibverbs)
dnl
dnl On success, it defines HAVE_RDMACM and sets RDMACM_LIBS
dnl to any libraries that are needed for linking
dnl rdmacm (e.g. -lrdmacm,...).
dnl
dnl It requires INFINIBAND.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a rdmacm
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_RDMACM.
dnl
dnl @version $Id: ac_rdmacm.m4 676 2006-05-16 20:44:08Z  pnorbert $
dnl @author Norbert Podhorszki <pnorbert@ornl.gov>

AC_DEFUN([AC_RDMACM], [
AC_REQUIRE([AC_INFINIBAND])
AC_LANG_SAVE
AC_LANG_C
ac_rdmacm_hdr_ok=no
ac_rdmacm_lib_ok=no
ac_with_rdmacm=no


RDMACM_CFLAGS=""
RDMACM_CPPFLAGS=""
RDMACM_LDFLAGS=""
RDMACM_LIBS=""

AC_ARG_WITH(rdmacm,
        [  --with-rdmacm=DIR      Location of rdmacm],
        [ ac_with_rdmacm=yes;])



if test x"$withval" = xno; then

        ac_with_rdmacm=no;

elif test x"$withval" = xyes -o x"$withval" = x; then

        RDMACM_CPPFLAGS="";
        RDMACM_LDFLAGS="";
        ac_with_rdmacm=yes;

else

        RDMACM_CPPFLAGS="-I$withval/include";
        RDMACM_LDFLAGS="-L$withval/lib64 -L$withval/lib";
        ac_with_rdmacm=yes;

fi

AM_CONDITIONAL(HAVE_RDMACM,test x$ac_with_rdmacm = xyes)


dnl Check for command-line disable
if test x"$ac_with_rdmacm" = xyes; then

        dnl Look for rdmacm header files
        save_CPPFLAGS=$CPPFLAGS;
        CPPFLAGS="$CPPFLAGS $RDMACM_CPPFLAGS"
        LDFLAGS="$LDFLAGS $RDMACM_LDFLAGS"


        if test x"$ac_rdmacm_hdr_ok" = xno; then
                AC_CHECK_HEADER(rdma/rdma_cma.h,
                         [AC_DEFINE(HAVE_RDMA_CMA_H, 1,
                                [Define to 1 if you have <rdma/rdma_cma.h>.])
                         ac_rdmacm_hdr_ok=yes;
                         RDMACM_CFLAGS="$RDMACM_CFLAGS $EXTRA_CFLAGS";
                         RDMACM_CPPFLAGS="$RDMACM_CPPFLAGS $EXTRA_CFLAGS"],
                         [ac_rdmacm_hdr_ok=no])
        fi

        if test x"$ac_rdmacm_hdr_ok" = xno; then
                CPPFLAGS=$save_CPPFLAGS
        fi

        dnl Look for -lrdmacm
        if test x"$ac_rdmacm_lib_ok" = xno -a x$ac_rdmacm_hdr_ok = xyes; then
            save_LIBS=$LIBS;
            LIBS=""
            AC_SEARCH_LIBS(rdma_connect,[rdmacm],
                    [ac_rdmacm_lib_ok=yes],
                    [ac_rdmacm_lib_ok=no],
                    [$save_LIBS $INFINIBAND_LIBS])
            if test -n $LIBS; then
                RDMACM_LIBS="$LIBS";
            fi
            LIBS=$save_LIBS;
        fi

        if test x"$ac_rdmacm_hdr_ok" = xno -o x$ac_rdmacm_lib_ok = xno; then
            AM_CONDITIONAL(HAVE_RDMACM,false)
            RDMACM_CFLAGS=""
            RDMACM_CPPFLAGS=""
            RDMACM_LDFLAGS=""
            RDMACM_LIBS=""
            ac_with_rdmacm=no
        else
            AM_CONDITIONAL(HAVE_RDMACM,true)
        fi


        AC_SUBST(RDMACM_CPPFLAGS)
        AC_SUBST(RDMACM_CFLAGS)
        AC_SUBST(RDMACM_LDFLAGS)
        AC_SUBST(RDMACM_LIBS)
        AC_SUBST(RDMACM_HEADER)


        # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
        if test x"$ac_rdmacm_lib_ok" = xyes; then
                ifelse([$1],,[AC_DEFINE(HAVE_RDMACM,1,[Define if you have the rdmacm.])],[$1])
                :
        else
                $2
                :
        fi

fi

AC_LANG_RESTORE
])dnl AC_RDMACM
