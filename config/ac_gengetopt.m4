#
#
# AC_PROG_GENGETOPT
#
# Test for gengetopt
# and set $GENGETOPT to the correct value.
#
#
dnl @synopsis AC_PROG_GENGETOPT
dnl
dnl This macro test if gengetopt is installed. If gengetopt
dnl is installed, it set $GENGETOPT to the right value
dnl
dnl @version 1.3
dnl @author Ron Oldfield raoldfi@sandia.gov
dnl
AC_DEFUN([AC_PROG_GENGETOPT],[
AC_ARG_WITH(gengetopt, AS_HELP_STRING([--with-gengetopt=<path>],[Location of gengetopt]),
        [GENGETOPT_PATH=$withval])
AC_PATH_PROGS(GENGETOPT,gengetopt,no,[$GENGETOPT_PATH:$PATH])
export GENGETOPT;
if test $GENGETOPT = "no" ;
then
        ac_gengetopt_ok=no;
        AC_MSG_WARN([Unable to find a gengetopt application]);
        AM_CONDITIONAL(HAVE_GENGETOPT, false)
else
        ac_gengetopt_ok=yes;
        AM_CONDITIONAL(HAVE_GENGETOPT, true)
fi;

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ac_gengetopt_ok" = xyes; then
        ifelse([$1],,[AC_DEFINE(HAVE_GENGETOPT,1,[Define if you have GENGETOPT.])],[$1])
        :
else
        $2
        :
fi

])


