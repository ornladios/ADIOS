AC_DEFUN([AX_NCSU_LIBTIMER], [

dnl Enable the --with-timer=path configure argument
AC_ARG_WITH(
  [timer],
  [AS_HELP_STRING(
    [--with-timer=DIR],
    [Location of the timer library]
  ),[],
  [with_timer=no]]dnl
)

dnl If the timer lib was specified, verify that it exists and can compile
if test "x$with_timer" != xno -a "x$with_timer" != x; then
    TIMER_CPPFLAGS="-I$with_timer/include"
    TIMER_LDFLAGS="-L$with_timer/lib"
    TIMER_LIBS="-ltimer"

    saveLIB="$LIB"
    saveLDFLAGS="$LDFLAGS"
    saveCPPFLAGS="$CPPFLAGS"
    LIB="$LIB $TIMER_LIBS"
    LDFLAGS="$LDFLAGS $TIMER_LDFLAGS"
    CPPFLAGS="$CPPFLAGS $TIMER_CPPFLAGS"

    AC_CHECK_HEADERS(
      [timer.h],
      [],
      [AC_MSG_FAILURE(
        [Cannot find timer.h from the timer lib. Make sure it has been properly installed at the path specified ($with_timer).]dnl
      )]dnl
    )

    AC_CHECK_LIB(
      [timer],
      [timer_init],
      [AC_DEFINE(
        [HAVE_LIBTIMER],
        [1],
        [Define if you have libtimer]
      )],
      [AC_MSG_FAILURE(
        [Cannot successfully link with the timer lib. Make sure it has been properly installed at the path specified ($with_timer).]dnl
      )]dnl
    )

    LIBS="$saveLIBS"
    LDFLAGS="$saveLDFLAGS"
    CPPFLAGS="$saveCPPFLAGS"

    AC_SUBST(TIMER_CPPFLAGS)
    AC_SUBST(TIMER_LDFLAGS)
    AC_SUBST(TIMER_LIBS)
fi

]) dnl End of DEFUN
