#!/usr/bin/env sh
aclocal -I config
case `uname` in Darwin*) glibtoolize --copy ;;
  *) libtoolize --copy ;; esac
autoconf
autoheader
automake --add-missing --copy
