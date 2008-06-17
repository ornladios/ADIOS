aclocal -I m4
libtoolize --force
autoconf
autoheader
automake -a
