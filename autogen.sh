aclocal -I config
libtoolize --force --copy
autoconf
autoheader
automake -a
