aclocal -I config
libtoolize --force --copy
autoconf
autoheader
automake --add-missing --copy
