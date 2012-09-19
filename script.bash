#!/bin/bash
echo "- libtoolize..."
libtoolize
echo "- aclocal..."
aclocal -I m4

echo "- autoconf..."
autoconf

echo "- autoheader..."
##autoheader

echo "- automake..."
#automake -a -v
automake -a -c

echo "- ready!"
exit
