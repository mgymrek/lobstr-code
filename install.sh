# Script to install lobSTR binaries, for use with Arvados

set -e

make clean
make distclean
./configure
make
make check
mkdir -p $1
mv src/lobSTR $1/
mv src/allelotype $1/
mv src/lobSTRIndex $1/