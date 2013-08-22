#!/bin/bash

# simple script to pack up a precompiled binary package, with the boost thread
# library statically linked in.
# Script provided graciously by Cole Trapnell


echo "packing up $1.tar.gz"
mkdir $1
make clean
make distclean
./configure
make
make check
mkdir $1/bin
mv src/lobSTR $1/bin/
mv src/allelotype $1/bin/
mv src/lobSTRIndex $1/bin/
cp config* $1
cp -r config* $1
cp -r m4 $1
cp install* $1
cp Make* $1
cp ax* $1
cp ac* $1
cp reconf $1
cp -r src $1
mkdir $1/scripts
cp scripts/lobstr_index.py $1/scripts/
cp scripts/*check*.py $1/scripts/
cp scripts/GetSTRInfo.py $1/scripts/GetSTRInfo.py
cp -r tests $1
cp -r models $1
cp README $1
cp AUTHORS $1
cp INSTALL $1
cp NEWS $1
cp ChangeLog $1
cp COPYING $1

tar cvfz $1.tar.gz $1
