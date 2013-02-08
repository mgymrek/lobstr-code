#!/bin/bash

# simple script to pack up lobstr source

echo "packing up $1.src.tar.gz"
mkdir $1
make check
make clean
make distclean
cp config* $1
cp -r config* $1
cp -r m4 $1
cp install* $1
cp Make* $1
cp ax* $1
cp ac* $1
cp reconf $1
cp -r src $1
cp -r data $1
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

tar cvfz $1.src.tar.gz $1
