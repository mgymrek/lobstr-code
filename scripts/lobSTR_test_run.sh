#!/bin/sh

# Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>
#
# This file is part of lobSTR.
#
# lobSTR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lobSTR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

#
# A quick test ro run lobstr with the test reference,
# Ensuring everything is installed properly.
#

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@">&2
    exit 1
}

##
## These directories are defined in Makefile.am,
## and created during "make install".
##
DIR=$(dirname -- "$0")
BINDIR="$DIR/../../../bin/"
## These directories are defined in tests/Makefile.am:
REFDIR="$DIR/../test-ref"
SAMPLEDIR="$DIR/../sample"

## Ensure directories exist
test -d "$BINDIR" \
    || die "installation error: bin-dir '$BINDIR' not found"
test -d "$REFDIR" \
    || die "installation error: ref-dir '$REFDIR' not found"
test -d "$SAMPLEDIR" \
    || die "installation error: sample-dir '$SAMPLEDIR' not found"

##
## Ensure expected files exists
##
LOBSTRBIN="$BINDIR/lobSTR"
REFPREFIX="$REFDIR/lobSTR_"
REFFILE="${REFPREFIX}chromsizes.tab"
SAMPLEFQ="$SAMPLEDIR/tiny.fq"

test -e "$LOBSTRBIN" \
    || die "installation error: lobSTR binary ($LOBSTRBIN) not found"
test -x "$LOBSTRBIN" \
    || die "installation error: lobSTR binary ($LOBSTRBIN) not executable"
test -e "$SAMPLEFQ" \
    || die "installation error: lobSTR sample FASTQ file ($SAMPLEFQ) not found"
test -e "$REFFILE" \
    || die "installation error: lobSTR reference file ($REFFILE) not found"

##
## All seems well, run lobSTR
##
TMP=$(mktemp -d) || die "failed to create temporary directory"
trap "rm -r '$TMP'" EXIT

"$LOBSTRBIN" \
    --verbose \
    --index-prefix "$REFPREFIX" \
    -f "$SAMPLEFQ" --fastq \
    --out "$TMP/test" \
    --rg-sample "test" \
    --rg-lib "test" \
    || die "lobSTR failed (check error messages above)"

