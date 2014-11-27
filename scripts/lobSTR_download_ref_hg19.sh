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
# A script to download and extract the lobSTR hg19 file
#

URL=http://files.teamerlich.org/lobstr/v3/ref/lobSTR_v3.0.2_hg19_resource_bundle.tar.gz
HOMEPAGE=http://lobstr.teamerlich.org/

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@">&2
    exit 1
}

usage()
{
    BASE=$(basename -- "$0")
    echo "lobSTR Reference Download

Usage:
   $BASE DEST-DIRECTORY

The script will create 'DEST-DIRECTORY', download the lobSTR
hg19 reference archive and extract it.

NOTE:
The lobSTR hg19 reference file is large (>2GB) - download might take a long
time, depending on your network connection.

For more information visit $HOMEPAGE
"
    exit 0
}

check_prog()
{
    which "$1" >/dev/null 2>/dev/null
}


DESTDIR="$1"
test -z "$DESTDIR" \
    && die "missing destination directory parameter. See --help for more information."
test "x$DESTDIR" = "x-h" || test "x$DESTDIR" = "x--help" \
    && usage

## Check if directory already exists
test -d "$DESTDIR" \
    && die "destination directory '$DESTDIR' already exists. aborting."

## Create the directory
mkdir -p -- "$DESTDIR" \
    || die "failed to create destination directory '$DESTDIR'"
cd "$DESTDIR" \
    || die "failed to change-directory to '$DESTDIR'"
DESTDIR=$PWD

## Download the archive (warning - can get big/slow)
BASEREF=$(basename -- "$URL")
if check_prog wget ; then
    wget -c -O "$BASEREF" -- "$URL" \
        || die "failed to download reference ($URL) using wget"
elif check_prog curl ; then
    curl --output "$BASEREF" --retry 10 -- "$URL" \
        || die "failed to download reference ($URL) using cURL"
else
    die "wget/curl not found - can not download file."
fi

## Extract files from active
tar -xzvf "$BASEREF" \
    || die "failed to extract archive '$BASEREF'"

## Remove archive to save space
rm -r "$BASEREF" \
    || die "failed to delete '$BASEREF'"

##
## Show the user how to use the newly downloaded reference
##

echo "

**
** lobSTR reference download completed
**

Reference stored in $DESTDIR/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref

To use the reference with lobSTR:

   lobSTR --index-prefix '$DESTDIR/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_' [...]

For more information visit $HOMEPAGE

"
