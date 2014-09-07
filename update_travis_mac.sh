#!/bin/sh
##
## Travis-CI for MAC-OS-X auto-update script.
##
## Copyright (C) 2014 Assaf Gordon <assafgordon@gmail.com>
##
## This script is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This script is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this script. If not, see <http://www.gnu.org/licenses/>.
##


##
## This script will:
##  1. branch from the current active branch,
##     to a new (fixed nam) Travis-CI-Mac-OS branch (re-creating it if needed).
##  2. Create Travis-CI configuration for MAC-OS-X
##  3. Push the new branch to GitHub (optionally deleting previous branch)
##  4. Return to previously active branch
##

## New of branch to use (locally and remotely)
TRAVIS_MAC_BRANCH=travis_mac_osx

## New Travis-CI configuration, adapted for Mac-OS
NEW_TRAVIS_YML=\
"language: objective-c
compiler:
  - clang
before_install:
  - brew update
  - brew install lcov
  - brew install gettext
  - brew install fftw
  - brew install boost
  - brew install cppunit
  - brew install gsl
script:
- \"( git fetch --unshallow ; true )\"
- ls -l
- ./reconf
- ls -l
- git describe --always
- ./configure
- make
- make check VERBOSE=yes
after_script:
- cat ./tests/test-suite.log ; true"




die()
{
BASE=$(basename "$0")
echo "$BASE: error: $@" >&2
exit 1
}

##
## Stop if there are uncommited changes.
##
UNCOMMITTED=$(git status --porcelain -uno) \
	|| die "failed to check for uncommited changes"
test -z "$UNCOMMITTED" \
	|| die "there are uncommited changes in the current directory. Aborting."

## Get current branch name
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD) \
	|| die "failed to find current branch name"
test -z "$CURRENT_BRANCH" \
	&& die "got empty current branch name"
test "$CURRENT_BRANCH" = "$TRAVIS_MAC_BRANCH" \
	&& die "current branch is same as travis-ci branch"


##
## Delete local Travis-Mac branch (if any)
##
LOCAL=$(git branch --list "$TRAVIS_MAC_BRANCH") \
	|| die "failed to check for local branch '$TRAVIS_MAC_BRANCH'"
if test -n "$LOCAL" ; then
	# There's a local branch, delete it
	git branch -q -D "$TRAVIS_MAC_BRANCH" \
		|| die "failed to delete local branch '$TRAVIS_MAC_BRANCH'"
fi

##
## Delete remote Travis-Mac branch (if any)
##
REMOTE=$(git branch -r | grep "origin/$TRAVIS_MAC_BRANCH\$")
if test -n "$REMOTE" ; then
	# Delete remote branch
	git branch -q -r -D "origin/$TRAVIS_MAC_BRANCH" \
		|| die "failed to delete remote branch 'origin/$TRAVIS_MAC_BRANCH'"

	# Push (delete) remote branch on temote server (e.g. github)
	git push -q origin ":$TRAVIS_MAC_BRANCH"
fi

##
## Create a new branch
##
git checkout -q -b "$TRAVIS_MAC_BRANCH" "$CURRENT_BRANCH" \
	|| die "failed to checkout new branch"

##
## Write a new Travis-CI configuration file
##

echo "$NEW_TRAVIS_YML" > .travis.yml \
	|| die "failed to update .travis.yml"

git add .travis.yml \
	|| die "failed to git-add .travis.yml"

git commit -q -m "Auto-Updated Travis-CI for Mac configuration" \
	|| die "failed to git-commit .travis.yml"

##
## Push new branch to remote server (e.g. github)
##
git push -q origin $TRAVIS_MAC_BRANCH:$TRAVIS_MAC_BRANCH \
	|| die "failed to push branch '$TRAVIS_MAC_BRANCH' to remote server"

##
## Go back to master branch
##
git checkout -q "$CURRENT_BRANCH" \
	|| die "failed to return to branch '$CURRENT_BRANCH'"

