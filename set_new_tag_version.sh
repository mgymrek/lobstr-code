#!/bin/sh

##
## This scripts lists the existing Git Tags, and asks to set a new git tag.
## Git tags are used to automatically set the program version.
##


##
## Show existing tags, revisions and dates
##
echo "Existing Git Tags:"
echo
for TAG in $(git tag | tac) ;
do
	echo "  Tag:	$TAG"
	git log --format="format:  Rev:	%H%n  Subj:	%s%n  Date:	%ad (%ar)" -n 1 "$TAG"
	echo
done

##
## Show the current git tag, which could be something like "2.0.3-112-abcd-dirty"
##
echo "Current Git ID:"
GITID=$(./config/git-version-gen .tarball-version) || exit 1
echo "  $GITID"
echo


##
## Ask if the user wants to create a new git tag
##
echo "Do you want to create a new git tag (i.e. mark a new version)?"
printf "Type 'yes' + enter to continue, or CTRL-C to stop: "
read NEWTAG
[ "$NEWTAG" = "yes" ] || { echo "Aborting." >&2 ; exit 1 ; }
echo
echo


##
## Ask for the git tag
##
printf "Enter new git tag, in the form of 'vX.Y.Z' (e.g. 'v2.1.13'): "
read NEWVERSION
echo "$NEWVERSION" | grep -E -q '^v[0-9]+\.[0-9]+\.[0-9]+$' || { echo "invalid version format ($NEWVERSION), expecting 'vX.Y.Z'" >&2 ; exit 1 ; }

##
## Add this tag to the git repository
##
git tag -a "$NEWVERSION" -m "New Version $NEWVERSION"
echo
echo
## Re-generate 'configure' with the new tagged version
echo "Running './reconf'..."
./reconf || exit 1
echo
echo


##
## Ask if the user wants to push (=publish) this version to GitHub.
##
echo "Do you want to publish this version (push it to github and make it public)?"
echo "NOTE: Modifying published tags is highly discouraged, so there's no going back."
echo "      If later on you want to change this tag (or re-tag another revision),"
echo "      The recommended way is to simply bump the verison and publish a new tag."
echo ""
echo "If you abort now, you can later delete the tag with 'git tag -d $NEWVERSION',"
echo "or publish the tag with 'git push --tags'"
echo ""
printf "Type 'yes' + enter to continue and publish the new tag, or CTRL-C to stop: "
read PUSHTAG
[ "$PUSHTAG" = "yes" ] || { echo "Aborting." >&2 ; exit 1 ; }
echo
echo
echo "Pushing tags"
git push --tags || exit
echo
echo
echo
echo
echo "Done!"
echo " To create new source distribution, run: ./make_src_dist.sh"
echo " To create new binary distribution, run: ./make_bin_dist.sh"
echo

