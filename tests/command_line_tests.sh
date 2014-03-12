#!/bin/sh

##
## This is a template to run tests through 'make check'.
## exit with zero code to signal success.
##

testcode() {
  if [ $? != "$@" ]; then exit 1; fi;
}
## Show Environment
echo "Running Dummy Test"
echo
echo "### Allelotype tests ###"

# Good bam file input, multiple bams
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.aligned.sorted.bam,test.aligned.sorted.bam \
  --out /tmp/lobtest \
  --noise_model ../models/illumina_v2.0.3
testcode 0
# Bad bam file input
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam bad.bam \
  --out /tmp/lobtest \
  --noise_model ../models/illumina_v2.0.3
testcode 1
# Bam file no index
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.noindex.bam \
  --out /tmp/lobtest \
  --noise_model ../models/illumina_v2.0.3
testcode 1
echo
echo
echo
env
echo
echo
echo
which lobSTR
echo
echo
echo
lobSTR --help
echo
echo
echo

exit 0
