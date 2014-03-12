B1;2c#!/bin/sh

##
## This is a template to run tests through 'make check'.
## exit with zero code to signal success.
##
OUTDIR=$(mktemp -d) || exit 1
testcode() {
  if [ $? != "$@" ]; then exit 1; fi;
}
## Show Environment

echo "### Allelotype tests ###"
# Good bam file input, multiple bams
echo "Testing good bam input..."
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.aligned.sorted.bam,test.aligned.sorted.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
# Bad bam file input
echo "Testing invalid path to bam input..."
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam bad.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
# Bam file no index
echo "Testing bam file with no index..."
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.noindex.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
# Bam file with no sample in read group
echo "Testing bam file with no sample in read group..."
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.nosample.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
# Bam file with no  read group
echo "Testing bam file with no read group..."
allelotype \
  --command classify \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --strinfo smallref/smallref_strinfo.tab \
  --bam test.norg.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1

exit 0
