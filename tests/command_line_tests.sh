#!/bin/sh

##
## This is a template to run tests through 'make check'.
## exit with zero code to signal success.
##
OUTDIR=$(mktemp -d) || exit 1
testcode() {
  if [ $? != "$@" ]; then exit 1; fi;
}

## Show Environment

echo "### lobSTR tests ###"
echo "Testing bam input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam  >/dev/null 2>&1
testcode 0 
echo "Testing fastq input single end..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fastq input paired end..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 ${LOBSTR_TEST_DIR=.}/tmp_1.fq --p2 ${LOBSTR_TEST_DIR=.}/tmp_2.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing bam paired end input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.nosample.bam \
  --bampair \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing zipped fastq input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.fq.gz \
  -q \
  --gzip \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fasta single end input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fa \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fasta paired end input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 ${LOBSTR_TEST_DIR=.}/tmp_1.fa --p2 ${LOBSTR_TEST_DIR=.}/tmp_2.fa \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing gzipped fasta input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.fa.gz \
  --gzip \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fastq input single end with multithread..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -p 2 \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "### Allelotype tests ###"
echo
echo "Testing good bam input..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam,${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
echo "Testing invalid path to bam input..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/bad.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
echo "Testing bam file with no index..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.noindex.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
echo "Testing bam file with no sample in read group..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.nosample.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
echo "Testing bam file with no read group..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.norg.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1

exit 0
