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
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam >/dev/null 2>&1
testcode 0 
echo "Testing fastq input single end..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f tmp_1.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fastq input paired end..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 tmp_1.fq --p2 tmp_2.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing bam paired end input..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f test.nosample.bam \
  --bampair \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing zipped fastq input..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f test.fq.gz \
  -q \
  --gzip \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fasta single end input..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f tmp_1.fa \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fasta paired end input..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 tmp_1.fa --p2 tmp_2.fa \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing gzipped fasta input..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f test.fa.gz \
  --gzip \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing fastq input single end with multithread..."
lobSTR \
  --index-prefix smallref/small_lobstr_ref/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f tmp_1.fq \
  -p 2 \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "### Allelotype tests ###"
echo
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
