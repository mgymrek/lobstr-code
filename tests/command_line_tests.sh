#!/bin/sh

##
## This is a template to run tests through 'make check'.
## exit with zero code to signal success.
##
OUTDIR=$(mktemp -d 2>/dev/null || mktemp -d -t'lobstrtmp') || exit 1
testcode() {
  if [ $? != "$@" ]; then exit 1; fi;
}

echo "### lobSTR index ###"
echo "Test lobSTRIndex usage..."
lobSTRIndex >/dev/null 2>&1
testcode 1

echo "Testing building small index..."
cp ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ref.fasta ${OUTDIR}/
lobSTRIndex index \
  -a is \
  ${OUTDIR}/lobSTR_ref.fasta >/dev/null 2>&1
testcode 0

echo "### lobSTR tests ###"
echo "Testing lobSTR usage and argument cases..."
lobSTR -h >/dev/null 2>&1
testcode 1
lobSTR -? >/dev/null 2>&1
testcode 1
lobSTR --badoption >/dev/null 2>&1
testcode 1
lobSTR >/dev/null 2>&1
testcode 1
lobSTR badarg >/dev/null 2>&1
testcode 1
lobSTR --version >/dev/null 2>&1
testcode 0
echo "Testing bam input..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam  >/dev/null 2>&1
testcode 0 
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam --gzip >/dev/null 2>&1
testcode 1
echo "Testing bam input with debug info..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam --verbose --align-debug --debug >/dev/null 2>&1
testcode 0
echo "Testing bam input with quiet mode and noweb..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam --quiet --noweb >/dev/null 2>&1
testcode 0 
echo "Testing bam input with unit mode..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --rg-lib test --rg-sample test \
  --bam -u >/dev/null 2>&1
testcode 0 
echo "Testing fastq input single end..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --rg-lib test --rg-sample test --debug >/dev/null 2>&1
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
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 ${LOBSTR_TEST_DIR=.}/tmp_1.fa --p2 ${LOBSTR_TEST_DIR=.}/tmp_2.fa \
  -p 2 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -p 0 \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
echo "Testing wrong file type specified..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fa -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq --bam \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
echo "Testing file does not exist..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${OUTDIR}/del_1.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 ${OUTDIR}/del_1.fq --p2 ${OUTDIR}/del_2.fq\
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${OUTDIR}/del_1.fq \
  -q \
  -p 2 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 ${OUTDIR}/del_1.fq --p2 ${OUTDIR}/del_2.fq\
  -q \
  -p 2 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
echo "Testing additional options..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -p 1 \
  -q \
  -m 1 \
  --multi \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -r 0.1 --extend 1000\
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --fft-window-size 12 --fft-window-step 24 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --minflank 50 --maxflank 10 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -p 1 \
  -q \
  -m -1 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -p 1 \
  -q \
  --fft-window-size 24 --fft-window-step 12 --entropy-threshold 0.4 \
  --maxflank 50 --minflank 10 --max-hits-quit-aln 1000 --max-diff-ref 50 \
  --min-read-length 36 --max-read-length 1000 --mapq 100 --bwaq 15 --oldillumina \
  -g 1 -e 1 --nw-score 1 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 0
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --minflank 0 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --maxflank 0 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --max-diff-ref 0 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  -f ${LOBSTR_TEST_DIR=.}/tmp_1.fq \
  -q \
  --extend 0 \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
echo "Testing different number of input files..."
lobSTR \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --p1 file1.fq,file2.fq --p2 file2.fq \
  -q \
  --rg-lib test --rg-sample test >/dev/null 2>&1
testcode 1
echo "### Allelotype tests ###"
echo
echo "Testing good bam input..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam,${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
echo "Testing good bam input with debug option..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam,${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --out ${OUTDIR}/lobtest \
  --verbose --debug \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
echo "Testing invalid path to bam input..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/bad.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
echo "Testing bam file with no index..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.noindex.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
echo "Testing bam file with no sample in read group..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.nosample.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 0
echo "Testing bam file with no read group..."
allelotype \
  --command classify \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.norg.bam \
  --out ${OUTDIR}/lobtest \
  --verbose \
  --noise_model ${LOBSTR_TEST_DIR=.}/../models/illumina_v2.0.3 >/dev/null 2>&1
testcode 1
echo "Testing training... not enough reads"
allelotype \
  --command train \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref/small_lobstr_ref_v2/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref/smallref_strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test.aligned.sorted.bam \
  --noise_model ${OUTDIR}/lobtest \
  --haploid chrY \
  --out ${OUTDIR}/lobtest \
  --verbose >/dev/null 2>&1
testcode 1
echo "Testing training... enough reads"
allelotype \
  --command train \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref_chrY/small_lobstr_ref_v2_chrY/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref_chrY/smallref_chrY.strinfo.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test_chrY.bam \
  --noise_model ${OUTDIR}/lobtest \
  --haploid chrY \
  --out ${OUTDIR}/lobtest \
  --min-bp-before-indel 0 --min-read-end-match 0 --maximal-end-match 0 --min-border 5 \
  --verbose --debug >/dev/null 2>&1
testcode 0
echo "Testing training... invalid strinfo"
allelotype \
  --command train \
  --index-prefix ${LOBSTR_TEST_DIR=.}/smallref_chrY/small_lobstr_ref_v2_chrY/lobSTR_ \
  --strinfo ${LOBSTR_TEST_DIR=.}/smallref_chrY/bad_strinfo_file.tab \
  --bam ${LOBSTR_TEST_DIR=.}/test_chrY.bam \
  --noise_model ${OUTDIR}/lobtest2 \
  --haploid chrY \
  --out ${OUTDIR}/lobtest \
  --verbose --debug >/dev/null 2>&1
testcode 1
exit 0
