# Copyright (C) 2014 Assaf Gordon <agordon@wi.mit.edu>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

##
## Add tests scripts here
##
TESTS = command_line_tests.sh

# (Copied from coreutils' tests/local.mk)
# Use automake's TESTS_ENVIRONMENT to add ./src/ to the search path.
TESTS_ENVIRONMENT =  \
  export \
  PATH='$(abs_top_builddir)/src$(PATH_SEPARATOR)'"$$PATH" \
  LOBSTR_TEST_DIR='$(abs_top_srcdir)/tests' \
  ;

EXTRA_DIST = \
    $(TESTS) \
    ./runtest.sh \
    ../models/illumina_v2.0.3.stepmodel \
    ../models/illumina_v2.0.3.stuttermodel \
    ../models/illumina_v2.0.3.stutterproblem \
    ./test.aligned.sorted.bam \
    ./test.aligned.sorted.bam.bai \
    ./test.fa.gz \
    ./test.fq.gz \
    ./test.noindex.bam \
    ./test.norg.bam \
    ./test.norg.bam.bai \
    ./test.nosample.bam \
    ./test.nosample.bam.bai \
    ./tiny.fq \
    ./tmp_1.fa \
    ./tmp_2.fa \
    ./tmp_1.fq \
    ./tmp_2.fq \
    ./small.fq \
    ./smallref/lobstr_test_ref.bed \
    ./smallref/readme.txt \
    ./smallref/smallref_strinfo.tab \
    ./smallref/small_lobstr_ref_v2/lobSTR_chromsizes.tab \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.ann \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.pac \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.rpac \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.sa \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.amb \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.bwt \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.rbwt \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref.fasta.rsa \
    ./smallref/small_lobstr_ref_v2/lobSTR_mergedref.bed \
    ./smallref/small_lobstr_ref_v2/lobSTR_ref_map.tab \
    ./smallref_chrY/bad_strinfo_file.tab \
    ./smallref_chrY/small_lobstr_ref_v2_chrY/lobSTR_chrY_ref.fasta \
    ./smallref_chrY/small_lobstr_ref_v2_chrY/lobSTR_chrY_ref_map.tab \
    ./smallref_chrY/smallref_chrY.strinfo.tab \
    ./test_chrY.bam \
    ./test_chrY.bam.bai
