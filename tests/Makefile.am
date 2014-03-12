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
    ./tiny.fq \
    ./test.aligned.sorted.bam \
    ./test.aligned.sorted.bam.bai \
    ./small.fq \
    ./smallref/lobstr_test_ref.bed \
    ./smallref/readme.txt \
    ./smallref/small_lobstr_ref/lobSTR_A.fa \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_A.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAAC.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_AAAAC.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_AAAC.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_AAAT.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_AAT.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.amb \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.ann \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.bwt \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.pac \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.rbwt \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.rpac \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.rsa \
    ./smallref/small_lobstr_ref/lobSTR_ACC.fa.sa \
    ./smallref/small_lobstr_ref/lobSTR_ref.fa \
    ./smallref/small_lobstr_ref/lobSTR_strdict.txt
