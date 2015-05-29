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
# A validation suite to check concordance with gold standard
# capillary electrophoresis calls.
#

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@">&2
    exit 1
}

usage()
{
    BASE=$(basename -- "$0")
    echo "allelotype validation suite
Usage:
    $BASE REFDIR SAMPLEDIR [ARGS]

REFDIR gives base directory of lobSTR hg19 reference bundle (e.g. hg19_v3.0.2/)
SAMPLEDIR gives base directory of validation data (e.g. lobSTR_validation_data_v3.0.2)
ARGS gives non-default options to pass to allelotype

The reference bundle and validation data can be downloaded from the lobSTR website:
http://lobstr.teamerlich.org/download.html
"
    exit 1
}

REFDIR="$1"
SAMPLEDIR="$2"
ARGS="$3"
test -z "${REFDIR}" || test -z "${SAMPLEDIR}" && usage

echo "[INPUT]: REFDIR=${REFDIR}"
echo "[INPUT]: SAMPLEDIR=${SAMPLEDIR}"
echo "[INPUT]: ARGS=${ARGS}"

##
## Check that required programs are installed
##
command -v vcf-sort >/dev/null 2>&1 || die "Could not find vcf-sort in the PATH"
command -v bgzip >/dev/null 2>&1 || die "Could not find bgzip in the PATH"
command -v tabix >/dev/null 2>&1 || die "Could not find tabix in the PATH"

##
## These directories are defined in Makefile.am,
## and created during "make install".
##
DIR=$(dirname -- "$0")
BINDIR="${DIR}/../src/"

## Ensure directories exist
test -d "${BINDIR}" \
    || die "installation error: bin-dir '${BINDIR}' not found"
test -d "${REFDIR}" \
    || die "input error: ref-dir '${REFDIR}' not found"
test -d "${SAMPLEDIR}" \
    || die "input error: sample-dir '${SAMPLEDIR}' not found"

##
## Ensure expected files exists
##
ALLELOTYPEBIN="${BINDIR}/allelotype"
SAMPLEBAM="${SAMPLEDIR}/marshfield.sorted.bam"
CORRECTIONS="${SAMPLEDIR}/marshfield_corrections.tab"
CAPILLARY="${SAMPLEDIR}/marshfield_capillary.stru"
CONVERSIONS="${SAMPLEDIR}/sample_conversions.tab"

test -e "${ALLELOTYPEBIN}" \
    || die "installation error: allelotype binary (${ALLELOTYPE}) not found"
test -x "${ALLELOTYPEBIN}" \
    || die "installation error: allelotype binary (${ALLELOTYPE}) not executable"
test -e "${SAMPLEBAM}" \
    || die "input error: sample bam file (${SAMPLEBAM}) not found"
test -e "${SAMPLEBAM}.bai" \
    || die "input error: bam index for ${SAMPLEBAM} not found"
test -e "${CORRECTIONS}" \
    || die "input error: locus corrections file (${CORRECTIONS}) not found"
test -e "${CAPILLARY}" \
    || die "input error: capillary genotypes file (${CAPILLARY}) not found"
test -e "${CONVERSIONS}" \
    || die "input error: sample conversions file (${CONVERSIONS}) not found"

TMPDIR=$(mktemp -d) || die "failed to create temporary directory"
trap "rm -r '${TMPDIR}'" EXIT

REFPREFIX="$(find ${REFDIR} -maxdepth 1 -mindepth 1 -type d)/lobSTR_"
STRINFO="$(ls ${REFDIR}/*strinfo*)"
NOISE_MODEL="${DIR}/../models/illumina_v2.0.3"

####################################
# Run allelotype on bwamem alignments
echo "[PROGRESS]: Running allelotype on bwa-mem alignments"

"${ALLELOTYPEBIN}" \
  --verbose \
  --command classify \
  --bam "${SAMPLEBAM}" \
  --noise_model "${NOISE_MODEL}" \
  --index-prefix "${REFPREFIX}" \
  --strinfo "${STRINFO}" \
  --out "${TMPDIR}/validation_bwamem" \
  --noweb ${ARGS} >/dev/null 2>&1 \
    || die "allelotype failed (check error messages above)"

vcf-sort "${TMPDIR}/validation_bwamem.vcf" | bgzip -c > "${TMPDIR}/validation_bwamem.vcf.gz"
tabix -p vcf "${TMPDIR}/validation_bwamem.vcf.gz"

####################################
# Compare to capillary - bwamem alignments
echo "[RESULTS]: Comparison of bwamem alignment calls to capillary"
"${DIR}/lobSTR_vcf_to_tab.py" "${TMPDIR}/validation_bwamem.vcf.gz" > "${TMPDIR}/validation_bwamem.tab" \
    || die "lobSTR_vcf_to_tab.py failed"
"${DIR}/lobSTR_capillary_comparator.py" \
    --lobSTR "${TMPDIR}/validation_bwamem.tab" \
    --cap "${CAPILLARY}" \
    --corrections "${CORRECTIONS}" \
    --sample-conversions "${CONVERSIONS}" || die "lobSTR_capillary_comparator.py failed"
