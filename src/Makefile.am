
# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = lobSTR \
	lobSTRIndex \
	allelotype

AM_LDFLAGS = $(BOOST_LDFLAGS) $(ZLIB_LDFLAGS)

AM_CPPFLAGS = $(FFTW_CFLAGS) $(PROTOBUF_CFLAGS) -I$(top_srcdir)/src/ -I$(top_srcdir)/

noinst_LIBRARIES = liblobstr.a libindex.a liballelotype.a libbamtools.a

libbamtools_a_SOURCES = \
	api/BamReader.cpp \
	api/BamMultiReader.cpp \
	api/BamWriter.cpp \
	api/BamAlignment.cpp \
	api/SamHeader.cpp \
	api/SamProgram.cpp \
	api/SamProgramChain.cpp \
	api/SamSequence.cpp \
	api/SamSequenceDictionary.cpp \
	api/SamReadGroup.cpp \
	api/SamReadGroupDictionary.cpp \
	api/internal/BamWriter_p.cpp \
	api/internal/SamFormatPrinter_p.cpp \
	api/internal/BamException_p.cpp \
	api/internal/BamDeviceFactory_p.cpp \
	api/internal/BamHeader_p.cpp \
	api/internal/BamRandomAccessController_p.cpp \
	api/internal/BamReader_p.cpp \
	api/internal/ILocalIODevice_p.cpp \
	api/internal/BgzfStream_p.cpp \
	api/internal/BamFile_p.cpp \
	api/internal/SamFormatParser_p.cpp \
	api/internal/SamHeaderValidator_p.cpp \
	api/internal/BamFtp_p.cpp \
	api/internal/BamHttp_p.cpp \
	api/internal/BamIndexFactory_p.cpp \
	api/internal/BamPipe_p.cpp \
	api/internal/BamStandardIndex_p.cpp \
	api/internal/BamToolsIndex_p.cpp \
	api/internal/BamMultiReader_p.cpp \
	utils/bamtools_utilities.cpp  \
	utils/bamtools_fasta.cpp \
	utils/bamtools_options.cpp \
	utils/bamtools_pileup_engine.cpp \
	api/BGZF.h \
	api/BamAlgorithms.h \
	api/algorithms/Sort.h \
	api/api_global.h \
	api/BamIndex.h \
	api/SamConstants.h \
	api/BamReader.h \
	api/BamMultiReader.h \
	api/BamConstants.h \
	api/BamWriter.h \
	api/BamAlignment.h \
	api/BamAux.h \
	api/IBamIODevice.h \
	api/SamHeader.h \
	api/SamProgram.h \
	api/SamSequence.h \
	api/SamSequenceDictionary.h \
	api/SamReadGroup.h \
	api/SamReadGroupDictionary.h \
	api/internal/BamWriter_p.h \
	api/SamProgramChain.h \
	api/internal/IRemoteIODevice_p.h \
	api/internal/IRemoteIODevice_p.cpp \
	api/internal/SamFormatPrinter_p.h \
	api/internal/BamException_p.h \
	api/internal/BamDeviceFactory_p.h \
	api/internal/BamHeader_p.h \
	api/internal/BamRandomAccessController_p.h \
	api/internal/BamReader_p.h \
	api/internal/ILocalIODevice_p.h \
	api/internal/BgzfStream_p.h \
	api/internal/BamFile_p.h \
	api/internal/SamFormatParser_p.h \
	api/internal/SamHeaderValidator_p.h \
	api/internal/BamFtp_p.h \
	api/internal/BamHttp_p.h \
	api/internal/BamIndexFactory_p.h \
	api/internal/BamMultiMerger_p.h \
	api/internal/BamPipe_p.h \
	api/internal/BamStandardIndex_p.h \
	api/internal/BamToolsIndex_p.h \
	api/internal/SamHeaderVersion_p.h \
	api/internal/BamMultiReader_p.h \
	shared/bamtools_global.h \
	utils/bamtools_fasta.h \
	utils/bamtools_filter_engine.h \
	utils/bamtools_filter_properties.h \
	utils/bamtools_filter_ruleparser.h \
	utils/bamtools_options.h \
	utils/bamtools_utilities.h \
	utils/bamtools_pileup_engine.h \
	utils/bamtools_variant.h \
	utils/utils_global.h

liblobstr_a_SOURCES = \
	BamFileReader.cpp BamFileReader.h \
	BamPairedFileReader.cpp BamPairedFileReader.h \
	BWAReadAligner.cpp BWAReadAligner.h \
	CachedWindowGenerator.cpp CachedWindowGenerator.h \
	common.cpp common.h \
	EntropyDetection.cpp EntropyDetection.h \
	FastaFileReader.cpp FastaFileReader.h \
	FastaPairedFileReader.cpp FastaPairedFileReader.h \
	FastqFileReader.cpp FastqFileReader.h \
	FastqPairedFileReader.cpp FastqPairedFileReader.h \
	FFT_four_nuc_vectors.cpp FFT_four_nuc_vectors.h \
	FFT_nuc_vectors.cpp FFT_nuc_vectors.h \
	FilterCounter.cpp FilterCounter.h \
	gzstream.cpp gzstream.h \
	STRRecord.h \
	HammingWindowGenerator.cpp HammingWindowGenerator.h \
	MultithreadData.cpp MultithreadData.h \
	nw.cpp nw.h \
	NWNoRefEndPenalty.cpp NWNoRefEndPenalty.h \
	cigar.h \
	runtime_parameters.cpp runtime_parameters.h \
	IFileReader.h RunInfo.h \
	IFileWriter.h ISatellite.h MSReadRecord.h \
	SamFileWriter.cpp SamFileWriter.h \
	STRDetector.cpp STRDetector.h \
	TextFileReader.cpp TextFileReader.h \
	TextFileWriter.cpp TextFileWriter.h \
	TukeyWindowGenerator.cpp TukeyWindowGenerator.h \
	ZippedFastaFileReader.cpp ZippedFastaFileReader.h \
	ZippedFastqFileReader.cpp ZippedFastqFileReader.h \
	ZippedTextFileReader.cpp ZippedTextFileReader.h \
	api/BamReader.h api/BamMultiReader.h \
	api/BamAlignment.h api/BamWriter.h \
	api/SamHeader.h api/SamReadGroup.h \
	api/SamReadGroupDictionary.h

libindex_a_SOURCES = \
	bwtaln.c bwtaln.h \
	bwt.c bwt.h \
	bamlite.c \
	bamlite.h \
	stdaln.c stdaln.h \
	bwtgap.c bwtgap.h \
	bwaseqio.c \
	bwaseqio.h \
	bwtindex.c \
	bwtio.c \
	utils.c \
	utils.h \
	simple_dp.c \
	kstring.c kstring.h \
	bwtsw2_core.c \
	khash.h kseq.h ksort.h kvec.h \
	bwtsw2_main.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwt_lite.c bwt_lite.h \
	bwape.c \
	bwase.c \
	bwase.h \
	bwt_gen.c bwt_gen.h \
	bwtsw2.h \
	QSufSort.c QSufSort.h \
	bntseq.c bntseq.h \
	bwtmisc.c \
	is.c \
	cs2nt.c

liballelotype_a_SOURCES = \
	AlignedRead.h \
	AlignmentFilters.h AlignmentFilters.cpp \
	common.cpp common.h \
	BamFileReader.cpp BamFileReader.h \
	BamPairedFileReader.cpp BamPairedFileReader.h \
	FastaFileReader.cpp FastaFileReader.h \
	FastqFileReader.cpp FastqFileReader.h \
	FastaPairedFileReader.cpp FastaPairedFileReader.h \
	FastqPairedFileReader.cpp FastqPairedFileReader.h \
	FilterCounter.cpp FilterCounter.h \
	Genotyper.cpp Genotyper.h \
	gzstream.cpp gzstream.h \
	NoiseModel.cpp NoiseModel.h \
	linear.h linear.cpp tron.cpp tron.h \
	ReadContainer.cpp ReadContainer.h \
	ReadPair.h \
	ReferenceSTR.cpp ReferenceSTR.h \
	RemoveDuplicates.cpp RemoveDuplicates.h \
	STRRecord.h \
	runtime_parameters.cpp runtime_parameters.h \
	IFileReader.h RunInfo.h \
	IFileWriter.h ISatellite.h MSReadRecord.h \
	TextFileReader.cpp TextFileReader.h \
	TextFileWriter.cpp TextFileWriter.h \
	VCFWriter.cpp VCFWriter.h \
	ZAlgorithm.cpp ZAlgorithm.h \
	ZippedFastaFileReader.cpp ZippedFastaFileReader.h \
	ZippedFastqFileReader.cpp ZippedFastqFileReader.h \
	ZippedTextFileReader.cpp ZippedTextFileReader.h \
	linear.cpp linear.h \
	predict.c \
	train.c	\
	tron.cpp tron.h \
	api/BamReader.h api/BamMultiReader.h \
	api/BamAlignment.h api/BamWriter.h \
	api/SamHeader.h api/SamReadGroup.h \
	api/SamReadGroupDictionary.h

lobSTRIndex_SOURCES = main.c main.h
lobSTRIndex_CPPFLAGS = $(AM_CPPFLAGS)
lobSTRIndex_LDADD = $(AM_LDADD) libindex.a libbamtools.a
lobSTRIndex_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)

lobSTR_SOURCES = main.cpp
lobSTR_CPPFLAGS = $(AM_CPPFLAGS)
lobSTR_LDADD = liblobstr.a libbamtools.a libindex.a $(FFTW_LIBS) $(PROTOBUF_LIBS)
lobSTR_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)

allelotype_SOURCES = main_allelotype.cpp
allelotype_CPPFLAGS = $(AM_CPPFLAGS)
allelotype_LDADD = liballelotype.a libbamtools.a $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PROTOBUF_LIBS)
allelotype_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)

# (Copied from coreutils' tests/local.mk)
# Use automake's TESTS_ENVIRONMENT to point to the "tests" directory
TESTS_ENVIRONMENT =  \
  export \
  LOBSTR_TEST_DIR='$(abs_top_srcdir)/tests' \
  ;

TESTS = lobSTRTests
check_PROGRAMS = $(TESTS)
lobSTRTests_SOURCES = \
	tests/main_test.cpp \
	tests/AlignmentFilters_test.h \
	tests/AlignmentFilters_test.cpp \
	tests/BWAReadAligner_test.h \
	tests/BWAReadAligner_test.cpp \
	tests/common_test.h \
	tests/common_test.cpp \
	tests/DNATools.h \
	tests/DNATools.cpp \
	tests/NWNoRefEndPenalty_test.h \
	tests/NWNoRefEndPenalty_test.cpp \
	tests/ReadContainer_test.h \
	tests/ReadContainer_test.cpp \
	tests/RemoveDuplicates_test.h \
	tests/RemoveDuplicates_test.cpp \
	tests/VCFWriter_test.h \
	tests/VCFWriter_test.cpp \
	tests/ZAlgorithm_test.h \
	tests/ZAlgorithm_test.cpp \
	AlignmentFilters.cpp \
	BamFileReader.cpp \
	BamPairedFileReader.cpp \
	BWAReadAligner.cpp \
	CachedWindowGenerator.cpp \
	common.cpp \
	EntropyDetection.cpp \
	FastaFileReader.cpp \
	FastaPairedFileReader.cpp \
	FastqFileReader.cpp \
	FastqPairedFileReader.cpp \
	FFT_four_nuc_vectors.cpp \
	FFT_nuc_vectors.cpp \
	FilterCounter.h FilterCounter.cpp \
	gzstream.cpp \
	HammingWindowGenerator.cpp \
	MultithreadData.cpp \
	nw.cpp nw.h \
	NWNoRefEndPenalty.cpp NWNoRefEndPenalty.h \
	cigar.h \
	runtime_parameters.cpp \
	SamFileWriter.cpp \
	STRDetector.cpp \
	TextFileReader.cpp \
	TextFileWriter.cpp \
	TukeyWindowGenerator.cpp \
	ReadContainer.cpp \
	RemoveDuplicates.cpp \
	VCFWriter.cpp \
	ZAlgorithm.h \
	ZAlgorithm.cpp \
	ZippedFastaFileReader.cpp \
	ZippedFastqFileReader.cpp \
	ZippedTextFileReader.cpp \
	api/BamReader.h api/BamMultiReader.h \
	api/BamAlignment.h api/BamWriter.h \
	api/SamHeader.h api/SamReadGroup.h \
	api/SamReadGroupDictionary.h

lobSTRTests_CPPFLAGS = $(CPPUNIT_CFLAGS) $(lobSTR_CPPFLAGS)
lobSTRTests_LDADD = libbamtools.a libindex.a $(FFTW_LIBS) $(CPPUNIT_LIBS) -ldl
lobSTRTests_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)
