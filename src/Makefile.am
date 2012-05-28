
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

AM_CPPFLAGS = $(FFTW_CFLAGS) $(RINSIDECPPFLAGS) $(RCPPCPPFLAGS) $(RCPPFLAGS)

noinst_headers = \
	BamFileReader.h \
	BamPairedFileReader.h \
	BWAReadAligner.h \
	CachedWindowGenerator.h \
	cigar.h \
	common.h \
	EntropyDetection.h \
	FastaFileReader.h \
	FastaPairedFileReader.h \
	FastqFileReader.h \
	FastqPairedFileReader.h \
	FFT_nuc_vectors.h \
	FFT_four_nuc_vectors.h \
	gzstream.h \
	IFileReader.h \
	ISatellite.h \
	MSReadRecord.h \
	MultithreadData.h \
	nw.h \
	HammingWindowGenerator.h \
	TukeyWindowGenerator.h \
	ReadPair.h \
	runtime_parameters.h \
	SamFileWriter.h \
	STRDetector.h \
	TabFileWriter.h \
	TextFileReader.h \
	TextFileWriter.h \
	ZippedFastaFileReader.h \
	ZippedFastqFileReader.h \
	ZippedTextFileReader.h \
	bamlite.h \
	bntseq.h \
	bwase.h \
	bwt.h \
	bwtaln.h \
	bwtgap.h \
	utils.h \
	bwt_gen.h \
	bwt_lite.h \
	kstring.h \
	stdaln.h \
	api/algorithms/Sort.h \
	api/api_global.h \
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
	bwtaln.h \
	main.h \
	bwt.h \
	bamlite.h \
	bntseq.h \
	bwtsw2.h \
	stdaln.h \
	bwtgap.h \
	utils.h \
	kstring.h \
	bwt_lite.h \
	bwase.h \
	bwt_gen.h \
	Genotyper.h \
	NoiseModel.h \
	ReadContainer.h 

noinst_LIBRARIES = liblobstr.a libindex.a liballelotype.a 

liblobstr_a_SOURCES = \
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
	gzstream.cpp \
	HammingWindowGenerator.cpp \
	MultithreadData.cpp \
	nw.cpp \
	runtime_parameters.cpp \
	SamFileWriter.cpp \
	STRDetector.cpp \
	TabFileWriter.cpp \
	TextFileReader.cpp \
	TextFileWriter.cpp \
	TukeyWindowGenerator.cpp \
	ZippedFastaFileReader.cpp \
	ZippedFastqFileReader.cpp \
	ZippedTextFileReader.cpp \
	bamlite.c \
	bntseq.c \
	bwaseqio.c \
	bwape.c \
	bwase.c \
	bwt.c \
	bwtaln.c \
	bwtgap.c \
	bwtindex.c \
	bwtio.c \
	bwtmisc.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwtsw2_core.c \
	bwtsw2_main.c \
	bwt_gen.c \
	bwt_lite.c \
	cs2nt.c \
	is.c \
	kstring.c \
	QSufSort.c \
	simple_dp.c \
	stdaln.c \
	utils.c  \
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
	utils/bamtools_pileup_engine.cpp

libindex_a_SOURCES = \
	bwtaln.c \
	bwt.c \
	bamlite.c \
	stdaln.c \
	bwtgap.c \
	bwaseqio.c \
	bwtindex.c \
	bwtio.c \
	utils.c \
	simple_dp.c \
	kstring.c \
	bwtsw2_core.c \
	bwtsw2_main.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwt_lite.c \
	bwape.c \
	bwase.c \
	bwt_gen.c \
	QSufSort.c \
	bntseq.c \
	bwtmisc.c \
	is.c \
	cs2nt.c

liballelotype_a_SOURCES = \
	common.cpp \
	Genotyper.cpp \
	NoiseModel.cpp \
	ReadContainer.cpp \
	runtime_parameters.cpp \
	TextFileReader.cpp \
	TextFileWriter.cpp \
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
	utils/bamtools_pileup_engine.cpp

#-- scripts to be installed in $prefix/bin
dist_bin_SCRIPTS = \
	../scripts/lobstr_index.py
	../scripts/convertBam.py

CLEANFILES = $(bin_SCRIPTS)

lobSTRIndex_SOURCES = main.c
lobSTRIndex_CPPFLAGS = $(AM_CPPFLAGS)
lobSTRIndex_LDADD = $(AM_LDADD) libindex.a
lobSTRIndex_LDFLAGS = -static  $(AM_LDFLAGS)

lobSTR_SOURCES = main.cpp
lobSTR_CPPFLAGS = $(AM_CPPFLAGS)
lobSTR_LDADD = liblobstr.a $(FFTW_LIBS)
lobSTR_LDFLAGS = -static $(AM_LDFLAGS)


allelotype_SOURCES = 	main_allelotype.cpp
allelotype_CPPFLAGS = $(AM_CPPFLAGS) $(RCPPFLAGS) $(RINSIDECPPFLAGS) $(RCPPCPPFLAGS) 
allelotype_LDADD = liballelotype.a $(RLDFLAGS)  $(VAR_STATIC) $(RINSIDELDFLAGS) $(RCPPLDFLAGS)  $(VAR_DYNAMIC)

