
# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = lobSTR lobSTRIndex allelotype



AM_CPPFLAGS = $(FFTW_CFLAGS)

AM_LDADD = $(FFTW_LIBS)
AM_LDFLAGS = $(BOOST_LDFLAGS) $(ZLIB_LDFLAGS)

lobSTR_CPPFLAGS = $(AM_CPPFLAGS)
lobSTR_LDADD = $(AM_LDADD)
lobSTR_LDFLAGS = $(AM_LDFLAGS)

lobSTR_SOURCES = \
	main.cpp \
	BamFileReader.h \
	BamFileReader.cpp \
	BWAReadAligner.h \
	BWAReadAligner.cpp \
	CachedWindowGenerator.h \
	CachedWindowGenerator.cpp \
	cigar.h \
	common.h \
	common.cpp \
	EntropyDetection.h \
	EntropyDetection.cpp \
	FastaFileReader.h \
	FastaFileReader.cpp \
	FastqFileReader.h \
	FastqFileReader.cpp \
	FFT_four_nuc_vectors.h \
	FFT_four_nuc_vectors.cpp \
	FFT_nuc_vectors.h \
	FFT_nuc_vectors.cpp \
	HammingWindowGenerator.h \
	HammingWindowGenerator.cpp \
	IFileReader.h \
	ISatellite.h \
	MSReadRecord.h \
	MultithreadData.h \
	MultithreadData.cpp \
	nw.h \
	nw.cpp \
	runtime_parameters.cpp \
	runtime_parameters.h \
	SamFileWriter.cpp \
	SamFileWriter.h \
	STRDetector.h \
	STRDetector.cpp \
	TabFileWriter.h \
	TabFileWriter.cpp \
	TextFileReader.h \
	TextFileReader.cpp \
	TextFileWriter.h \
	TextFileWriter.cpp \
	TukeyWindowGenerator.h \
	TukeyWindowGenerator.cpp \
	bamlite.h \
	bamlite.c \
	bntseq.h \
	bntseq.c \
	bwaseqio.c \
	bwape.c \
	bwase.c \
	bwase.h \
	bwt.h \
	bwt.c \
	bwtaln.h \
	bwtaln.c \
	bwtgap.h \
	bwtgap.c \
	bwtindex.c \
	bwtio.c \
	bwtmisc.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwtsw2_core.c \
	bwtsw2_main.c \
	bwt_gen.c \
	bwt_gen.h \
	bwt_lite.c \
	bwt_lite.h \
	cs2nt.c \
	is.c \
	kstring.c \
	kstring.h \
	QSufSort.c \
	simple_dp.c \
	stdaln.h \
	stdaln.c \
	utils.h \
	utils.c  \
	api/algorithms/Sort.h \
	api/api_global.h \
	api/BamReader.h \
	api/BamReader.cpp \
	api/BamMultiReader.h \
	api/BamMultiReader.cpp \
	api/BamConstants.h \
	api/BamWriter.h \
	api/BamWriter.cpp \
	api/BamAlignment.h \
	api/BamAlignment.cpp \
	api/BamAux.h \
	api/IBamIODevice.h \
	api/SamHeader.h \
	api/SamHeader.cpp \
	api/SamProgram.cpp \
	api/SamProgram.h \
	api/SamProgramChain.cpp \
	api/SamProgramChain.h \
	api/SamSequence.cpp \
	api/SamSequence.h \
	api/SamSequenceDictionary.h \
	api/SamSequenceDictionary.cpp \
	api/SamReadGroup.h \
	api/SamReadGroup.cpp \
	api/SamReadGroupDictionary.h \
	api/SamReadGroupDictionary.cpp \
	api/internal/BamWriter_p.h \
	api/internal/BamWriter_p.cpp \
	api/internal/SamFormatPrinter_p.h \
	api/internal/SamFormatPrinter_p.cpp \
	api/internal/BamException_p.h \
	api/internal/BamException_p.cpp \
	api/internal/BamDeviceFactory_p.h \
	api/internal/BamDeviceFactory_p.cpp \
	api/internal/BamHeader_p.h \
	api/internal/BamHeader_p.cpp \
	api/internal/BamRandomAccessController_p.cpp \
	api/internal/BamRandomAccessController_p.h \
	api/internal/BamReader_p.cpp \
	api/internal/BamReader_p.h \
	api/internal/ILocalIODevice_p.h \
	api/internal/ILocalIODevice_p.cpp \
	api/internal/BgzfStream_p.cpp \
	api/internal/BgzfStream_p.h \
	api/internal/BamFile_p.h \
	api/internal/BamFile_p.cpp \
	api/internal/SamFormatParser_p.h \
	api/internal/SamFormatParser_p.cpp \
	api/internal/SamHeaderValidator_p.h \
	api/internal/SamHeaderValidator_p.cpp \
	api/internal/BamFtp_p.h \
	api/internal/BamFtp_p.cpp \
	api/internal/BamHttp_p.h \
	api/internal/BamHttp_p.cpp \
	api/internal/BamIndexFactory_p.h \
	api/internal/BamIndexFactory_p.cpp \
	api/internal/BamMultiMerger_p.h \
	api/internal/BamPipe_p.h \
	api/internal/BamPipe_p.cpp \
	api/internal/BamStandardIndex_p.cpp \
	api/internal/BamStandardIndex_p.h \
	api/internal/BamToolsIndex_p.h \
	api/internal/BamToolsIndex_p.cpp \
	api/internal/SamHeaderVersion_p.h \
	api/internal/BamMultiReader_p.h \
	api/internal/BamMultiReader_p.cpp \
	shared/bamtools_global.h \
	utils/bamtools_utilities.h \
	utils/bamtools_utilities.cpp  \
	utils/bamtools_fasta.cpp \
	utils/bamtools_fasta.h \
	utils/bamtools_filter_engine.h \
	utils/bamtools_filter_properties.h \
	utils/bamtools_filter_ruleparser.h \
	utils/bamtools_options.h \
	utils/bamtools_options.cpp \
	utils/bamtools_pileup_engine.h \
	utils/bamtools_pileup_engine.cpp \
	utils/bamtools_variant.h \
	utils/utils_global.h

lobSTRIndex_CPPFLAGS = $(AM_CPPFLAGS)
lobSTRIndex_LDADD = $(AM_LDADD)
lobSTRIndex_LDFLAGS = $(AM_LDFLAGS)

lobSTRIndex_SOURCES = \
	main.c main.h \
	bwtaln.h \
	bwtaln.c \
	bwt.h bwt.c \
	bamlite.h \
	bamlite.c \
	bntseq.h \
	bwtsw2.h \
	stdaln.h \
	stdaln.c \
	bwtgap.h \
	bwtgap.c \
	bwaseqio.c \
	bwtindex.c \
	bwtio.c \
	utils.h \
	utils.c \
	simple_dp.c \
	kstring.c \
	kstring.h \
	bwtsw2_core.c \
	bwtsw2_main.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwt_lite.c \
	bwt_lite.h \
	bwape.c \
	bwase.c \
	bwase.h \
	bwt_gen.c \
	bwt_gen.h \
	QSufSort.c \
	bntseq.c \
	bwtmisc.c \
	is.c \
	cs2nt.c

# Required for R (R CMD config --cppflags, --ldflags)
#RCPPFLAGS = -I/usr/share/R/include 
#RLDFLAGS = -L/usr/lib64/R/lib -lR

# Required for Rcpp (echo 'Rcpp:::LdFlags()' | R --vanilla --slave)
#RCPPCPPFLAGS = -I/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/Rcpp/include
#RCPPLDFLAGS = -L/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/Rcpp/lib -lRcpp -Wl,-rpath,/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/Rcpp/lib

# Required for RInside (echo 'RInside:::CxxFlags()' | R --vanilla --slave)
#RINSIDECPPFLAGS = -I/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/RInside/include
#RINSIDELDFLAGS = -L/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/RInside/lib -lRInside -Wl,-rpath,/home/mgymrek/R/x86_64-pc-linux-gnu-library/2.14/RInside/lib

allelotype_CPPFLAGS = $(AM_CPPFLAGS) $(RCPPFLAGS) $(RCPPCPPFLAGS) $(RINSIDECPPFLAGS)
allelotype_LDADD = $(AM_LDADD) 
allelotype_LDFLAGS = $(AM_LDFLAGS) $(RLDFLAGS) $(RCPPLDFLAGS) $(RINSIDELDFLAGS)

allelotype_SOURCES = main_allelotype.cpp \
	Genotyper.h \
	Genotyper.cpp \
	NoiseModel.h \
	NoiseModel.cpp \
	ReadContainer.h \
	ReadContainer.cpp \
	runtime_parameters.h \
	runtime_parameters.cpp \
	TextFileReader.h \
	TextFileReader.cpp \
	TextFileWriter.h \
	TextFileWriter.cpp \
	api/algorithms/Sort.h \
	api/api_global.h \
	api/BamReader.h \
	api/BamReader.cpp \
	api/BamMultiReader.h \
	api/BamMultiReader.cpp \
	api/BamConstants.h \
	api/BamWriter.h \
	api/BamWriter.cpp \
	api/BamAlignment.h \
	api/BamAlignment.cpp \
	api/BamAux.h \
	api/IBamIODevice.h \
	api/SamHeader.h \
	api/SamHeader.cpp \
	api/SamProgram.cpp \
	api/SamProgram.h \
	api/SamProgramChain.cpp \
	api/SamProgramChain.h \
	api/SamSequence.cpp \
	api/SamSequence.h \
	api/SamSequenceDictionary.h \
	api/SamSequenceDictionary.cpp \
	api/SamReadGroup.h \
	api/SamReadGroup.cpp \
	api/SamReadGroupDictionary.h \
	api/SamReadGroupDictionary.cpp \
	api/internal/BamWriter_p.h \
	api/internal/BamWriter_p.cpp \
	api/internal/SamFormatPrinter_p.h \
	api/internal/SamFormatPrinter_p.cpp \
	api/internal/BamException_p.h \
	api/internal/BamException_p.cpp \
	api/internal/BamDeviceFactory_p.h \
	api/internal/BamDeviceFactory_p.cpp \
	api/internal/BamHeader_p.h \
	api/internal/BamHeader_p.cpp \
	api/internal/BamRandomAccessController_p.cpp \
	api/internal/BamRandomAccessController_p.h \
	api/internal/BamReader_p.cpp \
	api/internal/BamReader_p.h \
	api/internal/ILocalIODevice_p.h \
	api/internal/ILocalIODevice_p.cpp \
	api/internal/BgzfStream_p.cpp \
	api/internal/BgzfStream_p.h \
	api/internal/BamFile_p.h \
	api/internal/BamFile_p.cpp \
	api/internal/SamFormatParser_p.h \
	api/internal/SamFormatParser_p.cpp \
	api/internal/SamHeaderValidator_p.h \
	api/internal/SamHeaderValidator_p.cpp \
	api/internal/BamFtp_p.h \
	api/internal/BamFtp_p.cpp \
	api/internal/BamHttp_p.h \
	api/internal/BamHttp_p.cpp \
	api/internal/BamIndexFactory_p.h \
	api/internal/BamIndexFactory_p.cpp \
	api/internal/BamMultiMerger_p.h \
	api/internal/BamPipe_p.h \
	api/internal/BamPipe_p.cpp \
	api/internal/BamStandardIndex_p.cpp \
	api/internal/BamStandardIndex_p.h \
	api/internal/BamToolsIndex_p.h \
	api/internal/BamToolsIndex_p.cpp \
	api/internal/SamHeaderVersion_p.h \
	api/internal/BamMultiReader_p.h \
	api/internal/BamMultiReader_p.cpp \
	shared/bamtools_global.h \
	utils/bamtools_utilities.h \
	utils/bamtools_utilities.cpp  \
	utils/bamtools_fasta.cpp \
	utils/bamtools_fasta.h \
	utils/bamtools_filter_engine.h \
	utils/bamtools_filter_properties.h \
	utils/bamtools_filter_ruleparser.h \
	utils/bamtools_options.h \
	utils/bamtools_options.cpp \
	utils/bamtools_pileup_engine.h \
	utils/bamtools_pileup_engine.cpp \
	utils/bamtools_variant.h \
	utils/utils_global.h