# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = lobSTR strangers

lobSTR_CPPFLAGS = -I/usr/include/bamtools  -L/lib

AM_CPPFLAGS = $(FFTW_CFLAGS)

LDADD = $(FFTW_LIBS) -lbamtools

lobSTR_SOURCES = main.cpp runtime_parameters.cpp runtime_parameters.h MSReadRecord.h IFileReader.h FastqFileReader.h FastqFileReader.cpp FastaFileReader.h FastaFileReader.cpp TextFileReader.h TextFileReader.cpp common.h common.cpp TabFileWriter.h TabFileWriter.cpp TextFileWriter.h TextFileWriter.cpp STRDetector.h STRDetector.cpp MainLobeDetection.h MainLobeDetection.cpp TukeyWindowGenerator.h TukeyWindowGenerator.cpp FFT_four_nuc_vectors.h FFT_four_nuc_vectors.cpp MicroSatelliteDetection.h MicroSatelliteDetection.cpp HammingWindowGenerator.h HammingWindowGenerator.cpp CachedWindowGenerator.h CachedWindowGenerator.cpp FFT_nuc_vectors.h FFT_nuc_vectors.cpp GST.h GST.cpp GSTContainer.h GSTContainer.cpp SuffixTree.h SuffixTree.cpp NodeAlignment.h NodeAlignment.cpp MSRecord.h stree.h stree.cpp MSTableFileReader.h MSTableFileReader.cpp GenomeReader.h GenomeReader.cpp ReadAligner.h ReadAligner.cpp SamFileWriter.cpp SamFileWriter.h MultithreadData.h MultithreadData.cpp ISatellite.h Genotyper.h Genotyper.cpp AlignmentFileReader.cpp AlignmentFileReader.h EntropyDetection.h EntropyDetection.cpp

strangers_SOURCES = main_strangers.cpp FastaFileReader.h FastaFileReader.cpp FastqFileReader.h FastqFileReader.cpp IFileReader.h AnonymizerMultiThreadData.h AnonymizerMultiThreadData.cpp runtime_parameters.h runtime_parameters.cpp TextFileReader.h TextFileReader.cpp common.h common.cpp Anonymizer.h Anonymizer.cpp TextFileWriter.h TextFileWriter.cpp

