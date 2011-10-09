# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = lobSTR 

lobSTR_CPPFLAGS = -I/usr/include/bamtools  -L/lib

AM_CPPFLAGS = $(FFTW_CFLAGS)

LDADD = $(FFTW_LIBS) -lbamtools

lobSTR_SOURCES = main.cpp runtime_parameters.cpp runtime_parameters.h MSReadRecord.h IFileReader.h FastqFileReader.h FastqFileReader.cpp FastaFileReader.h FastaFileReader.cpp TextFileReader.h TextFileReader.cpp common.h common.cpp TabFileWriter.h TabFileWriter.cpp TextFileWriter.h TextFileWriter.cpp STRDetector.h STRDetector.cpp MainLobeDetection.h MainLobeDetection.cpp TukeyWindowGenerator.h TukeyWindowGenerator.cpp FFT_four_nuc_vectors.h FFT_four_nuc_vectors.cpp MicroSatelliteDetection.h MicroSatelliteDetection.cpp HammingWindowGenerator.h HammingWindowGenerator.cpp CachedWindowGenerator.h CachedWindowGenerator.cpp FFT_nuc_vectors.h FFT_nuc_vectors.cpp SamFileWriter.cpp SamFileWriter.h MultithreadData.h MultithreadData.cpp ISatellite.h Genotyper.h Genotyper.cpp AlignmentFileReader.cpp AlignmentFileReader.h EntropyDetection.h EntropyDetection.cpp BWAReadAligner.h BWAReadAligner.cpp bwtaln.h bwtaln.c bwt.h bwt.c bamlite.h bamlite.c bntseq.h   bwtsw2.h stdaln.h stdaln.c bwtgap.h bwtgap.c bwaseqio.c bwtindex.c bwtio.c utils.h utils.c simple_dp.h simple_dp.c kstring.c kstring.h bwtsw2_core.h bwtsw2_core.c bwtsw2_main.h bwtsw2_main.c bwtsw2_aux.h bwtsw2_aux.c bwtsw2_chain.h bwtsw2_chain.c bwt_lite.c bwt_lite.h bwape.c bwase.c bwase.h bwt_gen.c bwt_gen.h QSufSort.c SQufSort.h bntseq.c bwtmisc.c bwtmisc.h is.c cs2nt.c

