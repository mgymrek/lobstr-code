# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

DATA_FILES = \
	data/lobstr.allelepriors.hg19.min50.tab \
	data/trf.ucsc.hg19.filtered.bed \
	data/trf.ucsc.hg18.filtered.bed \
	data/strinfo.hg18.tab \
	data/strinfo.hg19.tab \
	data/vcf_positions_to_exclude_hg18.tab \
	data/vcf_positions_to_exclude_hg19.tab


MODEL_FILES = \
	models/illumina_v2.0.3.stuttermodel \
	models/illumina_v2.0.3.stutterproblem \
	models/illumina_v2.0.3.stepmodel


WWW_FILES = \
	www/contact-us.html \
	www/documentation.html \
	www/download.html \
	www/faq.html \
	www/file-formats.html \
	www/index.html \
	www/install.html \
	www/news.html \
	www/usage.html


EXTRA_DIST = reconf configure \
	$(WWW_FILES) \
	$(MODEL_FILES) \
	$(DATA_FILES)

SUBDIRS = m4 src scripts
