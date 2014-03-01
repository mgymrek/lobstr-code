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
	config/git-version-gen \
	lobSTR.rb \
	$(top_srcdir)/.version \
	$(WWW_FILES) \
	$(DATA_FILES)

SUBDIRS = m4 src scripts models
ACLOCAL_AMFLAGS = -I m4

BUILT_SOURCES = $(top_srcdir)/.version

$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
