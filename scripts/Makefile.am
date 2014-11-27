# Copyright (C) 2014 Assaf Gordon <agordon@wi.mit.edu>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



## Scripts listed here will be installed in the "bin" directory
## (e.g. "/usr/local/bin") and will be available on the $PATH.
## Make sure their name isn't arbitrary (which might conflict with other scripts).
##
bin_SCRIPTS =

##
## These scripts will be installed in the "/usr/share/APPNAME/scripts" directory.
## Available to users after installation, but not on the "$PATH".
## Suitible for support scripts that's not commonly used.
##
scriptdir = $(pkgdatadir)/scripts

script_SCRIPTS = \
	GetSTRInfo.py \
	lobstr_index.py \
	lobSTR_test_run.sh \
	lobSTR_download_ref_hg19.sh

EXTRA_DIST = $(script_SCRIPTS) \
	$(bin_SCRIPTS)

