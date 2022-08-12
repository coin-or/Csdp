# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

AUTOMAKE_OPTIONS = foreign

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = lib solver theta

# We don't want to compile the test subdirectory, unless the test target is
# specified. But we need to list it as a subdirectory to make sure that it
# is included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test example
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST =

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest

########################################################################
#                  Installation of the .pc file                        #
########################################################################

#pkgconfiglibdir = $(libdir)/pkgconfig
#pkgconfiglib_DATA = coinutils.pc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean

DISTCLEANFILES =
