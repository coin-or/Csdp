# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                        theta executable                              #
########################################################################

# Name of the executables compiled in this directory.
noinst_PROGRAMS = example

# List all source files for each executable
example_SOURCES = example.c

# List all additionally required libraries
example_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../include $(CSDPLIB_CFLAGS) $(OPENMP_CFLAGS)

#
# To clean up the directory.
#
clean-local:
	rm -f prob.dat-s
	rm -f prob.sol

