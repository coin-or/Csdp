# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                         csdp executable                              #
########################################################################

# Name of the executable compiled in this directory.
# We want it to be installed in the 'bin' directory
bin_PROGRAMS = csdp

# List all source files for this executable
csdp_SOURCES = csdp.c

# List all additionally required libraries
csdp_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../include $(CSDPLIB_CFLAGS) $(OPENMP_CFLAGS)
