# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                        theta executable                              #
########################################################################

# Name of the executables compiled in this directory.
# We want it to be installed in the 'bin' directory
bin_PROGRAMS = theta complement rand_graph graphtoprob

# List all source files for each executable
theta_SOURCES = theta.c
complement_SOURCES = complement.c
rand_graph_SOURCES = rand_graph.c
graphtoprob_SOURCES = graphtoprob.c

# List all additionally required libraries
theta_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)
complement_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)
rand_graph_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)
graphtoprob_LDADD = ../lib/libsdp.la $(CSDPLIB_LFLAGS) $(OPENMP_CFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../include $(CSDPLIB_CFLAGS) $(OPENMP_CFLAGS)
