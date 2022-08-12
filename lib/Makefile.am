# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                               libsdp                                 #
########################################################################

# Name of the library compiled in this directory.
noinst_LTLIBRARIES = libsdp.la

# List all source files for this library
libsdp_la_SOURCES = add_mat.c \
  addscaledmat.c \
  allocmat.c \
  calc_dobj.c \
  calc_pobj.c \
  chol.c \
  copy_mat.c \
  easysdp.c \
  Fnorm.c \
  freeprob.c \
  initparams.c \
  initsoln.c \
  linesearch.c \
  makefill.c \
  make_i.c \
  mat_mult.c \
  mat_multsp.c \
  matvec.c \
  norms.c \
  op_a.c \
  op_at.c \
  op_o.c \
  packed.c \
  psd_feas.c \
  qreig.c \
  readprob.c \
  readsol.c \
  sdp.c \
  solvesys.c \
  sortentries.c \
  sym_mat.c \
  trace_prod.c \
  tweakgap.c \
  user_exit.c \
  writeprob.c \
  writesol.c \
  zero_mat.c

# This is for libtool.
libsdp_la_LDFLAGS = $(LT_LDFLAGS)

# Add library dependencies to libCoinUtils.
libsdp_la_LIBADD = $(CSDPLIB_LFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../include $(CSDPLIB_CFLAGS) $(OPENMP_CFLAGS)
