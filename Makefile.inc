#
# Build settings shared by the subdirectory Makefiles.  Every variable here
# can be set in the environment or on the make command line, for example
#
#	make CC=clang CFLAGS="-O3 -g" BLAS_LIBS="-lopenblas"
#
# GNU make is required.  See INSTALL for the full list of variables.
#
# CSDP_CPPFLAGS and CSDP_CFLAGS hold the flags the build itself requires.
# They are kept separate from CPPFLAGS and CFLAGS so that setting the latter
# cannot silently drop -I../include or the feature macros.
#

CC ?= cc
INSTALL ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 755
INSTALL_DATA ?= $(INSTALL) -m 644

# theta, complement, graphtoprob and rand_graph have rather generic names, so
# distributors may want to install them under a prefix, e.g. TOOL_PREFIX=csdp-
TOOL_PREFIX ?=

prefix ?= /usr/local
exec_prefix ?= $(prefix)
bindir ?= $(exec_prefix)/bin
libdir ?= $(exec_prefix)/lib
includedir ?= $(prefix)/include
pkgconfigdir ?= $(libdir)/pkgconfig

CFLAGS ?= -O2 -Wall

# Set OPENMP_CFLAGS to nothing to build CSDP without OpenMP.
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LIBS ?=

# The BLAS and LAPACK libraries, plus any Fortran run time library they need.
# On macOS this is "-framework Accelerate".
BLAS_LIBS ?= -llapack -lblas

# LIBS is anything else to link against.  It is deliberately additive and has
# no default, so that a build system which exports an empty LIBS -- FreeBSD's
# bsd.port.mk does -- cannot knock out the BLAS.

# Where the top of the tree is, relative to whichever makefile included this
# one: "./" at the top level, "../" from a subdirectory, and so on.  Spelling
# the paths this way rather than hardcoding ../ keeps them right whatever the
# depth, and keeps them relative, which -I and -L in a reproducible build want
# to be.
CSDP_TOP := $(dir $(lastword $(MAKEFILE_LIST)))

CSDP_CPPFLAGS = -I$(CSDP_TOP)include
CSDP_CFLAGS = -std=gnu99
# The build's own -L has to precede the user's LDFLAGS, for the same reason
# CSDP_CPPFLAGS precedes CPPFLAGS: -L is searched in command line order, so a
# -L naming a directory that already holds an installed libsdp would otherwise
# win over the one just built here.
CSDP_LDFLAGS = -L$(CSDP_TOP)lib
CSDP_LIBS = -lsdp

# The probes below have to see the same compiler and flags as the build, or
# they will answer for the wrong target: a packager who selects 32 bit with
# CFLAGS=-m32 rather than CC="gcc -m32" would otherwise get -DBIT64.
CSDP_PROBE = $(CC) $(CPPFLAGS) $(CFLAGS)

# Set SHARED=yes to build a shared library alongside libsdp.a.  The programs
# then link against it, since -lsdp prefers the shared library when both are
# present.
SHARED ?= no

# CSDP_VERSION is the release, and is what goes in csdp.pc.
CSDP_VERSION = 6.2.0

ifeq ($(SHARED),yes)

# The shared library is versioned independently of the release, because its
# ABI has nothing to do with the release number.  These are libtool's
# -version-info numbers, and the rules for bumping them are libtool's:
#
#	if the code changed at all                       revision + 1
#	if any interface was added, removed or changed   current + 1, revision 0
#	if any interface was added and none removed
#	  or changed                                     age + 1
#	if any interface was removed or changed          age 0
#
CSDP_CURRENT = 0
CSDP_REVISION = 0
CSDP_AGE = 0

# As libtool names them: the soname is libsdp.so.(current - age), and the file
# it points at is libsdp.so.(current - age).(age).(revision).
CSDP_SOVERSION := $(shell expr $(CSDP_CURRENT) - $(CSDP_AGE))
CSDP_LIBVERSION := $(CSDP_SOVERSION).$(CSDP_AGE).$(CSDP_REVISION)

# Ask the compiler what it is targeting rather than asking uname what we are
# running on, so that this is right when cross compiling.
CSDP_TARGET_MACROS := $(shell : | $(CSDP_PROBE) -dM -E - 2>/dev/null)

ifneq ($(findstring __APPLE__,$(CSDP_TARGET_MACROS)),)
CSDP_CFLAGS += -fPIC
CSDP_SHLIB = libsdp.$(CSDP_LIBVERSION).dylib
CSDP_SONAME = libsdp.$(CSDP_SOVERSION).dylib
CSDP_SHLIB_DEV = libsdp.dylib
# Darwin's ld does not accept 0 for these, so libtool passes current + 1.
CSDP_COMPAT_VERSION := $(shell expr $(CSDP_CURRENT) + 1)
CSDP_SHLIB_FLAGS = -dynamiclib -install_name $(libdir)/$(CSDP_SONAME) \
	-compatibility_version $(CSDP_COMPAT_VERSION) \
	-current_version $(CSDP_COMPAT_VERSION).$(CSDP_REVISION)
CSDP_SHLIB_LINKS = ln -sf $(CSDP_SHLIB) $(CSDP_SONAME) && \
	ln -sf $(CSDP_SHLIB) $(CSDP_SHLIB_DEV)

else ifneq ($(findstring _WIN32,$(CSDP_TARGET_MACROS)),)
# Windows has no sonames and no symlinks.  The DLL is named for its ABI
# version and belongs next to the programs; what goes in $(libdir) is the
# import library.  Code is always position independent, so -fPIC is not used.
CSDP_SHLIB = libsdp-$(CSDP_SOVERSION).dll
CSDP_IMPLIB = libsdp.dll.a
CSDP_SHLIB_FLAGS = -shared -Wl,--out-implib,$(CSDP_IMPLIB)
CSDP_SHLIB_LINKS = :

else
CSDP_CFLAGS += -fPIC
CSDP_SHLIB = libsdp.so.$(CSDP_LIBVERSION)
CSDP_SONAME = libsdp.so.$(CSDP_SOVERSION)
CSDP_SHLIB_DEV = libsdp.so
CSDP_SHLIB_FLAGS = -shared -Wl,-soname,$(CSDP_SONAME)
CSDP_SHLIB_LINKS = ln -sf $(CSDP_SHLIB) $(CSDP_SONAME) && \
	ln -sf $(CSDP_SHLIB) $(CSDP_SHLIB_DEV)
endif

endif

# Both of these are POSIX interfaces that Windows lacks: USESIGTERM needs
# setrlimit(), and USEGETTIME needs gettimeofday().  Probe with -include
# rather than a here-doc, so that these lines contain no "#", which make 3.81
# would take as the start of a comment.
ifneq ($(shell : | $(CSDP_PROBE) -E -include sys/resource.h - >/dev/null 2>&1 && echo yes),)
CSDP_CPPFLAGS += -DUSESIGTERM
endif

ifneq ($(shell : | $(CSDP_PROBE) -E -include sys/time.h - >/dev/null 2>&1 && echo yes),)
CSDP_CPPFLAGS += -DUSEGETTIME
endif

LONG_BYTES := $(shell echo csdp_long __SIZEOF_LONG__ | $(CSDP_PROBE) -E -P - \
	2>/dev/null | sed -n '/csdp_long/s/[^0-9]//gp')
ifeq ($(LONG_BYTES),8)
CSDP_CPPFLAGS += -DBIT64
endif

# BIT64 switches the indexing macros in index.h, which is an installed header,
# so a consumer has to be given it too: its ijtok() must do the same arithmetic
# as the library's did.  The other feature macros appear only in the .c files,
# so they stay private.
CSDP_PC_CPPFLAGS = $(filter -DBIT64,$(CSDP_CPPFLAGS))

# USEOPENMP includes <omp.h> and calls the OpenMP runtime, so define it only
# when we are actually compiling with OpenMP.
ifneq ($(strip $(OPENMP_CFLAGS)),)
CSDP_CPPFLAGS += -DUSEOPENMP -DSETNUMTHREADS
endif

CSDP_COMPILE = $(CC) $(CSDP_CPPFLAGS) $(CPPFLAGS) $(CSDP_CFLAGS) \
	$(OPENMP_CFLAGS) $(CFLAGS)
CSDP_LINK = $(CC) $(CSDP_CFLAGS) $(OPENMP_CFLAGS) $(CFLAGS) \
	$(CSDP_LDFLAGS) $(LDFLAGS)

%.o: %.c
	$(CSDP_COMPILE) -c $< -o $@
