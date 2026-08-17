#
# This Makefile can be used to automatically build the entire package.  
#
# If you make changes in the Makefile or code under any subdirectory, you can
# rebuild the system with "make clean" followed by "make all".
#
#
# Build settings are in Makefile.inc.  Every variable there can be set in the
# environment or on the make command line; see INSTALL.
#
include Makefile.inc

.PHONY: all clean install unitTest

#
# On most systems, this should handle everything.
#
all:
	$(MAKE) -C lib all
	$(MAKE) -C solver all
	$(MAKE) -C theta all
	$(MAKE) -C example all

#
# Perform a unitTest
#

unitTest: all
	$(MAKE) -C test all

csdp.pc: csdp.pc.in
	sed -e 's|@prefix@|$(prefix)|' -e 's|@exec_prefix@|$(exec_prefix)|' \
	    -e 's|@libdir@|$(libdir)|' -e 's|@includedir@|$(includedir)|' \
	    -e 's|@VERSION@|$(CSDP_VERSION)|' -e 's|@BLAS_LIBS@|$(BLAS_LIBS)|' \
	    -e 's|@LIBS@|$(LIBS)|' -e 's|@OPENMP_LIBS@|$(OPENMP_LIBS)|' \
	    -e 's|@PC_CPPFLAGS@|$(CSDP_PC_CPPFLAGS)|' \
	    csdp.pc.in > $@

#
# Install the executables in $(bindir), and the library, headers and
# pkg-config file under $(libdir) and $(includedir).  install is among the
# .PHONY targets above so that the INSTALL file does not satisfy it on case
# insensitive filesystems.
#
install: all csdp.pc
	$(INSTALL) -d "$(DESTDIR)$(bindir)"
	$(INSTALL_PROGRAM) solver/csdp "$(DESTDIR)$(bindir)/csdp"
	$(INSTALL_PROGRAM) theta/theta "$(DESTDIR)$(bindir)/$(TOOL_PREFIX)theta"
	$(INSTALL_PROGRAM) theta/graphtoprob "$(DESTDIR)$(bindir)/$(TOOL_PREFIX)graphtoprob"
	$(INSTALL_PROGRAM) theta/complement "$(DESTDIR)$(bindir)/$(TOOL_PREFIX)complement"
	$(INSTALL_PROGRAM) theta/rand_graph "$(DESTDIR)$(bindir)/$(TOOL_PREFIX)rand_graph"
	$(INSTALL) -d "$(DESTDIR)$(libdir)"
	$(INSTALL_DATA) lib/libsdp.a "$(DESTDIR)$(libdir)/libsdp.a"
ifeq ($(SHARED),yes)
ifeq ($(CSDP_IMPLIB),)
	$(INSTALL_PROGRAM) lib/$(CSDP_SHLIB) "$(DESTDIR)$(libdir)/$(CSDP_SHLIB)"
	ln -sf $(CSDP_SHLIB) "$(DESTDIR)$(libdir)/$(CSDP_SONAME)"
	ln -sf $(CSDP_SHLIB) "$(DESTDIR)$(libdir)/$(CSDP_SHLIB_DEV)"
else
	$(INSTALL_PROGRAM) lib/$(CSDP_SHLIB) "$(DESTDIR)$(bindir)/$(CSDP_SHLIB)"
	$(INSTALL_DATA) lib/$(CSDP_IMPLIB) "$(DESTDIR)$(libdir)/$(CSDP_IMPLIB)"
endif
endif
	$(INSTALL) -d "$(DESTDIR)$(includedir)/csdp"
	for h in include/*.h; do \
	  $(INSTALL_DATA) $$h "$(DESTDIR)$(includedir)/csdp" || exit 1; \
	done
	$(INSTALL) -d "$(DESTDIR)$(pkgconfigdir)"
	$(INSTALL_DATA) csdp.pc "$(DESTDIR)$(pkgconfigdir)/csdp.pc"

#
# Clean out all of the directories.
# 

clean:
	rm -f csdp.pc
	$(MAKE) -C lib clean
	$(MAKE) -C solver clean
	$(MAKE) -C theta clean
	$(MAKE) -C test clean
	$(MAKE) -C example clean
