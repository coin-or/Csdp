#
#  Build the library.
#
include ../Makefile.inc

.PHONY: all clean

OBJS = readprob.o sdp.o op_o.o psd_feas.o op_a.o op_at.o  Fnorm.o calc_pobj.o calc_dobj.o trace_prod.o zero_mat.o mat_mult.o sym_mat.o copy_mat.o addscaledmat.o  user_exit.o make_i.o allocmat.o initsoln.o initparams.o add_mat.o writesol.o readsol.o easysdp.o writeprob.o solvesys.o makefill.o mat_multsp.o norms.o linesearch.o matvec.o chol.o  qreig.o tweakgap.o freeprob.o packed.o  sortentries.o

ifeq ($(SHARED),yes)
all: libsdp.a $(CSDP_SHLIB)
else
all: libsdp.a
endif

libsdp.a: $(OBJS)
	$(AR) crs libsdp.a $(OBJS)

$(CSDP_SHLIB): $(OBJS)
	$(CSDP_LINK) $(CSDP_SHLIB_FLAGS) -o $(CSDP_SHLIB) $(OBJS) \
		$(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm
	$(CSDP_SHLIB_LINKS)
#
# To clean things up.
#
clean:
	rm -f *.o
	rm -f libsdp.a
	rm -f libsdp.so* libsdp*.dylib libsdp*.dll libsdp.dll.a
