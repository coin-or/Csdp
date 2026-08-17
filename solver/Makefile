#
#  This builds the stand alone solver.
#
include ../Makefile.inc

.PHONY: all clean

all: csdp

csdp: csdp.o 
	$(CSDP_LINK) csdp.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o csdp
#
# To clean out the directory:
#
clean:
	rm -f *.o
	rm -f csdp
