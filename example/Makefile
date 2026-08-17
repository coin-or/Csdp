#
# Target all builds everything.
#
include ../Makefile.inc

.PHONY: all clean

all: example
#
# This builds the example code.
#
example: example.o 
	$(CSDP_LINK) example.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o example
#
# To clean up the directory.
#
clean:
	rm -f *.o
	rm -f example
	rm -f prob.dat-s
	rm -f prob.sol
