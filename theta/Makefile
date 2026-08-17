#
# Target all builds everything.
#
include ../Makefile.inc

.PHONY: all clean

all: theta complement rand_graph graphtoprob
#
# This builds the theta number code. 
#
theta: theta.o 
	$(CSDP_LINK) theta.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o theta
#
# Complement computes the complement of a graph.
#
complement: complement.o 
	$(CSDP_LINK) complement.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o complement
#
# rand_graph generates a random graph.  
#
rand_graph: rand_graph.o
	$(CSDP_LINK) rand_graph.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o rand_graph
#
# graphtoprob converts a file in the graph format to an SDP problem in our
# SDP format.
#
graphtoprob: graphtoprob.o 
	$(CSDP_LINK) graphtoprob.o $(CSDP_LIBS) $(BLAS_LIBS) $(LIBS) $(OPENMP_LIBS) -lm -o graphtoprob
#
# To clean up the directory.
#
clean:
	rm -f *.o
	rm -f theta
	rm -f complement
	rm -f rand_graph
	rm -f graphtoprob
