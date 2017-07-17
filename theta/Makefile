#
# Target all builds everything.
#
all: theta complement rand_graph graphtoprob
#
# This builds the theta number code. 
#
theta: theta.o 
	$(CC) $(CFLAGS) theta.o $(LIBS) -o theta
#
# Complement computes the complement of a graph.
#
complement: complement.o 
	$(CC) $(CFLAGS) complement.o $(LIBS) -o complement
#
# rand_graph generates a random graph.  
#
rand_graph: rand_graph.o
	$(CC) $(CFLAGS) rand_graph.o $(LIBS) -o rand_graph
#
# graphtoprob converts a file in the graph format to an SDP problem in our
# SDP format.
#
graphtoprob: graphtoprob.o 
	$(CC) $(CFLAGS) graphtoprob.o $(LIBS) -o graphtoprob
#
# To clean up the directory.
#
clean:
	rm -f *.o
	rm -f theta
	rm -f complement
	rm -f rand_graph
	rm -f graphtoprob


