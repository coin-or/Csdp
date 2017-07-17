#
# Target all builds everything.
#
all: example
#
# This builds the example code.
#
example: example.o 
	$(CC) $(CFLAGS) example.o $(LIBS) -o example
#
# To clean up the directory.
#
clean:
	rm -f *.o
	rm -f example
	rm -f prob.dat-s
	rm -f prob.sol



