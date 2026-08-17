#
# Currently, we have only two tests to complete.  The first is the solution
# of the SDPA format problem theta1.dat-s.  The makefile will run CSDP on 
# this problem and store the output in theta1.out.  The second test is the
# computation of the Lovasz theta number of the graph g50.  
#
# Once the tests have run, examine the .out files and compare them with
# theta1.correct and g50.correct.  The optimal values should agree to 
# six digits or more, and the DIMACS errors should all be smaller than
# 1.0e-6.  
#

.PHONY: all clean

all: theta1.out g50.out 
 
#
# The programs are run out of the build tree, so with SHARED=yes they have to
# find libsdp there rather than in $(libdir).  Exported rather than set on each
# recipe, so that adding a platform, or a test, does not mean editing every
# other line.
#
export LD_LIBRARY_PATH := ../lib
export DYLD_LIBRARY_PATH := ../lib
export PATH := ../lib:$(PATH)

#
# Test the solver on theta1.dat-s
#

theta1.out: 
	../solver/csdp theta1.dat-s >theta1.out

#
# Test theta on the g50 graph.
#

g50.out: 
	../theta/theta g50 >g50.out

clean:
	rm -f *.out
