all: run_eduk.out test_eduk.out run_ukp5.out test_ukp5.out run_bb.out test_bb.out run_per.out test_per.out 

CC=g++
# The boost linkage is only needed if the macro PROFILE is defined...
# The best would be only use it if the macro is defined. But this would need
# some control.
CFLAGS=-std=c++11 -Wall -Wpedantic -Wconversion -O3 -lboost_filesystem -lboost_system
H_DEPS=ukp_common.hpp test_common.hpp

run_eduk.out: run_eduk.cpp eduk.hpp $(H_DEPS)
	$(CC) -o $@ run_eduk.cpp $(CFLAGS)

test_eduk.out: test_eduk.cpp eduk.hpp $(H_DEPS)
	$(CC) -o $@ test_eduk.cpp $(CFLAGS)

run_ukp5.out: run_ukp5.cpp ukp5.hpp $(H_DEPS)
	$(CC) -o $@ run_ukp5.cpp $(CFLAGS)

test_ukp5.out: test_ukp5.cpp ukp5.hpp $(H_DEPS)
	$(CC) -o $@ test_ukp5.cpp $(CFLAGS)

run_per.out: run_per.cpp periodicity.hpp $(H_DEPS)
	$(CC) -o $@ run_per.cpp $(CFLAGS)

test_per.out: test_per.cpp periodicity.hpp $(H_DEPS)
	$(CC) -o $@ test_per.cpp $(CFLAGS)

run_bb.out: run_bb.cpp bb.hpp $(H_DEPS)
	$(CC) -o $@ run_bb.cpp $(CFLAGS)

test_bb.out: test_bb.cpp bb.hpp $(H_DEPS)
	$(CC) -o $@ test_bb.cpp $(CFLAGS)

.PHONY: clean

clean:
	rm *.out

