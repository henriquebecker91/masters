all: run_eduk.out run_ukp5.out test_ukp5.out run_per.out run_dom.out run_gardp.out test_gardp.out run_greendp.out test_greendp.out

CC=g++
# The -ftrapv makes the program abort if there's and integer overflow.
# This makes easies to debug overflow errors.
CFLAGS=-std=c++11 -Wall -Wpedantic -Wconversion -ftrapv -O3 #-g -O0 
LBOOST=-lboost_filesystem -lboost_system 
H_DEPS=ukp_common.hpp test_common.hpp
UKP5_DEPS=ukp5.hpp periodicity.hpp dominance.hpp
DEFS=-DHBM_CHECK_PERIODICITY -DHBM_XOR_SWAP

run_eduk.out: run_eduk.cpp eduk.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_eduk.cpp $(CFLAGS) $(DEFS)

run_ukp5.out: run_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ run_ukp5.cpp $(CFLAGS) $(LBOOST) $(DEFS)

test_ukp5.out: test_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ test_ukp5.cpp $(CFLAGS) $(LBOOST) $(DEFS)

run_per.out: run_per.cpp periodicity.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_per.cpp $(CFLAGS) $(DEFS)

run_dom.out: run_dom.cpp dominance.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_dom.cpp $(CFLAGS) $(DEFS)

run_gardp.out: run_gardp.cpp gardp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_gardp.cpp $(CFLAGS) $(DEFS)

test_gardp.out: test_gardp.cpp gardp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_gardp.cpp $(CFLAGS) $(DEFS)

run_greendp.out: run_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_greendp.cpp $(CFLAGS) $(DEFS)

test_greendp.out: test_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_greendp.cpp $(CFLAGS) $(DEFS)

.PHONY: clean

clean:
	rm *.out

