all: run_eduk.out test_eduk.out run_ukp5.out test_ukp5.out run_per.out test_per.out #run_bb.out test_bb.out

CC=g++
# The boost linkage is only needed if the macro HBM_PROFILE is defined...
# The best would be only use it if the macro is defined. But this would need
# some control.
CFLAGS=-std=c++11 -Wall -Wpedantic -Wconversion -O3 -lboost_filesystem -lboost_system
H_DEPS=ukp_common.hpp test_common.hpp
DEFS=-DHBM_CHECK_PERIODICITY -DHBM_XOR_SWAP

run_eduk.out: run_eduk.cpp eduk.hpp $(H_DEPS)
	$(CC) -o $@ run_eduk.cpp $(CFLAGS) $(DEFS)

test_eduk.out: test_eduk.cpp eduk.hpp $(H_DEPS)
	$(CC) -o $@ test_eduk.cpp $(CFLAGS) $(DEFS)

run_ukp5.out: run_ukp5.cpp ukp5.hpp $(H_DEPS)
	$(CC) -o $@ run_ukp5.cpp $(CFLAGS) $(DEFS)

test_ukp5.out: test_ukp5.cpp ukp5.hpp $(H_DEPS)
	$(CC) -o $@ test_ukp5.cpp $(CFLAGS) $(DEFS)

run_per.out: run_per.cpp periodicity.hpp $(H_DEPS)
	$(CC) -o $@ run_per.cpp $(CFLAGS) $(DEFS)

test_per.out: test_per.cpp periodicity.hpp $(H_DEPS)
	$(CC) -o $@ test_per.cpp $(CFLAGS) $(DEFS)

.PHONY: clean

clean:
	rm *.out

