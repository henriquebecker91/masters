all: run_eduk.out run_per.out run_dom.out ukp5 gardp greendp

CC=g++
# The -ftrapv makes the program abort if there's and integer overflow.
# This makes easies to debug overflow errors.
CFLAGS=-std=c++11 -Wall -Wpedantic -Wconversion -ftrapv -O3 #-g -O0 #
LBOOST=-lboost_filesystem -lboost_system 
H_DEPS=ukp_common.hpp test_common.hpp
UKP5_DEPS=ukp5.hpp periodicity.hpp dominance.hpp
DEFS=-DHBM_XOR_SWAP
UKP5_DEFS=$(DEFS) -DHBM_CHECK_PERIODICITY -DHBM_INIT_G_BY_CHUNKS

run_eduk.out: run_eduk.cpp eduk.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_eduk.cpp $(CFLAGS) $(DEFS)

run_per.out: run_per.cpp periodicity.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_per.cpp $(CFLAGS) $(DEFS)

run_dom.out: run_dom.cpp dominance.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_dom.cpp $(CFLAGS) $(DEFS)

run_ukp5.out: run_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ run_ukp5.cpp $(CFLAGS) $(LBOOST) $(UKP5_DEFS)

test_ukp5.out: test_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ test_ukp5.cpp $(CFLAGS) $(LBOOST) $(UKP5_DEFS)

ukp5: run_ukp5.out test_ukp5.out
	touch $@

run_gardp.out: run_gardp.cpp gardp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_gardp.cpp $(CFLAGS) $(DEFS)

test_gardp.out: test_gardp.cpp gardp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_gardp.cpp $(CFLAGS) $(DEFS)

gardp: run_gardp.out test_gardp.out
	touch $@

run_greendp.out: run_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_greendp.cpp $(CFLAGS) $(DEFS)

test_greendp.out: test_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_greendp.cpp $(CFLAGS) $(DEFS)

run_mgreendp.out: run_mgreendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_mgreendp.cpp $(CFLAGS) $(DEFS)

test_mgreendp.out: test_mgreendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_mgreendp.cpp $(CFLAGS) $(DEFS)

run_greendp1.out: run_greendp1.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_greendp1.cpp $(CFLAGS) $(DEFS)

test_greendp1.out: test_greendp1.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_greendp1.cpp $(CFLAGS) $(DEFS)

run_mgreendp1.out: run_mgreendp1.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_mgreendp1.cpp $(CFLAGS) $(DEFS)

test_mgreendp1.out: test_mgreendp1.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_mgreendp1.cpp $(CFLAGS) $(DEFS)

greendp: run_greendp.out test_greendp.out run_mgreendp.out test_mgreendp.out run_greendp1.out test_greendp1.out run_mgreendp1.out test_mgreendp1.out
	touch $@

p_run_greendp.out: run_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_greendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

p_test_greendp.out: test_greendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_greendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

p_run_mgreendp.out: run_mgreendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ run_mgreendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

p_test_mgreendp.out: test_mgreendp.cpp greendp.hpp wrapper.hpp $(H_DEPS)
	$(CC) -o $@ test_mgreendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

profile_greendp: p_run_greendp.out p_test_greendp.out p_run_mgreendp.out p_test_mgreendp.out
	touch $@

.PHONY: clean

clean:
	rm -f run_eduk.out run_per.out run_dom.out run_ukp5.out test_ukp5.out ukp5 run_gardp.out test_gardp.out gardp run_greendp.out test_greendp.out run_mgreendp.out test_mgreendp.out greendp p_run_greendp.out p_test_greendp.out p_run_mgreendp.out p_test_mgreendp.out profile_greendp run_mgreendp1.out test_mgreendp1.out run_mgreendp1.out test_mgreendp1.out

