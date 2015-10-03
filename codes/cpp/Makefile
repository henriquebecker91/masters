all: run_eduk.out run_ukp5.out test_eduk.out test_ukp5.out run_per.out test_per.out

CC=g++
CFLAGS=-std=c++11 -Wall -Wconversion -O3 
H_DEPS=test_common.hpp ukp_common.hpp
C_DEPS=test_common.cpp ukp_common.cpp

run_eduk.out: run_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ run_eduk.cpp eduk.cpp $(C_DEPS)

test_eduk.out: test_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ test_eduk.cpp eduk.cpp $(C_DEPS)

run_ukp5.out: run_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ run_ukp5.cpp ukp5.cpp $(C_DEPS)

test_ukp5.out: test_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ test_ukp5.cpp ukp5.cpp $(C_DEPS)

run_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ run_per.cpp periodicity.cpp $(C_DEPS)

test_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(C_DEPS)
	$(CC) $(CFLAGS) -o $@ test_per.cpp periodicity.cpp $(C_DEPS)

.PHONY: clean

clean:
	rm *.out

