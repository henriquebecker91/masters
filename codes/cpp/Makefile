all: run_eduk.out run_ukp5.out test_eduk.out test_ukp5.out

CC=g++
CFLAGS=-std=c++11 -Wall
DEPS=test_common.cpp test_common.hpp ukp_common.cpp ukp_common.hpp

run_eduk.out: run_eduk.cpp eduk.cpp eduk.hpp $(DEPS)
	$(CC) $(CFLAGS) -g -o $@ run_eduk.cpp test_common.cpp ukp_common.cpp

run_ukp5.out: run_ukp5.cpp ukp5.cpp ukp5.hpp $(DEPS)
	$(CC) $(CFLAGS) -o $@ run_ukp5.cpp ukp5.cpp test_common.cpp ukp_common.cpp

test_eduk.out: test_eduk.cpp eduk.cpp eduk.hpp $(DEPS)
	$(CC) $(CFLAGS) -o $@ test_eduk.cpp eduk.cpp test_common.cpp ukp_common.cpp

test_ukp5.out: test_ukp5.cpp ukp5.cpp ukp5.hpp $(DEPS)
	$(CC) $(CFLAGS) -o $@ test_ukp5.cpp ukp5.cpp test_common.cpp ukp_common.cpp

.PHONY: clean

clean:
	rm *.out

