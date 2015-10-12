all: run_eduk.out run_ukp5.out test_eduk.out test_ukp5.out run_per.out test_per.out

CC=g++
CFLAGS=-std=c++11 -Wall -Wconversion -O3 
H_DEPS=test_common.hpp ukp_common.hpp
O_DEPS=test_common.o ukp_common.o
DEFS=-DTWO_MULT_COMP

ukp_common.o: ukp_common.cpp ukp_common.hpp
	$(CC) $(DEFS) $(CFLAGS) -c -o $@ ukp_common.cpp

test_common.o: test_common.cpp test_common.hpp
	$(CC) $(DEFS) $(CFLAGS) -c -o $@ test_common.cpp

run_eduk.out: run_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ run_eduk.cpp eduk.cpp $(O_DEPS)

test_eduk.out: test_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ test_eduk.cpp eduk.cpp $(O_DEPS)

run_ukp5.out: run_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ run_ukp5.cpp ukp5.cpp $(O_DEPS)

test_ukp5.out: test_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ test_ukp5.cpp ukp5.cpp $(O_DEPS)

run_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ run_per.cpp periodicity.cpp $(O_DEPS)

test_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(O_DEPS)
	$(CC) $(DEFS) $(CFLAGS) -o $@ test_per.cpp periodicity.cpp $(O_DEPS)

.PHONY: clean

clean:
	rm *.out
	rm *.o

