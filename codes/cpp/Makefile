all: run_eduk.out run_ukp5.out test_eduk.out test_ukp5.out run_per.out test_per.out run_bb.out test_bb.out

CC=g++
# The boost linkage is only needed if the macro PROFILE is defined...
# The best would be only use it if the macro is defined. But this would need
# some control.
CFLAGS=-std=c++11 -Wall -Wpedantic -Wconversion -O3 -lboost_filesystem -lboost_system
H_DEPS=test_common.hpp ukp_common.hpp
O_DEPS=test_common.o ukp_common.o

ukp_common.o: ukp_common.cpp ukp_common.hpp
	$(CC) -c -o $@ ukp_common.cpp $(CFLAGS)

test_common.o: test_common.cpp test_common.hpp
	$(CC) -c -o $@ test_common.cpp $(CFLAGS)

run_eduk.out: run_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ run_eduk.cpp eduk.cpp $(O_DEPS) $(CFLAGS)

test_eduk.out: test_eduk.cpp eduk.cpp eduk.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ test_eduk.cpp eduk.cpp $(O_DEPS) $(CFLAGS)

run_ukp5.out: run_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ run_ukp5.cpp ukp5.cpp $(O_DEPS) $(CFLAGS)

test_ukp5.out: test_ukp5.cpp ukp5.cpp ukp5.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ test_ukp5.cpp ukp5.cpp $(O_DEPS) $(CFLAGS)

run_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ run_per.cpp periodicity.cpp $(O_DEPS) $(CFLAGS)

test_per.out: periodicity.cpp periodicity.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ test_per.cpp periodicity.cpp $(O_DEPS) $(CFLAGS)

run_bb.out: bb.cpp bb.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ run_bb.cpp bb.cpp $(O_DEPS) $(CFLAGS)

test_bb.out: bb.cpp bb.hpp $(H_DEPS) $(O_DEPS)
	$(CC) -o $@ test_bb.cpp bb.cpp $(O_DEPS) $(CFLAGS)

.PHONY: clean

clean:
	rm *.out
	rm *.o

