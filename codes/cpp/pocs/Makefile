# Created as part of the autodoc_code.cpp proof of concept, autodoc_code.cpp
# can't be compiled without passing the HBM_CURR_COMMIT_AT_COMPILATION flag
# by command line.

all: autodoc_code.out autodoc_code_fail.out get_clock_precision.out uninitialized_resize.out vector_index_vs_iterator.out

CC=g++
CFLAGS=-std=c++11 -Ilib -Wall -Wpedantic -Wconversion -O3
VERSION=$(shell git describe --abbrev=12 --dirty --always --tags)
VERSION_FLAG=-DHBM_GIT_HEAD_AT_COMPILATION=\"$(VERSION)\"

autodoc_code_fail.out: autodoc_code.cpp
	$(CC) -o $@ autodoc_code.cpp $(CFLAGS)
autodoc_code.out: autodoc_code.cpp
	$(CC) -o $@ autodoc_code.cpp $(CFLAGS) $(VERSION_FLAG)
get_clock_precision.out: get_clock_precision.cpp
	$(CC) -o $@ get_clock_precision.cpp $(CFLAGS)
uninitialized_resize.out: uninitialized_resize.cpp
	$(CC) -o $@ uninitialized_resize.cpp $(CFLAGS)
vector_index_vs_iterator.out: vector_index_vs_iterator.cpp
	$(CC) -o $@ vector_index_vs_iterator.cpp $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.out

