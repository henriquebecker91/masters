# CHECKLIST FOR ADDING NEW METHODS TO THIS MAKEFILE:
# 1º) Add the bin/run_*.out and bin/test_*.out compilation rules (follow the same
#     fashion already used, i.e. copy-paste and change relevant filenames),
#     give preference to "run_" before "test_";
# 2º) Create a rule that compiles both the run and the test binary using
#     the name of the algorithm as the rule name. Put the group rule
#     after the run_/test_ rules.
# 2.5º) Optionally group similar algorithms on a bigger group.
# 3º) Add the group name to the "all" rule below.
# Exceptions: If the method don't solve UKP (ex.: run_per.out and run_dom.out);
# or if the method is too slow to solve the instances on <repo_root>/data/ukp/
# on a reasonable time (ex.: eduk); then don't create a test_*.out rule, nor a
# group rule. Also, if it's a profile version don't add it to the 'all' rule
# (it's 'all except debug/profile build' in truth).

all: bin/ukp2sukp.out bin/sukp2ukp.out bin/run_eduk.out bin/run_per.out bin/run_dom.out bin/ukp5 bin/gardp bin/greendp bin/mtu1 bin/mtu2

CC=g++
# The -ftrapv makes the program abort if there's and integer overflow.
# This makes easies to debug overflow errors, but only works if the overflow
# is from signed numbers.
SFLAGS=-std=c++11 -Ilib -Wall -Wpedantic -Wconversion -ftrapv -O3 #-g -Og #
HEAD=$(shell git describe --abbrev=12 --dirty --always --tags)
CFLAGS=$(SFLAGS) -DHBM_GIT_HEAD_AT_COMPILATION=\"$(HEAD)\"
LBOOST=-lboost_filesystem -lboost_system 
# In 11/05/2016, with gcc 6.1.1 and boost 1.60.0-4, the multiprecision library
# gave compiler errors (seems like new gcc checks noticed a bug on boost).
# The only way to compile with boost is using -fpermissive, when boost
# releases its next version this can be removed.
# See: https://svn.boost.org/trac/boost/ticket/12181
LBOOST += -fpermissive
# To understand the utility of the pipe below ("| bin"), see
# https://www.gnu.org/software/make/manual/make.html#Prerequisite-Types
H_DEPS=lib/ukp_common.hpp lib/test_common.hpp | bin
DEFS=-DHBM_XOR_SWAP

# ukp5 specific
UKP5_DEPS=lib/ukp5.hpp lib/periodicity.hpp lib/dominance.hpp
UKP5_DEFS=$(DEFS) -DHBM_CHECK_PERIODICITY -DHBM_INIT_G_BY_CHUNKS

bin:
	mkdir bin

bin/ukp2sukp.out: mains/ukp2sukp.cpp lib/ukp_common.hpp | bin
	$(CC) -o $@ mains/ukp2sukp.cpp $(CFLAGS)

bin/sukp2ukp.out: mains/sukp2ukp.cpp lib/ukp_common.hpp | bin
	$(CC) -o $@ mains/sukp2ukp.cpp $(CFLAGS)

bin/run_eduk.out: mains/run_eduk.cpp lib/eduk.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_eduk.cpp $(CFLAGS) $(DEFS)

bin/run_per.out: mains/run_per.cpp lib/periodicity.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_per.cpp $(CFLAGS) $(LBOOST) $(DEFS)

bin/run_dom.out: mains/run_dom.cpp lib/dominance.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_dom.cpp $(CFLAGS) $(DEFS)

bin/run_ukp5.out: mains/run_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ mains/run_ukp5.cpp $(CFLAGS) $(LBOOST) $(UKP5_DEFS)

bin/test_ukp5.out: mains/test_ukp5.cpp $(UKP5_DEPS) $(H_DEPS)
	$(CC) -o $@ mains/test_ukp5.cpp $(CFLAGS) $(LBOOST) $(UKP5_DEFS)

bin/ukp5: bin/run_ukp5.out bin/test_ukp5.out
	touch $@

bin/run_gardp.out: mains/run_gardp.cpp lib/gardp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_gardp.cpp $(CFLAGS) $(LBOOST) $(DEFS)

bin/test_gardp.out: mains/test_gardp.cpp lib/gardp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_gardp.cpp $(CFLAGS) $(LBOOST) $(DEFS)

bin/gardp: bin/run_gardp.out bin/test_gardp.out
	touch $@

bin/run_greendp.out: mains/run_greendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp.cpp $(CFLAGS) $(DEFS)

bin/test_greendp.out: mains/test_greendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_greendp.cpp $(CFLAGS) $(DEFS)

bin/run_mgreendp.out: mains/run_mgreendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp.cpp $(CFLAGS) $(DEFS)

bin/test_mgreendp.out: mains/test_mgreendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_mgreendp.cpp $(CFLAGS) $(DEFS)

bin/run_greendp1.out: mains/run_greendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp1.cpp $(CFLAGS) $(DEFS)

bin/test_greendp1.out: mains/test_greendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_greendp1.cpp $(CFLAGS) $(DEFS)

bin/run_mgreendp1.out: mains/run_mgreendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp1.cpp $(CFLAGS) $(DEFS)

bin/test_mgreendp1.out: mains/test_mgreendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_mgreendp1.cpp $(CFLAGS) $(DEFS)

bin/run_greendp2.out: mains/run_greendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp2.cpp $(CFLAGS) $(DEFS)

bin/test_greendp2.out: mains/test_greendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_greendp2.cpp $(CFLAGS) $(DEFS)

bin/run_mgreendp2.out: mains/run_mgreendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp2.cpp $(CFLAGS) $(DEFS)

bin/test_mgreendp2.out: mains/test_mgreendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_mgreendp2.cpp $(CFLAGS) $(DEFS)

bin/greendp: bin/run_greendp.out bin/test_greendp.out bin/run_mgreendp.out bin/test_mgreendp.out bin/run_greendp1.out bin/test_greendp1.out bin/run_mgreendp1.out bin/test_mgreendp1.out bin/run_greendp2.out bin/test_greendp2.out bin/run_mgreendp2.out bin/test_mgreendp2.out
	touch $@

bin/p_run_greendp.out: mains/run_greendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/p_run_mgreendp.out: mains/run_mgreendp.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/p_run_greendp1.out: mains/run_greendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp1.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/p_run_mgreendp1.out: mains/run_mgreendp1.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp1.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/p_run_greendp2.out: mains/run_greendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_greendp2.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/p_run_mgreendp2.out: mains/run_mgreendp2.cpp lib/greendp.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mgreendp2.cpp $(CFLAGS) $(DEFS) -DHBM_PROFILE

bin/profile_greendp: bin/p_run_greendp.out bin/p_run_mgreendp.out bin/p_run_greendp1.out bin/p_run_mgreendp1.out bin/p_run_greendp2.out bin/p_run_mgreendp2.out
	touch $@

bin/run_mtu1.out: mains/run_mtu1.cpp lib/mtu.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mtu1.cpp $(CFLAGS) $(DEFS)

bin/test_mtu1.out: mains/test_mtu1.cpp lib/mtu.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_mtu1.cpp $(CFLAGS) $(DEFS)

bin/mtu1: bin/run_mtu1.out bin/test_mtu1.out
	touch $@

bin/run_mtu2.out: mains/run_mtu2.cpp lib/mtu.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_mtu2.cpp $(CFLAGS) $(DEFS)

bin/test_mtu2.out: mains/test_mtu2.cpp lib/mtu.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_mtu2.cpp $(CFLAGS) $(DEFS)

bin/mtu2: bin/run_mtu2.out bin/test_mtu2.out
	touch $@

.PHONY: clean

clean:
	rm -rf ./bin

