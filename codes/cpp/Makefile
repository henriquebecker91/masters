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
# (all means "all except debug/profile build" in truth).

all: bin/ukp2sukp.out bin/sukp2ukp.out bin/run_eduk.out bin/run_per.out bin/run_dom.out bin/ukp5 bin/gardp bin/greendp bin/mtu1 bin/mtu2 bin/cutstock bin/ordered_step_off bin/terminating_step_off bin/gurobi bin/cplex # bin/babayev 

CC=g++
# The -ftrapv makes the program abort if there's and integer overflow.
# This makes easies to debug overflow errors, but only works if the overflow
# is from signed numbers.
SFLAGS=-std=c++11 -Ilib -Wall -Wpedantic -Wconversion -ftrapv -O3 #-g -O0 #
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
H_DEPS=lib/ukp_common.hpp lib/test_common.hpp lib/type_name.hpp lib/workarounds.hpp | bin
DEFS=-DHBM_XOR_SWAP #-DHBM_PROFILE

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

bin/run_babayev.out: mains/run_babayev.cpp lib/babayev.hpp lib/wrapper.hpp lib/mtu.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_babayev.cpp $(CFLAGS) $(DEFS)

bin/test_babayev.out: mains/test_babayev.cpp lib/babayev.hpp lib/wrapper.hpp lib/mtu.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_babayev.cpp $(CFLAGS) $(DEFS)

bin/babayev: bin/run_babayev.out bin/test_babayev.out
	touch $@

bin/bench_babayev_hard.out: mains/bench_babayev_hard.cpp lib/chung.hpp lib/mtu.hpp lib/ukp5.hpp $(H_DEPS)
	$(CC) -o $@ mains/bench_babayev_hard.cpp $(CFLAGS) $(DEFS) $(LBOOST) -DHBM_PROFILE

bin/run_ordered_step_off.out: mains/run_ordered_step_off.cpp lib/stepoff.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_ordered_step_off.cpp $(CFLAGS) $(DEFS)

bin/test_ordered_step_off.out: mains/test_ordered_step_off.cpp lib/stepoff.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_ordered_step_off.cpp $(CFLAGS) $(DEFS)

bin/ordered_step_off: bin/run_ordered_step_off.out bin/test_ordered_step_off.out
	touch $@

bin/run_terminating_step_off.out: mains/run_terminating_step_off.cpp lib/stepoff.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/run_terminating_step_off.cpp $(CFLAGS) $(DEFS)

bin/test_terminating_step_off.out: mains/test_terminating_step_off.cpp lib/stepoff.hpp lib/wrapper.hpp $(H_DEPS)
	$(CC) -o $@ mains/test_terminating_step_off.cpp $(CFLAGS) $(DEFS)

bin/terminating_step_off: bin/run_terminating_step_off.out bin/test_terminating_step_off.out
	touch $@

GUROBI_HOME   =/opt/gurobi801/linux64
GUROBI_FLAGS  =-L${GUROBI_HOME}/lib -I${GUROBI_HOME}/include -lgurobi80 -lgurobi_g++5.2

bin/run_gurobi.out: mains/run_gurobi.cpp lib/gurobi_ukp_model.hpp $(H_DEPS)
	$(CC) -o $@ $< $(CFLAGS) $(GUROBI_FLAGS)

bin/test_gurobi.out: mains/test_gurobi.cpp lib/gurobi_ukp_model.hpp $(H_DEPS)
	$(CC) -o $@ $< $(CFLAGS) $(GUROBI_FLAGS)

bin/gurobi: bin/run_gurobi.out bin/test_gurobi.out
	touch $@

# Cutting-stock with multiple different knapsack solvers

# Options used by Leonardo Moura for compiling with CPlex
CPLEX_EFLAGS  = -m64 -fPIC -fno-strict-aliasing -fexceptions -DIL_STD #-DNDEBUG 

# Folders where CPlex is stored
SYSTEM        = x86-64_linux
LIBFORMAT     = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio128/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio128/concert

# linking cplex
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEX_LFLAGS  = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -ldl
# includes for cplex
CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEX_IFLAGS  = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

CPLEX_EFLAGS += $(CPLEX_LFLAGS) $(CPLEX_IFLAGS)

# KNAPSACK_SOLVER macro possible values
CPLEX=1
UKP5_FP=2
UKP5_FP_NS=3
UKP5_INT=4
UKP5_INT_NS=5
MTU1=6
MTU2=7
MGREENDP=8
MGREENDP1=9
# ORDSO = ORDered Step-Off
ORDSO_FP=10
ORDSO_FP_NS=11
ORDSO_INT=12
ORDSO_INT_NS=13
# TERSO = TERminating Step-Off
TERSO_FP=14
TERSO_INT=15

COMMON_DEPS = lib/ukp_common.hpp lib/type_name.hpp lib/workarounds.hpp | bin

# ukp5 specific
# REDEFINE UKP5_DEFS TO REMOVE HBM_XOR_SWAP (not used in fp builds)
# HBM_CHECK_PERIODICITY and HBM_INIT_G_BY_CHUNKS are also not very useful
# in this case and can be disabled
UKP5_DEPS = lib/ukp5.hpp lib/periodicity.hpp lib/dominance.hpp $(COMMON_DEPS)
UKP5_FP_FLAGS  = $(LBOOST)
UKP5_INT_FLAGS = -DHBM_XOR_SWAP $(LBOOST)

MTU_DEPS  = lib/mtu.hpp lib/wrapper.hpp $(COMMON_DEPS)
MTU_FLAGS = -DHBM_XOR_SWAP

GREENDP_DEPS  = lib/greendp.hpp lib/wrapper.hpp $(COMMON_DEPS)
GREENDP_FLAGS = -DHBM_XOR_SWAP

STEPOFF_DEPS  = lib/stepoff.hpp lib/wrapper.hpp $(COMMON_DEPS)
STEPOFF_FLAGS = -DHBM_XOR_SWAP

bin/cplex_cutstock.out: mains/cutstock.cpp $(UKP5_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(CPLEX) $(CPLEX_EFLAGS)

bin/ukp5_fp_cutstock.out: mains/cutstock.cpp $(UKP5_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(UKP5_FP) $(UKP5_FP_FLAGS) $(CPLEX_EFLAGS)

bin/ukp5_fp_ns_cutstock.out: mains/cutstock.cpp $(UKP5_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(UKP5_FP_NS) $(UKP5_FP_FLAGS) $(CPLEX_EFLAGS)

bin/ukp5_int_cutstock.out: mains/cutstock.cpp $(UKP5_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(UKP5_INT) $(UKP5_INT_FLAGS) $(CPLEX_EFLAGS)

bin/ukp5_int_ns_cutstock.out: mains/cutstock.cpp $(UKP5_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(UKP5_INT_NS) $(UKP5_INT_FLAGS) $(CPLEX_EFLAGS)

bin/mtu1_cutstock.out: mains/cutstock.cpp $(MTU_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(MTU1) $(MTU_FLAGS) $(CPLEX_EFLAGS)

bin/mgreendp_cutstock.out: mains/cutstock.cpp $(GREENDP_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(MGREENDP) $(GREENDP_FLAGS) $(CPLEX_EFLAGS)

bin/mgreendp1_cutstock.out: mains/cutstock.cpp $(GREENDP_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(MGREENDP1) $(GREENDP_FLAGS) $(CPLEX_EFLAGS)

bin/ordso_fp_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(ORDSO_FP) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/ordso_fp_ns_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(ORDSO_FP_NS) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/ordso_int_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(ORDSO_INT) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/ordso_int_ns_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(ORDSO_INT_NS) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/terso_fp_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(TERSO_FP) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/terso_int_cutstock.out: mains/cutstock.cpp $(STEPOFF_DEPS)
	$(CC) $(CFLAGS) $< -o $@ -DKNAPSACK_SOLVER=$(TERSO_INT) $(STEPOFF_FLAGS) $(CPLEX_EFLAGS)

bin/cutstock: bin/cplex_cutstock.out bin/ukp5_fp_cutstock.out bin/ukp5_fp_ns_cutstock.out bin/ukp5_int_cutstock.out bin/ukp5_int_ns_cutstock.out bin/mtu1_cutstock.out bin/mgreendp_cutstock.out bin/mgreendp1_cutstock.out bin/ordso_fp_cutstock.out bin/ordso_fp_ns_cutstock.out bin/ordso_int_cutstock.out bin/ordso_int_ns_cutstock.out bin/terso_fp_cutstock.out bin/terso_int_cutstock.out
	touch $@

bin/fast_cutstock: bin/ukp5_fp_cutstock.out bin/ukp5_fp_ns_cutstock.out bin/ukp5_int_cutstock.out bin/ukp5_int_ns_cutstock.out bin/mtu1_cutstock.out bin/ordso_fp_cutstock.out bin/ordso_fp_ns_cutstock.out bin/ordso_int_cutstock.out bin/ordso_int_ns_cutstock.out bin/terso_fp_cutstock.out bin/terso_int_cutstock.out
	touch $@

bin/run_cplex.out: mains/run_cplex.cpp lib/cplex_ukp_model.hpp $(H_DEPS)
	$(CC) $(CFLAGS) $< -o $@ $(CPLEX_EFLAGS)

bin/test_cplex.out: mains/test_cplex.cpp lib/cplex_ukp_model.hpp $(H_DEPS)
	$(CC) $(CFLAGS) $< -o $@ $(CPLEX_EFLAGS)

bin/cplex: bin/run_cplex.out bin/test_cplex.out
	touch $@

.PHONY: clean

clean:
	rm -rf ./bin

