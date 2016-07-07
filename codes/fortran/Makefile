# The '_f_' is used to not clash with the c++ binaries.
all: run_f_mtu.out

FC=gfortran
FFLAGS=-O3 -std=f2008 #-ftrapv #-g -Og # -Wall -Wpedantic -Wconversion

run_f_mtu.out: run_mtu.f
	$(FC) -o $@ $< $(FFLAGS)

.PHONY: clean

clean:
	rm *.out
