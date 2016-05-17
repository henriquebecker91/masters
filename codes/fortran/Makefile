# The '_f_' is used to not clash with the c++ binaries.
all: run_f_mtu2.out # run_f_mtu1.out

FC=gfortran
FFLAGS=-ftrapv -O3 #-g -Og # -Wall -Wpedantic -Wconversion

#run_f_mtu1.out: run_mtu1.f
#	$(FC) -o $@ $< $(FFLAGS)

run_f_mtu2.out: run_mtu2.f
	$(FC) -o $@ $< $(FFLAGS)

.PHONY: clean

clean:
	rm *.out
