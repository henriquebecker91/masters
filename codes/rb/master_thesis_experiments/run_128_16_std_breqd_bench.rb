#!/usr/bin/ruby

require 'batch_experiment'
require 'batch_experiment/sample_extractors'

# COMMon between COMMands
comm_comm = {
  pattern:    'INST_FILE',
  extractor:  BatchExperiment::UKP5Extractor,
}

# Commands and prefixes used.
#
# The codes will output their commit version, bit in general will be commit:
# commit 7826a0b88ade40a67879f321f40ae6df0378184d
# Date:   Mon Sep 12 22:45:36 2016 -0300
# Commit message: "Adding notes on the flaws of the paper testbed."
#
# The code was compiled using the respective Makefiles and the following
# GCC configuration:
# Using built-in specs.
# COLLECT_GCC=gcc
# COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-pc-linux-gnu/6.1.1/lto-wrapper
# Target: x86_64-pc-linux-gnu
# Configured with: /build/gcc/src/gcc/configure --prefix=/usr --libdir=/usr/lib --libexecdir=/usr/lib --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=https://bugs.archlinux.org/ --enable-languages=c,c++,ada,fortran,go,lto,objc,obj-c++ --enable-shared --enable-threads=posix --enable-libmpx --with-system-zlib --with-isl --enable-__cxa_atexit --disable-libunwind-exceptions --enable-clocale=gnu --disable-libstdcxx-pch --disable-libssp --enable-gnu-unique-object --enable-linker-build-id --enable-lto --enable-plugin --enable-install-libiberty --with-linker-hash-style=gnu --enable-gnu-indirect-function --disable-multilib --disable-werror --enable-checking=release
# Thread model: posix
# gcc version 6.1.1 20160802 (GCC)
eduk = {
  command:    'pyasukpt -nobb INST_FILE',
  prefix:     'eduk',
  pattern:    'INST_FILE',
  extractor:  BatchExperiment::PyaExtractor,
}
eduk2 = {
  command:    'pyasukpt INST_FILE',
  prefix:     'eduk2',
  pattern:    'INST_FILE',
  extractor:  BatchExperiment::PyaExtractor,
}
ukp5 = {
  command:    'run_ukp5.out INST_FILE',
  prefix:     'ukp5',
}.merge(comm_comm)
ukp5_sbw = { # this one has a specific commit: 490ed589f0d8
  command:    'run_ukp5_sbw.out INST_FILE',
  prefix:     'ukp5_sbw',
}.merge(comm_comm)
cpp_mtu1 = {
  command:    'run_mtu1.out INST_FILE',
  prefix:     'cpp-mtu1',
}.merge(comm_comm)
cpp_mtu2 = {
  command:    'run_mtu2.out INST_FILE',
  prefix:     'cpp-mtu2',
}.merge(comm_comm)
mgreendp1 = {
  command:    'run_mgreendp1.out INST_FILE',
  prefix:     'mgreendp1',
}.merge(comm_comm)
# The mgreendp has a specific commit too (the one below)
# commit f69e6052da01f01949ed403afe7526347495523c
# Author: Henrique Becker <henriquebecker91@gmail.com>
# Date:   Thu Sep 22 21:09:10 2016 -0300
# Added proof of concept for multiplying the doubles before truncating them and treating them as integers.
mgreendp = {
  command:    'run_mgreendp.out INST_FILE',
  prefix:     'mgreendp',
}.merge(comm_comm)

solvers = [eduk, eduk2, ukp5, ukp5_sbw, cpp_mtu1, cpp_mtu2, mgreendp, mgreendp1]

batch_info = {
  cpus_available: [1],  # only using one isolated core
  post_timeout:   10,   # ten second to killall -9
  cwd:            '/home/henrique/AreaDeTrabalho/128_16_std_breqd_benchmark/',
  output_dir:     '/home/henrique/AreaDeTrabalho/experiment_128_16_std_breqd/all_output/',
  timeout:        600, # ten minutes
}

file_list = []
10.times do | e |
  10.times do | s |
    file_list << "128_16_std_breqd-n#{2**(11+e)}-s#{s}.ukp"
  end
end

exp_conf = {
  csvfname:     '128_16_std_breqd_all.csv',
  qt_runs:      1,
  comms_order:  :by_file,
}

BatchExperiment.experiment(solvers, batch_info, exp_conf, file_list)

