#!/bin/sh

# This script assumes that it's inside the folder /codes/sh/ of a directory
# structure like the one in https://github.com/henriquebecker91/masters.
# Then it proceeds to make the working directory the directory where 
# the c++ code is inside (line bellow).
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd $dir
cd ../cpp/

# Before running this script you may want to isolate the CPU that will
# be used (to get more precise times). I (Henrique Becker) do that on
# linux by passing "isolcpus=3" to kernel by the bootloader. I also
# disable the hypertreading (can be done by the BIOS, not in my case,
# or by /sys/devices/system/cpu/cpuX/online, my case). Note that
# the CPU numbering starts at zero, despite what htop tells you.
CPU_FOR_USE=3

# The objective of this script is simply to test the performance of
# some different ways of sorting by efficiency.
for d in TWO_MULT_COMP INT_EFF FP_EFF RATIONAL_EFF
do
	make --always-make "DEFS=-D$d" test_per.out run_per.out
	taskset -c $CPU_FOR_USE ./test_per.out > "$d.txt"

	# If you downloaded this from the repository, comment the for
	# loop bellow, this instances aren't in the repo because they
	# are too big.
	for i in ../../data/ukp/big_instances/*
	do
		taskset -c $CPU_FOR_USE ./run_per.out "$i" >> "$d.txt"
	done

	# Remove some unnecessary text of the output, we are only interested
	# in the times. It's expected that the obtained values will
	# differ from the expect ones (we are not calculating the opt
	# but a periodicity bound).
	grep -v -E '(Obtai|Expe|The|opt|y_opt).*' "$d.txt" > "../../data/results/different_sorts/$d.txt"
	rm "$d.txt" 
done

echo 'If there were no errors, you can now check the updated files in the /data/results/different_sorts/ folder'

