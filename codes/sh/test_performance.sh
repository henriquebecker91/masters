#!/bin/bash

# Args:
# $1 Number of iterations
# $2 CPU core to execute computation
# $3 The relative path from the binaries directory (masters/codes/cpp)
#     to the folder that contain all the instances
# $4 A list of instance filenames that are in $3
function test_performance {
# This script assumes that it's inside the folder /codes/sh/ of a directory
# structure like the one in https://github.com/henriquebecker91/masters.
# Then it proceeds to make the working directory the directory where 
# the c++ code is inside (line bellow).
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd $dir
cd ../cpp/

# Number of iterations
N=$1

# Before running this script you may want to isolate the CPU that will
# be used (to get more precise times). I (Henrique Becker) do that on
# linux by passing "isolcpus=3" to kernel by the bootloader. I also
# disable the hypertreading (can be done by the BIOS, not in my case,
# or by /sys/devices/system/cpu/cpuX/online, my case). Note that
# the CPU numbering starts at zero, despite what htop tells you.
CPU_FOR_USE=$2

# relative path from where the binaries are to the instances
PATH_I=$3
# relative path from where the binaries are to where to put the results
PATH_R='../../data/results/tmp/'

TIME_FORMAT="Ext_time: %e\nExt_mem: %M\n"
CSV_HEADER="Filename;UKP5 Internal Time (seconds);UKP5 External Time (seconds);UKP5 Max Memory Use (Kb);UKP5 opt;Pyasukp Internal Time (seconds);Pyasukp External Time (seconds);Pyasukp Max Memory Use (Kb);PYAsUKP opt;"

declare -a files=("${!4}")
for f in "${files[@]}"
do
	ukp5_res="$PATH_R${f}5.res"
	pya_res="$PATH_R${f}_pya.res"
	for ((i=0; i < N; ++i))
	do
		taskset -c "$CPU_FOR_USE" time -f "$TIME_FORMAT" ./run_ukp5.out "$PATH_I$f" >> "$ukp5_res" 2>&1
		taskset -c "$CPU_FOR_USE" time -f "$TIME_FORMAT" pyasukp "$PATH_I$f" >> "$pya_res" 2>&1
	done

	grep -E '^Seconds: .*' "$ukp5_res" | cut -d ' ' -f2 > "${ukp5_res}.itime"
	grep -E '^Ext_time: .*' "$ukp5_res" | cut -d ' ' -f2 > "${ukp5_res}.etime"
	grep -E '^Ext_mem: .*' "$ukp5_res" | cut -d ' ' -f2 > "${ukp5_res}.emem"
	grep -E '^opt: .*' "$ukp5_res" | cut -d ':' -f2 | tr -d ' ' > "${ukp5_res}.opt"
	paste -d\; "${ukp5_res}.itime" "${ukp5_res}.etime" "${ukp5_res}.emem" "${ukp5_res}.opt" > "$ukp5_res"
	rm "${ukp5_res}.itime" "${ukp5_res}.etime" "${ukp5_res}.emem" "${ukp5_res}.opt"

	grep -E '^Total Time :.*' "$pya_res" | sed -e 's/Total Time :[[:blank:]]\+//' > "${pya_res}.itime"
	grep -E '^Ext_time: .*' "$pya_res" | cut -d ' ' -f2 > "${pya_res}.etime"
	grep -E '^Ext_mem: .*' "$pya_res" | cut -d ' ' -f2 > "${pya_res}.emem"
	pcregrep -M '^#The optimal value for the given capacity\n\d' "$pya_res" | tail -n 1 > "${pya_res}.opt"
	paste -d\; "${pya_res}.itime" "${pya_res}.etime" "${pya_res}.emem" "${pya_res}.opt" > "$pya_res"
	rm "${pya_res}.itime" "${pya_res}.etime" "${pya_res}.emem" "${pya_res}.opt"

	echo $CSV_HEADER > "$PATH_R${f}.all"
	tmp=`mktemp 'filename.tmp.XXX'`
	yes "$f" | head -n `wc -l "$ukp5_res" | cut -d\  -f1` > $tmp
	paste -d\; "$tmp" "$ukp5_res" "$pya_res" >> "$PATH_R${f}.all"
	rm "$tmp" "$ukp5_res" "$pya_res"
done

echo "If there were no errors the results should be at the $PATH_R folder (from the sh folder)"
}

# Almost gone insane trying to pass lists to functions in bash, see
# http://stackoverflow.com/questions/1063347/passing-arrays-as-parameters-in-bash
#pyasukp_bench_folder='../../data/ukp/'
#pyasukp_bench_files=(exnsd16.ukp exnsd18.ukp exnsd20.ukp exnsd26.ukp exnsdbis18.ukp exnsds12.ukp)
#test_performance 1000 3 "$pyasukp_bench_folder" pyasukp_bench_files[@]

#buriol_bench_folder='../../data/ukp/buriol/'
#buriol_bench_files=($(ls ${buriol_bench_folder}))
#test_performance 1 3 "$buriol_bench_folder" buriol_bench_files[@]

#myinst_bench_folder='../../data/ukp/myinst/'
#myinst_bench_files=($(ls ${myinst_bench_folder}))
#test_performance 1 3 "$myinst_bench_folder" myinst_bench_files[@]

ss_folder='../../data/ukp/ss/'
ss_files=($(ls ${ss_folder}))
test_performance 1 3 "$ss_folder" ss_files[@]

