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
local dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd $dir
cd ../cpp/

# Number of iterations
local N=$1

# Before running this script you may want to isolate the CPU that will
# be used (to get more precise times). I (Henrique Becker) do that on
# linux by passing "isolcpus=3" to kernel by the bootloader. I also
# disable the hypertreading (can be done by the BIOS, not in my case,
# or by /sys/devices/system/cpu/cpuX/online, my case). Note that
# the CPU numbering starts at zero, despite what htop tells you.
local CPU_FOR_USE=$2

# relative path from where the binaries are to the instances
local PATH_I=$3
# relative path from where the binaries are to where to put the results
local PATH_R='../../data/results/tmp/'

local TIME_FORMAT="Ext_time: %e\nExt_mem: %M\n"
local CSV_HEADER="Filename;UKP5 Internal Time (seconds);UKP5 External Time (seconds);UKP5 Max Memory Use (Kb);UKP5 opt;Pyasukp Internal Time (seconds);Pyasukp External Time (seconds);Pyasukp Max Memory Use (Kb);PYAsUKP opt;"

local files="$4"
for f in $files
do
	local ukp5_res="$PATH_R${f}5.res"
	local pya_res="$PATH_R${f}_pya.res"
	for ((i=0; i < N; ++i))
	do
		taskset -c "$CPU_FOR_USE" time -f "$TIME_FORMAT" ./run_ukp5.out "$PATH_I$f" >> "$ukp5_res" 2>&1
		taskset -c "$CPU_FOR_USE" time -f "$TIME_FORMAT" pyasukpt -nobb -src "$PATH_I$f" >> "$pya_res" 2>&1
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
	local tmp=`mktemp 'filename.tmp.XXX'`
	yes "$f" | head -n `wc -l "$ukp5_res" | cut -d\  -f1` > $tmp
	paste -d\; "$tmp" "$ukp5_res" "$pya_res" >> "$PATH_R${f}.all"
	rm "$tmp" "$ukp5_res" "$pya_res"
done

echo "If there were no errors the results should be at the $PATH_R folder (from the sh folder)"
}

function itest_performance {
	echo $1 $2 $3
	local name=$1
	local patt=$2
	local core=$3
	local folder="../../data/ukp/too_big_instances/instances_after_mail/$name/"
	local files=`ls ${folder} | grep "$patt"`
	test_performance 1 $core "$folder" "${files[@]}"
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

#sc_folder='../../data/ukp/sc/'
#sc_files=($(ls ${sc_folder}))
#test_performance 1 3 "$sc_folder" sc_files[@]

#wcd_folder='../../data/ukp/wcd/'
#wcd_files=($(ls ${wcd_folder} | grep 'hi_.*\.ukp'))
#test_performance 1 2 "$wcd_folder" wcd_files[@]

#echo "ss + ss2 + sc"
#itest_performance 'ss' 'ss_.*\.ukp' 1 &
##itest_performance 'ss2' 'ss2_.*\.ukp' 1 &
##itest_performance 'sc' 'sc_.*-[0-6]-.*\.ukp' 1 &
##test_performance 'sc' 'sc_.*-\([7-9]\|1[0-3]\)-.*\.ukp' 2 &
##itest_performance 'sc' 'sc_.*-1[4-9]-.*\.ukp' 3 &

##itest_performance 'wcd' 'hi_n5000-.*.ukp' 2 &
##itest_performance 'wcd' 'hi_n10000-.*.ukp' 2 &
itest_performance 'wcd' 'hi_n20000-.*.ukp' 2 &
#itest_performance 'wcd' 'hi_n50000-.*.ukp' 2
##itest_performance 'postponed_per' 'nsds2_n20000.*\.ukp' 3 &
#itest_performance 'postponed_per' 'nsds2_n50000.*\.ukp' 3
#itest_performance 'wcd' 'hi2_n5000-.*.ukp' 1
#itest_performance 'wcd' 'hi2_n10000-.*.ukp' 1
#itest_performance 'wcd' 'hi2_n20000-.*.ukp' 1
#itest_performance 'wcd' 'hi2_n50000-.*.ukp' 2
##itest_performance 'saw' 'saw_n50000wmin5000.*\.ukp' 3 &
#itest_performance 'saw' 'saw_n50000wmin10000.*\.ukp' 1
#itest_performance 'saw' 'saw_n10000w.*\.ukp' 1
#itest_performance 'saw' 'saw_n100000w.*\.ukp' 3
#itest_performance 'postponed_per' 'nsds2_n50000wmin50000-\([7-9]\|79\|8[0-9]\|9[0-9]\)-.*\.ukp' 3

