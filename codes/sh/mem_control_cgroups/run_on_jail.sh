#!/bin/bash

user=`id -un`
cgcreate -a "$user" -t "$user" -g "memory:test/$$"

# 2 * 1024**3 == 2GiB
mem_limit=$((2 * 1024**3))
echo $mem_limit > "/sys/fs/cgroup/memory/test/$$/memory.limit_in_bytes"

cgexec -g memory:test run_ukp5.out test

