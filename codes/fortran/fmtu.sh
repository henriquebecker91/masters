#!/bin/bash

echo "$2"
wd=`pwd`
sukp_file=`mktemp -p $wd`
ukp2sukp.out < "$2" > "$sukp_file"
run_f_mtu.out "$1" "$sukp_file"
rm "$sukp_file"

