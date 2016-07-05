#!/bin/bash

new_file="$2.sukp"
ukp2sukp.out < "$2" > "$new_file"
run_f_mtu.out "$1" "$new_file"
rm "$new_file"

