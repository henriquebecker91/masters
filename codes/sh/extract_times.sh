#!/bin/bash

# $1 Name of a file that contains one or more concatened run_ukp5.out outputs.

tmp_dir=`mktemp -d`

echo "Filename" > "$tmp_dir/filenames_$1"
grep "^.*\.ukp$" "$1" | sed 's/.*\///' >> "$tmp_dir/filenames_$1"
echo "Sort time" > "$tmp_dir/sort_time_$1"
grep "Sort time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/sort_time_$1"
echo "Dom  time" > "$tmp_dir/dom_time_$1"
grep "Dom  time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/dom_time_$1"
echo "Vect time" > "$tmp_dir/vect_time_$1"
grep "Vect time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/vect_time_$1"
echo "O(n) time" > "$tmp_dir/on_time_$1"
grep "O(n) time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/on_time_$1"
echo "pha1 time" > "$tmp_dir/pha1_time_$1"
grep "pha1 time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/pha1_time_$1"
echo "pha2 time" > "$tmp_dir/pha2_time_$1"
grep "pha2 time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/pha2_time_$1"
echo "Sum  time" > "$tmp_dir/sum_time_$1"
grep "Sum  time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/sum_time_$1"
echo "All  time" > "$tmp_dir/all_time_$1"
grep "All  time: " "$1" | cut -d: -f2 | cut -d\( -f1 | tr -d ' s' >> "$tmp_dir/all_time_$1"

paste -d\; "$tmp_dir/filenames_$1" "$tmp_dir/sort_time_$1" "$tmp_dir/dom_time_$1" "$tmp_dir/vect_time_$1" "$tmp_dir/on_time_$1" "$tmp_dir/pha1_time_$1" "$tmp_dir/pha2_time_$1" "$tmp_dir/sum_time_$1" "$tmp_dir/all_time_$1" > "times_$1"

