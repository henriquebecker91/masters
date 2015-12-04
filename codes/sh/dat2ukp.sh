#!/bin/bash
# This script takes only one parameter,
# the name of a file in the following format:
#
# <number of items>
#	1	<profit of first item>	<weight of first item>
#	2	<profit of second item>	<weight of second item>
#	...	...			...
#	<n>	<profit of the n item>	<weight of the n item>
# <knapsack capacity>
#
# and converts to the following format:
#
# n: <number of itens>
# c: <knapsack capacity>
# begin data
# <profit of first item>	<weight of first item>
# <profit of second item>	<weight of second item>
# ...	...			...
# <profit of the n item>	<weight of the n item>
# end data
#

f="$1"

itemqt=`head -n 1 "$f"`
capacity=`tail -n 1 "$f"`

tmp=`mktemp`
# get only the columns
grep '^[ \t]*[0-9]\+[ \t]\+\([0-9.]\+\)[ \t]\+\([0-9.]\+\)' "$f" | sed -e 's/^[ \t]*[0-9]\+[ \t]\+\([0-9.]\+\)[ \t]\+\([0-9.]\+\)/\2\t\1/g' > "$tmp"
# Remove the first column
out="${f%.*}.ukp"

echo "n: $itemqt"   >  "$out"
echo "c: $capacity" >> "$out"
echo "begin data"   >> "$out"
cat "$tmp"          >> "$out"
echo "end data"     >> "$out"

