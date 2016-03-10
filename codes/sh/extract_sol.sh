#!/bin/bash
# To be used on the codes/cpp/run_ukp5.out output, it removes the optimal
# solution items from the rest of the log putting them on a two column format
# (first weight, second profit).

grep "\(\.\./\|ix\)" "$1" | sed 's/.*w: \([0-9]\+\) p: \([0-9]\+\)/\1 \2/g'

