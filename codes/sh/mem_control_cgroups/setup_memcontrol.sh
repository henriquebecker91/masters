#!/bin/bash

if [[ $# -ne 1 ]] then
	echo "usage: $0 <user to own the cgroup>"
	echo "aborted: the first parameter is required"
	exit
fi
if [[ $EUID -ne 0 ]] then
	echo "aborted: needs to be run as root"
	exit
fi

user="$1"
swapoff -a
cgcreate -a "$user" -t "$user" -g 'memory:test'

