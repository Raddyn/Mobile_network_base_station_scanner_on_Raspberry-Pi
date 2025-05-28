#!/bin/bash

# Zjistí složku, kde leží tento skript
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Odkaz na python kód o úroveň níže
PYTHON_SCRIPT="$SCRIPT_DIR/cli/main.py" 


# Usage message
usage() {
    echo "Usage: $0 -f <frequency> [-o <open_file>] [-S <save_file>] [-s <sample_rate>] [-T <time>] [-d <debug>] [-n <num_of_scans>]"
    exit 1
}

# Default values
sample_rate="1.92e6"
capture_time="0.02"
debug=""
num_of_scans="3"
open_file=""
save_file=""
frequency=""

# Parse options
while getopts ":o:S:f:s:T:dn:" opt; do
  case $opt in
    o) open_file="$OPTARG" ;;
    S) save_file="$OPTARG" ;;
    f) frequency="$OPTARG" ;;
    s) sample_rate="$OPTARG" ;;
    T) capture_time="$OPTARG" ;;
    d) debug="--debug" ;;
    n) num_of_scans="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

if [ -z "$frequency" ]; then
    echo "Error: Frequency (-f) is required."
    usage
fi

cmd="python3 \"$PYTHON_SCRIPT\" -f $frequency -s $sample_rate -T $capture_time -n $num_of_scans"
[ -n "$open_file" ] && cmd+=" -o $open_file"
[ -n "$save_file" ] && cmd+=" -S $save_file"
[ -n "$debug" ] && cmd+=" $debug"

eval $cmd
if [ $? -ne 0 ]; then
    echo "Error: Command failed."
    exit 1
fi

