#!/bin/bash


if [[ "$1" == "--install" ]]; then
    chmod +x "$0"
    sudo ln -sf "$0" /usr/local/bin/lte-scan
    echo "Done, run as: lte-scan ..."
    exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/cli/main.py" 

usage() {
  echo "Usage: $0 -f <frequency> [-o <open_file>] [-S <save_file>] [-s <sample_rate>] [-T <time>] [-d <debug>] [-n <num_of_scans>] [-N <FFT_size>]"
  exit 1
}

sample_rate="1.92e6"
capture_time="0.02"
debug=""
num_of_scans="1"
FFT_size="128"
open_file=""
save_file=""
frequency=""

while getopts ":o:S:f:s:T:dn:N:" opt; do
  case $opt in
  o) open_file="$OPTARG" ;;
  S) save_file="$OPTARG" ;;
  f) frequency="$OPTARG" ;;
  s) sample_rate="$OPTARG" ;;
  T) capture_time="$OPTARG" ;;
  d) debug="--debug" ;;
  n) num_of_scans="$OPTARG" ;;
  N) FFT_size="$OPTARG" ;;
  \?) echo "Invalid option -$OPTARG" >&2; usage ;;
  :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done


cmd="python3 \"$PYTHON_SCRIPT\""
[ -n "$open_file" ] && cmd+=" -o \"$open_file\""
[ -n "$frequency" ] && cmd+=" -f $frequency"
[ -n "$save_file" ] && cmd+=" -S \"$save_file\""
[ -n "$sample_rate" ] && cmd+=" -s $sample_rate"
[ -n "$capture_time" ] && cmd+=" -T $capture_time"
[ -n "$debug" ] && cmd+=" $debug"
[ -n "$num_of_scans" ] && cmd+=" -n $num_of_scans"
[ -n "$FFT_size" ] && cmd+=" -N $FFT_size"

eval $cmd
if [ $? -ne 0 ]; then
  echo "Error: Command failed."
  exit 1
fi
