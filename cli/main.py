import argparse
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from core.capture import capture_samples
from core.cell_detect import lte_cell_scan


def main():
    parser = argparse.ArgumentParser(description="LTE Cell Scanner")
    parser.add_argument(
        "-o", "--open", type=str, required=False, help="Path to the captured waveform file")
    parser.add_argument(
        "-S", "--save", type=str, required=False, help="Path to save the waveform file")
    parser.add_argument(
        "-f", "--frequency", type=float, required=True, help="Frequency of the captured waveform")
    parser.add_argument(
        "-s", "--sample_rate", type=float, required=False, default=1.92e6, help="Sample rate of the captured waveform")
    parser.add_argument(
        "-T", "--time", type=float, required=False, default=0.02, help="Time to capture samples")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug mode")
    
    args = parser.parse_args()        

    #======Check if file exists and is a valid waveform file=========================================    
    if args.open is None:
        waveform = capture_samples(f_capture=args.frequency, sample_rate=int(args.sample_rate), num_samples=int(args.time*args.sample_rate))
    else:
        if not os.path.isfile(args.open):
            print(f"Error: File {args.open} does not exist.")
            sys.exit()
        else:
            waveform = lte_cell_scan.load_waveform(args.open)
            if waveform is None:
                print(f"Error: File {args.open} is not a valid waveform file.")
                sys.exit()
    #=====Check if save directory exists====================================================================    
    if args.save:
        # Save the waveform to a file
        if not os.path.isdir(os.path.dirname(args.save)):
            print(f"Error: Directory {os.path.dirname(args.save)} does not exist.")
            sys.exit()
        else:
            lte_cell_scan.save_waveform(waveform, args.save)

    #====Scan the waveform for LTE transmission========================================
    lte_cell_scan.detect_cells(waveform,sample_rate=args.sample_rate,debug=True)

if __name__ == "__main__":
    main()