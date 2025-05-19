import argparse
import os
import sys
import time

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
    parser.add_argument(
        "-n", "--num_of_scans", type=int, required=False, default=3, help="Number of iterations used to scan the current frequency, outputs most common NID1 and NID2")
    
    args = parser.parse_args()        

    #======Check if file exists and is a valid waveform file=========================================    
    if args.open is None:
        waveform = capture_samples(f_capture=args.frequency, sample_rate=int(args.sample_rate), num_samples=int(args.time*args.sample_rate))
    else:
        if not os.path.isfile(args.open):
            print(f"Error: File {args.open} does not exist.")
            sys.exit()
        else:
            if waveform is None:
                print(f"Error: File {args.open} is not a valid waveform file.")
                sys.exit()
            waveform = lte_cell_scan.load_waveform(args.open)
    #=====Check if save directory exists====================================================================    
    if args.save:
        # Save the waveform to a file
        if not os.path.isdir(os.path.dirname(args.save)):
            print(f"Error: Directory {os.path.dirname(args.save)} does not exist.")
            sys.exit()
        else:
            lte_cell_scan.save_waveform(waveform, args.save)

    #====Scan the waveform for LTE transmission========================================
    NID_2, NID_1 = lte_cell_scan(waveform,sample_rate=args.sample_rate,debug=True)
    print("==== Scan parameters ====")
    print("Frequency:", args.frequency)
    print("Time:", time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()))
    #TODO: Add a check to see if the NID_1 and NID_2 are valid,
    # and if not, print an error message and exit
    print("===== Found Cell ID =====")
    print("Cell ID:", NID_1 * 3 + NID_2)
    sys.exit()
    
if __name__ == "__main__":
    main()