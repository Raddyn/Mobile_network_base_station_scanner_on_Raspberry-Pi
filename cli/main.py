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
        "-s", "--save", type=str, required=False, help="Path to save the waveform file")
    parser.add_argument(
        "-f", "--frequency", type=float, required=True, help="Frequency of the captured waveform")
    parser.add_argument(
        "-T", "--time", type=float, required=False, default=0.06, help="Time to capture samples")
    parser.add_argument(
        "-t", "--tolerance", type=float, required=False, default=0.5, help="Tolerance for cell detection")
    # parser.add_argument(
        # "-n", "--num_cells", type=int, required=False, default=1, help="Number of cells to detect")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug mode")
    
    args = parser.parse_args()        
    #------------------------------------------------------------------------------------------------------------------
    
    # Check file argument, if not provided, capture using Adalm-Pluto
    if args.open is None:
        waveform = capture_samples(f_capture=args.frequency, sample_rate=int(30.72e6), num_samples=int(args.time*30.72e6))
    else:
    # Load waveform from file
        if not os.path.isfile(args.open):
            print(f"Error: File {args.open} does not exist.")
            sys.exit()
        else:
            waveform = lte_cell_scan.load_waveform(args.open)
            if waveform is None:
                print(f"Error: File {args.open} is not a valid waveform file.")
                sys.exit()
    #------------------------------------------------------------------------------------------------------------------
    
    # Save the waveform to a file if the save argument is provided    
    if args.save:
        # Save the waveform to a file
        if not os.path.isdir(os.path.dirname(args.save)):
            print(f"Error: Directory {os.path.dirname(args.save)} does not exist.")
            sys.exit()
        else:
            lte_cell_scan.save_waveform(waveform, args.save)
    #------------------------------------------------------------------------------------------------------------------
    # scan the waveform for LTE transmission and determine bandwidtht in order find the center frequency in order to
    # set the correct number of samples for the IFFT

    #------------------------------------------------------------------------------------------------------------------

    # LTE Cell detection
    lte_cell_scan.detect_cells(waveform, args.frequency, args.tolerance, args.debug)

    #------------------------------------------------------------------------------------------------------------------
    

if __name__ == "__main__":
    main()