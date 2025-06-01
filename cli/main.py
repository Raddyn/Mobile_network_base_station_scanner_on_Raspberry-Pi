import argparse
import os
import sys
import time

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from core.capture import capture_samples
from core.cell_detect import lte_cell_scan
from collections import Counter
import scipy.io as sio
import numpy as np


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="LTE Cell Scanner")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-o",
        "--open",
        type=str,
        help="Path to the captured waveform file",
        required=False,
    )
    group.add_argument(
        "-f",
        "--frequency",
        type=float,
        required=False,
        nargs="+",
        help="Frequency of the captured waveform",
    )
    parser.add_argument(
        "-S", "--save", type=str, required=False, help="Path to save the waveform file"
    )
    parser.add_argument(
        "-s",
        "--sample_rate",
        type=float,
        required=False,
        default=1.92e6,
        help="Sample rate of the captured waveform",
    )
    parser.add_argument(
        "-T",
        "--time",
        type=float,
        required=False,
        default=0.02,
        help="Time to capture samples",
    )
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    parser.add_argument(
        "-n",
        "--num_of_scans",
        type=int,
        required=False,
        default=1,
        help="Number of iterations used to scan the current frequency, outputs most common NID1 and NID2",
    )
    parser.add_argument(
        "-N",
        "--FFT_size",
        type=int,
        default=128,
        help="Size of the FFT to be used in the analysis",
    )
    args = parser.parse_args()

    first_run = True

    # Main program logic
    if args.open is None:  # Determine whther to run in capture mode or file mode
        for freq in args.frequency:  # Loop through each frequency provided
            NID_2 = []
            NID_1 = []
            SSS_flag = False
            if not first_run:
                print("\n")
            if freq < 0 or freq > 6e9:
                print(f"Error: Frequency {freq} Hz is out of range (0 - 6 GHz).")
                sys.exit()

            print(f"{'Frequency:':<20}{freq / 1e6:.2f} MHz")
            print(f"{'Sample Rate:':<20}{args.sample_rate / 1e6:.2f} MS/s")
            print(f"{'Capture Duration:':<20}{args.time:.2f} seconds")
            if args.debug:
                print(f"{'Debug Mode:':<20}{'Enabled'}")
            if args.num_of_scans > 1:
                print(f"{'Number of Scans:':<20}{args.num_of_scans} scans")
            print(
                f"{'Capture Time:':<20}{time.strftime('%d-%m-%Y %H:%M:%S', time.localtime())}"
            )
            print("=====================================================")
            # Capture samples for the given frequency
            for i in range(args.num_of_scans):
                timeout = 0
                while True:
                    waveform = capture_samples(
                        f_capture=freq,
                        sample_rate=int(args.sample_rate),
                        num_samples=int(args.time * args.sample_rate),
                    )
                    if waveform is not None and len(waveform) > 0:
                        break

                    timeout += 1
                    if timeout > 3:
                        print(
                            f"{'Error:':<20} No samples captured after multiple attempts."
                        )
                        sys.exit()
                    print(
                        f"{'Warning:':<20} No samples captured, retrying... (Attempt {timeout})"
                    )
                    time.sleep(0.5)

                # Scan the waveform, show debug on the last scan if enabled
                if i == args.num_of_scans - 1:
                    nid2, nid1 = lte_cell_scan(
                        waveform,
                        sample_rate=args.sample_rate,
                        debug=args.debug,
                        N=args.FFT_size,
                    )
                    if nid2 == -1 or nid1 == -1:
                        if i > 0:
                            i -= 1
                        break
                    NID_2.append(nid2)
                    NID_1.append(nid1)
                else:
                    nid2, nid1 = lte_cell_scan(
                        waveform, sample_rate=args.sample_rate, N=args.FFT_size
                    )
                    if nid2 == -1 or nid1 == -1:
                        if i > 0:
                            i -= 1
                        break
                    NID_2.append(nid2)
                    NID_1.append(nid1)
            # Determine whether the NID_2 and NID_1 are stable across scans, if not, deem them invalid
            most_common_nid2, count_nid2 = Counter(NID_2).most_common(1)[0]
            most_common_nid1, count_nid1 = Counter(NID_1).most_common(1)[0]
            if count_nid2 < args.num_of_scans / 2:
                print(f"{'Error:':<20} No valid PSS detected")
                sys.exit()
            if count_nid1 < args.num_of_scans / 2:
                SSS_flag = True

            print(f"{'NID_2:':<20} {most_common_nid2}")
            if SSS_flag:
                print(f"{'Warning:':<20} No fix on SSS")
            else:
                print(f"{'NID_1:':<20} {most_common_nid1}")
                print("=====================================================")
                print(f"{'Cell ID:':<20} {most_common_nid1 * 3 + most_common_nid2}")
            first_run = False
    else:
        # If the user provided a file, load the waveform from the file
        if not os.path.isfile(args.open):
            print(f"Error: File {args.open} does not exist.")
            sys.exit()
        else:
            print(f"{'Input File:':<20} {args.open}")
            print(f"{'Sample Rate:':<20}{args.sample_rate / 1e6:.2f} MS/s")
            print("=====================================================")
            if args.open.endswith(".mat"):
                data = sio.loadmat(args.open)
                iWave = data["iWave"]
                qWave = data["qWave"]
                waveform = iWave.squeeze() + 1j * qWave.squeeze()
            elif args.open.endswith(".npy"):
                waveform = np.load(args.open)
            else:
                print(f"Error: File {args.open} is not a valid waveform file.")
                sys.exit()
            # Scan the waveform, show debug on the last scan if enabled
            NID_2, NID_1 = lte_cell_scan(
                waveform, sample_rate=args.sample_rate, debug=args.debug
            )
            print(f"{'NID_2:':<20}{NID_2}")
            print(f"{'NID_1:':<20}{NID_1}")
            print("=====================================================")
            print(f"{'Cell ID:':<20} {NID_1 * 3 + NID_2}")

    # Save the waveform if the user provided a save path
    if args.save:
        # Save the waveform to a file
        if not os.path.isdir(os.path.dirname(args.save)):
            print(f"Error: Directory {os.path.dirname(args.save)} does not exist.")
            sys.exit()
        else:
            np.save(args.save, waveform)
    print(f"{'Waveform saved to:':<20} {args.save}") if args.save else None

    sys.exit()


if __name__ == "__main__":
    main()
