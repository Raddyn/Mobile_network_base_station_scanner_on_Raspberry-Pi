import argparse
import os
import sys
import time

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from core.capture import capture_samples
from core.cell_detect import lte_cell_scan
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description="LTE Cell Scanner")
    parser.add_argument(
        "-o",
        "--open",
        type=str,
        required=False,
        help="Path to the captured waveform file",
    )
    parser.add_argument(
        "-S", "--save", type=str, required=False, help="Path to save the waveform file"
    )
    parser.add_argument(
        "-f",
        "--frequency",
        type=float,
        required=True,
        help="Frequency of the captured waveform",
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
        default=3,
        help="Number of iterations used to scan the current frequency, outputs most common NID1 and NID2",
    )

    args = parser.parse_args()

    NID_2 = []
    NID_1 = []
    SSS_flag = False
    print("==== LTE Cell Scanner ====")
    print("==== Scan parameters: ===")
    print(f"{'Frequency:':<20}{args.frequency / 1e6:.2f} MHz")
    print(f"{'Sample Rate:':<20}{args.sample_rate / 1e6:.2f} MS/s")
    print(f"{'Capture Duration:':<20}{args.time:.2f} seconds")
    print(f"{'Debug Mode:':<20}{'Enabled' if args.debug else 'Disabled'}")
    print(f"{'Number of Scans:':<20}{args.num_of_scans} scans")
    print(f"{'Capture Time:':<20}{time.strftime('%d-%m-%Y %H:%M:%S', time.localtime())}")
    print("==========================")

    if args.open is None:
        for i in range(args.num_of_scans):
            if i != 0:
                time.sleep(1)  # Small delay to avoid overloading the device
            waveform = capture_samples(
                f_capture=args.frequency,
                sample_rate=int(args.sample_rate),
                num_samples=int(args.time * args.sample_rate),
            )
            if waveform is None:
                print(f"{'Error:':<20} No samples captured.")
                sys.exit()
            # Scan the waveform, show debug on the last scan if enabled
            if i == args.num_of_scans - 1:
                nid2, nid1 = lte_cell_scan(
                    waveform, sample_rate=args.sample_rate, debug=args.debug
                )
                NID_2.append(nid2)
                NID_1.append(nid1)
            else:
                nid2, nid1 = lte_cell_scan(waveform, sample_rate=args.sample_rate)
                NID_2.append(nid2)
                NID_1.append(nid1)
        most_common_nid2, count_nid2 = Counter(NID_2).most_common(1)[0]
        most_common_nid1, count_nid1 = Counter(NID_1).most_common(1)[0]
        if count_nid2 < args.num_of_scans / 2:
            print("Error: No valid PSS detected")
            sys.exit()
        if count_nid1 < args.num_of_scans / 2:
            SSS_flag = True

        print("NID_2:", most_common_nid2)
        if SSS_flag:
            print("Warning: Failed to detect SSS")
        else:
            print("NID_1:", most_common_nid1)
            print("==========================")
            print("Cell ID:", most_common_nid1 * 3 + most_common_nid2)

    else:
        # If the user provided a file, load the waveform from the file
        if not os.path.isfile(args.open):
            print(f"Error: File {args.open} does not exist.")
            sys.exit()
        else:
            waveform = lte_cell_scan.load_waveform(args.open)
            if waveform is None:
                print(f"Error: File {args.open} is not a valid waveform file.")
                sys.exit()
            # Scan the waveform, show debug on the last scan if enabled
            NID_2, NID_1 = lte_cell_scan(
                waveform, sample_rate=args.sample_rate, debug=args.debug
            )
            print("NID_2:", NID_2)
            print("NID_1:", NID_1)
            print("==========================")
            print("Cell ID:", NID_1 * 3 + NID_2)

    # =====Check if save directory exists====================================================================
    if args.save:
        # Save the waveform to a file
        if not os.path.isdir(os.path.dirname(args.save)):
            print(f"Error: Directory {os.path.dirname(args.save)} does not exist.")
            sys.exit()
        else:
            lte_cell_scan.save_waveform(waveform, args.save)
    sys.exit()


if __name__ == "__main__":
    main()
