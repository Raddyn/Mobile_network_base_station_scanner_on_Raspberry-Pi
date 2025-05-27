import argparse
import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from core.capture import capture_samples


def main():
    parser = argparse.ArgumentParser(description="Spectrograam capture")
    parser.add_argument(
        "-f", "--frequency", type=float, required=True, help="Frequency of the captured waveform")
    parser.add_argument(
        "-s", "--sample_rate", type=float, required=False, default=1.92e6, help="Sample rate of the captured waveform")
    parser.add_argument(
        "-T", "--time", type=float, required=False, default=0.02, help="Time to capture samples")
    
    
    args = parser.parse_args()
    # Capture samples
    waveform = capture_samples(f_capture=args.frequency, sample_rate=int(args.sample_rate), num_samples=int(args.time * args.sample_rate))
    if waveform is None:
        print("Error: No samples captured.")
        sys.exit()
    # Plot the waveform
    sample_rate = int(args.sample_rate)
    plt.figure()
    plt.specgram(waveform, Fs=sample_rate, NFFT=int(sample_rate//15000), noverlap=int((sample_rate//15000)/2))
    plt.title("Waveform Spectrogram")
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.colorbar(label="Magnitude")
    plt.grid()
    #plot spectrum
    plt.figure(figsize=(10, 6))
    plt.stem(np.fft.fftfreq(len(waveform), 1/sample_rate), np.abs(np.fft.fft(waveform)), use_line_collection=True)
    plt.show()
    
if __name__ == "__main__":
    main()
    print('Done!')
    sys.exit()