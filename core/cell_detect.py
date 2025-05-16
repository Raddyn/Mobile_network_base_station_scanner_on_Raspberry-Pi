import numpy as np
from sequence_generators import pss_gen, sss_gen

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.normalisation import normalise_signal
from utils.padding import pad_sync_sig


def lte_cell_scan(captured_waveform, frequency, N=128, tolerance=0.5, debug=False):

    
    # Normalize the waveform data
    captured_waveform = normalise_signal(captured_waveform)
    
    # Check wavform for PSS
    for NID2 in range(3):
        pss = pss_gen(NID2)
        
        # padd PSS for the correct N for IFFT
        pss =  pad_sync_sig(pss, N)
        
        
    #     correlation = np.correlate(captured_waveform, pss, mode='same')
    #     if np.max(np.abs(correlation)) > tolerance:
    #         print(f"Detected PSS with NID2: {NID2}")
    #         if debug:
    #             print(f"Correlation: {correlation}")
    #         break
    # else:
    #     print("No PSS detected.")
    #     return 1
        
    
    return 0



# Test script
if __name__ == "__main__":
    # Example usage
    captured_waveform = np.random.rand(1000)  # Replace with actual waveform data
    lte_cell_scan(captured_waveform)