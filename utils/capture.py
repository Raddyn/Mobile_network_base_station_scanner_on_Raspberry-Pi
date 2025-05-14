import adi
import numpy as np

def capture_samples(f_capture, sample_rate, num_samples,ip_address='ip:192.168.2.1',save=False):
    """
    Capture samples from the PlutoSDR.
    
    Parameters:
    f_capture (float): The center frequency for capturing samples.
    sample_rate (float): The sample rate for capturing samples.
    num_samples (int): The number of samples to capture.
    
    Returns:
    np.ndarray: The captured samples.
    """
    # Create a PlutoSDR object
    sdr = adi.Pluto(uri=ip_address)
    
    # Set the sample rate and center frequency
    sdr.sample_rate = int(sample_rate)
    sdr.rx_lo = int(f_capture)
    
    # Set the buffer size
    sdr.rx_buffer_size = num_samples
    
    # Capture samples
    samples = sdr.rx()
    
    if save:
        # Save to a file as a numpy array
        np.save('captured_samples.npy', samples)
    
    return samples