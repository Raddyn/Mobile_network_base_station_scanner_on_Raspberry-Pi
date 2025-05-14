import numpy as np

def normalise_signal(signal):
    """
    Normalise the signal to the range [-1, 1].

    Parameters:
        signal (np.ndarray): complex signal to be normalised.

    Returns:
        np.ndarray: The normalised signal.
    """
    # Calculate the maximum absolute value of the signal
    max_val = np.max(np.abs(signal))

    # Avoid division by zero
    if max_val == 0:
        return signal

    # Normalise the signal
    normalised_signal = signal / max_val

    return normalised_signal