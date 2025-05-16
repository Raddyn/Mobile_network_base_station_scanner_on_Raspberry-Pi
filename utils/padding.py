import numpy as np


def pad_sync_sig(signal, target_length):
    """
    Pad the synchronization signal to the target length by adding zeros in the middle.
    The DC component is set to zero.

    Parameters:
        signal (np.ndarray): The input signal to be padded.
        target_length (int): The desired length of the output signal.

    Returns:
        padded _signal: The padded signal.
    """
    padded_signal = np.zeros(target_length, dtype=complex)

    padded_signal[0] = 0
    padded_signal[target_length - 31 : target_length] = signal[0:31]
    padded_signal[1:32] = signal[31:62]

    return padded_signal
