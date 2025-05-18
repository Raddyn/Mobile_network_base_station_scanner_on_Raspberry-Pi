import numpy as np


def pad_sync_sig(signal, N):
    """
    Pad the synchronization signal to the target length by adding zeros in the middle.
    The DC component is set to zero.

    Parameters:
        signal (np.ndarray): The input signal to be padded.
        N (int): The desired length of the output signal.

    Returns:
        padded_signal: The padded signal.
    """
    padded_signal = np.zeros(N, dtype=complex)

    padded_signal[N-31:N] = signal[:31]
    padded_signal[1:32] = signal[31:62]
    padded_signal[0] = 0

    return padded_signal
