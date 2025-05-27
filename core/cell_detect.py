import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import scipy.io as sio

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def lte_cell_scan(waveform, sample_rate=int(1.92e6), debug=False):
    # Parameters
    N = 128
    NID_1 = None
    NID_2 = None

    # show the waveform using spectrogram
    if debug:	
        plt.figure()
        plt.specgram(waveform, Fs=sample_rate, NFFT=int(sample_rate//15000), noverlap=int((sample_rate//15000)/2))
        plt.title("Waveform Spectrogram")
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (Hz)")
        plt.colorbar(label="Magnitude")
        plt.grid()
    
    
    # generate PSS sequences
    pss = np.zeros((3, 62), dtype=complex)
    for i in range(3):
        pss[i, :] = pss_gen(i)

    # decimate the waveform
    if sample_rate > 1.92e6:
        waveform = sig.decimate(waveform, int(sample_rate // (15000 * N)))

    # normalise the waveform
    waveform = normalise_signal(waveform)

    # pad the PSS sequences
    padded_pss_sequences = np.zeros((3, N), dtype=complex)
    for i in range(3):
        padded_pss_sequences[i, N - 31 : N] = pss[i, :31]
        padded_pss_sequences[i, 1:32] = pss[i, 31:62]
        padded_pss_sequences[i, 0] = 0

    # Perform IFFT on PSS sequences
    ifft_pss_sequences = np.fft.ifft(padded_pss_sequences, axis=1)

    # Perform correlation
    corr = np.zeros(
        (3, len(waveform) + len(ifft_pss_sequences[0, :]) - 1), dtype=complex
    )

    # pad the waveform
    waveform = np.pad(waveform, (N // 2, -1 + N // 2), mode="constant")
    if debug:
        plt.figure()
    # Perform correlation
    for i in range(3):
        corr[i, :] = sig.correlate(waveform, ifft_pss_sequences[i, :], mode="same")
        if debug:
            plt.subplot(4, 1, i + 1)
            plt.plot(np.abs(corr[i, :]))
            plt.title(f"Correlation with PSS {i}")
            plt.xlabel("Samples")
            plt.ylabel("Magnitude")
            plt.ylim(0, 3)
            plt.grid()

    # find peaks in the correlation
    max_corr = np.zeros(3)
    for i in range(3):
        max_corr[i] = np.max(np.abs(corr[i, :]))
    NID_2 = np.argmax(max_corr)

    if debug:
        plt.subplot(4, 1, 4)
        plt.stem(max_corr)
        plt.title("Max correlation")
        plt.xlabel("NID2")
        plt.ylabel("Magnitude")
        plt.ylim(0, 3)
        plt.grid()
        plt.tight_layout()


    # Locate the PSS sequence in the waveform
    pss_center_in_waveform = np.argmax(np.abs(corr[NID_2, :]))


    # Locate the SSS sequences in the waveform
    sss_waveform_cp_normal = waveform[pss_center_in_waveform - ((9*(N//128)+N) + N // 2) : pss_center_in_waveform- ((9*(N//128)+N) + N // 2) + N]
    plt.figure()
    plt.plot((sss_waveform_cp_normal))
    plt.title("SSS in waveform normal CP")
    sss_waveform_cp_extended = waveform[pss_center_in_waveform - ((32*(N//128)+N) + N // 2) : pss_center_in_waveform- ((32*(N//128)+N) + N // 2) + N]
    plt.figure()
    plt.plot((sss_waveform_cp_extended))
    plt.title("SSS in waveform extended CP")
    
    



    # generate SSS sequences
    sss_sub0 = np.zeros((168, 62), dtype=complex)
    sss_sub5 = np.zeros((168, 62), dtype=complex)

    for i in range(168):
        sss_sub0[i], sss_sub5[i] = sss_gen(i, NID_2)

    # Pad the SSS sequences
    padded_sss_sub0 = np.zeros((168, N), dtype=complex)
    padded_sss_sub5 = np.zeros((168, N), dtype=complex)

    for i in range(168):
        padded_sss_sub0[i, N - 31 : N] = sss_sub0[i, :31]
        padded_sss_sub0[i, 1:32] = sss_sub0[i, 31:62]
        padded_sss_sub0[i, 0] = 0

        padded_sss_sub5[i, N - 31 : N] = sss_sub5[i, :31]
        padded_sss_sub5[i, 1:32] = sss_sub5[i, 31:62]
        padded_sss_sub5[i, 0] = 0

    # # Perform IFFT on SSS sequences
    ifft_sss_sub0 = np.zeros((168, N), dtype=complex)
    ifft_sss_sub5 = np.zeros((168, N), dtype=complex)

    for i in range(168):
        ifft_sss_sub0[i, :] = np.fft.ifft((padded_sss_sub0[i, :]), n=N)
        ifft_sss_sub5[i, :] = np.fft.ifft((padded_sss_sub5[i, :]), n=N)

    # Perform correlation
    corr_sub0 = np.zeros(
        (168, N + len(ifft_sss_sub0[0, :]) - 1), dtype=complex
    )
    corr_sub5 = np.zeros(
        (168, N + len(ifft_sss_sub5[0, :]) - 1), dtype=complex
    )

    peak_diff_sub0 = np.zeros(2)
    peak_diff_sub5 = np.zeros(2)

    max_corr_sub0 = np.zeros(168)
    max_corr_sub5 = np.zeros(168)
    
    SSS_corr_info = [None, -1, -1, None] # ["subframe", "peak_diff", "NID1", cp_type]
    for K in range(2):
        for i in range(168):
            if K == 0:
                corr_sub0[i, :] = sig.correlate(sss_waveform_cp_normal, ifft_sss_sub0[i, :], mode="full")
                corr_sub5[i, :] = sig.correlate(sss_waveform_cp_normal, ifft_sss_sub5[i, :], mode="full")
            else:
                corr_sub0[i, :] = sig.correlate(sss_waveform_cp_extended, ifft_sss_sub0[i, :], mode="full")
                corr_sub5[i, :] = sig.correlate(sss_waveform_cp_extended, ifft_sss_sub5[i, :], mode="full")
        
        # find peaks in the correlation
        for i in range(168):
            max_corr_sub0[i] = np.max(np.abs(corr_sub0[i, :]))
            max_corr_sub5[i] = np.max(np.abs(corr_sub5[i, :]))
        # Find the NID1 with the maximum correlation
        peak_diff_sub0[K] = max_corr_sub0[np.argmax(max_corr_sub0)] - np.mean(max_corr_sub0)
        peak_diff_sub5[K] = max_corr_sub5[np.argmax(max_corr_sub5)] - np.mean(max_corr_sub5)
        if debug:
            fig, axs = plt.subplots(2, 1, figsize=(8.27, 11.69 / 2))
            plt.rc('font', family='serif')

            axs[0].stem(max_corr_sub0, linefmt='k-', markerfmt='k_', basefmt='k')
            axs[0].plot(np.mean(max_corr_sub0) * np.ones(168), 'r--', label='Mean')
            axs[0].plot(np.max(max_corr_sub0) * np.ones(168), 'g--', label='Max')
            axs[0].stem(np.argmax(max_corr_sub0), max_corr_sub0[np.argmax(max_corr_sub0)], 'b', basefmt=' ', markerfmt='bx', label=f'NID_1 {np.argmax(max_corr_sub0)}')
            axs[0].legend(loc='lower right', framealpha=1)
            axs[0].set_title(f"SSS subframe: 0, Max-Mean difference = {peak_diff_sub0[K]:.2f}")
            axs[0].set_xlabel("NID1")
            axs[0].set_ylabel("Correlation coefficient")
            axs[0].grid()

            axs[1].stem(max_corr_sub5, linefmt='k-', markerfmt='k_', basefmt='k')
            axs[1].plot(np.mean(max_corr_sub5) * np.ones(168), 'r--', label='Mean')
            axs[1].plot(np.max(max_corr_sub5) * np.ones(168), 'g--', label='Max')
            axs[1].stem(np.argmax(max_corr_sub5), max_corr_sub5[np.argmax(max_corr_sub5)], 'b', basefmt=' ', markerfmt='bx', label=f'NID_1 {np.argmax(max_corr_sub5)}')
            axs[1].legend(loc='lower right', framealpha=1)
            axs[1].set_title(f"SSS subframe: 5, Max-Mean difference = {peak_diff_sub5[K]:.2f}")
            axs[1].set_xlabel("NID1")
            axs[1].set_ylabel("Correlation coefficient")
            axs[1].grid()

            cp_type = "normal" if K == 0 else "extended"
            fig.suptitle(f"Max correlation with SSS, Cyclic prefix: {cp_type}", fontweight="bold", y=0.975)
            plt.tight_layout()
        # Find the NID1 with the maximum correlation
        if SSS_corr_info[1] < peak_diff_sub0[K] or SSS_corr_info[1] < peak_diff_sub5[K]:
            SSS_corr_info[3] = "normal" if K == 0 else "extended"   
            if peak_diff_sub0[K] > peak_diff_sub5[K]:
                SSS_corr_info[0] = "sub0"
                SSS_corr_info[1] = peak_diff_sub0[K]
                SSS_corr_info[2] = np.argmax(max_corr_sub0)
            else:
                SSS_corr_info[0] = "sub5"
                SSS_corr_info[1] = peak_diff_sub5[K]
                SSS_corr_info[2] = np.argmax(max_corr_sub5)
                
    # Print the results
    if debug:
        print("SSS correlation info:")
        print("Subframe:", SSS_corr_info[0])
        print("Peak difference:", SSS_corr_info[1])
        print("NID1:", SSS_corr_info[2])
        print("Cyclic prefix type:", SSS_corr_info[3])
        print("NID2:", NID_2)
        
        print("Peak difference sub0:", peak_diff_sub0)
        print("Peak difference sub5:", peak_diff_sub5)
    NID_1 = SSS_corr_info[2]



    if debug:
        plt.show()
    return NID_2, NID_1


def pss_gen(NID2):
    """
    Generate the Primary Synchronization Signal (PSS) for a given NID2.
    This implementation is based on the LTE PSS generation algorithm as described in
    https://www.sharetechnote.com/html/Handbook_LTE_PSS.html accessed on 19.05.2025

    Parameters:
        NID2 (int): The NID2 value (0, 1, or 2).


    Returns:
        np.ndarray: The generated PSS sequence.
    """
    # -------------------------------------------------------

    if NID2 not in [0, 1, 2]:
        raise ValueError("NID2 must be 0, 1, or 2.")

    # -------------------------------------------------------
    i = NID2
    u_shift = [25, 29, 34]
    pss = np.array(
        [
            [
                np.exp(-1j * np.pi * u_shift[i] * n * (n + 1) / 63)
                if n <= 30
                else np.exp(-1j * np.pi * u_shift[i] * (n + 1) * (n + 2) / 63)
                for n in range(62)
            ]
        ]
    )
    return pss


def sss_gen(NID1, NID2):
    """
    Generate the Secondary Synchronization Signal (SSS) for given NID1 and NID2.

    This implementation is based on the LTE SSS generation algorithm as described in
    https://www.sharetechnote.com/html/Handbook_LTE_SSS.html accessed on 19.05.2025

    Parameters:
        NID1 (int): The NID1 value (0 to 167).
        NID2 (int): The NID2 value (0 to 2).

    Returns:
        tuple: Two numpy arrays representing the SSS sequences.

        sss0 (np.ndarray): The first SSS sequence.
        sss5 (np.ndarray): The second SSS sequence.
    """
    # -------------------------------------------------------
    if NID1 < 0 or NID1 > 167:
        raise ValueError("NID1 must be between 0 and 167.")
    if NID2 < 0 or NID2 > 2:
        raise ValueError("NID2 must be between 0 and 2.")
    # -------------------------------------------------------
    q_prime = np.floor(NID1 / 30)
    q = np.floor((NID1 + q_prime * (q_prime + 1) / 2) / 30)
    m_prime = NID1 + q * (q + 1) / 2
    m0 = m_prime % 31
    m1 = (m0 + np.floor(m_prime / 31) + 1) % 31
    x_s = np.zeros(31)
    x_s[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_s[i + 5] = (x_s[i + 2] + x_s[i]) % 2
    x_c = np.zeros(31)
    x_c[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_c[i + 5] = (x_c[i + 3] + x_c[i]) % 2
    s_tilda = 1 - 2 * x_s
    c_tilda = 1 - 2 * x_c
    s0_m0_even = np.zeros(31)
    s1_m1_even = np.zeros(31)
    for n in range(31):
        s0_m0_even[n] = s_tilda[int((n + m0) % 31)]
        s1_m1_even[n] = s_tilda[int((n + m1) % 31)]
    c0_even = np.zeros(31)
    for n in range(31):
        c0_even[n] = c_tilda[int((n + NID2) % 31)]
    d_even_sub0 = s0_m0_even * c0_even
    d_even_sub5 = s1_m1_even * c0_even
    x_z = np.zeros(31)
    x_z[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_z[i + 5] = (x_z[i + 4] + x_z[i + 2] + x_z[i + 1] + x_z[i]) % 2
    z_tilda = 1 - 2 * x_z
    s1_m1_odd = np.zeros(31)
    s0_m0_odd = np.zeros(31)
    for n in range(31):
        s1_m1_odd[n] = s_tilda[int((n + m1) % 31)]
        s0_m0_odd[n] = s_tilda[int((n + m0) % 31)]
    c1_odd = np.zeros(31)
    for n in range(31):
        c1_odd[n] = c_tilda[int((n + NID2 + 3) % 31)]
    z1_m0_odd = np.zeros(31)
    z1_m1_odd = np.zeros(31)
    for n in range(31):
        z1_m0_odd[n] = z_tilda[int((n + m0 % 8) % 31)]
        z1_m1_odd[n] = z_tilda[int((n + m1 % 8) % 31)]
    d_odd_sub0 = s1_m1_odd * c1_odd * z1_m0_odd
    d_odd_sub5 = s0_m0_odd * c1_odd * z1_m1_odd
    d_sub0 = np.zeros(62)
    d_sub0[::2] = d_even_sub0
    d_sub0[1::2] = d_odd_sub0
    sss0 = d_sub0
    d_sub5 = np.zeros(62)
    d_sub5[::2] = d_even_sub5
    d_sub5[1::2] = d_odd_sub5
    sss5 = d_sub5

    return sss0, sss5

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


# Test script
if __name__ == "__main__":
    data = sio.loadmat('data2.mat')
    iWave = data['iWave']
    qWave = data['qWave']
    waveform = iWave.squeeze() + 1j * qWave.squeeze()

    waveform = np.load("LTE_cell_192_128.npy")

    # Load the captured waveform
    NID_2, NID_1 = lte_cell_scan(waveform, sample_rate=1.92e6, debug=True)
    print("----- Cell IDentification -----")
    print("NCellID:", NID_1 * 3 + NID_2)
