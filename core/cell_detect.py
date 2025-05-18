import numpy as np
from sequence_generators import pss_gen, sss_gen
import scipy.signal as sig
import scipy.io as sio
import matplotlib.pyplot as plt

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.normalisation import normalise_signal


def lte_cell_scan(waveform, sample_rate = int(1.92e6), debug=False):
    
    # Parameters
    N=64
    NID_1 = None
    NID_2 = None
    
    
    #generate PSS sequences
    pss = np.zeros((3, 62), dtype=complex)
    for i in range(3):
        pss[i, :] = pss_gen(i)
    print(pss.shape)
    
    # decimate the waveform
    waveform = sig.decimate(waveform, int(sample_rate // (15000*64)))
    
    # normalise the waveform
    waveform = normalise_signal(waveform)

    # pad the PSS sequences
    padded_pss_sequences = np.zeros((3, N), dtype=complex)
    for i in range(3):
        padded_pss_sequences[i, N - 31 : N] = pss[i, :31]
        padded_pss_sequences[i, 1:32] = pss[i, 31:62]
        padded_pss_sequences[i, 0] = 0
        
    print(padded_pss_sequences.shape)
    
    # Perform IFFT on PSS sequences
    ifft_pss_sequences = np.fft.ifft(padded_pss_sequences, axis=1)
    print(ifft_pss_sequences.shape)
    
    # Perform correlation
    corr = np.zeros((3, len(waveform) + len(ifft_pss_sequences[0, :]) - 1), dtype=complex)

    for i in range(3):
        corr[i, :] = sig.correlate(waveform, ifft_pss_sequences[i, :], mode='full')
        if debug:
            plt.subplot(4, 1, i + 1)
            plt.plot(np.abs(corr[i, :]))
            plt.title(f'Correlation with PSS {i}')
            plt.xlabel('Samples')
            plt.ylabel('Magnitude')
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
        plt.title('Max correlation')
        plt.xlabel('NID2')
        plt.ylabel('Magnitude')
        plt.ylim(0, 3)
        plt.grid()
        plt.tight_layout()
        
    # find peaks
    peaks, _ = sig.find_peaks(corr[np.argmax(max_corr), :],distance = 4000 ,height=np.max(corr[np.argmax(max_corr), :])-np.average(np.abs(corr[np.argmax(max_corr), :])))
    if debug:
        plt.figure()
        plt.plot(np.abs(corr[np.argmax(max_corr), :]))
        plt.plot(peaks, corr[np.argmax(max_corr), peaks], "x")
        plt.title('PSS correlation NID2:{}'.format(NID_2))
        plt.xlabel('Samples')
        plt.ylabel('Magnitude')
        plt.grid()
        plt.tight_layout()
    
    # TODO: Equalize the waveform
    
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
    
    # Perform IFFT on SSS sequences
    ifft_sss_sub0 = np.fft.ifft(padded_sss_sub0, axis=1)
    ifft_sss_sub5 = np.fft.ifft(padded_sss_sub5, axis=1)
    
    # Perform correlation
    corr_sub0 = np.zeros((168, len(waveform) + len(ifft_sss_sub0[0, :]) - 1), dtype=complex)
    corr_sub5 = np.zeros((168, len(waveform) + len(ifft_sss_sub5[0, :]) - 1), dtype=complex)
    
    for i in range(168):
        corr_sub0[i, :] = sig.correlate(waveform, ifft_sss_sub0[i, :], mode='full')
        corr_sub5[i, :] = sig.correlate(waveform, ifft_sss_sub5[i, :], mode='full')

    
    # find peaks in the correlation
    max_corr_sub0 = np.zeros(168)
    max_corr_sub5 = np.zeros(168)
    for i in range(168):
        max_corr_sub0[i] = np.max(np.abs(corr_sub0[i, :]))
        max_corr_sub5[i] = np.max(np.abs(corr_sub5[i, :]))
    NID_1_sub0 = np.argmax(max_corr_sub0)
    NID_1_sub5 = np.argmax(max_corr_sub5)
    
    NID_1 = None
    if NID_1_sub0 == NID_1_sub5:
        NID_1 = NID_1_sub0
    else:
        NID_1 = -1
    
    if debug:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.stem(max_corr_sub0)
        plt.title('Max correlation with SSS sub0')
        plt.xlabel('NID1')
        plt.ylabel('Magnitude')
        plt.ylim(0, 3)
        plt.grid()
        plt.subplot(2, 1, 2)
        plt.stem(max_corr_sub5)
        plt.title('Max correlation with SSS sub5')
        plt.xlabel('NID1')
        plt.ylabel('Magnitude')
        plt.ylim(0, 3)
        plt.grid()
        plt.tight_layout()
            
    # find peaks
    peaks_sub0, _ = sig.find_peaks(corr_sub0[NID_1_sub0, :],distance = 4000 ,height=np.max(corr_sub0[NID_1_sub0, :])-np.average(np.abs(corr_sub0[NID_1_sub0, :])))
    peaks_sub5, _ = sig.find_peaks(corr_sub5[NID_1_sub5, :],distance = 4000 ,height=np.max(corr_sub5[NID_1_sub5, :])-np.average(np.abs(corr_sub5[NID_1_sub5, :])))
    if debug:
        plt.figure()
        plt.plot(np.abs(corr_sub0[NID_1_sub0, :]))
        plt.plot(peaks_sub0, corr_sub0[NID_1_sub0, peaks_sub0], "x")
        plt.title('SSS subframe0 NID1:{}'.format(NID_1_sub0))
        plt.xlabel('Samples')
        plt.ylabel('Magnitude')
        plt.grid()
        plt.tight_layout()
        
        plt.figure()
        plt.plot(np.abs(corr_sub5[NID_1_sub5, :]))
        plt.plot(peaks_sub5, corr_sub5[NID_1_sub5, peaks_sub5], "x")
        plt.title('SSS subframe5 NID1:{}'.format(NID_1_sub5))
        plt.xlabel('Samples')
        plt.ylabel('Magnitude')
        plt.grid()
        plt.tight_layout()
        plt.show()
    return NID_2, NID_1



# Test script
if __name__ == "__main__":

    data = sio.loadmat('data1_20Mhz.mat')
    iWave = data['iWave']
    qWave = data['qWave']
    waveform = iWave.squeeze() + 1j * qWave.squeeze()

    waveform = np.load('LTE_cell_1536_1024.npy')

    # Load the captured waveform
    lte_cell_scan(waveform,sample_rate=15.36e6,debug=True)