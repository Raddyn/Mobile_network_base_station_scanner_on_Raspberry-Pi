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
    
    # decimate the waveform
    waveform = sig.decimate(waveform, int(sample_rate // (15000*N)))
    
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
    corr = np.zeros((3, len(waveform) + len(ifft_pss_sequences[0, :]) - 1), dtype=complex)

    # pad the waveform
    waveform = np.pad(waveform, (N//2, -1+N//2), mode='constant')
    # Perform correlation
    for i in range(3):
        corr[i, :] = sig.correlate(waveform, ifft_pss_sequences[i, :], mode='same')
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
    peaks, _ = sig.find_peaks(corr[np.argmax(max_corr), :],distance = 4000 ,height=np.max(corr[np.argmax(max_corr), :]*0.5))
    if debug:
        plt.figure()
        plt.plot(np.abs(corr[np.argmax(max_corr), :]))
        plt.plot(peaks, corr[np.argmax(max_corr), peaks], "x")
        plt.title('PSS correlation NID2:{}'.format(NID_2))
        plt.xlabel('Samples')
        plt.ylabel('Magnitude')
        plt.grid()
        plt.tight_layout()
    
    print(waveform.shape)
    print(corr.shape)
    print(peaks)
    
    
    # Locate the PSS sequence in the waveform
    pss_center_in_waveform = peaks[0] - 1
    pss_waveform = waveform[pss_center_in_waveform-N//2:pss_center_in_waveform+N//2]
    
    # Equalize the waveform
    fft_pss_waveform = np.fft.fft(pss_waveform,n=N)
    fft_pss_waveform[:] = fft_pss_waveform[:]*(np.max(padded_pss_sequences[NID_2,:])/np.max(fft_pss_waveform[:]))
    eq_func = np.zeros(N, dtype=complex)
    for i in range(N):
        if padded_pss_sequences[NID_2,i] == 0:
            eq_func[i] = 0
        else:
            eq_func[i] = np.conj(padded_pss_sequences[NID_2,i]) / fft_pss_waveform[i]
            
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(fft_pss_waveform*eq_func)
    plt.title('FFT of PSS waveform')
    plt.subplot(3, 1, 2)
    plt.plot(padded_pss_sequences[NID_2,:])
    plt.title('Padded PSS sequence')
    plt.subplot(3, 1, 3)
    plt.plot(eq_func)
    plt.title('Equalization function')
    plt.tight_layout()
    
    # Locate the SSS sequences in the waveform
    
    sss_waveform = waveform[pss_center_in_waveform-(69+N//2):pss_center_in_waveform-(69+N//2)+62]
    
    
    #generate SSS sequences
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
    
    ifft_sss_sub0 = np.zeros((168, N), dtype=complex)
    ifft_sss_sub5 = np.zeros((168, N), dtype=complex)
    
    # # Perform IFFT on SSS sequences
    for i in range(168):
        ifft_sss_sub0[i, :] = np.fft.ifft(padded_sss_sub0[i, :], n=N)
        ifft_sss_sub5[i, :] = np.fft.ifft(padded_sss_sub5[i, :], n=N)

    # Perform correlation
    corr_sub0 = np.zeros((168, len(sss_waveform) + len(ifft_sss_sub0[0, :]) - 1), dtype=complex)
    corr_sub5 = np.zeros((168, len(sss_waveform) + len(ifft_sss_sub5[0, :]) - 1), dtype=complex)
    
    for i in range(168):
        corr_sub0[i, :] = sig.correlate(sss_waveform, ifft_sss_sub0[i, :], mode='full')
        corr_sub5[i, :] = sig.correlate(sss_waveform, ifft_sss_sub5[i, :], mode='full')

    # find peaks in the correlation
    max_corr_sub0 = np.zeros(168)
    max_corr_sub5 = np.zeros(168)
    for i in range(168):
        max_corr_sub0[i] = np.max(np.abs(corr_sub0[i, :]))
        max_corr_sub5[i] = np.max(np.abs(corr_sub5[i, :]))
        
    # pick larger correlation of two subframes
    if np.max(max_corr_sub0) > np.max(max_corr_sub5):
        NID_1 = np.argmax(max_corr_sub0)
    else:
        NID_1 = np.argmax(max_corr_sub5)
        




    # sss_waveform = waveform[pss_center_in_waveform-(69+N//2):pss_center_in_waveform-(69+N//2)+62]

    # plt.figure()
    # plt.subplot(2, 1, 1)
    # plt.plot(sss_waveform)
    # plt.title('SSS waveform')
    # plt.subplot(2, 1, 2)
    # plt.plot(ifft_sss_sub0[0,:])
    # plt.title('IFFT of SSS sub0 sequence')
    # plt.tight_layout()
    
    # plt.figure()
    # plt.plot(waveform)
    
    # plt.figure()
    # plt.subplot(2, 1, 1)
    # plt.plot(sss_waveform)
    # plt.title('SSS waveform')
    # plt.subplot(2, 1, 2)
    # plt.plot(ifft_sss_sub0[0,:])
    # plt.plot(ifft_sss_sub5[0,:])
    # plt.title('IFFT of SSS sub0 and sub5')
    
    
    
    
    
    
    
    
    
    
    
    
    
    # plt.figure()
    # plt.subplot(2, 1, 1)
    # plt.plot(waveform[pss_center_in_waveform-N//2 - 2*N:pss_center_in_waveform+N])
    # plt.title('Waveform')
    # plt.subplot(2, 1, 2)
    # plt.plot(ifft_sss_sub0[0,:])
    # plt.title('IFFT of SSS sub0')

    
    # # Perform correlation
    # corr_sub0 = np.zeros((168, len(waveform) + len(ifft_sss_sub0[0, :]) - 1), dtype=complex)
    # corr_sub5 = np.zeros((168, len(waveform) + len(ifft_sss_sub5[0, :]) - 1), dtype=complex)
    
    # for i in range(168):
    #     corr_sub0[i, :] = sig.correlate(waveform, ifft_sss_sub0[i, :], mode='full')
    #     corr_sub5[i, :] = sig.correlate(waveform, ifft_sss_sub5[i, :], mode='full')

    
    # # find peaks in the correlation
    # max_corr_sub0 = np.zeros(168)
    # max_corr_sub5 = np.zeros(168)
    # for i in range(168):
    #     max_corr_sub0[i] = np.max(np.abs(corr_sub0[i, :]))
    #     max_corr_sub5[i] = np.max(np.abs(corr_sub5[i, :]))
    # NID_1_sub0 = np.argmax(max_corr_sub0)
    # NID_1_sub5 = np.argmax(max_corr_sub5)
    
    # NID_1 = None
    # if NID_1_sub0 == NID_1_sub5:
    #     NID_1 = NID_1_sub0
    # else:
    #     NID_1 = -1
    
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
            
    # # find peaks
    # peaks_sub0, _ = sig.find_peaks(corr_sub0[NID_1_sub0, :],distance = 4000)
    # peaks_sub5, _ = sig.find_peaks(corr_sub5[NID_1_sub5, :],distance = 4000)
    # if debug:
    #     plt.figure()
    #     plt.plot(np.abs(corr_sub0[NID_1_sub0, :]))
    #     plt.plot(peaks_sub0, corr_sub0[NID_1_sub0, peaks_sub0], "x")
    #     plt.title('SSS subframe0 NID1:{}'.format(NID_1_sub0))
    #     plt.xlabel('Samples')
    #     plt.ylabel('Magnitude')
    #     plt.grid()
    #     plt.tight_layout()
        
    #     plt.figure()
    #     plt.plot(np.abs(corr_sub5[NID_1_sub5, :]))
    #     plt.plot(peaks_sub5, corr_sub5[NID_1_sub5, peaks_sub5], "x")
    #     plt.title('SSS subframe5 NID1:{}'.format(NID_1_sub5))
    #     plt.xlabel('Samples')
    #     plt.ylabel('Magnitude')
    #     plt.grid()
    #     plt.tight_layout()
    #     plt.show()
    
    if debug:
        plt.show()
    return NID_2, NID_1



# Test script
if __name__ == "__main__":

    # data = sio.loadmat('data1_20Mhz.mat')
    # iWave = data['iWave']
    # qWave = data['qWave']
    # waveform = iWave.squeeze() + 1j * qWave.squeeze()

    waveform = np.load('LTE_cell_192_128.npy')

    # Load the captured waveform
    NID_2, NID_1 = (lte_cell_scan(waveform,sample_rate=1.92e6, debug=True))
    print("----- Cell IDentification -----")
    print("NCellID:", NID_1*3 + NID_2)