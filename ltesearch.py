import numpy as np
import matplotlib.pyplot as plt

def find_pss(waveform, pss_sequences, N, chunk_size=int(5e-3*1.94e6), plot=False):
    '''
        # Description:
        Function to find the PSS sequences in a given waveform
        # Inputs:
            waveform: normalized np.array of complex values representing the waveform
            pss_sequences: np.array of complex values representing the PSS sequences
            N: int representing the number of samples in the IFFT
            chunk_size: int representing the number of samples processed at a time
        # Outputs:
            int representing the NID2 value found in the waveform
    '''
    ## Create paddded PSS sequences for IFFT
    padded_pss_sequences = np.zeros((3, N), dtype=complex)
    
    for i in range(3):
        padded_pss_sequences[i, N-31:N] = pss_sequences[i, :31]
        padded_pss_sequences[i, 1:32] = pss_sequences[i, 31:62]
        padded_pss_sequences[i, 0] = 0
        
    ## Perform IFFT on PSS sequences
    for i in range(3):
        padded_pss_sequences[i] = np.fft.ifft(padded_pss_sequences[i], N)
    
    ## Pregen corr array
    corr = np.zeros((3,(chunk_size + len(padded_pss_sequences[0,:]) - 1)), dtype=complex)
    
    ## Number of iq chunks
    iterations = int(len(waveform)//chunk_size)

    ## Loop through chunks and find max correlation
    for i in range(iterations):
        chunk = waveform[i*chunk_size:(i+1)*chunk_size]
        ## Perform correlation
        for i in range(3):
            corr[i,:] = np.correlate(padded_pss_sequences[i], chunk, mode='full')
        
        ## Find max correlation
        max_corr = np.zeros(3)
        for i in range(3):
            if max(abs(corr[i,:])) > max_corr[i]:
                max_corr[i] = max(abs(corr[i,:]))
    return np.argmax(max_corr)
                  
def normalize_waveform(data):
    '''
        # Description:
        Function to normalize the waveform data
        # Inputs:
            data: np.array of complex values representing the waveform
        # Outputs:
            np.array of normalized complex values representing the waveform
    '''
    return data / np.max(np.abs(data))

def find_sss(waveform, sss_sequences, N,chunksize=int(5e-3*1.94e6), plot=False):
    '''
        # Description:
        Function to find the SSS sequences in a given waveform
        # Inputs:
            waveform: normalized np.array of complex values representing the waveform
            sss_sequences: np.array of complex values representing the SSS sequences
            N: int representing the number of samples in the IFFT
        # Outputs:
            int representing the NID1 value found in the waveform
    '''
    
    ## return the index of the max correlation
    return np.argmax(max_corr)
