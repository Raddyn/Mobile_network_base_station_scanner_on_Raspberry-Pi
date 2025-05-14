import numpy as np
from matplotlib import pyplot as plt
import scipy.io as sio
from utils.normalisation import normalise_signal

## Variables
N = 128


## Load pregenerated PSS sequences
pss_sequences = np.load('pss_sequences.npy')


## Load waveform data
data = sio.loadmat('data0.mat')
# data = sio.loadmat('data1.mat')
# data = sio.loadmat('data2.mat')
# data = sio.loadmat('data3.mat')


iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.squeeze() + 1j * qWave.squeeze()

# waveform = np.load('data.npy')
# waveform = normalise_signal(waveform)


## Create paddded PSS sequences for IFFT
padded_pss_sequences = np.zeros((3, N), dtype=complex)

for i in range(3):
    padded_pss_sequences[i, 97:128] = pss_sequences[i, :31]
    padded_pss_sequences[i, 1:32] = pss_sequences[i, 31:62]
    padded_pss_sequences[i, 0] = 0

## Perform IFFT on PSS sequences
for i in range(3):

    padded_pss_sequences[i] = np.fft.ifft(padded_pss_sequences[i], n=128)

## pregenerate
corr = np.zeros((3,(len(waveform) + len(padded_pss_sequences[0,:]) - 1)), dtype=complex)

#pad the waveform to the same length as the PSS sequences
padded_waveform = np.pad(waveform, (N-1,N-1), 'constant', constant_values=(0, 0))

## Perform correlation
for i in range(3):
    corr[i,:] = np.correlate(padded_pss_sequences[i], padded_waveform, mode='valid')

## Find max correlation
max_corr = np.zeros(3)
for i in range(3):
    max_corr[i] = max(abs(corr[i,:]))



print("----- PSS_detection -----")
print("Found NID2:",np.argmax(max_corr))

print("Correlation values:")
print(max_corr)

plt.figure()
plt.plot((waveform))
plt.title('Waveform')


plt.figure()
for i in range(3):
    plt.subplot(3,1,i+1)
    plt.plot(np.abs(corr[i,:]))
    plt.title('NID2 = ' + str(i))
    #set fixed y-axis limits
    plt.ylim(0,2)
plt.tight_layout()
plt.show()

#extract the PSS sequence from the highest peak
# find the index of the highest peak
max_index = np.argmax(max_corr)
# find the index of the highest peak in the correlation
max_index = np.argmax(abs(corr[max_index,:]))
# extract the PSS sequence from the correlation
pss_sequence_data = waveform[max_index-N//2:max_index+N//2]

plt.figure()
plt.subplot(3,1,1)
plt.plot(np.fft.fft(pss_sequence_data, n=128))
plt.title('Extracted PSS sequence')
plt.subplot(3,1,2)
plt.plot(np.fft.fft(padded_pss_sequences[np.argmax(max_corr),:], n=128))
plt.title('Original PSS sequence')
plt.subplot(3,1,3)
plt.plot(pss_sequences[np.argmax(max_corr),:])
plt.show()


# equalize the waveform
