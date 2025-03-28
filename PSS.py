import numpy as np
from matplotlib import pyplot as plt
import scipy.io as sio

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

## Perform correlation
for i in range(3):
    corr[i,:] = np.correlate(padded_pss_sequences[i], waveform, mode='full')

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
plt.tight_layout()
plt.show()


