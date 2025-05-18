import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

data = sio.loadmat('data25.mat')

iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.squeeze() + 1j * qWave.squeeze()

# decimate the waveform
decimation_factor = 1
waveform = sig.decimate(waveform, decimation_factor)

# resample 
# waveform = sig.resample(waveform, len(waveform) * 16)

N = 128


window = np.blackman(N)
spectrum = np.fft.fftshift(np.fft.fft(waveform[:N]* window))
freqs = np.fft.fftshift(np.fft.fftfreq(N, d=1/int(1.92e6//decimation_factor)))

power = 20 * np.log10(np.abs(spectrum))
plt.figure(figsize=(10, 5))
plt.plot(freqs/1e6, power)
plt.xlabel("Frekvence [MHz]")
plt.ylabel("Síla signálu [dB]")
plt.title("Spektrum LTE signálu")
plt.grid()

plt.figure(figsize=(10, 5))
plt.specgram(waveform, NFFT=N, Fs=1.92e6/decimation_factor, noverlap=N//2, scale='dB')
plt.show()