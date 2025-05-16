import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

data = sio.loadmat('data1.mat')

iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.squeeze() + 1j * qWave.squeeze()

N = 128


window = np.hanning(N)
spectrum = np.fft.fftshift(np.fft.fft(waveform[:N]* window))
freqs = np.fft.fftshift(np.fft.fftfreq(N, d=1/int(1.92e6)))

power = 20 * np.log10(np.abs(spectrum))

plt.plot(freqs/1e6, power)
plt.xlabel("Frekvence [MHz]")
plt.ylabel("Síla signálu [dB]")
plt.title("Spektrum LTE signálu")
plt.grid()
plt.show()