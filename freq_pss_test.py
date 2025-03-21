import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy.fft import fft, ifft
from scipy import signal
import scipy.io as sio

sample_rate = 3.84e6

u_shift = [25, 29, 34]
pss_sequences = np.array([
    [
        np.exp(-1j * np.pi * u_shift[0] * n * (n + 1) / 63)
        if n <= 30 else
        np.exp(-1j * np.pi * u_shift[0] * (n + 1) * (n + 2) / 63)
        for n in range(62)
    ]
    ,
    [
        np.exp(-1j * np.pi * u_shift[1] * n * (n + 1) / 63)
        if n <= 30 else
        np.exp(-1j * np.pi * u_shift[1] * (n + 1) * (n + 2) / 63)
        for n in range(62)
    ]
    ,
    [
        np.exp(-1j * np.pi * u_shift[2] * n * (n + 1) / 63)
        if n <= 30 else
        np.exp(-1j * np.pi * u_shift[2] * (n + 1) * (n + 2) / 63)
        for n in range(62)
    ]
])



data = sio.loadmat('data1.mat')
iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.flatten() + 1j * qWave.flatten()

fft_waveform = fft(waveform)
corr0= signal.correlate(waveform, np.conjugate(pss_sequences[0]))
corr1= signal.correlate(waveform, np.conjugate(pss_sequences[1]))
corr2= signal.correlate(waveform, np.conjugate(pss_sequences[2]))
print(max(abs(corr0)),'\n',max(abs(corr1)),'\n',max(abs(corr2)))
plt.figure()
plt.plot(np.abs(corr0))
plt.plot(np.abs(corr1))
plt.plot(np.abs(corr2))
plt.legend(["PSS0", "PSS1", "PSS2"])
plt.xlabel("Sample")
plt.ylabel("Correlation")
plt.show()



plt.figure()
plt.subplot(3,1,1)
# remove lines and add a marker
plt.plot(np.real(pss_sequences[0]),np.imag(pss_sequences[0]), marker='o', color='b', linestyle='None')
plt.subplot(3,1,2)
plt.plot(np.real(pss_sequences[1]),np.imag(pss_sequences[1]), marker='o', color='r',linestyle='None')
plt.subplot(3,1,3)
plt.plot(np.real(pss_sequences[2]),np.imag(pss_sequences[2]), marker='o', color='g',linestyle='None')
plt.show()