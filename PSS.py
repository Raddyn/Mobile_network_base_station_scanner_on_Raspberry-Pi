import numpy as np
from matplotlib import pyplot as plt
import scipy.io as sio
from scipy import signal



samples_in_ofdm = 2457600
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


data = np.load('pss_sequences.npy')
print(data)


data = sio.loadmat('data2.mat')
print(data.keys())
# read  'iWave', 'qWave' from data and make a complex waveform
iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.flatten() + 1j * qWave.flatten()





# plt.figure()
# #subplots
# plt.subplot(2,1,1)
# plt.plot((np.fft.ifft(pss_sequences[0])))
# plt.subplot(2,1,2)
# plt.plot(signal.resample(np.fft.ifft(pss_sequences[0]),245760))
# plt.show()


# plt.figure()
# plt.plot(waveform)
# plt.xlabel("Sample")
# plt.ylabel("Amplitude")
# plt.show()

r_s0 = signal.resample(pss_sequences[0],int(sample_rate*62/16129))
r_s1 = signal.resample(pss_sequences[1],int(sample_rate*62/16129))
r_s2 = signal.resample(pss_sequences[2],int(sample_rate*62/16129))

# r_w = signal.resample(waveform,int(sample_rate*len(waveform)/3.84e6))
r_w = waveform
plt.figure()
plt.subplot(2,1,1)
plt.plot(r_w)
plt.subplot(2,1,2)
plt.plot(waveform)
plt.show()

corr0 = max(abs(np.correlate(np.fft.ifft(r_w),np.conjugate(r_s0))))
corr1 = max(abs(np.correlate(np.fft.ifft(r_w),np.conjugate(r_s1))))
corr2 = max(abs(np.correlate(np.fft.ifft(r_w),np.conjugate(r_s2))))

print('------NID2_0------')
print(corr0)
print(corr1)
print(corr2)


