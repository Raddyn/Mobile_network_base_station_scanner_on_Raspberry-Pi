import scipy.io as sio
import numpy as np
import ltesearch as lte

data = sio.loadmat('data0.mat')

iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.squeeze() + 1j * qWave.squeeze()

pss_sequences = np.load('pss_sequences.npy')
# waveform = np.load('data.npy')

print(lte.find_pss(waveform, pss_sequences, 128, plot=True))
