import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.signal as sig
from ltesearch import find_pss

def generate_sss(NID1,NID2):
    q_prime = np.floor(NID1/30)
    q = np.floor((NID1 + q_prime*(q_prime+1)/2)/30)
    m_prime = NID1 + q*(q+1)/2
    m0 = m_prime % 31
    m1 = (m0 + np.floor(m_prime/31) + 1) % 31
    x_s = np.zeros(31)
    x_s[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_s[i+5] = (x_s[i+2] + x_s[i]) % 2
    x_c = np.zeros(31)
    x_c[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_c[i+5] = (x_c[i+3] + x_c[i]) % 2
    s_tilda = 1 - 2*x_s
    c_tilda = 1 - 2*x_c
    s0_m0_even = np.zeros(31)
    s1_m1_even = np.zeros(31)
    for n in range(31):
        s0_m0_even[n] = s_tilda[int((n+m0)%31)]
        s1_m1_even[n] = s_tilda[int((n+m1)%31)]
    c0_even = np.zeros(31)
    for n in range(31):
        c0_even[n] = c_tilda[int((n+NID2)%31)]
    d_even_sub0 = s0_m0_even * c0_even
    d_even_sub5 = s1_m1_even * c0_even
    x_z = np.zeros(31)
    x_z[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_z[i+5] = (x_z[i+4] + x_z[i+2] + x_z[i+1] + x_z[i]) % 2
    z_tilda = 1 - 2*x_z
    s1_m1_odd = np.zeros(31)
    s0_m0_odd = np.zeros(31)
    for n in range(31):
        s1_m1_odd[n] = s_tilda[int((n+m1)%31)]
        s0_m0_odd[n] = s_tilda[int((n+m0)%31)]
    c1_odd = np.zeros(31)
    for n in range(31):
        c1_odd[n] = c_tilda[int((n+NID2+3)%31)]
    z1_m0_odd = np.zeros(31)
    z1_m1_odd = np.zeros(31)
    for n in range(31):
        z1_m0_odd[n] = z_tilda[int((n+m0%8)%31)]
        z1_m1_odd[n] = z_tilda[int((n+m1%8)%31)]
    d_odd_sub0 = s1_m1_odd * c1_odd * z1_m0_odd
    d_odd_sub5 = s0_m0_odd * c1_odd * z1_m1_odd
    d_sub0 = np.zeros(62)
    d_sub0[::2] = d_even_sub0
    d_sub0[1::2] = d_odd_sub0
    sss0 = d_sub0
    d_sub5 = np.zeros(62)
    d_sub5[::2] = d_even_sub5
    d_sub5[1::2] = d_odd_sub5
    sss5 = d_sub5
    return sss0, sss5

NID2 = None

data = sio.loadmat('data1.mat')
iWave = data['iWave']
qWave = data['qWave']
waveform = iWave.squeeze() + 1j * qWave.squeeze()

# waveform = np.load('LTE_cell_192_128.npy')
# waveform = np.load('LTE_cell_1536_1024.npy')
# waveform = np.load('data.npy')


# data = sio.loadmat('data1_20Mhz.mat')

# iWave = data['iWave']
# qWave = data['qWave']
# waveform = iWave.squeeze() + 1j * qWave.squeeze()

waveform = sig.decimate(waveform, 2)

N = 64
NID2 = find_pss(waveform, np.load('pss_sequences.npy'), N=N, plot=True)

print("----- PSS_detection -----")
print("Found NID2:", NID2)

sss_signals_sub0 = np.zeros((168, 62), dtype=complex)
sss_signals_sub5 = np.zeros((168, 62), dtype=complex)
for i in range(168):
    sss_signals_sub0[i], sss_signals_sub5[i] = generate_sss(i, NID2)

padded_sss_signals_sub0 = np.zeros((168, N), dtype=complex)
padded_sss_signals_sub5 = np.zeros((168, N), dtype=complex)  
for i in range(168):
    padded_sss_signals_sub0[i, N-31:N] = sss_signals_sub0[i, :31]
    padded_sss_signals_sub0[i, 1:32] = sss_signals_sub0[i, 31:62]
    padded_sss_signals_sub0[i, 0] = 0
    padded_sss_signals_sub5[i, N-31:N] = sss_signals_sub5[i, :31]
    padded_sss_signals_sub5[i, 1:32] = sss_signals_sub5[i, 31:62]
    padded_sss_signals_sub5[i, 0] = 0
    
for i in range(168):
    padded_sss_signals_sub0[i] = np.fft.ifft(padded_sss_signals_sub0[i], N)
    padded_sss_signals_sub5[i] = np.fft.ifft(padded_sss_signals_sub5[i], N)
    

corr_sub0 = np.zeros((168,(len(waveform) + len(padded_sss_signals_sub0[0,:]) - 1)), dtype=complex)
corr_sub5 = np.zeros((168,(len(waveform) + len(padded_sss_signals_sub5[0,:]) - 1)), dtype=complex)

for i in range(168):
    corr_sub0[i,:] = np.correlate(padded_sss_signals_sub0[i], waveform, mode='full')
    corr_sub5[i,:] = np.correlate(padded_sss_signals_sub5[i], waveform, mode='full')

max_corr_sub0 = np.zeros(168)
max_corr_sub5 = np.zeros(168)
for i in range(168):
    max_corr_sub0[i] = max(abs(corr_sub0[i,:]))
    max_corr_sub5[i] = max(abs(corr_sub5[i,:]))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(np.abs(corr_sub0[0,:]))
plt.title("Correlation with SSS sub0")
plt.xlabel("Samples")
plt.ylabel("Correlation")
plt.subplot(2, 1, 2)
plt.plot(np.abs(corr_sub5[0,:]))
plt.title("Correlation with SSS sub5")
plt.xlabel("Samples")
plt.ylabel("Correlation")
plt.tight_layout()
plt.show()



print("----- SSS_detection -----")
print("Found NID1_sub0:",np.argmax(max_corr_sub0))
print("Found NID1_sub1:",np.argmax(max_corr_sub5))


print("----- NcellID -----")
print("NcellID:",np.argmax(max_corr_sub0)*3 + NID2)