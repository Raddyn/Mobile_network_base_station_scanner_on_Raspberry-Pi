import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# 1. Načti IQ data
# data = sio.loadmat('data1_20Mhz.mat')
data = sio.loadmat('data2.mat')




iWave = data['iWave']
qWave = data['qWave']
iq = iWave.squeeze() + 1j * qWave.squeeze()


iq = np.load('data.npy')
print("Data1:", len(iq),"\n")

fs = 1.92e6                    # sample rate // sampling time

# 2. Parametry
fft_size = 2048
step = fft_size // 2
n_avg = 100 # kolik FFT oken zprůměrovat

# 3. Výpočet průměrného spektra
specs = []
for i in range(0, len(iq) - fft_size, step):
    window = iq[i:i+fft_size]
    spectrum = np.fft.fftshift(np.abs(np.fft.fft(window))**2)
    specs.append(spectrum)
    if len(specs) >= n_avg:
        break

avg_spectrum = np.mean(specs, axis=0)
freqs = np.fft.fftshift(np.fft.fftfreq(fft_size, d=1/fs))

# 4. Normalizace a prahování
avg_spectrum_db = 10 * np.log10(avg_spectrum + 1e-10)
threshold_db = np.max(avg_spectrum_db) - 10  # např. -10 dB od maxima
active_bins = freqs[(avg_spectrum_db > threshold_db)]

# 5. Výpočet šířky pásma
bw_detected = active_bins[-1] - active_bins[0]  # v Hz
bw_detected_mhz = np.round(bw_detected / 1e6, 1)

# 6. Zaokrouhlení na standardní LTE BW
standard_bw = np.array([1.4, 3, 5, 10, 15, 20])
closest_bw = standard_bw[np.argmin(np.abs(standard_bw - bw_detected_mhz))]

print(f"Detekovaná šířka pásma: {bw_detected_mhz} MHz (přiřazeno: {closest_bw} MHz)")
plt.plot(freqs / 1e6, avg_spectrum_db)
plt.axhline(threshold_db, color='red', linestyle='--')
plt.title("Průměrné spektrum")
plt.xlabel("Frekvence [MHz]")
plt.ylabel("Výkon [dB]")
plt.grid()
plt.show()