import adi 
import numpy as np
import matplotlib.pyplot as plt

# SDR nastavení
fs = int(30.72e6)  # vzorkovací frekvence
sample_time = 0.1  # čas vzorkování v sekundách

sdr = adi.Pluto(uri='ip:192.168.2.1')
sdr.rx_lo = int(796e6)
sdr.rx_buffer_size = int(fs * sample_time)
sdr.sample_rate = fs

# Načti data
iq = sdr.rx()

# Parametry FFT
fft_size = 2048
step = fft_size // 2
n_avg = 100  # kolik FFT oken zprůměrovat

# Výpočet průměrného spektra
specs = []
for i in range(0, len(iq) - fft_size, step):
    window = iq[i:i + fft_size]
    spectrum = np.fft.fftshift(np.abs(np.fft.fft(window))**2)
    specs.append(spectrum)
    if len(specs) >= n_avg:
        break

avg_spectrum = np.mean(specs, axis=0)
freqs = np.fft.fftshift(np.fft.fftfreq(fft_size, d=1 / fs))
avg_spectrum_db = 10 * np.log10(avg_spectrum + 1e-10)

# === Nový přístup: prahování podle středu pásma ===
center_bin = len(avg_spectrum_db) // 2
bw_ref_mhz = 1.0
bw_ref_hz = bw_ref_mhz * 1e6
bin_width = fs / fft_size
half_bin_count = int(bw_ref_hz / bin_width / 2)

# Výpočet referenčního výkonu ve středu
center_region = avg_spectrum_db[center_bin - half_bin_count : center_bin + half_bin_count]
reference_power = np.mean(center_region)  # nebo np.median(center_region)
threshold_db = reference_power - 10  # adaptivní práh

# Aktivní spektrální biny podle nového prahu
active_bins = freqs[(avg_spectrum_db > threshold_db)]
bw_detected = active_bins[-1] - active_bins[0] if len(active_bins) > 1 else 0
bw_detected_mhz = np.round(bw_detected / 1e6, 1)

# Přiřazení k nejbližší LTE šířce pásma
standard_bw = np.array([1.4, 3, 5, 10, 15, 20])
closest_bw = standard_bw[np.argmin(np.abs(standard_bw - bw_detected_mhz))]

print(f"Detekovaná šířka pásma: {bw_detected_mhz} MHz (přiřazeno: {closest_bw} MHz)")

# === Vizualizace ===
plt.figure(figsize=(10, 5))
plt.plot(freqs / 1e6, avg_spectrum_db)
plt.axhline(threshold_db, color='red', linestyle='--', label="Prahová úroveň")
plt.axvspan(freqs[center_bin - half_bin_count] / 1e6, freqs[center_bin + half_bin_count] / 1e6, 
            color='green', alpha=0.2, label="Referenční pásmo")
plt.title("Průměrné spektrum")
plt.xlabel("Frekvence [MHz]")
plt.ylabel("Výkon [dB]")
plt.legend()
plt.grid()

plt.figure(figsize=(10, 5))
plt.specgram(iq, NFFT=fft_size, Fs=fs, noverlap=fft_size//2, scale='dB')
plt.title("Spektrum LTE signálu")
plt.xlabel("Čas [s]")
plt.ylabel("Frekvence [Hz]")
plt.colorbar(label='Síla signálu [dB]')
plt.grid()
plt.tight_layout()

plt.show()