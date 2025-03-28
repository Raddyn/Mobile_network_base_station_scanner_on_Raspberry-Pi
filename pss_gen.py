import numpy as np

def lte_pss_gen(NID2):
    i = NID2
    u_shift = [25, 29, 34]
    pss_sequence = np.array([
        [
            np.exp(-1j * np.pi * u_shift[i] * n * (n + 1) / 63)
            if n <= 30 else
            np.exp(-1j * np.pi * u_shift[i] * (n + 1) * (n + 2) / 63)
            for n in range(62)
        ]
    ]
    )
    return pss_sequence
    