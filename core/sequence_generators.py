import numpy as np


def pss_gen(NID2):
    """
    Generate the Primary Synchronization Signal (PSS) for a given NID2.

    Parameters:
        NID2 (int): The NID2 value (0, 1, or 2).


    Returns:
        np.ndarray: The generated PSS sequence.
    """
    # -------------------------------------------------------

    if NID2 not in [0, 1, 2]:
        raise ValueError("NID2 must be 0, 1, or 2.")

    # -------------------------------------------------------
    i = NID2
    u_shift = [25, 29, 34]
    pss = np.array(
        [
            [
                np.exp(-1j * np.pi * u_shift[i] * n * (n + 1) / 63)
                if n <= 30
                else np.exp(-1j * np.pi * u_shift[i] * (n + 1) * (n + 2) / 63)
                for n in range(62)
            ]
        ]
    )
    return pss


def sss_gen(NID1, NID2):
    """
    Generate the Secondary Synchronization Signal (SSS) for given NID1 and NID2.

    Parameters:
        NID1 (int): The NID1 value (0 to 167).
        NID2 (int): The NID2 value (0 to 2).

    Returns:
        tuple: Two numpy arrays representing the SSS sequences.

        sss0 (np.ndarray): The first SSS sequence.
        sss5 (np.ndarray): The second SSS sequence.
    """
    # -------------------------------------------------------
    if NID1 < 0 or NID1 > 167:
        raise ValueError("NID1 must be between 0 and 167.")
    if NID2 < 0 or NID2 > 2:
        raise ValueError("NID2 must be between 0 and 2.")
    # -------------------------------------------------------
    q_prime = np.floor(NID1 / 30)
    q = np.floor((NID1 + q_prime * (q_prime + 1) / 2) / 30)
    m_prime = NID1 + q * (q + 1) / 2
    m0 = m_prime % 31
    m1 = (m0 + np.floor(m_prime / 31) + 1) % 31
    x_s = np.zeros(31)
    x_s[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_s[i + 5] = (x_s[i + 2] + x_s[i]) % 2
    x_c = np.zeros(31)
    x_c[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_c[i + 5] = (x_c[i + 3] + x_c[i]) % 2
    s_tilda = 1 - 2 * x_s
    c_tilda = 1 - 2 * x_c
    s0_m0_even = np.zeros(31)
    s1_m1_even = np.zeros(31)
    for n in range(31):
        s0_m0_even[n] = s_tilda[int((n + m0) % 31)]
        s1_m1_even[n] = s_tilda[int((n + m1) % 31)]
    c0_even = np.zeros(31)
    for n in range(31):
        c0_even[n] = c_tilda[int((n + NID2) % 31)]
    d_even_sub0 = s0_m0_even * c0_even
    d_even_sub5 = s1_m1_even * c0_even
    x_z = np.zeros(31)
    x_z[:5] = [0, 0, 0, 0, 1]
    for i in range(26):
        x_z[i + 5] = (x_z[i + 4] + x_z[i + 2] + x_z[i + 1] + x_z[i]) % 2
    z_tilda = 1 - 2 * x_z
    s1_m1_odd = np.zeros(31)
    s0_m0_odd = np.zeros(31)
    for n in range(31):
        s1_m1_odd[n] = s_tilda[int((n + m1) % 31)]
        s0_m0_odd[n] = s_tilda[int((n + m0) % 31)]
    c1_odd = np.zeros(31)
    for n in range(31):
        c1_odd[n] = c_tilda[int((n + NID2 + 3) % 31)]
    z1_m0_odd = np.zeros(31)
    z1_m1_odd = np.zeros(31)
    for n in range(31):
        z1_m0_odd[n] = z_tilda[int((n + m0 % 8) % 31)]
        z1_m1_odd[n] = z_tilda[int((n + m1 % 8) % 31)]
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
