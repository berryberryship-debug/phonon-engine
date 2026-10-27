#!/usr/bin/env python3
"""
V61 — Moteur de dispersion phononique 3D, réseau cubique type CsCl.
"""

from itertools import product
import numpy as np

# --- Paramètres Physiques ---
U = 1.66053906660e-27
M1 = 132.90545196 * U  # Cs
M2 = 35.45 * U         # Cl
A = 4.12e-10
K_PARALLEL = 10.0
K_PERP = 2.0

I3 = np.eye(3)
NEIGHBORS = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    rhat = R / np.linalg.norm(R)
    P = np.outer(rhat, rhat)
    K = K_PERP * I3 + (K_PARALLEL - K_PERP) * P
    NEIGHBORS.append((R, K))

def dynamical_matrix(q):
    q = np.asarray(q, dtype=float)
    Daa = np.zeros((3, 3), dtype=complex)
    Dbb = np.zeros((3, 3), dtype=complex)
    Dab = np.zeros((3, 3), dtype=complex)
    for R, K in NEIGHBORS:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K / M1
        Dbb += K / M2
        Dab += -K * phase / np.sqrt(M1 * M2)
    D = np.zeros((6, 6), dtype=complex)
    D[:3, :3] = Daa
    D[:3, 3:] = Dab
    D[3:, :3] = Dab.conj().T
    D[3:, 3:] = Dbb
    return D
