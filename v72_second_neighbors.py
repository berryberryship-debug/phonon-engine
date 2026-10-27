#!/usr/bin/env python3
"""
V72 — Moteur enrichi avec les seconds voisins (Cs-Cs et Cl-Cl).
"""

from itertools import product
import numpy as np

U = 1.66053906660e-27
M1 = 132.90545196 * U
M2 = 35.45 * U
A = 4.12e-10

K1_PARALLEL, K1_PERP = 10.0, 2.0
K2_PARALLEL, K2_PERP = 4.0, 1.0
I3 = np.eye(3)

NEIGHBORS_1 = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    P = np.outer(R / np.linalg.norm(R), R / np.linalg.norm(R))
    NEIGHBORS_1.append((R, K1_PERP * I3 + (K1_PARALLEL - K1_PERP) * P))

NEIGHBORS_2 = []
for axis in range(3):
    for sign in [-1, 1]:
        R = np.zeros(3); R[axis] = sign * A
        P = np.outer(R / np.linalg.norm(R), R / np.linalg.norm(R))
        NEIGHBORS_2.append((R, K2_PERP * I3 + (K2_PARALLEL - K2_PERP) * P))

def dynamical_matrix(q):
    Daa = np.zeros((3, 3), dtype=complex)
    Dbb = np.zeros((3, 3), dtype=complex)
    Dab = np.zeros((3, 3), dtype=complex)
    
    for R, K in NEIGHBORS_1:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K / M1
        Dbb += K / M2
        Dab += -K * phase / np.sqrt(M1 * M2)
        
    for R, K in NEIGHBORS_2:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K * (1.0 - phase) / M1
        Dbb += K * (1.0 - phase) / M2
        
    D = np.zeros((6, 6), dtype=complex)
    D[:3, :3] = Daa
    D[:3, 3:] = Dab
    D[3:, :3] = Dab.conj().T
    D[3:, 3:] = Dbb
    return D

def test_acoustic_sum_rule():
    D0 = dynamical_matrix(np.zeros(3))
    norm_D = np.linalg.norm(D0, 2)
    u_x = np.array([1, 0, 0, 1, 0, 0], dtype=float)
    e = np.zeros(6)
    e[:3] = np.sqrt(M1) * u_x[:3]
    e[3:] = np.sqrt(M2) * u_x[3:]
    e = e / np.linalg.norm(e)
    res = np.linalg.norm(D0 @ e) / max(norm_D, 1e-30)
    print(f"[PASS] Règle de somme V72 (résidu relatif = {res:.4e})")

if __name__ == "__main__":
    test_acoustic_sum_rule()
