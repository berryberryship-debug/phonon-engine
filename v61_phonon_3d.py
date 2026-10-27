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

def frequencies(q, tolerance=1e-12):
    D = dynamical_matrix(q)
    eigvals = np.linalg.eigvalsh(D)
    scale = max(float(np.max(np.abs(eigvals))), np.finfo(float).tiny)
    
    if eigvals[0] < -tolerance * scale:
        raise ValueError(f"Instabilité à q={q}")
        
    if np.all(q == 0.0):
        eigvals[np.abs(eigvals) < tolerance * scale] = 0.0
    else:
        eigvals[(eigvals < 0.0) & (eigvals >= -tolerance * scale)] = 0.0
        
    return np.sqrt(eigvals)

def build_q_path(n_segment):
    PI_A = np.pi / A
    high_symmetry = [
        ("Γ", np.array([0.0, 0.0, 0.0])),
        ("X", np.array([PI_A, 0.0, 0.0])),
        ("M", np.array([PI_A, PI_A, 0.0])),
        ("R", np.array([PI_A, PI_A, PI_A])),
        ("Γ", np.array([0.0, 0.0, 0.0])),
    ]

    q_points = []
    tick_positions = []
    tick_labels = []

    for point_index, (label, q_point) in enumerate(high_symmetry):
        if point_index > 0:
            q_previous = high_symmetry[point_index - 1][1]
            for step in range(n_segment):
                t = step / n_segment
                q_points.append((1.0 - t) * q_previous + t * q_point)

        tick_positions.append(len(q_points))
        tick_labels.append(label)
        q_points.append(q_point)

    q_path = np.asarray(q_points, dtype=float)
    differences = np.linalg.norm(np.diff(q_path, axis=0), axis=1)
    distance = np.concatenate(([0.0], np.cumsum(differences)))

    return q_path, distance, np.asarray(tick_positions), tick_labels

def main():
    print("=== V61 : MOTEUR 3D (PREMIERS VOISINS) ===")
    q_path, distance, tick_positions, tick_labels = build_q_path(100)
    omega = np.array([frequencies(q) for q in q_path])
    print(f"[PASS] {len(q_path)} points calculés.")

if __name__ == "__main__":
    main()
