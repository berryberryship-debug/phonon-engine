#!/usr/bin/env python3
"""
V61 — Moteur de dispersion phononique 3D pour le CsCl.
"""

from itertools import product
from pathlib import Path
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

# --- Paramètres Physiques Réels du CsCl ---
U = 1.66053906660e-27
M1 = 132.90545196 * U  # Masse du Cs
M2 = 35.45 * U         # Masse du Cl (Corrigée)
A = 4.12e-10           # Paramètre de maille
K_PARALLEL = 10.0      # Constante de force (à ajuster)
K_PERP = 2.0           # Constante de force (à ajuster)

N_SEGMENT = 100
CSV_FILE = Path("phonon_3d_cscl.csv")
PNG_FILE = Path("phonon_3d_cscl.png")

# --- Géométrie (8 voisins) ---
I3 = np.eye(3)
NEIGHBORS = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    rhat = R / np.linalg.norm(R)
    P = np.outer(rhat, rhat)
    K = K_PERP * I3 + (K_PARALLEL - K_PERP) * P
    NEIGHBORS.append((R, K))

# --- Matrice Dynamique ---
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

def frequencies(q, tolerance=1e-7):
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

# --- Trajet de haute symétrie ---
PI_A = np.pi / A
HIGH_SYMMETRY = [
    ("Γ", np.array([0.0, 0.0, 0.0])),
    ("X", np.array([PI_A, 0.0, 0.0])),
    ("M", np.array([PI_A, PI_A, 0.0])),
    ("R", np.array([PI_A, PI_A, PI_A])),
    ("Γ", np.array([0.0, 0.0, 0.0])),
]

def build_q_path(n_segment):
    points = []
    for seg in range(len(HIGH_SYMMETRY) - 1):
        q0 = HIGH_SYMMETRY[seg][1]
        q1 = HIGH_SYMMETRY[seg + 1][1]
        for j in range(n_segment):
            t = j / n_segment
            points.append((1.0 - t) * q0 + t * q1)
    points.append(HIGH_SYMMETRY[-1][1])
    
    points = np.array(points)
    distances = np.zeros(len(points))
    for i in range(1, len(points)):
        distances[i] = distances[i-1] + np.linalg.norm(points[i] - points[i-1])
        
    tick_positions = [distances[i * n_segment] for i in range(len(HIGH_SYMMETRY))]
    tick_labels = [x[0] for x in HIGH_SYMMETRY]
    return points, distances, tick_positions, tick_labels

# --- Génération CSV et Graphique ---
def main():
    q_path, distance, tick_positions, tick_labels = build_q_path(N_SEGMENT)
    omega = np.array([frequencies(q) for q in q_path])

    with CSV_FILE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["q_x", "q_y", "q_z", "dist", "w1", "w2", "w3", "w4", "w5", "w6"])
        for q, d, w in zip(q_path, distance, omega):
            writer.writerow([*q, d, *w])

    plt.figure(figsize=(10, 6))
    for branch in range(6):
        plt.plot(distance, omega[:, branch] / 1e12, linewidth=1.5, label=f"Br. {branch+1}")

    for x in tick_positions: plt.axvline(x, linewidth=0.6, alpha=0.4)
    plt.xticks(tick_positions, tick_labels)
    plt.ylabel(r"Fréquence \(\omega\) (\(10^{12}\) rad/s)")
    plt.title("Dispersion 3D — CsCl")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.savefig(PNG_FILE, dpi=180)
    print("[PASS] phonon_3d_cscl.png et .csv générés.")

if __name__ == "__main__":
    main()
