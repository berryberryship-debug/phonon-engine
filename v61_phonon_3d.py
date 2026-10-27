#!/usr/bin/env python3
"""
V61 — Moteur de dispersion phononique 3D, réseau cubique type CsCl.
"""

from itertools import product
from pathlib import Path
import csv
import sys

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1. Paramètres du modèle
# ------------------------------------------------------------
U = 1.66053906660e-27
M1 = 132.90545196 * U
M2 = 126.90447 * U
A = 4.12e-10
K_PARALLEL = 10.0
K_PERP = 2.0

N_SEGMENT = 100
OUT_DIR = Path(".")
CSV_FILE = OUT_DIR / "phonon_3d_cscl.csv"
PNG_FILE = OUT_DIR / "phonon_3d_cscl.png"

# ------------------------------------------------------------
# 2. Géométrie des voisins
# ------------------------------------------------------------
I3 = np.eye(3)
NEIGHBORS = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    rhat = R / np.linalg.norm(R)
    P = np.outer(rhat, rhat)
    K = K_PERP * I3 + (K_PARALLEL - K_PERP) * P
    NEIGHBORS.append((R, K))

# ------------------------------------------------------------
# 3. Matrice dynamique D(q), dimension 6 x 6
# ------------------------------------------------------------
def dynamical_matrix(q):
    q = np.asarray(q, dtype=float)
    if q.shape != (3,):
        raise ValueError("q doit contenir exactement trois composantes.")
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
        raise ValueError(f"Instabilité à q={q} : omega^2={eigvals[0]:.6e}")

    if np.all(q == 0.0):
        eigvals[np.abs(eigvals) < tolerance * scale] = 0.0
    else:
        eigvals[(eigvals < 0.0) & (eigvals >= -tolerance * scale)] = 0.0

    return np.sqrt(eigvals)

# ------------------------------------------------------------
# 4. Trajet de haute symétrie
# ------------------------------------------------------------
PI_A = np.pi / A
HIGH_SYMMETRY = [
    ("Γ", np.array([0.0, 0.0, 0.0])),
    ("X", np.array([PI_A, 0.0, 0.0])),
    ("M", np.array([PI_A, PI_A, 0.0])),
    ("R", np.array([PI_A, PI_A, PI_A])),
    ("Γ", np.array([0.0, 0.0, 0.0])),
]

def build_q_path(n_segment):
    if n_segment < 1:
        raise ValueError("n_segment doit être >= 1")
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

# ------------------------------------------------------------
# 5. Audits
# ------------------------------------------------------------
def run_audits():
    print("=== V61 : AUDIT DU MODÈLE PHONONIQUE 3D ===")
    D0 = dynamical_matrix(np.zeros(3))
    evals0 = np.linalg.eigvalsh(D0)
    scale0 = max(np.max(np.abs(evals0)), 1.0)
    n_zero = np.count_nonzero(np.abs(evals0) < 1e-10 * scale0)
    if n_zero != 3:
        raise AssertionError(f"Attendu : 3 modes nuls à Gamma ; obtenu : {n_zero}")
    print("[PASS] Trois modes acoustiques nuls à Γ")

    q_path, _, _, _ = build_q_path(N_SEGMENT)
    for q in q_path:
        evals = np.linalg.eigvalsh(dynamical_matrix(q))
        scale = max(float(np.max(np.abs(evals))), 1.0)
        if evals[0] < -1e-7 * scale:
            raise AssertionError(f"Valeur propre négative à q={q}: {evals[0]:.6e}")
    print("[PASS] Positivité numérique sur le trajet")

def audit_analytical_references():
    print("=== AUDIT ANALYTIQUE DES POINTS DE HAUTE SYMÉTRIE ===")
    alpha0 = (K_PARALLEL + 2.0 * K_PERP) / 3.0
    beta0 = (K_PARALLEL - K_PERP) / 3.0
    w_A = np.sqrt(8.0 * alpha0 / M1)
    w_B = np.sqrt(8.0 * alpha0 / M2)
    w_opt_G = np.sqrt(8.0 * alpha0 * (1.0 / M1 + 1.0 / M2))
    
    delta = 8.0 * beta0 / np.sqrt(M1 * M2)
    tA = 8.0 * alpha0 / M1
    tB = 8.0 * alpha0 / M2
    tr = 0.5 * (tA + tB)
    diff = 0.5 * (tA - tB)
    disc = np.sqrt(diff**2 + delta**2)
    w_M_minus = np.sqrt(tr - disc)
    w_M_plus = np.sqrt(tr + disc)

    expected = {
        "Γ": np.sort([0.0, 0.0, 0.0, w_opt_G, w_opt_G, w_opt_G]),
        "X": np.sort([w_A, w_A, w_A, w_B, w_B, w_B]),
        "M": np.sort([w_M_minus, w_M_minus, w_A, w_B, w_M_plus, w_M_plus]),
        "R": np.sort([w_A, w_A, w_A, w_B, w_B, w_B]),
    }
    q_dict = {
        "Γ": np.zeros(3),
        "X": np.array([PI_A, 0.0, 0.0]),
        "M": np.array([PI_A, PI_A, 0.0]),
        "R": np.array([PI_A, PI_A, PI_A]),
    }

    for name, q in q_dict.items():
        w_computed = frequencies(q)
        w_expected = expected[name]
        np.testing.assert_allclose(
            w_computed / 1e12, w_expected / 1e12,
            rtol=1e-6, atol=1e-8,
            err_msg=f"Échec concordance fréquences au point {name}",
        )

        evals = np.sort(np.linalg.eigvalsh(dynamical_matrix(q)))
        scale = max(float(np.max(np.abs(evals))), np.finfo(float).tiny)
        tol = 1e-8 * scale

        groups = []
        for val in evals:
            if not groups or abs(val - groups[-1][-1]) > tol:
                groups.append([val])
            else:
                groups[-1].append(val)
        multiplicities = tuple(len(group) for group in groups)
        
        expected_mult = {"Γ": (3, 3), "X": (3, 3), "M": (2, 1, 1, 2), "R": (3, 3)}[name]
        if multiplicities != expected_mult:
            raise AssertionError(f"Multiplicités incorrectes au point {name} : obtenu {multiplicities}")
        print(f"[PASS] Point {name:2s} : fréquences et multiplicités {multiplicities} validées.")

# ------------------------------------------------------------
# 6. Export
# ------------------------------------------------------------
def main():
    run_audits()
    audit_analytical_references()

    q_path, distance, tick_positions, tick_labels = build_q_path(N_SEGMENT)
    omega = np.array([frequencies(q) for q in q_path])

    with CSV_FILE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["q_x_m^-1", "q_y_m^-1", "q_z_m^-1", "distance_path_m^-1", "omega_1_rad_s^-1", "omega_2_rad_s^-1", "omega_3_rad_s^-1", "omega_4_rad_s^-1", "omega_5_rad_s^-1", "omega_6_rad_s^-1"])
        for q, d, w in zip(q_path, distance, omega):
            writer.writerow([*q, d, *w])

    plt.figure(figsize=(10, 6))
    for branch in range(6):
        plt.plot(distance, omega[:, branch] / 1e12, linewidth=1.5, label=f"Branche {branch + 1}")

    for x in tick_positions:
        plt.axvline(x, linewidth=0.6, alpha=0.4)
    plt.xticks(tick_positions, tick_labels)
    plt.xlabel(r"Chemin dans l'espace réciproque (\(m^{-1}\))")
    plt.ylabel(r"Fréquence angulaire \(\omega\) (\(10^{12}\) rad/s)")
    plt.title("Dispersion phononique 3D — réseau cubique type CsCl")
    plt.grid(True, alpha=0.25)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(PNG_FILE, dpi=180)
    plt.close()

    print(f"[PASS] CSV exporté : {CSV_FILE}")
    print(f"[PASS] Graphique exporté : {PNG_FILE}")
    print("V61 terminé.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ÉCHEC V61 : {exc}", file=sys.stderr)
        sys.exit(1)
