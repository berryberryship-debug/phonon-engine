#!/usr/bin/env python3
"""
V57 — Moteur de dispersion phononique 1D (Chaîne diatomique).
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Paramètres 1D ---
M1 = 2.0 * 1.66054e-27
M2 = 4.0 * 1.66054e-27
K = 10.0
A = 5.0e-10

def get_parameters(q):
    S = 1.0/M1 + 1.0/M2
    delta = (4.0 / (M1 * M2)) * np.sin(q * A / 2.0)**2
    B = np.sqrt(max(S**2 - delta, 0.0))
    w2_aco = K * delta / (S + B) if (S + B) > 0 else 0.0
    w2_opt = K * (S + B)
    return np.sqrt(max(w2_aco, 0.0)), np.sqrt(max(w2_opt, 0.0))

def main():
    print("=== V57 : MOTEUR 1D ===")
    q_vals = np.linspace(-np.pi / A, np.pi / A, 501)
    data = []
    for q in q_vals:
        w_aco, w_opt = get_parameters(q)
        data.append([q, w_aco, w_opt])
    
    data = np.array(data)
    csv_file = "phonon_1d.csv"
    with open(csv_file, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["q_m^-1", "omega_aco_rad_s^-1", "omega_opt_rad_s^-1"])
        writer.writerows(data)
    
    print(f"[PASS] Données 1D exportées : {csv_file}")

if __name__ == "__main__":
    main()
