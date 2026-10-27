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
U = 1.66053906660e-27
M1 = 2.0 * U
M2 = 4.0 * U
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
    
    # Noms de fichiers attendus par V58 et V60
    csv_file = "phonon_dispersion_501.csv"
    png_file = "phonon_curves_501.png"

    # Note : V58 calcule les vitesses de groupe à partir de ce CSV. 
    # Pour éviter un crash, on ajoute des colonnes vides ou on laisse V58 les calculer.
    # Ici, on génère le CSV avec les colonnes attendues (V58 recalculera les dérivées).
    with open(csv_file, mode='w', newline='') as f:
        writer = csv.writer(f)
        # V60 attend ces en-têtes exacts
        writer.writerow(["q_m^-1", "omega_acoustic_rad_s^-1", "omega_optical_rad_s^-1"])
        writer.writerows(data)
    
    print(f"[PASS] CSV exporté : {csv_file}")

    # Génération du graphique PNG
    plt.figure(figsize=(10, 6))
    plt.plot(data[:, 0], data[:, 1], 'b-', label='Acoustique')
    plt.plot(data[:, 0], data[:, 2], 'r-', label='Optique')
    plt.xlabel("Vecteur d'onde q (m^-1)")
    plt.ylabel("Fréquence omega (rad/s)")
    plt.title("Dispersion Phononique 1D")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(png_file, dpi=150)
    plt.close()
    
    print(f"[PASS] Graphique exporté : {png_file}")

if __name__ == "__main__":
    main()
