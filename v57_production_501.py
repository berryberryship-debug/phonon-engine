from __future__ import annotations
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

class LatticeModel:
    def __init__(self):
        self.m1 = 2.0 * 1.66054e-27
        self.m2 = 4.0 * 1.66054e-27
        self.k_spring = 10.0
        self.a = 5.0e-10

    def get_analytic_parameters(self, q):
        S = 1.0/self.m1 + 1.0/self.m2
        delta = (4.0 / (self.m1 * self.m2)) * np.sin(q * self.a / 2.0)**2
        B = np.sqrt(max(S**2 - delta, 0.0))
        w2_aco = self.k_spring * delta / (S + B) if (S + B) > 0 else 0.0
        w2_opt = self.k_spring * (S + B)
        return S, B, w2_aco, w2_opt

    def get_velocities(self, q):
        S, B, w2_aco, w2_opt = self.get_analytic_parameters(q)
        if B < 1e-15: return 0.0, 0.0
        w_aco = np.sqrt(max(w2_aco, 0.0))
        w_opt = np.sqrt(max(w2_opt, 0.0))
        prefactor = (self.k_spring * self.a * np.sin(q * self.a)) / (2.0 * self.m1 * self.m2 * B)
        vg_aco = prefactor / w_aco if w_aco > 0 else 0.0
        vg_opt = -prefactor / w_opt if w_opt > 0 else 0.0
        return vg_aco, vg_opt

def generate_science():
    print("=== MOTEUR DE DISPERSION PHONONIQUE (V57 - 501 POINTS) ===")
    csv_filename = "phonon_dispersion_501.csv"
    if os.path.exists(csv_filename):
        os.remove(csv_filename)
        print(f"[-] Ancien fichier {csv_filename} supprime.")

    lat = LatticeModel()
    q_vals = np.linspace(-np.pi / lat.a, np.pi / lat.a, 501)
    data = []
    
    for q in q_vals:
        S, B, w2_aco, w2_opt = lat.get_analytic_parameters(q)
        w_aco = np.sqrt(max(w2_aco, 0.0))
        w_opt = np.sqrt(max(w2_opt, 0.0))
        vg_aco, vg_opt = lat.get_velocities(q)
        data.append([q, w_aco, w_opt, vg_aco, vg_opt])
        
    data = np.array(data)

    with open(csv_filename, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["q_m^-1", "omega_acoustic_rad_s^-1", "omega_optical_rad_s^-1", "vg_aco_ms", "vg_opt_ms"])
        writer.writerows(data)
    print(f"[+] Donnees exportees : {csv_filename}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    ax1.plot(data[:, 0], data[:, 1], 'b-', label='Acoustique')
    ax1.plot(data[:, 0], data[:, 2], 'r-', label='Optique')
    ax1.set_title("Dispersion (501 points)")
    ax1.grid(True)
    ax2.plot(data[:, 0], data[:, 3], 'b-', label='v_g Aco')
    ax2.plot(data[:, 0], data[:, 4], 'r-', label='v_g Opt')
    ax2.set_title("Vitesses de groupe")
    ax2.grid(True)
    plt.tight_layout()
    plt.savefig("phonon_curves_501.png", dpi=300)
    print(f"[+] Graphiques generes : phonon_curves_501.png")

if __name__ == "__main__":
    generate_science()
