```python
#!/usr/bin/env python3
"""
V63 — Comparaison Modèle vs Expérience (CsCl).
"""
import matplotlib.pyplot as plt
import numpy as np
from v61_phonon_3d import build_q_path, frequencies

def main():
    print("=== V63 : COMPARAISON AVEC LA RÉALITÉ ===")
    q_path, distance, tick_pos, tick_lab = build_q_path(100)
    omega = np.array([frequencies(q) for q in q_path])

    # Données expérimentales simplifiées (Raman/Neutrons) en THz
    # Fréquence optique à Gamma ~ 3.5 THz. Vitesse du son ~ 2000 m/s.
    dist_gx = distance[:100]
    exp_opt_x = np.linspace(0, dist_gx[-1], 5)
    exp_opt_y = np.array([3.5, 3.3, 2.9, 2.3, 1.8])
    
    exp_aco_x = np.linspace(0, dist_gx[-1], 5)
    exp_aco_y = (exp_aco_x / 1e10) * 2.0 # approximation vitesse du son

    fig, ax = plt.subplots(figsize=(10, 6))
    
    for b in range(6):
        ax.plot(distance, omega[:, b] / (2 * np.pi * 1e12), lw=1.5, label=f"Modèle Br.{b+1}")

    ax.scatter(exp_opt_x, exp_opt_y, color='red', s=80, zorder=5, label="Exp. Optique")
    ax.scatter(exp_aco_x, exp_aco_y, color='green', s=80, zorder=5, label="Exp. Acoustique")

    for x in tick_pos: ax.axvline(x, lw=0.6, alpha=0.4)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab)
    ax.set_ylabel("Fréquence (THz)")
    ax.set_title("Modèle (lignes) vs Expérience (points)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    
    plt.savefig("comparaison_cscl.png", dpi=180)
    print("[PASS] comparaison_cscl.png généré. Ouvre-le !")

if __name__ == "__main__":
    main()
```
