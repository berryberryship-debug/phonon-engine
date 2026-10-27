#!/usr/bin/env python3
"""
V70 — Validation acoustique (V61) par régression, avec identification par polarisation.
"""

import numpy as np
import v61_phonon_3d as engine

def analyze_sound_velocities():
    print("=== ESTIMATION DES VITESSES ACOUSTIQUES (Γ → X) ===")
    q_vals = np.linspace(1e6, 1e8, 10)
    w_ta, w_la = [], []
    
    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        
        # Diagonalisation complète
        w2, V = np.linalg.eigh(D)
        
        # Vérifications de sécurité
        assert np.all(np.isfinite(w2)), "Valeurs propres non finies."
        assert w2[:3].min() > -1e-9 * np.abs(w2).max(), "Instabilité détectée."
        
        # On garde les 3 branches acoustiques
        w = np.sqrt(np.clip(w2[:3], 0, None))
        
        # Identification de LA par la polarisation selon x (car q est selon x)
        px = np.abs(V[0, :3])**2 + np.abs(V[3, :3])**2
        i_la = int(np.argmax(px))
        
        w_la.append(w[i_la])
        w_ta.append(np.mean(np.delete(w, i_la)))

    def fit_velocity(w):
        v, _, _, _ = np.linalg.lstsq(q_vals.reshape(-1, 1), w, rcond=None)
        return v[0]

    v_ta = fit_velocity(w_ta)
    v_la = fit_velocity(w_la)
    
    print(f"v_TA (moyenne) = {v_ta:.2f} m/s")
    print(f"v_LA (polarisée) = {v_la:.2f} m/s")
    
    ecart = abs(v_la - v_ta) / max(v_ta, 1e-10)
    if ecart < 1e-3:
        print(f"[INFO] Dégénérescence LA/TA le long de [100] (écart = {ecart:.1e}).")
    else:
        print(f"[INFO] Vitesses distinctes (écart = {ecart:.1%}).")

if __name__ == "__main__":
    analyze_sound_velocities()
