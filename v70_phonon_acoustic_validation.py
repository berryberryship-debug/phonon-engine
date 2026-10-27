#!/usr/bin/env python3
"""
V70 — Validation acoustique (V61) avec sélection dynamique et incertitudes.
"""

import numpy as np
import v61_phonon_3d as engine

def analyze_sound_velocities():
    print("=== ESTIMATION DES VITESSES ACOUSTIQUES (Γ → X) ===")
    
    q_vals = np.linspace(1e6, 1e8, 20)
    w_ta_list, w_la_list = [], []

    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        
        # Vérification hermitienne
        assert np.allclose(D, D.conj().T), "Matrice non hermitienne."
        
        w2, V = np.linalg.eigh(D)
        
        assert np.all(np.isfinite(w2)), "Valeurs propres non finies."
        scale = np.max(np.abs(w2))
        assert w2[0] > -1e-9 * scale, "Instabilité détectée."
        
        w = np.sqrt(np.clip(w2, 0, None))
        
        # Dé-massification et calcul du centre de masse (CM)
        m1, m2 = engine.M1, engine.M2
        u_cm = (np.sqrt(m1) * V[:3, :] + np.sqrt(m2) * V[3:, :]) / (m1 + m2)
        
        # Pureté longitudinale = projection sur x / norme
        px = np.abs(u_cm[0, :])**2 / (np.linalg.norm(u_cm, axis=0)**2 + 1e-30)
        
        idx_aco = np.argsort(w)[:3]
        i_la = idx_aco[np.argmax(px[idx_aco])]
        i_ta = [i for i in idx_aco if i != i_la]
        
        w_la_list.append(w[i_la])
        w_ta_list.append(np.mean(w[i_ta]))

    # Régression avec erreur standard (SEM)
    def fit_and_err(q, w):
        q = np.array(q)
        w = np.array(w)
        v = np.sum(q*w) / np.sum(q*q)
        n = len(q)
        sse = np.sum((w - v*q)**2)
        var_v = sse / (n - 1) / np.sum(q**2)
        return v, np.sqrt(var_v)

    v_ta, err_ta = fit_and_err(q_vals, w_ta_list)
    v_la, err_la = fit_and_err(q_vals, w_la_list)

    # Vérification de la linéarité (w/q)
    wl_la = np.array(w_la_list) / q_vals
    lin_la = (np.max(wl_la) / np.min(wl_la)) - 1

    print(f"v_TA = {v_ta:.2f} ± {err_ta:.2f} m/s")
    print(f"v_LA = {v_la:.2f} ± {err_la:.2f} m/s")
    print(f"Linéarité LA (w/q) = {lin_la:.2e}")
    print(f"Écart LA/TA = {abs(v_la - v_ta):.2f} m/s")

if __name__ == "__main__":
    analyze_sound_velocities()
