#!/usr/bin/env python3
"""
V70 — Validation acoustique (V61) avec sélection dynamique des modes et incertitude statistique.
"""

import numpy as np
import v61_phonon_3d as engine

def analyze_sound_velocities():
    print("=== ESTIMATION DES VITESSES ACOUSTIQUES (Γ → X) ===")
    q_vals = np.linspace(1e6, 1e8, 20)  # Plus de points pour une meilleure statistique
    w_ta, w_la = [], []

    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        w2, V = np.linalg.eigh(D)
        
        # Vérifications de sécurité
        assert np.all(np.isfinite(w2)), "Valeurs propres non finies."
        scale = np.max(np.abs(w2))
        assert scale > 0, "Spectre nul."
        assert w2[0] > -1e-9 * scale, "Instabilité détectée."
        
        w = np.sqrt(np.clip(w2, 0, None))
        
        # Calcul du déplacement du centre de masse (CM)
        m1, m2 = engine.M1, engine.M2
        # V[:3, :] est Cs, V[3:, :] est Cl
        u_cm = (m1 * V[:3, :] + m2 * V[3:, :]) / (m1 + m2)
        
        # Pureté longitudinale = projection sur x / norme totale du CM
        px = np.abs(u_cm[0, :])**2 / (np.linalg.norm(u_cm, axis=0)**2 + 1e-30)
        
        # Les 3 modes acoustiques sont les 3 plus basses fréquences
        idx_aco = np.argsort(w)[:3]
        
        # Parmi ces 3, LA est celui avec la plus grande pureté longitudinale
        i_la = idx_aco[np.argmax(px[idx_aco])]
        i_ta = [i for i in idx_aco if i != i_la]
        
        w_la.append(w[i_la])
        w_ta.append(np.mean(w[i_ta]))
        
        if q_val == q_vals[0]:
            print(f"Pureté LA (q={q_val:.1e}) = {px[i_la]:.4f}")

    # Régression avec erreur standard (forcée à passer par l'origine)
    def fit_and_err(q, w):
        q = np.array(q)
        w = np.array(w)
        v = np.sum(q*w) / np.sum(q*q)
        n = len(q)
        sse = np.sum((w - v*q)**2)
        var_v = sse / (n - 1) / np.sum(q**2)
        return v, np.sqrt(var_v)

    v_ta, err_ta = fit_and_err(q_vals, w_ta)
    v_la, err_la = fit_and_err(q_vals, w_la)

    print(f"v_TA (moyenne) = {v_ta:.2f} ± {err_ta:.2f} m/s")
    print(f"v_LA (polarisée) = {v_la:.2f} ± {err_la:.2f} m/s")
    print(f"Écart = {abs(v_la - v_ta):.2f} m/s")

if __name__ == "__main__":
    analyze_sound_velocities()
