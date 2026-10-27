#!/usr/bin/env python3
"""
V73 — Validation rigoureuse des vitesses acoustiques du modèle V72.
"""

import numpy as np
import v72_second_neighbors as engine

def compute_velocities():
    print("=== V73 : VITESSES ACOUSTIQUES (V72 - SECONDS VOISINS) ===")
    
    # 20 points pour une bonne statistique de régression
    q_vals = np.linspace(1e6, 1e8, 20)
    w_ta1_list, w_ta2_list, w_la_list = [], [], []

    m1, m2 = engine.M1, engine.M2

    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        
        # Vérification de l'herméticité
        if not np.allclose(D, D.conj().T):
            raise AssertionError("La matrice D(q) n'est pas hermitienne.")

        w2, V = np.linalg.eigh(D)
        
        # Vérifications de sécurité
        assert np.all(np.isfinite(w2)), "Valeurs propres non finies."
        scale = np.max(np.abs(w2))
        assert scale > 0, "Spectre nul."
        assert w2[0] > -1e-9 * scale, "Instabilité détectée."
        
        w = np.sqrt(np.clip(w2, 0, None))
        
        # --- Identification robuste des modes ---
        # On construit le vecteur de translation rigide selon x (massivement pondéré)
        # pour identifier le mode longitudinal acoustique (LA).
        v_rigid_x = np.zeros(6)
        v_rigid_x[:3] = np.array([1, 0, 0]) * np.sqrt(m1)
        v_rigid_x[3:] = np.array([1, 0, 0]) * np.sqrt(m2)
        v_rigid_x /= np.linalg.norm(v_rigid_x)
        
        # On cherche les 3 modes acoustiques (fréquences les plus basses)
        idx_aco = np.argsort(w)[:3]
        
        # Calcul du recouvrement avec la translation rigide x
        overlaps = [np.abs(np.vdot(v_rigid_x, V[:, i]))**2 for i in idx_aco]
        
        # Le mode avec le plus grand recouvrement est LA
        i_la = idx_aco[np.argmax(overlaps)]
        i_ta = [i for i in idx_aco if i != i_la]
        
        w_la_list.append(w[i_la])
        w_ta1_list.append(w[i_ta[0]])
        w_ta2_list.append(w[i_ta[1]])

    # Fonction de régression avec erreur standard (SEM)
    # Note : Il s'agit d'une erreur de régression, pas d'une incertitude expérimentale.
    def fit_and_err(q, w):
        q = np.array(q)
        w = np.array(w)
        v = np.sum(q*w) / np.sum(q*q) # Régression forcée à l'origine
        n = len(q)
        sse = np.sum((w - v*q)**2)
        var_v = sse / (n - 1) / np.sum(q**2)
        return v, np.sqrt(var_v)

    v_ta1, err_ta1 = fit_and_err(q_vals, w_ta1_list)
    v_ta2, err_ta2 = fit_and_err(q_vals, w_ta2_list)
    v_la, err_la = fit_and_err(q_vals, w_la_list)

    # Vérification de la linéarité (en divisant par q pour chaque point)
    wl_ta1 = np.array(w_ta1_list) / q_vals
    wl_la = np.array(w_la_list) / q_vals
    
    lin_ta = (np.max(wl_ta1) / np.min(wl_ta1)) - 1
    lin_la = (np.max(wl_la) / np.min(wl_la)) - 1

    print(f"v_TA1 = {v_ta1:.2f} ± {err_ta1:.2f} m/s")
    print(f"v_TA2 = {v_ta2:.2f} ± {err_ta2:.2f} m/s")
    print(f"v_LA  = {v_la:.2f} ± {err_la:.2f} m/s")
    print(f"Linéarité TA (w/q) = {lin_ta:.2e}")
    print(f"Linéarité LA (w/q) = {lin_la:.2e}")
    
    ecart = abs(v_la - v_ta1) / max(v_ta1, 1e-10)
    print(f"Écart relatif LA/TA = {ecart:.1%}")

if __name__ == "__main__":
    compute_velocities()
