#!/usr/bin/env python3
"""
V74 — Audit de stabilité sur la zone de Brillouin du réseau de Bravais sous-jacent au CsCl.
"""

import numpy as np
import v72_second_neighbors as engine

def audit_3d_grid(n_points=15):
    print("=== V74 : AUDIT DE STABILITÉ 3D COMPLET (SECONDS VOISINS) ===")
    
    # --- Test explicite des 3 modes nuls à Gamma ---
    D_gamma = engine.dynamical_matrix(np.zeros(3))
    evals_gamma = np.linalg.eigvalsh(D_gamma)
    scale_gamma = max(np.max(np.abs(evals_gamma)), 1.0)
    n_zero = np.sum(np.abs(evals_gamma) < 1e-12 * scale_gamma)
    if n_zero != 3:
        raise AssertionError(f"Attendu 3 modes nuls à Gamma, obtenu {n_zero}.")
    print("[PASS] 3 modes acoustiques nuls à Gamma.")

    # --- Grille sur la zone de Brillouin du réseau de Bravais sous-jacent ---
    PI_A = np.pi / engine.A
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    
    min_eval = np.inf
    min_eval_rel = np.inf
    min_q = None
    instabilites = 0
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                D = engine.dynamical_matrix(q)
                
                # Vérification de l'herméticité avec tolérance relative
                norm_D = np.linalg.norm(D, 2)
                if norm_D > 0:
                    diff = np.max(np.abs(D - D.conj().T))
                    if diff > 1e-12 * norm_D:
                        raise AssertionError(f"Matrice non hermitienne à q={q}")
                
                evals = np.linalg.eigvalsh(D)
                
                # Normalisation par l'échelle spectrale locale
                scale = max(np.max(np.abs(evals)), 1.0)
                current_min_rel = evals[0] / scale
                
                if evals[0] < min_eval:
                    min_eval = evals[0]
                    min_eval_rel = current_min_rel
                    min_q = q.copy()
                    
                # Seuil de diagnostic (choix numérique, pas une frontière physique universelle)
                if current_min_rel < -1e-5:
                    instabilites += 1
                    
    print(f"Grille testée : {n_points}^3 = {n_points**3} points")
    print(f"Valeur propre minimale absolue : {min_eval:.4e} s^-2")
    print(f"Valeur propre minimale relative : {min_eval_rel:.4e}")
    print(f"Au point q = [{min_q[0]:.4e}, {min_q[1]:.4e}, {min_q[2]:.4e}] m^-1")
    
    if instabilites > 0:
        print(f"[FAIL] {instabilites} instabilités physiques détectées !")
        return False
    else:
        print("[PASS] Aucune instabilité physique détectée sur la grille.")
        print("(Note : Ceci est un test numérique sur une grille, et ne constitue pas une preuve analytique sur tout le continuum).")
        return True

if __name__ == "__main__":
    audit_3d_grid()
