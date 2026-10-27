#!/usr/bin/env python3
"""
V74 — Audit de stabilité sur la zone de Brillouin du réseau de Bravais sous-jacent au CsCl.
"""

import numpy as np
import v72_second_neighbors as engine

def audit_3d_grid(n_points=15):
    print("=== V74 : AUDIT DE STABILITÉ 3D COMPLET (SECONDS VOISINS) ===")
    
    # --- Test explicite des 3 modes nuls à Gamma et du gap optique ---
    D_gamma = engine.dynamical_matrix(np.zeros(3))
    evals_gamma = np.linalg.eigvalsh(D_gamma)
    norm_D = np.linalg.norm(D_gamma, 2)
    
    if norm_D <= 0:
        raise AssertionError("La matrice D(0) est nulle.")
        
    tol = 1e-12 * norm_D
    acoustic = evals_gamma[:3]
    optic = evals_gamma[3:]
    
    if not np.all(np.abs(acoustic) < tol):
        raise AssertionError(f"Modes acoustiques non nuls à Gamma : {acoustic}")
    if not np.all(optic > tol):
        raise AssertionError(f"Modes optiques non strictement positifs à Gamma : {optic}")
        
    print("[PASS] 3 modes acoustiques nuls à Gamma et gap optique strictement positif.")

    # --- Grille sur la zone de Brillouin du réseau de Bravais sous-jacent ---
    PI_A = np.pi / engine.A
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    
    min_eval = np.inf
    min_eval_rel = np.inf
    min_q = None
    instabilites = 0
    SEUIL_RELATIF = -1e-5
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                D = engine.dynamical_matrix(q)
                
                # Vérification de l'herméticité
                norm_D_q = np.linalg.norm(D, 2)
                if norm_D_q > 0:
                    diff = np.max(np.abs(D - D.conj().T))
                    if diff > 1e-12 * norm_D_q:
                        raise AssertionError(f"Matrice non hermitienne à q={q}")
                
                evals = np.linalg.eigvalsh(D)
                
                # Normalisation par l'échelle spectrale locale
                scale = max(np.max(np.abs(evals)), 1.0)
                current_min_rel = evals[0] / scale
                
                if evals[0] < min_eval:
                    min_eval = evals[0]
                    min_eval_rel = current_min_rel
                    min_q = q.copy()
                    
                if current_min_rel < SEUIL_RELATIF:
                    instabilites += 1
                    
    print(f"Grille testée : {n_points}^3 = {n_points**3} points")
    print(f"Valeur propre minimale absolue : {min_eval:.4e} s^-2")
    print(f"Valeur propre minimale relative : {min_eval_rel:.4e}")
    print(f"Au point q = [{min_q[0]:.4e}, {min_q[1]:.4e}, {min_q[2]:.4e}] m^-1")
    
    if instabilites > 0:
        print(f"[FAIL] {instabilites} instabilités significatives détectées (seuil relatif < {SEUIL_RELATIF}).")
        return False
    else:
        print(f"[PASS] Aucune instabilité significative détectée sur la grille (seuil relatif < {SEUIL_RELATIF}).")
        print("(Note : Ceci est un test numérique sur une grille finie, et ne constitue pas une preuve analytique de stabilité sur tout le continuum).")
        return True

if __name__ == "__main__":
    audit_3d_grid()
