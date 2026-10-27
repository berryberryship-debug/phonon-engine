#!/usr/bin/env python3
"""
V74 — Audit de stabilité sur la zone de Brillouin 3D complète.
"""

import numpy as np
import v72_second_neighbors as engine

def audit_3d_grid(n_points=15):
    print("=== V74 : AUDIT DE STABILITÉ 3D COMPLET (SECONDS VOISINS) ===")
    PI_A = np.pi / engine.A
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    
    min_eval = np.inf
    min_q = None
    instabilites = 0
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                D = engine.dynamical_matrix(q)
                evals = np.linalg.eigvalsh(D)
                
                # On cherche la plus petite valeur propre (sans écrêtage)
                current_min = evals[0]
                if current_min < min_eval:
                    min_eval = current_min
                    min_q = q.copy()
                    
                # Tolérance relative pour le bruit numérique
                scale = max(np.max(np.abs(evals)), 1.0)
                if current_min < -1e-5 * scale:
                    instabilites += 1
                    
    print(f"Grille testée : {n_points}^3 = {n_points**3} points")
    print(f"Valeur propre minimale : {min_eval:.4e}")
    print(f"Au point q = {min_q}")
    
    if instabilites > 0:
        print(f"[FAIL] {instabilites} instabilités physiques détectées !")
    else:
        print("[PASS] Aucune instabilité physique sur la zone de Brillouin 3D.")

if __name__ == "__main__":
    audit_3d_grid()
