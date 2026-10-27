#!/usr/bin/env python3
"""
V62 (Corrigé) — Audit de stabilité sur une grille 3D incluant explicitement Gamma.
"""

import numpy as np
import sys

try:
    from v61_phonon_3d import dynamical_matrix, PI_A
except ImportError:
    print("Erreur : v61_phonon_3d.py introuvable.")
    sys.exit(1)

def audit_grille_3d(n_points=21):
    if n_points % 2 == 0:
        n_points += 1
        
    print(f"=== V62 : AUDIT DE STABILITÉ SUR GRILLE 3D ({n_points}x{n_points}x{n_points}) ===")
    print(f"Nombre total de points testés : {n_points**3}")
    
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    
    min_eigenvalue = np.inf
    min_q = None
    instabilites = 0
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                
                D = dynamical_matrix(q)
                eigvals = np.linalg.eigvalsh(D)
                
                current_min = np.min(eigvals)
                if current_min < min_eigenvalue:
                    min_eigenvalue = current_min
                    min_q = q.copy()
                    
                tolerance = -1e-5 * max(np.max(np.abs(eigvals)), 1.0)
                if current_min < tolerance:
                    instabilites += 1

    print(f"\nValeur propre minimale trouvée : {min_eigenvalue:.4e} s^-2")
    print(f"Au point q = [{min_q[0]:.4e}, {min_q[1]:.4e}, {min_q[2]:.4e}] m^-1")
    
    D_gamma = dynamical_matrix(np.zeros(3))
    eig_gamma = np.linalg.eigvalsh(D_gamma)
    n_zero_gamma = np.count_nonzero(np.abs(eig_gamma) < 1e-10 * max(np.max(np.abs(eig_gamma)), 1.0))
    print(f"\nModes nuls à Gamma (q=0) : {n_zero_gamma} (attendu : 3)")
    
    if instabilites > 0:
        print(f"\n[FAIL] {instabilites} points instables détectés !")
        sys.exit(1)
    else:
        print("\n[PASS] Aucune instabilité sur la grille.")
        print("Le modèle est dynamiquement stable sur le volume testé.")

if __name__ == "__main__":
    try:
        audit_grille_3d(n_points=21)
    except Exception as e:
        print(f"ÉCHEC V62 : {e}", file=sys.stderr)
        sys.exit(1)
