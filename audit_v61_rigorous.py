#!/usr/bin/env python3
"""
Audit strict de v61 : vecteurs rigides, hermiticité et dimension du noyau.
"""

import sys
import numpy as np
import v61_phonon_3d as engine

def test_acoustic_vectors():
    print("=== TEST : VECTEURS DE TRANSLATION RIGIDE À Γ ===")
    D0 = engine.dynamical_matrix(np.zeros(3))
    norm_D = np.linalg.norm(D0, 2)
    
    # 1. Vérifications structurelles de base
    assert norm_D > 0, "D(0) est nulle ou numériquement dégénérée."
    assert np.allclose(D0, D0.conj().T), "D(0) n'est pas hermitienne."
    print(f"[PASS] D(0) hermitienne et de norme non nulle ({norm_D:.2e})")

    # 2. Vérification du noyau (dimension 3) et des modes optiques
    w = np.linalg.eigvalsh(D0)
    tol = 1e-12 * norm_D
    n_zeros = np.sum(np.abs(w) < tol)
    assert n_zeros == 3, f"Le noyau n'est pas de dimension 3 (obtenu {n_zeros})."
    assert np.all(w[3:] > tol), "La branche optique contient une valeur propre non positive."
    print("[PASS] Noyau de dimension 3 et branches optiques strictement positives.")

    # 3. Test des vecteurs de translation rigide mass-normalisés
    m1, m2 = engine.M1, engine.M2
    u_x = np.array([1, 0, 0, 1, 0, 0], dtype=float)
    u_y = np.array([0, 1, 0, 0, 1, 0], dtype=float)
    u_z = np.array([0, 0, 1, 0, 0, 1], dtype=float)
    
    e_rigid = []
    for u in [u_x, u_y, u_z]:
        e = np.zeros(6)
        e[:3] = np.sqrt(m1) * u[:3]
        e[3:] = np.sqrt(m2) * u[3:]
        e_rigid.append(e / np.linalg.norm(e))
        
    for i, e in enumerate(e_rigid):
        res = np.linalg.norm(D0 @ e)
        rel_res = res / max(norm_D, 1e-30)
        print(f"Direction {['x','y','z'][i]} : résidu relatif = {rel_res:.4e}")
        if rel_res > 1e-12:
            sys.exit("[FAIL] Le vecteur rigide n'est pas dans le noyau de D(0).")
            
    print("[PASS] Les 3 vecteurs de translation rigide sont dans le noyau à la précision machine.")

if __name__ == "__main__":
    test_acoustic_vectors()
