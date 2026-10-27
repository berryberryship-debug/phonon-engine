#!/usr/bin/env python3
"""
V75 — Test de l'anisotropie de V61 le long de [110].
"""

import numpy as np
import v61_phonon_3d as engine

def test_direction_110():
    print("=== TEST [110] : MODE TRANSVERSE NUL ? ===")
    n = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    q_vals = np.linspace(1e6, 1e8, 10)
    r_nul, r_la = [], []
    pols = []

    for qv in q_vals:
        D = engine.dynamical_matrix(qv * n)
        assert np.allclose(D, D.conj().T), "D non hermitienne."
        w2, V = np.linalg.eigh(D)
        assert np.all(np.isfinite(w2))
        a = w2[:3]                      
        r_nul.append(a[0] / a[1])       
        r_la.append(a[2] / a[1])        

        e = V[:3, 0] + V[3:, 0]
        e = e / np.linalg.norm(e)
        pol = abs(e @ np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0))
        pols.append(pol)

    print(f"w²(bas)/w²(milieu) = {np.mean(r_nul):.3e}  (attendu ~ 0)")
    print(f"w²(haut)/w²(milieu) = {np.mean(r_la):.4f}  (attendu ~ 2)")
    print(f"Recouvrement min du mode bas avec [1,-1,0] = {min(pols):.4f}  (attendu ~ 1)")

    if abs(np.mean(r_nul)) < 1e-3 and abs(np.mean(r_la) - 2) < 1e-2:
        print("[PASS] Mode transverse [1,-1,0] nul : V61 est fortement anisotrope.")
    else:
        print("[ÉCART] La prédiction n'est pas confirmée.")

if __name__ == "__main__":
    test_direction_110()
