#!/usr/bin/env python3
"""
V75 — Anisotropie de V61 le long de [110].
Prédiction (K_PERP=2, K_PAR=10, q->0) : w^2 ∝ (3/7, 1, 11/7).
"""
import numpy as np
import v61_phonon_3d as engine

def test_direction_110():
    n = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    q_vals = np.linspace(1e6, 1e8, 10)
    r_bas, r_haut, pols = [], [], []
    ref = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    
    for qv in q_vals:
        D = engine.dynamical_matrix(qv * n)
        assert np.allclose(D, D.conj().T), "D non hermitienne."
        w2, V = np.linalg.eigh(D)
        assert np.all(np.isfinite(w2))
        a = w2[:3]
        
        r_bas.append(a[0] / a[1])
        r_haut.append(a[2] / a[1])
        
        e = V[:3, 0] + V[3:, 0]
        e = e / np.linalg.norm(e)
        pols.append(abs(e @ ref))
        
    rb, rh, pm = np.mean(r_bas), np.mean(r_haut), min(pols)
    print(f"w2(bas)/w2(milieu)  = {rb:.4f}  (attendu 3/7 = {3/7:.4f})")
    print(f"w2(haut)/w2(milieu) = {rh:.4f}  (attendu 11/7 = {11/7:.4f})")
    print(f"Recouvrement minimal du mode bas avec [1,-1,0] = {pm:.4f}  (attendu ~ 1)")
    
    ok = abs(rb - 3/7) < 5e-3 and abs(rh - 11/7) < 5e-3 and pm > 0.99
    print("[PASS] Anisotropie confirmée (Zener 7/3)." if ok else "[ECART] Prédiction non confirmée.")
    return ok

if __name__ == "__main__":
    raise SystemExit(0 if test_direction_110() else 1)
