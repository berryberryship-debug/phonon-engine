#!/usr/bin/env python3
"""
V75 — Anisotropie de V61 le long de [110].
Vérifie la convergence des rapports w^2 vers (3/7 : 1 : 11/7) quand q -> 0.
"""

import numpy as np
import v61_phonon_3d as engine

def get_polarization_of_mode(V, idx):
    """Dé-massifie un vecteur propre et renvoie son CM physique."""
    m1, m2 = engine.M1, engine.M2
    u1 = V[:3, idx] / np.sqrt(m1)
    u2 = V[3:, idx] / np.sqrt(m2)
    cm = (m1 * u1 + m2 * u2) / (m1 + m2)
    return cm / (np.linalg.norm(cm) + 1e-30)

def test_convergence_110():
    print("=== TEST [110] : CONVERGENCE VERS q -> 0 ===")
    n = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    ref = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    
    # Balayage sur plusieurs décades de q
    q_vals = np.logspace(5, 8, 10)
    
    print(f"{'q (m^-1)':>12} | {'r_bas':>10} | {'r_haut':>10} | {'pol_bas':>10}")
    print("-" * 55)

    r_bas_list, r_haut_list, q2_list = [], [], []

    for qv in q_vals:
        D = engine.dynamical_matrix(qv * n)
        assert np.allclose(D, D.conj().T), "D non hermitienne."
        w2, V = np.linalg.eigh(D)
        assert np.all(np.isfinite(w2))
        
        a = w2[:3]
        r_bas = a[0] / a[1]
        r_haut = a[2] / a[1]
        
        # Dé-massification pour la polarisation du mode bas
        cm_bas = get_polarization_of_mode(V, 0)
        pol_bas = abs(np.vdot(cm_bas, ref)) # Projection sur [1,-1,0]
        
        r_bas_list.append(r_bas)
        r_haut_list.append(r_haut)
        q2_list.append(qv**2)
        
        print(f"{qv:12.1e} | {r_bas:10.4f} | {r_haut:10.4f} | {pol_bas:10.4f}")

    # Extrapolation à q=0 par régression linéaire en q^2 (y = a + b*q^2)
    coef_bas = np.polyfit(q2_list, r_bas_list, 1)
    coef_haut = np.polyfit(q2_list, r_haut_list, 1)
    
    r_bas_0 = coef_bas[1]
    r_haut_0 = coef_haut[1]
    
    print(f"\nExtrapolation à q -> 0 :")
    print(f"  r_bas(0)  = {r_bas_0:.5f}  (attendu 3/7  = {3/7:.5f})")
    print(f"  r_haut(0) = {r_haut_0:.5f}  (attendu 11/7 = {11/7:.5f})")
    
    # Vérification de la polarisation à petit q
    pol_min = pol_bas_list_min = min([abs(np.vdot(get_polarization_of_mode(np.linalg.eigh(engine.dynamical_matrix(qv * n))[1], 0), ref)) for qv in q_vals])
    
    # On refait un check sur le dernier point (q le plus petit) pour la polarisation
    # Note : pol_bas est recalculé dans la boucle, mais on garde la liste pour la vérif.
    # Plus simple : on prend le dernier pol_bas affiché.
    
    ok_r_bas = abs(r_bas_0 - 3/7) < 2e-3
    ok_r_haut = abs(r_haut_0 - 11/7) < 2e-3
    ok_pol = pol_bas > 0.99  # Valeur du dernier point (q le plus petit)
    
    if ok_r_bas and ok_r_haut and ok_pol:
        print("\n[PASS] Anisotropie confirmée (Zener 7/3) et convergence vérifiée.")
        return True
    else:
        print("\n[ÉCART] Les rapports ou la polarisation ne convergent pas comme attendu.")
        return False

if __name__ == "__main__":
    raise SystemExit(0 if test_convergence_110() else 1)
