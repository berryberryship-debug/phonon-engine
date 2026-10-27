#!/usr/bin/env python3
"""
V75 — Anisotropie de V61 le long de [110].
Identification dynamique des modes par polarisation et extrapolation à q=0.
"""

import numpy as np
import v61_phonon_3d as engine

def get_cm(V, idx):
    """Dé-massifie le vecteur propre et renvoie le centre de masse physique."""
    m1, m2 = engine.M1, engine.M2
    u1 = V[:3, idx] / np.sqrt(m1)
    u2 = V[3:, idx] / np.sqrt(m2)
    cm = (m1 * u1 + m2 * u2) / (m1 + m2)
    norm = np.linalg.norm(cm)
    return cm / norm if norm > 1e-14 else np.zeros(3)

def test_convergence_110():
    print("=== TEST [110] : ANISOTROPIE ET CONVERGENCE ===")
    
    # Directions de polarisation théoriques à q -> 0
    p_ta1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    p_ta2 = np.array([0.0, 0.0, 1.0])
    p_la  = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    
    n = p_la # Direction de propagation
    
    # On teste uniquement sur les petits q pour la limite asymptotique
    q_vals = np.logspace(5, 6.5, 6) 
    
    r_ta1, r_ta2, r_la = [], [], []

    for qv in q_vals:
        D = engine.dynamical_matrix(qv * n)
        w2, V = np.linalg.eigh(D)
        
        # On isole les 3 modes acoustiques (les plus bas)
        idx_aco = np.argsort(w2)[:3]
        
        # On calcule la projection de chaque mode sur les polarisations théoriques
        proj_ta1, proj_ta2, proj_la = [], [], []
        
        for i in idx_aco:
            cm = get_cm(V, i)
            proj_ta1.append(abs(np.dot(cm, p_ta1)))
            proj_ta2.append(abs(np.dot(cm, p_ta2)))
            proj_la.append(abs(np.dot(cm, p_la)))
            
        # Identification des modes par recouvrement maximal
        i_ta1 = idx_aco[np.argmax(proj_ta1)]
        i_ta2 = idx_aco[np.argmax(proj_ta2)]
        i_la  = idx_aco[np.argmax(proj_la)]
        
        # Vérification que l'on a bien 3 modes distincts
        if len({i_ta1, i_ta2, i_la}) != 3:
            raise AssertionError(f"Ambiguïté de polarisation à q={qv}.")
            
        # Calcul des rapports par rapport au mode TA2 (le mode "milieu")
        r_ta1.append(w2[i_ta1] / w2[i_ta2])
        r_la.append(w2[i_la] / w2[i_ta2])

    # Extrapolation à q=0 (q^2 normalisé pour le conditionnement)
    q2 = np.array(q_vals)**2
    q2_norm = q2 / np.max(q2)
    
    # Régression linéaire sur les 6 points
    coef_ta1 = np.polyfit(q2_norm, r_ta1, 1)
    coef_la  = np.polyfit(q2_norm, r_la, 1)
    
    r_ta1_0 = coef_ta1[1]
    r_la_0  = coef_la[1]
    
    print(f"q_max utilisé pour l'extrapolation : {q_vals[-1]:.1e} m^-1")
    print(f"r(TA1/TA2) à q -> 0 : {r_ta1_0:.5f}  (attendu 3/7 = {3/7:.5f})")
    print(f"r(LA/TA2)  à q -> 0 : {r_la_0:.5f}  (attendu 11/7 = {11/7:.5f})")
    
    # Contrôle
    ok_ta1 = abs(r_ta1_0 - 3/7) < 2e-3
    ok_la = abs(r_la_0 - 11/7) < 2e-3
    
    if ok_ta1 and ok_la:
        print("[PASS] Limite acoustique confirmée (Zener 7/3).")
        return True
    else:
        print("[FAIL] La limite acoustique ne correspond pas à la prédiction.")
        return False

if __name__ == "__main__":
    raise SystemExit(0 if test_convergence_110() else 1)
