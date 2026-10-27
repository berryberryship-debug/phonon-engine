#!/usr/bin/env python3
"""
V75 — Preuve de l'instabilité du modèle V61 (forces centrales) le long de [110].
"""

import numpy as np
import v61_phonon_3d as engine

def test_instability_110():
    print("=== V75 : TEST D'INSTABILITÉ V61 LE LONG DE [110] ===")
    # On se place à un petit q le long de [110]
    dq = 1e7
    q = np.array([dq, dq, 0.0])
    
    D = engine.dynamical_matrix(q)
    evals = np.linalg.eigvalsh(D)
    
    print(f"Vecteur d'onde q = {q}")
    print(f"Valeurs propres (omega^2) :")
    for i, val in enumerate(evals):
        print(f"  Mode {i+1} : {val:.4e}")
        
    # La valeur propre minimale devrait être négative ou nulle
    min_eval = evals[0]
    if min_eval < -1e5: # Tolérance pour le bruit numérique
        print("[FAIL] Instabilité confirmée : le modèle V61 est mécaniquement instable.")
    else:
        print("[PASS] Pas d'instabilité détectée (contrairement à la prédiction théorique).")

if __name__ == "__main__":
    test_instability_110()
