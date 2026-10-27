#!/usr/bin/env python3
"""
V72 — Moteur enrichi avec les seconds voisins (Cs-Cs et Cl-Cl).
"""

from itertools import product
import numpy as np

U = 1.66053906660e-27
M1 = 132.90545196 * U
M2 = 35.45 * U
A = 4.12e-10

K1_PARALLEL, K1_PERP = 10.0, 2.0
K2_PARALLEL, K2_PERP = 4.0, 1.0
I3 = np.eye(3)

# --- Premiers voisins (8 atomes de Cl autour de Cs) ---
NEIGHBORS_1 = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    norm_R = np.linalg.norm(R)
    P = np.outer(R / norm_R, R / norm_R)
    K = K1_PERP * I3 + (K1_PARALLEL - K1_PERP) * P
    NEIGHBORS_1.append((R, K))

# --- Seconds voisins (6 atomes le long des axes) ---
NEIGHBORS_2 = []
for axis in range(3):
    for sign in [-1, 1]:
        R = np.zeros(3)
        R[axis] = sign * A
        norm_R = np.linalg.norm(R)
        P = np.outer(R / norm_R, R / norm_R)
        K = K2_PERP * I3 + (K2_PARALLEL - K2_PERP) * P
        NEIGHBORS_2.append((R, K))

def dynamical_matrix(q):
    q = np.asarray(q, dtype=float)
    Daa = np.zeros((3, 3), dtype=complex)
    Dbb = np.zeros((3, 3), dtype=complex)
    Dab = np.zeros((3, 3), dtype=complex)
    
    for R, K in NEIGHBORS_1:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K / M1
        Dbb += K / M2
        Dab += -K * phase / np.sqrt(M1 * M2)
        
    for R, K in NEIGHBORS_2:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K * (1.0 - phase) / M1
        Dbb += K * (1.0 - phase) / M2
        
    D = np.zeros((6, 6), dtype=complex)
    D[:3, :3] = Daa
    D[:3, 3:] = Dab
    D[3:, :3] = Dab.conj().T
    D[3:, 3:] = Dbb
    return D

def test_matrix_properties():
    print("=== TEST 1 : PROPRIÉTÉS DE LA MATRICE ===")
    q_test = np.array([0.1, 0.2, 0.3]) * 1e9
    D = dynamical_matrix(q_test)
    
    if not np.allclose(D, D.conj().T):
        raise AssertionError("La matrice globale D n'est pas hermitienne.")
    print("[PASS] La matrice globale D est hermitienne.")
    
    Daa = D[:3, :3]
    Dbb = D[3:, 3:]
    if not np.allclose(Daa, Daa.conj().T):
        raise AssertionError("Le bloc Daa n'est pas hermitien.")
    if not np.allclose(Dbb, Dbb.conj().T):
        raise AssertionError("Le bloc Dbb n'est pas hermitien.")
    print("[PASS] Les blocs diagonaux Daa et Dbb sont hermitiens.")
    
    Dab = D[:3, 3:]
    Dba = D[3:, :3]
    if not np.allclose(Dba, Dab.conj().T):
        raise AssertionError("La relation Dba = Dab^dagger n'est pas respectée.")
    print("[PASS] La relation Dba = Dab^dagger est respectée.")

def test_acoustic_sum_rule():
    print("\n=== TEST 2 : RÈGLE DE SOMME ACOUSTIQUE À Γ ===")
    D0 = dynamical_matrix(np.zeros(3))
    norm_D = np.linalg.norm(D0, 2)
    
    if norm_D <= 0:
        raise AssertionError("La matrice D(0) est nulle.")
    print(f"[PASS] Norme de D(0) non nulle ({norm_D:.2e})")
    
    u_x = np.array([1, 0, 0, 1, 0, 0], dtype=float)
    u_y = np.array([0, 1, 0, 0, 1, 0], dtype=float)
    u_z = np.array([0, 0, 1, 0, 0, 1], dtype=float)
    
    e_rigid = []
    for u in [u_x, u_y, u_z]:
        e = np.zeros(6)
        e[:3] = np.sqrt(M1) * u[:3]
        e[3:] = np.sqrt(M2) * u[3:]
        e_rigid.append(e / np.linalg.norm(e))
        
    for i, e in enumerate(e_rigid):
        res = np.linalg.norm(D0 @ e)
        rel_res = res / norm_D
        print(f"Direction {['x','y','z'][i]} : résidu relatif = {rel_res:.4e}")
        if rel_res > 1e-12:
            raise AssertionError(f"Le vecteur rigide selon {['x','y','z'][i]} n'est pas dans le noyau.")
            
    print("[PASS] Les 3 vecteurs de translation rigide sont dans le noyau.")
    
    # Vérification de la dimension du noyau
    evals = np.linalg.eigvalsh(D0)
    n_zero = np.sum(np.abs(evals) < 1e-12 * norm_D)
    if n_zero != 3:
        raise AssertionError(f"Le noyau n'est pas de dimension 3 (obtenu {n_zero}).")
    print("[PASS] Le noyau est de dimension 3.")

def test_eigenvalue_stability():
    print("\n=== TEST 3 : STABILITÉ SUR QUELQUES POINTS ===")
    q_points = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1e9, 0.0, 0.0]),
        np.array([1e9, 1e9, 0.0]),
        np.array([1e9, 1e9, 1e9])
    ]
    
    for q in q_points:
        D = dynamical_matrix(q)
        evals = np.linalg.eigvalsh(D)
        scale = np.max(np.abs(evals))
        
        if scale == 0:
            raise AssertionError(f"Spectre nul à q={q}.")
            
        if evals[0] < -1e-10 * scale:
            raise AssertionError(f"Instabilité détectée à q={q} : {evals[0]:.4e}")
        print(f"q={q} : eigenvalue min = {evals[0]:.4e} (stable)")
        
    print("[PASS] Aucune instabilité détectée sur les points testés.")
    print("(Note : la stabilité globale sur toute la zone de Brillouin est testée par V74).")

if __name__ == "__main__":
    test_matrix_properties()
    test_acoustic_sum_rule()
    test_eigenvalue_stability()
    print("\n=== FIN DES AUDITS V72 ===")
