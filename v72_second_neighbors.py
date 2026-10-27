#!/usr/bin/env python3
"""
V72 — Moteur enrichi avec les seconds voisins (Cs-Cs et Cl-Cl séparés).
Inclut les fonctions utilitaires (frequencies, build_q_path) pour remplacer V61.
"""

from itertools import product
import numpy as np

U = 1.66053906660e-27
M1 = 132.90545196 * U
M2 = 35.45 * U
A = 4.12e-10

# --- Paramètres des premiers voisins (Cs-Cl) ---
K1_PARALLEL, K1_PERP = 10.0, 2.0

# --- Paramètres des seconds voisins (Séparés) ---
K2_CS_CS_PARALLEL, K2_CS_CS_PERP = 4.0, 1.0
K2_CL_CL_PARALLEL, K2_CL_CL_PERP = 4.0, 1.0

I3 = np.eye(3)

# --- Premiers voisins ---
NEIGHBORS_1 = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    norm_R = np.linalg.norm(R)
    P = np.outer(R / norm_R, R / norm_R)
    K = K1_PERP * I3 + (K1_PARALLEL - K1_PERP) * P
    NEIGHBORS_1.append((R, K))

# --- Seconds voisins (Cs-Cs et Cl-Cl) ---
NEIGHBORS_2_CS_CS = []
NEIGHBORS_2_CL_CL = []
for axis in range(3):
    for sign in (-1, 1):
        R = np.zeros(3)
        R[axis] = sign * A
        norm_R = np.linalg.norm(R)
        P = np.outer(R / norm_R, R / norm_R)

        K_cs_cs = K2_CS_CS_PERP * I3 + (K2_CS_CS_PARALLEL - K2_CS_CS_PERP) * P
        K_cl_cl = K2_CL_CL_PERP * I3 + (K2_CL_CL_PARALLEL - K2_CL_CL_PERP) * P

        NEIGHBORS_2_CS_CS.append((R, K_cs_cs))
        NEIGHBORS_2_CL_CL.append((R, K_cl_cl))

def dynamical_matrix(q):
    q = np.asarray(q, dtype=float)
    Daa = np.zeros((3, 3), dtype=complex)
    Dbb = np.zeros((3, 3), dtype=complex)
    Dab = np.zeros((3, 3), dtype=complex)
    
    # Premiers voisins (Cs-Cl)
    for R, K in NEIGHBORS_1:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K / M1
        Dbb += K / M2
        Dab += -K * phase / np.sqrt(M1 * M2)
        
    # Seconds voisins (Cs-Cs)
    for R, K in NEIGHBORS_2_CS_CS:
        phase = np.exp(1j * np.dot(q, R))
        Daa += K * (1.0 - phase) / M1

    # Seconds voisins (Cl-Cl)
    for R, K in NEIGHBORS_2_CL_CL:
        phase = np.exp(1j * np.dot(q, R))
        Dbb += K * (1.0 - phase) / M2
        
    D = np.zeros((6, 6), dtype=complex)
    D[:3, :3] = Daa
    D[:3, 3:] = Dab
    D[3:, :3] = Dab.conj().T
    D[3:, 3:] = Dbb
    return D

def frequencies(q, tolerance=1e-12):
    D = dynamical_matrix(q)
    eigvals = np.linalg.eigvalsh(D)
    scale = max(float(np.max(np.abs(eigvals))), np.finfo(float).tiny)
    
    if eigvals[0] < -tolerance * scale:
        raise ValueError(f"Instabilité à q={q}")
        
    eigvals = np.where(
        (eigvals < 0.0) & (eigvals >= -tolerance * scale),
        0.0,
        eigvals,
    )
    return np.sqrt(eigvals)

def build_q_path(n_segment):
    PI_A = np.pi / A
    high_symmetry = [
        ("Γ", np.array([0.0, 0.0, 0.0])),
        ("X", np.array([PI_A, 0.0, 0.0])),
        ("M", np.array([PI_A, PI_A, 0.0])),
        ("R", np.array([PI_A, PI_A, PI_A])),
        ("Γ", np.array([0.0, 0.0, 0.0])),
    ]
    q_points = []
    tick_positions = []
    tick_labels = []

    for index, (label, q_point) in enumerate(high_symmetry):
        if index > 0:
            q_previous = high_symmetry[index - 1][1]
            for step in range(n_segment):
                t = step / n_segment
                q_points.append((1.0 - t) * q_previous + t * q_point)

        tick_positions.append(len(q_points))
        tick_labels.append(label)
        q_points.append(q_point)

    q_path = np.asarray(q_points, dtype=float)
    steps = np.linalg.norm(np.diff(q_path, axis=0), axis=1)
    distance = np.concatenate(([0.0], np.cumsum(steps)))

    return q_path, distance, np.asarray(tick_positions), tick_labels

# --- Audits ---
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
    
    eigenvalues = np.linalg.eigvalsh(D0)
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
    tolerance = 1.0e-12 * scale

    acoustic = eigenvalues[:3]
    optic = eigenvalues[3:]

    assert np.all(np.abs(acoustic) < tolerance), f"Modes acoustiques non nuls à Gamma : {acoustic}"
    assert np.all(optic > tolerance), f"Modes optiques non strictement positifs à Gamma : {optic}"

    print(f"[PASS] Trois modes acoustiques à Gamma : {np.array2string(acoustic, precision=3)}")
    print(f"[PASS] Trois modes optiques positifs : {np.array2string(optic, precision=3)}")
    print(f"[INFO] Gap optique / norme D : {optic[0] / norm_D:.6e}")

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
    print("(Note : la stabilité globale sur toute la zone est testée par V74).")

if __name__ == "__main__":
    test_matrix_properties()
    test_acoustic_sum_rule()
    test_eigenvalue_stability()
    print("\n=== FIN DES AUDITS V72 ===")
