#!/usr/bin/env python3
"""
V69 — Suivi 6x6, Procruste persistant, et vitesses du son par régression.
"""

from __future__ import annotations
import sys
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

try:
    import v61_phonon_3d as engine
    from v61_phonon_3d import build_q_path
except ImportError as exc:
    raise SystemExit("Erreur : v61_phonon_3d.py introuvable.")

@dataclass(frozen=True)
class Config:
    n_path_points: int = 100
    relative_eigenvalue_tolerance: float = 1.0e-10
    masses: tuple[float, float] = (132.90545196, 35.45)
    output_image: str = "dispersion_v69_final.png"

def diagonalize_dynamical_matrix(q: np.ndarray, cfg: Config):
    D = engine.dynamical_matrix(q)
    eigenvalues, eigenvectors = np.linalg.eigh(D)
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).eps)
    tolerance = cfg.relative_eigenvalue_tolerance * spectral_scale

    if eigenvalues[0] < -tolerance:
        raise ValueError(f"Instabilité à q={q}: lambda={eigenvalues[0]:.6e}")

    eigenvalues = np.where(eigenvalues < 0.0, 0.0, eigenvalues)
    return np.sqrt(eigenvalues), eigenvectors

def align_procruste(prev_vecs: np.ndarray, curr_vecs: np.ndarray) -> np.ndarray:
    """Aligne la base de curr_vecs sur prev_vecs via SVD."""
    M = prev_vecs.conj().T @ curr_vecs
    U, _, Vh = np.linalg.svd(M)
    R = U @ Vh
    return curr_vecs @ R

def extract_gamma_to_x(q_path, distance, tick_pos, tick_lab):
    if "Γ" not in tick_lab or "X" not in tick_lab:
        raise ValueError("Le chemin doit contenir Γ et X.")
    gamma_idx = int(tick_pos[tick_lab.index("Γ")])
    x_idx = int(tick_pos[tick_lab.index("X")])
    return q_path[gamma_idx : x_idx + 1], distance[gamma_idx : x_idx + 1]

def track_branches_6x6(cfg: Config):
    q_path, distance, tick_pos, tick_lab = build_q_path(cfg.n_path_points)
    q_segment, dist_segment = extract_gamma_to_x(q_path, distance, tick_pos, tick_lab)

    n_points, n_modes = len(q_segment), 6
    omega_tracked = np.zeros((n_points, n_modes))
    vecs_tracked = np.zeros((n_points, n_modes, 6), dtype=complex)

    # Init Gamma
    omega_g, vecs_g = diagonalize_dynamical_matrix(q_segment[0], cfg)
    omega_tracked[0] = omega_g
    vecs_tracked[0] = vecs_g.T

    prev_vecs = vecs_g

    for i in range(1, n_points):
        omega_c, vecs_c = diagonalize_dynamical_matrix(q_segment[i], cfg)
        
        # 1. Matrice de recouvrement 6x6 complète
        overlap = np.abs(prev_vecs.conj().T @ vecs_c)**2
        
        # 2. Algorithme Hongrois (maximise le recouvrement global)
        _, col_ind = linear_sum_assignment(overlap, maximize=True)
        
        # 3. Réordonnancement des vecteurs courants
        vecs_c_sorted = vecs_c[:, col_ind]
        
        # 4. Alignement de Procruste
        vecs_c_aligned = align_procruste(prev_vecs, vecs_c_sorted)
        
        # 5. Enregistrement
        omega_tracked[i] = omega_c[col_ind]
        vecs_tracked[i] = vecs_c_aligned.T
        
        # CRUCIAL : Mise à jour avec les vecteurs alignés pour l'itération suivante
        prev_vecs = vecs_c_aligned

    return q_segment, dist_segment, omega_tracked, vecs_tracked

def classify_polarization(q, vec, masses):
    """Détermine si le mode est L ou T en regardant le centre de masse."""
    q_norm = np.linalg.norm(q)
    if q_norm < 1e-12: return "N/A"
    
    q_hat = q / q_norm
    m_a, m_b = masses
    u_a, u_b = vec[:3], vec[3:]
    u_cm = (m_a * u_a + m_b * u_b) / (m_a + m_b)
    
    # Projection sur la direction de propagation
    dot = abs(np.dot(u_cm, q_hat)) / (np.linalg.norm(u_cm) + 1e-15)
    return "L" if dot > 0.7 else "T"

def compute_velocities(q_seg, omega, masses):
    """Calcule v_s par régression linéaire sur les 10 premiers points."""
    n_pts = 10
    q_norm = np.linalg.norm(q_seg[:n_pts], axis=1)
    
    print("\n=== VITESSES DU SON (RÉGRESSION LINÉAIRE) ===")
    for mode in range(3): # Les 3 branches acoustiques
        # Conversion rad/s
        w_rads = omega[:n_pts, mode] * 2 * np.pi * 1e12
        
        # Régression linéaire sans ordonnée à l'origine (w = v * q)
        v_s, _, _, _ = np.linalg.lstsq(q_norm.reshape(-1, 1), w_rads, rcond=None)
        v_s = v_s[0]
        
        # Polarisation au point milieu
        mid_idx = n_pts // 2
        pol = classify_polarization(q_seg[mid_idx], 
                                    np.real(np.linalg.eigh(engine.dynamical_matrix(q_seg[mid_idx]))[1])[:, mode], 
                                    masses)
        
        print(f"Mode {mode+1} ({pol}) : {v_s:.2f} m/s")

def main():
    cfg = Config()
    print("=== V69 : SUIVI OPTIMAL 6x6 ET VITESSES RIGOUREUSES ===")
    
    q_seg, dist_seg, freq_thz, vecs = track_branches_6x6(cfg)
    
    n_zero = np.sum(freq_thz[0] < 1e-6)
    print(f"Modes nuls à Gamma : {n_zero}/6")

    # Tracé
    fig, ax = plt.subplots(figsize=(10, 6))
    for mode in range(6):
        ax.plot(dist_seg, freq_thz[:, mode], lw=2, label=f"Branche {mode+1}")
    ax.set_title("V69 : Suivi 6x6 avec Procruste persistant")
    ax.set_ylabel("Fréquence (THz)")
    ax.set_xlabel("Γ → X")
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.savefig(cfg.output_image, dpi=180)
    print(f"[PASS] Figure : {cfg.output_image}")

    compute_velocities(q_seg, freq_thz, cfg.masses)
    
if __name__ == "__main__":
    main()
