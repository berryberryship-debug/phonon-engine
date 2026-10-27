```python
#!/usr/bin/env python3
"""
V68 — Suivi de branches par sous-espaces, alignement de Procruste et vitesses du son.
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
    degeneracy_tolerance: float = 1.0e-8
    minimum_valid_overlap: float = 0.90
    masses: tuple[float, float] = (132.90545196, 35.45)
    output_image: str = "dispersion_v68_final.png"

def diagonalize_dynamical_matrix(q: np.ndarray, cfg: Config):
    D = engine.dynamical_matrix(q)
    eigenvalues, eigenvectors = np.linalg.eigh(D)
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).eps)
    tolerance = cfg.relative_eigenvalue_tolerance * spectral_scale

    if eigenvalues[0] < -tolerance:
        raise ValueError(f"Instabilité à q={q}: lambda={eigenvalues[0]:.6e}")

    eigenvalues = np.where(eigenvalues < 0.0, 0.0, eigenvalues)
    return np.sqrt(eigenvalues), eigenvectors, eigenvalues

def degenerate_groups(eigenvalues: np.ndarray, cfg: Config) -> list[np.ndarray]:
    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).eps)
    tolerance = cfg.degeneracy_tolerance * spectral_scale
    groups: list[np.ndarray] = []
    start = 0
    for i in range(1, len(eigenvalues) + 1):
        if i == len(eigenvalues) or eigenvalues[i] - eigenvalues[i - 1] > tolerance:
            groups.append(np.arange(start, i))
            start = i
    return groups

def subspace_projectors(eigenvectors: np.ndarray, groups: list[np.ndarray]) -> list[np.ndarray]:
    return [eigenvectors[:, group] @ eigenvectors[:, group].conj().T for group in groups]

def subspace_overlap_matrix(prev_proj: list[np.ndarray], curr_proj: list[np.ndarray]) -> np.ndarray:
    overlap = np.empty((len(prev_proj), len(curr_proj)))
    for a, P_prev in enumerate(prev_proj):
        rank_a = max(int(round(np.trace(P_prev).real)), 1)
        for b, P_curr in enumerate(curr_proj):
            rank_b = max(int(round(np.trace(P_curr).real)), 1)
            raw = np.trace(P_prev @ P_curr).real
            overlap[a, b] = raw / np.sqrt(rank_a * rank_b)
    return overlap

def align_subspace_basis(prev_vectors: np.ndarray, curr_vectors: np.ndarray) -> np.ndarray:
    """Algorithme de Procruste pour aligner les phases dans un sous-espace dégénéré."""
    # Matrice de recouvrement
    M = prev_vectors.conj().T @ curr_vectors
    # SVD : M = U S V^dagger
    U, _, Vh = np.linalg.svd(M)
    # Matrice de rotation optimale
    R = U @ Vh
    return curr_vectors @ R

def polarization_fractions(q: np.ndarray, eigenvector: np.ndarray, masses: tuple[float, float]):
    q_norm = np.linalg.norm(q)
    if q_norm < 1.0e-12:
        return {"acoustic": np.nan, "longitudinal_acoustic": np.nan, "longitudinal_optic": np.nan}
    
    q_hat = q / q_norm
    m_a, m_b = masses
    u_a, u_b = eigenvector[:3], eigenvector[3:]
    
    cm = (m_a * u_a + m_b * u_b) / (m_a + m_b)
    rel = u_a - u_b
    
    norm_cm = np.vdot(cm, cm).real
    norm_rel = np.vdot(rel, rel).real
    total = norm_cm + norm_rel
    
    if total < 1.0e-30:
        return {"acoustic": np.nan, "longitudinal_acoustic": np.nan, "longitudinal_optic": np.nan}
        
    frac_aco = norm_cm / total
    lon_aco = abs(np.vdot(q_hat, cm))**2 / max(norm_cm, 1e-30)
    lon_opt = abs(np.vdot(q_hat, rel))**2 / max(norm_rel, 1e-30)
    
    return {"acoustic": frac_aco, "longitudinal_acoustic": lon_aco, "longitudinal_optic": lon_opt}

def extract_gamma_to_x(q_path, distance, tick_pos, tick_lab):
    """Correction 1 : Utilisation des VRAIS indices numériques."""
    if "Γ" not in tick_lab or "X" not in tick_lab:
        raise ValueError("Le chemin doit contenir Γ et X.")
    
    # tick_pos contient les indices réels (0, 100, 200...)
    gamma_idx = int(tick_pos[tick_lab.index("Γ")])
    x_idx = int(tick_pos[tick_lab.index("X")])
    
    return q_path[gamma_idx : x_idx + 1], distance[gamma_idx : x_idx + 1]

def track_branches(cfg: Config):
    q_path, distance, tick_pos, tick_lab = build_q_path(cfg.n_path_points)
    q_segment, dist_segment = extract_gamma_to_x(q_path, distance, tick_pos, tick_lab)

    n_points, n_modes = len(q_segment), 6
    omega_tracked = np.full((n_points, n_modes), np.nan)
    freq_thz = np.full((n_points, n_modes), np.nan)
    vecs_tracked = np.zeros((n_points, n_modes, 6), dtype=complex)

    # Init Gamma
    omega_g, vecs_g, evals_g = diagonalize_dynamical_matrix(q_segment[0], cfg)
    groups_g = degenerate_groups(evals_g, cfg)
    proj_g = subspace_projectors(vecs_g, groups_g)
    
    omega_tracked[0] = omega_g
    freq_thz[0] = omega_g / (2 * np.pi * 1e12)
    vecs_tracked[0] = vecs_g.T

    prev_proj, prev_groups, prev_vecs = proj_g, groups_g, vecs_g

    for i in range(1, n_points):
        omega_c, vecs_c, evals_c = diagonalize_dynamical_matrix(q_segment[i], cfg)
        groups_c = degenerate_groups(evals_c, cfg)
        proj_c = subspace_projectors(vecs_c, groups_c)

        overlap = subspace_overlap_matrix(prev_proj, proj_c)
        prev_idx, curr_idx = linear_sum_assignment(overlap, maximize=True)

        for p_idx, c_idx in zip(prev_idx, curr_idx):
            p_grp, c_grp = prev_groups[p_idx], groups_c[c_idx]
            
            # Correction 2 : Alignement de Procruste
            aligned_vecs = align_subspace_basis(prev_vecs[:, p_grp], vecs_c[:, c_grp])

            for local_i, mode_i in enumerate(p_grp):
                if local_i < len(c_grp):
                    c_mode = c_grp[local_i]
                    omega_tracked[i, mode_i] = omega_c[c_mode]
                    freq_thz[i, mode_i] = omega_c[c_mode] / (2 * np.pi * 1e12)
                    vecs_tracked[i, mode_i] = aligned_vecs[:, local_i]

        prev_proj, prev_groups, prev_vecs = proj_c, groups_c, vecs_c

    return q_segment, dist_segment, freq_thz, vecs_tracked

def compute_sound_velocities(q_segment, omega_tracked, masses):
    """Correction 3 : Calcul rigoureux de la pente à Gamma."""
    # On prend le point 1 (premier non nul)
    q0, q1 = q_segment[0], q_segment[1]
    dq = np.linalg.norm(q1 - q0)
    
    if dq < 1e-12:
        return np.zeros(3)
        
    dw = omega_tracked[1, :3] - omega_tracked[0, :3]
    velocities = dw / dq
    
    print("\n=== VITESSES DU SON CALCULÉES (m/s) ===")
    print(f"Branche Acoustique 1 : {velocities[0]:.2f} m/s")
    print(f"Branche Acoustique 2 : {velocities[1]:.2f} m/s")
    print(f"Branche Acoustique 3 : {velocities[2]:.2f} m/s")
    return velocities

def main():
    cfg = Config()
    print("=== V68 : SUIVI OPTIMAL, PROCRUSTE ET VITESSES DU SON ===")
    
    q_seg, dist_seg, freq_thz, vecs = track_branches(cfg)
    
    # Vérification
    n_zero = np.sum(freq_thz[0] < 1e-6)
    print(f"Modes nuls à Gamma : {n_zero}/6 (Attendu : 3)")

    # Tracé
    fig, ax = plt.subplots(figsize=(10, 6))
    for mode in range(6):
        ax.plot(dist_seg, freq_thz[:, mode], lw=2, label=f"Branche {mode+1}")
    ax.set_title("V68 : Suivi de branche optimal (CsCl)")
    ax.set_ylabel("Fréquence (THz)")
    ax.set_xlabel("Γ → X")
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.savefig(cfg.output_image, dpi=180)
    print(f"\n[PASS] Figure : {cfg.output_image}")

    # Vitesses du son
    # (On utilise les fréquences en rad/s pour le calcul de la vitesse)
    # Attention : on n'a que les THz dans le retour, on doit recalculer en rad/s
    # Pour l'exemple, on approxime :
    velocities = compute_sound_velocities(q_seg, freq_thz * 2 * np.pi * 1e12, cfg.masses)
    
if __name__ == "__main__":
    main()
```
