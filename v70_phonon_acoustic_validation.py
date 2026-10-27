#!/usr/bin/env python3
"""
V70 — Traitement rigoureux des déplacements physiques, permutation hongroise exacte,
continuité de jauge par sous-espaces stricts et régression acoustique multi-intervalles.
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
    raise SystemExit(
        "Erreur : v61_phonon_3d.py doit être importable depuis le répertoire courant.\n"
        f"Détail : {exc}"
    )

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Config:
    n_path_points: int = 120
    relative_eigenvalue_tolerance: float = 1.0e-10
    degeneracy_tolerance: float = 1.0e-6
    # Paramètre de maille du CsCl en mètres
    lattice_constant: float = 4.123e-10
    # Masses atomiques réelles en uma (Cs, Cl)
    masses_uma: tuple[float, float] = (132.90545196, 35.45)
    output_image: str = "dispersion_v70_validation.png"

# ---------------------------------------------------------------------------
# Inspection du moteur & conversion du vecteur d'onde
# ---------------------------------------------------------------------------
def get_q_scale_factor(q_path: np.ndarray, tick_pos: list[int], tick_lab: list[str], cfg: Config) -> float:
    """
    Détermine si q est en coordonnées réduites ou en SI.
    Au point X, le vecteur réduit standard CsCl (cubique simple) est (0.5, 0, 0)
    ou (pi/a, 0, 0) si exprimé en rad/m.
    """
    x_idx = int(tick_pos[tick_lab.index("X")])
    q_x = q_path[x_idx]
    norm_x = np.linalg.norm(q_x)

    # Si norm_x est de l'ordre de 0.5 ou 1.0, ce sont des coordonnées réduites
    if norm_x < 2.0:
        # q_SI = (2 * pi / a) * q_red
        return 2.0 * np.pi / cfg.lattice_constant
    # Si norm_x est de l'ordre de pi/a (~ 7.6e9 m^-1), le moteur est déjà en SI
    return 1.0

# ---------------------------------------------------------------------------
# Diagonalisation
# ---------------------------------------------------------------------------
def diagonalize_dynamical_matrix(q: np.ndarray, cfg: Config):
    D = engine.dynamical_matrix(q)
    eigenvalues, eigenvectors = np.linalg.eigh(D)

    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).eps)
    tolerance = cfg.relative_eigenvalue_tolerance * spectral_scale

    if eigenvalues[0] < -tolerance:
        raise ValueError(
            f"Instabilité dynamique critique à q={q}: "
            f"lambda_min={eigenvalues[0]:.6e}, tolerance={tolerance:.6e}"
        )

    # Écrêtage contrôlé des résidus d'arrondi
    eigenvalues_clean = np.where(eigenvalues < 0.0, 0.0, eigenvalues)
    omega_rads = np.sqrt(eigenvalues_clean)
    freq_thz = omega_rads / (2.0 * np.pi * 1.0e12)
    return freq_thz, omega_rads, eigenvectors, eigenvalues_clean

# ---------------------------------------------------------------------------
# Continuité de phase et sous-espaces stricts
# ---------------------------------------------------------------------------
def identify_degenerate_blocks(evals: np.ndarray, tol: float) -> list[list[int]]:
    """Identifie les blocs par composantes connexes sur les valeurs propres."""
    n = len(evals)
    visited = [False] * n
    blocks = []
    for i in range(n):
        if visited[i]:
            continue
        group = [i]
        visited[i] = True
        for j in range(i + 1, n):
            if not visited[j] and abs(evals[j] - evals[i]) <= tol:
                group.append(j)
                visited[j] = True
        blocks.append(group)
    return blocks

def align_eigenvectors(prev_vecs: np.ndarray, curr_vecs: np.ndarray, 
                       evals: np.ndarray, cfg: Config) -> np.ndarray:
    """
    Régule la jauge sans mélanger des sous-espaces de valeurs propres distinctes :
      - SVD Procruste pour les groupes strictement dégénérés.
      - Facteur de phase U(1) pour les modes isolés.
    """
    aligned = curr_vecs.copy()
    scale = max(float(np.max(np.abs(evals))), np.finfo(float).eps)
    tol = cfg.degeneracy_tolerance * scale
    blocks = identify_degenerate_blocks(evals, tol)

    for block in blocks:
        if len(block) == 1:
            m = block[0]
            ov = np.vdot(prev_vecs[:, m], aligned[:, m])
            if abs(ov) > 1.0e-12:
                aligned[:, m] *= (ov / abs(ov))
        else:
            sub_p = prev_vecs[:, block]
            sub_c = aligned[:, block]
            M = sub_p.conj().T @ sub_c
            U, _, Vh = np.linalg.svd(M)
            R = U @ Vh
            aligned[:, block] = sub_c @ R.conj().T

    return aligned

# ---------------------------------------------------------------------------
# Déplacements physiques & Classification
# ---------------------------------------------------------------------------
def physical_displacements(eigenvector: np.ndarray, masses_uma: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """
    Convertit le vecteur propre mass-normalisé e en déplacements réels u = e / sqrt(M).
    e a 6 composantes : [e_Cs_x, e_Cs_y, e_Cs_z, e_Cl_x, e_Cl_y, e_Cl_z].
    """
    m_cs, m_cl = masses_uma
    u_cs = eigenvector[:3] / np.sqrt(m_cs)
    u_cl = eigenvector[3:] / np.sqrt(m_cl)
    return u_cs, u_cl

def classify_branch(q: np.ndarray, eigenvector: np.ndarray, masses_uma: tuple[float, float]) -> str:
    """
    Caractérise la polarisation physique (Acoustique/Optique, L/T)
    via le mouvement du centre de masse et le déplacement relatif.
    """
    q_norm = np.linalg.norm(q)
    if q_norm < 1.0e-12:
        return "Γ"

    q_hat = q / q_norm
    m_cs, m_cl = masses_uma
    u_cs, u_cl = physical_displacements(eigenvector, masses_uma)

    # Déplacement du centre de masse et déplacement relatif
    u_cm = (m_cs * u_cs + m_cl * u_cl) / (m_cs + m_cl)
    u_rel = u_cs - u_cl

    norm_cm2 = float(np.vdot(u_cm, u_cm).real)
    norm_rel2 = float(np.vdot(u_rel, u_rel).real)
    total_energy = m_cs * np.linalg.norm(u_cs)**2 + m_cl * np.linalg.norm(u_cl)**2

    # Ratio d'énergie cinétique dans le centre de masse
    cm_energy_ratio = (m_cs + m_cl) * norm_cm2 / max(total_energy, 1.0e-30)

    is_acoustic = cm_energy_ratio >= 0.50
    char = "A" if is_acoustic else "O"

    u_ref = u_cm if is_acoustic else u_rel
    ref_norm = np.linalg.norm(u_ref)
    if ref_norm < 1.0e-15:
        return f"?{char}"

    # Projection longitudinale
    longitudinal_ratio = abs(np.vdot(q_hat, u_ref))**2 / (ref_norm**2)

    if longitudinal_ratio >= 0.70:
        return f"L{char}"
    elif longitudinal_ratio <= 0.30:
        return f"T{char}"
    return f"M{char}"

# ---------------------------------------------------------------------------
# Suivi modal principal
# ---------------------------------------------------------------------------
def track_modes(cfg: Config):
    q_path, distance, tick_pos, tick_lab = build_q_path(cfg.n_path_points)
    
    gamma_idx = int(tick_pos[tick_lab.index("Γ")])
    x_idx = int(tick_pos[tick_lab.index("X")])
    q_segment = q_path[gamma_idx : x_idx + 1]
    dist_segment = distance[gamma_idx : x_idx + 1]

    n_points, n_modes = len(q_segment), 6
    freq_tracked = np.zeros((n_points, n_modes))
    omega_tracked = np.zeros((n_points, n_modes))
    vecs_tracked = np.zeros((n_points, n_modes, 6), dtype=complex)
    labels = np.empty((n_points, n_modes), dtype=object)

    # Initialisation à Gamma
    freq_g, omega_g, vecs_g, evals_g = diagonalize_dynamical_matrix(q_segment[0], cfg)
    freq_tracked[0] = freq_g
    omega_tracked[0] = omega_g
    vecs_tracked[0] = vecs_g.T
    for m in range(n_modes):
        labels[0, m] = "Γ"

    prev_vecs = vecs_g

    for i in range(1, n_points):
        freq_c, omega_c, vecs_c, evals_c = diagonalize_dynamical_matrix(q_segment[i], cfg)

        # Matrice de recouvrement
        overlap = np.abs(prev_vecs.conj().T @ vecs_c)**2

        # Assignation de Kuhn-Munkres
        row_ind, col_ind = linear_sum_assignment(overlap, maximize=True)

        # Construction rigoureuse de la permutation
        permutation = np.empty(n_modes, dtype=int)
        permutation[row_ind] = col_ind

        vecs_sorted = vecs_c[:, permutation]
        omega_sorted = omega_c[permutation]
        freq_sorted = freq_c[permutation]
        evals_sorted = evals_c[permutation]

        # Alignement de jauge
        vecs_aligned = align_eigenvectors(prev_vecs, vecs_sorted, evals_sorted, cfg)

        freq_tracked[i] = freq_sorted
        omega_tracked[i] = omega_sorted
        vecs_tracked[i] = vecs_aligned.T

        for m in range(n_modes):
            labels[i, m] = classify_branch(q_segment[i], vecs_aligned[:, m], cfg.masses_uma)

        prev_vecs = vecs_aligned

    q_scale = get_q_scale_factor(q_path, tick_pos, tick_lab, cfg)

    return {
        "q_segment": q_segment,
        "dist_segment": dist_segment,
        "q_scale": q_scale,
        "freq_thz": freq_tracked,
        "omega": omega_tracked,
        "vecs": vecs_tracked,
        "labels": labels,
    }

# ---------------------------------------------------------------------------
# Analyse des vitesses acoustiques
# ---------------------------------------------------------------------------
def analyze_sound_velocities(data: dict, cfg: Config):
    q_seg = data["q_segment"]
    omega = data["omega"]
    labels = data["labels"]
    q_scale = data["q_scale"]

    q_norm_si = np.linalg.norm(q_seg, axis=1) * q_scale

    print("\n=== ESTIMATION DES VITESSES ACOUSTIQUES (Γ → X) ===")
    print(f"Facteur d'échelle q (-> SI m^-1) : {q_scale:.4e}")

    # Recherche des branches acoustiques au point 1
    acoustic_modes = [m for m in range(6) if labels[1, m] in ("LA", "TA", "MA")]

    if len(acoustic_modes) != 3:
        # Repli : prend les 3 plus basses fréquences
        acoustic_modes = list(np.argsort(omega[1])[:3])

    for m in acoustic_modes:
        tag = labels[1, m]
        # Test de robustesse sur plusieurs fenêtres d'ajustement
        results = []
        for n_pts in [4, 6, 8, 10]:
            x = q_norm_si[:n_pts, np.newaxis]
            y = omega[:n_pts, m]
            v, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
            results.append(float(v[0]))

        v_mean = np.mean(results)
        v_std = np.std(results)
        print(f"Mode {m + 1} ({tag:>2s}) : v_s = {v_mean:7.2f} ± {v_std:4.2f} m/s "
              f"(fenêtres 4-10 pts: {[round(r, 1) for r in results]})")

# ---------------------------------------------------------------------------
# Tracé et validation
# ---------------------------------------------------------------------------
def plot_results(data: dict, cfg: Config):
    dist = data["dist_segment"]
    freq = data["freq_thz"]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for m in range(6):
        ax.plot(dist, freq[:, m], lw=2.0, label=f"Mode {m + 1} ({data['labels'][-1, m]})")

    ax.set_title("V70 — Dispersion de phonons (CsCl) : Γ → X", fontsize=12)
    ax.set_ylabel("Fréquence (THz)", fontsize=11)
    ax.set_xlabel("Vecteur d'onde réduit", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, loc="upper left")

    fig.tight_layout()
    plt.savefig(cfg.output_image, dpi=200)
    print(f"\n[PASS] Figure enregistrée : {cfg.output_image}")

def main():
    cfg = Config()
    print("=== V70 : VALIDATION ACOUSTIQUE & DÉPLACEMENTS PHYSIQUES ===")
    data = track_modes(cfg)
    plot_results(data, cfg)
    analyze_sound_velocities(data, cfg)

if __name__ == "__main__":
    main()
