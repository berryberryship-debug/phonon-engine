```python
#!/usr/bin/env python3
"""
V67 — Suivi de branches phononiques par sous-espaces dégénérés.
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
    n_path_points: int = 100
    relative_eigenvalue_tolerance: float = 1.0e-10
    degeneracy_tolerance: float = 1.0e-8
    minimum_valid_overlap: float = 0.90
    masses: tuple[float, float] = (132.90545196, 35.45)  # Masses réelles CsCl
    output_image: str = "dispersion_v67_subspace_tracking.png"
    output_csv: str = "dispersion_v67_subspace_tracking.csv"

# ---------------------------------------------------------------------------
# Diagonalisation et diagnostic spectral
# ---------------------------------------------------------------------------
def diagonalize_dynamical_matrix(q: np.ndarray, cfg: Config):
    D = engine.dynamical_matrix(q)
    eigenvalues, eigenvectors = np.linalg.eigh(D)

    spectral_scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).eps)
    tolerance = cfg.relative_eigenvalue_tolerance * spectral_scale

    if eigenvalues[0] < -tolerance:
        raise ValueError(
            f"Instabilité dynamique à q={q}: "
            f"lambda_min={eigenvalues[0]:.6e}, tolerance={tolerance:.6e}"
        )

    eigenvalues = np.where(eigenvalues < 0.0, 0.0, eigenvalues)
    omega = np.sqrt(eigenvalues)
    return omega, eigenvectors, eigenvalues

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

# ---------------------------------------------------------------------------
# Suivi par sous-espaces
# ---------------------------------------------------------------------------
def subspace_projectors(eigenvectors: np.ndarray, groups: list[np.ndarray]) -> list[np.ndarray]:
    projectors = []
    for group in groups:
        vectors = eigenvectors[:, group]
        projectors.append(vectors @ vectors.conj().T)
    return projectors

def subspace_overlap_matrix(previous_projectors: list[np.ndarray], current_projectors: list[np.ndarray]) -> np.ndarray:
    overlap = np.empty((len(previous_projectors), len(current_projectors)))
    for a, previous_projection in enumerate(previous_projectors):
        rank_a = max(int(round(np.trace(previous_projection).real)), 1)
        for b, current_projection in enumerate(current_projectors):
            rank_b = max(int(round(np.trace(current_projection).real)), 1)
            raw = np.trace(previous_projection @ current_projection).real
            overlap[a, b] = raw / np.sqrt(rank_a * rank_b)
    return overlap

def assign_subspaces(overlap: np.ndarray):
    previous_indices, current_indices = linear_sum_assignment(overlap, maximize=True)
    return previous_indices, current_indices

# ---------------------------------------------------------------------------
# Classification physique
# ---------------------------------------------------------------------------
def polarization_fractions(q: np.ndarray, eigenvector: np.ndarray, masses: tuple[float, float]) -> dict[str, float]:
    q_norm = np.linalg.norm(q)
    if q_norm < 1.0e-12:
        return {k: np.nan for k in ["acoustic", "optic", "longitudinal_acoustic", "transverse_acoustic", "longitudinal_optic", "transverse_optic"]}

    q_hat = q / q_norm
    m_a, m_b = masses
    total_mass = m_a + m_b
    u_a = eigenvector[:3]
    u_b = eigenvector[3:]

    center_of_mass = (m_a * u_a + m_b * u_b) / total_mass
    relative = u_a - u_b

    acoustic_norm2 = float(np.vdot(center_of_mass, center_of_mass).real)
    optic_norm2 = float(np.vdot(relative, relative).real)
    total_norm2 = acoustic_norm2 + optic_norm2

    if total_norm2 < 1.0e-30:
        return {k: np.nan for k in ["acoustic", "optic", "longitudinal_acoustic", "transverse_acoustic", "longitudinal_optic", "transverse_optic"]}

    acoustic_fraction = acoustic_norm2 / total_norm2
    optic_fraction = optic_norm2 / total_norm2
    longitudinal_acoustic = abs(np.vdot(q_hat, center_of_mass)) ** 2 / max(acoustic_norm2, 1.0e-30)
    longitudinal_optic = abs(np.vdot(q_hat, relative)) ** 2 / max(optic_norm2, 1.0e-30)

    return {
        "acoustic": acoustic_fraction,
        "optic": optic_fraction,
        "longitudinal_acoustic": longitudinal_acoustic,
        "transverse_acoustic": 1.0 - longitudinal_acoustic,
        "longitudinal_optic": longitudinal_optic,
        "transverse_optic": 1.0 - longitudinal_optic,
    }

def diagnostic_label(fractions: dict[str, float]) -> str:
    if np.isnan(fractions["acoustic"]): return "Gamma"
    acoustic = fractions["acoustic"]
    if acoustic >= 0.75: character = "A"
    elif acoustic <= 0.25: character = "O"
    else: character = "M"

    if character == "A": longitudinal = fractions["longitudinal_acoustic"]
    elif character == "O": longitudinal = fractions["longitudinal_optic"]
    else: return "M"

    if longitudinal >= 0.75: return f"L{character}"
    if longitudinal <= 0.25: return f"T{character}"
    return f"M{character}"

# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------
def extract_gamma_to_x(q_path, distance, tick_pos, tick_lab, cfg):
    if "Γ" not in tick_lab or "X" not in tick_lab:
        raise ValueError("Le chemin doit contenir les points Gamma et X.")

    gamma_index = int(tick_lab.index("Γ"))
    x_index = int(tick_lab.index("X"))
    if x_index <= gamma_index:
        raise ValueError("X doit apparaître après Gamma dans le chemin.")

    q_segment = q_path[gamma_index : x_index + 1]
    distance_segment = distance[gamma_index : x_index + 1]
    return q_segment, distance_segment

def track_branches(cfg: Config):
    q_path, distance, tick_pos, tick_lab = build_q_path(cfg.n_path_points)
    q_segment, distance_segment = extract_gamma_to_x(q_path, distance, tick_pos, tick_lab, cfg)

    n_points = len(q_segment)
    n_modes = 6

    omega_tracked = np.full((n_points, n_modes), np.nan)
    frequency_thz = np.full((n_points, n_modes), np.nan)
    labels = np.empty((n_points, n_modes), dtype=object)
    acoustic_fraction = np.full((n_points, n_modes), np.nan)
    longitudinal_fraction = np.full((n_points, n_modes), np.nan)
    overlap_history = np.full((n_points, n_modes), np.nan)

    omega_gamma, vectors_gamma, eigenvalues_gamma = diagonalize_dynamical_matrix(q_segment[0], cfg)
    groups_gamma = degenerate_groups(eigenvalues_gamma, cfg)
    projectors_gamma = subspace_projectors(vectors_gamma, groups_gamma)

    omega_tracked[0, :] = omega_gamma
    frequency_thz[0, :] = omega_gamma / (2.0 * np.pi * 1.0e12)

    for mode in range(n_modes):
        labels[0, mode] = "Gamma"
        fractions = polarization_fractions(q_segment[0], vectors_gamma[:, mode], cfg.masses)
        acoustic_fraction[0, mode] = fractions["acoustic"]
        longitudinal_fraction[0, mode] = fractions["longitudinal_acoustic"]

    previous_projectors = projectors_gamma
    previous_groups = groups_gamma

    for i in range(1, n_points):
        omega_current, vectors_current, eigenvalues_current = diagonalize_dynamical_matrix(q_segment[i], cfg)
        groups_current = degenerate_groups(eigenvalues_current, cfg)
        projectors_current = subspace_projectors(vectors_current, groups_current)

        overlap = subspace_overlap_matrix(previous_projectors, projectors_current)
        previous_indices, current_indices = assign_subspaces(overlap)

        for previous_group_index, current_group_index in zip(previous_indices, current_indices):
            previous_group = previous_groups[previous_group_index]
            current_group = groups_current[current_group_index]
            overlap_value = overlap[previous_group_index, current_group_index]

            if overlap_value < cfg.minimum_valid_overlap:
                print(f"[WARN] q index {i}: sous-espace ambigu (overlap={overlap_value:.4f}).")

            for local_index, mode_index in enumerate(previous_group):
                if local_index < len(current_group):
                    current_mode = current_group[local_index]
                    omega_tracked[i, mode_index] = omega_current[current_mode]
                    frequency_thz[i, mode_index] = omega_current[current_mode] / (2.0 * np.pi * 1.0e12)
                    overlap_history[i, mode_index] = overlap_value

                    fractions = polarization_fractions(q_segment[i], vectors_current[:, current_mode], cfg.masses)
                    acoustic_fraction[i, mode_index] = fractions["acoustic"]
                    labels[i, mode_index] = diagnostic_label(fractions)

                    if fractions["acoustic"] >= 0.75:
                        longitudinal_fraction[i, mode_index] = fractions["longitudinal_acoustic"]
                    elif fractions["acoustic"] <= 0.25:
                        longitudinal_fraction[i, mode_index] = fractions["longitudinal_optic"]

        previous_projectors = projectors_current
        previous_groups = groups_current

    return {
        "q_segment": q_segment,
        "distance": distance_segment,
        "frequency_thz": frequency_thz,
        "labels": labels,
        "acoustic_fraction": acoustic_fraction,
        "longitudinal_fraction": longitudinal_fraction,
        "overlap_history": overlap_history,
    }

def plot_results(results: dict, cfg: Config):
    distance = results["distance"]
    frequency_thz = results["frequency_thz"]

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0.0, 0.9, frequency_thz.shape[1]))

    for mode in range(frequency_thz.shape[1]):
        valid = np.isfinite(frequency_thz[:, mode])
        ax.plot(distance[valid], frequency_thz[valid, mode], color=colors[mode], lw=2.0, label=f"Branche suivie {mode + 1}")

    ax.axvline(distance[-1], color="gray", ls="--", lw=1.0)
    ax.text(distance[-1], ax.get_ylim()[1] * 0.97, "X", ha="right", va="top", color="gray")

    ax.set_xlabel("Distance réduite le long de Γ–X")
    ax.set_ylabel("Fréquence (THz)")
    ax.set_title("V67 — Suivi modal par sous-espaces dégénérés (CsCl)")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=9)

    fig.tight_layout()
    fig.savefig(cfg.output_image, dpi=180)
    plt.close(fig)

def run_validation(results: dict, cfg: Config):
    frequency_thz = results["frequency_thz"]
    overlap_history = results["overlap_history"]

    print("\n=== Validation numérique V67 ===")
    gamma_frequencies = frequency_thz[0]
    n_zero_modes = int(np.sum(gamma_frequencies < 1.0e-6))
    print(f"Modes à fréquence quasi nulle à Gamma : {n_zero_modes}/6")

    if n_zero_modes != 3:
        print("[FAIL] Un cristal stable à deux atomes par maille doit avoir trois modes acoustiques.")
    else:
        print("[PASS] Trois modes acoustiques présents à Gamma.")

    finite_overlaps = overlap_history[np.isfinite(overlap_history)]
    if finite_overlaps.size:
        print(f"Chevauchement de sous-espace : min={finite_overlaps.min():.6f}, médiane={np.median(finite_overlaps):.6f}")
        if finite_overlaps.min() < cfg.minimum_valid_overlap:
            print("[WARN] Au moins un suivi de sous-espace est ambigu.")
        else:
            print("[PASS] Tous les chevauchements dépassent le seuil déclaré.")

def main():
    cfg = Config()
    print("=== V67 : suivi modal par sous-espaces dégénérés ===")
    results = track_branches(cfg)
    plot_results(results, cfg)
    run_validation(results, cfg)
    print(f"\n[PASS] Figure : {cfg.output_image}")

if __name__ == "__main__":
    main()
```
