#!/usr/bin/env python3
"""
V73 — Vitesses acoustiques de V72 le long de Γ → X.
"""

import numpy as np
import v72_second_neighbors as engine


def rigid_translation(direction: np.ndarray, m1: float, m2: float) -> np.ndarray:
    vector = np.zeros(6, dtype=float)
    vector[:3] = direction * np.sqrt(m1)
    vector[3:] = direction * np.sqrt(m2)
    return vector / np.linalg.norm(vector)


def acoustic_modes_along_x(q_value: float):
    q = np.array([q_value, 0.0, 0.0])
    D = engine.dynamical_matrix(q)

    eigenvalues, eigenvectors = np.linalg.eigh(D)
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)

    # Tolérance relative à l'échelle globale de la matrice (évite le faux positif à Gamma)
    stability_tolerance = 1.0e-10 * scale

    if eigenvalues[0] < -stability_tolerance:
        raise ValueError(
            f"Instabilité détectée à q={q}: "
            f"lambda_min={eigenvalues[0]:.6e}, "
            f"tolérance={stability_tolerance:.6e}"
        )

    frequencies = np.sqrt(np.clip(eigenvalues, 0.0, None))
    acoustic_indices = np.argsort(eigenvalues)[:3]

    q_hat = np.array([1.0, 0.0, 0.0])
    longitudinal_fractions = []

    for index in acoustic_indices:
        vector = eigenvectors[:, index]

        # Dé-massification pour obtenir les vrais déplacements
        displacement = np.empty(6, dtype=float)
        displacement[:3] = vector[:3] / np.sqrt(engine.M1)
        displacement[3:] = vector[3:] / np.sqrt(engine.M2)

        center_of_mass = (
            engine.M1 * displacement[:3] + engine.M2 * displacement[3:]
        ) / (engine.M1 + engine.M2)

        longitudinal_fraction = (
            abs(np.vdot(q_hat, center_of_mass)) ** 2
            / max(np.vdot(center_of_mass, center_of_mass).real, 1.0e-30)
        )
        longitudinal_fractions.append(longitudinal_fraction)

    # Le mode LA est celui qui a la plus grande fraction longitudinale
    la_idx_in_aco = int(np.argmax(longitudinal_fractions))
    la_index = acoustic_indices[la_idx_in_aco]
    ta_indices = [index for index in acoustic_indices if index != la_index]

    # Vérification que le mode identifié est bien longitudinal (fraction > 0.9)
    if longitudinal_fractions[la_idx_in_aco] < 0.90:
        raise ValueError(
            f"Le mode identifié comme LA n'est pas assez longitudinal "
            f"(fraction={longitudinal_fractions[la_idx_in_aco]:.2f}) à q={q_value}"
        )

    return {
        "q": q_value,
        "frequency_la": frequencies[la_index],
        "frequency_ta_degenerate": float(np.mean(frequencies[ta_indices])),
    }


def slope_through_origin(q: np.ndarray, omega: np.ndarray):
    q = np.asarray(q, dtype=float)
    omega = np.asarray(omega, dtype=float)

    slope = float(np.sum(q * omega) / np.sum(q * q))
    residuals = omega - slope * q

    n = len(q)
    degrees_of_freedom = max(n - 1, 1)
    residual_variance = float(np.sum(residuals**2) / degrees_of_freedom)
    slope_error = float(np.sqrt(residual_variance / np.sum(q * q)))

    relative_residual = float(
        np.max(np.abs(residuals) / np.maximum(np.abs(omega), np.finfo(float).tiny))
    )

    return slope, slope_error, relative_residual


def compute_velocities():
    print("=== V73 : VITESSES ACOUSTIQUES DE V72 LE LONG DE Γ–X ===")

    q_max_values = np.array([1.0e7, 2.5e7, 5.0e7, 1.0e8])
    q_reference = np.linspace(1.0e6, 1.0e8, 40)

    la_frequencies = []
    ta_frequencies = []

    for q_value in q_reference:
        result = acoustic_modes_along_x(q_value)
        la_frequencies.append(result["frequency_la"])
        ta_frequencies.append(result["frequency_ta_degenerate"])

    la_frequencies = np.asarray(la_frequencies)
    ta_frequencies = np.asarray(ta_frequencies)

    print("\n--- Pente globale sur [1e6, 1e8] m^-1 ---")

    v_la, v_la_error, la_residual = slope_through_origin(q_reference, la_frequencies)
    v_ta, v_ta_error, ta_residual = slope_through_origin(q_reference, ta_frequencies)

    print(f"v_LA = {v_la:.3f} ± {v_la_error:.3f} m/s")
    print(f"v_TA = {v_ta:.3f} ± {v_ta_error:.3f} m/s")
    print(f"Résidu relatif maximal LA : {la_residual:.3e}")
    print(f"Résidu relatif maximal TA : {ta_residual:.3e}")

    print("\n--- Convergence de la pente quand q_max diminue ---")
    print("q_max (m^-1) | v_LA (m/s) | v_TA (m/s) | résidu LA | résidu TA")

    for q_max in q_max_values:
        mask = q_reference <= q_max

        v_la_local, _, la_residual_local = slope_through_origin(
            q_reference[mask], la_frequencies[mask]
        )
        v_ta_local, _, ta_residual_local = slope_through_origin(
            q_reference[mask], ta_frequencies[mask]
        )

        print(
            f"{q_max:12.3e} | "
            f"{v_la_local:10.3f} | "
            f"{v_ta_local:10.3f} | "
            f"{la_residual_local:9.3e} | "
            f"{ta_residual_local:9.3e}"
        )

    anisotropy = abs(v_la - v_ta) / v_ta
    print(f"\nAnisotropie relative LA/TA : {anisotropy:.3%}")

if __name__ == "__main__":
    compute_velocities()
