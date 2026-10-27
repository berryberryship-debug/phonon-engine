#!/usr/bin/env python3
"""
Audit strict de v61 : hermiticité, noyau acoustique, gap optique,
vecteurs rigides et cohérence des masses.
"""

import sys
import numpy as np
import v61_phonon_3d as engine


def test_acoustic_vectors():
    print("=== TEST : VECTEURS DE TRANSLATION RIGIDE À Γ ===")

    q_gamma = np.zeros(3)
    D0 = engine.dynamical_matrix(q_gamma)
    norm_D = np.linalg.norm(D0, ord=2)

    assert norm_D > 0.0, "D(0) est nulle ou numériquement dégénérée."
    assert np.allclose(D0, D0.conj().T, atol=1e-12 * norm_D), "D(0) n'est pas hermitienne."
    print(f"[PASS] D(0) hermitienne, norme = {norm_D:.6e}")

    eigenvalues = np.linalg.eigvalsh(D0)
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
    tolerance = 1.0e-12 * scale

    acoustic = eigenvalues[:3]
    optic = eigenvalues[3:]

    assert np.all(np.abs(acoustic) < tolerance), f"Modes acoustiques non nuls à Gamma : {acoustic}"
    assert np.all(optic > tolerance), f"Gap optique nul ou négatif à Gamma : {optic}"

    print(f"[PASS] Trois modes acoustiques à Gamma : {np.array2string(acoustic, precision=3)}")
    print(f"[PASS] Trois modes optiques positifs : {np.array2string(optic, precision=3)}")
    print(f"[INFO] Gap optique / norme de D : {optic[0] / norm_D:.6e}")

    m1, m2 = engine.M1, engine.M2
    rigid_patterns = [
        np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0, 0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0]),
    ]

    for direction, pattern in zip(("x", "y", "z"), rigid_patterns):
        mass_weighted = np.empty(6, dtype=float)
        mass_weighted[:3] = np.sqrt(m1) * pattern[:3]
        mass_weighted[3:] = np.sqrt(m2) * pattern[3:]
        mass_weighted /= np.linalg.norm(mass_weighted)

        residual = np.linalg.norm(D0 @ mass_weighted)
        relative_residual = residual / norm_D

        print(f"Translation {direction} : résidu relatif = {relative_residual:.6e}")
        assert relative_residual < 1e-12, f"Le vecteur rigide {direction} n'appartient pas au noyau de D(0)."

    print("[PASS] Les trois translations rigides sont dans le noyau de D(0).")


def test_mass_consistency():
    print("\n=== TEST : COHÉRENCE DES MASSES ===")
    u = 1.66053906660e-27
    expected = np.array([132.90545196 * u, 35.45 * u])
    actual = np.array([engine.M1, engine.M2])
    assert np.allclose(actual, expected, rtol=1e-12), f"Masses incohérentes : moteur={actual}, attendu={expected}"
    print(f"[PASS] Masses moteur : M_Cs={engine.M1:.6e} kg, M_Cl={engine.M2:.6e} kg")


def test_hermiticity_away_from_gamma():
    print("\n=== TEST : HERMITICITÉ HORS DE Γ ===")
    for q in (np.array([0.1, 0.0, 0.0]), np.array([0.1, 0.2, 0.3]), np.array([np.pi / engine.A, 0.0, 0.0])):
        D = engine.dynamical_matrix(q)
        scale = max(float(np.linalg.norm(D, ord=2)), np.finfo(float).tiny)
        assert np.allclose(D, D.conj().T, atol=1e-12 * scale), f"D(q) n'est pas hermitienne à q={q}."
    print("[PASS] D(q) est hermitienne sur les points testés.")


def test_optical_modes_positive_on_gamma_x():
    print("\n=== TEST : STABILITÉ SUR Γ-X ===")
    n_points = 21
    q_gamma = np.zeros(3)
    q_x = np.array([np.pi / engine.A, 0.0, 0.0])
    for i in range(n_points):
        t = i / (n_points - 1)
        q = (1.0 - t) * q_gamma + t * q_x
        frequencies = engine.frequencies(q)
        assert np.all(np.isfinite(frequencies)), f"Fréquence non finie à q={q}."
        assert np.all(frequencies >= 0.0), f"Fréquence négative à q={q}."
    print(f"[PASS] {n_points} points de Γ à X sont dynamiquement stables.")


if __name__ == "__main__":
    test_acoustic_vectors()
    test_mass_consistency()
    test_hermiticity_away_from_gamma()
    test_optical_modes_positive_on_gamma_x()
    print("\n[PASS] Audit V61 terminé.")
