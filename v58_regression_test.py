#!/usr/bin/env python3
import csv
import sys
from pathlib import Path
import numpy as np

CSV_PATH = Path("phonon_dispersion_501.csv")
PNG_PATH = Path("phonon_curves_501.png")
SYMMETRY_RTOL = 1e-12
DERIVATIVE_RTOL = 2e-4
FREQUENCY_AT_ZERO_ATOL = 1e-6

def load_data(path):
    if not path.is_file(): raise FileNotFoundError(f"CSV introuvable : {path}")
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"q_m^-1", "omega_acoustic_rad_s^-1", "omega_optical_rad_s^-1", "vg_aco_ms", "vg_opt_ms"}
        rows = list(reader)
    data = {name: np.asarray([float(row[name]) for row in rows], dtype=float) for name in required}
    return data

def run_tests(data):
    q = data["q_m^-1"]; wa = data["omega_acoustic_rad_s^-1"]; wo = data["omega_optical_rad_s^-1"]
    va = data["vg_aco_ms"]; vo = data["vg_opt_ms"]
    results = []
    def check(name, condition, detail=""):
        results.append(bool(condition))
        print(f"[{'PASS' if condition else 'FAIL'}] {name} {detail}")

    check("Grille strictement croissante", np.all(np.diff(q) > 0))
    check("Grille symetrique", np.allclose(q, -q[::-1], rtol=SYMMETRY_RTOL, atol=1e-8))
    check("Point central q=0", q[len(q) // 2] == 0.0)
    check("Frequences finies", np.isfinite(wa).all() and np.isfinite(wo).all())
    check("Frequences non negatives", np.all(wa >= 0) and np.all(wo >= 0))
    check("Mode acoustique nul a q=0", abs(wa[len(q)//2]) <= FREQUENCY_AT_ZERO_ATOL)
    check("Gap optique positif", np.min(wo) > 0.0)
    check("Symetrie acoustique", np.allclose(wa, wa[::-1], rtol=SYMMETRY_RTOL, atol=1e-6))
    check("Symetrie optique", np.allclose(wo, wo[::-1], rtol=SYMMETRY_RTOL, atol=1e-6))

    pos = q >= 0
    check("Acoustique croissante pour q>=0", np.all(np.diff(wa[pos]) >= -1e-12 * np.max(wa)))
    check("Optique decroissante pour q>=0", np.all(np.diff(wo[pos]) <= 1e-12 * np.max(wo)))

    dwa = np.gradient(wa, q, edge_order=2); dwo = np.gradient(wo, q, edge_order=2)
    mask = np.ones(len(q), dtype=bool); mask[[0, 1, len(q)//2, len(q)-2, len(q)-1]] = False

    for name, analytic, numeric in (("Vitesse acoustique", va, dwa), ("Vitesse optique", vo, dwo)):
        scale = max(float(np.max(np.abs(numeric[mask]))), 1.0)
        rel = float(np.max(np.abs(analytic[mask] - numeric[mask])) / scale)
        check(f"Derivee numerique : {name}", rel < DERIVATIVE_RTOL, f"(err={rel:.3e})")

    check("Graphique PNG present", PNG_PATH.is_file() and PNG_PATH.stat().st_size > 0)
    passed = sum(results)
    print(f"\nResultat V58 : {passed}/{len(results)} tests reussis")
    return passed == len(results)

if __name__ == "__main__":
    try:
        data = load_data(CSV_PATH)
    except Exception as exc:
        print(f"[ERREUR] {exc}"); sys.exit(2)
    sys.exit(0 if run_tests(data) else 1)
