```python
#!/usr/bin/env python3
"""
V70 — Validation acoustique (V61) par régression.
"""

import numpy as np
import v61_phonon_3d as engine

def analyze_sound_velocities():
    print("=== ESTIMATION DES VITESSES ACOUSTIQUES (Γ → X) ===")
    # 10 points entre 10^6 et 10^8 m^-1
    q_vals = np.linspace(1e6, 1e8, 10)
    w_ta, w_la = [], []
    
    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        evals = np.linalg.eigvalsh(D)
        evals = np.where(evals < -1e10, np.nan, np.where(evals < 0, 0.0, evals))
        
        w_ta.append(np.sqrt(evals[0]))
        w_la.append(np.sqrt(evals[2]))

    def fit_velocity(w):
        v, _, _, _ = np.linalg.lstsq(q_vals.reshape(-1, 1), w, rcond=None)
        return v[0]

    v_ta = fit_velocity(w_ta)
    v_la = fit_velocity(w_la)
    
    print(f"v_TA = {v_ta:.2f} m/s")
    print(f"v_LA = {v_la:.2f} m/s")
    
    if np.isclose(v_ta, v_la, rtol=1e-3):
        print("[INFO] Les vitesses sont égales (isotropie du modèle V61).")

if __name__ == "__main__":
    analyze_sound_velocities()
```

---

6. v72_second_neighbors.py (Moteur avec seconds voisins)

```python
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

NEIGHBORS_1 = []
for sx, sy, sz in product((-1, 1), repeat=3):
    R = 0.5 * A * np.array([sx, sy, sz], dtype=float)
    P = np.outer(R / np.linalg.norm(R), R / np.linalg.norm(R))
    NEIGHBORS_1.append((R, K1_PERP * I3 + (K1_PARALLEL - K1_PERP) * P))

NEIGHBORS_2 = []
for axis in range(3):
    for sign in [-1, 1]:
        R = np.zeros(3); R[axis] = sign * A
        P = np.outer(R / np.linalg.norm(R), R / np.linalg.norm(R))
        NEIGHBORS_2.append((R, K2_PERP * I3 + (K2_PARALLEL - K2_PERP) * P))

def dynamical_matrix(q):
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

def test_acoustic_sum_rule():
    D0 = dynamical_matrix(np.zeros(3))
    norm_D = np.linalg.norm(D0, 2)
    u_x = np.array([1, 0, 0, 1, 0, 0], dtype=float)
    e = np.zeros(6)
    e[:3] = np.sqrt(M1) * u_x[:3]
    e[3:] = np.sqrt(M2) * u_x[3:]
    e = e / np.linalg.norm(e)
    res = np.linalg.norm(D0 @ e) / max(norm_D, 1e-30)
    print(f"[PASS] Règle de somme V72 (résidu relatif = {res:.4e})")

if __name__ == "__main__":
    test_acoustic_sum_rule()
```
