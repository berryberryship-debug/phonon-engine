```python
#!/usr/bin/env python3
"""
V73 — Validation rigoureuse du modèle à seconds voisins.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import v72_second_neighbors as engine

def compute_velocities():
    print("=== V73 : VITESSES ACOUSTIQUES (RÉGRESSION MULTI-Q) ===")
    q_vals = np.linspace(1e5, 1e7, 10)
    w_ta1, w_ta2, w_la = [], [], []

    for q_val in q_vals:
        q = np.array([q_val, 0.0, 0.0])
        D = engine.dynamical_matrix(q)
        evals = np.linalg.eigvalsh(D)
        evals = np.where(evals < -1e10, np.nan, np.where(evals < 0, 0.0, evals))
        w_ta1.append(np.sqrt(evals[0]))
        w_ta2.append(np.sqrt(evals[1]))
        w_la.append(np.sqrt(evals[2]))

    def fit_velocity(w):
        v, _, _, _ = np.linalg.lstsq(q_vals.reshape(-1, 1), w, rcond=None)
        return v[0]

    v_ta1 = fit_velocity(w_ta1)
    v_ta2 = fit_velocity(w_ta2)
    v_la = fit_velocity(w_la)

    print(f"v_TA1 = {v_ta1:.2f} m/s")
    print(f"v_TA2 = {v_ta2:.2f} m/s")
    print(f"v_LA  = {v_la:.2f} m/s")
    print(f"Écart relatif LA/TA = {v_la / v_ta1:.4f}")

if __name__ == "__main__":
    compute_velocities()
```

---

8. v74_full_bz_stability.py (Stabilité 3D)

```python
#!/usr/bin/env python3
"""
V74 — Audit de stabilité sur la zone de Brillouin 3D complète.
"""

import numpy as np
import v72_second_neighbors as engine

def audit_3d_grid(n_points=15):
    print("=== V74 : AUDIT DE STABILITÉ 3D ===")
    PI_A = np.pi / engine.A
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    
    min_eval = np.inf
    instabilites = 0
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                D = engine.dynamical_matrix(q)
                evals = np.linalg.eigvalsh(D)
                current_min = evals[0]
                if current_min < min_eval:
                    min_eval = current_min
                scale = max(np.max(np.abs(evals)), 1.0)
                if current_min < -1e-5 * scale:
                    instabilites += 1
                    
    print(f"Grille : {n_points}^3 = {n_points**3} points")
    print(f"Valeur propre min : {min_eval:.4e}")
    
    if instabilites > 0:
        print(f"[FAIL] {instabilites} instabilités détectées !")
    else:
        print("[PASS] Aucune instabilité physique détectée.")

if __name__ == "__main__":
    audit_3d_grid()
```
