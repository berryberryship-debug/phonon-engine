```python
#!/usr/bin/env python3
"""
V62 — Audit de stabilité sur une grille 3D (inclut Gamma).
"""
import numpy as np
import sys
from v61_phonon_3d import dynamical_matrix, PI_A

def audit_grille(n_points=21):
    print(f"=== V62 : AUDIT STABILITÉ ({n_points}x{n_points}x{n_points}) ===")
    q_vals = np.linspace(-PI_A, PI_A, n_points)
    min_eval, min_q, instabilites = np.inf, None, 0
    
    for qx in q_vals:
        for qy in q_vals:
            for qz in q_vals:
                q = np.array([qx, qy, qz])
                eigvals = np.linalg.eigvalsh(dynamical_matrix(q))
                
                current_min = np.min(eigvals)
                if current_min < min_eval:
                    min_eval, min_q = current_min, q.copy()
                    
                if current_min < -1e-5 * max(np.max(np.abs(eigvals)), 1.0):
                    instabilites += 1

    print(f"Min eigenvalue : {min_eval:.4e} s^-2 à q={min_q}")
    
    D_gamma = dynamical_matrix(np.zeros(3))
    n_zero = np.count_nonzero(np.abs(np.linalg.eigvalsh(D_gamma)) < 1e-10)
    print(f"Modes nuls à Gamma : {n_zero} (attendu : 3)")
    
    if instabilites > 0:
        print("[FAIL] Instabilités détectées.")
        sys.exit(1)
    print("[PASS] Stable sur la grille.")

if __name__ == "__main__":
    audit_grille(21)
```
