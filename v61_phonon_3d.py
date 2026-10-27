# Dans v61_phonon_3d.py, remplace la fonction frequencies par :

def frequencies(q, tolerance=1e-12):  # <-- 1e-7 remplacé par 1e-12
    D = dynamical_matrix(q)
    eigvals = np.linalg.eigvalsh(D)
    scale = max(float(np.max(np.abs(eigvals))), np.finfo(float).tiny)
    if eigvals[0] < -tolerance * scale:
        raise ValueError(f"Instabilité à q={q}")
    if np.all(q == 0.0):
        eigvals[np.abs(eigvals) < tolerance * scale] = 0.0
    else:
        eigvals[(eigvals < 0.0) & (eigvals >= -tolerance * scale)] = 0.0
    return np.sqrt(eigvals)
