import numpy as np

def generer_courbe_dispersion_fictive(n_points=100):
    q = np.linspace(-np.pi, np.pi, n_points)
    freq_acoustiques = np.abs(np.sin(q/2)) * 2.0 
    freq_optiques = 5.0 + np.cos(q) 
    dispersion = np.column_stack([
        freq_acoustiques * 0.8, freq_acoustiques * 1.0, freq_acoustiques * 1.2,
        freq_optiques * 0.9, freq_optiques * 1.0, freq_optiques * 1.1
    ])
    return dispersion

def encodage_positionnel(n_tokens, d_model):
    positions = np.arange(n_tokens)[:, np.newaxis]
    div_term = np.exp(np.arange(0, d_model, 2) * -(np.log(10000.0) / d_model))
    pe = np.zeros((n_tokens, d_model))
    pe[:, 0::2] = np.sin(positions * div_term)
    pe[:, 1::2] = np.cos(positions * div_term)
    return pe

def preparer_tenseur_pour_ia(dispersion, d_model=64):
    n_tokens, n_branches = dispersion.shape
    moyenne = np.mean(dispersion, axis=0)
    ecart_type = np.std(dispersion, axis=0) + 1e-8
    dispersion_normalisee = (dispersion - moyenne) / ecart_type
    poids_projection = np.random.randn(n_branches, d_model) * 0.1
    embedding = dispersion_normalisee @ poids_projection
    position_encoding = encodage_positionnel(n_tokens, d_model)
    tenseur_final = embedding + position_encoding
    return np.expand_dims(tenseur_final, axis=0)
