import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# 1. DONNÉES
# ============================================================

def generer_dispersion(n_points: int = 50) -> np.ndarray:
    """Génère une courbe de dispersion fictive (2 acoustiques + 2 optiques)."""
    q = np.linspace(-np.pi, np.pi, n_points)

    omega = np.stack([
        2.5 * np.abs(np.sin(q / 2)),          # LA
        1.8 * np.abs(np.sin(q / 2)),          # TA
        4.2 + 0.30 * np.cos(q),               # LO
        3.6 + 0.25 * np.cos(q),               # TO
    ], axis=1)

    omega += 0.02 * np.random.randn(*omega.shape)
    return np.maximum(omega, 0.05)

def tokenizer(dispersion: np.ndarray, d_model: int = 32) -> np.ndarray:
    """Transforme la dispersion en séquence d'embeddings (N, d_model)."""
    N = dispersion.shape[0]
    q = np.linspace(-np.pi, np.pi, N)

    features = [
        (q / np.pi)[:, None],
        dispersion / (dispersion.max() + 1e-8),
    ]

    for k in range(1, 5):
        features += [np.sin(k * q)[:, None], np.cos(k * q)[:, None]]

    grad = np.gradient(dispersion, axis=0)
    features.append(grad / (np.abs(grad).max() + 1e-8))

    X = np.concatenate(features, axis=1)

    rng = np.random.RandomState(42)
    W = rng.randn(X.shape[1], d_model) * 0.5
    T = X @ W
    return (T - T.mean(0)) / (T.std(0) + 1e-8)

# ============================================================
# 2. MINI-TRANSFORMER + ADAM
# ============================================================

class MiniTransformer:
    def __init__(self, d_model: int = 32):
        self.d_model = d_model
        scale = 1.0 / np.sqrt(d_model)

        self.W_Q = np.random.randn(d_model, d_model) * scale
        self.W_K = np.random.randn(d_model, d_model) * scale
        self.W_V = np.random.randn(d_model, d_model) * scale
        self.W_O = np.random.randn(d_model, d_model) * scale

        # États Adam
        self.m = {name: np.zeros_like(getattr(self, name)) for name in ["W_Q", "W_K", "W_V", "W_O"]}
        self.v = {name: np.zeros_like(getattr(self, name)) for name in ["W_Q", "W_K", "W_V", "W_O"]}
        self.t = 0
        self.beta1, self.beta2, self.eps = 0.9, 0.999, 1e-8

    def forward(self, X: np.ndarray):
        self.X = X
        Q = X @ self.W_Q
        K = X @ self.W_K
        V = X @ self.W_V

        scores = (Q @ K.T) / np.sqrt(self.d_model)
        scores -= scores.max(axis=-1, keepdims=True)          # stabilité numérique
        A = np.exp(scores)
        self.A = A / A.sum(axis=-1, keepdims=True)

        context = self.A @ V
        self.Z = context @ self.W_O
        return self.A, self.Z

    def _adam_update(self, name: str, grad: np.ndarray, lr: float):
        self.t += 1
        m = self.m[name] = self.beta1 * self.m[name] + (1 - self.beta1) * grad
        v = self.v[name] = self.beta2 * self.v[name] + (1 - self.beta2) * grad**2

        m_hat = m / (1 - self.beta1 ** self.t)
        v_hat = v / (1 - self.beta2 ** self.t)

        param = getattr(self, name)
        setattr(self, name, param - lr * m_hat / (np.sqrt(v_hat) + self.eps))

    def backward(self, dZ: np.ndarray, lr: float = 0.005):
        # Gradient à travers W_O
        d_context = dZ @ self.W_O.T
        dW_O = (self.A @ (self.X @ self.W_V)).T @ dZ

        # Gradient attention
        V = self.X @ self.W_V
        dA = d_context @ V.T
        dV = self.A.T @ d_context

        # Jacobien du softmax
        dS = self.A * (dA - (dA * self.A).sum(axis=-1, keepdims=True))

        dQ = (dS @ (self.X @ self.W_K)) / np.sqrt(self.d_model)
        dK = (dS.T @ (self.X @ self.W_Q)) / np.sqrt(self.d_model)

        dW_Q = self.X.T @ dQ
        dW_K = self.X.T @ dK
        dW_V = self.X.T @ dV

        # Mises à jour Adam
        self._adam_update("W_Q", dW_Q, lr)
        self._adam_update("W_K", dW_K, lr)
        self._adam_update("W_V", dW_V, lr)
        self._adam_update("W_O", dW_O, lr)

# ============================================================
# 3. ENTRAÎNEMENT & ANALYSE
# ============================================================

def main():
    np.random.seed(42)

    N, d_model = 50, 32
    dispersion = generer_dispersion(N)
    T = tokenizer(dispersion, d_model)

    idx_mask = N // 2          # q ≈ 0
    model = MiniTransformer(d_model)

    losses, entropies = [], []
    n_iters, lr = 1000, 0.005

    print(f"🚀 Entraînement auto-encodeur ({N} tokens, d_model={d_model})")
    print(f"   Masque à l'index {idx_mask} (q = 0)\n")

    for i in range(n_iters):
        X_masked = T.copy()
        X_masked[idx_mask] = 0.0

        A, Z = model.forward(X_masked)

        # Perte sur toute la séquence (force l'attention globale)
        loss = np.mean((Z - T) ** 2)
        losses.append(loss)

        dZ = 2.0 * (Z - T) / (N * d_model)
        model.backward(dZ, lr=lr)

        # Entropie moyenne
        ent = -np.mean(np.sum(A * np.log(A + 1e-12), axis=-1))
        entropies.append(ent)

        if i % 100 == 0 or i == n_iters - 1:
            print(f"  Iter {i:4d} │ Loss {loss:.5f} │ Entropy {ent:.4f}")

    # ---------- Analyse finale ----------
    print("\n--- Résultats ---")
    A_final, _ = model.forward(T)

    attn_row = A_final[idx_mask]
    q = np.linspace(-np.pi, np.pi, N)
    dist_q = np.abs(q - q[idx_mask])

    corr = np.corrcoef(attn_row, dist_q)[0, 1]
    print(f"Corrélation attention (q=0) ↔ distance spatiale : {corr:.4f}")

    if corr < -0.5:
        print("→ SUCCÈS : l'attention s'est concentrée sur le voisinage local.")
    else:
        print("→ L'attention a privilégié d'autres structures (fréquences, etc.).")

    # ---------- Visualisation ----------
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

    axes[0].plot(losses, color="crimson", lw=1.6)
    axes[0].set_title("Perte (MSE)")
    axes[0].grid(alpha=0.3)

    axes[1].plot(entropies, color="royalblue", lw=1.6)
    axes[1].axhline(np.log(N), color="k", ls="--", lw=1, label="max")
    axes[1].set_title("Entropie moyenne")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    im = axes[2].imshow(A_final, cmap="viridis", aspect="auto")
    axes[2].set_title("Matrice d'attention finale")
    fig.colorbar(im, ax=axes[2], fraction=0.046)

    plt.tight_layout()
    plt.savefig("resultats_attention_final.png", dpi=150)
    print("\n✅ Graphique sauvegardé : resultats_attention_final.png")

if __name__ == "__main__":
    main()
