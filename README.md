# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers des structures cristallines 3D de type CsCl (chlorure de césium).

## 🚀 Fonctionnalités

- **Moteurs physiques** : matrice dynamique pour une chaîne 1D (`v57`) et pour des réseaux cubiques 3D 6x6 (`v61_phonon_3d.py` et `v72_second_neighbors.py`).
- **Suivi de branches** : algorithmes d'affectation (Kuhn-Munkres) avec gestion des dégénérescences et alignement de phase (Procrustes).
- **Validation** : diagnostics et audits vérifiant la règle de somme acoustique et la stabilité 3D.
- **Masses atomiques** : masses réelles pour le CsCl ($M_{Cs} = 132.90$, $M_{Cl} = 35.45$ uma).

## 📂 Scripts principaux

- `v61_phonon_3d.py` : moteur 3D de base, 8 premiers voisins.
- `v72_second_neighbors.py` : moteur 3D enrichi des seconds voisins (Cs-Cs et Cl-Cl).
- `v73_rigorous_validation.py` : vitesses acoustiques (V72) par régression multi-q.
- `v74_full_bz_stability.py` : audit de stabilité sur une grille 3D de 3375 points (V72).
- `v75_direction_110_test.py` : test d'anisotropie de V61 le long de [110] (rapports attendus 3/7 : 1 : 11/7).
- `v70_phonon_acoustic_validation.py` : diagnostic acoustique du modèle V61.
- `audit_v61_rigorous.py` : audit de l'invariance par translation (V61).
- `v57_production_501.py` : moteur 1D historique.

Autres scripts présents dans le dépôt, non décrits ici : `v58_regression_test.py`, `v60_validate_artifacts.py`, `v61_cristal_3d.py`, `run_all.py`.

*(Les scripts V62 à V69 sont à la racine du dépôt à titre d'archives ; ils documentent l'évolution des algorithmes de suivi de branches.)*

## 📊 Résultats physiques et limites

1. **Modèle à premiers voisins (V61)**
   - **Dégénérescence exacte des vitesses acoustiques le long de [100]** : $v_{LA} = v_{TA} \approx 2380.5$ m/s. Elle vaut pour toutes les valeurs de $K_\parallel$ et $K_\perp$.
   - **Anisotropie** : le cristal n'est pas isotrope. Le cisaillement $C'$ est proportionnel à $K_\perp$ ; le rapport de Zener vaut $(K_\parallel + 2K_\perp)/(3K_\perp) = 7/3$.
     - [110] : trois branches distinctes (≈ 1558, 2381, 2984 m/s).
     - [111] : TA doublement dégénérée (≈ 1873 m/s), LA distincte (≈ 3160 m/s).
   - **Stabilité** : $C' > 0$, le cristal n'est pas marginalement instable.

2. **Modèle à seconds voisins (V72-V74)**
   - $K_{2\parallel} = 4.0$ N/m augmente $C_{11}$ ; $K_{2\perp} = 1.0$ N/m augmente $C_{44}$ et diminue $C_{12}$ d'autant ($\delta C_{12} = -\delta C_{44}$). La dégénérescence LA/TA le long de [100] est levée.
   - Vitesses le long de [100] :
     - TA : $2623.2 \pm 0.2$ m/s
     - LA : $3244.1 \pm 0.2$ m/s
   - L'onde longitudinale est ~23,7 % plus rapide que la transverse le long de [100]. Le rapport de Zener descend à environ 1,26.

**Limites physiques**
- **Ionicité** : le CsCl étant ionique, un modèle de ressorts ne peut pas capturer la séparation LO-TO (interaction de Coulomb à longue portée). Les vitesses sont paramétriques, pas des prédictions expérimentales.
- **Constantes de force** : les constantes des seconds voisins sont illustratives ; $K_{2\perp}$ produit un $C_{12}$ négatif.
- **Stabilité V72** : grille de 3375 points de la zone de Brillouin ; aucune valeur propre négative au-delà du bruit numérique (≈ −10⁻¹⁶ relatif).

## 🛠️ Installation & lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
pip install -r requirements.txt

# Validation du modèle enrichi (V73)
python v73_rigorous_validation.py

# Audit de stabilité (V74)
python v74_full_bz_stability.py

# Test d'anisotropie de V61 (V75) — attendu : [PASS]
python v75_direction_110_test.py
```

## 🔭 Perspectives

Pour transformer ce modèle paramétrique en outil quantitatif pour le CsCl réel, la voie standard est la **DFT** ; la précision dépend de la fonctionnelle et des pseudopotentiels.

1. **Constantes réelles** : un code DFT (Quantum ESPRESSO ou autre) couplé à `Phonopy` pour calculer les forces et ajuster $K_\parallel$ et $K_\perp$.
2. **Correction LO-TO** : intégrer les charges de Born effectives et le tenseur diélectrique. Les modèles à portée finie, comme MACE, ont des limites connues sur ces interactions coulombiennes.
3. **Validation croisée** : comparer aux courbes mesurées par diffusion inélastique de neutrons : A. A. Z. Ahmad, H. G. Smith, N. Wakabayashi, M. K. Wilkinson, *Phys. Rev. B* **6**, 3956 (1972), à 78 K.

*Modèles existants* : la littérature donne des paramètres de *shell model* (par ex. Jindal et Mahesh, 1975), dont la conversion vers ce modèle de ressorts n'est pas directe. Les constantes élastiques $C_{11}$, $C_{12}$, $C_{44}$ fournissent 3 contraintes pour 4 paramètres, et sans interaction coulombienne les branches optiques ne s'ajusteront pas.

---

Développé sur Termux et VS Code.
