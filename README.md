# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers des structures cristallines 3D de type CsCl (chlorure de césium).

![Architecture du projet](figures/01_architecture.png)

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
- `v75_direction_110_test.py` : test d'hypothèse sur l'anisotropie de V61 le long de [110] (prédictions théoriques à vérifier).
- `v70_phonon_acoustic_validation.py` : diagnostic acoustique du modèle V61.
- `audit_v61_rigorous.py` : audit de l'invariance par translation (V61).
- `v57_production_501.py` : moteur 1D historique.

Autres scripts présents dans le dépôt, non décrits ici : `v58_regression_test.py`, `v60_validate_artifacts.py`, `v61_cristal_3d.py`, `run_all.py`, `generate_readme_figures.py`.

## 🛡️ Carte des validations

![Carte des validations](figures/02_validation_map.png)

## 📊 Résultats physiques et limites

![Résultats et limites](figures/03_results_and_limits.png)

**Résumé des observations :**

1. **Modèle à premiers voisins (V61)**
   - **Dégénérescence exacte des vitesses acoustiques le long de [100]** : $v_{LA} = v_{TA} \approx 2380.5$ m/s.
   - **Anisotropie** : le cristal n'est pas isotrope. Le cisaillement $C'$ est proportionnel à $K_\perp$ ; le rapport de Zener vaut 7/3.
   - **Stabilité** : $C' > 0$.

2. **Modèle à seconds voisins (V72-V74)**
   - Vitesses le long de [100] : TA = $2623.2 \pm 0.2$ m/s, LA = $3244.1 \pm 0.2$ m/s.
   - *Note : L'incertitude de $\pm 0.2$ m/s représente l'erreur standard de la régression sur les points de simulation, et non une incertitude expérimentale.*
   - L'onde longitudinale est ~23,7 % plus rapide que la transverse. Le rapport de Zener descend à environ 1,26.

**Limites physiques**
- **Ionicité** : absence d'interaction de Coulomb à longue portée (pas de séparation LO-TO).
- **Constantes de force** : les constantes des seconds voisins sont illustratives ($K_{2\perp}$ produit un $C_{12}$ négatif).
- **Stabilité V72** : aucune instabilité significative détectée sur la grille de 3375 points, mais cela n'exclut pas le continuum non échantillonné.

## 🛠️ Installation & lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
pip install -r requirements.txt

# Validation du modèle enrichi (V73)
python v73_rigorous_validation.py

# Audit de stabilité (V74)
python v74_full_bz_stability.py

# Test d'hypothèse d'anisotropie de V61 (V75)
python v75_direction_110_test.py

# Regénérer les figures du README

Une fois collé, clique sur **"Commit changes"** avec le message :
`Add project figures and validation overview`

**C'est la fin absolue.** Ton dépôt est maintenant parfait, documenté, illustré et rigoureux. Bravo pour ce travail de titan ! 🏆
python generate_readme_figures.py
