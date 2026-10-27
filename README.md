# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers des structures cristallines 3D de type CsCl (chlorure de césium).

![Architecture du projet](figures/01_architecture.png)

## 🚀 Fonctionnalités

- **Moteurs physiques** : matrice dynamique pour une chaîne 1D (`v57`) et pour des réseaux cubiques 3D 6x6 (`v61_phonon_3d.py` et `v72_second_neighbors.py`).
- **Suivi de branches** : algorithmes d'affectation (Kuhn-Munkres) pour le suivi des modes.
- **Validation** : diagnostics et audits vérifiant la règle de somme acoustique et la stabilité 3D.
- **Masses du modèle** : $M_{\mathrm{Cs}}=132.90545196$ uma et $M_{\mathrm{Cl}}=35.45$ uma, converties en kg dans le moteur via $U=1.66053906660\times10^{-27}$ kg/uma.

## 📂 Scripts principaux

- `v61_phonon_3d.py` : moteur 3D de base, 8 premiers voisins.
- `v72_second_neighbors.py` : moteur 3D enrichi des seconds voisins (Cs-Cs et Cl-Cl).
- `v73_rigorous_validation.py` : vitesses acoustiques (V72) par régression multi-q.
- `v74_full_bz_stability.py` : audit de stabilité sur une grille 3D de 3375 points (V72).
- `v75_direction_110_test.py` : test d'hypothèse sur l'anisotropie de V61 le long de [110].
- `v70_phonon_acoustic_validation.py` : diagnostic acoustique du modèle V61.
- `audit_v61_rigorous.py` : audit de l'invariance par translation (V61).
- `v57_production_501.py` : moteur 1D historique.

Autres scripts présents dans le dépôt, non décrits ici : `v58_regression_test.py`, `v60_validate_artifacts.py`, `v61_cristal_3d.py`, `run_all.py`, `generate_readme_figures.py`.

## 🛡️ Carte des validations

![Carte des validations](figures/02_validation_map.png)

## 📊 Résultats du modèle et limites

![Résultats et limites](figures/03_results_and_limits.png)

### Modèle à premiers voisins (V61)

Le long de [100], le modèle prédit une dégénérescence des vitesses acoustiques longitudinales et transverses :
- $v_{LA} = v_{TA} \approx 2380.5$ m/s.

Cette dégénérescence est une conséquence de la forme choisie pour les ressorts premiers voisins et de la symétrie cubique du modèle. Elle ne constitue pas une prédiction expérimentale pour le CsCl réel.

Le rapport de Zener calculé pour ce jeu de paramètres (à partir des constantes élastiques effectives) est :
$$A_Z = \frac{2C_{44}}{C_{11}-C_{12}} \approx 7/3.$$

### Modèle à seconds voisins (V72–V74)

Le long de [100], la régression linéaire dans le régime acoustique (sur 20 points, $q \in [10^6, 10^8]$ m$^{-1}$) donne :
- $v_{TA} = 2623.2 \pm 0.2$ m/s ;
- $v_{LA} = 3244.1 \pm 0.2$ m/s.

Les incertitudes sont les erreurs standard de régression sur les points de simulation. Elles n'incluent ni la courbure systématique de $\omega(q)$, ni une incertitude expérimentale.

L'onde longitudinale est environ 23,7 % plus rapide que l'onde transverse. Le rapport de Zener du modèle vaut environ 1,26.

### Limites

- Absence d'interactions coulombiennes à longue portée : le modèle ne produit pas de séparation LO-TO.
- Les constantes de force des seconds voisins sont illustratives. Le jeu de paramètres actuel conduit à un $C_{12}$ négatif.
- La stabilité a été vérifiée sur une grille finie uniforme de $15^3=3375$ points (incluant $\Gamma$, X, M, R). Aucune valeur propre inférieure à $-10^{-5}\lambda_{\max}$ n'a été détectée. Cela ne constitue pas une preuve de stabilité sur tout le continuum de la zone de Brillouin.
- Aucune comparaison quantitative avec Ahmad et al. (1972) n'est revendiquée tant que les points expérimentaux ne sont pas numérisés avec leurs unités et incertitudes.

## 🛠️ Installation & lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
pip install -r requirements.txt

# Lancer les validations
python v73_rigorous_validation.py
python v74_full_bz_stability.py
python v75_direction_110_test.py
