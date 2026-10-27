# Phonon Engine — Dispersion phononique 1D et 3D

Ce projet est un moteur de simulation en Python pour calculer des relations de dispersion phononiques. Il est parti d’une chaîne diatomique 1D à ressorts, puis a été étendu à des réseaux cubiques 3D de type CsCl (chlorure de césium), décrits par une maille primitive à deux atomes.

![Architecture du projet](figures/01_architecture.png)

## 🚀 Fonctionnalités

- **Moteurs physiques** : matrice dynamique pour une chaîne diatomique 1D (`v57`) et pour des réseaux cubiques 3D à deux atomes par maille primitive (`v61_phonon_3d.py`, `v72_second_neighbors.py`).
- **Suivi de branches** : affectation optimale des modes par algorithme Kuhn-Munkres (Hongrois). Le traitement des dégénérescences dépend de l’implémentation du suivi modal utilisée ; il doit être vérifié dans le code et les tests avant toute affirmation de continuité branche par branche.
- **Validation numérique** : hermiticité, règle de somme acoustique, vecteurs de translation rigide, stabilité des valeurs propres et vitesses acoustiques.
- **Masses du modèle** : \(M_{\mathrm{Cs}} = 132.90545196\) uma et \(M_{\mathrm{Cl}} = 35.45\) uma, converties en kg dans le moteur via \(U = 1.66053906660 \times 10^{-27}\) kg/uma.

## 📂 Scripts principaux

| Script | Rôle |
|---|---|
| `v61_phonon_3d.py` | Moteur 3D de base, huit premiers voisins Cs–Cl |
| `v72_second_neighbors.py` | Moteur 3D enrichi : seconds voisins Cs–Cs et Cl–Cl |
| `v73_rigorous_validation.py` | Vitesses acoustiques de V72 par régression multi-\(q\) |
| `v74_full_bz_stability.py` | Audit de stabilité de V72 sur une grille 3D de la zone de Brillouin |
| `v75_direction_110_test.py` | Test d’anisotropie de V61 le long de \([110]\) |
| `v70_phonon_acoustic_validation.py` | Diagnostic acoustique de V61 |
| `audit_v61_rigorous.py` | Audit de l’invariance par translation de V61 |
| `v57_production_501.py` | Moteur 1D historique |

Autres scripts présents dans le dépôt, non décrits ici : `v58_regression_test.py`, `v60_validate_artifacts.py`, `v61_cristal_3d.py`, `run_all.py` et `generate_readme_figures.py`.

## 🛡️ Carte des validations

![Carte des validations](figures/02_validation_map.png)

## 📊 Résultats du modèle et limites

![Résultats et limites](figures/03_results_and_limits.png)

### Modèle à premiers voisins — V61

Le long de \([100]\), le modèle prédit une dégénérescence des vitesses acoustiques longitudinales et transverses :

\[
v_{\mathrm{LA}} = v_{\mathrm{TA}} \approx 2380.5\ \mathrm{m/s}.
\]

Cette dégénérescence provient de la forme choisie pour les ressorts premiers voisins et de la symétrie cubique du modèle. Elle ne constitue pas une prédiction expérimentale pour le CsCl réel.

Le rapport de Zener, calculé à partir des constantes élastiques effectives du modèle, vaut :

\[
A_Z = \frac{2C_{44}}{C_{11}-C_{12}} \approx \frac{7}{3}.
\]

### Modèle à seconds voisins — V72–V74

Le long de \([100]\), une régression linéaire dans le régime acoustique, sur 20 points avec \(q \in [10^6, 10^8]\ \mathrm{m^{-1}}\), donne :

\[
v_{\mathrm{TA}} = 2623.2 \pm 0.2\ \mathrm{m/s},
\]

\[
v_{\mathrm{LA}} = 3244.1 \pm 0.2\ \mathrm{m/s}.
\]

Ces valeurs et leurs erreurs standard doivent être reproduites par `v73_rigorous_validation.py`. Les incertitudes caractérisent la régression sur les points de simulation ; elles n’incluent ni la courbure systématique de \(\omega(q)\), ni une incertitude expérimentale.

Dans ce modèle, l’onde longitudinale est environ 23,7 % plus rapide que l’onde transverse. Le rapport de Zener vaut environ 1,26.

### Stabilité — V72–V74

La stabilité a été vérifiée sur une grille uniforme de \(15^3 = 3375\) points, incluant \(\Gamma\), X, M et R. Aucune valeur propre inférieure à \(-10^{-5}\lambda_{\max}\) n’a été détectée.

Cette affirmation suppose que :

- les valeurs propres analysées sont bien celles de la matrice dynamique massiquement pondérée ;
- le seuil \(-10^{-5}\lambda_{\max}\) est appliqué de façon homogène sur toute la grille ;
- la grille et le critère sont ceux implémentés dans `v74_full_bz_stability.py`.

Une grille finie ne prouve pas la stabilité sur tout le continuum de la zone de Brillouin.

### Limites physiques

- Le modèle n’inclut pas d’interaction coulombienne à longue portée ; il ne produit donc pas de séparation LO–TO.
- Les constantes de force des seconds voisins sont illustratives. Avec le jeu de paramètres actuel, le modèle conduit à un \(C_{12}\) négatif.
- Aucune comparaison quantitative avec Ahmad et al. (1972) n’est revendiquée tant que les points expérimentaux ne sont pas numérisés avec leurs unités, directions cristallographiques et incertitudes.

## 🛠️ Installation et lancement

Prérequis : Python 3.9 ou plus récent, NumPy, SciPy et Matplotlib.

```bash
pip install -r requirements.txt
```

### Audits structurels et acoustiques

```bash
python audit_v61_rigorous.py
python v70_phonon_acoustic_validation.py
```

### Vitesses acoustiques et stabilité — V72

```bash
python v73_rigorous_validation.py
python v74_full_bz_stability.py
```

### Anisotropie — V61

```bash
python v75_direction_110_test.py
```

### Régénération des figures

```bash
python generate_readme_figures.py
```

## 📌 Statut scientifique

Ce dépôt valide la cohérence interne et la stabilité numérique de modèles phononiques simplifiés. Il ne prétend pas reproduire quantitativement les dispersions expérimentales du CsCl sans données expérimentales numérisées, leurs unités, leurs directions cristallographiques et leurs incertitudes.

Les résultats chiffrés du README doivent être considérés comme validés uniquement s’ils sont régénérés par les scripts correspondants et conservés dans des artifacts vérifiables.
