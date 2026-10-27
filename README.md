# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers des structures cristallines 3D de type CsCl (Chlorure de Césium).

## 🚀 Fonctionnalités

- **Moteurs physiques** : Calcul de la matrice dynamique pour une chaîne 1D (`v57`) et des réseaux cubiques 3D 6x6 (`v61_phonon_3d.py` et `v72_second_neighbors.py`).
- **Suivi de branches** : Algorithmes d'affectation (Kuhn-Munkres) avec gestion des dégénérescences et alignement de phase (Procrustes).
- **Validation** : Pipelines de diagnostic et audits formels vérifiant la règle de somme acoustique et la stabilité 3D.
- **Masses atomiques** : Prise en compte des masses réelles pour le CsCl ($M_{Cs} = 132.90$, $M_{Cl} = 35.45$ uma).

## 📂 Scripts Principaux

- `v61_phonon_3d.py` : Moteur 3D de base avec les 8 premiers voisins.
- `v72_second_neighbors.py` : Moteur 3D enrichi avec les seconds voisins (Cs-Cs et Cl-Cl).
- `v73_rigorous_validation.py` : Calcul des vitesses acoustiques (V72) par régression multi-q.
- `v74_full_bz_stability.py` : Audit de stabilité sur une grille 3D de 3375 points pour le modèle V72.
- `v75_direction_110_test.py` : Test d'anisotropie de V61 le long de [110].
- `v70_phonon_acoustic_validation.py` : Pipeline de diagnostic pour le modèle V61.
- `audit_v61_rigorous.py` : Audit strict de l'invariance par translation (V61).
- `v57_production_501.py` : Moteur 1D historique.
- `v58_regression_test.py`, `v60_validate_artifacts.py`, `v61_cristal_3d.py`, `run_all.py` : Scripts de test et d'orchestration associés au moteur 1D et à la visualisation.

*(Les scripts V62 à V69 sont conservés à la racine du dépôt à titre d'archives, documentant l'évolution des algorithmes de suivi de branches).*

## 📊 Résultats Physiques et Limites

L'analyse numérique met en évidence deux comportements distincts selon le modèle :

1. **Modèle à premiers voisins (V61)** : 
   - L'analyse révèle une **dégénérescence exacte des vitesses acoustiques le long de la direction [100]** ($v_{LA} = v_{TA} \approx 2380.5$ m/s).
   - **Anisotropie** : Le cristal n'est pas isotrope. Le cisaillement $C'$ est proportionnel à $K_\perp$ ; le rapport de Zener vaut $(K_\parallel + 2K_\perp)/(3K_\perp) = 7/3$. Les directions [110] et [111] présentent des branches distinctes.
   - **Stabilité** : Le cristal n'est pas marginalement instable (le module de cisaillement $C' > 0$).
   
2. **Modèle à seconds voisins (V72-V74)** : 
   - L'ajout des interactions Cs-Cs et Cl-Cl brise la dégénérescence le long de [100]. L'introduction d'une raideur transverse pour les seconds voisins ($K_{2\parallel} = 4.0$ N/m et $K_{2\perp} = 1.0$ N/m) modifie le module $C_{44}$.
   - Les vitesses calculées le long de [100] deviennent différentes :
     - Vitesse Transverse (TA) : $2623.15 \pm 0.2$ m/s
     - Vitesse Longitudinale (LA) : $3243.98 \pm 0.2$ m/s
   - **L'onde longitudinale est ~23.7% plus rapide que l'onde transverse le long de cette direction.** Le rapport de Zener descend à environ 1.26.

**Limites physiques** :
- **Ionicité** : Le CsCl étant un cristal ionique, un modèle de ressorts phénoménologique ne peut pas capturer la séparation LO-TO (interaction de Coulomb à longue portée). Les vitesses calculées sont paramétriques et ne constituent pas des prédictions expérimentales.
- **Constantes de force** : Les constantes utilisées pour les seconds voisins sont illustratives et produisent un module $C_{12}$ négatif, un artefact du modèle paramétrique.
- **Stabilité V72** : Testée sur 3375 points de la zone de Brillouin 3D. Aucune instabilité physique significative n'a été détectée.

## 🛠️ Installation & Lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
# Installer les dépendances
pip install -r requirements.txt

# Lancer la validation du modèle enrichi (V73)
python v73_rigorous_validation.py

# Lancer l'audit de stabilité (V74)
python v74_full_bz_stability.py

# Lancer le test d'anisotropie (V75)
python v75_direction_110_test.py
