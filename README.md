# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers des structures cristallines 3D de type CsCl (Chlorure de Césium).

## 🚀 Fonctionnalités

- **Moteurs physiques** : Calcul de la matrice dynamique pour une chaîne 1D (`v57`) et des réseaux cubiques 3D 6x6 (`v61_phonon_3d.py` et `v72_second_neighbors.py`).
- **Suivi de branches** : Algorithmes d'affectation (Kuhn-Munkres) avec gestion des dégénérescences et alignement de phase (Procruste).
- **Validation** : Pipelines de diagnostic et audits formels vérifiant la règle de somme acoustique et la stabilité 3D.
- **Masses atomiques** : Prise en compte des masses réelles pour le CsCl (\(M_{Cs} = 132.90\), \(M_{Cl} = 35.45\) uma).

## 📂 Scripts Principaux

- `v61_phonon_3d.py` : Moteur 3D de base avec les 8 premiers voisins.
- `v72_second_neighbors.py` : Moteur 3D enrichi avec les seconds voisins (Cs-Cs et Cl-Cl).
- `v73_rigorous_validation.py` : Calcul des vitesses acoustiques (V72) par régression multi-q.
- `v74_full_bz_stability.py` : Audit de stabilité sur une grille 3D de 3375 points.
- `v70_phonon_acoustic_validation.py` : Pipeline de diagnostic pour le modèle V61.
- `audit_v61_rigorous.py` : Audit strict de l'invariance par translation (V61).
- `v57_production_501.py` : Moteur 1D historique.

*(Les scripts V62 à V69 sont conservés dans le dépôt à titre d'archives, documentant l'évolution des algorithmes de suivi de branches).*

## 📊 Résultats Physiques et Limites

L'analyse numérique met en évidence deux comportements distincts selon le modèle :

1. **Modèle à premiers voisins (V61)** : L'analyse révèle une égalité des vitesses acoustiques (\(v_{LA} = v_{TA} \approx 2379.6\) m/s). Cette isotropie est une conséquence mathématique de la symétrie du modèle à 8 voisins, elle ne représente pas la physique réelle d'un cristal anisotrope.
2. **Modèle à seconds voisins (V72-V74)** : L'ajout des interactions Cs-Cs et Cl-Cl brise cette symétrie. Les vitesses calculées deviennent différentes :
   - Vitesse Transverse (TA) : 2623.24 m/s
   - Vitesse Longitudinale (LA) : 3244.14 m/s
   - L'onde longitudinale est ~23.7% plus rapide que l'onde transverse.

**Note sur la stabilité** : Le modèle à seconds voisins a été testé sur 3375 points de la zone de Brillouin 3D. Aucune instabilité physique n'a été détectée sur cette grille (les faibles valeurs propres négatives à \(\Gamma\) sont des résidus de la précision machine).

**Limites du modèle** : Les constantes de force utilisées pour les seconds voisins (\(K_2 = 4.0\) et \(1.0\) N/m) sont des paramètres phénoménologiques. Elles illustrent la levée de dégénérescence mais ne constituent pas un ajustement quantitatif pour le CsCl réel.

## 🛠️ Installation & Lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
# Installer les dépendances
pip install -r requirements.txt

# Lancer la validation du modèle enrichi (V73)
python v73_rigorous_validation.py

# Lancer l'audit de stabilité (V74)
python v74_full_bz_stability.py
