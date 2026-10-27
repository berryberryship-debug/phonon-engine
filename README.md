# Phonon Engine — 1D & 3D Phonon Dispersion

Ce projet est un moteur de simulation en Python qui calcule la relation de dispersion des phonons. Il a débuté avec une chaîne diatomique 1D (modèle de ressorts) et a évolué vers une structure cristalline 3D de type CsCl (Chlorure de Césium).

## 🚀 Fonctionnalités

- **Moteurs physiques** : Calcul de la matrice dynamique pour une chaîne 1D (`v57`) et un réseau cubique 3D 6x6 (`v61_phonon_3d.py`).
- **Suivi de branches** : Algorithme d'affectation (Kuhn-Munkres) avec gestion des dégénérescences et alignement de phase (Procruste). Ces algorithmes permettent de suivre les modes de façon continue, sans garantir à eux seuls leur exactitude physique.
- **Validation** : Pipeline de diagnostic (V70) et audits formels (V71) vérifiant la règle de somme acoustique.
- **Masses atomiques** : Prise en compte des masses réelles pour le CsCl (\(M_{Cs} = 132.90\), \(M_{Cl} = 35.45\) uma).

## 📂 Scripts Principaux

- `v61_phonon_3d.py` : Le moteur de calcul de la matrice dynamique 3D pour le CsCl.
- `v70_phonon_acoustic_validation.py` : Pipeline de diagnostic extrayant les branches acoustiques et calculant les vitesses du son par régression.
- `audit_v61_rigorous.py` : Audit strict de l'invariance par translation.
- `v57_production_501.py` : Moteur 1D historique.

*(Les scripts V62 à V69 sont conservés dans le dépôt à titre d'archives, documentant l'évolution des algorithmes de suivi de branches).*

## 📊 Résultats Physiques et Limites

Les audits menés sur le moteur `v61` ont établi les points suivants :

1. **Règle de somme acoustique** : La vérification directe de l'action de la matrice dynamique sur les trois vecteurs de translation rigide à \(\Gamma\) donne un résidu relatif de l'ordre de \(10^{-17}\), compatible avec la précision machine (double précision). L'invariance par translation est respectée par construction numérique.
2. **Stabilité sur le trajet \(\Gamma \to X\)** : Aucune instabilité n'est détectée. L'analyse utilise une tolérance relative à l'échelle spectrale de la matrice pour distinguer les erreurs d'arrondi d'une véritable instabilité.
3. **Égalité des vitesses acoustiques** : L'analyse numérique du modèle à 8 premiers voisins avec des ressorts centraux révèle une égalité des vitesses longitudinales et transverses (\(v_{LA} = v_{TA} \approx 2379.6\) m/s). Ce résultat est un constat numérique issu de ce modèle spécifique. Il reste à démontrer analytiquement si cette égalité est une conséquence stricte de la géométrie à 8 voisins ou du jeu de paramètres choisi.

## 🛠️ Installation & Lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
# Installer les dépendances
pip install -r requirements.txt

# Lancer le diagnostic des vitesses acoustiques (V70)
python v70_phonon_acoustic_validation.py

# Lancer l'audit de la règle de somme (V71)
python audit_v61_rigorous.py
## 🧬 Modèle enrichi avec seconds voisins (V72-V74)

Pour dépasser la limite d'isotropie du modèle à 8 premiers voisins, les interactions de seconds voisins (Cs-Cs et Cl-Cl) ont été ajoutées. 

- **Stabilité** : Le modèle est stable sur toute la zone de Brillouin 3D (testé sur 3375 points).
- **Vitesses acoustiques** : La symétrie est brisée. La régression sur plusieurs vecteurs d'onde donne :
  - Vitesse Transverse (TA) : 2623.24 m/s
  - Vitesse Longitudinale (LA) : 3244.14 m/s
  - L'onde longitudinale est ~23.7% plus rapide que l'onde transverse.

*Note : Les constantes de force pour les seconds voisins (K2_parallel=4.0 N/m, K2_perp=1.0 N/m) sont des paramètres phénoménologiques. Elles illustrent la levée de dégénérescence mais ne constituent pas un ajustement quantitatif pour le CsCl réel.*
