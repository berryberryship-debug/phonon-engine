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

1. **Règle de somme acoustique** : La vérification directe de l'action de la matrice dynamique sur les trois vecteurs de translation rigide à \(\Gamma\) donne un résidu relatif de l'ordre de \(10^{-17}\), inférieur à la précision machine. L'invariance par translation est respectée par construction numérique.
2. **Stabilité sur le trajet \(\Gamma \to X\)** : Aucune instabilité n'est détectée. L'analyse utilise une tolérance relative à l'échelle spectrale de la matrice pour distinguer les erreurs d'arrondi d'une véritable instabilité.
3. **Égalité des vitesses acoustiques** : L'analyse numérique du modèle à 8 premiers voisins avec des ressorts centraux révèle une égalité des vitesses longitudinales et transverses (\(v_{LA} = v_{TA} \approx 2379.6\) m/s). Ce résultat, cohérent avec la structure de la matrice dynamique implémentée, n'est pas un théorème universel pour tous les cristaux, mais une propriété spécifique de ce modèle phénoménologique simplifié.

## 🛠️ Installation & Lancement

Nécessite Python 3.9+, NumPy, SciPy et Matplotlib.

```bash
# Installer les dépendances
pip install -r requirements.txt

# Lancer le diagnostic des vitesses acoustiques (V70)
python v70_phonon_acoustic_validation.py

# Lancer l'audit de la règle de somme (V71)
python audit_v61_rigorous.py
