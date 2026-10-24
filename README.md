# Phonon Engine — 1D Diatomic Chain Dispersion

Ce projet est un moteur de simulation physique écrit en Python qui calcule la relation de dispersion et les vitesses de groupe des phonons dans une chaîne diatomique unidimensionnelle. 

Il est né d'une exploration théorique sur les tachyons et les ordinateurs quantiques, et a évolué vers un outil de physique du solide rigoureusement testé (60 itérations de développement et de validation).

## 🚀 Fonctionnalités

- **Modèle physique** : Chaîne diatomique 1D avec masses m1, m2 et ressort de raideur k.
- **Calculs** : Formules analytiques rationalisées pour éviter l'annulation catastrophique près de q=0.
- **Robustesse** : Testé avec 14 assertions de régression (symétrie, monotonie, mode de Goldstone, cohérence des vitesses de groupe).
- **Exports** : Génération des données en `.csv` et des graphiques en `.png`.

## 📊 Paramètres Physiques

- Masses : m1 = 2.0 uma, m2 = 4.0 uma
- Constante de raideur : k = 10.0 N/m
- Paramètre de maille : a = 5.0 × 10⁻¹⁰ m
- Nombre de points : 501 (inclut exactement q=0)

*Note : Ces paramètres sont illustratifs et ne constituent pas un ajustement à un matériau expérimental particulier.*

## 🛠️ Installation & Lancement

Nécessite Python 3.9+, NumPy et Matplotlib.

```bash
# Cloner le dépôt
git clone https://github.com/berryberryship-debug/phonon-engine.git
cd phonon-engine

# Lancer la génération des données et la validation complète
python run_all.py
