# Phonon Engine — 1D Diatomic Chain Dispersion

Ce projet est un moteur de simulation physique écrit en Python qui calcule la relation de dispersion et les vitesses de groupe des phonons dans une chaîne diatomique unidimensionnelle. 

Il est né d'une exploration théorique sur les tachyons et les ordinateurs quantiques, et a évolué vers un outil de physique du solide rigoureusement testé.

## 🚀 Fonctionnalités

- **Modèle physique** : Chaîne diatomique 1D avec masses m1, m2 et ressort de raideur k.
- **Calculs** : Formules analytiques rationalisées pour éviter l'annulation catastrophique près de q=0.
- **Exports** : Génération des données en `.csv` et des graphiques en `.png`.
- **Pipeline de validation** : Intègre la génération, la validation des artefacts et les tests de régression.

## 📊 Paramètres Physiques

- Masses : m1 = 2.0 uma, m2 = 4.0 uma
- Constante de raideur : k = 10.0 N/m
- Paramètre de maille : a = 5.0 × 10⁻¹⁰ m
- Nombre de points : 501 (inclut exactement q=0)

*Note : Ces paramètres sont illustratifs et ne constituent pas un ajustement à un matériau expérimental particulier.*

## 🛠️ Installation & Lancement

Nécessite Python 3.9+ (testé avec NumPy et Matplotlib récents).

    # Cloner le dépôt
    git clone https://github.com/berryberryship-debug/phonon-engine.git
    cd phonon-engine

    # Lancer le pipeline complet (Génération + Validations)
    python run_all.py

Le pipeline exécute séquentiellement :
1. `v57_production_501.py` : Calcul, export CSV et génération du graphique PNG.
2. `v60_validate_artifacts.py` : Vérification structurelle des fichiers générés (taille, format, point critique).
3. `v58_regression_test.py` : 14 tests de régression validant la cohérence numérique (symétrie, monotonie, mode de Goldstone, vitesses de groupe).

## 📈 Résultats

- `phonon_dispersion_501.csv` : Tableau de 501 lignes avec q, ω_aco, ω_opt, v_g_aco, v_g_opt.
- `phonon_curves_501.png` : Graphique des branches de dispersion et des vitesses de groupe.

## ⚠️ Portée et limites

Les tests vérifient la cohérence numérique interne de l'implémentation pour les propriétés contrôlées du modèle idéal. Ils ne constituent pas une validation expérimentale ni une preuve que les paramètres représentent un matériau réel.

---
*Développé et validé sur Termux (Android) via un processus itératif de développement et de relecture critique.*
