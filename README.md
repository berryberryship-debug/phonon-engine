# Phonon Engine — 1D Diatomic Chain Dispersion

Ce projet est un moteur de simulation physique écrit en Python qui calcule la relation de dispersion et les vitesses de groupe des phonons dans une chaîne **diatomique** unidimensionnelle (deux atomes de masses différentes par maille élémentaire). 

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

## 🔬 Nature exacte du modèle (Diatomique)

Contrairement à une chaîne monoatomique qui ne possède qu'une seule branche de dispersion lisse, ce modèle simule une chaîne **diatomique** (deux atomes par maille). Cela se traduit par deux branches distinctes :
1. La branche **acoustique** (où les deux atomes vibrent en phase).
2. La branche **optique** (où les deux atomes vibrent en opposition de phase).

C'est cette nature diatomique qui explique la "cuspide" (le saut brusque) de la vitesse de groupe de la branche acoustique à q=0. La fréquence se comporte comme ω ∝ |q|, dont la dérivée (la vitesse) est discontinue en 0.

## 🛠️ Installation & Reproduction

Nécessite Python 3.9+ (testé avec NumPy et Matplotlib récents).

Pour cloner et lancer le pipeline complet (Génération V57 + Validation des artefacts V60 + Tests de régression V58) :

    git clone https://github.com/berryberryship-debug/phonon-engine.git
    cd phonon-engine
    python run_all.py

Pour reproduire spécifiquement les tests de validation (V58) et de robustesse (V59) :

    # Lancer uniquement les 14 tests de régression physique (V58)
    python v58_regression_test.py

    # Lancer le test de robustesse du pipeline (V59)
    python v59_test_pipeline.py

## 📈 Résultats

- `phonon_dispersion_501.csv` : Tableau de 501 lignes avec q, ω_aco, ω_opt, v_g_aco, v_g_opt.
- `phonon_curves_501.png` : Graphique des branches de dispersion et des vitesses de groupe.

## 📐 Définition de l'erreur relative (Vitesses de groupe)

Les 14 tests de régression comparent les vitesses de groupe analytiques exactes (v_g = dω/dq) aux vitesses calculées par différences finies sur la grille discrète. L'erreur relative maximale est calculée comme suit :
Erreur = max( |v_g_numérique - v_g_analytique| ) / max( |v_g_analytique| )
Cette normalisation par l'échelle globale (le maximum) évite les divisions par zéro près de q=0.

## ⚠️ Portée et limites

Les tests vérifient la cohérence numérique interne de l'implémentation pour les propriétés contrôlées du modèle idéal. Ils ne constituent pas une validation expérimentale ni une preuve que les paramètres représentent un matériau réel.

---
*Développé et validé sur Termux (Android) via un processus itératif de développement et de relecture critique.*
