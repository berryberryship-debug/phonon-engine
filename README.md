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

Cette affirmation suppose que les valeurs propres analysées sont bien celles de la matrice dynamique massivement pondérée, que le seuil est appliqué de façon homogène sur toute la grille, et que la grille et le critère sont ceux implémentés dans `v74_full_bz_stability.py`.

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

### 1. Audits structurels et acoustiques — V61

```bash
python audit_v61_rigorous.py
python v70_phonon_acoustic_validation.py
```

### 2. Vitesses acoustiques et stabilité — V72

```bash
python v73_rigorous_validation.py
python v74_full_bz_stability.py
```

### 3. Test d’anisotropie — V61

```bash
python v75_direction_110_test.py
```

### 4. Régénération des figures

```bash
python generate_readme_figures.py
```

## ❓ FAQ scientifique

### Fondements physiques

#### Q1. Pourquoi la matrice dynamique est-elle de taille 6 × 6 ?

**[Démontré]**

Le modèle CsCl utilise une maille primitive à deux atomes, Cs et Cl. Chaque atome possède trois degrés de liberté cartésiens, donc la matrice dynamique contient \(3 \times 2 = 6\) modes : trois branches acoustiques et trois branches optiques.

Cette dimension est correcte uniquement si les deux atomes sont traités comme des sous-réseaux indépendants dans la maille primitive.

#### Q2. Qu’est-ce que la règle de somme acoustique ?

**[Démontré]**

Lorsque tous les atomes sont translatés rigoureusement de la même quantité, aucun ressort n’est étiré. À \(\Gamma\), la matrice dynamique doit donc posséder trois vecteurs propres de translation rigide associés à des fréquences nulles.

Dans ce projet, le test vérifie que les trois vecteurs \((\sqrt{M_{\mathrm{Cs}}},0,0,\sqrt{M_{\mathrm{Cl}}},0,0)\), et ses analogues selon \(y\) et \(z\), appartiennent au noyau de \(D(\Gamma)\).

#### Q3. Pourquoi le modèle ne montre-t-il pas de séparation LO–TO ?

**[Démontré pour ce modèle]**

Le modèle actuel utilise uniquement des interactions de type ressort à courte portée. Il n’inclut pas les forces coulombiennes à longue portée ni les charges ioniques effectives.

Sans ces termes, la fréquence des modes optiques longitudinaux et transversaux peut rester dégénérée à \(\Gamma\). Une séparation LO–TO nécessiterait un terme non analytique dépendant de la direction de \(\mathbf q\) au voisinage de \(\Gamma\).

#### Q4. Pourquoi V61 prédit-il \(v_{\mathrm{LA}} = v_{\mathrm{TA}}\) le long de \([100]\) ?

**[Observé numériquement]**

Dans V61, la régression acoustique le long de \([100]\) donne une dégénérescence des vitesses longitudinales et transverses : \(v_{\mathrm{LA}} = v_{\mathrm{TA}} \approx 2380.5\ \mathrm{m/s}\).

Cette égalité est une conséquence de la forme choisie pour les ressorts premiers voisins et de la symétrie du modèle. Elle ne doit pas être interprétée comme une prédiction expérimentale pour le CsCl réel.

**[Hypothèse]**

Une démonstration analytique complète nécessiterait d’exprimer les constantes élastiques effectives \(C_{11}\), \(C_{12}\) et \(C_{44}\) en fonction de \(K_{\parallel}\), \(K_{\perp}\), \(A\) et des masses, puis de vérifier les relations correspondantes entre vitesses et constantes.

### Constantes élastiques

#### Q5. Quel est le rapport de Zener du CsCl réel ?

**[Hypothèse — à ne pas affirmer sans source]**

Le rapport de Zener est défini par :

\[
A_Z = \frac{2C_{44}}{C_{11}-C_{12}}.
\]

Une valeur pour le CsCl réel ne doit pas être publiée sans une référence expérimentale vérifiée, une température précisée et une indication claire des constantes élastiques utilisées.

#### Q6. Quel rapport de Zener donne le modèle V61 ?

**[Observé numériquement]**

Pour le jeu de paramètres de V61, le modèle donne :

\[
A_Z \approx \frac{7}{3}.
\]

Cette valeur est un résultat du modèle, calculé à partir de ses constantes élastiques effectives. Elle ne constitue pas une mesure du CsCl réel.

#### Q7. Le signe de \(C'\) suffit-il à garantir la stabilité ?

**[Démontré partiellement]**

Non. La positivité de :

\[
C' = \frac{C_{11}-C_{12}}{2}
\]

est une condition nécessaire de stabilité élastique cubique, mais elle n’est pas suffisante à elle seule.

**[Observé numériquement]**

Pour V72–V74, la stabilité dynamique a été testée sur une grille finie de \(15^3 = 3375\) points de la zone de Brillouin. Aucune valeur propre inférieure à \(-10^{-5}\lambda_{\max}\) n’a été détectée.

**[Limite]**

Cela ne prouve pas la stabilité sur tout le continuum de la zone de Brillouin. Une validation plus complète exigerait un balayage plus dense, une recherche adaptative autour des minima de \(\lambda_{\min}(\mathbf q)\), et une vérification des conditions élastiques pertinentes.

#### Q8. Quelle est l’erreur physique du modèle par rapport au CsCl réel ?

**[Hypothèse — aucune valeur publiée]**

Aucune erreur physique quantitative n’est revendiquée. Une telle estimation exigerait des données expérimentales numérisées à partir d’Ahmad et al. (1972), avec leurs unités, directions cristallographiques, branches et incertitudes.

En l’absence de ces données, toute erreur de « plusieurs centaines de m/s » serait spéculative.

#### Q9. Pourquoi utiliser les seconds voisins ?

**[Hypothèse]**

Les seconds voisins Cs–Cs et Cl–Cl enrichissent le modèle en permettant aux deux sous-réseaux d’interagir avec leurs propres voisins. Ils peuvent modifier les pentes acoustiques et les fréquences optiques.

À \(\Gamma\), leur contribution aux blocs diagonaux s’annule, de sorte que la règle de somme acoustique reste satisfaite si elle l’était déjà dans le modèle à premiers voisins.

**[Observé numériquement]**

Avec le jeu de paramètres actuel, V72 donne des vitesses acoustiques distinctes :

\[
v_{\mathrm{TA}} = 2623.2 \pm 0.2\ \mathrm{m/s},
\]

\[
v_{\mathrm{LA}} = 3244.1 \pm 0.2\ \mathrm{m/s}.
\]

**[Limite]**

Les constantes de force des seconds voisins sont illustratives. Le jeu actuel conduit à un \(C_{12}\) négatif ; ce n’est donc pas un paramétrage physiquement définitif du CsCl.

### Fiabilité numérique

#### Q10. Comment les vitesses acoustiques sont-elles calculées ?

**[Démontré]**

Le long de \([100]\), on calcule \(\omega(q)\) pour plusieurs petites valeurs de \(q\), puis on ajuste une droite passant par l’origine :

\[
\omega(q) \approx vq.
\]

La pente \(v\) est la vitesse acoustique, en m/s lorsque \(\omega\) est en rad/s et \(q\) en m\(^{-1}\).

#### Q11. Que signifie l’incertitude \(\pm 0.2\) m/s ?

**[Observé numériquement]**

Il s’agit de l’erreur standard de la régression sur les points de simulation. Elle caractérise la dispersion des points autour de l’ajustement linéaire.

**[Limite]**

Elle n’inclut ni la courbure systématique de \(\omega(q)\), ni les choix d’intervalle en \(q\), ni une incertitude expérimentale. Ce n’est donc pas une incertitude physique totale.

#### Q12. Pourquoi utiliser un algorithme hongrois pour suivre les branches ?

**[Démontré]**

L’algorithme hongrois, ou Kuhn-Munkres, résout le problème d’affectation globale optimal. Il évite l’affectation gloutonne, où un mode peut prendre le meilleur correspondant d’un autre mode.

**[Limite]**

L’algorithme optimise l’affectation relativement à la matrice de scores choisie, par exemple un recouvrement de vecteurs propres ou de sous-espaces. Il ne garantit pas, à lui seul, un suivi physique correct en présence de dégénérescences.

#### Q13. Pourquoi comparer des sous-espaces plutôt que des vecteurs propres ?

**[Démontré]**

Lorsque plusieurs modes sont dégénérés, leurs vecteurs propres individuels ne sont pas uniques : toute rotation unitaire à l’intérieur du sous-espace décrit la même physique.

Comparer des projecteurs de sous-espaces,

\[
P = \sum_i |v_i\rangle\langle v_i|,
\]

est donc plus robuste que comparer vecteur par vecteur.

**[Hypothèse]**

Un alignement de type Procrustes peut aider à construire une base continue dans un sous-espace dégénéré, mais il ne garantit pas à lui seul un suivi correct des branches. Il doit être combiné à une métrique de recouvrement et à des tests de continuité.

#### Q14. Les forces angulaires amélioreraient-elles le modèle ?

**[Hypothèse]**

Oui, potentiellement. Des termes à trois corps ou des forces angulaires peuvent enrichir le modèle au-delà des simples paires de ressorts, notamment pour représenter plus finement les déformations de la structure.

**[Limite]**

La relation de Cauchy entre constantes élastiques dépend d’hypothèses précises sur les interactions, la symétrie et l’équilibre du réseau. La violation de cette relation ne peut pas être présentée comme une conséquence automatique de l’ajout d’un terme à trois corps.

## 📌 Niveau de preuve

| Niveau | Signification |
|---|---|
| **[Démontré]** | Propriété mathématique ou test reproductible du code |
| **[Démontré partiellement]** | Condition nécessaire vérifiée, sans preuve de suffisance |
| **[Observé numériquement]** | Résultat obtenu pour un jeu de paramètres et une grille donnés |
| **[Hypothèse]** | Piste physique, extension ou interprétation non encore validée |
| **[Limite]** | Restriction explicite de la portée d’un résultat ou d’un test |

---

Développé sur Termux et VS Code.
