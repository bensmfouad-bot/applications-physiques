# Cahier de Physique I : Laboratoire Numérique
**Applications physiques du crible congruentiel : vide adélique, renormalisation dyadique, chaos quantique, matière condensée.**

**Auteur :** Fouad Bensmail  
**Date :** Septembre 2026  
**Preuve d'antériorité et d'intégrité :** [Zenodo DOI: 10.5281/zenodo.22950183](https://doi.org/10.5281/zenodo.22950183)

---

## 📖 À propos de ce dépôt
Ce dépôt contient les compagnons de mesure (scripts Python) et les documents associés au *Cahier de Physique I*. Il s'agit d'un laboratoire numérique rigoureux qui applique les mesures arithmétiques (volumes locaux, résidus spectraux, grammaire dyadique) à quatre portes physiques, selon une méthode stricte héritée du cahier mathématique : *toute affirmation est une mesure avec garde ; tout refus est publié*.

## 🔬 Les Compagnons de Mesure
- **`adelic_amplitude.py` (Rang I)** : Mesure du vide adélique, du flot de renormalisation (RG) et de la séparation géométrie/matière. Audit Vert.
- **`dyadic_renormalization.py` (Rang II)** : Flot de renormalisation sur l'arbre 2-adique, enchevêtrement et saturation hiérarchique (verres de spin, réseaux MERA).
- **`quantum_scars.py` (Rang III)** : Analyse spectrale des écarts dépliés et recherche de résonances transitoires (états de Gamow, chaos quantique).
- **`prime_percolation.py` (Rang IV)** : Test de percolation arithmétique et transitions de phase sur la chaîne 1D des premiers.

## 🏭 Note de Transfert I : Audit des PRNG et Sécurité Matérielle
Au-delà de la recherche fondamentale, ce dépôt inclut une **Note de Transfert I** détaillant un protocole opérationnel pour l'industrie de la sécurité matérielle (HSM, FPGA, cartes à puce, CI/CD) :
- Détection de biais structurels profonds dans les générateurs pseudo-aléatoires (PRNG/LCG).
- Méthodologie du **Null Stratifié** et calibration par **surrogats de phase** pour éliminer les faux positifs institutionnels.
- Identification des vulnérabilités de type "héritage de RANDU" (biais de parité, lissité).

*Pour toute demande d'utilisation commerciale, d'intégration industrielle ou de partenariat, veuillez contacter l'auteur. L'utilisation commerciale requiert une licence spécifique (ce dépôt est sous licence CC BY-NC 4.0).*

## 🛡️ Méthodologie (Articles II & III de la Charte)
1. **Gardes pré-enregistrées** : Les seuils et métriques sont fixés avant toute mesure.
2. **Calibration par surrogats** : Utilisation de surrogats de phase pour un *null stratifié* réaliste, propre à la dynamique du système.
3. **Publication des refus** : Les limites des instruments et des analogies physiques sont cartographiées et documentées avec la même rigueur que les succès.

## ⚖️ Licence et Propriété Intellectuelle
**Copyright (c) 2026 Fouad Bensmail — Tous droits réservés.**  
Licence : [Creative Commons Attribution - Pas d'Utilisation Commerciale 4.0 International (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

La paternité de cette œuvre est légalement et cryptographiquement attestée par les dépôts Zenodo et GitHub. Toute reproduction, utilisation commerciale ou appropriation sans l'autorisation écrite explicite de l'auteur est strictement interdite.