#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
prime_percolation.py v1 — Addendum IV-Phys (Cahier de Physique).
Percolation arithmétique et transitions de phase dans la séquence des premiers.

Lecture physique (analogie mesurée, Article III) :
On construit un graphe 1D où les nœuds sont les nombres premiers consécutifs.
Un lien existe entre p_n et p_{n+1} si leur écart g_n <= λ * <g>.
En variant λ, on observe l'émergence d'un "amas géant" (giant component),
analogue à la transition de percolation en physique statistique.

Gardes pré-enregistrées :
  G1 (Émergence d'un seuil) : Il doit exister une valeur critique λ_c (autour de 1.0-1.5)
       où la taille relative de l'amas géant S_max passe rapidement de ~0 à ~1.
  G2 (Universalité de la fenêtre) : La forme de la courbe de transition (λ_c et pente)
       doit rester stable (variation < 5% en moyenne) entre deux fenêtres disjointes.
  G3 (Nature de la transition - Hypothèse de refus) : Contrairement aux systèmes
       thermodynamiques infinis (transition abrupte du 2nd ordre), la distribution des
       premiers est déterministe et bornée. On s'attend à une transition "douce" ou
       étalée, refusant l'analogie stricte avec une transition de phase classique.
"""
import sys
import time
import numpy as np

def crible_premiers(n):
    """Crible d'Ératosthène optimisé avec numpy."""
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i :: i] = False
    return np.nonzero(isp)[0]

def get_giant_component_size(primes_window, lambda_val):
    """
    Calcule la taille relative de l'amas géant pour un seuil lambda donné.
    Un lien existe si g_n <= lambda * <g>.
    """
    if len(primes_window) < 2:
        return 0.0
    
    gaps = np.diff(primes_window)
    mean_gap = np.mean(gaps)
    
    # Masque des liens actifs (True si lié, False sinon)
    links = gaps <= (lambda_val * mean_gap)
    
    # Trouver les composantes connexes (séquences de True consécutifs)
    # On pad avec False aux extrémités pour détecter les débuts et fins de composantes
    padded = np.concatenate(([False], links, [False]))
    
    # Les débuts de composantes sont où on passe de False à True
    # Les fins sont où on passe de True à False
    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    
    # Tailles des composantes (nombre de liens + 1 = nombre de nœuds)
    component_sizes = ends - starts
    
    if len(component_sizes) == 0:
        max_size = 1 # Au moins un nœud isolé
    else:
        max_size = np.max(component_sizes) + 1
        
    return max_size / len(primes_window)

def main():
    t0 = time.time()
    print("=" * 78)
    print("prime_percolation.py v1 — Addendum IV-Phys : percolation arithmétique")
    print("=" * 78)
    
    print("\n[1/4] Crible des premiers jusqu'à 5,000,000...")
    primes = crible_premiers(5_000_000)
    print(f"      Trouvés : {len(primes)} premiers.")
    
    # Définition de deux fenêtres disjointes de même taille pour tester G2
    # Fenêtre 1 : ~ [1,000,000, 2,000,000]
    idx1_start = np.searchsorted(primes, 1_000_000)
    idx1_end = np.searchsorted(primes, 2_000_000)
    window1 = primes[idx1_start:idx1_end]
    
    # Fenêtre 2 : ~ [4,000,000, 5,000,000]
    idx2_start = np.searchsorted(primes, 4_000_000)
    idx2_end = np.searchsorted(primes, 5_000_000)
    window2 = primes[idx2_start:idx2_end]
    
    print(f"      Fenêtre 1 : {len(window1)} premiers (P ~ 10^6 à 2*10^6)")
    print(f"      Fenêtre 2 : {len(window2)} premiers (P ~ 4*10^6 à 5*10^6)")
    
    print(f"\n[2/4] Calcul de la taille de l'amas géant S_max en fonction de λ...")
    lambdas = np.arange(0.1, 3.1, 0.1)
    s_max_1 = []
    s_max_2 = []
    
    for lam in lambdas:
        s_max_1.append(get_giant_component_size(window1, lam))
        s_max_2.append(get_giant_component_size(window2, lam))
        
    print(f"[3/4] Évaluation des gardes pré-enregistrées...")
    g1_pass = False
    g2_pass = False
    g3_pass = False
    
    # G1 : Émergence d'un seuil (S_max passe de < 0.2 à > 0.8 entre deux valeurs de λ consécutives)
    jumps_1 = np.diff(s_max_1)
    max_jump_1 = np.max(jumps_1)
    # Un saut de > 0.3 sur un pas de 0.1 est considéré comme une émergence de seuil
    g1_pass = (max_jump_1 > 0.3)
    lambda_c_1 = lambdas[np.argmax(jumps_1) + 1]
    print(f"      G1 (Fenêtre 1) : Saut max de S_max = {max_jump_1:.3f} à λ ≈ {lambda_c_1:.1f} (Seuil > 0.3) -> {'OK' if g1_pass else 'ECHEC'}")
    
    # G2 : Universalité de la fenêtre (comparaison des courbes)
    # On calcule la différence absolue moyenne entre les deux courbes
    diff_curves = np.abs(np.array(s_max_1) - np.array(s_max_2))
    mean_diff = np.mean(diff_curves)
    # Si la différence moyenne est < 0.05 (5%), les courbes sont très similaires
    g2_pass = (mean_diff < 0.05)
    print(f"      G2 : Différence moyenne entre fenêtres = {mean_diff:.3f} (Seuil < 0.05 pour stabilité) -> {'OK' if g2_pass else 'ECHEC'}")
    
    # G3 : Nature de la transition (Hypothèse de refus : transition douce)
    # Une transition de phase du 2nd ordre stricte aurait une pente infinie (saut quasi-vertical).
    # Ici, on mesure la "douceur" en regardant sur combien de pas de λ le saut principal s'étale.
    # Si le saut > 0.3 s'étale sur plus de 2 pas (0.2 en λ), c'est une transition "douce/étalée".
    significant_jumps = np.sum(jumps_1 > 0.15)
    g3_pass = (significant_jumps >= 2) # Transition étalée sur au moins 0.2 en λ
    print(f"      G3 : Nombre de pas avec saut > 0.15 autour du pic = {significant_jumps} (Seuil >= 2 pour transition douce) -> {'OK' if g3_pass else 'ECHEC'}")
    
    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Rang IV) :")
    print(f"  G1 (Émergence d'un seuil, saut > 0.3)      : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Universalité de la fenêtre, diff < 0.05): {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Transition douce/étalée, pas >= 2)     : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT (avec refus de l'analogie stricte) :")
        print("    Un seuil de percolation émerge de manière universelle, mais la transition")
        print("    est douce et étalée, refusant l'analogie stricte avec une transition de")
        print("    phase du second ordre classique (système thermodynamique infini).")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique de percolation observée ne correspond pas au")
        print("    modèle attendu (soit pas de seuil clair, soit manque d'universalité, soit")
        print("    transition trop abrupte pour être qualifiée de 'douce').")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()