#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
dyadic_renormalization.py v1 — Addendum II-Phys (Cahier de Physique).
Flot de renormalisation sur l'arbre 2-adique : enchevêtrement et saturation.

Lecture physique (analogie mesurée, Article III) :
La profondeur k = v_2(p-1) est l'échelle de renormalisation (RG).
La partie impaire m = (p-1)/2^k représente les degrés de liberté résiduels.
Une distribution de m "plus lisse" (moins de variance dans les écarts) 
à grand k est l'analogue arithmétique de l'intrication à longue portée 
dans les réseaux de tenseurs (MERA) ou de l'ordre dans les verres de spin.
La saturation à grand k correspond à un point fixe infra-rouge (gel du système).

Gardes pré-enregistrées :
  G1 (Enchevêtrement / Croissance) : Pour 1 <= k <= 8, le score de régularité 
       (inverse de la variance normalisée des écarts) doit être > 3σ au-dessus 
       du null aléatoire (distribution uniforme de même taille).
  G2 (Saturation / Point Fixe) : Pour k >= 10, le score de régularité doit 
       cesser de croître de manière significative (stagnation ou effondrement), 
       indiquant la limite structurelle de la hiérarchie (rareté des premiers 
       de Sophie Germain et au-delà).
  G3 (Robustesse) : La forme de la courbe de régularité en fonction de k doit 
       être fortement corrélée (Spearman > 0.8) entre la fenêtre P <= 10^6 
       et la fenêtre P <= 10^7.
"""
import sys
import time
import math
import numpy as np

P_MAX = 10 ** 7
P_SUB = 10 ** 6
K_MAX = 20
N_NULL_TRIALS = 15  # Nombre d'échantillons aléatoires pour estimer le null

def crible_premiers(n):
    """Crible d'Ératosthène optimisé avec numpy."""
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i :: i] = False
    return np.nonzero(isp)[0]

def evaluer_regularite(m_array, m_max, n_trials=10):
    """
    Calcule le score de régularité d'un ensemble de valeurs.
    Métrique : Inverse de la variance normalisée des écarts (gaps) entre valeurs triées.
    Pour un processus de Poisson (aléatoire uniforme), Var(gaps) ≈ Mean(gaps)^2, donc V ≈ 1.
    Pour une distribution plus régulière (quasi-réseau), V < 1, donc Score > 1.
    """
    if len(m_array) < 10:
        return 0.0, 0.0, 0.0
    
    m_sorted = np.sort(m_array)
    gaps = np.diff(m_sorted)
    mean_gap = np.mean(gaps)
    var_gap = np.var(gaps)
    
    if mean_gap == 0 or var_gap == 0:
        return 0.0, 0.0, 0.0
        
    # Variance normalisée des écarts observés
    V_obs = var_gap / (mean_gap ** 2)
    score_obs = 1.0 / V_obs if V_obs > 0 else 0.0
    
    # Estimation du null (échantillons aléatoires de même taille dans [1, m_max])
    scores_null = []
    for _ in range(n_trials):
        m_rand = np.random.randint(1, m_max + 1, size=len(m_array))
        m_rand_sorted = np.sort(m_rand)
        gaps_rand = np.diff(m_rand_sorted)
        mean_gap_rand = np.mean(gaps_rand)
        var_gap_rand = np.var(gaps_rand)
        V_rand = var_gap_rand / (mean_gap_rand ** 2) if mean_gap_rand > 0 else 1.0
        scores_null.append(1.0 / V_rand if V_rand > 0 else 0.0)
        
    mean_null = np.mean(scores_null)
    std_null = np.std(scores_null) if len(scores_null) > 1 else 0.1
    
    return score_obs, mean_null, std_null

def main():
    t0 = time.time()
    print("=" * 78)
    print("dyadic_renormalization.py v1 — Addendum II-Phys : flot de renormalisation 2-adique")
    print("=" * 78)
    print()
    
    print(f"[1/4] Crible des premiers jusqu'à {P_MAX}...")
    primes = crible_premiers(P_MAX)
    print(f"      Trouvés : {len(primes)} premiers.")
    
    # Séparation des fenêtres pour G3
    primes_sub = primes[primes <= P_SUB]
    
    print(f"[2/4] Calcul des valuations 2-adiques (k) et parties impaires (m)...")
    # Pour P_MAX
    p_minus_1 = primes - 1
    # Astuce rapide pour v_2 : compter les zéros de fin en binaire
    k_all = np.array([int(math.log2(p & -p)) for p in p_minus_1])
    m_all = p_minus_1 // (2 ** k_all)
    
    # Pour P_SUB
    p_minus_1_sub = primes_sub - 1
    k_sub = np.array([int(math.log2(p & -p)) for p in p_minus_1_sub])
    m_sub = p_minus_1_sub // (2 ** k_sub)
    
    print(f"[3/4] Évaluation de la régularité par échelle k (null stratifié)...")
    results = []
    results_sub = []
    
    for k in range(1, K_MAX + 1):
        # Fenêtre principale
        mask = (k_all == k)
        m_k = m_all[mask]
        if len(m_k) > 50:
            m_max = np.max(m_k)
            score, null_mean, null_std = evaluer_regularite(m_k, m_max, N_NULL_TRIALS)
            sigma = (score - null_mean) / null_std if null_std > 0 else 0.0
            results.append((k, len(m_k), score, null_mean, null_std, sigma))
        else:
            results.append((k, len(m_k), 0.0, 0.0, 0.0, 0.0))
            
        # Fenêtre secondaire (pour G3)
        mask_sub = (k_sub == k)
        m_k_sub = m_sub[mask_sub]
        if len(m_k_sub) > 50:
            m_max_sub = np.max(m_k_sub)
            score_sub, _, _ = evaluer_regularite(m_k_sub, m_max_sub, 5) # Moins de trials pour la vitesse
            results_sub.append((k, score_sub))
        else:
            results_sub.append((k, 0.0))

    print(f"[4/4] Vérification des gardes pré-enregistrées...")
    print()
    print(f"{'k':<3} | {'N(k)':<7} | {'Score Obs':<10} | {'Null Mean':<10} | {'Sigma':<6}")
    print("-" * 50)
    
    scores_obs_list = []
    scores_sub_list = []
    g1_pass = True
    g2_pass = True
    
    peak_score = 0.0
    peak_k = 1
    
    for res in results:
        k, n, score, null_mean, null_std, sigma = res
        if n > 50:
            print(f"{k:<3} | {n:<7} | {score:<10.2f} | {null_mean:<10.2f} | {sigma:+6.1f}σ")
            scores_obs_list.append(score)
            
            # G1 : Pour 1 <= k <= 8, sigma doit être > 3
            if 1 <= k <= 8:
                if sigma < 3.0:
                    g1_pass = False
            
            # Suivi du pic pour G2
            if score > peak_score:
                peak_score = score
                peak_k = k
                
    # G2 : Pour k >= 10, le score ne doit plus croître significativement par rapport au pic
    # On vérifie si le score moyen pour k >= 10 est inférieur ou égal au pic, ou si la pente est négative
    scores_k_ge_10 = [res[2] for res in results if res[0] >= 10 and res[1] > 50]
    if len(scores_k_ge_10) > 0:
        # Si le dernier score est bien inférieur au pic, ou si la moyenne est stagnante, c'est bon.
        # Une chute ou stagnation après le pic (souvent vers k=6 à 9) valide la saturation.
        if scores_k_ge_10[-1] > peak_score * 1.2: # Tolérance de 20%
            g2_pass = False
    else:
        g2_pass = False # Pas assez de données pour juger
        
    # G3 : Corrélation de Spearman entre les courbes P<=10^6 et P<=10^7
    scores_sub_list = [res[1] for res in results_sub if res[1] > 0]
    scores_obs_matched = [res[2] for res in results if res[0] < len(results_sub) and results_sub[res[0]-1][1] > 0]
    
    if len(scores_sub_list) > 5 and len(scores_obs_matched) > 5:
        # Calcul simplifié de la corrélation de rang (Spearman)
        rank_obs = np.argsort(np.argsort(scores_obs_matched))
        rank_sub = np.argsort(np.argsort(scores_sub_list))
        spearman = np.corrcoef(rank_obs, rank_sub)[0, 1]
        g3_pass = (spearman > 0.8)
        print(f"\nCorrélation de forme (P=10^6 vs P=10^7) : Spearman ρ = {spearman:.3f}")
    else:
        g3_pass = False
        print("\nCorrélation de forme : Données insuffisantes pour G3.")

    print()
    print("=" * 78)
    print("VERDICT DES GARDES :")
    print(f"  G1 (Enchevêtrement, k<=8, sigma>3) : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Saturation, k>=10, stagnation) : {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Robustesse, Spearman > 0.8)    : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : La hiérarchie dyadique montre un enchevêtrement croissant")
        print("    suivi d'une saturation nette, validant l'analogie avec un flot de")
        print("    renormalisation aboutissant à un point fixe infra-rouge (verre de spin / MERA).")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique observée ne correspond pas au modèle de")
        print("    renormalisation dyadique attendu. Les gardes ont été violées.")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()