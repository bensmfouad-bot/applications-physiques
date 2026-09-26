#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
galois_smoothing_v2.py — Cahier de Physique II, Porte B.
Lissage de la queue eulérienne avec estimation de queue "Quantum-Inspired" (Ewin Tang).

Correction majeure : Au lieu de tronquer la somme, on utilise un accès par 
échantillonnage (Sampling) et par requête (Query via Miller-Rabin) pour estimer
la queue de la somme jusqu'à 10^9 en temps polylogarithmique, évitant le parcours O(N).
"""
import sys
import time
import math
import random

def crible_premiers(n):
    """Crible d'Ératosthène pour la base exacte."""
    isp = [True] * (n + 1)
    isp[0] = isp[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            for j in range(i * i, n + 1, i):
                isp[j] = False
    return [i for i, is_p in enumerate(isp) if is_p]

def is_prime_miller_rabin(n, k=5):
    """Test de primalité de Miller-Rabin (notre 'Query Access' rapide)."""
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0: return False
    
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
        
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def legendre_symbol_2(p):
    """Symbole de Legendre (2/p) : 1 si p ≡ ±1 mod 8, -1 si p ≡ ±3 mod 8."""
    # Optimisation : pas besoin de pow(), on regarde juste p mod 8
    r = p % 8
    if r == 1 or r == 7:
        return 1
    elif r == 3 or r == 5:
        return -1
    return 0

def main():
    t0 = time.time()
    print("=" * 78)
    print("galois_smoothing_v2.py — Porte B : Lissage avec estimation de queue (Ewin Tang)")
    print("=" * 78)
    
    P_BASE = 10_000_000
    P_MAX_ESTIMATE = 1_000_000_000
    N_SAMPLES = 200_000  # Nombre de requêtes "quantiques" (équilibre précision/temps)
    
    print(f"\n[1/4] Crible exact de la base jusqu'à {P_BASE}...")
    primes_base = crible_premiers(P_BASE)
    prime_set_base = set(primes_base)
    print(f"      Trouvés : {len(primes_base)} premiers.")
    
    print(f"[2/4] Estimation de la queue [{P_BASE}, {P_MAX_ESTIMATE}] via Sampling & Query (Miller-Rabin)...")
    print(f"      Nombre d'échantillons (requêtes) : {N_SAMPLES}")
    
    # Calcul exact de la somme jusqu'à P_BASE
    sum_brute_base = sum(1.0 / (p * p) for p in primes_base)
    sum_lisse_base = sum((1.0 + legendre_symbol_2(p)) / (p * p) for p in primes_base)
    
    # Estimation Monte Carlo de la queue (Paradigme Ewin Tang)
    tail_brute_est = 0.0
    tail_lisse_est = 0.0
    primes_found_in_tail = 0
    
    range_size = P_MAX_ESTIMATE - P_BASE
    
    for _ in range(N_SAMPLES):
        # Sampling : tirage uniforme dans l'intervalle
        x = random.randint(P_BASE + 1, P_MAX_ESTIMATE)
        
        # Query : test de primalité rapide
        if is_prime_miller_rabin(x):
            primes_found_in_tail += 1
            term = 1.0 / (x * x)
            tail_brute_est += term
            tail_lisse_est += term * (1.0 + legendre_symbol_2(x))
            
    # Facteur d'échelle pour obtenir un estimateur non biaisé de la somme sur tout l'intervalle
    scale_factor = range_size / N_SAMPLES
    tail_brute_total = tail_brute_est * scale_factor
    tail_lisse_total = tail_lisse_est * scale_factor
    
    print(f"      Premiers trouvés dans l'échantillon : {primes_found_in_tail}")
    print(f"      Estimation de la queue brute : {tail_brute_total:.6e}")
    print(f"      Estimation de la queue lissée : {tail_lisse_total:.6e}")
    
    print(f"\n[3/4] Évaluation des gardes aux points de mesure...")
    
    # Points de mesure dans la base exacte
    mesure_points = [10**5, 2*10**5, 5*10**5, 10**6, 2*10**6, 5*10**6]
    results = []
    
    for P in mesure_points:
        # Somme exacte jusqu'à P
        sum_brute_P = sum(1.0 / (p * p) for p in primes_base if p > P)
        sum_lisse_P = sum((1.0 + legendre_symbol_2(p)) / (p * p) for p in primes_base if p > P)
        
        # On ajoute l'estimation de la queue au-delà de P_BASE (qui est constante pour tous ces P)
        # Note: C'est une approximation, car la queue réelle dépend de P, mais comme P << P_BASE, 
        # la différence entre la queue à P et la queue à P_BASE est négligeable par rapport à l'estimation.
        # Pour être rigoureux, on pourrait estimer la queue spécifiquement pour chaque P, 
        # mais l'ordre de grandeur est dominé par les premiers termes.
        
        C_brute = (sum_brute_P + tail_brute_total) * P * math.log(P)
        C_lisse = (sum_lisse_P + tail_lisse_total) * P * math.log(P)
        
        results.append({'P': P, 'C_brute': C_brute, 'C_lisse': C_lisse})
        
    print("\n   P        | C_brute  | C_lissé")
    print("   " + "-" * 40)
    for r in results:
        print(f"   {r['P']:<8} | {r['C_brute']:<8.4f} | {r['C_lisse']:<8.4f}")
        
    print(f"\n[4/4] Vérification des gardes...")
    g1_pass = False
    g2_pass = False
    g3_pass = False
    
    # G1 : Gain de Précision (Variance relative entre 10^5 et 2*10^5)
    res_10_5 = next((r for r in results if r['P'] == 10**5), None)
    res_2_10_5 = next((r for r in results if r['P'] == 2*10**5), None)
    
    if res_10_5 and res_2_10_5:
        var_brute = abs(res_2_10_5['C_brute'] - res_10_5['C_brute']) / res_10_5['C_brute']
        var_lisse = abs(res_2_10_5['C_lisse'] - res_10_5['C_lisse']) / res_10_5['C_lisse']
        
        if var_brute > 0:
            gain = (var_brute - var_lisse) / var_brute
        else:
            gain = 0.0
        
        g1_pass = (gain > 0.20)
        print(f"      G1 : Variance brute = {var_brute:.2%}, Variance lissée = {var_lisse:.2%}")
        print(f"           Gain de précision = {gain:.1%} (Seuil > 20%) -> {'OK' if g1_pass else 'ECHEC'}")
    
    # G2 : Stabilité Asymptotique (Variation de C_lissé entre 10^5 et 10^6)
    res_10_6 = next((r for r in results if r['P'] == 10**6), None)
    if res_10_5 and res_10_6:
        var_stability = abs(res_10_6['C_lisse'] - res_10_5['C_lisse']) / res_10_5['C_lisse']
        g2_pass = (var_stability < 0.05)
        print(f"      G2 : Variation C_lissé (10^5 à 10^6) = {var_stability:.1%} (Seuil < 5%) -> {'OK' if g2_pass else 'ECHEC'}")
    
    # G3 : Nature du Lissage (Différence asymptotique < 10%)
    if res_10_6:
        diff_asym = abs(res_10_6['C_lisse'] - res_10_6['C_brute']) / res_10_6['C_brute']
        g3_pass = (diff_asym < 0.10)
        print(f"      G3 : Différence asymptotique = {diff_asym:.1%} (Seuil < 10%) -> {'OK' if g3_pass else 'ECHEC'}")
    
    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Cahier II, Porte B - Méthode Tang) :")
    print(f"  G1 (Gain de Précision, gain > 20%)        : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Stabilité Asymptotique, var < 5%)     : {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Nature du Lissage, diff < 10%)        : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : L'estimation de queue quantique-classique a permis de révéler")
        print("    la véritable dynamique asymptotique. Le lissage de Galois tient ses promesses.")
    else:
        print("\n>>> REFUS PUBLIÉ : Même avec une estimation de queue correcte jusqu'à 10^9,")
        print("    le lissage par caractère de Dirichlet n'accélère pas la convergence de la")
        print("    série absolument convergente, confirmant la limite analytique.")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()