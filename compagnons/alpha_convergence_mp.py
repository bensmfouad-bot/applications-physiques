#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0
# =============================================================================
"""
alpha_convergence_mp.py — Mesure de haute précision de la décroissance de la queue eulérienne.
Utilise mpmath (50 décimales) pour éliminer les artefacts d'arrondi de l'ancre.
"""
import mpmath
import time

# 1. Configurer la haute précision (50 décimales)
mpmath.mp.dps = 50

def sieve(n):
    """Crible d'Ératosthène standard pour obtenir la liste des premiers."""
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i, x in enumerate(s) if x]

def theta_partial_mp(pattern, primes, P_max):
    """Calcule le produit partiel Θ_P avec une précision arbitraire."""
    k = len(pattern)
    prod = mpmath.mpf(1.0)
    for p in primes:
        if p > P_max:
            break
        nu = len(set(h % p for h in pattern))
        if nu >= p:
            return mpmath.mpf(0.0)
        
        # Calcul en haute précision
        term1 = 1.0 - mpmath.mpf(nu) / p
        term2 = (1.0 - mpmath.mpf(1.0) / p) ** (-k)
        prod *= term1 * term2
    return prod

print("=" * 85)
print("alpha_convergence_mp.py — Mesure de haute précision de la queue eulérienne")
print("=" * 85)

print("\n[1/3] Génération des premiers jusqu'à 10^7...")
t0 = time.time()
primes = sieve(10_000_000)
print(f"      Trouvés : {len(primes)} premiers en {time.time()-t0:.2f} s")

patterns = {
    "jumeaux":    [0, 2],
    "triplet_A":  [0, 2, 6],
    "quadruplet": [0, 2, 6, 8],
}

# 2. Calcul de l'ancre de haute précision à P = 10^7
print("\n[2/3] Calcul des ancres de haute précision (P = 10^7)...")
Theta_inf_mp = {}
for name, pat in patterns.items():
    Theta_inf_mp[name] = theta_partial_mp(pat, primes, 10_000_000)
    print(f"      {name:12s} : Θ_10^7 = {Theta_inf_mp[name]}")

# 3. Mesure de la décroissance et calcul de α_eff
print("\n[3/3] Mesure de la décroissance de l'erreur et calcul de α_eff...")
print("-" * 85)

for name, pat in patterns.items():
    print(f"\n=== {name.upper()} (k={len(pat)}) ===")
    row = []
    
    # On teste les échelles inférieures
    for P in [10**3, 10**4, 10**5, 10**6]:
        tp = theta_partial_mp(pat, primes, P)
        # La différence est calculée en haute précision
        d = abs(Theta_inf_mp[name] - tp)
        row.append((P, tp, d))
    
    for i in range(1, len(row)):
        P0, _, d0 = row[i-1]
        P1, tp1, d1 = row[i]
        
        if d0 > 0 and d1 > 0:
            # α_eff = - ln(d1 / d0) / ln(P1 / P0)
            alpha = -mpmath.log(d1 / d0) / mpmath.log(P1 / P0)
            
            # Théorie : si d ~ 1/(P ln P), alors α_eff ≈ 1 + 1/ln(P1)
            alpha_theo = 1.0 + 1.0 / mpmath.log(P1)
            
            print(f"  P={P1:<7} | d={float(d1):<14.3e} | α_eff={float(alpha):<7.4f} | α_théo={float            (alpha_theo):<7.4f}")

print("\n" + "=" * 85)
print("[OK] Course terminée.")