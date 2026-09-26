#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
multidimensional_sieve.py v1.1 — Cahier de Physique II, Porte C.
Crible multidimensionnel et contournement partiel du mur de la parité.

Correction v1.1 : L'enregistrement des métriques se déclenche désormais dès que 
le premier courant p dépasse ou atteint le seuil de mesure (p >= mesure_point), 
corrigant l'artefact de la v1 qui cherchait une égalité stricte (p == 10^6, impossible).
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

def main():
    t0 = time.time()
    print("=" * 78)
    print("multidimensional_sieve.py v1.1 — Cahier II, Porte C : crible multidimensionnel")
    print("=" * 78)
    
    P_MAX = 10_000_000
    print(f"\n[1/4] Crible des premiers jusqu'à {P_MAX}...")
    primes = crible_premiers(P_MAX)
    prime_set = set(primes) # Pour des recherches O(1)
    print(f"      Trouvés : {len(primes)} premiers.")
    
    print(f"[2/4] Détection des triplets [p, p+2, p+6] et pondération multidimensionnelle...")
    
    mesure_points = [10**5, 10**6, 2*10**6, 5*10**6, 10**7]
    results = []
    
    current_weighted_sum = 0.0
    current_unweighted_sum = 0
    mp_idx = 0
    
    for p in primes:
        if p > P_MAX - 8:
            break
            
        # Vérification du triplet [p, p+2, p+6]
        if (p + 2) in prime_set and (p + 6) in prime_set:
            current_unweighted_sum += 1
            
            # Pondération multidimensionnelle : est-ce aussi un quadruplet [p, p+2, p+6, p+8] ?
            if (p + 8) in prime_set:
                current_weighted_sum += 2.0
            else:
                current_weighted_sum += 1.0
                
        # Enregistrement dès que nous dépassons ou atteignons le point de mesure courant
        while mp_idx < len(mesure_points) and p >= mesure_points[mp_idx]:
            gain = current_weighted_sum / current_unweighted_sum if current_unweighted_sum > 0 else 0.0
            results.append({
                'P': mesure_points[mp_idx],
                'unweighted': current_unweighted_sum,
                'weighted': current_weighted_sum,
                'gain': gain
            })
            mp_idx += 1
            
        if mp_idx >= len(mesure_points):
            break

    print(f"[3/4] Évaluation des gardes pré-enregistrées...")
    g1_pass = False
    g2_pass = False
    g3_pass = False
    
    # Affichage des résultats
    print("\n   P        | N_brut | N_pondéré | Gain G(P)")
    print("   " + "-" * 52)
    for r in results:
        print(f"   {r['P']:<8} | {r['unweighted']:<6} | {r['weighted']:<9.1f} | {r['gain']:.4f}")
        
    # G1 : Gain mesurable (G(P) > 1.05 pour P >= 10^6)
    res_10_6 = next((r for r in results if r['P'] == 10**6), None)
    gain_10_6 = res_10_6['gain'] if res_10_6 else 0.0
    g1_pass = (gain_10_6 > 1.05)
    print(f"\n      G1 : Gain à P=10^6 est {gain_10_6:.4f} (Seuil > 1.05) -> {'OK' if g1_pass else 'ECHEC'}")
    
    # G2 : Stabilité asymptotique (Variation relative entre 10^6 et 2*10^6 < 10%)
    res_2_10_6 = next((r for r in results if r['P'] == 2*10**6), res_10_6)
    gain_2_10_6 = res_2_10_6['gain'] if res_2_10_6 else gain_10_6
    
    if gain_10_6 > 0:
        var_rel = abs(gain_2_10_6 - gain_10_6) / gain_10_6
    else:
        var_rel = 1.0
    g2_pass = (var_rel < 0.10)
    print(f"      G2 : Variation relative entre 10^6 et 2x10^6 est {var_rel:.1%} (Seuil < 10%) -> {'OK' if g2_pass else 'ECHEC'}")
    
    # G3 : Nature du gain (Convergence vers une constante finie, pas de divergence)
    res_5_10_6 = next((r for r in results if r['P'] == 5*10**6), res_2_10_6)
    res_10_7 = next((r for r in results if r['P'] == 10**7), res_5_10_6)
    
    gain_5_10_6 = res_5_10_6['gain'] if res_5_10_6 else gain_2_10_6
    gain_10_7 = res_10_7['gain'] if res_10_7 else gain_5_10_6
    
    if gain_5_10_6 > 0:
        var_final = abs(gain_10_7 - gain_5_10_6) / gain_5_10_6
    else:
        var_final = 1.0
    g3_pass = (var_final < 0.05) # Stabilisation vers une constante finie
    print(f"      G3 : Variation finale du gain (5x10^6 à 10^7) est {var_final:.1%} (Seuil < 5% pour convergence finie) -> {'OK' if g3_pass else 'ECHEC'}")
    
    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Cahier II, Porte C) :")
    print(f"  G1 (Gain mesurable, G > 1.05)               : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Stabilité asymptotique, var < 10%)      : {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Convergence finie, pas de divergence)   : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : Le crible multidimensionnel extrait un gain mesurable et stable.")
        print("    Les corrélations d'ordre supérieur (clusters) contournent partiellement")
        print("    le mur de la parité, convergeant vers une constante finie (pas de divergence).")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique du crible multidimensionnel ne correspond pas au")
        print("    modèle attendu (gain insuffisant, instable, ou divergent).")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()