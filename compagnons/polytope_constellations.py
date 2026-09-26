#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
polytope_constellations.py v1 — Cahier de Physique II, Porte A.
Géométrie des polytopes croisés et invariance fractale des constellations.

Lecture physique (analogie mesurée) :
Une constellation de premiers n'est pas qu'une suite 1D, c'est la projection 
d'un objet géométrique dans l'espace des résidus modulo M_k. L'axiome d'invariance 
stipule que le nombre de positions d'ancrage valides V_k est donné exactement 
par le produit géométrique : V_k = ∏_{p|M_k} (p - ν_H(p)).

Gardes pré-enregistrées :
  G1 (Exactitude Combinatoire) : Le comptage brutal des résidus valides modulo M_k 
       doit égaler exactement la formule géométrique pour k=1 à 5.
  G2 (Convergence de la Densité) : L'erreur relative entre la densité empirique 
       (fenêtre [10^6, 2*10^6]) et la densité théorique D_5 (M_5 = 2310) doit être < 10%.
  G3 (Stabilité du Raffinement) : La variation relative de la densité théorique D_k 
       entre k=4 (M_4=210) et k=5 (M_5=2310) doit être < 5%.
"""
import sys
import time
import numpy as np
import math

def crible_premiers(n):
    """Crible d'Ératosthène optimisé."""
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i :: i] = False
    return np.nonzero(isp)[0]

def nu_H(pattern, p):
    """Nombre de résidus distincts occupés par le motif modulo p."""
    return len(set(h % p for h in pattern))

def geometric_V_k(pattern, primes_k):
    """Calcule V_k via la formule géométrique ∏ (p - ν_H(p))."""
    v_k = 1
    for p in primes_k:
        v_k *= (p - nu_H(pattern, p))
    return v_k

def brute_V_k(pattern, M_k):
    """Calcule V_k par énumération brute des résidus valides modulo M_k."""
    valid_count = 0
    for x in range(M_k):
        is_valid = True
        for h in pattern:
            if (x + h) % M_k == 0: # Simplification: on teste la divisibilité par les facteurs de M_k
                # En réalité, il faut tester si (x+h) est divisible par l'un des p divisant M_k
                pass
        # Méthode plus robuste et directe :
        is_valid = True
        for p in primes_k: # primes_k doit être accessible, on va le passer en argument
            pass
    return valid_count

# Réécriture de brute_V_k pour être autonome et correcte
def brute_V_k_correct(pattern, M_k, prime_factors):
    valid_count = 0
    for x in range(M_k):
        is_valid = True
        for p in prime_factors:
            # Si pour un p, deux éléments du motif ou un élément et 0 tombent sur le même résidu
            # C'est-à-dire si (x + h) % p == 0 pour un h dans pattern
            if any((x + h) % p == 0 for h in pattern):
                is_valid = False
                break
        if is_valid:
            valid_count += 1
    return valid_count

def main():
    t0 = time.time()
    print("=" * 78)
    print("polytope_constellations.py v1 — Cahier II, Porte A : géométrie des polytopes")
    print("=" * 78)
    
    patterns = {
        "jumeaux": [0, 2],
        "triplet_A": [0, 2, 6]
    }
    
    # Les 5 premiers nombres premiers pour M_1 à M_5
    primes_k = [2, 3, 5, 7, 11]
    M_values = [math.prod(primes_k[:k]) for k in range(1, 6)] # [2, 6, 30, 210, 2310]
    
    print(f"\n[1/4] Vérification de l'exactitude combinatoire (G1)...")
    g1_pass = True
    for name, pat in patterns.items():
        print(f"\n  Motif: {name} {pat}")
        for k in range(1, 6):
            M_k = M_values[k-1]
            v_geo = geometric_V_k(pat, primes_k[:k])
            v_brute = brute_V_k_correct(pat, M_k, primes_k[:k])
            match = (v_geo == v_brute)
            if not match:
                g1_pass = False
            status = "OK" if match else "ECHEC"
            print(f"    k={k} (M={M_k:<4}): V_geo={v_geo:<4}, V_brute={v_brute:<4} -> {status}")

    print(f"\n[2/4] Crible et mesure de la densité empirique dans [10^6, 2*10^6]...")
    P_MAX = 2_000_000
    primes = crible_premiers(P_MAX)
    prime_set = set(primes)
    
    window_start = 1_000_000
    window_end = 2_000_000
    window_primes = [p for p in primes if window_start <= p <= window_end]
    window_size = window_end - window_start
    
    print(f"      Fenêtre : {window_size} entiers, contenant {len(window_primes)} premiers.")

    print(f"\n[3/4] Évaluation des gardes G2 (Convergence) et G3 (Stabilité)...")
    g2_pass = True
    g3_pass = True
    
    for name, pat in patterns.items():
        print(f"\n  Motif: {name} {pat}")
        
        # Comptage empirique
        emp_count = 0
        for p in window_primes:
            if p > window_end - max(pat):
                break
            if all((p + h) in prime_set for h in pat):
                emp_count += 1
                
        emp_density = emp_count / window_size
        
        # Densités théoriques D_k = V_k / M_k
        D_k_values = []
        for k in range(1, 6):
            v_geo = geometric_V_k(pat, primes_k[:k])
            M_k = M_values[k-1]
            D_k = v_geo / M_k
            D_k_values.append(D_k)
            
        D_5 = D_k_values[-1]
        
        # G2 : Erreur relative par rapport à D_5
        if D_5 > 0:
            err_rel = abs(emp_density - D_5) / D_5
        else:
            err_rel = 1.0
        g2_match = (err_rel < 0.10)
        if not g2_match: g2_pass = False
        print(f"    Densité empirique : {emp_density:.6f}")
        print(f"    Densité théorique D_5 (M=2310) : {D_5:.6f}")
        print(f"    G2 : Erreur relative = {err_rel:.1%} (Seuil < 10%) -> {'OK' if g2_match else 'ECHEC'}")
        
        # G3 : Stabilité du raffinement entre k=4 et k=5
        D_4 = D_k_values[-2]
        if D_4 > 0:
            var_refinement = abs(D_5 - D_4) / D_4
        else:
            var_refinement = 1.0
        g3_match = (var_refinement < 0.05)
        if not g3_match: g3_pass = False
        print(f"    G3 : Variation D_4 vers D_5 = {var_refinement:.1%} (Seuil < 5%) -> {'OK' if g3_match else 'ECHEC'}")

    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Cahier II, Porte A) :")
    print(f"  G1 (Exactitude Combinatoire)            : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Convergence de la Densité, err <10%): {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Stabilité du Raffinement, var < 5%) : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : La géométrie des polytopes croisés est validée.")
        print("    Le comptage des classes d'ancrage est exactement multiplicatif,")
        print("    et la densité théorique converge rapidement vers l'observation empirique.")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique géométrique ne correspond pas au modèle")
        print("    d'invariance fractale attendu (écart combinatoire ou divergence de densité).")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()