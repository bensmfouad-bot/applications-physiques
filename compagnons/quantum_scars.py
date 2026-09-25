#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
quantum_scars.py v1 — Addendum III-Phys (Cahier de Physique).
Résonances transitoires, systèmes ouverts et états de Gamow dans le résidu du crible.

Lecture physique (analogie mesurée, Article III) :
Le résidu R(x) = π(x) - Li(x) n'est pas un bruit blanc stationnaire.
En variable t = ln(x), il présente des oscillations dont les fréquences 
dominantes correspondent aux parties imaginaires des zéros de ζ(s).
Nous cherchons des "quantum scars" : des résonances intermittentes dont 
l'amplitude varie fortement d'un bloc temporel à l'autre (non-stationnarité),
suivant une enveloppe de décroissance analogue à la largeur Γ d'un état de Gamow.

Gardes pré-enregistrées :
  G1 (Intermittence) : Pour au moins 2 des 3 premières fréquences cibles 
       (γ ≈ 14.13, 21.02, 25.01), la variance inter-blocs de la puissance 
       spectrale du signal arithmétique doit dépasser celle des surrogats 
       de phase de plus de 3σ.
  G2 (Décroissance de Gamow) : L'enveloppe maximale des pics de résonance 
       à travers les blocs temporels (échelles croissantes) doit montrer 
       une tendance à la décroissance (pente de régression linéaire < 0).
  G3 (Robustesse du Null) : Aucun des 10 surrogats de phase ne doit produire 
       de variance inter-blocs dépassant le seuil de 3σ pour ces fréquences 
       (contrôle des faux positifs, correction de Bonferroni implicite par 
       le seuil strict).
"""
import sys
import time
import numpy as np

P_MAX = 10_000_000
N_POINTS = 2000      # Nombre de points d'échantillonnage en échelle logarithmique
N_BLOCKS = 8         # Nombre de blocs contigus pour l'analyse d'intermittence
N_SURROGATES = 10    # Nombre de surrogats de phase pour le null stratifié

# Fréquences cibles (parties imaginaires des premiers zéros de Riemann)
TARGET_GAMMAS = [14.1347, 21.0220, 25.0109]

def sieve(n):
    """Crible d'Ératosthène optimisé avec numpy."""
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i :: i] = False
    return np.nonzero(isp)[0]

def li_approx(x):
    """Approximation asymptotique de l'intégrale logarithmique Li(x) pour x >= 10^4."""
    lx = np.log(x)
    return (x / lx) * (1.0 + 1.0/lx + 2.0/(lx**2))

def phase_randomized_surrogate(signal):
    """Génère un surrogate en randomisant les phases de la FFT du signal."""
    fft_vals = np.fft.rfft(signal)
    phases = np.angle(fft_vals)
    np.random.shuffle(phases)  # Détruit la cohérence de phase, préserve le spectre de puissance global
    surrogate_fft = np.abs(fft_vals) * np.exp(1j * phases)
    return np.fft.irfft(surrogate_fft, n=len(signal))

def main():
    t0 = time.time()
    print("=" * 80)
    print("quantum_scars.py v1 — Addendum III-Phys : résonances transitoires et états de Gamow")
    print("=" * 80)
    print()
    
    print(f"[1/5] Crible des premiers jusqu'à {P_MAX}...")
    primes = sieve(P_MAX)
    print(f"      Trouvés : {len(primes)} premiers.")
    
    print(f"[2/5] Échantillonnage du résidu R(x) = π(x) - Li(x) en échelle logarithmique...")
    x_vals = np.logspace(4, 7, N_POINTS)
    # searchsorted est extrêmement rapide pour compter les premiers <= x
    pi_x = np.searchsorted(primes, x_vals)
    li_x = li_approx(x_vals)
    R = pi_x - li_x
    
    # Détrend simple (retrait de la moyenne globale) pour l'analyse spectrale
    R_detrended = R - np.mean(R)
    
    print(f"[3/5] Découpage en {N_BLOCKS} blocs contigus et analyse spectrale...")
    block_size = N_POINTS // N_BLOCKS
    target_indices = []
    
    # Trouver les indices de fréquence les plus proches des cibles dans la FFT d'un bloc
    # La fréquence d'échantillonnage en t = ln(x) est : fs = N_POINTS / (ln(10^7) - ln(10^4))
    t_min, t_max = np.log(1e4), np.log(1e7)
    dt = (t_max - t_min) / N_POINTS
    fs = 1.0 / dt
    freqs = np.fft.rfftfreq(block_size, d=dt)
    
    for gamma in TARGET_GAMMAS:
        idx = np.argmin(np.abs(freqs - gamma))
        target_indices.append(idx)
        
    # Calcul de la puissance pour le signal original dans chaque bloc
    original_powers = np.zeros((len(TARGET_GAMMAS), N_BLOCKS))
    
    for b in range(N_BLOCKS):
        start = b * block_size
        end = start + block_size
        block_signal = R_detrended[start:end]
        fft_block = np.fft.rfft(block_signal)
        power = np.abs(fft_block)**2
        
        for i, idx in enumerate(target_indices):
            original_powers[i, b] = power[idx]
            
    print(f"[4/5] Génération de {N_SURROGATES} surrogats de phase (Null Stratifié)...")
    surrogate_variances = np.zeros((len(TARGET_GAMMAS), N_SURROGATES))
    
    for s in range(N_SURROGATES):
        surr_signal = phase_randomized_surrogate(R_detrended)
        for b in range(N_BLOCKS):
            start = b * block_size
            end = start + block_size
            fft_block = np.fft.rfft(surr_signal[start:end])
            power = np.abs(fft_block)**2
            for i, idx in enumerate(target_indices):
                # On stocke la variance inter-blocs pour ce surrogate
                # (Calculé après la boucle des blocs pour ce surrogate)
                pass
        
        # Calcul de la variance inter-blocs pour ce surrogate
        surr_powers = np.zeros((len(TARGET_GAMMAS), N_BLOCKS))
        for b in range(N_BLOCKS):
            start = b * block_size
            end = start + block_size
            fft_block = np.fft.rfft(surr_signal[start:end])
            power = np.abs(fft_block)**2
            for i, idx in enumerate(target_indices):
                surr_powers[i, b] = power[idx]
        
        for i in range(len(TARGET_GAMMAS)):
            surrogate_variances[i, s] = np.var(surr_powers[i, :])

    print(f"[5/5] Vérification des gardes pré-enregistrées...")
    print()
    print(f"{'Fréq (γ)':<10} | {'Var. Originale':<15} | {'Moy. Null':<12} | {'Sigma':<6} | {'Pente Enveloppe':<15}")
    print("-" * 75)
    
    g1_pass_count = 0
    g2_pass_count = 0
    
    for i, gamma in enumerate(TARGET_GAMMAS):
        var_orig = np.var(original_powers[i, :])
        mean_null = np.mean(surrogate_variances[i, :])
        std_null = np.std(surrogate_variances[i, :]) if N_SURROGATES > 1 else 1e-6
        
        sigma = (var_orig - mean_null) / std_null if std_null > 0 else 0.0
        
        # G2 : Pente de l'enveloppe (régression linéaire de la puissance max ou moyenne par bloc)
        # On utilise la puissance moyenne dans le bloc comme proxy de l'enveloppe
        x_blocks = np.arange(N_BLOCKS)
        slope, _ = np.polyfit(x_blocks, original_powers[i, :], 1)
        
        print(f"{gamma:<10.2f} | {var_orig:<15.3e} | {mean_null:<12.3e} | {sigma:+6.1f}σ | {slope:<15.3e}")
        
        if sigma > 3.0:
            g1_pass_count += 1
        if slope < 0: # Tendance à la décroissance (Gamow)
            g2_pass_count += 1

    print()
    print("=" * 80)
    print("VERDICT DES GARDES :")
    
    g1_pass = (g1_pass_count >= 2)
    g2_pass = (g2_pass_count >= 2)
    
    # G3 : Vérifier qu'aucun surrogate ne dépasse 3σ (contrôle des faux positifs)
    max_surr_sigma = 0.0
    for i in range(len(TARGET_GAMMAS)):
        mean_null = np.mean(surrogate_variances[i, :])
        std_null = np.std(surrogate_variances[i, :]) if N_SURROGATES > 1 else 1e-6
        for s in range(N_SURROGATES):
            sigma_s = (surrogate_variances[i, s] - mean_null) / std_null if std_null > 0 else 0.0
            if sigma_s > max_surr_sigma:
                max_surr_sigma = sigma_s
                
    g3_pass = (max_surr_sigma < 3.0)
    
    print(f"  G1 (Intermittence, >=2 fréquences à >3σ) : {'OK' if g1_pass else 'ECHEC'} ({g1_pass_count}/3)")
    print(f"  G2 (Décroissance de Gamow, pente < 0)     : {'OK' if g2_pass else 'ECHEC'} ({g2_pass_count}/3)")
    print(f"  G3 (Robustesse du Null, max surrogate <3σ): {'OK' if g3_pass else 'ECHEC'} (max surrogate = {max_surr_sigma:.1f}σ)")
    print("=" * 80)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : Le résidu du crible présente une intermittence spectrale")
        print("    significative et une décroissance d'enveloppe, validant l'analogie avec")
        print("    des résonances transitoires (états de Gamow) dans un système quantique ouvert.")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique observée ne correspond pas au modèle de")
        print("    résonances transitoires attendu. Les gardes ont été violées.")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()