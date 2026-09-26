#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
quantum_scars_v2.py — Addendum III-Phys (v2).
Cicatrices quantiques et résonances transitoires via la statistique des écarts (gaps).

Lecture physique (analogie mesurée, Article III) :
Contrairement à la v1 (densité locale bruitée), la v2 analyse les écarts 
dépliés s_n = g_n / <g>. Les fluctuations δ_n = s_n - 1 portent la signature 
spectrale fine (loi de Montgomery-Odlyzko, états de Gamow) sans la tendance 
globale du théorème des nombres premiers.

Gardes pré-enregistrées :
  G1 (Non-stationnarité) : CV de l'amplitude des pics spectraux dominants > 0.15.
  G2 (Profil de Résonance) : Ajustement Lorentzien du pic le plus significatif avec R² > 0.85.
  G3 (Absence de Cristal) : CV des espacements entre fréquences de pics significatifs > 0.3.
"""
import sys
import time
import numpy as np
from scipy.optimize import curve_fit
from scipy.fft import fft, ifft

def crible_premiers(n):
    """Crible d'Ératosthène optimisé avec numpy."""
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i :: i] = False
    return np.nonzero(isp)[0]

def lorentzian(x, amp, center, gamma, offset):
    """Profil de Lorentz pour l'ajustement des résonances (états de Gamow)."""
    return amp / (1.0 + ((x - center) / (gamma / 2.0))**2) + offset

def main():
    t0 = time.time()
    print("=" * 78)
    print("quantum_scars_v2.py — Addendum III-Phys (v2) : statistique des écarts dépliés")
    print("=" * 78)
    
    # Paramètres de la fenêtre glissante (en nombre de premiers, pas en valeur P)
    P_MIN_IDX = 78498  # Index du premier ~ 1,000,000
    P_MAX_IDX = 148933 # Index du premier ~ 2,000,000
    N_GAPS = 5000      # Nombre d'écarts par fenêtre (taille fixe pour FFT constante)
    STEP = 2000        # Pas de glissement (en nombre de premiers)
    N_SURROGATES = 15  # Nombre de surrogats de phase pour le null stratifié

    print(f"\n[1/4] Crible des premiers jusqu'à 2,000,000...")
    primes = crible_premiers(2_000_000)
    print(f"      Trouvés : {len(primes)} premiers.")

    print(f"[2/4] Dépliage (unfolding) et analyse spectrale des fluctuations d'écarts...")
    peak_amplitudes = []
    significant_peaks = [] # Stocke (window_idx, freq, amplitude, sigma)

    # Création des axes de fréquence (normalisés par rapport à N_GAPS)
    freqs = np.fft.fftfreq(N_GAPS)[:N_GAPS//2]
    
    # Bande de fréquence d'intérêt (évite le DC et le bruit haute fréquence)
    mask = (freqs > 0.02) & (freqs < 0.30)
    
    for i, start_idx in enumerate(range(P_MIN_IDX, P_MAX_IDX - N_GAPS, STEP)):
        end_idx = start_idx + N_GAPS
        window_primes = primes[start_idx:end_idx+1]
        
        if len(window_primes) < N_GAPS + 1:
            break
            
        # 1. Calcul des écarts bruts
        gaps = np.diff(window_primes)
        
        # 2. Dépliage (Unfolding) : normalisation par la moyenne locale
        mean_gap = np.mean(gaps)
        s_n = gaps / mean_gap
        
        # 3. Fluctuations pures
        delta_n = s_n - 1.0
        
        # 4. FFT des fluctuations
        spectrum = np.abs(fft(delta_n))[:N_GAPS//2]
        valid_spectrum = spectrum[mask]
        valid_freqs = freqs[mask]
        
        if len(valid_spectrum) == 0:
            continue
            
        max_amp = np.max(valid_spectrum)
        peak_freq = valid_freqs[np.argmax(valid_spectrum)]
        peak_amplitudes.append(max_amp)
        
        # 5. Calibration par surrogats de phase (Null Stratifié)
        surrogate_max_amps = []
        fft_full = fft(delta_n)
        for _ in range(N_SURROGATES):
            random_phases = np.exp(1j * 2 * np.pi * np.random.rand(N_GAPS))
            surrogate_signal = np.real(ifft(np.abs(fft_full) * random_phases))
            surrogate_spectrum = np.abs(fft(surrogate_signal))[:N_GAPS//2]
            surrogate_max_amps.append(np.max(surrogate_spectrum[mask]))
            
        null_mean = np.mean(surrogate_max_amps)
        null_std = np.std(surrogate_max_amps) if len(surrogate_max_amps) > 1 else 0.1
        sigma = (max_amp - null_mean) / null_std
        
        if sigma > 3.0:
            significant_peaks.append((i, peak_freq, max_amp, sigma))

    print(f"[3/4] Évaluation des gardes pré-enregistrées...")
    g1_pass = False
    g2_pass = False
    g3_pass = False
    
    # G1 : Non-stationnarité (Variance des amplitudes)
    if len(peak_amplitudes) > 5:
        cv = np.std(peak_amplitudes) / np.mean(peak_amplitudes) if np.mean(peak_amplitudes) > 0 else 0
        g1_pass = (cv > 0.15)
        print(f"      G1 : Coefficient de variation des amplitudes = {cv:.3f} (Seuil > 0.15) -> {'OK' if g1_pass else 'ECHEC'}")

    # G2 : Profil de Lorentz sur le pic le plus significatif
    if len(significant_peaks) > 0:
        best_idx, best_freq, best_amp, best_sigma = max(significant_peaks, key=lambda x: x[3])
        print(f"      Pic le plus significatif : Fenêtre {best_idx}, f={best_freq:.4f}, σ={best_sigma:.1f}")
        
        # Reconstruire le signal de cette fenêtre spécifique pour l'ajustement
        start_idx = P_MIN_IDX + best_idx * STEP
        window_primes = primes[start_idx:start_idx + N_GAPS + 1]
        gaps = np.diff(window_primes)
        delta_n = (gaps / np.mean(gaps)) - 1.0
        spectrum_full = np.abs(fft(delta_n))[:N_GAPS//2]
        
        # Trouver l'indice du pic dans le spectre complet
        peak_bin_idx = np.argmax(spectrum_full[mask]) + np.sum(~mask & (np.arange(len(spectrum_full)) < len(spectrum_full)//2))
        # Simplification : on cherche l'indice global correspondant à best_freq
        peak_bin_idx = np.argmin(np.abs(freqs - best_freq))
        
        # Fenêtre d'ajustement de ±8 bins autour du pic
        fit_start = max(1, peak_bin_idx - 8)
        fit_end = min(len(freqs)-1, peak_bin_idx + 9)
        
        x_data = freqs[fit_start:fit_end]
        y_data = spectrum_full[fit_start:fit_end]
        
        try:
            initial_guess = [np.max(y_data), best_freq, 0.05, np.min(y_data)]
            popt, pcov = curve_fit(lorentzian, x_data, y_data, p0=initial_guess, maxfev=5000)
            
            y_fit = lorentzian(x_data, *popt)
            ss_res = np.sum((y_data - y_fit)**2)
            ss_tot = np.sum((y_data - np.mean(y_data))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            
            print(f"      G2 : Ajustement Lorentzien R² = {r_squared:.3f} (Seuil > 0.85) -> {'OK' if r_squared > 0.85 else 'ECHEC'}")
            print(f"           Paramètres : Amp={popt[0]:.2f}, Centre={popt[1]:.4f}, Largeur(Γ)={popt[2]:.4f}")
            g2_pass = (r_squared > 0.85)
        except Exception as e:
            print(f"      G2 : Échec de l'ajustement Lorentzien ({e}) -> ECHEC")
            g2_pass = False
    else:
        print("      G2 : Aucun pic significatif (> 3σ) trouvé pour l'ajustement -> ECHEC")
        g2_pass = False

    # G3 : Absence de Cristal (Variance des espacements entre pics significatifs)
    if len(significant_peaks) >= 3:
        sig_freqs = np.sort([p[1] for p in significant_peaks])
        spacings = np.diff(sig_freqs)
        cv_spacing = np.std(spacings) / np.mean(spacings) if np.mean(spacings) > 0 else 0
        
        g3_pass = (cv_spacing > 0.3)
        print(f"      G3 : Variabilité des espacements de pics (CV) = {cv_spacing:.3f} (Seuil > 0.3) -> {'OK' if g3_pass else 'ECHEC'}")
    else:
        print("      G3 : Trop peu de pics significatifs pour tester la périodicité -> ECHEC")
        g3_pass = False

    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Rang III, v2) :")
    print(f"  G1 (Non-stationnarité, CV > 0.15)      : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Profil de Résonance, R² > 0.85)    : {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Absence de Cristal, CV_esp > 0.3)  : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : Les fluctuations d'écarts dépliés révèlent des résonances")
        print("    transitoires localisées, ajustables par un profil de Lorentz, sans former")
        print("    de réseau périodique rigide. L'analogie avec les états de Gamow tient.")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique spectrale des écarts ne correspond pas au modèle")
        print("    de résonances transitoires non-stationnaires attendu. Les gardes ont été violées.")
        
    print(f"\nTemps total d'exécution : {time.time() - t0:.2f} s")
    print("[OK]")

if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()