#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
quantum_scars.py v1 — Addendum III-Phys (Cahier de Physique).
Cicatrices quantiques et résonances transitoires dans le résidu spectral.

Lecture physique (analogie mesurée, Article III) :
Le résidu du crible ne forme pas un "cristal" spectral rigide, mais présente
des résonances transitoires non-stationnaires, analogues aux états de Gamow
(systèmes ouverts à durée de vie finie) ou aux cicatrices quantiques (quantum scars)
le long d'orbites instables dans les systèmes chaotiques.

Gardes pré-enregistrées :
  G1 (Non-stationnarité) : La variance de l'amplitude des pics spectraux détectés 
       à travers les fenêtres glissantes doit être élevée (le pic est localisé 
       dans l'échelle, pas global).
  G2 (Profil de Résonance) : Le pic spectral le plus significatif doit s'ajuster 
       à un profil de Lorentz (signature d'une résonance à durée de vie finie) 
       avec un coefficient de détermination R² > 0.85.
  G3 (Absence de Cristal) : La distribution des fréquences des pics significatifs 
       ne doit pas former un réseau périodique rigide (la variance des espacements 
       entre pics doit être élevée, rejetant l'hypothèse d'un peigne de Dirac).
"""
import sys
import time
import numpy as np
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq

def crible_premiers(n):
    """Crible d'Ératosthène optimisé."""
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
    print("quantum_scars.py v1 — Addendum III-Phys : résonances transitoires et cicatrices")
    print("=" * 78)
    
    # Paramètres de la fenêtre glissante
    P_MIN = 1_000_000
    P_MAX = 2_000_000
    W = 50_000       # Taille de la fenêtre
    STEP = 20_000    # Pas de glissement
    N_BINS = 100     # Résolution spectrale par fenêtre
    N_SURROGATES = 15 # Nombre de surrogats de phase pour le null stratifié

    print(f"\n[1/4] Crible des premiers jusqu'à {P_MAX}...")
    primes = crible_premiers(P_MAX)
    primes = primes[primes >= P_MIN]
    print(f"      Trouvés : {len(primes)} premiers dans la fenêtre d'étude.")

    print(f"[2/4] Analyse spectrale glissante et calibration par surrogats de phase...")
    peak_amplitudes = []
    significant_peaks = [] # Stocke (fenetre_index, freq, amplitude, sigma)

    # Création des axes de fréquence
    freqs = fftfreq(N_BINS, d=1.0)[:N_BINS//2]
    
    for i, start_p in enumerate(range(P_MIN, P_MAX - W, STEP)):
        end_p = start_p + W
        # Extraire les premiers dans la fenêtre
        window_primes = primes[(primes >= start_p) & (primes < end_p)]
        
        if len(window_primes) < 100:
            continue
            
        # Discrétisation de la densité locale (signal)
        hist, _ = np.histogram(window_primes, bins=N_BINS, range=(start_p, end_p))
        expected = len(window_primes) / N_BINS
        signal = (hist - expected) / np.sqrt(expected) # Fluctuation normalisée
        
        # FFT du signal
        spectrum = np.abs(fft(signal))[:N_BINS//2]
        
        # On cherche le pic dans une bande de fréquence intermédiaire (évite le DC et le bruit haute freq)
        mask = (freqs > 0.05) & (freqs < 0.45)
        valid_spectrum = spectrum[mask]
        valid_freqs = freqs[mask]
        
        if len(valid_spectrum) == 0:
            continue
            
        max_amp = np.max(valid_spectrum)
        peak_freq = valid_freqs[np.argmax(valid_spectrum)]
        peak_amplitudes.append(max_amp)
        
        # Calibration par surrogats de phase (Null Stratifié)
        surrogate_max_amps = []
        for _ in range(N_SURROGATES):
            # Randomisation des phases de la FFT complète
            random_phases = np.exp(1j * 2 * np.pi * np.random.rand(N_BINS))
            surrogate_signal = np.real(np.fft.ifft(np.abs(fft(signal)) * random_phases))
            surrogate_spectrum = np.abs(fft(surrogate_signal))[:N_BINS//2]
            surrogate_max_amps.append(np.max(surrogate_spectrum[mask]))
            
        null_mean = np.mean(surrogate_max_amps)
        null_std = np.std(surrogate_max_amps) if len(surrogate_max_amps) > 1 else 0.1
        sigma = (max_amp - null_mean) / null_std
        
        if sigma > 3.0:
            significant_peaks.append((i, peak_freq, max_amp, sigma))

    print(f"[3/4] Ajustement du profil de résonance (Lorentzien) sur le pic le plus significatif...")
    g1_pass = False
    g2_pass = False
    g3_pass = False
    
    # G1 : Non-stationnarité (Variance des amplitudes)
    if len(peak_amplitudes) > 5:
        var_amp = np.var(peak_amplitudes)
        mean_amp = np.mean(peak_amplitudes)
        # Si le coefficient de variation est > 0.15, c'est non-stationnaire
        cv = np.sqrt(var_amp) / mean_amp if mean_amp > 0 else 0
        g1_pass = (cv > 0.15)
        print(f"      Coefficient de variation des amplitudes : {cv:.3f} (Seuil > 0.15)")

    # G2 : Profil de Lorentz
    if len(significant_peaks) > 0:
        # Prendre le pic le plus significatif
        best_idx, best_freq, best_amp, best_sigma = max(significant_peaks, key=lambda x: x[3])
        print(f"      Pic le plus significatif : Fenêtre {best_idx}, f={best_freq:.3f}, σ={best_sigma:.1f}")
        
        # Extraire un petit voisinage autour du pic pour l'ajustement
        start_p = P_MIN + best_idx * STEP
        end_p = start_p + W
        window_primes = primes[(primes >= start_p) & (primes < end_p)]
        hist, bins = np.histogram(window_primes, bins=N_BINS, range=(start_p, end_p))
        expected = len(window_primes) / N_BINS
        signal = (hist - expected) / np.sqrt(expected)
        spectrum = np.abs(fft(signal))[:N_BINS//2]
        
        # Trouver l'indice du pic dans le spectre complet
        peak_bin_idx = np.argmax(spectrum[(freqs > 0.05) & (freqs < 0.45)]) + np.sum(freqs <= 0.05)
        
        # Fenêtre d'ajustement de ±5 bins autour du pic
        fit_start = max(1, peak_bin_idx - 5)
        fit_end = min(len(freqs)-1, peak_bin_idx + 6)
        
        x_data = freqs[fit_start:fit_end]
        y_data = spectrum[fit_start:fit_end]
        
        try:
            # Ajustement Lorentzien
            initial_guess = [np.max(y_data), best_freq, 0.1, np.min(y_data)]
            popt, pcov = curve_fit(lorentzian, x_data, y_data, p0=initial_guess, maxfev=5000)
            
            # Calcul du R²
            y_fit = lorentzian(x_data, *popt)
            ss_res = np.sum((y_data - y_fit)**2)
            ss_tot = np.sum((y_data - np.mean(y_data))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            
            print(f"      Ajustement Lorentzien : R² = {r_squared:.3f} (Seuil > 0.85)")
            print(f"      Paramètres : Amplitude={popt[0]:.2f}, Centre={popt[1]:.3f}, Largeur(Γ)={popt[2]:.3f}")
            g2_pass = (r_squared > 0.85)
        except Exception:
            print("      Échec de l'ajustement Lorentzien (données trop bruitées).")
            g2_pass = False
    else:
        print("      Aucun pic significatif (> 3σ) trouvé pour l'ajustement.")
        g2_pass = False

    # G3 : Absence de Cristal (Variance des espacements entre pics significatifs)
    if len(significant_peaks) >= 3:
        sig_freqs = np.sort([p[1] for p in significant_peaks])
        spacings = np.diff(sig_freqs)
        var_spacing = np.var(spacings)
        mean_spacing = np.mean(spacings)
        cv_spacing = np.sqrt(var_spacing) / mean_spacing if mean_spacing > 0 else 0
        
        # Si c'était un cristal, cv_spacing serait proche de 0. On attend une variance élevée.
        g3_pass = (cv_spacing > 0.3)
        print(f"      Variabilité des espacements de pics (CV) : {cv_spacing:.3f} (Seuil > 0.3 pour rejeter le cristal)")
    else:
        print("      Trop peu de pics significatifs pour tester la périodicité.")
        g3_pass = False

    print("\n" + "=" * 78)
    print("VERDICT DES GARDES (Rang III) :")
    print(f"  G1 (Non-stationnarité, CV > 0.15)      : {'OK' if g1_pass else 'ECHEC'}")
    print(f"  G2 (Profil de Résonance, R² > 0.85)    : {'OK' if g2_pass else 'ECHEC'}")
    print(f"  G3 (Absence de Cristal, CV_esp > 0.3)  : {'OK' if g3_pass else 'ECHEC'}")
    print("=" * 78)
    
    if g1_pass and g2_pass and g3_pass:
        print("\n>>> AUDIT VERT : Le résidu spectral présente des résonances transitoires")
        print("    localisées, ajustables par un profil de Lorentz (analogue états de Gamow),")
        print("    sans former de réseau périodique rigide (pas de cristal spectral).")
    else:
        print("\n>>> REFUS PUBLIÉ : La dynamique spectrale observée ne correspond pas au modèle")
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