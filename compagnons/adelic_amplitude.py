# adelic_amplitude.py — Addendum I-Phys
# Mesure de la convergence adélique des volumes locaux vers une constante
# portant la signature des zéros de zeta(s).

# 1. Pour H = {0,2} (jumeaux), calculer beta_p(H) pour p <= P_max = 10^5
# 2. Construire l'amplitude p-adique généralisée :
#    A_p(s) = beta_p(H) * |p|_p^s = beta_p(H) * p^{-s}
#    où |p|_p = 1/p est la norme p-adique standard.
# 3. Pour s complexe sur la ligne critique s = 1/2 + it, t in [0, 50],
#    calculer le produit partiel P(s, P) = prod_{p <= P} A_p(s)
# 4. Mesurer la convergence : |P(s, 10^5) / P(s, 10^4) - 1|
# 5. Chercher les pôles de P(s, P) (les "masses" des états de la corde)
#    et comparer avec les zéros de zeta(s) (14,13 ; 21,02 ; 25,01 ; ...)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) 2026 Fouad Bensmail — Tous droits réservés / All rights reserved.
# Licence / License : CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
# =============================================================================
"""
adelic_amplitude.py v1 — Addendum I-Phys (Cahier de Physique).
Le vide adélique : couplage de Tamagawa comme point fixe infrarouge,
flot de renormalisation de la série singulière, et séparation spectrale
géométrie / matière (les zéros de zêta chantent-ils dans le vide eulérien ?).

Refus publié (v0) : l'amplitude naïve A_p(s) = beta_p(H) * p^{-s} est mal
définie — le facteur prod_p p^{-s} = exp(-s * theta(P)) s'effondre vers 0 pour
Re(s) > 0 ; renormalisée, toute dépendance en s disparaît. La v1 ne mesure
que des objets exacts : produit eulérien, queue O(1/(P ln P)), résidu FFT.

Lecture physique (analogie mesurée, pas théorème — Article III) :
  - S(H) = point fixe infrarouge de la couplance adélique ;
  - C_H /(P ln P) = flot RG (fonction bêta du vide) ;
  - résidu spectral du vide = secteur géométrie ; le résidu du crible
    (Addenda XXXVI-XLVII) = secteur matière. La garde G2 teste la séparation.

Gardes pré-enregistrées :
  G0 (couplage)  : |S_{1e6}(jumeaux) - 1.3203236| < 1e-6 (poignée XIII/XXXI).
  G1 (flot RG)   : variation relative de C_H(P) = (S_ancre - S_P)*P*ln(P)
                   < 15 % sur [1e4, 1e6], pour chacun des six motifs.
  G3 (fonction bêta) : pente de ln|S_ancre - S_P| + ln(ln P) vs ln P
                   dans [-1.15, -0.85] (queue en 1/(P ln P)).
  G2 (séparation) : aucun pic FFT du résidu eulérien > 3x la médiane sur
                   t dans [10, 50] -> vide spectralement SILENCIEUX (prédit).
                   Un pic à < 0,45 d'un zéro de zêta -> STRUCTURE (ouverture).
"""
import sys
import time
import math
import numpy as np

P_ANCRE = 10 ** 7
P_FLOW = 10 ** 6
N_PAL = 300
X_MIN_FLOW = 10 ** 3
ZETAS = [14.134725142, 21.022039639, 25.010857580, 30.424876126,
         32.935061588, 37.586178159, 40.918719012, 43.327073281,
         48.005150881]
MOTIFS = {
    "jumeaux": [0, 2],
    "cousins": [0, 4],
    "sexy": [0, 6],
    "triplet_A": [0, 2, 6],
    "triplet_B": [0, 4, 6],
    "quadruplet": [0, 2, 6, 8],
}


def crible_premiers(n):
    isp = np.ones(n + 1, dtype=bool)
    isp[0:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if isp[i]:
            isp[i * i:: i] = False
    return np.nonzero(isp)[0]


def volumes_locaux(H, primes):
    k = len(H)
    nu = np.array([len({h % p for h in H}) for p in primes], dtype=float)
    return (1.0 - nu / primes) * (1.0 - 1.0 / primes) ** (-k)


def main():
    t0 = time.time()
    print("=" * 78)
    print("adelic_amplitude.py v1 — Addendum I-Phys : le vide adélique")
    print("=" * 78)
    print()
    primes = crible_premiers(P_ANCRE)
    print(f"Premiers <= {P_ANCRE} : {len(primes)}")

    lnP = np.linspace(math.log(X_MIN_FLOW), math.log(P_FLOW), N_PAL)
    P_pal = np.exp(lnP)
    idx_pal = np.searchsorted(primes, P_pal, side="right") - 1

    print()
    print("G0 / G1 / G3 : couplage, flot RG, fonction bêta")
    g0 = False
    g1 = True
    g3 = True
    rows = []
    residus = {}
    for nom, H in MOTIFS.items():
        beta = volumes_locaux(H, primes)
        S_cum = np.cumprod(beta)
        S_ancre = S_cum[-1]
        S_pal = S_cum[idx_pal]
        S_1e6 = S_cum[np.searchsorted(primes, P_FLOW, side="right") - 1]
        if nom == "jumeaux":
            g0 = abs(S_1e6 - 1.3203236) < 1e-6
        C_H = (S_ancre - S_pal) * P_pal * np.log(P_pal)
        m_flow = (P_pal >= 1e4)
        var = (C_H[m_flow].max() - C_H[m_flow].min()) / np.median(C_H[m_flow])
        if var >= 0.15:
            g1 = False
        y = np.log(np.maximum(np.abs(S_ancre - S_pal), 1e-18)) + np.log(lnP)
        pente = np.polyfit(lnP[m_flow], y[m_flow], 1)[0]
        if not (-1.15 <= pente <= -0.85):
            g3 = False
        residus[nom] = C_H / np.median(C_H[m_flow]) - 1.0
        rows.append((nom, S_1e6, np.median(C_H[m_flow]), var, pente))
        print(f"   {nom:11s} S(1e6)={S_1e6:.7f}  C_H={np.median(C_H[m_flow]):8.4f}"
              f"  var={var:6.1%}  pente={pente:+.3f}")
    print(f"G0 couplage jumeaux vs 2C2 : {'OK' if g0 else 'ECHEC'}")
    print(f"G1 flot RG (var < 15 % partout) : {'OK' if g1 else 'ECHEC'}")
    print(f"G3 fonction bêta (pente ~ -1 partout) : {'OK' if g3 else 'ECHEC'}")
    print()

    # ---------------- G2 : spectre du vide eulérien ----------------
    print("G2 séparation géométrie / matière : spectre du résidu eulérien")
    ds = lnP[1] - lnP[0]
    w = np.hanning(N_PAL)
    gam = 2 * np.pi * np.fft.rfftfreq(N_PAL, d=ds) / ds
    ratio_max = 0.0
    t_pic = 0.0
    for nom, r in residus.items():
        sp = np.abs(np.fft.rfft(r * w))
        bande = (gam >= 10) & (gam <= 50)
        med = np.median(sp[bande])
        mx = sp[bande].max()
        if med > 0 and mx / med > ratio_max:
            ratio_max = mx / med
            t_pic = gam[bande][int(np.argmax(sp[bande]))]
        print(f"   {nom:11s} max/mediane = {mx / med:5.2f}"
              + (f"  (pic t={gam[bande][int(np.argmax(sp[bande]))]:.2f})" if mx / med > 3 else ""))
    g2_silence = ratio_max < 3.0
    print(f"G2 vide silencieusement géométrique (max/med < 3) : "
          f"{'OK' if g2_silence else 'STRUCTURE'} (max={ratio_max:.2f} à t={t_pic:.2f})")
    print("Amplitudes du vide aux neuf zéros de zêta (jumeaux) :")
    r = residus["jumeaux"]
    sp = np.abs(np.fft.rfft(r * w))
    for g in ZETAS:
        iB = int(np.argmin(np.abs(gam - g)))
        print(f"   t = {g:9.5f}  amplitude {sp[iB]:8.3f}")
    print()

    if not (g0 and g1 and g3):
        print("VERDICT : ECHEC DE GARDE (vide adélique hors de ses ancres)")
    elif g2_silence:
        print("VERDICT : AUDIT VERT — SEPARATION MESUREE : le vide adélique est")
        print("  spectralement silencieux ; les résonances zêta de l'édifice vivent")
        print("  dans le résidu du crible (matière), pas dans la queue eulérienne")
        print("  (géométrie). Le point fixe infrarouge de la couplance est le")
        print("  nombre de Tamagawa ; son flot RG est en 1/(P ln P).")
    else:
        print("VERDICT : STRUCTURE — le vide eulérien chante à t =", f"{t_pic:.2f}",
              ": ouverture publiée ; la séparation géométrie/matière est refusée.")
    print(f"\nTemps total : {time.time() - t0:.2f} s")
    print("[OK]")


if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()