# Entropy Walls and Sphaleron Windows
Baryogenesis in modulus-dominated universes and the Large Volume Scenario (LVS).

Light moduli dominate the early universe and release entropy when they decay. This repository contains the manuscript and the code behind its numbers: the sphaleron-equilibrium window in an early matter-dominated era, the conversion efficiencies of baryon- and lepton-violating decays, the kinetic-energy bound on axion-driven mechanisms, the volume bounds in explicit Swiss-cheese and K3-fibred vacua, the Affleck–Dine analysis with the operator scale M_* as an input, a cross-model comparison, and an explicit amplitude calculation of the CP asymmetry in higgsino decay.

## Layout
```
paper/                     entropy_walls_revised.tex / .pdf, figures/
code/
  make_figures.py          fig_eff, fig_eps, fig_admap
  lvs_decay_widths.py      explicit N=1 supergravity vacua and modulus decay widths
  lvs_decay_widths_scan.py scan over volume (Gamma_aa, Gamma_gg, gaugino mass, Higgs coupling)
  ad_scan_Mstar.py         Affleck–Dine maximal volume vs operator scale M_* = zeta M_Pl
  cross_model_table.py     entropy-wall quantities for KKLT / heterotic / G2 / sequestered LVS
  higgsino_cp.py           Kane–Winkler rates generalised to the stop propagator matrix
  amplitude_check/         explicit two-component-spinor amplitude calculation (see below)
```

## Reproducing the paper's numbers
| Paper item | Script |
|---|---|
| Fig. 2 (efficiencies), Fig. 3 (required epsilon), Fig. 5 (Affleck–Dine plane) | `code/make_figures.py` |
| Table 3 (branching ratios in the vacua) | `code/lvs_decay_widths_scan.py` (mpmath, ~minutes) |
| Table 6 and the zeta scan (Affleck–Dine) | `code/ad_scan_Mstar.py` |
| Table 4 (other moduli stabilisation schemes), m_Phi threshold | `code/cross_model_table.py` |
| Appendix D (explicit amplitude check), Fig. 4 (higgsino) | `code/amplitude_check/` : `test_kw.py`, `test_cp.py`, `test_cp2.py`, `hs_controls.py`, `hs_mixed.py`, `hs_fig.py` |

Requirements: Python 3, numpy, scipy, matplotlib, mpmath. The paper compiles with `pdflatex` (twice).

## Status of the results (read before citing)
* **Re-derived and checked:** the steady-state expansion rate and the dilution coefficient C = 0.16 g_R/g_*; the susceptibility 12/79 and the ratio 28/79; the upper edge of the sphaleron window and the 2.6 MeV closure; the axion bound F <~ 4.2e10 GeV; the volume bounds; the Affleck–Dine coefficient 0.13.
* **Explicit supergravity:** Gamma_aa = m^3/(48 pi M_Pl^2) is recovered to four digits in every vacuum. The vacua are AdS (no uplift), W0 = A = 1, g_s = 0.3.
* **Affleck–Dine scan:** reproduces the paper's zeta = 1 values to within 6%; the dependence on M_* enters only through |phi_0|^2 ~ zeta^{3/2} (the charge efficiency epsilon does not depend on zeta).
* **Higgsino CP asymmetry:** the explicit engine reproduces the tree rate (39) of Kane & Winkler (arXiv:1909.04705) exactly and the f_2 term of their Eq. (41) for right-handed couplings, including the x-dependence. It **differs** from their Eq. (41) in the relative sign of the left-handed term and in the absence of the f_1 term, which implies a vanishing asymmetry for a pure higgsino. This is a claim against a published formula and **has not yet been checked independently** (e.g. with FeynCalc or by the authors). The relative sign of the f and fbar cut contributions is a fermion-ordering sign calibrated against Kane–Winkler for right-handed couplings; it multiplies both chirality structures and does not affect the cancellation.
* **Not regenerated here:** Fig. 1 (entropy dilution) is from the numerical solution of the Boltzmann system of the original code. The ten explicit Swiss-cheese vacua (Table 2 of the first manuscript), used as inputs by `make_figures.py`, `ad_scan_Mstar.py` and `cross_model_table.py`, are taken from the original scripts; `lvs_decay_widths.py` is an independent re-implementation, validated against Gamma_aa.
* **To verify before submission:** author lists and journal references of arXiv:2203.05010, arXiv:1207.2771, arXiv:1312.0954, arXiv:1909.04705; the sequestered-LVS soft-term scaling quoted from arXiv:0906.3297.

## License
MIT (see `LICENSE`).
