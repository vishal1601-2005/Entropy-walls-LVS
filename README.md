# Entropy Walls and Sphaleron Windows in the Large Volume Scenario

Code, data and LaTeX source for the paper *Entropy Walls and Sphaleron Windows in the Large Volume Scenario*
(Vishal Hiralal Thakur, ORCID [0009-0008-8406-0508](https://orcid.org/0009-0008-8406-0508)).
Every number, table and figure in the paper is produced by the scripts in `code/`.

## What the paper does

Light moduli dominate the early universe before they decay, and the entropy released by the decay dilutes any
earlier baryon asymmetry. The paper asks which baryogenesis mechanisms survive such a moduli-dominated era and
applies the answer to explicit vacua of the Large Volume Scenario (LVS):

1. solves the coupled modulus-radiation-asymmetry Boltzmann system with electroweak sphaleron washout;
2. proves that a rotating field-space basis cannot source a baryon-number chemical potential, and bounds
   axion-driven mechanisms by an effective anomalous decay constant F;
3. computes explicit Swiss-cheese and K3-fibred LVS vacua from the supergravity potential;
4. maps which channels survive as a function of the compactification volume (pre-decay electroweak baryogenesis,
   decay leptogenesis, axion velocity, B-violating decays, Affleck-Dine), including a numerically converged
   Affleck-Dine charge efficiency.

## Layout

    paper/   main.tex (+ macro files numbers*.tex, tables, figures) and main.pdf
    code/    all scripts and the JSON/NPY data they use
    requirements.txt

## Reproducing the results

    pip install -r requirements.txt
    cd code

| Script | Produces | Approx. runtime |
|---|---|---|
| `xi.py` | baryon susceptibility 12/79, check B=(28/79)(B-L) | seconds |
| `emde.py` | EMDE background, dilution D and its coefficient | ~1 min |
| `sph.py`, `washout.py` | sphaleron window, washout efficiencies eta_B, eta_L | minutes |
| `gen.py`, `gen2.py`, `gen3.py` | Sections 2-4 numbers and figures | minutes |
| `scan_lvs.py` (uses `swiss3.py`, `swiss4.py`) | explicit Swiss-cheese vacua (`lvs_explicit.json`) | ~1-2 min |
| `k3geom.py`, `k3explicit.py`, `k3scan.py` | K3-fibred Kahler geometry and vacua | minutes |
| `gen4.py`, `gen5.py` | LVS bounds tables and figures | minutes |
| `ad_sim3.py`, `ad_scan3.py` | converged Affleck-Dine charge efficiency (matter-dominated background) and scan | ~1 min |
| `ad.py`, `ad_lib.py`, `adgen*.py` | Affleck-Dine volume windows | seconds |
| `check_dilution_chain.py` | independent check of the dilution factor used in the AD result | ~1 min |

Scripts read and write in the current directory. To rebuild the paper, copy the regenerated `numbers*.tex`,
table files and `fig_*.pdf` into `paper/` and run `pdflatex main.tex` twice.

## Notes

* The explicit vacua use the potential of Balasubramanian et al. (Swiss cheese) and Cicoli, Burgess, Quevedo
  (K3-fibred, eq. 3.56 and the parameter sets LV, SV1, SV2 of that paper).
* `ad_sim.py`, `ad_sim2.py`, `ad_scan.py`, `ad_scan2.py` are earlier Affleck-Dine simulations kept for transparency:
  the first had unconverged initial conditions, the second used a radiation background. The paper uses `ad_sim3.py`.
* Results depend on stated assumptions (lattice sphaleron-rate fit, approximate g_*(T), parametrised decay width);
  see the Discussion section of the paper.

## License

MIT, see `LICENSE`.
