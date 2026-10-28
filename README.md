# Chemical-Electric Multimode Propulsion for Small Satellites with Nitrous Oxide as a Shared Propellant

Code for O. D. FitzGerald, *Chemical-Electric Multimode Propulsion for Small
Satellites with Nitrous Oxide as a Shared Propellant*, MSc thesis, Delft
University of Technology, 2026. A reduced-order model of an electrodeless (ECR)
thruster run on nitrous oxide, built on the framework of Petro and Sedwick
(2017) and anchored to the xenon operating point of Cannat et al. (2015), and a
spacecraft-level model that sizes chemical, electric and multimode propulsion
systems sharing N2O for a LEO and a GEO mission.

## Installation

Python 3.13.

```bash
pip install -r requirements.txt
```

## Usage

```bash
python run.py check_thruster                   # thruster-model checks (33)
python run.py check_spacecraft                 # spacecraft-model checks (181)
python run.py all                              # thesis tables (CSV) and figures (PDF)
python scripts/ps2017_reproduction.py          # Petro and Sedwick (2017) reproduction figures
python scripts/make_appendix_b.py --help       # Appendix B tables from the two check logs
```

Outputs are written to `output/`. `python run.py --list` lists the tasks.

## Layout

- `multimode/thruster/`: thruster model (cross sections, ionisation cost, exhaust, thrust efficiency)
- `multimode/spacecraft/`: spacecraft model (parameters, missions, architectures, closure, feasibility, sweeps)
- `multimode/analyses/`: the tasks run by `run.py`, including the two check suites
- `multimode/plotting/`: figure style and thesis figures
- `scripts/`: data converters, the Petro and Sedwick reproduction, the response-variant sweep, Appendix B and the CEA cross-check
- `data/`: cross sections and digitised reference data

## Data sources

- Xe and Ar cross sections: Hayashi database, www.lxcat.net, retrieved 20 April 2026; not included (see below).
- N2O cross sections: Song, M.-Y. et al. (2019), J. Phys. Chem. Ref. Data 48, 043104.
- H2O cross sections: Itikawa, Y. and Mason, N. (2005), J. Phys. Chem. Ref. Data 34, 1; electronic excitation from Thorn et al. (2007).
- Reproduction curves: Petro, E. M. and Sedwick, R. J. (2017), J. Propul. Power 33(6), 1410-1417, digitised; detachment cross sections from Fridman (2021), digitised.
- Exit-velocity validation: Sheppard, A. J. and Little, J. M. (2022), J. Propul. Power 38(6), 998-1008, Figs. 7a, 7c and 7d, digitised.
- Thruster anchor: Cannat et al. (2015), Table II.

The Xe and Ar cross sections are not included in this repository (LXCat does
not allow redistribution). Download them from the Hayashi database at
[www.lxcat.net](https://www.lxcat.net) (Data center, Hayashi database; target
species Xe and Ar, all processes, in LXCat text format) and save them as

```
data/lxcat/xe_hayashi.txt
data/lxcat/ar_hayashi.txt
```

## Citation

O. D. FitzGerald, *Chemical-Electric Multimode Propulsion for Small Satellites
with Nitrous Oxide as a Shared Propellant*, MSc thesis, Delft University of
Technology, 2026.

## Licence

MIT; see [LICENSE](LICENSE).
