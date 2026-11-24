#!/usr/bin/env python3
"""E134 step 1: the flies alone. No model is read.

The entry is item 2 of experiment_priorities.md ("the models as populations at every recorded stage"), design step 4
("Step 1, the flies alone"). For every recorded population of `data_prep.stage_population_data` (the lock first):
- each view (strength by odor, sharing on who responds and on how strongly, breadth among units, and for glomeruli the
  matrix itself) with its one-fly spread (the root mean squared distance of one fly's view from the other flies' mean
  view) and its floor (a resampled set of the same flies against the measured set), and the flies' own range (each fly
  against the others, in one-fly units);
- the energy distance between the flies' pooled tunings, as the guard;
- the recorded steps (ORNs to PN dendrites, PN dendrites to PN boutons, PN boutons to KC claws, KC claws to KC somata,
  and the same fly's claws to somata) as changes of each view, with intervals from resampling the flies of both stages;
- the no-step yardstick of every downstream stage: the error of the upstream recordings' view taken as the prediction;
- Remy's megamat KC somata against Yang's control flies on the 12 odors the panels share.

The measures (views, one-fly units, the population error, the energy distance) are hal's kcq/population_fit.py; this
script only calls them.

Outputs: <shared>/results/E134_population_baseline/flies_*.csv and flies_table.md.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_flies.py [--n-boot 1000] [--seed 0]
"""
from __future__ import annotations

import argparse
from typing import Optional, Sequence

import pandas as pd

EXPERIMENT = 'E134_population_baseline'
# the recorded steps, upstream -> downstream, as population keys of data_prep.stage_population_data
STEPS = (
    (('orn', 'megamat', 'wild_type'), ('pn_dendrite', 'megamat', 'wild_type')),
    (('pn_dendrite', 'megamat', 'wild_type'), ('pn_bouton', 'megamat', 'wild_type')),
    (('pn_bouton', 'megamat', 'wild_type'), ('kc_claw', 'megamat', 'wild_type')),
    (('kc_claw', 'megamat', 'wild_type'), ('kc_soma', 'megamat', 'wild_type')),
    (('kc_claw', 'megamat', 'same_fly'), ('kc_soma', 'megamat', 'same_fly')),
    (('orn', 'validation2', 'wild_type'), ('kc_soma', 'validation2', 'wild_type')),
    (('orn', 'mixtures', 'wild_type'), ('kc_soma', 'mixtures', 'wild_type')),
)
N_BOOT = 1000
SEED = 0


def flies_views(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per population and view: the one-fly spread, the floor, and each fly against the others (one row per
    population, view and fly; flies_views.csv)."""
    raise NotImplementedError


def flies_energy(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per population: the energy distance of each fly from the others and its floor, the guard
    (flies_energy.csv)."""
    raise NotImplementedError


def recorded_steps(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per step of `STEPS` and view: the change from the upstream to the downstream flies' mean view, with a 90 %
    interval from resampling the flies of both stages (flies_steps.csv)."""
    raise NotImplementedError


def no_step_yardstick(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per downstream stage of `STEPS` and view: the error, in the downstream stage's one-fly units, of the upstream
    flies' mean view taken as the prediction (flies_no_step.csv)."""
    raise NotImplementedError


def megamat_vs_control(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Remy's megamat KC somata against Yang's control flies on the 12 shared odors, each view in the control flies'
    one-fly units and the reverse (flies_megamat_vs_control.csv)."""
    raise NotImplementedError


def write_table() -> None:
    """flies_table.md from the CSVs."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.parse_args(argv)
    raise NotImplementedError


if __name__ == '__main__':
    main()
