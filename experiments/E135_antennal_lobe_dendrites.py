#!/usr/bin/env python3
"""E135, the antennal lobe held to the recorded PN dendrites: the step from ORNs to PNs read glomerulus by glomerulus,
the measured step and then global parameters and terms per glomerulus fitted to the dendrites and read on held-out
flies and odors, and the recorded dendrites put in the model's PN stage.

STUB (2026-11-24): the interface only; bodies raise NotImplementedError except constants and trivial helpers.

The entry is item E135 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04, revised on moving after Tom's words of 09:38, about
09:41 and 09:55, Yang's OK of 09:43 and hal's notes on Zocchi et al. 2022). It builds on E134's shared code and does not
copy it: viki's data_prep/stage_population_data.py and its committed lock (experiments/E134_lock.json), hal's
kcq/population_fit.py (the views and one-fly units), atlas's kcq/stage_readings.py (the model read as each stage is
recorded) and experiments/E134_population_baseline.py (lanes, claims, the lock guard), viki's experiments/E134_flies.py
(the flies' recorded steps and the no-step yardstick).

The readings, fixed before any run (from the entry's Definitions):
- The shared glomeruli are the 28 that the ORN recordings, the PN-dendrite recordings and the model's input all hold
  (E90's `shared_glomeruli`, checked against the lock's populations).
- The ORN input of every arm is the model's own megamat input: the lab's recordings at the ORN axon terminals (olfdb
  dataset 1, 9 flies), converted to rates (`megamat_orn_deltas`). It already holds the lateral signal among ORNs. No arm
  uses Hallem's spiking data. On the APL16 input, paa, euc, -aPine and B-myr come from the validation2 terminal
  recordings.
- The stage's curve is one monotone curve from PN rate to dendrite signal, for all glomeruli and odors (stage_readings'
  `fit_stage_curve`, set on pooled values only).
- Imposed: on the shared glomeruli each PN's odor-evoked rate is set from the recorded dendrites, rank for rank on the
  model's own pooled distribution of PN rates. The arrangement over glomeruli and odors is the flies'; the distribution
  of rates and the time course are the model's. The other glomeruli keep the model's PN stage.
- The global parameters are pn_inhsc, pn_inhadd, ln_inhsc, ln_inhadd and the scale of the PN's saturating curve
  (olfsysm's pn.tanhsc).

The measured step (step 0, recordings only): for each shared glomerulus, the slope of its dendrite response on its own
ORN-terminal response across odors, as a multiple of the median glomerulus's (its measured gain); how closely the
dendrites follow the terminals (Pearson r over odors); the residual by odor; and the other glomeruli whose terminal
responses that residual follows best (its candidate partners). Intervals from resampled flies (both stages' flies,
independently). Each is given on all flies and on both halves of the PN flies, with their agreement. DM1, DL5 and VA2
are shown first. paa and B-myr have no dendrite recording (the dendrites hold the megamat panel only), so of the three
weak-input odors only 2-but is read.

The glomerulus terms: a gain is one factor per shared glomerulus on the drive its PN reads from its own ORN terminals,
between 0.25 and 4; a lateral term adds the terminal response (the ORN delta) of one other glomerulus to that drive,
with a weight between -1 and 1, at most two per glomerulus. olfsysm's LN layer reads only the mean ORN rate over
glomeruli and its PN layer reads its own glomerulus's ORN delta, so the terms act on the PN's drive only and the
normalization keeps the recorded ORN input (`PNStage`, a port of olfsysm's `sim_PN_layer` checked against it). What a
term stands for: something between the recorded ORN terminals and the PN dendrites; no term is an identified
interaction.

Held to the dendrites: every parameter of this step is fitted to the recorded PN dendrites and to nothing downstream.
No parameter per odor is fitted, and no other stage gains one. The PN flies are split in two halves and the 17 odors
into folds, fixed in this experiment's lock (experiments/E135_lock.json) before any fit. A fit sees one half of the
flies and the odors of its folds and is read on the other half, on the odors it did not see.

The rungs: R0, the model's stage as it is. R1, the global parameters. R2, R1 with the gains: first set from the
measured gains with no gain fitted, then fitted from there. R3, R2 with lateral terms, tried only among each
glomerulus's two best measured candidates and added one at a time for as long as the held-out error falls. A rung is
kept if its held-out error on the matrix itself is lower than the rung below by more than twice the standard error over
resamplings of the flies; where fitted gains are no better than the measured ones by that margin, the measured gains
are the rung. A lateral term must also gain more than the same term fitted to dendrites whose odor labels are permuted.

The gates (the entry's Design 1 and Status): step 0 reads recordings only and runs on E134's committed lock. No fit
(step 1) and no model run (step 2) before E134's measure gate reads "usable" and its megamat runs are scored
(`require_e134_gate`); a planner then re-reads the thresholds.

Step 2 (up to 73 runs at fixed scales, both response rates re-held, on the family's eight wirings): on the megamat
input R1, the highest kept rung (each odor predicted by the fit that left it out), that rung with its terms' glomerulus
labels permuted, the dendrites imposed, imposed with glomerulus labels permuted, imposed with odor labels permuted; on
the APL16 input R1, the highest kept rung, and the dendrites imposed for its 12 megamat odors; one run repeats an
imposed run by the second route. The rung's 24 runs are dropped if none is kept. Nothing of the family is refitted.

The two routes into the model's PN stage: route 1, olfsysm's option for the PN stage's rates given directly (viki's
Work claim, Yang's OK of 09:43; `route1_pn_rates`, not implemented until its commit is pushed); route 2, the ORN input
solved glomerulus by glomerulus so that the unchanged PN stage returns the target rates to within SOLVE_TOL_HZ (the
imposed rates on the shared glomeruli and the model's own on the others; Tom, 09:38: an inversion fixed by the
recordings, not a search). A kept rung enters as its table of PN rates by the same route, with no rank mapping.

The classes (on the eight wirings, against E134's own runs; thresholds reset to three times the standard error over
wirings in E134's table before any arm runs): "upstream", "upstream to the boutons only", "not upstream", "the global
parameters carry", "the glomerulus terms carry"; no class for a run outside a band on the APL16 input or with more than
1 % of bouton responses above 300 Hz (`classify`).

Concentrations: megamat throughout, at 10^-3 on both sides. The APL16 input's paa, euc and -aPine are the ORN data
within 1 log10 unit of the flies' concentrations; the imposed arm does not touch them, a kept rung's terms are applied
to them, and every APL16 number involves them.

Outputs: <shared>/results/E135_antennal_lobe_dendrites/ (compact: step0_*.csv, lock.json, step1_*.csv, runs_2_*.csv,
claims/, score_*.csv, provenance_*.json); outputs/E135_antennal_lobe_dendrites/ (local, large: step 2's runs).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part step0
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part lock [--n-folds K]
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E135_antennal_lobe_dendrites.py --part step1|step2 [--lane <host>[-<n>]] [--reverse]
    ./venv/bin/python -u experiments/E135_antennal_lobe_dendrites.py --part score|commands
"""
from __future__ import annotations

import argparse
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E135_antennal_lobe_dendrites'
SCRIPT = 'experiments/E135_antennal_lobe_dendrites.py'

ORN_KEY = ('orn', 'megamat', 'wild_type')
PN_KEY = ('pn_dendrite', 'megamat', 'wild_type')
MEGAMAT_ODORS: Tuple[str, ...] = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms',
                                  'B-cit', 'Lin', 'aa', 'benz', 'va')
N_SHARED = 28
NAMED_GLOMERULI: Tuple[str, ...] = ('DM1', 'DL5', 'VA2')
WEAK_INPUT_ODORS: Tuple[str, ...] = ('paa', 'B-myr', '2-but')    # only 2-but has a dendrite recording
# the units the measured step is read in: dF/F as recorded (the first look's) and the lock's trial-noise units
STEP_UNITS: Tuple[str, ...] = ('dff', 'noise')

GLOBAL_PARAMS: Tuple[str, ...] = ('pn_inhsc', 'pn_inhadd', 'ln_inhsc', 'ln_inhadd', 'pn_tanhsc')
OLFSYSM_DEFAULTS: Dict[str, float] = {'pn_inhsc': 368.66, 'pn_inhadd': 31.41, 'ln_inhsc': 500.0, 'ln_inhadd': 200.0}
GAIN_BOUNDS: Tuple[float, float] = (0.25, 4.0)
LATERAL_BOUNDS: Tuple[float, float] = (-1.0, 1.0)
MAX_LATERAL = 2
N_CANDIDATES = 2
KEEP_SE_MULT = 2.0
RUNGS: Tuple[str, ...] = ('R0', 'R1', 'R2', 'R3')
SOLVE_TOL_HZ = 0.5
N_BOOT = 1000
N_PERMUTATIONS = 20
SEED = 0
CI = 0.9

# step 2: the arms (input panel, arm); 8 wirings each, plus one repeat of an imposed run by the second route
Arm = namedtuple('Arm', 'panel arm')
MEGAMAT_ARMS: Tuple[str, ...] = ('R1', 'rung', 'rung_glom_perm', 'imposed', 'imposed_glom_perm', 'imposed_odor_perm')
APL16_ARMS: Tuple[str, ...] = ('R1', 'rung', 'imposed')
RUNG_ARMS: Tuple[str, ...] = ('rung', 'rung_glom_perm')
REPEAT_ARM = Arm('megamat', 'imposed_route2_repeat')

# the classes (the entry's; the thresholds are reset from E134's table before any arm runs)
MIN_WIRINGS = 6
SOMA_DROP_UPSTREAM = 0.15
SOMA_WITHIN_BOUTONS_ONLY = 0.05
BOUTON_STRENGTH_HALVED = 0.5
CARRY_SHARE = 0.5
THRESHOLD_SE_MULT = 3.0
BOUTON_RATE_CAP_HZ = 300.0
BOUTON_RATE_CAP_SHARE = 0.01


# ---------------------------------------------------------------------------
# paths and guards

def out_dir() -> Path:
    """<shared>/results/E135_antennal_lobe_dendrites/ (created)."""
    raise NotImplementedError


def require_e134_lock() -> dict:
    """E134's committed lock (atlas's `E134_population_baseline.require_lock`); every part reads the recordings through
    it."""
    raise NotImplementedError


def require_e134_gate() -> dict:
    """E134's measure gate, read before any fit or model run: refuses (PermissionError, naming what it waits for) until
    E134 has recorded that its measure is "usable" and its megamat runs are scored."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 0: the recorded step (recordings only)

def shared_glomeruli() -> List[str]:
    """The 28 glomeruli the ORN recordings, the PN dendrites and the model's input share."""
    raise NotImplementedError


def recorded_stage(key: Tuple[str, str, str], units: str = 'dff') -> Dict[str, pd.DataFrame]:
    """Per fly, glomeruli x the 17 megamat odors of trial means, in dF/F (`units='dff'`) or in the lock's trial-noise
    units ('noise'), read through the data module against E134's lock."""
    raise NotImplementedError


def fly_stack(tables: Mapping[str, pd.DataFrame], glomeruli: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """flies x glomeruli x odors (NaN where a fly lacks a glomerulus or an odor)."""
    raise NotImplementedError


def measured_step(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str]) -> pd.DataFrame:
    """The measured step on fly means: per glomerulus the slope through the origin of its dendrite response on its own
    terminal response over `odors`, the slope over the median glomerulus's (the measured gain), Pearson r, and the
    regression on its own response and the odor's summed terminal response."""
    raise NotImplementedError


def residuals(orn: np.ndarray, pn: np.ndarray, step: pd.DataFrame, glomeruli: Sequence[str],
              odors: Sequence[str]) -> pd.DataFrame:
    """The residual by odor of each glomerulus's dendrites about its own slope (glomeruli x odors)."""
    raise NotImplementedError


def candidate_partners(resid: pd.DataFrame, orn_mean: pd.DataFrame, n: int = N_CANDIDATES) -> pd.DataFrame:
    """Per glomerulus, the other glomeruli whose terminal responses its residual follows best (rank, partner, r, the
    least-squares weight of the partner's terminal response on the residual)."""
    raise NotImplementedError


def step_intervals(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                   n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """90 % intervals of the measured gain and r from resampling the flies of both stages."""
    raise NotImplementedError


def halves(pn_flies: Sequence[str], seed: int = SEED) -> Tuple[List[str], List[str]]:
    """The two halves of the PN flies (sizes floor(n/2) and the rest; hal's `split_halves` with one split)."""
    raise NotImplementedError


def part_step0(units: Sequence[str] = STEP_UNITS, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Step 0: the measured step on all flies and on both halves with their agreement, the named glomeruli first, the
    glomeruli that do not follow their terminals, the no-step yardstick (E134's), the saved PN stages against the
    dendrites, and what cannot be read; step0_*.csv, step0_summary.json, provenance_step0.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# this experiment's lock: the halves and the folds, before any fit

def lock_contents(n_folds: Optional[int], seed: int = SEED) -> dict:
    """The halves of the PN flies and the odor folds (each odor left out once), with the commits."""
    raise NotImplementedError


def write_lock(n_folds: Optional[int], seed: int = SEED) -> dict:
    """Writes results/.../lock.json and its committed copy experiments/E135_lock.json."""
    raise NotImplementedError


def read_lock() -> dict:
    """The committed lock (PermissionError if absent or uncommitted)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 1: the PN stage alone

@dataclass
class Terms:
    """The glomerulus terms of a rung: gains per shared glomerulus (1 where absent) and lateral terms
    (glomerulus, partner, weight)."""
    gains: Dict[str, float] = field(default_factory=dict)
    laterals: List[Tuple[str, str, float]] = field(default_factory=list)


class PNStage:
    """olfsysm's ORN and LN layers on the model's input (the LN keeps the recorded input), and its PN layer ported to
    Python with the glomerulus terms on each PN's drive; identical to olfsysm's `run_PN_sims` with no terms
    (`check_identity`)."""

    def __init__(self, orn: pd.DataFrame, mp=None, rv=None):
        raise NotImplementedError

    def run(self, params: Optional[Mapping[str, float]] = None, terms: Optional[Terms] = None) -> pd.DataFrame:
        """The PN's odor-evoked rates (odor-window mean minus the settled pre-odor mean), glomeruli x odors (Hz)."""
        raise NotImplementedError

    def check_identity(self, tol: float = 1e-9) -> dict:
        """The port against olfsysm's PN layer, with no terms, at the defaults and at another setting."""
        raise NotImplementedError


def matrix_error(pred: pd.DataFrame, flies: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                 units_from: Optional[np.ndarray] = None) -> float:
    """The matrix view's error in one-fly units of a prediction (after the stage's curve, set on pooled values against
    the flies' pooled values) against the flies' mean."""
    raise NotImplementedError


def fit_rung(stage: PNStage, rung: str, flies: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
             start: Optional[dict] = None, candidates: Optional[pd.DataFrame] = None) -> dict:
    """One rung fitted to the dendrites of `flies` on `odors` (R1 the global parameters; R2 the gains; R3 lateral terms
    added one at a time among the candidates)."""
    raise NotImplementedError


def cross_fit(stage: PNStage, rung: str, lock: dict, flies_by_half: Mapping[str, np.ndarray],
              glomeruli: Sequence[str], **kwargs) -> dict:
    """A rung fitted fold by fold on one half and read on the other: the assembled out-of-fold prediction (each odor by
    the fit that left it out) and its held-out errors with bootstrap draws."""
    raise NotImplementedError


def keep_rung(err_rung: np.ndarray, err_below: np.ndarray, point_rung: float, point_below: float,
              mult: float = KEEP_SE_MULT) -> bool:
    """Kept if the rung's held-out error is lower than the rung below by more than `mult` standard errors of the
    paired difference over resamplings of the flies."""
    raise NotImplementedError


def permute_odor_labels(flies: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """The dendrites with their odor labels permuted (one permutation for all flies and glomeruli)."""
    raise NotImplementedError


def permute_glomerulus_labels(x: pd.DataFrame, glomeruli: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    """A glomeruli x odors table with its shared glomeruli's labels permuted."""
    raise NotImplementedError


def part_step1(n_perm: int = N_PERMUTATIONS) -> dict:
    """Step 1: the rungs fitted in order, the controls, the kept rung; refuses before E134's gate."""
    require_e134_gate()
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the imposed dendrites, both routes

def imposed_rates(model_rates: pd.DataFrame, dendrites: pd.DataFrame, glomeruli: Sequence[str]) -> pd.DataFrame:
    """The imposed table: on the shared glomeruli the dendrites' arrangement mapped rank for rank onto the model's own
    pooled PN rates there; the other glomeruli keep the model's."""
    raise NotImplementedError


def solve_input(stage: PNStage, target: pd.DataFrame, tol: float = SOLVE_TOL_HZ, max_iter: int = 50) -> Tuple[
        pd.DataFrame, dict]:
    """Route 2: the ORN input solved glomerulus by glomerulus so that the unchanged PN stage returns `target` (every
    glomerulus of the input) to within `tol` Hz; returns the input and the solve's record."""
    raise NotImplementedError


def solve_identity(stage: PNStage, tol: float = SOLVE_TOL_HZ) -> dict:
    """The solve's identity check: the model's own rates as the target return the input it was given."""
    raise NotImplementedError


def route1_pn_rates(*args, **kwargs):
    """Route 1, olfsysm's option for the PN stage's rates given directly (viki's Work claim): not implemented until its
    commit is pushed to olfsysm."""
    raise NotImplementedError("route 1 waits on viki's olfsysm option (the PN stage's rates given directly)")


# ---------------------------------------------------------------------------
# step 2: the runs, and the score

def planned_step2(reverse: bool = False, with_rung: bool = True) -> List[Tuple[Arm, str]]:
    """Step 2's points (arm, wiring): 48 on the megamat input, 24 on the APL16 input, one repeat (73); without the
    rung's 24 if none is kept."""
    raise NotImplementedError


def part_step2(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 2, resumable, one point claimed at a time; refuses before E134's gate."""
    require_e134_gate()
    raise NotImplementedError


def class_thresholds() -> dict:
    """The classes' thresholds reset to THRESHOLD_SE_MULT standard errors over wirings in E134's table."""
    raise NotImplementedError


def classify(scores: pd.DataFrame, thresholds: Mapping[str, float]) -> dict:
    """The entry's classes from step 2's scores against E134's runs."""
    raise NotImplementedError


def part_score() -> dict:
    """Step 2's runs scored and classified; refuses before E134's gate."""
    require_e134_gate()
    raise NotImplementedError


def print_commands() -> None:
    """The lanes' command lines."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'lock', 'step1', 'step2', 'score', 'commands'))
    ap.add_argument('--n-folds', type=int, default=None, help='lock: the number of odor folds')
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--lane', default=None)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    a = ap.parse_args(argv)
    if a.part == 'step0':
        part_step0(n_boot=a.n_boot)
    elif a.part == 'lock':
        write_lock(a.n_folds)
    elif a.part == 'step1':
        part_step1()
    elif a.part == 'step2':
        part_step2(a.lane, a.reverse, a.reclaim)
    elif a.part == 'score':
        part_score()
    elif a.part == 'commands':
        print_commands()


if __name__ == '__main__':
    main()
