#!/usr/bin/env python3
"""E89: the co-response read with a checked estimate, and the breadth weighting inside the models' APL.

The entry is item E89 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-13 03:20). E88's factors throughout (R recruitment, L co-response, G amplitude factor, O = R * L, s = R * L * G;
`E88_responder_sets`), C2's family (the breadth-2 APL-effect fit of E67) on E72 part B's eight wirings, S = real, P1,
P2, M1 and T = P3, P4, M2, M3, halves A / B of the recordings E77's (`recording_halves()`, seed 0) and C / D the second
split's (seed 1).

Parts (CLI `--part`):
- A (no simulations): A0 the candidate estimates of the co-response on known truth (E88's pseudo-recordings, ten a
  wiring, three trials, seeds SEEDS) at three noise levels: L1 the reading's own noise (E73's calibration), L2 the noise
  raised to the flies' repeatability, independent between odors (E88's reading 1b), and L3 the flies' repeatability with
  a share of each KC's noise common to the odors of a repeat, set per wiring and condition so that the pseudo-recordings'
  same-repeat gap is the flies' (E81's statistic), and at half and one and a half times that share (L3h, L3x). The
  candidates: (a) E88's two estimates of O ('tm' trial means, 'nr' trial noise removed), (b) their geometric mean
  ('gm'), (c) latent classes after deconvolution ('lc'), (d) a correction by each odor's agreement between repeats
  ('hf'); R always from the trial means and s from the noise-removed estimate. The choice; A1 E88's part A read again
  with the chosen candidate (factors, cells, the split, the groups, the base); A2 E88's B1 with it. estimator_check.csv,
  estimator_choice.json, factors.csv, cells.csv, changes.csv, verdict_A.json.
- B (no simulations): the saved families in which only the weighting differs (E66's part C, E67's exponent fits, E12's
  and E19's grids apart) by the main target and E88's archive scan; weighting_archive.csv, verdict_B.json.
- C0 / C1 / C2 / C3 (simulations): twelve weightings of APL's inhibition by a KC's breadth (WEIGHTINGS) on S's four
  wirings, on the fits as they are (C1) and on the base with the flies' recruitment (C2); identities first (C0); the
  choice per base and half; the confirmation on T (C3) where there is something to confirm, with E83's gate and the
  permuted control. Cscore: readings_C.csv, coresponse_C.csv, choice_C.csv, gate_C.csv, verdict_C.json.
- variants: variants.csv and model_pair_distances_long.csv of C's runs, for the archive.
- table: comparison_table.md.

Outputs: results/E89_apl_breadth_weighting/; spike counts and responses under outputs/E89_apl_breadth_weighting/C/.

Usage (from the repository root):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E89_apl_breadth_weighting.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E89_apl_breadth_weighting.py --part C1
"""
from __future__ import annotations

import argparse
import json
import os
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E87_inhibition_form as E87
from experiments import E88_responder_sets as E88

EXPERIMENT = 'E89_apl_breadth_weighting'

# ---------------------------------------------------------------------------
# part A: the settings and the thresholds (the entry's Design; fixed before the flies are read with any candidate)

SEEDS = (E88.SEED, E88.SEED + 1)        # E88's pseudo-recordings' seed and a second beside it
N_PSEUDO, N_TRIALS = E88.N_PSEUDO, E88.N_TRIALS
EST_TOL = E88.EST_TOL                   # 0.03: a factor's group mean off by at most this, in units of the similarity
CANDIDATES = ('tm', 'nr', 'gm', 'lc', 'hf')
CANDIDATE_NAMES = {'tm': 'trial means (E88 reading 1)', 'nr': 'trial noise removed (E88 reading 2)',
                   'gm': 'the geometric mean of the two', 'lc': 'latent classes after deconvolution',
                   'hf': "a correction by each odor's agreement between repeats"}
LEVELS = ('L1', 'L2', 'L3h', 'L3', 'L3x')
SHARE_SCALE = {'L3h': 0.5, 'L3': 1.0, 'L3x': 1.5}
BINDING_LEVELS = ('L1', 'L3h', 'L3', 'L3x')
WHO_FACTORS = ('R', 'L', 's')           # passes for who responds
FULL_FACTORS = ('R', 'L', 's', 'G')     # passes in full
SHARE_BOUNDS = (0.0, 0.95)              # L3: the search range of the shared share of the noise variance
# (c)'s free choice, the classes' edges in z beside the criterion, set on the firing-ceiling population then left alone
LC_EDGES = {'crit': (), 'coarse': (0.25, 0.6, 1.5, 2.5, 4.0), 'fine': (0.0, 0.35, 0.7, 1.3, 1.7, 2.3, 3.2, 4.5)}
LC_STEP, LC_NPMLE_ITERS, LC_EM_ITERS = 0.05, 300, 200
LC_SET_POPULATION = 'fc'
CHECK_POPULATIONS = ('b2', 'base')      # the check is read on these
A_POPULATIONS = ('b2', 'fc', 'S_unchanged', 'base')
# E88's thresholds, unchanged
G_MIN, G_RATIO, GROUP_RATIO, D2_BASE_GAIN, B1_FACTOR = E88.G_MIN, E88.G_RATIO, E88.GROUP_RATIO, E88.D2_BASE_GAIN, E88.B1_FACTOR

# ---------------------------------------------------------------------------
# part B

B_GAIN = 0.003                          # the member nearest the flies' change lower than the exponent-2 member by this

# ---------------------------------------------------------------------------
# part C: the weightings and the thresholds

MEASURES = ('F', 'S', 'N')
MEASURE_NAMES = {'F': "the family's: the KC's spike total without APL in the bouton base",
                 'S': 'the summed response: the KC\'s summed spike count over the odors it responds to (two or more spikes) '
                      'in the run the APL acts on, without APL',
                 'N': 'the number of odors it responds to in that run'}
EXPONENTS = {'F': (0.25, 0.5, 1.0), 'S': (0.25, 0.5, 1.0, 2.0), 'N': (0.25, 0.5, 1.0, 2.0)}
REFERENCE = ('F', 2.0)
NO_WEIGHTING = ('F', 0.0)               # exponent 0: every factor 1 (the measure does not matter)


def weightings() -> List[Tuple[str, float]]:
    """The twelve weightings (measure, exponent), no weighting first, then the reference (F, 2)."""
    return [NO_WEIGHTING] + [(m, k) for m in MEASURES for k in EXPONENTS[m]] + [REFERENCE]


def wtag(w: Tuple[str, float]) -> str:
    """A weighting's tag: 'none', or '<measure>_k<exponent>' ('F_k2' the reference)."""
    m, k = w
    return 'none' if k == 0 else f'{m}_k{k:g}'


ARMS = ('C1', 'C2')                     # the fits as they are; the base with the flies' recruitment
RAW_CRIT = E87.RAW_CRIT                 # a model response is two or more spikes
CHOICE_MIN_GAIN = 0.0015                # below this gain on S against half X there is nothing to confirm
GATE_GAIN, GATE_BORDER = 0.003, (0.0015, 0.0045)
GATE_SINGLES = 6                        # of T's eight wiring-by-direction cells
GATE_RAW_SLACK, GATE_MEAN_SLACK, FLIES_MEAN = 0.005, 0.01, 0.178
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075
CONTROL_RATIO = 0.5
CONTROL_SEEDS = (1, 2)
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
SPLITS = {0: ('A', 'B'), 1: ('C', 'D')}  # E77's recording_halves(seed)
UNCHANGED_T = {'A': 0.0565, 'B': 0.0538}  # T's unchanged fits against each half (E83; the entry's references)
UNCHANGED_S = {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
C2_UNITS = 'noise'                      # the base's per-odor recruitment, as E87's D2 and E88's D1

_CACHE: Dict[tuple, object] = {}


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


# ---------------------------------------------------------------------------
# part A: the candidates as functions of a recording's single trials

def lc_overlap(trials: np.ndarray, crit: float, edges: Sequence[float] = (), se: Optional[np.ndarray] = None
               ) -> Tuple[np.ndarray, np.ndarray]:
    """(c) latent classes after deconvolution: (O, f) of one recording's single trials (repeats, KCs, odors; graded z).
    Per odor the distribution of a KC's mean response with the trial noise deconvolved (NPMLE on a grid, each KC's own
    trial SD pooled over odors), cut into classes at `crit` and at `edges`; per pair the table of class weights by
    maximum likelihood (EM); O (triu order) the share of KCs above `crit` for both over sqrt(f_a f_b), f per odor."""
    raise NotImplementedError


def hf_overlap(binary: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    """(d) a correction by each odor's agreement between repeats: (O, f, params) of one recording's thresholded single
    trials (repeats, KCs, odors; binary). Per odor two classes (responder or not) with a hit rate h and a false-alarm
    rate e, fitted with their share f to the numbers of repeats responding (a binomial mixture, EM); per pair the share
    of KCs in both classes n11 solved from the mean product of the two odors over different repeats, which is
    n11 (h_a - e_a)(h_b - e_b) + f_a h_a e_b + f_b e_a h_b + (1 - f_a - f_b) e_a e_b; O = n11 / sqrt(f_a f_b)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A: the model's readings at the three noise levels

class ModelReading89(E88.ModelReading):
    """E88's `ModelReading` with a seed (SEEDS) and the third noise level: a per-KC, per-repeat term common to the odors
    of a repeat, its share of each KC's noise variance set per condition so that the pseudo-recordings' same-repeat gap
    is the flies'. `level_data(level, c)` gives the trial means, single trials and both criteria of a level;
    `candidate_factors(level, cand)` the population-ready per-pseudo-recording factors."""

    def __init__(self, o, label: str = '', seed: int = E88.SEED, **kw):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# part A: the check, the choice, the readings

def estimator_check(members: Sequence['ModelReading89'], levels: Sequence[str] = LEVELS) -> pd.DataFrame:
    """Each candidate's factors (the population's) minus the truth by group, condition and level: 'err' in units of the
    similarity (E88's convention), 'off' beyond EST_TOL."""
    raise NotImplementedError


def choose_candidate(ec: pd.DataFrame, population: str) -> dict:
    """The thresholds of the entry on one population: per candidate whether it passes for who responds (R, L and s
    inside the tolerance at every binding level, both seeds) and in full (G too); the chosen one (in full first, then
    the smallest largest error of L at L3)."""
    raise NotImplementedError


def part_a(out: Path, n_workers: int = 4) -> dict:
    """A0, the choice, A1 and A2; verdict_A.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B

def part_b(out: Path) -> dict:
    """The weighting in the saved models: weighting_archive.csv and verdict_B.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part C: the weightings

def measure_values(measure: str, off_counts: pd.DataFrame, family_breadth: pd.Series) -> pd.Series:
    """A KC's breadth by `measure` (indexed by KC id): 'F' the family's breadth (`family_breadth`), 'S' the summed spike
    count over the odors with two or more spikes in `off_counts`, 'N' the number of those odors."""
    raise NotImplementedError


def weighting_factors(m: pd.Series, k: float) -> pd.Series:
    """((m + 1) / mean over the KCs of (m + 1)) ** k per KC (as `E10_joint_fit.breadth_weights` forms it)."""
    raise NotImplementedError


def apl_kc_weights(wnorm_aplkc: pd.Series, m: pd.Series, k: float) -> pd.Series:
    """The APL->KC weights (per claw) of a weighting: the connectome's normalized weights times the KC's factor, by
    `E10_joint_fit.breadth_weights` (so F at exponent 2 gives the fit's weights exactly)."""
    raise NotImplementedError


def part_c0(out: Path, arm: str, wirings: Sequence[str], n_workers: int = 16) -> dict:
    """C0's identities: F at exponent 2 through the new code path reproduces the saved fit (C1) or E87's D2 run with
    per-odor APL scales (C2) exactly. identity_C.json."""
    raise NotImplementedError


def part_c1(out: Path, wirings: Sequence[str], weights: Optional[Sequence[Tuple[str, float]]] = None) -> pd.DataFrame:
    """C1: each weighting on each wiring at its fit's parameters, the APL scale (TIED_ON) bisected to the fit's APL-on
    count; the APL-off run the saved one. Resumable: runs_C1.csv (or runs_C3.csv for T)."""
    raise NotImplementedError


def part_c2(out: Path, wirings: Sequence[str], weights: Optional[Sequence[Tuple[str, float]]] = None,
            n_workers: int = 16) -> pd.DataFrame:
    """C2: each weighting on the base with the flies' recruitment (E88's per-odor thresholds), each odor's APL scale
    bisected to the control flies' count for that odor. Resumable: runs_C2.csv (or runs_C3.csv for T)."""
    raise NotImplementedError


def score_c(out: Path) -> dict:
    """The readings, the co-response, the choice, the gate and the verdicts: readings_C.csv, coresponse_C.csv,
    choice_C.csv, gate_C.csv, verdict_C.json."""
    raise NotImplementedError


def part_c3(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """C3: the confirmation on T for each base and direction with something to confirm, and the permuted controls."""
    raise NotImplementedError


def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv of C's runs."""
    raise NotImplementedError


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B', 'C0', 'C1', 'C2', 'C3', 'Cscore', 'variants', 'table'), required=True)
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = H.results_dir(EXPERIMENT)
    raise NotImplementedError(args.part)


if __name__ == '__main__':
    main()
