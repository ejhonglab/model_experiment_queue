#!/usr/bin/env python3
"""E88: the error common to all wirings split into who responds and how strongly.

With sub-threshold responses set to 0, the similarity of two odors (one minus their cosine distance) is the product of
three factors (experiment_priorities.md, E88; planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-13 00:45):
- R, recruitment: the geometric mean of the two odors' response fractions;
- L, co-response: how many times more often than chance a KC that responds to one odor responds to the other
  (L = O / R, O one minus the distance of the binarized responses);
- G, the amplitude factor: G = s / O, the mean product of the shared responders' two responses over the root mean
  squares of the two odors' responses over their responders.
So s = R * L * G exactly. For the flies each of s, O and R is the mean over a genotype's recordings, and L and G follow;
for a model population the means are over its wirings and pseudo-recordings. The per-pair APL effect is s with APL off
minus s with APL on. A cell is that effect with the named factors taken from the flies and the others from the model,
in both conditions.

Readings (both sides read like for like):
- 'r1': trial means. The flies' trial means in z at 1 z; the model under E73's per-KC reading (`H.session_zscore`,
  saturating map, calibrated on each run's APL-on counts by E74's `_calibrate_zcell`) on N_PSEUDO pseudo-recordings of
  N_TRIALS trials through its measurement noise (numpy default_rng(SEED), the APL-on draws first: E74's `pseudo_noise`),
  the criterion 1 z on the trial means.
- 'r1b': as 'r1' with each condition's measurement noise multiplied until an odor's responder sets agree across repeats
  (the mean binary cosine between two different repeats of one odor) as the flies' do in that condition, and each
  condition's criterion re-set to keep 'r1''s response rates.
- 'r2': the trial noise removed by the same estimator on both sides (`H.noise_free_pair_distances`'s, graded and
  binarized): the flies' single trials at their single-trial criterion (E87's `single_trial_criterion`, 1.057 z), the
  model's pseudo-recordings' single trials at the criterion at which its APL-on single-trial rate is its 'r1' APL-on
  rate (both conditions). 'r2b' (the estimator check only): the same at 'r1b''s noise, each condition's single-trial
  criterion at 'r1''s rate.
- 'r3raw' / 'r3cal' (E87's convention, beside): the deterministic model, raw spike counts (a response two or more
  spikes) or E65's calibrated readout, against the flies with the trial noise removed ('r2''s flies).

Parts (CLI `--part`):
- A (no simulations): A0 the estimators on known truth (the model's noise-free session z-scores); A1 the factors R, L,
  G, O and s per pair, condition and reading for the flies (all recordings, the six typical, each recording; 90 %
  intervals over N_BOOT resamples of the recordings; also in baseline-noise units), each wiring and each population;
  A2 the cells against the flies' effect in the same reading, with floors; A3 the same on E87's D2 runs (E80's per-odor
  recruitment on S, three arms) against the unchanged fits of S; verdict_A.json by the thresholds below.
- B (no simulations; readings 1 and 1b): B1 the changes of L and R under APL per group; B2 the side of a difference in
  the co-response (a per-odor cut of the responses without APL, and the flies' chosen form, against what APL leaves);
  B3 the shared responders' amplitudes; B4 the model's KCs that carry a group's overlap; verdict_B.json.
- C (no simulations): every archived variant (E87's `archive_variants`, E86's excepted), E85's 72 model flies and both
  families by the factors, under reading 1 (each variant's own calibration) and reading 3; the stand-out rule;
  verdict_C.json.
- D1 / D2 (simulations): not implemented here (stubs).
- table: comparison_table.md from the parts' tables.
Outputs: results/E88_responder_sets/ (comparison_table.md first); large outputs under outputs/E88_responder_sets/.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E88_responder_sets'

# the readings' settings (E74's pseudo-recordings)
N_PSEUDO, N_TRIALS, SEED = 10, 3, 43
N_BOOT = 1000                  # resamples of the recordings for the flies' intervals and the cells' floors
INTERVAL = (5.0, 95.0)         # the flies' 90 % intervals
B3_MIN_SHARED = 5              # B3: a pair's rank correlation needs at least this many shared responders

# the thresholds, fixed before any run (the entry's Design)
EST_TOL = 0.03                 # A0: a factor whose group mean is off by more than this in a reading is not read there
G_MIN = 0.015                  # A: g_sets (g_amp) at least this ...
G_RATIO = 0.5                  # ... and the other at most this share of it
GROUP_RATIO = 0.5              # A: a group's |centred error| with one factor from the flies at most this share of the own
D2_BASE_GAIN = 0.010           # A3: what is left on the base with the flies' recruitment is the co-response
D2_BASE_ARM = 'thrapl'         # the arm with the flies' recruitment in both conditions (per-odor APL scales)
B1_FACTOR = 1.10               # B1: a group's change differs by this factor, outside the flies' interval
B2_LARGE = 2.0 / 3.0           # B2: differ by at least two thirds of the difference under APL
B2_SMALL = 1.0 / 3.0           # B2: agree to within a third of it
B4_TOP = 0.25                  # B4: the top quarter of the responders' breadth weights ...
B4_SHARE = 2.0 / 3.0           # ... holds at least two thirds of the KCs a cut would keep and APL removes
C_BIN_GAIN = 0.008             # C: a variant's binarized pattern part lower than its baseline's by this
C_HALFWAY = 0.5                # C: or the alcohol group's co-response change at least halfway to the flies'

GROUPS = ('nearest4', 'ester', 'alcohol', 'six', 'other')
FACTORS = ('R', 'L', 'G', 'O', 's')
CELLS = {'own': (), 'sets': ('R', 'L'), 'amp': ('G',), 'R': ('R',), 'L': ('L',), 'flies': ('R', 'L', 'G')}
MODEL_READINGS = ('r1', 'r1b', 'r2', 'r2b', 'r3raw', 'r3cal')
FLIES_READING = {'r1': 'r1', 'r1b': 'r1', 'r2': 'r2', 'r2b': 'r2', 'r3raw': 'r2', 'r3cal': 'r2'}
VERDICT_READINGS = ('r1', 'r2')
POPULATIONS = ('b2', 'fc', 'S_unchanged', 'D2_thr', 'D2_throff', 'D2_thrapl')


# ---------------------------------------------------------------------------
# the factors

@dataclass
class Factors:
    """One condition's s, O and R per pair (120-vectors in triu order of APL16_ODORS); L = O / R and G = s / O."""
    s: np.ndarray
    O: np.ndarray
    R: np.ndarray

    @property
    def L(self) -> np.ndarray:
        raise NotImplementedError

    @property
    def G(self) -> np.ndarray:
        raise NotImplementedError

    def get(self, name: str) -> np.ndarray:
        raise NotImplementedError


@dataclass
class RecFactors:
    """One condition's per-recording s, O and R (recordings x 120; R NaN where O is, so that the three means are over
    the same recordings), with `mean(rows)` -> Factors."""
    s: np.ndarray
    O: np.ndarray
    R: np.ndarray

    def mean(self, rows: Optional[np.ndarray] = None) -> Factors:
        raise NotImplementedError


def group_masks() -> Dict[str, np.ndarray]:
    """Boolean masks over the 120 pairs in triu order: 'nearest4' (E77's NEAREST), 'ester' (15), 'alcohol' (6), 'six'
    (E78's SIX), 'other' (the 99 pairs outside the ester and alcohol groups), 'all'."""
    raise NotImplementedError


def factors_of_matrix(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(s, O, R) of one KC x odor matrix whose sub-threshold values are already 0."""
    raise NotImplementedError


def factors_of_trials(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(s, O, R, reliability) of one recording's single trials (repeats, KCs, odors) whose sub-threshold values are
    already 0, with the trial noise removed (`H.noise_free_pair_distances`'s estimator, graded and binarized); R from
    the single-trial response fractions; reliability per odor (the mean binary cosine between two different repeats)."""
    raise NotImplementedError


def cell_effect(model: Dict[str, Factors], flies: Dict[str, Factors], take: Sequence[str]) -> np.ndarray:
    """The per-pair effect (s with APL off minus s with APL on) of the product R * L * G with the factors in `take`
    from the flies and the others from the model, in both conditions."""
    raise NotImplementedError


def scores88(pred: np.ndarray, data: np.ndarray) -> dict:
    """E87's `scores` plus every group's mean centred error (err_<group>)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the flies

def flies_factors(unit: str, reading: str, tnt: Optional[Sequence[tuple]] = None) -> Dict[str, RecFactors]:
    """{'on': control flies, 'off': APL output blocked (TNT)} per-recording factors in `unit` ('z' / 'noise') and
    reading ('r1' trial means at the criterion, 'r2' single trials at the single-trial criterion with the trial noise
    removed); `tnt` restricts the TNT recordings (default all)."""
    raise NotImplementedError


def flies_reliability(unit: str = 'z') -> Dict[str, float]:
    """The flies' mean binary cosine between two different repeats of one odor at the single-trial criterion, per
    condition (the target of reading 1b)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the models

class ModelReading:
    """One run's APL-on and APL-off spike counts read in every reading of MODEL_READINGS, with the noise-free truth of
    each noisy reading (the session z-scores without measurement noise at the criterion that gives the reading's
    response rate in that condition)."""

    def __init__(self, o: Dict[str, tuple], label: str = '', cal: Optional[dict] = None, K: Optional[float] = None,
                 rel_target: Optional[Dict[str, float]] = None, readings: Sequence[str] = MODEL_READINGS):
        raise NotImplementedError

    def factors(self, reading: str) -> Dict[str, RecFactors]:
        raise NotImplementedError

    def truth(self, reading: str) -> Dict[str, Factors]:
        raise NotImplementedError


def population(members: Sequence[ModelReading], reading: str) -> Dict[str, Factors]:
    """A population's factors in a reading: s, O and R the means over its wirings and pseudo-recordings."""
    raise NotImplementedError


def model_populations() -> Dict[str, List[Tuple[str, Dict[str, tuple], Optional[float]]]]:
    """{population: [(wiring, outputs, K of the calibrated readout)]}: 'b2' (E67's fit on E72 part B's eight wirings),
    'fc' (the firing-ceiling family), 'S_unchanged' (b2 on S), 'D2_<arm>' (E87's D2 runs, S's wirings; K the fit's)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A

def part_a0(out: Path, readings_by_pop: dict) -> pd.DataFrame:
    """A0: each factor's estimate minus its truth by group and condition on the breadth-2 population, at the readings'
    noise ('r1', 'r2') and at reading 1b's ('r1b', 'r2b'). estimator_check.csv."""
    raise NotImplementedError


def part_a1(out: Path, readings_by_pop: dict) -> pd.DataFrame:
    """A1: the factors per pair, condition and reading (flies, wirings, populations). factors.csv, factors_pairs.csv."""
    raise NotImplementedError


def part_a2(out: Path, readings_by_pop: dict) -> pd.DataFrame:
    """A2 and A3: the cells against the flies' effect in the same reading, with floors. cells.csv."""
    raise NotImplementedError


def verdict_a(out: Path) -> dict:
    """The thresholds of A, from cells.csv and estimator_check.csv. verdict_A.json."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B

def part_b(out: Path) -> dict:
    """B1 to B4 and B's verdict: changes.csv, cut.csv, shared_amplitudes.csv, shared_kcs.csv, verdict_B.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part C

def part_c(out: Path, n_workers: int = 8) -> dict:
    """The archive by the factors: archive_scan.csv, verdict_C.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part D (simulations; not implemented here)

def part_d1(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """D1: the flies' chosen form as a static inhibition in place of APL on S's four wirings with the flies' per-odor
    recruitment (E87's D2 thresholds, re-solved with `B3Pool`). Not implemented: part D is a later step."""
    raise NotImplementedError('E88 part D1 is not implemented yet')


def part_d2(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """D2: up to three handles that stand out in C, on S's four wirings. Not implemented: part D is a later step."""
    raise NotImplementedError('E88 part D2 is not implemented yet')


def comparison_table(out: Path) -> str:
    """comparison_table.md covering A to C, from the parts' tables."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B', 'C', 'D1', 'D2', 'table'), required=True)
    ap.add_argument('--workers', type=int, default=8, help='C: worker processes')
    args = ap.parse_args()
    out = H.results_dir(EXPERIMENT)
    if args.part == 'A':
        part_a(out)
    elif args.part == 'B':
        part_b(out)
    elif args.part == 'C':
        part_c(out, args.workers)
    elif args.part == 'D1':
        part_d1(out)
    elif args.part == 'D2':
        part_d2(out)
    comparison_table(out)


if __name__ == '__main__':
    main()
