#!/usr/bin/env python3
"""E92: the step from claws to somata with APL on, read like for like: the flies' somata with their trial noise removed
against the models' somata as the recordings are normalized, the line along which every saved form of APL's inhibition
trades the ester pairs against the alcohol pairs, and then APL's inhibition of the KCs as a gain on the claw's drive.

The entry is item E92 of experiment_priorities.md (proposed by hal 2026-11-13 19:04; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-13 18:38). C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the
0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit`). S: real, P1, P2, M1. T: P3, P4,
M2, M3. Halves A / B of the recordings E77's (`recording_halves()`, seed 0). The 66 pairs of the panel's 12 megamat
odors in E90's groups (ester 15, alcohol 6, ester with alcohol 24, other 21, E78's six near pairs, all 66), the six near
pairs also as the four ester near pairs together and 1-5ol with 1-6ol and 6al with t2h singly (GROUPS).

This module holds part A (no simulations; it reads E90's archived runs of C2's family, Yang's per-trial recordings and
E89's archived runs, and nothing of E91's). Parts B to F are in experiments.E92_gain; their names here are stubs.

Parts (CLI `--part`):
- A1: the flies' somata (somata_flies.csv): per recording of Yang's control flies and flies with APL output blocked
  (TNT), in Yang's z and in baseline-noise units, the 66 pair distances as trial means (E90's reading) and with the
  trial noise removed (E81's estimate from products of different repeats, `H.noise_free_pair_distances`), as cosine and,
  on single trials centred over the KCs, as Pearson; the mean over recordings per pair and group with 90 % intervals from
  N_BOOT resamples of the recordings. Remy's wild-type somata (per-fly correlation matrices only) as measured, flagged.
- Acheck: the known-truth check of the noise removal (estimate_check.csv): E89's pseudo-recordings of the model
  (`ModelReading89`: the eight wirings' runs read under E73's per-KC session z-score, three trials, at E89's noise
  levels, with a share of the noise common to a repeat's odors) through the same estimate, against the model's
  noise-free group distances; a group passes in a measure and units if the population's estimate is within CHECK_TOL of
  the truth at every one of CHECK_LEVELS in both conditions. A group that fails is not read with the noise removed.
- A2: the model's somata as the recordings are normalized (somata_model.csv): the eight wirings, B5's point at the
  flies' 2 s pulse and D1's four firing-ceiling fits, read under E73's per-KC session z-score (calibrated on each run's
  APL-on counts, saturating and linear maps), in the reading's z (as Yang's z) and in baseline-noise units (the
  amplitude over the common noise), without noise (against A1's noise-removed values) and through the reading's
  measurement noise in expectation (against the flies' trial means), as Pearson and cosine; beside, raw counts and the
  calibrated readout (E90's soma readings); per run and as populations (the eight; S; T; D1's four; C2 on D1's four).
- A3: the verdict (verdict_A.json, verdict_rows_A.csv, steps_like_for_like.csv): by E90's rule (`E90.differs`: the
  population's group mean VERDICT_DIFF or more from the flies', outside their 90 % interval, at least 6 of the 8 wirings
  on that side) per group, condition and like-for-like reading, the somata first (no claw reading), then the two steps
  with E90's claws on both sides; per group "with APL on", "without APL", "both" or "neither", with the sign; a statement
  holds only if every like-for-like reading that passed the check agrees.
- A4: the line (line.csv, line.json): for every archived variant of C2's family at 6.8 % inside the bands, the ester and
  the alcohol pairs' distances with APL on, each minus the mean over the 120 pairs', on raw counts and under the
  calibrated readout, as populations where four wirings exist; the line fitted through E89's 13 weightings (C1, S's
  four wirings); a variant is off the line toward the flies if its residual is on the flies' side by more than the
  weightings' scatter about the line, with the flies' point from half A and from half B; such variants are named with
  their pattern part beside. No variant is chosen by this.
- A: A1, Acheck, A2, A3 and A4 in order; table: comparison_table.md and the figures.
- B, C, D, E, F: implemented in experiments.E92_gain.

Functions for experiments.E92_gain: `somata_readings(on, off)` (one run's somata in A2's readings, per group and
reading) and `line(readout)` / `place_on_line(est_rel, alc_rel, ln)` / `relative_distances(d_on)` (A4's line and one
variant's place against it).

Outputs: results/E92_steps_like_for_like/ (somata_flies.csv, estimate_check.csv, somata_model.csv, verdict_rows_A.csv,
steps_like_for_like.csv, verdict_A.json, line.csv, line.json, line_variants.csv, comparison_table.md, figures,
provenance_<part>.json); per-recording and per-run pair tables under outputs/E92_steps_like_for_like/.

Usage (from the repository root):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E92_steps_like_for_like.py --part A
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E92_steps_like_for_like.py --part A4 --workers 8
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E90_model_stages as E90

EXPERIMENT = 'E92_steps_like_for_like'

# ---------------------------------------------------------------------------
# groups

NEAR_ESTER = (('eb', 'ep'), ('IaA', 'pa'), ('2h', 'pa'), ('2h', 'IaA'))
NEAR_SINGLE = {'near_15ol_16ol': ('1-5ol', '1-6ol'), 'near_6al_t2h': ('6al', 't2h')}
GROUPS = E90.GROUPS + ('near_ester',) + tuple(NEAR_SINGLE)
GROUP_NAMES = dict(E90.GROUP_NAMES, near_ester='the four ester near pairs', near_15ol_16ol='1-5ol with 1-6ol',
                   near_6al_t2h='6al with t2h')

# ---------------------------------------------------------------------------
# A1 and the check

UNITS = ('z', 'noise')                  # Yang's whole-recording z; baseline-noise units
MEASURES = ('pearson', 'cosine')
FLY_READINGS = ('trial_mean', 'noise_removed')
N_BOOT = E90.N_BOOT
CI = E90.CI
BOOT_SEED = E90.BOOT_SEED
CHECK_TOL = 0.03

# ---------------------------------------------------------------------------
# A2

MAPS = ('saturating', 'linear')         # E73's two amplitude maps
NOISE = ('none', 'expected')            # without noise / through the reading's measurement noise in expectation

# ---------------------------------------------------------------------------
# A3

VERDICT_DIFF = E90.VERDICT_DIFF
VERDICT_SAME_SIDE = E90.VERDICT_SAME_SIDE


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def group_masks(pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """E90's group masks over `pairs` (canonical order) plus the near pairs split (GROUPS)."""
    raise NotImplementedError


def group_means(pairs: pd.Series) -> Dict[str, float]:
    """A pair Series' mean over each group's pairs (GROUPS; NaN pairs left out)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A1: the flies

def centred(trials: np.ndarray) -> np.ndarray:
    """Single trials (repeats, KCs, odors) with each trial's mean over the KCs removed, so that their cosine is Pearson's
    r."""
    raise NotImplementedError


def recording_distances(trials: np.ndarray, odors: Sequence[str]) -> Dict[Tuple[str, str], pd.Series]:
    """One recording's 66 megamat pair distances {(reading, measure): pair Series} for FLY_READINGS x MEASURES."""
    raise NotImplementedError


def flies_samples(units: str) -> Dict[Tuple[str, str, str], pd.DataFrame]:
    """{(genotype key 'ctrl' / 'tnt', reading, measure): recordings x the 66 pairs} in `units`; cached."""
    raise NotImplementedError


def part_a1() -> pd.DataFrame:
    """somata_flies.csv: per genotype, units, reading and measure, per pair and group the mean over recordings, its 90 %
    interval and SD from N_BOOT resamples of the recordings; Remy's wild type as measured (flagged)."""
    raise NotImplementedError


def flies_group(units: str, genotype: str, reading: str, measure: str) -> pd.DataFrame:
    """The flies' group rows of somata_flies.csv: group -> mean, lo, hi."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the known-truth check

def estimate_of_trials(trials: np.ndarray, measure: str) -> pd.Series:
    """The noise-removed 66 pair distances of one recording's single trials (repeats, KCs, odors in APL16 order)."""
    raise NotImplementedError


def check_wiring(label: str) -> List[dict]:
    """One wiring's known-truth rows: per condition, CHECK_LEVELS, units and measure the group means of the estimate
    (mean over the pseudo-recordings) and of the noise-free truth."""
    raise NotImplementedError


def part_check(workers: int = 8) -> pd.DataFrame:
    """estimate_check.csv (per wiring and the population) and check_passes() per group, measure and units."""
    raise NotImplementedError


def check_passes(ec: Optional[pd.DataFrame] = None) -> Dict[Tuple[str, str, str], bool]:
    """{(group, measure, units): passes} from estimate_check.csv's population rows."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A2: the model's somata

def model_calibrations(on: pd.DataFrame) -> Dict[str, dict]:
    """E73's calibration of the per-KC session z-score on one run's APL-on counts, per map (MAPS)."""
    raise NotImplementedError


def expected_distances(center: np.ndarray, sd: np.ndarray, measure: str, odors: Sequence[str]) -> pd.Series:
    """The 66 pair distances of units x odors values `center` through independent measurement noise of per-unit SD `sd`
    (a trial mean's), in expectation as a ratio of expectations (`zcell_pair_distances`' closed form; for Pearson the
    values centred over the units, whose noise adds (1 - 1 / n) sum sd^2 to each odor's squared length)."""
    raise NotImplementedError


def somata_pairs(on: pd.DataFrame, off: pd.DataFrame, cal: Optional[Dict[str, dict]] = None) -> pd.DataFrame:
    """One run's somata in every reading of A2: rows cond, reading, odor_a, odor_b, value (the 66 megamat pairs);
    readings '<map>_<units>_<noise>_<measure>' and 'raw_<measure>', 'cal_<measure>'."""
    raise NotImplementedError


def somata_readings(on: pd.DataFrame, off: pd.DataFrame, cal: Optional[Dict[str, dict]] = None) -> pd.DataFrame:
    """One run's somata group distances in A2's readings (for experiments.E92_gain): `on` / `off` spike counts (KCs x
    odors, any odor labels `H.odor_abbrev` reads); returns rows cond ('on' / 'off'), reading, then one column per group
    (GROUPS)."""
    raise NotImplementedError


def part_a2(workers: int = 8) -> pd.DataFrame:
    """somata_model.csv: per run (the eight wirings, B5's point, D1's four fits) and population, condition and reading
    the group means."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A3: the verdict

def like_for_like() -> List[Tuple[str, str, str, str]]:
    """The like-for-like comparisons: (model reading prefix '<map>_<units>_<noise>', flies' units, flies' reading,
    measure)."""
    raise NotImplementedError


def soma_verdict_rows(sm: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool]) -> pd.DataFrame:
    """Every (group, condition, reading) comparison of the model's somata with the flies' (`E90.differs`)."""
    raise NotImplementedError


def step_verdict_rows(sm: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool]) -> pd.DataFrame:
    """The two steps like for like (somata minus E90's claws, both sides), every reading; rows as `soma_verdict_rows`."""
    raise NotImplementedError


def condition_verdict(on: dict, off: dict) -> str:
    """'with APL on', 'without APL', 'both', 'neither' or 'not resolved' from two conditions' element results
    (`E90.element_result`)."""
    raise NotImplementedError


def part_a3() -> dict:
    """verdict_A.json: per group the somata's and the steps' verdicts, with the signs, over the readings that passed the
    check; which of the entry's outcomes (a) to (c) holds."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A4: the line

def relative_distances(d_on: pd.Series) -> Tuple[float, float]:
    """(ester, alcohol): the ester and the alcohol groups' mean APL-on distances, each minus the mean over all of the
    Series' pairs (the 120 pairs)."""
    raise NotImplementedError


def fit_line(est_rel: Sequence[float], alc_rel: Sequence[float]) -> dict:
    """The least-squares line of the ester on the alcohol relative distances: slope, intercept, r, n, the residuals'
    SD (`scatter`, n - 2 degrees of freedom) and largest absolute residual."""
    raise NotImplementedError


def line(readout: str = 'raw') -> dict:
    """A4's line through E89's 13 weightings (C1, S's four wirings as populations) on `readout` ('raw' / 'cal'), with the
    flies' points (control flies, all recordings and halves A and B, trial means and noise removed); from line.json if
    part A4 has run, else computed."""
    raise NotImplementedError


def place_on_line(est_rel: float, alc_rel: float, ln: dict) -> dict:
    """One variant's place against a line (`line`): its residual, the residual signed toward the flies per flies' point,
    and whether it is off the line toward the flies (by more than the scatter, against half A and against half B)."""
    raise NotImplementedError


def a4_variants() -> pd.DataFrame:
    """The archived variants A4 reads: experiment, folder, family key, wiring, tag, inside the bands."""
    raise NotImplementedError


def part_a4(workers: int = 8) -> pd.DataFrame:
    """line.csv, line.json, line_variants.csv."""
    raise NotImplementedError


def comparison_table() -> None:
    """comparison_table.md and the figures."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# B to F: implemented in experiments.E92_gain

def part_b(*args, **kwargs):
    """Part B, the gain at fixed scales on S: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_c(*args, **kwargs):
    """Part C, choice and confirmation: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_d(*args, **kwargs):
    """Part D, the family fitted again: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_e(*args, **kwargs):
    """Part E, on a pass only: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_f(*args, **kwargs):
    """Part F, the summed form on the same readings after E91's note: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'A1', 'Acheck', 'A2', 'A3', 'A4', 'table', 'B', 'C', 'D', 'E', 'F'])
    ap.add_argument('--workers', type=int, default=8, help='Acheck / A2 / A4: worker processes')
    a = ap.parse_args()
    if a.part == 'A':
        part_a1()
        part_check(a.workers)
        part_a2(a.workers)
        part_a3()
        part_a4(a.workers)
        comparison_table()
        return
    {'A1': part_a1, 'Acheck': lambda: part_check(a.workers), 'A2': lambda: part_a2(a.workers), 'A3': part_a3,
     'A4': lambda: part_a4(a.workers), 'table': comparison_table, 'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e,
     'F': part_f}[a.part]()


if __name__ == '__main__':
    main()
