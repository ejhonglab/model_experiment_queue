#!/usr/bin/env python3
"""E91: APL's inhibition on the KC's summed drive instead of claw by claw.

The entry is item E91 of experiment_priorities.md (proposed by atlas 2026-11-13 11:38; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-13 11:00). In olfsysm each claw's drive x (its PN input minus APL's
inhibition of the claw) is clipped at 0 before the claws are summed. olfsysm's `p.kc.claw_inh_pool` (s, the pooled
share; `fit_mb_model(claw_inh_pool=...)`) makes a KC's summed drive the sum over its claws of max(x, 0) + s * min(x, 0):
s = 0 is every model so far, s = 1 is APL's inhibition on the summed drive. The saved claw states and APL's drive from
the claws stay the clipped ones, so with APL off nothing changes, and at fixed scales APL's activity and the boutons are
those of s = 0: only the KCs' summed drive loses s times the inhibition its claws do not absorb.

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's
(`recording_halves()`, seed 0), C / D the second split's (seed 1; E89's `targets89`).

Parts (CLI `--part`):
- 0, the identities: s = 0 through the new path (the parameter set on each S wiring's simulator, and a simulator built
  with `claw_inh_pool=0.0`) reproduces the saved fits with no count differences; at s = 1 the APL-off runs equal the
  saved ones; E89's C1 point F at exponent 1 on P1 reproduced exactly (APL scale and counts); the parameter set on the
  simulator equals the simulator built with it (s = 1, APL on). identity.json.
- A (no new parameter): on S at the fits (APL on, one run a wiring with the bouton rates and APL's activity saved), per KC
  and odor the inhibition its claws do not absorb (the sum over its claws of the part of APL's inhibition above the
  claw's drive; window mean over the odor pulse), relative to its summed drive and to its margin over threshold; by
  class (the responders with APL on; E88's classes of the ester and the alcohol groups' shared responders; number of
  claws; quartile of the breadth weight; KC type); the prediction of the responders s = 1 removes at fixed scales (each
  KC's summed drive, with and without the unabsorbed inhibition, through the KC's leaky integration, no spikes), read as
  the change of the per-pair APL effect, checked against B (`--part Acheck`, after B1). unabsorbed.csv,
  prediction.csv; per KC and odor under outputs/E91_summed_drive_inhibition/A/.
- B1: the fits as they are, s of 0.5 and 1 at the family's measure F with exponents 2, 1 and 0.5 and with no weighting
  (8 settings, 32 points), the tied APL scale bisected to each wiring's APL-on count (E89's `weighting_reading`); plus
  s = 1 at F exponent 2 at the fit's scales (no bisection; 4 points, for part A's check). E89's saved C1 points are the
  references at s = 0.
- B2: the base with the flies' recruitment (E88's per-odor thresholds; each odor's APL scale set again, as E89's
  `part_c2`), s = 1 at exponents 2 and 1 and with no weighting (12 points). E89's saved C2 points are the references.
- B3, for the record: olfsysm's `allow_net_inh_per_claw=True` at exponent 2 on the fits as they are (4 points).
- co: the alcohol and the ester groups' co-response change under APL beside every S point (E89's `coresponse_c`,
  readings 1 and 1b and the geometric mean of E88's two estimates at L1 and L3).
- C: the choice per base and half (the lowest pattern part of S's population under the calibrated readout among the
  settings inside the bands); a choice with s above 0 that gains CHOICE_MIN_GAIN or more over its base's reference goes
  to T with E83's gate against the other half, on both splits (`gate91`); the control at s = 0 (a static inhibition at
  the soma, olfsysm's `ext_inh_kc`, equal for all odors to each KC's mean pooled unabsorbed inhibition at the chosen
  setting) on S and T; the equal rule (`equal_rule91`).
- score: readings.csv, coresponse.csv, choice.csv, gate.csv, equal.csv, verdict.json, figures; comparison_table.md.
- D: the family fitted again (`E10_joint_fit.py ... --claw-inh-pool 1`, 150 evaluations, from each wiring's own fit of
  C2's family, on S's four wirings; at exponent 1 and on T only by the entry's rules); `--part Dscore` scores them.
- E: only on a pass (not written).
- variants: variants.csv and model_pair_distances_long.csv, for the archive (`archive_tools/archive_model_outputs.py
  E91_summed_drive_inhibition`, then `archive_tools/migrate_model_output_params.py`).
- table: comparison_table.md.

Outputs: results/E91_summed_drive_inhibition/; spike counts and responses under outputs/E91_summed_drive_inhibition/
(B/, C/ for the runs, A/ for part A's per KC tables; D's fits at the folder's root, as E10_joint_fit.py writes them).

Usage (from the repository root; one 16-thread process at a time on viki):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E91_summed_drive_inhibition.py --part 0
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E87_inhibition_form as E87
from experiments import E88_responder_sets as E88
from experiments import E89_apl_breadth_weighting as E89

EXPERIMENT = 'E91_summed_drive_inhibition'

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

POOLS = (0.5, 1.0)                                  # B1's pooled shares
B1_WEIGHTINGS = (('F', 2.0), ('F', 1.0), ('F', 0.5), ('F', 0.0))
B2_WEIGHTINGS = (('F', 2.0), ('F', 1.0), ('F', 0.0))
B2_POOLS = (1.0,)
REFERENCE = E89.REFERENCE                           # ('F', 2.0): the family's weighting
ARMS = ('B1', 'B2')                                 # the fits as they are; the base with the flies' recruitment
RAW_CRIT = E87.RAW_CRIT                             # a model response is two or more spikes
PULSE = (0.0, 0.5)                                  # the odor pulse (s): the window of the unabsorbed inhibition
DRIVE_WINDOW = (0.0, 0.5)                           # E78's window for the claws carrying half of a responder's drive

# part A
PREDICT_RHO = 0.8                                   # the prediction stands at this correlation or more

# part C (the entry's rules)
CHOICE_MIN_GAIN = 0.0015
GATE_GAIN = 0.003                                   # (i) the mean gain over the directions
GATE_SINGLES_SHARE = 0.75                           # (ii) 12 of the 16 single-wiring cells (4 wirings x 4 directions)
GATE_RAW_SLACK, GATE_MEAN_SLACK = 0.005, 0.01       # (iii)
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075              # (iv)
CONTROL_RATIO = 0.5
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN, EQUAL_MAX_EXPONENT = 0.003, 0.005, 0.01, 1.0
HALVES = ('A', 'B', 'C', 'D')
B3_DIFF = 0.003                                     # outcome (i): B3 differs from s = 1 by more than this

# part D
D_POOL = 1.0
D_MAXFEV = 150
D_S_GAIN = 0.0015                                   # S's fitted population this far below UNCHANGED_S['all'] -> T's fits
D_PASS_GAIN = 0.003                                 # the population of eight this far below UNCHANGED_8
D_WIRINGS_LOWER = 6                                 # of the eight
D_RAW_REF, D_RAW_SLACK, D_MEAN_SLACK = 0.073, 0.005, 0.01
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
FLIES_MEAN = E89.FLIES_MEAN

_CACHE: Dict[tuple, object] = {}


# ---------------------------------------------------------------------------
# names

def variant_name(arm: str, w: Tuple[str, float], s: float = 0.0, net: bool = False, fixed: bool = False,
                 control: bool = False) -> str:
    """A setting's variant name: '<arm>_<weighting>' at s = 0 (E89's names with the arm B1 / B2), '_s<s>' for s above 0,
    '_net' for `allow_net_inh_per_claw`, '_fixed' at the fit's APL scales (no bisection), '_ctrl' for the control."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    """`variant_name`'s fields back: arm, measure, exponent, s, net, fixed, control."""
    raise NotImplementedError


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


# ---------------------------------------------------------------------------
# the simulators

def build_simulator91(label: str, **model_kws):
    """E83's `build_simulator` (the lab's input, the wiring's fit, its breadth) with extra `fit_mb_model` kwargs
    (`claw_inh_pool`, `ext_inh_kc`, ...); returns (simulate, info) as E83's does."""
    raise NotImplementedError


def set_pool(mp, s: float = 0.0, net: bool = False) -> None:
    """Set the pooled share (`mp.kc.claw_inh_pool`) and `allow_net_inh_per_claw` on a simulator's olfsysm parameters."""
    raise NotImplementedError


def claw_stats(rv, mp, W: np.ndarray, kc_of_claw: np.ndarray, wk_scaled: np.ndarray, odors: Sequence[int],
               responders: np.ndarray, thr: Optional[np.ndarray] = None, predict: bool = False) -> dict:
    """From one run's saved bouton rates and APL activity: per KC and odor the summed clipped drive, the unabsorbed
    inhibition, the number of claws carrying half of a responder's drive, and (with `predict`) the leaky-integrated
    peaks with and without the unabsorbed inhibition."""
    raise NotImplementedError


def unabsorbed_of(x: np.ndarray, kc_of_claw: np.ndarray, n_kc: int) -> Tuple[np.ndarray, np.ndarray]:
    """(summed clipped drive, unabsorbed inhibition) per KC and time step from the claws' net drive x (claws x time)."""
    raise NotImplementedError


def n_claws_half(drive: np.ndarray) -> int:
    """The number of a KC's claws (their window-mean clipped drives) that carry half of its summed drive."""
    raise NotImplementedError


def leaky_peak(drive: np.ndarray, dt: float, taum: float, v0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """(peak, the summed supra-threshold area is computed by the caller) of each row's leaky integration of `drive`
    (KCs x time; dV = (-V + drive) dt / taum, no spikes); returns (V, peak)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# parts

def part_0(out: Path) -> dict:
    raise NotImplementedError


def part_a(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def part_a_check(out: Path) -> dict:
    raise NotImplementedError


def part_b1(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    raise NotImplementedError


def part_b2(out: Path, wirings: Sequence[str], n_workers: int = 16) -> pd.DataFrame:
    raise NotImplementedError


def part_b3(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    raise NotImplementedError


def coresponse91(out: Path, n_workers: int = 4) -> pd.DataFrame:
    raise NotImplementedError


def readings91(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def choose91(rd: pd.DataFrame, arm: str, half: str) -> dict:
    raise NotImplementedError


def gate91(rows: Sequence[dict]) -> dict:
    raise NotImplementedError


def equal_rule91(rd: pd.DataFrame, arm: str) -> List[dict]:
    raise NotImplementedError


def control_matrix(unabsorbed: pd.DataFrame, s: float) -> pd.DataFrame:
    """The control's static inhibition at the soma: per KC s times its unabsorbed inhibition averaged over the odors,
    the same for every odor (KCs x odors)."""
    raise NotImplementedError


def part_c(out: Path, n_workers: int = 16) -> dict:
    raise NotImplementedError


def verdict91(out: Path) -> dict:
    raise NotImplementedError


def part_d(out: Path, stage: str = 'S2') -> List[str]:
    raise NotImplementedError


def score_d(out: Path) -> dict:
    raise NotImplementedError


def part_e(out: Path) -> None:
    """Part E (only on a pass by C on the fits as they are, or by D): the firing-ceiling family, S at the flies' 2 s
    pulse, S through the trial noise, the block and the natural-source alignment, for the passing setting. Not written:
    it runs only on a pass."""
    raise NotImplementedError('E91 part E runs only on a pass by part C or part D; it is not written yet')


def write_variants(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def comparison_table(out: Path) -> str:
    raise NotImplementedError


def provenance(out: Path, part: str, **extra) -> None:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('0', 'A', 'Acheck', 'B1', 'B2', 'B3', 'co', 'C', 'score', 'D', 'Dscore', 'E', 'variants', 'table'))
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--wirings', default='S')
    ap.add_argument('--stage', default='S2', help="D: 'S2' (S at exponent 2), 'S1' (exponent 1), 'T2' (T at exponent 2)")
    args = ap.parse_args()
    raise NotImplementedError(args.part)


if __name__ == '__main__':
    main()
