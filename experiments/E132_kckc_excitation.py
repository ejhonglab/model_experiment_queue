#!/usr/bin/env python3
"""E132: the KC-KC synapses with a cholinergic sign: recurrent excitation among KCs through the hemibrain's claw-to-claw
synapses, with a pooled inhibition and alone, read as a dose through zero. Steps 0 and 1 (no olfsysm change).

The entry is item E132 of experiment_priorities.md (moved into the list 2026-11-23; claimed by hal 21:00). "The family" is
C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S is real,
P1, P2, M1 and T is P3, P4, M2, M3. The family's bouton base pins multibouton_claws='first' (E10_joint_fit.base_kws).
"KC-KC at zero" is E110's part A (KC-KC 0 in both conditions, both response rates held at fixed scales).

The dose and the pooled neuron. A dose d is the KC-KC scale (`kc_kc_inh`, the simulator's 'wKCKC_scale') on the
loader's normalized hemibrain matrix HB (post x pre, mean row sum 1, `connectome_KCKC_weights`), in absolute units
(the family's fitted scales s are +23 to +46); a negative dose is excitation. The pooled neuron is E122's uniform neuron
U (every entry HB's total / N^2, the diagonal included), filtered with 50 ms as the matrix is; its scale b ('wU_scale')
is absolute, and b / s is E129's factor (the pooled neuron's weight as a fraction of the matrix's at the fitted scale).

How the pooled neuron goes in (the entry leaves the route to the code): as a rank-1 part of the KC-KC matrix, E122's
route. The simulator's KC-KC matrix is set before every call to d HB + b U (summed dense, then CSC, as fit_mb_model
builds it) with 'wKCKC_scale' 1 (`Runner`). At b 0 the matrix olfsysm receives is bit for bit the family's (s h + 0 is
s h, with HB's sparsity); at d 0 it is bit for bit E129's pooled arm's (b u + 0). The feedback-cell path of olfsysm
fa3543d (`mbc1_inh` with uniform cell weights) is the same in exact arithmetic but filters the pooled spikes once
instead of each KC's activity, so it is not E129's run in floating point; the check part measures it on 'real'
(`route_diagnostic`). The checks below must pass with 0 differing spike counts before any run.

Readings the entry leaves open, fixed here before any run:
1. The pooled hold (step 1's arms; E129's, extended to both conditions; `pooled_hold`): the threshold stays the fit's.
   With APL off (APL->KC and APL->bouton 0) the pooled scale b = s x f is bisected (E69's `bisect_factor` through
   pyolfsysm.tuning.hold_rates with thr_key 'wU_scale', from f 1, E129's SCALE_BRACKET 1/50 to 50 and SCALE_STEP 1.25)
   to the family's APL-off count; then the tied APL scales (APL->KC, KC->APL, bouton->APL; E83.TIED_ON) are bisected
   with APL on (E77's APL_BRACKET 0.1 to 10, APL_STEP 1.25) to the family's APL-on count, the pooled neuron at that b.
   U0's first stage is E129's step-1a bisection of arm P exactly, so its APL-off run must equal E129's S1_P_<w> run
   (0 differing spike counts, checked in the run; RuntimeError otherwise).
2. The threshold hold (step 2's X arms; E110's part A): no pooled neuron; E77_input_screen.hold_rates (the threshold
   bisected with APL off, then the tied APL scales with APL on).
3. A corner: a bisected scale at a bound of its bracket, or a response count more than RATE_TOL (0.2 %) from its target
   (E129's `held_ok`). An arm with corners on more than CORNER_WIRINGS_MAX (two) wirings gets no class.
4. The classes (E123's, E126's and E131's, unchanged; calibrated pattern part, S's and T's populations apart, against
   all recordings and each half A to D; `classify`): "below the family" (lower by WITHIN 0.003 or more on both
   populations and on 6 of 8 wirings on every target, the mean APL effect under the readout within MEAN_TOL 0.02 of the
   flies' +0.178, the raw-count loss not higher by RAW_LOSS_TOL 0.005 on both populations, every run in the bands);
   "equal" (E123's, one-sided); "costs" (higher by more than WITHIN on both populations against all recordings);
   "neither" otherwise; "corner (no class)" by reading 3. "Keeps the gain" (E123's flag): at or below the half-way point
   between the family and KC-KC at zero (0.0554 on S, 0.0577 on T against all recordings) on every target.
5. Two runs "differ" only by E126's rule (E129's `differ`): one lower on 7 of 8 wirings against all recordings and on
   both populations against all five targets. Step 1 reads each mirror dose against U0 (its dose 0) and against the
   family; the sign's own effect (a negative dose against its mirror and its dose 0) is step 2's.
6. The stop rule (`stop_rule`, read on S): a dose passes if (i) S's population against all recordings is not above the
   family's by more than STOP_SLACK (0.005), (ii) the largest spike count per KC and odor with APL off is at most
   MAX_COUNT_FACTOR (2) times the family's on each of the eight wirings, (iii) every run is in the four bands with
   corners on at most two wirings, and (iv, added on moving) it is lower than its own dose 0 (U0 for the pooled arms,
   E110's part A for the arm alone) on S's population against all recordings and on 3 of S's 4 wirings. The mirror's
   doses are read by it and reported; only step 2's negative doses are candidates for step 3.
7. Where the excitation would land (step 0, a description; `landing`): each KC's KC-KC input per odor is E129's step-0
   reader (s HB @ n / 0.5 s, Hz; n the family's saved spike counts of the condition), with APL on and off; responders
   are the saved responses (two or more spikes) of the same condition. For an ordered pair (a, b): the mean input under
   odor a onto the KCs that respond to b and not to a, over the mean input under a onto the KCs that respond to
   neither; a pair's ratio is the mean of its two directed ratios. Group means over E90's groups of the megamat pairs
   (ester 15, alcohol 6, ester with alcohol 24, other 21) and over all 120 pairs, per wiring and for S4, T4 and all8 as
   the mean of their wirings' group means.
8. The matrix's side of stability (step 0; `stability`): the largest real eigenvalue (the Perron root; HB is
   nonnegative) of HB and of each KC type's diagonal block (γ, α′β′, αβ), and how much of the leading right and left
   eigenvectors is uniform over the KCs (cos² of their angle with the all-ones vector), the part a pooled neuron can
   cancel; beside it each wiring's s times the root and the doses 2, 5 and 10 times it.
9. The slope toward zero (step 0; `slopes`): under the threshold hold, the family (s), the matrix at s / 2 (E131's
   H50: fixed scales, both rates re-held) and KC-KC at zero (E110's part A). Per wiring, each reading's change per unit
   of scale on the segments s -> s / 2 and s / 2 -> 0; per population the same over the mean of its wirings' s. Carried
   through zero (the expectation on record for the excitation alone) the half->zero segment's slope: r(0) - slope x d at
   d 2, 5, 10 (doses -2, -5, -10).
10. Late firing (`late_share`): the share of spikes at or after the pulse's end (0.5 s), with APL on and off. The
    family's from E121's step-0 rebuild (`E121_onset_race.run_wiring`, committed code; part late); each step-1 run's
    from its held point simulated again with spike times, which must reproduce the held counts.
11. Sizes (per run): b / s; the threshold; the mean over KCs and odors of each part's input (d HB @ n / 0.5 s and
    b U @ n / 0.5 s, Hz) with APL on and off, over the fit's threshold; the largest count per KC and odor and the KCs
    with 20 or more spikes to some odor (E110's run_extras).

Steps 2 to 4 (the excitation, the shuffle control and the male calyces, the refit) need the olfsysm option for negative
KC-KC weights, held until Yang's OK: their arms are named here and the runner refuses any negative dose
(`require_nonnegative`), whatever the matrix's entries.

Parts (CLI `--part`):
- check (simulations, `--wirings`): (i) the runner at the fit's parameters (matrix s, pooled 0) reproduces the
  wiring's saved family run with 0 differing spike counts and responses, APL on and off; (ii) at matrix 0 and the
  pooled neuron at E129's bisected factor (its runs table) the APL-off run equals E129's saved S1_P_<w> run; on 'real'
  (iii) the route diagnostic (the feedback-cell path at the same b, its differing counts against E129's run, recorded,
  not a gate). check_<host>_<lane>.json; the runs refuse until every wiring has passed (i) and (ii).
- step0 (no simulation): step0_references.csv (readings of the references; `references`), step0_slopes.csv,
  step0_landing_pairs.csv, step0_landing.csv, step0_stability.csv, step0.json, provenance_step0.json.
- late (simulations): the family's late firing from E121's rebuild, late_family.csv (spike times under
  outputs/E132_kckc_excitation/late/).
- runs (simulations, resumable, per lane `--wirings`; a point claimed by mkdir of results/E132_kckc_excitation/claims/
  <tag>/, `--reclaim` to take over a dead process's): step 1's U0, UI2, UI5, UI10. Rows in runs_<host>_<lane>.csv,
  provenance in provenance_<host>_<lane>.json, runs and spike times saved locally (outputs/E132_kckc_excitation/step1/)
  and in the model-output archive (data/model_outputs/E132_kckc_excitation/step1/).
- score (no simulation): readings.csv, runs.csv, verdict.json, comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts, tagged by the run tag, the archive's file
  name, so kcq/apl_effect.py scores the realism bands).
- plan: the points a lane would run.

The holds (CLAUDE.md, "Connectome data"): only the model's hemibrain wiring and E72's seven re-draws of it, and the
loader's hemibrain KC->KC counts and KC types; no olfsysm change; nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part check
    ./venv/bin/python -u experiments/E132_kckc_excitation.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part late
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part runs --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part runs --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E132_kckc_excitation.py --part score && ./venv/bin/python -u experiments/E132_kckc_excitation.py --part variants
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E132_kckc_excitation'
SCRIPT = 'experiments/E132_kckc_excitation.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
POPS = ('S4', 'T4')
POP_WIRINGS = {'S4': S_WIRINGS, 'T4': T_WIRINGS, 'all8': WIRINGS}
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
CONDS = ('on', 'off')
# arm -> step, matrix dose d (absolute KC-KC scale; None: the fit's own s), pooled neuron, hold
ARMS: Dict[str, dict] = {
    'U0': dict(step=1, dose=0.0, pooled=True, hold='pooled'),
    'UI2': dict(step=1, dose=2.0, pooled=True, hold='pooled'),
    'UI5': dict(step=1, dose=5.0, pooled=True, hold='pooled'),
    'UI10': dict(step=1, dose=10.0, pooled=True, hold='pooled'),
    'UX2': dict(step=2, dose=-2.0, pooled=True, hold='pooled'),
    'UX5': dict(step=2, dose=-5.0, pooled=True, hold='pooled'),
    'UX10': dict(step=2, dose=-10.0, pooled=True, hold='pooled'),
    'X2': dict(step=2, dose=-2.0, pooled=False, hold='threshold'),
    'X5': dict(step=2, dose=-5.0, pooled=False, hold='threshold'),
    'X10': dict(step=2, dose=-10.0, pooled=False, hold='threshold'),
}
STEP1_ARMS = ('U0', 'UI2', 'UI5', 'UI10')
STEP2_ARMS = ('UX2', 'UX5', 'UX10', 'X2', 'X5', 'X10')
MIRROR = {'UX2': 'UI2', 'UX5': 'UI5', 'UX10': 'UI10'}
DOSE0 = {'U': 'U0', 'X': 'E110_A'}
ARM_LABELS = {
    'U0': "the pooled neuron alone in the matrix's place (matrix 0), pooled hold",
    'UI2': 'the matrix at +2 beside the pooled neuron, pooled hold',
    'UI5': 'the matrix at +5 beside the pooled neuron, pooled hold',
    'UI10': 'the matrix at +10 beside the pooled neuron, pooled hold',
    'UX2': 'the matrix at -2 (excitation) beside the pooled neuron, pooled hold',
    'UX5': 'the matrix at -5 (excitation) beside the pooled neuron, pooled hold',
    'UX10': 'the matrix at -10 (excitation) beside the pooled neuron, pooled hold',
    'X2': 'the matrix at -2 (excitation) alone, threshold hold',
    'X5': 'the matrix at -5 (excitation) alone, threshold hold',
    'X10': 'the matrix at -10 (excitation) alone, threshold hold',
}
# references read in step 0 and score: kind -> what it is
REFS = {
    'family': "the family's saved fits (matrix at s)",
    'H50': "the matrix at s / 2, both rates re-held (E131's H50)",
    'E110_D': "the matrix at s / 2 held, the rest refitted (E110's part D; S's wirings)",
    'E110_A': 'KC-KC 0, both rates re-held (E110 part A)',
    'E110_B': "KC-KC 0, the rest refitted (E110's part B, counted start)",
    'E110_C': "the return fits from KC-KC 2 (E110's part C; real, P1, P3, M2)",
    'U': "the uniform neuron at the matrix's whole weight, both rates re-held (E122 step 1 on S, E123 on T)",
    'E129_P': "E129's pooled arm: the family's APL-on run, the pooled neuron in the matrix's place with APL off",
    'PO50': "E125's pooled survivor at 0.5 (PO50)",
}
PULSE_S = 0.5
STIM_END = 0.5
LANDING_DOSES = (2.0, 5.0, 10.0)
SCALE_BRACKET, SCALE_STEP = (1 / 50, 50.0), 1.25
RATE_TOL = 0.002
WITHIN = 0.003
N_WIRINGS_NEEDED = 6
MEAN_TOL = 0.02
RAW_LOSS_TOL = 0.005
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7
STOP_SLACK = 0.005
MAX_COUNT_FACTOR = 2.0
DOSE0_WIRINGS = 3
REALISM_THR_HZ = 60.0
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS + ('wU_scale',)
_CACHE: dict = {}


class SignedDoseError(NotImplementedError):
    """A negative KC-KC dose before olfsysm has the option for negative KC-KC weights (held until Yang's OK)."""


def require_nonnegative(d: float, b: float = 0.0) -> None:
    """Refuses a negative dose or pooled scale (SignedDoseError / ValueError)."""
    raise NotImplementedError


def uniform_neuron(ref: np.ndarray) -> np.ndarray:
    """E122's U: every entry ref's total / N^2 (the diagonal included)."""
    raise NotImplementedError


def combined_matrix(HB: np.ndarray, U: np.ndarray, d: float, b: float):
    """d HB + b U as a CSC matrix (the route: the simulator's KC-KC matrix with 'wKCKC_scale' 1)."""
    raise NotImplementedError


class Runner:
    """One wiring's family fit as a simulator whose KC-KC matrix is d HB + b U ('wKCKC_scale' d, 'wU_scale' b)."""

    def __init__(self, w: str):
        raise NotImplementedError

    def __call__(self, params: Mapping[str, float], spike_times: bool = False):
        raise NotImplementedError


def pooled_hold(runner: Runner, d: float, label: str) -> dict:
    """Reading 1: the pooled hold of one point."""
    raise NotImplementedError


def threshold_hold(runner: Runner, d: float, label: str) -> dict:
    """Reading 2: the threshold hold of one point (step 2's arms; refuses a negative dose)."""
    raise NotImplementedError


def held_ok(n: int, target: int) -> bool:
    """Reading 3: the response count within RATE_TOL of the target."""
    raise NotImplementedError


def late_share(spikes: pd.DataFrame, stim_end: float = STIM_END) -> float:
    """Reading 10: the share of spikes at or after the pulse's end."""
    raise NotImplementedError


def landing(counts: pd.DataFrame, responses: pd.DataFrame, M: np.ndarray, scale: float) -> pd.DataFrame:
    """Reading 7: per pair of odors, the directed ratios and the pair's ratio."""
    raise NotImplementedError


def stability(HB: np.ndarray, types: np.ndarray) -> pd.DataFrame:
    """Reading 8: the Perron root and the uniform share of the leading eigenvectors of HB and of each type's block."""
    raise NotImplementedError


def slopes(rd: pd.DataFrame, scales: Mapping[str, float]) -> pd.DataFrame:
    """Reading 9: the change per unit of scale from the family to H50 to E110_A, and carried through zero."""
    raise NotImplementedError


def references() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """The references' outputs and runs rows (REFS)."""
    raise NotImplementedError


def part_check(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> dict:
    raise NotImplementedError


def part_step0() -> dict:
    raise NotImplementedError


def part_late(wirings: Sequence[str] = WIRINGS) -> pd.DataFrame:
    raise NotImplementedError


def part_runs(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
              arms: Optional[Sequence[str]] = None) -> pd.DataFrame:
    raise NotImplementedError


def part_step2(*args, **kwargs):
    """Step 2 (UX2, UX5, UX10 with the pooled neuron; X2, X5, X10 alone): needs the olfsysm option."""
    raise SignedDoseError('step 2 needs the olfsysm option for negative KC-KC weights (held until Yang\'s OK)')


def part_step3(*args, **kwargs):
    """Step 3 (the shuffle control at a passing negative dose; the male calyces beside): needs the olfsysm option."""
    raise SignedDoseError('step 3 needs the olfsysm option for negative KC-KC weights (held until Yang\'s OK)')


def part_step4(*args, **kwargs):
    """Step 4 (the refit with a signed KC-KC scale in E10_joint_fit.py, the read on T, the sign control): needs the
    option and E10_joint_fit.py's signed scale, written only if the gate opens."""
    raise SignedDoseError('step 4 needs the olfsysm option and a signed KC-KC scale in E10_joint_fit.py')


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, flies_mean: float) -> dict:
    """Reading 4."""
    raise NotImplementedError


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    """Reading 5 (E129's `differ`)."""
    raise NotImplementedError


def stop_rule(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """Reading 6."""
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
