#!/usr/bin/env python3
"""E131: MB-C1's two cells beside the inhibition among KCs: a γ-directed feedback that the block spares, at its own
scale next to the KC-KC matrix, in half of the matrix's place, and with a time constant and a threshold of its own.

The entry is item E131 of experiment_priorities.md (proposed by viki 2026-11-23 16:10, moved into the list under the
idle-host rule; planned by a Fable subagent, claude-fable-5-2, spawned from viki 15:53). "The family" is C2's breadth-2
APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S is real, P1, P2, M1 and T
is P3, P4, M2, M3. "KC-KC at zero" is E110's part A. MB-C1 is olfsysm fa3543d's option: the hemibrain's two MB-C1_R
cells as feedback cells of their own (`fit_mb_model(mbc1_inh=, mbc1_tau=, mbc1_thr=, mbc1_pooled=)`, the simulator's
'wMBC1_scale'), alongside the KC-KC matrix, each with its own scale, present with APL on and with APL off.

The protocol (design item 1; nothing fitted). Each run is one of E83's simulators of the wiring at its fit, built with
the arm's MB-C1 options (`build_runner`), both response rates re-held to the family's by E77_input_screen.hold_rates
(the threshold bisected with APL off, then the tied APL scales with APL on). The KC-KC scale ('wKCKC_scale') and
MB-C1's ('wMBC1_scale') are simulator parameters, so one simulator serves every arm with the plain cells.

The arms (tags S1_<arm>_<wiring>; doses m are multiples of the wiring's fitted KC-KC scale s):
- B (24 runs), beside the matrix: the matrix at s, MB-C1 at m 0.25, 1 and 4 (B25, B100, B400).
- R (16 runs), in half of the matrix's place: the matrix at s / 2 with MB-C1 at m 0.5 (R50), and the matrix at s / 2
  alone as the reference (H50).
- X (8 runs), E125's split with a named survivor: the matrix at s / 2 with APL on and at 0 with APL off, MB-C1 at m 0.5
  in both (X50); read against E125's pooled neuron (PO50) and neuron per type (TY50) at the same d.
- The cell's own properties (16 runs, 8 more if gated), beside the matrix at m 1: late (L100: MB-C1's time constant
  1 s, its scale step 0's), strong drive only (D100: MB-C1's threshold step 0's), pooled (P100: the two cells as one,
  only if step 0's gate opens).
- References, read and not run: the family; KC-KC at zero (E110's part A, 'E110_A'); E125's PO50 and TY50.
- Not written: step 2 (the fit) and the permuted control; they run only if a gate opens (`--part score` reports it),
  and their code is written then.

Readings the entry leaves open, fixed here before any run:
1. The dose. m is MB-C1's scale as a multiple of the wiring's fitted KC-KC scale: wMBC1_scale = m x s (`mbc1_scale`),
   also in R and X (R50's MB-C1 at 0.5 s beside the matrix at 0.5 s). The loader's weights (pyolfsysm's
   `connectome_MBC1_cell_weights`, normalize=True) give each cell input weights summing to 1 and an equivalent KC x KC
   matrix (cells' output x input, summed over the cells) with mean row sum 1, as the hemibrain matrix has; so for
   uniform KC activity a the cells' total output is m s a N, the matrix's s a N, and m 1 delivers the same total. Step 0
   checks this on the model's KCs (`dose_check`: both totals N, each cell's input weights summing to 1, the equivalent
   matrix equal to `connectome_MBC1_KCKC_weights(pooled=False)`).
2. H50 is the plain cells' simulator at wMBC1_scale 0, which check (ii) shows equals the option off exactly.
3. X50's APL-off condition has wKCKC_scale 0 (hold_rates' off_extra, as E123's HBab), MB-C1 at 0.5 s in both.
4. A cell's activity (step 0, `cell_activity`): olfsysm's recursion m(t) = (1 - dt / tau) m(t - 1) + w_in @ spikes(t - 1)
   / tau, applied to the family's saved spike steps (E121 step 0's spike times, checked to add up to the saved counts),
   in Hz. Per cell and odor: its mean over the pulse (the steps with the odor on, 0 to 0.5 s) and its peak over the
   trial.
5. D100's threshold, per wiring: the largest activity of either cell over the 16 odors and every time step with APL on
   in the family's run (tau 50 ms), so that the cells are silent wherever the family's APL is intact. Check (iv) runs
   APL on at the family's parameters with MB-C1 at m 1 and this threshold and must reproduce the family's APL-on run
   exactly, else D100 is refused.
6. L100's scale, per wiring: s times the ratio of the cells' summed output over the pulse with APL off in the family's
   run at 50 ms to the same at 1 s (each cell's pulse mean weighted by its total output weight, summed over the cells
   and odors; `late_ratio`), so that the late cells' mean output over the pulse with APL off equals the plain cells'.
7. The pooled arm's gate: the Pearson r over the 16 odors of the two cells' pulse means, per wiring and condition (APL
   on and off); P100 runs if any of the 16 values is below POOLED_GATE_R (0.95).
8. A corner: the threshold or the APL scale at a bound of its bisection, or either response count more than RATE_TOL
   (0.2 %) from its target (E129's `held_ok`). An arm with corners on more than two wirings gets no class.
9. The classes (calibrated pattern part, S's and T's populations apart, against all recordings and each half A to D):
   - "below the family": lower than the family by WITHIN (0.003) or more on both populations and on 6 of 8 wirings, on
     every target; the mean APL effect under the readout within MEAN_TOL (0.02) of the flies' (E108's 'all' target,
     +0.178) on both populations; the raw-count loss against all recordings not higher than the family's by more than
     RAW_LOSS_TOL (0.005) on both populations; every run in the four bands;
   - "equal": E123's (`E123.classify`, one-sided);
   - "costs": the pattern part against all recordings higher than the family's by more than WITHIN on both
     populations;
   - "neither" otherwise; in that order of precedence; a corner arm "corner (no class)".
   - Beside the class: "keeps the gain" (E123's flag, read for R50 and X50), "stands in" (R50: E123's equal flag and
     R50 lower than H50 by E126's "differ" rule, E129's `differ`).
10. The outcomes are flags for the write-up: (a) every B dose costs, R50 does not stand in, X50 does not keep the gain;
    (b) B25 equal and every B dose above the largest tolerated one (the largest dose with every dose up to it equal)
    costs; (c) R50 stands in or X50 keeps the gain; (d) a B dose below the family (the fit's half, on T, is step 2's);
    (e) L100, D100 or P100 against B100 by the "differ" rule; (f) a corner, or a run outside the bands. (f) is a flag:
    bands enter the classes only where the classes name them ("equal", "below"), so a dose that leaves a band can
    still cost.
11. Step 2's gate (reported, not run): a B dose "below the family", R50 "stands in" or X50 "keeps the gain", none a
    corner.
12. Realism of MB-C1 per run (`cell_realism`, from the held point simulated again with its spike times and APL's
    trace, which must reproduce the held counts): each cell's activity in Hz (pulse mean over the odors, largest peak)
    with APL on and off; MB-C1's weight per output synapse against a KC->KC synapse's; and with APL on, the share of
    the inhibition delivered to γ KCs over the pulse that is MB-C1's against MB-C1's and APL's together (APL's: each
    KC's APL->claw weights times APL's trace), beside MB-C1's share of the inhibitory synapses onto γ against APL's
    (step 0's tables).

Parts (CLI `--part`):
- check (simulations, `--wirings`): on every wiring, at the family's saved parameters, (i) the option off and (ii) the
  option on at scale 0 reproduce the saved family run (APL on and off) with 0 differing spike counts; on 'real', (iii)
  the option at m 1 changes the counts; once step 0 is done, (iv) D100's threshold check. check_<host>_<lane>.json;
  the runs refuse until every wiring has passed (i) and (ii).
- step0 (no simulation): step0_synapses.csv (each cell's synapses in and out by KC type, with APL's and the KC->KC
  synapses beside), step0_activity.csv (each cell's activity per odor, wiring and condition), step0_gamma_input.csv
  (each γ KC's MB-C1 input at m 1 beside its KC-KC input), step0.json (the two fixed numbers per wiring, the cells'
  agreement and the pooled gate, the dose arithmetic, the dose check), provenance_step0.json.
- runs (simulations, resumable, per lane `--wirings`, `--arms`; a point claimed by mkdir of
  results/E131_mbc1_cells/claims/<tag>/, skipped when done or claimed, `--reclaim` to take over a dead process's): rows
  in runs_<host>_<lane>.csv, provenance (`H.model_provenance` per point) in provenance_<host>_<lane>.json, runs saved
  locally (outputs/E131_mbc1_cells/step1/) and in the model-output archive (data/model_outputs/E131_mbc1_cells/step1/).
- score (no simulation): readings.csv, runs.csv, verdict.json (classes, comparisons, outcomes, step 2's gate),
  comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- plan: the points a lane would run.

The holds (CLAUDE.md, "Connectome data"): only the model's hemibrain wiring and E72's seven re-draws of it, the
loader's hemibrain tables (KC->KC, APL<->KC, KC types) and Pratyush's hemibrain MB-C1 tables as pyolfsysm reads them;
the eight wirings keep the hemibrain's 1732 KCs, so MB-C1's weights are the same on each. No male CNS table, no olfsysm
change; nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part check
    ./venv/bin/python -u experiments/E131_mbc1_cells.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part runs --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part runs --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E131_mbc1_cells.py --part score && ./venv/bin/python -u experiments/E131_mbc1_cells.py --part variants
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E131_mbc1_cells'
SCRIPT = 'experiments/E131_mbc1_cells.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
POPS = ('S4', 'T4')
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
# E121 step 0 names the hemibrain wiring 'S'
E121_LABEL = {'real': 'S'}
# arm -> (group, KC-KC scale with APL on / off as multiples of s, MB-C1 dose m, MB-C1 form)
ARMS: Dict[str, dict] = {
    'B25': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=0.25, form='plain'),
    'B100': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=1.0, form='plain'),
    'B400': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=4.0, form='plain'),
    'R50': dict(group='R', kckc_on=0.5, kckc_off=0.5, m=0.5, form='plain'),
    'H50': dict(group='R', kckc_on=0.5, kckc_off=0.5, m=0.0, form='plain'),
    'X50': dict(group='X', kckc_on=0.5, kckc_off=0.0, m=0.5, form='plain'),
    'L100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='late'),
    'D100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='strong'),
    'P100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='pooled'),
}
GROUPS = ('B', 'R', 'X', 'forms')
B_ARMS = ('B25', 'B100', 'B400')
FORM_ARMS = ('L100', 'D100', 'P100')
E125_READ = ('PO50', 'TY50')
REFS = ('family', 'E110_A') + E125_READ
LATE_TAU_S = 1.0
PLAIN_TAU_S = 0.05
PULSE_S = 0.5
POOLED_GATE_R = 0.95
RATE_TOL = 0.002
WITHIN = 0.003
N_WIRINGS_NEEDED = 6
MEAN_TOL = 0.02
RAW_LOSS_TOL = 0.005
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7
REALISM_THR_HZ = 60.0
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS + ('wMBC1_scale',)
# (b, c) of the "differ" comparisons
COMPARISONS = (('R50', 'H50'), ('X50', 'PO50'), ('X50', 'TY50'), ('L100', 'B100'), ('D100', 'B100'), ('P100', 'B100'))
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E131_mbc1_cells/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    raise NotImplementedError


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def tag_of(arm: str, w: str) -> str:
    return f'S1_{arm}_{w}'


def mbc1_scale(arm: str, s: float, late_ratio: Optional[float] = None) -> float:
    """Reading 1: MB-C1's 'wMBC1_scale' of an arm on a wiring with fitted KC-KC scale `s` (L100: times `late_ratio`)."""
    raise NotImplementedError


def arm_params(arm: str, p: Mapping[str, float], late_ratio: Optional[float] = None) -> Tuple[Dict[str, float], dict]:
    """(the arm's starting parameters: the fit's with wKCKC_scale kckc_on x s and wMBC1_scale; hold_rates' off_extra:
    wKCKC_scale kckc_off x s where it differs from APL on, else {})."""
    raise NotImplementedError


def build_kws(arm: str, w: str, step0: Optional[Mapping] = None) -> dict:
    """The MB-C1 kwargs of an arm's simulator (`make_fixed_param_kc_simulator` model kwargs): mbc1_inh (overridden by
    'wMBC1_scale'), and mbc1_tau (L100), mbc1_thr (D100, step 0's per wiring) or mbc1_pooled (P100)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 0 (no simulation)

def cell_weights(kc_index: pd.Index, pooled: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """pyolfsysm's `connectome_MBC1_cell_weights` for the model's KCs (normalize=True), as fit_mb_model sets them."""
    raise NotImplementedError


def dose_check(kc_index: pd.Index) -> dict:
    """Reading 1's check on the model's KCs: the equivalent matrix's total and the hemibrain matrix's (both N), each
    cell's input weights' sum (1), and the equivalent matrix against `connectome_MBC1_KCKC_weights(pooled=False)`."""
    raise NotImplementedError


def synapse_table(kc_index: pd.Index) -> pd.DataFrame:
    """Each MB-C1 cell's synapses in and out by KC type on the model's KCs, with APL's calyx synapses onto them and the
    KC->KC synapses among them beside, per KC type (counts and per KC)."""
    raise NotImplementedError


def dose_arithmetic(syn: pd.DataFrame) -> dict:
    """The dose at which an MB-C1 output synapse weighs what a KC->KC synapse does (average and per cell), and MB-C1's
    share of the inhibitory synapses onto γ against APL's."""
    raise NotImplementedError


def spike_times(w: str, cond: str) -> pd.DataFrame:
    """E121 step 0's spike times of a wiring and condition (outputs/E121_onset_race/step0/), checked to add up to the
    family's saved spike counts."""
    raise NotImplementedError


def cell_activity(spikes: pd.DataFrame, w_in: np.ndarray, n_odors: int, tau: float, dt: float, n_steps: int,
                  pulse_steps: Tuple[int, int]) -> Dict[str, np.ndarray]:
    """Reading 4: each cell's activity (cells x odors): 'pulse_mean' and 'peak', in Hz."""
    raise NotImplementedError


def late_ratio(act50: Mapping[str, np.ndarray], act_late: Mapping[str, np.ndarray], out_w: np.ndarray) -> float:
    """Reading 6: the cells' summed output over the pulse at 50 ms over the same at 1 s (APL off)."""
    raise NotImplementedError


def cells_agreement(pulse_mean: np.ndarray) -> float:
    """Reading 7: the Pearson r over the odors of the two cells' pulse means (nan for one cell)."""
    raise NotImplementedError


def pooled_gate(agreement: Mapping[Tuple[str, str], float]) -> dict:
    """Reading 7: open if any (wiring, condition) agreement is below POOLED_GATE_R."""
    raise NotImplementedError


def gamma_input(counts: pd.DataFrame, s: float) -> pd.DataFrame:
    """Each γ KC's MB-C1 input at m 1 beside its KC-KC input (trial integrals, s x matrix @ counts), per odor."""
    raise NotImplementedError


def part_step0() -> dict:
    raise NotImplementedError


def step0_fixed() -> dict:
    """step0.json (raises if step 0 has not run)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the runner and the checks (simulations)

def build_runner(arm: str, w: str, params: Mapping[str, float], step0: Optional[Mapping] = None) -> Tuple[Callable, dict]:
    """(simulate, info): E83's simulator of wiring `w` at `params` with the arm's MB-C1 kwargs (None: the option off)."""
    raise NotImplementedError


def part_check(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> dict:
    raise NotImplementedError


def check_passed(need_threshold: bool = False) -> bool:
    raise NotImplementedError


def held_ok(n: int, target: int) -> bool:
    """Reading 8: the response count within RATE_TOL of the target."""
    raise NotImplementedError


def held(arm: str, w: str, label: str) -> dict:
    """One held point of an arm on a wiring (hold_rates on `build_runner`'s simulator)."""
    raise NotImplementedError


def cell_realism(sim: Callable, info: dict, res: Mapping, arm: str, w: str) -> dict:
    """Reading 12."""
    raise NotImplementedError


def planned(wirings: Sequence[str] = WIRINGS, arms: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """The (arm, wiring) points of a lane: group by group (B, R, X, the forms), each wiring's arms together."""
    raise NotImplementedError


def runnable_arms(step0: Optional[Mapping] = None, check: Optional[bool] = None) -> Tuple[str, ...]:
    """Every arm but P100 unless step 0's pooled gate is open, and D100 only once check (iv) has passed."""
    raise NotImplementedError


def part_runs(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
              arms: Optional[Sequence[str]] = None) -> pd.DataFrame:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# score

def flies_mean_effect() -> float:
    """The flies' mean APL effect against all recordings (E108's 'all' target)."""
    raise NotImplementedError


def corners(runs: pd.DataFrame, kind: str) -> dict:
    raise NotImplementedError


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, flies_mean: float) -> dict:
    """Reading 9."""
    raise NotImplementedError


def outcomes(classes: Mapping[str, dict], comps: Mapping[str, dict], runs: pd.DataFrame) -> dict:
    """Reading 10."""
    raise NotImplementedError


def step2_gate(classes: Mapping[str, dict]) -> dict:
    """Reading 11."""
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', 'step0', 'runs', 'score', 'variants', 'plan'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    ap.add_argument('--arms', nargs='+', default=None, choices=list(ARMS))
    ap.add_argument('--lane', default=None, help='the lane name in the tables (default the wirings joined by -)')
    ap.add_argument('--reclaim', action='store_true', help='take over claims left by a process that died')
    return ap


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
