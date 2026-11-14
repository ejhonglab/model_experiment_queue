#!/usr/bin/env python3
"""E95: what the best family's APL is for a KC: its inhibition at the PN boutons, and then all of it, replaced by a static
inhibition per KC that is present only with APL on (each KC's own loss of drive, or its number of claws alone); the same
offset made at the KCs in the feedback APL with the bouton path at a tenth; and the family fitted again that way.

The entry is item E95 of experiment_priorities.md (proposed by viki 2026-11-14 06:01; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 05:18).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E94's `WiringModel94`: E83's `WiringModel` through E91's `WiringModel91`, with a dose on the APL->bouton scale). S: real,
P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's, C / D the second split's (E89's `targets89`).

Definitions, per wiring at its fit, from E90's archived stage runs (window means over the 0.5 s pulse) and the simulator's
claw x bouton weights (E94's `wiring_W`, `kc_drive`): L^b(k, o) the loss of KC k's summed claw drive made at the boutons
for odor o (the drive from the boutons without APL minus that with APL), L^c(k, o) the loss made at the claws (the drive
from the boutons with APL minus the claws' states), b_k and c_k their means over the 16 odors, n_k the KC's number of claws
(`kc_losses`). A static inhibition is olfsysm's `ext_inh` (subtracted from a KC's input while the odor is on; KCs by
odors), set for the APL-on run only (`set_ext_inh_ko`, which maps the outputs' odor columns to olfsysm's odor order and
refuses a negative or non-finite value in Python: olfsysm's own check is inside its parallel region and aborts the
process). The APL-off run of every point is the saved fit's. The KC-KC scale and the somatic threshold stay at the fit's.

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B, C and D refuse to run until it passes): on every
  wiring, a static inhibition of zero (a KCs x odors matrix of zeros, set) at the fit's dose and weights reproduces the
  saved APL-on run, and APL removed (every APL scale 0) with a static inhibition of zero the saved APL-off run, with 0
  count differences; on the real wiring the KC-by-odor setter with the same value for every odor is E91's per-KC setter
  (the same olfsysm matrix and the same counts), a matrix non-zero for one odor changes only that odor's column, and part
  D's weights with no unweighted component at dose 0.1 reproduce E94's saved point B_sub_k2_d0.1 (its factor and its
  counts). Each wiring's per-KC thresholds (fixed_thr + 2 x spont_in), claws, connectome and family APL->KC weights are
  cached under outputs/E95_static_apl/W/ for part A.
- A (no simulations; part 0's cache and E94's wiring cache): A1 the split of L^b and L^c into a value per KC, per odor and
  the rest, by claw class and over the KCs that respond without APL (a1_split.csv); A2 b_k against n_k (the slope a claw,
  the scatter, the log-log slope, what carries the rest: the APL->bouton weights of the KC's boutons, its glomeruli
  against E75's map, its type) and c_k against the breadth factor and n_k (a2_per_kc.csv, a2_glomeruli.csv); A3 b_k, c_k
  and their sum against each KC's threshold, and the responses without APL a threshold raise of b_k, or b_k + c_k, would
  remove by the window-mean drive, by odor group and the co-responses by group of pairs (a3_sizes.csv, a3_removed.csv; a
  count-level description only); A4 APL's activity by odor and before the odor in the fits and in E94's saved runs, the
  connectome's per-claw APL->KC weights on every wiring (a4_apl.csv, a4_weights.csv); A5 where E87's and E88's static
  forms, E91's static control and E80's threshold difference sit (a5_archive.csv). a_summary.json.
- B, C, D: the point lists (`--arm`, `--wirings`, resumable; runs_B.csv, runs_C.csv, runs_D.csv; runs under
  outputs/E95_static_apl/runs/ with their stage matrices and each static inhibition saved per KC):
  - B, the bouton path as a static inhibition (dose 0; the claw path's feedback kept with the family's weights; the tied
    APL scale bisected in APL_BRACKET): B1 b_k; B2 beta n_k (beta the wiring's least-squares slope through the origin);
    B3 n_k^e for e in B3_EXPONENTS scaled to b's mean; B4 b_k permuted among the KCs with the same number of claws and
    over all KCs (seed PERM_SEED); B5 L^b(k, o).
  - C, the whole APL as a static inhibition (every APL scale 0; the static inhibition times a factor bisected in
    STATIC_BRACKET): C1 b_k + c_k; C2 L^b + L^c; C3 beta' n_k (beta' the slope of b + c); C4 beta n_k + c_k; C5 C1's
    values permuted over all KCs; C6 the claw path alone static (dose 1, the APL->KC weights 0, c_k, the tied scale
    bisected in APL_BRACKET).
  - D, the offset at the KCs in the feedback APL (dose D_DOSE; the APL->KC weight of a claw the family's plus an
    unweighted component, d in units of what the bouton path took, `component_unit`; the tied scale bisected): D1 the
    same weight on every claw, D2 the connectome's per-claw weights unweighted by breadth, both at D_SIZES; D3 D1 at
    d = 1 with no bouton path. D0 (dose 0.1, no component: E94's B_sub_k2_d0.1) is the control of part D's settings: on S
    E94's saved runs, on T run here (part E).
- score (`--no-co`, `--no-zcell`): readings.csv per variant and unit (each wiring, S and T where all four have a run) and
  readout ('cal' the calibrated readout, K set on each run's APL-on counts, the reading of the choice and the gate;
  'calK' the fit's K; 'raw'; 'zcell' and 'zcell_lin' E73's per-KC reading under the saturating and the linear map) against
  every target of `targets89` (E92's `score_row92`: loss, pattern part with its amplitude and shape parts, spread ratio,
  Pearson r and ρ on the per-pair APL effect, the two groups' errors, the mean APL effect); pairs.csv, line.csv (E92's
  part A line), coresponse.csv (E89's readings of the alcohol pairs' co-response change); then part E's tables.
- E: the choice per half (`choose95`), T's runs of a chosen setting and its control (runs_E.csv), the gate (E92's
  `gate_direction92` / `gate92`, E91's rules), the equal rule in three kinds (`equal_rule95`), verdict.json.
- F: the family fitted again at the connectome's bouton share with the unweighted component (`E10_joint_fit.py -k 2
  --aplpn-hold F --aplkc-flat uniform ...`, `f_command`), only on its condition (`f_condition`); `--part Fscore` scores
  them by E94's D1 rule (`f_scores`).
- G: documented stubs (`part_g`): run only by the entry's rules.
- steps (`--variants`): the step from the claws to the somata in wild-type flies for chosen or equal settings (E94's
  practice).
- variants: variants.csv and model_pair_distances_long.csv for the archive.
- table: comparison_table.md and the figure.

Deviations from the entry, forced by the code: (filled in as the parts are written)

Outputs: results/E95_static_apl/; spike counts, responses, stage matrices and static inhibitions under
outputs/E95_static_apl/ (runs/, W/, F/).

Usage (from the repository root; one 16-thread process at a time on viki; check `uptime` first):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E95_static_apl --part 0
    ./venv/bin/python -u -m experiments.E95_static_apl --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E95_static_apl --part B --arm B1
    ./venv/bin/python -u -m experiments.E95_static_apl --part score
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E89_apl_breadth_weighting as E89
from experiments import E91_summed_drive_inhibition as E91
from experiments import E92_gain as E92G
from experiments import E94_bouton_share as E94

EXPERIMENT = 'E95_static_apl'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

APL_BRACKET = E92G.APL_BRACKET                      # (0.01, 10): the tied APL scale (B, C6, D), E92's bracket
STATIC_BRACKET = (0.05, 20.0)                       # C: the static inhibition's factor (the probes: 1.24 to 1.48)
APL_STEP = E83.APL_STEP
RAW_CRIT = E92G.RAW_CRIT                            # a model response is two or more spikes
CLAW_GROUPS = E94.CLAW_GROUPS                       # 1-4, 5-8, 9+ claws
B3_EXPONENTS = (0.0, 0.5, 1.5, 2.0)                 # B3: n_k^e scaled to b's mean (1 is B2)
PERM_SEED = 1                                       # B4, C5 and the gate's controls: one seed
D_DOSE = 0.1                                        # D: E94's anatomical dose
D_SIZES = (0.5, 1.0, 2.0)                           # D1, D2: the component's size in units of the bouton path
D3_SIZE = 1.0
COMPONENTS = ('uniform', 'connectome')
CONNECTOME_SHARE_RANGE = (0.096, 0.118)             # the connectome's share of the claws' loss made at the boutons
E94_D0_TAG = 'B_sub_k2_d0.1'                        # E94's point at dose 0.1 (D's control; part 0's identity)

# the variants of each arm (88 points on S: B 36, C 24, D 28)
ARMS: Dict[str, Tuple[str, ...]] = {
    'B1': ('B1_bk',),
    'B2': ('B2_nk',),
    'B3': tuple(f'B3_nk{e:g}' for e in B3_EXPONENTS),
    'B4': (f'B4_bk_permclaw{PERM_SEED}', f'B4_bk_permall{PERM_SEED}'),
    'B5': ('B5_bko',),
    'C1': ('C1_bck',),
    'C2': ('C2_bcko',),
    'C3': ('C3_nk',),
    'C4': ('C4_nkck',),
    'C5': (f'C5_bck_permall{PERM_SEED}',),
    'C6': ('C6_ck',),
    'D1': tuple(f'D1_uniform_d{d:g}' for d in D_SIZES),
    'D2': tuple(f'D2_connectome_d{d:g}' for d in D_SIZES),
    'D3': (f'D3_uniform_d{D3_SIZE:g}_dose0',),
}
PART_OF_ARM = {a: a[0] for a in ARMS}
REF = 'ref'                                         # the unchanged fits (source 'saved')
D0 = 'D0'                                           # dose 0.1 without the component (E94's point): the control of D

# part E (the entry's rules: E91's, as in E94)
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
CONTROL_RATIO = E91.CONTROL_RATIO
HALVES = E91.HALVES
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
EQUAL_KINDS = {'bouton_static': ('B1', 'B2', 'B5'), 'whole_static': ('C1', 'C2', 'C4', 'C6'), 'at_kcs': ('D1', 'D2', 'D3')}
CHOICE_ARMS = ('B1', 'B2', 'B3', 'B5', 'C1', 'C2', 'C3', 'C4', 'C6', 'D1', 'D2', 'D3')   # B4 and C5 excepted

# part F
F_COND = 0.008                                      # S's best of D1 or D3 within this of the unchanged fits' 0.0500
F_FLAT_BOUNDS = (0.0, 8.0)                          # the component's size, a free parameter
F_MAXFEV = 200
F_TAG = '_f1'
F_START_FILE = 'f_start_params.json'
E94_REFIT_S = 0.0593                                # E94's refit at the same bouton share without the component (S)
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
FLIES_MEAN = E89.FLIES_MEAN

_CACHE: Dict[tuple, object] = {}


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# names

def parse_variant(v: str) -> dict:
    """A variant's fields: arm, part ('B' / 'C' / 'D' / 'ref'), mode ('bouton' the bouton path static, 'whole' every APL
    scale 0, 'claw' C6, 'kcs' part D, 'ref'), static (the static inhibition's form: 'b', 'nk', 'nk^e', 'b_o', 'bc',
    'bc_o', 'nk_bc', 'nk_c', 'c', or None), exponent, permute ('claw' / 'all' / None) and seed, component, size and dose
    (part D), control and control_of, screen (the arm, 'ref' or 'control') and the kind of equal it can be."""
    raise NotImplementedError


def control_of(variant: str) -> str:
    """The gate's control of a setting: a static setting's values permuted over the KCs ('<v>_ctrl', seed PERM_SEED);
    part D's settings: D0 (the same setting without the component)."""
    raise NotImplementedError


def arm_variants(arm: str) -> List[str]:
    return list(ARMS[arm])


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


# ---------------------------------------------------------------------------
# the losses per KC, the static matrices, the component

def kc_losses(wiring: str, kc_index: Optional[pd.Index] = None) -> dict:
    """L^b, L^c (KCs x odors, the outputs' row and column order), b, c (per KC), n (claws per KC), D_off (the summed
    claw drive without APL, KCs x odors) of a wiring's fit, from E90's stage runs and E94's wiring cache."""
    raise NotImplementedError


def slope_through_origin(x: np.ndarray, y: np.ndarray) -> float:
    """The least-squares slope of y on x through the origin."""
    raise NotImplementedError


def permute_values(v: np.ndarray, seed: int, groups: Optional[np.ndarray] = None) -> np.ndarray:
    """`v`'s rows permuted over all rows, or within each group of `groups` (e.g. the KCs' claw counts)."""
    raise NotImplementedError


def static_matrix(pv: dict, losses: dict) -> pd.DataFrame:
    """The static inhibition (KCs x odors, the outputs' order, not negative) of a parsed variant from `kc_losses`."""
    raise NotImplementedError


def component_unit(bclaw_mean: float, apl_scale: float, apl_mean: float, w_mean: float = 1.0) -> float:
    """The unscaled APL->KC weight, per unit of a component's mean weight `w_mean`, that at the fit's APL->KC scale and
    APL's mean activity takes from a claw the fit's mean loss a claw at the boutons (`bclaw_mean`): d = 1."""
    raise NotImplementedError


def d_weights(w_fam: pd.Series, wnorm: pd.Series, kind: str, d: float, unit: float) -> pd.Series:
    """Part D's APL->KC weights per claw: the family's plus d x unit x the component ('uniform': the connectome's mean
    weight on every claw; 'connectome': each claw's connectome weight, unweighted by breadth)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the simulator: the static inhibition KCs x odors

def set_ext_inh_ko(mp, rv, M: Optional[np.ndarray], order: Optional[np.ndarray] = None, tol: float = 1e-9) -> dict:
    """Set olfsysm's static inhibition from a KCs x odors matrix in the outputs' column order (`order[j]` olfsysm's odor
    of output column j), clipped at 0 (values below -tol x max refused), checked finite, the right size, ext_inh_scale 1
    and n_claws_active_to_spike <= 0, all in Python; None removes it. Returns the clipping's record."""
    raise NotImplementedError


class StaticModel(E94.WiringModel94):
    """E94's `WiringModel94` (a dose on the APL->bouton scale; olfsysm's `mp`, `rv`, `order`) with the runs of this
    experiment: `run_on(x)` (the tied APL scales times x, the dose), `run_static(f)` (every APL scale 0, the static
    inhibition times f), each with the static inhibition set for that run only."""


def reading95(model: StaticModel, pv: dict, M: Optional[pd.DataFrame], wk: Optional[pd.Series], label: str = '') -> dict:
    """One setting on a wiring at its fit: the weights, the dose, the static inhibition; the bisection (the tied APL scale
    in APL_BRACKET, or the static factor in STATIC_BRACKET), the run again with APL's activity saved, its stage matrices;
    everything reset on return, also on an error."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# realism beside the four bands

def static_realism(M: Optional[pd.DataFrame], thr: pd.Series, n_claws: pd.Series, d_off: pd.DataFrame,
                   on_counts: Optional[pd.DataFrame] = None, off_counts: Optional[pd.DataFrame] = None) -> dict:
    """The static inhibition per KC (its mean over the odors) as a multiple of the KC's threshold (10th percentile, median,
    90th, largest; by claw class), the share of KCs whose static inhibition exceeds their largest drive without APL, the
    share silent by the window-mean drive, and the shares of KCs silent in the run with APL on and off."""
    raise NotImplementedError


def d_realism(W: pd.DataFrame, b_off: np.ndarray, b_on: np.ndarray, c_on: np.ndarray, w_fam: np.ndarray,
              w_comp: np.ndarray, kc_of: np.ndarray, responders: Optional[np.ndarray] = None) -> dict:
    """Part D's: the share of the claws' loss made at the boutons, by the unweighted component and by the weighted one
    (each claw's loss at the claws split in proportion to its weights), and the spread of the per-KC summed APL->KC
    weights over the KCs that respond without APL."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part 0 and part A

def part_0(out: Path, wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json; raises if one fails."""
    raise NotImplementedError


def identity_passed() -> bool:
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """A1 to A5 and a_summary.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# parts B, C, D

def run_point(model: StaticModel, v: str, part: str = '', cc: Optional[dict] = None) -> Tuple[dict, dict]:
    """One variant on a wiring at its fit, saved with its stage matrices and its static inhibition; returns (row,
    reading)."""
    raise NotImplementedError


def part_points(out: Path, part: str, arms: Sequence[str], wirings: Sequence[str]) -> pd.DataFrame:
    """The arms' points on `wirings` (resumable; runs_<part>.csv). Needs part 0 to have passed."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the readings, part E

def runs_table(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def readings95(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = ()) -> pd.DataFrame:
    raise NotImplementedError


def choose95(rd: pd.DataFrame, half: str, runs: Optional[pd.DataFrame] = None) -> dict:
    """The choice on S against half `half`: the lowest pattern part of S's population under the readout among CHOICE_ARMS'
    settings (inside the bands on all four wirings, not at a bracket's bound) and the unchanged fits; something to confirm
    if a setting gains CHOICE_MIN_GAIN or more over the unchanged fits."""
    raise NotImplementedError


def equal_rule95(rd: pd.DataFrame) -> List[dict]:
    """E94's rule in the three kinds of EQUAL_KINDS."""
    raise NotImplementedError


def part_e(out: Path) -> dict:
    raise NotImplementedError


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part F, part G

def f_condition(rd: pd.DataFrame) -> dict:
    """Whether part F runs: S's best setting of D1 or D3 within F_COND of UNCHANGED_S['all'] against all recordings;
    D2's (`--aplkc-flat connectome`) only if a D2 setting is an equal at fixed scales."""
    raise NotImplementedError


def f_command(label: str, component: str = 'uniform', maxfev: int = F_MAXFEV) -> List[str]:
    raise NotImplementedError


def part_f(out: Path, pop: str = 'S', wirings: Optional[Sequence[str]] = None, component: str = 'uniform',
           force: bool = False, dry_run: bool = False) -> List[str]:
    raise NotImplementedError


def f_scores(out: Path, write: bool = True) -> dict:
    raise NotImplementedError


def part_g(out: Path) -> None:
    """Part G (only by the entry's rules): on a pass by (E) or as a fitted family of eight, the firing-ceiling family with
    the response rate with APL on re-held, S at the flies' 2 s pulse with both response rates re-held, and S through the
    trial noise as E72's part E ran it; on an equal of C1, C1 on S at the flies' 2 s pulse with both response rates
    re-held (E80's part A). Not written: it runs only by those rules."""
    raise NotImplementedError('E95 part G runs only by the entry\'s rules (a pass, or an equal of C1); it is not written')


# ---------------------------------------------------------------------------
# the archive, the table

def write_variants(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def comparison_table(out: Path) -> str:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'score', 'E', 'F', 'Fscore', 'G', 'steps',
                                                      'variants', 'table'))
    ap.add_argument('--arm', default=None, help='B: B1..B5, C: C1..C6, D: D1..D3 (comma-separated; default all of the part)')
    ap.add_argument('--wirings', default=None, help="'S', 'T', 'all8' or a comma-separated list (0: all8; B / C / D: S)")
    a = ap.parse_args()
    raise NotImplementedError(f'E95 part {a.part}: not written yet')


if __name__ == '__main__':
    main()
