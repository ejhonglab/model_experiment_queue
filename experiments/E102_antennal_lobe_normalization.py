#!/usr/bin/env python3
"""E102: the antennal lobe's normalization at the strength the flies' APL asks for. The flies' APL falls with concentration
where the best family's holds its strength; in the family APL follows its claws, and the claws follow PN totals that the
model's antennal lobe makes nearly equal. What a weaker PN normalization (olfsysm's pn_inhadd raised, pn_inhsc with it)
costs or gains on the per-pair APL effect, at fixed scales on the eight wirings and with the family fitted again, and what
it does to the KCs' responses at low concentrations.

The entry is item E102 of experiment_priorities.md (proposed by atlas 2026-11-15 19:50; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-15 19:33; claimed by hal 19:53). C2's family (the breadth-2 APL-effect fit
of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit` /
`build_simulator`; S: real, P1, P2, M1; T: P3, P4, M2, M3). Halves A and B of the recordings are E77's, C and D E86's.
The input is the lab's in every verdict run.

Reused, not duplicated:
- E77 (experiments/E77_input_screen.py): `inhsc_for` (pn_inhsc raised with pn_inhadd as E3 / E32 did), `handle_breadth`
  (the breadth behind the APL->KC weights for a handle, here on each wiring's claw x bouton matrix), `hold_rates` (both
  response rates re-held), and the S1 points of its screen (part 0's second identity, part A3).
- E83 (experiments/E83_cross_wiring_fit.py): the wirings (`saved_fit`, `wiring_parts`, `build_simulator`, which takes the
  PN kwargs as `extra_kws` since this experiment), the scores of a population (`score_runs`, `effect_scores`,
  `pair_distances`, `realism`).
- E86: the targets with halves C and D (`targets`).
- E93 / E95: E73's per-KC reading under the saturating and the linear map (`E95_static_apl.zcell_distances_map`).
- E98 (experiments/E98_apl_response_strength.py): the extended input (`b1_input`, 52 stimuli with the blank), its saved
  layout (B1_LAYOUT; `save_b1`, `load_b1`), APL's and the PN stage's reading (`read_apl_pn`, `relabel`), the pre-odor
  levels (`with_a0`: the blank's window carries, the settled level is the partner, E90's level is shown and decides
  nothing; E98's amendment of 2026-11-15 17:59), the amended identity (`a0_identity`), the indicator family, `classify`,
  `comparison_verdict` and the verdict itself (`verdict_core`, read with the flies' intervals of E98's saved
  a3_ratios.csv and its saved pulse test and solvent rates: no recording is read again).
- E100: W1's six weak-odor ratios (`read_ratio`, `family_on`, `match_b1`, `estimate_labels`, `held_class`,
  `w1_verdict`), against the flies' intervals of its saved b_ratios.csv; the flies' KC response rates at low
  concentration (`kc_data_rates`, E72 part F's data_fractions.csv).
- E101: the claw stage on the extended input (E90's `read_stages`, `family_table`, `unit_meta`), X2's model side
  (`model_x2`, `x2_model_classes`) against the flies' intervals of its saved c_x2_flies.csv (pa over ms near 10^-3, and
  APL's ratio over the claws', Q).
- E90: the PN stage against the PN dendrites (`pn_statistics`: the totals' CV over the 12 megamat odors, the rms
  difference from the dendrites' group distances).

Definitions (the entry's, fixed before any reading):
- The normalization: olfsysm's PN gain is pn_inhsc / (pn_inhadd + 0.25 inhA + 0.75 inhB) inside the PN tanh. The ladder
  (LADDER): pn_inhadd 100, 314, 1000, 3141, 10000 and the gain fixed (pn_inhadd FIXED_GAIN_INHADD, 10^6), pn_inhsc
  raised with it (`E77_input_screen.inhsc_for`) so that the gain at a typical inhibition is kept. The identity point
  (STEP_IDENTITY) passes the defaults (31.4088, 368.6631) explicitly.
- Rates re-held (E77's screen; `hold_point`): the breadth weights recomputed for the step on each wiring's own wiring
  (`point_breadth`), the somatic threshold bisected to the wiring's fit's count of responses with APL off, then the
  tied APL scale (APL->KC with KC->APL and bouton->APL) to its count with APL on; the KC-KC and APL->bouton scales as
  fitted.
- Each point is also run at its held parameters on E98's extended input (52 stimuli and the blank), APL on and acutely
  off (the APL-off run at hold_rates' params_off, as E98's B4), with PN rates, APL and the claws' summed drive saved
  (`run_extended`).
- The model's APL is read as E98 reads it: against the blank's window in the same run (the carrying level), with the
  settled pre-odor level as partner; E90's pre-odor mean decides nothing (unresolved_issues.md, the pre-odor start-up).
- Every pattern part is given with its amplitude part, its shape part, the spread ratio, ρ and Pearson r on the per-pair
  APL effect, the loss and the mean APL effect, under the calibrated readout (K set on each run's APL-on counts) and on
  raw counts. No level is chosen on the main target or on a correlation.

The levels (part C; `levels`, fixed here):
- L_apl: the smallest pn_inhadd at which V2 (paa at 10^-5 over 10^-4) is inside by E98's rule (`verdict_core`: at least 6
  of the 8 wirings, the carrying reading and level, held under its partners).
- L_lin: the smallest at which V2's linear reading (n 1, evoked; carrying reading and level) lies in the flies' interval
  on at least LEVEL_MIN_WIRINGS (6) of the 8 wirings.
- L_free: the largest that is an equal at fixed scales on S by E94's rule (`equal_rule`).
- A step with a flagged point (any of its eight runs outside the realism bands) is not a level. Any of the three may not
  exist. If L_apl does not exist, the refit runs at the last unflagged step (`refit_levels`).

The classes (part D; `refit_class`, fixed here), against the unchanged fits on S (BASE_PATTERN, recomputed in part 0):
- an equal (E94's rule): the pattern part within EQUAL_PATTERN (0.003) against all recordings and against halves A and B,
  the four fits inside the bands, the raw-count loss within EQUAL_RAW (0.005), the mean APL effect within EQUAL_MEAN
  (0.01);
- better: lower by BETTER_MARGIN (0.003) or more against all recordings and lower against each of halves A to D;
- a price: anything else.
- On an equal or better, T's four wirings are refitted at that level; better stands only if T's population is also lower
  than its fits by BETTER_MARGIN or more against all recordings and lower against each half, at least
  T_BETTER_MIN_WIRINGS (6) of the 8 single wirings are lower and the raw-count loss is not higher (`t_better_holds`).

The KCs at low concentration (outcome (e); `low_conc_class`, fixed here): the model's response rate at ma 10^-7 (the
kiwi-control session's estimate, E49's source, carrying; the other sessions beside), the median over the eight wirings,
against the flies' (E72 part F; control flies with APL on, APL output blocked (TNT) with APL off). The distance from
the flies is |log(model / flies)|; at a step it is 'further' if it exceeds C2's (the identity point's) by a factor of
LOW_CONC_FACTOR (1.25) or more, 'nearer' if C2's exceeds it by that factor, 'about the same' otherwise.

Parts (CLI `--part`):
- 0, the identities (identity.json; a failure stops the experiment, and every later part refuses to run past it,
  `require_identity`): the ladder's identity point reproduces each wiring's saved fit with 0 count differences (16-odor
  input, APL on and off; factors 1) and, on the real wiring, E98's saved B1 run on the extended input (counts and APL's
  means; `identity_point`); pn_inhadd 100 to 3141 on the real wiring reproduce E77's screen (thresholds, scales, counts and
  scores; `e77_identity`, read from part B's runs once they exist); the extended-input runs at 1000 and 3141 on the real
  wiring reproduce E98's B4 (counts and APL window means on their common stimuli; `b4_identity`), and every run passes
  E98's amended blank identity (`a0_identity`).
- A (no simulations): A1 the chain of totals from the ORN input to APL per stimulus on the real wiring at C2 (E98's B1
  run), 1000 and 3141 (E98's B4 runs): the input's total, its broad component (36 x the median glomerulus) and the rest,
  the PN total and APL's response at the carrying level, as ratios for E98's, E100's and E101's stimuli, with the exponent
  a power law from the input's total to APL would need, in the flies and in the model, and the log split of the model's
  flattening into the antennal lobe's part and the mushroom body's (a1_chain.csv, a1_ratios.csv); A2 the KCs at low
  concentration under E98's two B4 settings: response rates with APL on and off at ma 10^-7 and 2h 10^-5 (exact input),
  2h 10^-6 and farn 10^-2 (standing for 10^-7 and 10^-3), against the flies' and C2's (a2_low_conc.csv); A3 E77's S1
  points under the calibrated readout and on raw counts with the amplitude and shape parts, from its saved pair distances,
  the per-odor part beside as a lead (a3_s1_parts.csv, a3_s1_per_odor.csv); A4 the prediction for every ladder step on S
  (pattern part, mean APL effect on raw counts, the paa ratio), written before part B (a4_prediction.csv,
  a4_prediction.json). a_summary.json.
- B (simulations; resumable; `--wirings`, `--steps`): each step on each wiring, rates re-held (`hold_point`), its 16-odor
  outputs and parameters saved, then the extended run (`run_extended`), under outputs/E102_antennal_lobe_normalization/B/
  <step>/ in E98's B1 layout plus <wiring>_family.csv (E101's family table) and <wiring>_point.json.
- score (no simulations; from part B's saved runs): every step's readings (b_readings.csv per unit and readout and
  target; b_points.csv per wiring; b_apl.csv E98's, E100's and E101's comparisons per step; b_low_conc.csv; b_pn.csv) and
  part C's levels (levels.json).
- D (the refits): the commands of `E10_joint_fit.py` with C2's options at each level (`refit_commands`; printed, not run),
  and their classes once fitted (`--score-refits`).

Outputs: results/E102_antennal_lobe_normalization/ (the tables above, identity.json, provenance_<part>.json); counts,
responses, APL and PN means of part B under outputs/E102_antennal_lobe_normalization/, archived to
data/model_outputs/E102_antennal_lobe_normalization/ with pn_inhadd, pn_inhsc, the thresholds and every scale as absolute
values. Every table holding Yang's recordings' values (the flies' intervals and rates) opens with the quality caveat
(YANG_QUALITY_CAVEAT; her answers of 2026-11-15 18:08) as a '#' line.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E102_antennal_lobe_normalization --part A
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E102_antennal_lobe_normalization --part 0
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E102_antennal_lobe_normalization --part B --wirings S
    ./venv/bin/python -u -m experiments.E102_antennal_lobe_normalization --part score
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E102_antennal_lobe_normalization'
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
ALL_WIRINGS = S_WIRINGS + T_WIRINGS
WIRING_SETS = {'S': S_WIRINGS, 'T': T_WIRINGS, 'all8': ALL_WIRINGS}

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any reading)

PN_INHADD_DEFAULT, PN_INHSC_DEFAULT = 31.4088, 368.6631     # olfsysm's (E3's PN_INHADD_DEFAULT, PN_INHSC_DEFAULT)
FIXED_GAIN_INHADD = 1e6                                     # the ladder's last step: the gain fixed
LADDER = (100.0, 314.0, 1000.0, 3141.0, 10000.0, FIXED_GAIN_INHADD)
STEP_IDENTITY = 'identity'
E77_S1 = (100.0, 314.0, 1000.0, 3141.0)                     # E77's screen's S1 points (part 0's second identity)
E98_B4 = (1000.0, 3141.0)                                   # E98's B4 points (part 0's third identity)

# the levels (part C) and the classes (part D)
LEVEL_MIN_WIRINGS = 6                                       # of the 8 (E98's VERDICT_MIN_WIRINGS)
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01   # E94's equal rule
EQUAL_TARGETS = ('z8', 'halfA', 'halfB')                    # an equal: the pattern part against all recordings, A and B
BETTER_MARGIN = 0.003
BETTER_HALVES = ('halfA', 'halfB', 'halfC', 'halfD')
T_BETTER_MIN_WIRINGS = 6                                    # of the 8 single wirings (E83's gate)
# the unchanged fits' pattern parts under the calibrated readout against all recordings (the entry's; part 0 recomputes)
BASE_PATTERN = {'S': 0.0500, 'T': 0.0486, 'all8': 0.0482}
BASE_PATTERN_TOL = 5e-5                                     # part 0: the recomputed values agree to the entry's 4 decimals

# the KCs at low concentration (outcome (e))
LOW_CONC_FACTOR = 1.25
LOW_STIMULI = (                                             # (B1 label of the model's input, the flies' stimulus, role)
    ('ma [kiwi-control/diag] @ -7', ('ma', -7.0), 'exact'),
    ('2h [kiwi-control/control] @ -5', ('2h', -5.0), 'exact'),
    ('2h [kiwi-control/diag] @ -6', ('2h', -7.0), 'stands for 10^-7'),
    ('farn [kiwi-control/diag] @ -2', ('farn', -3.0), 'stands for 10^-3'))
LOW_CARRYING = 'ma [kiwi-control/diag] @ -7'                # outcome (e)'s stimulus (E49's source)
LOW_BESIDE = ('ma [megamat/diag] @ -7', 'ma [validation2/diag] @ -7')
FLIES_OFF_ON_MA7 = 2.29                                     # the flies' APL-off over APL-on response rate at ma 10^-7

# refits (part D): E10_joint_fit.py with C2's options
REFIT_ARGS = ('-k', '2', '--aplpn-max', '3', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both')
REFIT_MAXFEV = 200

YANG_QUALITY_CAVEAT = ("Yang has not checked these recordings' quality (her answers of 2026-11-15 18:08: all data can be "
                       "used, she is unsure of quality questions); the flies' values are read from E98's, E100's, E101's "
                       "and E72's saved tables")


# ---------------------------------------------------------------------------
# folders, records

def out_dir() -> Path:
    """results/E102_antennal_lobe_normalization/ in the shared directory."""
    raise NotImplementedError


def local_dir(sub: str = '') -> Path:
    """outputs/E102_antennal_lobe_normalization/<sub>/ in this repository."""
    raise NotImplementedError


def step_tag(step) -> str:
    """A step's folder and tag name: 'identity', 'inhadd100' ... 'inhadd10000', 'gain_fixed' (ASCII only)."""
    raise NotImplementedError


def step_kws(step) -> Dict[str, float]:
    """A step's PN kwargs: pn_inhadd and pn_inhsc (`E77_input_screen.inhsc_for`; the defaults at STEP_IDENTITY)."""
    raise NotImplementedError


def write_table(df: pd.DataFrame, path: Path, caveat: bool = False) -> None:
    """A results table; with `caveat` its first line is '# ' + YANG_QUALITY_CAVEAT (report_values skips it)."""
    raise NotImplementedError


def read_table(path: Path) -> pd.DataFrame:
    """A results table written by `write_table` (leading '#' lines skipped)."""
    raise NotImplementedError


def require_identity(which: str = 'all') -> None:
    """Refuses to go on past a failed (or missing) identity of part 0: raises unless identity.json says it passed."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B: one point (a step on a wiring)

def point_breadth(wiring: str, kws: Mapping[str, float]) -> pd.Series:
    """The breadth behind the APL->KC weights at a step on one wiring (`E77_input_screen.handle_breadth` with the
    wiring's claw x bouton matrix)."""
    raise NotImplementedError


def hold_point(wiring: str, step) -> dict:
    """Both response rates re-held at a step on one wiring (`E77_input_screen.hold_rates` on E83's simulator with the
    step's PN kwargs and breadth; the wiring's fit's targets and parameters). Returns hold_rates' dict with 'kws',
    'step', 'wiring'."""
    raise NotImplementedError


def run_extended(wiring: str, step, held: Mapping) -> Dict[str, dict]:
    """The held point on E98's extended input (`E98.b1_input`), APL on (held params) and acutely off (params_off): per
    condition E98's outputs (responses, spike counts, APL's means with the settled level, PN window and pre-odor means,
    APL's trace) and E90's claw stage (`read_stages`)."""
    raise NotImplementedError


def run_point(wiring: str, step, root: Optional[Path] = None, force: bool = False, extended: bool = True) -> dict:
    """One point of part B (resumable): `hold_point`, the 16-odor outputs and the parameters saved, then `run_extended`
    saved in E98's B1 layout with E101's family table; checks that the extended run's 16 panel columns equal the held
    run's counts, and E98's amended blank identity. Returns the point's record (<wiring>_point.json)."""
    raise NotImplementedError


def part_b(wirings: Sequence[str], steps: Sequence, root: Optional[Path] = None) -> None:
    """Part B's points on `wirings` x `steps` (resumable; every point's identity checked)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part 0

def identity_point(wiring: str, extended: bool = False) -> dict:
    """The ladder's identity point on one wiring: factors 1 and 0 count differences against the saved fit (16 odors,
    APL on and off), the breadth equal to the wiring's; with `extended`, the extended run against E98's saved B1 run
    (real wiring) or E101's A3 rows (counts of the 16 panel odors, APL's window and pre-odor means)."""
    raise NotImplementedError


def e77_identity(root: Optional[Path] = None) -> dict:
    """pn_inhadd 100 to 3141 on the real wiring (part B's runs) against E77's screen: thresholds, scales, response
    counts, spike counts and scores (screen.csv, screen_params.json, the archived screen outputs)."""
    raise NotImplementedError


def b4_identity(root: Optional[Path] = None) -> dict:
    """The extended runs at 1000 and 3141 on the real wiring against E98's B4 (counts on the common stimuli, APL window
    means) and every run's amended blank identity."""
    raise NotImplementedError


def part_0(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """identity.json: `identity_point` on `wirings` (the extended check on the real wiring), the unchanged fits'
    populations against BASE_PATTERN, and `e77_identity` / `b4_identity` where part B's runs exist."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A (no simulations)

def chain_table(runs: Mapping[str, Mapping[str, dict]], stimuli: pd.DataFrame) -> pd.DataFrame:
    """A1's chain of totals per setting and stimulus (APL on): the ORN input's total, positive total, median glomerulus,
    broad component and the rest; the PN total (all PN rows and the input's glomeruli) at the blank's level and E90's;
    APL's window mean at the blank's level, the settled level and E90's."""
    raise NotImplementedError


def chain_ratios(chain: pd.DataFrame) -> pd.DataFrame:
    """A1's ratios: per setting and fixed ratio (E98's V1, V2, V3 to V5; E100's W1 six; E101's pa over ms), the input's,
    the PN totals' and APL's ratios, the flies' APL ratio with its interval, the power-law exponents and the log split of
    the model's flattening (antennal lobe, mushroom body)."""
    raise NotImplementedError


def low_conc_table(runs: Mapping[str, Mapping[str, dict]]) -> pd.DataFrame:
    """A2 (and part B's): per setting, LOW_STIMULI's response rates with APL on and off, the off / on ratio, against the
    flies' (`E100.kc_data_rates`)."""
    raise NotImplementedError


def s1_parts() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: E77's S1 points (screen_pairs.csv) under the calibrated readout and on raw counts against all recordings: the
    loss, the mean APL effect, the pattern part with its amplitude and shape parts, the spread ratio, ρ and Pearson r on
    the per-pair APL effect; and the per-odor part of the pattern part's error (a lead)."""
    raise NotImplementedError


def prediction(a3: pd.DataFrame, ratios: pd.DataFrame) -> pd.DataFrame:
    """A4: the prediction for every ladder step on S, written before part B."""
    raise NotImplementedError


def part_a() -> dict:
    """Part A (no simulations): A1 to A4 (module docstring)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# score, part C's levels, part D's classes

def step_readings(step, root: Optional[Path] = None) -> Dict[str, pd.DataFrame]:
    """One step's readings from its saved runs (module docstring, score)."""
    raise NotImplementedError


def equal_rule(d_pattern: Mapping[str, float], d_raw_loss: float, d_mean: float, n_flagged: int, n_at_bound: int) -> dict:
    """E94's equal rule on S: 'equal' if the pattern part is within EQUAL_PATTERN against every EQUAL_TARGETS, no run is
    flagged or at a bracket's bound, the raw-count loss within EQUAL_RAW and the mean APL effect within EQUAL_MEAN."""
    raise NotImplementedError


def levels(steps: pd.DataFrame) -> dict:
    """Part C's levels from one row per step (value pn_inhadd; v2 E98's held V2 class; v2_linear the carrying reading's
    linear category; equal_S the step's equal on S; flagged): L_apl, L_lin, L_free (None where none exists)."""
    raise NotImplementedError


def refit_levels(lv: Mapping, steps: pd.DataFrame) -> List[float]:
    """The levels the refits run at: L_apl (or the last unflagged step when it does not exist), and L_lin if it exists
    and differs."""
    raise NotImplementedError


def refit_class(d_pattern: Mapping[str, float], d_raw_loss: float, d_mean: float, n_flagged: int) -> str:
    """'equal', 'better' or 'price' for S's refitted population against the unchanged fits."""
    raise NotImplementedError


def t_better_holds(d_pattern_t: Mapping[str, float], d_single: Mapping[str, float], d_raw_loss_t: float) -> bool:
    """Whether 'better' stands on T (E83's gate thresholds)."""
    raise NotImplementedError


def low_conc_class(model: float, flies: float, model_c2: float, factor: float = LOW_CONC_FACTOR) -> str:
    """'further', 'nearer' or 'about the same' (outcome (e))."""
    raise NotImplementedError


def refit_commands(level: float, wirings: Sequence[str] = S_WIRINGS) -> List[str]:
    """Part D's commands (E10_joint_fit.py with C2's options, `--pn-inhadd`, from each wiring's point of part B)."""
    raise NotImplementedError


def part_score(root: Optional[Path] = None) -> dict:
    """The score (module docstring) and part C's levels."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the CLI

def _wirings(arg: str) -> Tuple[str, ...]:
    if arg in WIRING_SETS:
        return WIRING_SETS[arg]
    w = tuple(x.strip() for x in arg.split(',') if x.strip())
    bad = [x for x in w if x not in ALL_WIRINGS]
    if bad:
        raise SystemExit(f'unknown wirings {bad}')
    return w


def _steps(arg: str) -> Tuple:
    if arg == 'ladder':
        return LADDER
    out = []
    for x in arg.split(','):
        x = x.strip()
        out.append(STEP_IDENTITY if x == STEP_IDENTITY else float(x))
    return tuple(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'score', 'D'])
    ap.add_argument('--wirings', default='all8', help="'S', 'T', 'all8' or a comma-separated list")
    ap.add_argument('--steps', default='ladder', help="'ladder' or comma-separated pn_inhadd values ('identity' allowed)")
    ap.add_argument('--root', default=None, help="part B: put the runs under another folder (default the local outputs)")
    ap.add_argument('--force', action='store_true', help='part B: run points that are saved again')
    ap.add_argument('--level', type=float, default=None, help='part D: the level (pn_inhadd) to print the refits for')
    a = ap.parse_args()
    root = Path(a.root) if a.root else None
    if a.part == '0':
        print(json.dumps(part_0(_wirings(a.wirings)), indent=1, default=str))
    elif a.part == 'A':
        print(json.dumps(part_a(), indent=1, default=str))
    elif a.part == 'B':
        part_b(_wirings(a.wirings), _steps(a.steps), root)
    elif a.part == 'score':
        print(json.dumps(part_score(root), indent=1, default=str))
    elif a.part == 'D':
        if a.level is None:
            raise SystemExit('--part D needs --level')
        print('\n'.join(refit_commands(a.level, _wirings(a.wirings))))


if __name__ == '__main__':
    main()
