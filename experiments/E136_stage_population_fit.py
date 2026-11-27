#!/usr/bin/env python3
"""E136, the best family fitted to all megamat stages at once: one loss over the PN, bouton, claw and soma populations,
five structural arms, a control fit on odor-permuted flies, and reads on held-out wirings, held-out flies and three
panels the fit never sees.

The entry is item E136 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04; its code claimed by disco 2026-11-24 13:29, after
atlas handed it on; its nine open definitions ruled by a Fable subagent, claude-fable-5-2, spawned from disco
2026-11-24 13:36, and written into the entry at 14:00). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part
B's eight wirings, each at its own fit (E83's `members`, `saved_fit`, `build_simulator`); S is real, P1, P2, M1 and T is
P3, P4, M2, M3. The measures, the readings and the stage table are those of the models as populations at every recorded
stage (E134): the model read as each stage is recorded (atlas's kcq/stage_readings.py), the views, one-fly units and
population errors (hal's kcq/population_fit.py), the flies' populations and the lock (viki's
data_prep/stage_population_data.py, experiments/E134_lock.json), the family built and run with every stage saved
(atlas's experiments/E134_population_baseline.py: `build_family`, `run_stages` with its `claw_reader`, `held_targets`).
This module calls them and copies none of them; a change those files need is asked of their owners.

The guard: nothing is fitted before E134's measure gate reads "usable" (results/E134_population_baseline/
measure_gate.json, written by viki from the gate planner's reading) and its stage table is in (score_classes.json,
score_stage_table.csv); `require_gate` refuses otherwise. A1 also waits for E135's class (`require_e135`); the other
arms do not.

The loss (the entry's, with the rulings): the megamat stage error over the stages that E134 found readable
(`readable_stages`), on the S wirings' mean views (hal's `wirings_as_flies`, row 'common'), against the fit half of
each stage's flies (`split_names`: hal's `split_halves` at seed 0, the first part; 3 of the 7 PN-dendrite flies, 2 of
the 5 bouton flies, 2 of the 4 claw flies, 2 of the 4 soma flies), each stage in its default reading with its curve set
anew at every evaluation; the half gives the targets (its mean views), the curve and the noise, and the one-fly units
are those of all the stage's flies (`P.subset_reference(..., keep_units=True)`, `Flies`). The PN dendrites on four
views, the other stages on three (E134's `megamat_stage_error`).

The operating point: both response rates are held on every wiring at every evaluation, by two spikes, to within
HOLD_TOL of E134's held megamat targets for that wiring (`B.held_targets`: APL on at the lock's 0.0915, APL off at 0.0915
times the wiring's own APL-off / APL-on ratio on its APL16 fit): the APL-off count sets the somatic threshold, then the
APL-on count the APL->KC scale with its two tied scales, each by a secant in log factor warm-started from the nearest
evaluated point (`secant_hold`, at most HOLD_MAX_SIMS simulations per rate, within HOLD_RANGE of the start). A hold
that misses adds one one-fly unit per HOLD_MISS_PER of relative miss. The reads are held exactly (E77's `hold_rates`).

The arms (`ARMS`): each searches shared log factors on each wiring's own parameters, the threshold and the APL->KC
scale being set by the operating point (so A0, A1, A2 have four parameters, A3 five, A4 three):
- A0, the family as it is: the APL->bouton and KC-KC factors;
- A1, with the antennal-lobe stage that E135 hands on (its table of PN rates by E135's route, `pn_rates`, fixed; if
  it hands on none, the normalization's strength and a gain on pn_inhsc searched too): waits on E135;
- A2, APL's inhibition of the claws as a gain (E92's form, `apl_inh_divisive`), its claws read as olfsysm divides
  them (`claw_inputs_divisive`, identical to olfsysm's claw_sims, `claw_identity`);
- A3, the per-claw threshold searched (E16's form, `claw_thr` in Hz, coordinate claw_thr / CLAW_THR_UNIT), set in place
  on the built simulator where that equals a fresh build exactly (`claw_thr_check`), else rebuilt at every evaluation;
  its claws read before the threshold, with the reading after it beside at the reads;
- A4, without the inhibition among KCs (the KC-KC scale at 0): the APL->bouton factor alone;
- the control (CONTROL): A0 against flies whose odor labels are permuted, one permutation for all flies and stages,
  its loss the mean over N_PERMS seeded permutations that leave no odor in place (`permutations`).
The search (`fit_arm`): Nelder-Mead from the family's held megamat points of E134 (every factor 1, no per-claw
threshold), first steps FIRST_STEP, XATOL / FATOL, MAXFEV by arm, one restart from the best point with half the steps
if it stops with RESTART_MIN_LEFT or more evaluations left (that point kept in fit_restart_<arm>.json, so that a resumed
fit restarts from it again; `restart_point`); a point outside BOUNDS is evaluated at the clipped point plus
one one-fly unit per BOUND_PENALTY_PER of overshoot. In A0, A2, A3 and A4 nothing searched acts before the boutons, so
the PN-dendrite error must be the same at every evaluation (checked).

Parts (CLI `--part`):
- identity: A0 at z = 0 against all flies: the hold changes nothing, every run reproduces E134's saved megamat_held
  counts and readings exactly, and each wiring's population errors reproduce E134's score_scores.csv where it exists.
- claw_identity: A2's divided claw drive against olfsysm's claw_sims; claw_thr_check: A3's in-place threshold against a
  fresh build at CLAW_THR_CHECK_HZ.
- fit: one arm's fit (`--arm`, CONTROL for the control), resumable from its history (on any host: `archive` copies
  the evaluations the earlier host saved, `--evals`, by default its best, so that the fit's best can be archived).
- read: the reads at the fitted parameters, held exactly at each panel's operating point (`planned_reads`, lanes and
  claims as E134's step 3): megamat on T against the held-out half; validation2 (the lock's 0.0504) and the mixtures
  (E134's held mixture rate) with the wiring's own ratio for APL off; the APL16 input at the wiring's own fit's two
  counts. The unrefitted family's reads are E134's megamat_held and mixtures_held runs and E90's part B (the APL16
  fit) where they exist; only its validation2 reads are run here. `read_at` reads an arm's model at any coordinates,
  with a tag, claims root and run folders of the caller's (E139's points; `read_arm` is it at an arm's fitted point).
- score: the reads scored and classed (`score_megamat`, `score_panels`, `apl16_rows`, `classify`, `part_score`); a
  variant whose runs live elsewhere is scored by the same functions once registered (`register_variant`, its loader
  consulted by `load_read`; read on every wiring, where a fitted arm's S wirings are its fit's best evaluation).
- variants: variants.csv, one row per run in the model-output archive with its parameters, for the archive's
  params.csv (`write_variants`; archive_tools/archive_model_outputs.py reads it).

The classes (the entry's, read as ruled): each wiring is one fly; a value is the mean over the four T wirings of each
wiring's error (the error of the T wirings' mean views beside); a change is the mean over T of the arm minus the
unrefitted family on the same wiring; a spread is the SD over the family's eight wirings of the unrefitted family's
value; the standard error over wirings is the larger of half the megamat stage error's spread and the SD of the
per-wiring changes over T divided by 2. "improves": the change is below -IMPROVE_SE standard errors and no stage's
change is above that stage's spread; "holds": it improves and no held-out panel (validation2, cmix, APL16 with APL
on, APL16 with APL off) is worse by more than its spread (kmix is scored beside, with E137's two readings of its
components as bounds: Tom, 2026-11-25 10:42); "megamat only"; "one stage at another's cost"; "adds" on
the per-wiring changes from A0; no class outside a band, at a bound, or no better than the control.

Realism (reported, not optimized): the APL16 reads' bands (kcq.apl_effect's, the KC data in baseline-noise units), the
megamat rates beside 0.0915, the largest counts, every fitted scale beside the family's; named in advance: the
APL->bouton factor against its bound, the bouton share of APL's inhibition, a per-claw threshold above 70 Hz, A1's
glomerulus terms, A2's divisor at the responders' claws and the KC-KC scale (a fitted term with no named carrier in the
calyx, Tom 2026-11-24 via hal). A fit outside a band gets a "Realism flag:" post and no class.

Concentrations: the fit sees megamat only, at 10^-3 on both sides. The validation2 read is one panel at matched
concentrations, and the mixtures are used as measured. The APL16 read involves paa, euc and -aPine, whose input is the
ORN data within 1 log10 unit of the flies' concentrations.

Holds (CLAUDE.md): only Pratyush's hemibrain (the family's wirings); no olfsysm change (A2's and A3's forms are options
at olfsysm b1e89d7; A1's stage enters by `pn_rates`); no input search; only Yang's processed products; Remy's processed
data without the bouton fly 2022-10-21 fly02 (viki's lock excludes it). Nothing is published.

Outputs:
- outputs/E136_stage_population_fit/ (local, large): fit_<arm>/ (every evaluation's stages, E134's step-3 layout),
  read/ (the reads' runs, also in the archive data/model_outputs/E136_stage_population_fit/read/), logs/.
- <shared>/results/E136_stage_population_fit/ (compact): fit_history_<arm>.csv, fit_<arm>.json, identity_A0.json,
  claw_identity_A2.json, claw_thr_check.json, reads_<host>_<lane>.csv, provenance_read_<host>_<lane>.json, claims/,
  score_*.csv, classes.json.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \
        experiments/E136_stage_population_fit.py --part fit --arm A0
"""
from __future__ import annotations

import argparse
import contextlib
import gc
import json
import os
import socket
import time
from collections import namedtuple
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import stage_readings as R

EXPERIMENT = 'E136_stage_population_fit'
SCRIPT = 'experiments/E136_stage_population_fit.py'
E134 = 'E134_population_baseline'
GATE_FILE = 'measure_gate.json'
# the fitted panel, and the E134 group the fits start from (its held megamat points)
PANEL = 'megamat'
START_GROUP = 'megamat_held'

# the family's wirings (E83's, as E134 takes them; `check_wirings` checks they are E134's)
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS

# the stages a loss can read (E134's MEGAMAT_ERROR_STAGES); which enter is E134's verdict (`readable_stages`)
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
HALVES: Tuple[str, ...] = ('fit', 'held')
HALF_SEED = 0
# the APL-off / APL-on ratio the entry named (the wirings' own ratios, 2.09 to 2.46, average it; reported only)
OFF_ON_RATIO = 2.2
MAX_FREE = 8
READ_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')

# the log factors: name -> the absolute parameters it multiplies (E72's TIED: the KC->APL and bouton->APL scales follow
# the APL->KC scale at each wiring's own ratios); 'thr' and 'apl' are set by the operating point's hold (HELD), the others
# searched; the scales APL off sets to 0
FACTORS: Dict[str, Tuple[str, ...]] = {
    'thr': ('fixed_thr',),
    'apl': ('wAPLKC_scale', 'wKCAPL_scale', 'wPNAPL_scale'),
    'aplpn': ('wAPLPN_scale',),
    'kckc': ('wKCKC_scale',),
}
HELD: Tuple[str, ...] = ('thr', 'apl')
OFF_ZERO: Tuple[str, ...] = ('wAPLKC_scale', 'wAPLPN_scale')

# the arms: name -> (description, the searched coordinates, the fixed model kwargs beside the family's)
Arm = namedtuple('Arm', 'name description searched kws')
ARMS: Dict[str, Arm] = {
    'A0': Arm('A0', 'the family as it is, refitted under the new loss', ('aplpn', 'kckc'), {}),
    'A1': Arm('A1', "with the antennal-lobe stage that E135 hands on, fixed", ('aplpn', 'kckc'), {'pn_rates': 'E135'}),
    'A2': Arm('A2', "APL's inhibition of the claws as a gain (E92's form)", ('aplpn', 'kckc'),
              {'apl_inh_divisive': True}),
    'A3': Arm('A3', "the per-claw threshold searched (E16's form)", ('aplpn', 'kckc', 'claw_thr'), {}),
    'A4': Arm('A4', "without the inhibition among KCs (E110's form)", ('aplpn',), {'kckc_fixed': 0.0}),
}
CONTROL = 'A0_perm'
# A1's form fitted to the permuted flies, built as A1 (E145's four-factor control; `fit_arm(..., control=True)` of A1)
CONTROL_A1 = 'A1_perm'
CONTROLS: Dict[str, str] = {CONTROL: 'A0', CONTROL_A1: 'A1'}
FAMILY = 'family'   # the unrefitted family in the reads: z = 0, each wiring at its held megamat point

# the search (the rulings): the coordinates (ln factors; the per-claw threshold in units of CLAW_THR_UNIT Hz), their
# bounds, the first steps, the tolerances, the budgets and the restart
CLAW_THR_UNIT = 50.0
BOUNDS: Dict[str, Tuple[float, float]] = {
    'aplpn': (float(np.log(1.0 / 30.0)), float(np.log(1.07))),
    'kckc': (float(np.log(0.1)), float(np.log(7.3))),
    'claw_thr': (0.0, 100.0 / CLAW_THR_UNIT),
    # A1's fallback (PN_COORDS): pn_inhadd from a tenth of olfsysm's default to about E102's top rung (10^4), and the
    # gain on pn_inhsc within a factor of three
    'inhadd': (float(np.log(0.1)), float(np.log(320.0))),
    'inhsc': (float(np.log(1.0 / 3.0)), float(np.log(3.0))),
}
# inhadd's first step is E102's first rung (pn_inhadd 100 over olfsysm's default 31.4088, E3's PN_INHADD_DEFAULT)
FIRST_STEP: Dict[str, float] = {'aplpn': -0.7, 'kckc': -0.7, 'claw_thr': 0.6, 'inhadd': float(np.log(100.0 / 31.4088)),
                                'inhsc': 0.4}
BOUND_PENALTY_PER = 0.1
XATOL, FATOL = 0.02, 0.005
MAXFEV: Dict[str, int] = {'A0': 80, 'A1': 80, 'A2': 80, 'A3': 110, 'A4': 40, CONTROL: 80}
# A1 with the fallback's two more coordinates (four searched, as A3's three take 110)
MAXFEV_A1_SEARCH = 140
RESTART_MIN_LEFT = 15
# "at a bound": within 2 % of a factor's bound, or 2 Hz of the threshold's
AT_BOUND_REL, AT_BOUND_HZ = 0.02, 2.0
# the hold (the rulings): tolerance, simulations per rate, ranges of the factors, the miss's penalty, the secant's
# first slope (d ln responses / d ln factor) where nothing is known yet
HOLD_TOL = 0.005
HOLD_MAX_SIMS = 8
# a wiring's first hold has no evaluated point to start from (an arm's change at the family's start, e.g. A4's KC-KC at
# 0, can move the threshold several-fold: about 4.5-fold in the smoke test of 2026-11-24), so it may take as many
# simulations as E69's bisect_factor
HOLD_FIRST_MAX_SIMS = 20
HOLD_RANGE: Tuple[float, float] = (0.1, 10.0)
HOLD_RANGE_DIVISIVE_APL: Tuple[float, float] = (0.01, 10.0)
HOLD_MISS_PER = 0.1
HOLD_SLOPE0: Dict[str, float] = {'thr': -3.0, 'apl': -1.0}
# the control's permutations
PERM_SEED, N_PERMS = 136, 20
# A3's in-place check
CLAW_THR_CHECK_HZ = 40.0
# the classes' margin, in standard errors over the wirings
IMPROVE_SE = 2.0
# the identity's tolerance on the population errors against E134's score ("to rounding")
IDENTITY_TOL = 1e-9
# the start check (the re-read on the adopted measure): A0's first evaluation (z = 0, E134's held points) on strength
# against the fit half, as the planner read it from E134's saved runs (4 decimals): the loss, then the four stages
START_CHECK: Dict[str, Dict[str, float]] = {
    'strength': dict(stage_error=0.9747, pn_dendrite_error=1.1284, pn_bouton_error=1.0516, kc_claw_error=0.7938,
                     kc_soma_error=0.8892),
}
START_CHECK_TOL = 6e-5
# a factor whose return to the family's value raises the loss by less than this is "not determined" by the measure
# (the re-read on the adopted measure: half the megamat stage error's spread on strength, 0.020)
NOT_DETERMINED = 0.010
# A1's PN stage (`require_e135`): E135's outcome (a) hands on R1's table, (b) the kept rung's (its step 1's
# 'step2_rung'), fixed and entering by E135's route 1 (olfsysm's `pn_rates`); (c) to (f) hand on none, and A1 searches
# the normalization instead (PN_COORDS: ln of pn_inhadd over olfsysm's default, pn_inhsc following as in E102's ladder,
# `E77_input_screen.inhsc_for`, and ln of a gain on that pn_inhsc), the simulator rebuilt at every evaluation
E135_TABLE_OUTCOMES: Dict[str, str] = {'a': 'R1', 'b': 'rung'}
E135_NONE_OUTCOMES: Tuple[str, ...] = ('c', 'd', 'e', 'f')
PN_COORDS: Tuple[str, ...] = ('inhadd', 'inhsc')

ReadPoint = namedtuple('ReadPoint', 'arm panel wiring half')

# the measure the loss and the classes read (E134's gate of 2026-11-24 15:19 read the three-view population error "not
# usable" and left the measure to Tom and Yang; the planner's cut is strength by odor alone, with the two-view error
# beside): 'all' every view (the PN dendrites on four, the matrix view included), or the views named. A measure other than
# 'all' writes every output in a folder of its own (`out_dir`, `local_dir`: measure_<name>/). Set by --measure
# (`set_measure`).
MEASURES: Dict[str, Optional[Tuple[str, ...]]] = {
    'all': None,
    'strength': ('strength',),
    'strength_sharing': ('strength', 'sharing'),
}
MEASURE = 'all'


def set_measure(name: str) -> None:
    """Sets the measure in force (MEASURES)."""
    global MEASURE
    if name not in MEASURES:
        raise ValueError(f'measure {name!r} not in {sorted(MEASURES)}')
    MEASURE = name


def measure_parts(views: Optional[Sequence[str]]) -> Optional[List[str]]:
    """hal's parts of the views named (None for every part)."""
    from kcq import population_fit as P
    return None if views is None else [p for p in P.ALL_PARTS if p.split('/')[0] in views]


def error_column(stage: str, views: Optional[Sequence[str]]) -> str:
    """The population error's column of a stage under a measure: the PN dendrites on four views where the matrix view
    is read, otherwise the error over the views read."""
    return 'error_with_matrix' if stage == 'pn_dendrite' and (views is None or 'matrix' in views) else 'error'


# ---------------------------------------------------------------------------
# paths and guards

def _measure_sub() -> str:
    return '' if MEASURE == 'all' else f'measure_{MEASURE}'


# folders another experiment's parts have E136's parts write under, in place of E136's own (`paths`): 'out' (the
# results: histories, fit records, restarts), 'local' (the outputs: each evaluation's stages) and 'archive' (the
# model-output archive: a fit's best evaluation); empty outside that context
PATHS: Dict[str, Path] = {}


@contextlib.contextmanager
def paths(out: Optional[Path] = None, local: Optional[Path] = None, archive: Optional[Path] = None):
    """For the block, `out_dir`, `local_dir` and `archive_dir` return the folders given (as they are: no measure
    subfolder; each one not given keeps E136's), so that E136's fit, its single-factor returns and its rescore
    (`fit_arm`, `determined`, `rescore`) write another experiment's records (E139's step 3). The reads' folders
    (`read_dirs`) and E134's are not changed."""
    old = dict(PATHS)
    PATHS.update({k: Path(v) for k, v in (('out', out), ('local', local), ('archive', archive)) if v is not None})
    try:
        yield
    finally:
        PATHS.clear()
        PATHS.update(old)


@contextlib.contextmanager
def extra_arms(*arms: Arm):
    """For the block, ARMS also holds `arms` (each an `Arm` whose name is none of E136's: an arm, CONTROL, FAMILY or a
    registered variant), so that E136's fit and reads build them as their own (`ArmModel`, `fit_arm`, `hold_at`,
    `read_at`; e.g. E139's A0 with the normalization fixed); they are removed after. Their budgets are passed to
    `fit_arm` (MAXFEV holds E136's arms only)."""
    added = []
    try:
        for a in arms:
            if a.name in ARMS or a.name in (*CONTROLS, FAMILY) or a.name in VARIANT_LOADERS:
                raise ValueError(f'{a.name!r} is already a name of E136 or a registered variant')
            ARMS[a.name] = a
            added.append(a.name)
        yield
    finally:
        for n in added:
            ARMS.pop(n, None)


# which of `split_names`' two parts is the fit half: 0 E136's (its first part), 1 the other (inside `cell`)
FIT_PART = 0


@contextlib.contextmanager
def cell(searched: Sequence[str], fit_part: int, start_check: Optional[Mapping[str, float]] = None):
    """For the block, E136's fit, single-factor returns, rescore, reads and score functions work in another cell of
    wirings by halves of the flies (E145, the fit of all megamat stages made again): S_WIRINGS (the searched wirings) is
    `searched`, four of WIRINGS, and T_WIRINGS (the read wirings) the other four in WIRINGS' order; the fit half 'fit' is
    part `fit_part` of `split_names`' split (0 E136's fit half, 1 its held-out half) and 'held' the other part, so that
    'fit' and 'held' name the cell's own halves. `start_check` (the measure in force's start values, START_CHECK's keys:
    the family at z = 0 on the searched wirings against the cell's fit half, read from saved runs) replaces START_CHECK
    for A0's and A1's first evaluation; without it, in a cell other than E136's, the first evaluation is not checked.
    WIRINGS, and so every wiring's seed index (`seed_index`), are unchanged. Restored after, also on an error.
    ValueError for wirings that are not four distinct names of WIRINGS, or a part other than 0 and 1."""
    global S_WIRINGS, T_WIRINGS, FIT_PART, START_CHECK
    searched = tuple(str(w) for w in searched)
    if len(searched) != len(WIRINGS) // 2 or len(set(searched)) != len(searched) or not set(searched) <= set(WIRINGS):
        raise ValueError(f'a cell searches four distinct wirings of {WIRINGS}: {searched}')
    if fit_part not in (0, 1):
        raise ValueError(f'{fit_part=}: 0 or 1')
    old = S_WIRINGS, T_WIRINGS, FIT_PART, START_CHECK
    own = (tuple(searched), int(fit_part)) == (tuple(old[0]), int(old[2]))
    S_WIRINGS = searched
    T_WIRINGS = tuple(w for w in WIRINGS if w not in searched)
    FIT_PART = int(fit_part)
    if start_check is not None:
        START_CHECK = {MEASURE: {k: float(v) for k, v in start_check.items()}}
    elif not own:
        START_CHECK = {}
    try:
        yield
    finally:
        S_WIRINGS, T_WIRINGS, FIT_PART, START_CHECK = old


def out_dir() -> Path:
    """<shared>/results/E136_stage_population_fit/[measure_<name>/] (created); PATHS['out'] inside `paths`."""
    if 'out' in PATHS:
        PATHS['out'].mkdir(parents=True, exist_ok=True)
        return PATHS['out']
    d = H.results_dir(EXPERIMENT)
    if _measure_sub():
        d = d / _measure_sub()
        d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    """outputs/E136_stage_population_fit/[measure_<name>/][<part>/] (created); PATHS['local']/[<part>/] inside
    `paths`."""
    if 'local' in PATHS:
        d = PATHS['local'] / part if part else PATHS['local']
        d.mkdir(parents=True, exist_ok=True)
        return d
    d = H.outputs_dir(EXPERIMENT)
    for x in (_measure_sub(), part):
        if x:
            d = d / x
    d.mkdir(parents=True, exist_ok=True)
    return d


def e134_dir() -> Path:
    """<shared>/results/E134_population_baseline/."""
    return H.results_dir(E134)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _commits() -> dict:
    return {k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}


def check_wirings() -> None:
    """S_WIRINGS and T_WIRINGS are E134's (and so E83's)."""
    from experiments import E134_population_baseline as B
    if (S_WIRINGS, T_WIRINGS) != (tuple(B.S_WIRINGS), tuple(B.T_WIRINGS)):
        raise RuntimeError(f"the wirings differ from E134's: {S_WIRINGS} {T_WIRINGS}")


def seed_index(wiring: str) -> int:
    """A wiring's seed index in atlas's readings (its position in WIRINGS, as E134's score reads the family), so that
    a wiring is pooled and noised the same whichever set of wirings is read with it."""
    return WIRINGS.index(wiring)


def require_gate() -> dict:
    """E134's measure gate (results/E134_population_baseline/measure_gate.json) and its classes (score_classes.json);
    RuntimeError unless the gate's verdict is 'usable', the stage table is in, and, where the gate names an adopted
    measure (its 'measure' block: 'views', 'adopted'; strength by odor alone since 2026-11-24 15:28), the measure in
    force reads exactly those views (the re-read on the adopted measure, planned by a Fable subagent, claude-fable-5-2,
    from disco 15:26). Returns {'gate', 'classes'}."""
    g = _read_json(e134_dir() / GATE_FILE)
    if not g:
        raise RuntimeError(f'E134 has no {GATE_FILE} yet: no fit before its measure gate reads "usable"')
    if str(g.get('verdict')) != 'usable':
        raise RuntimeError(f"E134's measure gate reads {g.get('verdict')!r}, not 'usable'")
    meas = g.get('measure') or {}
    if meas.get('views') is not None:
        want = MEASURES[MEASURE]
        if not meas.get('adopted'):
            raise RuntimeError(f"the gate's measure {meas['views']} is not adopted")
        if want is None or set(want) != set(meas['views']):
            names = [k for k, v in MEASURES.items() if v is not None and set(v) == set(meas['views'])]
            raise RuntimeError(f"the adopted measure is {meas['views']}, not {MEASURE!r}: run with --measure "
                               f"{names[0] if names else '?'}")
    c = _read_json(e134_dir() / 'score_classes.json')
    if not c or not (e134_dir() / 'score_stage_table.csv').exists():
        raise RuntimeError("E134's stage table is not in (score_classes.json, score_stage_table.csv)")
    return dict(gate=g, classes=c)


def e135_verdict_path() -> Path:
    """E135's score verdict (results/E135_antennal_lobe_dendrites/score_verdict.json, E135's `outcome` under 'outcome')."""
    from experiments import E135_antennal_lobe_dendrites as AL
    return AL.out_dir() / AL.SCORE_VERDICT_FILE


def require_e135() -> dict:
    """What E135 hands on to A1, from its score's verdict (`e135_verdict_path`): outcome (a) R1's table, (b) the kept
    rung's (step1_fits.json's 'step2_rung'), mode 'table' (A1's PN stage fixed as that table, `a1_stage`); (c) to (f)
    none, mode 'search' (the normalization searched, PN_COORDS). An outcome that hands on a stage is taken whatever its
    guard says; a raised guard ('improvement' False, 'at_a_cost_beside') is carried into the record. RuntimeError while
    the verdict is not written, is incomplete, or waits (the dF/F repeat; a unit- or reading-dependent outcome, which
    goes to Tom and Yang).

    Returns dict(outcome, mode, rung (None in 'search'), improvement, at_a_cost_beside, reading_dependent, cls,
    verdict_sha256, verdict_path)."""
    import hashlib
    from experiments import E135_antennal_lobe_dendrites as AL
    path = e135_verdict_path()
    v = _read_json(path)
    if not v:
        raise RuntimeError(f'A1 waits on E135: its score verdict is not written ({path})')
    oc = v.get('outcome') or {}
    o = oc.get('outcome')
    rec = dict(outcome=o, cls=v.get('cls'), improvement=oc.get('improvement'),
               at_a_cost_beside=list(oc.get('at_a_cost_beside') or []),
               reading_dependent=bool(oc.get('reading_dependent', False)),
               verdict_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), verdict_path=str(path))
    if o in E135_TABLE_OUTCOMES:
        rung = 'R1' if E135_TABLE_OUTCOMES[o] == 'R1' else AL.step2_rung_of(AL.read_step1_fits())
        if rung is None:
            raise RuntimeError(f"E135's outcome ({o}) hands on the kept rung, but step 1 kept none above R1")
        return dict(rec, mode='table', rung=str(rung))
    if o in E135_NONE_OUTCOMES:
        return dict(rec, mode='search', rung=None)
    raise RuntimeError(f'A1 waits on E135: its outcome is {o!r} ({path})')


_A1_STAGES: Dict[Tuple[str, str, str, str], Tuple[pd.DataFrame, object, dict]] = {}


def a1_stage(wiring: str, panel: str = PANEL, orn: Optional[pd.DataFrame] = None) -> Tuple[pd.DataFrame, object, dict]:
    """A1's PN stage on a wiring and a panel's input when E135 hands one on (`require_e135` mode 'table'): E135's
    `kept_stage_table` of its step-1 fits on all seven PN flies at the handed-on rung, on the wiring's PN stage
    (`family_stage` on megamat, then the input put on it, as the table's docstring allows) with the panel's input
    (`orn`, default E134's `panel_input(panel)`, the input this module builds on). Cached per wiring, panel, rung and
    verdict. Returns (the table: the stage's glomeruli x the input's odor labels, Hz; the stage on the input, for
    `run_pn_check`; the record: E135's table record and the verdict's)."""
    from experiments import E134_population_baseline as B
    from experiments import E135_antennal_lobe_dendrites as AL
    h = require_e135()
    if h['mode'] != 'table':
        raise RuntimeError(f"E135's outcome ({h['outcome']}) hands on no stage: A1 searches the normalization")
    key = (wiring, panel, h['rung'], h['verdict_sha256'])
    if key not in _A1_STAGES:
        orn = B.panel_input(panel) if orn is None else orn
        stage, _ = AL.family_stage('megamat', wiring)
        fits = AL.read_step1_fits()['all_flies']
        table, rec = AL.kept_stage_table(fits, h['rung'], stage, orn=orn, panel=panel)
        _A1_STAGES[key] = (table, AL.stage_on_input(stage, orn), dict(rec, e135=h, wiring=wiring))
    return _A1_STAGES[key]


def pn_norm_kws(z: Mapping[str, float]) -> Dict[str, float]:
    """A1's fallback normalization at coordinates `z` (PN_COORDS): pn_inhadd = olfsysm's default x exp(z['inhadd']),
    pn_inhsc = E77's `inhsc_for(pn_inhadd)` x exp(z['inhsc']); {} (olfsysm's defaults) where both are 0."""
    a, g = float(z.get('inhadd', 0.0)), float(z.get('inhsc', 0.0))
    if a == 0.0 and g == 0.0:
        return {}
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    from experiments.E77_input_screen import inhsc_for
    inhadd = float(PN_INHADD_DEFAULT * np.exp(a))
    return dict(pn_inhadd=inhadd, pn_inhsc=float(inhsc_for(inhadd) * np.exp(g)))


def a1_kws(wiring: str, panel: str = PANEL, orn: Optional[pd.DataFrame] = None,
           z: Optional[Mapping[str, float]] = None) -> dict:
    """A1's kwargs for one wiring on a panel's input at coordinates `z`: mode 'table', dict(pn_rates=<the kept table on
    the input's own columns>) (E135's `input_columns` of `route1_pn_rates`, every entry set); mode 'search',
    `pn_norm_kws(z)`."""
    from experiments import E134_population_baseline as B
    from experiments import E135_antennal_lobe_dendrites as AL
    if require_e135()['mode'] == 'search':
        return pn_norm_kws(z or {})
    orn = B.panel_input(panel) if orn is None else orn
    table, _, _ = a1_stage(wiring, panel, orn)
    return dict(pn_rates=AL.input_columns(AL.route1_pn_rates(table, list(table.index)), orn))


def a1_pn_check(built, wiring: str, panel: str = PANEL) -> Optional[dict]:
    """A1 in mode 'table': the built wiring's PN rates after its last simulation against the kept table (E135's
    `run_pn_check`, route 1, strict: RuntimeError beyond ROUTE1_TOL_HZ); None in mode 'search'."""
    from experiments import E135_antennal_lobe_dendrites as AL
    if require_e135()['mode'] != 'table':
        return None
    table, stage, _ = a1_stage(wiring, panel, built.orn)
    return AL.run_pn_check(built, stage, table, list(table.index))


def stage_classes_on_measure(views: Sequence[str], scores: Optional[pd.DataFrame] = None,
                             group: str = START_GROUP) -> Dict[str, dict]:
    """E134's stage class (`classify_stages`' rule) on a measure of fewer views, from E134's committed score table
    (score_scores.csv): per stage of LOSS_STAGES and each of its first two readings, each of the eight wirings' error on
    the measure's views (the root mean square of its view_<view> columns) against the stage's most atypical fly on the
    same views (hal's `leave_one_out` with the measure's parts, all the stage's megamat flies); 'inside' on
    E134's INSIDE_MIN_WIRINGS of 8, 'far' above its FAR_ERROR on as many, else 'between'; 'not readable' where the two
    readings differ."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    sc = scores if scores is not None else pd.read_csv(e134_dir() / 'score_scores.csv')
    out = {}
    for st, fl in all_flies(LOSS_STAGES).items():
        fp, have = B.fly_side_read(fl)
        od = [o for o in have]
        loo = P.leave_one_out(fp, od, B.breadth_edges(len(od)), parts=measure_parts(views))
        yard = float(loo['error'].max())
        by = {}
        for rd in R.readings_of(st)[:2]:
            x = sc[(sc.group == group) & (sc.stage == st) & (sc.panel == PANEL) & (sc.condition == 'on')
                   & (sc.reading == rd.name) & sc.wiring.isin(WIRINGS)]
            if len(x) < len(WIRINGS):
                continue
            err = np.sqrt((x[[f'view_{v}' for v in views]].astype(float) ** 2).mean(axis=1)).to_numpy()
            inside, far = int((err <= yard).sum()), int((err > B.FAR_ERROR).sum())
            by[rd.name] = dict(cls='inside' if inside >= B.INSIDE_MIN_WIRINGS else 'far' if far >= B.INSIDE_MIN_WIRINGS
                               else 'between', n_inside=inside, n_far=far, max=float(err.max()))
        vals = [v['cls'] for v in by.values()]
        cls = (vals[0] if len(vals) == 2 and len(set(vals)) == 1 else 'not readable' if len(set(vals)) > 1
               else f'{vals[0]} (one reading)' if vals else 'not scored')
        out[st] = dict(by_reading=by, cls=cls, yardstick=yard, views=list(views))
    return out


def readable_stages(classes: Optional[Mapping] = None, group: str = START_GROUP) -> Tuple[str, ...]:
    """The stages of LOSS_STAGES that E134 found readable, in stage order: on every view (MEASURE 'all') from
    score_classes.json's stage classes of `group` ('not readable' left out); on a measure of fewer views from
    `stage_classes_on_measure` (written to <out_dir>/stage_classes.json), which must find all four readable (the re-read
    on the adopted measure; RuntimeError otherwise)."""
    if MEASURE != 'all':
        cl = stage_classes_on_measure(MEASURES[MEASURE], group=group)
        _write_json(out_dir() / 'stage_classes.json', dict(measure=MEASURE, group=group, classes=cl, time=_now(),
                                                            host=_host(), **_commits()))
        bad = [s for s in LOSS_STAGES if cl[s]['cls'] in ('not readable', 'not scored')]
        if bad:
            raise RuntimeError(f'on the measure {MEASURE!r} the stages {bad} are not readable: {cl}')
        return LOSS_STAGES
    c = classes if classes is not None else _read_json(e134_dir() / 'score_classes.json')
    st = ((c or {}).get('stages') or {}).get(group)
    if not st:
        raise RuntimeError(f"E134's score_classes.json has no stage classes for {group}")
    return tuple(s for s in LOSS_STAGES if s in st and st[s].get('cls') != 'not readable')


# ---------------------------------------------------------------------------
# the flies: halves, references and the permuted control

def all_flies(stages: Sequence[str] = LOSS_STAGES, panel: str = PANEL) -> Dict[str, object]:
    """Per stage, viki's recorded population of the panel's wild-type flies (APL on), as atlas's `fly_side` loads it
    (refuses before the lock is committed)."""
    out = {}
    for st in stages:
        fl = R.fly_side(st, panel, 'on')
        if fl is None:
            raise RuntimeError(f'no recorded {st} on {panel}')
        out[st] = fl
    return out


def split_names(flies, seed: int = HALF_SEED) -> Dict[str, List[str]]:
    """One stage's flies split as ruled: hal's `split_halves(n, 1, seed)`, its first part the fit half (floor(n / 2)
    flies) and the rest the held-out half, in the population's order; inside a `cell` with FIT_PART 1 the two parts
    are exchanged (the fit half is the second part)."""
    from kcq import population_fit as P
    names = [str(f) for f in flies.flies]
    a, b = P.split_halves(len(names), 1, seed)[0]
    if FIT_PART:
        a, b = b, a
    return dict(fit=[names[i] for i in a], held=[names[i] for i in b])


def fly_halves(seed: int = HALF_SEED, stages: Sequence[str] = LOSS_STAGES) -> Dict[str, Dict[str, List[str]]]:
    """Per stage of `stages`, the flies of the fit half and of the held-out half (`split_names`)."""
    return {st: split_names(fl, seed) for st, fl in all_flies(stages).items()}


def permutations(odors: Sequence[str], n: int = N_PERMS, seed: int = PERM_SEED) -> List[Dict[str, str]]:
    """The control's permutations of the odor labels: the first `n` draws of `np.random.default_rng(seed)`'s
    permutations of `odors` that leave no odor in place, each as {an odor: the label it takes}."""
    odors = list(odors)
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        p = rng.permutation(len(odors))
        if (p == np.arange(len(odors))).any():
            continue
        out.append({odors[i]: odors[int(p[i])] for i in range(len(odors))})
    return out


def permute_population(pop, perm: Mapping[str, str]):
    """hal's `Population` with its odor labels permuted (`perm`: an odor -> the label it takes), columns in their
    original order."""
    from kcq import population_fit as P
    cols = list(pop.values.columns)
    return P.Population(pop.values.rename(columns=dict(perm))[cols], pop.responds.rename(columns=dict(perm))[cols],
                        name=pop.name, named_units=pop.named_units, graded=pop.graded)


class Flies:
    """The flies a loss reads against: for reading (curve, noise, criterion) a FlyCache of one half of each stage's
    flies (all of them for `half` None); for scoring, per stage the references of the half (its mean views) in all the
    stage's flies' one-fly units, with the odor labels as recorded and, for the control, under each permutation; on the
    views of `views` ('measure': the measure in force, MEASURES; None every view).

    Attributes: half, stages, panel, views, full (per stage, viki's population of all flies), names (per stage, the
    halves' fly names), perms (the control's permutations, or []), cache (the FlyCache)."""

    def __init__(self, half: Optional[str], stages: Sequence[str] = LOSS_STAGES,
                 perms: Optional[Sequence[Mapping[str, str]]] = None, panel: str = PANEL, seed: int = HALF_SEED,
                 full: Optional[Mapping[str, object]] = None, views: Optional[Sequence[str]] = 'measure'):
        if half is not None and half not in HALVES:
            raise ValueError(f'{half=}')
        self.half, self.stages, self.panel = half, tuple(stages), panel
        self.views = MEASURES[MEASURE] if views == 'measure' else (None if views is None else tuple(views))
        self.full = dict(full) if full is not None else all_flies(stages, panel)
        self.names = {st: split_names(fl, seed) for st, fl in self.full.items()}
        self.perms = list(perms or [])
        fixed = {}
        for st, fl in self.full.items():
            fixed[(st, panel, 'on')] = fl if half is None else R.subset_flies(fl, self.names[st][half])
            fixed[(st, panel, 'off')] = None
        self.cache = R.FlyCache(fixed=fixed, fallback=False)
        self._refs: Dict[tuple, object] = {}

    def reference(self, stage: str, odors: Sequence[str], perm: Optional[int] = None):
        """The reference of `stage` on the model's `odors` (those the flies have; E134's `fly_side_read`): all flies'
        parts, kept to the half's flies with all flies' one-fly units; `perm` the index of a permutation, or None for
        the recorded labels. Cached."""
        from experiments import E134_population_baseline as B
        from kcq import population_fit as P
        key = (stage, tuple(odors), perm)
        if key not in self._refs:
            fp, have = B.fly_side_read(self.full[stage])
            od = [o for o in odors if o in set(have)]
            if perm is not None:
                fp = [permute_population(f, self.perms[perm]) for f in fp]
            ref = P.fly_reference(fp, od, B.breadth_edges(len(od)), parts=measure_parts(self.views))
            if self.half is not None:
                ref = P.subset_reference(ref, self.names[stage][self.half], keep_units=True)
            self._refs[key] = ref
        return self._refs[key]


# ---------------------------------------------------------------------------
# the claws: A2's divisor and A3's threshold

def claw_inputs_divisive(W: np.ndarray, boutons: np.ndarray, inh: np.ndarray, wAPLKC: np.ndarray) -> np.ndarray:
    """The claws' net drive as olfsysm's sim_KC_layer forms it with `apl_inh_divisive` (E92's form; no claw dynamics):
    wiring x bouton rate at step t divided by 1 plus the claw's scaled APL weight times APL's activity at step t - 1
    (that product clipped at 0), then clipped at 0 (claws x steps). Arguments as E90's `claw_inputs` (atlas's
    `claw_reader` signature)."""
    x = W @ boutons
    i = np.outer(np.asarray(wAPLKC, dtype=float).ravel(), np.asarray(inh, dtype=float).ravel())
    np.maximum(i, 0.0, out=i)
    x /= 1.0 + i
    np.maximum(x, 0.0, out=x)
    return x


def claw_reader_after_threshold(claw_thr: float, divisive: bool = False) -> Callable:
    """A claw reader of the drive after a per-claw threshold (A3's reading beside): max(drive - claw_thr, 0), the
    drive E90's `claw_inputs` (or `claw_inputs_divisive`)."""
    from experiments import E90_model_stages as E90
    base = claw_inputs_divisive if divisive else E90.claw_inputs

    def reader(W, boutons, inh, wAPLKC):
        x = base(W, boutons, inh, wAPLKC)
        x -= float(claw_thr)
        np.maximum(x, 0.0, out=x)
        return x
    reader.__name__ = f'claw_after_threshold_{claw_thr:g}'
    return reader


def claw_reader_of(arm: str) -> Optional[Callable]:
    """The claw reader of an arm's runs (atlas's `run_stages` `claw_reader`): `claw_inputs_divisive` where the arm's
    kwargs make APL's inhibition of the claws divisive, else None (E90's subtractive `claw_inputs`)."""
    return claw_inputs_divisive if ARMS[base_arm(arm)].kws.get('apl_inh_divisive') else None


def claw_identity(arm: str = 'A2', odors: Sequence[str] = ('eb', '1-8ol'), wiring: str = 'real') -> dict:
    """The claw reading's identity with olfsysm (E90's `claw_identity` for an arm's form): `wiring` at its held megamat
    point with the arm's kwargs, the megamat input restricted to `odors` (the breadth held at the lab input's), built
    with olfsysm's dynamics recorded, APL on; the claw inputs formed by the arm's claw reader from its bouton_sims,
    inh_sims and APL->claw weights must equal its claw_sims at every step after the first."""
    from experiments import E10_joint_fit as J
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E90_model_stages as E90
    from experiments import E134_population_baseline as B
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator, get_sim_var
    Pp = E83.wiring_parts(wiring)
    orn = B.panel_input(PANEL)
    cols = [c for c in orn.columns if H.odor_abbrev(c) in odors]
    p, _ = start_params(wiring)
    wnorm, W = Pp['wnorm'], Pp['W']
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], Pp['breadth'], E83.BREADTH_K)
    sim, _ = make_fixed_param_kc_simulator(orn[cols], initial_params=p, **dict(J.base_kws(None), **arm_kws(arm)),
                                           _wPNKC=W, _wAPLKC=wAPLKC, _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'],
                                           _wPNAPL=wnorm['wPNAPL'], kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU,
                                           return_dynamics=True, delete_pretime=False)
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    rv = cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    sim(p)
    reader = claw_reader_of(arm) or E90.claw_inputs
    claws = get_sim_var(rv.kc, 'claw_sims')
    btn = get_sim_var(rv.pn, 'bouton_sims')
    inh = get_sim_var(rv.kc, 'inh_sims')
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    Wm = W.to_numpy(dtype=float)
    nd, mx, nsteps = 0, 0.0, 0
    for j in range(len(order)):
        i = int(order[j])
        b = np.asarray(btn[i])
        ok = np.flatnonzero(~np.isnan(b).any(axis=0))
        ok = ok[ok >= 1]
        x = reader(Wm, b[:, ok], np.asarray(inh[i]).ravel()[ok - 1], wk)
        c = np.asarray(claws[i])[:, ok]
        nd += int((x != c).sum())
        mx = max(mx, float(np.max(np.abs(x - c))))
        nsteps += int(len(ok))
    res = dict(arm=arm, wiring=wiring, odors=[H.odor_abbrev(c) for c in cols], reader=reader.__name__,
               n_claws=int(Wm.shape[0]), steps_compared=nsteps, claw_values_differ=nd, claw_max_abs_diff=mx,
               exact=bool(nd == 0))
    del sim, cells, rv, claws, btn
    return res


def _outs_differ(a: Mapping[str, dict], b: Mapping[str, dict]) -> Dict[str, int]:
    """Per condition and frame of two `run_stages` outputs, how many entries differ (-1 where shapes or indices do)."""
    from experiments import E134_population_baseline as B
    d = {}
    for c in R.CONDITIONS:
        for k in sorted(set(a[c]) | set(b[c])):
            if k not in a[c] or k not in b[c]:
                d[f'{c}/{k}'] = -1
                continue
            d[f'{c}/{k}'] = B._differ(a[c][k], b[c][k])
    return d


def claw_thr_check(wirings: Optional[Sequence[str]] = None, value: float = CLAW_THR_CHECK_HZ) -> dict:
    """A3's mechanics: on each wiring at its held megamat point, a run with the per-claw threshold set in place on the
    simulator built without it (olfsysm's mp.kc.claw_thr) equals a fresh build at that threshold exactly, counts and
    every stage; written to claw_thr_check.json ('exact' True only if every wiring's is)."""
    from experiments import E134_population_baseline as B
    wirings = S_WIRINGS if wirings is None else wirings
    rec = dict(value=value, wirings={}, time=_now(), host=_host())
    for w in wirings:
        p, poff = start_params(w)
        a = B.build_family(w, PANEL, dict(p))
        a.vars()[0].kc.claw_thr = float(value)
        oa = B.run_stages(B.Point(PANEL, 'held', w), a, p, poff)
        b = B.build_family(w, PANEL, dict(p), extra_kws={'claw_thr': float(value)})
        ob = B.run_stages(B.Point(PANEL, 'held', w), b, p, poff)
        d = _outs_differ(oa, ob)
        rec['wirings'][w] = dict(differ=d, exact=all(v == 0 for v in d.values()),
                                 n_on=n_responses(oa['on']['responses']))
        del a, b
    rec['exact'] = all(v['exact'] for v in rec['wirings'].values())
    rec.update(_commits())
    _write_json(out_dir() / 'claw_thr_check.json', rec)
    return rec


# ---------------------------------------------------------------------------
# the model: an arm on several wirings, its rates held

def base_arm(arm: str) -> str:
    """The arm an arm name, a control (CONTROLS: CONTROL as A0, CONTROL_A1 as A1) or the unrefitted family (FAMILY) is
    built as."""
    if arm in CONTROLS:
        return CONTROLS[arm]
    return 'A0' if arm == FAMILY else arm


# an arm's kwargs that set its APL-off condition only (`off_extra_of`), not model kwargs
OFF_ONLY_KWS = ('kckc_off_factor',)


def arm_kws(arm: str) -> dict:
    """The arm's fixed model kwargs beside the family's (make_fixed_param_kc_simulator kwargs; A4's KC-KC scale is not
    a kwarg and is left out, nor are the APL-off condition's own settings, OFF_ONLY_KWS, and A1's PN stage depends on
    the wiring, the panel and E135: `a1_kws`)."""
    return {k: v for k, v in ARMS[base_arm(arm)].kws.items() if k not in ('pn_rates', 'kckc_fixed', *OFF_ONLY_KWS)}


def off_extra_of(arm: str, p: Mapping[str, float]) -> Optional[Dict[str, float]]:
    """The APL-off condition's own parameters for an arm at APL-on parameters `p` (E77's `hold_rates` `off_extra`):
    with 'kckc_off_factor' in the arm's kwargs, the KC-KC scale at that factor of `p`'s (E141's inhibition among KCs
    present with APL on only); None for every other arm (the APL-off condition is the APL-on one with APL's scales 0)."""
    f = ARMS[base_arm(arm)].kws.get('kckc_off_factor')
    return None if f is None else dict(wKCKC_scale=float(p['wKCKC_scale']) * float(f))


def arm_searched(arm: str) -> Tuple[str, ...]:
    """An arm's searched coordinates: ARMS', and for A1 PN_COORDS beside them when E135 hands on no stage
    (`require_e135` mode 'search'; RuntimeError while A1 waits on E135)."""
    a = ARMS[base_arm(arm)]
    if a.name == 'A1' and require_e135()['mode'] == 'search':
        return tuple(a.searched) + PN_COORDS
    return tuple(a.searched)


def coord_value(k: str, v: float) -> float:
    """A coordinate's value as reported: the factor exp(v) for the log factors and PN_COORDS, the per-claw threshold in
    Hz."""
    return CLAW_THR_UNIT * float(v) if k == 'claw_thr' else float(np.exp(v))


def start_params(wiring: str, panel: str = PANEL) -> Tuple[Dict[str, float], Dict[str, float]]:
    """A wiring's starting APL-on and APL-off absolute parameters: its held megamat point of E134's step 3 (the run's
    saved params and params_off)."""
    from experiments import E134_population_baseline as B
    run = B.load_run(B.tag_of(B.Point(panel, 'held', wiring)))
    if run is None:
        raise FileNotFoundError(f"E134's {panel}_held_{wiring} is not saved")
    return ({k: float(v) for k, v in run.info['params'].items()},
            {k: float(v) for k, v in run.info['params_off'].items()})


def n_responses(responses: pd.DataFrame) -> int:
    """The number of responses (two or more spikes) of a responses table."""
    return int((responses.to_numpy(dtype=float) > 0.5).sum())


def secant_hold(run: Callable[[float], tuple], target: int, x0: float, lo: float, hi: float, slope0: float,
                max_sims: int = HOLD_MAX_SIMS, tol: float = HOLD_TOL) -> Tuple[float, tuple, dict]:
    """The log factor x in [lo, hi] at which `run(x)`'s number of responses (falling with x) is within `tol` of
    `target`: from `x0`, Newton steps on ln(responses) with the slope `slope0` (or the secant slope of the last two
    points where it falls) until the target is bracketed, then the secant of the bracket's two ends kept within its
    middle 90 %, or the bracket's midpoint after two runs in a row on the same side of the target (so that a curved
    count does not stall the secant); at most `max_sims` runs.

    Returns (x, outputs, info): the run nearest the target; info has 'n', 'target', 'n_sims', 'ok' (within `tol`),
    'at_end' (that run at an end of the range without reaching the target), 'miss' (|n - target| / target), 'slope'
    (the last slope estimate), 'x' and 'history'."""
    lt = float(np.log(max(target, 1)))
    pts: List[Tuple[float, int, tuple]] = []
    slope = float(slope0)
    x = float(min(max(x0, lo), hi))
    sides: List[bool] = []
    for _ in range(max_sims):
        out = run(x)
        n = n_responses(out[0])
        pts.append((x, n, out))
        if abs(n - target) <= tol * target:
            break
        above = [q for q in pts if q[1] > target]
        below = [q for q in pts if q[1] < target]
        if above and below:
            sides.append(n > target)
            a = max(above, key=lambda q: q[0])
            b = min(below, key=lambda q: q[0])
            if a[0] < b[0] and len(sides) >= 3 and sides[-1] == sides[-2]:
                xn = 0.5 * (a[0] + b[0])
            elif a[0] < b[0]:
                ya, yb = np.log(max(a[1], 1)), np.log(max(b[1], 1))
                if yb != ya:
                    slope = float((yb - ya) / (b[0] - a[0]))
                xn = a[0] + (lt - ya) / slope if slope < 0 else 0.5 * (a[0] + b[0])
                m = 0.05 * (b[0] - a[0])
                xn = min(max(xn, a[0] + m), b[0] - m)
            else:
                xn = 0.5 * (a[0] + b[0])
        else:
            if len(pts) >= 2:
                (x1, n1, _), (x2, n2, _) = pts[-2], pts[-1]
                if x2 != x1 and n1 > 0 and n2 > 0:
                    s = (np.log(n2) - np.log(n1)) / (x2 - x1)
                    if np.isfinite(s) and s < 0:
                        slope = float(s)
            xn = x + (lt - np.log(max(n, 1))) / slope
            xn = float(min(max(xn, lo), hi))
            if xn == x:
                break
        x = float(xn)
    best = min(pts, key=lambda q: abs(q[1] - target))
    ok = abs(best[1] - target) <= tol * target
    info = dict(n=int(best[1]), target=int(target), n_sims=len(pts), ok=bool(ok), slope=slope,
                at_end=bool((not ok) and (best[0] <= lo + 1e-12 or best[0] >= hi - 1e-12)),
                miss=float(abs(best[1] - target) / max(target, 1)), x=float(best[0]),
                history=[(float(q[0]), int(q[1])) for q in pts])
    return float(best[0]), best[2], info


class ArmModel:
    """One arm of the family on a set of wirings (the megamat panel): each wiring built once (E134's `build_family`,
    with the arm's kwargs, at its held megamat point), run at shared searched coordinates `z` with both rates held per
    wiring (`hold`) and every stage saved (E134's `run_stages`, the arm's claw reader).

    Attributes: arm, panel, wirings, searched (the arm's searched coordinates), base / base_off (per wiring, the start's
    absolute parameters of each condition), fixed (parameters the arm holds at a value, A4's KC-KC 0), kws (the arm's
    model kwargs), claw_reader, claw_mode (A3: 'in_place' or 'rebuild'), targets (per wiring, (APL-off, APL-on)
    response counts), built (per wiring, E134's `Built`), memory (per wiring, the evaluated points' coordinates and
    held log factors), slopes (per wiring, the holds' last slopes)."""

    def __init__(self, arm: str, wirings: Optional[Sequence[str]] = None, claw_mode: str = 'in_place',
                 lock: Optional[Mapping] = None):
        from experiments import E83_cross_wiring_fit as E83
        from experiments import E134_population_baseline as B
        self.arm, self.panel, self.wirings = arm, PANEL, tuple(S_WIRINGS if wirings is None else wirings)
        a = ARMS[base_arm(arm)]
        self.searched = arm_searched(arm)
        if len(self.searched) + len(HELD) > MAX_FREE:
            raise ValueError(f'{arm}: more than {MAX_FREE} parameters')
        self.fixed = {'wKCKC_scale': float(a.kws['kckc_fixed'])} if 'kckc_fixed' in a.kws else {}
        if any(k in a.kws for k in OFF_ONLY_KWS):
            raise NotImplementedError(f"{arm}: a fit of an arm with an APL-off condition of its own "
                                      f"({[k for k in OFF_ONLY_KWS if k in a.kws]}) is not written; reads only")
        self.kws = arm_kws(arm)
        self.pn_mode = require_e135()['mode'] if a.name == 'A1' else None
        self.pn_search = self.pn_mode == 'search'
        # per wiring, the kwargs every build of it has (A1's table in mode 'table'), and those its current build has
        self.wkws = {w: dict(self.kws, **(a1_kws(w, PANEL) if self.pn_mode == 'table' else {})) for w in self.wirings}
        self._built_kws: Dict[str, dict] = {}
        self.pn_checks: Dict[str, dict] = {}
        self.claw_reader = claw_reader_of(arm)
        self.claw_mode = claw_mode if 'claw_thr' in self.searched else None
        if self.claw_mode not in (None, 'in_place', 'rebuild'):
            raise ValueError(f'{claw_mode=}')
        self.apl_range = HOLD_RANGE_DIVISIVE_APL if self.kws.get('apl_inh_divisive') else HOLD_RANGE
        self.base, self.base_off, self.built, self.targets = {}, {}, {}, {}
        self._claw_thr: Dict[str, float] = {}
        lock = lock if lock is not None else R.require_lock()
        for w in self.wirings:
            p, poff = start_params(w)
            p.update(self.fixed)
            poff.update(self.fixed)
            self.base[w], self.base_off[w] = p, poff
            self.built[w] = B.build_family(w, PANEL, dict(p), extra_kws=self.wkws[w] or None)
            self.built[w].info.update(arm=arm)
            self._built_kws[w] = dict(self.wkws[w])
            self._claw_thr[w] = 0.0
            n_kc = E83.saved_fit(w)['saved']['on'][0].shape[0]
            t_off, t_on, _, _ = B.held_targets(B.Point(PANEL, 'held', w), n_kc * self.built[w].orn.shape[1], lock)
            self.targets[w] = (int(t_off), int(t_on))
        self.memory: Dict[str, List[Tuple[Tuple[float, ...], float, float]]] = {w: [] for w in self.wirings}
        self.slopes: Dict[str, Dict[str, float]] = {w: dict(HOLD_SLOPE0) for w in self.wirings}

    def coords(self, z: Mapping[str, float]) -> Tuple[float, ...]:
        return tuple(float(z.get(k, 0.0)) for k in self.searched)

    def params_for(self, wiring: str, z: Mapping[str, float], x_thr: float = 0.0, x_apl: float = 0.0
                   ) -> Tuple[Dict[str, float], Dict[str, float], Optional[float]]:
        """The APL-on and APL-off absolute parameters of a wiring at searched coordinates `z` (0 is the start) and held
        log factors `x_thr`, `x_apl`, and the per-claw threshold (Hz; None where the arm has none): every parameter a
        factor names is multiplied by exp(factor), so the KC->APL and bouton->APL scales keep the wiring's own ratios to
        the APL->KC scale; APL off keeps APL->KC and APL->bouton at 0 and takes the other parameters from APL on's."""
        bad = set(z) - set(self.searched)
        if bad:
            raise ValueError(f'{self.arm}: {sorted(bad)} are not searched')
        p, poff = dict(self.base[wiring]), dict(self.base_off[wiring])
        lf = {k: float(z[k]) for k in z if k in FACTORS}
        lf.update(thr=float(x_thr), apl=float(x_apl))
        for k, v in lf.items():
            f = float(np.exp(v))
            for kk in FACTORS[k]:
                p[kk] = p[kk] * f
                if kk not in OFF_ZERO:
                    poff[kk] = poff[kk] * f
        for kk in p:
            if kk not in OFF_ZERO and kk not in FACTORS['apl']:
                poff[kk] = p[kk]
        ct = CLAW_THR_UNIT * float(z['claw_thr']) if 'claw_thr' in z else (0.0 if 'claw_thr' in self.searched else None)
        return p, poff, ct

    def set_claw_thr(self, wiring: str, value: Optional[float]) -> None:
        """A3: the per-claw threshold of the wiring's simulator, in place on olfsysm's parameters or by a rebuild."""
        from experiments import E134_population_baseline as B
        if value is None or self.claw_mode is None or float(value) == self._claw_thr[wiring]:
            return
        if self.claw_mode == 'in_place':
            self.built[wiring].vars()[0].kc.claw_thr = float(value)
        else:
            kws = dict(self.wkws[wiring], **({'claw_thr': float(value)} if value > 0 else {}))
            self._rebuild(wiring, kws)
        self._claw_thr[wiring] = float(value)

    def _rebuild(self, wiring: str, kws: Mapping) -> None:
        """The wiring's simulator built again with `kws` (E134's `build_family` at the start's parameters), the old one
        dropped and collected first. pyolfsysm's fixed-parameter simulator is a reference cycle (its closure refers to
        itself), so a dropped one, with olfsysm's variables behind it (about 0.4 GB a wiring on megamat), is freed only by
        the cycle collector, which counts Python objects and not that memory: without the collection A1's search, which
        rebuilds every wiring at every evaluation, grew by about 1 GB per evaluation (2026-11-27, disco's probe)."""
        from experiments import E134_population_baseline as B
        self.built.pop(wiring, None)
        gc.collect()
        self.built[wiring] = B.build_family(wiring, PANEL, dict(self.base[wiring]), extra_kws=dict(kws) or None)
        self._built_kws[wiring] = dict(kws)

    def set_pn(self, wiring: str, z: Mapping[str, float]) -> None:
        """A1 when E135 hands on no stage: the wiring's simulator rebuilt with the normalization at `z`
        (`pn_norm_kws`) where its current build has other kwargs (`_rebuild`)."""
        if not self.pn_search:
            return
        kws = dict(self.wkws[wiring], **pn_norm_kws(z))
        if kws == self._built_kws[wiring]:
            return
        self._rebuild(wiring, kws)
        self.built[wiring].info.update(arm=self.arm)

    def hold(self, wiring: str, z: Mapping[str, float]) -> Tuple[Dict[str, float], Dict[str, float], dict]:
        """Both response counts of a wiring held at `z` (`secant_hold`, warm-started from the evaluated point nearest
        `z`, at most HOLD_MAX_SIMS simulations per rate, HOLD_FIRST_MAX_SIMS for the wiring's first hold): the APL-off
        count by the threshold, then the APL-on count by the APL->KC scale with its tied scales.
        Returns (params, params_off, info: x_thr, x_apl, each hold's info, the held outputs)."""
        c = self.coords(z)
        mem = self.memory[wiring]
        if mem:
            _, xt0, xa0 = min(mem, key=lambda m: sum((a - b) ** 2 for a, b in zip(m[0], c)))
        else:
            xt0 = xa0 = 0.0
        sim = self.built[wiring].sim
        t_off, t_on = self.targets[wiring]
        lt, ht = np.log(HOLD_RANGE[0]), np.log(HOLD_RANGE[1])
        la, ha = np.log(self.apl_range[0]), np.log(self.apl_range[1])
        n_max = HOLD_MAX_SIMS if mem else HOLD_FIRST_MAX_SIMS
        xt, out_off, h_off = secant_hold(lambda x: sim(self.params_for(wiring, z, x, 0.0)[1]), t_off, xt0, lt, ht,
                                         self.slopes[wiring]['thr'], n_max)
        xa, out_on, h_on = secant_hold(lambda x: sim(self.params_for(wiring, z, xt, x)[0]), t_on, xa0, la, ha,
                                       self.slopes[wiring]['apl'], n_max)
        self.slopes[wiring] = dict(thr=h_off['slope'], apl=h_on['slope'])
        mem.append((c, xt, xa))
        p, poff, _ = self.params_for(wiring, z, xt, xa)
        return p, poff, dict(x_thr=xt, x_apl=xa, off=h_off, on=h_on, outputs={'on': out_on, 'off': out_off})

    def run(self, z: Mapping[str, float], tag_prefix: str = '') -> Dict[str, 'R.SavedRun']:
        """Every wiring at `z`, both rates held, both conditions run again with every stage saved (which must give the
        hold's outputs exactly), as in-memory runs (`in_memory_run`; tags '<tag_prefix><wiring>'; the hold in
        run.info['hold'])."""
        from experiments import E134_population_baseline as B
        runs = {}
        for w in self.wirings:
            self.set_pn(w, z)
            self.set_claw_thr(w, self.params_for(w, z)[2])
            p, poff, hd = self.hold(w, z)
            ct = self.params_for(w, z)[2]
            out = B.run_stages(B.Point(PANEL, 'held', w), self.built[w], p, poff, claw_reader=self.claw_reader)
            for c in R.CONDITIONS:
                r, s = hd['outputs'][c]
                if B._differ(out[c]['spike_counts'], s) or B._differ(out[c]['responses'], r):
                    raise RuntimeError(f"{self.arm} {w} {c}: the saved run differs from the hold's outputs")
            info = dict(model=B.FAMILY, wiring=w, panel=PANEL, reading='fit', arm=self.arm, z=dict(z), params=p,
                        params_off=poff, claw_thr=ct, tag=f'{tag_prefix}{w}')
            if self.pn_mode == 'table' and w not in self.pn_checks:
                self.pn_checks[w] = a1_pn_check(self.built[w], w)
            if self.pn_mode is not None:
                info.update(pn_mode=self.pn_mode, pn_kws={k: v for k, v in self._built_kws[w].items() if k != 'pn_rates'
                                                          and k.startswith('pn_')}, pn_check=self.pn_checks.get(w))
            runs[w] = in_memory_run(f'{tag_prefix}{w}', out, self.built[w], info)
            runs[w].info['_out'] = out
            runs[w].info['hold'] = {k: v for k, v in hd.items() if k != 'outputs'}
        return runs


def in_memory_run(tag: str, out: Mapping[str, dict], built, info: Mapping[str, object]) -> 'R.SavedRun':
    """E134's `run_stages` output as atlas's `SavedRun` without writing it: the layout `load_step3_run` gives after
    `save_step3_run` (stage matrices with odors abbreviated, the settled pre-odor means in place of E90's window with
    E90's kept as '<stage>_pre_e90', the responses and spike counts as run), checked by `R.check_meta`; the test checks
    a saved and reloaded run reads the same."""
    conds = {}
    for c in R.CONDITIONS:
        st = {}
        for k, x in out[c].items():
            if k == 'apl':
                st[k] = x.copy()
            elif k in ('responses', 'spike_counts'):
                st[k] = x
            elif k.split('_')[0] in R.STAGE_FILES and k.endswith(('_win', '_pre', '_pre_settled')):
                st[k] = R._abbrev(x)
        for k in R.STAGE_FILES:
            if f'{k}_pre_settled' in st:
                st[f'{k}_pre_e90'] = st[f'{k}_pre']
                st[f'{k}_pre'] = st.pop(f'{k}_pre_settled')
            elif f'{k}_pre' in st:
                raise RuntimeError(f'{tag} {c} {k}: no settled pre-odor means')
        conds[c] = st
    run = R.SavedRun(source='E136', tag=tag, model=str(info['model']), wiring=str(info['wiring']),
                     input_panel=str(info['panel']), conditions=conds,
                     meta={k: np.asarray(v) for k, v in built.meta.items()}, orn=R._abbrev(built.orn), files=[],
                     settled=True, info=dict(info))
    R.check_meta(run)
    return run


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E136_stage_population_fit/[measure_<name>/]; PATHS['archive'] inside `paths`."""
    if 'archive' in PATHS:
        return PATHS['archive']
    d = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return d / _measure_sub() if _measure_sub() else d


def fit_dirs(arm: str) -> List[Path]:
    """Where a fit's saved evaluations are: outputs/.../fit_<arm>/ on the host that ran it (every evaluation), then the
    archive's fit_<arm>/ (the best evaluation, `archive_best`)."""
    return [local_dir(f'fit_{arm}'), archive_dir() / f'fit_{arm}']


def archive_best(arm: str, i: int) -> List[Path]:
    """Copies a fit's best evaluation (every S wiring's files, '<arm>_e<i>_<wiring>_*') from this host's
    outputs/.../fit_<arm>/ to the archive's fit_<arm>/, so the score can read it on any host. A wiring whose evaluation
    this host did not run (a fit resumed from another host's history) is taken as it is if the archive already holds it
    complete (its run.json), as that host's `--part archive` copies it."""
    import shutil
    src, dst = local_dir(f'fit_{arm}'), archive_dir() / f'fit_{arm}'
    dst.mkdir(parents=True, exist_ok=True)
    out = []
    for w in S_WIRINGS:
        files = sorted(src.glob(f'{arm}_e{i}_{w}_*'))
        if not any(f.name.endswith('_run.json') for f in files):
            have = sorted(dst.glob(f'{arm}_e{i}_{w}_*'))
            if any(f.name.endswith('_run.json') for f in have):
                out += have
                continue
            raise FileNotFoundError(f'{arm} eval {i} {w}: not saved under {src} nor in the archive {dst} (run '
                                    f'`--part archive --arm {arm} --evals {i}` on the host that ran it)')
        for f in sorted(files, key=lambda f: f.name.endswith('_run.json')):   # the run.json last, as it says complete
            tmp = dst / f'.{f.name}.{_host()}.{os.getpid()}.tmp'
            shutil.copyfile(f, tmp)
            os.replace(tmp, dst / f.name)
            out.append(dst / f.name)
    return out


def archive_evals(arm: str, evals: Optional[Sequence[int]] = None) -> List[int]:
    """Copies evaluations this host saved for a fit to the archive (`archive_best`), so that a fit resumed on another
    host from the history can archive its best; `evals` by default the best of the evaluations this host ran (by the
    history's 'host'). Returns the evaluations copied."""
    hist = pd.read_csv(out_dir() / f'fit_history_{arm}.csv', float_precision='round_trip')
    if evals is None:
        mine = hist[hist['host'].astype(str) == _host()]
        if mine.empty:
            raise ValueError(f'{arm}: no evaluation in the history was run on {_host()}')
        evals = [int(mine.loc[mine['loss_sim'].idxmin(), 'eval'])]
    for i in evals:
        archive_best(arm, int(i))
        print(f'E136 {arm}: eval {int(i)} archived', flush=True)
    return [int(i) for i in evals]


# the absolute parameters of a run, both conditions (APL off: 'off_<name>'), as E135's variants.csv names them
VARIANT_PARAMS = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale')


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv in the measure's out_dir, for the model-output archive's params.csv: one row per run the archive
    holds (`arch`, default archive_dir(): each fit's archived evaluations in fit_<arm>/, every read in read/), from its
    <tag>_run.json. A row has the tag (the files' variant name), the folder, arm, panel, wiring, the reading ('fit' or
    'read'), the evaluation (a fit's), a description, the searched coordinates (z_<name>), both conditions' absolute
    parameters (VARIANT_PARAMS; 'off_<name>' with APL off), A3's per-claw threshold (Hz), A1's PN normalization terms
    (pn_inhadd, pn_inhsc), and where and when it ran. The runs were written into the archive from every host that ran
    them, so the archive, not this host's outputs, is the list."""
    arch = archive_dir() if arch is None else Path(arch)
    rows = []
    for p in sorted(arch.glob('*/*_run.json')):
        d = json.loads(p.read_text())
        arm = str(d.get('arm'))
        what = 'the unrefitted family' if arm == FAMILY else f'arm {arm}' + (' (the control, A0 on odor-permuted flies)'
                                                                             if arm == CONTROL else '')
        reading = d.get('reading')
        kind = (f"{what}, fitted to the megamat stages: evaluation {d.get('eval')} on searched wiring {d.get('wiring')}"
                if reading == 'fit' else f"{what} at its classed point, read on the {d.get('panel')} input, wiring "
                                          f"{d.get('wiring')}, both response rates held")
        r = dict(tag=d['tag'], variant=d['tag'], part=p.parent.name, subdir=p.parent.name, arm=arm, panel=d.get('panel'),
                 wiring=d.get('wiring'), reading=reading, eval=d.get('eval'), model=d.get('model'), kind=kind)
        r.update({f'z_{k}': v for k, v in (d.get('z') or {}).items()})
        r.update({k: (d.get('params') or {}).get(k) for k in VARIANT_PARAMS})
        r.update({f'off_{k}': (d.get('params_off') or {}).get(k) for k in VARIANT_PARAMS})
        r['claw_thr'] = d.get('claw_thr')
        pn = d.get('pn_kws') if isinstance(d.get('pn_kws'), dict) else {}
        r['pn_inhadd'] = d.get('pn_inhadd', pn.get('pn_inhadd'))
        r['pn_inhsc'] = d.get('pn_inhsc', pn.get('pn_inhsc'))
        r.update(host=d.get('host'), lane=d.get('lane'), time=d.get('time'))
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        raise FileNotFoundError(f'no <tag>_run.json under {arch}')
    if df.tag.duplicated().any():
        raise ValueError(f'tags in two folders: {sorted(df.tag[df.tag.duplicated()])}')
    out = out_dir() / 'variants.csv'
    df.to_csv(out, index=False)
    print(f'variants.csv: {len(df)} runs -> {out}', flush=True)
    return df


def save_eval(arm: str, i, runs: Mapping[str, 'R.SavedRun']) -> List[Path]:
    """An evaluation's runs under outputs/.../fit_<arm>/ (E134's step-3 layout, `save_step3_run`; tag
    '<arm>_e<i>_<wiring>' for an evaluation's index, '<arm>_<i>_<wiring>' for a label, e.g. the single-factor
    returns' 'det_<factor>0')."""
    d = local_dir(f'fit_{arm}')
    written = []
    for w, run in runs.items():
        tag = f'{arm}_e{i}_{w}' if isinstance(i, (int, np.integer)) else f'{arm}_{i}_{w}'
        info = {k: v for k, v in run.info.items() if not str(k).startswith('_')}
        info.update(tag=tag, eval=i, host=_host(), time=_now())
        written += R.save_step3_run(tag, run.info['_out'], run.meta, run.orn, info, [d])
    return written


# ---------------------------------------------------------------------------
# the loss

def rates(runs: Mapping[str, 'R.SavedRun']) -> Dict[str, float]:
    """The wirings' mean response rates with APL on and off (two or more spikes), their ratio, and each wiring's."""
    out = {}
    for c in R.CONDITIONS:
        each = {w: float((run.conditions[c]['responses'].to_numpy(dtype=float) > 0.5).mean()) for w, run in runs.items()}
        out[f'rate_{c}'] = float(np.mean(list(each.values())))
        out.update({f'rate_{c}_{w}': v for w, v in each.items()})
    out['rate_ratio'] = out['rate_off'] / out['rate_on'] if out['rate_on'] > 0 else np.nan
    return out


def hold_penalty(runs: Mapping[str, 'R.SavedRun']) -> float:
    """The holds' penalty: one one-fly unit per HOLD_MISS_PER of relative miss of every hold that did not reach its
    target (0 for runs without a hold)."""
    pen = 0.0
    for run in runs.values():
        for c in ('off', 'on'):
            h = (run.info.get('hold') or {}).get(c)
            if h and not h['ok']:
                pen += h['miss'] / HOLD_MISS_PER
    return float(pen)


SPLIT_PARTS = ('strength/rate', 'strength/mean')
SPLIT_TERMS = ('level', 'spread', 'order', 'sd_ratio', 'r')


def strength_split(pops: Sequence['P.Population'], ref: 'P.FlyReference') -> pd.DataFrame:
    """The strength parts' errors split into level, spread and order (the amplitude and shape split of 2026-11-13, on
    strength by odor; viki's request of 2026-11-24 22:54): for each part of SPLIT_PARTS the reference has, with d the
    model's entries minus the flies' mean over the entries both define, sd the standard deviation over them (ddof 0)
    and r the Pearson correlation of the model's entries with the flies' mean, 'level' is |mean d|, 'spread' |sd_model -
    sd_flies| and 'order' sqrt(2 sd_model sd_flies (1 - r)), each in the part's one-fly unit, so that level^2 +
    spread^2 + order^2 is the part's error^2 (`P.score_arrays`); 'sd_ratio' is sd_model / sd_flies and 'r' the
    correlation. The strength parts are multiples of the panel's mean, so the level is near 0. Rows as hal's
    `wirings_as_flies`: one per wiring, 'mean' (their mean) and 'common' (the wirings' mean views); columns
    '<view>_<part>_<term>'."""
    from kcq import population_fit as P
    parts = [q for q in SPLIT_PARTS if q in ref.parts]
    views = [P.population_views(w, ref.odors, ref.breadth_edges, list(ref.parts)) for w in pops]
    names = [w.name or str(i) for i, w in enumerate(pops)]
    arrays = P.stack_views(views + [P.mean_views(views)], ref)
    rows = {name: {} for name in names + ['common']}
    for q in parts:
        m_all = ref.mean[q].to_numpy(dtype=float)
        unit = float(ref.one_fly[q])
        for name, a_all in zip(names + ['common'], arrays[q]):
            ok = np.isfinite(a_all) & np.isfinite(m_all)
            a, m = a_all[ok], m_all[ok]
            key = q.replace('/', '_')
            if len(a) < 2:
                rows[name].update({f'{key}_{k}': np.nan for k in SPLIT_TERMS})
                continue
            sa, sm = float(a.std()), float(m.std())
            r = float(np.corrcoef(a, m)[0, 1]) if sa > 0 and sm > 0 else np.nan
            cross = 2.0 * sa * sm - 2.0 * float(np.mean((a - a.mean()) * (m - m.mean())))
            rows[name].update({f'{key}_level': abs(float(np.mean(a - m))) / unit, f'{key}_spread': abs(sa - sm) / unit,
                               f'{key}_order': float(np.sqrt(max(cross, 0.0))) / unit,
                               f'{key}_sd_ratio': sa / sm if sm > 0 else np.nan, f'{key}_r': r})
    df = pd.DataFrame.from_dict(rows, orient='index')
    mean = df.loc[names].mean().to_frame('mean').T
    return pd.concat([df.loc[names], mean, df.loc[['common']]])


class StageLoss:
    """The megamat stage error of a set of runs: per stage of `stages` (default readings, the flies of `flies` for
    reading), the population errors of the wirings' mean views (hal's `wirings_as_flies`, row 'common') and of each
    wiring against `flies`' references (PN dendrites on four views, the other stages on three, as E134's
    `megamat_stage_error`), their root mean square over the stages; for the control, the same under each permutation and
    their mean. Returns a dict: 'stage_error' (the recorded labels), '<stage>_error' and its views, each wiring's
    '<stage>_error_<wiring>' and 'stage_error_<wiring>', the error over the stages after the PN dendrites
    ('stage_error_downstream', and each wiring's), 'stage_error_perm', '<stage>_error_perm' and each wiring's
    'stage_error_perm_<wiring>' (the permutations' means, for the control), 'loss_stages' (the permutations' mean for the control, else 'stage_error') and the rates. With
    `check_pn`, the PN-dendrite error must be the same at every call. With `split` (the score's), each stage's strength
    split (`strength_split`) as '<stage>_split_<column>' (the wirings' mean views) and '<stage>_split_<column>_<wiring>'."""

    def __init__(self, stages: Sequence[str], flies: Flies, lock: Optional[Mapping] = None,
                 subset_glomeruli: Optional[Sequence[str]] = None, check_pn: bool = False, split: bool = False):
        self.stages = tuple(stages)
        bad = [s for s in self.stages if s not in LOSS_STAGES]
        if bad or not self.stages:
            raise ValueError(f'stages {self.stages}')
        self.flies = flies
        self.lock = lock if lock is not None else R.require_lock()
        self.subset_glomeruli = subset_glomeruli
        self.check_pn = check_pn
        self.split = split
        self._pn: Optional[float] = None

    def read(self, runs: Mapping[str, 'R.SavedRun']) -> Dict[str, 'R.ModelReading']:
        """Each wiring's populations at the loss's stages (atlas's `read_run`, default readings, the wiring's seed
        index)."""
        return {w: R.read_run(run, stages=self.stages, lock=self.lock, flies=self.flies.cache, seed_index=seed_index(w),
                              subset_glomeruli=self.subset_glomeruli)
                for w, run in runs.items()}

    def __call__(self, runs: Mapping[str, 'R.SavedRun']) -> Dict[str, float]:
        from kcq import population_fit as P
        reads = self.read(runs)
        out: Dict[str, float] = {}
        common, each = [], {w: [] for w in runs}
        perm_common: List[List[float]] = [[] for _ in self.flies.perms]
        perm_each: Dict[str, List[List[float]]] = {w: [[] for _ in self.flies.perms] for w in runs}
        for st in self.stages:
            key = (st, self.flies.panel, 'on')
            pops = [reads[w].populations[key] for w in runs]
            odors = list(pops[0].values.columns)
            col = error_column(st, self.flies.views)
            tab = P.wirings_as_flies(pops, self.flies.reference(st, odors))
            out[f'{st}_error'] = float(tab.loc['common', col])
            out.update({f'{st}_{c}': float(tab.loc['common', c]) for c in tab.columns
                        if c.startswith(('view_', 'part_'))})
            common.append(out[f'{st}_error'])
            for w, pop in zip(runs, pops):
                out[f'{st}_error_{w}'] = float(tab.loc[pop.name, col])
                each[w].append(out[f'{st}_error_{w}'])
            if self.split:
                sp = strength_split(pops, self.flies.reference(st, odors))
                out.update({f'{st}_split_{c}': float(sp.loc['common', c]) for c in sp.columns})
                for w, pop in zip(runs, pops):
                    out.update({f'{st}_split_{c}_{w}': float(sp.loc[pop.name, c]) for c in sp.columns})
            for k in range(len(self.flies.perms)):
                tk = P.wirings_as_flies(pops, self.flies.reference(st, odors, k))
                perm_common[k].append(float(tk.loc['common', col]))
                for w, pop in zip(runs, pops):
                    perm_each[w][k].append(float(tk.loc[pop.name, col]))
            if self.flies.perms:
                out[f'{st}_error_perm'] = float(np.mean([pc[-1] for pc in perm_common]))

        def rms(xs):
            return float(np.sqrt(np.mean(np.square(xs))))
        out['stage_error'] = rms(common)
        down = [e for s, e in zip(self.stages, common) if s != 'pn_dendrite']
        out['stage_error_downstream'] = rms(down) if down else np.nan
        for w in runs:
            out[f'stage_error_{w}'] = rms(each[w])
        for w in runs:
            out[f'stage_error_downstream_{w}'] = rms([e for s, e in zip(self.stages, each[w]) if s != 'pn_dendrite']) \
                if any(s != 'pn_dendrite' for s in self.stages) else np.nan
        if self.flies.perms:
            out['stage_error_perm'] = float(np.mean([rms(pc) for pc in perm_common]))
            for w in runs:
                out[f'stage_error_perm_{w}'] = float(np.mean([rms(pe) for pe in perm_each[w]]))
        out['loss_stages'] = out['stage_error_perm'] if self.flies.perms else out['stage_error']
        out.update(rates(runs))
        if self.check_pn and 'pn_dendrite' in self.stages:
            v = out['pn_dendrite_error']
            if self._pn is None:
                self._pn = v
            elif v != self._pn:
                raise RuntimeError(f'the PN-dendrite error changed ({self._pn!r} -> {v!r}) though nothing searched acts '
                                   f'before the boutons')
        return out


# ---------------------------------------------------------------------------
# the parts: identity and the fit

def identity(arm: str = 'A0', wirings: Optional[Sequence[str]] = None) -> dict:
    """At z = 0 (each wiring at its held megamat point of E134's step 3): the hold changes nothing (one simulation per
    rate, no factor), every run reproduces E134's saved counts and responses exactly, the in-memory runs read exactly as
    the saved runs do (every population's values and responses), and, against all flies, each wiring's population
    errors reproduce E134's score_scores.csv ('megamat_held', the default readings) to IDENTITY_TOL where that file
    exists; written to identity_<arm>.json; RuntimeError otherwise."""
    from experiments import E134_population_baseline as B
    wirings = S_WIRINGS if wirings is None else wirings
    t0 = time.time()
    lock = R.require_lock()
    m = ArmModel(arm, wirings, lock=lock)
    runs = m.run({}, tag_prefix='identity_')
    flies = Flies(None, views=None)
    rec = dict(arm=arm, wirings=list(wirings), time=_now(), host=_host(), holds={}, counts={}, readings={}, scores={})
    ok = True
    for w, run in runs.items():
        h = run.info['hold']
        rec['holds'][w] = dict(x_thr=h['x_thr'], x_apl=h['x_apl'], sims=[h['off']['n_sims'], h['on']['n_sims']])
        ok &= h['x_thr'] == 0.0 and h['x_apl'] == 0.0 and h['off']['n_sims'] == 1 and h['on']['n_sims'] == 1
        saved = B.load_run(B.tag_of(B.Point(PANEL, 'held', w)))
        for c in R.CONDITIONS:
            for k in ('spike_counts', 'responses'):
                d = B._differ(R._abbrev(run.conditions[c][k]), R._abbrev(saved.conditions[c][k]))
                rec['counts'][f'{w}/{c}/{k}'] = d
                ok &= d == 0
        a = R.read_run(run, stages=LOSS_STAGES, lock=lock, flies=flies.cache, seed_index=seed_index(w))
        b = R.read_run(saved, stages=LOSS_STAGES, lock=lock, flies=flies.cache, seed_index=seed_index(w))
        for key, pa in a.populations.items():
            pb = b.populations[key]
            dv = B._differ(pa.values, pb.values)
            dr = B._differ(pa.responds.astype(float), pb.responds.astype(float))
            rec['readings']['/'.join((w,) + key)] = [dv, dr]
            ok &= dv == 0 and dr == 0
    res = StageLoss(LOSS_STAGES, flies, lock)(runs)
    rec['loss_all_flies'] = {k: v for k, v in res.items() if not k.endswith(tuple(f'_{w}' for w in wirings))}
    sc_path = e134_dir() / 'score_scores.csv'
    if sc_path.exists():
        sc = pd.read_csv(sc_path)
        for st in LOSS_STAGES:
            col = 'error_with_matrix' if st == 'pn_dendrite' else 'error'
            for w in wirings:
                x = sc[(sc.group == START_GROUP) & (sc.stage == st) & (sc.panel == PANEL) & (sc.condition == 'on')
                       & (sc.reading == R.DEFAULT_READINGS[st].name) & (sc.wiring == w)]
                if not len(x):
                    rec['scores'][f'{st}/{w}'] = 'not in score_scores.csv'
                    continue
                d = abs(float(x[col].iloc[0]) - res[f'{st}_error_{w}'])
                rec['scores'][f'{st}/{w}'] = d
                ok &= d <= IDENTITY_TOL
    else:
        rec['scores'] = "E134's score_scores.csv is not written yet"
    rec.update(exact=bool(ok), seconds=round(time.time() - t0, 1), **_commits())
    _write_json(out_dir() / f'identity_{arm}.json', rec)
    if not ok:
        raise RuntimeError(f'E136 identity fails: {rec}')
    return rec


def clip_coords(searched: Sequence[str], zv: Sequence[float]) -> Tuple[Tuple[float, ...], float]:
    """Coordinates clipped to BOUNDS and the bound penalty (one one-fly unit per BOUND_PENALTY_PER of overshoot,
    summed over the coordinates)."""
    zc, over = [], 0.0
    for k, v in zip(searched, zv):
        lo, hi = BOUNDS[k]
        c = float(min(max(float(v), lo), hi))
        over += abs(float(v) - c)
        zc.append(c)
    return tuple(zc), float(over / BOUND_PENALTY_PER)


def at_bound(searched: Sequence[str], zc: Sequence[float]) -> List[str]:
    """The coordinates at a bound as ruled: within AT_BOUND_REL of a factor's bound, AT_BOUND_HZ of the threshold's."""
    out = []
    for k, v in zip(searched, zc):
        lo, hi = BOUNDS[k]
        tol = AT_BOUND_HZ / CLAW_THR_UNIT if k == 'claw_thr' else float(np.log1p(AT_BOUND_REL))
        if v - lo <= tol:
            out.append(f'{k} low')
        if hi - v <= tol:
            out.append(f'{k} high')
    return out


def restore_memory(m: ArmModel, hist: pd.DataFrame) -> List[dict]:
    """The holds' warm starts and slopes of an arm model restored from a fit's history (its rows, in order), as the
    fit had them; returns the rows."""
    rows = hist.to_dict('records')
    for r in rows:
        c = tuple(float(r[f'z_{k}']) for k in m.searched)
        for w in m.wirings:
            m.memory[w].append((c, float(r[f'x_thr_{w}']), float(r[f'x_apl_{w}'])))
            m.slopes[w] = dict(thr=float(r[f'slope_thr_{w}']), apl=float(r[f'slope_apl_{w}']))
    return rows


def restart_point(name: str, hist: Sequence[dict], resumed: bool) -> dict:
    """The row the restart starts from: the best evaluation when the first run stopped, kept in fit_restart_<arm>.json
    so that a resumed fit starts its restart where it started it before. Without the record, the best of the history
    at a resume can be a point the restart itself found, and restarting from there begins the restart again with all
    of its budget left (A3's fit, 2026-11-25, began it from evals 22, 68 and 99 over its resumes). A fit that resumes
    with no record (one begun before the record existed) takes the best of the history, as those fits did."""
    path = out_dir() / f'fit_restart_{name}.json'
    rec = _read_json(path) if resumed else {}
    if rec:
        rows = [r for r in hist if int(r['eval']) == int(rec['from_eval'])]
        if not rows:
            raise RuntimeError(f"{path} names eval {rec['from_eval']}, which is not in the history")
        return rows[0]
    best = min(hist, key=lambda r: r['loss_sim'])
    _write_json(path, dict(arm=name, from_eval=int(best['eval']), loss=float(best['loss_sim']), n_evals=len(hist),
                           host=_host(), time=_now()))
    return best


def fit_arm(arm: str, control: bool = False, maxfev: Optional[int] = None, resume: bool = True,
            save_stages: bool = True, claw_mode: Optional[str] = None) -> dict:
    """One arm's fit on the S wirings against the fit half (`control`: the arm against the permuted fit half, as
    CONTROL for A0 and CONTROL_A1 for A1, E145's four-factor control), as the module docstring says. Every simulated evaluation's loss, parts, rates, holds and parameters are appended to
    fit_history_<arm>.csv and its stages saved (`save_eval`); resumable: an evaluation whose clipped coordinates are in
    the history is not simulated again, and the holds' warm starts and slopes are restored from it, so a restarted fit
    retraces its path and continues. A3's per-claw threshold is set in place only where claw_thr_check.json says it
    equals a fresh build (`claw_mode` overrides). Writes fit_<arm>.json."""
    from scipy.optimize import minimize
    from data_prep import stage_population_data as D
    gate = require_gate()
    stages = readable_stages(gate['classes'])
    lock = R.require_lock()
    name = (CONTROL_A1 if base_arm(arm) == 'A1' else CONTROL) if control else arm
    base = base_arm(name)
    if base == 'A1':
        require_e135()
    if 'claw_thr' in ARMS[base].searched and claw_mode is None:
        chk = _read_json(out_dir() / 'claw_thr_check.json')
        if not chk or set(chk.get('wirings', {})) < set(S_WIRINGS):
            chk = claw_thr_check()
        claw_mode = 'in_place' if chk['exact'] else 'rebuild'
    perms = permutations(list(D.panel_odors(PANEL))) if control else None
    flies = Flies('fit', stages, perms=perms)
    loss = StageLoss(stages, flies, lock, check_pn=base != 'A1')
    m = ArmModel(name, S_WIRINGS, lock=lock, claw_mode=claw_mode or 'in_place')
    searched = m.searched
    budget = int(maxfev or (MAXFEV_A1_SEARCH if m.pn_search else MAXFEV[name]))
    hist_path = out_dir() / f'fit_history_{name}.csv'
    memo: Dict[Tuple[float, ...], float] = {}
    hist: List[dict] = []
    if resume and hist_path.exists():
        for r in restore_memory(m, pd.read_csv(hist_path, float_precision='round_trip')):
            memo[tuple(float(r[f'z_{k}']) for k in searched)] = float(r['loss_sim'])
            hist.append(r)
        print(f'E136 {name}: resumed from {len(hist)} evaluations', flush=True)
    resumed = bool(hist)
    if not resumed:
        (out_dir() / f'fit_restart_{name}.json').unlink(missing_ok=True)
    calls = {'n': 0}
    t0 = time.time()

    def f(zv: np.ndarray) -> float:
        calls['n'] += 1
        zc, bpen = clip_coords(searched, zv)
        if zc not in memo:
            i = len(hist)
            z = dict(zip(searched, zc))
            t1 = time.time()
            runs = m.run(z, tag_prefix=f'{name}_e{i}_')
            res = loss(runs)
            res['hold_penalty'] = hold_penalty(runs)
            res['loss_sim'] = res['loss_stages'] + res['hold_penalty']
            if save_stages:
                save_eval(name, i, runs)
            row = dict(eval=i, **{f'z_{k}': v for k, v in z.items()}, **res, seconds=round(time.time() - t1, 1),
                       host=_host(), time=_now())
            for w, run in runs.items():
                hd = run.info['hold']
                row.update({f'x_thr_{w}': hd['x_thr'], f'x_apl_{w}': hd['x_apl'],
                            f'slope_thr_{w}': hd['off']['slope'], f'slope_apl_{w}': hd['on']['slope'],
                            f'hold_sims_{w}': hd['off']['n_sims'] + hd['on']['n_sims'],
                            f'hold_ok_{w}': bool(hd['off']['ok'] and hd['on']['ok']),
                            f'hold_at_end_{w}': bool(hd['off']['at_end'] or hd['on']['at_end'])})
                row.update({f'{k}_{w}': v for k, v in run.info['params'].items()})
                if run.info.get('claw_thr') is not None:
                    row[f'claw_thr_{w}'] = run.info['claw_thr']
                row.update({f'{k}_{w}': v for k, v in (run.info.get('pn_kws') or {}).items()})
            if i == 0 and (name == 'A0' or m.pn_search) and not any(zc) and MEASURE in START_CHECK:
                bad = {k: (res[k], v) for k, v in START_CHECK[MEASURE].items() if abs(res[k] - v) > START_CHECK_TOL}
                if bad:
                    raise RuntimeError(f'the start check fails (got, expected): {bad}')
            hist.append(row)
            memo[zc] = res['loss_sim']
            pd.DataFrame(hist).to_csv(hist_path, index=False)
            print(f"E136 {name} eval {i}: loss {res['loss_sim']:.4f} (stages {res['loss_stages']:.4f}, holds "
                  f"{res['hold_penalty']:.3f}; rates {res['rate_on']:.4f} / {res['rate_off']:.4f}), "
                  f"{row['seconds']:.0f} s", flush=True)
        return memo[zc] + bpen

    def run_nm(x0: np.ndarray, steps: np.ndarray, n: int):
        simplex = np.vstack([x0] + [x0 + steps[i] * np.eye(len(x0))[i] for i in range(len(x0))])
        return minimize(f, x0, method='Nelder-Mead',
                        options=dict(maxfev=n, initial_simplex=simplex, xatol=XATOL, fatol=FATOL))

    steps = np.array([FIRST_STEP[k] for k in searched], dtype=float)
    res = run_nm(np.zeros(len(searched)), steps, budget)
    messages = [str(res.message)]
    left = budget - calls['n']
    restarted, restart_from = False, None
    if left >= RESTART_MIN_LEFT:
        start = restart_point(name, hist, resumed)
        restart_from = int(start['eval'])
        x1 = np.array([float(start[f'z_{k}']) for k in searched])
        res2 = run_nm(x1, steps / 2.0, left)
        messages.append(str(res2.message))
        restarted = True
    best = min(hist, key=lambda r: r['loss_sim'])
    bz = {k: float(best[f'z_{k}']) for k in searched}
    flags = at_bound(searched, [bz[k] for k in searched])
    flags += [f'hold at an end on {w}' for w in S_WIRINGS if best.get(f'hold_at_end_{w}')]
    out = dict(arm=name, base_arm=base, searched=list(searched), stages=list(stages), start=START_GROUP, measure=MEASURE,
               gate_measure=gate['gate'].get('measure'), gate_verdict_as_read=gate['gate'].get('verdict_as_read'),
               n_evals=len(hist), n_calls=calls['n'], budget=budget, restarted=restarted, restart_from=restart_from,
               best_eval=int(best['eval']),
               best_z=bz, at_bound=flags, best=dict(best), messages=messages, claw_mode=m.claw_mode,
               seconds=round(time.time() - t0, 1), host=_host(), time=_now(), **_commits())
    out['factors'] = {k: coord_value(k, v) for k, v in bz.items()}
    if base == 'A1':
        out['e135'] = require_e135()
        out['pn_kws'] = pn_norm_kws(bz) if m.pn_search else None
        out['pn_checks'] = m.pn_checks or None
    if save_stages:
        out['archived'] = [str(p) for p in archive_best(name, int(best['eval']))]
    _write_json(out_dir() / f'fit_{name}.json', out)
    return out


def determined(arm: str, threshold: float = NOT_DETERMINED, save_stages: bool = True) -> dict:
    """The single-factor returns of the re-read on the adopted measure: from a fit's best point, each searched factor
    alone returned to the family's value (coordinate 0: factor 1, no per-claw threshold), both rates held again (warm
    starts from the fit's history), one evaluation each on S against the fit's half and measure. A factor whose return
    raises the loss by less than `threshold` is "not determined". If any is, the point with all of them returned (one
    more evaluation where more than one) becomes the arm's classed point when its loss is within `threshold` of the
    best's; otherwise the fit stands as found. fit_<arm>.json gains 'determined' (per factor: the loss returned, its
    rise, determined or not), 'classed_z', 'classed_loss' and 'at_bound_classed' (`at_bound` on the determined factors
    of the classed point); every evaluation is a row of determined_<arm>.csv and its stages are saved
    ('<arm>_det_<label>_<wiring>')."""
    from data_prep import stage_population_data as D
    gate = require_gate()
    stages = readable_stages(gate['classes'])
    lock = R.require_lock()
    fit = _read_json(out_dir() / f'fit_{arm}.json')
    if not fit:
        raise RuntimeError(f'{arm} is not fitted yet (fit_{arm}.json)')
    control = arm in CONTROLS
    base = base_arm(arm)
    perms = permutations(list(D.panel_odors(PANEL))) if control else None
    loss = StageLoss(stages, Flies('fit', stages, perms=perms), lock, check_pn=base != 'A1')
    m = ArmModel(arm, S_WIRINGS, lock=lock, claw_mode=fit.get('claw_mode') or 'in_place')
    restore_memory(m, pd.read_csv(out_dir() / f'fit_history_{arm}.csv', float_precision='round_trip'))
    best_z = {k: float(v) for k, v in fit['best_z'].items()}
    best_loss = float(fit['best']['loss_sim'])
    rows, det = [], {}

    def evaluate(label: str, z: Dict[str, float]) -> float:
        t1 = time.time()
        runs = m.run(z, tag_prefix=f'{arm}_det_{label}_')
        res = loss(runs)
        res['hold_penalty'] = hold_penalty(runs)
        res['loss_sim'] = res['loss_stages'] + res['hold_penalty']
        if save_stages:
            save_eval(arm, f'det_{label}', runs)
        rows.append(dict(label=label, **{f'z_{k}': v for k, v in z.items()}, **res, seconds=round(time.time() - t1, 1),
                         host=_host(), time=_now()))
        pd.DataFrame(rows).to_csv(out_dir() / f'determined_{arm}.csv', index=False)
        print(f"E136 {arm} {label}: loss {res['loss_sim']:.4f} (best {best_loss:.4f})", flush=True)
        return res['loss_sim']

    for k in m.searched:
        if best_z[k] == 0.0:
            det[k] = dict(returned_loss=best_loss, rise=0.0, determined=False, note='at the family value already')
            continue
        lk = evaluate(f'{k}0', dict(best_z, **{k: 0.0}))
        det[k] = dict(returned_loss=lk, rise=lk - best_loss, determined=bool(lk - best_loss >= threshold))
    free = [k for k in m.searched if not det[k]['determined'] and best_z[k] != 0.0]
    classed, classed_loss = dict(best_z), best_loss
    if free:
        z_all = dict(best_z, **{k: 0.0 for k in free})
        l_all = det[free[0]]['returned_loss'] if len(free) == 1 else evaluate('none', z_all)
        if l_all - best_loss <= threshold:
            classed, classed_loss = z_all, l_all
    dets = [k for k in m.searched if det[k]['determined']]
    fit.update(determined=det, classed_z=classed, classed_loss=classed_loss, not_determined_threshold=threshold,
               at_bound_classed=at_bound(dets, [classed[k] for k in dets]), determined_time=_now(),
               determined_commits=_commits())
    _write_json(out_dir() / f'fit_{arm}.json', fit)
    return fit


def history_evals(arm: str) -> List[int]:
    """The evaluations of a fit's history, in order."""
    h = pd.read_csv(out_dir() / f'fit_history_{arm}.csv')
    return [int(i) for i in h['eval']]


def rescore(arm: str, measures: Sequence[str] = tuple(MEASURES), halves: Sequence[Optional[str]] = ('fit', 'held'),
            evals: Optional[Sequence[int]] = None) -> pd.DataFrame:
    """A fit's saved evaluations (this host's outputs/.../fit_<arm>/, the history's evaluations) scored again under
    every measure of `measures` and against every half of `halves`, with no simulation (the re-read on the adopted
    measure: every view and part, the population error over all views, the two-view error). Writes
    rescore_<arm>.csv (one row per evaluation, measure and half) to the fit's out_dir, merged with the rows already
    there: a row this call scores replaces the one with the same evaluation, measure and half, and the others are kept
    (so a fit whose evaluations are saved on two hosts is rescored on each). An evaluation not saved on this host is
    skipped and listed. Returns the merged table."""
    lock = R.require_lock()
    evals = history_evals(arm) if evals is None else list(evals)
    runs_of: Dict[int, Dict[str, 'R.SavedRun']] = {}
    missing = []
    for i in evals:
        try:
            runs_of[i] = {w: R.load_step3_run(f'{arm}_e{i}_{w}', dirs=fit_dirs(arm)) for w in S_WIRINGS}
        except FileNotFoundError:
            missing.append(i)
    rows = []
    for meas in measures:
        for half in halves:
            loss = StageLoss(LOSS_STAGES, Flies(half, LOSS_STAGES, views=MEASURES[meas]), lock)
            for i, runs in runs_of.items():
                res = loss(runs)
                rows.append(dict(arm=arm, eval=i, measure=meas, half=half or 'all',
                                 **{k: v for k, v in res.items() if not k.startswith('rate_')}))
    df = pd.DataFrame(rows)
    path = out_dir() / f'rescore_{arm}.csv'
    if path.exists() and len(df):
        old = pd.read_csv(path)
        key = ['eval', 'measure', 'half']
        done = pd.MultiIndex.from_frame(df[key])
        old = old[~pd.MultiIndex.from_frame(old[key]).isin(done)]
        df = pd.concat([old, df], ignore_index=True).sort_values(['measure', 'half', 'eval'], kind='stable',
                                                                 ignore_index=True)
    elif path.exists():
        df = pd.read_csv(path)
    tmp = path.with_suffix('.csv.tmp')
    df.to_csv(tmp, index=False)
    tmp.replace(path)
    if missing:
        print(f'E136 rescore {arm}: {len(missing)} evaluations not saved on this host: {missing}', flush=True)
    return df


def consistency(arms: Optional[Sequence[str]] = None, measure_name: Optional[str] = None) -> pd.DataFrame:
    """The measure's consistency at the megamat stages, a description (the re-read on the adopted measure): over all
    the rescored evaluations of all arms (rescore_<arm>.csv), the Spearman ρ between the megamat stage error against the
    fit half and against the held-out half, with 0.8 as the yardstick, for the stage error and each stage."""
    from scipy import stats
    meas = measure_name or MEASURE
    fs = [out_dir() / f'rescore_{a}.csv' for a in (arms or [*ARMS, CONTROL])]
    df = pd.concat([pd.read_csv(f) for f in fs if f.exists()], ignore_index=True) if any(f.exists() for f in fs) \
        else pd.DataFrame()
    rows = []
    if not len(df):
        return pd.DataFrame(rows)
    d = df[df.measure == meas]
    a = d[d.half == 'fit'].set_index(['arm', 'eval'])
    b = d[d.half == 'held'].set_index(['arm', 'eval'])
    idx = a.index.intersection(b.index)
    for col in ['stage_error', *(f'{s}_error' for s in LOSS_STAGES)]:
        if col in a and col in b and len(idx) > 2:
            rho = float(stats.spearmanr(a.loc[idx, col], b.loc[idx, col])[0])
            rows.append(dict(measure=meas, statistic=col, n_evals=len(idx), rho=rho, yardstick=0.8))
    out = pd.DataFrame(rows)
    out.to_csv(out_dir() / f'consistency_{meas}.csv', index=False)
    return out


# ---------------------------------------------------------------------------
# the reads

def read_tag(pt: ReadPoint) -> str:
    """'<arm>_<panel>_<wiring>' (the arm 'family' for the unrefitted family)."""
    if pt.panel not in READ_PANELS or pt.wiring not in WIRINGS or pt.arm not in (*ARMS, FAMILY, CONTROL) \
            or (pt.arm == CONTROL and pt.panel != PANEL):
        raise ValueError(f'{pt}')
    return f'{pt.arm}_{pt.panel}_{pt.wiring}'


def planned_reads(arms: Optional[Sequence[str]] = None) -> List[ReadPoint]:
    """The reads at fixed parameters, in the order a forward lane takes them: the unrefitted family's validation2 reads
    (its other reads are E134's and E90's saved runs, `load_read`), then for each arm its fit on the megamat input on
    the T wirings (read against the held-out half; its S wirings are its fit's best evaluation), then on the eight
    wirings with the validation2, the mixture and the APL16 inputs; the control's best point on the T wirings (read
    against the held-out half under the permutations); `arms` (FAMILY and CONTROL among them for theirs) restricts
    the list."""
    want = [FAMILY, *ARMS, CONTROL] if arms is None else list(arms)
    pts = [ReadPoint(FAMILY, 'validation2', w, None) for w in WIRINGS] if FAMILY in want else []
    for a in want:
        if a == FAMILY:
            continue
        if a == CONTROL:
            pts += [ReadPoint(CONTROL, PANEL, w, 'held') for w in T_WIRINGS]
            continue
        pts += [ReadPoint(a, PANEL, w, 'held') for w in T_WIRINGS]
        pts += [ReadPoint(a, p, w, None) for p in READ_PANELS if p != PANEL for w in WIRINGS]
    return pts


def fitted_z(arm: str) -> Dict[str, float]:
    """An arm's fitted coordinates: its classed point after the single-factor returns (fit_<arm>.json 'classed_z',
    `determined`), else its best point ('best_z'); {} for the unrefitted family. RuntimeError if the fit is not
    written."""
    if arm == FAMILY:
        return {}
    f = _read_json(out_dir() / f'fit_{arm}.json')
    if not f:
        raise RuntimeError(f'{arm} is not fitted yet (fit_{arm}.json)')
    return {k: float(v) for k, v in (f.get('classed_z') or f['best_z']).items()}


def read_targets(wiring: str, panel: str, n_odors: int, lock: Optional[Mapping] = None
                 ) -> Tuple[int, int, float, float]:
    """(target_off, target_on, fraction_off, fraction_on) of a read as ruled: megamat E134's `held_targets`;
    validation2 the lock's rate of its flies and the mixtures E134's held mixture rate, each with the wiring's own
    APL-off / APL-on ratio; the APL16 input the wiring's own fit's two counts."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E134_population_baseline as B
    F = E83.saved_fit(wiring)
    size = F['saved']['on'][0].shape[0] * int(n_odors)
    lk = lock if lock is not None else R.require_lock()
    if panel == PANEL:
        return B.held_targets(B.Point(PANEL, 'held', wiring), size, lk)
    if panel == 'apl16':
        return int(F['target_off']), int(F['target_on']), F['target_off'] / F['size'], F['target_on'] / F['size']
    if panel == 'validation2':
        f_on = float(lk['populations']['kc_soma__validation2__wild_type']['response_rate']['mean_over_flies'])
    elif panel == 'mixtures':
        f_on = float(B.held_rate('mixtures', lk))
    else:
        raise ValueError(f'{panel=}')
    f_off = f_on * F['target_off'] / F['target_on']
    return int(round(f_off * size)), int(round(f_on * size)), float(f_off), float(f_on)


def read_dirs() -> List[Path]:
    """Where the reads' runs are: outputs/.../read/, then the archive's."""
    arch = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [local_dir('read'), (arch / _measure_sub() if _measure_sub() else arch) / 'read']


def _append_row(path: Path, row: Mapping) -> None:
    """Appends one row to a lane's table, rewriting it with the union of the columns (atomically) when they differ;
    only the lane's own process writes its file."""
    new = pd.DataFrame([dict(row)])
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        new.to_csv(path, index=False)
        return
    old = pd.read_csv(path)
    if list(old.columns) == list(new.columns):
        new.to_csv(path, mode='a', header=False, index=False)
        return
    cols = list(old.columns) + [c for c in new.columns if c not in old.columns]
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    pd.concat([old, new], ignore_index=True, sort=False)[cols].to_csv(tmp, index=False)
    os.replace(tmp, path)


def reads_table() -> pd.DataFrame:
    """Every lane's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def hold_at(arm: str, z: Mapping[str, float], wiring: str, panel: str, label: str = 'E136',
            lock: Optional[Mapping] = None, extra_kws: Optional[Mapping] = None) -> dict:
    """Both response rates of `arm`'s model (built as `base_arm(arm)`) on `wiring` and `panel` at coordinates `z`
    (applied as log factors to the wiring's held megamat point, `start_params`) held exactly at the panel's operating
    point, as a read holds them (`read_at`): E77's `hold_rates` on E83's simulator with the lab input's breadth, the
    threshold with APL off, then the tied APL scales, to `read_targets`' counts. Returns dict(hold: hold_rates' dict,
    kws: the model kwargs beside the family's (the arm's, A3's per-claw threshold, A1's PN stage), orn: the panel's
    input, breadth: the breadth behind the APL->KC weights, claw_thr: A3's threshold in Hz or None, targets:
    (target_off, target_on, fraction_off, fraction_on), p0: the parameters the hold started from, off_extra: the APL-off
    condition's own parameters or None). An arm whose kwargs carry 'kckc_off_factor' (OFF_ONLY_KWS; E141's inhibition
    among KCs present with APL on only) runs its APL-off condition with the KC-KC scale at that factor of its APL-on
    scale (E77's `hold_rates` `off_extra`), the threshold held there. `extra_kws` (default none) are further model
    kwargs beside the arm's, added last (e.g. E142's route-1 pn_rates on A1's model; they are in the returned 'kws')."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E134_population_baseline as B
    from experiments.E77_input_screen import hold_rates
    base = base_arm(arm)
    a = ARMS[base]
    lock = lock if lock is not None else R.require_lock()
    p0, _ = start_params(wiring)
    if 'kckc_fixed' in a.kws:
        p0['wKCKC_scale'] = float(a.kws['kckc_fixed'])
    for k, v in z.items():
        if k in FACTORS:
            for kk in FACTORS[k]:
                p0[kk] = p0[kk] * float(np.exp(v))
    ct = CLAW_THR_UNIT * float(z['claw_thr']) if 'claw_thr' in z else None
    orn = B.panel_input(panel)
    kws = dict(arm_kws(base), **({'claw_thr': ct} if ct else {}))
    if base == 'A1':
        kws.update(a1_kws(wiring, panel, orn, z))
    if extra_kws:
        kws.update(extra_kws)
    Pp = E83.wiring_parts(wiring)
    targets = read_targets(wiring, panel, orn.shape[1], lock)
    t_off, t_on = targets[0], targets[1]
    rng = HOLD_RANGE_DIVISIVE_APL if kws.get('apl_inh_divisive') else HOLD_RANGE
    off_extra = off_extra_of(base, p0)
    hold = hold_rates(lambda: E83.build_simulator(wiring, orn=orn, breadth=Pp['breadth'], params=dict(p0),
                                                  extra_kws=kws or None),
                      target_off=t_off, target_on=t_on, params=dict(p0), label=label, apl_bracket=rng,
                      off_extra=off_extra)
    return dict(hold=hold, kws=kws, orn=orn, breadth=Pp['breadth'], claw_thr=ct, targets=tuple(targets), p0=p0,
                off_extra=off_extra)


def read_at(tag: str, arm: str, z, wiring: str, panel: str, half: Optional[str], lane: str, root: Path,
            dirs: Sequence[Path], reclaim: bool = False, done: Optional[set] = None, extra: Optional[Mapping] = None,
            label: str = 'E136', script: str = SCRIPT) -> Optional[dict]:
    """One read of an arm's model at given coordinates, claimed as E134's step-3 points are (mkdir of
    <root>/claims/<tag>/; None if done or claimed elsewhere): `arm`'s model (built as `base_arm(arm)`) on `wiring` and
    `panel` at coordinates `z` (a mapping of coordinates, or a function returning it, called once the claim is held)
    on the wiring's held megamat point, both rates held exactly at the panel's operating point (E77's `hold_rates`, the
    threshold with APL off, then the tied APL scales; `read_targets`), built again at the held parameters and run with
    every stage saved (which must give the hold's outputs exactly); A3's claws also read after the threshold (saved as
    '<tag>_clawafter'). Saved in each of `dirs` (E134's step-3 layout); its row (with `extra`'s fields after the half,
    also recorded in the run's info) in <root>/reads_<host>_<lane>.csv and its provenance in
    <root>/provenance_read_<host>_<lane>.json (`script` named there). The claim is released on failure. `label` heads
    the printed lines."""
    from experiments import E134_population_baseline as B
    from experiments.E77_input_screen import frame_sha256
    root = Path(root)
    extra = dict(extra or {})
    if done is not None and tag in done:
        return None
    if not B.claim(tag, lane, reclaim, root=root):
        print(f'{label} {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        z = dict(z() if callable(z) else z)
        base = base_arm(arm)
        lock = R.require_lock()
        h = hold_at(arm, z, wiring, panel, label=f'{label} {tag}', lock=lock)
        hold, kws, orn, ct = h['hold'], h['kws'], h['orn'], h['claw_thr']
        t_off, t_on, f_off, f_on = h['targets']
        built = B.build_family(wiring, panel, hold['params'], orn=orn, extra_kws=kws or None)
        p, poff = dict(hold['params']), dict(hold['params_off'])
        reader = claw_reader_of(base)
        out = B.run_stages(B.Point(PANEL, 'held', wiring), built, p, poff, claw_reader=reader)
        for c in R.CONDITIONS:
            if B._differ(out[c]['spike_counts'], hold[c][1]) or B._differ(out[c]['responses'], hold[c][0]):
                raise RuntimeError(f"{tag} {c}: the saved run differs from the hold's outputs")
        info = dict(model=B.FAMILY, wiring=wiring, panel=panel, reading='read', arm=arm, z=z, params=p,
                    params_off=poff, claw_thr=ct, tag=tag, host=_host(), lane=lane, time=_now(), **extra)
        pn_rec = {}
        if base == 'A1':
            h = require_e135()
            pn_rec = dict(pn_mode=h['mode'], pn_rung=h['rung'], e135_outcome=h['outcome'],
                          **{k: float(v) for k, v in kws.items() if k in ('pn_inhadd', 'pn_inhsc')})
            chk = a1_pn_check(built, wiring, panel)
            info.update(pn_rec, pn_check=chk)
            if chk is not None:
                pn_rec['pn_check_max_abs_diff'] = max(chk['max_abs_diff_set'], chk['max_abs_diff_unset'])
        R.save_step3_run(tag, out, built.meta, built.orn, info, dirs)
        if ct:
            out2 = B.run_stages(B.Point(PANEL, 'held', wiring), built, p, poff,
                                claw_reader=claw_reader_after_threshold(ct, bool(kws.get('apl_inh_divisive'))))
            R.save_step3_run(f'{tag}_clawafter', out2, built.meta, built.orn,
                             dict(info, tag=f'{tag}_clawafter', reading='read, claws after the threshold'), dirs)
        rr = {}
        for c in R.CONDITIONS:
            rr.update({f'rate_{c}': float((out[c]['responses'].to_numpy(dtype=float) > 0.5).mean()),
                       f'max_count_{c}': float(out[c]['spike_counts'].to_numpy(dtype=float).max())})
        row = dict(tag=tag, arm=arm, panel=panel, wiring=wiring, half=half, **extra,
                   **{f'z_{k}': v for k, v in z.items()}, claw_thr=ct, **pn_rec,
                   **{k: float(v) for k, v in p.items()}, **{f'off_{k}': float(v) for k, v in poff.items()}, **rr,
                   target_on=t_on, target_off=t_off, frac_on_target=f_on, frac_off_target=f_off,
                   thr_factor=hold['thr_factor'], apl_factor=hold['apl_factor'], at_bound_thr=hold['at_bound_thr'],
                   at_bound_apl=hold['at_bound_apl'], hold_n_sims=hold['n_sims'], seconds=round(time.time() - t0, 1),
                   host=_host(), lane=lane, time=_now(), threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=frame_sha256(built.orn), **_commits())
        _append_row(root / f'reads_{_host()}_{lane}.csv', row)
        prov_path = root / f'provenance_read_{_host()}_{lane}.json'
        prov = _read_json(prov_path) or dict(part='read', lane=lane, script=script, points={})
        prov['points'][tag] = dict(info={k: v for k, v in built.info.items() if not str(k).startswith('_')},
                                   params=p, params_off=poff, model=H.model_provenance(orn_deltas=built.orn))
        _write_json(prov_path, prov)
        print(f"{label} {tag}: {row['seconds']:.0f} s, rate on {rr['rate_on']:.4f} / off {rr['rate_off']:.4f}",
              flush=True)
        return row
    except BaseException:
        B.release(tag, root=root)
        raise


def read_arm(pt: ReadPoint, lane: str, reclaim: bool = False, done: Optional[set] = None,
             archive: bool = True) -> Optional[dict]:
    """One read (`read_at`), claimed under this experiment's out_dir() as E134's step-3 points are (None if done or
    claimed elsewhere): the arm on the point's wiring and panel at its fitted coordinates (`fitted_z`) on the wiring's
    held megamat point, both rates held exactly at the panel's operating point; saved locally and in the archive
    (`read_dirs`; locally only without `archive`), its row in reads_<host>_<lane>.csv and its provenance in
    provenance_read_<host>_<lane>.json. The claim is released on failure."""
    dirs = read_dirs() if archive else read_dirs()[:1]
    return read_at(read_tag(pt), pt.arm, lambda: fitted_z(pt.arm), pt.wiring, pt.panel, pt.half, lane, out_dir(), dirs,
                   reclaim=reclaim, done=done)


def part_read(lane: str, arms: Optional[Sequence[str]] = None, reverse: bool = False, reclaim: bool = False) -> None:
    """A lane of reads: every planned read whose arm is fitted (the family's always), skipping what is done or claimed
    elsewhere."""
    pts = planned_reads(arms)
    pts = [pt for pt in (pts[::-1] if reverse else pts)
           if pt.arm == FAMILY or (out_dir() / f'fit_{pt.arm}.json').exists()]
    for pt in pts:
        read_arm(pt, lane, reclaim, set(reads_table().tag))


# ---------------------------------------------------------------------------
# the score

# variants scored as E136 scores its arms whose runs live elsewhere (e.g. E139's points): a name that is not an arm, the
# control or the family -> its loader, called as loader(pt, claw_after) and returning the saved run or None
# (`register_variant`); `load_read` consults it first, and `score_megamat` reads such a variant on every wiring from it
VARIANT_LOADERS: Dict[str, Callable[..., Optional['R.SavedRun']]] = {}


def register_variant(name: str, loader: Callable[..., Optional['R.SavedRun']]) -> None:
    """Registers a variant for E136's score functions (VARIANT_LOADERS): `loader(pt, claw_after)` gives the run of
    `ReadPoint(name, panel, wiring, half)` or None. ValueError for a name of E136's own (an arm, CONTROL, FAMILY)."""
    if name in (*ARMS, *CONTROLS, FAMILY):
        raise ValueError(f'{name!r} is one of E136\'s own names')
    VARIANT_LOADERS[name] = loader


def fitted_arm(arm: str) -> bool:
    """Whether `arm` is one of E136's fitted arms (an arm or the control), whose S wirings are its fit's best
    evaluation; the unrefitted family and registered variants are read on every wiring."""
    return arm in ARMS or arm == CONTROL


def load_read(pt: ReadPoint, claw_after: bool = False) -> Optional['R.SavedRun']:
    """A saved read (locally, else the archive; with `claw_after` A3's run with its claws read after the threshold,
    the same run elsewhere); the unrefitted family's megamat and mixture runs are E134's saved megamat_held and
    mixtures_held runs, and its APL16 runs E90's part B (the wirings' own fits); a registered variant's from its loader
    (VARIANT_LOADERS); None if absent."""
    from experiments import E134_population_baseline as B
    if pt.arm in VARIANT_LOADERS:
        return VARIANT_LOADERS[pt.arm](pt, claw_after)
    if claw_after and pt.arm in ARMS and 'claw_thr' in ARMS[pt.arm].searched:
        try:
            return R.load_step3_run(f'{read_tag(pt)}_clawafter', dirs=read_dirs())
        except FileNotFoundError:
            return None
    if pt.arm == FAMILY and pt.panel in (PANEL, 'mixtures'):
        return B.load_run(B.tag_of(B.Point(pt.panel, 'held', pt.wiring)))
    if pt.arm == FAMILY and pt.panel == 'apl16':
        try:
            return R.load_e90_run(pt.wiring, 'B')
        except FileNotFoundError:
            return None
    try:
        return R.load_step3_run(read_tag(pt), dirs=read_dirs())
    except FileNotFoundError:
        return None


def score_megamat(arm: str, stages: Sequence[str], half: Optional[str] = 'held',
                  wirings: Optional[Sequence[str]] = None, perms: bool = False,
                  claw_after: bool = False, views='measure') -> Optional[Dict[str, float]]:
    """An arm's (or the unrefitted family's, or a registered variant's) megamat stages on `wirings` against a half of
    the flies (all of them for None; under the control's permutations with `perms`): `StageLoss`'s dict. A fitted
    arm's S wirings are its fit's best evaluation (saved under outputs/.../fit_<arm>/; `fitted_arm`); every other
    name's runs are its reads on every wiring (`load_read`). None while a run is missing."""
    from data_prep import stage_population_data as D
    wirings = T_WIRINGS if wirings is None else wirings
    runs = megamat_runs(arm, wirings, half, claw_after)
    if runs is None:
        return None
    lock = R.require_lock()
    flies = Flies(half, stages, perms=permutations(list(D.panel_odors(PANEL))) if perms else None, views=views)
    return StageLoss(stages, flies, lock, split=True)(runs)


def megamat_runs(arm: str, wirings: Sequence[str], half: Optional[str] = 'held',
                 claw_after: bool = False) -> Optional[Dict[str, 'R.SavedRun']]:
    """The megamat runs `score_megamat` reads, per wiring: a fitted arm's best evaluation on the S wirings, every other
    run a read (`load_read`); None while one is missing."""
    runs = {}
    for w in wirings:
        if fitted_arm(arm) and w in S_WIRINGS:
            f = _read_json(out_dir() / f'fit_{arm}.json')
            try:
                r = R.load_step3_run(f"{arm}_e{f['best_eval']}_{w}", dirs=fit_dirs(arm)) if f else None
            except FileNotFoundError:
                r = None
        else:
            r = load_read(ReadPoint(arm, PANEL, w, half), claw_after)
        if r is None:
            return None
        runs[w] = r
    return runs


def stage_without_odor(arm: str, stage: str = 'kc_soma', odor: str = 'aa', half: Optional[str] = 'held',
                       wirings: Optional[Sequence[str]] = None) -> Optional[pd.DataFrame]:
    """A megamat stage's strength error with and without one odor's entry (E140's outcome (d): the KC somata without
    aa's), read as `score_megamat` reads it (the runs of `megamat_runs`, the half's references in all flies' units):
    E138's `odor_entries`, whose 'view_error_without' takes the odor's term out of the squared error and keeps the
    number of odors, nothing renormalized, which equals the odor's entry set to the flies' value. One row per wiring
    and 'common' (the wirings' mean views): arm, stage, odor, half, row, view_error, view_error_without, the odor's
    view_share. None while a run is missing."""
    from experiments import E138_strength_by_odor_entries as E138
    wirings = T_WIRINGS if wirings is None else wirings
    runs = megamat_runs(arm, wirings, half)
    if runs is None:
        return None
    flies = Flies(half, [stage])
    lock = R.require_lock()
    pops = [R.read_run(runs[w], stages=[stage], lock=lock, flies=flies.cache, seed_index=seed_index(w)
                       ).populations[(stage, PANEL, 'on')] for w in wirings]
    ref = flies.reference(stage, list(pops[0].values.columns))
    ent = E138.odor_entries(pops, ref, names=list(wirings))
    ent = ent[ent['odor'].astype(str) == odor]
    if ent.empty:
        raise KeyError(f'{odor!r} is not among the {stage} odors')
    out = ent[['row', 'view_error', 'view_error_without', 'view_share']].copy()
    for k, v in (('half', 'all' if half is None else half), ('odor', odor), ('stage', stage), ('arm', arm)):
        out.insert(0, k, v)
    return out.reset_index(drop=True)


def score_panels(arm: str, flies: Optional['R.FlyCache'] = None, views='measure') -> pd.DataFrame:
    """An arm's KC somata on the held-out panels (the validation2, mixture and APL16 inputs on the eight wirings, every
    panel and condition the runs are read on with flies to read them against, default readings), under the measure in
    force: per (input, panel, condition) each wiring's population error, their mean and the error of their mean views
    (hal's `wirings_as_flies` against all the panel's flies, E134's `fly_side_read`), with the views and the strength
    split ('split_...', `strength_split`)."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    flies = flies if flies is not None else R.FlyCache()
    lock = R.require_lock()
    views = MEASURES[MEASURE] if views == 'measure' else views
    rows = []
    for panel in READ_PANELS:
        if panel == PANEL:
            continue
        runs = [load_read(ReadPoint(arm, panel, w, None)) for w in WIRINGS]
        if any(r is None for r in runs):
            continue
        reads = {w: R.read_run(r, stages=['kc_soma'], lock=lock, flies=flies, seed_index=seed_index(w))
                 for w, r in zip(WIRINGS, runs)}
        for key in sorted({k for rd in reads.values() for k in rd.populations}):
            st, pan, cond = key
            fl = flies.get(key)
            if fl is None:
                continue
            pops = [reads[w].populations[key] for w in WIRINGS]
            fp, have = B.fly_side_read(fl)
            od = [o for o in pops[0].values.columns if o in set(have)]
            ref = P.fly_reference(fp, od, B.breadth_edges(len(od)), parts=measure_parts(views))
            tab = P.wirings_as_flies(pops, ref)
            sp = strength_split(pops, ref).add_prefix('split_')
            for name, row in tab.iterrows():
                wiring = str(name).partition('/')[2] or str(name)
                rows.append(dict(arm=arm, input=panel, stage=st, panel=pan, condition=cond, wiring=wiring,
                                 views='all' if views is None else '+'.join(views), n_odors=len(od), n_flies=len(fp),
                                 **{k: float(v) for k, v in row.items()},
                                 **({k: float(v) for k, v in sp.loc[name].items()} if name in sp.index else {})))
    return pd.DataFrame(rows)


def score_e137_beside(arm: str) -> pd.DataFrame:
    """An arm's mixture reads on the eight wirings scored beside the class under E137's readings of the mixture flies
    (E137_BESIDE; E137's committed lock, each reading's own units and flies, the flies' trial noise added to the model
    after the curve, as E137's re-score reads E134's runs: `read_model_panel`, seeds by the wirings' order, which is
    E137's), under the measure in force and under every view: per panel, reading, views and wiring (and the 'mean' and
    'common' rows) the KC somata's population error and its views. Empty if a mixture read is missing."""
    from data_prep import stage_population_data as D
    from experiments import E137_heldout_readings as E137
    from kcq import population_fit as P
    runs = [load_read(ReadPoint(arm, 'mixtures', w, None)) for w in WIRINGS]
    if any(r is None for r in runs):
        return pd.DataFrame()
    lock = E137.load_lock()
    rows = []
    for panel, readings in E137_BESIDE.items():
        od = list(D.panel_odors(panel))
        for reading in readings:
            units = E137.PRIMARY_UNITS[reading]
            pop, noise = E137.reading_population(panel, reading, units, lock)
            fps = P.fly_populations(pop)
            pops = E137.read_model_panel(runs, panel, pop, noise)
            for views in dict.fromkeys((MEASURES[MEASURE], None)):
                ref = P.fly_reference(fps, od, D.breadth_bins(len(od)), parts=measure_parts(views))
                for name, r in P.wirings_as_flies(pops, ref).iterrows():
                    rows.append(dict(arm=arm, input='mixtures', stage='kc_soma', panel=panel, condition='on',
                                     reading=reading, units=units, views='all' if views is None else '+'.join(views),
                                     wiring=str(name).partition('/')[2] or str(name), n_odors=len(od),
                                     n_flies=len(fps), **{k: float(v) for k, v in r.items()}))
    return pd.DataFrame(rows)


def e137_beside(rows: pd.DataFrame, arms: Sequence[str]) -> Dict[str, dict]:
    """`score_e137_beside`'s rows compared as the classes compare a held-out panel (`compare` on T against the unrefitted
    family, the spread over its eight wirings), per arm and '<panel>/<reading>/<views>'; gates nothing."""
    out: Dict[str, dict] = {}
    if not len(rows):
        return out
    fam = rows[rows.arm == FAMILY]
    for key, g in fam.groupby(['panel', 'reading', 'views']):
        k = '/'.join(key)
        base = {str(r.wiring): float(r.error) for r in g.itertuples() if r.wiring in WIRINGS}
        spread = _sd(base.values())
        for a in arms:
            x = rows[(rows.arm == a) & (rows.panel == key[0]) & (rows.reading == key[1]) & (rows.views == key[2])]
            out.setdefault(a, {})[k] = compare({str(r.wiring): float(r.error) for r in x.itertuples()}, base, spread)
    return out


def apl16_realism(on: pd.DataFrame, off: pd.DataFrame) -> dict:
    """The realism bands of one APL16 run (kcq.apl_effect's, as its by_variant scores the archive: the model's two-spike
    responses against the KC data in baseline-noise units): the rates, their ratio, the response rate ρ with APL on and
    off, the odors-per-KC distances and the flag."""
    from kcq import apl_effect as A
    row = {}
    for c, r in (('on', on), ('off', off)):
        d = A.odors_per_kc(r)
        row[f'rate_{c}'] = float((r.to_numpy(dtype=float) > 0.5).mean())
        row[f'opk_{c}_tv'] = float(0.5 * np.abs(d - A.opk_data(c, 'noise')).sum())
        per_odor = H.abbrev_columns(r.astype(float) > 0.5).mean().reindex(H.APL16_ODORS)
        row[f'rate_rho_{c}'] = H.spearman(per_odor, A.frac_data(c, 'noise'))
    row['rate_ratio'] = row['rate_off'] / row['rate_on'] if row['rate_on'] > 0 else np.nan
    row['realism_flag'] = A.realism_flag(row['rate_on'], row['rate_off'], row['rate_rho_on'], row['opk_on_tv'],
                                         row['opk_off_tv'], units='noise')
    return row


def apl16_rows(arm: str) -> pd.DataFrame:
    """An arm's APL16 reads, one row per wiring: E67's `full_row` (the APL-effect statistics, the two ρ's) of its APL-on
    and APL-off runs, the realism bands (`apl16_realism`) and the largest counts."""
    from experiments.E67_apl_effect_fits import full_row
    rows = []
    for w in WIRINGS:
        r = load_read(ReadPoint(arm, 'apl16', w, None))
        if r is None:
            continue
        o = {c: (r.conditions[c]['responses'], r.conditions[c]['spike_counts']) for c in R.CONDITIONS}
        row = dict(arm=arm, wiring=w, **full_row(o))
        row.update(apl16_realism(o['on'][0], o['off'][0]))
        for c in R.CONDITIONS:
            row[f'max_count_{c}'] = float(o[c][1].to_numpy(dtype=float).max())
        rows.append(row)
    return pd.DataFrame(rows)


# the held-out panels the class reads ("holds"); the kiwi mixtures left the class on Tom's word (2026-11-25 10:42, in
# hal's terminal, relayed by hal: "go with your read on the mixture panels") and are scored beside it (BESIDE_PANELS),
# with E137's two readings of the components as bounds; the control mixtures stay in, read as uploaded
HOLDS_PANELS: Tuple[Tuple[str, str, str], ...] = (('validation2', 'validation2', 'on'), ('mixtures', 'cmix', 'on'),
                                                   ('apl16', 'apl16', 'on'), ('apl16', 'apl16', 'off'))
BESIDE_PANELS: Tuple[Tuple[str, str, str], ...] = (('mixtures', 'kmix', 'on'),)
SCORED_PANELS: Tuple[Tuple[str, str, str], ...] = HOLDS_PANELS + BESIDE_PANELS
# E137's readings of the mixture flies scored beside (its lock, each reading's own units and flies' trial noise): the
# kiwi mixtures' two bounds and the control mixtures' second block
E137_BESIDE: Dict[str, Tuple[str, ...]] = {'kmix': ('components_first_block', 'components_second_block'),
                                           'cmix': ('components_second_block',)}


def _sd(x) -> float:
    x = np.asarray(list(x), dtype=float)
    return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan


def megamat_values(sc: Mapping[str, float], stages: Sequence[str], wirings: Sequence[str]) -> Dict[str, Dict[str, float]]:
    """From `StageLoss`'s dict, per wiring its stage error ('stage'), its error over the stages after the PN dendrites
    ('downstream') and each stage's population error, with the strength split ('<stage>_split_...') where scored."""
    return {w: dict(stage=sc[f'stage_error_{w}'], downstream=sc.get(f'stage_error_downstream_{w}', np.nan),
                    **{st: sc[f'{st}_error_{w}'] for st in stages},
                    **{k[:-len(w) - 1]: v for k, v in sc.items() if '_split_' in k and k.endswith(f'_{w}')})
            for w in wirings}


def panel_values(panels: pd.DataFrame, arm: str) -> Dict[Tuple[str, str, str], Dict[str, float]]:
    """From `score_panels`' rows, per (input, panel, condition) of SCORED_PANELS each wiring's KC population error."""
    out = {}
    for inp, pan, cond in SCORED_PANELS:
        x = panels[(panels.arm == arm) & (panels.input == inp) & (panels.panel == pan) & (panels.condition == cond)
                   & panels.wiring.isin(WIRINGS)] if len(panels) else panels
        out[(inp, pan, cond)] = {str(r.wiring): float(r.error) for r in x.itertuples()} if len(x) else {}
    return out


def compare(arm_v: Mapping[str, float], base_v: Mapping[str, float], spread: float,
            wirings: Optional[Sequence[str]] = None) -> dict:
    """The classes' comparison of one value as ruled: the per-wiring changes (arm minus the reference on the same
    wiring) over `wirings`, their mean, the standard error over wirings (the larger of half the spread and the SD of
    the changes over the wirings divided by sqrt(n)), and whether the change is an improvement (below -IMPROVE_SE
    standard errors) or a worsening beyond the spread."""
    wirings = T_WIRINGS if wirings is None else wirings
    ch = [arm_v[w] - base_v[w] for w in wirings if w in arm_v and w in base_v]
    if len(ch) < len(wirings):
        return dict(n=len(ch), complete=False)
    mean = float(np.mean(ch))
    se = float(max(spread / 2.0 if np.isfinite(spread) else 0.0, _sd(ch) / np.sqrt(len(ch))))
    return dict(n=len(ch), complete=True, changes=[float(c) for c in ch], change=mean, se=se, spread=spread,
                arm=float(np.mean([arm_v[w] for w in wirings])), ref=float(np.mean([base_v[w] for w in wirings])),
                improves=bool(mean < -IMPROVE_SE * se), worse_than_spread=bool(np.isfinite(spread) and mean > spread))


def classify(arm: str, stages: Sequence[str], mega: Mapping[str, Mapping[str, Dict[str, float]]],
             panels: Mapping[str, Mapping[Tuple[str, str, str], Dict[str, float]]], spreads: Mapping[str, float],
             realism: Mapping[str, Mapping[str, str]], fit: Mapping, control_fall: Optional[float],
             guard: Optional[Mapping] = None) -> dict:
    """One arm's class as the entry fixes it and the rulings read it, on the T wirings and the held-out half, against
    the unrefitted family there (`mega` / `panels`: per arm, FAMILY included, the values of `megamat_values` /
    `panel_values`; `spreads`: 'stage', each stage and each SCORED_PANELS key as 'input/panel/condition', the SD over the
    family's eight wirings; `realism`: per arm, per wiring its APL16 read's realism flag; `fit`: the arm's
    fit_<arm>.json; `control_fall`: the control fit's fall of the megamat stage error on T against the permuted
    held-out half; `guard`: rule 3's guard on a measure of fewer views, the same values under every view: 'mega',
    'panels' and 'spreads' as above). "Lower" in "one stage at another's cost" is read with the same margin as
    "improves". With the guard, a stage whose error over all its views rises by more than that error's spread makes an
    arm that lowers the strength error "strength at another view's cost" (it does not improve), and a held-out panel is
    worse when either error is. A fit's 'at_bound_classed' (the determined factors of its classed point, `determined`)
    stands in for 'at_bound' where it exists. BESIDE_PANELS are compared the same way, into 'panels_beside' (their guard
    into 'guard'), and gate nothing."""
    a, f = mega[arm], mega[FAMILY]
    m = compare({w: v['stage'] for w, v in a.items()}, {w: v['stage'] for w, v in f.items()}, spreads['stage'])
    st = {s: compare({w: v[s] for w, v in a.items()}, {w: v[s] for w, v in f.items()}, spreads[s]) for s in stages}
    pn, pb = {}, {}
    for keys, out in ((HOLDS_PANELS, pn), (BESIDE_PANELS, pb)):
        for key in keys:
            k = '/'.join(key)
            out[k] = compare(panels[arm].get(key, {}), panels[FAMILY].get(key, {}), spreads.get(k, np.nan))
    rec = dict(arm=arm, megamat=m, stages=st, panels=pn, panels_beside=pb)
    reasons = []
    if not m.get('complete') or not all(v.get('complete') for v in st.values()):
        reasons.append('reads missing')
    bounds = fit.get('at_bound_classed') if 'at_bound_classed' in fit else fit.get('at_bound')
    if bounds:
        reasons.append('at a bound: ' + ', '.join(bounds))
    left = [w for w in T_WIRINGS if realism.get(arm, {}).get(w) and not realism.get(FAMILY, {}).get(w)]
    if left:
        reasons.append('leaves a realism band on ' + ', '.join(left))
    fall = -m['change'] if m.get('complete') else np.nan
    rec.update(fall=fall, control_fall=control_fall)
    if control_fall is not None and np.isfinite(fall) and not fall > control_fall:
        reasons.append('no better than the permuted control')
    stage_worse = [s for s, v in st.items() if v.get('worse_than_spread')]
    panel_worse = [k for k, v in pn.items() if v.get('complete') and v.get('worse_than_spread')]
    panel_missing = [k for k, v in pn.items() if not v.get('complete')]
    guard_worse, guard_rec = [], {}
    if guard is not None:
        ga, gf = guard['mega'].get(arm, {}), guard['mega'].get(FAMILY, {})
        for s in stages:
            c = compare({w: v[s] for w, v in ga.items()}, {w: v[s] for w, v in gf.items()}, guard['spreads'].get(s, np.nan))
            guard_rec[s] = c
            if c.get('worse_than_spread'):
                guard_worse.append(s)
        for key in SCORED_PANELS:
            k = '/'.join(key)
            c = compare(guard['panels'].get(arm, {}).get(key, {}), guard['panels'].get(FAMILY, {}).get(key, {}),
                        guard['spreads'].get(k, np.nan))
            guard_rec[k] = c
            if key in HOLDS_PANELS and c.get('complete') and c.get('worse_than_spread') and k not in panel_worse:
                panel_worse.append(k)
    improves = bool(m.get('improves')) and not stage_worse and not guard_worse
    rec.update(improves=improves, stage_worse=stage_worse, guard_worse=guard_worse, guard=guard_rec,
               panel_worse=panel_worse, panel_missing=panel_missing)
    if reasons:
        cls = 'no class'
    elif improves and not panel_worse and not panel_missing:
        cls = 'holds'
    elif improves and panel_worse:
        cls = 'megamat only'
    elif improves:
        cls = 'improves (held-out panels not all read)'
    elif m.get('improves') and stage_worse:
        cls = "one stage at another's cost"
    elif m.get('improves') and guard_worse:
        cls = "strength at another view's cost"
    else:
        cls = 'does not improve'
    if arm != 'A0' and 'A0' in mega:
        ad = compare({w: v['stage'] for w, v in a.items()}, {w: v['stage'] for w, v in mega['A0'].items()},
                     spreads['stage'])
        adds = bool(ad.get('improves'))
        if arm == 'A1':
            dn = compare({w: v['downstream'] for w, v in a.items()}, {w: v['downstream'] for w, v in mega['A0'].items()},
                         spreads.get('downstream', np.nan))
            adds = adds and bool(dn.get('improves'))
            rec['adds_downstream'] = dn
        rec.update(adds_vs_A0=ad, adds=adds)
    rec.update(cls=cls, reasons=reasons)
    return rec


def realism_of(rows: pd.DataFrame) -> Dict[str, str]:
    """Per wiring, the APL16 read's realism flag ('' inside every band)."""
    return {str(r.wiring): (r.realism_flag if isinstance(r.realism_flag, str) else '') for r in rows.itertuples()} \
        if len(rows) else {}


def yardsticks(stages: Sequence[str] = LOSS_STAGES) -> Dict[str, dict]:
    """The yardsticks beside every error on the measure in force (the re-read on the adopted measure), against the
    held-out half in all flies' units: 'flat' (a profile with every odor at the panel's mean: each strength part 1 for
    every odor; strength views only) and 'half_to_half' (the fit half's mean views scored against the held-out half);
    E134's no-step yardstick is its flies_no_step.csv."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    odors = list(D.panel_odors(PANEL))
    fit, held = Flies('fit', stages), Flies('held', stages)
    out = {}
    for st in stages:
        rh, rf = held.reference(st, odors), fit.reference(st, odors)
        rec = dict(half_to_half=float(P.score_views(dict(rf.mean), rh).error))
        flat = {p: pd.Series(1.0, index=m.index) for p, m in rh.mean.items() if p.startswith('strength/')}
        if flat:
            rec['flat'] = float(P.score_views(flat, P.subset_reference(rh, rh.flies, keep_units=True)).views.get(
                'strength', np.nan))
        out[st] = rec
    return out


def part_score(arms: Optional[Sequence[str]] = None) -> dict:
    """The reads scored and classed: per arm (and the unrefitted family) the megamat stages on T against the held-out
    half (and, beside, against all flies and the fit half, and on S from the fit), A3's under both claw readings, the
    control's against the permuted held-out half, the held-out panels' KC population errors, the APL16 rows with the
    realism bands; the spreads over the family's eight wirings; the classes (`classify`); beside, the mixture panels
    under E137's readings (`score_e137_beside`). Writes score_megamat.csv, score_panels.csv, score_panels_e137.csv,
    score_apl16.csv, classes.json and provenance_score.json."""
    gate = require_gate()
    stages = readable_stages(gate['classes'])
    t0 = time.time()
    arms = [a for a in (arms or ARMS) if (out_dir() / f'fit_{a}.json').exists()]
    mrows, mega, mega_after, panels_rows = [], {}, {}, []
    for a in [FAMILY, *arms]:
        for half in ('held', None, 'fit'):
            ws = WIRINGS if (a == FAMILY and half == 'held') else T_WIRINGS
            sc = score_megamat(a, stages, half, ws)
            if sc is None:
                continue
            vals = megamat_values(sc, stages, ws)
            if half == 'held':
                mega[a] = vals
            for w, v in vals.items():
                mrows.append(dict(arm=a, half=half or 'all', wiring=w, **v))
            mrows.append(dict(arm=a, half=half or 'all', wiring='common', stage=sc['stage_error'],
                              downstream=sc['stage_error_downstream'], **{s: sc[f'{s}_error'] for s in stages},
                              **{k: v for k, v in sc.items() if '_split_' in k
                                 and not any(k.endswith(f'_{w}') for w in WIRINGS)}))
        if a in ARMS and 'claw_thr' in ARMS[a].searched:
            sc = score_megamat(a, stages, 'held', T_WIRINGS, claw_after=True)
            if sc is not None:
                mega_after[a] = megamat_values(sc, stages, T_WIRINGS)
        p = score_panels(a)
        if len(p):
            panels_rows.append(p)
    panels = pd.concat(panels_rows, ignore_index=True) if panels_rows else pd.DataFrame()
    pv = {a: panel_values(panels, a) for a in [FAMILY, *arms]}
    fam = mega.get(FAMILY, {})
    spreads = dict(stage=_sd(v['stage'] for v in fam.values()), downstream=_sd(v['downstream'] for v in fam.values()),
                   **{s: _sd(v[s] for v in fam.values()) for s in stages})
    for key in SCORED_PANELS:
        spreads['/'.join(key)] = _sd(pv[FAMILY].get(key, {}).values())
    apl_rows = [apl16_rows(a) for a in [FAMILY, *arms]]
    apl = pd.concat([r for r in apl_rows if len(r)], ignore_index=True) if any(len(r) for r in apl_rows) else \
        pd.DataFrame()
    realism = {a: realism_of(apl[apl.arm == a]) if len(apl) else {} for a in [FAMILY, *arms]}
    control_fall = None
    if (out_dir() / f'fit_{CONTROL}.json').exists():
        fp = score_megamat(FAMILY, stages, 'held', T_WIRINGS, perms=True)
        cp = score_megamat(CONTROL, stages, 'held', T_WIRINGS, perms=True)
        if fp is not None and cp is not None:
            control_fall = float(np.mean([fp[f'stage_error_perm_{w}'] - cp[f'stage_error_perm_{w}'] for w in T_WIRINGS]))
    guard, beside = None, {}
    if MEASURE != 'all':
        # rule 3's guard: the same reads under every view (four at the PN dendrites), with their own spreads
        g_mega, g_rows = {}, []
        for a in [FAMILY, *arms]:
            ws = WIRINGS if a == FAMILY else T_WIRINGS
            sc = score_megamat(a, stages, 'held', ws, views=None)
            if sc is not None:
                g_mega[a] = megamat_values(sc, stages, ws)
            p = score_panels(a, views=None)
            if len(p):
                g_rows.append(p)
        g_panels = pd.concat(g_rows, ignore_index=True) if g_rows else pd.DataFrame()
        g_pv = {a: panel_values(g_panels, a) for a in [FAMILY, *arms]}
        gf = g_mega.get(FAMILY, {})
        g_spreads = dict(stage=_sd(v['stage'] for v in gf.values()), **{s: _sd(v[s] for v in gf.values()) for s in stages})
        for key in SCORED_PANELS:
            g_spreads['/'.join(key)] = _sd(g_pv[FAMILY].get(key, {}).values())
        guard = dict(mega=g_mega, panels=g_pv, spreads=g_spreads)
        g_panels.to_csv(out_dir() / 'score_panels_all_views.csv', index=False)
        # beside: the two-view error (strength and sharing) and its change on T
        tv = {}
        for a in [FAMILY, *arms]:
            sc = score_megamat(a, stages, 'held', T_WIRINGS, views=MEASURES['strength_sharing'])
            if sc is not None:
                tv[a] = {w: sc[f'stage_error_{w}'] for w in T_WIRINGS}
        beside['two_view'] = {a: compare(tv[a], tv[FAMILY], np.nan) for a in arms if a in tv and FAMILY in tv}
        beside['yardsticks'] = yardsticks(stages)
        beside['spreads_all_views'] = g_spreads
    # beside: the mixture panels under E137's readings (the kiwi mixtures' two bounds, the control mixtures' second block)
    e_rows = [score_e137_beside(a) for a in [FAMILY, *arms]]
    e137 = pd.concat([r for r in e_rows if len(r)], ignore_index=True) if any(len(r) for r in e_rows) else pd.DataFrame()
    beside['e137_readings'] = e137_beside(e137, arms)
    beside['panels_beside'] = ['/'.join(k) for k in BESIDE_PANELS]
    classes = {}
    for a in arms:
        fit = _read_json(out_dir() / f'fit_{a}.json')
        if a not in mega or FAMILY not in mega:
            classes[a] = dict(cls='reads missing')
            continue
        classes[a] = classify(a, stages, mega, pv, spreads, realism, fit, control_fall, guard)
        if a in mega_after:
            alt = classify(a, stages, dict(mega, **{a: mega_after[a]}), pv, spreads, realism, fit, control_fall, guard)
            classes[a]['claws_after_threshold'] = dict(cls=alt['cls'], megamat=alt['megamat'])
            if alt['cls'] != classes[a]['cls']:
                classes[a]['reasons'].append('the class differs between the claws read before and after the threshold')
                classes[a]['cls'] = 'no class'
    pd.DataFrame(mrows).to_csv(out_dir() / 'score_megamat.csv', index=False)
    panels.to_csv(out_dir() / 'score_panels.csv', index=False)
    e137.to_csv(out_dir() / 'score_panels_e137.csv', index=False)
    apl.to_csv(out_dir() / 'score_apl16.csv', index=False)
    try:
        beside['consistency'] = consistency(arms + ([CONTROL] if (out_dir() / f'fit_{CONTROL}.json').exists() else []))\
            .to_dict('records')
    except Exception as e:   # a description only: a missing rescore does not stop the score
        beside['consistency'] = f'not read: {e}'
    _write_json(out_dir() / 'classes.json', dict(stages=list(stages), measure=MEASURE, spreads=spreads,
                                                 control_fall=control_fall, beside=beside, classes=classes))
    _write_json(out_dir() / 'provenance_score.json', dict(part='score', script=SCRIPT, lock=R.require_lock().get('sha256'),
                seconds=round(time.time() - t0, 1), host=_host(), time=_now(), arms=arms, **_commits()))
    return classes


# ---------------------------------------------------------------------------
# CLI

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('identity', 'claw_identity', 'claw_thr_check', 'fit', 'determined', 'rescore', 'read',
                             'score', 'archive', 'variants'))
    ap.add_argument('--arm', choices=tuple(ARMS) + (CONTROL, FAMILY))
    ap.add_argument('--evals', nargs='+', type=int, help="--part archive: the evaluations (default this host's best); --part rescore: the evaluations "
                         "(default the history's)")
    ap.add_argument('--maxfev', type=int)
    ap.add_argument('--wirings', nargs='+', choices=WIRINGS)
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--measure', choices=tuple(MEASURES), default='all',
                    help="the views the loss and the classes read (outputs of a measure other than 'all' in "
                         "measure_<name>/)")
    a = ap.parse_args(argv)
    check_wirings()
    set_measure(a.measure)
    if a.part == 'identity' and a.measure != 'all':
        raise SystemExit('the identity reads every view (against E134\'s score)')
    if a.part == 'identity':
        print(json.dumps(identity(a.arm or 'A0', a.wirings or S_WIRINGS), indent=1, default=str))
    elif a.part == 'claw_identity':
        r = claw_identity(a.arm or 'A2')
        _write_json(out_dir() / f'claw_identity_{r["arm"]}.json', dict(r, **_commits()))
        print(json.dumps(r, indent=1))
    elif a.part == 'claw_thr_check':
        print(json.dumps(claw_thr_check(a.wirings or S_WIRINGS), indent=1, default=str))
    elif a.part == 'fit':
        if a.arm in (None, FAMILY):
            raise SystemExit('--part fit needs --arm (an arm, or the control)')
        fit_arm(base_arm(a.arm), control=a.arm == CONTROL, maxfev=a.maxfev)
    elif a.part == 'read':
        require_gate()
        part_read(a.lane, None if a.arm is None else [a.arm], a.reverse, a.reclaim)
    elif a.part == 'determined':
        if a.arm in (None, FAMILY):
            raise SystemExit('--part determined needs --arm')
        print(json.dumps(determined(a.arm), indent=1, default=str))
    elif a.part == 'rescore':
        if a.arm in (None, FAMILY):
            raise SystemExit('--part rescore needs --arm')
        print(rescore(a.arm, evals=a.evals).groupby(['measure', 'half']).stage_error.describe().to_string())
    elif a.part == 'score':
        print(json.dumps(part_score(None if a.arm is None else [a.arm]), indent=1, default=str))
    elif a.part == 'archive':
        if a.arm in (None, FAMILY):
            raise SystemExit('--part archive needs --arm')
        archive_evals(a.arm, a.evals)
    elif a.part == 'variants':
        write_variants()


if __name__ == '__main__':
    main()
