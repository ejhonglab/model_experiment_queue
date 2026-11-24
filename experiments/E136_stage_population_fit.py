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
if it stops with RESTART_MIN_LEFT or more evaluations left; a point outside BOUNDS is evaluated at the clipped point plus
one one-fly unit per BOUND_PENALTY_PER of overshoot. In A0, A2, A3 and A4 nothing searched acts before the boutons, so
the PN-dendrite error must be the same at every evaluation (checked).

Parts (CLI `--part`):
- identity: A0 at z = 0 against all flies: the hold changes nothing, every run reproduces E134's saved megamat_held
  counts and readings exactly, and each wiring's population errors reproduce E134's score_scores.csv where it exists.
- claw_identity: A2's divided claw drive against olfsysm's claw_sims; claw_thr_check: A3's in-place threshold against a
  fresh build at CLAW_THR_CHECK_HZ.
- fit: one arm's fit (`--arm`, CONTROL for the control), resumable from its history.
- read: the reads at the fitted parameters, held exactly at each panel's operating point (`planned_reads`, lanes and
  claims as E134's step 3): megamat on T against the held-out half; validation2 (the lock's 0.0504) and the mixtures
  (E134's held mixture rate) with the wiring's own ratio for APL off; the APL16 input at the wiring's own fit's two
  counts. The unrefitted family's reads are E134's megamat_held and mixtures_held runs and E90's part B (the APL16
  fit) where they exist; only its validation2 reads are run here.
- score: the reads scored and classed (`score_megamat`, `score_panels`, `apl16_rows`, `classify`, `part_score`).

The classes (the entry's, read as ruled): each wiring is one fly; a value is the mean over the four T wirings of each
wiring's error (the error of the T wirings' mean views beside); a change is the mean over T of the arm minus the
unrefitted family on the same wiring; a spread is the SD over the family's eight wirings of the unrefitted family's
value; the standard error over wirings is the larger of half the megamat stage error's spread and the SD of the
per-wiring changes over T divided by 2. "improves": the change is below -IMPROVE_SE standard errors and no stage's
change is above that stage's spread; "holds": it improves and no held-out panel (validation2, kmix, cmix, APL16 with
APL on, APL16 with APL off) is worse by more than its spread; "megamat only"; "one stage at another's cost"; "adds" on
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
FAMILY = 'family'   # the unrefitted family in the reads: z = 0, each wiring at its held megamat point

# the search (the rulings): the coordinates (ln factors; the per-claw threshold in units of CLAW_THR_UNIT Hz), their
# bounds, the first steps, the tolerances, the budgets and the restart
CLAW_THR_UNIT = 50.0
BOUNDS: Dict[str, Tuple[float, float]] = {
    'aplpn': (float(np.log(1.0 / 30.0)), float(np.log(1.07))),
    'kckc': (float(np.log(0.1)), float(np.log(7.3))),
    'claw_thr': (0.0, 100.0 / CLAW_THR_UNIT),
}
FIRST_STEP: Dict[str, float] = {'aplpn': -0.7, 'kckc': -0.7, 'claw_thr': 0.6}
BOUND_PENALTY_PER = 0.1
XATOL, FATOL = 0.02, 0.005
MAXFEV: Dict[str, int] = {'A0': 80, 'A1': 80, 'A2': 80, 'A3': 110, 'A4': 40, CONTROL: 80}
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


def out_dir() -> Path:
    """<shared>/results/E136_stage_population_fit/[measure_<name>/] (created)."""
    d = H.results_dir(EXPERIMENT)
    if _measure_sub():
        d = d / _measure_sub()
        d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    """outputs/E136_stage_population_fit/[measure_<name>/][<part>/] (created)."""
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


def require_e135() -> dict:
    """E135's class and what it hands on to A1 (its table of PN rates by its route, as 'extra_kws'); RuntimeError while
    E135 has no class (A1 waits on it)."""
    raise RuntimeError('A1 waits on E135: what it hands on is not defined yet')


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
    flies) and the rest the held-out half, in the population's order."""
    from kcq import population_fit as P
    names = [str(f) for f in flies.flies]
    a, b = P.split_halves(len(names), 1, seed)[0]
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


def claw_thr_check(wirings: Sequence[str] = S_WIRINGS, value: float = CLAW_THR_CHECK_HZ) -> dict:
    """A3's mechanics: on each wiring at its held megamat point, a run with the per-claw threshold set in place on the
    simulator built without it (olfsysm's mp.kc.claw_thr) equals a fresh build at that threshold exactly, counts and
    every stage; written to claw_thr_check.json ('exact' True only if every wiring's is)."""
    from experiments import E134_population_baseline as B
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
    """The arm an arm name, the control (CONTROL) or the unrefitted family (FAMILY) is built as."""
    return 'A0' if arm in (CONTROL, FAMILY) else arm


def arm_kws(arm: str) -> dict:
    """The arm's fixed model kwargs beside the family's (make_fixed_param_kc_simulator kwargs; A1's table from E135
    and A4's KC-KC scale are not kwargs and are left out)."""
    kws = {k: v for k, v in ARMS[base_arm(arm)].kws.items() if k not in ('pn_rates', 'kckc_fixed')}
    if base_arm(arm) == 'A1':
        kws.update(require_e135()['extra_kws'])
    return kws


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

    def __init__(self, arm: str, wirings: Sequence[str] = S_WIRINGS, claw_mode: str = 'in_place',
                 lock: Optional[Mapping] = None):
        from experiments import E83_cross_wiring_fit as E83
        from experiments import E134_population_baseline as B
        self.arm, self.panel, self.wirings = arm, PANEL, tuple(wirings)
        a = ARMS[base_arm(arm)]
        self.searched = tuple(a.searched)
        if len(self.searched) + len(HELD) > MAX_FREE:
            raise ValueError(f'{arm}: more than {MAX_FREE} parameters')
        self.fixed = {'wKCKC_scale': float(a.kws['kckc_fixed'])} if 'kckc_fixed' in a.kws else {}
        self.kws = arm_kws(arm)
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
            self.built[w] = B.build_family(w, PANEL, dict(p), extra_kws=self.kws or None)
            self.built[w].info.update(arm=arm)
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
            kws = dict(self.kws, **({'claw_thr': float(value)} if value > 0 else {}))
            self.built[wiring] = B.build_family(wiring, PANEL, dict(self.base[wiring]), extra_kws=kws or None)
        self._claw_thr[wiring] = float(value)

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
    """<shared>/data/model_outputs/E136_stage_population_fit/[measure_<name>/]."""
    d = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return d / _measure_sub() if _measure_sub() else d


def fit_dirs(arm: str) -> List[Path]:
    """Where a fit's saved evaluations are: outputs/.../fit_<arm>/ on the host that ran it (every evaluation), then the
    archive's fit_<arm>/ (the best evaluation, `archive_best`)."""
    return [local_dir(f'fit_{arm}'), archive_dir() / f'fit_{arm}']


def archive_best(arm: str, i: int) -> List[Path]:
    """Copies a fit's best evaluation (every S wiring's files, '<arm>_e<i>_<wiring>_*') from this host's
    outputs/.../fit_<arm>/ to the archive's fit_<arm>/, so the score can read it on any host."""
    import shutil
    src, dst = local_dir(f'fit_{arm}'), archive_dir() / f'fit_{arm}'
    dst.mkdir(parents=True, exist_ok=True)
    out = []
    for w in S_WIRINGS:
        files = sorted(src.glob(f'{arm}_e{i}_{w}_*'))
        if not any(f.name.endswith('_run.json') for f in files):
            raise FileNotFoundError(f'{arm} eval {i} {w}: not saved under {src}')
        for f in sorted(files, key=lambda f: f.name.endswith('_run.json')):   # the run.json last, as it says complete
            tmp = dst / f'.{f.name}.{_host()}.{os.getpid()}.tmp'
            shutil.copyfile(f, tmp)
            os.replace(tmp, dst / f.name)
            out.append(dst / f.name)
    return out


def save_eval(arm: str, i: int, runs: Mapping[str, 'R.SavedRun']) -> List[Path]:
    """An evaluation's runs under outputs/.../fit_<arm>/ (E134's step-3 layout, `save_step3_run`; tag
    '<arm>_e<i>_<wiring>')."""
    d = local_dir(f'fit_{arm}')
    written = []
    for w, run in runs.items():
        info = {k: v for k, v in run.info.items() if not str(k).startswith('_')}
        info.update(tag=f'{arm}_e{i}_{w}', eval=i, host=_host(), time=_now())
        written += R.save_step3_run(f'{arm}_e{i}_{w}', run.info['_out'], run.meta, run.orn, info, [d])
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


class StageLoss:
    """The megamat stage error of a set of runs: per stage of `stages` (default readings, the flies of `flies` for
    reading), the population errors of the wirings' mean views (hal's `wirings_as_flies`, row 'common') and of each
    wiring against `flies`' references (PN dendrites on four views, the other stages on three, as E134's
    `megamat_stage_error`), their root mean square over the stages; for the control, the same under each permutation and
    their mean. Returns a dict: 'stage_error' (the recorded labels), '<stage>_error' and its views, each wiring's
    '<stage>_error_<wiring>' and 'stage_error_<wiring>', the error over the stages after the PN dendrites
    ('stage_error_downstream', and each wiring's), 'stage_error_perm', '<stage>_error_perm' and each wiring's
    'stage_error_perm_<wiring>' (the permutations' means, for the control), 'loss_stages' (the permutations' mean for the control, else 'stage_error') and the rates. With
    `check_pn`, the PN-dendrite error must be the same at every call."""

    def __init__(self, stages: Sequence[str], flies: Flies, lock: Optional[Mapping] = None,
                 subset_glomeruli: Optional[Sequence[str]] = None, check_pn: bool = False):
        self.stages = tuple(stages)
        bad = [s for s in self.stages if s not in LOSS_STAGES]
        if bad or not self.stages:
            raise ValueError(f'stages {self.stages}')
        self.flies = flies
        self.lock = lock if lock is not None else R.require_lock()
        self.subset_glomeruli = subset_glomeruli
        self.check_pn = check_pn
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
            out.update({f'{st}_{c}': float(tab.loc['common', c]) for c in tab.columns if c.startswith('view_')})
            common.append(out[f'{st}_error'])
            for w, pop in zip(runs, pops):
                out[f'{st}_error_{w}'] = float(tab.loc[pop.name, col])
                each[w].append(out[f'{st}_error_{w}'])
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

def identity(arm: str = 'A0', wirings: Sequence[str] = S_WIRINGS) -> dict:
    """At z = 0 (each wiring at its held megamat point of E134's step 3): the hold changes nothing (one simulation per
    rate, no factor), every run reproduces E134's saved counts and responses exactly, the in-memory runs read exactly as
    the saved runs do (every population's values and responses), and, against all flies, each wiring's population
    errors reproduce E134's score_scores.csv ('megamat_held', the default readings) to IDENTITY_TOL where that file
    exists; written to identity_<arm>.json; RuntimeError otherwise."""
    from experiments import E134_population_baseline as B
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


def fit_arm(arm: str, control: bool = False, maxfev: Optional[int] = None, resume: bool = True,
            save_stages: bool = True, claw_mode: Optional[str] = None) -> dict:
    """One arm's fit on the S wirings against the fit half (the control: A0 against the permuted fit half), as the
    module docstring says. Every simulated evaluation's loss, parts, rates, holds and parameters are appended to
    fit_history_<arm>.csv and its stages saved (`save_eval`); resumable: an evaluation whose clipped coordinates are in
    the history is not simulated again, and the holds' warm starts and slopes are restored from it, so a restarted fit
    retraces its path and continues. A3's per-claw threshold is set in place only where claw_thr_check.json says it
    equals a fresh build (`claw_mode` overrides). Writes fit_<arm>.json."""
    from scipy.optimize import minimize
    from data_prep import stage_population_data as D
    gate = require_gate()
    stages = readable_stages(gate['classes'])
    lock = R.require_lock()
    name = CONTROL if control else arm
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
    budget = int(maxfev or MAXFEV[name])
    hist_path = out_dir() / f'fit_history_{name}.csv'
    memo: Dict[Tuple[float, ...], float] = {}
    hist: List[dict] = []
    if resume and hist_path.exists():
        h = pd.read_csv(hist_path, float_precision='round_trip')
        for r in h.to_dict('records'):
            c = tuple(float(r[f'z_{k}']) for k in searched)
            memo[c] = float(r['loss_sim'])
            hist.append(r)
            for w in m.wirings:
                m.memory[w].append((c, float(r[f'x_thr_{w}']), float(r[f'x_apl_{w}'])))
                m.slopes[w] = dict(thr=float(r[f'slope_thr_{w}']), apl=float(r[f'slope_apl_{w}']))
        print(f'E136 {name}: resumed from {len(hist)} evaluations', flush=True)
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
            if i == 0 and name == 'A0' and not any(zc) and MEASURE in START_CHECK:
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
    restarted = False
    if left >= RESTART_MIN_LEFT:
        best = min(hist, key=lambda r: r['loss_sim'])
        x1 = np.array([float(best[f'z_{k}']) for k in searched])
        res2 = run_nm(x1, steps / 2.0, left)
        messages.append(str(res2.message))
        restarted = True
    best = min(hist, key=lambda r: r['loss_sim'])
    bz = {k: float(best[f'z_{k}']) for k in searched}
    flags = at_bound(searched, [bz[k] for k in searched])
    flags += [f'hold at an end on {w}' for w in S_WIRINGS if best.get(f'hold_at_end_{w}')]
    out = dict(arm=name, base_arm=base, searched=list(searched), stages=list(stages), start=START_GROUP, measure=MEASURE,
               gate_measure=gate['gate'].get('measure'), gate_verdict_as_read=gate['gate'].get('verdict_as_read'),
               n_evals=len(hist), n_calls=calls['n'], budget=budget, restarted=restarted, best_eval=int(best['eval']),
               best_z=bz, at_bound=flags, best=dict(best), messages=messages, claw_mode=m.claw_mode,
               seconds=round(time.time() - t0, 1), host=_host(), time=_now(), **_commits())
    out['factors'] = {k: (float(np.exp(v)) if k in FACTORS else CLAW_THR_UNIT * v) for k, v in bz.items()}
    if save_stages:
        out['archived'] = [str(p) for p in archive_best(name, int(best['eval']))]
    _write_json(out_dir() / f'fit_{name}.json', out)
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
    """An arm's fitted coordinates (fit_<arm>.json 'best_z'); {} for the unrefitted family. RuntimeError if the fit is
    not written."""
    if arm == FAMILY:
        return {}
    f = _read_json(out_dir() / f'fit_{arm}.json')
    if not f:
        raise RuntimeError(f'{arm} is not fitted yet (fit_{arm}.json)')
    return {k: float(v) for k, v in f['best_z'].items()}


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


def read_arm(pt: ReadPoint, lane: str, reclaim: bool = False, done: Optional[set] = None,
             archive: bool = True) -> Optional[dict]:
    """One read, claimed as E134's step-3 points are (mkdir of claims/<tag>/; None if done or claimed elsewhere): the
    arm on the point's wiring and panel at its fitted coordinates on the wiring's held megamat point, both rates held
    exactly at the panel's operating point (E77's `hold_rates`, the threshold with APL off, then the tied APL scales;
    `read_targets`), built again at the held parameters and run with every stage saved (which must give the hold's
    outputs exactly); A3's claws also read after the threshold (saved as '<tag>_clawafter'). Saved locally and in the
    archive (E134's step-3 layout); its row in reads_<host>_<lane>.csv and its provenance in
    provenance_read_<host>_<lane>.json. The claim is released on failure."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E134_population_baseline as B
    from experiments.E77_input_screen import frame_sha256, hold_rates
    tag = read_tag(pt)
    if done is not None and tag in done:
        return None
    if not B.claim(tag, lane, reclaim, root=out_dir()):
        print(f'E136 {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        z = fitted_z(pt.arm)
        arm = base_arm(pt.arm)
        a = ARMS[arm]
        lock = R.require_lock()
        p0, _ = start_params(pt.wiring)
        if 'kckc_fixed' in a.kws:
            p0['wKCKC_scale'] = float(a.kws['kckc_fixed'])
        for k, v in z.items():
            if k in FACTORS:
                for kk in FACTORS[k]:
                    p0[kk] = p0[kk] * float(np.exp(v))
        ct = CLAW_THR_UNIT * float(z['claw_thr']) if 'claw_thr' in z else None
        kws = dict(arm_kws(arm), **({'claw_thr': ct} if ct else {}))
        orn = B.panel_input(pt.panel)
        Pp = E83.wiring_parts(pt.wiring)
        t_off, t_on, f_off, f_on = read_targets(pt.wiring, pt.panel, orn.shape[1], lock)
        rng = HOLD_RANGE_DIVISIVE_APL if kws.get('apl_inh_divisive') else HOLD_RANGE
        hold = hold_rates(lambda: E83.build_simulator(pt.wiring, orn=orn, breadth=Pp['breadth'], params=dict(p0),
                                                      extra_kws=kws or None),
                          target_off=t_off, target_on=t_on, params=dict(p0), label=f'E136 {tag}', apl_bracket=rng)
        built = B.build_family(pt.wiring, pt.panel, hold['params'], orn=orn, extra_kws=kws or None)
        p, poff = dict(hold['params']), dict(hold['params_off'])
        reader = claw_reader_of(arm)
        out = B.run_stages(B.Point(PANEL, 'held', pt.wiring), built, p, poff, claw_reader=reader)
        for c in R.CONDITIONS:
            if B._differ(out[c]['spike_counts'], hold[c][1]) or B._differ(out[c]['responses'], hold[c][0]):
                raise RuntimeError(f"{tag} {c}: the saved run differs from the hold's outputs")
        info = dict(model=B.FAMILY, wiring=pt.wiring, panel=pt.panel, reading='read', arm=pt.arm, z=z, params=p,
                    params_off=poff, claw_thr=ct, tag=tag, host=_host(), lane=lane, time=_now())
        dirs = read_dirs() if archive else read_dirs()[:1]
        R.save_step3_run(tag, out, built.meta, built.orn, info, dirs)
        if ct:
            out2 = B.run_stages(B.Point(PANEL, 'held', pt.wiring), built, p, poff,
                                claw_reader=claw_reader_after_threshold(ct, bool(kws.get('apl_inh_divisive'))))
            R.save_step3_run(f'{tag}_clawafter', out2, built.meta, built.orn,
                             dict(info, tag=f'{tag}_clawafter', reading='read, claws after the threshold'), dirs)
        rr = {}
        for c in R.CONDITIONS:
            rr.update({f'rate_{c}': float((out[c]['responses'].to_numpy(dtype=float) > 0.5).mean()),
                       f'max_count_{c}': float(out[c]['spike_counts'].to_numpy(dtype=float).max())})
        row = dict(tag=tag, arm=pt.arm, panel=pt.panel, wiring=pt.wiring, half=pt.half,
                   **{f'z_{k}': v for k, v in z.items()}, claw_thr=ct,
                   **{k: float(v) for k, v in p.items()}, **{f'off_{k}': float(v) for k, v in poff.items()}, **rr,
                   target_on=t_on, target_off=t_off, frac_on_target=f_on, frac_off_target=f_off,
                   thr_factor=hold['thr_factor'], apl_factor=hold['apl_factor'], at_bound_thr=hold['at_bound_thr'],
                   at_bound_apl=hold['at_bound_apl'], hold_n_sims=hold['n_sims'], seconds=round(time.time() - t0, 1),
                   host=_host(), lane=lane, time=_now(), threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=frame_sha256(built.orn), **_commits())
        _append_row(out_dir() / f'reads_{_host()}_{lane}.csv', row)
        prov_path = out_dir() / f'provenance_read_{_host()}_{lane}.json'
        prov = _read_json(prov_path) or dict(part='read', lane=lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(info={k: v for k, v in built.info.items() if not str(k).startswith('_')},
                                   params=p, params_off=poff, model=H.model_provenance(orn_deltas=built.orn))
        _write_json(prov_path, prov)
        print(f"E136 {tag}: {row['seconds']:.0f} s, rate on {rr['rate_on']:.4f} / off {rr['rate_off']:.4f}", flush=True)
        return row
    except BaseException:
        B.release(tag, root=out_dir())
        raise


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

def load_read(pt: ReadPoint, claw_after: bool = False) -> Optional['R.SavedRun']:
    """A saved read (locally, else the archive; with `claw_after` A3's run with its claws read after the threshold,
    the same run elsewhere); the unrefitted family's megamat and mixture runs are E134's saved megamat_held and
    mixtures_held runs, and its APL16 runs E90's part B (the wirings' own fits); None if absent."""
    from experiments import E134_population_baseline as B
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
                  wirings: Sequence[str] = T_WIRINGS, perms: bool = False,
                  claw_after: bool = False) -> Optional[Dict[str, float]]:
    """An arm's (or the unrefitted family's) megamat stages on `wirings` against a half of the flies (all of them for
    None; under the control's permutations with `perms`): `StageLoss`'s dict. A fitted arm's S wirings are its fit's
    best evaluation (saved under outputs/.../fit_<arm>/). None while a run is missing."""
    from data_prep import stage_population_data as D
    runs = {}
    for w in wirings:
        if arm != FAMILY and w in S_WIRINGS:
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
    lock = R.require_lock()
    flies = Flies(half, stages, perms=permutations(list(D.panel_odors(PANEL))) if perms else None)
    return StageLoss(stages, flies, lock)(runs)


def score_panels(arm: str, flies: Optional['R.FlyCache'] = None) -> pd.DataFrame:
    """An arm's KC somata on the held-out panels (the validation2, mixture and APL16 inputs on the eight wirings, every
    panel and condition the runs are read on with flies to read them against, default readings), under the measure in
    force: per (input, panel, condition) each wiring's population error, their mean and the error of their mean views
    (hal's `wirings_as_flies` against all the panel's flies, E134's `fly_side_read`), with the views."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    flies = flies if flies is not None else R.FlyCache()
    lock = R.require_lock()
    views = MEASURES[MEASURE]
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
            for name, row in tab.iterrows():
                wiring = str(name).partition('/')[2] or str(name)
                rows.append(dict(arm=arm, input=panel, stage=st, panel=pan, condition=cond, wiring=wiring,
                                 measure=MEASURE, n_odors=len(od), n_flies=len(fp),
                                 **{k: float(v) for k, v in row.items()}))
    return pd.DataFrame(rows)


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


HOLDS_PANELS: Tuple[Tuple[str, str, str], ...] = (('validation2', 'validation2', 'on'), ('mixtures', 'kmix', 'on'),
                                                   ('mixtures', 'cmix', 'on'), ('apl16', 'apl16', 'on'),
                                                   ('apl16', 'apl16', 'off'))


def _sd(x) -> float:
    x = np.asarray(list(x), dtype=float)
    return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan


def megamat_values(sc: Mapping[str, float], stages: Sequence[str], wirings: Sequence[str]) -> Dict[str, Dict[str, float]]:
    """From `StageLoss`'s dict, per wiring its stage error ('stage'), its error over the stages after the PN dendrites
    ('downstream') and each stage's population error."""
    return {w: dict(stage=sc[f'stage_error_{w}'], downstream=sc.get(f'stage_error_downstream_{w}', np.nan),
                    **{st: sc[f'{st}_error_{w}'] for st in stages}) for w in wirings}


def panel_values(panels: pd.DataFrame, arm: str) -> Dict[Tuple[str, str, str], Dict[str, float]]:
    """From `score_panels`' rows, per (input, panel, condition) of HOLDS_PANELS each wiring's KC population error."""
    out = {}
    for inp, pan, cond in HOLDS_PANELS:
        x = panels[(panels.arm == arm) & (panels.input == inp) & (panels.panel == pan) & (panels.condition == cond)
                   & panels.wiring.isin(WIRINGS)] if len(panels) else panels
        out[(inp, pan, cond)] = {str(r.wiring): float(r.error) for r in x.itertuples()} if len(x) else {}
    return out


def compare(arm_v: Mapping[str, float], base_v: Mapping[str, float], spread: float,
            wirings: Sequence[str] = T_WIRINGS) -> dict:
    """The classes' comparison of one value as ruled: the per-wiring changes (arm minus the reference on the same
    wiring) over `wirings`, their mean, the standard error over wirings (the larger of half the spread and the SD of
    the changes over the wirings divided by sqrt(n)), and whether the change is an improvement (below -IMPROVE_SE
    standard errors) or a worsening beyond the spread."""
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
             mega_ctl: Optional[Mapping[str, Mapping[str, Dict[str, float]]]] = None) -> dict:
    """One arm's class as the entry fixes it and the rulings read it, on the T wirings and the held-out half, against
    the unrefitted family there (`mega` / `panels`: per arm, FAMILY included, the values of `megamat_values` /
    `panel_values`; `spreads`: 'stage', each stage and each HOLDS_PANELS key as 'input/panel/condition', the SD over the
    family's eight wirings; `realism`: per arm, per wiring its APL16 read's realism flag; `fit`: the arm's
    fit_<arm>.json; `control_fall`: the control fit's fall of the megamat stage error on T against the permuted
    held-out half). "Lower" in "one stage at another's cost" is read with the same margin as "improves"."""
    a, f = mega[arm], mega[FAMILY]
    m = compare({w: v['stage'] for w, v in a.items()}, {w: v['stage'] for w, v in f.items()}, spreads['stage'])
    st = {s: compare({w: v[s] for w, v in a.items()}, {w: v[s] for w, v in f.items()}, spreads[s]) for s in stages}
    pn = {}
    for key in HOLDS_PANELS:
        k = '/'.join(key)
        pn[k] = compare(panels[arm].get(key, {}), panels[FAMILY].get(key, {}), spreads.get(k, np.nan))
    rec = dict(arm=arm, megamat=m, stages=st, panels=pn)
    reasons = []
    if not m.get('complete') or not all(v.get('complete') for v in st.values()):
        reasons.append('reads missing')
    if fit.get('at_bound'):
        reasons.append('at a bound: ' + ', '.join(fit['at_bound']))
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
    improves = bool(m.get('improves')) and not stage_worse
    rec.update(improves=improves, stage_worse=stage_worse, panel_worse=panel_worse, panel_missing=panel_missing)
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


def part_score(arms: Optional[Sequence[str]] = None) -> dict:
    """The reads scored and classed: per arm (and the unrefitted family) the megamat stages on T against the held-out
    half (and, beside, against all flies and the fit half, and on S from the fit), A3's under both claw readings, the
    control's against the permuted held-out half, the held-out panels' KC population errors, the APL16 rows with the
    realism bands; the spreads over the family's eight wirings; the classes (`classify`). Writes score_megamat.csv,
    score_panels.csv, score_apl16.csv, classes.json and provenance_score.json."""
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
                              downstream=sc['stage_error_downstream'], **{s: sc[f'{s}_error'] for s in stages}))
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
    for key in HOLDS_PANELS:
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
    classes = {}
    for a in arms:
        fit = _read_json(out_dir() / f'fit_{a}.json')
        if a not in mega or FAMILY not in mega:
            classes[a] = dict(cls='reads missing')
            continue
        classes[a] = classify(a, stages, mega, pv, spreads, realism, fit, control_fall)
        if a in mega_after:
            alt = classify(a, stages, dict(mega, **{a: mega_after[a]}), pv, spreads, realism, fit, control_fall)
            classes[a]['claws_after_threshold'] = dict(cls=alt['cls'], megamat=alt['megamat'])
            if alt['cls'] != classes[a]['cls']:
                classes[a]['reasons'].append('the class differs between the claws read before and after the threshold')
                classes[a]['cls'] = 'no class'
    pd.DataFrame(mrows).to_csv(out_dir() / 'score_megamat.csv', index=False)
    panels.to_csv(out_dir() / 'score_panels.csv', index=False)
    apl.to_csv(out_dir() / 'score_apl16.csv', index=False)
    _write_json(out_dir() / 'classes.json', dict(stages=list(stages), spreads=spreads, control_fall=control_fall,
                                                 classes=classes))
    _write_json(out_dir() / 'provenance_score.json', dict(part='score', script=SCRIPT, lock=R.require_lock().get('sha256'),
                seconds=round(time.time() - t0, 1), host=_host(), time=_now(), arms=arms, **_commits()))
    return classes


# ---------------------------------------------------------------------------
# CLI

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('identity', 'claw_identity', 'claw_thr_check', 'fit', 'read', 'score'))
    ap.add_argument('--arm', choices=tuple(ARMS) + (CONTROL, FAMILY))
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
    elif a.part == 'score':
        print(json.dumps(part_score(None if a.arm is None else [a.arm]), indent=1, default=str))


if __name__ == '__main__':
    main()
