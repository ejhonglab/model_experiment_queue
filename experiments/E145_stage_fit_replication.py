#!/usr/bin/env python3
"""E145, the fit of all megamat stages made again on the other wirings and the other halves of the flies
(experiment_priorities.md, E145; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 19:30; its
code claimed by disco 2026-11-26 20:07).

E136's A1 (the best family with the antennal lobe's normalization searched beside the APL->bouton and KC-KC factors)
was fitted once: on the S wirings against half a of each stage's flies, read on the T wirings against half b. This fits
A1's form again, unchanged, in the three other cells of wirings by halves of the flies, each read where it was not
fitted, so that A1's fall, its "holds" and its fitted factors can be found a second time where nothing of them was
chosen.

The cells (CELLS): a cell is a set of four searched wirings with one half of each stage's flies.
- A1 (E136's, done): S against half a, read on T against half b;
- X, the exchange: T against half b, read on S against half a;
- F, the other flies: S against half b, read on T against half a;
- W, the other wirings: T against half a, read on S against half b.
"Half a" is E136's fit half (the first part of hal's `split_halves` at seed 0: 3 of 7 dendrite flies, 2 of 5 bouton
flies, 2 of 4 claw flies, 2 of 4 soma flies) and "half b" its held-out half. Inside a cell E136's 'fit' and 'held' name
the cell's own halves (E136's `cell`); this module's records name the halves 'a' and 'b'.

Everything else is E136's, unchanged: the measure (strength by odor alone, the adopted one), the loss (the megamat
stage error on the searched wirings' mean views against the fit half), the operating point and its holds, the start at
the family's held megamat points, A1's four coordinates, bounds, first steps, tolerances, its budget of 140 evaluations
with one restart, the single-factor returns with their 0.010 rule, the rescore, the reads and the classes. E136's own
machinery runs each cell (`e136_cell`: E136's `cell`, `paths` and measure), so that no fit is written twice.

Design (the entry's):
- step 0 (no simulation, `--part step0`): the lock (`--part lock`, before anything else); the family's held megamat
  runs and A1's saved runs scored where each cell is read (`cell_table`); each cell's start values from saved runs
  (`start_checks`: the family at z = 0 on the cell's searched wirings against its fit half, which the first evaluation
  must reproduce; E136's own cell is checked against E136's START_CHECK); A1's 124 saved evaluations against both
  halves, the ten best by each with their factors (`a1_evals`); the prediction per cell, written before any fit;
- step 1 (`--part fit --cell X|F|W`): per cell on one lane, E136's identity at z = 0 on the searched wirings (every run
  reproduces E134's held megamat runs with 0 differing counts), the fit (resumable from its history), the single-factor
  returns, the rescore against both halves under every measure with the order between the halves, then the cell's 28
  reads;
- the reads (`--part read`): per fitted cell, megamat on the four read wirings against the cell's read half, and the
  validation2, mixture and APL16 inputs with APL on and off on all eight wirings, at the classed point;
- the score, E143's resampled reading per cell and the four-factor control fit (`--part score`, `--part resampled`,
  `--part control`) are written once the fits run (stubs here).

Outputs: <shared>/results/E145_stage_fit_replication/ (E145_lock.json, step0.json, step0_*.csv, prediction.json, and per
cell <cell>/: E136's fit records, identity_A1.json, determined_A1.csv, rescore_A1.csv, consistency_*.csv, claims/,
reads_<host>_<lane>.csv, provenance_read_<host>_<lane>.json); outputs/E145_stage_fit_replication/<cell>/ (every
evaluation's stages, fit_A1/; the reads' runs, read/); data/model_outputs/E145_stage_fit_replication/<cell>/ (each
fit's best evaluation, fit_A1/, and the reads, read/).

Holds (CLAUDE.md, the entry): Pratyush's hemibrain only; no olfsysm change (pn_inhadd and pn_inhsc are options at
olfsysm b1e89d7; the wirings keep 'first' claws, pinned); four global factors searched, none per odor or glomerulus;
nothing fitted to a held-out panel; Remy's processed data without the bouton fly 2022-10-21 fly02 and without her older
megamat flies (E134's lock). Cells X and W search the T wirings, which every fit so far left unsearched: they are this
entry's cells only.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E145_stage_fit_replication --part lock | --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E145_stage_fit_replication --part fit --cell X --lane <host>[-<n>]
    ... --part read --cell X F W --lane <host>[-<n>] [--reverse]
"""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E145_stage_fit_replication'
SCRIPT = 'experiments/E145_stage_fit_replication.py'
LOCK_NAME = 'E145_lock.json'
PLANNER = 'a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 19:30'
# the measure every fit and class reads (E134's adopted one) and the arm whose form is fitted
MEASURE = 'strength'
ARM = 'A1'
PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
HALF_NAMES: Tuple[str, str] = ('a', 'b')
# the start check's tolerance (E136's START_CHECK_TOL; the values here are exact, E136's were rounded to 4 decimals)
START_TOL = 6e-5
# the ten best evaluations by each half (step 0)
N_BEST = 10


@dataclass(frozen=True)
class Cell:
    """One cell of wirings by halves of the flies.

    Attributes:
        name: 'A1' (E136's, done), 'X', 'F' or 'W'.
        label: what it is, in words.
        searched: 'S' or 'T', the wirings the fit searches (E136's S_WIRINGS or T_WIRINGS); the other four are read.
        fit_part: 0 if the fit half is half a (E136's fit half), 1 if half b (its held-out half); the read half is the
            other.
    """
    name: str
    label: str
    searched: str
    fit_part: int

    def searched_wirings(self) -> Tuple[str, ...]:
        from experiments import E136_stage_population_fit as F
        return tuple(F.WIRINGS[:4] if self.searched == 'S' else F.WIRINGS[4:])

    def read_wirings(self) -> Tuple[str, ...]:
        from experiments import E136_stage_population_fit as F
        return tuple(w for w in F.WIRINGS if w not in self.searched_wirings())

    @property
    def fit_half(self) -> str:
        return HALF_NAMES[self.fit_part]

    @property
    def read_half(self) -> str:
        return HALF_NAMES[1 - self.fit_part]

    @property
    def read_on(self) -> str:
        return 'T' if self.searched == 'S' else 'S'


CELLS: Dict[str, Cell] = {c.name: c for c in (
    Cell('A1', "E136's A1 (done): S against half a, read on T against half b", 'S', 0),
    Cell('X', 'the exchange: T against half b, read on S against half a', 'T', 1),
    Cell('F', 'the other flies: S against half b, read on T against half a', 'S', 1),
    Cell('W', 'the other wirings: T against half a, read on S against half b', 'T', 0),
)}
# the cells this experiment fits, in the entry's order (cell X first, then F, then W)
FIT_CELLS: Tuple[str, ...] = ('X', 'F', 'W')

# "a factor repeats" in a cell (the entry's definitions; the gain on pn_inhsc is not read)
REPEAT = dict(
    inhadd=dict(lo=562.0, hi=1778.0, unit='pn_inhadd',
                rule="pn_inhadd at the classed point between 562 and 1,778, nearer 1000 than E139's rungs at 314 and "
                     "3141"),
    aplpn=dict(lo=0.078, hi=0.31, unit='factor', determined=True,
               rule='the APL->bouton factor determined (its single-factor return raises the loss by 0.010 or more) and '
                    'between x0.078 and x0.31'),
    kckc=dict(lo=1.8, hi=None, unit='factor', determined=True,
              rule="the KC-KC factor determined and at x1.8 or more (the lower end of E141's range)"),
)

CLASSES: Dict[str, str] = {
    'per_cell': "E136's classes with E139's price, per cell, against the unrefitted family where the cell is read "
                "(its read wirings against its read half); A1's row beside as scored.",
    'fall_repeats': 'in all three cells the megamat stage error falls from the family\'s by more than the bar and by '
                    'more than 0.040 (the control fit\'s fall).',
    'megamat_fixes_a_factor': 'the factor repeats (REPEAT) in all three cells.',
    'hold_repeats': '"holds" as scored in cell X and in at least one of F and W; each with its share over the '
                    'resamples, a share under 0.8 named a near-tie.',
    'the_flies': 'cell F differs from A1 in its class or in whether KC-KC repeats, and cell W does not.',
    'the_wirings': 'cell W differs from A1 in its class or in whether KC-KC repeats, and cell F does not.',
    'no_class': 'for a cell whose fit ends at a bound or whose holds fail, or for a wiring outside a realism band; the '
                'number of cells is reported with any that holds.',
}
OUTCOMES: Dict[str, str] = {
    'a': 'The fall, the hold and all three factors repeat (about 1 in 5): A1 has been found twice where nothing of it '
         'was chosen; to Tom and Yang as a model that can be built on, with E141\'s question of the dose.',
    'b': 'The fall and the two factors before the calyx repeat, and KC-KC does not in two or more cells, which then read '
         '"megamat only" (about 1 in 3): megamat fixes the normalization and APL->bouton and nothing after them.',
    'c': 'The hold repeats at other factors (about 1 in 5): a larger additive term, or KC-KC elsewhere in E141\'s range; '
         'the fitted model is a range.',
    'd': 'The cells disagree (about 1 in 4): "the flies\'" or "the wirings\'"; A1\'s hold stays a near-tie.',
    'e': 'The fall does not repeat in a cell (about 1 in 20): to Tom and Yang at once.',
    'f': 'A fit ends at a bound (about 1 in 10, beside the others): that cell has no class; a planner rules on the bound '
         'before any wider search.',
}
PREDICTION: Dict[str, str] = {
    'fall_repeats_all_three': 'about 9 in 10',
    'holds_in_X': 'about 2 in 5',
    'holds_in_all_three': 'about 1 in 6',
    'kckc_determined_at_1.8_or_more_in_all_three': 'about 2 in 5',
    'kckc_determined_at_1.8_or_more_in_none': 'about 1 in 6',
    'source': "the planner's chances (the entry's), written before any fit",
}
JUDGMENTS: Tuple[str, ...] = (
    "A cell's fit is E136's `fit_arm` of A1 in its search mode inside E136's `cell` and `paths`, so that E136's code, "
    "not a copy, fits each cell; 'fit' and 'held' in its records are the cell's own halves.",
    "A cell's first evaluation (the family at z = 0) is checked against the start values step 0 reads from E134's "
    "saved held megamat runs on the cell's searched wirings against its fit half, to START_TOL (E136's START_CHECK_TOL); "
    "E136's own cell is read the same way and must reproduce E136's START_CHECK, which checks the reading.",
    "The 0 differing counts at z = 0 are E136's `identity` of A1 on the cell's searched wirings, run before the fit "
    "(identity_A1.json in the cell's folder): every run reproduces E134's saved counts and responses exactly.",
    "A fitted point is the cell's classed point (E136's `determined`: 'classed_z', else the best point), as E136 reads "
    "its arms.",
    "Reads are tagged '<cell>_A1_<panel>_<wiring>' and claimed in the cell's folder, so that no tag is shared between "
    "cells or with E136.",
    "The repeat of the additive term reads pn_inhadd at the classed point (E136's `pn_norm_kws`); the entry's 562 and "
    "1,778 are taken as written.",
    "The entry says the gain on pn_inhsc is not read; E136's returns found it determined in A1 (rise 0.064), and it is "
    "reported beside, reading nothing.",
)
# amendments to the lock after it is written, appended (none yet)
AMENDMENTS: Tuple[dict, ...] = ()


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir(cell: Optional[str] = None) -> Path:
    """<shared>/results/E145_stage_fit_replication/[<cell>/] (created)."""
    d = H.results_dir(EXPERIMENT)
    if cell:
        d = d / cell
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(cell: Optional[str] = None, part: str = '') -> Path:
    """outputs/E145_stage_fit_replication/[<cell>/][<part>/] (created)."""
    d = H.outputs_dir(EXPERIMENT)
    for x in (cell, part):
        if x:
            d = d / x
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir(cell: Optional[str] = None) -> Path:
    """<shared>/data/model_outputs/E145_stage_fit_replication/[<cell>/]."""
    d = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return d / cell if cell else d


def read_dirs(cell: str) -> List[Path]:
    """Where a cell's reads' runs are: outputs/.../<cell>/read/, then the archive's <cell>/read/."""
    return [local_dir(cell, 'read'), archive_dir(cell) / 'read']


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


def get_cell(name: str) -> Cell:
    if name not in CELLS:
        raise ValueError(f'cell {name!r} not in {sorted(CELLS)}')
    return CELLS[name]


@contextlib.contextmanager
def e136(measure: str = MEASURE) -> Iterator[object]:
    """E136's module with its measure set for the block (restored after)."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(measure)
    try:
        yield F
    finally:
        F.set_measure(old)


@contextlib.contextmanager
def e136_cell(name: str, start_check: Optional[Mapping[str, float]] = None) -> Iterator[object]:
    """E136 on the adopted measure inside a cell (E136's `cell`: the cell's searched and read wirings, its fit half as
    'fit') with its fit records, evaluations and archived best evaluation in this experiment's folders for the cell
    (E136's `paths`), so that E136's `identity`, `fit_arm`, `determined`, `rescore`, `consistency`, `fitted_z`,
    `hold_at` and `read_at` work in the cell. `start_check` as E136's `cell` takes it."""
    c = get_cell(name)
    with e136() as F, F.cell(c.searched_wirings(), c.fit_part, start_check), \
            F.paths(out=out_dir(name), local=local_dir(name), archive=archive_dir(name)):
        yield F


# ---------------------------------------------------------------------------
# the lock

def lock_dict() -> dict:
    """The lock: the cells, the fit as E136 fixes it (A1's coordinates, bounds, first steps, tolerances, budget, restart
    and the single-factor rule), the reads, the yardsticks, the repeats, the classes, the outcomes, the prediction, the
    judgments and the amendments, with the hashes of what it reads (E134's lock, E136's A1 fit and classes)."""
    import hashlib
    from kcq import stage_readings as R
    from experiments import E136_stage_population_fit as F
    with e136() as F:
        a1 = F.out_dir() / f'fit_{ARM}.json'
        cl = F.out_dir() / 'classes.json'
        control_fall = float(_read_json(cl)['control_fall'])
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    lk = R.require_lock()
    coords = ('aplpn', 'kckc', 'inhadd', 'inhsc')
    return dict(
        experiment=EXPERIMENT, script=SCRIPT, part='lock', planner=PLANNER,
        entry='experiment_priorities.md item E145 (proposed and moved by viki 2026-11-26 20:05)',
        cells={k: dict(label=c.label, searched=c.searched, searched_wirings=list(c.searched_wirings()),
                       read_on=c.read_on, read_wirings=list(c.read_wirings()), fit_half=c.fit_half,
                       read_half=c.read_half, fitted_here=k in FIT_CELLS) for k, c in CELLS.items()},
        halves=dict(a="E136's fit half: the first part of hal's split_halves(n, 1, seed 0) of each stage's flies",
                    b="E136's held-out half: the rest", seed=F.HALF_SEED),
        fit=dict(arm=ARM, measure=MEASURE, coordinates=list(coords),
                 bounds={k: [float(np.exp(v)) for v in F.BOUNDS[k]] for k in coords},
                 first_steps={k: float(F.FIRST_STEP[k]) for k in coords}, xatol=F.XATOL, fatol=F.FATOL,
                 budget=int(F.MAXFEV_A1_SEARCH), restart_min_left=int(F.RESTART_MIN_LEFT),
                 not_determined=float(F.NOT_DETERMINED), hold_tol=float(F.HOLD_TOL),
                 hold_max_sims=int(F.HOLD_MAX_SIMS), hold_first_max_sims=int(F.HOLD_FIRST_MAX_SIMS),
                 hold_miss_per=float(F.HOLD_MISS_PER), bound_penalty_per=float(F.BOUND_PENALTY_PER),
                 start="the family's held megamat points of the cell's searched wirings (E134's step 3), every factor 1",
                 start_check_tol=START_TOL),
        reads=dict(per_cell=len(planned_points(['X'])), panels=list(PANELS),
                   rule='megamat on the four read wirings against the read half; validation2, the mixtures and the '
                        'APL16 input with APL on and off on all eight wirings; at the classed point, both rates held '
                        "at each panel's operating point (E136's read_at)"),
        yardsticks=dict(bar_se=float(F.IMPROVE_SE), bar="twice the standard error of a change over the cell's four "
                        "read wirings", spreads="the family's over its eight wirings, as E136 scored them",
                        control_fall=control_fall),
        repeat=REPEAT, classes=CLASSES, outcomes=OUTCOMES, prediction=PREDICTION, judgments=list(JUDGMENTS),
        resampled=dict(n=1000, read_once_check=0.02, source="E143's reading, per cell"),
        control=dict(cell='X', arm="A1's form (its four coordinates)", permutations=int(F.N_PERMS),
                     perm_seed=int(F.PERM_SEED), budget=80, gates=False),
        e134_lock_sha256=lk.get('sha256'), e136_fit_a1_sha256=sha(a1), e136_classes_sha256=sha(cl),
        amendments=list(AMENDMENTS))


VOLATILE = ('written', 'host', 'commits', 'amended', 'amended_host', 'amended_commits')


def lock_path() -> Path:
    return out_dir() / LOCK_NAME


def committed_lock_path() -> Path:
    return H.HERE / 'experiments' / LOCK_NAME


def _strip(d: Mapping, also: Sequence[str] = ()) -> dict:
    return json.loads(json.dumps({k: v for k, v in d.items() if k not in VOLATILE + tuple(also)}, default=str))


def _same_lock(a: Mapping, b: Mapping) -> bool:
    return _strip(a) == _strip(b)


def _amended_only(cur: Mapping, new: Mapping) -> bool:
    """Whether `new` is the lock `cur` with amendments appended to its own and nothing else changed."""
    a = json.loads(json.dumps(list(cur.get('amendments') or []), default=str))
    b = json.loads(json.dumps(list(new.get('amendments') or []), default=str))
    return _strip(cur, ('amendments',)) == _strip(new, ('amendments',)) and len(b) > len(a) and b[:len(a)] == a


def any_fit_started() -> bool:
    """Whether any cell's fit has a history or a claim."""
    return any((out_dir(c) / f'fit_history_{ARM}.csv').exists() or (out_dir(c) / 'claims').exists() for c in FIT_CELLS)


def write_lock(force: bool = False) -> dict:
    """results/E145_stage_fit_replication/E145_lock.json and its committed copy experiments/E145_lock.json (`lock_dict`
    with the time, host and commits). A lock already written is kept when it is the same; a different one is refused
    once any fit has started, unless `force` or the only change is amendments appended (then written with its first
    writing's time, host and commits kept and the amendment's beside)."""
    cur = _read_json(lock_path())
    new = lock_dict()
    if cur and _same_lock(cur, new):
        return cur
    if cur and _amended_only(cur, new):
        new = dict(new, written=cur.get('written'), host=cur.get('host'), commits=cur.get('commits'), amended=_now(),
                   amended_host=_host(), amended_commits=H.code_commits())
    elif cur and not force and any_fit_started():
        raise RuntimeError(f"{lock_path()} differs from the code's lock and a fit has started: add an amendment")
    else:
        new = dict(new, written=_now(), host=_host(), commits=H.code_commits())
    _write_json(lock_path(), new)
    committed_lock_path().write_text(json.dumps(new, indent=1, default=str) + '\n')
    return new


def require_lock() -> dict:
    """The lock as written (results), checked against the code's (`lock_dict`) and its committed copy; RuntimeError
    otherwise."""
    cur = _read_json(lock_path())
    if not cur:
        raise RuntimeError(f'no lock yet ({lock_path()}): run --part lock first')
    com = _read_json(committed_lock_path())
    if not com or not _same_lock(cur, com):
        raise RuntimeError(f'the lock differs from its committed copy {committed_lock_path()}')
    if not _same_lock(cur, lock_dict()):
        raise RuntimeError("the lock differs from the code's definitions (lock_dict)")
    return cur


# ---------------------------------------------------------------------------
# step 0 (no simulation)

def half_of(name: str) -> str:
    """E136's name, outside any cell, of half 'a' ('fit') or 'b' ('held')."""
    return {'a': 'fit', 'b': 'held'}[name]


def wirings_of(which: str) -> Tuple[str, ...]:
    from experiments import E136_stage_population_fit as F
    return {'S': tuple(F.WIRINGS[:4]), 'T': tuple(F.WIRINGS[4:])}[which]


def cell_table() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Step 0's table from saved runs: the unrefitted family (E134's held megamat runs) and E136's A1 (its fit's best
    evaluation on S, its reads on T) on each set of wirings (S, T) against each half (a, b), on the adopted measure:
    per wiring its megamat stage error, its error after the PN dendrites and each stage's population error, and their
    mean over the four wirings ('mean', the classes' reading: E136's `compare` takes the mean of the per-wiring values),
    and beside it the error of the wirings' mean views ('common', the loss's reading: a fit's loss and START_CHECK).
    Returns (rows: arm, wirings, half, wiring (a wiring, 'mean' or 'common'), stage, downstream and each stage; the
    cells: per cell, the family's and A1's values where the cell is read and where it is fitted, in both readings,
    '_common' for the wirings' mean views)."""
    rows = []
    with e136() as F:
        stages = F.readable_stages(F.require_gate()['classes'])
        for arm in (F.FAMILY, ARM):
            for which in ('S', 'T'):
                ws = wirings_of(which)
                for h in HALF_NAMES:
                    sc = F.score_megamat(arm, stages, half_of(h), ws)
                    if sc is None:
                        raise FileNotFoundError(f'{arm} on {which}: a saved megamat run is missing')
                    vals = F.megamat_values(sc, stages, ws)
                    for w, v in vals.items():
                        rows.append(dict(arm=arm, wirings=which, half=h, wiring=w, stage=v['stage'],
                                         downstream=v['downstream'], **{s: v[s] for s in stages}))
                    rows.append(dict(arm=arm, wirings=which, half=h, wiring='mean',
                                     stage=float(np.mean([v['stage'] for v in vals.values()])),
                                     downstream=float(np.mean([v['downstream'] for v in vals.values()])),
                                     **{s: float(np.mean([v[s] for v in vals.values()])) for s in stages}))
                    rows.append(dict(arm=arm, wirings=which, half=h, wiring='common', stage=sc['stage_error'],
                                     downstream=sc['stage_error_downstream'],
                                     **{s: sc[f'{s}_error'] for s in stages}))
    df = pd.DataFrame(rows)
    m = df[df.wiring == 'mean'].set_index(['arm', 'wirings', 'half'])
    cm = df[df.wiring == 'common'].set_index(['arm', 'wirings', 'half'])
    cells = []
    for k, c in CELLS.items():
        r = dict(cell=k, label=c.label, searched=c.searched, fit_half=c.fit_half, read_on=c.read_on,
                 read_half=c.read_half)
        for tab, suf in ((m, ''), (cm, '_common')):
            r.update({f'family_where_read{suf}': float(tab.loc[('family', c.read_on, c.read_half), 'stage']),
                      f'a1_where_read{suf}': float(tab.loc[(ARM, c.read_on, c.read_half), 'stage']),
                      f'family_where_fitted{suf}': float(tab.loc[('family', c.searched, c.fit_half), 'stage']),
                      f'a1_where_fitted{suf}': float(tab.loc[(ARM, c.searched, c.fit_half), 'stage'])})
        r['a1_where_read_is_its_own_fit_wirings'] = bool(c.read_on == 'S')
        cells.append(r)
    return df, pd.DataFrame(cells)


def start_values(name: str) -> Dict[str, float]:
    """A cell's start values: the unrefitted family at z = 0 (E134's saved held megamat runs) on the cell's searched
    wirings against its fit half, as E136's fit reads its first evaluation (`StageLoss` on the readable stages, the
    measure in force): START_CHECK's keys (the loss and each stage's error)."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F0
    from kcq import stage_readings as R
    keys = list(F0.START_CHECK[MEASURE])
    with e136_cell(name) as F:
        stages = F.readable_stages(F.require_gate()['classes'])
        lock = R.require_lock()
        runs = {w: B.load_run(B.tag_of(B.Point(F.PANEL, 'held', w))) for w in F.S_WIRINGS}
        missing = [w for w, r in runs.items() if r is None]
        if missing:
            raise FileNotFoundError(f"E134's held megamat runs are missing on {missing}")
        res = F.StageLoss(stages, F.Flies('fit', stages), lock)(runs)
    return {k: float(res[k]) for k in keys}


def start_checks() -> Dict[str, dict]:
    """Every cell's start values (`start_values`); E136's own cell (A1) must reproduce E136's START_CHECK to START_TOL
    (it checks the reading). RuntimeError otherwise."""
    from experiments import E136_stage_population_fit as F
    out = {k: start_values(k) for k in CELLS}
    ref = F.START_CHECK[MEASURE]
    bad = {k: (out['A1'][k], v) for k, v in ref.items() if abs(out['A1'][k] - v) > START_TOL}
    if bad:
        raise RuntimeError(f"E136's own cell does not reproduce E136's START_CHECK (got, expected): {bad}")
    return out


def a1_evals() -> Tuple[pd.DataFrame, dict]:
    """A1's saved evaluations (E136's rescore_A1.csv on the adopted measure, against half a and half b, and its
    fit_history_A1.csv) with their factors: per evaluation its megamat stage error and each stage's against each half,
    its four factors (pn_inhadd as an absolute value), and whether both rates held on every S wiring. Returns (the
    table, a summary: the ten best by each half with their factors, and the Spearman ρ between the halves over the
    evaluations whose rates held on every wiring)."""
    from scipy import stats
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    with e136() as F:
        d = F.out_dir()
        rs = pd.read_csv(d / f'rescore_{ARM}.csv')
        hist = pd.read_csv(d / f'fit_history_{ARM}.csv', float_precision='round_trip')
        ws = list(F.S_WIRINGS)
    rs = rs[rs.measure == MEASURE]
    cols = ['stage_error', 'pn_dendrite_error', 'pn_bouton_error', 'kc_claw_error', 'kc_soma_error']
    a = rs[rs.half == 'fit'].set_index('eval')[cols].add_suffix('_a')
    b = rs[rs.half == 'held'].set_index('eval')[cols].add_suffix('_b')
    h = hist.set_index('eval')
    t = a.join(b, how='inner')
    t['aplpn'] = np.exp(h.loc[t.index, 'z_aplpn'])
    t['kckc'] = np.exp(h.loc[t.index, 'z_kckc'])
    t['inhadd_factor'] = np.exp(h.loc[t.index, 'z_inhadd'])
    t['pn_inhadd'] = PN_INHADD_DEFAULT * t['inhadd_factor']
    t['inhsc_gain'] = np.exp(h.loc[t.index, 'z_inhsc'])
    t['held_all'] = h.loc[t.index, [f'hold_ok_{w}' for w in ws]].astype(bool).all(axis=1)
    t['loss_sim'] = h.loc[t.index, 'loss_sim']
    t = t.reset_index()
    ok = t[t.held_all]
    rho = float(stats.spearmanr(ok['stage_error_a'], ok['stage_error_b'])[0]) if len(ok) > 2 else np.nan
    keep = ['eval', 'stage_error_a', 'stage_error_b', 'aplpn', 'kckc', 'inhadd_factor', 'pn_inhadd', 'inhsc_gain']
    summ = dict(n_evals=int(len(t)), n_held_all=int(len(ok)), rho_halves_held=rho,
                best_by_a=ok.nsmallest(N_BEST, 'stage_error_a')[keep].to_dict('records'),
                best_by_b=ok.nsmallest(N_BEST, 'stage_error_b')[keep].to_dict('records'))
    for k in ('best_by_a', 'best_by_b'):
        x = pd.DataFrame(summ[k])
        summ[f'{k}_ranges'] = {c: [float(x[c].min()), float(x[c].max())] for c in keep[1:]}
    return t, summ


def part_step0() -> dict:
    """Step 0 (no simulation): the lock checked, the cell table (step0_cells.csv, step0_table.csv), the start values
    (in step0.json), A1's evaluations (step0_a1_evals.csv) and the prediction (prediction.json), written before any
    fit; refuses once a fit has started, unless the files are already the same."""
    lock = require_lock()
    t0 = time.time()
    rows, cells = cell_table()
    checks = start_checks()
    ev, summ = a1_evals()
    if any_fit_started() and (out_dir() / 'step0.json').exists():
        old = _read_json(out_dir() / 'step0.json')
        if json.loads(json.dumps(old.get('start_checks'), default=str)) != json.loads(json.dumps(checks, default=str)):
            raise RuntimeError("a fit has started and step 0's start values would change")
    rows.to_csv(out_dir() / 'step0_table.csv', index=False)
    cells.to_csv(out_dir() / 'step0_cells.csv', index=False)
    ev.to_csv(out_dir() / 'step0_a1_evals.csv', index=False)
    rec = dict(part='step0', lock_written=lock.get('written'), start_checks=checks, start_check_tol=START_TOL,
               cells=cells.to_dict('records'), a1_evals=summ, seconds=round(time.time() - t0, 1), host=_host(),
               time=_now(), **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _write_json(out_dir() / 'step0.json', rec)
    pred = dict(PREDICTION, cells={r['cell']: dict(family_where_read=r['family_where_read'],
                                                   a1_where_read=r['a1_where_read']) for r in rec['cells']},
                written=_now(), host=_host(), before_any_fit=not any_fit_started())
    if not (out_dir() / 'prediction.json').exists():
        _write_json(out_dir() / 'prediction.json', pred)
    return rec


def require_step0() -> dict:
    """step0.json, with a start value for every cell; RuntimeError before step 0."""
    s = _read_json(out_dir() / 'step0.json')
    if not s or set(s.get('start_checks') or {}) < set(FIT_CELLS):
        raise RuntimeError('step 0 is not written (step0.json): run --part step0 first')
    return s


# ---------------------------------------------------------------------------
# step 1: the fits

@dataclass(frozen=True)
class Point:
    """A cell's read: A1 at the cell's classed point on one wiring and one panel's input."""
    cell: str
    panel: str
    wiring: str


def point_tag(pt: Point) -> str:
    """'<cell>_A1_<panel>_<wiring>', e.g. 'X_A1_megamat_real'."""
    return f'{pt.cell}_{ARM}_{pt.panel}_{pt.wiring}'


def parse_tag(tag: str) -> Point:
    cell, arm, panel, wiring = str(tag).split('_', 3)
    if cell not in FIT_CELLS or arm != ARM or panel not in PANELS:
        raise ValueError(f'{tag!r} is not a read of {EXPERIMENT}')
    return Point(cell, panel, wiring)


def planned_points(cells: Sequence[str] = FIT_CELLS) -> List[Point]:
    """A cell's 28 reads in a forward lane's order (E136's `planned_reads` of an arm): megamat on its four read
    wirings, then the validation2, mixture and APL16 inputs on the eight wirings."""
    from experiments import E136_stage_population_fit as F
    pts = []
    for k in cells:
        c = get_cell(k)
        pts += [Point(k, 'megamat', w) for w in c.read_wirings()]
        pts += [Point(k, p, w) for p in PANELS if p != 'megamat' for w in F.WIRINGS]
    return pts


def reads_table(cell: str) -> pd.DataFrame:
    """A cell's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir(cell).glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def fit_record(cell: str) -> dict:
    """A cell's fit_A1.json ({} before its fit ends)."""
    return _read_json(out_dir(cell) / f'fit_{ARM}.json')


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read, claimed in the cell's folder: E136's `read_at` of A1 at the cell's classed point (`fitted_z` inside
    the cell) on the point's wiring and panel, both rates held at the panel's operating point, saved in the cell's
    read folders; its row (with the cell and its halves) in the cell's reads_<host>_<lane>.csv. None if done or
    claimed elsewhere; RuntimeError before the cell's single-factor returns."""
    c = get_cell(pt.cell)
    tag = point_tag(pt)
    if done is not None and tag in done:
        return None
    if 'classed_z' not in fit_record(pt.cell):
        raise RuntimeError(f"cell {pt.cell}'s fit has no classed point yet (its single-factor returns)")
    half = c.read_half if pt.panel == 'megamat' else None
    extra = dict(cell=pt.cell, searched=c.searched, fit_half=c.fit_half, read_half=half, step='1')
    with e136_cell(pt.cell) as F:
        return F.read_at(tag, ARM, lambda: F.fitted_z(ARM), pt.wiring, pt.panel, 'held' if half else None, lane,
                         out_dir(pt.cell), read_dirs(pt.cell), reclaim=reclaim, done=done, extra=extra, label='E145',
                         script=SCRIPT)


def part_read(lane: str, cells: Sequence[str] = FIT_CELLS, reverse: bool = False, reclaim: bool = False) -> None:
    """A lane of reads: every planned read of the cells whose classed point is written, skipping what is done or
    claimed elsewhere."""
    ready = [k for k in cells if 'classed_z' in fit_record(k)]
    pts = planned_points(ready)
    for pt in (pts[::-1] if reverse else pts):
        read_point(pt, lane, reclaim, set(reads_table(pt.cell).tag))


FIT_PARTS: Tuple[str, ...] = ('identity', 'fit', 'determined', 'rescore', 'reads')


def part_fit(cell: str, lane: str, reclaim: bool = False, parts: Sequence[str] = FIT_PARTS) -> dict:
    """Step 1 for one cell, on one lane: claimed in the cell's folder ('fit_<cell>'; a fit resumed on another lane
    passes `reclaim`), then, each part skipped where its record is complete: E136's identity of A1 at z = 0 on the
    searched wirings (identity_A1.json: 0 differing counts against E134's held megamat runs), the fit (E136's `fit_arm`
    of A1, resumable from its history, the first evaluation checked against step 0's start values), the single-factor
    returns (`determined`), the rescore against both halves under every measure (`rescore`) and the order between the
    halves (`consistency`), then the cell's 28 reads. The claim is released when the lane stops on an error, so that
    another lane can resume the fit. Returns a summary."""
    from experiments import E134_population_baseline as B
    if cell not in FIT_CELLS:
        raise ValueError(f'cell {cell!r}: E145 fits {FIT_CELLS} (A1 is E136\'s)')
    require_lock()
    checks = require_step0()['start_checks']
    root = out_dir(cell)
    tag = f'fit_{cell}'
    if not B.claim(tag, lane, reclaim, root=root):
        raise RuntimeError(f'cell {cell} is claimed by another lane ({root}/claims/{tag}); pass --reclaim to resume it')
    out: Dict[str, object] = dict(cell=cell, lane=lane, host=_host())
    try:
        with e136_cell(cell, checks[cell]) as F:
            if 'identity' in parts and not (root / f'identity_{ARM}.json').exists():
                rec = F.identity(ARM, F.S_WIRINGS)
                out['identity'] = dict(exact=rec['exact'], seconds=rec['seconds'])
            if 'fit' in parts and not (root / f'fit_{ARM}.json').exists():
                if not (root / f'identity_{ARM}.json').exists():
                    raise RuntimeError(f'cell {cell}: no identity_{ARM}.json, the fit waits for the identity')
                F.fit_arm(ARM)
            fit = fit_record(cell)
            out['fit'] = {k: fit.get(k) for k in ('best_eval', 'factors', 'at_bound', 'n_evals', 'pn_kws')} \
                if fit else None
            if 'determined' in parts and fit and 'determined' not in fit:
                fit = F.determined(ARM)
            if 'rescore' in parts and fit:
                rp = root / f'rescore_{ARM}.csv'
                have = set(pd.read_csv(rp)['eval']) if rp.exists() else set()
                if set(F.history_evals(ARM)) - have:
                    F.rescore(ARM)
                    for meas in F.MEASURES:
                        F.consistency([ARM], meas)
        out['classed_z'] = fit_record(cell).get('classed_z')
        if 'reads' in parts and 'classed_z' in fit_record(cell):
            part_read(lane, [cell], reclaim=reclaim)
            out['reads'] = int(len(reads_table(cell)))
    except BaseException:
        B.release(tag, root=root)
        raise
    return out


# ---------------------------------------------------------------------------
# the score, the resampled reading and the control fit (written once the fits run)

def part_score() -> dict:
    """Step 5: per cell, E136's classes with E139's price against the unrefitted family where the cell is read, with
    A1's row beside as scored; the experiment's classes (CLASSES) and outcome. Not written yet."""
    raise NotImplementedError("E145's score is written while the fits run")


def part_resampled() -> dict:
    """Beside: E143's resampled reading of every cell (1,000 resamples, each class recomputed, after its read-once
    check to 0.02). Not written yet."""
    raise NotImplementedError("E145's resampled reading is written while the fits run")


def part_control(lane: str, reclaim: bool = False) -> dict:
    """Beside: A1's form fitted in cell X to flies with permuted odor labels (E136's 20 permutations, 80 evaluations),
    its fall where cell X is read against the permuted half; gates nothing. Not written yet (it needs E136's control
    on A1's form)."""
    raise NotImplementedError("E145's four-factor control fit is written while the fits run")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'step0', 'fit', 'read', 'score', 'resampled', 'control'))
    ap.add_argument('--cell', nargs='+', default=list(FIT_CELLS), choices=FIT_CELLS)
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--fit-parts', nargs='+', default=list(FIT_PARTS), choices=FIT_PARTS)
    a = ap.parse_args(argv)
    if a.part == 'lock':
        print(json.dumps(write_lock(), indent=1, default=str)[:3000])
    elif a.part == 'step0':
        print(json.dumps(part_step0(), indent=1, default=str))
    elif a.part == 'fit':
        if len(a.cell) != 1:
            raise SystemExit('--part fit takes one --cell')
        print(json.dumps(part_fit(a.cell[0], a.lane, a.reclaim, a.fit_parts), indent=1, default=str))
    elif a.part == 'read':
        part_read(a.lane, a.cell, a.reverse, a.reclaim)
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'resampled':
        print(json.dumps(part_resampled(), indent=1, default=str))
    elif a.part == 'control':
        print(json.dumps(part_control(a.lane, a.reclaim), indent=1, default=str))


if __name__ == '__main__':
    main()
