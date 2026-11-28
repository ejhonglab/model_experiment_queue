#!/usr/bin/env python3
"""E149, the megamat fits' picks under strength and sharing read where they were not fitted (experiment_priorities.md,
E149; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-28 04:46; disco's parts claimed
2026-11-28 05:06: the lock, the points' interface, step 0, the score and the write-up; hal's runner is
experiments/E149_reads.py).

Each of the four megamat fits (E136's A1, E145's X, F, W) is read, at fixed parameters, at the evaluation it would have
ended at under L1 (the megamat stage error on strength and sharing), on the wirings and the half of the flies it was not
fitted to, on the held-out panels and on the APL16 input, and scored by E136's classes with E139's price as E145 scored
a cell. Nothing is fitted, searched or weighted.

The points' interface (for hal's runner; every name below is in this module):
- The picks (`picks(lock=None) -> Dict[str, Pick]`, from the committed lock; `pick_table()` computes them from the saved
  tables, which the lock fixes). A pick's name is '<cell>e<eval>': step 1 'A1e46', 'Xe21', 'Fe86', 'We24' (L1's picks on
  E136's split, results/E148_megamat_other_views/losses_picks.csv); step 2's F L0 pick 'Fe116'; step 3's strength picks
  stopped early ('A1e45', 'Xe14', 'Fe80', 'We15' by the probe; the lock's `stopped`). `Pick` has: name, cell, eval, loss
  ('L1', 'L0' or 'L0_stopped'), z (E136's four coordinates aplpn, kckc, inhadd, inhsc: ln factors, from the cell's fit
  history), factors (exp z), pn_kws (E136's `pn_norm_kws(z)`: pn_inhadd, pn_inhsc), searched (the cell's four searched
  wirings), read (its four read wirings), fit_half and read_half ('a' or 'b'), arm ('A1', the E136 arm whose model it
  is), steps (the steps that read it).
- `Point(step, pick, panel, wiring)`: `.tag` (`point_tag`), `.cell`, `.eval`, `.half` (E136's 'held' for megamat, the
  cell's read half inside `E145.e136_cell(cell)`; None otherwise). Panels: 'megamat', 'validation2', 'mixtures', 'apl16'
  (E136's READ_PANELS) and 'extended' (E98's extended input, step 2).
- `point_tag(pick, panel, wiring)` = '<pick>_<panel>_<wiring>', e.g. 'Xe21_megamat_real'; `parse_tag(tag)` -> Point
  fields (step from the planned points).
- `planned_points(step, cells=None)` in a forward lane's order:
  - step 1 (112): megamat on each L1 pick's four read wirings (16; picks in A1, X, F, W order), then the APL16 input
    (32), validation2 (32) and the mixtures (32), each on all eight wirings (E136's WIRINGS);
  - step 2 (60): each L1 pick on the extended input on all eight wirings (32), then F's L0 pick 'Fe116' with E145's 28
    (megamat on F's read wirings, then APL16, validation2, mixtures on all eight);
  - step 3 (28 per cell): for each cell of `cells` (default the cells whose L1 pick holds, from the score's classes.json),
    its stopped strength pick with the 28 reads; RuntimeError before the score when `cells` is None.
- `read_extra(pt)`: the fields to pass as E136's `read_at(..., extra=...)` (step, pick, cell, eval, loss, searched,
  fit_half, read_half); E136's row then carries READ_COLUMNS.
- Folders: `results_dir()` (<shared>/results/E149_two_view_picks/; the rows reads_<host>_<lane>.csv and claims/ go here,
  as E136's `read_at(root=results_dir())` writes them), `local_dir(part='')` (outputs/E149_two_view_picks/),
  `archive_dir()` (<shared>/data/model_outputs/E149_two_view_picks/), `run_dirs(pt)` (outputs/.../step<k>/<pick>/<panel>/
  <wiring>/ and its archive twin, E134's step-3 layout, for E136's `read_at(dirs=...)`), `extended_dir(pick, where)`
  (step 2's extended runs in E102's part-B layout: outputs/.../step2/<pick>/extended/B/<E102 step tag of the pick's
  pn_inhadd>/, and the archive twin with where='archive', for E139's `read_beside` through a `BesideContext`),
  `identity_dir()` (results/.../identity/).
- The identity runs (step 2 of the entry's design): `identity_wiring(pick)` (the first of its searched wirings),
  `eval_params(pick, wiring)` (the saved evaluation's absolute parameters on that wiring from the fit history:
  fixed_thr, wAPLKC_scale, wKCAPL_scale, wPNAPL_scale, wAPLPN_scale, wKCKC_scale; APL off zeroes E136's OFF_ZERO),
  `eval_row(pick)` (the fit history's row: each stage's error per wiring and both response rates), `eval_tag(pick,
  wiring)` ('A1_e<eval>_<wiring>') and `eval_dirs(cell)` (where the evaluation's saved stages are).
- `reads_table()` (every lane's rows, the last per tag), `is_held(row)` (both counts within HOLD_TOL of their targets and
  neither factor at a bound), READ_COLUMNS (the columns the score reads), `require_lock()`.

Design (the entry's): step 0 (`--part step0`, no simulation: the identities, the saved settings' consistencies by the
fits' own reading, megamat's part of the L0 classed points' classes on E143's joint splits, a prediction per pick); the
identity runs, step 1's 112 reads and step 2's 60 (hal's runner, viki's and hal's lanes); the score (`--part score`):
each pick scored by E145's `score_cell` (E136's classes with E139's price, inside E145's `e136_cell`), the entry's
classes and outcome letters, the items beside; step 3 only for a cell whose L1 pick holds.

Outputs: <shared>/results/E149_two_view_picks/ (lock.json, step0.json, step0_*.csv, prediction.json, reads_<host>_<lane>
.csv, claims/, identity/, score_*.csv, classes.json, provenance_<part>.json); outputs/E149_two_view_picks/ (the runs);
data/model_outputs/E149_two_view_picks/ (the runs' archive twin).

Holds (CLAUDE.md, the entry): no olfsysm change (pn_inhadd and pn_inhsc are options at olfsysm b1e89d7; the wirings keep
'first' claws, pinned); nothing fitted, searched or weighted (the points are fixed by saved tables before any read); no
measure, class or yardstick changed; nothing fitted to a held-out panel, validation2 read as a held-out panel only;
Pratyush's hemibrain only; Remy's processed data without the bouton fly 2022-10-21 fly02 and her older megamat flies.

Usage (from the repository root, with -m):
    ./venv/bin/python -u -m experiments.E149_two_view_picks --part lock    # then commit experiments/E149_lock.json
    ./venv/bin/python -u -m experiments.E149_two_view_picks --part step0
    ./venv/bin/python -u -m experiments.E149_two_view_picks --part score
"""
from __future__ import annotations

import argparse
import json
import os
import re
import socket
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E149_two_view_picks'
SCRIPT = 'experiments/E149_two_view_picks.py'
RUNNER = 'experiments/E149_reads.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E149_lock.json'
PLANNER = 'a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-28 04:46'
E136 = 'E136_stage_population_fit'
E139 = 'E139_normalization_ladder'
E143 = 'E143_megamat_fly_sampling'
E145 = 'E145_stage_fit_replication'
E147 = 'E147_megamat_input_flies'
E148 = 'E148_megamat_other_views'
ARM = 'A1'
CELLS: Tuple[str, ...] = ('A1', 'X', 'F', 'W')
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
COORDS: Tuple[str, ...] = ('aplpn', 'kckc', 'inhadd', 'inhsc')
# E136's read panels and step 2's extended input (E98's, read by E102's rules as E139's `read_beside`)
READ_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
EXTENDED = 'extended'
PANELS: Tuple[str, ...] = READ_PANELS + (EXTENDED,)
# the panels after megamat in a pick's 28 reads, in a forward lane's order (the entry: megamat first, then the APL16
# input, validation2 and the mixtures)
HELD_OUT_ORDER: Tuple[str, ...] = ('apl16', 'validation2', 'mixtures')
STEPS: Tuple[int, ...] = (1, 2, 3)
# F's L0 pick as picked (step 2, beside): E145's classed point for F returned APL->bouton to the family's
F_L0_EVAL = 116
# "held" (E136's HOLD_TOL, E147's rule): both counts within this relative miss of their targets, no factor at a bound
HOLD_TOL = 0.005
# the yardsticks (E136's, as scored; the four-factor control and E147's spreads beside)
CONTROL_FALL = 0.040
IMPROVE_SE = 2.0
FOUR_FACTOR_CONTROL_FALL_FILE = ('control.json',)      # E145's, beside
# the columns of a read's row the score reads (E136's `read_at` writes them, with `read_extra`'s fields)
READ_COLUMNS: Tuple[str, ...] = (
    'tag', 'step', 'pick', 'cell', 'eval', 'loss', 'panel', 'wiring', 'half', 'read_half', 'searched', 'fit_half',
    'z_aplpn', 'z_kckc', 'z_inhadd', 'z_inhsc', 'pn_inhadd', 'pn_inhsc',
    'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale',
    'off_fixed_thr', 'off_wAPLKC_scale', 'off_wAPLPN_scale', 'off_wKCKC_scale',
    'rate_on', 'rate_off', 'max_count_on', 'max_count_off', 'target_on', 'target_off', 'frac_on_target',
    'frac_off_target', 'thr_factor', 'apl_factor', 'at_bound_thr', 'at_bound_apl',
    'seconds', 'host', 'lane', 'time', 'queue_commit', 'queue_dirty', 'pyolfsysm_commit', 'pyolfsysm_dirty')


# ---------------------------------------------------------------------------
# folders

DIRS: Dict[str, Path] = {}


def results_dir() -> Path:
    """<shared>/results/E149_two_view_picks/ (created); DIRS['results'] in tests."""
    d = Path(DIRS['results']) if 'results' in DIRS else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    """outputs/E149_two_view_picks/[<part>/] in this checkout (created)."""
    root = Path(DIRS['local']) if 'local' in DIRS else H.outputs_dir(EXPERIMENT)
    d = root / part if part else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E149_two_view_picks/ (not created)."""
    return Path(DIRS['archive']) if 'archive' in DIRS else H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def identity_dir() -> Path:
    """results/E149_two_view_picks/identity/ (created): the identity runs' records."""
    d = results_dir() / 'identity'
    d.mkdir(parents=True, exist_ok=True)
    return d


def lock_path() -> Path:
    return results_dir() / LOCK_NAME


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


# ---------------------------------------------------------------------------
# the picks

@dataclass(frozen=True)
class Pick:
    """One point read at fixed parameters: a saved evaluation of one cell's fit.

    Attributes:
        name: '<cell>e<eval>'.
        cell: 'A1', 'X', 'F' or 'W' (E145's CELLS).
        eval: the evaluation in the cell's fit history.
        loss: 'L1' (L1's pick on E136's split), 'L0' (F's L0 pick as picked) or 'L0_stopped' (the strength pick stopped
            early, step 3).
        z: E136's four coordinates (ln factors on each wiring's held megamat point).
        pn_kws: E136's `pn_norm_kws(z)` (pn_inhadd, pn_inhsc).
        searched, read: the cell's searched and read wirings.
        fit_half, read_half: 'a' or 'b'.
        steps: the steps that read it.
    """
    name: str
    cell: str
    eval: int
    loss: str
    z: Mapping[str, float]
    pn_kws: Mapping[str, float]
    searched: Tuple[str, ...]
    read: Tuple[str, ...]
    fit_half: str
    read_half: str
    steps: Tuple[int, ...]
    arm: str = ARM

    @property
    def factors(self) -> Dict[str, float]:
        return {k: float(np.exp(v)) for k, v in self.z.items()}


def pick_name(cell: str, i: int) -> str:
    if cell not in CELLS:
        raise ValueError(f'{cell=}')
    return f'{cell}e{int(i)}'


def cell_results(cell: str) -> Path:
    """A cell's fit records (E148's `cell_results`: E136's measure_strength/ for A1, E145's <cell>/ otherwise)."""
    from experiments import E148_megamat_other_views as M
    return M.cell_results(cell)


def history(cell: str) -> pd.DataFrame:
    """A cell's fit history (fit_history_A1.csv), exact floats."""
    return pd.read_csv(cell_results(cell) / f'fit_history_{ARM}.csv', float_precision='round_trip')


def make_pick(cell: str, i: int, loss: str, steps: Sequence[int], hist: Optional[pd.DataFrame] = None) -> Pick:
    """A pick from a cell's fit history: its z at evaluation `i` and E136's PN terms."""
    from experiments import E136_stage_population_fit as F
    from experiments import E145_stage_fit_replication as E
    h = history(cell) if hist is None else hist
    r = h[h['eval'] == int(i)]
    if len(r) != 1:
        raise ValueError(f'{cell}: evaluation {i} not in its fit history')
    z = {k: float(r.iloc[0][f'z_{k}']) for k in COORDS}
    c = E.get_cell(cell)
    return Pick(name=pick_name(cell, i), cell=cell, eval=int(i), loss=loss, z=z, pn_kws=F.pn_norm_kws(z),
                searched=tuple(c.searched_wirings()), read=tuple(c.read_wirings()), fit_half=c.fit_half,
                read_half=c.read_half, steps=tuple(int(s) for s in steps))


def stopped_eval(cell: str, k: int) -> int:
    """The strength pick stopped early: the held evaluation with the lowest L0 (the rescore's megamat stage error on
    strength) against the cell's fit half of E136's split among evaluations 0 to k."""
    from experiments import E148_megamat_other_views as M
    c = M.get_cell(cell)
    rs = M.rescore_table(cell)
    ev = M.cell_evals(cell).set_index('eval')
    s = M._tab(rs, 'strength', c.fit_half)
    idx = [int(i) for i in s.index if int(i) <= int(k) and bool(ev.loc[int(i), 'held'])]
    return int(s.loc[idx, 'stage_error'].idxmin())


def pick_table() -> Dict[str, Pick]:
    """The picks from the saved tables (what the lock fixes): L1's picks per cell on E136's split from E148's
    losses_picks.csv (step 1, and step 2's extended input), F's L0 pick as picked (step 2), and each cell's strength pick
    stopped early at its L1 pick's index (step 3)."""
    lp = pd.read_csv(H.results_dir(E148) / 'losses_picks.csv')
    out: Dict[str, Pick] = {}
    for c in CELLS:
        x = lp[(lp['cell'] == c) & (lp['loss'] == 'L1') & (lp['split'] == 'e136')]
        if len(x) != 1:
            raise RuntimeError(f"E148's losses_picks.csv has no single L1 pick for {c}")
        i = int(x['pick'].iloc[0])
        p = make_pick(c, i, 'L1', (1, 2))
        out[p.name] = p
    f0 = lp[(lp['cell'] == 'F') & (lp['loss'] == 'L0') & (lp['split'] == 'e136')]
    if int(f0['pick'].iloc[0]) != F_L0_EVAL:
        raise RuntimeError(f"F's L0 pick in E148's table is {int(f0['pick'].iloc[0])}, not {F_L0_EVAL}")
    p = make_pick('F', F_L0_EVAL, 'L0', (2,))
    out[p.name] = p
    for c in CELLS:
        k = next(q.eval for q in out.values() if q.cell == c and q.loss == 'L1')
        p = make_pick(c, stopped_eval(c, k), 'L0_stopped', (3,))
        out[p.name] = p
    return out


def _pick_record(p: Pick) -> dict:
    return dict(name=p.name, cell=p.cell, eval=p.eval, loss=p.loss, z=dict(p.z), factors=p.factors,
                pn_kws=dict(p.pn_kws), searched=list(p.searched), read=list(p.read), fit_half=p.fit_half,
                read_half=p.read_half, steps=list(p.steps), arm=p.arm)


def _pick_of(r: Mapping) -> Pick:
    return Pick(name=r['name'], cell=r['cell'], eval=int(r['eval']), loss=r['loss'],
                z={k: float(v) for k, v in r['z'].items()}, pn_kws={k: float(v) for k, v in r['pn_kws'].items()},
                searched=tuple(r['searched']), read=tuple(r['read']), fit_half=r['fit_half'], read_half=r['read_half'],
                steps=tuple(int(s) for s in r['steps']), arm=r.get('arm', ARM))


def picks(lock: Optional[Mapping] = None) -> Dict[str, Pick]:
    """The picks as the committed lock fixes them (`require_lock` unless `lock` is given)."""
    lock = require_lock() if lock is None else lock
    return {r['name']: _pick_of(r) for r in lock['picks']}


def l1_picks(ps: Mapping[str, Pick]) -> List[Pick]:
    """The four L1 picks, in CELLS' order."""
    return [next(p for p in ps.values() if p.cell == c and p.loss == 'L1') for c in CELLS]


# ---------------------------------------------------------------------------
# points and tags

@dataclass(frozen=True)
class Point:
    """One read: a pick on one panel's input and one wiring, in one step."""
    step: int
    pick: str
    panel: str
    wiring: str

    @property
    def tag(self) -> str:
        return point_tag(self.pick, self.panel, self.wiring)

    @property
    def cell(self) -> str:
        return parse_pick(self.pick)[0]

    @property
    def eval(self) -> int:
        return parse_pick(self.pick)[1]

    @property
    def half(self) -> Optional[str]:
        """E136's half name for `read_at` inside the cell: 'held' (the cell's read half) on megamat, None otherwise."""
        return 'held' if self.panel == 'megamat' else None


def parse_pick(name: str) -> Tuple[str, int]:
    m = re.fullmatch(r'(A1|X|F|W)e(\d+)', str(name))
    if not m:
        raise ValueError(f'{name!r} is not a pick name')
    return m.group(1), int(m.group(2))


def point_tag(pick: str, panel: str, wiring: str) -> str:
    """'<pick>_<panel>_<wiring>', e.g. 'Xe21_megamat_real' (ValueError for an unknown part)."""
    parse_pick(pick)
    if panel not in PANELS or wiring not in WIRINGS:
        raise ValueError(f'{pick=} {panel=} {wiring=}')
    return f'{pick}_{panel}_{wiring}'


def parse_tag(tag: str) -> Tuple[str, str, str]:
    """(pick, panel, wiring) of a tag (`point_tag`'s inverse)."""
    m = re.fullmatch(r'((?:A1|X|F|W)e\d+)_(' + '|'.join(PANELS) + r')_([A-Za-z0-9]+)', str(tag))
    if not m:
        raise ValueError(f'{tag!r} is not a tag of {EXPERIMENT}')
    pick, panel, wiring = m.groups()
    point_tag(pick, panel, wiring)
    return pick, panel, wiring


def pick_points(step: int, p: Pick) -> List[Point]:
    """A pick's 28 reads (E145's): megamat on its four read wirings, then the APL16 input, validation2 and the mixtures
    on all eight wirings."""
    pts = [Point(step, p.name, 'megamat', w) for w in p.read]
    pts += [Point(step, p.name, pan, w) for pan in HELD_OUT_ORDER for w in WIRINGS]
    return pts


def planned_points(step: int, cells: Optional[Sequence[str]] = None, lock: Optional[Mapping] = None) -> List[Point]:
    """A step's reads in a forward lane's order (the module docstring): step 1 the four L1 picks, megamat first (16),
    then APL16, validation2 and the mixtures (32 each), 112; step 2 the L1 picks on the extended input (32), then F's L0
    pick's 28, 60; step 3 for each cell of `cells` (default those whose L1 pick holds, from classes.json) its stopped
    strength pick's 28."""
    ps = picks(lock)
    step = int(step)
    if step == 1:
        l1 = l1_picks(ps)
        pts = [Point(1, p.name, 'megamat', w) for p in l1 for w in p.read]
        pts += [Point(1, p.name, pan, w) for pan in HELD_OUT_ORDER for p in l1 for w in WIRINGS]
        return pts
    if step == 2:
        pts = [Point(2, p.name, EXTENDED, w) for p in l1_picks(ps) for w in WIRINGS]
        f0 = next(p for p in ps.values() if p.loss == 'L0')
        return pts + pick_points(2, f0)
    if step == 3:
        if cells is None:
            cl = json.loads((results_dir() / 'classes.json').read_text()) if (results_dir() / 'classes.json').exists() \
                else None
            if cl is None:
                raise RuntimeError('step 3 waits for the score (classes.json): it reads only cells whose L1 pick holds')
            cells = [c for c in CELLS if ((cl.get('picks') or {}).get(next(p.name for p in l1_picks(ps)
                                                                           if p.cell == c)) or {}).get('cls') == 'holds']
        pts = []
        for c in cells:
            p = next(q for q in ps.values() if q.cell == c and q.loss == 'L0_stopped')
            pts += pick_points(3, p)
        return pts
    raise ValueError(f'{step=}: 1, 2 or 3')


def point_step(tag: str, lock: Optional[Mapping] = None) -> int:
    """The step a tag belongs to (step 3's for its stopped picks, which no other step reads)."""
    pick, panel, _ = parse_tag(tag)
    p = picks(lock)[pick]
    if panel == EXTENDED:
        return 2
    return int(p.steps[0]) if p.loss != 'L1' else 1


def read_extra(pt: Point, lock: Optional[Mapping] = None) -> dict:
    """The fields of a read's row beside E136's `read_at` own (its `extra`): step, pick, cell, eval, loss, searched ('S'
    or 'T'), fit_half and read_half (the cell's read half on megamat, None otherwise)."""
    p = picks(lock)[pt.pick]
    return dict(step=int(pt.step), pick=p.name, cell=p.cell, eval=p.eval, loss=p.loss,
                searched='S' if p.searched == S_WIRINGS else 'T', fit_half=p.fit_half,
                read_half=p.read_half if pt.panel == 'megamat' else None)


def run_dirs(pt: Point) -> List[Path]:
    """Where a read's runs are written and found: outputs/E149_two_view_picks/step<k>/<pick>/<panel>/<wiring>/ (created)
    and its archive twin (not created; E136's `read_at` writes into both)."""
    rel = Path(f'step{int(pt.step)}') / pt.pick / pt.panel / pt.wiring
    return [local_dir(str(rel)), archive_dir() / rel]


def extended_dir(pick: str, where: str = 'local', lock: Optional[Mapping] = None) -> Path:
    """Step 2's extended runs of a pick in E102's part-B layout (as E139's `extended_dir`):
    outputs/.../step2/<pick>/extended/B/<E102's step tag of the pick's pn_inhadd>/ ('local', created) or its archive
    twin ('archive')."""
    from experiments import E102_antennal_lobe_normalization as E102
    if where not in ('local', 'archive'):
        raise ValueError(f'{where=}')
    p = picks(lock)[pick]
    rel = Path('step2') / pick / EXTENDED
    root = local_dir(str(rel)) if where == 'local' else archive_dir() / rel
    return root / 'B' / E102.step_tag(float(p.pn_kws['pn_inhadd']))


def load_read(pt: Point):
    """A read's saved runs (E134's step-3 layout, `run_dirs`), or None."""
    from kcq import stage_readings as R
    try:
        return R.load_step3_run(pt.tag, dirs=run_dirs(pt))
    except FileNotFoundError:
        return None


# ---------------------------------------------------------------------------
# the identity runs' inputs (hal's runner)

def identity_wiring(p: Pick) -> str:
    """The searched wiring a pick's identity run is made on (the first of its cell's searched wirings)."""
    return p.searched[0]


def eval_row(p: Pick) -> pd.Series:
    """The pick's row of its cell's fit history (each stage's error per wiring, both response rates, the parameters)."""
    h = history(p.cell)
    return h[h['eval'] == p.eval].iloc[0]


def eval_params(p: Pick, wiring: str) -> Dict[str, float]:
    """The saved evaluation's absolute parameters on one of its searched wirings (APL on), from the fit history:
    fixed_thr, wAPLKC_scale, wKCAPL_scale, wPNAPL_scale, wAPLPN_scale, wKCKC_scale; APL off zeroes E136's OFF_ZERO."""
    from experiments import E136_stage_population_fit as F
    if wiring not in p.searched:
        raise ValueError(f'{p.name}: {wiring} is not one of its searched wirings')
    r = eval_row(p)
    return {k: float(r[f'{k}_{wiring}']) for k in F.VARIANT_PARAMS}


def eval_tag(p: Pick, wiring: str) -> str:
    """The saved evaluation's stages' tag ('A1_e<eval>_<wiring>')."""
    return f'{ARM}_e{p.eval}_{wiring}'


def eval_dirs(cell: str) -> List[Path]:
    """Where a cell's evaluations' stages are (this checkout's outputs, then the archive; E148's `_fit_dirs`)."""
    from experiments import E148_megamat_other_views as M
    return M._fit_dirs(cell)


# ---------------------------------------------------------------------------
# the reads' rows

def reads_table() -> pd.DataFrame:
    """Every lane's rows (results/.../reads_*.csv), the last per tag."""
    fs = [pd.read_csv(f) for f in sorted(results_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=list(READ_COLUMNS))
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def _truthy(x) -> bool:
    return str(x).strip().lower() in ('true', '1', '1.0')


def is_held(row: Mapping) -> bool:
    """A read is held: both response rates within HOLD_TOL of their targets (relative) and neither factor at a bound."""
    try:
        on = abs(float(row['rate_on']) / float(row['frac_on_target']) - 1.0) <= HOLD_TOL
        off = abs(float(row['rate_off']) / float(row['frac_off_target']) - 1.0) <= HOLD_TOL
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
    return bool(on and off and not _truthy(row.get('at_bound_thr')) and not _truthy(row.get('at_bound_apl')))


# ---------------------------------------------------------------------------
# the lock

def lock_contents() -> dict:
    """The lock: the picks (from the saved tables), the steps' points, the reading, the yardsticks, the classes, the
    outcomes, the prediction's events, the judgments and the files with their sha256."""
    raise NotImplementedError


def part_lock() -> dict:
    """Writes the lock to results/.../lock.json and LOCK_COPY (commit the copy before any other part)."""
    raise NotImplementedError


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (LOCK_COPY committed with no uncommitted change, the same bytes as the shared copy), its files
    checked; PermissionError otherwise (E147's)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 0 and the score

def part_step0() -> dict:
    """Step 0 (no simulation): the identities (E148's picks table, the probe's rows, E145's four classes reproduced by
    this score from E145's saved reads); the saved settings' consistencies by the fits' own reading on the saved splits
    and on every split's halves read anew, on E143's splits and its 200 joint splits; megamat's part of the four L0
    classed points' classes on those joint splits; the prediction per pick (prediction.json, once, before any read)."""
    raise NotImplementedError


def score_pick(p: Pick) -> dict:
    """One pick scored as E145 scored a cell (E145's `score_cell` on the pick's reads inside E145's `e136_cell`): E136's
    classes with E139's price against the unrefitted family where the cell is read, rule 3's guard, the held-out panels;
    with the entry's per-pick statistics."""
    raise NotImplementedError


def part_score() -> dict:
    """The score (designs 5 and 6, step 3's condition): every pick of steps 1 and 2 whose reads are in, the entry's
    classes and outcome letters (classes.json), the tables (score_*.csv) and the items beside."""
    raise NotImplementedError


PARTS: Tuple[str, ...] = ('lock', 'step0', 'score')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    a = ap.parse_args(argv)
    fn = globals()[f'part_{a.part}']
    print(json.dumps(fn(), indent=1, default=str)[:20000])


if __name__ == '__main__':
    main()
