#!/usr/bin/env python3
"""E148, the megamat fits' saved evaluations read on sharing and breadth and under other losses (experiment_priorities.md,
E148; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-27 19:55; code claimed by disco
2026-11-27 20:21).

Nothing is run, fitted, searched or chosen. The four fits of all megamat stages (E136's A1 and E145's X, F, W) saved every
evaluation's stages on their four searched wirings, and E143 read 23 saved settings on T. Their views are read once per
cell (`--part views --cell <cell>`), and every later part reads only those views: which views and which losses megamat's
flies can decide, what each loss would have picked and how far that pick is from the other half's best, whose claws the
claws' cost is, and whether the megamat stage error without the boutons still orders the settings.

Definitions (the entry's, fixed in the lock before any statistic beyond the planner's probe; experiments/E148_lock.json):
- the evaluations (CELLS): A1's 124 (E136, measure strength), X's 139, F's 126 and W's 139 (E145), each on its four
  searched wirings; an evaluation whose response rates were not held on every searched wiring (its fit history's
  hold_ok_<wiring>) is left out and counted; an evaluation's value at a stage is the error of its wirings' mean views
  (row 'common', the fit's loss reading, rescore's 'stage_error'), each wiring's beside;
- the saved settings: E143's 23 (role 'setting'), with its two groups that share a PN stage (CALYX_GROUPS: nine at
  pn_inhadd 1000, five at the family's normalization), read on T (a value the mean over T of each wiring's error, as
  E143 reads them); the unrefitted family on all eight wirings and E145's cells' classed points (E145_A1, _X, _F, _W:
  each fit's best evaluation on its searched wirings, its reads elsewhere) beside;
- the views: strength by odor (rate, mean), sharing between odors (who responds, how strongly) and breadth among units,
  with the matrix view at the PN dendrites where E134 locks four views, as hal's kcq/population_fit.py computes them;
- the reading: each run read once through E134's default reading with the curve, noise and criterion of all the stage's
  flies, seed index the wiring's position in E136's WIRINGS (E143's `read_runs`); a subset of the flies changes only the
  reference (its mean view, in all the stage's flies' one-fly units). Beside, E136's own reading (the half's flies set the
  curve, noise and criterion) on E136's split at every stage and on every pair of the claws' and somata's four flies
  (their three splits), saved with the views; where the read-once check (E143's CHECK_TOL, 0.02) misses in a cell, that
  reading is used there;
- the splits: every split of a stage's flies as E143's `stage_splits` (35, 10, 3, 3), E143's 200 joint splits (seed 0)
  with E136's beside, and each claw fly and each soma fly alone as the reference;
- the losses (LOSSES), each the root mean square of its terms, a term a stage's error over a set of views;
- a pick: the held evaluation with the lowest loss against the fit half; its regret: against the other half and under the
  same loss, its distance from that half's lowest held evaluation, as a share of the distance from the start (evaluation
  0, the family at z = 0) to that lowest;
- megamat's part of "holds" at an evaluation against a half, on the searched wirings (HOLDS_PART): a fall of the megamat
  stage error on strength from evaluation 0 by more than the control fit's 0.040, no stage's strength error up by more
  than its spread, and the claws' and somata's errors over all views not up by more than theirs (E136's spreads as
  scored); a description, since a class is scored on unsearched wirings with the held-out panels;
- the classes (CLASSES) and outcomes (OUTCOMES), fixed beforehand.

Design (CLI `--part`):
- lock: the lock, written to results/E148_megamat_other_views/lock.json and experiments/E148_lock.json (commit the copy
  before any other part);
- views --cell A1|X|F|W|settings: every evaluation's (or setting's) views per stage, wiring and reading, written to
  results/E148_megamat_other_views/views/<cell>.parquet (with the flies' parts, flies__<key>.parquet) on the host that
  holds the cell's stages; it says which evaluations it could not find (views/<cell>_summary.json);
- step0: the planner's probe reproduced from E136's and E145's rescore tables and fit histories, the prediction per stage
  and view (written once, before any views are read), and the identities on every cell whose views exist (the rescore
  tables reproduced from the views by E136's own reading, the consistencies, E143's consistency table on strength, the
  read-once check);
- room (design 3), views_decided (4), losses (5), claws (6), no_boutons (7), beside (8), score (classes.json with the
  classes and outcome letters).

Outputs: <shared>/results/E148_megamat_other_views/ (lock.json, views/, step0.json, step0_*.csv, prediction.json,
identities.csv, room.csv, views_decided.csv, losses*.csv, claws*.csv, no_boutons*.csv, beside_*.csv, classes.json,
provenance_<part>.json); outputs/E148_megamat_other_views/ (logs).

Holds (CLAUDE.md, the entry): no simulation, no lane, no fit, no olfsysm change; nothing fitted, searched or weighted; no
measure, class or yardstick changed; nothing fitted to a held-out panel, validation2's saved scores read beside only;
Pratyush's hemibrain only; Remy's processed data without the bouton fly 2022-10-21 fly02 and her older megamat flies.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E148_megamat_other_views --part lock    # then commit experiments/E148_lock.json
    OMP_NUM_THREADS=1 nice -n 10 ./venv/bin/python -u -m experiments.E148_megamat_other_views --part views --cell F
    ./venv/bin/python -u -m experiments.E148_megamat_other_views --part step0
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E148_megamat_other_views'
SCRIPT = 'experiments/E148_megamat_other_views.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E148_lock.json'
PLANNER = 'a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-27 19:55'
E134 = 'E134_population_baseline'
E136 = 'E136_stage_population_fit'
E139 = 'E139_normalization_ladder'
E141 = 'E141_kckc_rise_forms'
E143 = 'E143_megamat_fly_sampling'
E145 = 'E145_stage_fit_replication'
E147 = 'E147_megamat_input_flies'
PANEL = 'megamat'
ARM = 'A1'
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
CALYX_STAGES: Tuple[str, ...] = ('kc_claw', 'kc_soma')
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
# the views read (hal's VIEWS) and the PN dendrites' fourth (E134's lock)
VIEWS: Tuple[str, ...] = ('strength', 'sharing', 'breadth')
MATRIX = 'matrix'


@dataclass(frozen=True)
class Cell:
    """One fit of all megamat stages whose evaluations are read.

    Attributes:
        name: 'A1' (E136's), 'X', 'F' or 'W' (E145's).
        source: the experiment whose fit it is.
        searched: 'S' or 'T', the wirings it searched (its evaluations' wirings).
        fit_half: 'a' (E136's fit half) or 'b' (E136's held-out half), the half it was fitted to.
        n_evals: the evaluations of its fit history.
    """
    name: str
    source: str
    searched: str
    fit_half: str
    n_evals: int

    @property
    def wirings(self) -> Tuple[str, ...]:
        return S_WIRINGS if self.searched == 'S' else T_WIRINGS

    @property
    def other_half(self) -> str:
        return 'b' if self.fit_half == 'a' else 'a'


CELLS: Dict[str, Cell] = {c.name: c for c in (
    Cell('A1', E136, 'S', 'a', 124),
    Cell('X', E145, 'T', 'b', 139),
    Cell('F', E145, 'S', 'b', 126),
    Cell('W', E145, 'T', 'a', 139),
)}
# the saved settings' pseudo-cell of `--part views` (E143's settings on T, the family on all eight wirings, E145's
# classed points), and the cells `--part views` takes
SETTINGS_CELL = 'settings'
VIEW_CELLS: Tuple[str, ...] = (*CELLS, SETTINGS_CELL)
# the two pooled sets: the two fits on the same wirings
POOLS: Dict[str, Tuple[str, ...]] = {'S': ('A1', 'F'), 'T': ('X', 'W')}
# E145's cell variants read beside the settings (each fit's best evaluation on its searched wirings, its reads elsewhere)
CELL_VARIANTS: Tuple[str, ...] = ('E145_A1', 'E145_X', 'E145_F', 'E145_W')

# the losses: name -> (description, terms); a term (stage, views) is the stage's population error over those views (the
# root mean square of the views' errors; 'all' is E136's measure 'all': three views, four at the PN dendrites with the
# matrix); a loss is the root mean square of its terms
STRENGTH = ('strength',)
STRENGTH_SHARING = ('strength', 'sharing')
ALL = 'all'
LOSSES: Dict[str, Tuple[str, Tuple[Tuple[str, object], ...]]] = {
    'L0': ('strength at the four stages, the adopted measure and the fits\' loss',
           tuple((s, STRENGTH) for s in LOSS_STAGES)),
    'L1': ("strength and sharing at the four stages, E134's two-view error",
           tuple((s, STRENGTH_SHARING) for s in LOSS_STAGES)),
    'L2': ("all three views at the four stages (E136's measure 'all': the PN dendrites on four views with the matrix)",
           tuple((s, ALL) for s in LOSS_STAGES)),
    'L3': ("L0 with the claws' and somata's errors over all views, the classes' guard made a part of the loss",
           tuple((s, STRENGTH) for s in LOSS_STAGES) + tuple((s, ALL) for s in CALYX_STAGES)),
    'L4': ('L0 without the boutons', tuple((s, STRENGTH) for s in LOSS_STAGES if s != 'pn_bouton')),
    'L5': ('strength at the PN dendrites and boutons', (('pn_dendrite', STRENGTH), ('pn_bouton', STRENGTH))),
    'L6': ('strength at the claws and somata', (('kc_claw', STRENGTH), ('kc_soma', STRENGTH))),
}
# the yardsticks and thresholds (the entry's)
CONTROL_FALL_KEY = 'control_fall'          # E136's measure_strength classes.json (0.0403)
DECIDED_MEDIAN, DECIDED_MIN = 0.8, 0.5     # "is decided", "can be decided" (E143's tells-apart yardstick)
NO_READING_SPREADS = 2.0                   # "no reading": the models' range under this many of the view's spreads
TRANSFERS_MAX_REGRET = 0.25
CHECK_TOL = 0.02                           # E143's read-once check
IDENTITY_TOL = 1e-9                        # the rescore tables reproduced from the views by E136's own reading
N_JOINT, JOINT_SEED = 200, 0               # E143's joint splits
PROBE_TOL = 0.005                          # the planner's probe reproduced (its table rounds to 2 decimals)

# the planner's probe as the entry gives it (rounded as printed), reproduced in step 0
PROBE: Dict[str, Dict[str, object]] = {
    'consistency_L0': dict(A1=0.34, X=0.79, F=0.80, W=0.34),
    'consistency_L1': dict(A1=0.82, X=0.95, F=0.94, W=0.95),
    'consistency_L2': dict(A1=0.70, X=0.95, F=0.72, W=0.94),
    'soma_strength_sharing_breadth': dict(A1=(0.28, 0.97, 0.95), X=(0.50, 0.99, 0.95), F=(0.79, 0.97, 0.91),
                                          W=(-0.04, 0.99, 0.98)),
    'claw_strength_sharing_breadth': dict(A1=(-0.02, 0.54, -0.37), X=(0.42, 0.33, -0.09), F=(0.32, 0.08, -0.17),
                                          W=(-0.32, 0.48, -0.49)),
    'soma_range_strength_sharing': dict(A1=(0.17, 0.53), X=(0.29, 0.47), F=(0.36, 0.60), W=(0.18, 0.57)),
    'holds_part_half_a': dict(A1=(5, 111), X=(1, 116), F=(0, 106), W=(9, 123)),
    'holds_part_half_b': dict(A1=(102, 112), X=(121, 127), F=(102, 112), W=(19, 126)),
    'pick_L0_factors': dict(A1=(0.155, 3.24, 32), X=(0.175, 7.3, 68), F=(0.43, 2.53, 52), W=(0.082, 2.65, 35)),
    'pick_L1_factors': dict(A1=(0.23, 3.4, 43), X=(0.73, 1.6, 72), F=(0.61, 2.4, 72), W=(0.32, 2.4, 42)),
    'regret_L0_L1': dict(A1=(0.15, 0.04), X=(0.12, 0.04), F=(0.07, 0.00), W=(0.20, 0.11)),
    'pooled_L0': dict(S=0.52, T=-0.02),
    'pooled_L1': dict(S=0.97, T=0.97),
    'pooled_L2': dict(S=0.91, T=0.97),
    'claws_worse_among_lowering': dict(a=(429, 456), b=(7, 477)),
    'w_somata_all_views_fail_half_b': (101, 126),
}
# E145's consistency tables and step 0's (the identities): X, F, W from E145's consistency_strength.csv; A1 over its
# 120 evaluations with both rates held (E145's step 0) and over its 124 (E136's `consistency` on A1's rescore alone)
CONSISTENCY_IDENTITIES: Dict[str, float] = {'X': 0.786, 'F': 0.802, 'W': 0.341, 'A1_held': 0.338, 'A1_all': 0.357}

CLASSES: Dict[str, str] = {
    'decided': f'a view "is decided" at a stage, over a set of models: its consistency is {DECIDED_MEDIAN:g} or more at '
               f'the median split and {DECIDED_MIN:g} or more on every split (E143\'s yardstick); it has "no reading" if '
               f'the models\' errors range over less than {NO_READING_SPREADS:g} times the view\'s spread over the '
               "family's eight wirings",
    'reads_calyx': '"megamat reads the calyx beside strength": sharing or breadth is decided at the claws or the somata '
                   'over every cell\'s evaluations, over both pooled sets and over the nine settings at pn_inhadd 1000 '
                   '(the same stage and view on every set)',
    'loss_decided': f'a loss "can be decided": {DECIDED_MEDIAN:g} or more at the median joint split and {DECIDED_MIN:g} '
                    'or more on every joint split, in all four cells and on both pooled sets',
    'transfers': f'a pick "transfers": its regret is {TRANSFERS_MAX_REGRET:g} or less',
    'claws_read_half': '"the claws\' cost is the read half\'s": in all four cells the classed point\'s claws on strength '
                       'are worse than the spread against one half of E136\'s split and not against the other, the same '
                       'half in each; "one fly\'s" if one claw fly alone gives that sign (in all four cells, and no '
                       'other claw fly alone does), "two against two" otherwise',
    'no_boutons_orders': '"without the boutons the megamat stage error orders the settings": E143\'s yardstick on the '
                         '200 joint splits, over E143\'s settings without E135\'s four arms',
    'no_class': 'no class for a cell whose identities fail',
}
OUTCOMES: Dict[str, str] = {
    'a': 'Megamat reads the calyx beside strength, a loss with sharing can be decided, and its picks transfer (about 2 '
         'in 5): to Tom and Yang under both open decisions, with where that loss\'s picks put the calyx.',
    'b': 'Decided over the fits\' evaluations, but not among the settings that differ only in the calyx, or not on '
         'every split (about 2 in 5): E143\'s reading stands for the other views too.',
    'c': 'No view at the claws or the somata is decided on every split (about 1 in 6): megamat\'s calyx stages decide '
         'nothing under any reading.',
    'd': 'The claws\' cost is the read half\'s (about 9 in 10, beside the others): the classes stand as scored; how the '
         'claws enter a class goes to Tom and Yang.',
    'e': 'Without the boutons the megamat stage error still orders the settings (about 4 in 5, beside the others).',
    'f': 'An identity or a read-once check fails in a cell (about 1 in 5, beside the others): that cell is read by '
         'E136\'s own reading on fewer splits.',
}
PREDICTION: Dict[str, str] = {
    'reads_calyx_over_every_fit_on_all_three_splits': 'about 3 in 4',
    'reads_calyx_over_the_nine_settings_too': 'about 1 in 2',
    'a_loss_with_sharing_can_be_decided': 'about 3 in 5',
    'claws_cost_is_the_read_half': 'about 9 in 10',
    'claws_cost_is_one_fly': 'about 1 in 3',
    'without_boutons_orders_the_settings': 'about 4 in 5',
    'source': "the planner's chances (the entry's), written before any view is read",
}
# the judgments the entry leaves open (the narrowest reading, closest to E136's, E143's and E145's code)
JUDGMENTS: Tuple[str, ...] = (
    "An evaluation's value at a stage is the error of its four searched wirings' mean views (hal's wirings_as_flies, "
    "row 'common'), as the fits' loss and E136's rescore read it; the planner's probe reproduces on that reading "
    "(checked in step 0). Each wiring's value is kept beside. The settings are read as E143 reads them: the mean over T "
    "of each wiring's error.",
    "The probe counted every evaluation (A1's 124, F's 126); the lock's definitions leave out an evaluation whose rates "
    "were not held on every searched wiring and count it (A1: 4, F: 2), so step 0 reproduces the probe on every "
    "evaluation and every later part reads the held ones.",
    "L2 and 'over all views' at the PN dendrites read four views with the matrix (E134's lock, E136's measure 'all', "
    "error_with_matrix); the claws and somata have three. The probe's L2 consistencies reproduce on that reading.",
    "A loss's term is a stage's population error over the term's views (the root mean square of the views' errors); a "
    "loss is the root mean square of its terms, so L0, L1 and L2 are E136's measures 'strength', 'strength_sharing' "
    "and 'all'.",
    "A pick is among the held evaluations; its regret's start is evaluation 0 (the family at z = 0 on the searched "
    "wirings) and the other half's lowest is over the held evaluations.",
    "Megamat's part of \"holds\" reads changes from evaluation 0 on the 'common' reading against the half, with "
    "E136's spreads as scored (classes.json 'spreads' for strength, beside.spreads_all_views for all views) and the "
    "control fit's fall (classes.json control_fall) as the bar, as the probe reads it (reproduced in step 0).",
    "The read-once check in a cell is E143's: on E136's split, against each half, the read-once megamat stage error on "
    "strength at the cell's evaluation 0 and at its pick (the fit's best evaluation) within CHECK_TOL of E136's own "
    "reading of the same runs; the largest difference over every evaluation, per stage and view, is given beside.",
    "'Every split' of a view's consistency at a stage is that stage's splits (E143's stage_splits: 35, 10, 3, 3); a "
    "loss's consistency reads E143's 200 joint splits, which decide, with E136's split beside (as E143's tells-apart).",
    "A view's spread is the SD over the family's eight wirings of its per-wiring error on that view against all the "
    "stage's flies (read once); a set's range is the largest minus the smallest of its models' errors on that view "
    "against all the stage's flies.",
    "\"Megamat reads the calyx beside strength\" needs one stage and view (claws or somata, sharing or breadth) decided "
    "on every set it names: each of the four cells' held evaluations, both pooled sets and the nine settings at "
    "pn_inhadd 1000.",
    "The classed points' claws (design 6, the claws' class) are E136's per-wiring comparison (`compare`) of the cell's "
    "classed point with the family on the cell's four read wirings, against each half of E136's split, with E136's "
    "kc_claw spread on strength; a single claw fly as the reference likewise. \"One fly's\": one claw fly alone gives "
    "the sign in all four cells and no other claw fly alone does.",
    "Without the boutons (L4) the settings are E143's settings without E135's four arms (E135 set their PN stages "
    "from the dendrite flies); with them beside.",
    "No view's reading is weighted by its reliability; no entry is dropped (validation2 and E138's readings are not "
    "used).",
)
HOLDS: Tuple[str, ...] = (
    'no simulation, no model run, no olfsysm change, nothing fitted, searched or chosen',
    'no measure, class or yardstick changed: strength by odor alone still chooses and gates',
    'no held-out panel enters a fit; validation2\'s saved scores are read beside only, nothing chosen from E138\'s '
    'readings of it',
    "Pratyush's hemibrain only; Remy's processed data without the bouton fly 2022-10-21 fly02 and her older megamat "
    'flies (E134\'s lock)',
)


# ---------------------------------------------------------------------------
# paths and small helpers

DIRS: Dict[str, Path] = {}


def results_dir() -> Path:
    """<shared>/results/E148_megamat_other_views/ (created); DIRS['results'] in tests."""
    d = Path(DIRS['results']) if 'results' in DIRS else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def views_dir() -> Path:
    """results/E148_megamat_other_views/views/ (created): each cell's views and the flies' parts."""
    d = results_dir() / 'views'
    d.mkdir(parents=True, exist_ok=True)
    return d


def views_path(cell: str) -> Path:
    """A cell's views (`part_views`), one parquet file."""
    return views_dir() / f'{cell}.parquet'


def local_dir(part: str = '') -> Path:
    """outputs/E148_megamat_other_views/[<part>/] (created)."""
    root = Path(DIRS['local']) if 'local' in DIRS else H.outputs_dir(EXPERIMENT)
    d = root / part if part else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def lock_path() -> Path:
    """The lock's shared copy; the committed copy is LOCK_COPY."""
    return results_dir() / LOCK_NAME


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def get_cell(name: str) -> Cell:
    if name not in CELLS:
        raise ValueError(f'cell {name!r} not in {sorted(CELLS)}')
    return CELLS[name]


# ---------------------------------------------------------------------------
# the lock

def lock_files() -> Dict[str, Path]:
    """The files the lock records by sha256 (name -> path): E134's lock, E136's and E145's fit records, histories,
    rescore tables, consistency tables and classes, E143's lock and consistency tables."""
    raise NotImplementedError


def lock_contents() -> dict:
    """The lock's contents: the cells, pools, settings and their groups, views, readings, splits, losses, pick and
    regret, megamat's part of "holds", the classes' thresholds, the outcomes, the prediction, the judgments, the holds
    and the files with their sha256."""
    raise NotImplementedError


def part_lock() -> dict:
    """Writes the lock to results/.../lock.json and LOCK_COPY (commit the copy before any other part)."""
    raise NotImplementedError


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (LOCK_COPY committed with no uncommitted change, and the same bytes as the shared copy), its
    files checked; PermissionError otherwise."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design 2: the views, once per cell

def cell_evals(cell: str) -> pd.DataFrame:
    """A cell's evaluations from its fit history: eval, its four factors (aplpn, kckc, inhadd as factors and pn_inhadd
    absolute, inhsc gain), whether both rates held on every searched wiring ('held'), and its loss_sim."""
    raise NotImplementedError


def cell_runs(cell: str, i: int) -> Optional[Dict[str, object]]:
    """One evaluation's saved runs on the cell's searched wirings (E136's fit folders for A1: outputs/.../measure_strength
    /fit_A1/ then the archive; E145's for X, F, W: outputs/.../<cell>/fit_A1/ then the archive's <cell>/fit_A1/); None
    if any wiring is missing."""
    raise NotImplementedError


def readings_of_stage(stage: str) -> Dict[str, Tuple[int, ...]]:
    """The readings saved at a stage: 'all' (every fly) and, as 'sub:<indices>', E136's own reading by each half of
    E136's split and, at the claws and somata, by every pair of their four flies (their three splits)."""
    raise NotImplementedError


def read_views(runs: Mapping[str, object], name: str) -> pd.DataFrame:
    """One model's runs (per wiring) read at the four megamat stages under every reading of `readings_of_stage`:
    long rows (name, wiring, reading, stage, part, entry, value), the model's parts as E143's `views_rows` reads them
    against all the stage's flies' reference."""
    raise NotImplementedError


def part_views(cell: str, evals: Optional[Sequence[int]] = None) -> dict:
    """Design 2: every evaluation's views (or, for 'settings', every saved setting's) written to views/<cell>.parquet,
    resumable (evaluations already written are kept); the evaluations it could not find are listed in
    views/<cell>_summary.json and printed. Returns the summary."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 0 (design 1)

def probe_table() -> pd.DataFrame:
    """The planner's probe reproduced from E136's and E145's rescore tables and fit histories (every evaluation, as
    the planner read them), beside PROBE's values."""
    raise NotImplementedError


def prediction() -> dict:
    """The prediction per stage and view, from the probe's consistencies on E136's split, with PREDICTION; written once,
    before any view is read."""
    raise NotImplementedError


def identities(cells: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The identities on every cell whose views exist: the rescore tables reproduced from the views by E136's own
    reading (IDENTITY_TOL), the consistencies (CONSISTENCY_IDENTITIES), E143's consistency table on strength (from the
    settings' views), and the read-once check (CHECK_TOL)."""
    raise NotImplementedError


def part_step0() -> dict:
    """Step 0: the lock checked, the probe (step0_probe.csv), the prediction (prediction.json, once) and the identities
    of the cells whose views exist (identities.csv), step0.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# designs 3 to 8 and the score

def part_room() -> dict:
    """Design 3: per stage and view, the family, A1 and each classed point against a half, beside half to half, on every
    split."""
    raise NotImplementedError


def part_views_decided() -> dict:
    """Design 4: consistency per stage and view over each cell's evaluations, the two pooled sets, the saved settings and
    their two groups, on every split; the sharing errors split into level, spread and order."""
    raise NotImplementedError


def part_losses() -> dict:
    """Design 5: L0 to L6 on E136's split and the 200 joint splits: consistency, pick, regret, megamat's part of "holds"
    at the pick against the other half, and its four factors beside the strength pick's."""
    raise NotImplementedError


def part_claws() -> dict:
    """Design 6: the family, the four classed points and the evaluations at the claws against each half, each split and
    each claw fly; the same at the somata; the odors on which the halves of the claw flies differ, as leads."""
    raise NotImplementedError


def part_no_boutons() -> dict:
    """Design 7: the megamat stage error without the boutons over the saved settings on the 200 joint splits, with and
    without E135's four arms; E141's and E145's classes recomputed with the boutons left out."""
    raise NotImplementedError


def part_beside() -> dict:
    """Design 8, classing nothing: megamat's sharing error at the somata beside the held-out panels' strength errors and
    E139's price; each stage's views regressed on the four factors per cell; E145's classes with E147's spreads; E147's
    saved populations' sharing and breadth."""
    raise NotImplementedError


def part_score() -> dict:
    """The classes and outcome letters (classes.json)."""
    raise NotImplementedError


PARTS: Tuple[str, ...] = ('lock', 'views', 'step0', 'room', 'views_decided', 'losses', 'claws', 'no_boutons', 'beside',
                          'score')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--cell', choices=VIEW_CELLS, help='--part views: the cell whose views to read')
    ap.add_argument('--evals', nargs='+', type=int, help='--part views: only these evaluations')
    a = ap.parse_args(argv)
    if a.part == 'lock':
        rec = part_lock()
        print(json.dumps({k: rec[k] for k in ('experiment', 'time', 'host')}, default=str))
    elif a.part == 'views':
        if not a.cell:
            raise SystemExit('--part views takes --cell')
        print(json.dumps(part_views(a.cell, a.evals), indent=1, default=str))
    else:
        fn = globals()[f'part_{a.part}']
        print(json.dumps(fn(), indent=1, default=str)[:20000])


if __name__ == '__main__':
    main()
