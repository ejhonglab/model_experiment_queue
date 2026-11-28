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
import contextlib
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
CONTROL_FALL_KEY = 'control_fall'          # E136's measure_strength classes.json (0.0403), beside
FALL_BAR = 0.040                           # megamat's part of "holds": "the control fit's 0.040", as the entry and the
                                           # probe read it (E147's CONTROL_FALL)
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
    "E136's spreads as scored (classes.json 'spreads' for strength, beside.spreads_all_views for all views) and 0.040 "
    "as the bar (the entry's \"the control fit's 0.040\", E147's CONTROL_FALL; the control fit's fall is 0.0403, beside), "
    "as the probe reads it: with 0.0403 A1's evaluation 9 (a fall of 0.04026 against half b) would not count. The probe's "
    "429 claws worse of 456 against half a reproduces with the claws' spread rounded to 0.017 (W's evaluation 113 at "
    "a change of 0.01686); the lock reads the spread as scored (0.01669), 430.",
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
# revisions of the lock after it was written, appended (each with its time and why); none changes a definition above
REVISIONS: Tuple[dict, ...] = (dict(
    time='2026-11-27 20:46',
    what=("A cell enters designs 3 to 8 once its views hold every evaluation of its fit history; a pool once both cells do. "
    "A set read by E136's own reading (a failed read-once check) is read on the subsets its views saved: E136's split at "
    "the PN dendrites and boutons and the claws' and somata's three splits; its joint splits are those nine.",
    "Outcome (c) reads every set and every view at the claws and somata (strength included); outcome (b) holds where a "
    "calyx stage's sharing or breadth is decided, or at 0.8 at the median split, over every fit's held evaluations and "
    "both pools but not over the nine settings at pn_inhadd 1000."),
    why=("two readings the first lock left open, written when designs 3 to 8 were coded, before any of their statistics "
         "was read; step 0 had read only the probe and the read-once check of F's evaluations 0, 1 and 116")),)
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

def _e136_dir() -> Path:
    """E136's results on the adopted measure (its A1 cell's records)."""
    return H.results_dir(E136) / 'measure_strength'


def cell_results(cell: str) -> Path:
    """A cell's fit records (fit history, fit record, rescore and consistency tables): E136's measure_strength/ for A1,
    E145's <cell>/ for the others."""
    get_cell(cell)
    return _e136_dir() if cell == 'A1' else H.results_dir(E145) / cell


def lock_files() -> Dict[str, Path]:
    """The files the lock records by sha256 (name -> path): E134's lock, each cell's fit record, fit history, rescore
    and consistency tables, E136's classes, E145's classes and step 0, and E143's lock and consistency tables."""
    r = H.shared_dir() / 'results'
    files = {'e134_lock': H.HERE / 'experiments' / 'E134_lock.json',
             'e136_classes': _e136_dir() / 'classes.json',
             'e136_consistency_strength': _e136_dir() / 'consistency_strength.csv',
             'e145_classes': r / E145 / 'classes.json',
             'e145_step0': r / E145 / 'step0.json',
             'e143_lock': r / E143 / 'lock.json',
             'e143_splits_consistency': r / E143 / 'splits_consistency.csv',
             'e143_splits_joint_consistency': r / E143 / 'splits_joint_consistency.csv',
             'e143_calyx_consistency': r / E143 / 'calyx_consistency.csv',
             'e143_classes': r / E143 / 'classes.json',
             'e139_classes': r / E139 / 'classes.json',
             'e141_classes': r / E141 / 'classes.json'}
    for c in CELLS:
        d = cell_results(c)
        files.update({f'{c}_fit': d / f'fit_{ARM}.json', f'{c}_history': d / f'fit_history_{ARM}.csv',
                      f'{c}_rescore': d / f'rescore_{ARM}.csv'})
        if c != 'A1':
            files.update({f'{c}_consistency_{m}': d / f'consistency_{m}.csv'
                          for m in ('strength', 'strength_sharing', 'all')})
    return files


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def file_record(path: Path) -> dict:
    """{'path' ('<shared>/...' or '<repo>/...'), 'bytes', 'sha256'} of a file."""
    p = Path(path).resolve()
    rel = str(p)
    for root, tag in ((H.shared_dir().resolve(), '<shared>'), (H.HERE.resolve(), '<repo>')):
        try:
            rel = f'{tag}/{p.relative_to(root)}'
            break
        except ValueError:
            continue
    return dict(path=rel, bytes=p.stat().st_size, sha256=_sha256(p))


def resolve(path: str) -> Path:
    """A path the lock records."""
    if path.startswith('<shared>/'):
        return H.shared_dir() / path[len('<shared>/'):]
    if path.startswith('<repo>/'):
        return H.HERE / path[len('<repo>/'):]
    return Path(path)


def check_files(lock: Mapping) -> None:
    """Every file of the lock's 'files' has its recorded sha256; ValueError otherwise."""
    for name, e in lock['files'].items():
        got = _sha256(resolve(e['path']))
        if got != e['sha256']:
            raise ValueError(f"{name} ({e['path']}): sha256 {got[:12]} differs from the lock's {e['sha256'][:12]}")


def _jsonable(x):
    return json.loads(json.dumps(x, default=lambda v: v.item() if hasattr(v, 'item') else str(v)))


def lock_contents() -> dict:
    """The lock's contents: the cells, pools, settings and their groups, views, readings, splits, losses, pick and
    regret, megamat's part of "holds", the classes' thresholds, the outcomes, the prediction, the judgments, the holds
    and the files with their sha256."""
    from experiments import E143_megamat_fly_sampling as E143m
    cl = json.loads((_e136_dir() / 'classes.json').read_text())
    settings = [dict(name=s.name, source=s.source, arm=s.arm, label=s.label, role=s.role, calyx=s.calyx)
                for s in E143m.SETTINGS]
    return _jsonable(dict(
        experiment=EXPERIMENT, step='lock', script=SCRIPT, planner=PLANNER, host=_host(), time=_now(),
        commits=H.code_commits(),
        entry='experiment_priorities.md item E148 (proposed and moved by viki 2026-11-27 20:16)',
        cells={k: dict(source=c.source, searched=c.searched, wirings=list(c.wirings), fit_half=c.fit_half,
                       other_half=c.other_half, n_evals=c.n_evals) for k, c in CELLS.items()},
        halves=dict(a="E136's fit half: the first part of hal's split_halves(n, 1, seed 0) of each stage's flies",
                    b="E136's held-out half: the rest"),
        pools={k: list(v) for k, v in POOLS.items()},
        evaluations=dict(value="the error of the searched wirings' mean views (row 'common'), each wiring's beside",
                         held="both rates held on every searched wiring (the fit history's hold_ok_<wiring>); an "
                              "evaluation not held is left out and counted", start='evaluation 0 (the family at z = 0)'),
        settings=dict(list=settings, read_on=list(T_WIRINGS), value="the mean over T of each wiring's error (E143)",
                      groups={g: list(v) for g, v in E143m.CALYX_GROUPS.items()},
                      without_e135=[s.name for s in E143m.SETTINGS if s.role == 'setting' and s.source != 'E135'],
                      beside=dict(family="the unrefitted family on all eight wirings (E134's held megamat runs)",
                                  cell_variants=list(CELL_VARIANTS),
                                  cell_variants_rule="each fit's best evaluation on its searched wirings, its reads "
                                                     "elsewhere (E145_score.cell_loader)")),
        views=dict(views=list(VIEWS), matrix="the PN dendrites' fourth view (E134's lock)",
                   parts="hal's kcq/population_fit.py ALL_PARTS"),
        reading=dict(default="E134's default readings (kcq.stage_readings DEFAULT_READINGS), the curve, noise and "
                             "criterion of all the stage's flies, seed index the wiring's position in E136's WIRINGS "
                             "(E143's read_runs); a subset of flies changes only the reference, in all the stage's "
                             "flies' one-fly units",
                     beside="E136's own reading (the half's flies set the curve, noise and criterion) on E136's split "
                            "at every stage and on every pair of the claws' and somata's four flies",
                     check_tol=CHECK_TOL, identity_tol=IDENTITY_TOL),
        splits=dict(stage="E143's stage_splits (35, 10, 3, 3)", joint=dict(n=N_JOINT, seed=JOINT_SEED,
                    rule="E143's joint_splits, E136's own split beside"), single_fly=list(CALYX_STAGES)),
        losses={k: dict(description=d, terms=[[s, v if isinstance(v, str) else list(v)] for s, v in t])
                for k, (d, t) in LOSSES.items()},
        pick='the held evaluation with the lowest loss against the fit half',
        regret="against the other half under the same loss: (the pick's loss - that half's lowest) / (evaluation 0's "
               "loss - that half's lowest), over the held evaluations",
        holds_part=dict(rule="a fall of the megamat stage error on strength from evaluation 0 by more than the control "
                             "fit's fall, no stage's strength error up by more than its spread, the claws' and somata's "
                             "errors over all views not up by more than theirs; on the 'common' reading against the half",
                        bar=FALL_BAR, control_fall_beside=float(cl[CONTROL_FALL_KEY]), spreads=cl['spreads'],
                        spreads_all_views=cl['beside']['spreads_all_views']),
        thresholds=dict(decided_median=DECIDED_MEDIAN, decided_min=DECIDED_MIN, no_reading_spreads=NO_READING_SPREADS,
                        transfers_max_regret=TRANSFERS_MAX_REGRET),
        probe=PROBE, consistency_identities=CONSISTENCY_IDENTITIES, probe_tol=PROBE_TOL,
        classes=CLASSES, outcomes=OUTCOMES, prediction=PREDICTION, judgments=list(JUDGMENTS), holds=list(HOLDS),
        revisions=list(REVISIONS),
        files={k: file_record(v) for k, v in lock_files().items()}))


def part_lock() -> dict:
    """Writes the lock to results/.../lock.json and LOCK_COPY (commit the copy before any other part); refuses once any
    view is written, unless the lock is unchanged."""
    lock = lock_contents()
    if LOCK_COPY.exists() and any(views_dir().glob('*.parquet')):
        old = json.loads(LOCK_COPY.read_text())
        strip = lambda d, also=(): {k: v for k, v in d.items() if k not in ('time', 'host', 'commits') + tuple(also)}
        if strip(old) == strip(lock):
            return old
        r_old, r_new = list(old.get('revisions') or []), list(lock.get('revisions') or [])
        if not (strip(old, ('revisions',)) == strip(lock, ('revisions',)) and len(r_new) > len(r_old)
                and r_new[:len(r_old)] == r_old):
            raise RuntimeError('views exist and the lock would change: append a revision (REVISIONS) instead')
        lock = dict(lock, time=old.get('time'), host=old.get('host'), commits=old.get('commits'),
                    revised=_now(), revised_host=_host(), revised_commits=H.code_commits())
    txt = json.dumps(lock, indent=1) + '\n'
    LOCK_COPY.write_text(txt)
    lock_path().write_text(txt)
    print(f'wrote {LOCK_COPY} and {lock_path()}; commit {LOCK_COPY.name} before any other part', flush=True)
    return lock


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (LOCK_COPY committed with no uncommitted change, and the same bytes as the shared copy), its
    files checked; PermissionError otherwise."""
    from kcq import stage_readings as R
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it)')
    if not R._committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = lock_path() if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    check_files(lock)
    lock['sha256'] = _sha256(copy)
    return lock


def provenance(part: str, lock: Mapping, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time, the lock's sha256 and `extra`."""
    _write_json(results_dir() / f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
                     lock_sha256=lock.get('sha256'), **extra))


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(_jsonable(obj), indent=1) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    path = results_dir() / name
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)
    return path


# ---------------------------------------------------------------------------
# design 2: the views, once per cell

def cell_evals(cell: str) -> pd.DataFrame:
    """A cell's evaluations from its fit history: eval, its four factors (aplpn, kckc, inhadd as factors and pn_inhadd
    absolute, inhsc gain), whether both rates held on every searched wiring ('held'), and its loss_sim."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    c = get_cell(cell)
    h = pd.read_csv(cell_results(cell) / f'fit_history_{ARM}.csv', float_precision='round_trip')
    out = pd.DataFrame({'eval': h['eval'].astype(int)})
    for k in ('aplpn', 'kckc', 'inhadd', 'inhsc'):
        out[k] = np.exp(h[f'z_{k}'].astype(float))
    out['pn_inhadd'] = PN_INHADD_DEFAULT * out['inhadd']
    out['held'] = h[[f'hold_ok_{w}' for w in c.wirings]].astype(bool).all(axis=1).to_numpy()
    out['loss_sim'] = h['loss_sim'].to_numpy(dtype=float)
    return out


# other checkouts' outputs/ folders where a cell's saved stages may be (`--outputs-root`; e.g. a scratch worktree reading
# the stages its host's main checkout saved), searched after this checkout's outputs and before the archive
EXTRA_OUTPUTS: List[Path] = []


def _fit_dirs(cell: str) -> List[Path]:
    """Where a cell's evaluations' stages are: this checkout's outputs, those of EXTRA_OUTPUTS, then the archive."""
    from experiments import E145_stage_fit_replication as E
    if cell == 'A1':
        from experiments import E145_score as S
        d = S._e136_dirs()['fit']
        extra = [Path(r) / E136 / 'measure_strength' / f'fit_{ARM}' for r in EXTRA_OUTPUTS]
    else:
        d = [E.local_dir(cell, f'fit_{ARM}'), E.archive_dir(cell) / f'fit_{ARM}']
        extra = [Path(r) / E145 / cell / f'fit_{ARM}' for r in EXTRA_OUTPUTS]
    return d[:1] + [x for x in extra if x.is_dir()] + d[1:]


def cell_runs(cell: str, i: int) -> Optional[Dict[str, object]]:
    """One evaluation's saved runs on the cell's searched wirings (E136's fit folders for A1: outputs/.../measure_strength
    /fit_A1/ then the archive; E145's for X, F, W: outputs/.../<cell>/fit_A1/ then the archive's <cell>/fit_A1/); None
    if any wiring is missing."""
    from kcq import stage_readings as R
    dirs = _fit_dirs(cell)
    try:
        return {w: R.load_step3_run(f'{ARM}_e{int(i)}_{w}', dirs=dirs) for w in get_cell(cell).wirings}
    except FileNotFoundError:
        return None


def setting_runs() -> Dict[str, Dict[str, object]]:
    """The saved settings' runs (the 'settings' cell): E143's settings on T (E143's `setting_runs`), the unrefitted
    family on all eight wirings, and E145's cells' variants on all eight (E145_score's `cell_loader`); per name, its
    runs per wiring (missing wirings left out)."""
    from experiments import E136_stage_population_fit as F
    from experiments import E143_megamat_fly_sampling as E143m
    from experiments import E145_score as S
    out: Dict[str, Dict[str, object]] = {}
    for st in E143m.SETTINGS:
        out[st.name] = E143m.setting_runs(st, PANEL, T_WIRINGS)
    for w in S_WIRINGS:
        r = F.load_read(F.ReadPoint(F.FAMILY, PANEL, w, 'held'))
        if r is not None:
            out['family'][w] = r
    for v in CELL_VARIANTS:
        load = S.cell_loader(v.split('_', 1)[1])
        out[v] = {w: r for w in WIRINGS for r in [load(F.ReadPoint(v, PANEL, w, None))] if r is not None}
    return out


def sub_name(idx: Sequence[int]) -> str:
    """A reading by a subset of a stage's flies, by their indices in the stage's population: 'sub:0-2'."""
    return 'sub:' + '-'.join(str(int(i)) for i in sorted(idx))


def sub_indices(reading: str, n: int) -> Tuple[int, ...]:
    """The flies' indices of a reading ('all': every fly)."""
    if reading == 'all':
        return tuple(range(n))
    return tuple(int(x) for x in reading.split(':', 1)[1].split('-'))


_FLIES: Dict[str, object] = {}


def stage_flies(stage: str):
    """viki's recorded megamat population of a stage (all flies), cached."""
    from kcq import stage_readings as R
    if stage not in _FLIES:
        fl = R.fly_side(stage, PANEL, 'on')
        if fl is None:
            raise RuntimeError(f'no recorded {stage} on {PANEL}')
        _FLIES[stage] = fl
    return _FLIES[stage]


def readings_of_stage(stage: str, n: Optional[int] = None) -> Dict[str, Tuple[int, ...]]:
    """The readings saved at a stage: 'all' (every fly) and, as 'sub:<indices>', E136's own reading by each half of
    E136's split and, at the claws and somata, by every pair of their four flies (their three splits)."""
    from experiments import E143_megamat_fly_sampling as E143m
    n = len(stage_flies(stage).flies) if n is None else int(n)
    out = {'all': tuple(range(n))}
    for h in E143m.e136_split(n):
        out[sub_name(h)] = tuple(sorted(h))
    if stage in CALYX_STAGES:
        for a, b in E143m.stage_splits(n):
            out[sub_name(a)] = tuple(sorted(a))
            out[sub_name(b)] = tuple(sorted(b))
    return out


def reading_groups() -> List[Dict[str, str]]:
    """The reads of one run: per read, the reading of each stage it sets ('all' at every stage first, then the subsets,
    the j-th subset of every stage that has one in the j-th read)."""
    subs = {st: [r for r in readings_of_stage(st) if r != 'all'] for st in LOSS_STAGES}
    groups = [{st: 'all' for st in LOSS_STAGES}]
    for j in range(max(len(v) for v in subs.values())):
        groups.append({st: v[j] for st, v in subs.items() if j < len(v)})
    return groups


def _cache_of(group: Mapping[str, str]):
    """The FlyCache of one read: per stage of `group`, its reading's flies (megamat, APL on), nothing else."""
    from kcq import stage_readings as R
    fixed = {}
    for st, rd in group.items():
        fl = stage_flies(st)
        names = [str(f) for f in fl.flies]
        keep = [names[i] for i in sub_indices(rd, len(names))]
        fixed[(st, PANEL, 'on')] = fl if rd == 'all' else R.subset_flies(fl, keep)
        fixed[(st, PANEL, 'off')] = None
    return R.FlyCache(fixed=fixed, fallback=False)


def read_views(runs: Mapping[str, object], name: str, lock: Optional[Mapping] = None) -> pd.DataFrame:
    """One model's runs (per wiring) read at the four megamat stages under every reading of `readings_of_stage`:
    long rows (name, wiring, reading, stage, part, entry, value), the model's parts as E143's `views_rows` reads them
    against all the stage's flies' reference (E143's `_reference`, inside `views_at` on views_dir: the flies' parts
    saved there once)."""
    from experiments import E136_stage_population_fit as F
    from experiments import E143_megamat_fly_sampling as E143m
    from kcq import population_fit as P
    from kcq import stage_readings as R
    lock = lock if lock is not None else R.require_lock()
    groups = reading_groups()
    caches = [_cache_of(g) for g in groups]
    rows = []
    with _in_views():
        for w, run in runs.items():
            for g, cache in zip(groups, caches):
                rd = R.read_run(run, stages=list(g), lock=lock, flies=cache, seed_index=F.seed_index(w))
                for st, reading in g.items():
                    key = (st, PANEL, 'on')
                    if key not in rd.populations:
                        raise RuntimeError(f'{name} {w}: {st} not read ({rd.skipped})')
                    pop = rd.populations[key]
                    ref = E143m._reference('/'.join(key), list(pop.values.columns))
                    for q, x in P.reference_views(pop, ref).items():
                        if q in ref.parts:
                            rows.append(pd.DataFrame(dict(name=name, wiring=w, reading=reading, stage=st, part=q,
                                                          entry=list(x.index), value=x.to_numpy(dtype=float))))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


@contextlib.contextmanager
def _in_views():
    """E143's `views_at` on views_dir(), unless already inside it (so that its caches of the flies' references are kept
    over a cell's evaluations)."""
    from experiments import E143_megamat_fly_sampling as E143m
    if E143m.VIEWS_AT.get('dir') == views_dir():
        yield views_dir()
    else:
        with E143m.views_at(views_dir()) as d:
            yield d


VIEW_COLUMNS: Tuple[str, ...] = ('cell', 'name', 'eval', 'wiring', 'reading', 'stage', 'part', 'entry', 'value')


def _categorical(df: pd.DataFrame) -> pd.DataFrame:
    for c in ('cell', 'name', 'wiring', 'reading', 'stage', 'part', 'entry'):
        df[c] = df[c].astype(str).astype('category')
    return df


def part_views(cell: str, evals: Optional[Sequence[int]] = None) -> dict:
    """Design 2: every evaluation's views (or, for 'settings', every saved setting's) written to views/<cell>.parquet,
    resumable (each evaluation's or setting's rows are kept in outputs/.../views_chunks/<cell>/ as they are read, and
    the file is written from them); the evaluations it could not find are listed in views/<cell>_summary.json and
    printed. Returns the summary."""
    from kcq import stage_readings as R
    if cell not in VIEW_CELLS:
        raise ValueError(f'cell {cell!r} not in {VIEW_CELLS}')
    lock = require_lock()
    elock = R.require_lock()
    t0 = time.time()
    commits = {k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}   # the code that reads
    chunks = local_dir(f'views_chunks/{cell}')
    found, missing, read_here = [], [], []
    ctx = _in_views()
    ctx.__enter__()
    if cell == SETTINGS_CELL:
        todo = setting_runs()
        names = list(todo)
        for nm, runs in todo.items():
            path = chunks / f'{nm.replace(":", "__")}.parquet'
            if not runs:
                missing.append(nm)
                continue
            found.append(dict(name=nm, wirings=sorted(runs)))
            if path.exists() and sorted(set(pd.read_parquet(path, columns=['wiring'])['wiring'])) == sorted(runs):
                continue
            df = read_views(runs, nm, elock)
            df.insert(0, 'cell', cell)
            df.insert(2, 'eval', -1)
            df.to_parquet(path, index=False)
            read_here.append(nm)
            print(f'E148 views settings {nm}: {len(runs)} wirings, {time.time() - t0:.0f} s', flush=True)
    else:
        ev = cell_evals(cell)
        names = [int(i) for i in ev['eval']] if evals is None else [int(i) for i in evals]
        for i in names:
            path = chunks / f'e{i}.parquet'
            if path.exists():
                found.append(i)
                continue
            runs = cell_runs(cell, i)
            if runs is None:
                missing.append(i)
                continue
            df = read_views(runs, f'e{i}', elock)
            df.insert(0, 'cell', cell)
            df.insert(2, 'eval', i)
            df.to_parquet(path, index=False)
            found.append(i)
            read_here.append(i)
            if len(read_here) % 10 == 0:
                print(f'E148 views {cell}: {len(read_here)} read, {time.time() - t0:.0f} s', flush=True)
    ctx.__exit__(None, None, None)
    files = sorted(chunks.glob('*.parquet'))
    if files:
        df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)[list(VIEW_COLUMNS)]
        df = _categorical(df.sort_values(['eval', 'name', 'wiring', 'reading', 'stage', 'part'], kind='stable',
                                         ignore_index=True))
        out = views_path(cell)
        tmp = out.with_name(f'.{out.name}.{_host()}.{os.getpid()}.tmp')
        df.to_parquet(tmp, index=False, compression='zstd')
        os.replace(tmp, out)
    summ = dict(cell=cell, host=_host(), time=_now(), seconds=round(time.time() - t0, 1), n_wanted=len(names),
                n_found=len(found), n_read_here=len(read_here), missing=missing, complete=not missing,
                lock_sha256=lock.get('sha256'), e134_lock_sha256=elock.get('sha256'),
                readings={st: list(readings_of_stage(st)) for st in LOSS_STAGES},
                chunks_kept=len(found) - len(read_here) if cell != SETTINGS_CELL else None, **commits)
    if cell == SETTINGS_CELL:
        summ['found'] = found
    _write_json(views_dir() / f'{cell}_summary.json', summ)
    if missing:
        print(f'E148 views {cell}: {len(missing)} not found on {_host()}: {missing}', flush=True)
    return summ


def views_summary(cell: str) -> dict:
    """views/<cell>_summary.json ({} before the cell's views)."""
    return _read_json(views_dir() / f'{cell}_summary.json')


# ---------------------------------------------------------------------------
# the views as sides of E143's scoring

def flies_table(stage: str) -> pd.DataFrame:
    """The flies' parts of a megamat stage as E143's `_reference` saved them (views/flies__<stage>_megamat_on.parquet)."""
    from experiments import E143_megamat_fly_sampling as E143m
    with E143m.views_at(views_dir()):
        return pd.read_parquet(E143m.flies_path(f'{stage}/{PANEL}/on'))


_VIEWS: Dict[str, pd.DataFrame] = {}


def cell_views(cell: str) -> pd.DataFrame:
    """A cell's saved views (cached); FileNotFoundError before its `--part views`."""
    if cell not in _VIEWS:
        p = views_path(cell)
        if not p.exists():
            raise FileNotFoundError(f'{p}: run --part views --cell {cell} on the host that holds its stages')
        _VIEWS[cell] = pd.read_parquet(p)
    return _VIEWS[cell]


def have_views(cell: str) -> bool:
    return views_path(cell).exists()


def views_complete(cell: str) -> bool:
    """Whether a cell's views are written with every evaluation of its fit history (its summary's 'complete')."""
    return have_views(cell) and bool(views_summary(cell).get('complete'))


def model_rows(cell: str, stage: str, reading: str = 'all', names: Optional[Sequence[str]] = None,
               common: bool = True, each: bool = False, wirings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """A cell's model rows at one stage and reading, in E143's long layout (setting, wiring, part, entry, value):
    'setting' is '<cell>:<name>' (an evaluation's name 'e<i>'); with `common` a row per model of its wirings' mean views
    (wiring 'common', NaN-aware as hal's `mean_views`, over `wirings` or every wiring it has), with `each` every wiring's
    own rows."""
    v = cell_views(cell)
    d = v[(v['stage'] == stage) & (v['reading'] == reading)]
    if names is not None:
        d = d[d['name'].isin([str(n) for n in names])]
    if wirings is not None:
        d = d[d['wiring'].isin(list(wirings))]
    d = pd.DataFrame(dict(setting=cell + ':' + d['name'].astype(str), wiring=d['wiring'].astype(str),
                          part=d['part'].astype(str), entry=d['entry'].astype(str), value=d['value'].astype(float)))
    out = []
    if each:
        out.append(d)
    if common and len(d):
        g = d.groupby(['setting', 'part', 'entry'], sort=False, observed=True)['value']
        m = g.mean().reset_index()      # pandas' mean skips NaN, as _nanmean_rows does; all-NaN gives NaN
        m.insert(1, 'wiring', 'common')
        out.append(m[['setting', 'wiring', 'part', 'entry', 'value']])
    return pd.concat(out, ignore_index=True) if out else d


def make_side(stage: str, model: pd.DataFrame):
    """E143's Side of a megamat stage from model rows (`model_rows`) and the flies' saved parts."""
    from experiments import E143_megamat_fly_sampling as E143m
    return E143m.side_from_tables(f'{stage}/{PANEL}/on', flies_table(stage), model)


def view_errors(side, weights: np.ndarray, views) -> np.ndarray:
    """References x rows: every row's error on `views` ('all' every view, E136's error column at the PN dendrites; a
    tuple of views otherwise) against each weighted reference (E143's `side_errors`)."""
    from experiments import E143_megamat_fly_sampling as E143m
    vs = None if views == ALL else tuple(views)
    return E143m.side_errors(side, weights, vs, E143m.stage_error_col(side.key, vs))


def term_errors(stage_sides: Mapping[str, object], weights: Mapping[str, np.ndarray], terms) -> List[np.ndarray]:
    """Each term's references x rows errors (a term (stage, views); the same rows at every stage)."""
    return [view_errors(stage_sides[st], weights[st], v) for st, v in terms]


def rms(xs: Sequence[np.ndarray]) -> np.ndarray:
    return np.sqrt(np.mean(np.square(np.asarray(xs, dtype=float)), axis=0))


def single_weights(n: int, idx: Sequence[int]) -> np.ndarray:
    """1 x n weights of one subset of flies."""
    from experiments import E143_megamat_fly_sampling as E143m
    return E143m.half_weights(n, [list(idx)])


def e136_halves(stage: str) -> Dict[str, Tuple[int, ...]]:
    """E136's split of a stage's flies: {'a': its fit half's indices, 'b': its held-out half's}."""
    from experiments import E143_megamat_fly_sampling as E143m
    a, b = E143m.e136_split(len(stage_flies(stage).flies))
    return {'a': tuple(sorted(a)), 'b': tuple(sorted(b))}


# ---------------------------------------------------------------------------
# step 0 (design 1)

def rescore_table(cell: str) -> pd.DataFrame:
    """A cell's rescore table with its halves named 'a' / 'b' (E136's split; 'fit' and 'held' are the cell's own)."""
    c = get_cell(cell)
    rs = pd.read_csv(cell_results(cell) / f'rescore_{ARM}.csv', float_precision='round_trip')
    rs['half_ab'] = rs['half'].map({'fit': c.fit_half, 'held': c.other_half})
    return rs


def _rho(x, y) -> float:
    from scipy import stats
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(stats.spearmanr(x[ok], y[ok])[0]) if ok.sum() > 2 else np.nan


def _tab(rs: pd.DataFrame, measure: str, half: str) -> pd.DataFrame:
    return rs[(rs['measure'] == measure) & (rs['half_ab'] == half)].set_index('eval').sort_index()


MEASURE_OF: Dict[str, str] = {'L0': 'strength', 'L1': 'strength_sharing', 'L2': 'all'}


def holds_part_counts(rs: pd.DataFrame, half: str, evals: Sequence[int], cl: Mapping,
                      claw_spread: Optional[float] = None) -> dict:
    """Megamat's part of "holds" against a half over `evals` from a rescore table (the 'common' reading): the
    evaluations that lower the megamat stage error on strength by more than the control fit's fall, those of them with
    every condition of the rule, and those whose claws on strength are worse than their spread, and whose somata over
    all views are."""
    sp, spa, fall_bar = cl['spreads'], cl['beside']['spreads_all_views'], FALL_BAR
    s, a = _tab(rs, 'strength', half), _tab(rs, 'all', half)
    ev = [i for i in evals if i in s.index]
    fall = -(s['stage_error'] - s.loc[0, 'stage_error'])
    lower = fall > fall_bar
    worse = pd.Series(False, index=s.index)
    for st in LOSS_STAGES:
        worse |= (s[f'{st}_error'] - s.loc[0, f'{st}_error']) > sp[st]
    claw = (s['kc_claw_error'] - s.loc[0, 'kc_claw_error']) > (sp['kc_claw'] if claw_spread is None else claw_spread)
    soma_all = (a['kc_soma_error'] - a.loc[0, 'kc_soma_error']) > spa['kc_soma']
    for st in CALYX_STAGES:
        worse |= (a[f'{st}_error'] - a.loc[0, f'{st}_error']) > spa[st]
    lo = lower.loc[ev]
    return dict(n=len(ev), lowering=int(lo.sum()), holds_part=int((lower & ~worse).loc[ev].sum()),
                claws_worse_among_lowering=int((lower & claw).loc[ev].sum()),
                somata_all_views_worse_among_lowering=int((lower & soma_all).loc[ev].sum()))


def probe_cell(cell: str, evals: str = 'all') -> dict:
    """The planner's probe of one cell from its rescore table and fit history (`evals` 'all', as the planner read them,
    or 'held')."""
    c = get_cell(cell)
    rs = rescore_table(cell)
    ev = cell_evals(cell).set_index('eval')
    cl = json.loads((_e136_dir() / 'classes.json').read_text())
    use = list(ev.index) if evals == 'all' else list(ev.index[ev['held']])
    out: Dict[str, object] = dict(cell=cell, evals=evals, n=len(use))
    for L, meas in MEASURE_OF.items():
        a, b = _tab(rs, meas, 'a'), _tab(rs, meas, 'b')
        idx = [i for i in use if i in a.index and i in b.index]
        out[f'consistency_{L}'] = _rho(a.loc[idx, 'stage_error'], b.loc[idx, 'stage_error'])
    a, b = _tab(rs, 'all', 'a'), _tab(rs, 'all', 'b')
    idx = [i for i in use if i in a.index and i in b.index]
    for st, key in (('kc_soma', 'soma'), ('kc_claw', 'claw')):
        out[f'{key}_strength_sharing_breadth'] = tuple(_rho(a.loc[idx, f'{st}_view_{v}'], b.loc[idx, f'{st}_view_{v}'])
                                                       for v in VIEWS)
    f = _tab(rs, 'all', c.fit_half).loc[idx]
    out['soma_range_strength_sharing'] = tuple(float(f[f'kc_soma_view_{v}'].max() - f[f'kc_soma_view_{v}'].min())
                                               for v in ('strength', 'sharing'))
    for h in ('a', 'b'):
        hp = holds_part_counts(rs, h, idx, cl)
        out[f'holds_part_half_{h}'] = (hp['holds_part'], hp['lowering'])
        out[f'claws_worse_half_{h}'] = (hp['claws_worse_among_lowering'], hp['lowering'])
        out[f'somata_all_views_worse_half_{h}'] = (hp['somata_all_views_worse_among_lowering'], hp['lowering'])
    for L in ('L0', 'L1'):
        fit_t, oth = _tab(rs, MEASURE_OF[L], c.fit_half).loc[idx], _tab(rs, MEASURE_OF[L], c.other_half).loc[idx]
        p = int(fit_t['stage_error'].idxmin())
        out[f'pick_{L}'] = p
        out[f'pick_{L}_factors'] = (float(ev.loc[p, 'aplpn']), float(ev.loc[p, 'kckc']), float(ev.loc[p, 'inhadd']))
        lo = float(oth['stage_error'].min())
        out[f'regret_{L}'] = float((oth.loc[p, 'stage_error'] - lo) / (oth.loc[0, 'stage_error'] - lo))
    out['regret_L0_L1'] = (out.pop('regret_L0'), out.pop('regret_L1'))
    return out


def probe_pooled(evals: str = 'all') -> Dict[str, Dict[str, float]]:
    """The probe's pooled consistencies: over the two fits on the same wirings, each evaluation's megamat stage error
    against half a and against half b, under L0, L1 and L2."""
    out: Dict[str, Dict[str, float]] = {}
    for L, meas in MEASURE_OF.items():
        out[f'pooled_{L}'] = {}
        for pool, cells in POOLS.items():
            xa, xb = [], []
            for c in cells:
                rs = rescore_table(c)
                ev = cell_evals(c).set_index('eval')
                use = list(ev.index) if evals == 'all' else list(ev.index[ev['held']])
                a, b = _tab(rs, meas, 'a'), _tab(rs, meas, 'b')
                idx = [i for i in use if i in a.index and i in b.index]
                xa += list(a.loc[idx, 'stage_error'])
                xb += list(b.loc[idx, 'stage_error'])
            out[f'pooled_{L}'][pool] = _rho(xa, xb)
            out[f'pooled_{L}'][f'{pool}_n'] = len(xa)
    return out


def _close(got, want, tol: float) -> bool:
    if isinstance(want, (tuple, list)):
        return len(got) == len(want) and all(_close(g, w, tol) for g, w in zip(got, want))
    if isinstance(want, int) and not isinstance(want, bool):
        return int(got) == want
    return bool(np.isfinite(got) and abs(float(got) - float(want)) <= tol)


def probe_table() -> pd.DataFrame:
    """The planner's probe reproduced from E136's and E145's rescore tables and fit histories (every evaluation, as
    the planner read them), beside PROBE's values: one row per statistic and cell (or pool / half), with the value on
    the held evaluations beside and whether the printed value is reproduced (PROBE_TOL, counts exactly, the factors to
    their printed digits)."""
    rows = []
    allp = {c: probe_cell(c, 'all') for c in CELLS}
    heldp = {c: probe_cell(c, 'held') for c in CELLS}
    for stat, want in PROBE.items():
        if stat.startswith('pooled_') or stat in ('claws_worse_among_lowering', 'w_somata_all_views_fail_half_b'):
            continue
        for c, w in want.items():
            got, got_h = allp[c][stat], heldp[c][stat]
            if stat == 'consistency_L0' and c == 'A1':      # the probe's 0.34 is step 0's 0.338, A1's held evaluations
                got, got_h = got_h, got
            tol = PROBE_TOL
            if stat.startswith('pick_'):
                # to the printed digits, one unit of the last (A1's L1 pick prints 0.2245 as 0.23)
                ok = all(abs(g - x) <= 10 ** -_digits(x) + 1e-9 for g, x in zip(got, w))
            else:
                ok = _close(got, w, tol)
            rows.append(dict(statistic=stat, cell=c, printed=json.dumps(w), value=json.dumps(_jsonable(got)),
                             value_held=json.dumps(_jsonable(got_h)), reproduced=bool(ok)))
    pa, ph = probe_pooled('all'), probe_pooled('held')
    for stat in ('pooled_L0', 'pooled_L1', 'pooled_L2'):
        for pool, w in PROBE[stat].items():
            got = pa[stat][pool]
            rows.append(dict(statistic=stat, cell=pool, printed=json.dumps(w), value=json.dumps(got),
                             value_held=json.dumps(ph[stat][pool]), reproduced=_close(got, w, PROBE_TOL)))
    cl = json.loads((_e136_dir() / 'classes.json').read_text())
    for h in ('a', 'b'):
        got = tuple(int(sum(allp[c][f'claws_worse_half_{h}'][k] for c in CELLS)) for k in (0, 1))
        at17 = sum(holds_part_counts(rescore_table(c), h, list(cell_evals(c)['eval']), cl, 0.017)
                   ['claws_worse_among_lowering'] for c in CELLS)
        goth = tuple(int(sum(heldp[c][f'claws_worse_half_{h}'][k] for c in CELLS)) for k in (0, 1))
        w = PROBE['claws_worse_among_lowering'][h]
        rows.append(dict(statistic='claws_worse_among_lowering', cell=h, printed=json.dumps(w), value=json.dumps(got),
                         value_held=json.dumps(goth), reproduced=got == tuple(w),
                         note=f'{at17} with the claws\' spread rounded to 0.017'))
    got = allp['W']['somata_all_views_worse_half_b']
    w = PROBE['w_somata_all_views_fail_half_b']
    rows.append(dict(statistic='w_somata_all_views_fail_half_b', cell='W', printed=json.dumps(w),
                     value=json.dumps(got), value_held=json.dumps(heldp['W']['somata_all_views_worse_half_b']),
                     reproduced=tuple(got) == tuple(w)))
    return pd.DataFrame(rows)


def _digits(x: float) -> int:
    """The decimals a printed factor has (0.155 -> 3, 3.24 -> 2, 32 -> 0)."""
    s = repr(float(x)).rstrip('0').rstrip('.')
    return len(s.split('.')[1]) if '.' in s else 0


def prediction(probe: Optional[pd.DataFrame] = None) -> dict:
    """The prediction per stage and view, from the rescore tables on E136's split (over each cell's held evaluations,
    measure 'all'): a view is predicted decided at a stage over every cell if its consistency is DECIDED_MEDIAN or more
    in every cell on E136's split; with each cell's value and PREDICTION (the planner's chances); written once, before
    any view is read."""
    rows = {}
    for c in CELLS:
        rs = rescore_table(c)
        ev = cell_evals(c).set_index('eval')
        held = list(ev.index[ev['held']])
        a, b = _tab(rs, 'all', 'a'), _tab(rs, 'all', 'b')
        idx = [i for i in held if i in a.index and i in b.index]
        for st in LOSS_STAGES:
            for v in VIEWS + ((MATRIX,) if st == 'pn_dendrite' else ()):
                rows.setdefault(f'{st}/{v}', {})[c] = _rho(a.loc[idx, f'{st}_view_{v}'], b.loc[idx, f'{st}_view_{v}'])
    per = {k: dict(each=v, predicted_decided_over_every_cell=bool(all(x >= DECIDED_MEDIAN for x in v.values())))
           for k, v in rows.items()}
    return dict(per_stage_view=per, rule='predicted decided over every cell where its consistency on E136\'s split over '
                'the held evaluations (rescore, measure all) is at least DECIDED_MEDIAN in every cell', **PREDICTION)


def rescore_from_views(cell: str, evals: Sequence[int]) -> pd.DataFrame:
    """E136's rescore of a cell's evaluations recomputed from the saved views by E136's own reading (each half of
    E136's split read by its own flies, 'sub:' readings): per evaluation, measure and half, each stage's error, its
    views and each wiring's error, and the megamat stage error, in rescore's columns."""
    names = [f'e{int(i)}' for i in evals]
    c = get_cell(cell)
    rows = []
    for half in ('a', 'b'):
        per_stage = {}
        for st in LOSS_STAGES:
            rd = sub_name(e136_halves(st)[half])
            m = model_rows(cell, st, rd, names=names, common=True, each=True)
            sd = make_side(st, m)
            W = single_weights(len(sd.flies), e136_halves(st)[half])
            per_stage[st] = (sd, W)
        for meas, views in (('strength', STRENGTH), ('strength_sharing', STRENGTH_SHARING), ('all', ALL)):
            res = {}
            for st, (sd, W) in per_stage.items():
                res[st] = view_errors(sd, W, views)[0]
                if meas == 'all':
                    for v in VIEWS + ((MATRIX,) if st == 'pn_dendrite' else ()):
                        res[f'{st}_view_{v}'] = view_errors(sd, W, (v,))[0]
                else:
                    for v in views:
                        res[f'{st}_view_{v}'] = view_errors(sd, W, (v,))[0]
            rix = per_stage['pn_dendrite'][0].rows
            for nm in names:
                k = rix.get_loc((f'{cell}:{nm}', 'common'))
                r = dict(eval=int(nm[1:]), measure=meas, half_ab=half)
                for st in LOSS_STAGES:
                    r[f'{st}_error'] = float(res[st][k])
                    for col, x in res.items():
                        if col.startswith(f'{st}_view_'):
                            r[col] = float(x[k])
                    for w in c.wirings:
                        r[f'{st}_error_{w}'] = float(res[st][rix.get_loc((f'{cell}:{nm}', w))])
                r['stage_error'] = float(np.sqrt(np.mean([r[f'{st}_error'] ** 2 for st in LOSS_STAGES])))
                rows.append(r)
    return pd.DataFrame(rows)


def identity_rescore(cell: str) -> dict:
    """The cell's rescore table reproduced from its views by E136's own reading: the largest difference over every
    saved evaluation, measure, half and column, against IDENTITY_TOL."""
    rs = rescore_table(cell)
    have = sorted(int(x) for x in set(cell_views(cell)['eval']))
    got = rescore_from_views(cell, have)
    diffs = []
    for r in got.itertuples(index=False):
        ref = rs[(rs['eval'] == r.eval) & (rs['measure'] == r.measure) & (rs['half_ab'] == r.half_ab)]
        if ref.empty:
            continue
        ref = ref.iloc[0]
        for col, v in r._asdict().items():
            if col in ('eval', 'measure', 'half_ab') or col not in ref.index:
                continue
            diffs.append((abs(float(v) - float(ref[col])), r.eval, r.measure, r.half_ab, col))
    worst = max(diffs) if diffs else (np.nan, None, None, None, None)
    return dict(cell=cell, identity='rescore reproduced from the views (E136\'s own reading)', n_evals=len(have),
                n_values=len(diffs), max_abs_diff=float(worst[0]), worst=list(worst[1:]),
                passed=bool(diffs) and float(worst[0]) <= IDENTITY_TOL)


def consistency_from_views(cell: str, evals: Sequence[int], reading: str = 'own') -> float:
    """The megamat stage error on strength's consistency between the halves of E136's split over `evals`, from the
    views: by E136's own reading ('own') or read once ('all')."""
    names = [f'e{int(i)}' for i in evals]
    vals = {}
    for half in ('a', 'b'):
        terms = []
        for st in LOSS_STAGES:
            rd = sub_name(e136_halves(st)[half]) if reading == 'own' else 'all'
            sd = make_side(st, model_rows(cell, st, rd, names=names))
            terms.append(view_errors(sd, single_weights(len(sd.flies), e136_halves(st)[half]), STRENGTH)[0])
        vals[half] = rms(terms)
    return _rho(vals['a'], vals['b'])


def read_once_check(cell: str) -> dict:
    """E143's read-once check in a cell: on E136's split, against each half, the read-once megamat stage error on
    strength at evaluation 0 and at the fit's pick (its best evaluation) within CHECK_TOL of E136's own reading; the
    largest difference over every saved evaluation, per stage and view, beside."""
    fit = _read_json(cell_results(cell) / f'fit_{ARM}.json')
    have = sorted(int(x) for x in set(cell_views(cell)['eval']))
    names = [f'e{i}' for i in have]
    out = dict(cell=cell, best_eval=fit.get('best_eval'), n_evals=len(have))
    stage_vals = {}
    beside = {}
    for half in ('a', 'b'):
        own_t, once_t = [], []
        for st in LOSS_STAGES:
            W = None
            per = {}
            for rd in (sub_name(e136_halves(st)[half]), 'all'):
                sd = make_side(st, model_rows(cell, st, rd, names=names))
                W = single_weights(len(sd.flies), e136_halves(st)[half])
                per[rd] = {v: view_errors(sd, W, (v,))[0] for v in VIEWS + ((MATRIX,) if st == 'pn_dendrite' else ())}
                per[rd]['strength_err'] = view_errors(sd, W, STRENGTH)[0]
            own, once = per[sub_name(e136_halves(st)[half])], per['all']
            own_t.append(own['strength_err'])
            once_t.append(once['strength_err'])
            for v in own:
                if v != 'strength_err':
                    beside[f'{half}/{st}/{v}'] = float(np.nanmax(np.abs(own[v] - once[v])))
        stage_vals[half] = (rms(own_t), rms(once_t))
    idx = {f'e{i}': k for k, i in enumerate(have)}
    checks = []
    for half, (own, once) in stage_vals.items():
        for nm in ('e0', f"e{fit.get('best_eval')}"):
            if nm in idx:
                d = float(abs(own[idx[nm]] - once[idx[nm]]))
                checks.append(dict(half=half, eval=nm, own=float(own[idx[nm]]), once=float(once[idx[nm]]), diff=d))
        beside[f'{half}/stage_error_strength'] = float(np.nanmax(np.abs(own - once)))
    out.update(checks=checks, passed=bool(checks) and all(c['diff'] <= CHECK_TOL for c in checks),
               max_diff_every_eval=beside)
    return out


def identity_e143(lock: Mapping) -> dict:
    """E143's consistency table on strength (splits_consistency.csv: every stage's splits, over its settings and
    without E135's) reproduced from the settings' views (read once; the mean over T of each wiring's error)."""
    from experiments import E143_megamat_fly_sampling as E143m
    ref = pd.read_csv(H.results_dir(E143) / 'splits_consistency.csv')
    names = [s.name for s in E143m.SETTINGS if s.role == 'setting']
    wo = [s.name for s in E143m.SETTINGS if s.role == 'setting' and s.source != 'E135']
    rows = []
    for st in LOSS_STAGES:
        sd = make_side(st, model_rows(SETTINGS_CELL, st, 'all', names=[s.name for s in E143m.SETTINGS], common=False,
                                      each=True, wirings=T_WIRINGS))
        n = len(sd.flies)
        sp = E143m.stage_splits(n)
        Wa = E143m.half_weights(n, [a for a, _ in sp])
        Wb = E143m.half_weights(n, [b for _, b in sp])
        ea, eb = view_errors(sd, Wa, STRENGTH), view_errors(sd, Wb, STRENGTH)
        va, vb = E143m.values_of(ea, sd.rows), E143m.values_of(eb, sd.rows)
        va.columns = [c.split(':', 1)[1] for c in va.columns]
        vb.columns = [c.split(':', 1)[1] for c in vb.columns]
        for i in range(len(sp)):
            r = ref[(ref['stage'] == st) & (ref['split'] == i)]
            got = E143m.consistency(va.loc[i, names], vb.loc[i, names])
            got_wo = E143m.consistency(va.loc[i, wo], vb.loc[i, wo])
            rows.append(dict(stage=st, split=i, rho=got, rho_e143=float(r['rho'].iloc[0]) if len(r) else np.nan,
                             rho_without_e135=got_wo,
                             rho_without_e135_e143=float(r['rho_without_e135'].iloc[0]) if len(r) else np.nan))
    df = pd.DataFrame(rows)
    d = np.nanmax(np.abs(np.r_[df['rho'] - df['rho_e143'], df['rho_without_e135'] - df['rho_without_e135_e143']]))
    _write_csv(df, 'identity_e143_splits_consistency.csv')
    return dict(identity="E143's consistency table on strength from the settings' views", n_rows=len(df),
                max_abs_diff=float(d), passed=bool(np.isfinite(d) and d <= 1e-9))


def identities(cells: Optional[Sequence[str]] = None, lock: Optional[Mapping] = None) -> pd.DataFrame:
    """The identities on every cell whose views exist: the rescore tables reproduced from the views by E136's own
    reading (IDENTITY_TOL), the consistencies (CONSISTENCY_IDENTITIES) from the rescore tables and from the views, E143's
    consistency table on strength (from the settings' views), and the read-once check (CHECK_TOL)."""
    rows = []
    for c in (cells or CELLS):
        ev = cell_evals(c).set_index('eval')
        held = [int(i) for i in ev.index[ev['held']]]
        allv = [int(i) for i in ev.index]
        rs = rescore_table(c)
        a, b = _tab(rs, 'strength', 'a'), _tab(rs, 'strength', 'b')
        r_all = _rho(a.loc[allv, 'stage_error'], b.loc[allv, 'stage_error'])
        r_held = _rho(a.loc[held, 'stage_error'], b.loc[held, 'stage_error'])
        if c == 'A1':
            want = [('A1_held', r_held), ('A1_all', r_all)]
        else:
            tab = pd.read_csv(cell_results(c) / 'consistency_strength.csv')
            r_tab = float(tab.loc[tab['statistic'] == 'stage_error', 'rho'].iloc[0])
            want = [(c, r_all)]
            rows.append(dict(cell=c, identity="E145's consistency_strength.csv against the rescore table",
                             value=r_all, expected=r_tab, passed=abs(r_all - r_tab) <= 1e-12, needs_views=False))
        for key, val in want:
            rows.append(dict(cell=c, identity=f'consistency {key} (rescore table)', value=val,
                             expected=CONSISTENCY_IDENTITIES[key],
                             passed=abs(val - CONSISTENCY_IDENTITIES[key]) <= 0.0005, needs_views=False))
        if not have_views(c):
            rows.append(dict(cell=c, identity='views', value=np.nan, expected=np.nan, passed=None, needs_views=True,
                             note=f'no views yet: --part views --cell {c} on the host that holds its stages'))
            continue
        summ = views_summary(c)
        complete = views_complete(c)
        if not complete:
            rows.append(dict(cell=c, identity='views', value=np.nan, expected=np.nan, passed=None, needs_views=True,
                             note=f"views incomplete: {len(summ.get('missing') or [])} evaluations missing on "
                                  f"{summ.get('host')}; --part views --cell {c} on the host that holds them"))
        idr = identity_rescore(c)
        rows.append(dict(cell=c, identity=idr['identity'], value=idr['max_abs_diff'], expected=IDENTITY_TOL,
                         passed=idr['passed'], needs_views=True, note=json.dumps(dict(n_evals=idr['n_evals'],
                         worst=idr['worst'], missing=summ.get('missing')))))
        have = sorted(int(x) for x in set(cell_views(c)['eval']))
        hv = [i for i in held if i in have]
        key = 'A1_held' if c == 'A1' else c
        for rd in ('own', 'all') if complete else ():
            val = consistency_from_views(c, hv if c == 'A1' else [i for i in allv if i in have], rd)
            rows.append(dict(cell=c, identity=f'consistency {key} from the views ({rd} reading)', value=val,
                             expected=CONSISTENCY_IDENTITIES[key], needs_views=True,
                             passed=abs(val - CONSISTENCY_IDENTITIES[key]) <= 0.0005 if rd == 'own' else None,
                             note='' if rd == 'own' else 'beside: the read-once reading'))
        chk = read_once_check(c)
        rows.append(dict(cell=c, identity="E143's read-once check", value=max(x['diff'] for x in chk['checks']),
                         expected=CHECK_TOL, passed=chk['passed'], needs_views=True,
                         note=json.dumps(_jsonable(chk))))
    if have_views(SETTINGS_CELL):
        e = identity_e143(lock or {})
        rows.append(dict(cell=SETTINGS_CELL, identity=e['identity'], value=e['max_abs_diff'], expected=1e-9,
                         passed=e['passed'], needs_views=True))
    else:
        rows.append(dict(cell=SETTINGS_CELL, identity="E143's consistency table on strength", value=np.nan,
                         expected=np.nan, passed=None, needs_views=True,
                         note='no views yet: --part views --cell settings'))
    df = pd.DataFrame(rows)
    df['kind'] = ['check' if i == "E143's read-once check" else 'waiting' if i in ('views',) or (
        i.startswith("E143's consistency table") and p is None) else 'beside' if p is None else 'identity'
        for i, p in zip(df['identity'], df['passed'])]
    return df


def part_step0() -> dict:
    """Step 0: the lock checked, the probe (step0_probe.csv), the prediction (prediction.json, written once, before any
    view is read; refused afterwards if missing) and the identities of the cells whose views exist (identities.csv),
    step0.json."""
    lock = require_lock()
    t0 = time.time()
    probe = probe_table()
    _write_csv(probe, 'step0_probe.csv')
    pp = results_dir() / 'prediction.json'
    any_views = any(views_dir().glob('*.parquet'))
    if not pp.exists():
        if any_views:
            raise RuntimeError('views exist and no prediction was written before them')
        _write_json(pp, dict(prediction(), written=_now(), host=_host(), before_any_view=True))
    ids = identities(lock=lock)
    _write_csv(ids, 'identities.csv')
    failed = sorted(set(ids.loc[(ids['kind'] == 'identity') & (ids['passed'] == False), 'cell']))  # noqa: E712
    read_once_failed = sorted(set(ids.loc[(ids['kind'] == 'check') & (ids['passed'] == False), 'cell']))  # noqa: E712
    rec = dict(part='step0', probe_reproduced=bool(probe['reproduced'].all()),
               probe_not_reproduced=probe.loc[~probe['reproduced'], ['statistic', 'cell']].to_dict('records'),
               identities_checked=sorted(set(ids.loc[ids['kind'].isin(['identity', 'check']) & ids['needs_views'],
                                                     'cell'])),
               identities_waiting=sorted(set(ids.loc[ids['kind'] == 'waiting', 'cell'])),
               cells_failed=failed, read_once_failed=read_once_failed, seconds=round(time.time() - t0, 1),
               host=_host(), time=_now(),
               lock_sha256=lock.get('sha256'))
    _write_json(results_dir() / 'step0.json', rec)
    provenance('step0', lock)
    return rec


def require_identities(cells: Sequence[str] = tuple(CELLS)) -> Dict[str, Dict[str, bool]]:
    """Per cell, from identities.csv: 'identities' (its views exist and every identity passed: the rescore reproduced
    and the consistencies) and 'read_once' (E143's read-once check passed); RuntimeError before step 0."""
    p = results_dir() / 'identities.csv'
    if not p.exists():
        raise RuntimeError('identities.csv is not written: run --part step0')
    ids = pd.read_csv(p)
    truth = lambda x: str(x) == 'True'
    out = {}
    for c in cells:
        x = ids[ids['cell'] == c]
        idn = x[x['kind'] == 'identity']
        chk = x[x['kind'] == 'check']
        waiting = (x['kind'] == 'waiting').any()
        out[c] = dict(identities=bool(len(idn) and not waiting and all(truth(v) for v in idn['passed'])
                                      and (idn['needs_views'].map(truth)).any()),
                      read_once=bool(len(chk) and all(truth(v) for v in chk['passed'])))
    return out


# ---------------------------------------------------------------------------
# the model sets and their errors against any flies

def read_once_ok(cell: str) -> bool:
    """Whether a cell is read once (E143's read-once check passed in it, identities.csv); the saved settings are read
    once (E143's own check passed on T, results/E143_megamat_fly_sampling/check.json)."""
    if cell == SETTINGS_CELL:
        return bool(_read_json(H.results_dir(E143) / 'check.json').get('passed'))
    return require_identities([cell])[cell]['read_once']


@dataclass(frozen=True)
class ModelSet:
    """A set of models whose errors are compared between halves of the flies.

    Attributes:
        name: 'A1', 'X', 'F', 'W', 'pool_S', 'pool_T', 'settings', 'settings_without_e135', 'inhadd1000',
            'family_norm'.
        members: per member (a cell of VIEW_CELLS, its model names: 'e<i>' or a setting's name).
        read_once: whether every member is read once (else E136's own reading, on the subsets it saved).
    """
    name: str
    members: Tuple[Tuple[str, Tuple[str, ...]], ...]
    read_once: bool

    @property
    def labels(self) -> List[str]:
        return [f'{c}:{n}' for c, ns in self.members for n in ns]


def held_names(cell: str) -> List[str]:
    """A cell's held evaluations whose views are saved, as 'e<i>'."""
    ev = cell_evals(cell)
    have = set(int(x) for x in set(cell_views(cell)['eval']))
    return [f'e{int(i)}' for i, h in zip(ev['eval'], ev['held']) if h and int(i) in have]


def model_sets(cells: Optional[Sequence[str]] = None, settings: bool = True) -> Dict[str, ModelSet]:
    """The sets of design 4 whose views exist: each cell's held evaluations, the two pools (read once only if both
    members are), the saved settings (role 'setting'), without E135's four, and E143's two calyx groups."""
    from experiments import E143_megamat_fly_sampling as E143m
    out: Dict[str, ModelSet] = {}
    have = [c for c in (cells or CELLS) if views_complete(c)]
    for c in have:
        out[c] = ModelSet(c, ((c, tuple(held_names(c))),), read_once_ok(c))
    for p, cs in POOLS.items():
        if all(c in have for c in cs):
            out[f'pool_{p}'] = ModelSet(f'pool_{p}', tuple((c, tuple(held_names(c))) for c in cs),
                                        all(read_once_ok(c) for c in cs))
    if settings and have_views(SETTINGS_CELL):
        ro = read_once_ok(SETTINGS_CELL)
        sets = {'settings': [s.name for s in E143m.SETTINGS if s.role == 'setting'],
                'settings_without_e135': [s.name for s in E143m.SETTINGS if s.role == 'setting' and s.source != 'E135'],
                **{g: list(v) for g, v in E143m.CALYX_GROUPS.items()}}
        for k, ns in sets.items():
            out[k] = ModelSet(k, ((SETTINGS_CELL, tuple(ns)),), ro)
    return out


_SIDES: Dict[tuple, object] = {}


def _side(cell: str, stage: str, reading: str, kind: str):
    """A cached side: an evaluation cell's 'common' rows, or the settings' each-wiring rows on T ('T') or on every
    wiring ('each')."""
    key = (cell, stage, reading, kind)
    if key not in _SIDES:
        if kind == 'common':
            m = model_rows(cell, stage, reading, common=True, each=False)
        elif kind == 'T':
            m = model_rows(cell, stage, reading, common=False, each=True, wirings=T_WIRINGS)
        else:
            m = model_rows(cell, stage, reading, common=False, each=True)
        _SIDES[key] = make_side(stage, m)
    return _SIDES[key]


def _readings_at(stage: str) -> Dict[Tuple[int, ...], str]:
    return {v: k for k, v in readings_of_stage(stage).items() if k != 'all'}


def set_errors(ms: ModelSet, stage: str, subsets: Sequence[Sequence[int]], views) -> np.ndarray:
    """References x models: each model's error at a stage on `views` against each subset of the stage's flies (its
    mean in all the flies' one-fly units); an evaluation's value its wirings' mean views', a setting's the mean over T
    of each wiring's. Read once where the set is, else by E136's own reading on the subsets it saved (NaN for the
    others)."""
    from experiments import E143_megamat_fly_sampling as E143m
    n = len(stage_flies(stage).flies)
    subsets = [tuple(sorted(int(i) for i in s)) for s in subsets]
    cols = []
    for cell, names in ms.members:
        kind = 'T' if cell == SETTINGS_CELL else 'common'
        if ms.read_once or cell == SETTINGS_CELL and read_once_ok(cell):
            sd = _side(cell, stage, 'all', kind)
            e = view_errors(sd, E143m.half_weights(n, subsets), views)
            cols.append(_by_name(e, sd, cell, names, kind))
            continue
        rd = _readings_at(stage)
        block = np.full((len(subsets), len(names)), np.nan)
        for k, s in enumerate(subsets):
            if s == tuple(range(n)):
                r = 'all'
            elif s in rd:
                r = rd[s]
            else:
                continue
            sd = _side(cell, stage, r, kind)
            block[k] = _by_name(view_errors(sd, E143m.half_weights(n, [s]), views), sd, cell, names, kind)[0]
        cols.append(block)
    return np.hstack(cols) if cols else np.zeros((len(subsets), 0))


def _by_name(e: np.ndarray, sd, cell: str, names: Sequence[str], kind: str) -> np.ndarray:
    """References x names from a side's errors: an evaluation's 'common' row, a setting's mean over T (NaN unless all
    four are read)."""
    from experiments import E143_megamat_fly_sampling as E143m
    if kind == 'common':
        idx = sd.rows.get_indexer([(f'{cell}:{n}', 'common') for n in names])
        out = np.where(idx[None, :] >= 0, e[:, np.maximum(idx, 0)], np.nan)
        return out
    v = E143m.values_of(e, sd.rows, T_WIRINGS)
    return np.column_stack([v[f'{cell}:{n}'].to_numpy() if f'{cell}:{n}' in v.columns else np.full(len(e), np.nan)
                            for n in names])


def stage_splits_of(stage: str) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """E143's splits of a stage's flies (half a, half b)."""
    from experiments import E143_megamat_fly_sampling as E143m
    return E143m.stage_splits(len(stage_flies(stage).flies))


def e136_split_index(stage: str) -> int:
    from experiments import E143_megamat_fly_sampling as E143m
    n = len(stage_flies(stage).flies)
    return E143m.split_index(n, E143m.e136_split(n))


def decided(rho: Sequence[float]) -> bool:
    """"Is decided" / "can be decided": DECIDED_MEDIAN or more at the median split and DECIDED_MIN or more on every
    split (E143's tells-apart; splits without a value, NaN, left out)."""
    r = np.asarray(rho, dtype=float)
    r = r[np.isfinite(r)]
    return bool(len(r) and np.median(r) >= DECIDED_MEDIAN and np.min(r) >= DECIDED_MIN)


def family_spread(stage: str, views) -> float:
    """A view's spread: the SD (ddof 1) over the family's eight wirings of its per-wiring error against all the
    stage's flies (read once)."""
    from experiments import E143_megamat_fly_sampling as E143m
    sd = _side(SETTINGS_CELL, stage, 'all', 'each')
    n = len(sd.flies)
    e = view_errors(sd, E143m.half_weights(n, [tuple(range(n))]), views)[0]
    vals = [e[sd.rows.get_loc(('settings:family', w))] for w in WIRINGS if ('settings:family', w) in sd.rows]
    return float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan


def stage_views(stage: str) -> Tuple[str, ...]:
    return VIEWS + ((MATRIX,) if stage == 'pn_dendrite' else ())


# ---------------------------------------------------------------------------
# design 3: the room

def fly_parts(stage: str) -> Dict[str, Tuple[np.ndarray, float, List[str]]]:
    """The flies' saved parts of a megamat stage: per part, (flies x entries, all the flies' one-fly unit, the
    entries), the flies in the population's order."""
    ft = flies_table(stage)
    names = list(dict.fromkeys(ft['fly'].astype(str)))
    out = {}
    for q, f in ft.groupby('part', sort=False):
        entries = list(dict.fromkeys(f['entry'].astype(str)))
        X = f.pivot(index='fly', columns='entry', values='value').reindex(index=names, columns=entries)
        out[str(q)] = (X.to_numpy(dtype=float), float(f['unit'].iloc[0]), entries)
    return out


def _subset_mean(X: np.ndarray, idx: Sequence[int]) -> np.ndarray:
    Y = X[list(idx)]
    fin = np.isfinite(Y)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(fin.any(axis=0), np.where(fin, Y, 0.0).sum(axis=0) / fin.sum(axis=0), np.nan)


def half_to_half(stage: str, a: Sequence[int], b: Sequence[int], views) -> float:
    """One half's mean view as the prediction for the other: its error against the other half, in all the flies'
    one-fly units, on `views` ('all': every view, four at the PN dendrites; symmetric in the halves)."""
    from kcq import population_fit as P
    fp = fly_parts(stage)
    errs = {}
    for v in (*P.VIEWS, P.NAMED_VIEW):
        if views != ALL and v not in views:
            continue
        pe = []
        for x in P.PARTS[v]:
            q = f'{v}/{x}'
            if q not in fp:
                continue
            X, unit, _ = fp[q]
            d = _subset_mean(X, a) - _subset_mean(X, b)
            ok = np.isfinite(d)
            pe.append(np.sqrt(np.mean(d[ok] ** 2)) / unit if ok.any() else np.nan)
        if pe:
            errs[v] = float(np.sqrt(np.mean(np.square(pe))))
    use = [errs[v] for v in P.VIEWS if v in errs]
    if views == ALL and P.NAMED_VIEW in errs:
        use.append(errs[P.NAMED_VIEW])
    return float(np.sqrt(np.mean(np.square(use)))) if use else np.nan


ROOM_MODELS: Tuple[Tuple[str, str], ...] = (('family', 'T'), ('family', 'S'), ('E136:A1', 'T'), ('E145_A1', 'T'),
                                            ('E145_X', 'S'), ('E145_F', 'T'), ('E145_W', 'S'))


def part_room() -> dict:
    """Design 3: per stage and view, the family (on T and on S), A1 (E136's reads on T) and each classed point (E145's
    cells on their read wirings) against each half of every split (read once; the mean over the wirings of each
    wiring's error), beside half to half; room.csv and a summary per stage and view (medians over the splits and
    E136's split)."""
    from experiments import E143_megamat_fly_sampling as E143m
    lock = require_lock()
    rows = []
    for st in LOSS_STAGES:
        sp = stage_splits_of(st)
        n = len(stage_flies(st).flies)
        sd = _side(SETTINGS_CELL, st, 'all', 'each')
        subsets = [a for a, _ in sp] + [b for _, b in sp] + [tuple(range(n))]
        W = E143m.half_weights(n, subsets)
        for v in (*stage_views(st), ALL):
            vs = v if v == ALL else (v,)
            e = view_errors(sd, W, vs)
            hh = [half_to_half(st, a, b, vs) for a, b in sp]
            for name, wh in ROOM_MODELS:
                ws = S_WIRINGS if wh == 'S' else T_WIRINGS
                vals = E143m.values_of(e, sd.rows, ws)
                col = f'settings:{name}'
                if col not in vals.columns or vals[col].isna().all():
                    continue
                x = vals[col].to_numpy()
                for i, (a, b) in enumerate(sp):
                    rows.append(dict(stage=st, view=v, model=name, wirings=wh, split=i, e136=i == e136_split_index(st),
                                     against_a=float(x[i]), against_b=float(x[len(sp) + i]),
                                     against_all=float(x[-1]), half_to_half=hh[i]))
    df = pd.DataFrame(rows)
    _write_csv(df, 'room.csv')
    summ = (df.groupby(['stage', 'view', 'model', 'wirings'], sort=False)
            .agg(median_a=('against_a', 'median'), median_b=('against_b', 'median'),
                 against_all=('against_all', 'first'), median_half_to_half=('half_to_half', 'median'))
            .reset_index())
    e136 = df[df['e136']].set_index(['stage', 'view', 'model', 'wirings'])[['against_a', 'against_b', 'half_to_half']]
    summ = summ.join(e136.add_prefix('e136_'), on=['stage', 'view', 'model', 'wirings'])
    _write_csv(summ, 'room_summary.csv')
    provenance('room', lock)
    return dict(rows=len(df), summary=summ.to_dict('records'))


# ---------------------------------------------------------------------------
# design 4: which views the flies decide

def sharing_split(ms: ModelSet, stage: str, subset: Sequence[int], part: str) -> pd.DataFrame:
    """A sharing part's error of every model of a set against one subset of the flies split into level, spread and order
    (E136's `strength_split` terms on the part's entries: with d the model's entries minus the subset's mean, |mean d|,
    |sd_model - sd_flies| and sqrt(2 sd_model sd_flies (1 - r)), each over the part's one-fly unit), read once or by
    the subset's own reading as `set_errors` reads it. One row per model."""
    n = len(stage_flies(stage).flies)
    s = tuple(sorted(int(i) for i in subset))
    rows = []
    for cell, names in ms.members:
        kind = 'T' if cell == SETTINGS_CELL else 'common'
        if ms.read_once or cell == SETTINGS_CELL:
            r = 'all'
        else:
            r = 'all' if s == tuple(range(n)) else _readings_at(stage).get(s)
        if r is None:
            rows += [dict(model=f'{cell}:{nm}', level=np.nan, spread=np.nan, order=np.nan) for nm in names]
            continue
        sd = _side(cell, stage, r, kind)
        X = sd.fly_parts[part]
        fin = np.isfinite(X[list(s)])
        with np.errstate(invalid='ignore'):
            mean = np.where(fin.any(axis=0), np.nansum(X[list(s)], axis=0) / np.maximum(fin.sum(axis=0), 1), np.nan)
        M = sd.model[part]
        unit = sd.units[part]
        terms = np.full((M.shape[0], 3), np.nan)
        for i in range(M.shape[0]):
            ok = np.isfinite(M[i]) & np.isfinite(mean)
            if ok.sum() < 2:
                continue
            a, m = M[i, ok], mean[ok]
            sa, sm = a.std(), m.std()
            cross = 2.0 * sa * sm - 2.0 * float(np.mean((a - a.mean()) * (m - m.mean())))
            terms[i] = [abs(float(np.mean(a - m))) / unit, abs(sa - sm) / unit, np.sqrt(max(cross, 0.0)) / unit]
        df = pd.DataFrame(terms, columns=['level', 'spread', 'order'])
        df['setting'] = [r_[0] for r_ in sd.rows]
        df['wiring'] = [r_[1] for r_ in sd.rows]
        for nm in names:
            lab = f'{cell}:{nm}'
            if kind == 'common':
                x = df[(df.setting == lab) & (df.wiring == 'common')]
            else:
                x = df[(df.setting == lab) & df.wiring.isin(T_WIRINGS)]
                x = x if len(x) == len(T_WIRINGS) else x.iloc[0:0]
            rows.append(dict(model=lab, **({k: float(x[k].mean()) for k in ('level', 'spread', 'order')} if len(x)
                                          else dict(level=np.nan, spread=np.nan, order=np.nan))))
    return pd.DataFrame(rows)


def part_views_decided() -> dict:
    """Design 4: consistency per stage and view over every model set (each cell's held evaluations, the two pools, the
    saved settings, without E135's, and E143's two calyx groups) on every split of the stage (views_decided.csv), with
    "is decided", "no reading" (the set's range against all the flies against NO_READING_SPREADS times the family's
    spread) and E136's split's value (views_decided_summary.csv); the sharing parts' level, spread and order against
    each half of E136's split with their consistency (sharing_split.csv)."""
    lock = require_lock()
    sets = model_sets()
    rows, srows, shrows = [], [], []
    for st in LOSS_STAGES:
        sp = stage_splits_of(st)
        n = len(stage_flies(st).flies)
        e136i = e136_split_index(st)
        subsets = [a for a, _ in sp] + [b for _, b in sp] + [tuple(range(n))]
        for v in (*stage_views(st), ALL):
            vs = v if v == ALL else (v,)
            spread = family_spread(st, vs)
            for k, ms in sets.items():
                e = set_errors(ms, st, subsets, vs)
                rho = [_rho(e[i], e[len(sp) + i]) for i in range(len(sp))]
                for i, r in enumerate(rho):
                    rows.append(dict(stage=st, view=v, set=k, split=i, e136=i == e136i, rho=r,
                                     n_models=int(np.isfinite(e[i]).sum())))
                allf = e[-1]
                rng = float(np.nanmax(allf) - np.nanmin(allf)) if np.isfinite(allf).any() else np.nan
                rr = np.asarray(rho, dtype=float)
                srows.append(dict(stage=st, view=v, set=k, read_once=ms.read_once, n_models=len(ms.labels),
                                  n_splits=int(np.isfinite(rr).sum()), of_splits=len(sp),
                                  median=float(np.nanmedian(rr)) if np.isfinite(rr).any() else np.nan,
                                  min=float(np.nanmin(rr)) if np.isfinite(rr).any() else np.nan,
                                  e136=float(rr[e136i]), decided=decided(rr), range_all_flies=rng, spread=spread,
                                  no_reading=bool(np.isfinite(rng) and rng < NO_READING_SPREADS * spread)))
        for part in ('sharing/who', 'sharing/how_strongly'):
            for k, ms in sets.items():
                h = e136_halves(st)
                ta, tb = sharing_split(ms, st, h['a'], part), sharing_split(ms, st, h['b'], part)
                for term in ('level', 'spread', 'order'):
                    shrows.append(dict(stage=st, part=part, set=k, term=term, rho=_rho(ta[term], tb[term]),
                                       median_a=float(ta[term].median()), median_b=float(tb[term].median())))
    df, summ, sh = pd.DataFrame(rows), pd.DataFrame(srows), pd.DataFrame(shrows)
    _write_csv(df, 'views_decided.csv')
    _write_csv(summ, 'views_decided_summary.csv')
    _write_csv(sh, 'sharing_split.csv')
    provenance('views_decided', lock, sets={k: dict(read_once=m.read_once, n=len(m.labels)) for k, m in sets.items()})
    return dict(sets=sorted(sets), decided=summ[summ['decided']][['stage', 'view', 'set']].to_dict('records'))


# ---------------------------------------------------------------------------
# design 5: which loss the flies decide

def joint_splits_for(read_once: bool) -> pd.DataFrame:
    """The joint splits a set is read on: E143's 200 joint splits with E136's (read once), or, by E136's own reading,
    E136's split at the PN dendrites and boutons with each of the claws' and somata's three splits (nine, E136's
    among them). Columns: joint, <stage> (an index into the stage's splits), e136."""
    from experiments import E143_megamat_fly_sampling as E143m
    sizes = {st: len(stage_flies(st).flies) for st in LOSS_STAGES}
    if read_once:
        return E143m.joint_splits(sizes, N_JOINT, JOINT_SEED)
    rows = []
    e136 = {st: e136_split_index(st) for st in LOSS_STAGES}
    for i in range(len(stage_splits_of('kc_claw'))):
        for j in range(len(stage_splits_of('kc_soma'))):
            r = dict(pn_dendrite=e136['pn_dendrite'], pn_bouton=e136['pn_bouton'], kc_claw=i, kc_soma=j)
            rows.append(dict(joint=len(rows), **r, e136=(i, j) == (e136['kc_claw'], e136['kc_soma'])))
    return pd.DataFrame(rows)


def term_values(ms: ModelSet, js: pd.DataFrame) -> Dict[Tuple[str, object], Tuple[np.ndarray, np.ndarray]]:
    """Per term (stage, views) of every loss, the joint splits x models errors against half a and half b of each
    joint split's split of that stage."""
    terms = {(s, v) for _, (_, t) in LOSSES.items() for s, v in t}
    out = {}
    for st, v in terms:
        sp = stage_splits_of(st)
        idx = js[st].to_numpy(dtype=int)
        e = set_errors(ms, st, [a for a, _ in sp] + [b for _, b in sp], v)
        out[(st, v)] = (e[idx], e[len(sp) + idx])
    return out


def loss_values(tv: Mapping, loss: str) -> Tuple[np.ndarray, np.ndarray]:
    """A loss's joint splits x models values against half a and half b (the root mean square of its terms)."""
    terms = LOSSES[loss][1]
    return rms([tv[t][0] for t in terms]), rms([tv[t][1] for t in terms])


def pick_and_regret(fit: np.ndarray, other: np.ndarray, start: int = 0) -> Tuple[int, float]:
    """The pick (the lowest of `fit`) and its regret against `other`: (other[pick] - min other) / (other[start] - min
    other)."""
    p = int(np.nanargmin(fit))
    lo = float(np.nanmin(other))
    den = float(other[start] - lo)
    return p, float((other[p] - lo) / den) if den > 0 else np.nan


def holds_part_at(tv: Mapping, j: int, half: int, k: int, start: int, cl: Mapping) -> bool:
    """Megamat's part of "holds" of model k against half `half` (0 a, 1 b) of joint split j: the megamat stage error on
    strength falls from the start's by more than FALL_BAR, no stage's strength error rises by more than its spread, and
    the claws' and somata's errors over all views not by more than theirs."""
    sp, spa = cl['spreads'], cl['beside']['spreads_all_views']
    l0 = rms([tv[(s, STRENGTH)][half][j] for s in LOSS_STAGES])
    if not l0[start] - l0[k] > FALL_BAR:
        return False
    for s in LOSS_STAGES:
        x = tv[(s, STRENGTH)][half][j]
        if x[k] - x[start] > sp[s]:
            return False
    for s in CALYX_STAGES:
        x = tv[(s, ALL)][half][j]
        if x[k] - x[start] > spa[s]:
            return False
    return True


def part_losses() -> dict:
    """Design 5: L0 to L6 over each cell's held evaluations and the two pools, on the joint splits (`joint_splits_for`)
    and E136's: per joint split the consistency (losses_joint.csv); per cell and loss the pick on its fit half, its regret
    against the other half, megamat's part of "holds" at the pick against the other half, and its four factors, on
    E136's split and over the joint splits (losses_picks.csv); the summary with "can be decided" and "transfers"
    (losses_summary.csv)."""
    lock = require_lock()
    cl = json.loads((_e136_dir() / 'classes.json').read_text())
    sets = {k: v for k, v in model_sets(settings=False).items()}
    jrows, prows, srows = [], [], []
    for k, ms in sets.items():
        js = joint_splits_for(ms.read_once)
        tv = term_values(ms, js)
        drawn = ~js['e136'].to_numpy(dtype=bool) if ms.read_once else np.ones(len(js), dtype=bool)
        e136j = int(np.flatnonzero(js['e136'].to_numpy(dtype=bool))[0])
        for L in LOSSES:
            la, lb = loss_values(tv, L)
            rho = np.array([_rho(la[j], lb[j]) for j in range(len(js))])
            for j in range(len(js)):
                jrows.append(dict(set=k, loss=L, joint=int(js['joint'].iloc[j]), e136=bool(js['e136'].iloc[j]),
                                  rho=float(rho[j])))
            rec = dict(set=k, loss=L, read_once=ms.read_once, n_models=len(ms.labels), n_joint=int(drawn.sum()),
                       median=float(np.nanmedian(rho[drawn])), min=float(np.nanmin(rho[drawn])),
                       e136=float(rho[e136j]), decided=decided(rho[drawn]))
            if k in CELLS:
                c = get_cell(k)
                fit_i, oth_i = (0, 1) if c.fit_half == 'a' else (1, 0)
                names = ms.members[0][1]
                if names[0] != 'e0':
                    raise RuntimeError(f'{k}: evaluation 0 is not held or not saved')
                ev = cell_evals(k).set_index('eval')
                L_ = (la, lb)
                regs, holds = [], []
                for j in range(len(js)):
                    p, rg = pick_and_regret(L_[fit_i][j], L_[oth_i][j])
                    hp = holds_part_at(tv, j, oth_i, p, 0, cl)
                    regs.append(rg)
                    holds.append(hp)
                    if j == e136j:
                        i = int(names[p][1:])
                        prows.append(dict(cell=k, loss=L, split='e136', pick=i, regret=rg, transfers=rg <= TRANSFERS_MAX_REGRET,
                                          holds_part_other_half=hp, fit_loss=float(L_[fit_i][j][p]),
                                          other_loss=float(L_[oth_i][j][p]), other_min=float(np.nanmin(L_[oth_i][j])),
                                          start_other=float(L_[oth_i][j][0]),
                                          **{f: float(ev.loc[i, f]) for f in ('aplpn', 'kckc', 'inhadd', 'pn_inhadd',
                                                                              'inhsc', 'loss_sim')}))
                regs, holds = np.asarray(regs), np.asarray(holds)
                rec.update(pick_e136=prows[-1]['pick'], regret_e136=prows[-1]['regret'],
                           transfers_e136=prows[-1]['transfers'], holds_part_e136=prows[-1]['holds_part_other_half'],
                           regret_median=float(np.nanmedian(regs[drawn])),
                           share_transfers=float(np.mean(regs[drawn] <= TRANSFERS_MAX_REGRET)),
                           share_holds_part=float(np.mean(holds[drawn])))
            srows.append(rec)
    jdf, pdf, sdf = pd.DataFrame(jrows), pd.DataFrame(prows), pd.DataFrame(srows)
    _write_csv(jdf, 'losses_joint.csv')
    _write_csv(pdf, 'losses_picks.csv')
    _write_csv(sdf, 'losses_summary.csv')
    provenance('losses', lock, sets={k: dict(read_once=m.read_once, n=len(m.labels)) for k, m in sets.items()})
    return dict(summary=sdf.to_dict('records'))


# ---------------------------------------------------------------------------
# design 6: the claws' cost

CLASSED: Dict[str, Tuple[str, Tuple[str, ...]]] = {
    'A1': ('E145_A1', T_WIRINGS), 'X': ('E145_X', S_WIRINGS), 'F': ('E145_F', T_WIRINGS), 'W': ('E145_W', S_WIRINGS)}


def per_wiring_errors(stage: str, name: str, wirings: Sequence[str], subset: Sequence[int], views,
                      reading: str = 'all') -> Dict[str, float]:
    """A saved setting's per-wiring error at a stage against a subset of the flies (read once, or by `reading`)."""
    from experiments import E143_megamat_fly_sampling as E143m
    sd = _side(SETTINGS_CELL, stage, reading, 'each')
    n = len(sd.flies)
    e = view_errors(sd, E143m.half_weights(n, [list(subset)]), views)[0]
    return {w: float(e[sd.rows.get_loc((f'settings:{name}', w))]) for w in wirings
            if (f'settings:{name}', w) in sd.rows}


def claws_references(stage: str) -> List[Tuple[str, Tuple[int, ...]]]:
    """The references of design 6 at a stage: E136's halves, each half of every split, and each fly alone."""
    n = len(stage_flies(stage).flies)
    h = e136_halves(stage)
    out = [('e136_a', h['a']), ('e136_b', h['b'])]
    for i, (a, b) in enumerate(stage_splits_of(stage)):
        out += [(f'split{i}_a', tuple(a)), (f'split{i}_b', tuple(b))]
    out += [(f'fly{i}', (i,)) for i in range(n)]
    return out


def part_claws() -> dict:
    """Design 6: at the claws and at the somata, each classed point against the family on its cell's read wirings
    (E136's `compare` with the stage's strength spread as scored) against E136's halves (by E136's own reading as
    scored, and read once beside), each half of every split and each fly alone (read once); every cell's held
    evaluations' strength error changes from evaluation 0 against the same references (the share worse than the
    spread); the odors on which the halves of the claw flies differ most (leads). claws.csv, claws_evals.csv,
    claws_odor_leads.csv."""
    from experiments import E136_stage_population_fit as F
    lock = require_lock()
    cl = json.loads((_e136_dir() / 'classes.json').read_text())
    rows, erows = [], []
    for st in CALYX_STAGES:
        spread = float(cl['spreads'][st])
        for ref, sub in claws_references(st):
            readings = ['all']
            if ref.startswith('e136_'):
                readings.append(sub_name(sub))
            for rd in readings:
                for cell, (name, ws) in CLASSED.items():
                    a = per_wiring_errors(st, name, ws, sub, STRENGTH, rd)
                    f = per_wiring_errors(st, 'family', ws, sub, STRENGTH, rd)
                    cmp_ = F.compare(a, f, spread, ws)
                    rows.append(dict(stage=st, reference=ref, flies=sub_name(sub), reading='read once' if rd == 'all'
                                     else "E136's own", cell=cell, classed=name, change=cmp_.get('change'),
                                     se=cmp_.get('se'), spread=spread, worse_than_spread=cmp_.get('worse_than_spread'),
                                     improves=cmp_.get('improves')))
        for k, ms in model_sets(settings=False).items():
            if k not in CELLS:
                continue
            refs = claws_references(st)
            e = set_errors(ms, st, [s for _, s in refs], STRENGTH)
            for r, (ref, sub) in enumerate(refs):
                ch = e[r] - e[r][0]
                erows.append(dict(stage=st, cell=k, reference=ref, flies=sub_name(sub), read_once=ms.read_once,
                                  n=int(np.isfinite(ch).sum()), share_worse=float(np.nanmean(ch > spread))
                                  if np.isfinite(ch).any() else np.nan,
                                  median_change=float(np.nanmedian(ch)) if np.isfinite(ch).any() else np.nan))
    # the odors on which the halves of the claw flies differ (leads): the strength parts' halves' means
    fp = fly_parts('kc_claw')
    lrows = []
    for i, (a, b) in enumerate(stage_splits_of('kc_claw')):
        for q in ('strength/rate', 'strength/mean'):
            if q not in fp:
                continue
            X, unit, entries = fp[q]
            d = (_subset_mean(X, a) - _subset_mean(X, b)) / unit
            for rank, j in enumerate(np.argsort(-np.abs(np.nan_to_num(d)))[:5]):
                lrows.append(dict(split=i, e136=i == e136_split_index('kc_claw'), half_a=sub_name(a),
                                  half_b=sub_name(b), part=q, rank=rank + 1, odor=entries[j], a_minus_b=float(d[j])))
    df, edf, ldf = pd.DataFrame(rows), pd.DataFrame(erows), pd.DataFrame(lrows)
    _write_csv(df, 'claws.csv')
    _write_csv(edf, 'claws_evals.csv')
    _write_csv(ldf, 'claws_odor_leads.csv')
    provenance('claws', lock)
    return dict(claws=claws_class(df))


def claws_class(df: pd.DataFrame) -> dict:
    """"The claws' cost is the read half's" from claws.csv (the claws on strength; E136's halves by E136's own reading,
    as scored), and "one fly's" or "two against two" (each claw fly alone, read once)."""
    x = df[(df.stage == 'kc_claw')]
    own = x[x.reading == "E136's own"]
    per = {}
    for c in CELLS:
        y = own[own.cell == c].set_index('reference')['worse_than_spread']
        per[c] = dict(a=bool(y.get('e136_a', False)), b=bool(y.get('e136_b', False)))
    one_half = {c: (v['a'] != v['b']) for c, v in per.items()}
    which = {c: ('a' if v['a'] else 'b') for c, v in per.items() if one_half[c]}
    read_half = bool(all(one_half.values()) and len(set(which.values())) == 1)
    flies = sorted({r for r in x.reference if r.startswith('fly')})
    fly_sign = {f: all(bool(x[(x.cell == c) & (x.reference == f) & (x.reading == 'read once')]['worse_than_spread']
                            .any()) for c in CELLS) for f in flies}
    n_fly = sum(fly_sign.values())
    kind = "one fly's" if read_half and n_fly == 1 else 'two against two' if read_half else None
    return dict(per_cell=per, read_half=read_half, half=(next(iter(which.values())) if read_half else None),
                fly_sign=fly_sign, kind=kind)


# ---------------------------------------------------------------------------
# design 7: without the boutons

def part_no_boutons() -> dict:
    """Design 7: the megamat stage error without the boutons (L4) over the saved settings on the 200 joint splits and
    E136's, with and without E135's four arms, each setting's value the mean over T of each wiring's root mean square of
    its stages' strength errors (E143's reading), with L0 beside (E143's own splits_joint_consistency.csv reproduced);
    E141's and E145's classes recomputed with the boutons left out (E142's `reclass` and reading)."""
    from experiments import E143_megamat_fly_sampling as E143m
    lock = require_lock()
    js = joint_splits_for(True)
    names = [s.name for s in E143m.SETTINGS if s.role == 'setting']
    wo = [s.name for s in E143m.SETTINGS if s.role == 'setting' and s.source != 'E135']
    per_stage = {}
    for st in LOSS_STAGES:
        sp = stage_splits_of(st)
        n = len(stage_flies(st).flies)
        sd = _side(SETTINGS_CELL, st, 'all', 'T')
        e = view_errors(sd, E143m.half_weights(n, [a for a, _ in sp] + [b for _, b in sp]), STRENGTH)
        idx = js[st].to_numpy(dtype=int)
        per_stage[st] = (pd.DataFrame(e[idx], columns=sd.rows), pd.DataFrame(e[len(sp) + idx], columns=sd.rows))
    cols = None
    for a, _ in per_stage.values():
        cols = a.columns if cols is None else cols.intersection(a.columns, sort=False)
    rows = []
    for loss, stages in (('L4', [s for s in LOSS_STAGES if s != 'pn_bouton']), ('L0', list(LOSS_STAGES))):
        ea = rms([per_stage[s][0][cols].to_numpy() for s in stages])
        eb = rms([per_stage[s][1][cols].to_numpy() for s in stages])
        va = E143m.values_of(ea, cols)
        vb = E143m.values_of(eb, cols)
        va.columns = [c.split(':', 1)[1] for c in va.columns]
        vb.columns = [c.split(':', 1)[1] for c in vb.columns]
        for j in range(len(js)):
            rows.append(dict(loss=loss, joint=int(js['joint'].iloc[j]), e136=bool(js['e136'].iloc[j]),
                             rho=E143m.consistency(va.loc[j, names], vb.loc[j, names]),
                             rho_without_e135=E143m.consistency(va.loc[j, wo], vb.loc[j, wo])))
    df = pd.DataFrame(rows)
    _write_csv(df, 'no_boutons_joint.csv')
    ref = pd.read_csv(H.results_dir(E143) / 'splits_joint_consistency.csv')
    l0 = df[df.loss == 'L0'].reset_index(drop=True)
    e143_diff = float(np.nanmax(np.abs(l0['rho'].to_numpy() - ref['rho'].to_numpy()))) if len(ref) == len(l0) \
        else np.nan
    summ = {}
    for loss, g in df.groupby('loss', sort=False):
        d = g[~g.e136]
        summ[loss] = dict(median=float(d['rho'].median()), min=float(d['rho'].min()),
                          median_without_e135=float(d['rho_without_e135'].median()),
                          min_without_e135=float(d['rho_without_e135'].min()),
                          e136=float(g.loc[g.e136, 'rho'].iloc[0]),
                          e136_without_e135=float(g.loc[g.e136, 'rho_without_e135'].iloc[0]),
                          orders=decided(d['rho']), orders_without_e135=decided(d['rho_without_e135']))
    rec = dict(summary=summ, e143_l0_reproduced_max_abs_diff=e143_diff,
               e143_l0_reproduced=bool(np.isfinite(e143_diff) and e143_diff <= 1e-9),
               reclasses=reclass_without_boutons())
    _write_json(results_dir() / 'no_boutons.json', rec)
    provenance('no_boutons', lock)
    return rec


def reclass_without_boutons() -> dict:
    """E141's and E145's classes recomputed with the megamat stage error over the three other stages (E142's `reclass`
    with the boutons out of the stage error, the stage comparisons and the guard, the spreads from the unrefitted
    family's eight wirings against E136's held-out half, as E142 read E136's and E139's), with the default reading's
    re-read beside (which must give the classes as scored)."""
    from experiments import E142_step0 as S2
    three = [s for s in LOSS_STAGES if s != 'pn_bouton']
    e136m = pd.read_csv(_e136_dir() / 'score_megamat.csv')
    fam_rows = e136m[(e136m.half == 'held') & (e136m.arm == 'family') & e136m.wiring.isin(WIRINGS)]
    fam = {str(r.wiring): {s: float(getattr(r, s)) for s in LOSS_STAGES} for r in fam_rows.itertuples()}
    out = {}
    for reading, stages, drop in (('without_boutons', three, ('pn_bouton',)), ('default', list(LOSS_STAGES), ())):
        f = S2.with_stage(fam, stages)
        spreads = dict(stage=float(np.std([v['stage'] for v in f.values()], ddof=1)),
                       **{s: float(np.std([v[s] for v in f.values()], ddof=1)) for s in stages})
        rec = dict(stages=stages, spreads=spreads, classes={})
        e141 = json.loads((H.results_dir(E141) / 'classes.json').read_text())
        m141 = S2.megamat_rows(H.results_dir(E141) / 'score_megamat.csv', 'variant')
        for a, orig in e141['classes'].items():
            if a not in m141:       # E141's arms whose reads stopped after their holds failed: no class as scored
                rec['classes'][f'E141:{a}'] = dict(cls='no class' if orig.get('reasons') else 'reads missing',
                                                   original_cls=orig.get('cls'), note='no megamat reads')
                continue
            av = S2.with_stage({w: v for w, v in m141[a].items() if w in T_WIRINGS}, stages)
            r = S2.reclass(orig, av, f, spreads, stages, e141.get('control_fall'), drop)
            rec['classes'][f'E141:{a}'] = dict(cls=r['cls'], original_cls=r['original_cls'])
        e145 = json.loads((H.results_dir(E145) / 'classes.json').read_text())
        m145 = pd.read_csv(H.results_dir(E145) / 'score_megamat.csv')
        for c, crec in e145['cells'].items():
            orig = crec.get('cls') or {}
            ws = list(crec.get('read_wirings') or [])
            half = crec.get('read_half')
            x = m145[(m145.cell == c) & (m145.half == half)]
            vals = {}
            for nm in ('family', crec.get('variant')):
                y = x[(x.name == nm) & x.wiring.isin(ws)]
                vals[nm] = {str(r.wiring): {s: float(getattr(r, s)) for s in LOSS_STAGES} for r in y.itertuples()}
            if len(vals.get(crec.get('variant')) or {}) < len(ws) or len(vals['family']) < len(ws):
                rec['classes'][f'E145:{c}'] = dict(cls='reads missing', original_cls=orig.get('cls'))
                continue
            av = S2.with_stage(vals[crec.get('variant')], stages)
            fv = S2.with_stage(vals['family'], stages)
            r = _reclass_on(S2, orig, av, fv, f, spreads, stages, e145.get('control_fall'), drop, ws)
            rec['classes'][f'E145:{c}'] = dict(cls=r['cls'], original_cls=r['original_cls'])
        out[reading] = rec
    out['default_check'] = {k: dict(original=v['original_cls'], reread=v['cls'], same=v['original_cls'] == v['cls'])
                            for k, v in out['default']['classes'].items()}
    return out


def _reclass_on(S2, orig: Mapping, av: Mapping, fv: Mapping, fam8: Mapping, spreads: Mapping, stages: Sequence[str],
                control_fall, drop: Sequence[str], wirings: Sequence[str]) -> dict:
    """E142's `reclass` on a cell's own read wirings (E136's T_WIRINGS set to them for the comparison)."""
    from experiments import E136_stage_population_fit as F
    old = F.T_WIRINGS
    F.T_WIRINGS = tuple(wirings)
    try:
        return S2.reclass(orig, av, fv, spreads, stages, control_fall, drop)
    finally:
        F.T_WIRINGS = old


# ---------------------------------------------------------------------------
# design 8: beside, classing nothing

def part_beside() -> dict:
    """Design 8, classing nothing: (1) megamat's sharing (and breadth) error at the somata against all the flies over
    the saved settings, beside their held-out panels' strength errors (the class tables' comparisons, the mean over T)
    and E139's price (beside_settings.csv); (2) each stage's views regressed on the four log factors over each cell's
    held evaluations, against each half of E136's split (beside_regressions.csv); (3) E145's classes recomputed with
    E147's spreads over model flies (beside_e147_spreads.csv); (4) the sharing and breadth of E147's saved populations
    beside the reading at the lab's input (beside_e147_populations.csv)."""
    from experiments import E143_megamat_fly_sampling as E143m
    from experiments import E147_megamat_input_flies as E147m
    lock = require_lock()
    out = {}
    # (1)
    n = len(stage_flies('kc_soma').flies)
    sd = _side(SETTINGS_CELL, 'kc_soma', 'all', 'T')
    W = E143m.half_weights(n, [tuple(range(n))])
    vals = {v: E143m.values_of(view_errors(sd, W, (v,)), sd.rows).iloc[0] for v in VIEWS}
    tables = {E136: json.loads((_e136_dir() / 'classes.json').read_text()),
              E139: json.loads((H.results_dir(E139) / 'classes.json').read_text()),
              E141: json.loads((H.results_dir(E141) / 'classes.json').read_text())}
    rows = []
    for s in E143m.SETTINGS:
        r = dict(setting=s.name, role=s.role, calyx=s.calyx,
                 **{f'soma_{v}': float(vals[v].get(f'settings:{s.name}', np.nan)) for v in VIEWS})
        rec = (tables.get({'E136': E136, 'E139': E139, 'E141': E141}.get(s.source, ''), {}).get('classes') or {}).get(s.arm)
        if rec:
            for k, c in (rec.get('panels') or {}).items():
                r[f'panel_{k}'] = c.get('arm')
                r[f'panel_{k}_change'] = c.get('change')
            pr = rec.get('price_on_apl_effect') or {}
            r['e139_price_value'] = pr.get('price_value')
            r['e139_price'] = pr.get('price')
        rows.append(r)
    _write_csv(pd.DataFrame(rows), 'beside_settings.csv')
    # (2)
    reg = []
    for k, ms in model_sets(settings=False).items():
        if k not in CELLS:
            continue
        ev = cell_evals(k).set_index('eval')
        idx = [int(nm[1:]) for nm in ms.members[0][1]]
        Z = np.column_stack([np.ones(len(idx))] + [np.log(ev.loc[idx, f].to_numpy(dtype=float))
                                                    for f in ('aplpn', 'kckc', 'inhadd', 'inhsc')])
        for st in LOSS_STAGES:
            h = e136_halves(st)
            for v in stage_views(st):
                e = set_errors(ms, st, [h['a'], h['b']], (v,))
                for hi, half in enumerate('ab'):
                    y = e[hi]
                    ok = np.isfinite(y)
                    if ok.sum() < 6:
                        continue
                    beta, *_ = np.linalg.lstsq(Z[ok], y[ok], rcond=None)
                    res = y[ok] - Z[ok] @ beta
                    r2 = 1.0 - float(res @ res) / float(((y[ok] - y[ok].mean()) ** 2).sum())
                    reg.append(dict(cell=k, stage=st, view=v, half=half, n=int(ok.sum()), r2=r2,
                                    **{f'b_{f}': float(b) for f, b in zip(('const', 'aplpn', 'kckc', 'inhadd', 'inhsc'),
                                                                         beta)}))
    _write_csv(pd.DataFrame(reg), 'beside_regressions.csv')
    # (3)
    e147 = json.loads((H.results_dir(E147) / 'classes.json').read_text())
    spm = e147.get('spreads_model_flies_megamat') or {}
    e145 = json.loads((H.results_dir(E145) / 'classes.json').read_text())
    tab = {E145: {'classes': {c: v['cls'] for c, v in e145['cells'].items() if v.get('cls')}}}
    try:
        cw = E147m.classes_with_spreads(spm.get('strength') or {}, spm.get('guard') or {}, tables=tab)
        _write_csv(cw, 'beside_e147_spreads.csv')
        out['e147_spreads'] = cw.to_dict('records')
    except Exception as e:      # a description only: recorded, the other items go on
        out['e147_spreads'] = f'not recomputed: {type(e).__name__}: {e}'
    # (4)
    pp = H.results_dir(E147) / 'score_populations.csv'
    if pp.exists():
        p = pd.read_csv(pp)
        keep = [c for c in p.columns if c in ('population', 'reading', 'stage', 'n_model_flies') or
                any(c.endswith(x) for x in ('view_strength', 'view_sharing', 'view_breadth', 'view_matrix'))]
        _write_csv(p[keep], 'beside_e147_populations.csv')
        out['e147_populations'] = int(len(p))
    provenance('beside', lock, e147_classes_sha256=_sha256(H.results_dir(E147) / 'classes.json'),
               e147_populations_sha256=_sha256(pp) if pp.exists() else None)
    return out


# ---------------------------------------------------------------------------
# the score

def part_score() -> dict:
    """The classes and outcome letters (classes.json) from the parts' tables: per cell whether its identities passed
    and how it is read; "megamat reads the calyx beside strength"; whether a loss with sharing (L1, L2, L3) "can be
    decided" and whether its picks transfer; "the claws' cost is the read half's" ("one fly's" or "two against two");
    "without the boutons the megamat stage error orders the settings"; the outcome letters."""
    lock = require_lock()
    ids = require_identities()
    summ = pd.read_csv(results_dir() / 'views_decided_summary.csv')
    ls = pd.read_csv(results_dir() / 'losses_summary.csv')
    claws = claws_class(pd.read_csv(results_dir() / 'claws.csv'))
    nb = _read_json(results_dir() / 'no_boutons.json')
    ok_cells = [c for c in CELLS if ids[c]['identities']]
    sets_needed = [*ok_cells, 'pool_S', 'pool_T', 'inhadd1000']
    reads = []
    for st in CALYX_STAGES:
        for v in ('sharing', 'breadth'):
            x = summ[(summ.stage == st) & (summ.view == v)].set_index('set')
            have = [s for s in sets_needed if s in x.index]
            if len(have) == len(sets_needed) and all(str(x.loc[s, 'decided']) == 'True' for s in sets_needed):
                reads.append(f'{st}/{v}')
    over_fits = []
    for st in CALYX_STAGES:
        for v in ('sharing', 'breadth'):
            x = summ[(summ.stage == st) & (summ.view == v)].set_index('set')
            fs = [*ok_cells, 'pool_S', 'pool_T']
            if all(s in x.index and str(x.loc[s, 'decided']) == 'True' for s in fs):
                over_fits.append(f'{st}/{v}')
    any_calyx_decided = bool((summ.stage.isin(CALYX_STAGES) & summ.view.isin(VIEWS)
                              & (summ['decided'].astype(str) == 'True')).any())
    fit_sets = [*ok_cells, 'pool_S', 'pool_T']
    over_fits_median = []
    for st in CALYX_STAGES:
        for v in ('sharing', 'breadth'):
            x = summ[(summ.stage == st) & (summ.view == v)].set_index('set')
            if all(s in x.index and float(x.loc[s, 'median']) >= DECIDED_MEDIAN for s in fit_sets):
                over_fits_median.append(f'{st}/{v}')
    loss_dec = {}
    for L in LOSSES:
        x = ls[ls.loss == L].set_index('set')
        need = [*ok_cells, 'pool_S', 'pool_T']
        loss_dec[L] = dict(decided=bool(all(s in x.index and str(x.loc[s, 'decided']) == 'True' for s in need)),
                           transfers={c: bool(str(x.loc[c, 'transfers_e136']) == 'True') for c in ok_cells
                                      if c in x.index},
                           regret_e136={c: float(x.loc[c, 'regret_e136']) for c in ok_cells if c in x.index})
    sharing_losses = [L for L in ('L1', 'L2', 'L3') if loss_dec[L]['decided']]
    transfers = [L for L in sharing_losses if all(loss_dec[L]['transfers'].values())]
    no_b = (nb.get('summary') or {}).get('L4', {})
    orders = bool(no_b.get('orders_without_e135'))
    failed = [c for c in CELLS if not ids[c]['identities']]
    letters = []
    if reads and sharing_losses and transfers:
        letters.append('a')
    elif (over_fits or over_fits_median) and not reads:
        letters.append('b')
    if not any_calyx_decided:
        letters.append('c')
    if claws['read_half']:
        letters.append('d')
    if orders:
        letters.append('e')
    if failed or not all(ids[c]['read_once'] for c in CELLS):
        letters.append('f')
    rec = dict(cells={c: dict(ids[c], read_by="read once" if ids[c]['read_once'] else
                              "E136's own reading on fewer splits") for c in CELLS},
               no_class=failed, reads_calyx=dict(holds=bool(reads), stage_views=reads,
                                                 decided_over_the_fits=over_fits,
                                                 median_over_the_fits=over_fits_median),
               losses=loss_dec, loss_with_sharing_decided=sharing_losses, picks_transfer=transfers,
               claws=claws, no_boutons=dict(orders=orders, L4=no_b, L0=(nb.get('summary') or {}).get('L0')),
               outcomes=letters, outcome_text={k: OUTCOMES[k] for k in letters}, time=_now(), host=_host(),
               lock_sha256=lock.get('sha256'))
    complete = {c: views_complete(c) for c in VIEW_CELLS}
    rec['views_complete'] = complete
    rec['complete'] = all(complete.values())
    name = 'classes.json' if rec['complete'] else 'classes_partial.json'
    if not rec['complete']:
        rec['note'] = ('provisional: the views of ' + ', '.join(c for c, v in complete.items() if not v) + ' are not '
                       'complete, so these classes and letters read the other cells only; classes.json waits for them')
    _write_json(results_dir() / name, rec)
    provenance('score', lock, written=name)
    return rec


PARTS: Tuple[str, ...] = ('lock', 'views', 'step0', 'room', 'views_decided', 'losses', 'claws', 'no_boutons', 'beside',
                          'score')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--cell', choices=VIEW_CELLS, help='--part views: the cell whose views to read')
    ap.add_argument('--evals', nargs='+', type=int, help='--part views: only these evaluations')
    ap.add_argument('--outputs-root', nargs='+', default=[],
                    help="--part views: other checkouts' outputs/ folders to find the cell's saved stages in")
    a = ap.parse_args(argv)
    EXTRA_OUTPUTS[:] = [Path(x) for x in a.outputs_root]
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
