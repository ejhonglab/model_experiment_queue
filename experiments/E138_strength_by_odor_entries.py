#!/usr/bin/env python3
"""E138, the odors that carry the strength error and the glomeruli that feed them: the strength view's entries odor by
odor at the KC somata of every panel, on the flies' halves, along the megamat stages, and the glomeruli behind the
odors that carry it (no simulation).

The entry is item E138 of experiment_priorities.md (proposed by viki 2026-11-24 22:54 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-25 02:08 under the idle-host rule; claimed by atlas 02:11).
"The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings (E83's), each read as one fly; the
views, one-fly units, readings and flies are E134's (kcq/population_fit.py, kcq/stage_readings.py,
data_prep/stage_population_data.py and its lock experiments/E134_lock.json). The strength view ("strength by odor") has
two parts, each odor's response rate and its mean response as multiples of the panel's mean (part 'rate' and part
'mean'); a part's error is the root mean square over odors of the model's entry minus the flies' mean entry, in the
flies' one-fly units, and the view's error the root mean square of its two parts' errors. Nothing is run, fitted or
searched: E134's, E107's and E90's saved runs are read again.

The guard: no residual (model minus flies) is computed before this experiment's classes lock exists, is committed with
no uncommitted change and equals its copy in the shared results (`require_classes_lock`); and E134's lock guards every
reading of the flies (`stage_readings.require_lock`). Every part but 'classes' calls both first.

Parts (CLI `--part`):
- classes (design item 1; no residual): the three classes of glomeruli, written to results/E138_strength_by_odor_entries/
  classes_lock.json and its committed copy experiments/E138_classes_lock.json (then commit it), with the rules of the
  classes read here ("concentrated", (d)) and the feeders' rule, fixed before any residual:
  - 'ir': the glomeruli fed by ORNs whose ligand-specific receptor is an ionotropic receptor (IR), from Task et al. 2022
    (eLife 11:e72599), the table of receptors and co-receptors by glomerulus as drosolf ships it (data/task22_table3.csv;
    the same table as al_analysis's data/sent_to_nikola/2025-02-25/task22_table3.csv), whose IR glomeruli are those of
    Silbering et al. 2011 (J Neurosci 31:13357); the co-receptors (Ir8a, Ir76b, Ir25a) beside, not a class;
  - 'not_dendrite28': the model's input glomeruli absent from the 28 glomeruli of the PN-dendrite recordings (E135's
    lock, `glomeruli`; E90's `shared_glomeruli`);
  - 'few_hemibrain': the model's input glomeruli with fewer boutons, or fewer claws, than the median input glomerulus in
    the family's hemibrain wiring (E130's `step0_wiring_only('hemibrain-first')`, the claw x bouton matrix of E83's
    'real'), each sub-list beside the union.
- entries (design items 2 and 7): per panel and condition, each odor's entry in both strength parts for the family as
  fitted and as held and for the beside models: the model's entry (each wiring, the wirings' mean views 'common'),
  the flies' mean entry and range over flies, the residual in one-fly units, each odor's share of each part's squared
  error and of the view's, and the view's error without that odor's entry; the strength split (E136's
  `strength_split`: level, spread and order) beside every error. Panels (`PANEL_SPECS`): megamat (E134's step 3,
  17 odors), the APL16 panel with APL on against the control flies and APL off against APL output blocked (TNT)
  (E90's part B, the family's own fit panel, and E90's D1 fits beside), validation2 (E107's part B for the family, E134's
  step 3 for the beside models), and the kiwi and control mixtures (E134's step 3) only once the held-out panels read
  like for like (E137) has its lock committed (`require_e137_lock`); until then the mixture rows are skipped and the
  outputs say so.
- halves (design items 3 and 7): the same on the two seeded halves of each panel's flies (E136's `split_names`: hal's
  `split_halves(n, 1, 0)`, the first part the 'fit' half; the APL16 panel's control half paired with the TNT half of the
  same name, as the E134 lock's split_halves rule pairs them), each half read as the lock's split_halves rule says
  (curve, noise, reference and one-fly units from the half), on the odors of the all-flies reading; the correlation over
  odors between the halves of the residual entries and of the shares. Remy's 14 older megamat flies with partial panels
  are not in viki's data module (stage_population_data.py leaves them out), so nothing goes beside.
- stages (design items 4 and 7): megamat along the stages: each odor's strength entries at the flies' ORNs, PN
  dendrites, PN boutons, KC claws and KC somata (mean and range over flies) beside the model's (the family as fitted and
  held, the beside models), each stage in its default reading.
- feeders (design item 5): for the odors that carry each panel's error (`carrying_odors`), the glomeruli that carry their
  ORN input (`feeders`: the model's input on its 36 glomeruli, and the panel's ORN recordings on every recorded
  glomerulus beside), each flagged for being among the model's 36, among the dendrites' 28, in each class of the lock,
  labelled by the bouton line (`BOUTON_LINE`; its published expression, `bouton_line_labels`), and its boutons and
  claws, counts and shares, in the hemibrain and in the male CNS's two calyces (Pratyush's tables through olfsysm's
  'malecns-left' / 'malecns-right' loaders, E130's `step0_wiring_only`).
- score: the classes read beforehand: "concentrated" per panel (`concentrated`) and (d) (`d_reading`); classes.json,
  tables.md, provenance_score.json.
- regression (design item 6): HELD. The regression of the family's residuals on the share of each odor's input on a
  class, read on validation2's odors against random classes, fits a quantity to a held-out panel, which waits for Tom's
  or Yang's word (CLAUDE.md, the entry's Status). The command line refuses it unless `--authorized-by` names that word
  (who, when, where); `class_regression` is the function, never run on the real residuals here.

The classes read here (fixed beforehand, in the lock's 'rules'):
- "concentrated": a panel's two largest odors (by their share of the squared strength error on all flies, the family as
  fitted, the wirings' mean views) carry more than half of the squared strength error on both halves of the flies;
- (d): the male CNS gives those glomeruli (the feeders of the panel's two largest odors) more boutons or claws: on both
  male calyces their summed share of the calyx's boutons, or of its claws, is above the hemibrain's. Read on megamat as
  the entry says, and given for every panel beside.
(a) and (b) need item 6 and are not read.

Concentrations: megamat at 10^-3 on both sides; the APL16 entries involve paa, euc and -aPine, whose model input is the
ORN data within 1 log10 unit of the flies' concentrations; validation2 is one panel at the same concentrations; the
mixtures are used as measured.

Holds (CLAUDE.md): only Pratyush's hemibrain and male CNS (olfsysm's loaders of his tables); no olfsysm change; no fit to
KC recordings; item 6 held. Nothing is published.

Outputs (<shared>/results/E138_strength_by_odor_entries/): classes_lock.json, classes_glomeruli.csv; entries.csv,
entries_split.csv, entries_summary.csv; halves.csv, halves_split.csv, halves_agreement.csv; stages.csv,
stages_split.csv; feeders.csv, feeders_glomeruli.csv; classes.json, tables.md; provenance_<part>.json.

Usage (from the repository root; atlas is shared, so two threads under nice):
    ./venv/bin/python -u experiments/E138_strength_by_odor_entries.py --part classes   # then commit the lock copy
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E138_strength_by_odor_entries.py \\
        --part entries   # likewise halves, stages, feeders, score
"""
from __future__ import annotations

import argparse
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E138_strength_by_odor_entries'
SCRIPT = 'experiments/E138_strength_by_odor_entries.py'
CLASSES_LOCK = 'classes_lock.json'
CLASSES_LOCK_COPY = H.HERE / 'experiments' / 'E138_classes_lock.json'
E137_LOCK_COPY = H.HERE / 'experiments' / 'E137_lock.json'

S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
STRENGTH_PARTS: Tuple[str, ...] = ('strength/rate', 'strength/mean')
HALVES: Tuple[str, ...] = ('fit', 'held')
HALF_SEED = 0
CLASSES: Tuple[str, ...] = ('ir', 'not_dendrite28', 'few_hemibrain')
CALYCES: Tuple[str, ...] = ('hemibrain-first', 'malecns-left', 'malecns-right')
MALE_CALYCES: Tuple[str, ...] = ('malecns-left', 'malecns-right')
# the bouton recordings' driver line (Remy's PN boutons: GH146-Gal4 > 20xUAS-IVS-jGCaMP8m, the thesis's methods)
BOUTON_LINE = 'GH146-Gal4'
# the class rules (fixed in the lock before any residual)
CONCENTRATED_TOP = 2
CONCENTRATED_SHARE = 0.5
# the feeders of an odor: the smallest set of glomeruli, by descending positive input, carrying this share of the odor's
# summed positive input over the model's 36 glomeruli
FEEDER_SHARE = 2.0 / 3.0
# the odors that carry a panel's error: its CONCENTRATED_TOP largest, and any other with at least CARRY_MIN_SHARE of the
# squared view error, at most CARRY_MAX
CARRY_MIN_SHARE = 0.10
CARRY_MAX = 5
# item 6 (held): random classes of the same size
N_RANDOM_CLASSES = 1000
REGRESSION_SEED = 0

# the panels read: (panel, condition) -> the groups read on it, each (group, loader key)
PanelSpec = namedtuple('PanelSpec', 'panel condition groups')


def out_dir() -> Path:
    """<shared>/results/E138_strength_by_odor_entries/ (created)."""
    return H.results_dir(EXPERIMENT)


# ---------------------------------------------------------------------------
# design item 1: the classes lock

def ir_glomeruli() -> dict:
    """The glomeruli fed by IR-tuned ORNs (Task et al. 2022's table, receptors starting 'Ir'), with the table's file,
    sha256 and citation, and each glomerulus's co-receptors beside."""
    raise NotImplementedError


def dendrite_glomeruli() -> List[str]:
    """The 28 glomeruli of the PN-dendrite recordings (E135's lock 'glomeruli')."""
    raise NotImplementedError


def wiring_counts(calyx: str) -> pd.DataFrame:
    """Per glomerulus of a calyx's claw x bouton wiring: its boutons and claws, as counts and as shares of the calyx's."""
    raise NotImplementedError


def classes_contents() -> dict:
    """The classes lock's contents (no residual is read)."""
    raise NotImplementedError


def part_classes() -> dict:
    """Design item 1: classes_lock.json (shared results) and its committed copy; classes_glomeruli.csv."""
    raise NotImplementedError


def require_classes_lock() -> dict:
    """The committed classes lock (tracked, no uncommitted change, the same as the shared copy); PermissionError
    otherwise."""
    raise NotImplementedError


def require_e137_lock() -> dict:
    """E137's committed lock, before any mixture row is read; PermissionError otherwise."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the entries

def odor_entries(pops: Sequence, ref, names: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Per population (and the wirings' mean views, 'common'), odor and strength part: the model's entry, the flies'
    mean entry and range, the residual in one-fly units, its share of the part's and of the view's squared error, and
    the view's error without that odor's entry."""
    raise NotImplementedError


def concentrated(shares_all: pd.Series, shares_halves: Mapping[str, pd.Series], top: int = CONCENTRATED_TOP,
                 share: float = CONCENTRATED_SHARE) -> dict:
    """The class "concentrated" of one panel: the `top` odors with the largest share on all flies carry more than
    `share` of the squared strength error on both halves."""
    raise NotImplementedError


def feeders(inputs: pd.Series, share: float = FEEDER_SHARE) -> pd.Series:
    """The glomeruli that carry an odor's input: the smallest set, by descending positive input, carrying `share` of its
    summed positive input."""
    raise NotImplementedError


def d_reading(glomeruli: Sequence[str], counts: Mapping[str, pd.DataFrame]) -> dict:
    """(d): whether both male calyces give `glomeruli` a larger summed share of boutons, or of claws, than the
    hemibrain."""
    raise NotImplementedError


def class_regression(residuals: pd.Series, class_share: pd.Series, read_residuals: pd.Series,
                     read_share: pd.Series) -> dict:
    """Design item 6 (held): the residuals regressed on the share of each odor's input on a class, and the fit's
    correlation with the read panel's residuals."""
    raise NotImplementedError


def part_entries() -> pd.DataFrame:
    raise NotImplementedError


def part_halves() -> pd.DataFrame:
    raise NotImplementedError


def part_stages() -> pd.DataFrame:
    raise NotImplementedError


def part_feeders() -> pd.DataFrame:
    raise NotImplementedError


def part_score() -> dict:
    raise NotImplementedError


def part_regression(authorized_by: Optional[str]) -> dict:
    """Design item 6: refused unless `authorized_by` names Tom's or Yang's word."""
    if not authorized_by:
        raise SystemExit('E138 design item 6 (the regression on a class, read on validation2) is held for Tom\'s or '
                         'Yang\'s word; pass --authorized-by "<who, when, where>" only once it is given')
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('classes', 'entries', 'halves', 'stages', 'feeders', 'score', 'regression'))
    ap.add_argument('--authorized-by', default=None,
                    help="regression only: Tom's or Yang's word that lifts the hold on design item 6 (who, when, where)")
    a = ap.parse_args(argv)
    if a.part == 'regression':
        part_regression(a.authorized_by)
        return
    dict(classes=part_classes, entries=part_entries, halves=part_halves, stages=part_stages, feeders=part_feeders,
         score=part_score)[a.part]()


if __name__ == '__main__':
    main()
