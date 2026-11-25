#!/usr/bin/env python3
"""E137, the held-out panels read like for like: the mixture recordings by position in the session, and validation2 on
all its flies and odors (no simulation).

The entry is item E137 of experiment_priorities.md (proposed by viki 2026-11-24 22:54 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-25 02:08 under the idle-host rule). "The family", the views,
one-fly units and the floor are E134's (kcq/population_fit.py, kcq/stage_readings.py, experiments/
E134_population_baseline.py). No model is run and no model parameter changes: E134's saved populations and runs are
re-read and re-scored.

Parts (CLI `--part`, one per design item):
- lock (design items 1 and 3): for the eleven mixture flies of olfdb dataset 8, Remy's per-trial processed files
  (xrds_mean_peak.nc: her kc_soma_nls folders on tensor, or her suite2p 'combined' respvec files on matrix) with their
  sha256, the stimulus order of each session (the file's own trial list and Remy's olfactometer config), and the
  reproduction of the uploaded table from them (which cells, trials and blocks it averages); a fly whose files are not
  found is listed and stays read as uploaded. The five readings of the flies (`READINGS`) are fixed here, before any
  score, each with its own one-fly units and floor. Written to results/E137_heldout_readings/lock.json and
  lock_flies.csv; the committed copy is experiments/E137_lock.json.
- position (design item 2): per fly and stimulus, the response rate and mean response of each presentation, the third
  over the first within a block, the second block over the first, the probe odors and solvent at the session's start
  and end; whether the decline is general or specific to the mixture's components.
- readings (design items 3 and 4): the flies' populations under each reading (baseline-noise units where Remy's frames
  are on hand, her session z otherwise, at the lock's criterion), with their floor and their most atypical fly.
- rescore (design items 4 and 5): E134's saved mixture runs (the family as fitted and held, and the five beside models)
  read with the flies' trial noise on the model's mixture somata, scored under each reading, with the flat profile and
  the ORN recordings as the prediction beside; every strength error split into level, spread and order (E136's
  `strength_split`).
- rise (design item 6): the mixture's rise over its three dilutions, the flies' and the models'.
- validation2 (design item 7): validation2 on all its 20 flies and 22 odors under the two candidate rules for partial
  panels (`V2_RULES`), the six-fly reading of E134 and the two noise groups beside.
- order (design item 8): for megamat and validation2, how much of the one-fly unit an odor's place in its block carries;
  for the APL16 recordings, whether Yang's processed tables hold the order (otherwise a question for Yang).
- classes: the classes the entry fixes beforehand (`classify`) and the outcome they point to.

The public interface other experiments read (E138's mixture rows; keep these signatures stable):
- `load_lock()`: the committed lock, per fly its files with their sha256, its stimulus order, its reproduction and
  the readings it can be read under;
- `fly_mixture_view(panel, reading)`: the flies' populations of one mixture panel's KC somata under a named reading
  (one `kcq.population_fit.Population` per fly, with the trials where on hand), as the re-score reads them.

Holds (CLAUDE.md): only Remy's processed files (her per-trial respvec files and stimulus configs; no raw movies), Yang's
processed products only; no connectome table; no olfsysm change; nothing fitted to a held-out panel (the one fitted
quantity is a scale per fly among validation2's own flies, one of two candidate readings of the flies).

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part lock
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part position
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part readings
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part rescore
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part rise
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part validation2
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part order
    nice -n 10 ./venv/bin/python -u experiments/E137_heldout_readings.py --part classes
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E137_heldout_readings'
SCRIPT = 'experiments/E137_heldout_readings.py'
HERE = Path(__file__).resolve().parents[1]
LOCK_COPY = HERE / 'experiments' / 'E137_lock.json'

MIXTURE_PANELS: Tuple[str, ...] = ('kmix', 'cmix')
# the five readings of the flies (design item 3), fixed in the lock before any score
READINGS: Tuple[str, ...] = ('as_uploaded', 'first_presentations', 'components_first_block',
                             'components_second_block', 'components_both_blocks')
POSITION_READINGS: Tuple[str, ...] = READINGS[1:]
# validation2's two candidate rules for flies with partial panels (design item 7)
V2_RULES: Tuple[str, ...] = ('common_mean', 'fly_scale')
V2_COMMON_ODORS: Tuple[str, ...] = ('paa', 'PEA', 'guai')


def out_dir() -> Path:
    """<shared>/results/E137_heldout_readings/ (created)."""
    from kcq import harness as H
    return H.results_dir(EXPERIMENT)


# ---------------------------------------------------------------------------
# the public interface (E138 reads it)

def load_lock(path: Optional[Path] = None) -> dict:
    """The committed lock (experiments/E137_lock.json, or `path`): 'flies' maps each mixture fly ('<date>/<fly_num>')
    to its panel, its files (path, bytes, sha256), its session's stimulus order, whether the uploaded table is
    reproduced from them (and which cells, trials and blocks it averages), and the readings it is read under
    ('as_uploaded' only for a fly whose files are missing); 'readings' defines each of `READINGS`; 'criteria' the
    response criterion of each unit. Raises FileNotFoundError before the lock is written, and PermissionError if the
    copy has uncommitted changes."""
    raise NotImplementedError


def fly_mixture_view(panel: str, reading: str, trials: bool = True, lock: Optional[Mapping] = None):
    """The flies' KC somata of one mixture panel ('kmix' or 'cmix') under a named reading (one of `READINGS`), as the
    re-score reads them: a list of `kcq.population_fit.Population`, one per fly (named '<date>/<fly_num>'), cells x the
    panel's stimuli (data_prep.stage_population_data's KMIX_ODORS / CMIX_ODORS), in the lock's units and at its
    criterion; with `trials`, also the single trials (returned as the second item: {fly: DataFrame with columns
    (odor, repeat)}, None for a fly read as uploaded). Reads only the lock's files (sha256 checked)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design item 1: the lock

def mixture_flies() -> pd.DataFrame:
    """The eleven flies of olfdb dataset 8 (fly '<date>/<fly_num>', panel, number of cells)."""
    raise NotImplementedError


def find_fly_files(fly: str, panel: str) -> dict:
    """Remy's per-trial processed files of one fly's mixture session (respvec xrds_mean_peak.nc on tensor or matrix),
    its stimulus list and olfactometer config; {} entries where not found."""
    raise NotImplementedError


def reproduce_uploaded(fly: str, panel: str, files: Mapping) -> dict:
    """Whether the uploaded table (olfdb dataset 8) of one fly is reproduced from its per-trial file: the cell set,
    trials and blocks averaged, and the largest absolute difference."""
    raise NotImplementedError


def part_lock() -> dict:
    """Design items 1 and 3: lock.json and lock_flies.csv (results/), and the committed copy."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design item 2: the flies by position

def part_position() -> pd.DataFrame:
    """Per fly and stimulus, each presentation's response rate and mean response, the within-block and between-block
    ratios, the probes and solvent at the session's start and end; position_*.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design items 3 and 4: the readings of the flies

def part_readings() -> pd.DataFrame:
    """Each reading's flies: floor, most atypical fly and response rates per panel; readings_flies.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design items 4 and 5: the re-score

def flat_profile(ref) -> object:
    """The flat profile: every odor at the panel's mean, scored as a model's mean view."""
    raise NotImplementedError


def part_rescore() -> pd.DataFrame:
    """E134's saved mixture runs re-read with the flies' trial noise and scored under each reading, beside the flat
    profile and the ORN recordings; rescore_scores.csv, rescore_split.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design item 6: the rise

def part_rise() -> pd.DataFrame:
    """The mixture's rise over its three dilutions (response rate and mean response at 0 over -2), the flies' and the
    models'; rise.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design item 7: validation2 on all flies and odors

def v2_reference(flies, odors: Sequence[str], rule: str):
    """The flies' reference of validation2's KC somata on all flies, each fly on its own odors, under a partial-panel
    rule (`V2_RULES`)."""
    raise NotImplementedError


def part_validation2() -> pd.DataFrame:
    """Validation2 on 6 and on 20 flies under both rules, the two noise groups beside; validation2_*.csv."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# design item 8: the order at the other panels

def part_order() -> dict:
    """The share of the one-fly unit an odor's place in its block carries, at megamat and validation2; the APL16
    recordings' order, or the question for Yang; order_*.csv, order.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the classes

def classify(scores: pd.DataFrame) -> dict:
    """The entry's classes: a mixture panel 'readable on strength' if the family's class by E134's stage rule is the
    same under the four position readings; validation2 'readable on strength' if the class is the same on 6 and 20
    flies under both rules; a reading 'replaces' the lock's where readable and reproduced; the outcome (a) to (e)."""
    raise NotImplementedError


def part_classes() -> dict:
    """classes.json and provenance.json."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'position', 'readings', 'rescore', 'rise',
                                                      'validation2', 'order', 'classes'))
    a = ap.parse_args(argv)
    dict(lock=part_lock, position=part_position, readings=part_readings, rescore=part_rescore, rise=part_rise,
         validation2=part_validation2, order=part_order, classes=part_classes)[a.part]()


if __name__ == '__main__':
    main()
