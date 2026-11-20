#!/usr/bin/env python3
"""E113: the input of paa and B-myr without a search: each ORN fly's own recording, and the three entries that two
independent searches both named (DM1 for both odors, DP1m for paa), read on the family's eight wirings.

STUB (interface first): signatures and docstrings; bodies raise NotImplementedError until the implementation lands.

The entry is item E113 of experiment_priorities.md as revised 2026-11-20 10:11 (cut to its no-search parts after E116's
outcome (d); planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-20 10:08). Canonical name: "the
searched input of paa, B-myr and 2-but inside the ORN flies' spread, reduced and fitted again".
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

EXPERIMENT = 'E113_input_within_spread'
SCRIPT = 'experiments/E113_input_within_spread.py'
PAA_BMYR = ('paa', 'B-myr')
ARMS = ('D2', 'D1', 'D1x')
CLASS_ARM = 'D2'
CONTROL_GLOMERULI = {'paa': ('VL2a', 'DM6', 'DL2v', 'VA3'), 'B-myr': ('VM5v', 'DL2v', 'DC1', 'DM4')}
CONTROLS = ('C1', 'C2', 'C3', 'C4', 'D1inv')
N_FLIES = 6
FLY_KINDS = ('paa', 'pb')


def out_dir() -> Path:
    """results/E113_input_within_spread/ (E113_OUT_DIR in the environment instead, for smoke tests)."""
    raise NotImplementedError


def local_dir(sub: str = '') -> Path:
    """outputs/E113_input_within_spread/<sub> (E113_LOCAL_DIR in the environment instead)."""
    raise NotImplementedError


def entry_sizes() -> Dict[str, Dict[Tuple[str, str], float]]:
    """{arm: {(odor, glomerulus): log2 factor}} of arm D and its controls, checked against E109's and E116's tables."""
    raise NotImplementedError


def control_glomeruli(lab: Optional[pd.DataFrame] = None) -> Dict[str, Tuple[str, ...]]:
    """The four glomeruli nearest DM1 in the lab's response of each odor (the exclusions of the module docstring)."""
    raise NotImplementedError


def factors_of(name: str) -> pd.DataFrame:
    """An arm's or control's factors, glomeruli x abbreviated odors (1 where unchanged)."""
    raise NotImplementedError


def fly_deviations() -> dict:
    """The six validation2 ORN flies' paa and B-myr in the model's units, as deviations from the flies' mean."""
    raise NotImplementedError


def fly_input(f: int, kind: str) -> pd.DataFrame:
    """Fly `f`'s input: 'paa' (paa alone), 'pb' (paa and B-myr) or 'bmyr' (B-myr alone, E108's B<f>)."""
    raise NotImplementedError


def input_of(name: str) -> pd.DataFrame:
    """An input by name ('I0', an arm, a control, 'F<f>_<kind>')."""
    raise NotImplementedError


def write_inputs() -> dict:
    """inputs/<name>.csv, sha256.json and inputs.json (each input's description and its entries in the flies' SDs)."""
    raise NotImplementedError


def identity(wirings: Optional[Sequence[str]] = None) -> dict:
    """Step 1: the held reading of the lab's input reproduces each wiring's saved fit (0 differences, both factors 1);
    E108's B0 on the real wiring reproduces its saved run (whether part 0 reuses E108's B-myr-alone runs)."""
    raise NotImplementedError


def plan(wirings: Optional[Sequence[str]] = None) -> List[dict]:
    """The held readings of the simulation part, in order: arm D, the controls, part 0 (dict(change, wiring, tag, kind))."""
    raise NotImplementedError


def part_sims(lane: Optional[str] = None, wirings: Optional[Sequence[str]] = None,
              max_points: Optional[int] = None) -> pd.DataFrame:
    """The held readings (resumable; `lane` 'i/n' takes every n-th point from the i-th)."""
    raise NotImplementedError


def readings(runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings.csv: per change, unit and readout, E108's `unit_scores`."""
    raise NotImplementedError


def gate_arm(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """The gate of one arm of step 3 on T, with the five controls."""
    raise NotImplementedError


def gate_flies(rd: pd.DataFrame, runs: pd.DataFrame, kind: str) -> dict:
    """Part 0's rule for one kind of fly input: the fly picked on S against one half, read on T against the other."""
    raise NotImplementedError


def outcomes(arms: Mapping[str, dict], flies: Mapping[str, dict]) -> List[str]:
    """The outcome letters (a) to (e)."""
    raise NotImplementedError


def refit_choice(arms: Mapping[str, dict], flies: Mapping[str, dict], rd: pd.DataFrame) -> Optional[str]:
    """The one input the refit runs with, or None (no pass)."""
    raise NotImplementedError


def part_refit(dry_run: bool = False) -> List[str]:
    """Step 5: the family fitted again with the chosen input on T's four wirings (only after a pass)."""
    raise NotImplementedError


def part_remy() -> pd.DataFrame:
    """Step 6, beside: D1 and D2 on the validation2 ORN input, T as fitted, against Remy's validation2 somata."""
    raise NotImplementedError


def score() -> dict:
    """Part score (no simulations)."""
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv for the model-output archive."""
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
