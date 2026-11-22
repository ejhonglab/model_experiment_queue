#!/usr/bin/env python3
"""E126: the ladder of the pooled survivor: one neuron per KC type along the dose and KC groups without type.

Stub (interface first): signatures and docstrings; the bodies follow. The entry is item E126 of
experiment_priorities.md; it builds on E125 (experiments/E125_kckc_survivor.py), whose runner, check and classes it
reuses by import.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E126_kckc_survivor_ladder'
STEP1_ARMS = ('TY75', 'TY25')
GROUP_SEEDS = (0, 1)
STEP2_ARMS = tuple(f'RG75s{s}' for s in GROUP_SEEDS)


def form_of(arm: str) -> str:
    """'TY', 'RG' (or 'Z' for Z0)."""
    raise NotImplementedError


def dose_of(arm: str) -> float:
    """d of an arm (TY25 -> 0.25, RG75s0 -> 0.75; Z0 -> 0)."""
    raise NotImplementedError


def seed_of(arm: str) -> Optional[int]:
    """The draw's seed of a random-group arm (RG75s1 -> 1), else None."""
    raise NotImplementedError


def random_groups(types: np.ndarray, seed: int) -> np.ndarray:
    """Group labels: the KCs dealt at random into groups of the KC types' sizes (the type labels permuted)."""
    raise NotImplementedError


def group_survivor(ref: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """One neuron per group: within each group every entry the group's total incoming weight of `ref` over its size
    squared (with labels = the types, E122's TY exactly)."""
    raise NotImplementedError


def condition_arrays(arm: str, ref: np.ndarray, types: np.ndarray) -> Dict[str, np.ndarray]:
    """{'on': (1 - d) ref + d X, 'off': d X} of an arm (E125's for TY and Z0)."""
    raise NotImplementedError


def part_check() -> dict:
    """E125's check: d 0 (Z0) on the two-simulator runner reproduces E123's S1_HBab_real exactly; the matrices."""
    raise NotImplementedError


def part_runs(step: int, arms: Sequence[str], wirings: Sequence[str], lane: Optional[str] = None) -> pd.DataFrame:
    """Step 1 (TY25, TY75) or step 2 (RG75s0, RG75s1; refused unless the gate is open), resumable, per lane."""
    raise NotImplementedError


def gate(classes: Mapping[str, dict]) -> dict:
    """Step 2's gate: TY75 keeps the gain by E125's class (no corner) with no threshold under 60 Hz on any wiring."""
    raise NotImplementedError


def differ(rd: pd.DataFrame, a: str, b: str) -> dict:
    """The "differ" rule: one form lower on at least 7 of 8 wirings (all recordings) and on S and T against all five
    targets."""
    raise NotImplementedError


def score() -> dict:
    """readings.csv, runs.csv, verdict.json, comparison_table.md, provenance.json."""
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
