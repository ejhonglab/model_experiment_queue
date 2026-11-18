#!/usr/bin/env python3
"""E108: paa's input at the flies' concentration. paa's ORN input moved from the recorded 10^-4 to the flies' 10^-3 by
rules fixed before any run, read in the best family at fixed scales with both response rates re-held.

STUB (interface first): signatures and docstrings; the bodies follow.

The entry is item E108 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08 and 16:19; moved into the list 16:34; claimed by hal 16:37).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part B --lane a
    ./venv/bin/python -u experiments/E108_paa_extrapolation.py --part score
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E108_paa_extrapolation'


# ---------------------------------------------------------------------------
# the rules (made-up data in the tests)

def rule_x1(lo: np.ndarray, hi: np.ndarray, cap: float) -> np.ndarray:
    """X1, linear: hi plus its change from lo, capped at `cap`."""
    raise NotImplementedError


def rule_x2(lo: np.ndarray, hi: np.ndarray, ceiling: float) -> np.ndarray:
    """X2, saturating: hi moved toward `ceiling` by the fraction of the remaining distance covered from lo to hi."""
    raise NotImplementedError


def fit_x3(triples: pd.DataFrame) -> np.ndarray:
    """X3's coefficients (a, b, c) of top = a + b mid + c low, least squares over the triples ('low', 'mid', 'top')."""
    raise NotImplementedError


def rule_x3(lo: np.ndarray, hi: np.ndarray, coef: np.ndarray) -> np.ndarray:
    """X3: a + b hi + c lo."""
    raise NotImplementedError


def rule_scores(triples: pd.DataFrame, cap: float) -> pd.DataFrame:
    """The score that chooses the rule: per rule, the RMSE of the predicted top concentration over the triples (X1, X2
    on every triple; X3 with one odor's ramp left out at a time)."""
    raise NotImplementedError


def carrying_rule(scores: pd.DataFrame) -> str:
    """The rule with the lowest RMSE."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the data

def ramp_rates() -> pd.DataFrame:
    """olfdb dataset 21's ramps in estimated spike deltas (pyolfsysm's dF/F-to-spike fit): long table (odor, log10
    concentration, glomerulus, value = mean over flies)."""
    raise NotImplementedError


def ramp_triples(rates: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """The three triples of adjacent decades per single-odor ramp and glomerulus of the model's 36."""
    raise NotImplementedError


def paa_rates() -> pd.DataFrame:
    """paa at 10^-5 (diagnostic) and 10^-4 (the panel) in the model's units, the 36 glomeruli ('lo', 'hi')."""
    raise NotImplementedError


def input_of(name: str) -> pd.DataFrame:
    """An input by name: 'I0' (the lab's), 'X1', 'X2', 'X3', 'W5' (paa at 10^-5), 'P1'..'P4' (X1's changes permuted),
    'S1' (paa scaled to X1's total), 'B<f>' (validation2 fly f's B-myr)."""
    raise NotImplementedError


def write_inputs() -> dict:
    """inputs/<name>.csv with sha256, scores.csv, rule.json, the compact data copy in the shared data/."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# readings and the gate

def rho_pairs_with(model: pd.Series, data: pd.Series, odor: str) -> float:
    """ρ on odor-odor distances over the pairs with `odor` (temporary until E106's subset function exists)."""
    raise NotImplementedError


def identity() -> dict:
    """Part 0: every wiring's fit at the lab's input reproduces its saved counts (0 differences)."""
    raise NotImplementedError


def part_b(lane: Optional[str] = None) -> pd.DataFrame:
    """Part B: the held readings (resumable)."""
    raise NotImplementedError


def readings() -> pd.DataFrame:
    """Every input's readings per unit (wiring, population) and readout."""
    raise NotImplementedError


def gate(rd: pd.DataFrame, runs: pd.DataFrame, name: str) -> dict:
    """The gate of the entry for input `name` against the lab's input."""
    raise NotImplementedError


def class_of(g: dict) -> str:
    """'better', 'a price' or 'no gain'."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True)
    ap.add_argument('--lane', default=None)
    ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
