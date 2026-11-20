#!/usr/bin/env python3
"""E120: the bouton stage held to the recorded boutons. E119's per-ROI bouton statistics (the per-odor totals' CV, the
responders' CV, the top-1 % share, the share of entries below minus the threshold), at E119's criteria and readings, read
on saved runs (part A, no simulations), and the level of the antennal lobe's normalization for part B (E102's part D at
that level on the four searched wirings, run elsewhere).

STUB (interface first): signatures and docstrings; bodies raise NotImplementedError.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E120_bouton_stage_measured'


def indicator_change(win: pd.DataFrame, pre: pd.DataFrame, n: float, form: str) -> pd.DataFrame:
    """A stage's units x odors change seen through an indicator of exponent n (E98's family, forms 'baseline' and
    'evoked')."""
    raise NotImplementedError


def stats_rows(units: Dict[str, Tuple[pd.DataFrame, np.ndarray]], stage: str, readings: Sequence[str], pool: str,
               match_on: Optional[Dict[str, Tuple[pd.DataFrame, np.ndarray]]] = None, **keys) -> pd.DataFrame:
    """E119's per-ROI statistics of each wiring's units (E119's readings, pooling, matched threshold), with the flies'
    sides."""
    raise NotImplementedError


def part_a1() -> pd.DataFrame:
    """A1: the family's boutons and claws (E90's archive) through each member of the indicator family."""
    raise NotImplementedError


def part_a2() -> pd.DataFrame:
    """A2: E102's normalization ladder, the boutons with APL off from the saved PN stage, beside E102's cost."""
    raise NotImplementedError


def part_a3() -> pd.DataFrame:
    """A3: E94's saved bouton-dose settings on the four searched wirings, APL on, beside E94's cost."""
    raise NotImplementedError


def part_a4() -> pd.DataFrame:
    """A4: E114's maps (inputs, PN stage and bouton stage where saved); can be run again as E114 adds files."""
    raise NotImplementedError


def part_a5() -> pd.DataFrame:
    """A5: the bound: the per-odor totals' CV of the ORN input, the ORN recordings, the PN dendrites and each model
    stage."""
    raise NotImplementedError


def choose_level(a2: pd.DataFrame) -> dict:
    """Part B's level by the entry's rule (level.json)."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
