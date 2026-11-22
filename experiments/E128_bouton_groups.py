#!/usr/bin/env python3
"""E128: the flies' bouton totals by tuning group (no simulations, no olfsysm change).

Stub (interface first); bodies follow.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

EXPERIMENT = 'E128_bouton_groups'


def fly_groups(M: np.ndarray, criterion: float = 2.0) -> np.ndarray:
    """E119 A4's tuning groups of one fly's trial means M (ROIs x odors): group labels (0 the remainder of ROIs in no
    group, i.e. not responsive) per ROI."""
    raise NotImplementedError


def variance_shares(M: np.ndarray, labels: np.ndarray) -> pd.Series:
    """Each group's share of the variance across odors of the summed total: cov(group sum, total) / var(total)."""
    raise NotImplementedError


def readings(M: np.ndarray, labels: np.ndarray, order: Optional[Sequence[int]] = None) -> Dict[str, float]:
    """The totals' CV as measured, with every group at equal weight, and with the 1, 2, 3, 5 groups of largest shares
    left out (the remainder never left out)."""
    raise NotImplementedError


def part0() -> dict:
    """Reproduces E119's a2 statistics and a4 groups per fly; stops on any difference."""
    raise NotImplementedError


def part_a(n_deals: int = 200, seed: int = 0) -> dict:
    """The flies: readings, variance shares, the null (ROIs dealt into groups of the same sizes) and the repeat check."""
    raise NotImplementedError


def part_b(n_deals: int = 1000, seed: int = 0) -> dict:
    """The model: each fly's groups dealt onto the 28 shared glomeruli of E90's saved stages, eight wirings, APL on."""
    raise NotImplementedError


def score() -> dict:
    """The classes and the outcome letter (verdict.json, comparison_table.md)."""
    raise NotImplementedError
