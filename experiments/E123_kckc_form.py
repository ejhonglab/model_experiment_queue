#!/usr/bin/env python3
"""E123: the form of the inhibition among KCs: draws of the shuffle, random sparse matrices and the inhibition as a part
of APL (stub; the bodies follow).

Builds on E122 (experiments/E122_feedback_neuron.py) and E110 (experiments/E110_no_kckc.py, part A).
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

EXPERIMENT = 'E123_kckc_form'


def amount_kept_shuffle(ref: np.ndarray, types: np.ndarray, seed: int = 0) -> np.ndarray:
    """Each KC's row permuted within each presynaptic type (the diagonal kept 0): row and per-type amounts kept."""
    raise NotImplementedError


def random_k_partners(n: int, k: int, total: float, seed: int = 0) -> np.ndarray:
    """k random partners per KC (no self), equal weights, every row sum total / n."""
    raise NotImplementedError


def arm_array(arm: str, ref: np.ndarray, types: np.ndarray) -> np.ndarray:
    raise NotImplementedError


def part_check() -> dict:
    raise NotImplementedError


def part_1(arms: Sequence[str], wirings: Sequence[str], lane: Optional[str] = None) -> pd.DataFrame:
    raise NotImplementedError


def score() -> dict:
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    raise NotImplementedError
