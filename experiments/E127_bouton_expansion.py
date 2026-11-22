#!/usr/bin/env python3
"""E127 step 0: an expansive output at the PN boutons, read on the family's saved stages (no simulation, no olfsysm
change). Stub (interface first); bodies follow.

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u experiments/E127_bouton_expansion.py --part step0
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd

EXPERIMENT = 'E127_bouton_expansion'
QS = tuple(np.round(np.arange(1.0, 3.0 + 1e-9, 0.125), 3))
REF_FACTORS = (0.5, 1.0, 2.0)


def expand(r, r_ref: float, q: float):
    """The expansive output r_ref (max(r, 0) / r_ref)^q."""
    raise NotImplementedError


def r_ref_of(wiring: str) -> float:
    """The wiring's mean PN rate over the panel (E90's saved APL-on PN stage, window mean, all glomeruli, 16 odors)."""
    raise NotImplementedError


def expanded_boutons(pn_win, pn_pre, bouton_win, bouton_pre, gl, r_ref: float, q: float, cond: str):
    """Boutons x stimuli change (window minus pre-odor) after the expansion; APL on: the saved APL inhibition at each
    bouton held (first order), APL off: exact."""
    raise NotImplementedError


def bouton_stats(qs: Sequence[float] = QS) -> pd.DataFrame:
    """E119's per-ROI bouton statistics of every q, wiring, reading and criterion (APL on first order, off exact)."""
    raise NotImplementedError


def apl_ratios(qs: Sequence[float] = QS) -> pd.DataFrame:
    """E98's carrying APL ratios predicted per q, wiring and r_ref factor through E101's claw-to-APL line."""
    raise NotImplementedError


def gate(rows: pd.DataFrame) -> Dict:
    """The step-0 gate of the entry and the q the boutons ask for."""
    raise NotImplementedError


def part_step0() -> Dict:
    raise NotImplementedError
