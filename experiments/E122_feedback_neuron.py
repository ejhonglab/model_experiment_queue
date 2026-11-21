#!/usr/bin/env python3
"""E122: one feedback neuron that survives the block in place of KC-KC inhibition (STUB: interface only).

The fitted KC-KC matrix of the family (C2's breadth-2 APL-effect fit on E72's eight wirings) is replaced by a rank-1 (or
per-type block rank-1) matrix with the same total weight, passed through olfsysm's existing KC-KC path (`_wKCKC`; no
olfsysm change). Arms: U (one neuron, all KCs in and out, uniform), TY (one neuron per KC type), G (all KCs in, the gamma
KCs out); controls SH (the hemibrain matrix shuffled within each type pair) and GAB (G's output moved to the alpha-beta
KCs). Parts: 0 (no simulation), 1 (fixed scales on S, both response rates held), score1 (the step-1 gate), 2 (refits of
the passing arms on S, read on T), score.
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd

EXPERIMENT = 'E122_feedback_neuron'
ARMS = ('U', 'TY', 'G', 'SH', 'GAB')


def kc_types_of(kc_index: pd.Index) -> np.ndarray:
    """The KC type ('g', "a'b'", 'ab', ...) of each model KC."""
    raise NotImplementedError


def hemibrain_kckc(kc_ids: Sequence[int]) -> pd.DataFrame:
    """The normalized hemibrain KC x KC matrix (post x pre, mean row sum 1), as fit_mb_model builds it."""
    raise NotImplementedError


def arm_matrix(arm: str, ref: pd.DataFrame, types: np.ndarray, seed: int = 0) -> pd.DataFrame:
    """The arm's KC x KC matrix (post x pre) with the same total weight as `ref`."""
    raise NotImplementedError


def part_0() -> Dict[str, pd.DataFrame]:
    """Step 0, no simulation: the within-type share of the matrix's weight, and how well each KC's received KC-KC input in
    the family's saved runs is predicted by its type's pooled activity."""
    raise NotImplementedError


def part_1(arms: Sequence[str] = ARMS, wirings: Sequence[str] = ('real', 'P1', 'P2', 'M1')) -> pd.DataFrame:
    """Step 1: each arm at the fit's scales on S's wirings, both response rates held (E77_input_screen.hold_rates)."""
    raise NotImplementedError


if __name__ == '__main__':
    raise SystemExit('E122: stub')
