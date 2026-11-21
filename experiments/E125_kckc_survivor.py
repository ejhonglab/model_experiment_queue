#!/usr/bin/env python3
"""E125: the part of the inhibition among KCs that outlasts the block: a pooled neuron or a KC's own feedback in place
of the matrix.

STUB (interface first): signatures and docstrings; the bodies follow.

The entry is item E125 of experiment_priorities.md (proposed by viki 2026-11-21 14:27, moved into the list 14:33 under
the idle-host rule; planned by Fable subagents, claude-fable-5-2). It builds on E124 (experiments/E124_kckc_apl_borne.py:
its ladder, the hemibrain matrix's KC-KC scale with APL off at 0.25, 0.5 and 0.75 of the fit's) and E123 / E122 (the
runner at fixed scales, the classes, the pooled and per-type matrices).

Parts (CLI `--part`): check, arm1, arm2, arm3, score, variants, chain.
"""
from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence

import numpy as np

EXPERIMENT = 'E125_kckc_survivor'
SCRIPT = 'experiments/E125_kckc_survivor.py'
DOSES = (0.25, 0.5, 0.75)


def survivor_array(form: str, ref: np.ndarray, types: np.ndarray) -> np.ndarray:
    """The survivor X of a form ('PO' pooled, 'DG' a KC's own feedback, 'TY' one neuron per type) with the hemibrain
    matrix's total weight (times d by the caller)."""
    raise NotImplementedError


def condition_arrays(arm: str, ref: np.ndarray, types: np.ndarray) -> Dict[str, np.ndarray]:
    """{'on': (1 - d) ref + d X, 'off': d X} of an arm."""
    raise NotImplementedError


def two_condition_simulator(w: str, params: Mapping[str, float], M_on, M_off):
    """A simulate function that sends the APL-on condition to a simulator built with `M_on` and the APL-off condition
    to one built with `M_off`."""
    raise NotImplementedError


def part_check() -> dict:
    """d = 0 reproduces E123's S1_HBab_real; a diagonal matrix alone changes a KC's spike counts. check.json."""
    raise NotImplementedError


def part_runs(arms: Sequence[str], wirings: Sequence[str], lane: Optional[str] = None):
    """Arms 1 to 3 at fixed scales, both response rates re-held (resumable)."""
    raise NotImplementedError


def score() -> dict:
    """readings.csv, runs.csv, verdict.json, comparison_table.md."""
    raise NotImplementedError


def write_variants():
    """variants.csv and model_pair_distances_long.csv for the archive."""
    raise NotImplementedError


def chain(wirings: Sequence[str]) -> None:
    """check, arms 1 to 3, score, variants."""
    raise NotImplementedError


if __name__ == '__main__':
    raise SystemExit('E125: stub')
