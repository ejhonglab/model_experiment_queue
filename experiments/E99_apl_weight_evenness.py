#!/usr/bin/env python3
"""E99: the evenness of APL's weights among KCs and over a KC's claws: what loses the mean APL effect when APL's counted
synapses take the bouton path's place at the KCs (the weaker growth of a KC's amount with its number of claws, the
inequality among KCs with the same number of claws, or the uneven placement over a KC's claws, each changed alone at fixed
scales), and whether each KC's count of all of APL's synapses, spread evenly over its claws, holds the mean.

STUB (interface first): the signatures and docstrings below are the module's interface; the bodies follow.

The entry is item E99 of experiment_priorities.md (proposed by viki 2026-11-15 12:48; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-15 12:27). Built on E96 (`experiments/E96_offclaw_apl.py`: the counts, the
components, the weights, the runner, the readings, the choice, the gate, the equal rule, part D), E95 (`StaticModel`,
`kc_losses`, `d_unit`, `d_realism`), E94 (the stage readers) and E92 (the bracket).
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E99_apl_weight_evenness'

# the new settings (the entry's Definitions): each a per-claw APL->KC weight component with mean 1 over the claws
SHAPES = ('K1', 'K2', 'K3', 'K4', 'W1', 'W2', 'T1', 'T2')
GROUP_MIN = 20                                      # a claw-count group with fewer KCs is pooled with the next lower count


def claw_groups(n_claws: np.ndarray, min_kcs: int = GROUP_MIN) -> np.ndarray:
    """Each KC's claw-count group G (labelled by the group's lowest claw count): the KCs with the same number of claws,
    a group with fewer than `min_kcs` KCs pooled with the next lower count."""
    raise NotImplementedError


def shape99(name: str, kc: pd.DataFrame, claws: pd.DataFrame, counts: str = 'offclaw', frac: float = 0.5) -> np.ndarray:
    """The component per claw (mean 1 over the claws) of a setting `name` (SHAPES) from a wiring's per-KC and per-claw
    tables (E96's part 0 cache): K1 to K4 on `counts` ('offclaw' the calyx's off-claw counts a^c_k, 'offclaw_ped' with the
    peduncle's), K3 with `frac` of the inequality kept."""
    raise NotImplementedError


def variance_split(x: np.ndarray, kc_of_claw: np.ndarray, group_of_claw: np.ndarray) -> dict:
    """The variance of a per-claw component (mean 1) split into the part between claw-count groups, among the KCs of a
    group and within KCs (claw-weighted; they add up to the squared CV per claw)."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    """A variant's fields (this experiment's settings, E96's and E95's saved settings, the reference, D0)."""
    raise NotImplementedError


def part_0(out: Path, wirings: Sequence[str]) -> dict:
    """The identities (four runs on the real wiring, and the shapes' identities without simulations)."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """A1 the shapes, A2 the predictions (written before part B's first point), A3 who responds, A4 the bound's
    yardstick (no simulations)."""
    raise NotImplementedError


def part_b(out: Path, blocks: Optional[Sequence[str]] = None, wirings: Sequence[str] = ()) -> pd.DataFrame:
    """Part B's points (blocks I and II) on S, resumable (runs_B.csv)."""
    raise NotImplementedError


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """The readings, part C's verdicts (the shares and classes, the bound, presence), the choice, gate and equal rule,
    part D's rule, the table."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'score', 'C', 'Dcheck', 'D', 'Dscore', 'E', 'steps',
                                                      'variants', 'table'))
    ap.parse_args()
    raise NotImplementedError('E99: stub')


if __name__ == '__main__':
    main()
