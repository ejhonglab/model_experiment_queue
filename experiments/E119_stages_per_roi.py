#!/usr/bin/env python3
"""E119: the model's boutons, claws and claw-to-soma step against Remy's per-ROI megamat recordings (no simulations).

The entry is item E119 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-20 13:32). Part 0 is data_prep/remy_per_roi_data.py (its outputs in <shared>/data/remy_per_roi/). The model
side is E90's archive (data/model_outputs/E90_model_stages/B/: C2's family on E72 part B's eight wirings, APL on and
acutely off) at the settled pre-odor level (E90_settled_preodor's preodor_levels.parquet), the 12 megamat odors of the
panel; E90's readings of a stage (`stage_units`, `claw_signal`, `session_z`), its groups and its rule (`differs`,
`element_result`, `group_verdict`), and E81's noise removal (`H.noise_free_pair_distances`, Pearson on centred trials as
E92's `noise_removed`).

Parts (CLI `--part`):
- A, the boutons (6 flies): A1 per-odor totals; A2 per-ROI response rates, odors per ROI, amplitude spread; A3 cosine
  and Pearson distances, noise removed; A4 (a lead) bouton groups against the model's glomeruli.
- B, the claws (4 flies): B1 E90's nine claw readings scored on per-ROI statistics (kept / not); B2 E90's group
  verdicts under the kept readings; B3 the bouton-to-claw change.
- C (not run until part 0b is cleared): the same fly's claw-to-soma step.
- table: comparison_table.md.

Outputs: <shared>/results/E119_stages_per_roi/; large tables under outputs/E119_stages_per_roi/.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E90_model_stages as E90
from kcq import harness as H

EXPERIMENT = 'E119_stages_per_roi'
ODORS = E90.ODORS
CRITERIA = (1.0, 2.0, 3.0)
VERDICT_CRITERION = 2.0
CLAW_POOLS = (1, 7, 14, 28)
BOUTON_POOLS = (1, 3)
VERDICT_POOL = {'bouton': 1, 'claw': 14}


def fly_tables(stage: str) -> dict:
    raise NotImplementedError


def roi_stats(V: pd.DataFrame, thr: float) -> dict:
    raise NotImplementedError


def model_rois(stage: str, wiring: str, cond: str, reading: str, share: float, pool: int, seed: int) -> pd.DataFrame:
    raise NotImplementedError


def part_a() -> None:
    raise NotImplementedError


def part_b() -> None:
    raise NotImplementedError


def part_c() -> None:
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
