#!/usr/bin/env python3
"""E96: APL's synapses on the KCs outside the claws, in the calyx and in the peduncle, which the model's weights leave out,
as a component of the APL->KC weights per KC in the feedback APL with the bouton path at a tenth; whether weights that the
connectome has do what E95's uniform component did; and the family fitted again that way.

STUB (interface first): the signatures and docstrings below are the interface; the bodies are being written.

The entry is item E96 of experiment_priorities.md (proposed by viki 2026-11-14 10:21; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 10:02). Built on E95 (`experiments/E95_static_apl.py`, its part D and F).

Parts (CLI `--part`): 0 (the identities), A (A1 to A4, no simulations), B (`--arm` B1 to B6, 64 points on S), score, C
(the choice per half and the gate on T), Dcheck (E10_joint_fit.py's new --aplkc-flat modes at size 0 against E94's D1
fit's first evaluation), D (the family fitted again, two starts per S wiring: `--start own` / `--start e95`), Dscore, E
(stubs), steps, variants, table.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E96_offclaw_apl'

# the settings (the entry's Design; fixed before any run)
D_DOSE = 0.1                                        # the bouton path at a tenth (E94's anatomical dose)
B_SIZES = (0.5, 1.0, 2.0)                           # the component's size d in E95's unit
PERM_SEED = 1
MODES = ('offclaw', 'offclaw_ped', 'allsyn')        # B1, B2 (B3 to B5 too), B6


def synapse_counts(claw_index: pd.MultiIndex) -> Tuple[pd.DataFrame, pd.Series]:
    """Per KC of the wiring (and per kept claw) the APL->KC synapses of the table the model reads, with the loader's
    filters: on kept claws, on claws the wiring drops, on no claw in the calyx (a^c_k), in the peduncle (a^p_k)."""
    raise NotImplementedError


def component_shape(mode: str, counts: pd.DataFrame, kept_per_claw: pd.Series, claw_index: pd.MultiIndex,
                    permute: Optional[str] = None, seed: int = PERM_SEED) -> np.ndarray:
    """The component per claw (mean 1 over the claws): a_k / n_k on each of KC k's claws ('offclaw': a^c_k,
    'offclaw_ped': a^c_k + a^p_k; 'allsyn': the kept claw's own count plus the KC's off-claw part, in the proportion in
    which they are counted; 'uniform' 1 everywhere), the counts permuted within a claw count or over all KCs."""
    raise NotImplementedError


def e96_weights(w_fam: pd.Series, wnorm: pd.Series, shape: np.ndarray, d: float, unit: float) -> pd.Series:
    """The APL->KC weights per claw: the family's plus d x unit x mean(wnorm) x shape (E95's `d_weights` for the uniform
    shape); d = 0 the family's exactly."""
    raise NotImplementedError


def parse_variant(v: str) -> dict:
    """A variant's fields (arm, mode, size, pool, dose, permute, seed, control)."""
    raise NotImplementedError


def part_0(out: Path, wirings: Sequence[str]) -> dict:
    """The identities of part (0) (identity.json)."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """A1 to A4 (no simulations)."""
    raise NotImplementedError


def part_points(out: Path, arms: Optional[Sequence[str]], wirings: Sequence[str]) -> pd.DataFrame:
    """Part B's points (resumable; runs_B.csv)."""
    raise NotImplementedError


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """readings.csv and the tables of part C."""
    raise NotImplementedError


def choose96(rd: pd.DataFrame, half: str) -> dict:
    """The choice per half (the entry's (C))."""
    raise NotImplementedError


def d_condition(rd: pd.DataFrame) -> dict:
    """Whether part D runs, which component, the pooled share, B6."""
    raise NotImplementedError


def d_command(label: str, start: str, mode: str, pool: Optional[float]) -> List[str]:
    """One part D fit's E10_joint_fit.py command."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'A1', 'A2', 'A3', 'A4', 'B', 'score', 'C', 'Dcheck', 'D',
                                                      'Dscore', 'E', 'steps', 'variants', 'table'))
    ap.add_argument('--arm', default=None)
    ap.add_argument('--wirings', default=None)
    ap.add_argument('--pop', default='S', choices=('S', 'T'))
    ap.add_argument('--start', default='own', choices=('own', 'e95'))
    ap.add_argument('--dry-run', action='store_true')
    ap.parse_args()
    raise NotImplementedError('E96 is a stub so far')


if __name__ == '__main__':
    main()
