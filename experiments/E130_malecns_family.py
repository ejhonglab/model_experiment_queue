#!/usr/bin/env python3
"""E130: the best family on the male CNS's two calyces (a second fly's wiring, APL synapses and claw-to-claw synapses
among KCs in place of the hemibrain's).

STUB (interface first): signatures and docstrings; the bodies follow.

"The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's
`saved_fit`). The calyces: the hemibrain with multibouton_claws 'argmax' (olfsysm's default since 60cddc3) and the male
CNS's left and right calyces (olfsysm's pn2kc_connections 'malecns-left' / 'malecns-right'), through
E10_joint_fit.py's `--connectome`.

Parts (CLI `--part`): step0, check, 1, fit (`--targets`), score, variants, commands.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E130_malecns_family'
SCRIPT = 'experiments/E130_malecns_family.py'
CALYCES = ('hemibrain-argmax', 'malecns-left', 'malecns-right')
MALE = ('malecns-left', 'malecns-right')
TARGETS = ('malecns-left', 'malecns-left-kckc0', 'hemibrain-argmax', 'malecns-right', 'malecns-right-kckc0')


def malecns_kckc(connectome: str, kc_ids: Sequence[int]) -> pd.DataFrame:
    """The male calyx's KC x KC matrix (post x pre) among the model's KCs `kc_ids`, from Pratyush's KCclaw-to-KCclaw
    table, normalized as olfsysm's `connectome_KCKC_weights` (mean row sum 1)."""
    raise NotImplementedError


def step0() -> Dict[str, pd.DataFrame]:
    """The three calyces side by side (no simulation)."""
    raise NotImplementedError


def part_check() -> dict:
    """The identities: 'first' reproduces the family's saved hemibrain run; the male claw counts are olfsysm's README's."""
    raise NotImplementedError


def part_1() -> pd.DataFrame:
    """Step 1: six runs at fixed scales (the three calyces, KC-KC at the fit's scale and at 0), both rates re-held."""
    raise NotImplementedError


def run_fits(targets: Sequence[str], lane: Optional[str] = None, reclaim: bool = False) -> List[str]:
    """Step 2: the fits of `targets` (E72's part-B command), claimed by mkdir, resumable."""
    raise NotImplementedError


def score() -> dict:
    """The readings and the classes (no simulation)."""
    raise NotImplementedError


def write_variants() -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv for the model-output archive."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'check', '1', 'fit', 'score', 'variants', 'commands'))
    ap.add_argument('--targets', nargs='+', default=list(TARGETS), choices=list(TARGETS))
    a = ap.parse_args()
    raise NotImplementedError(a.part)


if __name__ == '__main__':
    main()
