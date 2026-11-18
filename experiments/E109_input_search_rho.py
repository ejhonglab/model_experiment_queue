#!/usr/bin/env python3
"""E109: the input of paa, B-myr and 2-but searched for the two ρ's on four wirings and one half of the recordings, read
on the other four wirings and the other half.

STUB (interface first): signatures and docstrings; the bodies raise NotImplementedError until the implementation lands.

The entry is item E109 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2; moved into the list 2026-11-17 23:56 under the idle-host rule; claimed by viki 23:56).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E109_input_search_rho'
SCRIPT = 'experiments/E109_input_search_rho.py'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
SEARCH_ODORS = ('paa', 'B-myr', '2-but')         # E64's order
CONTROL_ODORS = ('2h', '6al')                    # E64's controls
GATE_SIZE, FULL_SIZE = E83.GATE_SIZE, 1.0
INVERT_SIZE = -0.6                               # E86's
PERM_SEEDS = (1, 2, 3, 4)                        # E86's B1 and B3 seeds of E83's `permuted_factors`
PATTERN_MARGIN = 0.0015


class StartRefused(RuntimeError):
    """A search of the three odors asked to start before E108's class is written (names what it waits for)."""


def e108_class(results: Optional[Path] = None) -> dict:
    """E108's class as its results give it: 'written', 'class', 'carrying_rule', 'source', 'waiting_for'."""
    raise NotImplementedError


def require_e108_class(name: str, confirmed_in: Optional[str] = None) -> dict:
    """The start rule: raises StartRefused for a search of SEARCH_ODORS until E108's class is written."""
    raise NotImplementedError


def start_input(e108: Optional[dict] = None) -> Tuple[pd.DataFrame, dict]:
    """The input the searches start from: the lab's, or (E108 'better') the lab's with E108's carrying paa column."""
    raise NotImplementedError


def two_rhos_population(on: Mapping[str, pd.Series], off: Mapping[str, pd.Series], data_on: pd.Series,
                        data_off: pd.Series) -> Dict[str, float]:
    """The two ρ's of the mean matrices over the given wirings' pair distances against one half's KC matrices."""
    raise NotImplementedError


def pattern_penalty(rms: float, rms0: float, margin: float = PATTERN_MARGIN) -> float:
    """The penalty on a pattern part more than `margin` above the unchanged fits' (`rms0`)."""
    raise NotImplementedError


def search(name: str, resume: bool = True) -> pd.DataFrame:
    """One of the searches (SEARCHES): E83 part C's search with the two ρ's as its objective."""
    raise NotImplementedError


def part0() -> dict:
    """Part 0: the identities, the bound and the references."""
    raise NotImplementedError


def gate(rd: pd.DataFrame, runs: pd.DataFrame, ref: str = 'I0') -> dict:
    """The entry's gate on the held readings."""
    raise NotImplementedError


def class_of(g: Mapping, with_pending: bool = True) -> str:
    """'passes', 'carries below the gate', 'a price' or 'does not carry'."""
    raise NotImplementedError


def recurring_entries() -> pd.DataFrame:
    """Step 6's table of the glomerulus x odor entries that recur across the six searches."""
    raise NotImplementedError


def step7() -> dict:
    """Step 7, beside with no gate: paa against B-myr, and Remy's validation2 flies."""
    raise NotImplementedError


def main() -> None:
    raise NotImplementedError


if __name__ == '__main__':
    main()
