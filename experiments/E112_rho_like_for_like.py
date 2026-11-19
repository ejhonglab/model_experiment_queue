#!/usr/bin/env python3
"""E112: the two ρ's on odor-odor distances read as the recordings are made: the best family's saved runs through the
calibrated readout, the per-KC session z-score and the trial noise; the gap, its bound and the gates of the input
experiments read again (no new simulations in part A).

The entry is item E112 of experiment_priorities.md (proposed by viki 2026-11-19 13:58; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 13:50; claimed by hal 14:01). STUB (interface first): the bodies raise
NotImplementedError until the implementation is committed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E112_rho_like_for_like'
READINGS = ('R0', 'R1', 'R2', 'R3', 'R4')
N_PSEUDO = {'on': 11, 'off': 8}
N_TRIALS = 3


def reading_distances(o: Mapping[str, tuple], reading: str, seed: int = 0) -> Dict[str, pd.Series]:
    """{'on', 'off'} pair distances (canonical order) of a saved run `o` ({'on' / 'off': (responses, spike counts)})
    under reading R0 (raw counts), R1 (E65's calibrated readout, K on the run's APL-on counts), R2 (E73's per-KC session
    z-score, saturating map, graded, in expectation over the measurement noise) or R3 (R2 through sampled measurement
    noise on pseudo-recordings: 11 with APL on and 8 with APL off, 3 trials each, the distance averaged over them)."""
    raise NotImplementedError


def pseudo_recording_distances(z: Mapping[str, object], n_recordings: int, n_trials: int = N_TRIALS,
                               seed: int = 0) -> pd.Series:
    """The mean over `n_recordings` pseudo-recordings of the cosine distances of their trial means: a `session_zscore`
    reading's noise-free values plus independent Gaussian measurement noise of each KC's trial sd, per trial."""
    raise NotImplementedError


def trial_noise_distances(pairs_long: pd.DataFrame) -> Dict[str, Dict[str, pd.Series]]:
    """R4 from E43's pipeline table (partB/pairs_noise_long.csv): {model: {'on', 'off'}} the mean over the
    pseudo-recordings of the distances of their three-trial means (graded spike counts), APL on at the level matched to
    the control flies and APL off at the level matched to APL output blocked (TNT)."""
    raise NotImplementedError


def reading_scores(d: Mapping[str, pd.Series]) -> dict:
    """The scores of one unit's distances under a reading: the two ρ's over the 120, 66, 54 and 91 pairs, the distance
    of paa from B-myr, the APL-effect loss with its mean, pattern, amplitude and shape parts, ρ and Pearson r on the
    per-pair APL effect, the mean distances."""
    raise NotImplementedError


def bound_scores(d: Mapping[str, pd.Series]) -> dict:
    """E106's bound under a reading: the two ρ's with the pairs of paa and B-myr (and of paa, B-myr and 2-but) at the
    flies' values (raw and quantile), and the odor part and pair part hybrids."""
    raise NotImplementedError


def moves_verdict(rows: pd.DataFrame) -> dict:
    """The entry's "the reading moves ρ" verdict from the per-wiring and population ρ's under R0, R3 and R4."""
    raise NotImplementedError


def reread_gate(experiment: str, reading: str) -> dict:
    """A saved gate (E109's 184 held readings, E108's 156 runs, E97's 84 points) re-scored under `reading`."""
    raise NotImplementedError


def e109_holds(gates: Mapping[str, dict]) -> dict:
    """The entry's "E109's gain holds like for like" verdict per reading."""
    raise NotImplementedError


def part_a() -> None:
    """Part A (no simulations)."""
    raise NotImplementedError


def part_b_prepare() -> pd.DataFrame:
    """Part B's model table (fits.csv rows for E43's pipeline) for E109's inputs on the four held-out wirings."""
    raise NotImplementedError


def part_b_score() -> None:
    """Part B's scores (R4 on all eight wirings and E109's inputs)."""
    raise NotImplementedError
