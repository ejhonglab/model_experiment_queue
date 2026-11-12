#!/usr/bin/env python3
"""E86: the reduced searched input as a gated candidate. E77's input searched on the real wiring against one half of the
recordings, at 0.6 of its log size, read on the four held-out wirings (T: P3, P4, M2, M3) against the half its search
never saw: the label-permuted and the inverted controls, the search repeated with a second draw and on a second split of
the recordings, and a verdict fixed before any run.

The entry is item E86 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from hal
2026-11-12 11:08). C2's family throughout (the breadth-2 APL-effect fit of E67), each wiring at its own fit of E72's part
B, the 0.5 s pulse; S: real, P1, P2, M1; T: P3, P4, M2, M3 (E83's division). Halves A and B are E77's
(`recording_halves()`, seed 0); halves C and D are `recording_halves(1)`'s 'A' and 'B'. Every input is read at 0.6 of its
log factors (E83's `input_at`) with both response rates re-held (E83's `input_reading`: the breadth recomputed for the
input, the threshold bisected to the wiring's APL-off response count, then the APL scale, KC->APL and bouton->APL tied,
to its APL-on count), its breadth run cached under outputs/E86_reduced_input_candidate/. For an input, g is T's pattern
part under the calibrated readout, unchanged minus changed, against the half of its split that its search never saw
(positive is better).

Parts (CLI `--part`):
- A (no simulations): the four fitted inputs (E77's two searched inputs, E83 arm I's two) at 0.6 and at full size, by
  the pairs they move (size, projection on what the unseen half asks, the ester / alcohol / other groups' shares of the
  squared pattern error), wiring by wiring as paired differences, the agreement of their factors and of their per-pair
  changes, from E83's readings (pairs.csv, readings.csv); E85's part A4 tables are cited, not recomputed.
- identity: the lab's input reproduces T's saved fits, and E83's part A2 point (P3, half A's input at 0.6) is
  reproduced exactly.
- B1: E83's `permuted_factors` of the candidate's factors, seeds 1 and 2, at 0.6 on T, both directions (16 points).
- B2: the inverted input, the candidate's factors at -0.6 of their log size (8 points).
- B3: seeds 3 and 4, only if a B1 input has g above 0 (16 points).
- C-check: E77's part E search, first split, generator seed 0, the first generation of eb only: must reproduce E77's
  oracle_history_A.csv to the last digit.
- C (`--search C1_A|C1_B|C2_C|C2_D`): E77's part E search as it was, on the real wiring (`input_oracle`), C1 on halves
  A and B with generator seeds 10 and 11, C2 on halves C and D with generator seeds 20 and 21; resumable from a
  checkpoint per generation.
- C-read: each of the four new inputs at 0.6 on T, both response rates re-held (16 points).
- score: readings.csv, pairs.csv, two_by_two.csv, inputs.csv (g per input) from every run so far.
- D: the verdict (gate.csv, replication.csv, null.csv, verdict.json, comparison_table.md).
- E, F: only on "passes" / "carries" (not written).

Outputs: results/E86_reduced_input_candidate/ (shared directory) and outputs/E86_reduced_input_candidate/ (local: spike
counts and responses per run under <part>/, the searches' checkpoints under C/, the breadth runs).

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E86_reduced_input_candidate.py --part B
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E86_reduced_input_candidate.py --part C --search C1_A
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E86_reduced_input_candidate'
SIZE = 0.6
INVERT_SIZE = -0.6
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
B1_SEEDS = (1, 2)
B3_SEEDS = (3, 4)
# the halves of each split of the recordings: split seed -> {half name: the half of recording_halves(seed)}
SPLITS = {0: {'A': 'A', 'B': 'B'}, 1: {'C': 'A', 'D': 'B'}}
UNSEEN = {'A': 'B', 'B': 'A', 'C': 'D', 'D': 'C'}
# part C's searches: name -> the split, the half searched (its name in this experiment) and input_oracle's `seed`
# (input_oracle seeds its generator with seed for the split's first half and seed + 1 for its second)
SEARCHES = {'C1_A': dict(split_seed=0, half='A', seed=10), 'C1_B': dict(split_seed=0, half='B', seed=10),
            'C2_C': dict(split_seed=1, half='C', seed=20), 'C2_D': dict(split_seed=1, half='D', seed=20)}
GENERATOR_SEEDS = {'C1_A': 10, 'C1_B': 11, 'C2_C': 20, 'C2_D': 21}
# the verdict's thresholds (the entry's part D, fixed before any run)
CONTROL_SHARE, B3_MEAN_SHARE = 0.5, 0.25
GATE_GAIN, GATE_SINGLES = 0.003, 12
GATE_RAW_SLACK, GATE_MEAN_SLACK, GATE_SQ_ALL, GATE_SQ_SIX = 0.005, 0.01, 0.15, 0.075
PASS_MEAN, CARRY_MEAN, MIN_POSITIVE = 0.002, 0.001, 3


def halves() -> Dict[str, Dict[str, list]]:
    """The four halves of the recordings: {'A', 'B', 'C', 'D'}: {'ctrl' / 'tnt': [recording keys]}."""
    raise NotImplementedError


def targets() -> Dict[str, pd.Series]:
    """E83's targets (halfA, halfB, z8, z6, noise_free, noise8) and halfC, halfD (canonical pair order)."""
    raise NotImplementedError


def candidate_factors(half: str) -> pd.DataFrame:
    """E77's searched factors on half `half` ('A' or 'B'), glomeruli x odors (E83's `read_factors`)."""
    raise NotImplementedError


def control_factors(kind: str, half: str, seed: Optional[int] = None) -> Tuple[pd.DataFrame, float]:
    """A control input's factors and size: ('perm', seed) E83's `permuted_factors` of the candidate's factors at SIZE;
    ('inverted', None) the candidate's factors at INVERT_SIZE."""
    raise NotImplementedError


def held_reading(label: str, orn: pd.DataFrame, input_tag: str) -> dict:
    """E83's `input_reading` with this experiment's breadth cache (outputs/E86_reduced_input_candidate/)."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """Part A (no simulations): see the module docstring."""
    raise NotImplementedError


def identity(out: Path) -> dict:
    """The identity checks: T's simulators reproduce their saved fits; E83's A2 point (P3, half A at 0.6) exactly."""
    raise NotImplementedError


def plan_b(part: str) -> List[dict]:
    """The points of part B1, B2 or B3: [dict(change, half, unseen, kind, seed, size, wiring)]."""
    raise NotImplementedError


def run_points(out: Path, part: str, plan: Sequence[dict]) -> pd.DataFrame:
    """Held readings of a plan's points (resumable; runs_<part>.csv, provenance_<part>.json, outputs under
    outputs/E86_reduced_input_candidate/<part>/)."""
    raise NotImplementedError


def check_search(out: Path) -> dict:
    """Part C's check: input_oracle with the first split, generator seed 0, one generation of eb, against E77's
    oracle_history_A.csv to the last digit (check_search.json)."""
    raise NotImplementedError


def search(out: Path, name: str, resume: bool = True) -> pd.DataFrame:
    """One of part C's searches (SEARCHES), resumable from its checkpoint."""
    raise NotImplementedError


def score(out: Path) -> pd.DataFrame:
    """readings.csv, pairs.csv, two_by_two.csv and inputs.csv (one row per input: g and the gate's statistics)."""
    raise NotImplementedError


def verdict(inputs: pd.DataFrame, flies_mean: float) -> dict:
    """The entry's verdict (part D) from inputs.csv's table (see the module docstring)."""
    raise NotImplementedError


def write_verdict(out: Path) -> dict:
    """gate.csv, replication.csv, null.csv, verdict.json and the verdict's section of comparison_table.md."""
    raise NotImplementedError


def part_e(out: Path) -> None:
    """Part E, only on "passes" or "carries" (not written)."""
    raise NotImplementedError('part E runs only on "passes" or "carries"; not written yet')


def part_f(out: Path) -> None:
    """Part F, only on "passes" (not written; needs E83's build_simulator with the family and the pulse as arguments)."""
    raise NotImplementedError('part F runs only on "passes"; not written yet')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'identity', 'B1', 'B2', 'B3', 'B', 'C-check', 'C', 'C-read',
                                                      'score', 'D', 'E', 'F'])
    ap.add_argument('--search', choices=list(SEARCHES), help='part C: the search to run')
    ap.add_argument('--no-resume', action='store_true', help='part C: start the search over (ignore its checkpoint)')
    args = ap.parse_args()
    raise NotImplementedError(f'part {args.part}: stub')


if __name__ == '__main__':
    main()
