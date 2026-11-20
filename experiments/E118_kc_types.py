#!/usr/bin/env python3
"""E118: the KC types in the best family and in the recordings: each type's response rate, the two ρ's and the per-pair
APL effect by type and by the sample's type composition (no simulations), then each type's response rate held at
measured ratios (parts B and C, on the entry's rule only).

The entry is item E118 of experiment_priorities.md (proposed by viki 2026-11-19 19:16; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-19 19:07). "The family" is the breadth-2 APL-effect fit (C2, E67) on its eight
wirings (the hemibrain wiring and seven re-drawn within type; E72's part B, `E77_input_side_common_error.members`), read
per wiring and as the population (the mean over wirings of the pair distances per condition, E72's). "The two ρ's" are
ρ on odor-odor distances with APL on against the control flies and with APL off against APL output blocked (TNT) (E106's
`two_rhos`; the KC matrices the mean over recordings of the cosine distance), the model's the cosine distance on raw
spike counts. The KC types are the hemibrain's labels of the model's KCs (αβ 'ab', α′β′ "a'b'", γ 'g').

Parts (CLI `--part`):
- 0, identities: the all-KC reading of each wiring's saved counts reproduces E106's pair distances (raw counts) and its
  population ρ's (0.690 / 0.756 over the 120 pairs); the types' KCs and spikes add up to each wiring's; the weighted
  reading with every weight 1 is the plain reading; E108's 'all' targets are E106's data matrices.
- A (no simulations): A1 by type (each type's KCs alone: response rates and their ratio, share of spikes, the two ρ's
  over the 120, 66, 54 and 91 pairs and against halves A and B of E77's split, the distance of paa from B-myr, the
  APL-effect statistics on raw counts and under E65's calibrated readout with K set on all KCs); A2 the composition (γ
  KCs weighted 0.5, 1, 2, 4, 8 and α′β′ KCs 0 or 1, as if the recorded sample held that share; a weight w multiplies a
  KC's row by sqrt(w), so that cosine similarities count it w times), also against Remy's validation2 somata on E107's
  as-fitted runs (a description); A3 controls (the same weights with the type labels permuted among KCs of equal claw
  count, 20 permutations; random subsets of each type's size; each type's distance of paa from B-myr against its
  response rate over the 48 type x condition x wiring points); A4 what is measured (published per-type response rates,
  only as on record in the queue's notes; E80's saved per-type points read for the two ρ's from E106's archive table);
  the rule for part B and the outcome letter as far as part A decides it.
- B, C: the per-type threshold factors (E80_screen's `type_factors` through E77's `hold_rates`) and the refit; stubs
  only until E114's runs are done.

Outputs: <shared>/results/E118_kc_types/; nothing large (part A reads archives only).
Usage (from the repository root):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E118_kc_types.py --part 0
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E118_kc_types.py --part A
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E118_kc_types'
TYPES = ('ab', "a'b'", 'g')
TYPE_LABEL = {'ab': 'αβ', "a'b'": 'α′β′', 'g': 'γ', 'all': 'all KCs'}
G_WEIGHTS = (0.5, 1.0, 2.0, 4.0, 8.0)
APBP_WEIGHTS = (0.0, 1.0)
N_PERM = 20
N_SUBSETS = 20
HALVES = ('A', 'B')
PAA_BMYR = ('paa', 'B-myr')
# the rule for part B (fixed before any reading is scored)
RULE_GAIN, RULE_OTHER, RULE_MAX_W, RULE_WIRINGS, RULE_CONTROL = 0.010, 0.005, 4.0, 6, 0.5


def out_dir() -> Path:
    """<shared>/results/E118_kc_types/ (created)."""
    raise NotImplementedError


def members() -> List[Tuple[str, str, str]]:
    """The family's eight wirings, [(label, folder, tag)] (E77's `members`)."""
    raise NotImplementedError


def load_family() -> Dict[str, Dict[str, tuple]]:
    """{wiring label: {'on', 'off'}: (responses, spike counts)} of the eight saved fits (E72's `load_fit`)."""
    raise NotImplementedError


def claw_counts(index: pd.Index) -> pd.Series:
    """Claws per KC (the hemibrain bouton base's wiring, archived inputs/wPNKC_bouton_base.parquet), in `index`'s
    order. The re-draws keep each KC's claw count, so one count serves all eight wirings."""
    raise NotImplementedError


def kc_types(index: pd.Index) -> np.ndarray:
    """The KC type label of each row of `index` (level 'kc_type')."""
    raise NotImplementedError


def kc_weights(labels: np.ndarray, g: float = 1.0, apbp: float = 1.0, ab: float = 1.0) -> np.ndarray:
    """Per-KC composition weights: `g` for the γ KCs, `apbp` for the α′β′ KCs, `ab` for the αβ KCs."""
    raise NotImplementedError


def permuted_labels(labels: np.ndarray, claws: np.ndarray, seed: int) -> np.ndarray:
    """The type labels permuted among KCs of equal claw count (each claw count's labels shuffled among its KCs)."""
    raise NotImplementedError


def weighted_pair_distances(counts: pd.DataFrame, weights: np.ndarray, kind: str = 'raw', param=None) -> pd.Series:
    """Cosine odor-odor distances after E65's readout (`kind`, `param`) with each KC's row multiplied by sqrt(weight)
    (KCs of weight 0 dropped); weights all 1 give `kcq.harness.readout_pair_distances`."""
    raise NotImplementedError


def weighted_rate(responses: pd.DataFrame, weights: np.ndarray) -> float:
    """The response rate of the weighted sample: sum over KCs of weight x fraction of odors responded to, over the
    summed weight."""
    raise NotImplementedError


def wiring_K(o: Mapping[str, tuple]) -> float:
    """E65's calibrated readout's K of one fit, set on all its KCs' APL-on counts (E72's `distances` 'cal')."""
    raise NotImplementedError


def distances_for(o: Mapping[str, tuple], weights: np.ndarray, K: float) -> Dict[str, Dict[str, pd.Series]]:
    """{'raw', 'cal'}: {'on', 'off'} weighted pair distances of one fit."""
    raise NotImplementedError


def targets() -> dict:
    """The KC data: 'dist' {'all', 'A', 'B'}: {'on', 'off'} mean matrices (all recordings E106's `data_matrices`, halves
    E108's `dist_targets`); 'effect' {'all', 'A', 'B'}: the per-pair APL effect."""
    raise NotImplementedError


def score_distances(d: Mapping[str, Mapping[str, pd.Series]], T: Mapping) -> Dict[str, float]:
    """The readings of one set of distances: the two ρ's over the 120, 66, 54 and 91 pairs (all recordings) and over the
    120 against halves A and B; the distance of paa from B-myr per condition; the mean distance per condition; the
    APL-effect statistics (loss, mean, pattern, amplitude, shape, spread ratio, ρ and Pearson r on the per-pair APL
    effect) on raw counts and under the calibrated readout against all recordings, and the loss, mean and pattern part
    against each half."""
    raise NotImplementedError


def reading(fam: Mapping[str, Mapping[str, tuple]], weights: Mapping[str, np.ndarray], Ks: Mapping[str, float],
            T: Mapping) -> pd.DataFrame:
    """Rows per wiring and 'population' (the mean over wirings of the distances): response rates on / off and their
    ratio, the share of spikes per condition, and `score_distances`."""
    raise NotImplementedError


def part0(out: Path) -> dict:
    """Part 0's identities (part0_identity.json); raises on a failure."""
    raise NotImplementedError


def type_anatomy(fam: Mapping[str, Mapping[str, tuple]]) -> pd.DataFrame:
    """Per type and all KCs: the number of KCs, the mean claws per KC, and on the real wiring the mean APL→KC weight
    (the connectome's, normalized, and the family's breadth-weighted ones)."""
    raise NotImplementedError


def part_a1(out: Path, fam, T, Ks) -> pd.DataFrame:
    """A1 (by_type.csv)."""
    raise NotImplementedError


def part_a2(out: Path, fam, T, Ks) -> pd.DataFrame:
    """A2 (composition.csv)."""
    raise NotImplementedError


def validation2_composition(out: Path) -> pd.DataFrame:
    """A2 against Remy's validation2 somata (E107's as-fitted runs, ρ over her carrying set of pairs, APL on; a
    description; composition_validation2.csv)."""
    raise NotImplementedError


def part_a3(out: Path, fam, T, Ks) -> Dict[str, pd.DataFrame]:
    """A3: the permuted-label control (composition_permuted.csv), random subsets (random_subsets.csv), the sparseness
    points and their fit (sparseness_points.csv, sparseness_fit.csv)."""
    raise NotImplementedError


def published_type_rates() -> pd.DataFrame:
    """Published per-type response rates as on record in the queue's notes (type, rate, ratio to αβ, source); empty
    when none is on record."""
    raise NotImplementedError


def part_a4(out: Path, a1: pd.DataFrame) -> pd.DataFrame:
    """A4 (measured.csv, e80_type_points.csv)."""
    raise NotImplementedError


def rule_b(comp: pd.DataFrame, perm: pd.DataFrame, a1: pd.DataFrame, published: pd.DataFrame) -> dict:
    """The rule for part B: 'outside' (γ's APL-on response rate over αβ's outside the published ratios; None when none
    is on hand) and 'composition' (at a γ weight of RULE_MAX_W or less, one ρ higher by RULE_GAIN or more on
    RULE_WIRINGS of 8 wirings and against each half, the other not lower by more than RULE_OTHER, and the permuted-label
    control gaining less than RULE_CONTROL times as much); 'run_b' if either holds."""
    raise NotImplementedError


def outcome_a(rule: dict) -> dict:
    """The outcome letter as far as part A decides it."""
    raise NotImplementedError


def part_a(out: Path) -> dict:
    """Part A: A1 to A4, the rule, the outcome, comparison_table.md, realism_table.md, provenance_A.json."""
    raise NotImplementedError


def type_threshold_factors(index: pd.Index, ratios: Mapping[str, float]) -> pd.Series:
    """Part B: per-KC threshold factors (one per type) aiming the types' APL-on response rates at `ratios` to αβ's."""
    raise NotImplementedError


def part_b(out: Path) -> None:
    """Part B (only on the rule): held readings at measured, equal and half-way ratios, with controls."""
    raise NotImplementedError


def part_c(out: Path) -> None:
    """Part C (only on 'better' or 'an equal' in part B): the family fitted again with the types' ratios held."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C'))
    args = ap.parse_args()
    raise NotImplementedError(args.part)


if __name__ == '__main__':
    main()
