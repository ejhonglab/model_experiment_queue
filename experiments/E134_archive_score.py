#!/usr/bin/env python3
"""E134 step 2, the archive at the KC somata: every archived variant on the APL16 panel as a population.

The entry is item 2 of experiment_priorities.md ("the models as populations at every recorded stage"), design steps 6
and 7. For every variant of results/apl_effect/by_variant.csv whose APL-on and acute APL-off spike counts are in the
model-output archive (data/model_outputs/<folder>/ or a subfolder; `analyses.apl_effect_metrics.archived_counts`) and
whose odors are the APL16 panel's 16:
- its KC population with APL on is scored against the control flies and with APL off against APL output blocked (TNT),
  with kcq/population_fit.py: the population error, each view and part in one-fly units, and the energy-distance guard;
- the order of the variants is compared between two random halves of the recordings (each genotype's recordings split
  in two), 200 times: the ρ over variants between a measure's values against one half and against the other, for the
  population error, each view, the energy distance and the two ρ's (ρ on odor-odor distances and response rate ρ),
  over all distinct variants and over those inside the realism bands (the measure's gate: 0.8 or more);
- over the variants inside the bands, the population error against the calibrated pattern part (the APL-effect
  statistics, a held-out test under the adopted goal), and which variant and family is best by each.

The flies are the lock's populations (`data_prep.stage_population_data.load_population('kc_soma', 'apl16', 'control')`
and `'tnt'`): trial means in baseline-noise units, a response above the lock's criterion. The model is read as the lock
says (`--reading`): 'raw' is the spike counts with a response at two or more spikes (the planner's probe, no curve and no
noise); the other readings are atlas's kcq/stage_readings.py's for the KC somata, once pushed. Variants whose per-pair
effect duplicates an earlier one (the same fit saved by several experiments) count once in the consistency.

No model is read before the lock is committed: `main` refuses to run until the lock (experiments/E134_lock.json) is
in git.

Writes <shared>/results/E134_population_baseline/:
    archive_scores_<reading>.csv          one row per folder x tag x condition (on / off): the Score columns of
                                          kcq.population_fit, the energy guard, the realism flag, the duplicate flag
    archive_consistency_<reading>.csv     per condition, measure and variant set: the split-half ρ (mean, 5 and 95 %)
    archive_tradeoff_<reading>.csv        per variant inside the bands: population errors, calibrated pattern part,
                                          and the best variant and family by each
    archive_provenance_<reading>.json
and its per-variant views to outputs/E134_population_baseline/archive_views_<reading>.parquet.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_archive_score.py [--reading raw] [--n-splits 200] [--jobs 8] [--seed 0]
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import population_fit as P

EXPERIMENT = 'E134_population_baseline'
CONDITIONS = {'on': ('on', 'control'), 'off': ('acute_off', 'tnt')}   # model condition -> (archive name, fly cohort)
N_SPLITS = 200
SEED = 0
# the planner's probe: a model KC responds at two or more spikes, on raw counts
RAW_CRITERION = 2.0
# a reading maps one variant's spike counts (KCs x odors) to its population (values and response table)
Reading = Callable[[pd.DataFrame], P.Population]


def raw_reading(counts: pd.DataFrame) -> P.Population:
    """The 'raw' reading: spike counts as they are, a response at two or more spikes."""
    raise NotImplementedError


def readings() -> Dict[str, Reading]:
    """The readings `--reading` may name: 'raw', and those of kcq/stage_readings.py for the KC somata when present."""
    raise NotImplementedError


def lock_committed() -> dict:
    """The committed lock (experiments/E134_lock.json, tracked by git with no uncommitted change); raises RuntimeError
    otherwise (no model is read before the lock)."""
    raise NotImplementedError


def fly_populations(cohort: str) -> List[P.Population]:
    """The APL16 KC somata of one cohort ('control' or 'tnt') from the lock's data module, one population per
    recording."""
    raise NotImplementedError


def archived_variants(by_variant: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Every variant of results/apl_effect/by_variant.csv with archived APL-on and acute APL-off spike counts: columns
    exp, folder, tag, archive_dir, realism_flag, realism, and the APL-effect statistics by_variant.csv carries."""
    raise NotImplementedError


def read_counts(archive_dir: Path, tag: str, name: str) -> Optional[pd.DataFrame]:
    """A variant's spike counts (KCs x odors, odors abbreviated, in the APL16 order), or None if its odors are not
    the APL16 panel's."""
    raise NotImplementedError


def variant_views(variants: pd.DataFrame, reading: Reading, odors: Sequence[str], breadth_edges: Sequence[int],
        jobs: int = 1) -> Tuple[pd.DataFrame, Dict[Tuple[str, str, str], Dict[str, pd.Series]]]:
    """Each variant's views per condition, and its per-recording-independent summaries (response rate per odor,
    pair distances) for the two ρ's.

    Returns:
        (the variants that could be read, with a `duplicate` column; {(folder, tag, condition): parts}).
    """
    raise NotImplementedError


def score_variants(views: Dict[Tuple[str, str, str], Dict[str, pd.Series]],
        refs: Dict[str, P.FlyReference]) -> pd.DataFrame:
    """Every variant's Score against the flies of its condition (`refs`: condition -> FlyReference)."""
    raise NotImplementedError


def energy_scores(variants: pd.DataFrame, reading: Reading, flies: Dict[str, List[P.Population]],
        jobs: int = 1, seed: int = SEED) -> pd.DataFrame:
    """Each variant's energy distance to each recording of its condition (one column per recording), so that any
    subset of the recordings gives the guard without recomputing."""
    raise NotImplementedError


def consistency(views: Dict[Tuple[str, str, str], Dict[str, pd.Series]], flies: Dict[str, List[P.Population]],
        variants: pd.DataFrame, energy: Optional[pd.DataFrame] = None, n_splits: int = N_SPLITS,
        seed: int = SEED, units_from: str = 'half') -> pd.DataFrame:
    """The split-half consistency of every measure over the variants.

    Args:
        views: per (folder, tag, condition) the variant's parts.
        flies: condition -> the recordings.
        variants: the variants, with `duplicate` and `realism` columns.
        energy: `energy_scores`' table, for the energy distance's consistency.
        n_splits: random splits of each condition's recordings in two.
        seed: of the splits.
        units_from: 'half' to score against a half in that half's own one-fly units (the half read as the flies),
            'all' to keep all recordings' one-fly units.

    Returns:
        One row per condition, measure and variant set ('all', 'inside_bands'): split_half_rho (mean over splits), lo,
        hi (5 and 95 %), n_variants, n_splits.
    """
    raise NotImplementedError


def tradeoff(scores: pd.DataFrame, variants: pd.DataFrame, calibrated: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Over the variants inside the bands: the population errors with APL on and off beside the calibrated pattern
    part, with the best variant and the best family (experiment folder) by each."""
    raise NotImplementedError


def calibrated_pattern_parts(variants: pd.DataFrame, jobs: int = 1) -> pd.DataFrame:
    """Each variant's pattern part under the calibrated readout (as analyses/apl_effect_metrics.py --readout cal
    reads it: K calibrated on the APL-on counts to the control flies' responders-only participation ratio)."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--reading', default='raw')
    ap.add_argument('--n-splits', type=int, default=N_SPLITS)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--units-from', choices=('half', 'all'), default='half')
    ap.add_argument('--no-energy', action='store_true')
    ap.parse_args(argv)
    raise NotImplementedError


if __name__ == '__main__':
    main()
