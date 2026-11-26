"""E142's step 0 (no simulation): the flies' step by odor, every saved run's entries by odor, the count-weighted
yardsticks, what the model's bouton reading samples, the equal-weight reading of the boutons and the classes read again,
the planner's probe numbers reproduced, and step 0's classes. The definitions, the saved runs, the halves, the seeds and
the thresholds are those of experiments/E142_dendrite_bouton_step.py and its committed lock (experiments/E142_lock.json);
every piece calls `E142.require_all` first.

Pieces (`PIECES`, run in this order by `run`):
- flies: flies_steps.csv (each step's change by odor and part on all flies and each half, with its interval),
  flies_steps_summary.csv (the odors outside zero on both halves), no_step_entries.csv (the no-step yardstick's entries
  by odor against all bouton flies and each half, and the same for the step to the claws).
- models: models_entries.csv (every saved run's entries by odor at the dendrites, boutons and claws, each wiring and
  'common', on every fly set), models_concentrated.csv, models_step.csv (the model's step against the flies', r over the
  17 odors).
- yardsticks: counts.csv (each glomerulus's boutons and claws in each calyx, against E138's feeders.csv),
  yardsticks.csv (the count-weighted yardsticks, the no-step yardsticks and the null's percentiles), yardstick_null.csv.
- sampling: bouton_sampling.csv (the model's bouton units per glomerulus beside the counts), bouton_sampling_summary.csv.
- equal: equal_weight_groups.csv (the flies' groups), equal_weight_errors.csv (the flies and every saved run under the
  equal-weight reading), reclasses.json (E136's and E139's classes read again).
- probe: probe.csv (the planner's numbers against the reproduced ones).
- classes: classes_step0.json, tables.md.
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E142_dendrite_bouton_step as E142

PIECES: Tuple[str, ...] = ('flies', 'models', 'yardsticks', 'sampling', 'equal', 'probe', 'classes')


def weighted_strength(pop, weights: Mapping[str, float], odors: Sequence[str]) -> Dict[str, pd.Series]:
    """The strength parts of a named population with each unit (glomerulus) weighted by `weights` (0 where absent):
    equal to hal's `population_views` of the population with each row repeated its weight's number of times."""
    raise NotImplementedError


def flies_steps(lock: Mapping) -> Dict[str, pd.DataFrame]:
    """The flies' step by odor (STEPS) on every fly set, with intervals; the odors outside zero on both halves; the
    no-step yardstick's entries by odor."""
    raise NotImplementedError


def models_entries(lock: Mapping) -> Dict[str, pd.DataFrame]:
    """Every saved run's entries by odor at ENTRY_STAGES on every fly set, the concentrated class and the model's step
    against the flies'."""
    raise NotImplementedError


def yardsticks(lock: Mapping) -> Dict[str, pd.DataFrame]:
    """The counts, the count-weighted yardsticks with their null, and the carried glomeruli beside."""
    raise NotImplementedError


def bouton_sampling(lock: Mapping) -> Dict[str, pd.DataFrame]:
    """What the model's bouton reading samples, per wiring."""
    raise NotImplementedError


def equal_weight(lock: Mapping) -> dict:
    """The equal-weight reading of the flies and of every saved run, and the classes read again."""
    raise NotImplementedError


def probe(lock: Mapping) -> pd.DataFrame:
    """The planner's probe numbers (E142.PROBE) against step 0's reproductions."""
    raise NotImplementedError


def classes(lock: Mapping) -> dict:
    """Step 0's classes ('counts_make_step', 'readings') and tables.md."""
    raise NotImplementedError


def run(pieces: Optional[Sequence[str]] = None) -> dict:
    """Step 0's pieces in order (all of PIECES by default)."""
    raise NotImplementedError
