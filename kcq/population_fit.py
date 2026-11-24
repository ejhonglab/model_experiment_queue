"""The models as populations at every recorded stage: the three views, their errors in one-fly units and one
population error per stage (E134; the main goal adopted on 2026-11-24, CLAUDE.md "The main goal for new experiments").

A population is the table of trial-mean responses of all units of one stage (glomeruli, bouton regions, claw regions,
KCs) to the odors of one panel: a units x odors DataFrame, odors by their abbreviations, one per fly or per model wiring,
stage, panel and condition (the definition atlas's kcq/stage_readings.py and viki's data_prep/stage_population_data.py
use). Which units respond comes from the stage's response criterion, which the lock gives (the flies: a trial mean above
a threshold in baseline-noise units; the model at the KC somata: two or more spikes). `Population` holds the values and
that boolean table; `as_population` makes one from a bare DataFrame and a criterion, and every function here also
accepts a `Population` made by the caller with its own response table.

The three views, the same at every stage, each made of one or two parts:
- "strength by odor" (`strength`): each odor's response rate (the fraction of units that respond to it) as a multiple
  of the mean over the panel's odors (part `rate`), and each odor's mean response over all units likewise (part `mean`);
- "sharing between odors" (`sharing`): for each pair of odors, the cosine distance between their responses across the
  units, once on who responds (the boolean table, part `who`) and once on how strongly the responders respond (the
  values with the non-responders set to 0, part `how_strongly`);
- "breadth among units" (`breadth`): the distribution of the number of odors a unit responds to, as the fraction of
  units in each bin (part `bins`; bins given by their lower edges, the APL16 panel's 0, 1, 2-3, 4-7, 8 or more by
  default, other panels' bins from the lock).
Where the units are named glomeruli (ORNs, PN dendrites) a fourth view is read, "the matrix itself" (`matrix`): each
glomerulus's response to each odor (part `values`, glomeruli and odors matched by name).

One-fly units. A part's error is the root mean squared difference between the model's part and the flies' mean part,
over the part's entries (the odors, the pairs, the bins), divided by the root mean squared distance of one fly's part
from the mean part of the other flies (over flies and entries together). 1 says that the model is as far from the flies
as one fly is from the others. A view's error is the root mean square of its parts' errors, and a stage's population
error the root mean square of its views' errors: one number per stage, panel and condition, lower is better. (The
planner's probe of 2026-11-24 combined the two sharing parts so: its view errors 1.29 / 1.13 / 1.15 / 0.28 give a
population error of 1.01 only as three views with sharing the root mean square of its two parts.) The population error
is given over the three views (`error`) and, where the matrix view exists, over all four (`error_with_matrix`); which
one enters the megamat stage error is the lock's to say. Every part's and view's error is always returned beside it.

The floor is what a resampled set of the same flies gives against the measured set: the flies are drawn with
replacement, the resampled set's mean views are scored against the measured set's, in the measured set's one-fly
units, and the floor is the mean over draws (about 0.25 for the 11 control recordings). "One recording against the
others" scores each fly against the mean of the other flies, in the units of all flies. The model as a population of
flies: each wiring's error, their mean, and the error of the wirings' mean views (the error common to the wirings);
one wiring against the other wirings taken as the flies is the measure gate's check (1.2 or less).

The steps between stages: the change of each part from the upstream to the downstream stage (entry by entry, on the
odors both stages have), with an interval from resampling the flies of both stages; and the no-step yardstick, the
error of taking the upstream flies' mean views as the prediction for the downstream flies.

The energy distance between the sets of the units' tuning vectors (each unit's responses to the panel's odors) is kept
as a guard only: the model's mean distance to the flies as a multiple of the mean distance between two flies, with an
optional rank-for-rank map of each population's pooled values onto one reference distribution.

General functions here would belong in hong2p or pyolfsysm by the queue's rules; they stay in kcq until E134's measure
gate reads "usable".
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


VIEWS: Tuple[str, ...] = ('strength', 'sharing', 'breadth')
NAMED_VIEW = 'matrix'
PARTS: Dict[str, Tuple[str, ...]] = {
    'strength': ('rate', 'mean'),
    'sharing': ('who', 'how_strongly'),
    'breadth': ('bins',),
    'matrix': ('values',),
}
# the APL16 panel's bins of the number of odors a KC responds to (0, 1, 2-3, 4-7, 8 or more of 16), as the realism band
# on odors per KC (kcq.apl_effect.OPK_BINS); other panels' bins come from the lock
APL16_BREADTH_EDGES: Tuple[int, ...] = (0, 1, 2, 4, 8)
N_BOOT = 1000
N_SPLITS = 200
ENERGY_MAX_UNITS = 2000
SEED = 0

PopulationLike = Union['Population', pd.DataFrame]


@dataclass
class Population:
    """One stage's population of one fly or one model wiring, on one panel and condition.

    Attributes:
        values: units x odors trial-mean responses (odors by their abbreviations; for named stages, units are
            glomerulus names).
        responds: units x odors boolean table of which units respond (same index and columns as `values`).
        name: label of the fly or wiring.
        named_units: True where the units are named glomeruli (ORNs, PN dendrites), so the matrix view is read.
    """
    values: pd.DataFrame
    responds: pd.DataFrame
    name: str = ''
    named_units: bool = False

    @classmethod
    def from_criterion(cls, values: pd.DataFrame, criterion: float, at_or_above: bool = False, name: str = '',
            named_units: bool = False) -> 'Population':
        """Makes a population whose units respond where their value passes the stage's criterion.

        Args:
            values: units x odors trial-mean responses.
            criterion: the stage's response threshold (the lock's).
            at_or_above: True to count values equal to the criterion as responses (the model's "two or more spikes"),
                False for strictly above (the flies' noise-unit threshold).
            name: label of the fly or wiring.
            named_units: whether the units are named glomeruli.
        """
        raise NotImplementedError


def as_population(pop: PopulationLike, criterion: Optional[float] = None, at_or_above: bool = False,
        name: str = '', named_units: bool = False) -> Population:
    """Returns `pop` unchanged if it is a `Population`, else one made from the DataFrame with `criterion`."""
    raise NotImplementedError


def breadth_counts(responds: pd.DataFrame, edges: Sequence[int] = APL16_BREADTH_EDGES) -> pd.Series:
    """Fraction of units whose number of odors responded to falls in each bin.

    Args:
        responds: units x odors boolean table.
        edges: lower edges of the bins, increasing from 0; the last bin is open above.

    Returns:
        Series indexed by bin label ('0', '1', '2_3', '4_7', '8plus' for the default edges).
    """
    raise NotImplementedError


def population_views(pop: Population, odors: Optional[Sequence[str]] = None,
        breadth_edges: Sequence[int] = APL16_BREADTH_EDGES) -> Dict[str, pd.Series]:
    """The views' parts of one population.

    Args:
        pop: the population.
        odors: the odors to read (default all of `pop`'s, in their order); the strength parts are multiples of the
            mean over these odors.
        breadth_edges: lower edges of the breadth bins.

    Returns:
        Dict keyed by '<view>/<part>' ('strength/rate', 'strength/mean', 'sharing/who', 'sharing/how_strongly',
        'breadth/bins', and 'matrix/values' for named units), each a Series of the part's entries (indexed by odor,
        by 'a|b' odor pair, by bin, or by 'glomerulus|odor'). Undefined entries (an odor no unit responds to) are NaN.
    """
    raise NotImplementedError


@dataclass
class FlyReference:
    """The flies' side of one stage, panel and condition.

    Attributes:
        parts: per '<view>/<part>', a flies x entries DataFrame of each fly's part.
        mean: per part, the flies' mean part (NaN-aware over flies).
        one_fly: per part, the root mean squared distance of one fly's part from the mean of the other flies (over
            flies and entries): the part's one-fly unit.
        flies: the flies' names, in row order.
        odors: the odors read.
        breadth_edges: the breadth bins' lower edges.
    """
    parts: Dict[str, pd.DataFrame]
    mean: Dict[str, pd.Series]
    one_fly: Dict[str, float]
    flies: List[str]
    odors: List[str]
    breadth_edges: Tuple[int, ...] = APL16_BREADTH_EDGES


def fly_reference(flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        breadth_edges: Sequence[int] = APL16_BREADTH_EDGES) -> FlyReference:
    """Builds the flies' reference from their populations (at least two flies).

    Args:
        flies: one population per fly (or per wiring, where wirings stand in for flies).
        odors: the odors to read (default those every fly has, in the first fly's order).
        breadth_edges: lower edges of the breadth bins.
    """
    raise NotImplementedError


@dataclass
class Score:
    """Errors in one-fly units of one population (or one set of views) against a `FlyReference`.

    Attributes:
        parts: per '<view>/<part>', the part's error.
        views: per view, the root mean square of its parts' errors.
        error: the population error over the three views (strength, sharing, breadth).
        error_with_matrix: the population error over the four views where the matrix view exists, else NaN.
    """
    parts: Dict[str, float]
    views: Dict[str, float]
    error: float
    error_with_matrix: float = np.nan

    def as_series(self, prefix: str = '') -> pd.Series:
        """One flat Series: 'error', 'error_with_matrix', 'view_<view>' and 'part_<view>_<part>'."""
        raise NotImplementedError


def score_views(views: Mapping[str, pd.Series], ref: FlyReference) -> Score:
    """Scores a set of parts (a model's, a resampled set's mean, the upstream flies' mean) against the flies."""
    raise NotImplementedError


def score(pop: Population, ref: FlyReference) -> Score:
    """Scores one population (one model wiring) against the flies, in their one-fly units."""
    raise NotImplementedError


def mean_views(views: Sequence[Mapping[str, pd.Series]]) -> Dict[str, pd.Series]:
    """Entry-by-entry mean of several sets of parts (NaN-aware)."""
    raise NotImplementedError


def leave_one_out(flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        breadth_edges: Sequence[int] = APL16_BREADTH_EDGES) -> pd.DataFrame:
    """One recording against the others: each fly scored against the mean of the other flies, in all flies' units.

    Returns:
        One row per fly (index the flies' names) with the `Score.as_series` columns.
    """
    raise NotImplementedError


def floor(ref: FlyReference, n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """The floor: resampled sets of the same flies (drawn with replacement) scored against the measured set.

    Returns:
        One row per draw with the `Score.as_series` columns; the floor is their mean (`floor_summary`).
    """
    raise NotImplementedError


def floor_summary(draws: pd.DataFrame) -> pd.Series:
    """Mean of every column of `floor`'s draws."""
    raise NotImplementedError


def wirings_as_flies(wirings: Sequence[Population], ref: FlyReference) -> pd.DataFrame:
    """The model as a population of flies: each wiring's score, their mean, and the score of the wirings' mean views.

    Returns:
        One row per wiring (index the wirings' names), then a row 'mean' (the mean of the wirings' rows) and a row
        'common' (the wirings' mean views scored: the error common to the wirings), with `Score.as_series` columns.
    """
    raise NotImplementedError


def split_halves(n: int, n_splits: int = N_SPLITS, seed: int = SEED) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Random splits of n flies into two halves (sizes floor(n/2) and the rest), as index arrays."""
    raise NotImplementedError


def step_change(upstream: Mapping[str, pd.Series], downstream: Mapping[str, pd.Series]) -> Dict[str, pd.Series]:
    """The change of each part from the upstream to the downstream stage, entry by entry (downstream minus upstream),
    on the entries both have."""
    raise NotImplementedError


def step_interval(upstream: Sequence[Population], downstream: Sequence[Population],
        odors: Optional[Sequence[str]] = None, breadth_edges: Sequence[int] = APL16_BREADTH_EDGES,
        n_boot: int = N_BOOT, ci: float = 0.9, seed: int = SEED) -> Dict[str, pd.DataFrame]:
    """The flies' change of each part between two stages, with an interval from resampling the flies of both stages.

    Returns:
        Per part, a DataFrame indexed by entry with columns 'change' (the measured flies' mean change), 'lo' and 'hi'
        (the `ci` interval of the resampled change).
    """
    raise NotImplementedError


def no_step_score(upstream: Sequence[Population], downstream: FlyReference) -> Score:
    """The no-step yardstick: the upstream flies' mean views as the prediction for the downstream flies (on the odors
    both stages have; parts the upstream stage lacks are left out)."""
    raise NotImplementedError


def rank_map(values: pd.DataFrame, reference: np.ndarray) -> pd.DataFrame:
    """Maps a population's pooled values rank for rank onto a reference distribution (quantile matching)."""
    raise NotImplementedError


def energy_distance(a: np.ndarray, b: np.ndarray, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> float:
    """Energy distance between two sets of units' tuning vectors (rows), sqrt(2 E|X-Y| - E|X-X'| - E|Y-Y'|) with
    Euclidean norms (scipy's convention), on at most `max_units` units drawn from each set."""
    raise NotImplementedError


def energy_guard(pop: Population, flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        reference: Optional[np.ndarray] = None, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> Dict[str, float]:
    """The energy-distance guard: the population's mean energy distance to the flies, the flies' mean distance to each
    other, and their ratio.

    Args:
        pop: the model's (or one fly's) population.
        flies: the flies' populations (a fly identical to `pop`'s name is not compared with itself).
        odors: the odors read (default those all have).
        reference: if given, every population's pooled values are first mapped rank for rank onto it (`rank_map`).

    Returns:
        Dict with 'to_flies', 'fly_to_fly' and 'ratio'.
    """
    raise NotImplementedError
