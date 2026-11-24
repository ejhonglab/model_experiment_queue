"""The models as populations at every recorded stage: the three views, their errors in one-fly units and one
population error per stage (E134; the main goal adopted on 2026-11-24, CLAUDE.md "The main goal for new experiments").

A population is the table of trial-mean responses of all units of one stage (glomeruli, bouton regions, claw regions,
KCs) to the odors of one panel: a units x odors DataFrame, odors by their abbreviations, one per fly or per model wiring,
stage, panel and condition (the definition atlas's kcq/stage_readings.py and viki's data_prep/stage_population_data.py
use). Which units respond comes from the stage's response criterion, which the lock gives (the flies: a trial mean above
a threshold in baseline-noise units; the model at the KC somata: two or more spikes). `Population` holds one fly's or
one wiring's values and that response table; `as_population` makes one from a bare DataFrame and a criterion, and
`fly_populations` makes one per fly from the data module's recorded population (whose `Population` holds all flies of
a stage; the names differ by what they hold, and every function here takes this module's).

The three views, the same at every stage, each made of one or two parts:
- "strength by odor" (`strength`): each odor's response rate (the fraction of units that respond to it) as a multiple
  of the mean over the panel's odors (part `rate`), and each odor's mean response over all units likewise (part `mean`;
  NaN where the panel's mean is not positive);
- "sharing between odors" (`sharing`): for each pair of odors, the cosine distance between their responses across the
  units, once on who responds (the response table, part `who`) and once on how strongly the responders respond (the
  values with the non-responders set to 0, part `how_strongly`; the same as the cosine over the units that respond to
  either odor). A pair with an odor no unit responds to is NaN;
- "breadth among units" (`breadth`): the distribution of the number of odors a unit responds to, as the fraction of
  units in each bin (part `bins`). Bins are given by their lower edges (the APL16 panel's 0, 1, 2-3, 4-7, 8 or more by
  default) or as closed (lo, hi) pairs (the data module's `breadth_bins`, scaled to the panel's size by the lock).
Where the units are named glomeruli (ORNs, PN dendrites) a fourth view is read, "the matrix itself" (`matrix`): each
glomerulus's response to each odor (part `values`, entries 'glomerulus|odor', matched by name). A population with
binary responses only (`graded=False`) has no `strength/mean`, `sharing/how_strongly` or `matrix/values`.

One-fly units. A part's error is the root mean squared difference between the model's part and the flies' mean part,
over the part's entries (the odors, the pairs, the bins), divided by the root mean squared distance of one fly's part
from the mean part of the other flies (over flies and entries together). 1 says that the model is as far from the flies
as one fly is from the others. A view's error is the root mean square of its parts' errors, and a stage's population
error the root mean square of its views' errors: one number per stage, panel and condition, lower is better. (The
planner's probe of 2026-11-24 combined the two sharing parts so: its view errors 1.29 / 1.13 / 1.15 / 0.28 give a
population error of 1.01 only as three views with sharing the root mean square of its two parts.) The population error
is given over the three views (`error`) and, where the matrix view exists, over all four (`error_with_matrix`); which
one enters the megamat stage error is the lock's to say. Every part's and view's error is always returned beside it,
and `fly_reference(parts=...)` reads a subset of the parts (the probe's strength was the rate part only). Entries
undefined in the model (NaN) are left out of its part's root mean square and counted (`Score.missing`); a part with no
defined entry is NaN, and so are its view and the population error.

The floor is what a resampled set of the same flies gives against the measured set: the flies are drawn with
replacement, the resampled set's mean parts are scored against the measured set's, in the measured set's one-fly
units, and the floor is the mean over draws (about 0.25 for the 11 control recordings). "One recording against the
others" scores each fly against the mean of the other flies, in the units of all flies; the root mean square of its
part errors over the flies is 1 by construction. The model as a population of flies: each wiring's error, their mean,
and the error of the wirings' mean parts (the error common to the wirings); one wiring against the other wirings taken
as the flies (`leave_one_out` on the wirings) is the measure gate's check (1.2 or less). `subset_reference` reads the
flies of a half (split-half consistency) in that half's own units or in all flies'.

The steps between stages: the change of each part from the upstream to the downstream stage (entry by entry, on the
odors both stages have), with an interval from resampling the flies of both stages (the description); the step's
resampled distance, by which a step is classed (the planner's ruling of 2026-11-24 11:34): per part, the RMS distance
of a model's per-entry change from the flies' change, inside if at or below the 95th percentile of the same distance
for resampled flies (`step_resample`, `step_distance`; how many wirings must be inside is the runner's); and the
no-step yardstick, the error of taking the upstream flies' mean parts as the prediction for the downstream flies.

The energy distance between the sets of the units' tuning vectors (each unit's responses to the panel's odors) is kept
as a guard only: the model's mean distance to the flies as a multiple of the mean distance between two flies, with an
optional rank-for-rank map of each population's pooled values onto one reference distribution.

General functions here would belong in hong2p or pyolfsysm by the queue's rules; they stay in kcq until E134's measure
gate reads "usable".
"""
from __future__ import annotations

import dataclasses
import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist


VIEWS: Tuple[str, ...] = ('strength', 'sharing', 'breadth')
NAMED_VIEW = 'matrix'
PARTS: Dict[str, Tuple[str, ...]] = {
    'strength': ('rate', 'mean'),
    'sharing': ('who', 'how_strongly'),
    'breadth': ('bins',),
    'matrix': ('values',),
}
ALL_PARTS: Tuple[str, ...] = tuple(f'{v}/{p}' for v in (*VIEWS, NAMED_VIEW) for p in PARTS[v])
GRADED_PARTS = ('strength/mean', 'sharing/how_strongly', 'matrix/values')
# the APL16 panel's bins of the number of odors a KC responds to (0, 1, 2-3, 4-7, 8 or more of 16), as the realism band
# on odors per KC (kcq.apl_effect.OPK_BINS); other panels' bins come from the lock
APL16_BREADTH_EDGES: Tuple[int, ...] = (0, 1, 2, 4, 8)
N_BOOT = 1000
N_SPLITS = 200
ENERGY_MAX_UNITS = 2000
SEED = 0
PAIR_SEP = '|'

BinSpec = Sequence[Union[int, Tuple[int, int]]]


@dataclass
class Population:
    """One stage's population of one fly or one model wiring, on one panel and condition.

    Attributes:
        values: units x odors trial-mean responses (odors by their abbreviations; for named stages, units are
            glomerulus names). NaN where a unit was not given an odor.
        responds: units x odors table of which units respond (1.0 / 0.0, NaN where `values` is NaN; booleans are
            converted).
        name: label of the fly or wiring.
        named_units: True where the units are named glomeruli (ORNs, PN dendrites), so the matrix view is read.
        graded: False where only binary responses are on hand (the parts on how strongly are not read).
    """
    values: pd.DataFrame
    responds: pd.DataFrame
    name: str = ''
    named_units: bool = False
    graded: bool = True

    def __post_init__(self):
        v = self.values.astype(float)
        r = self.responds.reindex(index=v.index, columns=v.columns).astype(float)
        self.values = v
        self.responds = r.where(v.notna())

    @classmethod
    def from_criterion(cls, values: pd.DataFrame, criterion: float, at_or_above: bool = False, name: str = '',
            named_units: bool = False, graded: bool = True) -> 'Population':
        """Makes a population whose units respond where their value passes the stage's criterion.

        Args:
            values: units x odors trial-mean responses.
            criterion: the stage's response threshold (the lock's).
            at_or_above: True to count values equal to the criterion as responses (the model's "two or more spikes"),
                False for strictly above (the flies' noise-unit threshold).
            name: label of the fly or wiring.
            named_units: whether the units are named glomeruli.
            graded: whether the values are graded responses (False for 0 / 1 tables).
        """
        v = values.astype(float)
        r = (v >= criterion) if at_or_above else (v > criterion)
        return cls(v, r.astype(float).where(v.notna()), name=name, named_units=named_units, graded=graded)

    @property
    def odors(self) -> List[str]:
        """The odors the population was given (columns not all NaN)."""
        return [o for o in self.values.columns if self.values[o].notna().any()]


PopulationLike = Union[Population, pd.DataFrame]


def as_population(pop: PopulationLike, criterion: Optional[float] = None, at_or_above: bool = False,
        name: str = '', named_units: bool = False) -> Population:
    """Returns `pop` unchanged if it is a `Population`, else one made from the DataFrame with `criterion`."""
    if isinstance(pop, Population):
        return pop
    if criterion is None:
        raise ValueError('a bare DataFrame needs the response criterion')
    return Population.from_criterion(pop, criterion, at_or_above=at_or_above, name=name, named_units=named_units)


def fly_populations(recorded) -> List[Population]:
    """One `Population` per fly of a recorded population of data_prep/stage_population_data.py (its `tables`, `flies`,
    `criterion`, `named_units` and `graded`), a response being a trial mean above the criterion."""
    out = []
    for fly in recorded.flies:
        t = recorded.tables[fly]
        try:
            r = recorded.responses(fly)
        except NotImplementedError:
            r = (t.astype(float) > recorded.criterion).astype(float).where(t.notna())
        out.append(Population(t, r, name=str(fly), named_units=bool(recorded.named_units),
            graded=bool(getattr(recorded, 'graded', True))))
    return out


def breadth_bins(spec: BinSpec = APL16_BREADTH_EDGES) -> List[Tuple[int, float]]:
    """Closed (lo, hi) bins from lower edges (the last open above: hi infinite) or from (lo, hi) pairs."""
    spec = list(spec)
    if not spec:
        raise ValueError('no breadth bins')
    if all(isinstance(s, (tuple, list)) for s in spec):
        bins = [(int(lo), float(hi)) for lo, hi in spec]
    else:
        e = [int(x) for x in spec]
        if e[0] != 0 or any(b <= a for a, b in zip(e, e[1:])):
            raise ValueError(f'breadth edges must rise from 0: {e}')
        bins = [(lo, float(nxt - 1)) for lo, nxt in zip(e, e[1:])] + [(e[-1], math.inf)]
    return bins


def breadth_labels(spec: BinSpec = APL16_BREADTH_EDGES) -> List[str]:
    """Labels of the bins: '0', '1', '2_3', '4_7', '8plus' for the default edges; '8_16' for a closed pair."""
    out = []
    for lo, hi in breadth_bins(spec):
        out.append(f'{lo}plus' if math.isinf(hi) else (str(lo) if hi == lo else f'{lo}_{int(hi)}'))
    return out


def breadth_counts(responds: pd.DataFrame, edges: BinSpec = APL16_BREADTH_EDGES) -> pd.Series:
    """Fraction of units whose number of odors responded to falls in each bin.

    Args:
        responds: units x odors table of responses (1 / 0 or booleans; NaN not counted, so a unit missing an odor
            is binned by the odors it was given).
        edges: lower edges of the bins (the last open above), or closed (lo, hi) pairs.

    Returns:
        Series indexed by bin label ('0', '1', '2_3', '4_7', '8plus' for the default edges).
    """
    n = np.nansum(responds.to_numpy(dtype=float), axis=1)
    bins, labels = breadth_bins(edges), breadth_labels(edges)
    if len(n) == 0:
        return pd.Series(np.nan, index=labels)
    return pd.Series([float(np.mean((n >= lo) & (n <= hi))) for lo, hi in bins], index=labels)


def odor_pairs(odors: Sequence[str]) -> List[str]:
    """'a|b' labels of the pairs of `odors`, a before b in the given order."""
    return [f'{a}{PAIR_SEP}{b}' for a, b in itertools.combinations(odors, 2)]


def cosine_pairs(X: np.ndarray) -> np.ndarray:
    """Upper-triangle cosine distances between the columns of X (units x odors), in `odor_pairs`' order (NaN counted
    as 0; a column that is all NaN or all 0 gives NaN)."""
    absent = np.all(np.isnan(X), axis=0)
    X0 = np.nan_to_num(X, nan=0.0)
    norm = np.sqrt((X0 ** 2).sum(axis=0))
    bad = absent | (norm == 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        C = (X0.T @ X0) / np.outer(norm, norm)
    D = 1.0 - C
    D[bad, :] = np.nan
    D[:, bad] = np.nan
    iu = np.triu_indices(X.shape[1], k=1)
    return D[iu]


def _relative(x: pd.Series) -> pd.Series:
    m = float(np.nanmean(x.to_numpy(dtype=float))) if x.notna().any() else np.nan
    if not np.isfinite(m) or m <= 0:
        return pd.Series(np.nan, index=x.index)
    return x / m


def population_views(pop: Population, odors: Optional[Sequence[str]] = None,
        breadth_edges: BinSpec = APL16_BREADTH_EDGES, parts: Optional[Sequence[str]] = None) -> Dict[str, pd.Series]:
    """The views' parts of one population.

    Args:
        pop: the population.
        odors: the odors to read (default all of `pop`'s columns, in their order); the strength parts are multiples
            of the mean over these odors and breadth counts responses to these odors only.
        breadth_edges: the breadth bins (lower edges or closed pairs).
        parts: the parts to compute (default all that the population allows: `ALL_PARTS` less the graded parts when
            `pop.graded` is False and the matrix when its units are not named).

    Returns:
        Dict keyed by '<view>/<part>' ('strength/rate', 'strength/mean', 'sharing/who', 'sharing/how_strongly',
        'breadth/bins', and 'matrix/values' for named units), each a Series of the part's entries (indexed by odor,
        by 'a|b' odor pair, by bin, or by 'glomerulus|odor'). Undefined entries are NaN; an odor the population was
        not given makes its entries NaN and its breadth part NaN.
    """
    odors = list(pop.values.columns) if odors is None else list(odors)
    allowed = [p for p in ALL_PARTS if (pop.graded or p not in GRADED_PARTS) and (pop.named_units or not
        p.startswith(NAMED_VIEW + '/'))]
    want = allowed if parts is None else [p for p in parts if p in allowed]
    V = pop.values.reindex(columns=odors)
    R = pop.responds.reindex(columns=odors)
    out: Dict[str, pd.Series] = {}
    if 'strength/rate' in want:
        out['strength/rate'] = _relative(R.mean(axis=0))
    if 'strength/mean' in want:
        out['strength/mean'] = _relative(V.mean(axis=0))
    pairs = odor_pairs(odors)
    if 'sharing/who' in want:
        out['sharing/who'] = pd.Series(cosine_pairs(R.to_numpy(dtype=float)), index=pairs)
    if 'sharing/how_strongly' in want:
        W = V.where(R > 0.5, 0.0).where(V.notna())
        out['sharing/how_strongly'] = pd.Series(cosine_pairs(W.to_numpy(dtype=float)), index=pairs)
    if 'breadth/bins' in want:
        if V.isna().all(axis=0).any():
            out['breadth/bins'] = pd.Series(np.nan, index=breadth_labels(breadth_edges))
        else:
            out['breadth/bins'] = breadth_counts(R, breadth_edges)
    if 'matrix/values' in want:
        out['matrix/values'] = pd.Series(V.to_numpy(dtype=float).ravel(),
            index=[f'{g}{PAIR_SEP}{o}' for g in V.index for o in V.columns])
    return out


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
        breadth_edges: the breadth bins.
    """
    parts: Dict[str, pd.DataFrame]
    mean: Dict[str, pd.Series]
    one_fly: Dict[str, float]
    flies: List[str]
    odors: List[str]
    breadth_edges: Tuple = APL16_BREADTH_EDGES


def _loo_diffs(A: np.ndarray) -> np.ndarray:
    """Each row minus the NaN-aware mean of the other rows (NaN where undefined)."""
    finite = np.isfinite(A)
    S = np.where(finite, A, 0.0).sum(axis=0)
    C = finite.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        others = (S[None, :] - np.where(finite, A, 0.0)) / (C[None, :] - finite)
    return np.where(finite & ((C[None, :] - finite) > 0), A - others, np.nan)


def _one_fly_unit(A: np.ndarray) -> float:
    d = _loo_diffs(A)
    return float(np.sqrt(np.nanmean(d ** 2))) if np.isfinite(d).any() else np.nan


def _nanmean_rows(A: np.ndarray) -> np.ndarray:
    with np.errstate(invalid='ignore', divide='ignore'):
        finite = np.isfinite(A)
        return np.where(finite.any(axis=-2), np.where(finite, A, 0.0).sum(axis=-2) / finite.sum(axis=-2), np.nan)


def _reference_from_parts(parts: Dict[str, pd.DataFrame], flies: List[str], odors: List[str], edges,
        one_fly: Optional[Dict[str, float]] = None) -> FlyReference:
    mean = {p: pd.Series(_nanmean_rows(df.to_numpy(dtype=float)), index=df.columns) for p, df in parts.items()}
    if one_fly is None:
        one_fly = {p: _one_fly_unit(df.to_numpy(dtype=float)) for p, df in parts.items()}
    return FlyReference(parts=parts, mean=mean, one_fly=dict(one_fly), flies=list(flies), odors=list(odors),
        breadth_edges=tuple(edges))


def common_odors(pops: Sequence[Population]) -> List[str]:
    """The odors every population was given, in the first one's column order."""
    have = [set(p.odors) for p in pops]
    return [o for o in pops[0].values.columns if all(o in h for h in have)]


def fly_reference(flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        breadth_edges: BinSpec = APL16_BREADTH_EDGES, parts: Optional[Sequence[str]] = None) -> FlyReference:
    """Builds the flies' reference from their populations (at least two flies).

    Args:
        flies: one population per fly (or per wiring, where wirings stand in for flies).
        odors: the odors to read (default those every fly has, in the first fly's order).
        breadth_edges: the breadth bins (lower edges or closed pairs).
        parts: the parts to read (default every part all flies allow).
    """
    if len(flies) < 2:
        raise ValueError('a reference needs at least two flies')
    odors = common_odors(flies) if odors is None else list(odors)
    views = [population_views(f, odors, breadth_edges, parts) for f in flies]
    keep = [p for p in ALL_PARTS if all(p in v for v in views)]
    names = [f.name or str(i) for i, f in enumerate(flies)]
    table = {}
    for p in keep:
        cols = list(dict.fromkeys(itertools.chain.from_iterable(v[p].index for v in views)))
        table[p] = pd.DataFrame([v[p].reindex(cols).to_numpy(dtype=float) for v in views], index=names, columns=cols)
    return _reference_from_parts(table, names, odors, breadth_edges)


def subset_reference(ref: FlyReference, flies: Sequence[Union[int, str]], keep_units: bool = False) -> FlyReference:
    """The reference of a subset of its flies (by position or name), e.g. one half of the recordings.

    Args:
        ref: the full reference.
        flies: the flies kept.
        keep_units: True to keep `ref`'s one-fly units; False (default) to compute the subset's own (the subset read
            as the flies).
    """
    idx = [ref.flies.index(f) if isinstance(f, str) else int(f) for f in flies]
    parts = {p: df.iloc[idx] for p, df in ref.parts.items()}
    return _reference_from_parts(parts, [ref.flies[i] for i in idx], ref.odors, ref.breadth_edges,
        ref.one_fly if keep_units else None)


@dataclass
class Score:
    """Errors in one-fly units of one population (or one set of parts) against a `FlyReference`.

    Attributes:
        parts: per '<view>/<part>', the part's error.
        views: per view, the root mean square of its parts' errors.
        error: the population error over the three views (strength, sharing, breadth) the reference has.
        error_with_matrix: the population error over the four views where the matrix view exists, else NaN.
        missing: per part, the entries the flies define and the scored side does not.
    """
    parts: Dict[str, float]
    views: Dict[str, float]
    error: float
    error_with_matrix: float = np.nan
    missing: Dict[str, int] = field(default_factory=dict)

    def as_series(self, prefix: str = '') -> pd.Series:
        """One flat Series: 'error', 'error_with_matrix', 'view_<view>', 'part_<view>_<part>' and 'n_missing'."""
        d = {'error': self.error, 'error_with_matrix': self.error_with_matrix}
        d.update({f'view_{v}': e for v, e in self.views.items()})
        d.update({'part_' + p.replace('/', '_'): e for p, e in self.parts.items()})
        d['n_missing'] = float(sum(self.missing.values()))
        return pd.Series(d).add_prefix(prefix)


def _rms(xs: Sequence[float]) -> float:
    xs = list(xs)
    return float(np.sqrt(np.mean(np.square(xs)))) if xs else np.nan


def _combine(part_err: Dict[str, float]) -> Tuple[Dict[str, float], float, float]:
    views = {}
    for v in (*VIEWS, NAMED_VIEW):
        ps = [part_err[f'{v}/{p}'] for p in PARTS[v] if f'{v}/{p}' in part_err]
        if ps:
            views[v] = _rms(ps)
    three = _rms([views[v] for v in VIEWS if v in views])
    four = _rms([views[v] for v in (*VIEWS, NAMED_VIEW) if v in views]) if NAMED_VIEW in views else np.nan
    return views, three, four


def stack_views(views: Sequence[Mapping[str, pd.Series]], ref: FlyReference) -> Dict[str, np.ndarray]:
    """Per part of `ref`, an array (one row per set of parts in `views`, columns the reference's entries)."""
    out = {}
    for p, m in ref.mean.items():
        out[p] = np.array([np.asarray(v[p].reindex(m.index), dtype=float) if p in v else np.full(len(m), np.nan)
            for v in views]).reshape(len(views), len(m))
    return out


def score_arrays(arrays: Mapping[str, np.ndarray], ref: FlyReference, index: Optional[Sequence] = None) -> pd.DataFrame:
    """Scores many sets of parts at once (`stack_views`' arrays) against the flies; one row per set, with the
    `Score.as_series` columns."""
    return score_rows(arrays, {p: ref.mean[p].to_numpy(dtype=float) for p in arrays},
        {p: ref.one_fly[p] for p in arrays}, index)


def score_rows(arrays: Mapping[str, np.ndarray], means: Mapping[str, np.ndarray],
        units: Mapping[str, Union[float, np.ndarray]], index: Optional[Sequence] = None) -> pd.DataFrame:
    """`score_arrays` with a reference of its own for each row (e.g. one set of parts scored against each half of
    many splits of the recordings at once).

    Args:
        arrays: per part, sets x entries (`stack_views`' arrays, columns the reference's entries).
        means: per part, the flies' mean part: entries (one reference for every row) or sets x entries.
        units: per part, the one-fly unit: a number or one per set.
        index: the rows' labels.

    Returns:
        One row per set, with the `Score.as_series` columns.
    """
    n = len(next(iter(arrays.values())))
    cols: Dict[str, np.ndarray] = {}
    n_missing = np.zeros(n)
    for p, A in arrays.items():
        m = np.broadcast_to(np.asarray(means[p], dtype=float), A.shape)
        d = A - m
        defined = np.isfinite(m)
        n_missing += (defined & ~np.isfinite(A)).sum(axis=1)
        with np.errstate(invalid='ignore'):
            ok = np.isfinite(d)
            msd = np.where(ok.any(axis=1), np.where(ok, d ** 2, 0.0).sum(axis=1) / np.maximum(ok.sum(axis=1), 1), np.nan)
            cols['part_' + p.replace('/', '_')] = np.sqrt(msd) / np.asarray(units[p], dtype=float)
    out = pd.DataFrame(cols, index=index)
    for v in (*VIEWS, NAMED_VIEW):
        ps = ['part_' + f'{v}_{p}' for p in PARTS[v] if f'{v}/{p}' in arrays]
        if ps:
            out[f'view_{v}'] = np.sqrt((out[ps] ** 2).mean(axis=1, skipna=False))
    three = [f'view_{v}' for v in VIEWS if f'view_{v}' in out]
    out['error'] = np.sqrt((out[three] ** 2).mean(axis=1, skipna=False)) if three else np.nan
    four = three + ([f'view_{NAMED_VIEW}'] if f'view_{NAMED_VIEW}' in out else [])
    out['error_with_matrix'] = np.sqrt((out[four] ** 2).mean(axis=1, skipna=False)) if len(four) > len(three) else np.nan
    out['n_missing'] = n_missing
    first = ['error', 'error_with_matrix'] + [c for c in out if c.startswith('view_')]
    return out[first + [c for c in out if c.startswith('part_')] + ['n_missing']]


def _score_from_row(row: pd.Series, missing: Dict[str, int]) -> Score:
    parts = {p: float(row['part_' + p.replace('/', '_')]) for p in ALL_PARTS if 'part_' + p.replace('/', '_') in row}
    views = {v: float(row[f'view_{v}']) for v in (*VIEWS, NAMED_VIEW) if f'view_{v}' in row}
    return Score(parts=parts, views=views, error=float(row['error']), error_with_matrix=float(row['error_with_matrix']),
        missing=missing)


def score_views(views: Mapping[str, pd.Series], ref: FlyReference) -> Score:
    """Scores a set of parts (a model's, a resampled set's mean, the upstream flies' mean) against the flies."""
    arrays = stack_views([views], ref)
    row = score_arrays(arrays, ref).iloc[0]
    missing = {p: int((np.isfinite(ref.mean[p].to_numpy()) & ~np.isfinite(A[0])).sum()) for p, A in arrays.items()}
    return _score_from_row(row, missing)


def score(pop: Population, ref: FlyReference) -> Score:
    """Scores one population (one model wiring) against the flies, in their one-fly units."""
    return score_views(population_views(pop, ref.odors, ref.breadth_edges, list(ref.parts)), ref)


def score_many(pops: Sequence[Population], ref: FlyReference) -> pd.DataFrame:
    """`score` of every population, one row each (index the populations' names)."""
    views = [population_views(p, ref.odors, ref.breadth_edges, list(ref.parts)) for p in pops]
    return score_arrays(stack_views(views, ref), ref, index=[p.name or str(i) for i, p in enumerate(pops)])


def mean_views(views: Sequence[Mapping[str, pd.Series]]) -> Dict[str, pd.Series]:
    """Entry-by-entry mean of several sets of parts (NaN-aware; parts all sets have)."""
    keep = [p for p in ALL_PARTS if all(p in v for v in views)]
    out = {}
    for p in keep:
        idx = list(dict.fromkeys(itertools.chain.from_iterable(v[p].index for v in views)))
        A = np.array([v[p].reindex(idx).to_numpy(dtype=float) for v in views])
        out[p] = pd.Series(_nanmean_rows(A), index=idx)
    return out


def leave_one_out(flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        breadth_edges: BinSpec = APL16_BREADTH_EDGES, parts: Optional[Sequence[str]] = None,
        ref: Optional[FlyReference] = None) -> pd.DataFrame:
    """One recording against the others: each fly scored against the mean of the other flies, in all flies' units.

    Args:
        flies: the flies' (or wirings') populations; ignored when `ref` is given.
        odors, breadth_edges, parts: as `fly_reference`.
        ref: a reference already built from the flies.

    Returns:
        One row per fly (index the flies' names) with the `Score.as_series` columns.
    """
    if ref is None:
        ref = fly_reference(flies, odors, breadth_edges, parts)
    rows = []
    for i in range(len(ref.flies)):
        others = subset_reference(ref, [j for j in range(len(ref.flies)) if j != i], keep_units=True)
        arrays = {p: df.to_numpy(dtype=float)[i:i + 1] for p, df in ref.parts.items()}
        rows.append(score_arrays(arrays, others).iloc[0])
    return pd.DataFrame(rows, index=ref.flies)


def floor(ref: FlyReference, n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """The floor: resampled sets of the same flies (drawn with replacement) scored against the measured set.

    Returns:
        One row per draw with the `Score.as_series` columns; the floor is their mean (`floor_summary`).
    """
    rng = np.random.default_rng(seed)
    n = len(ref.flies)
    draws = rng.integers(0, n, size=(n_boot, n))
    arrays = {p: _nanmean_rows(df.to_numpy(dtype=float)[draws]) for p, df in ref.parts.items()}
    return score_arrays(arrays, ref)


def floor_summary(draws: pd.DataFrame) -> pd.Series:
    """Mean of every column of `floor`'s draws."""
    return draws.mean()


def wirings_as_flies(wirings: Sequence[Population], ref: FlyReference) -> pd.DataFrame:
    """The model as a population of flies: each wiring's score, their mean, and the score of the wirings' mean views.

    Returns:
        One row per wiring (index the wirings' names), then a row 'mean' (the mean of the wirings' rows) and a row
        'common' (the wirings' mean views scored: the error common to the wirings), with `Score.as_series` columns.
    """
    views = [population_views(w, ref.odors, ref.breadth_edges, list(ref.parts)) for w in wirings]
    names = [w.name or str(i) for i, w in enumerate(wirings)]
    each = score_arrays(stack_views(views, ref), ref, index=names)
    common = score_arrays(stack_views([mean_views(views)], ref), ref, index=['common'])
    mean = each.mean().to_frame('mean').T
    return pd.concat([each, mean, common])


def split_halves(n: int, n_splits: int = N_SPLITS, seed: int = SEED) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Random splits of n flies into two halves (sizes floor(n/2) and the rest), as index arrays."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_splits):
        p = rng.permutation(n)
        out.append((np.sort(p[:n // 2]), np.sort(p[n // 2:])))
    return out


def step_change(upstream: Mapping[str, pd.Series], downstream: Mapping[str, pd.Series]) -> Dict[str, pd.Series]:
    """The change of each part from the upstream to the downstream stage, entry by entry (downstream minus upstream),
    on the entries both have."""
    out = {}
    for p in ALL_PARTS:
        if p in upstream and p in downstream:
            u, d = upstream[p], downstream[p]
            idx = [e for e in d.index if e in set(u.index)]
            out[p] = d.reindex(idx) - u.reindex(idx)
    return out


def step_interval(upstream: Sequence[Population], downstream: Sequence[Population],
        odors: Optional[Sequence[str]] = None, breadth_edges: BinSpec = APL16_BREADTH_EDGES,
        n_boot: int = N_BOOT, ci: float = 0.9, seed: int = SEED,
        parts: Optional[Sequence[str]] = None) -> Dict[str, pd.DataFrame]:
    """The flies' change of each part between two stages, with an interval from resampling the flies of both stages.

    Args:
        upstream, downstream: the two stages' flies.
        odors: the odors read (default those every fly of both stages has).
        breadth_edges: the breadth bins.
        n_boot: resamples (each stage's flies drawn with replacement, independently).
        ci: the interval's coverage.
        seed: of the resampling.
        parts: the parts to read (default those both stages have).

    Returns:
        Per part, a DataFrame indexed by entry with columns 'change' (the measured flies' mean change), 'lo' and 'hi'
        (the `ci` interval of the resampled change).
    """
    rs = step_resample(upstream, downstream, odors, breadth_edges, n_boot, seed, parts)
    q = ((1 - ci) / 2, 1 - (1 - ci) / 2)
    out = {}
    for p, ch in rs.change.items():
        B = rs.resampled[p]
        lo, hi = np.nanquantile(B, q[0], axis=0), np.nanquantile(B, q[1], axis=0)
        out[p] = pd.DataFrame({'change': ch.to_numpy(dtype=float), 'lo': lo, 'hi': hi}, index=ch.index)
    return out


STEP_DISTANCE_Q = 0.95
STEP_DRAWS = ('sets', 'one_fly')


@dataclass
class StepResample:
    """The flies' side of one step (upstream stage -> downstream stage) for the step's resampled distance.

    Attributes:
        change: per part, the measured flies' mean change by entry (downstream minus upstream, `step_change`).
        resampled: per part, n_boot x entries: the change of each resampled set of flies (columns `change`'s entries).
        n_up, n_down: the flies of each stage.
        draw: what one resample is ('sets', 'one_fly'; `step_resample`).
        corrected: whether the resampled deviations were scaled by sqrt(n / (n - 1)) per stage.
        seed: of the resampling.
    """
    change: Dict[str, pd.Series]
    resampled: Dict[str, np.ndarray]
    n_up: int
    n_down: int
    draw: str = 'sets'
    corrected: bool = False
    seed: int = SEED

    def distances(self, part: str) -> np.ndarray:
        """The RMS distance of each resampled change from the measured change over the entries both define."""
        return _rms_rows(self.resampled[part] - self.change[part].to_numpy(dtype=float)[None, :])


def _rms_rows(D: np.ndarray) -> np.ndarray:
    ok = np.isfinite(D)
    with np.errstate(invalid='ignore', divide='ignore'):
        msd = np.where(ok, D ** 2, 0.0).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)
        return np.where(ok.any(axis=1), np.sqrt(msd), np.nan)


def step_resample(upstream: Sequence[Population], downstream: Sequence[Population],
        odors: Optional[Sequence[str]] = None, breadth_edges: BinSpec = APL16_BREADTH_EDGES, n_boot: int = N_BOOT,
        seed: int = SEED, parts: Optional[Sequence[str]] = None, draw: str = 'sets',
        correct: bool = False) -> StepResample:
    """The flies' change of each part between two stages and its resamples (`step_interval`'s draws).

    Args:
        upstream, downstream: the two stages' flies (at least two each).
        odors: the odors read (default those every fly of both stages has).
        breadth_edges: the breadth bins.
        n_boot: resamples.
        seed: of the resampling.
        parts: the parts to read (default those both stages have).
        draw: 'sets' (`step_interval`'s: each stage's flies drawn with replacement, independently, as many as it has,
            the change of their mean parts), or 'one_fly' (one fly drawn from each stage: the change of a single
            animal, the spread of one wiring's change if the model were a fly).
        correct: scale each stage's resampled deviation from its measured mean part by sqrt(n / (n - 1)) ('sets'
            only), which removes the bootstrap's shrinkage of a mean's variance by (n - 1) / n (a quarter of the
            variance for four flies).

    Returns:
        A `StepResample`.
    """
    if draw not in STEP_DRAWS:
        raise ValueError(f'{draw=}')
    if correct and draw != 'sets':
        raise ValueError("the correction applies to draw 'sets' only")
    odors = common_odors(list(upstream) + list(downstream)) if odors is None else list(odors)
    up = fly_reference(upstream, odors, breadth_edges, parts)
    down = fly_reference(downstream, odors, breadth_edges, parts)
    nu, nd = len(up.flies), len(down.flies)
    rng = np.random.default_rng(seed)
    k_u, k_d = (nu, nd) if draw == 'sets' else (1, 1)
    du = rng.integers(0, nu, size=(n_boot, k_u))
    dd = rng.integers(0, nd, size=(n_boot, k_d))
    change, resampled = {}, {}
    for p, ch in step_change(up.mean, down.mean).items():
        U = up.parts[p].reindex(columns=ch.index).to_numpy(dtype=float)
        D = down.parts[p].reindex(columns=ch.index).to_numpy(dtype=float)
        bu, bd = _nanmean_rows(U[du]), _nanmean_rows(D[dd])
        if correct:
            mu = up.mean[p].reindex(ch.index).to_numpy(dtype=float)
            md = down.mean[p].reindex(ch.index).to_numpy(dtype=float)
            bu = mu + (bu - mu) * np.sqrt(nu / (nu - 1))
            bd = md + (bd - md) * np.sqrt(nd / (nd - 1))
        change[p] = ch
        resampled[p] = bd - bu
    return StepResample(change=change, resampled=resampled, n_up=nu, n_down=nd, draw=draw, corrected=bool(correct),
        seed=seed)


def step_null_distance(upstream: Sequence[Population], downstream: Sequence[Population],
        odors: Optional[Sequence[str]] = None, breadth_edges: BinSpec = APL16_BREADTH_EDGES,
        q: float = STEP_DISTANCE_Q, n_boot: int = N_BOOT, seed: int = SEED, parts: Optional[Sequence[str]] = None,
        draw: str = 'sets', correct: bool = False) -> Dict[str, float]:
    """Per part, the threshold of the step class: the `q` quantile of the RMS distance of a resampled change from the
    flies' measured change, over all the entries the measured change defines (the lock's rules -> step_class).

    The draws are `step_resample`'s, the same as `step_interval`'s for the same `n_boot` and `seed`, so the measured
    change is the one `step_interval` gives. A model's change of a part is inside when its RMS distance from that
    change is at or below the threshold; for a model that lacks some of a part's entries, `step_distance` compares
    it with the null on the entries it has.

    Args:
        upstream, downstream: the two stages' flies.
        odors: the odors read (default those every fly of both stages has).
        breadth_edges: the breadth bins.
        q: the quantile (0.95 by the lock).
        n_boot, seed, parts, draw, correct: as `step_resample`.

    Returns:
        {part: threshold} (NaN for a part with no defined entry).
    """
    rs = step_resample(upstream, downstream, odors, breadth_edges, n_boot, seed, parts, draw, correct)
    out = {}
    for p in rs.change:
        d = rs.distances(p)
        out[p] = float(np.quantile(d[np.isfinite(d)], q)) if np.isfinite(d).any() else np.nan
    return out


def step_distance(model_changes: Union[Mapping[str, pd.Series], Mapping[str, Mapping[str, pd.Series]],
        Sequence[Mapping[str, pd.Series]]], flies: StepResample, q: float = STEP_DISTANCE_Q) -> pd.DataFrame:
    """The step's resampled distance (the step class's reading, 2026-11-24): per model (wiring) and part, the RMS
    distance of the model's per-entry change from the flies' change, against the same distance for resampled flies.

    A model's distance is over the entries that it and the flies' change define (the others are counted,
    `n_missing`); each resampled change's distance is taken over the same entries, so a model that lacks entries is
    compared with a null on those it has. "Inside" is the model's distance at or below the `q` quantile of the
    resampled distances. How many wirings must be inside is the caller's (six of eight, the runner's).

    Args:
        model_changes: one model's change (`step_change`'s dict, part -> Series by entry), several as a mapping
            name -> change, or a sequence of changes (named by position).
        flies: `step_resample`'s flies' side of the step.
        q: the quantile of the resampled distances a model's distance is compared with.

    Returns:
        One row per model and part (by part in the flies' order, then the models in the given order): model, part,
        distance, threshold (the `q` quantile of the resampled distances), inside (False where the distance is NaN),
        share_closer (the share of resampled changes closer to the flies' change than the model's), resampled_median,
        n_entries (scored), n_missing.
    """
    if isinstance(model_changes, Mapping) and all(isinstance(v, pd.Series) for v in model_changes.values()):
        models = {'model': model_changes}
    elif isinstance(model_changes, Mapping):
        models = dict(model_changes)
    else:
        models = {str(i): m for i, m in enumerate(model_changes)}
    names = list(models)
    rows = []
    for p, ch in flies.change.items():
        c = ch.to_numpy(dtype=float)
        M = np.array([np.asarray(models[n][p].reindex(ch.index), dtype=float) if p in models[n]
            else np.full(len(c), np.nan) for n in names]).reshape(len(names), len(c))
        dm = M - c[None, :]
        use = np.isfinite(dm)
        dist = _rms_rows(dm)
        D2 = (flies.resampled[p] - c[None, :]) ** 2
        okb = np.isfinite(D2)
        W = use.astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            null = np.sqrt((np.where(okb, D2, 0.0) @ W.T) / (okb.astype(float) @ W.T))   # n_boot x models
        thr = np.array([np.nanquantile(null[:, i], q) if np.isfinite(null[:, i]).any() else np.nan
            for i in range(len(names))])
        med = np.array([np.nanmedian(null[:, i]) if np.isfinite(null[:, i]).any() else np.nan
            for i in range(len(names))])
        fin = np.isfinite(null)
        closer = np.array([np.mean(null[fin[:, i], i] < dist[i]) if np.isfinite(dist[i]) and fin[:, i].any()
            else np.nan for i in range(len(names))])
        for i, n in enumerate(names):
            rows.append(dict(model=n, part=p, distance=dist[i], threshold=thr[i],
                inside=bool(np.isfinite(dist[i]) and dist[i] <= thr[i]), share_closer=closer[i],
                resampled_median=med[i], n_entries=int(use[i].sum()),
                n_missing=int((np.isfinite(c) & ~np.isfinite(M[i])).sum())))
    return pd.DataFrame(rows, columns=['model', 'part', 'distance', 'threshold', 'inside', 'share_closer',
        'resampled_median', 'n_entries', 'n_missing'])


def population_step_change(upstream: Population, downstream: Population, odors: Sequence[str],
        breadth_edges: BinSpec = APL16_BREADTH_EDGES, parts: Optional[Sequence[str]] = None) -> Dict[str, pd.Series]:
    """One model wiring's change of each part from its upstream to its downstream population (`step_change` of their
    views on `odors`), the input of `step_distance`."""
    return step_change(population_views(upstream, odors, breadth_edges, parts),
        population_views(downstream, odors, breadth_edges, parts))


def no_step_score(upstream: Sequence[Population], downstream: FlyReference) -> Score:
    """The no-step yardstick: the upstream flies' mean views as the prediction for the downstream flies.

    `downstream` must be built on odors every upstream fly has (so that both sides' strength parts are multiples of
    the mean over the same odors); parts the upstream stage lacks are left out of the score.
    """
    lack = [o for o in downstream.odors if o not in set(common_odors(upstream))]
    if lack:
        raise ValueError(f'the upstream flies lack odors of the downstream reference: {lack}')
    views = [population_views(u, downstream.odors, downstream.breadth_edges, list(downstream.parts)) for u in upstream]
    m = mean_views(views)
    ref = dataclasses.replace(downstream, parts={p: d for p, d in downstream.parts.items() if p in m},
        mean={p: s for p, s in downstream.mean.items() if p in m},
        one_fly={p: s for p, s in downstream.one_fly.items() if p in m})
    return score_views(m, ref)


def rank_map(values: pd.DataFrame, reference: np.ndarray) -> pd.DataFrame:
    """Maps a population's pooled values rank for rank onto a reference distribution (quantile matching: the value
    of rank r of n goes to the reference's (r - 0.5) / n quantile, ties sharing their mean rank); NaN stays NaN."""
    ref = np.asarray(reference, dtype=float)
    ref = ref[np.isfinite(ref)]
    v = values.to_numpy(dtype=float)
    ok = np.isfinite(v)
    out = np.full(v.shape, np.nan)
    if ok.any():
        r = stats.rankdata(v[ok])
        out[ok] = np.quantile(ref, (r - 0.5) / ok.sum())
    return pd.DataFrame(out, index=values.index, columns=values.columns)


def energy_sample(x: np.ndarray, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> np.ndarray:
    """The units (rows) of one set that `energy_distance` uses: rows with a NaN dropped, then at most `max_units`
    drawn without replacement by a generator of `seed` (so a set's draw depends only on the set and the seed)."""
    x = np.asarray(x, dtype=float)
    x = x.reshape(len(x), -1)
    x = x[np.isfinite(x).all(axis=1)]
    if len(x) > max_units:
        x = x[np.sort(np.random.default_rng(seed).choice(len(x), max_units, replace=False))]
    return x


def energy_from_samples(x: np.ndarray, y: np.ndarray, exx: Optional[float] = None,
        eyy: Optional[float] = None) -> float:
    """The energy distance of two drawn sets (`energy_sample`), with their mean within-set distances E|X-X'| and
    E|Y-Y'| passed in when already known (`mean_self_distance`)."""
    if not len(x) or not len(y):
        return np.nan
    exx = mean_self_distance(x) if exx is None else exx
    eyy = mean_self_distance(y) if eyy is None else eyy
    return float(np.sqrt(max(2 * cdist(x, y).mean() - exx - eyy, 0.0)))


def mean_self_distance(x: np.ndarray) -> float:
    """E|X-X'| over all pairs of rows, each row with itself included (scipy's convention)."""
    return float(cdist(x, x).mean()) if len(x) else np.nan


def energy_distance(a: np.ndarray, b: np.ndarray, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> float:
    """Energy distance between two sets of units' tuning vectors (rows), sqrt(2 E|X-Y| - E|X-X'| - E|Y-Y'|) with
    Euclidean norms and all pairs counted (scipy.stats.energy_distance's convention in one dimension), on at most
    `max_units` units of each set (`energy_sample`; rows with a NaN are dropped)."""
    return energy_from_samples(energy_sample(a, max_units, seed), energy_sample(b, max_units, seed))


def _tuning(pop: Population, odors: Sequence[str], reference: Optional[np.ndarray]) -> np.ndarray:
    v = pop.values.reindex(columns=list(odors))
    if reference is not None:
        v = rank_map(v, reference)
    return v.to_numpy(dtype=float)


def energy_to_each(pop: Population, flies: Sequence[Population], odors: Sequence[str],
        reference: Optional[np.ndarray] = None, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> pd.Series:
    """The population's energy distance to each fly (index the flies' names)."""
    x = _tuning(pop, odors, reference)
    return pd.Series([energy_distance(x, _tuning(f, odors, reference), max_units, seed) for f in flies],
        index=[f.name or str(i) for i, f in enumerate(flies)], dtype=float)


def fly_energy_matrix(flies: Sequence[Population], odors: Sequence[str], reference: Optional[np.ndarray] = None,
        max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> pd.DataFrame:
    """The energy distances between every two flies (symmetric, 0 on the diagonal)."""
    names = [f.name or str(i) for i, f in enumerate(flies)]
    X = [_tuning(f, odors, reference) for f in flies]
    M = np.zeros((len(X), len(X)))
    for i, j in itertools.combinations(range(len(X)), 2):
        M[i, j] = M[j, i] = energy_distance(X[i], X[j], max_units, seed)
    return pd.DataFrame(M, index=names, columns=names)


def energy_guard(pop: Population, flies: Sequence[Population], odors: Optional[Sequence[str]] = None,
        reference: Optional[np.ndarray] = None, max_units: int = ENERGY_MAX_UNITS, seed: int = SEED) -> Dict[str, float]:
    """The energy-distance guard: the population's mean energy distance to the flies, the flies' mean distance to each
    other, and their ratio.

    Args:
        pop: the model's (or one fly's) population.
        flies: the flies' populations (a fly with `pop`'s name is neither compared with it nor counted among the
            flies, so a fly is read against the others).
        odors: the odors read (default those all have).
        reference: if given, every population's pooled values are first mapped rank for rank onto it (`rank_map`).
        max_units: units drawn from each population.
        seed: of the draw.

    Returns:
        Dict with 'to_flies', 'fly_to_fly' and 'ratio'.
    """
    flies = [f for f in flies if not (pop.name and f.name == pop.name)]
    odors = common_odors([pop, *flies]) if odors is None else list(odors)
    to = energy_to_each(pop, flies, odors, reference, max_units, seed)
    M = fly_energy_matrix(flies, odors, reference, max_units, seed).to_numpy()
    ff = float(M[np.triu_indices(len(M), k=1)].mean()) if len(M) > 1 else np.nan
    return dict(to_flies=float(to.mean()), fly_to_fly=ff, ratio=float(to.mean()) / ff if ff else np.nan)
