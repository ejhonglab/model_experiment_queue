"""E142's step 0 (no simulation): the flies' step by odor, every saved run's entries by odor, the count-weighted
yardsticks, what the model's bouton reading samples, the equal-weight reading of the boutons and the classes read again,
the planner's probe numbers reproduced, and step 0's classes. The definitions, the saved runs, the halves, the seeds and
the thresholds are those of experiments/E142_dendrite_bouton_step.py and its committed lock (experiments/E142_lock.json);
`run` calls `E142.require_all` first, and every piece reads the context it builds (`Context`).

Pieces (`PIECES`, run in this order by `run`; a later piece reads the earlier pieces' files):
- flies: flies_steps.csv (each step's change by odor and part on all flies and each half, with its interval),
  flies_steps_summary.csv (per odor and part: outside zero on all flies, on each half, on both halves),
  no_step_entries.csv (the no-step yardsticks' entries by odor: the dendrites as the prediction for the boutons and for
  the claws, the boutons for the claws; against all flies and each half; row 'common', the upstream flies' mean views).
- models: models_entries.csv (every saved run's entries by odor at the dendrites, boutons and claws, each wiring and
  'common', on every fly set), models_concentrated.csv (E138's class per group, stage), models_step.csv (the model's
  step by odor against the flies', r over the 17 odors), models_equal_views.csv (each run's boutons with each
  glomerulus one unit: the strength parts), models_read_checks.csv (the glomerulus labels of the read boutons checked
  against the values before the curve), models_missing.json.
- yardsticks: counts.csv (each glomerulus's boutons and claws in each calyx, beside E138's feeders.csv), carried.csv
  (the carried glomeruli's values), yardsticks.csv (the count-weighted, no-step and permuted-null errors),
  yardstick_entries.csv (the count-weighted yardsticks' entries by odor, all flies), yardstick_null.csv.
- sampling: bouton_sampling.csv (per wiring and glomerulus: the wiring's boutons, the reading's units, each calyx's
  counts), bouton_sampling_summary.csv.
- equal: equal_weight_groups.csv (the flies' groups per fly and held-out repeat), equal_weight_groups_check.csv (the
  groups on the trial means against E119's A4), equal_weight_errors.csv (the flies' floor and leave-one-out, and every
  saved run under the equal-weight reading, beside the default reading), reclasses.json (E136's and E139's classes read
  again with the equal-weight reading and over the three other stages).
- probe: probe.csv (the planner's numbers against step 0's reproductions and the saved tables they came from).
- classes: classes_step0.json, tables.md.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E142_dendrite_bouton_step as E142
from kcq import harness as H

PIECES: Tuple[str, ...] = ('flies', 'models', 'yardsticks', 'sampling', 'equal', 'probe', 'classes')
PARTS = list(E142.STRENGTH_PARTS)
PART_KEY = {'strength/rate': 'rate', 'strength/mean': 'mean'}


# ---------------------------------------------------------------------------
# small helpers

def _write_csv(df: pd.DataFrame, name: str, out: Optional[Path] = None) -> Path:
    p = (E142.out_dir() if out is None else Path(out)) / name
    tmp = p.with_name(f'.{p.name}.{E142._host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, p)
    return p


def _write_json(name: str, obj, out: Optional[Path] = None) -> Path:
    p = (E142.out_dir() if out is None else Path(out)) / name
    E142._e140()._write_json(p, obj)
    return p


def _read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(E142.out_dir() / name)


def _tag(df: pd.DataFrame, **cols) -> pd.DataFrame:
    df = df.copy()
    for i, (k, v) in enumerate(cols.items()):
        df.insert(i, k, v)
    return df


def _sd(x) -> float:
    x = np.asarray(list(x), dtype=float)
    return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan


def _rms(xs) -> float:
    xs = np.asarray(list(xs), dtype=float)
    return float(np.sqrt(np.mean(xs ** 2))) if len(xs) else np.nan


def provenance(piece: str, ctx: 'Context', **extra) -> None:
    """provenance_step0_<piece>.json: the script, `H.code_commits()`, host, time, the lock's sha256 and `extra`."""
    _write_json(f'provenance_step0_{piece}.json',
                dict(piece=piece, script='experiments/E142_step0.py', main=E142.SCRIPT, **H.code_commits(),
                     host=E142._host(), time=E142._now(), lock_sha256=ctx.lock.get('sha256'),
                     e134_lock_sha256=ctx.e134.get('sha256'), **extra))


# ---------------------------------------------------------------------------
# the context: the flies of every set, read once

class Context:
    """What every piece reads: the lock and E134's, the odors, the flies of each set (E136's `Flies` on the strength
    parts, all flies and each half), the flies' populations per stage (E134's `fly_side_read`), the halves' names, and
    the soma flies only through E136's saved scores."""

    def __init__(self, lock: Mapping, e134: Mapping):
        from experiments import E134_population_baseline as B
        from experiments import E136_stage_population_fit as E136
        E136.set_measure(E142.E136_MEASURE)
        self.lock, self.e134 = dict(lock), dict(e134)
        self.odors = list(lock['odors'])
        stages = tuple(E142.ENTRY_STAGES)
        self.flies = {s: E136.Flies(None if s == 'all' else s, stages, views=('strength',)) for s in E142.FLY_SETS}
        full = self.flies['all'].full
        self.pops = {st: B.fly_side_read(full[st])[0] for st in stages}
        self.names = {st: self.flies['all'].names[st] for st in stages}
        for st in stages:
            want = list(lock['recordings'][st]['flies'])
            got = [p.name for p in self.pops[st]]
            if got != want:
                raise RuntimeError(f"{st}: flies {got} differ from the lock's {want}")
            if self.names[st] != lock['recordings'][st]['halves']:
                raise RuntimeError(f"{st}: halves {self.names[st]} differ from the lock's")
        self.criterion = {st: float(full[st].criterion) for st in stages}
        self.edges = B.breadth_edges(len(self.odors))

    def set_pops(self, stage: str, fly_set: str) -> List:
        """The flies' populations of a stage in one fly set."""
        if fly_set == 'all':
            return list(self.pops[stage])
        keep = set(self.names[stage][fly_set])
        return [p for p in self.pops[stage] if p.name in keep]

    def reference(self, stage: str, fly_set: str, odors: Optional[Sequence[str]] = None):
        """E136's reference of a stage for a fly set (a half's mean views in all the stage's flies' one-fly units), on
        the strength parts."""
        return self.flies[fly_set].reference(stage, list(odors or self.odors))


# ---------------------------------------------------------------------------
# pure functions

def weighted_arrays(W: np.ndarray, values: np.ndarray, responds: np.ndarray) -> Dict[str, np.ndarray]:
    """The strength parts of one named population for many weightings of its units at once.

    Args:
        W: weightings x units (each row one weighting; 0 leaves a unit out).
        values: units x odors (NaN where a unit was not given an odor).
        responds: units x odors of 1 / 0 (NaN where `values` is NaN).

    Returns:
        {'strength/rate', 'strength/mean'}: weightings x odors, each row relative to its mean over the odors (NaN where
        that mean is not positive), as hal's `population_views` gives for the population with each unit repeated its
        weight's number of times.
    """
    W = np.atleast_2d(np.asarray(W, dtype=float))
    V = np.asarray(values, dtype=float)
    Rr = np.asarray(responds, dtype=float)
    M = np.isfinite(V).astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        den = W @ M
        raw = {'strength/rate': (W @ np.nan_to_num(Rr * M)) / den, 'strength/mean': (W @ np.nan_to_num(V)) / den}
        out = {}
        for p, A in raw.items():
            m = np.nanmean(np.where(np.isfinite(A), A, np.nan), axis=1, keepdims=True)
            out[p] = np.where(np.isfinite(m) & (m > 0), A / m, np.nan)
    return out


def weighted_strength(pop, weights: Mapping[str, float], odors: Sequence[str]) -> Dict[str, pd.Series]:
    """The strength parts of a named population with each unit (glomerulus) weighted by `weights` (0 where absent):
    equal to hal's `population_views` of the population with each row repeated its weight's number of times."""
    od = list(odors)
    units = [str(u) for u in pop.values.index]
    w = np.array([float(weights.get(u, 0.0)) for u in units])
    a = weighted_arrays(w[None, :], pop.values.reindex(columns=od).to_numpy(dtype=float),
                        pop.responds.reindex(columns=od).to_numpy(dtype=float))
    return {p: pd.Series(A[0], index=od) for p, A in a.items()}


def repeated_population(pop, weights: Mapping[str, float]):
    """The named population with each unit's row repeated its (integer) weight's number of times (0 drops it)."""
    from kcq import population_fit as P
    idx = [u for u in pop.values.index for _ in range(int(round(float(weights.get(str(u), 0.0)))))]
    V = pop.values.loc[idx].reset_index(drop=True)
    Rr = pop.responds.loc[idx].reset_index(drop=True)
    return P.Population(V, Rr, name=pop.name, named_units=False, graded=pop.graded)


def fly_matrices(pops: Sequence, glomeruli: Sequence[str], odors: Sequence[str],
                 extra: Optional[pd.DataFrame] = None, criterion: Optional[float] = None
                 ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Per fly, (values, responds) on `glomeruli` x `odors` (NaN where the fly lacks a glomerulus); `extra` rows
    (glomeruli x odors, the carried glomeruli) set in every fly where the fly lacks them, responding above
    `criterion`."""
    out = []
    for p in pops:
        V = p.values.copy()
        V.index = V.index.astype(str)
        Rr = p.responds.copy()
        Rr.index = Rr.index.astype(str)
        V = V.reindex(index=list(glomeruli), columns=list(odors))
        Rr = Rr.reindex(index=list(glomeruli), columns=list(odors))
        if extra is not None:
            for g in extra.index:
                if g in V.index and V.loc[g].isna().all():
                    v = extra.loc[g].reindex(list(odors)).astype(float)
                    V.loc[g] = v
                    Rr.loc[g] = (v > float(criterion)).astype(float).where(v.notna())
        out.append((V.to_numpy(dtype=float), Rr.to_numpy(dtype=float)))
    return out


def yardstick_scores(W: np.ndarray, mats: Sequence[Tuple[np.ndarray, np.ndarray]], ref) -> pd.DataFrame:
    """The count-weighted yardsticks of many weightings: per weighting, the flies' (`mats`) mean weighted strength parts
    scored against `ref` (hal's `score_arrays`: one row per weighting, columns view_strength, part_strength_rate,
    part_strength_mean, ...)."""
    from kcq import population_fit as P
    W = np.atleast_2d(np.asarray(W, dtype=float))
    per = [weighted_arrays(W, V, Rr) for V, Rr in mats]
    arrays = {}
    for p in PARTS:
        S = np.stack([x[p] for x in per])                # flies x weightings x odors
        arrays[p] = P._nanmean_rows(np.moveaxis(S, 0, 1))  # weightings x odors (NaN-aware mean over flies)
    return P.score_arrays(arrays, ref)


def permutations(n: int, n_perm: int = E142.N_PERM, seed: int = E142.PERM_SEED) -> np.ndarray:
    """N_PERM permutations of n positions (rows), from `np.random.default_rng(seed)`."""
    rng = np.random.default_rng(seed)
    return np.array([rng.permutation(n) for _ in range(n_perm)])


def null_summary(actual: float, null: np.ndarray, q: float = E142.NULL_PERCENTILE) -> dict:
    """The null's q-th percentile, the share of permuted errors below the actual one (its percentile) and whether the
    actual error is below the q-th percentile."""
    x = np.asarray(null, dtype=float)
    x = x[np.isfinite(x)]
    p = float(np.percentile(x, q)) if len(x) else np.nan
    return dict(null_p05=p, null_median=float(np.median(x)) if len(x) else np.nan,
                null_percentile_of_actual=float(100.0 * np.mean(x < actual)) if len(x) else np.nan,
                below_p05=bool(np.isfinite(p) and actual < p), n_null=int(len(x)))


def outside_zero(lo: float, hi: float) -> int:
    """+1 when an interval is above 0, -1 when below, 0 when it holds 0 (or is undefined)."""
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return 0
    return 1 if lo > 0 else (-1 if hi < 0 else 0)


def bouton_region_glomeruli(run, odors: Sequence[str], condition: str = 'on', reading=None,
                            subset_glomeruli: Optional[Sequence[str]] = None) -> Tuple[np.ndarray, pd.DataFrame]:
    """The run's own glomerulus label of each bouton region its default bouton reading keeps (pool '1': each bouton
    one region; the subset's boutons, the silent ones dropped as E90's `_drop_silent` drops them), in the regions'
    order, with those boutons' values before the curve (for the check against `ModelReading.before_curve`). No region
    is given a glomerulus from the flies: the labels are the model's wiring's."""
    from experiments import E90_model_stages as E90
    from kcq import stage_readings as R
    reading = reading or R.DEFAULT_READINGS['pn_bouton']
    if str(reading.pool) != '1':
        raise ValueError(f'the labels follow pool 1 only, not {reading.pool!r}')
    X = R.stage_changes(run, 'pn_bouton', condition, list(odors), reading)
    g = np.asarray(run.meta['bouton_glomerulus']).astype(str)
    if reading.subset == 'all':
        m = np.ones(len(g), dtype=bool)
    else:
        sub = [str(x) for x in (subset_glomeruli if subset_glomeruli is not None else E90.shared_glomeruli())]
        m = np.isin(g, sub)
    Xm = X.iloc[np.flatnonzero(m)]
    keep = (Xm.abs() > 1e-12).any(axis=1).to_numpy()
    return g[m][keep], Xm.iloc[np.flatnonzero(keep)]


def glomerulus_units(pop, labels: Sequence[str], criterion: float):
    """A read bouton population with each glomerulus's regions averaged into one unit (the model's equal-weight
    reading), responding above `criterion`."""
    from kcq import population_fit as P
    lab = np.asarray(labels).astype(str)
    if len(lab) != len(pop.values):
        raise ValueError(f'{len(lab)} labels for {len(pop.values)} regions')
    V = pop.values.groupby(lab, sort=True).mean()
    return P.Population.from_criterion(V, criterion, name=pop.name)


def fly_group_units(trials: pd.DataFrame, odors: Sequence[str], group_odors: Sequence[str], held_out,
                    criterion: float = E142.GROUP_CRITERION, cut: float = E142.GROUP_CUT):
    """One bouton fly's equal-weight population for one held-out repeat: E128's tuning groups (`fly_groups`) made on
    the mean of the other repeats restricted to `group_odors`, each group (and the remainder, label 0) one unit valued
    at its ROIs' mean on the held-out repeat on `odors`, responding above `criterion`. Returns (population, labels)."""
    from experiments import E128_bouton_groups as E128
    from kcq import population_fit as P
    reps = sorted(set(trials.columns.get_level_values('repeat')))
    T = {r: trials.xs(r, axis=1, level='repeat').reindex(columns=list(odors)).astype(float) for r in reps}
    others = [r for r in reps if r != held_out]
    M2 = sum(T[r] for r in others) / len(others)
    lab = E128.fly_groups(M2[list(group_odors)].to_numpy(dtype=float), criterion, cut)
    V = T[held_out].groupby(lab, sort=True).mean()
    V.index = [f'group{int(i)}' for i in V.index]
    return P.Population.from_criterion(V, criterion, name=str(held_out)), lab


def reference_from_views(views: Mapping[str, Mapping[str, pd.Series]], odors: Sequence[str], edges):
    """A flies' reference (hal's FlyReference) from each fly's strength parts (fly -> part -> Series by odor), its
    one-fly units those of these flies."""
    from kcq import population_fit as P
    names = list(views)
    parts = {p: pd.DataFrame([views[n][p].reindex(list(odors)).to_numpy(dtype=float) for n in names], index=names,
                             columns=list(odors)) for p in PARTS}
    return P._reference_from_parts(parts, names, list(odors), edges)


def reclass(orig: Mapping, arm: Mapping[str, Mapping[str, float]], fam: Mapping[str, Mapping[str, float]],
            spreads: Mapping[str, float], stages: Sequence[str], control_fall: Optional[float],
            drop_guard: Sequence[str] = ()) -> dict:
    """An arm's class read again on new megamat values (E136's `classify` decision, on E136's `compare`): `arm` and
    `fam` per wiring {'stage': ..., <stage>: ...}; the held-out panels, the all-views guard, the realism, the bounds and
    the other reasons kept from `orig` (the arm's original record), the permuted control's test made again with the
    original control fall; `drop_guard` the stages left out of the guard (the boutons, over the three other stages)."""
    from experiments import E136_stage_population_fit as E136
    m = E136.compare({w: v['stage'] for w, v in arm.items()}, {w: v['stage'] for w, v in fam.items()}, spreads['stage'])
    st = {s: E136.compare({w: v[s] for w, v in arm.items()}, {w: v[s] for w, v in fam.items()}, spreads[s])
          for s in stages}
    reasons = [r for r in (orig.get('reasons') or []) if not str(r).startswith('no better than the permuted control')]
    if not m.get('complete') or not all(v.get('complete') for v in st.values()):
        reasons.append('reads missing')
    fall = -m['change'] if m.get('complete') else np.nan
    if control_fall is not None and np.isfinite(fall) and not fall > control_fall:
        reasons.append('no better than the permuted control')
    stage_worse = [s for s, v in st.items() if v.get('worse_than_spread')]
    guard_worse = [s for s in (orig.get('guard_worse') or []) if s not in set(drop_guard)]
    panel_worse = list(orig.get('panel_worse') or [])
    panel_missing = list(orig.get('panel_missing') or [])
    improves = bool(m.get('improves')) and not stage_worse and not guard_worse
    if reasons:
        cls = 'no class'
    elif improves and not panel_worse and not panel_missing:
        cls = 'holds'
    elif improves and panel_worse:
        cls = 'megamat only'
    elif improves:
        cls = 'improves (held-out panels not all read)'
    elif m.get('improves') and stage_worse:
        cls = "one stage at another's cost"
    elif m.get('improves') and guard_worse:
        cls = "strength at another view's cost"
    else:
        cls = 'does not improve'
    return dict(cls=cls, original_cls=orig.get('cls'), megamat=m, stages=st, fall=fall, control_fall=control_fall,
                improves=improves, stage_worse=stage_worse, guard_worse=guard_worse, panel_worse=panel_worse,
                panel_missing=panel_missing, reasons=reasons)


def readings_class(change: float, spread: float) -> dict:
    """"The reading's": A1's change at the boutons under the equal-weight reading is inside the spread or of the other
    sign (change >= -spread)."""
    ok = bool(np.isfinite(change) and np.isfinite(spread) and change >= -spread)
    return dict(change=float(change), spread=float(spread), inside_spread=bool(abs(change) <= spread),
                other_sign=bool(change > 0), holds=ok)


def counts_class(rows: pd.DataFrame, calyx: str, stage: str = 'pn_bouton', glomeruli: str = 'dendrites') -> dict:
    """"The counts make the step" for one calyx from yardsticks.csv rows: the count-weighted yardstick below the
    no-step yardstick by COUNTS_MARGIN or more on both halves, and below the null's 5th percentile on all flies."""
    x = rows[(rows.calyx == calyx) & (rows.stage == stage) & (rows.glomeruli == glomeruli)].set_index('fly_set')
    diffs = {h: float(x.loc[h, 'no_step'] - x.loc[h, 'counts']) for h in E142.HALVES if h in x.index}
    margin = {h: bool(d >= E142.COUNTS_MARGIN) for h, d in diffs.items()}
    below = bool(x.loc['all', 'below_p05']) if 'all' in x.index else False
    holds = bool(len(margin) == len(E142.HALVES) and all(margin.values()) and below)
    return dict(calyx=calyx, stage=stage, glomeruli=glomeruli,
                counts={s: float(x.loc[s, 'counts']) for s in x.index},
                no_step={s: float(x.loc[s, 'no_step']) for s in x.index}, below_no_step_by=diffs,
                margin_met=margin, null_p05=float(x.loc['all', 'null_p05']) if 'all' in x.index else np.nan,
                null_percentile_of_actual=float(x.loc['all', 'null_percentile_of_actual']) if 'all' in x.index
                else np.nan, below_p05=below, holds=holds)


# ---------------------------------------------------------------------------
# piece 1: the flies

def flies_steps(ctx: Context) -> Dict[str, pd.DataFrame]:
    """The flies' step by odor (STEPS) on every fly set, with intervals (hal's `step_interval`); per odor and part
    whether it is outside zero on all flies, on each half and on both halves; the no-step yardsticks' entries by odor
    (E138's `odor_entries` of the upstream flies against the downstream reference, row 'common')."""
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import population_fit as P
    t0 = time.time()
    rows = []
    for up, down in E142.STEPS:
        for s in E142.FLY_SETS:
            U, D_ = ctx.set_pops(up, s), ctx.set_pops(down, s)
            iv = P.step_interval(U, D_, ctx.odors, ctx.edges, E142.STEP_N_BOOT, E142.STEP_CI, E142.STEP_SEED,
                                 parts=PARTS)
            mu = P.fly_reference(U, ctx.odors, ctx.edges, PARTS).mean
            md = P.fly_reference(D_, ctx.odors, ctx.edges, PARTS).mean
            for p, df in iv.items():
                for o, r in df.iterrows():
                    rows.append(dict(step=f'{up}->{down}', upstream=up, downstream=down, fly_set=s,
                                     n_up=len(U), n_down=len(D_), part=PART_KEY[p], odor=o,
                                     up=float(mu[p][o]), down=float(md[p][o]), change=float(r.change),
                                     lo=float(r.lo), hi=float(r.hi), outside=outside_zero(r.lo, r.hi)))
    steps = pd.DataFrame(rows)
    summ = []
    for (step, part, odor), x in steps.groupby(['step', 'part', 'odor'], sort=False):
        o = x.set_index('fly_set')
        a, f, h = (int(o.loc[k, 'outside']) for k in ('all', 'fit', 'held'))
        summ.append(dict(step=step, part=part, odor=odor, change_all=float(o.loc['all', 'change']),
                         lo_all=float(o.loc['all', 'lo']), hi_all=float(o.loc['all', 'hi']),
                         change_fit=float(o.loc['fit', 'change']), change_held=float(o.loc['held', 'change']),
                         outside_all=a, outside_fit=f, outside_held=h,
                         outside_both_halves=int(f if (f != 0 and f == h) else 0)))
    summary = pd.DataFrame(summ)
    # the no-step yardsticks' entries by odor
    ent = []
    for up, down in (('pn_dendrite', 'pn_bouton'), ('pn_dendrite', 'kc_claw'), ('pn_bouton', 'kc_claw')):
        for s in E142.FLY_SETS:
            ref = ctx.reference(down, s)
            e = E138.odor_entries(ctx.pops[up], ref, names=[p.name for p in ctx.pops[up]])
            ent.append(_tag(e[e.row == 'common'], upstream=up, downstream=down, fly_set=s))
    no_step = pd.concat(ent, ignore_index=True)
    _write_csv(steps, 'flies_steps.csv')
    _write_csv(summary, 'flies_steps_summary.csv')
    _write_csv(no_step, 'no_step_entries.csv')
    provenance('flies', ctx, seconds=round(time.time() - t0, 1))
    print(f'flies: done ({time.time() - t0:.0f} s)', flush=True)
    return dict(steps=steps, summary=summary, no_step=no_step)


# ---------------------------------------------------------------------------
# piece 2: the saved runs

def run_dirs(source: str) -> Tuple[List[Path], str]:
    """(the folders a group's runs are read from, the tag pattern) of a MODEL_GROUPS source."""
    from kcq import stage_readings as R
    exp, pattern = source.split(':', 1)
    if exp == 'E134':
        return R.step3_dirs(), pattern
    if exp == 'E136':
        from experiments import E136_stage_population_fit as E136
        E136.set_measure(E142.E136_MEASURE)
        return E136.read_dirs(), pattern
    if exp == 'E139':
        return R._archive_dirs(E142.E139_EXPERIMENT, 'read'), pattern
    if exp == 'E135':
        from experiments import E135_antennal_lobe_dendrites as E135
        return E135.step2_dirs(), pattern
    raise ValueError(source)


def group_runs(g) -> Tuple[Dict[str, object], List[str]]:
    """A model group's saved runs ({wiring: SavedRun}, in WIRINGS order) and the tags not found."""
    from kcq import stage_readings as R
    dirs, pattern = run_dirs(g.source)
    runs, missing = {}, []
    for w in g.wirings:
        tag = pattern.format(w=w)
        try:
            runs[w] = R.load_step3_run(tag, dirs=dirs)
        except FileNotFoundError:
            missing.append(tag)
    return runs, missing


def read_group(runs: Mapping[str, object], ctx: Context, fly_set: str) -> Dict[str, object]:
    """Each run read as E136 reads a wiring against a fly set (atlas's `read_run`, default readings, ENTRY_STAGES, the
    wiring's seed index)."""
    from kcq import stage_readings as R
    return {w: R.read_run(r, stages=list(E142.ENTRY_STAGES), lock=ctx.e134, flies=ctx.flies[fly_set].cache,
                          seed_index=E142.WIRINGS.index(w)) for w, r in runs.items()}


def models_entries(ctx: Context) -> Dict[str, pd.DataFrame]:
    """Every saved run's entries by odor at ENTRY_STAGES on every fly set (E138's `odor_entries`), E138's
    "concentrated" per group and stage, the model's step against the flies' (r over the odors), and each run's
    equal-weight bouton views (each glomerulus's read boutons one unit), with the labels checked."""
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import population_fit as P
    t0 = time.time()
    ent, eqv, checks, missing, files = [], [], [], {}, {}
    for g in E142.MODEL_GROUPS:
        runs, miss = group_runs(g)
        if miss:
            missing[g.name] = miss
        if not runs:
            print(f'models: {g.name}: no runs', flush=True)
            continue
        files[g.name] = {w: [dict(path=f['path'], sha256=f.get('sha256')) for f in r.files] for w, r in runs.items()}
        for s in E142.FLY_SETS:
            reads = read_group(runs, ctx, s)
            for st in E142.ENTRY_STAGES:
                pops = [reads[w].populations[(st, E142.PANEL, 'on')] for w in runs]
                ref = ctx.reference(st, s, list(pops[0].values.columns))
                e = E138.odor_entries(pops, ref, names=list(runs))
                ent.append(_tag(e, group=g.name, fly_set=s, stage=st))
            key = ('pn_bouton', E142.PANEL, 'on')
            for w, run in runs.items():
                lab, Xk = bouton_region_glomeruli(run, ctx.odors)
                pop = reads[w].populations[key]
                bc = reads[w].before_curve[key]
                same = bool(len(lab) == len(pop.values) == len(bc)
                            and np.allclose(Xk.to_numpy(dtype=float), bc.to_numpy(dtype=float), rtol=0, atol=1e-12))
                checks.append(dict(group=g.name, fly_set=s, wiring=w, n_regions=int(len(pop.values)),
                                   n_labels=int(len(lab)), labels_match_values=same))
                if not same:
                    raise RuntimeError(f'{g.name} {w} {s}: the bouton labels do not follow the read regions')
                eq = glomerulus_units(pop, lab, ctx.criterion['pn_bouton'])
                v = P.population_views(eq, ctx.odors, ctx.edges, PARTS)
                for p in PARTS:
                    eqv += [dict(group=g.name, fly_set=s, wiring=w, part=PART_KEY[p], odor=o, value=float(x),
                                 n_units=int(len(eq.values))) for o, x in v[p].items()]
            del reads
        print(f'models: {g.name} ({len(runs)} wirings, {time.time() - t0:.0f} s)', flush=True)
    entries = pd.concat(ent, ignore_index=True)
    # E138's concentrated, per group and stage, on the row 'common'
    conc = []
    for (gname, st), x in entries.groupby(['group', 'stage'], sort=False):
        by = {s: y for s, y in x.groupby('fly_set')}
        if not all(s in by for s in E142.FLY_SETS):
            continue
        c = E138.concentrated(by['all'], {h: by[h] for h in E142.HALVES}, row='common', top=E142.CONCENTRATED_TOP,
                              share=E142.CONCENTRATED_SHARE)
        conc.append(dict(group=gname, stage=st, cls=c['cls'], odors=';'.join(c['odors']),
                         share_all=c['share_all'], share_fit=c['share_halves'].get('fit'),
                         share_held=c['share_halves'].get('held'),
                         own_fit=';'.join(c['halves_own']['fit']['odors']),
                         own_held=';'.join(c['halves_own']['held']['odors'])))
    # the model's step by odor against the flies' (all flies; each row and 'common')
    from scipy import stats
    steps = []
    a = entries[entries.fly_set == 'all']
    for gname, x in a.groupby('group', sort=False):
        for up, down in E142.STEPS:
            for row in x.row.unique():
                u = x[(x.stage == up) & (x.row == row)].set_index('odor')
                d = x[(x.stage == down) & (x.row == row)].set_index('odor')
                if u.empty or d.empty:
                    continue
                for k in ('rate', 'mean'):
                    ms = (d[f'{k}_model'] - u[f'{k}_model']).reindex(ctx.odors)
                    fs = (d[f'{k}_flies'] - u[f'{k}_flies']).reindex(ctx.odors)
                    ok = ms.notna() & fs.notna()
                    steps.append(dict(group=gname, step=f'{up}->{down}', row=row, part=k, n_odors=int(ok.sum()),
                                      r=float(np.corrcoef(ms[ok], fs[ok])[0, 1]) if ok.sum() > 2 else np.nan,
                                      rho=float(stats.spearmanr(ms[ok], fs[ok]).statistic) if ok.sum() > 2 else np.nan,
                                      model_step_sd=float(ms[ok].std(ddof=0)), flies_step_sd=float(fs[ok].std(ddof=0)),
                                      rms_diff=float(np.sqrt(np.mean((ms[ok] - fs[ok]) ** 2)))))
    _write_csv(entries, 'models_entries.csv')
    _write_csv(pd.DataFrame(conc), 'models_concentrated.csv')
    _write_csv(pd.DataFrame(steps), 'models_step.csv')
    _write_csv(pd.DataFrame(eqv), 'models_equal_views.csv')
    _write_csv(pd.DataFrame(checks), 'models_read_checks.csv')
    _write_json('models_missing.json', missing)
    provenance('models', ctx, groups=files, missing=missing, seconds=round(time.time() - t0, 1))
    print(f'models: done ({time.time() - t0:.0f} s); missing {missing}', flush=True)
    return dict(entries=entries, missing=missing)


# ---------------------------------------------------------------------------
# piece 3: the count-weighted yardsticks

def connectome_counts(lock: Mapping) -> Tuple[pd.DataFrame, dict]:
    """Each glomerulus's boutons and claws in each calyx (E138's `all_counts`), long, with the check against E138's
    feeders.csv (every glomerulus and calyx it lists)."""
    from experiments import E138_strength_by_odor_entries as E138
    c = E138.all_counts()
    rows = []
    for cal, df in c.items():
        for gl, r in df.iterrows():
            rows.append(dict(calyx=cal, glomerulus=str(gl), n_boutons=int(r.n_boutons), n_claws=int(r.n_claws),
                             share_boutons=float(r.share_boutons), share_claws=float(r.share_claws)))
    counts = pd.DataFrame(rows)
    fe = pd.read_csv(lock['inputs']['results']['e138_feeders']['path'])
    fe = fe.drop_duplicates('glomerulus').set_index('glomerulus')
    diffs = []
    for cal in E142.CALYCES:
        x = counts[counts.calyx == cal].set_index('glomerulus')
        for k in ('n_boutons', 'n_claws'):
            col = f'{cal}__{k}'
            if col not in fe:
                diffs.append(f'{col} not in feeders.csv')
                continue
            f = fe[col].dropna()
            mine = x[k].reindex(f.index).fillna(0)
            bad = [g for g in f.index if int(mine[g]) != int(f[g])]
            if bad:
                diffs.append(f'{col}: {bad}')
    check = dict(feeders_glomeruli=int(len(fe)), differences=diffs, agrees=not diffs)
    return counts, check


def carried_rows(ctx: Context) -> Tuple[pd.DataFrame, dict]:
    """The carried glomeruli (the lock's: ORN-recorded, not dendrite-recorded): the median over E135's 28 of its
    measured slope (units 'noise', all flies) times the ORN flies' mean response (trial-noise SD units, over the flies
    that have it), glomeruli x odors."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    ms = pd.read_csv(ctx.lock['inputs']['data']['e135_measured_step']['path'])
    ms = ms[(ms.units == E142.E135_GAIN_UNITS) & (ms.flies == E142.E135_GAIN_FLIES)].set_index('glomerulus')
    g28 = list(ctx.lock['glomeruli']['dendrite28'])
    gain = float(ms.slope.reindex(g28).median())
    orn = P.fly_populations(D.load_population('orn', E142.PANEL, 'wild_type'))
    car = list(ctx.lock['glomeruli']['carried'])
    tabs = []
    for o in orn:
        V = o.values.copy()
        V.index = V.index.astype(str)
        tabs.append(V.reindex(index=car, columns=ctx.odors).to_numpy(dtype=float))
    stack = np.stack(tabs)
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(stack, axis=0)
    out = pd.DataFrame(gain * mean, index=car, columns=ctx.odors)
    info = dict(median_slope=gain, n_slopes=int(ms.slope.reindex(g28).notna().sum()),
                orn_flies=[o.name for o in orn], n_orn_flies_per_glomerulus={g: int(np.isfinite(stack[:, i, 0]).sum())
                                                                              for i, g in enumerate(car)})
    return out, info


def yardsticks(ctx: Context) -> Dict[str, pd.DataFrame]:
    """The counts, the count-weighted yardsticks (the 29 dendrite glomeruli, and with the carried glomeruli beside)
    against the bouton and the claw flies on every fly set, the no-step yardsticks beside, and the null of permuted
    counts (on every set; the class reads all flies')."""
    from experiments import E138_strength_by_odor_entries as E138
    from kcq import population_fit as P
    t0 = time.time()
    counts, check = connectome_counts(ctx.lock)
    extra, cinfo = carried_rows(ctx)
    g29 = list(ctx.lock['glomeruli']['dendrites'])
    g38 = g29 + list(ctx.lock['glomeruli']['carried'])
    sets = {'dendrites': (g29, None), 'with_carried': (g38, extra)}
    perms = {k: permutations(len(v[0])) for k, v in sets.items()}
    crit = ctx.criterion['pn_dendrite']
    mats = {k: fly_matrices(ctx.pops['pn_dendrite'], gl, ctx.odors, ex, crit) for k, (gl, ex) in sets.items()}
    rows, nulls, ents = [], [], []
    for stage, cnt in E142.COUNT_OF.items():
        for s in E142.FLY_SETS:
            ref = ctx.reference(stage, s)
            bouton_no_step = (P.no_step_score(ctx.pops['pn_bouton'], ref).views['strength']
                              if stage == 'kc_claw' else np.nan)
            for gname, (gl, ex) in sets.items():
                ones = np.ones((1, len(gl)))
                ns = yardstick_scores(ones, mats[gname], ref).iloc[0]
                if gname == 'dendrites':
                    chk = P.no_step_score(ctx.pops['pn_dendrite'], ref).views['strength']
                    if not np.isclose(chk, ns.view_strength, rtol=0, atol=1e-10):
                        raise RuntimeError(f'the unweighted yardstick {ns.view_strength} is not the no-step {chk}')
                for cal in E142.CALYCES:
                    w = counts[counts.calyx == cal].set_index('glomerulus')[cnt].reindex(gl).fillna(0).to_numpy(float)
                    sc = yardstick_scores(w[None, :], mats[gname], ref).iloc[0]
                    nl = yardstick_scores(w[perms[gname]], mats[gname], ref)
                    su = null_summary(float(sc.view_strength), nl.view_strength.to_numpy())
                    rows.append(dict(stage=stage, count=cnt, glomeruli=gname, n_glomeruli=len(gl), calyx=cal,
                                     fly_set=s, n_flies=len(ref.flies), counts=float(sc.view_strength),
                                     counts_rate=float(sc.part_strength_rate), counts_mean=float(sc.part_strength_mean),
                                     no_step=float(ns.view_strength), no_step_rate=float(ns.part_strength_rate),
                                     no_step_mean=float(ns.part_strength_mean),
                                     no_step_from_boutons=float(bouton_no_step),
                                     below_no_step_by=float(ns.view_strength - sc.view_strength), **su,
                                     n_weighted_units=float(w.sum())))
                    nulls.append(pd.DataFrame(dict(stage=stage, glomeruli=gname, calyx=cal, fly_set=s,
                                                   perm=np.arange(len(nl)), view_strength=nl.view_strength.to_numpy(),
                                                   part_rate=nl.part_strength_rate.to_numpy(),
                                                   part_mean=nl.part_strength_mean.to_numpy())))
                    if s == 'all':
                        # entries by odor of the count-weighted yardstick, by repeating each glomerulus's row
                        # (checks weighted_arrays against hal's population_views as well)
                        pops = []
                        for p in ctx.pops['pn_dendrite']:
                            pp = p
                            if ex is not None:
                                V, Rr = fly_matrices([p], gl, ctx.odors, ex, crit)[0]
                                pp = type(p)(pd.DataFrame(V, index=gl, columns=ctx.odors),
                                             pd.DataFrame(Rr, index=gl, columns=ctx.odors), name=p.name,
                                             named_units=True, graded=p.graded)
                            pops.append(repeated_population(pp, dict(zip(gl, w))))
                        e = E138.odor_entries(pops, ref, names=[p.name for p in pops])
                        ec = e[e.row == 'common']
                        if not np.isclose(float(ec.view_error.iloc[0]), float(sc.view_strength), rtol=0, atol=1e-10):
                            raise RuntimeError(f'{cal} {stage} {gname}: the repeated rows give '
                                               f'{ec.view_error.iloc[0]}, the weights {sc.view_strength}')
                        ents.append(_tag(ec, stage=stage, glomeruli=gname, calyx=cal))
        print(f'yardsticks: {stage} ({time.time() - t0:.0f} s)', flush=True)
    out = pd.DataFrame(rows)
    _write_csv(counts, 'counts.csv')
    _write_csv(extra.rename_axis('glomerulus').reset_index(), 'carried.csv')
    _write_csv(out, 'yardsticks.csv')
    _write_csv(pd.concat(ents, ignore_index=True), 'yardstick_entries.csv')
    _write_csv(pd.concat(nulls, ignore_index=True), 'yardstick_null.csv')
    _write_json('yardsticks_info.json', dict(counts_check=check, carried=cinfo, glomeruli=sets_info(sets)))
    provenance('yardsticks', ctx, counts_check=check, carried=cinfo, seconds=round(time.time() - t0, 1))
    if not check['agrees']:
        raise RuntimeError(f"the counts differ from E138's feeders.csv: {check['differences']}")
    print(f'yardsticks: done ({time.time() - t0:.0f} s)', flush=True)
    return dict(yardsticks=out, counts=counts)


def sets_info(sets: Mapping[str, tuple]) -> dict:
    return {k: dict(glomeruli=list(v[0]), carried=(list(v[1].index) if v[1] is not None else [])) for k, v in
            sets.items()}


# ---------------------------------------------------------------------------
# piece 4: what the model's bouton reading samples

def bouton_sampling(ctx: Context) -> Dict[str, pd.DataFrame]:
    """Per wiring of the family held: each glomerulus's boutons in the wiring, the default reading's units (APL on,
    after the silent drop), each calyx's counts, and whether the glomerulus has model input or is among the reading's
    28; per wiring the shares and the correlations with the calyces' counts; the bouton flies' numbers of regions."""
    from data_prep import stage_population_data as D
    t0 = time.time()
    g = next(x for x in E142.MODEL_GROUPS if x.name == 'family_held')
    runs, _ = group_runs(g)
    counts = _read_csv('counts.csv')
    cw = {c: counts[counts.calyx == c].set_index('glomerulus') for c in E142.CALYCES}
    m36 = set(ctx.lock['glomeruli']['model36'])
    sub = set(ctx.lock['glomeruli']['bouton_reading_subset'])
    rows, summ = [], []
    for w, run in runs.items():
        gl = pd.Series(np.asarray(run.meta['bouton_glomerulus']).astype(str))
        nb = gl.value_counts()
        lab, _ = bouton_region_glomeruli(run, ctx.odors)
        nr = pd.Series(lab).value_counts()
        allg = sorted(set(nb.index) | set().union(*[set(x.index) for x in cw.values()]))
        for x in allg:
            rows.append(dict(wiring=w, glomerulus=x, model_input=x in m36, in_reading=x in sub,
                             wiring_boutons=int(nb.get(x, 0)), read_units=int(nr.get(x, 0)),
                             **{f'{c}__n_boutons': int(cw[c].n_boutons.get(x, 0)) for c in E142.CALYCES}))
        t = pd.DataFrame([r for r in rows if r['wiring'] == w]).set_index('glomerulus')
        rec = dict(wiring=w, n_boutons=int(nb.sum()), n_glomeruli=int((nb > 0).sum()),
                   n_read_units=int(nr.sum()), n_read_glomeruli=int((nr > 0).sum()),
                   share_boutons_no_input=float(t.wiring_boutons[~t.model_input].sum() / t.wiring_boutons.sum()),
                   n_glomeruli_no_input=int(((~t.model_input) & (t.wiring_boutons > 0)).sum()),
                   share_boutons_input_not_read=float(t.wiring_boutons[t.model_input & ~t.in_reading].sum()
                                                      / t.wiring_boutons.sum()),
                   share_boutons_in_reading=float(t.wiring_boutons[t.in_reading].sum() / t.wiring_boutons.sum()),
                   share_read_units_no_input=float(t.read_units[~t.model_input].sum() / max(t.read_units.sum(), 1)),
                   share_read_boutons_kept=float(nr.sum() / max(t.wiring_boutons[t.in_reading].sum(), 1)))
        r28 = t[t.in_reading]
        for c in E142.CALYCES:
            col = f'{c}__n_boutons'
            rec[f'r_read_units_vs_{c}'] = float(np.corrcoef(r28.read_units, r28[col])[0, 1])
            rec[f'r_wiring_boutons_vs_{c}'] = float(np.corrcoef(r28.wiring_boutons, r28[col])[0, 1])
            rec[f'share_{c}_boutons_on_reading28'] = float(cw[c].n_boutons.reindex(list(sub)).fillna(0).sum()
                                                          / cw[c].n_boutons.sum())
            rec[f'share_{c}_boutons_no_input'] = float(cw[c].n_boutons[[x not in m36 for x in cw[c].index]].sum()
                                                       / cw[c].n_boutons.sum())
        rec['wiring_equals_hemibrain'] = bool((t.wiring_boutons == t['hemibrain-first__n_boutons']).all())
        summ.append(rec)
    bp = D.load_population('pn_bouton', E142.PANEL, 'wild_type')
    flies = [dict(wiring=f'fly:{f}', n_read_units=int(len(bp.tables[f]))) for f in bp.flies]
    out, s = pd.DataFrame(rows), pd.concat([pd.DataFrame(summ), pd.DataFrame(flies)], ignore_index=True)
    _write_csv(out, 'bouton_sampling.csv')
    _write_csv(s, 'bouton_sampling_summary.csv')
    provenance('sampling', ctx, seconds=round(time.time() - t0, 1))
    print(f'sampling: done ({time.time() - t0:.0f} s)', flush=True)
    return dict(sampling=out, summary=s)


# ---------------------------------------------------------------------------
# piece 5: the equal-weight reading

def flies_equal_weight(ctx: Context) -> Tuple[Dict[str, object], pd.DataFrame, pd.DataFrame]:
    """The bouton flies under the equal-weight reading: per fly, the mean over its held-out repeats of the strength
    parts of its group units (`fly_group_units`); the references of every fly set (all five flies' one-fly units);
    the groups table; the check of E128's groups on the trial means against E119's a4_groups.csv / a4_summary.csv."""
    from data_prep import stage_population_data as D
    from experiments import E119_stages_per_roi as E119
    from experiments import E128_bouton_groups as E128
    from kcq import population_fit as P
    bp = D.load_population('pn_bouton', E142.PANEL, 'wild_type', trials=True)
    if list(bp.flies) != list(ctx.lock['recordings']['pn_bouton']['flies']):
        raise RuntimeError("the bouton flies differ from the lock's")
    g12 = list(E119.ODORS)
    views, grows, chk = {}, [], []
    a4 = pd.read_csv(ctx.lock['inputs']['results']['e119_a4_groups']['path'])
    s4 = pd.read_csv(ctx.lock['inputs']['results']['e119_a4_summary']['path']).set_index('acq')
    for f in bp.flies:
        tr = bp.trials[f]
        reps = sorted(set(tr.columns.get_level_values('repeat')))
        per = []
        for h in reps:
            pop, lab = fly_group_units(tr, ctx.odors, g12, h)
            v = P.population_views(pop, ctx.odors, ctx.edges, PARTS)
            per.append(v)
            ids, n = np.unique(lab, return_counts=True)
            grows += [dict(fly=f, held_out=int(h), group=int(i), n_rois=int(c), remainder=bool(i == 0))
                      for i, c in zip(ids, n)]
        views[f] = P.mean_views(per)
        # E128's groups on the trial means (all repeats), against E119's A4
        M = bp.tables[f].reindex(columns=g12).to_numpy(dtype=float)
        lab = E128.fly_groups(M)
        ids, cnt = np.unique(lab[lab > 0], return_counts=True)
        big = {int(i): int(c) for i, c in zip(ids, cnt) if c >= E128.MIN_GROUP}
        ref = {int(r.group): int(r.n_rois) for r in a4[a4.acq == f].itertuples()}
        row = s4.loc[f] if f in s4.index else None
        chk.append(dict(fly=f, n_rois=int(M.shape[0]), n_responsive=int((lab > 0).sum()), n_groups=int(len(ids)),
                        n_big=len(big), big_equal_e119=bool(big == ref),
                        e119_n_responsive=int(row.n_responsive) if row is not None else None,
                        e119_n_groups=int(row.n_tuning_groups) if row is not None else None,
                        agrees=bool(row is not None and big == ref and int((lab > 0).sum()) == int(row.n_responsive)
                                    and len(ids) == int(row.n_tuning_groups))))
    ref_all = reference_from_views(views, ctx.odors, ctx.edges)
    refs = {'all': ref_all}
    for h in E142.HALVES:
        refs[h] = P.subset_reference(ref_all, ctx.names['pn_bouton'][h], keep_units=True)
    return refs, pd.DataFrame(grows), pd.DataFrame(chk)


def score_views_row(views: Mapping[str, pd.Series], ref) -> dict:
    """The strength view's error and its two parts of a set of parts against a reference."""
    from kcq import population_fit as P
    s = P.score_views(dict(views), ref)
    return dict(view_strength=float(s.views['strength']), part_rate=float(s.parts['strength/rate']),
                part_mean=float(s.parts['strength/mean']))


def equal_weight_errors(ctx: Context, refs: Mapping[str, object]) -> pd.DataFrame:
    """The flies' floor and one fly against the others under the equal-weight reading, every saved run's equal-weight
    bouton errors (each wiring, 'common' the wirings' mean views, 'mean' their mean) on every fly set, beside its
    default-reading errors (models_entries.csv)."""
    from kcq import population_fit as P
    rows = []
    fl = P.floor_summary(P.floor(refs['all']))
    loo = P.leave_one_out(None, ref=refs['all'])
    rows.append(dict(group='flies', fly_set='all', row='floor', equal=float(fl.view_strength)))
    rows.append(dict(group='flies', fly_set='all', row='one_fly_rms', equal=float(_rms(loo.view_strength))))
    for h in E142.HALVES:
        other = 'held' if h == 'fit' else 'fit'
        rows.append(dict(group='flies', fly_set=h, row=f'half_to_half_from_{other}',
                         equal=score_views_row(refs[other].mean, refs[h])['view_strength']))
    ev = _read_csv('models_equal_views.csv')
    me = _read_csv('models_entries.csv')
    me = me[me.stage == 'pn_bouton'].drop_duplicates(['group', 'fly_set', 'row'])
    for (g, s), x in ev.groupby(['group', 'fly_set'], sort=False):
        ref = refs[s]
        per = {}
        for w, y in x.groupby('wiring', sort=False):
            per[w] = {f'strength/{k}': y[y.part == k].set_index('odor').value.reindex(ctx.odors) for k in ('rate', 'mean')}
        each = {w: score_views_row(v, ref) for w, v in per.items()}
        common = score_views_row(P.mean_views(list(per.values())), ref)
        d = me[(me.group == g) & (me.fly_set == s)].set_index('row').view_error
        for w, sc in each.items():
            rows.append(dict(group=g, fly_set=s, row=w, equal=sc['view_strength'], equal_rate=sc['part_rate'],
                             equal_mean=sc['part_mean'], default=float(d.get(w, np.nan)),
                             n_units=int(x[x.wiring == w].n_units.iloc[0])))
        rows.append(dict(group=g, fly_set=s, row='mean', equal=float(np.mean([v['view_strength'] for v in each.values()])),
                         default=float(np.mean([d.get(w, np.nan) for w in each]))))
        rows.append(dict(group=g, fly_set=s, row='common', equal=common['view_strength'], equal_rate=common['part_rate'],
                         equal_mean=common['part_mean'], default=float(d.get('common', np.nan))))
    return pd.DataFrame(rows)


def megamat_rows(path: Path, name_col: str) -> Dict[str, Dict[str, Dict[str, float]]]:
    """E136's or E139's score_megamat.csv, half 'held': per arm (or variant), per wiring its stage errors."""
    df = pd.read_csv(path)
    df = df[(df.half == 'held') & df.wiring.isin(E142.WIRINGS)]
    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    for r in df.itertuples():
        out.setdefault(str(getattr(r, name_col)), {})[str(r.wiring)] = {s: float(getattr(r, s)) for s in
                                                                         E142.LOSS_STAGES}
    return out


def with_stage(vals: Mapping[str, Dict[str, float]], stages: Sequence[str]) -> Dict[str, Dict[str, float]]:
    """Each wiring's values with 'stage' the RMS over `stages`."""
    return {w: dict(v, stage=_rms([v[s] for s in stages])) for w, v in vals.items()}


def reclasses(ctx: Context, errs: pd.DataFrame) -> dict:
    """E136's and E139's classes of the saved runs' arms read again: (1) with each wiring's bouton strength error under
    the equal-weight reading (against the held-out half) in place of the default's; (2) with the megamat stage error
    over the three other stages. The unrefitted family's values on its eight wirings give the spreads; E136's records
    give what is kept. "The reading's" from (1)."""
    res = ctx.lock['inputs']['results']
    e136c = json.loads(Path(res['e136_classes']['path']).read_text())
    e139c = json.loads(Path(res['e139_classes']['path']).read_text())
    m136 = megamat_rows(Path(res['e136_score_megamat']['path']), 'arm')
    m139 = megamat_rows(Path(res['e139_score_megamat']['path']), 'variant')
    fam = m136[E142.FAMILY]
    if len(fam) != len(E142.WIRINGS):
        raise RuntimeError(f"E136's family has {len(fam)} wirings")
    eq = errs[(errs.fly_set == 'held') & errs.row.isin(E142.WIRINGS)]
    eqv = {(r.group, r.row): float(r.equal) for r in eq.itertuples()}
    group_of = {'family': 'family_held'}
    arms = {'E136': [a for a in ('A0', 'A1', 'A4') if a in e136c['classes']],
            'E139': [g.name for g in E142.MODEL_GROUPS if g.source.startswith('E139:') and g.name in e139c['classes']]}

    def equal_vals(name: str, vals: Mapping[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
        grp = group_of.get(name, name)
        return {w: dict(v, pn_bouton=eqv[(grp, w)]) for w, v in vals.items() if (grp, w) in eqv}

    four, three = list(E142.LOSS_STAGES), [s for s in E142.LOSS_STAGES if s != 'pn_bouton']
    readings = {
        'equal_weight': (four, lambda n, v: equal_vals(n, v), ()),
        'without_boutons': (three, lambda n, v: dict(v), ('pn_bouton',)),
        'default': (four, lambda n, v: dict(v), ()),
    }
    out = {}
    for rname, (stages, fn, drop) in readings.items():
        f = with_stage(fn(E142.FAMILY, fam), stages)
        spreads = dict(stage=_sd(v['stage'] for v in f.values()), **{s: _sd(v[s] for v in f.values()) for s in stages})
        rec = dict(stages=stages, spreads=spreads, family_wirings=sorted(f), classes={})
        for src, names, cl, mv in (('E136', arms['E136'], e136c, m136), ('E139', arms['E139'], e139c, m139)):
            for a in names:
                if a not in mv:
                    rec['classes'][f'{src}:{a}'] = dict(cls='reads missing')
                    continue
                av = with_stage(fn(a, {w: v for w, v in mv[a].items() if w in E142.T_WIRINGS}), stages)
                rec['classes'][f'{src}:{a}'] = reclass(cl['classes'][a], av, f, spreads, stages,
                                                      cl.get('control_fall'), drop)
        out[rname] = rec
    # the default reading's re-read must give the original classes (a check of the decision tree and the values)
    out['default_check'] = {k: dict(original=v['original_cls'], reread=v['cls'], same=v['original_cls'] == v['cls'])
                            for k, v in out['default']['classes'].items()}
    # the reading's: A1's change at the boutons from the family under the equal-weight reading
    a1 = equal_vals('A1', {w: v for w, v in m136['A1'].items() if w in E142.T_WIRINGS})
    fe = equal_vals(E142.FAMILY, fam)
    ch = [a1[w]['pn_bouton'] - fe[w]['pn_bouton'] for w in E142.T_WIRINGS]
    spread = _sd(v['pn_bouton'] for v in fe.values())
    dch = [m136['A1'][w]['pn_bouton'] - fam[w]['pn_bouton'] for w in E142.T_WIRINGS]
    out['readings_class'] = dict(readings_class(float(np.mean(ch)), spread), changes=ch,
                                 a1_equal=float(np.mean([a1[w]['pn_bouton'] for w in E142.T_WIRINGS])),
                                 family_equal=float(np.mean([fe[w]['pn_bouton'] for w in E142.T_WIRINGS])),
                                 default_change=float(np.mean(dch)),
                                 default_spread=float(_sd(v['pn_bouton'] for v in fam.values())))
    return out


def equal_weight(ctx: Context) -> dict:
    """The equal-weight reading of the flies and of every saved run, and the classes read again."""
    t0 = time.time()
    refs, groups, chk = flies_equal_weight(ctx)
    errs = equal_weight_errors(ctx, refs)
    rc = reclasses(ctx, errs)
    _write_csv(groups, 'equal_weight_groups.csv')
    _write_csv(chk, 'equal_weight_groups_check.csv')
    _write_csv(errs, 'equal_weight_errors.csv')
    _write_json('reclasses.json', rc)
    provenance('equal', ctx, groups_check_agrees=bool(chk.agrees.all()), seconds=round(time.time() - t0, 1))
    print(f'equal: done ({time.time() - t0:.0f} s); groups agree with E119: {bool(chk.agrees.all())}', flush=True)
    return dict(errors=errs, reclasses=rc, groups_check=chk)


# ---------------------------------------------------------------------------
# piece 6: the planner's probe numbers

def _probe_row(name: str, quoted, got, source: str, tol: float = E142.PROBE_TOL, digits: Optional[int] = 2,
               **extra) -> dict:
    """One probe number: `got` rounded to `digits` (None: not rounded) within `tol` of `quoted`."""
    ok = False
    if quoted is not None and got is not None and np.isfinite(float(got)):
        g = round(float(got), digits) if digits is not None else float(got)
        ok = abs(g - float(quoted)) <= tol + 1e-12
    return dict(name=name, quoted=quoted, reproduced=got, difference=(float(got) - float(quoted)) if (
        quoted is not None and got is not None) else np.nan, reproduces=bool(ok), source=source, **extra)


def probe(ctx: Context) -> pd.DataFrame:
    """The planner's probe numbers (E142.PROBE) against step 0's own reads (models_entries.csv, all flies) and against
    the saved tables the planner read (E140's models_entries.csv / flies_entries.csv, E138's stages.csv, E136's and
    E139's scores, E134's and E135's)."""
    pr = E142.PROBE
    res = ctx.lock['inputs']['results']
    rows = []
    me = _read_csv('models_entries.csv')
    a = me[(me.fly_set == 'all')]
    cols = pr['lead_columns']

    def mine(col: str, odor: Optional[str] = None) -> pd.Series:
        if col == 'flies_dendrites':
            x = a[(a.group == 'family_held') & (a.stage == 'pn_dendrite') & (a.row == 'common')]
            return x.set_index('odor').rate_flies
        if col == 'flies_boutons':
            x = a[(a.group == 'family_held') & (a.stage == 'pn_bouton') & (a.row == 'common')]
            return x.set_index('odor').rate_flies
        x = a[(a.group == col) & (a.stage == 'pn_bouton') & (a.row == 'common')]
        return x.set_index('odor').rate_model
    e140m = pd.read_csv(res['e140_models_entries']['path'])
    e140f = pd.read_csv(res['e140_flies_entries']['path'])

    def e140(col: str) -> pd.Series:
        if col.startswith('flies_'):
            st = 'pn_dendrite' if col == 'flies_dendrites' else 'pn_bouton'
            x = e140f[(e140f.stage == st) & (e140f.cohort == 'wild_type') & (e140f.recording == 'mean')
                      & (e140f.reading == 'lock')]
            return x.set_index('odor').rate
        x = e140m[(e140m.group == col) & (e140m.condition == 'on') & (e140m.stage == 'pn_bouton')
                  & (e140m.reading == 'curve_noise') & (e140m.wiring == 'common')]
        return x.set_index('odor').rate
    series = {c: mine(c) for c in cols}
    series140 = {c: e140(c) for c in cols}
    for odor, vals in pr['lead_rate'].items():
        for c, q in zip(cols, vals):
            rows.append(_probe_row(f'lead {c} {odor}', q, float(series[c].get(odor, np.nan)), 'step 0 (all flies)',
                                   e140=float(series140[c].get(odor, np.nan))))
    fb = series['flies_boutons'].reindex(ctx.odors)
    for c, q, qr in zip(cols, pr['lead_sd'], pr['lead_r']):
        x = series[c].reindex(ctx.odors)
        # the planner's SDs are pandas' (ddof 1); ddof 0 beside
        rows.append(_probe_row(f'lead {c} SD over 17 odors (ddof 1)', q, float(x.std(ddof=1)), 'step 0 (all flies)',
                               ddof0=float(x.std(ddof=0)), e140=float(series140[c].reindex(ctx.odors).std(ddof=1))))
        if qr is not None:
            rows.append(_probe_row(f'lead {c} r with the flies\' boutons', qr, float(np.corrcoef(x, fb)[0, 1]),
                                   'step 0 (all flies)',
                                   e140=float(np.corrcoef(series140[c].reindex(ctx.odors),
                                                          series140['flies_boutons'].reindex(ctx.odors))[0, 1])))
    # E138's stages.csv for the family held (row common, the rate part at the boutons)
    st138 = pd.read_csv(res['e138_stages']['path'])
    x = st138[(st138.group == 'family_held') & (st138.stage == 'pn_bouton') & (st138.row == 'common')
              & (st138.condition == 'on')].set_index('odor').rate_model
    d = float(np.nanmax(np.abs(x.reindex(ctx.odors) - series['family_held'].reindex(ctx.odors))))
    rows.append(_probe_row("E138's stages.csv family held boutons (max |difference| over odors)", 0.0, d,
                           'E138 stages.csv against step 0', tol=1e-9, digits=None))
    # E139's stage table: half to half, the family, A1 on T against the held-out half (E136's scores and yardsticks)
    c136 = json.loads(Path(res['e136_classes']['path']).read_text())
    m136 = pd.read_csv(res['e136_score_megamat']['path'])
    held = m136[(m136.half == 'held') & m136.wiring.isin(E142.T_WIRINGS)]
    hh = c136.get('beside', {}).get('yardsticks', {})
    for st, (q_hh, q_f, q_a) in pr['stage_table'].items():
        f = float(held[held.arm == E142.FAMILY][st].mean())
        a1 = float(held[held.arm == 'A1'][st].mean())
        rows.append(_probe_row(f'stage table {st}: the family', q_f, f, "E136 score_megamat.csv (T, held)", tol=0.0005,
                               digits=3))
        rows.append(_probe_row(f'stage table {st}: A1', q_a, a1, "E136 score_megamat.csv (T, held)", tol=0.0005,
                               digits=3))
        if st in hh:
            rows.append(_probe_row(f'stage table {st}: half to half', q_hh, float(hh[st]['half_to_half']),
                                   "E136 classes.json beside.yardsticks", tol=0.0005, digits=3))
        elif st == 'stage' and all(s in hh for s in E142.LOSS_STAGES):
            rows.append(_probe_row('stage table stage: half to half', q_hh,
                                   _rms([hh[s]['half_to_half'] for s in E142.LOSS_STAGES]),
                                   'RMS of the stages\' half to half', tol=0.0005, digits=3))
    # the room's arithmetic (on the stage table's numbers, as the planner did)
    t = pr['stage_table']
    a1s = {s: t[s][2] for s in E142.LOSS_STAGES}
    rows.append(_probe_row('room: A1 with the boutons at half to half', pr['room']['boutons_at_half'],
                           _rms([a1s['pn_dendrite'], t['pn_bouton'][0], a1s['kc_claw'], a1s['kc_soma']]),
                           "arithmetic on the stage table", digits=3, tol=0.0005))
    rows.append(_probe_row('room: and the claws there too', pr['room']['boutons_and_claws_at_half'],
                           _rms([a1s['pn_dendrite'], t['pn_bouton'][0], t['kc_claw'][0], a1s['kc_soma']]),
                           "arithmetic on the stage table", digits=3, tol=0.0005))
    # E134: the family's boutons against the no-step yardstick (strength view, all flies)
    ns = pd.read_csv(res['e134_flies_no_step']['path'])
    ns = ns[ns.step.str.startswith('pn_dendrite') & ns.downstream.str.startswith('pn_bouton')
            & (ns.name == 'view_strength')]
    mine_ns = _read_csv('no_step_entries.csv')
    mine_ns = mine_ns[(mine_ns.upstream == 'pn_dendrite') & (mine_ns.downstream == 'pn_bouton')
                      & (mine_ns.fly_set == 'all')].view_error.iloc[0]
    rows.append(_probe_row("E134 no-step dendrites -> boutons (strength)", pr['e134_no_step']['no_step'],
                           float(ns.no_step.iloc[0]) if len(ns) else np.nan, 'E134 flies_no_step.csv',
                           step0=float(mine_ns)))
    fam_all = a[(a.group == 'family_held') & (a.stage == 'pn_bouton') & a.row.isin(E142.WIRINGS)]
    fam_each = fam_all.drop_duplicates('row').view_error.mean()
    fam_common = a[(a.group == 'family_held') & (a.stage == 'pn_bouton') & (a.row == 'common')].view_error.iloc[0]
    s134 = pd.read_csv(res['e134_score_scores']['path'])
    # the entry's "1.06" is E134's family as fitted (its write-up's stage table, strength, mean over wirings); the family
    # held (step 0's group) beside
    fit134 = s134[(s134.group == 'megamat_fitted') & (s134.stage == 'pn_bouton') & (s134.reading == 'curve_noise')
                  & s134.wiring.isin(E142.WIRINGS)]
    rows.append(_probe_row("E134 family boutons (strength, as fitted, mean over wirings, all flies)",
                           pr['e134_no_step']['family_boutons'], float(fit134.view_strength.mean()),
                           'E134 score_scores.csv (megamat_fitted)', held_step0=float(fam_each),
                           held_common_step0=float(fam_common)))
    # E135: the imposed dendrites' boutons against the family's (strength, all flies, eight wirings)
    s135 = pd.read_csv(res['e135_score_scores']['path'])
    x = s135[(s135.input_panel == 'megamat') & (s135.stage == 'pn_bouton') & (s135.reading == 'curve_noise')]
    imp = x[(x.arm == 'imposed') & x.wiring.isin(E142.WIRINGS)].view_strength.mean()
    imp_mine = a[(a.group == 'e135_imposed') & (a.stage == 'pn_bouton') & a.row.isin(E142.WIRINGS)].drop_duplicates(
        'row').view_error.mean()
    rows.append(_probe_row("E135 imposed boutons (strength, mean over eight wirings)", pr['e135_imposed']['imposed_boutons'],
                           float(imp), 'E135 score_scores.csv', step0=float(imp_mine)))
    y = s134[(s134.group == 'megamat_held') & (s134.stage == 'pn_bouton') & (s134.reading == 'curve_noise')
             & s134.wiring.isin(E142.WIRINGS)]
    rows.append(_probe_row("E135's family boutons (strength, mean over eight wirings; E134's megamat_held)",
                           pr['e135_imposed']['family_boutons'], float(y.view_strength.mean()),
                           'E134 score_scores.csv', step0=float(fam_each)))
    # step 0's default reads against E136's scores (the same reading: family on 8, A1 / A0 / A4 on T, held-out half)
    me_h = me[(me.fly_set == 'held') & (me.stage == 'pn_bouton')].drop_duplicates(['group', 'row'])
    m136h = m136[m136.half == 'held']
    for grp, arm in (('family_held', E142.FAMILY), ('A0', 'A0'), ('A1', 'A1'), ('A4', 'A4')):
        x = me_h[me_h.group == grp].set_index('row').view_error
        y = m136h[m136h.arm == arm].set_index('wiring').pn_bouton
        w = [k for k in x.index if k in y.index]
        d = float(np.max(np.abs(x[w] - y[w]))) if w else np.nan
        rows.append(_probe_row(f"step 0's {grp} bouton errors against E136's score (max |difference|, held)", 0.0, d,
                               'E136 score_megamat.csv', tol=1e-9, digits=None, n_wirings=len(w)))
    out = pd.DataFrame(rows)
    _write_csv(out, 'probe.csv')
    provenance('probe', ctx, n=len(out), n_reproduced=int(out.reproduces.sum()))
    print(f'probe: {int(out.reproduces.sum())} of {len(out)} reproduced', flush=True)
    return out


# ---------------------------------------------------------------------------
# piece 7: the classes and tables.md

def classes(ctx: Context) -> dict:
    """Step 0's classes ('counts_make_step' per calyx and connectome, with the beside readings; 'readings') and the
    outcomes they point to as far as step 0 can say; classes_step0.json and tables.md."""
    ys = _read_csv('yardsticks.csv')
    rc = json.loads((E142.out_dir() / 'reclasses.json').read_text())
    cm = {c: counts_class(ys, c) for c in E142.CALYCES}
    conn = {k: dict(calyces=list(v), holds=bool(all(cm[c]['holds'] for c in v))) for k, v in E142.CONNECTOMES.items()}
    beside = {f'{st}/{gl}/{c}': counts_class(ys, c, st, gl) for st in E142.COUNT_OF for gl in ('dendrites',
                                                                                              'with_carried')
              for c in E142.CALYCES if not (st == 'pn_bouton' and gl == 'dendrites')}
    any_counts = any(v['holds'] for v in conn.values())
    rd = rc['readings_class']
    out = dict(counts_make_step=dict(rule=E142.CLASSES['counts_make_step']['rule'], calyces=cm, connectomes=conn,
                                     holds=any_counts, beside=beside),
               readings=dict(rule=E142.CLASSES['readings']['rule'], **rd),
               reclasses={k: {a: dict(cls=v['cls'], original=v.get('original_cls')) for a, v in rc[k]['classes'].items()}
                          for k in ('equal_weight', 'without_boutons', 'default')},
               default_check=rc['default_check'],
               outcomes=dict(step0=(['a (if the reading does not sample the counts)'] if any_counts else
                                    ['not (a): the counts do not make the step']) + (['f'] if rd['holds'] else []),
                             step1='(b) to (e) wait for step 1'),
               lock_sha256=ctx.lock['sha256'], time=E142._now(), host=E142._host())
    _write_json('classes_step0.json', out)
    write_tables(ctx, out)
    provenance('classes', ctx)
    print(json.dumps(dict(counts={k: v['holds'] for k, v in conn.items()}, readings=rd['holds'])), flush=True)
    return out


def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '–'
    except (TypeError, ValueError):
        return str(x)


def write_tables(ctx: Context, cls: Mapping) -> None:
    """tables.md: the flies' step by odor, the no-step yardstick by odor, the saved runs' bouton and claw errors and
    concentration, the model's step, the yardsticks, the sampling, the equal-weight reading and the classes read
    again, the probe, and the classes."""
    d = E142.out_dir()
    L = ['# E142 step 0: the step from the PN dendrites to the PN boutons (no simulation)', '',
         "Entries are multiples of the mean over the 17 megamat odors; errors are the strength view's (rate and mean "
         "parts, RMS) in the downstream flies' one-fly units. Halves are E136's (boutons 2 / 3, dendrites 3 / 4, claws "
         '2 / 2).', '']
    st = _read_csv('flies_steps_summary.csv')
    for step in st.step.unique():
        L += [f"## The flies' step by odor, {step} (90 % interval on all flies; halves' changes)", '',
              '| odor | rate: change [interval] | fit | held | both halves | mean: change [interval] | fit | held | both '
              'halves |', '|---|---:|---:|---:|:---:|---:|---:|---:|:---:|']
        x = st[st.step == step]
        for o in ctx.odors:
            r = x[(x.odor == o) & (x.part == 'rate')].iloc[0]
            m = x[(x.odor == o) & (x.part == 'mean')].iloc[0]
            sym = {1: '+', -1: '-', 0: ''}
            L.append(f'| {o} | {_f(r.change_all, "+.2f")} [{_f(r.lo_all, "+.2f")}, {_f(r.hi_all, "+.2f")}] | '
                     f'{_f(r.change_fit, "+.2f")} | {_f(r.change_held, "+.2f")} | {sym[int(r.outside_both_halves)]} | '
                     f'{_f(m.change_all, "+.2f")} [{_f(m.lo_all, "+.2f")}, {_f(m.hi_all, "+.2f")}] | '
                     f'{_f(m.change_fit, "+.2f")} | {_f(m.change_held, "+.2f")} | {sym[int(m.outside_both_halves)]} |')
        L.append('')
    ns = _read_csv('no_step_entries.csv')
    L += ['## The no-step yardsticks (strength error; the upstream flies as the prediction)', '',
          '| upstream -> downstream | all | fit half | held half | largest shares (all flies) |', '|---|---:|---:|---:|---|']
    for (u, dn), x in ns.groupby(['upstream', 'downstream'], sort=False):
        e = {s: float(x[x.fly_set == s].view_error.iloc[0]) for s in E142.FLY_SETS}
        top = x[x.fly_set == 'all'].sort_values('view_share', ascending=False).head(3)
        L.append(f'| {u} -> {dn} | {_f(e["all"])} | {_f(e["fit"])} | {_f(e["held"])} | '
                 + ', '.join(f'{r.odor} {_f(r.view_share, ".2f")}' for r in top.itertuples()) + ' |')
    me = _read_csv('models_entries.csv')
    conc = _read_csv('models_concentrated.csv')
    L += ['', "## The saved runs: strength errors (wirings' mean views 'common'; held half on each run's wirings)", '',
          '| group | boutons all | boutons held | boutons: two largest (all) / concentrated | claws all | claws held | '
          'claws: two largest / concentrated |', '|---|---:|---:|---|---:|---:|---|']
    for g in [x.name for x in E142.MODEL_GROUPS if x.name in set(me.group)]:
        x = me[(me.group == g) & (me.row == 'common')].drop_duplicates(['fly_set', 'stage'])
        e = {(r.fly_set, r.stage): r.view_error for r in x.itertuples()}
        cb = conc[(conc.group == g) & (conc.stage == 'pn_bouton')]
        cc = conc[(conc.group == g) & (conc.stage == 'kc_claw')]
        L.append(f'| {g} | {_f(e.get(("all", "pn_bouton")))} | {_f(e.get(("held", "pn_bouton")))} | '
                 f'{cb.odors.iloc[0] if len(cb) else ""} ({_f(cb.share_all.iloc[0], ".2f") if len(cb) else ""}) / '
                 f'{cb.cls.iloc[0] if len(cb) else ""} | {_f(e.get(("all", "kc_claw")))} | '
                 f'{_f(e.get(("held", "kc_claw")))} | {cc.odors.iloc[0] if len(cc) else ""} / '
                 f'{cc.cls.iloc[0] if len(cc) else ""} |')
    ms = _read_csv('models_step.csv')
    L += ['', "## The model's step by odor against the flies' (r over the 17 odors, row 'common', all flies)", '',
          '| group | dendrites -> boutons, rate | mean | boutons -> claws, rate | mean |', '|---|---:|---:|---:|---:|']
    for g in [x.name for x in E142.MODEL_GROUPS if x.name in set(ms.group)]:
        x = ms[(ms.group == g) & (ms.row == 'common')].set_index(['step', 'part']).r
        L.append(f'| {g} | {_f(x.get(("pn_dendrite->pn_bouton", "rate")), ".2f")} | '
                 f'{_f(x.get(("pn_dendrite->pn_bouton", "mean")), ".2f")} | '
                 f'{_f(x.get(("pn_bouton->kc_claw", "rate")), ".2f")} | {_f(x.get(("pn_bouton->kc_claw", "mean")), ".2f")} |')
    ys = _read_csv('yardsticks.csv')
    L += ['', '## The count-weighted yardsticks (strength error; the null: 1,000 permutations of the counts)', '',
          '| stage | glomeruli | calyx | counts: all / fit / held | no-step: all / fit / held | null 5th pct (all) | '
          'percentile of the counts (all) |', '|---|---|---|---|---|---:|---:|']
    for (stg, gl, cal), x in ys.groupby(['stage', 'glomeruli', 'calyx'], sort=False):
        o = x.set_index('fly_set')
        L.append(f'| {stg} | {gl} | {cal} | ' + ' / '.join(_f(o.loc[s, 'counts']) for s in E142.FLY_SETS) + ' | '
                 + ' / '.join(_f(o.loc[s, 'no_step']) for s in E142.FLY_SETS) + f' | {_f(o.loc["all", "null_p05"])} | '
                 f'{_f(o.loc["all", "null_percentile_of_actual"], ".1f")} |')
    if (d / 'bouton_sampling_summary.csv').exists():
        bs = _read_csv('bouton_sampling_summary.csv')
        bw = bs[~bs.wiring.astype(str).str.startswith('fly:')]
        L += ['', "## What the model's bouton reading samples (the family held, APL on)", '',
              '| wiring | boutons | on glomeruli without input | read units (28 glomeruli) | r(read units, hemibrain '
              'counts) | r(read units, male CNS left) |', '|---|---:|---:|---:|---:|---:|']
        for r in bw.to_dict('records'):
            L.append(f"| {r['wiring']} | {int(r['n_boutons'])} | {_f(r['share_boutons_no_input'], '.2f')} | "
                     f"{int(r['n_read_units'])} | {_f(r['r_read_units_vs_hemibrain-first'], '.2f')} | "
                     f"{_f(r['r_read_units_vs_malecns-left'], '.2f')} |")
        fl = bs[bs.wiring.astype(str).str.startswith('fly:')]
        L.append('')
        L.append("The bouton flies' regions: " + ', '.join(f'{r.wiring[4:]} {int(r.n_read_units)}' for r in fl.itertuples()))
    if (d / 'equal_weight_errors.csv').exists():
        ew = _read_csv('equal_weight_errors.csv')
        L += ['', '## The equal-weight reading of the boutons (held half; T wirings for the arms)', '',
              '| group | equal weight, mean over wirings | default, mean over wirings |', '|---|---:|---:|']
        h = ew[(ew.fly_set == 'held') & ew.row.isin(E142.T_WIRINGS)]
        for g in [x.name for x in E142.MODEL_GROUPS if x.name in set(h.group)]:
            x = h[h.group == g]
            L.append(f'| {g} | {_f(x.equal.mean())} | {_f(x.default.mean())} |')
        fl = ew[ew.group == 'flies']
        L.append('')
        L.append("The flies under it: " + ', '.join(f'{r.fly_set} {r.row} {_f(r.equal)}' for r in fl.itertuples()))
        rc = json.loads((d / 'reclasses.json').read_text())
        L += ['', '## The classes read again', '', '| arm | original | equal-weight boutons | over the three other stages |',
              '|---|---|---|---|']
        for a in rc['default']['classes']:
            L.append(f"| {a} | {rc['default']['classes'][a].get('original_cls')} | "
                     f"{rc['equal_weight']['classes'][a]['cls']} | {rc['without_boutons']['classes'][a]['cls']} |")
    if (d / 'probe.csv').exists():
        pb = _read_csv('probe.csv')
        bad = pb[~pb.reproduces]
        L += ['', f"## The planner's probe numbers: {int(pb.reproduces.sum())} of {len(pb)} reproduced", '']
        if len(bad):
            L += ['| number | quoted | reproduced | source |', '|---|---:|---:|---|']
            for r in bad.itertuples():
                L.append(f'| {r.name} | {_f(r.quoted)} | {_f(r.reproduced)} | {r.source} |')
    L += ['', '## Classes', '', '```', json.dumps(cls, indent=1, default=E142._e140()._json_default), '```', '']
    (d / 'tables.md').write_text('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------

def run(pieces: Optional[Sequence[str]] = None) -> dict:
    """Step 0's pieces in order (all of PIECES by default), after the lock and E134's lock are required."""
    lock, e134 = E142.require_all()
    t0 = time.time()
    ctx = Context(lock, e134)
    fns = dict(flies=flies_steps, models=models_entries, yardsticks=yardsticks, sampling=bouton_sampling,
               equal=equal_weight, probe=probe, classes=classes)
    want = list(pieces or PIECES)
    bad = [p for p in want if p not in fns]
    if bad:
        raise ValueError(f'unknown pieces {bad}; {PIECES}')
    out = {}
    for p in PIECES:
        if p in want:
            out[p] = fns[p](ctx)
    print(f'step 0: {want} done ({time.time() - t0:.0f} s)', flush=True)
    return out
