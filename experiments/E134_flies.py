#!/usr/bin/env python3
"""E134 step 1: the flies alone. No model is read.

The entry is item 2 of experiment_priorities.md ("the models as populations at every recorded stage"), design step 4
("Step 1, the flies alone"). For every recorded population of `data_prep.stage_population_data` (read against the
committed lock):
- each view (strength by odor, sharing between odors on who responds and on how strongly, breadth among units, and for
  glomeruli the matrix itself), part by part: its one-fly spread (the root mean squared distance of one fly's part from
  the other flies' mean part, the part's unit), its floor (resampled sets of the same flies against the measured set,
  the mean and the 90 % range of the draws), and the flies' own range (each fly against the others, in one-fly units);
- the energy distance between each fly's units' tuning and the other flies', as the guard: as a multiple of the mean
  distance between two of the other flies, on the values as recorded and mapped rank for rank onto the population's
  pooled values;
- the recorded steps (`STEPS`) as changes of each part, with 90 % intervals from resampling the flies of both stages;
- the no-step yardstick of every downstream stage: the error, in the downstream flies' one-fly units, of the upstream
  flies' mean parts taken as the prediction, beside the downstream flies' own floor and range;
- Remy's megamat KC somata against Yang's control flies on the 12 odors the panels share, each in the other's one-fly
  units, at two criteria for the megamat flies (`MEGAMAT_CRITERIA`): the lock's (megamat's own level) and the APL16
  somata's 1.9684 (the lock's for every soma panel until the readings of 2026-11-24 11:11);
- in flies_table.md, the KC somata's criteria and response rates per set of recordings (the lock's soma_criteria).

The measures are hal's kcq/population_fit.py (`fly_reference`, `floor`, `leave_one_out`, `step_interval`,
`no_step_score`, `score_many`, `score_views`, `mean_views`, `fly_energy_matrix`, `rank_map`); this script only calls
them. The odors a population is read on are those every one of its flies was given (`P.common_odors`; validation2's
flies were each given part of the panel, so `V2_ODOR_SET` picks the largest set of odors enough of them share), and
the breadth bins are the data module's for that number of odors.

Outputs (<shared>/results/E134_population_baseline/): flies_views.csv, flies_loo.csv, flies_energy.csv,
flies_steps.csv, flies_steps_summary.csv, flies_no_step.csv, flies_megamat_vs_control.csv, flies_table.md,
provenance_flies.json.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_flies.py [--n-boot 1000] [--seed 0]
    ./venv/bin/python -u experiments/E134_flies.py --table      # flies_table.md from the CSVs
"""
from __future__ import annotations

import argparse
import dataclasses
import itertools
import json
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from data_prep import stage_population_data as D
from kcq import population_fit as P

EXPERIMENT = 'E134_population_baseline'
Key = Tuple[str, str, str]
# the recorded steps, upstream -> downstream, as population keys of data_prep.stage_population_data
STEPS: Tuple[Tuple[Key, Key], ...] = (
    (('orn', 'megamat', 'wild_type'), ('pn_dendrite', 'megamat', 'wild_type')),
    (('pn_dendrite', 'megamat', 'wild_type'), ('pn_bouton', 'megamat', 'wild_type')),
    (('pn_bouton', 'megamat', 'wild_type'), ('kc_claw', 'megamat', 'wild_type')),
    (('kc_claw', 'megamat', 'wild_type'), ('kc_soma', 'megamat', 'wild_type')),
    (('kc_claw', 'megamat', 'same_fly'), ('kc_soma', 'megamat', 'same_fly')),
    (('orn', 'validation2', 'wild_type'), ('kc_soma', 'validation2', 'wild_type')),
    (('orn', 'kmix', 'wild_type'), ('kc_soma', 'kmix', 'wild_type')),
    (('orn', 'cmix', 'wild_type'), ('kc_soma', 'cmix', 'wild_type')),
)
N_BOOT = 1000
SEED = 0
CI = 0.9
# validation2: the odor set a population is read on is the largest set of odors held by at least this many of its flies
V2_MIN_FLIES = 6
MEGAMAT12 = D.APL16_ODORS[:12]
# the megamat flies' criteria in the megamat-against-control reading: the lock's, and the APL16 somata's (the lock's
# before the readings of 2026-11-24 11:11)
MEGAMAT_CRITERIA: Tuple[str, ...] = ('lock', 'apl16')


def out_dir():
    d = D.RESULTS
    d.mkdir(parents=True, exist_ok=True)
    return d


def pname(key: Key) -> str:
    return D.pop_name(key)


# ---------------------------------------------------------------------------
# reading the populations

_REC: Dict[Key, D.Population] = {}


def recorded(key: Key) -> D.Population:
    """The recorded population, loaded against the committed lock (`D.require_lock`)."""
    if key not in _REC:
        D.require_lock()
        _REC[key] = D.load_population(*key)
    return _REC[key]


def flies_of(key: Key) -> List[P.Population]:
    """One hal `Population` per fly."""
    return P.fly_populations(recorded(key))


def v2_odor_set(flies: Sequence[P.Population], min_flies: int = V2_MIN_FLIES) -> Tuple[List[str], List[str]]:
    """The largest set of odors that at least `min_flies` of the flies were all given (ties: the set held by more flies,
    then the first in panel order), and those flies: the flies' own odor sets are candidates, and each candidate's
    holders are the flies given all of it."""
    sets = {f.name: set(f.odors) for f in flies}
    best = None
    for name, s in sets.items():
        cands = [s] + [s & t for t in sets.values()]
        for c in cands:
            holders = [n for n, t in sets.items() if c <= t]
            if len(holders) < min_flies:
                continue
            score = (len(c), len(holders))
            if best is None or score > best[0]:
                best = (score, c, holders)
    if best is None:
        raise ValueError(f'no odor set held by {min_flies} flies')
    order = [o for o in flies[0].values.columns if o in best[1]]
    return order, sorted(best[2])


def reading(key: Key) -> Tuple[List[P.Population], List[str]]:
    """The flies and the odors `key` is read on: all flies and the odors all of them were given, or for validation2's
    KC somata the `v2_odor_set` and its flies."""
    flies = flies_of(key)
    if key == ('kc_soma', 'validation2', 'wild_type'):
        odors, keep = v2_odor_set(flies)
        flies = [f for f in flies if f.name in keep]
        return flies, odors
    return flies, P.common_odors(flies)


def edges(n_odors: int):
    return D.breadth_bins(n_odors)


# ---------------------------------------------------------------------------
# each population alone

def _ci(x: pd.Series) -> Tuple[float, float]:
    q = ((1 - CI) / 2, 1 - (1 - CI) / 2)
    return float(np.nanquantile(x, q[0])), float(np.nanquantile(x, q[1]))


def population_rows(key: Key, n_boot: int = N_BOOT, seed: int = SEED) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """One population's views (one row per part, view and population error: the one-fly unit where it is a part's,
    the floor's mean and 90 % range, the flies' range against the others) and its leave-one-out rows (one per fly)."""
    flies, odors = reading(key)
    ref = P.fly_reference(flies, odors, edges(len(odors)))
    fl = P.floor(ref, n_boot=n_boot, seed=seed)
    loo = P.leave_one_out(flies, ref=ref)
    rows = []
    for col in fl.columns:
        if col == 'n_missing':
            continue
        level = 'part' if col.startswith('part_') else ('view' if col.startswith('view_') else 'population')
        name = col[5:] if level == 'part' else (col[5:] if level == 'view' else col)
        unit = np.nan
        if level == 'part':
            p = name.replace('_', '/', 1)
            unit = ref.one_fly.get(p, np.nan)
        lo, hi = _ci(fl[col])
        rows.append(dict(population=pname(key), stage=key[0], panel=key[1], cohort=key[2], level=level, name=name,
                         one_fly_unit=unit, floor=float(fl[col].mean()), floor_lo=lo, floor_hi=hi,
                         loo_mean=float(loo[col].mean()), loo_min=float(loo[col].min()),
                         loo_max=float(loo[col].max()), n_flies=len(flies), n_odors=len(odors),
                         odors=' '.join(odors), units=recorded(key).units, criterion=recorded(key).criterion))
    loo = loo.assign(population=pname(key)).rename_axis('fly').reset_index()
    return pd.DataFrame(rows), loo


def flies_views(n_boot: int = N_BOOT, seed: int = SEED) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per population with two or more flies, the views with spread, floor and range (flies_views.csv) and each fly
    against the others (flies_loo.csv)."""
    v, l = [], []
    for key in D.POPULATIONS:
        if len(recorded(key).flies) < 2:
            continue
        t0 = time.time()
        a, b = population_rows(key, n_boot, seed)
        v.append(a)
        l.append(b)
        print(f'views {pname(key)}: {time.time() - t0:.0f} s', flush=True)
    return pd.concat(v, ignore_index=True), pd.concat(l, ignore_index=True)


def flies_energy(max_units: int = P.ENERGY_MAX_UNITS, seed: int = SEED) -> pd.DataFrame:
    """Per population and fly, the energy guard against the other flies (as `P.energy_guard`, from one matrix of the
    flies' distances): the fly's mean distance to the others, the others' mean distance to each other and the ratio, on
    the values as recorded and mapped rank for rank onto the population's pooled values (flies_energy.csv)."""
    rows = []
    for key in D.POPULATIONS:
        if len(recorded(key).flies) < 3:
            continue
        flies, odors = reading(key)
        pooled = np.concatenate([f.values.reindex(columns=odors).to_numpy(dtype=float).ravel() for f in flies])
        for kind, refv in (('raw', None), ('ranked', pooled)):
            M = P.fly_energy_matrix(flies, odors, reference=refv, max_units=max_units, seed=seed).to_numpy()
            n = len(M)
            for i, f in enumerate(flies):
                oth = [j for j in range(n) if j != i]
                to = float(M[i, oth].mean())
                sub = M[np.ix_(oth, oth)]
                ff = float(sub[np.triu_indices(len(oth), k=1)].mean())
                rows.append(dict(population=pname(key), reading=kind, fly=f.name, to_flies=to, fly_to_fly=ff,
                                 ratio=to / ff if ff else np.nan, n_odors=len(odors)))
        print(f'energy {pname(key)}', flush=True)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the steps

def step_odors(up: Key, down: Key) -> List[str]:
    """The odors both stages' flies (as each is read) were all given, in the downstream panel's order."""
    fu, ou = reading(up)
    fd, od = reading(down)
    return [o for o in od if o in set(ou)]


def recorded_steps(n_boot: int = N_BOOT, seed: int = SEED) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per step of `STEPS` and part, the change of each entry from the upstream to the downstream flies' mean, with a
    90 % interval from resampling the flies of both stages (flies_steps.csv), and a summary per part: the root mean
    square of the change over the entries and the number of entries whose interval excludes 0
    (flies_steps_summary.csv). A stage of one fly (the same fly's claws and somata) gives the change with no interval."""
    rows, summ = [], []
    for up, down in STEPS:
        fu, _ = reading(up)
        fd, _ = reading(down)
        odors = step_odors(up, down)
        e = edges(len(odors))
        name = f'{pname(up)} -> {pname(down)}'
        if len(fu) >= 2 and len(fd) >= 2:
            si = P.step_interval(fu, fd, odors=odors, breadth_edges=e, n_boot=n_boot, ci=CI, seed=seed)
        else:
            vu = P.mean_views([P.population_views(f, odors, e) for f in fu])
            vd = P.mean_views([P.population_views(f, odors, e) for f in fd])
            si = {p: pd.DataFrame({'change': c, 'lo': np.nan, 'hi': np.nan}) for p, c in P.step_change(vu, vd).items()}
        for p, df in si.items():
            df = df.rename_axis('entry').reset_index()
            rows.append(df.assign(step=name, part=p, n_odors=len(odors)))
            excl = int(((df.lo > 0) | (df.hi < 0)).sum())
            summ.append(dict(step=name, part=p, n_entries=int(df.change.notna().sum()),
                             rms_change=float(np.sqrt(np.nanmean(df.change ** 2))),
                             n_interval_excludes_0=excl if df.lo.notna().any() else np.nan,
                             n_up_flies=len(fu), n_down_flies=len(fd), n_odors=len(odors)))
        print(f'step {name}', flush=True)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summ)


def no_step_yardstick(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per step of `STEPS` with two or more flies downstream, the no-step yardstick (`P.no_step_score`: the upstream
    flies' mean parts as the prediction of the downstream flies, in the downstream flies' one-fly units, on the odors
    both stages were given), beside the downstream flies' floor and their own range on the same odors
    (flies_no_step.csv)."""
    rows = []
    for up, down in STEPS:
        fu, _ = reading(up)
        fd, _ = reading(down)
        if len(fd) < 2:
            continue
        odors = step_odors(up, down)
        ref = P.fly_reference(fd, odors, edges(len(odors)))
        s = P.no_step_score(fu, ref).as_series()
        fl = P.floor(ref, n_boot=n_boot, seed=seed).mean()
        loo = P.leave_one_out(fd, ref=ref).mean()
        for col in s.index:
            if col == 'n_missing':
                continue
            rows.append(dict(step=f'{pname(up)} -> {pname(down)}', downstream=pname(down), name=col,
                             no_step=float(s[col]), down_floor=float(fl.get(col, np.nan)),
                             down_loo_mean=float(loo.get(col, np.nan)), n_odors=len(odors)))
        print(f'no step {pname(down)}', flush=True)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Remy's megamat somata against the control flies

def megamat_criterion(name: str, megamat: D.Population) -> float:
    """The megamat flies' criterion `name` of `MEGAMAT_CRITERIA`: 'lock' (the population's, megamat's own level) or
    'apl16' (the APL16 somata's, `D.CRITERIA['kc_soma']`)."""
    return float(megamat.criterion) if name == 'lock' else float(D.CRITERIA['kc_soma'])


def megamat_vs_control(n_boot: int = N_BOOT, seed: int = SEED,
                       criteria: Sequence[str] = MEGAMAT_CRITERIA) -> pd.DataFrame:
    """Remy's megamat KC somata against Yang's control flies on the 12 odors the panels share: each megamat fly and the
    megamat flies' mean parts in the control flies' one-fly units, and the reverse, each beside the reference's own
    floor and range on the same odors, with the megamat flies at each criterion of `criteria` (columns criterion,
    megamat_criterion, and the megamat flies' mean response rate over the 17 odors and the 12; the control flies keep
    the lock's 1.9684) (flies_megamat_vs_control.csv)."""
    odors = list(MEGAMAT12)
    e = edges(len(odors))
    mm = recorded(('kc_soma', 'megamat', 'wild_type'))
    control = flies_of(('kc_soma', 'apl16', 'control'))
    rows = []
    for cname in criteria:
        crit = megamat_criterion(cname, mm)
        mmc = dataclasses.replace(mm, criterion=crit)
        sides = {'megamat': P.fly_populations(mmc), 'control': control}
        rate17 = float(np.mean([D.fly_rate(mmc.tables[f], crit) for f in mmc.flies]))
        rate12 = float(np.mean([D.fly_rate(mmc.tables[f][odors], crit) for f in mmc.flies]))
        for ref_side, other in (('control', 'megamat'), ('megamat', 'control')):
            ref = P.fly_reference(sides[ref_side], odors, e)
            each = P.score_many(sides[other], ref)
            mv = P.mean_views([P.population_views(f, odors, e, list(ref.parts)) for f in sides[other]])
            common = P.score_views(mv, ref).as_series()
            fl = P.floor(ref, n_boot=n_boot, seed=seed).mean()
            loo = P.leave_one_out(sides[ref_side], ref=ref).mean()
            for col in common.index:
                if col == 'n_missing':
                    continue
                rows.append(dict(criterion=cname, megamat_criterion=crit, megamat_rate_17=rate17,
                                 megamat_rate_12=rate12, reference=ref_side, scored=other, name=col,
                                 each_mean=float(each[col].mean()), each_min=float(each[col].min()),
                                 each_max=float(each[col].max()), mean_views=float(common[col]),
                                 ref_floor=float(fl.get(col, np.nan)), ref_loo_mean=float(loo.get(col, np.nan)),
                                 n_scored=len(sides[other]), n_ref=len(sides[ref_side]), n_odors=len(odors)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the table

def _f(x, spec: str = '.2f') -> str:
    try:
        return '—' if x is None or not np.isfinite(float(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


STAGE_LABEL = {'orn': 'ORNs', 'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
               'kc_soma': 'KC somata'}
COHORT_LABEL = {'wild_type': '', 'same_fly': ' (the same fly)', 'control': ', control flies',
                'tnt': ', APL output blocked (TNT)'}


def _label(pop: str) -> str:
    s, p, c = pop.split('__')
    return f'{STAGE_LABEL[s]}, {p}{COHORT_LABEL[c]}'


def write_table() -> None:
    """flies_table.md from the CSVs: the populations' floors and ranges, the parts' floors, the energy guard, the steps
    and the no-step yardstick, Remy's megamat somata against the control flies."""
    od = out_dir()
    v = pd.read_csv(od / 'flies_views.csv')
    en = pd.read_csv(od / 'flies_energy.csv')
    ss = pd.read_csv(od / 'flies_steps_summary.csv')
    ns = pd.read_csv(od / 'flies_no_step.csv')
    mc = pd.read_csv(od / 'flies_megamat_vs_control.csv')
    L = ['# E134 step 1: the flies alone', '',
         f'Written by experiments/E134_flies.py on {D._host()}, {D._now()}. No model is read. Every number is in one-fly '
         'units (1: as far from the flies as one fly is from the others) unless the column says otherwise. The floor is '
         'the mean over resampled sets of the same flies; "one fly against the others" is the mean (and range) of each '
         'fly scored against the other flies.', '']
    L += ['## The population error per stage and panel', '',
          '| population | flies | odors | units | floor (90 % range) | one fly against the others (range) |',
          '|---|---:|---:|---|---:|---:|']
    for pop, g in v[v.level == 'population'].groupby('population', sort=False):
        r = g[g.name == 'error'].iloc[0]
        L.append(f"| {_label(pop)} | {r.n_flies} | {r.n_odors} | {r.units} | {_f(r.floor)} ({_f(r.floor_lo)} to "
                 f"{_f(r.floor_hi)}) | {_f(r.loo_mean)} ({_f(r.loo_min)} to {_f(r.loo_max)}) |")
    wm = v[(v.level == 'population') & (v.name == 'error_with_matrix') & v.floor.notna()]
    if len(wm):
        L += ['', 'With the matrix itself (glomeruli): ' + '; '.join(
            f"{_label(r.population)} floor {_f(r.floor)}, one fly against the others {_f(r.loo_mean)}"
            for r in wm.itertuples()) + '.']
    L += ['', '## The views: floor / one fly against the others', '']
    views = [x for x in ('strength', 'sharing', 'breadth', 'matrix')]
    L += ['| population | ' + ' | '.join(views) + ' |', '|---|' + '---:|' * len(views)]
    for pop, g in v[v.level == 'view'].groupby('population', sort=False):
        cells = []
        for w in views:
            r = g[g.name == w]
            cells.append('—' if r.empty else f'{_f(r.iloc[0].floor)} / {_f(r.iloc[0].loo_mean)}')
        L.append(f'| {_label(pop)} | ' + ' | '.join(cells) + ' |')
    L += ['', '## The energy guard (each fly against the others)', '',
          'The fly\'s mean energy distance to the other flies over the others\' mean distance to each other (1: as far '
          'as two of the others are apart), on the values as recorded and mapped rank for rank onto the population\'s '
          'pooled values.', '', '| population | as recorded: mean (range) | rank-mapped: mean (range) |', '|---|---:|---:|']
    for pop, g in en.groupby('population', sort=False):
        cells = []
        for kind in ('raw', 'ranked'):
            x = g[g.reading == kind].ratio
            cells.append(f'{_f(x.mean())} ({_f(x.min())} to {_f(x.max())})')
        L.append(f'| {_label(pop)} | ' + ' | '.join(cells) + ' |')
    L += ['', '## The recorded steps', '',
          'Per step and part: the root mean square over the entries of the change from the upstream flies\' mean to the '
          'downstream flies\' (in the part\'s own units: relative strengths, cosine distances, fractions of units), and '
          'the number of entries whose 90 % interval from resampling the flies of both stages excludes no change.', '']
    parts = list(dict.fromkeys(ss.part))
    L += ['| step | odors | ' + ' | '.join(parts) + ' |', '|---|---:|' + '---:|' * len(parts)]
    for step, g in ss.groupby('step', sort=False):
        cells = []
        for p in parts:
            r = g[g.part == p]
            if r.empty:
                cells.append('—')
            else:
                r = r.iloc[0]
                k = '' if not np.isfinite(r.n_interval_excludes_0) else f' ({int(r.n_interval_excludes_0)} of {r.n_entries})'
                cells.append(f'{_f(r.rms_change)}{k}')
        up, down = step.split(' -> ')
        L.append(f'| {_label(up)} to {_label(down)} | {g.n_odors.iloc[0]} | ' + ' | '.join(cells) + ' |')
    L += ['', '## The no-step yardstick', '',
          'The upstream flies\' mean views as the prediction of the downstream flies, in the downstream flies\' one-fly '
          'units, beside the downstream flies\' floor and one fly against the others on the same odors. A model explains '
          'a step only so far as it is nearer the downstream flies than this.', '',
          '| downstream (from) | odors | no step | downstream floor | downstream one fly against the others |',
          '|---|---:|---:|---:|---:|']
    for step, g in ns.groupby('step', sort=False):
        r = g[g.name == 'error'].iloc[0]
        up, down = step.split(' -> ')
        L.append(f'| {_label(down)} (from {_label(up)}) | {r.n_odors} | {_f(r.no_step)} | {_f(r.down_floor)} | '
                 f'{_f(r.down_loo_mean)} |')
    L += ['', 'By view (no step / downstream one fly against the others):', '']
    for step, g in ns.groupby('step', sort=False):
        up, down = step.split(' -> ')
        vv = g[g.name.str.startswith('view_')]
        L.append(f'- {_label(down)} from {_label(up)}: ' + '; '.join(
            f'{r.name[5:]} {_f(r.no_step)} / {_f(r.down_loo_mean)}' for r in vv.itertuples()))
    if 'criterion' not in mc:                                   # a CSV from before the two criteria
        mc = mc.assign(criterion='lock', megamat_criterion=np.nan, megamat_rate_17=np.nan, megamat_rate_12=np.nan)
    L += ['', "## Remy's megamat KC somata against the control flies (12 shared odors)", '',
          "The megamat flies at each criterion (the lock's, megamat's own level; and the APL16 somata's 1.9684, the "
          "lock's before the readings of 2026-11-24 11:11); the control flies at 1.9684 throughout.", '',
          '| megamat criterion | megamat flies\' response rate, 17 / 12 odors | reference flies | scored flies | '
          'each scored fly: mean (range) | their mean views | reference floor | reference one fly against the others |',
          '|---:|---:|---|---|---:|---:|---:|---:|']
    lab = {'control': 'control flies (APL16)', 'megamat': "Remy's megamat flies"}
    for (cn, rs, sc), g in mc.groupby(['criterion', 'reference', 'scored'], sort=False):
        r = g[g.name == 'error'].iloc[0]
        L.append(f'| {_f(r.megamat_criterion, ".4f")} ({cn}) | {_f(r.megamat_rate_17, ".4f")} / '
                 f'{_f(r.megamat_rate_12, ".4f")} | {lab[rs]} ({r.n_ref}) | {lab[sc]} ({r.n_scored}) | '
                 f'{_f(r.each_mean)} ({_f(r.each_min)} to {_f(r.each_max)}) | {_f(r.mean_views)} | {_f(r.ref_floor)} | '
                 f'{_f(r.ref_loo_mean)} |')
    L += ['', 'By view (their mean views):', '']
    for (cn, rs, sc), g in mc.groupby(['criterion', 'reference', 'scored'], sort=False):
        vv = g[g.name.str.startswith('view_')]
        L.append(f'- {sc} in {rs} units, megamat at the {cn} criterion: ' +
                 '; '.join(f'{r.name[5:]} {_f(r.mean_views)}' for r in vv.itertuples()))
    L += soma_criteria_lines(D.read_lock())
    L += ['', 'Concentrations: on the APL16 panel the model input of paa, euc and -aPine is ORN data within 1 log10 unit '
          'of the flies\' concentrations; the 12 shared odors here are all at 10^-3 in both panels. Validation2\'s KC '
          'flies were each given part of the panel; they are read on the largest odor set at least '
          f'{V2_MIN_FLIES} of them share (flies_views.csv, column odors).']
    (od / 'flies_table.md').write_text('\n'.join(L) + '\n')


def soma_criteria_lines(lock: Optional[dict]) -> List[str]:
    """flies_table.md's section on the KC somata's criteria per set of recordings, from the lock's soma_criteria (none
    when the lock has none)."""
    f = (lock or {}).get('soma_criteria')
    if not f:
        return []
    L = ['', "## The KC somata's response criteria, per set of recordings", '',
         "From the lock's soma_criteria (the planner's reading 14). Validation2's groups and the same fly are on megamat's "
         "footing: their values times megamat's median noise bridge over the group's "
         f"({_f(f['noise_bridge_reference']['median'], '.4f')}, over the four megamat flies). The last two columns are "
         "the lock's criterion for every soma panel before the readings of 2026-11-24 11:11.", '',
         '| KC somata | flies | criterion | values times | response rate, mean (range) | before: criterion | '
         'before: response rate |', '|---|---:|---:|---:|---:|---:|---:|']
    for r in f['table']:
        z = ' z' if 'session z' in r['units'] else ''
        L.append(f"| {r['cohort']} | {r['n_flies']} | {_f(r['criterion'], '.4f')}{z} | "
                 f"{_f(r['factor'], '.3f') if r['factor'] is not None else '—'} | {_f(r['rate_mean'], '.4f')} "
                 f"({_f(r['rate_min'], '.3f')} to {_f(r['rate_max'], '.3f')}) | {_f(r['apl16_criterion'], '.4f')}{z} | "
                 f"{_f(r['rate_mean_at_apl16_criterion'], '.4f')} |")
    g = f['validation2']['groups']
    L += ['', 'Noise bridge (median trial-noise SD in baseline-noise units, before the footing): megamat '
          f"{_f(f['noise_bridge_reference']['median'], '.3f')}; validation2 tsub03 {_f(g['tsub03']['noise_bridge_median'], '.3f')} "
          f"({_f(g['tsub03']['noise_bridge_min'], '.2f')} to {_f(g['tsub03']['noise_bridge_max'], '.2f')}), other "
          f"{_f(g['other']['noise_bridge_median'], '.3f')} ({_f(g['other']['noise_bridge_min'], '.2f')} to "
          f"{_f(g['other']['noise_bridge_max'], '.2f')}); the same fly "
          f"{_f(list(f['same_fly']['noise_bridge_per_fly'].values())[0], '.3f')}."]
    return L


def provenance(n_boot: int, seed: int, t0: float) -> None:
    from kcq import harness as H
    D._write_json(out_dir() / 'provenance_flies.json', dict(
        script=H.script_path(__file__), host=D._host(), time=D._now(), seconds=round(time.time() - t0),
        commits=H.code_commits(), n_boot=n_boot, seed=seed, ci=CI, v2_min_flies=V2_MIN_FLIES,
        megamat_criteria=list(MEGAMAT_CRITERIA), lock=str(D.LOCK_COPY), lock_time=(D.read_lock() or {}).get('time'),
        lock_sha256=D.sha256_file(D.LOCK_COPY) if D.LOCK_COPY.exists() else None))


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--table', action='store_true', help='only rewrite flies_table.md from the CSVs')
    a = ap.parse_args(argv)
    od = out_dir()
    if not a.table:
        t0 = time.time()
        D.require_lock()
        v, l = flies_views(a.n_boot, a.seed)
        v.to_csv(od / 'flies_views.csv', index=False)
        l.to_csv(od / 'flies_loo.csv', index=False)
        flies_energy(seed=a.seed).to_csv(od / 'flies_energy.csv', index=False)
        st, ss = recorded_steps(a.n_boot, a.seed)
        st.to_csv(od / 'flies_steps.csv', index=False)
        ss.to_csv(od / 'flies_steps_summary.csv', index=False)
        no_step_yardstick(a.n_boot, a.seed).to_csv(od / 'flies_no_step.csv', index=False)
        megamat_vs_control(a.n_boot, a.seed).to_csv(od / 'flies_megamat_vs_control.csv', index=False)
        provenance(a.n_boot, a.seed, t0)
    write_table()
    print(f'wrote {od}', flush=True)


if __name__ == '__main__':
    main()
