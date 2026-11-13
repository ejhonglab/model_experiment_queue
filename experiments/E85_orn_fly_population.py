#!/usr/bin/env python3
"""E85: the ORN input as a population of flies: each ORN fly's own responses as the model's input on the eight wirings,
read as the KC data are made, and the searched inputs set against the flies' spread.

The entry is item E85 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-12 10:43). C2's family throughout (the breadth-2 APL-effect fit of E67), each of the eight wirings at its own fit
of E72's part B, the 0.5 s pulse; S: real, P1, P2, M1; T: P3, P4, M2, M3 (E83's division). A model fly is one wiring with
one ORN fly's input; a population is the mean over its model flies of the pair distances per condition (E72).

Part A (no simulations):
- A1: the nine megamat ORN flies' trial-level responses (olfdb dataset 1, http://hal:8080, cached under
  outputs/E85_orn_fly_population/), the 12 megamat odors at the model input's concentrations, the model's 36 glomeruli.
  Each fly is scaled as the lab's pipeline scales it ('to-avg-max' in pyolfsysm's
  `scale_dff_to_est_spike_deltas_using_hallem`: its largest trial-mean response over all its panels and ROIs brought
  to the flies' average largest, AVG_FLYMAX from pyolfsysm's dff2spiking_model_choices.csv), then one slope to spike
  deltas, no constant: the least-squares slope through the origin of the lab's input (`H.apl16_orn_deltas()`) on the
  mean over flies; its R² is checked against 0.999 and reported (the entry: if the check fails, the regression's slope
  is used and the table says so).
- A2: a model fly's input: the lab's input plus s times (the fly's trial mean minus the mean over the flies that have
  the glomerulus), for the megamat odors and the glomeruli measured in that fly; the lab's value elsewhere. Variants:
  'repro' (s = r_g, the glomerulus's reliability: one minus the variance of a trial mean over the variance across flies
  of the trial means, pooled over odors, clipped to [0, 1], times the common factor of the A3 check); 'measured'
  (s = 1); 'gainfree' (the same as 'repro' on the fly's responses with its gain on each glomerulus divided out:
  least squares through the origin over the 12 odors against the mean over flies, for glomeruli whose largest mean
  response is at least GAIN_MIN_FRAC of the largest entry's, clipped to GAIN_CLIP; reliabilities recomputed on the
  divided responses); 'permuted' ('repro''s deviations of each fly permuted over its measured glomeruli within each
  odor, one generator per fly index, then re-centred over the flies per entry so that the mean over the nine flies is
  still the lab's input). The mean over the nine flies of every variant's input is the lab's input exactly.
- A3: fly_spread.csv (per entry: the SD across flies of log2(fly / mean) and in spike deltas, the reliability);
  orn_distances_by_fly.csv (per variant and fly, every pair's ORN distance, one minus Pearson over the 36 glomeruli,
  against the fly-mean vectors' and, for the 66 megamat pairs, the fly's own trial-mean and across-trial estimates on
  its measured glomeruli). The check: on the six near pairs (E78's SIX) the mean ratio of the 'repro' variant's
  distances (mean over flies) to the across-trial ones is within NEAR_TOL of 1; if not, one common factor on r_g is
  bisected to make it 1 and named in part_a.json.
- A4: searched_inputs_vs_fly_spread.csv: E77's two searched inputs and E83 arm I's two at SEARCH_SIZE of their log
  size, per entry against the flies' spread (log2 and spike deltas), with Hallem's data where E38 has it; and
  searched_inputs_span.csv: per input and odor, the share of the change's squared size (in spike deltas) inside the
  span of the nine flies' deviations, beside N_PERM label permutations of the factors (E83's `permuted_factors`).

Part B (simulations; `--part B --stage ...`, resumable): each model fly read as E77 / E83 read an input (breadth
recomputed for the input, the somatic threshold bisected to the wiring's APL-off response count, then the APL scale,
KC->APL and bouton->APL tied, to its APL-on count; `E77_input_screen.hold_rates`), except B5 (the fit's parameters,
nothing re-held). Stages (STAGES): B1S / B1T the 'repro' flies on S / T; B2 'measured', B3 'gainfree', B4 'permuted'
on S; B5 the 'repro' flies on all eight at the fits' parameters. `--part identity` first: the lab's input through the
same reading reproduces every wiring's saved fit exactly, and E83's unchanged readings (0.0500 / 0.0486 / 0.0482).

Part C (`--part C`, no simulations; runs on whatever part B has saved): readings.csv, pairs.csv, populations.csv,
draw_terms.csv, yardstick.csv, fly_mixture.csv, verdict.json, comparison_table.md, figures. Part D (the firing-ceiling
family) runs only if B1 moves the target (not written).

Outputs: results/E85_orn_fly_population/; spike counts and responses under outputs/E85_orn_fly_population/<part>/
(tag <stage>_<variant>_f<i>_<wiring>, conditions on / acute_off).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E85_orn_fly_population.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E85_orn_fly_population.py --part identity
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E85_orn_fly_population.py --part B
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E

EXPERIMENT = 'E85_orn_fly_population'
OLFDB_DATASET = 1
VARIANTS = ('repro', 'measured', 'gainfree', 'permuted')
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E.S_WIRINGS, E.T_WIRINGS, E.ALL_WIRINGS
# gains of the 'gainfree' variant: only glomeruli whose largest mean response is at least this share of the largest
# entry's, clipped to this range (1 elsewhere)
GAIN_MIN_FRAC = 0.1
GAIN_CLIP = (0.25, 4.0)
# the A3 check on the six near pairs: mean ratio of the 'repro' distances to the across-trial ones within this of 1
NEAR_TOL = 0.10
C_BRACKET = (0.1, 10.0)
SEARCH_SIZE = 0.6
N_PERM = 100
SEARCHED = {'E77_halfA': ('E77_input_side_common_error', 'oracle_factors_A.csv'),
            'E77_halfB': ('E77_input_side_common_error', 'oracle_factors_B.csv'),
            'E83_armI_halfA': ('E83_cross_wiring_fit', 'input_factors_A.csv'),
            'E83_armI_halfB': ('E83_cross_wiring_fit', 'input_factors_B.csv')}
# stage -> (variant, wirings, reading)
STAGES = {'B1S': ('repro', S_WIRINGS, 'held'), 'B1T': ('repro', T_WIRINGS, 'held'),
          'B2': ('measured', S_WIRINGS, 'held'), 'B3': ('gainfree', S_WIRINGS, 'held'),
          'B4': ('permuted', S_WIRINGS, 'held'), 'B5': ('repro', ALL_WIRINGS, 'free')}
STAGE_ORDER = tuple(STAGES)
# E83's unchanged readings at the lab's input (pattern part under the readout against all recordings)
UNCHANGED = {'S': 0.0500, 'T': 0.0486, 'all8': 0.0482}


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _local() -> Path:
    d = H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# A1: the ORN flies

def olfdb_long(refresh: bool = False) -> Tuple[pd.DataFrame, dict]:
    """olfdb dataset 1 (diagnostics and megamat, 9 flies) in its long format and its metadata, cached under
    outputs/E85_orn_fly_population/ (fetched from http://hal:8080 the first time)."""
    from data_prep import orn16_distances_from_olfdb as O
    p, pm = _local() / f'olfdb_dataset{OLFDB_DATASET}.parquet', _local() / f'olfdb_dataset{OLFDB_DATASET}.json'
    if refresh or not (p.exists() and pm.exists()):
        df, meta = O.fetch(OLFDB_DATASET)
        df.to_parquet(p)
        pm.write_text(json.dumps({k: meta.get(k) for k in ('id', 'title', 'data_sha256', 'n_flies', 'created_at')},
                                 indent=1) + '\n')
    return pd.read_parquet(p), json.loads(pm.read_text())


def avg_flymax() -> float:
    """The flies' average largest dF/F of the lab's dF/F -> spike delta model (pyolfsysm's
    dff2spiking_model_choices.csv)."""
    from pyolfsysm.mb_model import data_root
    c = pd.read_csv(Path(data_root) / 'internal' / 'dff2spiking_model_choices.csv', header=None, index_col=0)[1]
    if str(c['scaling_method']) != 'to-avg-max':
        raise RuntimeError(f"the lab's scaling is {c['scaling_method']}, not to-avg-max")
    return float(c['avg_flymax'])


_FLIES: Dict[str, object] = {}


def flies() -> dict:
    """The ORN flies, cached: 'labels' (f0 ... f8), 'keys' ((date, fly_num) per fly), 'odors' (MEGAMAT12), 'gloms' (the
    lab input's 36), 'x' (spike-delta trials, flies x odors x glomeruli x trials, NaN where not measured), 'slope',
    'r2', 'r2_origin', 'fly_max' (each fly's largest |trial-mean dF/F|), 'avg_flymax', 'meta' (olfdb's)."""
    if _FLIES:
        return _FLIES
    df, meta = olfdb_long()
    orn = H.abbrev_columns(H.apl16_orn_deltas())
    gloms, odors = list(orn.index), list(H.MEGAMAT12)
    keys = sorted({(str(d), int(f)) for d, f in zip(df.date, df.fly_num)})
    df = df.assign(fly=[keys.index((str(d), int(f))) for d, f in zip(df.date, df.fly_num)])
    # each fly's largest |trial mean| over all its panels and ROIs ('to-avg-max', applied to trial means)
    tm = df.groupby(['fly', 'panel', 'odor1', 'roi']).value.mean()
    fly_max = tm.abs().groupby(level='fly').max().reindex(range(len(keys)))
    want = {f'{o} @ {H.ORN_INPUT_CONC[o]:g}': o for o in odors}
    d = df[(df.panel == 'megamat') & df.odor1.isin(want) & df.roi.isin(gloms)]
    if set(d.odor1) != set(want):
        raise RuntimeError(f'olfdb dataset {OLFDB_DATASET}: missing {sorted(set(want) - set(d.odor1))}')
    reps = sorted(d['repeat'].unique())
    x = np.full((len(keys), len(odors), len(gloms), len(reps)), np.nan)
    fi = d.fly.to_numpy()
    oi = np.array([odors.index(want[o]) for o in d.odor1])
    gi = np.array([gloms.index(g) for g in d.roi])
    ti = np.array([reps.index(r) for r in d['repeat']])
    x[fi, oi, gi, ti] = d.value.to_numpy(dtype=float)
    A = avg_flymax()
    x = x / fly_max.to_numpy()[:, None, None, None] * A
    M = np.nanmean(np.nanmean(x, axis=3), axis=0)              # odors x glomeruli
    lab = orn[odors].T.to_numpy(dtype=float)                    # odors x glomeruli
    ok = ~np.isnan(M)
    slope = float((M[ok] @ lab[ok]) / (M[ok] @ M[ok]))
    res = lab[ok] - slope * M[ok]
    r2 = float(1.0 - (res @ res) / ((lab[ok] - lab[ok].mean()) @ (lab[ok] - lab[ok].mean())))
    r2o = float(1.0 - (res @ res) / (lab[ok] @ lab[ok]))
    _FLIES.update(labels=[f'f{i}' for i in range(len(keys))], keys=keys, odors=odors, gloms=gloms, x=x * slope,
                  slope=slope, r2=r2, r2_origin=r2o, fly_max=fly_max.to_numpy(), avg_flymax=A, meta=meta,
                  repeats=[int(r) for r in reps])
    return _FLIES


def trial_stats(x: np.ndarray) -> dict:
    """Trial means and reliabilities of flies x odors x glomeruli x trials responses: 'm' (trial means), 'M' (mean
    over the flies that have the entry), 'D' (m - M, 0 where not measured), 'measured' (bool), 'r_g' (per glomerulus,
    pooled over odors, clipped to [0, 1]), 'r_entry' (per odor and glomerulus), 'v_within' / 'v_across'."""
    import warnings
    with warnings.catch_warnings():
        # entries a fly did not measure (all-NaN slices)
        warnings.simplefilter('ignore', RuntimeWarning)
        m = np.nanmean(x, axis=3)
        vw = np.nanvar(x, axis=3, ddof=1) / np.sum(~np.isnan(x), axis=3)   # variance of a trial mean
    with np.errstate(invalid='ignore', divide='ignore'):
        M = np.nanmean(m, axis=0)
        va = np.nanvar(m, axis=0, ddof=1)                         # across flies, per odor x glomerulus
        vw_m = np.nanmean(vw, axis=0)
        r_g = 1.0 - np.nansum(vw_m, axis=0) / np.nansum(va, axis=0)
        r_entry = 1.0 - vw_m / va
    r_g = np.clip(np.nan_to_num(r_g, nan=0.0), 0.0, 1.0)
    measured = ~np.isnan(m)
    D = np.where(measured, m - M[None], 0.0)
    return dict(m=m, M=M, D=D, measured=measured, r_g=r_g, r_entry=np.clip(r_entry, 0.0, 1.0), v_within=vw_m,
                v_across=va)


def fly_gains(x: np.ndarray) -> Tuple[np.ndarray, dict]:
    """Each fly's gain on each glomerulus (flies x glomeruli): least squares through the origin of its trial means over
    the 12 odors on the mean over flies, for glomeruli whose largest mean response is at least GAIN_MIN_FRAC of the
    largest entry's; clipped to GAIN_CLIP; 1 elsewhere (and where not measured or not positive)."""
    st = trial_stats(x)
    m, M = st['m'], st['M']
    big = np.nanmax(M, axis=0) >= GAIN_MIN_FRAC * np.nanmax(M)
    num = np.nansum(m * M[None], axis=1)
    den = np.nansum(np.where(np.isnan(m), np.nan, M[None]) ** 2, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        a = num / den
    use = big[None, :] & np.isfinite(a) & (a > 0)
    clipped = use & ((a < GAIN_CLIP[0]) | (a > GAIN_CLIP[1]))
    a = np.where(use, np.clip(a, *GAIN_CLIP), 1.0)
    return a, dict(n_used=int(use.sum()), n_clipped=int(clipped.sum()), n_gloms_big=int(big.sum()))


def deviations(variant: str, c: float = 1.0) -> np.ndarray:
    """flies x odors x glomeruli deviations of `variant` (VARIANTS) to add to the lab's input of the megamat odors; `c`
    the common factor on the reliabilities (A3's check; not used by 'measured')."""
    F = flies()
    if variant == 'measured':
        return trial_stats(F['x'])['D']
    if variant in ('repro', 'permuted'):
        st = trial_stats(F['x'])
        d = np.clip(c * st['r_g'], 0.0, 1.0)[None, None, :] * st['D']
        if variant == 'repro':
            return d
        return recenter(permute_deviations(d, st['measured']), st['measured'])
    if variant == 'gainfree':
        a, _ = fly_gains(F['x'])
        st = trial_stats(F['x'] / a[:, None, :, None])
        return np.clip(c * st['r_g'], 0.0, 1.0)[None, None, :] * st['D']
    raise ValueError(variant)


def permute_deviations(d: np.ndarray, measured: np.ndarray) -> np.ndarray:
    """Each fly's deviations permuted over its measured glomeruli within each odor (generator seed = the fly's index,
    odors in order)."""
    out = d.copy()
    for f in range(d.shape[0]):
        rng = np.random.default_rng(f)
        for o in range(d.shape[1]):
            idx = np.flatnonzero(measured[f, o])
            out[f, o, idx] = d[f, o, idx[rng.permutation(len(idx))]]
    return out


def recenter(d: np.ndarray, measured: np.ndarray) -> np.ndarray:
    """Deviations re-centred per odor and glomerulus over the flies that have the entry (0 elsewhere), so that their sum
    over flies is 0."""
    cnt = measured.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = np.where(cnt > 0, np.where(measured, d, 0.0).sum(axis=0) / cnt, 0.0)
    return np.where(measured, d - mu[None], 0.0)


def model_input(dev: np.ndarray) -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()` layout) with odors x glomeruli deviations `dev` added to the megamat
    odors."""
    F = flies()
    orn = H.apl16_orn_deltas().astype(float).copy()
    ab = [H.odor_abbrev(c) for c in orn.columns]
    if list(orn.index) != F['gloms']:
        raise RuntimeError('glomerulus order differs from the flies\'')
    for j, o in enumerate(F['odors']):
        col = orn.columns[ab.index(o)]
        orn[col] = orn[col].to_numpy() + dev[j]
    return orn


_C: Dict[str, float] = {}


def common_factor() -> float:
    """The common factor on the reliabilities set by part A's check (part_a.json), 1 before part A has run."""
    if 'c' not in _C:
        p = out_dir() / 'part_a.json'
        _C['c'] = float(json.loads(p.read_text())['near_check']['c']) if p.exists() else 1.0
    return _C['c']


def fly_input(variant: str, f: int, c: Optional[float] = None) -> pd.DataFrame:
    """Model fly `f`'s input of `variant` (the common factor of part A unless `c` is given)."""
    return model_input(deviations(variant, common_factor() if c is None else c)[f])


# ---------------------------------------------------------------------------
# A3: the flies' spread and the ORN distances within a fly

def spread_table() -> pd.DataFrame:
    """Per odor and glomerulus: the mean over flies (spike deltas), its share of the largest entry, the number of flies,
    the SD across flies of log2(fly / mean) (flies and mean positive, at least 3 flies), the SD across flies in spike
    deltas, the share of the across-fly variance that repeats across trials (the entry's reliability), r_g, and the
    lab's input."""
    F = flies()
    st = trial_stats(F['x'])
    m, M = st['m'], st['M']
    lab = H.abbrev_columns(H.apl16_orn_deltas())[F['odors']].T.to_numpy(dtype=float)
    rows = []
    top = np.nanmax(M)
    for j, o in enumerate(F['odors']):
        for k, g in enumerate(F['gloms']):
            v = m[:, j, k]
            v = v[~np.isnan(v)]
            pos = v[v > 0]
            sdl = float(np.std(np.log2(pos / M[j, k]), ddof=1)) if (M[j, k] > 0 and len(pos) >= 3) else np.nan
            rows.append(dict(odor=o, glomerulus=g, mean=float(M[j, k]), share_of_largest=float(M[j, k] / top),
                             n_flies=int(len(v)), n_positive=int(len(pos)), sd_log2_ratio=sdl,
                             sd_spike_deltas=float(np.std(v, ddof=1)) if len(v) >= 2 else np.nan,
                             reliability=float(st['r_entry'][j, k]) if np.isfinite(st['r_entry'][j, k]) else np.nan,
                             r_g=float(st['r_g'][k]), lab_input=float(lab[j, k])))
    return pd.DataFrame(rows)


def _pearson_dist(a: np.ndarray, b: np.ndarray) -> float:
    return float(1.0 - np.corrcoef(a, b)[0, 1])


def within_fly_distances(f: int) -> pd.DataFrame:
    """Fly `f`'s own ORN distances of the 66 megamat pairs on its measured glomeruli: from its trial means and across
    trials (each trial vector centred over the glomeruli; covariances between different trials only, so that trial
    noise drops out)."""
    F = flies()
    x = F['x'][f]                                      # odors x glomeruli x trials
    meas = ~np.isnan(np.nanmean(x, axis=2)).any(axis=0)
    x = x[:, meas, :]
    m = np.nanmean(x, axis=2)
    xc = x - np.nanmean(x, axis=1, keepdims=True)
    nt = x.shape[2]
    pairs_t = [(t, u) for t in range(nt) for u in range(nt) if t != u]

    def cross(a, b):
        return float(np.mean([np.nansum(xc[a, :, t] * xc[b, :, u]) for t, u in pairs_t]))

    rows = []
    odors = F['odors']
    for a, b in combinations(range(len(odors)), 2):
        caa, cbb, cab = cross(a, a), cross(b, b), cross(a, b)
        r = cab / np.sqrt(caa * cbb) if caa > 0 and cbb > 0 else np.nan
        rows.append(dict(fly=F['labels'][f], odor_a=odors[a], odor_b=odors[b], n_gloms=int(meas.sum()),
                         trial_means=_pearson_dist(m[a], m[b]), across_trials=float(1.0 - r) if np.isfinite(r) else np.nan))
    return pd.DataFrame(rows)


def input_distances(orn: pd.DataFrame) -> pd.Series:
    """1 - Pearson over the 36 glomeruli of every pair of an input (canonical order)."""
    return E._cs(H.orn_input_pairs(orn))


def distances_table(c: float) -> pd.DataFrame:
    """orn_distances_by_fly.csv: per variant and fly, every pair's input distance, with the lab input's (the fly-mean
    vectors') and, for the megamat pairs, the fly's own trial-mean and across-trial distances on its measured
    glomeruli."""
    F = flies()
    lab = input_distances(H.apl16_orn_deltas())
    own = pd.concat([within_fly_distances(f) for f in range(len(F['labels']))], ignore_index=True)
    own.index = pd.MultiIndex.from_tuples([(fl, *sorted((a, b))) for fl, a, b in zip(own.fly, own.odor_a, own.odor_b)])
    parts = []
    for v in VARIANTS:
        dev = deviations(v, c)
        for f, fl in enumerate(F['labels']):
            d = input_distances(model_input(dev[f]))
            t = pd.DataFrame(dict(variant=v, fly=fl, odor_a=d.index.get_level_values(0), odor_b=d.index.get_level_values(1),
                                  distance=d.to_numpy(), lab_distance=lab.reindex(d.index).to_numpy()))
            key = [(fl, a, b) for a, b in d.index]
            t['own_trial_means'] = own['trial_means'].reindex(key).to_numpy()
            t['own_across_trials'] = own['across_trials'].reindex(key).to_numpy()
            parts.append(t)
    return pd.concat(parts, ignore_index=True)


def near_ratio(dist: pd.DataFrame, variant: str = 'repro') -> Tuple[float, pd.DataFrame]:
    """The mean over the six near pairs (E78's SIX) of the variant's distance (mean over flies) over the across-trial
    distance (mean over flies), and the per-pair table."""
    from experiments.E78_claw_knee import SIX
    d = dist[dist.variant == variant]
    g = d.groupby(['odor_a', 'odor_b'])[['distance', 'own_across_trials', 'own_trial_means', 'lab_distance']].mean()
    rows = g.reindex([tuple(sorted(p)) for p in SIX])
    rows = rows.assign(ratio=rows.distance / rows.own_across_trials)
    return float(rows.ratio.mean()), rows.reset_index()


# ---------------------------------------------------------------------------
# A4: the searched inputs against the flies' spread

def searched_factors(name: str, tables: Optional[Mapping[str, Tuple[str, str]]] = None) -> pd.DataFrame:
    """A searched input's factors as glomeruli x odors (E83's `read_factors`); `tables` maps names to (results folder,
    factor table) (SEARCHED by default)."""
    folder, f = (SEARCHED if tables is None else tables)[name]
    return E.read_factors(H.results_dir(folder) / f)


def hallem_log2() -> pd.Series:
    """log2(Hallem / lab) for the (odor, glomerulus) entries E38's shared table has, where both are positive."""
    p = H.results_dir('E38_hallem_input') / 'inputs' / 'hallem_shared_glomeruli.csv'
    h = pd.read_csv(p, index_col=0)
    lab = H.abbrev_columns(H.apl16_orn_deltas())
    out = {}
    for o in h.columns:
        for g in h.index:
            if o in lab.columns and g in lab.index and h.loc[g, o] > 0 and lab.loc[g, o] > 0:
                out[(o, g)] = float(np.log2(h.loc[g, o] / lab.loc[g, o]))
    return pd.Series(out, dtype=float)


def span_share(dev_o: np.ndarray, v: np.ndarray) -> float:
    """The share of |v|^2 inside the span of the rows of dev_o (flies x glomeruli)."""
    if not np.any(v):
        return np.nan
    q, s, _ = np.linalg.svd(dev_o.T, full_matrices=False)
    q = q[:, s > 1e-9 * max(s.max(), 1e-300)]
    p = q.T @ v
    return float(p @ p / (v @ v))


def searched_tables(tables: Optional[Mapping[str, Tuple[str, str]]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per entry, per input and odor) tables of part A4 (the module docstring), for the inputs of `tables` (name ->
    (results folder, factor table); SEARCHED by default; E86's part E passes its new searches)."""
    tables = SEARCHED if tables is None else tables
    F = flies()
    sp = spread_table().set_index(['odor', 'glomerulus'])
    st = trial_stats(F['x'])
    lab = H.abbrev_columns(H.apl16_orn_deltas())
    hal = hallem_log2()
    rows, srows = [], []
    for name in tables:
        fac = searched_factors(name, tables)
        searched = [o for o in fac.columns if not np.allclose(fac[o].to_numpy(), 1.0)]
        for o in searched:
            if o not in F['odors']:
                continue
            j = F['odors'].index(o)
            for g in F['gloms']:
                fv = float(fac.loc[g, o])
                l2 = SEARCH_SIZE * np.log2(fv)
                s = sp.loc[(o, g)]
                lv = float(lab.loc[g, o])
                dl = lv * (fv ** SEARCH_SIZE - 1.0)
                rows.append(dict(input=name, odor=o, glomerulus=g, factor=fv, log2_change=l2, lab_input=lv,
                                 sd_log2_ratio=s.sd_log2_ratio, z_log2=l2 / s.sd_log2_ratio if s.sd_log2_ratio > 0 else np.nan,
                                 change_spike_deltas=dl, sd_spike_deltas=s.sd_spike_deltas,
                                 z_spike_deltas=dl / s.sd_spike_deltas if s.sd_spike_deltas > 0 else np.nan,
                                 share_of_largest=s.share_of_largest, hallem_log2_vs_lab=hal.get((o, g), np.nan)))
            # the span of the nine flies' deviations (spike deltas), against N_PERM label permutations of the factors
            dev_o = st['D'][:, j, :]
            v = lab[o].reindex(F['gloms']).to_numpy(dtype=float) * (fac[o].reindex(F['gloms']).to_numpy() ** SEARCH_SIZE - 1.0)
            share = span_share(dev_o, v)
            perm = []
            for seed in range(1, N_PERM + 1):
                pf = E.permuted_factors(fac, seed)
                vp = lab[o].reindex(F['gloms']).to_numpy(dtype=float) * (pf[o].reindex(F['gloms']).to_numpy() ** SEARCH_SIZE - 1.0)
                perm.append(span_share(dev_o, vp))
            perm = np.asarray(perm)
            srows.append(dict(input=name, odor=o, share_in_span=share, perm_mean=float(np.nanmean(perm)),
                              perm_p95=float(np.nanpercentile(perm, 95)), perm_frac_ge=float(np.mean(perm >= share)),
                              sq_size=float(v @ v), rank=int(np.linalg.matrix_rank(dev_o))))
    return pd.DataFrame(rows), pd.DataFrame(srows)


def searched_summary(entries: pd.DataFrame, span: pd.DataFrame) -> pd.DataFrame:
    """Per searched input: shares of entries within one and two SDs of the flies' spread (log2 and spike deltas), and
    the squared-size-weighted share in the flies' span against the permutations'."""
    rows = []
    for name, e in entries.groupby('input', sort=False):
        zl = e.z_log2.abs().dropna()
        zs = e.z_spike_deltas.abs().dropna()
        sp = span[span.input == name]
        w = sp.sq_size / sp.sq_size.sum()
        rows.append(dict(input=name, n_entries=len(e), n_log2=len(zl), within1_log2=float((zl <= 1).mean()),
                         within2_log2=float((zl <= 2).mean()), rms_log2_change=float(np.sqrt((e.log2_change ** 2).mean())),
                         n_spike_deltas=len(zs), within1_spike_deltas=float((zs <= 1).mean()),
                         within2_spike_deltas=float((zs <= 2).mean()), span_share=float((w * sp.share_in_span).sum()),
                         span_share_perm_mean=float((w * sp.perm_mean).sum()),
                         n_odors_above_perm_p95=int((sp.share_in_span > sp.perm_p95).sum()), n_odors=len(sp)))
    return pd.DataFrame(rows)


def _write_input(name: str, orn: pd.DataFrame) -> Path:
    d = out_dir() / 'input_flies'
    d.mkdir(exist_ok=True)
    lay = H.abbrev_columns(orn)
    lay.index.name = 'glomerulus'
    lay.to_csv(d / f'{name}.csv')
    return d / f'{name}.csv'


def part_a() -> dict:
    """Part A: A1 to A4 (the module docstring). Writes flies.csv, fly_spread.csv, orn_distances_by_fly.csv,
    near_pairs_check.csv, searched_inputs_vs_fly_spread.csv, searched_inputs_span.csv, searched_inputs_summary.csv,
    input_flies/, part_a.json."""
    out = out_dir()
    F = flies()
    lab = H.apl16_orn_deltas()
    st = trial_stats(F['x'])
    pd.DataFrame(dict(fly=F['labels'], date=[k[0] for k in F['keys']], fly_num=[k[1] for k in F['keys']],
                      fly_max_dff=F['fly_max'], n_gloms=[int(st['measured'][f].any(axis=0).sum()) for f in range(len(F['labels']))])
                 ).to_csv(out / 'flies.csv', index=False)
    spread = spread_table()
    spread.to_csv(out / 'fly_spread.csv', index=False)
    # A3's check and the common factor on the reliabilities
    dist = distances_table(1.0)
    ratio1, near1 = near_ratio(dist)
    c, ratio = 1.0, ratio1
    if abs(ratio1 - 1.0) > NEAR_TOL:
        lo, hi = C_BRACKET
        f_lo = near_ratio(distances_table(lo))[0] - 1.0
        f_hi = near_ratio(distances_table(hi))[0] - 1.0
        if f_lo * f_hi > 0:
            c = lo if abs(f_lo) < abs(f_hi) else hi
        else:
            for _ in range(30):
                mid = float(np.sqrt(lo * hi))
                fm = near_ratio(distances_table(mid))[0] - 1.0
                if fm * f_lo > 0:
                    lo, f_lo = mid, fm
                else:
                    hi, f_hi = mid, fm
                if hi / lo < 1.001:
                    break
            c = float(np.sqrt(lo * hi))
        dist = distances_table(c)
        ratio = near_ratio(dist)[0]
    _C['c'] = c
    dist.to_csv(out / 'orn_distances_by_fly.csv', index=False)
    near = near_ratio(dist)[1].assign(c=c)
    near.to_csv(out / 'near_pairs_check.csv', index=False)
    # the model flies' inputs
    neg = {}
    for v in VARIANTS:
        dev = deviations(v, c)
        if not np.allclose(dev.sum(axis=0), 0.0, atol=1e-8):
            raise RuntimeError(f'{v}: the deviations do not sum to 0 over the flies')
        for f, fl in enumerate(F['labels']):
            orn_f = model_input(dev[f])
            _write_input(f'{v}_{fl}', orn_f)
            neg[f'{v}_{fl}'] = int((orn_f.to_numpy() < lab.to_numpy().min()).sum())
    # A4
    entries, span = searched_tables()
    entries.to_csv(out / 'searched_inputs_vs_fly_spread.csv', index=False)
    span.to_csv(out / 'searched_inputs_span.csv', index=False)
    summ = searched_summary(entries, span)
    summ.to_csv(out / 'searched_inputs_summary.csv', index=False)
    _, ginfo = fly_gains(F['x'])
    meg = spread[spread['mean'] > 0]
    res = dict(**H.code_commits(), olfdb=dict(dataset=OLFDB_DATASET, **{k: F['meta'].get(k) for k in ('title', 'data_sha256', 'n_flies')}),
               scaling=dict(method='to-avg-max on trial means over all panels and ROIs', avg_flymax=F['avg_flymax'],
                            slope=F['slope'], r2=F['r2'], r2_origin=F['r2_origin'], r2_check_0999=F['r2'] >= 0.999,
                            note=('the mean over flies is proportional to the lab input with R2 >= 0.999' if F['r2'] >= 0.999 else
                                  "R2 below 0.999: the slope is the regression's (least squares through the origin)")),
               reliability=dict(r_g_median=float(np.median(st['r_g'])), r_g_min=float(st['r_g'].min()),
                                r_g_max=float(st['r_g'].max())),
               gains=ginfo,
               near_check=dict(ratio_at_c1=ratio1, c=c, ratio=ratio, tol=NEAR_TOL, passed_at_c1=abs(ratio1 - 1.0) <= NEAR_TOL),
               spread=dict(median_sd_log2_all=float(spread.sd_log2_ratio.median()),
                           median_sd_log2_above10=float(meg[meg.share_of_largest >= 0.1].sd_log2_ratio.median()),
                           median_sd_log2_above30=float(meg[meg.share_of_largest >= 0.3].sd_log2_ratio.median())),
               inputs_below_lab_min=neg)
    (out / 'part_a.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(json.dumps({k: res[k] for k in ('scaling', 'reliability', 'gains', 'near_check', 'spread')}, indent=1, default=float), flush=True)
    print(summ.to_string(index=False), flush=True)
    return res


# ---------------------------------------------------------------------------
# B: the model flies

def build(label: str, orn: pd.DataFrame, input_tag: str):
    """`E83.build_simulator` with the breadth recomputed for `orn`, its APL-off breadth run cached under this
    experiment's outputs (by `input_tag`); the wiring and the normalized APL weights must be the lab input's."""
    from experiments import E10_joint_fit as J
    P = E.wiring_parts(label)
    _, wnorm, W, breadth = J.bouton_base(orn, _local(), None, input_tag=input_tag, **E.wiring_kws(label))
    if not W.equals(P['W']):
        raise RuntimeError(f'{label}: the wiring built with input {input_tag} differs from the lab input\'s')
    for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
        if not P['wnorm'][k].equals(wnorm[k]):
            raise RuntimeError(f'{label}: {k} built with input {input_tag} differs from the lab input\'s')
    return E.build_simulator(label, orn, breadth=breadth)


def held_reading(label: str, orn: pd.DataFrame, input_tag: str) -> dict:
    """E83's `input_reading` (both response rates held) with this experiment's cache."""
    from experiments import E77_input_screen as S
    P = E.saved_fit(label)
    res = S.hold_rates(lambda: build(label, orn, input_tag), target_off=P['target_off'], target_on=P['target_on'],
                       params=P['p'], label=f'{label} {input_tag}')
    res['input_tag'] = input_tag
    return res


def free_reading(label: str, orn: pd.DataFrame, input_tag: str) -> dict:
    """The wiring's fit's parameters with the input (breadth recomputed), nothing re-held: APL on, and APL off with the
    APL->KC and APL->bouton scales 0 (as hold_rates' APL-off runs)."""
    t0 = time.time()
    P = E.saved_fit(label)
    sim, info = build(label, orn, input_tag)
    p = dict(P['p'])
    on = sim(p)
    off = sim(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))
    del sim, info
    return dict(on=on, off=off, params=p, params_off=dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0), thr_factor=1.0,
                apl_factor=1.0, n_on=E._n(on[0]), n_off=E._n(off[0]), target_on=P['target_on'], target_off=P['target_off'],
                at_bound_thr=False, at_bound_apl=False, n_sims=2, history=[], breadth_sha256=None, input_tag=input_tag,
                seconds=time.time() - t0)


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = _local() / part
    d.mkdir(parents=True, exist_ok=True)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}): the local outputs, else the model-output archive; or None."""
    for d in (_local() / part, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E._load_pair(d, tag)
        if o is not None:
            return o
    return None


def stage_plan(stages: Sequence[str]) -> List[Tuple[str, str, int, str]]:
    """[(stage, variant, fly index, wiring)], S's wirings before T's, fly by fly within a wiring."""
    n = len(flies()['labels'])
    return [(s, STAGES[s][0], f, w) for s in stages for w in STAGES[s][1] for f in range(n)]


def run_tag(stage: str, variant: str, f: int, w: str) -> str:
    return f'{stage}_{variant}_f{f}_{w}'


def identity() -> dict:
    """The lab's input through this experiment's held reading on every wiring reproduces the saved fit exactly, and the
    unchanged populations read E83's values (UNCHANGED, four decimals). Writes identity.json; raises on a mismatch."""
    out = out_dir()
    lab = H.apl16_orn_deltas()
    res = dict(**H.code_commits(), wirings={}, unchanged={})
    for w in ALL_WIRINGS:
        r = held_reading(w, lab, 'E85_lab')
        ns, nr = E._same({'on': r['on'], 'off': r['off']}, E.saved_fit(w)['saved'])
        res['wirings'][w] = dict(thr_factor=r['thr_factor'], apl_factor=r['apl_factor'], count_diffs=ns, response_diffs=nr,
                                 n_sims=r['n_sims'], seconds=r['seconds'])
        print(w, res['wirings'][w], flush=True)
        if ns or nr:
            raise RuntimeError(f'{w}: the lab input through the reading differs from the saved fit in {ns} counts, {nr} responses')
        del r
        gc.collect()
    tg = E.targets()
    d = {w: E.pair_distances(E.saved_fit(w)['saved'], 'cal') for w in ALL_WIRINGS}
    for name, ws in (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS)):
        v = E._rms(E.pattern_error(E.population([d[w] for w in ws]), tg['z8']))
        res['unchanged'][name] = dict(value=v, expected=UNCHANGED[name], ok=round(v, 4) == UNCHANGED[name])
        if round(v, 4) != UNCHANGED[name]:
            raise RuntimeError(f'{name}: unchanged pattern part {v:.5f}, E83 {UNCHANGED[name]}')
    (out / 'identity.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(json.dumps(res['unchanged'], indent=1), flush=True)
    return res


def part_b(stages: Sequence[str] = STAGE_ORDER) -> pd.DataFrame:
    """Part B's stages in order; resumable (runs under outputs/E85_orn_fly_population/B/, rows in runs.csv, provenance
    per point in provenance_B.json). Needs part A (part_a.json, the common factor) and identity.json."""
    out = out_dir()
    if not (out / 'part_a.json').exists() or not (out / 'identity.json').exists():
        raise RuntimeError('run --part A and --part identity first')
    rp, pp = out / 'runs.csv', out / 'provenance_B.json'
    done = set(E._read_rows(rp).get('tag', pd.Series(dtype=str)))
    F = flies()
    tg = E.targets()
    c = common_factor()
    devs = {}
    for stage, variant, f, w in stage_plan(stages):
        tag = run_tag(stage, variant, f, w)
        if tag in done:
            continue
        if variant not in devs:
            devs[variant] = deviations(variant, c)
        orn_f = model_input(devs[variant][f])
        name = f'{variant}_{F["labels"][f]}'
        path = _write_input(name, orn_f)
        itag = E.input_tag_for(orn_f, f'E85_{name}')
        reading = STAGES[stage][2]
        res = held_reading(w, orn_f, itag) if reading == 'held' else free_reading(w, orn_f, itag)
        o = {'on': res['on'], 'off': res['off']}
        save_run('B', tag, o)
        row = E._run_row('B', f'{stage}_{variant}', w, tag, res, o, stage=stage, variant=variant, fly=F['labels'][f],
                         fly_index=f, reading=reading, c=c, input_tag=itag, orn_input=str(path.relative_to(H.shared_dir())))
        E._append_row(rp, row)
        E._provenance_add(pp, 'B', tag, H.model_provenance(orn_deltas=orn_f))
        d0, d1 = E.pair_distances(E.saved_fit(w)['saved'], 'cal'), E.pair_distances(o, 'cal')
        print(f"{tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, "
              f"{res['n_off']}/{res['target_off']} off); pattern under the readout against all recordings "
              f"{E._rms(E.pattern_error(d0, tg['z8'])):.4f} -> {E._rms(E.pattern_error(d1, tg['z8'])):.4f}; "
              f"flag '{row.get('realism_flag') or ''}'; {res['seconds']:.0f}s", flush=True)
        del res, o
        gc.collect()
    return E._read_rows(rp)


# ---------------------------------------------------------------------------
# C: the readings

READOUTS = ('cal', 'raw')
# thresholds of the entry, fixed before any run (on B1, the 'repro' flies held on all eight wirings)
MOVE_GAIN, NULL_GAIN, MIN_WIRINGS, RAW_SLACK, CONTROL_SHARE = 0.003, 0.0015, 6, 0.005, 0.5
MIXTURE_GAIN = 0.003
N_DRAWS = 100


def group_of(stage: str) -> str:
    """The reading group of a stage: B1S and B1T together are 'B1'."""
    return 'B1' if stage.startswith('B1') else stage


def _pairs_masks(index: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """Boolean masks over pairs (canonical order): 'nearest' (E77's four), 'six' (E78's SIX), 'megamat' (both odors
    megamat)."""
    from experiments.E78_claw_knee import NEAREST, SIX
    keys = [tuple(x) for x in index]
    a, b = index.get_level_values(0), index.get_level_values(1)
    return dict(nearest=np.array([k in {tuple(sorted(p)) for p in NEAREST} for k in keys]),
                six=np.array([k in {tuple(sorted(p)) for p in SIX} for k in keys]),
                megamat=np.asarray(a.isin(H.MEGAMAT12) & b.isin(H.MEGAMAT12)))


_DATA: Dict[str, Dict[str, pd.Series]] = {}


def data_distances() -> Dict[str, pd.Series]:
    """The KC data's mean pair distances, control flies ('on') and APL output blocked (TNT) ('off'), all recordings
    (E77's 'z8')."""
    if not _DATA:
        from experiments import E77_input_side_common_error as E77
        _DATA.update(E77.data_distances('z8'))
    return dict(_DATA)


def scores(d: Dict[str, pd.Series], tg: Dict[str, pd.Series]) -> dict:
    """E83's `effect_scores` plus the mean distances and ρ on odor-odor distances with APL on (against the control
    flies) and APL off (against APL output blocked (TNT)), the near pairs' distances and APL effect, and the pattern part
    over the 66 megamat pairs against all recordings."""
    row = E.effect_scores(d, tg)
    dd = data_distances()
    on, off = E._cs(d['on']), E._cs(d['off'])
    me = E._cs(H.apl_effect(on, off))
    masks = _pairs_masks(me.index)
    row.update(mean_on=float(on.mean()), mean_off=float(off.mean()),
               rho_dist_on=H.spearman(on, dd['on'].reindex(on.index)), rho_dist_off=H.spearman(off, dd['off'].reindex(off.index)))
    for k in ('nearest', 'six'):
        row[f'{k}_on'] = float(on[masks[k]].mean())
        row[f'{k}_off'] = float(off[masks[k]].mean())
        row[f'{k}_effect'] = float(me[masks[k]].mean())
    de = E._cs(tg['z8']).reindex(me.index)
    mm = masks['megamat']
    err = (me[mm] - me[mm].mean()) - (de[mm] - de[mm].mean())
    row['z8_pattern_megamat'] = E._rms(err.to_numpy())
    row['z8_loss_megamat'] = E._rms((me[mm] - de[mm]).to_numpy())
    return row


def data_scores() -> dict:
    """The KC data's own values of `scores`' descriptive columns (mean distances, near pairs)."""
    dd = data_distances()
    on, off = E._cs(dd['on']), E._cs(dd['off'])
    me = E._cs(H.apl_effect(on, off))
    masks = _pairs_masks(me.index)
    row = dict(mean_effect=float(me.mean()), mean_on=float(on.mean()), mean_off=float(off.mean()))
    for k in ('nearest', 'six'):
        row.update({f'{k}_on': float(on[masks[k]].mean()), f'{k}_off': float(off[masks[k]].mean()),
                    f'{k}_effect': float(me[masks[k]].mean())})
    return row


def load_runs() -> Tuple[pd.DataFrame, Dict[str, Dict[str, tuple]]]:
    """runs.csv and {tag: outputs} of every saved run."""
    rt = E._read_rows(out_dir() / 'runs.csv')
    runs = {}
    for r in rt.itertuples():
        o = load_run('B', r.tag)
        if o is None:
            print(f'{r.tag}: outputs missing', flush=True)
            continue
        runs[r.tag] = o
    return rt[rt.tag.isin(runs)].reset_index(drop=True), runs


def mixture_distances(fly_d: Sequence[Dict[str, pd.Series]], weights: Sequence[float]) -> Dict[str, pd.Series]:
    """The weighted mean over input flies of their populations' distances per condition (weights summing to 1); with
    equal weights, the population of every model fly in `fly_d`'s populations."""
    w = np.asarray(weights, dtype=float)
    return {c: E._cs(sum(wi * E._cs(d[c]) for wi, d in zip(w, fly_d))) for c in ('on', 'off')}


def fit_mixture(fly_d: Sequence[Dict[str, pd.Series]], target: pd.Series) -> Tuple[np.ndarray, float]:
    """Non-negative weights over the input flies, summing to 1, that minimize the pattern part of the mixture's
    per-pair APL effect against `target` (SLSQP from equal weights); (weights, the fitted pattern part)."""
    from scipy.optimize import minimize
    n = len(fly_d)
    on = np.stack([E._cs(d['on']).to_numpy() for d in fly_d])
    off = np.stack([E._cs(d['off']).to_numpy() for d in fly_d])
    idx = E._cs(fly_d[0]['on']).index
    t = E._cs(target).reindex(idx).to_numpy()
    t = t - t.mean()

    def f(w):
        e = w @ on - w @ off
        e = e - e.mean()
        return float(np.mean((e - t) ** 2))

    r = minimize(f, np.full(n, 1.0 / n), method='SLSQP', bounds=[(0.0, 1.0)] * n,
                 constraints=[dict(type='eq', fun=lambda w: w.sum() - 1.0)], options=dict(maxiter=500, ftol=1e-14))
    w = np.clip(r.x, 0.0, None)
    w = w / w.sum()
    return w, float(np.sqrt(f(w)))


def _dist(o: Dict[str, tuple], ro: str) -> Dict[str, pd.Series]:
    return E.pair_distances(o, ro)


def part_c() -> dict:
    """Part C on the runs saved so far: readings.csv, pairs.csv, populations.csv, draw_terms.csv, yardstick.csv,
    fly_mixture.csv, verdict.json, comparison_table.md, the figure, and the archive tables (model_pair_distances_long.csv,
    variants.csv)."""
    out = out_dir()
    tg = E.targets()
    rt, runs = load_runs()
    if rt.empty:
        raise RuntimeError('no runs yet')
    rt['group'] = [group_of(s) for s in rt.stage]
    flies_n = len(flies()['labels'])
    # distances of every run and of the wirings at the lab's input
    dist = {ro: {t: _dist(runs[t], ro) for t in rt.tag} for ro in READOUTS}
    lab = {ro: {w: _dist(E.saved_fit(w)['saved'], ro) for w in ALL_WIRINGS} for ro in READOUTS}
    rows, prs = [], []

    def add(group, unit, kind, n, ro, d, ref):
        sc = scores(d, tg)
        row = dict(group=group, unit=unit, kind=kind, n_model_flies=n, readout=ro, **sc)
        if ref is not None:
            sr = scores(ref, tg)
            for k in ('z8_pattern', 'z8_loss', 'mean_effect', 'z8_rho', 'z6_pattern', 'noise_free_pattern', 'noise8_pattern',
                      'halfA_pattern', 'halfB_pattern', 'mean_on', 'mean_off', 'rho_dist_on', 'rho_dist_off', 'nearest_effect',
                      'six_effect', 'z8_pattern_megamat'):
                row[f'lab_{k}'] = sr[k]
                row[f'gain_{k}'] = sr[k] - sc[k]
        rows.append(row)
        me = E._cs(H.apl_effect(d['on'], d['off']))
        prs.append(pd.DataFrame(dict(group=group, unit=unit, readout=ro, odor_a=me.index.get_level_values(0),
                                     odor_b=me.index.get_level_values(1), on=E._cs(d['on']).reindex(me.index).to_numpy(),
                                     off=E._cs(d['off']).reindex(me.index).to_numpy(), effect=me.to_numpy())))

    pops = {'S': S_WIRINGS, 'T': T_WIRINGS, 'all': ALL_WIRINGS}
    for ro in READOUTS:
        for w in ALL_WIRINGS:
            add('lab', w, 'wiring', 1, ro, lab[ro][w], None)
        for name, ws in pops.items():
            add('lab', name, 'population', len(ws), ro, E.population([lab[ro][w] for w in ws]), None)
        for g, gt in rt.groupby('group', sort=False):
            for r in gt.itertuples():
                add(g, r.tag, 'model_fly', 1, ro, dist[ro][r.tag], lab[ro][r.wiring])
            for w, wt in gt.groupby('wiring', sort=False):
                if len(wt) == flies_n:
                    add(g, w, 'wiring_nine_flies', flies_n, ro, E.population([dist[ro][t] for t in wt.tag]), lab[ro][w])
            for name, ws in pops.items():
                sub = gt[gt.wiring.isin(ws)]
                if len(sub) == flies_n * len(ws):
                    add(g, name, 'population', len(sub), ro, E.population([dist[ro][t] for t in sub.tag]),
                        E.population([lab[ro][w] for w in ws]))
    readings = pd.DataFrame(rows)
    readings.to_csv(out / 'readings.csv', index=False)
    pd.concat(prs, ignore_index=True).to_csv(out / 'pairs.csv', index=False)
    pop = readings[readings.kind.isin(['population', 'wiring_nine_flies'])]
    pop.to_csv(out / 'populations.csv', index=False)

    res = dict(**H.code_commits(), n_runs=int(len(rt)), groups={g: int(n) for g, n in rt.group.value_counts().items()})
    b1 = rt[rt.group == 'B1']
    full_b1 = len(b1) == flies_n * len(ALL_WIRINGS)
    # draws of eight model flies, eight wirings and eight input flies
    if full_b1:
        key = {(r.wiring, r.fly_index): r.tag for r in b1.itertuples()}
        rng = np.random.default_rng(0)
        drows = []
        for k in range(N_DRAWS):
            perm = rng.permutation(flies_n)[:len(ALL_WIRINGS)]
            tags = [key[(w, int(f))] for w, f in zip(ALL_WIRINGS, perm)]
            for ro in READOUTS:
                sc = scores(E.population([dist[ro][t] for t in tags]), tg)
                drows.append(dict(draw=k, readout=ro, flies=' '.join(f'f{int(f)}' for f in perm),
                                  **{c: sc[c] for c in ('z8_pattern', 'z8_loss', 'mean_effect', 'z8_rho', 'mean_on', 'mean_off',
                                                        'rho_dist_on', 'rho_dist_off', 'six_effect', 'nearest_effect')}))
        draws = pd.DataFrame(drows)
        draws.to_csv(out / 'draws.csv', index=False)
        # the draw terms of the per-pair APL effect (under the readout): wiring, input fly, interaction
        dt = []
        for ro in READOUTS:
            eff = np.stack([[E._cs(H.apl_effect(dist[ro][key[(w, f)]]['on'], dist[ro][key[(w, f)]]['off'])).to_numpy()
                             for f in range(flies_n)] for w in ALL_WIRINGS])           # wirings x flies x pairs
            for centred in (False, True):
                e = eff - eff.mean(axis=2, keepdims=True) if centred else eff
                gm = e.mean(axis=(0, 1))
                wm = e.mean(axis=1) - gm
                fm = e.mean(axis=0) - gm
                inter = e - gm - wm[:, None, :] - fm[None, :, :]
                single = readings[(readings.group == 'B1') & (readings.kind == 'model_fly') & (readings.readout == ro)]
                single = single.merge(b1[['tag', 'wiring', 'fly_index']], left_on='unit', right_on='tag')
                labw = readings[(readings.group == 'lab') & (readings.kind == 'wiring') & (readings.readout == ro)].set_index('unit')
                sd_in = float(np.sqrt(single.groupby('wiring').z8_pattern.var(ddof=1).mean()))
                dt.append(dict(readout=ro, effect=('pattern (each model fly mean removed)' if centred else 'per-pair APL effect'),
                               var_wiring=float((wm ** 2).mean()), var_fly=float((fm ** 2).mean()),
                               var_interaction=float((inter ** 2).mean()),
                               sd_single_pattern_within_wiring=sd_in,
                               sd_single_pattern_over_wirings_lab=float(labw.loc[list(ALL_WIRINGS), 'z8_pattern'].std(ddof=1)),
                               mean_single_minus_lab=float((single.z8_pattern - single.wiring.map(labw.z8_pattern)).mean())))
        pd.DataFrame(dt).to_csv(out / 'draw_terms.csv', index=False)
        # the yardstick: APL on from one model fly, APL off from another of a different wiring and input fly
        yrows = []
        de = E._cs(tg['z8'])
        for ro in READOUTS:
            items = [(w, f, dist[ro][key[(w, f)]]) for w in ALL_WIRINGS for f in range(flies_n)]
            vals = []
            for w1, f1, d1 in items:
                for w2, f2, d2 in items:
                    if w1 == w2 or f1 == f2:
                        continue
                    me = E._cs(H.apl_effect(d1['on'], d2['off']))
                    dd = de.reindex(me.index)
                    vals.append((E._rms((me - dd).to_numpy()), E._rms(((me - me.mean()) - (dd - dd.mean())).to_numpy())))
            lv = []
            for w1 in ALL_WIRINGS:
                for w2 in ALL_WIRINGS:
                    if w1 != w2:
                        me = E._cs(H.apl_effect(lab[ro][w1]['on'], lab[ro][w2]['off']))
                        dd = de.reindex(me.index)
                        lv.append(E._rms((me - dd).to_numpy()))
            v = np.asarray(vals)
            yrows.append(dict(readout=ro, n=len(v), loss_median=float(np.median(v[:, 0])), loss_mean=float(v[:, 0].mean()),
                              pattern_median=float(np.median(v[:, 1])), lab_two_wirings_loss_median=float(np.median(lv))))
        pd.DataFrame(yrows).to_csv(out / 'yardstick.csv', index=False)
        # the mixture of input flies: fitted on S against one half, scored on T against the other
        mrows = []
        for X, Y in (('A', 'B'), ('B', 'A')):
            fly_S = [E.population([dist['cal'][key[(w, f)]] for w in S_WIRINGS]) for f in range(flies_n)]
            fly_T = [E.population([dist['cal'][key[(w, f)]] for w in T_WIRINGS]) for f in range(flies_n)]
            w_fit, fitted = fit_mixture(fly_S, tg[f'half{X}'])
            eq = np.full(flies_n, 1.0 / flies_n)
            sT_fit = E._rms(E.pattern_error(mixture_distances(fly_T, w_fit), tg[f'half{Y}']))
            sT_eq = E._rms(E.pattern_error(mixture_distances(fly_T, eq), tg[f'half{Y}']))
            sT_lab = E._rms(E.pattern_error(E.population([lab['cal'][w] for w in T_WIRINGS]), tg[f'half{Y}']))
            sS_eq = E._rms(E.pattern_error(mixture_distances(fly_S, eq), tg[f'half{X}']))
            mrows.append(dict(fitted_on=X, scored_on=Y, S_fitted=fitted, S_equal=sS_eq, T_unseen_fitted=sT_fit,
                              T_unseen_equal=sT_eq, T_unseen_lab=sT_lab, gain_vs_equal=sT_eq - sT_fit, gain_vs_lab=sT_lab - sT_fit,
                              **{f'w_f{i}': float(x) for i, x in enumerate(w_fit)}))
        mix = pd.DataFrame(mrows)
        mix.to_csv(out / 'fly_mixture.csv', index=False)
        res['mixture'] = dict(passes=bool((mix.gain_vs_equal >= MIXTURE_GAIN).all()),
                              gain_vs_equal=[float(x) for x in mix.gain_vs_equal], gain_vs_lab=[float(x) for x in mix.gain_vs_lab])
    # the verdict on B1, the control and the matrices
    r = readings[readings.readout == 'cal'].set_index(['group', 'unit'])
    rr = readings[readings.readout == 'raw'].set_index(['group', 'unit'])
    if ('B1', 'all') in r.index:
        g_all = float(r.loc[('B1', 'all'), 'gain_z8_pattern'])
        wires = r.loc['B1'].reindex(list(ALL_WIRINGS))
        n_lower = int((wires.gain_z8_pattern > 0).sum())
        pr = population_realism_table(rt, runs)
        pb = pr[(pr.group == 'B1')]
        flags = int((pb.realism_flag.fillna('').astype(str).str.strip() != '').sum())
        run_flags = int((b1.realism_flag.fillna('').astype(str).str.strip() != '').sum())
        raw_change = float(rr.loc[('B1', 'all'), 'z8_loss'] - rr.loc[('lab', 'all'), 'z8_loss'])
        moves = bool(g_all >= MOVE_GAIN and n_lower >= MIN_WIRINGS and flags == 0 and raw_change <= RAW_SLACK)
        res['B1'] = dict(pattern=float(r.loc[('B1', 'all'), 'z8_pattern']), lab=float(r.loc[('lab', 'all'), 'z8_pattern']),
                         gain=g_all, n_wirings_lower=n_lower, n_flagged_populations=flags, n_flagged_model_flies=run_flags,
                         raw_loss_change=raw_change, moves=moves,
                         no_move=bool(abs(g_all) < NULL_GAIN), worse=bool(g_all <= -MOVE_GAIN),
                         rho_dist_on_change=float(r.loc[('B1', 'all'), 'rho_dist_on'] - r.loc[('lab', 'all'), 'rho_dist_on']),
                         mean_on_change=float(r.loc[('B1', 'all'), 'mean_on'] - r.loc[('lab', 'all'), 'mean_on']),
                         mean_off_change=float(r.loc[('B1', 'all'), 'mean_off'] - r.loc[('lab', 'all'), 'mean_off']))
        ds = data_scores()
        res['B1']['matrices'] = bool(res['B1']['rho_dist_on_change'] >= 0.03 or (
            res['B1']['mean_on_change'] >= (ds['mean_on'] - float(r.loc[('lab', 'all'), 'mean_on'])) / 3.0
            and res['B1']['mean_off_change'] >= (ds['mean_off'] - float(r.loc[('lab', 'all'), 'mean_off'])) / 3.0))
    if ('B1', 'S') in r.index and ('B4', 'S') in r.index:
        gS, gC = float(r.loc[('B1', 'S'), 'gain_z8_pattern']), float(r.loc[('B4', 'S'), 'gain_z8_pattern'])
        res['control'] = dict(B1_S_gain=gS, B4_S_gain=gC, size_not_structure=bool(gS > 0 and gC >= CONTROL_SHARE * gS))
    if not (out / 'population_realism.csv').exists() or ('B1', 'all') not in r.index:
        population_realism_table(rt, runs)
    (out / 'verdict.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    write_archive_tables(rt, runs)
    comparison_table(readings, rt, res)
    figure(readings, rt)
    print(json.dumps({k: v for k, v in res.items() if k in ('B1', 'control', 'mixture', 'groups')}, indent=1, default=float), flush=True)
    return res


def run_realism_parts(o: Dict[str, tuple]) -> dict:
    """A run's response rates, per-odor response fractions (APL16 order) and odors-per-KC distributions, both
    conditions."""
    from kcq import apl_effect as A
    row = {}
    for c in ('on', 'off'):
        r = o[c][0]
        row[f'rate_{c}'] = float(r.to_numpy(dtype=float).mean())
        row[f'frac_{c}'] = H.abbrev_columns(r.astype(float) > 0.5).mean().reindex(H.APL16_ODORS).to_numpy(dtype=float)
        row[f'opk_{c}'] = A.odors_per_kc(r)
    return row


def population_realism(parts: Sequence[dict]) -> dict:
    """The realism bands read on a population (the rule of 16:10: two spikes against the KC data in baseline-noise
    units): the mean response rates over its model flies and their ratio, the response rate ρ of the mean per-odor
    response fractions with APL on, the odors-per-KC distances of the mean distributions, and the flag."""
    from kcq import apl_effect as A
    row = {}
    for c in ('on', 'off'):
        row[f'rate_{c}'] = float(np.mean([q[f'rate_{c}'] for q in parts]))
        frac = pd.Series(np.mean([q[f'frac_{c}'] for q in parts], axis=0), index=H.APL16_ODORS)
        row[f'rate_rho_{c}'] = H.spearman(frac, A.frac_data(c, 'noise'))
        d = np.mean([q[f'opk_{c}'] for q in parts], axis=0)
        row[f'opk_{c}_tv'] = float(0.5 * np.abs(d - A.opk_data(c, 'noise')).sum())
    row['ratio'] = row['rate_off'] / row['rate_on']
    row['realism_flag'] = A.realism_flag(row['rate_on'], row['rate_off'], row['rate_rho_on'], row['opk_on_tv'],
                                         row['opk_off_tv'], units='noise')
    return row


def population_realism_table(rt: pd.DataFrame, runs: Dict[str, Dict[str, tuple]]) -> pd.DataFrame:
    """population_realism.csv: `population_realism` of the wirings at the lab's input (each, and S, T, all) and, per
    group, of each wiring's nine model flies and of S, T and all."""
    lab = {w: run_realism_parts(E.saved_fit(w)['saved']) for w in ALL_WIRINGS}
    parts = {t: run_realism_parts(runs[t]) for t in rt.tag}
    pops = {'S': S_WIRINGS, 'T': T_WIRINGS, 'all': ALL_WIRINGS}
    rows = [dict(group='lab', unit=w, kind='wiring', n_model_flies=1, **population_realism([lab[w]])) for w in ALL_WIRINGS]
    rows += [dict(group='lab', unit=n, kind='population', n_model_flies=len(ws), **population_realism([lab[w] for w in ws]))
             for n, ws in pops.items()]
    n9 = len(flies()['labels'])
    for g, gt in rt.assign(group=[group_of(s) for s in rt.stage]).groupby('group', sort=False):
        for w, wt in gt.groupby('wiring', sort=False):
            if len(wt) == n9:
                rows.append(dict(group=g, unit=w, kind='wiring_nine_flies', n_model_flies=n9,
                                 **population_realism([parts[t] for t in wt.tag])))
        for n, ws in pops.items():
            sub = gt[gt.wiring.isin(ws)]
            if len(sub) == n9 * len(ws):
                rows.append(dict(group=g, unit=n, kind='population', n_model_flies=len(sub),
                                 **population_realism([parts[t] for t in sub.tag])))
    out = pd.DataFrame(rows)
    out.to_csv(out_dir() / 'population_realism.csv', index=False)
    return out


def write_archive_tables(rt: pd.DataFrame, runs: Dict[str, Dict[str, tuple]]) -> None:
    """model_pair_distances_long.csv (raw spike counts, the harness's primary model RDM) and variants.csv (absolute
    parameters) of every run, for kcq/apl_effect.py and the model-output archive."""
    out = out_dir()
    parts, vrows = [], []
    for r in rt.itertuples():
        o = runs[r.tag]
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
        vrows.append(dict(tag=r.tag, variant=r.tag, part='B', stage=r.stage, input_variant=r.variant, fly=r.fly,
                          wiring=r.wiring, reading=r.reading, base='bouton', k=E.BREADTH_K, fit_folder=r.fit_folder,
                          fit_tag=r.fit_tag, kc_kc_inh=r.wKCKC_scale,
                          **{k: getattr(r, k) for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                                        'thr_factor', 'apl_factor', 'rate_on', 'rate_off')},
                          orn_input=r.orn_input, subdir='B'))
    pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    pd.DataFrame(vrows).to_csv(out / 'variants.csv', index=False)


def _f(x, spec='.4f') -> str:
    return E._fmt(x if x is None else float(x), spec)


def comparison_table(readings: pd.DataFrame, rt: pd.DataFrame, res: dict) -> None:
    """comparison_table.md: the populations against the wirings at the lab's input and the KC data, the per-wiring paired
    gains, realism per group, the verdict."""
    out = out_dir()
    r = readings.set_index(['group', 'unit', 'readout'])
    ds = data_scores()
    cols = [('lab', 'S'), ('lab', 'T'), ('lab', 'all')] + [(g, u) for g in ('B1', 'B2', 'B3', 'B4', 'B5')
                                                           for u in ('S', 'T', 'all') if (g, u, 'cal') in r.index]
    names = {'lab': 'lab input', 'B1': 'repro, held', 'B2': 'as measured, held', 'B3': 'gain-free, held',
             'B4': 'permuted, held', 'B5': 'repro, free'}
    stats = [('pattern part, readout, all recordings (lower is better)', 'cal', 'z8_pattern', '.4f', None),
             ('APL-effect loss, readout', 'cal', 'z8_loss', '.4f', None),
             ('mean APL effect, readout', 'cal', 'mean_effect', '+.3f', ds['mean_effect']),
             ('ρ on the per-pair APL effect, readout', 'cal', 'z8_rho', '.2f', None),
             ('pattern part, readout, six typical TNT recordings', 'cal', 'z6_pattern', '.4f', None),
             ('pattern part, readout, noise-free effect', 'cal', 'noise_free_pattern', '.4f', None),
             ('pattern part, readout, baseline-noise units', 'cal', 'noise8_pattern', '.4f', None),
             ('pattern part, readout, the 66 megamat pairs', 'cal', 'z8_pattern_megamat', '.4f', None),
             ('APL-effect loss, raw counts', 'raw', 'z8_loss', '.4f', None),
             ('mean APL effect, raw counts', 'raw', 'mean_effect', '+.3f', ds['mean_effect']),
             ('mean distance, APL on / control flies, readout', 'cal', 'mean_on', '.3f', ds['mean_on']),
             ('mean distance, APL off / APL output blocked (TNT), readout', 'cal', 'mean_off', '.3f', ds['mean_off']),
             ('ρ on odor-odor distances, APL on, readout', 'cal', 'rho_dist_on', '.2f', None),
             ('ρ on odor-odor distances, APL off, readout', 'cal', 'rho_dist_off', '.2f', None),
             ('APL effect, four nearest ester pairs, readout', 'cal', 'nearest_effect', '+.3f', ds['nearest_effect']),
             ('APL effect, six near pairs, readout', 'cal', 'six_effect', '+.3f', ds['six_effect'])]
    lines = ['# E85: the ORN input as a population of flies', '',
             f"Generated by experiments/E85_orn_fly_population.py at {res['queue_commit'][:7]} from {res['n_runs']} runs "
             f"({', '.join(f'{g} {n}' for g, n in res['groups'].items())}).", '',
             '## Populations', '',
             '| statistic | KC data | ' + ' | '.join(f'{names[g]}<br>{u}' for g, u in cols) + ' |',
             '|---|---:|' + '---:|' * len(cols)]
    for label, ro, c, spec, dv in stats:
        cells = [_f(r.loc[(g, u, ro), c], spec) if (g, u, ro) in r.index else '–' for g, u in cols]
        lines.append(f'| {label} | {_f(dv, spec) if dv is not None else "–"} | ' + ' | '.join(cells) + ' |')
    lines += ['', '## The repro flies on each wiring (B1, readout, all recordings)', '',
              '| wiring | lab input | nine flies | gain | single flies, min to max |', '|---|---:|---:|---:|---:|']
    mf = readings[(readings.group == 'B1') & (readings.kind == 'model_fly') & (readings.readout == 'cal')]
    mf = mf.merge(rt[['tag', 'wiring']], left_on='unit', right_on='tag')
    for w in ALL_WIRINGS:
        if ('B1', w, 'cal') in r.index:
            s = mf[mf.wiring == w].z8_pattern
            lines.append(f"| {w} | {_f(r.loc[('lab', w, 'cal'), 'z8_pattern'])} | {_f(r.loc[('B1', w, 'cal'), 'z8_pattern'])} | "
                         f"{_f(r.loc[('B1', w, 'cal'), 'gain_z8_pattern'], '+.4f')} | {_f(s.min())} to {_f(s.max())} |")
    lines += ['', '## Realism per group (every run)', '',
              '| group | runs | response rate on | off / on | response rate ρ, on | odors per KC TV, on / off | largest count on / off | thr factor | APL factor | flagged |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for g, gt in rt.assign(group=[group_of(s) for s in rt.stage]).groupby('group', sort=False):
        fl = int((gt.realism_flag.fillna('').astype(str).str.strip() != '').sum())
        lines.append(f'| {names.get(g, g)} | {len(gt)} | {gt.rate_on.min():.3f} to {gt.rate_on.max():.3f} | '
                     f'{gt.ratio.min():.2f} to {gt.ratio.max():.2f} | {gt.rate_rho_on.min():.2f} to {gt.rate_rho_on.max():.2f} | '
                     f'{gt.opk_on_tv.max():.3f} / {gt.opk_off_tv.max():.3f} | {gt.max_count_on.max():.0f} / {gt.max_count_off.max():.0f} | '
                     f'{gt.thr_factor.min():.3f} to {gt.thr_factor.max():.3f} | {gt.apl_factor.min():.3f} to {gt.apl_factor.max():.3f} | {fl} |')
    pr = pd.read_csv(out / 'population_realism.csv')
    pr = pr[pr.kind != 'wiring']
    lines += ['', '## Realism per population (the bands read on the mean over its model flies)', '',
              '| group | unit | model flies | response rate on | off / on | response rate ρ, on | odors per KC TV, on / off | flag |',
              '|---|---|---:|---:|---:|---:|---:|---|']
    for q in pr.itertuples():
        fl = q.realism_flag if isinstance(q.realism_flag, str) and q.realism_flag.strip() else '–'
        lines.append(f'| {names.get(q.group, q.group)} | {q.unit} | {q.n_model_flies} | {q.rate_on:.3f} | {q.ratio:.2f} | '
                     f'{q.rate_rho_on:.2f} | {q.opk_on_tv:.3f} / {q.opk_off_tv:.3f} | {fl} |')
    pa = json.loads((out / 'part_a.json').read_text())
    lines += ['', '## Part A: the input flies', '',
              f"Scaling: {pa['scaling']['method']}, slope {pa['scaling']['slope']:.2f}, R² of the lab input on the mean over flies "
              f"{pa['scaling']['r2']:.3f} ({pa['scaling']['note']}). Reliability r_g: median {pa['reliability']['r_g_median']:.2f}, "
              f"{pa['reliability']['r_g_min']:.2f} to {pa['reliability']['r_g_max']:.2f}. Near-pair check: ratio "
              f"{pa['near_check']['ratio_at_c1']:.3f} at factor 1, common factor {pa['near_check']['c']:g}. SD across flies of "
              f"log2(fly / mean), median: {pa['spread']['median_sd_log2_all']:.2f} (all entries), "
              f"{pa['spread']['median_sd_log2_above10']:.2f} (mean above 10 % of the largest), "
              f"{pa['spread']['median_sd_log2_above30']:.2f} (above 30 %).", '',
              '| searched input at 0.6 | entries within 1 SD (log2) | within 2 SD | within 1 SD (spike deltas) | within 2 SD | share in the flies\' span | label permutations | odors above the permutations\' 95th percentile |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
    ss = pd.read_csv(out / 'searched_inputs_summary.csv')
    for q in ss.itertuples():
        lines.append(f'| {q.input} | {q.within1_log2:.1%} | {q.within2_log2:.1%} | {q.within1_spike_deltas:.1%} | '
                     f'{q.within2_spike_deltas:.1%} | {q.span_share:.1%} | {q.span_share_perm_mean:.1%} | {q.n_odors_above_perm_p95} of {q.n_odors} |')
    lines += ['', '## Verdict', '', '```', json.dumps({k: v for k, v in res.items() if k in ('B1', 'control', 'mixture')},
                                                      indent=1, default=float), '```', '']
    (out / 'comparison_table.md').write_text('\n'.join(lines))


def figure(readings: pd.DataFrame, rt: pd.DataFrame) -> None:
    """E85_populations.png: per wiring, the pattern part under the readout against all recordings at the lab's input,
    for each repro model fly and for the nine flies' population (B1); and the near pairs' APL effect."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r = readings[readings.readout == 'cal']
    mf = r[(r.group == 'B1') & (r.kind == 'model_fly')].merge(rt[['tag', 'wiring']], left_on='unit', right_on='tag')
    if mf.empty:
        return
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    ax = axs[0]
    for i, w in enumerate(ALL_WIRINGS):
        s = mf[mf.wiring == w].z8_pattern
        ax.scatter(np.full(len(s), i) + np.linspace(-0.2, 0.2, len(s)), s, s=10, color='0.6',
                   label='single model flies' if i == 0 else None)
        lab = r[(r.group == 'lab') & (r.unit == w)].z8_pattern
        ax.scatter([i], lab, marker='_', s=300, color='k', label='lab input' if i == 0 else None)
        nine = r[(r.group == 'B1') & (r.unit == w) & (r.kind == 'wiring_nine_flies')].z8_pattern
        if len(nine):
            ax.scatter([i], nine, marker='o', s=40, color='C0', label='nine input flies' if i == 0 else None)
    ax.set(xticks=range(len(ALL_WIRINGS)), xticklabels=ALL_WIRINGS,
           ylabel='pattern part of the APL-effect loss\n(calibrated readout, all recordings)', title='each wiring')
    ax.legend(fontsize=8)
    ax = axs[1]
    pr = pd.read_csv(out_dir() / 'pairs.csv')
    dd = data_distances()
    de = E._cs(H.apl_effect(dd['on'], dd['off']))
    pops = [('lab', 'all'), ('B1', 'all')]
    x = np.arange(len(de))
    order = np.argsort(de.to_numpy())
    ax.plot(x, de.to_numpy()[order], color='k', label='KC data')
    for g, u in pops:
        p = pr[(pr.group == g) & (pr.unit == u) & (pr.readout == 'cal')].set_index(['odor_a', 'odor_b']).effect
        if len(p):
            ax.plot(x, E._cs(p).reindex(de.index).to_numpy()[order], lw=0.8,
                    label='eight wirings, lab input' if g == 'lab' else 'eight wirings x nine input flies')
    ax.set(xlabel='odor pair (ordered by the KC data)', ylabel='APL effect (APL on minus off)\ncosine distance',
           title='per-pair APL effect, populations')
    ax.legend(fontsize=8)
    fig.savefig(out_dir() / 'E85_populations.png', dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('A', 'identity', 'B', 'C'))
    ap.add_argument('--stage', default=','.join(STAGE_ORDER), help='part B: comma-separated stages (STAGES)')
    a = ap.parse_args()
    t0 = time.time()
    if a.part == 'A':
        part_a()
    elif a.part == 'identity':
        identity()
    elif a.part == 'B':
        stages = [s.strip() for s in a.stage.split(',') if s.strip()]
        bad = set(stages) - set(STAGES)
        if bad:
            raise SystemExit(f'unknown stages {sorted(bad)}')
        part_b(stages)
    elif a.part == 'C':
        part_c()
    print(f'done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
