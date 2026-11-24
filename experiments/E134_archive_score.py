#!/usr/bin/env python3
"""E134 step 2, the archive at the KC somata: every archived variant on the APL16 panel as a population.

The entry is item 2 of experiment_priorities.md ("the models as populations at every recorded stage"), design steps 6
and 7. For every variant of results/apl_effect/by_variant.csv whose APL-on and acute APL-off spike counts are in the
model-output archive (data/model_outputs/<folder>/ or a subfolder, not a withheld_* one, and not a variant of
`kcq.apl_effect.HELD_VARIANTS`) and whose odors include the APL16 panel's 16:
- its KC population with APL on is scored against the control flies and with APL off against APL output blocked (TNT),
  with kcq/population_fit.py: the population error, each view and part in one-fly units, and the energy-distance guard
  (the mean energy distance to the recordings as a multiple of the recordings' mean distance from each other);
- the order of the variants is compared between two random halves of the recordings (each genotype's recordings split
  in two), 200 times: the ρ over variants between a measure's values against one half and against the other, for the
  population error, each view and part, the energy distance and the two ρ's (ρ on odor-odor distances and response
  rate ρ), over all distinct variants and over the distinct variants inside the realism bands (the measure's gate:
  0.8 or more). Against a half, the half is read as the flies (its own one-fly units; `--units-from all` keeps all
  recordings' units);
- over the distinct variants inside the bands, the population error against the calibrated pattern part (the
  APL-effect statistics, a held-out test under the adopted goal), and the best variant and the best family (here an
  experiment folder, by its best variant) by each.

The flies are the lock's populations (`data_prep.stage_population_data.load_population('kc_soma', 'apl16', 'control')`
and `'tnt'`): trial means in baseline-noise units, a response above the lock's criterion. The two ρ's are computed
against the same recordings: ρ on odor-odor distances against the mean over a set's recordings of each recording's
cosine distances across its KCs, response rate ρ against the mean of their per-odor response rates (so they are not
the archive's `dist_rho_*` columns, which are against Yang's z-scored matrices). The model is read as the lock says
(`--reading`), with kcq/stage_readings.py's KC-soma readings (`SomaReading`, `soma_readings()`), as the family's runs
are read: its default ('curve_noise': one monotone curve set on the APL-on counts and applied unchanged to APL off, and
the flies' measurement noise) or an alternative ('raw': the spike counts with a response at two or more spikes, the
planner's probe; ...). Which reading is the stage's own is the lock's to fix. A reading's curve is set against all
recordings, also where a half is scored. `--parts probe` reads the probe's parts only (strength by odor on the
response-rate part). The energy distance maps every population's pooled values rank for rank onto the pooled values of
the condition's recordings (`--energy-reference pooled`; `none` reads the values as they are). Variants with identical
spike counts (the same fit saved by several experiments) count once in the consistency and the trade-off
(`duplicate`).

No model is read before the lock is committed: `main` refuses to run until the lock (experiments/E134_lock.json) is
in git.

Writes <shared>/results/E134_population_baseline/ (the `<reading>` suffix also carries a non-default `--parts`):
    archive_scores_<reading>.csv          one row per folder x tag x condition (on / off): the Score columns of
                                          kcq.population_fit, the energy guard, the two ρ's against all recordings,
                                          inside_bands, duplicate
    archive_consistency_<reading>.csv     per condition, measure and variant set: the split-half ρ (mean, 5 and 95 %)
    archive_tradeoff_<reading>.csv        per distinct variant inside the bands: population errors with APL on and off,
                                          the calibrated and raw pattern parts, and the ranks by each
    archive_tradeoff_best_<reading>.csv   per measure: the best variant and the best family (folder) by it
    archive_provenance_<reading>.json
and the per-variant parts to outputs/E134_population_baseline/archive_views_<reading>.parquet.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_archive_score.py [--reading curve_noise|raw|...] [--parts all|probe]
        [--n-splits 200]
        [--jobs 8] [--seed 0] [--units-from half|all] [--energy-reference pooled|none] [--no-energy] [--no-calibrated]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from kcq import harness as H
from kcq import apl_effect as AE
from kcq import population_fit as P
from kcq import stage_readings as SR

EXPERIMENT = 'E134_population_baseline'
RESULTS = H.shared_dir() / 'results' / EXPERIMENT
OUTPUTS = H.HERE / 'outputs' / EXPERIMENT
ODORS: Tuple[str, ...] = tuple(H.APL16_ODORS)
CONDITIONS = {'on': ('on', 'control'), 'off': ('acute_off', 'tnt')}   # model condition -> (archive name, fly cohort)
N_SPLITS = 200
SEED = 0
PROBE_PARTS = ('strength/rate', 'sharing/who', 'sharing/how_strongly', 'breadth/bins')
EFFECT_COLS = tuple(f'{AE.V}__apl_effect_{k}' for k in ('rmse', 'rmse_centered', 'mean', 'spearman', 'pearson'))
ENERGY_MAX_UNITS = P.ENERGY_MAX_UNITS
BY_VARIANT = AE.OUT / 'by_variant.csv'
# a reading maps one condition's spike counts (KCs x odors) to its population: `reading(counts_on)` for APL on,
# `reading(counts_off, condition='off', curve_counts=counts_on)` for APL off (stage_readings.SomaReading)
Reading = Callable[..., P.Population]
DEFAULT_READING = SR.DEFAULT_READINGS['kc_soma'].name


def readings() -> Dict[str, Reading]:
    """The readings `--reading` may name: stage_readings' KC-soma readings (`SR.soma_readings()`)."""
    return SR.soma_readings()


def lock_committed() -> dict:
    """The committed lock (experiments/E134_lock.json, tracked by git with no uncommitted change); raises RuntimeError
    otherwise (no model is read before the lock). Uses the data module's `require_lock` once it is written."""
    from data_prep import stage_population_data as S
    try:
        return S.require_lock()
    except NotImplementedError:
        pass
    path = Path(S.LOCK_COPY)
    rel = str(path.relative_to(H.HERE)) if path.is_relative_to(H.HERE) else str(path)
    git = lambda *a: subprocess.run(['git', '-C', str(H.HERE), *a], capture_output=True).returncode == 0
    if not (path.exists() and git('ls-files', '--error-unmatch', rel) and git('diff', '--quiet', 'HEAD', '--', rel)):
        raise RuntimeError(f'the E134 lock ({path}) is not committed: no model is read before it')
    return json.loads(path.read_text())


def fly_populations(cohort: str) -> List[P.Population]:
    """The APL16 KC somata of one cohort ('control' or 'tnt') from the lock's data module, one population per
    recording, columns in the APL16 order."""
    from data_prep import stage_population_data as S
    rec = S.load_population('kc_soma', 'apl16', cohort)
    if set(rec.odors) != set(ODORS):
        raise ValueError(f'the APL16 recordings have odors {rec.odors}, not the panel {ODORS}')
    return [P.Population(p.values[list(ODORS)], p.responds[list(ODORS)], name=p.name, graded=p.graded)
        for p in P.fly_populations(rec)]


def archived_counts(folder: str, archive: Optional[Path] = None) -> Dict[str, Path]:
    """{tag: directory} of a folder's archived variants with APL-on and acute APL-off spike counts, the folder itself
    first and then its subfolders in sorted order (as analyses/apl_effect_metrics.py), withheld_* folders left out."""
    root = Path(archive or AE.ARCHIVE) / folder
    out: Dict[str, Path] = {}
    if not root.is_dir():
        return out
    for q in sorted(root.rglob('*_on_spike_counts.parquet'), key=lambda q: (q.parent != root, str(q.parent), q.name)):
        if any(AE.WITHHELD_DIR.match(x) for x in q.relative_to(root).parts[:-1]):
            continue
        tag = q.name[:-len('_on_spike_counts.parquet')]
        if tag not in out and (q.parent / f'{tag}_acute_off_spike_counts.parquet').exists():
            out[tag] = q.parent
    return out


def archived_variants(by_variant: Optional[pd.DataFrame] = None, archive: Optional[Path] = None) -> pd.DataFrame:
    """Every variant of results/apl_effect/by_variant.csv with archived APL-on and acute APL-off spike counts: columns
    exp, folder, tag, archive_dir, realism_flag, inside_bands, rate_on, rate_off, and the APL-effect statistics
    by_variant.csv carries (`EFFECT_COLS`)."""
    bv = pd.read_csv(BY_VARIANT, dtype={'tag': str}) if by_variant is None else by_variant
    rows = []
    for folder, g in bv.groupby('folder', sort=False):
        where = archived_counts(folder, archive)
        for _, r in g.iterrows():
            tag = str(r['tag'])
            if AE.held(folder, tag) or tag not in where:
                continue
            flag = r.get('realism_flag', np.nan)
            row = dict(exp=r.get('exp', ''), folder=folder, tag=tag, archive_dir=str(where[tag]),
                realism_flag='' if pd.isna(flag) else str(flag), rate_on=r.get('rate_on', np.nan),
                rate_off=r.get('rate_off', np.nan))
            row['inside_bands'] = row['realism_flag'] == ''
            row.update({c: r.get(c, np.nan) for c in EFFECT_COLS})
            rows.append(row)
    cols = ['exp', 'folder', 'tag', 'archive_dir', 'realism_flag', 'inside_bands', 'rate_on', 'rate_off', *EFFECT_COLS]
    return pd.DataFrame(rows, columns=cols)


def read_counts(archive_dir: Path, tag: str, name: str) -> Optional[pd.DataFrame]:
    """A variant's spike counts (KCs x odors, odors abbreviated, in the APL16 order), or None if the file is missing
    or lacks one of the APL16 panel's odors."""
    p = Path(archive_dir) / f'{tag}_{name}_spike_counts.parquet'
    if not p.exists():
        return None
    df = H.abbrev_columns(pd.read_parquet(p))
    if not set(ODORS) <= set(df.columns):
        return None
    return df[list(ODORS)].astype(float)


# worker state of `variant_views` (set by `_init`)
_W: dict = {}


def _init(reading: str, parts: Optional[Sequence[str]], edges, energy: Optional[dict], readers: Optional[dict]) -> None:
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    _W.clear()
    _W.update(reading=(readers or readings())[reading], parts=parts, edges=edges, energy=energy)


def _variant_job(job: Tuple[str, str, str]) -> dict:
    folder, tag, where = job
    out = dict(folder=folder, tag=tag, conds={})
    h = hashlib.sha256()
    counts = {}
    for cond, (name, _) in CONDITIONS.items():
        counts[cond] = read_counts(Path(where), tag, name)
        if counts[cond] is None:
            return dict(folder=folder, tag=tag, error=f'no APL16 {name} counts')
        h.update(np.ascontiguousarray(counts[cond].to_numpy()).tobytes())
    for cond in CONDITIONS:
        if cond == 'on':
            pop = _W['reading'](counts['on'])
        else:
            pop = _W['reading'](counts[cond], condition=cond, curve_counts=counts['on'])
        rec = dict(views=P.population_views(pop, ODORS, _W['edges'], _W['parts']),
            rates=pop.responds.mean(axis=0).to_numpy(dtype=float),
            dists=P.cosine_pairs(pop.values.to_numpy(dtype=float)))
        e = (_W['energy'] or {}).get(cond)
        if e is not None:
            v = pop.values if e['reference'] is None else P.rank_map(pop.values, e['reference'])
            x = P.energy_sample(v.to_numpy(dtype=float), e['max_units'], e['seed'])
            exx = P.mean_self_distance(x)
            rec['energy'] = np.array([P.energy_from_samples(x, y, exx, eyy) for y, eyy in zip(e['samples'], e['self'])])
        out['conds'][cond] = rec
    out['hash'] = h.hexdigest()
    return out


def energy_setup(flies: Dict[str, List[P.Population]], reference: str = 'pooled', max_units: int = ENERGY_MAX_UNITS,
        seed: int = SEED) -> Dict[str, dict]:
    """Per condition, what the energy guard needs from the recordings: the reference of the rank map (`reference`
    'pooled': the pooled values of the condition's recordings; 'none': no map), each recording's drawn units and their
    mean within-set distance, and the recordings' distances from each other."""
    out = {}
    for cond, pops in flies.items():
        ref = np.concatenate([p.values.to_numpy(dtype=float).ravel() for p in pops]) if reference == 'pooled' else None
        samples = []
        for p in pops:
            v = p.values if ref is None else P.rank_map(p.values, ref)
            samples.append(P.energy_sample(v.to_numpy(dtype=float), max_units, seed))
        selfd = [P.mean_self_distance(x) for x in samples]
        n = len(samples)
        M = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                M[i, j] = M[j, i] = P.energy_from_samples(samples[i], samples[j], selfd[i], selfd[j])
        out[cond] = dict(reference=ref, samples=samples, self=selfd, max_units=max_units, seed=seed,
            fly_matrix=pd.DataFrame(M, index=[p.name for p in pops], columns=[p.name for p in pops]))
    return out


def variant_views(variants: pd.DataFrame, reading: str = DEFAULT_READING, parts: Optional[Sequence[str]] = None,
        breadth_edges=P.APL16_BREADTH_EDGES, energy: Optional[Dict[str, dict]] = None, jobs: int = 1,
        readers: Optional[Dict[str, Reading]] = None) -> Tuple[pd.DataFrame, Dict[str, dict]]:
    """Each variant's parts per condition, its per-odor response rates and pair distances (for the two ρ's), and its
    energy distance to each recording.

    Args:
        variants: `archived_variants`' table.
        reading: a name of `readings()` (or of `readers`).
        parts: the parts read (default all).
        breadth_edges: the breadth bins.
        energy: `energy_setup`'s table, or None for no energy guard.
        jobs: worker processes.
        readers: the readings by name (default `readings()`; must be picklable for jobs > 1).

    Returns:
        (the variants that could be read, with `duplicate` and `key` columns, and those left out with their reason in
        `skipped`; {key: {condition: dict(views, rates, dists[, energy])}}, key 'folder/tag').
    """
    work = [(r.folder, r.tag, r.archive_dir) for r in variants.itertuples()]
    e = None if energy is None else {c: {k: v for k, v in d.items() if k != 'fly_matrix'} for c, d in energy.items()}
    init = (reading, parts, breadth_edges, e, readers)
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs, initializer=_init, initargs=init) as ex:
            res = list(ex.map(_variant_job, work, chunksize=4))
    else:
        _init(*init)
        res = [_variant_job(w) for w in work]
    seen, keep, data, skipped = set(), [], {}, []
    for (_, r), out in zip(variants.iterrows(), res):
        if 'error' in out:
            skipped.append(dict(folder=r.folder, tag=r.tag, why=out['error']))
            continue
        key = f'{r.folder}/{r.tag}'
        keep.append(dict(r, key=key, duplicate=out['hash'] in seen))
        seen.add(out['hash'])
        data[key] = out['conds']
    meta = pd.DataFrame(keep, columns=[*variants.columns, 'key', 'duplicate'])
    meta.attrs['skipped'] = skipped
    return meta, data


def _rowwise_spearman(M: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Spearman's ρ of each row of M with d, over the entries both define."""
    out = np.full(len(M), np.nan)
    okd = np.isfinite(d)
    clean = np.isfinite(M).all(axis=1) & okd.all()
    if clean.any():
        rm = stats.rankdata(M[clean], axis=1)
        rd = stats.rankdata(d)
        rm = rm - rm.mean(axis=1, keepdims=True)
        rd = rd - rd.mean()
        with np.errstate(invalid='ignore', divide='ignore'):
            out[clean] = (rm @ rd) / (np.sqrt((rm ** 2).sum(axis=1)) * np.sqrt((rd ** 2).sum()))
    for i in np.flatnonzero(~clean):
        ok = np.isfinite(M[i]) & okd
        if ok.sum() > 2:
            out[i] = stats.spearmanr(M[i, ok], d[ok]).statistic
    return out


def fly_summaries(flies: Sequence[P.Population]) -> Tuple[np.ndarray, np.ndarray]:
    """Per recording, its per-odor response rates (recordings x odors) and its pair distances (recordings x pairs,
    cosine across its units), on the APL16 odors."""
    rates = np.array([p.responds[list(ODORS)].mean(axis=0).to_numpy(dtype=float) for p in flies])
    dists = np.array([P.cosine_pairs(p.values[list(ODORS)].to_numpy(dtype=float)) for p in flies])
    return rates, dists


def _measures(arrays: Dict[str, np.ndarray], ref: P.FlyReference, rates: np.ndarray, dists: np.ndarray,
        fly_rates: np.ndarray, fly_dists: np.ndarray, energy: Optional[np.ndarray], fly_energy: Optional[np.ndarray],
        idx: np.ndarray) -> pd.DataFrame:
    """Every measure of every variant against the recordings `idx` (their reference `ref`)."""
    out = P.score_arrays(arrays, ref)
    out['dist_rho'] = _rowwise_spearman(dists, np.nanmean(fly_dists[idx], axis=0))
    out['rate_rho'] = _rowwise_spearman(rates, np.nanmean(fly_rates[idx], axis=0))
    if energy is not None:
        to = np.nanmean(energy[:, idx], axis=1)
        sub = fly_energy[np.ix_(idx, idx)]
        ff = sub[np.triu_indices(len(idx), k=1)].mean() if len(idx) > 1 else np.nan
        out['energy_to_flies'] = to
        out['energy_ratio'] = to / ff
    return out


def _condition_arrays(meta: pd.DataFrame, data: Dict[str, dict], cond: str, ref: P.FlyReference):
    keys = list(meta['key'])
    arrays = P.stack_views([data[k][cond]['views'] for k in keys], ref)
    rates = np.array([data[k][cond]['rates'] for k in keys]).reshape(len(keys), len(ODORS))
    dists = np.array([data[k][cond]['dists'] for k in keys]).reshape(len(keys), -1)
    energy = np.array([data[k][cond]['energy'] for k in keys]) if keys and 'energy' in data[keys[0]][cond] else None
    return arrays, rates, dists, energy


def score_variants(meta: pd.DataFrame, data: Dict[str, dict], flies: Dict[str, List[P.Population]],
        refs: Dict[str, P.FlyReference], energy: Optional[Dict[str, dict]] = None) -> pd.DataFrame:
    """Every variant's measures against all recordings of its condition: one row per variant and condition."""
    out = []
    for cond, ref in refs.items():
        arrays, rates, dists, en = _condition_arrays(meta, data, cond, ref)
        fr, fd = fly_summaries(flies[cond])
        fe = None if energy is None or en is None else energy[cond]['fly_matrix'].to_numpy()
        m = _measures(arrays, ref, rates, dists, fr, fd, en if fe is not None else None, fe, np.arange(len(ref.flies)))
        m.index = meta.index
        base = meta[['exp', 'folder', 'tag', 'inside_bands', 'duplicate', 'realism_flag']]
        out.append(pd.concat([base.assign(condition=cond), m], axis=1))
    return pd.concat(out, ignore_index=True)


def measure_columns(df: pd.DataFrame) -> List[str]:
    """The measures compared in the consistency: the population errors, views, parts, the two ρ's and the energy
    distance (columns `_measures` writes that are not all NaN)."""
    skip = {'n_missing', 'energy_to_flies'}
    cols = [c for c in df.columns if c not in skip and (c.startswith(('error', 'view_', 'part_'))
        or c in ('dist_rho', 'rate_rho', 'energy_ratio'))]
    return [c for c in cols if df[c].notna().any()]


def consistency(meta: pd.DataFrame, data: Dict[str, dict], flies: Dict[str, List[P.Population]],
        refs: Dict[str, P.FlyReference], energy: Optional[Dict[str, dict]] = None, n_splits: int = N_SPLITS,
        seed: int = SEED, units_from: str = 'half') -> pd.DataFrame:
    """The split-half consistency of every measure over the variants.

    Args:
        meta, data: `variant_views`' output.
        flies: condition -> the recordings.
        refs: condition -> the recordings' reference (all of them).
        energy: `energy_setup`'s table, for the energy distance's consistency.
        n_splits: random splits of each condition's recordings in two.
        seed: of the splits.
        units_from: 'half' to score against a half in that half's own one-fly units (the half read as the flies),
            'all' to keep all recordings' one-fly units.

    Returns:
        One row per condition, measure and variant set ('all' and 'inside_bands', distinct variants only):
        split_half_rho (mean over splits), lo, hi (5 and 95 %), n_variants, n_splits.
    """
    sets = {'all': ~meta['duplicate'].to_numpy(bool),
            'inside_bands': ~meta['duplicate'].to_numpy(bool) & meta['inside_bands'].to_numpy(bool)}
    rows = []
    for ci, (cond, ref) in enumerate(refs.items()):
        arrays, rates, dists, en = _condition_arrays(meta, data, cond, ref)
        fr, fd = fly_summaries(flies[cond])
        fe = None if energy is None or en is None else energy[cond]['fly_matrix'].to_numpy()
        res: Dict[Tuple[str, str], List[float]] = {}
        for a, b in P.split_halves(len(ref.flies), n_splits, seed + ci):
            halves = [_measures(arrays, P.subset_reference(ref, h, keep_units=units_from == 'all'), rates, dists, fr, fd,
                en if fe is not None else None, fe, h) for h in (a, b)]
            for name, mask in sets.items():
                for col in measure_columns(halves[0]):
                    x, y = halves[0][col].to_numpy()[mask], halves[1][col].to_numpy()[mask]
                    ok = np.isfinite(x) & np.isfinite(y)
                    r = stats.spearmanr(x[ok], y[ok]).statistic if ok.sum() > 2 else np.nan
                    res.setdefault((name, col), []).append(r)
        for (name, col), v in res.items():
            v = np.asarray(v, dtype=float)
            rows.append(dict(condition=cond, measure=col, variant_set=name, split_half_rho=np.nanmean(v),
                lo=np.nanpercentile(v, 5), hi=np.nanpercentile(v, 95), n_variants=int(sets[name].sum()),
                n_splits=n_splits))
    return pd.DataFrame(rows)


TRADEOFF_MEASURES = ('error_on', 'error_off', 'cal_pattern_part', 'raw_pattern_part')


def tradeoff(scores: pd.DataFrame, variants: pd.DataFrame, calibrated: Optional[pd.DataFrame] = None
        ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Over the distinct variants inside the bands: the population errors with APL on and off beside the calibrated
    pattern part (and the raw one), with the best variant and the best family (experiment folder) by each.

    Args:
        scores: `score_variants`' table.
        variants: the read variants (`variant_views`' meta).
        calibrated: `calibrated_pattern_parts`' table (folder, tag, cal_pattern_part), or None.

    Returns:
        (per variant: folder, tag, the measures and their ranks, 1 the lowest; per measure: the best variant, its
        value, the best folder by its lowest variant, and the Spearman ρ of the measure with the calibrated pattern
        part over the variants).
    """
    keep = variants[variants['inside_bands'].astype(bool) & ~variants['duplicate'].astype(bool)]
    t = keep[['exp', 'folder', 'tag']].copy()
    for cond in CONDITIONS:
        s = scores[scores.condition == cond].set_index(['folder', 'tag'])['error']
        t[f'error_{cond}'] = [s.get((f, g), np.nan) for f, g in zip(t.folder, t.tag)]
    t['raw_pattern_part'] = keep[f'{AE.V}__apl_effect_rmse_centered'].to_numpy(dtype=float)
    t['cal_pattern_part'] = np.nan
    if calibrated is not None and len(calibrated):
        c = calibrated.set_index(['folder', 'tag'])['cal_pattern_part']
        t['cal_pattern_part'] = [c.get((f, g), np.nan) for f, g in zip(t.folder, t.tag)]
    for m in TRADEOFF_MEASURES:
        t[f'rank_{m}'] = t[m].rank(method='min')
    best = []
    for m in TRADEOFF_MEASURES:
        if not t[m].notna().any():
            continue
        i = t[m].idxmin()
        fam = t.groupby('folder')[m].min().dropna().sort_values()
        ok = t[m].notna() & t['cal_pattern_part'].notna()
        rho = stats.spearmanr(t.loc[ok, m], t.loc[ok, 'cal_pattern_part']).statistic if ok.sum() > 2 else np.nan
        best.append(dict(measure=m, best_folder_tag=f"{t.at[i, 'folder']}/{t.at[i, 'tag']}", best_value=t.at[i, m],
            best_family=fam.index[0], best_family_value=fam.iloc[0], second_family=fam.index[1] if len(fam) > 1 else '',
            rho_with_cal_pattern_part=rho, n_variants=int(t[m].notna().sum())))
    return t.reset_index(drop=True), pd.DataFrame(best)


def _cal_job(job: Tuple[str, str, str, float]) -> dict:
    folder, tag, where, target = job
    try:
        on = pd.read_parquet(Path(where) / f'{tag}_on_spike_counts.parquet')
        off = pd.read_parquet(Path(where) / f'{tag}_acute_off_spike_counts.parquet')
        K = H.calibrate_readout(on, 'saturating', target, stat='participation_ratio', responders_min=2.0)
        s = H.apl_effect_stats_readout(on, off, 'saturating', K)
        return dict(folder=folder, tag=tag, K=float(K), cal_pattern_part=s['apl_effect_rmse_centered'],
            cal_loss=s['apl_effect_rmse'], cal_mean=s['apl_effect_mean'])
    except Exception as ex:   # e.g. a variant the readout cannot calibrate
        return dict(folder=folder, tag=tag, error=repr(ex)[:200])


def calibrated_pattern_parts(variants: pd.DataFrame, jobs: int = 1, target: Optional[float] = None) -> pd.DataFrame:
    """Each variant's pattern part under the calibrated readout (as analyses/apl_effect_metrics.py --readout cal reads
    it: the saturating readout's K calibrated on the APL-on counts to the control flies' responders-only participation
    ratio, registry ctrl_participation_ratio_resp, and used unchanged with APL off): folder, tag, K, cal_pattern_part,
    cal_loss, cal_mean (error where it failed)."""
    target = H.data_target('ctrl_participation_ratio_resp') if target is None else target
    work = [(r.folder, r.tag, r.archive_dir, target) for r in variants.itertuples()]
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs, initializer=os.environ.setdefault,
                initargs=('OMP_NUM_THREADS', '1')) as ex:
            res = list(ex.map(_cal_job, work, chunksize=4))
    else:
        res = [_cal_job(w) for w in work]
    return pd.DataFrame(res)


def views_long(meta: pd.DataFrame, data: Dict[str, dict]) -> pd.DataFrame:
    """The per-variant parts as one long table: folder, tag, condition, part, entry, value."""
    rows = []
    for r in meta.itertuples():
        for cond, d in data[r.key].items():
            for part, s in d['views'].items():
                rows.append(pd.DataFrame(dict(folder=r.folder, tag=r.tag, condition=cond, part=part,
                    entry=s.index.astype(str), value=s.to_numpy(dtype=float))))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=['folder', 'tag', 'condition', 'part', 'entry', 'value'])


def _suffix(reading: str, parts: str) -> str:
    return reading if parts == 'all' else f'{reading}_{parts}'


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--reading', default=DEFAULT_READING, help="a KC-soma reading of kcq/stage_readings.py (default "
        "its default, the lock's to fix)")
    ap.add_argument('--parts', choices=('all', 'probe'), default='all',
        help="'probe': strength by odor on its response-rate part only, as the planner's probe read it")
    ap.add_argument('--n-splits', type=int, default=N_SPLITS)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--units-from', choices=('half', 'all'), default='half')
    ap.add_argument('--energy-reference', choices=('pooled', 'none'), default='pooled')
    ap.add_argument('--energy-max-units', type=int, default=ENERGY_MAX_UNITS)
    ap.add_argument('--no-energy', action='store_true')
    ap.add_argument('--no-calibrated', action='store_true')
    a = ap.parse_args(argv)

    lock = lock_committed()
    readers = readings()
    if a.reading not in readers:
        raise SystemExit(f'unknown reading {a.reading!r}; known: {sorted(readers)}')
    parts = None if a.parts == 'all' else list(PROBE_PARTS)
    from data_prep import stage_population_data as S
    try:
        edges = S.breadth_bins(len(ODORS))
    except NotImplementedError:
        edges = P.APL16_BREADTH_EDGES
    flies = {cond: fly_populations(cohort) for cond, (_, cohort) in CONDITIONS.items()}
    refs = {cond: P.fly_reference(f, list(ODORS), edges, parts) for cond, f in flies.items()}
    energy = None if a.no_energy else energy_setup(flies, a.energy_reference, a.energy_max_units, a.seed)

    variants = archived_variants()
    print(f'{len(variants)} archived variants; reading them ({a.jobs} jobs)', flush=True)
    meta, data = variant_views(variants, a.reading, parts, edges, energy, a.jobs)
    print(f"{len(meta)} read ({int((~meta.duplicate).sum())} distinct, "
          f"{int((~meta.duplicate & meta.inside_bands).sum())} inside the bands); {len(meta.attrs['skipped'])} left out",
          flush=True)
    sfx = _suffix(a.reading, a.parts)
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    scores = score_variants(meta, data, flies, refs, energy)
    scores.to_csv(RESULTS / f'archive_scores_{sfx}.csv', index=False)
    views_long(meta, data).to_parquet(OUTPUTS / f'archive_views_{sfx}.parquet')
    cons = consistency(meta, data, flies, refs, energy, a.n_splits, a.seed, a.units_from)
    cons.to_csv(RESULTS / f'archive_consistency_{sfx}.csv', index=False)
    cal = None
    if not a.no_calibrated:
        sel = meta[meta.inside_bands.astype(bool) & ~meta.duplicate.astype(bool)]
        cal = calibrated_pattern_parts(sel, a.jobs)
    per, best = tradeoff(scores, meta, cal)
    per.to_csv(RESULTS / f'archive_tradeoff_{sfx}.csv', index=False)
    best.to_csv(RESULTS / f'archive_tradeoff_best_{sfx}.csv', index=False)
    prov = dict(script='experiments/E134_archive_score.py', argv=list(argv) if argv is not None else None,
        args=vars(a), generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
        lock_sha256=hashlib.sha256(json.dumps(lock, sort_keys=True).encode()).hexdigest(),
        by_variant_sha256=hashlib.sha256(Path(BY_VARIANT).read_bytes()).hexdigest(),
        breadth_bins=[list(map(float, b)) for b in P.breadth_bins(edges)], parts=list(refs['on'].parts),
        n_variants=len(variants), n_read=len(meta), n_distinct=int((~meta.duplicate).sum()),
        n_inside_bands_distinct=int((~meta.duplicate & meta.inside_bands).sum()), skipped=meta.attrs['skipped'][:50],
        n_skipped=len(meta.attrs['skipped']),
        calibrated_errors=[] if cal is None or 'error' not in cal else cal['error'].dropna().tolist()[:20],
        flies={c: [p.name for p in f] for c, f in flies.items()})
    (RESULTS / f'archive_provenance_{sfx}.json').write_text(json.dumps(prov, indent=1, default=str))
    gate = cons[(cons.measure == 'error') & (cons.variant_set == 'inside_bands')]
    print(gate.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
