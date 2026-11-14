#!/usr/bin/env python3
"""Other measures of agreement between a model's per-pair APL effect and the KC data's,
compared on every saved model variant.

The main target (CLAUDE.md, "The main goal for new experiments") is scored by
`kcq.harness.apl_effect_stats`: the APL-effect loss (root mean squared difference over the 120
pairs), its mean part and its pattern part (the loss with both means removed), with ρ on the
per-pair APL effect (Spearman) and Pearson's r beside. This script adds, for the same per-pair
effects (cosine distance on raw spike counts, APL on minus acute APL off; the KC data's control
flies minus APL output blocked (TNT), each the mean over recordings), and for every variant that
kcq/apl_effect.py scores:

- the pattern part split exactly into an amplitude part and a shape part (the Taylor-diagram
  identity, pattern^2 = (sd_m - sd_d)^2 + 2 sd_m sd_d (1 - r), with sd the spread of each side's
  per-pair effect over the pairs and r Pearson's): whether a model's pattern error is a spread of
  its per-pair effects that is too small or too large, or a different shape;
- Kendall's tau beside Spearman's ρ and Pearson's r;
- the pattern part and Pearson's r weighted by each pair's reliability in the data (inverse of
  the squared standard error of the data's per-pair effect over recordings, both genotypes), so
  that pairs the flies measure poorly count less;
- the pattern part split into its per-odor part (the projection of the error onto the additive
  space over odors, a value per odor summed over a pair's two odors) and its pair-specific part
  (the residual), as `kcq.harness.apl_effect_pair_parts` does for one variant.

It then compares the measures: how they correlate over the variants (all, and those inside the
realism bands), how reliably each ranks the variants against two independent halves of the
recordings (random splits of each genotype's recordings; the ρ over variants between a measure's
values against one half and against the other, mean over splits), which variant each calls best
(overall and per experiment) and how much their top 25 overlap. Variants whose per-pair effect is
identical to an earlier one's (the same fit saved by several experiments) count once in the
correlations, the reliabilities and the overlaps (`duplicate` in by_variant_metrics.csv).

Writes <shared>/results/apl_effect_metrics/: by_variant_metrics.csv, metric_correlations.csv
(all variants) and metric_correlations_realistic.csv, split_half_reliability.csv,
best_by_metric.csv, best_per_experiment.csv, top25_overlap.csv, shares.csv (per variant inside
the bands, the mean part's share of the loss squared, the amplitude part's and the per-odor part's
shares of the pattern part squared, and the ratio of the spreads), odor_share_reliability.csv (how
reproducible single odors' shares are between random halves of the recordings: the data's own per-odor
profile and the error of the two models quoted on 2026-11-13, ODOR_SHARE_MODELS), data_pairs.csv (the data's
per-pair effect and its standard error) and provenance.json. Three further tables come from the same
splits as split_half_reliability.csv: split_half_subsets.csv (the same consistency over all distinct
variants, over the distinct variants inside the bands, and within experiments: the mean over experiment
folders with 8 or more distinct variants inside the bands), best_pick_stability.csv (per experiment
folder with 2 or more variants inside the bands, how often a measure's best variant differs between the
halves and from the pattern part's best; of the distinct variants inside the bands, the rank by the
other half of the variant one half calls best) and odor_share_agreement.csv (the r over the 16 odors
between the halves' single-odor shares of each variant's error, median over the distinct variants inside
the bands and over the 25 of them with the lowest pattern part). These follow the read-only check of
2026-11-13 by hal's planner (notes/metrics_rules_scripts/check2b_picks.py and check3_readout.py in the
shared directory).

`--readout cal` reads the same variants through the calibrated readout instead of raw counts (the reading
the gates use since the APL-effect fit on re-drawn wirings, E72): every variant of
results/apl_effect/by_variant.csv whose APL-on and acute APL-off spike counts are in the model-output
archive (data/model_outputs/<folder>/, or a subfolder of it), its saturating readout's K calibrated on its
APL-on counts to the control flies' responders-only participation ratio (registry
ctrl_participation_ratio_resp; responders two or more spikes; as E65, E67 and E83 do) and used unchanged
with APL off, the per-pair APL effect from the readout's cosine distances. Each variant's raw counts from
the archive are checked against by_variant.csv's loss and pattern part (a variant whose archived counts do
not give its saved distances is left out). It writes the same tables with `_cal` before the extension
(by_variant_metrics_cal.csv, with K, K_at_bound and archive_dir; ...; provenance_cal.json), and no
data_pairs.csv (the data side does not change). The realism bands are read on the counts as before.

Usage: ./venv/bin/python analyses/apl_effect_metrics.py [--n-splits 200] [--readout raw|cal] [--jobs 4]
       [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from kcq import harness as H
from kcq import apl_effect as AE

OUT = H.shared_dir() / 'results' / 'apl_effect_metrics'
SEED = 0
TOP_K = 25

# (column, lower is better)
METRICS: Dict[str, bool] = {
    'loss': True, 'mean_part': True, 'pattern_part': True,
    'amplitude_part': True, 'shape_part': True, 'sd_ratio': None,
    'spearman': False, 'pearson': False, 'kendall': False,
    'pattern_part_w': True, 'pearson_w': False,
    'per_odor_part': True, 'pair_specific_part': True,
}
# the measures compared (sd_ratio is described, not ranked: 1 is best, both sides of it are worse)
RANKED = [k for k, v in METRICS.items() if v is not None]


def data_recordings() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """The per-recording pair distances of each genotype (rows recordings, columns the 120 pairs in
    the data effect's order) and the data's per-pair effect (`H.yang_apl_effect()`)."""
    d = H.yang_apl_effect()
    ctrl = H.yang_recording_rdms('TNT3c', 'cosine')[d.index]
    tnt = H.yang_recording_rdms('TNT', 'cosine')[d.index]
    check = ctrl.mean() - tnt.mean()
    if not np.allclose(check.to_numpy(), d.to_numpy(), atol=1e-12):
        raise RuntimeError('the per-recording distances do not give harness.yang_apl_effect()')
    return ctrl, tnt, d


def pair_se(ctrl: pd.DataFrame, tnt: pd.DataFrame) -> np.ndarray:
    """Standard error of the data's per-pair effect (difference of two means over recordings)."""
    return np.sqrt(ctrl.var(ddof=1).to_numpy() / len(ctrl) + tnt.var(ddof=1).to_numpy() / len(tnt))


def model_effects(index: pd.Index) -> Tuple[pd.DataFrame, np.ndarray]:
    """Every scored variant's per-pair effect (rows exp / folder / tag, columns `index`), read as
    kcq/apl_effect.py reads them (the spike-count cosine rows of each results folder's
    model_pair_distances_long.csv, APL on minus acute APL off), variants with all 120 pairs only."""
    rows, mats = [], []
    for folder in sorted(AE.RESULTS.glob('E*_*')):
        m = re.match(r'(E\d+)_', folder.name)
        if not m or not folder.is_dir():
            continue
        for tag, (on, off) in AE.model_pairs(folder).items():
            e = H.apl_effect(on, off).reindex(index)
            if e.isna().any():
                continue
            rows.append(dict(exp=m.group(1), folder=folder.name, tag=tag))
            mats.append(e.to_numpy(dtype=float))
    return pd.DataFrame(rows), np.array(mats)


def _ranks(a: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(stats.rankdata, -1, a)


def _rowcorr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson's r of each row of a (n x p) with b (p,) or with each row of b (n x p)."""
    ac = a - a.mean(-1, keepdims=True)
    bc = b - b.mean(-1, keepdims=True)
    return (ac * bc).sum(-1) / np.sqrt((ac ** 2).sum(-1) * (bc ** 2).sum(-1))


def metrics(M: np.ndarray, d: np.ndarray, w: np.ndarray, X: np.ndarray, kendall: bool = True) -> pd.DataFrame:
    """The measures for every row of M (variants x pairs) against the data effect d.

    Args:
        M: model per-pair effects.
        d: data per-pair effect.
        w: pair weights (normalized to sum 1 here).
        X: the additive design over odors of the pairs (`H.pair_design`).
        kendall: compute Kendall's tau (the slowest; skipped in the split halves).
    """
    n = M.shape[1]
    err = M - d
    mean_err = err.mean(1)
    Mc = M - M.mean(1, keepdims=True)
    dc = d - d.mean()
    ec = Mc - dc
    pattern = np.sqrt((ec ** 2).mean(1))
    sd_m = np.sqrt((Mc ** 2).mean(1))
    sd_d = float(np.sqrt((dc ** 2).mean()))
    r = (Mc @ dc) / (n * sd_m * sd_d)
    amp = (sd_m - sd_d) ** 2
    shape = 2 * sd_m * sd_d * (1 - r)
    # (exact: pattern^2 = amplitude + shape)
    assert np.allclose(pattern ** 2, amp + shape, atol=1e-12)
    rho = _rowcorr(_ranks(M), stats.rankdata(d))
    w = w / w.sum()
    Mw = M - (M @ w)[:, None]
    dw = d - d @ w
    pattern_w = np.sqrt(((Mw - dw) ** 2) @ w)
    r_w = (Mw * dw) @ w / np.sqrt(((Mw ** 2) @ w) * ((dw ** 2) @ w))
    P = X @ np.linalg.pinv(X)
    per_odor = ec @ P.T
    out = dict(loss=np.sqrt((err ** 2).mean(1)), mean_part=np.abs(mean_err), mean_error=mean_err,
        pattern_part=pattern, amplitude_part=amp, shape_part=shape, sd_model=sd_m, sd_ratio=sd_m / sd_d,
        spearman=rho, pearson=r, pattern_part_w=pattern_w, pearson_w=r_w,
        per_odor_part=np.sqrt((per_odor ** 2).mean(1)), pair_specific_part=np.sqrt(((ec - per_odor) ** 2).mean(1)))
    if kendall:
        out['kendall'] = np.array([stats.kendalltau(m, d).statistic for m in M])
    return pd.DataFrame(out)


def _oriented(df: pd.DataFrame) -> pd.DataFrame:
    """Each ranked measure with higher = better (so correlations between measures read as agreement)."""
    return pd.DataFrame({k: (-df[k] if METRICS[k] else df[k]) for k in RANKED})


def split_half_reliability(M: np.ndarray, ctrl: pd.DataFrame, tnt: pd.DataFrame, X: np.ndarray, n_splits: int,
                           rng: np.random.Generator) -> pd.DataFrame:
    """For each measure, the ρ over variants between its values against one random half of the
    recordings (each genotype split in two) and against the other half; mean and 5 / 95 % over
    splits. The weights of the weighted measures come from each half's own recordings."""
    C, T = ctrl.to_numpy(), tnt.to_numpy()
    res = {k: [] for k in RANKED if k != 'kendall'}
    for _ in range(n_splits):
        pc, pt = rng.permutation(len(C)), rng.permutation(len(T))
        halves = []
        for sc, st in ((pc[:len(C) // 2], pt[:len(T) // 2]), (pc[len(C) // 2:], pt[len(T) // 2:])):
            dh = C[sc].mean(0) - T[st].mean(0)
            se2 = C[sc].var(0, ddof=1) / len(sc) + T[st].var(0, ddof=1) / len(st)
            halves.append(metrics(M, dh, 1 / se2, X, kendall=False))
        for k in res:
            res[k].append(stats.spearmanr(halves[0][k], halves[1][k]).statistic)
    return pd.DataFrame([dict(metric=k, split_half_rho=np.mean(v), lo=np.percentile(v, 5), hi=np.percentile(v, 95),
        n_splits=n_splits) for k, v in res.items()])


# the models whose single-odor shares of the error were quoted (viki's reply of 2026-11-13 13:00)
ODOR_SHARE_MODELS = {'C2 (breadth-2 APL-effect fit, real wiring; E67)': ('E67_apl_effect_fits', 'rate0.068_aplboth_k2_cap3'),
                     'E86 C1_B_P3': ('E86_reduced_input_candidate', 'C1_B_P3')}


def odor_share_reliability(df: pd.DataFrame, M: np.ndarray, ctrl: pd.DataFrame, tnt: pd.DataFrame, X: np.ndarray,
                           n_splits: int, seed: int = SEED) -> pd.DataFrame:
    """How reproducible single odors' shares are between random halves of the recordings: per split, the r over
    the 16 odors between the shares (`X`'s least-squares projection) computed against each half, for the data's own
    effect, for the error of each of ODOR_SHARE_MODELS (each side's mean removed) and the r over pairs of its
    pair-specific part; mean and 5 / 95 % over splits."""
    rng = np.random.default_rng(seed)
    C, T = ctrl.to_numpy(), tnt.to_numpy()
    Xp = np.linalg.pinv(X)
    models = {}
    for lab, (folder, tag) in ODOR_SHARE_MODELS.items():
        i = df.index[(df.folder == folder) & (df.tag == tag)]
        if len(i):
            models[lab] = M[i[0]]
    res = {'data effect: odor shares': []}
    for lab in models:
        res[f'{lab}: error odor shares'] = []
        res[f'{lab}: error pair-specific part'] = []
    for _ in range(n_splits):
        pc, pt = rng.permutation(len(C)), rng.permutation(len(T))
        ds = [C[pc[:len(C) // 2]].mean(0) - T[pt[:len(T) // 2]].mean(0), C[pc[len(C) // 2:]].mean(0) - T[pt[len(T) // 2:]].mean(0)]
        res['data effect: odor shares'].append(np.corrcoef(*[Xp @ (dh - dh.mean()) for dh in ds])[0, 1])
        for lab, m in models.items():
            es = [(m - m.mean()) - (dh - dh.mean()) for dh in ds]
            sh = [Xp @ e for e in es]
            res[f'{lab}: error odor shares'].append(np.corrcoef(*sh)[0, 1])
            res[f'{lab}: error pair-specific part'].append(np.corrcoef(*[e - X @ s for e, s in zip(es, sh)])[0, 1])
    return pd.DataFrame([dict(what=k, r_between_halves=np.mean(v), lo=np.percentile(v, 5), hi=np.percentile(v, 95),
        n_splits=n_splits) for k, v in res.items()])


# the measures whose best variant per experiment is compared (best_pick_stability.csv)
PICK = ['loss', 'pattern_part', 'shape_part', 'spearman', 'pearson', 'pattern_part_w']
# experiment folders with at least this many distinct variants inside the bands enter split_half_subsets.csv's
# 'within experiments'
MIN_WITHIN = 8


def _best(tab: pd.DataFrame, k: str, g: np.ndarray) -> int:
    """The row of `g` that measure `k` of `tab` calls best (the first of ties)."""
    v = tab[k].to_numpy()[g]
    return int(g[np.argmin(v) if METRICS[k] else np.argmax(v)])


def split_half_extras(df: pd.DataFrame, M: np.ndarray, ctrl: pd.DataFrame, tnt: pd.DataFrame, X: np.ndarray,
                      n_splits: int, seed: int = SEED) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The comparisons of hal's read-only check (2026-11-13) on the splits of `split_half_reliability` (a generator
    seeded alike, so the same halves).

    Args:
        df: by_variant_metrics' rows (`folder`, `duplicate`, `realistic` and the measures against all recordings),
            row i the variant of M's row i.
        M: the variants' per-pair effects.

    Returns:
        subsets: per ranked measure and subset ('all distinct', 'inside the bands': the distinct variants inside the
            bands, 'within experiments': the mean over experiment folders with MIN_WITHIN or more of those of the ρ
            within the folder), the ρ over variants between the measure against one half and against the other (mean
            over splits, 5 / 95 %).
        picks: per measure of PICK, over the experiment folders with 2 or more variants inside the bands (duplicates
            kept, as best_per_experiment.csv): the number whose best variant differs from the pattern part's against
            all recordings and against one same half, and between the two halves; and of the distinct variants inside
            the bands, the rank by the second half of the variant the first half calls best (median over splits) and
            how often it is in the second half's best TOP_K.
        odors: per subset (the distinct variants inside the bands, and the TOP_K of them with the lowest pattern part
            against all recordings), the r over the 16 odors between the single-odor shares of a variant's error
            (each side's mean removed) against the two halves, median over the subset's variants (mean over splits,
            5 / 95 %).
    """
    rng = np.random.default_rng(seed)
    C, T = ctrl.to_numpy(), tnt.to_numpy()
    uniq, real = ~df.duplicate.to_numpy(bool), df.realistic.to_numpy(bool)
    ub = np.flatnonzero(uniq & real)
    folders = df.folder.to_numpy()
    groups = [np.flatnonzero(real & (folders == f)) for f in pd.unique(folders)]
    groups = [g for g in groups if len(g) >= 2]
    within = [g[uniq[g]] for g in groups]
    within = [g for g in within if len(g) >= MIN_WITHIN]
    lead = ub[np.argsort(df.pattern_part.to_numpy()[ub], kind='stable')[:TOP_K]]
    keys = [k for k in RANKED if k != 'kendall']
    subsets = {'all distinct': np.flatnonzero(uniq), 'inside the bands': ub}
    rel = {(s, k): [] for s in list(subsets) + ['within experiments'] for k in keys}
    full_best = {k: [_best(df, k, g) for g in groups] for k in PICK}
    differs, cross, rank, top = ({k: [] for k in PICK} for _ in range(4))
    Pinv = np.linalg.pinv(X)
    odor = {'inside the bands': [], f'the {TOP_K} with the lowest pattern part': []}
    Mc = M - M.mean(1, keepdims=True)
    for _ in range(n_splits):
        pc, pt = rng.permutation(len(C)), rng.permutation(len(T))
        halves, ds = [], []
        for sc, st in ((pc[:len(C) // 2], pt[:len(T) // 2]), (pc[len(C) // 2:], pt[len(T) // 2:])):
            dh = C[sc].mean(0) - T[st].mean(0)
            se2 = C[sc].var(0, ddof=1) / len(sc) + T[st].var(0, ddof=1) / len(st)
            halves.append(metrics(M, dh, 1 / se2, X, kendall=False))
            ds.append(dh)
        for k in keys:
            a, b = halves[0][k].to_numpy(), halves[1][k].to_numpy()
            for s, ix in subsets.items():
                rel[(s, k)].append(stats.spearmanr(a[ix], b[ix]).statistic)
            rel[('within experiments', k)].append(np.nanmean([stats.spearmanr(a[g], b[g]).statistic for g in within])
                if within else np.nan)
        pat0 = [_best(halves[0], 'pattern_part', g) for g in groups]
        for k in PICK:
            b0 = [_best(halves[0], k, g) for g in groups]
            b1 = [_best(halves[1], k, g) for g in groups]
            differs[k].append(sum(x != y for x, y in zip(b0, b1)))
            cross[k].append(sum(x != y for x, y in zip(b0, pat0)))
            sgn = 1 if METRICS[k] else -1
            a, b = sgn * halves[0][k].to_numpy()[ub], sgn * halves[1][k].to_numpy()[ub]
            r = stats.rankdata(b)[np.argmin(a)]
            rank[k].append(r)
            top[k].append(r <= TOP_K)
        rr = _rowcorr((Mc - (ds[0] - ds[0].mean())) @ Pinv.T, (Mc - (ds[1] - ds[1].mean())) @ Pinv.T)
        odor['inside the bands'].append(np.median(rr[ub]))
        odor[f'the {TOP_K} with the lowest pattern part'].append(np.median(rr[lead]))
    pct = lambda v: (float(np.nanpercentile(v, 5)), float(np.nanpercentile(v, 95)))
    n_of = {'all distinct': len(subsets['all distinct']), 'inside the bands': len(ub), 'within experiments': len(within)}
    sub = pd.DataFrame([dict(metric=k, subset=s, split_half_rho=float(np.nanmean(v)), lo=pct(v)[0], hi=pct(v)[1],
        n=n_of[s], n_splits=n_splits) for (s, k), v in rel.items()])
    picks = pd.DataFrame([dict(metric=k, n_experiments=len(groups),
        differs_from_pattern_part_all_recordings=int(sum(x != y for x, y in zip(full_best[k], full_best['pattern_part']))),
        differs_from_pattern_part_one_half=float(np.mean(cross[k])),
        differs_between_halves=float(np.mean(differs[k])), differs_between_halves_lo=pct(differs[k])[0],
        differs_between_halves_hi=pct(differs[k])[1], n_inside_bands=len(ub),
        best_rank_by_other_half_median=float(np.median(rank[k])), best_in_other_top25_share=float(np.mean(top[k])),
        n_splits=n_splits) for k in PICK])
    odors = pd.DataFrame([dict(subset=s, n_variants=len(ub) if s == 'inside the bands' else len(lead),
        median_r_between_halves=float(np.mean(v)), lo=pct(v)[0], hi=pct(v)[1], n_splits=n_splits) for s, v in odor.items()])
    return sub, picks, odors


def archived_counts(folder: str) -> Dict[str, Path]:
    """{tag: directory} of the variants of an experiment's model-output archive folder with both APL-on and acute
    APL-off spike counts (`<tag>_on_spike_counts.parquet`, `<tag>_acute_off_spike_counts.parquet`): the folder itself
    first, then its subfolders in sorted order."""
    root = AE.ARCHIVE / folder
    out: Dict[str, Path] = {}
    if not root.is_dir():
        return out
    for q in sorted(root.rglob('*_on_spike_counts.parquet'), key=lambda q: (q.parent != root, str(q.parent), q.name)):
        tag = q.name[:-len('_on_spike_counts.parquet')]
        if tag not in out and (q.parent / f'{tag}_acute_off_spike_counts.parquet').exists():
            out[tag] = q.parent
    return out


def _readout_effect(job: tuple) -> dict:
    """One archived variant under the calibrated readout (a worker of `readout_effects`)."""
    folder, tag, where, target, index = job
    try:
        on = pd.read_parquet(Path(where) / f'{tag}_on_spike_counts.parquet')
        off = pd.read_parquet(Path(where) / f'{tag}_acute_off_spike_counts.parquet')
        K = H.calibrate_readout(on, 'saturating', target, stat='participation_ratio', responders_min=2.0)
        at_bound = abs(H.readout_stat(on, 'saturating', K, 'participation_ratio', 2.0) - target) > 0.02 * target
        cal = H.apl_effect(H.readout_pair_distances(on, 'saturating', K), H.readout_pair_distances(off, 'saturating', K))
        raw = H.apl_effect(H.readout_pair_distances(on, 'raw'), H.readout_pair_distances(off, 'raw'))
        return dict(folder=folder, tag=tag, K=float(K), K_at_bound=bool(at_bound), where=str(where),
                    cal=cal.reindex(index).to_numpy(dtype=float), raw=raw.reindex(index).to_numpy(dtype=float))
    except Exception as ex:   # e.g. a variant of other odors
        return dict(folder=folder, tag=tag, error=repr(ex)[:200])


def readout_effects(index: pd.Index, target: float, jobs: int = 4) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, dict]:
    """Every variant of results/apl_effect/by_variant.csv with archived spike counts, read through the calibrated
    readout (`--readout cal`, see the module docstring).

    Returns:
        (meta: exp, folder, tag, K, K_at_bound, archive_dir; the per-pair effects under the readout and on raw counts,
        rows as meta, columns `index`; counts of the variants scored, archived and left out, and why).
    """
    bv = pd.read_csv(AE.OUT / 'by_variant.csv', dtype={'tag': str})
    work = []
    for folder, g in bv.groupby('folder', sort=False):
        where = archived_counts(folder)
        work += [(folder, t, str(where[t]), target, index) for t in g.tag if t in where]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        res = list(ex.map(_readout_effect, work, chunksize=8))
    rows, cal, raw, errors, partial = [], [], [], [], 0
    for r in res:
        if 'error' in r:
            errors.append(f"{r['folder']}/{r['tag']}: {r['error']}")
            continue
        if np.isnan(r['cal']).any() or np.isnan(r['raw']).any():
            partial += 1
            continue
        m = re.match(r'(E\d+)_', r['folder'])
        rows.append(dict(exp=m.group(1) if m else '', folder=r['folder'], tag=r['tag'], K=r['K'], K_at_bound=r['K_at_bound'],
                         archive_dir=str(Path(r['where']).relative_to(AE.ARCHIVE))))
        cal.append(r['cal'])
        raw.append(r['raw'])
    counts = dict(n_scored=len(bv), n_archived=len(work), n_unreadable=len(errors), n_not_all_pairs=partial, errors=errors[:20])
    return pd.DataFrame(rows), np.array(cal), np.array(raw), counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--n-splits', type=int, default=200)
    ap.add_argument('--readout', choices=('raw', 'cal'), default='raw', help="'raw': the per-pair effects of the "
        "spike counts' cosine distances (the saved distances); 'cal': the calibrated readout's, from the model-output "
        "archive's spike counts (tables named *_cal)")
    ap.add_argument('--jobs', type=int, default=4, help='worker processes reading the archive (--readout cal)')
    ap.add_argument('--out', default=str(OUT), help='output folder (default the shared results/apl_effect_metrics/)')
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sfx = '' if a.readout == 'raw' else f'_{a.readout}'
    rng = np.random.default_rng(SEED)

    ctrl, tnt, dser = data_recordings()
    d = dser.to_numpy(dtype=float)
    se = pair_se(ctrl, tnt)
    X, odors = H.pair_design(dser.index, H.APL16_ODORS)
    by_variant = pd.read_csv(AE.OUT / 'by_variant.csv', dtype={'tag': str})
    extra_prov = {}
    if a.readout == 'raw':
        pd.DataFrame(dict(odor_a=[p[0] for p in dser.index], odor_b=[p[1] for p in dser.index], effect=d, se=se,
            weight=(1 / se ** 2) / (1 / se ** 2).sum())).to_csv(out / 'data_pairs.csv', index=False)
        meta, M = model_effects(dser.index)
    else:
        os.environ.setdefault('OMP_NUM_THREADS', '1')
        target = H.data_target('ctrl_participation_ratio_resp')
        meta, M, Mraw, counts = readout_effects(dser.index, target, a.jobs)
        # the archived counts must give the saved distances (raw counts: the loss and pattern part of by_variant.csv)
        raw = metrics(Mraw, d, 1 / se ** 2, X, kendall=False)
        chk = meta[['folder', 'tag']].assign(loss=raw.loss.to_numpy(), pattern=raw.pattern_part.to_numpy()).merge(
            by_variant[['folder', 'tag', f'{AE.V}__apl_effect_rmse', f'{AE.V}__apl_effect_rmse_centered']], on=['folder', 'tag'],
            how='left')
        dev = np.maximum((chk.loss - chk[f'{AE.V}__apl_effect_rmse']).abs(),
                         (chk.pattern - chk[f'{AE.V}__apl_effect_rmse_centered']).abs()).to_numpy()
        keep = dev <= 1e-9
        meta, M = meta[keep].reset_index(drop=True), M[keep]
        extra_prov = dict(readout='saturating, K calibrated per variant on its APL-on counts to the participation ratio of '
            'the responders (two or more spikes), registry ctrl_participation_ratio_resp', cal_target=target, **counts,
            n_raw_mismatch=int((~keep).sum()), raw_mismatch=[f'{f}/{t}' for f, t in zip(chk.folder[~keep], chk.tag[~keep])][:20],
            raw_check_max_abs_diff=float(np.nanmax(dev[keep])) if keep.any() else None,
            n_K_at_bound=int(meta.K_at_bound.sum()))
        print(f"{counts['n_scored']} scored variants, {counts['n_archived']} with archived spike counts; left out: "
              f"{counts['n_unreadable']} unreadable, {counts['n_not_all_pairs']} without all 120 pairs, {int((~keep).sum())} "
              'whose archived counts do not give the saved distances', flush=True)
    _, first = np.unique(np.round(M, 12), axis=0, return_index=True)
    meta['duplicate'] = ~np.isin(np.arange(len(M)), first)
    df = pd.concat([meta, metrics(M, d, 1 / se ** 2, X)], axis=1)
    bv = by_variant[['folder', 'tag', 'realism_flag', 'rate_on', 'rate_off']]
    df = df.merge(bv, on=['folder', 'tag'], how='left')
    df['realistic'] = df.realism_flag.isna() | (df.realism_flag.astype(str) == '')
    if a.readout == 'raw':
        # the same loss as kcq/apl_effect.py (a check that the variants are read alike)
        chk = by_variant.merge(df[['folder', 'tag', 'loss', 'spearman']], on=['folder', 'tag'])
        if len(chk) and not np.allclose(chk['loss'], chk[f'{AE.V}__apl_effect_rmse'], atol=1e-9):
            raise RuntimeError('the loss differs from results/apl_effect/by_variant.csv')
    df.to_csv(out / f'by_variant_metrics{sfx}.csv', index=False)

    uniq = ~df.duplicate.to_numpy()
    o = _oriented(df)
    o[uniq].corr(method='spearman').to_csv(out / f'metric_correlations{sfx}.csv')
    o[uniq & df.realistic.to_numpy()].corr(method='spearman').to_csv(out / f'metric_correlations_realistic{sfx}.csv')
    g = df[df.realistic]
    pd.DataFrame(dict(folder=g.folder, tag=g.tag, mean_share=g.mean_part ** 2 / g.loss ** 2,
        amplitude_share=g.amplitude_part / g.pattern_part ** 2, per_odor_share=g.per_odor_part ** 2 / g.pattern_part ** 2,
        sd_ratio=g.sd_ratio)).to_csv(out / f'shares{sfx}.csv', index=False)

    rel = split_half_reliability(M[uniq], ctrl, tnt, X, a.n_splits, rng)
    rel.to_csv(out / f'split_half_reliability{sfx}.csv', index=False)
    odor_share_reliability(df, M, ctrl, tnt, X, a.n_splits).to_csv(out / f'odor_share_reliability{sfx}.csv', index=False)
    sub, picks, odor_agree = split_half_extras(df, M, ctrl, tnt, X, a.n_splits)
    sub.to_csv(out / f'split_half_subsets{sfx}.csv', index=False)
    picks.to_csv(out / f'best_pick_stability{sfx}.csv', index=False)
    odor_agree.to_csv(out / f'odor_share_agreement{sfx}.csv', index=False)

    best_rows, top = [], {}
    for k in RANKED:
        for subset, g in (('all', df), ('realistic', df[df.realistic])):
            i = g[k].idxmin() if METRICS[k] else g[k].idxmax()
            best_rows.append(dict(metric=k, subset=subset, exp=g.at[i, 'exp'], folder=g.at[i, 'folder'], tag=g.at[i, 'tag'],
                value=g.at[i, k], loss=g.at[i, 'loss'], pattern_part=g.at[i, 'pattern_part'], spearman=g.at[i, 'spearman']))
        g = df[df.realistic & ~df.duplicate]
        top[k] = set((g.nsmallest(TOP_K, k) if METRICS[k] else g.nlargest(TOP_K, k)).index)
    pd.DataFrame(best_rows).to_csv(out / f'best_by_metric{sfx}.csv', index=False)
    pd.DataFrame({k: {j: len(top[k] & top[j]) for j in RANKED} for k in RANKED}).to_csv(out / f'top25_overlap{sfx}.csv')

    per_exp = []
    for (exp, folder), g in df[df.realistic].groupby(['exp', 'folder'], sort=False):
        if len(g) < 2:
            continue
        row = dict(exp=exp, folder=folder, n_variants=len(g))
        for k in ('loss', 'pattern_part', 'spearman', 'pearson', 'pattern_part_w', 'shape_part'):
            row[f'best_{k}'] = g.at[g[k].idxmin() if METRICS[k] else g[k].idxmax(), 'tag']
        per_exp.append(row)
    pe = pd.DataFrame(per_exp)
    pe.to_csv(out / f'best_per_experiment{sfx}.csv', index=False)

    (out / f'provenance{sfx}.json').write_text(json.dumps(dict(script='analyses/apl_effect_metrics.py',
        run_at=datetime.now(timezone.utc).isoformat(), n_variants=len(df), n_realistic=int(df.realistic.sum()),
        n_unique=int(uniq.sum()), n_unique_realistic=int((uniq & df.realistic.to_numpy()).sum()),
        n_splits=a.n_splits, seed=SEED, top_k=TOP_K, **extra_prov, **H.code_commits()), indent=1, default=str) + '\n')
    print(f'{len(df)} variants ({int(df.realistic.sum())} inside the bands; {int(uniq.sum())} distinct); outputs in {out}')


if __name__ == '__main__':
    main()
