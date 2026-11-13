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
shares of the pattern part squared, and the ratio of the spreads), data_pairs.csv (the data's
per-pair effect and its standard error) and provenance.json.

Usage: ./venv/bin/python analyses/apl_effect_metrics.py [--n-splits 200]
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from typing import Dict, Tuple

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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--n-splits', type=int, default=200)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    ctrl, tnt, dser = data_recordings()
    d = dser.to_numpy(dtype=float)
    se = pair_se(ctrl, tnt)
    X, odors = H.pair_design(dser.index, H.APL16_ODORS)
    pd.DataFrame(dict(odor_a=[p[0] for p in dser.index], odor_b=[p[1] for p in dser.index], effect=d, se=se,
        weight=(1 / se ** 2) / (1 / se ** 2).sum())).to_csv(OUT / 'data_pairs.csv', index=False)

    meta, M = model_effects(dser.index)
    _, first = np.unique(np.round(M, 12), axis=0, return_index=True)
    meta['duplicate'] = ~np.isin(np.arange(len(M)), first)
    df = pd.concat([meta, metrics(M, d, 1 / se ** 2, X)], axis=1)
    bv = pd.read_csv(AE.OUT / 'by_variant.csv', dtype={'tag': str})[['folder', 'tag', 'realism_flag', 'rate_on', 'rate_off']]
    df = df.merge(bv, on=['folder', 'tag'], how='left')
    df['realistic'] = df.realism_flag.isna() | (df.realism_flag.astype(str) == '')
    # the same loss as kcq/apl_effect.py (a check that the variants are read alike)
    chk = pd.read_csv(AE.OUT / 'by_variant.csv', dtype={'tag': str}).merge(df[['folder', 'tag', 'loss', 'spearman']],
        on=['folder', 'tag'])
    if len(chk) and not np.allclose(chk['loss'], chk[f'{AE.V}__apl_effect_rmse'], atol=1e-9):
        raise RuntimeError('the loss differs from results/apl_effect/by_variant.csv')
    df.to_csv(OUT / 'by_variant_metrics.csv', index=False)

    uniq = ~df.duplicate.to_numpy()
    o = _oriented(df)
    o[uniq].corr(method='spearman').to_csv(OUT / 'metric_correlations.csv')
    o[uniq & df.realistic.to_numpy()].corr(method='spearman').to_csv(OUT / 'metric_correlations_realistic.csv')
    g = df[df.realistic]
    pd.DataFrame(dict(folder=g.folder, tag=g.tag, mean_share=g.mean_part ** 2 / g.loss ** 2,
        amplitude_share=g.amplitude_part / g.pattern_part ** 2, per_odor_share=g.per_odor_part ** 2 / g.pattern_part ** 2,
        sd_ratio=g.sd_ratio)).to_csv(OUT / 'shares.csv', index=False)

    rel = split_half_reliability(M[uniq], ctrl, tnt, X, a.n_splits, rng)
    rel.to_csv(OUT / 'split_half_reliability.csv', index=False)

    best_rows, top = [], {}
    for k in RANKED:
        for subset, g in (('all', df), ('realistic', df[df.realistic])):
            i = g[k].idxmin() if METRICS[k] else g[k].idxmax()
            best_rows.append(dict(metric=k, subset=subset, exp=g.at[i, 'exp'], folder=g.at[i, 'folder'], tag=g.at[i, 'tag'],
                value=g.at[i, k], loss=g.at[i, 'loss'], pattern_part=g.at[i, 'pattern_part'], spearman=g.at[i, 'spearman']))
        g = df[df.realistic & ~df.duplicate]
        top[k] = set((g.nsmallest(TOP_K, k) if METRICS[k] else g.nlargest(TOP_K, k)).index)
    pd.DataFrame(best_rows).to_csv(OUT / 'best_by_metric.csv', index=False)
    pd.DataFrame({k: {j: len(top[k] & top[j]) for j in RANKED} for k in RANKED}).to_csv(OUT / 'top25_overlap.csv')

    per_exp = []
    for (exp, folder), g in df[df.realistic].groupby(['exp', 'folder'], sort=False):
        if len(g) < 2:
            continue
        row = dict(exp=exp, folder=folder, n_variants=len(g))
        for k in ('loss', 'pattern_part', 'spearman', 'pearson', 'pattern_part_w', 'shape_part'):
            row[f'best_{k}'] = g.at[g[k].idxmin() if METRICS[k] else g[k].idxmax(), 'tag']
        per_exp.append(row)
    pe = pd.DataFrame(per_exp)
    pe.to_csv(OUT / 'best_per_experiment.csv', index=False)

    (OUT / 'provenance.json').write_text(json.dumps(dict(script='analyses/apl_effect_metrics.py',
        run_at=datetime.now(timezone.utc).isoformat(), n_variants=len(df), n_realistic=int(df.realistic.sum()),
        n_unique=int(uniq.sum()), n_unique_realistic=int((uniq & df.realistic.to_numpy()).sum()),
        n_splits=a.n_splits, seed=SEED, top_k=TOP_K, **H.code_commits()), indent=1, default=str) + '\n')
    print(f'{len(df)} variants ({int(df.realistic.sum())} inside the bands; {int(uniq.sum())} distinct); outputs in {OUT}')


if __name__ == '__main__':
    main()
