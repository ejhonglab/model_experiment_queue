#!/usr/bin/env python3
"""How much the 16-odor ORN reference matrix's mixed recipe matters: per-fly correlations for the pairs
within a panel, fly-mean vectors for the pairs across the two panels.

`kcq.harness.orn16_consensus_pairs()` (data_prep/orn16_distances_from_olfdb.py) gives each of the 72 pairs within
a panel (66 megamat x megamat, 6 validation2 x validation2) as the mean over flies of each fly's 1 - Pearson r
across the 36 glomeruli, and each of the 48 pairs across the panels as 1 - r between the two odors' fly-mean
vectors, since no fly was imaged with both panels. A per-fly correlation carries that fly's measurement noise, so
the two kinds of pairs sit on different scales. Yang asked (2026-11-13) why the matrix is not built one way for
all 120 pairs, from fly-mean glomerulus vectors stitched across the panels (the column
`distance_flymean_vectors` of orn16_corrdist_olfdb_36glom_pairs.csv, which the same script writes).

This compares the two recipes on everything the ORN matrix enters:
- the matrix itself: the per-fly minus fly-mean difference of the 72 within-panel pairs, the mean distance of
  the within-panel and cross-panel pairs under each recipe, and the registry's ORN entries (mean, min, max,
  pairs more than 0.8 apart, within-group and fruity x small-ester means);
- the KC data against the ORNs: ρ on odor-odor distances of the control flies' and of the APL output blocked
  (TNT) matrix with the ORNs (registry rho_ctrl_orn, rho_tnt_orn) and `compare_on_off`'s ORN statistics
  (spearman_orn_on / _off, gap_closed, apl_share_resid_var, spearman_resid_orn_on_off);
- every saved model variant's ρ on odor-odor distances with the ORNs, APL on and APL off (the stats tables'
  spearman_orn_on / _off of the model side), from each results folder's model_pair_distances_long.csv (spike-count
  cosine rows, as kcq/apl_effect.py reads them).

Writes <shared>/results/orn_recipe_check/: matrix_summary.csv, data_stats.csv, model_rho_orn.csv (per variant,
both recipes), model_rho_orn_summary.csv, pairs.csv and provenance.json.
Usage: ./venv/bin/python analyses/orn_recipe_check.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import apl_effect as AE

INPUTS = H.shared_dir() / 'data' / 'model_outputs' / 'inputs'
OUT = H.shared_dir() / 'results' / 'orn_recipe_check'


def recipes() -> tuple:
    """(pairs table, {recipe: pair Series}): 'mixed' (the reference until 2026-11-13,
    `H.orn16_consensus_pairs_mixed()`) and 'flymean' (fly-mean vectors for all 120 pairs, the reference since,
    `H.orn16_consensus_pairs()`), each checked against the harness. The pairs table's columns are
    distance_mixed and distance_flymean_vectors since 2026-11-13 (before, `distance` was the mixed recipe)."""
    t = pd.read_csv(INPUTS / 'orn16_corrdist_olfdb_36glom_pairs.csv')
    idx = pd.MultiIndex.from_arrays([t.odor_a, t.odor_b], names=['odor_a', 'odor_b'])
    mixed_col = 'distance_mixed' if 'distance_mixed' in t else 'distance'
    mixed = H.canonical_pairs(pd.Series(t[mixed_col].to_numpy(), index=idx))
    flymean = H.canonical_pairs(pd.Series(t.distance_flymean_vectors.to_numpy(), index=idx))
    for name, cur, rec in (('orn16_consensus_pairs_mixed', H.orn16_consensus_pairs_mixed(), mixed),
                           ('orn16_consensus_pairs', H.orn16_consensus_pairs(), flymean)):
        a, b = H.align_pairs(cur, rec)
        if len(a) != 120 or not np.allclose(a.to_numpy(), b.to_numpy(), atol=1e-12):
            raise RuntimeError(f'the pairs table does not give harness.{name}()')
    t['within_panel'] = t.panel_a == t.panel_b
    return t, {'mixed': mixed, 'flymean': flymean}


def matrix_summary(t: pd.DataFrame, rec: dict) -> pd.DataFrame:
    rows = []
    kc = H.yang_recording_rdms('TNT3c', 'cosine').mean()
    for name, m in rec.items():
        m, _ = H.align_pairs(m, kc)
        idx = m.index
        within = np.array([H.ODOR_GROUPS.get(x, x) == H.ODOR_GROUPS.get(y, y) for x, y in idx])
        block = np.array([H.block_label(x, y) == H.FRUITY_SMALL_BLOCK for x, y in idx])
        panel = {(r.odor_a, r.odor_b): r.within_panel for r in t.itertuples()}
        wp = np.array([panel.get((x, y), panel.get((y, x))) for x, y in idx], dtype=bool)
        rows.append(dict(recipe=name, mean=m.mean(), min=m.min(), max=m.max(), n_gt08=int((m > 0.8).sum()),
            within_group_mean=m[within].mean(), fruity_small_mean=m[block].mean(),
            within_panel_mean=m[wp].mean(), cross_panel_mean=m[~wp].mean()))
    return pd.DataFrame(rows)


def data_stats(rec: dict) -> pd.DataFrame:
    kc = {g: H.yang_recording_rdms(g, 'cosine').mean() for g in ('TNT3c', 'TNT')}
    rows = []
    for name, orn in rec.items():
        o, c = H.align_pairs(orn, kc['TNT3c'])
        _, t = H.align_pairs(c, kc['TNT'])
        st, _ = H.compare_on_off(c, t, orn=o)
        rows.append(dict(recipe=name, rho_ctrl_orn=H.spearman(c, o), rho_tnt_orn=H.spearman(t, o),
            **{k: st[k] for k in ('spearman_orn_on', 'spearman_orn_off', 'gap_closed', 'apl_share_resid_var',
                                  'spearman_resid_orn_on_off') if k in st}))
    return pd.DataFrame(rows)


def model_rho(rec: dict) -> pd.DataFrame:
    rows = []
    for folder in sorted(AE.RESULTS.glob('E*_*')):
        if not folder.is_dir():
            continue
        for tag, (on, off) in AE.model_pairs(folder).items():
            row = dict(folder=folder.name, tag=tag)
            for name, orn in rec.items():
                for cond, d in (('on', on), ('off', off)):
                    a, b = H.align_pairs(orn, d)
                    row[f'rho_{cond}_{name}'] = H.spearman(b, a) if len(a) >= 100 else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    for cond in ('on', 'off'):
        df[f'diff_{cond}'] = df[f'rho_{cond}_flymean'] - df[f'rho_{cond}_mixed']
    return df


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t, rec = recipes()
    wp = t[t.within_panel]
    d = (wp.distance_mixed if 'distance_mixed' in wp else wp.distance) - wp.distance_flymean_vectors
    t.to_csv(OUT / 'pairs.csv', index=False)
    ms = matrix_summary(t, rec)
    ms.to_csv(OUT / 'matrix_summary.csv', index=False)
    ds = data_stats(rec)
    ds.to_csv(OUT / 'data_stats.csv', index=False)
    mr = model_rho(rec)
    mr.to_csv(OUT / 'model_rho_orn.csv', index=False)
    summ = []
    for cond in ('on', 'off'):
        x = mr[f'diff_{cond}'].dropna()
        summ.append(dict(condition=cond, n_variants=len(x), mean_change=x.mean(), median_abs_change=x.abs().median(),
            p95_abs_change=x.abs().quantile(0.95), max_abs_change=x.abs().max(),
            rank_corr_over_variants=mr[[f'rho_{cond}_mixed', f'rho_{cond}_flymean']].corr('spearman').iloc[0, 1]))
    pd.DataFrame(summ).to_csv(OUT / 'model_rho_orn_summary.csv', index=False)
    (OUT / 'provenance.json').write_text(json.dumps(dict(script='analyses/orn_recipe_check.py',
        run_at=datetime.now(timezone.utc).isoformat(), within_panel_perfly_minus_flymean=dict(mean=float(d.mean()),
            min=float(d.min()), max=float(d.max()), n=int(len(d))),
        inputs={str(INPUTS / 'orn16_corrdist_olfdb_36glom_pairs.csv'): AE._sha(INPUTS / 'orn16_corrdist_olfdb_36glom_pairs.csv')},
        **H.code_commits()), indent=1, default=str) + '\n')
    pd.set_option('display.width', 200)
    print(f'within-panel pairs, per-fly minus fly-mean: mean {d.mean():+.3f}, min {d.min():+.3f}, max {d.max():+.3f} (n {len(d)})')
    print(ms.round(3).to_string(index=False))
    print(ds.round(3).to_string(index=False))
    print(pd.DataFrame(summ).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
