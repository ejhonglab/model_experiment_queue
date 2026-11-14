#!/usr/bin/env python3
"""Rebuild the 16-odor ORN odor-odor distance matrix (1 - Pearson r across glomeruli) from
olfdb, over the 36 glomeruli the model's ORN input uses.

The matrix the program has used as its 16-odor ORN reference (the shared directory's apl_pair_analysis/
orn16_consensus_corrdist.csv, now `kcq.harness.orn16_consensus_pairs_zero_filled()`) was computed by
the original report's analysis from al_analysis/data/sent_to_vlad/2024-06-21/consensus_mean.csv
over the union of the megamat and validation2 panels' glomeruli, with 0 in place of glomeruli a
panel did not image (DL4 and VM4 for the megamat odors; VA4 and VL1 for the validation2 odors),
so unmeasured glomeruli counted as zero responses (unresolved_issues.md). This script recomputes
it from Tom's ORN data in olfdb (dataset 1, megamat panel, 9 flies; dataset 4, validation2 panel,
6 flies; http://hal:8080) over the 36 glomeruli both panels imaged, the same set as the model's
input (`kcq.harness.apl16_orn_deltas()`):

- every pair (since 2026-11-13, Tom's decision on Yang's question): 1 - r between the two odors' fly-mean
  glomerulus vectors (trial means per fly, then the mean over the flies of the odor's panel), the same
  recipe for the 72 pairs within a panel and the 48 across the panels, and the same vectors the model's
  input is stitched from.
Until 2026-11-13 the reference used a mixed recipe, kept for comparison: the pairs within a panel (66
megamat x megamat, 6 validation2 x validation2) as the mean over flies of each fly's 1 - r, with
al_analysis's `al_util.mean_of_fly_corrs` (trial means, one correlation per fly across its ROIs, averaged
over flies ignoring NaN), and the pairs across panels, where no fly was imaged with both panels, from the
fly-mean vectors. A per-fly correlation carries that fly's noise, so those 72 pairs sat about 0.03 further
apart than the 48 across the panels (unresolved_issues.md, Resolved).

Outputs (shared directory, data/model_outputs/inputs/; the existing ORN input
orn_deltas_apl16*.csv files there are not touched):
  orn16_corrdist_olfdb_36glom_flymean.csv  16 x 16, the matrix described above (diagonal 0); the reference,
                                           `kcq.harness.orn16_consensus_pairs()`
  orn16_corrdist_olfdb_36glom.csv          16 x 16, the mixed recipe used until 2026-11-13,
                                           `kcq.harness.orn16_consensus_pairs_mixed()`; rewritten only if
                                           its values change, so that a process reading it keeps them
  orn16_corrdist_olfdb_36glom_pairs.csv    120 pairs: the reference, the mixed recipe, the method each used,
                                           and the zero-filled value of the original report
  orn16_corrdist_olfdb_36glom.json       provenance (olfdb dataset ids and sha256, glomeruli, odors,
                                         this repository's and al_analysis's commits)
Usage: ./venv/bin/python data_prep/orn16_distances_from_olfdb.py
"""
from __future__ import annotations

import io
import json
import subprocess
from datetime import datetime, timezone
from itertools import combinations
from urllib.request import urlopen

import numpy as np
import pandas as pd

from kcq import harness as H

OLFDB = 'http://hal:8080'
DATASETS = {'megamat': 1, 'validation2': 4}
OUT = H.shared_dir() / 'data' / 'model_outputs' / 'inputs'
STEM = 'orn16_corrdist_olfdb_36glom'


def fetch(ds_id: int) -> tuple:
    with urlopen(f'{OLFDB}/api/download?ids={ds_id}&format=long') as r:
        df = pd.read_csv(io.BytesIO(r.read()))
    with urlopen(f'{OLFDB}/api/datasets/{ds_id}') as r:
        meta = json.loads(r.read())
    return df, meta


def panel_trials(df: pd.DataFrame, panel: str, odors, gloms) -> pd.DataFrame:
    """Trial-level responses of `odors` (abbreviations, at the model input's concentrations):
    rows (odor1, repeat), columns (date, fly_num, roi) restricted to `gloms`."""
    want = {f'{o} @ {H.ORN_INPUT_CONC[o]:g}' for o in odors}
    d = df[(df.panel == panel) & df.odor1.isin(want) & df.roi.isin(gloms)]
    missing = want - set(d.odor1)
    assert not missing, (panel, missing)
    return d.pivot_table(index=['odor1', 'repeat'], columns=['date', 'fly_num', 'roi'], values='value')


def main() -> None:
    from al_analysis.al_util import mean_of_fly_corrs
    gloms = list(H.apl16_orn_deltas().index)
    assert len(gloms) == 36
    trials, metas = {}, {}
    panel_odors = {'megamat': list(H.MEGAMAT12), 'validation2': list(H.VALIDATION2_4)}
    for panel, ds in DATASETS.items():
        df, metas[panel] = fetch(ds)
        absent = set(gloms) - set(df[df.panel == panel].roi)
        assert not absent, (panel, absent)
        trials[panel] = panel_trials(df, panel, panel_odors[panel], gloms)
    abbrev = lambda s: s.split(' @ ')[0]

    # within-panel pairs: mean over flies of per-fly correlations (al_analysis's function)
    within = {}
    n_flies = {}
    for panel, t in trials.items():
        corr = mean_of_fly_corrs(t, square=True)
        corr.index = [abbrev(o) for o in corr.index]; corr.columns = [abbrev(o) for o in corr.columns]
        n_flies[panel] = t.columns.droplevel('roi').nunique()
        for a, b in combinations(panel_odors[panel], 2):
            within[(a, b)] = 1.0 - float(corr.loc[a, b])

    # fly-mean glomerulus vectors (trial mean per fly, then mean over flies), all 16 odors
    vec = {}
    for panel, t in trials.items():
        fm = t.groupby(level='odor1').mean().T.groupby(level='roi').mean()  # roi x odor
        for c in fm.columns:
            vec[abbrev(c)] = fm[c].reindex(gloms)
    V = pd.DataFrame(vec)[list(H.APL16_ODORS)]
    assert not V.isna().any().any()

    old = H.orn16_consensus_pairs_zero_filled()
    get = lambda s, a, b: s.get((a, b), s.get((b, a), np.nan))
    rows = []
    for a, b in combinations(H.APL16_ODORS, 2):
        pa = 'megamat' if a in H.MEGAMAT12 else 'validation2'
        pb = 'megamat' if b in H.MEGAMAT12 else 'validation2'
        flymean = 1.0 - float(np.corrcoef(V[a], V[b])[0, 1])
        if pa == pb:
            mixed, mixed_method = within[(a, b)], 'mean over flies of per-fly 1 - r (al_util.mean_of_fly_corrs)'
        else:
            mixed, mixed_method = flymean, '1 - r of fly-mean vectors (no fly imaged with both panels)'
        rows.append(dict(odor_a=a, odor_b=b, panel_a=pa, panel_b=pb, distance=flymean,
            method='1 - r of fly-mean vectors', distance_flymean_vectors=flymean,
            distance_mixed=mixed, method_mixed=mixed_method, distance_old_zero_filled=get(old, a, b)))
    pairs = pd.DataFrame(rows)

    def square(col):
        sq = pd.DataFrame(0.0, index=H.APL16_ODORS, columns=H.APL16_ODORS)
        for r in pairs.itertuples():
            sq.loc[r.odor_a, r.odor_b] = sq.loc[r.odor_b, r.odor_a] = getattr(r, col)
        return sq
    OUT.mkdir(parents=True, exist_ok=True)
    square('distance').to_csv(OUT / f'{STEM}_flymean.csv')
    mixed_sq, mixed_path = square('distance_mixed'), OUT / f'{STEM}.csv'
    if not (mixed_path.exists() and np.allclose(pd.read_csv(mixed_path, index_col=0).loc[H.APL16_ODORS, H.APL16_ODORS].to_numpy(),
                                                mixed_sq.to_numpy(), rtol=0, atol=1e-12)):
        mixed_sq.to_csv(mixed_path)
    pairs.to_csv(OUT / f'{STEM}_pairs.csv', index=False)
    prov = dict(generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), script='data_prep/orn16_distances_from_olfdb.py',
        **H.code_commits(),
        olfdb=OLFDB, datasets={p: dict(id=DATASETS[p], title=metas[p].get('title'), data_sha256=metas[p].get('data_sha256'), n_flies=int(n_flies[p]))
            for p in DATASETS},
        glomeruli=gloms, odors={o: f'{o} @ {H.ORN_INPUT_CONC[o]:g}' for o in H.APL16_ODORS},
        reference=f'{STEM}_flymean.csv: 1 - Pearson r between fly-mean glomerulus vectors (36 glomeruli; trial means per fly, then the mean over the panel\'s flies), all 120 pairs (since 2026-11-13)',
        mixed=f'{STEM}.csv (until 2026-11-13): within a panel al_analysis.al_util.mean_of_fly_corrs (trial means; per-fly Pearson r across that fly\'s ROIs among the 36; mean over flies ignoring NaN), distance = 1 - mean r; across panels the fly-mean vectors',
        approx_conc_note=H.APPROX_CONC_NOTE)
    (OUT / f'{STEM}.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')

    # comparison with the old matrix and the KC data
    pd.set_option('display.width', 200)
    pairs['kind'] = np.where(pairs.panel_a == pairs.panel_b, pairs.panel_a + ' x ' + pairs.panel_b, 'megamat x validation2')
    s = pairs.assign(d_ref_zero_filled=(pairs.distance - pairs.distance_old_zero_filled).abs(),
        d_mixed_ref=(pairs.distance_mixed - pairs.distance))
    print(s.groupby('kind')[['d_ref_zero_filled', 'd_mixed_ref']].agg(['mean', 'max']).round(4).to_string())
    ref = pairs.set_index(['odor_a', 'odor_b']).distance
    mixed = pairs.set_index(['odor_a', 'odor_b']).distance_mixed
    for g, lab in (('TNT3c', 'KCs, control flies'), ('TNT', 'KCs, APL output blocked (TNT)')):
        kc = H.yang_mean_rdm(g, 'cosine')
        out = []
        for name, o in (('zero-filled (original report)', old), ('mixed (until 2026-11-13)', mixed), ('fly-mean vectors (the reference)', ref)):
            a, b = H.align_pairs(kc, o)
            out.append(f'{name} ρ={H.spearman(a, b):.3f}')
        print(f'{lab} vs ORNs: ' + '; '.join(out))
    print('wrote', OUT / f'{STEM}_flymean.csv')


if __name__ == '__main__':
    main()
