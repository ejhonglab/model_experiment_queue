#!/usr/bin/env python3
"""Yang's request (2026-11-18, viki's terminal): the best family's odor-odor cosine similarity with APL on, with APL off
and the change between the two, as heatmaps side by side, with the flies' three beside them.

"The best family" is the page's (E106's definition): the breadth-2 APL-effect fit (C2, E67) on its eight wirings (the
hemibrain wiring and seven re-drawn; E72's part B), read as the mean over the wirings of the cosine distance on raw spike
counts per condition (`E106.population`); APL off is the acute knockout (APL removed, nothing re-tuned). The flies' are
the mean over recordings of the cosine distance, control flies (TNT3c) for APL on and APL output blocked (TNT) for APL
off (`E106.data_matrices`). Similarity = 1 − distance; the change is APL off minus APL on, the page's per-pair APL effect
(positive where APL decorrelates the pair). Nothing is simulated.

Outputs: <shared>/results/yang_family_similarity/ (similarity.png, pairs.csv: the six values per pair, provenance.json).
Usage: ./venv/bin/python analyses/yang_family_similarity.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from kcq import harness as H
from experiments import E106_rho_gap_parts as E106

NAME = 'yang_family_similarity'


def similarity(dist: pd.Series) -> pd.Series:
    """1 − cosine distance, in canonical pair order."""
    return 1.0 - H.canonical_pairs(dist)


def pairs_table(model: Dict[str, pd.Series], data: Dict[str, pd.Series]) -> pd.DataFrame:
    """One row per pair: similarity with APL on / off and the change (off − on), for the family and for the flies."""
    cols = {}
    for who, d in (('model', model), ('flies', data)):
        on, off = similarity(d['on']), similarity(d['off'])
        cols[f'{who}_on'], cols[f'{who}_off'], cols[f'{who}_change'] = on, off, off - on
    return pd.DataFrame(cols).dropna()


def square(s: pd.Series, order: List[str]) -> np.ndarray:
    sq = np.array(H.square_from_pairs(s).reindex(index=order, columns=order).to_numpy(dtype=float))
    np.fill_diagonal(sq, np.nan)
    return sq


def group_bounds(order: List[str]) -> List[float]:
    """Positions between two odors of different chemical groups (`H.ODOR_GROUPS`)."""
    g = [H.ODOR_GROUPS.get(o, o) for o in order]
    return [i - 0.5 for i in range(1, len(order)) if g[i] != g[i - 1]]


def draw(t: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    order = H.odor_order(H.APL16_ODORS)
    seq = matplotlib.colormaps['Purples'].copy()
    seq.set_bad('#c8c8c8')
    div = H.oklab_diverging_cmap()
    div.set_bad('#c8c8c8')
    lim = float(np.ceil(t[['model_change', 'flies_change']].abs().to_numpy().max() * 20) / 20)
    rho = lambda a, b: spearmanr(t[a], t[b]).statistic
    titles = {
        'model_on': f'best family, APL on\nρ {rho("model_on", "flies_on"):.2f} with the control flies',
        'model_off': f'best family, APL off (acute)\nρ {rho("model_off", "flies_off"):.2f} with APL output blocked (TNT)',
        'model_change': f'best family: APL off minus APL on\nmean {t.model_change.mean():+.3f}, '
                        f'ρ {rho("model_change", "flies_change"):.2f} with the flies\' change',
        'flies_on': 'KCs, control flies (APL on)',
        'flies_off': 'KCs, APL output blocked (TNT) (APL off)',
        'flies_change': f'KCs: APL off minus APL on\nmean {t.flies_change.mean():+.3f}',
    }
    bounds = group_bounds(order)
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.6), gridspec_kw=dict(hspace=0.36, wspace=0.34))
    for i, who in enumerate(('model', 'flies')):
        for j, part in enumerate(('on', 'off', 'change')):
            ax, col = axes[i, j], f'{who}_{part}'
            mat = square(t[col], order)
            if part == 'change':
                im = ax.imshow(mat, cmap=div, vmin=-lim, vmax=lim)
                label = 'change in cosine similarity, APL off minus APL on'
            else:
                im = ax.imshow(mat, cmap=seq, vmin=0.0, vmax=1.0)
                label = 'cosine similarity'
            for b in bounds:
                ax.axhline(b, color='white', lw=1.2)
                ax.axvline(b, color='white', lw=1.2)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, rotation=90, fontsize=8)
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels(order, fontsize=8)
            ax.tick_params(length=0)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.set_title(titles[col], fontsize=10)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.set_label(label, fontsize=8)
            cb.ax.tick_params(labelsize=7)
            cb.outline.set_visible(False)
    fig.suptitle('Odor-odor cosine similarity: the best family (breadth-2 APL-effect fit, mean over its 8 wirings, raw '
                 'spike counts) and the KC data\nρ: Spearman correlation over the 120 pairs; odors in chemical groups '
                 '(white lines); grey diagonal', fontsize=11)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    t0 = time.time()
    out = H.shared_dir() / 'results' / NAME
    members = E106.family_members('family')
    model = E106.population(members)
    data = E106.data_matrices()
    t = pairs_table(model, data)
    if len(t) != 120:
        raise RuntimeError(f'{len(t)} pairs with all six values, not 120')
    out.mkdir(parents=True, exist_ok=True)
    t.rename_axis(['odor_a', 'odor_b']).to_csv(out / 'pairs.csv')
    draw(t, out / 'similarity.png')
    prov = H.model_provenance()
    prov.update(script=H.script_path(__file__), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), seconds=time.time() - t0,
                family=[dict(wiring=w, folder=f, tag=g) for w, f, g in members],
                data='E106.data_matrices (kcq.harness.yang_recording_rdms, cosine; TNT3c for APL on, TNT for APL off)')
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    rho = lambda a, b: spearmanr(t[a], t[b]).statistic
    print(f'{len(members)} wirings; ρ on {rho("model_on", "flies_on"):.3f}, off {rho("model_off", "flies_off"):.3f}, '
          f'change {rho("model_change", "flies_change"):.3f} (Pearson {t.model_change.corr(t.flies_change):.3f}); '
          f'mean change model {t.model_change.mean():+.3f}, flies {t.flies_change.mean():+.3f}')
    print(f'-> {out}')


if __name__ == '__main__':
    main()
