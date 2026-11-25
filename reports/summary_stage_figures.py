#!/usr/bin/env python3
"""The Summary's figures of the stage target (reports/summary_template.md), drawn from saved tables and the flies'
recordings; no model is run. build_summary.py runs it before filling the template.

Writes, in the shared directory's results/summary_stages/:
- stage_recordings.png (the background figure): the megamat flies at each recorded stage (ORNs, PN dendrites, PN
  boutons, KC claws, KC somata), one column each, in the three views: strength by odor (each odor's response rate as a
  multiple of the panel's mean; the flies' mean and range), sharing between odors (the flies' mean cosine distance
  between the sets of units that respond to two odors) and breadth among units (the share of units responding to a
  number of odors; the flies' mean and each fly). Read as E134 reads the flies (experiments/E134_flies.py, `reading`
  and `edges`; kcq/population_fit.py, `fly_reference`); stage_recordings_values.csv holds the plotted values;
- stage_strength_by_model.png: strength by odor at each megamat stage after the ORNs and as the megamat stage error, in
  one-fly units: the flies' floor and most atypical fly, the recordings upstream taken as the prediction, the best
  family held at the flies' response rate (each wiring and the wirings' mean), and the three PN stages of E135 (the
  recorded dendrites imposed, R1 and R3). Sources: E134's writeup_stage_table.csv, writeup_steps.csv and
  writeup_stage_error.csv; E135's writeup_values.csv;
- strength_by_odor_somata.png: strength by odor at the megamat KC somata, odor by odor, the flies against the family
  as fitted and held (E138's entries.csv), and aa along the four stages (E138's stages.csv).

Usage: ./venv/bin/python reports/summary_stage_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from kcq import harness as H

RESULTS = H.shared_dir() / 'results'
OUT = RESULTS / 'summary_stages'
STAGES = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
STAGE_NAMES = {'orn': 'ORNs', 'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
               'kc_soma': 'KC somata'}
ODORS = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms', 'B-cit', 'Lin', 'benz',
         'aa', 'va')
WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
# the marks of E137's strength_by_reading.png where an element is shared: dark floor, orange most atypical fly, blue family
FLOOR, ATYPICAL, FAMILY, FLIES = '#3a3a38', '#eb6834', '#2a78d6', '#3a3a38'
ARMS = (('imposed', 'the recorded dendrites imposed', '#6250d6'),
        ('R1', 'fitted to the dendrites: the global stage', '#d55181'),
        ('rung', 'fitted to the dendrites: the stage with glomerulus terms', '#1baf7a'))
BLUES = ('#f4f8fe', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b')


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    return plt


def _bin_label(b: str) -> str:
    lo, _, hi = b.partition('_')
    return f'{lo}-{hi}' if hi else lo


def stage_recordings() -> pd.DataFrame:
    """The background figure and its values (module docstring)."""
    from experiments import E134_flies as F
    from kcq import population_fit as P
    from matplotlib.colors import LinearSegmentedColormap
    plt = _plt()
    refs = {}
    for st in STAGES:
        flies, odors = F.reading((st, 'megamat', 'wild_type'))
        refs[st] = P.fly_reference(flies, odors, F.edges(len(odors)))
    fig, axes = plt.subplots(3, 5, figsize=(17, 9.6), gridspec_kw=dict(height_ratios=[1, 1.25, 0.9]))
    cmap = LinearSegmentedColormap.from_list('blues', BLUES)
    rows, ymax, im = [], 0.0, None
    for j, st in enumerate(STAGES):
        ref = refs[st]
        n = len(ref.flies)
        # strength by odor
        ax = axes[0, j]
        mean = ref.mean['strength/rate'].reindex(ODORS)
        parts = ref.parts['strength/rate'].reindex(columns=ODORS)
        lo, hi = parts.min(), parts.max()
        x = np.arange(len(ODORS))
        ax.vlines(x, lo, hi, color=FLIES, alpha=0.35, lw=3)
        ax.plot(x, mean, 'o', color=FLIES, ms=5)
        ax.axhline(1, color='0.55', lw=1, ls='--')
        ax.set_xticks(x, ODORS, rotation=90, fontsize=8)
        ax.set_title(f'{STAGE_NAMES[st]} ({n} flies)', fontsize=11)
        ymax = max(ymax, float(np.nanmax(hi.to_numpy())))
        for o in ODORS:
            rows.append(dict(stage=st, part='strength/rate', entry=o, flies_mean=mean.get(o), flies_min=lo.get(o),
                             flies_max=hi.get(o), n_flies=n))
        # sharing between odors, who responds
        ax = axes[1, j]
        sh = ref.mean['sharing/who']
        M = pd.DataFrame(np.nan, index=list(ODORS), columns=list(ODORS))
        for k, v in sh.items():
            a, b = k.split('|')
            if a in M.index and b in M.index:
                M.loc[a, b] = M.loc[b, a] = v
            rows.append(dict(stage=st, part='sharing/who', entry=k, flies_mean=v,
                             flies_min=ref.parts['sharing/who'][k].min(), flies_max=ref.parts['sharing/who'][k].max(),
                             n_flies=n))
        im = ax.imshow(M.to_numpy(dtype=float), cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
        ax.set_xticks(x, ODORS, rotation=90, fontsize=7)
        ax.set_yticks(x, ODORS, fontsize=7)
        ax.spines[:].set_visible(False)
        # breadth among units
        ax = axes[2, j]
        br = ref.mean['breadth/bins']
        bp = ref.parts['breadth/bins']
        xb = np.arange(len(br))
        ax.bar(xb, br.to_numpy(dtype=float), color='0.82', width=0.7)
        for f in bp.index:
            ax.plot(xb + np.linspace(-0.2, 0.2, len(bp))[list(bp.index).index(f)], bp.loc[f].to_numpy(dtype=float),
                    'o', color=FLIES, ms=3.5)
        ax.set_xticks(xb, [_bin_label(b) for b in br.index])
        ax.set_xlabel('number of odors a unit responds to', fontsize=9)
        ax.set_ylim(0, 1)
        for b in br.index:
            rows.append(dict(stage=st, part='breadth/bins', entry=b, flies_mean=br[b], flies_min=bp[b].min(),
                             flies_max=bp[b].max(), n_flies=n))
    for j in range(5):
        axes[0, j].set_ylim(0, ymax * 1.05)
        axes[0, j].grid(axis='y', color='0.92', lw=0.8)
        axes[0, j].set_axisbelow(True)
        if j:
            axes[0, j].set_yticklabels([])
            axes[1, j].set_yticklabels([])
            axes[2, j].set_yticklabels([])
    axes[0, 0].set_ylabel("response rate,\nmultiple of the panel's mean")
    axes[2, 0].set_ylabel('share of units')
    fig.tight_layout()
    cb = fig.colorbar(im, ax=axes[1, :].tolist(), shrink=0.85, pad=0.01)
    cb.set_label('cosine distance across units\n(which units respond)')
    fig.savefig(OUT / 'stage_recordings.png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    V = pd.DataFrame(rows)
    V.to_csv(OUT / 'stage_recordings_values.csv', index=False)
    return V


def stage_strength_by_model() -> None:
    plt = _plt()
    e134 = RESULTS / 'E134_population_baseline'
    st_tab = pd.read_csv(e134 / 'writeup_stage_table.csv')
    st_tab = st_tab[(st_tab.group == 'megamat_fitted') & (st_tab.panel == 'megamat') & (st_tab.condition == 'on')]
    steps = pd.read_csv(e134 / 'writeup_steps.csv')
    steps = steps[steps.group == 'megamat_fitted'].set_index('step')
    err = pd.read_csv(e134 / 'writeup_stage_error.csv')
    held = err[err.group == 'megamat_held'].set_index('wiring')
    vals = pd.read_csv(RESULTS / 'E135_antennal_lobe_dendrites' / 'writeup_values.csv')
    vals = vals[(vals.panel == 'megamat') & (vals.condition == 'on')]
    groups = list(STAGES[1:]) + ['stage_error']
    names = [STAGE_NAMES[s] for s in STAGES[1:]] + ['megamat\nstage error']
    fig, ax = plt.subplots(figsize=(11, 4.6))
    seen = set()

    def lab(s):
        if s in seen:
            return None
        seen.add(s)
        return s
    for i, g in enumerate(groups):
        if g != 'stage_error':
            r = st_tab[st_tab.stage == g].iloc[0]
            ax.hlines(r.strength_floor, i - 0.42, i + 0.42, color=FLOOR, lw=3, label=lab("the flies' floor"))
            ax.hlines(r.strength_loo_max, i - 0.42, i + 0.42, color=ATYPICAL, lw=3, label=lab('most atypical fly'))
            up = STAGES[STAGES.index(g) - 1]
            ns = steps.loc[f'{up}->{g}', 'no_step_strength']
            ax.plot(i - 0.3, ns, 's', color='0.45', ms=7, mec='white', label=lab('the recordings upstream taken as the prediction'))
            col = f'{g}_strength'
            stat, stage = 'strength', g
        else:
            col = 'stage_error_strength'
            stat, stage = 'stage_error_strength', 'four_stages'
        w = held.reindex(list(WIRINGS))[col].to_numpy(dtype=float)
        ax.plot(np.full(len(w), i - 0.12) + np.linspace(-0.05, 0.05, len(w)), w, 'o', color=FAMILY, ms=3.5, alpha=0.6,
                label=lab('the best family, each wiring'))
        ax.plot(i - 0.12, held.loc['model', col], 'o', color=FAMILY, ms=9, mec='white', mew=1.2,
                label=lab("the best family, wirings' mean"))
        for k, (arm, name, c) in enumerate(ARMS):
            v = vals[(vals.arm == arm) & (vals.stage == stage)].set_index('statistic').value
            m, lo, hi = (v.get(f'{stat}_mean'), v.get(f'{stat}_min'), v.get(f'{stat}_max'))
            if m is None or not np.isfinite(m):
                continue
            xx = i + 0.1 + 0.1 * k
            ax.vlines(xx, lo, hi, color=c, lw=1.2)
            ax.plot(xx, m, '^', color=c, ms=8, mec='white', mew=1, label=lab(name))
    ax.axhline(1, color='0.55', lw=1, ls=':')
    ax.set_xticks(range(len(groups)), names)
    ax.set_ylim(0, None)
    ax.set_ylabel('strength by odor,\none-fly units')
    ax.grid(axis='y', color='0.92', lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'stage_strength_by_model.png', dpi=140, bbox_inches='tight')
    plt.close(fig)


def strength_by_odor_somata() -> None:
    plt = _plt()
    e138 = RESULTS / 'E138_strength_by_odor_entries'
    ent = pd.read_csv(e138 / 'entries.csv')
    ent = ent[(ent.panel == 'megamat') & (ent.condition == 'on') & (ent.reading == 'e134') & (ent.flies == 'all')
              & (ent.row == 'common')]
    fit = ent[ent.group == 'family_fitted'].set_index('odor')
    hel = ent[ent.group == 'family_held'].set_index('odor')
    order = fit.rate_flies.sort_values(ascending=False).index.tolist()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw=dict(width_ratios=[3, 1]), sharey=True)
    x = np.arange(len(order))
    a1.vlines(x, fit.rate_flies_min[order], fit.rate_flies_max[order], color=FLIES, alpha=0.35, lw=3)
    a1.plot(x, fit.rate_flies[order], 'o', color=FLIES, ms=6, label='the flies (mean and range)')
    a1.plot(x - 0.15, fit.rate_model[order], 'D', color=FAMILY, ms=6, mec='white', label='the best family, as fitted')
    a1.plot(x + 0.15, hel.rate_model.reindex(order), 'D', mfc='white', mec=FAMILY, mew=1.5, ms=6,
            label="the best family, held at the flies' response rate")
    a1.axhline(1, color='0.55', lw=1, ls='--')
    a1.set_xticks(x, order, rotation=90)
    a1.set_title('KC somata, megamat')
    a1.set_ylabel("response rate,\nmultiple of the panel's mean")
    a1.legend(frameon=False, fontsize=9)
    st = pd.read_csv(e138 / 'stages.csv')
    st = st[(st.panel == 'megamat') & (st.condition == 'on') & (st.row == 'common') & (st.odor == 'aa')]
    stages = [s for s in STAGES[1:] if s in set(st.stage)]
    xs = np.arange(len(stages))
    f = st[st.group == 'family_fitted'].set_index('stage').reindex(stages)
    h = st[st.group == 'family_held'].set_index('stage').reindex(stages)
    a2.vlines(xs, f.rate_flies_min, f.rate_flies_max, color=FLIES, alpha=0.35, lw=3)
    a2.plot(xs, f.rate_flies, 'o-', color=FLIES, ms=6, lw=1)
    a2.plot(xs, f.rate_model, '-', color=FAMILY, lw=2)
    if h.rate_model.notna().any():
        a2.plot(xs, h.rate_model, '--', color=FAMILY, lw=2)
    a2.axhline(1, color='0.55', lw=1, ls='--')
    a2.set_xticks(xs, [STAGE_NAMES[s] for s in stages], rotation=90)
    a2.set_title('aa along the stages')
    for a in (a1, a2):
        a.grid(axis='y', color='0.92', lw=0.8)
        a.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(OUT / 'strength_by_odor_somata.png', dpi=140, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    stage_recordings()
    stage_strength_by_model()
    strength_by_odor_somata()
    prov = dict(script='reports/summary_stage_figures.py', **H.code_commits(),
                figures=['stage_recordings.png', 'stage_strength_by_model.png', 'strength_by_odor_somata.png'])
    (OUT / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print(f'-> {OUT}')


if __name__ == '__main__':
    main()
