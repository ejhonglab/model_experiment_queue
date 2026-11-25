#!/usr/bin/env python3
"""The components' fall across Remy's mixture sessions (a read-only follow-up of E137, at Tom's request of 2026-11-25).

Tom asked for a page on how much a component's second block (given after the mixtures) evokes against its first, before
he decides whether the mixture panels enter E136's "holds". No model is run. Everything is read from E137's locked
per-trial files (experiments/E137_lock.json, sha256 checked) in its baseline-noise units and at megamat's soma level,
and from Remy's four megamat flies (E137's `megamat_order_trials`), whose order is drawn anew per block.

Parts, written to results/E137_heldout_readings/component_fall/:
- components.csv: per mixture fly and component with a second block, each block's response rate and mean response, the
  second block over the first with a 90 % interval from resampling the fly's cells, the component's place in the second
  block, and the cross-adaptation check (below);
- presentations.csv: each of a component's six presentations, its response rate and mean response, and both over the
  first presentation's;
- megamat_repeats.csv: per megamat fly, each odor's repeats over its first, and a fit of each trial's log response rate
  (and log mean response) on odor effects plus a factor per repeat of the odor plus a factor per place in the block, so
  that a fall with the odor's own repetition and a fall with time are told apart ('fall_with': 'repeats' where a
  repeat costs more response rate than the 17 places between repeats, else 'time');
- mixture_ratio.csv: per mixture fly, the mixture at 0 over the components' mean response rate under each of E137's
  readings (from E137's rise.csv), beside the family's raw counts;
- fall_by_component.png, fall_presentations.png; provenance.json.

The cross-adaptation check: over a fly's cells, the second block's mean response to a component regressed (least squares)
on its first block's and on the cell's response to the mixture at 0; a negative coefficient on the mixture says that the
cells the mixture drove lost more of their response to the component than other cells with the same first-block response.

Usage: ./venv/bin/python analyses/E137_component_fall.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E137_heldout_readings as E

N_BOOT = 1000
SEED = 137
COLORS = {'kmix': '#2a78d6', 'cmix': '#d55181'}
PANEL_NAMES = {'kmix': 'kiwi flies', 'cmix': 'control flies'}


def out_dir() -> Path:
    d = E.out_dir() / 'component_fall'
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fly_frames(fly: str, lock: dict):
    """One fly's cells x (label, block, rep) single-trial values in baseline-noise units (E137's footing), its session's
    trial order (one row per trial), or None for a fly whose files are missing."""
    red = E.fly_reduced(fly, lock)
    if red is None:
        return None
    red = red.assign(v=E.noise_values(red, lock['footing']['factor']))
    W = red.pivot_table(index='cell', columns=['label', 'block', 'rep'], values='v')
    seq = red.drop_duplicates('trial').sort_values('trial')[['trial', 'label', 'block', 'rep', 'in_panel']]
    return W, seq


def _boot_ratio(m0: np.ndarray, m1: np.ndarray, crit: float, rng: np.random.Generator):
    """90 % intervals of the second block over the first, for the response rate and the mean, from resampling cells."""
    n = len(m0)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    a0, a1 = m0[idx], m1[idx]
    r0, r1 = (a0 > crit).mean(1), (a1 > crit).mean(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        rr = np.where(r0 > 0, r1 / r0, np.nan)
        mr = a1.mean(1) / a0.mean(1)
    q = lambda x: (float(np.nanquantile(x, 0.05)), float(np.nanquantile(x, 0.95)))
    return q(rr), q(mr)


def part_components(lock: dict):
    crit = float(lock['criteria']['noise'])
    rng = np.random.default_rng(SEED)
    rows, prow = [], []
    for fly, rec in lock['flies'].items():
        fr = _fly_frames(fly, lock)
        if fr is None:
            continue
        W, seq = fr
        panel = rec['panel']
        mix = f'{panel}0'
        second = seq[(seq.block == 1) & seq.in_panel & ~seq.label.str.contains('mix')]
        comps = list(dict.fromkeys(second.label))
        mixm = W[mix][0].mean(axis=1).to_numpy()
        first_block_order = list(dict.fromkeys(seq[(seq.block == 0) & seq.in_panel & ~seq.label.str.contains('mix')].label))
        for c in (comps or first_block_order):
            b0 = W[c][0]
            m0 = b0.mean(axis=1).to_numpy()
            pres = [b0[k].to_numpy() for k in sorted(b0.columns)]
            if c in comps:
                b1 = W[c][1]
                pres += [b1[k].to_numpy() for k in sorted(b1.columns)]
            for k, x in enumerate(pres):
                prow.append(dict(fly=fly, panel=panel, component=c, presentation=k + 1, block=k // 3 + 1,
                                 rate=float((x > crit).mean()), mean=float(x.mean())))
            if c not in comps:
                continue
            m1 = b1.mean(axis=1).to_numpy()
            r0, r1 = float((m0 > crit).mean()), float((m1 > crit).mean())
            (rlo, rhi), (mlo, mhi) = _boot_ratio(m0, m1, crit, rng)
            X = np.column_stack([np.ones_like(m0), m0, mixm])
            beta = np.linalg.lstsq(X, m1, rcond=None)[0]
            t_first = int(seq[(seq.label == c) & (seq.block == 1)].trial.min())
            t_mix_end = int(seq[seq.label.str.contains('mix') | seq.label.str.contains('_no_')].trial.max())
            rows.append(dict(fly=fly, panel=panel, component=c, place_in_second_block=comps.index(c) + 1,
                             trials_after_last_mixture=t_first - t_mix_end - 1,
                             trials_between_blocks=t_first - int(seq[(seq.label == c) & (seq.block == 0)].trial.min()),
                             n_cells=int(len(m0)), rate_block1=r0, rate_block2=r1,
                             rate_ratio=r1 / r0 if r0 > 0 else np.nan, rate_ratio_lo=rlo, rate_ratio_hi=rhi,
                             mean_block1=float(m0.mean()), mean_block2=float(m1.mean()),
                             mean_ratio=float(m1.mean() / m0.mean()), mean_ratio_lo=mlo, mean_ratio_hi=mhi,
                             mixture_coef=float(beta[2]), first_block_coef=float(beta[1]),
                             corr_blocks=float(np.corrcoef(m0, m1)[0, 1])))
    C = pd.DataFrame(rows)
    P = pd.DataFrame(prow)
    for q in ('rate', 'mean'):
        P[f'{q}_over_first'] = P[q] / P.groupby(['fly', 'component'])[q].transform('first')
    C.to_csv(out_dir() / 'components.csv', index=False)
    P.to_csv(out_dir() / 'presentations.csv', index=False)
    return C, P


def part_megamat():
    t, crit = E.megamat_order_trials()
    tr = t.drop_duplicates(['fly', 'trial']).sort_values(['fly', 'trial'])
    tr = tr.assign(rep=tr.groupby(['fly', 'odor']).cumcount())
    t = t.merge(tr[['fly', 'trial', 'rep']], on=['fly', 'trial'])
    per = t.groupby(['fly', 'trial', 'odor', 'rep', 'pos']).value.agg(
        rate=lambda s: float((s > crit).mean()), mean='mean').reset_index()
    rows = []
    for fly, g in per.groupby('fly'):
        rec = dict(fly=fly, n_trials=int(len(g)), n_odors=int(g.odor.nunique()), trials_per_block=int(g.pos.max()) + 1)
        for q in ('rate', 'mean'):
            P = g.pivot_table(index='odor', columns='rep', values=q)
            tot = P.mean()
            rec[f'{q}_repeat2_over_1'] = float(tot[1] / tot[0])
            rec[f'{q}_repeat3_over_1'] = float(tot[2] / tot[0])
            y = np.log(np.clip(g[q].to_numpy(dtype=float), 1e-3, None))
            X = np.column_stack([pd.get_dummies(g.odor).to_numpy(dtype=float), g.rep.to_numpy(dtype=float),
                                 g.pos.to_numpy(dtype=float)])
            b = np.linalg.lstsq(X, y, rcond=None)[0]
            rec[f'{q}_factor_per_repeat'] = float(np.exp(b[-2]))
            rec[f'{q}_factor_per_block_of_places'] = float(np.exp(b[-1] * rec['trials_per_block']))
        # what the fall follows: the odor's repeats if a repeat costs more than as many trials of other odors
        rec['fall_with'] = ('repeats' if rec['rate_factor_per_repeat'] < rec['rate_factor_per_block_of_places']
                            else 'time')
        rows.append(rec)
    M = pd.DataFrame(rows)
    M.to_csv(out_dir() / 'megamat_repeats.csv', index=False)
    return M


def part_mixture_ratio():
    r = pd.read_csv(E.out_dir() / 'rise.csv')
    f = r[(r.side == 'flies') & (r.units == 'noise')]
    W = f.pivot_table(index=['panel', 'wiring'], columns='reading', values='rate_mix0_over_components').reset_index()
    W = W.rename(columns={'wiring': 'fly'})
    fam = r[(r.reading == 'raw_counts') & (r.group == 'fitted')].groupby('panel').rate_mix0_over_components.median()
    W['family_raw_counts_median'] = W.panel.map(fam)
    W.to_csv(out_dir() / 'mixture_ratio.csv', index=False)
    return W


def figures(C: pd.DataFrame, P: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    # 1: the second block over the first, per component and fly
    order = [('kmix', c) for c in ('EtOH', 'IAol', 'IaA', 'ea', 'eb')] + [('cmix', c) for c in ('ms', 'va', 'fur', '2h', '1o3ol')]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, q, title in ((axes[0], 'rate', 'response rate'), (axes[1], 'mean', 'mean response')):
        for i, (panel, c) in enumerate(order):
            x = C[(C.panel == panel) & (C.component == c)]
            if x.empty:
                continue
            off = np.linspace(-0.25, 0.25, len(x)) if len(x) > 1 else np.zeros(1)
            xs = i + (1 if panel == 'cmix' else 0) + off
            y = x[f'{q}_ratio'].to_numpy()
            ax.vlines(xs, x[f'{q}_ratio_lo'], x[f'{q}_ratio_hi'], color=COLORS[panel], lw=1.2, alpha=0.6)
            ax.plot(xs, y, 'o', ms=6, color=COLORS[panel], mec='white', mew=1.2)
        ax.axhline(1, color='0.6', lw=1, ls=':')
        ticks = [i + (1 if p == 'cmix' else 0) for i, (p, _) in enumerate(order)]
        ax.set_xticks(ticks, [c for _, c in order])
        ax.set_title(f'{title}: second block over first')
        ax.grid(axis='y', color='0.92', lw=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('second block / first block')
    axes[0].set_ylim(0, 1.3)
    for panel, xm in (('kmix', 2), ('cmix', 8)):
        axes[0].text(xm, 1.22, PANEL_NAMES[panel], ha='center', color='0.25')
        axes[1].text(xm, 1.22, PANEL_NAMES[panel], ha='center', color='0.25')
    fig.text(0.5, -0.02, 'components in their order in each block; one dot per fly, bars 90 % intervals from resampling '
             'its cells', ha='center', color='0.35', fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir() / 'fall_by_component.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    # 2: the six presentations, over the first
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    lead = {'kmix': 'EtOH', 'cmix': 'ms'}
    for ax, panel in zip(axes, ('kmix', 'cmix')):
        p = P[(P.panel == panel) & P.component.isin(C[C.panel == panel].component.unique())
              & P.fly.isin(C[C.panel == panel].fly.unique())]
        med = p.groupby(['component', 'presentation']).rate_over_first.median().unstack('presentation')
        xs = np.array([1, 2, 3, 5, 6, 7])
        ends = med.iloc[:, -1].sort_values()
        ypos, last = {}, -1.0
        for c, y in ends.items():           # end labels at least 0.045 apart
            last = max(float(y), last + 0.045)
            ypos[c] = last
        for c, row in med.iterrows():
            hi = c == lead[panel]
            ax.plot(xs[:row.size], row.to_numpy(), '-o', ms=6 if hi else 4, lw=2,
                    color=COLORS[panel] if hi else '0.65', mec='white', mew=1, zorder=3 if hi else 2)
            ax.text(7.25, ypos[c], c, va='center', fontsize=9, color='0.2' if hi else '0.45')
        ax.axvspan(3.5, 4.5, color='0.94', lw=0)
        ax.text(4, 1.12, 'the mixtures', ha='center', fontsize=9, color='0.35')
        ax.set_xticks(xs, ['1', '2', '3', '1', '2', '3'])
        ax.set_xlabel('presentation: first block | second block')
        ax.set_xlim(0.6, 7.9)
        ax.set_title(f"{PANEL_NAMES[panel]} (median over the {p.fly.nunique()} flies with a second block)")
        ax.grid(axis='y', color='0.92', lw=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("response rate / the first presentation's")
    axes[0].set_ylim(0, 1.2)
    fig.tight_layout()
    fig.savefig(out_dir() / 'fall_presentations.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    lock = E.load_lock()
    C, P = part_components(lock)
    M = part_megamat()
    R = part_mixture_ratio()
    figures(C, P)
    prov = dict(script='analyses/E137_component_fall.py', time=E._now(), host=E._host(), **H.code_commits(),
                lock=E._rel(E.LOCK_COPY), criterion_noise=lock['criteria']['noise'], footing=lock['footing']['factor'],
                n_boot=N_BOOT, seed=SEED, n_flies=int(C.fly.nunique()), n_megamat_flies=int(len(M)))
    (out_dir() / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    pd.set_option('display.width', 200)
    print(C.round(3).to_string())
    print(M.round(3).to_string())
    print(R.round(3).to_string())


if __name__ == '__main__':
    main()
