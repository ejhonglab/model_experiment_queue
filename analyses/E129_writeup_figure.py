"""Write-up figure for E129 (the inhibition that outlasts the block reduced to a static level), from verdict.json only.

Two panels, no simulation:
- left, the class decision: each arm's calibrated pattern part minus the half-way point between the family and KC-KC
  at zero (E110 part A), on S's and T's populations, against all recordings and each of the four halves A to D. An arm
  "keeps the gain" when all ten points are at or below zero (`class_margins`);
- right, wiring by wiring: each arm's calibrated pattern part minus the family's on that wiring, against all recordings,
  with the 0.003 margin of the class "equal" (`wiring_diffs`).

Writes E129_pattern_part_classes.png, writeup_class_margins.csv and writeup_wiring_diffs.csv to the results folder,
with provenance_writeup_figure.json.

Usage (from the repository root): ./venv/bin/python analyses/E129_writeup_figure.py [--out DIR]
"""
import argparse
import json
import time
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E129_block_survivor_static'
ARMS = ('Gpre', 'G', 'T', 'P', 'S')
TARGETS = ('all', 'A', 'B', 'C', 'D')
POPULATIONS = ('S4', 'T4')
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
LABELS = {'family': 'the family', 'E110_A': 'KC–KC at zero', 'Gpre': 'G, preview\n(saved runs)',
          'G': 'G, one level', 'T': 'T, time-mean\nper KC\nand odor', 'P': 'P, pooled\nneuron\n(50 ms)',
          'S': 'S, matrix\nat 1 s'}
POP_STYLE = {'S4': dict(color='#2a78d6', marker='o', label='S (four searched wirings)'),
             'T4': dict(color='#d55181', marker='s', label='T (four held-out wirings)')}


def class_margins(verdict: dict, arms: Sequence[str] = ARMS, with_references: bool = True) -> pd.DataFrame:
    """Each arm's pattern part minus the half-way point, per target and population.

    Args:
        verdict: E129's verdict.json, read.
        arms: the arms to read (their 'classes' entries).
        with_references: also give the family and KC-KC at zero (E110 part A), read from the first arm's entry.

    Returns:
        One row per (kind, target, population): the arm's value, the half-way point, the margin (value minus half-way;
        at or below 0 on all ten rows is "keeps the gain") and the class.
    """
    rows = []
    for arm in arms:
        c = verdict['classes'][arm]
        for t in TARGETS:
            for pop in POPULATIONS:
                d = c['targets'][t]['populations'][pop]
                rows.append(dict(kind=arm, target=t, population=pop, value=d['arm'], half_way=d['half_way'],
                                 margin=d['arm'] - d['half_way'], cls=c['cls']))
    if with_references:
        c = verdict['classes'][arms[0]]
        for ref in ('family', 'E110_A'):
            for t in TARGETS:
                for pop in POPULATIONS:
                    d = c['targets'][t]['populations'][pop]
                    rows.append(dict(kind=ref, target=t, population=pop, value=d[ref], half_way=d['half_way'],
                                     margin=d[ref] - d['half_way'], cls=''))
    return pd.DataFrame(rows)


def wiring_diffs(verdict: dict, arms: Sequence[str] = ARMS) -> pd.DataFrame:
    """Each arm's pattern part minus the family's, wiring by wiring, against all recordings.

    Returns:
        One row per (kind, wiring) with the population the wiring belongs to (S4 for real, P1, P2, M1; T4 otherwise).
    """
    rows = []
    for arm in arms:
        for w, d in verdict['classes'][arm]['targets']['all']['d_family_per_wiring'].items():
            rows.append(dict(kind=arm, wiring=w, population='S4' if w in S_WIRINGS else 'T4', d_family=d))
    return pd.DataFrame(rows)


def draw(margins: pd.DataFrame, diffs: pd.DataFrame, path: Path, within: float = 0.003) -> None:
    """Draws the two panels to `path` (PNG)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[7, 5]))
    kinds = ['family', 'E110_A'] + [k for k in ARMS if k in set(margins.kind)]
    for i, k in enumerate(kinds):
        for j, pop in enumerate(POPULATIONS):
            st = POP_STYLE[pop]
            m = margins[(margins.kind == k) & (margins.population == pop)]
            x = i + (j - 0.5) * 0.3
            halves = m[m.target != 'all']
            a.scatter(np.full(len(halves), x), halves.margin, s=22, facecolors='none', edgecolors=st['color'],
                      marker=st['marker'], linewidths=1.2, zorder=3)
            al = m[m.target == 'all']
            a.scatter([x], al.margin, s=70, color=st['color'], marker=st['marker'], edgecolors='white',
                      linewidths=1.5, zorder=4, label=st['label'] if i == 0 else None)
    a.axhline(0, color='#52514e', lw=1, ls='--', zorder=1)
    a.text(-0.5, 0.0004, 'half-way point', ha='left', va='bottom', fontsize=8.5,
           color='#52514e')
    a.set_xticks(range(len(kinds)), [LABELS[k] for k in kinds], fontsize=8.5)
    a.set_xlim(-0.6, len(kinds) - 0.4)
    a.set_ylabel('pattern part minus the half-way point')
    a.set_title('The classes: at or below the line on all ten targets keeps the gain\n'
                '(filled: all recordings; open: the four halves)', fontsize=10)
    a.legend(loc='upper left', fontsize=8.5, frameon=False)

    arms = [k for k in ARMS if k in set(diffs.kind)]
    rng = np.random.default_rng(0)
    for i, k in enumerate(arms):
        for pop in POPULATIONS:
            st = POP_STYLE[pop]
            d = diffs[(diffs.kind == k) & (diffs.population == pop)]
            x = i + rng.uniform(-0.14, 0.14, len(d))
            b.scatter(x, d.d_family, s=34, color=st['color'], marker=st['marker'], edgecolors='white',
                      linewidths=1, zorder=3)
    b.axhline(0, color='#52514e', lw=1, zorder=1)
    b.axhline(within, color='#52514e', lw=1, ls=':', zorder=1)
    b.text(-0.55, within - 0.0004, f'"equal" margin (+{within})', ha='left', va='top', fontsize=8.5,
           color='#52514e')
    b.set_xticks(range(len(arms)), [LABELS[k] for k in arms], fontsize=8.5)
    b.set_xlim(-0.6, len(arms) - 0.4)
    b.set_ylabel('pattern part minus the family\'s, per wiring')
    b.set_title('Wiring by wiring, against all recordings', fontsize=10)
    for ax in (a, b):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e4e3df', lw=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='folder for the figure and tables (default: the results '
                                                            'folder, which also holds verdict.json)')
    args = ap.parse_args(argv)
    res = H.results_dir(EXPERIMENT)
    out = args.out or res
    out.mkdir(parents=True, exist_ok=True)
    verdict = json.loads((res / 'verdict.json').read_text())
    M, D = class_margins(verdict), wiring_diffs(verdict)
    M.to_csv(out / 'writeup_class_margins.csv', index=False)
    D.to_csv(out / 'writeup_wiring_diffs.csv', index=False)
    fig = out / 'E129_pattern_part_classes.png'
    draw(M, D, fig, within=verdict.get('rules', {}).get('within', 0.003))
    prov = dict(script=H.script_path(__file__), **H.code_commits(), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
                verdict_generated=verdict.get('generated'), outputs=[fig.name, 'writeup_class_margins.csv',
                                                                     'writeup_wiring_diffs.csv'])
    (out / 'provenance_writeup_figure.json').write_text(json.dumps(prov, indent=1, default=str))
    print(M[M.target == 'all'].pivot(index='kind', columns='population', values='margin').round(4))
    print(f'wrote {fig}')


if __name__ == '__main__':
    main()
