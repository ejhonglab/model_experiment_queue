"""Write-up tables and figure for E131 (MB-C1's two cells beside the inhibition among KCs), from E131's results folder
and the data-target registry only; no simulation.

The numbers E131's note gives that the score's tables do not give directly:
- `class_margins`: each arm's calibrated pattern part, the family's, KC-KC at zero's (E110 part A) and the half-way point
  between them, per target (all recordings, halves A to D) and population (S4, T4), with the arm minus the family and
  minus the half-way point (verdict.json's class targets);
- `wiring_summary`: per arm and target, the arm minus the family wiring by wiring (verdict.json's d_family_per_wiring),
  with the number of wirings at or under the 0.003 margin, under the family, and under it by 0.003 or more, and the
  range over the eight wirings;
- `cost_split`: for each pair (a, b) of the "differ" comparisons and for L100, D100 and P100 against B100, b minus the
  family, a minus the family, b minus a and the share of b's rise over the family that a removes, per target and
  population;
- `run_ranges`: per arm and statistic of runs.csv, the mean over the eight wirings, over S's and T's four, the range
  over the eight (with the wirings at the ends) and, for the threshold, the wirings under 60 Hz (E126's line);
- `yardsticks`: step 0's dose arithmetic from the output table (reading 13) beside the same arithmetic from the input
  table (what the entry used), MB-C1's synapses against APL's on the claws and in the calyx, the per-synapse weight of
  each arm against a KC->KC synapse, step 0's fixed numbers (D100's threshold, L100's ratio, the cells' agreement) and
  activities, each gamma KC's MB-C1 input at m 1 against its KC-KC input, and the KC data's values (the registry's) the
  note sets beside the model's;
- `differ_extra`: comparisons the entry did not fix, by the experiment's own "differ" rule (E129's `differ` through
  experiments/E131_mbc1_cells.py, as the score's comparisons use it, on readings.csv and runs.csv): each MB-C1 arm
  against the family, and the slow and the thresholded cells (L100, D100) against the plain cells at a quarter of the
  dose (B25). Per comparison: the verdict, the wirings on which each is lower (against all recordings), whether each
  is lower on all ten targets, the corner arms, a minus b over all eight wirings, on S and T, and wiring by wiring.

The figure (E131_pattern_part_classes.png): left, each arm's calibrated pattern part minus the family's on S's and T's
populations, against all recordings (filled) and each half (open), with the 0.003 margin; right, wiring by wiring
against all recordings.

Checks (raise on failure): the mean over the eight wirings of each runs.csv statistic equals readings.csv's all8 run_
column; verdict.json's d_family equals arm minus family; the per-synapse weights computed here equal runs.csv's
mbc1_per_synapse_vs_kckc; the step-0 totals equal step0_synapses.csv's; the score's own comparisons, recomputed by the
same rule from readings.csv and runs.csv, equal verdict.json's.

Writes writeup_class_margins.csv, writeup_wiring_summary.csv, writeup_cost_split.csv, writeup_run_ranges.csv,
writeup_yardsticks.csv, writeup_differ_extra.csv, E131_pattern_part_classes.png and provenance_writeup_tables.json to
the results folder.

Usage (from the repository root): ./venv/bin/python analyses/E131_writeup_tables.py [--out DIR]
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E131_mbc1_cells'
WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
TARGETS = ('all', 'A', 'B', 'C', 'D')
POPULATIONS = ('S4', 'T4')
REFS = ('family', 'E110_A', 'PO50', 'TY50')
ARMS = ('B25', 'B100', 'B400', 'R50', 'H50', 'X50', 'L100', 'D100', 'P100')
CLASSED = ('PO50', 'TY50') + ARMS
FORM_PAIRS = (('L100', 'B100'), ('D100', 'B100'), ('P100', 'B100'), ('R50', 'H50'), ('X50', 'PO50'), ('X50', 'TY50'),
              ('B25', 'B100'), ('B100', 'B400'), ('R50', 'B100'))
MBC1_ARMS = ('B25', 'B100', 'B400', 'R50', 'X50', 'L100', 'D100', 'P100')
# the "differ" comparisons the entry did not fix, (a, b)
EXTRA_PAIRS = tuple((a, 'family') for a in MBC1_ARMS) + (('L100', 'B25'), ('D100', 'B25'))
DIFFER_KEYS = ('verdict', 'n_wirings_a_lower', 'n_wirings_b_lower', 'a_lower_on_all_targets', 'b_lower_on_all_targets')
RUN_STATS = ('fixed_thr', 'wAPLKC_scale', 'wKCKC_scale', 'wMBC1_scale', 'rate_on', 'rate_off', 'ratio', 'rate_rho_on',
             'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on', 'max_count_off', 'rate_on_g', 'rate_off_g',
             'ratio_g', 'rate_on_apbp', 'rate_off_apbp', 'rate_on_ab', 'rate_off_ab', 'mbc1_act_on', 'mbc1_act_off',
             'mbc1_peak_on', 'mbc1_peak_off', 'share_mbc1_vs_apl_g', 'share_mbc1_vs_apl_all', 'share3_mbc1_g',
             'share3_apl_g', 'share3_kckc_g', 'mbc1_per_synapse_vs_kckc', 'thr_factor', 'apl_factor')
THR_LINE = 60.0
TOL = 1e-9
REGISTRY = ('noise_ctrl_frac', 'noise_tnt_frac', 'noise_rate_ratio_fraction', 'noise_ctrl_odors_per_kc_tv_max',
            'noise_tnt_odors_per_kc_tv_max', 'apl_effect_mean', 'apl_effect_pattern_floor', 'apl_effect_rmse_floor',
            'apl_effect_rho_ceiling', 'ctrl_rate_rho_ceiling', 'tnt_rate_rho_ceiling', 'ctrl_ceiling', 'tnt_ceiling',
            'ctrl_mean_dist', 'tnt_mean_dist')
LABELS = {'E110_A': 'KC–KC\nat zero', 'PO50': 'PO50\nE125', 'TY50': 'TY50\nE125', 'B25': 'B25', 'B100': 'B100',
          'B400': 'B400\ncorner', 'R50': 'R50', 'H50': 'H50', 'X50': 'X50', 'L100': 'L100\n1 s',
          'D100': 'D100\nthresh-\nold', 'P100': 'P100\npooled'}
POP_STYLE = {'S4': dict(color='#2a78d6', marker='o', label='S (real, P1, P2, M1)'),
             'T4': dict(color='#d55181', marker='s', label='T (P3, P4, M2, M3)')}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def class_margins(verdict: dict) -> pd.DataFrame:
    """Each classed arm's pattern part against the family, KC-KC at zero and the half-way point, per target and
    population (calibrated readout), with the references as rows of their own."""
    rows = []
    for arm in CLASSED:
        c = verdict['classes'][arm]
        for t in TARGETS:
            for pop in POPULATIONS:
                d = c['targets'][t]['populations'][pop]
                if abs(d['arm'] - d['family'] - d['d_family']) > TOL:
                    raise ValueError(f'{arm} {t} {pop}: d_family is not arm minus family')
                rows.append(dict(kind=arm, target=t, population=pop, value=d['arm'], family=d['family'],
                                 E110_A=d['E110_A'], half_way=d['half_way'], minus_family=d['d_family'],
                                 minus_half_way=d['arm'] - d['half_way'], cls=c['cls'],
                                 keeps_gain=c.get('keeps_gain')))
    c = verdict['classes'][CLASSED[0]]
    for ref in ('family', 'E110_A'):
        for t in TARGETS:
            for pop in POPULATIONS:
                d = c['targets'][t]['populations'][pop]
                rows.append(dict(kind=ref, target=t, population=pop, value=d[ref], family=d['family'],
                                 E110_A=d['E110_A'], half_way=d['half_way'], minus_family=d[ref] - d['family'],
                                 minus_half_way=d[ref] - d['half_way'], cls='', keeps_gain=None))
    return pd.DataFrame(rows)


def wiring_summary(verdict: dict, within: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per arm and target, the arm minus the family wiring by wiring, and its counts and range over the wirings.

    Returns:
        (per-wiring rows: kind, target, wiring, population, minus_family; summary rows: kind, target, n_within (at or
        under +within), n_lower (under 0), n_below (at or under -within), min, wiring_min, max, wiring_max).
    """
    per, summ = [], []
    for arm in CLASSED:
        for t in TARGETS:
            dw = verdict['classes'][arm]['targets'][t]['d_family_per_wiring']
            s = pd.Series({w: dw[w] for w in WIRINGS})
            for w, v in s.items():
                per.append(dict(kind=arm, target=t, wiring=w, population='S4' if w in S_WIRINGS else 'T4',
                                minus_family=v))
            summ.append(dict(kind=arm, target=t, n_within=int((s <= within).sum()), n_lower=int((s < 0).sum()),
                             n_below=int((s <= -within).sum()), min=s.min(), wiring_min=s.idxmin(), max=s.max(),
                             wiring_max=s.idxmax()))
    return pd.DataFrame(per), pd.DataFrame(summ)


def cost_split(verdict: dict, pairs: Sequence[Tuple[str, str]] = FORM_PAIRS) -> pd.DataFrame:
    """For each pair (a, b): b minus the family, a minus the family, b minus a, and the share of b's rise over the
    family that a removes ((b - a) / (b - family)), per target and population, with the "differ" verdict where the
    score made one."""
    rows = []
    comps = verdict.get('comparisons', {})
    for a, b in pairs:
        verdict_ab = comps.get(f'{a}|{b}', {}).get('verdict', '')
        for t in TARGETS:
            for pop in POPULATIONS:
                da = verdict['classes'][a]['targets'][t]['populations'][pop]['d_family']
                db = verdict['classes'][b]['targets'][t]['populations'][pop]['d_family']
                rows.append(dict(a=a, b=b, target=t, population=pop, b_minus_family=db, a_minus_family=da,
                                 b_minus_a=db - da, share_removed=(db - da) / db if db != 0 else np.nan,
                                 differ_verdict=verdict_ab))
    return pd.DataFrame(rows)


def run_ranges(runs: pd.DataFrame, readings: pd.DataFrame) -> pd.DataFrame:
    """Per arm and statistic of runs.csv: the means over the eight wirings and over S's and T's four, the range over
    the eight with the wirings at its ends, and for the threshold the wirings under 60 Hz; the mean over eight is
    checked against readings.csv's all8 run_ column."""
    rows = []
    all8 = readings[(readings.unit == 'all8') & (readings.readout == 'cal')].set_index('kind')
    for kind in REFS + ARMS:
        r = runs[runs.kind == kind].set_index('wiring')
        if sorted(r.index) != sorted(WIRINGS):
            raise ValueError(f'{kind}: wirings {sorted(r.index)}')
        for stat in RUN_STATS:
            if stat not in r.columns:
                continue
            s = pd.to_numeric(r[stat], errors='coerce').reindex(WIRINGS)
            if s.isna().all():
                continue
            ref = f'run_{stat}'
            if ref in all8.columns and kind in all8.index and pd.notna(all8.loc[kind, ref]):
                if abs(s.mean() - all8.loc[kind, ref]) > 1e-6 * max(1.0, abs(s.mean())):
                    raise ValueError(f'{kind} {stat}: mean {s.mean()} != readings {all8.loc[kind, ref]}')
            row = dict(kind=kind, stat=stat, mean8=s.mean(), mean_S4=s[list(S_WIRINGS)].mean(),
                       mean_T4=s.drop(list(S_WIRINGS)).mean(), min8=s.min(), wiring_min=s.idxmin(), max8=s.max(),
                       wiring_max=s.idxmax(), n_wirings=int(s.notna().sum()))
            if stat == 'fixed_thr':
                under = [w for w in WIRINGS if s[w] < THR_LINE]
                row['wirings_under_60'] = ' '.join(under)
                row['n_under_60'] = len(under)
            rows.append(row)
    return pd.DataFrame(rows)


def yardsticks(step0: dict, synapses: pd.DataFrame, runs: pd.DataFrame, registry: dict) -> pd.DataFrame:
    """Step 0's arithmetic and fixed numbers, the per-synapse weight of each arm and the KC data's values (rows:
    group, name, value, note)."""
    rows: List[dict] = []

    def add(group, name, value, note=''):
        rows.append(dict(group=group, name=name, value=float(value), note=note))

    da = step0['dose_arithmetic']
    syn = synapses.set_index('kc_type')
    tot = syn.loc['all']
    for col, key in (('kckc_out', 'kckc_synapses'), ('mbc1_out', 'mbc1_out_synapses'), ('mbc1_in', 'mbc1_in_synapses')):
        if abs(tot[col] - da[key]) > TOL:
            raise ValueError(f'{col} {tot[col]} != step0 {key} {da[key]}')
    g = syn.loc['g']
    n_g = g['n_kcs']
    # The output table (reading 13; MB-C1 presynaptic), as step 0 computes it.
    add('synapses', 'kckc_synapses', tot['kckc_out'], 'KC->KC synapses among the 1732 modelled KCs')
    add('synapses', 'mbc1_out_synapses', tot['mbc1_out'], 'MB-C1 -> modelled KCs (MBC1>KC table), both cells')
    add('synapses', 'mbc1_in_synapses', tot['mbc1_in'], 'modelled KCs -> MB-C1 (KC>MBC1 table), both cells')
    add('synapses', 'apl_out_claws', tot['apl_out_claws'], "APL -> modelled KCs on the wiring's claws")
    add('synapses', 'apl_out_calyx', tot['apl_out_calyx'], 'APL -> modelled KCs in the calyx')
    add('output_table', 'm_equal_synapse', tot['mbc1_out'] / tot['kckc_out'],
        'the dose m at which an MB-C1 output synapse weighs a KC->KC synapse')
    add('output_table', 'mbc1_out_over_apl_claws', tot['mbc1_out'] / tot['apl_out_claws'], '')
    add('output_table', 'mbc1_out_over_apl_calyx', tot['mbc1_out'] / tot['apl_out_calyx'], '')
    add('output_table', 'mbc1_out_per_g_kc', g['mbc1_out'] / n_g, '')
    add('output_table', 'apl_claws_per_g_kc', g['apl_out_claws'] / n_g, '')
    add('output_table', 'apl_calyx_per_g_kc', g['apl_out_calyx'] / n_g, '')
    add('output_table', 'g_share_mbc1_vs_apl_claws', g['mbc1_out'] / (g['mbc1_out'] + g['apl_out_claws']),
        "MB-C1's share of the inhibitory synapses onto gamma KCs, against APL's on the claws")
    add('output_table', 'g_share_mbc1_vs_apl_calyx', g['mbc1_out'] / (g['mbc1_out'] + g['apl_out_calyx']),
        "the same against APL's in the calyx")
    for t in ('g', "a'b'", 'ab'):
        add('output_table', f'mbc1_out_share_{t}', syn.loc[t, 'mbc1_out'] / tot['mbc1_out'], '')
        add('output_table', f'mbc1_in_share_{t}', syn.loc[t, 'mbc1_in'] / tot['mbc1_in'], '')
        add('output_table', f'kc_share_{t}', syn.loc[t, 'n_kcs'] / tot['n_kcs'], '')
    # The same arithmetic from the input table (KC -> MB-C1), which the entry's numbers came from.
    add('input_table', 'm_equal_synapse', tot['mbc1_in'] / tot['kckc_out'], "the entry's 'about 0.16'")
    add('input_table', 'mbc1_over_apl_claws', tot['mbc1_in'] / tot['apl_out_claws'],
        "the entry's 'about a third of APL's onto KCs'")
    add('input_table', 'mbc1_per_g_kc', g['mbc1_in'] / n_g, "the entry's 'about 6.5'")
    add('input_table', 'g_share_mbc1_vs_apl_claws', g['mbc1_in'] / (g['mbc1_in'] + g['apl_out_claws']),
        "the entry's 'about 0.3'")
    for k in ('m_equal_synapse', 'g_share_mbc1_vs_apl_claws', 'g_share_mbc1_vs_apl_calyx', 'mbc1_out_per_g_kc'):
        mine = [r['value'] for r in rows if r['group'] == 'output_table' and r['name'] == k][0]
        if abs(mine - da[k]) > 1e-9:
            raise ValueError(f'{k}: {mine} != step0 {da[k]}')
    # Each arm's weight per MB-C1 output synapse against a KC->KC synapse: its MB-C1 scale over the wiring's fitted
    # KC-KC scale (the family's), over m_equal_synapse.
    m_eq = tot['mbc1_out'] / tot['kckc_out']
    fam = runs[runs.kind == 'family'].set_index('wiring')['wKCKC_scale']
    for kind in ARMS:
        r = runs[runs.kind == kind].set_index('wiring')
        m = (r['wMBC1_scale'] / fam.reindex(r.index)).astype(float)
        w = m / m_eq
        if (w - r['mbc1_per_synapse_vs_kckc'].astype(float)).abs().max() > 1e-6:
            raise ValueError(f'{kind}: per-synapse weight differs from runs.csv')
        add('per_synapse', kind, w.mean(), f'dose m {m.mean():.4g} of the fitted KC-KC scale')
        add('per_synapse_input_table', kind, (m / (tot['mbc1_in'] / tot['kckc_out'])).mean(),
            "the same by the entry's arithmetic (the input table's count in place of the output table's)")
    # Step 0's fixed numbers and the cells' activities in the family's runs.
    for name in ('threshold', 'late_ratio'):
        v = pd.Series(step0[name])
        add('step0', f'{name}_min', v.min(), v.idxmin())
        add('step0', f'{name}_max', v.max(), v.idxmax())
    ag = pd.Series(step0['agreement'])
    add('step0', 'agreement_min', ag.min(), ag.idxmin())
    add('step0', 'agreement_max', ag.max(), ag.idxmax())
    add('step0', 'agreement_n_below_gate', (ag < 0.95).sum(), f'of {len(ag)}')
    gi = step0['gamma_input']
    for cond in ('on', 'off'):
        keys = [k for k in gi if k.endswith('|' + cond)]
        pm = [v for k in keys for v in gi[k]['pulse_mean_hz'].values()]
        pk = [v for k in keys for v in gi[k]['peak_hz'].values()]
        tot_ratio = [gi[k]['mbc1_over_kckc_total'] for k in keys]
        med_ratio = [gi[k]['mbc1_over_kckc_median'] for k in keys]
        add('step0', f'cell_pulse_mean_hz_{cond}_min', min(pm), 'either cell, eight wirings')
        add('step0', f'cell_pulse_mean_hz_{cond}_max', max(pm), 'either cell, eight wirings')
        add('step0', f'cell_peak_hz_{cond}_min', min(pk), '')
        add('step0', f'cell_peak_hz_{cond}_max', max(pk), '')
        add('step0', f'g_mbc1_over_kckc_total_{cond}_min', min(tot_ratio), "gamma KCs' MB-C1 input at m 1 over KC-KC")
        add('step0', f'g_mbc1_over_kckc_total_{cond}_max', max(tot_ratio), '')
        add('step0', f'g_mbc1_over_kckc_median_{cond}_min', min(med_ratio), 'the median gamma KC')
        add('step0', f'g_mbc1_over_kckc_median_{cond}_max', max(med_ratio), '')
    # The KC data's values the note sets beside the model's.
    for k in REGISTRY:
        e = registry[k]
        add('kc_data', k, e['value'] if isinstance(e, dict) else e, 'registry')
        ci = e.get('ci') if isinstance(e, dict) else None
        if ci:
            add('kc_data', f'{k}.lo', ci[0], 'registry, 90 % interval')
            add('kc_data', f'{k}.hi', ci[1], 'registry, 90 % interval')
    return pd.DataFrame(rows)


def _cal(readings: pd.DataFrame, kind: str, unit: str, col: str = 'all_pattern') -> float:
    r = readings[(readings.kind == kind) & (readings.unit == unit) & (readings.readout == 'cal')]
    if len(r) != 1:
        raise ValueError(f'{kind} {unit}: {len(r)} calibrated rows')
    return float(r[col].iloc[0])


def differ_extra(readings: pd.DataFrame, runs: pd.DataFrame, verdict: dict,
                 pairs: Sequence[Tuple[str, str]] = EXTRA_PAIRS) -> pd.DataFrame:
    """The "differ" rule (experiments/E131_mbc1_cells.py's `differ`, E129's) on comparisons the entry did not fix,
    after checking that it reproduces verdict.json's comparisons from the same readings.csv and runs.csv.

    Returns:
        one row per (a, b): group, a, b, verdict, the wirings on which each is lower (counts and names), whether each
        is lower on all ten targets, the targets on which each is lower, the corner arms, a minus b in the calibrated
        pattern part over all eight wirings, on S and T against all recordings, its range over the ten targets, and
        wiring by wiring (d_<wiring>).
    """
    from experiments import E131_mbc1_cells as E131
    for key, c in verdict.get('comparisons', {}).items():
        a, b = key.split('|')
        r = E131.differ(readings, runs, a, b)
        bad = {k: (r.get(k), c.get(k)) for k in DIFFER_KEYS if r.get(k) != c.get(k)}
        if bad:
            raise ValueError(f'{key}: the rule on readings.csv differs from verdict.json: {bad}')
    rows = []
    for a, b in pairs:
        r = E131.differ(readings, runs, a, b)
        if r['verdict'] == 'incomplete':
            raise ValueError(f'{a}|{b}: incomplete')
        pw = r['d_per_wiring']
        per_t = [v for t in r['d_per_target'].values() for v in t.values()]
        row = dict(group='arm_vs_family' if b == 'family' else 'form_vs_B25', a=a, b=b, verdict=r['verdict'],
                   n_wirings_a_lower=r['n_wirings_a_lower'], n_wirings_b_lower=r['n_wirings_b_lower'],
                   wirings_a_lower=' '.join(w for w in WIRINGS if pw[w] < 0),
                   wirings_b_lower=' '.join(w for w in WIRINGS if pw[w] > 0),
                   a_lower_on_all_targets=r['a_lower_on_all_targets'],
                   b_lower_on_all_targets=r['b_lower_on_all_targets'], n_targets=len(per_t),
                   n_targets_a_lower=sum(v < 0 for v in per_t), n_targets_b_lower=sum(v > 0 for v in per_t),
                   corner_arms=' '.join(r['corner_arms']),
                   a_minus_b_all8=_cal(readings, a, 'all8') - _cal(readings, b, 'all8'),
                   a_minus_b_S4=r['d_per_target']['all']['S4'], a_minus_b_T4=r['d_per_target']['all']['T4'],
                   a_minus_b_target_min=min(per_t), a_minus_b_target_max=max(per_t))
        row.update({f'd_{w}': pw[w] for w in WIRINGS})
        rows.append(row)
    return pd.DataFrame(rows)


def draw(margins: pd.DataFrame, per_wiring: pd.DataFrame, path: Path, within: float) -> None:
    """Draws the two panels to `path` (PNG)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw=dict(width_ratios=[1, 1]))
    kinds = ['E110_A'] + list(CLASSED)
    for i, k in enumerate(kinds):
        for j, pop in enumerate(POPULATIONS):
            st = POP_STYLE[pop]
            m = margins[(margins.kind == k) & (margins.population == pop)]
            x = i + (j - 0.5) * 0.3
            halves = m[m.target != 'all']
            a.scatter(np.full(len(halves), x), halves.minus_family, s=20, facecolors='none', edgecolors=st['color'],
                      marker=st['marker'], linewidths=1.1, zorder=3)
            al = m[m.target == 'all']
            a.scatter([x], al.minus_family, s=64, color=st['color'], marker=st['marker'], edgecolors='white',
                      linewidths=1.4, zorder=4, label=st['label'] if i == 0 else None)
    for ax in (a, b):
        ax.axhline(0, color='#52514e', lw=1, zorder=1)
        ax.axhline(within, color='#52514e', lw=1, ls=':', zorder=1)
    a.axvline(2.5, color='#c9c7c1', lw=0.8, zorder=0)
    a.text(len(kinds) - 0.45, within + 0.0003, f'+{within}', ha='right', va='bottom', fontsize=8.5, color='#52514e')
    a.set_xticks(range(len(kinds)), [LABELS[k] for k in kinds], fontsize=8)
    a.set_xlim(-0.6, len(kinds) - 0.4)
    a.set_ylabel("pattern part minus the family's")
    a.set_title('The classes: "costs" is over the dotted line on S and T against all recordings\n'
                '(filled: all recordings; open: the four halves)', fontsize=10)
    a.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=8.5, frameon=False)

    rng = np.random.default_rng(0)
    d_all = per_wiring[per_wiring.target == 'all']
    kinds_b = list(CLASSED)
    for i, k in enumerate(kinds_b):
        for pop in POPULATIONS:
            st = POP_STYLE[pop]
            d = d_all[(d_all.kind == k) & (d_all.population == pop)]
            x = i + rng.uniform(-0.15, 0.15, len(d))
            b.scatter(x, d.minus_family, s=30, color=st['color'], marker=st['marker'], edgecolors='white',
                      linewidths=0.9, zorder=3)
    b.axvline(1.5, color='#c9c7c1', lw=0.8, zorder=0)
    b.set_xticks(range(len(kinds_b)), [LABELS[k] for k in kinds_b], fontsize=8)
    b.set_xlim(-0.6, len(kinds_b) - 0.4)
    b.set_ylabel("pattern part minus the family's, per wiring")
    b.set_title('Wiring by wiring, against all recordings', fontsize=10)
    for ax in (a, b):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e4e3df', lw=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='folder for the tables and the figure (default: the '
                                                            'results folder, which holds the inputs)')
    args = ap.parse_args(argv)
    res = H.results_dir(EXPERIMENT)
    out = args.out or res
    out.mkdir(parents=True, exist_ok=True)
    paths = dict(verdict=res / 'verdict.json', readings=res / 'readings.csv', runs=res / 'runs.csv',
                 step0=res / 'step0.json', synapses=res / 'step0_synapses.csv',
                 registry=H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json')
    verdict = json.loads(paths['verdict'].read_text())
    readings = pd.read_csv(paths['readings'])
    runs = pd.read_csv(paths['runs'])
    step0 = json.loads(paths['step0'].read_text())
    synapses = pd.read_csv(paths['synapses'])
    registry = json.loads(paths['registry'].read_text())
    within = verdict.get('rules', {}).get('within', 0.003)

    M = class_margins(verdict)
    W, WS = wiring_summary(verdict, within)
    C = cost_split(verdict)
    R = run_ranges(runs, readings)
    Y = yardsticks(step0, synapses, runs, registry)
    D = differ_extra(readings, runs, verdict)
    M.to_csv(out / 'writeup_class_margins.csv', index=False)
    W.merge(WS, on=['kind', 'target']).to_csv(out / 'writeup_wiring_summary.csv', index=False)
    C.to_csv(out / 'writeup_cost_split.csv', index=False)
    R.to_csv(out / 'writeup_run_ranges.csv', index=False)
    Y.to_csv(out / 'writeup_yardsticks.csv', index=False)
    D.to_csv(out / 'writeup_differ_extra.csv', index=False)
    fig = out / 'E131_pattern_part_classes.png'
    draw(M, W, fig, within)
    prov = dict(script=H.script_path(__file__), **H.code_commits(), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
                verdict_generated=verdict.get('generated'),
                inputs={str(p.relative_to(H.shared_dir())): _sha256(p) for p in paths.values()},
                outputs=[fig.name, 'writeup_class_margins.csv', 'writeup_wiring_summary.csv', 'writeup_cost_split.csv',
                         'writeup_run_ranges.csv', 'writeup_yardsticks.csv', 'writeup_differ_extra.csv'],
                differ_extra_pairs=[f'{a}|{b}' for a, b in EXTRA_PAIRS],
                checks=['class targets: d_family = arm - family', 'runs.csv means over 8 = readings.csv all8 run_',
                        'per-synapse weights = runs.csv mbc1_per_synapse_vs_kckc',
                        'step-0 synapse totals = step0.json dose_arithmetic',
                        "the differ rule on readings.csv and runs.csv = verdict.json's comparisons"])
    (out / 'provenance_writeup_tables.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    pd.set_option('display.width', 200)
    print(M[M.target == 'all'].pivot(index='kind', columns='population', values='minus_family').round(4))
    print(WS[WS.target == 'all'].round(4).to_string(index=False))
    print(C[C.target == 'all'].round(4).to_string(index=False))
    print(Y.to_string(index=False))
    print(D.drop(columns=[f'd_{w}' for w in WIRINGS]).round(4).to_string(index=False))
    print(f'wrote {fig}')


if __name__ == '__main__':
    main()
