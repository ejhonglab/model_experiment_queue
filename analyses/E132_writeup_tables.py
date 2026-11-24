"""Write-up tables and figure for E132 (recurrent excitation among KCs with a pooled inhibition), from E132's results
folder, the saved runs its score reads and the data-target registry; no simulation.

The numbers E132's note gives that the score's tables do not give directly:
- `run_ranges`: per kind and statistic of runs.csv, the mean over the eight wirings, over S's and T's four, the range
  over the eight (with the wirings at the ends), and the corners and runs outside the bands; the family's late shares
  from late_family.csv (E121's rebuild), which runs.csv does not carry;
- `mean_effects`: per kind and population (S4, T4, all8), the mean APL effect on raw counts and under the readout, the
  raw-count loss and the calibrated pattern part, each also minus the family's and minus the kind's own dose 0 (U0 for
  the pooled arms, KC-KC at zero, E110's part A, for the arm alone): outcome (b)'s baseline;
- `differ_extra`: comparisons the entry did not fix, by the experiment's own "differ" rule (`E132.differ` on
  readings.csv and runs.csv, as the score uses it): each arm of steps 1 and 2 against the family;
- `group_errors`: E92's groups of the megamat pairs (ester 15, alcohol 6; E90's `group_masks`) and all 120 pairs, per
  kind, population, readout and condition: the model's mean pair distance, the flies' (all recordings and the halves A
  to D, E108's `dist_targets`, the targets of the classes) and the error (model minus flies), each also relative to the
  mean over the 120 pairs (E92's `relative_distances`, the coordinates of its line);
- `group_compare`: outcome (c)'s reading for the negative doses at -2 (and the mirror's doses beside) against the
  family, U0 and KC-KC at zero, with APL on: per population, readout and target, on the groups' mean distances and on
  the distances relative to the mean pair (E92's line's coordinates), the change of the alcohol and ester groups (a
  minus b), whether the alcohol pairs are closer (a lower and a smaller error against that target) and the ester pairs
  not closer (a not lower), and the shift of E92's line residual (on its S line of the readout, line.json; signed toward
  the control flies' side); a row 'every half' says whether each holds on all four halves;
- `expectations`: the entry's expected ranges beside the measured values, and step 0's slope carried through zero
  beside the arm alone (X2, X5, X10).

The figure (E132_dose_curve.png): the calibrated pattern part against all recordings (S's and T's populations) and the
largest spike count per KC and odor with APL off (mean over the eight wirings), against the dose, for the pooled arm
(U0, UI and UX) and the arm alone (KC-KC at zero, X), with the family and step 0's expectation.

Checks (raise on failure): the mean over the eight wirings of each runs.csv statistic equals readings.csv's all8 run_
column; the population distances recomputed here give readings.csv's mean distances; the raw-count mean APL effect
minus the family's and minus the dose 0's equal verdict.json's step-2 values; the score's own comparisons, recomputed by
the same rule, equal verdict.json's.

Writes writeup_run_ranges.csv, writeup_mean_effects.csv, writeup_differ_extra.csv, writeup_group_errors.csv,
writeup_group_compare.csv, writeup_expectations.csv, E132_dose_curve.png and provenance_writeup_tables.json to the
results folder.

Usage (from the repository root): ./venv/bin/python analyses/E132_writeup_tables.py [--out DIR]
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E132_kckc_excitation'
WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
POPS = {'S4': ('real', 'P1', 'P2', 'M1'), 'T4': ('P3', 'P4', 'M2', 'M3'), 'all8': WIRINGS}
TARGETS = ('all', 'A', 'B', 'C', 'D')
HALVES = ('A', 'B', 'C', 'D')
KINDS = ('family', 'H50', 'E110_A', 'E110_B', 'U', 'E129_P', 'PO50', 'U0', 'UI2', 'UI5', 'UI10', 'UX2', 'UX5', 'UX10',
         'X2', 'X5', 'X10')
ARMS = ('U0', 'UI2', 'UI5', 'UI10', 'UX2', 'UX5', 'UX10', 'X2', 'X5', 'X10')
DOSE0 = {'UI2': 'U0', 'UI5': 'U0', 'UI10': 'U0', 'UX2': 'U0', 'UX5': 'U0', 'UX10': 'U0', 'X2': 'E110_A',
         'X5': 'E110_A', 'X10': 'E110_A'}
DOSE = {'U0': 0, 'UI2': 2, 'UI5': 5, 'UI10': 10, 'UX2': -2, 'UX5': -5, 'UX10': -10, 'X2': -2, 'X5': -5, 'X10': -10,
        'E110_A': 0}
RUN_STATS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
             'max_count_off', 'n_kc_ge20_on', 'n_kc_ge20_off', 'silent_on', 'silent_off', 'late_share_on',
             'late_share_off', 'fixed_thr', 'thr_factor', 'apl_factor', 'wAPLKC_scale', 'wKCKC_scale', 'wU_scale',
             'pooled_factor', 'kckc_in_over_thr_on', 'kckc_in_over_thr_off', 'pooled_in_over_thr_on',
             'pooled_in_over_thr_off', 'rate_on_g', 'rate_off_g', 'ratio_g', 'n_on', 'n_off', 'target_on',
             'target_off')
GROUP_KINDS = ('family', 'H50', 'E110_A', 'U0', 'UI2', 'UI5', 'UI10', 'UX2', 'UX5', 'UX10', 'X2', 'X5', 'X10')
GROUPS = ('ester', 'alcohol', 'all120')
COMPARE = (('UX2', 'family'), ('UX2', 'U0'), ('UX2', 'E110_A'), ('UX2', 'UI2'), ('X2', 'family'), ('X2', 'U0'),
           ('X2', 'E110_A'), ('UI10', 'family'), ('UI10', 'U0'))
EXTRA_PAIRS = tuple((a, 'family') for a in ARMS)
# The entry's expected ranges (experiment_priorities.md, E132, "How it moves the APL-effect loss"; all eight wirings)
ENTRY_EXPECTED = {
    'U0': ('pattern', 0.054, 0.060),
    'UI10': ('pattern', 0.050, 0.056),
    'X2': ('pattern', 0.063, 0.072),
    'X5': ('pattern', 0.063, 0.072),
    'X10': ('pattern', 0.063, 0.072),
}
ENTRY_EXPECTED_ABOVE_DOSE0 = (0.000, 0.008)     # the pooled arm at -2, -5, -10: dose 0 plus 0.000 to 0.008
ENTRY_EXPECTED_MEAN_CAL = (0.17, 0.21)          # the mean APL effect under the readout
ENTRY_EXPECTED_MEAN_RAW_CHANGE = (-0.02, 0.04)  # the raw-count mean APL effect's change from the family's
TOL = 1e-9
COLORS = {'pooled': '#2a78d6', 'alone': '#d55181'}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def run_ranges(runs: pd.DataFrame, rd: pd.DataFrame, late_family: pd.DataFrame) -> pd.DataFrame:
    """Per kind and RUN_STATS column: means over the eight wirings and over S's and T's four, the range with its
    wirings, and the counts of corners and runs outside the bands (statistic 'n_corner', 'n_outside_bands'); the
    family's late shares from late_family.csv."""
    rows = []
    for kind in KINDS:
        sub = runs[runs['kind'] == kind].set_index('wiring')
        if not len(sub):
            continue
        for st in RUN_STATS:
            if st not in sub.columns:
                continue
            s = pd.to_numeric(sub[st], errors='coerce')
            source = 'runs.csv'
            if kind == 'family' and st.startswith('late_share_') and s.isna().all():
                c = st.rsplit('_', 1)[1]
                lf = late_family[late_family['condition'] == c].set_index('wiring')['late_share']
                s, source = lf.reindex(sub.index).astype(float), 'late_family.csv (E121 rebuild)'
            if s.isna().all():
                continue
            have = [w for w in WIRINGS if w in s.index and np.isfinite(s.get(w, np.nan))]
            row = dict(kind=kind, statistic=st, source=source, n_wirings=len(have),
                       mean8=float(s.reindex(list(WIRINGS)).mean()) if len(have) == 8 else np.nan,
                       mean_S4=float(s.reindex(list(POPS['S4'])).mean()),
                       mean_T4=float(s.reindex(list(POPS['T4'])).mean()),
                       min=float(s.min()), max=float(s.max()), wiring_min=str(s.idxmin()), wiring_max=str(s.idxmax()))
            rows.append(row)
            if source == 'runs.csv' and len(have) == 8:
                col = f'run_{st}'
                ref = _v(rd, kind, 'all8', col)
                if np.isfinite(ref) and abs(ref - row['mean8']) > 1e-6 * max(1.0, abs(ref)):
                    raise RuntimeError(f'{kind} {st}: mean over 8 {row["mean8"]} != readings.csv all8 {ref}')
        for st, col in (('n_corner', 'corner'), ('n_outside_bands', 'in_bands')):
            v = sub[col].astype(str).str.lower().isin(['true', '1', '1.0'])
            n = int(v.sum()) if col == 'corner' else int((~v).sum())
            ws = sorted(sub.index[v] if col == 'corner' else sub.index[~v])
            rows.append(dict(kind=kind, statistic=st, source='runs.csv', n_wirings=len(sub), mean8=float(n),
                             wiring_min=' '.join(ws)))
    return pd.DataFrame(rows)


def mean_effects(rd: pd.DataFrame, verdict: dict) -> pd.DataFrame:
    """Per kind and population: mean APL effect (raw, cal), raw-count loss and calibrated pattern part, also minus the
    family's and minus the dose 0's; checked against verdict.json's step-2 values."""
    rows = []
    for kind in KINDS:
        for u in POPS:
            vals = dict(mean_effect_raw=_v(rd, kind, u, 'mean_effect', 'raw'),
                        mean_effect_cal=_v(rd, kind, u, 'mean_effect', 'cal'),
                        loss_raw=_v(rd, kind, u, 'all_loss', 'raw'),
                        pattern_cal=_v(rd, kind, u, 'all_pattern', 'cal'),
                        pattern_raw=_v(rd, kind, u, 'all_pattern', 'raw'))
            if not np.isfinite(vals['pattern_cal']):
                continue
            row = dict(kind=kind, unit=u, dose=DOSE.get(kind, np.nan), dose0=DOSE0.get(kind, ''), **vals)
            for k, v in vals.items():
                row[f'{k}_minus_family'] = v - {'mean_effect_raw': _v(rd, 'family', u, 'mean_effect', 'raw'),
                                                'mean_effect_cal': _v(rd, 'family', u, 'mean_effect', 'cal'),
                                                'loss_raw': _v(rd, 'family', u, 'all_loss', 'raw'),
                                                'pattern_cal': _v(rd, 'family', u, 'all_pattern', 'cal'),
                                                'pattern_raw': _v(rd, 'family', u, 'all_pattern', 'raw')}[k]
                d0 = DOSE0.get(kind)
                if d0:
                    ref = {'mean_effect_raw': ('mean_effect', 'raw'), 'mean_effect_cal': ('mean_effect', 'cal'),
                           'loss_raw': ('all_loss', 'raw'), 'pattern_cal': ('all_pattern', 'cal'),
                           'pattern_raw': ('all_pattern', 'raw')}[k]
                    row[f'{k}_minus_dose0'] = v - _v(rd, d0, u, ref[0], ref[1])
            rows.append(row)
    out = pd.DataFrame(rows)
    for arm, a in verdict.get('step2', {}).get('arms', {}).items():
        for u in POPS:
            r = out[(out['kind'] == arm) & (out['unit'] == u)].iloc[0]
            for mine, theirs in (('mean_effect_raw_minus_family', 'mean_effect_raw_minus_family'),
                                 ('mean_effect_raw_minus_dose0', 'mean_effect_raw_minus_dose0'),
                                 ('loss_raw', 'loss_raw'), ('mean_effect_cal', 'mean_effect_cal')):
                if abs(float(r[mine]) - float(a[theirs][u])) > TOL:
                    raise RuntimeError(f'{arm} {u} {mine}: {r[mine]} != verdict.json {a[theirs][u]}')
    return out


def differ_extra(rd: pd.DataFrame, runs: pd.DataFrame, verdict: dict, E132) -> pd.DataFrame:
    """EXTRA_PAIRS by `E132.differ`; first recomputes the score's own comparisons and checks them against verdict.json."""
    for key, c in verdict['comparisons'].items():
        a, b = key.split('|')
        mine = E132.differ(rd, runs, a, b)
        if mine['verdict'] != c['verdict'] or mine.get('n_wirings_a_lower') != c.get('n_wirings_a_lower'):
            raise RuntimeError(f'{key}: {mine["verdict"]} != verdict.json {c["verdict"]}')
    rows = []
    for a, b in EXTRA_PAIRS:
        d = E132.differ(rd, runs, a, b)
        row = dict(a=a, b=b, verdict=d['verdict'], n_wirings_a_lower=d.get('n_wirings_a_lower'),
                   n_wirings_b_lower=d.get('n_wirings_b_lower'), a_lower_on_all_targets=d.get('a_lower_on_all_targets'),
                   b_lower_on_all_targets=d.get('b_lower_on_all_targets'),
                   corner_arms=' '.join(d.get('corner_arms', [])),
                   wirings_a_lower=' '.join(w for w, v in d.get('d_per_wiring', {}).items() if v < 0))
        for u in ('S4', 'T4', 'all8'):
            row[f'a_minus_b_{u}'] = _v(rd, a, u, 'all_pattern') - _v(rd, b, u, 'all_pattern')
        for t in TARGETS:
            for u in ('S4', 'T4'):
                row[f'a_minus_b_{t}_{u}'] = d.get('d_per_target', {}).get(t, {}).get(u, np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def group_errors(E132, rd: pd.DataFrame) -> pd.DataFrame:
    """E92's groups per kind, population, readout, condition and target (module docstring); the population distances
    are checked against readings.csv's mean distances."""
    from experiments import E108_paa_extrapolation as E108
    from experiments import E90_model_stages as E90
    from experiments import E83_cross_wiring_fit as E83
    outs, _ = E132.collect()
    flies = E108.dist_targets()
    rows = []
    for kind in GROUP_KINDS:
        if not all((kind, w) in outs for w in WIRINGS):
            continue
        for ro in ('raw', 'cal'):
            dd = {w: E132._dist(outs[(kind, w)], ro) for w in WIRINGS}
            for u, ws in POPS.items():
                d = E83.population([dd[w] for w in ws])
                for c in ('on', 'off'):
                    m = E90._cs(d[c])
                    ref = _v(rd, kind, u, f'mean_{c}', ro)
                    if abs(float(m.mean()) - ref) > 1e-9:
                        raise RuntimeError(f'{kind} {u} {ro} {c}: mean distance {m.mean()} != readings.csv {ref}')
                    gm = E90.group_masks(m.index)
                    mv = m.to_numpy(dtype=float)
                    mmu = float(np.nanmean(mv))
                    for t in TARGETS:
                        f = E90._cs(flies[t][c]).reindex(m.index)
                        fv = f.to_numpy(dtype=float)
                        fmu = float(np.nanmean(fv))
                        for g in GROUPS:
                            mask = np.ones(len(mv), dtype=bool) if g == 'all120' else np.asarray(gm[g])
                            mg, fg = float(np.nanmean(mv[mask])), float(np.nanmean(fv[mask]))
                            rows.append(dict(kind=kind, unit=u, readout=ro, condition=c, target=t, group=g,
                                             n_pairs=int(mask.sum()), model=mg, flies=fg, error=mg - fg,
                                             model_rel=mg - mmu, flies_rel=fg - fmu, error_rel=(mg - mmu) - (fg - fmu)))
    E132._save_dist_cache()
    return pd.DataFrame(rows)


def group_compare(ge: pd.DataFrame, line: dict) -> pd.DataFrame:
    """Outcome (c)'s reading, APL on (module docstring), in two readings: on the groups' mean distances ('abs_') and on
    the distances relative to the mean over the 120 pairs ('rel_', the coordinates of E92's line). Per reading:
    the alcohol and ester groups' change a minus b, 'alcohol_closer' (a lower than b and a smaller absolute error against
    the target), 'ester_not_closer' (a not lower than b) and 'leaves_toward_flies' (both); beside, the shift of E92's
    line residual (its S line of the readout), signed toward the control flies' side."""
    on = ge[ge['condition'] == 'on'].set_index(['kind', 'unit', 'readout', 'target', 'group'])
    rows = []
    for a, b in COMPARE:
        if not ((a, 'S4', 'raw', 'all', 'ester') in on.index and (b, 'S4', 'raw', 'all', 'ester') in on.index):
            continue
        for u in ('S4', 'T4'):
            for ro in ('raw', 'cal'):
                ln = line['lines'][ro]['S']
                side = float(np.sign(ln['flies']['trial_mean_all']['resid']))
                per_t = []
                for t in TARGETS:
                    A = {g: on.loc[(a, u, ro, t, g)] for g in GROUPS}
                    B = {g: on.loc[(b, u, ro, t, g)] for g in GROUPS}
                    r = dict(a=a, b=b, unit=u, readout=ro, target=t)
                    for pre, mcol, ecol in (('abs', 'model', 'error'), ('rel', 'model_rel', 'error_rel')):
                        d_alc = float(A['alcohol'][mcol] - B['alcohol'][mcol])
                        d_est = float(A['ester'][mcol] - B['ester'][mcol])
                        alc = bool(d_alc < 0 and abs(A['alcohol'][ecol]) < abs(B['alcohol'][ecol]))
                        est = bool(d_est >= 0)
                        r.update({f'{pre}_alcohol_error_a': float(A['alcohol'][ecol]),
                                  f'{pre}_alcohol_error_b': float(B['alcohol'][ecol]),
                                  f'{pre}_ester_error_a': float(A['ester'][ecol]),
                                  f'{pre}_ester_error_b': float(B['ester'][ecol]),
                                  f'{pre}_alcohol_a_minus_b': d_alc, f'{pre}_ester_a_minus_b': d_est,
                                  f'{pre}_alcohol_closer': alc, f'{pre}_ester_not_closer': est,
                                  f'{pre}_leaves_toward_flies': bool(alc and est)})
                    d_resid = float((A['ester']['model_rel'] - B['ester']['model_rel'])
                                    - ln['slope'] * (A['alcohol']['model_rel'] - B['alcohol']['model_rel']))
                    r.update(line_resid_shift_toward_flies=d_resid * side, line_max_abs_resid=float(ln['max_abs_resid']))
                    rows.append(r)
                    if t in HALVES:
                        per_t.append(r)
                e = dict(a=a, b=b, unit=u, readout=ro, target='every half')
                for pre in ('abs', 'rel'):
                    for k in ('alcohol_closer', 'ester_not_closer', 'leaves_toward_flies'):
                        e[f'{pre}_{k}'] = all(r[f'{pre}_{k}'] for r in per_t)
                    for k in ('alcohol_a_minus_b', 'ester_a_minus_b'):
                        e[f'{pre}_{k}'] = per_t[0][f'{pre}_{k}']
                e.update(line_resid_shift_toward_flies=per_t[0]['line_resid_shift_toward_flies'],
                         line_max_abs_resid=per_t[0]['line_max_abs_resid'])
                rows.append(e)
    return pd.DataFrame(rows)


def expectations(rd: pd.DataFrame, slopes: pd.DataFrame) -> pd.DataFrame:
    """The entry's expected ranges and step 0's slope carried through zero, beside the measured values."""
    rows = []
    for kind, (stat, lo, hi) in ENTRY_EXPECTED.items():
        v = _v(rd, kind, 'all8', 'all_pattern')
        rows.append(dict(kind=kind, statistic='pattern part, readout, all eight wirings', source='the entry',
                         expected_lo=lo, expected_hi=hi, measured=v, inside=bool(lo <= v <= hi)))
    for kind in ('UX2', 'UX5', 'UX10'):
        v = _v(rd, kind, 'all8', 'all_pattern') - _v(rd, 'U0', 'all8', 'all_pattern')
        lo, hi = ENTRY_EXPECTED_ABOVE_DOSE0
        rows.append(dict(kind=kind, statistic='pattern part minus U0, readout, all eight wirings', source='the entry',
                         expected_lo=lo, expected_hi=hi, measured=v, inside=bool(lo <= v <= hi)))
    for kind in ARMS:
        v = _v(rd, kind, 'all8', 'mean_effect', 'cal')
        lo, hi = ENTRY_EXPECTED_MEAN_CAL
        rows.append(dict(kind=kind, statistic='mean APL effect, readout, all eight wirings', source='the entry',
                         expected_lo=lo, expected_hi=hi, measured=v, inside=bool(lo <= v <= hi)))
        v = _v(rd, kind, 'all8', 'mean_effect', 'raw') - _v(rd, 'family', 'all8', 'mean_effect', 'raw')
        lo, hi = ENTRY_EXPECTED_MEAN_RAW_CHANGE
        rows.append(dict(kind=kind, statistic='mean APL effect on raw counts minus the family\'s, all eight wirings',
                         source='the entry', expected_lo=lo, expected_hi=hi, measured=v,
                         inside=bool(lo <= v <= hi)))
    sl = slopes[(slopes['readout'] == 'cal') & (slopes['reading'] == 'pattern')].set_index('unit')
    for kind, dcol in (('X2', 'at_m2'), ('X5', 'at_m5'), ('X10', 'at_m10')):
        for u in POPS:
            rows.append(dict(kind=kind, statistic=f'pattern part, readout, {u}',
                             source='step 0: the slope from half to zero carried through zero',
                             expected_lo=float(sl.loc[u, dcol]), expected_hi=float(sl.loc[u, dcol]),
                             measured=_v(rd, kind, u, 'all_pattern'),
                             measured_minus_expected=_v(rd, kind, u, 'all_pattern') - float(sl.loc[u, dcol])))
    return pd.DataFrame(rows)


def draw(rd: pd.DataFrame, slopes: pd.DataFrame, runs: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pooled = ('UX10', 'UX5', 'UX2', 'U0', 'UI2', 'UI5', 'UI10')
    alone = ('X10', 'X5', 'X2', 'E110_A')
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    for arm_set, key, label in ((pooled, 'pooled', 'with the pooled neuron (pooled hold)'),
                                (alone, 'alone', 'the matrix alone (threshold hold)')):
        x = [DOSE[k] for k in arm_set]
        for u, ls, mfc, mk in (('S4', '-', None, 'o'), ('T4', '--', 'white', 's')):
            y = [_v(rd, k, u, 'all_pattern') for k in arm_set]
            ax.plot(x, y, ls=ls, color=COLORS[key], lw=2, marker=mk, ms=7,
                    mfc=COLORS[key] if mfc is None else mfc, mec=COLORS[key], zorder=3,
                    label=f'{label}, {"S" if u == "S4" else "T"}')
    sl = slopes[(slopes['readout'] == 'cal') & (slopes['reading'] == 'pattern')].set_index('unit')
    ax.plot([0, -2, -5, -10], [sl.loc['S4', 'at_zero'], sl.loc['S4', 'at_m2'], sl.loc['S4', 'at_m5'],
                               sl.loc['S4', 'at_m10']], ls=':', color='#57606a', lw=1.5,
            label="step 0's slope carried through zero, S")
    for u, ls in (('S4', '-'), ('T4', '--')):
        ax.axhline(_v(rd, 'family', u, 'all_pattern'), color='#8c959f', ls=ls, lw=1, zorder=1)
    ax.text(10.4, _v(rd, 'family', 'S4', 'all_pattern'), 'the family', va='center', ha='left', fontsize=8,
            color='#57606a')
    ax.axvspan(-10.8, -3.5, color='#eef0f2', zorder=0)
    ax.text(-7.1, 0.057, 'corners, runaway\n(flagged)', ha='center', va='center', fontsize=8, color='#57606a')
    ax.set_yscale('log')
    ax.set_yticks([0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25])
    ax.set_yticklabels(['0.05', '0.06', '0.08', '0.10', '0.15', '0.20', '0.25'])
    ax.minorticks_off()
    ax.set_xlim(-10.8, 12.6)
    ax.set_xticks([-10, -5, -2, 0, 2, 5, 10])
    ax.set_xlabel('KC–KC dose (negative: excitation)')
    ax.set_ylabel('pattern part of the APL-effect loss\n(calibrated readout, all recordings; lower is better)')
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    ax = axes[1]
    for arm_set, key, label in ((pooled, 'pooled', 'with the pooled neuron'), (alone, 'alone', 'the matrix alone')):
        x = [DOSE[k] for k in arm_set]
        y = [float(pd.to_numeric(runs.loc[runs['kind'] == k, 'max_count_off']).mean()) for k in arm_set]
        ax.plot(x, y, color=COLORS[key], lw=2, marker='o', ms=7, zorder=3, label=label)
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    fam = float(pd.to_numeric(runs.loc[runs['kind'] == 'family', 'max_count_off']).mean())
    ax.axhline(fam, color='#8c959f', lw=1)
    ax.axhline(2 * fam, color='#8c959f', lw=1, ls=':')
    ax.text(10.4, fam, 'the family', va='center', ha='left', fontsize=8, color='#57606a')
    ax.text(10.4, 2 * fam, "2× the family\n(stop rule)", va='center', ha='left', fontsize=8, color='#57606a')
    ax.axvspan(-10.8, -3.5, color='#eef0f2', zorder=0)
    ax.set_yscale('log')
    ax.set_xlim(-10.8, 12.6)
    ax.set_xticks([-10, -5, -2, 0, 2, 5, 10])
    ax.set_xlabel('KC–KC dose (negative: excitation)')
    ax.set_ylabel('largest spike count per KC and odor, APL off\n(mean over the eight wirings)')
    for a in axes:
        a.spines[['top', 'right']].set_visible(False)
        a.grid(axis='y', color='#e4e3df', lw=0.6, zorder=0)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='folder for the tables and the figure (default: the '
                                                            'results folder, which holds the inputs)')
    args = ap.parse_args(argv)
    from experiments import E132_kckc_excitation as E132
    res = H.results_dir(EXPERIMENT)
    out = args.out or res
    out.mkdir(parents=True, exist_ok=True)
    line_path = H.results_dir('E92_steps_like_for_like') / 'line.json'
    paths = dict(verdict=res / 'verdict.json', readings=res / 'readings.csv', runs=res / 'runs.csv',
                 slopes=res / 'step0_slopes.csv', late_family=res / 'late_family.csv', e92_line=line_path)
    verdict = json.loads(paths['verdict'].read_text())
    rd = pd.read_csv(paths['readings'])
    runs = pd.read_csv(paths['runs'])
    runs['corner'] = runs['corner'].astype(str).str.lower().isin(['true', '1', '1.0'])
    runs['in_bands'] = runs['in_bands'].astype(str).str.lower().isin(['true', '1', '1.0'])
    runs['thr_below_60'] = pd.to_numeric(runs['fixed_thr'], errors='coerce') < 60.0
    slopes = pd.read_csv(paths['slopes'])
    late_family = pd.read_csv(paths['late_family'])
    line = json.loads(line_path.read_text())

    R = run_ranges(runs, rd, late_family)
    M = mean_effects(rd, verdict)
    D = differ_extra(rd, runs, verdict, E132)
    X = expectations(rd, slopes)
    G = group_errors(E132, rd)
    C = group_compare(G, line)
    tables = dict(writeup_run_ranges=R, writeup_mean_effects=M, writeup_differ_extra=D, writeup_group_errors=G,
                  writeup_group_compare=C, writeup_expectations=X)
    for name, t in tables.items():
        t.to_csv(out / f'{name}.csv', index=False)
    fig = out / 'E132_dose_curve.png'
    draw(rd, slopes, runs, fig)
    prov = dict(script=H.script_path(__file__), **H.code_commits(), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
                verdict_generated=verdict.get('generated'),
                inputs={str(p.relative_to(H.shared_dir())): _sha256(p) for p in paths.values()},
                also_read='the saved runs E132_kckc_excitation.collect() loads (references and steps 1 and 2), through '
                          'E132._dist (its distances cache); E108.dist_targets (the flies)',
                outputs=[fig.name] + [f'{n}.csv' for n in tables],
                compare_pairs=[f'{a}|{b}' for a, b in COMPARE],
                differ_extra_pairs=[f'{a}|{b}' for a, b in EXTRA_PAIRS],
                checks=['runs.csv means over 8 = readings.csv all8 run_', 'population mean distances = readings.csv '
                        'mean_on / mean_off', "raw mean APL effect minus family / dose 0, raw loss and mean effect "
                        "under the readout = verdict.json step2", "the differ rule on readings.csv and runs.csv = "
                        "verdict.json's comparisons"])
    (out / 'provenance_writeup_tables.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 40)
    print(M[M.unit != 'all8'].round(4).to_string(index=False))
    print(D.drop(columns=[c for c in D.columns if c.startswith('a_minus_b_') and c.count('_') > 4]).round(4)
          .to_string(index=False))
    print(X.round(4).to_string(index=False))
    print(C[(C.target == 'every half')].round(4).to_string(index=False))
    print(f'wrote {fig}')


if __name__ == '__main__':
    main()
