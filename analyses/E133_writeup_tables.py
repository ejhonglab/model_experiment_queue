"""Write-up tables for E133 (the differences between the hemibrain and the male CNS's left calyx moved one at a time),
steps 0 to 2, from E133's results folder, the model's tables and the data-target registry; no simulation.

The numbers E133's note gives that the score's own tables (step0_table.md, step1_table.md, step2_table.md) do not give
directly, each checked against verdict1.json, verdict2.json, step0.json, d2.csv, builder_sizes.csv or the step-0
tables (a difference raises):
- `step0_groups`: the calibrated pattern error by pair group (the 42 pairs that hold 2-but, paa or B-myr, the 78 that
  hold none, the 39 with exactly one of the three, the 3 among them, and the 81 of the counts as first written) for
  E130's three fixed-scale runs and three fits and the family's eight wirings (range); the left minus the hemibrain on
  every group and target; and the split of the squared gap over the groups (the pattern part squared is the
  pair-count-weighted mean of its groups' squares, so the split is exact);
- `step1_redraws`: per calyx, the seven re-draws' mean, range and SD of the calibrated pattern part and its groups,
  the own wiring's value and its rank among the eight; L - H and R - H on every target; the same on raw counts; the
  squared gap's split;
- `step2_arms`: per arm and host, d (calibrated pattern part, all recordings) with its four wirings, SD and standard
  error, its halves, its groups, its amplitude and shape parts and its raw-count d, recomputed here from readings1.csv
  and readings2.csv; the classes re-read by E133_score2's own rule functions on d2.csv; G, the local gap and its groups,
  and the derived readings (ALL's margin over G/2 and its leave-one-out means, the part of the local gap no arm moves,
  how the arms add, the right's and the left-minus-right gaps on the four wirings for step 3); and the classes re-read
  with the local gap in G's place (a reading beside the rule);
- `realism`: the four bands' inputs (runs1.csv, runs2.csv) for the three calyces' own runs, the re-draws, the
  transplants, B and ALL on the left host and the runs lowest on the calibrated pattern part, beside the KC data's
  values (the registry's noise_* keys, the bands of kcq/apl_effect.py); and per arm and host the gamma KCs' APL-off /
  APL-on ratio, the held threshold and APL->KC scale over the host wiring's own, the largest spike counts over the own
  run's, and the mean APL effect under the readout (the runs outside the flies' value +-0.02 listed);
- `b_factors`: arm B's factor per glomerulus (E133_builders.build('B') on the model's tables, `light_tables`: APL's
  counted synapses, which the runner's normalized weights differ from by one factor per vector, and no factor depends
  on), the boutons per glomerulus and each glomerulus's mean APL->bouton and bouton->APL weight per bouton as a multiple
  of the calyx's mean, on the three calyces; whether the glomerulus has ORN input in the model (the 36 of
  `apl16_orn_deltas`); the share of APL's summed weight at the boutons on the glomeruli with input, before and after B
  on each host;
- `planner_check`: every number of the three gate readings (the Fable planners of 2026-11-24 06:41, 07:43 and 08:57,
  chat.md and E133's Status line) that these tables hold, beside the script's value at the precision written, and
  whether it reproduces (a number that does not is reported, not raised).

Writes writeup_step0_groups.csv, writeup_step0_gap_split.csv, writeup_step1_redraws.csv, writeup_step1_gaps.csv,
writeup_step2_arms.csv, writeup_step2_derived.csv, writeup_step2_local_yardstick.csv, writeup_realism.csv, writeup_realism_arms.csv,
writeup_mean_effect_outside.csv, writeup_sizes.csv, writeup_b_factors.csv, writeup_b_shares.csv,
writeup_planner_check.csv, writeup_tables.md (all of them as markdown) and provenance_writeup_tables.json to the
results folder.

Usage (from the repository root; light, no simulation):
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 19 ./venv/bin/python analyses/E133_writeup_tables.py [--out DIR]
"""
import argparse
import hashlib
import json
import math
import re
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E133_calyx_transplants'
SCRIPT = 'analyses/E133_writeup_tables.py'
HEMI, LEFT, RIGHT = 'hemibrain-argmax', 'malecns-left', 'malecns-right'
CALYCES = (HEMI, LEFT, RIGHT)
SHORT = {HEMI: 'hemibrain (argmax)', LEFT: 'male CNS left', RIGHT: 'male CNS right'}
HOSTS = (LEFT, HEMI)
DONOR = {LEFT: HEMI, HEMI: LEFT}
ARMS = ('C', 'S', 'A', 'B', 'K', 'ALL')
ARM_NAME = {'C': 'C, claws per KC by type', 'S': "S, the glomeruli's claw shares", 'A': "A, APL's synapses at the KCs",
            'B': "B, APL's synapses at the boutons", 'K': "K, the matrix's type-by-type shares",
            'ALL': 'ALL, the five together'}
SEEDS1 = tuple(range(1, 8))
EIGHT = (0,) + SEEDS1
HOST_SEEDS = (0, 1, 2, 3)
TARGETS = ('all', 'A', 'B', 'C', 'D')
HALVES = TARGETS[1:]
GROUPS = ('p120', 'p42', 'p78', 'p39', 'p81', 'p3')
NPAIRS = {'p120': 120, 'p42': 42, 'p78': 78, 'p39': 39, 'p81': 81, 'p3': 3}
GROUP_NAME = {'p120': 'all 120', 'p42': 'the 42 that hold 2-but, paa or B-myr', 'p78': 'the 78 that hold none',
              'p39': 'the 39 with exactly one of the three', 'p81': 'the other 81 (the counts as first written)',
              'p3': 'the 3 among the three'}
UNITS0 = ('S1:hemibrain-argmax', 'S1:malecns-left', 'S1:malecns-right', 'fit:hemibrain-argmax', 'fit:malecns-left',
          'fit:malecns-right')
FAMILY = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')
B_GLOMERULI = ('DM1', 'DM4', 'VM7d', 'DP1m', 'DL5', 'VP1m', 'DC4')
TOL = 1e-12
MEAN_EFFECT_TOL = 0.02     # the entry's expectation: the mean APL effect under the readout within +-0.02 of the flies'


def _sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _col(t: str, g: str) -> str:
    """readings column of the pattern error of group g against target t."""
    return f'{t}_pattern' if g == 'p120' else f'{t}_pattern_{g}'


def _mean(x) -> float:
    x = [float(v) for v in x]
    return math.fsum(x) / len(x)


def _sd(x) -> float:
    return float(np.std(np.asarray(x, dtype=float), ddof=1))


class Checks:
    """Every comparison with a saved number; `eq` raises on a difference."""

    def __init__(self):
        self.n = 0

    def eq(self, what: str, a, b, tol: float = TOL) -> None:
        self.n += 1
        if isinstance(a, str) or isinstance(b, str) or isinstance(a, (bool, np.bool_)):
            if a != b:
                raise AssertionError(f'{what}: {a!r} != {b!r}')
            return
        a, b = float(a), float(b)
        if not (np.isfinite(a) and np.isfinite(b)) or abs(a - b) > tol * max(1.0, abs(a), abs(b)):
            raise AssertionError(f'{what}: {a!r} != {b!r}')


# ---------------------------------------------------------------------------
# inputs

def load(res: Path) -> dict:
    j = lambda n: json.loads((res / n).read_text())
    c = lambda n: pd.read_csv(res / n, low_memory=False)
    reg = json.loads((H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json').read_text())
    return dict(rd1=c('readings1.csv'), runs1=c('runs1.csv'), rd2=c('readings2.csv'), runs2=c('runs2.csv'),
                d2=c('d2.csv'), v1=j('verdict1.json'), v2=j('verdict2.json'), s0=j('step0.json'),
                sp=c('step0_patterns.csv'), sg=c('step0_glomeruli.csv'), sgg=c('step0_glomeruli_by_glomerulus.csv'),
                sizes=c('builder_sizes.csv'), reg=reg)


def _rd(rd: pd.DataFrame, calyx: str, seed: int, readout: str = 'cal', arm: Optional[str] = None) -> pd.Series:
    m = (rd['calyx'] == calyx) & (rd['seed'] == seed) & (rd['readout'] == readout)
    if arm is not None:
        m &= rd['arm'] == arm
    r = rd[m]
    if len(r) != 1:
        raise AssertionError(f'{len(r)} rows for {calyx} seed {seed} {readout} {arm}')
    return r.iloc[0]


def _run(runs: pd.DataFrame, calyx: str, seed: int, arm: Optional[str] = None) -> pd.Series:
    m = (runs['calyx'] == calyx) & (runs['seed'] == seed)
    if arm is not None:
        m &= runs['arm'] == arm
    r = runs[m]
    if len(r) != 1:
        raise AssertionError(f'{len(r)} run rows for {calyx} seed {seed} {arm}')
    return r.iloc[0]


# ---------------------------------------------------------------------------
# step 0

def step0_groups(I: dict, ck: Checks):
    sp = I['sp'][I['sp']['readout'] == 'cal'].set_index('unit')
    rows = []
    for u in UNITS0:
        r = dict(run=u)
        for g in GROUPS:
            r[g] = float(sp.loc[u, _col('all', g)])
        for k in ('amplitude', 'shape', 'sd_ratio', 'pearson_effect'):
            r[k] = float(sp.loc[u, f'all_{k}'])
        rows.append(r)
    fam = [f'family:{w}' for w in FAMILY]
    for agg in ('min', 'max'):
        r = dict(run=f'family eight wirings, {agg}')
        for g in GROUPS:
            r[g] = float(getattr(sp.loc[fam, _col('all', g)], agg)())
        rows.append(r)
    groups = pd.DataFrame(rows)
    for t in TARGETS:                      # E130's left minus hemibrain, against step0.json
        for g in ('p42', 'p78', 'p39', 'p81'):
            d = float(sp.loc['S1:malecns-left', _col(t, g)] - sp.loc['S1:hemibrain-argmax', _col(t, g)])
            ck.eq(f'step0 L-H {t} {g}', d, I['s0']['fixed_scales_left_minus_hemi'][t][f'left_minus_hemi_{g}'])
    ck.eq('step0 n_pairs', sum(I['s0']['n_pairs'][g] for g in ('p42', 'p78')), 120)
    for g in GROUPS:
        ck.eq(f'step0 n_pairs {g}', I['s0']['n_pairs'][g], NPAIRS[g])
    rows = []
    for kind in ('S1', 'fit'):
        L, Hh = f'{kind}:malecns-left', f'{kind}:hemibrain-argmax'
        for t in TARGETS:
            r = dict(runs=kind, target=t)
            for g in GROUPS:
                r[f'L_minus_H_{g}'] = float(sp.loc[L, _col(t, g)] - sp.loc[Hh, _col(t, g)])
            tot = float(sp.loc[L, _col(t, 'p120')] ** 2 - sp.loc[Hh, _col(t, 'p120')] ** 2)
            for g in ('p42', 'p78', 'p39', 'p81', 'p3'):
                r[f'share_sq_{g}'] = NPAIRS[g] / 120 * float(sp.loc[L, _col(t, g)] ** 2 - sp.loc[Hh, _col(t, g)] ** 2) / tot
            ck.eq(f'step0 split {kind} {t} 42+78', r['share_sq_p42'] + r['share_sq_p78'], 1.0, 1e-9)
            ck.eq(f'step0 split {kind} {t} 39+81', r['share_sq_p39'] + r['share_sq_p81'], 1.0, 1e-9)
            ck.eq(f'step0 split {kind} {t} 39+3', r['share_sq_p39'] + r['share_sq_p3'], r['share_sq_p42'], 1e-9)
            rows.append(r)
    split = pd.DataFrame(rows)
    return groups, split


# ---------------------------------------------------------------------------
# step 1

def step1(I: dict, ck: Checks):
    rd = I['rd1']
    v1 = I['v1']
    rows = []
    for c in CALYCES:
        for ro in ('cal', 'raw'):
            for t in TARGETS:
                for g in GROUPS:
                    col = _col(t, g)
                    own = float(_rd(rd, c, 0, ro)[col])
                    red = [float(_rd(rd, c, s, ro)[col]) for s in SEEDS1]
                    eight = [own] + red
                    rows.append(dict(calyx=c, readout=ro, target=t, group=g, own=own, mean7=_mean(red),
                                     min7=min(red), max7=max(red), sd7=_sd(red), mean8=_mean(eight),
                                     own_rank_of_8=int(1 + sum(x < own for x in red))))
    red = pd.DataFrame(rows)

    def m(c, ro, t, g, k='mean7'):
        x = red[(red['calyx'] == c) & (red['readout'] == ro) & (red['target'] == t) & (red['group'] == g)]
        return float(x[k].iloc[0])
    gaps = []
    for ro in ('cal', 'raw'):
        for t in TARGETS:
            r = dict(readout=ro, target=t)
            hv = [float(_rd(rd, HEMI, s, ro)[_col(t, 'p120')]) for s in SEEDS1]
            lv = [float(_rd(rd, LEFT, s, ro)[_col(t, 'p120')]) for s in SEEDS1]
            rv = [float(_rd(rd, RIGHT, s, ro)[_col(t, 'p120')]) for s in SEEDS1]
            med = float(np.median(hv))
            r.update(L=m(LEFT, ro, t, 'p120'), H=m(HEMI, ro, t, 'p120'), R=m(RIGHT, ro, t, 'p120'),
                     hemi_median=med, left_above_hemi_median=int(sum(x > med for x in lv)),
                     right_above_hemi_median=int(sum(x > med for x in rv)))
            r['L_minus_H'], r['R_minus_H'], r['L_minus_R'] = r['L'] - r['H'], r['R'] - r['H'], r['L'] - r['R']
            for g in ('p42', 'p78', 'p3', 'p39'):
                r[f'L_minus_H_{g}'] = m(LEFT, ro, t, g) - m(HEMI, ro, t, g)
                r[f'R_minus_H_{g}'] = m(RIGHT, ro, t, g) - m(HEMI, ro, t, g)
            # the squared gap's split over the groups: on the re-draws' means (squares of the means), and exactly
            # (the mean over the re-draws of each run's square)
            tot = r['L'] ** 2 - r['H'] ** 2
            sq = lambda c, g: _mean([float(_rd(rd, c, s, ro)[_col(t, g)]) ** 2 for s in SEEDS1])
            tot_exact = sq(LEFT, 'p120') - sq(HEMI, 'p120')
            for g in ('p42', 'p78', 'p3', 'p39'):
                r[f'share_sq_{g}'] = NPAIRS[g] / 120 * (m(LEFT, ro, t, g) ** 2 - m(HEMI, ro, t, g) ** 2) / tot
                r[f'share_sq_exact_{g}'] = NPAIRS[g] / 120 * (sq(LEFT, g) - sq(HEMI, g)) / tot_exact
            ck.eq(f'step1 exact split {ro} {t}', r['share_sq_exact_p42'] + r['share_sq_exact_p78'], 1.0, 1e-9)
            if ro == 'cal':
                v = v1[t]
                ck.eq(f'step1 L {t}', r['L'], v['L'])
                ck.eq(f'step1 H {t}', r['H'], v['H'])
                ck.eq(f'step1 gap {t}', r['L_minus_H'], v['gap'])
                ck.eq(f'step1 median {t}', med, v['hemi_median'])
                ck.eq(f'step1 above {t}', r['left_above_hemi_median'], v['left_above_hemi_median'])
                ck.eq(f'step1 R {t}', r['R'], v['right']['R'])
                ck.eq(f'step1 R-H {t}', r['R_minus_H'], v['right']['R_minus_H'])
                ck.eq(f'step1 R-L {t}', -r['L_minus_R'], v['right']['R_minus_L'])
                ck.eq(f'step1 right above {t}', r['right_above_hemi_median'], v['right']['right_above_hemi_median'])
                ck.eq(f'step1 class {t}', v['cls'], 'kept')
                for c in CALYCES:
                    ck.eq(f'step1 own {t} {c}', m(c, 'cal', t, 'p120', 'own'), v['own'][c])
                    for g in ('p42', 'p78', 'p39', 'p81', 'p3'):
                        ck.eq(f'step1 group {t} {g} {c}', m(c, 'cal', t, g), v['groups'][g][c])
            gaps.append(r)
    ck.eq('step1 no-class runs', len(v1['no_class_runs']), 0)
    return red, pd.DataFrame(gaps)


# ---------------------------------------------------------------------------
# step 2

D_STATS = {'pattern': 'all_pattern', 'amplitude': 'all_amplitude', 'shape': 'all_shape', 'mean_effect': 'mean_effect',
           **{g: f'all_pattern_{g}' for g in ('p42', 'p78', 'p39', 'p81', 'p3')},
           **{f'half_{t}': f'{t}_pattern' for t in HALVES}}


def step2(I: dict, ck: Checks):
    from experiments import E133_score2 as S
    rd1, rd2, d2, v2 = I['rd1'], I['rd2'], I['d2'], I['v2']
    # G and the local gap, from step 1's readings
    G = {}
    for t in TARGETS:
        col = _col(t, 'p120')
        L8 = _mean([_rd(rd1, LEFT, s)[col] for s in EIGHT])
        H8 = _mean([_rd(rd1, HEMI, s)[col] for s in EIGHT])
        L4 = _mean([_rd(rd1, LEFT, s)[col] for s in HOST_SEEDS])
        H4 = _mean([_rd(rd1, HEMI, s)[col] for s in HOST_SEEDS])
        G[t] = dict(G=L8 - H8, local=L4 - H4, L8=L8, H8=H8, L4=L4, H4=H4)
        for k in G[t]:
            ck.eq(f'G {t} {k}', G[t][k], v2['G'][t][k])
    g = G['all']['G']
    ck.eq('G/2', g / 2, v2['thresholds']['G_half'])
    ck.eq('G/4', g / 4, v2['thresholds']['G_quarter'])
    rows = []
    for arm in ARMS:
        for h in HOSTS:
            r = dict(arm=arm, host=h)
            for ro in ('cal', 'raw'):
                for name, col in D_STATS.items():
                    if ro == 'raw' and name not in ('pattern', 'amplitude', 'shape', 'mean_effect'):
                        continue
                    dw = [float(_rd(rd2, h, s, ro, arm)[col]) - float(_rd(rd1, h, s, ro)[col]) for s in HOST_SEEDS]
                    key = name if ro == 'cal' else f'raw_{name}'
                    r[f'd_{key}'] = _mean(dw)
                    x = d2[(d2['arm'] == arm) & (d2['host'] == h) & (d2['readout'] == ro) & (d2['stat'] == col)]
                    if len(x) == 1:
                        ck.eq(f'd2 {arm} {h} {ro} {col}', r[f'd_{key}'], x['d'].iloc[0])
                        for i, s in enumerate(HOST_SEEDS):
                            ck.eq(f'd2 {arm} {h} {ro} {col} seed{s}', dw[i], x[f'd_seed{s}'].iloc[0])
                    if name == 'pattern' and ro == 'cal':
                        for i, s in enumerate(HOST_SEEDS):
                            r[f'd_seed{s}'] = dw[i]
                        r['sd_wirings'] = _sd(dw)
                        r['se'] = _sd(dw) / 2.0
                        r['n_toward'] = int(sum((x < 0) if h == LEFT else (x > 0) for x in dw))
                        r['d_min'], r['d_max'] = min(dw), max(dw)
                        for k, s in enumerate(HOST_SEEDS):
                            r[f'd_without_seed{s}'] = _mean([dw[i] for i in range(4) if i != k])
                        x7 = d2[(d2['arm'] == arm) & (d2['host'] == h) & (d2['readout'] == 'cal')
                                & (d2['stat'] == 'all_pattern')]
                        ck.eq(f'd2 sd {arm} {h}', r['sd_wirings'], x7['sd_wirings'].iloc[0], 1e-9)
                        r['sd_redraws'] = float(x7['sd_redraws'].iloc[0])
            va = v2['arms'][arm]
            hv = va['hosts'][h]
            ck.eq(f'v2 d {arm} {h}', r['d_pattern'], hv['d'])
            ck.eq(f'v2 n_toward {arm} {h}', r['n_toward'], hv['n_toward'])
            for t in HALVES:
                ck.eq(f'v2 half {arm} {h} {t}', r[f'd_half_{t}'], hv['d_halves'][t])
            for gg in ('p42', 'p78', 'p39', 'p81', 'p3'):
                ck.eq(f'v2 group {arm} {h} {gg}', r[f'd_{gg}'], va['groups'][gg][h])
            ck.eq(f'v2 raw {arm} {h}', r['d_raw_pattern'], va['raw_pattern'][h])
            r['every_half'] = bool(hv['every_half'])
            r['cls'], r['side'] = va['cls'], va['side'] or ''
            rows.append(r)
    arms = pd.DataFrame(rows)
    # the classes re-read by the score's own rule functions on d2.csv
    rr = S.arm_readings(d2, g)
    for arm in ARMS:
        ck.eq(f'class {arm}', rr[arm]['rule_class'], v2['arms'][arm]['cls'])
        ck.eq(f'side {arm}', str(rr[arm]['side']), str(v2['arms'][arm]['side']))
        ck.eq(f'directional {arm}', 'cls_directional' in rr[arm], False)
    ck.eq('step3 selected', ' '.join(S.step3_selection(rr_cls(rr))), ' '.join(v2['step3']['selected']))
    ck.eq('runs in', v2['n_runs_in'], 48)
    ck.eq('no-class runs', len(v2['no_class_runs']), 0)
    ck.eq('realism flags', len(v2['realism_flags']), 0)
    # derived readings
    A = arms.set_index(['arm', 'host'])
    d = lambda arm, h, k='d_pattern': float(A.loc[(arm, h), k])
    der = {}
    der['G'], der['G_half'], der['G_quarter'] = g, g / 2, g / 4
    der['local'] = G['all']['local']
    for t in HALVES:
        der[f'G_{t}'] = G[t]['G']
        der[f'local_{t}'] = G[t]['local']
    for gg in ('p42', 'p78', 'p3', 'p39'):
        der[f'local_{gg}'] = (_mean([_rd(rd1, LEFT, s)[_col('all', gg)] for s in HOST_SEEDS])
                              - _mean([_rd(rd1, HEMI, s)[_col('all', gg)] for s in HOST_SEEDS]))
    for c in CALYCES:
        for gg in ('p120', 'p42', 'p78', 'p3'):
            der[f'own4_{c}_{gg}'] = _mean([_rd(rd1, c, s)[_col('all', gg)] for s in HOST_SEEDS])
    for gg in ('p120', 'p42', 'p78'):
        der[f'R4_minus_H4_{gg}'] = der[f'own4_{RIGHT}_{gg}'] - der[f'own4_{HEMI}_{gg}']
        der[f'L4_minus_R4_{gg}'] = der[f'own4_{LEFT}_{gg}'] - der[f'own4_{RIGHT}_{gg}']
    der['ALL_left_over_G_half'] = -d('ALL', LEFT) - g / 2
    der['ALL_hemi_over_G_half'] = d('ALL', HEMI) - g / 2
    der['ALL_left_without_seed0'] = d('ALL', LEFT, 'd_without_seed0')
    der['ALL_hemi_without_seed2'] = d('ALL', HEMI, 'd_without_seed2')
    der['unmoved_left'] = der['local'] + d('ALL', LEFT)
    der['unmoved_hemi'] = der['local'] - d('ALL', HEMI)
    der['S_left_under_G_quarter'] = g / 4 + d('S', LEFT)
    der['B_hemi_under_G_quarter'] = g / 4 - d('B', HEMI)
    der['B_plus_S_hemi'] = d('B', HEMI) + d('S', HEMI)
    der['B_left'], der['ALL_left'], der['ALL_hemi'] = d('B', LEFT), d('ALL', LEFT), d('ALL', HEMI)
    der['B_left_shape'] = d('B', LEFT, 'd_shape')
    der['B_left_p3'] = d('B', LEFT, 'd_p3')
    der['S_left_se'] = d('S', LEFT, 'se')
    der['left_fit_minus_G_half'] = float(I['sp'][(I['sp']['readout'] == 'cal')
                                                  & (I['sp']['unit'] == 'fit:malecns-left')]['all_pattern'].iloc[0]) - g / 2
    lines = [g / 2, g / 4, 0.004, 0.003]
    der['nearest_d_to_a_line'] = min(abs(abs(x) - ln) for x in arms['d_pattern'] for ln in lines)
    derived = pd.DataFrame([dict(quantity=k, value=v) for k, v in der.items()])
    # the classes with the local gap on step 2's four wirings in G's place (a reading beside the rule, which reads G)
    rl = S.arm_readings(d2, G['all']['local'])
    yard = pd.DataFrame([dict(arm=arm, cls_G=rr[arm]['rule_class'], side_G=rr[arm]['side'] or '',
                              cls_local=rl[arm]['rule_class'], side_local=rl[arm]['side'] or '',
                              **{f'part_local_{h}': bool(rl[arm]['hosts'][h]['part']) for h in HOSTS},
                              **{f'full_local_{h}': bool(rl[arm]['hosts'][h]['full']) for h in HOSTS})
                         for arm in ARMS])
    return arms, derived, der, yard


def rr_cls(rr: dict) -> dict:
    """E133_score2.arm_readings' output with 'cls' set to the rule class (no run is without a class here)."""
    return {a: dict(r, cls=r['rule_class']) for a, r in rr.items()}


# ---------------------------------------------------------------------------
# realism

R_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv')


def _true(x: pd.Series) -> pd.Series:
    """A column of booleans as written to CSV (True / False, or their strings)."""
    return x.astype(str).str.strip().str.lower().eq('true')


def _reg(I: dict, k: str) -> float:
    return float(I['reg'][k]['value'])


def realism(I: dict, ck: Checks):
    from kcq import apl_effect as AE
    runs1, runs2, rd1, rd2 = I['runs1'], I['runs2'], I['rd1'], I['rd2']
    pat1 = rd1[rd1['readout'] == 'cal'].set_index('tag')['all_pattern']
    pat2 = rd2[rd2['readout'] == 'cal'].set_index('tag')['all_pattern']
    allruns = pd.concat([runs1.assign(arm='own', step=1), runs2.assign(step=2)], ignore_index=True)
    allruns['pattern'] = [float(pat1[t]) if s == 1 else float(pat2[t]) for t, s in zip(allruns['tag'], allruns['step'])]
    ck.eq('runs', len(allruns), 72)

    def row(label: str, sub: pd.DataFrame) -> dict:
        r = dict(runs=label, n=len(sub), tags=' '.join(sub['tag']) if len(sub) <= 4 else '')
        for k in ('pattern',) + R_COLS:
            x = pd.to_numeric(sub[k], errors='coerce')
            r[f'{k}_min'], r[f'{k}_max'] = float(x.min()), float(x.max())
        r['n_in_bands'] = int(_true(sub['in_bands']).sum())
        r['n_flagged'] = int(sub['realism_flag'].fillna('').astype(str).str.strip().ne('').sum())
        r['n_corner'] = int(_true(sub['corner']).sum())
        r['n_no_class'] = int(_true(sub['no_class']).sum())
        return r
    rows = []
    for c in CALYCES:
        rows.append(row(f'{SHORT[c]}, own wiring (fixed scales)', allruns[(allruns['step'] == 1) & (allruns['calyx'] == c)
                                                                         & (allruns['seed'] == 0)]))
    for c in CALYCES:
        rows.append(row(f'{SHORT[c]}, seven re-draws', allruns[(allruns['step'] == 1) & (allruns['calyx'] == c)
                                                             & (allruns['seed'] > 0)]))
    for h in HOSTS:
        rows.append(row(f'the 24 transplants on the {SHORT[h]}', allruns[(allruns['step'] == 2) & (allruns['calyx'] == h)]))
    for arm in ('B', 'ALL'):
        rows.append(row(f'{arm} on the male CNS left', allruns[(allruns['step'] == 2) & (allruns['calyx'] == LEFT)
                                                              & (allruns['arm'] == arm)]))
    best = allruns.loc[[allruns['pattern'].idxmin()]]
    rows.append(row(f'lowest calibrated pattern part of the 72 ({best["tag"].iloc[0]})', best))
    bl = allruns[allruns['calyx'] == LEFT]
    bl = bl.loc[[bl['pattern'].idxmin()]]
    rows.append(row(f'lowest on the male CNS left ({bl["tag"].iloc[0]})', bl))
    ck.eq('all 72 in bands', int(_true(allruns['in_bands']).sum()), 72)
    ck.eq('no corner', int(_true(allruns['corner']).sum()), 0)
    ck.eq('no no-class', int(_true(allruns['no_class']).sum()), 0)
    ck.eq('no realism flag', int(allruns['realism_flag'].fillna('').astype(str).str.strip().ne('').sum()), 0)
    data = dict(runs='KC data', n=np.nan, tags='', pattern_min=_reg(I, 'apl_effect_pattern_floor'),
                rate_on_min=_reg(I, 'noise_ctrl_frac'), rate_off_min=_reg(I, 'noise_tnt_frac'),
                ratio_min=_reg(I, 'noise_rate_ratio_fraction'), rate_rho_on_min=_reg(I, 'ctrl_rate_rho_ceiling'),
                rate_rho_off_min=_reg(I, 'tnt_rate_rho_ceiling'), opk_on_tv_min=_reg(I, 'noise_ctrl_odors_per_kc_tv_max'),
                opk_off_tv_min=_reg(I, 'noise_tnt_odors_per_kc_tv_max'))
    data['ratio_ci'] = ' to '.join(f'{x:.2f}' for x in I['reg']['noise_rate_ratio_fraction']['ci'])
    data['band_rate_on'] = f'{AE.RATE_ON_BAND[0]} to {AE.RATE_ON_BAND[1]}'
    data['band_ratio'] = f'{AE.RATE_RATIO_BAND[0]} to {AE.RATE_RATIO_BAND[1]}'
    data['band_rate_rho_on'] = f'at least {AE.RATE_RHO_MIN:.2f}'
    rows.append(data)
    main = pd.DataFrame(rows)
    # per arm and host: gamma's ratio, the held parameters over the own run's, the largest counts, the mean APL effect
    me1 = rd1[rd1['readout'] == 'cal'].set_index('tag')['mean_effect']
    me2 = rd2[rd2['readout'] == 'cal'].set_index('tag')['mean_effect']
    arows = []
    for c in CALYCES:
        sub = runs1[runs1['calyx'] == c]
        me = [float(me1[t]) for t in sub['tag']]
        arows.append(dict(arm='own and re-draws (step 1)', host=c, n=len(sub), ratio_g_min=sub['ratio_g'].min(),
                          ratio_g_max=sub['ratio_g'].max(), thr_factor_min=sub['thr_factor'].min(),
                          thr_factor_max=sub['thr_factor'].max(), apl_factor_min=sub['apl_factor'].min(),
                          apl_factor_max=sub['apl_factor'].max(), mean_effect_min=min(me), mean_effect_max=max(me)))
    for arm in ARMS:
        for h in HOSTS:
            sub = runs2[(runs2['arm'] == arm) & (runs2['calyx'] == h)]
            ck.eq(f'runs2 {arm} {h}', len(sub), 4)
            cr = np.maximum(sub['max_count_on'] / sub['max_count_on_own'], sub['max_count_off'] / sub['max_count_off_own'])
            me = [float(me2[t]) for t in sub['tag']]
            arows.append(dict(arm=arm, host=h, n=len(sub), ratio_g_min=sub['ratio_g'].min(), ratio_g_max=sub['ratio_g'].max(),
                              ratio_g_own_min=sub['ratio_g_own'].min(), ratio_g_own_max=sub['ratio_g_own'].max(),
                              thr_vs_own_min=sub['fixed_thr_vs_own'].min(), thr_vs_own_max=sub['fixed_thr_vs_own'].max(),
                              aplkc_vs_own_min=sub['wAPLKC_scale_vs_own'].min(),
                              aplkc_vs_own_max=sub['wAPLKC_scale_vs_own'].max(),
                              aplpn_vs_own_min=sub['wAPLPN_scale_vs_own'].min(),
                              aplpn_vs_own_max=sub['wAPLPN_scale_vs_own'].max(),
                              max_count_vs_own_max=float(cr.max()), rate_rho_on_min=sub['rate_rho_on'].min(),
                              rate_rho_on_max=sub['rate_rho_on'].max(), mean_effect_min=min(me), mean_effect_max=max(me)))
            ck.eq(f'ratio_g_own {arm} {h}', sub['ratio_g_own'].min(),
                  runs1[(runs1['calyx'] == h) & (runs1['seed'].isin(HOST_SEEDS))]['ratio_g'].min(), 1e-9)
    arms = pd.DataFrame(arows)
    # the mean APL effect under the readout against the flies' value +- 0.02
    flies = _reg(I, 'apl_effect_mean')
    out = []
    for step, rd in ((1, rd1), (2, rd2)):
        c = rd[rd['readout'] == 'cal']
        for _, r in c.iterrows():
            dv = float(r['mean_effect']) - flies
            out.append(dict(step=step, tag=r['tag'], calyx=r['calyx'], arm=r.get('arm', 'own') if step == 2 else 'own',
                            seed=int(r['seed']), mean_effect=float(r['mean_effect']), minus_flies=dv,
                            outside=bool(abs(dv) > MEAN_EFFECT_TOL)))
    me = pd.DataFrame(out)
    return main, arms, me, allruns


# ---------------------------------------------------------------------------
# the transplants' sizes, and B's factors per glomerulus

def sizes(I: dict) -> pd.DataFrame:
    s = I['sizes']
    keep = [c for c in s.columns if c.startswith(('size_A_', 'size_B_', 'size_C_off', 'size_C_copied', 'size_S_moved',
                                                  'size_K_within', 'size_C_claws_per_kc'))]
    out = s[(s['host_seed'] == 0) & (s['arm'].isin(('A', 'B', 'C', 'S', 'K')))][['host', 'donor', 'arm'] + keep]
    return out.dropna(axis=1, how='all')


def b_factors(I: dict, ck: Checks):
    from experiments import E133_builders as B
    orn = H.apl16_orn_deltas()
    with_input = set(str(g) for g in orn.index)
    ck.eq('glomeruli with ORN input', len(with_input), 36)
    T = {c: B.light_tables(c, 0) for c in CALYCES}
    out = {h: B.build('B', T[h], T[DONOR[h]], seed=0) for h in HOSTS}
    rows = []
    gl = {c: np.asarray([str(x) for x in T[c]['W'].columns.get_level_values('glomerulus')], dtype=object) for c in CALYCES}
    allg = sorted(set().union(*[set(gl[c]) for c in CALYCES]))
    for gname in allg:
        r = dict(glomerulus=gname, orn_input=gname in with_input)
        for c in CALYCES:
            m = gl[c] == gname
            r[f'boutons_{c}'] = int(m.sum())
            for k, key in (('wAPLPN', 'aplpn'), ('wPNAPL', 'pnapl')):
                w = T[c]['wnorm'][k].to_numpy(dtype=float)
                r[f'{key}_count_{c}'] = float(w[m].sum())
                r[f'{key}_rel_{c}'] = float(w[m].mean() / w.mean()) if m.any() else np.nan
        for h in HOSTS:
            m = gl[h] == gname
            for k, key in (('wAPLPN', 'aplpn'), ('wPNAPL', 'pnapl')):
                old = T[h]['wnorm'][k].to_numpy(dtype=float)[m]
                new = out[h]['wnorm'][k].to_numpy(dtype=float)[m]
                pos = old > 0
                if not m.any():
                    f = np.nan
                elif not pos.any():
                    f = np.nan        # no APL synapse in the host: kept
                else:
                    fs = new[pos] / old[pos]
                    if np.ptp(fs) > 1e-12 * max(1.0, abs(fs).max()):
                        raise AssertionError(f'B {h} {k} {gname}: factors differ over the boutons')
                    f = float(fs[0])
                r[f'{key}_factor_{h}'] = f
        rows.append(r)
    bf = pd.DataFrame(rows)
    # checks against the builder's recorded sizes and step 0's glomerulus tables
    sz = I['sizes']
    for h in HOSTS:
        srow = sz[(sz['host'] == h) & (sz['donor'] == DONOR[h]) & (sz['host_seed'] == 0) & (sz['arm'] == 'B')].iloc[0]
        meta = out[h]['meta']['transplant']['parts']['B']
        for key in ('aplpn', 'pnapl'):
            ck.eq(f'B rescale {h} {key}', meta[f'{key}_rescale'], srow[f'size_B_{key}_rescale'], 1e-9)
            k = 'wAPLPN' if key == 'aplpn' else 'wPNAPL'
            w = T[h]['wnorm'][k].to_numpy(dtype=float)
            fac = out[h]['wnorm'][k].to_numpy(dtype=float)
            elig = np.zeros(len(w), dtype=bool)
            for gname in set(gl[h]):
                m = gl[h] == gname
                if (gl[DONOR[h]] == gname).any() and w[m].sum() > 0:
                    elig |= m
            f = fac[elig] / np.where(w[elig] > 0, w[elig], np.nan)
            f = f[np.isfinite(f)]
            ck.eq(f'B factor min {h} {key}', f.min(), srow[f'size_B_{key}_factor_min'], 1e-9)
            ck.eq(f'B factor max {h} {key}', f.max(), srow[f'size_B_{key}_factor_max'], 1e-9)
            # share of the summed weight on the glomeruli with ORN input, before and after
    shares = []
    for h in HOSTS:
        for k, key in (('wAPLPN', 'aplpn'), ('wPNAPL', 'pnapl')):
            w = T[h]['wnorm'][k].to_numpy(dtype=float)
            wn = out[h]['wnorm'][k].to_numpy(dtype=float)
            m = np.isin(gl[h], list(with_input))
            ck.eq(f'B total kept {h} {key}', wn.sum(), w.sum(), 1e-9)
            fin = bf[bf['orn_input']][f'{key}_factor_{h}'].dropna()
            fall = bf[f'{key}_factor_{h}'].dropna()
            shares.append(dict(host=h, donor=DONOR[h], direction=key, share_with_input_before=float(w[m].sum() / w.sum()),
                               share_with_input_after=float(wn[m].sum() / wn.sum()),
                               factor_min_with_input=float(fin.min()),
                               factor_min_with_input_glomerulus=str(bf.loc[fin.idxmin(), 'glomerulus']),
                               factor_max_with_input=float(fin.max()),
                               factor_max_with_input_glomerulus=str(bf.loc[fin.idxmax(), 'glomerulus']),
                               factor_max_all=float(fall.max()),
                               factor_max_all_glomerulus=str(bf.loc[fall.idxmax(), 'glomerulus'])))
    sgg = I['sgg'].drop_duplicates(['calyx', 'glomerulus']).set_index(['calyx', 'glomerulus'])
    B_ = bf.set_index('glomerulus')
    for (c, gname), r in sgg.iterrows():
        ck.eq(f'step0 boutons {c} {gname}', B_.loc[gname, f'boutons_{c}'], r['n_boutons'])
        ck.eq(f'step0 apl_pn {c} {gname}', B_.loc[gname, f'aplpn_count_{c}'], r['apl_pn'])
        ck.eq(f'step0 pn_apl {c} {gname}', B_.loc[gname, f'pnapl_count_{c}'], r['pn_apl'])
    sg = I['sg'][I['sg']['odor'] == 'paa'].set_index('calyx')       # paa's glomeruli: DM1 alone
    for c in CALYCES:
        ck.eq(f'step0 DM1 APL->bouton per bouton {c}', B_.loc['DM1', f'aplpn_rel_{c}'], sg.loc[c, 'apl_pn_per_bouton_rel'], 1e-9)
        ck.eq(f'step0 DM1 bouton->APL per bouton {c}', B_.loc['DM1', f'pnapl_rel_{c}'], sg.loc[c, 'pn_apl_per_bouton_rel'], 1e-9)
    return bf, pd.DataFrame(shares)


# ---------------------------------------------------------------------------
# the planners' numbers

def _decimals(s: str) -> int:
    m = re.search(r'\.(\d+)', s)
    return len(m.group(1)) if m else 0


def planner_check(g0: pd.DataFrame, sp0: pd.DataFrame, red: pd.DataFrame, gaps: pd.DataFrame, arms: pd.DataFrame,
                  der: dict, rmain: pd.DataFrame, rarms: pd.DataFrame, bf: pd.DataFrame, bsh: pd.DataFrame,
                  I: dict) -> pd.DataFrame:
    """Each number of the three gate readings that the tables hold: (reading, what, written, script value); a percent
    is written with '%'. Reproduced: the script's value rounded to the written precision is the written value."""
    G0 = g0.set_index('run')
    R = red.set_index(['calyx', 'readout', 'target', 'group'])
    GP = gaps.set_index(['readout', 'target'])
    A = arms.set_index(['arm', 'host'])
    BF = bf.set_index('glomerulus')
    BS = bsh.set_index(['host', 'direction'])
    RA = rarms.set_index(['arm', 'host'])
    sp = I['sp'][I['sp']['readout'] == 'cal'].set_index('unit')
    sz = I['sizes']
    S0 = sp0.set_index(['runs', 'target'])
    hl = lambda kind, g: [float(S0.loc[(kind, t), f'L_minus_H_{g}']) for t in HALVES]
    r7 = lambda c, g, k='mean7': float(R.loc[(c, 'cal', 'all', g), k])
    items = []

    def add(reading, what, written, value):
        items.append((reading, what, written, float(value)))
    s0, s1, s2 = 'step 0 (06:41)', 'step 1 (07:43)', 'step 2 (08:57)'
    # step 0's table and its halves
    for u, lab in (('S1:hemibrain-argmax', 'hemibrain'), ('S1:malecns-left', 'left'), ('S1:malecns-right', 'right')):
        for g, w in zip(('p120', 'p42', 'p39', 'p3', 'p78'),
                        {'hemibrain': ('0.0576', '0.051', '0.052', '0.038', '0.061'),
                         'left': ('0.0719', '0.090', '0.082', '0.162', '0.060'),
                         'right': ('0.0625', '0.063', '0.059', '0.101', '0.062')}[lab]):
            add(s0, f'{lab}, fixed scales, {GROUP_NAME[g]}', w, G0.loc[u, g])
    for g, lo, hi in (('p120', '0.047', '0.064'), ('p42', '0.041', '0.061'), ('p39', '0.040', '0.058'),
                      ('p3', '0.021', '0.091'), ('p78', '0.048', '0.070')):
        add(s0, f'family range low, {GROUP_NAME[g]}', lo, G0.loc['family eight wirings, min', g])
        add(s0, f'family range high, {GROUP_NAME[g]}', hi, G0.loc['family eight wirings, max', g])
    for g, lo, hi in (('p42', '+0.032', '+0.041'), ('p39', '+0.024', '+0.032'), ('p78', '-0.006', '+0.006')):
        add(s0, f'left minus hemibrain on the halves, {GROUP_NAME[g]}, low', lo, min(hl('S1', g)))
        add(s0, f'left minus hemibrain on the halves, {GROUP_NAME[g]}, high', hi, max(hl('S1', g)))
    add(s0, 'left fit, the 42', '0.092', G0.loc['fit:malecns-left', 'p42'])
    add(s0, 'hemibrain fit, the 42', '0.048', G0.loc['fit:hemibrain-argmax', 'p42'])
    add(s0, 'left fit, the 78', '0.057', G0.loc['fit:malecns-left', 'p78'])
    add(s0, 'hemibrain fit, the 78', '0.064', G0.loc['fit:hemibrain-argmax', 'p78'])
    add(s0, 'left, the other 81', '0.067', G0.loc['S1:malecns-left', 'p81'])
    add(s0, 'hemibrain, the other 81', '0.060', G0.loc['S1:hemibrain-argmax', 'p81'])
    add(s0, "the 3 pairs' share of the squared gap", 'about a third', 100 * S0.loc[('S1', 'all'), 'share_sq_p3'])
    add(s0, "the 78's share of the squared gap", '-4 %', 100 * S0.loc[('S1', 'all'), 'share_sq_p78'])
    # step 1
    c1 = GP.loc[('cal', 'all')]
    add(s1, 'L', '0.0750', c1['L'])
    add(s1, 'H', '0.0583', c1['H'])
    add(s1, 'L - H', '+0.0167', c1['L_minus_H'])
    add(s1, 'left above the hemibrain median (of 7)', '7', c1['left_above_hemi_median'])
    add(s1, 'L - H on the halves, low', '+0.0126', min(GP.loc[('cal', t), 'L_minus_H'] for t in HALVES))
    add(s1, 'L - H on the halves, high', '+0.0169', max(GP.loc[('cal', t), 'L_minus_H'] for t in HALVES))
    add(s1, 'R - H', '+0.0095', c1['R_minus_H'])
    add(s1, 'R - H on the halves, low', '+0.0048', min(GP.loc[('cal', t), 'R_minus_H'] for t in HALVES))
    add(s1, 'R - H on the halves, high', '+0.0127', max(GP.loc[('cal', t), 'R_minus_H'] for t in HALVES))
    add(s1, 'a re-draw scatters (SD), low', '0.006', min(r7(c, 'p120', 'sd7') for c in (HEMI, LEFT)))
    add(s1, 'a re-draw scatters (SD), high', '0.007', max(r7(c, 'p120', 'sd7') for c in (HEMI, LEFT)))
    tab = {HEMI: (('p120', '0.0583', '0.050', '0.065'), ('p42', '0.048', '0.041', '0.054'), ('p3', '0.044', '0.015', '0.079'),
                  ('p78', '0.063', '0.054', '0.072')),
           LEFT: (('p120', '0.0750', '0.065', '0.086'), ('p42', '0.091', '0.079', '0.108'), ('p3', '0.202', '0.158', '0.249'),
                  ('p78', '0.065', '0.056', '0.072')),
           RIGHT: (('p120', '0.0678', '0.062', '0.074'), ('p42', '0.067', '0.056', '0.076'), ('p3', '0.112', '0.072', '0.159'),
                   ('p78', '0.068', '0.062', '0.073'))}
    for c, rows_ in tab.items():
        for g, mn, lo, hi in rows_:
            add(s1, f'{SHORT[c]} re-draws, {GROUP_NAME[g]}, mean', mn, r7(c, g))
            add(s1, f'{SHORT[c]} re-draws, {GROUP_NAME[g]}, low', lo, r7(c, g, 'min7'))
            add(s1, f'{SHORT[c]} re-draws, {GROUP_NAME[g]}, high', hi, r7(c, g, 'max7'))
    add(s1, "the 42's share of the squared gap (squares of the re-draws' means)", '94 %', 100 * c1['share_sq_p42'])
    add(s1, "the 3's share of the squared gap (squares of the re-draws' means)", '43 %', 100 * c1['share_sq_p3'])
    add(s1, "the 42's share of the squared gap (exact: the re-draws' mean squares)", '94 %', 100 * c1['share_sq_exact_p42'])
    add(s1, "the 3's share of the squared gap (exact: the re-draws' mean squares)", '43 %', 100 * c1['share_sq_exact_p3'])
    add(s1, "the right's own wiring, rank of its eight (1 the lowest)", '2', r7(RIGHT, 'p120', 'own_rank_of_8'))
    add(s1, 'left-only part (L - R)', '+0.007', c1['L_minus_R'])
    rw = GP.loc[('raw', 'all')]
    add(s1, 'raw counts: L - H', '+0.0096', rw['L_minus_H'])
    add(s1, 'raw counts: left above the hemibrain median (of 7)', '6', rw['left_above_hemi_median'])
    rsd = [float(R.loc[(c, 'raw', 'all', 'p120'), 'sd7']) for c in (HEMI, LEFT)]
    add(s1, 'raw counts: re-draw SD, low', '0.009', min(rsd))
    add(s1, 'raw counts: re-draw SD, high', '0.010', max(rsd))
    add(s1, 'raw counts: R - H', '+0.002', rw['R_minus_H'])
    add(s1, 'raw counts: L - H on the 42', '+0.040', rw['L_minus_H_p42'])
    add(s1, 'raw counts: R - H on the 42', '+0.018', rw['R_minus_H_p42'])
    s_moved = sz[(sz['arm'] == 'S') & (sz['host_seed'] == 0) & sz['host'].isin(HOSTS)
                 & (sz['donor'] == sz['host'].map(DONOR))]['size_S_moved_share']
    c_sz = sz[(sz['arm'] == 'C') & (sz['host_seed'] == 0) & sz['host'].isin(HOSTS) & (sz['donor'] == sz['host'].map(DONOR))]
    c_moved = c_sz['size_C_off_share'] + c_sz['size_C_copied_share']
    # written "C and S re-deal claws at random (6.5 % and 13 to 14 % of a host's claws)": read as written (C 6.5 %,
    # S 13 to 14 %) and swapped (S 6.5 %, C 13 to 14 %)
    add(s1, "C's claws switched off or copied, as written (C 6.5 %), the smaller host", '6.5 %', 100 * c_moved.min())
    add(s1, "S's claws moved, as written (S 13 to 14 %), the larger host", '14 %', 100 * s_moved.max())
    add(s1, "S's claws moved, swapped (S 6.5 %), the smaller host", '6.5 %', 100 * s_moved.min())
    add(s1, "S's claws moved, swapped (S 6.5 %), the larger host", '6.5 %', 100 * s_moved.max())
    add(s1, "C's claws switched off or copied, swapped (C 13 to 14 %), low", '13 %', 100 * c_moved.min())
    add(s1, "C's claws switched off or copied, swapped (C 13 to 14 %), high", '14 %', 100 * c_moved.max())
    # step 2
    add(s2, 'G', '+0.0164', der['G'])
    add(s2, 'G/2', '0.0082', der['G_half'])
    add(s2, 'G/4', '0.0041', der['G_quarter'])
    add(s2, 'the local gap', '+0.0143', der['local'])
    add(s2, 'the local gap, the 42', '+0.038', der['local_p42'])
    add(s2, 'the local gap, the 78', '+0.0005', der['local_p78'])
    for arm, (dl, dh, l42, h42, l78, h78) in {
            'C': ('+0.0003', '+0.0015', '-0.0018', '+0.0026', '+0.0014', '+0.0011'),
            'S': ('-0.0041', '+0.0049', '-0.0037', '+0.0070', '-0.0046', '+0.0041'),
            'A': ('-0.0010', '+0.0003', '-0.0024', '+0.0012', '+0.0000', '-0.0001'),
            'B': ('-0.0097', '+0.0036', '-0.0171', '+0.0126', '-0.0047', '-0.0008'),
            'K': ('-0.0003', '+0.0014', '+0.0001', '+0.0015', '-0.0006', '+0.0014'),
            'ALL': ('-0.0098', '+0.0094', '-0.0241', '+0.0185', '-0.0009', '+0.0050')}.items():
        add(s2, f'{arm} d, left host', dl, A.loc[(arm, LEFT), 'd_pattern'])
        add(s2, f'{arm} d, hemibrain host', dh, A.loc[(arm, HEMI), 'd_pattern'])
        add(s2, f'{arm} d on the 42, left host', l42, A.loc[(arm, LEFT), 'd_p42'])
        add(s2, f'{arm} d on the 42, hemibrain host', h42, A.loc[(arm, HEMI), 'd_p42'])
        add(s2, f'{arm} d on the 78, left host', l78, A.loc[(arm, LEFT), 'd_p78'])
        add(s2, f'{arm} d on the 78, hemibrain host', h78, A.loc[(arm, HEMI), 'd_p78'])
    add(s2, 'B on the left: the four wirings, low', '-0.0111', A.loc[('B', LEFT), 'd_min'])
    add(s2, 'B on the left: the four wirings, high', '-0.0091', A.loc[('B', LEFT), 'd_max'])
    add(s2, "B on the left: d of the shape part", '-0.0095', der['B_left_shape'])
    add(s2, 'B on the left: d on the 3 pairs', '-0.064', der['B_left_p3'])
    add(s2, "the left's own error on the 3 pairs (four wirings)", '0.178', der[f'own4_{LEFT}_p3'])
    add(s2, "the hemibrain's own error on the 42 (four wirings)", '0.048', der[f'own4_{HEMI}_p42'])
    add(s2, "the hemibrain's own error on the 78 (four wirings)", '0.060', der[f'own4_{HEMI}_p78'])
    add(s2, "B's hemibrain d under the G/4 line", '0.0005', der['B_hemi_under_G_quarter'])
    add(s2, "S's d on the left", '-0.004093', A.loc[('S', LEFT), 'd_pattern'])
    add(s2, "the 'a part' line, G/4", '0.004098', der['G_quarter'])
    add(s2, "S's d under the line", '0.000005', der['S_left_under_G_quarter'])
    add(s2, "S's standard error", '0.0015', der['S_left_se'])
    add(s2, 'S raw-count d, left host', '+0.0004', A.loc[('S', LEFT), 'd_raw_pattern'])
    add(s2, 'S raw-count d, hemibrain host', '-0.0012', A.loc[('S', HEMI), 'd_raw_pattern'])
    add(s2, 'ALL over G/2, left', '0.0016', der['ALL_left_over_G_half'])
    add(s2, 'ALL over G/2, hemibrain', '0.0012', der['ALL_hemi_over_G_half'])
    add(s2, 'ALL on the left without seed 0', '-0.0081', der['ALL_left_without_seed0'])
    add(s2, 'ALL on the hemibrain without seed 2', '+0.0064', der['ALL_hemi_without_seed2'])
    add(s2, 'the local gap moved by no arm, left host', '0.0044', der['unmoved_left'])
    add(s2, 'the local gap moved by no arm, hemibrain host', '0.0049', der['unmoved_hemi'])
    add(s2, 'B plus S on the hemibrain', '+0.0085', der['B_plus_S_hemi'])
    add(s2, "the right's own gap from the hemibrain, four wirings", '+0.0105', der['R4_minus_H4_p120'])
    add(s2, 'the left minus the right, four wirings, all', '+0.0038', der['L4_minus_R4_p120'])
    add(s2, 'the left minus the right, four wirings, the 42', '+0.020', der['L4_minus_R4_p42'])
    add(s2, 'the left minus the right, four wirings, the 78', '-0.0065', der['L4_minus_R4_p78'])
    add(s2, "step 4's line: the left's fit minus G/2", '0.0633', der['left_fit_minus_G_half'])
    add(s2, 'a quarter of B on the left (the right-host line)', '0.0024', -A.loc[('B', LEFT), 'd_pattern'] / 4)
    add(s2, 'a quarter of ALL on the left (the right-host line)', '0.0025', -A.loc[('ALL', LEFT), 'd_pattern'] / 4)
    # B's factors (the planner's read-only table)
    btab = {'DM1': ('14', '9', '10', '2.36', '5.38', '3.33', '0.38', '0.64'),
            'DM4': ('10', '9', '6', '1.11', '1.65', '1.76', '0.58', '0.71'),
            'VM7d': ('6', '10', '11', '0.48', '0.35', '0.70', '1.19', '1.60'),
            'DP1m': ('16', '10', '9', '0.87', '1.95', '2.04', '0.39', '0.32'),
            'VP1m': ('1', '14', '10', '7.32', '0.60', '0.98', '10.6', '14.8')}
    for gname, (bh, bl, br, mh, ml, mr, fa, fp) in btab.items():
        for c, w in ((HEMI, bh), (LEFT, bl), (RIGHT, br)):
            add(s2, f'{gname} boutons, {SHORT[c]}', w, BF.loc[gname, f'boutons_{c}'])
        for c, w in ((HEMI, mh), (LEFT, ml), (RIGHT, mr)):
            add(s2, f"{gname} APL->bouton per bouton over the calyx's mean, {SHORT[c]}", w, BF.loc[gname, f'aplpn_rel_{c}'])
        add(s2, f"{gname} B's factor on the left host, APL->bouton", fa, BF.loc[gname, f'aplpn_factor_{LEFT}'])
        add(s2, f"{gname} B's factor on the left host, bouton->APL", fp, BF.loc[gname, f'pnapl_factor_{LEFT}'])
    for gname in ('DM1', 'DM4', 'VM7d', 'DP1m'):
        add(s2, f'{gname} has ORN input in the model (1 yes)', '1', float(BF.loc[gname, 'orn_input']))
    add(s2, 'VP1m has ORN input in the model (0 no)', '0', float(BF.loc['VP1m', 'orn_input']))
    bl_ = BS.loc[(LEFT, 'aplpn')]
    bp_ = BS.loc[(LEFT, 'pnapl')]
    add(s2, "among the 36, B's lowest APL->bouton factor on the left (DC4)", '0', bl_['factor_min_with_input'])
    add(s2, "among the 36, B's highest APL->bouton factor on the left", '2.6', bl_['factor_max_with_input'])
    add(s2, "among the 36, B's lowest bouton->APL factor on the left", '0.32', bp_['factor_min_with_input'])
    add(s2, "among the 36, B's highest bouton->APL factor on the left", '1.7', bp_['factor_max_with_input'])
    for gname in ('DM1', 'DP1m', 'DL5'):
        add(s2, f"{gname}'s APL->bouton factor on the left ('about 0.4')", '0.4', BF.loc[gname, f'aplpn_factor_{LEFT}'])
    add(s2, "share of APL->bouton weight on the glomeruli with input, left host, before", '0.86', bl_['share_with_input_before'])
    add(s2, "share of APL->bouton weight on the glomeruli with input, left host, after", '0.69', bl_['share_with_input_after'])
    bh_ = BS.loc[(HEMI, 'aplpn')]
    add(s2, "share of APL->bouton weight on the glomeruli with input, hemibrain host, before", '0.865', bh_['share_with_input_before'])
    add(s2, "share of APL->bouton weight on the glomeruli with input, hemibrain host, after", '0.905', bh_['share_with_input_after'])
    add(s2, 'B on the left: APL->KC over the own run, low', '1.36', RA.loc[('B', LEFT), 'aplkc_vs_own_min'])
    add(s2, 'B on the left: APL->KC over the own run, high', '1.40', RA.loc[('B', LEFT), 'aplkc_vs_own_max'])
    # realism of step 2
    RM = rmain.set_index('runs')
    t2 = [f'the 24 transplants on the {SHORT[h]}' for h in HOSTS]
    add(s2, 'response rate with APL on (48 runs), low', '0.0668', min(RM.loc[t, 'rate_on_min'] for t in t2))
    add(s2, 'response rate with APL on (48 runs), high', '0.0668', max(RM.loc[t, 'rate_on_max'] for t in t2))
    add(s2, 'response rate with APL off (48 runs), low', '0.1468', min(RM.loc[t, 'rate_off_min'] for t in t2))
    add(s2, 'response rate with APL off (48 runs), high', '0.1468', max(RM.loc[t, 'rate_off_max'] for t in t2))
    add(s2, 'their ratio (48 runs), low', '2.20', min(RM.loc[t, 'ratio_min'] for t in t2))
    add(s2, 'their ratio (48 runs), high', '2.20', max(RM.loc[t, 'ratio_max'] for t in t2))
    add(s2, "the flies' ratio", '2.28', RM.loc['KC data', 'ratio_min'])
    add(s2, 'response rate ρ with APL on (48 runs), low', '0.48', min(RM.loc[t, 'rate_rho_on_min'] for t in t2))
    add(s2, 'response rate ρ with APL on (48 runs), high', '0.80', max(RM.loc[t, 'rate_rho_on_max'] for t in t2))
    add(s2, "the largest spike count over the host wiring's own, highest of the 48", '1.63',
        rarms[rarms['arm'].isin(ARMS)]['max_count_vs_own_max'].max())
    add(s2, "γ's ratio, ALL on the left, low", '3.1', RA.loc[('ALL', LEFT), 'ratio_g_min'])
    add(s2, "γ's ratio, ALL on the left, high", '3.4', RA.loc[('ALL', LEFT), 'ratio_g_max'])
    add(s2, "γ's ratio, the left's own four wirings, low", '2.2', RA.loc[('ALL', LEFT), 'ratio_g_own_min'])
    add(s2, "γ's ratio, the left's own four wirings, high", '2.4', RA.loc[('ALL', LEFT), 'ratio_g_own_max'])
    add(s2, "γ's ratio, the hemibrain's own four wirings, low", '3.5', RA.loc[('ALL', HEMI), 'ratio_g_own_min'])
    add(s2, "γ's ratio, the hemibrain's own four wirings, high", '4.0', RA.loc[('ALL', HEMI), 'ratio_g_own_max'])
    add(s2, "γ's ratio, C on the left, low", '2.4', RA.loc[('C', LEFT), 'ratio_g_min'])
    add(s2, "γ's ratio, C on the left, high", '2.6', RA.loc[('C', LEFT), 'ratio_g_max'])
    a_sz = sz[(sz['arm'] == 'A') & (sz['host'] == LEFT) & (sz['donor'] == HEMI) & (sz['host_seed'] == 0)].iloc[0]
    add(s2, "A's left KCs left without APL->KC weight", '24', a_sz['size_A_aplkc_kcs_to_zero'])
    rows = []
    for reading, what, written, value in items:
        if written == 'about a third':
            # an approximate reading: within a tenth of a third
            num, d = 100 / 3, None
            rep, trunc = bool(abs(value - num) <= num / 10), None
        else:
            d = _decimals(written.replace('%', '').strip())
            num = float(written.replace('%', '').replace('+', '').strip())
            rep = bool(abs(value - num) <= 0.5 * 10 ** (-d) + 1e-12)
            # whether the written value is the script's cut off (toward zero) at the written decimals
            trunc = bool(abs(math.trunc(round(value * 10 ** d, 9)) / 10 ** d - num) <= 1e-12)
        rows.append(dict(reading=reading, what=what, written=written, script=value, reproduced=rep,
                         truncated_matches=trunc))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# markdown

def _fmt(x, sp='.4f') -> str:
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '–'
    return format(float(x), sp).replace('-', '−')


def _rng(a, b, sp='.4f') -> str:
    return _fmt(a, sp) if format(float(a), sp) == format(float(b), sp) else f'{_fmt(a, sp)} to {_fmt(b, sp)}'


def markdown(g0, sp0, red, gaps, arms, der, yard, rmain, rarms, me, bf, bsh, pc) -> str:
    L = ['# E133 write-up tables (generated by analyses/E133_writeup_tables.py; no simulation)', '']
    L += ['## Step 0: the pattern error by pair group (calibrated readout, all recordings)', '',
          '| run | 120 | 42 with 2-but, paa or B-myr | 39 with one | 3 among the three | 78 with none | 81 (as first written) |',
          '|---|---:|---:|---:|---:|---:|---:|']
    for _, r in g0.iterrows():
        L.append(f"| {r['run']} | {_fmt(r['p120'])} | {_fmt(r['p42'], '.3f')} | {_fmt(r['p39'], '.3f')} | "
                 f"{_fmt(r['p3'], '.3f')} | {_fmt(r['p78'], '.3f')} | {_fmt(r['p81'], '.3f')} |")
    L += ['', '### The left minus the hemibrain, and the squared gap by pair group', '',
          '| runs | target | 120 | 42 | 39 | 3 | 78 | share of squared gap: 42 / 39 / 3 / 78 |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in sp0.iterrows():
        L.append(f"| {r['runs']} | {r['target']} | {_fmt(r['L_minus_H_p120'], '+.4f')} | {_fmt(r['L_minus_H_p42'], '+.3f')} | "
                 f"{_fmt(r['L_minus_H_p39'], '+.3f')} | {_fmt(r['L_minus_H_p3'], '+.3f')} | {_fmt(r['L_minus_H_p78'], '+.3f')} | "
                 f"{_fmt(r['share_sq_p42'], '.2f')} / {_fmt(r['share_sq_p39'], '.2f')} / {_fmt(r['share_sq_p3'], '.2f')} / "
                 f"{_fmt(r['share_sq_p78'], '.2f')} |")
    L += ['', '## Step 1: the seven re-draws (calibrated readout, all recordings)', '',
          '| calyx | group | own wiring | re-draws, mean | range | SD | own rank of 8 |', '|---|---|---:|---:|---:|---:|---:|']
    for c in CALYCES:
        for g in ('p120', 'p42', 'p3', 'p78'):
            r = red[(red['calyx'] == c) & (red['readout'] == 'cal') & (red['target'] == 'all') & (red['group'] == g)].iloc[0]
            L.append(f"| {SHORT[c]} | {GROUP_NAME[g]} | {_fmt(r['own'])} | {_fmt(r['mean7'])} | "
                     f"{_rng(r['min7'], r['max7'], '.3f')} | {_fmt(r['sd7'])} | {int(r['own_rank_of_8'])} |")
    L += ['', '| readout | target | L | H | R | L − H | left above H median | R − H | L − R | L − H on 42 / 78 / 3 | '
          'share of squared gap 42 / 3 (exact 42 / 3) |', '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in gaps.iterrows():
        L.append(f"| {r['readout']} | {r['target']} | {_fmt(r['L'])} | {_fmt(r['H'])} | {_fmt(r['R'])} | "
                 f"{_fmt(r['L_minus_H'], '+.4f')} | {int(r['left_above_hemi_median'])} | {_fmt(r['R_minus_H'], '+.4f')} | "
                 f"{_fmt(r['L_minus_R'], '+.4f')} | {_fmt(r['L_minus_H_p42'], '+.3f')} / {_fmt(r['L_minus_H_p78'], '+.3f')} / "
                 f"{_fmt(r['L_minus_H_p3'], '+.3f')} | {_fmt(r['share_sq_p42'], '.2f')} / {_fmt(r['share_sq_p3'], '.2f')} "
                 f"({_fmt(r['share_sq_exact_p42'], '.2f')} / {_fmt(r['share_sq_exact_p3'], '.2f')}) |")
    L += ['', '## Step 2: the arms (calibrated readout, all recordings; d = changed run minus own run, mean of four wirings)', '',
          '| arm | host | class | d | the four wirings | SD / SE | halves A / B / C / D | 42 / 78 / 3 | amplitude / shape | raw |',
          '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in arms.iterrows():
        cls = r['cls'] + (f" ({r['side']})" if r['side'] else '')
        L.append(f"| {r['arm']} | {SHORT[r['host']]} | {cls} | {_fmt(r['d_pattern'], '+.4f')} | "
                 + ' / '.join(_fmt(r[f'd_seed{s}'], '+.4f') for s in HOST_SEEDS) + ' | '
                 f"{_fmt(r['sd_wirings'])} / {_fmt(r['se'])} | " + ' / '.join(_fmt(r[f'd_half_{t}'], '+.4f') for t in HALVES)
                 + f" | {_fmt(r['d_p42'], '+.4f')} / {_fmt(r['d_p78'], '+.4f')} / {_fmt(r['d_p3'], '+.4f')} | "
                 f"{_fmt(r['d_amplitude'], '+.4f')} / {_fmt(r['d_shape'], '+.4f')} | {_fmt(r['d_raw_pattern'], '+.4f')} |")
    L += ['', '### The classes with the local gap in G\'s place (a reading beside the rule)', '',
          '| arm | class (G) | class (the local gap) | a part on the left / hemibrain (local) |', '|---|---|---|---|']
    for _, r in yard.iterrows():
        L.append(f"| {r['arm']} | {r['cls_G']}{' (' + r['side_G'] + ')' if r['side_G'] else ''} | {r['cls_local']}"
                 f"{' (' + r['side_local'] + ')' if r['side_local'] else ''} | "
                 f"{'yes' if r['part_local_' + LEFT] else 'no'} / {'yes' if r['part_local_' + HEMI] else 'no'} |")
    L += ['', '### Derived readings', '', '| quantity | value |', '|---|---:|']
    for k, v in der.items():
        L.append(f'| {k} | {_fmt(v, "+.6f")} |')
    L += ['', '## Realism', '', '| runs | n | pattern part | rate on | rate off | off / on | rate ρ on | rate ρ off | '
          'odors per KC on / off | in bands |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _, r in rmain.iterrows():
        if r['runs'] == 'KC data':
            L.append(f"| KC data | | floor {_fmt(r['pattern_min'], '.3f')} | {_fmt(r['rate_on_min'], '.3f')} (band "
                     f"{r['band_rate_on']}) | {_fmt(r['rate_off_min'], '.3f')} | {_fmt(r['ratio_min'], '.2f')} ({r['ratio_ci']}; "
                     f"band {r['band_ratio']}) | ceiling {_fmt(r['rate_rho_on_min'], '.2f')}; band {r['band_rate_rho_on']} | "
                     f"ceiling {_fmt(r['rate_rho_off_min'], '.2f')} | bands {_fmt(r['opk_on_tv_min'], '.3f')} / "
                     f"{_fmt(r['opk_off_tv_min'], '.3f')} | |")
            continue
        L.append(f"| {r['runs']} | {int(r['n'])} | {_rng(r['pattern_min'], r['pattern_max'])} | "
                 f"{_rng(r['rate_on_min'], r['rate_on_max'])} | {_rng(r['rate_off_min'], r['rate_off_max'])} | "
                 f"{_rng(r['ratio_min'], r['ratio_max'], '.2f')} | {_rng(r['rate_rho_on_min'], r['rate_rho_on_max'], '.2f')} | "
                 f"{_rng(r['rate_rho_off_min'], r['rate_rho_off_max'], '.2f')} | "
                 f"{_rng(r['opk_on_tv_min'], r['opk_on_tv_max'], '.3f')} / {_rng(r['opk_off_tv_min'], r['opk_off_tv_max'], '.3f')} | "
                 f"{int(r['n_in_bands'])} of {int(r['n'])} |")
    L += ['', '| arm | host | γ off / on (own) | threshold × own | APL→KC × own | largest count × own, max | '
          'rate ρ on | mean APL effect |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in rarms.iterrows():
        if r['arm'] not in ARMS:
            L.append(f"| {r['arm']} | {SHORT[r['host']]} | {_rng(r['ratio_g_min'], r['ratio_g_max'], '.2f')} | "
                     f"hold factor {_rng(r['thr_factor_min'], r['thr_factor_max'], '.2f')} | hold factor "
                     f"{_rng(r['apl_factor_min'], r['apl_factor_max'], '.2f')} | – | – | "
                     f"{_rng(r['mean_effect_min'], r['mean_effect_max'], '+.3f')} |")
            continue
        L.append(f"| {r['arm']} | {SHORT[r['host']]} | {_rng(r['ratio_g_min'], r['ratio_g_max'], '.2f')} "
                 f"({_rng(r['ratio_g_own_min'], r['ratio_g_own_max'], '.2f')}) | {_rng(r['thr_vs_own_min'], r['thr_vs_own_max'], '.3f')} | "
                 f"{_rng(r['aplkc_vs_own_min'], r['aplkc_vs_own_max'], '.3f')} | {_fmt(r['max_count_vs_own_max'], '.2f')} | "
                 f"{_rng(r['rate_rho_on_min'], r['rate_rho_on_max'], '.2f')} | {_rng(r['mean_effect_min'], r['mean_effect_max'], '+.3f')} |")
    o = me[me['outside']]
    L += ['', f'Mean APL effect under the readout outside the flies\' value ±{MEAN_EFFECT_TOL}: {len(o)} of {len(me)} runs '
          f'({int((o["calyx"] != RIGHT).sum())} on the two hosts, {int((o["calyx"] == RIGHT).sum())} on the male CNS right).', '',
          '| step | run | mean APL effect | minus the flies\' |', '|---:|---|---:|---:|']
    for _, r in o.iterrows():
        L.append(f"| {r['step']} | {r['tag']} | {_fmt(r['mean_effect'], '+.4f')} | {_fmt(r['minus_flies'], '+.4f')} |")
    L += ['', "## B's factors per glomerulus", '',
          '| glomerulus | ORN input | boutons hb / left / right | APL→bouton per bouton × mean, hb / left / right | '
          'bouton→APL per bouton × mean, hb / left / right | factor on the left host, APL→bouton / bouton→APL | '
          'factor on the hemibrain host |', '|---|---|---:|---:|---:|---:|---:|']
    B_ = bf.set_index('glomerulus')
    for gname in B_GLOMERULI:
        r = B_.loc[gname]
        L.append(f"| {gname} | {'yes' if r['orn_input'] else 'no'} | "
                 + ' / '.join(str(int(r[f'boutons_{c}'])) for c in CALYCES) + ' | '
                 + ' / '.join(_fmt(r[f'aplpn_rel_{c}'], '.2f') for c in CALYCES) + ' | '
                 + ' / '.join(_fmt(r[f'pnapl_rel_{c}'], '.2f') for c in CALYCES) + ' | '
                 f"{_fmt(r[f'aplpn_factor_{LEFT}'], '.2f')} / {_fmt(r[f'pnapl_factor_{LEFT}'], '.2f')} | "
                 f"{_fmt(r[f'aplpn_factor_{HEMI}'], '.2f')} / {_fmt(r[f'pnapl_factor_{HEMI}'], '.2f')} |")
    L += ['', '| host | direction | share on the glomeruli with input, before / after | factor range among the 36 | '
          'largest factor |', '|---|---|---:|---:|---:|']
    for _, r in bsh.iterrows():
        L.append(f"| {SHORT[r['host']]} | {r['direction']} | {_fmt(r['share_with_input_before'], '.3f')} / "
                 f"{_fmt(r['share_with_input_after'], '.3f')} | {_fmt(r['factor_min_with_input'], '.2f')} "
                 f"({r['factor_min_with_input_glomerulus']}) to {_fmt(r['factor_max_with_input'], '.2f')} "
                 f"({r['factor_max_with_input_glomerulus']}) | {_fmt(r['factor_max_all'], '.2f')} ({r['factor_max_all_glomerulus']}) |")
    L += ['', "## The planners' numbers", '', f"{int(pc['reproduced'].sum())} of {len(pc)} reproduce at the precision "
          'written ("cut off" where the written value is the script\'s cut off at that precision instead of rounded).', '',
          '| reading | what | written | script | reproduced |', '|---|---|---:|---:|---|']
    for _, r in pc.iterrows():
        yn = 'yes' if r['reproduced'] else ('no (cut off)' if r['truncated_matches'] else 'no')
        L.append(f"| {r['reading']} | {r['what']} | {r['written']} | {_fmt(r['script'], '+.6f')} | {yn} |")
    return '\n'.join(L) + '\n'


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', type=Path, default=None, help='output folder (default: the results folder)')
    a = ap.parse_args(argv)
    res = H.results_dir(EXPERIMENT)
    out = a.out or res
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    I = load(res)
    ck = Checks()
    g0, sp0 = step0_groups(I, ck)
    red, gaps = step1(I, ck)
    arms, derived, der, yard = step2(I, ck)
    rmain, rarms, me, _ = realism(I, ck)
    sz = sizes(I)
    bf, bsh = b_factors(I, ck)
    pc = planner_check(g0, sp0, red, gaps, arms, der, rmain, rarms, bf, bsh, I)
    tables = dict(step0_groups=g0, step0_gap_split=sp0, step1_redraws=red, step1_gaps=gaps, step2_arms=arms,
                  step2_derived=derived, step2_local_yardstick=yard, realism=rmain, realism_arms=rarms, mean_effect_outside=me, sizes=sz,
                  b_factors=bf, b_shares=bsh, planner_check=pc)
    for k, df in tables.items():
        df.to_csv(out / f'writeup_{k}.csv', index=False)
    (out / 'writeup_tables.md').write_text(markdown(g0, sp0, red, gaps, arms, der, yard, rmain, rarms, me, bf, bsh, pc))
    inputs = ['readings1.csv', 'runs1.csv', 'readings2.csv', 'runs2.csv', 'd2.csv', 'verdict1.json', 'verdict2.json',
              'step0.json', 'step0_patterns.csv', 'step0_glomeruli.csv', 'step0_glomeruli_by_glomerulus.csv',
              'builder_sizes.csv']
    prov = dict(script=SCRIPT, generated=time.strftime('%Y-%m-%d %H:%M'), seconds=round(time.time() - t0, 1),
                checks_passed=ck.n, planner_numbers=int(len(pc)), planner_reproduced=int(pc['reproduced'].sum()),
                not_reproduced=pc[~pc['reproduced']][['reading', 'what', 'written', 'script',
                                                      'truncated_matches']].to_dict('records'),
                commits=H.code_commits(),
                inputs={**{f'results/{EXPERIMENT}/{n}': _sha256(res / n) for n in inputs},
                        'results/data_targets/data_targets.json':
                            _sha256(H.shared_dir() / 'results' / 'data_targets' / 'data_targets.json')},
                note=("No simulation. B's factors are E133_builders.build('B') on light_tables (the model's tables with "
                      "APL's counted synapses), checked against builder_sizes.csv and step 0's glomerulus tables."))
    (out / 'provenance_writeup_tables.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    print(f'{ck.n} checks passed; planner numbers reproduced {int(pc["reproduced"].sum())} of {len(pc)}; wrote {out}',
          flush=True)
    for _, r in pc[~pc['reproduced']].iterrows():
        cut = ' (cut off, not rounded)' if r['truncated_matches'] else ''
        print(f"  not reproduced: {r['reading']}: {r['what']}: written {r['written']}, script {r['script']:.6f}{cut}",
              flush=True)


if __name__ == '__main__':
    main()
