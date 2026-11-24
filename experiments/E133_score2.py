#!/usr/bin/env python3
"""E133 step 2's score: the transplants' d's, the arms' classes and step 3's arms (no simulation).

The entry is item E133 of experiment_priorities.md ("The classes", step 2; the Design's step 3; the Status line's step-1
gate reading of 2026-11-24 07:43, which asks for each arm's d on the 42 and the 78 pairs beside its class). This module
is separate from experiments/E133_calyx_transplants.py so that it can be written and committed while step-2 lanes run
from the same checkout: no lane loads it (its test checks the import graph), so E108's `_still_clean` guard never sees
it change. It imports the runner's readings and never re-implements them: a run's pair distances are
E133_calyx_transplants._dist (E83's `pair_distances`, raw counts or the calibrated readout with K re-calibrated on the
run's APL-on counts) and its readings E133_calyx_transplants.unit_readings (E110's `scores_of`: the loss, the pattern
part, its amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair APL effect, against all
recordings and the halves A to D; and `group_patterns`, the pattern error on the 42 / 78 / 39 / 81 / 3 pairs).

What it reads:
- step 2's 48 points (E133_calyx_transplants.planned_step2: the six arms C, S, A, B, K, ALL on each host's four
  wirings), their rows in every host's and lane's runs_2_*.csv and their saved runs (outputs/E133_calyx_transplants/
  step2/ or the model-output archive's step2/). Rows of other arms (step 3's options, step 2′'s PA and CI) are listed
  and not scored.
- step 1's readings and runs (readings1.csv, runs1.csv): the hosts' own runs at seeds 0 to 3 (seed 0 E130's
  fixed-scale run S1_<calyx>, seeds 1 to 3 step 1's re-draws) and the eight wirings of each calyx behind G. The own
  runs' readings are recomputed here once and compared with readings1.csv (`own_check`), so that both sides of every d
  come from the same code.

The quantities (the entry's, on the calibrated pattern part against all recordings; each half beside):
- G: the mean over the left's eight wirings (its own and the seven re-draws) minus the hemibrain's, from step 1's
  readings. The local gap, the same over the four wirings of step 2 (seeds 0 to 3), is reported beside it.
- d, per arm and host: the mean over the host's four wirings of the changed run minus that wiring's own run.

The classes (the entry's rules; `host_reading`, `arm_class`):
- "closes" (host the left): d at most -G/2 and at most -0.004, the changed run lower on 3 of the 4 wirings, and d
  negative against every half. "Opens" (host the hemibrain) is the mirror.
- "carries": it opens and closes. "a part": at least G/4 and at least 0.003, in those directions on both hosts, by the
  same rules (3 of the 4 wirings, every half). "one-sided": it closes or opens, and is under G/4 on the other host.
  "none": under 0.003 on both hosts. Anything else "unclear".
- no class (cls None, the would-be class recorded) for an arm any of whose runs, or of its hosts' own runs at seeds 0
  to 3, is outside a band, at a corner or over the guard (the runs table's `no_class`, set by the runner).
- step 3: the arms that carry, are a part or are one-sided, at most the three largest (`step3_selection`).

Readings the entry leaves open, fixed here before step 2's runs are scored:
- "Under G/4 on the other host" and "under 0.003 on both hosts" are read on |d|, the size of the move whichever its
  direction (the two phrases are worded alike, and "none" can only mean a small move). The class under the other
  reading of "under G/4" (the move toward the other calyx under G/4, any move away counting as under) is recorded
  beside as 'cls_directional' where it differs.
- "Lower on 3 of the 4 wirings": at least 3 of the 4 d's strictly below 0 (above 0 for the hemibrain host). "Negative
  against every half": strictly, on each of the halves A to D. The thresholds are inclusive ("at most", "at least").
- d's mean is an exactly rounded sum (math.fsum) over the four wirings, divided by 4.
- "The three largest" are ranked by the arm's larger |d| of its two hosts (calibrated pattern part, all recordings);
  exact ties by its other |d|, then carries before a part before one-sided, then the order C, S, A, B, K, ALL. ALL is
  an arm of step 2 and is eligible; the list without it is recorded beside ('selected_without_ALL'), since step 3's
  splits name single arms only.
- An arm with fewer than its eight runs is 'incomplete' (no class). Without all 48 runs the script writes nothing
  unless `--partial`, which writes a provisional table (files *_provisional.*, every heading marked, no class and no
  step-3 list: complete arms show only a 'provisional reading').

Outputs (results/E133_calyx_transplants/; with --partial the same names with '_provisional', removed by the final
scoring): readings2.csv (per step-2 run and readout, `unit_readings`), runs2.csv (the rows of the 48 points, with the
host wiring's own run's realism columns and parameters beside, '<col>_own' and '<param>_vs_own'), d2.csv (per arm,
host, readout and statistic: the four wirings' d's, their mean and SD, the SD of the host's seven re-draws and of its
eight wirings), verdict2.json (G, the classes, step 3's list, the runs without a class, the realism flags, the own
check, provenance) and step2_table.md.

Usage (from the repository root; no simulation, light):
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 19 ./venv/bin/python -u experiments/E133_score2.py
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 19 ./venv/bin/python -u experiments/E133_score2.py --partial
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E133_calyx_transplants as E
from kcq import harness as H

SCRIPT = 'experiments/E133_score2.py'
ARMS = ('C', 'S', 'A', 'B', 'K', 'ALL')
LEFT, HEMI = E.LEFT, E.HEMI
HOSTS = (LEFT, HEMI)                       # the class reads the left first (closes), then the hemibrain (opens)
SIGN = {LEFT: -1.0, HEMI: 1.0}             # the direction of a host's d toward the other calyx
VERB = {LEFT: 'closes', HEMI: 'opens'}
SEEDS = E.HOST_SEEDS                       # (0, 1, 2, 3)
TARGETS = E.TARGETS_READ                   # ('all', 'A', 'B', 'C', 'D')
HALVES = TARGETS[1:]
# the entry's thresholds: (fraction of G, floor)
FULL = (0.5, 0.004)                        # closes / opens
PART = (0.25, 0.003)                       # a part
ONE_SIDED_OTHER = 0.25                     # one-sided: under G/4 on the other host
NONE_MAX = 0.003                           # none: under 0.003 on both hosts
MIN_TOWARD = 3                             # lower (higher) on 3 of the 4 wirings
STEP3_MAX = 3
ELIGIBLE = ('carries', 'a part', 'one-sided')
CLASS_ORDER = {c: i for i, c in enumerate(ELIGIBLE)}
OWN_TOL = 1e-12
CAL_STATS = (tuple(f'{t}_pattern' for t in TARGETS) + tuple(f'{t}_pattern_{g}' for t in TARGETS for g in E.GROUPS)
             + ('all_amplitude', 'all_shape', 'all_sd_ratio', 'all_pearson_effect', 'all_rho_effect', 'mean_effect',
                'all_loss'))
RAW_STATS = ('all_pattern', 'all_amplitude', 'all_shape', 'all_sd_ratio', 'all_pearson_effect', 'mean_effect', 'all_loss')
OWN_CHECK_COLS = tuple(f'{t}_pattern' for t in TARGETS) + tuple(f'all_pattern_{g}' for g in E.GROUPS) + (
    'all_amplitude', 'all_shape', 'mean_effect')
OWN_BESIDE = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'ratio_g', 'rate_on_g', 'rate_on_apbp', 'rate_on_ab',
              'rate_off_g', 'rate_off_apbp', 'rate_off_ab', 'max_count_on', 'max_count_off', 'n_kc_ge20_on', 'n_kc_ge20_off',
              'thr_factor', 'apl_factor') + E.PARAM_KEYS
OUTPUTS = ('readings2.csv', 'runs2.csv', 'd2.csv', 'verdict2.json', 'step2_table.md')
SIZE_COLS = {'C': ('size_C_off_share', 'size_C_copied_share', 'size_C_claws_per_kc_g_host', 'size_C_claws_per_kc_g_after',
                   'size_C_claws_per_kc_ab_host', 'size_C_claws_per_kc_ab_after'),
             'S': ('size_S_moved_share', 'size_S_n_glomeruli_changed', 'size_S_r_min', 'size_S_r_max'),
             'A': ('size_A_aplkc_factor_min', 'size_A_aplkc_factor_max', 'size_A_kcapl_factor_min', 'size_A_kcapl_factor_max'),
             'B': ('size_B_aplpn_factor_min', 'size_B_aplpn_factor_max', 'size_B_pnapl_factor_min', 'size_B_pnapl_factor_max'),
             'K': ('size_K_within_type_share_host', 'size_K_within_type_share_after', 'size_K_rescale')}
SIZE_COLS['ALL'] = tuple(c for a in ('S', 'C', 'A', 'B', 'K') for c in SIZE_COLS[a])


def out_name(name: str, provisional: bool) -> str:
    """'d2.csv' -> 'd2_provisional.csv' when `provisional`."""
    if not provisional:
        return name
    stem, ext = name.rsplit('.', 1)
    return f'{stem}_provisional.{ext}'


def _sha256_file(p: Path) -> Optional[str]:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def _mean_exact(x: Sequence[float]) -> float:
    x = [float(v) for v in x]
    return math.fsum(x) / len(x) if x else float('nan')


def _sd(x: Sequence[float]) -> float:
    x = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    return float(x.std(ddof=1)) if x.size >= 2 else float('nan')


# ---------------------------------------------------------------------------
# G and d

def gap(rd1: pd.DataFrame, targets: Sequence[str] = TARGETS) -> Dict[str, dict]:
    """G per target t (the left's eight wirings, seeds 0 to 7, minus the hemibrain's, mean over each; calibrated
    '<t>_pattern' of step 1's readings) and the local gap (the same over seeds 0 to 3). NaN without all eight (four)."""
    cal = rd1[rd1['readout'] == 'cal']
    out = {}
    eight = (0,) + tuple(E.SEEDS1)
    for t in targets:
        col = f'{t}_pattern'

        def vals(c, seeds):
            s = cal[cal['calyx'] == c].drop_duplicates('seed', keep='last').set_index('seed')[col]
            return [float(s[k]) for k in seeds if k in s.index and np.isfinite(s[k])]
        L8, H8, L4, H4 = vals(LEFT, eight), vals(HEMI, eight), vals(LEFT, SEEDS), vals(HEMI, SEEDS)
        ok8, ok4 = len(L8) == len(eight) and len(H8) == len(eight), len(L4) == len(SEEDS) and len(H4) == len(SEEDS)
        out[t] = dict(G=float(np.mean(L8) - np.mean(H8)) if ok8 else float('nan'),
                      L8=float(np.mean(L8)) if L8 else float('nan'), H8=float(np.mean(H8)) if H8 else float('nan'),
                      local=float(np.mean(L4) - np.mean(H4)) if ok4 else float('nan'),
                      L4=float(np.mean(L4)) if L4 else float('nan'), H4=float(np.mean(H4)) if H4 else float('nan'),
                      n_left=len(L8), n_hemi=len(H8))
    return out


def d_table(rd2: pd.DataFrame, rd1: pd.DataFrame, arms: Sequence[str] = ARMS, hosts: Sequence[str] = HOSTS,
            seeds: Sequence[int] = SEEDS) -> pd.DataFrame:
    """Per arm, host, readout and statistic (CAL_STATS, RAW_STATS): the changed run minus its wiring's own run at each
    seed ('d_seed<n>', NaN where either is missing), 'n' (wirings with both), 'd' (their mean), 'sd_wirings' (SD of the
    d's), and the host's own scatter of the statistic in step 1: 'sd_redraws' (its seven re-draws, seeds 1 to 7) and
    'sd_eight' (with its own wiring). `rd2` has columns arm, calyx, seed, readout; `rd1` calyx, seed, readout."""
    a = rd2.drop_duplicates(['arm', 'calyx', 'seed', 'readout'], keep='last').set_index(['arm', 'calyx', 'seed', 'readout'])
    b = rd1.drop_duplicates(['calyx', 'seed', 'readout'], keep='last').set_index(['calyx', 'seed', 'readout'])
    rows = []
    for arm in arms:
        for h in hosts:
            for ro, cols in (('cal', CAL_STATS), ('raw', RAW_STATS)):
                for col in cols:
                    dw = {}
                    for s in seeds:
                        ka, kb = (arm, h, int(s), ro), (h, int(s), ro)
                        if ka in a.index and kb in b.index and col in a.columns and col in b.columns:
                            x, y = float(a.loc[ka, col]), float(b.loc[kb, col])
                            if np.isfinite(x) and np.isfinite(y):
                                dw[int(s)] = x - y
                    red = [float(b.loc[(h, s, ro), col]) for s in E.SEEDS1 if (h, s, ro) in b.index and col in b.columns]
                    own = [float(b.loc[(h, 0, ro), col])] if (h, 0, ro) in b.index and col in b.columns else []
                    rows.append(dict(arm=arm, host=h, readout=ro, stat=col, n=len(dw),
                                     **{f'd_seed{int(s)}': dw.get(int(s), float('nan')) for s in seeds},
                                     d=_mean_exact(list(dw.values())), sd_wirings=_sd(list(dw.values())),
                                     sd_redraws=_sd(red), sd_eight=_sd(own + red)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the classes

def host_reading(d_w: Sequence[float], d_halves: Mapping[str, float], sign: float, G: float,
                 n_wirings: int = len(SEEDS)) -> dict:
    """One arm on one host: d (the exact mean of the wirings' d's `d_w`), 'd_toward' (sign x d: positive toward the
    other calyx; sign -1 for the left, whose closing lowers the pattern part), 'n_toward' (wirings moved that way), the
    halves' d's toward ('every_half': all strictly positive), and the rules: 'full' (closes / opens: toward by at least
    G/2 and 0.004, on 3 of the 4 wirings, every half), 'part' (the same at G/4 and 0.003), 'under_quarter' (|d| < G/4),
    'under_quarter_directional' (d_toward < G/4), 'under_none' (|d| < 0.003). 'complete' needs all `n_wirings` d's and
    every half's; the rules are False otherwise."""
    d_w = [float(x) for x in d_w if np.isfinite(x)]
    hv = {k: float(v) for k, v in d_halves.items()}
    complete = bool(len(d_w) == n_wirings and set(hv) >= set(HALVES) and all(np.isfinite(hv[k]) for k in HALVES)
                    and np.isfinite(G))
    d = _mean_exact(d_w)
    dt = sign * d
    n_toward = int(sum(sign * x > 0 for x in d_w))
    every_half = bool(complete and all(sign * hv[k] > 0 for k in HALVES))

    def meets(frac: float, floor: float) -> bool:
        return bool(complete and dt >= frac * G and dt >= floor and n_toward >= MIN_TOWARD and every_half)
    return dict(d=d, d_toward=dt, n_wirings=len(d_w), n_toward=n_toward, d_halves=hv, every_half=every_half,
                complete=complete, full=meets(*FULL), part=meets(*PART),
                under_quarter=bool(complete and abs(d) < ONE_SIDED_OTHER * G),
                under_quarter_directional=bool(complete and dt < ONE_SIDED_OTHER * G),
                under_none=bool(complete and abs(d) < NONE_MAX))


def arm_class(left: Mapping, hemi: Mapping, directional: bool = False) -> Tuple[str, Optional[str]]:
    """(class, side) of an arm from its two host readings (`host_reading`; the left's and the hemibrain's): 'carries'
    (closes and opens), 'a part' (both at the part rule), 'one-sided' (closes or opens, under G/4 on the other host:
    |d|, or with `directional` the move toward), 'none' (|d| under 0.003 on both), 'unclear' otherwise; 'incomplete'
    without both hosts' four wirings and halves. side: 'closes (left host)' / 'opens (hemibrain host)' for one-sided."""
    if not (left['complete'] and hemi['complete']):
        return 'incomplete', None
    uq = 'under_quarter_directional' if directional else 'under_quarter'
    if left['full'] and hemi['full']:
        return 'carries', None
    if left['part'] and hemi['part']:
        return 'a part', None
    if left['full'] and hemi[uq]:
        return 'one-sided', 'closes (left host)'
    if hemi['full'] and left[uq]:
        return 'one-sided', 'opens (hemibrain host)'
    if left['under_none'] and hemi['under_none']:
        return 'none', None
    return 'unclear', None


def _halves_d(dt: pd.DataFrame, arm: str, h: str) -> Dict[str, float]:
    x = dt[(dt['arm'] == arm) & (dt['host'] == h) & (dt['readout'] == 'cal')].set_index('stat')
    return {t: float(x.loc[f'{t}_pattern', 'd']) if f'{t}_pattern' in x.index and x.loc[f'{t}_pattern', 'n'] == len(SEEDS)
            else float('nan') for t in HALVES}


def arm_readings(dt: pd.DataFrame, G: float, arms: Sequence[str] = ARMS) -> Dict[str, dict]:
    """Per arm: each host's `host_reading` (calibrated pattern part, all recordings, the halves beside) and the class
    (`arm_class`, magnitude reading; 'cls_directional' the other reading of "under G/4")."""
    out = {}
    for arm in arms:
        hr = {}
        for h in HOSTS:
            x = dt[(dt['arm'] == arm) & (dt['host'] == h) & (dt['readout'] == 'cal') & (dt['stat'] == 'all_pattern')]
            d_w = [float(x[f'd_seed{s}'].iloc[0]) for s in SEEDS] if len(x) else []
            hr[h] = host_reading(d_w, _halves_d(dt, arm, h), SIGN[h], G)
        cls, side = arm_class(hr[LEFT], hr[HEMI])
        cd, sd = arm_class(hr[LEFT], hr[HEMI], directional=True)
        out[arm] = dict(rule_class=cls, side=side, hosts=hr,
                        **({'cls_directional': cd, 'side_directional': sd} if (cd, sd) != (cls, side) else {}))
    return out


def classes(readings: Mapping[str, dict], no_class_tags: Mapping[str, Sequence[str]], G: float,
            provisional: bool = False) -> Dict[str, dict]:
    """The arms' classes: 'cls' the rule class, or None with 'would_be' when any of the arm's runs or its hosts' own
    runs has no class (`no_class_tags`: {arm: tags}, the own runs' under every arm), when G is not positive, or when
    `provisional` (then 'provisional_reading' holds the rule class)."""
    out = {}
    for arm, r in readings.items():
        rec = dict(r)
        cls = r['rule_class']
        nc = list(no_class_tags.get(arm, ()))
        if not (np.isfinite(G) and G > 0):
            rec.update(cls=None, would_be=cls, reason=f'G is {G} (not positive)')
        elif provisional:
            rec.update(cls=None, provisional_reading=cls, reason='provisional: not all of step 2\'s runs are in')
        elif cls == 'incomplete':
            rec.update(cls=None, would_be=cls, reason='not all eight runs are in')
        elif nc:
            rec.update(cls=None, would_be=cls, reason='a run without a class (bands, corner or guard)', no_class_runs=nc)
        else:
            rec['cls'] = cls
        out[arm] = rec
    return out


def step3_selection(arms: Mapping[str, Mapping], k: int = STEP3_MAX, include_all: bool = True,
                    key: str = 'cls') -> List[str]:
    """Step 3's arms: those whose `key` ('cls', or 'would_be' for the list a no-class run took away) carries, is a
    part or is one-sided, at most `k`, the largest first: by the larger |d| of the two hosts (calibrated pattern part,
    all recordings), exact ties by the other |d|, then by class (carries, a part, one-sided), then by ARMS order."""
    cand = [a for a, r in arms.items() if r.get(key) in ELIGIBLE and (include_all or a != 'ALL')]

    def order(a):
        r = arms[a]
        m = sorted((abs(float(r['hosts'][h]['d'])) for h in HOSTS), reverse=True)
        return (-m[0], -m[1], CLASS_ORDER[r[key]], ARMS.index(a) if a in ARMS else len(ARMS))
    return sorted(cand, key=order)[:k]


# ---------------------------------------------------------------------------
# collecting the runs

def _runs_table2(tries: int = 3, wait: float = 5.0) -> pd.DataFrame:
    """E133_calyx_transplants.runs_table(2), read again after `wait` s if a lane's append was caught half written."""
    import time
    for i in range(tries):
        try:
            return E.runs_table(2)
        except (pd.errors.ParserError, ValueError, OSError):
            if i == tries - 1:
                raise
            time.sleep(wait)
    raise AssertionError('unreachable')


def collect2(points: Optional[Sequence] = None) -> Tuple[Dict[str, Dict[str, tuple]], pd.DataFrame, List[dict], List[str]]:
    """Step 2's planned points that are in: ({tag: outputs}, their runs-table rows, the missing points ({tag, why}),
    the tags of other step-2 rows (step 3's options, step 2′), not scored here)."""
    pts = list(E.planned_step2() if points is None else points)
    tags = [E.tag_of(p) for p in pts]
    runs = _runs_table2()
    ri = runs.set_index('tag') if len(runs) else pd.DataFrame()
    outs, rows, missing = {}, [], []
    for pt, tag in zip(pts, tags):
        if tag not in ri.index:
            missing.append(dict(tag=tag, why='no row in runs_2_*.csv'))
            continue
        r = ri.loc[tag].to_dict()
        if (r.get('calyx'), int(r.get('seed')), str(r.get('arm'))) != (pt.calyx, int(pt.seed), pt.arm):
            raise RuntimeError(f'{tag}: the row says {r.get("calyx")}, seed {r.get("seed")}, arm {r.get("arm")}')
        o = E.load_run(tag, 'step2')
        if o is None:
            missing.append(dict(tag=tag, why='a row but no outputs here or in the archive'))
            continue
        outs[tag] = o
        rows.append(dict(r, tag=tag))
    extra = sorted(set(ri.index) - set(tags)) if len(ri) else []
    return outs, pd.DataFrame(rows), missing, extra


def readings2(outs: Mapping[str, Dict[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per step-2 run and readout: the arm, host (calyx) and seed, and E133_calyx_transplants.unit_readings."""
    idx = runs.set_index('tag')
    rows = []
    for ro in E.READOUTS:
        for tag, o in outs.items():
            rows.append(dict(tag=tag, arm=idx.loc[tag, 'arm'], calyx=idx.loc[tag, 'calyx'], seed=int(idx.loc[tag, 'seed']),
                             readout=ro, **E.unit_readings(E._dist(o, ro))))
    E._save_dist_cache()
    return pd.DataFrame(rows)


def own_outputs(c: str, seed: int) -> Dict[str, tuple]:
    """A host wiring's own run: E130's S1_<calyx> at seed 0, step 1's re-draw otherwise."""
    if seed == 0:
        return E.e130_s1(c)
    o = E.load_run(E.tag_of(E.Point(1, c, int(seed), None)), 'step1')
    if o is None:
        raise RuntimeError(f'S1_{c}_seed{seed}: no outputs here or in the archive')
    return o


def own_check(rd1: pd.DataFrame, tol: float = OWN_TOL) -> dict:
    """The hosts' own runs at seeds 0 to 3 read again (both readouts) against readings1.csv: the largest absolute
    difference over OWN_CHECK_COLS; 'passes' at most `tol`."""
    b = rd1.drop_duplicates(['calyx', 'seed', 'readout'], keep='last').set_index(['calyx', 'seed', 'readout'])
    worst, where, n = 0.0, None, 0
    for c in HOSTS:
        for s in SEEDS:
            o = own_outputs(c, s)
            for ro in E.READOUTS:
                r = E.unit_readings(E._dist(o, ro))
                n += 1
                for col in OWN_CHECK_COLS:
                    x = abs(float(r[col]) - float(b.loc[(c, int(s), ro), col]))
                    x = x if np.isfinite(x) else float('inf')
                    if x > worst:
                        worst, where = x, f'{c} seed {s} {ro} {col}'
    E._save_dist_cache()
    return dict(max_abs_diff=worst, where=where, n_readings=n, tol=tol, passes=bool(worst <= tol))


def runs2_table(runs: pd.DataFrame, runs1: pd.DataFrame) -> pd.DataFrame:
    """The step-2 rows with their host wiring's own run beside: '<col>_own' for OWN_BESIDE, '<param>_vs_own' (the
    changed run's held parameter over the own run's) for E133's PARAM_KEYS, the hold's factors likewise."""
    own = runs1.drop_duplicates(['calyx', 'seed'], keep='last').set_index(['calyx', 'seed'])
    rows = []
    for r in runs.to_dict('records'):
        k = (r['calyx'], int(r['seed']))
        o = own.loc[k] if k in own.index else None
        r = dict(r, own_tag=E.tag_of(E.Point(1, r['calyx'], int(r['seed']), None)), own_no_class=(
            bool(o['no_class']) if o is not None else None))
        for col in OWN_BESIDE:
            r[f'{col}_own'] = o[col] if o is not None and col in o.index else np.nan
        for col in E.PARAM_KEYS + ('thr_factor', 'apl_factor'):
            ov = float(o[col]) if o is not None and col in o.index else np.nan
            r[f'{col}_vs_own'] = float(r[col]) / ov if np.isfinite(ov) and ov != 0 else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def no_class_tags(runs2: pd.DataFrame, runs1: pd.DataFrame, arms: Sequence[str] = ARMS) -> Dict[str, List[str]]:
    """{arm: the tags without a class among its runs and its hosts' own runs at seeds 0 to 3}."""
    own = runs1.drop_duplicates(['calyx', 'seed'], keep='last').set_index(['calyx', 'seed'])
    own_nc = [E.tag_of(E.Point(1, c, int(s), None)) for c in HOSTS for s in SEEDS
              if (c, int(s)) in own.index and bool(own.loc[(c, int(s)), 'no_class'])]
    out = {}
    for arm in arms:
        x = runs2[runs2['arm'] == arm] if len(runs2) else runs2
        out[arm] = [t for t, nc in zip(x.get('tag', []), x.get('no_class', [])) if bool(nc)] + own_nc
    return out


def realism_flags(runs2: pd.DataFrame) -> List[dict]:
    """The step-2 runs without a class or with a realism flag (E83's `realism`): their bands, corner and guard."""
    out = []
    for r in runs2.to_dict('records'):
        flag = r.get('realism_flag')
        flag = '' if flag is None or (isinstance(flag, float) and not np.isfinite(flag)) else str(flag)
        if bool(r.get('no_class')) or flag:
            out.append(dict(tag=r['tag'], flag=flag, no_class=bool(r.get('no_class')),
                            bands={b: bool(r.get(b)) for b in E.BAND_COLS}, corner=bool(r.get('corner')),
                            guard=bool(r.get('guard_max_count')), max_count_on=r.get('max_count_on'),
                            max_count_off=r.get('max_count_off'), ref_max_count_on=r.get('ref_max_count_on'),
                            ref_max_count_off=r.get('ref_max_count_off')))
    return out


# ---------------------------------------------------------------------------
# the verdict and the table

def verdict2(dt: pd.DataFrame, G: Mapping[str, dict], nct: Mapping[str, Sequence[str]], provisional: bool) -> dict:
    """G, the arms' classes and step 3's list (none while `provisional`)."""
    g = float(G['all']['G'])
    arms = classes(arm_readings(dt, g), nct, g, provisional)
    for arm, rec in arms.items():
        grp = {}
        for gname in E.GROUPS:
            x = dt[(dt['arm'] == arm) & (dt['readout'] == 'cal') & (dt['stat'] == f'all_pattern_{gname}')].set_index('host')['d']
            grp[gname] = {h: float(x.get(h, np.nan)) for h in HOSTS}
        raw = dt[(dt['arm'] == arm) & (dt['readout'] == 'raw') & (dt['stat'] == 'all_pattern')].set_index('host')['d']
        rec['groups'] = grp
        rec['raw_pattern'] = {h: float(raw.get(h, np.nan)) for h in HOSTS}
    v = dict(provisional=bool(provisional),
             rules=dict(readout='cal', column='all_pattern', halves=list(HALVES), closes_opens=dict(G_fraction=FULL[0], floor=FULL[1]),
                        part=dict(G_fraction=PART[0], floor=PART[1]), one_sided_other=f'|d| < {ONE_SIDED_OTHER} G',
                        none=f'|d| < {NONE_MAX} on both hosts', wirings_toward=f'{MIN_TOWARD} of {len(SEEDS)}',
                        sign=dict(left='closes: d < 0', hemibrain='opens: d > 0')),
             G={t: dict(G[t]) for t in TARGETS}, thresholds=dict(G_half=g / 2, G_quarter=g / 4, floor_full=FULL[1],
                                                                  floor_part=PART[1], none=NONE_MAX),
             arms=arms)
    if provisional:
        v['step3'] = dict(selected=None, note="withheld until all 48 of step 2's runs are in")
    else:
        v['step3'] = dict(
            selected=step3_selection(arms), selected_without_ALL=step3_selection(arms, include_all=False),
            would_be_without_no_class=step3_selection({a: dict(r, would_be=r.get('cls') or r.get('would_be')) for a, r in arms.items()},
                                                      key='would_be'),
            rule=(f'at most {STEP3_MAX} of the arms that carry, are a part or are one-sided, the largest first by the larger |d| '
                  'of the two hosts; exact ties by the other |d|, then carries, a part, one-sided, then C, S, A, B, K, ALL'))
        v['carries'] = [a for a, r in arms.items() if r.get('cls') == 'carries']
    return v


def _f(x, sp: str = '.4f') -> str:
    return E._f(x, sp)


def _yn(x) -> str:
    return 'yes' if bool(x) else 'no'


def _stat(dt: pd.DataFrame, arm: str, h: str, stat: str, ro: str = 'cal') -> pd.Series:
    x = dt[(dt['arm'] == arm) & (dt['host'] == h) & (dt['readout'] == ro) & (dt['stat'] == stat)]
    return x.iloc[0] if len(x) else pd.Series(dtype=object)


def _size_cell(rows: pd.DataFrame, arm: str) -> str:
    cells = []
    for col in SIZE_COLS.get(arm, ()):
        if col not in rows.columns:
            continue
        v = pd.to_numeric(rows[col], errors='coerce').dropna()
        if not len(v):
            continue
        lo, hi = float(v.min()), float(v.max())
        name = col.replace('size_', '')
        cells.append(f'{name} {lo:.3g}' if math.isclose(lo, hi, rel_tol=1e-9, abs_tol=1e-12) else f'{name} {lo:.3g} to {hi:.3g}')
    return '; '.join(cells)


def step2_table(v: Mapping, dt: pd.DataFrame, rd2: pd.DataFrame, rd1: pd.DataFrame, runs2: pd.DataFrame,
                missing: Sequence[Mapping], n_in: int, n_planned: int) -> str:
    """Step 2's tables as markdown."""
    prov = bool(v['provisional'])
    P = ' (PROVISIONAL)' if prov else ''
    G = v['G']
    L = [f'# E133 step 2{": PROVISIONAL, " + str(n_in) + " of " + str(n_planned) + " runs in, no class" if prov else " (generated)"}', '']
    if prov:
        L += ['**Provisional.** Not all of step 2\'s runs are in. No arm has a class and no step-3 list is given. '
              '"Provisional reading" is the rule\'s reading on the runs in so far, not a class.', '',
              f'Missing ({len(missing)}): ' + ', '.join(m['tag'] for m in missing) + '.', '']
    L += [f"G (the left's eight wirings minus the hemibrain's; step 1's readings, calibrated pattern part): {_f(G['all']['G'], '+.4f')} "
          f"(L {_f(G['all']['L8'])}, H {_f(G['all']['H8'])}); halves A / B / C / D "
          + ' / '.join(_f(G[t]['G'], '+.4f') for t in HALVES) + '.', '',
          f"The local gap on the four wirings of step 2 (seeds 0 to 3): {_f(G['all']['local'], '+.4f')} "
          f"(L {_f(G['all']['L4'])}, H {_f(G['all']['H4'])}); halves " + ' / '.join(_f(G[t]['local'], '+.4f') for t in HALVES) + '.', '',
          f"Thresholds: closes / opens at G/2 = {_f(v['thresholds']['G_half'])} and {FULL[1]}; a part at G/4 = "
          f"{_f(v['thresholds']['G_quarter'])} and {PART[1]}; one-sided: |d| under G/4 on the other host; none: |d| under "
          f"{NONE_MAX} on both. Both with 3 of the 4 wirings and every half in the direction.", '',
          f'## Classes{P}', '',
          '| arm | class | d, left host | d, hemibrain host | left: lower on / every half | hemibrain: higher on / every half '
          '| closes / part, left | opens / part, hemibrain |',
          '|---|---|---:|---:|---:|---:|---|---|']
    for arm, r in v['arms'].items():
        hl, hh = r['hosts'][LEFT], r['hosts'][HEMI]
        if prov:
            c = f"provisional reading: {r.get('provisional_reading')}"
        elif r.get('cls'):
            c = r['cls'] + (f" ({r['side']})" if r.get('side') else '')
        else:
            c = f"no class (would be {r.get('would_be')})"
        if r.get('cls_directional'):
            c += f"; {r['cls_directional']} under the directional reading"
        L.append(f"| {arm} | {c} | {_f(hl['d'], '+.4f')} | {_f(hh['d'], '+.4f')} | "
                 f"{hl['n_toward']} of {hl['n_wirings']} / {_yn(hl['every_half'])} | {hh['n_toward']} of {hh['n_wirings']} / "
                 f"{_yn(hh['every_half'])} | {_yn(hl['full'])} / {_yn(hl['part'])} | {_yn(hh['full'])} / {_yn(hh['part'])} |")
    s3 = v['step3']
    L += ['', f"Step 3: {s3['note'] if prov else (', '.join(s3['selected']) or 'no arm')}"
          + ('' if prov else f"; without ALL: {', '.join(s3['selected_without_ALL']) or 'no arm'}. Rule: {s3['rule']}."), '']
    if not prov:
        nc = [(a, r.get('no_class_runs')) for a, r in v['arms'].items() if r.get('no_class_runs')]
        L += ['Arms without a class for a run outside the bands, at a corner or over the guard: '
              + ('; '.join(f'{a} ({", ".join(t)})' for a, t in nc) if nc else 'none') + '.', '']
    rf = v.get('realism_flags') or []
    L += ['Step-2 runs without a class or with a realism flag: '
          + ('; '.join(f"{x['tag']} ({x['flag'] or 'no class'})" for x in rf) if rf else 'none') + '.', '']
    L += [f'## d by pair group, calibrated, all recordings{P}', '',
          'Each cell: the left host / the hemibrain host. The last column is the pattern part on raw counts.', '',
          '| arm | 120 pairs | 42 with 2-but, paa or B-myr | other 78 | 3 among the three | 39 (exactly one) | 81 | raw, 120 |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for arm in v['arms']:
        cells = [' / '.join(_f(_stat(dt, arm, h, st).get('d'), '+.4f') for h in HOSTS)
                 for st in ('all_pattern', 'all_pattern_p42', 'all_pattern_p78', 'all_pattern_p3', 'all_pattern_p39', 'all_pattern_p81')]
        cells.append(' / '.join(_f(_stat(dt, arm, h, 'all_pattern', 'raw').get('d'), '+.4f') for h in HOSTS))
        L.append(f'| {arm} | ' + ' | '.join(cells) + ' |')
    L += ['', f'## The four wirings\' d, calibrated pattern part{P}', '',
          'SD of the four: the d\'s over the wirings. SD of the re-draws: the host\'s seven re-draws\' pattern parts in step 1.', '',
          '| arm | host | seed 0 | seed 1 | seed 2 | seed 3 | d | SD of the four | SD of the re-draws | d on halves A / B / C / D |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for arm in v['arms']:
        for h in HOSTS:
            x = _stat(dt, arm, h, 'all_pattern')
            hv = v['arms'][arm]['hosts'][h]['d_halves']
            L.append(f"| {arm} | {E.SHORT[h]} | " + ' | '.join(_f(x.get(f'd_seed{s}'), '+.4f') for s in SEEDS)
                     + f" | {_f(x.get('d'), '+.4f')} | {_f(x.get('sd_wirings'))} | {_f(x.get('sd_redraws'))} | "
                     + ' / '.join(_f(hv.get(t), '+.4f') for t in HALVES) + ' |')
    a1 = rd1.drop_duplicates(['calyx', 'seed', 'readout'], keep='last').set_index(['calyx', 'seed', 'readout'])
    a2 = rd2.set_index(['tag', 'readout']) if len(rd2) else None
    ri = runs2.set_index('tag') if len(runs2) else None
    for h in HOSTS:
        L += ['', f'## Runs on the {E.SHORT[h]}: readings, calibrated, all recordings{P}', '',
              '| run | pattern part | d | amplitude | shape | spread ratio | Pearson r | mean APL effect | 42 | 78 | 3 | raw pattern part |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for s in SEEDS:
            if (h, s, 'cal') in a1.index:
                r, rr = a1.loc[(h, s, 'cal')], a1.loc[(h, s, 'raw')]
                L.append(f"| own, seed {s} | {_f(r['all_pattern'])} | | {_f(r['all_amplitude'])} | {_f(r['all_shape'])} | "
                         f"{_f(r['all_sd_ratio'], '.2f')} | {_f(r['all_pearson_effect'], '.3f')} | {_f(r['mean_effect'], '+.3f')} | "
                         f"{_f(r['all_pattern_p42'], '.3f')} | {_f(r['all_pattern_p78'], '.3f')} | {_f(r['all_pattern_p3'], '.3f')} | "
                         f"{_f(rr['all_pattern'])} |")
        for arm in v['arms']:
            for s in SEEDS:
                tag = E.tag_of(E.Point(2, h, s, arm))
                if a2 is None or (tag, 'cal') not in a2.index:
                    continue
                r, rr = a2.loc[(tag, 'cal')], a2.loc[(tag, 'raw')]
                own = a1.loc[(h, s, 'cal')]['all_pattern'] if (h, s, 'cal') in a1.index else np.nan
                L.append(f"| {arm}, seed {s} | {_f(r['all_pattern'])} | {_f(float(r['all_pattern']) - float(own), '+.4f')} | "
                         f"{_f(r['all_amplitude'])} | {_f(r['all_shape'])} | {_f(r['all_sd_ratio'], '.2f')} | "
                         f"{_f(r['all_pearson_effect'], '.3f')} | {_f(r['mean_effect'], '+.3f')} | {_f(r['all_pattern_p42'], '.3f')} | "
                         f"{_f(r['all_pattern_p78'], '.3f')} | {_f(r['all_pattern_p3'], '.3f')} | {_f(rr['all_pattern'])} |")
        L += ['', f'## Runs on the {E.SHORT[h]}: realism{P}', '',
              'Own: the host wiring\'s own run. "x own": the held parameter over the own run\'s.', '',
              '| run | response rate on / off | off / on | γ off / on (own) | rate on γ / α′β′ / αβ | largest count on / off (own) '
              '| KCs with 20+ spikes on / off | threshold x own | APL->KC x own | APL->bouton x own | in bands / corner / guard |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
        for arm in v['arms']:
            for s in SEEDS:
                tag = E.tag_of(E.Point(2, h, s, arm))
                if ri is None or tag not in ri.index:
                    continue
                r = ri.loc[tag]
                L.append(f"| {arm}, seed {s} | {_f(r['rate_on'])} / {_f(r['rate_off'])} | {_f(r['ratio'], '.2f')} | "
                         f"{_f(r['ratio_g'], '.2f')} ({_f(r['ratio_g_own'], '.2f')}) | "
                         f"{_f(r['rate_on_g'], '.3f')} / {_f(r['rate_on_apbp'], '.3f')} / {_f(r['rate_on_ab'], '.3f')} | "
                         f"{_f(r['max_count_on'], '.0f')} / {_f(r['max_count_off'], '.0f')} ({_f(r['max_count_on_own'], '.0f')} / "
                         f"{_f(r['max_count_off_own'], '.0f')}) | {_f(r['n_kc_ge20_on'], '.0f')} / {_f(r['n_kc_ge20_off'], '.0f')} | "
                         f"{_f(r['fixed_thr_vs_own'], '.3f')} | {_f(r['wAPLKC_scale_vs_own'], '.3f')} | {_f(r['wAPLPN_scale_vs_own'], '.3f')} | "
                         f"{'yes' if r['in_bands'] else 'NO'} / {'YES' if r['corner'] else 'no'} / {'YES' if r['guard_max_count'] else 'no'} |")
    L += ['', f'## The transplants\' sizes (the runs\' rows, range over the wirings in){P}', '', '| arm | host | size |', '|---|---|---|']
    for arm in v['arms']:
        for h in HOSTS:
            x = runs2[(runs2['arm'] == arm) & (runs2['calyx'] == h)] if len(runs2) else runs2
            if len(x):
                L.append(f'| {arm} | {E.SHORT[h]} | {_size_cell(x, arm)} |')
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------------------
# the part

def score2(partial: bool = False, check_own: bool = True, out: Optional[Path] = None) -> Optional[dict]:
    """Step 2's score (the module docstring). Without all 48 runs: the missing points are printed and nothing is
    written (returns None), unless `partial`, which writes the provisional files."""
    out = E.out_dir() if out is None else Path(out)
    pts = E.planned_step2()
    outs, runs, missing, extra = collect2(pts)
    n_in, n_planned = len(outs), len(pts)
    print(f'E133 score2: {n_in} of {n_planned} step-2 runs in; other step-2 rows not scored: {extra or "none"}', flush=True)
    if missing:
        print('E133 score2: missing: ' + ', '.join(f"{m['tag']} ({m['why']})" for m in missing), flush=True)
        if not partial:
            print('E133 score2: no verdict until all 48 are in (--partial writes a provisional table)', flush=True)
            return None
    prov = bool(missing)
    rd1 = pd.read_csv(out / 'readings1.csv')
    runs1 = pd.read_csv(out / 'runs1.csv')
    oc = own_check(rd1) if check_own else dict(skipped=True)
    if check_own and not oc['passes'] and not prov:
        raise RuntimeError(f"E133 score2: the own runs' readings differ from readings1.csv ({oc}); run --part score1 again")
    rd2 = readings2(outs, runs) if len(runs) else pd.DataFrame(columns=['tag', 'arm', 'calyx', 'seed', 'readout'])
    r2 = runs2_table(runs, runs1) if len(runs) else runs
    dt = d_table(rd2, rd1)
    G = gap(rd1)
    nct = no_class_tags(r2, runs1)
    v = verdict2(dt, G, nct, prov)
    v.update(n_runs_in=n_in, n_planned=n_planned, missing=missing, other_step2_rows=extra,
             no_class_runs=sorted({t for ts in nct.values() for t in ts}), realism_flags=realism_flags(r2), own_check=oc,
             inputs=dict(readings1_sha256=_sha256_file(out / 'readings1.csv'), runs1_sha256=_sha256_file(out / 'runs1.csv'),
                         runs2_files=sorted(p.name for p in out.glob('runs_2_*.csv'))),
             generated=E._now(), host=E._host(), script=SCRIPT, commits=H.code_commits())
    rd2.to_csv(out / out_name('readings2.csv', prov), index=False)
    r2.to_csv(out / out_name('runs2.csv', prov), index=False)
    dt.to_csv(out / out_name('d2.csv', prov), index=False)
    E._write_json(out / out_name('verdict2.json', prov), v)
    (out / out_name('step2_table.md', prov)).write_text(step2_table(v, dt, rd2, rd1, r2, missing, n_in, n_planned))
    if not prov:
        for name in OUTPUTS:
            p = out / out_name(name, True)
            if p.exists():
                p.unlink()
    print(json.dumps(dict(provisional=prov, G=v['G']['all']['G'], step3=v['step3'],
                          arms={a: r.get('cls') or r.get('provisional_reading') or f"no class, would be {r.get('would_be')}"
                                for a, r in v['arms'].items()}), default=str), flush=True)
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--partial', action='store_true',
                    help="with runs missing, write a provisional table (*_provisional.*; no class, no step-3 list)")
    ap.add_argument('--no-own-check', action='store_true', help="skip reading the hosts' own runs again (own_check)")
    a = ap.parse_args()
    v = score2(partial=a.partial, check_own=not a.no_own_check)
    if v is None:
        sys.exit(2)


if __name__ == '__main__':
    main()
