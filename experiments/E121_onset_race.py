#!/usr/bin/env python3
"""E121 step 0: the onset race, read from the family's spikes and APL's trace (no new model code).

The family is the breadth-2 APL-effect fit on its eight wirings (E72): S = the hemibrain wiring (C2, E67's
rate0.068_aplboth_k2_cap3), P1, P2, M1; T = P3, P4, M2, M3. For each wiring the fitted model is rebuilt from its
rebuild archive (pyolfsysm.rebuild, data/model_outputs/<exp>/rebuild/<tag>/) at the default 0.5 s pulse, APL on and
APL off (acute), with olfsysm's spike recordings and APL trace (`inh_sims`) saved: one `fit_mb_model` call per
condition with `return_olfsysm_vars=True`, then the KC layer run again with `mp.kc.save_spike_recordings` and
`mp.kc.save_inh_sims` set (committed olfsysm, no option added). The rerun must reproduce the archive's saved spike
counts and responses exactly, else the run stops.

Per odor and wiring: APL's rise, the time from odor onset to 90 % of its plateau (`apl_rise`; the plateau is the
mean of the trace over the last 100 ms of the pulse; a second reading takes the trace's maximum over the pulse).
Per responding KC and odor (2 or more spikes) the time of its second spike, the response's completion under the
2-spike criterion. A response "completes before APL's rise" when its second spike is at or before APL's 90 % point
of that odor and wiring (APL on's trace in both conditions: APL off is the reference, read against the same clock).

The gate (E121's entry): fewer than 20 % of APL-on responses complete before APL's 90 % point: no race, and the
experiment closes without code.  The family's APL trace overshoots (it peaks 50 to 100 ms after onset and sags to a lower
sustained level by the pulse's end), so "its plateau" has two readings, and the gate is read on both: when they
disagree the verdict is "split" (race None in step0_gate.json) and goes to the planner. Beside them, a reading that
needs no plateau: APL's level at each response's completion as a fraction of its peak and of its sustained level
(`apl_level_at`).

Groups (E77_input_screen's): ester group 2h, IaA, pa, eb, ep, ms; alcohol group 1-8ol, 1-5ol, 1-6ol, 6al; the rest
(2-but, t2h, paa, -aPine, euc, B-myr) "other". Breadth: the number of the 16 odors the KC responds to in that
condition (bins 1, 2-3, 4-7, 8+).

Outputs (results/E121_onset_race/): step0_shares.csv (share complete before APL's 90 % point, overall, by group,
by breadth, by odor, per wiring and pooled, APL on and off, both plateau readings), step0_apl_rise.csv (per wiring
and odor: plateau, time to 50 % and 90 %, both readings), step0_gate.json, provenance.json. Per-wiring spike times
and APL traces under outputs/E121_onset_race/step0/.

Part `reads` (no simulation; added 2026-11-23 for the note): the reads a Fable planner made inline from step 0's saved
traces and spike times when it read the split gate (claude-fable-5-2, spawned from viki 2026-11-23 21:37), by its
stated definitions, so the note quotes committed code. APL's levels count its activity before the odor (a KC meets
APL's whole activity, not its rise above rest): per wiring and odor, on APL on's trace, rest (the trace at the last step
at or before 0 s), peak (its maximum over the pulse) and its time, and sustained (its mean over the last LATE_WINDOW s
of the pulse); a response's level at completion is the trace at the last step at or before its second spike, read as
a fraction of the peak and of the sustained level. Survival: among the KC x odor entries responding with APL off, the
survivors respond with APL on too; the APL-off second spike's time and spike count are scored as tellers of survival by
their AUC (Mann-Whitney). Windows: per condition, responders with two or more spikes in [0, 0.1), [0.25, 0.5) and after
the pulse (half-open windows), the shares of spikes in the first 100 ms and after the pulse, and spike counts. Beside
them the family's per-odor response rates relative to the mean odor's against the control flies', the family's realism
table (results/apl_effect/by_variant.csv), and the values the planner's reasons cite from other experiments' tables.
Writes reads_*.csv, reads_summary.json (with the planner's check values against the code's) and provenance_reads.json.
The saved spike times are checked against the archives' spike counts and step 0's hashes first.

Two reads added after the placement (2026-11-23; the placement's Fable planner, claude-fable-5-2, spawned from viki
2026-11-23 22:23, made them inline from the model-output archive), in reads_additions.csv (one row per wiring) and
reads_summary.json's 'additions': (1) the APL-on responders that do not respond with APL off, in the family's rebuild
archives and in E110's arm A (KC-KC at 0 at fixed scales, both response counts held; data/model_outputs/E110_no_kckc/A/,
variant A_real for S); (2) among the APL-off responders, a KC's summed APL->KC weight (wAPLKC.parquet of the rebuild
archive, summed over the KC's claws) and the archive's breadth value (breadth.parquet) as tellers of survival under
APL, lower favouring survival, by their AUC, with the survivors' and the removed's median weight over the mean KC's.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E121_onset_race.py step0
    ./venv/bin/python experiments/E121_onset_race.py reads
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

EXPERIMENT = 'E121_onset_race'
GATE_SHARE = 0.20
RISE_FRACTION = 0.9
# (the plateau of APL's trace: 'late', its mean over the last LATE_WINDOW s of the pulse, the primary reading; 'max',
# its maximum over the pulse)
PLATEAUS = ('late', 'max')
LATE_WINDOW = 0.1
STIM_END = 0.5
ESTER = ('2h', 'IaA', 'pa', 'eb', 'ep', 'ms')
ALCOHOL = ('1-8ol', '1-5ol', '1-6ol', '6al')
BREADTH_BINS = ((1, 1, '1'), (2, 3, '2-3'), (4, 7, '4-7'), (8, 16, '8+'))
S_SET = ('S', 'P1', 'P2', 'M1')
T_SET = ('P3', 'P4', 'M2', 'M3')
TAG = 'rate0.068_aplboth_k2_cap3'
# wiring label -> (archive experiment folder, rebuild tag)
WIRINGS = {'S': ('E67_apl_effect_fits', TAG)}
WIRINGS.update({f'P{i}': ('E72_wiring_population', f'wiring-claw_permutation_within_type-seed{i}_{TAG}')
                for i in range(1, 5)})
WIRINGS.update({f'M{i}': ('E72_wiring_population', f'wiring-curveball_within_type_matched-seed{i}_{TAG}')
                for i in range(1, 4)})
CONDITIONS = ('on', 'off')
# (the reads: half-open windows [lo, hi) in s from onset; a response completing before EARLY_T2 counts as early)
READ_WINDOWS = (('first100', 0.0, 0.1), ('last250', 0.25, STIM_END), ('after', STIM_END, float('inf')))
EARLY_T2 = 0.045
# (absorbs the float difference between a spike's time and the trace's step times, which are rounded to 0.1 ms)
TIME_TOL = 1e-9
LEVEL_BYS = (('overall', None), ('group', 'group'), ('breadth', 'breadth_bin'), ('odor', 'odor_name'))


def odor_name(col) -> str:
    """'2h @ -3' -> '2h'."""
    return str(col).split(' @ ')[0]


def odor_group(odor: str) -> str:
    """The odor's group: 'ester' (E77's ester group), 'alcohol' (E77's alcohol group) or 'other'."""
    return 'ester' if odor in ESTER else 'alcohol' if odor in ALCOHOL else 'other'


def breadth_bin(n: int) -> str:
    """The bin of a KC's number of odors responded to (1, 2-3, 4-7, 8+)."""
    for lo, hi, lab in BREADTH_BINS:
        if lo <= n <= hi:
            return lab
    raise ValueError(f'breadth {n} outside 1-16')


def apl_rise(trace: np.ndarray, time_s: np.ndarray, stim_end: float = STIM_END, plateau: str = 'late',
             fraction: float = RISE_FRACTION) -> Tuple[float, float]:
    """Returns (plateau, time from onset to `fraction` of it) of one APL trace.

    The trace is measured from its value at odor onset (the last step at or before 0 s), so a resting APL activity
    does not count as rise. The plateau is the trace's mean over the last LATE_WINDOW s of the pulse ('late') or
    its maximum over the pulse ('max'); the time is the first step after onset at which the trace reaches
    `fraction` of the plateau (NaN if the plateau is not above the onset value).

    Args:
        trace: APL activity per time step.
        time_s: the steps' times, seconds from odor onset.
        stim_end: the pulse's end, seconds from onset.
        plateau: 'late' or 'max'.
        fraction: of the plateau.
    """
    trace = np.asarray(trace, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    if trace.shape != time_s.shape:
        raise ValueError('trace and time_s differ in shape')
    pre = time_s <= 0
    base = trace[pre][-1] if pre.any() else 0.0
    during = (time_s > 0) & (time_s <= stim_end)
    if plateau == 'late':
        top = trace[during & (time_s > stim_end - LATE_WINDOW)].mean()
    elif plateau == 'max':
        top = trace[during].max()
    else:
        raise ValueError(f'plateau {plateau!r} not in {PLATEAUS}')
    amp = top - base
    if not amp > 0:
        return float(amp), float('nan')
    hit = np.nonzero(during & (trace - base >= fraction * amp))[0]
    return float(amp), (float(time_s[hit[0]]) if len(hit) else float('nan'))


def second_spike_times(spikes: pd.DataFrame, n: int = 2) -> pd.DataFrame:
    """Returns the `n`-th spike's time per (kc, odor) with `n` or more spikes.

    Args:
        spikes: one row per spike, columns 'kc', 'odor', 'time_s' (as `simulate(..., spike_times=True)`'s).

    Returns DataFrame with columns 'kc', 'odor', 't_complete' (the n-th spike's time) and 'n_spikes'.
    """
    s = spikes.sort_values(['odor', 'kc', 'time_s'])
    g = s.groupby(['odor', 'kc'], sort=False)
    rank = g.cumcount()
    nsp = g['time_s'].transform('size')
    out = s.loc[rank.to_numpy() == n - 1, ['kc', 'odor', 'time_s']].rename(columns={'time_s': 't_complete'})
    out['n_spikes'] = nsp[rank.to_numpy() == n - 1].to_numpy()
    return out.reset_index(drop=True)


def completion_table(spikes: pd.DataFrame, odors: Sequence[str], n_kcs: int) -> pd.DataFrame:
    """Per response (KC x odor with 2 or more spikes): its completion time, odor, group and the KC's breadth."""
    c = second_spike_times(spikes)
    c['odor_name'] = [odors[i] for i in c['odor']]
    c['group'] = c['odor_name'].map(odor_group)
    breadth = np.bincount(c['kc'].to_numpy(), minlength=n_kcs)
    c['breadth'] = breadth[c['kc'].to_numpy()]
    c['breadth_bin'] = [breadth_bin(int(b)) for b in c['breadth']]
    return c


def mark_before(completion: pd.DataFrame, rise: pd.DataFrame, plateau: str) -> pd.DataFrame:
    """Adds 't_apl90' (APL's 90 % point of the response's odor, from `rise`, rows of one wiring with columns
    'odor_name', 'plateau', 't90') and 'before' (the response completes at or before it)."""
    r = rise[rise['plateau'] == plateau].set_index('odor_name')['t90']
    c = completion.copy()
    c['t_apl90'] = c['odor_name'].map(r).to_numpy()
    c['before'] = c['t_complete'] <= c['t_apl90']
    return c


def shares(marked: pd.DataFrame) -> pd.DataFrame:
    """Returns the shares of responses complete before APL's 90 % point: overall, by group, by breadth, by odor.

    Args:
        marked: `mark_before`'s rows (one wiring, condition and plateau reading, or several pooled).

    Returns DataFrame with columns 'by', 'level', 'n_responses', 'n_before', 'share', 'median_t_complete',
    'median_t_apl90'.
    """
    rows = []
    for by, col in (('overall', None), ('group', 'group'), ('breadth', 'breadth_bin'), ('odor', 'odor_name')):
        groups = [('all', marked)] if col is None else list(marked.groupby(col, sort=True))
        for level, g in groups:
            n = len(g)
            rows.append(dict(by=by, level=str(level), n_responses=n, n_before=int(g['before'].sum()),
                             share=(float(g['before'].mean()) if n else float('nan')),
                             median_t_complete=float(g['t_complete'].median()) if n else float('nan'),
                             median_t_apl90=float(g['t_apl90'].median()) if n else float('nan'),
                             **{f'median_{k}': (float(g[k].median()) if (n and k in g) else float('nan'))
                                for k in ('apl_frac_peak', 'apl_frac_late')}))
    return pd.DataFrame(rows)


def apl_level_at(completion: pd.DataFrame, apl: pd.DataFrame) -> pd.DataFrame:
    """Adds APL's level at each response's completion: 'apl_frac_peak' (APL's trace above its onset value at
    `t_complete`, as a fraction of its maximum over the pulse) and 'apl_frac_late' (the same, as a fraction of its
    mean over the last LATE_WINDOW s of the pulse); a reading of the race that needs no choice of plateau.

    Args:
        completion: `completion_table`'s rows (one wiring and condition).
        apl: odors x time steps, the APL trace of the same wiring (APL on), columns the times from onset.
    """
    t = apl.columns.to_numpy(dtype=float)
    pre, during = t <= 0, (t > 0) & (t <= STIM_END)
    late = during & (t > STIM_END - LATE_WINDOW)
    c = completion.copy()
    fp = np.full(len(c), np.nan)
    fl = np.full(len(c), np.nan)
    for odor, idx in c.groupby('odor_name').indices.items():
        tr = apl.loc[odor].to_numpy(dtype=float)
        base = tr[pre][-1]
        lvl = np.interp(c['t_complete'].to_numpy()[idx], t, tr) - base
        fp[idx] = lvl / (tr[during].max() - base)
        fl[idx] = lvl / (tr[late].mean() - base)
    c['apl_frac_peak'], c['apl_frac_late'] = fp, fl
    return c


def run_wiring(label: str, archive: Path, out: Path) -> Dict[str, object]:
    """Rebuilds one wiring's fit, APL on and off, and returns its spikes and APL traces.

    One `fit_mb_model` call per condition at the archive's parameters (`pyolfsysm.rebuild.simulate_at`, with
    `multibouton_claws='first'` pinned as the family was fitted, and `return_olfsysm_vars`), then olfsysm's KC
    layer run again with the spike recordings and APL's trace saved. Both runs must reproduce the archive's
    saved spike counts and responses exactly (RuntimeError otherwise).

    Returns dict with, per condition c in CONDITIONS: 'spikes_<c>' (one row per spike: kc, odor, step, time_s),
    'apl_<c>' (odors x time steps, APL's trace, columns the times from onset), 'counts_<c>'; and 'odors' (names, in
    output order), 'n_kcs', 'check' (the reproduction checks).
    """
    import olfsysm as osm
    from pyolfsysm import rebuild as R
    from pyolfsysm.mb_model import _spike_times_frame, get_sim_var

    a = R.load_archive(archive)
    res: Dict[str, object] = {'check': {}}
    for c in CONDITIONS:
        kws = dict(a['model_kws'], return_olfsysm_vars=True, multibouton_claws='first')
        r, s, _, pdict = R.simulate_at(a['params'], a['wiring'], a['weights'], a['orn_deltas'], model_kws=kws,
                                       apl_off=(c == 'off'))
        mp, rv = pdict['mp'], pdict['rv']
        rs, ss = a['outputs'][c]
        model_odors = a['orn_deltas'].columns
        if 'panel' in model_odors.names:
            model_odors = model_odors.droplevel('panel')
        odor_order = model_odors.get_indexer(r.columns)
        assert (odor_order >= 0).all()
        mp.kc.save_spike_recordings = True
        mp.kc.save_inh_sims = True
        osm.run_KC_sims(mp, rv, False)
        s2 = pd.DataFrame(np.asarray(rv.kc.spike_counts)[:, odor_order], index=s.index, columns=s.columns)
        ok = dict(fit_counts=bool(s.equals(ss)), fit_responses=bool(r.equals(rs)),
                  rerun_counts=bool(np.array_equal(s2.to_numpy(), ss.to_numpy())))
        res['check'][c] = ok
        if not all(ok.values()):
            raise RuntimeError(f'{label} {c}: the rebuild does not reproduce the archive ({ok})')
        spikes = _spike_times_frame(mp, rv, odor_order, s2.to_numpy())
        inh = get_sim_var(rv.kc, 'inh_sims')
        inh = np.asarray(inh, dtype=float).reshape(inh.shape[0], -1)[odor_order]
        t = mp.time_pre_start + np.arange(inh.shape[1]) * mp.time_dt
        if not np.isclose(mp.time_stim_end, STIM_END):
            raise RuntimeError(f'pulse ends at {mp.time_stim_end}, not the default {STIM_END}')
        res[f'spikes_{c}'] = spikes
        res[f'apl_{c}'] = pd.DataFrame(inh, index=[odor_name(o) for o in r.columns], columns=np.round(t, 6))
        res[f'counts_{c}'] = s2
        res['odors'] = [odor_name(o) for o in r.columns]
        res['n_kcs'] = int(s.shape[0])
        res['apl_taum'] = float(mp.kc.apl_taum)
        del mp, rv, pdict
    return res


def rise_table(label: str, apl: pd.DataFrame) -> pd.DataFrame:
    """Per odor and plateau reading: APL's rise (plateau above onset, times to 50 % and 90 %)."""
    t = apl.columns.to_numpy(dtype=float)
    rows = []
    for odor, tr in apl.iterrows():
        for pl in PLATEAUS:
            amp, t90 = apl_rise(tr.to_numpy(), t, plateau=pl)
            _, t50 = apl_rise(tr.to_numpy(), t, plateau=pl, fraction=0.5)
            rows.append(dict(wiring=label, odor_name=odor, group=odor_group(odor), plateau=pl, amplitude=amp,
                             t50=t50, t90=t90))
    return pd.DataFrame(rows)


def _sha(df: pd.DataFrame) -> str:
    return hashlib.sha256(pd.util.hash_pandas_object(df, index=True).to_numpy().tobytes()).hexdigest()


def step0(out_dir: Optional[Path] = None, res_dir: Optional[Path] = None,
          wirings: Optional[Sequence[str]] = None) -> dict:
    """Runs step 0 over the family's wirings and writes its tables; returns the gate dict."""
    from kcq import harness as H
    out_dir = out_dir or (H.outputs_dir(EXPERIMENT) / 'step0')
    res_dir = res_dir or H.results_dir(EXPERIMENT)
    out_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)
    archive_root = H.shared_dir() / 'data' / 'model_outputs'
    labels = list(wirings or WIRINGS)
    marked: List[pd.DataFrame] = []
    rises: List[pd.DataFrame] = []
    checks, hashes = {}, {}
    taum = None
    for lab in labels:
        exp, tag = WIRINGS[lab]
        arch = archive_root / exp / 'rebuild' / tag
        print(f'{lab}: {arch}', flush=True)
        o = run_wiring(lab, arch, out_dir)
        checks[lab] = o['check']
        taum = o['apl_taum']
        rise = rise_table(lab, o['apl_on'])
        rises.append(rise)
        for c in CONDITIONS:
            o[f'spikes_{c}'].to_parquet(out_dir / f'{lab}_{c}_spike_times.parquet')
            ap = o[f'apl_{c}'].copy()
            ap.columns = [f'{x:.4f}' for x in ap.columns]
            ap.to_parquet(out_dir / f'{lab}_{c}_apl_trace.parquet')
            hashes[f'{lab}_{c}_spike_times'] = _sha(o[f'spikes_{c}'])
            comp = apl_level_at(completion_table(o[f'spikes_{c}'], o['odors'], o['n_kcs']), o['apl_on'])
            for pl in PLATEAUS:
                # (APL on's clock in both conditions: APL off is read against the rise of the intact APL)
                m = mark_before(comp, rise, pl)
                m['wiring'], m['condition'], m['plateau'] = lab, c, pl
                marked.append(m)
        print(f'  done, checks {o["check"]}', flush=True)
    M = pd.concat(marked, ignore_index=True)
    rise_all = pd.concat(rises, ignore_index=True)
    rows = []
    sets = {lab: [lab] for lab in labels}
    sets.update({'S_set': [x for x in S_SET if x in labels], 'T_set': [x for x in T_SET if x in labels],
                 'pooled': labels})
    for name, members in sets.items():
        if not members:
            continue
        for (c, pl), g in M[M['wiring'].isin(members)].groupby(['condition', 'plateau'], sort=False):
            sh = shares(g)
            sh.insert(0, 'plateau', pl)
            sh.insert(0, 'condition', c)
            sh.insert(0, 'wiring', name)
            rows.append(sh)
    S = pd.concat(rows, ignore_index=True)
    S.to_csv(res_dir / 'step0_shares.csv', index=False)
    rise_all.to_csv(res_dir / 'step0_apl_rise.csv', index=False)

    def get(w, c, pl, by, level):
        r = S[(S.wiring == w) & (S.condition == c) & (S.plateau == pl) & (S.by == by) & (S.level == level)]
        return float(r['share'].iloc[0]) if len(r) else float('nan')

    per_w = {lab: get(lab, 'on', 'late', 'overall', 'all') for lab in labels}
    pooled = get('pooled', 'on', 'late', 'overall', 'all')
    pooled_max = get('pooled', 'on', 'max', 'overall', 'all')
    race_late, race_max = pooled >= GATE_SHARE, pooled_max >= GATE_SHARE
    race = race_late if race_late == race_max else None
    on_l = M[(M.condition == 'on') & (M.plateau == 'late')]
    rl = rise_all[rise_all.plateau == 'late']
    gate = dict(
        experiment=EXPERIMENT, step='step0',
        rule=(f'fewer than {GATE_SHARE:.0%} of APL-on responses complete (second spike) before APL reaches '
              f'{RISE_FRACTION:.0%} of its plateau: no race, the experiment closes without code'),
        verdict=('race: step 0 passes (the share is at or above the gate); steps 1-3 wait on the olfsysm option'
                 if race else 'no race: the experiment closes without code' if race is not None else
                 ("split: the two readings of APL's plateau disagree (APL's trace overshoots, peaking early in the "
                  "pulse and sagging to a lower sustained level); "
                  f"by the sustained level {'race' if race_late else 'no race'}, by the peak "
                  f"{'race' if race_max else 'no race'}: for the planner")),
        race=race, race_plateau_late=bool(race_late), race_plateau_max=bool(race_max),
        apl_frac_peak_at_completion_on=dict(median=float(on_l['apl_frac_peak'].median()),
            q25=float(on_l['apl_frac_peak'].quantile(0.25)), q75=float(on_l['apl_frac_peak'].quantile(0.75))),
        apl_frac_late_at_completion_on=dict(median=float(on_l['apl_frac_late'].median()),
            q25=float(on_l['apl_frac_late'].quantile(0.25)), q75=float(on_l['apl_frac_late'].quantile(0.75))),
        apl_peak_over_late=float(np.median(
            rise_all[rise_all.plateau == 'max'].set_index(['wiring', 'odor_name'])['amplitude'] /
            rise_all[rise_all.plateau == 'late'].set_index(['wiring', 'odor_name'])['amplitude'])),
        apl_t90_peak_s=dict(median=float(rise_all[rise_all.plateau == 'max']['t90'].median()),
            min=float(rise_all[rise_all.plateau == 'max']['t90'].min()),
            max=float(rise_all[rise_all.plateau == 'max']['t90'].max())),
        share_on_by_group_plateau_max={g: get('pooled', 'on', 'max', 'group', g) for g in ('ester', 'alcohol', 'other')},
        share_on_pooled=pooled, share_on_pooled_plateau_max=pooled_max,

        share_on_per_wiring=per_w,
        share_on_range_over_wirings=[min(per_w.values()), max(per_w.values())],
        share_on_by_group={g: get('pooled', 'on', 'late', 'group', g) for g in ('ester', 'alcohol', 'other')},
        share_on_by_breadth={b: get('pooled', 'on', 'late', 'breadth', b) for _, _, b in BREADTH_BINS},
        share_off_pooled=get('pooled', 'off', 'late', 'overall', 'all'),
        share_off_by_group={g: get('pooled', 'off', 'late', 'group', g) for g in ('ester', 'alcohol', 'other')},
        apl_t90_s=dict(median=float(rl['t90'].median()), min=float(rl['t90'].min()), max=float(rl['t90'].max()),
                       by_group={g: float(rl[rl.group == g]['t90'].median()) for g in ('ester', 'alcohol', 'other')}),
        median_t_complete_on_s=float(M[(M.condition == 'on') & (M.plateau == 'late')]['t_complete'].median()),
        n_responses_on=int(((M.condition == 'on') & (M.plateau == 'late')).sum()),
        apl_taum_s=taum, plateau_definition=(f"'late': the mean of APL's trace (inh_sims, measured from its value at "
            f"onset) over the last {LATE_WINDOW} s of the {STIM_END} s pulse (primary); 'max': its maximum over the pulse"),
        groups=dict(ester=list(ESTER), alcohol=list(ALCOHOL),
                    other='the rest of the 16 APL16 odors (2-but, t2h, paa, -aPine, euc, B-myr)'),
        breadth='the number of the 16 odors the KC responds to (2 or more spikes) in that condition',
        wirings={lab: '/'.join(WIRINGS[lab]) for lab in labels},
        reproduction_checks=checks,
        approx_conc_note=getattr(H, 'APPROX_CONC_NOTE', None),
        time=datetime.now(timezone.utc).isoformat(timespec='seconds'),
    )
    (res_dir / 'step0_gate.json').write_text(json.dumps(gate, indent=1, default=str) + '\n')
    prov = dict(script='experiments/E121_onset_race.py', part='step0', time=gate['time'],
                model_provenance=H.model_provenance(), code_commits=H.code_commits(),
                archives={lab: str(archive_root / WIRINGS[lab][0] / 'rebuild' / WIRINGS[lab][1]) for lab in labels},
                spike_times_sha256=hashes, outputs_dir=str(out_dir),
                how=('each archive rebuilt with pyolfsysm.rebuild.simulate_at (multibouton_claws first, '
                     'return_olfsysm_vars), then olfsysm.run_KC_sims with save_spike_recordings and save_inh_sims; '
                     'spike counts and responses checked equal to the archive'))
    (res_dir / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    return gate


# --- part reads: the planner's reads of step 0's saved files (no simulation) ---

def apl_levels(trace: np.ndarray, time_s: np.ndarray, stim_end: float = STIM_END,
               late_window: float = LATE_WINDOW) -> Dict[str, float]:
    """Returns APL's levels of one trace, each including its level before the odor (not measured from it).

    Args:
        trace: APL activity per time step.
        time_s: the steps' times, seconds from odor onset (increasing).
        stim_end: the pulse's end, seconds from onset.
        late_window: the sustained level is the mean over the pulse's last `late_window` s.

    Returns dict with 'rest' (the trace at the last step at or before 0 s), 'peak' (its maximum over the pulse,
    0 < t <= stim_end, as `apl_rise`'s), 't_peak' (the time of the first such maximum), 'sustained' (its mean over
    stim_end - late_window < t <= stim_end, `apl_rise`'s 'late' plateau), 'rest_over_peak' and 'peak_over_sustained'.
    """
    trace = np.asarray(trace, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    if trace.shape != time_s.shape:
        raise ValueError('trace and time_s differ in shape')
    pre = np.nonzero(time_s <= 0)[0]
    if not len(pre):
        raise ValueError('no step at or before onset')
    during = np.nonzero((time_s > 0) & (time_s <= stim_end))[0]
    late = (time_s > stim_end - late_window) & (time_s <= stim_end)
    i_pk = during[int(np.argmax(trace[during]))]
    rest, peak, sus = float(trace[pre[-1]]), float(trace[i_pk]), float(trace[late].mean())
    return dict(rest=rest, peak=peak, t_peak=float(time_s[i_pk]), sustained=sus, rest_over_peak=rest / peak,
                peak_over_sustained=peak / sus)


def level_at(trace: np.ndarray, time_s: np.ndarray, t: Sequence[float]) -> np.ndarray:
    """Returns the trace at the last step at or before each time in `t` (NaN before the first step).

    A spike's time equals a step's time up to rounding; TIME_TOL absorbs that, so a spike at a step reads that step.
    """
    trace = np.asarray(trace, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    idx = np.searchsorted(time_s, np.asarray(t, dtype=float) + TIME_TOL, side='right') - 1
    out = np.full(idx.shape, np.nan)
    ok = idx >= 0
    out[ok] = trace[idx[ok]]
    return out


def auc(pos: Sequence[float], neg: Sequence[float]) -> float:
    """Returns the AUC of telling `pos` from `neg` by larger values: P(a pos value > a neg value) plus half the
    probability of a tie (the Mann-Whitney U over the product of the sizes). Negate both for 'smaller favours pos'."""
    pos, neg = np.asarray(pos, dtype=float), np.asarray(neg, dtype=float)
    if not (len(pos) and len(neg)):
        return float('nan')
    r = pd.Series(np.concatenate([pos, neg])).rank(method='average').to_numpy()
    n1 = len(pos)
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * len(neg)))


def window_counts(spikes: pd.DataFrame, lo: float, hi: float) -> pd.Series:
    """Returns the number of spikes in [lo, hi) per (kc, odor) with any there (index kc, odor)."""
    m = (spikes['time_s'] >= lo - TIME_TOL) & (spikes['time_s'] < hi - TIME_TOL)
    return spikes[m].groupby(['kc', 'odor']).size()


def window_reads(spikes: pd.DataFrame) -> Dict[str, float]:
    """Returns one wiring and condition's firing in READ_WINDOWS (from its spikes, one row per spike).

    A responder has two or more spikes over the trial. Returns dict with 'n_responders', 'n_<window>' (KC x odor with
    two or more spikes in the window), 'share_responders_first100' and '_last250' (n_<window> / n_responders), 'share_spikes_<window>'
    (of all spikes), 'share_exactly_two' (of responders), 'median_spikes' and 'max_spikes' (per responder; the
    largest count per KC and odor), 'n_spikes', and the second spikes' 't2_p05', 't2_median', 't2_p95' (s).
    """
    n = spikes.groupby(['kc', 'odor']).size()
    resp = n[n >= 2]
    out: Dict[str, float] = dict(n_responders=int(len(resp)), n_spikes=int(len(spikes)))
    t = spikes['time_s']
    for name, lo, hi in READ_WINDOWS:
        out[f'n_{name}'] = int((window_counts(spikes, lo, hi) >= 2).sum())
        out[f'share_spikes_{name}'] = float(((t >= lo - TIME_TOL) & (t < hi - TIME_TOL)).mean())
    for name in ('first100', 'last250'):
        out[f'share_responders_{name}'] = out[f'n_{name}'] / out['n_responders'] if len(resp) else float('nan')
    out['share_exactly_two'] = float((resp == 2).mean())
    out['median_spikes'] = float(resp.median())
    out['max_spikes'] = int(n.max())
    t2 = second_spike_times(spikes)['t_complete']
    out.update(t2_p05=float(t2.quantile(0.05)), t2_median=float(t2.median()), t2_p95=float(t2.quantile(0.95)))
    return out


def levels_at_completion(completion: pd.DataFrame, apl: pd.DataFrame, levels: pd.DataFrame) -> pd.DataFrame:
    """Adds APL's level at each response's completion (its second spike) to `completion_table`'s rows.

    Args:
        completion: one wiring and condition's `completion_table` rows.
        apl: odors x time steps, APL on's trace of the wiring (columns the times from onset).
        levels: `apl_levels` per odor of the same trace (index odor_name).

    Adds 'apl_level' (the trace at the last step at or before the second spike), 'frac_peak', 'frac_sustained' (that
    level over the odor's peak and sustained level, both counting APL's level before the odor), 'below_sustained'
    and 'before_peak' (the second spike before the peak's time).
    """
    t = apl.columns.to_numpy(dtype=float)
    c = completion.copy()
    lvl = np.full(len(c), np.nan)
    for odor, idx in c.groupby('odor_name').indices.items():
        lvl[idx] = level_at(apl.loc[odor].to_numpy(dtype=float), t, c['t_complete'].to_numpy()[idx])
    lv = levels.loc[c['odor_name']]
    c['apl_level'] = lvl
    c['frac_peak'] = lvl / lv['peak'].to_numpy()
    c['frac_sustained'] = lvl / lv['sustained'].to_numpy()
    c['below_sustained'] = c['frac_sustained'] < 1
    c['before_peak'] = c['t_complete'].to_numpy() < lv['t_peak'].to_numpy() - TIME_TOL
    return c


def level_summary(rows: pd.DataFrame) -> pd.DataFrame:
    """Returns the level at completion over `levels_at_completion` rows: overall, by group, breadth and odor.

    Columns 'by', 'level', 'n', 'frac_peak_median', '_q25', '_q75', 'frac_sustained_median', '_q25', '_q75',
    'share_below_sustained', 'share_before_peak'.
    """
    out = []
    for by, col in LEVEL_BYS:
        groups = [('all', rows)] if col is None else list(rows.groupby(col, sort=True))
        for level, g in groups:
            d = dict(by=by, level=str(level), n=int(len(g)))
            for k in ('frac_peak', 'frac_sustained'):
                d[f'{k}_median'] = float(g[k].median())
                d[f'{k}_q25'] = float(g[k].quantile(0.25))
                d[f'{k}_q75'] = float(g[k].quantile(0.75))
            d['share_below_sustained'] = float(g['below_sustained'].mean())
            d['share_before_peak'] = float(g['before_peak'].mean())
            out.append(d)
    return pd.DataFrame(out)


def survival_reads(comp_on: pd.DataFrame, comp_off: pd.DataFrame) -> Dict[str, float]:
    """Returns who survives APL among the APL-off responders of one wiring (`second_spike_times` rows per condition).

    A survivor responds to the odor with APL off and with APL on; the removed respond with APL off only. The APL-off
    second spike's time is scored as a teller of survival with earlier favouring survival ('auc_t2'), the APL-off spike
    count with more favouring it ('auc_count'), 0.5 being chance. Beside them, the APL-on responders that do not
    respond with APL off ('n_on_not_off', 'share_on_not_off' of the APL-on responders).
    """
    on = set(zip(comp_on['kc'].to_numpy(), comp_on['odor'].to_numpy()))
    off = set(zip(comp_off['kc'].to_numpy(), comp_off['odor'].to_numpy()))
    surv = np.array([k in on for k in zip(comp_off['kc'].to_numpy(), comp_off['odor'].to_numpy())], dtype=bool)
    t2, nsp = comp_off['t_complete'].to_numpy(dtype=float), comp_off['n_spikes'].to_numpy(dtype=float)
    n_on_not_off = sum(k not in off for k in on)
    return dict(n_off_responders=int(len(comp_off)), n_survivors=int(surv.sum()), n_removed=int((~surv).sum()),
                auc_t2=auc(-t2[surv], -t2[~surv]), auc_count=auc(nsp[surv], nsp[~surv]),
                median_t2_survivors=float(np.median(t2[surv])), median_t2_removed=float(np.median(t2[~surv])),
                share_early_survivors=float((t2[surv] < EARLY_T2 - TIME_TOL).mean()),
                share_early_removed=float((t2[~surv] < EARLY_T2 - TIME_TOL).mean()),
                n_on_responders=int(len(on)), n_on_not_off=int(n_on_not_off),
                share_on_not_off=float(n_on_not_off / len(on)) if on else float('nan'))


E110_ARM_A = ('E110_no_kckc', 'A')


def e110_arm_a_variant(label: str) -> str:
    """E110's arm-A variant of a wiring of the family ('A_real' for S, 'A_<label>' otherwise)."""
    return 'A_real' if label == 'S' else f'A_{label}'


def on_not_off_counts(resp_on: pd.DataFrame, resp_off: pd.DataFrame) -> Dict[str, float]:
    """Returns the APL-on responders that do not respond with APL off, from two KCs x odors response tables.

    Args:
        resp_on, resp_off: KCs x odors, 1 (or True) where the KC responds to the odor, same index and columns.

    Returns dict with 'n_on_responders', 'n_off_responders', 'n_on_not_off' and 'share_on_not_off' (of the APL-on
    responders).
    """
    if not (resp_on.index.equals(resp_off.index) and list(resp_on.columns) == list(resp_off.columns)):
        raise ValueError('the two response tables differ in their KCs or odors')
    on, off = resp_on.to_numpy() > 0, resp_off.to_numpy() > 0
    n_on, n = int(on.sum()), int((on & ~off).sum())
    return dict(n_on_responders=n_on, n_off_responders=int(off.sum()), n_on_not_off=n,
                share_on_not_off=float(n / n_on) if n_on else float('nan'))


def kc_apl_weight(w_apl_kc: pd.DataFrame, kc_ids: Sequence) -> np.ndarray:
    """Returns each KC's summed APL->KC weight, in the order of `kc_ids`.

    Args:
        w_apl_kc: the rebuild archive's wAPLKC (one row per claw, index levels including 'kc_id', one column).
        kc_ids: the KCs' IDs in the order of the spike counts' rows (a KC's position is its 'kc' in the spike times).
    """
    s = w_apl_kc.iloc[:, 0].groupby(level='kc_id').sum()
    out = s.reindex(pd.Index(kc_ids))
    if out.isna().any():
        raise ValueError(f'{int(out.isna().sum())} KCs without an APL->KC weight')
    return out.to_numpy(dtype=float)


def survival_by_weight(comp_on: pd.DataFrame, comp_off: pd.DataFrame, weight: np.ndarray,
                       breadth: np.ndarray) -> Dict[str, float]:
    """Returns how well a KC's APL->KC weight and its breadth value tell survival under APL (one wiring).

    As `survival_reads`: among the APL-off responders (`second_spike_times` rows), a survivor responds with APL on too.
    Each entry carries its KC's value (`weight`, `breadth`: per KC position); lower favours survival in both AUCs
    ('auc_apl_weight', 'auc_breadth'; 0.5 is chance). Beside them the survivors' and the removed's median weight over
    the mean KC's (all KCs of the wiring): 'weight_survivors_over_mean', 'weight_removed_over_mean'.
    """
    weight, breadth = np.asarray(weight, dtype=float), np.asarray(breadth, dtype=float)
    on = set(zip(comp_on['kc'].to_numpy(), comp_on['odor'].to_numpy()))
    kc = comp_off['kc'].to_numpy()
    surv = np.array([k in on for k in zip(kc, comp_off['odor'].to_numpy())], dtype=bool)
    w, b = weight[kc], breadth[kc]
    mean_w = float(weight.mean())
    return dict(auc_apl_weight=auc(-w[surv], -w[~surv]), auc_breadth=auc(-b[surv], -b[~surv]),
                weight_survivors_over_mean=float(np.median(w[surv]) / mean_w),
                weight_removed_over_mean=float(np.median(w[~surv]) / mean_w), mean_kc_weight=mean_w)


# The placement planner's values (claude-fable-5-2, spawned from viki 2026-11-23 22:23; read inline from the archive),
# checked as PLANNER_CHECKS are. The family's counts are per wiring in WIRINGS order.
PLACEMENT_CHECKS: Dict[str, Tuple[Tuple[float, ...], float, int]] = {
    'family_n_on_not_off_archive': ((373, 559, 357, 621, 510, 365, 312, 345), 1, 0),
    'kckc0_n_on_not_off_range': ((0, 0), 1, 0),
    'kckc0_n_on_responders_range': ((1852, 1931), 1, 0),
    'auc_apl_weight_range': ((0.71, 0.76), 1, 2),
    'auc_breadth_range': ((0.67, 0.71), 1, 2),
    'weight_survivors_over_mean_range': ((0.05, 0.22), 1, 2),
    'weight_removed_over_mean_range': ((0.86, 1.80), 1, 2),
}


# The planner's check values (2026-11-23 21:52; e121g_notebrief.txt), in the units given: name -> (values, scale,
# decimals). The code's values (times in s, shares as fractions) are multiplied by `scale` and agree when within half
# a unit of the last decimal.
PLANNER_CHECKS: Dict[str, Tuple[Tuple[float, ...], float, int]] = {
    'n_on_responses_pooled': ((15175,), 1, 0),
    'rest_over_peak_range_wirings': ((0.32, 0.33), 1, 2),
    'peak_over_sustained_range_wirings': ((1.55, 1.59), 1, 2),
    't_peak_range_ms': ((67, 93), 1000, 0),
    'frac_peak_median': ((0.80,), 1, 2),
    'frac_peak_iqr': ((0.74, 0.86), 1, 2),
    'frac_peak_median_range_wirings': ((0.79, 0.82), 1, 2),
    'frac_sustained_median': ((1.27,), 1, 2),
    'frac_sustained_iqr': ((1.16, 1.38), 1, 2),
    'below_sustained_pct': ((3.7,), 100, 1),
    'before_peak_pct': ((99.6,), 100, 1),
    'before_peak_pct_range_wirings': ((99.0, 100.0), 100, 1),
    'frac_peak_ester_alcohol_other': ((0.79, 0.81, 0.81), 1, 2),
    'frac_peak_breadth_1_8plus': ((0.84, 0.77), 1, 2),
    'frac_peak_odor_pa_bmyr_ms': ((0.77, 0.87, 0.87), 1, 2),
    'n_off_responders_range': ((3960, 4633), 1, 0),
    'n_survivors_range': ((1310, 1580), 1, 0),
    'auc_t2_range': ((0.60, 0.63), 1, 2),
    'auc_count_range': ((0.64, 0.68), 1, 2),
    'median_t2_survivors_ms_range': ((34, 36), 1000, 0),
    'median_t2_removed_ms_range': ((38, 40), 1000, 0),
    'early_removed_pct_range': ((69, 72), 100, 0),
    'early_survivors_pct_range': ((87, 91), 100, 0),
    'n_on_not_off_range': ((312, 621), 1, 0),
    'on_not_off_pct_range': ((16, 32), 100, 0),
    'on_n_responders_mean': ((1897,), 1, 0),
    'on_n_first100_mean': ((1891,), 1, 0),
    'on_n_last250_mean': ((96,), 1, 0),
    'on_n_last250_range': ((58, 152), 1, 0),
    'on_n_after_mean': ((0,), 1, 0),
    'on_spikes_first100_pct_mean': ((80,), 100, 0),
    'on_spikes_first100_pct_range': ((74, 86), 100, 0),
    'on_exactly_two_pct_mean': ((44,), 100, 0),
    'on_exactly_two_pct_range': ((40, 48), 100, 0),
    'on_median_spikes': ((3,), 1, 0),
    'off_n_responders_mean': ((4167,), 1, 0),
    'off_n_last250_mean': ((1113,), 1, 0),
    'off_spikes_first100_pct_mean': ((50,), 100, 0),
    't2_on_S_p05_p95_median_ms': ((31, 56, 41), 1000, 0),
    'rel_rate_ms_family_flies': ((0.45, 0.81), 1, 2),
    'rel_rate_1-8ol_family_flies': ((0.87, 0.68), 1, 2),
}


def check_planner(code: Dict[str, Sequence[float]],
                  checks: Dict[str, Tuple[Tuple[float, ...], float, int]] = PLANNER_CHECKS) -> List[dict]:
    """Returns one row per check: the planner's values, the code's (scaled), and whether they agree within rounding."""
    rows = []
    for name, (want, scale, nd) in checks.items():
        got = tuple(float(x) * scale for x in code[name])
        if len(got) != len(want):
            raise ValueError(f'{name}: {len(got)} values, the check has {len(want)}')
        ok = all(abs(g - w) <= 0.5 * 10 ** -nd + 1e-9 for g, w in zip(got, want))
        rows.append(dict(name=name, planner=list(want), code=[round(g, nd + 2) for g in got], decimals=nd,
                         agrees=bool(ok)))
    return rows


def _file_sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _range(x) -> Tuple[float, float]:
    x = list(x)
    return float(min(x)), float(max(x))


def reason_values(shared: Path) -> pd.DataFrame:
    """Returns the values the planner's reasons cite from other experiments, each read from its committed table.

    Columns 'reason', 'experiment', 'file', 'selector', 'column', 'value'. The family's calibrated pattern parts and the
    spread of its calibrated loss over the eight wirings come from E129's readings (the family's own rows), static APL
    from E95's readings, APL's time constant at the 2 s pulse from E82's runs (each arm minus the fit's 2 s point with
    the default time constant, part B1), the searched input's gain from E86's gate and the male calyces from E130's
    readings.
    """
    res = shared / 'results'
    rows = []

    def add(reason, exp, file, sel, col, value):
        rows.append(dict(reason=reason, experiment=exp, file=str(Path('results') / file), selector=sel, column=col,
                         value=float(value)))

    f129 = 'E129_block_survivor_static/readings.csv'
    r = pd.read_csv(res / f129)
    fam = r[(r['kind'] == 'family') & (r['readout'] == 'cal')].set_index('unit')
    for u, lab in (('all8', 'all eight wirings'), ('S4', 'S'), ('T4', 'T'), ('real', 'the hemibrain wiring')):
        add(f"the family's calibrated pattern part, {lab}", 'E129', f129, f'kind family, readout cal, unit {u}',
            'all_pattern', fam.loc[u, 'all_pattern'])
    one = fam.loc[[u for u in fam.index if u in WIRINGS or u == 'real']]
    add("the family's calibrated pattern part, lowest of the eight wirings", 'E129', f129,
        'kind family, readout cal, single wirings', 'all_pattern', one['all_pattern'].min())
    add("the family's calibrated pattern part, highest of the eight wirings", 'E129', f129,
        'kind family, readout cal, single wirings', 'all_pattern', one['all_pattern'].max())
    add("the family's calibrated loss, SD over the eight wirings (ddof 1; E72's draw term)", 'E129', f129,
        'kind family, readout cal, single wirings', 'all_loss', one['all_loss'].std(ddof=1))

    f95 = 'E95_static_apl/readings.csv'
    r = pd.read_csv(res / f95)
    for v, lab in (('C2_bcko', 'the whole APL static per KC and odor (C2)'), ('ref', 'the fits (reference)')):
        x = r[(r['variant'] == v) & (r['unit'] == 'S') & (r['readout'] == 'cal')]
        add(f'static APL, S: {lab}', 'E95', f95, f'variant {v}, unit S, readout cal', 'z8_pattern',
            x['z8_pattern'].iloc[0])

    f82 = 'E82_flies_protocol/runs.csv'
    r = pd.read_csv(res / f82).set_index('point_id')
    for model in ('fit0', 'C2'):
        base = r.loc[f'B1_{model}', 'cal_pattern']
        for p in r.index[(r['part'] == 'B5') & (r['model'] == model) & ~r['identity'].astype(bool)]:
            add(f"APL's time constant {r.loc[p, 'apl_taum']:g} s at the 2 s pulse, {model}, minus its B1 point",
                'E82', f82, f'{p} minus B1_{model}', 'cal_pattern', r.loc[p, 'cal_pattern'] - base)

    f86 = 'E86_reduced_input_candidate/gate.csv'
    r = pd.read_csv(res / f86)
    add('the reduced searched input on held-out wirings against unseen recordings, mean gain of the two splits',
        'E86', f86, 'change gate', 'g', r.loc[r['change'] == 'gate', 'g'].iloc[0])

    f130 = 'E130_malecns_family/readings.csv'
    r = pd.read_csv(res / f130)
    r = r[r['readout'] == 'cal'].set_index('unit')
    for u, lab in (('fit:malecns-left', 'the male CNS left calyx'), ('fit:malecns-right', 'the male CNS right calyx'),
                   ('family:real', 'the hemibrain wiring')):
        add(f'the family fitted on {lab}, calibrated pattern part', 'E130', f130, f'unit {u}, readout cal',
            'all_pattern', r.loc[u, 'all_pattern'])
    return pd.DataFrame(rows)


def reads(out_dir: Optional[Path] = None, res_dir: Optional[Path] = None) -> dict:
    """Computes the planner's reads from step 0's saved traces and spike times (no simulation) and writes them.

    Writes to `res_dir` (results/E121_onset_race/): reads_apl_levels.csv (per wiring and odor), reads_levels.csv (the
    level at completion, per wiring and pooled, overall and by group, breadth and odor), reads_survival.csv,
    reads_windows.csv (per wiring and condition), reads_relative_rates.csv, reads_realism_family.csv, reads_reasons.csv,
    reads_additions.csv (the two reads added after the placement, per wiring), reads_summary.json (the headline values
    and the planner's checks; the additions and their checks under 'additions') and provenance_reads.json. Returns the
    summary.
    """
    from kcq import harness as H
    out_dir = out_dir or (H.outputs_dir(EXPERIMENT) / 'step0')
    res_dir = res_dir or H.results_dir(EXPERIMENT)
    shared = H.shared_dir()
    step0_prov = json.loads((res_dir / 'provenance.json').read_text())
    hashes = step0_prov['spike_times_sha256']
    inputs: Dict[str, str] = {}
    lev_rows, comp_rows, surv_rows, win_rows, add_rows, checks = [], [], [], [], [], {}
    pooled_spikes: Dict[str, List[pd.DataFrame]] = {}
    n_kcs = None
    for lab in WIRINGS:
        exp, tag = WIRINGS[lab]
        arch = shared / 'data' / 'model_outputs' / exp / 'rebuild' / tag
        sp, cmp_ = {}, {}
        for c in CONDITIONS:
            p = out_dir / f'{lab}_{c}_spike_times.parquet'
            inputs[str(p)] = _file_sha(p)
            sp[c] = pd.read_parquet(p)
            if _sha(sp[c]) != hashes[f'{lab}_{c}_spike_times']:
                raise RuntimeError(f'{p}: not the spike times step 0 saved (hash differs)')
            ap_ = arch / ('expected_spike_counts_on.parquet' if c == 'on' else
                          'expected_spike_counts_acute_off.parquet')
            inputs[str(ap_)] = _file_sha(ap_)
            expect = pd.read_parquet(ap_).to_numpy()
            got = np.zeros_like(expect)
            np.add.at(got, (sp[c]['kc'].to_numpy(), sp[c]['odor'].to_numpy()), 1)
            if not np.array_equal(got, expect):
                raise RuntimeError(f'{lab} {c}: the saved spike times do not add up to the archive\'s spike counts')
            n_kcs = expect.shape[0] if n_kcs is None else n_kcs
            if expect.shape[0] != n_kcs:
                raise RuntimeError('the wirings differ in their number of KCs')
            checks[f'{lab}_{c}'] = dict(hash_matches_step0=True, counts_match_archive=True)
        p = out_dir / f'{lab}_on_apl_trace.parquet'
        inputs[str(p)] = _file_sha(p)
        apl = pd.read_parquet(p)
        apl.columns = apl.columns.astype(float)
        t = apl.columns.to_numpy(dtype=float)
        lv = pd.DataFrame({o: apl_levels(apl.loc[o].to_numpy(dtype=float), t) for o in apl.index}).T
        lv.index.name = 'odor_name'
        lev_rows.append(lv.reset_index().assign(wiring=lab, group=lambda d: d['odor_name'].map(odor_group)))
        if lev_rows[0]['odor_name'].tolist() != list(apl.index):
            raise RuntimeError(f'{lab}: odors in another order than the first wiring\'s')
        odors = list(apl.index)
        for c in CONDITIONS:
            cmp_[c] = second_spike_times(sp[c])
            w = window_reads(sp[c])
            win_rows.append(dict(wiring=lab, condition=c, **w))
        comp = levels_at_completion(completion_table(sp['on'], odors, n_kcs), apl, lv)
        comp['wiring'] = lab
        comp_rows.append(comp)
        surv_rows.append(dict(wiring=lab, **survival_reads(cmp_['on'], cmp_['off'])))

        # (the additions after placement: the archives' on-not-off counts, the family's and E110 arm A's, and the
        # KC's APL->KC weight and breadth value as tellers of survival)
        resp = {}
        for c, f in (('on', 'expected_responses_on.parquet'), ('off', 'expected_responses_acute_off.parquet')):
            inputs[str(arch / f)] = _file_sha(arch / f)
            resp[c] = pd.read_parquet(arch / f)
        fam = on_not_off_counts(resp['on'], resp['off'])
        a_dir = shared / 'data' / 'model_outputs' / E110_ARM_A[0] / E110_ARM_A[1]
        var = e110_arm_a_variant(lab)
        a_resp = {}
        for c, cond in (('on', 'on'), ('off', 'acute_off')):
            fr, fc = a_dir / f'{var}_{cond}_responses.parquet', a_dir / f'{var}_{cond}_spike_counts.parquet'
            inputs[str(fr)], inputs[str(fc)] = _file_sha(fr), _file_sha(fc)
            a_resp[c] = pd.read_parquet(fr)
            if not np.array_equal(a_resp[c].to_numpy() > 0, pd.read_parquet(fc).to_numpy() >= 2):
                raise RuntimeError(f'{fr}: responses are not two or more spikes of {fc.name}')
            if not a_resp[c].index.equals(resp[c].index):
                raise RuntimeError(f'{fr}: KCs differ from the rebuild archive\'s')
        kc0 = on_not_off_counts(a_resp['on'], a_resp['off'])
        kc_ids = resp['on'].index.get_level_values('kc_id')
        if kc_ids.has_duplicates:
            raise RuntimeError(f'{lab}: duplicate KC IDs in the rebuild archive')
        for f in ('wAPLKC.parquet', 'breadth.parquet'):
            inputs[str(arch / f)] = _file_sha(arch / f)
        wt = kc_apl_weight(pd.read_parquet(arch / 'wAPLKC.parquet'), kc_ids)
        br = pd.read_parquet(arch / 'breadth.parquet').iloc[:, 0].reindex(kc_ids).to_numpy(dtype=float)
        if np.isnan(br).any():
            raise RuntimeError(f'{lab}: KCs without a breadth value')
        add_rows.append(dict(wiring=lab, family_n_on_responders=fam['n_on_responders'],
                             family_n_on_not_off=fam['n_on_not_off'], family_share_on_not_off=fam['share_on_not_off'],
                             e110_variant=var, kckc0_n_on_responders=kc0['n_on_responders'],
                             kckc0_n_off_responders=kc0['n_off_responders'], kckc0_n_on_not_off=kc0['n_on_not_off'],
                             kckc0_share_on_not_off=kc0['share_on_not_off'],
                             **survival_by_weight(cmp_['on'], cmp_['off'], wt, br)))
        for c in CONDITIONS:
            done = pooled_spikes.setdefault(c, [])
            done.append(sp[c].assign(kc=sp[c]['kc'] + len(done) * n_kcs))
    L = pd.concat(lev_rows, ignore_index=True)
    C = pd.concat(comp_rows, ignore_index=True)
    SV = pd.DataFrame(surv_rows)
    AD = pd.DataFrame(add_rows)
    # (the family's counts from the archives' responses must be the spike times' own)
    if not (AD['family_n_on_not_off'].to_numpy() == SV['n_on_not_off'].to_numpy()).all():
        raise RuntimeError('the rebuild archives\' on-not-off counts differ from the saved spike times\'')
    AD.insert(AD.columns.get_loc('auc_apl_weight'), 'auc_count', SV['auc_count'].to_numpy())
    add_code = {
        'family_n_on_not_off_archive': tuple(AD['family_n_on_not_off']),
        'kckc0_n_on_not_off_range': _range(AD['kckc0_n_on_not_off']),
        'kckc0_n_on_responders_range': _range(AD['kckc0_n_on_responders']),
        'auc_apl_weight_range': _range(AD['auc_apl_weight']),
        'auc_breadth_range': _range(AD['auc_breadth']),
        'weight_survivors_over_mean_range': _range(AD['weight_survivors_over_mean']),
        'weight_removed_over_mean_range': _range(AD['weight_removed_over_mean']),
    }
    add_pc = check_planner(add_code, PLACEMENT_CHECKS)
    # (pooled over the wirings: each wiring's KCs counted apart)
    for c in CONDITIONS:
        win_rows.append(dict(wiring='pooled', condition=c,
                             **window_reads(pd.concat(pooled_spikes[c], ignore_index=True))))
    W = pd.DataFrame(win_rows)
    lev_sum = []
    for name, members in [(lab, [lab]) for lab in WIRINGS] + [('pooled', list(WIRINGS))]:
        s = level_summary(C[C['wiring'].isin(members)])
        s.insert(0, 'wiring', name)
        lev_sum.append(s)
    LS = pd.concat(lev_sum, ignore_index=True)

    # the family's per-odor response rates with APL on against the control flies' (relative to the mean odor's)
    cnt = C.groupby('odor_name').size().reindex(odors)
    fam_rate = cnt / (len(WIRINGS) * n_kcs)
    f_yang = shared / 'data' / 'yang_apl16_kc_olfdb' / 'computed' / 'response_fraction_recording_mean.csv'
    inputs[str(f_yang)] = _file_sha(f_yang)
    y = pd.read_csv(f_yang)
    y = y[y['genotype'] == 'TNT3c'].set_index('odor')['frac_responding_mean'].reindex(odors)
    if y.isna().any():
        raise RuntimeError(f'{f_yang}: odors missing')
    noise = pd.Series({o: H.data_target(f'noise_ctrl_frac_{o}') for o in odors})
    RR = pd.DataFrame({'family_rate_on': fam_rate, 'family_rel': fam_rate / fam_rate.mean(),
                       'flies_frac_1z': y, 'flies_rel_1z': y / y.mean(),
                       'flies_frac_noise': noise, 'flies_rel_noise': noise / noise.mean()})
    RR.index.name = 'odor_name'
    RR.insert(0, 'group', [odor_group(o) for o in RR.index])

    # the family's standard realism table (results/apl_effect/by_variant.csv)
    from kcq import apl_effect as AE
    f_bv = shared / 'results' / 'apl_effect' / 'by_variant.csv'
    inputs[str(f_bv)] = _file_sha(f_bv)
    bv = pd.read_csv(f_bv)
    rrows = []
    for lab, (exp, tag) in WIRINGS.items():
        x = bv[(bv['folder'] == exp) & (bv['tag'] == tag)]
        if len(x) != 1:
            raise RuntimeError(f'{lab}: {len(x)} rows in by_variant.csv for {exp}/{tag}')
        x = x.iloc[0]
        rrows.append(dict(wiring=lab, rate_on=x['rate_on'], rate_off=x['rate_off'], ratio=x['rate_off'] / x['rate_on'],
                          rate_rho_on=x['rate_rho_on'], rate_rho_off=x['rate_rho_off'], opk_on_tv=x['opk_on_tv'],
                          opk_off_tv=x['opk_off_tv'], realism=x['realism'],
                          raw_loss=x['spike_counts_cosine__apl_effect_rmse'],
                          raw_pattern=x['spike_counts_cosine__apl_effect_rmse_centered'],
                          raw_mean_effect=x['spike_counts_cosine__apl_effect_mean']))
    RF = pd.DataFrame(rrows)
    mean = RF.drop(columns=['wiring', 'realism']).mean()
    RF = pd.concat([RF, pd.DataFrame([dict(wiring='mean', realism='', **mean.to_dict())])], ignore_index=True)
    kc_data = dict(rate_on=H.data_target('noise_ctrl_frac'), rate_off=H.data_target('noise_tnt_frac'),
                   ratio=H.data_target('noise_rate_ratio_fraction'),
                   rate_rho_ceiling_on=H.data_target('ctrl_rate_rho_ceiling'),
                   rate_rho_ceiling_off=H.data_target('tnt_rate_rho_ceiling'),
                   opk_on_band=H.data_target('noise_ctrl_odors_per_kc_tv_max'),
                   opk_off_band=H.data_target('noise_tnt_odors_per_kc_tv_max'),
                   rate_on_band=list(AE.RATE_ON_BAND), ratio_band=list(AE.RATE_RATIO_BAND),
                   rate_rho_on_min=AE.RATE_RHO_MIN, apl_effect_mean=H.data_target('apl_effect_mean'),
                   apl_effect_rmse_floor=H.data_target('apl_effect_rmse_floor'),
                   apl_effect_pattern_floor=H.data_target('apl_effect_pattern_floor'))

    RS = reason_values(shared)

    # the headline values, in the planner's terms
    lw = L.groupby('wiring')[['rest_over_peak', 'peak_over_sustained']].mean()
    pooled = LS[LS['wiring'] == 'pooled'].set_index(['by', 'level'])
    perw = LS[(LS['wiring'] != 'pooled') & (LS['by'] == 'overall')].set_index('wiring')
    on = W[(W['condition'] == 'on') & (W['wiring'] != 'pooled')].set_index('wiring')
    off = W[(W['condition'] == 'off') & (W['wiring'] != 'pooled')].set_index('wiring')
    po = pooled.loc[('overall', 'all')]
    code = {
        'n_on_responses_pooled': (len(C),),
        'rest_over_peak_range_wirings': _range(lw['rest_over_peak']),
        'peak_over_sustained_range_wirings': _range(lw['peak_over_sustained']),
        't_peak_range_ms': _range(L['t_peak']),
        'frac_peak_median': (po['frac_peak_median'],),
        'frac_peak_iqr': (po['frac_peak_q25'], po['frac_peak_q75']),
        'frac_peak_median_range_wirings': _range(perw['frac_peak_median']),
        'frac_sustained_median': (po['frac_sustained_median'],),
        'frac_sustained_iqr': (po['frac_sustained_q25'], po['frac_sustained_q75']),
        'below_sustained_pct': (po['share_below_sustained'],),
        'before_peak_pct': (po['share_before_peak'],),
        'before_peak_pct_range_wirings': _range(perw['share_before_peak']),
        'frac_peak_ester_alcohol_other': tuple(pooled.loc[('group', g), 'frac_peak_median']
                                              for g in ('ester', 'alcohol', 'other')),
        'frac_peak_breadth_1_8plus': tuple(pooled.loc[('breadth', b), 'frac_peak_median'] for b in ('1', '8+')),
        'frac_peak_odor_pa_bmyr_ms': tuple(pooled.loc[('odor', o), 'frac_peak_median'] for o in ('pa', 'B-myr', 'ms')),
        'n_off_responders_range': _range(SV['n_off_responders']),
        'n_survivors_range': _range(SV['n_survivors']),
        'auc_t2_range': _range(SV['auc_t2']),
        'auc_count_range': _range(SV['auc_count']),
        'median_t2_survivors_ms_range': _range(SV['median_t2_survivors']),
        'median_t2_removed_ms_range': _range(SV['median_t2_removed']),
        'early_removed_pct_range': _range(SV['share_early_removed']),
        'early_survivors_pct_range': _range(SV['share_early_survivors']),
        'n_on_not_off_range': _range(SV['n_on_not_off']),
        'on_not_off_pct_range': _range(SV['share_on_not_off']),
        'on_n_responders_mean': (on['n_responders'].mean(),),
        'on_n_first100_mean': (on['n_first100'].mean(),),
        'on_n_last250_mean': (on['n_last250'].mean(),),
        'on_n_last250_range': _range(on['n_last250']),
        'on_n_after_mean': (on['n_after'].mean(),),
        'on_spikes_first100_pct_mean': (on['share_spikes_first100'].mean(),),
        'on_spikes_first100_pct_range': _range(on['share_spikes_first100']),
        'on_exactly_two_pct_mean': (on['share_exactly_two'].mean(),),
        'on_exactly_two_pct_range': _range(on['share_exactly_two']),
        'on_median_spikes': (on['median_spikes'].median(),),
        'off_n_responders_mean': (off['n_responders'].mean(),),
        'off_n_last250_mean': (off['n_last250'].mean(),),
        'off_spikes_first100_pct_mean': (off['share_spikes_first100'].mean(),),
        't2_on_S_p05_p95_median_ms': (on.loc['S', 't2_p05'], on.loc['S', 't2_p95'], on.loc['S', 't2_median']),
        'rel_rate_ms_family_flies': (RR.loc['ms', 'family_rel'], RR.loc['ms', 'flies_rel_1z']),
        'rel_rate_1-8ol_family_flies': (RR.loc['1-8ol', 'family_rel'], RR.loc['1-8ol', 'flies_rel_1z']),
    }
    pc = check_planner(code)
    now = datetime.now(timezone.utc).isoformat(timespec='seconds')
    summary = dict(
        experiment=EXPERIMENT, part='reads', time=now, n_kcs=n_kcs, wirings=list(WIRINGS),
        code_values={k: [float(x) for x in v] for k, v in code.items()},
        planner_checks=pc, planner_checks_disagree=[r['name'] for r in pc if not r['agrees']],
        kc_data=kc_data, input_checks=checks,
        additions=dict(
            added='2026-11-23, after the placement (the placement planner claude-fable-5-2, spawned from viki '
                  '2026-11-23 22:23, read these inline from the model-output archive)',
            table='reads_additions.csv', code_values={k: [float(x) for x in v] for k, v in add_code.items()},
            placement_checks=add_pc, placement_checks_disagree=[r['name'] for r in add_pc if not r['agrees']],
            definitions=dict(
                on_not_off=('KC x odor entries responding with APL on and not with APL off (responses: two or more '
                            'spikes), from the rebuild archives\' expected responses (the family) and from E110\'s '
                            'arm A, KC-KC at 0 at fixed scales with both response counts held '
                            '(data/model_outputs/E110_no_kckc/A/, A_real for S)'),
                apl_weight=('among the APL-off responders, survivors respond with APL on too; each entry carries its '
                            'KC\'s summed APL->KC weight (wAPLKC.parquet of the rebuild archive, summed over the KC\'s '
                            'claws, unscaled) and its breadth value (breadth.parquet); AUC with lower favouring '
                            'survival; the survivors\' and the removed\'s median weight over the mean of all KCs\''))),
        definitions=dict(
            levels=("APL on's trace (inh_sims), each level including APL's activity before the odor: rest the trace at "
                    "the last step at or before 0 s, peak its maximum over 0 < t <= 0.5 s and its time, sustained its "
                    f"mean over the last {LATE_WINDOW} s of the pulse; rest / peak and peak / sustained as the mean "
                    "over the 16 odors per wiring"),
            completion='the second spike (two or more spikes: a response); the level at completion is the trace at '
                       'the last step at or before it, as a fraction of the odor\'s peak and sustained level',
            breadth='the number of the 16 odors the KC responds to with APL on',
            survival=('among the KC x odor entries responding with APL off, survivors respond with APL on too; AUC of '
                      'the APL-off second spike time (earlier favours survival) and spike count (more favours it); '
                      f'early: second spike before {EARLY_T2 * 1000:g} ms'),
            windows=f'half-open windows {[(n, lo, hi) for n, lo, hi in READ_WINDOWS]} s from onset; a window '
                    'responder has two or more spikes in it',
            relative_rates=("APL on, the family's responders per odor pooled over the eight wirings over the mean "
                            "odor's; the control flies' fraction responding (trial mean above 1 z, "
                            "response_fraction_recording_mean.csv, TNT3c, 16 odors) and in baseline-noise units "
                            "(registry noise_ctrl_frac_<odor>), each over the mean of the 16 odors")),
    )
    L.to_csv(res_dir / 'reads_apl_levels.csv', index=False)
    LS.to_csv(res_dir / 'reads_levels.csv', index=False)
    SV.to_csv(res_dir / 'reads_survival.csv', index=False)
    W.to_csv(res_dir / 'reads_windows.csv', index=False)
    RR.to_csv(res_dir / 'reads_relative_rates.csv')
    RF.to_csv(res_dir / 'reads_realism_family.csv', index=False)
    RS.to_csv(res_dir / 'reads_reasons.csv', index=False)
    AD.to_csv(res_dir / 'reads_additions.csv', index=False)
    (res_dir / 'reads_summary.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
    prov = dict(script='experiments/E121_onset_race.py', part='reads', time=now, code_commits=H.code_commits(),
                inputs_sha256=inputs, step0_provenance=str(res_dir / 'provenance.json'),
                reasons_sources=sorted(set(RS['file'])), registry=str(H.data_targets_path()),
                how='no simulation: reads of step 0\'s saved spike times and APL traces (checked against step 0\'s '
                    'hashes and the archives\' spike counts), results/apl_effect/by_variant.csv, Yang\'s response '
                    'fractions, the data-target registry and the other experiments\' tables in reads_reasons.csv; '
                    'since the placement (2026-11-23) also the rebuild archives\' responses, wAPLKC and breadth and '
                    'E110\'s arm-A responses and spike counts (reads_additions.csv)')
    (res_dir / 'provenance_reads.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('part', choices=['step0', 'reads'])
    ap.add_argument('--wirings', nargs='+', choices=list(WIRINGS), default=None)
    args = ap.parse_args()
    if args.part == 'reads':
        s = reads()
        print(json.dumps(dict(planner_checks_disagree=s['planner_checks_disagree'],
                              checks=[r for r in s['planner_checks'] if not r['agrees']],
                              placement_checks=s['additions']['placement_checks']), indent=1))
        return
    gate = step0(wirings=args.wirings)
    print(json.dumps({k: gate[k] for k in ('verdict', 'share_on_pooled', 'share_on_pooled_plateau_max',
                                            'apl_frac_peak_at_completion_on', 'apl_peak_over_late',
                                            'share_on_range_over_wirings',
                                            'share_on_by_group', 'share_on_by_breadth', 'share_off_pooled',
                                            'apl_t90_s')}, indent=1))


if __name__ == '__main__':
    main()
