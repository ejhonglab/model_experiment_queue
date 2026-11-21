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
experiment closes without code.

Groups (E77_input_screen's): ester group 2h, IaA, pa, eb, ep, ms; alcohol group 1-8ol, 1-5ol, 1-6ol, 6al; the rest
(2-but, t2h, paa, -aPine, euc, B-myr) "other". Breadth: the number of the 16 odors the KC responds to in that
condition (bins 1, 2-3, 4-7, 8+).

Outputs (results/E121_onset_race/): step0_shares.csv (share complete before APL's 90 % point, overall, by group,
by breadth, by odor, per wiring and pooled, APL on and off, both plateau readings), step0_apl_rise.csv (per wiring
and odor: plateau, time to 50 % and 90 %, both readings), step0_gate.json, provenance.json. Per-wiring spike times
and APL traces under outputs/E121_onset_race/step0/.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E121_onset_race.py step0
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
                             median_t_apl90=float(g['t_apl90'].median()) if n else float('nan')))
    return pd.DataFrame(rows)


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
            comp = completion_table(o[f'spikes_{c}'], o['odors'], o['n_kcs'])
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
    race = pooled >= GATE_SHARE
    rl = rise_all[rise_all.plateau == 'late']
    gate = dict(
        experiment=EXPERIMENT, step='step0',
        rule=(f'fewer than {GATE_SHARE:.0%} of APL-on responses complete (second spike) before APL reaches '
              f'{RISE_FRACTION:.0%} of its plateau: no race, the experiment closes without code'),
        verdict=('race: step 0 passes (the share is at or above the gate); steps 1-3 wait on the olfsysm option'
                 if race else 'no race: the experiment closes without code'),
        race=bool(race),
        share_on_pooled=pooled, share_on_pooled_plateau_max=pooled_max,
        verdict_same_with_plateau_max=bool((pooled_max >= GATE_SHARE) == race),
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('part', choices=['step0'])
    ap.add_argument('--wirings', nargs='+', choices=list(WIRINGS), default=None)
    args = ap.parse_args()
    gate = step0(wirings=args.wirings)
    print(json.dumps({k: gate[k] for k in ('verdict', 'share_on_pooled', 'share_on_range_over_wirings',
                                            'share_on_by_group', 'share_on_by_breadth', 'share_off_pooled',
                                            'apl_t90_s')}, indent=1))


if __name__ == '__main__':
    main()
