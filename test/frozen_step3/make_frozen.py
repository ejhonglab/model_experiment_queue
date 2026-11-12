#!/usr/bin/env python3
"""Makes the frozen files of the step-3 move (the queue's generally useful functions moved to hong2p and pyolfsysm,
2026-11-12): small inputs cut from real data (C2's saved spike counts and responses, E82's spike times, Yang's per-trial
KC data and per-recording distance matrices, the report's matrix, the model's ORN input) and synthetic traces, and every
case's output computed by the queue's code (cases.py). Run once, on the code before the move, from the repository root:
`./venv/bin/python test/frozen_step3/make_frozen.py`. Writes frozen_hong2p.npz, frozen_pyolfsysm.npz and frozen_queue.npz
here; the first two are copied to hong2p's test/data/ and olfsysm's test/test_data/.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from frozen_step3 import frozen_outputs as F  # noqa: E402
from frozen_step3 import cases as C  # noqa: E402
from kcq import harness as H  # noqa: E402

ARCHIVE = H.shared_dir() / 'data' / 'model_outputs'
C2 = ARCHIVE / 'E67_apl_effect_fits' / 'rate0.068_aplboth_k2_cap3'
SPIKE_TIMES = ARCHIVE / 'E82_flies_protocol' / 'rate0.068_aplboth_k2_pulse2_cap3_on_spike_times.parquet'


def _frame(d: dict, name: str, df: pd.DataFrame) -> None:
    d[f'{name}/values'] = df.to_numpy(dtype=float)
    d[f'{name}/index'] = np.asarray(list(df.index))
    d[f'{name}/columns'] = np.array([str(c) for c in df.columns])


def _pair_series(d: dict, name: str, s: pd.Series) -> None:
    d[f'{name}/values'] = s.to_numpy(dtype=float)
    d[f'{name}/a'] = np.array([a for a, _ in s.index])
    d[f'{name}/b'] = np.array([b for _, b in s.index])
    d[f'{name}/name'] = np.array(['' if s.name is None else str(s.name)])


def build_inputs() -> dict:
    d = {}
    for cond, f in (('on', 'on'), ('off', 'acute_off')):
        for kind, short in (('spike_counts', 'counts'), ('responses', 'resp')):
            df = pd.read_parquet(f'{C2}_{f}_{kind}.parquet').iloc[::4]
            df.index = df.index.get_level_values('kc_id')
            _frame(d, f'{short}_{cond}', df)
    sq = pd.read_csv(H.apl_pair_dir() / 'TNT3c_dist_apl16.csv', index_col=0)
    _frame(d, 'sq_ctrl', sq)
    d['sq_ctrl/index'] = np.array([str(i) for i in sq.index])
    orn = H.abbrev_columns(H.apl16_orn_deltas())
    _frame(d, 'orn', orn)
    d['orn/index'] = np.array([str(i) for i in orn.index])
    mats = H.yang_recording_matrices('TNT3c')
    tm = mats[sorted(mats)[0]]
    _frame(d, 'tm_mat', tm.iloc[:300])
    d['ref_amp'] = np.sort(tm.to_numpy()[tm.to_numpy() > 1.0])
    for name, g in (('rec_ctrl', 'TNT3c'), ('rec_tnt', 'TNT')):
        r = H.yang_recording_rdms(g)
        d[f'{name}/values'] = r.to_numpy(dtype=float)
        for lvl in ('date', 'fly_num', 'hemibrain'):
            v = r.index.get_level_values(lvl)
            d[f'{name}/{lvl}'] = np.array([int(x) for x in v]) if lvl == 'fly_num' else np.array([str(x) for x in v])
        d[f'{name}/a'] = np.array([a for a, _ in r.columns])
        d[f'{name}/b'] = np.array([b for _, b in r.columns])
    _pair_series(d, 'data_effect', H.yang_apl_effect())
    tr = H.yang_recording_trials('TNT3c')
    keys = sorted(tr)
    d['trials_ctrl'] = tr[keys[0]][:, :200]
    d['trials_ctrl2'] = tr[keys[1]][:, :200]
    tt = H.yang_recording_trials('TNT')
    d['trials_tnt'] = tt[sorted(tt)[0]][:, :200]
    d['odors'] = np.array(H.APL16_ODORS)
    st = pd.read_parquet(SPIKE_TIMES)
    st = st[st['kc'] < 300]
    d['spikes/kc'] = st['kc'].to_numpy(dtype=np.int64)
    d['spikes/odor'] = st['odor'].to_numpy(dtype=np.int64)
    d['spikes/time_s'] = st['time_s'].to_numpy(dtype=float)
    d['spikes/n_kc'] = np.array([300], dtype=np.int64)
    d['cal_target'] = np.array([H.data_target('ctrl_participation_ratio_resp')])
    d['volume_rate'] = np.array([C.RATE])
    # synthetic traces (as test_e84_trace_windows's, smaller) and a noisy mean trace for the kernel fits
    rng = np.random.default_rng(0)
    n_units, n_trials, gap, n_before, n_after = 10, 4, 139, 15, 119
    T = n_before + n_trials * gap + n_after
    onsets = n_before + 5 + gap * np.arange(n_trials)
    t = np.arange(T)
    x = 100 + 20 * rng.random((n_units, 1)) + 0.002 * t[None] * rng.standard_normal((n_units, 1))
    x = x + rng.standard_normal((n_units, T))
    for n in onsets:
        amp = rng.gamma(1.0, 2.0, (n_units, 1)) * (rng.random((n_units, 1)) < 0.5)
        x = x + amp * H.indicator_kernel((t - n) / C.RATE, 1.0, 10.0)[None]
    d['traces'], d['onsets'] = x, onsets.astype(np.int64)
    ts = np.arange(-15, 119) / C.RATE
    d['synth_t'] = ts
    d['synth_y'] = 2.0 * H.kernel_signal(ts, 0.9, 7.0, 'steady', pulse=2.0) + 0.05 * np.random.default_rng(4).standard_normal(len(ts))
    return d


# the inputs each group's cases use (the queue's file has none: its cases read the queue's data)
INPUTS = {'hong2p': {'sq_ctrl', 'counts_on', 'counts_off', 'orn', 'tm_mat', 'rec_ctrl', 'rec_tnt', 'data_effect',
    'trials_ctrl', 'trials_tnt', 'odors', 'traces', 'onsets', 'synth_t', 'synth_y', 'volume_rate'},
    'pyolfsysm': {'counts_on', 'counts_off', 'resp_on', 'resp_off', 'tm_mat', 'trials_ctrl', 'trials_ctrl2',
    'trials_tnt', 'odors', 'spikes', 'ref_amp', 'cal_target', 'data_effect'},
    'queue': set()}


def expected(cases: C.Cases, O: dict) -> dict:
    out = {}
    for name, (fn, _) in cases.items():
        print(name, flush=True)
        F.flatten(fn(O), name, out)
    return out


def main() -> None:
    inp = build_inputs()
    for name, cases in (('hong2p', C.hong2p_cases()), ('pyolfsysm', C.pyolfsysm_cases()), ('queue', C.queue_cases())):
        keep = INPUTS.get(name)
        sub = {k: v for k, v in inp.items() if keep is None or k.split('/')[0] in keep}
        F.save(HERE / f'frozen_{name}.npz', sub, expected(cases, F.objects(sub)))
        print(f'wrote frozen_{name}.npz', flush=True)


if __name__ == '__main__':
    main()
