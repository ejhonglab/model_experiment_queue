#!/usr/bin/env python3
"""E90 under the settled pre-odor level: every reading of the model's stages against the flies' stages (E90,
experiments/E90_model_stages.py) whose stage starts from zero at the run's start, or inherits that start-up, read again
with the pre-odor mean taken over SETTLED_PRE_WINDOW (-0.25 to 0 s) instead of E90's PRE_WINDOW (-0.5 to 0 s).

olfsysm's KC layer, and with it APL, starts from zero at time_start (-0.5 s) and reaches its spontaneous level about 0.2
to 0.3 s later (unresolved_issues.md, 2026-11-15 17:48), so E90's pre-odor means of the claws and of APL read low, and the
boutons inherit APL's start-up with APL on. The PNs start at time_pre_start (-2 s) and drift by at most 0.018 Hz between
the two windows, so the PN stage, and the boutons and claws with APL off, move only by that drift. The somata (spike counts
in the odor window) are not affected. The settled window is E98's SETTLED_PRE_WINDOW (closed on the left, open on the
right), the level the queue uses since E98's amendment of 2026-11-15 17:59.

E90's saved stage runs hold window and pre-odor means only, not traces, so the settled levels come from short runs at the
saved parameters (no fitting) of every run E90 read: part B's eight wirings (C2's family; E83's `build_simulator`), B5's
point at the flies' 2 s pulse and D1's four firing-ceiling fits (E82's `build_simulator`), APL on and acutely off, each on
LEVEL_ODORS only with the breadth held at the full input's. Before the odor every odor's run is the same, so a few odors
give every unit's pre-odor level. Each short run is checked before its levels are used: its spike counts and responses for
those odors equal the saved run's, its pre-odor means are the same for every odor it ran, and its means over E90's window
equal E90's saved pre-odor means (within LEVEL_TOL) for the PNs, the boutons and the claws.

The settled levels then replace E90's pre-odor means in the saved stage matrices, and E90's own reading, scoring and
verdict functions run on them unchanged, with E90's seeds (`read_wiring`, `pn_rows`, `soma_rows`, `noise_check`,
`model_steps_of`, `verdict_rows`, `group_verdict`). The same pipeline on the saved means ('e90') must reproduce E90's saved
readings exactly (`reproduces`), so every difference is the pre-odor level's. viki's reconstruction of the same levels
(results/preodor_startup_check_viki/e90_i0_settled.py and claw_steps_fits.py in the shared directory: APL's saved traces of
E98 B1 / E101 A3, the PN trace of check_layers_real.py and the family's weights) is compared where its files are present.

Parts (CLI `--part`):
- levels (simulations): the short runs; preodor_levels.parquet (every unit's level over both windows) and
  preodor_levels.csv (per run and stage: the means, E90's window against E90's saved means, the identity checks).
- score (no simulations): the readings of part B, B5 and D1 at both levels, E90's tables at the settled level
  (stages_model.csv, noise_check.csv, steps_groups.csv, verdict.json, verdict_rows.csv, steps_by_reading.csv,
  stages_b5.csv, b5_steps.csv, d1_steps.csv) and the comparisons with E90's (compare_stages.csv, compare_steps.csv,
  compare_verdict.csv, silent_claws.csv).
- table: comparison_table.md.

Outputs: results/E90_model_stages/settled_preodor/ (E90's own outputs are not touched); the readings' pair tables under
outputs/E90_model_stages/settled_preodor/.

Usage (from the repository root):
    env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E90_settled_preodor.py --part levels
    ./venv/bin/python -u experiments/E90_settled_preodor.py --part score
    ./venv/bin/python -u experiments/E90_settled_preodor.py --part table
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E90_model_stages as E90
from kcq import harness as H

SETTLED_PRE_WINDOW = (-0.25, 0.0)       # E98's SETTLED_PRE_WINDOW: the pre-odor level after the KC layer's start-up (s)
LEVEL_ODORS = ('2h', 'eb', '1-8ol', 'paa')   # the short runs' odors (any odors give the pre-odor levels; four = 4 threads)
LEVEL_TOL = 1e-9                        # E90's window in a short run against E90's saved pre-odor means
STAGES = ('pn', 'bouton', 'claw')
VERSIONS = ('e90', 'settled')
VIKI_DIR = ('results', 'preodor_startup_check_viki', 'e90i0')
B5_TAG = 'real_pulse2'


def out_dir() -> Path:
    d = E90.out_dir() / 'settled_preodor'
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir() -> Path:
    return E90.local_dir('settled_preodor')


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=str) + '\n')


# ---------------------------------------------------------------------------
# the windows and the levels

def pre_window_steps(t: np.ndarray, window: Tuple[float, float], closed_left: bool) -> np.ndarray:
    """The indices of the steps at times `t` inside `window` (lo, hi): open on the right (the odor starts at hi), closed
    on the left if `closed_left` (E98's settled window), open otherwise (E90's PRE_WINDOW: olfsysm's step at time_start
    has no bouton value)."""
    t = np.asarray(t, dtype=float)
    lo, hi = window
    left = (t >= lo - 1e-9) if closed_left else (t > lo + 1e-9)
    return np.flatnonzero(left & (t < hi - 1e-9))


def window_steps(t: np.ndarray) -> Dict[str, np.ndarray]:
    """The pre-odor steps of both levels: 'e90' (E90's own `_windows`, PRE_WINDOW) and 'settled' (SETTLED_PRE_WINDOW)."""
    return dict(e90=E90._windows(t, E90.MODEL_WINDOW)['pre'], settled=pre_window_steps(t, SETTLED_PRE_WINDOW, True))


def settled_level(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    """The settled pre-odor level of traces `x` (..., steps) at step times `t`: their mean over SETTLED_PRE_WINDOW."""
    return np.asarray(x, dtype=float)[..., pre_window_steps(t, SETTLED_PRE_WINDOW, True)].mean(axis=-1)


def read_levels(mp, rv, W: pd.DataFrame, order: np.ndarray, odors: Sequence[str]) -> dict:
    """The pre-odor levels of the last simulation (olfsysm's variables `mp`, `rv`; output column j is olfsysm's odor
    order[j]) over both windows (`window_steps`), as E90's `read_stages` reads its stages: the PNs (pn_sims, glomeruli in
    the wiring's order), the boutons (bouton_sims) and the claws (E90's `claw_inputs` from the boutons and APL's activity
    one step earlier), and APL (inh_sims). Returns {version: {'pn', 'bouton', 'claw': unit arrays, 'apl': float}} from the
    first odor, and 'spread': the largest difference of any odor's levels from the first's (0 when the pre-odor period
    is the same for every odor)."""
    from pyolfsysm.mb_model import get_sim_var
    t = E90.step_times(mp, rv)
    ws = window_steps(t)
    pn = get_sim_var(rv.pn, 'pn_sims')
    btn = get_sim_var(rv.pn, 'bouton_sims')
    inh = get_sim_var(rv.kc, 'inh_sims')
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    Wm = W.to_numpy(dtype=float)
    if btn.shape[1] != Wm.shape[1] or len(wk) != Wm.shape[0]:
        raise RuntimeError(f'bouton_sims {btn.shape}, {len(wk)} APL->claw weights, wiring {Wm.shape}')
    per_odor = []
    for j in range(len(odors)):
        i = int(order[j])
        lv = {}
        for v, idx in ws.items():
            b = np.asarray(btn[i])[:, idx]
            if np.isnan(b).any():
                raise RuntimeError(f'NaN bouton rates in the {v} window')
            a = np.asarray(inh[i]).ravel()
            lv[v] = dict(pn=np.asarray(pn[i])[:, idx].mean(axis=1), bouton=b.mean(axis=1),
                         claw=E90.claw_inputs(Wm, b, a[idx - 1], wk).mean(axis=1), apl=float(a[idx].mean()))
        per_odor.append(lv)
    first = per_odor[0]
    spread = 0.0
    for lv in per_odor[1:]:
        for v in VERSIONS:
            for k in STAGES + ('apl',):
                spread = max(spread, float(np.max(np.abs(np.asarray(lv[v][k]) - np.asarray(first[v][k])))))
    return dict(first, spread=spread, n_steps={v: int(len(ix)) for v, ix in ws.items()},
                window_s={v: [float(t[ix[0]]), float(t[ix[-1]])] for v, ix in ws.items()})


def _identity(r: pd.DataFrame, s: pd.DataFrame, saved: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    """A short run's responses and counts `r`, `s` against the saved run's (all odors) for the odors it ran."""
    rs, ss = (H.abbrev_columns(x) for x in saved)
    ra, sa = H.abbrev_columns(r), H.abbrev_columns(s)
    cols = list(sa.columns)
    nd = int((sa.to_numpy(dtype=float) != ss[cols].to_numpy(dtype=float)).sum())
    nr = int((ra.to_numpy(dtype=float) != rs[cols].to_numpy(dtype=float)).sum())
    same_index = bool(sa.index.equals(ss.index))
    return dict(odors=cols, counts_differ=nd, responses_differ=nr, index_equal=same_index,
                exact=bool(nd == 0 and nr == 0 and same_index))


def _short_runs(sim, info: dict, W: pd.DataFrame, params: Dict[str, Dict[str, float]],
                saved: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]) -> Dict[str, dict]:
    """Both conditions of a short-run simulator (`params` and `saved` keyed 'on' / 'off'): identity and levels."""
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    mp.kc.save_inh_sims = True
    wp = np.asarray(rv.kc.wPNKC, dtype=float)
    if not (wp.shape == W.shape and np.array_equal(wp, W.to_numpy(dtype=float))):
        raise RuntimeError("olfsysm's wPNKC differs from the wiring's claw x bouton matrix")
    out = {}
    for c, _ in E90.CONDS:
        r, s = sim(params[c])
        ident = _identity(r, s, saved[c])
        if not ident['exact']:
            raise RuntimeError(f'the short run differs from the saved run ({c}): {ident}')
        out[c] = dict(identity=ident, levels=read_levels(mp, rv, W, order, [H.odor_abbrev(x) for x in r.columns]))
    return out


def _odor_columns(orn: pd.DataFrame, odors: Sequence[str] = LEVEL_ODORS) -> List:
    cols = [c for c in orn.columns if H.odor_abbrev(c) in odors]
    if len(cols) != len(odors):
        raise RuntimeError(f'{len(cols)} of the odors {odors} in the ORN input')
    return cols


def c2_short_runs(label: str, odors: Sequence[str] = LEVEL_ODORS) -> Dict[str, dict]:
    """Part B's wiring `label` (C2's family at its own fit; E83's `build_simulator` on `odors`, the breadth held)."""
    from experiments import E83_cross_wiring_fit as E83
    P = E83.wiring_parts(label)
    sim, info = E83.build_simulator(label, orn=P['orn'][_odor_columns(P['orn'], odors)], breadth=P['breadth'])
    params = dict(on=dict(P['p']), off=dict(P['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
    res = _short_runs(sim, info, P['W'], params, P['saved'])
    del sim, info
    return res


def b5_short_runs(odors: Sequence[str] = LEVEL_ODORS) -> Dict[str, dict]:
    """B5's point (E82's re-held B1_C2, the real wiring at the 2 s pulse) as E90's `part_b5` builds it, on `odors`."""
    from experiments import E80_pulse as E80P
    from experiments import E82_protocol_runs as R
    from experiments import E77_input_screen as S
    b = E90.b5_spec()
    spec = b['spec']
    m = R.model_parts('C2')
    breadth = R.point_breadth(spec, m)
    if E80P._sha(breadth.to_numpy()) != b['rec']['breadth_sha256']:
        raise RuntimeError("B5: the breadth pre-pass at 2 s differs from E82's B1_C2 point's")
    orn = S.c2_parts()['orn']
    sim, info = R.build_simulator(R.model_kws(spec, m), wPNKC=None, breadth=breadth, params=b['params'],
                                  kc_kc_tau=spec['kc_kc_tau'], orn=orn[_odor_columns(orn, odors)])
    res = _short_runs(sim, info, info['wPNKC'], dict(on=b['params'], off=b['params_off']), E90.b5_saved())
    del sim, info
    return res


def d1_short_runs(model: str, odors: Sequence[str] = LEVEL_ODORS) -> Dict[str, dict]:
    """D1's firing-ceiling fit `model` as E90's `d1_run` builds it, on `odors`."""
    from experiments import E82_protocol_runs as R
    from experiments import E77_input_screen as S
    m = R.model_parts(model)
    W = None if m['wiring'] is None else R.wiring_parts(m['wiring'])['W']
    orn = S.c2_parts()['orn']
    sim, info = R.build_simulator(m['kc'], wPNKC=W, breadth=R.own_breadth(m), params=m['p'],
                                  orn=orn[_odor_columns(orn, odors)])
    params = dict(on=dict(m['p']), off=dict(m['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
    res = _short_runs(sim, info, info['wPNKC'], params, m['saved'])
    del sim, info
    return res


def runs() -> List[Tuple[str, str, str]]:
    """Every run of E90's whose stages were read: (part, tag, wiring of its units)."""
    out = [('B', w, w) for w in E90.ALL_WIRINGS] + [('B5', B5_TAG, 'real')]
    d1 = json.loads((E90.out_dir() / 'd1.json').read_text())
    out += [('D1', m, d1[m]['wiring']) for m in E90.D1_MODELS]
    return out


def _saved_stage(part: str, tag: str, cc: str) -> dict:
    st = E90.load_stage(part, f'{tag}_{cc}')
    if st is None:
        raise RuntimeError(f"E90's saved stage run {part} {tag}_{cc} not found")
    return st


def part_levels(only: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The short runs of every run (`runs`; resumable: a run whose levels are saved under the local folder is skipped),
    each checked against E90's saved pre-odor means over E90's window; preodor_levels.parquet / .csv."""
    import gc
    t0 = time.time()
    d = local_dir() / 'levels'
    d.mkdir(exist_ok=True)
    for part, tag, _ in runs():
        if only and tag not in only:
            continue
        path = d / f'{part}_{tag}.json'
        if path.exists():
            continue
        t1 = time.time()
        res = (c2_short_runs(tag) if part == 'B' else b5_short_runs() if part == 'B5' else d1_short_runs(tag))
        rec = {}
        for c, cc in E90.CONDS:
            st = _saved_stage(part, tag, cc)
            lv = res[c]['levels']
            chk = {}
            for k in STAGES:
                x = st[f'{k}_pre'].to_numpy(dtype=float)
                if x.shape[0] != len(lv['e90'][k]):
                    raise RuntimeError(f'{part} {tag} {c} {k}: {x.shape[0]} saved units, {len(lv["e90"][k])} in the run')
                chk[k] = float(np.max(np.abs(x - np.asarray(lv['e90'][k])[:, None])))
            chk['apl'] = float(np.max(np.abs(st['apl']['pre'].to_numpy(dtype=float) - lv['e90']['apl'])))
            if max(chk.values()) > LEVEL_TOL:
                raise RuntimeError(f"{part} {tag} {c}: E90's window in the short run differs from E90's saved pre-odor "
                                   f'means: {chk}')
            rec[c] = dict(identity=res[c]['identity'], e90_window_vs_saved=chk, spread_over_odors=lv['spread'],
                          n_steps=lv['n_steps'], window_s=lv['window_s'],
                          levels={v: {k: np.asarray(lv[v][k]).tolist() if k != 'apl' else lv[v][k]
                                      for k in STAGES + ('apl',)} for v in VERSIONS})
        path.write_text(json.dumps(rec) + '\n')
        print(f'{part} {tag}: exact; E90 window vs saved {[rec[c]["e90_window_vs_saved"] for c in rec]}; '
              f'spread over odors {[rec[c]["spread_over_odors"] for c in rec]}; {time.time() - t1:.0f} s', flush=True)
        del res
        gc.collect()
    df = collect_levels()
    _write_json(out_dir() / 'provenance_levels.json', dict(part='levels', script=H.script_path(__file__),
        model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), runs=[list(r) for r in runs()],
        settings=dict(SETTLED_PRE_WINDOW=SETTLED_PRE_WINDOW, PRE_WINDOW=E90.PRE_WINDOW, LEVEL_ODORS=LEVEL_ODORS,
                      LEVEL_TOL=LEVEL_TOL), seconds=round(time.time() - t0, 1)))
    return df


def load_levels(part: str, tag: str) -> Optional[dict]:
    """A run's saved levels (`part_levels`): {cond 'on' / 'off': {'levels': {version: {stage: array, 'apl': float}}, ...}}."""
    path = local_dir() / 'levels' / f'{part}_{tag}.json'
    if not path.exists():
        return None
    rec = json.loads(path.read_text())
    for c in rec:
        for v in VERSIONS:
            for k in STAGES:
                rec[c]['levels'][v][k] = np.asarray(rec[c]['levels'][v][k], dtype=float)
    return rec


def collect_levels() -> pd.DataFrame:
    """preodor_levels.parquet (every unit's level, long: part, tag, cond, stage, unit, e90, settled, saved) and
    preodor_levels.csv (per run, condition and stage: the means over units at both levels, the largest change, the
    largest difference of E90's window from E90's saved means, the identity of the short run and the number of claws
    whose level falls to 0)."""
    long, rows = [], []
    for part, tag, _ in runs():
        rec = load_levels(part, tag)
        if rec is None:
            continue
        for c, cc in E90.CONDS:
            r = rec[c]
            for k in STAGES + ('apl',):
                e, s = np.atleast_1d(r['levels']['e90'][k]), np.atleast_1d(r['levels']['settled'][k])
                if k != 'apl':
                    long.append(pd.DataFrame(dict(part=part, tag=tag, cond=cc, stage=k, unit=np.arange(len(e)), e90=e,
                                                  settled=s)))
                rows.append(dict(part=part, tag=tag, cond=cc, stage=k, n_units=int(len(e)), e90_mean=float(e.mean()),
                                 settled_mean=float(s.mean()), change_mean=float((s - e).mean()),
                                 change_max_abs=float(np.abs(s - e).max()),
                                 n_zero_e90=int((e <= 0).sum()), n_zero_settled=int((s <= 0).sum()),
                                 e90_window_vs_saved=float(r['e90_window_vs_saved'][k]),
                                 spread_over_odors=float(r['spread_over_odors']), short_run_exact=bool(r['identity']['exact']),
                                 n_steps_e90=r['n_steps']['e90'], n_steps_settled=r['n_steps']['settled']))
    if not rows:
        return pd.DataFrame()
    pd.concat(long, ignore_index=True).to_parquet(out_dir() / 'preodor_levels.parquet')
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'preodor_levels.csv', index=False)
    return df


def viki_check() -> pd.DataFrame:
    """viki's reconstructed settled levels of part B's runs (results/preodor_startup_check_viki/e90i0/, made from APL's
    saved traces of E98 B1 / E101 A3; her I0_<wiring>_<cond>_<stage>_pre_settled.parquet) against the short runs': the
    largest absolute difference per run and stage (where her files exist)."""
    d = H.shared_dir().joinpath(*VIKI_DIR)
    rows = []
    for w in E90.ALL_WIRINGS:
        rec = load_levels('B', w)
        if rec is None:
            continue
        for c, cc in E90.CONDS:
            for k in STAGES:
                f = d / f'I0_{w}_{cc}_{k}_pre_settled.parquet'
                if not f.exists():
                    continue
                v = pd.read_parquet(f)[k].to_numpy(dtype=float)
                s = rec[c]['levels']['settled'][k]
                rows.append(dict(wiring=w, cond=cc, stage=k, n=int(len(v)),
                                 max_abs_diff=float(np.max(np.abs(v - s))) if len(v) == len(s) else np.nan))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the readings at both levels

def stage_at(part: str, tag: str, cc: str, version: str) -> dict:
    """E90's saved stage dict of a run, with its pre-odor means replaced by the run's settled levels ('settled') or as
    saved ('e90'); APL's pre-odor mean likewise."""
    st = _saved_stage(part, tag, cc)
    if version == 'e90':
        return st
    rec = load_levels(part, tag)
    if rec is None:
        raise RuntimeError(f'{part} {tag}: no levels (run --part levels)')
    c = dict((cc_, c_) for c_, cc_ in E90.CONDS)[cc]
    lv = rec[c]['levels'][version]
    st = dict(st)
    for k in STAGES:
        x = st[f'{k}_pre']
        st[f'{k}_pre'] = pd.DataFrame(np.repeat(np.asarray(lv[k])[:, None], x.shape[1], axis=1), index=x.index,
                                      columns=x.columns)
    a = st['apl'].copy()
    a['pre'] = float(lv['apl'])
    st['apl'] = a
    return st


def _read(label: str, st: Dict[str, dict], units: str, seed: int) -> List[dict]:
    """E90's readings of one run pair (APL on / off): bouton and claw stages, the PN stage and the somata, as E90's
    `model_pairs`, `part_b5score` and `part_d1` read them (the calibrated readout's K on the run's own APL-on counts)."""
    from experiments import E83_cross_wiring_fit as E83
    rows = E90.read_wiring(label, st['on'], st['off'], E90.unit_meta(units), np.random.default_rng(seed))
    for c in ('on', 'off'):
        rows += E90.pn_rows(label, st[c], c)
    K = E83.readout_K({'on': (st['on']['responses'], st['on']['spike_counts'])})
    rows += E90.soma_rows(label, {c: st[c]['spike_counts'] for c in ('on', 'off')}, K)
    return rows


def pairs(kind: str, version: str, force: bool = False) -> pd.DataFrame:
    """The readings' pair table of part B ('B', E90's model_pairs), B5 ('B5', E90's b5_pairs without the 0.5 s rows) or
    D1 ('D1', E90's d1_pairs) at a level, cached under the local folder, with E90's seeds."""
    path = local_dir() / f'pairs_{kind}_{version}.parquet'
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = time.time()
    rows = []
    if kind == 'B':
        for i, w in enumerate(E90.ALL_WIRINGS):
            st = {c: stage_at('B', w, cc, version) for c, cc in E90.CONDS}
            rows += _read(w, st, w, E90.NOISE_SEED + i)
    elif kind == 'B5':
        st = {c: stage_at('B5', B5_TAG, cc, version) for c, cc in E90.CONDS}
        rows += _read(B5_TAG, st, 'real', E90.NOISE_SEED + 100)
    elif kind == 'D1':
        d1 = json.loads((E90.out_dir() / 'd1.json').read_text())
        for i, m in enumerate(E90.D1_MODELS):
            st = {c: stage_at('D1', m, cc, version) for c, cc in E90.CONDS}
            rows += _read(m, st, d1[m]['wiring'], E90.NOISE_SEED + 200 + i)
    else:
        raise ValueError(kind)
    df = pd.DataFrame(rows)
    df.to_parquet(path)
    print(f'pairs {kind} {version}: {len(df)} rows in {time.time() - t0:.0f} s', flush=True)
    return df


E90_SAVED_PAIRS = {'B': 'model_pairs.parquet', 'B5': 'b5_pairs.parquet', 'D1': 'd1_pairs.parquet'}


def reproduces(kind: str, mp_e90: pd.DataFrame) -> Optional[bool]:
    """Whether the 'e90' pairs equal E90's saved pair table of `kind` exactly (None when it is not on this host)."""
    path = E90.local_dir() / E90_SAVED_PAIRS[kind]
    if not path.exists():
        return None
    old = pd.read_parquet(path)
    if kind == 'B5':
        old = old[old.wiring == B5_TAG]
    old = old.reset_index(drop=True)
    new = mp_e90[list(old.columns)].reset_index(drop=True) if set(old.columns) <= set(mp_e90.columns) else None
    if new is None or len(new) != len(old):
        return False
    try:
        pd.testing.assert_frame_equal(new, old, check_dtype=False, check_exact=True)
    except AssertionError:
        return False
    return True


def stage_table(mp: pd.DataFrame) -> pd.DataFrame:
    """E90's stages_model.csv from a pair table: the population's and each wiring's group means (`part_bscore`)."""
    gw = E90._group_table(mp, E90.KEY + ['wiring'])
    gp = E90._group_table(E90.populations(mp), E90.KEY)
    gp['wiring'] = 'population'
    return pd.concat([gp, gw], ignore_index=True)


def verdict_of(sm: pd.DataFrame, steps: pd.DataFrame, nc: pd.DataFrame, effect_error: Dict[str, float]
               ) -> Tuple[dict, pd.DataFrame]:
    """E90's `part_c` on given tables (the flies' tables read from E90's results): the verdict per group on the primary
    subsets, on both subsets as readings and on the alternatives; the PN stage beside; whether D2 would be due."""
    passed = E90.check_passes(nc)
    vr = E90.verdict_rows(sm, steps, passed)
    both = {'bouton': ('shared', 'all'), 'claw': ('shared', 'all')}
    vb = E90.verdict_rows(sm, steps, passed, both)
    alt = {'claws_input': {'bouton': ('shared',), 'claw': ('input',)},
           'claws_shared': {'bouton': ('shared',), 'claw': ('shared',)},
           'all_three': {'bouton': E90.SUBSETS, 'claw': E90.SUBSETS}}
    va = {k: E90.verdict_rows(sm, steps, passed, v) for k, v in alt.items()}
    pr = E90.pn_verdict_rows(sm)
    rows = pd.concat([vr.assign(verdict_subsets='primary'), vb.assign(verdict_subsets='both'),
                      pr.assign(verdict_subsets='pn')] + [x.assign(verdict_subsets=k) for k, x in va.items()],
                     ignore_index=True)
    els = ('bouton', 'claw', 'step_on', 'step_off')
    out = dict(noise_check_passed=passed, effect_error=effect_error, groups={}, both_subsets={},
               alternatives={k: dict(groups={}) for k in alt})
    pn_res = {g: E90.element_result(pr[pr.group == g]) for g in E90.GROUPS}
    for g in E90.GROUPS:
        res = {el: E90.element_result(vr[(vr.group == g) & (vr.element == el)]) for el in els}
        out['groups'][g] = dict(verdict=E90.group_verdict(res, effect_error[g]), elements=res, pn=pn_res[g])
        rb = {el: E90.element_result(vb[(vb.group == g) & (vb.element == el)]) for el in els}
        out['both_subsets'][g] = dict(verdict=E90.group_verdict(rb, effect_error[g]), elements=rb)
        for k, x in va.items():
            ra = {el: E90.element_result(x[(x.group == g) & (x.element == el)]) for el in els}
            out['alternatives'][k]['groups'][g] = dict(verdict=E90.group_verdict(ra, effect_error[g]), elements=ra)
    out['d2_due'] = bool(any(v['verdict'] == 'upstream' and v['elements']['bouton']['differs'] for v in out['groups'].values())
                         or any(v['differs'] for v in pn_res.values()))
    return out, rows


def steps_by_reading_of(steps: pd.DataFrame) -> pd.DataFrame:
    """E90's `steps_by_reading` (the addendum's steps by soma reading, session z claws, drive alone) on a steps table."""
    gs = E90._group_table(steps[(steps.share == 0.0) & (steps.reading == 'session_z')], ['step', 'soma', 'subset', 'wiring'])
    rows = []
    for step in ('on', 'off'):
        for soma in E90.STEP_SOMA_READINGS:
            fl = E90.flies_step_group(E90.FLY_STEP_OF[(step, soma.split('_')[1])], 'noise_removed')
            for sub in ('input', 'all', 'shared'):
                x = gs[(gs.step == step) & (gs.soma == soma) & (gs.subset == sub)]
                pop, wir = x[x.wiring == 'population'].iloc[0], x[x.wiring != 'population']
                for g in E90.GROUPS:
                    r = E90.differs(float(pop[g]), wir[g].to_numpy(), float(fl.loc[g, 'mean']), float(fl.loc[g, 'lo']),
                                    float(fl.loc[g, 'hi']))
                    rows.append(dict(step=step, soma=soma, claws=sub, group=g, model=float(pop[g]),
                                     model_min=float(wir[g].min()), model_max=float(wir[g].max()),
                                     flies=float(fl.loc[g, 'mean']), flies_lo=float(fl.loc[g, 'lo']),
                                     flies_hi=float(fl.loc[g, 'hi']), diff=r['diff'], n_same_side=r['n_same_side'],
                                     differs=r['differs']))
    return pd.DataFrame(rows)


def d1_table(mp_b: pd.DataFrame, fc: pd.DataFrame) -> pd.DataFrame:
    """E90's d1_steps.csv (part_d1's table) from the C2 pairs (four wirings) and the firing-ceiling pairs."""
    c2 = mp_b[mp_b.wiring.isin(E90.D1_C2_WIRINGS)]
    out = []
    for fam, df in (('C2_S4', c2), ('fc_S4', fc)):
        steps = E90.model_steps_of(df)
        gt = E90._group_table(steps[steps.wiring == 'population'], ['step', 'soma', 'subset', 'share', 'reading'])
        gt.insert(0, 'family', fam)
        out.append(gt)
        stg = E90._group_table(E90.populations(df), E90.KEY)
        stg.insert(0, 'family', fam)
        stg.insert(1, 'step', 'stage')
        out.append(stg)
    return pd.concat(out, ignore_index=True)


def b5_tables(b5: pd.DataFrame, mp_b: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """E90's stages_b5.csv (B5's group means beside the real wiring's at the 0.5 s pulse) and the steps of both, from all,
    the input and the shared glomeruli's claws (the addendum's 2 s table: somata in every reading minus the APL-on claws
    under the session z-score, drive alone)."""
    both = pd.concat([b5, mp_b[mp_b.wiring == 'real']], ignore_index=True)
    g = E90._group_table(both, E90.KEY + ['wiring'])
    st = E90.model_steps_of(both)
    st = st[(st.wiring != 'population') & (st.share == 0.0) & (st.reading == 'session_z')]
    return g, E90._group_table(st, ['step', 'soma', 'subset', 'wiring'])


def pn_totals_cv(version: str) -> float:
    """Outcome (g)'s statistic of E90's `pn_stage_table` at a level: the CV over the 12 megamat odors of the real wiring's
    summed PN rate changes (APL on) on the glomeruli of the PN dendrites' table (al_transfer.csv)."""
    st = stage_at('B', 'real', 'on', version)
    a3 = pd.read_csv(E90.out_dir() / 'al_transfer.csv')
    ent = a3[(a3.section == 'entries') & (a3['map'] == 'linear') & (a3.stat == 'pn')]
    gl = sorted(set(ent.key.str.split('|').str[0]))
    X = E90.stage_units(st, 'pn').loc[gl, list(E90.ODORS)]
    return E90._cv(X.sum(axis=0))


def _silent_claws(mp: pd.DataFrame) -> pd.DataFrame:
    x = mp[(mp.stage == 'claw') & (mp.share == 0.0) & (mp.reading == 'session_z')]
    return x.groupby(['wiring', 'subset', 'cond'])['n_units'].first().reset_index()


def part_score() -> dict:
    """The readings at both levels and E90's tables at the settled level, with the comparisons (module docstring)."""
    t0 = time.time()
    o = out_dir()
    P = {(k, v): pairs(k, v) for k in ('B', 'B5', 'D1') for v in VERSIONS}
    repro = {k: reproduces(k, P[(k, 'e90')]) for k in ('B', 'B5', 'D1')}
    print(f"the 'e90' readings reproduce E90's saved pairs: {repro}", flush=True)
    old_verdict = json.loads((E90.out_dir() / 'verdict.json').read_text())
    ee = old_verdict['effect_error']          # the somata's errors in the effect: spike counts, not affected
    res, tabs = {}, {}
    for v in VERSIONS:
        mp = P[('B', v)]
        nc = E90.noise_check(mp, path=local_dir() / f'noise_check_{v}.csv')
        steps = E90.model_steps_of(mp)
        sm = stage_table(mp)
        verdict, vrows = verdict_of(sm, steps, nc, ee)
        sbr = steps_by_reading_of(steps)
        gsteps = E90._group_table(steps, ['step', 'soma', 'subset', 'share', 'reading', 'wiring'])
        b5g, b5s = b5_tables(P[('B5', v)], mp)
        d1 = d1_table(mp, P[('D1', v)])
        tabs[v] = dict(nc=nc, sm=sm, verdict=verdict, vrows=vrows, sbr=sbr, gsteps=gsteps, b5g=b5g, b5s=b5s, d1=d1)
        res[v] = verdict
    s = tabs['settled']
    s['sm'].to_csv(o / 'stages_model.csv', index=False)
    s['nc'].to_csv(o / 'noise_check.csv', index=False)
    s['gsteps'].to_csv(o / 'steps_groups.csv', index=False)
    s['vrows'].to_csv(o / 'verdict_rows.csv', index=False)
    s['sbr'].to_csv(o / 'steps_by_reading.csv', index=False)
    s['b5g'].to_csv(o / 'stages_b5.csv', index=False)
    s['b5s'].to_csv(o / 'b5_steps.csv', index=False)
    s['d1'].to_csv(o / 'd1_steps.csv', index=False)
    _write_json(o / 'verdict.json', dict(s['verdict'], pre_window=SETTLED_PRE_WINDOW))
    # the 'e90' pipeline against E90's saved tables (the verdict)
    e = tabs['e90']
    same_verdict = all(e['verdict']['groups'][g]['verdict'] == old_verdict['groups'][g]['verdict'] for g in E90.GROUPS)
    # comparisons
    G = list(E90.GROUPS)
    cs = e['sm'].merge(s['sm'], on=E90.KEY + ['wiring'], suffixes=('_e90', '_settled'))
    for g in G:
        cs[f'{g}_change'] = cs[f'{g}_settled'] - cs[f'{g}_e90']
    cs.to_csv(o / 'compare_stages.csv', index=False)
    cst = e['gsteps'].merge(s['gsteps'], on=['step', 'soma', 'subset', 'share', 'reading', 'wiring'],
                            suffixes=('_e90', '_settled'))
    for g in G:
        cst[f'{g}_change'] = cst[f'{g}_settled'] - cst[f'{g}_e90']
    cst.to_csv(o / 'compare_steps.csv', index=False)
    for name, key, on in (('compare_b5.csv', 'b5g', E90.KEY + ['wiring']), ('compare_b5_steps.csv', 'b5s',
                          ['step', 'soma', 'subset', 'wiring']),
                          ('compare_d1.csv', 'd1', ['family', 'step', 'soma', 'subset', 'share', 'reading', 'stage', 'cond'])):
        a, b = e[key], s[key]
        on = [c for c in on if c in a.columns]
        m = a.merge(b, on=on, suffixes=('_e90', '_settled'))
        for g in G:
            m[f'{g}_change'] = m[f'{g}_settled'] - m[f'{g}_e90']
        m.to_csv(o / name, index=False)
    vrow = []
    for sel, get in (('primary', lambda d: d['groups']), ('both', lambda d: d['both_subsets'])) + tuple(
            (k, (lambda k_: (lambda d: d['alternatives'][k_]['groups']))(k)) for k in ('claws_input', 'claws_shared', 'all_three')):
        for g in G:
            a, b = get(e['verdict'])[g], get(s['verdict'])[g]
            row = dict(subsets=sel, group=g, verdict_e90=a['verdict'], verdict_settled=b['verdict'],
                       changes=a['verdict'] != b['verdict'], effect_error=ee[g])
            for el in ('bouton', 'claw', 'step_on', 'step_off'):
                row[f'{el}_e90'] = f"{a['elements'][el]['n_differ']}/{a['elements'][el]['n_readings']}"
                row[f'{el}_settled'] = f"{b['elements'][el]['n_differ']}/{b['elements'][el]['n_readings']}"
            vrow.append(row)
    cv = pd.DataFrame(vrow)
    cv.to_csv(o / 'compare_verdict.csv', index=False)
    sc = _silent_claws(P[('B', 'e90')]).merge(_silent_claws(P[('B', 'settled')]), on=['wiring', 'subset', 'cond'],
                                              suffixes=('_e90', '_settled'))
    sc['became_silent'] = sc['n_units_e90'] - sc['n_units_settled']
    sc.to_csv(o / 'silent_claws.csv', index=False)
    vk = viki_check()
    if len(vk):
        vk.to_csv(o / 'viki_reconstruction_check.csv', index=False)
    summary = dict(reproduces_e90_pairs=repro, e90_pipeline_reproduces_e90_verdict=same_verdict,
                   noise_check_passed={v: res[v]['noise_check_passed'] for v in VERSIONS},
                   verdict={v: {g: res[v]['groups'][g]['verdict'] for g in G} for v in VERSIONS},
                   d2_due={v: res[v]['d2_due'] for v in VERSIONS}, n_verdicts_changed=int(cv.changes.sum()),
                   pn_totals_cv={v: pn_totals_cv(v) for v in VERSIONS},
                   viki_max_abs_diff=float(vk.max_abs_diff.max()) if len(vk) else None)
    _write_json(o / 'summary.json', summary)
    _write_json(o / 'provenance_score.json', dict(part='score', script=H.script_path(__file__), **H.code_commits(),
        inputs=dict(levels=str(local_dir() / 'levels'), e90_stage_runs='outputs/E90_model_stages/{B,B5,D1}/',
                    flies=str(E90.out_dir())), settings=dict(SETTLED_PRE_WINDOW=SETTLED_PRE_WINDOW, PRE_WINDOW=E90.PRE_WINDOW,
                    NOISE_SEED=E90.NOISE_SEED, N_NOISE_REPS=E90.N_NOISE_REPS), seconds=round(time.time() - t0, 1)))
    print(json.dumps(summary, indent=1), flush=True)
    return summary


# ---------------------------------------------------------------------------
# the comparison table

TG = E90.TABLE_GROUPS


def _cells(df: pd.DataFrame, fmt: str = '.3f', sign: bool = False) -> List[str]:
    spec = ('+' if sign else '') + fmt
    return [f"{E90._f(df[f'{g}_e90'], spec)} / {E90._f(df[f'{g}_settled'], spec)}" for g in TG]


def part_table() -> None:
    """comparison_table.md: E90's numbers that move, at E90's level / the settled level, beside the flies'."""
    o = out_dir()
    cs = pd.read_csv(o / 'compare_stages.csv')
    cst = pd.read_csv(o / 'compare_steps.csv')
    cv = pd.read_csv(o / 'compare_verdict.csv')
    lv = pd.read_csv(o / 'preodor_levels.csv')
    sc = pd.read_csv(o / 'silent_claws.csv')
    sm_flies = {k: E90.flies_group(k) for k in ('bouton_noise_removed', 'claw_noise_removed')}
    L = ['# E90 at the settled pre-odor level', '',
         f'Pre-odor level: E90\'s mean over {E90.PRE_WINDOW} s (before) / the settled mean over {SETTLED_PRE_WINDOW} s '
         '(after). Cells "a / b" are E90\'s value / the settled value. Distances are 1 - Pearson r across units, group '
         'means over the pairs; the model is the population of eight (C2\'s family, APL on unless stated).', '']
    L += ['## Pre-odor levels (mean over units)', '',
          '| run | cond | PNs | boutons | claws | APL | claws at 0 (E90 / settled) |', '|---|---|---:|---:|---:|---:|---:|']
    for (part, tag, cond), x in lv.groupby(['part', 'tag', 'cond'], sort=False):
        x = x.set_index('stage')
        cl = x.loc['claw']
        L.append(f"| {part} {tag} | {cond} | " + ' | '.join(
            f"{x.loc[k, 'e90_mean']:.3f} / {x.loc[k, 'settled_mean']:.3f}" for k in ('pn', 'bouton', 'claw', 'apl'))
            + f" | {int(cl['n_zero_e90'])} / {int(cl['n_zero_settled'])} |")
    L += ['', '## The verdict', '', '| subsets | group | verdict (E90 / settled) | boutons | claws | step with APL on | '
          'step without APL |', '|---|---|---|---:|---:|---:|---:|']
    for _, r in cv[cv.group.isin(TG)].iterrows():
        L.append(f"| {r.subsets} | {r.group} | {r.verdict_e90} / {r.verdict_settled}{' **changes**' if r.changes else ''} | "
                 + ' | '.join(f"{r[f'{el}_e90']} / {r[f'{el}_settled']}" for el in ('bouton', 'claw', 'step_on', 'step_off'))
                 + ' |')
    L += ['', '## Stages (population, APL on)', '', '| stage | ' + ' | '.join(E90.TABLE_GROUP_HEAD[g] for g in TG) + ' |',
          '|---|' + '---:|' * len(TG)]
    for stage, flies in (('bouton', 'bouton_noise_removed'), ('claw', 'claw_noise_removed')):
        fl = sm_flies[flies]
        L.append(f'| {stage}s, flies, noise removed | ' + ' | '.join(E90._flies_cell(fl, g) for g in TG) + ' |')
        x = cs[(cs.stage == stage) & (cs.cond == 'on') & (cs.wiring == 'population')]
        for sub in ('shared', 'all', 'input'):
            for share in ((0.0, 0.5) if stage == 'claw' else (0.0,)):
                for rd in ('session_z', 'rate', 'noisy_nr'):
                    y = x[(x.subset == sub) & (x.share == share) & (x.reading == rd)]
                    if len(y):
                        L.append(f'| model {stage}s, {sub}, share {share}, {rd} | ' + ' | '.join(_cells(y.iloc[0])) + ' |')
    L += ['', '## Steps (somata minus APL-on claws; population)', '',
          '| step | soma reading | claws | claw reading | ' + ' | '.join(E90.TABLE_GROUP_HEAD[g] for g in TG) + ' |',
          '|---|---|---|---|' + '---:|' * len(TG)]
    sbr = pd.read_csv(E90.out_dir() / 'steps_by_reading.csv')
    for step in ('on', 'off'):
        for soma in E90.STEP_SOMA_READINGS:
            fl = sbr[(sbr.step == step) & (sbr.soma == soma) & (sbr.claws == 'input')].set_index('group')
            L.append(f'| {step}, flies | {soma.split("_")[1]} | noise removed | | '
                     + ' | '.join(f"{fl.loc[g, 'flies']:+.3f} [{fl.loc[g, 'flies_lo']:.2f}, {fl.loc[g, 'flies_hi']:.2f}]"
                                  for g in TG) + ' |')
            for sub in ('all', 'input'):
                y = cst[(cst.step == step) & (cst.soma == soma) & (cst.subset == sub) & (cst.share == 0.0)
                        & (cst.reading == 'session_z') & (cst.wiring == 'population')]
                L.append(f'| {step}, model | {soma} | {sub} | session_z | ' + ' | '.join(_cells(y.iloc[0], sign=True)) + ' |')
    b5 = pd.read_csv(o / 'compare_b5.csv')
    b5s = pd.read_csv(o / 'compare_b5_steps.csv')
    L += ['', "## B5: the real wiring at the 0.5 s pulse and at the flies' 2 s pulse (APL on)", '',
          '| reading | pulse | ' + ' | '.join(E90.TABLE_GROUP_HEAD[g] for g in TG) + ' |', '|---|---|' + '---:|' * len(TG)]
    for stage, sub, rd in (('bouton', 'shared', 'session_z'), ('bouton', 'shared', 'rate'), ('claw', 'all', 'session_z'),
                           ('claw', 'input', 'session_z'), ('claw', 'shared', 'session_z')):
        for w, pulse in (('real', '0.5 s'), (B5_TAG, '2 s')):
            y = b5[(b5.stage == stage) & (b5.subset == sub) & (b5.share == 0.0) & (b5.cond == 'on') & (b5.reading == rd)
                   & (b5.wiring == w)]
            L.append(f'| {stage}s, {sub}, {rd} | {pulse} | ' + ' | '.join(_cells(y.iloc[0])) + ' |')
    for soma in E90.STEP_SOMA_READINGS:
        for w, pulse in (('real', '0.5 s'), (B5_TAG, '2 s')):
            y = b5s[(b5s.step == 'on') & (b5s.soma == soma) & (b5s.subset == 'input') & (b5s.wiring == w)]
            L.append(f'| step with APL on, input claws, {soma} | {pulse} | ' + ' | '.join(_cells(y.iloc[0], sign=True)) + ' |')
    d1 = pd.read_csv(o / 'compare_d1.csv')
    L += ['', "## D1: the firing-ceiling family beside C2's on the same four wirings (steps from all / the input "
          "glomeruli's claws; somata under the calibrated readout, Pearson; claws APL on, session z, drive alone)", '',
          '| family | step | claws | ' + ' | '.join(E90.TABLE_GROUP_HEAD[g] for g in TG) + ' |', '|---|---|---|' + '---:|' * len(TG)]
    for fam in ('C2_S4', 'fc_S4'):
        for step in ('on', 'off'):
            for sub in ('all', 'input'):
                y = d1[(d1.family == fam) & (d1.step == step) & (d1.soma == 'cal_pearson') & (d1.subset == sub)
                       & (d1.share == 0.0) & (d1.reading == 'session_z')]
                L.append(f'| {fam} | {step} | {sub} | ' + ' | '.join(_cells(y.iloc[0], sign=True)) + ' |')
    L += ['', '## Claws that become silent at the settled level (claws of all glomeruli, drive alone)', '',
          '| wiring | cond | units read (E90 / settled) |', '|---|---|---:|']
    for _, r in sc[sc.subset == 'all'].iterrows():
        L.append(f'| {r.wiring} | {r.cond} | {r.n_units_e90} / {r.n_units_settled} |')
    (o / 'comparison_table.md').write_text('\n'.join(L) + '\n')
    print(f"wrote {o / 'comparison_table.md'}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['levels', 'score', 'table'])
    ap.add_argument('--only', nargs='*', default=None, help='levels: only these tags (wirings or D1 models, real_pulse2)')
    a = ap.parse_args()
    {'levels': lambda: part_levels(a.only), 'score': part_score, 'table': part_table}[a.part]()


if __name__ == '__main__':
    main()
