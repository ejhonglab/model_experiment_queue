#!/usr/bin/env python3
"""E125: the part of the inhibition among KCs that outlasts the block: a pooled neuron or a KC's own feedback in place
of the matrix.

The entry is item E125 of experiment_priorities.md (proposed by viki 2026-11-21 14:27, moved into the list 14:33 under
the idle-host rule; planned by Fable subagents, claude-fable-5-2). It builds on E124 (experiments/E124_kckc_apl_borne.py:
its ladder, the hemibrain matrix's KC-KC scale with APL off at 0.25, 0.5 and 0.75 of the fit's), E123 (the runner at
fixed scales, the classes) and E122 (the pooled and per-type matrices). "The family" is C2's breadth-2 APL-effect fit
(E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`), the bouton base pinning
multibouton_claws='first' (E3's BOUTON through E10_joint_fit.base_kws); S is real, P1, P2, M1 and T is P3, P4, M2, M3.

The design (fixed scales, nothing fitted). In every arm the KC-KC scale is the wiring fit's, s, in both conditions, and
the two conditions get matrices of their own: with APL on (1 - d) HB + d X, with APL off d X (`condition_arrays`), HB the
hemibrain matrix (mean row sum 1, total weight N) and X a survivor with HB's total weight, so the hemibrain matrix at
(1 - d) of the fitted scale is present with APL on only and the survivor at d of the total weight in both. The two
conditions go to two simulators built with their own matrices (`two_condition_simulator`: the APL-off condition, which
`E77_input_screen.hold_rates` marks through `off_extra`, to the one built with the APL-off matrix). Both response rates
are re-held (the threshold bisected with APL off, then the tied APL scales with APL on). Tags S1_<arm>_<w>.
- Arm 1, a pooled survivor: X is E122's U (one neuron pooling all KCs, every entry total / N^2, the diagonal included):
  PO25, PO50, PO75 (d 0.25, 0.5, 0.75). 24 runs.
- Arm 2, a KC's own feedback: X is diagonal, each KC inhibited by its own spikes at its type's mean amount (the
  hemibrain matrix's mean row sum over the KC's type): DG25, DG50, DG75. 24 runs.
- Arm 3, one neuron per KC type: X is E122's TY (block-diagonal, each type's total incoming weight the hemibrain's): TY50
  (d 0.5). 8 runs.
- References, read and not run again: the family, KC-KC at zero (E110's part A; part B, the fitted reference, beside),
  E124's ladder at the same d (HBoff25, HBoff50, HBoff75; HBoff0 is E123's HBab), and E123's uniform neuron U (d 1, both
  conditions U at the fit's scale).

Parts (CLI `--part`):
- check (before any run): (i) the hemibrain matrix is E123's (sha256 against E123's HBab runs); (ii) with d 0 (Z0: on
  HB, off a matrix of zeros, at the fit's scale) the two-simulator runner reproduces E123's held point S1_HBab_real
  exactly (spike counts with APL on and off; the factors to 1e-12 relative); (iii) a diagonal matrix alone changes a
  KC's spike counts: on 'real' at the family's parameters with APL off, DG's survivor at d 0.5 against the zero matrix
  (the number of KC x odor spike counts that differ). check.json; (ii) failing stops everything, (iii) failing drops
  arm 2 (reported as not testable without a model change; `arm2_testable` false) and leaves arms 1 and 3.
- arm1, arm2, arm3 (simulations, resumable, per wiring; `--wirings`): rows in runs_1_<host>_<lane>.csv, held points in
  starts_1_<host>_<lane>.json, provenance (`H.model_provenance` per point) in provenance_1_<host>_<lane>.json; runs saved
  locally (outputs/E125_kckc_survivor/step1/) and in the model-output archive (data/model_outputs/E125_kckc_survivor/
  step1/).
- score (no simulation): readings.csv, runs.csv, verdict.json, comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- chain: check, arm1, arm2 (if testable), arm3 on the given wirings, then score and variants when every point is done.

The classes (fixed before any run; E123's, unchanged; calibrated pattern part, S's and T's populations apart, against
all recordings and the four halves A to D):
- "keeps the gain" (`E123.classify`): at or below the half-way point between the family and E110's part A on both
  populations, on every target (half-way points against all recordings 0.0554 on S and 0.0577 on T, checked in
  verdict.json).
- "does what the matrix does": keeps the gain at d 0.5 or 0.75 and is not more than WITHIN (0.003) above E124's rung at
  the same d (HBoff50, HBoff75) on both populations, on every target (the reading against all recordings only is
  recorded beside it). A form (pooled, a KC's own feedback) does what the matrix does if one of its arms at d 0.5 or
  0.75 does.
- A run with the threshold or APL at a bound of its bisection is a corner; an arm with corners on more than two
  wirings gets no class ('corner (no class)').
- The outcomes of the entry, read from the forms' flags (pooled PO, a KC's own feedback DG; TY beside): (a) one of the
  two does what the matrix does; (b) neither keeps the gain at any d where the ladder does; (c) only one of the two does
  what the matrix does, or only one keeps the gain at a d where the ladder does; (d) an arm is a corner. They are
  recorded as flags for the write-up, which reads them.
- No correlation chooses anything, no pair is weighted, and no threshold moves after a run is read.
Realism: every run's four bands, γ's response rate and off / on ratio, the largest spike counts, the threshold and the
APL->KC scale; a threshold under REALISM_THR_HZ (60 Hz) is recorded as a realism cost even inside the bands
(`thr_below_60`).

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it, and the model loader's hemibrain KC->KC counts and KC types; no other connectome table; no raw-only recording; no
olfsysm change. Nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E125_kckc_survivor.py --part chain
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E125_kckc_survivor.py --part arm1 --wirings P3 P4
    ./venv/bin/python -u experiments/E125_kckc_survivor.py --part score && ./venv/bin/python -u experiments/E125_kckc_survivor.py --part variants
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import time
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108
from experiments import E122_feedback_neuron as E122
from experiments import E123_kckc_form as E123
from experiments import E124_kckc_apl_borne as E124

EXPERIMENT = 'E125_kckc_survivor'
SCRIPT = 'experiments/E125_kckc_survivor.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E123.S_WIRINGS, E123.T_WIRINGS, E123.WIRINGS
TYPES = E123.TYPES
DOSES = (0.25, 0.5, 0.75)
FORMS = {'PO': 'one neuron pooling all KCs (E122\'s U)', 'DG': "a KC's own feedback (the diagonal, its type's mean amount)",
         'TY': "one neuron per KC type (E122's TY)"}
PARTS = {'arm1': tuple(f'PO{round(d * 100)}' for d in DOSES), 'arm2': tuple(f'DG{round(d * 100)}' for d in DOSES),
         'arm3': ('TY50',)}
ARMS = PARTS['arm1'] + PARTS['arm2'] + PARTS['arm3']
RUN_ORDER = ARMS
ARM_LABELS = {a: f"the hemibrain matrix at {1 - int(a[2:]) / 100:g} of the fit's scale with APL on only, and "
                 f"{FORMS[a[:2]]} at {int(a[2:]) / 100:g} of its total weight in both conditions" for a in ARMS}
ARM_LABELS['Z0'] = "d 0: the hemibrain matrix with APL on only, a matrix of zeros with APL off (E123's HBab)"
LADDER = {0.25: 'HBoff25', 0.5: 'HBoff50', 0.75: 'HBoff75'}          # E124's rungs (read)
READ_KINDS = ('family', 'E110_A', 'E110_B', 'HBoff0') + tuple(LADDER.values())
E123_U = 'U1'                                                       # E123's uniform neuron, d 1 (read; E123 kind 'U')
OFF_MARK = '_E125_apl_off'
WITHIN = 0.003
MATRIX_LIKE_DOSES = (0.5, 0.75)
CORNER_WIRINGS_MAX = 2
REALISM_THR_HZ = 60.0
DIAG_CHECK_ARM = 'DG50'
RUNGS_ENTRY = {'HBoff25': (0.0536, 0.0570), 'HBoff50': (0.0537, 0.0552), 'HBoff75': (0.0529, 0.0504)}
RUNGS_TOL = 5e-5
TARGETS, POPS, READOUTS = E123.TARGETS, E123.POPS, E123.READOUTS
PARAM_KEYS = E83.PARAM_KEYS
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths

def out_dir() -> Path:
    """<shared>/results/E125_kckc_survivor/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return E122._host()


def _now() -> str:
    return E122._now()


def lane_name(wirings: Sequence[str]) -> str:
    return '-'.join(wirings)


# ---------------------------------------------------------------------------
# the matrices

def form_of(arm: str) -> str:
    """'PO', 'DG', 'TY' (or 'Z' for Z0)."""
    return 'Z' if arm == 'Z0' else arm[:2]


def dose_of(arm: str) -> float:
    """d of an arm (PO25 -> 0.25; Z0 -> 0)."""
    return 0.0 if arm == 'Z0' else int(arm[2:]) / 100.0


def survivor_array(form: str, ref: np.ndarray, types: np.ndarray) -> np.ndarray:
    """The survivor X of a form with the hemibrain matrix's total weight (the module docstring): 'PO' E122's U, 'TY'
    E122's TY, 'DG' the diagonal with each KC's entry its type's mean row sum of `ref`; 'Z' zeros."""
    n = len(types)
    if form == 'PO':
        return E122.arm_array('U', ref, types)
    if form == 'TY':
        return E122.arm_array('TY', ref, types)
    if form == 'DG':
        diag = np.zeros(n)
        for t in TYPES:
            m = types == t
            if m.any():
                diag[m] = ref[m].sum() / m.sum()
        return np.diag(diag)
    if form == 'Z':
        return np.zeros((n, n))
    raise ValueError(f'unknown form {form!r}')


def condition_arrays(arm: str, ref: np.ndarray, types: np.ndarray) -> Dict[str, np.ndarray]:
    """{'on': (1 - d) ref + d X, 'off': d X} of an arm (Z0: {'on': ref, 'off': zeros})."""
    d = dose_of(arm)
    X = survivor_array(form_of(arm), ref, types)
    if arm == 'Z0':
        return {'on': ref.copy(), 'off': np.zeros_like(ref)}
    return {'on': (1.0 - d) * ref + d * X, 'off': d * X}


def matrix_path(arm: str, cond: str) -> Path:
    return local_dir('matrices') / f'kckc_{arm}_{cond}.parquet'


def write_matrices(arm: str) -> Dict[str, Tuple[Path, pd.DataFrame, str]]:
    """{'on' | 'off': (path, matrix, sha256 of its values)} for the model's KCs, written to
    outputs/.../matrices/kckc_<arm>_<cond>.parquet unless there with the same values (columns as strings)."""
    from experiments.E77_input_screen import frame_sha256
    if arm in _CACHE.get('mat', {}):
        return _CACHE['mat'][arm]
    idx = E122.model_kc_index()
    ref = E122.hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL))
    arrs = condition_arrays(arm, ref.to_numpy(dtype=float), E122.kc_types_of(idx))
    out = {}
    for c, a in arrs.items():
        M = pd.DataFrame(a, index=ref.index, columns=ref.columns)
        sha = frame_sha256(M)
        p = matrix_path(arm, c)
        ok = False
        if p.exists():
            old = pd.read_parquet(p)
            old.columns = pd.Index([int(x) for x in old.columns], name=old.columns.name)
            ok = frame_sha256(old) == sha
        if not ok:
            W = M.copy()
            W.columns = W.columns.astype(str)
            tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
            W.to_parquet(tmp)
            os.replace(tmp, p)
        out[c] = (p, M, sha)
    _CACHE.setdefault('mat', {})[arm] = out
    return out


# ---------------------------------------------------------------------------
# the runner with a matrix for each condition

def is_off(params: Mapping[str, float]) -> bool:
    """Whether a call is the APL-off condition (hold_rates' off_extra carries OFF_MARK); checked against the APL scales."""
    off = bool(params.get(OFF_MARK, 0.0))
    if off and (float(params['wAPLKC_scale']) != 0.0 or float(params['wAPLPN_scale']) != 0.0):
        raise RuntimeError('E125: an APL-off call with APL scales not 0')
    if not off and float(params['wAPLKC_scale']) == 0.0:
        raise RuntimeError('E125: a call with APL->KC 0 not marked as the APL-off condition')
    return off


def dispatch(sim_on: Callable, sim_off: Callable) -> Callable:
    """params -> (responses, spike counts): the APL-off condition (OFF_MARK) to `sim_off`, the rest to `sim_on`, the
    mark removed."""
    def simulate(params):
        q = {k: v for k, v in params.items() if k != OFF_MARK}
        return (sim_off if is_off(params) else sim_on)(q)
    return simulate


def two_condition_simulator(w: str, params: Mapping[str, float], M_on: pd.DataFrame, M_off: pd.DataFrame
                            ) -> Tuple[Callable, dict]:
    """(simulate, info): two of E83's simulators of wiring `w` at `params`, one built with `M_on` as its KC-KC matrix
    and one with `M_off` (both at the KC-KC scale the params give), behind `dispatch`; info is the APL-on simulator's,
    with 'kws_off' the other's scalar kwargs."""
    sim_on, info = E83.build_simulator(w, params=dict(params), extra_kws=dict(_wKCKC=M_on))
    sim_off, info_off = E83.build_simulator(w, params=dict(params), extra_kws=dict(_wKCKC=M_off))
    info = dict(info, kws_off=info_off['kws'])
    return dispatch(sim_on, sim_off), info


def held(arm: str, w: str, params: Mapping[str, float], label: str) -> dict:
    """One held point: E77_input_screen.hold_rates on the two-condition simulator of the arm's matrices (the APL-off
    condition marked by off_extra; the KC-KC scale the fit's in both)."""
    from experiments.E77_input_screen import hold_rates
    P = E83.wiring_parts(w)
    m = write_matrices(arm)
    res = hold_rates(lambda: two_condition_simulator(w, params, m['on'][1], m['off'][1]),
                     target_off=P['target_off'], target_on=P['target_on'], params=dict(params), label=label,
                     off_extra={OFF_MARK: 1.0})
    res['params_off'] = {k: v for k, v in res['params_off'].items() if k != OFF_MARK}
    return res


# ---------------------------------------------------------------------------
# runs

def _tables(prefix: str) -> List[Path]:
    return sorted(out_dir().glob(f'{prefix}_*.csv'))


def runs_table() -> pd.DataFrame:
    """Every row of every host and lane (runs_1_*.csv)."""
    fs = [pd.read_csv(p) for p in _tables('runs_1')]
    return pd.concat(fs, ignore_index=True, sort=False) if fs else pd.DataFrame(columns=['tag'])


def save_run(tag: str, o: Dict[str, tuple]) -> None:
    """The run's matrices, locally and in the model-output archive (which another host's score part reads)."""
    for d in (local_dir('step1'), archive_dir() / 'step1'):
        d.mkdir(parents=True, exist_ok=True)
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1])):
                p = d / f'{tag}_{cc}_{k}.parquet'
                tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
                x.to_parquet(tmp)
                os.replace(tmp, p)


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    for d in (local_dir('step1'), archive_dir() / 'step1'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def part_of(arm: str) -> str:
    return next(p for p, a in PARTS.items() if arm in a)


def run_row(tag: str, arm: str, w: str, res: dict, o: Dict[str, tuple], cc: dict, lane: str) -> dict:
    P = E83.wiring_parts(w)
    m = write_matrices(arm)
    row = dict(part=part_of(arm), tag=tag, arm=arm, form=form_of(arm), dose=dose_of(arm), wiring=w,
               kckc_scale=float(res['params']['wKCKC_scale']), kckc_scale_off=float(res['params_off']['wKCKC_scale']),
               fit_folder=P['folder'], fit_tag=P['tag'], matrix_on=m['on'][0].name, matrix_on_sha256=m['on'][2],
               matrix_off=m['off'][0].name, matrix_off_sha256=m['off'][2],
               weight_on=float(m['on'][1].to_numpy().sum()), weight_off=float(m['off'][1].to_numpy().sum()),
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
               target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), host=_host(), lane=lane)
    row.update(E83.realism(o))
    row.update(E122.type_rates(o))
    row.update(E123.bands(row))
    row['thr_below_60'] = bool(float(row['fixed_thr']) < REALISM_THR_HZ)
    row.update(cc)
    return row


def planned(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS) -> List[Tuple[str, str]]:
    """The (arm, wiring) points run here, in RUN_ORDER (arm 1, then arms 2 and 3)."""
    return [(a, w) for a in RUN_ORDER if a in arms for w in WIRINGS if w in wirings]


# ---------------------------------------------------------------------------
# check

def diag_effect(w: str = 'real', arm: str = DIAG_CHECK_ARM) -> dict:
    """Design item 5's second half: at the family's parameters on `w` with APL off, the spike counts with the arm's
    APL-off matrix (its survivor at d alone) against those with a matrix of zeros (the number that differ)."""
    P = E83.wiring_parts(w)
    p = dict(P['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    m = write_matrices(arm)['off'][1]
    counts = {}
    for name, M in (('survivor', m), ('zeros', pd.DataFrame(0.0, index=m.index, columns=m.columns))):
        sim, _ = E83.build_simulator(w, params=dict(P['p']), extra_kws=dict(_wKCKC=M))
        counts[name] = sim(p)[1].to_numpy(dtype=float)
        del sim
        gc.collect()
    diff = counts['survivor'] != counts['zeros']
    return dict(wiring=w, arm=arm, n_differing=int(diff.sum()), n_kcs_differing=int(diff.any(axis=1).sum()),
                total_spikes_survivor=float(counts['survivor'].sum()), total_spikes_zeros=float(counts['zeros'].sum()),
                diag_only=bool(np.count_nonzero(m.to_numpy() - np.diag(np.diag(m.to_numpy()))) == 0),
                effect=bool(diff.any()))


def part_check() -> dict:
    """(i) the hemibrain matrix is E123's; (ii) Z0 (d 0) on 'real' reproduces E123's S1_HBab_real exactly; (iii) a
    diagonal matrix alone changes a KC's spike counts. check.json (a failed one is kept beside as
    check_failed_<time>.json). `passes` is (i) and (ii); `arm2_testable` is (iii)."""
    path = out_dir() / 'check.json'
    res = E122._read_json(path)
    if res.get('passes'):
        return res
    if res:
        E122._write_json(out_dir() / f"check_failed_{res.get('checked', '').replace(' ', '_').replace(':', '')}.json", res)
    rt = E123.runs_table()
    s123 = set(rt.loc[rt['arm'] == 'HBab', 'matrix_sha256'])
    m = write_matrices('Z0')
    mats = dict(sha256=m['on'][2], e123=sorted(s123), same=s123 == {m['on'][2]},
                off_weight=float(m['off'][1].to_numpy().sum()))
    cc = E108._commits_clean()
    P = E83.wiring_parts('real')
    t0 = time.time()
    hr = held('Z0', 'real', P['p'], 'E125 check Z0 real')
    saved = E123.load_run('S1_HBab_real')
    if saved is None:
        raise RuntimeError("E125 check: E123's S1_HBab_real is neither local nor in the archive")
    rec = rt[rt['tag'] == 'S1_HBab_real'].iloc[-1]
    diff = {c: int((hr[c][1].to_numpy(dtype=float) != saved[c][1].to_numpy(dtype=float)).sum()) for c in ('on', 'off')}
    fac = dict(thr=(hr['thr_factor'], float(rec.thr_factor)), apl=(hr['apl_factor'], float(rec.apl_factor)))
    exact = all(v == 0 for v in diff.values()) and all(abs(a - b) <= 1e-12 * abs(b) for a, b in fac.values())
    held_point = dict(differing_counts=diff, factors=fac, exact=exact, seconds=time.time() - t0,
                      kckc_scale_off=float(hr['params_off']['wKCKC_scale']),
                      what="Z0 (d 0: on the hemibrain matrix, off a matrix of zeros, two simulators) on 'real' at the "
                           "family fit, against E123's S1_HBab_real (one simulator, KC-KC scale 0 with APL off)")
    del hr
    gc.collect()
    dg = diag_effect()
    res = dict(checked=_now(), host=_host(), matrices=mats, held_point=held_point, diagonal=dg, commits=cc)
    res['passes'] = bool(exact and mats['same'])
    res['arm2_testable'] = bool(dg['effect'] and dg['diag_only'])
    E122._write_json(path, res)
    print(f"E125 check: passes {res['passes']}; arm 2 testable {res['arm2_testable']}; matrix {mats['same']}; held "
          f"point {held_point}; diagonal {dg}", flush=True)
    if not res['passes']:
        raise RuntimeError('E125 check failed (check.json)')
    return res


def wait_for_check(wait_s: float = 3 * 3600, every: float = 60) -> dict:
    t0 = time.time()
    while True:
        res = E122._read_json(out_dir() / 'check.json')
        if res.get('passes'):
            return res
        if 'passes' in res and not res['passes']:
            raise RuntimeError('check.json: the check failed')
        if time.time() - t0 > wait_s:
            raise RuntimeError('no passing check.json yet')
        print(f'E125: waiting for check.json ({_now()})', flush=True)
        time.sleep(every)


def runnable_arms(check: Mapping) -> Tuple[str, ...]:
    """ARMS without arm 2 when the check found the diagonal without effect."""
    return tuple(a for a in ARMS if check.get('arm2_testable') or a not in PARTS['arm2'])


# ---------------------------------------------------------------------------
# arms 1 to 3

def part_runs(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> pd.DataFrame:
    """Arms 1 to 3 (the module docstring), resumable, one process per lane; arm 2 refused if not testable."""
    chk = wait_for_check()
    ok = runnable_arms(chk)
    dropped = [a for a in arms if a not in ok]
    if dropped:
        print(f'E125: {dropped} not run (the check found the diagonal without effect: arm 2 not testable)', flush=True)
    arms = [a for a in arms if a in ok]
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    starts = out_dir() / f'starts_1_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in planned(arms, wirings):
        tag = f'S1_{arm}_{w}'
        if tag in set(runs_table()['tag']):
            print(f'E125 {tag}: done', flush=True)
            continue
        E108._still_clean(cc)
        P = E83.wiring_parts(w)
        res = held(arm, w, P['p'], f'E125 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = run_row(tag, arm, w, res, o, cc, lane)
        E83._append_row(table, row)
        prov = E122._read_json(prov_path) or dict(part='step1', script=SCRIPT, points={})
        m = write_matrices(arm)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
            change=(f"arm {arm} ({ARM_LABELS[arm]}): KC-KC matrix (1 - d) HB + d X with APL on and d X with APL off "
                    f"(d {dose_of(arm):g}, X {FORMS[form_of(arm)]}), each its own simulator, at the fit's KC-KC scale; "
                    "fixed_thr re-held with APL off to the fit's APL-off count, then the tied APL scales to its APL-on "
                    "count (E77_input_screen.hold_rates)"),
            matrices={c: dict(path=str(m[c][0].relative_to(H.HERE)), sha256=m[c][2]) for c in ('on', 'off')},
            fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
        E122._write_json(prov_path, prov)
        st = E122._read_json(starts)
        st[tag] = dict({k: float(res['params'][k]) for k in PARAM_KEYS}, source=f'held point {tag}')
        E122._write_json(starts, st)
        print(f"E125 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/{res['target_off']} off); thr {row['fixed_thr']:.1f}; bands {row['in_bands']}; flag "
              f"'{row.get('realism_flag') or ''}'; {res.get('seconds', 0):.0f}s", flush=True)
        del res, o
        gc.collect()
    return runs_table()


def all_done(arms: Optional[Sequence[str]] = None) -> bool:
    have = set(runs_table()['tag'])
    return all(f'S1_{a}_{w}' in have for a, w in planned(ARMS if arms is None else arms))


# ---------------------------------------------------------------------------
# score

def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83.pair_distances), cached locally by the sha256 of its spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist' not in _CACHE:
        try:
            _CACHE['dist'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist'] = {}
    st = _CACHE['dist']
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in st:
        st[key] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    return st[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run: ({(kind, wiring): outputs}, runs table). Kinds: READ_KINDS from E124's collect (the family, E110's A
    and B, E123's HBab as HBoff0, E124's ladder), E123's U as 'U1', and this experiment's arms."""
    o124, r124 = E124.collect()
    outs = {k: v for k, v in o124.items() if k[0] in READ_KINDS}
    rows = [r for r in r124.to_dict('records') if r['kind'] in READ_KINDS]
    o123, r123 = E123.collect()
    for r in r123.to_dict('records'):
        if r['kind'] == 'U':
            outs[(E123_U, r['wiring'])] = o123[('U', r['wiring'])]
            rows.append(dict(r, kind=E123_U, dose=1.0, form='PO'))
    del o123, o124

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E122.type_rates(o))
        r.update(E123.bands(r))
        rows.append(r)

    for r in runs_table().drop_duplicates('tag', keep='last').to_dict('records'):
        o = load_run(r['tag'])
        if o is not None:
            add(r['arm'], r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/step1', source=f"E125 ({r['host']})",
                form=r['form'], dose=r['dose'], thr_factor=r['thr_factor'], apl_factor=r['apl_factor'],
                at_bound_thr=bool(r['at_bound_thr']), at_bound_apl=bool(r['at_bound_apl']),
                **{k: float(r[k]) for k in PARAM_KEYS})
    runs = pd.DataFrame(rows)
    runs['thr_below_60'] = runs['fixed_thr'].astype(float) < REALISM_THR_HZ
    return outs, runs


def readings(outs, runs) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E122's
    `scores_of`, the runs' means of E122.RUN_COLS ('run_<col>'), the number of runs outside the bands."""
    rows = []
    sets = (('S4', S_WIRINGS), ('T4', T_WIRINGS), ('all8', WIRINGS))
    for kind in dict.fromkeys(runs['kind']):
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].drop_duplicates('wiring', keep='last').set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in sets if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                d = dd[unit] if len(ws) == 1 else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, n_wirings=len(ws), readout=ro,
                                 n_outside_bands=int((~s['in_bands'].astype(bool)).sum()),
                                 **{f'run_{c}': float(s[c].astype(float).mean()) for c in E122.RUN_COLS if c in s.columns},
                                 **E122.scores_of(d)))
    _save_dist_cache()
    return pd.DataFrame(rows)


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    return E122._v(rd, kind, unit, col, readout)


def corners(runs: pd.DataFrame, arm: str) -> dict:
    """The wirings on which an arm's run has the threshold or APL at its bisection's bound."""
    sub = runs[runs['kind'] == arm]
    at = sub[sub['at_bound_thr'].eq(True) | sub['at_bound_apl'].eq(True)]
    ws = sorted(at['wiring'])
    return dict(corner_wirings=ws, corner=len(ws) > CORNER_WIRINGS_MAX)


def matrix_like(rd: pd.DataFrame, arm: str, keeps_gain: bool) -> dict:
    """"Does what the matrix does": d 0.5 or 0.75, keeps the gain, and not more than WITHIN above E124's rung at the same
    d on both populations, on every target (and against all recordings only, beside)."""
    d = dose_of(arm)
    if d not in MATRIX_LIKE_DOSES:
        return dict(matrix_like=None, matrix_like_all_only=None, rung=None)
    rung = LADDER[d]
    per = {t: {u: _v(rd, arm, u, f'{t}_pattern') - _v(rd, rung, u, f'{t}_pattern') for u in POPS} for t in TARGETS}
    within = all(v <= WITHIN for t in TARGETS for v in per[t].values())
    within_all = all(v <= WITHIN for v in per['all'].values())
    return dict(matrix_like=bool(keeps_gain and within), matrix_like_all_only=bool(keeps_gain and within_all), rung=rung,
                d_rung=per)


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """An arm's class: E123's (keeps the gain), "does what the matrix does", and the corner rule."""
    c = E123.classify(rd, runs, arm)
    if not c.get('complete'):
        return c
    c.update(corners(runs, arm))
    c.update(matrix_like(rd, arm, bool(c['keeps_gain'])))
    cls = 'does what the matrix does' if c['matrix_like'] else ('keeps the gain' if c['keeps_gain'] else 'neither')
    if c['corner']:
        c['cls_before_corner'] = cls
        cls = 'corner (no class)'
    c['cls'] = cls
    sub = runs[runs['kind'] == arm]
    c.update(form=form_of(arm), dose=dose_of(arm), thr_below_60=sorted(sub.loc[sub['thr_below_60'], 'wiring']),
             fixed_thr_mean=float(sub['fixed_thr'].mean()), wAPLKC_scale_mean=float(sub['wAPLKC_scale'].mean()),
             apl_factor_mean=float(sub['apl_factor'].mean()))
    return c


def rungs_check(rd: pd.DataFrame) -> dict:
    """E124's rungs (calibrated pattern part, all recordings) against the entry's numbers."""
    out = {}
    for k, (s, t) in RUNGS_ENTRY.items():
        got = (_v(rd, k, 'S4', 'all_pattern'), _v(rd, k, 'T4', 'all_pattern'))
        out[k] = dict(got=got, entry=(s, t), same=bool(abs(got[0] - s) <= RUNGS_TOL and abs(got[1] - t) <= RUNGS_TOL))
    return out


def outcomes(classes: Mapping[str, dict], ladder: Mapping[str, dict]) -> dict:
    """The entry's outcomes as flags from the forms' classes (the module docstring)."""
    def arms_of(form):
        return [a for a in ARMS if form_of(a) == form]

    ladder_keeps = {d: bool(ladder.get(k, {}).get('keeps_gain')) for d, k in LADDER.items()}
    forms = {}
    for f in ('PO', 'DG', 'TY'):
        cs = {a: classes.get(a, {}) for a in arms_of(f)}
        complete = all(c.get('complete') for c in cs.values())
        ok = {a: c for a, c in cs.items() if c.get('complete') and not c.get('corner')}
        forms[f] = dict(complete=complete, arms={a: c.get('cls') for a, c in cs.items()},
                        matrix_like=any(bool(c.get('matrix_like')) for c in ok.values()),
                        keeps_at_ladder_d=any(bool(c.get('keeps_gain')) and ladder_keeps.get(dose_of(a), False)
                                              for a, c in ok.items()),
                        corner_arms=[a for a, c in cs.items() if c.get('corner')])
    po, dg = forms['PO'], forms['DG']
    two = [x for x in (po, dg) if x['complete']]
    return dict(forms=forms, ladder_keeps_gain=ladder_keeps,
                a=any(x['matrix_like'] for x in two),
                b=bool(two) and not any(x['keeps_at_ladder_d'] for x in two),
                c=(len(two) == 2 and (po['matrix_like'] != dg['matrix_like']
                                      or po['keeps_at_ladder_d'] != dg['keeps_at_ladder_d'])),
                d=any(f['corner_arms'] for f in forms.values()),
                complete_forms=[f for f, x in forms.items() if x['complete']])


def arm_table(rd: pd.DataFrame, runs: pd.DataFrame, classes: Mapping[str, dict]) -> List[dict]:
    """Per arm and reference: the calibrated pattern part on S and T (all recordings), the rung beside, the share of the
    gain recovered, the threshold, APL->KC, the wirings at a bound and under 60 Hz, the class."""
    out = []
    kinds = ['family', 'E110_A', 'HBoff0'] + list(LADDER.values()) + [E123_U] + list(ARMS)
    for k in kinds:
        sub = runs[runs['kind'] == k]
        if len(sub) < len(WIRINGS):
            out.append(dict(kind=k, complete=False))
            continue
        r = dict(kind=k, complete=True, cls=classes.get(k, {}).get('cls'))
        for u in POPS + ('all8',):
            r[f'pattern_cal_{u}'] = _v(rd, k, u, 'all_pattern')
            r[f'gain_share_{u}'] = E124.gain_share(rd, k, u)
        r['mean_effect_cal_all8'] = _v(rd, k, 'all8', 'mean_effect')
        r['mean_effect_raw_all8'] = _v(rd, k, 'all8', 'mean_effect', 'raw')
        r['fixed_thr_mean'] = float(sub['fixed_thr'].astype(float).mean())
        r['wAPLKC_scale_mean'] = float(sub['wAPLKC_scale'].astype(float).mean())
        r['n_thr_below_60'] = int(sub['thr_below_60'].astype(bool).sum())
        if 'at_bound_apl' in sub:
            r['n_at_bound'] = int((sub['at_bound_apl'].eq(True)
                                   | sub.get('at_bound_thr', pd.Series(False, index=sub.index)).eq(True)).sum())
        if k in ARMS and dose_of(k) in LADDER:
            r['rung'] = LADDER[dose_of(k)]
        out.append(r)
    return out


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    present = set(runs['kind'])
    classes = {a: classify(rd, runs, a) for a in ARMS if a in present}
    ladder = {k: E123.classify(rd, runs, k) for k in ('HBoff0',) + tuple(LADDER.values()) + (E123_U,) if k in present}
    rt = runs_table()
    bad = {}
    for a in set(rt.get('arm', [])):
        m = write_matrices(a)
        for c in ('on', 'off'):
            s = sorted(set(rt.loc[rt['arm'] == a, f'matrix_{c}_sha256']))
            if s != [m[c][2]]:
                bad[f'{a}_{c}'] = s
    if bad:
        raise RuntimeError(f"E125 score: an arm's matrix differs between hosts or from this host's: {bad}")
    chk = E122._read_json(out / 'check.json')
    at = arm_table(rd, runs, classes)
    v = dict(generated=_now(), host=_host(), check=dict(passes=chk.get('passes'), arm2_testable=chk.get('arm2_testable')),
             half_way=E124.half_way_check(rd), rungs=rungs_check(rd), classes=classes, references=ladder,
             outcomes=outcomes(classes, ladder), arm_table=at,
             fitted_reference_E110_B={u: _v(rd, 'E110_B', u, 'all_pattern') for u in POPS + ('all8',)}
             if 'E110_B' in present else None,
             rules=dict(within=WITHIN, matrix_like_doses=list(MATRIX_LIKE_DOSES), corner_wirings_max=CORNER_WIRINGS_MAX,
                        targets=list(TARGETS), populations=list(POPS), realism_thr_hz=REALISM_THR_HZ),
             n_runs=int(len(runs)), complete=all(classes.get(a, {}).get('complete') for a in runnable_arms(chk)))
    E122._write_json(out / 'verdict.json', v)
    kinds = ['family', 'E110_A', 'E110_B', 'HBoff0'] + list(LADDER.values()) + [E123_U] + list(ARMS)
    md = ['# E125 comparison tables (generated)', '',
          'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
          f"Outcomes: {json.dumps({k: v['outcomes'][k] for k in 'abcd'})}", '',
          f"Half-way points: {json.dumps(v['half_way'])}", '', f"E124's rungs: {json.dumps(v['rungs'])}", '',
          '## The arms beside the ladder', '',
          '| kind | rung | pattern S4 | pattern T4 | gain share S4 | gain share T4 | threshold (Hz) | APL->KC | '
          'at bound | thr < 60 Hz | class |', '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in at:
        if not r.get('complete'):
            md.append(f"| {r['kind']} | | incomplete |||||||| |")
            continue
        md.append(f"| {r['kind']} | {r.get('rung', '')} | {r['pattern_cal_S4']:.4f} | {r['pattern_cal_T4']:.4f} | "
                  f"{r['gain_share_S4']:.2f} | {r['gain_share_T4']:.2f} | {r['fixed_thr_mean']:.1f} | "
                  f"{r['wAPLKC_scale_mean']:.4g} | {r.get('n_at_bound', '')} | {r['n_thr_below_60']} | {r.get('cls') or ''} |")
    md.append('')
    for u in ('S4', 'T4', 'all8'):
        md += [f'## {u}', '', E123.comparison_table(rd, kinds, u), '']
    md += ['## Runs and bands', '', E123.runs_bands_table(runs), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    E122._write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT,
                                                         commits=H.code_commits()))
    print(json.dumps({**{a: c.get('cls') for a, c in classes.items()},
                      'outcomes': {k: v['outcomes'][k] for k in 'abcd'}}), flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every E125 run, absolute parameters) and model_pair_distances_long.csv (raw spike counts)."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    mine = runs[runs['folder'].astype(str).str.startswith(EXPERIMENT)]
    rt = runs_table().drop_duplicates('tag', keep='last').set_index('tag')
    for r in mine.itertuples():
        o = outs[(r.kind, r.wiring)]
        t = rt.loc[r.tag]
        rows.append(dict(tag=r.tag, variant=r.tag, part=part_of(r.kind), arm=r.kind, arm_label=ARM_LABELS[r.kind],
                         form=form_of(r.kind), dose=dose_of(r.kind), base='bouton', k=E124.K, wiring=r.wiring,
                         kind=f"C2's family's {r.wiring}: {ARM_LABELS[r.kind]}, at the fit's scales, both response "
                              'rates held',
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         kckc_matrix_on=f'matrices/{t["matrix_on"]}', kckc_matrix_off=f'matrices/{t["matrix_off"]}',
                         kckc_matrix_on_sha256=t['matrix_on_sha256'], kckc_matrix_off_sha256=t['matrix_off_sha256'],
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, kckc_scale_off=r.wKCKC_scale, rate_on=r.rate_on, rate_off=r.rate_off,
                         max_count_on=r.max_count_on, max_count_off=r.max_count_off, realism_flag=r.realism_flag,
                         thr_below_60=bool(r.thr_below_60), subdir='step1'))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def chain(wirings: Sequence[str] = WIRINGS, wait_others_s: float = 6 * 3600) -> None:
    """check, arms 1 to 3 on `wirings` (arm 2 only if testable), then (waiting up to `wait_others_s` for the other
    lanes) score and variants."""
    chk = part_check()
    for part in ('arm1', 'arm2', 'arm3'):
        part_runs(PARTS[part], wirings, lane=f'{part}-{lane_name(wirings)}')
    arms = runnable_arms(chk)
    t0 = time.time()
    while not all_done(arms):
        if time.time() - t0 > wait_others_s:
            print('E125: other lanes not done; scoring what there is', flush=True)
            break
        print(f'E125: waiting for the other lanes ({_now()})', flush=True)
        time.sleep(300)
    score()
    write_variants()
    print(f'E125 chain done {_now()}', flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', 'arm1', 'arm2', 'arm3', 'score', 'variants', 'chain', 'plan'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    a = ap.parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'check':
        part_check()
    elif a.part in PARTS:
        part_runs(PARTS[a.part], ws, lane=f'{a.part}-{lane_name(ws)}')
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'chain':
        chain(ws)
    elif a.part == 'plan':
        for arm, w in planned(ARMS, ws):
            print(arm, w)


if __name__ == '__main__':
    main()
