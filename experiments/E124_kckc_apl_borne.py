#!/usr/bin/env python3
"""E124: the inhibition among KCs as a part of APL, fitted with a KC-KC scale for each condition, and what sparse random
inhibition lacks.

The entry is item E124 of experiment_priorities.md (proposed by viki 2026-11-21 08:32, moved into the list 08:42 under
the idle-host rule; planned by Fable subagents, claude-fable-5-2). It builds on E123 (experiments/E123_kckc_form.py: the
runner at fixed scales, the classes, the matrices) and on E110 (KC-KC at zero, at fixed scales in part A and fitted in
part B). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's
`saved_fit`), the bouton base pinning multibouton_claws='first' (E3's BOUTON through E10_joint_fit.base_kws); S is real,
P1, P2, M1 and T is P3, P4, M2, M3.

The arms at fixed scales (E123's runner: the wiring's family fit, the arm's matrix at the fit's KC-KC scale, both response
rates re-held by `E77_input_screen.hold_rates`, the threshold bisected with APL off, then the tied APL scales with APL
on; tag S1_<arm>_<w>):
- Part A1, the dose with APL off: HBoff25, HBoff50, HBoff75: the hemibrain matrix, the KC-KC scale of the APL-off
  condition at 0.25, 0.5 and 0.75 of the fit's (`off_extra=dict(wKCKC_scale=f x the fit's)`), the APL-on one at the
  fit's. Dose 0 is E123's HBab and dose 1 the family (read, not run again).
- Part B, what the ladder lacks: TS16s0, TS16s1: 16 partners per KC drawn within the hemibrain's totals between KC types
  (`type_structured_partners`, seeds 0 and 1); K4s1, K4s2, K16s1, K16s2: E123's equal k-partner matrices
  (`E123.random_k_partners`) with seeds 1 and 2 (E123's K4 and K16 are seed 0, read as K4s0 and K16s0).
The type-structured matrix (fixed here before any run): for each postsynaptic type a (n_a KCs, the hemibrain's total
input R_a and its total T[a, b] from each presynaptic type b), the 16 n_a entries are split among the presynaptic types
in proportion to T[a, b] (largest remainders), dealt to type a's KCs 16 each in a random order, and each KC's partners of
type b drawn without replacement among type b's KCs (no self). Every entry of type a's rows is R_a / (16 n_a): each KC
has 16 partners with equal weights, every KC of a type the type's mean amount, the total weight the hemibrain's exactly,
and each block between types its total up to the rounding of the counts. (The ladder ignores the types and gives every
KC the same amount; the shuffles keep both and each KC's own amount; this keeps the types' totals and amounts only.)

Part A2, the fit (part A3 the same on E123's seed-0 shuffle, SH0, only if A2 keeps the gain): E10_joint_fit.py with
E110's / E72's part-B command (`-k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl --apl-readout both`, 150
evaluations) and `--off-kckc-scale 0` (the KC-KC scale of the blocked condition 0, the APL-on one free with the
threshold, APL->KC with KC->APL and bouton->APL tied, and APL->bouton under its cap), started at E123's held point of the
arm (S1_HBab_<w>; A3 S1_SH0ab_<w>, with `--kckc-matrix` of SH0's matrix) on each of the eight wirings, each its own
process and folder (results and outputs of E124_kckc_apl_borne/fits/<A2|A3>_<w>), claimed by mkdir so that wirings can be
split over hosts (`--wirings`).

Parts (CLI `--part`):
- check: (i) the hemibrain matrix and the seed-0 K4 / K16 matrices of this script are E123's (sha256 against E123's
  runs); (ii) the runner reproduces E123's held point S1_HBab_real exactly (spike counts with APL on and off, the factors
  to 1e-12 relative) as HBoff0 (dose 0). check.json. a1 / b wait for a passing check.json.
- a1, b (simulations, resumable, per wiring): rows in runs_1_<host>_<lane>.csv, held points in starts_1_<host>_<lane>.json,
  provenance (`H.model_provenance` per point) in provenance_1_<host>_<lane>.json; runs saved locally
  (outputs/E124_kckc_apl_borne/step1/) and in the model-output archive (data/model_outputs/E124_kckc_apl_borne/step1/).
- score (no simulation): readings.csv, runs.csv, verdict.json, comparison_table.md, a2_start_check.csv,
  provenance_score.json (below).
- a2, a3 (fits): `--wirings` takes a subset; a3 refuses unless verdict.json's A2 keeps the gain (or `--force`).
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- chain: check, a1 and b on the given wirings, then score and variants when every point is done.
- commands: prints the a2 / a3 lane commands per wiring and their E10_joint_fit.py commands.

The classes (fixed before any run; E123's, unchanged; `E123.classify`): calibrated pattern part, S's and T's populations
apart, against all recordings and the four halves; the reference "KC-KC at zero" is E110's part A at fixed scales
(half-way points against all recordings 0.0554 on S and 0.0577 on T, checked in verdict.json); the fitted reference
(E110's part B, the counted start) is read and reported beside every class and moves none.
- "keeps the gain": at or below the half-way point between the family and E110's part A on both populations, on every
  target; "equal": not more than 0.003 above the family on both populations and on 6 of 8 wirings, on every target, with
  every run in the four bands.
- A fit (A2, A3) that ends with APL->KC or the threshold at a bound of its fit (E10_joint_fit.fit_one's: the threshold
  within [1/4, 4] x the bouton base's, APL->KC within [1/30, 30] x the start's; at a bound = within BOUND_TOL, relative,
  on 1 or more wirings) is reported as a corner and gets no class.
- Part B: "the types make up the rest" if both type-structured draws keep the gain and no equal draw of 16 partners
  (K16s0, K16s1, K16s2) does; "a gap of draws" if the equal draws of one k (K4s0 to K4s2, or K16s0 to K16s2) straddle
  the half-way point (some keep the gain and some do not); "neither" otherwise (the KCs' unequal amounts remain).
- Part A1's reading (no class beyond each dose's): per dose 0, 0.25, 0.5, 0.75, 1 the pattern part on S and T, the share
  of the gain recovered ((E110 A - dose) / (E110 A - family)), the wirings with APL at the bracket's bound, the APL factor
  and the threshold.
- No correlation chooses anything, no pair is weighted, and no threshold moves after a run is read.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it, and the model loader's hemibrain KC->KC counts and KC types; no other connectome table; no raw-only recording; no
olfsysm change. Nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E124_kckc_apl_borne.py --part chain
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E124_kckc_apl_borne.py --part a2 --wirings real P1
    ./venv/bin/python -u experiments/E124_kckc_apl_borne.py --part score && ./venv/bin/python -u experiments/E124_kckc_apl_borne.py --part variants
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108
from experiments import E122_feedback_neuron as E122
from experiments import E123_kckc_form as E123

EXPERIMENT = 'E124_kckc_apl_borne'
SCRIPT = 'experiments/E124_kckc_apl_borne.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E123.S_WIRINGS, E123.T_WIRINGS, E123.WIRINGS
TYPES = E123.TYPES
DOSES = (0.25, 0.5, 0.75)
A1_ARMS = tuple(f'HBoff{round(f * 100)}' for f in DOSES)
TS_K, TS_SEEDS = 16, (0, 1)
B_ARMS = tuple(f'TS{TS_K}s{s}' for s in TS_SEEDS) + ('K4s1', 'K4s2', 'K16s1', 'K16s2')
ARMS = A1_ARMS + B_ARMS
RUN_ORDER = ARMS
# read from E123 (its arm names): dose 0, the seed-0 rungs, and the APL-borne shuffle (A3's start)
E123_READ = {'HBoff0': 'HBab', 'K4s0': 'K4', 'K16s0': 'K16', 'SH0ab': 'SH0ab'}
K_DRAWS = {4: ('K4s0', 'K4s1', 'K4s2'), 16: ('K16s0', 'K16s1', 'K16s2')}
ARM_LABELS = {**{f'HBoff{round(f * 100)}': f'the hemibrain matrix, KC-KC with APL off at {f:g} of the fit\'s scale'
                 for f in (0.0,) + DOSES},
              **{f'TS{TS_K}s{s}': f'{TS_K} partners per KC within the hemibrain\'s totals between KC types, equal weights '
                 f'per type, seed {s}' for s in TS_SEEDS},
              **{f'K{k}s{s}': f'{k} random partners per KC, equal weights and amounts, seed {s}' for k in (4, 16)
                 for s in (0, 1, 2)},
              'SH0ab': "E123's seed-0 shuffle, KC-KC 0 with APL off"}
FITS = ('A2', 'A3')
FIT_ARM = {'A2': 'HB', 'A3': 'SH0'}
FIT_START = {'A2': 'S1_HBab', 'A3': 'S1_SH0ab'}
FIT_PREFIX = f'{EXPERIMENT}/fits'
MAXFEV = 150
K = 2.0
OFF_KCKC = 0.0
BOUND_TOL = 0.01
THR_BOUNDS = (0.25, 4.0)            # E10_joint_fit.fit_one: fixed_thr within [1/4, 4] x the base's
APLKC_BOUNDS = (1 / 30, 30.0)       # wAPLKC_scale within [1/30, 30] x the start's
HALF_WAY_ENTRY = {'S4': 0.0554, 'T4': 0.0577}
HALF_WAY_TOL = 5e-5
TARGETS, POPS, READOUTS = E123.TARGETS, E123.POPS, E123.READOUTS
PARAM_KEYS = E83.PARAM_KEYS
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths

def out_dir() -> Path:
    """<shared>/results/E124_kckc_apl_borne/ (created)."""
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

def base_of(arm: str) -> str:
    """The matrix an arm uses: 'HB' for the dose arms (HBoff<d>) and A2, the arm itself otherwise."""
    return 'HB' if arm.startswith('HBoff') or arm == 'HB' else arm


def off_kckc_fraction(arm: str) -> Optional[float]:
    """The APL-off KC-KC scale as a fraction of the fit's for the dose arms (HBoff<d>: d / 100), None for the others."""
    if arm.startswith('HBoff'):
        return int(arm[len('HBoff'):]) / 100.0
    return None


def type_structured_partners(ref: np.ndarray, types: np.ndarray, k: int = TS_K, seed: int = 0) -> np.ndarray:
    """k partners per KC drawn within the hemibrain's totals between KC types (the module docstring): for each
    postsynaptic type a, the k n_a entries split among the presynaptic types in proportion to the block totals T[a, b]
    (largest remainders), dealt to type a's KCs k each in a random order, each KC's partners of a type drawn without
    replacement (no self); every entry of type a's rows R_a / (k n_a), R_a the type's total input."""
    rng = np.random.default_rng(seed)
    n = len(types)
    out = np.zeros((n, n))
    tl = list(dict.fromkeys(types))
    cols = {b: np.flatnonzero(types == b) for b in tl}
    for a in tl:
        rows = cols[a]
        T = np.array([ref[np.ix_(rows, cols[b])].sum() for b in tl])
        R = T.sum()
        if R <= 0:
            continue
        total = k * len(rows)
        q = total * T / R
        c = np.floor(q).astype(int)
        for i in np.argsort(-(q - c), kind='stable')[:total - c.sum()]:
            c[i] += 1
        labels = rng.permutation(np.repeat(np.arange(len(tl)), c))
        w = R / total
        for j, i in enumerate(rows):
            mine = labels[j * k:(j + 1) * k]
            for bi in np.unique(mine):
                m = int((mine == bi).sum())
                cand = cols[tl[bi]]
                cand = cand[cand != i]
                if m > len(cand):
                    raise ValueError(f'{m} partners of type {tl[bi]} for a KC of type {a}: only {len(cand)}')
                out[i, rng.choice(cand, size=m, replace=False)] = w
    return out


def arm_array(arm: str, ref: np.ndarray, types: np.ndarray) -> np.ndarray:
    """The arm's KC x KC matrix (post x pre) with the hemibrain matrix's total weight."""
    b = base_of(arm)
    if b == 'HB':
        return ref.copy()
    if b == 'SH0':
        return E123.arm_array('SH0', ref, types)
    if b.startswith('TS'):
        k, s = b[2:].split('s')
        return type_structured_partners(ref, types, int(k), int(s))
    if b.startswith('K') and 's' in b:
        k, s = b[1:].split('s')
        return E123.random_k_partners(len(types), int(k), float(ref.sum()), int(s))
    raise ValueError(f'unknown arm {arm!r}')


def matrix_path(arm: str) -> Path:
    return local_dir('matrices') / f'kckc_{base_of(arm)}.parquet'


def write_matrix(arm: str) -> Tuple[Path, pd.DataFrame, str]:
    """The arm's matrix for the model's KCs, written to outputs/.../matrices/kckc_<base>.parquet unless there with the
    same values (columns as strings, as E10_joint_fit.py's --kckc-matrix reads them); (path, matrix, sha256 of its
    values, E77's `frame_sha256`)."""
    from experiments.E77_input_screen import frame_sha256
    key = base_of(arm)
    if key in _CACHE.get('mat', {}):
        return _CACHE['mat'][key]
    idx = E122.model_kc_index()
    ref = E122.hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL))
    M = pd.DataFrame(arm_array(key, ref.to_numpy(dtype=float), E122.kc_types_of(idx)), index=ref.index,
                     columns=ref.columns)
    sha = frame_sha256(M)
    p = matrix_path(arm)
    ok = False
    if p.exists():
        old = pd.read_parquet(p)
        old.columns = pd.Index([int(c) for c in old.columns], name=old.columns.name)
        ok = frame_sha256(old) == sha
    if not ok:
        W = M.copy()
        W.columns = W.columns.astype(str)
        tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
        W.to_parquet(tmp)
        os.replace(tmp, p)
    _CACHE.setdefault('mat', {})[key] = (p, M, sha)
    return p, M, sha


def block_shares(M: np.ndarray, types: np.ndarray) -> dict:
    """The share of the total weight within type, and each (post, pre) block's share."""
    tot = M.sum()
    out = {}
    for a in TYPES:
        for b in TYPES:
            out[f'share_{E122.TYPE_COL[a]}_from_{E122.TYPE_COL[b]}'] = float(M[np.ix_(types == a, types == b)].sum() / tot)
    out['within_type_share'] = float(sum(M[np.ix_(types == t, types == t)].sum() for t in TYPES) / tot)
    return out


# ---------------------------------------------------------------------------
# runs

def _tables(prefix: str) -> List[Path]:
    return sorted(out_dir().glob(f'{prefix}_*.csv'))


def runs_table() -> pd.DataFrame:
    """Every fixed-scale row of every host and lane (runs_1_*.csv)."""
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


def off_extra_of(arm: str, params: Mapping[str, float]) -> Optional[Dict[str, float]]:
    """hold_rates' off_extra for an arm: the dose arms' APL-off KC-KC scale (f x the fit's), None for the others."""
    f = off_kckc_fraction(arm)
    return None if f is None else dict(wKCKC_scale=f * float(params['wKCKC_scale']))


def held(arm: str, w: str, params: Mapping[str, float], label: str) -> dict:
    """One held point (E123's runner: E77_input_screen.hold_rates on E83's simulator with the arm's matrix)."""
    from experiments.E77_input_screen import hold_rates
    P = E83.wiring_parts(w)
    _, M, _ = write_matrix(arm)
    return hold_rates(lambda: E83.build_simulator(w, params=dict(params), extra_kws=dict(_wKCKC=M)),
                      target_off=P['target_off'], target_on=P['target_on'], params=dict(params), label=label,
                      off_extra=off_extra_of(arm, params))


def run_row(tag: str, arm: str, w: str, res: dict, o: Dict[str, tuple], cc: dict, lane: str) -> dict:
    P = E83.wiring_parts(w)
    path, _, sha = write_matrix(arm)
    f = off_kckc_fraction(arm)
    row = dict(part='A1' if arm in A1_ARMS else 'B', tag=tag, arm=arm, wiring=w, off_fraction=f,
               kckc_scale_off=float(res['params_off']['wKCKC_scale']),
               fit_folder=P['folder'], fit_tag=P['tag'], matrix=path.name, matrix_sha256=sha,
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
               target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), host=_host(), lane=lane)
    row.update(E83.realism(o))
    row.update(E122.type_rates(o))
    row.update(E123.bands(row))
    row.update(cc)
    return row


def planned(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS) -> List[Tuple[str, str]]:
    """The (arm, wiring) points run here, in RUN_ORDER (A1's doses first, then part B)."""
    return [(a, w) for a in RUN_ORDER if a in arms for w in WIRINGS if w in wirings]


# ---------------------------------------------------------------------------
# check

def part_check() -> dict:
    """(i) the HB and seed-0 K4 / K16 matrices are E123's; (ii) HBoff0 on 'real' reproduces E123's S1_HBab_real exactly.
    check.json (a failed one is kept beside as check_failed_<time>.json)."""
    path = out_dir() / 'check.json'
    res = E122._read_json(path)
    if res.get('passes'):
        return res
    if res:
        E122._write_json(out_dir() / f"check_failed_{res.get('checked', '').replace(' ', '_').replace(':', '')}.json", res)
    rt = E123.runs_table()
    mats = {}
    for mine, theirs in (('HB', 'HBab'), ('K4s0', 'K4'), ('K16s0', 'K16')):
        s123 = set(rt.loc[rt['arm'] == theirs, 'matrix_sha256'])
        sha = write_matrix(mine)[2]
        mats[mine] = dict(sha256=sha, e123=sorted(s123), same=s123 == {sha})
    cc = E108._commits_clean()
    P = E83.wiring_parts('real')
    t0 = time.time()
    hr = held('HBoff0', 'real', P['p'], 'E124 check HBoff0 real')
    saved = E123.load_run('S1_HBab_real')
    if saved is None:
        raise RuntimeError("E124 check: E123's S1_HBab_real is neither local nor in the archive")
    rec = rt[rt['tag'] == 'S1_HBab_real'].iloc[-1]
    diff = {c: int((hr[c][1].to_numpy(dtype=float) != saved[c][1].to_numpy(dtype=float)).sum()) for c in ('on', 'off')}
    fac = dict(thr=(hr['thr_factor'], float(rec.thr_factor)), apl=(hr['apl_factor'], float(rec.apl_factor)))
    exact = all(v == 0 for v in diff.values()) and all(abs(a - b) <= 1e-12 * abs(b) for a, b in fac.values())
    res = dict(checked=_now(), host=_host(), matrices=mats,
               held_point=dict(differing_counts=diff, factors=fac, exact=exact, seconds=time.time() - t0,
                               off_kckc=float(hr['params_off']['wKCKC_scale']),
                               what="HBoff0 (dose 0) on 'real' at the family fit, against E123's S1_HBab_real"),
               commits=cc)
    res['passes'] = bool(exact and all(m['same'] for m in mats.values()))
    E122._write_json(path, res)
    print(f"E124 check: passes {res['passes']}; matrices {[(k, v['same']) for k, v in mats.items()]}; held point "
          f"{res['held_point']}", flush=True)
    if not res['passes']:
        raise RuntimeError('E124 check failed (check.json)')
    del hr
    gc.collect()
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
        print(f'E124: waiting for check.json ({_now()})', flush=True)
        time.sleep(every)


# ---------------------------------------------------------------------------
# parts a1 and b

def part_runs(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> pd.DataFrame:
    """Parts a1 and b (the module docstring), resumable, one process per lane."""
    wait_for_check()
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    starts = out_dir() / f'starts_1_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in planned(arms, wirings):
        tag = f'S1_{arm}_{w}'
        if tag in set(runs_table()['tag']):
            print(f'E124 {tag}: done', flush=True)
            continue
        E108._still_clean(cc)
        P = E83.wiring_parts(w)
        res = held(arm, w, P['p'], f'E124 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = run_row(tag, arm, w, res, o, cc, lane)
        E83._append_row(table, row)
        prov = E122._read_json(prov_path) or dict(part='step1', script=SCRIPT, points={})
        path, _, sha = write_matrix(arm)
        f = off_kckc_fraction(arm)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
            change=(f"KC-KC matrix of arm {arm} ({ARM_LABELS[arm]}, same total weight) at the fit's KC-KC scale"
                    + ('' if f is None else f"; KC-KC scale {f:g} x the fit's in the APL-off condition only")
                    + "; fixed_thr re-held with APL off to the fit's APL-off count, then the tied APL scales to its "
                      "APL-on count (E77_input_screen.hold_rates)"),
            matrix=dict(arm=arm, base=base_of(arm), path=str(path.relative_to(H.HERE)), sha256=sha),
            off_extra=off_extra_of(arm, P['p']), fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
        E122._write_json(prov_path, prov)
        st = E122._read_json(starts)
        st[tag] = dict({k: float(res['params'][k]) for k in PARAM_KEYS}, source=f'held point {tag}')
        E122._write_json(starts, st)
        print(f"E124 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/{res['target_off']} off); bands {row['in_bands']}; flag '{row.get('realism_flag') or ''}'"
              f"; {res.get('seconds', 0):.0f}s", flush=True)
        del res, o
        gc.collect()
    return runs_table()


def all_done() -> bool:
    have = set(runs_table()['tag'])
    return all(f'S1_{a}_{w}' in have for a, w in planned())


# ---------------------------------------------------------------------------
# parts a2 and a3 (fits)

def fit_name(part: str, w: str) -> str:
    return f'{part}_{w}'


def fit_dirs(part: str, w: str) -> Tuple[Path, Path]:
    """(results, outputs) folders of a fit (not created)."""
    n = fit_name(part, w)
    return H.shared_dir() / 'results' / FIT_PREFIX / n, H.HERE / 'outputs' / FIT_PREFIX / n


def fit_row(part: str, w: str) -> Optional[pd.Series]:
    p = fit_dirs(part, w)[0] / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    return f.iloc[-1] if len(f) else None


def load_fit(part: str, w: str) -> Optional[Tuple[pd.Series, Dict[str, tuple]]]:
    """(fits.csv row, saved runs) of a fit: this host's outputs, else the archive's; or None."""
    r = fit_row(part, w)
    if r is None:
        return None
    for d in (fit_dirs(part, w)[1], archive_dir() / 'fits' / fit_name(part, w)):
        o = E83._load_pair(d, str(r['tag']))
        if o is not None:
            return r, o
    return None


def starts_path() -> Path:
    return out_dir() / 'starts.json'


def e123_start(key: str) -> dict:
    """E123's held point `key` (S1_<arm>_<w>) from its starts_*.json files."""
    for p in sorted((H.shared_dir() / 'results' / E123.EXPERIMENT).glob('starts_*.json')):
        d = E122._read_json(p)
        if key in d:
            return dict(d[key], file=p.name)
    raise RuntimeError(f"E123 has no held point {key}")


def write_starts(part: str, wirings: Sequence[str] = WIRINGS) -> dict:
    """starts.json entries <part>_<w>: E123's held point of the arm (S1_HBab_<w> for A2, S1_SH0ab_<w> for A3), its APL-on
    parameters (KC-KC at the fit's), absolute."""
    from experiments import E110_no_kckc as E110
    entries = {}
    for w in wirings:
        e = e123_start(f'{FIT_START[part]}_{w}')
        entries[fit_name(part, w)] = E110.start_entry({k: float(e[k]) for k in PARAM_KEYS},
                                                      source=f"E123's held point {FIT_START[part]}_{w} ({e['file']})")
    return E110.merge_starts(starts_path(), entries)


def fit_command(part: str, w: str, python: Optional[str] = None) -> List[str]:
    """The E10_joint_fit.py command of one fit (E110's / E72's part-B command, `--off-kckc-scale 0`, the start; A3 also
    SH0's matrix by `--kckc-matrix`)."""
    from experiments import E110_no_kckc as E110
    py = str(H.HERE / 'venv' / 'bin' / 'python') if python is None else python
    cmd = [py, '-u', 'experiments/E10_joint_fit.py', '-k', f'{K:g}', '--aplpn-max', '3', '--rate-target', 'ctrl_frac',
           '--loss', 'apl', '--apl-readout', 'both', '--tag', f'_cap3_{part}', '--maxfev', str(MAXFEV),
           '--experiment', f'{FIT_PREFIX}/{fit_name(part, w)}'] + E110.wiring_args(w)
    if part == 'A3':
        cmd += ['--kckc-matrix', str(matrix_path('SH0')), '--kckc-matrix-tag', 'E123SH0']
    cmd += ['--off-kckc-scale', repr(OFF_KCKC), '--init-params', f'results/{EXPERIMENT}/starts.json:{fit_name(part, w)}']
    return cmd


def claim(name: str, reclaim: bool = False) -> bool:
    d = out_dir() / 'claims' / name
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    E122._write_json(d / 'claim.json', dict(name=name, host=_host(), pid=os.getpid(), time=_now()))
    return True


def release(name: str) -> None:
    shutil.rmtree(out_dir() / 'claims' / name, ignore_errors=True)


def a3_allowed() -> Tuple[bool, str]:
    """Part A3's gate: verdict.json's A2 complete, not a corner, and keeping the gain."""
    return a3_allowed_from(E122._read_json(out_dir() / 'verdict.json').get('fits', {}))


def part_fits(part: str, wirings: Optional[Sequence[str]] = None, dry_run: bool = False, reclaim: bool = False,
              force: bool = False) -> List[str]:
    """Parts a2 and a3: one E10_joint_fit.py process after another on `wirings` (default all eight), skipping those done
    or claimed; a row per fit in runs_fits_<host>.csv."""
    from experiments import E110_no_kckc as E110
    if part not in FITS:
        raise ValueError(part)
    ws = list(wirings or WIRINGS)
    if part == 'A3' and not dry_run:
        ok, why = a3_allowed()
        if not ok and not force:
            print(f'E124 A3: not run ({why})', flush=True)
            return []
    if dry_run:
        for w in ws:
            print(' '.join(fit_command(part, w)), flush=True)
        return []
    write_starts(part, ws)
    if part == 'A3':
        write_matrix('SH0')
    cc = E122.fit_clean()
    done = []
    for w in ws:
        name = fit_name(part, w)
        rd, od = fit_dirs(part, w)
        if fit_row(part, w) is not None:
            print(f'E124 {name}: done', flush=True)
            done.append(w)
            continue
        if not claim(name, reclaim):
            print(f'E124 {name}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            cmd = fit_command(part, w)
            print(f'E124 {name}: ' + ' '.join(cmd), flush=True)
            rec = E110.launch(cmd, rd, local_dir('logs') / f'{name}.log', cc, name)
            rec['launcher'] = SCRIPT
            E122._write_json(rd / 'launch.json', rec)
            if rec['returncode']:
                raise RuntimeError(f"E124 {name}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
            r = fit_row(part, w)
            if r is None:
                raise RuntimeError(f"E124 {name}: no row in its fits.csv (log {rec['log']})")
        except BaseException:
            release(name)
            raise
        # the fit's outputs in the archive, so that another host's score part reads them
        dst = archive_dir() / 'fits' / name
        dst.mkdir(parents=True, exist_ok=True)
        for p in od.glob(f"{r['tag']}_*.parquet"):
            shutil.copy2(p, dst / p.name)
        E83._append_row(out_dir() / f'runs_fits_{_host()}.csv', dict(
            name=name, part=part, wiring=w, tag=r['tag'], folder=f'{FIT_PREFIX}/{name}', host=_host(),
            seconds=rec['seconds'], command=' '.join(cmd), **E110.params_of_row(r), loss=float(r['loss']),
            n_evals=int(r['n_evals']), off_kckc_scale=float(r.get('off_kckc_scale', np.nan)), **cc))
        print(f"E124 {name}: loss {float(r['loss']):.4f}, KC-KC (APL on) {float(r['kc_kc_inh']):.3g}, "
              f"{rec['seconds'] / 60:.0f} min", flush=True)
        done.append(w)
    return done


def corner_of(row: Mapping, start: Mapping) -> dict:
    """Whether a fit ended with the threshold or APL->KC at a bound of its fit (BOUND_TOL, relative)."""
    thr = float(row['thr_factor'])
    apl = float(row['wAPLKC_scale']) / float(start['wAPLKC_scale'])
    at = lambda x, lo, hi: bool(x <= lo * (1 + BOUND_TOL) or x >= hi * (1 - BOUND_TOL))
    aplpn = float(row['aplpn_factor_vs_e3']) if 'aplpn_factor_vs_e3' in row else np.nan
    return dict(thr_factor_vs_base=thr, aplkc_factor_vs_start=apl, thr_at_bound=at(thr, *THR_BOUNDS),
                aplkc_at_bound=at(apl, *APLKC_BOUNDS), aplpn_factor_vs_base=aplpn,
                aplpn_at_cap=bool(np.isfinite(aplpn) and aplpn >= 3 * (1 - BOUND_TOL)))


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
    """Every run: ({(kind, wiring): outputs}, runs table). Kinds: 'family', 'E110_A', 'E110_B' (the counted fit), E123's
    runs read (HBoff0, K4s0, K16s0, SH0ab), this experiment's arms, and the fits 'A2', 'A3'."""
    from experiments import E110_no_kckc as E110
    outs, rows = {}, []
    ra = pd.concat([pd.read_csv(p) for p in sorted((H.shared_dir() / 'results' / 'E110_no_kckc').glob('runs_A_*.csv'))],
                   ignore_index=True).drop_duplicates('tag', keep='last').set_index('tag')

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E122.type_rates(o))
        r.update(E123.bands(r))
        rows.append(r)

    for w in WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], source='saved fit', **F['p'])
        a = E110.load_a(w)
        if a is not None:
            add('E110_A', w, a, tag=f'A_{w}', folder='E110_no_kckc/A', source='E110 part A (read)',
                **{k: float(ra.loc[f'A_{w}', k]) for k in PARAM_KEYS})
        b = E110.counted_b(w)
        if b is not None:
            o = E110.load_fit_outputs(b['name'], b['tag'])
            if o is not None:
                add('E110_B', w, o, tag=b['tag'], folder=f"E110_no_kckc/fits/{b['name']}",
                    source='E110 part B, counted start (read)', **E110.params_of_row(b['row']))
    r123 = E123.runs_table().drop_duplicates('tag', keep='last').set_index('tag')
    for kind, e in E123_READ.items():
        for w in WIRINGS:
            t = f'S1_{e}_{w}'
            o = E123.load_run(t) if t in r123.index else None
            if o is not None:
                r = r123.loc[t]
                add(kind, w, o, tag=t, folder=f'{E123.EXPERIMENT}/step1', source='E123 (read)',
                    off_fraction=0.0 if kind in ('HBoff0', 'SH0ab') else np.nan, thr_factor=float(r['thr_factor']),
                    apl_factor=float(r['apl_factor']), at_bound_apl=bool(r['at_bound_apl']),
                    **{k: float(r[k]) for k in PARAM_KEYS})
    for r in runs_table().drop_duplicates('tag', keep='last').to_dict('records'):
        o = load_run(r['tag'])
        if o is not None:
            add(r['arm'], r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/step1', source=f"E124 ({r['host']})",
                off_fraction=r.get('off_fraction'), thr_factor=r['thr_factor'], apl_factor=r['apl_factor'],
                at_bound_apl=bool(r['at_bound_apl']), **{k: float(r[k]) for k in PARAM_KEYS})
    for part in FITS:
        for w in WIRINGS:
            f = load_fit(part, w)
            if f is not None:
                r, o = f
                add(part, w, o, tag=str(r['tag']), folder=f'{FIT_PREFIX}/{fit_name(part, w)}', source=f'E124 fit {part}',
                    loss=float(r['loss']), n_evals=int(r['n_evals']), off_kckc_scale=float(r.get('off_kckc_scale', np.nan)),
                    **E110.params_of_row(r))
    return outs, pd.DataFrame(rows)


def readings(outs, runs) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E122's
    `scores_of`, the runs' means of E122.RUN_COLS ('run_<col>'), the number of runs outside the bands (E123's)."""
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


def half_way_check(rd: pd.DataFrame) -> dict:
    """The half-way points against all recordings (family and E110 part A) against the entry's 0.0554 / 0.0577."""
    out = {}
    for u, want in HALF_WAY_ENTRY.items():
        hw = (_v(rd, 'family', u, 'all_pattern') + _v(rd, 'E110_A', u, 'all_pattern')) / 2
        out[u] = dict(half_way=hw, entry=want, same=bool(abs(hw - want) <= HALF_WAY_TOL))
    return out


def gain_share(rd: pd.DataFrame, kind: str, unit: str) -> float:
    """The share of the gain recovered on `unit` (all recordings, readout): (E110 A - kind) / (E110 A - family)."""
    f, z, a = (_v(rd, k, unit, 'all_pattern') for k in ('family', 'E110_A', kind))
    return float((z - a) / (z - f))


def a1_reading(rd: pd.DataFrame, runs: pd.DataFrame, classes: Mapping[str, dict]) -> List[dict]:
    """Part A1's reading per dose (0: E123's HBab, 1: the family)."""
    out = []
    for f, kind in ((0.0, 'HBoff0'),) + tuple(zip(DOSES, A1_ARMS)) + ((1.0, 'family'),):
        sub = runs[runs['kind'] == kind]
        if len(sub) < len(WIRINGS):
            out.append(dict(dose=f, kind=kind, complete=False))
            continue
        row = dict(dose=f, kind=kind, complete=True, cls=classes.get(kind, {}).get('cls'),
                   keeps_gain=classes.get(kind, {}).get('keeps_gain'))
        for u in POPS + ('all8',):
            row[f'pattern_cal_{u}'] = _v(rd, kind, u, 'all_pattern')
            row[f'gain_share_{u}'] = gain_share(rd, kind, u)
        row['mean_effect_raw_all8'] = _v(rd, kind, 'all8', 'mean_effect', 'raw')
        row['loss_raw_all8'] = _v(rd, kind, 'all8', 'all_loss', 'raw')
        if kind != 'family':
            row['n_apl_at_bound'] = int(sub['at_bound_apl'].astype(bool).sum())
            row['apl_factor_mean'] = float(sub['apl_factor'].astype(float).mean())
            row['thr_factor_mean'] = float(sub['thr_factor'].astype(float).mean())
        row['fixed_thr_mean'] = float(sub['fixed_thr'].astype(float).mean())
        row['wAPLKC_scale_mean'] = float(sub['wAPLKC_scale'].astype(float).mean())
        row['rate_on_mean'] = float(sub['rate_on'].mean())
        out.append(row)
    return out


def b_reading(classes: Mapping[str, dict], rd: pd.DataFrame) -> dict:
    """Part B's reading (the module docstring)."""
    ts = [f'TS{TS_K}s{s}' for s in TS_SEEDS]
    need = ts + list(K_DRAWS[4]) + list(K_DRAWS[16])
    if not all(classes.get(a, {}).get('complete') for a in need):
        return dict(complete=False, missing=[a for a in need if not classes.get(a, {}).get('complete')])
    kg = {a: bool(classes[a]['keeps_gain']) for a in need}
    types_rest = all(kg[a] for a in ts) and not any(kg[a] for a in K_DRAWS[16])
    straddle = {k: len({kg[a] for a in ds}) == 2 for k, ds in K_DRAWS.items()}
    gap = any(straddle.values())
    reading = 'the types make up the rest' if types_rest else ('a gap of draws' if gap else 'neither')
    if types_rest and gap:
        reading = 'the types make up the rest; a gap of draws'
    vals = {a: {u: _v(rd, a, u, 'all_pattern') for u in POPS} for a in need}
    return dict(complete=True, reading=reading, types_make_up_the_rest=bool(types_rest), gap_of_draws=bool(gap),
                straddle=straddle, keeps_gain=kg, classes={a: classes[a]['cls'] for a in need}, pattern_cal_all=vals)


def fit_verdict(rd: pd.DataFrame, runs: pd.DataFrame, part: str) -> dict:
    """A fit part's class (E123.classify) with the corner rule, and the fitted parameters beside the family's."""
    c = E123.classify(rd, runs, part)
    if not c.get('complete'):
        return c
    st = E122._read_json(starts_path())
    corners = {}
    for w in WIRINGS:
        r = fit_row(part, w)
        corners[w] = corner_of(r, st[fit_name(part, w)])
    cw = sorted(w for w, x in corners.items() if x['thr_at_bound'] or x['aplkc_at_bound'])
    c.update(corner=bool(cw), corner_wirings=cw, corners=corners)
    if cw:
        c['cls_before_corner'] = c['cls']
        c['cls'] = 'corner (no class)'
    sub = runs[runs['kind'] == part].set_index('wiring')
    fam = runs[runs['kind'] == 'family'].set_index('wiring')
    c['params'] = {w: dict(fixed_thr=float(sub.loc[w, 'fixed_thr']), family_fixed_thr=float(fam.loc[w, 'fixed_thr']),
                           wAPLKC_scale=float(sub.loc[w, 'wAPLKC_scale']),
                           family_wAPLKC_scale=float(fam.loc[w, 'wAPLKC_scale']),
                           kckc_on=float(sub.loc[w, 'wKCKC_scale']), family_kckc=float(fam.loc[w, 'wKCKC_scale']))
                   for w in WIRINGS}
    c['vs_E110_B'] = {u: dict(fit=_v(rd, part, u, 'all_pattern'), E110_B=_v(rd, 'E110_B', u, 'all_pattern'))
                      for u in POPS + ('all8',)} if ((runs['kind'] == 'E110_B').sum() == len(WIRINGS)) else None
    return c


def a2_start_check(runs: pd.DataFrame) -> pd.DataFrame:
    """Each fit's first evaluation (its history) against the held run it started from (E123's): the response rates."""
    rows = []
    for part in FITS:
        start_kind = {'A2': 'HBoff0', 'A3': 'SH0ab'}[part]
        for w in WIRINGS:
            r = fit_row(part, w)
            if r is None:
                continue
            hp = fit_dirs(part, w)[0] / f"history_{r['tag']}.csv"
            if not hp.exists():
                continue
            h0 = pd.read_csv(hp).iloc[0]
            s = runs[(runs['kind'] == start_kind) & (runs['wiring'] == w)]
            if not len(s):
                continue
            s = s.iloc[0]
            rows.append(dict(part=part, wiring=w, start=s['tag'], rate_on_first=float(h0['rate_on']),
                             rate_on_held=float(s['rate_on']), rate_off_first=float(h0['rate_off']),
                             rate_off_held=float(s['rate_off']),
                             same=bool(abs(float(h0['rate_on']) - float(s['rate_on'])) < 1e-9
                                       and abs(float(h0['rate_off']) - float(s['rate_off'])) < 1e-9)))
    return pd.DataFrame(rows)


def comparison_table(rd: pd.DataFrame, kinds: Sequence[str], unit: str) -> str:
    ks = [k for k in kinds if ((rd['kind'] == k) & (rd['unit'] == unit)).any()]
    lines = [f'| statistic ({unit}) | ' + ' | '.join(ks) + ' |', '|---|' + '---:|' * len(ks)]
    for lab, col, ro in E123.STATS:
        if col is None:
            cells = [' / '.join(f'{_v(rd, k, unit, f"{h}_pattern", ro):.3f}' for h in E108.HALVES) for k in ks]
        else:
            cells = [f'{_v(rd, k, unit, col, ro):.3f}' for k in ks]
        lines.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    kinds_present = set(runs['kind'])
    arms = [a for a in ('HBoff0',) + A1_ARMS + B_ARMS + ('K4s0', 'K16s0', 'SH0ab') if a in kinds_present]
    classes = {a: E123.classify(rd, runs, a) for a in dict.fromkeys(arms)}
    fits = {p: fit_verdict(rd, runs, p) for p in FITS if p in kinds_present}
    sc = a2_start_check(runs)
    sc.to_csv(out / 'a2_start_check.csv', index=False)
    idx = E122.model_kc_index()
    types = E122.kc_types_of(idx)
    mstats = {}
    ref = write_matrix('HB')[1].to_numpy(dtype=float)
    for b in ['HB'] + [a for a in B_ARMS] + ['K4s0', 'K16s0']:
        M = write_matrix(b)[1].to_numpy(dtype=float)
        mstats[b] = dict(E123.matrix_stats(M, ref), **block_shares(M, types))
    rt = runs_table()
    shas = {a: sorted(set(rt.loc[rt['arm'] == a, 'matrix_sha256'])) for a in set(rt.get('arm', []))}
    bad = {a: s for a, s in shas.items() if len(s) != 1 or s[0] != write_matrix(a)[2]}
    if bad:
        raise RuntimeError(f"E124 score: an arm's matrix differs between hosts or from this host's: {bad}")
    v = dict(generated=_now(), host=_host(), half_way=half_way_check(rd), classes=classes,
             a1=a1_reading(rd, runs, classes), b=b_reading(classes, rd), fits=fits,
             fitted_reference_E110_B={u: _v(rd, 'E110_B', u, 'all_pattern') for u in POPS + ('all8',)}
             if 'E110_B' in kinds_present else None,
             a3_gate=dict(zip(('allowed', 'why'), a3_allowed_from(fits))),
             start_check_all_same=bool(sc['same'].all()) if len(sc) else None,
             matrices=mstats, matrix_sha256=shas,
             rules=dict(within=E123.WITHIN, wirings_needed=E123.N_WIRINGS_NEEDED, targets=list(TARGETS),
                        populations=list(POPS), bound_tol=BOUND_TOL),
             n_runs=int(len(runs)),
             complete=all(classes.get(a, {}).get('complete') for a in ARMS))
    E122._write_json(out / 'verdict.json', v)
    kinds = ['family', 'E110_A', 'E110_B', 'HBoff0'] + list(A1_ARMS) + ['A2', 'A3'] + list(B_ARMS) + ['K4s0', 'K16s0', 'SH0ab']
    md = ['# E124 comparison tables (generated)', '',
          'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
          'Fits: ' + ', '.join(f"{p} {c.get('cls', 'incomplete')}" for p, c in fits.items()), '',
          f"Part B: {v['b'].get('reading', 'incomplete')}", '',
          f"Half-way points: {json.dumps(v['half_way'])}", '']
    md += ['## Part A1, the dose with APL off', '', '| dose | kind | pattern S4 | pattern T4 | gain share S4 | gain share T4 '
           '| APL at bound | APL factor | threshold (Hz) | class |', '|---:|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in v['a1']:
        if not r.get('complete'):
            md.append(f"| {r['dose']:g} | {r['kind']} | incomplete |||||||| |")
            continue
        md.append(f"| {r['dose']:g} | {r['kind']} | {r['pattern_cal_S4']:.4f} | {r['pattern_cal_T4']:.4f} | "
                  f"{r['gain_share_S4']:.2f} | {r['gain_share_T4']:.2f} | {r.get('n_apl_at_bound', '')} | "
                  f"{r.get('apl_factor_mean', float('nan')):.3f} | {r['fixed_thr_mean']:.1f} | {r.get('cls')} |")
    md.append('')
    for u in ('S4', 'T4', 'all8'):
        md += [f'## {u}', '', comparison_table(rd, kinds, u), '']
    md += ['## The matrices', '', '| matrix | ' + ' | '.join(next(iter(mstats.values())).keys()) + ' |',
           '|---|' + '---:|' * len(next(iter(mstats.values())))]
    md += [f'| {b} | ' + ' | '.join(f'{x:.3f}' for x in s.values()) + ' |' for b, s in mstats.items()]
    md += ['', '## Runs and bands', '', E123.runs_bands_table(runs), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    E122._write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT,
                                                         commits=H.code_commits()))
    print(json.dumps({**{a: c.get('cls') for a, c in classes.items()}, **{p: c.get('cls') for p, c in fits.items()},
                      'B': v['b'].get('reading')}), flush=True)
    return v


def a3_allowed_from(fits: Mapping[str, dict]) -> Tuple[bool, str]:
    """(allowed, why) of part A3's gate from the fits' verdicts."""
    v = fits.get('A2', {})
    if not v.get('complete'):
        return False, 'A2 not complete'
    if v.get('corner'):
        return False, f"A2 a corner ({v.get('corner_wirings')})"
    if not v.get('keeps_gain'):
        return False, f"A2 {v.get('cls')!r}"
    return True, f"A2 {v.get('cls')}"


def write_variants() -> pd.DataFrame:
    """variants.csv (every E124 run and fit, absolute parameters) and model_pair_distances_long.csv (raw spike counts)."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    mine = runs[runs['folder'].astype(str).str.startswith(EXPERIMENT)]
    for r in mine.itertuples():
        o = outs[(r.kind, r.wiring)]
        is_fit = r.kind in FITS
        label = (f"fit with KC-KC 0 with APL off ({'hemibrain matrix' if r.kind == 'A2' else 'seed-0 shuffle'})" if is_fit
                 else ARM_LABELS[r.kind])
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.kind if is_fit else ('A1' if r.kind in A1_ARMS else 'B'),
                         arm=r.kind, arm_label=label, base='bouton', k=K, wiring=r.wiring,
                         kckc_scale_off=(OFF_KCKC if is_fit else
                                         (r.wKCKC_scale * off_kckc_fraction(r.kind) if off_kckc_fraction(r.kind) is not None
                                          else r.wKCKC_scale)),
                         kind=(f"C2's family's {r.wiring}: {label}" + ('' if is_fit else ", at the fit's scales, both "
                                                                        'response rates held')),
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         kckc_matrix=f"matrices/kckc_{base_of(FIT_ARM[r.kind] if is_fit else r.kind)}.parquet",
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag,
                         subdir=f'fits/{fit_name(r.kind, r.wiring)}' if is_fit else 'step1'))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def chain(wirings: Sequence[str] = WIRINGS, arms: Sequence[str] = ARMS, wait_others_s: float = 6 * 3600) -> None:
    """check, parts a1 and b on `wirings`, then (waiting up to `wait_others_s` for the other lanes) score and variants."""
    part_check()
    part_runs(arms, wirings)
    t0 = time.time()
    while not all_done():
        if time.time() - t0 > wait_others_s:
            print('E124: other lanes not done; scoring what there is', flush=True)
            break
        print(f'E124: waiting for the other lanes ({_now()})', flush=True)
        time.sleep(300)
    score()
    write_variants()
    print(f'E124 chain done {_now()}', flush=True)


def print_commands() -> None:
    lane = 'OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    for part in FITS:
        for w in WIRINGS:
            print(f'{lane} --part {part.lower()} --wirings {w}')
            print('    ' + ' '.join(fit_command(part, w, python='./venv/bin/python')))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('check', 'a1', 'b', 'score', 'a2', 'a3', 'variants', 'chain', 'plan', 'commands'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    ap.add_argument('--dry-run', action='store_true', help='(a2, a3) print the E10_joint_fit.py commands only')
    ap.add_argument('--reclaim', action='store_true', help="(a2, a3) take over a claim left by a process that died")
    ap.add_argument('--force', action='store_true', help="(a3) run without A2's gate")
    a = ap.parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'check':
        part_check()
    elif a.part == 'a1':
        part_runs(A1_ARMS, ws, lane=f"A1-{lane_name(ws)}")
    elif a.part == 'b':
        part_runs(B_ARMS, ws, lane=f"B-{lane_name(ws)}")
    elif a.part == 'score':
        score()
    elif a.part in ('a2', 'a3'):
        part_fits(a.part.upper(), ws, dry_run=a.dry_run, reclaim=a.reclaim, force=a.force)
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'chain':
        chain(ws)
    elif a.part == 'plan':
        for arm, w in planned(ARMS, ws):
            print(arm, w)
    elif a.part == 'commands':
        print_commands()


if __name__ == '__main__':
    main()
