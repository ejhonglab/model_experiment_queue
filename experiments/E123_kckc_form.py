#!/usr/bin/env python3
"""E123: the form of the inhibition among KCs: draws of the shuffle, random sparse matrices and the inhibition as a part
of APL.

The entry is item E123 of experiment_priorities.md (proposed by viki 2026-11-21 04:30, planned by a Fable subagent,
claude-fable-5-2; moved into the list 04:35 under the idle-host rule). It builds on E122's runner at fixed scales
(experiments/E122_feedback_neuron.py: the family, S and T, the matched total weight, the held point) and on E110's
part A (KC-KC 0 at fixed scales on all eight wirings). "The family" is C2's breadth-2 APL-effect fit (E67) on E72 part
B's eight wirings, each at its own fit (E83's `saved_fit`), the bouton base pinning multibouton_claws='first' as E122's
does; S is real, P1, P2, M1 and T is P3, P4, M2, M3.

The arms (KC x KC matrices, post x pre, the model's KCs, every one with the hemibrain matrix's total weight; `arm_array`):
- SH1, SH2, SH3: E122's shuffle within (post type, pre type) blocks with seeds 1, 2, 3 (SH0 is E122's seed 0).
- AK (the amount kept): each KC's row permuted among the KCs of each presynaptic type (the diagonal kept 0), so every KC
  keeps its row sum and its amount from each type and loses its partners; seed 0.
- K4, K16, K64, K256 (sparseness alone): k distinct random partners per KC (no self), equal weights, every KC's row
  sum the hemibrain's mean (total / N); seed 0. E122's U is the end of the ladder.
- HBab, SH0ab (the inhibition as a part of APL): the hemibrain matrix and the seed-0 shuffle with the KC-KC scale 0 in
  the APL-off condition only (`E77_input_screen.hold_rates(off_extra=dict(wKCKC_scale=0.0))`): the threshold is held
  with APL off and KC-KC off, then the tied APL scales with APL on and KC-KC at the fit's scale.
- References: the family (saved fits, all eight), KC-KC 0 at fixed scales (E110's part A, all eight, read), U and SH0
  on S (E122's step-1 runs, read) and on T (run here, as E122's step 1 on T's wirings).
Each run: the wiring's family fit, the arm's matrix at the fit's KC-KC scale, both response rates re-held
(`E77_input_screen.hold_rates`: the threshold bisected with APL off, then the tied APL scales with APL on).

Parts (CLI `--part`):
- check (one held point, before any new run): (i) the matrices of SH0 and U are E122's (sha256 of the values as in
  E122's runs_1.csv); (ii) the references read: each of E110's part A runs and E122's S1_U / S1_SH runs has the
  response counts its table records, and every family run the counts E110 part A held to; (iii) the APL-borne check:
  the held point with the option off (SH0 on 'real', no off_extra) reproduces E122's saved S1_SH_real run exactly
  (spike counts with APL on and off, threshold and APL factors). check.json. Part 1 waits for a passing check.json
  (another host's lane may start first).
- 1 (simulations, resumable, per wiring): `--wirings` and `--arms` (default all, in RUN_ORDER: the draws and the
  APL-borne arms first, then AK, the ladder from k 4 up, then the references U and SH0 on T); rows in
  runs_1_<host>_<lane>.csv (one file per process, so lanes on several hosts never write the same file), held points in
  starts_<host>_<lane>.json, provenance in provenance_1_<host>_<lane>.json, the runs' matrices saved locally
  (outputs/E123_kckc_form/step1/) and in the model-output archive (data/model_outputs/E123_kckc_form/step1/), which the
  score part reads, so a run made on another host is collected without copying. A point is done when its tag is in
  any runs_1_*.csv.
- score (no simulation): readings.csv (every kind; each wiring, S4, T4 and all8; raw and calibrated; E108's
  `unit_scores` against all recordings and halves A to D, the two ρ's, run means), runs.csv (one row per run with
  realism and each band), verdict.json (the classes below), comparison_table.md (pattern part with its amplitude and
  shape parts and the spread ratio under the calibrated readout, raw loss and mean effect, the rates, bands per run),
  step2_plan.json (the step-2 gate), provenance_score.json.
- 2: the gate only (step2_plan.json); the refits are not implemented yet (they need, for the APL-borne arm, the same
  option in E10_joint_fit.py, checked unchanged when off).
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- chain: check, 1 (the given wirings), and score + variants when every point of every wiring is done.

The classes (fixed before any run, from the entry; calibrated pattern part, S's and T's populations apart, against all
recordings and each half A to D):
- "keeps the gain": the arm's pattern part at or below the half-way point between the family's and E110 part A's on
  both populations, on every target.
- "equal": the arm's pattern part not higher than the family's by more than WITHIN (0.003) on both populations and on
  6 of 8 wirings, on every target, with every run in the four bands. ("Within 0.003" read one-sided, as E122's
  "replaces" was: an arm lower than the family by more is also within; the two-sided reading is recorded beside it as
  `equal_two_sided`.)
- the class of an arm: "equal" if equal, else "keeps the gain" if so, else "neither" (both flags recorded).
- the draws (SH0 to SH3) agree if all four have one class. The shuffle "costs response rate ρ" if the mean over the
  four seeds of the response rate ρ with APL on (each seed's mean over the eight wirings) is below the family's by more
  than RHO_COST (0.05), and below it on 3 of 4 seeds.
- No correlation chooses anything, no pair is weighted, no threshold moves after a run is read.
- Step 2's gate: at most two refit sets, only for arms classed "equal" and not flagged: the simplest random arm
  (smallest k first, then AK) and HBab.

Holds: only the model's hemibrain wiring and E72's seven re-draws, and the model loader's hemibrain KC->KC counts and
KC types (no other connectome table); no raw-only recording; no olfsysm change. Nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E123_kckc_form.py --part chain --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E123_kckc_form.py --part 1 --wirings P3 P4
    ./venv/bin/python -u experiments/E123_kckc_form.py --part score && ./venv/bin/python -u experiments/E123_kckc_form.py --part variants
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108
from experiments import E122_feedback_neuron as E122

EXPERIMENT = 'E123_kckc_form'
S_WIRINGS, T_WIRINGS, WIRINGS = E122.S_WIRINGS, E122.T_WIRINGS, E122.WIRINGS
TYPES = E122.TYPES
LADDER = (4, 16, 64, 256)
DRAWS = ('SH0', 'SH1', 'SH2', 'SH3')
APL_BORNE = ('HBab', 'SH0ab')
RUN_ORDER = ('SH1', 'SH2', 'SH3', 'HBab', 'SH0ab', 'AK') + tuple(f'K{k}' for k in LADDER) + ('U', 'SH0')
ARMS = RUN_ORDER
E122_READ = {'U': 'U', 'SH0': 'SH'}          # arms read from E122's step 1 on S (E122's arm names)
ARM_LABELS = {**{f'SH{s}': f'the hemibrain matrix shuffled within each type pair, seed {s}' for s in range(4)},
              'AK': "each KC's row permuted within each presynaptic type (its amounts kept, its partners lost)",
              **{f'K{k}': f'{k} random partners per KC, equal weights and amounts' for k in LADDER},
              'HBab': 'the hemibrain matrix, KC-KC 0 with APL off (inhibition carried by APL)',
              'SH0ab': 'the seed-0 shuffle, KC-KC 0 with APL off (inhibition carried by APL)',
              'U': "E122's uniform neuron (every entry 1/N)"}
WITHIN = 0.003
N_WIRINGS_NEEDED = 6
RHO_COST, RHO_SEEDS_NEEDED = 0.05, 3
TARGETS = ('all',) + tuple(E108.HALVES)
POPS = ('S4', 'T4')
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS
CHECK_WAIT_S = 3 * 3600
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths

def out_dir() -> Path:
    """<shared>/results/E123_kckc_form/ (created)."""
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
    """The matrix an arm uses ('HB', 'SH<s>', 'AK', 'K<k>', 'U'): the APL-borne arms use their base's."""
    return {'HBab': 'HB', 'SH0ab': 'SH0'}.get(arm, arm)


def is_apl_borne(arm: str) -> bool:
    return arm in APL_BORNE


def amount_kept_shuffle(ref: np.ndarray, types: np.ndarray, seed: int = 0) -> np.ndarray:
    """Each KC's row permuted within each presynaptic type: row i's entries from the KCs of type t permuted among those
    KCs (KC i itself excluded, its diagonal kept 0), so each KC keeps its row sum and its amount from each type and loses
    its partners."""
    rng = np.random.default_rng(seed)
    out = ref.copy()
    for t in TYPES:
        cols = np.flatnonzero(types == t)
        for i in range(ref.shape[0]):
            c = cols[cols != i]
            out[i, c] = ref[i, c][rng.permutation(len(c))]
    return out


def random_k_partners(n: int, k: int, total: float, seed: int = 0) -> np.ndarray:
    """k distinct random partners per KC (no self), equal weights total / (n k), so every KC's row sum is total / n."""
    if not 0 < k < n:
        raise ValueError(f'{k=} for {n=} KCs')
    rng = np.random.default_rng(seed)
    out = np.zeros((n, n))
    w = total / n / k
    for i in range(n):
        j = rng.choice(n - 1, size=k, replace=False)
        j[j >= i] += 1
        out[i, j] = w
    return out


def arm_array(arm: str, ref: np.ndarray, types: np.ndarray) -> np.ndarray:
    """The arm's KC x KC matrix (post x pre) as an array with the same total weight as `ref` (the module docstring)."""
    b = base_of(arm)
    if b == 'HB':
        return ref.copy()
    if b == 'U':
        return E122.arm_array('U', ref, types)
    if b.startswith('SH'):
        return E122.shuffle_within_blocks(ref, types, int(b[2:]))
    if b == 'AK':
        return amount_kept_shuffle(ref, types, 0)
    if b.startswith('K'):
        return random_k_partners(len(types), int(b[1:]), float(ref.sum()), 0)
    raise ValueError(f'unknown arm {arm!r}')


def matrix_path(arm: str) -> Path:
    return local_dir('matrices') / f'kckc_{base_of(arm)}.parquet'


def write_matrix(arm: str) -> Tuple[Path, pd.DataFrame, str]:
    """The arm's matrix for the model's KCs (index / columns of the normalized hemibrain matrix), written to
    outputs/.../matrices/kckc_<base>.parquet unless there with the same values; (path, matrix, sha256 of its values,
    E77's `frame_sha256`, as E122 records)."""
    from experiments.E77_input_screen import frame_sha256
    key = base_of(arm)
    if key in _CACHE.get('mat', {}):
        return _CACHE['mat'][key]
    idx = E122.model_kc_index()
    ref = E122.hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL))
    M = pd.DataFrame(arm_array(arm, ref.to_numpy(dtype=float), E122.kc_types_of(idx)), index=ref.index,
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
        W.to_parquet(p)
    _CACHE.setdefault('mat', {})[key] = (p, M, sha)
    return p, M, sha


def matrix_stats(M: np.ndarray, ref: np.ndarray) -> dict:
    """The matrix's total, inputs per KC, row-sum CV, reciprocal share of its nonzero pairs, and the row sums' Pearson r
    with the hemibrain's."""
    nz = M > 0
    rs, rr = M.sum(axis=1), ref.sum(axis=1)
    recip = float((nz & nz.T).sum() / nz.sum()) if nz.sum() else np.nan
    return dict(total=float(M.sum()), inputs_per_kc=float(nz.sum(axis=1).mean()), row_sum_cv=float(rs.std() / rs.mean()),
                reciprocal_share=recip,
                row_sum_r_hemibrain=float(np.corrcoef(rs, rr)[0, 1]) if rs.std() > 0 else np.nan)


# ---------------------------------------------------------------------------
# runs

def _tables(prefix: str) -> List[Path]:
    return sorted(out_dir().glob(f'{prefix}_*.csv'))


def runs_table() -> pd.DataFrame:
    """Every step-1 row of every host and lane (runs_1_*.csv)."""
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


def held(arm: str, w: str, params: Mapping[str, float], label: str) -> dict:
    """One held point (E77_input_screen.hold_rates on E83's simulator with the arm's matrix; for the APL-borne arms the
    KC-KC scale 0 in the APL-off condition)."""
    from experiments.E77_input_screen import hold_rates
    P = E83.wiring_parts(w)
    _, M, _ = write_matrix(arm)
    return hold_rates(lambda: E83.build_simulator(w, params=dict(params), extra_kws=dict(_wKCKC=M)),
                      target_off=P['target_off'], target_on=P['target_on'], params=dict(params), label=label,
                      off_extra=dict(wKCKC_scale=0.0) if is_apl_borne(arm) else None)


def run_row(tag: str, arm: str, w: str, res: dict, o: Dict[str, tuple], cc: dict, lane: str) -> dict:
    P = E83.wiring_parts(w)
    path, _, sha = write_matrix(arm)
    row = dict(part='step1', tag=tag, arm=arm, wiring=w, apl_borne=is_apl_borne(arm),
               kckc_scale_off=0.0 if is_apl_borne(arm) else float(res['params']['wKCKC_scale']),
               fit_folder=P['folder'], fit_tag=P['tag'], matrix=path.name, matrix_sha256=sha,
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
               target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}), host=_host(), lane=lane)
    row.update(E83.realism(o))
    row.update(E122.type_rates(o))
    row.update(bands(row))
    row.update(cc)
    return row


def bands(r: Mapping) -> dict:
    """The four bands of a run (E110's `bands_of`) and whether it is inside all."""
    from experiments import E110_no_kckc as E110
    b = E110.bands_of(float(r['rate_on']), float(r['rate_off']), float(r['rate_rho_on']), float(r['opk_on_tv']),
                      float(r['opk_off_tv']))
    return dict(b, in_bands=all(b.values()))


def planned(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS) -> List[Tuple[str, str]]:
    """The (arm, wiring) points run here, in RUN_ORDER (U and SH0 on S are E122's)."""
    return [(a, w) for a in RUN_ORDER if a in arms for w in wirings if w in WIRINGS
            and not (a in E122_READ and w in S_WIRINGS)]


# ---------------------------------------------------------------------------
# check

def part_check() -> dict:
    """(i) SH0's and U's matrices are E122's; (ii) the references read have their recorded counts; (iii) the held point
    with the APL-borne option off reproduces E122's S1_SH_real exactly (the module docstring). check.json."""
    from experiments import E110_no_kckc as E110
    path = out_dir() / 'check.json'
    res = E122._read_json(path)
    if res.get('passes'):
        return res
    r1 = pd.read_csv(H.shared_dir() / 'results' / E122.EXPERIMENT / 'runs_1.csv')
    ra = pd.concat([pd.read_csv(p) for p in sorted((H.shared_dir() / 'results' / 'E110_no_kckc').glob('runs_A_*.csv'))],
                   ignore_index=True)
    # (i)
    mats = {}
    for arm, e in E122_READ.items():
        sha122 = set(r1.loc[r1['arm'] == e, 'matrix_sha256'])
        sha = write_matrix(arm)[2]
        mats[arm] = dict(sha256=sha, e122=sorted(sha122), same=sha122 == {sha})
    # (ii)
    def n(df):
        return int((df.to_numpy(dtype=float) > 0.5).sum())
    refs = {}
    for w in WIRINGS:
        a = E110.load_a(w)
        row = ra[ra['tag'] == f'A_{w}'].iloc[-1]
        refs[f'E110_A_{w}'] = dict(n=(n(a['on'][0]), n(a['off'][0])), recorded=(int(row.n_on), int(row.n_off)))
        F = E83.saved_fit(w)['saved']
        refs[f'family_{w}'] = dict(n=(n(F['on'][0]), n(F['off'][0])),
                                   recorded=(int(row.target_on), int(row.target_off)))
    for r in r1[r1['arm'].isin(list(E122_READ.values()))].itertuples():
        o = E122.load_run('step1', r.tag)
        refs[f'E122_{r.tag}'] = dict(n=(n(o['on'][0]), n(o['off'][0])), recorded=(int(r.n_on), int(r.n_off)))
    for v in refs.values():
        v['same'] = tuple(v['n']) == tuple(v['recorded'])
    # (iii)
    cc = E108._commits_clean()
    P = E83.wiring_parts('real')
    t0 = time.time()
    hr = held('SH0', 'real', P['p'], 'E123 check SH0 real')
    saved = E122.load_run('step1', 'S1_SH_real')
    rec = r1[r1['tag'] == 'S1_SH_real'].iloc[0]
    diff = {c: int((hr[c][1].to_numpy(dtype=float) != saved[c][1].to_numpy(dtype=float)).sum()) for c in ('on', 'off')}
    fac = dict(thr=(hr['thr_factor'], float(rec.thr_factor)), apl=(hr['apl_factor'], float(rec.apl_factor)))
    exact = all(v == 0 for v in diff.values()) and all(a == b for a, b in fac.values())
    res = dict(checked=_now(), host=_host(), matrices=mats, references=refs,
               apl_borne_off=dict(differing_counts=diff, factors=fac, exact=exact, seconds=time.time() - t0,
                                  what="SH0 on 'real' at the family fit with off_extra None, against E122's S1_SH_real"),
               commits=cc)
    res['passes'] = bool(exact and all(m['same'] for m in mats.values()) and all(v['same'] for v in refs.values()))
    E122._write_json(path, res)
    print(f"E123 check: passes {res['passes']}; matrices {[(k, v['same']) for k, v in mats.items()]}; references "
          f"{sum(v['same'] for v in refs.values())}/{len(refs)}; APL-borne off {res['apl_borne_off']}", flush=True)
    if not res['passes']:
        raise RuntimeError('E123 check failed (check.json)')
    del hr
    gc.collect()
    return res


def wait_for_check(wait_s: float = CHECK_WAIT_S, every: float = 60) -> dict:
    """A passing check.json (written by the host that runs `check`), waiting up to `wait_s`."""
    t0 = time.time()
    while True:
        res = E122._read_json(out_dir() / 'check.json')
        if res.get('passes'):
            return res
        if 'passes' in res and not res['passes']:
            raise RuntimeError('check.json: the check failed')
        if time.time() - t0 > wait_s:
            raise RuntimeError('no passing check.json yet')
        print(f'E123: waiting for check.json ({_now()})', flush=True)
        time.sleep(every)


# ---------------------------------------------------------------------------
# step 1

def part_1(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> pd.DataFrame:
    """Step 1 (the module docstring), resumable, one process per lane."""
    wait_for_check()
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    starts = out_dir() / f'starts_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in planned(arms, wirings):
        tag = f'S1_{arm}_{w}'
        if tag in set(runs_table()['tag']):
            print(f'E123 {tag}: done', flush=True)
            continue
        E108._still_clean(cc)
        P = E83.wiring_parts(w)
        res = held(arm, w, P['p'], f'E123 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = run_row(tag, arm, w, res, o, cc, lane)
        E83._append_row(table, row)
        prov = E122._read_json(prov_path) or dict(part='step1', script='experiments/E123_kckc_form.py', points={})
        path, _, sha = write_matrix(arm)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
            change=(f"KC-KC matrix of arm {arm} ({ARM_LABELS[arm]}, same total weight) at the fit's KC-KC scale"
                    + ("; KC-KC scale 0 in the APL-off condition only" if is_apl_borne(arm) else "")
                    + "; fixed_thr re-held with APL off to the fit's APL-off count, then the tied APL scales to its "
                      "APL-on count (E77_input_screen.hold_rates)"),
            matrix=dict(arm=arm, base=base_of(arm), path=str(path.relative_to(H.HERE)), sha256=sha),
            fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
        E122._write_json(prov_path, prov)
        st = E122._read_json(starts)
        st[tag] = dict({k: float(res['params'][k]) for k in PARAM_KEYS}, source=f'step 1 held point {tag}')
        E122._write_json(starts, st)
        print(f"E123 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/{res['target_off']} off); bands {row['in_bands']}; flag '{row.get('realism_flag') or ''}'"
              f"; {res.get('seconds', 0):.0f}s", flush=True)
        del res, o
        gc.collect()
    return runs_table()


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
        try:
            path.write_bytes(pickle.dumps(st))
        except OSError:
            pass
    return st[key]


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run: ({(kind, wiring): outputs}, runs table). Kinds: 'family', 'E110_A', and each arm (its step-1 runs
    here, U and SH0 on S from E122's step 1)."""
    from experiments import E110_no_kckc as E110
    outs, rows = {}, []
    ra = pd.concat([pd.read_csv(p) for p in sorted((H.shared_dir() / 'results' / 'E110_no_kckc').glob('runs_A_*.csv'))],
                   ignore_index=True).drop_duplicates('tag', keep='last').set_index('tag')

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E122.type_rates(o))
        r.update(bands(r))
        rows.append(r)

    for w in WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], source='saved fit', **F['p'])
        a = E110.load_a(w)
        if a is not None:
            add('E110_A', w, a, tag=f'A_{w}', folder='E110_no_kckc/A', source='E110 part A (read)',
                **{k: float(ra.loc[f'A_{w}', k]) for k in PARAM_KEYS})
    r1 = pd.read_csv(H.shared_dir() / 'results' / E122.EXPERIMENT / 'runs_1.csv').set_index('tag')
    for arm, e in E122_READ.items():
        for w in S_WIRINGS:
            o = E122.load_run('step1', f'S1_{e}_{w}')
            if o is not None:
                add(arm, w, o, tag=f'S1_{e}_{w}', folder=f'{E122.EXPERIMENT}/step1', source='E122 step 1 (read)',
                    apl_borne=False, **{k: float(r1.loc[f'S1_{e}_{w}', k]) for k in PARAM_KEYS})
    for r in runs_table().drop_duplicates('tag', keep='last').to_dict('records'):
        o = load_run(r['tag'])
        if o is not None:
            add(r['arm'], r['wiring'], o, tag=r['tag'], folder=f'{EXPERIMENT}/step1', source=f"E123 ({r['host']})",
                apl_borne=bool(r['apl_borne']), thr_factor=r['thr_factor'], apl_factor=r['apl_factor'],
                at_bound_thr=r['at_bound_thr'], **{k: float(r[k]) for k in PARAM_KEYS})
    return outs, pd.DataFrame(rows)


RUN_COLS = E122.RUN_COLS


def readings(outs, runs) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E110's
    `scores_of`, the runs' means of RUN_COLS ('run_<col>'), the number of runs outside the bands."""
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
                                 **{f'run_{c}': float(s[c].astype(float).mean()) for c in RUN_COLS if c in s.columns},
                                 **E122.scores_of(d)))
    return pd.DataFrame(rows)


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    return E122._v(rd, kind, unit, col, readout)


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """The class of an arm (the module docstring); incomplete unless it has all eight wirings."""
    if not all(((runs['kind'] == arm) & (runs['wiring'] == w)).any() for w in WIRINGS):
        return dict(arm=arm, complete=False)
    d, keeps, eq, eq2 = {}, True, True, True
    for t in TARGETS:
        col = f'{t}_pattern'
        pops = {}
        for u in POPS:
            a, f, z = _v(rd, arm, u, col), _v(rd, 'family', u, col), _v(rd, 'E110_A', u, col)
            pops[u] = dict(arm=a, family=f, E110_A=z, half_way=(f + z) / 2, d_family=a - f)
            keeps &= a <= (f + z) / 2
            eq &= a - f <= WITHIN
            eq2 &= abs(a - f) <= WITHIN
        per_w = {w: _v(rd, arm, w, col) - _v(rd, 'family', w, col) for w in WIRINGS}
        eq &= sum(v <= WITHIN for v in per_w.values()) >= N_WIRINGS_NEEDED
        eq2 &= sum(abs(v) <= WITHIN for v in per_w.values()) >= N_WIRINGS_NEEDED
        d[t] = dict(populations=pops, d_family_per_wiring=per_w)
    sub = runs[runs['kind'] == arm]
    in_bands = bool(sub['in_bands'].all())
    flagged = sorted(sub.loc[~sub['in_bands'].astype(bool), 'wiring'])
    equal = bool(eq and in_bands)
    cls = 'equal' if equal else ('keeps the gain' if keeps else 'neither')
    return dict(arm=arm, complete=True, cls=cls, keeps_gain=bool(keeps), equal=equal,
                equal_pattern_only=bool(eq), equal_two_sided=bool(eq2 and in_bands), all_in_bands=in_bands,
                runs_outside_bands=flagged, targets=d,
                rate_rho_on_mean8=float(sub['rate_rho_on'].mean()),
                mean_effect_raw_all8=_v(rd, arm, 'all8', 'mean_effect', 'raw'),
                loss_raw_all8=_v(rd, arm, 'all8', 'all_loss', 'raw'))


def draws_verdict(classes: Mapping[str, dict], runs: pd.DataFrame) -> dict:
    """Whether the four seeds of the shuffle agree, and whether the shuffle costs response rate ρ with APL on."""
    if not all(classes.get(s, {}).get('complete') for s in DRAWS):
        return dict(complete=False)
    fam = float(runs.loc[runs['kind'] == 'family', 'rate_rho_on'].mean())
    per = {s: float(runs.loc[runs['kind'] == s, 'rate_rho_on'].mean()) for s in DRAWS}
    cls = {s: classes[s]['cls'] for s in DRAWS}
    mean = float(np.mean(list(per.values())))
    below = sum(v < fam for v in per.values())
    return dict(complete=True, classes=cls, agree=len(set(cls.values())) == 1, family_rate_rho_on=fam,
                rate_rho_on_per_seed=per, rate_rho_on_mean=mean, seeds_below_family=int(below),
                costs_rate_rho=bool(mean < fam - RHO_COST and below >= RHO_SEEDS_NEEDED))


def step2_plan(classes: Mapping[str, dict]) -> dict:
    """The step-2 gate: the simplest random arm classed 'equal' (smallest k, then AK) and HBab if 'equal'."""
    def ok(a):
        return classes.get(a, {}).get('cls') == 'equal'
    random_arm = next((a for a in [f'K{k}' for k in LADDER] + ['AK'] if ok(a)), None)
    sets = [a for a in (random_arm, 'HBab' if ok('HBab') else None) if a]
    complete = all(classes.get(a, {}).get('complete') for a in [f'K{k}' for k in LADDER] + ['AK', 'HBab'])
    return dict(complete=complete, refit_sets=sets, random_arm=random_arm, apl_borne=ok('HBab'),
                note=('each set: refits on all eight wirings (threshold, APL->KC with KC->APL and bouton->APL tied, '
                      'APL->bouton cap 3, KC-KC scale free; 150 evaluations, from the step-1 held point), and S\'s mean '
                      'parameters carried to T as in E122. Not implemented yet: the APL-borne refit needs a KC-KC scale '
                      'for the APL-off condition in E10_joint_fit.py, checked unchanged when off.'))


STATS = (('pattern part, readout', 'all_pattern', 'cal'), ('amplitude part, readout', 'all_amplitude', 'cal'),
         ('shape part, readout', 'all_shape', 'cal'), ('spread ratio, readout', 'all_sd_ratio', 'cal'),
         ('Pearson r on the per-pair APL effect, readout', 'all_pearson_effect', 'cal'),
         ('ρ on the per-pair APL effect, readout', 'all_rho_effect', 'cal'),
         ('mean APL effect, readout', 'mean_effect', 'cal'), ('pattern part, half A / B / C / D', None, 'cal'),
         ('APL-effect loss, raw', 'all_loss', 'raw'), ('mean APL effect, raw', 'mean_effect', 'raw'),
         ('pattern part, raw', 'all_pattern', 'raw'),
         ('ρ on odor-odor distances, APL on', 'rho_on', 'raw'), ('ρ on odor-odor distances, APL off', 'rho_off', 'raw'),
         ('response rate, APL on', 'run_rate_on', 'raw'), ('response rate, APL off', 'run_rate_off', 'raw'),
         ('APL-off / APL-on response rate', 'run_ratio', 'raw'), ('response rate ρ, APL on', 'run_rate_rho_on', 'raw'),
         ('response rate ρ, APL off', 'run_rate_rho_off', 'raw'), ('odors per KC (TV), APL on', 'run_opk_on_tv', 'raw'),
         ('odors per KC (TV), APL off', 'run_opk_off_tv', 'raw'), ('γ response rate, APL on', 'run_rate_on_g', 'raw'),
         ('γ off / on ratio', 'run_ratio_g', 'raw'), ('largest spike count, APL on', 'run_max_count_on', 'raw'),
         ('largest spike count, APL off', 'run_max_count_off', 'raw'), ('somatic threshold', 'run_fixed_thr', 'raw'),
         ('KC-KC scale (APL on)', 'run_wKCKC_scale', 'raw'), ('runs outside the bands', 'n_outside_bands', 'raw'))


def comparison_table(rd: pd.DataFrame, kinds: Sequence[str], unit: str) -> str:
    ks = [k for k in kinds if ((rd['kind'] == k) & (rd['unit'] == unit)).any()]
    lines = [f'| statistic ({unit}) | ' + ' | '.join(ks) + ' |', '|---|' + '---:|' * len(ks)]
    for lab, col, ro in STATS:
        if col is None:
            cells = [' / '.join(f'{_v(rd, k, unit, f"{h}_pattern", ro):.3f}' for h in E108.HALVES) for k in ks]
        else:
            cells = [f'{_v(rd, k, unit, col, ro):.3f}' for k in ks]
        lines.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def runs_bands_table(runs: pd.DataFrame) -> str:
    cols = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'opk_on_tv', 'opk_off_tv', 'max_count_on', 'max_count_off',
            'fixed_thr')
    lines = ['| kind | wiring | ' + ' | '.join(cols) + ' | in bands | flag |', '|---|---|' + '---:|' * len(cols) + '---|---|']
    for r in runs.itertuples():
        lines.append(f'| {r.kind} | {r.wiring} | ' + ' | '.join(f'{float(getattr(r, c)):.3f}' for c in cols)
                     + f" | {'yes' if r.in_bands else 'no'} | {r.realism_flag if isinstance(r.realism_flag, str) else ''} |")
    return '\n'.join(lines)


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    arms = [a for a in ARMS if a in set(runs['kind'])]
    classes = {a: classify(rd, runs, a) for a in arms}
    dv = draws_verdict(classes, runs)
    plan = step2_plan(classes)
    ref = E122.hemibrain_kckc(E122.model_kc_index().get_level_values(H.KC_ID_LEVEL)).to_numpy(dtype=float)
    mstats = {b: matrix_stats(write_matrix(b)[1].to_numpy(dtype=float), ref)
              for b in dict.fromkeys(['HB'] + [base_of(a) for a in ARMS])}
    rt = runs_table()
    shas = {a: sorted(set(rt.loc[rt['arm'] == a, 'matrix_sha256'])) for a in set(rt.get('arm', []))}
    bad = {a: s for a, s in shas.items() if len(s) != 1 or s[0] != write_matrix(a)[2]}
    if bad:
        raise RuntimeError(f"E123 score: an arm's matrix differs between hosts or from this host's: {bad}")
    v = dict(generated=_now(), host=_host(), classes=classes, draws=dv, step2=plan, matrices=mstats, matrix_sha256=shas,
             rules=dict(within=WITHIN, wirings_needed=N_WIRINGS_NEEDED, rho_cost=RHO_COST,
                        rho_seeds_needed=RHO_SEEDS_NEEDED, targets=list(TARGETS), populations=list(POPS)),
             n_runs=int(len(runs)), complete=all(c.get('complete') for c in classes.values()) and len(arms) == len(ARMS))
    E122._write_json(out / 'verdict.json', v)
    E122._write_json(out / 'step2_plan.json', dict(plan, generated=_now()))
    kinds = ['family', 'E110_A'] + list(ARMS)
    md = ['# E123 comparison tables (generated)', '',
          'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
          f"Draws: {json.dumps({k: dv.get(k) for k in ('agree', 'costs_rate_rho', 'rate_rho_on_mean', 'family_rate_rho_on')})}",
          '', f"Step 2: {plan['refit_sets']}", '']
    for u in ('S4', 'T4', 'all8'):
        md += [f'## {u}', '', comparison_table(rd, kinds, u), '']
    md += ['## The matrices', '', '| matrix | ' + ' | '.join(next(iter(mstats.values())).keys()) + ' |',
           '|---|' + '---:|' * len(next(iter(mstats.values())))]
    md += [f'| {b} | ' + ' | '.join(f'{x:.3f}' for x in s.values()) + ' |' for b, s in mstats.items()]
    md += ['', '## Runs and bands', '', runs_bands_table(runs), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    E122._write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(),
                                                         script='experiments/E123_kckc_form.py', commits=H.code_commits()))
    print(json.dumps({a: c.get('cls') for a, c in classes.items()}), flush=True)
    return v


def part_2() -> dict:
    """Step 2's gate (step2_plan.json from `score`); the refits themselves are not implemented yet."""
    plan = E122._read_json(out_dir() / 'step2_plan.json')
    print(json.dumps(plan, indent=1), flush=True)
    if plan.get('refit_sets'):
        raise NotImplementedError(f"E123 step 2 refits of {plan['refit_sets']}: not implemented yet ({plan['note']})")
    return plan


def write_variants() -> pd.DataFrame:
    """variants.csv (every E123 run, absolute parameters) and model_pair_distances_long.csv (raw spike counts)."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    for r in runs[runs['folder'] == f'{EXPERIMENT}/step1'].itertuples():
        o = outs[(r.kind, r.wiring)]
        rows.append(dict(tag=r.tag, variant=r.tag, part='S1', arm=r.kind, arm_label=ARM_LABELS[r.kind], base='bouton', k=2.0,
                         wiring=r.wiring, apl_borne=bool(r.apl_borne), kckc_scale_off=0.0 if r.apl_borne else r.wKCKC_scale,
                         kind=(f"C2's family's {r.wiring}: KC-KC matrix of arm {r.kind} ({ARM_LABELS[r.kind]}), at the "
                               "fit's scales, both response rates held"),
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         kckc_matrix=f'matrices/kckc_{base_of(r.kind)}.parquet', fixed_thr=r.fixed_thr,
                         wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale,
                         wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale,
                         rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag, subdir='step1'))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def all_done() -> bool:
    have = set(runs_table()['tag'])
    return all(f'S1_{a}_{w}' in have for a, w in planned())


def chain(wirings: Sequence[str], arms: Sequence[str] = ARMS, wait_others_s: float = 6 * 3600) -> None:
    """check, step 1 on `wirings`, then (waiting up to `wait_others_s` for the other lanes) score and variants."""
    part_check()
    part_1(arms, wirings)
    t0 = time.time()
    while not all_done():
        if time.time() - t0 > wait_others_s:
            print('E123: other lanes not done; scoring what there is', flush=True)
            break
        print(f'E123: waiting for the other lanes ({_now()})', flush=True)
        time.sleep(300)
    score()
    write_variants()
    print(f'E123 chain done {_now()}', flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', '1', 'score', '2', 'variants', 'chain', 'plan'))
    ap.add_argument('--wirings', nargs='+', default=list(WIRINGS), choices=list(WIRINGS))
    ap.add_argument('--arms', nargs='+', default=list(ARMS), choices=list(ARMS))
    a = ap.parse_args()
    if a.part == 'check':
        part_check()
    elif a.part == '1':
        part_1(a.arms, a.wirings)
    elif a.part == 'score':
        score()
    elif a.part == '2':
        part_2()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'chain':
        chain(a.wirings, a.arms)
    elif a.part == 'plan':
        for arm, w in planned(a.arms, a.wirings):
            print(arm, w)


if __name__ == '__main__':
    main()
