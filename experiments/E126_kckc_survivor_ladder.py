#!/usr/bin/env python3
"""E126: the ladder of the pooled survivor: one neuron per KC type along the dose and KC groups without type.

The entry is item E126 of experiment_priorities.md (proposed by viki 2026-11-21 18:13 from E125's next steps, moved into
the list 2026-11-22 09:59 under the idle-host rule; planned by Fable subagents, claude-fable-5-2). It follows E125
(experiments/E125_kckc_survivor.py), whose runner (two simulators, one per condition's KC-KC matrix, both response rates
re-held by E77_input_screen.hold_rates), check and classes it reuses by import. "The family" is C2's breadth-2
APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit; S is real, P1, P2, M1 and T is P3, P4, M2, M3.

The design (fixed scales, nothing fitted; E125's). The KC-KC scale is the wiring fit's, s, in both conditions; with APL
on the matrix is (1 - d) HB + d X, with APL off d X, HB the hemibrain matrix (mean row sum 1) and X a survivor with HB's
total weight. So the hemibrain matrix at (1 - d) of the fitted scale is present with APL on only and the survivor at d of
its total weight in both.
- Step 1, the ladder (`--part step1`): X one neuron per KC type (E122's TY, as E125's TY50), at d 0.75 and 0.25: TY75,
  TY25 (16 runs; tags S1_<arm>_<wiring>). TY75 runs first, since step 2's gate is on it.
- Step 2, groups without type (`--part step2`, gated): X one neuron per group, the KCs dealt at random into groups of
  the KC types' sizes (`random_groups`: the model KCs' type labels permuted with numpy's default_rng(seed)), each
  group's block entries its members' total incoming hemibrain weight over its size squared (`group_survivor`, which
  with the groups the types is E122's TY exactly), two draws (seeds 0, 1) at d 0.75: RG75s0, RG75s1 (16 runs; tags
  S2_<arm>_<wiring>). Refused unless verdict.json (from `--part score` after step 1) has the gate open: TY75 keeps the
  gain by E125's class, is no corner, and has no threshold under 60 Hz on any wiring (`gate`).
- References, read and not run again (E125's `collect`): the family, KC-KC at zero (E110's part A; part B beside),
  E124's ladder (HBoff0, HBoff25, HBoff50, HBoff75), E123's uniform neuron U1, and E125's arms (PO25, PO50, PO75, DG25,
  DG50, DG75, TY50).

Parts (CLI `--part`):
- check: (i) the hemibrain matrix is E123's (sha256 against E123's HBab runs); (ii) with d 0 (Z0: on HB, off zeros) the
  two-simulator runner reproduces E123's held point S1_HBab_real exactly (spike counts; factors to 1e-12 relative);
  (iii) the survivor matrices: E126's TY at d 0.5 is E125's TY50 (sha256 of both conditions' matrices against E125's
  runs), `group_survivor` on the types is E122's TY exactly, and each draw's groups have the types' sizes. check.json;
  failing stops everything.
- step1, step2 (simulations, resumable, per lane; `--wirings`, `--arms`): rows in runs_<step>_<host>_<lane>.csv (a tag
  in any host's runs_1_*.csv / runs_2_*.csv is skipped), held points in starts_<step>_<host>_<lane>.json, provenance
  (`H.model_provenance` per point) in provenance_<step>_<host>_<lane>.json; the runs saved locally
  (outputs/E126_kckc_survivor_ladder/step<step>/) and in the model-output archive (data/model_outputs/<EXPERIMENT>/).
- score (no simulation): readings.csv, runs.csv, verdict.json (classes, gate, "differ" comparisons, gain kept,
  outcomes), comparison_table.md, provenance.json (`H.model_provenance`) and provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- plan: the points a lane would run.

The classes (fixed before any run; E125's, unchanged; calibrated pattern part, S's and T's populations apart, against
all recordings and the four halves A to D):
- "keeps the gain" (E123.classify) and "does what the matrix does" (keeps the gain at d 0.5 or 0.75 and not more than
  0.003 above E124's rung at the same d on both populations, every target), as in E125; the arm at 0.25 can only keep
  the gain.
- Two forms "differ" at a d (`differ`) only if one is lower on at least DIFFER_WIRINGS (7) of the 8 wirings against all
  recordings, and lower on both populations against every one of the five targets; otherwise "not told apart".
- A run with the threshold or APL at a bound of its bisection is a corner; an arm with corners on more than two wirings
  gets no class ('corner (no class)'), and a comparison with it none either.
- Beside every class, the gain kept (E124.gain_share: (E110 A - arm) / (E110 A - family), all recordings) on all eight
  wirings ('all8'), on S and T, and per wiring.
- The outcomes (flags for the write-up, which reads them): (a) TY75 does what the matrix does with no threshold under
  60 Hz; (b) TY75 fails to keep the gain, or needs a threshold under 60 Hz; 'between_a_b' TY75 keeps the gain with no
  threshold under 60 Hz but is more than 0.003 above E124's rung; with step 2: 'types_needed' (a) and a draw fails to
  keep the gain or TY75 is lower than a draw by the "differ" rule; (c) TY75 keeps the gain and both draws keep it,
  neither with TY75 lower; (d) an arm of this experiment is a corner.
- No correlation chooses anything, no pair is weighted, and no threshold moves after a run is read.
Realism: every run's four bands, γ's response rate and off / on ratio, the largest spike counts, the threshold and the
APL->KC scale; a threshold under 60 Hz is recorded as a realism cost even inside the bands (`thr_below_60`).

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it, and the model loader's hemibrain KC->KC counts and KC types; no other connectome table; no olfsysm change. Nothing is
published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E126_kckc_survivor_ladder.py --part check
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E126_kckc_survivor_ladder.py --part step1 --wirings P3 P4
    ./venv/bin/python -u experiments/E126_kckc_survivor_ladder.py --part score
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
from experiments import E123_kckc_form as E123
from experiments import E124_kckc_apl_borne as E124
from experiments import E125_kckc_survivor as E125

EXPERIMENT = 'E126_kckc_survivor_ladder'
SCRIPT = 'experiments/E126_kckc_survivor_ladder.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E125.S_WIRINGS, E125.T_WIRINGS, E125.WIRINGS
TYPES = E125.TYPES
STEP1_ARMS = ('TY75', 'TY25')
GROUP_SEEDS = (0, 1)
STEP2_ARMS = tuple(f'RG75s{s}' for s in GROUP_SEEDS)
STEPS = {1: STEP1_ARMS, 2: STEP2_ARMS}
ARMS = STEP1_ARMS + STEP2_ARMS
FORMS = dict(E125.FORMS, RG='one neuron per group of KCs dealt at random into groups of the KC types\' sizes')
E125_ARMS = E125.ARMS                                               # read (PO25..PO75, DG25..DG75, TY50)
LADDER = E125.LADDER                                                # E124's rungs (read)
E123_U = E125.E123_U
WITHIN = E125.WITHIN
MATRIX_LIKE_DOSES = E125.MATRIX_LIKE_DOSES
CORNER_WIRINGS_MAX = E125.CORNER_WIRINGS_MAX
REALISM_THR_HZ = E125.REALISM_THR_HZ
DIFFER_WIRINGS = 7
TARGETS, POPS, READOUTS = E125.TARGETS, E125.POPS, E125.READOUTS
PARAM_KEYS = E125.PARAM_KEYS
OFF_MARK = E125.OFF_MARK
# the "differ" comparisons, (a, b) at the same d
COMPARISONS = (('TY25', 'PO25'), ('TY25', 'HBoff25'), ('TY25', 'DG25'),
               ('TY50', 'PO50'), ('TY50', 'HBoff50'), ('TY50', 'DG50'),
               ('TY75', 'PO75'), ('TY75', 'HBoff75'), ('TY75', 'DG75')) + tuple(
    (r, b) for r in STEP2_ARMS for b in ('TY75', 'PO75', 'HBoff75'))
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths

def out_dir() -> Path:
    """<shared>/results/E126_kckc_survivor_ladder/ (created)."""
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
# the arms and their matrices

def form_of(arm: str) -> str:
    """'TY', 'RG', E125's forms 'PO', 'DG' (or 'Z' for Z0)."""
    return 'Z' if arm == 'Z0' else arm[:2]


def dose_of(arm: str) -> float:
    """d of an arm (TY25 -> 0.25, RG75s0 -> 0.75; Z0 -> 0)."""
    return 0.0 if arm == 'Z0' else int(arm[2:4]) / 100.0


def seed_of(arm: str) -> Optional[int]:
    """The draw's seed of a random-group arm (RG75s1 -> 1), else None."""
    return int(arm.split('s', 1)[1]) if form_of(arm) == 'RG' else None


def step_of(arm: str) -> int:
    return next(s for s, a in STEPS.items() if arm in a)


def arm_label(arm: str) -> str:
    d = dose_of(arm)
    if arm == 'Z0':
        return E125.ARM_LABELS['Z0']
    s = f', draw seed {seed_of(arm)}' if form_of(arm) == 'RG' else ''
    return (f"the hemibrain matrix at {1 - d:g} of the fit's scale with APL on only, and {FORMS[form_of(arm)]}{s} at "
            f"{d:g} of its total weight in both conditions")


def random_groups(types: np.ndarray, seed: int) -> np.ndarray:
    """Group labels: the KCs dealt at random into groups of the KC types' sizes (the type labels permuted by
    numpy's default_rng(seed)); a KC's group is named by the type whose size it has."""
    types = np.asarray(types, dtype=object)
    return types[np.random.default_rng(seed).permutation(len(types))]


def group_survivor(ref: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """One neuron per group (the groups named by TYPES; other labels get none): within each group every entry its
    members' total incoming weight of `ref` (rows are post) over its size squared. With labels = the types this is
    E122's TY exactly; the total is `ref`'s when every KC is in a group."""
    n = len(labels)
    out = np.zeros((n, n))
    for t in TYPES:
        m = labels == t
        if m.any():
            out[np.ix_(m, m)] = ref[m].sum() / m.sum() / m.sum()
    return out


def condition_arrays(arm: str, ref: np.ndarray, types: np.ndarray) -> Dict[str, np.ndarray]:
    """{'on': (1 - d) ref + d X, 'off': d X} of an arm: E125's for TY, PO, DG and Z0; for RG, X the random groups'."""
    if form_of(arm) != 'RG':
        return E125.condition_arrays(arm, ref, types)
    d = dose_of(arm)
    X = group_survivor(ref, random_groups(types, seed_of(arm)))
    return {'on': (1.0 - d) * ref + d * X, 'off': d * X}


def group_info(arm: str, types: np.ndarray) -> dict:
    """A random-group arm's draw: each group's size and type composition, and the share of KCs in the group named by
    their own type."""
    if form_of(arm) != 'RG':
        return {}
    g = random_groups(types, seed_of(arm))
    comp = {t: {u: int(((g == t) & (types == u)).sum()) for u in TYPES} for t in TYPES}
    return dict(group_seed=seed_of(arm), group_sizes=json.dumps({t: int((g == t).sum()) for t in TYPES}),
                group_composition=json.dumps(comp), group_same_as_type=float((g == types).mean()))


def matrix_path(arm: str, cond: str) -> Path:
    return local_dir('matrices') / f'kckc_{arm}_{cond}.parquet'


def _model_ref() -> Tuple[pd.DataFrame, np.ndarray]:
    if 'ref' not in _CACHE:
        idx = E122.model_kc_index()
        _CACHE['ref'] = (E122.hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL)), E122.kc_types_of(idx))
    return _CACHE['ref']


def write_matrices(arm: str) -> Dict[str, Tuple[Path, pd.DataFrame, str]]:
    """{'on' | 'off': (path, matrix, sha256 of its values)} for the model's KCs, written to
    outputs/.../matrices/kckc_<arm>_<cond>.parquet unless there with the same values (columns as strings)."""
    from experiments.E77_input_screen import frame_sha256
    if arm in _CACHE.get('mat', {}):
        return _CACHE['mat'][arm]
    ref, types = _model_ref()
    arrs = condition_arrays(arm, ref.to_numpy(dtype=float), types)
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


def held(arm: str, w: str, params: Mapping[str, float], label: str) -> dict:
    """One held point: E77_input_screen.hold_rates on E125's two-condition simulator of the arm's matrices (the APL-off
    condition marked by off_extra; the KC-KC scale the fit's in both)."""
    from experiments.E77_input_screen import hold_rates
    P = E83.wiring_parts(w)
    m = write_matrices(arm)
    res = hold_rates(lambda: E125.two_condition_simulator(w, params, m['on'][1], m['off'][1]),
                     target_off=P['target_off'], target_on=P['target_on'], params=dict(params), label=label,
                     off_extra={OFF_MARK: 1.0})
    res['params_off'] = {k: v for k, v in res['params_off'].items() if k != OFF_MARK}
    return res


# ---------------------------------------------------------------------------
# runs

def _tables(prefix: str) -> List[Path]:
    return sorted(out_dir().glob(f'{prefix}_*.csv'))


def runs_table() -> pd.DataFrame:
    """Every row of every host and lane (runs_1_*.csv and runs_2_*.csv)."""
    fs = [pd.read_csv(p) for s in STEPS for p in _tables(f'runs_{s}')]
    return pd.concat(fs, ignore_index=True, sort=False) if fs else pd.DataFrame(columns=['tag'])


def tag_of(arm: str, w: str) -> str:
    return f'S{step_of(arm)}_{arm}_{w}'


def save_run(tag: str, o: Dict[str, tuple], step: int) -> None:
    """The run's matrices, locally and in the model-output archive (which another host's score part reads)."""
    for d in (local_dir(f'step{step}'), archive_dir() / f'step{step}'):
        d.mkdir(parents=True, exist_ok=True)
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1])):
                p = d / f'{tag}_{cc}_{k}.parquet'
                tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
                x.to_parquet(tmp)
                os.replace(tmp, p)


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    step = int(tag[1])
    for d in (local_dir(f'step{step}'), archive_dir() / f'step{step}'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_row(tag: str, arm: str, w: str, res: dict, o: Dict[str, tuple], cc: dict, lane: str) -> dict:
    P = E83.wiring_parts(w)
    m = write_matrices(arm)
    row = dict(step=step_of(arm), tag=tag, arm=arm, form=form_of(arm), dose=dose_of(arm), wiring=w,
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
    row.update(group_info(arm, _model_ref()[1]))
    row.update(E83.realism(o))
    row.update(E122.type_rates(o))
    row.update(E123.bands(row))
    row['thr_below_60'] = bool(float(row['fixed_thr']) < REALISM_THR_HZ)
    row.update(cc)
    return row


def planned(arms: Sequence[str] = ARMS, wirings: Sequence[str] = WIRINGS) -> List[Tuple[str, str]]:
    """The (arm, wiring) points run here: arms in ARMS order (TY75 first), wirings in WIRINGS order."""
    return [(a, w) for a in ARMS if a in arms for w in WIRINGS if w in wirings]


def all_done(step: int) -> bool:
    have = set(runs_table()['tag'])
    return all(tag_of(a, w) in have for a, w in planned(STEPS[step]))


# ---------------------------------------------------------------------------
# check

def part_check() -> dict:
    """(i) the hemibrain matrix is E123's; (ii) Z0 (d 0) on 'real' reproduces E123's S1_HBab_real exactly; (iii) the
    survivors: TY at d 0.5 is E125's TY50, `group_survivor` on the types is E122's TY, the draws keep the types' sizes.
    check.json (a failed one kept beside as check_failed_<time>.json); `passes` all three."""
    from experiments.E77_input_screen import frame_sha256
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
    # (iii) the survivors
    ref, types = _model_ref()
    r = ref.to_numpy(dtype=float)
    r125 = E125.runs_table()
    ty50 = condition_arrays('TY50', r, types)
    sha_ty50 = {c: frame_sha256(pd.DataFrame(ty50[c], index=ref.index, columns=ref.columns)) for c in ('on', 'off')}
    e125_ty50 = {c: sorted(set(r125.loc[r125['arm'] == 'TY50', f'matrix_{c}_sha256'])) for c in ('on', 'off')}
    sizes = {t: int((types == t).sum()) for t in TYPES}
    draws = {}
    for a in STEP2_ARMS:
        g = random_groups(types, seed_of(a))
        X = condition_arrays(a, r, types)['off'] / dose_of(a)
        draws[a] = dict(sizes={t: int((g == t).sum()) for t in TYPES}, same_as_type=float((g == types).mean()),
                        total=float(X.sum()), ref_total=float(r.sum()))
    surv = dict(ty50_same_as_E125={c: e125_ty50[c] == [sha_ty50[c]] for c in ('on', 'off')},
                ty50_sha256=sha_ty50, e125_ty50_sha256=e125_ty50,
                groups_on_types_is_TY=bool(np.array_equal(group_survivor(r, types), E122.arm_array('TY', r, types))),
                type_sizes=sizes, draws=draws,
                draws_keep_sizes=all(d['sizes'] == sizes for d in draws.values()),
                draws_keep_total=all(abs(d['total'] - d['ref_total']) <= 1e-9 * d['ref_total'] for d in draws.values()),
                all_kcs_typed=bool(np.isin(types, TYPES).all()))
    surv_ok = bool(all(surv['ty50_same_as_E125'].values()) and surv['groups_on_types_is_TY'] and surv['draws_keep_sizes']
                   and surv['draws_keep_total'])
    cc = E108._commits_clean()
    P = E83.wiring_parts('real')
    t0 = time.time()
    hr = held('Z0', 'real', P['p'], 'E126 check Z0 real')
    saved = E123.load_run('S1_HBab_real')
    if saved is None:
        raise RuntimeError("E126 check: E123's S1_HBab_real is neither local nor in the archive")
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
    res = dict(checked=_now(), host=_host(), matrices=mats, held_point=held_point, survivors=surv, commits=cc)
    res['passes'] = bool(exact and mats['same'] and surv_ok)
    E122._write_json(path, res)
    print(f"E126 check: passes {res['passes']}; matrix {mats['same']}; survivors {surv_ok}; held point {held_point}",
          flush=True)
    if not res['passes']:
        raise RuntimeError('E126 check failed (check.json)')
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
        print(f'E126: waiting for check.json ({_now()})', flush=True)
        time.sleep(every)


# ---------------------------------------------------------------------------
# steps 1 and 2

def step2_allowed() -> Tuple[bool, str]:
    """Whether step 2 may run: verdict.json (written by score after every step-1 run) has the gate open, and it was
    scored on every step-1 run that exists now."""
    v = E122._read_json(out_dir() / 'verdict.json')
    if not v:
        return False, 'no verdict.json: run --part score after step 1'
    g = v.get('gate', {})
    have = sorted(t for t in runs_table()['tag'] if str(t).startswith('S1_'))
    if sorted(v.get('step1_tags', [])) != have:
        return False, 'verdict.json was not scored on every step-1 run there is: run --part score again'
    if not g.get('complete'):
        return False, f"step 1 not complete ({g.get('why')})"
    if not g.get('open'):
        return False, f"the gate is closed: {g.get('why')}"
    return True, g.get('why', '')


def part_runs(step: int, arms: Optional[Sequence[str]] = None, wirings: Sequence[str] = WIRINGS,
              lane: Optional[str] = None) -> pd.DataFrame:
    """Step 1 (TY75, TY25) or step 2 (RG75s0, RG75s1; refused unless `step2_allowed`), resumable, one process per
    lane."""
    wait_for_check()
    arms = [a for a in (arms or STEPS[step]) if a in STEPS[step]]
    if step == 2:
        ok, why = step2_allowed()
        if not ok:
            raise RuntimeError(f'E126 step 2 refused: {why}')
        print(f'E126 step 2: gate open ({why})', flush=True)
    lane = lane or f'step{step}-{lane_name(wirings)}'
    table = out_dir() / f'runs_{step}_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_{step}_{_host()}_{lane}.json'
    starts = out_dir() / f'starts_{step}_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in planned(arms, wirings):
        tag = tag_of(arm, w)
        if tag in set(runs_table()['tag']):
            print(f'E126 {tag}: done', flush=True)
            continue
        E108._still_clean(cc)
        P = E83.wiring_parts(w)
        res = held(arm, w, P['p'], f'E126 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o, step)
        row = run_row(tag, arm, w, res, o, cc, lane)
        E83._append_row(table, row)
        prov = E122._read_json(prov_path) or dict(part=f'step{step}', script=SCRIPT, points={})
        m = write_matrices(arm)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
            change=(f"arm {arm} ({arm_label(arm)}): KC-KC matrix (1 - d) HB + d X with APL on and d X with APL off "
                    f"(d {dose_of(arm):g}, X {FORMS[form_of(arm)]}), each its own simulator, at the fit's KC-KC scale; "
                    "fixed_thr re-held with APL off to the fit's APL-off count, then the tied APL scales to its APL-on "
                    "count (E77_input_screen.hold_rates)"),
            matrices={c: dict(path=str(m[c][0].relative_to(H.HERE)), sha256=m[c][2]) for c in ('on', 'off')},
            group_seed=seed_of(arm), fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws')))
        E122._write_json(prov_path, prov)
        st = E122._read_json(starts)
        st[tag] = dict({k: float(res['params'][k]) for k in PARAM_KEYS}, source=f'held point {tag}')
        E122._write_json(starts, st)
        print(f"E126 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/{res['target_off']} off); thr {row['fixed_thr']:.1f}; bands {row['in_bands']}; flag "
              f"'{row.get('realism_flag') or ''}'; {res.get('seconds', 0):.0f}s", flush=True)
        del res, o
        gc.collect()
    return runs_table()


# ---------------------------------------------------------------------------
# score

def _dist(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83.pair_distances), cached locally by the sha256 of its spike counts and the readout;
    the cache starts from E125's (read only), which holds the references'."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    if 'dist' not in _CACHE:
        st = {}
        for p in (E125.local_dir() / 'distances_cache.pkl', local_dir() / 'distances_cache.pkl'):
            try:
                if p.exists():
                    st.update(pickle.loads(p.read_bytes()))
            except Exception:
                pass
        _CACHE['dist'] = st
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
    """Every run: ({(kind, wiring): outputs}, runs table). E125's `collect` (the family, E110's A and B, E124's ladder
    with HBoff0, E123's U1, E125's arms), then this experiment's arms."""
    outs, r125 = E125.collect()
    rows = r125.to_dict('records')
    for r in runs_table().drop_duplicates('tag', keep='last').to_dict('records'):
        o = load_run(r['tag'])
        if o is None:
            continue
        outs[(r['arm'], r['wiring'])] = o
        row = dict(kind=r['arm'], wiring=r['wiring'], tag=r['tag'], folder=f"{EXPERIMENT}/step{int(r['step'])}",
                   source=f"E126 ({r['host']})", form=r['form'], dose=r['dose'], thr_factor=r['thr_factor'],
                   apl_factor=r['apl_factor'], at_bound_thr=bool(r['at_bound_thr']),
                   at_bound_apl=bool(r['at_bound_apl']), **{k: float(r[k]) for k in PARAM_KEYS},
                   **E83.realism(o), **E122.type_rates(o))
        row.update(E123.bands(row))
        rows.append(row)
    runs = pd.DataFrame(rows)
    runs['thr_below_60'] = runs['fixed_thr'].astype(float) < REALISM_THR_HZ
    return outs, runs


def readings(outs, runs) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E122's
    `scores_of`, the runs' means of E122.RUN_COLS ('run_<col>'), the number of runs outside the bands (E125's)."""
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


def complete(runs: pd.DataFrame, kind: str) -> bool:
    return all(((runs['kind'] == kind) & (runs['wiring'] == w)).any() for w in WIRINGS)


def matrix_like(rd: pd.DataFrame, arm: str, keeps_gain: bool) -> dict:
    """E125's "does what the matrix does" (d 0.5 or 0.75, keeps the gain, not more than WITHIN above E124's rung at the
    same d on both populations, every target), with this experiment's arm names."""
    d = dose_of(arm)
    if d not in MATRIX_LIKE_DOSES:
        return dict(matrix_like=None, matrix_like_all_only=None, rung=None)
    rung = LADDER[d]
    per = {t: {u: _v(rd, arm, u, f'{t}_pattern') - _v(rd, rung, u, f'{t}_pattern') for u in POPS} for t in TARGETS}
    within = all(v <= WITHIN for t in TARGETS for v in per[t].values())
    within_all = all(v <= WITHIN for v in per['all'].values())
    return dict(matrix_like=bool(keeps_gain and within), matrix_like_all_only=bool(keeps_gain and within_all), rung=rung,
                d_rung=per)


def gain_kept(rd: pd.DataFrame, kind: str) -> dict:
    """The share of the gain kept (E124.gain_share, all recordings, readout) on all eight wirings, S, T and each
    wiring."""
    out = {u: E124.gain_share(rd, kind, u) for u in ('all8',) + POPS}
    out['per_wiring'] = {w: E124.gain_share(rd, kind, w) for w in WIRINGS}
    return out


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """An arm's class (E125's: E123's keeps the gain, "does what the matrix does", the corner rule), with the gain kept
    and the realism beside."""
    c = E123.classify(rd, runs, arm)
    if not c.get('complete'):
        return c
    c.update(E125.corners(runs, arm))
    c.update(matrix_like(rd, arm, bool(c['keeps_gain'])))
    cls = 'does what the matrix does' if c['matrix_like'] else ('keeps the gain' if c['keeps_gain'] else 'neither')
    if c['corner']:
        c['cls_before_corner'] = cls
        cls = 'corner (no class)'
    c['cls'] = cls
    sub = runs[runs['kind'] == arm]
    c.update(form=form_of(arm), dose=dose_of(arm), thr_below_60=sorted(sub.loc[sub['thr_below_60'], 'wiring']),
             fixed_thr_mean=float(sub['fixed_thr'].mean()), fixed_thr_min=float(sub['fixed_thr'].min()),
             wAPLKC_scale_mean=float(sub['wAPLKC_scale'].mean()), gain_kept=gain_kept(rd, arm))
    if 'apl_factor' in sub:
        c['apl_factor_mean'] = float(sub['apl_factor'].astype(float).mean())
    return c


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    """The "differ" rule at a d: `a` and `b` differ only if one is lower (calibrated pattern part) on at least
    DIFFER_WIRINGS of the eight wirings against all recordings, and lower on both populations against every one of the
    five targets. 'verdict': '<a> lower', '<b> lower', 'not told apart', 'no class (corner)' (either arm with corners
    on more than two wirings) or 'incomplete'."""
    if not (complete(runs, a) and complete(runs, b)):
        return dict(a=a, b=b, verdict='incomplete')
    per_w = {w: _v(rd, a, w, 'all_pattern') - _v(rd, b, w, 'all_pattern') for w in WIRINGS}
    per_t = {t: {u: _v(rd, a, u, f'{t}_pattern') - _v(rd, b, u, f'{t}_pattern') for u in POPS} for t in TARGETS}
    n_a, n_b = sum(v < 0 for v in per_w.values()), sum(v > 0 for v in per_w.values())
    a_all = all(v < 0 for t in per_t.values() for v in t.values())
    b_all = all(v > 0 for t in per_t.values() for v in t.values())
    corner = [x for x in (a, b) if E125.corners(runs, x)['corner']]
    if corner:
        verdict = 'no class (corner)'
    elif n_a >= DIFFER_WIRINGS and a_all:
        verdict = f'{a} lower'
    elif n_b >= DIFFER_WIRINGS and b_all:
        verdict = f'{b} lower'
    else:
        verdict = 'not told apart'
    return dict(a=a, b=b, dose=dose_of(a), verdict=verdict, n_wirings_a_lower=int(n_a), n_wirings_b_lower=int(n_b),
                a_lower_on_all_targets=bool(a_all), b_lower_on_all_targets=bool(b_all), corner_arms=corner,
                d_per_wiring=per_w, d_per_target=per_t,
                d_S4_all=per_t['all']['S4'], d_T4_all=per_t['all']['T4'])


def gate(classes: Mapping[str, dict]) -> dict:
    """Step 2's gate: TY75 complete, keeps the gain by E125's class, no corner, and no threshold under 60 Hz on any
    wiring. A TY75 that keeps the gain but is more than 0.003 above E124's rung still opens it (the entry's note)."""
    c = classes.get('TY75', {})
    if not c.get('complete'):
        return dict(complete=False, open=False, why='TY75 not run on all eight wirings')
    if c.get('corner'):
        return dict(complete=True, open=False, why=f"TY75 is a corner (wirings {c.get('corner_wirings')})")
    if not c.get('keeps_gain'):
        return dict(complete=True, open=False, why=f"TY75 does not keep the gain (class {c.get('cls')})")
    if c.get('thr_below_60'):
        return dict(complete=True, open=False, why=f"TY75 needs a threshold under 60 Hz on {c['thr_below_60']}")
    return dict(complete=True, open=True, why=f"TY75 keeps the gain (class {c.get('cls')}), no threshold under 60 Hz "
                                              f"(lowest {c.get('fixed_thr_min'):.1f} Hz)")


def outcomes(classes: Mapping[str, dict], comps: Mapping[str, dict], g: Mapping) -> dict:
    """The entry's outcomes as flags (the module docstring)."""
    ty = classes.get('TY75', {})
    ok = bool(ty.get('complete')) and not ty.get('corner')
    no_low = not ty.get('thr_below_60')
    a = bool(ok and ty.get('matrix_like') and no_low)
    b = bool(ok and (not ty.get('keeps_gain') or not no_low))
    between = bool(ok and ty.get('keeps_gain') and no_low and not ty.get('matrix_like'))
    rg = {r: classes.get(r, {}) for r in STEP2_ARMS}
    step2_done = all(c.get('complete') for c in rg.values())
    ty_lower = [r for r in STEP2_ARMS if comps.get(f'{r}|TY75', {}).get('verdict') == 'TY75 lower']
    rg_keep = {r: (bool(c.get('keeps_gain')) and not c.get('corner')) if c.get('complete') else None
               for r, c in rg.items()}
    out = dict(a=a, b=b, between_a_b=between, gate_open=bool(g.get('open')), step2_complete=step2_done,
               draws_keep_gain=rg_keep, ty75_lower_than=ty_lower,
               types_needed=bool(a and step2_done and (not all(rg_keep.values()) or bool(ty_lower))) if step2_done else None,
               c=bool(g.get('open') and step2_done and ty.get('keeps_gain') and all(rg_keep.values()) and not ty_lower)
               if step2_done else None,
               d=any(classes.get(x, {}).get('corner') for x in ARMS))
    return out


def arm_table(rd: pd.DataFrame, runs: pd.DataFrame, classes: Mapping[str, dict]) -> List[dict]:
    """Per arm and reference: the calibrated pattern part on S and T (all recordings), the gain kept on all8 / S / T,
    the mean APL effect (readout and raw), the threshold, APL->KC, the wirings at a bound and under 60 Hz, the class."""
    out = []
    kinds = ['family', 'E110_A', 'HBoff0'] + list(LADDER.values()) + [E123_U] + list(E125_ARMS) + list(ARMS)
    for k in kinds:
        sub = runs[runs['kind'] == k]
        if not complete(runs, k):
            out.append(dict(kind=k, complete=False))
            continue
        r = dict(kind=k, complete=True, cls=classes.get(k, {}).get('cls'))
        for u in POPS + ('all8',):
            r[f'pattern_cal_{u}'] = _v(rd, k, u, 'all_pattern')
            r[f'gain_share_{u}'] = E124.gain_share(rd, k, u)
        r['mean_effect_cal_all8'] = _v(rd, k, 'all8', 'mean_effect')
        r['mean_effect_raw_all8'] = _v(rd, k, 'all8', 'mean_effect', 'raw')
        r['fixed_thr_mean'] = float(sub['fixed_thr'].astype(float).mean())
        r['fixed_thr_min'] = float(sub['fixed_thr'].astype(float).min())
        r['wAPLKC_scale_mean'] = float(sub['wAPLKC_scale'].astype(float).mean())
        r['n_thr_below_60'] = int(sub['thr_below_60'].astype(bool).sum())
        r['n_outside_bands'] = int((~sub['in_bands'].astype(bool)).sum())
        if 'at_bound_apl' in sub:
            r['n_at_bound'] = int((sub['at_bound_apl'].eq(True)
                                   | sub.get('at_bound_thr', pd.Series(False, index=sub.index)).eq(True)).sum())
        if k in E125_ARMS + ARMS and dose_of(k) in LADDER:
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
    classes = {a: classify(rd, runs, a) for a in E125_ARMS + ARMS if a in present}
    refs = {k: E123.classify(rd, runs, k) for k in ('HBoff0',) + tuple(LADDER.values()) + (E123_U,) if k in present}
    for k, c in refs.items():
        if c.get('complete'):
            c['gain_kept'] = gain_kept(rd, k)
    rt = runs_table()
    bad = {}
    for a in set(rt.get('arm', [])):
        m = write_matrices(a)
        for c in ('on', 'off'):
            s = sorted(set(rt.loc[rt['arm'] == a, f'matrix_{c}_sha256']))
            if s != [m[c][2]]:
                bad[f'{a}_{c}'] = s
    if bad:
        raise RuntimeError(f"E126 score: an arm's matrix differs between hosts or from this host's: {bad}")
    comps = {f'{a}|{b}': differ(rd, runs, a, b) for a, b in COMPARISONS if a in present and b in present}
    g = gate(classes)
    chk = E122._read_json(out / 'check.json')
    at = arm_table(rd, runs, classes)
    v = dict(generated=_now(), host=_host(), check=dict(passes=chk.get('passes')),
             half_way=E124.half_way_check(rd), rungs=E125.rungs_check(rd), classes=classes, references=refs,
             gate=g, comparisons=comps, outcomes=outcomes(classes, comps, g), arm_table=at,
             fitted_reference_E110_B={u: _v(rd, 'E110_B', u, 'all_pattern') for u in POPS + ('all8',)}
             if 'E110_B' in present else None,
             step1_tags=sorted(t for t in rt['tag'] if str(t).startswith('S1_')),
             step2_tags=sorted(t for t in rt['tag'] if str(t).startswith('S2_')),
             rules=dict(within=WITHIN, matrix_like_doses=list(MATRIX_LIKE_DOSES), corner_wirings_max=CORNER_WIRINGS_MAX,
                        differ_wirings=DIFFER_WIRINGS, targets=list(TARGETS), populations=list(POPS),
                        realism_thr_hz=REALISM_THR_HZ, group_seeds=list(GROUP_SEEDS)),
             n_runs=int(len(runs)),
             step1_complete=all(classes.get(a, {}).get('complete') for a in STEP1_ARMS),
             step2_complete=all(classes.get(a, {}).get('complete') for a in STEP2_ARMS))
    E122._write_json(out / 'verdict.json', v)
    kinds = (['family', 'E110_A', 'E110_B', 'HBoff0'] + list(LADDER.values()) + [E123_U] + list(E125_ARMS)
             + [a for a in ARMS if a in present])
    md = ['# E126 comparison tables (generated)', '',
          'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
          f"Gate for step 2: {json.dumps(g)}", '',
          f"Outcomes: {json.dumps(v['outcomes'])}", '',
          f"Half-way points: {json.dumps(v['half_way'])}", '', f"E124's rungs: {json.dumps(v['rungs'])}", '',
          '## The arms beside the ladder', '',
          '| kind | rung | pattern S4 | pattern T4 | gain kept all8 | gain kept S4 | gain kept T4 | mean effect, readout | '
          'threshold (Hz) | lowest threshold (Hz) | APL->KC | at bound | thr < 60 Hz | outside bands | class |',
          '|---|---|' + '---:|' * 12 + '---|']
    for r in at:
        if not r.get('complete'):
            md.append(f"| {r['kind']} | | incomplete |" + ' |' * 12)
            continue
        md.append(f"| {r['kind']} | {r.get('rung', '')} | {r['pattern_cal_S4']:.4f} | {r['pattern_cal_T4']:.4f} | "
                  f"{r['gain_share_all8']:.2f} | {r['gain_share_S4']:.2f} | {r['gain_share_T4']:.2f} | "
                  f"{r['mean_effect_cal_all8']:+.3f} | {r['fixed_thr_mean']:.1f} | {r['fixed_thr_min']:.1f} | "
                  f"{r['wAPLKC_scale_mean']:.4g} | {r.get('n_at_bound', '')} | {r['n_thr_below_60']} | "
                  f"{r['n_outside_bands']} | {r.get('cls') or ''} |")
    md += ['', '## The "differ" comparisons', '',
           '| a | b | d | wirings a lower | wirings b lower | a lower on all targets | b lower on all targets | '
           'a - b, S4 | a - b, T4 | verdict |', '|---|---|---:|---:|---:|---|---|---:|---:|---|']
    for c in comps.values():
        if c['verdict'] == 'incomplete':
            md.append(f"| {c['a']} | {c['b']} | | | | | | | | incomplete |")
            continue
        md.append(f"| {c['a']} | {c['b']} | {c['dose']:g} | {c['n_wirings_a_lower']} | {c['n_wirings_b_lower']} | "
                  f"{c['a_lower_on_all_targets']} | {c['b_lower_on_all_targets']} | {c['d_S4_all']:+.4f} | "
                  f"{c['d_T4_all']:+.4f} | {c['verdict']} |")
    md.append('')
    for u in ('S4', 'T4', 'all8'):
        md += [f'## {u}', '', E123.comparison_table(rd, kinds, u), '']
    md += ['## Runs and bands', '', E123.runs_bands_table(runs), '']
    (out / 'comparison_table.md').write_text('\n'.join(md))
    E122._write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT,
                                                         commits=H.code_commits()))
    E122._write_json(out / 'provenance.json', dict(
        generated=_now(), host=_host(), script=SCRIPT, experiment=EXPERIMENT,
        model=H.model_provenance(orn_deltas=E83.wiring_parts('real')['orn']),
        design=dict(step1=list(STEP1_ARMS), step2=list(STEP2_ARMS), group_seeds=list(GROUP_SEEDS),
                    arms={a: arm_label(a) for a in ARMS}, references_read=['family', 'E110_A', 'E110_B', 'HBoff0']
                    + list(LADDER.values()) + [E123_U] + list(E125_ARMS)),
        run_provenance=sorted(p.name for p in out.glob('provenance_[12]_*.json'))))
    print(json.dumps({**{a: c.get('cls') for a, c in classes.items()}, 'gate': g, 'outcomes': v['outcomes']}),
          flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every E126 run, absolute parameters) and model_pair_distances_long.csv (raw spike counts)."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    mine = runs[runs['folder'].astype(str).str.startswith(EXPERIMENT)]
    rt = runs_table().drop_duplicates('tag', keep='last').set_index('tag')
    for r in mine.itertuples():
        o = outs[(r.kind, r.wiring)]
        t = rt.loc[r.tag]
        rows.append(dict(tag=r.tag, variant=r.tag, step=step_of(r.kind), arm=r.kind, arm_label=arm_label(r.kind),
                         form=form_of(r.kind), dose=dose_of(r.kind), group_seed=seed_of(r.kind), base='bouton',
                         k=E124.K, wiring=r.wiring,
                         kind=f"C2's family's {r.wiring}: {arm_label(r.kind)}, at the fit's scales, both response "
                              'rates held',
                         ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         kckc_matrix_on=f'matrices/{t["matrix_on"]}', kckc_matrix_off=f'matrices/{t["matrix_off"]}',
                         kckc_matrix_on_sha256=t['matrix_on_sha256'], kckc_matrix_off_sha256=t['matrix_off_sha256'],
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, kckc_scale_off=r.wKCKC_scale, rate_on=r.rate_on, rate_off=r.rate_off,
                         max_count_on=r.max_count_on, max_count_off=r.max_count_off, realism_flag=r.realism_flag,
                         thr_below_60=bool(r.thr_below_60), subdir=f'step{step_of(r.kind)}'))
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', 'step1', 'step2', 'score', 'variants', 'plan'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    ap.add_argument('--arms', nargs='+', default=None, choices=list(ARMS))
    a = ap.parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'check':
        part_check()
    elif a.part in ('step1', 'step2'):
        step = int(a.part[-1])
        part_runs(step, a.arms, ws, lane=f'{a.part}-{lane_name(ws)}')
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'plan':
        for arm, w in planned(a.arms or ARMS, ws):
            print(arm, w, tag_of(arm, w))


if __name__ == '__main__':
    main()
