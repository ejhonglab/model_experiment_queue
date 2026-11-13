#!/usr/bin/env python3
"""E89: the co-response read with a checked estimate, and the breadth weighting inside the models' APL.

The entry is item E89 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-13 03:20). E88's factors throughout (R recruitment, L co-response, G amplitude factor, O = R * L, s = R * L * G;
`E88_responder_sets`), C2's family (the breadth-2 APL-effect fit of E67) on E72 part B's eight wirings, S = real, P1,
P2, M1 and T = P3, P4, M2, M3, halves A / B of the recordings E77's (`recording_halves()`, seed 0) and C / D the second
split's (seed 1).

Parts (CLI `--part`):
- A (no simulations): A0 the candidate estimates of the co-response on known truth (E88's pseudo-recordings, ten a
  wiring, three trials, seeds SEEDS) at three noise levels: L1 the reading's own noise (E73's calibration), L2 the noise
  raised to the flies' repeatability, independent between odors (E88's reading 1b), and L3 the flies' repeatability with
  a share of each KC's noise common to the odors of a repeat, set per wiring and condition so that the pseudo-recordings'
  same-repeat gap is the flies' (E81's statistic), and at half and one and a half times that share (L3h, L3x). The
  candidates: (a) E88's two estimates of O ('tm' trial means, 'nr' trial noise removed), (b) their geometric mean
  ('gm'), (c) latent classes after deconvolution ('lc'), (d) a correction by each odor's agreement between repeats
  ('hf'); R always from the trial means and s from the noise-removed estimate. The choice; A1 E88's part A read again
  with the chosen candidate (factors, cells, the split, the groups, the base); A2 E88's B1 with it. estimator_check.csv,
  estimator_choice.json, factors.csv, cells.csv, changes.csv, verdict_A.json.
- B (no simulations): the saved families in which only the weighting differs (E66's part C, E67's exponent fits, E12's
  and E19's grids apart) by the main target and E88's archive scan; weighting_archive.csv, verdict_B.json.
- C0 / C1 / C2 / C3 (simulations): twelve weightings of APL's inhibition by a KC's breadth (WEIGHTINGS) on S's four
  wirings, on the fits as they are (C1) and on the base with the flies' recruitment (C2); identities first (C0); the
  choice per base and half; the confirmation on T (C3) where there is something to confirm, with E83's gate and the
  permuted control. Cscore: readings_C.csv, coresponse_C.csv, choice_C.csv, gate_C.csv, verdict_C.json.
- variants: variants.csv and model_pair_distances_long.csv of C's runs, for the archive.
- table: comparison_table.md.

Outputs: results/E89_apl_breadth_weighting/; spike counts and responses under outputs/E89_apl_breadth_weighting/C/.

Usage (from the repository root):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E89_apl_breadth_weighting.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E89_apl_breadth_weighting.py --part C1
"""
from __future__ import annotations

import argparse
import json
import os
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E87_inhibition_form as E87
from experiments import E88_responder_sets as E88

EXPERIMENT = 'E89_apl_breadth_weighting'

# ---------------------------------------------------------------------------
# part A: the settings and the thresholds (the entry's Design; fixed before the flies are read with any candidate)

SEEDS = (E88.SEED, E88.SEED + 1)        # E88's pseudo-recordings' seed and a second beside it
N_PSEUDO, N_TRIALS = E88.N_PSEUDO, E88.N_TRIALS
EST_TOL = E88.EST_TOL                   # 0.03: a factor's group mean off by at most this, in units of the similarity
CANDIDATES = ('tm', 'nr', 'gm', 'lc', 'hf')
CANDIDATE_NAMES = {'tm': 'trial means (E88 reading 1)', 'nr': 'trial noise removed (E88 reading 2)',
                   'gm': 'the geometric mean of the two', 'lc': 'latent classes after deconvolution',
                   'hf': "a correction by each odor's agreement between repeats"}
LEVELS = ('L1', 'L2', 'L3h', 'L3', 'L3x')
SHARE_SCALE = {'L3h': 0.5, 'L3': 1.0, 'L3x': 1.5}
BINDING_LEVELS = ('L1', 'L3h', 'L3', 'L3x')
WHO_FACTORS = ('R', 'L', 's')           # passes for who responds
FULL_FACTORS = ('R', 'L', 's', 'G')     # passes in full
SHARE_BOUNDS = (0.0, 0.95)              # L3: the search range of the shared share of the noise variance
# (c)'s free choice, the classes' edges in z beside the criterion, set on the firing-ceiling population then left alone
LC_EDGES = {'crit': (), 'coarse': (0.25, 0.6, 1.5, 2.5, 4.0), 'fine': (0.0, 0.35, 0.7, 1.3, 1.7, 2.3, 3.2, 4.5)}
LC_STEP, LC_NPMLE_ITERS, LC_EM_ITERS = 0.05, 300, 200
LC_SET_POPULATION = 'fc'
CHECK_POPULATIONS = ('b2', 'base')      # the check is read on these
A_POPULATIONS = ('b2', 'fc', 'S_unchanged', 'base')
# E88's thresholds, unchanged
G_MIN, G_RATIO, GROUP_RATIO, D2_BASE_GAIN, B1_FACTOR = E88.G_MIN, E88.G_RATIO, E88.GROUP_RATIO, E88.D2_BASE_GAIN, E88.B1_FACTOR

# ---------------------------------------------------------------------------
# part B

B_GAIN = 0.003                          # the member nearest the flies' change lower than the exponent-2 member by this

# ---------------------------------------------------------------------------
# part C: the weightings and the thresholds

MEASURES = ('F', 'S', 'N')
MEASURE_NAMES = {'F': "the family's: the KC's spike total without APL in the bouton base",
                 'S': 'the summed response: the KC\'s summed spike count over the odors it responds to (two or more spikes) '
                      'in the run the APL acts on, without APL',
                 'N': 'the number of odors it responds to in that run'}
EXPONENTS = {'F': (0.25, 0.5, 1.0), 'S': (0.25, 0.5, 1.0, 2.0), 'N': (0.25, 0.5, 1.0, 2.0)}
REFERENCE = ('F', 2.0)
NO_WEIGHTING = ('F', 0.0)               # exponent 0: every factor 1 (the measure does not matter)


def weightings() -> List[Tuple[str, float]]:
    """The twelve weightings (measure, exponent), no weighting first, then the reference (F, 2)."""
    return [NO_WEIGHTING] + [(m, k) for m in MEASURES for k in EXPONENTS[m]] + [REFERENCE]


def wtag(w: Tuple[str, float]) -> str:
    """A weighting's tag: 'none', or '<measure>_k<exponent>' ('F_k2' the reference)."""
    m, k = w
    return 'none' if k == 0 else f'{m}_k{k:g}'


ARMS = ('C1', 'C2')                     # the fits as they are; the base with the flies' recruitment
RAW_CRIT = E87.RAW_CRIT                 # a model response is two or more spikes
CHOICE_MIN_GAIN = 0.0015                # below this gain on S against half X there is nothing to confirm
GATE_GAIN, GATE_BORDER = 0.003, (0.0015, 0.0045)
GATE_SINGLES = 6                        # of T's eight wiring-by-direction cells
GATE_RAW_SLACK, GATE_MEAN_SLACK, FLIES_MEAN = 0.005, 0.01, 0.178
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075
CONTROL_RATIO = 0.5
CONTROL_SEEDS = (1, 2)
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
SPLITS = {0: ('A', 'B'), 1: ('C', 'D')}  # E77's recording_halves(seed)
UNCHANGED_T = {'A': 0.0565, 'B': 0.0538}  # T's unchanged fits against each half (E83; the entry's references)
UNCHANGED_S = {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
C2_UNITS = 'noise'                      # the base's per-odor recruitment, as E87's D2 and E88's D1

_CACHE: Dict[tuple, object] = {}


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


# ---------------------------------------------------------------------------
# part A: the candidates as functions of a recording's single trials

def lc_overlap(trials: np.ndarray, crit: float, edges: Sequence[float] = (), se: Optional[np.ndarray] = None
               ) -> Tuple[np.ndarray, np.ndarray]:
    """(c) latent classes after deconvolution: (O, f) of one recording's single trials (repeats, KCs, odors; graded z).
    Per odor the distribution of a KC's mean response with the trial noise deconvolved (NPMLE on a grid, each KC's own
    trial SD pooled over odors), cut into classes at `crit` and at `edges`; per pair the table of class weights by
    maximum likelihood (EM); O (triu order) the share of KCs above `crit` for both over sqrt(f_a f_b), f per odor."""
    raise NotImplementedError


def hf_overlap(binary: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    """(d) a correction by each odor's agreement between repeats: (O, f, params) of one recording's thresholded single
    trials (repeats, KCs, odors; binary). Per odor two classes (responder or not) with a hit rate h and a false-alarm
    rate e, fitted with their share f to the numbers of repeats responding (a binomial mixture, EM); per pair the share
    of KCs in both classes n11 solved from the mean product of the two odors over different repeats, which is
    n11 (h_a - e_a)(h_b - e_b) + f_a h_a e_b + f_b e_a h_b + (1 - f_a - f_b) e_a e_b; O = n11 / sqrt(f_a f_b)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A: the model's readings at the three noise levels

class ModelReading89(E88.ModelReading):
    """E88's `ModelReading` with a seed (SEEDS) and the third noise level: a per-KC, per-repeat term common to the odors
    of a repeat, its share of each KC's noise variance set per condition so that the pseudo-recordings' same-repeat gap
    is the flies'. `level_data(level, c)` gives the trial means, single trials and both criteria of a level;
    `candidate_factors(level, cand)` the population-ready per-pseudo-recording factors."""

    def __init__(self, o, label: str = '', seed: int = E88.SEED, **kw):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# part A: the check, the choice, the readings

def estimator_check(members: Sequence['ModelReading89'], levels: Sequence[str] = LEVELS) -> pd.DataFrame:
    """Each candidate's factors (the population's) minus the truth by group, condition and level: 'err' in units of the
    similarity (E88's convention), 'off' beyond EST_TOL."""
    raise NotImplementedError


def choose_candidate(ec: pd.DataFrame, population: str) -> dict:
    """The thresholds of the entry on one population: per candidate whether it passes for who responds (R, L and s
    inside the tolerance at every binding level, both seeds) and in full (G too); the chosen one (in full first, then
    the smallest largest error of L at L3)."""
    raise NotImplementedError


def part_a(out: Path, n_workers: int = 4) -> dict:
    """A0, the choice, A1 and A2; verdict_A.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B

def part_b(out: Path) -> dict:
    """The weighting in the saved models: weighting_archive.csv and verdict_B.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part C: the weightings

def _kc_ids(df: pd.DataFrame) -> pd.Index:
    return pd.Index(df.index.get_level_values(H.KC_ID_LEVEL))


def measure_values(measure: str, off_counts: pd.DataFrame, family_breadth: pd.Series) -> pd.Series:
    """A KC's breadth by `measure` (indexed by KC id, the outputs' KC order for S and N): 'F' the family's breadth
    (`family_breadth`, E83's `wiring_parts(label)['breadth']`: the KC's spike total without APL in the bouton base), 'S'
    the summed spike count over the odors with two or more spikes in `off_counts` (the run the APL acts on, without APL),
    'N' the number of those odors; S and N are 0 for a KC without a response."""
    if measure == 'F':
        return family_breadth.astype(float)
    c = off_counts.to_numpy(dtype=float)
    resp = c > RAW_CRIT
    if measure == 'S':
        v = np.where(resp, c, 0.0).sum(axis=1)
    elif measure == 'N':
        v = resp.sum(axis=1).astype(float)
    else:
        raise ValueError(measure)
    return pd.Series(v, index=_kc_ids(off_counts), name=measure)


def weighting_factors(m: pd.Series, k: float) -> pd.Series:
    """((m + 1) / mean over the KCs of (m + 1)) ** k per KC (as `E10_joint_fit.breadth_weights` forms it)."""
    bp1 = m.astype(float) + 1.0
    return (bp1 / bp1.mean()) ** float(k)


def apl_kc_weights(wnorm_aplkc: pd.Series, m: pd.Series, k: float) -> pd.Series:
    """The APL->KC weights (per claw) of a weighting: the connectome's normalized weights times the KC's factor, by
    `E10_joint_fit.breadth_weights` itself (so F at exponent 2 gives the fit's weights exactly)."""
    from experiments import E10_joint_fit as J
    return J.breadth_weights(wnorm_aplkc, m.astype(float), float(k))


def factor_summary(f: pd.Series, responders: pd.Index) -> dict:
    """The factors' range: percentiles over all KCs and over the KCs that respond without APL (`responders`, KC ids)."""
    row = {}
    for lab, v in (('all', f.to_numpy(dtype=float)), ('resp', f.reindex(responders).to_numpy(dtype=float))):
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        q = np.quantile(v, [0.0, 0.1, 0.5, 0.9, 0.99, 1.0])
        row.update({f'fac_{lab}_{n}': float(x) for n, x in zip(('min', 'p10', 'p50', 'p90', 'p99', 'max'), q)})
        row[f'fac_{lab}_p90_over_p10'] = float(q[3] / q[1]) if q[1] > 0 else np.nan
    return row


def _responders(off_counts: pd.DataFrame) -> pd.Index:
    resp = (off_counts.to_numpy(dtype=float) > RAW_CRIT).any(axis=1)
    return _kc_ids(off_counts)[resp]


def _c_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT) / 'C'
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_tag(arm: str, w: Tuple[str, float], wiring: str, extra: str = '') -> str:
    return f'{arm}_{wtag(w)}{extra}_{wiring}'


def save_run(tag: str, o: Dict[str, tuple]) -> None:
    d = _c_dir()
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    """A saved C run ({'on', 'off'}: (responses, spike counts)) from the local outputs, else the archive; or None."""
    for d in (_c_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'C'):
        out = {}
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            ps = [d / f'{tag}_{cc}_{k}.parquet' for k in ('responses', 'spike_counts')]
            if not all(p.exists() for p in ps):
                out = None
                break
            out[c] = tuple(pd.read_parquet(p) for p in ps)
        if out:
            return out
    return None


def _commits_clean() -> dict:
    """H.code_commits() at the start of a simulation part, which must be clean (the rows carry this snapshot: the code
    the process loaded)."""
    cc = H.code_commits()
    dirty = [k for k, v in cc.items() if k.endswith('_dirty') and v]
    if dirty:
        raise RuntimeError(f'commit first: {dirty} (rows record the commits of the code that ran)')
    return cc


def _append(path: Path, row: dict) -> None:
    from experiments import E83_cross_wiring_fit as E
    E._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def _realism_row(o: Dict[str, tuple]) -> dict:
    from experiments import E83_cross_wiring_fit as E
    return E.realism(o)


# C1: the fits as they are (E83's WiringModel, the APL->KC weights set in place)

def set_weights(model, wk: Optional[pd.Series]) -> None:
    """Set a WiringModel's APL->KC weights (unscaled, per claw) in place: `wk` (indexed as the model's), or the fit's
    for None."""
    if wk is None:
        model._kc[...] = model._kc0
        return
    if not wk.index.equals(model.wAPLKC.index):
        raise RuntimeError(f'{model.label}: the weights are not indexed as the model\'s APL->KC weights')
    model._kc[...] = wk.to_numpy(dtype=float).reshape(model._kc0.shape)


def weighting_reading(model, wk: pd.Series, label: str = '') -> dict:
    """One weighting on a wiring at its fit (C1): the APL->KC weights set to `wk`, the APL scale (APL->KC, KC->APL,
    bouton->APL tied, E83's TIED_ON) bisected (E69's `bisect_factor`, E83's APL_BRACKET and APL_STEP, from 1) to the
    fit's APL-on response count; the APL-off run is the saved fit's. The weights are reset on return."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E83_cross_wiring_fit as E
    t0 = time.time()
    set_weights(model, wk)
    try:
        x, on, hist, at_bound = E69.bisect_factor(model.run_on, model.target_on, model.size, E.APL_BRACKET, E.APL_STEP,
                                                  f'{label} APL scale')
    finally:
        set_weights(model, None)
    p = model.p
    return dict(on=on, off=model.saved['off'], apl_factor=float(x), at_bound_apl=bool(at_bound), n_sims=len(hist),
                params=dict(p, **{k: p[k] * x for k in E.TIED_ON}), seconds=time.time() - t0, history=hist)


def _c1_measures(model) -> Dict[str, pd.Series]:
    from experiments import E83_cross_wiring_fit as E
    off = model.saved['off'][1]
    br = E.wiring_parts(model.label)['breadth']
    return {m: measure_values(m, off, br) for m in MEASURES}


def part_c1(out: Path, wirings: Sequence[str], weights: Optional[Sequence[Tuple[str, float]]] = None,
            runs_name: str = 'runs_C1.csv', extra: str = '', permute_seed: Optional[int] = None) -> pd.DataFrame:
    """C1 (and C3's runs of C1's kind): each weighting on each wiring at its fit's parameters, the APL scale bisected to
    the fit's APL-on count; the APL-off run the saved one. The reference (F at exponent 2) is C0's identity (the saved
    fit, the scale 1). With `permute_seed`, the weighting's per-KC factors permuted among the KCs of each KC type (E80's
    `permute_within_type`; the control of C3). Resumable (rows in `runs_name`, runs under outputs/.../C/)."""
    import gc
    from experiments import E80_screen as S80
    from experiments import E83_cross_wiring_fit as E
    cc = _commits_clean()
    weights = list(weights or weightings())
    path = out / runs_name
    for wiring in wirings:
        todo = [w for w in weights if run_tag('C1', w, wiring, extra) not in _done(path)]
        if not todo:
            continue
        model = E.WiringModel(wiring, check=False)
        ms = _c1_measures(model)
        wn = E.wiring_parts(wiring)['wnorm']['wAPLKC']
        resp = _responders(model.saved['off'][1])
        kc_type = model.saved['off'][1].index.to_frame(index=False).drop_duplicates(H.KC_ID_LEVEL).set_index(H.KC_ID_LEVEL)[H.KC_TYPE_LEVEL]
        for w in todo:
            tag = run_tag('C1', w, wiring, extra)
            m, k = w
            mv = ms[m]
            fac = weighting_factors(mv, k)
            if permute_seed is not None:
                fi = pd.MultiIndex.from_arrays([fac.index, kc_type.reindex(fac.index).to_numpy()], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
                fac = pd.Series(S80.permute_within_type(pd.Series(fac.to_numpy(), index=fi), permute_seed).to_numpy(), index=fac.index)
                # the permuted factors as the measure behind breadth_weights: m' + 1 = fac' * mean(m + 1) ** (1 / k)
                wk = wn * fac.reindex(wn.index.get_level_values(H.KC_ID_LEVEL)).fillna(1.0).to_numpy()
            else:
                wk = apl_kc_weights(wn, mv, k)
            res = weighting_reading(model, wk, label=tag)
            o = {'on': res['on'], 'off': res['off']}
            save_run(tag, o)
            row = dict(part='C1', arm='C1', tag=tag, weighting=wtag(w), measure=m if k else 'none', exponent=k, wiring=wiring,
                       permute_seed=permute_seed, fit_folder=model.folder, fit_tag=model.tag, apl_factor=res['apl_factor'],
                       apl_min=res['apl_factor'], apl_max=res['apl_factor'], at_bound=res['at_bound_apl'],
                       n_on=E._n(res['on'][0]), target_on=model.target_on, n_off=E._n(res['off'][0]), target_off=model.target_off,
                       n_sims=res['n_sims'], seconds=res['seconds'], **{kk: res['params'][kk] for kk in E.PARAM_KEYS},
                       wAPLKC_sha256=E80_sha(wk), K_fit=model.K, **factor_summary(fac, resp), **cc)
            row.update(_realism_row(o))
            _append(path, row)
            print(f'E89 {tag}: APL scale x {res["apl_factor"]:.4g}, {row["n_on"]} / {model.target_on} responses, '
                  f'{res["n_sims"]} sims, {res["seconds"]:.0f}s', flush=True)
        del model
        gc.collect()
    return pd.read_csv(path)


def E80_sha(a) -> str:
    from experiments import E80_screen as S80
    return S80._arr_sha256(np.asarray(a, dtype=float))


# C2: the base with the flies' recruitment (E87's B3Pool, the APL->KC weights set in its workers)

def _c2_worker(conn, label: str, parts: dict, fit: dict, js: List[int]) -> None:
    """A C2 worker: wiring `label`'s simulator on the ORN columns `js` (as E87's `_d2_worker`), then for each task its
    APL->KC weights set in place ('wAPLKC', per claw; None the fit's) and E80's B3 one-odor task (`_odor_task`) run,
    until None."""
    try:
        from experiments import E83_cross_wiring_fit as E
        from experiments import E80_screen as S80
        E._PARTS[label] = parts
        E._FITS[label] = fit
        sim, info = E.build_simulator(label, orn=parts['orn'].iloc[:, js], breadth=parts['breadth'])
        cells = info['cells']
        mp, order = cells['mp'].cell_contents, np.asarray(cells['odor_order'].cell_contents, dtype=int)
        names = [str(c) for c in cells['odor_index'].cell_contents]
        kc = cells['wAPLKC_unscaled'].cell_contents
        kc0 = kc.copy()
        n_kc = int(fit['saved']['on'][0].shape[0])
        conn.send(('ready', names))
    except Exception as e:
        conn.send(('error', repr(e)))
        return
    while True:
        msg = conn.recv()
        if msg is None:
            break
        try:
            wk = msg.pop('wAPLKC', None)
            kc[...] = kc0 if wk is None else np.asarray(wk, dtype=float).reshape(kc0.shape)
            try:
                res = S80._odor_task(sim, mp, order, names, dict(info['params']), n_kc, msg)
            finally:
                kc[...] = kc0
            conn.send(('ok', res))
        except Exception as e:
            conn.send(('error', repr(e)))


class C2Pool(E87.B3Pool):
    """E87's B3Pool with this script's worker: `run_task(task, wk)` runs E80's B3 task with the APL->KC weights `wk`
    (a per-claw Series indexed as the fit's, or None for the fit's) and returns ({'off', 'on1', 'on'}: (responses, spike
    counts)), per-odor records)."""

    _WORKER = staticmethod(_c2_worker)

    def run_task(self, task: dict, wk: Optional[pd.Series] = None):
        keys = ('thr', 'target_off', 'apl', 'target_on')
        arr = None if wk is None else wk.to_numpy(dtype=float)
        for conn, js in zip(self.conns, self.groups):
            conn.send(dict({k: (None if task.get(k) is None else [task[k][j] for j in js]) for k in keys},
                           label=task.get('label', ''), wAPLKC=arr))
        recs = [None] * len(self.cols)
        for conn, js in zip(self.conns, self.groups):
            status, v = conn.recv()
            if status != 'ok':
                raise RuntimeError(f'{self.label}: a worker failed: {v}')
            for j, rec in zip(js, v):
                recs[j] = rec
        from experiments import E83_cross_wiring_fit as E
        saved = E.saved_fit(self.label)['saved']
        outs = {}
        for key in ('off', 'on1', 'on'):
            src = saved['off' if key == 'off' else 'on']
            R, S = src[0].copy(), src[1].copy()
            for j, rec in enumerate(recs):
                if rec['odor'] != self.cols[j]:
                    raise RuntimeError(f"odor {j}: {rec['odor']} where {self.cols[j]} was expected")
                R.iloc[:, j] = rec[key][0][:, 0]
                S.iloc[:, j] = rec[key][1][:, 0]
            outs[key] = (R, S)
        return outs, recs


def base_thresholds(w: str, n_workers: int = 16) -> dict:
    """The base's per-odor threshold factors on wiring `w` ({'thr': {odor column: factor}, 'target_off', 'target_on',
    'cols'}): E88's (outputs/E88_responder_sets/D1/base_thr_<w>.json) for S's wirings; for the others solved as E88
    solved S's (E88's `base_thresholds`: E80's one-odor task with APL off, each odor's threshold bisected to the flies'
    count with APL output blocked (TNT), baseline-noise units) and saved to outputs/E89_apl_breadth_weighting/C/
    base_thr_<w>.json, after the identity of the one-odor path on that wiring (factors 1 reproduce the saved fit;
    identity_C.json)."""
    from experiments import E80_screen as S80
    from experiments import E83_cross_wiring_fit as E
    p88 = H.outputs_dir(E88.EXPERIMENT) / 'D1' / f'base_thr_{w}.json'
    if p88.exists():
        return json.loads(p88.read_text())
    path = _c_dir() / f'base_thr_{w}.json'
    if path.exists():
        return json.loads(path.read_text())
    fit = E.saved_fit(w)
    pool = C2Pool(w, n_workers)
    try:
        n = len(pool.cols)
        outs, _ = pool.run_task(dict(thr=[1.0] * n, apl=[1.0] * n, label=f'E89 one-odor identity {w}'))
        ns, nr = E._same({'on': outs['on1'], 'off': outs['off']}, fit['saved'])
        _identity_add(H.results_dir(EXPERIMENT), f'one_odor_path_{w}', dict(count_diffs=ns, response_diffs=nr))
        if ns or nr:
            raise RuntimeError(f'{w}: the one-odor path does not reproduce the saved fit')
        tg = S80.b3_targets(C2_UNITS, pool.cols, int(fit['saved']['on'][0].shape[0]))
        outs, recs = pool.run_task(dict(target_off=tg['target_off'], label=f'E89 base {w}'))
    finally:
        pool.close()
    thr = {r['odor']: float(r['thr']) for r in recs}
    path.write_text(json.dumps(dict(thr=thr, target_off=tg['target_off'], target_on=tg['target_on'], cols=list(pool.cols),
                                    at_bound=[bool(r['at_bound_thr']) for r in recs], **H.code_commits()), indent=1) + '\n')
    outs['off'][0].to_parquet(_c_dir() / f'base_{w}_acute_off_responses.parquet')
    outs['off'][1].to_parquet(_c_dir() / f'base_{w}_acute_off_spike_counts.parquet')
    return json.loads(path.read_text())


def base_off(w: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The base's APL-off run on wiring `w` (E88's D1 base for S's wirings, this script's for the others)."""
    for d in (H.outputs_dir(E88.EXPERIMENT) / 'D1', _c_dir(), H.shared_dir() / 'data' / 'model_outputs' / E88.EXPERIMENT / 'D1'):
        ps = [d / f'base_{w}_acute_off_{k}.parquet' for k in ('responses', 'spike_counts')]
        if all(p.exists() for p in ps):
            return tuple(pd.read_parquet(p) for p in ps)
    raise RuntimeError(f'no base APL-off run for {w}: run base_thresholds first')


def _identity_add(out: Path, key: str, d: dict) -> None:
    path = out / 'identity_C.json'
    v = json.loads(path.read_text()) if path.exists() else {}
    v[key] = dict(d, **H.code_commits())
    _write_json(path, v)


def part_c2(out: Path, wirings: Sequence[str], weights: Optional[Sequence[Tuple[str, float]]] = None,
            n_workers: int = 16, runs_name: str = 'runs_C2.csv', extra: str = '', permute_seed: Optional[int] = None,
            identity: bool = True) -> pd.DataFrame:
    """C2 (and C3's runs of C2's kind): on each wiring the base (per-odor threshold factors, `base_thresholds`), then each
    weighting with each odor's APL scale (tied, E80's B3) bisected to the control flies' count for that odor
    (`E80_screen.b3_targets`, baseline-noise units; E80's B3 bracket). The measures S and N from the base's APL-off run.
    The reference (F at exponent 2) first: C0's identity on S's wirings (with `identity`), whose APL-off and APL-on runs
    must equal E87's D2 run with per-odor APL scales exactly (E87 did not save the per-odor scales; the same bisection
    from the same start at the same thresholds recovers them, and they are saved here, base_apl_<wiring>.json). With
    `permute_seed`, the factors permuted within KC type (C3's control). Resumable (rows in `runs_name`)."""
    from experiments import E80_screen as S80
    from experiments import E83_cross_wiring_fit as E
    cc = _commits_clean()
    weights = list(weights or weightings())
    if REFERENCE in weights and permute_seed is None:
        weights = [REFERENCE] + [w for w in weights if w != REFERENCE]
    path = out / runs_name
    for wiring in wirings:
        todo = [w for w in weights if run_tag('C2', w, wiring, extra) not in _done(path)]
        if not todo:
            continue
        b = base_thresholds(wiring, n_workers)
        fit = E.saved_fit(wiring)
        parts = E.wiring_parts(wiring)
        off0 = base_off(wiring)
        ms = {m: measure_values(m, off0[1], parts['breadth']) for m in MEASURES}
        wn = parts['wnorm']['wAPLKC']
        resp = _responders(off0[1])
        kc_type = off0[1].index.to_frame(index=False).drop_duplicates(H.KC_ID_LEVEL).set_index(H.KC_ID_LEVEL)[H.KC_TYPE_LEVEL]
        pool = C2Pool(wiring, n_workers)
        try:
            thr = [float(b['thr'][c]) for c in pool.cols]
            tg = S80.b3_targets(C2_UNITS, pool.cols, int(fit['saved']['on'][0].shape[0]))
            for w in todo:
                tag = run_tag('C2', w, wiring, extra)
                m, k = w
                fac = weighting_factors(ms[m], k)
                if permute_seed is not None:
                    fi = pd.MultiIndex.from_arrays([fac.index, kc_type.reindex(fac.index).to_numpy()], names=[H.KC_ID_LEVEL, H.KC_TYPE_LEVEL])
                    fac = pd.Series(S80.permute_within_type(pd.Series(fac.to_numpy(), index=fi), permute_seed).to_numpy(), index=fac.index)
                    wk = wn * fac.reindex(wn.index.get_level_values(H.KC_ID_LEVEL)).fillna(1.0).to_numpy()
                else:
                    wk = apl_kc_weights(wn, ms[m], k)
                t0 = time.time()
                outs, recs = pool.run_task(dict(thr=thr, target_on=tg['target_on'], label=f'E89 {tag}'), wk)
                secs = time.time() - t0
                o = {'on': outs['on'], 'off': outs['off']}
                ns_off, nr_off = E._same({'on': outs['off'], 'off': outs['off']}, {'on': off0, 'off': off0})
                if ns_off or nr_off:
                    raise RuntimeError(f'{tag}: the APL-off run differs from the base\'s ({ns_off} counts)')
                if w == REFERENCE and permute_seed is None and identity and wiring in E.S_WIRINGS:
                    d2 = E87._load_d2(f'D2_b3_thrapl_{E87.D2_UNITS}_{wiring}')
                    ns, nr = E._same(o, d2)
                    _identity_add(out, f'C2_{wiring}', dict(count_diffs=ns, response_diffs=nr,
                                                            against=f'E87 D2_b3_thrapl_{E87.D2_UNITS}_{wiring}'))
                    print(f'E89 C0 (C2) identity {wiring}: {ns} count / {nr} response differences', flush=True)
                    if ns or nr:
                        raise RuntimeError(f"{wiring}: the base with F at exponent 2 does not reproduce E87's D2 run")
                if w == REFERENCE and permute_seed is None:
                    _write_json(_c_dir() / f'base_apl_{wiring}.json', dict(apl={r['odor']: r['apl'] for r in recs},
                                                                           thr={r['odor']: r['thr'] for r in recs}, **cc))
                save_run(tag, o)
                apl = np.array([r['apl'] for r in recs])
                row = dict(part='C2', arm='C2', tag=tag, weighting=wtag(w), measure=m if k else 'none', exponent=k, wiring=wiring,
                           permute_seed=permute_seed, fit_folder=fit['folder'], fit_tag=fit['tag'], apl_factor=float(np.exp(np.log(apl).mean())),
                           apl_min=float(apl.min()), apl_max=float(apl.max()),
                           at_bound=bool(any(r['at_bound_apl'] for r in recs)), n_on=E._n(outs['on'][0]), target_on=int(sum(tg['target_on'])),
                           n_off=E._n(outs['off'][0]), target_off=int(sum(b['target_off'])), n_sims=int(sum(r['n_sims'] for r in recs)),
                           seconds=secs, thr_min=float(min(thr)), thr_max=float(max(thr)),
                           apl_per_odor=json.dumps({r['odor']: r['apl'] for r in recs}),
                           thr_per_odor=json.dumps({r['odor']: r['thr'] for r in recs}),
                           **{kk: fit['p'][kk] for kk in E.PARAM_KEYS}, wAPLKC_sha256=E80_sha(wk), K_fit=fit['K'],
                           **factor_summary(fac, resp), **cc)
                row.update(_realism_row(o))
                _append(path, row)
                print(f'E89 {tag}: APL scales x {row["apl_min"]:.3g} to {row["apl_max"]:.3g}, {row["n_on"]} / {row["target_on"]} '
                      f'responses, {row["n_sims"]} sims, {secs:.0f}s', flush=True)
        finally:
            pool.close()
    return pd.read_csv(path)


def part_c0(out: Path, arm: str, wirings: Sequence[str], n_workers: int = 16) -> dict:
    """C0's identities. C1: on each wiring the simulator at its fit reproduces the saved APL-on and APL-off runs
    (E83's `WiringModel(check=True)`), and F at exponent 2 through the new code path (`apl_kc_weights`, set in place)
    reproduces the saved APL-on run exactly with the APL scale at 1; that run is C1's reference point. C2: the reference
    of `part_c2` (run first there) against E87's D2 run with per-odor APL scales. identity_C.json."""
    from experiments import E83_cross_wiring_fit as E
    if arm == 'C2':
        part_c2(out, wirings, [REFERENCE], n_workers)
        return json.loads((out / 'identity_C.json').read_text())
    cc = _commits_clean()
    path = out / 'runs_C1.csv'
    for wiring in wirings:
        tag = run_tag('C1', REFERENCE, wiring)
        if tag in _done(path):
            continue
        t0 = time.time()
        model = E.WiringModel(wiring, check=True)
        parts = E.wiring_parts(wiring)
        wk = apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], REFERENCE[1])
        same_w = bool(np.array_equal(wk.to_numpy(), model._kc0.ravel()))
        set_weights(model, wk)
        try:
            on = model.run_on(1.0)
        finally:
            set_weights(model, None)
        ns, nr = E._same({'on': on, 'off': model.saved['off']}, model.saved)
        _identity_add(out, f'C1_{wiring}', dict(weights_equal=same_w, count_diffs=ns, response_diffs=nr,
                                                against=f'{model.folder}:{model.tag}'))
        print(f'E89 C0 (C1) identity {wiring}: weights equal {same_w}, {ns} count / {nr} response differences', flush=True)
        if ns or nr or not same_w:
            raise RuntimeError(f'{wiring}: F at exponent 2 through the new path does not reproduce the saved fit')
        o = {'on': on, 'off': model.saved['off']}
        save_run(tag, o)
        fac = weighting_factors(parts['breadth'], REFERENCE[1])
        row = dict(part='C1', arm='C1', tag=tag, weighting=wtag(REFERENCE), measure='F', exponent=REFERENCE[1], wiring=wiring,
                   permute_seed=None, fit_folder=model.folder, fit_tag=model.tag, apl_factor=1.0, apl_min=1.0, apl_max=1.0,
                   at_bound=False, n_on=E._n(on[0]), target_on=model.target_on, n_off=E._n(model.saved['off'][0]),
                   target_off=model.target_off, n_sims=1, seconds=time.time() - t0, **{kk: model.p[kk] for kk in E.PARAM_KEYS},
                   wAPLKC_sha256=E80_sha(wk), K_fit=model.K, **factor_summary(fac, _responders(model.saved['off'][1])), **cc)
        row.update(_realism_row(o))
        _append(path, row)
        del model
    return json.loads((out / 'identity_C.json').read_text())


def score_c(out: Path) -> dict:
    """The readings, the co-response, the choice, the gate and the verdicts: readings_C.csv, coresponse_C.csv,
    choice_C.csv, gate_C.csv, verdict_C.json."""
    raise NotImplementedError


def part_c3(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """C3: the confirmation on T for each base and direction with something to confirm, and the permuted controls."""
    raise NotImplementedError


def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv and model_pair_distances_long.csv of C's runs."""
    raise NotImplementedError


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    raise NotImplementedError


def _wirings(arg: str) -> Tuple[str, ...]:
    from experiments import E83_cross_wiring_fit as E
    return {'S': E.S_WIRINGS, 'T': E.T_WIRINGS}.get(arg) or tuple(x for x in arg.split(',') if x)


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json: the script, the commits (H.code_commits()), the model's provenance with the lab's input
    (H.model_provenance) for the simulation parts, the settings and `extra`; a part run again adds a dated record."""
    path = out / f'provenance_{part}.json'
    old = json.loads(path.read_text()) if path.exists() else {}
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), **H.code_commits(), **extra)
    if part.startswith('C'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    _write_json(path, dict(script='experiments/E89_apl_breadth_weighting.py', part=part,
                           settings=dict(SEEDS=SEEDS, N_PSEUDO=N_PSEUDO, N_TRIALS=N_TRIALS, EST_TOL=EST_TOL, LEVELS=LEVELS,
                                         BINDING_LEVELS=BINDING_LEVELS, SHARE_SCALE=SHARE_SCALE, LC_EDGES=LC_EDGES,
                                         weightings=[wtag(w) for w in weightings()], C2_UNITS=C2_UNITS),
                           runs=runs))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B', 'C0', 'C1', 'C2', 'Cruns', 'C3', 'Cscore', 'variants', 'table'), required=True)
    ap.add_argument('--workers', type=int, default=4, help='A: processes; C2 / C3: one-odor workers (16)')
    ap.add_argument('--wirings', default='S', help="C0 / C1 / C2: 'S', 'T' or a comma-separated list")
    args = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = H.results_dir(EXPERIMENT)
    ws = _wirings(args.wirings)
    if args.part in ('C0', 'C1', 'C2', 'Cruns', 'C3'):
        provenance(out, args.part, wirings=list(ws), workers=args.workers)
    if args.part == 'A':
        part_a(out, args.workers)
    elif args.part == 'B':
        part_b(out)
    elif args.part == 'C0':
        part_c0(out, 'C1', ws)
        part_c0(out, 'C2', ws, args.workers)
    elif args.part == 'C1':
        part_c0(out, 'C1', ws)
        part_c1(out, ws)
    elif args.part == 'C2':
        part_c2(out, ws, n_workers=args.workers)
    elif args.part == 'Cruns':
        part_c0(out, 'C1', ws)
        part_c1(out, ws)
        part_c2(out, ws, n_workers=args.workers)
    elif args.part == 'C3':
        part_c3(out, args.workers)
    elif args.part == 'Cscore':
        score_c(out)
    elif args.part == 'variants':
        write_variants(out)
    if args.part in ('A', 'B', 'Cscore', 'C3', 'variants', 'table'):
        comparison_table(out)


if __name__ == '__main__':
    main()
