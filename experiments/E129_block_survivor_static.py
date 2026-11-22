#!/usr/bin/env python3
"""E129: the inhibition that outlasts the block reduced to a static level: in the blocked condition only, the family's
inhibition among KCs replaced by its time average, a constant per KC, a constant per KC type, or one constant for all
KCs, against a pooled neuron and a slow matrix.

The entry is item E129 of experiment_priorities.md (proposed by hal 2026-11-22 13:29 at Tom's request of 12:57, moved
into the list under the idle-host rule; planned by a Fable subagent, claude-fable-5-2, spawned from hal 13:05). "The
family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S
is real, P1, P2, M1 and T is P3, P4, M2, M3. "KC-KC at zero" is E110's part A (KC-KC 0 in both conditions, both response
rates held at fixed scales).

The protocol (design item 1; nothing fitted). On each wiring the APL-on run is the family's saved run, untouched, and
the common somatic threshold (fixed_thr) and every other parameter stay the fit's. Only the APL-off condition is
simulated again (`OffRunner`, one of E83's simulators of the wiring at its fit, APL->KC and APL->bouton 0): the KC-KC
matrix is replaced by the arm's inhibition, with one scale bisected (E69's `bisect_factor`, from factor 1) to hold the
family's APL-off response count (`target_off`, the response rate 0.150). Any change of the per-pair APL effect then comes
from the distances with APL off.

The arms (tags S1_<arm>_<wiring> for steps 1a and 1b, S1c_Z_<wiring>, S2_<arm>_<wiring>):
- T (step 1a), the time course removed: a static inhibition (olfsysm's `ext_inh`, KCs x odors, subtracted from each
  KC's input while the odor is on, where KC-KC enters and in the same units, `E95_static_apl.set_ext_inh_ko`) equal to
  each KC's and odor's time-mean of the family's own KC-KC input with APL off, and the KC-KC scale 0.
- G (1a), one level: the same constant for every KC and odor.
- P (1a), a pooled neuron: E122's uniform neuron U (every entry total / N^2, the diagonal included, the hemibrain
  matrix's total weight) in the matrix's place, filtered with 50 ms as the matrix is.
- S (1a), slow: the hemibrain matrix filtered with 1 s (olfsysm's `kc_kc_tau`, set on the APL-off simulator only).
- K (1b, gated), a constant per KC: its mean over the 16 odors of T's inhibition. Y (1b), a constant per KC type: the
  mean over the type's KCs and the odors. KP (1b, gated again), K's constants permuted among the KCs of each type.
- Z (1c, gated), G's run with the level removed: the APL-off run at the fit's threshold with no inhibition among KCs and
  no static inhibition, not bisected (the acute reading).
- C50, C0 (step 2, gated; APL-on runs): the hemibrain matrix at 0.5 and 0 of the fit's KC-KC scale with APL on, the tied
  APL scales (APL->KC, KC->APL, bouton->APL) bisected to the family's APL-on count, the threshold the fit's; each paired
  with the APL-off held run of G and of Y where that arm at least keeps the gain (kinds 'C50|G', 'C0|G', 'C50|Y',
  'C0|Y').
- MB (step 3): MB-C1's counted synapses in the matrix's place with APL off. Refused without --tom-approved-mbc1 (Tom's
  word on Pratyush's hemibrain MB-C1 tables, which the loader does not read); not implemented until then. Its reading,
  fixed now: one neuron, its input from each KC and its output to each KC MB-C1's synapse counts summed over the KC's
  claws, so the matrix in HB's place is the outer product (output x input) scaled to HB's total weight, filtered with
  50 ms, its scale bisected as P's.
- References, read and not run: the family; KC-KC at zero (E110's part A, 'E110_A'); the step-0 preview 'Gpre'.

Readings the entry leaves open, fixed here before any run:
1. T's time-mean. olfsysm filters each KC's spikes (one vesicle per spike: nves 1 while `mp.kc.ves_p` is 0, checked)
   into an activity whose integral over the trial is the KC's spike count, so the family's KC-KC input to KC i for odor o
   has the trial integral s x sum_j HB[i, j] n[j, o] (s the fit's KC-KC scale, HB the hemibrain matrix as the simulator
   holds it, checked equal to the loader's, n the family's saved APL-off spike counts, counted over the whole simulated
   trial). T's inhibition is that integral divided by the pulse length (0.5 s, the simulator's stim window, checked):
   the family's total KC-KC input delivered evenly while the odor is on (`time_mean_input`). It is the time-mean over
   the pulse wherever the spikes and their 50 ms tails fall inside it, and conserves the trial's total otherwise.
2. The static arms' base levels at factor 1 all have T's grand mean: G T's mean over KCs and odors, K each KC's mean
   over the odors, Y the mean over the KCs of each type and the odors, KP K's values permuted within each type
   (numpy's default_rng(KP_SEED), one permutation per type in TYPES order) (`static_base`). Their factors are therefore
   comparable, and each arm's held level is its factor times its base.
3. The one scale. Static arms: the factor multiplies the base matrix (ext_inh_scale stays 1), with the KC-KC scale 0.
   P and S: the factor multiplies the APL-off simulator's KC-KC scale (s x f). Every arm: E69's `bisect_factor` from 1,
   bracket SCALE_BRACKET (1/50 to 50), first step SCALE_STEP (1.25), to `target_off` (two or more spikes). Step 2's APL
   factor: E83's APL_BRACKET (0.1 to 10) and APL_STEP (1.25), to `target_on`.
4. A corner: the scale at a bound of its bracket, or the response rate not held, i.e. the count more than RATE_TOL
   (0.2 %) of the target away from it (about eight responses, a rate within 0.0003). An arm with corners on more than
   CORNER_WIRINGS_MAX (two) wirings gets no class.
5. P's matrix is E122's U; S's 1 s is set as `mp.kc.kc_kc_tau` on the APL-off simulator after it is built (the
   simulator does not reset it per call; checked before every run), not through E10_joint_fit.py.
6. The classes (E123's and E126's, unchanged; calibrated pattern part, S's and T's populations apart, against all
   recordings and each of the four halves A to D):
   - "equal" (`E123.classify`, one-sided): not higher than the family by more than WITHIN (0.003) on both populations
     and on 6 of 8 wirings, on every target, with every run in the four bands;
   - "keeps the gain" (`E123.classify`): at or below the half-way point between the family and KC-KC at zero on both
     populations, every target (0.0554 on S and 0.0577 on T against all recordings, checked in verdict.json);
   - "neither" otherwise; "corner (no class)" by reading 4;
   - two arms "differ" only by E126's rule (`differ`: one lower on 7 of 8 wirings against all recordings and on both
     populations against all five targets), else "not told apart".
   "At least keeps the gain" is the class 'equal' or 'keeps the gain'. "T and G fall in different classes" (step 1b's
   gate) and "K and Y fall in different classes" (KP's) compare the class strings, a corner's included.
7. Step 1c's gate: G at least keeps the gain. Step 2's gate: G or Y at least keeps the gain; its 16 APL-on runs are the
   same whichever static arm they are paired with, so they are run once and paired with each arm that opens it.
8. The step-0 preview 'Gpre': the family's saved APL-on runs paired with E110 part A's APL-off runs (KC-KC 0, the
   threshold bisected to the same count): a uniform level while the odor is on is the threshold raised while the odor
   is on, up to Vm's lag at the step. Scored and classed as an arm (its class is a preview, not G's); score reports G
   minus Gpre per population once G is run.
9. The split (step 0): T's base matrix (KCs x odors, Hz) per wiring, over all KCs and over the KCs that respond with APL
   off in the family: the sum of squares about the grand mean in shares per KC type, per KC within its type, per odor
   and the rest (adding to 1; 'share_kc' the first two together, the static per-KC part), the per type x odor part of
   the rest beside, the CV over the odors at the median KC, and the mean level per type in Hz and against fixed_thr
   (`split_kckc`); beside it, E110 part A's threshold raise (the level G should need) over T's grand mean (G's expected
   factor).
10. Realism of the static arms: the held level (mean, per type, largest) in Hz and against each KC's threshold as the
    simulator sets it (fixed_thr + 2 x spontaneous input, read from olfsysm) and against fixed_thr; the share of KC x
    odor entries whose level is above the KC's threshold ('level_above_thr_share'), flagged as hard to defend when
    positive ('level_above_thr').

Parts (CLI `--part`):
- step0 (no simulation): the preview and the split; step0_preview.csv, step0_runs.csv, step0_preview.json,
  step0_split.csv, provenance_step0.json.
- identity (simulations, `--wirings`): on every wiring, (i) the hemibrain matrix at scale 1 through `OffRunner` gives
  the saved APL-off run exactly, (ii) with a static inhibition of zeros set as well, also exactly; on 'real', (iii) T's
  inhibition at factor 1 with KC-KC 0 changes the counts (the static inhibition is read) and (iv) S's runner at factor
  1 changes them (the 1 s filter is read). identity_<host>_<lane>.json; step 1a refuses to run until every wiring has
  passed (`identity_passed`).
- step1a, step1b, step1c, step2 (simulations, resumable, per lane `--wirings`; a point is claimed by mkdir of
  results/E129_block_survivor_static/claims/<tag>/ and skipped when done or claimed, `--reclaim` taking over a claim left
  by a process that died): rows in runs_<step>_<host>_<lane>.csv, provenance (`H.model_provenance` per point) in
  provenance_<step>_<host>_<lane>.json; runs saved locally (outputs/E129_block_survivor_static/<step>/) and in the
  model-output archive (data/model_outputs/E129_block_survivor_static/<step>/), static inhibitions beside
  (<tag>_ext_inh.parquet). Steps 1b, 1c and 2 refuse unless verdict.json (from `score`) opens their gate.
- step3: refused without --tom-approved-mbc1.
- score (no simulation): readings.csv, runs.csv, verdict.json (classes, gates, comparisons, outcomes), comparison_table.md,
  provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- plan: the points a lane would run.

The outcomes (flags for the write-up, which reads them): (a) G equal; (b) Y or K equal and G not; (c) T equal and K not;
(d) T not at least keeping the gain while P or S does; (e) S against P by the "differ" rule; (f) corners, or no arm
told apart from another. No correlation chooses anything, no pair is weighted, and no threshold moves after a run is
read.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): only the model's hemibrain wiring and E72's seven re-draws of
it, and the model loader's hemibrain KC->KC counts and KC types (step 3's MB-C1 tables wait on Tom's word); no olfsysm
change; nothing is published. The simulation parts import no module another session may be editing (E122 to E126 are
read only by score, step0 and variants).

Usage (from the repository root; 16 threads, under nice; one lane per population):
    ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part identity
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step1a --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E129_block_survivor_static.py --part step1a --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E129_block_survivor_static.py --part score
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import socket
import time
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E129_block_survivor_static'
SCRIPT = 'experiments/E129_block_survivor_static.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
POP_WIRINGS = {'S4': S_WIRINGS, 'T4': T_WIRINGS, 'all8': WIRINGS}
POPS = ('S4', 'T4')
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
STATIC_ARMS = ('T', 'G', 'K', 'Y', 'KP')
STEPS = {'step1a': ('G', 'T', 'P', 'S'), 'step1b': ('K', 'Y', 'KP'), 'step1c': ('Z',), 'step2': ('C50', 'C0')}
TAG_PREFIX = {'step1a': 'S1', 'step1b': 'S1', 'step1c': 'S1c', 'step2': 'S2'}
STEP_OF = {a: s for s, arms in STEPS.items() for a in arms}
ON_DOSE = {'C50': 0.5, 'C0': 0.0}
PAIR_WITH = ('G', 'Y')
ARM_LABELS = {
    'T': "each KC's and odor's time-mean of the family's KC-KC input with APL off, as a static inhibition; KC-KC 0",
    'G': 'one static level for every KC and odor; KC-KC 0',
    'K': "a static level per KC (T's mean over the odors); KC-KC 0",
    'Y': "a static level per KC type (T's mean over the type's KCs and the odors); KC-KC 0",
    'KP': "K's levels permuted among the KCs of each type; KC-KC 0",
    'P': "E122's uniform neuron in the matrix's place with APL off, 50 ms",
    'S': 'the hemibrain matrix filtered with 1 s with APL off',
    'Z': 'no inhibition among KCs and no static level with APL off, not bisected (the acute reading)',
    'C50': "APL on with the hemibrain matrix at 0.5 of the fit's KC-KC scale, the tied APL scales re-held",
    'C0': 'APL on with KC-KC 0, the tied APL scales re-held',
    'Gpre': "preview: the family's APL-on runs with E110 part A's APL-off runs (KC-KC 0, threshold re-held)",
}
PULSE_S = 0.5
SLOW_TAU_S = 1.0
KP_SEED = 0
SCALE_BRACKET, SCALE_STEP = (1 / 50, 50.0), 1.25
RATE_TOL = 0.002
WITHIN = 0.003
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7
HALF_WAY_ENTRY = {'S4': 0.0554, 'T4': 0.0577}
HALF_WAY_TOL = 5e-5
PREVIEW_TOL = 0.001
REFS = ('family', 'E110_A', 'Gpre')
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS
COMPARISONS = (('G', 'T'), ('P', 'S'), ('T', 'P'), ('T', 'S'), ('G', 'P'), ('G', 'S'), ('K', 'Y'), ('K', 'KP'),
               ('Y', 'G'), ('K', 'T'))
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E129_block_survivor_static/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=str) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def lane_name(wirings: Sequence[str]) -> str:
    return '-'.join(wirings)


def n_resp(responses: pd.DataFrame) -> int:
    return int((responses.to_numpy(dtype=float) > 0.5).sum())


def tag_of(arm: str, w: str) -> str:
    return f'{TAG_PREFIX[STEP_OF[arm]]}_{arm}_{w}'


def kc_types_of(kc_index: pd.Index) -> np.ndarray:
    """The KC type ('g', "a'b'", 'ab') of each model KC (the index's kc_type level)."""
    return np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)


def hemibrain_matrix(kc_index: pd.Index) -> pd.DataFrame:
    """The model loader's normalized hemibrain KC x KC matrix (post x pre, mean row sum 1) for the model's KCs, in their
    order (pyolfsysm's `connectome_KCKC_weights`, as fit_mb_model builds it)."""
    from pyolfsysm.mb_model import connectome_KCKC_weights
    ids = pd.Index(kc_index.get_level_values(H.KC_ID_LEVEL))
    M = connectome_KCKC_weights(ids)
    if not (np.array_equal(M.index.to_numpy(), ids.to_numpy()) and np.array_equal(M.columns.to_numpy(), ids.to_numpy())):
        raise RuntimeError("the hemibrain KC-KC matrix is not in the model KCs' order")
    return M


def uniform_neuron(ref: np.ndarray) -> np.ndarray:
    """E122's U: every entry total / N^2 (the diagonal included), `ref`'s total weight."""
    n = ref.shape[0]
    return np.full((n, n), float(ref.sum()) / n / n)


def sha256_array(a: np.ndarray) -> str:
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=float), dtype='<f8').tobytes()).hexdigest()


# ---------------------------------------------------------------------------
# the static inhibitions

def time_mean_input(counts: np.ndarray, M: np.ndarray, scale: float, pulse_s: float = PULSE_S) -> np.ndarray:
    """Reading 1: scale x M @ counts / pulse_s (KCs x odors, Hz), M post x pre, counts KCs x odors (a spike's filtered
    activity integrates to one vesicle)."""
    counts, M = np.asarray(counts, dtype=float), np.asarray(M, dtype=float)
    if M.shape != (counts.shape[0], counts.shape[0]):
        raise ValueError(f'matrix {M.shape} for {counts.shape[0]} KCs')
    if not (np.isfinite(scale) and scale >= 0 and pulse_s > 0):
        raise ValueError(f'scale {scale}, pulse {pulse_s}')
    return float(scale) * (M @ counts) / float(pulse_s)


def family_time_mean(w: str) -> pd.DataFrame:
    """T's base on wiring `w` (reading 1): KCs x odors in the saved outputs' order, from the family's saved APL-off spike
    counts, the loader's hemibrain matrix and the fit's KC-KC scale."""
    F = E83.saved_fit(w)
    C = F['saved']['off'][1]
    M = hemibrain_matrix(C.index)
    return pd.DataFrame(time_mean_input(C.to_numpy(dtype=float), M.to_numpy(dtype=float), float(F['p']['wKCKC_scale'])),
                        index=C.index, columns=C.columns)


def static_base(arm: str, Tm: np.ndarray, types: np.ndarray, seed: int = KP_SEED) -> np.ndarray:
    """Reading 2: the base matrix (KCs x odors) of a static arm from T's `Tm`: 'T' Tm, 'G' its grand mean everywhere,
    'K' each KC's mean over the odors, 'Y' each type's mean over its KCs and the odors, 'KP' K's values permuted among
    the KCs of each type (default_rng(seed), types in TYPES order)."""
    Tm = np.asarray(Tm, dtype=float)
    types = np.asarray(types, dtype=object)
    if types.shape != (Tm.shape[0],):
        raise ValueError(f'{types.shape} types for {Tm.shape[0]} KCs')
    n_od = Tm.shape[1]
    if arm == 'T':
        return Tm.copy()
    if arm == 'G':
        return np.full(Tm.shape, Tm.mean())
    k = Tm.mean(axis=1)
    if arm == 'K':
        return np.repeat(k[:, None], n_od, axis=1)
    if arm == 'Y':
        v = np.zeros(len(k))
        for t in TYPES:
            m = types == t
            if m.any():
                v[m] = Tm[m].mean()
        return np.repeat(v[:, None], n_od, axis=1)
    if arm == 'KP':
        rng = np.random.default_rng(int(seed))
        v = k.copy()
        for t in TYPES:
            idx = np.flatnonzero(types == t)
            v[idx] = k[idx[rng.permutation(idx.size)]]
        return np.repeat(v[:, None], n_od, axis=1)
    raise ValueError(f'not a static arm: {arm!r}')


def split_kckc(A: np.ndarray, types: np.ndarray) -> dict:
    """Reading 9: the sum of squares of a KCs x odors matrix about its grand mean in shares per KC type, per KC within
    its type, per odor and the rest (adding to 1), 'share_kc' the per-KC part (type and within type together), the
    per type x odor part of the rest ('share_type_x_odor', inside 'share_rest'), the CV over the odors at the median KC
    (KCs with a positive mean), and the grand mean and each type's mean."""
    A = np.asarray(A, dtype=float)
    types = np.asarray(types, dtype=object)
    gm = A.mean()
    row_m = A.mean(axis=1)
    col = A.mean(axis=0) - gm
    type_m = np.zeros(len(row_m))
    # the type x odor interaction: each type's odor profile minus the type's mean and the odor's overall effect (the
    # projection of the rest on the type x odor cells)
    tx = np.zeros_like(A)
    out = dict(grand_mean=float(gm))
    for t in TYPES:
        m = types == t
        if m.any():
            type_m[m] = row_m[m].mean()
            tx[m] = (A[m].mean(axis=0) - A[m].mean() - col)[None, :]
        out[f'mean_{TYPE_COL[t]}'] = float(A[m].mean()) if m.any() else np.nan
        out[f'n_kcs_{TYPE_COL[t]}'] = int(m.sum())
    typ = type_m - gm
    within = row_m - type_m
    rest = A - gm - typ[:, None] - within[:, None] - col[None, :]
    tot = float(((A - gm) ** 2).sum())
    n_kc, n_od = A.shape
    mu = row_m
    ok = mu > 0
    cv = A[ok].std(axis=1) / mu[ok]
    if not tot > 0:
        out.update(share_type=np.nan, share_kc_within_type=np.nan, share_kc=np.nan, share_odor=np.nan, share_rest=np.nan,
                   share_type_x_odor=np.nan, cv_odors_median_kc=np.nan)
        return out
    s_type = float((typ ** 2).sum() * n_od / tot)
    s_within = float((within ** 2).sum() * n_od / tot)
    out.update(share_type=s_type, share_kc_within_type=s_within, share_kc=s_type + s_within,
               share_odor=float((col ** 2).sum() * n_kc / tot), share_rest=float((rest ** 2).sum() / tot),
               share_type_x_odor=float((tx ** 2).sum() / tot),
               cv_odors_median_kc=float(np.median(cv)) if cv.size else np.nan)
    return out


def level_realism(M: Optional[np.ndarray], thr: np.ndarray, types: np.ndarray, fixed_thr: float) -> dict:
    """Reading 10: a held static inhibition (KCs x odors, Hz; None: none) against each KC's threshold `thr` (fixed_thr +
    2 x spontaneous input) and fixed_thr."""
    keys = ('level_mean', 'level_max', 'level_over_thr_p50', 'level_over_thr_p90', 'level_over_thr_max',
            'level_over_fixed_thr_mean', 'level_above_thr_share') + tuple(f'level_mean_{TYPE_COL[t]}' for t in TYPES)
    if M is None:
        return {**{k: np.nan for k in keys}, 'level_above_thr': False}
    A = np.asarray(M, dtype=float)
    thr = np.asarray(thr, dtype=float).ravel()
    if thr.shape != (A.shape[0],):
        raise ValueError(f'{thr.shape} thresholds for {A.shape[0]} KCs')
    r = A / thr[:, None]
    q = np.quantile(r, [0.5, 0.9, 1.0])
    row = dict(level_mean=float(A.mean()), level_max=float(A.max()), level_over_thr_p50=float(q[0]),
               level_over_thr_p90=float(q[1]), level_over_thr_max=float(q[2]),
               level_over_fixed_thr_mean=float(A.mean() / fixed_thr), level_above_thr_share=float((r > 1).mean()))
    for t in TYPES:
        m = np.asarray(types, dtype=object) == t
        row[f'level_mean_{TYPE_COL[t]}'] = float(A[m].mean()) if m.any() else np.nan
    row['level_above_thr'] = bool(row['level_above_thr_share'] > 0)
    return row


def type_rates(o: Mapping[str, tuple]) -> dict:
    """Per KC type and condition present: the response rate (two or more spikes) and the mean spike count per KC and
    odor (the summed response); the off / on ratio of the rates where both are present (E122's keys)."""
    row = {}
    for c in ('on', 'off'):
        if c not in o:
            continue
        r, s = o[c]
        t = kc_types_of(r.index)
        R = r.to_numpy(dtype=float) > 0.5
        S = s.to_numpy(dtype=float)
        for tp in TYPES:
            row[f'rate_{c}_{TYPE_COL[tp]}'] = float(R[t == tp].mean())
            row[f'count_{c}_{TYPE_COL[tp]}'] = float(S[t == tp].mean())
    if 'on' in o and 'off' in o:
        for tp in TYPES:
            on = row[f'rate_on_{TYPE_COL[tp]}']
            row[f'ratio_{TYPE_COL[tp]}'] = row[f'rate_off_{TYPE_COL[tp]}'] / on if on > 0 else np.nan
    return row


def held_ok(n: int, target: int) -> bool:
    """Reading 4: the response count within RATE_TOL of the target."""
    return bool(abs(int(n) - int(target)) <= RATE_TOL * int(target))


# ---------------------------------------------------------------------------
# the runners

def _preload() -> None:
    """Imports every module the simulation parts use before E108's commit check (none another session may be editing)."""
    from experiments import E10_joint_fit  # noqa: F401
    from experiments import E69_kc_firing_ceiling  # noqa: F401
    from experiments import E95_static_apl  # noqa: F401
    from experiments import E108_paa_extrapolation  # noqa: F401
    from experiments import E110_no_kckc  # noqa: F401


class OffRunner:
    """The APL-off condition of one wiring at its fit (the module docstring), on one of E83's simulators: kind 'HB' (the
    hemibrain matrix as the family has it, 50 ms), 'U' (E122's uniform neuron in its place, 50 ms) or 'slow' (the
    hemibrain matrix, `kc_kc_tau` SLOW_TAU_S). `run(kckc_factor, ext)`: APL->KC and APL->bouton 0, the KC-KC scale the
    fit's times `kckc_factor`, the static inhibition `ext` (KCs x odors in the outputs' order; None for none) set for
    that run only; every other parameter the fit's.

    Attributes: w, kind, p, saved, target_off, size, kc_index, types, Tm (T's base, reading 1), tau, matrix_sha256,
    matrix_max_diff (the simulator's KC-KC matrix against the expected one), thr (each KC's threshold, after a run)."""

    KINDS = ('HB', 'U', 'slow')

    def __init__(self, w: str, kind: str = 'HB'):
        from experiments import E10_joint_fit as J
        if kind not in self.KINDS:
            raise ValueError(f'kind {kind!r}')
        _preload()
        P = E83.wiring_parts(w)
        self.w, self.kind = w, kind
        self.p = dict(P['p'])
        self.saved = P['saved']
        self.target_off, self.size = int(P['target_off']), int(P['size'])
        self.kc_index = self.saved['off'][1].index
        self.types = kc_types_of(self.kc_index)
        HB = hemibrain_matrix(self.kc_index)
        ref = HB.to_numpy(dtype=float)
        expected, extra = ref, None
        if kind == 'U':
            expected = uniform_neuron(ref)
            extra = dict(_wKCKC=pd.DataFrame(expected, index=HB.index, columns=HB.columns))
        self.sim, self.info = E83.build_simulator(w, extra_kws=extra)
        cells = self.info['cells']
        self.mp, self.rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        self.order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        Wc = cells['wKCKC_unscaled'].cell_contents
        Wd = Wc.toarray() if hasattr(Wc, 'toarray') else np.asarray(Wc, dtype=float)
        self.matrix_max_diff = float(np.abs(Wd - expected).max())
        if not np.allclose(Wd, expected, rtol=1e-12, atol=1e-15):
            raise RuntimeError(f'{w} {kind}: the simulator\'s KC-KC matrix differs from the expected one by up to '
                               f'{self.matrix_max_diff:.3g}')
        self.matrix_sha256 = sha256_array(expected)
        self.tau = SLOW_TAU_S if kind == 'slow' else float(J.TAU)
        self.mp.kc.kc_kc_tau = self.tau
        if float(self.mp.kc.ves_p) != 0.0:
            raise RuntimeError(f'{w}: KC vesicle depletion is on (ves_p {self.mp.kc.ves_p}); reading 1 assumes one vesicle')
        pulse = float(self.mp.time_stim_end - self.mp.time_stim_start)
        if abs(pulse - PULSE_S) > 1e-12:
            raise RuntimeError(f'{w}: the pulse is {pulse} s, not {PULSE_S}')
        self.Tm = time_mean_input(self.saved['off'][1].to_numpy(dtype=float), ref, float(self.p['wKCKC_scale']))
        self.thr = None

    def params(self, kckc_factor: float) -> Dict[str, float]:
        return dict(self.p, wAPLKC_scale=0.0, wAPLPN_scale=0.0, wKCKC_scale=float(self.p['wKCKC_scale']) * float(kckc_factor))

    def run(self, kckc_factor: float = 1.0, ext: Optional[np.ndarray] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        from experiments.E95_static_apl import set_ext_inh_ko
        if float(self.mp.kc.kc_kc_tau) != self.tau:
            raise RuntimeError(f'{self.w}: kc_kc_tau is {self.mp.kc.kc_kc_tau}, not {self.tau}')
        if ext is not None and np.asarray(ext).shape != (len(self.kc_index), len(self.order)):
            raise ValueError(f'ext_inh {np.asarray(ext).shape}')
        set_ext_inh_ko(self.mp, self.rv, ext, None if ext is None else self.order)
        try:
            out = self.sim(self.params(kckc_factor))
        finally:
            set_ext_inh_ko(self.mp, self.rv, None)
        self.thr = np.asarray(self.rv.kc.thr, dtype=float).ravel().copy()
        return out

    def arm_fn(self, arm: str) -> Tuple[Callable[[float], Tuple[pd.DataFrame, pd.DataFrame]], Optional[np.ndarray]]:
        """(f -> outputs, the static base or None) of an arm on this runner."""
        need = {'P': 'U', 'S': 'slow'}.get(arm, 'HB')
        if need != self.kind:
            raise ValueError(f'arm {arm} needs a {need} runner, not {self.kind}')
        if arm in STATIC_ARMS:
            base = static_base(arm, self.Tm, self.types)
            return (lambda f: self.run(0.0, float(f) * base)), base
        if arm in ('P', 'S'):
            return (lambda f: self.run(float(f), None)), None
        if arm == 'Z':
            return (lambda f: self.run(0.0, None)), None
        raise ValueError(f'arm {arm!r} is not an APL-off arm')


def runner_kind(arm: str) -> str:
    return {'P': 'U', 'S': 'slow', 'C50': 'on', 'C0': 'on'}.get(arm, 'HB')


def held_off(runner: OffRunner, arm: str, label: str) -> dict:
    """One APL-off point (reading 3): the arm's scale bisected to the family's APL-off count (Z: no bisection)."""
    from experiments import E69_kc_firing_ceiling as E69
    t0 = time.time()
    fn, base = runner.arm_fn(arm)
    if arm == 'Z':
        out = fn(1.0)
        f, hist, at_bound = np.nan, [], False
    else:
        f, out, hist, at_bound = E69.bisect_factor(fn, runner.target_off, runner.size, SCALE_BRACKET, SCALE_STEP, label)
    n = n_resp(out[0])
    ext = None if base is None else float(f) * base
    held = held_ok(n, runner.target_off) if arm != 'Z' else None
    return dict(factor=float(f), off=out, n_off=n, target_off=runner.target_off, at_bound=bool(at_bound), held=held,
                corner=bool(at_bound or (held is False)), history=hist, n_sims=max(len(hist), 1), ext=ext,
                params_off=runner.params(0.0 if arm in STATIC_ARMS or arm == 'Z' else f), thr=runner.thr,
                seconds=time.time() - t0)


def held_on(w: str, dose: float, label: str) -> dict:
    """One APL-on point of step 2: the hemibrain matrix at `dose` of the fit's KC-KC scale with APL on, the tied APL
    scales (E83.TIED_ON) bisected to the family's APL-on count, the threshold the fit's."""
    from experiments import E69_kc_firing_ceiling as E69
    _preload()
    t0 = time.time()
    P = E83.wiring_parts(w)
    p = dict(P['p'])
    sim, _ = E83.build_simulator(w)
    q = dict(p, wKCKC_scale=float(p['wKCKC_scale']) * float(dose))

    def fn(x):
        return sim(dict(q, **{k: p[k] * float(x) for k in E83.TIED_ON}))

    x, out, hist, at_bound = E69.bisect_factor(fn, int(P['target_on']), int(P['size']), E83.APL_BRACKET, E83.APL_STEP,
                                               label)
    n = n_resp(out[0])
    held = held_ok(n, int(P['target_on']))
    params = dict(q, **{k: p[k] * float(x) for k in E83.TIED_ON})
    del sim
    gc.collect()
    return dict(factor=float(x), on=out, n_on=n, target_on=int(P['target_on']), at_bound=bool(at_bound), held=held,
                corner=bool(at_bound or not held), history=hist, n_sims=len(hist), params_on=params,
                seconds=time.time() - t0)


# ---------------------------------------------------------------------------
# saving and loading runs

def save_run(step: str, tag: str, o: Mapping[str, tuple], ext: Optional[pd.DataFrame] = None) -> None:
    """The run's conditions (and the static inhibition), locally and in the model-output archive."""
    for d in (local_dir(step), archive_dir() / step):
        d.mkdir(parents=True, exist_ok=True)
        items = [(f'{tag}_{cc}_{k}.parquet', x) for c, cc in (('on', 'on'), ('off', 'acute_off')) if c in o
                 for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1]))]
        if ext is not None:
            items.append((f'{tag}_ext_inh.parquet', ext))
        for name, x in items:
            p = d / name
            tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
            x.to_parquet(tmp)
            os.replace(tmp, p)


def load_run(step: str, tag: str, conds: Sequence[str] = ('on', 'off')) -> Optional[Dict[str, tuple]]:
    for d in (local_dir(step), archive_dir() / step):
        out = {}
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            if c not in conds:
                continue
            ps = [d / f'{tag}_{cc}_{k}.parquet' for k in ('responses', 'spike_counts')]
            if not all(p.exists() for p in ps):
                break
            out[c] = tuple(pd.read_parquet(p) for p in ps)
        if len(out) == len(conds):
            return out
    return None


def runs_table(step: Optional[str] = None) -> pd.DataFrame:
    """Every row of every host and lane (runs_<step>_*.csv; all steps if None), the last row per tag."""
    steps = [step] if step else list(STEPS)
    fs = [pd.read_csv(p) for s in steps for p in sorted(out_dir().glob(f'runs_{s}_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag', 'arm', 'wiring', 'step'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def claim(tag: str, lane: str, reclaim: bool = False) -> bool:
    """Claims a point by mkdir of results/<EXPERIMENT>/claims/<tag>/ (atomic on NFS); False if another process holds
    it (unless `reclaim`)."""
    d = out_dir() / 'claims' / tag
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(tag=tag, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(tag: str) -> None:
    shutil.rmtree(out_dir() / 'claims' / tag, ignore_errors=True)


# ---------------------------------------------------------------------------
# step 0 (no simulation)

def preview_outputs() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """The references: ({(kind, wiring): outputs}, runs rows) for 'family', 'E110_A' and 'Gpre' (reading 8)."""
    from experiments import E110_no_kckc as E110
    outs, rows = {}, []
    ra = E110.runs_table_of('A').set_index('tag')
    for w in WIRINGS:
        F = E83.saved_fit(w)
        a = E110.load_a(w)
        if a is None:
            raise RuntimeError(f"E110 part A's run of {w} is neither local nor in the archive")
        ra_w = ra.loc[f'A_{w}']
        if (n_resp(a['off'][0]), n_resp(a['on'][0])) != (int(ra_w.n_off), int(ra_w.n_on)):
            raise RuntimeError(f"E110 part A's run of {w} has not the counts its table records")
        for kind, o, meta in (
                ('family', F['saved'], dict(F['p'], tag=F['tag'], folder=F['folder'])),
                ('E110_A', a, dict({k: float(ra_w[k]) for k in PARAM_KEYS}, tag=f'A_{w}', folder='E110_no_kckc/A')),
                ('Gpre', {'on': F['saved']['on'], 'off': a['off']},
                 dict(F['p'], fixed_thr_off=float(ra_w['fixed_thr']), tag=f'Gpre_{w}',
                      folder=f"{F['folder']} (on) + E110_no_kckc/A (off)"))):
            outs[(kind, w)] = o
            rows.append(run_meta_row(kind, w, o, **meta, corner=False))
    return outs, pd.DataFrame(rows)


def run_meta_row(kind: str, w: str, o: Mapping[str, tuple], **meta) -> dict:
    from experiments import E110_no_kckc as E110
    r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **type_rates(o))
    b = E110.bands_of(float(r['rate_on']), float(r['rate_off']), float(r['rate_rho_on']), float(r['opk_on_tv']),
                      float(r['opk_off_tv']))
    r.update(b, in_bands=all(b.values()))
    return r


def split_table() -> pd.DataFrame:
    """Reading 9 per wiring and KC set ('all', 'responders' the KCs responding to an odor with APL off in the family),
    with E110 part A's threshold raise and the factor it implies for G."""
    from experiments import E110_no_kckc as E110
    ra = E110.runs_table_of('A').set_index('tag')
    rows = []
    for w in WIRINGS:
        F = E83.saved_fit(w)
        Tm = family_time_mean(w)
        types = kc_types_of(Tm.index)
        A = Tm.to_numpy(dtype=float)
        resp = (F['saved']['off'][0].to_numpy(dtype=float) > 0.5).any(axis=1)
        fthr = float(F['p']['fixed_thr'])
        d_thr = float(ra.loc[f'A_{w}', 'fixed_thr']) - fthr
        for kset, m in (('all', np.ones(len(A), dtype=bool)), ('responders', resp)):
            s = split_kckc(A[m], types[m])
            rows.append(dict(wiring=w, kcs=kset, n_kcs=int(m.sum()), fixed_thr=fthr, kckc_scale=float(F['p']['wKCKC_scale']),
                             level_over_fixed_thr=s['grand_mean'] / fthr, e110a_thr_raise=d_thr,
                             g_factor_expected=d_thr / A.mean(), **s))
    return pd.DataFrame(rows)


def part_step0() -> dict:
    """The preview and the split (no simulation): step0_preview.csv (readings), step0_runs.csv, step0_preview.json,
    step0_split.csv, provenance_step0.json."""
    out = out_dir()
    outs, runs = preview_outputs()
    rd = readings(outs, runs)
    rd.to_csv(out / 'step0_preview.csv', index=False)
    runs.to_csv(out / 'step0_runs.csv', index=False)
    sp = split_table()
    sp.to_csv(out / 'step0_split.csv', index=False)
    cls = classify(rd, runs, 'Gpre')
    v = dict(generated=_now(), host=_host(), half_way=half_way_check(rd), preview_class=cls,
             pattern={k: {u: {t: _v(rd, k, u, f'{t}_pattern') for t in targets()} for u in POPS + ('all8',)} for k in REFS},
             parts={k: {u: {c: _v(rd, k, u, f'all_{c}') for c in ('amplitude', 'shape', 'sd_ratio', 'pearson_effect',
                                                                    'rho_effect')}
                         for u in POPS} for k in REFS},
             mean_effect={k: {ro: _v(rd, k, 'all8', 'mean_effect', ro) for ro in READOUTS} for k in REFS},
             loss_raw={k: _v(rd, k, 'all8', 'all_loss', 'raw') for k in REFS},
             gpre_minus_family={u: _v(rd, 'Gpre', u, 'all_pattern') - _v(rd, 'family', u, 'all_pattern') for u in POPS},
             gpre_minus_e110a={u: _v(rd, 'Gpre', u, 'all_pattern') - _v(rd, 'E110_A', u, 'all_pattern') for u in POPS},
             split_mean={kset: sp[sp['kcs'] == kset].drop(columns=['wiring', 'kcs']).mean().to_dict()
                         for kset in ('all', 'responders')})
    _write_json(out / 'step0_preview.json', v)
    _write_json(out / 'provenance_step0.json', dict(generated=_now(), host=_host(), script=SCRIPT, part='step0',
                                                    commits=H.code_commits(),
                                                    inputs=dict(family={w: list(E83.members()[w]) for w in WIRINGS},
                                                                E110_A='data/model_outputs/E110_no_kckc/A/A_<w>_*')))
    print(json.dumps(dict(pattern={k: {u: round(v['pattern'][k][u]['all'], 4) for u in POPS} for k in REFS},
                          preview_class=cls.get('cls'), split=v['split_mean']['all']), indent=1, default=str), flush=True)
    return v


# ---------------------------------------------------------------------------
# the identities (simulations)

def _ndiff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> Tuple[int, int]:
    return (int((a[1].to_numpy(dtype=float) != b[1].to_numpy(dtype=float)).sum()),
            int((a[0].to_numpy(dtype=float) != b[0].to_numpy(dtype=float)).sum()))


def part_identity(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> dict:
    """The identities (the module docstring): identity_<host>_<lane>.json."""
    from experiments import E108_paa_extrapolation as E108
    _preload()
    lane = lane or lane_name(wirings)
    path = out_dir() / f'identity_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    res = dict(_read_json(path), host=_host(), lane=lane, script=SCRIPT, commits=cc)
    res.setdefault('wirings', {})
    for w in wirings:
        if res['wirings'].get(w, {}).get('passes'):
            print(f'E129 identity {w}: done', flush=True)
            continue
        t0 = time.time()
        R = OffRunner(w, 'HB')
        saved = R.saved['off']
        a = R.run(1.0, None)
        b = R.run(1.0, np.zeros((len(R.kc_index), len(R.order))))
        r = dict(hb_scale1=dict(zip(('counts', 'responses'), _ndiff(a, saved))),
                 zeros_static=dict(zip(('counts', 'responses'), _ndiff(b, saved))),
                 matrix_max_diff=R.matrix_max_diff, matrix_sha256=R.matrix_sha256)
        r['passes'] = all(v == 0 for k in ('hb_scale1', 'zeros_static') for v in r[k].values())
        if w == 'real':
            t = R.run(0.0, R.Tm)
            r['t_factor1'] = dict(n_differing=_ndiff(t, saved)[0], n_responses=n_resp(t[0]), target_off=R.target_off)
            del R
            gc.collect()
            Sl = OffRunner(w, 'slow')
            s = Sl.run(1.0, None)
            r['slow_factor1'] = dict(n_differing=_ndiff(s, saved)[0], n_responses=n_resp(s[0]), target_off=Sl.target_off,
                                     tau=float(Sl.mp.kc.kc_kc_tau))
            del Sl
            r['passes'] = bool(r['passes'] and r['t_factor1']['n_differing'] > 0 and r['slow_factor1']['n_differing'] > 0)
        else:
            del R
        gc.collect()
        r.update(seconds=time.time() - t0, checked=_now())
        res['wirings'][w] = r
        _write_json(path, res)
        print(f'E129 identity {w}: {json.dumps(r, default=str)}', flush=True)
    res['passes'] = all(res['wirings'].get(w, {}).get('passes') for w in wirings)
    _write_json(path, res)
    return res


def identity_passed() -> bool:
    """Every wiring has passed the identities in some host's identity_*.json."""
    ok = set()
    for p in out_dir().glob('identity_*.json'):
        for w, r in _read_json(p).get('wirings', {}).items():
            if r.get('passes'):
                ok.add(w)
    return ok >= set(WIRINGS)


# ---------------------------------------------------------------------------
# steps 1 and 2 (simulations)

def gates_from_verdict() -> dict:
    return _read_json(out_dir() / 'verdict.json').get('gates', {})


def step_arms(step: str, g: Optional[Mapping] = None) -> Tuple[str, ...]:
    """The arms a step runs; step 1b's KP only if its gate is open (`g`, verdict.json's gates)."""
    if step == 'step1b':
        g = gates_from_verdict() if g is None else g
        return ('K', 'Y') + (('KP',) if g.get('kp', {}).get('open') else ())
    return STEPS[step]


def planned(step: str, wirings: Sequence[str] = WIRINGS, arms: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """The (arm, wiring) points of a lane in order: step 1a G and T per wiring (one simulator), then P, then S; the
    other steps arm by arm."""
    arms = tuple(STEPS[step] if arms is None else arms)
    ws = [w for w in WIRINGS if w in wirings]
    if step == 'step1a':
        first = [(a, w) for w in ws for a in ('G', 'T') if a in arms]
        return first + [(a, w) for a in ('P', 'S') if a in arms for w in ws]
    if step == 'step1b':
        return [(a, w) for w in ws for a in ('K', 'Y') if a in arms] + [('KP', w) for w in ws if 'KP' in arms]
    return [(a, w) for a in arms for w in ws]


def step_allowed(step: str) -> Tuple[bool, str]:
    if step == 'step1a':
        return identity_passed(), 'the identities have not passed on every wiring (run --part identity)'
    g = gates_from_verdict()
    key = {'step1b': 'step1b', 'step1c': 'step1c', 'step2': 'step2'}[step]
    return bool(g.get(key, {}).get('open')), f"verdict.json's gate {key}: {g.get(key, {}).get('why', 'no verdict yet')}"


def run_row(step: str, tag: str, arm: str, w: str, res: dict, o: Mapping[str, tuple], cc: dict, lane: str,
            runner: Optional[OffRunner] = None) -> dict:
    P = E83.wiring_parts(w)
    row = dict(step=step, tag=tag, arm=arm, wiring=w, factor=res['factor'], at_bound=res['at_bound'], held=res['held'],
               corner=res['corner'], n_sims=res['n_sims'], seconds=res['seconds'], fit_folder=P['folder'],
               fit_tag=P['tag'], **{k: float(P['p'][k]) for k in PARAM_KEYS}, host=_host(), lane=lane)
    if step == 'step2':
        row.update(n_on=res['n_on'], target_on=res['target_on'], dose_on=ON_DOSE[arm],
                   params_on=json.dumps({k: float(v) for k, v in res['params_on'].items()}),
                   **{f'{k}_on': float(res['params_on'][k]) for k in PARAM_KEYS})
        row.update({k: v for k, v in type_rates(o).items()})
        row['rate_on'] = float(o['on'][0].to_numpy(dtype=float).mean())
        row['max_count_on'] = int(o['on'][1].to_numpy(dtype=float).max())
    else:
        po = res['params_off']
        row.update(n_off=res['n_off'], target_off=res['target_off'], kckc_scale_off=float(po['wKCKC_scale']),
                   kc_kc_tau_off=float(runner.tau), runner_kind=runner.kind, kckc_matrix_sha256=runner.matrix_sha256,
                   params_off=json.dumps({k: float(v) for k, v in po.items()}))
        row.update(run_meta_row(arm, w, o))
        row.pop('kind', None)
        row.update(level_realism(res['ext'], res['thr'], runner.types, float(P['p']['fixed_thr'])))
        row['thr_p50'] = float(np.median(res['thr']))
    row.update(cc)
    return row


def part_runs(step: str, wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
              arms: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Steps 1a, 1b, 1c and 2 (the module docstring), resumable, one point claimed at a time."""
    from experiments import E108_paa_extrapolation as E108
    ok, why = step_allowed(step)
    if not ok:
        raise SystemExit(f'E129 {step} refused: {why}')
    arms = tuple(step_arms(step) if arms is None else arms)
    _preload()
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_{step}_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_{step}_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    runner: Optional[OffRunner] = None
    for arm, w in planned(step, wirings, arms):
        tag = tag_of(arm, w)
        if tag in set(runs_table()['tag']):
            print(f'E129 {tag}: done', flush=True)
            continue
        if not claim(tag, lane, reclaim):
            print(f'E129 {tag}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            P = E83.wiring_parts(w)
            ext_df = None
            if step == 'step2':
                if runner is not None:
                    runner = None
                    gc.collect()
                res = held_on(w, ON_DOSE[arm], f'E129 {tag}')
                o = {'on': res['on']}
                change = (f"arm {arm} ({ARM_LABELS[arm]}): APL on, wKCKC_scale {ON_DOSE[arm]:g} x the fit's, the tied APL "
                          f"scales {E83.TIED_ON} bisected to the family's APL-on count; threshold the fit's")
            else:
                kind = runner_kind(arm)
                if runner is None or runner.w != w or runner.kind != kind:
                    runner = None
                    gc.collect()
                    runner = OffRunner(w, kind)
                res = held_off(runner, arm, f'E129 {tag}')
                o = {'on': P['saved']['on'], 'off': res['off']}
                if res['ext'] is not None:
                    ext_df = pd.DataFrame(res['ext'], index=runner.kc_index, columns=P['saved']['off'][1].columns)
                change = (f"arm {arm} ({ARM_LABELS[arm]}): the APL-on run the family's saved one; with APL off "
                          + ('the static inhibition (ext_inh) factor x its base, KC-KC 0' if arm in STATIC_ARMS else
                             'KC-KC 0 and no static inhibition, not bisected' if arm == 'Z' else
                             f'the KC-KC scale factor x the fit\'s ({runner.kind} matrix, tau {runner.tau:g} s)')
                          + ", bisected to the family's APL-off count; threshold and the rest the fit's")
            save_run(step, tag, o, ext_df)
            row = run_row(step, tag, arm, w, res, o, cc, lane, runner)
            E83._append_row(table, row)
            prov = _read_json(prov_path) or dict(part=step, lane=lane, script=SCRIPT, points={})
            prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
                change=change, factor=res['factor'], fit=f"{P['folder']}:{P['tag']}",
                kckc_matrix_sha256=None if runner is None or step == 'step2' else runner.matrix_sha256,
                kc_kc_tau_off=None if runner is None or step == 'step2' else runner.tau,
                bracket=list(E83.APL_BRACKET if step == 'step2' else SCALE_BRACKET), history=res['history']))
            _write_json(prov_path, prov)
            print(f"E129 {tag}: factor {res['factor']:.4f} ({res.get('n_off', res.get('n_on'))}/"
                  f"{res.get('target_off', res.get('target_on'))}); corner {res['corner']}; bands {row.get('in_bands')}; "
                  f"flag '{row.get('realism_flag') or ''}'; {res['seconds']:.0f}s", flush=True)
            del res, o
            gc.collect()
        except BaseException:
            release(tag)
            raise
    return runs_table(step)


def part_step3(tom_approved_mbc1: bool = False) -> None:
    """Step 3, refused without Tom's word on the MB-C1 tables (the module docstring fixes its reading)."""
    if not tom_approved_mbc1:
        raise SystemExit("E129 step 3 waits on Tom's word on Pratyush's hemibrain MB-C1 tables (--tom-approved-mbc1)")
    raise NotImplementedError('E129 step 3: written once Tom has approved the MB-C1 tables (KCclaw-to-MBC1unit and '
                              'MBC1unit-to-KCclaw of run 2025-12-05); their reading is fixed in the module docstring')


# ---------------------------------------------------------------------------
# score

def targets() -> Tuple[str, ...]:
    from experiments import E108_paa_extrapolation as E108
    return ('all',) + tuple(E108.HALVES)


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _dist(o: Mapping[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83.pair_distances), cached by the sha256 of its spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist' not in _CACHE:
        try:
            _CACHE['dist'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist'] = {}
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in _CACHE['dist']:
        _CACHE['dist'][key] = E83.pair_distances(dict(o), readout)
        _CACHE['dist_dirty'] = True
    return _CACHE['dist'][key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale',
            'factor', 'level_mean', 'level_over_thr_p50', 'level_over_thr_max', 'level_above_thr_share') + tuple(
    f'{k}_{TYPE_COL[t]}' for t in TYPES for k in ('rate_on', 'rate_off', 'ratio', 'count_on', 'count_off'))


def readings(outs: Mapping[Tuple[str, str], Mapping[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E110's
    `scores_of`, the runs' means of RUN_COLS ('run_<col>'), the number of runs outside the bands."""
    from experiments import E110_no_kckc as E110
    rows = []
    for kind in dict.fromkeys(runs['kind']):
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].drop_duplicates('wiring', keep='last').set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in POP_WIRINGS.items() if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                d = dd[unit] if len(ws) == 1 else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, n_wirings=len(ws), readout=ro,
                                 n_outside_bands=int((~s['in_bands'].astype(bool)).sum()),
                                 **{f'run_{c}': float(pd.to_numeric(s[c], errors='coerce').mean())
                                    for c in RUN_COLS if c in s.columns},
                                 **E110.scores_of(d)))
    _save_dist_cache()
    return pd.DataFrame(rows)


def complete(runs: pd.DataFrame, kind: str) -> bool:
    return all(((runs['kind'] == kind) & (runs['wiring'] == w)).any() for w in WIRINGS)


def corners(runs: pd.DataFrame, kind: str) -> dict:
    """Reading 4: the wirings on which a kind's run is a corner, and whether there are more than CORNER_WIRINGS_MAX."""
    sub = runs[runs['kind'] == kind]
    c = sub['corner'].astype(bool) if 'corner' in sub else pd.Series(False, index=sub.index)
    ws = sorted(sub.loc[c, 'wiring'])
    return dict(corner_wirings=ws, corner=len(ws) > CORNER_WIRINGS_MAX)


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """An arm's class (reading 6): E123's ("equal", "keeps the gain", "neither") and the corner rule, with the gain kept
    (E124's share: (E110 A - arm) / (E110 A - family), all recordings) and the arm's run means beside."""
    from experiments import E123_kckc_form as E123
    c = E123.classify(rd, runs, arm)
    if not c.get('complete'):
        return c
    c.update(corners(runs, arm))
    cls = c['cls']
    if c['corner']:
        c['cls_before_corner'] = cls
        cls = 'corner (no class)'
    c['cls'] = cls
    c['at_least_keeps_gain'] = cls in ('equal', 'keeps the gain')
    gain = {}
    for u in POPS + ('all8',):
        f, z, a = (_v(rd, k, u, 'all_pattern') for k in ('family', 'E110_A', arm))
        gain[u] = float((z - a) / (z - f))
    c['gain_kept'] = gain
    sub = runs[runs['kind'] == arm]
    for col in ('factor', 'level_mean', 'level_over_thr_max', 'level_above_thr_share', 'fixed_thr'):
        if col in sub:
            c[f'{col}_mean'] = float(pd.to_numeric(sub[col], errors='coerce').mean())
    if 'level_above_thr' in sub:
        c['level_above_thr_wirings'] = sorted(sub.loc[sub['level_above_thr'].astype(bool), 'wiring'])
    return c


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    """E126's "differ" rule (reading 6) with this experiment's corners: 'verdict' '<a> lower', '<b> lower', 'not told
    apart', 'no class (corner)' or 'incomplete'."""
    if not (complete(runs, a) and complete(runs, b)):
        return dict(a=a, b=b, verdict='incomplete')
    per_w = {w: _v(rd, a, w, 'all_pattern') - _v(rd, b, w, 'all_pattern') for w in WIRINGS}
    per_t = {t: {u: _v(rd, a, u, f'{t}_pattern') - _v(rd, b, u, f'{t}_pattern') for u in POPS} for t in targets()}
    n_a, n_b = sum(v < 0 for v in per_w.values()), sum(v > 0 for v in per_w.values())
    a_all = all(v < 0 for t in per_t.values() for v in t.values())
    b_all = all(v > 0 for t in per_t.values() for v in t.values())
    corner = [x for x in (a, b) if corners(runs, x)['corner']]
    if corner:
        verdict = 'no class (corner)'
    elif n_a >= DIFFER_WIRINGS and a_all:
        verdict = f'{a} lower'
    elif n_b >= DIFFER_WIRINGS and b_all:
        verdict = f'{b} lower'
    else:
        verdict = 'not told apart'
    return dict(a=a, b=b, verdict=verdict, n_wirings_a_lower=int(n_a), n_wirings_b_lower=int(n_b),
                a_lower_on_all_targets=bool(a_all), b_lower_on_all_targets=bool(b_all), corner_arms=corner,
                d_per_wiring=per_w, d_per_target=per_t)


def gates(classes: Mapping[str, dict]) -> dict:
    """The gates of steps 1b (T and G in different classes), KP (K and Y in different classes), 1c (G at least keeps
    the gain) and 2 (G or Y at least keeps the gain; 'pair_with' the arms that do)."""
    def cls(a):
        c = classes.get(a, {})
        return c.get('cls') if c.get('complete') else None

    def keeps(a):
        return bool(classes.get(a, {}).get('at_least_keeps_gain'))

    out = {}
    t, g = cls('T'), cls('G')
    out['step1b'] = dict(open=bool(t and g and t != g), why=f'T {t}, G {g}')
    k, y = cls('K'), cls('Y')
    out['kp'] = dict(open=bool(k and y and k != y), why=f'K {k}, Y {y}')
    out['step1c'] = dict(open=keeps('G'), why=f'G {g}')
    pw = [a for a in PAIR_WITH if keeps(a)]
    out['step2'] = dict(open=bool(pw), pair_with=pw, why=f'G {g}, Y {y}')
    return out


def outcomes(classes: Mapping[str, dict], comps: Mapping[str, dict]) -> dict:
    """The entry's outcomes as flags (the module docstring); None where an arm they need is not run."""
    def done(a):
        return bool(classes.get(a, {}).get('complete'))

    def eq(a):
        return classes.get(a, {}).get('cls') == 'equal'

    def keeps(a):
        return bool(classes.get(a, {}).get('at_least_keeps_gain'))

    arms = [a for a in classes if a not in REFS]
    told = [k for k, c in comps.items() if c.get('verdict', '').endswith('lower')]
    return dict(
        a=eq('G') if done('G') else None,
        b=(bool((eq('Y') or eq('K')) and not eq('G')) if done('G') and (done('Y') or done('K')) else None),
        c=(bool(eq('T') and not eq('K')) if done('T') and done('K') else None),
        d=(bool(not keeps('T') and (keeps('P') or keeps('S'))) if done('T') and (done('P') or done('S')) else None),
        e=comps.get('S|P', {}).get('verdict'),
        f=bool(any(classes[a].get('corner') for a in arms if done(a)) or (bool(comps) and not told)),
        told_apart=told)


def half_way_check(rd: pd.DataFrame) -> dict:
    out = {}
    for u, want in HALF_WAY_ENTRY.items():
        hw = (_v(rd, 'family', u, 'all_pattern') + _v(rd, 'E110_A', u, 'all_pattern')) / 2
        out[u] = dict(half_way=hw, entry=want, same=bool(abs(hw - want) <= HALF_WAY_TOL))
    return out


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run: the references (`preview_outputs`), steps 1a to 1c (the family's APL-on run with each APL-off run), and
    step 2's pairings ('C50|G', ...; the APL-on run with the paired arm's APL-off run)."""
    outs, runs = preview_outputs()
    rows = runs.to_dict('records')
    rt = runs_table()
    off = {}
    for r in rt[rt['step'] != 'step2'].to_dict('records') if len(rt) else []:
        o = load_run(r['step'], r['tag'])
        if o is None:
            continue
        outs[(r['arm'], r['wiring'])] = o
        off[(r['arm'], r['wiring'])] = (o, r)
        rows.append(dict(r, kind=r['arm']))
    for r in rt[rt['step'] == 'step2'].to_dict('records') if len(rt) else []:
        o_on = load_run('step2', r['tag'], ('on',))
        if o_on is None:
            continue
        for a in PAIR_WITH:
            if (a, r['wiring']) not in off:
                continue
            o_off, r_off = off[(a, r['wiring'])]
            kind = f"{r['arm']}|{a}"
            o = {'on': o_on['on'], 'off': o_off['off']}
            outs[(kind, r['wiring'])] = o
            meta = {k: r[k] for k in ('tag', 'step', 'factor', 'at_bound', 'held', 'dose_on', 'host', 'lane')}
            meta.update({k: float(r[f'{k}_on']) for k in PARAM_KEYS}, tag_off=r_off['tag'],
                        corner=bool(r['corner']) or bool(r_off['corner']))
            rows.append(run_meta_row(kind, r['wiring'], o, **meta))
    out = pd.DataFrame(rows)
    for col in ('corner', 'level_above_thr'):
        if col in out:
            out[col] = out[col].fillna(False).astype(bool)
    return outs, out


def comparison_md(rd: pd.DataFrame, runs: pd.DataFrame, kinds: Sequence[str], classes: Mapping[str, dict]) -> str:
    from experiments import E123_kckc_form as E123
    lines = ['# E129 comparison tables (generated)', '',
             'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
             '## The arms (calibrated readout, all recordings)', '',
             '| kind | pattern S4 | pattern T4 | amplitude S4 | shape S4 | spread ratio S4 | Pearson r S4 | ρ effect S4 | '
             'mean effect (readout) | mean effect (raw) | loss (raw) | factor | level (Hz) | level / thr max | corners | '
             'class |', '|---|' + '---:|' * 14 + '---|']
    for k in kinds:
        if not ((rd['kind'] == k) & (rd['unit'] == 'all8')).any():
            continue
        c = classes.get(k, {})
        sub = runs[runs['kind'] == k]
        f = pd.to_numeric(sub.get('factor', pd.Series(dtype=float)), errors='coerce').mean()
        lv = pd.to_numeric(sub.get('level_mean', pd.Series(dtype=float)), errors='coerce').mean()
        lt = pd.to_numeric(sub.get('level_over_thr_max', pd.Series(dtype=float)), errors='coerce').max()
        lines.append(f"| {k} | {_v(rd, k, 'S4', 'all_pattern'):.4f} | {_v(rd, k, 'T4', 'all_pattern'):.4f} | "
                     f"{_v(rd, k, 'S4', 'all_amplitude'):.4f} | {_v(rd, k, 'S4', 'all_shape'):.4f} | "
                     f"{_v(rd, k, 'S4', 'all_sd_ratio'):.3f} | {_v(rd, k, 'S4', 'all_pearson_effect'):.3f} | "
                     f"{_v(rd, k, 'S4', 'all_rho_effect'):.3f} | {_v(rd, k, 'all8', 'mean_effect'):.3f} | "
                     f"{_v(rd, k, 'all8', 'mean_effect', 'raw'):.3f} | {_v(rd, k, 'all8', 'all_loss', 'raw'):.4f} | "
                     f"{f:.3f} | {lv:.1f} | {lt:.2f} | {len(corners(runs, k)['corner_wirings'])} | {c.get('cls', '')} |")
    lines.append('')
    for u in ('S4', 'T4', 'all8'):
        lines += [f'## {u}', '', E123.comparison_table(rd, list(kinds), u), '']
    lines += ['## Runs and bands', '', E123.runs_bands_table(runs), '']
    return '\n'.join(lines)


def acute_table(runs: pd.DataFrame) -> List[dict]:
    """Step 1c's reading per wiring: Z (no level) against G's held run, the response rate with APL off and the summed
    response (mean spike count per KC and odor) by KC type."""
    out = []
    for w in WIRINGS:
        z, g = (runs[(runs['kind'] == k) & (runs['wiring'] == w)] for k in ('Z', 'G'))
        if not (len(z) and len(g)):
            continue
        z, g = z.iloc[-1], g.iloc[-1]
        r = dict(wiring=w, rate_off_Z=float(z['rate_off']), rate_off_G=float(g['rate_off']))
        for t in TYPES:
            c = TYPE_COL[t]
            r.update({f'rate_off_{c}_Z': float(z[f'rate_off_{c}']), f'rate_off_{c}_G': float(g[f'rate_off_{c}']),
                      f'count_off_{c}_Z': float(z[f'count_off_{c}']), f'count_off_{c}_G': float(g[f'count_off_{c}'])})
        out.append(r)
    return out


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    kinds = list(dict.fromkeys(runs['kind']))
    arms = [k for k in kinds if k not in ('family', 'E110_A')]
    classes = {a: classify(rd, runs, a) for a in arms}
    comps = {f'{a}|{b}': differ(rd, runs, a, b) for a, b in COMPARISONS + (('S', 'P'),) if a in kinds and b in kinds}
    g = gates(classes)
    g_minus = {u: _v(rd, 'G', u, 'all_pattern') - _v(rd, 'Gpre', u, 'all_pattern') for u in POPS} if 'G' in kinds else None
    v = dict(generated=_now(), host=_host(), identity_passed=identity_passed(), half_way=half_way_check(rd),
             classes=classes, comparisons=comps, gates=g, outcomes=outcomes(classes, comps),
             g_minus_gpre=g_minus, g_matches_preview=(None if g_minus is None else
                                                      bool(all(abs(x) <= PREVIEW_TOL for x in g_minus.values()))),
             acute=acute_table(runs),
             rules=dict(within=WITHIN, corner_wirings_max=CORNER_WIRINGS_MAX, differ_wirings=DIFFER_WIRINGS,
                        rate_tol=RATE_TOL, scale_bracket=list(SCALE_BRACKET), scale_step=SCALE_STEP,
                        targets=list(targets()), populations=list(POPS)),
             n_runs=int(len(runs)))
    _write_json(out / 'verdict.json', v)
    (out / 'comparison_table.md').write_text(comparison_md(rd, runs, kinds, classes))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT, commits=H.code_commits()))
    print(json.dumps(dict(classes={a: c.get('cls') for a, c in classes.items()}, gates=g,
                          outcomes=v['outcomes']), default=str), flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every E129 run and pairing, absolute parameters) and model_pair_distances_long.csv (raw spike
    counts) for the archive."""
    out = out_dir()
    outs, runs = collect()
    rows, parts = [], []
    mine = runs[~runs['kind'].isin(REFS)]
    for r in mine.itertuples():
        o = outs[(r.kind, r.wiring)]
        arm = r.kind.split('|')[0]
        rows.append(dict(tag=r.tag if '|' not in r.kind else f'{r.tag}|{getattr(r, "tag_off")}', variant=f'{r.kind}_{r.wiring}',
                         arm=r.kind, arm_label=ARM_LABELS.get(arm, arm), step=STEP_OF.get(arm), base='bouton', k=2.0,
                         wiring=r.wiring, ref_experiment=E83.members()[r.wiring][0], ref_fit=E83.members()[r.wiring][1],
                         factor=r.factor, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         kc_kc_inh=r.wKCKC_scale, kckc_scale_off=getattr(r, 'kckc_scale_off', np.nan),
                         kc_kc_tau_off=getattr(r, 'kc_kc_tau_off', np.nan), level_mean=getattr(r, 'level_mean', np.nan),
                         rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag, corner=bool(r.corner),
                         subdir=STEP_OF.get(arm)))
        parts.append(pd.concat({(rows[-1]['variant'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (rows[-1]['variant'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'identity') + tuple(STEPS) + ('step3', 'score', 'variants',
                                                                                             'plan'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    ap.add_argument('--lane', default=None, help='the lane name in the tables (default the wirings joined by -)')
    ap.add_argument('--reclaim', action='store_true', help='take over claims left by a process that died')
    ap.add_argument('--tom-approved-mbc1', action='store_true',
                    help="(step3) Tom has approved Pratyush's hemibrain MB-C1 tables for this experiment")
    a = ap.parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'step0':
        part_step0()
    elif a.part == 'identity':
        part_identity(ws, a.lane)
    elif a.part in STEPS:
        part_runs(a.part, ws, a.lane, a.reclaim)
    elif a.part == 'step3':
        part_step3(a.tom_approved_mbc1)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'plan':
        for step in STEPS:
            for arm, w in planned(step, ws):
                print(step, tag_of(arm, w))


if __name__ == '__main__':
    main()
