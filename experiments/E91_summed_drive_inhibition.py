#!/usr/bin/env python3
"""E91: APL's inhibition on the KC's summed drive instead of claw by claw.

The entry is item E91 of experiment_priorities.md (proposed by atlas 2026-11-13 11:38; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-13 11:00). In olfsysm each claw's net drive x (its PN input minus APL's
inhibition of the claw) is clipped at 0 before the claws are summed. olfsysm's `p.kc.claw_inh_pool` (s, the pooled
share; `fit_mb_model(claw_inh_pool=...)`) makes a KC's summed drive the sum over its claws of max(x, 0) + s * min(x, 0):
s = 0 is every model so far, s = 1 is APL's inhibition on the summed drive. The saved claw states and APL's drive from
the claws stay the clipped ones, so with APL off nothing changes, and at fixed scales APL's activity and the boutons are
those of s = 0: only each KC's summed drive loses s times U, the inhibition its claws do not absorb (U = -sum min(x, 0)).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's
(`recording_halves()`, seed 0), C / D the second split's (seed 1; E89's `targets89`).

Parts (CLI `--part`):
- 0, the identities: s = 0 through the new path (the parameter set on each S wiring's simulator, and on the real wiring a
  simulator built with `claw_inh_pool=0.0`) reproduces the saved fits with no count differences; at s = 1 the APL-off
  runs equal the saved ones; the parameter set on a simulator equals the simulator built with it (s = 1, APL on, the
  real wiring); E89's C1 point F at exponent 1 on P1 reproduced exactly (APL scale and counts). identity.json.
- A (no new parameter): on S at the fits (APL on, one run a wiring with the bouton rates and APL's activity saved), per
  KC and odor the inhibition its claws do not absorb (U, window mean over the odor pulse), relative to its summed drive D
  and to its margin over threshold (the share of the margin it takes: (peak0 - peak1) / (peak0 - thr), peaks of the
  KC's leaky integration of D and of D - U, no spikes); by class (the responders with APL on; E88's classes of the
  ester and the alcohol groups' shared responders; number of claws; quartile of the breadth weight; KC type); the
  prediction of which responders s = 1 removes at fixed scales (a count-level reading that ignores the spikes' resets
  and the KC-KC inhibition: each KC's count times the ratio of its areas above threshold of the leaky integration of
  D - U and of D, rounded; a responder predicted below two spikes is removed), read as the change of the per-pair APL
  effect with the removed responders' counts set to 0 ('removed', the entry's check) and with every count scaled
  ('scaled', beside). unabsorbed.csv; per KC and odor under outputs/E91_summed_drive_inhibition/A/. Acheck (after
  B1): the predicted against the simulated change of the per-pair effect (B1's s = 1 at F exponent 2 re-held, the
  entry's check; its run at the fit's scales beside); prediction.csv, prediction.json.
- B1: the fits as they are, s of 0.5 and 1 at the family's measure F with exponents 2, 1 and 0.5 and with no weighting
  (8 settings, 32 points), the tied APL scale bisected to each wiring's APL-on count (E89's `weighting_reading`); plus
  s = 1 at F exponent 2 at the fit's scales (no bisection; 4 points, for part A's check, not a setting of the choice).
  E89's saved C1 points are the references at s = 0 (variants B1_F_k2 = the unchanged fits, B1_F_k1, B1_F_k0.5,
  B1_none).
- B2: the base with the flies' recruitment (E88's per-odor thresholds; each odor's APL scale bisected to the control
  flies' count for that odor, as E89's `part_c2`), s = 1 at exponents 2 and 1 and with no weighting (12 points). E89's
  saved C2 points are the references (B2_F_k2, B2_F_k1, B2_none); first, on the real wiring, this path at s = 0 must
  reproduce E89's C2 run at exponent 2 exactly (identity.json).
- B3, for the record: olfsysm's `allow_net_inh_per_claw=True` at exponent 2 on the fits as they are (4 points).
- co: the co-response change of the four nearest ester pairs, the ester and the alcohol groups under APL beside every S
  point (E89's `coresponse_c`: readings 1 and 1b, and the geometric mean of E88's two estimates, 'gm', at L1 and L3).
- C: the choice per base and half (`choose91`: the lowest pattern part of S's population under the calibrated readout
  among the settings inside the bands and not at a bracket's bound); a choice with s above 0 that gains CHOICE_MIN_GAIN
  or more over its base's reference (B1: the unchanged fits; B2: the base at exponent 2, s = 0) is run on T and read
  with E83's gate against the other half, on both splits (`gate91`, four directions); the control at s = 0 on S and T
  (`control_matrix`: a static inhibition at the soma, olfsysm's `ext_inh_kc`, the same for every odor and equal to each
  KC's mean over the odors of s times its unabsorbed inhibition at the chosen setting; the APL scale re-held); then
  `score`.
- score: readings.csv, coresponse.csv (unless --no-co), choice.csv, gate.csv, equal.csv, verdict.json, figures,
  comparison_table.md.
- D: the family fitted again (`E10_joint_fit.py -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl --apl-readout both
  --tag _cap3 --claw-inh-pool 1`, 150 evaluations, from each wiring's own fit of C2's family (exactly, `--init-from-exact`;
  at exponent 1 from its fitted factors, since an exact start needs the same exponent), one process at a time;
  `--stage S2` on S, `S1` at exponent 1 only if B1's s = 1 point at exponent 1 is an equal or better, `T2` / `T1` on T
  only if S's fitted population is D_S_GAIN or more below the unchanged fits'). `--part Dscore` scores them (d_scores.csv,
  the verdict's 'D').
- E: only on a pass (not written).
- variants: variants.csv and model_pair_distances_long.csv (every run, the references excepted, and D's fits), for the
  archive: `./venv/bin/python archive_tools/archive_model_outputs.py E91_summed_drive_inhibition`, then
  `./venv/bin/python archive_tools/migrate_model_output_params.py --apply`.
- table: comparison_table.md.

Every run row carries the absolute parameters (the tied APL scale's factor and the scales as run), s, the weighting, the
realism (E83's `realism`: the four bands under the rule of 16:10, `realism_flag`, with the earlier rule beside,
`realism_flag_zcrit`), the largest spike counts, the number of claws carrying half of a responder's drive (median over
the responders, APL on; and APL off), the unabsorbed inhibition over the responders, and `H.code_commits()` of the
process that ran it (which must be clean: E89's `_commits_clean`).

Outputs: results/E91_summed_drive_inhibition/; spike counts and responses under outputs/E91_summed_drive_inhibition/
(B/ for the screen and C's runs, A/ for part A's per KC tables; D's fits at the folder's root, as E10_joint_fit.py
writes them), each run's unabsorbed inhibition per KC and odor beside its counts (<tag>_unabsorbed.parquet).

Usage (from the repository root; one 16-thread process at a time on viki):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E91_summed_drive_inhibition.py --part 0
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E87_inhibition_form as E87
from experiments import E88_responder_sets as E88
from experiments import E89_apl_breadth_weighting as E89

EXPERIMENT = 'E91_summed_drive_inhibition'

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

POOLS = (0.5, 1.0)                                  # B1's pooled shares
B1_WEIGHTINGS = (('F', 2.0), ('F', 1.0), ('F', 0.5), ('F', 0.0))
B2_WEIGHTINGS = (('F', 2.0), ('F', 1.0), ('F', 0.0))
B2_POOLS = (1.0,)
REFERENCE = E89.REFERENCE                           # ('F', 2.0): the family's weighting
ARMS = ('B1', 'B2')                                 # the fits as they are; the base with the flies' recruitment
E89_ARM = {'B1': 'C1', 'B2': 'C2'}                  # E89's runs of the same base at s = 0
RAW_CRIT = E87.RAW_CRIT                             # a model response is two or more spikes
PULSE = (0.0, 0.5)                                  # the odor pulse (s; [start, end), the steps of ext_inh_kc)
DRIVE_WINDOW = (0.0, 0.5)                           # E78's window (closed) of the claws carrying half of a drive
PREDICT_FROM = -0.05                                # part A's leaky integration starts here (s), at the drive's level
CO_CANDIDATE = 'gm'                                 # the geometric mean of E88's two estimates of the co-response

# part A
PREDICT_RHO = 0.8                                   # the prediction stands at this correlation or more
PRED_KINDS = ('removed', 'scaled')                  # part A's predictions: the responders removed (primary); the counts

# part C (the entry's rules)
CHOICE_MIN_GAIN = 0.0015
GATE_GAIN = 0.003                                   # (i) the mean gain over the directions
GATE_SINGLES_SHARE = 0.75                           # (ii) 12 of the 16 single-wiring cells (4 wirings x 4 directions)
GATE_RAW_SLACK, GATE_MEAN_SLACK = 0.005, 0.01       # (iii)
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075              # (iv)
CONTROL_RATIO = 0.5
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN, EQUAL_MAX_EXPONENT = 0.003, 0.005, 0.01, 1.0
HALVES = ('A', 'B', 'C', 'D')
SPLIT_OF = {'A': 0, 'B': 0, 'C': 1, 'D': 1}
B3_DIFF = 0.003                                     # outcome (i): B3 differs from s = 1 by more than this

# part D
D_POOL = 1.0
D_MAXFEV = 150
D_STAGES = {'S2': ('S', 2.0), 'S1': ('S', 1.0), 'T2': ('T', 2.0), 'T1': ('T', 1.0)}
D_S_GAIN = 0.0015                                   # S's fitted population this far below UNCHANGED_S['all'] -> T's fits
D_PASS_GAIN = 0.003                                 # the population of eight this far below UNCHANGED_8
D_WIRINGS_LOWER = 6                                 # of the eight
D_RAW_REF, D_RAW_SLACK, D_MEAN_SLACK = 0.073, 0.005, 0.01
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
FLIES_MEAN = E89.FLIES_MEAN

_CACHE: Dict[tuple, object] = {}


def _write_json(path: Path, obj) -> None:
    E89._write_json(path, obj)


# ---------------------------------------------------------------------------
# names

def variant_name(arm: str, w: Tuple[str, float], s: float = 0.0, net: bool = False, fixed: bool = False,
                 control: bool = False) -> str:
    """A setting's variant name: '<arm>_<weighting>' (E89's `wtag`) at s = 0, so that E89's references read as
    'B1_F_k2' etc.; '_s<s>' for s above 0, '_net' for `allow_net_inh_per_claw`, '_fixed' at the fit's APL scales (no
    bisection), '_ctrl' for the control of the setting named before it (which itself runs at s = 0)."""
    v = f'{arm}_{E89.wtag((w[0], float(w[1])))}'
    if float(s) != 0.0:
        v += f'_s{float(s):g}'
    if net:
        v += '_net'
    if fixed:
        v += '_fixed'
    if control:
        v += '_ctrl'
    return v


def parse_variant(v: str) -> dict:
    """`variant_name`'s fields back: arm, measure ('none' for no weighting), exponent, weighting ((measure, exponent),
    E89's NO_WEIGHTING for none), wtag, s (as run: 0 for a control), s_of (the s in the name), net, fixed, control."""
    arm, rest = v.split('_', 1)
    parts = rest.split('_')
    if parts[0] == 'none':
        w, i = E89.NO_WEIGHTING, 1
    else:
        if len(parts) < 2 or not re.fullmatch(r'k[0-9.]+', parts[1]):
            raise ValueError(f'variant {v!r}: no weighting')
        w, i = (parts[0], float(parts[1][1:])), 2
    s, net, fixed, control = 0.0, False, False, False
    for p in parts[i:]:
        if p == 'net':
            net = True
        elif p == 'fixed':
            fixed = True
        elif p == 'ctrl':
            control = True
        elif re.fullmatch(r's[0-9.]+', p):
            s = float(p[1:])
        else:
            raise ValueError(f'variant {v!r}: unknown part {p!r}')
    return dict(arm=arm, measure='none' if w[1] == 0 else w[0], exponent=float(w[1]), weighting=w, wtag=E89.wtag(w),
                s=0.0 if control else s, s_of=s, net=net, fixed=fixed, control=control)


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


def ref_variant(arm: str) -> str:
    """The reference of a base: its family's weighting at s = 0 (B1_F_k2, the unchanged fits; B2_F_k2)."""
    return variant_name(arm, REFERENCE)


def _wirings(arg: str) -> Tuple[str, ...]:
    return {'S': E83.S_WIRINGS, 'T': E83.T_WIRINGS, 'all8': E83.ALL_WIRINGS}.get(arg) or tuple(x for x in arg.split(',') if x)


# ---------------------------------------------------------------------------
# the simulators

def build_simulator91(label: str, **model_kws):
    """E83's `build_simulator` for the lab's input (the wiring's fit, its breadth, the family's weights) with extra
    `fit_mb_model` kwargs (`claw_inh_pool`, ...) in every simulation; returns (simulate, info) with E83's info keys
    'label', 'params', 'kws', 'breadth', 'wAPLKC', 'wAPLPN', 'W', 'orn', 'cells'."""
    from experiments import E10_joint_fit as J
    from experiments.E77_input_screen import _scalar_kws
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator
    P = E83.wiring_parts(label)
    p = dict(P['p'])
    wnorm, W, breadth, orn = P['wnorm'], P['W'], P['breadth'], P['orn']
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], breadth, E83.BREADTH_K)
    sim, _ = make_fixed_param_kc_simulator(orn, initial_params=p, **J.base_kws(None), _wPNKC=W, _wAPLKC=wAPLKC,
                                           _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'], _wPNAPL=wnorm['wPNAPL'],
                                           kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU, **model_kws)
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    if not (np.array_equal(cells['wAPLKC_unscaled'].cell_contents.ravel(), wAPLKC.to_numpy())
            and np.array_equal(cells['wAPLPN_unscaled'].cell_contents.ravel(), wnorm['wAPLPN'].to_numpy())):
        raise RuntimeError("the simulator's unscaled APL weights are not the ones it was built with")
    kws = _scalar_kws(dict(J.base_kws(None), kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU,
                           **{k: v for k, v in model_kws.items() if np.ndim(v) == 0}))
    return sim, dict(label=label, params=p, kws=kws, breadth=breadth, wAPLKC=wAPLKC, wAPLPN=wnorm['wAPLPN'], W=W, orn=orn,
                     cells=cells)


class WiringModel91(E83.WiringModel):
    """E83's `WiringModel` (one wiring at its fit with the lab's input; `run_on(x)` with the tied APL scales times x),
    built with extra `fit_mb_model` kwargs when given (`build_simulator91`), and with the simulator's olfsysm variables
    at hand: `mp`, `rv` and `order` (the olfsysm odor index of each output column)."""

    def __init__(self, label: str, check: bool = True, **model_kws):
        P = E83.wiring_parts(label)
        self.label = label
        for k in ('folder', 'tag', 'p', 'saved', 'target_on', 'target_off', 'size', 'K', 'W'):
            setattr(self, k, P[k])
        self.model_kws = dict(model_kws)
        self.sim, self.info = build_simulator91(label, **model_kws) if model_kws else E83.build_simulator(label)
        self.wAPLKC, self.wAPLPN = self.info['wAPLKC'], self.info['wAPLPN']
        cells = self.info['cells']
        self._kc = cells['wAPLKC_unscaled'].cell_contents
        self._pn = cells['wAPLPN_unscaled'].cell_contents
        self._kc0, self._pn0 = self._kc.copy(), self._pn.copy()
        self.site, self.scores = None, {}
        self.mp, self.rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        self.order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        if check:
            E83._check_saved(label, self.sim, self.p, self.saved)

    def p_off(self) -> dict:
        return dict(self.p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)


def set_pool(mp, s: float = 0.0, net: bool = False) -> None:
    """Set the pooled share (`mp.kc.claw_inh_pool`) and `allow_net_inh_per_claw` on a simulator's olfsysm parameters
    (they stay set for its later simulations; reset with the defaults). An olfsysm without the parameter is accepted
    for s = 0 only."""
    s = float(s)
    if not (np.isfinite(s) and 0.0 <= s <= 1.0):
        raise ValueError(f'claw_inh_pool {s} must be finite and in [0, 1]')
    if hasattr(mp.kc, 'claw_inh_pool'):
        mp.kc.claw_inh_pool = s
    elif s != 0.0:
        raise RuntimeError('the installed olfsysm has no claw_inh_pool (pip install -e ../olfsysm)')
    mp.kc.allow_net_inh_per_claw = bool(net)


def set_ext_inh(mp, rv, per_kc: Optional[np.ndarray], n_odors: int = 0) -> None:
    """The control's static inhibition at the soma (olfsysm `rv.kc.ext_inh`, subtracted from each KC's input while the
    odor is on, times `mp.kc.ext_inh_scale`, which must be 1): `per_kc` (one value per KC, in olfsysm's KC order = the
    outputs' row order) the same for each of the simulator's `n_odors` odors; None removes it."""
    if per_kc is None:
        rv.kc.ext_inh = np.zeros((0, 0))
        return
    v = np.asarray(per_kc, dtype=float).reshape(-1)
    if v.size != int(mp.kc.N) or not np.isfinite(v).all() or (v < 0).any():
        raise ValueError(f'ext_inh: {v.size} values for {int(mp.kc.N)} KCs, all finite and not negative')
    if float(mp.kc.ext_inh_scale) != 1.0:
        raise RuntimeError(f'ext_inh_scale is {mp.kc.ext_inh_scale}, not 1')
    rv.kc.ext_inh = np.repeat(v[:, None], int(n_odors), axis=1)


def apl_on_only_ext(sim, mp, rv, ext: Optional[np.ndarray], n_odors: int):
    """`sim` with the control's static inhibition `ext` present only in runs with APL on: a run with every APL scale at
    0 (the APL-off run of E80's one-odor task, `wAPLKC_scale` and `wAPLPN_scale` 0) runs without it and restores it
    after. The control stands in for what APL's pooled inhibition does, which goes with APL, so the control's APL-off
    run is the base's (checked in `b2_points`). Returns `sim` itself when `ext` is None."""
    if ext is None:
        return sim

    def wrapped(params, *args, **kwargs):
        off = all(float(params.get(k, 1.0)) == 0.0 for k in ('wAPLKC_scale', 'wAPLPN_scale'))
        if not off:
            return sim(params, *args, **kwargs)
        set_ext_inh(mp, rv, None)
        try:
            return sim(params, *args, **kwargs)
        finally:
            set_ext_inh(mp, rv, ext, n_odors)
    return wrapped


# ---------------------------------------------------------------------------
# the claws' net drives, from a run's saved bouton rates and APL activity

def _kc_sum(kc_of_claw: np.ndarray, n_kc: int):
    from scipy.sparse import csr_matrix
    c = np.asarray(kc_of_claw, dtype=int).ravel()
    return csr_matrix((np.ones(c.size), (c, np.arange(c.size))), shape=(int(n_kc), c.size))


def unabsorbed_of(x: np.ndarray, kc_of_claw: np.ndarray, n_kc: int) -> Tuple[np.ndarray, np.ndarray]:
    """(D, U) per KC and time step from the claws' net drive x (claws x time; a claw's PN drive minus APL's inhibition
    of it): D the summed clipped drive, sum over the KC's claws of max(x, 0), and U the inhibition its claws do not
    absorb, sum of max(-x, 0) (olfsysm's summed drive is D - s U)."""
    A = _kc_sum(kc_of_claw, n_kc)
    x = np.asarray(x, dtype=float)
    return np.asarray(A @ np.maximum(x, 0.0)), np.asarray(A @ np.maximum(-x, 0.0))


def n_claws_half(drive: np.ndarray) -> int:
    """The number of a KC's claws that carry half of its summed drive (its claws' window-mean clipped drives, largest
    first; E78's `who_selects`); 0 for no drive."""
    v = np.sort(np.asarray(drive, dtype=float).ravel())[::-1]
    tot = v.sum()
    if not tot > 0:
        return 0
    return int(np.searchsorted(np.cumsum(v), 0.5 * tot) + 1)


def leaky(drive: np.ndarray, dt: float, taum: float, v0: Optional[np.ndarray] = None) -> np.ndarray:
    """Each row's leaky integration of `drive` (KCs x time) as olfsysm integrates a KC's potential, without spikes:
    V(t) = V(t-1) + (-V(t-1) + drive(t)) dt / taum, from `v0` (default the first column, the drive's level)."""
    d = np.asarray(drive, dtype=float)
    v = d[:, 0].copy() if v0 is None else np.asarray(v0, dtype=float).copy()
    out = np.empty_like(d)
    a = float(dt) / float(taum)
    for j in range(d.shape[1]):
        v = v + (-v + d[:, j]) * a
        out[:, j] = v
    return out


def claw_stats(rv, mp, odors: Sequence[int], responders: np.ndarray, predict: bool = False, s_pred: float = 1.0,
               apl: bool = True) -> Dict[str, np.ndarray]:
    """From the last run of a simulator (its saved bouton rates and, with `apl`, APL's activity: `mp.kc.save_inh_sims`
    must have been on), per KC and odor (`odors`: olfsysm odor indices, one per output column): 'D' and 'U' (window means
    over PULSE of `unabsorbed_of`'s D and U, with x = wPNKC @ bouton rates(t) - wAPLKC * APL(t - 1), as olfsysm forms
    it), 'n_half' (`n_claws_half` of the claws' clipped drives averaged over DRIVE_WINDOW, for the `responders`, KCs x
    odors; NaN for the others), 'n_driven' (claws with a positive mean drive); with `predict`, 'peak0' / 'peak1' and
    'area0' / 'area1' (the peak after onset, and the area above the KC's threshold, of `leaky` of D and of D - s_pred U
    from PREDICT_FROM)."""
    from pyolfsysm.mb_model import get_time_index
    t = np.asarray(get_time_index(mp, rv), dtype=float)
    Wp = np.asarray(rv.kc.wPNKC, dtype=float)
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    c2k = np.asarray(rv.kc.claw_to_kc, dtype=int).ravel()
    n_kc = int(mp.kc.N)
    if Wp.shape[0] != c2k.size or wk.size != c2k.size:
        raise RuntimeError(f'{Wp.shape[0]} claws in wPNKC, {wk.size} APL>claw weights, {c2k.size} claw KCs')
    responders = np.asarray(responders, dtype=bool).reshape(n_kc, len(odors))
    t0 = PREDICT_FROM if predict else min(PULSE[0], DRIVE_WINDOW[0])
    t1 = t[-1] if predict else max(PULSE[1], DRIVE_WINDOW[1])
    cols = np.flatnonzero((t >= t0 - 1e-12) & (t <= t1 + 1e-12))
    tc = t[cols]
    m_pulse = (tc >= PULSE[0] - 1e-12) & (tc < PULSE[1] - 1e-12)
    m_drive = (tc >= DRIVE_WINDOW[0] - 1e-12) & (tc <= DRIVE_WINDOW[1] + 1e-12)
    post = tc >= -1e-12
    bs = rv.pn.bouton_sims
    ih = rv.kc.inh_sims if apl else None
    if apl and (len(ih) == 0 or np.asarray(ih[int(odors[0])]).size == 0):
        raise RuntimeError("APL's activity was not saved (mp.kc.save_inh_sims)")
    A = _kc_sum(c2k, n_kc)
    order = np.argsort(c2k, kind='stable')
    claws_of = np.split(order, np.flatnonzero(np.diff(c2k[order])) + 1)
    kc_with = c2k[order][np.r_[0, np.flatnonzero(np.diff(c2k[order])) + 1]]
    claws = {int(k): c for k, c in zip(kc_with, claws_of)}
    shape = (n_kc, len(odors))
    out = {k: np.full(shape, np.nan) for k in ('D', 'U', 'n_half', 'n_driven')}
    if predict:
        out.update({k: np.full(shape, np.nan) for k in ('peak0', 'peak1', 'area0', 'area1')})
        thr = np.asarray(rv.kc.thr, dtype=float).ravel()
        dt, taum = float(mp.time_dt), float(mp.kc.taum)
    prev = np.maximum(cols - 1, 0)
    for j, oi in enumerate(odors):
        bt = np.asarray(bs[int(oi)], dtype=float)
        x = Wp @ bt[:, cols]
        if apl:
            a = np.asarray(ih[int(oi)], dtype=float).reshape(-1, t.size)[0]
            x -= wk[:, None] * a[prev][None, :]
        D = np.asarray(A @ np.maximum(x, 0.0))
        U = np.asarray(A @ np.maximum(-x, 0.0))
        out['D'][:, j] = D[:, m_pulse].mean(axis=1)
        out['U'][:, j] = U[:, m_pulse].mean(axis=1)
        xm = np.maximum(x[:, m_drive], 0.0).mean(axis=1)
        out['n_driven'][:, j] = np.asarray(A @ (xm > 0).astype(float)).ravel()
        for k in np.flatnonzero(responders[:, j]):
            out['n_half'][k, j] = n_claws_half(xm[claws[int(k)]]) if int(k) in claws else 0
        if predict:
            V0 = leaky(D, dt, taum)[:, post]
            V1 = leaky(D - float(s_pred) * U, dt, taum)[:, post]
            out['peak0'][:, j], out['peak1'][:, j] = V0.max(axis=1), V1.max(axis=1)
            out['area0'][:, j] = np.maximum(V0 - thr[:, None], 0.0).sum(axis=1) * dt
            out['area1'][:, j] = np.maximum(V1 - thr[:, None], 0.0).sum(axis=1) * dt
        del x, D, U
    if predict:
        out['thr'] = thr
    return out


def stats_summary(st: Dict[str, np.ndarray], counts: np.ndarray, s: float, prefix: str = 'on') -> dict:
    """A run's row columns from `claw_stats` (APL on: prefix 'on'): the claws carrying half of a responder's drive
    (median and mean over the responder entries), the unabsorbed inhibition over the responders (mean, and its mean
    share of their summed drive) and over all entries, the pooled part (s U over D) over the responders."""
    resp = np.asarray(counts, dtype=float) > RAW_CRIT
    nh = st['n_half'][resp]
    row = {f'claws_half_{prefix}_median': float(np.nanmedian(nh)) if nh.size else np.nan,
           f'claws_half_{prefix}_mean': float(np.nanmean(nh)) if nh.size else np.nan}
    if prefix == 'on':
        U, D = st['U'], st['D']
        with np.errstate(divide='ignore', invalid='ignore'):
            ud = np.where(D > 0, U / D, np.nan)
        row.update(unabsorbed_resp_mean=float(np.nanmean(U[resp])) if resp.any() else np.nan,
                   unabsorbed_all_mean=float(np.nanmean(U)), drive_resp_mean=float(np.nanmean(D[resp])) if resp.any() else np.nan,
                   unabsorbed_over_drive_resp=float(np.nanmean(ud[resp])) if resp.any() else np.nan,
                   pooled_over_drive_resp=float(s * np.nanmean(ud[resp])) if resp.any() else np.nan)
    return row


# ---------------------------------------------------------------------------
# saving and loading runs

def _b_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT) / 'B'
    d.mkdir(parents=True, exist_ok=True)
    return d


def _a_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT) / 'A'
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_run(tag: str, o: Dict[str, tuple], unabsorbed: Optional[pd.DataFrame] = None) -> None:
    d = _b_dir()
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')
    if unabsorbed is not None:
        unabsorbed.to_parquet(d / f'{tag}_unabsorbed.parquet')


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment ({'on', 'off'}: (responses, spike counts)): the local outputs, else the archive."""
    for d in (_b_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'B'):
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


def load_unabsorbed(tag: str) -> pd.DataFrame:
    for d in (_b_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'B'):
        p = d / f'{tag}_unabsorbed.parquet'
        if p.exists():
            return pd.read_parquet(p)
    raise RuntimeError(f'no unabsorbed inhibition saved for {tag}')


def load_any(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E91' this experiment's, 'E89' E89's C runs (the references), 'saved' the wiring's fit."""
    if source == 'saved':
        return E83.saved_fit(wiring)['saved']
    o = E89.load_run(tag) if source == 'E89' else load_run(tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


def _append(path: Path, row: dict) -> None:
    E83._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def _frame(st: Dict[str, np.ndarray], key: str, like: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(st[key], index=like.index, columns=like.columns)


# ---------------------------------------------------------------------------
# B1 / B3 / C on the fits as they are (E83's WiringModel, the APL->KC weights set in place)

def run_stats(model: WiringModel91, x: float, check: Optional[tuple] = None, predict: bool = False) -> Tuple[tuple, dict]:
    """The APL-on run at the tied APL scale factor `x` with APL's activity saved, and its `claw_stats` (all odors; the
    responders two or more spikes); with `check` ((responses, spike counts) of the same point) the counts must equal."""
    model.mp.kc.save_inh_sims = True
    try:
        on = model.run_on(x)
    finally:
        model.mp.kc.save_inh_sims = False
    if check is not None:
        ns = int((on[1].to_numpy(dtype=float) != check[1].to_numpy(dtype=float)).sum())
        if ns:
            raise RuntimeError(f'{model.label}: the run again with APL\'s activity saved differs in {ns} counts')
    st = claw_stats(model.rv, model.mp, model.order, on[1].to_numpy(dtype=float) > RAW_CRIT, predict=predict)
    return on, st


def off_stats(model: WiringModel91) -> dict:
    """The claws carrying half of a responder's drive with APL off (the saved APL-off run, identical for every point of
    a wiring; checked), cached in outputs/.../B/claws_off_<wiring>.json."""
    path = _b_dir() / f'claws_off_{model.label}.json'
    if path.exists():
        return json.loads(path.read_text())
    off = model.sim(model.p_off())
    ns = int((off[1].to_numpy(dtype=float) != model.saved['off'][1].to_numpy(dtype=float)).sum())
    if ns:
        raise RuntimeError(f'{model.label}: the APL-off run differs from the saved one in {ns} counts')
    st = claw_stats(model.rv, model.mp, model.order, off[1].to_numpy(dtype=float) > RAW_CRIT, apl=False)
    row = stats_summary(st, off[1].to_numpy(dtype=float), 0.0, 'off')
    path.write_text(json.dumps(row, indent=1) + '\n')
    return row


def b_reading(model: WiringModel91, wk: Optional[pd.Series], s: float = 0.0, net: bool = False, fixed: bool = False,
              ext: Optional[np.ndarray] = None, label: str = '') -> dict:
    """One setting on a wiring at its fit: the APL->KC weights `wk` (None: the fit's), the pooled share `s`, olfsysm's
    `allow_net_inh_per_claw` with `net`, the control's static inhibition `ext` (per KC); the tied APL scale (E83's
    TIED_ON) bisected to the fit's APL-on count as E89's `weighting_reading` (E69's `bisect_factor`, E83's APL_BRACKET and
    APL_STEP, from 1), or at the fit's scales with `fixed`; then the run again with APL's activity saved for
    `claw_stats`. The APL-off run is the saved fit's. Everything is reset on return."""
    from experiments import E69_kc_firing_ceiling as E69
    t0 = time.time()
    E89.set_weights(model, wk)
    set_pool(model.mp, s, net)
    set_ext_inh(model.mp, model.rv, ext, len(model.order))
    try:
        if fixed:
            x, at_bound = 1.0, False
            on = model.run_on(1.0)
            hist = [dict(stage=label, factor=1.0, n_responses=E83._n(on[0]))]
        else:
            x, on, hist, at_bound = E69.bisect_factor(model.run_on, model.target_on, model.size, E83.APL_BRACKET, E83.APL_STEP,
                                                      f'{label} APL scale')
        on2, st = run_stats(model, x, check=on)
    finally:
        set_ext_inh(model.mp, model.rv, None)
        set_pool(model.mp, 0.0, False)
        E89.set_weights(model, None)
    p = model.p
    return dict(on=on, off=model.saved['off'], apl_factor=float(x), at_bound_apl=bool(at_bound), n_sims=len(hist) + 1,
                params=dict(p, **{k: p[k] * x for k in E83.TIED_ON}), seconds=time.time() - t0, history=hist, stats=st)


def _setting_weights(model: WiringModel91, w: Tuple[str, float]) -> Tuple[pd.Series, pd.Series]:
    """(APL->KC weights, per-KC factors) of a weighting at the family's measure F on a wiring (E89's
    `apl_kc_weights`; F at exponent 2 gives the fit's weights exactly)."""
    if w[0] != 'F':
        raise NotImplementedError(f'E91 runs the family\'s measure F only, not {w[0]}')
    parts = E83.wiring_parts(model.label)
    return E89.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], w[1]), E89.weighting_factors(parts['breadth'], w[1])


def b1_points(out: Path, wirings: Sequence[str], settings: Sequence[dict], runs_name: str, part: str) -> pd.DataFrame:
    """Settings on the fits as they are (B1, B3 and C's runs of that base): each setting a dict(w, s, net, fixed,
    control) with `control` the variant whose control this is (None otherwise; the control runs at s = 0 with the
    static inhibition `control_matrix` of that variant's run on the same wiring). Resumable: rows in `runs_name`, runs
    under outputs/.../B/."""
    import gc
    cc = E89._commits_clean()
    path = out / runs_name
    for wiring in wirings:
        todo = []
        for st in settings:
            v = (variant_name('B1' if not st.get('net') else 'B3', st['w'], st.get('s_of', st.get('s', 0.0)), st.get('net', False),
                              st.get('fixed', False), control=st.get('control') is not None))
            if run_tag(v, wiring) not in _done(path):
                todo.append((v, st))
        if not todo:
            continue
        model = WiringModel91(wiring, check=False)
        offrow = off_stats(model)
        resp = E89._responders(model.saved['off'][1])
        for v, st in todo:
            tag = run_tag(v, wiring)
            pv = parse_variant(v)
            wk, fac = _setting_weights(model, st['w'])
            ext = None
            if st.get('control') is not None:
                ctag = run_tag(st['control'], wiring)
                cm = control_matrix(load_unabsorbed(ctag), parse_variant(st['control'])['s'])
                cm = cm.reindex(model.saved['on'][1].index)
                if cm.isna().any().any():
                    raise RuntimeError(f'{tag}: the control matrix does not cover the model\'s KCs')
                ext = cm.iloc[:, 0].to_numpy(dtype=float)
            res = b_reading(model, wk, s=pv['s'], net=pv['net'], fixed=pv['fixed'], ext=ext, label=tag)
            o = {'on': res['on'], 'off': res['off']}
            U = _frame(res['stats'], 'U', res['on'][1])
            save_run(tag, o, U)
            row = dict(part=part, arm=pv['arm'], tag=tag, variant=v, source='E91', weighting=pv['wtag'], measure=pv['measure'],
                       exponent=pv['exponent'], claw_inh_pool=pv['s'], s_of=pv['s_of'], allow_net_inh_per_claw=pv['net'],
                       fixed_scale=pv['fixed'], control_of=st.get('control') or '', wiring=wiring, permute_seed=np.nan,
                       fit_folder=model.folder, fit_tag=model.tag, apl_factor=res['apl_factor'], apl_min=res['apl_factor'],
                       apl_max=res['apl_factor'], at_bound=res['at_bound_apl'], n_on=E83._n(res['on'][0]),
                       target_on=model.target_on, n_off=E83._n(res['off'][0]), target_off=model.target_off,
                       n_sims=res['n_sims'], seconds=res['seconds'], **{k: res['params'][k] for k in E83.PARAM_KEYS},
                       wAPLKC_sha256=E89.E80_sha(wk), K_fit=model.K, ext_inh_kc_max=float(ext.max()) if ext is not None else np.nan,
                       ext_inh_kc_mean=float(ext.mean()) if ext is not None else np.nan, **E89.factor_summary(fac, resp),
                       **stats_summary(res['stats'], res['on'][1].to_numpy(dtype=float), pv['s']), **offrow, **cc)
            row.update(E89._realism_row(o))
            _append(path, row)
            print(f'E91 {tag}: APL scale x {res["apl_factor"]:.4g}, {row["n_on"]} / {model.target_on} responses, '
                  f'{res["n_sims"]} sims, {res["seconds"]:.0f}s; claws carrying half: {row["claws_half_on_median"]:g}; '
                  f'flag {row.get("realism_flag")!r}', flush=True)
        del model
        gc.collect()
    return pd.read_csv(path)


def b1_settings() -> List[dict]:
    """B1's 32 points' settings (exponents 2, 1, 0.5 and none at s of 0.5 and 1) and the 4 at the fit's scales (s = 1,
    exponent 2)."""
    out = [dict(w=w, s=s) for w in B1_WEIGHTINGS for s in POOLS]
    out.append(dict(w=REFERENCE, s=1.0, fixed=True))
    return out


def part_b1(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    return b1_points(out, wirings, b1_settings(), 'runs_B1.csv', 'B1')


def part_b3(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B3: `allow_net_inh_per_claw=True` at exponent 2 (s = 0), APL scale re-held (4 points)."""
    return b1_points(out, wirings, [dict(w=REFERENCE, s=0.0, net=True)], 'runs_B3.csv', 'B3')


# ---------------------------------------------------------------------------
# B2 / C on the base with the flies' recruitment (E89's C2Pool, with this script's worker)

def _odor_claw_stats(sim, mp, rv, order: np.ndarray, jj: int, params: dict, rec: dict, n_kc: int) -> dict:
    """In a B2 worker, after E80's one-odor task: odor `jj` again at its threshold and APL scale with APL's activity
    saved (checked against the task's counts) and with APL off, and their `claw_stats`."""
    from experiments import E80_screen as S80
    out = {}
    for cond, q in (('on', S80._tied(dict(params), float(rec['apl']))), ('off', dict(params, wAPLKC_scale=0.0, wAPLPN_scale=0.0))):
        mp.kc.save_inh_sims = cond == 'on'
        mp.sim_only = [int(order[jj])]
        try:
            r, s = sim(q, kc_thr_factor=np.full(n_kc, float(rec['thr'])))
        finally:
            mp.sim_only = []
            mp.kc.save_inh_sims = False
        c = s.iloc[:, jj].to_numpy(dtype=float)
        if not np.array_equal(c, np.asarray(rec[cond][1], dtype=float).ravel()):
            raise RuntimeError(f"odor {jj} ({cond}): the run again differs from the task's")
        st = claw_stats(rv, mp, [int(order[jj])], (c > RAW_CRIT)[:, None], apl=cond == 'on')
        out[cond] = {k: st[k][:, 0] for k in ('D', 'U', 'n_half')}
    return out


def _c91_worker(conn, label: str, parts: dict, fit: dict, js: List[int]) -> None:
    """A B2 worker: wiring `label`'s simulator on the ORN columns `js` (as E89's `_c2_worker`), then for each task its
    APL->KC weights ('wAPLKC', per claw; None the fit's), pooled share ('claw_inh_pool'), `allow_net_inh_per_claw`
    ('net') and the control's static inhibition ('ext_inh', per KC; None for none) set, E80's B3 one-odor task run and,
    with 'claw', each odor's `claw_stats`; everything reset after each task; until None."""
    try:
        from experiments import E80_screen as S80
        E83._PARTS[label] = parts
        E83._FITS[label] = fit
        sim, info = E83.build_simulator(label, orn=parts['orn'].iloc[:, js], breadth=parts['breadth'])
        cells = info['cells']
        mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
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
            s = float(msg.pop('claw_inh_pool', 0.0))
            net = bool(msg.pop('net', False))
            ext = msg.pop('ext_inh', None)
            claw = bool(msg.pop('claw', False))
            kc[...] = kc0 if wk is None else np.asarray(wk, dtype=float).reshape(kc0.shape)
            try:
                set_pool(mp, s, net)
                set_ext_inh(mp, rv, ext, len(order))
                # (the control's static inhibition only with APL on; its APL-off run is the base's)
                sim_t = apl_on_only_ext(sim, mp, rv, ext, len(order))
                res = S80._odor_task(sim_t, mp, order, names, dict(info['params']), n_kc, msg)
                if claw:
                    for jj, rec in enumerate(res):
                        rec['claw'] = _odor_claw_stats(sim_t, mp, rv, order, jj, info['params'], rec, n_kc)
            finally:
                kc[...] = kc0
                set_ext_inh(mp, rv, None)
                set_pool(mp, 0.0, False)
            conn.send(('ok', res))
        except Exception as e:
            conn.send(('error', repr(e)))


class C91Pool(E89.C2Pool):
    """E89's C2Pool with this script's worker: `run_task(task, wk, s, net, ext, claw)` runs E80's B3 task with the APL->KC
    weights `wk` (a per-claw Series indexed as the fit's, or None), the pooled share `s`, `allow_net_inh_per_claw` with
    `net`, the control's static inhibition `ext` (per KC, the outputs' row order) and, with `claw`, each odor's claw
    statistics; returns ({'off', 'on1', 'on'}: (responses, spike counts)), per-odor records)."""

    _WORKER = staticmethod(_c91_worker)

    def run_task(self, task: dict, wk: Optional[pd.Series] = None, s: float = 0.0, net: bool = False,
                 ext: Optional[np.ndarray] = None, claw: bool = True):
        keys = ('thr', 'target_off', 'apl', 'target_on')
        arr = None if wk is None else wk.to_numpy(dtype=float)
        ex = None if ext is None else np.asarray(ext, dtype=float)
        for conn, js in zip(self.conns, self.groups):
            conn.send(dict({k: (None if task.get(k) is None else [task[k][j] for j in js]) for k in keys},
                           label=task.get('label', ''), wAPLKC=arr, claw_inh_pool=float(s), net=bool(net), ext_inh=ex,
                           claw=bool(claw)))
        recs = [None] * len(self.cols)
        for conn, js in zip(self.conns, self.groups):
            status, v = conn.recv()
            if status != 'ok':
                raise RuntimeError(f'{self.label}: a worker failed: {v}')
            for j, rec in zip(js, v):
                recs[j] = rec
        saved = E83.saved_fit(self.label)['saved']
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


def _stats_from_recs(recs: List[dict], cond: str) -> Dict[str, np.ndarray]:
    return {k: np.stack([r['claw'][cond][k] for r in recs], axis=1) for k in ('D', 'U', 'n_half')}


def base_thresholds91(w: str, n_workers: int = 16) -> dict:
    """The base's per-odor threshold factors on wiring `w` ({'thr', 'target_off', 'target_on', 'cols'}): E88's (S's
    wirings, outputs/E88_responder_sets/D1/), E89's if it solved them (local or archived), else this experiment's
    (outputs/.../B/base_thr_<w>.json), solved as E89's `base_thresholds` solves them (the identity of the one-odor path
    first, in identity.json; then each odor's threshold bisected to the flies' count with APL output blocked (TNT),
    baseline-noise units), with the base's APL-off run saved beside."""
    from experiments import E80_screen as S80
    for p in (H.outputs_dir(E88.EXPERIMENT) / 'D1' / f'base_thr_{w}.json', H.outputs_dir(E89.EXPERIMENT) / 'C' / f'base_thr_{w}.json',
              H.shared_dir() / 'data' / 'model_outputs' / E88.EXPERIMENT / 'D1' / f'base_thr_{w}.json',
              H.shared_dir() / 'data' / 'model_outputs' / E89.EXPERIMENT / 'C' / f'base_thr_{w}.json',
              _b_dir() / f'base_thr_{w}.json'):
        if p.exists():
            return json.loads(p.read_text())
    cc = E89._commits_clean()
    fit = E83.saved_fit(w)
    pool = C91Pool(w, n_workers)
    try:
        n = len(pool.cols)
        outs, _ = pool.run_task(dict(thr=[1.0] * n, apl=[1.0] * n, label=f'E91 one-odor identity {w}'), claw=False)
        ns, nr = E83._same({'on': outs['on1'], 'off': outs['off']}, fit['saved'])
        _identity_add(H.results_dir(EXPERIMENT), f'one_odor_path_{w}', dict(count_diffs=ns, response_diffs=nr), cc)
        if ns or nr:
            raise RuntimeError(f'{w}: the one-odor path does not reproduce the saved fit')
        tg = S80.b3_targets(E89.C2_UNITS, pool.cols, int(fit['saved']['on'][0].shape[0]))
        outs, recs = pool.run_task(dict(target_off=tg['target_off'], label=f'E91 base {w}'), claw=False)
    finally:
        pool.close()
    path = _b_dir() / f'base_thr_{w}.json'
    path.write_text(json.dumps(dict(thr={r['odor']: float(r['thr']) for r in recs}, target_off=tg['target_off'],
                                    target_on=tg['target_on'], cols=list(pool.cols),
                                    at_bound=[bool(r['at_bound_thr']) for r in recs], **cc), indent=1) + '\n')
    outs['off'][0].to_parquet(_b_dir() / f'base_{w}_acute_off_responses.parquet')
    outs['off'][1].to_parquet(_b_dir() / f'base_{w}_acute_off_spike_counts.parquet')
    return json.loads(path.read_text())


def base_off91(w: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The base's APL-off run on wiring `w` (E88's D1 base for S's wirings, E89's or this experiment's for the others)."""
    for d in (H.outputs_dir(E88.EXPERIMENT) / 'D1', H.outputs_dir(E89.EXPERIMENT) / 'C', _b_dir(),
              H.shared_dir() / 'data' / 'model_outputs' / E88.EXPERIMENT / 'D1',
              H.shared_dir() / 'data' / 'model_outputs' / E89.EXPERIMENT / 'C',
              H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'B'):
        ps = [d / f'base_{w}_acute_off_{k}.parquet' for k in ('responses', 'spike_counts')]
        if all(p.exists() for p in ps):
            return tuple(pd.read_parquet(p) for p in ps)
    raise RuntimeError(f'no base APL-off run for {w}: run base_thresholds91 first')


def b2_points(out: Path, wirings: Sequence[str], settings: Sequence[dict], runs_name: str, part: str,
              n_workers: int = 16, identity: bool = False) -> pd.DataFrame:
    """Settings on the base with the flies' recruitment (B2 and C's runs of that base; dict(w, s, control) as
    `b1_points`): per wiring the base's per-odor thresholds (`base_thresholds91`), then each setting with each odor's
    APL scale (tied) bisected to the control flies' count for that odor (E80's B3, baseline-noise units). The APL-off
    run must equal the base's. With `identity`, first on S's first wiring this path at s = 0 and F at exponent 2 must
    reproduce E89's C2 run exactly (identity.json, 'B2_s0_<wiring>'). Resumable (rows in `runs_name`)."""
    from experiments import E80_screen as S80
    cc = E89._commits_clean()
    path = out / runs_name
    for wiring in wirings:
        todo = []
        for st in settings:
            v = variant_name('B2', st['w'], st.get('s_of', st.get('s', 0.0)), control=st.get('control') is not None)
            if run_tag(v, wiring) not in _done(path):
                todo.append((v, st))
        if not todo:
            continue
        b = base_thresholds91(wiring, n_workers)
        fit = E83.saved_fit(wiring)
        parts = E83.wiring_parts(wiring)
        off0 = base_off91(wiring)
        resp = E89._responders(off0[1])
        pool = C91Pool(wiring, n_workers)
        try:
            thr = [float(b['thr'][c]) for c in pool.cols]
            tg = S80.b3_targets(E89.C2_UNITS, pool.cols, int(fit['saved']['on'][0].shape[0]))
            idp = out / 'identity.json'
            key = f'B2_s0_{wiring}'
            if identity and wiring == E83.S_WIRINGS[0] and key not in (json.loads(idp.read_text()) if idp.exists() else {}):
                wk0 = E89.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], REFERENCE[1])
                outs, recs = pool.run_task(dict(thr=thr, target_on=tg['target_on'], label=f'E91 identity {key}'), wk0, s=0.0,
                                           claw=False)
                ref = E89.load_run(E89.run_tag('C2', REFERENCE, wiring))
                ns, nr = E83._same({'on': outs['on'], 'off': outs['off']}, ref) if ref is not None else (-1, -1)
                _identity_add(out, key, dict(count_diffs=ns, response_diffs=nr, against=f"E89 {E89.run_tag('C2', REFERENCE, wiring)}"), cc)
                print(f'E91 identity {key}: {ns} count / {nr} response differences', flush=True)
                if ns or nr:
                    raise RuntimeError(f"{key}: the B2 path at s = 0 does not reproduce E89's C2 run")
            for v, st in todo:
                tag = run_tag(v, wiring)
                pv = parse_variant(v)
                wk = E89.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], pv['exponent'])
                fac = E89.weighting_factors(parts['breadth'], pv['exponent'])
                ext = None
                if st.get('control') is not None:
                    cm = control_matrix(load_unabsorbed(run_tag(st['control'], wiring)), parse_variant(st['control'])['s'])
                    cm = cm.reindex(fit['saved']['on'][1].index)
                    if cm.isna().any().any():
                        raise RuntimeError(f'{tag}: the control matrix does not cover the model\'s KCs')
                    ext = cm.iloc[:, 0].to_numpy(dtype=float)
                t0 = time.time()
                outs, recs = pool.run_task(dict(thr=thr, target_on=tg['target_on'], label=f'E91 {tag}'), wk, s=pv['s'], ext=ext)
                secs = time.time() - t0
                o = {'on': outs['on'], 'off': outs['off']}
                ns_off, nr_off = E83._same({'on': outs['off'], 'off': outs['off']}, {'on': off0, 'off': off0})
                if ns_off or nr_off:
                    raise RuntimeError(f"{tag}: the APL-off run differs from the base's ({ns_off} counts)")
                st_on, st_off = _stats_from_recs(recs, 'on'), _stats_from_recs(recs, 'off')
                save_run(tag, o, _frame(st_on, 'U', outs['on'][1]))
                apl = np.array([r['apl'] for r in recs])
                row = dict(part=part, arm='B2', tag=tag, variant=v, source='E91', weighting=pv['wtag'], measure=pv['measure'],
                           exponent=pv['exponent'], claw_inh_pool=pv['s'], s_of=pv['s_of'], allow_net_inh_per_claw=False,
                           fixed_scale=False, control_of=st.get('control') or '', wiring=wiring, permute_seed=np.nan,
                           fit_folder=fit['folder'], fit_tag=fit['tag'], apl_factor=float(np.exp(np.log(apl).mean())),
                           apl_min=float(apl.min()), apl_max=float(apl.max()), at_bound=bool(any(r['at_bound_apl'] for r in recs)),
                           n_on=E83._n(outs['on'][0]), target_on=int(sum(tg['target_on'])), n_off=E83._n(outs['off'][0]),
                           target_off=int(sum(b['target_off'])), n_sims=int(sum(r['n_sims'] for r in recs)), seconds=secs,
                           thr_min=float(min(thr)), thr_max=float(max(thr)),
                           apl_per_odor=json.dumps({r['odor']: r['apl'] for r in recs}),
                           thr_per_odor=json.dumps({r['odor']: r['thr'] for r in recs}),
                           **{kk: fit['p'][kk] for kk in E83.PARAM_KEYS}, wAPLKC_sha256=E89.E80_sha(wk), K_fit=fit['K'],
                           ext_inh_kc_max=float(ext.max()) if ext is not None else np.nan,
                           ext_inh_kc_mean=float(ext.mean()) if ext is not None else np.nan,
                           **E89.factor_summary(fac, resp), **stats_summary(st_on, outs['on'][1].to_numpy(dtype=float), pv['s']),
                           **stats_summary(st_off, outs['off'][1].to_numpy(dtype=float), 0.0, 'off'), **cc)
                row.update(E89._realism_row(o))
                _append(path, row)
                print(f'E91 {tag}: APL scales x {row["apl_min"]:.3g} to {row["apl_max"]:.3g}, {row["n_on"]} / {row["target_on"]} '
                      f'responses, {row["n_sims"]} sims, {secs:.0f}s; flag {row.get("realism_flag")!r}', flush=True)
        finally:
            pool.close()
    return pd.read_csv(path)


def part_b2(out: Path, wirings: Sequence[str], n_workers: int = 16) -> pd.DataFrame:
    """B2: s = 1 at exponents 2, 1 and no weighting on the base with the flies' recruitment (12 points)."""
    return b2_points(out, wirings, [dict(w=w, s=s) for w in B2_WEIGHTINGS for s in B2_POOLS], 'runs_B2.csv', 'B2', n_workers,
                     identity=True)


# ---------------------------------------------------------------------------
# part 0: the identities

def _identity_add(out: Path, key: str, d: dict, cc: Optional[dict] = None) -> None:
    path = out / 'identity.json'
    v = json.loads(path.read_text()) if path.exists() else {}
    v[key] = dict(d, **(cc if cc is not None else H.code_commits()))
    _write_json(path, v)


def part_0(out: Path, wirings: Sequence[str] = E83.S_WIRINGS, kwarg_wiring: str = 'real',
           e89_point: Tuple[str, Tuple[str, float]] = ('P1', ('F', 1.0))) -> dict:
    """The identities (the module docstring's part 0); identity.json; raises on any failure."""
    import gc
    import olfsysm
    cc = E89._commits_clean()
    bi = olfsysm.build_info()
    _identity_add(out, 'olfsysm', {k: bi.get(k) for k in ('git_hash', 'git_dirty', 'cpp_dirty', 'source_sha256', 'build_time')}, cc)
    bad = []
    for w in wirings:
        m = WiringModel91(w, check=True)              # the simulator reproduces the saved APL-on and APL-off runs
        set_pool(m.mp, 0.0)
        on0, off0 = m.run_on(1.0), m.sim(m.p_off())
        ns0, nr0 = E83._same({'on': on0, 'off': off0}, m.saved)
        set_pool(m.mp, 1.0)
        try:
            off1, on1 = m.sim(m.p_off()), m.run_on(1.0)
        finally:
            set_pool(m.mp, 0.0)
        ns1, nr1 = E83._same({'on': off1, 'off': off1}, {'on': m.saved['off'], 'off': m.saved['off']})
        n_on1 = int((on1[1].to_numpy(dtype=float) != m.saved['on'][1].to_numpy(dtype=float)).sum())
        d = dict(s0_count_diffs=ns0, s0_response_diffs=nr0, s1_off_count_diffs=ns1 // 2, s1_off_response_diffs=nr1 // 2,
                 s1_on_count_diffs_at_fit_scales=n_on1, s1_on_responses_at_fit_scales=E83._n(on1[0]),
                 target_on=m.target_on, against=f'{m.folder}:{m.tag}')
        if w == kwarg_wiring:
            for s in (0.0, 1.0):
                mk = WiringModel91(w, check=False, claw_inh_pool=s)
                on_k, off_k = mk.run_on(1.0), mk.sim(mk.p_off())
                ref_on = m.saved['on'] if s == 0.0 else on1
                nsk, nrk = E83._same({'on': on_k, 'off': off_k}, {'on': ref_on, 'off': m.saved['off']})
                d[f'kwarg_s{s:g}_count_diffs'], d[f'kwarg_s{s:g}_response_diffs'] = nsk, nrk
                d[f'kwarg_s{s:g}_mp_value'] = float(mk.mp.kc.claw_inh_pool)
                del mk
                gc.collect()
        _identity_add(out, f'S_{w}', d, cc)
        print(f'E91 identity {w}: {json.dumps(d)}', flush=True)
        if ns0 or nr0 or ns1 or nr1 or any(v for k, v in d.items() if k.startswith('kwarg_') and k.endswith('_diffs')):
            bad.append(w)
        if n_on1 == 0:
            bad.append(f'{w}: s = 1 changes nothing with APL on')
        del m
        gc.collect()
    # E89's C1 point through this script's path at s = 0
    w, wt = e89_point
    m = WiringModel91(w, check=False)
    wk, _ = _setting_weights(m, wt)
    res = b_reading(m, wk, s=0.0, label=f'E91 identity E89 C1 {E89.wtag(wt)} {w}')
    ref = E89.load_run(E89.run_tag('C1', wt, w))
    r89 = pd.read_csv(H.results_dir(E89.EXPERIMENT) / 'runs_C1.csv', float_precision='round_trip').set_index('tag')
    r89 = r89.loc[E89.run_tag('C1', wt, w)]
    ns, nr = E83._same({'on': res['on'], 'off': res['off']}, ref) if ref is not None else (-1, -1)
    d = dict(count_diffs=ns, response_diffs=nr, apl_factor=res['apl_factor'], apl_factor_e89=float(r89['apl_factor']),
             against=f'E89 {E89.run_tag("C1", wt, w)}')
    _identity_add(out, f'E89_C1_{E89.wtag(wt)}_{w}', d, cc)
    print(f'E91 identity E89 point: {json.dumps(d)}', flush=True)
    if ns or nr or res['apl_factor'] != float(r89['apl_factor']):
        bad.append(f'E89 C1 {E89.wtag(wt)} {w}')
    if bad:
        raise RuntimeError(f'E91 identities failed: {bad} (identity.json)')
    return json.loads((out / 'identity.json').read_text())


# ---------------------------------------------------------------------------
# part A: what the clipping holds back in the fits as they are

def shared_classes(on: np.ndarray, off: np.ndarray) -> Dict[str, List[Tuple[int, int, int, str]]]:
    """E88's B4 classes (`E88.shared_kcs`, the same rule) of each group's shared responders with APL off: per group a
    list of (pair index, odor a, odor b, KC row, class) entries with class 'survive' (responders to both odors with APL
    on too), 'cut_kept_apl_removed' (kept for both by a per-odor cut of the APL-off counts to the run's APL-on
    fractions, not surviving) or 'rest'. `on` / `off`: KCs x the 16 odors (H.APL16_ODORS order) spike counts."""
    on, off = np.asarray(on, dtype=float), np.asarray(off, dtype=float)
    ron, roff = on > RAW_CRIT, off > RAW_CRIT
    cut = E87.per_odor_oracle(E87.thresholded(off, RAW_CRIT), ron.mean(axis=0), 'divisive', RAW_CRIT) > 0
    i, j = E88._iu()
    masks = E88.group_masks()
    out = {}
    for g in E88.GROUPS + ('all',):
        ent = []
        for p in np.flatnonzero(masks[g]):
            a, b = int(i[p]), int(j[p])
            sh = roff[:, a] & roff[:, b]
            surv = sh & ron[:, a] & ron[:, b]
            ckr = sh & cut[:, a] & cut[:, b] & ~surv
            rest = sh & ~surv & ~ckr
            for cl, msk in (('survive', surv), ('cut_kept_apl_removed', ckr), ('rest', rest)):
                ent += [(int(p), a, b, int(k), cl) for k in np.flatnonzero(msk)]
        out[g] = ent
    return out


def _kc_odor_table(w: str, o: Dict[str, tuple], st: dict, st_off: dict) -> pd.DataFrame:
    """Part A's per KC and odor table of one wiring (the outputs' KC and odor order; odors abbreviated)."""
    on, off = o['on'][1], o['off'][1]
    ids = on.index.get_level_values(H.KC_ID_LEVEL)
    kc = E88._kc_table(w).reindex(ids)
    odors = [H.odor_abbrev(c) for c in on.columns]
    n_kc, n_o = on.shape
    c_on, c_off = on.to_numpy(dtype=float), off.to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(st['area0'] > 0, st['area1'] / st['area0'], np.nan)
        margin_share = np.where(st['peak0'] > st['thr'][:, None], (st['peak0'] - st['peak1']) / (st['peak0'] - st['thr'][:, None]), np.nan)
    # each KC's count times the ratio of its areas above threshold, rounded to whole spikes (no area: unchanged)
    pred = np.where(c_on > 0, np.round(c_on * np.where(np.isfinite(ratio), ratio, 1.0)), c_on)
    rep = lambda a: np.repeat(np.asarray(a), n_o)
    df = pd.DataFrame(dict(wiring=w, kc_id=rep(ids), kc_type=rep(kc['kc_type'].to_numpy()), odor=np.tile(odors, n_kc),
                           n_claws=rep(kc['n_claws'].to_numpy()), breadth=rep(kc['breadth'].to_numpy()),
                           breadth_factor=rep(kc['breadth_factor'].to_numpy()), apl_weight=rep(kc['apl_weight'].to_numpy()),
                           count_on=c_on.ravel(), count_off=c_off.ravel(), resp_on=(c_on > RAW_CRIT).ravel(),
                           resp_off=(c_off > RAW_CRIT).ravel(), D=st['D'].ravel(), U=st['U'].ravel(),
                           thr=rep(st['thr']), peak0=st['peak0'].ravel(), peak1=st['peak1'].ravel(), area0=st['area0'].ravel(),
                           area1=st['area1'].ravel(), margin_share=margin_share.ravel(), n_half_on=st['n_half'].ravel(),
                           n_half_off=st_off['n_half'].ravel(), n_driven_on=st['n_driven'].ravel(), pred_count_scaled=pred.ravel()))
    with np.errstate(divide='ignore', invalid='ignore'):
        df['U_over_D'] = np.where(df.D > 0, df.U / df.D, np.nan)
    df['pred_removed'] = df.resp_on & (df.pred_count_scaled <= RAW_CRIT)
    df['pred_count_removed'] = np.where(df.pred_removed, 0.0, df.count_on)
    return df


def _class_stats(g: pd.DataFrame, **meta) -> dict:
    return dict(meta, n=int(len(g)), U_mean=float(g.U.mean()), U_median=float(g.U.median()), D_mean=float(g.D.mean()),
                U_over_D_mean=float(g.U_over_D.mean()), margin_share_mean=float(g.margin_share.mean()),
                margin_share_median=float(g.margin_share.median()),
                share_pred_removed=float(g.pred_removed.mean()) if 'pred_removed' in g and len(g) else np.nan,
                n_half_on_median=float(g.n_half_on.median()) if 'n_half_on' in g else np.nan)


def class_table(df: pd.DataFrame) -> pd.DataFrame:
    """unabsorbed.csv's rows from part A's per KC and odor table (one wiring or several): the responders with APL on
    (with E88's classes of the ester and the alcohol groups' shared responders, each entry a KC and a pair, read as the
    mean of its two odors), the responders with APL off that APL removes, the rest; the APL-on responders by number of
    claws, by quartile of the breadth weight (over the KCs that respond with APL off) and by KC type."""
    rows = []
    units = [('S' if df.wiring.nunique() > 1 else df.wiring.iloc[0], df)]
    if df.wiring.nunique() > 1:
        units += list(df.groupby('wiring', sort=False))
    for unit, d in units:
        ron = d[d.resp_on]
        rows.append(_class_stats(ron, unit=unit, kind='responders', cls='with APL on'))
        rows.append(_class_stats(d[d.resp_off & ~d.resp_on], unit=unit, kind='responders', cls='with APL off, removed by APL'))
        rows.append(_class_stats(d[~d.resp_off & ~d.resp_on], unit=unit, kind='responders', cls='neither'))
        for nc, g in ron.groupby(np.minimum(ron.n_claws, 12)):
            rows.append(_class_stats(g, unit=unit, kind='n_claws', cls=f'{int(nc)}' + ('+' if nc >= 12 else '')))
        kcs_off = d[d.resp_off].drop_duplicates(['wiring', 'kc_id'])
        q = np.quantile(kcs_off.breadth_factor, [0.25, 0.5, 0.75]) if len(kcs_off) else [np.nan] * 3
        quart = np.searchsorted(q, ron.breadth_factor.to_numpy(), side='right') + 1
        for qq, g in ron.groupby(quart):
            rows.append(_class_stats(g, unit=unit, kind='breadth_quartile', cls=f'Q{int(qq)}'))
        for t, g in ron.groupby('kc_type'):
            rows.append(_class_stats(g, unit=unit, kind='kc_type', cls=str(t)))
        # E88's classes of the groups' shared responders with APL off, per wiring then pooled
        ents = []
        for w, dw in d.groupby('wiring', sort=False):
            piv = {c: dw.pivot(index='kc_id', columns='odor', values=c).reindex(columns=H.APL16_ODORS)
                   for c in ('count_on', 'count_off', 'U', 'D', 'U_over_D', 'margin_share', 'n_half_on')}
            classes = shared_classes(piv['count_on'].to_numpy(), piv['count_off'].to_numpy())
            for grp in ('ester', 'alcohol'):
                for p, a, b, k, cl in classes[grp]:
                    ents.append(dict(group=grp, cls=cl, **{c: 0.5 * (piv[c].iat[k, a] + piv[c].iat[k, b])
                                                            for c in ('U', 'D', 'U_over_D', 'margin_share', 'n_half_on')}))
        if ents:
            e = pd.DataFrame(ents)
            for (grp, cl), g in e.groupby(['group', 'cls']):
                rows.append(_class_stats(g.assign(pred_removed=np.nan), unit=unit, kind=f'{grp}_shared', cls=cl))
    return pd.DataFrame(rows)


def predicted_run(df: pd.DataFrame, o: Dict[str, tuple], kind: str = 'removed') -> Dict[str, tuple]:
    """Part A's predicted run at s = 1 and fixed scales (`kind`, PRED_KINDS): 'removed', the responders predicted to
    fall below two spikes set to 0 and every other count as it is (the entry's prediction of which responders s = 1
    removes); 'scaled', every count times its KC's ratio of supra-threshold areas, rounded. Responses two or more
    spikes; the APL-off run unchanged."""
    on_r, on_s = o['on']
    ids = on_s.index.get_level_values(H.KC_ID_LEVEL)
    odors = [H.odor_abbrev(c) for c in on_s.columns]
    piv = df.pivot(index='kc_id', columns='odor', values=f'pred_count_{kind}').reindex(index=ids, columns=odors)
    s = pd.DataFrame(piv.to_numpy(dtype=float), index=on_s.index, columns=on_s.columns)
    r = (s > RAW_CRIT).astype(on_r.dtypes.iloc[0])
    return {'on': (r, s), 'off': o['off']}


def part_a(out: Path, wirings: Sequence[str] = E83.S_WIRINGS) -> pd.DataFrame:
    """Part A (the module docstring); unabsorbed.csv (class_table over S and each wiring) and, per wiring,
    outputs/.../A/unabsorbed_<wiring>.parquet (per KC and odor) and the predicted runs (A/pred_s1_<kind>_<wiring>_*.parquet)."""
    import gc
    cc = E89._commits_clean()
    frames = []
    for w in wirings:
        path = _a_dir() / f'unabsorbed_{w}.parquet'
        if path.exists():
            frames.append(pd.read_parquet(path))
            continue
        m = WiringModel91(w, check=False)
        set_pool(m.mp, 0.0)
        on, st = run_stats(m, 1.0, check=m.saved['on'], predict=True)
        off = m.sim(m.p_off())
        if int((off[1].to_numpy(dtype=float) != m.saved['off'][1].to_numpy(dtype=float)).sum()):
            raise RuntimeError(f'{w}: the APL-off run differs from the saved one')
        st_off = claw_stats(m.rv, m.mp, m.order, off[1].to_numpy(dtype=float) > RAW_CRIT, apl=False)
        df = _kc_odor_table(w, m.saved, st, st_off)
        df.to_parquet(path)
        for kind in PRED_KINDS:
            po = predicted_run(df, m.saved, kind)
            po['on'][0].to_parquet(_a_dir() / f'pred_s1_{kind}_{w}_on_responses.parquet')
            po['on'][1].to_parquet(_a_dir() / f'pred_s1_{kind}_{w}_on_spike_counts.parquet')
        frames.append(df)
        print(f'E91 A {w}: responders with APL on {int(df.resp_on.sum())}, unabsorbed over them {df[df.resp_on].U.mean():.3g} '
              f'(drive {df[df.resp_on].D.mean():.3g}); predicted removed at s = 1 {int(df.pred_removed.sum())}', flush=True)
        del m
        gc.collect()
    df = pd.concat(frames, ignore_index=True)
    ct = class_table(df)
    ct.to_csv(out / 'unabsorbed.csv', index=False)
    _write_json(out / 'part_A.json', dict(cc, wirings=list(wirings), pulse=PULSE, drive_window=DRIVE_WINDOW,
                                          predict_from=PREDICT_FROM))
    return ct


def _pred_run(w: str, kind: str = 'removed') -> Dict[str, tuple]:
    o = E83.saved_fit(w)['saved']
    ps = [_a_dir() / f'pred_s1_{kind}_{w}_on_{k}.parquet' for k in ('responses', 'spike_counts')]
    if not all(p.exists() for p in ps):
        raise RuntimeError(f'no prediction for {w}: run --part A')
    return {'on': tuple(pd.read_parquet(p) for p in ps), 'off': o['off']}


def part_a_check(out: Path) -> dict:
    """Part A's prediction against B1's simulations: per pair the change of the APL effect (s = 1 at F exponent 2 minus
    the unchanged fits) under each readout, on S's population and each wiring, predicted (A) and simulated (B1's
    re-held point, the entry's check; the run at the fit's scales beside), for both predictions (PRED_KINDS); the
    prediction stands if the correlation of the 'removed' prediction on S's population under the calibrated readout
    against the re-held point is PREDICT_RHO or more ('scaled' and raw counts beside). prediction.csv, prediction.json.
    Also the overlap of the responders predicted and simulated (fixed scales) to be removed."""
    rows, res = [], dict(threshold=PREDICT_RHO, units={})
    sims = {'reheld': variant_name('B1', REFERENCE, 1.0), 'fixed': variant_name('B1', REFERENCE, 1.0, fixed=True)}
    W = E83.S_WIRINGS
    for ro in ('cal', 'raw'):
        d0 = {w: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for w in W}
        ds = {}
        for k, v in sims.items():
            got = {w: load_run(run_tag(v, w)) for w in W}
            if all(x is not None for x in got.values()):
                ds[k] = {w: E83.pair_distances(got[w], ro) for w in W}
        for kind in PRED_KINDS:
            dp = {w: E83.pair_distances(_pred_run(w, kind), ro) for w in W}
            for unit, ws in [('S', W)] + [(w, (w,)) for w in W]:
                pop = lambda dd: E83.population([dd[w] for w in ws])
                e0 = E83._cs(H.apl_effect(*[pop(d0)[c] for c in ('on', 'off')]))
                ep = E83._cs(H.apl_effect(*[pop(dp)[c] for c in ('on', 'off')]))
                pred = ep - e0
                r = dict(readout=ro, unit=unit, prediction=kind)
                for k, dd in ds.items():
                    es = E83._cs(H.apl_effect(*[pop(dd)[c] for c in ('on', 'off')]))
                    simc = (es - e0).reindex(pred.index)
                    r[f'r_{k}'] = float(np.corrcoef(pred.to_numpy(), simc.to_numpy())[0, 1])
                    r[f'mean_change_{k}'] = float(simc.mean())
                    for (a, b), x, y in zip(pred.index, pred.to_numpy(), simc.to_numpy()):
                        rows.append(dict(readout=ro, unit=unit, prediction=kind, sim=k, odor_a=a, odor_b=b, predicted=x, simulated=y))
                r['mean_change_predicted'] = float(pred.mean())
                res['units'][f'{kind}:{unit}:{ro}'] = r
    # the responders removed: predicted against the fixed-scale simulation
    ov = {}
    for w in W:
        o = load_run(run_tag(sims['fixed'], w))
        if o is None:
            continue
        s0 = E83.saved_fit(w)['saved']['on'][1].to_numpy(dtype=float) > RAW_CRIT
        s1 = o['on'][1].to_numpy(dtype=float) > RAW_CRIT
        pr = _pred_run(w, 'removed')['on'][1].to_numpy(dtype=float) > RAW_CRIT
        rem_sim, rem_pred = s0 & ~s1, s0 & ~pr
        inter = int((rem_sim & rem_pred).sum())
        ov[w] = dict(removed_simulated=int(rem_sim.sum()), removed_predicted=int(rem_pred.sum()), both=inter,
                     added_simulated=int((~s0 & s1).sum()),
                     jaccard=float(inter / max(1, int((rem_sim | rem_pred).sum()))))
    res['removed_overlap_fixed_scales'] = ov
    key = res['units'].get('removed:S:cal', {})
    res['stands'] = bool(key.get('r_reheld', np.nan) >= PREDICT_RHO) if 'r_reheld' in key else None
    pd.DataFrame(rows).to_csv(out / 'prediction.csv', index=False)
    _write_json(out / 'prediction.json', res)
    print(f"E91 A check: r (the responders removed, S, calibrated readout) re-held {key.get('r_reheld', np.nan):.3f}, fixed scales "
          f"{key.get('r_fixed', np.nan):.3f}; stands: {res['stands']}", flush=True)
    return res


# ---------------------------------------------------------------------------
# the runs table and the readings

def _e89_reference_rows() -> pd.DataFrame:
    """E89's C1 / C2 rows at the weightings this experiment screens, as this experiment's s = 0 references (B1_<w>,
    B2_<w>; source 'E89'), and the unchanged fits on T (B1_F_k2; source 'saved', their realism computed here)."""
    rows = []
    r89 = H.results_dir(E89.EXPERIMENT)
    for arm, weights in (('B1', B1_WEIGHTINGS), ('B2', B2_WEIGHTINGS)):
        p = r89 / f'runs_{E89_ARM[arm]}.csv'
        if not p.exists():
            continue
        d = pd.read_csv(p)
        d = d[d.permute_seed.isna() & d.weighting.isin([E89.wtag(w) for w in weights])].drop_duplicates('tag', keep='last')
        for r in d.to_dict('records'):
            r.update(part=f'E89_{E89_ARM[arm]}', arm=arm, variant=variant_name(arm, E89._wt_of(f"{E89_ARM[arm]}_{r['weighting']}")),
                     source='E89', claw_inh_pool=0.0, s_of=0.0, allow_net_inh_per_claw=False, fixed_scale=False, control_of='')
            rows.append(r)
    have = {(r['variant'], r['wiring']) for r in rows}
    for w in E83.ALL_WIRINGS:
        if (ref_variant('B1'), w) in have:
            continue
        key = ('saved_row', w)
        if key not in _CACHE:
            f = E83.saved_fit(w)
            _CACHE[key] = dict(part='saved', arm='B1', tag=f'{f["folder"]}:{f["tag"]}', variant=ref_variant('B1'), source='saved',
                               weighting=E89.wtag(REFERENCE), measure='F', exponent=2.0, claw_inh_pool=0.0, s_of=0.0,
                               allow_net_inh_per_claw=False, fixed_scale=False, control_of='', wiring=w, permute_seed=np.nan,
                               fit_folder=f['folder'], fit_tag=f['tag'], apl_factor=1.0, at_bound=False,
                               **{k: f['p'][k] for k in E83.PARAM_KEYS}, **E83.realism(f['saved']))
        rows.append(_CACHE[key])
    return pd.DataFrame(rows)


def runs_table(out: Path) -> pd.DataFrame:
    """Every run of the readings: this experiment's (runs_B1 / B2 / B3 / C.csv), E89's references and the unchanged fits;
    with 'variant', 'flagged' and the parsed setting (s, net, fixed, control)."""
    parts = [pd.read_csv(out / n) for n in ('runs_B1.csv', 'runs_B2.csv', 'runs_B3.csv', 'runs_C.csv') if (out / n).exists()]
    parts.append(_e89_reference_rows())
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    if not len(df):
        return df
    df = df.drop_duplicates(['variant', 'wiring'], keep='first')
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('s', 'net', 'fixed', 'control'):
        df[c] = pv[c].to_numpy()
    df['permute_seed'] = np.nan
    df['flagged'] = df.realism_flag.fillna('').astype(str).str.strip().replace('nan', '') != ''
    df['at_bound'] = df.at_bound.fillna(False).astype(bool)
    return df


def _distances(source: str, tag: str, wiring: str) -> Dict[str, Dict[str, pd.Series]]:
    key = ('dist', source, tag, wiring)
    if key not in _CACHE:
        o = load_any(source, tag, wiring)
        K = E83.saved_fit(wiring)['K']
        _CACHE[key] = {'cal': E83.pair_distances(o, 'cal'), 'calK': E83.pair_distances(o, 'calK', K),
                       'raw': E83.pair_distances(o, 'raw')}
    return _CACHE[key]


def readings91(out: Path, runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings.csv: per variant and unit (each wiring; S and T where all four wirings have a run) and readout ('cal' the
    calibrated readout with K set on each run's APL-on counts, the readings and the gate; 'calK' the fit's K; 'raw' raw
    counts) E83's `effect_scores` against every target of E89's `targets89` (both splits' halves, all recordings, the six
    typical, the noise-free effect, baseline-noise units) and the block of the per-pair effect (E89's `_score_row`); per
    unit the number of runs outside the bands and at a bracket's bound; the setting's fields."""
    runs = runs_table(out) if runs is None else runs
    rows = []
    for (arm, variant), g in runs.groupby(['arm', 'variant'], sort=False):
        ds = {r.wiring: _distances(r.source, r.tag, r.wiring) for r in g.itertuples()}
        meta = {r.wiring: dict(flagged=bool(r.flagged), at_bound=bool(r.at_bound)) for r in g.itertuples()}
        pv = parse_variant(variant)
        units = [(w, (w,)) for w in ds] + [(p, ws) for p, ws in (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS)) if all(w in ds for w in ws)]
        for unit, ws in units:
            for ro in ('cal', 'calK', 'raw'):
                d = ds[unit][ro] if len(ws) == 1 else E83.population([ds[w][ro] for w in ws])
                rows.append(dict(arm=arm, variant=variant, weighting=pv['wtag'], exponent=pv['exponent'], s=pv['s'], s_of=pv['s_of'],
                                 net=pv['net'], fixed=pv['fixed'], control=pv['control'], permute_seed=np.nan, unit=unit,
                                 n_wirings=len(ws), readout=ro, n_flagged=int(sum(meta[w]['flagged'] for w in ws)),
                                 n_at_bound=int(sum(meta[w]['at_bound'] for w in ws)), **E89._score_row(d)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'readings.csv', index=False)
    return df


def _pat(rd: pd.DataFrame, variant: str, unit: str, target: str, col: str = 'pattern', readout: str = 'cal') -> float:
    return E89._pat(rd, variant, unit, target, col, readout)


SPREAD_KEYS = ('apl_effect_amplitude', 'apl_effect_shape', 'apl_effect_sd_model', 'apl_effect_sd_data', 'apl_effect_sd_ratio')


def _spread_rows(d: Dict[str, pd.Series], tg: Dict[str, pd.Series], **meta) -> List[dict]:
    """Per target of `tg`, the pattern part of a run's or population's distances `d` and its amplitude and shape parts
    (`H.apl_effect_spread_shape`; the rules of 2026-11-13 16:22: every pattern part a conclusion rests on is given with
    both parts and the spread ratio), with Pearson's r and ρ on the per-pair APL effect beside."""
    rows = []
    me = E83._cs(H.apl_effect(d['on'], d['off']))
    for t, de0 in tg.items():
        de = E83._cs(de0).reindex(me.index)
        # (H's functions pair the data by H's canonical orientation of the pairs, E83's `_cs` by another: the target
        # goes in H's, or pairs whose orientations differ would be dropped silently)
        sp = H.apl_effect_spread_shape(d['on'], d['off'], H.canonical_pairs(de))
        err = (me - me.mean()) - (de - de.mean())
        pattern = float(np.sqrt(np.mean(err.to_numpy() ** 2)))
        if not np.isclose(pattern ** 2, sp['apl_effect_amplitude'] ** 2 + sp['apl_effect_shape'] ** 2, rtol=1e-9, atol=1e-12):
            raise RuntimeError(f'{t}: the amplitude and shape parts do not add up to the pattern part (pairs misaligned?)')
        rows.append(dict(meta, target=t, pattern=pattern,
                         amplitude=sp['apl_effect_amplitude'], shape=sp['apl_effect_shape'], sd_model=sp['apl_effect_sd_model'],
                         sd_data=sp['apl_effect_sd_data'], sd_ratio=sp['apl_effect_sd_ratio'],
                         pearson=float(np.corrcoef(me.to_numpy(), de.to_numpy())[0, 1]), rho=H.spearman(me, de)))
    return rows


def spread91(out: Path) -> pd.DataFrame:
    """spread.csv: for every variant of runs_*.csv on S and T (populations, where all four wirings have a run) and for D's
    fitted populations (each stage's wirings, and their unchanged fits), under the calibrated readout (K set on each
    run's APL-on counts, the readings' 'cal') and on raw counts, against every target of E89's `targets89`: the pattern
    part, its amplitude and shape parts, the spread ratio, Pearson's r and ρ on the per-pair APL effect."""
    from experiments import E72_wiring_population as E72
    tg = E89.targets89()
    rows = []
    runs = runs_table(out)
    for (arm, variant), g in runs.groupby(['arm', 'variant'], sort=False):
        ds = {r.wiring: _distances(r.source, r.tag, r.wiring) for r in g.itertuples()}
        for unit, ws in (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS)):
            if not all(w in ds for w in ws):
                continue
            for ro in ('cal', 'raw'):
                d = E83.population([ds[w][ro] for w in ws])
                rows += _spread_rows(d, tg, kind='run', arm=arm, variant=variant, unit=unit, readout=ro)
    fits = _fits91(out)
    for stage, (pop, k) in D_STAGES.items():
        got = {}
        for w in _wirings(pop):
            tag = d_tag(w, k)
            if not fits.empty and tag in set(fits.tag):
                o = E72.load_fit(EXPERIMENT, tag)
                if o is not None:
                    got[w] = o
        if len(got) < len(_wirings(pop)):
            continue
        for ro in ('cal', 'raw'):
            d = E83.population([E83.pair_distances(got[w], ro) for w in got])
            d0 = E83.population([E83.pair_distances(E83.saved_fit(w)['saved'], ro) for w in got])
            rows += _spread_rows(d, tg, kind='fit', arm='D', variant=f'D_{stage}', unit=pop, readout=ro)
            rows += _spread_rows(d0, tg, kind='fit', arm='D', variant=f'D_{stage}_unchanged', unit=pop, readout=ro)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'spread.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# part C: the choice, the gate, the control, the equal rule

def _settings(rd: pd.DataFrame, arm: str) -> pd.DataFrame:
    """The arm's settings of the choice on S under the readout: its screen's points and the s = 0 references, not the
    runs at the fit's scales, the controls or B3."""
    return rd[(rd.arm == arm) & (rd.unit == 'S') & (rd.readout == 'cal') & ~rd.fixed.astype(bool) & ~rd.control.astype(bool)
              & ~rd.net.astype(bool)]


def choose91(rd: pd.DataFrame, arm: str, half: str) -> dict:
    """The choice on S against half `half`: the setting with the lowest pattern part of S's population under the
    calibrated readout among the base's reference and the settings whose four runs are inside the bands and not at a
    bracket's bound; it goes to T if its s is above 0 and it gains CHOICE_MIN_GAIN or more over the reference."""
    g = _settings(rd, arm)
    ref = ref_variant(arm)
    ok = g[(g.variant == ref) | ((g.n_flagged == 0) & (g.n_at_bound == 0))]
    if ok.empty or ref not in set(g.variant):
        return dict(arm=arm, half=half, chosen=None, confirm=False, note='runs missing')
    best = ok.loc[ok[f'half{half}_pattern'].idxmin()]
    p_ref = float(g[g.variant == ref][f'half{half}_pattern'].iloc[0])
    gain = p_ref - float(best[f'half{half}_pattern'])
    confirm = bool(float(best.s) > 0 and gain >= CHOICE_MIN_GAIN)
    return dict(arm=arm, half=half, chosen=best.variant, s=float(best.s), weighting=best.weighting,
                pattern=float(best[f'half{half}_pattern']), ref_pattern=p_ref, gain=gain, confirm=confirm, n_eligible=int(len(ok)),
                n_settings=int(len(g)), excluded=json.dumps(sorted(set(g.variant) - set(ok.variant))))


def gate_direction(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, X: str, chosen: Optional[str], ref: str,
                   controls: Sequence[str]) -> dict:
    """One direction of the gate on T (E89's `_gate_direction`: chosen on half X, read against the other half; gain g,
    the single wirings lower, the raw-count loss change, the mean APL effect under the readout, the squared pattern
    error's reduction against all recordings and the six typical, the runs outside the bands, the controls' gains)."""
    r = E89._gate_direction(rd, runs, arm, X, chosen, ref, controls)
    r['split'] = SPLIT_OF[X]
    return r


def gate91(rows: Sequence[dict]) -> dict:
    """E83's gate on the directions given (the entry's four: A to B, B to A, C to D, D to C; a direction whose half
    chose nothing to confirm counts with g 0, no single wiring lower and no reduction): (i) the mean gain over the
    directions GATE_GAIN or more; (ii) at least GATE_SINGLES_SHARE of the single-wiring cells (4 x directions) lower;
    (iii) every run of the confirmed directions inside the bands, T's raw-count loss not above the reference's by more
    than GATE_RAW_SLACK, T's mean APL effect under the readout within GATE_MEAN_SLACK of the flies'; (iv) the squared
    pattern error GATE_SQ_ALL lower against all recordings and GATE_SQ_SIX against the six typical (means over the
    directions); the control: in every confirmed direction with g CHOICE_MIN_GAIN or more, the control gains less than
    CONTROL_RATIO as much. 'passes' all five."""
    rows = list(rows)
    n = len(rows)
    conf = [r for r in rows if r.get('confirmed')]
    gs = [float(r.get('g', 0.0)) if r.get('confirmed') else 0.0 for r in rows]
    out = dict(directions=n, n_confirmed=len(conf), g_mean=float(np.mean(gs)) if n else 0.0, g_min=float(min(gs)) if n else 0.0)
    out['i'] = bool(n and out['g_mean'] >= GATE_GAIN)
    out['n_singles_lower'] = int(sum(int(r.get('n_singles_lower', 0)) for r in conf))
    out['n_singles'] = int(4 * n)
    out['ii'] = bool(n and out['n_singles_lower'] >= GATE_SINGLES_SHARE * out['n_singles'])
    if conf:
        fm = float(conf[0]['flies_mean'])
        worst = max((float(r['mean_effect_cal']) for r in conf), key=lambda x: abs(x - fm))
        out.update(raw_loss_change=float(max(float(r['raw_loss_change']) for r in conf)), mean_effect_cal=worst, flies_mean=fm,
                   n_flagged=int(sum(int(r.get('n_flagged', 0)) for r in conf)))
        out['iii'] = bool(out['n_flagged'] == 0 and out['raw_loss_change'] <= GATE_RAW_SLACK and abs(worst - fm) <= GATE_MEAN_SLACK)
    else:
        out['iii'] = False
    out['sq_reduction_all'] = float(np.mean([float(r.get('sq_reduction_all', 0.0)) if r.get('confirmed') else 0.0 for r in rows])) if n else 0.0
    out['sq_reduction_six'] = float(np.mean([float(r.get('sq_reduction_six', 0.0)) if r.get('confirmed') else 0.0 for r in rows])) if n else 0.0
    out['iv'] = bool(out['sq_reduction_all'] >= GATE_SQ_ALL and out['sq_reduction_six'] >= GATE_SQ_SIX)
    need = [r for r in conf if float(r.get('g', 0.0)) >= CHOICE_MIN_GAIN]
    out['control'] = bool(need and all(r.get('v') is True for r in need))
    out['passes'] = bool(out['i'] and out['ii'] and out['iii'] and out['iv'] and out['control'])
    out['passes_i_to_iv'] = bool(out['i'] and out['ii'] and out['iii'] and out['iv'])
    return out


def equal_rule91(rd: pd.DataFrame, arm: str) -> List[dict]:
    """The equal rule: a setting with s above 0 and an exponent of EQUAL_MAX_EXPONENT or less (no weighting included) is
    an equal of the family's (the base's reference, F at exponent 2 at s = 0) if, on S under the readout, its pattern
    part is within EQUAL_PATTERN of the reference's against all recordings and both halves, its four runs are inside
    the bands, its raw-count loss within EQUAL_RAW and its mean APL effect within EQUAL_MEAN of the reference's;
    'better' if its pattern part is lower on all three."""
    ref = ref_variant(arm)
    g = _settings(rd, arm)
    out = []
    if ref not in set(g.variant):
        return out
    for r in g.itertuples():
        if r.variant == ref or not float(r.s) > 0 or float(r.exponent) > EQUAL_MAX_EXPONENT:
            continue
        d = {t: _pat(rd, r.variant, 'S', t) - _pat(rd, ref, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = _pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - _pat(rd, ref, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - float(g[g.variant == ref].mean_effect.iloc[0])
        eq = bool(all(abs(x) <= EQUAL_PATTERN for x in d.values()) and int(r.n_flagged) == 0 and abs(raw) <= EQUAL_RAW
                  and abs(me) <= EQUAL_MEAN)
        out.append(dict(arm=arm, variant=r.variant, s=float(r.s), exponent=float(r.exponent), **{f'd_pattern_{t}': x for t, x in d.items()},
                        d_raw_loss=raw, d_mean_effect=me, n_flagged=int(r.n_flagged), equal=eq,
                        better=bool(all(x < 0 for x in d.values()))))
    return out


def control_matrix(unabsorbed: pd.DataFrame, s: float) -> pd.DataFrame:
    """The control's static inhibition at the soma (olfsysm `ext_inh_kc`, a KCs x odors matrix): per KC `s` times its
    unabsorbed inhibition averaged over the odors (`unabsorbed`: KCs x odors, window means over the pulse, the chosen
    setting's run), the same for every odor; it keeps how much each KC is inhibited and removes the dependence on which
    claws the odor drives."""
    u = unabsorbed.astype(float)
    v = float(s) * u.mean(axis=1)
    return pd.DataFrame(np.repeat(v.to_numpy()[:, None], u.shape[1], axis=1), index=u.index, columns=u.columns)


def plan_c(rd: pd.DataFrame) -> Dict[Tuple[str, str], dict]:
    return {(arm, h): choose91(rd, arm, h) for arm in ARMS for h in HALVES}


def part_c(out: Path, n_workers: int = 16) -> dict:
    """C: the choice per base and half (choice.csv); for every chosen setting with something to confirm, its runs on T
    (B1: at each wiring's fit; B2: the base on T, its per-odor thresholds solved first, the reference at s = 0 and the
    setting), then its control on S and T. Resumable (runs_C.csv). Then `score`."""
    rd = readings91(out)
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    for arm in ARMS:
        chosen = sorted({c['chosen'] for (a, h), c in plan.items() if a == arm and c['confirm']})
        if not chosen:
            print(f'E91 C {arm}: nothing to confirm', flush=True)
            continue
        T = list(E83.T_WIRINGS)
        for v in chosen:
            pv = parse_variant(v)
            if arm == 'B1':
                b1_points(out, T, [dict(w=pv['weighting'], s=pv['s'])], 'runs_C.csv', 'C')
            else:
                b2_points(out, T, [dict(w=REFERENCE, s=0.0), dict(w=pv['weighting'], s=pv['s'])], 'runs_C.csv', 'C', n_workers)
        for v in chosen:
            pv = parse_variant(v)
            ctl = dict(w=pv['weighting'], s=0.0, s_of=pv['s'], control=v)
            for ws in (list(E83.S_WIRINGS), T):
                if arm == 'B1':
                    b1_points(out, ws, [ctl], 'runs_C.csv', 'C')
                else:
                    b2_points(out, ws, [ctl], 'runs_C.csv', 'C', n_workers)
    return score(out, coresponse=False)


# ---------------------------------------------------------------------------
# the co-response beside every S point (E89's coresponse_c, this experiment's runs)

def _co_job91(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3."""
    source, tag, wiring, K, cand, rel, gap = job
    warnings.simplefilter('ignore')
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    o = load_any(source, tag, wiring)
    mr = E89.ModelReading89(o, label=tag, K=K, rel_target=rel, gap_target=gap)
    res = dict(tag=tag, wiring=wiring, readings={})
    for rdn in ('r1', 'r1b'):
        F = mr.mean_factors(rdn)
        res['readings'][rdn] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    for level in ('L1', 'L3'):
        v = {c: mr.level_values(level, c, None)[0] for c in E89.CONDS}
        F = {c: E89.cand_factors(v[c], cand) for c in E89.CONDS}
        res['readings'][f'cand_{level}'] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    return res


def coresponse91(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse.csv: per S setting (every variant with all four S runs, the controls, the fixed-scale runs and B3
    included) the changes of L (co-response) and R (recruitment) under APL of the four nearest ester pairs, the ester
    and the alcohol groups on S's population (E88's `log_changes`) in readings 1 and 1b and with the geometric mean
    (CO_CANDIDATE) at L1 and L3, against the flies' (90 % intervals over resampled recordings; E88's B1 rule:
    'differs' outside the flies' interval and by a factor of B1_FACTOR or more); 'reaches_flies' if it differs in none
    of the readings. Jobs cached in coresponse_cache.pkl."""
    import multiprocessing
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    runs = runs_table(out)
    runs = runs[runs.wiring.isin(E83.S_WIRINGS)]
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    cache = out / 'coresponse_cache.pkl'
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    jobs = [(r.source, r.tag, r.wiring, E83.saved_fit(r.wiring)['K'], CO_CANDIDATE, rel, gap) for r in runs.itertuples()
            if (r.source, r.tag) not in done]
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            for res, j in zip(ex.map(_co_job91, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'coresponse {j[1]}', flush=True)
    fl = {}
    fv = E89.flies_values({})
    for subset in E88.SUBSETS:
        tnt = E88.flies_subset(subset)
        rf = E88.flies_factors('z', 'r1', tnt)
        F = {c: rf[c].mean() for c in E89.CONDS}
        fl[('r1', subset)] = (E88.log_changes(F, ('L', 'R')), [E88.log_changes(b, ('L', 'R')) for b in E88.flies_boot(rf)])
        fl[('cand', subset)] = (E88.log_changes(E89.flies_cand(fv, CO_CANDIDATE, None, tnt), ('L', 'R')),
                                [E88.log_changes(b, ('L', 'R')) for b in E89.flies_cand_boot(fv, CO_CANDIDATE, None, tnt)])
    rows = []
    for (arm, variant), g in runs.groupby(['arm', 'variant'], sort=False):
        if set(g.wiring) != set(E83.S_WIRINGS):
            continue
        for rdn in ('r1', 'r1b', 'cand_L1', 'cand_L3'):
            Fs = [done[(r.source, r.tag)]['readings'][rdn] for r in g.itertuples()]
            Fp = {c: E88.Factors(*(E88._mean_rows(np.stack([f[c][k] for f in Fs])) for k in ('s', 'O', 'R'))) for c in E89.CONDS}
            ch = E88.log_changes(Fp, ('L', 'R'))
            for subset in E88.SUBSETS:
                fch, fb = fl[('r1' if rdn.startswith('r1') else 'cand', subset)]
                for f in ('L', 'R'):
                    for grp in ('nearest4', 'ester', 'alcohol'):
                        k = f'chg_{f}_{grp}'
                        lo, hi = np.nanpercentile([b[k] for b in fb], E88.INTERVAL)
                        x = ch[k]
                        differs = bool(not (lo <= x <= hi) and abs(x - fch[k]) >= np.log(E89.B1_FACTOR))
                        rows.append(dict(arm=arm, variant=variant, weighting=g.weighting.iloc[0], reading=rdn,
                                         candidate=CO_CANDIDATE if rdn.startswith('cand') else '', subset=subset, factor=f, group=grp,
                                         ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])), flies_lo=float(np.exp(lo)),
                                         flies_hi=float(np.exp(hi)), differs=differs))
    df = pd.DataFrame(rows)
    if len(df):
        reach = df.groupby(['arm', 'variant', 'subset', 'factor', 'group']).differs.any().rename('differs_somewhere').reset_index()
        reach['reaches_flies'] = ~reach.differs_somewhere
        df = df.merge(reach[['arm', 'variant', 'subset', 'factor', 'group', 'reaches_flies']],
                      on=['arm', 'variant', 'subset', 'factor', 'group'])
    df.to_csv(out / 'coresponse.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# scoring and the verdict

def score(out: Path, n_workers: int = 4, coresponse: bool = True) -> dict:
    """readings.csv, coresponse.csv (with `coresponse`), choice.csv, gate.csv, equal.csv, verdict.json, the figure."""
    rd = readings91(out)
    if coresponse:
        coresponse91(out, n_workers)
    v = verdict91(out, rd)
    figure(out)
    comparison_table(out)
    return v


def verdict91(out: Path, rd: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: per base the choices (four halves), the gate against the base's reference (and, for B2, against
    the unchanged fits beside), 'passes' / 'gains on S, not through the gate' / 'level', the equals, the flagged runs;
    part A's prediction (prediction.json, if run); B3 against s = 1 (outcome (i)); D (`score_d`, if fitted); the
    groups' errors per setting on S."""
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table(out)
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    grows, eq_rows = [], []
    v = dict(thresholds=dict(CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=GATE_GAIN, GATE_SINGLES_SHARE=GATE_SINGLES_SHARE,
                             GATE_RAW_SLACK=GATE_RAW_SLACK, GATE_MEAN_SLACK=GATE_MEAN_SLACK, GATE_SQ_ALL=GATE_SQ_ALL,
                             GATE_SQ_SIX=GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO, EQUAL_PATTERN=EQUAL_PATTERN, EQUAL_RAW=EQUAL_RAW,
                             EQUAL_MEAN=EQUAL_MEAN, B3_DIFF=B3_DIFF, PREDICT_RHO=PREDICT_RHO),
             readout="cal (the calibrated readout, K set on each run's APL-on counts; E83's readings and gate)", arms={})
    for arm in ARMS:
        refs = {'own': ref_variant(arm)}
        if arm == 'B2':
            refs['unchanged'] = ref_variant('B1')
        av = dict(choice={h: plan[(arm, h)] for h in HALVES}, gate={})
        for rname, ref in refs.items():
            dirs = []
            for X in HALVES:
                c = plan[(arm, X)]
                chosen = c['chosen'] if c['confirm'] else None
                ctrls = [f'{chosen}_ctrl'] if chosen else []
                r = gate_direction(rd, runs, arm, X, chosen, ref, ctrls)
                r.update(reference_kind=rname)
                dirs.append(r)
                grows.append(r)
            av['gate'][rname] = gate91(dirs)
        own = av['gate']['own']
        gains_s = any(c['confirm'] for c in av['choice'].values())
        av['verdict'] = 'passes' if own['passes'] else ('gains on S, not through the gate' if gains_s else 'level')
        av['equals'] = equal_rule91(rd, arm)
        eq_rows += av['equals']
        g = rd[(rd.arm == arm) & (rd.unit == 'S') & (rd.readout == 'cal')]
        av['S'] = {r.variant: dict(halfA=r.halfA_pattern, halfB=r.halfB_pattern, halfC=r.halfC_pattern, halfD=r.halfD_pattern,
                                   all=r.z8_pattern, loss=r.z8_loss, mean_effect=r.mean_effect, rho=r.z8_rho,
                                   err_ester=r.z8_err_ester, err_alcohol=r.z8_err_alcohol, n_flagged=int(r.n_flagged),
                                   n_at_bound=int(r.n_at_bound)) for r in g.itertuples()}
        rr = runs[(runs.arm == arm) & runs.flagged]
        av['flagged_runs'] = rr[['tag', 'realism_flag']].values.tolist()
        v['arms'][arm] = av
    pd.DataFrame(grows).to_csv(out / 'gate.csv', index=False)
    pd.DataFrame(eq_rows).to_csv(out / 'equal.csv', index=False)
    # B3 against s = 1 (outcome (i))
    b3, s1 = variant_name('B3', REFERENCE, net=True), variant_name('B1', REFERENCE, 1.0)
    if ((rd.variant == b3) & (rd.unit == 'S')).any() and ((rd.variant == s1) & (rd.unit == 'S')).any():
        d = {t: _pat(rd, b3, 'S', t) - _pat(rd, s1, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        v['B3'] = dict(d_pattern=d, differs=bool(any(abs(x) > B3_DIFF for x in d.values())))
    if (out / 'prediction.json').exists():
        p = json.loads((out / 'prediction.json').read_text())
        v['A_prediction'] = dict(stands=p.get('stands'), **{k: p['units'].get(k) for k in ('removed:S:cal', 'scaled:S:cal', 'removed:S:raw')})
    if (out / 'd_scores.json').exists():
        v['D'] = json.loads((out / 'd_scores.json').read_text())
    if (out / 'coresponse.csv').exists():
        co = pd.read_csv(out / 'coresponse.csv')
        if len(co):
            reach = co[(co.factor == 'L') & (co.group == 'alcohol')].drop_duplicates(['arm', 'variant', 'subset'])
            v['alcohol_coresponse_reaches_flies'] = {f'{r.arm}:{r.variant}:{r.subset}': bool(r.reaches_flies) for r in reach.itertuples()}
    _write_json(out / 'verdict.json', v)
    for arm in ARMS:
        a = v['arms'][arm]
        print(f"E91 C {arm}: {a['verdict']}; " + '; '.join(f"{h} {a['choice'][h].get('chosen')} (gain {a['choice'][h].get('gain', np.nan):.4f})"
                                                           for h in HALVES), flush=True)
    return v


# ---------------------------------------------------------------------------
# part D: the family fitted again

def d_member(label: str) -> dict:
    """A wiring's fit of C2's family (E83's `members`) with its E10_joint_fit.py wiring arguments (from its fits.csv
    row: --wiring / --wiring-seed, or --wiring-file / --wiring-tag)."""
    folder, tag = E83.members()[label]
    row = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[tag]
    args: List[str] = []
    wiring = row.get('wiring')
    if isinstance(wiring, str) and wiring:
        args += ['--wiring', wiring, '--wiring-seed', str(int(row['wiring_seed']))]
    wf = row.get('wiring_file')
    if isinstance(wf, str) and wf:
        args += ['--wiring-file', wf, '--wiring-tag', str(row['wiring_tag'])]
    return dict(label=label, folder=folder, tag=tag, wiring_args=args, k=float(row['k']))


def d_command(label: str, k: float, pool: float = D_POOL, maxfev: int = D_MAXFEV) -> List[str]:
    """The E10_joint_fit.py command of one D fit: C2's family (`-k <k> --aplpn-max 3 --rate-target ctrl_frac --loss apl
    --apl-readout both --tag _cap3`) with `--claw-inh-pool`, into this experiment's folders, from the wiring's own fit
    (exactly at the same exponent, `--init-from-exact`; from its fitted factors otherwise)."""
    m = d_member(label)
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', f'{k:g}', '--aplpn-max', '3',
           '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', '_cap3', '--maxfev', str(int(maxfev)),
           '--claw-inh-pool', f'{pool:g}', '--experiment', EXPERIMENT, '--init-from', m['tag'],
           '--init-from-experiment', m['folder']] + m['wiring_args']
    if float(k) == m['k']:
        cmd.append('--init-from-exact')
    return cmd


def d_tag(label: str, k: float, pool: float = D_POOL) -> str:
    """The tag E10_joint_fit.py gives a D fit (its prefix rules: the wiring, rate0.068_, aplboth_, pool<s>_, then
    k<k>_cap3)."""
    m = d_member(label)
    a = m['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    from experiments import E10_joint_fit as J
    return f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.claw_inh_pool_tag(pool)}k{k:g}_cap3"


def _fits91(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def d_allowed(out: Path, stage: str) -> Tuple[bool, str]:
    """The entry's rules for the conditional fits: S1 only if B1's s = 1 point at exponent 1 is an equal or better;
    T2 / T1 only if S's fitted population (S2) is D_S_GAIN or more below the unchanged fits' pattern part against all
    recordings under the readout."""
    if stage == 'S2':
        return True, 'always'
    if stage == 'S1':
        rd = pd.read_csv(out / 'readings.csv') if (out / 'readings.csv').exists() else None
        if rd is None:
            return False, 'no readings (run --part score after B1)'
        v = variant_name('B1', ('F', 1.0), 1.0)
        e = [x for x in equal_rule91(rd, 'B1') if x['variant'] == v]
        ok = bool(e and (e[0]['equal'] or e[0]['better']))
        return ok, f"{v}: {'equal or better' if ok else 'neither equal nor better'} ({e[0] if e else 'not read'})"
    src = 'S' + stage[1:]
    sd = score_d(out, write=False)
    p = sd.get(src, {}).get('S', {}).get('pattern_all')
    if p is None:
        return False, f'{src} not scored'
    ok = bool(p <= UNCHANGED_S['all'] - D_S_GAIN)
    return ok, f'{src} population pattern part {p:.4f} against {UNCHANGED_S["all"]} - {D_S_GAIN}'


def part_d(out: Path, stage: str = 'S2', force: bool = False, dry_run: bool = False) -> List[str]:
    """Runs the stage's fits one after another (`d_command`; a process each, its log under outputs/.../D/), skipping
    those already in this experiment's fits.csv; refuses a conditional stage the entry's rules do not allow (`d_allowed`)
    unless `force`. Returns the stage's tags."""
    pop, k = D_STAGES[stage]
    ok, why = d_allowed(out, stage)
    print(f'E91 D {stage}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
    if not ok and not force:
        return []
    E89._commits_clean()
    logs = H.outputs_dir(EXPERIMENT) / 'D'
    logs.mkdir(parents=True, exist_ok=True)
    tags = []
    for w in _wirings(pop):
        tag = d_tag(w, k)
        tags.append(tag)
        if tag in set(_fits91(out).tag):
            print(f'E91 D {tag}: done', flush=True)
            continue
        cmd = d_command(w, k)
        print('E91 D:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E91 D {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits91(out).tag):
            raise RuntimeError(f'E91 D: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def score_d(out: Path, write: bool = True) -> dict:
    """D's scores: per stage, wiring and population (S, and the eight when T's fits exist) E83's `effect_scores` under
    the calibrated readout and raw counts, against the unchanged fits of the same wirings (E72's part B and C2), with the
    realism; the verdict of the family as a population of eight (the entry's five rules). d_scores.csv, d_scores.json."""
    from experiments import E72_wiring_population as E72
    fits = _fits91(out)
    res, rows = {}, []
    tg = E89.targets89()
    for stage, (pop, k) in D_STAGES.items():
        got = {}
        for w in _wirings(pop):
            tag = d_tag(w, k)
            if tag in set(fits.tag):
                o = E72.load_fit(EXPERIMENT, tag)
                if o is not None:
                    got[w] = (tag, o)
        if not got:
            continue
        st = {}
        for w, (tag, o) in got.items():
            d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
            d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
            sc, sc0 = E83.effect_scores(d['cal'], tg), E83.effect_scores(d0['cal'], tg)
            raw, raw0 = E83.effect_scores(d['raw'], tg), E83.effect_scores(d0['raw'], tg)
            fr = fits[fits.tag == tag].iloc[0]
            rl = E83.realism(o)
            row = dict(stage=stage, unit=w, tag=tag, pattern_all=sc['z8_pattern'], pattern_all_unchanged=sc0['z8_pattern'],
                       halfA=sc['halfA_pattern'], halfB=sc['halfB_pattern'], loss=sc['z8_loss'], mean_effect=sc['mean_effect'],
                       raw_loss=raw['z8_loss'], raw_loss_unchanged=raw0['z8_loss'], err_ester=sc['z8_err_ester'],
                       err_alcohol=sc['z8_err_alcohol'], lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                       **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                                 'kc_kc_inh', 'claw_inh_pool', 'n_evals')},
                       fit_loss=fr.get('loss'), **rl)
            rows.append(row)
            st[w] = (d, d0, row)
        rs = {}
        for unit, ws in (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS)):
            if all(w in st for w in ws):
                dp = {ro: E83.population([st[w][0][ro] for w in ws]) for ro in ('cal', 'raw')}
                d0 = {ro: E83.population([st[w][1][ro] for w in ws]) for ro in ('cal', 'raw')}
                sc, sc0 = E83.effect_scores(dp['cal'], tg), E83.effect_scores(d0['cal'], tg)
                rs[unit] = dict(pattern_all=sc['z8_pattern'], pattern_all_unchanged=sc0['z8_pattern'], halfA=sc['halfA_pattern'],
                                halfB=sc['halfB_pattern'], mean_effect=sc['mean_effect'], loss=sc['z8_loss'],
                                raw_loss=E83.effect_scores(dp['raw'], tg)['z8_loss'],
                                raw_loss_unchanged=E83.effect_scores(d0['raw'], tg)['z8_loss'], n_lower=int(sum(st[w][2]['lower'] for w in ws)),
                                n_flagged=int(sum(bool(str(st[w][2].get('realism_flag') or '').strip() not in ('', 'nan')) for w in ws)))
                rows.append(dict(stage=stage, unit=unit, **rs[unit]))
        res[stage] = rs
    # the family as a population of eight (S's and T's fits at the same exponent)
    for k in (2.0, 1.0):
        sS, sT = f'S{k:g}', f'T{k:g}'
        if sS in res and sT in res and 'S' in res[sS] and 'T' in res[sT]:
            ds = []
            for stage in (sS, sT):
                pop, kk = D_STAGES[stage]
                for w in _wirings(pop):
                    ds.append(E83.pair_distances(E72.load_fit(EXPERIMENT, d_tag(w, kk)), 'cal'))
            dr = []
            for stage in (sS, sT):
                pop, kk = D_STAGES[stage]
                for w in _wirings(pop):
                    dr.append(E83.pair_distances(E72.load_fit(EXPERIMENT, d_tag(w, kk)), 'raw'))
            sc = E83.effect_scores(E83.population(ds), tg)
            raw = E83.effect_scores(E83.population(dr), tg)['z8_loss']
            n_lower = res[sS]['S']['n_lower'] + res[sT]['T']['n_lower']
            n_flag = res[sS]['S']['n_flagged'] + res[sT]['T']['n_flagged']
            eight = dict(pattern_all=sc['z8_pattern'], mean_effect=sc['mean_effect'], raw_loss=raw, n_lower=n_lower, n_flagged=n_flag)
            eight['passes'] = bool(sc['z8_pattern'] <= UNCHANGED_8 - D_PASS_GAIN and n_lower >= D_WIRINGS_LOWER and n_flag == 0
                                   and raw <= D_RAW_REF + D_RAW_SLACK and abs(sc['mean_effect'] - FLIES_MEAN) <= D_MEAN_SLACK)
            res[f'eight_k{k:g}'] = eight
    if write:
        pd.DataFrame(rows).to_csv(out / 'd_scores.csv', index=False)
        _write_json(out / 'd_scores.json', res)
    return res


# ---------------------------------------------------------------------------
# part E (only on a pass)

def part_e(out: Path) -> None:
    """Part E (only on a pass by C on the fits as they are, or by D): for the passing setting, the firing-ceiling family
    (E69's fit 1 and E72 part G's fits) with the response rate with APL on re-held; S at the flies' 2 s pulse with both
    response rates re-held; S through the trial noise as E72's part E ran it; the block under the readout and the
    natural-source alignment with APL on, off and their difference; each holds if the population keeps at least half of
    its gain at 0.5 s without noise under the same reading. Not written: it runs only on a pass."""
    raise NotImplementedError('E91 part E runs only on a pass by part C or part D; it is not written yet')


# ---------------------------------------------------------------------------
# the archive, the figure and the table

def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): one row per run of this experiment (B and C) with absolute parameters (the
    APL scales as run; on the base the per-odor threshold factors and APL scales in JSON), the pooled share, the
    weighting, the control's static inhibition, the factors' range and realism; and D's fits (fits.csv rows), whose
    distances E10_joint_fit.py wrote are kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E91']
    rows, parts = [], []
    for _, r in runs.iterrows():
        base = ('the fit as it is' if r.arm in ('B1', 'B3') else "the base with the flies' recruitment (E88's per-odor thresholds, each "
                "odor's APL scale re-held at the control flies' count for that odor)")
        kind = (f"{base}; APL->KC weights x ((F + 1) / mean(F + 1))^{r.exponent:g}; each KC's summed drive sum max(x, 0) + "
                f"{r.claw_inh_pool:g} sum min(x, 0) (olfsysm claw_inh_pool)"
                + ('; allow_net_inh_per_claw' if r.allow_net_inh_per_claw else '')
                + (f"; the control of {r.control_of}: a static inhibition at the soma (ext_inh_kc) per KC" if str(r.control_of) not in ('', 'nan') else '')
                + ('; at the fit\'s APL scales' if r.fixed_scale else '; the APL scale re-held') + "; APL-off run the fit's / the base's")
        row = dict(tag=r.tag, variant=r.tag, part=r.part, arm=r.arm, kind=kind, config=r.variant, breadth_def='F', breadth_k=r.exponent,
                   k=r.exponent, claw_inh_pool=r.claw_inh_pool, allow_net_inh_per_claw=r.allow_net_inh_per_claw, wiring=r.wiring,
                   ref_experiment=r.fit_folder, ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                   wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                   apl_factor=r.apl_factor, scale_min=r.apl_min, scale_max=r.apl_max, rate_on=r.rate_on, rate_off=r.rate_off,
                   max_count_on=r.max_count_on, max_count_off=r.max_count_off, wAPLKC_sha256=r.wAPLKC_sha256,
                   ext_inh_kc_max=r.get('ext_inh_kc_max'), control_of=r.control_of, realism_flag=r.realism_flag)
        if r.arm == 'B2':
            row.update(scaled_odors='per odor', multiplier_min=r.thr_min, multiplier_max=r.thr_max,
                       thr_per_odor=r.thr_per_odor, apl_per_odor=r.apl_per_odor)
        rows.append(row)
        o = load_run(r.tag)
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    long_path = out / 'model_pair_distances_long.csv'
    keep = []
    if long_path.exists():
        old = pd.read_csv(long_path)
        keep = [old[old.tag.isin(set(_fits91(out).tag))]]           # D's fits, written by E10_joint_fit.py
    pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def figure(out: Path) -> Optional[Path]:
    """E91_screen.png: per base, S's pattern part under the readout against half A, half B and all recordings by setting
    (the s = 0 references and the settings at s above 0; the base's reference dashed; hollow: a run outside the bands);
    and part A's predicted against simulated change of the per-pair APL effect, if checked."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not (out / 'readings.csv').exists():
        return None
    rd = pd.read_csv(out / 'readings.csv')
    has_pred = (out / 'prediction.csv').exists()
    fig, axes = plt.subplots(1, 3 if has_pred else 2, figsize=(17 if has_pred else 12, 4.8))
    for ax, arm in zip(axes, ARMS):
        g = _settings(rd, arm).sort_values(['s', 'exponent'], ascending=[True, False])
        x = np.arange(len(g))
        for t, col, lab in (('halfA', '#d1495b', 'half A'), ('halfB', '#2e86ab', 'half B'), ('z8', '0.2', 'all recordings')):
            y = g[f'{t}_pattern'].to_numpy()
            fl = (g.n_flagged > 0).to_numpy()
            ax.plot(x, y, '-', color=col, lw=0.8, label=lab)
            ax.scatter(x, y, s=28, facecolors=['none' if f else col for f in fl], edgecolors=col)
            ref = g[g.variant == ref_variant(arm)][f'{t}_pattern']
            if len(ref):
                ax.axhline(float(ref.iloc[0]), color=col, ls='--', lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([v.split('_', 1)[1] for v in g.variant], rotation=60, fontsize=8)
        ax.set_ylabel('pattern part of the APL-effect loss')
        ax.set_title(("the fits as they are" if arm == 'B1' else "the base with the flies' recruitment")
                     + ": S's population\n(calibrated readout; hollow: a run outside the bands)", fontsize=9)
        ax.legend(fontsize=8)
    if has_pred:
        p = pd.read_csv(out / 'prediction.csv')
        p = p[(p.unit == 'S') & (p.readout == 'cal') & (p.prediction == 'removed')]
        ax = axes[2]
        for k, col in (('reheld', '#2e86ab'), ('fixed', '#d1495b')):
            q = p[p.sim == k]
            if len(q):
                r = np.corrcoef(q.predicted, q.simulated)[0, 1]
                ax.scatter(q.predicted, q.simulated, s=10, color=col,
                           label=f"{'APL scale re-held' if k == 'reheld' else 'at the fit scales'} (r {r:.2f})")
        ax.axline((0, 0), slope=1, color='0.5', lw=0.7)
        ax.set_xlabel('predicted change of the per-pair APL effect (part A, the responders removed)')
        ax.set_ylabel('simulated change (s = 1, F at exponent 2)')
        ax.set_title("S's population, calibrated readout", fontsize=9)
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / 'E91_screen.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _f(x, spec: str = '.4f') -> str:
    return E87._f(x, spec)


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    L = ["# E91: APL's inhibition on the KC's summed drive instead of claw by claw", '',
         f'Generated by experiments/E91_summed_drive_inhibition.py at {H.code_commits().get("queue_commit", "")[:7]}.', '']
    if (out / 'identity.json').exists():
        L += ['## 0: identities', '', '```', (out / 'identity.json').read_text()[:6000], '```', '']
    if (out / 'unabsorbed.csv').exists():
        u = pd.read_csv(out / 'unabsorbed.csv')
        u = u[u.unit == 'S']
        L += ["## A: the inhibition the claws do not absorb (S, the fits as they are, APL on; window means over the pulse)", '',
              '| kind | class | entries | U mean | U median | D mean | U / D | share of the margin | predicted removed at s = 1 | claws carrying half |',
              '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in u.itertuples():
            L.append(f'| {r.kind} | {r.cls} | {r.n} | {_f(r.U_mean, ".3g")} | {_f(r.U_median, ".3g")} | {_f(r.D_mean, ".3g")} | '
                     f'{_f(r.U_over_D_mean, ".3f")} | {_f(r.margin_share_mean, ".3f")} | {_f(r.share_pred_removed, ".3f")} | '
                     f'{_f(r.n_half_on_median, ".1f")} |')
        L.append('')
    if (out / 'prediction.json').exists():
        p = json.loads((out / 'prediction.json').read_text())
        L += ["## A: the prediction against B1's simulations (correlation of the per-pair change of the APL effect)", '',
              '| prediction:unit:readout | re-held (the check) | at the fit\'s scales | mean change predicted | re-held | fixed |', '|---|---:|---:|---:|---:|---:|']
        for k, r in p['units'].items():
            L.append(f"| {k} | {_f(r.get('r_reheld'), '.3f')} | {_f(r.get('r_fixed'), '.3f')} | {_f(r.get('mean_change_predicted'), '+.4f')} | "
                     f"{_f(r.get('mean_change_reheld'), '+.4f')} | {_f(r.get('mean_change_fixed'), '+.4f')} |")
        L += ['', f"Stands (the responders removed, S, calibrated readout, re-held, r >= {PREDICT_RHO}): {p.get('stands')}.", '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        runs = runs_table(out)
        for arm in ARMS + ('B3',):
            g = rd[(rd.arm == arm) & (rd.unit == 'S') & (rd.readout == 'cal')]
            if g.empty:
                continue
            graw = rd[(rd.arm == arm) & (rd.unit == 'S') & (rd.readout == 'raw')].set_index('variant')
            rr = runs[runs.arm == arm].groupby('variant')
            L += [f"## B, {dict(B1='the fits as they are (B1)', B2='the base with the flies recruitment (B2)', B3='allow_net_inh_per_claw (B3)')[arm]}: S's population", '',
                  '| setting | pattern, half A | half B | all | loss | mean APL effect | ρ | raw-count loss | ester / alcohol error | block | '
                  'APL scale factor (range) | claws carrying half, on / off | flagged / at bound |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            for r in g.itertuples():
                q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()
                fac = (f"{q.apl_factor.min():.3g}-{q.apl_factor.max():.3g}" if len(q) and 'apl_factor' in q else '–')
                ch = (f"{_f(q.claws_half_on_median.mean(), '.1f')} / {_f(q.claws_half_off_median.mean(), '.1f')}"
                      if len(q) and 'claws_half_on_median' in q else '–')
                L.append(f'| {r.variant} | {r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | {r.z8_loss:.4f} | '
                         f'{r.mean_effect:+.3f} | {r.z8_rho:.2f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | '
                         f'{r.z8_err_ester:+.3f} / {r.z8_err_alcohol:+.3f} | {r.block:+.3f} | {fac} | {ch} | {r.n_flagged} / {r.n_at_bound} |')
            L.append('')
        rr = runs[runs.source == 'E91']
        if len(rr):
            L += ['## Realism of the runs (rule of 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | flag | flag (earlier rule) |',
                  '|---|---:|---:|---:|---:|---:|---:|---|---|']
            for r in rr.itertuples():
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | {r.opk_on_tv:.3f} / '
                         f'{r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | {r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
    for name in ('choice.csv', 'gate.csv', 'equal.csv', 'd_scores.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                L += [f'## {name}', '', '(none)', '']
                continue
            try:
                txt = d.to_markdown(index=False, floatfmt='.4f')
            except ImportError:
                txt = '```\n' + d.to_string() + '\n```'
            L += [f'## {name}', '', txt, '']
    if (out / 'verdict.json').exists():
        vc = json.loads((out / 'verdict.json').read_text())
        short = {a: dict(verdict=x['verdict'], gate=x['gate'], equals=[e for e in x['equals'] if e['equal']])
                 for a, x in vc.get('arms', {}).items()}
        short.update({k: vc[k] for k in ('B3', 'A_prediction', 'D') if k in vc})
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json: the script, the commits (H.code_commits()), olfsysm's build, the model's provenance with
    the lab's input (H.model_provenance) for the simulation parts, the settings and `extra`; a part run again adds a
    dated record."""
    import olfsysm
    path = out / f'provenance_{part}.json'
    old = json.loads(path.read_text()) if path.exists() else {}
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('0', 'A', 'B1', 'B2', 'B3', 'C', 'D'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    _write_json(path, dict(script='experiments/E91_summed_drive_inhibition.py', part=part,
                           settings=dict(POOLS=POOLS, B1_WEIGHTINGS=B1_WEIGHTINGS, B2_WEIGHTINGS=B2_WEIGHTINGS, B2_POOLS=B2_POOLS,
                                         PULSE=PULSE, DRIVE_WINDOW=DRIVE_WINDOW, PREDICT_FROM=PREDICT_FROM, C2_UNITS=E89.C2_UNITS,
                                         CO_CANDIDATE=CO_CANDIDATE, D_POOL=D_POOL, D_MAXFEV=D_MAXFEV),
                           runs=runs))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('0', 'A', 'Acheck', 'B1', 'B2', 'B3', 'co', 'C', 'score', 'D', 'Dscore', 'E', 'variants', 'table',
                             'spread'))
    ap.add_argument('--workers', type=int, default=16, help='B2 / C: one-odor workers (16); co / score: processes (4)')
    ap.add_argument('--wirings', default='S', help="0 / A / B1 / B2 / B3: 'S', 'T' or a comma-separated list")
    ap.add_argument('--stage', default='S2', choices=list(D_STAGES), help='D: which fits')
    ap.add_argument('--force', action='store_true', help="D: run a conditional stage the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='D: print the commands only')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    args = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = H.results_dir(EXPERIMENT)
    ws = _wirings(args.wirings)
    if args.part in ('0', 'A', 'B1', 'B2', 'B3', 'C', 'D'):
        provenance(out, args.part, wirings=list(ws), workers=args.workers, stage=args.stage if args.part == 'D' else None)
    if args.part == '0':
        part_0(out, ws)
    elif args.part == 'A':
        part_a(out, ws)
    elif args.part == 'Acheck':
        part_a_check(out)
    elif args.part == 'B1':
        part_b1(out, ws)
    elif args.part == 'B2':
        part_b2(out, ws, args.workers)
    elif args.part == 'B3':
        part_b3(out, ws)
    elif args.part == 'co':
        coresponse91(out, min(args.workers, 4))
    elif args.part == 'C':
        part_c(out, args.workers)
    elif args.part == 'score':
        score(out, min(args.workers, 4), coresponse=not args.no_co)
    elif args.part == 'D':
        part_d(out, args.stage, force=args.force, dry_run=args.dry_run)
    elif args.part == 'Dscore':
        score_d(out)
        verdict91(out)
        comparison_table(out)
    elif args.part == 'E':
        part_e(out)
    elif args.part == 'variants':
        write_variants(out)
    elif args.part == 'spread':
        spread91(out)
    elif args.part == 'table':
        figure(out)
        comparison_table(out)


if __name__ == '__main__':
    main()
