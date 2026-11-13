#!/usr/bin/env python3
"""E88: the error common to all wirings split into who responds and how strongly.

With sub-threshold responses set to 0, the similarity of two odors (one minus their cosine distance) is the product of
three factors (experiment_priorities.md, E88; planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-13 00:45):
- R, recruitment: the geometric mean of the two odors' response fractions;
- L, co-response: how many times more often than chance a KC that responds to one odor responds to the other
  (L = O / R, O one minus the distance of the binarized responses);
- G, the amplitude factor: G = s / O, the mean product of the shared responders' two responses over the root mean
  squares of the two odors' responses over their responders.
So s = R * L * G exactly. For the flies each of s, O and R is the mean over a genotype's recordings, and L and G follow;
for a model population the means are over its wirings and pseudo-recordings. The per-pair APL effect is s with APL off
minus s with APL on. A cell is that effect with the named factors taken from the flies and the others from the model,
in both conditions.

Readings (both sides read like for like):
- 'r1': trial means. The flies' trial means in z at 1 z; the model under E73's per-KC reading (`H.session_zscore`,
  saturating map, calibrated on each run's APL-on counts by E74's `_calibrate_zcell`) on N_PSEUDO pseudo-recordings of
  N_TRIALS trials through its measurement noise (numpy default_rng(SEED), the APL-on draws first: E74's `pseudo_noise`),
  the criterion 1 z on the trial means.
- 'r1b': as 'r1' with each condition's measurement noise multiplied until an odor's responder sets agree across repeats
  (the mean binary cosine between two different repeats of one odor, at the single-trial criterion that keeps the
  condition's 'r1' response rate) as the flies' do in that condition, and each condition's trial-mean criterion re-set
  to keep 'r1''s response rates.
- 'r2': the trial noise removed by the same estimator on both sides (`H.noise_free_pair_distances`'s, graded and
  binarized): the flies' single trials at their single-trial criterion (E87's `single_trial_criterion`, 1.057 z), the
  model's pseudo-recordings' single trials at the criterion at which its APL-on single-trial rate is its 'r1' APL-on
  rate (both conditions). R from the single-trial response fractions. 'r2b' (beside, and the estimator check): the
  same at 'r1b''s noise, each condition's single-trial criterion at 'r1''s rate in that condition.
- 'r3raw' / 'r3cal' (E87's convention, beside): the deterministic model, raw spike counts (a response two or more
  spikes) or E65's calibrated readout, against the flies with the trial noise removed ('r2''s flies).

Parts (CLI `--part`):
- A (no simulations): A0 the estimators on known truth (the model's noise-free session z-scores at the criterion that
  gives the reading's response rate in each condition); A1 the factors R, L, G, O and s per pair, condition and reading
  for the flies (all recordings, the six typical, each recording; 90 % intervals over N_BOOT resamples of the
  recordings; also in baseline-noise units), each wiring and each population; A2 the cells against the flies' effect
  in the same reading, with floors; A3 the same on E87's D2 runs (E80's per-odor recruitment on S, three arms) against
  the unchanged fits of S; verdict_A.json by the thresholds below.
- B (no simulations; readings 1 and 1b): B1 the changes of L and R under APL per group; B2 the side of a difference in
  the co-response (a per-odor cut of the responses without APL, and the flies' chosen form, against what APL leaves);
  B3 the shared responders' amplitudes; B4 the model's KCs that carry a group's overlap; verdict_B.json.
- C (no simulations): every archived variant (E87's `archive_variants`, E86's excepted), E85's 72 model flies and both
  families by the factors, under reading 1 (each variant's own calibration) and reading 3; the stand-out rule;
  verdict_C.json.
- D1 / D2 (simulations): not implemented here (stubs).
- table: comparison_table.md from the parts' tables.

Choices where the entry leaves the form open (stated here, and in the functions):
- A0 measures the error of s and O as the plain difference (estimate minus truth, both similarities) and the error of
  the multiplicative factors R, L and G in the units of the similarity they multiply into, s_true * ln(est / true) (the
  first-order error the factor's misestimate puts into s = R * L * G), so that EST_TOL (0.03) means the same for every
  factor; the plain differences are written beside.
- B3's top share uses, in both conditions, each odor's share of its APL-off responders that a per-odor cut to its APL-on
  fraction keeps (f_on / f_off).
- B4's "responders" are a wiring's KCs that respond (two or more spikes) to at least one odor with APL off; "the KCs a
  cut would keep and APL removes" are counted once per wiring and group.
- C's "handle the present code has" is read as: the variant's experiment has a script in experiments/ (the judgement of
  whether a handle composes with D2 is the session's).
Outputs: results/E88_responder_sets/ (comparison_table.md first); large outputs under outputs/E88_responder_sets/.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E87_inhibition_form as E87

EXPERIMENT = 'E88_responder_sets'

# the readings' settings (E74's pseudo-recordings)
N_PSEUDO, N_TRIALS, SEED = 10, 3, 43
N_BOOT = 1000                  # resamples of the recordings for the flies' intervals and the cells' floors
INTERVAL = (5.0, 95.0)         # the flies' 90 % intervals
B3_MIN_SHARED = 5              # B3: a pair's rank correlation needs at least this many shared responders
RAW_CRIT = E87.RAW_CRIT        # a model response is two or more spikes
MULT_BOUNDS = (0.1, 64.0)      # 1b: the search range of the noise multiplier

# the thresholds, fixed before any run (the entry's Design)
EST_TOL = 0.03                 # A0: a factor whose group mean is off by more than this in a reading is not read there
G_MIN = 0.015                  # A: g_sets (g_amp) at least this ...
G_RATIO = 0.5                  # ... and the other at most this share of it
GROUP_RATIO = 0.5              # A: a group's |centred error| with one factor from the flies at most this share of the own
D2_BASE_GAIN = 0.010           # A3: what is left on the base with the flies' recruitment is the co-response
D2_BASE_ARM = 'thrapl'         # the arm with the flies' recruitment in both conditions (per-odor APL scales)
B1_FACTOR = 1.10               # B1: a group's change differs by this factor, outside the flies' interval
B2_LARGE = 2.0 / 3.0           # B2: differ by at least two thirds of the difference under APL
B2_SMALL = 1.0 / 3.0           # B2: agree to within a third of it
B4_TOP = 0.25                  # B4: the top quarter of the responders' breadth weights ...
B4_SHARE = 2.0 / 3.0           # ... holds at least two thirds of the KCs a cut would keep and APL removes
C_BIN_GAIN = 0.008             # C: a variant's binarized pattern part lower than its baseline's by this
C_HALFWAY = 0.5                # C: or the alcohol group's co-response change at least halfway to the flies'

CONDS = ('on', 'off')
GROUPS = ('nearest4', 'ester', 'alcohol', 'six', 'other')
VERDICT_GROUPS = ('ester', 'alcohol')
FACTORS = ('R', 'L', 'G', 'O', 's')
MULT_FACTORS = ('R', 'L', 'G')
CELLS = {'own': (), 'sets': ('R', 'L'), 'amp': ('G',), 'R': ('R',), 'L': ('L',), 'flies': ('R', 'L', 'G')}
MODEL_READINGS = ('r1', 'r1b', 'r2', 'r2b', 'r3raw', 'r3cal')
NOISY_READINGS = ('r1', 'r1b', 'r2', 'r2b')
CELL_READINGS = ('r1', 'r1b', 'r2', 'r3raw', 'r3cal')
FLIES_READING = {'r1': 'r1', 'r1b': 'r1', 'r2': 'r2', 'r2b': 'r2', 'r3raw': 'r2', 'r3cal': 'r2'}
VERDICT_READINGS = ('r1', 'r2')
B_READINGS = ('r1', 'r1b')
POPULATIONS = ('b2', 'fc', 'S_unchanged', 'D2_thr', 'D2_throff', 'D2_thrapl')
SUBSETS = ('all', 'six')

_CACHE: Dict[tuple, object] = {}


def _iu() -> Tuple[np.ndarray, np.ndarray]:
    return E87._iu()


def n_pairs() -> int:
    return len(_iu()[0])


def pair_keys() -> List[Tuple[str, str]]:
    """The 120 pairs in triu order of APL16_ODORS, each as a canonical (sorted) pair."""
    i, j = _iu()
    o = H.APL16_ODORS
    return [tuple(sorted((o[a], o[b]))) for a, b in zip(i, j)]


def _mean_rows(a) -> np.ndarray:
    return E87._nanmean_rows(np.asarray(a, dtype=float))


def _nanmean(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(x.mean()) if x.size else float('nan')


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


# ---------------------------------------------------------------------------
# the factors

@dataclass
class Factors:
    """One condition's s, O and R per pair (120-vectors in triu order of APL16_ODORS); L = O / R and G = s / O."""
    s: np.ndarray
    O: np.ndarray
    R: np.ndarray

    @property
    def L(self) -> np.ndarray:
        with np.errstate(divide='ignore', invalid='ignore'):
            return self.O / self.R

    @property
    def G(self) -> np.ndarray:
        with np.errstate(divide='ignore', invalid='ignore'):
            return self.s / self.O

    def get(self, name: str) -> np.ndarray:
        return getattr(self, name)


@dataclass
class RecFactors:
    """One condition's per-recording s, O and R (recordings x 120; R NaN where O is, so that the three means are over
    the same recordings), with `mean(rows)` -> Factors."""
    s: np.ndarray
    O: np.ndarray
    R: np.ndarray

    def mean(self, rows: Optional[np.ndarray] = None) -> Factors:
        if rows is None:
            return Factors(_mean_rows(self.s), _mean_rows(self.O), _mean_rows(self.R))
        return Factors(_mean_rows(self.s[rows]), _mean_rows(self.O[rows]), _mean_rows(self.R[rows]))

    @property
    def n(self) -> int:
        return int(self.s.shape[0])

    def subset(self, rows) -> 'RecFactors':
        rows = np.asarray(rows, dtype=int)
        return RecFactors(self.s[rows], self.O[rows], self.R[rows])


def group_masks() -> Dict[str, np.ndarray]:
    """Boolean masks over the 120 pairs in triu order: 'nearest4' (E77's NEAREST), 'ester' (15), 'alcohol' (6), 'six'
    (E78's SIX), 'other' (the 99 pairs outside the ester and alcohol groups), 'all'."""
    if 'masks' not in _CACHE:
        from experiments import E77_input_side_common_error as E77
        from experiments import E78_claw_knee as E78
        keys = pair_keys()
        est = np.array([a in E77.ESTER and b in E77.ESTER for a, b in keys])
        alc = np.array([a in E77.ALCOHOL and b in E77.ALCOHOL for a, b in keys])
        near = {tuple(sorted(q)) for q in E77.NEAREST}
        six = {tuple(sorted(q)) for q in E78.SIX}
        _CACHE['masks'] = dict(nearest4=np.array([k in near for k in keys]), ester=est, alcohol=alc,
                               six=np.array([k in six for k in keys]), other=~(est | alc), all=np.ones(len(keys), dtype=bool))
    return _CACHE['masks']


def _sOR(s: np.ndarray, O: np.ndarray, R: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The three set to NaN together where any is (an odor without responders)."""
    bad = ~(np.isfinite(s) & np.isfinite(O) & np.isfinite(R))
    return tuple(np.where(bad, np.nan, v) for v in (s, O, R))


def factors_of_matrix(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(s, O, R) of one KC x odor matrix whose sub-threshold values are already 0 (triu order): s one minus the cosine
    distance, O one minus the distance of the binarized matrix, R the geometric mean of the two response fractions."""
    x = np.asarray(x, dtype=float)
    s = 1.0 - E87.cos_dist(x)
    b = (x > 0).astype(float)
    O = 1.0 - E87.cos_dist(b)
    f = b.mean(axis=0)
    R = np.sqrt(np.outer(f, f))[_iu()]
    return _sOR(s, O, R)


def factors_of_trials(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(s, O, R, reliability) of one recording's single trials (repeats, KCs, odors) whose sub-threshold values are
    already 0, with the trial noise removed (`H.noise_free_pair_distances`'s estimator: the mean single-trial cosine
    between different repeats over the square root of the two odors' own, graded and binarized); R from the single-trial
    response fractions (mean over repeats); reliability per odor (the mean binary cosine between two different
    repeats of the odor)."""
    x = np.asarray(x, dtype=float)
    mask = H._repeat_masks(x.shape[0])['different_repeat']
    iu = _iu()
    s = 1.0 - H._noise_free_matrix(H.trial_cosines(x), mask)[iu]
    b = (x > 0).astype(float)
    cb = H.trial_cosines(b)
    O = 1.0 - H._noise_free_matrix(cb, mask)[iu]
    f = b.mean(axis=(0, 1))
    R = np.sqrt(np.outer(f, f))[iu]
    rel = np.diag(H._kind_means(cb, mask)).copy()
    return (*_sOR(s, O, R), rel)


def rec_factors(mats: Sequence[np.ndarray]) -> RecFactors:
    """RecFactors of several thresholded matrices (one row each)."""
    fs = [factors_of_matrix(m) for m in mats]
    return RecFactors(*(np.stack([f[k] for f in fs]) for k in range(3)))


def rec_factors_trials(trials: Sequence[np.ndarray]) -> Tuple[RecFactors, np.ndarray]:
    """(RecFactors, reliability: recordings x odors) of several recordings' thresholded single trials."""
    fs = [factors_of_trials(t) for t in trials]
    return RecFactors(*(np.stack([f[k] for f in fs]) for k in range(3))), np.stack([f[3] for f in fs])


def cell_effect(model: Dict[str, Factors], flies: Dict[str, Factors], take: Sequence[str]) -> np.ndarray:
    """The per-pair effect (s with APL off minus s with APL on) of the product R * L * G with the factors in `take`
    from the flies and the others from the model, in both conditions."""
    def prod(c: str) -> np.ndarray:
        v = {f: (flies if f in take else model)[c].get(f) for f in MULT_FACTORS}
        # a pair with no shared responder (O = 0, so L = 0 and G = 0 / 0) has similarity 0 whatever G is
        with np.errstate(invalid='ignore'):
            return np.where((v['L'] == 0) | (v['R'] == 0), 0.0, v['R'] * v['L'] * v['G'])
    return prod('off') - prod('on')


def effect(F: Dict[str, Factors], factor: str = 's') -> np.ndarray:
    """A side's own per-pair effect on `factor` ('s': the graded effect; 'O': the binarized one): APL off minus on."""
    return F['off'].get(factor) - F['on'].get(factor)


def scores88(pred: np.ndarray, data: np.ndarray) -> dict:
    """E87's `scores` (loss, mean, pattern part, ρ, r² and slope, the ester and alcohol groups' centred errors, the
    four nearest and six near pairs' mean effect, the block) plus every group's mean centred error (err_<group>)."""
    s = E87.scores(pred, data)
    p, d = np.asarray(pred, dtype=float), np.asarray(data, dtype=float)
    ok = np.isfinite(p) & np.isfinite(d)
    ec = np.full(len(p), np.nan)
    ec[ok] = (p[ok] - p[ok].mean()) - (d[ok] - d[ok].mean())
    for g, m in group_masks().items():
        if g != 'all':
            s[f'err_{g}'] = _nanmean(ec[m])
    return s


def log_changes(F: Dict[str, Factors], factors: Sequence[str] = ('L', 'R', 'G', 'O', 's')) -> Dict[str, float]:
    """Per factor and group the mean over the group's pairs of log(factor on / factor off) ('chg_<factor>_<group>';
    pairs where the ratio is not positive and finite left out)."""
    out = {}
    for f in factors:
        with np.errstate(divide='ignore', invalid='ignore'):
            lr = np.log(F['on'].get(f) / F['off'].get(f))
        lr = np.where(np.isfinite(lr), lr, np.nan)
        for g, m in group_masks().items():
            out[f'chg_{f}_{g}'] = _nanmean(lr[m])
    return out


# ---------------------------------------------------------------------------
# the flies

def six_typical(unit: str = 'z') -> List[tuple]:
    """The six typical TNT recordings: the eight without E77's STRONG_TNT."""
    from experiments import E77_input_side_common_error as E77
    strong = {tuple(k) for k in E77.STRONG_TNT}
    keep = [k for k in E87.trialmean_mats(unit)['TNT'] if tuple(k) not in strong]
    if len(keep) != 6:
        raise RuntimeError(f'expected six typical TNT recordings, found {len(keep)}')
    return keep


def _flies_all(unit: str, reading: str) -> Dict[str, dict]:
    """{'on' / 'off': dict(keys, rf: RecFactors, rel: recordings x odors or None)} of every recording, cached."""
    key = ('flies', unit, reading)
    if key not in _CACHE:
        out = {}
        for c, g in (('on', 'TNT3c'), ('off', 'TNT')):
            if reading == 'r1':
                tm = E87.trialmean_mats(unit)[g]
                crit = E87.criterion(unit)
                keys = list(tm)
                out[c] = dict(keys=keys, rf=rec_factors([E87.thresholded(tm[k], crit) for k in keys]), rel=None)
            elif reading == 'r2':
                tr = E87.trial_arrays(unit)[g]
                c1 = E87.single_trial_criterion(unit)
                keys = list(tr)
                rf, rel = rec_factors_trials([E87.thresholded(tr[k], c1) for k in keys])
                out[c] = dict(keys=keys, rf=rf, rel=rel)
            else:
                raise ValueError(reading)
        _CACHE[key] = out
    return _CACHE[key]


def flies_factors(unit: str, reading: str, tnt: Optional[Sequence[tuple]] = None) -> Dict[str, RecFactors]:
    """{'on': control flies, 'off': APL output blocked (TNT)} per-recording factors in `unit` ('z' / 'noise') and
    reading ('r1' trial means at the criterion, 'r2' single trials at the single-trial criterion with the trial noise
    removed); `tnt` restricts the TNT recordings (default all)."""
    a = _flies_all(unit, reading)
    off = a['off']['rf']
    if tnt is not None:
        want = [tuple(k) for k in tnt]
        off = off.subset([a['off']['keys'].index(k) for k in want])
    return {'on': a['on']['rf'], 'off': off}


def flies_subset(subset: str, unit: str = 'z') -> Optional[List[tuple]]:
    return None if subset == 'all' else six_typical(unit)


def flies_reliability(unit: str = 'z') -> Dict[str, float]:
    """The flies' mean binary cosine between two different repeats of one odor at the single-trial criterion (mean
    over odors, then over recordings), per condition: the target of reading 1b."""
    a = _flies_all(unit, 'r2')
    out = {}
    for c in CONDS:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            out[c] = float(np.mean(np.nanmean(a[c]['rel'], axis=1)))
    return out


def boot_rows(n_on: int, n_off: int, n: int = N_BOOT, seed: int = 0) -> List[Tuple[np.ndarray, np.ndarray]]:
    """`n` resamples of the recordings with replacement within each genotype: [(control rows, TNT rows)]."""
    rng = np.random.default_rng(seed)
    return [(rng.integers(0, n_on, n_on), rng.integers(0, n_off, n_off)) for _ in range(n)]


def flies_boot(rf: Dict[str, RecFactors], n: int = N_BOOT, seed: int = 0) -> List[Dict[str, Factors]]:
    """The flies' factors over `n` resamples of the recordings."""
    return [{'on': rf['on'].mean(a), 'off': rf['off'].mean(b)} for a, b in boot_rows(rf['on'].n, rf['off'].n, n, seed)]


# ---------------------------------------------------------------------------
# the models

def _quantile_crit(x: np.ndarray, rate: float) -> float:
    """The criterion at which the share of `x` above it is `rate` (np.quantile at 1 - rate)."""
    return float(np.quantile(np.asarray(x, dtype=float).ravel(), 1.0 - rate))


class ModelReading:
    """One run's APL-on and APL-off spike counts read in every reading of MODEL_READINGS (lazily), with the
    noise-free truth of each noisy reading: the session z-scores without measurement noise at the criterion that gives
    the reading's response rate in that condition.

    Args:
        o: {'on', 'off'}: (responses or None, spike counts).
        label: a name.
        cal: E73's calibration of the reading (K, noise, noise_ratio, q, drift); default E74's `_calibrate_zcell` on the
            run's APL-on counts.
        K: the calibrated readout's K for 'r3cal' (default E83's `readout_K` on the run's APL-on counts).
        rel_target: the flies' reliability per condition (`flies_reliability`), needed for 'r1b' / 'r2b'.
        noise_scale: a factor on the measurement noise (1; 0 for the estimators' check at zero noise).
    """

    def __init__(self, o: Dict[str, tuple], label: str = '', cal: Optional[dict] = None, K: Optional[float] = None,
                 rel_target: Optional[Dict[str, float]] = None, noise_scale: float = 1.0):
        from experiments import E74_apl_effect_by_odors_per_kc as E74
        self.label, self.o, self.K, self.rel_target = label, o, K, rel_target
        self.counts = {c: E87._odor_matrix(o[c][1]) for c in CONDS}
        self.cal = dict(cal) if cal is not None else E74._calibrate_zcell(self.counts['on'])
        Kz = self.cal.get('K')
        kw = dict(q=float(self.cal['q']), noise=float(self.cal['noise']), amp=self.cal.get('amp', 'saturating'),
                  K=None if Kz is None or not np.isfinite(Kz) else float(Kz), noise_ratio=float(self.cal['noise_ratio']),
                  drift=float(self.cal['drift']))
        rng = np.random.default_rng(SEED)
        self.C, self.sd, self.xi, self.xim = {}, {}, {}, {}
        for c in CONDS:
            z = H.session_zscore(self.counts[c], **kw)
            self.C[c] = z['center']                                    # (1, KCs, odors)
            self.sd[c] = z['trial_sd'] if noise_scale == 1.0 else z['trial_sd'] * noise_scale   # (1, KCs)
            self.xi[c] = rng.standard_normal((N_PSEUDO, N_TRIALS) + self.C[c].shape[1:])
            self.xim[c] = self.xi[c].mean(axis=1)                      # E74's pseudo_noise, the same draws
        self._f: Dict[str, Dict[str, RecFactors]] = {}
        self._tm: Dict[str, Dict[str, np.ndarray]] = {}
        self._truth: Dict[str, Dict[str, Factors]] = {}
        self.crit: Dict[str, Dict[str, float]] = {}
        self.rate: Dict[str, Dict[str, float]] = {}
        self.rel: Dict[str, Dict[str, float]] = {}
        self.mult: Optional[Dict[str, float]] = None
        self.mult_info: Dict[str, dict] = {}

    # the pseudo-recordings
    def trial_means(self, c: str, m: float = 1.0) -> np.ndarray:
        """(N_PSEUDO, KCs, odors) trial means in z through the measurement noise times `m` (E74's `read_zcell` at m 1)."""
        sd = self.sd[c] if m == 1.0 else m * self.sd[c]
        return self.C[c] + sd[:, :, None] * self.xim[c]

    def trials(self, c: str, m: float = 1.0) -> np.ndarray:
        """(N_PSEUDO, N_TRIALS, KCs, odors) single trials in z through the measurement noise times `m`."""
        sd = self.sd[c] if m == 1.0 else m * self.sd[c]
        return self.C[c][:, None] + sd[:, None, :, None] * self.xi[c]

    def reliability(self, c: str, m: float, rate: float) -> float:
        """The mean binary cosine between two different repeats of one odor (mean over odors, then pseudo-recordings)
        at noise `m` and the single-trial criterion that gives `rate`."""
        tr = self.trials(c, m)
        crit = _quantile_crit(tr, rate)
        mask = H._repeat_masks(N_TRIALS)['different_repeat']
        vals = []
        for r in range(tr.shape[0]):
            d = np.diag(H._kind_means(H.trial_cosines((tr[r] > crit).astype(float)), mask))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                vals.append(float(np.nanmean(d)))
        return float(np.mean(vals))

    def multipliers(self) -> Dict[str, float]:
        """Reading 1b's noise multiplier per condition: geometric bisection in MULT_BOUNDS on the reliability at the
        condition's 'r1' rate against the flies' (`rel_target`)."""
        if self.mult is None:
            if self.rel_target is None:
                raise ValueError(f'{self.label}: reading 1b needs rel_target')
            self.factors('r1')
            self.mult = {}
            for c in CONDS:
                rate, tgt = self.rate['r1'][c], float(self.rel_target[c])
                lo, hi = MULT_BOUNDS
                r_lo, r_hi = self.reliability(c, lo, rate), self.reliability(c, hi, rate)
                at = ''
                if r_lo <= tgt:
                    m, at = lo, 'low'
                elif r_hi >= tgt:
                    m, at = hi, 'high'
                else:
                    for _ in range(30):
                        mid = float(np.sqrt(lo * hi))
                        if self.reliability(c, mid, rate) > tgt:
                            lo = mid
                        else:
                            hi = mid
                    m = float(np.sqrt(lo * hi))
                self.mult[c] = m
                self.mult_info[c] = dict(mult=m, at_bound=at, rel_target=tgt, rel_r1=self.reliability(c, 1.0, rate),
                                         rel_r1b=self.reliability(c, m, rate))
        return self.mult

    def factors(self, reading: str) -> Dict[str, RecFactors]:
        """{'on', 'off'}: RecFactors (pseudo-recordings as rows; one row for the deterministic readings)."""
        if reading in self._f:
            return self._f[reading]
        crit, rate, f, tms = {}, {}, {}, {}
        if reading == 'r1':
            for c in CONDS:
                tm = self.trial_means(c)
                crit[c] = float(H.KC_RESPONSE_Z)
                rate[c] = float((tm > crit[c]).mean())
                tms[c] = E87.thresholded(tm, crit[c])
                f[c] = rec_factors(list(tms[c]))
        elif reading == 'r2':
            self.factors('r1')
            c1 = _quantile_crit(self.trials('on'), self.rate['r1']['on'])
            rel = {}
            for c in CONDS:
                tr = self.trials(c)
                crit[c] = c1
                rate[c] = float((tr > c1).mean())
                f[c], rr = rec_factors_trials([E87.thresholded(t, c1) for t in tr])
                rel[c] = float(np.mean(np.nanmean(rr, axis=1)))
            self.rel['r2'] = rel
        elif reading in ('r1b', 'r2b'):
            self.factors('r1')
            m = self.multipliers()
            rel = {}
            for c in CONDS:
                r1rate = self.rate['r1'][c]
                if reading == 'r1b':
                    tm = self.trial_means(c, m[c])
                    crit[c] = _quantile_crit(tm, r1rate)
                    rate[c] = float((tm > crit[c]).mean())
                    tms[c] = E87.thresholded(tm, crit[c])
                    f[c] = rec_factors(list(tms[c]))
                else:
                    tr = self.trials(c, m[c])
                    crit[c] = _quantile_crit(tr, r1rate)
                    rate[c] = float((tr > crit[c]).mean())
                    f[c], rr = rec_factors_trials([E87.thresholded(t, crit[c]) for t in tr])
                    rel[c] = float(np.mean(np.nanmean(rr, axis=1)))
            if rel:
                self.rel[reading] = rel
        elif reading in ('r3raw', 'r3cal'):
            K = None
            if reading == 'r3cal':
                if self.K is None:
                    from experiments import E83_cross_wiring_fit as E
                    self.K = E.readout_K(self.o)
                K = float(self.K)
            for c in CONDS:
                x = self.counts[c] if K is None else self.counts[c] / (self.counts[c] + K)
                crit[c] = RAW_CRIT if K is None else RAW_CRIT / (RAW_CRIT + K)
                rate[c] = float((x > crit[c]).mean())
                tms[c] = E87.thresholded(x, crit[c])[None]
                f[c] = rec_factors(list(tms[c]))
        else:
            raise ValueError(reading)
        self.crit[reading], self.rate[reading], self._f[reading] = crit, rate, f
        if tms:
            self._tm[reading] = tms
        return f

    def thresholded_means(self, reading: str) -> Dict[str, np.ndarray]:
        """{'on', 'off'}: (pseudo-recordings, KCs, odors) thresholded trial means of reading 'r1' / 'r1b' (one row of
        thresholded values for the deterministic readings)."""
        self.factors(reading)
        return self._tm[reading]

    def mean_factors(self, reading: str) -> Dict[str, Factors]:
        return {c: f.mean() for c, f in self.factors(reading).items()}

    def truth(self, reading: str) -> Dict[str, Factors]:
        """The noise-free truth of a noisy reading: the factors of the session z-scores without measurement noise at
        the criterion that gives the reading's response rate in each condition."""
        if reading not in self._truth:
            self.factors(reading)
            out = {}
            for c in CONDS:
                C0 = self.C[c][0]
                ct = _quantile_crit(C0, self.rate[reading][c])
                out[c] = rec_factors([E87.thresholded(C0, ct)]).mean()
            self._truth[reading] = out
        return self._truth[reading]


def _pop_mean(fs: Sequence[Dict[str, Factors]]) -> Dict[str, Factors]:
    return {c: Factors(*(_mean_rows(np.stack([getattr(f[c], k) for f in fs])) for k in ('s', 'O', 'R'))) for c in CONDS}


def population(members: Sequence[ModelReading], reading: str) -> Dict[str, Factors]:
    """A population's factors in a reading: s, O and R the means over its wirings (each the mean over its
    pseudo-recordings)."""
    return _pop_mean([m.mean_factors(reading) for m in members])


def population_truth(members: Sequence[ModelReading], reading: str) -> Dict[str, Factors]:
    return _pop_mean([m.truth(reading) for m in members])


def model_populations() -> Dict[str, List[Tuple[str, tuple, Dict[str, tuple], Optional[float]]]]:
    """{population: [(wiring, cache key, outputs, K of the calibrated readout or None)]}: 'b2' (E67's fit on E72 part
    B's eight wirings), 'fc' (the firing-ceiling family), 'S_unchanged' (b2 on S), 'D2_<arm>' (E87's D2 runs on S's
    wirings; K the fit's, as E87's two-by-two on D2)."""
    from experiments import E83_cross_wiring_fit as E
    fams = E87.model_families()
    out = {}
    for fam in ('b2', 'fc'):
        out[fam] = [(w, (f, t), E87._load(f, t), None) for w, f, t in fams[fam]]
    b2 = {w: (k, o) for w, k, o, _ in out['b2']}
    out['S_unchanged'] = [(w, b2[w][0], b2[w][1], None) for w in E.S_WIRINGS]
    for arm in E87.D2_ARMS:
        mem = []
        for w in E.S_WIRINGS:
            tag = f'D2_b3_{arm}_{E87.D2_UNITS}_{w}'
            o = E87._load_d2(tag)
            if o is None:
                raise RuntimeError(f'no D2 run {tag} under {E87._d_dir()}')
            mem.append((w, ('E87_inhibition_form', tag), o, float(E.saved_fit(w)['K'])))
        out[f'D2_{arm}'] = mem
    return out


def build_members(rel_target: Dict[str, float]) -> Tuple[Dict[str, List[ModelReading]], Dict[str, list]]:
    """({population: [ModelReading]}, {population: [(wiring, folder, tag)]}); runs shared by populations are read once."""
    pops = model_populations()
    cache: Dict[tuple, ModelReading] = {}
    out, ids = {}, {}
    for pop, mem in pops.items():
        lst = []
        for w, key, o, K in mem:
            if key not in cache:
                cache[key] = ModelReading(o, label=w, K=K, rel_target=rel_target)
            lst.append(cache[key])
        out[pop] = lst
        ids[pop] = [(w, key[0], key[1]) for w, key, _, _ in mem]
    return out, ids


def _member_info(members: Dict[str, List[ModelReading]], ids: Dict[str, list]) -> pd.DataFrame:
    """Per population and wiring the calibration, the readings' criteria, rates, reliabilities and 1b's multipliers."""
    rows = []
    for pop, lst in members.items():
        for mr, (w, folder, tag) in zip(lst, ids[pop]):
            row = dict(population=pop, wiring=w, folder=folder, tag=tag, **{f'cal_{k}': v for k, v in mr.cal.items()},
                       K_readout=mr.K)
            for rd in mr.crit:
                for c in CONDS:
                    row[f'{rd}_crit_{c}'] = mr.crit[rd][c]
                    row[f'{rd}_rate_{c}'] = mr.rate[rd][c]
            for rd, v in mr.rel.items():
                for c in CONDS:
                    row[f'{rd}_rel_{c}'] = v[c]
            for c, v in mr.mult_info.items():
                row.update({f'mult_{c}_{k}': x for k, x in v.items()})
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# part A

def _group_row(vals: Dict[str, np.ndarray]) -> Dict[str, float]:
    return {f'{k}_{g}': _nanmean(v[m]) for k, v in vals.items() for g, m in group_masks().items()}


def estimator_check(members: Sequence[ModelReading], readings: Sequence[str] = NOISY_READINGS) -> pd.DataFrame:
    """Each factor's estimate (the population's) minus its truth by group and condition: 'diff' (estimate minus truth,
    group mean over pairs) and 'err' (s and O: diff; R, L and G: s_true * ln(est / true), the first-order error in s),
    'off' if |err| > EST_TOL."""
    rows = []
    masks = group_masks()
    for rd in readings:
        est, tru = population(members, rd), population_truth(members, rd)
        for c in CONDS:
            st = tru[c].s
            for f in FACTORS:
                e, t = est[c].get(f), tru[c].get(f)
                diff = e - t
                if f in MULT_FACTORS:
                    with np.errstate(divide='ignore', invalid='ignore'):
                        err = st * np.log(e / t)
                else:
                    err = diff
                for g in GROUPS + ('all',):
                    m = masks[g]
                    er = _nanmean(err[m])
                    rows.append(dict(reading=rd, condition=c, factor=f, group=g, n_pairs=int(m.sum()), estimate=_nanmean(e[m]),
                                     truth=_nanmean(t[m]), diff=_nanmean(diff[m]), err=er,
                                     off=bool(g != 'all' and np.isfinite(er) and abs(er) > EST_TOL)))
    return pd.DataFrame(rows)


def factors_read(ec: pd.DataFrame) -> Dict[str, Dict[str, bool]]:
    """{reading: {factor: read}} from an estimator check: a factor is not read in a reading where any group's mean
    error in either condition is beyond EST_TOL."""
    out = {}
    for rd, g in ec.groupby('reading'):
        out[rd] = {f: not bool(g[(g.factor == f)].off.any()) for f in FACTORS}
    return out


def part_a0(out: Path, members: Dict[str, List[ModelReading]]) -> pd.DataFrame:
    """A0: the estimators on known truth for every population (the thresholds read the breadth-2 population),
    estimator_check.csv; and the check at zero noise (each noisy reading of the real wiring's breadth-2 fit with no
    measurement noise against its truth, the largest absolute difference: estimator_zero_noise.json)."""
    parts = []
    for pop, lst in members.items():
        ec = estimator_check(lst)
        ec.insert(0, 'population', pop)
        parts.append(ec)
    df = pd.concat(parts, ignore_index=True)
    df.to_csv(out / 'estimator_check.csv', index=False)
    mr0 = members['b2'][0]
    z = ModelReading(mr0.o, label='zero noise', cal=mr0.cal, rel_target=mr0.rel_target, noise_scale=0.0)
    z.mult = {'on': 1.0, 'off': 1.0}
    zero = {}
    for rd in NOISY_READINGS:
        est, tru = z.mean_factors(rd), z.truth(rd)
        zero[rd] = {f: float(np.nanmax([np.nanmax(np.abs(est[c].get(f) - tru[c].get(f))) for c in CONDS])) for f in FACTORS}
    _write_json(out / 'estimator_zero_noise.json', zero)
    print('A0 zero noise, largest |estimate - truth|:', {k: max(v.values()) for k, v in zero.items()}, flush=True)
    return df


def _factor_rows(F: Dict[str, Factors], boot: Optional[List[Dict[str, Factors]]] = None) -> List[dict]:
    """Per condition and group the mean of each factor over the group's pairs (and the 90 % interval over `boot`)."""
    rows = []
    masks = group_masks()
    for c in CONDS:
        for g in GROUPS + ('all',):
            m = masks[g]
            row = dict(condition=c, group=g)
            for f in FACTORS:
                row[f] = _nanmean(F[c].get(f)[m])
                if boot is not None:
                    v = [_nanmean(b[c].get(f)[m]) for b in boot]
                    row[f'{f}_lo'], row[f'{f}_hi'] = (float(x) for x in np.nanpercentile(v, INTERVAL))
            rows.append(row)
    return rows


def _pair_rows(F: Dict[str, Factors], boot: Optional[List[Dict[str, Factors]]] = None, **meta) -> pd.DataFrame:
    keys = pair_keys()
    masks = group_masks()
    parts = []
    for c in CONDS:
        d = dict(**meta, condition=c, odor_a=[k[0] for k in keys], odor_b=[k[1] for k in keys],
                 group=['ester' if masks['ester'][i] else 'alcohol' if masks['alcohol'][i] else 'other' for i in range(len(keys))],
                 nearest4=masks['nearest4'], six=masks['six'])
        for f in FACTORS:
            d[f] = F[c].get(f)
            if boot is not None:
                arr = np.stack([b[c].get(f) for b in boot])
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    d[f'{f}_lo'], d[f'{f}_hi'] = np.nanpercentile(arr, INTERVAL, axis=0)
        parts.append(pd.DataFrame(d))
    return pd.concat(parts, ignore_index=True)


def part_a1(out: Path, members: Dict[str, List[ModelReading]], ids: Dict[str, list]) -> pd.DataFrame:
    """A1: the factors per condition and group (factors.csv) and per pair (factors_pairs.csv): the flies in z and in
    baseline-noise units, readings 1 and 2, all recordings and the six typical (90 % intervals over N_BOOT resamples),
    and each recording; each population's wirings and the population in every model reading."""
    rows, prs = [], []
    for unit in E87.UNITS:
        for rd in ('r1', 'r2'):
            for subset in SUBSETS:
                rf = flies_factors(unit, rd, flies_subset(subset, unit))
                F = {c: rf[c].mean() for c in CONDS}
                boot = flies_boot(rf)
                meta = dict(source='flies', unit=unit, reading=rd, population='flies', member=f'mean_{subset}', subset=subset)
                rows += [dict(meta, **r) for r in _factor_rows(F, boot)]
                prs.append(_pair_rows(F, boot, **meta))
            a = _flies_all(unit, rd)
            for c, g in (('on', 'TNT3c'), ('off', 'TNT')):
                for n, k in enumerate(a[c]['keys']):
                    rf1 = a[c]['rf'].subset([n])
                    F1 = {cc: rf1.mean() for cc in CONDS}   # one recording: its own condition only
                    meta = dict(source='flies', unit=unit, reading=rd, population='flies', member=f'{g}_{k[0]}_{k[1]}_{k[2]}',
                                subset='recording')
                    rows += [dict(meta, **r) for r in _factor_rows(F1) if r['condition'] == c]
                    if unit == 'z':
                        p = _pair_rows(F1, **meta)
                        prs.append(p[p.condition == c])
    for pop, lst in members.items():
        for rd in MODEL_READINGS:
            meta = dict(source='model', unit='z' if rd in NOISY_READINGS else rd[2:], reading=rd, population=pop, subset='')
            for mr, (w, _, _) in zip(lst, ids[pop]):
                F = mr.mean_factors(rd)
                rows += [dict(meta, member=w, **r) for r in _factor_rows(F)]
            Fp = population(lst, rd)
            rows += [dict(meta, member='population', **r) for r in _factor_rows(Fp)]
            prs.append(_pair_rows(Fp, member='population', **meta))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'factors.csv', index=False)
    pd.concat(prs, ignore_index=True).to_csv(out / 'factors_pairs.csv', index=False)
    print(f'A1: {len(df)} rows', flush=True)
    return df


def _cell_rows(Fm: Dict[str, Factors], Ff: Dict[str, Factors], boot: List[Dict[str, Factors]], **meta) -> List[dict]:
    """The cells of one model population against the flies' effect in one reading, with floors (the flies' factors
    taken from resampled recordings, the rest as measured: the mean over resamples of the pattern part / loss of the
    cell minus the cell as measured); 'binary' the model's binarized effect against the flies' binarized effect (its
    floor the flies' binarized effect resampled)."""
    data = effect(Ff)
    rows = []
    for name, take in CELLS.items():
        pred = cell_effect(Fm, Ff, take)
        fl_p, fl_l = [], []
        if take:
            for b in boot:
                e = cell_effect(Fm, b, take) - pred
                ok = np.isfinite(e)
                fl_p.append(float(np.std(e[ok])))
                fl_l.append(float(np.sqrt(np.mean(e[ok] ** 2))))
        rows.append(dict(meta, cell=name, take='+'.join(take) if take else 'none', **scores88(pred, data),
                         floor_pattern=float(np.mean(fl_p)) if fl_p else 0.0, floor_loss=float(np.mean(fl_l)) if fl_l else 0.0))
    db = effect(Ff, 'O')
    pb = effect(Fm, 'O')
    fl_p = []
    for b in boot:
        e = effect(b, 'O') - db
        ok = np.isfinite(e)
        fl_p.append(float(np.std(e[ok])))
    rows.append(dict(meta, cell='binary', take='binarized', **scores88(pb, db), floor_pattern=float(np.mean(fl_p)), floor_loss=np.nan))
    return rows


def part_a2(out: Path, members: Dict[str, List[ModelReading]]) -> pd.DataFrame:
    """A2 and A3: every population's cells against the flies' effect in the same reading (CELL_READINGS; the flies in
    z, all recordings and the six typical), with floors. cells.csv."""
    rows = []
    for subset in SUBSETS:
        tnt = flies_subset(subset)
        fl = {}
        for frd in ('r1', 'r2'):
            rf = flies_factors('z', frd, tnt)
            fl[frd] = ({c: rf[c].mean() for c in CONDS}, flies_boot(rf))
        for pop, lst in members.items():
            for rd in CELL_READINGS:
                Ff, boot = fl[FLIES_READING[rd]]
                rows += _cell_rows(population(lst, rd), Ff, boot, population=pop, reading=rd, subset=subset)
        print(f'A2 {subset}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'cells.csv', index=False)
    return df


def main_target(out: Path, members: Dict[str, List[ModelReading]]) -> pd.DataFrame:
    """Beside: each population's own scores on the main target (E83's `effect_scores`, raw counts and the calibrated
    readout). main_target.csv."""
    from experiments import E83_cross_wiring_fit as E
    rows = []
    for pop, lst in members.items():
        for ro in ('raw', 'cal'):
            ds = []
            for mr in lst:
                if ro == 'cal' and mr.K is None:
                    mr.K = E.readout_K(mr.o)
                ds.append(E.pair_distances(mr.o, ro, mr.K if ro == 'cal' else None))
            rows.append(dict(population=pop, readout=ro, **E.effect_scores(E.population(ds))))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'main_target.csv', index=False)
    return df


def _split_verdict(own: float, sets: float, amp: float) -> dict:
    g_sets, g_amp = own - sets, own - amp
    if g_sets >= G_MIN and g_amp <= G_RATIO * g_sets:
        v = 'who responds'
    elif g_amp >= G_MIN and g_sets <= G_RATIO * g_amp:
        v = 'how strongly'
    else:
        v = 'in both'
    return dict(own=own, sets=sets, amp=amp, g_sets=g_sets, g_amp=g_amp, verdict=v)


def _group_verdict(own: float, r_err: float, l_err: float) -> dict:
    r_ok = abs(r_err) <= GROUP_RATIO * abs(own)
    l_ok = abs(l_err) <= GROUP_RATIO * abs(own)
    v = 'recruitment' if r_ok and not l_ok else 'co-response' if l_ok and not r_ok else 'both'
    return dict(own_err=own, R_err=r_err, L_err=l_err, verdict=v)


def _agree(vals: Sequence[str]) -> str:
    vals = list(vals)
    return vals[0] if all(v == vals[0] for v in vals) else 'the readings disagree'


def verdict_a(out: Path) -> dict:
    """The thresholds of A on the breadth-2 population of eight wirings in readings 1 and 2 (which must agree), from
    cells.csv and estimator_check.csv, every verdict repeated with the six typical recordings; reading 3 beside.
    verdict_A.json."""
    cells = pd.read_csv(out / 'cells.csv')
    ec = pd.read_csv(out / 'estimator_check.csv')
    read = factors_read(ec[ec.population == 'b2'])

    def cell(pop, rd, subset, name):
        r = cells[(cells.population == pop) & (cells.reading == rd) & (cells.subset == subset) & (cells.cell == name)]
        return r.iloc[0]

    v = dict(estimator=dict(tolerance=EST_TOL, factors_read=read,
                            worst={rd: {f: float(g[g.factor == f].err.abs().max()) for f in FACTORS}
                                   for rd, g in ec[(ec.population == 'b2') & (ec.group != 'all')].groupby('reading')}))
    split, groups, d2 = {}, {}, {}
    for subset in SUBSETS:
        sp = {}
        for rd in VERDICT_READINGS + ('r1b', 'r3raw', 'r3cal'):
            r = _split_verdict(*(float(cell('b2', rd, subset, n).pattern) for n in ('own', 'sets', 'amp')))
            if rd in read:
                r['factors_read'] = all(read[rd][f] for f in MULT_FACTORS)
            r['rho_own'] = float(cell('b2', rd, subset, 'own').rho)
            r['rho_sets'] = float(cell('b2', rd, subset, 'sets').rho)
            r['floor_sets'] = float(cell('b2', rd, subset, 'sets').floor_pattern)
            r['binary'] = float(cell('b2', rd, subset, 'binary').pattern)
            sp[rd] = r
        sp['verdict'] = _agree([sp[rd]['verdict'] for rd in VERDICT_READINGS])
        ok = all(sp[rd].get('factors_read', True) for rd in VERDICT_READINGS)
        sp['verdict_with_estimator_check'] = sp['verdict'] if ok else (
            'not read: ' + '; '.join(f'{rd}: ' + ', '.join(f for f in MULT_FACTORS if not read[rd][f])
                                     for rd in VERDICT_READINGS if not all(read[rd][f] for f in MULT_FACTORS)))
        split[subset] = sp
        gr = {}
        for g in GROUPS:
            gg = {}
            for rd in VERDICT_READINGS + ('r1b', 'r3raw', 'r3cal'):
                gg[rd] = _group_verdict(*(float(cell('b2', rd, subset, n)[f'err_{g}']) for n in ('own', 'R', 'L')))
                if rd in read:
                    gg[rd]['factors_read'] = bool(read[rd]['R'] and read[rd]['L'])
            gg['verdict'] = _agree([gg[rd]['verdict'] for rd in VERDICT_READINGS])
            gr[g] = gg
        groups[subset] = gr
        dd = {}
        for rd in VERDICT_READINGS + ('r1b', 'r3raw', 'r3cal'):
            own = float(cell(f'D2_{D2_BASE_ARM}', rd, subset, 'own').pattern)
            lc = float(cell(f'D2_{D2_BASE_ARM}', rd, subset, 'L').pattern)
            dd[rd] = dict(own=own, L_cell=lc, gain=own - lc, coresponse_left=bool(own - lc >= D2_BASE_GAIN),
                          sets_cell=float(cell(f'D2_{D2_BASE_ARM}', rd, subset, 'sets').pattern),
                          unchanged_own=float(cell('S_unchanged', rd, subset, 'own').pattern))
        dd['verdict'] = _agree(['the co-response is left' if dd[rd]['coresponse_left'] else 'not the co-response alone'
                                for rd in VERDICT_READINGS])
        d2[subset] = dd
    v.update(split=split, groups=groups, d2_base=dict(arm=D2_BASE_ARM, **d2))
    changed = []
    if split['six']['verdict'] != split['all']['verdict']:
        changed.append(f"split: {split['all']['verdict']} with all recordings, {split['six']['verdict']} with the six typical")
    for g in GROUPS:
        if groups['six'][g]['verdict'] != groups['all'][g]['verdict']:
            changed.append(f"{g}: {groups['all'][g]['verdict']} with all, {groups['six'][g]['verdict']} with the six typical")
    if d2['six']['verdict'] != d2['all']['verdict']:
        changed.append(f"D2 base: {d2['all']['verdict']} with all, {d2['six']['verdict']} with the six typical")
    v['depends_on_the_two_strongest_recordings'] = changed
    v['thresholds'] = dict(G_MIN=G_MIN, G_RATIO=G_RATIO, GROUP_RATIO=GROUP_RATIO, D2_BASE_GAIN=D2_BASE_GAIN, EST_TOL=EST_TOL)
    _write_json(out / 'verdict_A.json', v)
    return v


def figure_a(out: Path) -> Path:
    """E88_cells.png: the breadth-2 population's per-pair effect against the flies' (z, all recordings, reading 1 and
    reading 2): the model's own, the cell with the flies' responder sets and the cell with the flies' amplitudes."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fp = pd.read_csv(out / 'factors_pairs.csv')
    masks = group_masks()
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.4), sharex=True, sharey=True)
    for i, rd in enumerate(('r1', 'r2')):
        def get(source, pop, member, reading, unit, subset=''):
            g = fp[(fp.source == source) & (fp.population == pop) & (fp.member == member) & (fp.reading == reading) & (fp.unit == unit)]
            if subset:
                g = g[g.subset == subset]
            return {c: Factors(*(g[g.condition == c][k].to_numpy() for k in ('s', 'O', 'R'))) for c in CONDS}
        Ff = get('flies', 'flies', 'mean_all', rd, 'z', 'all')
        Fm = get('model', 'b2', 'population', rd, 'z')
        data = effect(Ff)
        for j, (name, take) in enumerate((("the model's own", ()), ("the flies' responder sets", ('R', 'L')),
                                          ("the flies' amplitudes", ('G',)))):
            ax = axes[i, j]
            pred = cell_effect(Fm, Ff, take)
            s = scores88(pred, data)
            lo, hi = -0.05, 0.45
            ax.plot([lo, hi], [lo, hi], color='0.6', lw=1)
            for g, col in (('other', '0.45'), ('ester', '#d1495b'), ('alcohol', '#2e86ab')):
                m = masks[g]
                ax.scatter(data[m], pred[m], s=14, color=col, label=f'{g} pairs')
            ax.set_title(f'{name}, {"trial means" if rd == "r1" else "trial noise removed"}\n'
                         f'(pattern part {s["pattern"]:.3f}, mean {s["mean"]:+.3f})', fontsize=9)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            if i == 1:
                ax.set_xlabel('KC data: per-pair APL effect (z, sub-threshold set to 0)')
        axes[i, 0].set_ylabel('breadth-2 population: per-pair APL effect')
    axes[0, 0].legend(fontsize=8, loc='upper left')
    fig.tight_layout()
    path = out / 'E88_cells.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _provenance(out: Path, part: str, **extra) -> None:
    _write_json(out / f'provenance_{part}.json', dict(script='experiments/E88_responder_sets.py', part=part, **H.code_commits(),
                                                     settings=dict(N_PSEUDO=N_PSEUDO, N_TRIALS=N_TRIALS, SEED=SEED, N_BOOT=N_BOOT,
                                                                   INTERVAL=INTERVAL, MULT_BOUNDS=MULT_BOUNDS),
                                                     **extra))


def _setup_members(out: Path) -> Tuple[Dict[str, List[ModelReading]], Dict[str, list], Dict[str, float]]:
    rel = flies_reliability('z')
    print('flies reliability (z, single-trial criterion):', rel, flush=True)
    t0 = time.time()
    members, ids = build_members(rel)
    for pop, lst in members.items():
        for mr in lst:
            for rd in MODEL_READINGS:
                mr.factors(rd)
        print(f'{pop}: read in {time.time() - t0:.0f}s', flush=True)
    return members, ids, rel


def part_a(out: Path) -> dict:
    """A0 to A3, the main target beside, A's verdict and the figure."""
    members, ids, rel = _setup_members(out)
    info = _member_info(members, ids)
    info.to_csv(out / 'members.csv', index=False)
    _provenance(out, 'A', flies_reliability=rel, criteria={u: E87.criterion(u) for u in E87.UNITS},
                single_trial_criteria={u: E87.single_trial_criterion(u) for u in E87.UNITS},
                populations={k: [list(x) for x in v] for k, v in ids.items()}, chosen_form=_chosen_form().name)
    part_a0(out, members)
    part_a1(out, members, ids)
    part_a2(out, members)
    main_target(out, members)
    v = verdict_a(out)
    figure_a(out)
    print('A verdict:', v['split']['all']['verdict'], '| six:', v['split']['six']['verdict'], flush=True)
    return v


# ---------------------------------------------------------------------------
# part B

def _chosen_form() -> 'E87.Form':
    return E87.chosen_form(H.results_dir('E87_inhibition_form'))


def _flies_cut_form(tnt: Optional[Sequence[tuple]]) -> Dict[str, Dict[str, Factors]]:
    """The flies (z, trial means): {'apl': {'on': control, 'off': TNT}, 'cut': {'on': the TNT recordings cut per odor to
    the control flies' fractions (E87's per-odor oracle, divisive), 'off': TNT}, 'form': {'on': the flies' chosen form on
    the TNT recordings at the control flies' rate, 'off': TNT}}."""
    rd = E87.Reading('z', 'trial_means', tnt=tnt)
    crit = rd.crit
    fc = np.mean([(m > crit).mean(axis=0) for m in rd.tm_ctrl], axis=0)
    off = rec_factors([E87.thresholded(m, crit) for m in rd.tm_tnt]).mean()
    on = rec_factors([E87.thresholded(m, crit) for m in rd.tm_ctrl]).mean()
    cut = rec_factors([E87.per_odor_oracle(m, fc, 'divisive', crit) for m in rd.tm_tnt]).mean()
    f_ch = _chosen_form()
    _, c = rd.predict(f_ch)
    form = rec_factors(rd.reduced(f_ch, c)).mean()
    return {'apl': {'on': on, 'off': off}, 'cut': {'on': cut, 'off': off}, 'form': {'on': form, 'off': off}}


def flies_fractions() -> np.ndarray:
    """The control flies' per-odor response fractions (z, trial means at 1 z; mean over recordings)."""
    rd = E87.Reading('z', 'trial_means')
    return np.mean([(m > rd.crit).mean(axis=0) for m in rd.tm_ctrl], axis=0)


def model_cut_form(mr: ModelReading, reading: str, fc: np.ndarray) -> Dict[str, Dict[str, Factors]]:
    """One wiring in reading 'r1' / 'r1b': {'apl', 'cut_own' (each odor's APL-off responders cut by amplitude to its
    APL-on fraction, the mean over the pseudo-recordings), 'cut_flies' (to the control flies' fractions), 'form' (the
    flies' chosen form on the APL-off responses, c over the pseudo-recordings at the wiring's APL-on rate)}: each
    {'on': the result, 'off': the APL-off responses}."""
    tms = mr.thresholded_means(reading)
    crit = mr.crit[reading]
    off_list = list(tms['off'])
    off = mr.mean_factors(reading)['off']
    f_on = np.mean([(m > 0).mean(axis=0) for m in tms['on']], axis=0)
    cut_own = rec_factors([E87.per_odor_oracle(m, f_on, 'divisive', crit['off']) for m in off_list]).mean()
    cut_fl = rec_factors([E87.per_odor_oracle(m, fc, 'divisive', crit['off']) for m in off_list]).mean()
    f_ch = _chosen_form()
    inhs = [E87.inhibition(m, f_ch, crit['off']) for m in off_list]
    c = E87.fit_c(off_list, inhs, crit['off'], mr.rate[reading]['on'])
    form = rec_factors([E87.apply_form(m, c, i, f_ch.kind, crit['off']) for m, i in zip(off_list, inhs)]).mean()
    return {'apl': mr.mean_factors(reading), 'cut_own': {'on': cut_own, 'off': off}, 'cut_flies': {'on': cut_fl, 'off': off},
            'form': {'on': form, 'off': off}}


def part_b1(out: Path, members: Dict[str, List[ModelReading]], ids: Dict[str, list]) -> pd.DataFrame:
    """B1: per group the mean over its pairs of log(L on / L off) and of log(R on / R off) (and of G, O, s) for the
    flies (z, readings 1 and 2; 90 % intervals over N_BOOT resamples; all and the six typical) and each population and
    wiring in every model reading. changes.csv (long: one row per source, reading, factor and group)."""
    rows = []
    for subset in SUBSETS:
        for frd in ('r1', 'r2'):
            rf = flies_factors('z', frd, flies_subset(subset))
            F = {c: rf[c].mean() for c in CONDS}
            ch = log_changes(F)
            bch = [log_changes(b) for b in flies_boot(rf)]
            for k, v in ch.items():
                _, f, g = k.split('_', 2)
                lo, hi = np.nanpercentile([b[k] for b in bch], INTERVAL)
                rows.append(dict(source='flies', population='flies', member='mean', subset=subset, reading=frd, factor=f, group=g,
                                 log_change=v, ratio=float(np.exp(v)), lo=float(lo), hi=float(hi)))
    for pop, lst in members.items():
        for rd in MODEL_READINGS:
            per = [log_changes(mr.mean_factors(rd)) for mr in lst]
            ch = log_changes(population(lst, rd))
            for k, v in ch.items():
                _, f, g = k.split('_', 2)
                ws = [p[k] for p in per]
                rows.append(dict(source='model', population=pop, member='population', subset='', reading=rd, factor=f, group=g,
                                 log_change=v, ratio=float(np.exp(v)), lo=float(np.nanmin(ws)), hi=float(np.nanmax(ws))))
                for (w, _, _), x in zip(ids[pop], ws):
                    rows.append(dict(source='model', population=pop, member=w, subset='', reading=rd, factor=f, group=g,
                                     log_change=x, ratio=float(np.exp(x)), lo=np.nan, hi=np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'changes.csv', index=False)
    return df


def b1_differs(ch: pd.DataFrame, pop: str, group: str, factor: str, subset: str = 'all') -> dict:
    """B1's rule for one population, group and factor: at each of readings 1 and 1b the population's change against
    the flies' (reading 1) interval; differs if outside it and by a factor of B1_FACTOR or more, at both."""
    fl = ch[(ch.source == 'flies') & (ch.subset == subset) & (ch.reading == 'r1') & (ch.factor == factor) & (ch.group == group)].iloc[0]
    out = dict(flies=float(fl.log_change), flies_lo=float(fl.lo), flies_hi=float(fl.hi))
    both = True
    for rd in B_READINGS:
        m = ch[(ch.source == 'model') & (ch.population == pop) & (ch.member == 'population') & (ch.reading == rd)
               & (ch.factor == factor) & (ch.group == group)].iloc[0]
        x = float(m.log_change)
        outside = not (fl.lo <= x <= fl.hi)
        big = abs(x - fl.log_change) >= np.log(B1_FACTOR)
        out[rd] = dict(model=x, wiring_min=float(m.lo), wiring_max=float(m.hi), outside=bool(outside), by_factor=bool(big),
                       differs=bool(outside and big))
        both &= outside and big
    out['differs'] = bool(both)
    return out


def _b2_side(d_apl: float, d_cut: float) -> str:
    """B2's rule on the differences (model minus flies) of the changes under APL and under the cut."""
    if not np.isfinite(d_apl) or d_apl == 0:
        return 'no difference under APL'
    sg = np.sign(d_apl)
    if d_cut * sg >= B2_LARGE * abs(d_apl):
        return 'in the responses without APL'
    if abs(d_cut) <= B2_SMALL * abs(d_apl) and (d_apl - d_cut) * sg >= B2_LARGE * abs(d_apl):
        return 'in what APL does to them'
    return 'both'


def part_b2(out: Path, members: Dict[str, List[ModelReading]], ids: Dict[str, list]) -> pd.DataFrame:
    """B2: per group the change of L (and R) under APL, under a per-odor cut of the responses without APL to the APL-on
    fractions (the model's own and the flies'), and under the flies' chosen form; flies (z, trial means; all and the six
    typical) and each population and wiring in readings 1 and 1b. cut.csv."""
    rows = []
    for subset in SUBSETS:
        fl = _flies_cut_form(flies_subset(subset))
        for what, F in fl.items():
            for k, v in log_changes(F, ('L', 'R')).items():
                _, f, g = k.split('_', 2)
                rows.append(dict(source='flies', population='flies', member='mean', subset=subset, reading='r1', what=what,
                                 factor=f, group=g, log_change=v, ratio=float(np.exp(v))))
    fc = flies_fractions()
    for pop, lst in members.items():
        for rd in B_READINGS:
            per = [model_cut_form(mr, rd, fc) for mr in lst]
            for what in ('apl', 'cut_own', 'cut_flies', 'form'):
                Fp = _pop_mean([p[what] for p in per])
                for k, v in log_changes(Fp, ('L', 'R')).items():
                    _, f, g = k.split('_', 2)
                    rows.append(dict(source='model', population=pop, member='population', subset='', reading=rd, what=what,
                                     factor=f, group=g, log_change=v, ratio=float(np.exp(v))))
                for (w, _, _), p in zip(ids[pop], per):
                    for k, v in log_changes(p[what], ('L', 'R')).items():
                        _, f, g = k.split('_', 2)
                        rows.append(dict(source='model', population=pop, member=w, subset='', reading=rd, what=what,
                                         factor=f, group=g, log_change=v, ratio=float(np.exp(v))))
        print(f'B2 {pop}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'cut.csv', index=False)
    return df


def _rank_corr(u: np.ndarray, v: np.ndarray) -> float:
    from scipy.stats import rankdata
    if len(u) < 2:
        return float('nan')
    a, b = rankdata(u), rankdata(v)
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def shared_amplitudes(x: np.ndarray, q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """For one thresholded matrix (KCs x odors) per pair (triu order): the rank correlation of the two odors' responses
    over their shared responders (NaN with fewer than B3_MIN_SHARED), and the share of the shared responders that lie in
    the top `q[o]` of each odor's responders (ties at the edge kept)."""
    i, j = _iu()
    r = x > 0
    top = np.zeros_like(r)
    for o in range(x.shape[1]):
        n = int(r[:, o].sum())
        k = int(round(float(q[o]) * n))
        if k > 0:
            th = np.sort(x[r[:, o], o])[-k]
            top[:, o] = r[:, o] & (x[:, o] >= th)
    rc, sh = np.full(len(i), np.nan), np.full(len(i), np.nan)
    for p, (a, b) in enumerate(zip(i, j)):
        s = r[:, a] & r[:, b]
        n = int(s.sum())
        if n:
            sh[p] = float((s & top[:, a] & top[:, b]).sum()) / n
        if n >= B3_MIN_SHARED:
            rc[p] = _rank_corr(x[s, a], x[s, b])
    return rc, sh


def part_b3(out: Path, members: Dict[str, List[ModelReading]]) -> pd.DataFrame:
    """B3: per group and condition the rank correlation of the two responses over a pair's shared responders and the
    share of the shared responders in the top fraction of both odors' responders (each odor's f_on / f_off, the share
    of its APL-off responders a per-odor cut keeps; the same fraction in both conditions), each the mean over the
    recordings (pseudo-recordings), then over the group's pairs: flies (z, trial means; all and the six typical) and each
    population in readings 1 and 1b (means over its wirings). shared_amplitudes.csv."""
    rows = []
    masks = group_masks()
    for subset in SUBSETS:
        rd = E87.Reading('z', 'trial_means', tnt=flies_subset(subset))
        mats = {'on': [E87.thresholded(m, rd.crit) for m in rd.tm_ctrl], 'off': [E87.thresholded(m, rd.crit) for m in rd.tm_tnt]}
        f = {c: np.mean([(m > 0).mean(axis=0) for m in mats[c]], axis=0) for c in CONDS}
        with np.errstate(divide='ignore', invalid='ignore'):
            q = np.clip(np.where(f['off'] > 0, f['on'] / f['off'], 0.0), 0.0, 1.0)
        for c in CONDS:
            res = [shared_amplitudes(m, q) for m in mats[c]]
            rc, sh = _mean_rows(np.stack([x[0] for x in res])), _mean_rows(np.stack([x[1] for x in res]))
            for g in GROUPS + ('all',):
                rows.append(dict(source='flies', population='flies', subset=subset, reading='r1', condition=c, group=g,
                                 rank_corr=_nanmean(rc[masks[g]]), top_share=_nanmean(sh[masks[g]])))
    for pop, lst in members.items():
        for rdn in B_READINGS:
            acc = {c: ([], []) for c in CONDS}
            for mr in lst:
                tms = mr.thresholded_means(rdn)
                f = {c: np.mean([(m > 0).mean(axis=0) for m in tms[c]], axis=0) for c in CONDS}
                with np.errstate(divide='ignore', invalid='ignore'):
                    q = np.clip(np.where(f['off'] > 0, f['on'] / f['off'], 0.0), 0.0, 1.0)
                for c in CONDS:
                    res = [shared_amplitudes(m, q) for m in tms[c]]
                    acc[c][0].append(_mean_rows(np.stack([x[0] for x in res])))
                    acc[c][1].append(_mean_rows(np.stack([x[1] for x in res])))
            for c in CONDS:
                rc, sh = _mean_rows(np.stack(acc[c][0])), _mean_rows(np.stack(acc[c][1]))
                for g in GROUPS + ('all',):
                    rows.append(dict(source='model', population=pop, subset='', reading=rdn, condition=c, group=g,
                                     rank_corr=_nanmean(rc[masks[g]]), top_share=_nanmean(sh[masks[g]])))
        print(f'B3 {pop}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'shared_amplitudes.csv', index=False)
    return df


def _kc_table(w: str) -> pd.DataFrame:
    """Per KC of wiring `w` (indexed by KC id): E70's `wiring_features` of the wiring's claw x bouton matrix, the base
    breadth behind the family's APL weights (E83's `wiring_parts`, E70's `base_breadth` of the bouton base), its breadth
    factor ((b + 1) / mean(b + 1)) ** 2, the resulting APL->KC weight (normalized to mean 1), and the claws per
    glomerulus (columns 'glom:<name>')."""
    key = ('kc_table', w)
    if key not in _CACHE:
        from experiments import E10_joint_fit as J
        from experiments import E70_high_firing_kcs as E70
        from experiments import E83_cross_wiring_fit as E
        p = E.wiring_parts(w)
        W = p['W']
        f = E70.wiring_features(W)
        b = p['breadth'].reindex(f.index).astype(float)
        bp1 = b + 1.0
        f['breadth'] = b.to_numpy()
        f['breadth_factor'] = ((bp1 / bp1.mean()) ** E.BREADTH_K).to_numpy()
        wk = J.breadth_weights(p['wnorm']['wAPLKC'], p['breadth'], E.BREADTH_K)
        wk = pd.Series(np.asarray(wk, dtype=float), index=wk.index.get_level_values(H.KC_ID_LEVEL))
        wk = wk.groupby(level=0).sum().reindex(f.index)
        f['apl_weight'] = (wk / wk.mean()).to_numpy()
        glom = W.columns.get_level_values('glomerulus')
        bg = W.T.groupby(glom).sum().T.groupby(level=H.KC_ID_LEVEL, sort=False).sum().reindex(f.index)
        bg.columns = [f'glom:{c}' for c in bg.columns]
        _CACHE[key] = pd.concat([f, bg], axis=1)
    return _CACHE[key]


def shared_kcs(o: Dict[str, tuple], kc: pd.DataFrame) -> Tuple[List[dict], Dict[str, set], pd.Series]:
    """B4 on one run's deterministic counts (a response two or more spikes): for each group's pairs the shared
    responders with APL off divided into 'survive' (shared responders with APL on too), 'cut_kept_apl_removed' (kept for
    both odors by a per-odor cut of the APL-off counts to the run's APL-on fractions, not surviving) and 'rest'; per
    group and class the KC-pair entries' features. Returns (rows, {group: KC ids of class 'cut_kept_apl_removed'},
    the responders' breadth factors)."""
    on_df, off_df = H.abbrev_columns(o['on'][1]), H.abbrev_columns(o['off'][1])
    ids = pd.Index(off_df.index.get_level_values(H.KC_ID_LEVEL))
    on = on_df[list(H.APL16_ODORS)].to_numpy(dtype=float)
    off = off_df[list(H.APL16_ODORS)].to_numpy(dtype=float)
    ron, roff = on > RAW_CRIT, off > RAW_CRIT
    f_on = ron.mean(axis=0)
    cut = E87.per_odor_oracle(E87.thresholded(off, RAW_CRIT), f_on, 'divisive', RAW_CRIT) > 0
    n_odors_off = roff.sum(axis=1)
    # each odor's APL-off responders ranked by amplitude (percentile among them; 1 the largest)
    rank = np.full(off.shape, np.nan)
    for j in range(off.shape[1]):
        r = roff[:, j]
        if r.any():
            from scipy.stats import rankdata
            rank[r, j] = rankdata(off[r, j]) / r.sum()
    k = kc.reindex(ids)
    responders = n_odors_off > 0
    bf = pd.Series(k['breadth_factor'].to_numpy(), index=ids)[responders]
    glcols = [c for c in k.columns if c.startswith('glom:')]
    rows, kept_removed = [], {}
    i, j = _iu()
    masks = group_masks()
    for g in GROUPS + ('all',):
        ent = {cl: [] for cl in ('survive', 'cut_kept_apl_removed', 'rest')}
        n_shared, dot_total = 0, 0.0
        for p in np.flatnonzero(masks[g]):
            a, b = i[p], j[p]
            sh = roff[:, a] & roff[:, b]
            n_shared += int(sh.sum())
            dot = off[:, a] * off[:, b]
            dot_total += float(dot[sh].sum())
            surv = sh & ron[:, a] & ron[:, b]
            ckr = sh & cut[:, a] & cut[:, b] & ~surv
            rest = sh & ~surv & ~ckr
            for cl, msk in (('survive', surv), ('cut_kept_apl_removed', ckr), ('rest', rest)):
                for kk in np.flatnonzero(msk):
                    ent[cl].append((kk, float(dot[kk]), 0.5 * (rank[kk, a] + rank[kk, b])))
        kept_removed[g] = {ids[e[0]] for e in ent['cut_kept_apl_removed']}
        for cl, es in ent.items():
            row = dict(group=g, kc_class=cl, n_entries=len(es), n_kcs=len({e[0] for e in es}),
                       share_of_overlap=len(es) / n_shared if n_shared else np.nan,
                       share_of_dot=sum(e[1] for e in es) / dot_total if dot_total else np.nan)
            if es:
                rr = np.array([e[0] for e in es])
                kk_ = k.iloc[rr]
                row.update(n_odors_off=float(n_odors_off[rr].mean()), breadth=float(kk_.breadth.mean()),
                           breadth_factor=float(kk_.breadth_factor.mean()), apl_weight=float(kk_.apl_weight.mean()),
                           amp_rank=float(np.nanmean([e[2] for e in es])), n_claws=float(kk_.n_claws.mean()),
                           top_glomerulus_share=float((kk_.max_claws_one_glomerulus / kk_.n_claws).mean()),
                           top_quarter_share=float((kk_.breadth_factor.to_numpy() >= np.quantile(bf, 1 - B4_TOP)).mean()))
                for t, v in kk_.kc_type.value_counts(normalize=True).items():
                    row[f'type_{t}'] = float(v)
                gs = kk_[glcols].sum()
                gs = gs[gs > 0].sort_values(ascending=False)
                row['top_glomeruli'] = ', '.join(f'{c[5:]} {v / gs.sum():.2f}' for c, v in gs.head(4).items()) if len(gs) else ''
            rows.append(row)
    return rows, kept_removed, bf


def part_b4(out: Path, members: Dict[str, List[ModelReading]], ids: Dict[str, list]) -> Tuple[pd.DataFrame, dict]:
    """B4: the model's KCs that carry each group's overlap with APL off, on the breadth-2 fits' eight wirings and D2's
    arm with per-odor APL scales (deterministic counts): per wiring and pooled. shared_kcs.csv; the breadth-weighting
    rule per population and group (returned)."""
    rows, rule = [], {}
    for pop in ('b2', f'D2_{D2_BASE_ARM}'):
        pooled_kr: Dict[str, List[bool]] = {g: [] for g in GROUPS + ('all',)}
        prow = []
        for mr, (w, folder, tag) in zip(members[pop], ids[pop]):
            kc = _kc_table(w)
            rs, kr, bf = shared_kcs(mr.o, kc)
            thr = float(np.quantile(bf, 1 - B4_TOP))
            for r in rs:
                rows.append(dict(population=pop, wiring=w, **r))
                prow.append(r)
            for g, s in kr.items():
                pooled_kr[g] += [bool(kc.loc[x, 'breadth_factor'] >= thr) for x in s]
            print(f'B4 {pop} {w}', flush=True)
        df_p = pd.DataFrame(prow)
        for (g, cl), sub in df_p.groupby(['group', 'kc_class']):
            tot = sub.n_entries.sum()
            row = dict(population=pop, wiring='pooled', group=g, kc_class=cl, n_entries=int(tot), n_kcs=int(sub.n_kcs.sum()))
            for col in ('share_of_overlap', 'share_of_dot'):
                row[col] = float(sub[col].mean())
            for col in ('n_odors_off', 'breadth', 'breadth_factor', 'apl_weight', 'amp_rank', 'n_claws', 'top_glomerulus_share',
                        'top_quarter_share'):
                if col in sub and tot:
                    row[col] = float(np.nansum(sub[col] * sub.n_entries) / tot)
            rows.append(row)
        rule[pop] = {g: dict(n_kcs=len(v), share_top_quarter=float(np.mean(v)) if v else float('nan'),
                             breadth_weighting_named=bool(v and np.mean(v) >= B4_SHARE)) for g, v in pooled_kr.items()}
    df = pd.DataFrame(rows)
    df.to_csv(out / 'shared_kcs.csv', index=False)
    return df, rule


def verdict_b(out: Path, b4_rule: dict) -> dict:
    """B's verdict: B1 (per population and group, L and R, all recordings and the six typical), B2's side for the
    groups whose co-response change differs (at both noise levels), B3's summary and B4's rule. verdict_B.json."""
    ch = pd.read_csv(out / 'changes.csv')
    cut = pd.read_csv(out / 'cut.csv')
    sa = pd.read_csv(out / 'shared_amplitudes.csv')
    v = dict(B1={}, B2={}, B3={}, B4=b4_rule)
    for subset in SUBSETS:
        b1, b2 = {}, {}
        for pop in POPULATIONS:
            b1[pop] = {g: {f: b1_differs(ch, pop, g, f, subset) for f in ('L', 'R')} for g in GROUPS}
            gsides = {}
            for g in GROUPS:
                fl = cut[(cut.source == 'flies') & (cut.subset == subset) & (cut.factor == 'L') & (cut.group == g)].set_index('what').log_change
                sides = {}
                for rd in B_READINGS:
                    m = cut[(cut.source == 'model') & (cut.population == pop) & (cut.member == 'population') & (cut.reading == rd)
                            & (cut.factor == 'L') & (cut.group == g)].set_index('what').log_change
                    d_apl = float(m['apl'] - fl['apl'])
                    d_cut = float(m['cut_own'] - fl['cut'])
                    d_form = float(m['form'] - fl['form'])
                    sides[rd] = dict(model_apl=float(m['apl']), flies_apl=float(fl['apl']), model_cut=float(m['cut_own']),
                                     flies_cut=float(fl['cut']), model_cut_flies_fractions=float(m['cut_flies']),
                                     model_form=float(m['form']), flies_form=float(fl['form']), d_apl=d_apl, d_cut=d_cut,
                                     d_departure=d_apl - d_cut, d_form=d_form, side=_b2_side(d_apl, d_cut))
                s = [sides[rd]['side'] for rd in B_READINGS]
                sides['side'] = s[0] if all(x == s[0] for x in s) else 'undecided (the noise levels disagree)'
                sides['co_response_change_differs'] = b1[pop][g]['L']['differs']
                gsides[g] = sides
            b2[pop] = {g: x for g, x in gsides.items() if x['co_response_change_differs']}
            b2[pop]['all_groups'] = {g: x['side'] for g, x in gsides.items()}
        v['B1'][subset] = b1
        v['B2'][subset] = b2
    for pop in ('flies',) + POPULATIONS:
        g = sa[(sa.population == pop) & ((sa.subset == 'all') | (sa.source == 'model'))]
        v['B3'][pop] = {f'{r.reading}_{r.condition}_{r.group}': dict(rank_corr=r.rank_corr, top_share=r.top_share)
                        for r in g.itertuples() if r.group in ('ester', 'alcohol', 'nearest4')}
    changed = []
    for pop in POPULATIONS:
        for g in GROUPS:
            for f in ('L', 'R'):
                if v['B1']['all'][pop][g][f]['differs'] != v['B1']['six'][pop][g][f]['differs']:
                    changed.append(f'B1 {pop} {g} {f}')
            sa_, sb_ = v['B2']['all'][pop]['all_groups'][g], v['B2']['six'][pop]['all_groups'][g]
            if sa_ != sb_:
                changed.append(f'B2 {pop} {g}: {sa_} with all, {sb_} with the six typical')
    v['depends_on_the_two_strongest_recordings'] = changed
    v['thresholds'] = dict(B1_FACTOR=B1_FACTOR, B2_LARGE=B2_LARGE, B2_SMALL=B2_SMALL, B4_TOP=B4_TOP, B4_SHARE=B4_SHARE)
    _write_json(out / 'verdict_B.json', v)
    return v


def part_b(out: Path) -> dict:
    """B1 to B4 and B's verdict: changes.csv, cut.csv, shared_amplitudes.csv, shared_kcs.csv, verdict_B.json."""
    members, ids, rel = _setup_members(out)
    _provenance(out, 'B', flies_reliability=rel, populations={k: [list(x) for x in v] for k, v in ids.items()},
                chosen_form=_chosen_form().name)
    part_b1(out, members, ids)
    print('B1 done', flush=True)
    part_b2(out, members, ids)
    part_b3(out, members)
    _, rule = part_b4(out, members, ids)
    v = verdict_b(out, rule)
    print('B verdict written', flush=True)
    return v


# ---------------------------------------------------------------------------
# part C

C_READINGS = ('r1', 'r3raw', 'r3cal')
_W: dict = {}


def _c_init() -> None:
    warnings.simplefilter('ignore')
    _W['flies'] = {}
    for frd in ('r1', 'r2'):
        rf = flies_factors('z', frd)
        _W['flies'][frd] = {c: rf[c].mean() for c in CONDS}


def _c_job(job: tuple) -> List[dict]:
    """One variant (folder, variant, source, on path, off path) under C_READINGS: the cells of A2 (pattern part, mean,
    ρ, ester / alcohol errors), the binarized pattern part and the changes of B1, against the flies (all recordings)."""
    folder, variant, source, on, off = job
    base = dict(folder=folder, variant=variant, source=source)
    try:
        if 'flies' not in _W:
            _c_init()
        o = {'on': (None, pd.read_parquet(on)), 'off': (None, pd.read_parquet(off))}
        if not set(H.APL16_ODORS) <= set(H.abbrev_columns(o['on'][1]).columns):
            return [dict(base, error='not the APL16 odors')]
        mr = ModelReading(o, label=variant)
        Fs = {rd: mr.mean_factors(rd) for rd in C_READINGS}
        rows = []
        for rd in C_READINGS:
            Fm = Fs[rd]
            Ff = _W['flies'][FLIES_READING[rd]]
            data = effect(Ff)
            row = dict(base, reading=rd, n_kc=int(mr.counts['on'].shape[0]), rate_on=mr.rate[rd]['on'], rate_off=mr.rate[rd]['off'],
                       z_K=mr.cal['K'], z_noise=mr.cal['noise'], K_readout=mr.K)
            for name, take in CELLS.items():
                if name == 'flies':
                    continue
                s = scores88(cell_effect(Fm, Ff, take), data)
                row.update({f'{name}_{k}': s[k] for k in ('pattern', 'loss', 'mean', 'rho', 'err_ester', 'err_alcohol')})
            sb = scores88(effect(Fm, 'O'), effect(Ff, 'O'))
            row.update(bin_pattern=sb['pattern'], bin_mean=sb['mean'], bin_rho=sb['rho'])
            row.update({k: v for k, v in log_changes(Fm, ('L', 'R')).items() if k.split('_', 2)[2] in ('ester', 'alcohol', 'nearest4')})
            rows.append(row)
        return rows
    except Exception as e:   # a malformed file is listed, not fatal
        return [dict(base, error=f'{type(e).__name__}: {e}'[:300])]


def _wiring_of(tag: str) -> str:
    """The wiring a tag gives (E72's draw prefixes anywhere in it, or an E85 model fly's last field), else 'real'."""
    from experiments import E72_wiring_population as E72
    from experiments import E83_cross_wiring_fit as E
    for lab, prefix in E72.DRAWS:
        if prefix in tag:
            return lab
    last = tag.rsplit('_', 1)[-1]
    if last in E.ALL_WIRINGS:
        return last
    return 'real'


def _handle_in_code(folder: str) -> bool:
    m = re.match(r'E(\d+)_', folder)
    here = Path(__file__).resolve().parent
    return bool(m and any(here.glob(f'E{m.group(1)}_*.py')))


def _realism_own(folder: str, variant: str, on: str, off: str, bv: pd.DataFrame) -> str:
    """'inside the bands', the flag, or 'not scored': results/apl_effect/by_variant.csv where it has the response
    rates, else E77's `realism_of` on the run's spike counts (as E87's `rescore_archive`)."""
    from experiments import E77_input_side_common_error as E77
    row = bv[(bv.folder == folder) & (bv.tag == variant)]
    if len(row) and pd.notna(row.rate_on.iloc[0]):
        f = row.realism_flag.iloc[0]
        return f if isinstance(f, str) and f else 'inside the bands'
    try:
        s_on, s_off = pd.read_parquet(on), pd.read_parquet(off)
        o = {'on': ((s_on >= 2).astype(float), s_on), 'off': ((s_off >= 2).astype(float), s_off)}
        r = E77.realism_of(o)
        return r['realism_flag'] if isinstance(r['realism_flag'], str) and r['realism_flag'] else 'inside the bands'
    except Exception as e:
        return f'not scored ({type(e).__name__})'


def c_jobs() -> List[tuple]:
    """Every archived variant (E87's `archive_variants`, E86's excepted) and E85's 72 model flies (viki's
    outputs/E85_orn_fly_population/B/, where the archive does not already hold them)."""
    jobs = [(f, v, 'archive', str(on), str(off)) for f, v, on, off in E87.archive_variants()]
    have = {(j[0], j[1]) for j in jobs}
    d = H.outputs_dir('E85_orn_fly_population') / 'B'
    for p in sorted(d.glob('B1[ST]_repro_f*_*_on_spike_counts.parquet')):
        tag = p.name[:-len('_on_spike_counts.parquet')]
        q = d / f'{tag}_acute_off_spike_counts.parquet'
        if q.exists() and ('E85_orn_fly_population', tag) not in have:   # archived ones are scanned once
            jobs.append(('E85_orn_fly_population', tag, 'E85_model_flies', str(p), str(q)))
    return jobs


def part_c(out: Path, n_workers: int = 8) -> dict:
    """The archive by the factors (no simulations): each variant under reading 1 (its own calibration) and reading 3
    (raw counts and the calibrated readout) against the flies (z, all recordings; reading 3 against the flies with the
    trial noise removed); joined with by_variant.csv's realism flag; the stand-out rule under readings 1 and 3 against
    the variant's baseline (the breadth-2 fit on the same wiring, else the real wiring's), with realism on its own run
    and the handle in the present code. archive_scan.csv, verdict_C.json."""
    from concurrent.futures import ProcessPoolExecutor
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_THREAD_LIMIT'):
        os.environ[k] = '1'
    jobs = c_jobs()
    _provenance(out, 'C', n_variants=len(jobs), n_workers=n_workers)
    print(f'C: {len(jobs)} variants on {n_workers} workers', flush=True)
    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_c_init) as ex:
        for n, r in enumerate(ex.map(_c_job, jobs, chunksize=4)):
            rows += r
            if n % 100 == 0:
                print(f'C {n} of {len(jobs)} ({time.time() - t0:.0f}s)', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'archive_scan_raw.csv', index=False)
    return score_c(out, df)


def score_c(out: Path, df: Optional[pd.DataFrame] = None) -> dict:
    """C's stand-out rule and verdict on the scan's rows (archive_scan_raw.csv). archive_scan.csv, verdict_C.json."""
    from experiments import E83_cross_wiring_fit as E
    if df is None:
        df = pd.read_csv(out / 'archive_scan_raw.csv')
    # a run scanned from two places (the archive and an experiment's outputs) is kept once
    df = df.drop_duplicates(['folder', 'variant', 'reading'], keep='first').reset_index(drop=True)
    fams = E87.model_families()
    fam_of = {(f, t): fam for fam, mem in fams.items() for _, f, t in mem}
    df['family'] = [fam_of.get((f, v), 'E85 model flies' if s == 'E85_model_flies' else '')
                    for f, v, s in zip(df.folder, df.variant, df.source)]
    df['wiring'] = df.variant.astype(str).map(_wiring_of)
    df['handle_in_code'] = df.folder.astype(str).map(_handle_in_code)
    bv = pd.read_csv(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv')
    flag = bv.drop_duplicates(['folder', 'tag']).set_index(['folder', 'tag']).realism_flag
    df['realism_flag_by_variant'] = [flag.get((f, v), np.nan) for f, v in zip(df.folder, df.variant)]
    # the flies' changes per reading (z, all recordings)
    fl = {}
    for frd in ('r1', 'r2'):
        rf = flies_factors('z', frd)
        fl[frd] = log_changes({c: rf[c].mean() for c in CONDS}, ('L', 'R'))
    base = {}
    for w, (f, t) in E.members().items():
        for rd in C_READINGS:
            m = df[(df.folder == f) & (df.variant == t) & (df.reading == rd)]
            if len(m):
                base[(w, rd)] = m.iloc[0]
    cols = ('bin_pattern', 'chg_L_alcohol', 'chg_L_ester', 'chg_R_ester')
    for c in cols:
        df[f'base_{c}'] = [base[(w, rd)][c] if (w, rd) in base else np.nan for w, rd in zip(df.wiring, df.reading)]
    df['flies_chg_L_alcohol'] = [fl[FLIES_READING[rd]]['chg_L_alcohol'] if isinstance(rd, str) else np.nan for rd in df.reading]
    df['flies_chg_L_ester'] = [fl[FLIES_READING[rd]]['chg_L_ester'] if isinstance(rd, str) else np.nan for rd in df.reading]
    df['flies_chg_R_ester'] = [fl[FLIES_READING[rd]]['chg_R_ester'] if isinstance(rd, str) else np.nan for rd in df.reading]
    df['bin_gain'] = df.base_bin_pattern - df.bin_pattern
    with np.errstate(divide='ignore', invalid='ignore'):
        df['alcohol_L_progress'] = (df.chg_L_alcohol - df.base_chg_L_alcohol) / (df.flies_chg_L_alcohol - df.base_chg_L_alcohol)
    df['ester_L_no_further'] = (df.chg_L_ester - df.flies_chg_L_ester).abs() <= (df.base_chg_L_ester - df.flies_chg_L_ester).abs()
    for k, f in (('ester_R', 'chg_R_ester'), ('ester_L', 'chg_L_ester'), ('alcohol_L', 'chg_L_alcohol')):
        df[f'{k}_towards_flies'] = (df[f] - df[f'flies_{f}']).abs() < (df[f'base_{f}'] - df[f'flies_{f}']).abs()
    # "at least halfway from the baseline's to the flies'": the variant's change ends within (1 - C_HALFWAY) of the
    # baseline's distance from the flies' (progress 0.5 to 1.5; an overshoot past that is not counted)
    df['alcohol_L_halfway'] = ((df.chg_L_alcohol - df.flies_chg_L_alcohol).abs()
                               <= (1 - C_HALFWAY) * (df.base_chg_L_alcohol - df.flies_chg_L_alcohol).abs())
    df['cand'] = df.reading.isin(['r1', 'r3raw']) & ((df.bin_gain >= C_BIN_GAIN) | (df.alcohol_L_halfway & df.ester_L_no_further))
    df['cand_by'] = np.where(df.cand & (df.bin_gain >= C_BIN_GAIN), 'binarized pattern part',
                             np.where(df.cand, 'alcohol co-response change', ''))
    roots = H.shared_dir() / 'data' / 'model_outputs'
    jobs = {(f, v): (on, off) for f, v, s, on, off in c_jobs()}
    real = {}
    for f, v in df.loc[df.cand, ['folder', 'variant']].drop_duplicates().itertuples(index=False):
        on, off = jobs.get((f, v), (None, None))
        real[(f, v)] = _realism_own(f, v, on, off, bv) if on else 'not scored'
    df['realism_own'] = [real.get((f, v), np.nan) for f, v in zip(df.folder, df.variant)]
    df['standout'] = df.cand & (df.realism_own == 'inside the bands') & df.handle_in_code
    df.to_csv(out / 'archive_scan.csv', index=False)
    ok = df[df.error.isna()] if 'error' in df else df
    by_handle = []
    for (rd, folder), g in ok[ok.reading.isin(['r1', 'r3raw'])].groupby(['reading', 'folder']):
        by_handle.append(dict(reading=rd, folder=folder, n=len(g), best_bin_gain=float(g.bin_gain.max()),
                              best_variant=str(g.loc[g.bin_gain.idxmax(), 'variant']) if g.bin_gain.notna().any() else '',
                              median_bin_pattern=float(g.bin_pattern.median()),
                              share_ester_R_towards=float(g.ester_R_towards_flies.mean()),
                              share_ester_L_towards=float(g.ester_L_towards_flies.mean()),
                              share_alcohol_L_towards=float(g.alcohol_L_towards_flies.mean()),
                              best_alcohol_L_progress=float(g.alcohol_L_progress.max()), n_cand=int(g.cand.sum()),
                              n_standout=int(g.standout.sum())))
    bh = pd.DataFrame(by_handle).sort_values(['reading', 'best_bin_gain'], ascending=[True, False])
    bh.to_csv(out / 'archive_by_handle.csv', index=False)
    fam_rank = {}
    for rd in ('r1', 'r3raw'):
        g = ok[(ok.reading == rd) & (ok.family != '')]
        fam_rank[rd] = {f'{r.family}:{r.wiring}': dict(bin_pattern=r.bin_pattern, own_pattern=r.own_pattern, sets_pattern=r.sets_pattern,
                                                        chg_L_alcohol=r.chg_L_alcohol, chg_L_ester=r.chg_L_ester, chg_R_ester=r.chg_R_ester)
                        for r in g.itertuples()}
    st = df[df.standout]
    v = dict(n_variants=int(df[['folder', 'variant']].drop_duplicates().shape[0]), n_errors=int(df.error.notna().sum()) if 'error' in df else 0,
             n_candidates=int(df.cand.sum()), n_standout=int(len(st)),
             n_standout_by={f'{k[0]}: {k[1]}': int(v) for k, v in st.groupby(['reading', 'cand_by']).size().items()} if len(st) else {},
             standout_variants=int(st[['folder', 'variant']].drop_duplicates().shape[0]),
             standout_handles={f: int(n) for f, n in st.drop_duplicates(['folder', 'variant']).folder.value_counts().items()},
             standouts=st[['folder', 'variant', 'reading', 'wiring', 'cand_by', 'bin_pattern', 'base_bin_pattern', 'alcohol_L_progress',
                           'chg_L_alcohol', 'base_chg_L_alcohol', 'realism_own']].values.tolist(),
             candidates_outside_bands=df[df.cand & (df.realism_own != 'inside the bands')][['folder', 'variant', 'reading', 'realism_own']].values.tolist()[:50],
             flies_changes={rd: fl[rd] for rd in fl}, families=fam_rank,
             handles_top=bh.head(15).to_dict(orient='records'),
             thresholds=dict(C_BIN_GAIN=C_BIN_GAIN, C_HALFWAY=C_HALFWAY),
             handle_rule='a handle the present code has: the variant\'s experiment has a script in experiments/ (whether it composes '
                         'with D2 is for the session to judge)')
    _write_json(out / 'verdict_C.json', v)
    print(f'C: {len(df)} rows, {int(df.cand.sum())} candidates, {len(st)} standouts', flush=True)
    return v


# ---------------------------------------------------------------------------
# part D (simulations; not implemented here)

def part_d1(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """D1: the flies' chosen form as a static inhibition in place of APL on S's four wirings with the flies' per-odor
    recruitment (E87's D2 thresholds, re-solved with `B3Pool`). Not implemented: part D is a later step."""
    raise NotImplementedError('E88 part D1 is not implemented yet')


def part_d2(out: Path, n_workers: int = 16) -> pd.DataFrame:
    """D2: up to three handles that stand out in C, on S's four wirings. Not implemented: part D is a later step."""
    raise NotImplementedError('E88 part D2 is not implemented yet')


# ---------------------------------------------------------------------------
# the comparison table

def _f(x, spec: str = '.4f') -> str:
    return E87._f(x, spec)


def comparison_table(out: Path) -> str:
    """comparison_table.md covering A to C, from the parts' tables."""
    lines = ['# E88: the common error split into who responds and how strongly', '',
             f'Generated by experiments/E88_responder_sets.py at {H.code_commits().get("queue_commit", "")[:7]}. z, sub-threshold '
             'responses set to 0; the pattern part (lower is better) against the flies\' per-pair effect in the same reading. '
             'Readings: r1 trial means (the model under E73\'s per-KC reading through its measurement noise), r1b the same with the '
             'noise raised to the flies\' repeatability, r2 the trial noise removed by the same estimator on both sides, r3raw / '
             'r3cal the deterministic model (raw counts / calibrated readout) against the flies with the noise removed.', '']
    if (out / 'members.csv').exists():
        mi = pd.read_csv(out / 'members.csv')
        lines += ['## Readings of the populations (mean over members)', '',
                  '| population | r1 rate on / off | r2 criterion | r2 reliability on / off | 1b multiplier on / off | r1b reliability on / off |',
                  '|---|---:|---:|---:|---:|---:|']
        for pop, g in mi.groupby('population', sort=False):
            lines.append(f'| {pop} | {_f(g.r1_rate_on.mean(), ".3f")} / {_f(g.r1_rate_off.mean(), ".3f")} | {_f(g.r2_crit_on.mean(), ".3f")} | '
                         f'{_f(g.r2_rel_on.mean(), ".2f")} / {_f(g.r2_rel_off.mean(), ".2f")} | {_f(g.mult_on_mult.mean(), ".2f")} / '
                         f'{_f(g.mult_off_mult.mean(), ".2f")} | {_f(g.mult_on_rel_r1b.mean(), ".2f")} / {_f(g.mult_off_rel_r1b.mean(), ".2f")} |')
    if (out / 'estimator_check.csv').exists():
        ec = pd.read_csv(out / 'estimator_check.csv')
        b = ec[(ec.population == 'b2') & (ec.group != 'all')]
        lines += ['', '## A0: the estimators on known truth (breadth-2 population; the largest |group mean error| over the five groups '
                  'and both conditions; s and O plain, R, L and G in units of s; not read beyond 0.03)', '',
                  '| reading | R | L | G | O | s |', '|---|---:|---:|---:|---:|---:|']
        for rd, g in b.groupby('reading'):
            lines.append(f'| {rd} | ' + ' | '.join(_f(g[g.factor == f].err.abs().max()) + (' (not read)' if g[g.factor == f].off.any() else '')
                                                  for f in FACTORS) + ' |')
    if (out / 'cells.csv').exists():
        cells = pd.read_csv(out / 'cells.csv')
        lines += ['', '## A2 / A3: the cells (pattern part; floor in brackets; ester / alcohol centred error)', '']
        for subset in SUBSETS:
            lines += [f'### {"all recordings" if subset == "all" else "the six typical recordings with APL output blocked"}', '',
                      '| population | reading | own | flies\' sets | flies\' amplitudes | flies\' R | flies\' L | binarized | ρ own / sets | '
                      'own ester / alcohol | sets ester / alcohol | R ester / alcohol | L ester / alcohol |',
                      '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            g = cells[cells.subset == subset].set_index(['population', 'reading', 'cell'])
            for pop in POPULATIONS:
                for rd in CELL_READINGS:
                    if (pop, rd, 'own') not in g.index:
                        continue
                    r = {c: g.loc[(pop, rd, c)] for c in ('own', 'sets', 'amp', 'R', 'L', 'binary')}
                    ea = lambda x: f"{_f(x.err_ester, '+.3f')} / {_f(x.err_alcohol, '+.3f')}"
                    lines.append(f"| {pop} | {rd} | {_f(r['own'].pattern)} | {_f(r['sets'].pattern)} ({_f(r['sets'].floor_pattern, '.3f')}) | "
                                 f"{_f(r['amp'].pattern)} ({_f(r['amp'].floor_pattern, '.3f')}) | {_f(r['R'].pattern)} | {_f(r['L'].pattern)} | "
                                 f"{_f(r['binary'].pattern)} | {_f(r['own'].rho, '.2f')} / {_f(r['sets'].rho, '.2f')} | {ea(r['own'])} | "
                                 f"{ea(r['sets'])} | {ea(r['R'])} | {ea(r['L'])} |")
            lines.append('')
            fr = cells[(cells.subset == subset) & (cells.cell == 'flies') & (cells.population == 'b2')]
            for _, r in fr.iterrows():
                lines.append(f'Flies ({r.reading}): mean effect {_f(r["mean"], "+.3f")}, the floor of the pattern part {_f(r.floor_pattern)}.')
            lines.append('')
    if (out / 'main_target.csv').exists():
        mt = pd.read_csv(out / 'main_target.csv')
        lines += ['## The main target beside (E83\'s effect_scores; all recordings)', '',
                  '| population | readout | pattern part | loss | mean APL effect | ρ | ester / alcohol |', '|---|---|---:|---:|---:|---:|---:|']
        for _, r in mt.iterrows():
            lines.append(f'| {r.population} | {r.readout} | {_f(r.z8_pattern)} | {_f(r.z8_loss)} | {_f(r.mean_effect, "+.3f")} | '
                         f'{_f(r.z8_rho, ".2f")} | {_f(r.z8_err_ester, "+.3f")} / {_f(r.z8_err_alcohol, "+.3f")} |')
        lines.append('')
    if (out / 'changes.csv').exists():
        ch = pd.read_csv(out / 'changes.csv')
        lines += ['## B1: what APL does to each factor (on over off, geometric mean over the group\'s pairs; four nearest / ester / alcohol)', '',
                  '| source | reading | recruitment R | co-response L | amplitude G |', '|---|---|---:|---:|---:|']

        def trio(sub, f):
            vals = []
            for gname in ('nearest4', 'ester', 'alcohol'):
                r = sub[(sub.factor == f) & (sub.group == gname)]
                if not len(r):
                    vals.append('–')
                    continue
                r = r.iloc[0]
                iv = f' [{_f(np.exp(r.lo), ".2f")}, {_f(np.exp(r.hi), ".2f")}]' if np.isfinite(r.lo) else ''
                vals.append(f'{_f(r.ratio, ".2f")}{iv}')
            return ' / '.join(vals)
        for subset in SUBSETS:
            for frd in ('r1', 'r2'):
                sub = ch[(ch.source == 'flies') & (ch.subset == subset) & (ch.reading == frd)]
                lines.append(f'| KC data ({subset}) | {frd} | {trio(sub, "R")} | {trio(sub, "L")} | {trio(sub, "G")} |')
        for pop in POPULATIONS:
            for rd in ('r1', 'r1b', 'r2', 'r3raw'):
                sub = ch[(ch.source == 'model') & (ch.population == pop) & (ch.member == 'population') & (ch.reading == rd)]
                lines.append(f'| {pop} | {rd} | {trio(sub, "R")} | {trio(sub, "L")} | {trio(sub, "G")} |')
        lines += ['', 'Model intervals are the range over wirings; the flies\' are 90 % over resampled recordings.', '']
    if (out / 'cut.csv').exists():
        cut = pd.read_csv(out / 'cut.csv')
        lines += ['## B2: the co-response under APL, under a per-odor cut and under the flies\' chosen form (L on over off; four nearest / ester / alcohol)', '',
                  '| source | reading | APL | cut (own fractions) | cut (flies\' fractions) | flies\' chosen form |', '|---|---|---:|---:|---:|---:|']

        def trio2(sub, what):
            return ' / '.join(_f(sub[(sub.what == what) & (sub.factor == 'L') & (sub.group == gname)].ratio.iloc[0], '.2f')
                              if len(sub[(sub.what == what) & (sub.factor == 'L') & (sub.group == gname)]) else '–'
                              for gname in ('nearest4', 'ester', 'alcohol'))
        for subset in SUBSETS:
            sub = cut[(cut.source == 'flies') & (cut.subset == subset)]
            lines.append(f'| KC data ({subset}) | r1 | {trio2(sub, "apl")} | {trio2(sub, "cut")} | – | {trio2(sub, "form")} |')
        for pop in POPULATIONS:
            for rd in B_READINGS:
                sub = cut[(cut.source == 'model') & (cut.population == pop) & (cut.member == 'population') & (cut.reading == rd)]
                lines.append(f'| {pop} | {rd} | {trio2(sub, "apl")} | {trio2(sub, "cut_own")} | {trio2(sub, "cut_flies")} | {trio2(sub, "form")} |')
        lines.append('')
    if (out / 'shared_amplitudes.csv').exists():
        sa = pd.read_csv(out / 'shared_amplitudes.csv')
        lines += ['## B3: the shared responders\' amplitudes (rank correlation / top share; ester / alcohol)', '',
                  '| source | reading | APL on | APL off |', '|---|---|---:|---:|']
        for (src, pop, subset, rd), g in sa.groupby(['source', 'population', 'subset', 'reading'], sort=False, dropna=False):
            def pr(c):
                return ' ; '.join(f"{_f(g[(g.condition == c) & (g.group == gg)].rank_corr.iloc[0], '.2f')} / "
                                  f"{_f(g[(g.condition == c) & (g.group == gg)].top_share.iloc[0], '.2f')}" for gg in ('ester', 'alcohol'))
            lines.append(f'| {pop}{" (" + subset + ")" if isinstance(subset, str) and subset else ""} | {rd} | {pr("on")} | {pr("off")} |')
        lines.append('')
    if (out / 'shared_kcs.csv').exists():
        sk = pd.read_csv(out / 'shared_kcs.csv')
        p = sk[sk.wiring == 'pooled']
        lines += ['## B4: the model\'s shared responders with APL off (pooled over wirings)', '',
                  '| population | group | class | entries | share of overlap | odors (APL off) | breadth factor | top quarter | amplitude rank | claws | top glomerulus share |',
                  '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for _, r in p[p.group.isin(['nearest4', 'ester', 'alcohol'])].iterrows():
            lines.append(f'| {r.population} | {r.group} | {r.kc_class} | {r.n_entries} | {_f(r.share_of_overlap, ".2f")} | {_f(r.get("n_odors_off"), ".1f")} | '
                         f'{_f(r.get("breadth_factor"), ".2f")} | {_f(r.get("top_quarter_share"), ".2f")} | {_f(r.get("amp_rank"), ".2f")} | '
                         f'{_f(r.get("n_claws"), ".1f")} | {_f(r.get("top_glomerulus_share"), ".2f")} |')
        lines.append('')
    if (out / 'archive_by_handle.csv').exists():
        bh = pd.read_csv(out / 'archive_by_handle.csv')
        lines += ['## C: the archive by handle (best binarized gain over the baseline; share of variants moving each change towards the flies)', '',
                  '| reading | folder | n | best gain | best variant | ester R | ester L | alcohol L | best alcohol progress | candidates | standouts |',
                  '|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|']
        for _, r in bh.groupby('reading').head(12).iterrows():
            lines.append(f'| {r.reading} | {r.folder} | {r.n} | {_f(r.best_bin_gain)} | {r.best_variant} | {_f(r.share_ester_R_towards, ".2f")} | '
                         f'{_f(r.share_ester_L_towards, ".2f")} | {_f(r.share_alcohol_L_towards, ".2f")} | {_f(r.best_alcohol_L_progress, ".2f")} | '
                         f'{r.n_cand} | {r.n_standout} |')
        lines.append('')
    lines += ['## Verdicts', '']
    for part in ('A', 'B', 'C'):
        p = out / f'verdict_{part}.json'
        if p.exists():
            v = json.loads(p.read_text())
            if part == 'B':
                v = {k: v[k] for k in ('B2', 'B4', 'depends_on_the_two_strongest_recordings') if k in v}
            if part == 'C':
                v = {k: v[k] for k in v if k not in ('families', 'handles_top')}
            lines += [f'### {part}', '', '```', json.dumps(v, indent=1)[:20000], '```', '']
    text = '\n'.join(lines) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'verdictA', 'B', 'C', 'Cscore', 'D1', 'D2', 'table'), required=True)
    ap.add_argument('--workers', type=int, default=8, help='C: worker processes')
    args = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = H.results_dir(EXPERIMENT)
    if args.part == 'A':
        part_a(out)
    elif args.part == 'verdictA':
        verdict_a(out)
    elif args.part == 'B':
        part_b(out)
    elif args.part == 'C':
        part_c(out, args.workers)
    elif args.part == 'Cscore':
        score_c(out)
    elif args.part == 'D1':
        part_d1(out)
    elif args.part == 'D2':
        part_d2(out)
    comparison_table(out)


if __name__ == '__main__':
    main()
