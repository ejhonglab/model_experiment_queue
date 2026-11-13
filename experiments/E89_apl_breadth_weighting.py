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

CONDS = E88.CONDS
FACTORS = E88.FACTORS
MULT_FACTORS = E88.MULT_FACTORS
GROUPS = E88.GROUPS


def _iu() -> Tuple[np.ndarray, np.ndarray]:
    return E88._iu()


def kc_trial_se(trials: np.ndarray) -> np.ndarray:
    """Per KC the standard error of its trial mean: the SD of its single trials around their per-odor means, pooled over
    the odors, over sqrt(repeats). trials: (repeats, KCs, odors)."""
    n = trials.shape[0]
    res = trials - trials.mean(axis=0, keepdims=True)
    var = (res ** 2).sum(axis=(0, 2)) / ((n - 1) * trials.shape[2])
    return np.sqrt(var / n)


def lc_grid(x_min: float, x_max: float) -> np.ndarray:
    """The deconvolution's grid: LC_STEP from -2 to 6 (offset by half a step, so that no grid point lies on a class edge
    given to two decimals), 0.25 beyond, covering [x_min - 0.5, x_max + 0.5]."""
    fine = np.arange(-2.0 + LC_STEP / 2, 6.0, LC_STEP)
    lo, hi = min(x_min - 0.5, -2.0), max(x_max + 0.5, 6.0)
    left = np.arange(-2.0 - 0.125, lo - 0.25, -0.25)[::-1]
    right = np.arange(6.0 + 0.125, hi + 0.25, 0.25)
    return np.concatenate([left, fine, right])


def lc_prior(xbar: np.ndarray, se: np.ndarray, iters: int = LC_NPMLE_ITERS) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per odor the NPMLE of the distribution of a KC's mean response (EM on `lc_grid`; each KC's trial-mean `se`,
    floored at half a grid step), and each KC's likelihood at each grid point. Returns (grid, g: odors x grid,
    lik: odors x KCs x grid; float32)."""
    grid = lc_grid(float(np.nanmin(xbar)), float(np.nanmax(xbar)))
    se = np.maximum(np.asarray(se, dtype=float), LC_STEP / 2)
    z = ((xbar.T[:, :, None] - grid[None, None, :]) / se[None, :, None]).astype(np.float32)   # odors, KCs, grid
    lik = np.exp(-0.5 * z * z) / se[None, :, None].astype(np.float32)
    del z
    g = np.full((xbar.shape[1], len(grid)), 1.0 / len(grid), dtype=np.float32)
    for _ in range(iters):
        num = lik * g[:, None, :]
        num /= num.sum(axis=2, keepdims=True) + np.float32(1e-30)
        g = num.mean(axis=1)
    return grid, g, lik


def lc_classes(prior: Tuple[np.ndarray, np.ndarray, np.ndarray], cuts: Sequence[float]) -> Tuple[np.ndarray, np.ndarray]:
    """From `lc_prior`'s output, the classes between consecutive `cuts`: (mass: odors x classes, A: odors x KCs x
    classes, each KC's likelihood of its observed mean under each class)."""
    grid, g, lik = prior
    bins = np.searchsorted(np.sort(np.asarray(cuts, dtype=float)), grid)
    M = len(cuts) + 1
    mass = np.stack([g[:, bins == m].sum(axis=1) for m in range(M)], axis=1).astype(float)
    A = np.stack([(lik[:, :, bins == m] * g[:, None, bins == m]).sum(axis=2) for m in range(M)], axis=2).astype(float)
    A /= np.maximum(mass[:, None, :], 1e-12)
    return mass, A


def lc_overlap(trials: np.ndarray, crit: float, edges: Sequence[float] = (), se: Optional[np.ndarray] = None,
               em_iters: int = LC_EM_ITERS, prior=None) -> Tuple[np.ndarray, np.ndarray]:
    """(c) latent classes after deconvolution: (O, f) of one recording's single trials (repeats, KCs, odors; graded z).
    Per odor the distribution of a KC's mean response with the trial noise deconvolved (`lc_prior`: NPMLE on a grid,
    each KC's own trial SD pooled over odors, `kc_trial_se`), cut into classes at `crit` and at `edges`; per pair the
    table of class weights by maximum likelihood (EM, all pairs at once, from the independent table); O (triu order) the
    share of KCs above `crit` for both over sqrt(f_a f_b), f per odor the deconvolved share above `crit`. `prior`:
    `lc_prior`'s output for these trials, when several edges are read on one recording."""
    trials = np.asarray(trials, dtype=float)
    if prior is None:
        prior = lc_prior(trials.mean(axis=0), kc_trial_se(trials) if se is None else se)
    cuts = sorted(set(float(e) for e in edges if abs(float(e) - float(crit)) > 1e-9) | {float(crit)})
    mass, A = lc_classes(prior, cuts)
    M = mass.shape[1]
    top = np.arange(M) > cuts.index(float(crit))
    f = mass[:, top].sum(axis=1)
    i, j = _iu()
    Aa, Ab = A[i], A[j]                                      # pairs, KCs, classes
    W = mass[i][:, :, None] * mass[j][:, None, :]            # pairs, classes, classes
    K = A.shape[1]
    AaT = np.transpose(Aa, (0, 2, 1))
    for _ in range(em_iters):
        d = np.einsum('pkm,pkm->pk', Aa @ W, Ab) + 1e-300
        W = W * (AaT @ (Ab / d[:, :, None])) / K
    n11 = W[:, top][:, :, top].sum(axis=(1, 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        O = n11 / np.sqrt(f[i] * f[j])
    return O, f


def hf_fit(nresp: np.ndarray, n_rep: int, iters: int = 1000, tol: float = 1e-10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per odor (columns of `nresp`, the number of a KC's `n_rep` repeats that responded) the binomial mixture of
    responders (share f, hit rate h) and the rest (false-alarm rate e) by EM, vectorized over odors; h >= e. Returns
    (f, h, e)."""
    n = np.asarray(nresp, dtype=float)
    n_o = n.shape[1]
    f = np.clip((n >= n_rep - 0.5).mean(axis=0) + 0.01, 0.01, 0.9)
    h = np.full(n_o, 0.8)
    e = np.full(n_o, 0.02)
    for _ in range(iters):
        with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
            l1 = f * h ** n * (1 - h) ** (n_rep - n)
            l0 = (1 - f) * e ** n * (1 - e) ** (n_rep - n)
            p = np.where(l1 + l0 > 0, l1 / (l1 + l0), (n > n_rep / 2).astype(float))
        sp, sq = p.sum(axis=0), (1 - p).sum(axis=0)
        f_new = p.mean(axis=0)
        h_new = np.where(sp > 0, (p * n).sum(axis=0) / (n_rep * np.maximum(sp, 1e-300)), h)
        e_new = np.where(sq > 0, ((1 - p) * n).sum(axis=0) / (n_rep * np.maximum(sq, 1e-300)), e)
        done = max(np.abs(f_new - f).max(), np.abs(h_new - h).max(), np.abs(e_new - e).max()) < tol
        f, h, e = f_new, h_new, e_new
        if done:
            break
    swap = h < e
    f = np.where(swap, 1 - f, f)
    h, e = np.where(swap, e, h), np.where(swap, h, e)
    return f, h, e


def hf_overlap(binary: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    """(d) a correction by each odor's agreement between repeats: (O, f, params) of one recording's thresholded single
    trials (repeats, KCs, odors; any value above 0 a response). Per odor two classes (responder or not) with a hit rate
    h and a false-alarm rate e, fitted with their share f to the numbers of repeats responding (`hf_fit`); per pair the
    share of KCs in both classes n11 solved from the mean product of the two odors over different repeats (the mean
    over KCs and over repeat pairs i != j of x_a,i x_b,j), which is n11 (h_a - e_a)(h_b - e_b) + f_a h_a e_b
    + f_b e_a h_b + (1 - f_a - f_b) e_a e_b; n11 kept within its bounds [max(0, f_a + f_b - 1), min(f_a, f_b)];
    O = n11 / sqrt(f_a f_b)."""
    x = (np.asarray(binary, dtype=float) > 0).astype(float)
    n_rep = x.shape[0]
    f, h, e = hf_fit(x.sum(axis=0), n_rep)
    S = x.sum(axis=0)                                        # KCs, odors
    K = x.shape[1]
    same = np.einsum('rka,rkb->ab', x, x)                    # sum over KCs and repeats of x_a,i x_b,i
    E_ab = (S.T @ S - same) / (K * n_rep * (n_rep - 1))
    i, j = _iu()
    base = f[i] * h[i] * e[j] + f[j] * e[i] * h[j] + (1 - f[i] - f[j]) * e[i] * e[j]
    with np.errstate(divide='ignore', invalid='ignore'):
        n11 = (E_ab[i, j] - base) / ((h[i] - e[i]) * (h[j] - e[j]))
        n11 = np.clip(n11, np.maximum(0.0, f[i] + f[j] - 1.0), np.minimum(f[i], f[j]))
        O = n11 / np.sqrt(f[i] * f[j])
    return O, f, dict(f=f, h=h, e=e)


def same_repeat_gap(trials: np.ndarray) -> float:
    """E81's same-repeat gap of one recording's single trials (repeats, KCs, odors) as given: the mean over the two-odor
    pairs of the single-trial cosine distance between different repeats minus that within the same repeat
    (`H.same_repeat_gap` of `H.distances_by_kind`, computed directly)."""
    C = H.trial_cosines(np.asarray(trials, dtype=float))
    m = H._repeat_masks(trials.shape[0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        dr = 1.0 - H._kind_means(C, m['different_repeat'])
        sr = 1.0 - H._kind_means(C, m['same_repeat'])
    i, j = _iu()
    return float(np.nanmean(dr[i, j] - sr[i, j]))


def recording_values(tm: np.ndarray, trials: np.ndarray, crit_tm: float, crit_st: float,
                     lc_edges: Optional[Dict[str, Sequence[float]]] = None) -> Dict[str, np.ndarray]:
    """One recording's per-pair values for every candidate: from the thresholded trial means (s_tm, O_tm, R_tm), from
    the thresholded single trials with the trial noise removed (s_nr, O_nr, R_nr), (c) O_lc:<edges> for each of
    `lc_edges` (graded single trials, cut at `crit_tm`), (d) O_hf (single trials at `crit_st`)."""
    out = {}
    s, O, R = E88.factors_of_matrix(E87.thresholded(tm, crit_tm))
    out.update(s_tm=s, O_tm=O, R_tm=R)
    th = E87.thresholded(trials, crit_st)
    s, O, R, _ = E88.factors_of_trials(th)
    out.update(s_nr=s, O_nr=O, R_nr=R)
    out['O_hf'] = hf_overlap(th)[0]
    if lc_edges:
        prior = lc_prior(np.asarray(trials, dtype=float).mean(axis=0), kc_trial_se(np.asarray(trials, dtype=float)))
        for name, edges in lc_edges.items():
            out[f'O_lc:{name}'] = lc_overlap(trials, crit_tm, edges, prior=prior)[0]
    return out


def cand_factors(v: Dict[str, np.ndarray], cand: str, lc_name: Optional[str] = None) -> E88.Factors:
    """A candidate's Factors from per-pair means (`v`: the means over recordings, pseudo-recordings or wirings of
    `recording_values`' keys): R from the trial means, s from the noise-removed estimate, O the candidate's ('gm' the
    geometric mean of the two means)."""
    if cand == 'tm':
        O = v['O_tm']
    elif cand == 'nr':
        O = v['O_nr']
    elif cand == 'gm':
        O = np.sqrt(v['O_tm'] * v['O_nr'])
    elif cand == 'lc':
        O = v[f'O_lc:{lc_name}']
    elif cand == 'hf':
        O = v['O_hf']
    else:
        raise ValueError(cand)
    return E88.Factors(*E88._sOR(np.asarray(v['s_nr'], dtype=float), np.asarray(O, dtype=float), np.asarray(v['R_tm'], dtype=float)))


def _mean_dicts(ds: Sequence[Dict[str, np.ndarray]], rows: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
    """The per-pair NaN-mean over `ds` (or the rows `rows` of them) of each key."""
    ds = list(ds) if rows is None else [ds[r] for r in rows]
    return {k: E88._mean_rows(np.stack([d[k] for d in ds])) for k in ds[0]}


# ---------------------------------------------------------------------------
# part A: the model's readings at the three noise levels

class ModelReading89(E88.ModelReading):
    """E88's `ModelReading` with a seed (SEEDS; E88's draws at E88.SEED) and the third noise level: a per-KC, per-repeat
    term common to the odors of a repeat (`eta`, its own stream), mixed with the odor-specific draws as sqrt(1 - share)
    xi + sqrt(share) eta, which keeps each entry's noise variance and the agreement between repeats of one odor. The share
    is set per condition, at reading 1b's noise multiplier, so that the pseudo-recordings' same-repeat gap on their
    single trials thresholded at the level's single-trial criterion (`same_repeat_gap`, mean over pseudo-recordings) is
    the flies' (`gap_target`, their single trials thresholded at E87's single-trial criterion); geometric bisection.

    Levels (`level_params`): L1 the reading's own noise (E88's readings 1 and 2: the trial means at 1 z, the single
    trials at E88's reading-2 criterion); L2 reading 1b's multiplier, independent between odors (E88's 1b and 2b); L3h,
    L3, L3x the multiplier with the share times 0.5, 1 and 1.5. At L2 and L3 each condition's trial-mean and single-trial
    criteria keep reading 1's response rate in that condition. The truth of every level is E88's truth of reading 1
    (the noise-free z at reading 1's response rates)."""

    def __init__(self, o, label: str = '', seed: int = E88.SEED, gap_target: Optional[Dict[str, float]] = None, **kw):
        super().__init__(o, label=label, **kw)
        self.seed = int(seed)
        if self.seed != E88.SEED:
            rng = np.random.default_rng(self.seed)
            for c in CONDS:
                self.xi[c] = rng.standard_normal((N_PSEUDO, N_TRIALS) + self.C[c].shape[1:])
                self.xim[c] = self.xi[c].mean(axis=1)
        rng2 = np.random.default_rng([self.seed, 89])
        self.eta = {c: rng2.standard_normal((N_PSEUDO, N_TRIALS, self.C[c].shape[1], 1)) for c in CONDS}
        self.gap_target = gap_target
        self.share: Optional[Dict[str, float]] = None
        self.share_info: Dict[str, dict] = {}

    def trials_l(self, c: str, m: float, share: float) -> np.ndarray:
        """(N_PSEUDO, N_TRIALS, KCs, odors) single trials in z at noise multiplier `m` with `share` of the noise variance
        common to the odors of a repeat (E88's `trials` at share 0)."""
        if share == 0:
            return self.trials(c, m)
        xi = np.sqrt(1.0 - share) * self.xi[c] + np.sqrt(share) * self.eta[c]
        return self.C[c][:, None] + (m * self.sd[c])[:, None, :, None] * xi

    def means_l(self, c: str, m: float, share: float) -> np.ndarray:
        """(N_PSEUDO, KCs, odors) trial means of `trials_l` (E88's `trial_means` at share 0)."""
        if share == 0:
            return self.trial_means(c, m)
        xim = np.sqrt(1.0 - share) * self.xim[c] + np.sqrt(share) * self.eta[c].mean(axis=1)
        return self.C[c] + (m * self.sd[c])[:, :, None] * xim

    def gap(self, c: str, m: float, share: float, rate: float) -> float:
        tr = self.trials_l(c, m, share)
        crit = E88._quantile_crit(tr, rate)
        return float(np.mean([same_repeat_gap(E87.thresholded(t, crit)) for t in tr]))

    def shares(self) -> Dict[str, float]:
        """L3's share per condition (bisection on [SHARE_BOUNDS], 30 steps, on the gap at reading 1b's multiplier and
        reading 1's rate)."""
        if self.share is None:
            if self.gap_target is None:
                raise ValueError(f'{self.label}: level 3 needs gap_target')
            mult = self.multipliers()
            self.share = {}
            for c in CONDS:
                rate, tgt = self.rate['r1'][c], float(self.gap_target[c])
                lo, hi = SHARE_BOUNDS
                g_lo, g_hi = self.gap(c, mult[c], lo, rate), self.gap(c, mult[c], hi, rate)
                at = ''
                if g_lo >= tgt:
                    s, at = lo, 'low'
                elif g_hi <= tgt:
                    s, at = hi, 'high'
                else:
                    for _ in range(30):
                        mid = 0.5 * (lo + hi)
                        if self.gap(c, mult[c], mid, rate) < tgt:
                            lo = mid
                        else:
                            hi = mid
                    s = 0.5 * (lo + hi)
                self.share[c] = float(s)
                self.share_info[c] = dict(share=float(s), at_bound=at, gap_target=tgt, gap_share0=g_lo,
                                          gap_at_share=self.gap(c, mult[c], s, rate))
        return self.share

    def level_params(self, level: str, c: str) -> dict:
        """{m, share, crit_tm, crit_st} of a level in condition `c`."""
        self.factors('r1')
        if level == 'L1':
            self.factors('r2')
            return dict(m=1.0, share=0.0, crit_tm=self.crit['r1'][c], crit_st=self.crit['r2'][c])
        m = self.multipliers()[c]
        share = 0.0 if level == 'L2' else min(SHARE_BOUNDS[1], SHARE_SCALE[level] * self.shares()[c])
        rate = self.rate['r1'][c]
        return dict(m=m, share=share, crit_tm=E88._quantile_crit(self.means_l(c, m, share), rate),
                    crit_st=E88._quantile_crit(self.trials_l(c, m, share), rate))

    def level_values(self, level: str, c: str, lc_edges: Optional[Dict[str, Sequence[float]]] = None
                     ) -> Tuple[Dict[str, np.ndarray], dict]:
        """(the means over the pseudo-recordings of `recording_values`, the level's parameters with the achieved rates,
        gap and reliability) in condition `c`."""
        pr = self.level_params(level, c)
        tr = self.trials_l(c, pr['m'], pr['share'])
        tm = self.means_l(c, pr['m'], pr['share'])
        vals = [recording_values(tm[r], tr[r], pr['crit_tm'], pr['crit_st'], lc_edges) for r in range(N_PSEUDO)]
        info = dict(pr, rate_tm=float((tm > pr['crit_tm']).mean()), rate_st=float((tr > pr['crit_st']).mean()),
                    gap=float(np.mean([same_repeat_gap(E87.thresholded(t, pr['crit_st'])) for t in tr])),
                    gap_graded=float(np.mean([same_repeat_gap(t) for t in tr])))
        mask = H._repeat_masks(N_TRIALS)['different_repeat']
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            info['reliability'] = float(np.mean([np.nanmean(np.diag(H._kind_means(H.trial_cosines((t > pr['crit_st']).astype(float)), mask)))
                                                 for t in tr]))
        return _mean_dicts(vals), info

    def truth_values(self) -> Dict[str, E88.Factors]:
        return self.truth('r1')


def flies_gap_target(unit: str = 'z') -> Dict[str, dict]:
    """The flies' same-repeat gap per condition (mean over recordings): on their single trials thresholded at E87's
    single-trial criterion (the target of L3) and graded (E81's 0.025 / 0.070, beside)."""
    key = ('gap', unit)
    if key not in _CACHE:
        tr = E87.trial_arrays(unit)
        c1 = E87.single_trial_criterion(unit)
        out = {}
        for c, g in (('on', 'TNT3c'), ('off', 'TNT')):
            out[c] = dict(thresholded=float(np.mean([same_repeat_gap(E87.thresholded(a, c1)) for a in tr[g].values()])),
                          graded=float(np.mean([same_repeat_gap(a) for a in tr[g].values()])))
        _CACHE[key] = out
    return _CACHE[key]


# ---------------------------------------------------------------------------
# part A: the flies

def flies_values(lc_edges: Dict[str, Sequence[float]], unit: str = 'z') -> Dict[str, dict]:
    """{'on' / 'off': dict(keys, vals: [recording_values per recording])} of the flies (control flies / APL output
    blocked (TNT)): trial means at the trial-mean criterion, single trials at the single-trial criterion (E87's)."""
    key = ('flies_vals', unit, tuple(sorted(lc_edges)))
    if key not in _CACHE:
        tm, tr = E87.trialmean_mats(unit), E87.trial_arrays(unit)
        crit, c1 = E87.criterion(unit), E87.single_trial_criterion(unit)
        out = {}
        for c, g in (('on', 'TNT3c'), ('off', 'TNT')):
            keys = list(tr[g])
            if list(tm[g]) != keys:
                raise RuntimeError('the trial means and the single trials list the recordings differently')
            out[c] = dict(keys=keys, vals=[recording_values(tm[g][k], tr[g][k], crit, c1, lc_edges) for k in keys])
        _CACHE[key] = out
    return _CACHE[key]


def flies_cand(fv: Dict[str, dict], cand: str, lc_name: Optional[str] = None, tnt: Optional[Sequence[tuple]] = None,
               rows: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> Dict[str, E88.Factors]:
    """The flies' factors with a candidate: all recordings (or the TNT recordings `tnt`), or the resampled `rows`."""
    on, off = fv['on']['vals'], fv['off']['vals']
    if tnt is not None:
        want = [tuple(k) for k in tnt]
        off = [off[fv['off']['keys'].index(k)] for k in want]
    if rows is None:
        return {'on': cand_factors(_mean_dicts(on), cand, lc_name), 'off': cand_factors(_mean_dicts(off), cand, lc_name)}
    return {'on': cand_factors(_mean_dicts(on, rows[0]), cand, lc_name), 'off': cand_factors(_mean_dicts(off, rows[1]), cand, lc_name)}


def flies_cand_boot(fv: Dict[str, dict], cand: str, lc_name: Optional[str] = None, tnt=None, n: int = E88.N_BOOT,
                    seed: int = 0) -> List[Dict[str, E88.Factors]]:
    n_off = len(fv['off']['vals']) if tnt is None else len(tnt)
    return [flies_cand(fv, cand, lc_name, tnt, r) for r in E88.boot_rows(len(fv['on']['vals']), n_off, n, seed)]


# ---------------------------------------------------------------------------
# part A: the populations, read in a worker per run and seed

def a_runs() -> Dict[str, List[Tuple[str, tuple, Optional[float]]]]:
    """{population: [(wiring, (folder, tag), K or None)]} of A_POPULATIONS (E88's `model_populations`: 'b2' E67's fit on
    E72 part B's eight wirings, 'fc' the firing-ceiling family, 'S_unchanged' b2 on S, 'base' E87's D2 runs with
    per-odor APL scales on S)."""
    pops = E88.model_populations()
    out = {}
    for pop in A_POPULATIONS:
        src = pops['D2_thrapl' if pop == 'base' else pop]
        out[pop] = [(w, key, K) for w, key, _, K in src]
    return out


def _load_run(key: tuple) -> Dict[str, tuple]:
    folder, tag = key
    if folder == 'E87_inhibition_form':
        o = E87._load_d2(tag)
        if o is None:
            raise RuntimeError(f'no D2 run {tag}')
        return o
    return E87._load(folder, tag)


def _a_job(job: tuple) -> dict:
    """One run at one seed: every level's per-pair means of every candidate, the truth, the level parameters."""
    key, label, K, seed, rel, gap, lc_edges = job
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    warnings.simplefilter('ignore')
    t0 = time.time()
    mr = ModelReading89(_load_run(key), label=label, seed=seed, K=K, rel_target=rel, gap_target=gap)
    res = dict(key=key, label=label, seed=seed, levels={}, info={}, cal=dict(mr.cal))
    for level in LEVELS:
        res['levels'][level] = {}
        for c in CONDS:
            v, info = mr.level_values(level, c, lc_edges)
            res['levels'][level][c] = v
            res['info'][(level, c)] = dict(info, rate_r1=mr.rate['r1'][c])
    tr = mr.truth_values()
    res['truth'] = {c: dict(s=tr[c].s, O=tr[c].O, R=tr[c].R) for c in CONDS}
    res['mult'] = dict(mr.mult_info)
    res['share'] = dict(mr.share_info)
    res['seconds'] = time.time() - t0
    return res


def run_a_jobs(runs: Sequence[tuple], seeds: Sequence[int], lc_edges: Dict[str, Sequence[float]], n_workers: int
               ) -> Dict[Tuple[tuple, int], dict]:
    """`_a_job` for every (run, seed), in `n_workers` processes; {(key, seed): result}."""
    from concurrent.futures import ProcessPoolExecutor
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in flies_gap_target('z').items()}
    jobs = [(key, label, K, s, rel, gap, lc_edges) for label, key, K in runs for s in seeds]
    out = {}
    t0 = time.time()
    if n_workers <= 1:
        for j in jobs:
            r = _a_job(j)
            out[(r['key'], r['seed'])] = r
            print(f"A {r['label']} seed {r['seed']}: {r['seconds']:.0f}s", flush=True)
        return out
    import multiprocessing
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
        for r in ex.map(_a_job, jobs):
            out[(r['key'], r['seed'])] = r
            print(f"A {r['label']} {r['key'][1][:40]} seed {r['seed']}: {r['seconds']:.0f}s ({time.time() - t0:.0f}s)", flush=True)
    return out


def pop_values(res: Dict[Tuple[tuple, int], dict], members: Sequence[tuple], seed: int, level: str
               ) -> Dict[str, Dict[str, np.ndarray]]:
    """A population's per-pair means (over its wirings, each the mean over its pseudo-recordings) per condition."""
    return {c: _mean_dicts([res[(key, seed)]['levels'][level][c] for _, key, _ in members]) for c in CONDS}


def pop_truth(res, members, seed: int) -> Dict[str, E88.Factors]:
    return {c: E88.Factors(*(E88._mean_rows(np.stack([res[(key, seed)]['truth'][c][k] for _, key, _ in members]))
                             for k in ('s', 'O', 'R'))) for c in CONDS}


def check_rows(est: Dict[str, E88.Factors], tru: Dict[str, E88.Factors], **meta) -> List[dict]:
    """E88's `estimator_check` rows for one candidate (estimate minus truth per factor, group and condition; R, L and G in
    units of the similarity, s_true * ln(est / true))."""
    rows = []
    masks = E88.group_masks()
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
                er = E88._nanmean(err[m])
                rows.append(dict(meta, condition=c, factor=f, group=g, n_pairs=int(m.sum()), estimate=E88._nanmean(e[m]),
                                 truth=E88._nanmean(t[m]), diff=E88._nanmean(diff[m]), err=er,
                                 off=bool(g != 'all' and np.isfinite(er) and abs(er) > EST_TOL)))
    return rows


def estimator_check(res: Dict[Tuple[tuple, int], dict], runs: Dict[str, list], lc_names: Sequence[str],
                    populations: Sequence[str]) -> pd.DataFrame:
    """Each candidate's factors (the population's) minus the truth, by population, seed, level, condition, factor and
    group: 'err' in units of the similarity (E88's convention), 'off' beyond EST_TOL."""
    rows = []
    for pop in populations:
        mem = runs[pop]
        for seed in SEEDS:
            if not all((key, seed) in res for _, key, _ in mem):
                continue
            tru = pop_truth(res, mem, seed)
            for level in LEVELS:
                v = pop_values(res, mem, seed, level)
                for cand in CANDIDATES:
                    for lc in (lc_names if cand == 'lc' else (None,)):
                        est = {c: cand_factors(v[c], cand, lc) for c in CONDS}
                        rows += check_rows(est, tru, population=pop, seed=seed, level=level, candidate=cand,
                                           lc_edges=lc or '', binding=level in BINDING_LEVELS)
    return pd.DataFrame(rows)


def worst_errors(ec: pd.DataFrame) -> pd.DataFrame:
    """Per population, candidate (and edges), level and factor the largest |group-mean error| over the five groups, both
    conditions and both seeds."""
    g = ec[ec.group != 'all'].assign(abs_err=lambda d: d.err.abs())
    return g.groupby(['population', 'candidate', 'lc_edges', 'level', 'factor']).abs_err.max().reset_index()


def lc_choice(ec: pd.DataFrame) -> dict:
    """(c)'s free choice on the firing-ceiling population: the edges with the smallest largest error of L over the binding
    levels (then of R, L, s and G together)."""
    w = worst_errors(ec[(ec.population == LC_SET_POPULATION) & (ec.candidate == 'lc')])
    w = w[w.level.isin(BINDING_LEVELS)]
    score = {}
    for name, g in w.groupby('lc_edges'):
        score[name] = (float(g[g.factor == 'L'].abs_err.max()), float(g[g.factor.isin(FULL_FACTORS)].abs_err.max()))
    best = min(score, key=lambda k: score[k])
    return dict(chosen=best, population=LC_SET_POPULATION, worst_L_binding={k: v[0] for k, v in score.items()},
                worst_RLsG_binding={k: v[1] for k, v in score.items()}, edges={k: list(LC_EDGES[k]) for k in score})


def choose_candidate(ec: pd.DataFrame, population: str, lc_name: str) -> dict:
    """The entry's thresholds on one population: a candidate passes for who responds if R, L and s are inside EST_TOL
    for every group, both conditions and both seeds at every binding level, in full if G is as well; the chosen one is a
    candidate that passes in full if any, else one that passes for who responds, with the smallest largest error of L at
    L3 (both seeds). Level 2 is given beside and does not bind."""
    w = worst_errors(ec[(ec.population == population) & ((ec.candidate != 'lc') | (ec.lc_edges == lc_name))])
    out = dict(population=population, tolerance=EST_TOL, candidates={})
    for cand in CANDIDATES:
        g = w[w.candidate == cand]
        if g.empty:
            continue
        worst = {lv: {f: float(g[(g.level == lv) & (g.factor == f)].abs_err.max()) for f in FACTORS} for lv in LEVELS}
        bind = g[g.level.isin(BINDING_LEVELS)]
        who = bool((bind[bind.factor.isin(WHO_FACTORS)].abs_err <= EST_TOL).all())
        full = bool(who and (bind[bind.factor == 'G'].abs_err <= EST_TOL).all())
        l2 = g[(g.level == 'L2') & g.factor.isin(WHO_FACTORS)]
        out['candidates'][cand] = dict(name=CANDIDATE_NAMES[cand], passes_who=who, passes_full=full,
                                       passes_level2_who=bool((l2.abs_err <= EST_TOL).all()),
                                       worst_L_L3=worst['L3']['L'], worst=worst,
                                       worst_binding={f: float(bind[bind.factor == f].abs_err.max()) for f in FACTORS})
    c = out['candidates']
    full = [k for k in c if c[k]['passes_full']]
    who = [k for k in c if c[k]['passes_who']]
    pool = full or who
    out['chosen'] = min(pool, key=lambda k: c[k]['worst_L_L3']) if pool else None
    out['chosen_passes'] = 'in full' if full else ('for who responds' if who else 'none')
    return out


# ---------------------------------------------------------------------------
# part A: A1 and A2 with the chosen candidate

def _pop_cand(res, mem, seed, level, cand, lc) -> Dict[str, E88.Factors]:
    v = pop_values(res, mem, seed, level)
    return {c: cand_factors(v[c], cand, lc) for c in CONDS}


def part_a1(out: Path, res, runs, fv, cands: Sequence[str], lc: str, seed: int = SEEDS[0]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A1: factors.csv (the flies with every candidate in `cands`, all recordings and the six typical, 90 % intervals over
    E88.N_BOOT resamples; each population and wiring at every level) and cells.csv (each population's cells at every
    level against the flies' factors with the same candidate, all and the six typical, with E88's floors)."""
    frows, crows = [], []
    for cand in cands:
        for subset in E88.SUBSETS:
            tnt = E88.flies_subset(subset)
            Ff = flies_cand(fv, cand, lc, tnt)
            boot = flies_cand_boot(fv, cand, lc, tnt)
            frows += [dict(source='flies', candidate=cand, subset=subset, population='flies', member=f'mean_{subset}', level='',
                           **r) for r in E88._factor_rows(Ff, boot)]
            for pop, mem in runs.items():
                for level in LEVELS:
                    Fm = _pop_cand(res, mem, seed, level, cand, lc)
                    crows += E88._cell_rows(Fm, Ff, boot, population=pop, level=level, subset=subset, candidate=cand, seed=seed)
        for pop, mem in runs.items():
            for level in LEVELS:
                for w, key, _ in mem:
                    v = res[(key, seed)]['levels'][level]
                    F = {c: cand_factors(v[c], cand, lc) for c in CONDS}
                    frows += [dict(source='model', candidate=cand, subset='', population=pop, member=w, level=level, **r)
                              for r in E88._factor_rows(F)]
                frows += [dict(source='model', candidate=cand, subset='', population=pop, member='population', level=level, **r)
                          for r in E88._factor_rows(_pop_cand(res, mem, seed, level, cand, lc))]
                frows += [dict(source='truth', candidate='', subset='', population=pop, member='population', level=level, **r)
                          for r in E88._factor_rows(pop_truth(res, mem, seed))] if (cand == cands[0] and level == 'L1') else []
        print(f'A1 {cand}', flush=True)
    fdf, cdf = pd.DataFrame(frows), pd.DataFrame(crows)
    fdf.to_csv(out / 'factors.csv', index=False)
    cdf.to_csv(out / 'cells.csv', index=False)
    return fdf, cdf


def part_a2(out: Path, res, runs, fv, cands: Sequence[str], lc: str, seed: int = SEEDS[0]) -> pd.DataFrame:
    """A2: E88's B1 per candidate: per group the change of each factor under APL (the mean over the group's pairs of
    log(on / off)) for the flies (90 % intervals over E88.N_BOOT resamples; all and the six typical) and each population
    at every level, with the known truth and the candidate's error on it. changes.csv."""
    rows = []
    for cand in cands:
        for subset in E88.SUBSETS:
            tnt = E88.flies_subset(subset)
            ch = E88.log_changes(flies_cand(fv, cand, lc, tnt))
            bch = [E88.log_changes(b) for b in flies_cand_boot(fv, cand, lc, tnt)]
            for k, v in ch.items():
                _, f, g = k.split('_', 2)
                lo, hi = np.nanpercentile([b[k] for b in bch], E88.INTERVAL)
                rows.append(dict(source='flies', candidate=cand, population='flies', member='mean', subset=subset, level='',
                                 factor=f, group=g, log_change=v, ratio=float(np.exp(v)), lo=float(lo), hi=float(hi)))
        for pop, mem in runs.items():
            tch = E88.log_changes(pop_truth(res, mem, seed))
            for level in LEVELS:
                ch = E88.log_changes(_pop_cand(res, mem, seed, level, cand, lc))
                per = [E88.log_changes({c: cand_factors(res[(key, seed)]['levels'][level][c], cand, lc) for c in CONDS})
                       for _, key, _ in mem]
                for k, v in ch.items():
                    _, f, g = k.split('_', 2)
                    ws = [p[k] for p in per]
                    rows.append(dict(source='model', candidate=cand, population=pop, member='population', subset='', level=level,
                                     factor=f, group=g, log_change=v, ratio=float(np.exp(v)), lo=float(np.nanmin(ws)),
                                     hi=float(np.nanmax(ws)), truth=tch[k], truth_ratio=float(np.exp(tch[k])), error=v - tch[k]))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'changes.csv', index=False)
    return df


def b1_differs(ch: pd.DataFrame, cand: str, pop: str, group: str, factor: str, subset: str = 'all',
               levels: Sequence[str] = BINDING_LEVELS) -> dict:
    """E88's B1 rule with one candidate: at each level of `levels` the population's change against the flies' interval;
    differs if outside it and by a factor of B1_FACTOR or more at every level."""
    fl = ch[(ch.source == 'flies') & (ch.candidate == cand) & (ch.subset == subset) & (ch.factor == factor) & (ch.group == group)].iloc[0]
    out = dict(flies=float(fl.log_change), flies_ratio=float(np.exp(fl.log_change)), flies_lo=float(fl.lo), flies_hi=float(fl.hi))
    every = True
    for lv in levels:
        m = ch[(ch.source == 'model') & (ch.candidate == cand) & (ch.population == pop) & (ch.level == lv)
               & (ch.factor == factor) & (ch.group == group)].iloc[0]
        x = float(m.log_change)
        outside = not (fl.lo <= x <= fl.hi)
        big = abs(x - fl.log_change) >= np.log(B1_FACTOR)
        out[lv] = dict(model=x, model_ratio=float(np.exp(x)), truth_ratio=float(m.truth_ratio), error=float(m.error),
                       outside=bool(outside), by_factor=bool(big), differs=bool(outside and big))
        every &= outside and big
    out['differs'] = bool(every)
    return out


def verdict_a(out: Path, choice: dict, cands_read: Sequence[str], lc: str) -> dict:
    """The thresholds of A on the chosen candidate (estimator_choice.json): the split on the breadth-2 population read
    only with a candidate that passes in full there and only if the verdict is the same at every binding level; the
    groups (recruitment / co-response) with a candidate that passes for who responds; the base (the co-response is what
    is left if the cell with the flies' L alone is lower by D2_BASE_GAIN or more) with a candidate that passes for who
    responds on the base; every verdict repeated with the six typical recordings; A2's changes by E88's B1 rule.
    verdict_A.json."""
    cells = pd.read_csv(out / 'cells.csv')
    ch = pd.read_csv(out / 'changes.csv')
    v = dict(choice={p: dict(chosen=c['chosen'], passes=c['chosen_passes']) for p, c in choice['populations'].items()},
             lc_edges=lc, thresholds=dict(EST_TOL=EST_TOL, G_MIN=G_MIN, G_RATIO=G_RATIO, GROUP_RATIO=GROUP_RATIO,
                                          D2_BASE_GAIN=D2_BASE_GAIN, B1_FACTOR=B1_FACTOR))

    def cell(cand, pop, level, subset, name):
        return cells[(cells.candidate == cand) & (cells.population == pop) & (cells.level == level) & (cells.subset == subset)
                     & (cells.cell == name)].iloc[0]

    b2 = choice['populations']['b2']
    cb2 = b2['chosen']
    v['outcome_d_no_candidate'] = cb2 is None
    split, groups, base = {}, {}, {}
    base_c = choice['populations']['base']['chosen']
    for subset in E88.SUBSETS:
        if cb2 is not None:
            sp = {lv: E88._split_verdict(*(float(cell(cb2, 'b2', lv, subset, n).pattern) for n in ('own', 'sets', 'amp')))
                  for lv in LEVELS}
            vals = [sp[lv]['verdict'] for lv in BINDING_LEVELS]
            agree = all(x == vals[0] for x in vals)
            sp['verdict'] = (vals[0] if agree else 'the binding levels disagree') if b2['chosen_passes'] == 'in full' else (
                'not read: the chosen candidate does not pass in full' if b2['chosen_passes'] == 'for who responds' else 'not read')
            sp['verdict_if_read'] = vals[0] if agree else 'the binding levels disagree'
            split[subset] = sp
            gr = {}
            for g in GROUPS:
                gg = {lv: E88._group_verdict(*(float(cell(cb2, 'b2', lv, subset, n)[f'err_{g}']) for n in ('own', 'R', 'L')))
                      for lv in LEVELS}
                vals = [gg[lv]['verdict'] for lv in BINDING_LEVELS]
                gg['verdict'] = vals[0] if all(x == vals[0] for x in vals) else 'the binding levels disagree'
                gr[g] = gg
            groups[subset] = gr
        if base_c is not None:
            dd = {}
            for lv in LEVELS:
                own = float(cell(base_c, 'base', lv, subset, 'own').pattern)
                lcl = float(cell(base_c, 'base', lv, subset, 'L').pattern)
                dd[lv] = dict(own=own, L_cell=lcl, gain=own - lcl, coresponse_left=bool(own - lcl >= D2_BASE_GAIN),
                              sets_cell=float(cell(base_c, 'base', lv, subset, 'sets').pattern),
                              unchanged_own=float(cell(base_c, 'S_unchanged', lv, subset, 'own').pattern))
            vals = ['the co-response is left' if dd[lv]['coresponse_left'] else 'not the co-response alone' for lv in BINDING_LEVELS]
            dd['verdict'] = vals[0] if all(x == vals[0] for x in vals) else 'the binding levels disagree'
            base[subset] = dd
    v.update(split=split, groups=groups, base=dict(candidate=base_c, **base))
    changed = []
    if split and split['six']['verdict'] != split['all']['verdict']:
        changed.append(f"split: {split['all']['verdict']} with all, {split['six']['verdict']} with the six typical")
    for g in (GROUPS if groups else ()):
        if groups['six'][g]['verdict'] != groups['all'][g]['verdict']:
            changed.append(f"{g}: {groups['all'][g]['verdict']} with all, {groups['six'][g]['verdict']} with the six typical")
    if base and base['six']['verdict'] != base['all']['verdict']:
        changed.append(f"base: {base['all']['verdict']} with all, {base['six']['verdict']} with the six typical")
    v['depends_on_the_two_strongest_recordings'] = changed
    a2 = {}
    for cand in cands_read:
        a2[cand] = {}
        for subset in E88.SUBSETS:
            a2[cand][subset] = {pop: {g: {f: b1_differs(ch, cand, pop, g, f, subset) for f in ('L', 'R')}
                                      for g in ('nearest4', 'ester', 'alcohol')} for pop in A_POPULATIONS}
    v['A2'] = a2
    v['A2_candidates'] = list(cands_read)
    _write_json(out / 'verdict_A.json', v)
    return v


def part_a(out: Path, n_workers: int = 4) -> dict:
    """A0 (the check on known truth, both seeds, every level; (c)'s edges set on the firing-ceiling population first),
    the choice per population, A1 and A2 with the chosen candidate (and E88's two readings and the geometric mean beside),
    the verdict. estimator_check.csv, estimator_choice.json, levels.csv, factors.csv, cells.csv, changes.csv,
    verdict_A.json."""
    provenance(out, 'A', workers=n_workers, flies_gap=flies_gap_target('z'), flies_reliability=E88.flies_reliability('z'))
    runs = a_runs()
    t0 = time.time()
    # (c)'s edges on the firing-ceiling population, every choice
    res = run_a_jobs([(w, key, K) for w, key, K in runs[LC_SET_POPULATION]], SEEDS, LC_EDGES, n_workers)
    ec_fc = estimator_check(res, runs, list(LC_EDGES), [LC_SET_POPULATION])
    lcc = lc_choice(ec_fc)
    lc = lcc['chosen']
    print(f'A0: (c) edges set on {LC_SET_POPULATION}: {lc} ({time.time() - t0:.0f}s)', lcc, flush=True)
    seen = {key for _, key, _ in runs[LC_SET_POPULATION]}
    rest = []
    for pop, mem in runs.items():
        for w, key, K in mem:
            if key not in seen:
                seen.add(key)
                rest.append((w, key, K))
    res2 = run_a_jobs(rest, SEEDS, {lc: LC_EDGES[lc]}, n_workers)
    for k, r in res.items():                 # the firing-ceiling runs keep only the chosen edges' (c)
        for level in r['levels']:
            for c in CONDS:
                r['levels'][level][c] = {kk: vv for kk, vv in r['levels'][level][c].items()
                                         if not kk.startswith('O_lc:') or kk == f'O_lc:{lc}'}
    res.update(res2)
    ec = estimator_check(res, runs, [lc], list(A_POPULATIONS))
    ec_fc_other = ec_fc[(ec_fc.candidate == 'lc') & (ec_fc.lc_edges != lc)]
    pd.concat([ec, ec_fc_other], ignore_index=True).to_csv(out / 'estimator_check.csv', index=False)
    lv_rows = []
    for (key, seed), r in res.items():
        for (level, c), info in r['info'].items():
            lv_rows.append(dict(folder=key[0], tag=key[1], wiring=r['label'], seed=seed, level=level, condition=c,
                                **{k: x for k, x in info.items() if np.isscalar(x)}, seconds=r['seconds']))
    pd.DataFrame(lv_rows).to_csv(out / 'levels.csv', index=False)
    choice = dict(lc_edges=lcc, populations={p: choose_candidate(ec, p, lc) for p in A_POPULATIONS},
                  binding_levels=list(BINDING_LEVELS), seeds=list(SEEDS), flies_gap=flies_gap_target('z'),
                  flies_reliability=E88.flies_reliability('z'))
    _write_json(out / 'estimator_choice.json', choice)
    for p in CHECK_POPULATIONS:
        print(f"A0 {p}: chosen {choice['populations'][p]['chosen']} ({choice['populations'][p]['chosen_passes']})", flush=True)
    fv = flies_values({lc: LC_EDGES[lc]})
    chosen = [choice['populations'][p]['chosen'] for p in CHECK_POPULATIONS]
    cands = list(dict.fromkeys([c for c in chosen if c] + ['tm', 'nr', 'gm']))
    part_a1(out, res, runs, fv, cands, lc)
    part_a2(out, res, runs, fv, cands, lc)
    v = verdict_a(out, choice, cands, lc)
    _CACHE['a_res'] = res
    print(f'A done ({time.time() - t0:.0f}s)', flush=True)
    return v


# ---------------------------------------------------------------------------
# part B: the weighting in the saved models (no simulations)

def b_families() -> Dict[str, dict]:
    """The families in which only the weighting differs: {family: dict(folder, members: {tag: (measure, exponent)},
    ref: the exponent-2 member's tag, operating point)}. E66's part C per reference model (spike total at exponents 0 to
    3, number of odors at 1 to 3, the APL-on rate re-held; the reference is the fitted spike total at 2); E67's fits of
    C2's family on the real wiring re-fitted under the calibrated readout at exponents 1, 1.5 and 3, with C2 itself
    (rate0.068_aplboth_k2_cap3) as the exponent-2 member and E67's readout-only fit at 2 beside; E12's grid (spike total
    and number of odors, the 5 % operating point) and E19's per threshold and KC-KC scale (5 %), given apart."""
    fam = {}
    for ref in ('thr40_gain1.5', 'thr70_gain1.5', 'thr70_gain1.25'):
        mem = {f'C_{ref}_spikes_k{k}': ('spikes', float(k)) for k in (0, 1, 2, 3)}
        mem.update({f'C_{ref}_n_odors_k{k}': ('n_odors', float(k)) for k in (1, 2, 3)})
        fam[f'E66_{ref}'] = dict(folder='E66_per_odor_apl_profile', members=mem, ref=f'C_{ref}_spikes_k2', operating_point=0.068)
    fam['E67_C2'] = dict(folder='E67_apl_effect_fits', ref='rate0.068_aplboth_k2_cap3', operating_point=0.068,
                         members={'rate0.068_aplcal_k1_cap3_c2': ('spikes', 1.0), 'rate0.068_aplcal_k1.5_cap3_c2': ('spikes', 1.5),
                                  'rate0.068_aplboth_k2_cap3': ('spikes', 2.0), 'rate0.068_aplcal_k2_cap3': ('spikes', 2.0),
                                  'rate0.068_aplcal_k3_cap3': ('spikes', 3.0)})
    mem = {'k0': ('none', 0.0)}
    mem.update({f'spikes_k{k:g}': ('spikes', float(k)) for k in (-1, 0.5, 1, 2, 3, 4)})
    mem.update({f'n_odors_k{k:g}': ('n_odors', float(k)) for k in (-1, 0.5, 1, 2)})
    fam['E12_grid'] = dict(folder='E12_breadth_apl', members=mem, ref='spikes_k2', operating_point=0.05)
    for thr in (40, 70):
        for kk in (0, 30, 100):
            fam[f'E19_thr{thr}_kckc{kk}'] = dict(folder='E19_thr_breadth_kckc', ref=f'thr{thr}_k2_kckc{kk}', operating_point=0.05,
                                                  members={f'thr{thr}_k{k}_kckc{kk}': ('spikes', float(k)) for k in (0, 1, 2, 3)})
    return fam


def _archive_paths() -> Dict[Tuple[str, str], Tuple[str, str]]:
    if 'archive' not in _CACHE:
        _CACHE['archive'] = {(f, v): (str(on), str(off)) for f, v, on, off in E87.archive_variants()}
    return _CACHE['archive']


def _b_reading(on: str, off: str) -> dict:
    """The main target of one archived variant under the calibrated readout (K set on its APL-on counts, E83's 'cal')
    and on raw counts (E83's `effect_scores`)."""
    from experiments import E83_cross_wiring_fit as E
    o = {'on': (None, pd.read_parquet(on)), 'off': (None, pd.read_parquet(off))}
    row = {}
    for ro in ('cal', 'raw'):
        sc = E.effect_scores(E.pair_distances(o, ro))
        row.update({f'{ro}_{k}': sc[k] for k in ('z8_pattern', 'z8_loss', 'mean_effect', 'z8_rho', 'z8_err_ester', 'z8_err_alcohol',
                                                  'z6_pattern', 'noise_free_pattern')})
    return row


def _b_candidate_job(job: tuple) -> dict:
    """One archived variant read with part A's candidate at L1 (its own calibration) and L3 (the middle share): the
    changes of L and R under APL per group (E88's `log_changes`)."""
    folder, variant, on, off, cand, lc, rel, gap = job
    warnings.simplefilter('ignore')
    try:
        o = {'on': (None, pd.read_parquet(on)), 'off': (None, pd.read_parquet(off))}
        mr = ModelReading89(o, label=variant, rel_target=rel, gap_target=gap)
        row = dict(folder=folder, variant=variant)
        for level in ('L1', 'L3'):
            v = {c: mr.level_values(level, c, {lc: LC_EDGES[lc]} if cand == 'lc' else None)[0] for c in CONDS}
            ch = E88.log_changes({c: cand_factors(v[c], cand, lc) for c in CONDS}, ('L', 'R'))
            row.update({f'{level}_{k}': x for k, x in ch.items() if k.split('_', 2)[2] in ('nearest4', 'ester', 'alcohol')})
        return row
    except Exception as e:
        return dict(folder=folder, variant=variant, error=f'{type(e).__name__}: {e}'[:300])


def part_b(out: Path, n_workers: int = 4) -> dict:
    """The weighting in the saved models: per family member the main target under the readout and on raw counts (from
    the archived counts), results/apl_effect/by_variant.csv's scores and realism flags, E88's archive_scan.csv rows (the
    cells, the group errors and the changes of the co-response and of the recruitment, as trial means 'r1' and on the
    deterministic counts 'r3raw'), and, with a candidate of part A that passes for who responds on the breadth-2
    population, the changes read again with it (L1 and L3). The rule per family, fixed beforehand: the pattern part
    follows the alcohols' co-response change if the member whose change is nearest the flies' (trial means, like for
    like; on the deterministic counts and with the candidate beside) has a pattern part lower than the family's
    exponent-2 member by B_GAIN or more, under the readout and on raw counts. weighting_archive.csv, verdict_B.json."""
    from concurrent.futures import ProcessPoolExecutor
    provenance(out, 'B', workers=n_workers)
    fams = b_families()
    paths = _archive_paths()
    scan = pd.read_csv(H.results_dir(E88.EXPERIMENT) / 'archive_scan.csv', low_memory=False)
    bv = pd.read_csv(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv')
    fl = {}
    for frd in ('r1', 'r2'):
        rf = E88.flies_factors('z', frd)
        fl[frd] = E88.log_changes({c: rf[c].mean() for c in CONDS}, ('L', 'R'))
    choice_path = out / 'estimator_choice.json'
    cand, lc = None, None
    if choice_path.exists():
        ch = json.loads(choice_path.read_text())
        b2 = ch['populations']['b2']
        if b2['chosen'] is not None and b2['candidates'][b2['chosen']]['passes_who']:
            cand, lc = b2['chosen'], ch['lc_edges']['chosen']
    rows = []
    for fam, d in fams.items():
        for tag, (measure, k) in d['members'].items():
            row = dict(family=fam, folder=d['folder'], variant=tag, measure=measure, exponent=k, ref=tag == d['ref'],
                       operating_point=d['operating_point'])
            if (d['folder'], tag) not in paths:
                row['error'] = 'not in the archive'
                rows.append(row)
                continue
            on, off = paths[(d['folder'], tag)]
            row.update(_b_reading(on, off))
            b = bv[(bv.folder == d['folder']) & (bv.tag == tag)]
            if len(b):
                b = b.iloc[0]
                row.update(bv_rmse=b['spike_counts_cosine__apl_effect_rmse'], bv_pattern=b['spike_counts_cosine__apl_effect_rmse_centered'],
                           bv_mean=b['spike_counts_cosine__apl_effect_mean'], bv_rho=b['spike_counts_cosine__apl_effect_spearman'],
                           rate_on=b['rate_on'], rate_off=b['rate_off'], realism_flag=b['realism_flag'])
            for rd in ('r1', 'r3raw', 'r3cal'):
                s = scan[(scan.folder == d['folder']) & (scan.variant == tag) & (scan.reading == rd)]
                if len(s):
                    s = s.iloc[0]
                    for c in ('own_pattern', 'own_err_alcohol', 'own_err_ester', 'sets_pattern', 'amp_pattern', 'bin_pattern',
                              'chg_L_alcohol', 'chg_L_ester', 'chg_L_nearest4', 'chg_R_alcohol', 'chg_R_ester', 'chg_R_nearest4'):
                        if c in s:
                            row[f'{rd}_{c}'] = s[c]
            rows.append(row)
    df = pd.DataFrame(rows)
    if cand is not None:
        rel = E88.flies_reliability('z')
        gap = {c: v['thresholded'] for c, v in flies_gap_target('z').items()}
        jobs = [(r.folder, r.variant, *paths[(r.folder, r.variant)], cand, lc, rel, gap)
                for r in df.itertuples() if (r.folder, r.variant) in paths]
        import multiprocessing
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            cr = pd.DataFrame(list(ex.map(_b_candidate_job, jobs)))
        cr = cr.rename(columns={c: f'cand_{c}' for c in cr.columns if c not in ('folder', 'variant')})
        df = df.merge(cr, on=['folder', 'variant'], how='left')
        fv = flies_values({lc: LC_EDGES[lc]} if cand == 'lc' else {})
        fl['cand'] = E88.log_changes(flies_cand(fv, cand, lc), ('L', 'R'))
    df.to_csv(out / 'weighting_archive.csv', index=False)
    v = dict(rule=dict(B_GAIN=B_GAIN, nearest='the member whose alcohol co-response change is nearest the flies\' (log scale)'),
             candidate=cand, lc_edges=lc, flies_changes={k: {kk: float(np.exp(x)) for kk, x in d.items()
                                                             if kk.split('_', 2)[2] in ('nearest4', 'ester', 'alcohol')}
                                                         for k, d in fl.items()}, families={})
    for fam, d in fams.items():
        g = df[(df.family == fam) & df.get('error', pd.Series(np.nan, index=df.index)).isna()]
        if g.empty or not g.ref.any():
            v['families'][fam] = dict(error='members missing')
            continue
        ref = g[g.ref].iloc[0]
        fv_ = dict(operating_point=d['operating_point'], given_apart=d['operating_point'] != 0.068, ref=d['ref'],
                   ref_cal_pattern=float(ref.cal_z8_pattern), ref_raw_pattern=float(ref.raw_z8_pattern), nearest={})
        for lab, col, flies in (('trial_means', 'r1_chg_L_alcohol', fl['r1']['chg_L_alcohol']),
                                ('deterministic_counts', 'r3raw_chg_L_alcohol', fl['r2']['chg_L_alcohol']),
                                ('candidate_L1', 'cand_L1_chg_L_alcohol', fl.get('cand', {}).get('chg_L_alcohol')),
                                ('candidate_L3', 'cand_L3_chg_L_alcohol', fl.get('cand', {}).get('chg_L_alcohol'))):
            if col not in g or flies is None or g[col].isna().all():
                continue
            dist = (g[col] - flies).abs()
            n = g.loc[dist.idxmin()]
            gain_cal = float(ref.cal_z8_pattern - n.cal_z8_pattern)
            gain_raw = float(ref.raw_z8_pattern - n.raw_z8_pattern)
            fv_['nearest'][lab] = dict(member=n.variant, measure=n.measure, exponent=float(n.exponent),
                                       change=float(np.exp(n[col])), flies=float(np.exp(flies)), ref_change=float(np.exp(ref[col])),
                                       cal_pattern=float(n.cal_z8_pattern), raw_pattern=float(n.raw_z8_pattern),
                                       gain_cal=gain_cal, gain_raw=gain_raw, realism_flag=n.get('realism_flag'),
                                       follows=bool(gain_cal >= B_GAIN and gain_raw >= B_GAIN))
        fv_['follows'] = bool(fv_['nearest'].get('trial_means', {}).get('follows', False))
        fv_['members'] = {r.variant: dict(measure=r.measure, exponent=float(r.exponent), cal_pattern=float(r.cal_z8_pattern),
                                          raw_pattern=float(r.raw_z8_pattern), cal_mean=float(r.cal_mean_effect),
                                          r1_alcohol_L_change=float(np.exp(r.r1_chg_L_alcohol)) if 'r1_chg_L_alcohol' in g and pd.notna(r.r1_chg_L_alcohol) else None,
                                          realism_flag=r.realism_flag if 'realism_flag' in g and isinstance(r.realism_flag, str) else '')
                          for r in g.itertuples()}
        v['families'][fam] = fv_
    v['follows_any_at_the_flies_operating_point'] = [f for f, x in v['families'].items()
                                                     if x.get('follows') and not x.get('given_apart')]
    _write_json(out / 'verdict_B.json', v)
    print('B: families that follow:', [f for f, x in v['families'].items() if x.get('follows')], flush=True)
    return v


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
        part_b(out, args.workers)
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
