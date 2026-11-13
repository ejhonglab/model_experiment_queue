#!/usr/bin/env python3
"""E87: the form of APL's inhibition, read from the flies' own responses.

Which inhibition, applied to the KC responses with APL output blocked (TNT), gives the control flies' distances; what
the same inhibition makes of the models' responses without APL; and a target for those responses that needs no APL
(experiment_priorities.md, E87; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-12 20:10).

A recording's responses are its KC x 16-odor trial means (the solvent left out). A response is a value above the
criterion (1 z; the registry's noise_threshold in baseline-noise units); values not above it are set to 0 on both sides
(the registry's thresholded effect). A form reduces KC k's response to odor o by c * A_o^q * B_k^p, A_o the odor's mean
(thresholded) response over the recording's KCs with APL output blocked relative to its mean over odors, B_k the KC's
breadth with APL output blocked relative to its mean over responding KCs (the number of odors it responds to, or its
summed response). Subtractive: the reduced value is kept if it stays above the criterion; divisive: the criterion is
raised by the same amount and surviving values are unchanged. c is set so that the mean response rate over the
recordings is the control flies'. The predicted control matrix is the mean over the recordings with APL output blocked
of the cosine distances of the reduced responses; the predicted effect is that minus the same recordings' own.

Parts (CLI `--part`):
- A: the forms on the flies' responses (no simulations): A1 the grid, the per-odor oracle and the form from KCs of
  similar tuning, as trial means, with the trial noise removed and binarized, in z and in baseline-noise units, with
  floors and permuted controls; A2 the form chosen on one half of the recordings and scored on the other.
- B: the models' APL-off responses through the same forms (no simulations); the two-by-two.
- C: the effect of a form within the recordings with APL output blocked as a target that needs no APL; the models and
  the archive against it.
- D1: the forms as a static inhibition in place of APL on the eight wirings (simulations); D2 only if C3 has a variant
  that stands out.
Outputs: results/E87_inhibition_form/ (comparison_table.md first); large outputs under outputs/E87_inhibition_form/.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E87_inhibition_form'

# the grid of A1 and the thresholds, fixed before any run (the entry's Design)
QS = (0.0, 0.5, 1.0, 1.5, 2.0)
PS = (0.0, 0.25, 0.5, 1.0, 2.0)
KINDS = ('subtractive', 'divisive')
BREADTHS = ('count', 'summed')
SIMILAR_POWERS = (1, 2, 4)
UNITS = ('z', 'noise')
READINGS = ('trial_means', 'noise_free', 'binary')
N_PERM = 20
N_RANDOM_SPLITS = 200
HOLD_GAIN = 0.005          # A2: a form beyond a uniform inhibition holds
BREADTH_GAIN = 0.005       # A2: the breadth dependence is needed in a unit
SIMILAR_GAIN = 0.003       # A2: the similar-tuning form counts
B_RESPONSES_NOT = 0.004    # B: the flies' form on the model's responses is not lower than the model's APL by this
B_FLIES_GAIN = 0.010       # B: ... while on the flies' responses it is lower than the model's APL by this
B_FORM_GAIN = 0.006        # B: the error is in the form of the model's APL (one wiring draw)
B_FORM_WIRINGS = 6
C_FLOOR_RATIO = 0.5        # C1: usable if its pattern floor is at most half the main target's ...
C_CEILING = 0.90           # ... and its ceiling on ρ is at least this
C_STANDOUT = 0.008         # C3: a variant stands out by this against the breadth-2 fit on the same wiring
C_MEAN_SHARE = 2.0 / 3.0
D_GAIN = 0.003             # D1: a static form does better / worse than the family's APL
D_WIRINGS = 6
D_RAW_SLACK = 0.005
D_UNCHANGED = 0.0482       # the eight wirings' pattern part under the readout against all recordings (E83)
D_RAW_UNCHANGED = 0.0729
D_PROXY_R = 0.8
D_PROXY_DIFF = 0.006

NEAR4 = (('2h', 'IaA'), ('2h', 'pa'), ('IaA', 'pa'), ('eb', 'ep'))


@dataclass(frozen=True)
class Form:
    """One form of inhibition: `kind` 'subtractive' or 'divisive', exponents `q` (on the odor's response over all KCs)
    and `p` (on the KC's breadth), `breadth` 'count' (odors responded to) or 'summed' (summed response). The two kinds
    remove the same responses (v - c I > crit and v > crit + c I are one condition), so they share c and differ only in
    the surviving values."""
    kind: str = 'subtractive'
    q: float = 0.0
    p: float = 0.0
    breadth: str = 'count'

    @property
    def name(self) -> str:
        b = '' if self.p == 0 else f'_{self.breadth}'
        return f'{self.kind[:3]}_q{self.q:g}_p{self.p:g}{b}'


UNIFORM = Form('subtractive', 0.0, 0.0, 'count')
UNIFORM_DIV = Form('divisive', 0.0, 0.0, 'count')

_CACHE: Dict[tuple, object] = {}


def _iu() -> Tuple[np.ndarray, np.ndarray]:
    return np.triu_indices(len(H.APL16_ODORS), 1)


def pair_index() -> pd.MultiIndex:
    """The 120 pairs in canonical (sorted-name) order, as E83's `_cs`."""
    i, j = _iu()
    o = H.APL16_ODORS
    return pd.MultiIndex.from_tuples([tuple(sorted((o[a], o[b]))) for a, b in zip(i, j)], names=['odor_a', 'odor_b'])


def _pairs(values: np.ndarray) -> pd.Series:
    """A 120-vector in triu order of APL16_ODORS as a canonical pair Series."""
    return pd.Series(np.asarray(values, dtype=float), index=pair_index()).sort_index()


def unit_scale(unit: str) -> Optional[pd.Series]:
    """Each KC's factor from Yang's z to `unit` (None for 'z'; `data_targets.kc_scale_ratios` for 'noise')."""
    if unit == 'z':
        return None
    if unit == 'noise':
        from kcq import data_targets as DT
        return DT.kc_scale_ratios()
    raise ValueError(unit)


def criterion(unit: str) -> float:
    """The trial-mean response criterion of `unit` ('z': 1; 'noise': the registry's noise_threshold)."""
    return float(H.KC_RESPONSE_Z) if unit == 'z' else float(H.data_target('noise_threshold'))


def trialmean_mats(unit: str) -> Dict[str, Dict[tuple, np.ndarray]]:
    """{'TNT3c' / 'TNT': {recording: KCs x 16 trial means (APL16 order)}} in `unit`, recordings sorted, KCs sorted by
    roi (`data_targets.kc_response_matrices`)."""
    key = ('tm', unit)
    if key not in _CACHE:
        from kcq import data_targets as DT
        mats = DT.kc_response_matrices(unit_scale(unit))
        _CACHE[key] = {g: {(str(r[0]), int(r[1]), str(r[2])): m.sort_index().to_numpy(dtype=float) for r, m in d.items()}
                       for g, d in mats.items()}
    return _CACHE[key]


def trial_arrays(unit: str) -> Dict[str, Dict[tuple, np.ndarray]]:
    """{'TNT3c' / 'TNT': {recording: (repeats, KCs, 16) single trials}} in `unit` (`H.yang_recording_trials`; KCs sorted
    by roi, as `trialmean_mats`)."""
    key = ('trials', unit)
    if key not in _CACHE:
        _CACHE[key] = {g: H.yang_recording_trials(g, unit) for g in ('TNT3c', 'TNT')}
    return _CACHE[key]


def _rate(mats: Sequence[np.ndarray], crit: float) -> float:
    return float(np.mean([float((m > crit).mean()) for m in mats]))


def control_rate(unit: str) -> float:
    """The control flies' mean response rate over recordings at `criterion(unit)` (0.068 by construction)."""
    return _rate(list(trialmean_mats(unit)['TNT3c'].values()), criterion(unit))


def single_trial_criterion(unit: str) -> float:
    """The single-trial criterion at which the control flies' single-trial response rate (mean over recordings of the
    share of repeat x KC x odor entries above it) is their trial-mean rate: geometric bisection on [0.1, 100]."""
    key = ('crit1', unit)
    if key not in _CACHE:
        target = control_rate(unit)
        arrs = list(trial_arrays(unit)['TNT3c'].values())
        lo, hi = 0.1, 100.0
        for _ in range(60):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if _rate(arrs, mid) > target else (lo, mid)
        _CACHE[key] = float(np.sqrt(lo * hi))
    return _CACHE[key]


def thresholded(m: np.ndarray, crit: float) -> np.ndarray:
    return np.where(m > crit, m, 0.0)


def strengths(m: np.ndarray, crit: float) -> np.ndarray:
    """A_o: each odor's mean thresholded response over the KCs (APL output blocked), relative to its mean over odors."""
    a = thresholded(m, crit).mean(axis=0)
    return a / a.mean() if a.mean() > 0 else np.ones_like(a)


def breadths(m: np.ndarray, crit: float, measure: str = 'count') -> np.ndarray:
    """B_k: each KC's breadth (the number of odors above the criterion, or its summed thresholded response), relative
    to its mean over responding KCs; 0 for KCs that respond to nothing."""
    if measure == 'count':
        b = (m > crit).sum(axis=1).astype(float)
    elif measure == 'summed':
        b = thresholded(m, crit).sum(axis=1)
    else:
        raise ValueError(measure)
    pos = b > 0
    return b / b[pos].mean() if pos.any() else b


def inhibition(m: np.ndarray, form: Form, crit: float, A: Optional[np.ndarray] = None,
               B: Optional[np.ndarray] = None) -> np.ndarray:
    """A_o^q * B_k^p (KCs x odors), the form's inhibition before its constant c. `A` / `B` replace the recording's own
    strengths / breadths (the permuted controls)."""
    A = strengths(m, crit) if A is None else A
    rowf = np.ones(m.shape[0])
    if form.p != 0:
        B = breadths(m, crit, form.breadth) if B is None else B
        rowf = np.where(B > 0, B, 0.0) ** form.p
    return np.outer(rowf, A ** form.q)


def apply_form(m: np.ndarray, c: float, inh: np.ndarray, kind: str, crit: float) -> np.ndarray:
    """The reduced responses (thresholded convention: values not above the criterion are 0). Subtractive: m - c inh
    where that stays above the criterion; divisive: m where m > crit + c inh. Any array shape (single trials too)."""
    keep = m - c * inh > crit
    if kind == 'subtractive':
        return np.where(keep, m - c * inh, 0.0)
    if kind == 'divisive':
        return np.where(keep, m, 0.0)
    raise ValueError(kind)


def fit_c(mats: Sequence[np.ndarray], inhs: Sequence[np.ndarray], crit: float, rate: float) -> float:
    """The constant c at which the mean over `mats` of the share of entries that survive (v - c inh > crit) is `rate`,
    exactly: an entry survives while c is below (v - crit) / inh, so c is read off the weighted, sorted thresholds."""
    ts, ws = [], []
    n_rec = len(mats)
    for m, inh in zip(mats, inhs):
        resp = m > crit
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.where(inh[resp] > 0, (m[resp] - crit) / inh[resp], np.inf)
        ts.append(t)
        ws.append(np.full(t.size, 1.0 / (m.size * n_rec)))
    t = np.concatenate(ts)
    w = np.concatenate(ws)
    always = float(w[np.isinf(t)].sum())
    if always >= rate:
        return 0.0 if always == rate else float('nan')
    fin = ~np.isinf(t)
    t, w = t[fin], w[fin]
    order = np.argsort(-t)
    t, cum = t[order], always + np.cumsum(w[order])
    k = int(np.searchsorted(cum, rate))   # the first k entries survive at c just below t[k - 1]
    if k <= 0:
        return float(t[0])
    if k >= len(t):
        return 0.0
    return float(0.5 * (t[k - 1] + t[k]))


def similar_inhibition(m: np.ndarray, crit: float, power: float) -> np.ndarray:
    """Each KC's inhibition by the KCs of similar tuning: for odor o, the weighted mean of the other responding KCs'
    thresholded responses to o, weights the positive part of the KC's correlation with each over the other 15 odors,
    to `power` (0 for a KC with no weight)."""
    r = thresholded(m, crit)
    rows = np.flatnonzero((r > 0).any(axis=1))
    out = np.zeros_like(r)
    if len(rows) < 2:
        return out
    x = r[rows]
    n_od = x.shape[1]
    for o in range(n_od):
        y = np.delete(x, o, axis=1)
        y = y - y.mean(axis=1, keepdims=True)
        sd = np.sqrt((y ** 2).sum(axis=1))
        ok = sd > 0
        z = np.zeros_like(y)
        z[ok] = y[ok] / sd[ok, None]
        cor = z @ z.T
        np.fill_diagonal(cor, 0.0)
        w = np.clip(cor, 0.0, None) ** power
        s = w.sum(axis=1)
        out[rows, o] = np.where(s > 0, (w @ x[:, o]) / np.where(s > 0, s, 1.0), 0.0)
    return out


def per_odor_oracle(m: np.ndarray, fracs: np.ndarray, kind: str, crit: float) -> np.ndarray:
    """Each odor's response fraction brought to `fracs` (the control flies' per-odor fractions): the k = round(frac N)
    largest values kept (and above the criterion); 'divisive' keeps their values, 'subtractive' subtracts the raised
    criterion's excess."""
    out = np.zeros_like(m)
    n = m.shape[0]
    for j in range(m.shape[1]):
        col = m[:, j]
        k = int(round(fracs[j] * n))
        if k <= 0:
            continue
        th = max(float(np.sort(col)[-k]), crit)
        keep = col >= th
        keep &= col > crit
        out[:, j] = np.where(keep, col if kind == 'divisive' else col - (th - crit), 0.0)
        if kind == 'subtractive':
            out[:, j] = np.where(out[:, j] > crit, out[:, j], np.where(keep, np.nextafter(crit, np.inf), 0.0))
    return out


def cos_dist(x: np.ndarray) -> np.ndarray:
    """1 - cosine similarity of the odors' vectors across units (a 120-vector in triu order; NaN where an odor is
    silent)."""
    g = x.T @ x
    n = np.sqrt(np.diag(g))
    with np.errstate(divide='ignore', invalid='ignore'):
        c = g / np.outer(n, n)
    return (1.0 - c)[_iu()]


def mean_dist(mats: Sequence[np.ndarray]) -> np.ndarray:
    """The mean over recordings of `cos_dist` (NaN-aware: a recording with a silent odor leaves its pairs out)."""
    return np.nanmean(np.stack([cos_dist(m) for m in mats]), axis=0)


def mean_dist_noise_free(trials: Sequence[np.ndarray]) -> np.ndarray:
    """The mean over recordings of the noise-free distances (`H.noise_free_pair_distances`, products of different
    repeats) of thresholded-convention single trials, in triu order."""
    i, j = _iu()
    o = H.APL16_ODORS
    keys = [(o[a], o[b]) for a, b in zip(i, j)]
    out = []
    for tr in trials:
        s = H.noise_free_pair_distances(tr, H.APL16_ODORS, values='graded')
        s.index = [tuple(x) for x in s.index]
        out.append(np.array([s.get(k, s.get(k[::-1], np.nan)) for k in keys], dtype=float))
    return np.nanmean(np.stack(out), axis=0)


def binarize(mats: Sequence[np.ndarray]) -> List[np.ndarray]:
    return [(m > 0).astype(float) for m in mats]


# ---------------------------------------------------------------------------
# the data and the predictions in one reading

class Reading:
    """The flies' recordings of one unit and one reading ('trial_means', 'noise_free' or 'binary'), restricted to the
    given recordings (default all), with the data's effect in the thresholded convention and `predict(form, ...)`."""

    def __init__(self, unit: str = 'z', reading: str = 'trial_means', ctrl: Optional[Sequence[tuple]] = None,
                 tnt: Optional[Sequence[tuple]] = None):
        self.unit, self.reading = unit, reading
        self.crit = criterion(unit)
        tm = trialmean_mats(unit)
        self.ctrl_keys = list(tm['TNT3c']) if ctrl is None else [tuple(k) for k in ctrl]
        self.tnt_keys = list(tm['TNT']) if tnt is None else [tuple(k) for k in tnt]
        self.tm_ctrl = [tm['TNT3c'][k] for k in self.ctrl_keys]
        self.tm_tnt = [tm['TNT'][k] for k in self.tnt_keys]
        self.rate = control_rate(unit)     # the control flies' rate over all recordings: the forms' operating point
        if reading == 'noise_free':
            tr = trial_arrays(unit)
            self.crit1 = single_trial_criterion(unit)
            self.tr_ctrl = [tr['TNT3c'][k] for k in self.ctrl_keys]
            self.tr_tnt = [tr['TNT'][k] for k in self.tnt_keys]
            self.rate1 = _rate(list(trial_arrays(unit)['TNT3c'].values()), self.crit1)
            d_on = mean_dist_noise_free([thresholded(x, self.crit1) for x in self.tr_ctrl])
            d_off = mean_dist_noise_free([thresholded(x, self.crit1) for x in self.tr_tnt])
        elif reading in ('trial_means', 'binary'):
            on = [thresholded(m, self.crit) for m in self.tm_ctrl]
            off = [thresholded(m, self.crit) for m in self.tm_tnt]
            if reading == 'binary':
                on, off = binarize(on), binarize(off)
            d_on, d_off = mean_dist(on), mean_dist(off)
        else:
            raise ValueError(reading)
        self.d_on, self.d_off = d_on, d_off
        self.data = d_on - d_off

    def inhibitions(self, form, A_perm=None, B_perm=None) -> List[np.ndarray]:
        """The inhibition (before c) of each TNT recording: a Form, ('similar', power), or with permuted strengths /
        breadths (lists of per-recording arrays)."""
        out = []
        for n, m in enumerate(self.tm_tnt):
            if isinstance(form, tuple) and form[0] == 'similar':
                key = ('sim', self.unit, self.tnt_keys[n], form[1])
                if key not in _CACHE:
                    _CACHE[key] = similar_inhibition(m, self.crit, form[1])
                out.append(_CACHE[key])
            else:
                out.append(inhibition(m, form, self.crit, None if A_perm is None else A_perm[n],
                                      None if B_perm is None else B_perm[n]))
        return out

    def predict(self, form, kind: Optional[str] = None, A_perm=None, B_perm=None, c: Optional[float] = None
                ) -> Tuple[np.ndarray, float]:
        """(the predicted per-pair effect in triu order, c): the form on the TNT recordings at the control flies' rate
        (or at the given c), its distances minus the same recordings' own."""
        kind = form.kind if isinstance(form, Form) else (kind or 'subtractive')
        inhs = self.inhibitions(form, A_perm, B_perm)
        if self.reading == 'noise_free':
            if c is None:
                c = fit_c(self.tr_tnt, [np.broadcast_to(i, x.shape) for i, x in zip(inhs, self.tr_tnt)], self.crit1, self.rate1)
            red = [apply_form(x, c, i[None], kind, self.crit1) for x, i in zip(self.tr_tnt, inhs)]
            d = mean_dist_noise_free(red)
        else:
            if c is None:
                c = fit_c(self.tm_tnt, inhs, self.crit, self.rate)
            red = [apply_form(m, c, i, kind, self.crit) for m, i in zip(self.tm_tnt, inhs)]
            if self.reading == 'binary':
                red = binarize(red)
            d = mean_dist(red)
        return d - self.d_off, c

    def predict_oracle(self, kind: str) -> np.ndarray:
        """The per-odor oracle's predicted effect (each odor's fraction set to the control flies' of these recordings;
        trial means and binary only)."""
        fr = np.mean([(m > self.crit).mean(axis=0) for m in self.tm_ctrl], axis=0)
        red = [per_odor_oracle(m, fr, kind, self.crit) for m in self.tm_tnt]
        if self.reading == 'binary':
            red = binarize(red)
        return mean_dist(red) - self.d_off

    def reduced(self, form, c: float) -> List[np.ndarray]:
        """The reduced trial-mean matrices of the TNT recordings (for the realism reading)."""
        kind = form.kind if isinstance(form, Form) else 'subtractive'
        return [apply_form(m, c, i, kind, self.crit) for m, i in zip(self.tm_tnt, self.inhibitions(form))]


def scores(pred: np.ndarray, data: np.ndarray) -> dict:
    """loss, mean, pattern part, ρ, the share of the data's variance explained by a line on the prediction (r²) and its
    slope, the ester and alcohol groups' mean centred errors, the four nearest ester pairs and the six near pairs, and
    the fruity x small-ester block of the predicted effect (all in canonical pair order)."""
    from experiments import E77_input_side_common_error as E77
    from experiments import E78_claw_knee as E78
    ps, ds = _pairs(pred), _pairs(data)
    ok = ps.notna() & ds.notna()
    e = (ps - ds)[ok]
    ec = (ps[ok] - ps[ok].mean()) - (ds[ok] - ds[ok].mean())
    masks = E77.kind_masks(ps.index)
    six = np.array([tuple(x) in {tuple(sorted(q)) for q in E78.SIX} for x in ps.index])
    near = masks['nearest']
    r = float(np.corrcoef(ps[ok], ds[ok])[0, 1])
    slope = float(np.polyfit(ps[ok], ds[ok], 1)[0])
    blk = H.block_means(ps.dropna())
    return dict(loss=float(np.sqrt((e ** 2).mean())), mean=float(ps.mean()), pattern=float(np.sqrt((ec ** 2).mean())),
                rho=float(H.spearman(ps[ok], ds[ok])), r2=r * r, slope=slope,
                err_ester=float(ec[masks['ester'][ok.to_numpy()]].mean()), err_alcohol=float(ec[masks['alcohol'][ok.to_numpy()]].mean()),
                near4=float(ps[near].mean()), six=float(ps[six].mean()),
                block=float(blk.get('fruity-small_ester', np.nan)), n_pairs=int(ok.sum()))


def realism_of_form(rd: Reading, red: Sequence[np.ndarray]) -> dict:
    """What a form predicts for the control flies beside what they have: the response rate, ρ of the per-odor response
    fractions, the distance of the odors-per-KC distribution from the control flies' (total variation) and the
    responders' mean amplitude over the control flies'."""
    from kcq import data_targets as DT
    fr_p = np.mean([(m > 0).mean(axis=0) for m in red], axis=0)
    fr_c = np.mean([(m > rd.crit).mean(axis=0) for m in rd.tm_ctrl], axis=0)
    opk_p = np.mean([DT.opk_bin_fractions((m > 0).sum(axis=1)) for m in red], axis=0)
    opk_c = np.mean([DT.opk_bin_fractions((m > rd.crit).sum(axis=1)) for m in rd.tm_ctrl], axis=0)
    amp_p = np.mean([m[m > 0].mean() for m in red if (m > 0).any()])
    amp_c = np.mean([m[m > rd.crit].mean() for m in rd.tm_ctrl])
    return dict(rate_pred=float(np.mean([(m > 0).mean() for m in red])), rate_ctrl=float(fr_c.mean()),
                rate_rho=float(H.spearman(pd.Series(fr_p), pd.Series(fr_c))),
                opk_tv=float(0.5 * np.abs(opk_p - opk_c).sum()), amp_ratio=float(amp_p / amp_c))


def grid() -> List[Form]:
    """A1's forms: q x p x kind x breadth (the breadth measure only where p > 0)."""
    out = []
    for kind in KINDS:
        for q in QS:
            for p in PS:
                for b in (BREADTHS if p > 0 else ('count',)):
                    out.append(Form(kind, q, p, b))
    return out


def _perm_strengths(rd: Reading, seed: int) -> List[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [rng.permutation(strengths(m, rd.crit)) for m in rd.tm_tnt]


def _perm_breadths(rd: Reading, form: Form, seed: int) -> List[np.ndarray]:
    """Each TNT recording's breadths permuted over its responding KCs (non-responders keep 0)."""
    rng = np.random.default_rng(seed)
    out = []
    for m in rd.tm_tnt:
        b = breadths(m, rd.crit, form.breadth)
        pos = np.flatnonzero(b > 0)
        b2 = b.copy()
        b2[pos] = b[rng.permutation(pos)]
        out.append(b2)
    return out


def _choose(rd: Reading, forms: Sequence) -> Tuple[object, dict]:
    """The form of `forms` with the lowest pattern part on `rd` (ties: the first), with its scores."""
    best, best_s = None, None
    for f in forms:
        pred, c = rd.predict(f)
        s = scores(pred, rd.data)
        if best_s is None or s['pattern'] < best_s['pattern']:
            best, best_s = f, dict(s, c=c)
    return best, best_s


def _form_cols(f) -> dict:
    if isinstance(f, Form):
        return dict(form=f.name, family='global', kind=f.kind, q=f.q, p=f.p, breadth=f.breadth if f.p else '')
    if isinstance(f, tuple) and f[0] == 'similar':
        return dict(form=f'similar_pow{f[1]:g}', family='similar', kind='subtractive', q=np.nan, p=np.nan, breadth='')
    raise ValueError(f)


def part_a1(out: Path, units: Sequence[str] = UNITS, readings: Sequence[str] = READINGS) -> pd.DataFrame:
    """A1: every form of the grid, the per-odor oracle (both kinds) and the similar-tuning form (three powers) on all
    recordings, per unit and reading: forms.csv (scores, c, the realism reading for trial means) and forms_pairs.csv
    (the predicted effect per pair, and the data's)."""
    rows, prs = [], []
    for unit in units:
        for reading in readings:
            rd = Reading(unit, reading)
            prs.append(pd.DataFrame(dict(unit=unit, reading=reading, form='data', effect=_pairs(rd.data))).reset_index())
            forms = grid() + [('similar', pw) for pw in SIMILAR_POWERS]
            for f in forms:
                pred, c = rd.predict(f)
                row = dict(unit=unit, reading=reading, **_form_cols(f), c=c, **scores(pred, rd.data))
                if reading == 'trial_means':
                    row.update(realism_of_form(rd, rd.reduced(f, c)))
                rows.append(row)
                prs.append(pd.DataFrame(dict(unit=unit, reading=reading, form=_form_cols(f)['form'], effect=_pairs(pred))).reset_index())
            if reading != 'noise_free':
                for kind in KINDS:
                    pred = rd.predict_oracle(kind)
                    rows.append(dict(unit=unit, reading=reading, form=f'per_odor_oracle_{kind[:3]}', family='oracle', kind=kind,
                                     q=np.nan, p=np.nan, breadth='', c=np.nan, **scores(pred, rd.data)))
                    prs.append(pd.DataFrame(dict(unit=unit, reading=reading, form=f'per_odor_oracle_{kind[:3]}', effect=_pairs(pred))).reset_index())
            rows.append(dict(unit=unit, reading=reading, form='none', family='none', kind='', q=np.nan, p=np.nan, breadth='',
                             c=0.0, **scores(np.zeros_like(rd.data), rd.data)))
            print(f'A1 {unit} {reading}: {len(forms)} forms', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'forms.csv', index=False)
    pd.concat(prs, ignore_index=True).to_csv(out / 'forms_pairs.csv', index=False)
    return df


def part_a1_controls(out: Path, forms_df: pd.DataFrame, units: Sequence[str] = UNITS) -> pd.DataFrame:
    """A1's controls on all recordings (trial means and the trial noise removed): for the best global form of each unit
    and reading by pattern part and for q 1 / p 0.25, the odors' strengths permuted over the odors and the breadths
    permuted over the responding KCs (N_PERM seeds each); and the floors: the recordings of both genotypes resampled
    (300 draws), the data's own floor and the spread of the held form's prediction. controls.csv, floors.csv."""
    rows, frows = [], []
    for unit in units:
        for reading in ('trial_means', 'noise_free'):
            rd = Reading(unit, reading)
            sub = forms_df[(forms_df.unit == unit) & (forms_df.reading == reading) & (forms_df.family == 'global')]
            best = sub.sort_values('pattern').iloc[0]
            cand = {'best': Form(best.kind, float(best.q), float(best.p), best.breadth or 'count'),
                    'q1_p0.25': Form('subtractive', 1.0, 0.25, 'count')}
            for lab, f in cand.items():
                base = scores(rd.predict(f)[0], rd.data)
                for ctl in ('strengths', 'breadths'):
                    if ctl == 'breadths' and f.p == 0:
                        continue
                    for s in range(N_PERM):
                        if ctl == 'strengths':
                            pred, c = rd.predict(f, A_perm=_perm_strengths(rd, s))
                        else:
                            pred, c = rd.predict(f, B_perm=_perm_breadths(rd, f, s))
                        rows.append(dict(unit=unit, reading=reading, which=lab, form=f.name, control=ctl, seed=s, c=c,
                                         base_pattern=base['pattern'], base_loss=base['loss'], **scores(pred, rd.data)))
            if reading == 'trial_means':
                rng = np.random.default_rng(0)
                fu = UNIFORM
                fb = cand['best']
                cu = rd.predict(fu)[1]
                cb = rd.predict(fb)[1]
                pu0, pb0 = rd.predict(fu)[0], rd.predict(fb)[0]
                for k in range(300):
                    ic = rng.integers(0, len(rd.ctrl_keys), len(rd.ctrl_keys))
                    it = rng.integers(0, len(rd.tnt_keys), len(rd.tnt_keys))
                    r2 = Reading(unit, reading, [rd.ctrl_keys[i] for i in ic], [rd.tnt_keys[i] for i in it])
                    e = r2.data - rd.data
                    ok = ~np.isnan(e)
                    ec = e[ok] - e[ok].mean()
                    pu, pb = r2.predict(fu, c=cu)[0], r2.predict(fb, c=cb)[0]
                    frows.append(dict(unit=unit, reading=reading, draw=k,
                                      data_loss=float(np.sqrt(np.mean(e[ok] ** 2))), data_pattern=float(np.sqrt(np.mean(ec ** 2))),
                                      uniform_spread_pattern=float(np.nanstd(pu - pu0)), best_spread_pattern=float(np.nanstd(pb - pb0))))
            print(f'A1 controls {unit} {reading}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'controls.csv', index=False)
    pd.DataFrame(frows).to_csv(out / 'floors.csv', index=False)
    return df


def part_a1_recordings(out: Path, units: Sequence[str] = UNITS) -> pd.DataFrame:
    """A1's subsets (trial means): the six typical TNT recordings, the two strongest, and each TNT recording alone
    predicting all control recordings: the best global form by pattern part and the uniform form's scores.
    recordings.csv."""
    from experiments import E77_input_side_common_error as E77
    rows = []
    for unit in units:
        tm = trialmean_mats(unit)
        allt = list(tm['TNT'])
        strong = [tuple(k) for k in E77.STRONG_TNT]
        subsets = {'six_typical': [k for k in allt if k not in strong], 'two_strongest': strong}
        subsets.update({f'alone_{k[0]}_{k[1]}_{k[2]}': [k] for k in allt})
        for lab, keys in subsets.items():
            rd = Reading(unit, 'trial_means', tnt=keys)
            f, s = _choose(rd, grid())
            u = scores(rd.predict(UNIFORM)[0], rd.data)
            rows.append(dict(unit=unit, subset=lab, n_tnt=len(keys), **{f'best_{k}': v for k, v in _form_cols(f).items()},
                             **{f'best_{k}': v for k, v in s.items()}, uniform_pattern=u['pattern'], uniform_loss=u['loss'],
                             uniform_mean=u['mean'], data_mean=float(np.nanmean(rd.data))))
        print(f'A1 recordings {unit}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'recordings.csv', index=False)
    return df


def _halves(seed: int) -> Dict[str, Dict[str, list]]:
    from experiments import E77_input_side_common_error as E77
    h = E77.recording_halves(seed)
    return {k: {g: [(str(r[0]), int(r[1]), str(r[2])) for r in v[g]] for g in ('ctrl', 'tnt')} for k, v in h.items()}


def held_out_one(unit: str, reading: str, split: int, direction: str, forms: Sequence[Form],
                 controls: bool = True) -> dict:
    """One held-out reading: the global form (of `forms`) and the similar-tuning power chosen by the pattern part on
    half X's recordings (its TNT recordings predicting its control recordings), each scored with the uniform form, the
    best p = 0 form and their permuted controls on half Y."""
    h = _halves(split)
    X, Y = direction[0], direction[-1]
    rx = Reading(unit, reading, h[X]['ctrl'], h[X]['tnt'])
    ry = Reading(unit, reading, h[Y]['ctrl'], h[Y]['tnt'])
    f, sx = _choose(rx, forms)
    f0, _ = _choose(rx, [g for g in forms if g.p == 0])
    sim, _ = _choose(rx, [('similar', pw) for pw in SIMILAR_POWERS])
    row = dict(unit=unit, reading=reading, split=split, direction=direction, **{f'chosen_{k}': v for k, v in _form_cols(f).items()},
               chosen_pattern_X=sx['pattern'], p0_form=f0.name, similar_form=_form_cols(sim)['form'])
    for lab, g in (('chosen', f), ('uniform', UNIFORM), ('p0', f0), ('similar', sim)):
        s = scores(ry.predict(g)[0], ry.data)
        row.update({f'{lab}_{k}': s[k] for k in ('pattern', 'loss', 'mean', 'rho')})
    if controls:
        perm = [scores(ry.predict(f, A_perm=_perm_strengths(ry, s))[0], ry.data) for s in range(N_PERM)]
        row['chosen_strengths_perm_pattern'] = float(np.mean([s['pattern'] for s in perm]))
        if f.p > 0:
            pb = [scores(ry.predict(f, B_perm=_perm_breadths(ry, f, s))[0], ry.data) for s in range(N_PERM)]
            row['chosen_breadths_perm_loss'] = float(np.mean([s['loss'] for s in pb]))
            row['chosen_breadths_perm_pattern'] = float(np.mean([s['pattern'] for s in pb]))
    return row


def part_a2(out: Path, units: Sequence[str] = UNITS, n_random: int = N_RANDOM_SPLITS) -> pd.DataFrame:
    """A2: the four held-out readings (E77's `recording_halves()` and `recording_halves(1)`, both directions) per unit,
    as trial means and with the trial noise removed (held_out.csv); then n_random further splits (seeds 2 onward) in z,
    trial means, chosen against uniform only (held_out_random.csv)."""
    rows = []
    forms = grid()
    for unit in units:
        for reading in ('trial_means', 'noise_free'):
            for split in (0, 1):
                for direction in ('A->B', 'B->A'):
                    rows.append(held_out_one(unit, reading, split, direction, forms))
                    print(f'A2 {unit} {reading} split {split} {direction}: chosen {rows[-1]["chosen_form"]} '
                          f'{rows[-1]["chosen_pattern"]:.4f} uniform {rows[-1]["uniform_pattern"]:.4f}', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'held_out.csv', index=False)
    rrows = []
    for seed in range(2, 2 + n_random):
        for direction in ('A->B', 'B->A'):
            r = held_out_one('z', 'trial_means', seed, direction, forms, controls=False)
            rrows.append({k: r[k] for k in ('split', 'direction', 'chosen_form', 'chosen_pattern', 'uniform_pattern',
                                            'p0_pattern', 'similar_pattern', 'chosen_loss', 'uniform_loss')})
        if seed % 20 == 0:
            pd.DataFrame(rrows).to_csv(out / 'held_out_random.csv', index=False)
    pd.DataFrame(rrows).to_csv(out / 'held_out_random.csv', index=False)
    return df


def verdict_a(out: Path) -> dict:
    """The thresholds of A2, fixed before any run (HOLD_GAIN, BREADTH_GAIN, SIMILAR_GAIN), read from held_out.csv."""
    ho = pd.read_csv(out / 'held_out.csv')
    v = {}
    zz = ho[ho.unit == 'z']
    holds = {}
    for reading, g in zz.groupby('reading'):
        gain = g.uniform_pattern - g.chosen_pattern
        perm_gain = g.uniform_pattern - g.chosen_strengths_perm_pattern
        holds[reading] = dict(gains=[float(x) for x in gain], n_ok=int((gain >= HOLD_GAIN).sum()), mean_gain=float(gain.mean()),
                              perm_mean_gain=float(perm_gain.mean()),
                              holds=bool((gain >= HOLD_GAIN).sum() >= 3 and gain.mean() >= HOLD_GAIN
                                         and perm_gain.mean() < 0.5 * gain.mean()))
    v['form_holds'] = dict(by_reading=holds, holds=all(x['holds'] for x in holds.values()))
    br = {}
    for unit, g in ho[ho.reading == 'trial_means'].groupby('unit'):
        gain = g.p0_loss - g.chosen_loss
        has_p = g.chosen_p > 0
        perm = (g.p0_loss - g.chosen_breadths_perm_loss) if 'chosen_breadths_perm_loss' in g else pd.Series(np.nan, index=g.index)
        br[unit] = dict(gains=[float(x) for x in gain], chosen_p=[float(x) for x in g.chosen_p], mean_gain=float(gain.mean()),
                        perm_mean_gain=float(perm[has_p].mean()) if has_p.any() else float('nan'),
                        needed=bool(has_p.all() and gain.mean() >= BREADTH_GAIN and (perm[has_p].mean() < 0.5 * gain.mean())))
    v['breadth_needed'] = br
    sim_gain = zz.chosen_pattern - zz.similar_pattern
    v['similar_counts'] = dict(mean_gain=float(sim_gain.mean()), counts=bool(sim_gain.mean() >= SIMILAR_GAIN))
    (out / 'verdict_A.json').write_text(__import__('json').dumps(v, indent=1) + '\n')
    return v


def part_a(out: Path) -> dict:
    """A1 and A2 with their controls, then A's verdict."""
    import json
    (out / 'provenance_A.json').write_text(json.dumps(dict(script='experiments/E87_inhibition_form.py', part='A',
        **H.code_commits(), criteria={u: criterion(u) for u in UNITS}, single_trial_criteria={u: single_trial_criterion(u) for u in UNITS},
        control_rate={u: control_rate(u) for u in UNITS}), indent=1, default=str) + '\n')
    f = part_a1(out)
    part_a1_controls(out, f)
    part_a1_recordings(out)
    part_a2(out)
    return verdict_a(out)


def part_b(out: Path) -> dict:
    """The models' APL-off responses through the forms: model_forms.csv, two_by_two.csv."""
    raise NotImplementedError


def part_c(out: Path) -> dict:
    """The target without APL (apl_free_target.csv), the models against it, the archive scan (archive_scan.csv)."""
    raise NotImplementedError


def part_d1(out: Path) -> dict:
    """The forms as a static inhibition in place of APL on the eight wirings (runs_D.csv, readings_D.csv)."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'B', 'C', 'D1', 'D2', 'verdict'), required=True)
    args = ap.parse_args()
    out = H.results_dir(EXPERIMENT)
    {'A': part_a, 'B': part_b, 'C': part_c, 'D1': part_d1}[args.part](out)


if __name__ == '__main__':
    main()
