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


def rec_dist(mats: Sequence[np.ndarray]) -> np.ndarray:
    """`cos_dist` of each recording (recordings x 120, triu order)."""
    return np.stack([cos_dist(m) for m in mats]) if len(mats) else np.full((0, len(_iu()[0])), np.nan)


def rec_dist_noise_free(trials: Sequence[np.ndarray]) -> np.ndarray:
    """The noise-free distances (`H.noise_free_pair_distances`, products of different repeats) of each recording's
    thresholded-convention single trials (recordings x 120, triu order)."""
    i, j = _iu()
    o = H.APL16_ODORS
    keys = [(o[a], o[b]) for a, b in zip(i, j)]
    out = []
    for tr in trials:
        s = H.noise_free_pair_distances(tr, H.APL16_ODORS, values='graded')
        s.index = [tuple(x) for x in s.index]
        out.append(np.array([s.get(k, s.get(k[::-1], np.nan)) for k in keys], dtype=float))
    return np.stack(out) if out else np.full((0, len(keys)), np.nan)


def _nanmean_rows(r: np.ndarray) -> np.ndarray:
    """The mean over recordings (rows), NaN-aware: a recording with a silent odor leaves its pairs out."""
    if r.shape[0] == 0:
        return np.full(r.shape[1], np.nan)
    with np.errstate(invalid='ignore'):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return np.nanmean(r, axis=0)


def mean_dist(mats: Sequence[np.ndarray]) -> np.ndarray:
    """The mean over recordings of `cos_dist`."""
    return _nanmean_rows(rec_dist(mats))


def mean_dist_noise_free(trials: Sequence[np.ndarray]) -> np.ndarray:
    """The mean over recordings of `rec_dist_noise_free`."""
    return _nanmean_rows(rec_dist_noise_free(trials))


def binarize(mats: Sequence[np.ndarray]) -> List[np.ndarray]:
    return [(m > 0).astype(float) for m in mats]


# ---------------------------------------------------------------------------
# the data and the predictions in one reading

class Reading:
    """The flies' recordings of one unit and one reading ('trial_means', 'noise_free' or 'binary'), restricted to the
    given recordings (default all; an empty list of control recordings leaves the data's side out), with the data's
    effect in the thresholded convention, the per-recording distances (`rec_on`, `rec_off`: recordings x 120) and
    `predict(form, ...)` / `predict_rows(form, c)`."""

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
            self.rec_on = rec_dist_noise_free([thresholded(x, self.crit1) for x in self.tr_ctrl])
            self.rec_off = rec_dist_noise_free([thresholded(x, self.crit1) for x in self.tr_tnt])
        elif reading in ('trial_means', 'binary'):
            on = [thresholded(m, self.crit) for m in self.tm_ctrl]
            off = [thresholded(m, self.crit) for m in self.tm_tnt]
            if reading == 'binary':
                on, off = binarize(on), binarize(off)
            self.rec_on, self.rec_off = rec_dist(on), rec_dist(off)
        else:
            raise ValueError(reading)
        self.d_on, self.d_off = _nanmean_rows(self.rec_on), _nanmean_rows(self.rec_off)
        self.data = self.d_on - self.d_off

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

    def predict_rows(self, form, kind: Optional[str] = None, A_perm=None, B_perm=None, c: Optional[float] = None
                     ) -> Tuple[np.ndarray, float]:
        """(the reduced TNT recordings' distances, recordings x 120, c): the form on the TNT recordings at the control
        flies' rate (or at the given c)."""
        kind = form.kind if isinstance(form, Form) else (kind or 'subtractive')
        inhs = self.inhibitions(form, A_perm, B_perm)
        if self.reading == 'noise_free':
            if c is None:
                c = fit_c(self.tr_tnt, [np.broadcast_to(i, x.shape) for i, x in zip(inhs, self.tr_tnt)], self.crit1, self.rate1)
            red = [apply_form(x, c, i[None], kind, self.crit1) for x, i in zip(self.tr_tnt, inhs)]
            return rec_dist_noise_free(red), c
        if c is None:
            c = fit_c(self.tm_tnt, inhs, self.crit, self.rate)
        red = [apply_form(m, c, i, kind, self.crit) for m, i in zip(self.tm_tnt, inhs)]
        if self.reading == 'binary':
            red = binarize(red)
        return rec_dist(red), c

    def predict(self, form, kind: Optional[str] = None, A_perm=None, B_perm=None, c: Optional[float] = None
                ) -> Tuple[np.ndarray, float]:
        """(the predicted per-pair effect in triu order, c): the form on the TNT recordings at the control flies' rate
        (or at the given c), its distances minus the same recordings' own."""
        rows, c = self.predict_rows(form, kind, A_perm, B_perm, c)
        return _nanmean_rows(rows) - self.d_off, c

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
    if float(np.std(ps[ok])) > 0:
        r = float(np.corrcoef(ps[ok], ds[ok])[0, 1])
        slope = float(np.polyfit(ps[ok], ds[ok], 1)[0])
    else:   # a constant prediction (no inhibition): no line
        r, slope = float('nan'), float('nan')
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


# ---------------------------------------------------------------------------
# part B: the models' responses through the same forms

RAW_CRIT = 1.5      # a model response is two or more spikes (counts are integers: above 1.5)


def _odor_matrix(counts: pd.DataFrame) -> np.ndarray:
    """KCs x the 16 APL odors (APL16 order) of a model's spike counts (columns abbreviated)."""
    c = H.abbrev_columns(counts)
    return c[list(H.APL16_ODORS)].to_numpy(dtype=float)


class ModelResponses:
    """One model's APL-on and APL-off responses read as the flies' are in part A: 'raw' (spike counts, a response two
    or more spikes) or 'cal' (E65's calibrated saturating readout n / (n + K), K set on the run's APL-on counts, the
    criterion the readout of 1.5 spikes); values not above the criterion set to 0."""

    def __init__(self, o: Dict[str, tuple], readout: str, K: Optional[float] = None, label: str = ''):
        from experiments import E83_cross_wiring_fit as E
        self.label, self.readout = label, readout
        on, off = _odor_matrix(o['on'][1]), _odor_matrix(o['off'][1])
        if readout == 'raw':
            self.crit = RAW_CRIT
        elif readout == 'cal':
            K = E.readout_K(o) if K is None else K
            self.K = K
            on, off = on / (on + K), off / (off + K)
            self.crit = RAW_CRIT / (RAW_CRIT + K)
        else:
            raise ValueError(readout)
        self.on, self.off = on, off
        self.rate_on = float((on > self.crit).mean())
        self.d_on = cos_dist(thresholded(on, self.crit))
        self.d_off = cos_dist(thresholded(off, self.crit))

    def form(self, form, c: Optional[float] = None) -> Tuple[np.ndarray, float]:
        """(the APL-off distances after the form, c): c set so that the APL-off response rate after the form is the
        model's own APL-on rate."""
        kind = form.kind if isinstance(form, Form) else 'subtractive'
        if isinstance(form, tuple) and form[0] == 'similar':
            inh = similar_inhibition(self.off, self.crit, form[1])
        else:
            inh = inhibition(self.off, form, self.crit)
        if c is None:
            c = fit_c([self.off], [inh], self.crit, self.rate_on)
        return cos_dist(apply_form(self.off, c, inh, kind, self.crit)), c


def model_families() -> Dict[str, List[Tuple[str, str, str]]]:
    """{'b2': [(wiring, folder, tag)] the breadth-2 APL-effect fits of the eight wirings (E72 part B, E83's members),
    'fc': the firing-ceiling family (E69's fit at the literature setting on the real wiring, E72 part G's on P1, P2,
    M1)}."""
    from experiments import E72_wiring_population as E72
    from experiments import E83_cross_wiring_fit as E
    b2 = [(w, f, t) for w, (f, t) in E.members().items()]
    fc = E72.new_fit_members(prefix_tag=E72.G_SUFFIX, real_folder=E72.G_REAL[0], real_tag=E72.G_REAL[1], draws=E72.G_DRAWS)
    return {'b2': b2, 'fc': [tuple(x) for x in fc]}


def _load(folder: str, tag: str) -> Dict[str, tuple]:
    from experiments import E72_wiring_population as E72
    o = E72.load_fit(folder, tag)
    if o is None:
        raise RuntimeError(f'no saved outputs for {folder}:{tag}')
    return o


def chosen_form(out: Path) -> Form:
    """The form A2 chose: of the global forms chosen on half X in the four held-out readings in z with the trial noise
    removed (the reading the deterministic models are set against), the one chosen most often, ties broken by the lower
    mean pattern part on the unseen halves."""
    ho = pd.read_csv(out / 'held_out.csv')
    g = ho[(ho.unit == 'z') & (ho.reading == 'noise_free')]
    agg = g.groupby('chosen_form').agg(n=('chosen_form', 'size'), pat=('chosen_pattern', 'mean')).reset_index()
    best = agg.sort_values(['n', 'pat'], ascending=[False, True]).iloc[0].chosen_form
    r = g[g.chosen_form == best].iloc[0]
    return Form(r.chosen_kind, float(r.chosen_q), float(r.chosen_p), r.chosen_breadth if isinstance(r.chosen_breadth, str) and r.chosen_breadth else 'count')


def _pop(ds: Sequence[np.ndarray]) -> np.ndarray:
    return np.nanmean(np.stack(ds), axis=0)


def part_b(out: Path) -> dict:
    """The models through the forms (no simulations): each family's wirings and population, raw and under the
    calibrated readout, with the model's own APL, the uniform form, the form A2 chose, the family's own exponent (q 0,
    p 2) and the best of the grid for that model (by pattern part against the data), scored against the data's effect
    with the trial noise removed and as trial means (z); the form nearest the model's own APL (the grid form whose APL-off
    result is closest, rms over pairs, to the model's APL-on distances) applied to the flies; the two-by-two and B's
    verdict. model_forms.csv, two_by_two.csv, verdict_B.json."""
    import json
    f_ch = chosen_form(out)
    forms = grid()
    data = {r: Reading('z', r) for r in ('noise_free', 'trial_means')}
    rows, tb = [], []
    vb = {}
    for fam, members in model_families().items():
        for readout in ('raw', 'cal'):
            mrs = [ModelResponses(_load(f, t), readout, label=w) for w, f, t in members]
            res = {}   # (name, unit) -> effect vector
            per = {'own': [], 'uniform': [], 'chosen': [], 'q0_p2': []}
            for mr in mrs:
                per['own'].append((mr.d_on, mr.d_off))
                for lab, f in (('uniform', UNIFORM), ('chosen', f_ch), ('q0_p2', Form('subtractive', 0.0, 2.0, 'count'))):
                    per[lab].append((mr.form(f)[0], mr.d_off))
            for lab, pairs in per.items():
                for n, (d1, d0) in enumerate(pairs):
                    res[(lab, mrs[n].label)] = d1 - d0
                res[(lab, 'population')] = _pop([a for a, _ in pairs]) - _pop([b for _, b in pairs])
            # the best of the grid for the population, and the form nearest the model's own APL
            pop_on = _pop([mr.d_on for mr in mrs])
            pop_off = _pop([mr.d_off for mr in mrs])
            best, best_pat, near, near_rms = None, np.inf, None, np.inf
            grid_d = {}
            for f in forms:
                d = _pop([mr.form(f)[0] for mr in mrs])
                grid_d[f] = d
                pat = scores(d - pop_off, data['noise_free'].data)['pattern']
                if pat < best_pat:
                    best, best_pat = f, pat
                rms = float(np.sqrt(np.nanmean((d - pop_on) ** 2)))
                if rms < near_rms:
                    near, near_rms = f, rms
            res[('best_grid', 'population')] = grid_d[best] - pop_off
            for (lab, unit), eff in res.items():
                for rname, rd in data.items():
                    rows.append(dict(family=fam, readout=readout, inhibition=lab, unit=unit, data_reading=rname,
                                     form={'own': "the model's own APL", 'uniform': UNIFORM.name, 'chosen': f_ch.name,
                                           'q0_p2': 'sub_q0_p2_count', 'best_grid': best.name}[lab], **scores(eff, rd.data)))
            # the flies' side of the two-by-two
            flies = {}
            for rname, rd in data.items():
                flies[('uniform', rname)] = scores(rd.predict(UNIFORM)[0], rd.data)
                flies[('chosen', rname)] = scores(rd.predict(f_ch)[0], rd.data)
                flies[('nearest_to_model', rname)] = scores(rd.predict(near)[0], rd.data)
            for rname, rd in data.items():
                for who, lab, sc in (
                        ('flies', 'uniform', flies[('uniform', rname)]), ('flies', 'chosen', flies[('chosen', rname)]),
                        ('flies', 'nearest_to_model', flies[('nearest_to_model', rname)]),
                        ('model', 'uniform', scores(res[('uniform', 'population')], rd.data)),
                        ('model', 'chosen', scores(res[('chosen', 'population')], rd.data)),
                        ('model', 'own', scores(res[('own', 'population')], rd.data)),
                        ('model', 'best_grid', scores(res[('best_grid', 'population')], rd.data))):
                    tb.append(dict(family=fam, readout=readout, data_reading=rname, responses=who, inhibition=lab,
                                   form={'uniform': UNIFORM.name, 'chosen': f_ch.name, 'nearest_to_model': near.name,
                                         'own': "the model's own APL", 'best_grid': best.name}[lab],
                                   nearest_rms=near_rms if lab == 'nearest_to_model' else np.nan, **sc))
            if fam == 'b2':
                nf = data['noise_free'].data
                own = scores(res[('own', 'population')], nf)['pattern']
                mch = scores(res[('chosen', 'population')], nf)['pattern']
                fch = flies[('chosen', 'noise_free')]['pattern']
                lower = sum(scores(res[('chosen', mr.label)], nf)['pattern'] <= scores(res[('own', mr.label)], nf)['pattern'] - B_FORM_GAIN
                            for mr in mrs)
                vb[readout] = dict(own=own, model_chosen=mch, flies_chosen=fch, n_wirings_chosen_lower=int(lower),
                                   in_responses=bool(own - mch < B_RESPONSES_NOT and fch <= own - B_FLIES_GAIN),
                                   in_form=bool(mch <= own - B_FORM_GAIN and lower >= B_FORM_WIRINGS),
                                   nearest_form=near.name, nearest_rms=near_rms, best_grid=best.name)
            print(f'B {fam} {readout}: own {scores(res[("own", "population")], data["noise_free"].data)["pattern"]:.4f} '
                  f'chosen {scores(res[("chosen", "population")], data["noise_free"].data)["pattern"]:.4f} '
                  f'best {best.name} {best_pat:.4f} nearest {near.name} {near_rms:.4f}', flush=True)
    pd.DataFrame(rows).to_csv(out / 'model_forms.csv', index=False)
    pd.DataFrame(tb).to_csv(out / 'two_by_two.csv', index=False)
    v = dict(chosen_form=f_ch.name, by_readout=vb,
             verdict=('in the responses' if any(x['in_responses'] for x in vb.values()) and not any(x['in_form'] for x in vb.values())
                      else 'in the form' if any(x['in_form'] for x in vb.values()) and not any(x['in_responses'] for x in vb.values())
                      else 'both'))
    (out / 'verdict_B.json').write_text(json.dumps(v, indent=1) + '\n')
    (out / 'provenance_B.json').write_text(json.dumps(dict(script='experiments/E87_inhibition_form.py', part='B', **H.code_commits(),
        families={k: [list(x) for x in v2] for k, v2 in model_families().items()}), indent=1, default=str) + '\n')
    return v


# ---------------------------------------------------------------------------
# part C: a target that needs no APL

def _tnt_flies(keys: Sequence[tuple]) -> List[tuple]:
    return sorted({(k[0], k[1]) for k in keys})


def _halves_of_flies(keys: Sequence[tuple], rng) -> Tuple[np.ndarray, np.ndarray]:
    """Row indices of `keys` split by fly (a fly's hemispheres together) into two halves at random."""
    flies = _tnt_flies(keys)
    perm = [flies[i] for i in rng.permutation(len(flies))]
    a = set(perm[:len(perm) // 2])
    ia = np.array([n for n, k in enumerate(keys) if (k[0], k[1]) in a])
    ib = np.array([n for n, k in enumerate(keys) if (k[0], k[1]) not in a])
    return ia, ib


def _sb_ceiling(r: float) -> float:
    return float(np.sqrt(2 * r / (1 + r))) if r > 0 else float('nan')


def within_target_stats(unit: str, reading: str, form: Form, tnt: Optional[Sequence[tuple]] = None,
                        n_boot: int = 1000, n_split: int = 200, seed: int = 0) -> dict:
    """The effect of `form` within the recordings with APL output blocked (its distances minus the recordings' own,
    the form's c at the control flies' rate): the target (the mean over recordings), its floor (the pattern part and
    loss of the mean over resampled recordings against it, c held) and the agreement of two halves of the flies (ρ,
    Spearman-Brown-corrected to the ceiling). Resampling is over the per-recording rows."""
    rd = Reading(unit, reading, ctrl=[], tnt=tnt)
    rows, c = rd.predict_rows(form)
    eff = rows - rd.rec_off
    t = _nanmean_rows(eff)
    rng = np.random.default_rng(seed)
    n = eff.shape[0]
    fl_p, fl_l = [], []
    for _ in range(n_boot):
        e = _nanmean_rows(eff[rng.integers(0, n, n)]) - t
        ok = ~np.isnan(e)
        fl_l.append(float(np.sqrt(np.mean(e[ok] ** 2))))
        fl_p.append(float(np.std(e[ok])))
    rs = []
    for _ in range(n_split):
        ia, ib = _halves_of_flies(rd.tnt_keys, rng)
        ta, tb = _nanmean_rows(eff[ia]), _nanmean_rows(eff[ib])
        ok = ~(np.isnan(ta) | np.isnan(tb))
        rs.append(float(H.spearman(pd.Series(ta[ok]), pd.Series(tb[ok]))))
    r = float(np.nanmean(rs))
    return dict(target=t, c=c, floor_pattern=float(np.mean(fl_p)), floor_loss=float(np.mean(fl_l)), half_rho=r,
                ceiling=_sb_ceiling(r), mean=float(np.nanmean(t)), sd=float(np.nanstd(t)), n_tnt=n)


def main_target_stats(unit: str, reading: str, n_boot: int = 1000, n_split: int = 200, seed: int = 0) -> dict:
    """The main target in the same convention (the thresholded data effect): the floor of its pattern part (the
    recordings of both genotypes resampled) and the Spearman-Brown ceiling of two halves of each genotype's flies."""
    rd = Reading(unit, reading)
    rng = np.random.default_rng(seed)
    nc, nt = rd.rec_on.shape[0], rd.rec_off.shape[0]
    fl = []
    for _ in range(n_boot):
        e = (_nanmean_rows(rd.rec_on[rng.integers(0, nc, nc)]) - _nanmean_rows(rd.rec_off[rng.integers(0, nt, nt)])) - rd.data
        ok = ~np.isnan(e)
        fl.append(float(np.std(e[ok])))
    rs = []
    for _ in range(n_split):
        ca, cb = _halves_of_flies(rd.ctrl_keys, rng)
        ta, tb = _halves_of_flies(rd.tnt_keys, rng)
        da = _nanmean_rows(rd.rec_on[ca]) - _nanmean_rows(rd.rec_off[ta])
        db = _nanmean_rows(rd.rec_on[cb]) - _nanmean_rows(rd.rec_off[tb])
        ok = ~(np.isnan(da) | np.isnan(db))
        rs.append(float(H.spearman(pd.Series(da[ok]), pd.Series(db[ok]))))
    r = float(np.nanmean(rs))
    return dict(floor_pattern=float(np.mean(fl)), half_rho=r, ceiling=_sb_ceiling(r), mean=float(np.nanmean(rd.data)),
                sd=float(np.nanstd(rd.data)))


def _model_within(mr: 'ModelResponses', form) -> np.ndarray:
    return mr.form(form)[0] - mr.d_off


def _e85_model_flies() -> List[Tuple[str, Dict[str, tuple]]]:
    """E85's 72 model flies (B1: each wiring with each ORN fly's reproducible input, both response rates held), from
    outputs/E85_orn_fly_population/B/ (viki's)."""
    d = H.outputs_dir('E85_orn_fly_population') / 'B'
    out = []
    for p in sorted(d.glob('B1[ST]_repro_f*_*_on_spike_counts.parquet')):
        tag = p.name[:-len('_on_spike_counts.parquet')]
        q = d / f'{tag}_acute_off_spike_counts.parquet'
        if q.exists():
            out.append((tag, {'on': (None, pd.read_parquet(p)), 'off': (None, pd.read_parquet(q))}))
    return out


def archive_variants() -> List[Tuple[str, str, Path, Path]]:
    """Every archived variant with APL-on and acute APL-off spike counts on the 16 APL odors: [(folder, variant, on, off)],
    E86's excepted."""
    root = H.shared_dir() / 'data' / 'model_outputs'
    out = []
    for off in sorted(root.rglob('*_acute_off_spike_counts.parquet')):
        rel = off.relative_to(root)
        folder = rel.parts[0]
        if folder in ('inputs', 'E86_reduced_input_candidate'):
            continue
        on = off.with_name(off.name.replace('_acute_off_spike_counts.parquet', '_on_spike_counts.parquet'))
        if on.exists():
            out.append((folder, off.name[:-len('_acute_off_spike_counts.parquet')], on, off))
    return out


def part_c(out: Path) -> dict:
    """C1 (apl_free_target.csv: the within-fly effect of the uniform form, the form A2 chose and q 1 / p 1, as trial
    means and with the trial noise removed, in both units, with the six typical recordings beside, its floor, the halves'
    agreement and ceiling, beside the main target's in the same convention); C2 (apl_free_models.csv: both families'
    wirings and populations and E85's 72 model flies against it); C3 (archive_scan.csv: every archived variant); C's
    verdict (verdict_C.json)."""
    import json
    from experiments import E77_input_side_common_error as E77
    f_ch = chosen_form(out)
    forms = {'uniform': UNIFORM, 'chosen': f_ch, 'q1_p1': Form('subtractive', 1.0, 1.0, 'count')}
    rows, tg = [], {}
    strong = [tuple(k) for k in E77.STRONG_TNT]
    for unit in UNITS:
        for reading in ('trial_means', 'noise_free'):
            mt = main_target_stats(unit, reading)
            for lab, f in forms.items():
                for subset in ('all', 'six_typical'):
                    tnt = None if subset == 'all' else [k for k in trialmean_mats(unit)['TNT'] if k not in strong]
                    st = within_target_stats(unit, reading, f, tnt=tnt)
                    if subset == 'all':
                        tg[(unit, reading, lab)] = st['target']
                    rows.append(dict(unit=unit, reading=reading, form_label=lab, form=f.name, subset=subset,
                                     **{k: v for k, v in st.items() if k != 'target'},
                                     main_floor_pattern=mt['floor_pattern'], main_ceiling=mt['ceiling'], main_half_rho=mt['half_rho'],
                                     usable=bool(st['floor_pattern'] <= C_FLOOR_RATIO * mt['floor_pattern'] and st['ceiling'] >= C_CEILING)))
                    print(f'C1 {unit} {reading} {lab} {subset}: floor {st["floor_pattern"]:.4f} (main {mt["floor_pattern"]:.4f}) '
                          f'ceiling {st["ceiling"]:.3f} (main {mt["ceiling"]:.3f})', flush=True)
    c1 = pd.DataFrame(rows)
    c1.to_csv(out / 'apl_free_target.csv', index=False)
    pd.concat([pd.DataFrame(dict(unit=u, reading=r, form_label=l, effect=_pairs(v))).reset_index() for (u, r, l), v in tg.items()],
              ignore_index=True).to_csv(out / 'apl_free_target_pairs.csv', index=False)
    # C2: the models against it (z, both readings, the chosen form)
    mrows = []
    fams = model_families()
    for readout in ('raw', 'cal'):
        for fam, members in fams.items():
            mrs = [ModelResponses(_load(f, t), readout, label=w) for w, f, t in members]
            for lab, f in forms.items():
                effs = [(mr.label, _model_within(mr, f)) for mr in mrs]
                effs.append(('population', _pop([mr.form(f)[0] for mr in mrs]) - _pop([mr.d_off for mr in mrs])))
                for unit_lab, e in effs:
                    for reading in ('trial_means', 'noise_free'):
                        mrows.append(dict(source=fam, readout=readout, unit=unit_lab, form_label=lab, reading=reading,
                                          **scores(e, tg[('z', reading, lab)])))
        flies = _e85_model_flies()
        if flies:
            mrs = [ModelResponses(o, readout, label=tag) for tag, o in flies]
            for lab, f in forms.items():
                pop = _pop([mr.form(f)[0] for mr in mrs]) - _pop([mr.d_off for mr in mrs])
                for reading in ('trial_means', 'noise_free'):
                    mrows.append(dict(source='E85_model_flies', readout=readout, unit='population', form_label=lab, reading=reading,
                                      **scores(pop, tg[('z', reading, lab)])))
                    for mr in mrs:
                        mrows.append(dict(source='E85_model_flies', readout=readout, unit=mr.label, form_label=lab, reading=reading,
                                          **scores(_model_within(mr, f), tg[('z', reading, lab)])))
        print(f'C2 {readout}', flush=True)
    c2 = pd.DataFrame(mrows)
    c2.to_csv(out / 'apl_free_models.csv', index=False)
    c3 = archive_scan(out, f_ch, tg)
    v = dict(chosen_form=f_ch.name,
             usable={f'{u}_{lab}': bool(c1[(c1.unit == u) & (c1.form_label == lab) & (c1.subset == 'all')].usable.all())
                     for u in UNITS for lab in forms},
             n_standout=int(c3.standout.sum()) if len(c3) else 0,
             standouts=c3[c3.standout][['folder', 'variant']].values.tolist() if len(c3) else [])
    (out / 'verdict_C.json').write_text(json.dumps(v, indent=1) + '\n')
    (out / 'provenance_C.json').write_text(json.dumps(dict(script='experiments/E87_inhibition_form.py', part='C', **H.code_commits()),
                                                      indent=1, default=str) + '\n')
    return v


def archive_scan(out: Path, f_ch: Form, tg: dict) -> pd.DataFrame:
    """C3: every archived variant (`archive_variants`) through the chosen form against C1's target (z, the trial noise
    removed, and trial means), raw and under the readout; joined with results/apl_effect/by_variant.csv (its own
    APL-effect loss and realism); a variant stands out if its pattern part against the target (raw, noise-free) is lower
    than that of the breadth-2 fit on the same wiring by C_STANDOUT or more, with a mean effect of at least C_MEAN_SHARE of
    the flies', inside the realism bands. archive_scan.csv."""
    from experiments import E72_wiring_population as E72
    rows = []
    vs = archive_variants()
    for n, (folder, variant, on, off) in enumerate(vs):
        try:
            o = {'on': (None, pd.read_parquet(on)), 'off': (None, pd.read_parquet(off))}
            cols = list(H.abbrev_columns(o['on'][1]).columns)
            if not set(H.APL16_ODORS) <= set(cols):
                continue
            for readout in ('raw', 'cal'):
                mr = ModelResponses(o, readout, label=variant)
                if mr.rate_on <= 0 or not (mr.off > mr.crit).any():
                    continue
                e = _model_within(mr, f_ch)
                row = dict(folder=folder, variant=variant, readout=readout, n_kc=int(mr.off.shape[0]), rate_on=mr.rate_on,
                           rate_off=float((mr.off > mr.crit).mean()))
                for reading in ('noise_free', 'trial_means'):
                    s = scores(e, tg[('z', reading, 'chosen')])
                    row.update({f'{reading}_{k}': s[k] for k in ('pattern', 'loss', 'mean', 'rho')})
                rows.append(row)
        except Exception as ex:   # a malformed file is listed, not fatal
            rows.append(dict(folder=folder, variant=variant, error=repr(ex)[:200]))
        if n % 200 == 0:
            print(f'C3 {n} of {len(vs)}', flush=True)
    df = pd.DataFrame(rows)
    bv = pd.read_csv(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv')
    bv = bv.rename(columns={'tag': 'variant', 'spike_counts_cosine__apl_effect_rmse': 'own_apl_loss',
                            'spike_counts_cosine__apl_effect_rmse_centered': 'own_apl_pattern'})
    df = df.merge(bv[['folder', 'variant', 'own_apl_loss', 'own_apl_pattern', 'realism_flag']].drop_duplicates(['folder', 'variant']),
                  on=['folder', 'variant'], how='left')
    # the breadth-2 fit on the same wiring
    ref = {}
    for w, f, t in model_families()['b2']:
        ref[w] = (f, t)
    prefixes = {lab: prefix for lab, prefix in E72.DRAWS}

    def wiring_of(v: str) -> str:
        for lab, prefix in prefixes.items():
            if v.startswith(prefix):
                return lab
        return 'real' if not v.startswith('wiring-') else ''
    df['wiring'] = df.variant.astype(str).map(wiring_of)
    base = {}
    for w, (f, t) in ref.items():
        m = df[(df.folder == f) & (df.variant == t) & (df.readout == 'raw')]
        base[w] = float(m.noise_free_pattern.iloc[0]) if len(m) else np.nan
    df['b2_same_wiring_pattern'] = df.wiring.map(base)
    flies_mean = float(np.nanmean(tg[('z', 'noise_free', 'chosen')]))
    df['standout'] = ((df.readout == 'raw') & (df.noise_free_pattern <= df.b2_same_wiring_pattern - C_STANDOUT)
                      & (df.noise_free_mean >= C_MEAN_SHARE * flies_mean) & df.realism_flag.isna())
    df.to_csv(out / 'archive_scan.csv', index=False)
    ok = df[(df.readout == 'raw') & df.own_apl_loss.notna() & df.noise_free_pattern.notna()]
    if len(ok) > 2:
        print(f'C3: {len(df)} rows; ρ of the pattern against the target with the own APL-effect loss '
              f'{H.spearman(ok.noise_free_pattern, ok.own_apl_loss):.2f}; standouts {int(df.standout.sum())}', flush=True)
    return df


# ---------------------------------------------------------------------------
# part D: the forms as a static inhibition in place of APL (simulations)

D_BRACKET, D_STEP = (1.0, 200.0), 1.3    # bisect_factor on v = 1 + c (c = 0 at v = 1, responses fall with v)


def d_forms(out: Path) -> Dict[str, Form]:
    """D1's five forms: uniform; q 1, p 0; the form A2 chose; q 1, p 1; q 0, p 2 (breadth by count unless the chosen
    form's)."""
    return {'uniform': UNIFORM, 'q1_p0': Form('subtractive', 1.0, 0.0), 'chosen': chosen_form(out),
            'q1_p1': Form('subtractive', 1.0, 1.0), 'q0_p2': Form('subtractive', 0.0, 2.0)}


def _d_worker(conn, label: str, parts: dict, fit: dict, js: List[int]) -> None:
    """A part-D worker: wiring `label`'s simulator (E83's `build_simulator`, the fit's parameters, the lab's input and
    breadth) on the ORN columns `js`, then for each task (a list of per-odor KC threshold factors, or None) the APL-off
    one-odor runs (`mp.sim_only`) of its odors, until None."""
    try:
        from experiments import E83_cross_wiring_fit as E
        E._PARTS[label] = parts
        E._FITS[label] = fit
        sim, info = E.build_simulator(label, orn=parts['orn'].iloc[:, js], breadth=parts['breadth'])
        cells = info['cells']
        mp, order = cells['mp'].cell_contents, np.asarray(cells['odor_order'].cell_contents, dtype=int)
        names = [str(c) for c in cells['odor_index'].cell_contents]
        p_off = dict(info['params'], wAPLKC_scale=0.0, wAPLPN_scale=0.0)
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
            res = []
            for jj in range(len(order)):
                fac = np.ones(n_kc) if msg[jj] is None else np.asarray(msg[jj], dtype=float)
                mp.sim_only = [int(order[jj])]
                try:
                    r, s = sim(p_off, kc_thr_factor=fac)
                finally:
                    mp.sim_only = []
                res.append((names[jj], r.iloc[:, jj].to_numpy(), s.iloc[:, jj].to_numpy()))
            conn.send(('ok', res))
        except Exception as e:
            conn.send(('error', repr(e)))


class WiringPool:
    """One wiring's one-odor APL-off runs spread over worker processes (spawned, one OpenMP thread each), as E80's
    OdorPool does for C2; `run(factors)` takes per-odor KC threshold factors (KCs x odors, the saved outputs' column
    order; None for all ones) and returns (responses, spike counts) in the saved outputs' layout."""

    def __init__(self, label: str, n_workers: int = 16):
        import multiprocessing
        import os
        from experiments import E83_cross_wiring_fit as E
        self.label = label
        parts = E.wiring_parts(label)
        fit = E.saved_fit(label)
        self.template = fit['saved']['off']
        self.cols = [str(c) for c in self.template[0].columns]
        n = len(self.cols)
        self.groups = [list(map(int, g)) for g in np.array_split(np.arange(n), max(1, min(n_workers, n)))]
        send = {k: parts[k] for k in parts}
        ctx = multiprocessing.get_context('spawn')
        keep = {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT')}
        os.environ['OMP_NUM_THREADS'] = os.environ['OMP_THREAD_LIMIT'] = '1'
        self.procs, self.conns = [], []
        try:
            for js in self.groups:
                a, b = ctx.Pipe()
                pr = ctx.Process(target=_d_worker, args=(b, label, send, fit, js), daemon=True)
                pr.start()
                self.procs.append(pr)
                self.conns.append(a)
        finally:
            for k, v in keep.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        for conn, js in zip(self.conns, self.groups):
            status, v = conn.recv()
            if status != 'ready':
                self.close()
                raise RuntimeError(f'{label}: a worker failed to start: {v}')
            if v != [self.cols[j] for j in js]:
                self.close()
                raise RuntimeError(f'{label}: a worker has odors {v}, expected {[self.cols[j] for j in js]}')

    def run(self, factors: Optional[np.ndarray]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        for conn, js in zip(self.conns, self.groups):
            conn.send([None if factors is None else factors[:, j] for j in js])
        R, S = self.template[0].copy(), self.template[1].copy()
        for conn, js in zip(self.conns, self.groups):
            status, v = conn.recv()
            if status != 'ok':
                raise RuntimeError(f'{self.label}: a worker failed: {v}')
            for j, (name, r, s) in zip(js, v):
                if name != self.cols[j]:
                    raise RuntimeError(f'{self.label}: odor {name} where {self.cols[j]} was expected')
                R.iloc[:, j] = r
                S.iloc[:, j] = s
        return R, S

    def close(self) -> None:
        for conn in self.conns:
            try:
                conn.send(None)
            except Exception:
                pass
        for pr in self.procs:
            pr.join(timeout=30)
            if pr.is_alive():
                pr.terminate()
        self.procs, self.conns = [], []


def _d_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT) / 'D'
    d.mkdir(parents=True, exist_ok=True)
    return d


def _save_d(tag: str, r: pd.DataFrame, s: pd.DataFrame) -> None:
    d = _d_dir()
    r.to_parquet(d / f'{tag}_on_responses.parquet')
    s.to_parquet(d / f'{tag}_on_spike_counts.parquet')


def _load_d(tag: str) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    d = _d_dir()
    ps = [d / f'{tag}_on_{k}.parquet' for k in ('responses', 'spike_counts')]
    return tuple(pd.read_parquet(p) for p in ps) if all(p.exists() for p in ps) else None


def static_factors(off_counts: pd.DataFrame, form: Form, c: float) -> np.ndarray:
    """1 + c * A_o^q * B_k^p (KCs x odors, the outputs' column order) from the wiring's own APL-off spike counts (a
    response two or more spikes)."""
    m = off_counts.to_numpy(dtype=float)
    return 1.0 + c * inhibition(m, form, RAW_CRIT)


def part_d1(out: Path, wirings: Optional[Sequence[str]] = None, n_workers: int = 16) -> pd.DataFrame:
    """D1: on each wiring of C2's family, first the identity (the saved APL-off run reproduced through the one-odor
    path), then for each of D1's forms the per-odor KC thresholds raised by 1 + c A^q B^p with APL removed, c bisected
    to the wiring's APL-on response count; the APL-off run the wiring's own. Resumable (runs under outputs/.../D/, rows in
    runs_D.csv). Then the readings (readings_D.csv) and D1's verdict (verdict_D.json)."""
    import json
    import time
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E83_cross_wiring_fit as E
    wirings = list(E.ALL_WIRINGS if wirings is None else wirings)
    forms = d_forms(out)
    runs_path = out / 'runs_D.csv'
    done = set(pd.read_csv(runs_path).tag) if runs_path.exists() else set()
    ident_path = out / 'identity_D.json'
    ident = json.loads(ident_path.read_text()) if ident_path.exists() else {}
    for w in wirings:
        todo = [k for k in forms if f'D1_{k}_{w}' not in done]
        if not todo and w in ident:
            continue
        fit = E.saved_fit(w)
        pool = WiringPool(w, n_workers)
        try:
            if w not in ident:
                r, s = pool.run(None)
                ns = int((s.to_numpy(dtype=float) != fit['saved']['off'][1].to_numpy(dtype=float)).sum())
                nr = int((r.to_numpy(dtype=float) != fit['saved']['off'][0].to_numpy(dtype=float)).sum())
                ident[w] = dict(count_diffs=ns, response_diffs=nr, **H.code_commits())
                ident_path.write_text(json.dumps(ident, indent=1, default=str) + '\n')
                print(f'D1 identity {w}: {ns} count / {nr} response differences', flush=True)
                if ns or nr:
                    raise RuntimeError(f'{w}: the one-odor path does not reproduce the saved APL-off run')
            off = fit['saved']['off']
            for k in todo:
                f = forms[k]
                tag = f'D1_{k}_{w}'
                t0 = time.time()
                v, (r, s), hist, at_b = E69.bisect_factor(lambda v: pool.run(static_factors(off[1], f, v - 1.0)),
                                                          fit['target_on'], fit['size'], D_BRACKET, D_STEP, f'{tag} c')
                _save_d(tag, r, s)
                o = {'on': (r, s), 'off': off}
                row = dict(part='D1', form_label=k, form=f.name, wiring=w, tag=tag, c=v - 1.0, n_on=int((r.to_numpy() > 0.5).sum()),
                           target_on=fit['target_on'], at_bound=bool(at_b), n_sims=len(hist), seconds=time.time() - t0,
                           fit_folder=fit['folder'], fit_tag=fit['tag'], **{kk: fit['p'][kk] for kk in E.PARAM_KEYS},
                           K_fit=fit['K'], **H.code_commits())
                row.update(E.realism(o))
                pd.DataFrame([row]).to_csv(runs_path, mode='a', header=not runs_path.exists(), index=False)
                print(f'D1 {tag}: c {v - 1.0:.4g}, {row["n_on"]} / {fit["target_on"]} responses, {len(hist)} sims, '
                      f'{time.time() - t0:.0f}s', flush=True)
        finally:
            pool.close()
    df = score_d1(out)
    comparison_table(out)
    return df


def score_d1(out: Path) -> pd.DataFrame:
    """D1's readings: per form, each wiring and the population (S, T, all eight), E83's `effect_scores` under the
    readout and on raw counts against all of E83's targets; the agreement with part B's count-level effect of the same
    form (raw counts, z with the trial noise removed in B's convention is not used here: the correlation of the two
    per-pair effects and the difference of their pattern parts against the main target, both on raw counts); the
    verdict. readings_D.csv, verdict_D.json."""
    import json
    from experiments import E83_cross_wiring_fit as E
    runs = pd.read_csv(out / 'runs_D.csv')
    forms = d_forms(out)
    rows, vd = [], {}
    pops = {'S': E.S_WIRINGS, 'T': E.T_WIRINGS, 'all8': E.ALL_WIRINGS}
    for k, f in forms.items():
        ws = [w for w in E.ALL_WIRINGS if f'D1_{k}_{w}' in set(runs.tag)]
        dist = {ro: {} for ro in ('cal', 'raw')}
        base = {ro: {} for ro in ('cal', 'raw')}
        proxy = {}
        for w in ws:
            fit = E.saved_fit(w)
            r, s = _load_d(f'D1_{k}_{w}')
            o = {'on': (r, s), 'off': fit['saved']['off']}
            for ro in ('cal', 'raw'):
                dist[ro][w] = E.pair_distances(o, ro, fit['K'] if ro == 'cal' else None)
                base[ro][w] = E.pair_distances(fit['saved'], ro, fit['K'] if ro == 'cal' else None)
            mr = ModelResponses(fit['saved'], 'raw', label=w)
            proxy[w] = mr.form(f)[0]
        for ro in ('cal', 'raw'):
            units = {w: [w] for w in ws}
            units.update({p: list(m) for p, m in pops.items() if all(w in ws for w in m)})
            for u, members in units.items():
                d = E.population([dist[ro][w] for w in members])
                d0 = E.population([base[ro][w] for w in members])
                sc, sc0 = E.effect_scores(d), E.effect_scores(d0)
                rows.append(dict(form_label=k, form=f.name, readout=ro, unit=u, n_wirings=len(members),
                                 **sc, **{f'unchanged_{kk}': vv for kk, vv in sc0.items()},
                                 mean_on=float(d['on'].mean()), mean_off=float(d['off'].mean())))
        if all(w in ws for w in E.ALL_WIRINGS):
            r8 = [x for x in rows if x['form_label'] == k and x['unit'] == 'all8']
            cal, raw = next(x for x in r8 if x['readout'] == 'cal'), next(x for x in r8 if x['readout'] == 'raw')
            singles = [x for x in rows if x['form_label'] == k and x['readout'] == 'cal' and x['unit'] in E.ALL_WIRINGS]
            n_lower = sum(x['z8_pattern'] < x['unchanged_z8_pattern'] for x in singles)
            # the proxy: the simulated against the count-level (part B, raw counts, thresholded convention) population effect
            sim_eff = E.population([dist['raw'][w] for w in E.ALL_WIRINGS])
            sim_e = _cs_effect(sim_eff)
            cnt = _pop([proxy[w] for w in E.ALL_WIRINGS]) - _pop([ModelResponses(E.saved_fit(w)['saved'], 'raw').d_off for w in E.ALL_WIRINGS])
            cnt_e = _pairs(cnt)
            ok = sim_e.notna() & cnt_e.notna()
            r_proxy = float(np.corrcoef(sim_e[ok], cnt_e[ok])[0, 1])
            tg = E.targets()['z8']
            pat_sim = scores(_pairs_to_triu(sim_e), _pairs_to_triu(tg))['pattern']
            pat_cnt = scores(cnt, _pairs_to_triu(tg))['pattern']
            vd[k] = dict(form=f.name, pattern_cal=cal['z8_pattern'], unchanged=D_UNCHANGED, n_wirings_lower=int(n_lower),
                         raw_loss=raw['z8_loss'], better=bool(cal['z8_pattern'] <= D_UNCHANGED - D_GAIN and n_lower >= D_WIRINGS
                                                              and raw['z8_loss'] <= D_RAW_UNCHANGED + D_RAW_SLACK),
                         worse=bool(cal['z8_pattern'] >= D_UNCHANGED + D_GAIN), proxy_r=r_proxy,
                         proxy_pattern_diff=float(abs(pat_sim - pat_cnt)),
                         proxy_stands=bool(r_proxy >= D_PROXY_R and abs(pat_sim - pat_cnt) < D_PROXY_DIFF))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'readings_D.csv', index=False)
    (out / 'verdict_D.json').write_text(json.dumps(vd, indent=1) + '\n')
    return df


def _cs_effect(d: Dict[str, pd.Series]) -> pd.Series:
    from experiments import E83_cross_wiring_fit as E
    return E._cs(H.apl_effect(d['on'], d['off']))


def _pairs_to_triu(s: pd.Series) -> np.ndarray:
    """A canonical pair Series as a 120-vector in triu order of APL16_ODORS (the inverse of `_pairs`)."""
    i, j = _iu()
    o = H.APL16_ODORS
    keys = [tuple(sorted((o[a], o[b]))) for a, b in zip(i, j)]
    s = s.copy()
    s.index = [tuple(sorted(x)) for x in s.index]
    return np.array([s.get(k, np.nan) for k in keys], dtype=float)


# ---------------------------------------------------------------------------
# the comparison table

def _f(x, spec: str = '.4f') -> str:
    try:
        return '–' if x is None or (isinstance(x, float) and np.isnan(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def comparison_table(out: Path) -> str:
    """comparison_table.md: A's forms (z, all recordings), the held-out readings, B's two-by-two (breadth-2 population,
    the trial noise removed), C's target, D1's populations and the verdicts, from the parts' tables."""
    import json
    lines = ['# E87: the form of APL\'s inhibition in the flies\' own responses', '',
             f'Generated by experiments/E87_inhibition_form.py at {H.code_commits().get("queue_commit", "")[:7]}. All values '
             'z with sub-threshold responses set to 0 unless stated; lower is better for the pattern part.', '']
    f = pd.read_csv(out / 'forms.csv')
    keep = ['none', 'sub_q0_p0', 'div_q0_p0', 'per_odor_oracle_sub', 'sub_q1_p0', 'sub_q1_p0.25_count', 'sub_q1_p0.25_summed',
            'sub_q1_p1_count', 'sub_q0_p2_count', 'similar_pow1']
    tm = f[(f.unit == 'z') & (f.reading == 'trial_means')].set_index('form')
    nf = f[(f.unit == 'z') & (f.reading == 'noise_free')].set_index('form')
    lines += ['## A: the forms on all recordings', '',
              '| form | mean | pattern part | ρ | r² / slope | ester / alcohol error | four nearest | pattern, trial noise removed | odors per KC TV |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for k in keep:
        if k not in tm.index:
            continue
        r = tm.loc[k]
        lines.append(f"| {k} | {_f(r['mean'], '+.3f')} | {_f(r.pattern)} | {_f(r.rho, '.2f')} | {_f(r.r2, '.2f')} / {_f(r.slope, '.2f')} | "
                     f"{_f(r.err_ester, '+.3f')} / {_f(r.err_alcohol, '+.3f')} | {_f(r.near4, '+.3f')} | "
                     f"{_f(nf.loc[k].pattern) if k in nf.index else '–'} | {_f(r.get('opk_tv'), '.3f')} |")
    ho = pd.read_csv(out / 'held_out.csv')
    lines += ['', '## A2: held out (pattern part on the unseen half)', '',
              '| unit | reading | split | direction | chosen | chosen | uniform | strengths permuted | p 0 (loss) | chosen (loss) | similar tuning |',
              '|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|']
    for _, r in ho.iterrows():
        lines.append(f'| {r.unit} | {r.reading} | {r.split} | {r.direction} | {r.chosen_form} | {_f(r.chosen_pattern)} | '
                     f'{_f(r.uniform_pattern)} | {_f(r.chosen_strengths_perm_pattern)} | {_f(r.p0_loss)} | {_f(r.chosen_loss)} | '
                     f'{_f(r.similar_pattern)} |')
    if (out / 'held_out_random.csv').exists():
        rr = pd.read_csv(out / 'held_out_random.csv')
        lines += ['', f'Random splits (z, trial means, {len(rr)} readings): chosen {_f(rr.chosen_pattern.mean())}, uniform '
                  f'{_f(rr.uniform_pattern.mean())}; the chosen form lower by 0.005 or more in '
                  f'{_f(100 * ((rr.uniform_pattern - rr.chosen_pattern) >= HOLD_GAIN).mean(), ".0f")} %.']
    if (out / 'two_by_two.csv').exists():
        tb = pd.read_csv(out / 'two_by_two.csv')
        lines += ['', '## B: the two-by-two (populations, the data\'s effect with the trial noise removed)', '',
                  '| family | readout | responses | inhibition | form | mean | pattern part | ρ | ester error |',
                  '|---|---|---|---|---|---:|---:|---:|---:|']
        for _, r in tb[tb.data_reading == 'noise_free'].iterrows():
            lines.append(f'| {r.family} | {r.readout} | {r.responses} | {r.inhibition} | {r.form} | {_f(r["mean"], "+.3f")} | '
                         f'{_f(r.pattern)} | {_f(r.rho, ".2f")} | {_f(r.err_ester, "+.3f")} |')
    if (out / 'apl_free_target.csv').exists():
        c1 = pd.read_csv(out / 'apl_free_target.csv')
        lines += ['', '## C: the target without APL', '',
                  '| unit | reading | form | subset | mean | floor | halves ρ | ceiling | main floor | main ceiling | usable |',
                  '|---|---|---|---|---:|---:|---:|---:|---:|---:|---|']
        for _, r in c1.iterrows():
            lines.append(f'| {r.unit} | {r.reading} | {r.form} | {r.subset} | {_f(r["mean"], "+.3f")} | {_f(r.floor_pattern)} | '
                         f'{_f(r.half_rho, ".2f")} | {_f(r.ceiling, ".2f")} | {_f(r.main_floor_pattern)} | {_f(r.main_ceiling, ".2f")} | {r.usable} |')
    if (out / 'readings_D.csv').exists():
        rd = pd.read_csv(out / 'readings_D.csv')
        lines += ['', '## D1: the forms as a static inhibition in place of APL (all recordings, the main target)', '',
                  '| form | readout | unit | pattern part | unchanged | loss | mean APL effect | ρ | mean distance on / off |',
                  '|---|---|---|---:|---:|---:|---:|---:|---:|']
        for _, r in rd[rd.unit.isin(['S', 'T', 'all8'])].iterrows():
            lines.append(f'| {r.form} | {r.readout} | {r.unit} | {_f(r.z8_pattern)} | {_f(r.unchanged_z8_pattern)} | {_f(r.z8_loss)} | '
                         f'{_f(r.mean_effect, "+.3f")} | {_f(r.z8_rho, ".2f")} | {_f(r.mean_on, ".3f")} / {_f(r.mean_off, ".3f")} |')
    lines += ['', '## Verdicts', '']
    for part in ('A', 'B', 'C', 'D'):
        p = out / f'verdict_{part}.json'
        if p.exists():
            lines += [f'### {part}', '', '```', json.dumps(json.loads(p.read_text()), indent=1), '```', '']
    text = '\n'.join(lines) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', choices=('A', 'A2', 'verdictA', 'B', 'C', 'D1', 'D1score', 'table'), required=True)
    ap.add_argument('--wirings', default=None, help='D1: a comma-separated list of wirings (default all eight)')
    ap.add_argument('--workers', type=int, default=16, help='D1: worker processes (one thread each)')
    args = ap.parse_args()
    out = H.results_dir(EXPERIMENT)
    if args.part == 'A':
        part_a(out)
    elif args.part == 'A2':
        part_a2(out)
        verdict_a(out)
    elif args.part == 'verdictA':
        print(verdict_a(out))
    elif args.part == 'B':
        part_b(out)
    elif args.part == 'C':
        part_c(out)
    elif args.part == 'D1':
        part_d1(out, None if args.wirings is None else args.wirings.split(','), args.workers)
    elif args.part == 'D1score':
        score_d1(out)
        comparison_table(out)
    elif args.part == 'table':
        comparison_table(out)


if __name__ == '__main__':
    main()
