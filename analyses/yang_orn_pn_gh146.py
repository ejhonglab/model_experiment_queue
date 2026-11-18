#!/usr/bin/env python3
"""The models' ORN -> PN step against the GH146 PN recordings (Yang's question, 2026-11-17, viki's terminal): how many
versions of the step there are, the model's PN responses beside the GH146 PN dendrites and the ORNs (heatmaps side by
side), their correlations on the measured glomeruli only, and whether another version or setting fits GH146 better.

Data (read only): the GH146 PN dendrites and the pebbled ORNs of al_analysis's sent_to_remy/2023-10-29 tables
(`E90_model_stages.al_table`; 7 GH146 flies, 9 pebbled flies; trial means of dF/F, the megamat panel), the fly mean of
each glomerulus over the flies that have it. The model's ORN input is `kcq.harness.apl16_orn_deltas()` (spike deltas,
36 glomeruli). Everything is read on the 28 glomeruli that GH146, the pebbled ORNs and the model's input share
(`E90_model_stages.shared_glomeruli`, the measured glomeruli) and the 12 megamat odors of the APL panel (GH146 has no
recording of paa, -aPine, euc or B-myr). The model's 18 glomeruli without ORN data and its 8 input glomeruli without
GH146 data enter no correlation (the 18 still feed olfsysm's LNs at their spontaneous rates, as in every model run).

Versions, each a glomerulus x odor matrix of evoked responses:
- orn_input: the model's ORN input itself (Hz), i.e. no transformation;
- orn_dff: the pebbled ORNs' fly-mean dF/F, the measured ORNs in the same indicator as GH146;
- olfsysm_default: olfsysm's PN layer (Kennedy's dynamic ORN -> LN -> PN model) at its defaults, the step of every
  model on the summary pages: the mean PN rate over the odor (0 to 0.5 s) minus the settled pre-odor mean (-0.25 to
  0 s, E98's SETTLED_PRE_WINDOW). Re-simulated here on the ORN, LN and PN layers only (`osm.run_ORN_LN_sims`,
  `osm.run_PN_sims`) of one fit_mb_model run of the bouton base (no KC tuning); the odor-window means are checked
  against E32's saved glomerulus_pn_rates_default.csv and _inhadd1000.csv (IDENTITY_TOL);
- olfsysm_inhadd<x>: the same with the normalization weakened, pn_inhadd raised and pn_inhsc raised with it so that
  the gain at a typical inhibition is kept (E3 / E32 / E77 / E102's ladder: 100, 314, 1000, 3141, 10000, and 10^6, the
  gain fixed); E102's refit at 1000 is the one fitted family that used another setting;
- kennedy_ss: the steady-state form of Kennedy's PN equation used by Yang's own Claude project (Report 07, Equation 2,
  of the page "APL Synapses, Tuned Claws": PN = 200 tanh((D + c) h / 200 * k1 / (k2 + sum_g (nu0 + D_g))), rectified,
  minus its value at D = 0; c 2.941, h 5.3395, k1 368.66, k2 31.41, nu0 13.3 Hz), here on the model's ORN input in Hz
  (no dF/F scale needed), the sum over the 36 input glomeruli (hers was over the 38 of her table);
- olsen_fit: the divisive normalization of Olsen, Bhandawat and Wilson (2010), PN = E^1.5 / (E^1.5 + sigma^1.5 +
  (m S)^1.5) with E the ORN input rectified at 0 and S its sum over the 36 input glomeruli; sigma and m fitted to
  GH146 (its output scale does not enter a correlation);
- orn_power_fit: the ORN input to a power p, signed (sign(D) |D|^p), p fitted in POWER_BOUNDS: below 1 a compression,
  above 1 an expansion, with no normalization across glomeruli;
- olfsysm_fit: olfsysm's PN layer with pn_inhadd, the PN gain and pn_offset fitted to GH146 inside OLFSYSM_BOUNDS
  (pn_inhsc = the gain factor times `inhsc_for(pn_inhadd)`; pn_tanhsc only multiplies pn_inhsc inside the tanh, so it is
  held at its default).
The fits maximize r_all (below): olsen_fit by Nelder-Mead on log sigma, log m; orn_power_fit by a bounded scalar search;
olfsysm_fit by Nelder-Mead on log pn_inhadd, log gain factor, pn_offset from FIT_STARTS (the fit's own point and the
default for the out-of-sample folds).

Out of sample (`loo`), every version alike: for each odor, a fitted version is refitted on the other 11 odors, and the
straight line from the version's 11 training columns to GH146's maps the held-out column onto GH146's units (a
version's own unit, Hz or dF/F or a.u., can change between a fit's folds, so the columns of different folds are
comparable only through these maps, each learnt without its held-out odor). The metrics of the assembled matrix are
the `oos` columns; for a version with nothing fitted they differ from its plain values only through the maps.

Metrics of a version X against GH146 G (both 28 x 12):
- r_all: Pearson r over the 336 entries; rho_all: Spearman rho over them (blind to a monotone map such as GCaMP's
  saturation, applied to all entries alike);
- r_pattern: Pearson r over glomeruli within each odor, mean over odors (the glomerular pattern);
- r_tuning: Pearson r over odors within each glomerulus, mean over glomeruli (each glomerulus's tuning);
- r_totals: Pearson r over odors of each odor's summed response (what the normalization changes most); cv_totals: the
  coefficient of variation of those sums (ddof 1);
- r_similarity: Pearson r between the 66 odor-pair correlations (Pearson over glomeruli) of X and of G (the PN
  representation's odor-odor similarity).
The ceiling of each r metric: G's flies split 3 / 4 in all 35 ways, the metric between the two half means,
Spearman-Brown corrected to the 7 flies (rel = 2 r / (1 + r)), ceiling sqrt(rel) (what a noiseless true model would
reach against the 7 flies' mean; approximate for the mean-of-correlations metrics). Intervals: the 90 % spread of N_BOOT
resamples of the GH146 flies with replacement around the point value (`centered_interval`; seed BOOT_SEED; a glomerulus
missing from all drawn flies is dropped from that resample for every version), of the plain values (lo, hi) and of the
out-of-sample ones (oos_lo, oos_hi);
p_vs_<ref> / oos_p_vs_<ref>: the fraction of resamples in which the version's metric is not above the reference's
(paired). A fitted version's plain value (`value`) is in sample.

Outputs (<shared>/results/yang_orn_pn_gh146/): versions.csv (long: version, sample, glomerulus, odor, value; each
version as computed and out of sample in GH146's units), metrics.csv, ceiling.csv, fits.json, identity.json, totals.csv, heatmaps.png,
metrics.png, provenance.json.

Usage (from the repository root): ./venv/bin/python -u analyses/yang_orn_pn_gh146.py
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcq import harness as H  # noqa: E402
from experiments import E90_model_stages as E90  # noqa: E402
from experiments.E3_drive_dependent_gain import BOUTON, PN_INHADD_DEFAULT, PN_INHSC_DEFAULT, TYPICAL_INH  # noqa: E402

NAME = 'yang_orn_pn_gh146'
ODORS: List[str] = list(E90.ODORS)
ODOR_WINDOW = (0.0, 0.5)          # closed at both ends, as E32's glomerulus_pn_rates
PRE_WINDOW = (-0.25, 0.0)         # E98's SETTLED_PRE_WINDOW: closed on the left, open on the right
IDENTITY_TOL = 0.05               # Hz, odor-window means against E32's saved tables
LADDER = (100.0, 314.0, 1000.0, 3141.0, 10000.0, 1e6)
PN_OFFSET_DEFAULT = 2.9410
KENNEDY = dict(c=2.941, h=5.3395, k1=368.66, k2=31.41, nu0=13.3)
OLSEN_N = 1.5
# olfsysm_fit's bounds: pn_inhadd from below the default to the ladder's top (10^6: the gain fixed), the PN gain at a
# typical inhibition 1/20 to 20 times the default's (the tanh from nearly linear to fully saturated), the offset +-30 Hz;
# without bounds the fit runs off to pn_inhadd 1e43 and pn_inhsc 1e-35 (a vanishing gain: the tanh linear, every
# rate ~1e-75 Hz), the same shape as gain_fixed at a gain no rate can be read at
OLFSYSM_BOUNDS = dict(pn_inhadd=(1.0, 1e6), gain_factor=(0.05, 20.0), pn_offset=(-30.0, 30.0))
FIT_STARTS = ((PN_INHADD_DEFAULT, 1.0, PN_OFFSET_DEFAULT), (1000.0, 1.0, PN_OFFSET_DEFAULT),
              (1e5, 1.0, PN_OFFSET_DEFAULT), (1e5, 0.2, 0.0))
POWER_BOUNDS = (0.2, 4.0)
N_BOOT = 1000
BOOT_SEED = 11
R_METRICS = ('r_all', 'rho_all', 'r_pattern', 'r_tuning', 'r_totals', 'r_similarity')


def inhsc_for(inhadd: float) -> float:
    """pn_inhsc raised with pn_inhadd so that the gain at a typical inhibition is kept (E3 / E32 / E77 / E102)."""
    return float(PN_INHSC_DEFAULT * (inhadd + TYPICAL_INH) / (PN_INHADD_DEFAULT + TYPICAL_INH))


# ---------------------------------------------------------------------------------------------------------- the data

def fly_mean(table: pd.DataFrame, glomeruli: Sequence[str]) -> pd.DataFrame:
    """Glomeruli x odors: the mean over the flies that have each glomerulus (rows of `table`: fly, odor)."""
    m = table.groupby(level='odor').mean()
    return m.loc[ODORS, list(glomeruli)].T


def gh146_flies(gh: pd.DataFrame, glomeruli: Sequence[str]) -> np.ndarray:
    """flies x glomeruli x odors (NaN where a fly lacks a glomerulus)."""
    flies = list(dict.fromkeys(gh.index.get_level_values('fly')))
    out = np.full((len(flies), len(glomeruli), len(ODORS)), np.nan)
    for i, f in enumerate(flies):
        x = gh.xs(f, level='fly').reindex(index=ODORS, columns=list(glomeruli))
        out[i] = x.to_numpy(dtype=float).T
    return out


# ------------------------------------------------------------------------------------------------- the model's PN layer

class PNLayer:
    """olfsysm's ORN, LN and PN layers re-run on one fit_mb_model run's ModelParams / RunVars."""

    def __init__(self, orn: pd.DataFrame):
        import olfsysm as osm
        from pyolfsysm.mb_model import fit_mb_model, get_time_index
        self.osm = osm
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _, _, wPNKC, pd0 = fit_mb_model(orn, fixed_thr=1e4, wAPLKC=1.0, wKCAPL=1.0, wAPLPN=1.0, wPNAPL=1.0,
                                            **BOUTON, make_plots=False, silent=True, return_olfsysm_vars=True)
        self.mp, self.rv = pd0['mp'], pd0['rv']
        self.glomeruli = list(dict.fromkeys(wPNKC.columns.get_level_values('glomerulus')))
        self.odors = [H.odor_abbrev(c) for c in orn.columns]
        t = np.asarray(get_time_index(self.mp, self.rv))
        self.win = (t >= ODOR_WINDOW[0]) & (t <= ODOR_WINDOW[1])
        self.pre = (t >= PRE_WINDOW[0]) & (t < PRE_WINDOW[1])
        self.defaults = dict(inhadd=self.mp.pn.inhadd, inhsc=self.mp.pn.inhsc, offset=self.mp.pn.offset,
                             tanhsc=self.mp.pn.tanhsc)
        assert abs(self.defaults['inhadd'] - PN_INHADD_DEFAULT) < 1e-9 and abs(self.defaults['inhsc'] - PN_INHSC_DEFAULT) < 1e-9
        assert abs(self.defaults['offset'] - PN_OFFSET_DEFAULT) < 1e-9

    def run(self, inhadd: float, inhsc: float, offset: float = PN_OFFSET_DEFAULT) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """(odor-window means, evoked = odor-window minus settled pre-odor means), glomeruli x the 16 odors (Hz)."""
        from pyolfsysm.mb_model import get_sim_var
        mp = self.mp
        mp.pn.inhadd, mp.pn.inhsc, mp.pn.offset = float(inhadd), float(inhsc), float(offset)
        try:
            self.osm.run_ORN_LN_sims(mp, self.rv)
            self.osm.run_PN_sims(mp, self.rv)
            pn = np.asarray(get_sim_var(self.rv.pn, 'pn_sims'))
        finally:
            mp.pn.inhadd, mp.pn.inhsc, mp.pn.offset = (self.defaults[k] for k in ('inhadd', 'inhsc', 'offset'))
        on = pn[:, :, self.win].mean(axis=2).T
        pre = pn[:, :, self.pre].mean(axis=2).T
        mk = lambda a: pd.DataFrame(a, index=self.glomeruli, columns=self.odors)
        return mk(on), mk(on - pre)


def kennedy_ss(orn: pd.DataFrame) -> pd.DataFrame:
    """Yang's Claude's steady-state Kennedy step on the ORN input (glomeruli x odors, Hz), minus its value at D = 0."""
    k = KENNEDY
    d = orn.to_numpy(dtype=float)
    den = k['k2'] + (k['nu0'] + d).sum(axis=0, keepdims=True)
    pn = np.maximum(200.0 * np.tanh((d + k['c']) * k['h'] / 200.0 * k['k1'] / den), 0.0)
    den0 = k['k2'] + k['nu0'] * d.shape[0]
    pn0 = max(200.0 * np.tanh(k['c'] * k['h'] / 200.0 * k['k1'] / den0), 0.0)
    return pd.DataFrame(pn - pn0, index=orn.index, columns=orn.columns)


def olsen(orn: pd.DataFrame, sigma: float, m: float) -> pd.DataFrame:
    e = np.maximum(orn.to_numpy(dtype=float), 0.0)
    s = e.sum(axis=0, keepdims=True)
    n = OLSEN_N
    return pd.DataFrame(e ** n / (e ** n + sigma ** n + (m * s) ** n), index=orn.index, columns=orn.columns)


# --------------------------------------------------------------------------------------------------------- the metrics

def _r(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def _pair_corrs(x: np.ndarray) -> np.ndarray:
    """The 66 odor-pair Pearson correlations over glomeruli (x: glomeruli x odors, no NaN)."""
    c = np.corrcoef(x.T)
    iu = np.triu_indices(c.shape[0], 1)
    return c[iu]


def metrics(x: np.ndarray, g: np.ndarray) -> Dict[str, float]:
    """x, g: glomeruli x odors; rows of g with any NaN are dropped from both."""
    keep = np.isfinite(g).all(axis=1) & np.isfinite(x).all(axis=1)
    x, g = x[keep], g[keep]
    tx, tg = x.sum(axis=0), g.sum(axis=0)
    return dict(
        r_all=_r(x.ravel(), g.ravel()),
        rho_all=float(stats.spearmanr(x.ravel(), g.ravel()).correlation),
        r_pattern=float(np.nanmean([_r(x[:, j], g[:, j]) for j in range(x.shape[1])])),
        r_tuning=float(np.nanmean([_r(x[i], g[i]) for i in range(x.shape[0])])),
        r_totals=_r(tx, tg),
        r_similarity=_r(_pair_corrs(x), _pair_corrs(g)),
        cv_totals=float(tx.std(ddof=1) / tx.mean()) if tx.mean() != 0 else float('nan'),
        n_glomeruli=int(keep.sum()))


def ceiling(flies: np.ndarray) -> pd.DataFrame:
    """Split-half ceilings of the r metrics: all 3 / 4 splits of the flies."""
    n = flies.shape[0]
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for a in combinations(range(n), 3):
            b = [i for i in range(n) if i not in a]
            ga, gb = np.nanmean(flies[list(a)], axis=0), np.nanmean(flies[b], axis=0)
            rows.append(metrics(ga, gb))
    df = pd.DataFrame(rows)
    out = []
    for k in R_METRICS:
        r = float(df[k].mean())
        rel = 2 * r / (1 + r) if r > -1 else float('nan')
        out.append(dict(metric=k, split_half_r=r, spearman_brown=rel, ceiling=float(np.sqrt(rel)) if rel > 0 else float('nan'),
                        n_splits=len(df)))
    return pd.DataFrame(out)


def centered_interval(value: float, boot: pd.Series) -> Tuple[float, float]:
    """The 90 % spread of the bootstrap around the point value: value - (mean - 5th percentile), value + (95th
    percentile - mean). A resample with repeated flies has a noisier fly mean, so every r against it is lower than
    against the 7 flies' mean, and the raw percentiles sit below the point value; the spread is what they carry."""
    b = np.asarray(boot, dtype=float)
    m = float(np.nanmean(b))
    return float(value - (m - np.nanpercentile(b, 5))), float(value + (np.nanpercentile(b, 95) - m))


# -------------------------------------------------------------------------------------------------------------- fits

def _neg_r(pred: pd.DataFrame, g: pd.DataFrame, odors: Sequence[str]) -> float:
    r = _r(pred.loc[g.index, list(odors)].to_numpy().ravel(), g.loc[:, list(odors)].to_numpy().ravel())
    return 1.0 if not np.isfinite(r) else -r


def fit_olsen(orn: pd.DataFrame, g: pd.DataFrame, odors: Sequence[str]) -> Dict[str, float]:
    best = None
    for s0, m0 in ((12.0, 0.05), (50.0, 0.02), (5.0, 0.2)):
        f = lambda p: _neg_r(olsen(orn, float(np.exp(p[0])), float(np.exp(p[1]))), g, odors)
        res = optimize.minimize(f, x0=[np.log(s0), np.log(m0)], method='Nelder-Mead',
                                options=dict(xatol=1e-4, fatol=1e-7, maxiter=2000))
        if best is None or res.fun < best.fun:
            best = res
    return dict(sigma=float(np.exp(best.x[0])), m=float(np.exp(best.x[1])), r_train=float(-best.fun), nfev=int(best.nfev))


def signed_power(orn: pd.DataFrame, p: float) -> pd.DataFrame:
    x = orn.to_numpy(dtype=float)
    return pd.DataFrame(np.sign(x) * np.abs(x) ** p, index=orn.index, columns=orn.columns)


def fit_power(orn: pd.DataFrame, g: pd.DataFrame, odors: Sequence[str]) -> Dict[str, float]:
    res = optimize.minimize_scalar(lambda p: _neg_r(signed_power(orn, p), g, odors), bounds=POWER_BOUNDS, method='bounded',
                                   options=dict(xatol=1e-4))
    return dict(power=float(res.x), r_train=float(-res.fun), at_bound=bool(min(abs(res.x - b) for b in POWER_BOUNDS) < 1e-3))


def olfsysm_point(p: Sequence[float]) -> Tuple[float, float, float]:
    """(pn_inhadd, pn_inhsc, pn_offset) of a fit's point (log pn_inhadd, log gain factor, pn_offset): pn_inhsc is the
    gain factor times `inhsc_for(pn_inhadd)`, so the factor is the PN gain at a typical inhibition over the default's."""
    inhadd = float(np.exp(p[0]))
    return inhadd, float(np.exp(p[1]) * inhsc_for(inhadd)), float(p[2])


def fit_olfsysm(layer: PNLayer, g: pd.DataFrame, odors: Sequence[str],
                starts: Sequence[Tuple[float, float, float]] = FIT_STARTS) -> Dict[str, float]:
    """Nelder-Mead inside OLFSYSM_BOUNDS (a point outside scores as r = -1): pn_inhadd, the gain factor, pn_offset."""
    lo = np.array([np.log(OLFSYSM_BOUNDS['pn_inhadd'][0]), np.log(OLFSYSM_BOUNDS['gain_factor'][0]), OLFSYSM_BOUNDS['pn_offset'][0]])
    hi = np.array([np.log(OLFSYSM_BOUNDS['pn_inhadd'][1]), np.log(OLFSYSM_BOUNDS['gain_factor'][1]), OLFSYSM_BOUNDS['pn_offset'][1]])

    def f(p):
        if np.any(p < lo) or np.any(p > hi):
            return 1.0
        return _neg_r(layer.run(*olfsysm_point(p))[1], g, odors)
    best = None
    for a0, f0, o0 in starts:
        res = optimize.minimize(f, x0=[np.log(a0), np.log(f0), o0], method='Nelder-Mead',
                                options=dict(xatol=1e-3, fatol=1e-6, maxiter=600))
        if best is None or res.fun < best.fun:
            best = res
    inhadd, inhsc, offset = olfsysm_point(best.x)
    p = np.asarray(best.x)
    at = [k for k, i in zip(('pn_inhadd', 'gain_factor', 'pn_offset'), range(3)) if min(p[i] - lo[i], hi[i] - p[i]) < 0.02]
    return dict(pn_inhadd=inhadd, gain_factor=float(np.exp(best.x[1])), pn_inhsc=inhsc, pn_offset=offset,
                r_train=float(-best.fun), nfev=int(best.nfev), at_bounds=at)


def loo(predict: Callable[[Sequence[str]], Tuple[pd.DataFrame, dict]], label: str, g: pd.DataFrame,
        verbose: bool = True) -> Tuple[pd.DataFrame, List[dict]]:
    """Out of sample, in GH146's units: for each odor, `predict(train)` on the other odors gives a matrix (glomeruli x
    odors) and its parameters; the straight line from the matrix's training columns to `g`'s (least squares over their
    entries) maps the held-out column onto g's scale. A matrix is in a version's own units (Hz, dF/F, a.u.), and a fit
    of a fold can change that unit, so the held-out columns of different folds are comparable only through each fold's
    map, learnt without the held-out odor. Returns the assembled predictions and one row per fold."""
    cols, rows = {}, []
    for o in ODORS:
        train = [x for x in ODORS if x != o]
        pred, params = predict(train)
        x = pred.loc[g.index, train].to_numpy(dtype=float).ravel()
        y = g.loc[:, train].to_numpy(dtype=float).ravel()
        ok = np.isfinite(x) & np.isfinite(y)
        a, b = np.polyfit(x[ok], y[ok], 1)
        cols[o] = a * pred.loc[g.index, o] + b
        rows.append(dict(version=label, held_out=o, slope=float(a), intercept=float(b), **params))
        if verbose and params:
            print(f'  {label} without {o}: {params}', flush=True)
    return pd.DataFrame(cols)[ODORS], rows


# --------------------------------------------------------------------------------------------------------- figures

DIVERGING = ['#104281', '#3987e5', '#9ec5f4', '#f0efec', '#f3a59c', '#e34948', '#8f1f1f']
INK, INK2, MUTED, SURFACE = '#0b0b0b', '#52514e', '#8a8984', '#fcfcfb'
SERIES = {'GH146 PN dendrites': INK, 'ORN input (no step)': '#2a78d6', 'model PN, olfsysm default': '#eb6834',
          'best fit': '#1baf7a'}


def heatmaps(panels: List[Tuple[str, str, pd.DataFrame]], path: Path, note: str) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    cmap = LinearSegmentedColormap.from_list('div', DIVERGING)
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(3.05 * n + 0.6, 8.2), facecolor=SURFACE)
    for i, (ax, (title, unit, m)) in enumerate(zip(axes, panels)):
        a = m.to_numpy(dtype=float)
        lim = float(np.nanpercentile(np.abs(a), 99)) or 1.0
        im = ax.imshow(a, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim),
                       interpolation='nearest')
        ax.set_facecolor(SURFACE)
        ax.set_xticks(range(m.shape[1]))
        ax.set_xticklabels(m.columns, rotation=90, fontsize=8, color=INK2)
        ax.set_yticks(range(m.shape[0]))
        ax.set_yticklabels(m.index if i == 0 else [], fontsize=8, color=INK2)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title, fontsize=9.5, color=INK, loc='left', pad=8)
        cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.11, aspect=25)
        cb.set_label(unit, fontsize=8, color=INK2)
        cb.ax.tick_params(labelsize=7, colors=INK2, length=2)
        cb.outline.set_visible(False)
    fig.text(0.01, 0.005, note, fontsize=7.5, color=MUTED, ha='left', va='bottom', wrap=True)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def metrics_figure(met: pd.DataFrame, ceil: pd.DataFrame, totals: pd.DataFrame, order: List[str], labels: Dict[str, str],
                   highlight: Dict[str, str], path: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    show = [('r_all', 'all entries, r'), ('rho_all', 'all entries, ρ'), ('r_pattern', 'pattern within an odor, r'),
            ('r_totals', 'per-odor totals, r'), ('r_similarity', 'odor-odor similarity, r')]
    fig = plt.figure(figsize=(15.5, 4.6 + 0.26 * len(order)), facecolor=SURFACE)
    gs = fig.add_gridspec(2, len(show), height_ratios=[0.26 * len(order) + 0.6, 3.1], hspace=0.55, wspace=0.08)
    y = np.arange(len(order))[::-1]
    for j, (k, title) in enumerate(show):
        ax = fig.add_subplot(gs[0, j])
        ax.set_facecolor(SURFACE)
        for yi, v in zip(y, order):
            row = met[(met.version == v) & (met.metric == k)].iloc[0]
            col = highlight.get(v, INK2)
            ax.plot([row.oos_lo, row.oos_hi], [yi, yi], color=col, lw=2, solid_capstyle='round')
            ax.plot(row.oos, yi, 'o', ms=8, color=col, mec=SURFACE, mew=2)
            if row.fitted:
                ax.plot(row.value, yi, 'o', ms=8, mfc='none', mec=col, mew=1.5)
        c = ceil.set_index('metric').loc[k, 'ceiling']
        if np.isfinite(c):
            ax.axvline(c, color=MUTED, lw=1, ls=(0, (3, 3)))
            ax.text(c, y.max() + 0.7, 'ceiling', fontsize=7, color=MUTED, ha='center', va='bottom')
        ax.set_yticks(y)
        ax.set_yticklabels([labels[v] for v in order] if j == 0 else [], fontsize=8, color=INK2)
        ax.set_ylim(-0.6, y.max() + 1.3)
        ax.set_title(title, fontsize=9, color=INK, loc='left')
        ax.tick_params(labelsize=7.5, colors=INK2, length=2)
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', color='#e6e5e1', lw=0.6)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_color('#c3c2b7')
    ax = fig.add_subplot(gs[1, :])
    ax.set_facecolor(SURFACE)
    x = np.arange(len(totals.index))
    labs = [lab for lab in SERIES if lab in totals.columns]
    for i, lab in enumerate(labs):
        ax.plot(x, totals[lab], '-o', color=SERIES[lab], lw=2, ms=6.5, mec=SURFACE, mew=1.5, label=lab,
                zorder=3 if i == 0 else 2)
    # direct labels at the right end, spread apart where the series end close together
    ends = sorted(((float(totals[lab].iloc[-1]), lab) for lab in labs))
    span = float(np.nanmax(totals.to_numpy()) - np.nanmin(totals.to_numpy()))
    ys = []
    for y0, _ in ends:
        ys.append(max(y0, ys[-1] + 0.075 * span) if ys else y0)
    for (y0, lab), y1 in zip(ends, ys):
        ax.text(x[-1] + 0.25, y1, lab, fontsize=7.5, color=INK2, va='center')
    ax.axhline(1.0, color='#c3c2b7', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(totals.index, fontsize=8, color=INK2)
    ax.set_xlim(-0.4, len(x) + 1.9)
    ax.set_ylabel("odor's summed response / mean over odors", fontsize=8, color=INK2)
    ax.set_title('Per-odor totals over the 28 glomeruli (odors ordered by the GH146 total)', fontsize=9, color=INK, loc='left',
                 pad=26)
    fig.text(0.01, 0.995, 'Filled dots: out of sample (leave one odor out, mapped onto GH146 units); hollow: a fitted step '
             'in sample. Lines: 90 % spread over resampled GH146 flies. Dashed: ceiling (split halves of the 7 flies).',
             fontsize=7.5, color=MUTED, ha='left', va='top')
    ax.tick_params(labelsize=7.5, colors=INK2, length=2)
    ax.legend(fontsize=7.5, frameon=False, loc='lower left', bbox_to_anchor=(0.0, 1.06), ncol=4)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.savefig(path, dpi=160, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------ main

def main() -> None:
    t_start = time.time()
    out = H.shared_dir() / 'results' / NAME
    out.mkdir(parents=True, exist_ok=True)

    gh_tab, pb_tab = E90.al_table('GH146'), E90.al_table('pebbled')
    glom = E90.shared_glomeruli()
    assert len(glom) == 28, len(glom)
    G = fly_mean(gh_tab, glom)
    flies = gh146_flies(gh_tab, glom)
    assert np.allclose(np.nanmean(flies, axis=0), G.to_numpy(), equal_nan=True)
    orn_all = H.apl16_orn_deltas()
    orn = orn_all.copy()
    orn.columns = [H.odor_abbrev(c) for c in orn.columns]
    print(f'GH146: {flies.shape[0]} flies; {len(glom)} glomeruli x {len(ODORS)} odors', flush=True)

    # versions
    V: Dict[str, pd.DataFrame] = {}
    V['orn_input'] = orn.loc[glom, ODORS]
    V['orn_dff'] = fly_mean(pb_tab, glom)
    layer = PNLayer(orn_all)
    on_def, ev_def = layer.run(PN_INHADD_DEFAULT, PN_INHSC_DEFAULT)
    ident = {}
    for tag, (a, s) in {'default': (PN_INHADD_DEFAULT, PN_INHSC_DEFAULT), 'inhadd1000': (1000.0, inhsc_for(1000.0))}.items():
        on = on_def if tag == 'default' else layer.run(a, s)[0]
        saved = pd.read_csv(H.shared_dir() / 'results' / 'E32_pn_stage' / f'glomerulus_pn_rates_{tag}.csv', index_col=0)
        d = float((on.loc[saved.index, saved.columns] - saved).abs().max().max())
        ident[tag] = dict(max_abs_diff_hz=d, tol=IDENTITY_TOL, ok=d < IDENTITY_TOL, pn_inhadd=a, pn_inhsc=s)
        print(f'identity {tag}: max |diff| against E32 {d:.4f} Hz', flush=True)
    (out / 'identity.json').write_text(json.dumps(ident, indent=1) + '\n')
    if not all(v['ok'] for v in ident.values()):
        raise RuntimeError(f'PN layer does not reproduce E32: {ident}')
    V['olfsysm_default'] = ev_def.loc[glom, ODORS]
    for a in LADDER:
        V[f'olfsysm_inhadd{a:g}'] = layer.run(a, inhsc_for(a))[1].loc[glom, ODORS]
    V['kennedy_ss'] = kennedy_ss(orn.loc[:, ODORS]).loc[glom, ODORS]

    # fits, in sample and out of sample
    fits = {}
    t0 = time.time()
    fo = fit_olsen(orn.loc[:, ODORS], G, ODORS)
    fits['olsen_fit'] = fo
    V['olsen_fit'] = olsen(orn.loc[:, ODORS], fo['sigma'], fo['m']).loc[glom, ODORS]
    print(f'olsen_fit: {fo} ({time.time() - t0:.0f} s)', flush=True)
    t0 = time.time()
    ff = fit_olfsysm(layer, G, ODORS)
    fits['olfsysm_fit'] = ff
    V['olfsysm_fit'] = layer.run(ff['pn_inhadd'], ff['pn_inhsc'], ff['pn_offset'])[1].loc[glom, ODORS]
    print(f'olfsysm_fit: {ff} ({time.time() - t0:.0f} s)', flush=True)
    fp = fit_power(orn.loc[:, ODORS], G, ODORS)
    fits['orn_power_fit'] = fp
    V['orn_power_fit'] = signed_power(orn.loc[:, ODORS], fp['power']).loc[glom, ODORS]
    print(f'orn_power_fit: {fp}', flush=True)

    def pred_olsen(train):
        p = fit_olsen(orn.loc[:, ODORS], G, train)
        return olsen(orn.loc[:, ODORS], p['sigma'], p['m']).loc[glom, ODORS], p

    full_start = ((ff['pn_inhadd'], ff['gain_factor'], ff['pn_offset']), (PN_INHADD_DEFAULT, 1.0, PN_OFFSET_DEFAULT))

    def pred_olfsysm(train):
        p = fit_olfsysm(layer, G, train, starts=full_start)
        return layer.run(p['pn_inhadd'], p['pn_inhsc'], p['pn_offset'])[1].loc[glom, ODORS], p

    def pred_power(train):
        p = fit_power(orn.loc[:, ODORS], G, train)
        return signed_power(orn.loc[:, ODORS], p['power']).loc[glom, ODORS], p

    # out of sample, every version through the same leave-one-odor-out map onto GH146's units (the fitted ones also
    # refitted without the held-out odor)
    loo_rows = []
    L = {}
    t0 = time.time()
    fitted = {'olsen_fit': pred_olsen, 'orn_power_fit': pred_power, 'olfsysm_fit': pred_olfsysm}
    for v in V:
        pred = fitted.get(v, lambda train, m=V[v]: (m, {}))
        L[v], r = loo(pred, v, G)
        loo_rows += r
    print(f'leave-one-odor-out: {time.time() - t0:.0f} s', flush=True)
    fits['loo'] = loo_rows

    # metrics, intervals, ceiling
    order = list(V)
    g = G.to_numpy(dtype=float)
    X_in = {v: V[v].loc[glom, ODORS].to_numpy(dtype=float) for v in order}
    X_oos = {v: L[v].loc[glom, ODORS].to_numpy(dtype=float) for v in order}
    point_in = {v: metrics(X_in[v], g) for v in order}
    point_oos = {v: metrics(X_oos[v], g) for v in order}
    rng = np.random.default_rng(BOOT_SEED)
    boot_in = {v: [] for v in order}
    boot_oos = {v: [] for v in order}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for _ in range(N_BOOT):
            idx = rng.integers(0, flies.shape[0], flies.shape[0])
            gb = np.nanmean(flies[idx], axis=0)
            for v in order:
                boot_in[v].append(metrics(X_in[v], gb))
                boot_oos[v].append(metrics(X_oos[v], gb))
    boot_in = {v: pd.DataFrame(b) for v, b in boot_in.items()}
    boot_oos = {v: pd.DataFrame(b) for v, b in boot_oos.items()}
    rows = []
    for v in order:
        for k in R_METRICS + ('cv_totals',):
            bi, bo = boot_in[v][k], boot_oos[v][k]
            lo, hi = centered_interval(point_in[v][k], bi)
            oos_lo, oos_hi = centered_interval(point_oos[v][k], bo)
            row = dict(version=v, metric=k, value=point_in[v][k], lo=lo, hi=hi, oos=point_oos[v][k], oos_lo=oos_lo,
                       oos_hi=oos_hi, boot_mean=float(np.nanmean(bi)), n_glomeruli=point_in[v]['n_glomeruli'],
                       fitted=v in fitted)
            for ref in ('olfsysm_default', 'orn_input', 'orn_dff'):
                if k in R_METRICS and v != ref:
                    row[f'p_vs_{ref}'] = float(np.mean((bi - boot_in[ref][k]).to_numpy() <= 0))
                    row[f'oos_p_vs_{ref}'] = float(np.mean((bo - boot_oos[ref][k]).to_numpy() <= 0))
            rows.append(row)
    met = pd.DataFrame(rows)
    gh_cv = metrics(g, g)['cv_totals']
    gb_cv = [metrics(np.nanmean(flies[rng.integers(0, flies.shape[0], flies.shape[0])], axis=0),
                     np.nanmean(flies, axis=0))['cv_totals'] for _ in range(N_BOOT)]
    met = pd.concat([met, pd.DataFrame([dict(version='gh146', metric='cv_totals', value=gh_cv, oos=gh_cv,
                                             lo=float(np.nanpercentile(gb_cv, 5)), hi=float(np.nanpercentile(gb_cv, 95)),
                                             n_glomeruli=len(glom))])], ignore_index=True)
    ceil = ceiling(flies)
    met.to_csv(out / 'metrics.csv', index=False)
    ceil.to_csv(out / 'ceiling.csv', index=False)

    long = []
    for v, m in V.items():
        s = m.stack().rename('value').reset_index()
        s.columns = ['glomerulus', 'odor', 'value']
        s.insert(0, 'sample', 'in_sample')
        s.insert(0, 'version', v)
        long.append(s)
    for v, m in L.items():
        s = m.loc[glom, ODORS].stack().rename('value').reset_index()
        s.columns = ['glomerulus', 'odor', 'value']
        s.insert(0, 'sample', 'leave_one_odor_out_in_gh146_units')
        s.insert(0, 'version', v)
        long.append(s)
    s = G.stack().rename('value').reset_index()
    s.columns = ['glomerulus', 'odor', 'value']
    s.insert(0, 'sample', 'fly_mean')
    s.insert(0, 'version', 'gh146')
    long.append(s)
    pd.concat(long, ignore_index=True).to_csv(out / 'versions.csv', index=False)
    (out / 'fits.json').write_text(json.dumps(fits, indent=1) + '\n')

    # the best version by r_all (out of sample for the fitted ones), and the figures
    rall = met[met.metric == 'r_all'].set_index('version')['oos'].drop('gh146', errors='ignore')
    transforms = [v for v in order if v.startswith('olfsysm') or v in ('kennedy_ss', 'olsen_fit', 'orn_power_fit')]
    best = str(rall.loc[transforms].idxmax())
    labels = {'orn_input': 'ORN input, no step (Hz)', 'orn_dff': 'ORNs measured, pebbled (dF/F)',
              'olfsysm_default': 'olfsysm PN, default (all models)', 'kennedy_ss': "Kennedy steady state (Yang's Claude)",
              'olsen_fit': 'Olsen 2010, fitted (out of sample)', 'olfsysm_fit': 'olfsysm PN, fitted (out of sample)',
              'orn_power_fit': 'ORN input to a power, fitted (out of sample)'}
    for a in LADDER:
        labels[f'olfsysm_inhadd{a:g}'] = (f'olfsysm PN, pn_inhadd {a:g}' if a < 1e6 else 'olfsysm PN, gain fixed (inhadd 10^6)')
    ord_g = list(pd.Series(np.argmax(G.to_numpy(), axis=1), index=G.index).sort_values(kind='stable').index)
    r_of = lambda v, col='value': met[(met.version == v) & (met.metric == 'r_all')][col].iloc[0]
    panels = [(f'ORN input to the model\n(no step; r = {r_of("orn_input"):.2f})', 'spike delta (Hz)', V['orn_input'].loc[ord_g]),
              (f'ORNs measured (pebbled)\n(r = {r_of("orn_dff"):.2f})', 'dF/F', V['orn_dff'].loc[ord_g]),
              (f'Model PN, olfsysm default\n(r = {r_of("olfsysm_default"):.2f})', 'evoked rate (Hz)', V['olfsysm_default'].loc[ord_g]),
              (f'Best step: {labels[best].split(" (")[0]}\n(r = {r_of(best):.2f}; {r_of(best, "oos"):.2f} out of sample)',
               'evoked rate (Hz)' if best.startswith('olfsysm') or best == 'kennedy_ss' else 'response (a.u.)',
               V[best].loc[ord_g]),
              ('GH146 PN dendrites\n(7 flies, fly mean)', 'dF/F', G.loc[ord_g])]
    note = ('28 glomeruli measured in GH146, the pebbled ORNs and the model input; 12 megamat odors of the APL panel. '
            'r: Pearson r with GH146 over all 336 entries. Colors: each panel on its own scale, 0 at gray. '
            'Glomeruli ordered by the odor of their largest GH146 response.')
    heatmaps(panels, out / 'heatmaps.png', note)
    tot = pd.DataFrame({'GH146 PN dendrites': G.sum(), 'ORN input (no step)': V['orn_input'].sum(),
                        'model PN, olfsysm default': V['olfsysm_default'].sum(), 'best fit': V[best].sum()})
    tot = tot / tot.mean()
    tot = tot.loc[tot['GH146 PN dendrites'].sort_values(ascending=False).index]
    tot.rename_axis('odor').to_csv(out / 'totals.csv')
    fig_order = ['orn_input', 'orn_dff', 'olfsysm_default'] + [f'olfsysm_inhadd{a:g}' for a in LADDER] + \
                ['kennedy_ss', 'olsen_fit', 'orn_power_fit', 'olfsysm_fit']
    metrics_figure(met, ceil, tot, fig_order, labels, {'olfsysm_default': '#eb6834', best: '#1baf7a', 'orn_input': '#2a78d6'},
                   out / 'metrics.png')

    prov = H.model_provenance()
    prov.update(script=H.script_path(__file__), generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), seconds=time.time() - t_start,
                inputs=dict(gh146=str(E90.gh146_dir() / 'GH146_ij_certain-roi_stats.p'),
                            pebbled=str(E90.gh146_dir() / 'pebbled_ij_certain-roi_stats.p'),
                            orn_input='kcq.harness.apl16_orn_deltas()', e32=str(H.shared_dir() / 'results' / 'E32_pn_stage')),
                glomeruli=glom, odors=ODORS, best_by_r_all=best, n_boot=N_BOOT, boot_seed=BOOT_SEED,
                windows=dict(odor=ODOR_WINDOW, pre=PRE_WINDOW))
    (out / 'provenance.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        piv = met[met.metric.isin(R_METRICS + ('cv_totals',))].pivot(index='version', columns='metric', values='value')
        print(piv.loc[[v for v in fig_order] + ['gh146']].round(3).to_string())
        print(ceil.round(3).to_string())
    print(f'best by r_all: {best}; wrote {out} ({time.time() - t_start:.0f} s)')


if __name__ == '__main__':
    main()
