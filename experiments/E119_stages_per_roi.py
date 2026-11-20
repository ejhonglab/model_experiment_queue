#!/usr/bin/env python3
"""E119: the model's boutons, claws and claw-to-soma step against Remy's per-ROI megamat recordings (no simulations).

The entry is item E119 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-20 13:32). Part 0 is data_prep/remy_per_roi_data.py: Remy's per-ROI files reduced to ROI x trial window means
and each ROI's pre-odor noise SD (<shared>/data/remy_per_roi/), 5 bouton flies (2022-10-21 fly02 excluded, as Remy's
dashboards leave it out) and 4 claw flies, the 12 megamat odors of the APL16 panel (66 pairs).

The model side is E90's archive (data/model_outputs/E90_model_stages/B/: C2's family, the breadth-2 APL-effect fit of
E67, on E72 part B's eight wirings, APL on and acutely off; window 0 to 0.5 s) with the pre-odor means at the settled
level (E90_settled_preodor's preodor_levels.parquet, -0.25 to 0 s). E90's readings of a stage are reused: the change
over the window (`stage_units`), the claw signal with a share of the KC's spike count (`claw_signal`, CLAW_SPIKE_SHARES),
the session z-score (`session_z`), its groups and its rule (`differs`, `element_result`, `group_verdict`). E81's noise
removal from products of different repeats (`H.noise_free_pair_distances`; Pearson on trials centred over the units, as
E92's `noise_removed`).

The model's units are pooled into ROIs before any per-ROI statistic (CLAW_POOLS, BOUTON_POOLS): claws in groups of 14
of one glomerulus carry the verdict ('g14'): a claw ROI averages about 14 claws of different KCs (Remy's notes), and
those claws surround one PN bouton (a microglomerulus), so they share their PN's input; groups of 7 and 28, and random
groups ('r<n>'), beside. Boutons are not pooled for the verdict (Remy's bouton ROIs outnumber the model's 389 boutons
2 to 3 times; groups of 3 beside). A ROI's value is the mean of its
units' values, read afterwards (session z per ROI, then measurement noise).

Readings of a model ROI (E90's): 'rate' (the change itself), 'session_z' and 'noisy' (session z through measurement noise
at the flies' per-odor trial-to-trial reliability, three trials, as E90's `noisy_reading`; its distances with the noise
removed by E81's estimator, its per-ROI statistics on the trial mean); claws also at each CLAW_SPIKE_SHARES (E90's nine
claw readings).

Per-ROI statistics (`roi_stats`), at a criterion: the flies' trial means in baseline-noise units (the window's change
over the ROI's pre-odor frame SD) above c SDs, c in CRITERIA; the model's at the threshold that gives the flies' mean
overall response rate at c (`matched_threshold`), so the overall rate is matched by construction and the comparison is of
how the responses are spread: the per-odor response rates, the number of odors per ROI (bins 0, 1, 2-3, 4-7, 8-12 of the
12), the responders' coefficient of variation, the top-1 % share, the share of entries below -threshold (inhibited) and
the per-odor totals' coefficient of variation. A wiring lies within the flies' spread on a scalar if it is inside the
flies' range, on a profile (per-odor rates, odors per ROI) if its distance from the flies' mean profile (RMS, total
variation) is at most the most atypical fly's from the other flies' mean (`profile_band`). Every "differs" follows E90's
rule: on distances `E90.differs` (the population VERDICT_DIFF or more from the flies' mean, outside their 90 % interval
over resampled flies, 6 of 8 wirings on that side); on per-ROI statistics 6 of 8 wirings outside the flies' spread on
the same side.

Parts (CLI `--part`):
- A, the boutons: A1 per-odor totals (a1_totals.csv), beside E114's exponent; A2 per-ROI statistics at the three
  criteria (a2_roi_stats.csv, a2_within.csv); A3 distances per fly, cosine and Pearson, trial means and noise removed,
  E90's groups, and the bouton element of E90's verdict read on them (a3_flies.csv, a3_model.csv, a3_element.csv);
  A4 (a lead) bouton groups by tuning and by co-variation across repeats, each tuning group matched to the nearest model
  glomerulus (a4_groups.csv, a4_summary.csv).
- B, the claws: B1 E90's nine claw readings scored on the per-ROI statistics, kept if inside the flies' spread on 6 of
  8 wirings (b1_readings.csv); B2 the claw element on the per-ROI flies under the kept readings, and E90's group
  verdicts read again under the kept readings only (b2_element.csv, b2_verdicts.csv); B3 the bouton-to-claw change
  (b3_change.csv).
- C, the same fly's claw-to-soma step (part 0b's trials; `--part C`): C1 the step per pair from the same trials (c1_step.csv),
  C2 the somata's response rate against the claws' summed response (c2_curve.csv).
- table: comparison_table.md and outcome.json.

Outputs: <shared>/results/E119_stages_per_roi/ (the CSVs, comparison_table.md, outcome.json, provenance_<part>.json);
caches under outputs/E119_stages_per_roi/.

Usage (from the repository root; light, no simulations):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E119_stages_per_roi.py --part A
    ... --part B; --part C; --part table
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E90_model_stages as E90
from kcq import harness as H

EXPERIMENT = 'E119_stages_per_roi'
ODORS = tuple(E90.ODORS)
GROUPS = E90.GROUPS
STAGES = ('bouton', 'claw')
WINDOW = 'e90'                          # part 0's window read: E90's STAGE_WINDOW, 0.5 to 3 s
NOISE_VAR = {'bouton': 'F', 'claw': 'Fc'}           # read in baseline-noise units
DIST_VAR = {'bouton': 'F_zscore', 'claw': 'Fc_zscore'}   # Remy's variables, for the distances
CRITERIA = (1.0, 2.0, 3.0)              # baseline-noise SDs
VERDICT_CRITERION = 2.0
# pooling of the model's units into ROIs, as labels: '1' (none), 'r<n>' (random groups of n), 'g<n>' (groups of n units of
# one glomerulus: a claw ROI surrounds one PN bouton, a microglomerulus, so its claws share their PN's input)
CLAW_POOLS = ('1', 'g7', 'g14', 'g28', 'r7', 'r14', 'r28')
BOUTON_POOLS = ('1', 'g3', 'r3')
VERDICT_POOL = {'bouton': '1', 'claw': 'g14'}
POOL_SEED = 0
NOISE_SEED = 0
READINGS = ('rate', 'session_z', 'noisy')
SHARES = {'bouton': (0.0,), 'claw': tuple(E90.CLAW_SPIKE_SHARES)}
SUBSET = {'bouton': 'shared', 'claw': 'all'}         # E90's VERDICT_SUBSETS
CONDS = (('on', 'on'), ('off', 'acute_off'))
WIRINGS = tuple(E90.ALL_WIRINGS)
MIN_SAME = E90.VERDICT_SAME_SIDE        # 6 of 8
OPR_BINS = ((0, 0), (1, 1), (2, 3), (4, 7), (8, 12))
OPR_LABELS = ('0', '1', '2-3', '4-7', '8-12')
MEASURES = ('cosine', 'pearson')
N_BOOT = E90.N_BOOT
A4_TUNING_CUT = 0.3                     # A4: average-linkage cut, 1 - Pearson r of the 12-odor tuning
A4_COVAR_CUT = 0.7                      # A4: average-linkage cut, 1 - Pearson r of the repeat residuals
A4_MIN_ROIS = 5
A4_MATCH_R = 0.5                        # A4: a tuning group is unmatched if no model glomerulus reaches this r ...
A4_NULL_N = 200                         # ... or its best r is not above the 95th percentile of the best r over the
                                        # glomeruli with the group's odors permuted (the chance level of a best of 54)
ARCHIVE = ('data', 'model_outputs', 'E90_model_stages', 'B')
LEVELS = ('results', 'E90_model_stages', 'settled_preodor', 'preodor_levels.parquet')
LEVEL_TOL = 1e-9


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir() -> Path:
    return H.outputs_dir(EXPERIMENT)


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def provenance(part: str, **settings) -> None:
    _write_json(out_dir() / f'provenance_{part}.json', dict(part=part, script=H.script_path(__file__), **H.code_commits(),
        generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), settings=settings))


# ---------------------------------------------------------------------------
# the flies (part 0's tables)

def roi_dir() -> Path:
    return H.shared_dir() / 'data' / 'remy_per_roi'


@lru_cache(maxsize=None)
def _part0() -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = roi_dir()
    lock = json.loads((d / 'lock.json').read_text())
    if not lock['checks']['passed']:
        raise RuntimeError("part 0's checks did not pass")
    return pd.read_parquet(d / 'window_means.parquet'), pd.read_parquet(d / 'noise_sd.parquet')


def trials_array(wm: pd.DataFrame, col: str, odors: Sequence[str] = ODORS) -> Tuple[np.ndarray, np.ndarray]:
    """(repeats x ROIs x odors) of one acq's column `col`, and the ROI ids."""
    sub = wm[wm.odor.isin(list(odors))]
    p = sub.pivot_table(index=['repeat', 'roi'], columns='odor', values=col, aggfunc='first')
    rois = np.asarray(sorted(sub.roi.unique()))
    reps = sorted(sub.repeat.unique())
    X = np.stack([p.loc[r].reindex(rois)[list(odors)].to_numpy(dtype=float) for r in reps])
    return X, rois


def fly_trials(stage: str, odors: Sequence[str] = ODORS) -> Dict[str, Dict[str, np.ndarray]]:
    """Per fly (acq) of `stage`: 'noise' (repeats x ROIs x odors: the window's change, window mean minus the pre-odor
    mean of NOISE_VAR, over the ROI's pre-odor noise SD), 'z' (the same change of Remy's z variable, DIST_VAR), 'z_raw'
    (the window mean of DIST_VAR as Remy reads it, no pre-odor mean removed), 'rois'. ROIs with a zero noise SD are left
    out."""
    wm, ns = _part0()
    out = {}
    for acq, sub in wm[wm.stage == stage].groupby('acq'):
        v, z = NOISE_VAR[stage], DIST_VAR[stage]
        Xw, rois = trials_array(sub, f'{v}__{WINDOW}', odors)
        Xp, _ = trials_array(sub, f'{v}__pre', odors)
        Zw, _ = trials_array(sub, f'{z}__{WINDOW}', odors)
        Zp, _ = trials_array(sub, f'{z}__pre', odors)
        sd = ns[(ns.stage == stage) & (ns.acq == acq)].set_index('roi')[v].reindex(rois).to_numpy(dtype=float)
        keep = np.isfinite(sd) & (sd > 0)
        out[acq] = dict(noise=((Xw - Xp) / sd[None, :, None])[:, keep], z=(Zw - Zp)[:, keep], z_raw=Zw[:, keep],
                        rois=rois[keep])
    return out


def fly_reliability(stage: str) -> pd.Series:
    """The flies' per-odor trial-to-trial reliability: per fly and odor the mean over the pairs of different repeats of
    the Pearson r across ROIs (z change), averaged over flies (the noise level of the model's 'noisy' reading)."""
    rows = []
    for acq, d in fly_trials(stage).items():
        X = d['z']
        R = X.shape[0]
        r = []
        for k in range(X.shape[2]):
            r.append(np.mean([np.corrcoef(X[i, :, k], X[j, :, k])[0, 1] for i in range(R) for j in range(i + 1, R)]))
        rows.append(r)
    return pd.Series(np.mean(rows, axis=0), index=list(ODORS))


# ---------------------------------------------------------------------------
# per-ROI statistics

def opr_distribution(n_odors: np.ndarray) -> np.ndarray:
    """The fraction of ROIs in each OPR_BINS bin of the number of odors they respond to."""
    n = np.asarray(n_odors)
    return np.array([np.mean((n >= lo) & (n <= hi)) for lo, hi in OPR_BINS])


def top_share(V: np.ndarray, frac: float = 0.01) -> float:
    """Per odor, the share of the summed positive response carried by the top `frac` of ROIs (at least one), averaged
    over odors (NaN odors with no positive response left out)."""
    V = np.asarray(V, dtype=float)
    k = max(1, int(np.ceil(frac * V.shape[0])))
    out = []
    for j in range(V.shape[1]):
        p = np.clip(V[:, j], 0, None)
        s = p.sum()
        if s > 0:
            out.append(np.sort(p)[::-1][:k].sum() / s)
    return float(np.mean(out)) if out else np.nan


def roi_stats(V: pd.DataFrame, thr: float) -> dict:
    """The per-ROI statistics of ROIs x odors values V at threshold `thr`: 'rate' (overall response rate), 'per_odor'
    (Series), 'opr' (odors-per-ROI distribution, OPR_BINS), 'opr_mean', 'inhibited' (share of entries below -thr),
    'resp_cv' (the responders' coefficient of variation), 'top1' (`top_share`), 'totals_cv' (the coefficient of
    variation over odors of the mean response per odor, all ROIs, signed), 'n_rois'."""
    X = V[list(ODORS)].to_numpy(dtype=float)
    R = X > thr
    resp = X[R]
    tot = X.mean(axis=0)
    return dict(rate=float(R.mean()), per_odor=pd.Series(R.mean(axis=0), index=list(ODORS)),
                opr=opr_distribution(R.sum(axis=1)), opr_mean=float(R.sum(axis=1).mean()),
                inhibited=float((X < -thr).mean()),
                resp_cv=float(resp.std() / resp.mean()) if resp.size > 1 and resp.mean() > 0 else np.nan,
                top1=top_share(X), totals_cv=float(tot.std() / tot.mean()) if tot.mean() != 0 else np.nan,
                n_rois=int(X.shape[0]))


def matched_threshold(V: pd.DataFrame, rate: float) -> float:
    """The threshold at which the overall response rate of V (ROIs x ODORS entries above it) is `rate`."""
    X = V[list(ODORS)].to_numpy(dtype=float).ravel()
    return float(np.quantile(X, 1.0 - rate))


SCALARS = ('rate', 'opr_mean', 'inhibited', 'resp_cv', 'top1', 'totals_cv')


def profile_distance(kind: str, a, b) -> float:
    """'per_odor': RMS difference of two per-odor rate profiles; 'opr': total variation distance of two odors-per-ROI
    distributions."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean((a - b) ** 2))) if kind == 'per_odor' else float(0.5 * np.abs(a - b).sum())


def profile_band(profiles: Sequence[np.ndarray], kind: str) -> Tuple[np.ndarray, float]:
    """The flies' mean profile and the band: the largest distance of a fly's profile from the mean of the other flies'."""
    P = np.asarray([np.asarray(p, dtype=float) for p in profiles])
    band = max(profile_distance(kind, P[i], np.delete(P, i, axis=0).mean(axis=0)) for i in range(len(P)))
    return P.mean(axis=0), float(band)


def fly_roi_stats(stage: str) -> Dict[float, Dict[str, dict]]:
    """{criterion: {acq: roi_stats}} of the flies' trial means in baseline-noise units."""
    ft = fly_trials(stage)
    out = {}
    for c in CRITERIA:
        out[c] = {a: roi_stats(pd.DataFrame(d['noise'].mean(axis=0), columns=list(ODORS)), c) for a, d in ft.items()}
    return out


def fly_spread(stats: Dict[str, dict]) -> dict:
    """The flies' spread of the per-ROI statistics: per scalar (min, max, mean); per profile (mean, band)."""
    out = {}
    for s in SCALARS:
        v = np.array([x[s] for x in stats.values()], dtype=float)
        out[s] = dict(lo=float(np.nanmin(v)), hi=float(np.nanmax(v)), mean=float(np.nanmean(v)))
    for k in ('per_odor', 'opr'):
        m, b = profile_band([np.asarray(x[k], dtype=float) for x in stats.values()], k)
        out[k] = dict(mean=m, band=b)
    return out


def within(stats: dict, spread: dict) -> Dict[str, object]:
    """Where one wiring's statistics lie against the flies' spread: per scalar the side (-1 below the range, 0 inside,
    +1 above), per profile its distance and whether inside the band."""
    out = {}
    for s in SCALARS:
        v = stats[s]
        out[f'{s}_side'] = 0 if not np.isfinite(v) else (-1 if v < spread[s]['lo'] else (1 if v > spread[s]['hi'] else 0))
    for k in ('per_odor', 'opr'):
        d = profile_distance(k, np.asarray(stats[k], dtype=float), spread[k]['mean'])
        out[f'{k}_dist'] = d
        out[f'{k}_in'] = bool(d <= spread[k]['band'])
    return out


def scalar_differs(sides: Sequence[int], min_same: int = MIN_SAME) -> Tuple[bool, int]:
    """The per-ROI rule: at least `min_same` wirings outside the flies' range on the same side; returns (differs,
    the side)."""
    s = np.asarray(sides)
    for side in (-1, 1):
        if (s == side).sum() >= min_same:
            return True, side
    return False, 0


# ---------------------------------------------------------------------------
# the model

def archive_dir() -> Path:
    return H.shared_dir().joinpath(*ARCHIVE)


@lru_cache(maxsize=None)
def _levels() -> pd.DataFrame:
    return pd.read_parquet(H.shared_dir().joinpath(*LEVELS))


@lru_cache(maxsize=None)
def unit_meta(wiring: str) -> Dict[str, np.ndarray]:
    """E90's units' metadata of a wiring (bouton_glomerulus, claw_glomerulus, kc_of_claw), from the archive."""
    d = json.loads((archive_dir() / f'{wiring}_units.json').read_text())
    return {k: np.asarray(v) for k, v in d.items()}


def model_stage(wiring: str, cc: str) -> dict:
    """E90's saved stage run of a wiring and condition ('on' / 'acute_off') with every unit's pre-odor mean at the
    settled level. Checks that the levels' E90-window column equals the saved pre-odor means."""
    st = E90.load_stage('B', f'{wiring}_{cc}')
    if st is None:
        raise RuntimeError(f"E90's stage run B {wiring}_{cc} not found")
    lv = _levels()
    lv = lv[(lv.part == 'B') & (lv.tag == wiring) & (lv.cond == cc)]
    st = dict(st)
    for k in ('pn', 'bouton', 'claw'):
        x = st[f'{k}_pre']
        l = lv[lv.stage == k].sort_values('unit')
        if len(l) != x.shape[0]:
            raise RuntimeError(f'{wiring} {cc} {k}: {len(l)} levels for {x.shape[0]} units')
        if np.max(np.abs(l.e90.to_numpy() - x.iloc[:, 0].to_numpy())) > 1e-6 * max(1.0, float(np.abs(x.to_numpy()).max())):
            raise RuntimeError(f"{wiring} {cc} {k}: the levels' E90 window is not the saved pre-odor mean")
        st[f'{k}_pre'] = pd.DataFrame(np.repeat(l.settled.to_numpy()[:, None], x.shape[1], axis=1), index=x.index,
                                      columns=x.columns)
    return st


@lru_cache(maxsize=None)
def model_units(stage: str, wiring: str, share: float) -> Dict[str, pd.DataFrame]:
    """A wiring's stage as units x the panel's odors (E90's `stage_units`, the settled pre-odor level), APL on and off,
    restricted to E90's verdict subset (SUBSET); claws as E90's `claw_signal` at `share` (scales from the APL-on run)."""
    st = {c: model_stage(wiring, cc) for c, cc in CONDS}
    meta = unit_meta(wiring)
    out = {}
    if stage == 'claw':
        ch_on = E90.stage_units(st['on'], 'claw')
        scales = (E90._rms(ch_on.to_numpy()), E90._rms(st['on']['spike_counts'].to_numpy(dtype=float)))
    for c in ('on', 'off'):
        U = E90.stage_units(st[c], stage)
        m = E90.subset_mask(stage, U.index, meta, SUBSET[stage])
        if stage == 'claw':
            U = E90.claw_signal(U, st[c]['spike_counts'], np.asarray(meta['kc_of_claw']), share, scales)
        out[c] = U.loc[m]
    gl = meta['claw_glomerulus'] if stage == 'claw' else meta['bouton_glomerulus']
    out['groups'] = np.asarray(gl)[E90.subset_mask(stage, E90.stage_units(st['on'], stage).index, meta, SUBSET[stage])]
    return out


def pool_units(U: pd.DataFrame, pool: str, seed: int, groups: Optional[np.ndarray] = None) -> pd.DataFrame:
    """Units pooled into ROIs, each ROI the mean of its units: `pool` '1' returns U; 'r<n>' a random partition
    (generator `seed`) into groups of n (the remainder dropped); 'g<n>' within each value of `groups` (the units'
    glomeruli) a random partition into groups of n, a remainder of at least n / 2 units kept as a smaller ROI."""
    pool = str(pool)
    if pool == '1':
        return U.reset_index(drop=True)
    kind, size = pool[0], int(pool[1:])
    rng = np.random.default_rng(seed)
    X = U.to_numpy(dtype=float)
    if kind == 'r':
        idx = rng.permutation(U.shape[0])
        n = U.shape[0] // size
        Y = X[idx[:n * size]].reshape(n, size, U.shape[1]).mean(axis=1)
    elif kind == 'g':
        if groups is None or len(groups) != U.shape[0]:
            raise ValueError("pooling by glomerulus needs the units' glomeruli")
        rows = []
        g = np.asarray(groups)
        for k in sorted(set(g.tolist())):
            idx = rng.permutation(np.flatnonzero(g == k))
            for i in range(0, len(idx), size):
                chunk = idx[i:i + size]
                if len(chunk) == size or len(chunk) >= size / 2:
                    rows.append(X[chunk].mean(axis=0))
        Y = np.asarray(rows)
    else:
        raise ValueError(pool)
    return pd.DataFrame(Y, columns=U.columns)


def noisy_trials(Z: pd.DataFrame, reliability: pd.Series, rng: np.random.Generator, n_trials: int = 3) -> np.ndarray:
    """Three trials (trials x ROIs x ODORS) of session-z values Z through measurement noise at the flies' per-odor
    reliability r (sigma^2 = variance over ROIs x (1 - r) / r, E90's `noisy_reading`)."""
    X = Z[list(ODORS)].to_numpy(dtype=float)
    r = reliability.reindex(list(ODORS)).to_numpy(dtype=float)
    sig = np.sqrt(X.var(axis=0) * (1.0 - r) / r)
    return X[None] + sig[None, None, :] * rng.standard_normal((n_trials,) + X.shape)


def model_rois(stage: str, wiring: str, cond: str, reading: str, share: float, pool: str,
               reliability: Optional[pd.Series] = None) -> np.ndarray:
    """A wiring's ROIs in a reading: trials x ROIs x ODORS (one 'trial' for the noise-free readings). Pooling seed
    POOL_SEED + the wiring's index; noise seed NOISE_SEED + the wiring's index (the same for both conditions)."""
    wi = WIRINGS.index(wiring)
    mu = model_units(stage, wiring, float(share))
    return rois_of_units(mu[cond], mu['groups'], stage, reading, pool, wi, reliability)


def rois_of_units(U: pd.DataFrame, groups: Optional[np.ndarray], stage: str, reading: str, pool: str, wi: int,
                  reliability: Optional[pd.Series] = None) -> np.ndarray:
    """`model_rois` on given units x ODORS values U (and their glomeruli `groups`, for 'g<n>' pooling): pooled with seed
    POOL_SEED + wi, silent ROIs dropped, read in `reading` (noise seed NOISE_SEED + wi); trials x ROIs x ODORS. E120
    reads other runs' stages with it."""
    V = E90._drop_silent(pool_units(U, pool, POOL_SEED + wi, groups))
    if reading == 'rate':
        return V[list(ODORS)].to_numpy(dtype=float)[None]
    Z = E90.session_z(V)
    if reading == 'session_z':
        return Z[list(ODORS)].to_numpy(dtype=float)[None]
    if reading == 'noisy':
        rel = fly_reliability(stage) if reliability is None else reliability
        return noisy_trials(Z, rel, np.random.default_rng(NOISE_SEED + wi))
    raise ValueError(reading)


# ---------------------------------------------------------------------------
# distances

def centred(trials: np.ndarray) -> np.ndarray:
    x = np.asarray(trials, dtype=float)
    return x - x.mean(axis=1, keepdims=True)


def trial_mean_distances(trials: np.ndarray, measure: str) -> pd.Series:
    """The 66 pair distances of the trial means (trials x units x ODORS), 1 - cosine or 1 - Pearson across units."""
    X = pd.DataFrame(np.asarray(trials, dtype=float).mean(axis=0), columns=list(ODORS))
    return E90.odor_distances(X, measure)


def noise_removed(trials: np.ndarray, measure: str) -> pd.Series:
    """E81's noise-removed 66 pair distances of trials (repeats x units x ODORS), Pearson on centred trials (E92)."""
    x = centred(trials) if measure == 'pearson' else np.asarray(trials, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        d = H.noise_free_pair_distances(x, list(ODORS))
    return E90._cs(d).reindex(E90.megamat_pairs())


def model_distances(trials: np.ndarray, measure: str) -> pd.Series:
    """A model reading's distances: noise-free readings (one trial) directly, the noisy reading with the noise removed."""
    return trial_mean_distances(trials, measure) if trials.shape[0] == 1 else noise_removed(trials, measure)


def fly_distances(stage: str) -> pd.DataFrame:
    """Per fly, reading ('trial_mean', 'noise_removed', 'remy_raw': Remy's trial-mean Pearson on the window means
    without the pre-odor mean removed) and measure: the 66 pair distances (rows acq, reading, measure; columns pairs)."""
    rows = {}
    for acq, d in fly_trials(stage).items():
        for ms in MEASURES:
            rows[(acq, 'trial_mean', ms)] = trial_mean_distances(d['z'], ms)
            rows[(acq, 'noise_removed', ms)] = noise_removed(d['z'], ms)
        rows[(acq, 'remy_raw', 'pearson')] = trial_mean_distances(d['z_raw'], 'pearson')
    df = pd.DataFrame(rows).T
    df.index.names = ['acq', 'reading', 'measure']
    return df


def group_table_flies(fd: pd.DataFrame, reading: str, measure: str) -> pd.DataFrame:
    """The flies' group means of one reading and measure: per group the mean over flies, 90 % interval over resampled
    flies (E90's `bootstrap_mean`), and each fly's value."""
    S = fd.xs((reading, measure), level=('reading', 'measure'))
    pairs = E90.megamat_pairs()
    boot = E90.bootstrap_mean(S[pairs])
    gb = E90._group_boot(boot, pairs)
    mean = S[pairs].mean(axis=0)
    gm = E90.group_means(mean)
    rows = []
    for g in GROUPS:
        lo, hi = np.nanpercentile(gb[g], E90.CI)
        per = {f'fly_{a}': E90.group_means(S.loc[a, pairs])[g] for a in S.index}
        rows.append(dict(group=g, mean=gm[g], lo=float(lo), hi=float(hi), **per))
    return pd.DataFrame(rows).set_index('group')


def stage_element(stage: str, readings: Sequence[Tuple[str, float]], pool: str, fd: pd.DataFrame,
                  fly_reading: str = 'noise_removed') -> pd.DataFrame:
    """E90's element of a stage read on the per-ROI flies: per reading (name, share), measure and group, the model's
    population (mean over wirings of the group means, APL on) against the flies' (`group_table_flies`) by E90's
    `differs`. Rows: reading, share, measure, group, model, flies, lo, hi, diff, n_same_side, differs, and each wiring."""
    rows = []
    rel = fly_reliability(stage)
    for rd, sh in readings:
        for ms in MEASURES:
            fl = group_table_flies(fd, fly_reading, ms)
            per_w = {}
            for w in WIRINGS:
                per_w[w] = E90.group_means(model_distances(model_rois(stage, w, 'on', rd, sh, pool, rel), ms))
            for g in GROUPS:
                vals = [per_w[w][g] for w in WIRINGS]
                pop = float(np.mean(vals))
                r = E90.differs(pop, vals, fl.loc[g, 'mean'], fl.loc[g, 'lo'], fl.loc[g, 'hi'])
                rows.append(dict(stage=stage, reading=rd, share=sh, pool=pool, measure=ms, group=g, model=pop,
                                 flies=fl.loc[g, 'mean'], lo=fl.loc[g, 'lo'], hi=fl.loc[g, 'hi'], **r,
                                 **{f'w_{w}': v for w, v in zip(WIRINGS, vals)}))
    return pd.DataFrame(rows)


def element_by_group(el: pd.DataFrame) -> pd.DataFrame:
    """E90's `element_result` per group over the rows (readings x measures)."""
    rows = []
    for g in GROUPS:
        sub = el[el.group == g]
        rows.append(dict(group=g, **E90.element_result(sub)))
    return pd.DataFrame(rows).set_index('group')


# ---------------------------------------------------------------------------
# per-ROI tables (A2, B1)

def roi_table(stage: str, pools: Sequence[str], conds: Sequence[str] = ('on', 'off')) -> pd.DataFrame:
    """Every wiring's per-ROI statistics at the flies' matched rates: rows stage, criterion, reading, share, pool, cond,
    wiring, threshold, the scalars, the profile distances and within-ness (`within`), and the flies' spread; the flies'
    own rows (wiring 'fly:<acq>')."""
    fs = fly_roi_stats(stage)
    rel = fly_reliability(stage)
    rows = []
    for c in CRITERIA:
        sp = fly_spread(fs[c])
        rate = sp['rate']['mean']
        for acq, st in fs[c].items():
            rows.append(dict(stage=stage, criterion=c, reading='flies', share=np.nan, pool='flies', cond='flies',
                             wiring=f'fly:{acq}', threshold=c, **{s: st[s] for s in SCALARS},
                             opr=json.dumps([round(float(x), 6) for x in st['opr']]), **within(st, sp)))
        for rd in READINGS:
            for sh in SHARES[stage]:
                for pool in pools:
                    for cond in conds:
                        for w in WIRINGS:
                            T = model_rois(stage, w, cond, rd, sh, pool, rel)
                            V = pd.DataFrame(T.mean(axis=0), columns=list(ODORS))
                            # the threshold is matched on the APL-on ROIs and used for both conditions
                            Von = V if cond == 'on' else pd.DataFrame(
                                model_rois(stage, w, 'on', rd, sh, pool, rel).mean(axis=0), columns=list(ODORS))
                            thr = matched_threshold(Von, rate)
                            st = roi_stats(V, thr)
                            rows.append(dict(stage=stage, criterion=c, reading=rd, share=sh, pool=pool, cond=cond, wiring=w,
                                             threshold=thr, **{s: st[s] for s in SCALARS},
                                             opr=json.dumps([round(float(x), 6) for x in st['opr']]), **within(st, sp)))
    df = pd.DataFrame(rows)
    sp_rows = []
    for c in CRITERIA:
        sp = fly_spread(fs[c])
        sp_rows.append(dict(stage=stage, criterion=c, **{f'{s}_{k}': sp[s][k] for s in SCALARS for k in ('lo', 'hi', 'mean')},
                            per_odor_band=sp['per_odor']['band'], opr_band=sp['opr']['band'],
                            opr_mean_profile=json.dumps([round(float(x), 6) for x in sp['opr']['mean']]),
                            per_odor_mean=json.dumps({o: round(float(v), 6) for o, v in zip(ODORS, sp['per_odor']['mean'])})))
    df.attrs['spread'] = pd.DataFrame(sp_rows)
    return df


def within_summary(rt: pd.DataFrame) -> pd.DataFrame:
    """Per stage, criterion, reading, share, pool and condition: per scalar the number of wirings below / inside / above
    the flies' range and whether it differs (`scalar_differs`), per profile the number inside the band."""
    m = rt[rt.cond != 'flies']
    rows = []
    for k, sub in m.groupby(['stage', 'criterion', 'reading', 'share', 'pool', 'cond']):
        r = dict(zip(['stage', 'criterion', 'reading', 'share', 'pool', 'cond'], k), n_wirings=len(sub))
        for s in SCALARS:
            sides = sub[s + '_side']
            d, side = scalar_differs(sides.to_numpy())
            r.update({f'{s}_below': int((sides == -1).sum()), f'{s}_above': int((sides == 1).sum()), f'{s}_differs': d,
                      f'{s}_side': side})
        for p in ('per_odor', 'opr'):
            r[f'{p}_n_in'] = int(sub[f'{p}_in'].sum())
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# part A

def a1_totals(rt: pd.DataFrame) -> pd.DataFrame:
    """A1: the per-odor totals' coefficient of variation over the 12 odors, flies and model (APL on), with E90's PN
    stage CV and E114's exponent beside (a1_totals.csv)."""
    sub = rt[(rt.stage == 'bouton') & (rt.criterion == VERDICT_CRITERION)]
    rows = []
    for k, s in sub.groupby(['reading', 'share', 'pool', 'cond'], dropna=False):
        rows.append(dict(zip(['reading', 'share', 'pool', 'cond'], k), mean=float(s.totals_cv.mean()),
                         lo=float(s.totals_cv.min()), hi=float(s.totals_cv.max()), n=len(s)))
    return pd.DataFrame(rows)


def e114_exponent() -> Optional[dict]:
    p = H.shared_dir() / 'results' / 'E114_input_compression' / 'choice.json'
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    o = H.shared_dir() / 'results' / 'E114_input_compression' / 'outcome.json'
    return dict(exponent=d.get('exponent'), reason=d.get('reason'),
                outcome=json.loads(o.read_text()) if o.exists() else None)


def a4_groups(link_cut_tuning: float = A4_TUNING_CUT, link_cut_covar: float = A4_COVAR_CUT) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A4 (a lead): per bouton fly, the responsive ROIs (any odor above VERDICT_CRITERION noise SDs) grouped by their
    12-odor tuning (trial means, baseline-noise units; average linkage on 1 - Pearson r, cut at `link_cut_tuning`) and by
    their co-variation across repeats (each trial minus its odor's mean; cut at `link_cut_covar`); every tuning group of
    A4_MIN_ROIS or more ROIs matched to the model glomerulus (real wiring, APL on, mean bouton change of its boutons, all
    54) whose tuning has the highest Pearson r. Returns (groups, summary)."""
    from scipy.cluster.hierarchy import fcluster, linkage
    st = model_stage('real', 'on')
    U = E90.stage_units(st, 'bouton')
    gl = unit_meta('real')['bouton_glomerulus']
    G = U[list(ODORS)].groupby(np.asarray(gl)).mean()
    G = G[G.std(axis=1) > 0]
    inputs = set(H.apl16_orn_deltas().index)
    rows, summ = [], []
    for acq, d in fly_trials('bouton').items():
        X = d['noise']
        M = X.mean(axis=0)
        keep = (M > VERDICT_CRITERION).any(axis=1)
        Mk, Xk = M[keep], X[:, keep]
        if Mk.shape[0] < 2:
            continue
        lt = fcluster(linkage(Mk, 'average', metric='correlation'), link_cut_tuning, 'distance')
        res = (Xk - Mk[None]).transpose(1, 0, 2).reshape(Mk.shape[0], -1)
        lc = fcluster(linkage(res, 'average', metric='correlation'), link_cut_covar, 'distance')
        n_unm = 0
        n_in = 0
        for g in np.unique(lt):
            idx = lt == g
            if idx.sum() < A4_MIN_ROIS:
                continue
            tun = Mk[idx].mean(axis=0)
            Gm = G.to_numpy(dtype=float)
            Gz = (Gm - Gm.mean(axis=1, keepdims=True)) / Gm.std(axis=1, keepdims=True)

            def best_r(v):
                vz = (v - v.mean()) / v.std()
                return Gz @ vz / len(v)
            r = pd.Series(best_r(tun), index=G.index)
            best = r.idxmax()
            prng = np.random.default_rng(int(g))
            null95 = float(np.percentile([best_r(prng.permutation(tun)).max() for _ in range(A4_NULL_N)], 95))
            cv = pd.Series(lc[idx]).value_counts()
            rows.append(dict(acq=acq, group=int(g), n_rois=int(idx.sum()), best_glomerulus=best, r_best=float(r.max()),
                             best_has_input=best in inputs, null95=null95,
                             matched=bool(r.max() >= A4_MATCH_R and r.max() > null95),
                             covar_groups=int(len(cv)), covar_top_share=float(cv.iloc[0] / idx.sum())))
            n_in += int(idx.sum())
            n_unm += int(idx.sum()) if not (r.max() >= A4_MATCH_R and r.max() > null95) else 0
        summ.append(dict(acq=acq, n_rois=int(M.shape[0]), n_responsive=int(keep.sum()), n_tuning_groups=int(len(np.unique(lt))),
                         n_big_groups=int(sum((lt == g).sum() >= A4_MIN_ROIS for g in np.unique(lt))),
                         n_covar_groups=int(len(np.unique(lc))), share_in_unmatched=n_unm / n_in if n_in else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(summ)


def part_a() -> dict:
    """Part A's tables (a1_totals.csv, a2_roi_stats.csv, a2_spread.csv, a2_within.csv, a3_flies.csv, a3_element.csv,
    a4_groups.csv, a4_summary.csv)."""
    t0 = time.time()
    od = out_dir()
    rt = roi_table('bouton', BOUTON_POOLS)
    rt = _side_cols(rt)
    rt.to_csv(od / 'a2_roi_stats.csv', index=False)
    rt.attrs['spread'].to_csv(od / 'a2_spread.csv', index=False)
    ws = within_summary(rt)
    ws.to_csv(od / 'a2_within.csv', index=False)
    a1 = a1_totals(rt)
    a1.to_csv(od / 'a1_totals.csv', index=False)
    fd = fly_distances('bouton')
    fd.to_csv(local_dir() / 'a3_fly_pairs.csv')
    ft = []
    for rd in ('trial_mean', 'noise_removed'):
        for ms in MEASURES:
            t = group_table_flies(fd, rd, ms).reset_index()
            t.insert(0, 'measure', ms)
            t.insert(0, 'reading', rd)
            ft.append(t)
    t = group_table_flies(fd, 'remy_raw', 'pearson').reset_index()
    t.insert(0, 'measure', 'pearson')
    t.insert(0, 'reading', 'remy_raw')
    ft.append(t)
    pd.concat(ft, ignore_index=True).to_csv(od / 'a3_flies.csv', index=False)
    el = pd.concat([stage_element('bouton', [(rd, 0.0) for rd in READINGS], p, fd) for p in BOUTON_POOLS], ignore_index=True)
    el.to_csv(od / 'a3_element.csv', index=False)
    g, s = a4_groups()
    g.to_csv(od / 'a4_groups.csv', index=False)
    s.to_csv(od / 'a4_summary.csv', index=False)
    provenance('A', CRITERIA=CRITERIA, BOUTON_POOLS=BOUTON_POOLS, VERDICT_POOL=VERDICT_POOL, WINDOW=WINDOW, SUBSET=SUBSET,
               POOL_SEED=POOL_SEED, NOISE_SEED=NOISE_SEED, A4=dict(tuning=A4_TUNING_CUT, covar=A4_COVAR_CUT,
               min_rois=A4_MIN_ROIS, match_r=A4_MATCH_R, null_n=A4_NULL_N), seconds=round(time.time() - t0, 1),
               e114=e114_exponent(), part0_lock=str(roi_dir() / 'lock.json'))
    print(f'part A in {time.time() - t0:.0f} s', flush=True)
    return dict(rt=rt, el=el)


def _side_cols(rt: pd.DataFrame) -> pd.DataFrame:
    """Recomputes each row's per-scalar side against the flies' spread of its criterion as '<scalar>_side' (the value
    columns keep the statistic itself)."""
    sp = rt.attrs['spread'].set_index('criterion')
    rt = rt.copy()
    for s in SCALARS:
        lo = rt.criterion.map(sp[f'{s}_lo'])
        hi = rt.criterion.map(sp[f'{s}_hi'])
        v = rt[s]
        rt[f'{s}_side'] = np.where(~np.isfinite(v), 0, np.where(v < lo, -1, np.where(v > hi, 1, 0)))
    rt.attrs['spread'] = sp.reset_index()
    return rt


# ---------------------------------------------------------------------------
# part B

def b1_readings(rt: pd.DataFrame) -> pd.DataFrame:
    """B1: E90's nine claw readings (READINGS x shares; APL on, VERDICT_POOL, VERDICT_CRITERION) scored on the per-ROI
    statistics: per wiring inside the flies' spread on all three (per-odor response rates within the band, odors per ROI
    within the band, the responders' CV inside the range); kept if on MIN_SAME of 8 wirings."""
    sub = rt[(rt.stage == 'claw') & (rt.cond == 'on') & (rt.criterion == VERDICT_CRITERION)]
    rows = []
    for (rd, sh, pool), s in sub.groupby(['reading', 'share', 'pool']):
        ok = s.per_odor_in & s.opr_in & (s.resp_cv_side == 0)
        rows.append(dict(reading=rd, share=sh, pool=pool, n_per_odor_in=int(s.per_odor_in.sum()),
                         n_opr_in=int(s.opr_in.sum()), n_resp_cv_in=int((s.resp_cv_side == 0).sum()), n_all_in=int(ok.sum()),
                         kept=bool(ok.sum() >= MIN_SAME), per_odor_dist=float(s.per_odor_dist.mean()),
                         opr_dist=float(s.opr_dist.mean()), resp_cv=float(s.resp_cv.mean()), opr_mean=float(s.opr_mean.mean()),
                         inhibited=float(s.inhibited.mean())))
    return pd.DataFrame(rows)


def e90_verdict_rows() -> pd.DataFrame:
    return pd.read_csv(H.shared_dir() / 'results' / 'E90_model_stages' / 'settled_preodor' / 'verdict_rows.csv')


def e90_verdicts_kept(kept: Sequence[Tuple[str, float]]) -> pd.DataFrame:
    """E90's group verdicts (settled pre-odor level, its primary subsets, its flies from the decoded matrices) read again
    with the claw and step elements restricted to the kept claw readings (E90's 'noisy_nr' is this script's 'noisy');
    the bouton element as E90 read it. Rows per group: verdict, each element's result, E90's own verdict."""
    v = e90_verdict_rows()
    v = v[v.verdict_subsets == 'primary']
    eff = json.loads((H.shared_dir() / 'results' / 'E90_model_stages' / 'settled_preodor' / 'verdict.json').read_text())
    name = {'noisy': 'noisy_nr'}
    keep = {(name.get(r, r), float(s)) for r, s in kept}
    rows = []
    for g in GROUPS:
        res = {}
        for e in ('bouton', 'claw', 'step_on', 'step_off'):
            sub = v[(v.element == e) & (v.group == g)]
            if e != 'bouton':
                sub = sub[[(r, float(s)) in keep for r, s in zip(sub.model_reading, sub.share)]]
            res[e] = E90.element_result(sub)
        verdict = E90.group_verdict(res, eff['effect_error'].get(g, 0.0)) if all(r['n_readings'] for r in res.values()) else 'no kept reading'
        rows.append(dict(group=g, verdict=verdict, e90_verdict=eff['groups'].get(g, {}).get('verdict'),
                         **{f'{e}_{k}': res[e][k] for e in res for k in ('differs', 'mixed', 'n_differ', 'n_readings')}))
    return pd.DataFrame(rows)


def per_roi_verdicts(bouton_el: pd.DataFrame, claw_el: pd.DataFrame, kept: Sequence[Tuple[str, float]]) -> pd.DataFrame:
    """The stage elements read on the per-ROI flies (the boutons' three readings, the claws' kept readings; both measures;
    verdict pools) with E90's steps (kept claw readings; E90's flies), by E90's `group_verdict`."""
    v = e90_verdict_rows()
    v = v[v.verdict_subsets == 'primary']
    eff = json.loads((H.shared_dir() / 'results' / 'E90_model_stages' / 'settled_preodor' / 'verdict.json').read_text())
    name = {'noisy': 'noisy_nr'}
    keep = {(name.get(r, r), float(s)) for r, s in kept}
    kset = {(r, float(s)) for r, s in kept}
    rows = []
    for g in GROUPS:
        b = bouton_el[(bouton_el.group == g) & (bouton_el.pool.astype(str) == VERDICT_POOL['bouton'])]
        c = claw_el[(claw_el.group == g) & (claw_el.pool.astype(str) == VERDICT_POOL['claw'])]
        c = c[[(r, float(s)) in kset for r, s in zip(c.reading, c.share)]]
        res = dict(bouton=E90.element_result(b), claw=E90.element_result(c))
        for e in ('step_on', 'step_off'):
            sub = v[(v.element == e) & (v.group == g)]
            res[e] = E90.element_result(sub[[(r, float(s)) in keep for r, s in zip(sub.model_reading, sub.share)]])
        verdict = E90.group_verdict(res, eff['effect_error'].get(g, 0.0)) if res['claw']['n_readings'] else 'no kept reading'
        rows.append(dict(group=g, verdict=verdict, **{f'{e}_{k}': res[e][k] for e in res for k in
                                                      ('differs', 'mixed', 'n_differ', 'n_readings', 'signs')}))
    return pd.DataFrame(rows)


def b3_change(rt_b: pd.DataFrame, rt_c: pd.DataFrame, kept: Sequence[Tuple[str, float]], n_boot: int = N_BOOT) -> pd.DataFrame:
    """B3: the bouton-to-claw change of the per-ROI statistics (VERDICT_CRITERION, APL on): the flies' (claw flies' mean
    minus bouton flies' mean, 90 % interval resampling each stage's flies) and the model's per wiring (the claws at
    VERDICT_POOL in each reading named, minus the boutons at VERDICT_POOL in the same reading with share 0)."""
    stats = ('opr_mean', 'resp_cv', 'top1', 'inhibited', 'totals_cv')
    fb = rt_b[(rt_b.cond == 'flies') & (rt_b.criterion == VERDICT_CRITERION)]
    fc = rt_c[(rt_c.cond == 'flies') & (rt_c.criterion == VERDICT_CRITERION)]
    rng = np.random.default_rng(E90.BOOT_SEED)
    rows = []
    for s in stats + ('rate',):
        a, b = fc[s].to_numpy(dtype=float), fb[s].to_numpy(dtype=float)
        bt = [rng.choice(a, len(a)).mean() - rng.choice(b, len(b)).mean() for _ in range(n_boot)]
        lo, hi = np.percentile(bt, E90.CI)
        rows.append(dict(source='flies', reading='', share=np.nan, stat=s, change=float(a.mean() - b.mean()), lo=float(lo),
                         hi=float(hi), n_inside=np.nan))
    readings = list(dict.fromkeys(list(kept) + [('session_z', 0.0)]))
    mb = rt_b[(rt_b.cond == 'on') & (rt_b.criterion == VERDICT_CRITERION) & (rt_b.pool.astype(str) == VERDICT_POOL['bouton'])]
    mc = rt_c[(rt_c.cond == 'on') & (rt_c.criterion == VERDICT_CRITERION) & (rt_c.pool.astype(str) == VERDICT_POOL['claw'])]
    fl = {r['stat']: r for r in rows}
    for rd, sh in readings:
        cb = mb[(mb.reading == rd)].set_index('wiring')
        cc = mc[(mc.reading == rd) & (mc.share == sh)].set_index('wiring')
        for s in stats:
            d = (cc[s] - cb[s]).reindex(list(WIRINGS)).to_numpy(dtype=float)
            inside = int(((d >= fl[s]['lo']) & (d <= fl[s]['hi'])).sum())
            rows.append(dict(source='model', reading=rd, share=sh, stat=s, change=float(np.nanmean(d)), lo=float(np.nanmin(d)),
                             hi=float(np.nanmax(d)), n_inside=inside))
    return pd.DataFrame(rows)


def part_b() -> dict:
    """Part B's tables (b1_roi_stats.csv, b1_spread.csv, b1_within.csv, b1_readings.csv, b2_flies.csv, b2_element.csv,
    b2_verdicts_e90.csv, b2_verdicts_per_roi.csv, b3_change.csv). Reads part A's a2_roi_stats.csv and a3_element.csv."""
    t0 = time.time()
    od = out_dir()
    rt = _side_cols(roi_table('claw', CLAW_POOLS, conds=('on',)))
    rt.to_csv(od / 'b1_roi_stats.csv', index=False)
    rt.attrs['spread'].to_csv(od / 'b1_spread.csv', index=False)
    within_summary(rt).to_csv(od / 'b1_within.csv', index=False)
    b1 = b1_readings(rt)
    b1.to_csv(od / 'b1_readings.csv', index=False)
    kept = [(r.reading, float(r.share)) for r in b1[(b1.pool.astype(str) == VERDICT_POOL['claw']) & b1.kept].itertuples()]
    fd = fly_distances('claw')
    fd.to_csv(local_dir() / 'b2_fly_pairs.csv')
    ft = []
    for rd in ('trial_mean', 'noise_removed', 'remy_raw'):
        for ms in (MEASURES if rd != 'remy_raw' else ('pearson',)):
            t = group_table_flies(fd, rd, ms).reset_index()
            t.insert(0, 'measure', ms)
            t.insert(0, 'reading', rd)
            ft.append(t)
    pd.concat(ft, ignore_index=True).to_csv(od / 'b2_flies.csv', index=False)
    nine = [(rd, sh) for rd in READINGS for sh in SHARES['claw']]
    el = stage_element('claw', nine, VERDICT_POOL['claw'], fd)
    el.to_csv(od / 'b2_element.csv', index=False)
    e90k = e90_verdicts_kept(kept)
    e90k.to_csv(od / 'b2_verdicts_e90.csv', index=False)
    a3 = pd.read_csv(od / 'a3_element.csv')
    pr = per_roi_verdicts(a3, el, kept)
    pr.to_csv(od / 'b2_verdicts_per_roi.csv', index=False)
    # beside, not decisive: every one of the nine readings, the per-ROI flies' and E90's own (which must give E90's
    # settled verdicts back: the check of the restriction)
    per_roi_verdicts(a3, el, nine).to_csv(od / 'b2_verdicts_per_roi_all9.csv', index=False)
    e9 = e90_verdicts_kept(nine)
    e9.to_csv(od / 'b2_verdicts_e90_all9.csv', index=False)
    if not (e9.verdict == e9.e90_verdict).all():
        raise RuntimeError("E90's verdicts under all nine readings do not reproduce E90's settled verdicts")
    rtb = pd.read_csv(od / 'a2_roi_stats.csv')
    b3 = b3_change(rtb, rt, kept)
    b3.to_csv(od / 'b3_change.csv', index=False)
    provenance('B', CRITERIA=CRITERIA, CLAW_POOLS=CLAW_POOLS, VERDICT_POOL=VERDICT_POOL, WINDOW=WINDOW, SUBSET=SUBSET,
               POOL_SEED=POOL_SEED, NOISE_SEED=NOISE_SEED, kept=kept, seconds=round(time.time() - t0, 1),
               part0_lock=str(roi_dir() / 'lock.json'))
    print(f'part B in {time.time() - t0:.0f} s; kept {kept}', flush=True)
    return dict(rt=rt, b1=b1, kept=kept)


# ---------------------------------------------------------------------------
# part C (the same fly; part 0b's trials)

def same_fly_path() -> Path:
    return roi_dir() / 'same_fly_trials.npz'


def same_fly_trials() -> Dict[str, np.ndarray]:
    """Part 0b's same-fly trials: 'claw_noise', 'soma_noise', 'claw_z', 'soma_z' (repeats x ROIs x ODORS)."""
    d = np.load(same_fly_path(), allow_pickle=False)
    odors = [str(x) for x in d['odors']]
    ix = [odors.index(o) for o in ODORS]
    return {k: d[k][:, :, ix] for k in ('claw_noise', 'soma_noise', 'claw_z', 'soma_z')}


def c1_step(tr: Dict[str, np.ndarray]) -> pd.DataFrame:
    """C1: the same fly's step per pair (somata minus claws), trial means and noise removed, cosine and Pearson (z
    change), with the group means. Rows reading, measure, odor_a, odor_b, claw, soma, step."""
    rows = []
    for rd in ('trial_mean', 'noise_removed'):
        for ms in MEASURES:
            f = trial_mean_distances if rd == 'trial_mean' else noise_removed
            c, s = f(tr['claw_z'], ms), f(tr['soma_z'], ms)
            for (a, b) in c.index:
                rows.append(dict(reading=rd, measure=ms, odor_a=a, odor_b=b, claw=float(c[(a, b)]), soma=float(s[(a, b)]),
                                 step=float(s[(a, b)] - c[(a, b)])))
    return pd.DataFrame(rows)


def c1_groups(c1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (rd, ms), sub in c1.groupby(['reading', 'measure']):
        for col in ('claw', 'soma', 'step'):
            s = pd.Series(sub[col].to_numpy(), index=pd.MultiIndex.from_arrays([sub.odor_a, sub.odor_b]))
            rows.append(dict(reading=rd, measure=ms, what=col, **E90.group_means(s)))
    return pd.DataFrame(rows)


def c2_curve(tr: Dict[str, np.ndarray], thr: float = VERDICT_CRITERION) -> pd.DataFrame:
    """C2: per odor the somata's response rate (trial mean above `thr` noise SDs) against the claws' summed response
    (the mean over claw ROIs of the trial mean in noise units)."""
    s = tr['soma_noise'].mean(axis=0)
    c = tr['claw_noise'].mean(axis=0)
    return pd.DataFrame(dict(odor=list(ODORS), soma_rate=(s > thr).mean(axis=0), claw_total=c.mean(axis=0)))


def model_c2(wiring: str, share: float = 0.0) -> pd.DataFrame:
    """C2's model curve of a wiring (APL on): per odor the KCs' response rate (two or more spikes) against the claws'
    summed drive change (all claws, the mean)."""
    st = model_stage(wiring, 'on')
    U = E90.stage_units(st, 'claw')
    sc = H.abbrev_columns(st['spike_counts'].astype(float))
    return pd.DataFrame(dict(odor=list(ODORS), kc_rate=(sc[list(ODORS)] >= 2).mean(axis=0).to_numpy(),
                             claw_total=U[list(ODORS)].mean(axis=0).to_numpy(), wiring=wiring))


E90_STEPS_FLIES = {'pearson': ('on_remy', 'on_ctrl'), 'cosine': ('on_ctrl_cos',)}
C_MODEL_CLAW = ('session_z', 0.0, 'all')        # the model's step: calibrated somata minus E90's claws in this reading


def c1_compare(c1: pd.DataFrame) -> pd.DataFrame:
    """C1 against E90's steps (settled pre-odor level): per measure and group, the same fly's noise-removed step against
    the between-fly step of the flies (E90's steps_flies.csv, claws noise removed; Remy's somata and Yang's control
    flies for Pearson, Yang's control flies for cosine) and against the model's step with APL on (E90 settled
    steps_groups.csv: the calibrated somata minus the claws in C_MODEL_CLAW, the population and the eight wirings). The
    same fly differs from the between-fly step if it lies outside the between-fly 90 % interval and VERDICT_DIFF or more
    from its mean; the model differs from the same fly by E90's rule without an interval (one fly): VERDICT_DIFF or more
    and 6 of 8 wirings on that side."""
    g = c1_groups(c1)
    sf = g[(g.reading == 'noise_removed') & (g.what == 'step')].set_index('measure')
    fl = pd.read_csv(H.shared_dir() / 'results' / 'E90_model_stages' / 'steps_flies.csv')
    fl = fl[(fl.kind == 'group') & (fl.claw == 'claw_noise_removed')]
    sm = pd.read_csv(H.shared_dir() / 'results' / 'E90_model_stages' / 'settled_preodor' / 'steps_groups.csv')
    rd, sh, sub = C_MODEL_CLAW
    rows = []
    for ms in MEASURES:
        m = sm[(sm.step == 'on') & (sm.soma == f'cal_{ms}') & (sm.reading == rd) & (sm.share == sh) & (sm.subset == sub)]
        mw = m[m.wiring.isin(WIRINGS)].set_index('wiring')
        for grp in GROUPS:
            v = float(sf.loc[ms, grp])
            vals = mw[grp].reindex(list(WIRINGS)).to_numpy(dtype=float)
            pop = float(np.mean(vals))
            d = pop - v
            n_same = int((np.sign(vals - v) == np.sign(d)).sum())
            r = dict(measure=ms, group=grp, same_fly=v, model=pop, model_diff=d, model_n_same_side=n_same,
                     model_differs=bool(abs(d) >= E90.VERDICT_DIFF and n_same >= MIN_SAME))
            for step in E90_STEPS_FLIES[ms]:
                f = fl[(fl.step == step) & (fl.group == grp)].iloc[0]
                out = bool((v < f['lo'] or v > f['hi']) and abs(v - f['mean']) >= E90.VERDICT_DIFF)
                r.update({f'{step}': float(f['mean']), f'{step}_lo': float(f['lo']), f'{step}_hi': float(f['hi']),
                          f'{step}_differs': out})
            rows.append(r)
    return pd.DataFrame(rows)


def c2_summary(c2: pd.DataFrame, cm: pd.DataFrame) -> pd.DataFrame:
    """C2 in numbers: per source (the fly; each wiring) the coefficient of variation over the 12 odors of the somata's
    response rate and of the claws' summed response, and the slope of the response rate on the claws' summed response
    relative to their means (the elasticity, a least-squares fit on the log-free relative values), with the Spearman
    correlation beside as a description."""
    rows = []

    def one(src, rate, tot):
        rate, tot = np.asarray(rate, dtype=float), np.asarray(tot, dtype=float)
        x, y = tot / tot.mean() - 1, rate / rate.mean() - 1
        slope = float(np.polyfit(x, y, 1)[0]) if x.std() > 0 else np.nan
        return dict(source=src, rate_cv=float(rate.std() / rate.mean()), total_cv=float(tot.std() / tot.mean()),
                    elasticity=slope, spearman=float(pd.Series(rate).corr(pd.Series(tot), method='spearman')))
    rows.append(one('fly 2023-04-11 fly03 (megamat1_redo)', c2.soma_rate, c2.claw_total))
    for w, d in cm.groupby('wiring'):
        rows.append(one(w, d.kc_rate, d.claw_total))
    return pd.DataFrame(rows)


def part_c() -> dict:
    """Part C's tables (c1_step.csv, c1_groups.csv, c2_curve.csv, c2_model.csv), from part 0b's same_fly_trials.npz."""
    t0 = time.time()
    od = out_dir()
    tr = same_fly_trials()
    c1 = c1_step(tr)
    c1.to_csv(od / 'c1_step.csv', index=False)
    c1_groups(c1).to_csv(od / 'c1_groups.csv', index=False)
    c1_compare(c1).to_csv(od / 'c1_compare.csv', index=False)
    c2 = c2_curve(tr)
    c2.to_csv(od / 'c2_curve.csv', index=False)
    cm = pd.concat([model_c2(w) for w in WIRINGS], ignore_index=True)
    cm.to_csv(od / 'c2_model.csv', index=False)
    c2_summary(c2, cm).to_csv(od / 'c2_summary.csv', index=False)
    provenance('C', seconds=round(time.time() - t0, 1), same_fly=str(same_fly_path()))
    return dict(c1=c1)


# ---------------------------------------------------------------------------
# the table and the outcome

def _f(x, spec: str = '.3f') -> str:
    try:
        return '–' if x is None or not np.isfinite(float(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def _rng(lo, hi, spec: str = '.3f') -> str:
    return f'{_f(lo, spec)} to {_f(hi, spec)}'


TG = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near', 'all66')
TG_HEAD = ('ester (15)', 'alcohol (6)', 'ester with alcohol (24)', 'other (21)', 'near (6)', 'all 66')


def outcome(od: Path) -> dict:
    """The outcome letters the design allows from the tables (parts A, B and, where part 0b passed, C)."""
    ws = pd.read_csv(od / 'a2_within.csv', dtype={'pool': str})
    w = ws[(ws.criterion == VERDICT_CRITERION) & (ws.pool == VERDICT_POOL['bouton'])]
    on, off = w[w.cond == 'on'].set_index('reading'), w[w.cond == 'off'].set_index('reading')
    b_tot = bool(all(on.loc[r, 'totals_cv_differs'] and on.loc[r, 'totals_cv_side'] == -1 for r in READINGS))
    inh_on = {r: bool(on.loc[r, 'inhibited_differs'] and on.loc[r, 'inhibited_side'] == 1) for r in READINGS}
    inh_off = {r: bool(off.loc[r, 'inhibited_differs'] and off.loc[r, 'inhibited_side'] == 1) for r in READINGS}
    cv_on = {r: bool(on.loc[r, 'resp_cv_differs']) for r in READINGS}
    a3 = pd.read_csv(od / 'a3_element.csv', dtype={'pool': str})
    a3 = a3[a3.pool == VERDICT_POOL['bouton']]
    a3_any = bool(a3.differs.any())
    b1 = pd.read_csv(od / 'b1_readings.csv', dtype={'pool': str})
    kept = b1[(b1.pool == VERDICT_POOL['claw']) & b1.kept]
    letters = {}
    letters['a'] = dict(holds=bool(not b_tot and not a3_any and not any(inh_on.values())),
                        text="the boutons match in strength, rates and distances")
    letters['b'] = dict(holds=b_tot, text="the boutons' per-odor totals are more spread in the flies than in the model, in "
                        "every reading (below the flies' range on 6 or more of 8 wirings)")
    c_all = all(inh_on.values()) and not any(inh_off.values())
    c_some = any(inh_on.values()) and not any(inh_off.values())
    letters['c'] = dict(holds=bool(c_all), partial=bool(c_some and not c_all),
                        readings_inhibited_on=[r for r, v in inh_on.items() if v],
                        text="the model's boutons with APL on are more inhibited below baseline than the flies' (and not "
                        "with APL off): " + ('in every reading' if c_all else 'in ' + ', '.join(r for r, v in inh_on.items() if v)
                                              + ' only' if c_some else 'in no reading'))
    letters['d'] = dict(holds=bool(len(kept) == 1), kept=[(r.reading, float(r.share)) for r in kept.itertuples()],
                        text='one claw reading is kept')
    letters['e'] = dict(holds=bool(len(kept) == 0), text='no claw reading is kept (pooled ROIs cannot settle the claws)')
    p0 = json.loads((roi_dir() / 'lock.json').read_text())['checks']['passed']
    letters['g'] = dict(holds=not p0, text="part 0's checks fail")
    c1p = od / 'c1_compare.csv'
    sfc = roi_dir() / 'same_fly_checks.json'
    if c1p.exists() and sfc.exists() and json.loads(sfc.read_text()).get('passed'):
        c1 = pd.read_csv(c1p)
        cols = [c for c in c1.columns if c.endswith('_differs') and c != 'model_differs']
        diff_groups = sorted({f'{r.group} ({r.measure}, {c[:-8]})' for r in c1.itertuples() for c in cols
                              if str(getattr(r, c)) == 'True'})
        letters['f'] = dict(holds=bool(len(diff_groups) > 0), groups=diff_groups,
                            text='the same-fly step differs from the between-fly step (one fly; a lead until more)')
    else:
        letters['f'] = dict(holds=None, text='part C not run')
    return dict(letters=letters, holds=[k for k, v in letters.items() if v.get('holds')],
                cv_narrower_on=cv_on, inhibited_off=inh_off)


def write_table() -> None:
    """comparison_table.md and outcome.json from the parts' CSVs."""
    od = out_dir()
    L: List[str] = ['# E119: the model\'s boutons, claws and claw-to-soma step against the per-ROI recordings', '']
    lock = json.loads((roi_dir() / 'lock.json').read_text())
    ck = pd.read_csv(roi_dir() / 'checks.csv')
    L += ['## Part 0: the data and the checks', '',
          '| check | stage | flies | passed | worst value |', '|---|---|---:|---:|---:|']
    for (c, st), sub in ck.groupby(['check', 'stage']):
        worst = sub.value.max() if c != '2_dashboard' else sub[~sub.note.str.startswith('excluded')].value.min()
        L.append(f'| {c} | {st} | {sub.acq.nunique()} | {int(sub.passed.sum())} of {len(sub)} | {_f(worst, ".2e" if c != "2_dashboard" else ".4f")} |')
    ex = ck[(ck.check == '2_dashboard') & ck.note.str.startswith('excluded')]
    L += ['', f"Excluded from every analysis: {', '.join(ex.acq)} (best r with any dashboard fly "
          f"{_f(ex.value.max())}; Remy's dashboards leave it out). Check 2 maps the other five bouton flies to dashboard "
          'flies 1 to 5 in date order.', '']
    sfc = roi_dir() / 'same_fly_checks.json'
    if sfc.exists():
        s0 = json.loads(sfc.read_text())
        L += [f"Part 0b (the same fly, 2023-04-11 fly 3, megamat1_redo; {s0['n_claw']} claw and {s0['n_soma']} soma ROIs): "
              f"the decoded same-fly series reproduced at r {_f(s0['claws_r'])} (claws) and {_f(s0['somata_r'])} (somata), "
              f"the claws-against-somata agreement over 1 to 2 s {_f(s0['agreement_1_2s'])} (Remy: 0.86); "
              f"{'passed' if s0['passed'] else 'FAILED'}.", '']
    # A1 / A2
    sp = pd.read_csv(od / 'a2_spread.csv').set_index('criterion').loc[VERDICT_CRITERION]
    rt = pd.read_csv(od / 'a2_roi_stats.csv', dtype={'pool': str})
    m = rt[(rt.criterion == VERDICT_CRITERION) & (rt.cond != 'flies') & (rt.pool == VERDICT_POOL['bouton'])]
    ws = pd.read_csv(od / 'a2_within.csv', dtype={'pool': str})
    ws = ws[(ws.criterion == VERDICT_CRITERION) & (ws.pool == VERDICT_POOL['bouton'])].set_index(['reading', 'cond'])
    e114 = e114_exponent()
    L += ['## Part A: the boutons (5 flies)', '',
          f'Per-ROI statistics at {VERDICT_CRITERION:g} baseline-noise SDs (the flies) and at the model threshold that gives '
          f"the flies' mean overall response rate ({_f(sp['rate_mean'])}), model boutons not pooled, E90's shared-glomerulus "
          "subset. Model cells: mean over the 8 wirings, and in brackets the number of wirings below / above the flies' "
          'range (or inside the band, for the profiles).', '',
          '| statistic | flies (range) | ' + ' | '.join(f'{r}<br>APL {c}' for r in READINGS for c in ('on', 'off')) + ' |',
          '|---|---:|' + '---:|' * (2 * len(READINGS))]
    for stat, name in (('totals_cv', "per-odor totals' CV over the 12 odors (A1)"), ('resp_cv', "responders' CV"),
                       ('top1', 'top-1 % share'), ('inhibited', 'share of entries below -threshold'),
                       ('opr_mean', 'odors per ROI, mean')):
        cells = []
        for r in READINGS:
            for c in ('on', 'off'):
                v = m[(m.reading == r) & (m.cond == c)][stat].mean()
                w = ws.loc[(r, c)]
                cells.append(f"{_f(v)} ({int(w[f'{stat}_below'])} / {int(w[f'{stat}_above'])})")
        L.append(f"| {name} | {_f(sp[f'{stat}_mean'])} ({_rng(sp[f'{stat}_lo'], sp[f'{stat}_hi'])}) | " + ' | '.join(cells) + ' |')
    for prof, name, band in (('per_odor', 'per-odor response rates, RMS from the flies', 'per_odor_band'),
                             ('opr', 'odors per ROI, total variation from the flies', 'opr_band')):
        cells = []
        for r in READINGS:
            for c in ('on', 'off'):
                v = m[(m.reading == r) & (m.cond == c)][f'{prof}_dist'].mean()
                cells.append(f"{_f(v)} ({int(ws.loc[(r, c)][f'{prof}_n_in'])} in)")
        L.append(f'| {name} | band {_f(sp[band])} | ' + ' | '.join(cells) + ' |')
    if e114:
        L += ['', f"Beside A1: E114's exponent of the compressive map from the ORN recordings to the input rates is "
              f"{e114.get('exponent')} ({e114.get('reason')}; E114 still running). E90's PN-dendrite stage: the flies' CV 0.312 "
              'against the model\'s 0.074.']
    # A3
    fl = pd.read_csv(od / 'a3_flies.csv').set_index(['reading', 'measure', 'group'])
    el = pd.read_csv(od / 'a3_element.csv', dtype={'pool': str})
    el = el[el.pool == VERDICT_POOL['bouton']]
    L += ['', '### A3: distances across bouton ROIs, per fly, E90\'s groups', '',
          "Flies: the mean over the 5 flies of each fly's group mean (90 % interval over resampled flies), noise removed by "
          "E81's estimator; model: the population over 8 wirings (mean of the wirings' group means), APL on; an asterisk "
          "marks \"differs\" by E90's rule.", '',
          '| measure, reading | ' + ' | '.join(TG_HEAD) + ' |', '|---|' + '---:|' * len(TG)]
    for ms in MEASURES:
        L.append(f'| {ms}, flies | ' + ' | '.join(f"{_f(fl.loc[('noise_removed', ms, g), 'mean'])} "
                 f"({_rng(fl.loc[('noise_removed', ms, g), 'lo'], fl.loc[('noise_removed', ms, g), 'hi'], '.2f')})" for g in TG) + ' |')
        for rd in READINGS:
            sub = el[(el.measure == ms) & (el.reading == rd)].set_index('group')
            L.append(f'| {ms}, model {rd} | ' + ' | '.join(f"{_f(sub.loc[g, 'model'])}{'*' if sub.loc[g, 'differs'] else ''}"
                                                            for g in TG) + ' |')
    eg = element_by_group(el)
    L += ['', 'Bouton element over the 3 readings x 2 measures: ' + '; '.join(
        f"{g} {'differs' if eg.loc[g, 'differs'] else ('mixed' if eg.loc[g, 'mixed'] else 'same')} "
        f"({eg.loc[g, 'n_differ']} of {eg.loc[g, 'n_readings']})" for g in TG) + '.']
    a4 = pd.read_csv(od / 'a4_summary.csv')
    L += ['', '### A4 (a lead): bouton groups against the model\'s glomeruli', '',
          '| fly | ROIs | responsive | tuning groups (>= 5 ROIs) | co-variation groups | share of grouped ROIs whose best model glomerulus is not above chance |',
          '|---|---:|---:|---:|---:|---:|']
    for r in a4.itertuples():
        L.append(f'| {r.acq} | {r.n_rois} | {r.n_responsive} | {r.n_tuning_groups} ({r.n_big_groups}) | {r.n_covar_groups} | {_f(r.share_in_unmatched, ".2f")} |')
    L += ['', f"Chance: the best Pearson r of 54 model glomeruli with a group's 12-odor tuning, odors permuted, is about "
          f"{_f(pd.read_csv(od / 'a4_groups.csv').null95.median(), '.2f')} (median 95th percentile), so a match on 12 odors is weak "
          'evidence; every best match was a glomerulus with ORN input.']
    # B
    b1 = pd.read_csv(od / 'b1_readings.csv', dtype={'pool': str})
    spc = pd.read_csv(od / 'b1_spread.csv').set_index('criterion').loc[VERDICT_CRITERION]
    L += ['', '## Part B: the claws (4 flies)', '', '### B1: E90\'s nine claw readings on the per-ROI statistics', '',
          f"Model claws pooled in groups of 14 of one glomerulus ({VERDICT_POOL['claw']}), APL on, at the flies' mean rate at "
          f"{VERDICT_CRITERION:g} SDs ({_f(spc['rate_mean'])}). Flies: per-odor response rates band {_f(spc['per_odor_band'])}, "
          f"odors per ROI band {_f(spc['opr_band'])}, responders' CV {_rng(spc['resp_cv_lo'], spc['resp_cv_hi'])}, share "
          f"below -threshold {_rng(spc['inhibited_lo'], spc['inhibited_hi'])}. A reading is kept if all three statistics are "
          'inside on 6 of 8 wirings.', '',
          "| reading | share | per-odor rates inside | odors per ROI inside | responders' CV inside | responders' CV | below -threshold | kept |",
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in b1[b1.pool == VERDICT_POOL['claw']].itertuples():
        L.append(f'| {r.reading} | {r.share:g} | {r.n_per_odor_in} | {r.n_opr_in} | {r.n_resp_cv_in} | {_f(r.resp_cv)} | '
                 f"{_f(r.inhibited)} | {'yes' if r.kept else 'no'} |")
    other = b1[b1.pool != VERDICT_POOL['claw']]
    L += ['', f"Other poolings (1, g7, g28, r7, r14, r28): kept in {int(other.kept.sum())} of {len(other)} reading-pool pairs; "
          f"responders' CV {_rng(other.resp_cv.min(), other.resp_cv.max())} against the flies' {_rng(spc['resp_cv_lo'], spc['resp_cv_hi'])}."]
    fc = pd.read_csv(od / 'b2_flies.csv').set_index(['reading', 'measure', 'group'])
    e2 = pd.read_csv(od / 'b2_element.csv')
    L += ['', '### B2: the claw distances on the per-ROI flies (all nine readings; none kept, so none decides)', '',
          '| measure, reading | ' + ' | '.join(TG_HEAD) + ' |', '|---|' + '---:|' * len(TG)]
    for ms in MEASURES:
        L.append(f'| {ms}, flies | ' + ' | '.join(f"{_f(fc.loc[('noise_removed', ms, g), 'mean'])} "
                 f"({_rng(fc.loc[('noise_removed', ms, g), 'lo'], fc.loc[('noise_removed', ms, g), 'hi'], '.2f')})" for g in TG) + ' |')
        for rd in READINGS:
            for sh in SHARES['claw']:
                sub = e2[(e2.measure == ms) & (e2.reading == rd) & (e2.share == sh)].set_index('group')
                L.append(f'| {ms}, model {rd} {sh:g} | ' + ' | '.join(
                    f"{_f(sub.loc[g, 'model'])}{'*' if sub.loc[g, 'differs'] else ''}" for g in TG) + ' |')
    v9 = pd.read_csv(od / 'b2_verdicts_per_roi_all9.csv').set_index('group')
    vk = pd.read_csv(od / 'b2_verdicts_e90.csv').set_index('group')
    L += ['', "Group verdicts. Under the kept readings only: " + '; '.join(f'{g} {vk.loc[g, "verdict"]}' for g in TG) +
          ". Beside, not deciding: with all nine readings on the per-ROI flies (E90's steps): " +
          '; '.join(f'{g} {v9.loc[g, "verdict"]}' for g in TG) + ". E90's settled verdicts (all nine readings, decoded "
          'matrices), reproduced by this script: ' + '; '.join(f'{g} {vk.loc[g, "e90_verdict"]}' for g in TG) + '.']
    b3 = pd.read_csv(od / 'b3_change.csv')
    L += ['', '### B3: the bouton-to-claw change', '',
          "| statistic | flies: claws minus boutons (90 % interval) | model (session z, claws g14 minus boutons): mean (range over wirings) |",
          '|---|---:|---:|']
    f3 = b3[b3.source == 'flies'].set_index('stat')
    m3 = b3[(b3.source == 'model') & (b3.reading == 'session_z')].set_index('stat')
    for st, name in (('totals_cv', "per-odor totals' CV"), ('resp_cv', "responders' CV"), ('top1', 'top-1 % share'),
                     ('inhibited', 'share below -threshold'), ('opr_mean', 'odors per ROI, mean'), ('rate', 'overall response rate')):
        mc = (f"{_f(m3.loc[st, 'change'])} ({_rng(m3.loc[st, 'lo'], m3.loc[st, 'hi'])})" if st in m3.index
              else 'matched by construction')
        L.append(f"| {name} | {_f(f3.loc[st, 'change'])} ({_rng(f3.loc[st, 'lo'], f3.loc[st, 'hi'])}) | {mc} |")
    # C
    if (od / 'c1_compare.csv').exists():
        c1 = pd.read_csv(od / 'c1_compare.csv')
        L += ['', '## Part C: the same fly\'s claw-to-soma step (one fly)', '',
              'C1: the step (somata minus claws), noise removed, per group; between-fly steps from E90 (claws noise removed); '
              "the model's step: calibrated somata minus E90's claws (session z, all claws), APL on, population of 8 wirings. "
              'An asterisk marks a between-fly step the same fly lies outside of (by VERDICT_DIFF or more), or a model step that '
              'differs from the same fly on 6 of 8 wirings.', '',
              '| measure, source | ' + ' | '.join(TG_HEAD) + ' |', '|---|' + '---:|' * len(TG)]
        for ms in MEASURES:
            sub = c1[c1.measure == ms].set_index('group')
            L.append(f'| {ms}, same fly | ' + ' | '.join(_f(sub.loc[g, 'same_fly']) for g in TG) + ' |')
            for step in E90_STEPS_FLIES[ms]:
                L.append(f'| {ms}, between flies ({step}) | ' + ' | '.join(
                    f"{_f(sub.loc[g, step])} ({_rng(sub.loc[g, step + '_lo'], sub.loc[g, step + '_hi'], '.2f')}){'*' if sub.loc[g, step + '_differs'] else ''}"
                    for g in TG) + ' |')
            L.append(f'| {ms}, model | ' + ' | '.join(f"{_f(sub.loc[g, 'model'])}{'*' if sub.loc[g, 'model_differs'] else ''}"
                                                     for g in TG) + ' |')
        c2 = pd.read_csv(od / 'c2_summary.csv')
        L += ['', "C2: across the 12 odors, the somata's response rate against the claws' summed response (the fly: somata "
              f"above {VERDICT_CRITERION:g} noise SDs, claws in noise units; the model: KCs with 2 or more spikes, the claws' summed "
              'drive change). Spearman beside, as a description.', '',
              "| source | response rate's CV | claws' summed response, CV | elasticity | Spearman |", '|---|---:|---:|---:|---:|']
        for r in c2.itertuples():
            L.append(f'| {r.source} | {_f(r.rate_cv)} | {_f(r.total_cv)} | {_f(r.elasticity, ".2f")} | {_f(r.spearman, ".2f")} |')
    oc = outcome(od)
    _write_json(od / 'outcome.json', oc)
    L += ['', '## Outcome', '']
    for k in sorted(oc['letters']):
        v = oc['letters'][k]
        state = 'holds' if v.get('holds') else ('partly' if v.get('partial') else ('not run' if v.get('holds') is None else 'no'))
        L.append(f'- ({k}) {state}: {v["text"]}' + (f" ({', '.join(v['groups'])})" if v.get('groups') else ''))
    (od / 'comparison_table.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L), flush=True)


# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('A', 'B', 'C', 'table'))
    a = ap.parse_args(argv)
    if a.part == 'A':
        part_a()
    elif a.part == 'B':
        part_b()
    elif a.part == 'C':
        part_c()
    else:
        write_table()


if __name__ == '__main__':
    main()
