#!/usr/bin/env python3
"""E92: the step from claws to somata with APL on, read like for like: the flies' somata with their trial noise removed
against the models' somata as the recordings are normalized, the line along which every saved form of APL's inhibition
trades the ester pairs against the alcohol pairs, and then APL's inhibition of the KCs as a gain on the claw's drive.

The entry is item E92 of experiment_priorities.md (proposed by hal 2026-11-13 19:04; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-13 18:38). C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the
0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit`). S: real, P1, P2, M1. T: P3, P4,
M2, M3. Halves A / B of the recordings E77's (`recording_halves()`, seed 0). The 66 pairs of the panel's 12 megamat
odors in E90's groups (ester 15, alcohol 6, ester with alcohol 24, other 21, E78's six near pairs, all 66), the six near
pairs also as the four ester near pairs together and 1-5ol with 1-6ol and 6al with t2h singly (GROUPS).

This module holds part A (no simulations; it reads E90's archived runs of C2's family, Yang's per-trial recordings and
the archived runs of E77 to E89, and nothing of E91's). Parts B to F are in experiments.E92_gain; their names here are
stubs.

Parts (CLI `--part`):
- A1: the flies' somata (somata_flies.csv): per recording of Yang's control flies and flies with APL output blocked
  (TNT), in Yang's z and in baseline-noise units, the pair distances as trial means (E90's reading) and with the trial
  noise removed (E81's estimate from products of different repeats, `H.noise_free_pair_distances`), as cosine and, on
  single trials centred over the KCs, as Pearson; the mean over recordings per pair and group with 90 % intervals from
  N_BOOT resamples of the recordings. Remy's wild-type somata (per-fly correlation matrices only) as measured, flagged.
- Acheck: the known-truth check of the noise removal (estimate_check.csv): E89's pseudo-recordings of the model
  (`E89.ModelReading89`: each of the eight wirings' runs under E73's per-KC session z-score, saturating map, ten
  pseudo-recordings of three trials, at E89's binding noise levels CHECK_LEVELS, the flies' repeatability and a share of
  the noise common to the odors of a repeat) through the same estimate, against the model's noise-free distances; a
  group passes in a measure and units if the population's estimate is within CHECK_TOL of the truth at every level in
  both conditions. A group that fails is not read with the noise removed.
- A2: the model's somata as the recordings are normalized (somata_model.csv): the eight wirings, B5's point at the
  flies' 2 s pulse and D1's four firing-ceiling fits (E90's archived counts), read under E73's per-KC session z-score
  (calibrated on each run's APL-on counts, saturating and linear maps; `model_calibrations`), in the reading's z (as
  Yang's z) and in baseline-noise units (the amplitude over the common frame noise), without noise (against A1's
  noise-removed values) and through the reading's measurement noise in expectation (against the flies' trial means;
  `distance_matrix` with `sd`), as Pearson and cosine; beside, raw counts and the calibrated readout (E90's
  `soma_rows`). The linear map in baseline-noise units without noise is the raw counts exactly. Per run and as
  populations (all8, S, T; fc_S4 for D1's four).
- A3: the verdict (verdict_rows_A.csv, steps_like_for_like.csv, verdict_A.json): by E90's rule (`E90.differs`: the
  population's group mean VERDICT_DIFF or more from the flies', outside their 90 % interval, at least 6 of the 8 wirings
  on that side) per group, condition and like-for-like reading (`like_for_like`: 2 maps x 2 units x 2 noise readings x
  2 measures), the somata first (no claw reading), then the two steps (somata minus E90's claws: the model's claws as E90
  read them, session z without noise or through the noise as measured, all claws, drive alone; Remy's claws noise
  removed or as measured). Per group "with APL on", "without APL", "both", "neither" or "not resolved", with the sign; a
  statement holds only if every like-for-like reading that passed the check agrees. Beside, not in the verdict: raw counts
  and the calibrated readout against the flies in baseline-noise units, B5 against the real wiring at 0.5 s, D1's four.
- A4: the line (line.json, line.csv): for every archived variant of C2's family (`a4_variants`) the ester and the
  alcohol pairs' mean distance with APL on, each minus the mean over the 120 pairs' (`relative_distances`), on raw counts
  and under the calibrated readout (K set on each run's APL-on counts), as populations where four wirings exist (S, T),
  single runs otherwise; the line fitted through E89's 13 weightings (C1) as S populations (`line`), and one line per S
  wiring through the same weightings for single runs; a variant is off the line toward the flies (`place_on_line`) if
  its residual, signed toward the side of the line the control flies' point lies on, exceeds the weightings' scatter
  about the line (the residuals' SD), with the flies' point from half A and from half B (cosine trial means, the main
  target's reading; noise removed beside). On T the residual is taken relative to C2's unchanged population on T (E90's
  runs), against S's line. Every variant carries its pattern part (calibrated readout, all recordings) and its realism
  (E83's `realism`: inside the bands if every member run is); one with an odor silent in both conditions (NaN
  distances, `n_pairs` below 120) is placed but never named. No variant is chosen by this.
- A: A1, Acheck, A2, A3 and A4 in order, then table; table: comparison_table.md and the figures.
- B, C, D, E, F: implemented in experiments.E92_gain.

Functions for experiments.E92_gain: `somata_readings(on, off)` (one run's somata in A2's readings, per group and
reading), `line(readout)` (A4's lines, from line.json once part A4 has run), `place_on_line(est_rel, alc_rel, ln)` and
`relative_distances(d_on)` (one variant's place against a line).

Outputs: results/E92_steps_like_for_like/ (somata_flies.csv, estimate_check.csv, somata_model.csv,
somata_model_calibration.csv, verdict_rows_A.csv, steps_like_for_like.csv, verdict_A.json, line.json, line.csv,
comparison_table.md, E92_somata.png, E92_line.png, provenance_<part>.json); per-recording and per-run pair tables under
outputs/E92_steps_like_for_like/.

Usage (from the repository root):
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E92_steps_like_for_like.py --part A --workers 8
"""
from __future__ import annotations

import os

# (only when run as a script: experiments.E92_gain imports this module and keeps its own thread count; the worker
# processes inherit the setting)
if __name__ == '__main__':
    for _v in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(_v, '1')

import argparse
import json
import re
import time
import warnings
from itertools import combinations
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E90_model_stages as E90

EXPERIMENT = 'E92_steps_like_for_like'
ODORS16 = tuple(H.APL16_ODORS)
S_WIRINGS = ('real', 'P1', 'P2', 'M1')
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
ALL_WIRINGS = E90.ALL_WIRINGS
CONDS = ('on', 'off')
GENO = {'on': 'ctrl', 'off': 'tnt'}
ARCHIVE = H.shared_dir() / 'data' / 'model_outputs'

# ---------------------------------------------------------------------------
# groups

NEAR_ESTER = (('eb', 'ep'), ('IaA', 'pa'), ('2h', 'pa'), ('2h', 'IaA'))
NEAR_SINGLE = {'near_15ol_16ol': ('1-5ol', '1-6ol'), 'near_6al_t2h': ('6al', 't2h')}
GROUPS = E90.GROUPS + ('near_ester',) + tuple(NEAR_SINGLE)
GROUP_NAMES = dict(E90.GROUP_NAMES, near_ester='the four ester near pairs', near_15ol_16ol='1-5ol with 1-6ol',
                   near_6al_t2h='6al with t2h')
TABLE_GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near', 'near_ester', 'near_15ol_16ol', 'near_6al_t2h')
TABLE_HEAD = dict(E90.TABLE_GROUP_HEAD, near_ester='four ester near', near_15ol_16ol='1-5ol, 1-6ol',
                  near_6al_t2h='6al, t2h')

# ---------------------------------------------------------------------------
# A1 and the check

UNITS = ('z', 'noise')                  # Yang's whole-recording z; baseline-noise units
MEASURES = ('pearson', 'cosine')
FLY_READINGS = ('trial_mean', 'noise_removed')
N_BOOT = E90.N_BOOT
CI = E90.CI
BOOT_SEED = E90.BOOT_SEED
CHECK_TOL = 0.03
CHECK_LEVELS = ('L1', 'L3h', 'L3', 'L3x')   # E89's BINDING_LEVELS
REMY_FLAG = "as measured: Remy's wild-type somata exist only as per-fly correlation matrices"

# ---------------------------------------------------------------------------
# A2

MAPS = ('saturating', 'linear')         # E73's two amplitude maps
NOISE = ('none', 'expected')            # without noise / through the reading's measurement noise in expectation
D1_WIRING = {'fit1': 'real', 'P1_fc': 'P1', 'P2_fc': 'P2', 'M1_fc': 'M1'}

# ---------------------------------------------------------------------------
# A3

# the groups the outcomes are read on (the two single near pairs are leads)
OUTCOME_GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near_ester')
VERDICT_DIFF = E90.VERDICT_DIFF
VERDICT_SAME_SIDE = E90.VERDICT_SAME_SIDE
CLAW_READING = {'none': ('session_z', 'noise_removed'), 'expected': ('noisy_measured', 'measured')}   # model, flies

# ---------------------------------------------------------------------------
# A4

A4_FOLDERS = ('E77_input_side_common_error', 'E78_claw_knee', 'E80_remaining_candidates', 'E82_flies_protocol',
              'E83_cross_wiring_fit', 'E85_orn_fly_population', 'E86_reduced_input_candidate', 'E87_inhibition_form',
              'E88_responder_sets', 'E89_apl_breadth_weighting')
# variants of other families in those folders: the firing-ceiling fits (E69's fits 0 and 1, E72 part G's fc fits, E82's
# refit of that family and its D1 runs)
NOT_C2 = {'E78_claw_knee': r'^fit1_', 'E80_remaining_candidates': r'^fit0_',
          'E82_flies_protocol': r'fit0|fit1|_fc|^nu_|D1'}
SKIP_DIRS = ('smoke', 'e10_checks', 'backup')
WIRING_RE = r'^(.*)_(real|P1|P2|P3|P4|M1|M2|M3)$'
LINE_EXPERIMENT, LINE_PART = 'E89_apl_breadth_weighting', 'C1'
READOUTS = ('raw', 'cal')
FLY_POINTS = ('all', 'A', 'B')


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def provenance(part: str, **settings) -> None:
    """provenance_<part>.json: the script, the commits (`H.code_commits()`), the settings."""
    _write_json(out_dir() / f'provenance_{part}.json', dict(part=part, script=H.script_path(__file__),
                started=time.strftime('%Y-%m-%d %H:%M:%S%z'), **H.code_commits(), settings=settings))


def _pool(workers: int):
    from concurrent.futures import ProcessPoolExecutor
    return ProcessPoolExecutor(max_workers=max(1, int(workers)))


def _map(fn: Callable, jobs: Sequence, workers: int) -> List:
    """fn over jobs, in `workers` processes (in this process if workers <= 1), results in job order."""
    if workers <= 1 or len(jobs) <= 1:
        return [fn(j) for j in jobs]
    with _pool(min(workers, len(jobs))) as ex:
        return list(ex.map(fn, jobs))


# ---------------------------------------------------------------------------
# pairs and groups

PAIRS120 = pd.MultiIndex.from_tuples(sorted({tuple(sorted(p)) for p in combinations(ODORS16, 2)}),
                                     names=['odor_a', 'odor_b'])
PAIRS66 = E90.megamat_pairs()


def group_masks(pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """E90's group masks over `pairs` (canonical order) plus the near pairs split (GROUPS)."""
    m = E90.group_masks(pairs)
    keys = [tuple(sorted(map(str, x))) for x in pairs]
    ne = {tuple(sorted(p)) for p in NEAR_ESTER}
    m['near_ester'] = np.array([k in ne for k in keys])
    for g, p in NEAR_SINGLE.items():
        m[g] = np.array([k == tuple(sorted(p)) for k in keys])
    return m


def group_means(pairs: pd.Series) -> Dict[str, float]:
    """A pair Series' mean over each group's pairs (GROUPS; NaN pairs left out)."""
    s = E90._cs(pairs)
    m = group_masks(s.index)
    v = s.to_numpy(dtype=float)
    out = {}
    for g in GROUPS:
        x = v[m[g]]
        out[g] = float(np.nanmean(x)) if np.isfinite(x).any() else np.nan
    return out


def _group_boot(boot: np.ndarray, pairs: pd.MultiIndex) -> Dict[str, np.ndarray]:
    m = group_masks(pairs)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return {g: np.nanmean(boot[:, m[g]], axis=1) for g in GROUPS}


def summary_rows(key_cols: dict, mean: pd.Series, boot: np.ndarray, n: int) -> List[dict]:
    """Per pair and per group (GROUPS): the mean, its 90 % interval and SD over the resamples `boot` (resamples x the
    pairs of `mean`, in its order)."""
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        lo, hi = np.nanpercentile(boot, CI, axis=0)
        sd = np.nanstd(boot, axis=0)
    for i, (a, b) in enumerate(mean.index):
        rows.append(dict(key_cols, n=n, kind='pair', odor_a=a, odor_b=b, group='', mean=float(mean.iloc[i]),
                         lo=float(lo[i]), hi=float(hi[i]), sd_boot=float(sd[i])))
    gm = group_means(mean)
    gb = _group_boot(boot, mean.index)
    for g in GROUPS:
        glo, ghi = np.nanpercentile(gb[g], CI)
        rows.append(dict(key_cols, n=n, kind='group', odor_a='', odor_b='', group=g, mean=gm[g], lo=float(glo),
                         hi=float(ghi), sd_boot=float(np.nanstd(gb[g]))))
    return rows


def _square_pairs(D: np.ndarray, odors: Sequence[str], index: pd.MultiIndex) -> pd.Series:
    pos = {o: i for i, o in enumerate(odors)}
    return pd.Series([float(D[pos[a], pos[b]]) for a, b in index], index=index)


def distance_matrix(X: np.ndarray, measure: str, sd=None) -> np.ndarray:
    """odors x odors distances between the columns (odors) of X (units x odors): 1 - cosine, or 1 - Pearson r (the
    columns centred over the units). With `sd` (per unit, or one value), in expectation over independent measurement
    noise of that SD added to every entry, as a ratio of expectations (`zcell_pair_distances`' closed form): the noise
    adds sum(sd^2) to each odor's squared length (centred: (1 - 1 / n_units) sum(sd^2)) and nothing to the products of
    two odors."""
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    if measure == 'pearson':
        Xc = X - X.mean(axis=0, keepdims=True)
        f = 1.0 - 1.0 / n
    elif measure == 'cosine':
        Xc = X
        f = 1.0
    else:
        raise ValueError(measure)
    G = Xc.T @ Xc
    n2 = np.diag(G).copy()
    if sd is not None:
        n2 = n2 + f * float(np.sum(np.broadcast_to(np.asarray(sd, dtype=float), (n,)) ** 2))
    with np.errstate(invalid='ignore', divide='ignore'):
        return 1.0 - G / np.sqrt(np.outer(n2, n2))


def _abbrev16(counts: pd.DataFrame) -> pd.DataFrame:
    return H.abbrev_columns(counts.astype(float))[list(ODORS16)]


# ---------------------------------------------------------------------------
# A1: the flies

def centred(trials: np.ndarray) -> np.ndarray:
    """Single trials (repeats, KCs, odors) with each trial's mean over the KCs removed, so that their cosine is Pearson's
    r."""
    x = np.asarray(trials, dtype=float)
    return x - x.mean(axis=1, keepdims=True)


def noise_removed(trials: np.ndarray, measure: str, odors: Sequence[str] = ODORS16) -> pd.Series:
    """The noise-removed 120 pair distances of one recording's single trials (repeats, units, odors): E81's estimate
    from products of different repeats (`H.noise_free_pair_distances`), on the trials as they are (cosine) or centred
    over the units (Pearson)."""
    x = centred(trials) if measure == 'pearson' else np.asarray(trials, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        d = H.noise_free_pair_distances(x, list(odors))
    return E90._cs(d).reindex(PAIRS120)


def recording_distances(trials: np.ndarray, odors: Sequence[str] = ODORS16) -> Dict[Tuple[str, str], pd.Series]:
    """One recording's 120 pair distances {(reading, measure): pair Series} for FLY_READINGS x MEASURES: 'trial_mean' the
    distances of the trial means (E90's reading), 'noise_removed' (`noise_removed`)."""
    xbar = np.asarray(trials, dtype=float).mean(axis=0)
    out = {}
    for ms in MEASURES:
        out[('trial_mean', ms)] = _square_pairs(distance_matrix(xbar, ms), odors, PAIRS120)
        out[('noise_removed', ms)] = noise_removed(trials, ms, odors)
    return out


_FLIES: Dict[str, Dict[Tuple[str, str, str], pd.DataFrame]] = {}


def recording_label(rec) -> str:
    return '/'.join(map(str, rec))


def flies_samples(units: str) -> Dict[Tuple[str, str, str], pd.DataFrame]:
    """{(genotype key 'ctrl' / 'tnt', reading, measure): recordings x the 120 pairs} in `units` (Yang's per-trial values,
    `H.yang_recording_trials`; rows labelled 'date/fly/hemibrain'); cached."""
    if units not in _FLIES:
        out = {}
        for g, gk in (('TNT3c', 'ctrl'), ('TNT', 'tnt')):
            rows: Dict[tuple, dict] = {}
            for rec, X in H.yang_recording_trials(g, units).items():
                for k, s in recording_distances(X).items():
                    rows.setdefault(k, {})[recording_label(rec)] = s
            for (rd, ms), d in rows.items():
                out[(gk, rd, ms)] = pd.DataFrame(d).T
        _FLIES[units] = out
    return _FLIES[units]


def part_a1() -> pd.DataFrame:
    """somata_flies.csv: per genotype, units, reading and measure, per pair (all 120) and group the mean over
    recordings, its 90 % interval and SD from N_BOOT resamples of the recordings; Remy's wild type as measured (the 66
    megamat pairs, Pearson, flagged). The per-recording tables go to outputs/E92_steps_like_for_like/A1/."""
    t0 = time.time()
    provenance('A1', N_BOOT=N_BOOT, CI=CI, BOOT_SEED=BOOT_SEED, UNITS=UNITS, MEASURES=MEASURES)
    rows = []
    long = []
    for units in UNITS:
        for (gk, rd, ms), X in flies_samples(units).items():
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                mean = pd.Series(np.nanmean(X.to_numpy(dtype=float), axis=0), index=X.columns)
            rows += summary_rows(dict(genotype=gk, units=units, reading=rd, measure=ms, flag=''), mean,
                                 E90.bootstrap_mean(X, N_BOOT, BOOT_SEED), int(X.shape[0]))
            x = X.copy()
            x.columns = [f'{a}|{b}' for a, b in X.columns]
            long.append(x.reset_index(names='recording').melt(id_vars='recording', var_name='pair', value_name='value')
                        .assign(genotype=gk, units=units, reading=rd, measure=ms))
    R = E90.remy_soma_fly_distances()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        mean = pd.Series(np.nanmean(R.to_numpy(dtype=float), axis=0), index=R.columns)
    rows += summary_rows(dict(genotype='remy', units='remy', reading='trial_mean', measure='pearson', flag=REMY_FLAG),
                         mean, E90.bootstrap_mean(R, N_BOOT, BOOT_SEED), int(R.shape[0]))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'somata_flies.csv', index=False)
    pd.concat(long, ignore_index=True).to_parquet(local_dir('A1') / 'flies_recordings.parquet')
    _FLIES_GROUP.clear()
    print(f'A1 in {time.time() - t0:.0f} s', flush=True)
    return df


_FLIES_GROUP: Dict[str, pd.DataFrame] = {}


def flies_group(units: str, genotype: str, reading: str, measure: str) -> pd.DataFrame:
    """The flies' group rows of somata_flies.csv: group -> mean, lo, hi."""
    if 'df' not in _FLIES_GROUP:
        _FLIES_GROUP['df'] = pd.read_csv(out_dir() / 'somata_flies.csv')
    d = _FLIES_GROUP['df']
    x = d[(d.kind == 'group') & (d.units == units) & (d.genotype == genotype) & (d.reading == reading) & (d.measure == measure)]
    if not len(x):
        raise KeyError((units, genotype, reading, measure))
    return x.set_index('group')[['mean', 'lo', 'hi']]


# ---------------------------------------------------------------------------
# the known-truth check

def _e90_dir(part: str) -> List[Path]:
    return [H.outputs_dir('E90_model_stages') / part, ARCHIVE / 'E90_model_stages' / part]


def e90_run(part: str, label: str) -> Dict[str, tuple]:
    """One of E90's archived runs ({'on', 'off'}: (responses, spike counts)): part 'B' (a wiring of C2's family), 'B5'
    ('real_pulse2') or 'D1' (a firing-ceiling fit)."""
    for d in _e90_dir(part):
        if (d / f'{label}_on_spike_counts.parquet').exists():
            return {c: tuple(pd.read_parquet(d / f'{label}_{cc}_{k}.parquet') for k in ('responses', 'spike_counts'))
                    for c, cc in (('on', 'on'), ('off', 'acute_off'))}
    raise FileNotFoundError(f'E90 {part} {label}')


def _zkw(cal: dict) -> dict:
    K = cal.get('K')
    return dict(q=float(cal['q']), noise=float(cal['noise']), amp=cal.get('amp', 'saturating'),
                K=None if K is None or not np.isfinite(K) else float(K), noise_ratio=float(cal['noise_ratio']),
                drift=float(cal['drift']))


def check_wiring(job: tuple) -> List[dict]:
    """One wiring's known-truth rows: per condition, level of CHECK_LEVELS, units and measure the group means of the
    estimate (`noise_removed`, mean over E89's pseudo-recordings), of the trial means (beside) and of the noise-free
    truth. job: (label, the flies' reliability, the flies' same-repeat gap; E89's targets)."""
    label, rel, gap = job
    warnings.simplefilter('ignore')
    from experiments import E88_responder_sets as E88
    from experiments import E89_apl_breadth_weighting as E89
    if tuple(E89.BINDING_LEVELS) != CHECK_LEVELS:
        raise RuntimeError(f"E89's binding levels {E89.BINDING_LEVELS} are not CHECK_LEVELS {CHECK_LEVELS}")
    t0 = time.time()
    o = e90_run('B', label)
    mr = E89.ModelReading89(o, label=label, rel_target=rel, gap_target=gap)
    kw = _zkw(mr.cal)
    rows = []
    for c in E88.CONDS:
        sr = H.session_zscore(mr.counts[c], **kw)['scale_ratio'][0]
        for level in CHECK_LEVELS:
            pr = mr.level_params(level, c)
            tr = mr.trials_l(c, pr['m'], pr['share'])
            for units in UNITS:
                f = np.ones_like(sr) if units == 'z' else sr
                truth_x = mr.C[c][0] * f[:, None]
                T_ = tr * f[None, None, :, None]
                for ms in MEASURES:
                    truth = group_means(_square_pairs(distance_matrix(truth_x, ms), ODORS16, PAIRS120))
                    est = pd.concat([noise_removed(T_[r], ms) for r in range(T_.shape[0])], axis=1).mean(axis=1)
                    tm = pd.concat([_square_pairs(distance_matrix(T_[r].mean(axis=0), ms), ODORS16, PAIRS120)
                                    for r in range(T_.shape[0])], axis=1).mean(axis=1)
                    ge, gt = group_means(est), group_means(tm)
                    for g in GROUPS:
                        rows.append(dict(wiring=label, cond=c, level=level, units=units, measure=ms, group=g,
                                         truth=truth[g], estimate=ge[g], trial_mean=gt[g], error=ge[g] - truth[g],
                                         error_trial_mean=gt[g] - truth[g], multiplier=float(pr['m']),
                                         share=float(pr['share'])))
    print(f'check {label}: {time.time() - t0:.0f} s', flush=True)
    return rows


def check_table(rows: pd.DataFrame) -> pd.DataFrame:
    """The check's rows with the population's (the mean over the wirings) added and `within` (|error| <= CHECK_TOL)."""
    key = ['cond', 'level', 'units', 'measure', 'group']
    pop = rows.groupby(key)[['truth', 'estimate', 'trial_mean', 'multiplier', 'share']].mean().reset_index()
    pop['error'] = pop.estimate - pop.truth
    pop['error_trial_mean'] = pop.trial_mean - pop.truth
    pop['wiring'] = 'population'
    df = pd.concat([rows, pop], ignore_index=True)
    df['within'] = df.error.abs() <= CHECK_TOL
    return df


def part_check(workers: int = 8) -> pd.DataFrame:
    """estimate_check.csv (per wiring and the population); the passes per group, measure and units in
    estimate_check.json."""
    from experiments import E88_responder_sets as E88
    from experiments import E89_apl_breadth_weighting as E89
    t0 = time.time()
    provenance('Acheck', CHECK_TOL=CHECK_TOL, CHECK_LEVELS=CHECK_LEVELS, wirings=ALL_WIRINGS, seed=E88.SEED,
               n_pseudo=E88.N_PSEUDO, n_trials=E88.N_TRIALS, workers=workers)
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    res = _map(check_wiring, [(w, rel, gap) for w in ALL_WIRINGS], workers)
    df = check_table(pd.DataFrame([r for rr in res for r in rr]))
    df.to_csv(out_dir() / 'estimate_check.csv', index=False)
    passes = check_passes(df)
    _write_json(out_dir() / 'estimate_check.json', dict(
        rule=f'a group passes in a measure and units if the population\'s (eight wirings) noise-removed estimate is within '
             f'{CHECK_TOL} of the model\'s noise-free group distance at every level of {list(CHECK_LEVELS)} in both '
             f'conditions', flies_reliability=rel, flies_gap=gap,
        passes={f'{g}|{m}|{u}': v for (g, m, u), v in passes.items()},
        largest_population_error={f'{m}|{u}': float(df[(df.wiring == 'population') & (df.measure == m) & (df.units == u)
                                                       & df.group.isin(TABLE_GROUPS)].error.abs().max())
                                  for m in MEASURES for u in UNITS}, seconds=round(time.time() - t0, 1)))
    print(f'check in {time.time() - t0:.0f} s; fails: {[k for k, v in passes.items() if not v]}', flush=True)
    return df


def check_passes(ec: Optional[pd.DataFrame] = None) -> Dict[Tuple[str, str, str], bool]:
    """{(group, measure, units): passes} from the check's population rows (estimate_check.csv by default)."""
    ec = pd.read_csv(out_dir() / 'estimate_check.csv') if ec is None else ec
    p = ec[ec.wiring == 'population']
    return {(g, m, u): bool(x.within.all()) for (g, m, u), x in p.groupby(['group', 'measure', 'units'])}


# ---------------------------------------------------------------------------
# A2: the model's somata

def model_calibrations(on: pd.DataFrame) -> Dict[str, dict]:
    """E73's calibration of the per-KC session z-score on one run's APL-on counts (KCs x the 16 odors), per map
    (MAPS): the noise so that the APL-on response rate under the data's criterion is the control flies', K of the
    saturating map so that the responders' mean amplitude is theirs, the noise ratio from their silent KCs' trial SD."""
    from experiments import E73_per_kc_normalization as E73
    t = E73.refs()['targets']
    q = E73.q_values()[0]
    drift = E73.refs()['drift']
    on = _abbrev16(on)
    return {m: H.calibrate_session_zscore(on, q=q, amp=m, target_rate=t['rate'],
                                          target_amp=t['amp'] if m == 'saturating' else None,
                                          target_silent_sd=t['silent_sd'], drift=drift) for m in MAPS}


def somata_pairs(on: pd.DataFrame, off: pd.DataFrame, cal: Optional[Dict[str, dict]] = None,
                 K: Optional[float] = None) -> pd.DataFrame:
    """One run's somata in every reading of A2 (the 66 megamat pairs): rows cond, reading, odor_a, odor_b, value;
    readings '<map>_<units>_<noise>_<measure>' (the per-KC session z-score: `model_calibrations` unless `cal` is given;
    units 'z' the reading's z, 'noise' the amplitude over the frame noise; noise 'none' or 'expected' through the
    reading's measurement noise of a three-trial mean) and E90's 'raw_<measure>' and 'cal_<measure>' (`E90.soma_rows`,
    K set on the APL-on counts unless given)."""
    from experiments import E83_cross_wiring_fit as E83
    on16, off16 = _abbrev16(on), _abbrev16(off)
    cal = model_calibrations(on16) if cal is None else cal
    parts = []
    for m in MAPS:
        kw = _zkw(dict(cal[m], amp=m))
        for cond, X in (('on', on16), ('off', off16)):
            z = H.session_zscore(X, **kw)
            C = z['center'][0]
            sd = z['trial_sd'][0] / np.sqrt(z['n_trials'])
            sr = z['scale_ratio'][0]
            for units in UNITS:
                V, s = (C, sd) if units == 'z' else (C * sr[:, None], sd * sr)
                for noise in NOISE:
                    for ms in MEASURES:
                        d = _square_pairs(distance_matrix(V, ms, None if noise == 'none' else s), ODORS16, PAIRS66)
                        parts.append(pd.DataFrame(dict(cond=cond, reading=f'{m}_{units}_{noise}_{ms}',
                                                       odor_a=d.index.get_level_values(0), odor_b=d.index.get_level_values(1),
                                                       value=d.to_numpy())))
    K = E83.readout_K({'on': (None, on)}) if K is None else K
    sr_ = pd.DataFrame(E90.soma_rows('run', {'on': on, 'off': off}, K))
    parts.append(sr_[['cond', 'reading', 'odor_a', 'odor_b', 'value']])
    return pd.concat(parts, ignore_index=True)


def _groups_of_pairs(df: pd.DataFrame, by: Sequence[str]) -> pd.DataFrame:
    """Group means (GROUPS) of pair rows (odor_a, odor_b, value) grouped by `by`."""
    out = []
    for k, sub in df.groupby(list(by), sort=False):
        s = pd.Series(sub.value.to_numpy(), index=pd.MultiIndex.from_arrays([sub.odor_a, sub.odor_b]))
        out.append(dict(zip(by, k if isinstance(k, tuple) else (k,)), **group_means(s)))
    return pd.DataFrame(out)


def somata_readings(on: pd.DataFrame, off: pd.DataFrame, cal: Optional[Dict[str, dict]] = None,
                    K: Optional[float] = None) -> pd.DataFrame:
    """One run's somata group distances in A2's readings (for experiments.E92_gain): `on` / `off` spike counts (KCs x
    odors, odor labels as olfsysm writes them or abbreviated); returns rows cond ('on' / 'off'), reading, then one column
    per group (GROUPS). `cal` (`model_calibrations`) and `K` (the calibrated readout's) are set on `on` unless given."""
    return _groups_of_pairs(somata_pairs(on, off, cal, K), ['cond', 'reading'])


def a2_runs() -> List[Tuple[str, str, str, str]]:
    """A2's runs: (E90's part, label, family, wiring)."""
    return ([('B', w, 'C2', w) for w in ALL_WIRINGS] + [('B5', 'real_pulse2', 'C2_pulse2', 'real')] +
            [('D1', m, 'fc', w) for m, w in D1_WIRING.items()])


def _a2_job(job: tuple) -> Tuple[pd.DataFrame, dict]:
    part, label, family, wiring = job
    warnings.simplefilter('ignore')
    t0 = time.time()
    o = e90_run(part, label)
    on, off = o['on'][1], o['off'][1]
    cal = model_calibrations(on)
    pairs = somata_pairs(on, off, cal)
    pairs.to_parquet(local_dir('A2') / f'{part}_{label}_pairs.parquet')
    info = dict(part=part, run=label, family=family, wiring=wiring, seconds=round(time.time() - t0, 1))
    for m, c in cal.items():
        info.update({f'{m}_{k}': c.get(k) for k in ('noise', 'K', 'noise_ratio', 'rate', 'resp_amp_mean',
                                                     'nonresp_trial_sd', 'noise_at_bound', 'K_at_bound', 'converged')})
    print(f'A2 {part} {label}: {info["seconds"]:.0f} s', flush=True)
    return pairs.assign(part=part, run=label, family=family, wiring=wiring), info


POPULATIONS = {'all8': ('C2', ALL_WIRINGS), 'S': ('C2', S_WIRINGS), 'T': ('C2', T_WIRINGS), 'fc_S4': ('fc', S_WIRINGS)}


def somata_model_table(pairs: pd.DataFrame) -> pd.DataFrame:
    """somata_model.csv's rows from A2's pair rows: per run and per population (POPULATIONS: the mean over its runs of
    the pair distances), condition and reading, the group means."""
    runs = _groups_of_pairs(pairs, ['part', 'run', 'family', 'wiring', 'cond', 'reading'])
    runs.insert(0, 'kind', 'run')
    pops = []
    for name, (fam, ws) in POPULATIONS.items():
        sub = pairs[(pairs.family == fam) & pairs.wiring.isin(ws)]
        if sub.run.nunique() != len(ws):
            continue
        p = sub.groupby(['cond', 'reading', 'odor_a', 'odor_b'], sort=False)['value'].mean().reset_index()
        g = _groups_of_pairs(p, ['cond', 'reading'])
        g.insert(0, 'kind', 'population')
        g.insert(1, 'part', '')
        g.insert(2, 'run', name)
        g.insert(3, 'family', fam)
        g.insert(4, 'wiring', name)
        pops.append(g)
    return pd.concat([runs] + pops, ignore_index=True)


def part_a2(workers: int = 8) -> pd.DataFrame:
    """somata_model.csv (per run and population, condition and reading, the group means) and
    somata_model_calibration.csv (E73's calibration of each run)."""
    t0 = time.time()
    provenance('A2', MAPS=MAPS, UNITS=UNITS, NOISE=NOISE, runs=a2_runs(), workers=workers)
    res = _map(_a2_job, a2_runs(), workers)
    pairs = pd.concat([r[0] for r in res], ignore_index=True)
    pairs.to_parquet(local_dir('A2') / 'somata_model_pairs.parquet')
    pd.DataFrame([r[1] for r in res]).to_csv(out_dir() / 'somata_model_calibration.csv', index=False)
    sm = somata_model_table(pairs)
    sm.to_csv(out_dir() / 'somata_model.csv', index=False)
    print(f'A2 in {time.time() - t0:.0f} s', flush=True)
    return sm


# ---------------------------------------------------------------------------
# A3: the verdict

def like_for_like() -> List[Tuple[str, str, str, str]]:
    """The like-for-like comparisons: (the model's reading prefix '<map>_<units>_<noise>', the flies' units, the flies'
    reading, measure): the per-KC session z-score in its z against Yang's z, in baseline-noise units against the flies'
    baseline-noise units; without noise against the flies' noise-removed values, through the noise against their trial
    means."""
    return [(f'{m}_{u}_{n}', u, 'noise_removed' if n == 'none' else 'trial_mean', ms)
            for m in MAPS for u in UNITS for n in NOISE for ms in MEASURES]


def compare(pop: Dict[str, float], wirings: pd.DataFrame, flies: pd.DataFrame, n_same: int = VERDICT_SAME_SIDE,
            groups: Sequence[str] = GROUPS) -> List[dict]:
    """E90's rule per group: the population's group means `pop` against the flies' (`flies`: group -> mean, lo, hi), the
    wirings' (rows of `wirings`, one column per group) on its side; `n_same` the number of wirings required on that
    side (E90's 6 of 8)."""
    rows = []
    for g in groups:
        r = E90.differs(float(pop[g]), wirings[g].to_numpy(dtype=float), float(flies.loc[g, 'mean']),
                        float(flies.loc[g, 'lo']), float(flies.loc[g, 'hi']))
        if r['n_same_side'] < n_same:
            r['same_side'] = False
            r['differs'] = False
        elif n_same != VERDICT_SAME_SIDE:
            r['same_side'] = True
            r['differs'] = bool(r['big'] and r['outside'])
        rows.append(dict(group=g, model=float(pop[g]), flies=float(flies.loc[g, 'mean']), lo=float(flies.loc[g, 'lo']),
                         hi=float(flies.loc[g, 'hi']), n_wirings=int(len(wirings)), **r))
    return rows


def _pop_and_wirings(sm: pd.DataFrame, reading: str, cond: str, pop: str = 'all8', family: str = 'C2',
                     wirings: Sequence[str] = ALL_WIRINGS) -> Tuple[pd.Series, pd.DataFrame]:
    x = sm[(sm.reading == reading) & (sm.cond == cond)]
    p = x[(x.kind == 'population') & (x.run == pop)]
    w = x[(x.kind == 'run') & (x.family == family) & x.wiring.isin(wirings)]
    if len(p) != 1 or len(w) != len(wirings):
        raise RuntimeError(f'{reading} {cond}: {len(p)} population rows, {len(w)} wiring rows')
    return p.iloc[0], w


def soma_verdict_rows(sm: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool],
                      flies: Callable[[str, str, str, str], pd.DataFrame] = None) -> pd.DataFrame:
    """Every (group, condition, like-for-like reading) comparison of the model's somata (the population of eight) with
    the flies' (`compare`); a noise-removed comparison of a group that failed the check is kept with read False. Beside
    (set 'beside'): raw counts and the calibrated readout against the flies in baseline-noise units."""
    flies = flies_group if flies is None else flies
    rows = []
    comps = [(p, u, fr, ms, 'like_for_like') for p, u, fr, ms in like_for_like()]
    comps += [(ro, 'noise', fr, ms, 'beside') for ro in ('raw', 'cal') for fr in FLY_READINGS for ms in MEASURES]
    for prefix, u, fr, ms, cset in comps:
        reading = f'{prefix}_{ms}'
        for cond in CONDS:
            pop, wir = _pop_and_wirings(sm, reading, cond)
            fl = flies(u, GENO[cond], fr, ms)
            for r in compare(pop, wir, fl):
                read = bool(fr != 'noise_removed' or passes.get((r['group'], ms, u), False))
                rows.append(dict(element='soma', set=cset, cond=cond, reading=reading, flies_units=u, flies_reading=fr,
                                 measure=ms, read=read, **r))
    return pd.DataFrame(rows)


def flies_claws() -> Dict[str, pd.DataFrame]:
    """Remy's claws per fly (flies x the 66 pairs, 1 - Pearson): 'noise_removed' (block averages) and 'measured' (E90's
    readers, over E90's STAGE_WINDOW)."""
    _, _, odors = E90.remy_matrices('claw', 'stim')
    return {'noise_removed': E90._flies_pairs(E90.remy_noise_removed('claw'), odors),
            'measured': E90._flies_pairs(E90.remy_window_means('claw'), odors)}


def model_claws() -> pd.DataFrame:
    """E90's model claws (all claws, the drive alone, APL on) per wiring: rows wiring, reading ('session_z' /
    'noisy_measured'), then one column per group (GROUPS)."""
    mp = None
    for d in (H.outputs_dir('E90_model_stages'), ARCHIVE / 'E90_model_stages'):
        if (d / 'model_pairs.parquet').exists():
            mp = pd.read_parquet(d / 'model_pairs.parquet')
            break
    if mp is None:
        raise FileNotFoundError("E90's model_pairs.parquet")
    c = mp[(mp.stage == 'claw') & (mp.subset == 'all') & (mp.share == 0.0) & (mp.cond == 'on')
           & mp.reading.isin(['session_z', 'noisy_measured'])]
    return _groups_of_pairs(c, ['wiring', 'reading'])


def flies_step_groups(soma: pd.DataFrame, claw: pd.DataFrame) -> pd.DataFrame:
    """The flies' step (somata minus claws) per group: mean, lo, hi (both stages' flies resampled, E90's `_step_boot`)."""
    soma = soma.reindex(columns=PAIRS66)
    claw = claw.reindex(columns=PAIRS66)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        st = pd.Series(np.nanmean(soma.to_numpy(dtype=float), axis=0) - np.nanmean(claw.to_numpy(dtype=float), axis=0),
                       index=PAIRS66)
    boot = E90._step_boot(soma, claw)
    gm, gb = group_means(st), _group_boot(boot, PAIRS66)
    return pd.DataFrame({g: dict(mean=gm[g], lo=float(np.nanpercentile(gb[g], CI[0])),
                                 hi=float(np.nanpercentile(gb[g], CI[1]))) for g in GROUPS}).T


def step_verdict_rows(sm: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool]) -> pd.DataFrame:
    """The two steps like for like (somata minus E90's claws on both sides, `CLAW_READING`), every like-for-like reading,
    by `compare`; rows as `soma_verdict_rows` with element 'step' (cond 'on' the step with APL on, 'off' the step without
    APL, both from the claws with APL on, as E90's). A noise-removed comparison is read if the group passed A1's check and
    E90's check of the claws' noise removal passed."""
    e90v = json.loads((H.results_dir('E90_model_stages') / 'verdict.json').read_text())
    claw_ok = bool(e90v['noise_check_passed']['claw'])
    mc = model_claws()
    fc = flies_claws()
    rows = []
    for prefix, u, fr, ms in like_for_like():
        noise = prefix.split('_')[-1]
        mcr, fcr = CLAW_READING[noise]
        reading = f'{prefix}_{ms}'
        for cond in CONDS:
            pop, wir = _pop_and_wirings(sm, reading, cond)
            cl = mc[mc.reading == mcr].set_index('wiring')
            w_step = wir.set_index('wiring')[list(GROUPS)] - cl.loc[list(ALL_WIRINGS), list(GROUPS)]
            p_step = {g: float(pop[g]) - float(cl.loc[list(ALL_WIRINGS), g].mean()) for g in GROUPS}
            fl = flies_step_groups(flies_samples(u)[(GENO[cond], fr, ms)], fc[fcr])
            for r in compare(p_step, w_step, fl):
                read = bool(fr != 'noise_removed' or (passes.get((r['group'], ms, u), False) and claw_ok))
                rows.append(dict(element='step', set='like_for_like', cond=cond, reading=reading, flies_units=u,
                                 flies_reading=fr, measure=ms, read=read, claws_model=mcr, claws_flies=fcr, **r))
    return pd.DataFrame(rows)


def element_result(rows: pd.DataFrame) -> dict:
    """E90's `element_result` over the rows that are read: 'differs' (every reading differs), 'mixed', 'n_differ',
    'n_readings', 'signs'."""
    return E90.element_result(rows[rows.read])


def condition_verdict(on: dict, off: dict) -> str:
    """'with APL on', 'without APL', 'both', 'neither', 'not resolved' (a condition's readings disagree) or 'not read'
    from the two conditions' element results (`element_result`)."""
    if on['n_readings'] == 0 or off['n_readings'] == 0:
        return 'not read'
    if on['mixed'] or off['mixed']:
        return 'not resolved'
    if on['differs'] and off['differs']:
        return 'both'
    if on['differs']:
        return 'with APL on'
    if off['differs']:
        return 'without APL'
    return 'neither'


def _breakdown(rows: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Per factor of the reading (measure, map or readout, units, the flies' reading), the readings that differ out of
    those read."""
    r = rows[rows.read].copy()
    if not len(r):
        return {}
    r['map'] = r.reading.str.split('_').str[0]
    out = {}
    for f in ('measure', 'map', 'flies_units', 'flies_reading'):
        out[f] = {str(k): f"{int(x.differs.sum())}/{len(x)}" + (f" ({'+' if x[x.differs]['diff'].mean() > 0 else '-'})"
                                                                if x.differs.any() else '')
                  for k, x in r.groupby(f)}
    return out


def verdicts(vr: pd.DataFrame) -> Dict[str, dict]:
    """Per group: each condition's element result and the condition verdict, with the signs and the breakdown."""
    out = {}
    for g in GROUPS:
        x = vr[vr.group == g]
        res = {c: element_result(x[x.cond == c]) for c in CONDS}
        out[g] = dict(verdict=condition_verdict(res['on'], res['off']), on=res['on'], off=res['off'],
                      mean_diff_on=float(x[(x.cond == 'on') & x.read]['diff'].mean()),
                      mean_diff_off=float(x[(x.cond == 'off') & x.read]['diff'].mean()),
                      breakdown={c: _breakdown(x[x.cond == c]) for c in CONDS})
    return out


def outcomes(soma: Dict[str, dict], step: Dict[str, dict]) -> dict:
    """Which of the entry's outcomes (a) to (c) the like-for-like verdicts support, on OUTCOME_GROUPS: (a) no group's
    somata differ without APL in any reading, the ester group's somata differ with APL on in every reading too alike
    (sign -) and the alcohol group's too distinct (+); (b) some group's somata or step without APL differ in every
    reading, too alike (the step without APL stays short); (c) some group's readings disagree in a condition; and whether
    the ester group's error with APL on has one sign in every reading that differs (part B's condition under (c))."""
    groups = OUTCOME_GROUPS
    off_clean = all(soma[g]['off']['n_differ'] == 0 for g in groups)
    est, alc = soma['ester']['on'], soma['alcohol']['on']
    a = bool(off_clean and est['differs'] and est['signs'] == [-1] and alc['differs'] and alc['signs'] == [1])
    b = bool(any((v[g]['off']['differs'] and v[g]['off']['signs'] == [-1]) for v in (soma, step) for g in groups))
    c = bool(any(v[g][cn]['mixed'] for v in (soma, step) for g in groups for cn in CONDS))
    return dict(a=a, b=b, c=c, somata_off_clean=off_clean, ester_on_signs=est['signs'],
                ester_on_one_sign=bool(len(est['signs']) <= 1 and est['n_differ'] > 0),
                ester_on=f"{est['n_differ']}/{est['n_readings']}", alcohol_on=f"{alc['n_differ']}/{alc['n_readings']}")


def beside_rows(sm: pd.DataFrame) -> pd.DataFrame:
    """For information, not in the verdict: D1's four firing-ceiling fits as a population of four (3 of 4 on the side),
    and B5 (the real wiring at the 2 s pulse) and the real wiring at 0.5 s alone (the difference from the flies only), in
    the like-for-like readings."""
    rows = []
    for prefix, u, fr, ms in like_for_like():
        reading = f'{prefix}_{ms}'
        for cond in CONDS:
            fl = flies_group(u, GENO[cond], fr, ms)
            x = sm[(sm.reading == reading) & (sm.cond == cond)]
            p = x[(x.kind == 'population') & (x.run == 'fc_S4')].iloc[0]
            w = x[(x.kind == 'run') & (x.family == 'fc')]
            for r in compare(p, w, fl, n_same=3):
                rows.append(dict(element='soma', set='fc_S4', cond=cond, reading=reading, flies_units=u, flies_reading=fr,
                                 measure=ms, **r))
            for run in ('real', 'real_pulse2'):
                rr = x[(x.kind == 'run') & (x.run == run)].iloc[0]
                for g in GROUPS:
                    rows.append(dict(element='soma', set=f'single_{run}', cond=cond, reading=reading, flies_units=u,
                                     flies_reading=fr, measure=ms, group=g, model=float(rr[g]), flies=float(fl.loc[g, 'mean']),
                                     lo=float(fl.loc[g, 'lo']), hi=float(fl.loc[g, 'hi']),
                                     diff=float(rr[g]) - float(fl.loc[g, 'mean'])))
    return pd.DataFrame(rows)


def part_a3() -> dict:
    """verdict_rows_A.csv (every comparison), steps_like_for_like.csv (the step rows), verdict_A.json: per group the
    somata's and the steps' verdicts with the signs and the breakdown by reading, over the readings that passed the
    check; which of the entry's outcomes (a) to (c) holds; beside, raw counts and the calibrated readout, D1, B5."""
    t0 = time.time()
    provenance('A3', VERDICT_DIFF=VERDICT_DIFF, VERDICT_SAME_SIDE=VERDICT_SAME_SIDE, CLAW_READING=CLAW_READING)
    sm = pd.read_csv(out_dir() / 'somata_model.csv')
    passes = check_passes()
    sv = soma_verdict_rows(sm, passes)
    st = step_verdict_rows(sm, passes)
    bs = beside_rows(sm)
    pd.concat([sv, st, bs], ignore_index=True).to_csv(out_dir() / 'verdict_rows_A.csv', index=False)
    st.to_csv(out_dir() / 'steps_like_for_like.csv', index=False)
    soma = verdicts(sv[sv.set == 'like_for_like'])
    step = verdicts(st)
    beside = verdicts(sv[sv.set == 'beside'].assign(read=True))
    by_map = {}
    for m in MAPS:
        sm_ = verdicts(sv[(sv.set == 'like_for_like') & sv.reading.str.startswith(m + '_')])
        st_ = verdicts(st[st.reading.str.startswith(m + '_')])
        by_map[m] = dict(somata={g: v['verdict'] for g, v in sm_.items()}, steps={g: v['verdict'] for g, v in st_.items()},
                         somata_on={g: v['on'] for g, v in sm_.items()}, somata_off={g: v['off'] for g, v in sm_.items()},
                         outcomes=outcomes(sm_, st_))
    out = dict(somata=soma, steps=step, beside_raw_cal=beside, outcomes=outcomes(soma, step), by_map=by_map,
               check_passes={f'{g}|{m}|{u}': v for (g, m, u), v in passes.items()},
               rules=dict(VERDICT_DIFF=VERDICT_DIFF, VERDICT_SAME_SIDE=VERDICT_SAME_SIDE, readings=like_for_like(),
                          statement='a condition differs only if every like-for-like reading read differs; it does not '
                                    'differ only if none does; otherwise not resolved'))
    _write_json(out_dir() / 'verdict_A.json', out)
    print(json.dumps({g: (soma[g]['verdict'], step[g]['verdict']) for g in TABLE_GROUPS}, indent=1), flush=True)
    print(json.dumps(out['outcomes'], indent=1), f'A3 in {time.time() - t0:.0f} s', flush=True)
    return out


# ---------------------------------------------------------------------------
# A4: the line

def relative_distances(d_on: pd.Series) -> Tuple[float, float]:
    """(ester, alcohol): the ester and the alcohol groups' mean APL-on distances, each minus the mean over all of the
    Series' pairs (the 120 pairs)."""
    s = E90._cs(d_on)
    m = group_masks(s.index)
    v = s.to_numpy(dtype=float)
    mu = float(np.nanmean(v))
    return float(np.nanmean(v[m['ester']])) - mu, float(np.nanmean(v[m['alcohol']])) - mu


def fit_line(est_rel: Sequence[float], alc_rel: Sequence[float]) -> dict:
    """The least-squares line of the ester on the alcohol relative distances: slope, intercept, r, n, the residuals' SD
    (`scatter`, n - 2 degrees of freedom) and largest absolute residual."""
    y, x = np.asarray(est_rel, dtype=float), np.asarray(alc_rel, dtype=float)
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    n = len(x)
    return dict(slope=float(b), intercept=float(a), r=float(np.corrcoef(x, y)[0, 1]), n=int(n),
                scatter=float(np.sqrt((res ** 2).sum() / (n - 2))), max_abs_resid=float(np.abs(res).max()),
                residuals=[float(v) for v in res])


def flies_points() -> Dict[str, Dict[str, float]]:
    """The control flies' points (ester, alcohol relative distances with APL on, cosine, Yang's z): per reading
    ('trial_mean', the main target's; 'noise_removed') and recordings ('all', half 'A', half 'B': E77's halves)."""
    from experiments import E77_input_side_common_error as E77
    halves = E77.recording_halves()
    S = flies_samples('z')
    out = {}
    for rd in FLY_READINGS:
        X = S[('ctrl', rd, 'cosine')]
        for k in FLY_POINTS:
            rows = list(X.index) if k == 'all' else [recording_label(r) for r in halves[k]['ctrl']]
            if not set(rows) <= set(X.index):
                raise RuntimeError(f'half {k}: recordings {sorted(set(rows) - set(X.index))} not in the flies\' tables')
            e, a = relative_distances(X.loc[rows].mean(axis=0))
            out[f'{rd}_{k}'] = dict(ester=e, alcohol=a)
    return out


def _with_flies(ln: dict, fp: Dict[str, Dict[str, float]]) -> dict:
    ln = dict(ln)
    ln['flies'] = {k: dict(v, resid=float(v['ester'] - (ln['intercept'] + ln['slope'] * v['alcohol'])))
                   for k, v in fp.items()}
    return ln


def place_on_line(est_rel: float, alc_rel: float, ln: dict, ref_resid: float = 0.0) -> dict:
    """One variant's place against a line (one of `line(readout)`'s, with its flies' points): its residual (minus
    `ref_resid`, the residual of the variant's reference where the line was not fitted on its wirings), the residual
    signed toward each flies' point's side of the line ('toward_<point>'), and whether it is off the line toward the flies
    by more than the scatter against half A and against half B, in the trial means ('off_toward_flies', the main
    target's reading) and noise removed ('off_toward_flies_nr'); 'beyond_every_weighting' if also beyond the largest
    residual of the line's own points."""
    r = float(est_rel - (ln['intercept'] + ln['slope'] * alc_rel)) - float(ref_resid)
    out = dict(resid=r)
    for k, v in ln['flies'].items():
        out[f'toward_{k}'] = r * float(np.sign(v['resid']))
    for suf, rd in (('', 'trial_mean'), ('_nr', 'noise_removed')):
        tw = [out[f'toward_{rd}_{h}'] for h in ('A', 'B')]
        out[f'off_toward_flies{suf}'] = bool(all(t > ln['scatter'] for t in tw))
        out[f'beyond_every_weighting{suf}'] = bool(all(t > ln['max_abs_resid'] for t in tw))
    return out


_LINE: Dict[str, dict] = {}


def line(readout: str = 'raw') -> dict:
    """A4's lines on `readout` ('raw' / 'cal'): {'S': the line through E89's 13 weightings (C1) as populations of S's four
    wirings, 'real' / 'P1' / 'P2' / 'M1': through the same weightings on that wiring alone}, each `fit_line`'s keys, the
    weightings' points ('points') and the control flies' points with their residuals ('flies'); from line.json once part
    A4 has run, else computed (which reads E89's 52 C1 runs)."""
    if readout not in _LINE:
        path = out_dir() / 'line.json'
        if path.exists():
            _LINE.update(json.loads(path.read_text())['lines'])
        else:
            _LINE.update(compute_lines(_line_runs(1))[0])
    return _LINE[readout]


def _variant_key(folder: str, tag: str, wiring_col) -> Tuple[str, str]:
    """(family, wiring) of an archived tag."""
    w = str(wiring_col) if isinstance(wiring_col, str) and wiring_col in ALL_WIRINGS else None
    if folder == 'E78_claw_knee':
        m = re.match(r'^(C2|P1)_(.*)$', tag)
        return (m.group(2), 'real' if m.group(1) == 'C2' else 'P1') if m else (tag, 'real')
    if folder == 'E82_flies_protocol':
        return re.sub(r'_(C2|(?:P1|P2|M1)_b2)$', '', tag), w or 'real'
    m = re.match(WIRING_RE, tag)
    if m:
        return m.group(1), m.group(2)
    return tag, w or 'real'


def _find_files(folder: str) -> Dict[str, List[Path]]:
    """{tag: [paths of <tag>_on_spike_counts.parquet]} under an archive folder, SKIP_DIRS left out."""
    out: Dict[str, List[Path]] = {}
    root = ARCHIVE / folder
    for p in root.rglob('*_on_spike_counts.parquet'):
        if any(s in part for part in p.relative_to(root).parts[:-1] for s in SKIP_DIRS):
            continue
        out.setdefault(p.name[:-len('_on_spike_counts.parquet')], []).append(p)
    return out


def a4_variants() -> pd.DataFrame:
    """The archived runs A4 reads: every tag of the folders' params.csv (A4_FOLDERS) with both conditions' counts and
    responses, and E90's eight runs of C2 (the references); folder, tag, family, wiring, kind, path of the APL-on counts,
    included (False with the reason for the other families' runs)."""
    rows = []
    for folder in A4_FOLDERS:
        p = pd.read_csv(ARCHIVE / folder / 'params.csv', low_memory=False).drop_duplicates('tag')
        files = _find_files(folder)
        wcol = p['wiring'] if 'wiring' in p.columns else pd.Series([None] * len(p), index=p.index)
        for tag, wc in zip(p.tag, wcol):
            fam, w = _variant_key(folder, tag, wc)
            reason = ''
            if folder in NOT_C2 and re.search(NOT_C2[folder], tag):
                reason = 'another family (firing ceiling)'
            paths = files.get(tag, [])
            path, off, off_from = None, None, ''
            for f in sorted(paths, key=lambda x: len(x.parts)):
                if not f.with_name(f'{tag}_on_responses.parquet').exists():
                    continue
                off, off_from = f.with_name(f'{tag}_acute_off_spike_counts.parquet'), 'own'
                if not off.exists():
                    off, off_from = _off_reference(folder, f, w)
                if off is not None and off.exists() and _resp(off).exists():
                    path = f
                    break
            if path is None and not reason:
                reason = 'no saved APL-on counts with an APL-off run'
            rows.append(dict(folder=folder, tag=tag, family=fam, wiring=w,
                             kind='oracle' if 'oracle' in tag else 'variant', path=str(path) if path else '',
                             off_path=str(off) if path else '', off_from=off_from if path else '', n_paths=len(paths),
                             included=not reason, reason=reason))
    for w in ALL_WIRINGS:
        d = next(x for x in _e90_dir('B') if (x / f'{w}_on_spike_counts.parquet').exists())
        rows.append(dict(folder='E90_model_stages', tag=w, family='C2_unchanged', wiring=w, kind='reference',
                         path=str(d / f'{w}_on_spike_counts.parquet'), off_path=str(d / f'{w}_acute_off_spike_counts.parquet'),
                         off_from='own', n_paths=1, included=True, reason=''))
    return pd.DataFrame(rows)


def _resp(counts_path: Path) -> Path:
    return counts_path.with_name(counts_path.name.replace('_spike_counts.parquet', '_responses.parquet'))


def _off_reference(folder: str, on_path: Path, wiring: str) -> Tuple[Optional[Path], str]:
    """The APL-off run of a variant saved with APL on only, as its experiment scored it: E87's static forms (D1) the
    wiring's own fit's (E90's run of it), E88's the base's with the flies' recruitment (base_<wiring>)."""
    if folder == 'E87_inhibition_form':
        d = next((x for x in _e90_dir('B') if (x / f'{wiring}_acute_off_spike_counts.parquet').exists()), None)
        return (d / f'{wiring}_acute_off_spike_counts.parquet', 'the wiring\'s fit (E90)') if d else (None, '')
    if folder == 'E88_responder_sets':
        return on_path.with_name(f'base_{wiring}_acute_off_spike_counts.parquet'), f'base_{wiring}'
    return None, ''


def _load_pair(path: str, off_path: Optional[str] = None) -> Dict[str, tuple]:
    p = Path(path)
    q = p.with_name(p.name.replace('_on_spike_counts.parquet', '_acute_off_spike_counts.parquet')) if not off_path \
        else Path(off_path)
    return {c: (pd.read_parquet(_resp(f)), pd.read_parquet(f)) for c, f in (('on', p), ('off', q))}


def _a4_job(job: tuple) -> dict:
    """One run: its realism (E83's `realism`), its 120 pair distances with APL on and off on raw counts and under the
    calibrated readout (as arrays in PAIRS120 order) and the readout's K."""
    idx, path, off_path = job
    warnings.simplefilter('ignore')
    from experiments import E83_cross_wiring_fit as E83
    try:
        o = _load_pair(path, off_path)
        out = dict(idx=idx, error='')
        rz = E83.realism(o)
        out.update({k: rz.get(k) for k in ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'opk_on_tv', 'opk_off_tv',
                                            'max_count_on', 'max_count_off', 'realism_flag', 'realism_flag_zcrit')})
        K = E83.readout_K(o)
        out['K'] = K
        for ro in READOUTS:
            d = E83.pair_distances(o, ro, K if ro == 'cal' else None)
            for c in CONDS:
                out[f'{ro}_{c}'] = E90._cs(d[c]).reindex(PAIRS120).to_numpy(dtype=float)
        return out
    except Exception as e:
        return dict(idx=idx, error=f'{type(e).__name__}: {e}')


def _line_runs(workers: int) -> Tuple[pd.DataFrame, List[dict]]:
    """E89's C1 runs (the line's 13 weightings on S's four wirings) as a4_variants rows and their `_a4_job` results."""
    runs = pd.read_csv(H.results_dir(LINE_EXPERIMENT) / f'runs_{LINE_PART}.csv')
    v = a4_variants()
    v = v[(v.folder == LINE_EXPERIMENT) & v.tag.isin(runs.tag)].reset_index(drop=True)
    res = _map(_a4_job, list(zip(v.index, v.path, v.off_path)), workers)
    return v.assign(weighting=v.tag.map(dict(zip(runs.tag, runs.weighting)))), res


def _pop(res_rows: Sequence[dict], key: str) -> pd.Series:
    return pd.Series(np.mean(np.stack([r[key] for r in res_rows]), axis=0), index=PAIRS120)


def compute_lines(line_runs: Tuple[pd.DataFrame, List[dict]]) -> Tuple[Dict[str, dict], pd.DataFrame]:
    """The lines per readout ({'raw' / 'cal': {'S', 'real', 'P1', 'P2', 'M1'}}) from E89's C1 runs, with the flies'
    points, and the weightings' points (a table)."""
    v, res = line_runs
    by = {r['idx']: r for r in res}
    if any(by[i]['error'] for i in v.index):
        raise RuntimeError([by[i]['error'] for i in v.index if by[i]['error']])
    fp = flies_points()
    lines, pts = {}, []
    for ro in READOUTS:
        lines[ro] = {}
        for key, ws in [('S', S_WIRINGS)] + [(w, (w,)) for w in S_WIRINGS]:
            e, a, names = [], [], []
            for wt, g in v.groupby('weighting'):
                rr = [by[i] for i in g.index if g.loc[i, 'wiring'] in ws]
                if len(rr) != len(ws):
                    raise RuntimeError(f'{wt}: {len(rr)} runs on {ws}')
                ee, aa = relative_distances(_pop(rr, f'{ro}_on'))
                e.append(ee)
                a.append(aa)
                names.append(wt)
                pts.append(dict(readout=ro, line=key, weighting=wt, ester=ee, alcohol=aa))
            if len(names) != 13:
                raise RuntimeError(f'{len(names)} weightings, not 13')
            ln = _with_flies(fit_line(e, a), fp)
            ln['points'] = dict(zip(names, zip(a, e)))
            lines[ro][key] = ln
    return lines, pd.DataFrame(pts)


def effect_parts(on: pd.Series, off: pd.Series, data: pd.Series) -> Dict[str, float]:
    """`H.apl_effect_stats` and `H.apl_effect_spread_shape` of a model's pair distances against the data's per-pair
    effect, every series first put in the harness's pair orientation (`H.canonical_pairs`: the harness aligns pairs by
    orientation and drops those that differ, and E83's / E90's series come in `_cs`'s sorted orientation); raises unless
    every pair of `on` with a finite distance in both conditions and in the data survives the alignment (a run with an
    odor silent in both conditions has NaN distances for its pairs)."""
    on, off, data = (H.canonical_pairs(x) for x in (on, off, data))
    if len(off) != len(on) or len(data) != len(on):
        raise RuntimeError(f'{len(on)}, {len(off)} and {len(data)} pairs')
    n_finite = int((np.isfinite(on.to_numpy(dtype=float)) & np.isfinite(off.reindex(on.index).to_numpy(dtype=float))
                    & np.isfinite(data.reindex(on.index).to_numpy(dtype=float))).sum())
    st = H.apl_effect_stats(on, off, data)
    if int(st['apl_effect_n_pairs']) != n_finite:
        raise RuntimeError(f"{int(st['apl_effect_n_pairs'])} of {n_finite} pairs survive the alignment")
    return dict(st, **H.apl_effect_spread_shape(on, off, data))


def variant_rows(v: pd.DataFrame, res: Dict[int, dict], lines: Dict[str, Dict[str, dict]]) -> pd.DataFrame:
    """line.csv's rows: per (folder, family) the populations of S and of T where all four wirings exist, single runs
    otherwise; per readout the relative distances, the place against the line (`place_on_line`: S's line for S
    populations and for T, there relative to C2's unchanged population on T or the wiring's run; the wiring's own line for
    single runs on S), the pattern part and the other scores under that readout against all recordings and the halves,
    the realism of the members."""
    from experiments import E83_cross_wiring_fit as E83
    tg = E83.targets()
    ok = v[v.included.to_numpy() & np.array([i in res and not res[i]['error'] for i in v.index])]
    units = []
    for (folder, fam), g in ok.groupby(['folder', 'family'], sort=False):
        ws = dict(zip(g.wiring, g.index))
        done = set()
        for name, need in (('S', S_WIRINGS), ('T', T_WIRINGS)):
            if set(need) <= set(ws):
                units.append((folder, fam, name, [ws[w] for w in need]))
                done |= set(need)
        for w, i in ws.items():
            if w not in done:
                units.append((folder, fam, w, [i]))
    ref = {}
    refs = ok[ok.family == 'C2_unchanged']
    for name, need in (('T', T_WIRINGS),) + tuple((w, (w,)) for w in T_WIRINGS):
        idx = [refs.index[refs.wiring == w][0] for w in need]
        for ro in READOUTS:
            e, a = relative_distances(_pop([res[i] for i in idx], f'{ro}_on'))
            ln = lines[ro]['S']
            ref[(ro, name)] = e - (ln['intercept'] + ln['slope'] * a)
    rows = []
    for folder, fam, name, idx in units:
        rr = [res[i] for i in idx]
        flags = [r['realism_flag'] or '' for r in rr]
        base = dict(folder=folder, family=fam, population=name, n_runs=len(idx), tags=';'.join(ok.loc[idx, 'tag']),
                    off_from=';'.join(sorted(set(ok.loc[idx, 'off_from']))),
                    kind='line' if folder == LINE_EXPERIMENT and fam.startswith(LINE_PART + '_') else ok.loc[idx[0], 'kind'],
                    in_bands=all(f == '' for f in flags),
                    realism_flags=' | '.join(sorted({f for f in flags if f})),
                    rate_on=float(np.mean([r['rate_on'] for r in rr])), rate_off=float(np.mean([r['rate_off'] for r in rr])),
                    rate_rho_on=float(np.mean([r['rate_rho_on'] for r in rr])),
                    max_count_on=int(max(r['max_count_on'] for r in rr)))
        for ro in READOUTS:
            d = {c: _pop(rr, f'{ro}_{c}') for c in CONDS}
            e, a = relative_distances(d['on'])
            if name in S_WIRINGS:
                ln, rref, basis = lines[ro][name], 0.0, f'the line on {name}'
            elif name == 'S':
                ln, rref, basis = lines[ro]['S'], 0.0, "S's line"
            else:
                ln, rref, basis = lines[ro]['S'], ref[(ro, name)], f"S's line, relative to C2 unchanged on {name}"
            pl = place_on_line(e, a, ln, rref)
            n_pairs = int(np.isfinite(d['on'].to_numpy()).sum())
            if n_pairs < len(PAIRS120):
                pl.update({k: False for k in pl if k.startswith(('off_toward', 'beyond'))})
            sc = E83.effect_scores(d, {k: tg[k] for k in ('z8', 'halfA', 'halfB')})
            ef = effect_parts(d['on'], d['off'], tg['z8'])
            rows.append(dict(base, readout=ro, n_pairs=n_pairs, ester_rel=e, alcohol_rel=a, basis=basis, **pl,
                             loss=sc['z8_loss'], pattern=sc['z8_pattern'], pattern_halfA=sc['halfA_pattern'],
                             pattern_halfB=sc['halfB_pattern'], mean_effect=sc['mean_effect'],
                             amplitude=ef['apl_effect_amplitude'], shape=ef['apl_effect_shape'],
                             sd_ratio=ef['apl_effect_sd_ratio'], pearson_r=ef['apl_effect_pearson'],
                             pattern_check=ef['apl_effect_rmse_centered'], rho=ef['apl_effect_spearman'],
                             err_ester=sc['z8_err_ester'], err_alcohol=sc['z8_err_alcohol']))
    return pd.DataFrame(rows)


def part_a4(workers: int = 8) -> pd.DataFrame:
    """line.json (the lines, the flies' points, the variants off the line), line.csv (every variant), line_points.csv
    (the weightings' points), a4_variants.csv (every archived run considered, with the reason for those left out)."""
    t0 = time.time()
    provenance('A4', A4_FOLDERS=A4_FOLDERS, NOT_C2=NOT_C2, SKIP_DIRS=SKIP_DIRS, LINE=(LINE_EXPERIMENT, LINE_PART),
               READOUTS=READOUTS, workers=workers)
    v = a4_variants()
    v.to_csv(out_dir() / 'a4_variants.csv', index=False)
    inc = v[v.included]
    print(f'A4: {len(inc)} runs of {len(v)} ({(~v.included).sum()} left out)', flush=True)
    out = _map(_a4_job, list(zip(inc.index, inc.path, inc.off_path)), workers)
    res = {r['idx']: r for r in out}
    errs = {int(k): r['error'] for k, r in res.items() if r['error']}
    lv = v[(v.folder == LINE_EXPERIMENT) & v.tag.isin(pd.read_csv(H.results_dir(LINE_EXPERIMENT) / f'runs_{LINE_PART}.csv').tag)]
    runs = pd.read_csv(H.results_dir(LINE_EXPERIMENT) / f'runs_{LINE_PART}.csv')
    lines, pts = compute_lines((lv.assign(weighting=lv.tag.map(dict(zip(runs.tag, runs.weighting)))),
                                [res[i] for i in lv.index]))
    pts.to_csv(out_dir() / 'line_points.csv', index=False)
    vr = variant_rows(v, res, lines)
    vr.to_csv(out_dir() / 'line.csv', index=False)
    off = vr[vr.in_bands & vr.off_toward_flies & (vr.kind != 'line')]
    named = vr[vr.in_bands & (vr.kind != 'line')]
    _write_json(out_dir() / 'line.json', dict(
        lines=lines, errors=errs, n_runs=int(len(inc)), n_variants=int(len(vr) // len(READOUTS)),
        n_in_bands=int(vr[vr.readout == 'raw'].in_bands.sum()),
        counts={ro: {f'{fo}|{k}': dict(n=int(len(x)), off_toward_flies=int(x.off_toward_flies.sum()),
                                       beyond_every_weighting=int(x.beyond_every_weighting.sum()))
                     for (fo, k), x in named[named.readout == ro].groupby(['folder', 'kind'])} for ro in READOUTS},
        off_toward_flies={ro: off[off.readout == ro][['folder', 'family', 'population', 'kind', 'resid',
                                                       'toward_trial_mean_A', 'toward_trial_mean_B',
                                                       'beyond_every_weighting', 'pattern', 'loss', 'mean_effect']]
                          .sort_values('resid', ascending=False).to_dict(orient='records') for ro in READOUTS},
        rule='off the line toward the flies: the residual (ester relative distance minus the line at the variant\'s '
             'alcohol relative distance), signed toward the side of the line the control flies\' point lies on, above '
             'the line\'s scatter (the residuals\' SD of its 13 weightings) with the flies\' point from half A and from '
             'half B (cosine trial means)', seconds=round(time.time() - t0, 1)))
    _LINE.clear()
    _LINE.update(lines)
    print(f'A4 in {time.time() - t0:.0f} s; errors {len(errs)}; off the line toward the flies: '
          f'{ {ro: len(off[off.readout == ro]) for ro in READOUTS} }', flush=True)
    return vr


# ---------------------------------------------------------------------------
# the table and the figures

def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


def comparison_table() -> None:
    """comparison_table.md: A1's flies, the check, A2's model, A3's verdict, A4's line; and the figures."""
    import subprocess
    o = out_dir()
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=H.HERE, capture_output=True, text=True).stdout.strip()
    G = TABLE_GROUPS
    L = ['# E92 part A: the step from claws to somata read like for like', '',
         f'Generated by experiments/E92_steps_like_for_like.py at {commit}. Distances between two odors\' KC responses '
         '(1 - Pearson r or 1 - cosine across KCs), group means over the 66 megamat pairs; the flies\' values with 90 % '
         f'intervals from {N_BOOT} resamples of the recordings. The model is C2\'s family on E72 part B\'s eight wirings; '
         '"population" is the mean over the eight of the pair distances.', '',
         '## A1: the flies\' somata', '',
         '| genotype | units | reading | measure | ' + ' | '.join(TABLE_HEAD[g] for g in G) + ' |',
         '|---|---|---|---|' + '---:|' * len(G)]
    fl = pd.read_csv(o / 'somata_flies.csv')
    fg = fl[fl.kind == 'group']
    lab = {'ctrl': 'control flies', 'tnt': 'APL output blocked (TNT)', 'remy': "Remy's wild type (as measured)"}
    for (gk, u, rd, ms), x in fg.groupby(['genotype', 'units', 'reading', 'measure'], sort=False):
        x = x.set_index('group')
        L.append(f'| {lab[gk]} | {u} | {rd} | {ms} | ' + ' | '.join(
            f"{_f(x.loc[g, 'mean'])} [{_f(x.loc[g, 'lo'], '.2f')}, {_f(x.loc[g, 'hi'], '.2f')}]" for g in G) + ' |')
    if (o / 'estimate_check.csv').exists():
        ec = pd.read_csv(o / 'estimate_check.csv')
        p = ec[(ec.wiring == 'population') & ec.group.isin(G)]
        L += ['', '## The known-truth check of the noise removal', '',
              f'Largest |estimate - truth| of the population over both conditions and the levels {list(CHECK_LEVELS)}; '
              f'tolerance {CHECK_TOL}. Trial means beside (their error is the inflation the estimate removes).', '',
              '| measure | units | ' + ' | '.join(TABLE_HEAD[g] for g in G) + ' | trial means, largest error |',
              '|---|---|' + '---:|' * (len(G) + 1)]
        for (ms, u), x in p.groupby(['measure', 'units']):
            L.append(f'| {ms} | {u} | ' + ' | '.join(
                f"{_f(x[x.group == g].error.abs().max(), '.4f')}{'' if x[x.group == g].within.all() else ' FAIL'}" for g in G)
                + f" | {_f(x.error_trial_mean.abs().max(), '.3f')} |")
    if (o / 'somata_model.csv').exists():
        sm = pd.read_csv(o / 'somata_model.csv')
        L += ['', '## A2: the model\'s somata (population of eight; the wirings\' range in parentheses)', '',
              '| APL | reading | ' + ' | '.join(TABLE_HEAD[g] for g in G) + ' |', '|---|---|' + '---:|' * len(G)]
        for cond in CONDS:
            for rd in sorted(sm.reading.unique()):
                x = sm[(sm.cond == cond) & (sm.reading == rd)]
                p = x[(x.kind == 'population') & (x.run == 'all8')]
                w = x[(x.kind == 'run') & (x.family == 'C2')]
                if len(p):
                    L.append(f'| {cond} | {rd} | ' + ' | '.join(
                        f"{_f(p.iloc[0][g])} ({_f(w[g].min(), '.2f')}-{_f(w[g].max(), '.2f')})" for g in G) + ' |')
    if (o / 'verdict_A.json').exists():
        v = json.loads((o / 'verdict_A.json').read_text())
        L += ['', '## A3: the verdict (like-for-like readings that passed the check)', '',
              '| group | somata: verdict | APL on: differ / read | APL off: differ / read | step: verdict | step on | '
              'step off |', '|---|---|---:|---:|---|---:|---:|']
        for g in G:
            s, t = v['somata'][g], v['steps'][g]

            def cell(e):
                sg = ''.join('+' if x > 0 else '-' for x in e['signs'])
                return f"{e['n_differ']}/{e['n_readings']}" + (f' ({sg})' if sg else '')
            L.append(f"| {GROUP_NAMES[g]} | {s['verdict']} | {cell(s['on'])} | {cell(s['off'])} | {t['verdict']} | "
                     f"{cell(t['on'])} | {cell(t['off'])} |")
        L += ['', f"Outcomes: {json.dumps(v['outcomes'])}", '']
    if (o / 'line.json').exists():
        lj = json.loads((o / 'line.json').read_text())
        L += ['## A4: the line', '', '| readout | line | slope | intercept | r | scatter | largest residual | control flies, '
              'trial means: residual (all / A / B) | noise removed |', '|---|---|---:|---:|---:|---:|---:|---|---|']
        for ro, lns in lj['lines'].items():
            for k, ln in lns.items():
                f_ = ln['flies']
                L.append(f"| {ro} | {k} | {_f(ln['slope'])} | {_f(ln['intercept'])} | {_f(ln['r'])} | {_f(ln['scatter'], '.4f')} | "
                         f"{_f(ln['max_abs_resid'], '.4f')} | " +
                         ' / '.join(_f(f_[f'trial_mean_{h}']['resid']) for h in FLY_POINTS) + ' | ' +
                         ' / '.join(_f(f_[f'noise_removed_{h}']['resid']) for h in FLY_POINTS) + ' |')
        lc = pd.read_csv(o / 'line.csv')
        off = lc[lc.in_bands & lc.off_toward_flies & (lc.kind != 'line')].sort_values(['readout', 'resid'], ascending=[True, False])
        L += ['', f"Variants: {lj['n_variants']} ({lj['n_in_bands']} inside the bands); off the line toward the flies "
              '(inside the bands):', '',
              '| readout | experiment | family | wirings | residual | toward the flies (A / B) | pattern part (cal) | loss | '
              'mean APL effect |', '|---|---|---|---|---:|---|---:|---:|---:|']
        for _, r in off.iterrows():
            L.append(f'| {r.readout} | {r.folder} | {r.family} | {r.population} | {_f(r.resid, "+.3f")} | '
                     f'{_f(r.toward_trial_mean_A, "+.3f")} / {_f(r.toward_trial_mean_B, "+.3f")} | {_f(r.pattern, ".4f")} | '
                     f'{_f(r.loss, ".4f")} | {_f(r.mean_effect, "+.3f")} |')
    (o / 'comparison_table.md').write_text('\n'.join(L) + '\n')
    figures()


def figures() -> None:
    """E92_somata.png (per group the flies' somata, trial means and noise removed, against the model's population in the
    per-KC readings, APL on and off) and E92_line.png (the line, its weightings, every variant and the flies' points)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    o = out_dir()
    # model readings in fixed categorical order (validated: dataviz validate_palette.js, light surface); markers differ too
    ink, ink2, grid = '#0b0b0b', '#52514e', '#e6e5e1'
    blue, magenta, purple, aqua = '#2a78d6', '#c2449f', '#5b4cc4', '#138a60'
    if (o / 'somata_model.csv').exists():
        sm = pd.read_csv(o / 'somata_model.csv')
        G = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near_ester', 'near_15ol_16ol', 'near_6al_t2h')
        fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharey=False)
        for i, cond in enumerate(CONDS):
            for j, ms in enumerate(MEASURES):
                ax = axes[i, j]
                x = np.arange(len(G))
                for k, (rd, mk, col, lab) in enumerate((('trial_mean', 'o', ink, 'flies, trial means'),
                                                        ('noise_removed', 's', ink2, 'flies, trial noise removed'))):
                    f = flies_group('z', GENO[cond], rd, ms)
                    ax.errorbar(x - 0.25 + 0.1 * k, [f.loc[g, 'mean'] for g in G],
                                yerr=[[f.loc[g, 'mean'] - f.loc[g, 'lo'] for g in G], [f.loc[g, 'hi'] - f.loc[g, 'mean'] for g in G]],
                                fmt=mk, color=col, ms=5, capsize=2, lw=1, label=lab, mfc='white' if k else col)
                for k, (rd, mk, col, lab) in enumerate((
                        (f'saturating_z_none_{ms}', 'D', blue, 'model, per-KC z, no noise'),
                        (f'saturating_z_expected_{ms}', 'd', magenta, 'model, per-KC z, through the noise'),
                        (f'raw_{ms}', 'v', purple, 'model, raw counts'),
                        (f'cal_{ms}', '^', aqua, 'model, calibrated readout'))):
                    p = sm[(sm.kind == 'population') & (sm.run == 'all8') & (sm.cond == cond) & (sm.reading == rd)]
                    if len(p):
                        ax.plot(x + 0.05 + 0.08 * k, [p.iloc[0][g] for g in G], mk, color=col, ms=5, label=lab)
                ax.set_xticks(x)
                ax.set_xticklabels([TABLE_HEAD[g] for g in G], fontsize=8, rotation=20)
                ax.set_title(f"{'APL on: control flies' if cond == 'on' else 'APL off: APL output blocked (TNT)'}, "
                             f"{'1 - Pearson r' if ms == 'pearson' else 'cosine distance'} (Yang's z)", fontsize=9)
                ax.grid(True, color=grid, lw=0.6)
                ax.spines[['top', 'right']].set_visible(False)
        axes[0, 0].set_ylabel('distance between two odors, group mean')
        axes[1, 0].set_ylabel('distance between two odors, group mean')
        axes[0, 0].legend(fontsize=7, frameon=False, ncol=2)
        fig.suptitle('KC somata by group: the flies (90 % intervals over recordings) and the model (population of eight)',
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(o / 'E92_somata.png', dpi=150)
        plt.close(fig)
    if (o / 'line.csv').exists():
        lc = pd.read_csv(o / 'line.csv')
        lj = json.loads((o / 'line.json').read_text())
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for ax, ro in zip(axes, READOUTS):
            ln = lj['lines'][ro]['S']
            x = lc[(lc.readout == ro) & (lc.population == 'S')]
            ax.scatter(x.alcohol_rel, x.ester_rel, s=14, c=np.where(x.in_bands, '#9aa7b4', '#e0d6c8'), label='S populations',
                       zorder=2)
            off = x[x.in_bands & x.off_toward_flies & (x.kind != 'line')]
            ax.scatter(off.alcohol_rel, off.ester_rel, s=30, facecolors='none', edgecolors=magenta, lw=1.2,
                       label='off the line toward the flies', zorder=3)
            pts = np.array(list(ln['points'].values()))
            ax.scatter(pts[:, 0], pts[:, 1], s=26, c=blue, label="E89's 13 weightings", zorder=4)
            xs = np.linspace(min(pts[:, 0].min(), x.alcohol_rel.min()) - 0.02, max(pts[:, 0].max(), x.alcohol_rel.max()) + 0.02, 10)
            ys = ln['intercept'] + ln['slope'] * xs
            ax.plot(xs, ys, '-', color=blue, lw=1)
            ax.fill_between(xs, ys - ln['scatter'], ys + ln['scatter'], color=blue, alpha=0.08, lw=0)
            for k, mk in (('trial_mean', '*'), ('noise_removed', 'P')):
                for h in FLY_POINTS:
                    f_ = ln['flies'][f'{k}_{h}']
                    ax.plot(f_['alcohol'], f_['ester'], mk, color=ink if h == 'all' else ink2, ms=11 if h == 'all' else 7,
                            label=f"control flies, {k.replace('_', ' ')}" if h == 'all' else None)
            ax.set_xlabel('alcohol pairs\' APL-on distance minus the mean pair\'s (cosine)', fontsize=9)
            ax.set_ylabel('ester pairs\' APL-on distance minus the mean pair\'s (cosine)', fontsize=9)
            ax.set_title(f"{'raw counts' if ro == 'raw' else 'calibrated readout'}: slope {ln['slope']:.2f}, "
                         f"scatter {ln['scatter']:.3f}", fontsize=10)
            ax.grid(True, color=grid, lw=0.6)
            ax.spines[['top', 'right']].set_visible(False)
        axes[0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(o / 'E92_line.png', dpi=150)
        plt.close(fig)


# ---------------------------------------------------------------------------
# B to F: implemented in experiments.E92_gain

def part_b(*args, **kwargs):
    """Part B, the gain at fixed scales on S: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_c(*args, **kwargs):
    """Part C, choice and confirmation: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_d(*args, **kwargs):
    """Part D, the family fitted again: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_e(*args, **kwargs):
    """Part E, on a pass only: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def part_f(*args, **kwargs):
    """Part F, the summed form on the same readings after E91's note: implemented in experiments.E92_gain."""
    raise NotImplementedError('implemented in experiments.E92_gain')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'A1', 'Acheck', 'A2', 'A3', 'A4', 'table', 'B', 'C', 'D', 'E', 'F'])
    ap.add_argument('--workers', type=int, default=8, help='Acheck / A2 / A4: worker processes')
    a = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    if a.part == 'A':
        part_a1()
        part_check(a.workers)
        part_a2(a.workers)
        part_a3()
        part_a4(a.workers)
        comparison_table()
        return
    {'A1': part_a1, 'Acheck': lambda: part_check(a.workers), 'A2': lambda: part_a2(a.workers), 'A3': part_a3,
     'A4': lambda: part_a4(a.workers), 'table': comparison_table, 'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e,
     'F': part_f}[a.part]()


if __name__ == '__main__':
    main()
