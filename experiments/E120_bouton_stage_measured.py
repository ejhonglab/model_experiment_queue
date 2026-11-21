#!/usr/bin/env python3
"""E120: the bouton stage held to the recorded boutons. E119's per-ROI bouton statistics read on saved runs (part A, no
simulations), the level of the antennal lobe's normalization for part B fixed by the entry's rule, and part B's helpers
(E102's part D at that level on the four searched wirings runs elsewhere; each refit then once with its stage matrices
saved, and E119's statistics read on it).

The entry is item E120 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from viki
2026-11-20 13:32; revised on moving 14:06). Reused, not duplicated:
- E119 (experiments/E119_stages_per_roi.py): the flies (part 0's per-ROI tables, `fly_roi_stats`, `fly_spread`; the
  excluded bouton fly stays excluded by part 0's lock), the model's ROIs in a reading (`rois_of_units`: pooling, silent
  ROIs dropped, 'rate', 'session_z', 'noisy'), the statistics (`roi_stats`, `matched_threshold`, `within`), the default's
  stage with the settled pre-odor level (`model_stage`) and its units' metadata (`unit_meta`).
- E90: `stage_units`, `subset_mask` (the boutons of the 28 shared glomeruli, E119's subset), `load_stage`, `read_stages`,
  `_save_stage`, al_transfer.csv (the ORN recordings' and PN dendrites' totals) and d2_pn_settings.csv.
- E98: the indicator family (INDICATOR_NS, FORMS; `indicator_ratio`'s two forms applied per unit and odor).
- E102: the ladder's saved PN stage (data/model_outputs/E102_antennal_lobe_normalization/B/<step>/), its costs
  (b_readings.csv), step tags and kwargs, part D's refit tags and parameters (`refit_tag`, `refit_params`).
- E94: the dose ladder's saved bouton stage (data/model_outputs/E94_bouton_share/B/) and its costs (equal.csv).
- E114: its inputs (results/E114_input_compression/inputs/), its PN and bouton stage tables when written.

Definitions (fixed before reading):
- The statistics are E119's: per wiring (or fly) the per-odor totals' CV over the 12 megamat odors (the mean change per
  odor over all ROIs, signed), the responders' CV, the top-1 % share, the share of entries below minus the threshold,
  at criteria 1, 2 and 3 baseline-noise SDs, the model's threshold matched to the flies' mean overall response rate at
  that criterion (on the run's APL-on ROIs where it has them). Boutons are not pooled (E119's verdict pool), claws in
  groups of 14 of one glomerulus. "Inside" is inside the flies' range (5 bouton flies, 4 claw flies). The decisions are
  read at E119's verdict criterion (2 SDs).
- The indicator (A1): each unit's window mean x1 and pre-odor level x0 (the settled level, as E119) read as
  x1^n - x0^n ('baseline'; rates clipped at 0) or sign(x1 - x0) |x1 - x0|^n ('evoked'), n in E98's INDICATOR_NS
  (1 to 2.5), at the default (E90's archive, the eight wirings, APL on). n 1 is the linear reading (E119's).
- A member passes A1 if, under one of E119's readings, the totals' CV, the responders' CV and the top-1 % share are all
  inside the flies' ranges on at least 6 of the 8 wirings; part B does not run if a member passes.
- The normalization (A2): the boutons with APL off are their PNs' rates (each bouton its glomerulus's PN; the model has
  no other input to a bouton with APL->bouton at 0), read from E102's saved PN stage of each step and wiring (window
  minus pre-odor mean over E90's windows; the settled pre-odor level differs from E90's by at most 0.018 Hz on the PNs).
  The threshold is matched on the APL-off ROIs themselves (no APL-on boutons are saved). E119's on-minus-off difference
  of each statistic at the default (per wiring) is given beside, with the APL-off value plus it as an estimate.
- The level for part B (the entry's rule): the smallest step from 100 to 10000 whose totals' CV ('rate' reading, APL off)
  is inside the flies' range on at least 6 of the 8 wirings; if none, 10000. The fixed gain is the bound, never a level.
  The other readings and the APL-on estimate are reported beside, and do not decide. If the level is 1000, E102's part D
  is the answer and part B does not run.
- The dose (A3): E94's saved settings (B1: the subtraction at exponent 2, doses 0 to 4; B2: at exponent 1; B3: the gain),
  APL on, its pre-odor mean as saved (E90's window); dose 1 is E90's archive of the same fits read the same way, and
  with the settled level beside. Four wirings: a setting is inside on a statistic at 3 of 4. E94's cost per setting
  (equal.csv) beside.
- A4 (E114's maps): the per-odor totals' CV of each input (E114 rescales each odor to its own sum, so the totals do not
  change, which the table checks), its PN stage (pn_stage.csv) and bouton stage (bouton files) when E114 has written
  them; E90's D2 exponents (no per-odor rescale) beside. Run again (`--part A4`) as E114 adds files.
- A5 (the bound): the totals' CV of the model's ORN input (36 and the 28 shared glomeruli), of the megamat ORN recordings
  and the PN dendrites (E90's al_transfer.csv, 28 shared glomeruli, linear and through the saturating map, 90 %
  intervals), of the model's PN stage (E102's identity step), its boutons and claws (E119), and of the flies' boutons and
  claws.

Outputs: <shared>/results/E120_bouton_stage_measured/ (a1_*.csv, a2_*.csv, a3_*.csv, a4_*.csv, a5_bound.csv,
level.json, outcome.json, comparison_table.md, provenance_<part>.json; after part B, b_stage_*_<step>.csv and
variants.csv, the four refits' parameters, the one table the archive's params.csv is built from); part B's stage runs
under outputs/E120_bouton_stage_measured/B/<step>/. outcome.json's `letter` is the entry's outcome (a) to (d); (e), the
dose, is `dose`; `part_b` holds E102's class of the refit at the level and the bouton totals' CV, once part B is read.

Usage (from the repository root; part A is light, no simulations):
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E120_bouton_stage_measured.py --part A
    ... --part A4 (again later); --part table
Part B (on hal, after E102's part D at the level; see `part_b_commands`):
    ... --part Bstage --level <L> --wirings real,P1 (16 threads); --part Bread --level <L>
"""
from __future__ import annotations

import argparse
import json
import time
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from experiments import E119_stages_per_roi as E119
from experiments import E90_model_stages as E90
from kcq import harness as H

EXPERIMENT = 'E120_bouton_stage_measured'
ODORS = tuple(E119.ODORS)                     # the 12 megamat odors
WIRINGS = tuple(E119.WIRINGS)                 # E72 part B's eight (S then T)
S_WIRINGS = ('real', 'P1', 'P2', 'M1')        # the four searched wirings (E102's S)
READINGS = tuple(E119.READINGS)
CRITERIA = tuple(E119.CRITERIA)
VERDICT_CRITERION = E119.VERDICT_CRITERION
MIN_IN = 6                                    # of the 8 wirings
MIN_IN_S = 3                                  # of the 4 searched wirings (A3)
THREE = ('totals_cv', 'resp_cv', 'top1')      # the three bouton statistics of the A1 rule
STATS = THREE + ('inhibited',)
CLAW_STATS = STATS + ('per_odor', 'opr')
POOL = {'bouton': '1', 'claw': 'g14'}         # E119's verdict pools
SUBSET = {'bouton': 'shared', 'claw': 'all'}  # E119's (E90's verdict subsets)
LEVEL_STEPS = (100.0, 314.0, 1000.0, 3141.0, 10000.0)
LEVEL_FALLBACK = 10000.0
E102_STEPS = ('identity', 100.0, 314.0, 1000.0, 3141.0, 10000.0, 1e6)
E102_DONE_LEVEL = 1000.0                      # E102's part D ran here
E94_WIRINGS = S_WIRINGS


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive(experiment: str) -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / experiment


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def provenance(part: str, **settings) -> None:
    _write_json(out_dir() / f'provenance_{part}.json', dict(part=part, script=H.script_path(__file__), **H.code_commits(),
                generated=time.strftime('%Y-%m-%d %H:%M:%S %Z'), settings=settings))


def _cv(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.std(x) / np.mean(x))


# ---------------------------------------------------------------------------
# the flies and the reading of a model's units

@lru_cache(maxsize=None)
def fly_stats(stage: str) -> Dict[float, Dict[str, dict]]:
    return E119.fly_roi_stats(stage)


@lru_cache(maxsize=None)
def fly_spreads(stage: str) -> Dict[float, dict]:
    return {c: E119.fly_spread(fly_stats(stage)[c]) for c in CRITERIA}


@lru_cache(maxsize=None)
def fly_reliability(stage: str) -> pd.Series:
    return E119.fly_reliability(stage)


def indicator_change(win: pd.DataFrame, pre: pd.DataFrame, n: float, form: str) -> pd.DataFrame:
    """A stage's units x odors change seen through an indicator of exponent n (E98's family): 'baseline'
    max(win, 0)^n - max(pre, 0)^n, 'evoked' sign(win - pre) |win - pre|^n. n 1 is the change itself (rates are >= 0)."""
    if form == 'baseline':
        x = np.clip(win.to_numpy(dtype=float), 0, None) ** n - np.clip(pre.to_numpy(dtype=float), 0, None) ** n
    elif form == 'evoked':
        d = win.to_numpy(dtype=float) - pre.to_numpy(dtype=float)
        x = np.sign(d) * np.abs(d) ** n
    else:
        raise ValueError(f'unknown form {form!r}')
    out = pd.DataFrame(x, index=win.index, columns=[H.odor_abbrev(str(c)) for c in win.columns])
    return out


def stage_subset(U: pd.DataFrame, stage: str, wiring: str) -> Tuple[pd.DataFrame, np.ndarray]:
    """The units of E119's subset (bouton 'shared', claw 'all') and their glomeruli (the wiring's E90 metadata; the
    units must be in its order)."""
    meta = E119.unit_meta(wiring)
    gl = np.asarray(meta['claw_glomerulus'] if stage == 'claw' else meta['bouton_glomerulus'])
    if len(gl) != U.shape[0]:
        raise RuntimeError(f'{wiring} {stage}: {U.shape[0]} units, {len(gl)} in the metadata')
    m = E90.subset_mask(stage, U.index, meta, SUBSET[stage])
    return U.loc[m], gl[m]


def wiring_reading(U: pd.DataFrame, groups: np.ndarray, stage: str, reading: str, wiring: str) -> pd.DataFrame:
    """One wiring's ROIs x ODORS in a reading (E119's `rois_of_units`, trial mean)."""
    wi = WIRINGS.index(wiring)
    T = E119.rois_of_units(U, groups, stage, reading, POOL[stage], wi, fly_reliability(stage))
    return pd.DataFrame(T.mean(axis=0), columns=list(ODORS))


def stats_of(V: pd.DataFrame, V_match: pd.DataFrame, stage: str) -> List[dict]:
    """E119's statistics of ROIs x ODORS V at each criterion, the threshold matched on V_match to the flies' mean rate;
    each with its sides against the flies' range (`E119.within`)."""
    out = []
    for c in CRITERIA:
        sp = fly_spreads(stage)[c]
        thr = E119.matched_threshold(V_match, sp['rate']['mean'])
        st = E119.roi_stats(V, thr)
        w = E119.within(st, sp)
        out.append(dict(criterion=c, threshold=thr, n_rois=st['n_rois'], **{s: st[s] for s in E119.SCALARS},
                        **{k: v for k, v in w.items()}))
    return out


def read_units(units: Mapping[str, Tuple[pd.DataFrame, np.ndarray]], stage: str, readings: Sequence[str] = READINGS,
               match_on: Optional[Mapping[str, Tuple[pd.DataFrame, np.ndarray]]] = None, **keys) -> pd.DataFrame:
    """E119's per-ROI statistics of each wiring's units ({wiring: (units x panel odors, glomeruli)}), in each reading and
    criterion; the threshold matched on `match_on`'s units of the same wiring (default the units themselves). Rows carry
    `keys`."""
    rows = []
    for rd in readings:
        for w, (U, g) in units.items():
            V = wiring_reading(U, g, stage, rd, w)
            Vm = V if match_on is None else wiring_reading(match_on[w][0], match_on[w][1], stage, rd, w)
            for r in stats_of(V, Vm, stage):
                rows.append(dict(**keys, stage=stage, reading=rd, wiring=w, **r))
    return pd.DataFrame(rows)


def fly_rows(stage: str) -> pd.DataFrame:
    """The flies' own statistics per criterion (rows wiring 'fly:<acq>')."""
    rows = []
    for c in CRITERIA:
        for acq, st in fly_stats(stage)[c].items():
            rows.append(dict(stage=stage, criterion=c, wiring=f'fly:{acq}', **{s: st[s] for s in E119.SCALARS}))
    return pd.DataFrame(rows)


def spread_table(stage: str) -> pd.DataFrame:
    """The flies' range of each statistic per criterion (lo, hi, mean) and the profile bands."""
    rows = []
    for c in CRITERIA:
        sp = fly_spreads(stage)[c]
        rows.append(dict(stage=stage, criterion=c, **{f'{s}_{k}': sp[s][k] for s in E119.SCALARS for k in ('lo', 'hi', 'mean')},
                         per_odor_band=sp['per_odor']['band'], opr_band=sp['opr']['band']))
    return pd.DataFrame(rows)


def summarize(rows: pd.DataFrame, by: Sequence[str], stats: Sequence[str] = STATS) -> pd.DataFrame:
    """Per group: each statistic's mean over the wirings and the number of wirings below / inside / above the flies'
    range (profiles: inside the band)."""
    out = []
    for k, sub in rows.groupby(list(by), dropna=False, sort=False):
        r = dict(zip(by, k if isinstance(k, tuple) else (k,)), n_wirings=int(len(sub)))
        for s in stats:
            if s in ('per_odor', 'opr'):
                r[f'{s}_dist'] = float(sub[f'{s}_dist'].mean())
                r[f'{s}_n_in'] = int(sub[f'{s}_in'].sum())
                continue
            side = sub[f'{s}_side']
            r[s] = float(sub[s].mean())
            r[f'{s}_lo_w'], r[f'{s}_hi_w'] = float(sub[s].min()), float(sub[s].max())
            r[f'{s}_below'], r[f'{s}_in'], r[f'{s}_above'] = int((side == -1).sum()), int((side == 0).sum()), int((side == 1).sum())
        out.append(r)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# A1: the indicator family at the default

def indicator_members() -> List[Tuple[float, str]]:
    from experiments import E98_apl_response_strength as E98
    return [(float(n), f) for f in E98.FORMS for n in E98.INDICATOR_NS]


def default_units(stage: str, wiring: str, n: float = 1.0, form: str = 'evoked', cond: str = 'on'
                  ) -> Tuple[pd.DataFrame, np.ndarray]:
    """The default family's stage (E90's archive, the settled pre-odor level, E119's `model_stage`) through the
    indicator, E119's subset."""
    st = E119.model_stage(wiring, cond)
    U = indicator_change(st[f'{stage}_win'], st[f'{stage}_pre'], n, form)
    return stage_subset(U, stage, wiring)


def part_a1(wirings: Sequence[str] = WIRINGS) -> Dict[str, pd.DataFrame]:
    """A1: every member of the indicator family on the default's boutons and claws (APL on, eight wirings). Writes
    a1_rows.csv (per member, reading, wiring, criterion), a1_summary.csv and a1_pass.csv (the rule, per member and
    reading)."""
    t0 = time.time()
    rows = []
    for n, f in indicator_members():
        for stage in ('bouton', 'claw'):
            units = {w: default_units(stage, w, n, f) for w in wirings}
            rows.append(read_units(units, stage, n=n, form=f))
        print(f'A1 n {n:g} {f}: {time.time() - t0:.0f} s', flush=True)
    rows = pd.concat(rows, ignore_index=True)
    rows.to_csv(out_dir() / 'a1_rows.csv', index=False)
    summ = summarize(rows, ('stage', 'n', 'form', 'reading', 'criterion'), CLAW_STATS)
    summ.to_csv(out_dir() / 'a1_summary.csv', index=False)
    ps = a1_pass(summ)
    ps.to_csv(out_dir() / 'a1_pass.csv', index=False)
    provenance('A1', members=indicator_members(), wirings=list(wirings), POOL=POOL, SUBSET=SUBSET, CRITERIA=CRITERIA,
               seconds=round(time.time() - t0, 1))
    return dict(rows=rows, summary=summ, passes=ps)


def a1_pass(summ: pd.DataFrame, min_in: int = MIN_IN) -> pd.DataFrame:
    """The A1 rule per member and reading (bouton, VERDICT_CRITERION): all of THREE inside on `min_in` wirings."""
    b = summ[(summ.stage == 'bouton') & (summ.criterion == VERDICT_CRITERION)].copy()
    b['passes'] = np.all([b[f'{s}_in'] >= min_in for s in THREE], axis=0)
    return b[['n', 'form', 'reading', *[f'{s}_in' for s in THREE], *THREE, 'passes']].reset_index(drop=True)


# ---------------------------------------------------------------------------
# A2: the normalization ladder (E102's saved PN stage)

def e102_step_tag(step) -> str:
    from experiments import E102_antennal_lobe_normalization as E102
    return E102.step_tag(step)


def e102_pn(step, wiring: str, cond: str = 'on') -> pd.DataFrame:
    """E102's saved PN stage of a step and wiring: glomeruli x the panel's 16 odors, window minus pre-odor mean."""
    tag = e102_step_tag(step)
    d = archive('E102_antennal_lobe_normalization') / 'B' / tag
    w = pd.read_parquet(d / f'B_{tag}_{wiring}_{cond}_pn_win.parquet')
    p = pd.read_parquet(d / f'B_{tag}_{wiring}_{cond}_pn_pre.parquet')
    cols = list(H.APL16_ODORS)
    return (w - p)[cols]


def pn_boutons(pn: pd.DataFrame, wiring: str) -> Tuple[pd.DataFrame, np.ndarray]:
    """The boutons with APL off as their PNs' changes (each bouton its glomerulus's row), E119's subset."""
    meta = E119.unit_meta(wiring)
    gl = np.asarray(meta['bouton_glomerulus'])
    missing = sorted(set(gl) - set(pn.index))
    if missing:
        raise RuntimeError(f'{wiring}: boutons of glomeruli without a PN row: {missing[:5]}')
    U = pd.DataFrame(pn.loc[gl].to_numpy(dtype=float), columns=list(pn.columns))
    m = E90.subset_mask('bouton', U.index, meta, SUBSET['bouton'])
    return U.loc[m].reset_index(drop=True), gl[m]


def e102_costs() -> pd.DataFrame:
    """E102's cost per step at fixed scales (b_readings.csv): the population of the eight wirings and of S, under the
    calibrated readout (pattern part against all recordings, its change from the default, the amplitude and shape parts,
    the spread ratio, the mean APL effect) and on raw counts (the loss's change, the mean APL effect)."""
    from experiments import E102_antennal_lobe_normalization as E102
    b = pd.read_csv(E102.out_dir() / 'b_readings.csv', comment='#')
    rows = []
    for unit in ('all8', 'S'):
        cal = b[(b.unit == unit) & (b.readout == 'cal')].set_index('change')
        raw = b[(b.unit == unit) & (b.readout == 'raw')].set_index('change')
        for tag in cal.index:
            rows.append(dict(step=tag, unit=unit, pattern_cal=cal.loc[tag, 'z8_pattern'], d_pattern_cal=cal.loc[tag, 'd_z8_pattern'],
                             amplitude_cal=cal.loc[tag, 'z8_amplitude'], shape_cal=cal.loc[tag, 'z8_shape'],
                             sd_ratio_cal=cal.loc[tag, 'z8_sd_ratio'], mean_effect_cal=cal.loc[tag, 'mean_effect'],
                             d_loss_raw=raw.loc[tag, 'd_z8_loss'], mean_effect_raw=raw.loc[tag, 'mean_effect']))
    return pd.DataFrame(rows)


def e119_on_minus_off() -> pd.DataFrame:
    """E119's difference between APL on and off at the default, per wiring, reading and criterion (a2_roi_stats.csv,
    boutons not pooled)."""
    a = pd.read_csv(E119.out_dir() / 'a2_roi_stats.csv')
    a = a[(a.pool.astype(str) == '1') & (a.cond.isin(['on', 'off']))]
    p = a.pivot_table(index=['reading', 'criterion', 'wiring'], columns='cond', values=list(STATS))
    out = pd.DataFrame({s: p[(s, 'on')] - p[(s, 'off')] for s in STATS}).reset_index()
    return out


def part_a2(wirings: Sequence[str] = WIRINGS, steps: Sequence = E102_STEPS) -> Dict[str, pd.DataFrame]:
    """A2: E102's steps on the eight wirings, the boutons with APL off from the saved PN stage; the identity check against
    E119's APL-off boutons; E119's on-minus-off beside, the APL-off value plus it as an estimate; E102's cost per step.
    Writes a2_rows.csv, a2_summary.csv, a2_identity.csv, a2_costs.csv."""
    t0 = time.time()
    rows = []
    ident = []
    for step in steps:
        tag = e102_step_tag(step)
        units = {w: pn_boutons(e102_pn(step, w), w) for w in wirings}
        if tag == 'identity':
            for w in wirings:
                ref = E119.model_units('bouton', w, 0.0)['off']
                a = units[w][0][list(H.APL16_ODORS)].to_numpy(dtype=float)
                b = ref[list(H.APL16_ODORS)].to_numpy(dtype=float)
                ident.append(dict(wiring=w, n_boutons=int(a.shape[0]), n_ref=int(b.shape[0]),
                                  max_abs_diff=float(np.abs(a - b).max()) if a.shape == b.shape else np.nan,
                                  max_abs_ref=float(np.abs(b).max())))
        r = read_units(units, 'bouton', step=tag, pn_inhadd=np.nan if tag == 'identity' else float(step), cond='off')
        rows.append(r)
        print(f'A2 {tag}: {time.time() - t0:.0f} s', flush=True)
    rows = pd.concat(rows, ignore_index=True)
    om = e119_on_minus_off().rename(columns={s: f'd_on_{s}' for s in STATS})
    rows = rows.merge(om, on=['reading', 'criterion', 'wiring'], how='left')
    for s in STATS:
        rows[f'{s}_on_est'] = rows[s] + rows[f'd_on_{s}']
    # sides of the estimate against the flies' range
    for s in STATS:
        lo = rows.criterion.map({c: fly_spreads('bouton')[c][s]['lo'] for c in CRITERIA})
        hi = rows.criterion.map({c: fly_spreads('bouton')[c][s]['hi'] for c in CRITERIA})
        v = rows[f'{s}_on_est']
        rows[f'{s}_on_est_side'] = np.where(v < lo, -1, np.where(v > hi, 1, 0))
    rows.to_csv(out_dir() / 'a2_rows.csv', index=False)
    summ = summarize(rows, ('step', 'reading', 'criterion'))
    est = rows.groupby(['step', 'reading', 'criterion'], sort=False).agg(
        **{f'{s}_on_est': (f'{s}_on_est', 'mean') for s in STATS},
        **{f'{s}_on_est_in': (f'{s}_on_est_side', lambda x: int((x == 0).sum())) for s in STATS}).reset_index()
    summ = summ.merge(est, on=['step', 'reading', 'criterion'], how='left')
    summ.to_csv(out_dir() / 'a2_summary.csv', index=False)
    pd.DataFrame(ident).to_csv(out_dir() / 'a2_identity.csv', index=False)
    costs = e102_costs()
    costs.to_csv(out_dir() / 'a2_costs.csv', index=False)
    provenance('A2', steps=[e102_step_tag(s) for s in steps], wirings=list(wirings), seconds=round(time.time() - t0, 1),
               source=str(archive('E102_antennal_lobe_normalization') / 'B'))
    return dict(rows=rows, summary=summ, identity=pd.DataFrame(ident), costs=costs)


def choose_level(summ: pd.DataFrame, min_in: int = MIN_IN) -> dict:
    """Part B's level by the entry's rule from a2_summary (the totals' CV, 'rate' reading, APL off; it does not depend on
    the criterion, read at VERDICT_CRITERION): the smallest LEVEL_STEPS step inside the flies' range on `min_in` of the
    8 wirings, else LEVEL_FALLBACK; whether part B runs (not at E102's done level). The other readings and the APL-on
    estimate beside."""
    s = summ[summ.criterion == VERDICT_CRITERION].copy()
    s['step'] = s['step'].astype(str)

    def first_inside(reading: str, col: str) -> Optional[float]:
        for v in LEVEL_STEPS:
            r = s[(s.step == e102_step_tag(v)) & (s.reading == reading)]
            if len(r) and int(r[col].iloc[0]) >= min_in:
                return float(v)
        return None
    found = first_inside('rate', 'totals_cv_in')
    level = found if found is not None else LEVEL_FALLBACK
    beside = {f'{rd}_{k}': first_inside(rd, col) for rd in READINGS for k, col in (('off', 'totals_cv_in'),
                                                                                   ('on_est', 'totals_cv_on_est_in'))}
    return dict(level=level, level_tag=e102_step_tag(level), found_inside=found is not None,
                part_b_runs=bool(level != E102_DONE_LEVEL), rule=(
                    "the smallest pn_inhadd of 100, 314, 1000, 3141, 10000 whose per-odor totals' CV (E119's 'rate' "
                    "reading, boutons with APL off from E102's saved PN stage) is inside the flies' range on 6 of 8 "
                    "wirings; 10000 if none; part B does not run at 1000 (E102's part D)"),
                beside_first_inside=beside)


# ---------------------------------------------------------------------------
# A3: the dose ladder (E94's saved bouton stage)

def e94_variants() -> List[str]:
    d = archive('E94_bouton_share') / 'B'
    vs = sorted({p.name.split('_real_on_bouton_win')[0] for p in d.glob('B_*_real_on_bouton_win.parquet')})
    return vs


def e94_units(variant: str, wiring: str) -> Tuple[pd.DataFrame, np.ndarray]:
    """E94's saved boutons of a setting and wiring, APL on (its pre-odor mean as saved), E119's subset; 'ref' is E90's
    archive of the fits (dose 1), read the same way, and 'ref_settled' with the settled level (E119's reading)."""
    if variant == 'ref':
        st = E90.load_stage('B', f'{wiring}_on')
        w, p = st['bouton_win'], st['bouton_pre']
    elif variant == 'ref_settled':
        st = E119.model_stage(wiring, 'on')
        w, p = st['bouton_win'], st['bouton_pre']
    else:
        d = archive('E94_bouton_share') / 'B'
        w = pd.read_parquet(d / f'{variant}_{wiring}_on_bouton_win.parquet')
        p = pd.read_parquet(d / f'{variant}_{wiring}_on_bouton_pre.parquet')
    U = indicator_change(w, p, 1.0, 'evoked')
    return stage_subset(U, 'bouton', wiring)


def e94_costs() -> pd.DataFrame:
    from experiments import E94_bouton_share as E94
    e = pd.read_csv(E94.out_dir() / 'equal.csv', comment='#')
    keep = ['variant', 'screen', 'form', 'exponent', 'dose', 'd_pattern_z8', 'd_pattern_halfA', 'd_pattern_halfB',
            'd_raw_loss', 'd_mean_effect', 'n_flagged', 'equal', 'better']
    return e[[c for c in keep if c in e.columns]]


def dose_of(variant: str) -> Tuple[str, float, float]:
    if variant.startswith('ref'):
        return 'sub', 2.0, 1.0
    from experiments import E94_bouton_share as E94
    pv = E94.parse_variant(variant)
    return pv['form'], float(pv['exponent']), float(pv['dose'])


def part_a3(wirings: Sequence[str] = E94_WIRINGS) -> Dict[str, pd.DataFrame]:
    """A3: E94's saved settings on the four searched wirings, APL on; dose 1 from E90's archive (as saved and settled).
    Writes a3_rows.csv, a3_summary.csv (with E94's cost), a3_dose.json (outcome (e))."""
    t0 = time.time()
    rows = []
    for v in ['ref', 'ref_settled', *e94_variants()]:
        form, k, dose = dose_of(v)
        units = {w: e94_units(v, w) for w in wirings}
        rows.append(read_units(units, 'bouton', variant=v, form=form, exponent=k, dose=dose))
    rows = pd.concat(rows, ignore_index=True)
    rows.to_csv(out_dir() / 'a3_rows.csv', index=False)
    summ = summarize(rows, ('variant', 'form', 'exponent', 'dose', 'reading', 'criterion'))
    summ = summ.merge(e94_costs().drop(columns=['form', 'exponent', 'dose'], errors='ignore'), on='variant', how='left')
    summ.to_csv(out_dir() / 'a3_summary.csv', index=False)
    od = dose_outcome(summ)
    _write_json(out_dir() / 'a3_dose.json', od)
    provenance('A3', variants=sorted(rows.variant.unique()), wirings=list(wirings), MIN_IN_S=MIN_IN_S,
               seconds=round(time.time() - t0, 1))
    return dict(rows=rows, summary=summ, dose=od)


def dose_outcome(summ: pd.DataFrame, min_in: int = MIN_IN_S) -> dict:
    """Outcome (e) per reading (B1's subtraction at exponent 2 and the fits, VERDICT_CRITERION): the doses whose share of
    entries below minus the threshold is inside the flies' range on `min_in` of 4 wirings; 'below_fitted_only' if they
    are all below dose 1 (and some are), 'every_dose' if all are; across readings 'reading_decides' if the sets differ."""
    s = summ[(summ.criterion == VERDICT_CRITERION) & (summ.form == 'sub') & (summ.exponent == 2.0)
             & (summ.variant != 'ref_settled')]
    per = {}
    for rd, sub in s.groupby('reading'):
        doses = sorted(float(x) for x in sub.dose.unique())
        inside = sorted(float(x) for x in sub[sub.inhibited_in >= min_in].dose)
        cls = ('every_dose' if inside == doses else 'below_fitted_only' if inside and max(inside) < 1.0 else
               'none' if not inside else 'other')
        per[rd] = dict(doses=doses, inside=inside, cls=cls,
                       values={f'{d:g}': float(sub[sub.dose == d].inhibited.iloc[0]) for d in doses})
    sets = {tuple(v['inside']) for v in per.values()}
    overall = 'reading_decides' if len(sets) > 1 else next(iter(per.values()))['cls'] if per else 'none'
    return dict(per_reading=per, overall=overall, min_in=min_in, criterion=VERDICT_CRITERION)


# ---------------------------------------------------------------------------
# A4: E114's maps (run again as E114 adds files)

def lab_input() -> pd.DataFrame:
    """The model's ORN input (glomeruli x the 16 panel odors, abbreviated)."""
    lab = H.apl16_orn_deltas().copy()
    lab.columns = [H.odor_abbrev(c[-1] if isinstance(c, tuple) else str(c)) for c in lab.columns]
    return lab


def e114_dir() -> Path:
    return H.shared_dir() / 'results' / 'E114_input_compression'


def part_a4() -> Dict[str, pd.DataFrame]:
    """A4: per E114 input (inputs/<name>.csv; I0 the lab's input) the per-odor totals over the glomeruli (signed and
    positive parts) and their CV over the 12 megamat odors, the largest change of a total from I0's; its PN stage's
    totals' CV (pn_stage.csv) and bouton stage (any '<...>_on_bouton_win.parquet' E114 saved, read with E119's
    statistics on the real wiring) when written; E90's D2 PN settings (S3 exponents, no per-odor rescale) beside.
    Writes a4_inputs.csv, a4_stages.csv."""
    rows = []
    lab = lab_input()
    ins = {'I0': lab}
    for p in sorted((e114_dir() / 'inputs').glob('*.csv')):
        x = pd.read_csv(p, index_col=0)
        x.columns = [H.odor_abbrev(str(c)) for c in x.columns]
        ins[p.stem] = x
    t0 = lab[list(ODORS)].sum(axis=0)
    for name, x in ins.items():
        tot = x[list(ODORS)].sum(axis=0)
        pos = x[list(ODORS)].clip(lower=0).sum(axis=0)
        rows.append(dict(input=name, n_glomeruli=int(x.shape[0]), totals_cv=_cv(tot), pos_totals_cv=_cv(pos),
                         max_rel_change_from_I0=float(((tot - t0).abs() / t0.abs()).max()) if x.shape[0] == lab.shape[0] else np.nan))
    a = pd.DataFrame(rows)
    a.to_csv(out_dir() / 'a4_inputs.csv', index=False)
    st = []
    pp = e114_dir() / 'pn_stage.csv'
    if pp.exists():
        for _, r in pd.read_csv(pp).iterrows():
            st.append(dict(source='E114 pn_stage.csv', input=r['input'], stage='pn', totals_cv=r['cv_totals'],
                           totals_cv_flies=r.get('cv_totals_flies', np.nan)))
    d2 = E90.out_dir() / 'd2_pn_settings.csv'
    if d2.exists():
        for _, r in pd.read_csv(d2).iterrows():
            if r['handle'] == 'S3':
                st.append(dict(source='E90 D2 (S3: positive deltas to a power, no per-odor rescale)', input=r['setting'],
                               stage='pn', totals_cv=r['cv_totals'], totals_cv_flies=r['cv_totals_flies'],
                               cal_pattern=r.get('cal_pattern', np.nan)))
    for root in (e114_dir(), H.outputs_dir('E114_input_compression'), archive('E114_input_compression')):
        for f in sorted(Path(root).rglob('*_on_bouton_win.parquet')) if Path(root).exists() else []:
            pre = f.with_name(f.name.replace('_bouton_win', '_bouton_pre'))
            if not pre.exists():
                continue
            U = indicator_change(pd.read_parquet(f), pd.read_parquet(pre), 1.0, 'evoked')
            units = {'real': stage_subset(U, 'bouton', 'real')}
            r = read_units(units, 'bouton', readings=('rate',))
            r = r[r.criterion == VERDICT_CRITERION].iloc[0]
            st.append(dict(source=str(f), input=f.name.replace('_on_bouton_win.parquet', ''), stage='bouton',
                           totals_cv=r['totals_cv'], resp_cv=r['resp_cv'], top1=r['top1'], inhibited=r['inhibited']))
    s = pd.DataFrame(st)
    s.to_csv(out_dir() / 'a4_stages.csv', index=False)
    provenance('A4', inputs=list(ins), pn_stage=pp.exists(), n_stage_rows=len(s))
    return dict(inputs=a, stages=s)


# ---------------------------------------------------------------------------
# A5: the bound

def part_a5() -> pd.DataFrame:
    """A5: the per-odor totals' CV over the 12 megamat odors from the input to the boutons and claws, model and flies.
    Writes a5_bound.csv (stage, source, value, lo, hi)."""
    rows = []
    lab = lab_input()
    shared = E90.shared_glomeruli()
    rows.append(dict(stage='ORN input', source='model input, 36 glomeruli', value=_cv(lab[list(ODORS)].sum(axis=0))))
    rows.append(dict(stage='ORN input', source='model input, 28 shared glomeruli',
                     value=_cv(lab.loc[shared, list(ODORS)].sum(axis=0))))
    al = pd.read_csv(E90.out_dir() / 'al_transfer.csv')
    sm = al[al.section == 'summary']
    for mp in ('linear', 'saturating'):
        for k, stg in (('cv_total_orn', 'ORNs, megamat recordings (pebbled)'), ('cv_total_pn', 'PN dendrites (GH146)')):
            r = sm[(sm['map'] == mp) & (sm.stat == k)].iloc[0]
            rows.append(dict(stage=stg, source=f'flies, 28 shared glomeruli, {mp}', value=r['value'], lo=r['lo'], hi=r['hi']))
    from experiments import E102_antennal_lobe_normalization as E102
    bp = pd.read_csv(E102.out_dir() / 'b_pn.csv', comment='#')
    for _, r in bp[bp.wiring == 'real'].iterrows():
        rows.append(dict(stage='PN stage', source=f"model, E102 step {r['step']} (28 shared glomeruli)", value=r['cv_totals']))
    a1 = pd.read_csv(E119.out_dir() / 'a1_totals.csv')
    for _, r in a1[(a1.pool.astype(str) == '1') & (a1.reading == 'rate')].iterrows():
        rows.append(dict(stage='boutons', source=f"model default, APL {r['cond']}, rate (8 wirings)", value=r['mean'], lo=r['lo'], hi=r['hi']))
    for stage in ('bouton', 'claw'):
        sp = fly_spreads(stage)[VERDICT_CRITERION]['totals_cv']
        rows.append(dict(stage=f'{stage}s', source=f'flies ({len(fly_stats(stage)[VERDICT_CRITERION])})',
                         value=sp['mean'], lo=sp['lo'], hi=sp['hi']))
    b1 = pd.read_csv(E119.out_dir() / 'b1_roi_stats.csv')
    c = b1[(b1.stage == 'claw') & (b1.reading == 'rate') & (b1.share == 0) & (b1.pool.astype(str) == 'g14')
           & (b1.cond == 'on') & (b1.criterion == VERDICT_CRITERION)]
    if len(c):
        rows.append(dict(stage='claws', source='model default, APL on, rate, g14 (8 wirings)', value=float(c.totals_cv.mean()),
                         lo=float(c.totals_cv.min()), hi=float(c.totals_cv.max())))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'a5_bound.csv', index=False)
    provenance('A5')
    return df


# ---------------------------------------------------------------------------
# part B: the commands, the stage runs at the refits, their reading

def part_b_commands(level: float, lanes: Sequence[Sequence[str]] = (('real', 'P1'), ('P2', 'M1')), threads: int = 16
                    ) -> List[str]:
    """The exact commands for part B on hal (two lanes, from committed code): E102's part D at the level (its start
    checks first), then this script's stage runs at each lane's refits; the scoring after both lanes."""
    env = f'OMP_NUM_THREADS={threads} OMP_THREAD_LIMIT={threads} PYTHONUNBUFFERED=1 nice -n 10'
    py = './venv/bin/python -u'
    out = []
    for lane in lanes:
        w = ','.join(lane)
        out.append(f'{env} {py} experiments/E102_antennal_lobe_normalization.py --part D --level {level:g} --wirings {w} '
                   f'--e120-level && {env} {py} experiments/E120_bouton_stage_measured.py --part Bstage --level {level:g} '
                   f'--wirings {w}')
    out.append(f'{py} experiments/E102_antennal_lobe_normalization.py --part D --level {level:g} --score-refits && '
               f'{py} experiments/E120_bouton_stage_measured.py --part Bread --level {level:g}')
    return out


def stage_dir(level: float) -> Path:
    return local_dir(f'B/{e102_step_tag(level)}')


def refit_stage_run(level: float, wiring: str, force: bool = False) -> dict:
    """One refit of E102's part D at `level` run once more on the lab's 16-odor input, APL on and acutely off (APL->KC
    and APL->bouton 0), with its stage matrices saved (E90's `read_stages` and `_save_stage`, tag '<wiring>_on' /
    '<wiring>_acute_off' under `stage_dir`); its counts and responses must equal E10's saved outputs of the refit
    exactly, otherwise it raises. Returns the record (<wiring>_stage.json)."""
    from experiments import E102_antennal_lobe_normalization as E102
    from experiments import E72_wiring_population as E72
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E98_apl_response_strength as E98
    d = stage_dir(level)
    rec_path = d / f'{wiring}_stage.json'
    rec = _read_json(rec_path)
    if rec is not None and rec.get('exact') and not force:
        return rec
    tag = E102.refit_tag(level, wiring)
    saved = E72.load_fit(E102.EXPERIMENT, tag)
    if saved is None:
        raise RuntimeError(f'{tag}: no saved refit outputs (run E102 part D at {level:g} first)')
    p = E102.refit_params(tag)
    kws = E102.step_kws(level)
    breadth = E102.point_breadth(wiring, kws)
    sim, info = E83.build_simulator(wiring, breadth=breadth, params=dict(p), extra_kws=kws)
    mp, rv, order = E98._sim_cells(sim)
    W = info['W']
    res = dict(wiring=wiring, level=float(level), tag=tag, params=p, pn_kws=kws, commits=H.code_commits(),
               time=time.strftime('%Y-%m-%d %H:%M'))
    exact = True
    for c, ct in (('on', 'on'), ('off', 'acute_off')):
        q = dict(p) if c == 'on' else dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
        r, s = sim(q)
        rs, ss = saved[c]
        nd = int((s.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum())
        nr = int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum())
        res[f'{c}_counts_differ'], res[f'{c}_responses_differ'] = nd, nr
        exact &= nd == 0 and nr == 0 and s.index.equals(ss.index)
        st = E90.read_stages(mp, rv, W, order, [H.odor_abbrev(str(x)) for x in r.columns], E90.MODEL_WINDOW)
        E90._save_stage(d, f'{wiring}_{ct}', st, (r, s))
    res['exact'] = bool(exact)
    _write_json(rec_path, res)
    if not exact:
        raise RuntimeError(f'{wiring} at {level:g}: the stage run differs from the refit (see {rec_path})')
    return res


def part_bstage(level: float, wirings: Sequence[str]) -> None:
    for w in wirings:
        t0 = time.time()
        refit_stage_run(level, w)
        print(f'E120 Bstage {w} at {level:g}: {time.time() - t0:.0f} s', flush=True)


def load_refit_stage(level: float, wiring: str, ct: str) -> Optional[dict]:
    """A part-B stage run (`refit_stage_run`) from the local outputs or the archive."""
    tag = e102_step_tag(level)
    for d in (H.outputs_dir(EXPERIMENT) / 'B' / tag, archive(EXPERIMENT) / 'B' / tag):
        if (d / f'{wiring}_{ct}_bouton_win.parquet').exists():
            return {k: pd.read_parquet(d / f'{wiring}_{ct}_{k}.parquet') for k in
                    ('bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')}
    return None


def part_bread(level: float, wirings: Sequence[str] = S_WIRINGS) -> pd.DataFrame:
    """E119's statistics on part B's stage runs (boutons and claws, APL on and off, the threshold matched on APL on),
    beside the default's on the same wirings. Writes b_stage_rows.csv and b_stage_summary.csv."""
    rows = []
    for stage in ('bouton', 'claw'):
        for c, ct in (('on', 'on'), ('off', 'acute_off')):
            units, units_on = {}, {}
            for w in wirings:
                st, so = load_refit_stage(level, w, ct), load_refit_stage(level, w, 'on')
                if st is None or so is None:
                    continue
                units[w] = stage_subset(indicator_change(st[f'{stage}_win'], st[f'{stage}_pre'], 1.0, 'evoked'), stage, w)
                units_on[w] = stage_subset(indicator_change(so[f'{stage}_win'], so[f'{stage}_pre'], 1.0, 'evoked'), stage, w)
            if units:
                rows.append(read_units(units, stage, match_on=units_on, model=f'refit {e102_step_tag(level)}', cond=c))
            dflt = {w: default_units(stage, w, cond=ct) for w in wirings}
            dflt_on = {w: default_units(stage, w, cond='on') for w in wirings}
            rows.append(read_units(dflt, stage, match_on=dflt_on, model='default', cond=c))
    rows = pd.concat(rows, ignore_index=True)
    rows.to_csv(out_dir() / f'b_stage_rows_{e102_step_tag(level)}.csv', index=False)
    summ = summarize(rows, ('model', 'stage', 'cond', 'reading', 'criterion'), CLAW_STATS)
    summ.to_csv(out_dir() / f'b_stage_summary_{e102_step_tag(level)}.csv', index=False)
    provenance('Bread', level=level, wirings=list(wirings))
    write_variants(level, wirings)
    return summ


# part B ran at olfsysm bde95cb, whose default was multibouton_claws='first'; since the queue's 7af5a27 E3's BOUTON (and
# so E83's simulator for these refits) pins 'first', so a rerun from current code builds the same model
RUN_MULTIBOUTON = 'first'


def stage_record(level: float, wiring: str) -> Optional[dict]:
    """`refit_stage_run`'s record (<wiring>_stage.json) from the local outputs or the archive."""
    tag = e102_step_tag(level)
    for d in (H.outputs_dir(EXPERIMENT) / 'B' / tag, archive(EXPERIMENT) / 'B' / tag):
        if (d / f'{wiring}_stage.json').exists():
            return dict(_read_json(d / f'{wiring}_stage.json'), _dir=str(d))
    return None


def write_variants(level: float, wirings: Sequence[str] = S_WIRINGS) -> Optional[pd.DataFrame]:
    """variants.csv: one row per part-B stage run (tag = its file names' stem under B/<step>/, the wiring), with E102's
    refit parameters as absolute values, the PN normalization, the refit's E102 tag and the overall response rates with
    APL on and acutely off (from the saved responses). The archive's params.csv is built from this table only
    (archive_model_outputs.PARAMS_TABLES); the readings of part A are not model variants."""
    rows = []
    for w in wirings:
        rec = stage_record(level, w)
        if rec is None:
            continue
        d = Path(rec['_dir'])
        r = dict(tag=w, variant=w, part='B', wiring=w, folder=f"B/{e102_step_tag(level)}", **rec.get('pn_kws', {}),
                 **rec['params'], multibouton_claws=RUN_MULTIBOUTON, target_rate=0.068,
                 ref_experiment='E102_antennal_lobe_normalization', ref_fit=rec['tag'], exact=rec.get('exact'))
        for c, ct in (('on', 'on'), ('acute_off', 'acute_off')):
            f = d / f'{w}_{ct}_responses.parquet'
            r[f'rate_{c}'] = float(pd.read_parquet(f).to_numpy(dtype=float).mean()) if f.exists() else np.nan
        rows.append(r)
    if not rows:
        return None
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'variants.csv', index=False)
    return df


def outcome_letter(a1_outcome_a: bool, level: Mapping) -> str:
    """The entry's outcome among (a) to (d): (a) an indicator member passes A1; (b) a step up to 1000 reaches the flies'
    totals; (c) a step above 1000 does; (d) none does. (e), the dose, is read separately (`dose_outcome`)."""
    if a1_outcome_a:
        return 'a'
    if level.get('found_inside'):
        return 'b' if float(level['level']) <= E102_DONE_LEVEL else 'c'
    return 'd'


def part_b_reading(level: float) -> Optional[dict]:
    """Part B's reading for outcome.json: E102's class of the refit at `level` (d_class_<step>.json) and the bouton
    totals' CV (rate reading, verdict criterion) of the default and the refit with APL on, with the number of wirings
    inside the flies' range. None before part B is read."""
    tag = e102_step_tag(level)
    bs = out_dir() / f'b_stage_summary_{tag}.csv'
    cls_path = H.shared_dir() / 'results' / 'E102_antennal_lobe_normalization' / f'd_class_{tag}.json'
    if not bs.exists():
        return None
    b = pd.read_csv(bs)
    b = b[(b.stage == 'bouton') & (b.cond == 'on') & (b.reading == 'rate') & (b.criterion == VERDICT_CRITERION)]
    out = dict(step=tag)
    for m, k in ((f'refit {tag}', 'refit'), ('default', 'default')):
        r = b[b.model == m]
        if len(r):
            out[f'totals_cv_on_{k}'] = float(r.totals_cv.iloc[0])
            out[f'totals_cv_on_{k}_in'] = int(r.totals_cv_in.iloc[0])
            out['n_wirings'] = int(r.n_wirings.iloc[0])
    cls = _read_json(cls_path)
    if cls is not None:
        S = cls.get('S', {})
        out.update(refit_class=S.get('class'), pattern_unchanged=S.get('pattern_unchanged'), pattern_refit=S.get('pattern'),
                   d_pattern=S.get('d_pattern'), d_raw_loss=S.get('d_raw_loss'), n_flagged=S.get('n_flagged'))
    return out


# ---------------------------------------------------------------------------
# the table

def _f(x, spec: str = '.3f') -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return '-'
    return format(x, spec) if np.isfinite(x) else '-'


def _md(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] + ['---:'] * (len(header) - 1)) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return out


STAT_NAMES = {'totals_cv': "per-odor totals' CV", 'resp_cv': "responders' CV", 'top1': 'top-1 % share',
              'inhibited': 'share below -threshold'}


def _range(stage: str, s: str) -> str:
    sp = fly_spreads(stage)[VERDICT_CRITERION][s]
    return f"{_f(sp['mean'])} ({_f(sp['lo'])} to {_f(sp['hi'])})"


def write_table() -> dict:
    """comparison_table.md and outcome.json from part A's tables."""
    od = out_dir()
    L = ['# E120: the bouton stage held to the recorded boutons', '',
         "E119's per-ROI statistics (boutons not pooled, E90's shared-glomerulus subset; claws in groups of 14 of one "
         "glomerulus) at 2 baseline-noise SDs, the model's threshold matched to the flies' mean overall response rate. "
         "Model cells: mean over the wirings [number inside the flies' range]. Flies: mean (range), 5 bouton flies "
         "(2022-10-21 fly02 excluded) and 4 claw flies.", '']
    out = {}
    # A1
    a1 = pd.read_csv(od / 'a1_summary.csv')
    ps = pd.read_csv(od / 'a1_pass.csv')
    L += ['## A1: the indicator family at the default (8 wirings, APL on, boutons)', '']
    hdr = ['member', 'reading'] + [STAT_NAMES[s] for s in STATS] + ['passes']
    rws = [['flies', '-'] + [_range('bouton', s) for s in STATS] + ['']]
    b = a1[(a1.stage == 'bouton') & (a1.criterion == VERDICT_CRITERION)]
    for _, r in b.iterrows():
        p = ps[(ps.n == r.n) & (ps.form == r.form) & (ps.reading == r.reading)].passes.iloc[0]
        rws.append([f"n {r.n:g} {r.form}", r.reading] + [f"{_f(r[s])} [{int(r[f'{s}_in'])}]" for s in STATS]
                   + ['yes' if p else 'no'])
    L += _md(hdr, rws) + ['']
    c = a1[(a1.stage == 'claw') & (a1.criterion == VERDICT_CRITERION) & (a1.reading == 'rate')]
    L += ['Claws (rate reading; per-odor rates and odors per ROI: number of wirings inside the band):', '']
    hdr = ['member'] + [STAT_NAMES[s] for s in STATS] + ['per-odor rates in', 'odors per ROI in']
    rws = [['flies'] + [_range('claw', s) for s in STATS] + ['', '']]
    for _, r in c.iterrows():
        rws.append([f"n {r.n:g} {r.form}"] + [f"{_f(r[s])} [{int(r[f'{s}_in'])}]" for s in STATS]
                   + [int(r['per_odor_n_in']), int(r['opr_n_in'])])
    L += _md(hdr, rws) + ['']
    passing = ps[ps.passes]
    out['a1_passing'] = [f"n {r.n:g} {r.form} {r.reading}" for _, r in passing.iterrows()]
    out['a1_outcome_a'] = bool(len(passing))
    # A2
    a2 = pd.read_csv(od / 'a2_summary.csv')
    costs = pd.read_csv(od / 'a2_costs.csv')
    idn = pd.read_csv(od / 'a2_identity.csv')
    L += ["## A2: the normalization ladder (E102's saved PN stage; boutons with APL off, 8 wirings)", '',
          f"Identity: at E102's identity step the PN-derived boutons equal E119's APL-off boutons within "
          f"{_f(idn.max_abs_diff.max(), '.4f')} Hz (E90's pre-odor window against the settled level; largest change "
          f"{_f(idn.max_abs_ref.max(), '.1f')} Hz). E102's cost: pattern part under the calibrated readout against all "
          "recordings, its change from the default, population of the eight wirings at fixed scales; raw-count loss change.", '']
    hdr = ['step', 'reading'] + [STAT_NAMES[s] for s in STATS] + ["totals' CV, APL-on estimate", 'pattern part (cal)',
                                                                    'change', 'mean APL effect (cal)', 'raw-loss change']
    rws = [['flies', '-'] + [_range('bouton', s) for s in STATS] + [_range('bouton', 'totals_cv'), '', '', _f(H.data_target('apl_effect_mean')), '']]
    s2 = a2[a2.criterion == VERDICT_CRITERION]
    for _, r in s2.iterrows():
        cs = costs[(costs.step == r.step) & (costs.unit == 'all8')]
        cc = cs.iloc[0] if len(cs) else {}
        rws.append([r.step, r.reading] + [f"{_f(r[s])} [{int(r[f'{s}_in'])}]" for s in STATS]
                   + [f"{_f(r['totals_cv_on_est'])} [{int(r['totals_cv_on_est_in'])}]",
                      _f(cc.get('pattern_cal')), _f(cc.get('d_pattern_cal'), '+.4f'), _f(cc.get('mean_effect_cal')),
                      _f(cc.get('d_loss_raw'), '+.4f')])
    L += _md(hdr, rws) + ['']
    lv = choose_level(a2)
    _write_json(od / 'level.json', dict(**lv, time=time.strftime('%Y-%m-%d %H:%M'), commits=H.code_commits()))
    out['level'] = lv
    L += [f"Level for part B: pn_inhadd {lv['level']:g} ({'inside on 6 of 8' if lv['found_inside'] else 'no step inside; the fallback'}); "
          f"part B {'runs' if lv['part_b_runs'] else 'does not run'}. Beside (first step inside, by reading): "
          + ', '.join(f"{k} {v if v is None else format(v, 'g')}" for k, v in lv['beside_first_inside'].items()) + '.', '']
    # A3
    a3 = pd.read_csv(od / 'a3_summary.csv')
    dose = _read_json(od / 'a3_dose.json')
    L += ["## A3: the bouton-dose ladder (E94's saved settings, 4 searched wirings, APL on)", '',
          "Inside: of 4 wirings. 'ref' is the fits (dose 1, E90's archive) with the pre-odor mean as saved, as E94's runs; "
          "'ref_settled' with the settled level (E119's reading). E94's cost: the pattern part's change under the "
          "calibrated readout against all recordings (S), the raw-count loss change, the mean APL effect's change.", '']
    hdr = ['setting', 'reading'] + [STAT_NAMES[s] for s in STATS] + ['pattern change', 'raw-loss change', 'mean-effect change']
    rws = [['flies', '-'] + [_range('bouton', s) for s in STATS] + ['', '', '']]
    s3 = a3[a3.criterion == VERDICT_CRITERION].copy()
    s3['_o'] = s3.form.map({'sub': 0, 'gain': 1}) * 100 + s3.exponent * 10 + s3.dose
    for _, r in s3.sort_values(['_o', 'variant', 'reading']).iterrows():
        rws.append([r.variant, r.reading] + [f"{_f(r[s])} [{int(r[f'{s}_in'])}]" for s in STATS]
                   + [_f(r.get('d_pattern_z8'), '+.4f'), _f(r.get('d_raw_loss'), '+.4f'), _f(r.get('d_mean_effect'), '+.3f')])
    L += _md(hdr, rws) + ['']
    L += [f"Outcome (e): {dose['overall']} (" + '; '.join(f"{k}: inside at doses {v['inside']}" for k, v in
                                                       dose['per_reading'].items()) + ').', '']
    out['dose'] = dose['overall']
    # A4
    L += ["## A4: E114's maps", '']
    a4i = pd.read_csv(od / 'a4_inputs.csv') if (od / 'a4_inputs.csv').exists() else pd.DataFrame()
    a4s = pd.read_csv(od / 'a4_stages.csv') if (od / 'a4_stages.csv').exists() and (od / 'a4_stages.csv').stat().st_size > 1 else pd.DataFrame()
    if len(a4i):
        L += _md(['input', "totals' CV (signed)", "positive totals' CV", 'largest change of a total from I0'],
                 [[r.input, _f(r.totals_cv), _f(r.pos_totals_cv), _f(r.max_rel_change_from_I0, '.2e')] for _, r in a4i.iterrows()]) + ['']
    if len(a4s):
        L += _md(['source', 'input', 'stage', "totals' CV"], [[r.source, r.input, r.stage, _f(r.totals_cv)] for _, r in a4s.iterrows()]) + ['']
    has_e114_stage = bool(len(a4s) and a4s.source.astype(str).str.contains('E114').any())
    L += [("E114's PN or bouton stage is read above." if has_e114_stage else
           "E114 has not written its PN stage (pn_stage.csv) or saved bouton stage files yet; run `--part A4` again "
           "when it has."), '']
    out['a4_e114_stage'] = has_e114_stage
    # A5
    a5 = pd.read_csv(od / 'a5_bound.csv') if (od / 'a5_bound.csv').exists() else pd.DataFrame(columns=['stage', 'source', 'value'])
    L += ["## A5: the bound (per-odor totals' CV over the 12 megamat odors)", '']
    L += _md(['stage', 'source', 'value', '90 % interval or range'],
             [[r.stage, r.source, _f(r.value), f"{_f(r.get('lo'))} to {_f(r.get('hi'))}" if np.isfinite(r.get('lo', np.nan)) else '-']
              for _, r in a5.iterrows()]) + ['']
    out['part_b_runs'] = bool(lv['part_b_runs'] and not out['a1_outcome_a'])
    out['part_b_reason'] = ('A1: an indicator member brings the boutons inside' if out['a1_outcome_a'] else
                            'the level is 1000, E102 part D is the answer' if not lv['part_b_runs'] else
                            f"level {lv['level']:g}")
    out['part_b_commands'] = part_b_commands(lv['level']) if out['part_b_runs'] else []
    L += ['## Part B', '', f"Runs: {'yes' if out['part_b_runs'] else 'no'} ({out['part_b_reason']}).", '']
    if out['part_b_runs']:
        L += ['```'] + out['part_b_commands'] + ['```', '']
    pb = part_b_reading(lv['level']) if out['part_b_runs'] else None
    out['part_b'] = pb
    if pb is not None:
        tag = pb['step']
        bs = pd.read_csv(od / f'b_stage_summary_{tag}.csv')
        bs = bs[(bs.criterion == VERDICT_CRITERION)]
        L += [f"E119's statistics on the refits at {tag} ({pb.get('n_wirings', '-')} searched wirings) and on the default "
              "fits of the same wirings; threshold matched on APL on.", '']
        hdr = ['model', 'stage', 'APL', 'reading'] + [STAT_NAMES[s] for s in STATS]
        rws = []
        for stage in ('bouton', 'claw'):
            rws.append(['flies', stage, 'on', '-'] + [_range(stage, s) for s in STATS])
            for _, r in bs[bs.stage == stage].iterrows():
                rws.append([r.model, stage, r.cond, r.reading] + [f"{_f(r[s])} [{int(r[f'{s}_in'])}]" for s in STATS])
        L += _md(hdr, rws) + ['']
        if pb.get('refit_class') is not None:
            L += [f"E102's class of the refit at {tag} (pattern part under the calibrated readout, the four searched "
                  f"wirings): {pb['refit_class']} ({_f(pb['pattern_unchanged'], '.4f')} unchanged, "
                  f"{_f(pb['pattern_refit'], '.4f')} refitted; {pb.get('n_flagged')} flagged).", '']
    out['letter'] = outcome_letter(out['a1_outcome_a'], lv)
    L += [f"Outcome: ({out['letter']}); the dose, (e): {out['dose']}.", '']
    (od / 'comparison_table.md').write_text('\n'.join(L) + '\n')
    _write_json(od / 'outcome.json', dict(**out, time=time.strftime('%Y-%m-%d %H:%M'), commits=H.code_commits()))
    return out


def part_a() -> dict:
    part_a1()
    part_a2()
    part_a3()
    part_a4()
    part_a5()
    return write_table()


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('A', 'A1', 'A2', 'A3', 'A4', 'A5', 'table', 'Bstage', 'Bread',
                                                      'commands'))
    ap.add_argument('--level', type=float, default=None)
    ap.add_argument('--wirings', default=','.join(S_WIRINGS))
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args(argv)
    ws = [w.strip() for w in a.wirings.split(',') if w.strip()]
    if a.part == 'A':
        print(json.dumps(part_a(), indent=1, default=str))
    elif a.part in ('A1', 'A2', 'A3', 'A4', 'A5'):
        {'A1': part_a1, 'A2': part_a2, 'A3': part_a3, 'A4': part_a4, 'A5': part_a5}[a.part]()
        if (out_dir() / 'a1_summary.csv').exists() and (out_dir() / 'a2_summary.csv').exists() and (out_dir() / 'a3_summary.csv').exists():
            write_table()
    elif a.part == 'table':
        print(json.dumps(write_table(), indent=1, default=str))
    elif a.part == 'commands':
        lv = _read_json(out_dir() / 'level.json') or {}
        print('\n'.join(part_b_commands(a.level if a.level is not None else lv['level'])))
    else:
        if a.level is None:
            raise SystemExit(f'--part {a.part} needs --level')
        if a.part == 'Bstage':
            part_bstage(a.level, ws)
        else:
            print(part_bread(a.level, ws).to_string())


if __name__ == '__main__':
    main()
