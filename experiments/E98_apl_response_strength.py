#!/usr/bin/env python3
"""E98: the strength of APL's response in the flies against the model's APL. In the best family APL's activity is
nearly the same for every odor of the panel, and two fifths of it is there before the odor; whether Yang's 2021 calyx
maps, read as each fly's APL response by odor and concentration, agree with that APL wherever the two can be compared,
and which recordings would tell a flat APL from one that follows the KCs' response. A diagnostic: it moves no loss.

The entry is item E98 of experiment_priorities.md (proposed by atlas 2026-11-15 12:45; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-15 12:26). C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the
0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's `saved_fit` / `build_simulator`), nothing
re-tuned. The flies: Yang's 2021 single-plane APL calyx maps (data/yang_apl_calyx_maps/2021_APL_18um/ in the shared
directory; 12 recordings, 8 flies; GCaMP6f in APL; pa, paa and ms as concentration series; 7 recordings with a 1 s
pulse, 5 with 0.5 s), read with E79's reader, baselines, blank and fly bootstrap (experiments/E79_apl_calyx_constraints.py).

Definitions (the entry's, fixed before any reading):
- A recording's response to a stimulus: the calyx-mean peak dF/F inside the outline under E79's order-matched baseline,
  the mean of the two presentations, with the recording's position-matched blank subtracted (MAIN_READING). A recording
  enters a ratio only if both responses exceed MIN_RESPONSE (0.02).
- A concentration that was not imaged: linear interpolation of the response on log10 concentration between its two
  imaged neighbours, inside the imaged range only; outside it, the nearest imaged concentration stands in when it is
  within MATCH_LOG_UNITS (1, inclusive; Tom's rule) of the one asked for (`value_at`).
- The flies' ratio for two stimuli: the geometric mean over the 8 flies (a fly's recordings averaged on the log scale),
  with a 90 % interval from N_BOOT (10,000) bootstrap draws of flies (E79's `fly_bootstrap`, seed SEED).
- The model's APL for a stimulus: the mean of APL's activity in the odor window (a; E90's MODEL_WINDOW, 0-0.5 s) and its
  pre-odor mean (a0), as E90's stage runs save them (`<wiring>_on_apl.csv`).
- The model's ratio through an indicator of exponent n (INDICATOR_NS, 1 to 2.5 in steps of 0.25), in two forms
  (`indicator_ratio`): 'baseline', (a1^n - a0^n) / (a2^n - a0^n), the pre-odor level in the indicator's baseline; and
  'evoked', ((a1 - a0) / (a2 - a0))^n, only the evoked part seen. GCaMP6f's Hill coefficient is about 2.3 (a literature
  value: Chen et al. 2013, Nature 499:295).
- A wiring is inside a comparison if any of its ratios over n and the two forms lies in the flies' interval; flatter if
  all of them are nearer 1 than the interval or on the other side of 1; steeper if all are beyond the interval's far end
  (`classify`). A comparison is inside, flatter or steeper if at least VERDICT_MIN_WIRINGS (6) of the 8 wirings are;
  otherwise not decided (`comparison_verdict`). The linear reading (n 1, evoked) is given beside every verdict.

Parts (CLI `--part`):
- 0, the identities (identity.json; part A refuses to run until the reader's pass, part B until all pass): the reader
  reproduces E79's per-stimulus calyx means under its three baselines (E79_maps_across_baselines.csv) and its blank
  table (E79_baseline_noise.csv, E79_blank_vs_odor_maps.csv); the extended-input run reproduces each wiring's saved run
  on the 16 panel odors with 0 count differences and E90's saved APL window and pre-odor means (real wiring: mean over
  odors 19.548, pre-odor share 0.423); a blank stimulus gives an evoked APL of zero (`--no-runs`: the reader only).
- A (no simulations): A1 the flies' responses per recording and stimulus under E79's three baselines, with and without
  the blank subtracted, the two presentations' agreement, the block order and the drift of the first-in-block paraffin
  oil response (a1_responses.csv, a1_agreement.csv, a1_pfo_drift.csv); A2 the dose-response per odor, recording and fly:
  the fold change per log unit, the lowest concentration above the blank, the values at the model inputs' concentrations
  (a2_dose_response.csv, a2_at_model_concs.csv, a2_medians.csv); A3 the fixed list of ratios (RATIOS: V1 to V5 and the
  entry's probe rows) with their intervals, under every reading and the sensitivity subsets (SUBSETS: the 1 s and the
  0.5 s recordings apart; without the two recordings with a large blank and 20210805_3FRfun1) (a3_ratios.csv,
  a3_per_fly.csv); A4 the family's saved runs (E90's stage runs): APL's window mean, pre-odor mean and evoked activity
  per odor on the eight wirings, the profile over the 16 odors against each odor's input total and PN total, and the
  model's ratios of the fixed list that the panel's odors give, through the indicator family (a4_apl.csv,
  a4_profile.csv, a4_model_ratios.csv, a4_preview.csv); A5 what other APLs would show for the same ratios: one in
  proportion to the KCs' response rates in each genotype, one in proportion to the PN dendrites' totals (the megamat
  odors, E90's tables), E37's saved table for the spike-driven fits (a5_alternatives.csv); A6 the power table: for every
  pair of panel odors the ratio by the family's APL and by a proportional APL, and the number of flies that separates the
  two (a6_power.csv, the five that need the fewest in a_summary.json). a_summary.json.
- B (simulations at fixed parameters, APL on and acutely off, each saved with its counts, responses and APL means; phase
  2): B1 the extended input on the eight wirings (the real wiring first), restricted to the input's 36 glomeruli, 40
  stimuli (B1_STIMULI): the 16 panel odors; paa at 10^-5 from the validation2 session's diagnostics, with the megamat and
  natural-mixture sessions' beside; the other sessions' ms at 10^-3; 2h, t2h, 2-but and 1-6ol at 10^-6 and 2h at 10^-5;
  ma at 10^-7 and ea at 10^-8 (E49's stimuli); Hallem's pa and ms at 10^-2, -4, -6 and -8 and paa at 10^-2 (drosolf's
  table, E38's conversion); a blank. With B1 the KCs' response rate to the blank, APL on and off, against the flies'
  solvent. B2 (real wiring): the 1 s pulse; the fuller input of E93; E95's refit at the connectome's bouton share. B3
  (real wiring): E36's KC-spike-driven fit and E37's spike-only fit at 40 Hz with the gain. B4 (real wiring): E77's
  screen settings with weaker PN normalization (pn_inhadd 1000 and 3141) at their held readings.
- C (no simulations; phase 2): the verdict, fixed here: V1 (paa over ms near 10^-3, the lab's input, the hull of the two
  readings of paa's concentration, the model under each estimate of ms's input; `v1_verdict`) and V2 (paa at 10^-5 over
  paa at 10^-4); V3 to V5 (Hallem's input), B2 to B4 and the pulse-matched readings described beside (verdict.json).
- table: comparison_table.md and the figures.

Deviations from the entry (each also in the note):
- MAIN_READING subtracts the position-matched blank from the order-matched responses, as the entry defines. That blank is
  the paraffin oil block's 2nd and 3rd presentations minus the first-in-block ones, and the order-matched baseline
  already subtracts the paraffin oil block's presentation at the same position, so the first-in-block baseline with the
  blank subtracted equals the order-matched baseline without it exactly (`blank_identity`): the main reading removes the
  position effect twice. Every ratio is therefore also given under the other five readings (a3_ratios.csv), the
  order-matched one without the blank (the entry's probes) first among them.
- The one-log-unit rule (`value_at`) stands in a recording's nearest imaged concentration only for a concentration
  outside the series' design range (the union of the recordings' ranges); one inside it that a recording did not image
  (ms at 10^-2 in 20210805_3FRfun1 and 20210823_4FRfun2) leaves that recording out of readings that need it, rather than
  reading those two recordings at another concentration than the other ten. A3's subset 'all_recording_ranges' gives the
  other reading (the rule on each recording's own range; ms at 10^-3 is then read at 10^-3.431 in those two).
- E79's reader leaves out the doubled pa 10^-6.771 block of 20210805_3FRfun1 and 20210823_4FRfun2 (Yang, 2026-11-13: both
  are pa 10^-6.771 repeats); kept out here so that the reader reproduces E79, and no ratio of the fixed list reads it.
- E37's saved APL table has no pre-odor level and is the mean (and peak) over the whole time after odor onset, not over
  0-0.5 s, so its ratios in A5 are of the window values, not of the evoked activity.

Outputs: results/E98_apl_response_strength/ (the tables above, identity.json, provenance_<part>.json, comparison_table.md,
figures); counts, responses and APL means of part B under outputs/E98_apl_response_strength/, archived to
data/model_outputs/E98_apl_response_strength/.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part 0 --no-runs
    ./venv/bin/python -u -m experiments.E98_apl_response_strength --part A
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E98_apl_response_strength --part B --arm B1 --wirings real --threads 4
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E98_apl_response_strength'
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's ALL_WIRINGS (S then T)

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any reading)

MIN_RESPONSE = 0.02                    # a recording enters a ratio only if both responses exceed this (dF/F)
MATCH_LOG_UNITS = 1.0                  # Tom's rule: a concentration within one log10 unit (inclusive) stands in
N_BOOT = 10000                         # bootstrap draws of flies
SEED = 0
CI = (5.0, 95.0)                       # 90 % intervals
INDICATOR_NS = tuple(float(x) for x in np.arange(1.0, 2.5 + 1e-9, 0.25))
FORMS = ('baseline', 'evoked')
VERDICT_MIN_WIRINGS = 6                # of the 8 wirings
BASELINES = ('order_matched', 'first_in_block', 'all_pfo')          # E79's three
MAIN_READING = ('order_matched', True)                             # (baseline, blank subtracted)
READINGS = (MAIN_READING, ('order_matched', False), ('first_in_block', True), ('first_in_block', False),
            ('all_pfo', True), ('all_pfo', False))
LARGE_BLANK = 0.05                     # a recording's blank calyx mean above this is large (E79: 0.09 and 0.18, the rest ~0)
NOT_IN_18UM_LIST = ('20210805_3FRfun1',)
SUBSETS = ('all', 'pulse_1s', 'pulse_0.5s', 'without_large_blank_and_3FRfun1', 'all_recording_ranges')
POWER_MIN_CONC = -5.0                  # A6's between-fly SD: A3's across-odor ratios with both stimuli at this or above
POWER_Z = (1.6449, 0.8416)             # A6: a two-sided 90 % interval excluding the other value with 80 % power

# the fixed list of ratios: (id, numerator, denominator, group); a stimulus is (odor, log10 concentration)
Stim = Tuple[str, float]
RATIOS: Tuple[Tuple[str, Stim, Stim, str], ...] = (
    ('V1_paa-4', ('paa', -4.0), ('ms', -3.0), 'V1'),
    ('V1_paa-2.954', ('paa', -2.954), ('ms', -3.0), 'V1'),
    ('V2', ('paa', -5.0), ('paa', -4.0), 'V2'),
    ('V3_ms-4', ('ms', -4.0), ('ms', -2.0), 'V3'),
    ('V3_ms-6', ('ms', -6.0), ('ms', -2.0), 'V3'),
    ('V3_ms-8', ('ms', -8.0), ('ms', -2.0), 'V3'),
    ('V4', ('pa', -6.0), ('pa', -4.0), 'V4'),
    ('V5_-4', ('pa', -4.0), ('ms', -4.0), 'V5'),
    ('V5_-6', ('pa', -6.0), ('ms', -6.0), 'V5'),
    ('V5_paa-2', ('paa', -2.0), ('ms', -2.0), 'V5'),
    ('probe_pa-4.386', ('pa', -4.386), ('ms', -3.0), 'probe'),
    ('probe_-4.863', ('pa', -4.863), ('ms', -4.863), 'probe'),
)
# the planner's probes (the entry's table; E79's order-matched means, no blank subtracted)
PROBES: Dict[str, Tuple[float, float, float]] = {
    'V1_paa-4': (0.85, 0.70, 0.97), 'V1_paa-2.954': (1.14, 1.03, 1.25), 'probe_pa-4.386': (1.14, 1.05, 1.24),
    'V2': (0.53, 0.46, 0.61), 'probe_-4.863': (1.96, 1.61, 2.39)}
# the model inputs' concentrations (log10) of the three odors: the lab's input (megamat pa and ms at -3, validation2 paa
# at -4), the KC panel's (paa at -2.954), B1's added stimuli and Hallem's series
MODEL_CONCS: Dict[str, Tuple[float, ...]] = {'pa': (-3.0, -4.0, -6.0, -8.0, -2.0), 'ms': (-3.0, -2.0, -4.0, -6.0, -8.0),
                                             'paa': (-4.0, -2.954, -5.0, -2.0)}
# the panel's inputs of the three odors (odor -> log10 concentration of the lab's input)
PANEL_INPUT_CONC = {'pa': -3.0, 'ms': -3.0, 'paa': -4.0}


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(path.read_text()) if path.exists() else None


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                                         **extra))
    _write_json(path, d)


# ---------------------------------------------------------------------------
# the flies: responses, concentrations, ratios

CONC_TOL = 1e-6                        # two log10 concentrations closer than this are the same


def reading_label(reading: Tuple[str, bool]) -> str:
    return f"{reading[0]}{'_blank' if reading[1] else ''}"


def recording_info() -> pd.DataFrame:
    """recordings.csv indexed by recording name, with 'fly' (E79's `<date>_fly<n>`) and 'pulse_s'."""
    from experiments import E79_apl_calyx_constraints as E79
    r = pd.read_csv(E79.DATA / 'recordings.csv')
    r.index = [f.split('_moco')[0] for f in r['file']]
    r['fly'] = [f'{d}_fly{f}' for d, f in zip(r['date'], r['fly'])]
    r['pulse_s'] = r['odor_pulse_s'].astype(float)
    return r


def load_recordings():
    """E79's `load_recordings` (the 12 recordings, their stacks, frames and calyx masks)."""
    from experiments import E79_apl_calyx_constraints as E79
    return E79.load_recordings()


def stimulus_responses(recs, pulses: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """Per recording, stimulus, baseline (BASELINES) and blank (subtracted or not): the calyx means of the two
    presentations ('rep1', 'rep2', in frame order) and of their mean ('response'; E79's `stimulus_maps`, doubtful
    presentations left out), the recording's blank calyx mean ('blank'), its fly, side, pulse (`pulses`, default
    recordings.csv's) and the stimulus's block."""
    from experiments import E79_apl_calyx_constraints as E79
    if pulses is None:
        pulses = recording_info()['pulse_s'].to_dict()
    rows = []
    for rec in recs:
        fr, st, mask = rec.frames, rec.stack, rec.mask
        blank = float(E79.position_matched_blank(st, fr)[mask].mean())
        ok = fr[(fr['odor'] != 'pfo') & ~fr['doubtful']]
        block = {(o, float(c)): int(g['block'].min()) for (o, c), g in ok.groupby(['odor', 'log10_conc'])}
        for bl in BASELINES:
            for (o, c), (m, reps) in E79.stimulus_maps(st, fr, bl).items():
                if len(reps) != 2:
                    raise RuntimeError(f'{rec.name} {o} {c:g}: {len(reps)} presentations')
                vals = dict(response=float(m[mask].mean()), rep1=float(reps[0][mask].mean()), rep2=float(reps[1][mask].mean()))
                for sub in (False, True):
                    rows.append(dict(recording=rec.name, fly=rec.fly, side=rec.side, pulse_s=float(pulses[rec.name]),
                                     odor=o, conc=float(c), stimulus=E79.stim_label((o, c)), block=block[(o, float(c))],
                                     baseline=bl, blank_subtracted=sub, blank=blank,
                                     **{k: v - blank if sub else v for k, v in vals.items()}))
    return pd.DataFrame(rows)


def _reading(resp: pd.DataFrame, reading: Tuple[str, bool]) -> pd.DataFrame:
    return resp[(resp['baseline'] == reading[0]) & (resp['blank_subtracted'].astype(bool) == bool(reading[1]))]


def blank_identity(resp: pd.DataFrame) -> float:
    """The largest |first-in-block with the blank subtracted - order-matched without it| over recordings and stimuli
    (zero up to rounding: the two are the same reading)."""
    a = _reading(resp, ('first_in_block', True)).set_index(['recording', 'stimulus'])['response']
    b = _reading(resp, ('order_matched', False)).set_index(['recording', 'stimulus'])['response']
    if not a.index.sort_values().equals(b.index.sort_values()):
        raise RuntimeError('the two readings have different stimuli')
    return float((a - b.reindex(a.index)).abs().max())


def value_at(concs: Sequence[float], values: Sequence[float], c: float,
             design_range: Optional[Tuple[float, float]] = None) -> Tuple[float, str]:
    """A recording's response at log10 concentration `c` from its imaged concentrations `concs` and responses `values`:
    ('imaged') if imaged; ('interpolated') linear on log10 concentration between its two imaged neighbours inside its
    imaged range; ('substituted') the nearest imaged end when `c` is outside the series' `design_range` (default the
    recording's own range) and within MATCH_LOG_UNITS of that end, inclusive; otherwise NaN ('outside', or 'missing'
    when `c` is inside the design range but outside the recording's own)."""
    x = np.asarray(concs, dtype=float)
    y = np.asarray(values, dtype=float)
    if x.size == 0:
        return np.nan, 'missing'
    o = np.argsort(x)
    x, y = x[o], y[o]
    if np.any(np.diff(x) < CONC_TOL):
        raise ValueError(f'repeated concentrations {list(x)}')
    c = float(c)
    hit = np.flatnonzero(np.abs(x - c) <= CONC_TOL)
    if hit.size:
        return float(y[hit[0]]), 'imaged'
    lo, hi = float(x[0]), float(x[-1])
    if lo < c < hi:
        return float(np.interp(c, x, y)), 'interpolated'
    dlo, dhi = design_range if design_range is not None else (lo, hi)
    if dlo - CONC_TOL <= c <= dhi + CONC_TOL:
        return np.nan, 'missing'
    m = concentration_match(c, [lo, hi])
    if m is None:
        return np.nan, 'outside'
    return float(y[0] if m == lo else y[-1]), 'substituted'


def concentration_match(c: float, imaged: Sequence[float], max_units: float = MATCH_LOG_UNITS) -> Optional[float]:
    """The imaged concentration nearest to `c` if it is within `max_units` log10 units (inclusive), else None."""
    x = np.asarray(imaged, dtype=float)
    if x.size == 0:
        return None
    d = np.abs(x - float(c))
    i = int(np.argmin(d))
    return float(x[i]) if d[i] <= max_units + CONC_TOL else None


def fly_ratio(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The flies' ratio from per-recording log ratios (columns 'fly', 'log_ratio'; natural log): per fly the mean on
    the log scale, the geometric mean over flies, its 90 % interval from `n_boot` bootstrap draws of flies (E79's
    `fly_bootstrap`), the between-fly SD of the log2 ratio and the counts."""
    from experiments import E79_apl_calyx_constraints as E79
    d = per_recording.dropna(subset=['log_ratio'])
    per_fly = d.groupby('fly')['log_ratio'].mean()
    out = dict(n_flies=int(len(per_fly)), n_recordings=int(len(d)), ratio=np.nan, lo90=np.nan, hi90=np.nan,
               sd_log2_flies=np.nan, n_flies_above_1=int((per_fly > 0).sum()))
    if len(per_fly) == 0:
        return out
    mean, lo, hi = E79.fly_bootstrap(per_fly.to_numpy(), n_boot=n_boot, seed=seed)
    out.update(ratio=float(np.exp(mean[0])), lo90=float(np.exp(lo[0])), hi90=float(np.exp(hi[0])))
    if len(per_fly) > 1:
        out['sd_log2_flies'] = float(per_fly.std(ddof=1) / np.log(2))
    return out


def design_ranges(resp: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Each odor's concentration range over all recordings (the series' design range)."""
    g = resp.groupby('odor')['conc']
    return {o: (float(g.min()[o]), float(g.max()[o])) for o in g.groups}


def subset_recordings(resp: pd.DataFrame, subset: str) -> List[str]:
    """The recordings of a subset (SUBSETS): 'all'; 'pulse_1s' / 'pulse_0.5s'; 'without_large_blank_and_3FRfun1' (the
    recordings whose blank's calyx mean exceeds LARGE_BLANK, two, and NOT_IN_18UM_LIST left out); 'all_recording_ranges'
    (all, with the one-log rule applied to each recording's own range: `flies_ratios`)."""
    r = resp.drop_duplicates('recording').set_index('recording')
    if subset in ('all', 'all_recording_ranges'):
        keep = r.index
    elif subset == 'pulse_1s':
        keep = r.index[np.isclose(r['pulse_s'], 1.0)]
    elif subset == 'pulse_0.5s':
        keep = r.index[np.isclose(r['pulse_s'], 0.5)]
    elif subset == 'without_large_blank_and_3FRfun1':
        large = r.index[r['blank'] > LARGE_BLANK]
        if len(large) != 2:
            raise RuntimeError(f'expected two recordings with a large blank, found {list(large)}')
        keep = [x for x in r.index if x not in set(large) | set(NOT_IN_18UM_LIST)]
    else:
        raise ValueError(f'unknown subset {subset!r}')
    return sorted(keep)


def _stim(s: Stim) -> str:
    return f'{s[0]} {s[1]:g}'


def flies_ratios(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING, subset: str = 'all',
                 ratios=RATIOS) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: every ratio of `ratios` under one reading and subset of the recordings: (one row per ratio, per recording
    rows with each response, how it was read and whether the recording entered). The one-log rule stands in a
    recording's end for a concentration outside the series' design range (`value_at`), or, in subset
    'all_recording_ranges', outside the recording's own range (ms at 10^-3 then reads ms at 10^-3.431 in the two
    recordings without ms at 10^-2)."""
    dr = design_ranges(resp) if subset != 'all_recording_ranges' else {}
    keep = set(subset_recordings(resp, subset))
    r = _reading(resp, reading)
    r = r[r['recording'].isin(keep)]
    rec_rows, rows = [], []
    for rid, num, den, group in ratios:
        per = []
        for rec, g in r.groupby('recording'):
            row = dict(ratio_id=rid, reading=reading_label(reading), subset=subset, recording=rec, fly=g['fly'].iloc[0],
                       pulse_s=float(g['pulse_s'].iloc[0]))
            for side, (o, c) in (('num', num), ('den', den)):
                go = g[g['odor'] == o]
                v, how = value_at(go['conc'], go['response'], c, dr.get(o))
                row[f'{side}_value'], row[f'{side}_how'] = v, how
            row['entered'] = bool(np.isfinite(row['num_value']) and np.isfinite(row['den_value'])
                                  and row['num_value'] > MIN_RESPONSE and row['den_value'] > MIN_RESPONSE)
            row['log_ratio'] = float(np.log(row['num_value'] / row['den_value'])) if row['entered'] else np.nan
            per.append(row)
        pr = pd.DataFrame(per)
        rec_rows.append(pr)
        fr = fly_ratio(pr[pr['entered']])
        probe = PROBES.get(rid, (np.nan,) * 3)
        rows.append(dict(ratio_id=rid, group=group, numerator=_stim(num), denominator=_stim(den),
                         reading=reading_label(reading), subset=subset, **fr,
                         n_recordings_in_subset=int(len(pr)), n_recordings_left_out=int((~pr['entered']).sum()),
                         num_how=';'.join(sorted(set(pr['num_how']))), den_how=';'.join(sorted(set(pr['den_how']))),
                         probe_ratio=probe[0], probe_lo90=probe[1], probe_hi90=probe[2]))
    return pd.DataFrame(rows), pd.concat(rec_rows, ignore_index=True)


def presentation_se(g: pd.DataFrame) -> float:
    """A recording's typical standard error of a stimulus's response, the median over its stimuli of |rep1 - rep2| / 2."""
    return float(np.median(np.abs(g['rep1'] - g['rep2']) / 2.0))


def dose_response(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A2: per recording and odor the fold change of the response per log10 unit (exp of the slope of ln response on
    log10 concentration over the concentrations above MIN_RESPONSE, two or more needed), the linear slope (dF/F per log10
    unit, all concentrations), the lowest concentration above the blank (the lowest whose response, the blank
    subtracted, exceeds both MIN_RESPONSE and twice the recording's `presentation_se`); and the values at MODEL_CONCS
    (`value_at`)."""
    from scipy.stats import spearmanr
    dr = design_ranges(resp)
    r = _reading(resp, reading)
    rows, at_rows = [], []
    for rec, g in r.groupby('recording'):
        se = presentation_se(g)
        for o, go in g.groupby('odor'):
            go = go.sort_values('conc')
            x, y = go['conc'].to_numpy(float), go['response'].to_numpy(float)
            above = y > MIN_RESPONSE
            fold = float(np.exp(np.polyfit(x[above], np.log(y[above]), 1)[0])) if above.sum() >= 2 else np.nan
            thr = max(MIN_RESPONSE, 2.0 * se)
            over = np.flatnonzero(y > thr)
            rows.append(dict(recording=rec, fly=go['fly'].iloc[0], pulse_s=float(go['pulse_s'].iloc[0]), odor=o,
                             reading=reading_label(reading), n_concs=int(len(x)), n_above=int(above.sum()),
                             fold_per_log_unit=fold, slope_dff_per_log_unit=float(np.polyfit(x, y, 1)[0]),
                             spearman_conc=float(spearmanr(x, y)[0]), presentation_se=se, threshold=thr,
                             lowest_above_blank=float(x[over[0]]) if over.size else np.nan,
                             all_above_from_lowest=bool(over.size and (y[over[0]:] > thr).all()),
                             response_min=float(y.min()), response_max=float(y.max())))
            for c in MODEL_CONCS.get(o, ()):
                v, how = value_at(x, y, c, dr[o])
                at_rows.append(dict(recording=rec, fly=go['fly'].iloc[0], pulse_s=float(go['pulse_s'].iloc[0]), odor=o,
                                    conc=c, reading=reading_label(reading), value=v, how=how))
    return pd.DataFrame(rows), pd.DataFrame(at_rows)


def fly_geomean(per_recording: pd.DataFrame, col: str) -> dict:
    """`fly_ratio` of a positive per-recording quantity (its log averaged within flies)."""
    d = per_recording[['fly', col]].rename(columns={col: 'v'})
    d = d[np.isfinite(d['v']) & (d['v'] > 0)]
    return fly_ratio(pd.DataFrame(dict(fly=d['fly'], log_ratio=np.log(d['v']))))


# ---------------------------------------------------------------------------
# the model: APL's activity, the indicator, the comparison

E90_EXPERIMENT = 'E90_model_stages'


def e90_dir(tag: str) -> Path:
    """The folder holding E90's stage run `tag` (`E90_model_stages.load_stage`'s order: the local outputs first, then
    the model-output archive)."""
    for d in (H.HERE / 'outputs' / E90_EXPERIMENT / 'B', H.shared_dir() / 'data' / 'model_outputs' / E90_EXPERIMENT / 'B'):
        if (d / f'{tag}_apl.csv').exists():
            return d
    raise FileNotFoundError(f"E90's stage run {tag} not found")


def input_totals(orn: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per odor the lab input's total over its glomeruli: the sum of the ORN deltas ('input_total') and of their
    positive parts ('input_total_pos')."""
    orn = H.apl16_orn_deltas() if orn is None else orn
    t = pd.DataFrame(dict(input_total=orn.sum(axis=0), input_total_pos=orn.clip(lower=0).sum(axis=0)))
    t.index = [H.odor_abbrev(c) for c in orn.columns]
    return t


def family_apl(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """A4: per wiring and odor, from E90's stage runs (APL on): APL's window mean (a), pre-odor mean (a0), evoked
    activity (a - a0), the PN stage's evoked total (sum over glomeruli of the window minus the pre-odor means), the
    input's total (sum over the 36 glomeruli of the lab's ORN deltas) and the KCs' response rates with APL on and
    acutely off (two or more spikes, the runs' responses)."""
    tot = input_totals()
    rows = []
    for w in wirings:
        d = e90_dir(f'{w}_on')
        apl = pd.read_csv(d / f'{w}_on_apl.csv', index_col=0)
        pn = pd.read_parquet(d / f'{w}_on_pn_win.parquet') - pd.read_parquet(d / f'{w}_on_pn_pre.parquet')
        rate = {c: H.abbrev_columns(pd.read_parquet(e90_dir(f'{w}_{cc}') / f'{w}_{cc}_responses.parquet').astype(float)).mean()
                for c, cc in (('on', 'on'), ('off', 'acute_off'))}
        for o in apl.index:
            rows.append(dict(wiring=w, odor=str(o), a=float(apl.loc[o, 'win']), a0=float(apl.loc[o, 'pre']),
                             evoked=float(apl.loc[o, 'win'] - apl.loc[o, 'pre']), pn_total=float(pn[str(o)].sum()),
                             input_total=float(tot.loc[o, 'input_total']), input_total_pos=float(tot.loc[o, 'input_total_pos']),
                             kc_rate_on=float(rate['on'][o]), kc_rate_off=float(rate['off'][o]), source=str(d)))
    return pd.DataFrame(rows)


def profile_summary(apl: pd.DataFrame) -> pd.DataFrame:
    """A4's profile per wiring: APL's mean window activity, its CV over the odors (sample SD, as E95's A4), the pre-odor
    share, the evoked activity's CV and range, ms's evoked activity over the odors' mean, and the evoked activity's
    Spearman ρ and Pearson r with the input's total, the PN stage's total and the KCs' response rates."""
    from scipy.stats import pearsonr, spearmanr
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        e = g['evoked']
        row = dict(wiring=w, a_mean=float(g['a'].mean()), a_cv=float(g['a'].std() / g['a'].mean()),
                   a0=float(g['a0'].mean()), a0_share=float(g['a0'].mean() / g['a'].mean()), evoked_mean=float(e.mean()),
                   evoked_cv=float(e.std() / e.mean()), evoked_min=float(e.min()), evoked_min_odor=str(e.idxmin()),
                   evoked_max=float(e.max()), evoked_max_odor=str(e.idxmax()), ms_over_mean=float(e['ms'] / e.mean()),
                   ms_rank_from_lowest=int(e.rank()['ms']))
        for col in ('input_total', 'input_total_pos', 'pn_total', 'kc_rate_on', 'kc_rate_off'):
            row[f'spearman_evoked_{col}'] = float(spearmanr(e, g[col])[0])
            row[f'pearson_evoked_{col}'] = float(pearsonr(e, g[col])[0])
        row['pn_total_cv'] = float(g['pn_total'].std(ddof=0) / g['pn_total'].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def indicator_ratio(a1, a2, a0, n: float, form: str):
    """The model's ratio of stimulus 1 over stimulus 2 seen through an indicator of exponent `n`: 'baseline',
    (a1^n - a0^n) / (a2^n - a0^n); 'evoked', ((a1 - a0) / (a2 - a0))^n. Equal to the evoked ratio at n 1."""
    a1, a2, a0 = (np.asarray(x, dtype=float) for x in (a1, a2, a0))
    if form == 'baseline':
        return (a1 ** n - a0 ** n) / (a2 ** n - a0 ** n)
    if form == 'evoked':
        return ((a1 - a0) / (a2 - a0)) ** n
    raise ValueError(f'unknown form {form!r}')


def model_ratio_set(a1: float, a2: float, a0: float, ns: Sequence[float] = INDICATOR_NS,
                    forms: Sequence[str] = FORMS) -> pd.DataFrame:
    """One wiring's ratios over `ns` and `forms` (columns n, form, ratio)."""
    return pd.DataFrame([dict(n=float(n), form=f, ratio=float(indicator_ratio(a1, a2, a0, n, f))) for f in forms for n in ns])


def classify(ratios: Sequence[float], lo: float, hi: float) -> str:
    """One wiring against the flies' interval [lo, hi]: 'inside' if any ratio lies in it; 'flatter' if all are nearer 1
    than the interval or on the other side of 1; 'steeper' if all are beyond the interval's far end (from 1); 'mixed'
    otherwise (an interval that holds 1 has no 'flatter': nothing is nearer 1 than it)."""
    r = np.asarray(ratios, dtype=float)
    if not np.all(np.isfinite(r) & (r > 0)):
        raise ValueError(f'ratios must be positive and finite: {ratios}')
    x = np.log(r)
    L, U = np.log(float(lo)), np.log(float(hi))
    if not L <= U:
        raise ValueError(f'empty interval [{lo}, {hi}]')
    if np.any((x >= L) & (x <= U)):
        return 'inside'
    if L > 0:
        return 'flatter' if np.all(x < L) else 'steeper' if np.all(x > U) else 'mixed'
    if U < 0:
        return 'flatter' if np.all(x > U) else 'steeper' if np.all(x < L) else 'mixed'
    return 'steeper' if (np.all(x > U) or np.all(x < L)) else 'mixed'


def comparison_verdict(classes: Sequence[str], min_wirings: int = VERDICT_MIN_WIRINGS) -> str:
    """'inside', 'flatter' or 'steeper' if at least `min_wirings` of the wirings are; otherwise 'not decided'."""
    c = list(classes)
    for v in ('inside', 'flatter', 'steeper'):
        if sum(x == v for x in c) >= min_wirings:
            return v
    return 'not decided'


def hull(intervals: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    """The smallest interval holding all of `intervals` (V1: the two readings of paa's concentration)."""
    iv = [(float(a), float(b)) for a, b in intervals]
    if not iv or any(not (np.isfinite(a) and np.isfinite(b) and a <= b) for a, b in iv):
        raise ValueError(f'bad intervals {intervals}')
    return min(a for a, _ in iv), max(b for _, b in iv)


def v1_verdict(per_estimate: Dict[str, str]) -> str:
    """V1 from its comparison verdict under each estimate of ms's input: 'inside' if inside under any; 'flatter' or
    'steeper' only if that under all; otherwise 'not decided'."""
    v = list(per_estimate.values())
    if not v:
        raise ValueError('no estimate of ms\'s input')
    if 'inside' in v:
        return 'inside'
    for d in ('flatter', 'steeper'):
        if all(x == d for x in v):
            return d
    return 'not decided'


def wiring_classes(apl: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring the model's ratio of odor `num` over odor `den` (the panel's runs) through the indicator family, its
    range, the linear reading (n 1, evoked) and its class against [lo, hi]."""
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        a0 = float(g['a0'].iloc[0])
        if not np.allclose(g['a0'], a0):
            raise RuntimeError(f'{w}: the pre-odor level differs between odors')
        rs = model_ratio_set(g.loc[num, 'a'], g.loc[den, 'a'], a0)
        rows.append(dict(wiring=w, ratio_linear=float(indicator_ratio(g.loc[num, 'a'], g.loc[den, 'a'], a0, 1.0, 'evoked')),
                         ratio_min=float(rs['ratio'].min()), ratio_max=float(rs['ratio'].max()),
                         ratio_n2p25_baseline=float(rs[(rs.n == 2.25) & (rs.form == 'baseline')]['ratio'].iloc[0]),
                         ratio_n2p25_evoked=float(rs[(rs.n == 2.25) & (rs.form == 'evoked')]['ratio'].iloc[0]),
                         cls=classify(rs['ratio'], lo, hi), cls_linear=classify([rs['ratio'].iloc[0]], lo, hi)))
    return pd.DataFrame(rows)


def kc_rates(criterion: str = 'z1') -> pd.DataFrame:
    """The KC data's per-odor response rates (columns 'ctrl' and 'tnt', rows the 16 odors and 'pfo' where available):
    'z1' the mean over recordings of the fraction above Yang's 1 z; 'noise' the registry's rates in baseline-noise
    units (noise_ctrl_frac_<odor>, noise_tnt_frac_<odor>; the bands' criterion)."""
    if criterion == 'z1':
        return pd.DataFrame({k: H.yang_response_fractions(g, thresh=1.0).mean() for k, g in (('ctrl', 'TNT3c'), ('tnt', 'TNT'))})
    if criterion == 'noise':
        return pd.DataFrame({k: pd.Series({o: H.data_target(f'noise_{k}_frac_{o}') for o in H.APL16_ODORS})
                             for k in ('ctrl', 'tnt')})
    raise ValueError(f'unknown criterion {criterion!r}')


ACROSS_ODORS = (('paa_over_ms', 'paa', 'ms'), ('pa_over_ms', 'pa', 'ms'), ('pa_over_paa', 'pa', 'paa'))
E37_TABLE = ('E37_spike_only_apl', 'apl_activity.csv')


def alternative_apls(apl: pd.DataFrame) -> pd.DataFrame:
    """A5: the ratios of the across-odor comparisons (ACROSS_ODORS, the panel's concentrations) by other APLs: the
    family's (evoked, and the window mean; median and range over the wirings of `apl`); in proportion to the KCs'
    response rates in each genotype (Yang's criterion of 1 z, and the noise units of the bands beside); in proportion to
    the PN dendrites' totals and to the pebbled ORNs' totals (the megamat odors, the 28 glomeruli they share with the
    model; E90's `al_tables` / `al_statistics`; paa is not megamat); E37's saved table (window mean and peak of every
    variant, no pre-odor level)."""
    rows = []

    def add(source: str, kind: str, values: Dict[str, float], note: str = '', lo=None, hi=None):
        for rid, a, b in ACROSS_ODORS:
            va, vb = values.get(a, np.nan), values.get(b, np.nan)
            r = va / vb if np.isfinite(va) and np.isfinite(vb) and vb > 0 else np.nan
            rows.append(dict(source=source, kind=kind, ratio_id=rid, ratio=float(r),
                             ratio_min=float(lo[rid]) if lo is not None else np.nan,
                             ratio_max=float(hi[rid]) if hi is not None else np.nan, note=note))

    for col, kind in (('evoked', 'family APL, evoked (a - a0), median over wirings'),
                      ('a', 'family APL, window mean (a), median over wirings')):
        per = {}
        for w, g in apl.groupby('wiring', sort=False):
            s = g.set_index('odor')[col]
            per[w] = {rid: s[a] / s[b] for rid, a, b in ACROSS_ODORS}
        per = pd.DataFrame(per).T
        med = per.median()
        rows += [dict(source='family', kind=kind, ratio_id=rid, ratio=float(med[rid]), ratio_min=float(per[rid].min()),
                      ratio_max=float(per[rid].max()), note='8 wirings, E90 stage runs') for rid, _, _ in ACROSS_ODORS]
    for crit in ('z1', 'noise'):
        k = kc_rates(crit)
        for g, lab in (('ctrl', 'control flies'), ('tnt', 'APL output blocked (TNT)')):
            add(f'kc_rate_{crit}', f'in proportion to the KCs\' response rates, {lab} ({crit})', k[g].to_dict(),
                note='KC panel concentrations (pa, ms, paa at 10^-3)')
    from experiments import E90_model_stages as E90
    st = E90.al_statistics(**{k: v for k, v in E90.al_tables().items()})
    add('pn_dendrites', 'in proportion to the PN dendrites\' summed response (GH146, 28 glomeruli)', st['total_pn'].to_dict(),
        note='megamat odors only (paa not imaged)')
    add('orn_pebbled', 'in proportion to the ORNs\' summed response (pebbled, 28 glomeruli)', st['total_orn'].to_dict(),
        note='megamat odors only')
    tot = input_totals()
    add('input', "in proportion to the model's input total (36 glomeruli)", tot['input_total'].to_dict(),
        note='the lab input; paa at 10^-4')
    e37 = pd.read_csv(H.shared_dir() / 'results' / E37_TABLE[0] / E37_TABLE[1], index_col=0)
    for c in e37.columns:
        add(f'E37:{c}', 'E37 saved APL (time >= 0; mean or peak, no pre-odor level)', e37[c].to_dict(),
            note='older bases fitted at 5 % under the structure loss')
    return pd.DataFrame(rows)


def power_sds(a3_main: pd.DataFrame) -> dict:
    """A6's between-fly SDs of a log2 ratio: the median and the largest over A3's across-odor ratios (numerator and
    denominator different odors) whose two concentrations are both POWER_MIN_CONC or above (the panel's odors are at
    about 10^-3; ratios of weak responses spread more), under one reading and subset (`a3_main`), with the ratios used."""
    ok = []
    for rid, num, den, _ in RATIOS:
        if num[0] != den[0] and min(num[1], den[1]) >= POWER_MIN_CONC - CONC_TOL:
            ok.append(rid)
    x = a3_main.set_index('ratio_id').loc[ok, 'sd_log2_flies'].astype(float)
    return dict(median=float(x.median()), max=float(x.max()), ratios={k: float(v) for k, v in x.items()})


def flies_needed(delta_log2: float, sd_log2: float, z: Tuple[float, float] = POWER_Z) -> float:
    """A6: the number of flies whose geometric mean's 90 % interval excludes a value `delta_log2` (log2) away with 80 %
    power, at a between-fly SD `sd_log2` (log2): ceil(((z_a + z_b) sd / delta)^2), at least 2; inf when delta is 0."""
    d = abs(float(delta_log2))
    if d == 0:
        return float('inf')
    return float(max(2, int(np.ceil(((z[0] + z[1]) * float(sd_log2) / d) ** 2 - 1e-9))))


SUBSTITUTED_INPUT = ('paa', 'euc', '-aPine')    # the model's ORN input at another concentration than the KC panel's


def power_table(apl: pd.DataFrame, sd_log2: Dict[str, float]) -> pd.DataFrame:
    """A6: for every pair of panel odors (oriented so that the proportional APL's ratio is at least 1) the family's
    evoked ratio (median and range over the wirings), the proportional APL's (the KCs' response rates with APL output
    blocked (TNT), control flies beside; Yang's 1 z), their log2 difference (and that from the wiring nearest the
    proportional APL) and the flies needed at each between-fly SD of `sd_log2`; sorted by the size of the difference,
    so the pairs that need the fewest flies come first. The ratios are linear readings (n 1, evoked): through an
    indicator of n above 1 the family's stays near 1 and the proportional APL's moves further away, so the counts are
    upper bounds in that respect."""
    from itertools import combinations
    k = kc_rates('z1')
    ev = apl.pivot(index='odor', columns='wiring', values='evoked')
    rows = []
    for a, b in combinations(list(H.APL16_ODORS), 2):
        num, den = (a, b) if k.loc[a, 'tnt'] >= k.loc[b, 'tnt'] else (b, a)
        fam = ev.loc[num] / ev.loc[den]
        prop = float(k.loc[num, 'tnt'] / k.loc[den, 'tnt'])
        fmed = float(fam.median())
        d = float(np.log2(prop) - np.log2(fmed))
        row = dict(numerator=num, denominator=den, substituted_input=bool({num, den} & set(SUBSTITUTED_INPUT)), family_ratio=fmed, family_min=float(fam.min()), family_max=float(fam.max()),
                   proportional_tnt=prop, proportional_ctrl=float(k.loc[num, 'ctrl'] / k.loc[den, 'ctrl']),
                   delta_log2=d, delta_log2_nearest_wiring=float(np.min(np.abs(np.log2(prop) - np.log2(fam.to_numpy(float))))))
        for name, sd in sd_log2.items():
            row[f'flies_needed_sd_{name}'] = flies_needed(d, sd)
            row[f'flies_needed_nearest_wiring_sd_{name}'] = flies_needed(row['delta_log2_nearest_wiring'], sd)
        rows.append(row)
    df = pd.DataFrame(rows)
    return df.reindex(df['delta_log2'].abs().sort_values(ascending=False, kind='stable').index).reset_index(drop=True)


# ---------------------------------------------------------------------------
# part B: the extended input and its runs (phase 2)

B_ARMS = ('B1', 'B2', 'B3', 'B4')


def b1_input(label: str = 'real') -> Tuple[pd.DataFrame, dict]:
    """B1's extended input (the lab's 16 odors, then the added stimuli, on the input's 36 glomeruli) and its record
    (every added stimulus's source session, its column and the glomeruli filled with 0)."""
    raise NotImplementedError


def hallem_columns(glomeruli: Sequence[str]) -> pd.DataFrame:
    """Hallem's pa and ms at 10^-2, -4, -6 and -8 and paa at 10^-2 from drosolf's concentration rows, onto the input's
    glomeruli by E38's conversion (glomeruli neither has filled with 0)."""
    raise NotImplementedError


class ExtendedRun:
    """One wiring at its fit (E83's `build_simulator`, breadth held at the lab input's) on B1's extended input, APL's
    activity saved; `run(cond)` simulates APL on or acutely off and returns the outputs and APL's window and pre-odor
    means per stimulus."""

    def __init__(self, label: str):
        raise NotImplementedError

    def run(self, cond: str) -> dict:
        raise NotImplementedError


def run_identity(label: str) -> dict:
    """Part 0's run identities on one wiring: the 16 panel odors' counts equal the saved fit's (APL on and off), APL's
    window and pre-odor means E90's saved ones, the blank's evoked APL zero."""
    raise NotImplementedError


def part_b(arms: Sequence[str], wirings: Sequence[str]) -> None:
    """Part B's runs (resumable), saved under outputs/E98_apl_response_strength/B/."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

E79_RESULTS = ('results', 'E79_apl_calyx_constraints')
IDENTITY_TOL = 1e-12


def _max_diff(a: pd.Series, b: pd.Series) -> Tuple[float, int]:
    """The largest |a - b| over their common index and the number of index entries in only one of them."""
    common = a.index.intersection(b.index)
    unmatched = int(len(a.index.symmetric_difference(b.index)))
    d = (a.loc[common].astype(float) - b.loc[common].astype(float)).abs()
    return (float(d.max()) if len(d) else np.nan), unmatched


def reader_identity(recs=None, resp: Optional[pd.DataFrame] = None) -> dict:
    """Part 0's reader identities: E79's per-stimulus calyx means under its three baselines
    (E79_maps_across_baselines.csv: mean_order_matched, mean_first_in_block, mean_all_pfo) and its blank table
    (E79_baseline_noise.csv: the blank's and the paraffin oil presentations' calyx means per recording;
    E79_blank_vs_odor_maps.csv: the blank and each stimulus's first-in-block mean per recording and stimulus)."""
    from experiments import E79_apl_calyx_constraints as E79
    recs = load_recordings() if recs is None else recs
    resp = stimulus_responses(recs) if resp is None else resp
    e79 = H.shared_dir().joinpath(*E79_RESULTS)
    out: Dict[str, object] = {}
    maps = pd.read_csv(e79 / 'E79_maps_across_baselines.csv').set_index(['recording', 'stimulus'])
    raw = resp[~resp['blank_subtracted'].astype(bool)]
    for bl in BASELINES:
        ours = raw[raw['baseline'] == bl].set_index(['recording', 'stimulus'])['response']
        d, u = _max_diff(ours, maps[f'mean_{bl}'])
        out[f'maps_{bl}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(ours)))
    noise = pd.read_csv(e79 / 'E79_baseline_noise.csv').set_index('recording')
    rows = {}
    for rec in recs:
        bf = E79.blank_frames(rec.frames)
        m, st = rec.mask, rec.stack
        rows[rec.name] = dict(blank_mean=float(E79.position_matched_blank(st, rec.frames)[m].mean()),
                              pfo_first_in_block_mean=float(st[bf['first_in_block']].mean(axis=0)[m].mean()),
                              pfo_block_pos2_mean=float(st[bf['pfo_pos2'][0]][m].mean()),
                              pfo_block_pos3_mean=float(st[bf['pfo_pos3'][0]][m].mean()))
    mine = pd.DataFrame(rows).T
    for col in mine.columns:
        d, u = _max_diff(mine[col], noise[col])
        out[f'blank_table_{col}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(mine)))
    bvo = pd.read_csv(e79 / 'E79_blank_vs_odor_maps.csv').set_index(['recording', 'stimulus'])
    fib = raw[raw['baseline'] == 'first_in_block'].drop_duplicates(['recording', 'stimulus']).set_index(['recording', 'stimulus'])
    for ours, theirs in (('blank', 'blank_mean'), ('response', 'odor_mean_first_in_block')):
        d, u = _max_diff(fib[ours], bvo[theirs])
        out[f'blank_vs_odor_{theirs}'] = dict(max_abs_diff=d, unmatched=u, n=int(len(fib)))
    out['exact'] = bool(all(v['unmatched'] == 0 and np.isfinite(v['max_abs_diff']) and v['max_abs_diff'] <= IDENTITY_TOL
                            for v in out.values() if isinstance(v, dict)))
    out['tolerance'] = IDENTITY_TOL
    return out


def identity_passed(which: str = 'reader') -> bool:
    """Whether identity.json records the reader's ('reader') or every ('all') identity as passed."""
    d = _read_json(out_dir() / 'identity.json') or {}
    if which == 'reader':
        return bool(d.get('reader', {}).get('exact'))
    if which == 'all':
        return bool(d.get('all_exact'))
    raise ValueError(f'unknown identity {which!r}')


def part_0(wirings: Sequence[str] = ALL_WIRINGS, runs: bool = True) -> dict:
    """Part 0 (module docstring): identity.json. The reader's identities always; the run identities on `wirings` when
    `runs` (part B's builder)."""
    path = out_dir() / 'identity.json'
    res = _read_json(path) or {}
    t0 = time.time()
    res['reader'] = dict(reader_identity(), time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(path, res)
    print(f"reader identity: exact {res['reader']['exact']} ({time.time() - t0:.0f} s)", flush=True)
    if not res['reader']['exact']:
        raise RuntimeError(f"the reader does not reproduce E79: {res['reader']}")
    if runs:
        for w in wirings:
            res[f'run_{w}'] = run_identity(w)
            _write_json(path, res)
    res['all_exact'] = bool(res['reader']['exact'] and all(res.get(f'run_{w}', {}).get('exact') for w in ALL_WIRINGS))
    _write_json(path, res)
    return res


def a1_tables(recs, resp: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """A1: the presentations' agreement per recording (under every reading: Pearson r over stimuli of the two
    presentations' calyx means, the median |rep1 - rep2| / 2, the median ratio rep2 / rep1 of responses above
    MIN_RESPONSE) and the drift of the first-in-block paraffin oil response (per recording and block its calyx mean, the
    block's stimulus; per recording the slope over blocks, the first and last; the Spearman ρ of block with concentration
    within each odor, the order's confound with concentration)."""
    from scipy.stats import spearmanr
    agree = []
    for (rec, bl, sub), g in resp.groupby(['recording', 'baseline', 'blank_subtracted']):
        ok = (g['rep1'] > MIN_RESPONSE) & (g['rep2'] > MIN_RESPONSE)
        agree.append(dict(recording=rec, fly=g['fly'].iloc[0], pulse_s=float(g['pulse_s'].iloc[0]),
                          reading=reading_label((bl, bool(sub))), n_stimuli=int(len(g)),
                          r_rep1_rep2=float(np.corrcoef(g['rep1'], g['rep2'])[0, 1]), presentation_se=presentation_se(g),
                          median_rep2_over_rep1=float(np.median(g.loc[ok, 'rep2'] / g.loc[ok, 'rep1'])) if ok.any() else np.nan,
                          mean_rep2_minus_rep1=float((g['rep2'] - g['rep1']).mean()), blank=float(g['blank'].iloc[0])))
    drift_rows, drift_sum = [], []
    for rec in recs:
        fr = rec.frames
        pfo1 = fr[(fr['odor'] == 'pfo') & (fr['position'] == 1)]
        vals = np.array([float(rec.stack[i][rec.mask].mean()) for i in pfo1.index])
        blocks = pfo1['block'].to_numpy(int)
        stim = {int(b): '; '.join(sorted({f"{o} {c:g}" for o, c in zip(g['odor'], g['log10_conc']) if o != 'pfo'}) or ['pfo'])
                for b, g in fr.groupby('block')}
        for b, v in zip(blocks, vals):
            drift_rows.append(dict(recording=rec.name, block=int(b), stimulus=stim[int(b)], pfo_first_in_block=v))
        r = resp[(resp['recording'] == rec.name) & (resp['baseline'] == 'order_matched') & ~resp['blank_subtracted'].astype(bool)]
        order = {f'spearman_block_conc_{o}': float(spearmanr(g['block'], g['conc'])[0]) for o, g in r.groupby('odor')}
        drift_sum.append(dict(recording=rec.name, n_blocks=int(len(blocks)), pfo_first_mean=float(vals.mean()),
                              pfo_first_sd=float(vals.std(ddof=1)), pfo_first_first=float(vals[0]), pfo_first_last=float(vals[-1]),
                              slope_per_block=float(np.polyfit(blocks, vals, 1)[0]),
                              spearman_block=float(spearmanr(blocks, vals)[0]), **order))
    return dict(agreement=pd.DataFrame(agree), drift=pd.DataFrame(drift_rows), drift_summary=pd.DataFrame(drift_sum))


def a2_medians(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING) -> pd.DataFrame:
    """Per odor and imaged concentration the median, minimum and maximum over recordings of the response, and the
    geometric mean over flies (responses above MIN_RESPONSE) with its 90 % interval."""
    rows = []
    for (o, c), g in _reading(resp, reading).groupby(['odor', 'conc']):
        gm = fly_geomean(g, 'response')
        rows.append(dict(odor=o, conc=float(c), reading=reading_label(reading), n_recordings=int(len(g)),
                         median=float(g['response'].median()), min=float(g['response'].min()), max=float(g['response'].max()),
                         n_above=int((g['response'] > MIN_RESPONSE).sum()), fly_geomean=gm['ratio'], fly_lo90=gm['lo90'],
                         fly_hi90=gm['hi90'], n_flies=gm['n_flies']))
    return pd.DataFrame(rows)


def a2_fly_summary(dose: pd.DataFrame) -> pd.DataFrame:
    """Per odor the fold change per log unit over flies (geometric mean of the recordings' folds within flies, then over
    flies, with the 90 % interval) and the lowest concentration above the blank (median and range over recordings)."""
    rows = []
    for (o, rd), g in dose.groupby(['odor', 'reading']):
        gm = fly_geomean(g, 'fold_per_log_unit')
        rows.append(dict(odor=o, reading=rd, fold_per_log_unit=gm['ratio'], lo90=gm['lo90'], hi90=gm['hi90'],
                         n_flies=gm['n_flies'], fold_min=float(g['fold_per_log_unit'].min()),
                         fold_max=float(g['fold_per_log_unit'].max()),
                         lowest_above_blank_median=float(g['lowest_above_blank'].median()),
                         lowest_above_blank_min=float(g['lowest_above_blank'].min()),
                         lowest_above_blank_max=float(g['lowest_above_blank'].max()),
                         n_recordings_none_above=int(g['lowest_above_blank'].isna().sum())))
    return pd.DataFrame(rows)


def a4_model_ratios(apl: pd.DataFrame) -> pd.DataFrame:
    """The model's ratios that the panel's runs give, per wiring, n and form: paa over ms (the lab's input, paa at
    10^-4 and ms at 10^-3: V1 under the megamat estimate of ms's input) and pa over ms (pa at 10^-3; beside the flies'
    pa at 10^-4.386, outside the one-log rule)."""
    rows = []
    for w, g in apl.groupby('wiring', sort=False):
        g = g.set_index('odor')
        a0 = float(g['a0'].iloc[0])
        for rid, a, b in ACROSS_ODORS:
            rs = model_ratio_set(g.loc[a, 'a'], g.loc[b, 'a'], a0)
            for _, r in rs.iterrows():
                rows.append(dict(wiring=w, ratio_id=rid, n=r['n'], form=r['form'], ratio=r['ratio']))
    return pd.DataFrame(rows)


PREVIEW = (('V1 under the megamat estimate of ms\'s input', 'paa', 'ms', ('V1_paa-4', 'V1_paa-2.954')),
           ('pa over ms (the model\'s pa at 10^-3, the flies\' at 10^-4.386, outside the one-log rule)', 'pa', 'ms',
            ('probe_pa-4.386',)))


PREVIEW_SUBSETS = ('all', 'all_recording_ranges')


def a4_preview(apl: pd.DataFrame, a3: pd.DataFrame) -> pd.DataFrame:
    """The comparisons that the panel's runs already allow, under every reading of the flies (subsets PREVIEW_SUBSETS): per wiring
    the class (`classify`, n and both forms) and the linear reading's, and the comparison's verdict. V1 here is under one
    of ms's estimates only: inside here means V1 is inside; outside here leaves V1 to part C."""
    rows = []
    a3 = a3[a3['subset'].isin(PREVIEW_SUBSETS)]
    for label, num, den, ids in PREVIEW:
        for (sub, rd), g in a3.groupby(['subset', 'reading'], sort=False):
            g = g.set_index('ratio_id')
            iv = [(g.loc[i, 'lo90'], g.loc[i, 'hi90']) for i in ids]
            if any(not np.isfinite(x) for v in iv for x in v):
                continue
            lo, hi = hull(iv)
            wc = wiring_classes(apl, num, den, lo, hi)
            rows.append(dict(comparison=label, reading=rd, subset=sub, flies_lo90=lo, flies_hi90=hi,
                             verdict=comparison_verdict(wc['cls']), verdict_linear=comparison_verdict(wc['cls_linear']),
                             **{f'n_{c}': int((wc['cls'] == c).sum()) for c in ('inside', 'flatter', 'steeper', 'mixed')},
                             model_linear_min=float(wc['ratio_linear'].min()), model_linear_max=float(wc['ratio_linear'].max()),
                             model_min=float(wc['ratio_min'].min()), model_max=float(wc['ratio_max'].max()),
                             per_wiring='; '.join(f"{r.wiring} {r.cls} ({r.ratio_linear:.3f})" for r in wc.itertuples())))
    return pd.DataFrame(rows)


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def input_files(wirings: Sequence[str] = ALL_WIRINGS) -> Dict[str, str]:
    """Part A's inputs (path relative to the shared directory or the repository root -> sha256): the calyx maps (TIFFs,
    ROIs, recordings.csv, frames.csv), E79's tables the reader is checked against, E90's stage files read, E37's APL
    table and the data-target registry."""
    from experiments import E79_apl_calyx_constraints as E79
    sh = H.shared_dir()
    files = sorted(E79.DATA.glob('*.tif')) + sorted(E79.DATA.glob('*.roi')) + [E79.DATA / 'recordings.csv', E79.DATA / 'frames.csv']
    e79 = sh.joinpath(*E79_RESULTS)
    files += [e79 / f for f in ('E79_maps_across_baselines.csv', 'E79_baseline_noise.csv', 'E79_blank_vs_odor_maps.csv')]
    for w in wirings:
        files += [e90_dir(f'{w}_on') / f'{w}_on_{k}' for k in ('apl.csv', 'pn_win.parquet', 'pn_pre.parquet', 'responses.parquet')]
        files.append(e90_dir(f'{w}_acute_off') / f'{w}_acute_off_responses.parquet')
    files += [sh / 'results' / E37_TABLE[0] / E37_TABLE[1], H.data_targets_path()]
    out = {}
    for f in files:
        f = Path(f)
        key = str(f.relative_to(sh)) if f.is_relative_to(sh) else str(f.relative_to(H.HERE)) if f.is_relative_to(H.HERE) else str(f)
        out[key] = _sha256(f)
    return out


def _records(df: pd.DataFrame, cols: Sequence[str]) -> List[dict]:
    return [{c: (r[c] if not isinstance(r[c], (np.floating, np.integer)) else r[c].item()) for c in cols} for _, r in df.iterrows()]


def part_a() -> dict:
    """Part A (module docstring): A1 to A6 and a_summary.json."""
    if not identity_passed('reader'):
        raise RuntimeError("part A needs part 0's reader identities (run --part 0 --no-runs first)")
    out = out_dir()
    t0 = time.time()
    recs = load_recordings()
    resp = stimulus_responses(recs)
    resp.to_csv(out / 'a1_responses.csv', index=False)
    a1 = a1_tables(recs, resp)
    a1['agreement'].to_csv(out / 'a1_agreement.csv', index=False)
    a1['drift'].to_csv(out / 'a1_pfo_drift.csv', index=False)
    a1['drift_summary'].to_csv(out / 'a1_pfo_drift_summary.csv', index=False)
    print(f'A1 done ({time.time() - t0:.0f} s)', flush=True)
    dose, at, med = [], [], []
    for rd in READINGS:
        d, a = dose_response(resp, rd)
        dose.append(d)
        at.append(a)
        med.append(a2_medians(resp, rd))
    dose, at, med = (pd.concat(x, ignore_index=True) for x in (dose, at, med))
    dose.to_csv(out / 'a2_dose_response.csv', index=False)
    at.to_csv(out / 'a2_at_model_concs.csv', index=False)
    med.to_csv(out / 'a2_medians.csv', index=False)
    a2s = a2_fly_summary(dose)
    a2s.to_csv(out / 'a2_fly_summary.csv', index=False)
    print(f'A2 done ({time.time() - t0:.0f} s)', flush=True)
    a3, a3r = [], []
    for rd in READINGS:
        for sub in SUBSETS:
            s, r = flies_ratios(resp, rd, sub)
            a3.append(s)
            a3r.append(r)
    a3, a3r = pd.concat(a3, ignore_index=True), pd.concat(a3r, ignore_index=True)
    a3.to_csv(out / 'a3_ratios.csv', index=False)
    a3r.to_csv(out / 'a3_per_recording.csv', index=False)
    print(f'A3 done ({time.time() - t0:.0f} s)', flush=True)
    apl = family_apl()
    apl.to_csv(out / 'a4_apl.csv', index=False)
    prof = profile_summary(apl)
    prof.to_csv(out / 'a4_profile.csv', index=False)
    mr = a4_model_ratios(apl)
    mr.to_csv(out / 'a4_model_ratios.csv', index=False)
    prev = a4_preview(apl, a3)
    prev.to_csv(out / 'a4_preview.csv', index=False)
    print(f'A4 done ({time.time() - t0:.0f} s)', flush=True)
    alt = alternative_apls(apl)
    alt.to_csv(out / 'a5_alternatives.csv', index=False)
    print(f'A5 done ({time.time() - t0:.0f} s)', flush=True)
    main = a3[(a3['reading'] == reading_label(MAIN_READING)) & (a3['subset'] == 'all')]
    probe_rd = a3[(a3['reading'] == reading_label(('order_matched', False))) & (a3['subset'] == 'all')]
    sds = power_sds(main)
    power = power_table(apl, {k: v for k, v in sds.items() if k in ('median', 'max')})
    power.to_csv(out / 'a6_power.csv', index=False)
    print(f'A6 done ({time.time() - t0:.0f} s)', flush=True)
    ratio_cols = ['ratio_id', 'numerator', 'denominator', 'ratio', 'lo90', 'hi90', 'sd_log2_flies', 'n_flies',
                  'n_recordings', 'n_recordings_left_out', 'num_how', 'den_how', 'probe_ratio', 'probe_lo90', 'probe_hi90']
    summ = dict(
        reader_identity=(_read_json(out / 'identity.json') or {}).get('reader', {}).get('exact'),
        blank_identity_max_abs_diff=blank_identity(resp),
        large_blank_recordings=sorted(resp.loc[resp['blank'] > LARGE_BLANK, 'recording'].unique()),
        a3_main=_records(main if len(main) else a3[a3['reading'] == reading_label(MAIN_READING)], ratio_cols),
        a3_order_matched_no_blank=_records(probe_rd, ratio_cols),
        a4_profile=_records(prof, ['wiring', 'a_mean', 'a_cv', 'a0_share', 'evoked_mean', 'evoked_cv', 'evoked_min_odor',
                                   'evoked_max_odor', 'ms_over_mean', 'spearman_evoked_input_total',
                                   'spearman_evoked_pn_total', 'spearman_evoked_kc_rate_off']),
        a4_preview=_records(prev, ['comparison', 'reading', 'subset', 'flies_lo90', 'flies_hi90', 'verdict', 'verdict_linear',
                                   'n_inside', 'n_flatter', 'n_steeper', 'model_linear_min', 'model_linear_max']),
        a5=_records(alt[~alt['source'].str.startswith('E37:') | alt['source'].str.contains('_mean')],
                    ['source', 'kind', 'ratio_id', 'ratio', 'ratio_min', 'ratio_max']),
        a6_sd_log2=sds, a6_best5=_records(power.head(5), list(power.columns)),
        a6_best5_without_substituted_input=_records(power[~power['substituted_input']].head(5), list(power.columns)),
        a6_n_pairs_flies_needed_at_most_3=dict(
            median=int((power['flies_needed_sd_median'] <= 3).sum()), max=int((power['flies_needed_sd_max'] <= 3).sum())),
        runtime_s=time.time() - t0)
    _write_json(out / 'a_summary.json', summ)
    _write_json(out / 'provenance.json', dict(experiment=EXPERIMENT, script=H.script_path(__file__), parts=['0', 'A'],
        host=_host(), time=time.strftime('%Y-%m-%d %H:%M'), commits=dict(H.code_commits()), inputs=input_files(),
        settings=dict(MIN_RESPONSE=MIN_RESPONSE, MATCH_LOG_UNITS=MATCH_LOG_UNITS, N_BOOT=N_BOOT, SEED=SEED, CI=CI,
                      INDICATOR_NS=INDICATOR_NS, FORMS=FORMS, VERDICT_MIN_WIRINGS=VERDICT_MIN_WIRINGS,
                      MAIN_READING=MAIN_READING, LARGE_BLANK=LARGE_BLANK, POWER_Z=POWER_Z, POWER_MIN_CONC=POWER_MIN_CONC,
                      RATIOS=[list(r) for r in RATIOS])))
    return summ


def part_c() -> dict:
    """Part C (module docstring): verdict.json."""
    raise NotImplementedError


def comparison_table() -> str:
    """comparison_table.md from the parts' tables."""
    raise NotImplementedError


def _wirings(arg: str) -> Tuple[str, ...]:
    if arg == 'all8':
        return ALL_WIRINGS
    ws = tuple(x.strip() for x in arg.split(',') if x.strip())
    bad = [w for w in ws if w not in ALL_WIRINGS]
    if bad:
        raise SystemExit(f'unknown wirings {bad}')
    return ws


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'C', 'table'])
    ap.add_argument('--wirings', default='all8', help="all8 or a comma-separated list of wiring labels")
    ap.add_argument('--arm', default='B1', help="part B's arms, comma-separated (B1..B4)")
    ap.add_argument('--threads', type=int, default=4,
                    help='part B (and part 0 with runs): OMP threads (OMP_NUM_THREADS and OMP_THREAD_LIMIT must match)')
    ap.add_argument('--no-runs', action='store_true', help="part 0: the reader's identities only")
    a = ap.parse_args()
    ws = _wirings(a.wirings)
    if a.part in ('B',) or (a.part == '0' and not a.no_runs):
        for v in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT'):
            if os.environ.get(v) != str(a.threads):
                ap.error(f'part {a.part}: run with OMP_NUM_THREADS={a.threads} OMP_THREAD_LIMIT={a.threads} in the '
                         'environment')
    if a.part == '0':
        provenance('0', wirings=list(ws), runs=not a.no_runs)
        print(json.dumps(part_0(ws, runs=not a.no_runs), indent=1, default=float), flush=True)
    elif a.part == 'A':
        provenance('A')
        print(json.dumps(part_a(), indent=1, default=float), flush=True)
    elif a.part == 'B':
        arms = tuple(x for x in a.arm.split(',') if x)
        provenance('B', arms=list(arms), wirings=list(ws), threads=a.threads)
        part_b(arms, ws)
    elif a.part == 'C':
        provenance('C')
        print(json.dumps(part_c(), indent=1, default=float), flush=True)
    elif a.part == 'table':
        print(comparison_table(), flush=True)


if __name__ == '__main__':
    main()
