#!/usr/bin/env python3
"""E105: the KCs' response rates at the weak odors read like for like. Whether the best family's KCs respond to weak
input less than the flies' once the flies' three sets of recordings are kept apart, the solvent's responders are removed,
both panels are read in one unit and each session's estimate of the input is read by itself; and what the antennal
lobe's weaker normalization (E102) costs at low concentration on that reading. A diagnostic: no model is changed, nothing
is fitted or scored on the main target.

The entry is item E105 of experiment_priorities.md (proposed by viki 2026-11-17 01:36; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 01:16). The flies: Yang's processed per-trial values in the shared
directory (data/yang_weak_odor_kc/, data/yang_apl16_kc_olfdb/, data/yang_apl16_kc_noise_units/); the family: C2's (the
breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse) on E72 part B's eight wirings at their fits, nothing re-tuned,
read from the saved runs (E98's B1 for the real wiring, E101's archived A3 for the other seven, through E104's
`saved_run`; E102's archived ladder and refits for part B; E104's saved b_model.csv for part C).

The flies, in three sets that are never pooled (`SETS`, `assign_set`):
- D1: the diagnostic-mixture recordings to 2025-02-24 (9 control recordings of 7 flies; the 3 recordings of 2 flies with
  APL output blocked (TNT)); six trials, the solvent three.
- D2: the same panel from 2025-03-20 (7 control recordings of 5 flies); six trials (one recording four).
- N: the dilute natural-mixture recordings of 2026-01-16 to 02-10 (9 control recordings of 5 flies); two to four trials.
- APL16 (11 control recordings of 8 flies, 8 with APL output blocked (TNT) of 7 flies): for the unit and for the
  panel's values.
The reader (`read_trials`) refuses any path outside the shared directory (`check_path`), any recording outside the
panel's processed list (its computed/recordings.csv) and any recording data/yang_apl_recordings/recordings.csv marks raw
only (`check_rows`).

Definitions (the entry's; all in `definitions()` and hashed into the lock, LOCK, with part 0's c* and unit check, before
any class is computed; `require_lock` refuses parts A to C until the lock is committed and matches):
- A response of a fly's KC: its trial mean above the criterion. The trial-noise unit: the trial mean over the pooled SD of
  the KC's single trials around their own stimulus means (every stimulus the recording presents, the solvent included,
  from the same trials as the reading; `pooled_sd`). Its criterion c* is the value at which the APL16 control flies'
  response rate (the mean over their 11 recordings of each recording's share of KC x odor entries above it over the 16
  odors, the registry's definition of 0.068) equals theirs at 1 z; fixed in part 0 (`find_cstar`) and applied unchanged.
- The unit check (part 0, `unit_check`): the APL16 control flies' 16 per-odor response rates (means over recordings) in
  the trial-noise unit at c* within UNIT_RMS_TOL (RMS) of those in baseline-noise units (data/yang_apl16_kc_noise_units/,
  the registry's threshold 1.9684). If it fails, 1 z carries and the trial-noise unit becomes partner (1).
- The evoked response rate: a recording's response rate to a stimulus minus its response rate to the solvent, in the
  same unit and trials (`evoked`).
- A set's value: the mean over flies (a fly is a date and a fly number; its recordings averaged), with a 90 % interval
  from 10,000 bootstrap draws of flies (`set_value`, E79's `fly_bootstrap`).
- The flies with APL output blocked (TNT): their value times and over exp(h), h the half-width on the log scale of set
  D1's control flies' interval at the same stimulus and reading ((ln hi - ln lo) / 2; `blocked_interval`); the three
  recordings' range beside.
- The family's response: two spikes or more; its evoked response rate is its response rate minus the blank's (0).
- The class of a cell (stimulus, condition, estimate, set; `cell_class`): 'short' if at least 6 of the 8 wirings lie
  below the flies' interval, 'over' if at least 6 above, 'inside' if at least 6 inside (ends included) or if neither of
  the first two holds and the eight wirings' median lies inside; otherwise 'not decided'. APL on against control flies,
  APL off against APL output blocked (TNT) (set D1 only).
- W per estimate, condition and set (`w_class`): 'short' if ma is short, or 2h and farn are both short; 'not short' if ma
  is inside or over and 2h and farn are not both short; otherwise 'not decided'.
- The carrying reading: evoked response rates, the trial-noise unit, all trials. Held (`held_w`): W the same under the
  carrying reading and under each partner (PARTNERS): (1) 1 z in place of the trial-noise unit; (2) the first three
  trials alone (the unit's SD from them too); (3) Yang's definition, above 1 z on at least half the trials; (4) one spike
  or more as the family's response; (5) each side's evoked response rate over its own panel's (the flies' over the
  APL16 flies' of the same group, the set's interval divided by the APL16 group's value; the family's over the mean of
  its 16 panel odors' rates in the same run and condition); (6) each control fly of the set left out in turn (for APL
  off, the left-out fly changes set D1's control interval and with it the blocked flies'), each of which must agree.
- Beside in every table, as its first row (reading 'carried_so_far'): the solvent not removed, 1 z, all trials, D1 and
  D2 pooled as recording means (E49's values) with an interval from bootstrap draws of their 12 control flies (a drawn
  fly brings its recordings, E104's `flies_rate_ratio`), the three blocked recordings' mean, kiwi-control's estimate.
- Part B's class (`e102_class`): at a step a stimulus is 'further' if the median wiring's distance from the flies' value
  on the log scale (|ln(median of the 8 wirings' evoked rates / the set's evoked value)|) is 1.25 times the unchanged
  fits' (E102's identity step) or more, 'nearer' if 0.8 times or less, 'level' otherwise; E102's own rule (a difference of
  ln 1.25, E102's `low_conc_class`) beside. Outcome (e) of E102 "stands like for like" if ma at pn_inhadd 1000 is further
  in both conditions under every estimate and against every set.

Parts (CLI `--part`):
- 0 (p0_identity.json, p0_counts.csv, p0_cstar.json, p0_unit_check.csv; a failure stops the experiment, outcome (h)): the
  reader reproduces E49's data_fractions.csv and E104's b_flies.csv (carrying, first three trials, solvent subtracted);
  the family's response rates from the archived spike counts reproduce E104's b_model.csv at m = 1 and m = 0, E102's
  b_low_conc.csv and its refits' d_low_conc_inhadd1000.csv with 0 differences; the sets' counts; c* and the unit check.
  `--write-lock` then writes the lock (only if every identity passed).
- A (a1_flies.csv, a1_solvent.csv, a2_sets.csv, a3_cells.csv, a3_w.csv, a3_held.csv, a4_beside.csv, a_family.csv,
  a_realism.csv, a_summary.json with the outcomes (a) to (f)).
- B (b_cells.csv, b_e102_class.csv, b_refits.csv, b_summary.json).
- C (c_ratios.csv, c_dstar.csv, c_summary.json).
- D (simulations; no class rests on these; `--sub`): D1 the input's size (FACTORS), D2 the family read as the
  recordings were made (E73's per-KC session reading), D3 the trial noise on the real wiring (C2's matched threshold
  jitter, E72's part F); `--sub score` writes the d*_ tables.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E105_weak_odor_rates --part 0
    ./venv/bin/python -u -m experiments.E105_weak_odor_rates --part 0 --write-lock    # then commit the lock
    ./venv/bin/python -u -m experiments.E105_weak_odor_rates --part A
    ./venv/bin/python -u -m experiments.E105_weak_odor_rates --part B
    ./venv/bin/python -u -m experiments.E105_weak_odor_rates --part C
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E105_weak_odor_rates --part D --sub D1

Outputs: results/E105_weak_odor_rates/ in the shared directory; part D's counts and responses under
outputs/E105_weak_odor_rates/D/, archived to data/model_outputs/E105_weak_odor_rates/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import socket
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E105_weak_odor_rates'
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's / E98's
S_WIRINGS = ('real', 'P1', 'P2', 'M1')                              # E102's refits

# ---------------------------------------------------------------------------
# the flies' sets (the entry's; part 0 checks the counts)

WEAK_DIR = ('data', 'yang_weak_odor_kc')
APL16_DIR = ('data', 'yang_apl16_kc_olfdb')
NOISE_UNITS_DIR = ('data', 'yang_apl16_kc_noise_units')
APL_RECORDINGS_CSV = ('data', 'yang_apl_recordings', 'recordings.csv')    # the raw-only marks (use_until_further_notice)
APL16_FILES = {'control': ('TNT3c', 'dataset_20_long.csv'), 'apl_output_blocked': ('TNT', 'dataset_19_long.csv')}
APL16_GENOTYPE = {'control': 'TNT3c', 'apl_output_blocked': 'TNT'}
SETS = ('D1', 'D2', 'N')
SET_RULES = {'D1': ('diag_mix', None, '2025-02-24'), 'D2': ('diag_mix', '2025-03-20', None),
             'N': ('natmix_dilute', '2026-01-16', '2026-02-10')}    # (panel, first date, last date), inclusive
SET_COUNTS = {('D1', 'control'): (9, 7), ('D1', 'apl_output_blocked'): (3, 2), ('D2', 'control'): (7, 5),
              ('N', 'control'): (9, 5), ('APL16', 'control'): (11, 8), ('APL16', 'apl_output_blocked'): (8, 7)}
GROUPS = ('control', 'apl_output_blocked')
CONDITION_GROUP = {'on': 'control', 'off': 'apl_output_blocked'}
CONDITION_RUN = {'on': 'on', 'off': 'acute_off'}
SOLVENT = 'pfo @ 0'
POOLED_SET = 'D1+D2'          # the reading carried so far

# ---------------------------------------------------------------------------
# the stimuli

ESTIMATES = ('kiwi-control', 'megamat', 'validation2')
CARRIED_ESTIMATE = 'kiwi-control'
WEAK_ODORS = ('ma', '2h', 'farn')
FLIES_STIM = {'ma': 'ma @ -7', '2h': '2h @ -7', 'farn': 'farn @ -3'}
MODEL_STIM = {'ma': ('ma', -7.0), '2h': ('2h', -6.0), 'farn': ('farn', -2.0)}   # the input's (E98's B1 diagnostics)
INPUT_ROLE = {'ma': 'exact input', '2h': 'input at 10^-6 for the flies\' 10^-7',
              'farn': 'input at 10^-2 for the flies\' 10^-3'}
LOWER_STIM = {'ma': 'ma @ -8', '2h': '2h @ -8', 'farn': 'farn @ -4'}         # set D1, beside (A1)
BESIDE_STIM = (('2h', '2h @ -5', '2h [kiwi-control/control] @ -5', 'N', 'exact input, one estimate'),
               ('ea', 'ea @ -7', 'ea [kiwi-control/diag] @ -8', 'N', 'input a log unit lower: only inside or over is safe'))
A1_STIMULI = ('ma @ -7', '2h @ -7', 'farn @ -3', 'ma @ -8', '2h @ -8', 'farn @ -4', '2h @ -5', 'ea @ -7', SOLVENT)
BLANK = 'blank @ 0'

# ---------------------------------------------------------------------------
# the definitions (hashed into the lock)

Z1 = 1.0                      # Yang's criterion, 1 z
FIRST_TRIALS = 3
YANG_TRIAL_SHARE = 0.5        # Yang's definition: above 1 z on at least half the trials
UNIT_RMS_TOL = 0.010          # part 0: the trial-noise unit against baseline-noise units (16 per-odor rates, RMS)
NOISE_UNITS_THRESHOLD = 1.9684    # the registry's noise_threshold (data/yang_apl16_kc_noise_units/README.md)
CSTAR_BRACKET = (0.5, 6.0)
CSTAR_TOL = 1e-6
N_BOOT = 10000
SEED = 0
CI = (5.0, 95.0)
MIN_WIRINGS = 6               # of the 8
SPIKES = 2                    # the family's response
SPIKES_PARTNER = 1
CELL_CLASSES = ('short', 'over', 'inside', 'not decided')
W_CLASSES = ('short', 'not short', 'not decided')
PARTNERS = ('unit', 'first_three_trials', 'yang_definition', 'one_spike', 'relative_to_panel', 'fly_left_out')
COUNTING_PARTNERS = ('unit', 'first_three_trials', 'yang_definition', 'one_spike', 'relative_to_panel')   # outcome (e)
SPLIT_SHARE = 0.1             # E98's split rule (A2)
E102_FURTHER, E102_NEARER = 1.25, 0.8   # part B: the median wiring's log distance over the unchanged fits'
E102_STEP_OUTCOME = 'inhadd1000'
E102_STEPS = ('identity', 'inhadd100', 'inhadd314', 'inhadd1000', 'inhadd3141', 'inhadd10000', 'gain_fixed')
E102_UNCHANGED = 'identity'
FACTORS = (0.5, 0.71, 1.41, 2.0, 2.83, 4.0)   # part D1
D3_RECORDINGS, D3_TRIALS = 10, 6

LOCK = Path(__file__).resolve().parent / 'E105_lock.json'
IDENTITY_TOL = 1e-12          # the flies' reader against the saved fractions (floating-point summation order)

E104_RESULTS = ('results', 'E104_kc_calyx_vs_apl')
E102_RESULTS = ('results', 'E102_antennal_lobe_normalization')
E102_ARCHIVE = ('data', 'model_outputs', 'E102_antennal_lobe_normalization')
E49_FRACTIONS = ('results', 'E49_low_conc_kc', 'data_fractions.csv')

QUALITY_CAVEAT = ("Yang's processed per-trial values (data/yang_weak_odor_kc/, data/yang_apl16_kc_olfdb/); the family's "
                  "saved runs; nothing re-tuned")
UNMEASURED = ("unmeasured: z-scored nuclear GCaMP6s against spike counts; a chronic block (TNT) against the model's acute "
              "APL off; two flies with APL output blocked; other flies and sessions for the ORN input than for the KCs; "
              "the weak-odor recordings in baseline-noise units (their traces are not on hand)")


def out_dir() -> Path:
    """results/E105_weak_odor_rates/ in the shared directory."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    """outputs/E105_weak_odor_rates/<sub>/ in this repository."""
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _json_default(x):
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def write_table(df: pd.DataFrame, path: Path) -> None:
    """A results table: two '#' lines (QUALITY_CAVEAT, UNMEASURED), then the csv (`read_table` reads it)."""
    with open(path, 'w') as f:
        f.write(f'# {QUALITY_CAVEAT}\n# {UNMEASURED}\n')
        df.to_csv(f, index=False)


def read_table(path: Path) -> pd.DataFrame:
    """A table written by `write_table` (or any csv with leading '#' lines)."""
    n = 0
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line.startswith('#'):
                break
            n += 1
    return pd.read_csv(path, skiprows=n, float_precision='round_trip')


def _commits() -> dict:
    return {k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}


def provenance(part: str, **extra) -> None:
    """provenance.json in the results folder: one entry per run of a part (commits, host, time, settings)."""
    path = out_dir() / 'provenance.json'
    d = _read_json(path) or dict(experiment=EXPERIMENT, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(part=part, commits=_commits(), host=_host(),
                                         time=time.strftime('%Y-%m-%d %H:%M'), **extra))
    _write_json(path, d)


# ---------------------------------------------------------------------------
# the reader

def check_path(path: Path) -> Path:
    """`path` resolved, if it lies under the shared directory; PermissionError otherwise."""
    p = Path(path).resolve()
    root = H.shared_dir().resolve()
    if p != root and root not in p.parents:
        raise PermissionError(f'{path}: outside the shared directory ({root}); the reader reads nothing else')
    return p


def raw_only_stems(path: Optional[Path] = None) -> set:
    """The recordings data/yang_apl_recordings/recordings.csv marks raw only (use_until_further_notice 'raw only ...');
    empty when the file is absent (none of the KC recordings read here is in it)."""
    p = H.shared_dir().joinpath(*APL_RECORDINGS_CSV) if path is None else check_path(path)
    if not Path(p).exists():
        return set()
    r = pd.read_csv(p, usecols=['stem', 'use_until_further_notice'])
    return set(r.loc[r['use_until_further_notice'].astype(str).str.startswith('raw only'), 'stem'].astype(str))


def check_rows(df: pd.DataFrame, allowed: Sequence[str], raw_only: Optional[set] = None) -> pd.DataFrame:
    """`df` unchanged if every row's recording is in `allowed` (the panel's processed list) and none is raw only (in
    `raw_only`, or a 'use_until_further_notice' column saying 'raw only'); PermissionError otherwise."""
    rec = df['recording'].astype(str)
    if 'use_until_further_notice' in df.columns:
        bad = df['use_until_further_notice'].astype(str).str.startswith('raw only')
        if bad.any():
            raise PermissionError(f'raw-only rows: {sorted(set(rec[bad]))[:5]} (not used until further notice)')
    out = sorted(set(rec) - set(map(str, allowed)))
    if out:
        raise PermissionError(f"recordings outside the panel's processed list: {out[:5]}")
    raw = sorted(set(rec) & set(raw_only or ()))
    if raw:
        raise PermissionError(f'raw-only recordings: {raw[:5]} (not used until further notice)')
    return df


def assign_set(dates: pd.Series, panel: str) -> pd.Series:
    """Each row's set (SETS) from its panel and date ('YYYY-MM-DD', SET_RULES; inclusive); a row in no set raises."""
    d = dates.astype(str)
    out = pd.Series(np.nan, index=dates.index, dtype=object)
    for s, (p, lo, hi) in SET_RULES.items():
        if p != panel:
            continue
        m = np.ones(len(d), bool)
        if lo is not None:
            m &= (d >= lo).to_numpy()
        if hi is not None:
            m &= (d <= hi).to_numpy()
        out[m] = s
    if out.isna().any():
        raise ValueError(f'{panel}: rows in no set: {sorted(set(d[out.isna()]))[:5]}')
    return out


_TRIALS: Dict[str, pd.DataFrame] = {}


def read_trials(panel: str, path: Optional[Path] = None, recordings: Optional[Path] = None) -> pd.DataFrame:
    """A panel's per-trial values ('diag_mix', 'natmix_dilute' or 'APL16'), through `check_path` and `check_rows`:
    set, recording, genotype_group, date, fly_num, fly ('<date>_fly<n>'), roi, odor, repeat, value (Yang's z). `path` /
    `recordings`: another long table / processed list (tests), under the shared directory too. Cached per process."""
    key = f'{panel}|{path}|{recordings}'
    if key in _TRIALS:
        return _TRIALS[key]
    raw = raw_only_stems()
    if panel in ('diag_mix', 'natmix_dilute'):
        base = H.shared_dir().joinpath(*WEAK_DIR, panel)
        p = check_path(base / f'{panel}_long.csv' if path is None else path)
        rp = check_path(base / 'computed' / 'recordings.csv' if recordings is None else recordings)
        allowed = pd.read_csv(rp)['recording'].astype(str)
        t = pd.read_csv(p, usecols=['date', 'fly_num', 'hemibrain', 'roi', 'odor', 'repeat', 'value', 'recording',
                                    'genotype_group'])
        t = check_rows(t, allowed, raw)
        t['set'] = assign_set(t['date'], panel).to_numpy()
    elif panel == 'APL16':
        base = H.shared_dir().joinpath(*APL16_DIR)
        rp = check_path(base / 'computed' / 'recordings.csv' if recordings is None else recordings)
        rr = pd.read_csv(rp)
        parts = []
        for g, (sub, fn) in APL16_FILES.items():
            p = check_path(base / sub / fn if path is None else Path(path) / sub / fn)
            x = pd.read_csv(p, usecols=['date', 'fly_num', 'hemibrain', 'roi', 'odor', 'repeat', 'value'])
            x['genotype_group'] = g
            x['recording'] = x['date'].astype(str) + '_fly' + x['fly_num'].astype(str) + '_' + x['hemibrain'].astype(str)
            r = rr[rr['genotype'] == APL16_GENOTYPE[g]]
            allowed = (r['date'].astype(str) + '_fly' + r['fly_num'].astype(str) + '_' + r['hemibrain'].astype(str))
            parts.append(check_rows(x, allowed, raw))
        t = pd.concat(parts, ignore_index=True)
        t['set'] = 'APL16'
    else:
        raise ValueError(f'unknown panel {panel!r}')
    t['fly'] = t['date'].astype(str) + '_fly' + t['fly_num'].astype(str)
    _TRIALS[key] = t
    return t


def read_all_trials() -> pd.DataFrame:
    """The three weak-odor sets and the APL16 recordings, one frame (a 'set' column keeps them apart)."""
    return pd.concat([read_trials('diag_mix'), read_trials('natmix_dilute'), read_trials('APL16')], ignore_index=True)


def set_counts(trials: pd.DataFrame) -> Dict[Tuple[str, str], Tuple[int, int]]:
    """{(set, group): (recordings, flies)}."""
    g = trials.drop_duplicates(['set', 'genotype_group', 'recording'])
    return {(s, grp): (int(len(x)), int(x['fly'].nunique())) for (s, grp), x in g.groupby(['set', 'genotype_group'])}


# ---------------------------------------------------------------------------
# the units and the rates

KC_KEYS = ['set', 'recording', 'genotype_group', 'fly', 'roi']


def pooled_sd(trials: pd.DataFrame) -> pd.Series:
    """Per (set, recording, group, fly, roi), the pooled SD of its single trials around their own stimulus means: the
    sum over stimuli of the squared deviations over the trials minus the number of stimuli (stimuli with one trial add
    no degree of freedom)."""
    g = trials.groupby(KC_KEYS + ['odor'], sort=False)['value']
    n = g.size()
    ss = (g.var(ddof=1) * (n - 1)).fillna(0.0)
    k = pd.DataFrame(dict(ss=ss, n=n, one=1)).groupby(level=KC_KEYS, sort=False).sum()
    dof = k['n'] - k['one']
    return np.sqrt(k['ss'] / dof.where(dof > 0)).rename('pooled_sd')


def kc_table(trials: pd.DataFrame, first_n: Optional[int] = None) -> pd.DataFrame:
    """Per recording, KC and stimulus (from the same trials: all, or the first `first_n`): the trial mean in z ('tm'),
    the number of trials ('n'), the share of trials above 1 z ('share_gt1'), the KC's pooled SD ('pooled_sd') and the
    trial mean in the trial-noise unit ('tn' = tm / pooled_sd)."""
    t = trials if first_n is None else trials[trials['repeat'] < first_n]
    keys = KC_KEYS + ['odor']
    g = t.groupby(keys, sort=False)['value']
    out = pd.DataFrame(dict(tm=g.mean(), n=g.size(),
                            share_gt1=(t['value'] > Z1).groupby([t[k] for k in keys], sort=False).mean()))
    sd = pooled_sd(t)
    out = out.join(sd, on=KC_KEYS)
    out['tn'] = out['tm'] / out['pooled_sd']
    return out.reset_index()


UNITS = ('trial_noise', 'z1', 'yang')


def responses(kc: pd.DataFrame, unit: str, cstar: Optional[float] = None) -> pd.Series:
    """Per KC row of `kc`, whether it responds under `unit`: 'trial_noise' (tn above `cstar`), 'z1' (trial mean above 1
    z) or 'yang' (above 1 z on at least half the trials)."""
    if unit == 'trial_noise':
        if cstar is None:
            raise ValueError('the trial-noise unit needs c*')
        return kc['tn'] > float(cstar)
    if unit == 'z1':
        return kc['tm'] > Z1
    if unit == 'yang':
        return kc['share_gt1'] >= YANG_TRIAL_SHARE
    raise ValueError(f'unknown unit {unit!r}')


def response_rates(kc: pd.DataFrame, unit: str, cstar: Optional[float] = None) -> pd.DataFrame:
    """Per recording and stimulus, the response rate under `unit` (`responses`): set, recording, genotype_group, fly,
    odor, rate, n_kcs."""
    r = responses(kc, unit, cstar)
    keys = ['set', 'recording', 'genotype_group', 'fly', 'odor']
    g = r.groupby([kc[k] for k in keys], sort=False)
    return pd.DataFrame(dict(rate=g.mean(), n_kcs=g.size())).reset_index()


def evoked(rates: pd.DataFrame) -> pd.DataFrame:
    """`rates` with 'solvent' (the recording's rate to the solvent, same unit and trials) and 'evoked' (rate minus it)."""
    s = rates[rates['odor'] == SOLVENT].set_index('recording')['rate']
    x = rates.copy()
    x['solvent'] = x['recording'].map(s).astype(float)
    if x['solvent'].isna().any():
        raise RuntimeError(f"recordings without the solvent: {sorted(set(x.loc[x['solvent'].isna(), 'recording']))}")
    x['evoked'] = x['rate'] - x['solvent']
    return x


def panel_rate(rates: pd.DataFrame, odors: Sequence[str]) -> pd.DataFrame:
    """Per recording the mean over `odors` of its per-odor response rates (each recording has the same KCs for every
    odor, so this is its share of KC x odor entries above the criterion), and of its evoked rates when present."""
    x = rates[rates['odor'].isin(list(odors))]
    n = x.groupby('recording')['odor'].nunique()
    if (n != len(odors)).any():
        raise RuntimeError(f'recordings without all {len(odors)} odors: {list(n[n != len(odors)].index)[:5]}')
    cols = [c for c in ('rate', 'evoked') if c in x.columns]
    return x.groupby(['set', 'recording', 'genotype_group', 'fly'], sort=False)[cols].mean().reset_index()


def apl16_odor_labels(trials: pd.DataFrame) -> List[str]:
    """The APL16 recordings' 16 odor strings ('<abbrev> @ <conc>'), the solvent left out, in H.APL16_ODORS's order."""
    by = {o.split(' @ ')[0]: o for o in trials['odor'].unique()}
    return [by[o] for o in H.APL16_ODORS]


def _mean_panel_rate(tn: np.ndarray, rec_idx: np.ndarray, n_rec: int, c: float) -> float:
    hits = np.bincount(rec_idx, weights=(tn > c).astype(float), minlength=n_rec)
    tot = np.bincount(rec_idx, minlength=n_rec).astype(float)
    return float((hits / tot).mean())


def find_cstar(kc: pd.DataFrame, target: float, odors: Sequence[str], bracket: Tuple[float, float] = CSTAR_BRACKET,
               tol: float = CSTAR_TOL) -> dict:
    """c*: the smallest trial-noise criterion (to `tol`, by bisection over `bracket`) at which the mean over recordings of
    each recording's share of KC x odor entries (over `odors`) above it is at most `target`; with the rate there and
    just below. `kc`: one group's `kc_table` rows."""
    x = kc[kc['odor'].isin(list(odors))]
    codes, recs = pd.factorize(x['recording'])
    tn = x['tn'].to_numpy(float)
    f = lambda c: _mean_panel_rate(tn, codes, len(recs), c)
    lo, hi = bracket
    if not (f(lo) > target >= f(hi)):
        raise RuntimeError(f'c* is not inside {bracket}: rates {f(lo)}, {f(hi)} against {target}')
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) > target:
            lo = mid
        else:
            hi = mid
    return dict(cstar=float(hi), rate_at_cstar=f(hi), rate_below=f(lo), target=float(target), tol=tol,
                n_recordings=int(len(recs)))


def unit_check(kc_ctrl: pd.DataFrame, cstar: float, odors: Sequence[str], noise_path: Optional[Path] = None,
               group: str = 'control') -> dict:
    """The APL16 control flies' 16 per-odor response rates (means over their recordings) in the trial-noise unit at c*
    against those in baseline-noise units (data/yang_apl16_kc_noise_units/response_fraction_per_recording.csv,
    frac_noise: above the registry's 1.9684), and at 1 z beside: the table, the RMS and whether it is within
    UNIT_RMS_TOL. `group` 'apl_output_blocked' (with that group's rows as `kc_ctrl`): the same for the flies with APL
    output blocked (TNT), beside (part A; the entry's check is the control flies')."""
    rr = response_rates(kc_ctrl[kc_ctrl['odor'].isin(list(odors))], 'trial_noise', cstar)
    rz = response_rates(kc_ctrl[kc_ctrl['odor'].isin(list(odors))], 'z1')
    p = H.shared_dir().joinpath(*NOISE_UNITS_DIR, 'response_fraction_per_recording.csv') if noise_path is None \
        else check_path(noise_path)
    nu = pd.read_csv(p)
    nu = nu[nu['genotype'] == APL16_GENOTYPE[group]]
    rows = []
    for o in odors:
        a = o.split(' @ ')[0]
        rows.append(dict(odor=o, trial_noise=float(rr.loc[rr['odor'] == o, 'rate'].mean()),
                         baseline_noise=float(nu.loc[nu['odor'] == a, 'frac_noise'].mean()),
                         z1=float(rz.loc[rz['odor'] == o, 'rate'].mean()),
                         n_recordings_trial_noise=int((rr['odor'] == o).sum()),
                         n_recordings_baseline_noise=int((nu['odor'] == a).sum())))
    t = pd.DataFrame(rows)
    t['diff'] = t['trial_noise'] - t['baseline_noise']
    rms = float(np.sqrt((t['diff'] ** 2).mean()))
    rms_z = float(np.sqrt(((t['z1'] - t['baseline_noise']) ** 2).mean()))
    return dict(table=t, rms=rms, rms_z1_beside=rms_z, tol=UNIT_RMS_TOL, passed=bool(rms <= UNIT_RMS_TOL),
                noise_file=str(p), noise_file_sha256=_sha256(Path(p)))


# ---------------------------------------------------------------------------
# the flies' values

def fly_values(rec_values: pd.DataFrame, col: str) -> pd.Series:
    """Per fly, the mean of its recordings' `col`."""
    return rec_values.groupby('fly', sort=True)[col].mean()


def _ci(x: np.ndarray) -> Tuple[float, float]:
    return float(np.percentile(x, CI[0])), float(np.percentile(x, CI[1]))


def set_value(per_fly: pd.Series, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The mean over flies with its 90 % interval from `n_boot` bootstrap draws of flies (E79's `fly_bootstrap`)."""
    x = per_fly.dropna().to_numpy(float)
    if len(x) == 0:
        return dict(value=np.nan, lo90=np.nan, hi90=np.nan, n_flies=0)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    lo, hi = _ci(x[idx].mean(axis=1))
    return dict(value=float(x.mean()), lo90=lo, hi90=hi, n_flies=int(len(x)))


def recording_mean_value(rec: pd.DataFrame, col: str, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The reading carried so far: the mean over recordings (E49's), its 90 % interval from bootstrap draws of flies (a
    drawn fly brings its recordings; E104's `flies_rate_ratio`)."""
    flies = sorted(rec['fly'].unique())
    by_fly = [rec.loc[rec['fly'] == f, col].to_numpy(float) for f in flies]
    sums = np.array([v.sum() for v in by_fly])
    ns = np.array([len(v) for v in by_fly])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(flies), size=(n_boot, len(flies)))
    lo, hi = _ci(sums[idx].sum(axis=1) / ns[idx].sum(axis=1))
    return dict(value=float(rec[col].mean()), lo90=lo, hi90=hi, n_flies=len(flies))


def blocked_interval(value: float, ctrl_value: float, ctrl_lo: float, ctrl_hi: float) -> Tuple[float, float]:
    """The blocked flies' interval: `value` over and times exp(h), h the half-width of set D1's control flies' interval
    on the log scale ((ln hi - ln lo) / 2). NaN where the control interval's ends or the value are not positive."""
    if not (np.isfinite(value) and value > 0 and ctrl_lo > 0 and ctrl_hi > 0):
        return np.nan, np.nan
    h = 0.5 * (np.log(ctrl_hi) - np.log(ctrl_lo))
    return float(value * np.exp(-h)), float(value * np.exp(h))


# ---------------------------------------------------------------------------
# the classes

def cell_class(model: Sequence[float], lo: float, hi: float, min_wirings: int = MIN_WIRINGS) -> str:
    """A cell's class (CELL_CLASSES) from the wirings' values and the flies' interval [lo, hi] (ends inside): 'short' if
    at least `min_wirings` lie below, 'over' if at least that many above, 'inside' if that many inside or (neither of
    the two) the median lies inside; otherwise 'not decided'. 'not decided' when the interval is not defined."""
    v = np.asarray(model, float)
    if not (np.isfinite(lo) and np.isfinite(hi)) or len(v) == 0 or not np.isfinite(v).all():
        return 'not decided'
    below, above = int((v < lo).sum()), int((v > hi).sum())
    inside = len(v) - below - above
    if below >= min_wirings:
        return 'short'
    if above >= min_wirings:
        return 'over'
    if inside >= min_wirings:
        return 'inside'
    med = float(np.median(v))
    if lo <= med <= hi:
        return 'inside'
    return 'not decided'


def w_class(classes: Mapping[str, str]) -> str:
    """W (W_CLASSES) from the cells of ma, 2h and farn: 'short' if ma is short or 2h and farn both are; 'not short' if ma
    is inside or over and 2h and farn are not both short; otherwise 'not decided' (farn absent in set N: 'both short'
    cannot hold)."""
    ma = classes.get('ma')
    both = classes.get('2h') == 'short' and classes.get('farn') == 'short'
    if ma == 'short' or both:
        return 'short'
    if ma in ('inside', 'over'):
        return 'not short'
    return 'not decided'


def held_w(carrying: str, partners: Mapping[str, object]) -> dict:
    """Whether W is held: the same under the carrying reading and every partner (a partner's value may be a list, one W
    per left-out fly, each of which must agree). Returns 'held' (the carrying W, or 'not held'), the partners that
    differ, and those that differ with both W decided (short against not short)."""
    differ, flips = [], []
    for name, v in partners.items():
        vals = list(v) if isinstance(v, (list, tuple)) else [v]
        if any(x != carrying for x in vals):
            differ.append(name)
        if carrying in ('short', 'not short') and any(x in ('short', 'not short') and x != carrying for x in vals):
            flips.append(name)
    return dict(held=carrying if not differ else 'not held', differ=differ, flips=flips)


def _turns(values: Sequence[str]) -> Tuple[bool, bool]:
    """(strict, loose): both 'short' and 'not short' among `values`; any two different."""
    v = list(values)
    return ('short' in v and 'not short' in v), (len(set(v)) > 1)


def outcomes(held: pd.DataFrame) -> dict:
    """Outcomes (a) to (f) from part A's held table (one row per estimate, condition and set: 'carrying', 'held',
    'differ', 'flips'). "Turns on" is read strictly (a decided W that changes, short against not short) and loosely (any
    difference) beside; the strict reading gives the outcome."""
    h = held.copy()
    combos = len(h)
    res = {}
    res['a'] = dict(holds=bool(combos and (h['held'] == 'short').all()), n=combos,
                    n_held_short=int((h['held'] == 'short').sum()))
    res['b'] = dict(holds=bool(combos and (h['held'] == 'not short').all()), n=combos,
                    n_held_not_short=int((h['held'] == 'not short').sum()))
    c_strict, c_loose = [], []
    for (cond, s), g in h.groupby(['condition', 'set'], sort=False):
        st, lo = _turns(g['carrying'])
        if st:
            c_strict.append(f'{cond}/{s}')
        if lo:
            c_loose.append(f'{cond}/{s}')
    res['c'] = dict(holds=bool(c_strict), where=c_strict, where_loose=c_loose)
    d_strict, d_loose = [], []
    for (cond, e), g in h.groupby(['condition', 'estimate'], sort=False):
        if g['set'].nunique() < 2:
            continue
        st, lo = _turns(g['carrying'])
        if st:
            d_strict.append(f'{cond}/{e}')
        if lo:
            d_loose.append(f'{cond}/{e}')
    res['d'] = dict(holds=bool(d_strict), where=d_strict, where_loose=d_loose)
    e_strict, e_loose = [], []
    for r in h.itertuples(index=False):
        fl = [p for p in json.loads(r.flips) if p in COUNTING_PARTNERS]
        df_ = [p for p in json.loads(r.differ) if p in COUNTING_PARTNERS]
        if fl:
            e_strict.append(f'{r.condition}/{r.estimate}/{r.set}: {", ".join(fl)}')
        if df_:
            e_loose.append(f'{r.condition}/{r.estimate}/{r.set}: {", ".join(df_)}')
    res['e'] = dict(holds=bool(e_strict), where=e_strict, where_loose=e_loose)
    f_strict, f_loose = [], []
    d1 = h[h['set'] == 'D1']
    for e, g in d1.groupby('estimate', sort=False):
        if set(g['condition']) != {'on', 'off'}:
            continue
        st, lo = _turns(g['carrying'])
        if st:
            f_strict.append(e)
        if lo:
            f_loose.append(e)
    res['f'] = dict(holds=bool(f_strict), where=f_strict, where_loose=f_loose)
    fly = [f'{r.condition}/{r.estimate}/{r.set}' for r in h.itertuples(index=False)
           if 'fly_left_out' in json.loads(r.differ)]
    res['beside_fly_left_out_differs'] = fly
    return res


def log_distance(model: float, flies: float) -> float:
    """|ln(model / flies)|; inf for a model rate at or below 0; NaN for a flies' value at or below 0."""
    if not (np.isfinite(flies) and flies > 0):
        return np.nan
    if not (np.isfinite(model) and model > 0):
        return float('inf')
    return float(abs(np.log(model / flies)))


def e102_class(model_med: float, unchanged_med: float, flies: float, further: float = E102_FURTHER,
               nearer: float = E102_NEARER) -> str:
    """Part B's class: 'further' if the step's log distance from the flies (`log_distance` of the wirings' median) is at
    least `further` times the unchanged fits', 'nearer' if at most `nearer` times, 'level' otherwise ('not read' where
    the flies' value is not positive)."""
    d, d0 = log_distance(model_med, flies), log_distance(unchanged_med, flies)
    if not (np.isfinite(d) or np.isinf(d)) or np.isnan(d) or np.isnan(d0):
        return 'not read'
    if np.isinf(d) and np.isinf(d0):
        return 'level'
    if np.isinf(d):
        return 'further'
    if np.isinf(d0):
        return 'nearer'
    if d0 == 0:
        return 'level' if d == 0 else 'further'
    q = d / d0
    if q >= further - 1e-12:
        return 'further'
    if q <= nearer + 1e-12:
        return 'nearer'
    return 'level'


def e102_own_class(model_med: float, unchanged_med: float, flies: float) -> str:
    """E102's own rule beside (`low_conc_class`: a difference of the log distances of ln 1.25)."""
    if not (np.isfinite(flies) and flies > 0):
        return 'not read'
    from experiments import E102_antennal_lobe_normalization as E102
    return E102.low_conc_class(model_med, flies, unchanged_med)


# ---------------------------------------------------------------------------
# the lock

def definitions() -> dict:
    """The definitions, lists and thresholds (module constants only) as one JSON-able dict."""
    return dict(
        experiment=EXPERIMENT,
        sets=dict(rules={k: list(v) for k, v in SET_RULES.items()},
                  counts={f'{s}|{g}': list(v) for (s, g), v in SET_COUNTS.items()}, never_pooled=list(SETS),
                  fly='a date and a fly number; its recordings averaged',
                  apl16=dict(files={k: list(v) for k, v in APL16_FILES.items()}),
                  reader="refuses paths outside the shared directory, recordings outside the panel's computed/"
                         "recordings.csv and recordings data/yang_apl_recordings/recordings.csv marks raw only"),
        stimuli=dict(estimates=list(ESTIMATES), odors=list(WEAK_ODORS), flies=FLIES_STIM,
                     model={k: list(v) for k, v in MODEL_STIM.items()}, input_role=INPUT_ROLE, lower=LOWER_STIM,
                     beside=[list(x) for x in BESIDE_STIM], a1=list(A1_STIMULI), solvent=SOLVENT, blank=BLANK,
                     model_label="E98's session_label(odor, conc, estimate, glomeruli_diagnostics)"),
        response=dict(
            unit="trial-noise unit: the trial mean over the pooled SD of the KC's single trials around their own "
                 "stimulus means (every stimulus of the recording, the solvent included; sum of squares over the "
                 "trials minus the stimuli), from the reading's own trials",
            cstar="the smallest criterion (bisection to cstar_tol over cstar_bracket) at which the APL16 control flies' "
                  "mean over recordings of each recording's share of KC x odor entries above it over the 16 odors is "
                  "at most their value at 1 z",
            cstar_bracket=list(CSTAR_BRACKET), cstar_tol=CSTAR_TOL, z1=Z1,
            unit_check=dict(tol_rms=UNIT_RMS_TOL, baseline_noise_threshold=NOISE_UNITS_THRESHOLD,
                            rule="the 16 per-odor rates (means over the 11 control recordings) at c* against "
                                 "frac_noise of data/yang_apl16_kc_noise_units/; failing it, 1 z carries and the "
                                 "trial-noise unit is partner (1)"),
            evoked='a recording\'s rate minus its solvent\'s, same unit and trials',
            first_trials=FIRST_TRIALS, yang_trial_share=YANG_TRIAL_SHARE,
            family=f'{SPIKES} spikes or more; evoked = rate minus the blank\'s', family_partner=SPIKES_PARTNER),
        values=dict(set='mean over flies; 90 % interval from n_boot bootstrap draws of flies', n_boot=N_BOOT, seed=SEED,
                    ci=list(CI), blocked="value over and times exp(h), h = (ln hi - ln lo) / 2 of set D1's control "
                                         "flies' interval at the same stimulus and reading; the recordings' range beside",
                    carried_so_far="the solvent not removed, 1 z, all trials, D1 and D2 pooled as recording means "
                                   "(interval: bootstrap of their flies, a fly bringing its recordings), the three "
                                   "blocked recordings' mean, kiwi-control's estimate"),
        classes=dict(min_wirings=MIN_WIRINGS, labels=list(CELL_CLASSES),
                     cell="short if >= min_wirings of 8 below the interval; over if >= min_wirings above; inside if "
                          ">= min_wirings inside (ends included) or, neither, the median inside; else not decided; "
                          "APL on against control flies, APL off against APL output blocked (TNT), set D1 only",
                     w="short if ma short or 2h and farn both short; not short if ma inside or over and 2h and farn not "
                       "both short; else not decided", w_labels=list(W_CLASSES)),
        held=dict(carrying='evoked, the trial-noise unit, all trials', partners=list(PARTNERS),
                  rule="W the same under the carrying reading and every partner; partner fly_left_out: one reading "
                       "per control fly of the set left out (for APL off, set D1's), each must agree",
                  relative_to_panel="the set's evoked value and interval over the APL16 flies' evoked panel value of "
                                    "the same group (mean over flies of each recording's mean over the 16 odors of "
                                    "its evoked rates); the family's evoked rate over its 16 panel odors' mean rate "
                                    "in the same run and condition"),
        outcomes=dict(turns="strict: both short and not short among the values (decides the outcome); loose: any "
                            "difference (beside)", counting_partners=list(COUNTING_PARTNERS)),
        split_rule=dict(share=SPLIT_SHARE, values='per-fly evoked rates (arithmetic), E103\'s split_share',
                        ratio='D1 over the other set, joint bootstrap of both sets\' flies (E104\'s q_of_fall, '
                              'arithmetic)'),
        part_b=dict(steps=list(E102_STEPS), unchanged=E102_UNCHANGED, further=E102_FURTHER, nearer=E102_NEARER,
                    outcome_step=E102_STEP_OUTCOME,
                    rule="|ln(median of the wirings' evoked rates / the set's evoked value)| over the unchanged "
                         "fits' (E102's identity step): further at >= 1.25, nearer at <= 0.8; E102's own rule beside",
                    stands="ma at pn_inhadd 1000 further in both conditions under every estimate and against every "
                           "set", refits=list(S_WIRINGS)),
        part_c=dict(ratio="APL output blocked (TNT) over set D1's control flies alone (E104's flies_rate_ratio), "
                          "evoked and as recorded", d_star="E104's d_star on its saved b_model.csv",
                    y2="E104's y_class, y_over_wirings, y_over_estimates and y_verdict"),
        part_d=dict(factors=list(FACTORS), d3_recordings=D3_RECORDINGS, d3_trials=D3_TRIALS),
        identity_tol=IDENTITY_TOL)


def _canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), default=_json_default)


def definitions_sha(defs: dict) -> str:
    """sha256 of the canonical JSON (sorted keys, no spaces)."""
    return hashlib.sha256(_canonical(defs).encode()).hexdigest()


def _part0_locked_values(p0: dict) -> dict:
    return dict(cstar=float(p0['cstar']['cstar']), carrying_unit=p0['carrying_unit'], partner_unit=p0['partner_unit'],
                unit_check_rms=float(p0['unit_check']['rms']), unit_check_passed=bool(p0['unit_check']['passed']),
                cstar_target=float(p0['cstar']['target']))


def write_lock(p0: dict, path: Optional[Path] = None) -> dict:
    """The lock: the definitions, part 0's c* and unit check (which unit carries), their sha256s, host, time and
    commits. Refuses to write unless part 0 passed, and to overwrite a lock with other definitions or values."""
    path = LOCK if path is None else Path(path)
    if not p0.get('passed'):
        raise RuntimeError('part 0 did not pass: no lock')
    defs = definitions()
    vals = _part0_locked_values(p0)
    sha, vsha = definitions_sha(defs), definitions_sha(vals)
    old = _read_json(path)
    if old is not None and (old.get('sha256') != sha or old.get('values_sha256') != vsha):
        raise RuntimeError(f'{path.name} exists with other definitions or values; not overwritten')
    lock = dict(experiment=EXPERIMENT, sha256=sha, values_sha256=vsha, definitions=defs, values=vals,
                host=_host(), time=time.strftime('%Y-%m-%d %H:%M'), commits=_commits(),
                note="Written by part 0 before any class was computed. `require_lock` refuses parts A to C until this "
                     "file is committed, its sha256s are its contents' and its definitions equal the module's.")
    _write_json(path, lock)
    return lock


def _committed(path: Path) -> bool:
    repo = Path(path).resolve().parent
    try:
        r = subprocess.run(['git', '-C', str(repo), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(repo), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def locked(path: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether the lock exists, its sha256s are its contents', its definitions equal `definitions()` and (`check_git`)
    it is committed."""
    path = LOCK if path is None else Path(path)
    if not path.exists():
        return False
    try:
        d = json.loads(path.read_text())
    except ValueError:
        return False
    defs, vals = d.get('definitions'), d.get('values')
    if not isinstance(defs, dict) or not isinstance(vals, dict):
        return False
    if d.get('sha256') != definitions_sha(defs) or d.get('values_sha256') != definitions_sha(vals):
        return False
    if _canonical(defs) != _canonical(definitions()):
        return False
    return _committed(path) if check_git else True


def require_lock(path: Optional[Path] = None, check_git: bool = True) -> dict:
    """The lock's values (c*, the carrying and partner units); PermissionError unless `locked`."""
    path = LOCK if path is None else Path(path)
    if not locked(path, check_git):
        raise PermissionError(f'{path.name}: no class is computed before the lock is committed and matches the '
                              'definitions')
    return json.loads(path.read_text())['values']


# ---------------------------------------------------------------------------
# the family's saved runs

def _e104():
    from experiments import E104_kc_calyx_vs_apl as E104
    return E104


def model_label(odor: str, estimate: str) -> str:
    """The B1 label of a weak odor's input under an estimate's diagnostics."""
    from experiments import E98_apl_response_strength as E98
    o, c = MODEL_STIM[odor]
    return E98.session_label(o, c, estimate, E98.DIAG)


def run_rates(counts: pd.DataFrame, labels: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Per stimulus of a run's spike counts (KCs x stimuli): the response rate at two spikes or more ('rate'), at one or
    more ('rate_1spike'), the largest count."""
    c = counts.astype(float)
    if labels is not None:
        c = c[list(labels)]
    return pd.DataFrame(dict(rate=(c >= SPIKES).mean(), rate_1spike=(c >= SPIKES_PARTNER).mean(), max_count=c.max()))


def family_table(runs: Mapping[str, Mapping[str, dict]], source: str = 'saved') -> pd.DataFrame:
    """One row per wiring, condition ('on', 'off') and stimulus of `runs` ({wiring: {'on': {...}, 'acute_off': {...}}},
    each with 'spike_counts' and 'responses'): rate, rate_1spike, max_count, the blank's rates, evoked rates (rate minus
    the blank's), the panel's mean rates (the 16 APL16 odors), and whether the saved responses are the counts at two
    spikes (responses_equal)."""
    rows = []
    for w, run in runs.items():
        for cond, ct in CONDITION_RUN.items():
            o = run[ct]
            s = o['spike_counts']
            rt = run_rates(s)
            resp_eq = bool(np.array_equal((s.astype(float).to_numpy() >= SPIKES), o['responses'].astype(float).to_numpy() > 0))
            bl = rt.loc[BLANK] if BLANK in rt.index else pd.Series(dict(rate=np.nan, rate_1spike=np.nan))
            p16 = [x for x in H.APL16_ODORS if x in rt.index]
            pr, pr1 = (float(rt.loc[p16, 'rate'].mean()), float(rt.loc[p16, 'rate_1spike'].mean())) if len(p16) == 16 \
                else (np.nan, np.nan)
            for lab, r in rt.iterrows():
                rows.append(dict(source=source, wiring=w, condition=cond, stimulus=str(lab), rate=float(r['rate']),
                                 rate_1spike=float(r['rate_1spike']), max_count=float(r['max_count']),
                                 blank_rate=float(bl['rate']), blank_rate_1spike=float(bl['rate_1spike']),
                                 panel_rate=pr, panel_rate_1spike=pr1, n_kcs=int(len(s)), responses_equal=resp_eq,
                                 dir=str(o.get('dir', ''))))
    t = pd.DataFrame(rows)
    t['evoked'] = t['rate'] - t['blank_rate']
    t['evoked_1spike'] = t['rate_1spike'] - t['blank_rate_1spike']
    return t


_FAMILY: Dict[str, pd.DataFrame] = {}


def saved_family(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """`family_table` of the eight wirings' saved runs (E104's `saved_run`: E98's B1 for the real wiring, E101's archived
    A3 for the others). Cached per process."""
    key = ','.join(wirings)
    if key not in _FAMILY:
        E104 = _e104()
        _FAMILY[key] = family_table({w: E104.saved_run(w) for w in wirings})
    return _FAMILY[key]


def e102_run(step: str, wiring: str, refit: bool = False) -> Optional[Dict[str, dict]]:
    """An E102 extended run from its archive: part B's step folder (B/<step>/B_<step>_<wiring>_...), or with `refit`
    part D's refit at the level (D/runs/B/<step>/D_<step>_<wiring>_...)."""
    from experiments import E98_apl_response_strength as E98
    arc = H.shared_dir().joinpath(*E102_ARCHIVE)
    d = arc / 'D' / 'runs' / 'B' / step if refit else arc / 'B' / step
    label = f"{'D' if refit else 'B'}_{step}_{wiring}"
    return E98.load_b1_run(label, roots=[d])


def e102_family(steps: Sequence[str] = E102_STEPS, wirings: Sequence[str] = ALL_WIRINGS,
                refit_step: str = E102_STEP_OUTCOME) -> pd.DataFrame:
    """`family_table` of E102's ladder (every step, eight wirings) and its refits (S's four wirings at `refit_step`), with
    'step' and 'set' ('ladder' or 'refit')."""
    parts = []
    for st in steps:
        runs = {w: r for w in wirings if (r := e102_run(st, w)) is not None}
        if len(runs) != len(wirings):
            raise RuntimeError(f'E102 step {st}: {sorted(set(wirings) - set(runs))} missing in the archive')
        parts.append(family_table(runs, source=f'E102:B/{st}').assign(step=st, kind='ladder'))
    runs = {w: r for w in S_WIRINGS if (r := e102_run(refit_step, w, refit=True)) is not None}
    if len(runs) != len(S_WIRINGS):
        raise RuntimeError(f'E102 refits at {refit_step}: {sorted(set(S_WIRINGS) - set(runs))} missing')
    parts.append(family_table(runs, source=f'E102:D/runs/B/{refit_step}').assign(step=refit_step, kind='refit'))
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# part 0

def _rates_cache(trials_weak: pd.DataFrame, first_n: Optional[int]) -> pd.DataFrame:
    return kc_table(trials_weak, first_n)


def e49_identity(kc_all: pd.DataFrame) -> Tuple[dict, pd.DataFrame]:
    """The reader against E49's data_fractions.csv: per panel, group and odor the number of recordings, the recording
    mean of the response rates at 1 z ('frac'), its SEM and, for mixtures, the mean of the union of the components'
    responders ('union')."""
    ref = pd.read_csv(H.shared_dir().joinpath(*E49_FRACTIONS))
    x = kc_all[kc_all['set'].isin(SETS)].copy()
    x['panel'] = np.where(x['set'] == 'N', 'natmix_dilute', 'diag_mix')
    x['resp'] = responses(x, 'z1')
    rr = x.groupby(['panel', 'genotype_group', 'recording', 'odor'], sort=False)['resp'].mean().rename('rate').reset_index()
    resp = x.set_index(['panel', 'genotype_group', 'recording', 'roi', 'odor'])['resp']
    rows = []
    for r in ref.itertuples(index=False):
        m = rr[(rr['panel'] == r.panel) & (rr['genotype_group'] == r.group) & (rr['odor'] == r.odor)]
        n = int(len(m))
        frac = float(m['rate'].mean()) if n else np.nan
        sem = float(m['rate'].std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan
        union = np.nan
        if ' + ' in r.odor and n:
            comps = [c.strip() for c in r.odor.split(' + ')]
            u = []
            for rec in m['recording']:
                parts = [resp.xs((r.panel, r.group, rec), level=[0, 1, 2]).xs(c, level='odor') for c in comps]
                u.append(float((parts[0] | parts[1]).mean()))
            union = float(np.mean(u))
        rows.append(dict(panel=r.panel, group=r.group, odor=r.odor, n_ref=int(r.n_recordings), n=n, frac_ref=r.frac,
                         frac=frac, sem_ref=r.sem, sem=sem, union_ref=r.union, union=union))
    t = pd.DataFrame(rows)
    d = {c: float(np.nanmax(np.abs(t[c] - t[f'{c}_ref']))) if t[[c, f'{c}_ref']].notna().all(axis=1).any() else 0.0
         for c in ('frac', 'sem', 'union')}
    nan_mismatch = int(((t['union'].isna()) != (t['union_ref'].isna())).sum())
    passed = bool((t['n'] == t['n_ref']).all() and all(v <= IDENTITY_TOL for v in d.values()) and nan_mismatch == 0)
    return dict(passed=passed, max_abs_diff=d, n_rows=int(len(t)), n_mismatch=int((t['n'] != t['n_ref']).sum()),
                union_nan_mismatch=nan_mismatch), t


def _e104_rates(kc: pd.DataFrame) -> pd.DataFrame:
    """The reader's rates at 1 z for diag_mix in E104's `weak_kc_rates` layout (recording, genotype_group, date,
    fly_num, odor, rate, fly)."""
    x = kc[kc['set'].isin(('D1', 'D2'))]
    r = response_rates(x, 'z1')
    r['fly'] = r['fly'].astype(str)
    # E104's groupby sorts by recording (then odor): the recording means are sums in that order
    return r.sort_values(['recording', 'odor'], kind='stable').reset_index(drop=True)


def e104_identity(kc_all: pd.DataFrame, kc_first3: pd.DataFrame) -> Tuple[dict, pd.DataFrame]:
    """The reader against E104's b_flies.csv (carrying, first_three_trials, solvent_subtracted): E104's
    `flies_rate_ratio` on the reader's rates reproduces every value exactly."""
    E104 = _e104()
    ref = E104.read_table(H.shared_dir().joinpath(*E104_RESULTS, 'b_flies.csv'))
    all_t, first = _e104_rates(kc_all), _e104_rates(kc_first3)
    rows, worst = [], 0.0
    cols = ('ctrl', 'tnt', 'ratio', 'lo90', 'hi90', 'n_ctrl_recordings', 'n_ctrl_flies', 'n_tnt_recordings', 'n_tnt_flies')
    for rd, rates, solv in (('carrying', all_t, False), ('first_three_trials', first, False),
                            ('solvent_subtracted', all_t, True)):
        for o in E104.FLIES_WEAK:
            lab = E104.flies_stim_label(o)
            mine = E104.flies_rate_ratio(rates, lab, solv)
            r = ref[(ref['odor'] == o) & (ref['reading'] == rd)].iloc[0]
            diff = {c: float(abs(float(mine[c]) - float(r[c]))) for c in cols}
            worst = max(worst, max(diff.values()))
            rows.append(dict(reading=rd, odor=o, **{f'{c}': mine[c] for c in cols}, **{f'{c}_ref': r[c] for c in cols},
                             max_abs_diff=max(diff.values())))
    return dict(passed=bool(worst == 0.0), max_abs_diff=worst, n_rows=len(rows)), pd.DataFrame(rows)


def model_identity(fam: pd.DataFrame, e102: pd.DataFrame) -> Tuple[dict, pd.DataFrame]:
    """The family's rates from the archived spike counts against E104's b_model.csv (m = 1 and m = 0: rate and
    rate_1spike, every wiring and stimulus) and E102's b_low_conc.csv (rate_on, rate_off per step, wiring and stimulus)
    and d_low_conc_inhadd1000.csv ('refit' rows); the saved responses against the counts at two spikes. 0 differences
    to pass."""
    E104 = _e104()
    rows = []
    bm = E104.read_table(H.shared_dir().joinpath(*E104_RESULTS, 'b_model.csv'))
    f = fam.set_index(['wiring', 'condition', 'stimulus'])
    for r in bm[bm['m'].isin([0.0, 1.0])].itertuples(index=False):
        cond = 'on' if r.m == 1.0 else 'off'
        mine = f.loc[(r.wiring, cond, r.stimulus)]
        rows.append(dict(check='E104 b_model', step='', wiring=r.wiring, condition=cond, stimulus=r.stimulus,
                         ref=float(r.rate), mine=float(mine['rate']), ref_1spike=float(r.rate_1spike),
                         mine_1spike=float(mine['rate_1spike'])))
    e = e102.set_index(['kind', 'step', 'wiring', 'condition', 'stimulus'])
    bl = read_table(H.shared_dir().joinpath(*E102_RESULTS, 'b_low_conc.csv'))
    for r in bl.itertuples(index=False):
        for cond, col in (('on', 'rate_on'), ('off', 'rate_off')):
            mine = e.loc[('ladder', r.setting, r.wiring, cond, r.stimulus)]
            rows.append(dict(check='E102 b_low_conc', step=r.setting, wiring=r.wiring, condition=cond,
                             stimulus=r.stimulus, ref=float(getattr(r, col)), mine=float(mine['rate'])))
    dl = read_table(H.shared_dir().joinpath(*E102_RESULTS, f'd_low_conc_{E102_STEP_OUTCOME}.csv'))
    for r in dl[dl['set'] == 'refit'].itertuples(index=False):
        for cond, col in (('on', 'rate_on'), ('off', 'rate_off')):
            mine = e.loc[('refit', E102_STEP_OUTCOME, r.wiring, cond, r.stimulus)]
            rows.append(dict(check='E102 d_low_conc refit', step=r.setting, wiring=r.wiring, condition=cond,
                             stimulus=r.stimulus, ref=float(getattr(r, col)), mine=float(mine['rate'])))
    t = pd.DataFrame(rows)
    t['differs'] = (t['ref'] != t['mine'])
    if 'ref_1spike' in t.columns:
        t['differs_1spike'] = t['ref_1spike'].notna() & (t['ref_1spike'] != t['mine_1spike'])
    n = int(t['differs'].sum() + t.get('differs_1spike', pd.Series(dtype=bool)).sum())
    resp = bool(fam['responses_equal'].all() and e102['responses_equal'].all())
    by = t.groupby('check')['differs'].agg(['size', 'sum']).to_dict('index')
    return dict(passed=bool(n == 0 and resp), n_differ=n, responses_equal_counts_at_2=resp,
                by_check={k: dict(n=int(v['size']), differ=int(v['sum'])) for k, v in by.items()}), t


def part_0(write: bool = False) -> dict:
    """Part 0: the identities, the sets' counts, c* and the unit check (p0_*.json / .csv); with `write`, then the lock
    (only if every identity passed). Raises if an identity fails (outcome (h))."""
    out = out_dir()
    trials = read_all_trials()
    counts = set_counts(trials)
    crow = [dict(set=s, group=g, n_recordings=counts.get((s, g), (0, 0))[0], n_flies=counts.get((s, g), (0, 0))[1],
                 entry_recordings=v[0], entry_flies=v[1]) for (s, g), v in SET_COUNTS.items()]
    ct = pd.DataFrame(crow)
    extra = sorted(set(counts) - set(SET_COUNTS))
    counts_ok = bool((ct['n_recordings'] == ct['entry_recordings']).all() and (ct['n_flies'] == ct['entry_flies']).all()
                     and not extra)
    write_table(ct, out / 'p0_counts.csv')
    weak = trials[trials['set'].isin(SETS)]
    kc_all, kc3 = kc_table(weak), kc_table(weak, FIRST_TRIALS)
    e49, e49t = e49_identity(kc_all)
    write_table(e49t, out / 'p0_e49_identity.csv')
    e104, e104t = e104_identity(kc_all, kc3)
    write_table(e104t, out / 'p0_e104_flies_identity.csv')
    fam = saved_family()
    e102 = e102_family()
    mid, midt = model_identity(fam, e102)
    write_table(midt, out / 'p0_model_identity.csv')
    # c* and the unit check
    apl = trials[trials['set'] == 'APL16']
    kc_apl = kc_table(apl)
    odors16 = apl16_odor_labels(apl)
    ctrl = kc_apl[kc_apl['genotype_group'] == 'control']
    z_panel = panel_rate(response_rates(ctrl, 'z1'), odors16)
    target = float(z_panel['rate'].mean())
    cs = find_cstar(ctrl, target, odors16)
    uc = unit_check(ctrl, cs['cstar'], odors16)
    write_table(uc['table'], out / 'p0_unit_check.csv')
    carrying_unit = 'trial_noise' if uc['passed'] else 'z1'
    partner_unit = 'z1' if uc['passed'] else 'trial_noise'
    tnt = kc_apl[kc_apl['genotype_group'] == 'apl_output_blocked']
    tnt_rates = {u: float(panel_rate(response_rates(tnt, u, cs['cstar']), odors16)['rate'].mean())
                 for u in ('trial_noise', 'z1')}
    cstar = dict(cs, odors=odors16, target_definition="the mean over the 11 control recordings of each recording's share "
                 "of KC x odor entries above 1 z over the 16 odors", tnt_rate_trial_noise=tnt_rates['trial_noise'],
                 tnt_rate_z1=tnt_rates['z1'], tnt_ratio_trial_noise=tnt_rates['trial_noise'] / cs['rate_at_cstar'],
                 tnt_ratio_z1=tnt_rates['z1'] / target)
    _write_json(out / 'p0_cstar.json', cstar)
    res = dict(counts=dict(passed=counts_ok, extra=[list(x) for x in extra]), e49=e49, e104_flies=e104, model=mid,
               cstar=cstar, unit_check={k: v for k, v in uc.items() if k != 'table'},
               carrying_unit=carrying_unit, partner_unit=partner_unit, host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
               definitions_sha256=definitions_sha(definitions()))
    res['passed'] = bool(counts_ok and e49['passed'] and e104['passed'] and mid['passed'])
    _write_json(out / 'p0_identity.json', res)
    files = {}
    for p in (H.shared_dir().joinpath(*WEAK_DIR, 'diag_mix', 'diag_mix_long.csv'),
              H.shared_dir().joinpath(*WEAK_DIR, 'natmix_dilute', 'natmix_dilute_long.csv'),
              H.shared_dir().joinpath(*APL16_DIR, 'TNT3c', 'dataset_20_long.csv'),
              H.shared_dir().joinpath(*APL16_DIR, 'TNT', 'dataset_19_long.csv'),
              H.shared_dir().joinpath(*E49_FRACTIONS), H.shared_dir().joinpath(*E104_RESULTS, 'b_flies.csv'),
              H.shared_dir().joinpath(*E104_RESULTS, 'b_model.csv'),
              H.shared_dir().joinpath(*E102_RESULTS, 'b_low_conc.csv')):
        files[str(p)] = _sha256(p)
    provenance('0', passed=res['passed'], write_lock=write, files=files)
    if not res['passed']:
        raise RuntimeError(f"part 0: an identity fails (outcome (h)): counts {counts_ok}, E49 {e49['passed']}, "
                           f"E104 flies {e104['passed']}, model {mid['passed']} (p0_identity.json)")
    if write:
        res['lock'] = dict(path=str(LOCK), sha256=write_lock(res)['sha256'])
    return res


# ---------------------------------------------------------------------------
# part A: the flies' side

READINGS = ('carrying',) + PARTNERS[:-1]       # the fly-left-out readings are 'fly_left_out:<fly>'


def set_cells(r: pd.DataFrame, spec: dict, stimuli: Sequence[str], div: Optional[Mapping[str, float]] = None
              ) -> pd.DataFrame:
    """Per set, group and stimulus of `stimuli` the flies' value and interval from per-recording rates `r` (`evoked`'s
    columns, with 'set') under `spec` (`reading_spec`): the control flies by `set_value` over their fly means (or, for
    the pooled reading, the recording mean of D1 and D2 with `recording_mean_value`), the blocked flies by
    `blocked_interval` from the same set's control interval, the recordings' range beside. Each set is read by itself;
    only `spec['pooled']` pools D1 and D2. `div`: per group, a divisor of every value (partner (5))."""
    col = 'evoked' if spec['evoked'] else 'rate'
    div = dict(div) if div is not None else {g: 1.0 for g in GROUPS}
    rows = []
    sets = (POOLED_SET,) if spec['pooled'] else [s for s in SETS if s in set(r['set'])]
    for s in sets:
        in_set = r['set'].isin(('D1', 'D2')) if s == POOLED_SET else (r['set'] == s)
        ctrl_iv = {}
        for g in GROUPS:
            for stim in stimuli:
                x = r[in_set & (r['genotype_group'] == g) & (r['odor'] == stim)]
                if g == 'control' and spec.get('drop_fly') is not None:
                    x = x[x['fly'] != spec['drop_fly']]
                if not len(x):
                    continue
                y = x.assign(v=x[col] / div[g])
                if g == 'control':
                    v = recording_mean_value(y, 'v') if spec['pooled'] else set_value(fly_values(y, 'v'))
                    lo, hi = v['lo90'], v['hi90']
                    ctrl_iv[stim] = (v['value'], lo, hi)
                else:
                    v = dict(value=float(y['v'].mean()) if spec['pooled'] else float(fly_values(y, 'v').mean()),
                             n_flies=int(y['fly'].nunique()))
                    c = ctrl_iv.get(stim)
                    lo, hi = blocked_interval(v['value'], *c) if c is not None else (np.nan, np.nan)
                rows.append(dict(set=s, group=g, stimulus=stim, value=v['value'], lo90=lo, hi90=hi,
                                 n_flies=v['n_flies'], n_recordings=int(len(y)),
                                 rec_min=float(y['v'].min()), rec_max=float(y['v'].max()),
                                 divisor=div[g], measure=col))
    return pd.DataFrame(rows)


def reading_spec(reading: str, lock: Mapping[str, object]) -> dict:
    """The flies' and the family's side of a reading: unit, first_n, evoked, relative, spikes, drop_fly, pooled."""
    cu, pu = lock['carrying_unit'], lock['partner_unit']
    base = dict(unit=cu, first_n=None, evoked=True, relative=False, spikes=SPIKES, drop_fly=None, pooled=False)
    if reading == 'carrying':
        return base
    if reading == 'unit':
        return dict(base, unit=pu)
    if reading == 'first_three_trials':
        return dict(base, first_n=FIRST_TRIALS)
    if reading == 'yang_definition':
        return dict(base, unit='yang')
    if reading == 'one_spike':
        return dict(base, spikes=SPIKES_PARTNER)
    if reading == 'relative_to_panel':
        return dict(base, relative=True)
    if reading.startswith('fly_left_out:'):
        return dict(base, drop_fly=reading.split(':', 1)[1])
    if reading == 'carried_so_far':
        return dict(base, unit='z1', evoked=False, pooled=True)
    raise ValueError(f'unknown reading {reading!r}')


class Flies:
    """The flies' per-recording rates under every unit and trial set, computed once (`rates`), and their set values
    (`cells`)."""

    def __init__(self, cstar: float):
        self.cstar = float(cstar)
        t = read_all_trials()
        self.kc = {None: kc_table(t), FIRST_TRIALS: kc_table(t, FIRST_TRIALS)}
        apl = t[t['set'] == 'APL16']
        self.odors16 = apl16_odor_labels(apl)
        self._rates: Dict[Tuple[str, Optional[int]], pd.DataFrame] = {}

    def rates(self, unit: str, first_n: Optional[int] = None) -> pd.DataFrame:
        k = (unit, first_n)
        if k not in self._rates:
            self._rates[k] = evoked(response_rates(self.kc[first_n], unit, self.cstar))
        return self._rates[k]

    def panel_values(self, unit: str, first_n: Optional[int] = None) -> Dict[str, dict]:
        """The APL16 flies' evoked panel value per group (mean over flies of each recording's mean over the 16 odors of
        its evoked rates) and the response rate beside."""
        r = self.rates(unit, first_n)
        p = panel_rate(r[r['set'] == 'APL16'], self.odors16)
        out = {}
        for g in GROUPS:
            x = p[p['genotype_group'] == g]
            out[g] = dict(evoked=set_value(fly_values(x, 'evoked')), rate=set_value(fly_values(x, 'rate')))
        return out

    def cells(self, spec: dict, stimuli: Sequence[str]) -> pd.DataFrame:
        """`set_cells` of this unit's and trials' rates under `spec`, divided by the APL16 groups' evoked panel values
        when `spec['relative']`."""
        r = self.rates(spec['unit'], spec['first_n'])
        div = None
        if spec['relative']:
            pv = self.panel_values(spec['unit'], spec['first_n'])
            div = {g: pv[g]['evoked']['value'] for g in GROUPS}
        return set_cells(r, spec, stimuli, div)

    def control_flies(self, s: str) -> List[str]:
        r = self.rates('z1')
        return sorted(r.loc[(r['set'] == s) & (r['genotype_group'] == 'control'), 'fly'].unique())


def model_values(fam: pd.DataFrame, spec: dict, condition: str, label: str) -> pd.Series:
    """The eight wirings' values at a stimulus label under `spec`: evoked rate at its spike criterion (or the rate as
    run for the carried reading), over the wiring's panel rate when `relative`."""
    x = fam[(fam['condition'] == condition) & (fam['stimulus'] == label)].set_index('wiring')
    one = spec['spikes'] == SPIKES_PARTNER
    if spec['evoked']:
        v = x['evoked_1spike' if one else 'evoked']
    else:
        v = x['rate_1spike' if one else 'rate']
    if spec['relative']:
        v = v / x['panel_rate_1spike' if one else 'panel_rate']
    return v.astype(float)


def cell_rows(fam: pd.DataFrame, flies: pd.DataFrame, spec: dict, reading: str,
              estimates: Sequence[str] = ESTIMATES, wirings: Sequence[str] = ALL_WIRINGS) -> List[dict]:
    """Part A's cells under one reading: per set, condition, estimate and weak odor the flies' interval, the wirings'
    values, the counts below, inside and above, the median and the class."""
    rows = []
    for s in flies['set'].unique():
        for cond, g in CONDITION_GROUP.items():
            for e in estimates:
                for o in WEAK_ODORS:
                    f = flies[(flies['set'] == s) & (flies['group'] == g) & (flies['stimulus'] == FLIES_STIM[o])]
                    if not len(f):
                        continue
                    f = f.iloc[0]
                    lab = model_label(o, e)
                    mv = model_values(fam, spec, cond, lab).reindex(list(wirings))
                    rows.append(_cell(reading, s, cond, e, o, FLIES_STIM[o], lab, f, mv))
    return rows


def _cell(reading, s, cond, e, o, fstim, lab, f, mv) -> dict:
    v = mv.to_numpy(float)
    lo, hi = float(f['lo90']), float(f['hi90'])
    return dict(reading=reading, set=s, condition=cond, estimate=e, odor=o, flies_stimulus=fstim, model_stimulus=lab,
                flies_value=float(f['value']), flies_lo90=lo, flies_hi90=hi, flies_n_flies=int(f['n_flies']),
                flies_n_recordings=int(f['n_recordings']), flies_rec_min=float(f['rec_min']),
                flies_rec_max=float(f['rec_max']),
                n_below=int((v < lo).sum()) if np.isfinite(lo) else -1,
                n_inside=int(((v >= lo) & (v <= hi)).sum()) if np.isfinite(lo) else -1,
                n_above=int((v > hi).sum()) if np.isfinite(hi) else -1, median=float(np.median(v)),
                cls=cell_class(v, lo, hi), **{f'w_{w}': float(x) for w, x in mv.items()})


def w_rows(cells: pd.DataFrame) -> pd.DataFrame:
    """W per reading, estimate, condition and set from the cells (`w_class`), with the three cells' classes."""
    rows = []
    for (rd, e, cond, s), g in cells.groupby(['reading', 'estimate', 'condition', 'set'], sort=False):
        cl = dict(zip(g['odor'], g['cls']))
        rows.append(dict(reading=rd, estimate=e, condition=cond, set=s, W=w_class(cl),
                         **{f'cls_{o}': cl.get(o, '') for o in WEAK_ODORS}))
    return pd.DataFrame(rows)


def held_rows(w: pd.DataFrame) -> pd.DataFrame:
    """Per estimate, condition and set: the carrying W, every partner's (fly_left_out as the list over left-out flies),
    the held W, the partners that differ and those that flip (short against not short)."""
    rows = []
    base = w[w['reading'] == 'carrying']
    for r in base.itertuples(index=False):
        x = w[(w['estimate'] == r.estimate) & (w['condition'] == r.condition) & (w['set'] == r.set)]
        partners = {}
        for p in PARTNERS[:-1]:
            v = x.loc[x['reading'] == p, 'W']
            if len(v):
                partners[p] = str(v.iloc[0])
        lo = x.loc[x['reading'].str.startswith('fly_left_out:'), 'W'].tolist()
        if lo:
            partners['fly_left_out'] = lo
        h = held_w(r.W, partners)
        rows.append(dict(estimate=r.estimate, condition=r.condition, set=r.set, carrying=r.W, held=h['held'],
                         differ=json.dumps(h['differ']), flips=json.dumps(h['flips']),
                         **{f'W_{p}': (json.dumps(v) if isinstance(v, list) else v) for p, v in partners.items()}))
    return pd.DataFrame(rows)


def a2_sets(fl: Flies, unit: str) -> pd.DataFrame:
    """A2: D1 against D2 and against N on each shared stimulus (control flies): per-fly rates (evoked and as recorded),
    E98's split rule (E103's `split_share`), the ratio of the set values with its interval from joint bootstrap draws of
    both sets' flies (E104's `q_of_fall`, arithmetic)."""
    from experiments import E103_apl_across_odors_heldout as E103
    E104 = _e104()
    r = fl.rates(unit)
    rows = []
    for other in ('D2', 'N'):
        for stim in ('ma @ -7', '2h @ -7', 'farn @ -3', SOLVENT):
            for col in ('evoked', 'rate'):
                if col == 'evoked' and stim == SOLVENT:
                    continue
                a = r[(r['set'] == 'D1') & (r['genotype_group'] == 'control') & (r['odor'] == stim)]
                b = r[(r['set'] == other) & (r['genotype_group'] == 'control') & (r['odor'] == stim)]
                if not len(a) or not len(b):
                    continue
                fa, fb = fly_values(a, col), fly_values(b, col)
                share, n = E103.split_share(fa.to_numpy(float), fb.to_numpy(float))
                q = E104.q_of_fall(fa, fb, arithmetic=True)
                rows.append(dict(unit=unit, measure=col, stimulus=stim, set_a='D1', set_b=other, value_a=float(fa.mean()),
                                 value_b=float(fb.mean()), n_flies_a=len(fa), n_flies_b=len(fb), split_share=share,
                                 n_splits=n, differ=bool(share < SPLIT_SHARE), ratio=q['q'], ratio_lo90=q['q_lo90'],
                                 ratio_hi90=q['q_hi90']))
    return pd.DataFrame(rows)


def a1_tables(fl: Flies, lock: Mapping[str, object]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A1: per set, group and stimulus (A1_STIMULI the set has) the response rate and the evoked rate with intervals,
    in the trial-noise unit, at 1 z and by Yang's definition, all trials and the first three; and the solvent's rate per
    recording."""
    parts = []
    for unit in ('trial_noise', 'z1', 'yang'):
        for first_n in (None, FIRST_TRIALS):
            for ev in (False, True):
                spec = dict(unit=unit, first_n=first_n, evoked=ev, relative=False, spikes=SPIKES, drop_fly=None,
                            pooled=False)
                stims = [s for s in A1_STIMULI if not (ev and s == SOLVENT)]
                c = fl.cells(spec, stims)
                parts.append(c.assign(unit=unit, trials='all' if first_n is None else f'first{first_n}',
                                      carrying=bool(unit == lock['carrying_unit'] and first_n is None and ev)))
    a1 = pd.concat(parts, ignore_index=True)
    a1['beside'] = a1['stimulus'].isin(list(LOWER_STIM.values()) + [b[1] for b in BESIDE_STIM])
    sol = []
    for unit in ('trial_noise', 'z1', 'yang'):
        for first_n in (None, FIRST_TRIALS):
            r = fl.rates(unit, first_n)
            x = r[r['odor'] == SOLVENT]
            sol.append(x[['set', 'recording', 'genotype_group', 'fly', 'rate', 'n_kcs']].assign(
                unit=unit, trials='all' if first_n is None else f'first{first_n}'))
    return a1, pd.concat(sol, ignore_index=True)


def beside_rows(fam: pd.DataFrame, fl: Flies, spec: dict) -> List[dict]:
    """A4: 2h at 10^-5 (exact input, kiwi-control's control panel) and ea at 10^-7 (input at 10^-8), set N, APL on."""
    rows = []
    f = fl.cells(spec, [b[1] for b in BESIDE_STIM])
    for o, fstim, lab, s, note in BESIDE_STIM:
        x = f[(f['set'] == s) & (f['group'] == 'control') & (f['stimulus'] == fstim)]
        if not len(x):
            continue
        mv = model_values(fam, spec, 'on', lab).reindex(list(ALL_WIRINGS))
        row = _cell('carrying', s, 'on', CARRIED_ESTIMATE, o, fstim, lab, x.iloc[0], mv)
        row['note'] = note
        row['safe'] = (row['cls'] in ('inside', 'over')) if o == 'ea' else True
        rows.append(row)
    return rows


def realism_rows(fam: pd.DataFrame) -> pd.DataFrame:
    """The eight fits' panel response rates (16 odors, two spikes) with APL on and acutely off, their ratio, against
    the KC data's (the registry's ctrl_frac and ratio), and the largest spike counts at the weak stimuli."""
    rows = []
    weak = [model_label(o, e) for o in WEAK_ODORS for e in ESTIMATES] + [b[2] for b in BESIDE_STIM]
    for w, g in fam.groupby('wiring', sort=False):
        on = g[g['condition'] == 'on'].iloc[0]['panel_rate']
        off = g[g['condition'] == 'off'].iloc[0]['panel_rate']
        rows.append(dict(wiring=w, rate_on=on, rate_off=off, ratio=off / on,
                         max_count_weak_on=float(g[(g['condition'] == 'on') & g['stimulus'].isin(weak)]['max_count'].max()),
                         max_count_weak_off=float(g[(g['condition'] == 'off') & g['stimulus'].isin(weak)]['max_count'].max()),
                         max_count_all_on=float(g[g['condition'] == 'on']['max_count'].max()),
                         max_count_all_off=float(g[g['condition'] == 'off']['max_count'].max())))
    t = pd.DataFrame(rows)
    try:
        t['kc_data_rate_on'] = float(H.data_target('noise_ctrl_frac'))
        t['kc_data_ratio'] = float(H.data_target('noise_tnt_frac')) / float(H.data_target('noise_ctrl_frac'))
    except Exception:                                    # the registry's keys (CLAUDE.md: 0.068 and 2.28)
        t['kc_data_rate_on'], t['kc_data_ratio'] = 0.068, 2.28
    return t


def part_a() -> dict:
    """Part A (no simulation; only once the lock is committed): a1_flies.csv, a1_solvent.csv, a2_sets.csv,
    a3_cells.csv, a3_w.csv, a3_held.csv, a4_beside.csv, a_family.csv, a_realism.csv and a_summary.json (the outcomes)."""
    lock = require_lock()
    out = out_dir()
    fl = Flies(lock['cstar'])
    fam = saved_family()
    write_table(fam, out / 'a_family.csv')
    a1, sol = a1_tables(fl, lock)
    write_table(a1, out / 'a1_flies.csv')
    write_table(sol, out / 'a1_solvent.csv')
    write_table(pd.concat([a2_sets(fl, lock['carrying_unit']), a2_sets(fl, lock['partner_unit'])], ignore_index=True),
                out / 'a2_sets.csv')
    cells = []
    stims = [FLIES_STIM[o] for o in WEAK_ODORS]
    # the reading carried so far, first
    spec = reading_spec('carried_so_far', lock)
    cells += cell_rows(fam, fl.cells(spec, stims), spec, 'carried_so_far', estimates=(CARRIED_ESTIMATE,))
    for rd in READINGS:
        spec = reading_spec(rd, lock)
        cells += cell_rows(fam, fl.cells(spec, stims), spec, rd)
    for s in SETS:
        for fly in fl.control_flies(s):
            rd = f'fly_left_out:{fly}'
            spec = reading_spec(rd, lock)
            c = fl.cells(spec, stims)
            # only this set's cells (and set D1's APL-off cells, whose blocked interval rests on D1's control flies)
            c = c[c['set'] == s]
            cells += cell_rows(fam, c, spec, rd)
    cells = pd.DataFrame(cells)
    write_table(cells, out / 'a3_cells.csv')
    w = w_rows(cells)
    write_table(w, out / 'a3_w.csv')
    held = held_rows(w[w['reading'] != 'carried_so_far'])
    write_table(held, out / 'a3_held.csv')
    spec = reading_spec('carrying', lock)
    write_table(pd.DataFrame(beside_rows(fam, fl, spec)), out / 'a4_beside.csv')
    real = realism_rows(fam)
    write_table(real, out / 'a_realism.csv')
    # beside: the unit check on the APL16 flies with APL output blocked (TNT) (the entry's check is the control flies')
    kc16 = fl.kc[None]
    uc_tnt = unit_check(kc16[(kc16['set'] == 'APL16') & (kc16['genotype_group'] == 'apl_output_blocked')],
                        lock['cstar'], fl.odors16, group='apl_output_blocked')
    write_table(uc_tnt['table'], out / 'a0_unit_check_tnt_beside.csv')
    oc = outcomes(held)
    pv = {u: fl.panel_values(u) for u in ('trial_noise', 'z1')}
    summ = dict(lock=lock, outcomes=oc, w_carried_so_far=w[w['reading'] == 'carried_so_far'][
        ['estimate', 'condition', 'set', 'W']].to_dict('records'),
        held=held[['estimate', 'condition', 'set', 'carrying', 'held', 'differ', 'flips']].to_dict('records'),
        apl16_panel_values=pv, realism=real.to_dict('records'),
        unit_check_tnt_beside={k: v for k, v in uc_tnt.items() if k != 'table'})
    _write_json(out / 'a_summary.json', summ)
    provenance('A', lock_sha256=json.loads(LOCK.read_text())['sha256'])
    return dict(outcomes=oc, held=summ['held'])


# ---------------------------------------------------------------------------
# part B: E102's cost at low concentration, read again

def part_b() -> dict:
    """Part B (no simulation): per step of E102's ladder (eight wirings) the cells of part A under the carrying reading
    (b_cells.csv), E102's class on evoked response rates per step, stimulus, condition, estimate and set against the
    unchanged fits (b_e102_class.csv), the refits described (b_refits.csv), and the reading of E102's outcome (e)
    (b_summary.json)."""
    lock = require_lock()
    out = out_dir()
    fl = Flies(lock['cstar'])
    spec = reading_spec('carrying', lock)
    flies = fl.cells(spec, [FLIES_STIM[o] for o in WEAK_ODORS])
    e102 = e102_family()
    rows = []
    for st in E102_STEPS:
        fam = e102[(e102['kind'] == 'ladder') & (e102['step'] == st)]
        rows += [dict(r, step=st) for r in cell_rows(fam, flies, spec, 'carrying')]
    cells = pd.DataFrame(rows)
    write_table(cells, out / 'b_cells.csv')
    ref = cells[cells['step'] == E102_UNCHANGED].set_index(['set', 'condition', 'estimate', 'odor'])
    cls = []
    for r in cells.itertuples(index=False):
        k = (r.set, r.condition, r.estimate, r.odor)
        u = float(ref.loc[k, 'median'])
        cls.append(dict(step=r.step, set=r.set, condition=r.condition, estimate=r.estimate, odor=r.odor,
                        flies_value=r.flies_value, median=r.median, unchanged_median=u,
                        log_distance=log_distance(r.median, r.flies_value),
                        unchanged_log_distance=log_distance(u, r.flies_value),
                        side=('above' if r.median > r.flies_value else 'below' if r.median < r.flies_value else 'at'),
                        e102_class=e102_class(r.median, u, r.flies_value),
                        e102_own_rule=e102_own_class(r.median, u, r.flies_value), cell_class=r.cls))
    cls = pd.DataFrame(cls)
    cls['distance_ratio'] = cls['log_distance'] / cls['unchanged_log_distance']
    write_table(cls, out / 'b_e102_class.csv')
    # the refits (described: four wirings), against the unchanged fits on the same four
    rf = e102[e102['kind'] == 'refit']
    unch = e102[(e102['kind'] == 'ladder') & (e102['step'] == E102_UNCHANGED)]
    fixed = e102[(e102['kind'] == 'ladder') & (e102['step'] == E102_STEP_OUTCOME)]
    rrows = []
    for r in flies.itertuples(index=False):
        for cond, g in CONDITION_GROUP.items():
            if r.group != g:
                continue
            for e in ESTIMATES:
                for o in WEAK_ODORS:
                    if FLIES_STIM[o] != r.stimulus:
                        continue
                    lab = model_label(o, e)
                    vr = model_values(rf, spec, cond, lab).reindex(list(S_WIRINGS)).to_numpy(float)
                    vu = model_values(unch, spec, cond, lab).reindex(list(S_WIRINGS)).to_numpy(float)
                    vf = model_values(fixed, spec, cond, lab).reindex(list(S_WIRINGS)).to_numpy(float)
                    mr, mu, mf = (float(np.median(v)) for v in (vr, vu, vf))
                    rrows.append(dict(set=r.set, condition=cond, estimate=e, odor=o, flies_value=r.value,
                                      flies_lo90=r.lo90, flies_hi90=r.hi90, n_below=int((vr < r.lo90).sum()),
                                      n_inside=int(((vr >= r.lo90) & (vr <= r.hi90)).sum()),
                                      n_above=int((vr > r.hi90).sum()), median_refit=mr, median_unchanged_S=mu,
                                      median_fixed_S=mf, e102_class_refit=e102_class(mr, mu, r.value),
                                      e102_class_fixed_S=e102_class(mf, mu, r.value),
                                      e102_own_rule_refit=e102_own_class(mr, mu, r.value),
                                      **{f'w_{w}': float(x) for w, x in zip(S_WIRINGS, vr)}))
    rr = pd.DataFrame(rrows)
    write_table(rr, out / 'b_refits.csv')
    ma = cls[(cls['step'] == E102_STEP_OUTCOME) & (cls['odor'] == 'ma')]
    stands = bool(len(ma) and (ma['e102_class'] == 'further').all())
    where = ma[ma['e102_class'] == 'further'][['condition', 'estimate', 'set']].to_dict('records')
    not_where = ma[ma['e102_class'] != 'further'][['condition', 'estimate', 'set', 'e102_class', 'side']].to_dict('records')
    summ = dict(stands_like_for_like=stands, further_under=where, not_further_under=not_where,
                n_combinations=int(len(ma)),
                ma_inhadd1000=ma[['condition', 'estimate', 'set', 'flies_value', 'median', 'unchanged_median',
                                  'distance_ratio', 'side', 'e102_class', 'e102_own_rule', 'cell_class']].to_dict('records'),
                refits_ma=rr[rr['odor'] == 'ma'][['condition', 'estimate', 'set', 'median_refit', 'median_unchanged_S',
                                                  'e102_class_refit', 'e102_own_rule_refit']].to_dict('records'),
                max_count=float(e102['max_count'].max()))
    _write_json(out / 'b_summary.json', summ)
    provenance('B', lock_sha256=json.loads(LOCK.read_text())['sha256'])
    return summ


# ---------------------------------------------------------------------------
# part C: E104's ratios within set D1

def part_c() -> dict:
    """Part C (no simulation; described): the flies' ratio of response rates, APL output blocked (TNT) over set D1's
    control flies alone (E104's `flies_rate_ratio`), as recorded and evoked, at 1 z (E104's unit) and in the carrying
    unit; d* from E104's saved ladder (b_model.csv) at those ratios and Y2's class by E104's own rules, beside E104's
    verdict (c_ratios.csv, c_dstar.csv, c_summary.json)."""
    lock = require_lock()
    E104 = _e104()
    out = out_dir()
    fl = Flies(lock['cstar'])
    ratios = []
    for unit in dict.fromkeys(('z1', lock['carrying_unit'])):
        r = fl.rates(unit)
        r = r[r['set'] == 'D1'][['recording', 'genotype_group', 'fly', 'odor', 'rate']].copy()
        for o in E104.FLIES_WEAK:
            lab = E104.flies_stim_label(o)
            for solv in (False, True):
                x = E104.flies_rate_ratio(r, lab, solv)
                ratios.append(dict(unit=unit, odor=o, reading='evoked' if solv else 'as_recorded', **x))
    ratios = pd.DataFrame(ratios)
    write_table(ratios, out / 'c_ratios.csv')
    bm = E104.read_table(H.shared_dir().joinpath(*E104_RESULTS, 'b_model.csv'))
    rows = []
    for t in ratios.itertuples(index=False):
        for e in E104.X2_ESTIMATE_SESSIONS:
            lab = E104.model_label(t.odor, e)
            for w, g in bm[bm['stimulus'] == lab].groupby('wiring', sort=False):
                ds = E104.d_star(g['d'], g['ratio'], t.ratio, t.lo90, t.hi90)
                d1 = float(g.loc[g['m'] == 1.0, 'd'].iloc[0])
                rows.append(dict(unit=t.unit, reading=t.reading, odor=t.odor, estimate=e, wiring=w, flies_ratio=t.ratio,
                                 flies_lo90=t.lo90, flies_hi90=t.hi90, **ds, d_m1=d1,
                                 family_ratio_m1=float(g.loc[g['m'] == 1.0, 'ratio'].iloc[0]),
                                 y2=E104.y_class((d1, d1), (ds['d_star_lo'], ds['d_star_hi']))))
    dt = pd.DataFrame(rows)
    write_table(dt, out / 'c_dstar.csv')
    verdicts = {}
    for (unit, rd), g in dt.groupby(['unit', 'reading'], sort=False):
        per_odor = {}
        for o, go in g.groupby('odor', sort=False):
            pe = {e: E104.y_over_wirings(ge['y2']) for e, ge in go.groupby('estimate', sort=False)}
            per_odor[o] = dict(per_estimate=pe, over_estimates=E104.y_over_estimates(pe),
                               d_star_median=float(go['d_star'].median()))
        verdicts[f'{unit}|{rd}'] = dict(per_odor=per_odor,
                                        verdict=E104.y_verdict({o: v['over_estimates'] for o, v in per_odor.items()}))
    e104 = _read_json(H.shared_dir().joinpath(*E104_RESULTS, 'b_summary.json')) or {}
    summ = dict(y2_within_d1=verdicts, e104_y2=dict(carrying=e104.get('y2', {}).get('carrying'),
                                                    held=e104.get('y2', {}).get('held')),
                ratios=ratios[['unit', 'odor', 'reading', 'ctrl', 'tnt', 'ratio', 'lo90', 'hi90', 'n_ctrl_flies',
                               'n_tnt_flies']].to_dict('records'),
                note="described: E104's verdicts stand as its entry fixed them")
    _write_json(out / 'c_summary.json', summ)
    provenance('C', lock_sha256=json.loads(LOCK.read_text())['sha256'])
    return summ


# ---------------------------------------------------------------------------
# part D (simulations at fixed parameters; beside: no class rests on these)

D_RUNS = 'D/runs'
D3_DIR = 'D/D3'
D_TAG = 'D_{wiring}'                 # the part-D runs' variant names (D1 and D2 read the same runs)
D3_TAG = 'D3_real'
D_STEM = '{tag}_{cond}'              # E98's B1 layout under outputs/E105_weak_odor_rates/D/runs/
ADDED_PANEL = 'E105'
MIXES = (('2h', 'ma'), ('2h', 'farn'), ('ma', 'farn'))
EA4_LABEL = 'ea [kiwi-control/kiwi] @ -4.2'
D2_MAPS = ('linear', 'saturating')   # E73's kept map for C2 is linear (map_choice.csv); the saturating beside
D2_RULES = ('as_blank', 'as_above')  # 2h at 10^-8 and farn at 10^-4: as the blank, or as the stimulus a log unit above
FACTORS_ALL = tuple(sorted(set(FACTORS) | {1.0}))


def _e98():
    return _e104()._e98()


def scaled_label(odor: str, estimate: str, factor: float) -> str:
    """A D1 column: a weak stimulus's input under an estimate times `factor` on every glomerulus (one ' @ ')."""
    o, c = MODEL_STIM[odor]
    return f'{o} [{estimate}/diag x{factor:g}] @ {c:g}'


def mix_label(a: str, b: str, estimate: str) -> str:
    """A mixture column (the per-glomerulus maximum of the two components' inputs under an estimate)."""
    (oa, ca), (ob, cb) = MODEL_STIM[a], MODEL_STIM[b]
    return f'{oa}{ca:g}+{ob}{cb:g} [{estimate}/diag max] @ 0'


def ea8_label(estimate: str) -> str:
    E98 = _e98()
    return E98.session_label('ea', -8, estimate, E98.DIAG)


def d_input() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Part D's input (the lab's 36 glomeruli): the 16 panel odors (unchanged, first), the nine weak stimuli under the
    three estimates (B1's columns), each times FACTORS (D1), the mixtures by the per-glomerulus maximum under each
    estimate (D2), ea at 10^-8 under each estimate's diagnostics, 2h at 10^-5 (kiwi-control's control panel), ea at
    10^-4.2 (kiwi-control's kiwi panel) and the blank. Returns (input, table: label, role, odor, estimate, factor,
    source, the input total and its positive part)."""
    E98 = _e98()
    b1, _ = E98.b1_input()
    labels = E98.input_labels(b1)
    col = {lab: b1.iloc[:, i].astype(float) for i, lab in enumerate(labels)}
    panel = H.apl16_orn_deltas()
    gl = list(panel.index)
    sess = E98.session_inputs()
    rows, extra = [], {}

    def add(label, v, **kw):
        extra[label] = pd.Series(np.asarray(v, float), index=panel.index)
        rows.append(dict(label=label, input_total=float(np.sum(v)), input_total_pos=float(np.clip(v, 0, None).sum()), **kw))

    for o in H.APL16_ODORS:
        rows.append(dict(label=o, role='panel', odor=o, estimate='', factor=1.0, source='panel',
                         input_total=float(col[o].sum()), input_total_pos=float(col[o].clip(lower=0).sum())))
    for e in ESTIMATES:
        for o in WEAK_ODORS:
            lab = model_label(o, e)
            add(lab, col[lab].to_numpy(), role='weak', odor=o, estimate=e, factor=1.0, source=f'B1:{lab}')
            for f in FACTORS:
                add(scaled_label(o, e, f), f * col[lab].to_numpy(), role='scaled', odor=o, estimate=e, factor=float(f),
                    source=f'{f:g} x B1:{lab}')
        for a, b in MIXES:
            add(mix_label(a, b, e), np.maximum(col[model_label(a, e)].to_numpy(), col[model_label(b, e)].to_numpy()),
                role='mixture', odor=f'{a}+{b}', estimate=e, factor=1.0, source=f'max(B1:{model_label(a, e)}, '
                f'B1:{model_label(b, e)})')
        lab = ea8_label(e)
        if lab in col:
            v = col[lab].to_numpy()
        else:
            raw = sess[e][(E98.DIAG, 'ea @ -8')]
            v = raw.reindex(gl).fillna(0.0).astype(float).to_numpy()
        add(lab, v, role='ea8', odor='ea', estimate=e, factor=1.0, source=f'{e}/diag ea @ -8')
    lab = BESIDE_STIM[0][2]
    add(lab, col[lab].to_numpy(), role='beside', odor='2h', estimate='kiwi-control', factor=1.0, source=f'B1:{lab}')
    raw = sess['kiwi-control'][('kiwi', 'ea @ -4.2')]
    add(EA4_LABEL, raw.reindex(gl).fillna(0.0).astype(float).to_numpy(), role='ea4.2', odor='ea',
        estimate='kiwi-control', factor=1.0, source='kiwi-control/kiwi ea @ -4.2')
    add(BLANK, np.zeros(len(gl)), role='blank', odor='blank', estimate='', factor=1.0, source='blank')
    ex = pd.DataFrame(extra, index=panel.index)
    ex.columns = pd.MultiIndex.from_tuples([(ADDED_PANEL, k) for k in ex.columns], names=panel.columns.names)
    orn = pd.concat([panel.astype(float), ex], axis='columns')
    orn.index.name = panel.index.name
    return orn, pd.DataFrame(rows)


def _clean_commits() -> dict:
    c = _commits()
    if c.get('queue_dirty') or c.get('pyolfsysm_dirty'):
        raise RuntimeError(f'part D runs from clean checkouts only: {c}')
    return c


def _save_run(d: Path, stem: str, o: dict) -> None:
    for k in ('responses', 'spike_counts'):
        x = o[k].copy()
        x.columns = [str(c) for c in x.columns]
        x.to_parquet(d / f'{stem}_{k}.parquet')
    o['apl'].to_csv(d / f'{stem}_apl.csv')


def run_d(wirings: Sequence[str] = ALL_WIRINGS, force: bool = False) -> pd.DataFrame:
    """D1's and D2's runs (resumable): per wiring, `d_input` at its fit, APL on and acutely off (E98's `ExtendedRun`),
    saved under outputs/E105_weak_odor_rates/D/runs/ (D_<wiring>_<on|acute_off>_{spike_counts,responses}.parquet,
    _apl.csv) with a runs.csv row each (every parameter, absolute, the commits, threads). The columns the saved B1 run
    has must reproduce it exactly (0 count differences), else the run stops."""
    import gc
    E98 = _e98()
    d = local_dir(D_RUNS)
    commits = _clean_commits()
    orn, tab = d_input()
    tab.to_csv(d / 'd_stimuli.csv', index=False)
    for w in wirings:
        tag = D_TAG.format(wiring=w)
        if not force and all((d / f'{D_STEM.format(tag=tag, cond=ct)}_spike_counts.parquet').exists()
                             for ct in CONDITION_RUN.values()):
            continue
        t0 = time.time()
        er = E98.ExtendedRun(w, orn)
        saved = _e104().saved_run(w)
        for cond, ct in CONDITION_RUN.items():
            t1 = time.time()
            o = er.run('on' if cond == 'on' else 'off')
            s_new = o['spike_counts'].astype(float)
            s_old = saved[ct]['spike_counts'].astype(float)
            common = [c for c in s_old.columns if c in s_new.columns]
            same_idx = bool(np.array_equal(np.asarray(s_new.index), np.asarray(s_old.index)))
            nd = int((s_new[common].to_numpy() != s_old[common].to_numpy()).sum()) if same_idx else -1
            ident = dict(wiring=w, condition=ct, n_common=len(common), index_equal=same_idx, counts_differ=nd)
            _write_json(d / f'{D_STEM.format(tag=tag, cond=ct)}_identity.json', ident)
            if nd != 0:
                raise RuntimeError(f'{w} {ct}: the part-D run does not reproduce the saved B1 run ({ident})')
            stem = D_STEM.format(tag=tag, cond=ct)
            _save_run(d, stem, o)
            p = er.params('on' if cond == 'on' else 'off')
            row = dict(tag=tag, wiring=w, condition=ct, n_stimuli=int(s_new.shape[1]), n_common_with_b1=len(common),
                       max_spike_count=float(s_new.to_numpy().max()),
                       params_json=json.dumps({k: float(v) for k, v in p.items()}, sort_keys=True),
                       **{k: float(v) for k, v in p.items()}, runtime_s=round(time.time() - t1, 1), host=_host(),
                       time=time.strftime('%Y-%m-%d %H:%M'), threads=os.environ.get('OMP_NUM_THREADS'),
                       thread_limit=os.environ.get('OMP_THREAD_LIMIT'), **commits)
            pd.DataFrame([row]).to_csv(d / 'runs.csv', mode='a', header=not (d / 'runs.csv').exists(), index=False)
            print(f"{time.strftime('%H:%M:%S')} {w} {ct}: {time.time() - t1:.0f} s, {len(common)} B1 columns reproduced",
                  flush=True)
            del o
        del er
        gc.collect()
        print(f"{time.strftime('%H:%M:%S')} {w}: done in {time.time() - t0:.0f} s", flush=True)
    return pd.read_csv(d / 'runs.csv')


def d_runs(wirings: Sequence[str] = ALL_WIRINGS) -> Dict[str, Dict[str, dict]]:
    """The saved part-D runs ({wiring: {'on', 'acute_off': {'spike_counts', 'responses', 'apl'}}})."""
    d = local_dir(D_RUNS)
    out = {}
    for w in wirings:
        tag = D_TAG.format(wiring=w)
        out[w] = {}
        for ct in CONDITION_RUN.values():
            stem = D_STEM.format(tag=tag, cond=ct)
            out[w][ct] = dict(spike_counts=pd.read_parquet(d / f'{stem}_spike_counts.parquet'),
                              responses=pd.read_parquet(d / f'{stem}_responses.parquet'),
                              apl=pd.read_csv(d / f'{stem}_apl.csv', index_col=0), dir=str(d))
    return out


def _log_inverse(x: np.ndarray, y: np.ndarray, target: float) -> Tuple[float, str]:
    """The factor at which y (made non-decreasing in x, interpolated linearly in log x) reaches `target`: the midpoint of
    the first and last such x (E104's `_inverse`); censored below or above the ladder."""
    E104 = _e104()
    o = np.argsort(x)
    lx, f = np.log(np.asarray(x, float)[o]), E104.isotonic(np.asarray(y, float)[o], True)
    if not np.isfinite(target):
        return np.nan, 'no target'
    a, b = E104._inverse(lx, f, target, 'left'), E104._inverse(lx, f, target, 'right')
    if np.isinf(a) and a < 0:
        return float(np.exp(lx[0])), f'below {np.exp(lx[0]):g}'
    if np.isinf(b) and b > 0:
        return float(np.exp(lx[-1])), f'above {np.exp(lx[-1]):g}'
    return float(np.exp(0.5 * (a + b))), ''


def d1_tables(lock: Mapping[str, object], fam: pd.DataFrame, tab: pd.DataFrame, fl: 'Flies') -> Tuple[pd.DataFrame,
                                                                                                         pd.DataFrame]:
    """D1: per wiring, condition, estimate, odor and factor the evoked response rate (two spikes) and the largest count
    (d1_rates.csv); per set, condition, estimate and odor the factor at which each wiring's rate meets the flies'
    evoked value (carrying reading), their median, and the estimates' input totals and the factor between them
    (d1_factors.csv)."""
    spec = reading_spec('carrying', lock)
    flies = fl.cells(spec, [FLIES_STIM[o] for o in WEAK_ODORS])
    rows = []
    for e in ESTIMATES:
        for o in WEAK_ODORS:
            for f in FACTORS_ALL:
                lab = model_label(o, e) if f == 1.0 else scaled_label(o, e, f)
                x = fam[fam['stimulus'] == lab]
                for r in x.itertuples(index=False):
                    rows.append(dict(estimate=e, odor=o, factor=f, stimulus=lab, wiring=r.wiring, condition=r.condition,
                                     evoked=r.evoked, rate=r.rate, max_count=r.max_count))
    rates = pd.DataFrame(rows)
    tot = tab[tab['role'] == 'weak'].set_index(['odor', 'estimate'])
    frows = []
    for fr in flies.itertuples(index=False):
        cond = 'on' if fr.group == 'control' else 'off'
        for e in ESTIMATES:
            for o in WEAK_ODORS:
                if FLIES_STIM[o] != fr.stimulus:
                    continue
                g = rates[(rates['estimate'] == e) & (rates['odor'] == o) & (rates['condition'] == cond)]
                facs, cens = {}, {}
                for w, gw in g.groupby('wiring', sort=False):
                    facs[w], cens[w] = _log_inverse(gw['factor'].to_numpy(), gw['evoked'].to_numpy(), fr.value)
                v = np.array([facs[w] for w in ALL_WIRINGS])
                tots = {s: float(tot.loc[(o, s), 'input_total']) for s in ESTIMATES}
                frows.append(dict(set=fr.set, condition=cond, estimate=e, odor=o, flies_value=fr.value,
                                  flies_lo90=fr.lo90, flies_hi90=fr.hi90, median_factor=float(np.median(v)),
                                  n_censored_below=sum(c.startswith('below') for c in cens.values()),
                                  n_censored_above=sum(c.startswith('above') for c in cens.values()),
                                  input_total=tots[e], input_total_over_lowest=tots[e] / min(tots.values()),
                                  sessions_factor_max=max(tots.values()) / min(tots.values()),
                                  max_count=float(g['max_count'].max()),
                                  **{f'factor_{w}': facs[w] for w in ALL_WIRINGS},
                                  **{f'censored_{w}': cens[w] for w in ALL_WIRINGS}))
    return rates, pd.DataFrame(frows)


# D2: the sessions as the recordings were made

def stim_input(stim: str, estimate: str, rule: str) -> Optional[List[str]]:
    """The model columns whose per-glomerulus maximum is a flies' stimulus's input under an estimate, or None (no
    input: the stimulus enters the session as a blank). Tom's rule (within one log unit) where the concentration has no
    ORN data: ma at 10^-8 takes ma at 10^-7; 2h at 10^-7 takes 2h at 10^-6; farn at 10^-3 and 10^-2.5 take farn at
    10^-2; ea at 10^-7 takes ea at 10^-8 (the estimate's diagnostics); ea at 10^-5 takes ea at 10^-4.2
    (kiwi-control's kiwi panel); IaA at 10^-2 takes the panel's IaA (10^-3). 2h at 10^-8 and farn at 10^-4 (no ORN data
    within one log unit): the blank under rule 'as_blank', the stimulus a log unit above under 'as_above'. 1o3ol at
    10^-5, eb at 10^-7 and 10^-5, ban, Bmix28, air and the solvent have none. A mixture: its components'."""
    if ' + ' in stim:
        parts = [stim_input(c.strip(), estimate, rule) for c in stim.split(' + ')]
        cols = [c for p in parts if p for c in p]
        return cols or None
    single = {'2h @ -7': model_label('2h', estimate), 'ma @ -7': model_label('ma', estimate),
              'ma @ -8': model_label('ma', estimate), 'farn @ -3': model_label('farn', estimate),
              'farn @ -2.5': model_label('farn', estimate), '2h @ -5': BESIDE_STIM[0][2], 'ea @ -7': ea8_label(estimate),
              'ea @ -5': EA4_LABEL, 'IaA @ -2': 'IaA'}
    if stim in ('2h @ -8', 'farn @ -4'):
        if rule == 'as_blank':
            return None
        return [model_label('2h' if stim.startswith('2h') else 'farn', estimate)]
    if stim in single:
        return [single[stim]]
    return None


def recording_designs(trials: pd.DataFrame) -> pd.DataFrame:
    """Per recording of the three sets: its set, group, fly and design (the stimuli it presents with their trial
    numbers, as a JSON list of [stimulus, n_trials] pairs, sorted)."""
    t = trials[trials['set'].isin(SETS)]
    n = t.groupby(['set', 'recording', 'genotype_group', 'fly', 'odor'], sort=True)['repeat'].nunique().reset_index()
    rows = []
    for (s, rec, g, fly), x in n.groupby(['set', 'recording', 'genotype_group', 'fly'], sort=True):
        rows.append(dict(set=s, recording=rec, genotype_group=g, fly=fly,
                         design=json.dumps(sorted([[o, int(k)] for o, k in zip(x['odor'], x['repeat'])]))))
    return pd.DataFrame(rows)


def _session_columns(design: Sequence[Sequence], estimate: str, rule: str, counts: pd.DataFrame
                     ) -> Tuple[pd.DataFrame, List[dict], int, int]:
    """A design's session for `session_zscore`: the stimuli with input as columns (mixtures as the maximum of the
    components' inputs: the D2 mixture columns of the run where both components have input), each repeated by its trial
    number over the design's greatest common divisor of trial numbers (so that the session's means weight the stimuli
    by their trials), and the blank stimuli's count in those units. Returns (counts, per-stimulus records with the
    column and trial number, the unit of trials, the blank's multiplicity)."""
    from math import gcd
    from functools import reduce
    g = reduce(gcd, [int(k) for _, k in design])
    cols, recs, n_blank = [], [], 0
    for stim, k in design:
        inp = stim_input(stim, estimate, rule)
        m = int(k) // g
        if inp is None:
            n_blank += m
            recs.append(dict(stimulus=stim, n_trials=int(k), column=None, input='none (blank)'))
            continue
        if len(inp) == 1:
            col = inp[0]
        else:
            pair = tuple(next(o for o in WEAK_ODORS if model_label(o, estimate) == c) for c in inp)
            col = mix_label(*pair, estimate) if pair in MIXES else mix_label(*pair[::-1], estimate)
        recs.append(dict(stimulus=stim, n_trials=int(k), column=col, input=' max '.join(inp)))
        cols += [col] * m
    c = counts[cols].copy()
    c.columns = [f'{x}#{i}' for i, x in enumerate(cols)]
    return c, recs, g, n_blank


def d2_reading(counts: pd.DataFrame, design: Sequence[Sequence], estimate: str, rule: str, cal: dict,
               q: float, drift: float) -> List[dict]:
    """One design's session read as E73 reads a recording (`session_zscore` with the wiring's panel calibration, carried
    unchanged): per stimulus the expected share of KCs whose trial mean (over the stimulus's trials) is above 1 z."""
    c, recs, g, n_blank = _session_columns(design, estimate, rule, counts)
    z = H.session_zscore(c, q=q, noise=cal['noise'], amp=cal['amp'], K=cal['K'], noise_ratio=cal['noise_ratio'],
                         n_trials=g, n_blank=n_blank, drift=drift)
    out = []
    first = {}
    for j, name in enumerate(c.columns):
        first.setdefault(name.split('#')[0], j)
    for r in recs:
        sd = (z['trial_sd'] / np.sqrt(r['n_trials']))[:, :, None]
        if r['column'] is None:
            C = z['blank'][:, :, None]
        else:
            C = z['center'][:, :, [first[r['column']]]]
        P = H._tail_moments(C, sd, Z1)[0]
        out.append(dict(r, rate_z=float(P.mean()),
                        rate_two_spikes=float((counts[r['column']] >= SPIKES).mean()) if r['column'] else 0.0))
    return out


def d2_tables(lock: Mapping[str, object], runs: Mapping[str, Mapping[str, dict]], fl: 'Flies'
              ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """D2: per set, design, estimate, rule, map, wiring, condition and stimulus the model's rate under E73's per-KC
    session reading (d2_rates.csv); the calibration per wiring and map (d2_calibration.csv); per set, condition,
    estimate, rule, map and stimulus the set's value (mean over its flies of their recordings' designs) as recorded
    and evoked against the flies' at 1 z, with the wirings' values and the class (described; d2_cells.csv)."""
    from experiments import E73_per_kc_normalization as E73
    dt = H.data_target
    q = E73.q_values()[0]
    drift = E73.drift_ratio()
    targets = dict(rate=dt('ctrl_frac'), amp=dt('ctrl_resp_amp_mean_z'), silent_sd=dt('ctrl_trial_noise_sd'))
    designs = recording_designs(read_all_trials())
    cal_rows, rows = [], []
    for w, run in runs.items():
        on16 = run['on']['spike_counts'][list(H.APL16_ODORS)].astype(float)
        for a_map in D2_MAPS:
            cal = H.calibrate_session_zscore(on16, q=q, amp=a_map, target_rate=targets['rate'],
                                             target_amp=targets['amp'] if a_map == 'saturating' else None,
                                             target_silent_sd=targets['silent_sd'], drift=drift)
            cal = dict(cal, amp=a_map)
            cal_rows.append(dict(wiring=w, map=a_map, q=q, drift=drift, **{k: cal.get(k) for k in (
                'noise', 'K', 'noise_ratio', 'rate', 'resp_amp_mean', 'nonresp_trial_sd', 'noise_at_bound', 'K_at_bound',
                'converged')}))
            for cond, ct in CONDITION_RUN.items():
                counts = run[ct]['spike_counts'].astype(float)
                g = CONDITION_GROUP[cond]
                for design, x in designs[designs['genotype_group'] == g].groupby('design', sort=False):
                    dsg = json.loads(design)
                    for e in ESTIMATES:
                        for rule in D2_RULES:
                            for r in d2_reading(counts, dsg, e, rule, cal, q, drift):
                                rows.append(dict(r, wiring=w, map=a_map, condition=cond, estimate=e, rule=rule,
                                                 design=design, sets=';'.join(sorted(set(x['set']))),
                                                 recordings=';'.join(x['recording'])))
    rates = pd.DataFrame(rows)
    # the set values: each recording takes its design's value; a fly the mean of its recordings; the set the flies' mean
    sol = rates[rates['stimulus'] == SOLVENT][['wiring', 'map', 'condition', 'estimate', 'rule', 'design', 'rate_z']]
    rates = rates.merge(sol.rename(columns={'rate_z': 'solvent_z'}), on=['wiring', 'map', 'condition', 'estimate', 'rule',
                                                                        'design'], how='left')
    rates['evoked_z'] = rates['rate_z'] - rates['solvent_z']
    spec_r = dict(unit='z1', first_n=None, evoked=False, relative=False, spikes=SPIKES, drop_fly=None, pooled=False)
    spec_e = dict(spec_r, evoked=True)
    stims = sorted(set(rates['stimulus']) - {SOLVENT})
    fr, fe = fl.cells(spec_r, stims + [SOLVENT]), fl.cells(spec_e, stims)
    crows = []
    for (s, cond), d in designs.groupby(['set', 'genotype_group'], sort=False):
        c = 'on' if cond == 'control' else 'off'
        for (e, rule, a_map), rr in rates[rates['condition'] == c].groupby(['estimate', 'rule', 'map'], sort=False):
            for stim in sorted(set(rr['stimulus'])):
                vals = {}
                for meas in ('rate_z', 'evoked_z'):
                    per_rec = []
                    for rec in d.itertuples(index=False):
                        v = rr[(rr['design'] == rec.design) & (rr['stimulus'] == stim)]
                        if len(v):
                            per_rec.append(dict(fly=rec.fly, **v.set_index('wiring')[meas].reindex(list(ALL_WIRINGS))
                                                .rename(lambda w: f'w_{w}').to_dict()))
                    if not per_rec:
                        continue
                    pr = pd.DataFrame(per_rec).groupby('fly').mean()
                    vals[meas] = pr.mean()
                if not vals:
                    continue
                for meas, ff in (('rate_z', fr), ('evoked_z', fe)):
                    if meas not in vals or (meas == 'evoked_z' and stim == SOLVENT):
                        continue
                    f = ff[(ff['set'] == s) & (ff['group'] == cond) & (ff['stimulus'] == stim)]
                    v = vals[meas].reindex([f'w_{w}' for w in ALL_WIRINGS]).to_numpy(float)
                    lo = float(f['lo90'].iloc[0]) if len(f) else np.nan
                    hi = float(f['hi90'].iloc[0]) if len(f) else np.nan
                    crows.append(dict(set=s, condition=c, estimate=e, rule=rule, map=a_map, stimulus=stim, measure=meas,
                                      flies_value=float(f['value'].iloc[0]) if len(f) else np.nan, flies_lo90=lo,
                                      flies_hi90=hi, median=float(np.median(v)), cls_described=cell_class(v, lo, hi),
                                      **{f'w_{w}': x for w, x in zip(ALL_WIRINGS, v)}))
    return rates, pd.DataFrame(cal_rows), pd.DataFrame(crows)


# D3: the trial noise on the real wiring

D3_STIMULI = tuple(H.APL16_ODORS) + tuple(model_label(o, e) for e in ESTIMATES for o in WEAK_ODORS) + (
    BESIDE_STIM[0][2], 'ea [kiwi-control/diag] @ -8', BLANK)


def run_d3(force: bool = False) -> dict:
    """D3: the real wiring at its fit through C2's matched threshold jitter (E72's part F: E49's `matched_noise`
    ('c2_68'), E49's `noisy_trials` scheme: thr = max(thr0 (1 + sigma xi), 0.01 thr0), thr0 = fixed_thr x ratio + 2
    spont_in, one xi per KC and trial, the same thresholds with APL on and acutely off), D3_RECORDINGS pseudo-recordings
    of D3_TRIALS trials, on the panel's 16 odors, the nine weak stimuli under the three estimates, 2h at 10^-5, ea at
    10^-8 and the blank (the part-D input's columns). Saves outputs/E105_weak_odor_rates/D/D3/D3_real_trials.npz (uint16
    spike counts, conditions x recordings x trials x KCs x stimuli) and D3_real.json (sigma, the ratio, the absolute
    thresholds and scales, the seed)."""
    from experiments import E49_low_conc_kc as E49
    from experiments import E43_model_trial_noise as E43
    E98 = _e98()
    d = local_dir(D3_DIR)
    npz = d / f'{D3_TAG}_trials.npz'
    if npz.exists() and not force:
        return _read_json(d / f'{D3_TAG}.json')
    commits = _clean_commits()
    orn, _ = d_input()
    labs = E98.input_labels(orn)
    keep = [i for i, lab in enumerate(labs) if lab in set(D3_STIMULI)]
    orn = orn.iloc[:, keep]
    er = E98.ExtendedRun('real', orn)
    labels = list(er.labels)
    mn = E49.matched_noise('c2_68')
    sigma, ratio = float(mn['sigma']), float(mn['fixed_thr_ratio'])
    spont = np.asarray(er.rv.kc.spont_in, dtype=float).ravel()
    f0 = float(er.p['fixed_thr']) * ratio
    thr0 = f0 + 2.0 * spont
    rng = np.random.default_rng(SEED)
    n_kcs = spont.size
    spk = np.zeros((2, D3_RECORDINGS, D3_TRIALS, n_kcs, len(labels)), dtype=np.uint16)
    t0 = time.time()
    for rep in range(D3_RECORDINGS):
        for tr in range(D3_TRIALS):
            xi = rng.standard_normal(n_kcs)
            thr = np.maximum(thr0 * (1.0 + sigma * xi), 0.01 * thr0)
            with E43.kc_threshold_override(thr):
                for i, cond in enumerate(('on', 'off')):
                    p = dict(er.params(cond), fixed_thr=f0)
                    _, s = er.sim(p)
                    o = E98.relabel(s, s, labels)
                    spk[i, rep, tr] = o['spike_counts'].astype(float).to_numpy().astype(np.uint16)
        print(f"{time.strftime('%H:%M:%S')} D3 pseudo-recording {rep + 1} of {D3_RECORDINGS} ({time.time() - t0:.0f} s)",
              flush=True)
    np.savez_compressed(npz, spike_counts=spk, conditions=np.array(['on', 'acute_off']), stimuli=np.array(labels),
                        kcs=np.array([str(x) for x in o['spike_counts'].index]))
    rec = dict(tag=D3_TAG, wiring='real', sigma=sigma, fixed_thr_ratio=ratio, fixed_thr_fit=float(er.p['fixed_thr']),
               fixed_thr=f0, thr_min_frac=0.01, noise_source=mn['source'], seed=SEED, n_recordings=D3_RECORDINGS,
               n_trials=D3_TRIALS, stimuli=labels, params_on=dict(er.params('on'), fixed_thr=f0),
               params_off=dict(er.params('off'), fixed_thr=f0), max_spike_count=int(spk.max()),
               runtime_s=round(time.time() - t0, 1), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
               threads=os.environ.get('OMP_NUM_THREADS'), thread_limit=os.environ.get('OMP_THREAD_LIMIT'), **commits)
    _write_json(d / f'{D3_TAG}.json', rec)
    return rec


def d3_tables(lock: Mapping[str, object], fl: 'Flies') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """D3's readings: per condition, pseudo-recording, trials (all six, the first three) and stimulus the response rate
    by the trial-mean spike count at one or more ('count1', E49's data rule) and at two or more ('count2'), and per
    trial at two spikes ('trial2'), each also evoked (minus the blank's), and over the panel's 16 odors
    (d3_recordings.csv); per condition, estimate, odor and reading their mean and SD over pseudo-recordings beside the
    flies' values per set (carrying reading, and 1 z as recorded; d3_summary.csv)."""
    d = local_dir(D3_DIR)
    z = np.load(d / f'{D3_TAG}_trials.npz')
    spk = z['spike_counts'].astype(float)
    labels = [str(x) for x in z['stimuli']]
    bi = labels.index(BLANK)
    p16 = [labels.index(o) for o in H.APL16_ODORS]
    rows = []
    for i, cond in enumerate(('on', 'off')):
        for rep in range(spk.shape[1]):
            for nt in (D3_TRIALS, FIRST_TRIALS):
                x = spk[i, rep, :nt]
                m = x.mean(axis=0)
                rd = {'count1': (m >= 1).mean(axis=0), 'count2': (m >= 2).mean(axis=0),
                      'trial2': (x >= SPIKES).mean(axis=(0, 1))}
                for name, v in rd.items():
                    for j, lab in enumerate(labels):
                        rows.append(dict(condition=cond, recording=rep, n_trials=nt, reading=name, stimulus=lab,
                                         rate=float(v[j]), evoked=float(v[j] - v[bi]), panel_rate=float(v[p16].mean())))
    t = pd.DataFrame(rows)
    spec = reading_spec('carrying', lock)
    fc = fl.cells(spec, [FLIES_STIM[o] for o in WEAK_ODORS])
    fz = fl.cells(dict(spec, unit='z1', evoked=False), [FLIES_STIM[o] for o in WEAK_ODORS])
    srows = []
    for (cond, nt, name), g in t.groupby(['condition', 'n_trials', 'reading'], sort=False):
        for e in ESTIMATES:
            for o in WEAK_ODORS:
                x = g[g['stimulus'] == model_label(o, e)]
                grp = CONDITION_GROUP[cond]
                row = dict(condition=cond, n_trials=nt, reading=name, estimate=e, odor=o, stimulus=model_label(o, e),
                           rate_mean=float(x['rate'].mean()), rate_sd=float(x['rate'].std(ddof=1)),
                           evoked_mean=float(x['evoked'].mean()), evoked_sd=float(x['evoked'].std(ddof=1)),
                           panel_rate_mean=float(x['panel_rate'].mean()))
                for s in SETS:
                    a = fc[(fc['set'] == s) & (fc['group'] == grp) & (fc['stimulus'] == FLIES_STIM[o])]
                    b = fz[(fz['set'] == s) & (fz['group'] == grp) & (fz['stimulus'] == FLIES_STIM[o])]
                    row[f'flies_{s}_carrying'] = float(a['value'].iloc[0]) if len(a) else np.nan
                    row[f'flies_{s}_z1_as_recorded'] = float(b['value'].iloc[0]) if len(b) else np.nan
                srows.append(row)
    return t, pd.DataFrame(srows)


def write_variants(runs: pd.DataFrame, d3: Optional[dict]) -> pd.DataFrame:
    """results/E105_weak_odor_rates/variants.csv, which the archive builds params.csv from: one row per saved run (tag
    = the files' variant), with every factor, the absolute threshold and scales, all parameters (base_params, JSON),
    the stimuli and the largest spike counts."""
    rows = []
    for (tag, w), g in runs.drop_duplicates(['tag', 'condition'], keep='last').groupby(['tag', 'wiring'], sort=False):
        on = g[g['condition'] == 'on'].iloc[0]
        off = g[g['condition'] == 'acute_off'].iloc[0]
        rows.append(dict(tag=tag, variant=tag, part='D1/D2', wiring=w, ref_experiment='E72_wiring_population (part B '
                         'fits; E83 wiring labels), E98 ExtendedRun', factor=';'.join(f'{f:g}' for f in FACTORS_ALL),
                         scaled_odors=';'.join(f'{model_label(o, e)}' for e in ESTIMATES for o in WEAK_ODORS),
                         fixed_thr=on['fixed_thr'], wAPLKC_scale=on['wAPLKC_scale'], wKCAPL_scale=on['wKCAPL_scale'],
                         wAPLPN_scale=on['wAPLPN_scale'], wPNAPL_scale=on['wPNAPL_scale'],
                         kc_kc_inh=on.get('wKCKC_scale', np.nan), base_params=on['params_json'],
                         stimuli=f"{int(on['n_stimuli'])} columns: outputs/{EXPERIMENT}/{D_RUNS}/d_stimuli.csv",
                         max_count_on=on['max_spike_count'], max_count_off=off['max_spike_count'], host=on['host'],
                         queue_commit=on.get('queue_commit'), pyolfsysm_commit=on.get('pyolfsysm_commit')))
    if d3:
        p = d3['params_on']
        rows.append(dict(tag=d3['tag'], variant=d3['tag'], part='D3', wiring='real', ref_experiment='E72 part F matched '
                         f"threshold jitter ({d3['noise_source']})", fixed_thr=d3['fixed_thr'],
                         thr_factor=d3['fixed_thr_ratio'], wAPLKC_scale=p['wAPLKC_scale'], wKCAPL_scale=p['wKCAPL_scale'],
                         wAPLPN_scale=p['wAPLPN_scale'], wPNAPL_scale=p['wPNAPL_scale'],
                         kc_kc_inh=p.get('wKCKC_scale', np.nan), seed=d3['seed'],
                         base_params=json.dumps(dict(sigma=d3['sigma'], fixed_thr_ratio=d3['fixed_thr_ratio'],
                                                     fixed_thr_fit=d3['fixed_thr_fit'], thr_min_frac=d3['thr_min_frac'],
                                                     n_recordings=d3['n_recordings'], n_trials=d3['n_trials'],
                                                     params_on=d3['params_on'], params_off=d3['params_off']),
                                                sort_keys=True),
                         stimuli=';'.join(d3['stimuli']), max_count_on=d3['max_spike_count'], host=d3['host'],
                         queue_commit=d3.get('queue_commit'), pyolfsysm_commit=d3.get('pyolfsysm_commit')))
    v = pd.DataFrame(rows)
    v.to_csv(out_dir() / 'variants.csv', index=False)
    return v


def part_d(sub: str) -> dict:
    """Part D (`sub`): 'D1' (the runs D1 and D2 read: `run_d`), 'D3' (`run_d3`), 'score' (no simulation: d1_rates.csv,
    d1_factors.csv, d2_rates.csv, d2_calibration.csv, d2_cells.csv, d3_recordings.csv, d3_summary.csv, d_family.csv,
    variants.csv and d_summary.json), 'all' (the three in turn)."""
    lock = require_lock()
    out = out_dir()
    res = {}
    if sub in ('D1', 'D2', 'all'):
        runs = run_d()
        res['runs'] = int(len(runs))
        provenance('D_runs', model=H.model_provenance(orn_deltas=d_input()[0]), threads=os.environ.get('OMP_NUM_THREADS'))
    if sub in ('D3', 'all'):
        res['d3'] = {k: v for k, v in run_d3().items() if k not in ('stimuli',)}
        provenance('D3_runs', threads=os.environ.get('OMP_NUM_THREADS'))
    if sub in ('score', 'all'):
        fl = Flies(lock['cstar'])
        runs = d_runs()
        fam = family_table(runs, source='E105:D/runs')
        write_table(fam, out / 'd_family.csv')
        _, tab = d_input()
        write_table(tab, out / 'd_stimuli.csv')
        rates, factors = d1_tables(lock, fam, tab, fl)
        write_table(rates, out / 'd1_rates.csv')
        write_table(factors, out / 'd1_factors.csv')
        r2, cal, c2 = d2_tables(lock, runs, fl)
        write_table(r2, out / 'd2_rates.csv')
        write_table(cal, out / 'd2_calibration.csv')
        write_table(c2, out / 'd2_cells.csv')
        d3 = _read_json(local_dir(D3_DIR) / f'{D3_TAG}.json')
        if d3 is not None:
            t3, s3 = d3_tables(lock, fl)
            write_table(t3, out / 'd3_recordings.csv')
            write_table(s3, out / 'd3_summary.csv')
            _write_json(out / 'd3_run.json', {k: v for k, v in d3.items()})
        runs_csv = pd.read_csv(local_dir(D_RUNS) / 'runs.csv')
        write_variants(runs_csv, d3)
        summ = dict(factors_median=factors[['set', 'condition', 'estimate', 'odor', 'flies_value', 'median_factor',
                                            'n_censored_below', 'n_censored_above', 'sessions_factor_max',
                                            'max_count']].to_dict('records'),
                    max_count_d1=float(rates['max_count'].max()), max_count_all=float(fam['max_count'].max()),
                    d2_calibration=cal.to_dict('records'), d3=None if d3 is None else
                    dict(sigma=d3['sigma'], fixed_thr_ratio=d3['fixed_thr_ratio'], max_spike_count=d3['max_spike_count']))
        _write_json(out / 'd_summary.json', summ)
        provenance('D_score', lock_sha256=json.loads(LOCK.read_text())['sha256'])
        res['score'] = dict(n_factor_rows=int(len(factors)), n_d2_cells=int(len(c2)))
    return res


# ---------------------------------------------------------------------------
# figure (from part A's saved table)

FIG_ESTIMATE_COLORS = {'kiwi-control': '#2a78d6', 'megamat': '#eb6834', 'validation2': '#1baf7a'}   # categorical 1-3
FIG_PANELS = (('carried_so_far', POOLED_SET, 'on', 'the reading carried so far, control flies D1+D2 / APL on'),
              ('carried_so_far', POOLED_SET, 'off', 'the reading carried so far, APL output blocked (TNT) / APL off'),
              ('carrying', 'D1', 'on', 'set D1, control flies / APL on'),
              ('carrying', 'D1', 'off', 'set D1, APL output blocked (TNT) / APL off'),
              ('carrying', 'D2', 'on', 'set D2, control flies / APL on'),
              ('carrying', 'N', 'on', 'set N, control flies / APL on'))


def figure_cells(path: Optional[Path] = None) -> Path:
    """E105_cells.png: per panel (FIG_PANELS) and weak odor the flies' value (line) and 90 % interval (band) against
    the eight wirings' values (dots) and their median (bar) under each estimate; the carried reading first. From
    a3_cells.csv."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    c = read_table(out_dir() / 'a3_cells.csv')
    path = out_dir() / 'E105_cells.png' if path is None else Path(path)
    fig, axes = plt.subplots(3, 2, figsize=(11, 10.5), sharey=False)
    ink, muted = '#2b2b29', '#6b6a64'
    for ax, (rd, s, cond, title) in zip(axes.ravel(), FIG_PANELS):
        x = c[(c['reading'] == rd) & (c['set'] == s) & (c['condition'] == cond)]
        ests = [e for e in ESTIMATES if e in set(x['estimate'])]
        for i, o in enumerate(WEAK_ODORS):
            g = x[x['odor'] == o]
            if not len(g):
                continue
            f = g.iloc[0]
            x0, x1 = i - 0.42, i + 0.42
            ax.fill_between([x0, x1], f['flies_lo90'], f['flies_hi90'], color='#d9d8d2', lw=0, zorder=1)
            ax.plot([x0, x1], [f['flies_value']] * 2, color=muted, lw=1.5, zorder=2)
            for j, e in enumerate(ests):
                r = g[g['estimate'] == e]
                if not len(r):
                    continue
                r = r.iloc[0]
                xe = i + (j - (len(ests) - 1) / 2) * 0.26
                v = [r[f'w_{w}'] for w in ALL_WIRINGS]
                jit = np.linspace(-0.06, 0.06, len(v))
                ax.scatter(xe + jit, v, s=22, color=FIG_ESTIMATE_COLORS[e], edgecolor='white', linewidth=0.8, zorder=3)
                ax.plot([xe - 0.09, xe + 0.09], [r['median']] * 2, color=ink, lw=2, zorder=4)
                ax.annotate(r['cls'], (xe, 0), xycoords=('data', 'axes fraction'), xytext=(0, 3),
                            textcoords='offset points', ha='center', va='bottom', fontsize=7, color=muted, rotation=90)
        ax.set_xticks(range(len(WEAK_ODORS)))
        ax.set_xticklabels([f'{o}\n(flies {FLIES_STIM[o].split(" @ ")[1]})' for o in WEAK_ODORS], fontsize=9, color=ink)
        ax.set_title(title, fontsize=9.5, color=ink, loc='left')
        meas = 'response rate' if rd == 'carried_so_far' else 'evoked response rate'
        ax.set_ylabel(meas, fontsize=9, color=ink)
        ax.axhline(0, color='#bdbcb5', lw=0.8, zorder=0)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=muted, labelsize=8)
        ax.set_xlim(-0.6, len(WEAK_ODORS) - 0.4)
    handles = [Line2D([], [], marker='o', ls='', color=FIG_ESTIMATE_COLORS[e], label=f"{e}'s estimate of the input")
               for e in ESTIMATES]
    handles += [Line2D([], [], color=ink, lw=2, label='median of the eight wirings'),
                Line2D([], [], color=muted, lw=1.5, label="the flies' value (band: 90 % interval)")]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=8.5, frameon=False)
    fig.suptitle("The family's KCs at the weak odors against the flies', set by set (carrying reading: evoked, "
                 "trial-noise unit, all trials;\nfamily: two spikes or more). Labels: the cell's class", fontsize=10,
                 color=ink, x=0.02, ha='left')
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'C', 'D', 'figures'])
    ap.add_argument('--write-lock', action='store_true')
    ap.add_argument('--sub', default='all', choices=['D1', 'D2', 'D3', 'score', 'all'])
    a = ap.parse_args()
    if a.part == '0':
        r = part_0(write=a.write_lock)
    elif a.part == 'A':
        r = part_a()
    elif a.part == 'B':
        r = part_b()
    elif a.part == 'C':
        r = part_c()
    elif a.part == 'figures':
        r = dict(cells=str(figure_cells()))
    else:
        r = part_d(a.sub)
    print(json.dumps(r, indent=1, default=_json_default))


if __name__ == '__main__':
    main()
