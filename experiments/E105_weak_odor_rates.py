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
read from the saved runs (E98's B1 for the real wiring, E101's archived A3 for the other seven; E102's archived ladder and
refits for part B; E104's saved ladder table for part C).

The flies, in three sets that are never pooled (`SETS`):
- D1: the diagnostic-mixture recordings to 2025-02-24 (9 control recordings of 7 flies; the 3 recordings of 2 flies with
  APL output blocked (TNT)); six trials, the solvent three.
- D2: the same panel from 2025-03-20 (7 control recordings of 5 flies); six trials (one recording four).
- N: the dilute natural-mixture recordings of 2026-01-16 to 02-10 (9 control recordings of 5 flies); two to four trials.
- APL16 (11 control recordings, 8 with APL output blocked (TNT)): for the unit and for the panel's values.

Definitions (the entry's; all in `definitions()` and hashed into the lock, LOCK, before any class is computed):
- A response of a fly's KC: its trial mean above the criterion. The trial-noise unit: the trial mean over the pooled SD of
  the KC's single trials around their own stimulus means (every stimulus the recording presents, the solvent included;
  `pooled_sd`). Its criterion c* is the value at which the APL16 control flies' response rate (the mean over their 11
  recordings of each recording's share of KC x odor entries above it, over the 16 odors: the registry's definition of
  0.068) equals theirs at 1 z; fixed in part 0 (`find_cstar`) and applied unchanged.
- The evoked response rate: a recording's response rate to a stimulus minus its response rate to the solvent, in the
  same unit and trials.
- A set's value: the mean over flies (a fly is a date and a fly number; its recordings averaged), with a 90 % interval
  from 10,000 bootstrap draws of flies.
- The flies with APL output blocked (TNT): their value times and over the relative half-width, on the log scale, of set
  D1's control flies at the same stimulus (`blocked_interval`); the three recordings' range beside.
- The family's response: two spikes or more; its response rate minus the blank's is its evoked response rate (the blank
  gives none).
- The class of a cell (stimulus, condition, estimate, set; `cell_class`): 'short' if at least 6 of 8 wirings lie below
  the flies' interval, 'over' if at least 6 above, 'inside' if at least 6 inside or if neither holds and the eight
  wirings' median lies inside; otherwise 'not decided'. APL on against control flies, APL off against APL output blocked
  (TNT) (set D1 only).
- W per estimate, condition and set (`w_class`): 'short' if ma is short, or 2h and farn are both short; 'not short' if ma
  is inside or over and 2h and farn are not both short; otherwise 'not decided'.
- The carrying reading: evoked response rates, the trial-noise unit, all trials. Held (`held_w`): W the same under the
  carrying reading and under each partner (PARTNERS): (1) 1 z in place of the trial-noise unit; (2) the first three
  trials alone; (3) Yang's definition, above 1 z on at least half the trials; (4) one spike or more as the family's
  response; (5) each side's evoked response rate over its own panel's; (6) each control fly of the set left out in turn.
  If part 0's unit check fails, 1 z carries and the trial-noise unit is partner (1).
- Beside in every table, as its first row: the reading carried so far (the solvent not removed, 1 z, all trials, D1 and
  D2 pooled as recording means, kiwi-control's estimate).

Parts (CLI `--part`):
- 0: the identities (p0_identity.json; a failure stops the experiment, outcome (h)); c* and the unit check; with
  `--write-lock` the lock (LOCK), which must be committed before part A computes a class (`require_lock`).
- A: the flies against the family by set and by estimate (no simulation): A1 the flies' side, A2 the sets against each
  other (E98's split rule over flies), A3 the classes and W under the carrying reading and each partner, A4 beside (2h at
  10^-5, ea at 10^-7).
- B: E102's cost at low concentration read again (no simulation): the cells of part A per step of E102's ladder and for
  its four refits, E102's class on evoked response rates (`e102_class`).
- C: E104's ratios within set D1 (no simulation; described): d* from E104's saved ladder, Y2's class by E104's rules.
- D: beside, at fixed parameters (simulations; no class rests on these): D1 the input's size (FACTORS), D2 the family read
  as the recordings were made (E73's per-KC session reading), D3 the trial noise on the real wiring (C2's matched
  threshold jitter, E72's part F).

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
SETS = ('D1', 'D2', 'N')
SET_RULES = {'D1': ('diag_mix', None, '2025-02-24'), 'D2': ('diag_mix', '2025-03-20', None),
             'N': ('natmix_dilute', '2026-01-16', '2026-02-10')}    # (panel, first date, last date), inclusive
SET_COUNTS = {('D1', 'control'): (9, 7), ('D1', 'apl_output_blocked'): (3, 2), ('D2', 'control'): (7, 5),
              ('N', 'control'): (9, 5), ('APL16', 'control'): (11, 8), ('APL16', 'apl_output_blocked'): (8, 7)}
GROUPS = ('control', 'apl_output_blocked')
CONDITION_GROUP = {'on': 'control', 'off': 'apl_output_blocked'}
SOLVENT = 'pfo @ 0'

# ---------------------------------------------------------------------------
# the stimuli

ESTIMATES = ('kiwi-control', 'megamat', 'validation2')
WEAK_ODORS = ('ma', '2h', 'farn')
FLIES_STIM = {'ma': 'ma @ -7', '2h': '2h @ -7', 'farn': 'farn @ -3'}
MODEL_STIM = {'ma': ('ma', -7.0), '2h': ('2h', -6.0), 'farn': ('farn', -2.0)}   # the input's (E98's B1 diagnostics)
LOWER_STIM = {'ma': 'ma @ -8', '2h': '2h @ -8', 'farn': 'farn @ -4'}         # set D1, beside (A1)
BESIDE_STIM = (('2h @ -5', '2h [kiwi-control/control] @ -5', 'N', 'exact input, one estimate'),
               ('ea @ -7', 'ea [kiwi-control/diag] @ -8', 'N', 'input a log unit lower: only inside or over is safe'))
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
FACTORS = (0.5, 0.71, 1.41, 2.0, 2.83, 4.0)   # part D1
D3_RECORDINGS, D3_TRIALS = 10, 6

LOCK = Path(__file__).resolve().parent / 'E105_lock.json'
IDENTITY_TOL = 1e-12          # the flies' reader against the saved fractions (floating-point summation order)

QUALITY_CAVEAT = ("Yang's processed per-trial values (data/yang_weak_odor_kc/, data/yang_apl16_kc_olfdb/); she has not "
                  "been asked about these recordings' quality beyond her answers of 2026-11-15 18:08")
UNMEASURED = ("unmeasured: z-scored nuclear GCaMP6s against spike counts; a chronic block (TNT) against the model's acute "
              "APL off; two flies with APL output blocked; other flies and sessions for the ORN input than for the KCs; "
              "the weak-odor recordings' baseline-noise units (traces not on hand)")


def out_dir() -> Path:
    """results/E105_weak_odor_rates/ in the shared directory."""
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    """outputs/E105_weak_odor_rates/<sub>/ in this repository."""
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# the reader

def check_path(path: Path) -> Path:
    """`path` resolved, if it lies under the shared directory; PermissionError otherwise."""
    raise NotImplementedError


def raw_only_stems() -> set:
    """The recordings data/yang_apl_recordings/recordings.csv marks raw only (use_until_further_notice)."""
    raise NotImplementedError


def check_rows(df: pd.DataFrame, allowed: Sequence[str], raw_only: Optional[set] = None) -> pd.DataFrame:
    """`df` unchanged if every row's recording is in `allowed` (the panel's processed list) and none is raw only
    (`raw_only`, or a 'use_until_further_notice' column saying so); PermissionError otherwise."""
    raise NotImplementedError


def assign_set(df: pd.DataFrame, panel: str) -> pd.Series:
    """Each row's set (SETS) from its panel and date (SET_RULES); a row in no set raises."""
    raise NotImplementedError


def read_trials(panel: str, path: Optional[Path] = None) -> pd.DataFrame:
    """A panel's per-trial values (diag_mix, natmix_dilute or APL16): recording, set, genotype_group, date, fly_num,
    fly, roi, odor, repeat, value; through `check_path` and `check_rows`."""
    raise NotImplementedError


def set_counts(trials: pd.DataFrame) -> Dict[Tuple[str, str], Tuple[int, int]]:
    """{(set, group): (recordings, flies)}."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the units and the rates

def pooled_sd(trials: pd.DataFrame) -> pd.Series:
    """Per recording and KC, the pooled SD of its single trials around their own stimulus means (sum of squares over
    stimuli over the trials minus the stimuli)."""
    raise NotImplementedError


def kc_table(trials: pd.DataFrame, first_n: Optional[int] = None) -> pd.DataFrame:
    """Per recording, KC and stimulus: the trial mean (z), the number of trials, the share of trials above 1 z, the KC's
    pooled SD and the trial mean in the trial-noise unit (all from the same trials: all, or the first `first_n`)."""
    raise NotImplementedError


def response_rates(kc: pd.DataFrame, unit: str, cstar: Optional[float] = None) -> pd.DataFrame:
    """Per recording and stimulus, the response rate under `unit`: 'trial_noise' (trial mean in the trial-noise unit
    above `cstar`), 'z1' (trial mean above 1 z) or 'yang' (above 1 z on at least half the trials)."""
    raise NotImplementedError


def evoked(rates: pd.DataFrame) -> pd.DataFrame:
    """`rates` with 'evoked': each recording's rate minus its solvent's rate."""
    raise NotImplementedError


def panel_rate(rates: pd.DataFrame, odors: Sequence[str]) -> pd.Series:
    """Per recording the share of KC x odor entries above the criterion over `odors` (mean of the per-odor rates)."""
    raise NotImplementedError


def find_cstar(kc: pd.DataFrame, target: float, odors: Sequence[str]) -> float:
    """c*: the trial-noise criterion at which the mean over recordings of `panel_rate` equals `target` (bisection over
    CSTAR_BRACKET to CSTAR_TOL)."""
    raise NotImplementedError


def unit_check(kc_ctrl: pd.DataFrame, cstar: float) -> dict:
    """The APL16 control flies' 16 per-odor response rates (mean over recordings) in the trial-noise unit at c* against
    those in baseline-noise units (data/yang_apl16_kc_noise_units/response_fraction_per_recording.csv): RMS and passed."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the flies' values

def fly_values(rec_values: pd.DataFrame, col: str) -> pd.Series:
    """Per fly, the mean of its recordings' `col`."""
    raise NotImplementedError


def set_value(per_fly: pd.Series, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The mean over flies with its 90 % interval from bootstrap draws of flies."""
    raise NotImplementedError


def blocked_interval(value: float, ctrl_value: float, ctrl_lo: float, ctrl_hi: float) -> Tuple[float, float]:
    """The blocked flies' interval: `value` over and times exp(h), h the half-width of set D1's control flies' interval
    on the log scale ((ln hi - ln lo) / 2)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the classes

def cell_class(model: Sequence[float], lo: float, hi: float, min_wirings: int = MIN_WIRINGS) -> str:
    """A cell's class (CELL_CLASSES) from the eight wirings' values and the flies' interval (inclusive)."""
    raise NotImplementedError


def w_class(classes: Mapping[str, str]) -> str:
    """W (W_CLASSES) from the cells of ma, 2h and farn (farn absent in set N: 'both short' cannot hold)."""
    raise NotImplementedError


def held_w(carrying: str, partners: Mapping[str, object]) -> dict:
    """Whether W is held: the same under the carrying reading and every partner (a partner's value may be a list, one W
    per left-out fly, each of which must agree). Returns the held W ('not held' otherwise) and the partners that
    differ."""
    raise NotImplementedError


def outcomes(w: pd.DataFrame) -> dict:
    """Outcomes (a) to (g) from part A's W table (and part B's reading for (g))."""
    raise NotImplementedError


def e102_class(model_med: float, unchanged_med: float, flies: float, further: float = E102_FURTHER,
               nearer: float = E102_NEARER) -> str:
    """Part B's class: 'further' if |ln(model / flies)| is at least `further` times the unchanged fits', 'nearer' if at
    most `nearer` times, 'level' otherwise."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the lock

def definitions() -> dict:
    """The definitions, lists and thresholds (module constants only) as one JSON-able dict."""
    raise NotImplementedError


def definitions_sha(defs: dict) -> str:
    raise NotImplementedError


def write_lock(part0: dict, path: Optional[Path] = None) -> dict:
    """The lock: the definitions, part 0's c* and unit check (which unit carries), their sha256; refuses to overwrite a
    lock with other definitions."""
    raise NotImplementedError


def locked(path: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether the lock exists, matches `definitions()` and its sha256, and is committed."""
    raise NotImplementedError


def require_lock(**kw) -> dict:
    """The lock's contents; PermissionError unless `locked(**kw)`."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

def part_0(write: bool = False) -> dict:
    raise NotImplementedError


def part_a() -> dict:
    raise NotImplementedError


def part_b() -> dict:
    raise NotImplementedError


def part_c() -> dict:
    raise NotImplementedError


def part_d(sub: str) -> dict:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['0', 'A', 'B', 'C', 'D'])
    ap.add_argument('--write-lock', action='store_true')
    ap.add_argument('--sub', default='all', choices=['D1', 'D2', 'D3', 'score', 'all'])
    a = ap.parse_args()
    if a.part == '0':
        print(json.dumps(part_0(write=a.write_lock), indent=1, default=str))
    elif a.part == 'A':
        print(json.dumps(part_a(), indent=1, default=str))
    elif a.part == 'B':
        print(json.dumps(part_b(), indent=1, default=str))
    elif a.part == 'C':
        print(json.dumps(part_c(), indent=1, default=str))
    else:
        print(json.dumps(part_d(a.sub), indent=1, default=str))


if __name__ == '__main__':
    main()
