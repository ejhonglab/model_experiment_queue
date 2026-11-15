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
  reading those two recordings at another concentration than the other ten.
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
SUBSETS = ('all', 'pulse_1s', 'pulse_0.5s', 'without_large_blank_and_3FRfun1')
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

def stimulus_responses(recs) -> pd.DataFrame:
    """Per recording, stimulus, baseline (BASELINES) and blank (subtracted or not): the calyx means of the two
    presentations and of their mean (E79's `stimulus_maps`, doubtful presentations left out), the recording's blank
    calyx mean, its fly, side, pulse and the stimulus's block."""
    raise NotImplementedError


def blank_identity(resp: pd.DataFrame) -> float:
    """The largest |first-in-block with the blank subtracted - order-matched without it| over recordings and stimuli
    (zero up to rounding: the two are the same reading)."""
    raise NotImplementedError


def value_at(concs: Sequence[float], values: Sequence[float], c: float,
             design_range: Optional[Tuple[float, float]] = None) -> Tuple[float, str]:
    """A recording's response at log10 concentration `c` from its imaged concentrations `concs` and responses `values`:
    ('imaged') if imaged; ('interpolated') linear on log10 concentration between its two imaged neighbours inside its
    imaged range; ('substituted') the nearest imaged end when `c` is outside the series' `design_range` (default the
    recording's own range) and within MATCH_LOG_UNITS of that end, inclusive; otherwise NaN ('outside', or 'missing'
    when `c` is inside the design range but outside the recording's own)."""
    raise NotImplementedError


def concentration_match(c: float, imaged: Sequence[float], max_units: float = MATCH_LOG_UNITS) -> Optional[float]:
    """The imaged concentration nearest to `c` if it is within `max_units` log10 units (inclusive), else None."""
    raise NotImplementedError


def fly_ratio(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """The flies' ratio from per-recording log ratios (columns 'fly', 'log_ratio'): per fly the mean on the log scale,
    the geometric mean over flies, its 90 % interval from `n_boot` bootstrap draws of flies (E79's `fly_bootstrap`), the
    between-fly SD of the log2 ratio and the counts."""
    raise NotImplementedError


def flies_ratios(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING, subset: str = 'all',
                 ratios=RATIOS) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: every ratio of `ratios` under one reading and subset of the recordings: (one row per ratio, per recording
    rows with each response, how it was read and whether the recording entered)."""
    raise NotImplementedError


def dose_response(resp: pd.DataFrame, reading: Tuple[str, bool] = MAIN_READING) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A2: per recording and odor the fold change of the response per log10 unit (the slope of ln response on log10
    concentration over the concentrations above MIN_RESPONSE), the linear slope, the lowest concentration above the
    blank; and the values at MODEL_CONCS (`value_at`)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the model: APL's activity, the indicator, the comparison

def family_apl(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """A4: per wiring and odor, from E90's stage runs (APL on): APL's window mean (a), pre-odor mean (a0), evoked
    activity (a - a0), the PN stage's evoked total (sum over glomeruli of the window minus the pre-odor means) and the
    input's total (sum over the 36 glomeruli of the lab's ORN deltas)."""
    raise NotImplementedError


def indicator_ratio(a1, a2, a0, n: float, form: str):
    """The model's ratio of stimulus 1 over stimulus 2 seen through an indicator of exponent `n`: 'baseline',
    (a1^n - a0^n) / (a2^n - a0^n); 'evoked', ((a1 - a0) / (a2 - a0))^n. Equal to the evoked ratio at n 1."""
    raise NotImplementedError


def model_ratio_set(a1: float, a2: float, a0: float, ns: Sequence[float] = INDICATOR_NS,
                    forms: Sequence[str] = FORMS) -> pd.DataFrame:
    """One wiring's ratios over `ns` and `forms` (columns n, form, ratio)."""
    raise NotImplementedError


def classify(ratios: Sequence[float], lo: float, hi: float) -> str:
    """One wiring against the flies' interval [lo, hi]: 'inside' if any ratio lies in it; 'flatter' if all are nearer 1
    than the interval or on the other side of 1; 'steeper' if all are beyond the interval's far end (from 1); 'mixed'
    otherwise."""
    raise NotImplementedError


def comparison_verdict(classes: Sequence[str], min_wirings: int = VERDICT_MIN_WIRINGS) -> str:
    """'inside', 'flatter' or 'steeper' if at least `min_wirings` of the wirings are; otherwise 'not decided'."""
    raise NotImplementedError


def hull(intervals: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    """The smallest interval holding all of `intervals` (V1: the two readings of paa's concentration)."""
    raise NotImplementedError


def v1_verdict(per_estimate: Dict[str, str]) -> str:
    """V1 from its comparison verdict under each estimate of ms's input: 'inside' if inside under any; 'flatter' or
    'steeper' only if that under all; otherwise 'not decided'."""
    raise NotImplementedError


def alternative_apls() -> pd.DataFrame:
    """A5: the ratios of the across-odor comparisons by other APLs: in proportion to the KCs' response rates in each
    genotype (Yang's criterion of 1 z, and the noise units of the bands beside), in proportion to the PN dendrites'
    totals (the megamat odors; E90's `al_tables`), E37's saved table (window mean and peak of every variant)."""
    raise NotImplementedError


def flies_needed(delta_log2: float, sd_log2: float, z: Tuple[float, float] = POWER_Z) -> int:
    """A6: the number of flies whose geometric mean's 90 % interval excludes a value `delta_log2` (log2) away with 80 %
    power, at a between-fly SD `sd_log2` (log2): ceil(((z_a + z_b) sd / delta)^2)."""
    raise NotImplementedError


def power_table(apl: pd.DataFrame, sd_log2: Dict[str, float]) -> pd.DataFrame:
    """A6: for every pair of panel odors (oriented so that the proportional APL's ratio is at least 1) the family's
    evoked ratio (median and range over the wirings), the proportional APL's (the KCs' response rates with APL output
    blocked (TNT), control flies beside), their log2 difference and the flies needed at each between-fly SD of
    `sd_log2`."""
    raise NotImplementedError


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

def reader_identity(recs=None) -> dict:
    """Part 0's reader identities: E79's per-stimulus calyx means under its three baselines and its blank table."""
    raise NotImplementedError


def identity_passed(which: str = 'reader') -> bool:
    """Whether identity.json records the reader's ('reader') or every ('all') identity as passed."""
    raise NotImplementedError


def part_0(wirings: Sequence[str] = ALL_WIRINGS, runs: bool = True) -> dict:
    """Part 0 (module docstring): identity.json."""
    raise NotImplementedError


def part_a() -> dict:
    """Part A (module docstring): A1 to A6 and a_summary.json."""
    raise NotImplementedError


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
