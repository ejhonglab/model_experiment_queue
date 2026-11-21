#!/usr/bin/env python3
"""E114: a compressive map from the ORN recordings to the model's input rates. Every odor's positive input to an exponent
below 1 (or through a saturating map), on the family's eight wirings at fixed scales with both response rates re-held,
then the family fitted again at the exponent the gate prefers.

The entry is item E114 of experiment_priorities.md (proposed by viki 2026-11-19 13:58 at Tom's request of 13:49 relayed
by hal; planned by a Fable subagent, claude-fable-5-2, spawned from viki 13:50; claimed by viki 2026-11-19 19:04 on hal's
yes). It comes from handle S3 of the input side of the error common to all wirings (E77; outcome (a3)). "The family" is
C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit
(E83's `saved_fit`; S: real, P1, P2, M1; T: P3, P4, M2, M3), read as the mean matrix over the eight wirings (E72's
population, 'all8'), and over S and T beside. "The two ρ's" are ρ on odor-odor distances with APL on against the
control flies and with APL off against APL output blocked (TNT), on raw counts (cosine distance on spike counts), each
KC matrix the mean over all recordings or over one of E108's halves A to D (E108's `dist_targets`); the APL-effect
statistics are E108's `unit_scores` (raw counts and the calibrated readout, against all recordings and the halves).

Concentrations. paa's input is the lab's recording at 10^-4 for the flies' 10^-3; euc and -aPine stay at the recorded
10^-1.75 and 10^-1.5 for the flies' 10^-2; B-myr and 2-but are at the flies' concentrations. Every output involves paa,
euc and -aPine (the map acts on every odor). The holds stand: the model's hemibrain wiring and E72's re-draws of it only;
no other connectome table; no raw-only recording. Nothing is published to claude.ai.

The inputs (glomeruli x the 16 odors, the lab's 36 rows; inputs/<name>.csv with sha256.json; fixed before any run):
- pow<e>: E77's S3 rule (`E77_input_screen.orn_power`): each odor's positive deltas to the power e, rescaled so that
  their sum over glomeruli is the odor's own sum of positive deltas; zero and negative deltas unchanged. Part A:
  e = 0.6, 0.7, 0.8, 0.9 (EXPONENTS); the wrong-way control e = 1.15 (WRONG_WAY); pow1 is the identity (part 0).
- sat<c>: the saturating map x / (1 + x / c) on each positive delta (x in the model's input units, estimated spike
  deltas in Hz), c = 100, 200, 400 Hz (SAT_C), rescaled by the same rule (the odor's sum of positive deltas kept); zero
  and negative deltas unchanged. Part B.
- I0: the lab's input, `H.apl16_orn_deltas()` (the reference: the saved fits).

Parts (CLI `--part`; every simulation part is resumable, a point done when its row is in the part's table):
- inputs (no simulations): the inputs' files, sha256.json and inputs_summary.csv (per input and odor: the largest
  delta, the sum of positive deltas, the participation ratio over glomeruli).
- 0 (simulations): the identity. pow1 equals the lab's input bit for bit; on each wiring (`--wirings` restricts), the
  simulator built with pow1 and its breadth recomputed (E85's `build`, this experiment's cache) gives the wiring's
  breadth and reproduces the saved APL-on and APL-off runs to the last count and response, and the held reading of pow1
  (E85's `held_reading`) reproduces them with both factors exactly 1 (identity.json; 'all_exact' needs all eight).
  If anything differs the experiment stops there.
- 0c (one E10_joint_fit.py smoke run, 5 evaluations): the refit's path. Part C's command on the real wiring with
  `--orn-input results/E114_input_compression/inputs/pow1.csv` and `--init-params
  results/E114_input_compression/starts.json:check_real` (the fit's own absolute parameters; `--init-from-exact` refuses
  a run whose `--orn-input` differs from the fit's row, so the start goes through `--init-params`, as part C's does):
  its first evaluation equals the fit's best history row and its fits.csv row on every shared column (E110's
  `exact_check`; check.json).
- A (simulations): each of pow0.6 to pow0.9 and pow1.15 on the eight wirings at fixed scales with both response rates
  re-held (the breadth recomputed for the input, the somatic threshold bisected to the fit's APL-off response count,
  then the tied APL scale to its APL-on count; E85's `held_reading`). Rows in runs_A.csv (absolute parameters on and
  off, the factors, realism under both response criteria, `H.code_commits()`); runs under outputs/.../A/.
- B (simulations): sat100, sat200, sat400 on the eight wirings, the same way (runs_B.csv; outputs/.../B/).
- score (no simulations): readings.csv (per input, unit and readout, E108's `unit_scores` and the runs' means; 'd_'
  columns the change from I0's same unit and readout), runs.csv, gate.csv and gate.json (`gate` per input of parts A,
  B and C), choice.json (part C's exponent, `choose_c`), outcome.json (`outcome_of`), beside tables when their inputs
  exist (e112.csv, pn_stage.csv, bouton_stage.csv, d1.csv, d2.csv), comparison_table.md, realism_table.md.
- Cstart (simulations only when the chosen exponent has no part-A runs, i.e. 0.75): the held readings of the chosen
  exponent on the eight wirings (runs_Cs.csv; outputs/.../Cs/), then starts.json (C_<wiring>: the held point's APL-on
  parameters, absolute).
- C (fits, one E10_joint_fit.py process each): the family fitted again with the chosen map on each of the eight
  wirings: E72's part-B command (`-k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl --apl-readout both --tag
  _cap3_C --maxfev 150`, the wiring's options as E110's `wiring_args`) with `--orn-input
  results/E114_input_compression/inputs/pow<e>.csv --init-params results/E114_input_compression/starts.json:C_<w>`,
  into results/ and outputs/ of E114_input_compression/fits/C_<w>/ (launch.json beside; logs under outputs/.../logs/).
- D (simulations; only when part C's class is "better", `--force` overrides): D1, the chosen map on the validation2
  ORN input (E107's `panel_input('validation2')`, the same rule per odor), T's four wirings at their refits'
  parameters with the breadth of the mapped megamat input on that wiring (the refit's own), nothing re-held, APL on
  and acutely off; read against Remy's validation2 somata as E109's `remy_readings` (E107's carrying set of pairs and
  the 31 of paa and B-myr, her copy under E107's lock), with E107's saved runs of the family as fitted beside
  (E109's `e107_reference`). D2, the map on top of E113's best input when E113 has one (`--e113-input <csv>`, glomeruli
  x abbreviated odors as E10's `--orn-input` reads it; without it D2 records that E113 has none): E113's input and the
  map on it, each on the eight wirings held as in part A.
- beside (simulations, no class): the PN stage of every input (E90's `al_stage` at the bouton base's parameters,
  `pn_statistics` against the PN dendrites, E90's al_transfer.csv); the bouton stage on the real wiring, APL on, at each
  input's held point (I0 at the fit, part C's input at its refit): E90's `read_stages` over its MODEL_WINDOW, the
  boutons of the 28 glomeruli the flies share, 1 - Pearson across boutons per megamat pair, group means over E90's
  GROUPS for the 'rate' and 'session_z' readings, beside the flies' boutons as measured and with the noise removed
  (E90's stages_flies.csv). These need no KC data.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive and kcq/apl_effect.py.

The gate (the entry's, fixed before any run; `gate`; each input's population of eight against I0's, each wiring paired
with itself; E108's constants):
- Gain: one ρ higher by GAIN (0.010) or more for the population against all recordings, on at least 6 of the 8 wirings
  (against all recordings) and for the population against each of halves A to D; the other ρ not lower by more than
  OTHER_SLACK (0.005) for the population against all recordings.
- The wrong-way control loses: on each ρ that gains, pow1.15's population change against all recordings is below zero.
  The same part-A control serves parts B and C (it is never refitted). A failure is named among the fails.
- The APL effect: the pattern part under the calibrated readout not higher by more than 0.003 against all recordings
  and each of halves A to D; the raw-count loss (all recordings) not higher by more than 0.005; the mean APL effect
  (calibrated) not lower by more than 0.01. The amplitude part, the shape part, the spread ratio, ρ and Pearson r on the
  per-pair APL effect go beside.
- Realism: every one of the input's eight runs inside the four bands (E83's `realism`, the bands of 2026-11-10; the
  flags under Yang's z rule beside).
- Response rates: response rate ρ with APL on and with APL off (the bands' units, the mean over the eight runs) not
  lower by more than 0.05.
- The block, the natural-source alignment with APL on, off and on minus off, and the structure against the wild-type
  somata are columns only: they enter no class (the entry: "The block is a column only").
- Classes (`class_of`): "better" (the gain, the control and every constraint hold); "a price" (the gain holds and a
  constraint or the control fails; each failure named with its size); "an equal" (no gain, every constraint holds and
  neither ρ lower than I0's by more than OTHER_SLACK); "no gain" (otherwise).

Choices the entry left open (fixed here before any run, with the reason):
- "Loses" for the wrong-way control: its population change below zero on each ρ the candidate gains on (the entry's
  word; E77's screen had 1.33 lower on both ρ's, so a sign test is the plainest reading of 1.15).
- Part C's exponent ("the largest gated gain"): among pow0.6 to pow0.9, those whose class is "better"; of them the one
  with the largest population change against all recordings on a ρ that gains; on a tie the exponent nearer 1. None
  "better": 0.75, as the entry says ("if none passes"), whose held points are run first (part Cstart) as the refit's
  start and enter no choice; they are gated like part A's inputs once run (`gate_tables`), for outcome (b)'s "at fixed
  scales" (added 2026-11-20 after the chain: with the fallback 0.75 not among part A's exponents, (b) had no
  fixed-scale gate to read and the letter fell to '-'). The saturating form is not refitted (the entry refits "the exponent").
- The refit's start: each wiring's held point of the chosen exponent (the fit's own parameters with the somatic
  threshold and the tied APL scale re-held for the new input; E110 part B's precedent), one start, 150 evaluations
  (E72's), KC-KC free. The bounds are E10's (relative to E3's base parameters, which do not depend on the input).
- Part C's gate: the refits' population against the family's saved fits (I0), with part A's pow1.15 as its control.
- Part D runs on part C's class "better" only ("only on a pass after the refit"); D1's wirings are T's at their refits;
  D2's reference is E113's input held the same way, and its pass is not gated (beside).
- Outcome letters (`outcome_of`, after part C; a provisional letter from part A before it):
  (a) part C's class is "better";
  (b) part C's class is "a price" with the calibrated pattern part among its fails, and part A's chosen exponent's also
      (the price at fixed scales and after the refit);
  (c) no exponent of part A has the population gain, but on the real wiring alone one ρ is higher by GAIN or more with
      the other not lower by more than OTHER_SLACK (the hemibrain wiring's own point, E77's);
  (d) no exponent of part A has the population gain, and the saturating form differs from the exponent: a sat<c> has
      the population gain, or the best population change of a ρ over the sat<c> differs from that over the exponents
      by GAIN or more;
  '-' none of these (reported with the classes).
- The reading E112 names: E112 had not ended when this was written. score gives, beside, the two ρ's of every input's
  population under E112's readings R0 to R3 (`E112_rho_like_for_like.run_readings` / `reading_scores`) with E112's
  index row state, and names none as carrying; that is E112's.
- E108's halves A to D are the halves (E77's `recording_halves(0)` and `(1)`).

Departures from the entry: none so far.

Provenance: the other checkouts must be clean; this repository clean in every file a process loads (E108's
`_commits_clean`; for the fits E110's `_fit_clean`, nothing modified under experiments/ or kcq/); every row records
`H.code_commits()`; every simulation part's provenance has `H.model_provenance(orn_deltas=...)` per point.

Outputs: <shared>/results/E114_input_compression/; spike counts and responses under outputs/E114_input_compression/
(A, B, Cs, D, beside, fits/), then data/model_outputs/E114_input_compression/.
Part stages (added after the chain, for E120's part A4): beside's bouton stage again with its bouton windows saved
(outputs/E114_input_compression/stages/), checked against bouton_stage.csv (stage_check.csv).

Usage (from the repository root; one 16-thread lane on viki, under systemd-run --user; outputs/E114_input_compression/
chain.sh runs the parts in order):
    ./venv/bin/python -u experiments/E114_input_compression.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part 0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part 0c
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part B
    ./venv/bin/python -u experiments/E114_input_compression.py --part score
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part Cstart
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part C
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part D
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E114_input_compression.py --part beside
    ./venv/bin/python -u experiments/E114_input_compression.py --part score
    ./venv/bin/python -u experiments/E114_input_compression.py --part variants
"""
from __future__ import annotations

import argparse
import gc
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
from experiments import E83_cross_wiring_fit as E83
from experiments import E108_paa_extrapolation as E108

EXPERIMENT = 'E114_input_compression'
SCRIPT = 'experiments/E114_input_compression.py'
EXPONENTS = (0.6, 0.7, 0.8, 0.9)
WRONG_WAY = 1.15
FALLBACK_EXPONENT = 0.75
SAT_C = (100.0, 200.0, 400.0)
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
UNITS = (('all8', ALL_WIRINGS), ('S', S_WIRINGS), ('T', T_WIRINGS))
READOUTS = ('raw', 'cal')
HALVES = E108.HALVES
TARGETS = E108.TARGETS
# the gate (the entry's; E108's constants)
GAIN, GAIN_WIRINGS, OTHER_SLACK = E108.GAIN, E108.GAIN_WIRINGS, E108.OTHER_SLACK
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = E108.EQUAL_PATTERN, E108.EQUAL_RAW, E108.EQUAL_MEAN
RATE_RHO_SLACK = E108.RATE_RHO_SLACK
CLASSES = ('better', 'a price', 'an equal', 'no gain')
# part C (E72's part-B command; E110's helpers)
FIT_PREFIX = f'{EXPERIMENT}/fits'
CHECK_PREFIX = f'{EXPERIMENT}/check'
MAXFEV = 150
SMOKE_TAG = '_cap3_check0'
FIT_SUFFIX = 'C'
REFIT = 'refit'                                    # the input name of part C's runs in the readings
E113_REF, E113_MAPPED = 'e113', 'e113_map'         # part D2's input names

_CACHE: Dict[object, object] = {}


# ---------------------------------------------------------------------------
# folders and small helpers

def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(f'{EXPERIMENT}/{part}' if part else EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M:%S %Z')


def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


# ---------------------------------------------------------------------------
# the maps and the inputs

def power_map(orn: pd.DataFrame, e: float) -> pd.DataFrame:
    """pow<e>: E77's S3 rule (`E77_input_screen.orn_power`): each odor's positive deltas to the power `e`, rescaled so
    that their sum over glomeruli is the odor's own sum of positive deltas; zero and negative deltas unchanged. e = 1
    gives the input back bit for bit."""
    from experiments.E77_input_screen import orn_power
    return orn_power(orn, float(e))


def saturating_map(orn: pd.DataFrame, c: float) -> pd.DataFrame:
    """sat<c>: x / (1 + x / c) on each positive delta x (the model's input units), rescaled as the S3 rule does (the
    odor's sum of positive deltas kept); zero and negative deltas unchanged."""
    out = orn.astype(float).copy()
    x = out.to_numpy()
    pos = x > 0
    xp = np.where(pos, x, 0.0)
    y = np.where(pos, xp / (1.0 + xp / float(c)), x)
    s0 = xp.sum(axis=0)
    s1 = np.where(pos, y, 0.0).sum(axis=0)
    scale = np.divide(s0, s1, out=np.ones_like(s0), where=s1 > 0)
    out[:] = np.where(pos, y * scale[None, :], x)
    return out


def map_of(name: str):
    """The map of an input name ('pow<e>' or 'sat<c>'): a function of an ORN frame."""
    if name.startswith('pow'):
        e = float(name[3:])
        return lambda orn: power_map(orn, e)
    if name.startswith('sat'):
        c = float(name[3:])
        return lambda orn: saturating_map(orn, c)
    raise ValueError(f'no map {name!r}')


def pow_name(e: float) -> str:
    return f'pow{float(e):g}'


def sat_name(c: float) -> str:
    return f'sat{float(c):g}'


A_INPUTS = tuple(pow_name(e) for e in EXPONENTS) + (pow_name(WRONG_WAY),)
B_INPUTS = tuple(sat_name(c) for c in SAT_C)
CONTROL = pow_name(WRONG_WAY)
CANDIDATES = tuple(pow_name(e) for e in EXPONENTS) + B_INPUTS


def lab_input() -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()`), as E108 reads it."""
    return E108.lab_input()


def input_of(name: str) -> pd.DataFrame:
    """An input by name: 'I0' (the lab's), 'pow<e>' or 'sat<c>' (the map of the lab's input); D2's 'e113' and
    'e113_map' come from `e113_input`."""
    c = _CACHE.setdefault('inputs', {})
    if name not in c:
        if name == 'I0':
            c[name] = lab_input()
        elif name.startswith(('pow', 'sat')):
            c[name] = map_of(name)(lab_input())
        elif name in (E113_REF, E113_MAPPED):
            raise RuntimeError(f'{name}: set by part D2 (`e113_input`)')
        else:
            raise ValueError(f'unknown input {name!r}')
    return c[name].copy()


def input_path(name: str) -> Path:
    return out_dir() / 'inputs' / f'{name}.csv'


def input_spec(name: str) -> str:
    """The `--orn-input` of E10_joint_fit.py (relative to the shared directory)."""
    return str(input_path(name).relative_to(H.shared_dir()))


def write_input(name: str, orn: Optional[pd.DataFrame] = None) -> dict:
    """inputs/<name>.csv (glomeruli x abbreviated odors) and its record in inputs/sha256.json."""
    from experiments.E77_input_screen import frame_sha256
    orn = input_of(name) if orn is None else orn
    lay = H.abbrev_columns(orn)
    lay.index.name = 'glomerulus'
    p = input_path(name)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(lay.to_csv())
    rec = dict(file=str(p.relative_to(H.shared_dir())), file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
               frame_sha256=frame_sha256(orn))
    shas = _read_json(p.parent / 'sha256.json')
    shas[name] = rec
    _write_json(p.parent / 'sha256.json', shas)
    return rec


def input_record(name: str) -> dict:
    """The record part inputs wrote, checked against the input built now (its frame's sha256; within 1e-9 of the file
    otherwise, E108's rule)."""
    from experiments.E77_input_screen import frame_sha256
    rec = _read_json(input_path(name).parent / 'sha256.json').get(name)
    if rec is None:
        raise RuntimeError(f'run --part inputs first (no record of input {name})')
    built = _CACHE.get('inputs', {}).get(name)
    built = input_of(name) if built is None else built
    sha = frame_sha256(built)
    if rec['frame_sha256'] == sha:
        return rec
    diff = E108.recorded_input_diff(built, H.shared_dir() / rec['file'])
    if not diff <= E108.INPUT_ATOL:
        raise RuntimeError(f'input {name} differs from the one part inputs wrote (largest difference {diff:.3g})')
    return dict(rec, frame_sha256=sha, recorded_frame_sha256=rec['frame_sha256'], frame_max_abs_diff=diff)


def participation_ratio(v: np.ndarray) -> float:
    """(sum v)^2 / sum v^2 over the positive entries (the number of glomeruli an odor's input is spread over)."""
    v = np.clip(np.asarray(v, dtype=float), 0.0, None)
    s2 = float((v ** 2).sum())
    return float(v.sum() ** 2 / s2) if s2 > 0 else np.nan


def inputs_summary(names: Sequence[str]) -> pd.DataFrame:
    """Per input and odor: the largest delta, the sum of positive deltas, the participation ratio over glomeruli."""
    rows = []
    for n in names:
        x = H.abbrev_columns(input_of(n))
        for o in x.columns:
            v = x[o].to_numpy(dtype=float)
            rows.append(dict(input=n, odor=o, max=float(v.max()), sum_pos=float(v[v > 0].sum()), pr=participation_ratio(v)))
    return pd.DataFrame(rows)


def write_inputs() -> dict:
    """Part inputs (no simulations): pow1 (part 0c's), part A's and part B's inputs and pow0.75 (part C's fallback),
    inputs_summary.csv."""
    names = (pow_name(1.0),) + A_INPUTS + B_INPUTS + (pow_name(FALLBACK_EXPONENT),)
    recs = {n: write_input(n) for n in names}
    # the round trip through E10_joint_fit.py's reader: what a fit reads is the input built here
    from experiments.E10_joint_fit import read_orn_input
    lab = H.apl16_orn_deltas()
    for n in names:
        back = read_orn_input(input_spec(n), lab)
        d = float(np.max(np.abs(back.to_numpy(dtype=float) - input_of(n).to_numpy(dtype=float))))
        recs[n]['e10_round_trip_max_abs_diff'] = d
        if not d <= E108.INPUT_ATOL:
            raise RuntimeError(f'{n}: E10_joint_fit.read_orn_input differs from the input by {d:.3g}')
    inputs_summary(('I0',) + names).to_csv(out_dir() / 'inputs_summary.csv', index=False)
    _write_json(out_dir() / 'inputs.json', dict(inputs=recs, **H.code_commits()))
    print(json.dumps({n: r['frame_sha256'][:12] for n, r in recs.items()}), flush=True)
    return recs


# ---------------------------------------------------------------------------
# part 0: the identity

def cache_tag(orn: pd.DataFrame, name: str) -> str:
    return E83.input_tag_for(orn, f'E114_{name}')


def identity(units: Optional[Sequence[str]] = None) -> dict:
    """Part 0 (the module docstring); resumable per wiring; identity.json; raises if anything differs."""
    from experiments import E85_orn_fly_population as E85
    path = out_dir() / 'identity.json'
    res = _read_json(path)
    res.update(E108._commits_clean())
    lab = lab_input()
    p1 = power_map(lab, 1.0)
    res['pow1_equals_lab'] = bool(p1.index.equals(lab.index) and p1.columns.equals(lab.columns)
                                  and np.array_equal(p1.to_numpy(), lab.to_numpy()))
    if not res['pow1_equals_lab']:
        _write_json(path, res)
        raise RuntimeError('E114 part 0: pow1 is not the lab input bit for bit')
    res.setdefault('units', {})
    for u in (units or ALL_WIRINGS):
        if u in res['units'] and 'held_thr_factor' in res['units'][u]:
            continue
        t0 = time.time()
        F = E83.saved_fit(u)
        tag = cache_tag(p1, 'pow1')
        sim, info = E85.build(u, p1, tag, cache_dir=local_dir())
        own = E83.wiring_parts(u)['breadth']
        b = info['breadth']
        beq = bool(b.index.equals(own.index) and np.array_equal(b.to_numpy(dtype=float), own.to_numpy(dtype=float)))
        on = sim(dict(F['p']))
        off = sim(dict(F['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
        ns, nr = E83._same({'on': on, 'off': off}, F['saved'])
        del sim, info, on, off
        gc.collect()
        hr = E85.held_reading(u, p1, tag, cache_dir=local_dir())
        hs, hn = E83._same({'on': hr['on'], 'off': hr['off']}, F['saved'])
        res['units'][u] = dict(fit=f"{F['folder']}:{F['tag']}", breadth_equal=beq, count_diffs=ns, response_diffs=nr,
                               held_count_diffs=hs, held_response_diffs=hn, held_thr_factor=float(hr['thr_factor']),
                               held_apl_factor=float(hr['apl_factor']), held_n_sims=int(hr['n_sims']),
                               seconds=round(time.time() - t0, 1), host=_host(),
                               threads=os.environ.get('OMP_NUM_THREADS'))
        _write_json(path, res)
        print(f'E114 0 {u}: {res["units"][u]}', flush=True)
        del hr
        gc.collect()
    ok = lambda d: (d.get('breadth_equal') and d['count_diffs'] == 0 and d['response_diffs'] == 0
                    and d['held_count_diffs'] == 0 and d['held_response_diffs'] == 0
                    and d['held_thr_factor'] == 1.0 and d['held_apl_factor'] == 1.0)
    done = [res['units'][u] for u in ALL_WIRINGS if u in res['units']]
    res['all_exact'] = bool(len(done) == len(ALL_WIRINGS) and all(ok(d) for d in done))
    res['checked'] = _now()
    _write_json(path, res)
    bad = {u: d for u, d in res['units'].items() if not ok(d)}
    if bad:
        raise RuntimeError(f'E114 part 0: pow1 does not reproduce the saved fits: {bad} (the experiment stops here)')
    return res


def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('all_exact'))


def fit_base_command(name: str, prefix: str = FIT_PREFIX, suffix: str = FIT_SUFFIX) -> List[str]:
    """E72's part-B command (E110's `fit_command` without a start) into results/ and outputs/ of <prefix>/<name>."""
    py = str(H.HERE / 'venv' / 'bin' / 'python')
    return [py, '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3', '--rate-target', 'ctrl_frac', '--loss',
            'apl', '--apl-readout', 'both', '--tag', f'_cap3_{suffix}', '--maxfev', str(MAXFEV), '--experiment',
            f'{prefix}/{name}']


def check_command() -> List[str]:
    """Part 0c's command: pow1 on the real wiring from the fit's own parameters (starts.json:check_real; `--smoke`)."""
    return (fit_base_command('exact_real_pow1', CHECK_PREFIX, 'check0') + ['--orn-input', input_spec(pow_name(1.0)),
            '--init-params', f'results/{EXPERIMENT}/starts.json:check_real', '--smoke'])


def write_check_start() -> dict:
    """starts.json:check_real, the real wiring's fit's absolute parameters (part 0c's start)."""
    from experiments import E110_no_kckc as E110
    F = E83.saved_fit('real')
    e = dict({k: float(F['p'][k]) for k in E83.PARAM_KEYS}, source=f"the fit {F['folder']}:{F['tag']}")
    return E110.merge_starts(out_dir() / 'starts.json', {'check_real': e})


def _fits_row(res_dir: Path) -> Optional[pd.Series]:
    """The last row of a fit folder's fits.csv (one fit per folder), or None."""
    p = res_dir / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    return f.iloc[-1] if len(f) else None


def launch(cmd: List[str], res: Path, log: Path, cc: dict, name: str) -> dict:
    """Runs one E10_joint_fit.py process (E110's `launch`, with this script as the launcher)."""
    from experiments import E110_no_kckc as E110
    rec = dict(name=name, command=cmd, cwd=str(H.HERE), host=_host(), started=_now(), commits=cc,
               queue_dirty_paths=str(cc.get('queue_dirty_paths') or '').split(), log=str(log.relative_to(H.HERE)),
               threads=dict(OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS'),
                            OMP_THREAD_LIMIT=E110._child_env().get('OMP_THREAD_LIMIT')), launcher=SCRIPT)
    res.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    _write_json(res / 'launch.json', rec)
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT, env=E110._child_env())
    rec.update(finished=_now(), returncode=r.returncode, seconds=time.time() - t0)
    _write_json(res / 'launch.json', rec)
    return rec


def fit_path_check(rerun: bool = False) -> dict:
    """Part 0c (the module docstring); check.json; raises if the first evaluation differs from the fit."""
    from experiments import E110_no_kckc as E110
    if not identity_passed():
        raise RuntimeError('run --part 0 first')
    cmd = check_command()
    res_dir = H.shared_dir() / 'results' / CHECK_PREFIX / 'exact_real_pow1'
    rec = _read_json(res_dir / 'launch.json')
    hist = sorted(res_dir.glob('history_*.csv'))
    if rerun or not hist or rec.get('returncode') is None or not E110.same_command(rec.get('command') or [], cmd):
        cc = E110._fit_clean()
        write_check_start()
        print('E114 0c: ' + ' '.join(cmd), flush=True)
        rec = launch(cmd, res_dir, local_dir('logs') / 'check_exact_real_pow1.log', cc, 'exact_real_pow1')
        hist = sorted(res_dir.glob('history_*.csv'))
    out = dict(command=' '.join(cmd), launch={k: rec.get(k) for k in ('host', 'started', 'finished', 'returncode', 'seconds',
                                                                       'log', 'threads')}, commits=rec.get('commits'))
    if rec.get('returncode') != 0 or not hist:
        out['passes'] = False
    else:
        h = pd.read_csv(hist[-1])
        folder, ftag = E83.members()['real']
        fh = pd.read_csv(H.results_dir(folder) / f'history_{ftag}.csv')
        best = fh.loc[fh['loss'].idxmin()]
        row = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[ftag]
        ex = E110.exact_check(h.iloc[0], best, row)
        out.update(history=hist[-1].name, n_evaluations=int(len(h)), reference=f'{folder}:{ftag} (the best history row; '
                   'fits.csv row)', **ex)
        out['passes'] = bool(ex['exact'])
    out['checked'] = _now()
    _write_json(out_dir() / 'check.json', out)
    print(json.dumps({k: out.get(k) for k in ('passes', 'n_evaluations', 'history_differs', 'row_differs')}, default=str),
          flush=True)
    if not out['passes']:
        raise RuntimeError('E114 part 0c fails (check.json); part C does not run until it passes')
    return out


def check_passed() -> bool:
    return bool(_read_json(out_dir() / 'check.json').get('passes'))


# ---------------------------------------------------------------------------
# parts A, B and Cstart: held readings at fixed scales

def run_tag(part: str, name: str, unit: str) -> str:
    return f'{part}_{name}_{unit}'


def plan(part: str) -> List[Tuple[str, str]]:
    """[(input, wiring)] of a simulation part, in order."""
    if part == 'A':
        return [(n, w) for n in A_INPUTS for w in ALL_WIRINGS]
    if part == 'B':
        return [(n, w) for n in B_INPUTS for w in ALL_WIRINGS]
    if part == 'Cs':
        ch = _read_json(out_dir() / 'choice.json')
        if not ch.get('exponent'):
            raise RuntimeError('no choice.json: run --part score after parts A and B')
        n = pow_name(ch['exponent'])
        return [] if n in A_INPUTS else [(n, w) for w in ALL_WIRINGS]
    if part == 'D2':
        return [(n, w) for n in (E113_REF, E113_MAPPED) for w in ALL_WIRINGS]
    raise ValueError(part)


def runs_of(part: str) -> pd.DataFrame:
    p = out_dir() / f'runs_{part}.csv'
    return E108._read_rows(p) if p.exists() else pd.DataFrame(columns=['tag', 'input', 'unit'])


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run: this host's outputs, else the model-output archive; or None."""
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_row(part: str, name: str, unit: str, tag: str, res: dict, o: Dict[str, tuple], rec: dict, cc: dict) -> dict:
    F = E83.saved_fit(unit)
    row = dict(part=part, input=name, unit=unit, tag=tag, fit_folder=F['folder'], fit_tag=F['tag'],
               input_file_sha256=rec.get('file_sha256'), input_frame_sha256=rec['frame_sha256'], input_tag=res['input_tag'],
               orn_input=rec.get('file'), thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'],
               n_off=res['n_off'], target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
               at_bound_apl=res['at_bound_apl'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: res['params'][k] for k in E83.PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}),
               K_run=E83.readout_K(o), breadth_sha256=res.get('breadth_sha256'))
    row.update(E83.realism(o))
    row.update(cc)
    return row


def part_sims(part: str, wirings: Optional[Sequence[str]] = None, max_points: Optional[int] = None) -> pd.DataFrame:
    """Part A, B, Cs (Cstart's held readings) or D2: the plan's points in order, resumable (a point is done when its
    row is in runs_<part>.csv): E85's `held_reading` (both response rates held), runs saved under outputs/.../<part>/."""
    from experiments import E85_orn_fly_population as E85
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identity must hold before any reading)')
    cc = E108._commits_clean()
    path = out_dir() / f'runs_{part}.csv'
    prov_path = out_dir() / f'provenance_{part}.json'
    n = 0
    for name, unit in plan(part):
        if wirings is not None and unit not in wirings:
            continue
        tag = run_tag(part, name, unit)
        if tag in set(runs_of(part)['tag']) or (max_points is not None and n >= max_points):
            continue
        E108._still_clean(cc)
        orn = input_of(name) if part != 'D2' else _CACHE['inputs'][name]
        rec = input_record(name)
        res = E85.held_reading(unit, orn, cache_tag(orn, name), cache_dir=local_dir())
        o = {'on': res['on'], 'off': res['off']}
        save_run(part, tag, o)
        row = run_row(part, name, unit, tag, res, o, rec, cc)
        E83._append_row(path, row)
        prov = _read_json(prov_path) or dict(part=part, script=SCRIPT, points={})
        prov['points'][tag] = H.model_provenance(orn_deltas=orn)
        _write_json(prov_path, prov)
        print(f"E114 {part} {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/{res['target_off']} off); flag '{row.get('realism_flag') or ''}'; {res['seconds']:.0f}s",
              flush=True)
        n += 1
        del res, o
        gc.collect()
    return runs_of(part)


# ---------------------------------------------------------------------------
# part C: the refits

def chosen_name() -> str:
    ch = _read_json(out_dir() / 'choice.json')
    if not ch.get('exponent'):
        raise RuntimeError('no choice.json: run --part score after parts A and B')
    return pow_name(ch['exponent'])


def held_runs_of(name: str) -> pd.DataFrame:
    """The held runs of an input (part A's, or Cstart's)."""
    r = pd.concat([runs_of('A'), runs_of('Cs')], ignore_index=True, sort=False)
    return r[r['input'] == name] if len(r) else r


def write_starts() -> dict:
    """starts.json: C_<w>, the chosen exponent's held point of each wiring (its APL-on parameters, absolute); a start is
    fixed once written (E110's `merge_starts`)."""
    from experiments import E110_no_kckc as E110
    name = chosen_name()
    entries = {}
    for _, r in held_runs_of(name).iterrows():
        p = json.loads(r['params_on'])
        entries[f"C_{r['unit']}"] = dict({k: float(p[k]) for k in E83.PARAM_KEYS}, input=name,
                                         source=f"the held point {r['tag']} (results/{EXPERIMENT}/runs_*.csv)")
    if len(entries) != len(ALL_WIRINGS):
        raise RuntimeError(f'{name}: held points on {sorted(entries)} only; run --part Cstart')
    return E110.merge_starts(out_dir() / 'starts.json', entries)


def fit_command(w: str, name: Optional[str] = None) -> List[str]:
    """Part C's command on wiring `w` (the module docstring)."""
    from experiments import E110_no_kckc as E110
    name = chosen_name() if name is None else name
    return (fit_base_command(f'C_{w}') + E110.wiring_args(w) + ['--orn-input', input_spec(name), '--init-params',
            f'results/{EXPERIMENT}/starts.json:C_{w}'])


def fit_dirs(w: str) -> Tuple[Path, Path]:
    """(results, outputs) folders of part C's fit on `w`."""
    return H.shared_dir() / 'results' / FIT_PREFIX / f'C_{w}', H.HERE / 'outputs' / FIT_PREFIX / f'C_{w}'


def load_fit(w: str) -> Optional[Tuple[pd.Series, Dict[str, tuple]]]:
    """(fits.csv row, saved runs) of part C's fit on `w`: this host's outputs, else the archive's; or None."""
    rd, od = fit_dirs(w)
    row = _fits_row(rd)
    if row is None:
        return None
    for d in (od, H.shared_dir() / 'data' / 'model_outputs' / FIT_PREFIX / f'C_{w}'):
        o = E83._load_pair(d, str(row['tag']))
        if o is not None:
            return row, o
    return None


def part_c(wirings: Optional[Sequence[str]] = None) -> List[str]:
    """Part C: the refits, one E10_joint_fit.py process after another (skipping those done); runs_C.csv."""
    from experiments import E110_no_kckc as E110
    if not check_passed():
        raise RuntimeError('run --part 0c first')
    write_starts()
    cc = E110._fit_clean()
    done = []
    for w in (wirings or ALL_WIRINGS):
        rd, _ = fit_dirs(w)
        if _fits_row(rd) is not None:
            done.append(w)
            continue
        E108._still_clean(cc)
        cmd = fit_command(w)
        print(f'E114 C {w}: ' + ' '.join(cmd), flush=True)
        rec = launch(cmd, rd, local_dir('logs') / f'C_{w}.log', cc, f'C_{w}')
        if rec['returncode']:
            raise RuntimeError(f"E114 C {w}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
        r = _fits_row(rd)
        if r is None:
            raise RuntimeError(f"E114 C {w}: no row in its fits.csv (log {rec['log']})")
        E83._append_row(out_dir() / 'runs_C.csv', dict(name=f'C_{w}', unit=w, input=chosen_name(), tag=r['tag'],
                        folder=f'{FIT_PREFIX}/C_{w}', host=_host(), seconds=rec['seconds'], command=' '.join(cmd),
                        **E110.params_of_row(r), loss=float(r['loss']), n_evals=int(r['n_evals']), **cc))
        print(f"E114 C {w}: loss {float(r['loss']):.4f}, {rec['seconds'] / 60:.0f} min", flush=True)
        done.append(w)
    return done


# ---------------------------------------------------------------------------
# the readings

def run_outputs(name: str, unit: str) -> Dict[str, tuple]:
    """A run of the readings: I0 the unit's saved fit; REFIT part C's fit; a part's held run otherwise."""
    if name == 'I0':
        return E83.saved_fit(unit)['saved']
    if name == REFIT:
        f = load_fit(unit)
        if f is None:
            raise RuntimeError(f'no refit on {unit}')
        return f[1]
    for part in ('A', 'B', 'Cs', 'D2'):
        o = load_run(part, run_tag(part, name, unit))
        if o is not None:
            return o
    raise RuntimeError(f'no saved run of {name} on {unit}')


def _dist(name: str, unit: str, readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83's `pair_distances`), cached by the sha256 of its counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    mem = _CACHE.setdefault('dist_mem', {})
    key = (name, unit, readout)
    if key in mem:
        return mem[key]
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
    o = run_outputs(name, unit)
    skey = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if skey not in _CACHE['dist_store']:
        _CACHE['dist_store'][skey] = E83.pair_distances(o, readout)
        _CACHE['dist_dirty'] = True
    mem[key] = _CACHE['dist_store'][skey]
    return mem[key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist_store']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


def runs_table() -> pd.DataFrame:
    """Every run: I0 (the saved fits), parts A, B, Cs and D2 (held) and part C (refits), with 'flagged'."""
    rows = []
    for u in ALL_WIRINGS:
        F = E83.saved_fit(u)
        rows.append(dict(part='I0', tag=f'I0_{u}', input='I0', unit=u, thr_factor=1.0, apl_factor=1.0, fit_folder=F['folder'],
                         fit_tag=F['tag'], **{k: F['p'][k] for k in E83.PARAM_KEYS}, **E83.realism(F['saved'])))
    parts = [runs_of(p) for p in ('A', 'B', 'Cs', 'D2')]
    for u in ALL_WIRINGS:
        f = load_fit(u)
        if f is None:
            continue
        from experiments import E110_no_kckc as E110
        r, o = f
        rows.append(dict(part='C', tag=f'C_{u}', input=REFIT, unit=u, fit_folder=f'{FIT_PREFIX}/C_{u}', fit_tag=r['tag'],
                         loss=float(r['loss']), n_evals=int(r['n_evals']), **E110.params_of_row(r), **E83.realism(o)))
    df = pd.concat([pd.DataFrame(rows)] + [p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates(['input', 'unit'], keep='last').reset_index(drop=True)
    df['flagged'] = [E108._flagged(x) for x in df['realism_flag']]
    return df


RUN_MEANS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'rate_rho_on_zcrit', 'rate_rho_off_zcrit',
             'opk_on_tv', 'opk_off_tv', 'max_count_on', 'max_count_off', 'thr_factor', 'apl_factor')


def readings() -> pd.DataFrame:
    """readings.csv (per input, unit and readout: E108's `unit_scores` and the runs' means; 'd_' columns the change from
    I0's same unit and readout) and runs.csv."""
    out = out_dir()
    runs = runs_table()
    runs.to_csv(out / 'runs.csv', index=False)
    rows = []
    for name, g in runs.groupby('input', sort=False):
        have = set(g['unit'])
        units = [(u, (u,)) for u in g['unit']] + [(u, us) for u, us in UNITS if all(x in have for x in us)]
        for ro in READOUTS:
            dd = {u: _dist(name, u, ro) for u in g['unit']}
            for unit, us in units:
                single = len(us) == 1 and unit in dd
                d = dd[unit] if single else E83.population([dd[u] for u in us])
                sub = g[g['unit'].isin(us)]
                rows.append(dict(input=name, unit=unit, kind='wiring' if single else 'population', n_wirings=len(us),
                                 readout=ro, n_flagged=int(sub['flagged'].sum()),
                                 **{f'run_{k}': float(sub[k].astype(float).mean()) for k in RUN_MEANS if k in sub},
                                 **E108.unit_scores(d)))
    _save_dist_cache()
    rd = pd.DataFrame(rows)
    rd['vcf_diff'] = rd['vcf_on'] - rd['vcf_off']
    base = rd[rd['input'] == 'I0'].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c not in ('input', 'unit', 'kind', 'n_wirings', 'readout')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, r), np.nan) for u, r, v in zip(rd['unit'], rd['readout'], rd[c])]
    rd.to_csv(out / 'readings.csv', index=False)
    return rd


# ---------------------------------------------------------------------------
# the gate, the choice and the outcome

def _v(rd: pd.DataFrame, name: str, unit: str, col: str, readout: str = 'raw') -> float:
    r = rd[(rd['input'] == name) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _have(runs: pd.DataFrame, name: str, units: Sequence[str] = ALL_WIRINGS) -> bool:
    got = set(zip(runs['input'], runs['unit']))
    return all((name, u) in got for u in units)


def gate(rd: pd.DataFrame, runs: pd.DataFrame, name: str, ref: str = 'I0', control: Optional[str] = CONTROL) -> dict:
    """The entry's gate (the module docstring) for input `name` against `ref` on the population of the eight wirings.
    `rd`: `readings`' layout; `runs`: `runs_table`'s ('input', 'unit', 'flagged', 'rate_rho_<on|off>'). `control`: the
    wrong-way control's input (None: no control condition)."""
    out = dict(input=name, reference=ref, control=control, complete=bool(_have(runs, name) and _have(runs, ref)))
    if not out['complete']:
        out.update(gain=False, note='runs missing')
        return out
    D = lambda col, ro='raw', u='all8', n=name: _v(rd, n, u, col, ro) - _v(rd, ref, u, col, ro)
    other = {'on': 'off', 'off': 'on'}
    gains = {}
    for c in ('on', 'off'):
        d_all = D(f'rho_{c}_all')
        n_w = int(sum(D(f'rho_{c}_all', 'raw', w) >= GAIN for w in ALL_WIRINGS))
        d_h = {h: D(f'rho_{c}_{h}') for h in HALVES}
        gains[c] = bool(d_all >= GAIN and n_w >= GAIN_WIRINGS and all(v >= GAIN for v in d_h.values()))
        out.update({f'd_rho_{c}_all': d_all, f'n_wirings_gain_{c}': n_w, **{f'd_rho_{c}_{h}': v for h, v in d_h.items()},
                    **{f'd_rho_{c}_{u}': D(f'rho_{c}_all', 'raw', u) for u in ('S', 'T')}})
    for c in ('on', 'off'):
        out[f'gain_{c}'] = bool(gains[c] and out[f'd_rho_{other[c]}_all'] >= -OTHER_SLACK)
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    fails = {}
    # the wrong-way control
    if control is not None:
        ctl = {}
        for c in ('on', 'off'):
            if out[f'gain_{c}']:
                v = D(f'rho_{c}_all', n=control) if _have(runs, control) else np.nan
                ctl[c] = v
                if not (np.isfinite(v) and v < 0):
                    fails[f'wrong-way control ({control}), APL {c}'] = v
        out['control_d_rho'] = {c: D(f'rho_{c}_all', n=control) if _have(runs, control) else np.nan for c in ('on', 'off')}
        out['control_ok'] = (not any(k.startswith('wrong-way') for k in fails)) if ctl else None
    # the APL effect, under the calibrated readout; the raw-count loss
    for t in TARGETS:
        out[f'd_pattern_{t}'] = D(f'{t}_pattern', 'cal')
        if out[f'd_pattern_{t}'] > EQUAL_PATTERN:
            fails[f'pattern part ({t})'] = out[f'd_pattern_{t}']
    out['d_raw_loss'] = D('all_loss', 'raw')
    if out['d_raw_loss'] > EQUAL_RAW:
        fails['raw-count loss'] = out['d_raw_loss']
    out['d_mean_effect'] = D('mean_effect', 'cal')
    if out['d_mean_effect'] < -EQUAL_MEAN:
        fails['mean APL effect'] = out['d_mean_effect']
    for k in ('amplitude', 'shape', 'sd_ratio', 'rho_effect', 'pearson_effect', 'loss', 'pattern'):
        out[f'd_{k}_all_cal'] = D(f'all_{k}', 'cal')
    out['d_pattern_all_raw'] = D('all_pattern', 'raw')
    out['apl_effect_ok'] = not any(k.startswith(('pattern', 'raw', 'mean')) for k in fails)
    # realism
    r = runs[(runs['input'] == name) & runs['unit'].isin(ALL_WIRINGS)]
    out['n_flagged'] = int(r['flagged'].sum())
    out['n_flagged_zcrit'] = int(sum(E108._flagged(x) for x in r.get('realism_flag_zcrit', pd.Series(dtype=object))))
    if out['n_flagged']:
        fails['realism (runs outside the bands)'] = out['n_flagged']
    out['realism_ok'] = out['n_flagged'] == 0
    # response rates (the mean over the eight runs)
    rr = runs[(runs['input'] == ref) & runs['unit'].isin(ALL_WIRINGS)]
    for c in ('on', 'off'):
        v = float(r[f'rate_rho_{c}'].astype(float).mean()) - float(rr[f'rate_rho_{c}'].astype(float).mean())
        out[f'd_rate_rho_{c}'] = v
        if v < -RATE_RHO_SLACK:
            fails[f'response rate ρ, APL {c}'] = v
    out['rates_ok'] = not any(k.startswith('response rate') for k in fails)
    # columns only
    for k in ('block', 'vcf_on', 'vcf_off', 'vcf_diff', 'struct_wt_on', 'struct_wt_off'):
        out[f'd_{k}'] = D(k)
    out['fails'] = fails
    out['constraints_ok'] = bool(out['apl_effect_ok'] and out['realism_ok'] and out['rates_ok'])
    out['rhos_not_lower'] = bool(min(out['d_rho_on_all'], out['d_rho_off_all']) >= -OTHER_SLACK)
    return out


def class_of(g: Mapping) -> str:
    """'better' (the gain, the control and every constraint hold), 'a price' (the gain holds and a constraint or the
    control fails), 'an equal' (no gain, every constraint holds and neither ρ lower by more than OTHER_SLACK), 'no gain'
    (otherwise); '' for an incomplete gate."""
    if not g.get('complete'):
        return ''
    if g.get('gain'):
        return 'a price' if g['fails'] else 'better'
    return 'an equal' if (g['constraints_ok'] and g['rhos_not_lower']) else 'no gain'


def best_gain(g: Mapping) -> float:
    """The largest population change against all recordings on a ρ that gains (NaN when neither gains)."""
    v = [g[f'd_rho_{c}_all'] for c in ('on', 'off') if g.get(f'gain_{c}')]
    return float(max(v)) if v else np.nan


def choose_c(gates: Mapping[str, Mapping]) -> dict:
    """Part C's exponent (the module docstring): of pow0.6 to pow0.9 classed 'better', the largest `best_gain` (on a
    tie the exponent nearer 1); none: FALLBACK_EXPONENT. `gates`: {input name: gate dict}."""
    cands = []
    for e in EXPONENTS:
        g = gates.get(pow_name(e))
        if g is not None and class_of(g) == 'better':
            cands.append((best_gain(g), e))
    if not cands:
        return dict(exponent=FALLBACK_EXPONENT, reason='no exponent of part A is "better": the entry\'s 0.75',
                    candidates=[])
    cands.sort(key=lambda t: (-t[0], abs(1.0 - t[1])))
    return dict(exponent=cands[0][1], reason='the "better" exponent with the largest gated gain', candidates=cands)


def wiring_gain(rd: pd.DataFrame, name: str, unit: str = 'real', ref: str = 'I0') -> bool:
    """One ρ higher by GAIN or more on one wiring against all recordings, the other not lower by more than OTHER_SLACK."""
    d = {c: _v(rd, name, unit, f'rho_{c}_all') - _v(rd, ref, unit, f'rho_{c}_all') for c in ('on', 'off')}
    return bool((d['on'] >= GAIN and d['off'] >= -OTHER_SLACK) or (d['off'] >= GAIN and d['on'] >= -OTHER_SLACK))


def outcome_of(gates: Mapping[str, Mapping], rd: pd.DataFrame, choice: Optional[Mapping] = None) -> dict:
    """The entry's outcome letter (the module docstring): {'letter', 'text', 'provisional'}. `gates` has part A's and
    part B's inputs and, once part C is done, REFIT."""
    exps = [pow_name(e) for e in EXPONENTS]
    pa = {n: gates[n] for n in exps if n in gates and gates[n].get('complete')}
    pb = {n: gates[n] for n in B_INPUTS if n in gates and gates[n].get('complete')}
    gc_ = gates.get(REFIT)
    done_c = bool(gc_ and gc_.get('complete'))
    out = dict(provisional=not done_c)
    pattern_failed = lambda g: any(k.startswith('pattern part') for k in g.get('fails', {}))
    if done_c:
        cl = class_of(gc_)
        if cl == 'better':
            return dict(out, letter='a', text='"better" after the refit')
        chosen = pow_name((choice or {}).get('exponent', FALLBACK_EXPONENT))
        ga = gates.get(chosen)
        if cl == 'a price' and pattern_failed(gc_) and ga is not None and class_of(ga) == 'a price' and pattern_failed(ga):
            return dict(out, letter='b', text='"a price" on the calibrated pattern part, at fixed scales and after the refit')
    if pa and not any(g.get('gain') for g in pa.values()):
        if any(wiring_gain(rd, n) for n in pa):
            return dict(out, letter='c', text='a gain on the hemibrain wiring only')
        best = lambda gs: {c: max([g[f'd_rho_{c}_all'] for g in gs.values()]) for c in ('on', 'off')}
        if pb:
            be, bs = best(pa), best(pb)
            differs = any(g.get('gain') for g in pb.values()) or any(abs(bs[c] - be[c]) >= GAIN for c in ('on', 'off'))
            if differs:
                return dict(out, letter='d', text='no gain, and the saturating form differs from the exponent',
                            best_exponent=be, best_saturating=bs)
    if not done_c and pa and any(g.get('gain') for g in pa.values()):
        return dict(out, letter='-', text='part A has a gain; the letter waits for the refit')
    return dict(out, letter='-', text='none of the entry\'s outcomes')


def gate_tables(rd: pd.DataFrame, runs: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, dict], dict, dict]:
    """gate.csv, gate.json, choice.json (written once parts A and B are complete, never changed after part C starts)
    and outcome.json."""
    gates = {n: gate(rd, runs, n) for n in A_INPUTS + B_INPUTS if n != CONTROL}
    gates[CONTROL] = gate(rd, runs, CONTROL, control=None)
    if _have(runs, REFIT):
        gates[REFIT] = gate(rd, runs, REFIT)
    cs = _read_json(out_dir() / 'choice.json')
    # the chosen exponent's held points (part Cstart) when it is not one of part A's: outcome (b)'s fixed-scale reading
    chosen = pow_name(float(cs['exponent'])) if cs.get('exponent') else None
    if chosen and chosen not in gates and _have(runs, chosen):
        gates[chosen] = gate(rd, runs, chosen)
    if all(gates[n].get('complete') for n in A_INPUTS + B_INPUTS):
        new = choose_c(gates)
        if cs.get('exponent') and float(cs['exponent']) != float(new['exponent']):
            raise RuntimeError(f"choice.json says {cs['exponent']}, the gate now {new['exponent']} (a choice never moves)")
        if not cs.get('exponent'):
            cs = dict(new, written=_now(), **H.code_commits())
            _write_json(out_dir() / 'choice.json', cs)
    rows = []
    for n, g in gates.items():
        flat = {k: v for k, v in g.items() if not isinstance(v, (dict, list))}
        flat['fails'] = '; '.join(f'{k} {v:+.4f}' if isinstance(v, float) else f'{k} {v}' for k, v in g.get('fails', {}).items())
        for c, v in (g.get('control_d_rho') or {}).items():
            flat[f'control_d_rho_{c}'] = v
        flat['class'] = class_of(g)
        rows.append(flat)
    gt = pd.DataFrame(rows)
    gt.to_csv(out_dir() / 'gate.csv', index=False)
    oc = outcome_of(gates, rd, cs)
    _write_json(out_dir() / 'gate.json', dict(gates=gates, classes={n: class_of(g) for n, g in gates.items()}, **H.code_commits()))
    _write_json(out_dir() / 'outcome.json', dict(oc, choice=cs.get('exponent'), written=_now()))
    return gt, gates, cs, oc


# ---------------------------------------------------------------------------
# part D (only on a pass after the refit)

def e113_input(path: str) -> pd.DataFrame:
    """E113's best input (a CSV as E10_joint_fit.py's `--orn-input` reads it), in the lab input's layout."""
    from experiments.E10_joint_fit import read_orn_input
    return read_orn_input(path, H.apl16_orn_deltas())


def d1_runs(name: str) -> Dict[str, Dict[str, tuple]]:
    """D1's runs: on each of T's wirings, the map `name` on the validation2 input at the wiring's refit parameters, the
    breadth of the mapped megamat input (the refit's), APL on and acutely off, nothing re-held (saved under D1/)."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E107_validation2_panel as E107
    from experiments import E110_no_kckc as E110
    v2 = map_of(name)(E107.panel_input('validation2'))
    orn = input_of(name)
    out = {}
    for w in T_WIRINGS:
        tag = f'D1_{name}_{w}'
        o = load_run('D1', tag)
        if o is None:
            f = load_fit(w)
            if f is None:
                raise RuntimeError(f'D1: no refit on {w}')
            p = E110.params_of_row(f[0])
            _, info = E85.build(w, orn, cache_tag(orn, name), cache_dir=local_dir())
            b = info['breadth']
            del info
            sim, info = E83.build_simulator(w, orn=v2, breadth=b, params=p)
            o = {'on': E107._abbrev(sim(p)), 'off': E107._abbrev(sim(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)))}
            del sim, info
            gc.collect()
            save_run('D1', tag, o)
            prov = _read_json(out_dir() / 'provenance_D1.json') or dict(part='D1', script=SCRIPT, points={})
            prov['points'][tag] = dict(H.model_provenance(orn_deltas=v2), params=p)
            _write_json(out_dir() / 'provenance_D1.json', prov)
        out[w] = o
    return out


def d1_table(name: str) -> pd.DataFrame:
    """D1's readings: the mapped validation2 runs' population (T) against Remy's somata (E109's `remy_readings`), with
    E107's saved runs of the family as fitted beside."""
    from experiments import E106_rho_gap_parts as E106
    from experiments import E109_input_search_rho as E109
    runs = d1_runs(name)
    per = {w: {c: H.canonical_pairs(H.pair_distances(o[c][1], 'cosine')) for c in ('on', 'off')} for w, o in runs.items()}
    pop = E106.mean_pairs(list(per.values()))
    ref = E109.e107_reference(T_WIRINGS)
    rows = [dict(reading='map at the refits', input=name, **E109.remy_readings(pop['on'])),
            dict(reading="E107's family as fitted", input='I0', **E109.remy_readings(ref['pop']['on']))]
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'd1.csv', index=False)
    return df


def part_d(force: bool = False, e113: Optional[str] = None) -> dict:
    """Part D (the module docstring); d.json. Runs nothing unless part C's class is 'better' (or `force`)."""
    g = _read_json(out_dir() / 'gate.json')
    cl = g.get('classes', {}).get(REFIT)
    res = dict(refit_class=cl, forced=bool(force), checked=_now())
    if cl != 'better' and not force:
        res['ran'] = False
        res['note'] = f"part C's class is {cl!r}: part D runs only on a pass after the refit"
        _write_json(out_dir() / 'd.json', res)
        print(res['note'], flush=True)
        return res
    name = chosen_name()
    res['d1'] = d1_table(name).to_dict('records')
    if e113:
        ref = e113_input(e113)
        _CACHE.setdefault('inputs', {})[E113_REF] = ref
        _CACHE['inputs'][E113_MAPPED] = map_of(name)(ref)
        for n in (E113_REF, E113_MAPPED):
            write_input(n, _CACHE['inputs'][n])
        part_sims('D2')
        res['d2'] = dict(e113_input=e113, map=name)
    else:
        res['d2'] = 'E113 has no best input (no --e113-input): D2 not run'
    res['ran'] = True
    _write_json(out_dir() / 'd.json', res)
    return res


# ---------------------------------------------------------------------------
# beside: E112's readings and the stages

def e112_ended() -> Tuple[bool, str]:
    """Whether E112 has ended (its index row says 'completed'), and the row's state."""
    p = H.shared_dir() / 'experiment_priorities.md'
    for line in p.read_text().splitlines():
        if line.startswith('| E112 |'):
            st = line.split('|')[3].strip()
            return st == 'completed', st
    return False, 'not found'


def e112_table(runs: pd.DataFrame) -> pd.DataFrame:
    """Beside: the population of every complete input under E112's readings R0 to R3 (`run_readings`,
    `reading_scores`), with E112's state; e112.csv. The carrying reading is E112's to name."""
    from experiments import E112_rho_like_for_like as E112
    ended, state = e112_ended()
    rows = []
    for name in runs['input'].unique():
        if not _have(runs, name):
            continue
        per = [E112.run_readings(run_outputs(name, w), E112.COUNT_READINGS, seed=E112.wiring_seed(w)) for w in ALL_WIRINGS]
        for r in E112.COUNT_READINGS:
            pop = E112.population([p[r] for p in per])
            rows.append(dict(input=name, reading=r, label=E112.READING_LABEL.get(r, r), e112_ended=ended, e112_state=state,
                             **E112.reading_scores(pop)))
    df = pd.DataFrame(rows)
    if len(df):
        base = df[df['input'] == 'I0'].set_index('reading')
        for c in ('rho_on_all120', 'rho_off_all120', 'pattern', 'loss', 'mean_effect'):
            df[f'd_{c}'] = [v - base[c].get(r, np.nan) for r, v in zip(df['reading'], df[c])]
    df.to_csv(out_dir() / 'e112.csv', index=False)
    return df


def stage_inputs(runs: pd.DataFrame) -> List[Tuple[str, pd.DataFrame, Dict[str, float]]]:
    """(name, input, the real wiring's parameters) of every input with a run on the real wiring: I0 at the fit, a held
    input at its held point, REFIT (the chosen map) at its refit."""
    out = []
    for _, r in runs[runs['unit'] == 'real'].iterrows():
        n = r['input']
        if n == 'I0':
            out.append((n, lab_input(), dict(E83.saved_fit('real')['p'])))
        elif n == REFIT:
            from experiments import E110_no_kckc as E110
            out.append((n, input_of(chosen_name()), E110.params_of_row(load_fit('real')[0])))
        elif isinstance(r.get('params_on'), str) and n not in (E113_REF, E113_MAPPED):
            out.append((n, input_of(n), json.loads(r['params_on'])))
    return out


def pn_stage_row(name: str, orn: pd.DataFrame) -> dict:
    """The PN stage of an input (E90's `al_stage` at the bouton base's parameters, `pn_statistics`)."""
    from experiments import E90_model_stages as E90
    st = E90.al_stage(None, orn=orn)
    return dict(input=name, **E90.pn_statistics(st['win'], st['pre']))


def bouton_stage_rows(name: str, orn: pd.DataFrame, params: Mapping[str, float], save_dir: Optional[Path] = None
                      ) -> List[dict]:
    """The bouton stage on the real wiring, APL on, at `params` (the module docstring): group means of 1 - Pearson across
    the boutons of the shared glomeruli, 'rate' and 'session_z', beside the flies'. `save_dir`: also save the stage's
    bouton windows as `<name>_on_bouton_win.parquet` and `<name>_on_bouton_pre.parquet` there (E90's `read_stages`
    frames; E120's part A4 reads them)."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E90_model_stages as E90
    sim, info = E85.build('real', orn, cache_tag(orn, name), cache_dir=local_dir())
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    mp.kc.save_inh_sims = True
    r, s = sim(dict(params))
    W = E83.wiring_parts('real')['W']
    st = E90.read_stages(mp, rv, W, order, [H.odor_abbrev(c) for c in r.columns], E90.MODEL_WINDOW)
    if save_dir is not None:
        for k in ('bouton_win', 'bouton_pre'):
            st[k].to_parquet(Path(save_dir) / f'{name}_on_{k}.parquet')
    del sim, info
    gc.collect()
    U = E90.stage_units(st, 'bouton')
    meta = E90.unit_meta('real', W)
    V = E90._drop_silent(U.loc[E90.subset_mask('bouton', U.index, meta, 'shared')])
    fl = {k: E90.flies_group(k) for k in ('bouton_measured', 'bouton_noise_removed')}
    rows = []
    for reading, X in (('rate', V), ('session_z', E90.session_z(V))):
        gm = E90.group_means(E90.odor_distances(X))
        for g, v in gm.items():
            rows.append(dict(input=name, reading=reading, group=g, model=v, n_boutons=int(V.shape[0]),
                             **{f'flies_{k.split("_", 1)[1]}': float(fl[k].loc[g, 'mean']) if g in fl[k].index else np.nan
                                for k in fl}))
    return rows


def part_beside() -> dict:
    """Beside (simulations, no class): pn_stage.csv and bouton_stage.csv, resumable per input."""
    runs = runs_table()
    pp, bp = out_dir() / 'pn_stage.csv', out_dir() / 'bouton_stage.csv'
    have_p = set(pd.read_csv(pp)['input']) if pp.exists() else set()
    have_b = set(pd.read_csv(bp)['input']) if bp.exists() else set()
    for name, orn, p in stage_inputs(runs):
        if name not in have_p:
            E83._append_row(pp, pn_stage_row(name, orn))
            print(f'E114 beside PN {name}', flush=True)
        if name not in have_b:
            for row in bouton_stage_rows(name, orn, p):
                E83._append_row(bp, row)
            print(f'E114 beside boutons {name}', flush=True)
    _write_json(out_dir() / 'provenance_beside.json', dict(part='beside', script=SCRIPT, **H.code_commits()))
    return dict(pn=str(pp), bouton=str(bp))


def part_stages() -> pd.DataFrame:
    """Part stages (added 2026-11-20 for E120's part A4, after the chain): beside's bouton stage run again for every
    input with its bouton windows saved under outputs/E114_input_compression/stages/ (resumable per input), and its rows
    checked against bouton_stage.csv (stage_check.csv: the largest absolute difference per input)."""
    runs = runs_table()
    d = local_dir('stages')
    bp = out_dir() / 'bouton_stage.csv'
    ref = pd.read_csv(bp) if bp.exists() else pd.DataFrame(columns=['input', 'reading', 'group', 'model'])
    rows = []
    for name, orn, p in stage_inputs(runs):
        if (d / f'{name}_on_bouton_win.parquet').exists() and (d / f'{name}_on_bouton_pre.parquet').exists():
            continue
        new = pd.DataFrame(bouton_stage_rows(name, orn, p, save_dir=d))
        old = ref[ref['input'] == name]
        m = new.merge(old, on=['input', 'reading', 'group'], suffixes=('', '_saved'))
        diff = float((m['model'] - m['model_saved']).abs().max()) if len(m) else np.nan
        rows.append(dict(input=name, n_rows=len(new), n_matched=len(m), max_abs_diff=diff))
        print(f'E114 stages {name}: max |difference| from bouton_stage.csv {diff:.3g}', flush=True)
    sc = out_dir() / 'stage_check.csv'
    out = pd.concat([pd.read_csv(sc), pd.DataFrame(rows)], ignore_index=True) if sc.exists() and rows else pd.DataFrame(rows)
    if len(out):
        out.drop_duplicates('input', keep='last').to_csv(sc, index=False)
    _write_json(out_dir() / 'provenance_stages.json', dict(part='stages', script=SCRIPT, saved=str(d), **H.code_commits()))
    return out


# ---------------------------------------------------------------------------
# the tables

COMPARISON_COLS = (
    ('ρ on odor-odor distances, APL on (raw)', 'rho_on_all', 'raw'), ('ρ on odor-odor distances, APL off (raw)', 'rho_off_all', 'raw'),
    ('APL-effect loss (readout)', 'all_loss', 'cal'), ('pattern part (readout)', 'all_pattern', 'cal'),
    ('amplitude part (readout)', 'all_amplitude', 'cal'), ('shape part (readout)', 'all_shape', 'cal'),
    ('spread ratio (readout)', 'all_sd_ratio', 'cal'), ('mean APL effect (readout)', 'mean_effect', 'cal'),
    ('ρ on the per-pair APL effect (readout)', 'all_rho_effect', 'cal'), ('Pearson r on the per-pair APL effect (readout)', 'all_pearson_effect', 'cal'),
    ('APL-effect loss (raw)', 'all_loss', 'raw'), ('pattern part (raw)', 'all_pattern', 'raw'), ('mean APL effect (raw)', 'mean_effect', 'raw'),
    ('response rate, APL on', 'run_rate_on', 'raw'), ('response rate, APL off', 'run_rate_off', 'raw'),
    ('APL-off / APL-on response rate', 'run_ratio', 'raw'),
    ('response rate ρ, APL on', 'run_rate_rho_on', 'raw'), ('response rate ρ, APL off', 'run_rate_rho_off', 'raw'),
    ('odors per KC (TV), APL on', 'run_opk_on_tv', 'raw'), ('odors per KC (TV), APL off', 'run_opk_off_tv', 'raw'),
    ('largest spike count, APL on', 'run_max_count_on', 'raw'),
    ('the block (raw)', 'block', 'raw'), ('natural-source alignment, APL on (raw)', 'vcf_on', 'raw'),
    ('natural-source alignment, APL off (raw)', 'vcf_off', 'raw'), ('natural-source alignment, APL on minus off (raw)', 'vcf_diff', 'raw'),
    ('structure vs wild-type somata (raw)', 'struct_wt_on', 'raw'), ('runs outside the bands', 'n_flagged', 'raw'))


def comparison_table(rd: pd.DataFrame, gt: pd.DataFrame, choice: Mapping, oc: Mapping) -> str:
    """comparison_table.md: every input's population of eight beside I0, the gate per input, the choice and the outcome."""
    names = [n for n in ('I0',) + A_INPUTS + B_INPUTS + (pow_name(FALLBACK_EXPONENT), REFIT, E113_REF, E113_MAPPED)
             if len(rd[(rd['input'] == n) & (rd['unit'] == 'all8')])]
    L = ['# E114: a compressive map from the ORN recordings to the input rates', '',
         "paa's input is the 10^-4 recording for the flies' 10^-3; euc and -aPine stay at 10^-1.75 and 10^-1.5 for the "
         "flies' 10^-2; the map acts on every odor.", '', '## The family, mean matrix over its eight wirings', '',
         '| statistic | ' + ' | '.join(names) + ' |', '|---|' + '---:|' * len(names)]
    for lab, col, ro in COMPARISON_COLS:
        spec = '.0f' if col == 'n_flagged' else '.3f'
        L.append(f'| {lab} | ' + ' | '.join(_f(_v(rd, n, 'all8', col, ro), spec) for n in names) + ' |')
    if len(gt):
        ok = lambda x: '-' if x is None or (isinstance(x, float) and np.isnan(x)) else ('ok' if x else 'fails')
        L += ['', '## The gate (the eight wirings, against the family at the lab\'s input)', '',
              '| input | ρ change on / off | wirings ≥ 0.010 on / off | halves A-D on (min) / off (min) | control | APL effect | '
              'realism | response rates | class | fails |', '|---|---:|---:|---:|---|---|---|---|---|---|']
        for r in gt.to_dict('records'):
            if not r.get('complete'):
                L.append(f'| {r["input"]} | runs missing | | | | | | | | |')
                continue
            hmin = {c: min(r[f'd_rho_{c}_{h}'] for h in HALVES) for c in ('on', 'off')}
            L.append(f'| {r["input"]} | {_f(r["d_rho_on_all"], "+.3f")} / {_f(r["d_rho_off_all"], "+.3f")} | '
                     f'{int(r["n_wirings_gain_on"])} / {int(r["n_wirings_gain_off"])} | {_f(hmin["on"], "+.3f")} / '
                     f'{_f(hmin["off"], "+.3f")} | {ok(r.get("control_ok"))} | {ok(r["apl_effect_ok"])} | {ok(r["realism_ok"])} | '
                     f'{ok(r["rates_ok"])} | {r.get("class", "")} | {r.get("fails", "")} |')
    L += ['', f"Part C's exponent: {choice.get('exponent', '-')} ({choice.get('reason', 'not chosen yet')}).",
          f"Outcome: ({oc.get('letter')}) {oc.get('text')}{' (provisional)' if oc.get('provisional') else ''}.", '']
    s = '\n'.join(L) + '\n'
    (out_dir() / 'comparison_table.md').write_text(s)
    return s


def realism_table(rd: pd.DataFrame, choice: Mapping) -> str:
    """realism_table.md: the seed (I0), the best variant on the APL-effect loss (readout) among the inputs, part C's
    input and the refit, beside the KC data's values and bands."""
    from kcq import apl_effect as A
    pop = rd[(rd['unit'] == 'all8') & (rd['readout'] == 'cal')]
    best = pop.loc[pop['all_loss'].idxmin(), 'input'] if len(pop) else None
    names = [n for n in dict.fromkeys(['I0', best, pow_name(choice.get('exponent', FALLBACK_EXPONENT)), REFIT]) if n
             and len(pop[pop['input'] == n])]
    data = dict(rate_on=H.data_target('ctrl_frac'), rate_off=H.data_target('tnt_frac'))
    rows = [('response rate, APL on (band 0.045 to 0.085)', 'run_rate_on', f"{data['rate_on']:.3f}"),
            ('response rate, APL off', 'run_rate_off', f"{data['rate_off']:.3f}"),
            ('APL-off / APL-on (band 1.6 to 2.9)', 'run_ratio', f"{data['rate_off'] / data['rate_on']:.2f}"),
            ('response rate ρ, APL on', 'run_rate_rho_on', '-'), ('response rate ρ, APL off', 'run_rate_rho_off', '-'),
            ('odors per KC (TV), APL on', 'run_opk_on_tv', '-'), ('odors per KC (TV), APL off', 'run_opk_off_tv', '-'),
            ('response rate ρ, APL on (Yang\'s z)', 'run_rate_rho_on_zcrit', '-'),
            ('runs outside the bands (of 8)', 'n_flagged', '0')]
    L = ['| statistic | ' + ' | '.join(names) + ' | KC data |', '|---|' + '---:|' * (len(names) + 1)]
    for lab, col, dv in rows:
        L.append(f'| {lab} | ' + ' | '.join(_f(_v(rd, n, 'all8', col, 'raw'), '.0f' if col == 'n_flagged' else '.3f')
                                            for n in names) + f' | {dv} |')
    s = '\n'.join(L) + '\n'
    (out_dir() / 'realism_table.md').write_text(s)
    return s


def score() -> dict:
    """Part score (no simulations; E112's readings beside are computed from the saved counts)."""
    rd = readings()
    runs = runs_table()
    gt, gates, choice, oc = gate_tables(rd, runs)
    try:
        e112_table(runs)
    except Exception as e:      # beside, no class: a failure here is reported, not fatal
        print(f'E114 score: e112 beside failed: {e!r}', flush=True)
    comparison_table(rd, gt, choice, oc)
    realism_table(rd, choice)
    _write_json(out_dir() / 'provenance_score.json', dict(part='score', script=SCRIPT, **H.code_commits()))
    print(json.dumps(dict(classes={n: class_of(g) for n, g in gates.items()}, choice=choice.get('exponent'), outcome=oc),
                     indent=1, default=str), flush=True)
    return dict(classes={n: class_of(g) for n, g in gates.items()}, choice=choice, outcome=oc)


def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts), one row per held run and refit."""
    runs = runs_table()
    rows, parts = [], []
    for _, r in runs[runs['input'] != 'I0'].iterrows():
        o = run_outputs(r['input'], r['unit'])
        # a refit's tag is its fit's own stem, the name of its archived spike counts and responses (E110's layout), so
        # kcq/apl_effect.py finds its responses for the bands; `variant` keeps C_<wiring>
        tag = r['fit_tag'] if r['input'] == REFIT else r['tag']
        kind = (f"C2's family's {r['unit']} " + ('fitted again (E10_joint_fit.py, E72\'s part-B command) on ' if r['input'] == REFIT
                else 'fit, the somatic threshold and the tied APL scale re-held, on ') +
                f"the ORN input {r['input'] if r['input'] != REFIT else chosen_name()} (results/{EXPERIMENT}/inputs/)")
        rows.append(dict(tag=tag, variant=r['tag'], part=r['part'], arm=r['input'], base='bouton', k=2.0, breadth_k=2.0, kind=kind,
                         config=r['input'], wiring=r['unit'], ref_experiment=r.get('fit_folder'), ref_fit=r.get('fit_tag'),
                         **{k: r.get(k) for k in E83.PARAM_KEYS}, kc_kc_inh=r.get('wKCKC_scale'),
                         thr_factor=r.get('thr_factor'), apl_factor=r.get('apl_factor'), rate_on=r['rate_on'],
                         rate_off=r['rate_off'], max_count_on=r['max_count_on'], max_count_off=r['max_count_off'],
                         realism_flag=r['realism_flag'], subdir=(f'fits/C_{r["unit"]}' if r['input'] == REFIT else r['part'])))
        parts.append(pd.concat({(tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out_dir() / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('inputs', '0', '0c', 'A', 'B', 'score', 'Cstart', 'C', 'D', 'beside',
                                                      'variants', 'stages'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=ALL_WIRINGS)
    ap.add_argument('--max-points', type=int, default=None)
    ap.add_argument('--force', action='store_true', help='part D: run without a pass after the refit')
    ap.add_argument('--e113-input', default=None, help="part D2: E113's best input (CSV, relative to the shared directory)")
    ap.add_argument('--rerun', action='store_true', help='part 0c: run the check again')
    a = ap.parse_args()
    if a.part == 'inputs':
        write_inputs()
    elif a.part == '0':
        identity(a.wirings)
    elif a.part == '0c':
        fit_path_check(a.rerun)
    elif a.part in ('A', 'B'):
        part_sims(a.part, a.wirings, a.max_points)
    elif a.part == 'score':
        score()
    elif a.part == 'Cstart':
        part_sims('Cs', a.wirings, a.max_points)
        print(json.dumps(write_starts(), indent=1), flush=True)
    elif a.part == 'C':
        part_c(a.wirings)
    elif a.part == 'D':
        part_d(a.force, a.e113_input)
    elif a.part == 'beside':
        part_beside()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'stages':
        part_stages()


if __name__ == '__main__':
    main()
