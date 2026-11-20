#!/usr/bin/env python3
"""E113: the input of paa and B-myr without a search: each ORN fly's own recording, and the three entries that two
independent searches both named (DM1 for both odors, DP1m for paa), read on the family's eight wirings.

The entry is item E113 of experiment_priorities.md (proposed by viki 2026-11-19 13:58; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 13:50; moved into the list 13:58 on Tom's request of 13:49 relayed by hal), as revised
2026-11-20 10:11: cut to its no-search parts after E116's outcome (d) (planned by a Fable subagent, claude-fable-5-2,
spawned from viki 2026-11-20 10:08). Held, not run (reopened only on Tom's word): arms W and K, the penalty arm, their
controls and the refit with E109's input. Canonical name: "the searched input of paa, B-myr and 2-but inside the ORN
flies' spread, reduced and fitted again" (kept for the record; the page's title may say "without a search").

"The family" is C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at
its own fit (E83's `saved_fit`; S: real, P1, P2, M1; T: P3, P4, M2, M3). "The two ρ's" are ρ on odor-odor distances with
APL on against the control flies and with APL off against APL output blocked (TNT), on raw counts (cosine distance on
spike counts), each KC matrix the mean over all recordings ('all') or over one of E108's halves A to D (E108's
`dist_targets`; A and B are E77's `recording_halves(0)`). Every reading of an input is E85's `held_reading` (the breadth
recomputed for the input, the somatic threshold bisected to the fit's APL-off response count, then the tied APL scale to
its APL-on count), as E116's part H ran it, scored by E108's `unit_scores` and gated with E109's thresholds, read as
E116's part R (E116's `change_values`, `_measures`; E109's `class_of`).

Concentrations. paa's input is the lab's recording at 10^-4 for the flies' 10^-3, so every number involves it; euc and
-aPine stay at the recorded 10^-1.75 and 10^-1.5 for the flies' 10^-2 and are in every 120-pair statistic; B-myr is at
the flies' 10^-2 (`H.APPROX_CONC_NOTE` in every provenance file).

The inputs (glomeruli x the 16 odors, the lab's 36 rows; inputs/<name>.csv with sha256.json and inputs.json; fixed here
before any run; only paa's and B-myr's columns ever differ from the lab's input, `H.apl16_orn_deltas()`):
- Arm D, the entries, no search; every other entry the lab's. Sizes are log2 factors on the lab's response:
  - D1: DM1 lowered for paa and for B-myr by E109's mean factors read at 0.6 (results/E109_input_search_rho/
    step6_entries.csv, mean_log2_read): paa -0.5970922190910033, B-myr -0.4023030617464084 (the entry's -0.60, -0.40).
    Reported with the caveat that the sizes came from searches on Yang's recordings.
  - D1x: DM1 lowered for both by the factor-of-two bound (log2 -1.0). Same caveat.
  - D2: E116's mean factors read at 0.6 (results/E116_input_second_flies/entries.csv, mean_log2_read): DM1 for paa
    -0.2657535311916492, DM1 for B-myr -0.4021171825728125, DP1m for paa -0.4871020051360745 (the entry's -0.27, -0.40,
    -0.49). These owe nothing to Yang's recordings. D2 carries the class.
  `entry_sizes` checks the constants against both tables (to 1e-12) before inputs are written.
- The controls (D1's two factors elsewhere):
  - C1 to C4: D1's paa factor on the k-th of the four glomeruli nearest DM1 in the lab's paa response, and D1's B-myr
    factor on the k-th nearest DM1 in the lab's B-myr response. "Nearest": the smallest |lab(g) - lab(DM1)| for that
    odor, ties by the lab input's glomerulus order; excluded: DM1, DP1m and E109's eight recurring entries of both odors
    (DL5, VM7d, DP1l, VC3, VM5d). The list, fixed before any run (`control_glomeruli` on the lab's input; asserted):
      paa:   VL2a (lab 102.6 Hz; DM1 154.6), DM6 (70.2), DL2v (61.4), VA3 (55.1)
      B-myr: VM5v (57.9; DM1 72.7), DL2v (39.3), DC1 (37.9), DM4 (35.9)
    so C1 = (paa VL2a, B-myr VM5v), C2 = (DM6, DL2v), C3 = (DL2v, DC1), C4 = (VA3, DM4). DM6 is the glomerulus of E116's
    Hallem arm H_dm6 for paa; it is in the list by the rule.
  - D1inv: D1 inverted (DM1 raised by the same factors).
- Part 0, each fly's own input: the six validation2 ORN flies (olfdb dataset 4; E108's compact copy and units: each fly
  scaled to-avg-max, then E108's slope through the origin of the lab's four validation odors on the flies' mean; E109's
  `validation2_fly_values`), as E85's 'measured' model flies and E108's step 4: the lab's column plus the fly's deviation
  from the mean over the flies that have the glomerulus (0 where the fly lacks it). F<f>_paa: paa alone; F<f>_pb: paa and
  B-myr together (f = 0 to 5, the flies in E108's `bmyr_flies` order). B-myr alone is E108's step 4 (its B0 to B5, the
  same construction; `fly_deviations` asserts the deviations equal E108's): its runs are reused as F<f>_bmyr when step 1's
  check reproduces E108's B0 on the real wiring bit for bit, and run again in part sims otherwise. 2-but is dropped.

Parts (CLI `--part`; simulation parts resumable, a point done when its row is in any runs_sims_*.csv):
- inputs (no simulations): the inputs' files, sha256.json, inputs.json (each input's changed entries, as log2 factors
  and in the validation2 ORN flies' SDs, E109's `flies_spread`).
- identity (simulations, about 30 min on 16 threads, about 1 h on 8): step 1. On each wiring (`--wirings` restricts), the
  held reading of the lab's input reproduces the saved fit (0 count and 0 response differences, both factors exactly 1),
  as E109's and E116's part 0 found (their part0.json beside); on the real wiring, E108's B0 held again reproduces
  E108's saved run B_B0_real (bmyr_reuse). identity.json, merged across invocations; 'all_exact' needs all eight. If
  anything differs nothing else runs.
- sims (simulations; about 160 held readings, about 2 min each on 16 threads, 3 to 4 on 8): the plan (`plan`): arm D,
  then the controls, then part 0's 12 fly inputs (and the six B-myr-alone inputs if E108's are not reused), each on the
  eight wirings. `--lane i/n` takes every n-th point of the plan from the i-th (0-based), `--wirings` restricts; rows in
  runs_sims_<i>of<n>.csv (absolute parameters on and off, the factors, realism, the five odors' response fractions,
  `H.code_commits()`); runs under outputs/E113_input_within_spread/R/; each point's `H.model_provenance` in
  provenance_sims_<i>of<n>.json.
- remy (simulations, beside, no class; 8 runs): step 6, D1 and D2 on the validation2 ORN input (E107's `panel_input`,
  paa's and B-myr's columns times the arm's factors: E116's `v2_with` at size 1), T's four wirings as fitted (E109's
  `fitted_on`), read against Remy's validation2 somata (E116's `remy_targets`: her mean table and her halves a to d, over
  the 136 pairs and the 31 of paa and B-myr) with E107's as-fitted runs of the family beside; runs_V.csv, outputs/.../V/.
- score (no simulations): readings.csv (E108's `unit_scores` per change, unit and readout: each wiring, all8, S, T),
  runs.csv, gate.csv and gate.json (each arm, part 0's rule per kind, the refit's reading beside), the outcome letters,
  refit.json (the refit's input, or none), controls.csv, flies.csv, remy_readings.csv, e112 beside, realism.csv,
  comparison_table.md, provenance_score.json.
- refit (fits; only when score's refit.json names an input, i.e. after a pass): step 5, the family fitted again with that
  input on T's four wirings: E72's part-B command as E114's part C (`-k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl
  --apl-readout both --tag _cap3_refit --maxfev 150`, the wiring's options as E110's `wiring_args`) with `--orn-input
  results/E113_input_within_spread/inputs/<name>.csv --init-params results/E113_input_within_spread/starts.json:R_<w>`,
  into results/ and outputs/ of E113_input_within_spread/fits/refit_<w>/; runs_refit.csv. Then score again. `--dry-run`
  prints the commands.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive and kcq/apl_effect.py.

The gate (the entry's: E109's thresholds, read as E116's part R; `gate_arm`):
- Each arm on T (population of P3, P4, M2, M3) against the family at the lab's input (I0): ρ with APL on higher by 0.010
  or more, or with APL off by 0.015 or more, against all recordings and against each of halves A and B, on T's population
  and on at least 3 of T's 4 wirings against each; the other ρ not lower by more than 0.005 against any of the three.
- Controls: each of the five (C1 to C4, D1inv), on each ρ that gains, below half the arm's gain against all recordings on
  T (every control, not their mean: the entry says "each control"). The same five serve D1, D1x and D2.
- The other measures (E116's part R CHECKS, against all recordings, on T): the pattern part under the calibrated readout
  not higher by more than 0.003, the raw-count loss by 0.005, the mean APL effect (calibrated) not lower by 0.01, every
  run of the arm on T inside the four bands, response rate ρ with APL on and off (the mean over T's runs) not lower by
  0.05. The block, the natural-source alignment and the structure vs wild-type somata are columns only (the entry: "The
  block is a column only"): their failures are listed but no class reads them (the class is `class_without_pending`).
- Classes (E109's `class_of`): "passes", "a price", "carries below the gate", "does not carry". All eight wirings are
  reported beside (6 of 8), with no class.
- Part 0's rule (`gate_flies`, per kind: paa alone, paa and B-myr): no fly is picked on the readings it is judged by. The
  fly whose input gives S's population the highest ρ with APL on against half A is read on T against half B, and the
  fly with the highest against half B is read on T against half A. Each direction: a gain on T's population against the
  held-out half and on at least 3 of T's 4 wirings against it, the other ρ against that half not lower by more than
  0.005; a pass needs both directions on the same ρ. No controls (the fly-to-fly spread is the comparison); the other
  measures as above, the worst over the two directions. The B-myr-alone flies are read the same way, beside, and enter
  no outcome.
- The table carries the amplitude part, the shape part, the spread ratio and Pearson r on the per-pair APL effect (each
  calibrated, as the change from I0), and paa's distance from B-myr.

Outcomes (`outcomes`): (a) D2 passes; (b) D1 or D1x passes and D2 does not; (c) a kind of part 0 passes under the rule;
(d) none passes and something carries below the gate; (e) none passes and nothing carries, or every arm's controls gain
as much as the arm on each ρ. Several letters can hold ((a) with (c)).

Choices the entry left open (fixed here before any run):
- The sizes are the tables' exact means read at 0.6, not the entry's rounded values.
- The refit's input (`refit_choice`): the first that passes of D2, D1, D1x; else the passing kind of part 0 with the
  larger mean gain (summed over the ρ's that pass), and of its two picks the same fly if both directions picked it,
  otherwise the one with the higher ρ with APL on of S's population against all recordings (S only: T is what it is
  judged on). The refit starts from each T wiring's held point of that input (the fit's own parameters with the somatic
  threshold and the tied APL scale re-held: E114 part C's and E110 part B's precedent), 150 evaluations, KC-KC free.
  It is reported beside with E116's values on T against I0 and a parameter table against the family's fits (any of the
  threshold and the APL scales that moves by more than 20 % named); it is not gated.
- "Reported in the flies' SDs": each changed entry's log2 factor over the validation2 ORN flies' SD of log2(fly / mean)
  for that odor and glomerulus (E109's `flies_spread`).
- E112's reading (beside): E112 had not ended when this was written (its index row 'running'). score gives, when E112 has
  ended, the two ρ's of each arm's and control's population (all8, T) under E112's readings R0 to R3 (E112's
  `run_readings`, `reading_scores`) and names none as carrying; otherwise a line saying it had not ended.

Departures from the entry: none so far.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no connectome table is read beyond the model's own hemibrain
wiring and E72's re-draws of it (the eight wirings); nothing of FlyWire, FAFB or BANC; no table of Yang's own Claude's
project; no raw-only recording (the KC data are Yang's processed APL16 matrices, the ORN data olfdb dataset 4 through
E108's compact copy). Nothing is published to claude.ai.

Provenance: E116's `commits_clean` after importing this experiment's modules (E108's check: the other checkouts clean,
this repository clean in every file a process loads); every row records `H.code_commits()`; every run's provenance has
`H.model_provenance(orn_deltas=...)`; the refit's processes E110's `_fit_clean` (nothing modified under experiments/ or
kcq/).

Outputs: <shared>/results/E113_input_within_spread/; spike counts and responses under outputs/E113_input_within_spread/
(R, V, fits/), then data/model_outputs/E113_input_within_spread/. The scoring reads a run from this host's outputs, else
from the archive, so score runs on any host once the runs are archived.

Usage (from the repository root). On atlas (8 threads a lane; OMP_THREAD_LIMIT as well, CLAUDE.md "Gotchas"; under
systemd-run --user or nohup, with -u):
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part inputs
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E113_input_within_spread.py --part identity
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E113_input_within_spread.py --part sims --lane 0/1
  (two lanes at once when the machine has 16 cores free: `--lane 0/2` and `--lane 1/2` in two processes)
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E113_input_within_spread.py --part remy
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part score
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part variants
    ./venv/bin/python archive_tools/archive_model_outputs.py E113_input_within_spread --exclude 'input-E113_*,wiring-*'
    ./venv/bin/python archive_tools/migrate_model_output_params.py E113_input_within_spread --apply
  Only when score's refit.json names an input (a pass):
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part refit --dry-run
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E113_input_within_spread.py --part refit [--wirings P3 P4]
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part score
    ./venv/bin/python -u experiments/E113_input_within_spread.py --part variants
  then archive and migrate again. ('input-E113_*' and 'wiring-*' at the top of outputs/E113_input_within_spread/ are
  the breadth caches of E85's `build` (the lab wiring's and the re-draws'), which any host recomputes.)
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E109_input_search_rho as E109
from experiments import E116_input_second_flies as E116

EXPERIMENT = 'E113_input_within_spread'
SCRIPT = 'experiments/E113_input_within_spread.py'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
PAA_BMYR = ('paa', 'B-myr')
YANG_TARGETS = ('all', 'A', 'B')
UNSEEN = {'A': 'B', 'B': 'A'}
# arm D (log2 factors on the lab's response; the module docstring)
D1_LOG2 = {('paa', 'DM1'): -0.5970922190910033, ('B-myr', 'DM1'): -0.4023030617464084}
D2_LOG2 = {('paa', 'DM1'): -0.2657535311916492, ('B-myr', 'DM1'): -0.4021171825728125, ('paa', 'DP1m'): -0.4871020051360745}
BOUND_LOG2 = -1.0
ARMS = ('D2', 'D1', 'D1x')
CLASS_ARM = 'D2'
SIZE_SOURCES = {'D1': ('E109_input_search_rho', 'step6_entries.csv'), 'D2': ('E116_input_second_flies', 'entries.csv')}
EXCLUDED = ('DM1', 'DP1m') + tuple(sorted({g for gs in E116.E109_RECURRING.values() for g in gs} - {'DM1', 'DP1m'}))
CONTROL_GLOMERULI = {'paa': ('VL2a', 'DM6', 'DL2v', 'VA3'), 'B-myr': ('VM5v', 'DL2v', 'DC1', 'DM4')}
N_NEAREST = 4
CONTROLS = ('C1', 'C2', 'C3', 'C4', 'D1inv')
# part 0
N_FLIES = 6
FLY_KINDS = ('paa', 'pb')                  # paa alone; paa and B-myr together
BESIDE_KIND = 'bmyr'                       # B-myr alone (E108's step 4)
FLY_ODORS = {'paa': ('paa',), 'pb': ('paa', 'B-myr'), 'bmyr': ('B-myr',)}
# the gate (E109's thresholds, E116's part R reading)
GAIN, OTHER_SLACK, CONTROL_RATIO = E109.GAIN, E109.OTHER_SLACK, E109.CONTROL_RATIO
T_GAIN_WIRINGS, ALL8_GAIN_WIRINGS = 3, 6
REFIT = 'refit'
FIT_PREFIX = f'{EXPERIMENT}/fits'
MAXFEV = 150
MOVE_NAMED = 0.20
PRELOAD = E116.PRELOAD + ('experiments.E116_input_second_flies', 'experiments.E110_no_kckc',
                          'experiments.E112_rho_like_for_like', 'experiments.E77_input_screen')

_CACHE: Dict[object, object] = {}


# ---------------------------------------------------------------------------
# folders, files, provenance

def out_dir() -> Path:
    """results/E113_input_within_spread/ (E113_OUT_DIR in the environment instead, for smoke tests)."""
    d = Path(os.environ['E113_OUT_DIR']) if os.environ.get('E113_OUT_DIR') else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    """outputs/E113_input_within_spread/<sub> (E113_LOCAL_DIR in the environment instead)."""
    root = Path(os.environ['E113_LOCAL_DIR']) if os.environ.get('E113_LOCAL_DIR') else H.outputs_dir(EXPERIMENT)
    d = root / sub if sub else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _host() -> str:
    import socket
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    E109._write_json(path, obj)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def commits_clean() -> dict:
    """E116's `commits_clean` (E108's check) after importing this experiment's modules (PRELOAD)."""
    if 'commits' not in _CACHE:
        for m in PRELOAD:
            importlib.import_module(m)
        _CACHE['commits'] = E116.commits_clean()
    return dict(_CACHE['commits'])


def _provenance_add(path: Path, tag: str, prov: dict) -> None:
    d = _read_json(path) or dict(script=SCRIPT, points={})
    d['points'][tag] = prov
    _write_json(path, d)


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment: this host's outputs, else the model-output archive; or None."""
    for d in (local_dir(part), archive_dir() / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


# ---------------------------------------------------------------------------
# the inputs

def lab_input() -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()`), E108's reader."""
    from experiments import E108_paa_extrapolation as E108
    return E108.lab_input()


def entry_sizes(check: bool = True) -> Dict[str, Dict[Tuple[str, str], float]]:
    """{arm or control: {(odor, glomerulus): log2 factor}} (the module docstring). With `check`, D1's and D2's constants are
    compared with E109's and E116's tables (mean_log2_read, to 1e-12) and the control list with `control_glomeruli`."""
    if check:
        for arm, const in (('D1', D1_LOG2), ('D2', D2_LOG2)):
            exp, f = SIZE_SOURCES[arm]
            t = pd.read_csv(H.results_dir(exp) / f).set_index(['odor', 'glomerulus'])
            for k, v in const.items():
                got = float(t.loc[k, 'mean_log2_read'])
                if abs(got - v) > 1e-12:
                    raise RuntimeError(f'{arm} {k}: the constant {v} is not {exp}/{f} ({got})')
        cg = control_glomeruli()
        if {o: tuple(v) for o, v in cg.items()} != CONTROL_GLOMERULI:
            raise RuntimeError(f'the control glomeruli by the rule are {cg}, not {CONTROL_GLOMERULI}')
    out = {'D1': dict(D1_LOG2), 'D1x': {k: BOUND_LOG2 for k in D1_LOG2}, 'D2': dict(D2_LOG2)}
    for i in range(N_NEAREST):
        out[f'C{i + 1}'] = {(o, CONTROL_GLOMERULI[o][i]): D1_LOG2[(o, 'DM1')] for o in PAA_BMYR}
    out['D1inv'] = {k: -v for k, v in D1_LOG2.items()}
    return out


def control_glomeruli(lab: Optional[pd.DataFrame] = None, n: int = N_NEAREST,
                      excluded: Sequence[str] = EXCLUDED) -> Dict[str, Tuple[str, ...]]:
    """For each of paa and B-myr, the `n` glomeruli nearest DM1 in the lab's response to it (smallest |lab(g) -
    lab(DM1)|, ties by the input's glomerulus order), `excluded` left out."""
    lab = H.abbrev_columns(lab_input() if lab is None else lab).astype(float)
    out = {}
    for o in PAA_BMYR:
        c = lab[o]
        d = (c - c['DM1']).abs().drop([g for g in excluded if g in c.index])
        out[o] = tuple(d.sort_values(kind='stable').index[:n])
    return out


def factors_of(name: str, lab: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """An arm's or control's factors, glomeruli x abbreviated odors in the lab input's order (1 where unchanged)."""
    lab = lab_input() if lab is None else lab
    F = pd.DataFrame(1.0, index=lab.index, columns=[H.odor_abbrev(c) for c in lab.columns])
    for (o, g), l2 in entry_sizes(check=False)[name].items():
        F.loc[g, o] = float(2.0 ** l2)
    return F


def fly_deviations() -> dict:
    """The six validation2 ORN flies' paa and B-myr in the model's units (E109's `validation2_fly_values`): 'keys' (E108's
    `bmyr_flies` order), 'glomeruli', and per odor 'dev' (flies x glomeruli: the fly's value minus the mean over the
    flies that have the glomerulus, 0 where the fly lacks it) and 'values'. B-myr's deviations must equal E108's."""
    if 'flies' in _CACHE:
        return _CACHE['flies']
    from experiments import E108_paa_extrapolation as E108
    B = E108.bmyr_flies()
    V = E109.validation2_fly_values(PAA_BMYR)
    out = dict(keys=list(B['keys']), glomeruli=list(B['glomeruli']), dev={}, values={})
    for o, v in V.items():
        if v.shape[0] != N_FLIES:
            raise RuntimeError(f'{o}: {v.shape[0]} flies, not {N_FLIES}')
        with np.errstate(invalid='ignore'):
            mu = np.nanmean(v, axis=0)
        out['dev'][o] = np.where(np.isfinite(v), v - mu[None], 0.0)
        out['values'][o] = v
    if not np.array_equal(out['dev']['B-myr'], B['dev']):
        raise RuntimeError("B-myr's per-fly deviations differ from E108's step 4")
    _CACHE['flies'] = out
    return out


def fly_input(f: int, kind: str) -> pd.DataFrame:
    """Fly `f`'s input: the lab's with the fly's deviation added to the columns of FLY_ODORS[kind] ('paa' paa alone, 'pb'
    paa and B-myr, 'bmyr' B-myr alone, which is E108's B<f>)."""
    if kind == BESIDE_KIND:
        from experiments import E108_paa_extrapolation as E108
        return E108.input_of(f'B{f}')
    D = fly_deviations()
    orn = lab_input()
    if list(orn.index) != D['glomeruli']:
        raise RuntimeError("the lab input's glomeruli are not the flies'")
    for o in FLY_ODORS[kind]:
        col = [c for c in orn.columns if H.odor_abbrev(c) == o]
        if len(col) != 1:
            raise RuntimeError(f'{o}: {col}')
        orn[col[0]] = orn[col[0]].to_numpy(float) + D['dev'][o][f]
    return orn


def fly_name(f: int, kind: str) -> str:
    return f'F{f}_{kind}'


def parse_fly(name: str) -> Optional[Tuple[int, str]]:
    if name.startswith('F') and '_' in name:
        a, b = name[1:].split('_', 1)
        if a.isdigit() and b in FLY_KINDS + (BESIDE_KIND,):
            return int(a), b
    return None


def input_names(bmyr: bool = False) -> List[str]:
    """Every input of the simulation part in plan order: the arms, the controls, part 0's flies (and the B-myr-alone
    flies when `bmyr`)."""
    kinds = FLY_KINDS + ((BESIDE_KIND,) if bmyr else ())
    return list(ARMS) + list(CONTROLS) + [fly_name(f, k) for k in kinds for f in range(N_FLIES)]


def input_of(name: str) -> pd.DataFrame:
    """An input by name ('I0', an arm, a control, 'F<f>_<kind>'), the lab input's layout."""
    c = _CACHE.setdefault('inputs', {})
    if name not in c:
        if name == 'I0':
            c[name] = lab_input()
        elif name in ARMS or name in CONTROLS:
            c[name] = E83.input_at(factors_of(name), 1.0, lab_input())
        elif parse_fly(name):
            c[name] = fly_input(*parse_fly(name))
        else:
            raise ValueError(f'unknown input {name!r}')
    return c[name].copy()


def input_path(name: str) -> Path:
    return out_dir() / 'inputs' / f'{name}.csv'


def input_spec(name: str) -> str:
    """The `--orn-input` of E10_joint_fit.py (relative to the shared directory)."""
    return str(input_path(name).relative_to(H.shared_dir()))


def changed_entries(name: str, sp: Optional[pd.DataFrame] = None) -> List[dict]:
    """The entries an input changes: odor, glomerulus, lab value, new value, log2 factor (where both are positive) and
    the log2 factor in the validation2 ORN flies' SDs (E109's `flies_spread`)."""
    lab = H.abbrev_columns(lab_input())
    new = H.abbrev_columns(input_of(name))
    sp = E109.flies_spread() if sp is None else sp
    sd = sp.set_index(['odor', 'glomerulus'])['flies_sd_log2']
    rows = []
    for o in lab.columns:
        diff = lab[o].to_numpy(float) != new[o].to_numpy(float)
        for g in lab.index[diff]:
            a, b = float(lab.loc[g, o]), float(new.loc[g, o])
            l2 = float(np.log2(b / a)) if (a > 0 and b > 0) else np.nan
            s = float(sd.get((o, g), np.nan))
            rows.append(dict(input=name, odor=o, glomerulus=g, lab=a, new=b, log2_factor=l2, flies_sd_log2=s,
                             in_flies_sd=l2 / s if np.isfinite(s) and s > 0 else np.nan))
    return rows


def write_inputs() -> dict:
    """Part inputs: inputs/<name>.csv (E38's layout) for I0 and every input, sha256.json, inputs.json (descriptions) and
    entries.csv (`changed_entries` of the arms and controls; per fly the entries of DM1 and DP1m)."""
    from experiments.E77_input_screen import frame_sha256
    sizes = entry_sizes(check=True)
    fly_deviations()
    lab = lab_input()
    if not np.array_equal(E83.input_at(factors_of('D1') * 0 + 1, 1.0, lab).to_numpy(), lab.to_numpy()):
        raise RuntimeError('factor 1 does not give the lab input bit for bit')
    shas, desc, rows = {}, {}, []
    sp = E109.flies_spread()
    for name in ['I0'] + input_names(bmyr=True):
        orn = input_of(name)
        p = input_path(name)
        p.parent.mkdir(parents=True, exist_ok=True)
        lay = H.abbrev_columns(orn)
        lay.index.name = 'glomerulus'
        p.write_text(lay.to_csv())
        shas[name] = dict(file=str(p.relative_to(H.shared_dir())) if H.shared_dir() in p.parents else str(p),
                          file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(), frame_sha256=frame_sha256(orn))
        ch = changed_entries(name, sp) if name != 'I0' else []
        if name in sizes:
            desc[name] = dict(entries={f'{o} {g}': v for (o, g), v in sizes[name].items()},
                              kind='arm' if name in ARMS else 'control')
            rows += ch
        elif parse_fly(name):
            f, k = parse_fly(name)
            desc[name] = dict(fly=fly_deviations()['keys'][f], kind=f'fly_{k}', odors=list(FLY_ODORS[k]),
                              n_changed=len(ch), source=("E108's B input (step 4)" if k == BESIDE_KIND else
                                                         "the lab's column plus the fly's deviation (E109's units)"))
            rows += [r for r in ch if r['glomerulus'] in ('DM1', 'DP1m')]
        else:
            desc[name] = dict(kind='reference')
    _write_json(out_dir() / 'inputs' / 'sha256.json', shas)
    _write_json(out_dir() / 'inputs.json', dict(generated=_now(), script=SCRIPT, inputs=desc, control_glomeruli=CONTROL_GLOMERULI,
                                                excluded=list(EXCLUDED), approx_conc_note=H.APPROX_CONC_NOTE,
                                                **H.code_commits()))
    pd.DataFrame(rows).to_csv(out_dir() / 'entries.csv', index=False)
    print(f'E113 inputs: {len(shas)} written; controls {CONTROL_GLOMERULI}', flush=True)
    return shas


# ---------------------------------------------------------------------------
# step 1: identity

def held(label: str, orn: pd.DataFrame, name: str) -> dict:
    """A held reading (E85's `held_reading` with this experiment's breadth cache)."""
    from experiments import E85_orn_fly_population as E85
    return E85.held_reading(label, orn, E83.input_tag_for(orn, f'E113_{name}'), cache_dir=local_dir())


def _prior_identity() -> Dict[str, dict]:
    """E109's and E116's part-0 held paths (part0.json), for the record."""
    out = {}
    for exp in ('E109_input_search_rho', 'E116_input_second_flies'):
        p = H.shared_dir() / 'results' / exp / 'part0.json'
        out[exp] = (_read_json(p).get('held_path') or _read_json(p).get('reading_path') or {}) if p.exists() else {}
    return out


def identity(wirings: Optional[Sequence[str]] = None) -> dict:
    """Step 1 (the module docstring); identity.json, merged across invocations. Raises if a held path differs."""
    from experiments import E108_paa_extrapolation as E108
    cc = commits_clean()
    p = out_dir() / 'identity.json'
    res = _read_json(p) or dict(held_path={})
    res.setdefault('held_path', {})
    lab = lab_input()
    ws = list(wirings or ALL_WIRINGS)
    for w in ws:
        t0 = time.time()
        hr = held(w, lab, 'I0')
        ns, nr = E83._same({'on': hr['on'], 'off': hr['off']}, E83.saved_fit(w)['saved'])
        res['held_path'][w] = dict(count_diffs=ns, response_diffs=nr, thr_factor=float(hr['thr_factor']),
                                   apl_factor=float(hr['apl_factor']), seconds=time.time() - t0, host=_host(), at=_now(), **cc)
        print(f"E113 identity {w}: {ns} count and {nr} response differences, factors {hr['thr_factor']} / "
              f"{hr['apl_factor']} ({time.time() - t0:.0f}s)", flush=True)
        del hr
        gc.collect()
    if 'real' in ws:
        t0 = time.time()
        ref = E108.load_run(E108.run_tag('B0', 'real'))
        if ref is None:
            res['bmyr_check'] = dict(found=False, exact=False, note="E108's B_B0_real not found locally or in the archive")
        else:
            hb = held('real', fly_input(0, BESIDE_KIND), 'F0_bmyr')
            ns, nr = E83._same({'on': hb['on'], 'off': hb['off']}, ref)
            res['bmyr_check'] = dict(found=True, count_diffs=ns, response_diffs=nr, exact=bool(ns == 0 and nr == 0),
                                     seconds=time.time() - t0, host=_host(), at=_now())
            del hb
        print(f"E113 identity: E108's B0 on real {res['bmyr_check']}", flush=True)
    res['prior'] = _prior_identity()
    hp = res['held_path']
    res['all_exact'] = bool(set(hp) >= set(ALL_WIRINGS) and all(
        d['count_diffs'] == 0 and d['response_diffs'] == 0 and d['thr_factor'] == 1.0 and d['apl_factor'] == 1.0
        for d in hp.values()))
    res['bmyr_reuse'] = bool((res.get('bmyr_check') or {}).get('exact'))
    res['script'], res['updated'] = SCRIPT, _now()
    _write_json(p, res)
    bad = [w for w in ws if hp[w]['count_diffs'] or hp[w]['response_diffs'] or hp[w]['thr_factor'] != 1.0
           or hp[w]['apl_factor'] != 1.0]
    if bad:
        raise RuntimeError(f'E113 identity fails on {bad} (identity.json); nothing else runs')
    return res


def identity_state() -> dict:
    return _read_json(out_dir() / 'identity.json')


def bmyr_reuse() -> bool:
    return bool(identity_state().get('bmyr_reuse'))


# ---------------------------------------------------------------------------
# the held readings

def plan(wirings: Optional[Sequence[str]] = None, bmyr: Optional[bool] = None) -> List[dict]:
    """The held readings of part sims, in order: the arms, the controls, part 0's flies (the B-myr-alone flies too when
    E108's are not reused), each on the eight wirings (or `wirings`)."""
    if bmyr is None:
        bmyr = not bmyr_reuse()
    ws = list(wirings or ALL_WIRINGS)
    out = []
    for ch in input_names(bmyr=bmyr):
        kind = 'arm' if ch in ARMS else 'control' if ch in CONTROLS else f'fly_{parse_fly(ch)[1]}'
        out += [dict(change=ch, wiring=w, tag=f'{ch}_{w}', kind=kind) for w in ws]
    return out


def lane_points(points: List[dict], lane: Optional[str]) -> List[dict]:
    """The points of lane 'i/n' (every n-th from the i-th, 0-based); all for None."""
    if not lane:
        return points
    i, n = (int(x) for x in lane.split('/'))
    if not 0 <= i < n:
        raise ValueError(f'lane {lane}')
    return points[i::n]


def _lane_label(lane: Optional[str]) -> str:
    return '0of1' if not lane else '{}of{}'.format(*lane.split('/'))


def _sim_rows() -> pd.DataFrame:
    fs = sorted(out_dir().glob('runs_sims_*.csv'))
    parts = [E83._read_rows(f) for f in fs]
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True, sort=False) if parts else pd.DataFrame()


def part_sims(lane: Optional[str] = None, wirings: Optional[Sequence[str]] = None,
              max_points: Optional[int] = None) -> pd.DataFrame:
    """Part sims (the module docstring): the lane's points not yet in any runs_sims_*.csv, wiring by wiring, each wiring's
    simulator first checked against its saved fit."""
    from experiments.E77_input_screen import frame_sha256
    if not identity_state().get('all_exact'):
        raise RuntimeError('run --part identity on all eight wirings first (identity.json all_exact)')
    cc = commits_clean()
    lab = _lane_label(lane)
    rp, pp = out_dir() / f'runs_sims_{lab}.csv', out_dir() / f'provenance_sims_{lab}.json'
    done = set(_sim_rows().get('tag', pd.Series(dtype=str)))
    todo = [p for p in lane_points(plan(wirings), lane) if p['tag'] not in done]
    if max_points is not None:
        todo = todo[:max_points]
    print(f'E113 sims lane {lab}: {len(todo)} points to run', flush=True)
    for w in [w for w in ALL_WIRINGS if any(p['wiring'] == w for p in todo)]:
        sim, info = E83.build_simulator(w)
        E83._check_saved(w, sim, info['params'], E83.saved_fit(w)['saved'])
        del sim, info
        gc.collect()
        for pt in [p for p in todo if p['wiring'] == w]:
            orn = input_of(pt['change'])
            path = input_path(pt['change'])
            if not path.exists():
                raise RuntimeError(f'{path} missing: run --part inputs first')
            res = held(w, orn, pt['change'])
            o = {'on': res['on'], 'off': res['off']}
            save_run('R', pt['tag'], o)
            row = E83._run_row('sims', pt['change'], w, pt['tag'], res, o, kind=pt['kind'], lane=lab, host=_host(),
                               input_tag=res['input_tag'],
                               orn_input=str(path.relative_to(H.shared_dir())) if H.shared_dir() in path.parents else str(path),
                               input_frame_sha256=frame_sha256(orn),
                               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}))
            row.update(E109.run_extras(o))
            row.update({k: v for k, v in cc.items() if k == 'queue_dirty_paths'})
            E83._append_row(rp, row)
            _provenance_add(pp, pt['tag'], H.model_provenance(orn_deltas=orn))
            print(f"E113 sims {pt['tag']}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/"
                  f"{res['target_on']} on, {res['n_off']}/{res['target_off']} off); flag '{row.get('realism_flag') or ''}'; "
                  f"{res['seconds']:.0f}s", flush=True)
            del res, o
            gc.collect()
    return E83._read_rows(rp)


# ---------------------------------------------------------------------------
# step 6 beside: Remy's validation2 somata

def part_remy() -> pd.DataFrame:
    """Step 6 (beside, no class): D1 and D2 on the validation2 ORN input (E116's `v2_with` at size 1), T as fitted (E109's
    `fitted_on`); runs_V.csv."""
    from experiments.E77_input_screen import frame_sha256
    cc = commits_clean()
    rp, pp = out_dir() / 'runs_V.csv', out_dir() / 'provenance_V.json'
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    for arm in ('D1', 'D2'):
        orn = E116.v2_with(factors_of(arm), 1.0)
        for w in T_WIRINGS:
            tag = f'V_{arm}_{w}'
            if tag in done:
                continue
            t0 = time.time()
            o = E109.fitted_on(w, orn)
            save_run('V', tag, {'on': o['on'], 'off': o['off']})
            E83._append_row(rp, dict(part='V', change=arm, wiring=w, tag=tag, panel='validation2', reading='fitted',
                                     input_frame_sha256=frame_sha256(orn), fit_folder=E83.saved_fit(w)['folder'],
                                     fit_tag=E83.saved_fit(w)['tag'], **{k: o['params'][k] for k in E83.PARAM_KEYS},
                                     rate_on=float(o['on'][0].to_numpy(dtype=float).mean()),
                                     rate_off=float(o['off'][0].to_numpy(dtype=float).mean()),
                                     seconds=time.time() - t0, host=_host(), **cc))
            _provenance_add(pp, tag, H.model_provenance(orn_deltas=orn))
            print(f'E113 V {tag}: {time.time() - t0:.0f}s', flush=True)
    return E83._read_rows(rp)


def remy_readings() -> pd.DataFrame:
    """remy_readings.csv: the family (E107's as-fitted validation2 runs, 'I0') and D1, D2 (runs_V.csv) on T's population
    and each wiring of T, APL on, against Remy's mean table and her halves a to d: ρ over the 136 and the 31 of paa and
    B-myr, the distance of paa from B-myr."""
    rv = E83._read_rows(out_dir() / 'runs_V.csv')
    tg = E116.remy_targets()
    changes = {'I0': E116.v2_fitted(T_WIRINGS)}
    for arm in ('D1', 'D2'):
        rs = {w: load_run('V', f'V_{arm}_{w}') for w in T_WIRINGS}
        if len(rv) and all(o is not None for o in rs.values()):
            changes[arm] = rs
    rows = []
    for ch, runs in changes.items():
        for unit, ws in [('T', T_WIRINGS)] + [(w, (w,)) for w in T_WIRINGS]:
            pop = E116.population_on({w: runs[w]['on'][1] for w in ws})
            for t in ('all',) + E116.REMY_HALVES:
                rows.append(dict(change=ch, unit=unit, target=t, rho=E116.remy_rho(pop, tg[t]),
                                 rho_pb=E116.remy_rho(pop, tg[t], 'pb'), paa_bmyr=E116.pb_distance(pop)))
    df = pd.DataFrame(rows)
    if len(df):
        ref = df[df.change == 'I0'].set_index(['unit', 'target'])
        df['d_rho'] = [r.rho - ref['rho'].get((r.unit, r.target), np.nan) for r in df.itertuples()]
        df['d_rho_pb'] = [r.rho_pb - ref['rho_pb'].get((r.unit, r.target), np.nan) for r in df.itertuples()]
    df.to_csv(out_dir() / 'remy_readings.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# step 5: the refit (only on a pass)

def fit_dirs(w: str) -> Tuple[Path, Path]:
    """(results, outputs) folders of the refit on `w`."""
    return H.shared_dir() / 'results' / FIT_PREFIX / f'refit_{w}', H.HERE / 'outputs' / FIT_PREFIX / f'refit_{w}'


def _fits_row(res_dir: Path) -> Optional[pd.Series]:
    p = res_dir / 'fits.csv'
    if not p.exists():
        return None
    f = pd.read_csv(p)
    return f.iloc[-1] if len(f) else None


def load_fit(w: str) -> Optional[Tuple[pd.Series, Dict[str, tuple]]]:
    """(fits.csv row, saved runs) of the refit on `w`: this host's outputs, else the archive's; or None."""
    rd, od = fit_dirs(w)
    row = _fits_row(rd)
    if row is None:
        return None
    for d in (od, archive_dir() / 'fits' / f'refit_{w}'):
        o = E83._load_pair(d, str(row['tag']))
        if o is not None:
            return row, o
    return None


def refit_name() -> Optional[str]:
    """The refit's input from score's refit.json (None: no pass)."""
    return _read_json(out_dir() / 'refit.json').get('input')


def write_starts(name: str) -> dict:
    """starts.json: R_<w>, the held point of `name` on each of T's wirings (its APL-on parameters, absolute)."""
    from experiments import E110_no_kckc as E110
    runs = _sim_rows()
    entries = {}
    for w in T_WIRINGS:
        r = runs[(runs['change'] == name) & (runs['wiring'] == w)]
        if not len(r):
            raise RuntimeError(f'{name} has no held point on {w}')
        r = r.iloc[-1]
        entries[f'R_{w}'] = E110.start_entry(json.loads(r['params_on']), input=name,
                                             source=f"the held point {r['tag']} (results/{EXPERIMENT}/runs_sims_*.csv)")
    return E110.merge_starts(out_dir() / 'starts.json', entries)


def fit_command(w: str, name: str) -> List[str]:
    """The refit's command on `w` (the module docstring)."""
    from experiments import E110_no_kckc as E110
    py = str(H.HERE / 'venv' / 'bin' / 'python')
    return ([py, '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3', '--rate-target', 'ctrl_frac', '--loss',
             'apl', '--apl-readout', 'both', '--tag', '_cap3_refit', '--maxfev', str(MAXFEV), '--experiment',
             f'{FIT_PREFIX}/refit_{w}'] + E110.wiring_args(w)
            + ['--orn-input', input_spec(name), '--init-params', f'results/{EXPERIMENT}/starts.json:R_{w}'])


def part_refit(dry_run: bool = False, wirings: Optional[Sequence[str]] = None) -> List[str]:
    """Step 5 (only when refit.json names an input): one E10_joint_fit.py process after another on T's wirings (skipping
    those done); runs_refit.csv."""
    from experiments import E110_no_kckc as E110
    name = refit_name()
    if not name:
        print('E113 refit: no input passed the gate (refit.json); no refit', flush=True)
        return []
    if dry_run:
        for w in (wirings or T_WIRINGS):
            print(' '.join(fit_command(w, name)), flush=True)
        return []
    write_starts(name)
    cc = E110._fit_clean()
    done = []
    for w in (wirings or T_WIRINGS):
        rd, _ = fit_dirs(w)
        if _fits_row(rd) is not None:
            done.append(w)
            continue
        cmd = fit_command(w, name)
        print(f'E113 refit {w}: ' + ' '.join(cmd), flush=True)
        rec = E110.launch(cmd, rd, local_dir('logs') / f'refit_{w}.log', cc, f'refit_{w}')
        rec['launcher'] = SCRIPT
        _write_json(rd / 'launch.json', rec)
        if rec['returncode']:
            raise RuntimeError(f"E113 refit {w}: E10_joint_fit.py exited {rec['returncode']} (log {rec['log']})")
        r = _fits_row(rd)
        if r is None:
            raise RuntimeError(f"E113 refit {w}: no row in its fits.csv (log {rec['log']})")
        E83._append_row(out_dir() / 'runs_refit.csv', dict(name=f'refit_{w}', wiring=w, input=name, tag=r['tag'],
                        folder=f'{FIT_PREFIX}/refit_{w}', host=_host(), seconds=rec['seconds'], command=' '.join(cmd),
                        **E110.params_of_row(r), loss=float(r['loss']), n_evals=int(r['n_evals']), **cc))
        print(f"E113 refit {w}: loss {float(r['loss']):.4f}, {rec['seconds'] / 60:.0f} min", flush=True)
        done.append(w)
    return done


# ---------------------------------------------------------------------------
# scoring (no simulations)

def _e108_bmyr_rows() -> List[dict]:
    """E108's B-myr-alone runs as F<f>_bmyr rows (when step 1 reproduced E108's B0 and the runs are there)."""
    from experiments import E108_paa_extrapolation as E108
    if not bmyr_reuse():
        return []
    rows = []
    for f in range(N_FLIES):
        os_ = {w: E108.load_run(E108.run_tag(f'B{f}', w)) for w in ALL_WIRINGS}
        if any(o is None for o in os_.values()):
            continue
        for w, o in os_.items():
            rows.append(dict(change=fly_name(f, BESIDE_KIND), wiring=w, tag=f'{fly_name(f, BESIDE_KIND)}_{w}',
                             kind='fly_bmyr', source=f"E108's {E108.run_tag(f'B{f}', w)}", **E83.realism(o),
                             **E109.run_extras(o)))
    return rows


def runs_table() -> pd.DataFrame:
    """Every run read: the saved fits as 'I0' (all eight), part sims' rows, E108's B-myr-alone runs when reused, the
    refits as 'refit' (T); with 'flagged'."""
    rows = [dict(change='I0', wiring=w, tag=f'I0_{w}', kind='reference', **E83.realism(E83.saved_fit(w)['saved']),
                 **E109.run_extras(E83.saved_fit(w)['saved'])) for w in ALL_WIRINGS]
    rows += _e108_bmyr_rows()
    for w in T_WIRINGS:
        f = load_fit(w)
        if f is not None:
            rows.append(dict(change=REFIT, wiring=w, tag=f'{REFIT}_{w}', kind='refit', fit_tag=str(f[0]['tag']),
                             **E83.realism(f[1]), **E109.run_extras(f[1])))
    parts = [pd.DataFrame(rows), _sim_rows()]
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates('tag', keep='last').reset_index(drop=True)
    df['flagged'] = [E109._flagged(x) for x in df['realism_flag']]
    return df


def run_outputs(change: str, wiring: str) -> Dict[str, tuple]:
    if change == 'I0':
        return E83.saved_fit(wiring)['saved']
    if change == REFIT:
        f = load_fit(wiring)
        if f is None:
            raise RuntimeError(f'no refit on {wiring}')
        return f[1]
    o = load_run('R', f'{change}_{wiring}')
    if o is None and parse_fly(change) and parse_fly(change)[1] == BESIDE_KIND:
        from experiments import E108_paa_extrapolation as E108
        o = E108.load_run(E108.run_tag(f'B{parse_fly(change)[0]}', wiring))
    if o is None:
        raise RuntimeError(f'no saved run {change}_{wiring}')
    return o


def _dist(change: str, wiring: str, readout: str) -> Dict[str, pd.Series]:
    key = ('dist', change, wiring, readout)
    if key not in _CACHE:
        _CACHE[key] = E83.pair_distances(run_outputs(change, wiring), readout)
    return _CACHE[key]


def readings(runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings.csv: per change, unit (each wiring; 'all8', 'S' and 'T' for a change on all eight; 'T' for one on T) and
    readout ('raw', 'cal'), E108's `unit_scores`, the mean realism quantities over the unit's runs, the runs flagged and
    the distance of paa from B-myr."""
    from experiments import E108_paa_extrapolation as E108
    runs = runs_table() if runs is None else runs
    rows = []
    for ch, g in runs.groupby('change', sort=False):
        ws = set(g['wiring'])
        if ws >= set(ALL_WIRINGS):
            wirings, pops = ALL_WIRINGS, (('all8', ALL_WIRINGS), ('S', S_WIRINGS), ('T', T_WIRINGS))
        elif ws >= set(T_WIRINGS):
            wirings, pops = T_WIRINGS, (('T', T_WIRINGS),)
        else:
            continue
        meta = {r.wiring: r for r in g.itertuples()}
        for ro in ('raw', 'cal'):
            ds = {w: _dist(ch, w, ro) for w in wirings}
            for unit, uw in [(w, (w,)) for w in wirings] + list(pops):
                d = ds[unit] if unit in ds else E83.population([ds[w] for w in uw])
                sub = g[g['wiring'].isin(uw)]
                pb = E83._cs(d['on']).index.get_loc(tuple(sorted(PAA_BMYR)))
                rows.append(dict(change=ch, unit=unit, kind='wiring' if unit in ds else 'population', readout=ro,
                                 n_flagged=int(sum(meta[w].flagged for w in uw)),
                                 **{f'run_{k}': float(sub[k].astype(float).mean()) for k in E109.RUN_MEANS if k in sub},
                                 paa_bmyr_on=float(E83._cs(d['on']).iloc[pb]), paa_bmyr_off=float(E83._cs(d['off']).iloc[pb]),
                                 **E108.unit_scores(d)))
    rd = pd.DataFrame(rows)
    rd.to_csv(out_dir() / 'readings.csv', index=False)
    return rd


def _have(rd: pd.DataFrame, change: str, unit: str) -> bool:
    return bool(((rd['change'] == change) & (rd['unit'] == unit)).any())


def gate_arm(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, controls: Sequence[str] = CONTROLS) -> dict:
    """One arm's gate on T (the module docstring): the gain against all recordings, A and B on T and on 3 of its 4
    wirings, the other ρ, each control below half the gain, the other measures (E116's `_measures`); the eight wirings
    beside. 'cls' is E109's `class_of` without the pending items (columns only here); 'class_with_columns' with them."""
    out = dict(arm=arm, complete=False)
    v = E116.change_values(rd, runs, arm, 'T', T_WIRINGS, T_GAIN_WIRINGS)
    cv = {c: E116.change_values(rd, runs, c, 'T', T_WIRINGS, T_GAIN_WIRINGS) for c in controls}
    if v is None or any(x is None for x in cv.values()):
        out['note'] = 'readings missing: ' + ' '.join([arm] * (v is None) + [c for c, x in cv.items() if x is None])
        return out
    out['complete'] = True
    out['directions'] = {arm: v}
    for c in ('on', 'off'):
        ok, other = E116._gain([v], c)
        out[f'gain_{c}_parts'] = dict(gain=ok, other=other)
        out[f'gain_{c}'] = bool(ok and other)
        out[f'mean_gain_{c}'] = float(v[f'd_rho_{c}_all'])
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    ctl = {}
    for c in ('on', 'off'):
        g0 = out[f'mean_gain_{c}']
        each = {k: float(x[f'd_rho_{c}_all']) for k, x in cv.items()}
        ctl[c] = dict(candidate=g0, **each, max=float(max(each.values())),
                      ok=bool(g0 > 0 and all(x < CONTROL_RATIO * g0 for x in each.values())))
    out['controls'] = ctl
    E116._measures([v], out)
    out['cls'] = E109.class_of(out, with_pending=False)
    out['class_with_columns'] = E109.class_of(out)
    b = E116.change_values(rd, runs, arm, 'all8', ALL_WIRINGS, ALL8_GAIN_WIRINGS) if _have(rd, arm, 'all8') else None
    if b is not None:
        out['all8'] = {k: b[k] for k in ('d_rho_on_all', 'd_rho_off_all', 'd_rho_on_A', 'd_rho_off_A', 'd_rho_on_B',
                                         'd_rho_off_B', 'n_wirings_gain_on_all', 'n_wirings_gain_off_all',
                                         'gain_wirings_ok_on', 'gain_wirings_ok_off', 'd_pattern_cal', 'n_flagged')}
    return out


def fly_picks(rd: pd.DataFrame, kind: str) -> Optional[Dict[str, str]]:
    """{'A': the fly input with the highest ρ with APL on of S's population against half A, 'B': against half B}, or None
    if a fly's S reading is missing."""
    names = [fly_name(f, kind) for f in range(N_FLIES)]
    if not all(_have(rd, n, 'S') for n in names):
        return None
    out = {}
    for X in UNSEEN:
        vals = [E116._v(rd, n, 'S', f'rho_on_{X}') for n in names]
        out[X] = names[int(np.nanargmax(vals))]
    return out


def gate_flies(rd: pd.DataFrame, runs: pd.DataFrame, kind: str) -> dict:
    """Part 0's rule for one kind (the module docstring): each direction's pick read on T against the other half."""
    out = dict(kind=kind, complete=False)
    picks = fly_picks(rd, kind)
    if picks is None:
        out['note'] = 'S readings missing'
        return out
    dirs = {}
    for X, ch in picks.items():
        Y = UNSEEN[X]
        v = E116.change_values(rd, runs, ch, 'T', T_WIRINGS, T_GAIN_WIRINGS)
        if v is None:
            out['note'] = f'T readings of {ch} missing'
            return out
        d = dict(v, picked_on=X, read_on=Y, pick=ch, rho_on_S_picked=E116._v(rd, ch, 'S', f'rho_on_{X}'))
        for c in ('on', 'off'):
            o = 'off' if c == 'on' else 'on'
            d[f'd_rho_{c}_held'] = float(v[f'd_rho_{c}_{Y}'])
            d[f'n_wirings_gain_{c}_held'] = int(sum(x >= GAIN[c] for x in v[f'd_rho_{c}_{Y}_wirings'].values()))
            d[f'gain_{c}_held'] = bool(d[f'd_rho_{c}_held'] >= GAIN[c] and d[f'n_wirings_gain_{c}_held'] >= T_GAIN_WIRINGS
                                       and v[f'd_rho_{o}_{Y}'] >= -OTHER_SLACK)
        dirs[f'{X}to{Y}'] = d
    out['complete'] = True
    out['picks'] = picks
    out['directions'] = dirs
    for c in ('on', 'off'):
        out[f'gain_{c}'] = bool(all(d[f'gain_{c}_held'] for d in dirs.values()))
        out[f'mean_gain_{c}'] = float(np.mean([d[f'd_rho_{c}_held'] for d in dirs.values()]))
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    out['controls'] = None
    E116._measures(list(dirs.values()), out)
    out['cls'] = E109.class_of(out, with_pending=False)
    out['class_with_columns'] = E109.class_of(out)
    return out


def _passes(g: Mapping) -> bool:
    return bool(g.get('complete') and g.get('cls') == 'passes')


def outcomes(arms: Mapping[str, dict], flies: Mapping[str, dict]) -> List[str]:
    """The outcome letters (the module docstring) from the arms' gates and part 0's (FLY_KINDS only)."""
    out = []
    if _passes(arms.get(CLASS_ARM, {})):
        out.append('a')
    elif any(_passes(arms.get(a, {})) for a in ARMS if a != CLASS_ARM):
        out.append('b')
    if any(_passes(flies.get(k, {})) for k in FLY_KINDS):
        out.append('c')
    if out:
        return out
    gs = [g for g in list(arms.values()) + [flies.get(k, {}) for k in FLY_KINDS] if g.get('complete')]
    if not gs:
        return out
    carries = any(g.get('cls') == 'carries below the gate' for g in gs)
    ctl_as_well = all(all(g['controls'][c]['max'] >= g['controls'][c]['candidate'] for c in ('on', 'off'))
                      for g in arms.values() if g.get('complete'))
    out.append('d' if carries and not ctl_as_well else 'e')
    return out


def refit_choice(arms: Mapping[str, dict], flies: Mapping[str, dict], rd: pd.DataFrame) -> Optional[str]:
    """The refit's input (the module docstring), or None."""
    for a in ('D2', 'D1', 'D1x'):
        if _passes(arms.get(a, {})):
            return a
    ok = [k for k in FLY_KINDS if _passes(flies.get(k, {}))]
    if not ok:
        return None
    gain = lambda g: sum(g[f'mean_gain_{c}'] for c in ('on', 'off') if g.get(f'gain_{c}'))
    k = max(ok, key=lambda x: gain(flies[x]))
    picks = list(dict.fromkeys(flies[k]['picks'].values()))
    if len(picks) == 1:
        return picks[0]
    return max(picks, key=lambda n: E116._v(rd, n, 'S', 'rho_on_all'))


def refit_reading(rd: pd.DataFrame, runs: pd.DataFrame) -> dict:
    """The refit, beside (no gate): E116's values on T against I0, its class read the same way without controls, and the
    parameters against the family's fits (moves beyond MOVE_NAMED named)."""
    out = dict(complete=False, input=refit_name())
    v = E116.change_values(rd, runs, REFIT, 'T', T_WIRINGS, T_GAIN_WIRINGS) if _have(rd, REFIT, 'T') else None
    if v is None:
        return out
    out.update(complete=True, directions={REFIT: v})
    for c in ('on', 'off'):
        ok, other = E116._gain([v], c)
        out[f'gain_{c}'] = bool(ok and other)
        out[f'mean_gain_{c}'] = float(v[f'd_rho_{c}_all'])
    out['controls'] = None
    E116._measures([v], out)
    out['cls'] = E109.class_of(out, with_pending=False)
    rows, named = [], []
    for w in T_WIRINGS:
        f = load_fit(w)
        if f is None:
            continue
        from experiments import E110_no_kckc as E110
        p, p0 = E110.params_of_row(f[0]), E83.saved_fit(w)['p']
        for k in E83.PARAM_KEYS:
            rel = float(p[k] / p0[k] - 1) if p0[k] else np.nan
            rows.append(dict(wiring=w, param=k, family=float(p0[k]), refit=float(p[k]), rel_change=rel))
            if k != 'wKCKC_scale' and np.isfinite(rel) and abs(rel) > MOVE_NAMED:
                named.append(f'{w} {k} {rel:+.0%}')
    pd.DataFrame(rows).to_csv(out_dir() / 'refit_params.csv', index=False)
    out['moved_beyond_20pct'] = named
    return out


def e112_ended() -> Tuple[bool, str]:
    return E116.e112_ended()


def e112_table(runs: pd.DataFrame) -> pd.DataFrame:
    """Beside, only once E112 has ended: each arm's and control's population (all8, T) under E112's readings R0 to R3."""
    from experiments import E112_rho_like_for_like as E112
    ended, state = e112_ended()
    if not ended:
        return pd.DataFrame()
    rows = []
    for ch in ['I0'] + list(ARMS) + list(CONTROLS):
        if set(runs[runs['change'] == ch]['wiring']) < set(ALL_WIRINGS):
            continue
        per = {w: E112.run_readings(run_outputs(ch, w), E112.COUNT_READINGS, seed=E112.wiring_seed(w)) for w in ALL_WIRINGS}
        for unit, ws in (('all8', ALL_WIRINGS), ('T', T_WIRINGS)):
            for r in E112.COUNT_READINGS:
                pop = E112.population([per[w][r] for w in ws])
                rows.append(dict(change=ch, unit=unit, reading=r, label=E112.READING_LABEL.get(r, r), e112_state=state,
                                 **E112.reading_scores(pop)))
    df = pd.DataFrame(rows)
    if len(df):
        base = df[df['change'] == 'I0'].set_index(['unit', 'reading'])
        for c in ('rho_on_all120', 'rho_off_all120'):
            if c in df:
                df[f'd_{c}'] = [v - base[c].get((u, r), np.nan) for u, r, v in zip(df['unit'], df['reading'], df[c])]
    df.to_csv(out_dir() / 'e112.csv', index=False)
    return df


def flies_table(rd: pd.DataFrame) -> pd.DataFrame:
    """flies.csv: per fly input (part 0 and the B-myr-alone flies), ρ with APL on and off of S's population against A, B
    and all recordings, and T's change from I0 against all, A and B, with paa's distance from B-myr."""
    rows = []
    for k in FLY_KINDS + (BESIDE_KIND,):
        for f in range(N_FLIES):
            n = fly_name(f, k)
            if not _have(rd, n, 'T'):
                continue
            row = dict(input=n, kind=k, fly=f)
            for c in ('on', 'off'):
                for t in YANG_TARGETS:
                    row[f'S_rho_{c}_{t}'] = E116._v(rd, n, 'S', f'rho_{c}_{t}')
                    row[f'T_d_rho_{c}_{t}'] = E116._v(rd, n, 'T', f'rho_{c}_{t}') - E116._v(rd, 'I0', 'T', f'rho_{c}_{t}')
            row['T_d_paa_bmyr_on'] = E116._v(rd, n, 'T', 'paa_bmyr_on') - E116._v(rd, 'I0', 'T', 'paa_bmyr_on')
            row['T_d_pattern_cal'] = E116._v(rd, n, 'T', 'all_pattern', 'cal') - E116._v(rd, 'I0', 'T', 'all_pattern', 'cal')
            row['n_flagged_all8'] = int(E116._v(rd, n, 'all8', 'n_flagged')) if _have(rd, n, 'all8') else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'flies.csv', index=False)
    return df


def _f(x, spec: str = '.3f') -> str:
    return E109._f(x, spec)


def _md(df: pd.DataFrame, spec: str = '.4f') -> str:
    return E109._md_table(df, spec)


def _value_row(label: str, v: Mapping, cls: str = '') -> str:
    t = lambda k, s='+.4f': _f(v.get(k), s)
    return (f"| {label} | {' / '.join(t(f'd_rho_on_{x}') for x in YANG_TARGETS)} | "
            f"{' / '.join(t(f'd_rho_off_{x}') for x in YANG_TARGETS)} | "
            f"{v.get('n_wirings_gain_on_all', '')} / {v.get('n_wirings_gain_off_all', '')} | {t('d_pattern_cal')} | "
            f"{t('d_amplitude_cal')} | {t('d_shape_cal')} | {t('d_sd_ratio_cal', '+.3f')} | {t('d_pearson_effect_cal', '+.3f')} | "
            f"{t('d_loss_raw')} | {t('d_mean_effect_cal')} | {t('d_rate_rho_on', '+.3f')} / {t('d_rate_rho_off', '+.3f')} | "
            f"{v.get('n_flagged', '')} | {t('d_paa_bmyr_on', '+.3f')} / {t('d_paa_bmyr_off', '+.3f')} | "
            f"{t('d_block', '+.3f')} | {cls} |")


VALUE_HEADER = ('| input | Δρ on, all / A / B | Δρ off, all / A / B | T wirings gaining on / off (all) | pattern part (cal) | '
                'amplitude part | shape part | spread ratio | Pearson r on the per-pair APL effect | raw-count loss | '
                'mean APL effect (cal) | response rate ρ on / off | runs flagged | paa from B-myr on / off | the block (column) | '
                'class |\n|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---:|---|')


def comparison_table(rd: pd.DataFrame, runs: pd.DataFrame, arms: Dict[str, dict], flies: Dict[str, dict],
                     ctl_vals: Dict[str, dict], letters: List[str], refit: Optional[str], rr: dict,
                     remy: pd.DataFrame, e112: object) -> str:
    """comparison_table.md: the tables the note reads first."""
    idn = identity_state()
    L = ['# E113: the input of paa and B-myr without a search', '',
         "Changes are from the family at the lab's input (I0) on the same wirings, raw counts unless said; 'all', A, B: "
         "Yang's recordings (all, and E77's halves). Population of T (P3, P4, M2, M3) unless said.", '',
         f"Step 1: held reading of the lab's input reproduces the saved fits: all_exact {idn.get('all_exact')}; E108's B0 "
         f"on real reproduced: {idn.get('bmyr_reuse')} (B-myr-alone runs {'reused from E108' if idn.get('bmyr_reuse') else 'run here'}).", '',
         f"Family on T against all recordings: ρ {_f(E116._v(rd, 'I0', 'T', 'rho_on_all'))} / "
         f"{_f(E116._v(rd, 'I0', 'T', 'rho_off_all'))}; against A {_f(E116._v(rd, 'I0', 'T', 'rho_on_A'))} / "
         f"{_f(E116._v(rd, 'I0', 'T', 'rho_off_A'))}; against B {_f(E116._v(rd, 'I0', 'T', 'rho_on_B'))} / "
         f"{_f(E116._v(rd, 'I0', 'T', 'rho_off_B'))}; pattern part (cal) {_f(E116._v(rd, 'I0', 'T', 'all_pattern', 'cal'), '.4f')}; "
         f"paa from B-myr {_f(E116._v(rd, 'I0', 'T', 'paa_bmyr_on'))} / {_f(E116._v(rd, 'I0', 'T', 'paa_bmyr_off'))}.", '',
         '## Arm D and its controls (T)', '', VALUE_HEADER]
    for name in ARMS + CONTROLS:
        v = arms[name]['directions'][name] if arms.get(name, {}).get('complete') else ctl_vals.get(name)
        if v is not None:
            L.append(_value_row(name, v, arms.get(name, {}).get('cls', 'control' if name in CONTROLS else '')))
    L += ['', f"{CLASS_ARM} carries the class; D1 and D1x are reported with the caveat that their sizes came from searches "
          "on Yang's recordings.", '', '| arm | ρ | candidate (all recordings) | ' + ' | '.join(CONTROLS)
          + ' | controls hold |', '|---|---|---:|' + '---:|' * len(CONTROLS) + '---|']
    for a in ARMS:
        g = arms.get(a, {})
        if not g.get('complete'):
            L.append(f"| {a} | | {g.get('note', 'incomplete')} |" + ' |' * (len(CONTROLS) + 1))
            continue
        for c in ('on', 'off'):
            k = g['controls'][c]
            L.append(f"| {a} | APL {c} | {_f(k['candidate'], '+.4f')} | " + ' | '.join(_f(k[x], '+.4f') for x in CONTROLS)
                     + f" | {'yes' if k['ok'] else 'no'} |")
    L += ['', '| arm | class | class with the columns | gain on / off | fails | columns failing | all eight: Δρ on / off (all), '
          'wirings gaining |', '|---|---|---|---|---|---|---|']
    for a in ARMS:
        g = arms.get(a, {})
        if g.get('complete'):
            b = g.get('all8') or {}
            L.append(f"| {a} | {g['cls']} | {g['class_with_columns']} | {g['gain_on']} / {g['gain_off']} | "
                     f"{g['fails'] or 'none'} | {g['pending_fails'] or 'none'} | {_f(b.get('d_rho_on_all'), '+.4f')} / "
                     f"{_f(b.get('d_rho_off_all'), '+.4f')}, {b.get('n_wirings_gain_on_all', '')} / "
                     f"{b.get('n_wirings_gain_off_all', '')} of 8 |")
    ent = out_dir() / 'entries.csv'
    if ent.exists():
        e = pd.read_csv(ent)
        e = e[e['input'].isin(ARMS + CONTROLS)]
        if len(e):
            L += ['', "## The entries changed, in the validation2 ORN flies' SDs", '',
                  _md(e[['input', 'odor', 'glomerulus', 'lab', 'new', 'log2_factor', 'flies_sd_log2', 'in_flies_sd']], '.3f')]
    L += ['', "## Part 0: each ORN fly's own input", '']
    ft = flies_table(rd)
    if len(ft):
        L.append(_md(ft[['input', 'S_rho_on_A', 'S_rho_on_B', 'S_rho_on_all', 'T_d_rho_on_all', 'T_d_rho_off_all',
                         'T_d_rho_on_A', 'T_d_rho_on_B', 'T_d_rho_off_A', 'T_d_rho_off_B', 'T_d_paa_bmyr_on',
                         'T_d_pattern_cal', 'n_flagged_all8']], '+.4f'))
    L += ['', '| kind | direction | pick | Δρ on / off against the held-out half (T) | T wirings gaining on / off | '
          'passes on / off |', '|---|---|---|---|---|---|']
    for k in FLY_KINDS + (BESIDE_KIND,):
        g = flies.get(k, {})
        if not g.get('complete'):
            L.append(f"| {k} | | {g.get('note', 'incomplete')} | | | |")
            continue
        for dn, d in g['directions'].items():
            L.append(f"| {k} | picked on {d['picked_on']}, read on {d['read_on']} | {d['pick']} | "
                     f"{_f(d['d_rho_on_held'], '+.4f')} / {_f(d['d_rho_off_held'], '+.4f')} | {d['n_wirings_gain_on_held']} / "
                     f"{d['n_wirings_gain_off_held']} | {d['gain_on_held']} / {d['gain_off_held']} |")
        L.append(f"| {k} | class | {g['cls']}{' (beside: E108 step 4, no outcome)' if k == BESIDE_KIND else ''} | "
                 f"fails {g['fails'] or 'none'} | | |")
    L += ['', f"## Outcome: {', '.join('(' + x + ')' for x in letters) if letters else 'incomplete'}; refit input: "
          f"{refit or 'none (no pass)'}", '']
    if rr.get('complete'):
        L += ['## The refit (beside, no gate)', '', VALUE_HEADER, _value_row(f"refit ({rr.get('input')})", rr['directions'][REFIT],
              rr.get('cls', '')), '', f"Parameters moved by more than 20 %: {', '.join(rr['moved_beyond_20pct']) or 'none'} "
              "(refit_params.csv).", '']
    if len(remy):
        sub = remy[(remy['unit'] == 'T')]
        L += ["## Beside: D1 and D2 on the validation2 ORN input against Remy's somata (T as fitted, APL on)", '',
              _md(sub[['change', 'target', 'rho', 'd_rho', 'rho_pb', 'd_rho_pb', 'paa_bmyr']], '.4f'), '']
    if isinstance(e112, pd.DataFrame) and len(e112):
        L += ["## Beside: E112's readings", '', _md(e112[[c for c in e112.columns if c in (
            'change', 'unit', 'reading', 'rho_on_all120', 'rho_off_all120', 'd_rho_on_all120', 'd_rho_off_all120')]], '.4f')]
    else:
        ended, state = e112_ended()
        L.append(f"E112's reading: E112 had not ended at scoring (its index row: '{state}'), so the two ρ's are not given "
                 'under its reading.' if not ended else "E112's reading: see e112.csv.")
    return '\n'.join(L) + '\n'


def realism_table(runs: pd.DataFrame) -> pd.DataFrame:
    """realism.csv: per change (the family on all eight and on T, every input), the mean and range over its runs of the
    response rates, their ratio, response rate ρ, the odors-per-KC distances, the largest counts and the factors, and
    the runs flagged."""
    rows = []
    for ch, g in runs.groupby('change', sort=False):
        units = [('all8', g), ('T', g[g['wiring'].isin(T_WIRINGS)])] if ch == 'I0' else [('runs', g)]
        for unit, s in units:
            row = dict(change=ch, unit=unit, n_runs=len(s), n_flagged=int(s['flagged'].sum()),
                       flags='; '.join(sorted({str(x) for x in s['realism_flag'] if E109._flagged(x)})))
            for k in ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv',
                      'max_count_on', 'max_count_off', 'thr_factor', 'apl_factor'):
                if k in s:
                    x = s[k].astype(float)
                    row.update({k: float(x.mean()), f'{k}_min': float(x.min()), f'{k}_max': float(x.max())})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'realism.csv', index=False)
    return df


def score() -> dict:
    """Part score (no simulations; the module docstring)."""
    out = out_dir()
    runs = runs_table()
    runs.to_csv(out / 'runs.csv', index=False)
    rd = readings(runs)
    arms = {a: gate_arm(rd, runs, a) for a in ARMS}
    ctl_vals = {c: E116.change_values(rd, runs, c, 'T', T_WIRINGS, T_GAIN_WIRINGS) for c in CONTROLS}
    pd.DataFrame([dict(control=c, **{k: v for k, v in x.items() if not isinstance(v, (dict, list))})
                  for c, x in ctl_vals.items() if x is not None]).to_csv(out / 'controls.csv', index=False)
    flies = {k: gate_flies(rd, runs, k) for k in FLY_KINDS + (BESIDE_KIND,)}
    letters = outcomes(arms, flies)
    complete = all(g.get('complete') for g in arms.values()) and all(flies[k].get('complete') for k in FLY_KINDS)
    refit = refit_choice(arms, flies, rd) if complete else None
    _write_json(out / 'refit.json', dict(input=refit, complete=complete, decided=_now(),
                                         rule='D2, D1, D1x, then the passing kind of part 0 (refit_choice)'))
    rr = refit_reading(rd, runs)
    rows = [E109._flat(g, which='arm', name=a, carrying=bool(a == CLASS_ARM)) for a, g in arms.items()]
    rows += [E109._flat(g, which='flies', name=k, carrying=bool(k in FLY_KINDS)) for k, g in flies.items()]
    if rr.get('complete'):
        rows.append(E109._flat(rr, which='refit', name=REFIT, carrying=False))
    pd.DataFrame(rows).to_csv(out / 'gate.csv', index=False)
    remy = remy_readings() if (out / 'runs_V.csv').exists() else pd.DataFrame()
    e112 = e112_table(runs)
    verdict = dict(arms={a: {k: g.get(k) for k in ('complete', 'cls', 'class_with_columns', 'gain_on', 'gain_off',
                                                   'mean_gain_on', 'mean_gain_off', 'controls', 'fails', 'pending_fails',
                                                   'note', 'all8')} for a, g in arms.items()},
                   flies={k: {x: g.get(x) for x in ('complete', 'cls', 'picks', 'gain_on', 'gain_off', 'mean_gain_on',
                                                    'mean_gain_off', 'fails', 'note')} for k, g in flies.items()},
                   outcomes=letters, complete=complete, refit_input=refit,
                   refit={k: rr.get(k) for k in ('complete', 'input', 'cls', 'gain_on', 'gain_off', 'mean_gain_on',
                                                 'mean_gain_off', 'fails', 'moved_beyond_20pct')},
                   e112=dict(zip(('ended', 'state'), e112_ended())))
    _write_json(out / 'gate.json', dict(verdict=verdict, gates=dict(arms=arms, flies=flies, refit=rr), **H.code_commits()))
    realism_table(runs)
    (out / 'comparison_table.md').write_text(comparison_table(rd, runs, arms, flies, ctl_vals, letters, refit, rr, remy, e112))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), script=SCRIPT, host=_host(), **H.code_commits(),
                                                    approx_conc_note=H.APPROX_CONC_NOTE))
    print(json.dumps(verdict, indent=1, default=str), flush=True)
    return verdict


def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive (part sims' runs, R; Remy's panel, V; the refits) and
    model_pair_distances_long.csv (raw spike counts, Yang's panel: R and the refits) for kcq/apl_effect.py."""
    out = out_dir()
    rows, parts = [], []
    for r in _sim_rows().to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=r['tag'], part='R', change=r['change'], kind=r.get('kind'), size=1.0,
                         wiring=r['wiring'], base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'], orn_input=r.get('orn_input'),
                         **{k: r[k] for k in E83.PARAM_KEYS}, kc_kc_inh=r['wKCKC_scale'], thr_factor=r['thr_factor'],
                         apl_factor=r['apl_factor'], rate_on=r['rate_on'], rate_off=r['rate_off'],
                         realism_flag=r.get('realism_flag'), subdir='R', panel='APL16',
                         kind_text=(f"C2's family's {r['wiring']} fit with input {r['change']}; breadth recomputed for the "
                                    "input; the somatic threshold and the tied APL scale re-held at the fit's APL-off and "
                                    "APL-on counts")))
    for r in E83._read_rows(out / 'runs_V.csv').to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=r['tag'], part='V', change=r['change'], kind='remy_panel', size=1.0,
                         wiring=r['wiring'], base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'], **{k: r[k] for k in E83.PARAM_KEYS},
                         kc_kc_inh=r['wKCKC_scale'], thr_factor=1.0, apl_factor=1.0, rate_on=r['rate_on'],
                         rate_off=r['rate_off'], subdir='V', panel='validation2',
                         kind_text=(f"C2's family's {r['wiring']} fit on the validation2 panel as fitted (nothing re-held), "
                                    f"arm {r['change']}'s factors on paa and B-myr")))
    for r in E83._read_rows(out / 'runs_refit.csv').to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=f"refit_{r['wiring']}", part='fits', change=r['input'], kind='refit',
                         wiring=r['wiring'], base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         ref_experiment=r['folder'], ref_fit=r['tag'], fixed_thr=r['fixed_thr'],
                         wAPLKC_scale=r['wAPLKC_scale'], wKCAPL_scale=r['wKCAPL_scale'], wAPLPN_scale=r['wAPLPN_scale'],
                         wPNAPL_scale=r['wPNAPL_scale'], wKCKC_scale=r['wKCKC_scale'], kc_kc_inh=r['wKCKC_scale'],
                         subdir=f"fits/refit_{r['wiring']}", panel='APL16',
                         kind_text=f"C2's family fitted again on {r['wiring']} with input {r['input']} (E72's part-B options)"))
    for r in rows:
        if r['panel'] != 'APL16':
            continue
        o = load_fit(r['wiring'])[1] if r['part'] == 'fits' and load_fit(r['wiring']) else load_run(r['subdir'], r['tag'])
        if o is None:
            continue
        parts.append(pd.concat({(r['variant'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['variant'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('inputs', 'identity', 'sims', 'remy', 'score', 'refit', 'variants'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(ALL_WIRINGS))
    ap.add_argument('--lane', default=None, help="part sims: 'i/n', every n-th point of the plan from the i-th (0-based)")
    ap.add_argument('--max-points', type=int, default=None)
    ap.add_argument('--dry-run', action='store_true', help='part refit: print the commands')
    a = ap.parse_args()
    if a.part == 'inputs':
        write_inputs()
    elif a.part == 'identity':
        identity(a.wirings)
    elif a.part == 'sims':
        part_sims(a.lane, a.wirings, a.max_points)
    elif a.part == 'remy':
        part_remy()
    elif a.part == 'score':
        score()
    elif a.part == 'refit':
        part_refit(a.dry_run, a.wirings)
    elif a.part == 'variants':
        write_variants()


if __name__ == '__main__':
    main()
