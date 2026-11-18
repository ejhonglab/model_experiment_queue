#!/usr/bin/env python3
"""E109: the input of paa, B-myr and 2-but searched for the two ρ's on four wirings and one half of the recordings, read
on the other four wirings and the other half.

The entry is item E109 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08; moved into the list 2026-11-17 23:56 under the idle-host rule,
planned by a Fable subagent, claude-fable-5-2, spawned from viki 23:38; claimed by viki 23:56). "The family" is C2's (the
breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of E72 part B's eight wirings at its own fit (E83's
`saved_fit`; S, searched: real, P1, P2, M1; T, held out: P3, P4, M2, M3), read as the mean matrix over the wirings named
(E72's population) on raw counts (cosine distance on spike counts). "The two ρ's" are ρ on odor-odor distances with APL
on against the control flies and with APL off against APL output blocked (TNT), each KC matrix the mean over the
recordings of a half (cosine; E108's `dist_targets`: halves A and B are E77's `recording_halves(0)`, C and D its
`recording_halves(1)`, E86's second split), through E106's `two_rhos` / `rho_pairs`.

Concentrations. paa's input is the lab's recording at 10^-4 for the flies' 10^-3 (or, if E108's class is "better",
E108's extrapolation to 10^-3, which is not a recording); euc and -aPine stay at the recorded 10^-1.75 and 10^-1.5
against the flies' 10^-2; B-myr's and 2-but's input are at the flies' concentrations. Every output that involves paa,
euc or -aPine is affected by this substitution (`H.APPROX_CONC_NOTE` in every provenance file).

Parts (CLI `--part`), in the entry's order of work:
- 0 (simulations, about 30 min on one 16-thread lane): the identities, the bound, the references (`part0`):
  - the objective at identity: S's population at the saved fits, the two ρ's against each half A to D and the pattern
    part under the calibrated readout (each wiring's K at its fit's) against each half; against A and B the pattern part
    must be E83's IDENTITY_OBJECTIVE (0.0569, 0.0552);
  - the search's path: on each of S's wirings, the pseudo-odor run (E83's InputPool, E64's `candidate_orn`) at factor 1
    reproduces the saved APL-on and APL-off counts and responses of each of the five odors searched (0 differences);
  - the reading's path: on each of T's wirings, the lab's input with its breadth recomputed (E85's `build`, this
    experiment's cache) gives the wiring's breadth and reproduces the saved runs, and the held reading of the lab's
    input (E85's `held_reading`) reproduces them with both factors 1 (0 differences);
  - step 7's path: on each of T's wirings, the validation2 panel as fitted (`fitted_on`, E107's `run_fitted` with the
    input as an argument) reproduces E107's saved runs validation2_fitted_<wiring> (0 differences);
  - the bound (bound.csv): the family's two ρ's (all eight, S and T; against all recordings and each half) with every
    pair of paa, B-myr and 2-but (42 pairs) put at the flies' values (E106's `at_data_values`, 'raw' and 'quantile'),
    and the same for 2h and 6al (29 pairs, the controls' room);
  - the references (references.csv): the family at the lab's input (E108's `unit_scores` on all eight, S, T and each
    wiring, raw and calibrated), E107's validation2 reading of the family (eight wirings and T, as fitted, ρ against
    Remy's somata over the 136 pairs and the 31 of paa and B-myr), and E108's carrying input on T when E108's runs are
    there (beside);
  - the yardsticks (yardsticks.csv: E106's, as E108 reads them, with the gate's thresholds beside) and probes.csv (the
    move planner's probes beside the values computed here).
  Writes part0.json ('all_exact'); raises if anything differs. No search starts before part 0 has passed.
- search (`--search NAME`; resumable from a checkpoint per generation; one lane of four worker processes, one per wiring
  of S, `--threads` OpenMP threads each, 4 by default): E83 part C's search (`search`) with the objective changed:
  - SEARCHES: A1, A2 (half A, two draws), B1, B2 (half B), C1, D1 (E86's second split, one draw per half) for paa,
    B-myr and 2-but in E64's order; ctlA, ctlB (halves A and B) for 2h and 6al (E64's controls);
  - per odor in turn, factors in [0.5, 2] per glomerulus on the start input's response, a (1 + 32) evolution strategy for
    eight generations, sigma from 0.3 (E83 part C's constants), candidates simulated as pseudo-odors on the four
    wirings of S at their fits' parameters, the lab input's breadth held;
  - the objective, maximized: the sum of the two ρ's of S's population (raw counts) against the searched half, minus
    E83's rate penalty (10 per unit of relative deviation beyond 15 % of each wiring's own response fraction of the odor,
    APL on and APL off, summed over the four wirings) and minus the pattern-part penalty (`pattern_penalty`: 10 per unit
    of relative excess of S's pattern part under the calibrated readout against the searched half beyond 0.0015 above
    the start's);
  - the start rule: the six searches of paa, B-myr and 2-but refuse to start (StartRefused) until E108's class is
    written (`e108_class`, `require_e108_class`); if it is "better", paa starts from E108's carrying input in all six and
    paa's factors multiply that input; on any other class every odor starts from the lab's input. The control searches
    do not depend on E108 and start from the lab's input.
  Writes search_history_<name>.csv, search_candidates_<name>.csv, input_factors_<name>.csv, inputs/searched_<name>.csv,
  search_runs_<name>.csv (each wiring's final point with its realism), provenance_search_<name>.json; the final runs
  under outputs/E109_input_search_rho/S/ and the checkpoints under outputs/E109_input_search_rho/S/state_<name>/.
- read (simulations; resumable): the held readings on T (`part_read`): every finished search's input at 0.6 and at 1 of
  its log size; for each of the six, its factors permuted over each odor's responding glomeruli (E83's
  `permuted_factors`, seeds 1 to 4, E86's) at 0.6 and inverted (-0.6, E86's); the breadth recomputed for the input, the
  somatic threshold bisected to the wiring's APL-off response count, then the tied APL scale to its APL-on count (E85's
  `held_reading`). Rows in runs_R.csv (absolute parameters, factors, realism), provenance_R.json, the inputs in inputs/;
  runs under outputs/E109_input_search_rho/R/.
- step7 (simulations, beside, no gate): paa's and B-myr's factors of each of the six searches at 0.6 applied to the
  validation2 ORN input's own paa (10^-4) and B-myr, the four wirings of T run as fitted (`fitted_on`); runs under
  outputs/E109_input_search_rho/V/, rows in runs_V.csv.
- score (no simulations): readings.csv (E108's `unit_scores` per input, unit and readout), runs.csv, searches.csv, the
  gate (gate.csv, gate.json: `gate`, `class_of`, against the lab's input and, under "better", against E108's input
  beside), step6_entries.csv / step6_null.csv / step6_flies.csv (`recurring_entries`, once the six are finished),
  step7_paa_bmyr.csv / step7_remy.csv (`step7_tables`, when step 7's runs are there), comparison_table.md.
- variants: variants.csv and model_pair_distances_long.csv for the model-output archive.

The gate (the entry's, fixed before any run; `gate`): means over the four held-out wirings (T's population) against the
half each search never saw, in both directions of E77's split; a direction's value is the mean over its two draws.
- Gain: ρ with APL on higher by GAIN_ON (0.010) or more, or with APL off by GAIN_OFF (0.015) or more, on at least
  GAIN_WIRINGS (3) of the four wirings (each wiring's paired change, the mean over the direction's draws), in both
  directions; the other ρ not lower by more than OTHER_SLACK (0.005) in either direction; the same sign on the second
  split (the mean over C1 against D and D1 against C of T's change).
- Controls: the permuted factors (mean over the four seeds and the direction's draws), the inverted factors and the
  control search for 2h and 6al (each read against its unseen half), each gaining less than CONTROL_RATIO (0.5) times
  the candidate's gain (the mean over the two directions), on each ρ that gains.
- The other measures, as in E108's gate (each direction): the pattern part under the calibrated readout against the
  unseen half not higher by more than EQUAL_PATTERN (0.003), the raw-count loss against the unseen half not higher by
  more than EQUAL_RAW (0.005), the mean APL effect (calibrated) not lower by more than EQUAL_MEAN (0.01), every run inside
  the four bands, response rate ρ with APL on and with APL off (mean over the runs) not lower by more than
  RATE_RHO_SLACK (0.05). The amplitude part, the shape part, the spread ratio and Pearson r on the per-pair APL effect
  are given beside (no thresholds on them).
- Pending Tom's answers of 2026-11-17 14:44 (applied, and the class also given without them): the block not lower by
  more than BLOCK_SLACK (0.010); the natural-source alignment with APL on and with APL off, and the structure vs
  wild-type somata, not lower by more than STRUCTURE_SLACK (0.05); the alignment's on minus off given beside.
- Classes (`class_of`): "passes" (the gain and the controls hold and every other measure); "a price" (the gain and the
  controls hold and a measure fails, named with its size); "carries below the gate" (no gain, but a mean gain over the
  directions of CARRY_MEAN (0.005) or more on a ρ whose controls hold); "does not carry". Outcomes (a) to (d) beside.

Readings chosen where the entry is open (fixed here before any run):
- The rate penalty's reference ("the family's own" response fractions) and the pattern penalty's ("the unchanged
  fits'") are the start's: the saved fits when every odor starts from the lab's input (then both are exactly the
  family's own fits); under E108 "better", the four wirings at E108's paa with the fits' parameters, so that the start
  itself is not penalized.
- The pattern part in the search's penalty is S's population under the calibrated readout with each wiring's K held at
  its fit's (E83's PopulationObjective) against the searched half: the only half the search may see. The penalty is 10
  (E83's rate penalty constant) per unit of relative excess, a near-hard bound at the entry's margin of 0.0015.
- The gate's raw-count loss is against the unseen half (the gate reads every measure against the unseen half); the
  loss against all recordings is given beside.
- "Means over the four held-out wirings" is T's population (the mean matrix over the four wirings); each wiring's
  paired change gives the count of wirings for "at least 3 of the 4".
- A gain whose controls do not hold is "does not carry" (the controls hold is a condition of "carries below the gate"
  in the entry, and the classes "passes" and "a price" read the gain as the input's own).
- The runs that must lie inside the bands are those of the six searched inputs at 0.6 on T (the gate's reading).
- Under E108 "better" the class is carried against the family at the lab's input (I0) and given beside against the
  family at E108's input (E108's held runs of its carrying rule on T).
- The control searches are read with the same gate (gate rows "ctl"), and their gains enter the candidate's controls.
- Step 6's flies: paa's and B-myr's spread over the six validation2 ORN flies (olfdb dataset 4, E108's compact copy, in
  E108's units: to-avg-max, then the slope through the origin of the lab's four validation odors), 2-but's over the nine
  megamat ORN flies (E85's `spread_table`); the rule is E86's (five of six agree in sign, mean log2 factor beyond 0.5),
  with its label-permuted counts beside (E86's `five_of_six_null`).
- Step 7 reads Remy's matrices with E107's readers under E107's lock, from the copy in data/remy_validation2_kc/
  (`remy_copy`: E107's `remy_dir` pointed at its `data_dir` while reading, so the lock checks the copies' sha256).

Departures from the entry: none so far.

The holds (CLAUDE.md, "Connectome data (Tom, 2026-11-17)"): no connectome table is read beyond the model's own hemibrain
wiring and E72's re-draws of it (the eight wirings); nothing of FlyWire, FAFB or BANC; no table of Yang's own Claude's
project; no raw-only recording. Nothing is published to claude.ai.

Provenance: the other checkouts must be clean and this repository clean in every file a process loads (E108's
`_commits_clean`, every module imported first, PRELOAD); every row records `H.code_commits()`; every run's provenance
has `H.model_provenance(orn_deltas=...)`.

Outputs: <shared>/results/E109_input_search_rho/; spike counts and responses under outputs/E109_input_search_rho/
(S: the searches' final points; R: the held readings; V: step 7), then data/model_outputs/E109_input_search_rho/.

Usage (from the repository root; one 16-thread lane on viki, under systemd-run --user):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E109_input_search_rho.py --part 0
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E109_input_search_rho.py --part search --search ctlA
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u experiments/E109_input_search_rho.py --part search --search A1 \\
        [--e108-class <class> --e108-class-source '<where it is written>']
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E109_input_search_rho.py --part read
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E109_input_search_rho.py --part step7
    ./venv/bin/python -u experiments/E109_input_search_rho.py --part score
"""
from __future__ import annotations

import argparse
import contextlib
import gc
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
from experiments import E83_searches as SR

EXPERIMENT = 'E109_input_search_rho'
SCRIPT = 'experiments/E109_input_search_rho.py'
E108_EXPERIMENT = 'E108_paa_extrapolation'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
UNITS = (('all8', ALL_WIRINGS), ('S', S_WIRINGS), ('T', T_WIRINGS))
SEARCH_ODORS = ('paa', 'B-myr', '2-but')         # E64's order
CONTROL_ODORS = ('2h', '6al')                    # E64's controls
HALVES = ('A', 'B', 'C', 'D')
TARGETS = ('all',) + HALVES                      # E108's: 'all' = all recordings
UNSEEN = {'A': 'B', 'B': 'A', 'C': 'D', 'D': 'C'}
SPLIT_OF = {'A': 0, 'B': 0, 'C': 1, 'D': 1}
# the searches: half searched, draw, generator seed (new, one per search), odors, group
SEARCHES = {
    'A1': dict(half='A', draw=1, seed=1091, odors=SEARCH_ODORS, group='main'),
    'A2': dict(half='A', draw=2, seed=1092, odors=SEARCH_ODORS, group='main'),
    'B1': dict(half='B', draw=1, seed=1093, odors=SEARCH_ODORS, group='main'),
    'B2': dict(half='B', draw=2, seed=1094, odors=SEARCH_ODORS, group='main'),
    'C1': dict(half='C', draw=1, seed=1095, odors=SEARCH_ODORS, group='main'),
    'D1': dict(half='D', draw=1, seed=1096, odors=SEARCH_ODORS, group='main'),
    'ctlA': dict(half='A', draw=1, seed=1097, odors=CONTROL_ODORS, group='control'),
    'ctlB': dict(half='B', draw=1, seed=1098, odors=CONTROL_ODORS, group='control'),
}
MAIN = ('A1', 'A2', 'B1', 'B2', 'C1', 'D1')
FIRST_SPLIT = ('A1', 'A2', 'B1', 'B2')
SECOND_SPLIT = ('C1', 'D1')
CONTROL_SEARCHES = ('ctlA', 'ctlB')
# E83 part C's search constants
POP, GENS, SIGMA0 = SR.I_POP, SR.I_GENS, SR.I_SIGMA0
BOUNDS = SR.I_BOUNDS
RATE_TOL, PENALTY = SR.I_RATE_TOL, SR.I_PENALTY
WORKER_THREADS = 4                               # four workers, one 16-thread lane
PATTERN_MARGIN = 0.0015
# the readings
GATE_SIZE, FULL_SIZE = E83.GATE_SIZE, 1.0
INVERT_SIZE = -0.6                               # E86's
PERM_SEEDS = (1, 2, 3, 4)                        # E86's B1 and B3 seeds of E83's `permuted_factors`
# the gate (fixed in the entry)
GAIN = {'on': 0.010, 'off': 0.015}
GAIN_WIRINGS = 3
OTHER_SLACK = 0.005
CONTROL_RATIO = 0.5
CARRY_MEAN = 0.005
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
RATE_RHO_SLACK = 0.05
BLOCK_SLACK, STRUCTURE_SLACK = 0.010, 0.05
CLASSES = ('passes', 'a price', 'carries below the gate', 'does not carry')
OUTCOME = {'passes': 'a', 'a price': 'b', 'carries below the gate': 'c', 'does not carry': 'd'}
# the move planner's probes (the entry, 2026-11-17 23:40 to 23:50), for probes.csv
PROBES = {
    ('bound', 'all8', 'three', 'raw', 'all'): (0.807, 0.880), ('bound', 'all8', 'three', 'quantile', 'all'): (0.840, 0.899),
    ('bound', 'all8', 'three', 'as_is', 'all'): (0.690, 0.756), ('bound', 'T', 'three', 'raw', 'all'): (0.805, 0.886),
    ('bound', 'T', 'three', 'quantile', 'all'): (0.838, 0.907), ('bound', 'T', 'three', 'as_is', 'all'): (0.685, 0.762),
    ('bound', 'all8', 'controls', 'raw', 'all'): (0.753, 0.794),
    ('e107', 'all8', '136', 'fitted', 'remy'): (0.676, np.nan), ('e107', 'all8', '31', 'fitted', 'remy'): (0.542, np.nan),
}
# step 7's yardsticks (E107; cosine distance of paa from B-myr)
PAA_BMYR = ('paa', 'B-myr')
PRELOAD = ('experiments.E10_joint_fit', 'experiments.E64_input_oracle', 'experiments.E67_apl_effect_fits',
           'experiments.E72_wiring_population', 'experiments.E77_input_screen',
           'experiments.E77_input_side_common_error', 'experiments.E85_orn_fly_population',
           'experiments.E86_reduced_input_candidate', 'experiments.E89_apl_breadth_weighting',
           'experiments.E106_rho_gap_parts', 'experiments.E107_validation2_panel', 'experiments.E108_paa_extrapolation',
           'kcq.apl_effect', 'pyolfsysm.mb_model', 'pyolfsysm.tuning')

_CACHE: Dict[object, object] = {}


class StartRefused(RuntimeError):
    """A search of the three odors asked to start before E108's class is written (the message names what it waits
    for)."""


# ---------------------------------------------------------------------------
# folders, files, provenance

def out_dir() -> Path:
    """results/E109_input_search_rho/ in the shared directory (E109_OUT_DIR in the environment instead, for smoke
    tests)."""
    d = Path(os.environ['E109_OUT_DIR']) if os.environ.get('E109_OUT_DIR') else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    """outputs/E109_input_search_rho/<sub> (E109_LOCAL_DIR in the environment instead, for smoke tests)."""
    root = Path(os.environ['E109_LOCAL_DIR']) if os.environ.get('E109_LOCAL_DIR') else H.outputs_dir(EXPERIMENT)
    d = root / sub if sub else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_json(path: Path, obj) -> None:
    SR._write_json(path, json.loads(json.dumps(obj, default=lambda x: x.item() if hasattr(x, 'item') else str(x))))


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def _flagged(x) -> bool:
    return isinstance(x, str) and x.strip() not in ('', 'nan')


def commits_clean() -> dict:
    """`H.code_commits()` at the start of a simulation part, after importing every module the parts use (PRELOAD):
    E108's `_commits_clean` (the other checkouts clean; this repository clean in every file this process loads)."""
    from experiments import E108_paa_extrapolation as E108
    for m in PRELOAD:
        importlib.import_module(m)
    return E108._commits_clean()


def _provenance_add(path: Path, part: str, tag: str, prov: dict) -> None:
    d = _read_json(path) or dict(part=part, script=SCRIPT, points={})
    d['points'][tag] = prov
    _write_json(path, d)


# ---------------------------------------------------------------------------
# E108's class and the start rule

def e108_class(results: Optional[Path] = None) -> dict:
    """E108's class as its results give it (results/E108_paa_extrapolation/classes.csv, the carrying rule's row).

    Returns:
        dict with 'written' (the class is written: the carrying rule's runs complete, its class there and not marked
        provisional), 'cls' (the class, with Tom's pending items applied, as E108's `class_of`), 'class_without_pending',
        'carrying_rule', 'complete', 'provisional' (the name of the column that marks it provisional, if it does),
        'source' (the file), 'waiting_for' (what is missing, None when written).
    """
    p = Path(results or H.results_dir(E108_EXPERIMENT)) / 'classes.csv'
    out = dict(written=False, cls=None, class_without_pending=None, carrying_rule=None, complete=False, provisional=None,
               source=str(p), waiting_for=None)
    if not p.exists():
        out['waiting_for'] = (f"E108's classes.csv ({p}): E108 (paa's input at the flies' concentration, hal's) has not "
                              "scored its runs")
        return out
    ct = pd.read_csv(p)
    true = lambda v: str(v).strip().lower() in ('true', '1', '1.0')
    car = ct[[true(v) for v in ct.get('carrying', pd.Series(False, index=ct.index))]]
    if len(car) != 1:
        out['waiting_for'] = f"E108's carrying rule in {p} ({len(car)} rows marked carrying)"
        return out
    r = car.iloc[0]
    out['carrying_rule'] = str(r['input'])
    out['complete'] = true(r.get('complete', False))
    cls = r.get('class') if 'class' in ct.columns else None
    out['cls'] = None if cls is None or pd.isna(cls) or not str(cls).strip() else str(cls)
    cwp = r.get('class_without_pending') if 'class_without_pending' in ct.columns else None
    out['class_without_pending'] = None if cwp is None or pd.isna(cwp) else str(cwp)
    prov = [c for c in ct.columns if str(c).startswith('provisional') and true(r[c])]
    out['provisional'] = prov[0] if prov else None
    if not out['complete'] or out['cls'] is None:
        out['waiting_for'] = (f"E108's runs of its carrying rule {out['carrying_rule']} on the family's eight wirings and "
                              f"its score (classes.csv has no class yet: complete={out['complete']})")
    elif out['provisional']:
        out['waiting_for'] = (f"E108's class to be written in its note or its Status (classes.csv marks the class "
                              f"'{out['cls']}' as '{out['provisional']}'); once it is, start with --e108-class <the class "
                              "as written> --e108-class-source '<where it is written>'")
    else:
        out['written'] = True
    return out


def require_e108_class(name: str, confirmed: Optional[Mapping[str, str]] = None,
                       results: Optional[Path] = None) -> dict:
    """The start rule for search `name`: the searches of SEARCH_ODORS start only when E108's class is written; the
    control searches do not depend on it.

    Args:
        name: a key of SEARCHES.
        confirmed: {'cls': the class as written in E108's note or Status, 'source': where} when classes.csv still marks
            the class provisional (E108's code marks it so until its write-up); accepted only if E108's runs are complete
            and `cls` is classes.csv's class (with or without the pending items).
        results: E108's results folder (tests).

    Returns:
        `e108_class`'s dict with 'applies' (whether the rule applies to this search) and, for the six, 'start'
        ('e108' when the class is "better", else 'lab').

    Raises:
        StartRefused: naming what the search waits for.
    """
    if SEARCHES[name]['group'] != 'main':
        return dict(applies=False, start='lab')
    st = dict(e108_class(results), applies=True)
    if confirmed:
        c = str(confirmed.get('cls', '')).strip()
        if not (st['complete'] and st['cls'] is not None):
            raise StartRefused(f"E109 {name}: --e108-class given, but {st['waiting_for']}")
        if c not in {st['cls'], st['class_without_pending']}:
            raise StartRefused(f"E109 {name}: --e108-class '{c}' is not E108's class in {st['source']} ('{st['cls']}', "
                               f"without the pending items '{st['class_without_pending']}')")
        st.update(written=True, cls=c, written_in=str(confirmed.get('source') or ''), waiting_for=None)
    if not st['written']:
        raise StartRefused(f"E109 {name} (a search of paa, B-myr and 2-but) waits for E108's class: {st['waiting_for']}")
    st['start'] = 'e108' if st['cls'] == 'better' else 'lab'
    return st


def start_input(e108: Optional[Mapping] = None) -> Tuple[pd.DataFrame, dict]:
    """The input a search starts from: the lab's (`H.apl16_orn_deltas()`), or, when `e108` says 'start' 'e108' (E108's
    class "better"), the lab's with paa's column E108's carrying rule's (E108's `input_of`, checked against the file its
    part inputs recorded, `input_record`). Returns (input, info: 'start', 'rule', 'record', 'frame_sha256')."""
    from experiments.E77_input_screen import frame_sha256
    lab = H.apl16_orn_deltas()
    if not e108 or e108.get('start') != 'e108':
        return lab, dict(start='lab', rule=None, record=None, frame_sha256=frame_sha256(lab))
    from experiments import E108_paa_extrapolation as E108
    rule = str(e108['carrying_rule'])
    orn = E108.input_of(rule)
    rec = E108.input_record(rule)
    if not (orn.index.equals(lab.index) and orn.columns.equals(lab.columns)):
        raise RuntimeError(f"E108's input {rule} is not in the lab input's layout")
    differ = [H.odor_abbrev(c) for c in lab.columns if not np.array_equal(orn[c].to_numpy(float), lab[c].to_numpy(float))]
    if differ != ['paa']:
        raise RuntimeError(f"E108's input {rule} differs from the lab's in {differ}, not in paa's column alone")
    return orn, dict(start='e108', rule=rule, record=rec, frame_sha256=frame_sha256(orn))


def base_of(start: Mapping) -> pd.DataFrame:
    """The start input named by a search's start info (`start_input`'s)."""
    return start_input(dict(start=start['start'], carrying_rule=start.get('rule')))[0]


# ---------------------------------------------------------------------------
# the targets and the objective

def dist_targets() -> Dict[str, Dict[str, pd.Series]]:
    """The KC data's mean distances per condition for 'all' and the halves A to D (E108's `dist_targets`)."""
    from experiments import E108_paa_extrapolation as E108
    return E108.dist_targets()


def effect_targets() -> Dict[str, pd.Series]:
    """The per-pair APL effect for 'all' and the halves A to D (E108's `effect_targets`)."""
    from experiments import E108_paa_extrapolation as E108
    return E108.effect_targets()


def _counts(x) -> pd.DataFrame:
    return x[1] if isinstance(x, tuple) else x


def population_distances(counts: Mapping[str, object]) -> pd.Series:
    """The mean over the given wirings of their raw pair distances (cosine on spike counts, E83's `pair_distances`
    'raw'), canonical order (E83's)."""
    return E83._cs(pd.concat([E83._cs(H.readout_pair_distances(_counts(c), 'raw')) for c in counts.values()],
                             axis=1).mean(axis=1))


def two_rhos_population(on: Mapping[str, object], off: Mapping[str, object], data_on: pd.Series,
                        data_off: pd.Series) -> Dict[str, float]:
    """The two ρ's of the mean matrices over the given wirings (`on` / `off`: {wiring: spike counts or (responses,
    spike counts)}, raw counts) against one half's KC matrices (E106's `two_rhos`): {'on', 'off', 'sum'}."""
    from experiments import E106_rho_gap_parts as E106
    r = E106.two_rhos(population_distances(on), population_distances(off), data_on, data_off)
    return dict(on=float(r['on']), off=float(r['off']), sum=float(r['on'] + r['off']))


def pattern_penalty(rms: float, rms0: float, margin: float = PATTERN_MARGIN, penalty: float = PENALTY) -> float:
    """The penalty on a pattern part `rms` more than `margin` above the start's `rms0`: `penalty` per unit of relative
    excess, penalty x max(0, rms - rms0 - margin) / rms0."""
    return float(penalty * max(0.0, float(rms) - float(rms0) - float(margin)) / float(rms0))


class RhoObjective:
    """S's population read against one half: the two ρ's on raw counts (`two_rhos_population`) and the pattern part under
    the calibrated readout with each wiring's K held at its fit's (E83's PopulationObjective, the half's per-pair APL
    effect). `obj(on, off)` with {wiring: spike counts} returns dict(rho_on, rho_off, rho_sum, pattern)."""

    def __init__(self, half: str, wirings: Sequence[str] = S_WIRINGS):
        dt = dist_targets()[half]
        self.half = half
        self.data_on, self.data_off = dt['on'], dt['off']
        self.wirings = tuple(wirings)
        self.pattern = E83.PopulationObjective(effect_targets()[half], self.wirings)

    def __call__(self, on: Mapping[str, object], off: Mapping[str, object]) -> dict:
        r = two_rhos_population(on, off, self.data_on, self.data_off)
        p = self.pattern(on={w: _counts(v) for w, v in on.items()}, off={w: _counts(v) for w, v in off.items()})
        return dict(rho_on=r['on'], rho_off=r['off'], rho_sum=r['sum'], pattern=float(p['rms']))


def score_of(v: Mapping[str, float], fracs: Mapping[str, Tuple[float, float]], base_fracs: Mapping[str, Tuple[float, float]],
             rms0: float) -> dict:
    """A point's search score: the sum of the two ρ's minus E83's rate penalty on the odor's fractions (`fracs` against
    `base_fracs`, per wiring (APL on, APL off)) minus the pattern-part penalty (against the start's `rms0`)."""
    pr = SR.rate_penalty(fracs, base_fracs)
    pp = pattern_penalty(v['pattern'], rms0)
    return dict(v, rate_penalty=pr, pattern_penalty=pp, score=float(v['rho_sum'] - pr - pp))


def _odor_fracs(inc: Mapping[str, Mapping[str, tuple]], odor: str) -> Dict[str, Tuple[float, float]]:
    """Each wiring's response fraction of `odor`, APL on and APL off, in the incumbents `inc`."""
    oc = SR._odor_col(next(iter(inc.values()))['on'][1], odor)
    return {w: (float(inc[w]['on'][0][oc].astype(float).gt(0.5).mean()), float(inc[w]['off'][0][oc].astype(float).gt(0.5).mean()))
            for w in inc}


def incumbent_score(inc: Mapping[str, Mapping[str, tuple]], odor: str, obj: RhoObjective,
                    base_fracs: Mapping[str, Tuple[float, float]], rms0: float) -> dict:
    """The incumbents' score with `odor` the odor searched."""
    v = obj(on={w: inc[w]['on'][1] for w in inc}, off={w: inc[w]['off'][1] for w in inc})
    return score_of(v, _odor_fracs(inc, odor), base_fracs, rms0)


def candidate_scores(inc: Mapping[str, Mapping[str, tuple]], odor: str, res: Mapping[str, Tuple[np.ndarray, ...]],
                     obj: RhoObjective, base_fracs: Mapping[str, Tuple[float, float]], rms0: float) -> List[dict]:
    """Each candidate's score: the odor's column of every wiring's incumbent replaced by the candidate's (column i of
    `res[w]`, E83's InputPool layout: responses on, counts on, responses off, counts off), scored by `score_of`."""
    oc = SR._odor_col(next(iter(inc.values()))['on'][1], odor)
    n_cand = next(iter(res.values()))[1].shape[1]
    rows = []
    for i in range(n_cand):
        on = {w: SR._with_column(inc[w]['on'][1], oc, res[w][1][:, i]) for w in inc}
        off = {w: SR._with_column(inc[w]['off'][1], oc, res[w][3][:, i]) for w in inc}
        fr = {w: (float((res[w][0][:, i] > 0.5).mean()), float((res[w][2][:, i] > 0.5).mean())) for w in inc}
        rows.append(dict(score_of(obj(on, off), fr, base_fracs, rms0), fracs=fr))
    return rows


# ---------------------------------------------------------------------------
# the searches

def _sdir(name: str) -> Path:
    return local_dir(f'S/state_{name}')


def factors_path(name: str) -> Path:
    return out_dir() / f'input_factors_{name}.csv'


def write_factors(name: str, base: pd.DataFrame, done: Mapping[str, Sequence[float]], start: Mapping) -> pd.DataFrame:
    """input_factors_<name>.csv (search, half, odor, glomerulus, base_delta (the start input's), lab_delta, factor,
    log2_factor, searched_delta; the odors searched so far, in the search's order) and inputs/searched_<name>.csv (the
    searched input, glomeruli x abbreviated odors, E38's layout)."""
    cfg = SEARCHES[name]
    lab = H.apl16_orn_deltas()
    ab = [H.odor_abbrev(c) for c in base.columns]
    now = base.astype(float).copy()
    rows = []
    for odor in cfg['odors']:
        if odor not in done:
            continue
        f = np.asarray(done[odor], dtype=float)
        j = ab.index(odor)
        now.iloc[:, j] = base.iloc[:, j].to_numpy(dtype=float) * f
        for gi, g in enumerate(base.index):
            rows.append(dict(search=name, half=cfg['half'], start=start['start'], odor=odor, glomerulus=g,
                             base_delta=float(base.iloc[gi, j]), lab_delta=float(lab.iloc[gi, j]), factor=float(f[gi]),
                             log2_factor=float(np.log2(f[gi])), searched_delta=float(base.iloc[gi, j] * f[gi])))
    fac = pd.DataFrame(rows)
    fac.to_csv(factors_path(name), index=False)
    E83._write_input(out_dir(), f'searched_{name}', now)
    return fac


def search_finished(name: str) -> bool:
    """Whether search `name` has searched all its odors (its factor table lists them all)."""
    p = factors_path(name)
    return p.exists() and set(pd.read_csv(p)['odor']) == set(SEARCHES[name]['odors'])


def search_factors(name: str) -> pd.DataFrame:
    """A finished search's factors, glomeruli x odors (abbreviated, the lab input's order; 1 where not searched)."""
    if not search_finished(name):
        raise RuntimeError(f'search {name} is not finished ({factors_path(name)})')
    return E83.read_factors(factors_path(name))


def search_start(name: str) -> dict:
    """The start info a search recorded (provenance_search_<name>.json)."""
    p = out_dir() / f'provenance_search_{name}.json'
    st = _read_json(p).get('start')
    if not st:
        raise RuntimeError(f'{p}: no start info')
    return st


def part0_passed() -> bool:
    return bool(_read_json(out_dir() / 'part0.json').get('all_exact'))


def search(name: str, resume: bool = True, threads: int = WORKER_THREADS, confirmed: Optional[Mapping[str, str]] = None,
           gens: int = GENS, pop: int = POP) -> pd.DataFrame:
    """Search `name` (SEARCHES; see the module docstring): the start rule, the start input, the identity of the pseudo-odor
    path, then per odor the (1 + `pop`) evolution strategy for `gens` generations on S's four wirings, maximizing the
    score (`score_of`) against the searched half; checkpointed every generation, resumed by default. Returns the
    history."""
    from experiments.E77_input_screen import frame_sha256
    cfg = SEARCHES[name]
    half, unseen, odors = cfg['half'], UNSEEN[cfg['half']], tuple(cfg['odors'])
    e108 = require_e108_class(name, confirmed)
    if not part0_passed():
        raise RuntimeError('run --part 0 first (part0.json: the identities must hold before any search)')
    cc = commits_clean()
    base, start = start_input(e108 if e108.get('applies') else None)
    out = out_dir()
    sdir = _sdir(name)
    ck = sdir / 'checkpoint.json'
    hist_path, cand_path = out / f'search_history_{name}.csv', out / f'search_candidates_{name}.csv'
    lab = H.apl16_orn_deltas()
    glomeruli = list(lab.index)
    ab = [H.odor_abbrev(c) for c in lab.columns]
    differs = [o for o in ab if not np.array_equal(base[lab.columns[ab.index(o)]].to_numpy(float),
                                                     lab[lab.columns[ab.index(o)]].to_numpy(float))]
    saved = {w: E83.saved_fit(w)['saved'] for w in S_WIRINGS}
    obj, other = RhoObjective(half), RhoObjective(unseen)
    lo, hi = np.log(BOUNDS[0]), np.log(BOUNDS[1])
    t_start = time.time()
    pool = SR.InputPool(S_WIRINGS, threads, orn=None if start['start'] == 'lab' else base)
    try:
        # the pseudo-odor path at factor 1 of an odor the start does not change reproduces the saved runs
        check = next(o for o in odors + SEARCH_ODORS + CONTROL_ODORS if o not in differs)
        res = pool.run(check, np.ones((1, len(glomeruli))))
        oc = SR._odor_col(saved['real']['on'][1], check)
        for w in S_WIRINGS:
            for c, kr, kc in (('on', 0, 1), ('off', 2, 3)):
                if not (np.array_equal(res[w][kc][:, 0].astype(float), saved[w][c][1][oc].to_numpy(dtype=float))
                        and np.array_equal(res[w][kr][:, 0].astype(float), saved[w][c][0][oc].to_numpy(dtype=float))):
                    raise RuntimeError(f'{w}: the pseudo-odor run of {check} at factor 1 differs from the saved {c} run')
        print(f'{name}: the pseudo-odor runs of {check} at factor 1 reproduce the four saved fits', flush=True)
        if resume and ck.exists():
            st = json.loads(ck.read_text())
            if st.get('start', {}).get('frame_sha256') != start['frame_sha256']:
                raise RuntimeError(f'{name}: the checkpoint started from another input; run with --no-resume')
            done = {o: list(v) for o, v in st['done'].items()}
            rng = SR._rng_from(st['rng'])
            oi, g = int(st['odor_index']), int(st['generation'])
            frames = st['frames']
            inc = {w: {c: (SR._load_frame(sdir, f'{w}__{frames}', f'{c}_responses'),
                           SR._load_frame(sdir, f'{w}__{frames}', f'{c}_spike_counts')) for c in ('on', 'off')}
                   for w in S_WIRINGS}
            base_fracs = {o: {w: tuple(v) for w, v in d.items()} for o, d in st['base_fracs'].items()}
            rms0 = float(st['rms0'])
            parent, parent_score, sigma = None, None, SIGMA0
            if g + 1 < gens:
                start_odor, start_gen = oi, g + 1
                parent, parent_score, sigma = np.array(st['parent']), float(st['parent_score']), float(st['sigma'])
                now = incumbent_score(inc, odors[oi], obj, base_fracs[odors[oi]], rms0)['score']
                if abs(now - parent_score) > 1e-12:
                    raise RuntimeError(f'{name}: the checkpoint parent scores {now} now, {parent_score} when saved')
            else:
                start_odor, start_gen = oi + 1, 0
            pos = oi * gens + g
            order = {o: k for k, o in enumerate(odors)}
            SR._trim_csv(hist_path, lambda d: (d.odor == 'start') | (d.odor.map(order) * gens + d.generation <= pos))
            SR._trim_csv(cand_path, lambda d: d.odor.map(order) * gens + d.generation <= pos)
            print(f'{name}: resuming at odor {start_odor} generation {start_gen}', flush=True)
        else:
            for p in (hist_path, cand_path, factors_path(name)):
                if p.exists():
                    p.unlink()
            SR._reset_state(sdir)
            inc = {w: {c: tuple(x.copy() for x in saved[w][c]) for c in ('on', 'off')} for w in S_WIRINGS}
            # an odor whose start column is not the lab's starts from its own column at factor 1
            for odor in [o for o in odors if o in differs]:
                r1 = pool.run(odor, np.ones((1, len(glomeruli))))
                ocd = SR._odor_col(inc['real']['on'][1], odor)
                for w in S_WIRINGS:
                    inc[w] = {'on': (SR._with_column(inc[w]['on'][0], ocd, r1[w][0][:, 0]),
                                     SR._with_column(inc[w]['on'][1], ocd, r1[w][1][:, 0])),
                              'off': (SR._with_column(inc[w]['off'][0], ocd, r1[w][2][:, 0]),
                                      SR._with_column(inc[w]['off'][1], ocd, r1[w][3][:, 0]))}
            base_fracs = {o: _odor_fracs(inc, o) for o in odors}
            rms0 = float(obj(on={w: inc[w]['on'][1] for w in S_WIRINGS}, off={w: inc[w]['off'][1] for w in S_WIRINGS})['pattern'])
            done, rng = {}, np.random.default_rng(cfg['seed'])
            start_odor, start_gen, parent, parent_score, sigma = 0, 0, None, None, SIGMA0
            frames = 'start'
            for w in S_WIRINGS:
                SR._save_frames(sdir, f'{w}__{frames}', {'on_responses': inc[w]['on'][0], 'on_spike_counts': inc[w]['on'][1],
                                                         'off_responses': inc[w]['off'][0], 'off_spike_counts': inc[w]['off'][1]})
            v0 = obj(on={w: inc[w]['on'][1] for w in S_WIRINGS}, off={w: inc[w]['off'][1] for w in S_WIRINGS})
            u0 = other(on={w: inc[w]['on'][1] for w in S_WIRINGS}, off={w: inc[w]['off'][1] for w in S_WIRINGS})
            SR._append_csv(hist_path, [dict(search=name, half=half, odor='start', generation=-1, rho_on_half=v0['rho_on'],
                                            rho_off_half=v0['rho_off'], pattern_half=v0['pattern'], rho_on_unseen=u0['rho_on'],
                                            rho_off_unseen=u0['rho_off'], pattern_unseen=u0['pattern'],
                                            seconds=time.time() - t_start)])
            print(f"{name}: start ({start['start']}) S's two ρ's against half {half} {v0['rho_on']:.4f} / {v0['rho_off']:.4f}, "
                  f"pattern part {v0['pattern']:.4f}; against the unseen half {unseen} {u0['rho_on']:.4f} / "
                  f"{u0['rho_off']:.4f}", flush=True)
        old = _read_json(out / f'provenance_search_{name}.json')
        runs_so_far = list(old.get('invocations', [])) + [dict(generated=_now(), resumed=bool(resume and old), **cc)]
        prov = dict(generated=_now(), script=SCRIPT, search=name, **cfg, unseen=unseen, **cc, invocations=runs_so_far,
                    model=H.model_provenance(orn_deltas=base, _wPNKC=E83.wiring_parts('real')['W']),
                    start=dict(start, base_fracs=base_fracs, rms0=rms0, differs_from_lab=differs), e108=e108,
                    wirings={w: dict(folder=E83.saved_fit(w)['folder'], tag=E83.saved_fit(w)['tag'], params=E83.saved_fit(w)['p'],
                                     K=E83.saved_fit(w)['K']) for w in S_WIRINGS},
                    halves={h: {g: [list(map(str, r)) for r in v] for g, v in hh.items()}
                            for h, hh in _halves().items() if h in (half, unseen)},
                    pop=pop, generations=gens, bounds=BOUNDS, sigma0=SIGMA0, sigma_min=SR.SIGMA_MIN, sigma_max=SR.SIGMA_MAX,
                    success_up=SR.SUCCESS_UP, failure_down=SR.FAILURE_DOWN, rate_tol=RATE_TOL, penalty=PENALTY,
                    pattern_margin=PATTERN_MARGIN, workers=len(S_WIRINGS), threads_per_worker=threads,
                    objective='sum of the two ρ on odor-odor distances of S\'s population (raw counts) against the searched '
                              'half, minus the rate penalty and the pattern-part penalty',
                    breadth="each wiring's lab-input breadth (E83's wiring_parts) held for every candidate",
                    approx_conc_note=H.APPROX_CONC_NOTE)
        _write_json(out / f'provenance_search_{name}.json', prov)
        for oi in range(start_odor, len(odors)):
            odor = odors[oi]
            oc = SR._odor_col(inc['real']['on'][1], odor)
            if oi == start_odor and parent is not None:
                g_first = start_gen
            else:
                parent = np.ones(len(glomeruli))
                parent_score = incumbent_score(inc, odor, obj, base_fracs[odor], rms0)['score']
                sigma, g_first = SIGMA0, 0
            for g in range(g_first, gens):
                t0 = time.time()
                offspring = np.exp(SR.es_offspring(rng, np.log(parent), sigma, pop, lo, hi))
                res = pool.run(odor, offspring)
                rows = candidate_scores(inc, odor, res, obj, base_fracs[odor], rms0)
                scores = np.array([r['score'] for r in rows])
                b = int(np.argmax(scores))
                success = bool(scores[b] > parent_score)
                if success:
                    parent, parent_score = offspring[b].copy(), float(scores[b])
                    for w in S_WIRINGS:
                        r_on, s_on = inc[w]['on']
                        r_off, s_off = inc[w]['off']
                        inc[w] = {'on': (SR._with_column(r_on, oc, res[w][0][:, b]), SR._with_column(s_on, oc, res[w][1][:, b])),
                                  'off': (SR._with_column(r_off, oc, res[w][2][:, b]), SR._with_column(s_off, oc, res[w][3][:, b]))}
                sigma = SR.es_sigma(sigma, success)
                on = {w: inc[w]['on'][1] for w in S_WIRINGS}
                off = {w: inc[w]['off'][1] for w in S_WIRINGS}
                ph, pu = obj(on, off), other(on, off)
                fr = _odor_fracs(inc, odor)
                hrow = dict(search=name, half=half, odor=odor, generation=g, best_score=float(scores[b]),
                            parent_score=parent_score, success=success, sigma=sigma, best_rho_on=rows[b]['rho_on'],
                            best_rho_off=rows[b]['rho_off'], best_pattern=rows[b]['pattern'],
                            best_rate_penalty=rows[b]['rate_penalty'], best_pattern_penalty=rows[b]['pattern_penalty'],
                            rho_on_half=ph['rho_on'], rho_off_half=ph['rho_off'], pattern_half=ph['pattern'],
                            rho_on_unseen=pu['rho_on'], rho_off_unseen=pu['rho_off'], pattern_unseen=pu['pattern'],
                            **{f'frac_on_{w}': fr[w][0] for w in S_WIRINGS}, **{f'frac_off_{w}': fr[w][1] for w in S_WIRINGS},
                            **{f'base_frac_on_{w}': base_fracs[odor][w][0] for w in S_WIRINGS},
                            **{f'base_frac_off_{w}': base_fracs[odor][w][1] for w in S_WIRINGS},
                            log2_factor_rms=float(np.sqrt(np.mean(np.log2(parent) ** 2))),
                            n_at_bound=int((np.abs(np.abs(np.log(parent)) - hi) < 1e-12).sum()), seconds=time.time() - t0)
                SR._append_csv(cand_path, [dict(search=name, half=half, odor=odor, generation=g, candidate=i, score=r['score'],
                                                rho_on=r['rho_on'], rho_off=r['rho_off'], pattern=r['pattern'],
                                                rate_penalty=r['rate_penalty'], pattern_penalty=r['pattern_penalty'])
                                           for i, r in enumerate(rows)])
                if success:
                    frames = f'o{oi}g{g}'
                    for w in S_WIRINGS:
                        SR._save_frames(sdir, f'{w}__{frames}', {'on_responses': inc[w]['on'][0], 'on_spike_counts': inc[w]['on'][1],
                                                                 'off_responses': inc[w]['off'][0], 'off_spike_counts': inc[w]['off'][1]})
                if g == gens - 1:
                    done[odor] = parent.tolist()
                SR._write_json(ck, dict(search=name, odor_index=oi, generation=g, frames=frames, parent=parent.tolist(),
                                        parent_score=parent_score, sigma=sigma, rng=SR._rng_state(rng), done=done,
                                        base_fracs={o: {w: list(v) for w, v in d.items()} for o, d in base_fracs.items()},
                                        rms0=rms0, start=dict(start=start['start'], rule=start['rule'],
                                                              frame_sha256=start['frame_sha256'])))
                SR._drop_other_frames(sdir, frames)
                SR._append_csv(hist_path, [hrow])
                print(f"{name} {odor} gen {g}: best {scores[b]:.4f} (parent {parent_score:.4f}); S against half {half} "
                      f"{ph['rho_on']:.4f} / {ph['rho_off']:.4f} (pattern {ph['pattern']:.4f}), against {unseen} "
                      f"{pu['rho_on']:.4f} / {pu['rho_off']:.4f}; sigma {sigma:.2f}; {time.time() - t0:.0f}s", flush=True)
            parent = None
            write_factors(name, base, done, start)
    finally:
        pool.close()
    # the final point: each wiring's runs, absolute parameters and realism
    rows = []
    for w in S_WIRINGS:
        o = {c: inc[w][c] for c in ('on', 'off')}
        tag = f'{name}_{w}'
        save_run('S', tag, o)
        P = E83.saved_fit(w)
        rows.append(dict(search=name, half=half, wiring=w, tag=tag, fit_folder=P['folder'], fit_tag=P['tag'], apl_factor=1.0,
                         thr_factor=1.0, n_on=E83._n(o['on'][0]), target_on=P['target_on'], n_off=E83._n(o['off'][0]),
                         target_off=P['target_off'], K=P['K'], start=start['start'], **P['p'], **E83.realism(o), **cc))
    pd.DataFrame(rows).to_csv(out / f'search_runs_{name}.csv', index=False)
    print(pd.DataFrame(rows)[['wiring', 'n_on', 'target_on', 'n_off', 'target_off', 'rate_rho_on', 'realism_flag']].to_string(index=False),
          flush=True)
    return pd.read_csv(hist_path)


def _halves() -> Dict[str, Dict[str, list]]:
    from experiments import E86_reduced_input_candidate as E86
    return E86.halves()


# ---------------------------------------------------------------------------
# runs

def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment ({'on', 'off'}): this host's outputs, else the model-output archive; or None."""
    for d in (local_dir(part), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def cache_tag(orn: pd.DataFrame, name: str) -> str:
    return E83.input_tag_for(orn, f'E109_{name}')


def held(label: str, orn: pd.DataFrame, name: str) -> dict:
    """A held reading (E85's `held_reading` with this experiment's breadth cache): hold_rates' dict."""
    from experiments import E85_orn_fly_population as E85
    return E85.held_reading(label, orn, cache_tag(orn, name), cache_dir=local_dir())


def fitted_on(label: str, orn: pd.DataFrame) -> Dict[str, object]:
    """One wiring at its fit on another panel's input `orn`, as E107's `run_fitted` runs it (E83's `build_simulator`, the
    lab input's breadth held, nothing re-tuned), APL on and acutely off, columns abbreviated: {'on', 'off', 'params',
    'params_off'}."""
    from experiments import E107_validation2_panel as E107
    P = E83.wiring_parts(label)
    p = dict(P['p'])
    p_off = dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    sim, info = E83.build_simulator(label, orn=orn, breadth=P['breadth'], params=p)
    on, off = E107._abbrev(sim(p)), E107._abbrev(sim(p_off))
    del sim, info
    return dict(on=on, off=off, params=p, params_off=p_off)


# ---------------------------------------------------------------------------
# part 0

def _saved_distances(readout: str) -> Dict[str, Dict[str, pd.Series]]:
    """{wiring: {'on', 'off'}} pair distances of the saved fits under `readout` (E83's `pair_distances`)."""
    key = ('saved_dist', readout)
    if key not in _CACHE:
        _CACHE[key] = {w: E83.pair_distances(E83.saved_fit(w)['saved'], readout) for w in ALL_WIRINGS}
    return _CACHE[key]


def family_population(unit: str, readout: str = 'raw') -> Dict[str, pd.Series]:
    """The saved fits' population over a unit's wirings (UNITS) or one wiring."""
    ws = dict(UNITS).get(unit, (unit,))
    d = _saved_distances(readout)
    return E83.population([d[w] for w in ws])


def bound_table() -> pd.DataFrame:
    """The bound (part 0): the family's two ρ's on all eight, S and T against all recordings and each half, as it is
    and with every pair of the odors (the three searched: 42 pairs; the controls 2h, 6al: 29) put at the flies' values
    (E106's `at_data_values`, 'raw' and 'quantile')."""
    from experiments import E106_rho_gap_parts as E106
    dt = dist_targets()
    rows = []
    for unit, _ in UNITS:
        fam = family_population(unit, 'raw')
        for which, odors in (('three', SEARCH_ODORS), ('controls', CONTROL_ODORS)):
            n_pairs = int(sum(1 for p in fam['on'].index if p[0] in odors or p[1] in odors))
            for t in TARGETS:
                for how in ('as_is', 'raw', 'quantile'):
                    row = dict(unit=unit, odors=which, n_pairs=n_pairs, target=t, how=how)
                    for c in ('on', 'off'):
                        m = fam[c] if how == 'as_is' else E106.at_data_values(fam[c], dt[t][c], odors, how)
                        row[f'rho_{c}'] = float(E106.rho_pairs(m, dt[t][c]))
                    rows.append(row)
    return pd.DataFrame(rows)


@contextlib.contextmanager
def remy_copy():
    """E107's readers reading Remy's files from the copy in data/remy_validation2_kc/ (E107's `data_dir`) instead of her
    folder: E107's `remy_dir` points there while the block runs (its lock checks the copies' sha256); E107's caches are
    cleared on entry and exit."""
    from experiments import E107_validation2_panel as E107
    if not (E107.data_dir() / 'README.md').exists():
        raise RuntimeError(f"no copy of Remy's files in {E107.data_dir()} (E107's part data)")
    keep = E107.remy_dir
    E107.remy_dir = E107.data_dir
    E107._LOCK_CACHE.clear()
    E107._CACHE.clear()
    try:
        yield E107
    finally:
        E107.remy_dir = keep
        E107._LOCK_CACHE.clear()
        E107._CACHE.clear()


def e107_reference(wirings: Sequence[str]) -> Dict[str, Dict[str, pd.Series]]:
    """E107's saved validation2 runs as fitted on `wirings` ({wiring: {'on', 'off'}} cosine pair distances on spike
    counts, canonical pairs), and 'pop' (their mean)."""
    from experiments import E106_rho_gap_parts as E106
    from experiments import E107_validation2_panel as E107
    runs = E107.family_runs('validation2', 'fitted', list(wirings))
    per = {w: {c: H.canonical_pairs(H.pair_distances(o[c][1], 'cosine')) for c in ('on', 'off')} for w, o in runs.items()}
    return dict(per, pop=E106.mean_pairs(list(per.values())))


def remy_readings(pop_on: pd.Series) -> dict:
    """A validation2 population's APL-on distances against Remy's validation2 somata as E107 reads the family (her mean
    table, cosine; the carrying set's pairs and its 31 of paa and B-myr), read from the copy (`remy_copy`)."""
    with remy_copy() as E107:
        data = E107.remy_v2('table_cos')
        key = E107.carrying_set()
        sp = E107.set_pairs(key)
        pb = E107.subset_of(sp, 'pb')
        r_all, n_all = E107.rho(pop_on, data, sp)
        r_pb, n_pb = E107.rho(pop_on, data, pb)
        key_pb = E107.pkey(*PAA_BMYR)
        cp = H.canonical_pairs(pop_on)
        d_pb = float(cp.loc[key_pb]) if key_pb in cp.index else np.nan
        remy_pb = float(data.loc[key_pb])
    return dict(set=key, rho_all=float(r_all), n_all=int(n_all), rho_pb=float(r_pb), n_pb=int(n_pb), d_paa_bmyr_on=d_pb,
                remy_d_paa_bmyr=remy_pb)


def references_table() -> pd.DataFrame:
    """The references (part 0): the family at the lab's input (E108's `unit_scores` on all eight, S, T and each wiring,
    raw and calibrated: 'I0'); E108's carrying input on T when its runs are there ('E108ref', beside)."""
    from experiments import E108_paa_extrapolation as E108
    rows = []
    for ro in ('raw', 'cal'):
        for unit in [u for u, _ in UNITS] + list(ALL_WIRINGS):
            rows.append(dict(input='I0', unit=unit, kind='population' if unit in dict(UNITS) else 'wiring', readout=ro,
                             **E108.unit_scores(family_population(unit, ro))))
    ref = e108_reference_runs()
    if ref:
        for ro in ('raw', 'cal'):
            ds = {w: E83.pair_distances(o, ro) for w, o in ref['runs'].items()}
            for unit in ['T'] + list(T_WIRINGS):
                d = ds[unit] if unit in ds else E83.population(list(ds.values()))
                rows.append(dict(input='E108ref', unit=unit, kind='wiring' if unit in ds else 'population', readout=ro,
                                 rule=ref['rule'], **E108.unit_scores(d)))
    return pd.DataFrame(rows)


def e108_reference_runs() -> Optional[dict]:
    """E108's held runs of its carrying rule on T (B_<rule>_<wiring>), when all four are there: {'rule', 'runs'}."""
    from experiments import E108_paa_extrapolation as E108
    st = e108_class()
    rule = st.get('carrying_rule')
    if not rule:
        return None
    runs = {w: E108.load_run(E108.run_tag(rule, w)) for w in T_WIRINGS}
    if any(o is None for o in runs.values()):
        return None
    return dict(rule=rule, runs=runs)


def part0(threads: int = WORKER_THREADS) -> dict:
    """Part 0 (the module docstring): the identities, the bound, the references, the yardsticks; part0.json, bound.csv,
    references.csv, yardsticks.csv, probes.csv, e107_reference.csv. Raises if an identity fails."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E107_validation2_panel as E107
    from experiments import E108_paa_extrapolation as E108
    cc = commits_clean()
    out = out_dir()
    res = dict(script=SCRIPT, generated=_now(), host=E107._host(), **cc, objective={}, search_path={}, reading_path={},
               step7_path={}, e108=e108_class())
    lab = H.apl16_orn_deltas()
    # the objective at identity
    saved = {w: E83.saved_fit(w)['saved'] for w in ALL_WIRINGS}
    for h in HALVES:
        v = RhoObjective(h)(on={w: saved[w]['on'][1] for w in S_WIRINGS}, off={w: saved[w]['off'][1] for w in S_WIRINGS})
        want = E83.IDENTITY_OBJECTIVE.get(h)
        res['objective'][h] = dict(v, pattern_probe=want, ok=True if want is None else bool(abs(v['pattern'] - want) <= 5e-5))
    print('E109 0 objective at identity: ' + json.dumps(res['objective']), flush=True)
    # the search's path: the pseudo-odor runs at factor 1
    pool = SR.InputPool(S_WIRINGS, threads)
    try:
        for odor in SEARCH_ODORS + CONTROL_ODORS:
            r = pool.run(odor, np.ones((1, len(lab.index))))
            oc = SR._odor_col(saved['real']['on'][1], odor)
            diffs = {}
            for w in S_WIRINGS:
                n = 0
                for c, kr, kc in (('on', 0, 1), ('off', 2, 3)):
                    n += int((r[w][kc][:, 0].astype(float) != saved[w][c][1][oc].to_numpy(dtype=float)).sum())
                    n += int((r[w][kr][:, 0].astype(float) != saved[w][c][0][oc].to_numpy(dtype=float)).sum())
                diffs[w] = n
            res['search_path'][odor] = diffs
            print(f'E109 0 search path {odor}: differences {diffs}', flush=True)
    finally:
        pool.close()
    # the reading's path on T
    for w in T_WIRINGS:
        t0 = time.time()
        sim, info = E85.build(w, lab, cache_tag(lab, 'I0'), cache_dir=local_dir())
        own = E83.wiring_parts(w)['breadth']
        b = info['breadth']
        beq = bool(b.index.equals(own.index) and np.array_equal(b.to_numpy(dtype=float), own.to_numpy(dtype=float)))
        p = E83.saved_fit(w)['p']
        o = {'on': sim(dict(p)), 'off': sim(dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))}
        ns, nr = E83._same(o, saved[w])
        del sim, info
        gc.collect()
        hr = held(w, lab, 'I0')
        hs, hn = E83._same({'on': hr['on'], 'off': hr['off']}, saved[w])
        res['reading_path'][w] = dict(breadth_equal=beq, count_diffs=ns, response_diffs=nr, held_count_diffs=hs,
                                      held_response_diffs=hn, held_thr_factor=hr['thr_factor'], held_apl_factor=hr['apl_factor'],
                                      seconds=round(time.time() - t0, 1))
        print(f"E109 0 reading path {w}: {res['reading_path'][w]}", flush=True)
        del hr
        gc.collect()
    # step 7's path on T: E107's as-fitted validation2 runs
    v2 = E107.panel_input('validation2')
    for w in T_WIRINGS:
        ref = E107.load_run(E107.run_tag('validation2', 'fitted', w))
        if ref is None:
            res['step7_path'][w] = dict(error=f"E107's run {E107.run_tag('validation2', 'fitted', w)} not found", count_diffs=-1,
                                        response_diffs=-1)
            continue
        o = fitted_on(w, v2)
        ns, nr = E83._same({'on': o['on'], 'off': o['off']}, ref)
        res['step7_path'][w] = dict(count_diffs=ns, response_diffs=nr)
        print(f"E109 0 step 7 path {w}: {res['step7_path'][w]}", flush=True)
    # the bound, the references, E107's reference, the yardsticks, the probes
    bt = bound_table()
    bt.to_csv(out / 'bound.csv', index=False)
    rt = references_table()
    rt.to_csv(out / 'references.csv', index=False)
    e107 = []
    for unit, ws in (('all8', ALL_WIRINGS), ('T', T_WIRINGS)):
        ref = e107_reference(ws)
        e107.append(dict(unit=unit, reading='fitted', **remy_readings(ref['pop']['on'])))
    e107 = pd.DataFrame(e107)
    e107.to_csv(out / 'e107_reference.csv', index=False)
    y = E108.e106_yardsticks()
    gates = pd.DataFrame([dict(name=f'gate_{k}', value=v, source=SCRIPT, definition=d) for k, v, d in (
        ('gain_on', GAIN['on'], 'the gain in ρ on odor-odor distances with APL on that the gate asks of T against the unseen half'),
        ('gain_off', GAIN['off'], 'the same with APL off'),
        ('other_slack', OTHER_SLACK, 'the other ρ not lower by more than this'),
        ('carry_mean', CARRY_MEAN, 'the mean gain of "carries below the gate"'))])
    pd.concat([y, gates], ignore_index=True).to_csv(out / 'yardsticks.csv', index=False)
    probes = []
    for (kind, unit, which, how, t), (pon, poff) in PROBES.items():
        if kind == 'bound':
            r = bt[(bt.unit == unit) & (bt.odors == which) & (bt.how == how) & (bt.target == t)].iloc[0]
            got_on, got_off = float(r.rho_on), float(r.rho_off)
        else:
            r = e107[e107.unit == unit].iloc[0]
            got_on, got_off = float(r.rho_all if which == '136' else r.rho_pb), np.nan
        probes.append(dict(kind=kind, unit=unit, odors_or_pairs=which, how=how, target=t, probe_on=pon, value_on=got_on,
                           probe_off=poff, value_off=got_off))
    pr = pd.DataFrame(probes)
    pr['agrees_to_0.001'] = [bool(abs(a - b) <= 0.0015 and (np.isnan(c) or abs(c - d) <= 0.0015))
                             for a, b, c, d in zip(pr.probe_on, pr.value_on, pr.probe_off, pr.value_off)]
    pr.to_csv(out / 'probes.csv', index=False)
    res['probes_agree'] = bool(pr['agrees_to_0.001'].all())
    ok = (all(v['ok'] for v in res['objective'].values())
          and all(n == 0 for d in res['search_path'].values() for n in d.values())
          and all(d['breadth_equal'] and d['count_diffs'] == 0 and d['response_diffs'] == 0 and d['held_count_diffs'] == 0
                  and d['held_response_diffs'] == 0 and d['held_thr_factor'] == 1.0 and d['held_apl_factor'] == 1.0
                  for d in res['reading_path'].values())
          and all(d['count_diffs'] == 0 and d['response_diffs'] == 0 for d in res['step7_path'].values()))
    res['all_exact'] = bool(ok)
    _write_json(out / 'part0.json', res)
    print(json.dumps({k: v for k, v in res.items() if k in ('objective', 'all_exact', 'probes_agree')}, indent=1, default=str),
          flush=True)
    print(pr.to_string(index=False), flush=True)
    if not ok:
        raise RuntimeError('E109 part 0: an identity does not hold (part0.json); no search starts')
    return res


# ---------------------------------------------------------------------------
# the held readings (part read)

def change_meta(names: Optional[Sequence[str]] = None) -> Dict[str, dict]:
    """{change: dict(search, kind ('search', 'perm', 'inverted'), size, seed, half, unseen, split, group)} of the
    readings of the finished searches (`names`, default all): each at GATE_SIZE and FULL_SIZE; for the six, the
    permuted factors at GATE_SIZE (PERM_SEEDS) and the inverted at INVERT_SIZE."""
    out = {}
    for n in (names or SEARCHES):
        if not search_finished(n):
            continue
        cfg = SEARCHES[n]
        m = dict(search=n, half=cfg['half'], unseen=UNSEEN[cfg['half']], split=SPLIT_OF[cfg['half']], group=cfg['group'])
        out[f'{n}_s{GATE_SIZE:g}'] = dict(m, kind='search', size=GATE_SIZE, seed=None)
        out[f'{n}_s{FULL_SIZE:g}'] = dict(m, kind='search', size=FULL_SIZE, seed=None)
        if n in MAIN:
            for s in PERM_SEEDS:
                out[f'{n}_perm{s}'] = dict(m, kind='perm', size=GATE_SIZE, seed=s)
            out[f'{n}_inv'] = dict(m, kind='inverted', size=INVERT_SIZE, seed=None)
    return out


def change_input(change: str, meta: Optional[Mapping] = None) -> pd.DataFrame:
    """The ORN input of a reading: the search's start input times its factors (permuted for 'perm') at the size."""
    m = meta or change_meta([change.split('_')[0]])[change]
    F = search_factors(m['search'])
    if m['kind'] == 'perm':
        F = E83.permuted_factors(F, int(m['seed']))
    return E83.input_at(F, m['size'], orn=base_of(search_start(m['search'])))


def read_plan(names: Optional[Sequence[str]] = None) -> List[dict]:
    """The held readings to run, wiring by wiring: [dict(change, wiring, tag, **meta)]; the gate's inputs (0.6) first,
    then the controls, then the full sizes."""
    meta = change_meta(names)
    rank = lambda ch: (0 if meta[ch]['kind'] == 'search' and meta[ch]['size'] == GATE_SIZE else
                       1 if meta[ch]['kind'] != 'search' else 2)
    order = sorted(meta, key=lambda ch: (rank(ch), list(meta).index(ch)))
    return [dict(change=ch, wiring=w, tag=f'{ch}_{w}', **meta[ch]) for w in T_WIRINGS for ch in order]


def part_read(names: Optional[Sequence[str]] = None, max_points: Optional[int] = None) -> pd.DataFrame:
    """The held readings (resumable: points in runs_R.csv are skipped); each wiring's simulator is first checked against
    its saved fit."""
    from experiments.E77_input_screen import frame_sha256
    if not part0_passed():
        raise RuntimeError('run --part 0 first')
    cc = commits_clean()
    out = out_dir()
    rp, pp = out / 'runs_R.csv', out / 'provenance_R.json'
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    todo = [p for p in read_plan(names) if p['tag'] not in done]
    if max_points is not None:
        todo = todo[:max_points]
    print(f'E109 read: {len(todo)} points to run', flush=True)
    for w in [w for w in T_WIRINGS if any(p['wiring'] == w for p in todo)]:
        sim, info = E83.build_simulator(w)
        E83._check_saved(w, sim, info['params'], E83.saved_fit(w)['saved'])
        del sim, info
        gc.collect()
        for pt in [p for p in todo if p['wiring'] == w]:
            orn = change_input(pt['change'], pt)
            path = E83._write_input(out, pt['change'], orn)
            res = held(w, orn, pt['change'])
            o = {'on': res['on'], 'off': res['off']}
            save_run('R', pt['tag'], o)
            row = E83._run_row('R', pt['change'], w, pt['tag'], res, o, search=pt['search'], kind=pt['kind'],
                               seed=pt['seed'], half=pt['half'], unseen=pt['unseen'], size=pt['size'], group=pt['group'],
                               start=search_start(pt['search'])['start'], input_tag=res['input_tag'],
                               orn_input=str(path.relative_to(H.shared_dir())), input_frame_sha256=frame_sha256(orn),
                               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                               params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}))
            row.update(run_extras(o))
            row.update({k: v for k, v in cc.items() if k == 'queue_dirty_paths'})
            E83._append_row(rp, row)
            _provenance_add(pp, 'R', pt['tag'], H.model_provenance(orn_deltas=orn))
            print(f"E109 R {pt['tag']}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/"
                  f"{res['target_on']} on, {res['n_off']}/{res['target_off']} off); flag '{row.get('realism_flag') or ''}'; "
                  f"{res['seconds']:.0f}s", flush=True)
            del res, o
            gc.collect()
    return E83._read_rows(rp)


def run_extras(o: Dict[str, tuple]) -> dict:
    """The response fractions of the five odors searched, APL on and off (two or more spikes)."""
    row = {}
    for c in ('on', 'off'):
        r = H.abbrev_columns(o[c][0].astype(float) > 0.5)
        for od in SEARCH_ODORS + CONTROL_ODORS:
            row[f'frac_{od}_{c}'] = float(r[od].mean())
    return row


# ---------------------------------------------------------------------------
# step 7 (simulations)

def step7_input(name: str, size: float = GATE_SIZE) -> pd.DataFrame:
    """The validation2 ORN input (E107's `panel_input`) with its own paa (10^-4) and B-myr times search `name`'s factors
    for those odors at `size` of their log size."""
    from experiments import E107_validation2_panel as E107
    F = search_factors(name)
    v2 = E107.panel_input('validation2').astype(float).copy()
    for c in v2.columns:
        o = H.odor_abbrev(c)
        if o in PAA_BMYR:
            f = F[o].reindex(v2.index).to_numpy(dtype=float)
            v2[c] = v2[c].to_numpy(dtype=float) * np.power(f, float(size))
    return v2


def part_step7(names: Sequence[str] = MAIN) -> pd.DataFrame:
    """Step 7's runs (resumable): each finished search of `names` on the four wirings of T as fitted (`fitted_on`) with
    `step7_input`; rows in runs_V.csv."""
    from experiments.E77_input_screen import frame_sha256
    if not part0_passed():
        raise RuntimeError('run --part 0 first')
    cc = commits_clean()
    out = out_dir()
    rp, pp = out / 'runs_V.csv', out / 'provenance_V.json'
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    for n in [n for n in names if search_finished(n)]:
        orn = step7_input(n)
        for w in T_WIRINGS:
            tag = f'V_{n}_{w}'
            if tag in done:
                continue
            t0 = time.time()
            o = fitted_on(w, orn)
            save_run('V', tag, {'on': o['on'], 'off': o['off']})
            row = dict(part='V', search=n, wiring=w, tag=tag, size=GATE_SIZE, panel='validation2', reading='fitted',
                       input_frame_sha256=frame_sha256(orn), fit_folder=E83.saved_fit(w)['folder'],
                       fit_tag=E83.saved_fit(w)['tag'], **{k: o['params'][k] for k in E83.PARAM_KEYS},
                       rate_on=float(o['on'][0].to_numpy(dtype=float).mean()),
                       rate_off=float(o['off'][0].to_numpy(dtype=float).mean()), seconds=time.time() - t0, **cc)
            E83._append_row(rp, row)
            _provenance_add(pp, 'V', tag, H.model_provenance(orn_deltas=orn))
            print(f'E109 V {tag}: {time.time() - t0:.0f}s', flush=True)
    return E83._read_rows(rp)


# ---------------------------------------------------------------------------
# scoring (no simulations)

def runs_table() -> pd.DataFrame:
    """Every reading run on T: the saved fits as 'I0', E108's carrying input as 'E108ref' when its runs are there, and
    runs_R.csv's; with 'flagged'."""
    rows = [dict(change='I0', wiring=w, tag=f'I0_{w}', kind='reference', **E83.realism(E83.saved_fit(w)['saved']),
                 **run_extras(E83.saved_fit(w)['saved'])) for w in T_WIRINGS]
    ref = e108_reference_runs()
    if ref:
        rows += [dict(change='E108ref', wiring=w, tag=f'E108ref_{w}', kind='reference', **E83.realism(o), **run_extras(o))
                 for w, o in ref['runs'].items()]
    rr = E83._read_rows(out_dir() / 'runs_R.csv')
    df = pd.concat([pd.DataFrame(rows), rr], ignore_index=True, sort=False) if len(rr) else pd.DataFrame(rows)
    df = df.drop_duplicates('tag', keep='last').reset_index(drop=True)
    df['flagged'] = [_flagged(x) for x in df['realism_flag']]
    return df


def run_outputs(change: str, wiring: str) -> Dict[str, tuple]:
    if change == 'I0':
        return E83.saved_fit(wiring)['saved']
    if change == 'E108ref':
        return e108_reference_runs()['runs'][wiring]
    o = load_run('R', f'{change}_{wiring}')
    if o is None:
        raise RuntimeError(f'no saved run {change}_{wiring}')
    return o


RUN_MEANS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
             'max_count_off', 'thr_factor', 'apl_factor') + tuple(f'frac_{o}_{c}' for o in SEARCH_ODORS + CONTROL_ODORS
                                                                   for c in ('on', 'off'))


def readings(runs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """readings.csv: per input (change), unit (each wiring of T and T's population) and readout ('raw', 'cal'), E108's
    `unit_scores` (the two ρ's against all recordings and each half, the APL-effect statistics, the block, the
    natural-source alignment and the structure vs wild-type somata), the mean over the unit's runs of the realism
    quantities and the runs flagged; 'paa_bmyr_on' / 'paa_bmyr_off' (the distance of paa from B-myr)."""
    from experiments import E108_paa_extrapolation as E108
    runs = runs_table() if runs is None else runs
    rows = []
    for ch, g in runs.groupby('change', sort=False):
        if set(g['wiring']) != set(T_WIRINGS):
            continue
        meta = {r.wiring: r for r in g.itertuples()}
        for ro in ('raw', 'cal'):
            ds = {w: E83.pair_distances(run_outputs(ch, w), ro) for w in T_WIRINGS}
            for unit in list(T_WIRINGS) + ['T']:
                ws = T_WIRINGS if unit == 'T' else (unit,)
                d = ds[unit] if unit in ds else E83.population([ds[w] for w in ws])
                sub = g[g['wiring'].isin(ws)]
                pb = E83._cs(d['on']).index.get_loc(tuple(sorted(PAA_BMYR)))
                rows.append(dict(change=ch, unit=unit, kind='wiring' if unit in ds else 'population', readout=ro,
                                 n_flagged=int(sum(meta[w].flagged for w in ws)),
                                 **{f'run_{k}': float(sub[k].astype(float).mean()) for k in RUN_MEANS if k in sub},
                                 paa_bmyr_on=float(E83._cs(d['on']).iloc[pb]), paa_bmyr_off=float(E83._cs(d['off']).iloc[pb]),
                                 **E108.unit_scores(d)))
    rd = pd.DataFrame(rows)
    out = out_dir()
    rd.to_csv(out / 'readings.csv', index=False)
    return rd


def _v(rd: pd.DataFrame, change: str, unit: str, col: str, readout: str = 'raw') -> float:
    r = rd[(rd['change'] == change) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _d(rd: pd.DataFrame, change: str, ref: str, unit: str, col: str, readout: str = 'raw') -> float:
    return _v(rd, change, unit, col, readout) - _v(rd, ref, unit, col, readout)


def _reading(name: str, kind: str = 'search', size: float = GATE_SIZE, seed: Optional[int] = None) -> str:
    if kind == 'search':
        return f'{name}_s{size:g}'
    if kind == 'perm':
        return f'{name}_perm{seed}'
    if kind == 'inverted':
        return f'{name}_inv'
    raise ValueError(kind)


def direction_values(rd: pd.DataFrame, runs: pd.DataFrame, names: Sequence[str], ref: str = 'I0', kind: str = 'search',
                     size: float = GATE_SIZE, seeds: Sequence[Optional[int]] = (None,)) -> Optional[dict]:
    """The values of one direction (the searches `names` of one half, each read against its unseen half), each the mean
    over the searches (and `seeds`): T's change of each ρ, each wiring's change, the APL-effect statistics' changes, the
    runs flagged, the response rate ρ's change, the pending items' changes. None if a reading is missing."""
    chs = [_reading(n, kind, size, s) for n in names for s in seeds]
    have = set(rd['change'])
    if not chs or any(c not in have for c in chs):
        return None
    un = {c: UNSEEN[SEARCHES[n]['half']] for n in names for c in [_reading(n, kind, size, s) for s in seeds]}
    m = lambda f: float(np.mean([f(c, un[c]) for c in chs]))
    out = dict(changes=chs, unseen=sorted(set(un.values())))
    for c in ('on', 'off'):
        out[f'd_rho_{c}'] = m(lambda ch, y: _d(rd, ch, ref, 'T', f'rho_{c}_{y}'))
        out[f'd_rho_{c}_wirings'] = {w: m(lambda ch, y: _d(rd, ch, ref, w, f'rho_{c}_{y}')) for w in T_WIRINGS}
        out[f'n_wirings_gain_{c}'] = int(sum(v >= GAIN[c] for v in out[f'd_rho_{c}_wirings'].values()))
        out[f'd_rho_{c}_all'] = m(lambda ch, y: _d(rd, ch, ref, 'T', f'rho_{c}_all'))
    for k in ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson_effect', 'rho_effect', 'loss'):
        out[f'd_{k}_cal'] = m(lambda ch, y: _d(rd, ch, ref, 'T', f'{y}_{k}', 'cal'))
    out['d_pattern_all_cal'] = m(lambda ch, y: _d(rd, ch, ref, 'T', 'all_pattern', 'cal'))
    out['d_loss_raw'] = m(lambda ch, y: _d(rd, ch, ref, 'T', f'{y}_loss', 'raw'))
    out['d_loss_raw_all'] = m(lambda ch, y: _d(rd, ch, ref, 'T', 'all_loss', 'raw'))
    out['d_mean_effect_cal'] = m(lambda ch, y: _d(rd, ch, ref, 'T', 'mean_effect', 'cal'))
    for k in ('block', 'vcf_on', 'vcf_off', 'struct_wt_on', 'struct_wt_off', 'paa_bmyr_on', 'paa_bmyr_off'):
        out[f'd_{k}'] = m(lambda ch, y: _d(rd, ch, ref, 'T', k))
    out['d_vcf_on_minus_off'] = out['d_vcf_on'] - out['d_vcf_off']
    r = runs[runs['change'].isin(chs) & runs['wiring'].isin(T_WIRINGS)]
    r0 = runs[(runs['change'] == ref) & runs['wiring'].isin(T_WIRINGS)]
    out['n_flagged'] = int(r['flagged'].sum())
    out['n_runs'] = int(len(r))
    for c in ('on', 'off'):
        out[f'd_rate_rho_{c}'] = float(r[f'rate_rho_{c}'].astype(float).mean() - r0[f'rate_rho_{c}'].astype(float).mean())
    return out


def gate(rd: pd.DataFrame, runs: pd.DataFrame, ref: str = 'I0', main: Sequence[str] = FIRST_SPLIT,
         second: Sequence[str] = SECOND_SPLIT, size: float = GATE_SIZE, with_controls: bool = True) -> dict:
    """The entry's gate (the module docstring) on the readings `rd` (`readings`' layout) and the runs `runs`
    (`runs_table`'s: change, wiring, flagged, rate_rho_on / _off) against the reference `ref`, for the searches `main`
    (two directions of E77's split, by their halves) and `second` (E86's split). Returns a dict with the directions'
    values, the gain per ρ, the controls, the failures ('fails', 'pending_fails'), 'complete', 'cls',
    'class_without_pending', 'outcome'."""
    out = dict(reference=ref, size=size, searches=list(main), complete=False)
    dirs = {}
    for X in sorted({SEARCHES[n]['half'] for n in main}):
        dv = direction_values(rd, runs, [n for n in main if SEARCHES[n]['half'] == X], ref, 'search', size)
        if dv is None:
            out['note'] = f'readings missing (direction {X})'
            return out
        dirs[X] = dv
    if len(dirs) != 2:
        out['note'] = f'{len(dirs)} directions'
        return out
    sv = direction_values(rd, runs, list(second), ref, 'search', size) if second else None
    out['complete'] = True
    out['directions'] = dirs
    out['second'] = sv
    other = {'on': 'off', 'off': 'on'}
    for c in ('on', 'off'):
        ok_dirs = all(d[f'd_rho_{c}'] >= GAIN[c] and d[f'n_wirings_gain_{c}'] >= GAIN_WIRINGS for d in dirs.values())
        ok_other = all(d[f'd_rho_{other[c]}'] >= -OTHER_SLACK for d in dirs.values())
        # no second split (the control searches' own reading): not applicable
        ok_second = True if not second else bool(sv is not None and sv[f'd_rho_{c}'] > 0)
        out[f'mean_gain_{c}'] = float(np.mean([d[f'd_rho_{c}'] for d in dirs.values()]))
        out[f'second_d_rho_{c}'] = None if sv is None else float(sv[f'd_rho_{c}'])
        out[f'gain_{c}_parts'] = dict(directions=bool(ok_dirs), other=bool(ok_other), second_split=ok_second)
        out[f'gain_{c}'] = bool(ok_dirs and ok_other and ok_second)
    out['gain'] = bool(out['gain_on'] or out['gain_off'])
    # the controls, on each ρ, against the candidate's mean gain (None: not read, the control searches' own rows)
    ctl = {} if with_controls else None
    if with_controls:
        for c in ('on', 'off'):
            g0 = out[f'mean_gain_{c}']
            vals = {}
            for label, kind, seeds, names in (('permuted', 'perm', PERM_SEEDS, main), ('inverted', 'inverted', (None,), main),
                                              ('control search', 'search', (None,), CONTROL_SEARCHES)):
                per = []
                for X in sorted({SEARCHES[n]['half'] for n in names}):
                    dv = direction_values(rd, runs, [n for n in names if SEARCHES[n]['half'] == X], ref, kind, size, seeds)
                    per.append(np.nan if dv is None else dv[f'd_rho_{c}'])
                vals[label] = float(np.mean(per)) if per and np.isfinite(per).all() else np.nan
            ok = bool(g0 > 0 and all(np.isfinite(v) and v < CONTROL_RATIO * g0 for v in vals.values()))
            ctl[c] = dict(candidate=g0, **vals, ok=ok)
    out['controls'] = ctl
    # the other measures (each direction; the worst named)
    fails, pend = {}, {}

    def worst(key, sign):
        vs = [d[key] for d in dirs.values()]
        return float(max(vs) if sign > 0 else min(vs))
    checks = (('pattern part (calibrated, unseen half)', 'd_pattern_cal', +1, EQUAL_PATTERN, fails),
              ('raw-count loss (unseen half)', 'd_loss_raw', +1, EQUAL_RAW, fails),
              ('mean APL effect (calibrated)', 'd_mean_effect_cal', -1, EQUAL_MEAN, fails),
              ('response rate ρ, APL on', 'd_rate_rho_on', -1, RATE_RHO_SLACK, fails),
              ('response rate ρ, APL off', 'd_rate_rho_off', -1, RATE_RHO_SLACK, fails),
              ('the block', 'd_block', -1, BLOCK_SLACK, pend),
              ('natural-source alignment, APL on', 'd_vcf_on', -1, STRUCTURE_SLACK, pend),
              ('natural-source alignment, APL off', 'd_vcf_off', -1, STRUCTURE_SLACK, pend),
              ('structure vs wild-type somata', 'd_struct_wt_on', -1, STRUCTURE_SLACK, pend))
    for label, key, sign, slack, into in checks:
        w = worst(key, sign)
        out[f'worst_{key}'] = w
        if (sign > 0 and w > slack) or (sign < 0 and w < -slack):
            into[label] = w
    nf = int(sum(d['n_flagged'] for d in dirs.values()))
    out['n_flagged'] = nf
    if nf:
        fails['realism (runs outside the bands)'] = nf
    out['fails'], out['pending_fails'] = fails, pend
    out['cls'] = class_of(out)
    out['class_without_pending'] = class_of(out, with_pending=False)
    out['outcome'] = OUTCOME[out['cls']]
    return out


def class_of(g: Mapping, with_pending: bool = True) -> str:
    """The class of a gate's dict: 'passes' (the gain and its controls hold and nothing fails), 'a price' (the gain and its
    controls hold, a measure fails), 'carries below the gate' (otherwise, a mean gain of CARRY_MEAN or more on a ρ whose
    controls hold), 'does not carry'. The pending items count unless `with_pending` is False."""
    if not g.get('complete'):
        raise ValueError('the gate is not complete')
    ctl = g.get('controls')
    ctl_ok = lambda c: True if ctl is None else bool(ctl.get(c, {}).get('ok'))
    gaining = [c for c in ('on', 'off') if g.get(f'gain_{c}')]
    if gaining and all(ctl_ok(c) for c in gaining):
        failed = bool(g.get('fails')) or (with_pending and bool(g.get('pending_fails')))
        return 'a price' if failed else 'passes'
    if any(g.get(f'mean_gain_{c}', -np.inf) >= CARRY_MEAN and ctl_ok(c) for c in ('on', 'off')):
        return 'carries below the gate'
    return 'does not carry'


def _flat(g: Mapping, **extra) -> dict:
    """A gate's dict as one table row (nested values as JSON)."""
    row = dict(extra)
    for k, v in g.items():
        if k in ('directions', 'second'):
            continue
        if isinstance(v, dict):
            row[k] = '; '.join(f'{a} {b:+.4f}' if isinstance(b, float) else f'{a} {b}' for a, b in v.items()) \
                if k in ('fails', 'pending_fails') else json.dumps(v, default=str)
        elif isinstance(v, (list, tuple)):
            row[k] = ' '.join(map(str, v))
        else:
            row[k] = v
    for X, d in (g.get('directions') or {}).items():
        for k, v in d.items():
            if not isinstance(v, (dict, list)):
                row[f'dir{X}_{k}'] = v
    if g.get('second'):
        for k, v in g['second'].items():
            if not isinstance(v, (dict, list)):
                row[f'second_{k}'] = v
    return row


def gate_tables(rd: pd.DataFrame, runs: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """gate.csv (the six against the lab's input, at 0.6 (carrying) and at 1 (beside); against E108's input beside under
    "better"; the control searches read the same way) and gate.json (the carrying verdict)."""
    rows, full = [], {}
    refs = ['I0'] + (['E108ref'] if (rd['change'] == 'E108ref').any() else [])
    for ref in refs:
        for size in (GATE_SIZE, FULL_SIZE):
            g = gate(rd, runs, ref=ref, size=size)
            key = f'six:{ref}:s{size:g}'
            full[key] = g
            rows.append(_flat(g, which='six', carrying=bool(ref == 'I0' and size == GATE_SIZE)))
        g = gate(rd, runs, ref=ref, main=CONTROL_SEARCHES, second=(), size=GATE_SIZE, with_controls=False)
        full[f'ctl:{ref}:s{GATE_SIZE:g}'] = g
        rows.append(_flat(g, which='ctl', carrying=False))
    gt = pd.DataFrame(rows)
    gt.to_csv(out_dir() / 'gate.csv', index=False)
    v = full.get(f'six:I0:s{GATE_SIZE:g}', {})
    verdict = dict(complete=bool(v.get('complete')))
    if v.get('complete'):
        verdict.update(cls=v['cls'], class_without_pending=v['class_without_pending'], outcome=v['outcome'], gain_on=v['gain_on'],
                       gain_off=v['gain_off'], mean_gain_on=v['mean_gain_on'], mean_gain_off=v['mean_gain_off'],
                       controls=v['controls'], fails=v['fails'], pending_fails=v['pending_fails'])
    _write_json(out_dir() / 'gate.json', dict(verdict=verdict, gates=full, **H.code_commits()))
    return gt, verdict


def searches_table() -> pd.DataFrame:
    """searches.csv: per finished (or running) search, its start and last parent on S (the two ρ's and the pattern part
    against the searched and the unseen half), per odor the rms of its log2 factors and the factors at a bound."""
    rows = []
    for n, cfg in SEARCHES.items():
        p = out_dir() / f'search_history_{n}.csv'
        if not p.exists():
            continue
        h = pd.read_csv(p)
        s0 = h[h.odor == 'start']
        last = h[h.odor != 'start'].tail(1)
        row = dict(search=n, half=cfg['half'], unseen=UNSEEN[cfg['half']], odors=' '.join(cfg['odors']),
                   finished=search_finished(n), generations_done=int((h.odor != 'start').sum()))
        for k in ('rho_on_half', 'rho_off_half', 'pattern_half', 'rho_on_unseen', 'rho_off_unseen', 'pattern_unseen'):
            row[f'start_{k}'] = float(s0[k].iloc[0]) if len(s0) else np.nan
            row[f'last_{k}'] = float(last[k].iloc[0]) if len(last) else np.nan
        if search_finished(n):
            f = pd.read_csv(factors_path(n))
            for o, g in f.groupby('odor', sort=False):
                row[f'{o}_log2_rms'] = float(np.sqrt(np.mean(g.log2_factor ** 2)))
                row[f'{o}_n_at_bound'] = int((np.abs(np.abs(g.log2_factor) - 1.0) < 1e-9).sum())
        rows.append(row)
    st = pd.DataFrame(rows)
    st.to_csv(out_dir() / 'searches.csv', index=False)
    return st


# ---------------------------------------------------------------------------
# step 6: the entries that recur across the six searches (no simulations)

def validation2_fly_values(odors: Sequence[str] = PAA_BMYR) -> Dict[str, np.ndarray]:
    """The six validation2 ORN flies' responses (olfdb dataset 4, E108's compact copy) to `odors` at the lab input's
    concentrations, in the model's units as E108's `bmyr_flies` puts them (each fly scaled to-avg-max, then E108's slope
    through the origin of the lab's four validation odors on the flies' mean): {odor: flies x the lab input's 36
    glomeruli, NaN where a fly lacks the glomerulus}."""
    from experiments import E108_paa_extrapolation as E108
    B = E108.bmyr_flies()
    t = E108.compact_tables()
    gl = list(B['glomeruli'])
    fm = t['validation'][t['validation']['panel'] == 'validation2']
    m = fm.merge(t['validation_maxima'], on=['date', 'fly_num'], validate='many_to_one')
    A = E108.dff_fit()['avg_flymax']
    m = m.assign(scaled=m['value'] / m['fly_max'] * A, odor=[o.split(' @ ')[0] for o in m['odor1']])
    out = {}
    for od in odors:
        want = f'{od} @ {H.ORN_INPUT_CONC[od]:g}'
        b = m[(m['odor1'] == want) & m['roi'].isin(gl)]
        V = np.full((len(B['keys']), len(gl)), np.nan)
        for (d, f, g), v in b.groupby(['date', 'fly_num', 'roi'])['scaled'].mean().items():
            V[B['keys'].index((str(d), int(f))), gl.index(g)] = v * B['slope']
        out[od] = V
    return out


def spread_rows(V: np.ndarray, glomeruli: Sequence[str], odor: str, source: str) -> List[dict]:
    """Per glomerulus: the flies' mean, the number of flies, and the SD across flies of log2(fly / mean) (flies and mean
    positive, at least 3 flies; E85's `spread_table` definition)."""
    rows = []
    with np.errstate(invalid='ignore'):
        M = np.nanmean(V, axis=0)
    for k, g in enumerate(glomeruli):
        v = V[:, k]
        v = v[np.isfinite(v)]
        pos = v[v > 0]
        sdl = float(np.std(np.log2(pos / M[k]), ddof=1)) if (np.isfinite(M[k]) and M[k] > 0 and len(pos) >= 3) else np.nan
        rows.append(dict(odor=odor, glomerulus=g, flies_mean=float(M[k]) if np.isfinite(M[k]) else np.nan,
                         n_flies=int(len(v)), n_positive=int(len(pos)), flies_sd_log2=sdl, flies_source=source))
    return rows


def flies_spread() -> pd.DataFrame:
    """The ORN flies' spread for the three odors: paa and B-myr over the six validation2 flies (`validation2_fly_values`),
    2-but over the nine megamat flies (E85's `spread_table`)."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E108_paa_extrapolation as E108
    gl = list(E108.bmyr_flies()['glomeruli'])
    rows = []
    for od, V in validation2_fly_values().items():
        rows += spread_rows(V, gl, od, 'validation2 ORN flies (olfdb dataset 4, six flies; E108\'s units)')
    sp = E85.spread_table()
    for r in sp[sp['odor'].isin([o for o in SEARCH_ODORS if o not in PAA_BMYR])].itertuples():
        rows.append(dict(odor=r.odor, glomerulus=r.glomerulus, flies_mean=r.mean, n_flies=r.n_flies,
                         n_positive=r.n_positive, flies_sd_log2=r.sd_log2_ratio,
                         flies_source="megamat ORN flies (olfdb dataset 1, nine flies; E85's spread_table)"))
    return pd.DataFrame(rows)


def recurring_entries(F: Optional[Mapping[str, pd.DataFrame]] = None, n_null: Optional[int] = None,
                      write: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Step 6 (no simulations): E86's rule (`five_of_six`: five of the six searches agree in sign and the mean log2
    factor is beyond +-0.5; the change as read at 0.6 beside) on the six searches' factors for the three odors, with the
    lab's and the start input, the ORN flies' mean and spread (`flies_spread`) and the change as read in units of that
    spread, and Hallem's log2 ratio to the lab's input where E38 has it (E85's `hallem_log2`); and the rule's counts with
    the factors label-permuted (E86's `five_of_six_null`). Returns (entries, null); with `write`, step6_entries.csv,
    step6_null.csv and step6_flies.csv."""
    from experiments import E85_orn_fly_population as E85
    from experiments import E86_reduced_input_candidate as E86
    F = F or {n: search_factors(n) for n in MAIN}
    if set(F) != set(MAIN):
        raise RuntimeError(f'step 6 needs the six searches ({sorted(set(MAIN) - set(F))} missing)')
    L = {n: np.log2(F[n][list(SEARCH_ODORS)]) for n in MAIN}
    ent = E86.five_of_six(L, scale_read=GATE_SIZE)
    null = E86.five_of_six_null(dict(F), n=E86.N_NULL_E if n_null is None else n_null, odors=SEARCH_ODORS,
                                scale_read=GATE_SIZE)
    lab = H.abbrev_columns(H.apl16_orn_deltas())
    ent['lab_input'] = [float(lab.loc[g, o]) for o, g in zip(ent.odor, ent.glomerulus)]
    fs = flies_spread()
    ent = ent.merge(fs, on=['odor', 'glomerulus'], how='left')
    ent['read_in_flies_sd'] = ent['mean_log2_read'] / ent['flies_sd_log2']
    hal = E85.hallem_log2()
    ent['hallem_log2_vs_lab'] = [float(hal.get((o, g), np.nan)) for o, g in zip(ent.odor, ent.glomerulus)]
    if write:
        ent.to_csv(out_dir() / 'step6_entries.csv', index=False)
        null.to_csv(out_dir() / 'step6_null.csv', index=False)
        fs.to_csv(out_dir() / 'step6_flies.csv', index=False)
        # beside: E108's reading of B-myr within its flies' spread (its step 4, bmyr.csv), when written
        bm = H.results_dir(E108_EXPERIMENT) / 'bmyr.csv'
        if bm.exists():
            pd.read_csv(bm).to_csv(out_dir() / 'step6_e108_bmyr.csv', index=False)
    return ent, null


# ---------------------------------------------------------------------------
# step 7 (tables)

def step7_tables(rd: Optional[pd.DataFrame] = None) -> Dict[str, pd.DataFrame]:
    """Step 7 (beside, no gate): (a) the distance of paa from B-myr in each searched input (ORN input, cosine across the
    glomeruli) and in the family with it on T (raw counts, APL on and off), with the flies' values (control flies, APL
    output blocked (TNT), Remy's wild-type somata) and the family's at the lab's input (all eight and T); (b) each search's
    paa and B-myr factors at 0.6 on the validation2 input, T's population as fitted against Remy's somata (ρ over the
    carrying set's pairs and over the 31 of paa and B-myr), beside E107's T and eight-wiring references. Writes
    step7_paa_bmyr.csv and step7_remy.csv."""
    from experiments import E106_rho_gap_parts as E106
    rd = readings() if rd is None else rd
    pk = tuple(sorted(PAA_BMYR))
    dt = dist_targets()['all']
    flies = dict(ctrl=float(E83._cs(dt['on']).loc[pk]), tnt=float(E83._cs(dt['off']).loc[pk]))
    with remy_copy() as E107:
        flies['remy_wt'] = float(E107.remy_v2('table_cos').loc[E107.pkey(*PAA_BMYR)])
    fam8 = family_population('all8', 'raw')
    rows = [dict(input='I0', size=0.0, orn_cosine=float(H.canonical_pairs(H.pair_distances(H.apl16_orn_deltas(), 'cosine')).loc[PAA_BMYR]),
                 family_all8_on=float(E83._cs(fam8['on']).loc[pk]), family_all8_off=float(E83._cs(fam8['off']).loc[pk]),
                 family_T_on=_v(rd, 'I0', 'T', 'paa_bmyr_on'), family_T_off=_v(rd, 'I0', 'T', 'paa_bmyr_off'), **flies)]
    for n in MAIN:
        if not search_finished(n):
            continue
        for size in (GATE_SIZE, FULL_SIZE):
            ch = _reading(n, 'search', size)
            orn = change_input(ch, change_meta([n])[ch])
            rows.append(dict(input=n, size=size, orn_cosine=float(H.canonical_pairs(H.pair_distances(orn, 'cosine')).loc[PAA_BMYR]),
                             family_T_on=_v(rd, ch, 'T', 'paa_bmyr_on'), family_T_off=_v(rd, ch, 'T', 'paa_bmyr_off'), **flies))
    a = pd.DataFrame(rows)
    a.to_csv(out_dir() / 'step7_paa_bmyr.csv', index=False)
    rows = []
    for unit, ws in (('all8', ALL_WIRINGS), ('T', T_WIRINGS)):
        ref = e107_reference(ws)
        rows.append(dict(input='I0 (E107)', unit=unit, **remy_readings(ref['pop']['on'])))
    for n in MAIN:
        runs = {w: load_run('V', f'V_{n}_{w}') for w in T_WIRINGS}
        if any(o is None for o in runs.values()):
            continue
        per = [{c: H.canonical_pairs(H.pair_distances(o[c][1], 'cosine')) for c in ('on', 'off')} for o in runs.values()]
        pop = E106.mean_pairs(per)
        rows.append(dict(input=n, unit='T', **remy_readings(pop['on'])))
    b = pd.DataFrame(rows)
    b.to_csv(out_dir() / 'step7_remy.csv', index=False)
    return dict(paa_bmyr=a, remy=b)


# ---------------------------------------------------------------------------
# score, the comparison table, the archive tables

def _f(x, spec: str = '.3f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '-'
    except (TypeError, ValueError):
        return str(x)


def _md_table(df: pd.DataFrame, spec: str = '.4f') -> str:
    """A plain markdown table of `df` (numbers formatted with `spec`)."""
    cells = lambda r: [_f(v, spec) if isinstance(v, (float, np.floating)) else str(v) for v in r]
    lines = ['| ' + ' | '.join(map(str, df.columns)) + ' |', '|' + '---|' * len(df.columns)]
    lines += ['| ' + ' | '.join(cells(r)) + ' |' for r in df.itertuples(index=False)]
    return '\n'.join(lines)


def comparison_table(rd: pd.DataFrame, verdict: dict, gates: Mapping[str, dict]) -> str:
    """comparison_table.md: the references, the searches on S, the gate per direction (with the amplitude part, the shape
    part, the spread ratio and Pearson r beside the pattern part), the controls, the class; step 6 and step 7 when
    written."""
    out = out_dir()
    L = ['# E109: the input of paa, B-myr and 2-but searched for the two ρ\'s on held-out wirings and recordings', '',
         'paa\'s input is the lab\'s recording at 10^-4 (the flies\' 10^-3); euc and -aPine stay at the recorded 10^-1.75 and '
         '10^-1.5 (the flies\' 10^-2). Raw counts for the two ρ\'s; the APL-effect statistics under the calibrated readout '
         'unless named. T: the four held-out wirings (P3, P4, M2, M3), their mean matrix.', '']
    st = searches_table()
    if len(st):
        L += ['## The searches on S (the four searched wirings, their mean matrix)', '',
              '| search | half | odors | generations | ρ APL on, searched half (start / last) | ρ APL off, searched half | '
              'ρ APL on, unseen half | ρ APL off, unseen half | pattern part, searched half |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|']
        for r in st.itertuples():
            L.append(f'| {r.search} | {r.half} | {r.odors} | {r.generations_done} | {_f(r.start_rho_on_half)} / {_f(r.last_rho_on_half)} | '
                     f'{_f(r.start_rho_off_half)} / {_f(r.last_rho_off_half)} | {_f(r.start_rho_on_unseen)} / {_f(r.last_rho_on_unseen)} | '
                     f'{_f(r.start_rho_off_unseen)} / {_f(r.last_rho_off_unseen)} | {_f(r.start_pattern_half, ".4f")} / '
                     f'{_f(r.last_pattern_half, ".4f")} |')
        L.append('')
    g = gates.get(f'six:I0:s{GATE_SIZE:g}')
    if g and g.get('complete'):
        L += ['## The gate (T against the unseen half, changes from the family at the lab\'s input; each direction the mean '
              'over its two draws)', '',
              '| statistic | A searched, read on B | B searched, read on A | second split (C1 on D, D1 on C) | threshold |',
              '|---|---:|---:|---:|---|']
        dA, dB, sv = g['directions']['A'], g['directions']['B'], g.get('second') or {}
        for lab, k, spec, th in (('ρ on odor-odor distances, APL on', 'd_rho_on', '+.4f', f'≥ {GAIN["on"]} (or the other ρ)'),
                                 ('ρ on odor-odor distances, APL off', 'd_rho_off', '+.4f', f'≥ {GAIN["off"]}'),
                                 ('wirings gaining ≥ 0.010, APL on', 'n_wirings_gain_on', '.0f', f'≥ {GAIN_WIRINGS} of 4'),
                                 ('wirings gaining ≥ 0.015, APL off', 'n_wirings_gain_off', '.0f', f'≥ {GAIN_WIRINGS} of 4'),
                                 ('pattern part', 'd_pattern_cal', '+.4f', f'≤ +{EQUAL_PATTERN}'),
                                 ('amplitude part', 'd_amplitude_cal', '+.4f', 'beside'),
                                 ('shape part', 'd_shape_cal', '+.4f', 'beside'),
                                 ('spread ratio', 'd_sd_ratio_cal', '+.3f', 'beside'),
                                 ('Pearson r on the per-pair APL effect', 'd_pearson_effect_cal', '+.3f', 'beside'),
                                 ('ρ on the per-pair APL effect', 'd_rho_effect_cal', '+.3f', 'beside'),
                                 ('raw-count loss (unseen half)', 'd_loss_raw', '+.4f', f'≤ +{EQUAL_RAW}'),
                                 ('raw-count loss (all recordings)', 'd_loss_raw_all', '+.4f', 'beside'),
                                 ('mean APL effect', 'd_mean_effect_cal', '+.4f', f'≥ -{EQUAL_MEAN}'),
                                 ('response rate ρ, APL on', 'd_rate_rho_on', '+.3f', f'≥ -{RATE_RHO_SLACK}'),
                                 ('response rate ρ, APL off', 'd_rate_rho_off', '+.3f', f'≥ -{RATE_RHO_SLACK}'),
                                 ('runs outside the bands', 'n_flagged', '.0f', '0'),
                                 ('the block', 'd_block', '+.3f', f'≥ -{BLOCK_SLACK} (pending)'),
                                 ('natural-source alignment, APL on', 'd_vcf_on', '+.3f', f'≥ -{STRUCTURE_SLACK} (pending)'),
                                 ('natural-source alignment, APL off', 'd_vcf_off', '+.3f', f'≥ -{STRUCTURE_SLACK} (pending)'),
                                 ('natural-source alignment, APL on minus off', 'd_vcf_on_minus_off', '+.3f', 'beside'),
                                 ('structure vs wild-type somata', 'd_struct_wt_on', '+.3f', f'≥ -{STRUCTURE_SLACK} (pending)'),
                                 ('distance of paa from B-myr, APL on', 'd_paa_bmyr_on', '+.3f', 'beside')):
            L.append(f'| {lab} | {_f(dA.get(k), spec)} | {_f(dB.get(k), spec)} | {_f(sv.get(k), spec)} | {th} |')
        L += ['', '| ρ | candidate (mean gain) | permuted | inverted | control search | controls hold |', '|---|---:|---:|---:|---:|---|']
        for c, v in g['controls'].items():
            L.append(f'| APL {c} | {_f(v["candidate"], "+.4f")} | {_f(v["permuted"], "+.4f")} | {_f(v["inverted"], "+.4f")} | '
                     f'{_f(v["control search"], "+.4f")} | {"yes" if v["ok"] else "no"} |')
        L += ['', f'Class: **{verdict.get("cls")}** (without the pending items: {verdict.get("class_without_pending")}; outcome '
              f'({verdict.get("outcome")})). Fails: {json.dumps(verdict.get("fails"), default=str)}; pending: '
              f'{json.dumps(verdict.get("pending_fails"), default=str)}.', '']
    p = out / 'step6_entries.csv'
    if p.exists():
        e = pd.read_csv(p)
        m = e[e.meets]
        L += ['## Step 6: the entries five of the six searches agree on (mean log2 factor beyond ±0.5)', '',
              f'{len(m)} entries meet the rule ({int(e.meets_read.sum())} as read at 0.6).' + ('' if len(m) else ''), '',
              '| odor | glomerulus | lab input | mean log2 factor | as read at 0.6 | flies\' SD (log2) | as read, in flies\' SD | Hallem log2 vs lab |',
              '|---|---|---:|---:|---:|---:|---:|---:|']
        for r in m.itertuples():
            L.append(f'| {r.odor} | {r.glomerulus} | {_f(r.lab_input, ".1f")} | {_f(r.mean_log2, "+.2f")} | {_f(r.mean_log2_read, "+.2f")} | '
                     f'{_f(r.flies_sd_log2, ".2f")} | {_f(r.read_in_flies_sd, "+.2f")} | {_f(r.hallem_log2_vs_lab, "+.2f")} |')
        L.append('')
    p = out / 'step6_e108_bmyr.csv'
    if p.exists():
        L += ["## Beside step 6: E108's reading of B-myr within its flies' spread (E108's bmyr.csv: each validation2 fly's "
              "B-myr in the family, the two ρ's of the eight wirings)", '', _md_table(pd.read_csv(p)), '']
    for name, title in (('step7_paa_bmyr.csv', 'Step 7: the distance of paa from B-myr'),
                        ('step7_remy.csv', "Step 7: Remy's validation2 somata (APL on, as fitted)")):
        p = out / name
        if p.exists():
            L += [f'## {title}', '', _md_table(pd.read_csv(p)), '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def score() -> dict:
    """Part score (no simulations; the module docstring)."""
    runs = runs_table()
    runs.to_csv(out_dir() / 'runs.csv', index=False)
    rd = readings(runs)
    searches_table()
    gt, verdict = gate_tables(rd, runs)
    gates = _read_json(out_dir() / 'gate.json').get('gates', {})
    if all(search_finished(n) for n in MAIN):
        recurring_entries()
    if any(load_run('V', f'V_{n}_{T_WIRINGS[0]}') is not None for n in MAIN):
        step7_tables(rd)
    comparison_table(rd, verdict, gates)
    _write_json(out_dir() / 'provenance_score.json', dict(part='score', script=SCRIPT, generated=_now(), **H.code_commits()))
    print(json.dumps(verdict, indent=1, default=str), flush=True)
    return verdict


def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): the held readings (R) and the searches' final points on S (S)."""
    out = out_dir()
    rows, parts = [], []
    rr = E83._read_rows(out / 'runs_R.csv')
    for r in rr.to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=r['tag'], part='R', change=r['change'], search=r['search'], kind=r['kind'],
                         size=r['size'], wiring=r['wiring'], base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'], orn_input=r.get('orn_input'),
                         **{k: r[k] for k in E83.PARAM_KEYS}, kc_kc_inh=r['wKCKC_scale'], thr_factor=r['thr_factor'],
                         apl_factor=r['apl_factor'], rate_on=r['rate_on'], rate_off=r['rate_off'],
                         realism_flag=r.get('realism_flag'), subdir='R',
                         kind_text=(f"C2's family's {r['wiring']} fit with search {r['search']}'s input ({r['kind']}, size "
                                    f"{r['size']}); breadth recomputed for the input; the somatic threshold and the tied APL "
                                    "scale re-held at the fit's APL-off and APL-on counts")))
    for n in SEARCHES:
        p = out / f'search_runs_{n}.csv'
        if not p.exists():
            continue
        for r in pd.read_csv(p).to_dict('records'):
            rows.append(dict(tag=r['tag'], variant=r['tag'], part='S', change=f'{n}_search', search=n, kind='search_final',
                             size=1.0, wiring=r['wiring'], base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                             ref_experiment=r['fit_folder'], ref_fit=r['fit_tag'],
                             orn_input=f'results/{EXPERIMENT}/inputs/searched_{n}.csv',
                             **{k: r[k] for k in E83.PARAM_KEYS}, kc_kc_inh=r['wKCKC_scale'], thr_factor=1.0, apl_factor=1.0,
                             rate_on=r['rate_on'], rate_off=r['rate_off'], realism_flag=r.get('realism_flag'), subdir='S',
                             kind_text=(f"C2's family's {r['wiring']} fit with search {n}'s final input, at the fit's "
                                        "parameters with the lab input's breadth held (the search's own reading)")))
    for r in rows:
        o = load_run(r['subdir'], r['tag'])
        if o is None:
            continue
        parts.append(pd.concat({(r['tag'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['tag'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'search', 'read', 'step7', 'score', 'step6', 'variants', 'e108'))
    ap.add_argument('--search', nargs='+', default=None, choices=list(SEARCHES),
                    help='part search: the searches to run, in turn; part read: the searches whose readings to run')
    ap.add_argument('--threads', type=int, default=WORKER_THREADS, help="part search / 0: each worker's OpenMP threads")
    ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--e108-class', default=None, help="E108's class as written in its note or Status (when classes.csv "
                    'still marks it provisional)')
    ap.add_argument('--e108-class-source', default=None, help='where that class is written')
    ap.add_argument('--max-points', type=int, default=None)
    ap.add_argument('--generations', type=int, default=GENS, help='part search: generations per odor (smoke tests)')
    ap.add_argument('--pop', type=int, default=POP, help='part search: offspring per generation (smoke tests)')
    a = ap.parse_args()
    if a.part == '0':
        part0(a.threads)
    elif a.part == 'search':
        if not a.search:
            ap.error('--part search needs --search')
        confirmed = dict(cls=a.e108_class, source=a.e108_class_source) if a.e108_class else None
        for n in a.search:
            search(n, resume=not a.no_resume, threads=a.threads, confirmed=confirmed, gens=a.generations, pop=a.pop)
    elif a.part == 'read':
        part_read(a.search, a.max_points)
    elif a.part == 'step7':
        part_step7(a.search or MAIN)
    elif a.part == 'score':
        score()
    elif a.part == 'step6':
        recurring_entries()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'e108':
        print(json.dumps(e108_class(), indent=1, default=str))


if __name__ == '__main__':
    main()
