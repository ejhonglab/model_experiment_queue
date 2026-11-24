#!/usr/bin/env python3
"""E134, the models as populations at every recorded stage: the runner of step 2's saved runs and step 3's short runs
at fixed parameters, and their score.

STUB (phase 1, 2026-11-24): the interface only; bodies raise NotImplementedError except trivial helpers (the points,
the claims, the lock guard, the lanes' command lines).

The entry is item E134 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04), design steps 5 and 8. "The family" is C2's
breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `members`, `saved_fit`,
`build_simulator`); S is real, P1, P2, M1 and T is P3, P4, M2, M3. The code is split over three hosts, one owner per
file: viki's data_prep/stage_population_data.py (the flies' populations, the lock) and experiments/E134_flies.py (step
1); hal's kcq/population_fit.py (the views, one-fly units, the population error, the energy guard) and
experiments/E134_archive_score.py (step 2's archive at the KC somata, design steps 6 and 7); atlas's
kcq/stage_readings.py (the model read as each stage is recorded) and this runner. Nothing is fitted or searched, and no
model is changed.

The guard: no model is read until viki's lock (experiments/E134_lock.json) exists, is committed with no uncommitted
change, and passes viki's own check (`require_lock`); every part but 'commands' calls it first.

Parts (CLI `--part`):
- step2 (no simulation; design step 5): the saved runs read as populations (kcq/stage_readings.py) and scored against
  the flies (kcq/population_fit.py), in the stage's reading and its alternatives (`STEP2_SOURCES`):
  - E90's part B: the family's eight wirings on the APL16 input, APL on and off, with every stage saved: the PN
    dendrites, the PN boutons and the KC claws on the panel's 12 megamat odors against the megamat flies, the somata
    on the APL16 panel against the control flies and APL output blocked (TNT) and on the 12 megamat odors against
    Remy's megamat somata; the ORN stage is the lab's input on the 12 odors;
  - E107's part B: the family's somata on the megamat (17 odors) and validation2 (22) inputs, as fitted and re-held at
    6.8 %, against Remy's megamat and validation2 somata;
  - E90's D1 (the four firing-ceiling fits) and B5 (the hemibrain wiring at the flies' 2 s pulse), at the stages they
    saved; other archives with E90's stage layout where present (phase 2 lists them).
  Each wiring is read as one fly: every wiring's error, their mean and the error of the wirings' mean views (hal's
  `wirings_as_flies`), and the measure gate's second check, one wiring against the other seven as flies (1.2 or less).
- step3 (simulations, resumable, lanes; design step 8): about 50 short runs at fixed parameters (`planned_step3`, 47):
  - the family's eight wirings on the megamat input (17 odors) with every stage saved, as fitted (each must reproduce
    E107's saved megamat_fitted_<wiring> counts exactly, and the saved fit on the odors the panels share) and with the
    response rate held at the megamat flies' (0.0915 in Remy's four megamat flies by her criterion; the APL-off rate at
    the wiring's own APL-off / APL-on ratio on its fit panel, `HELD_OFF_RULE`): 16 runs;
  - the same on the kiwi / control mixture input (pyolfsysm's `natmix_orn_deltas`, the measured mixtures used as
    measured), held at the mixture flies' response rate (viki's lock): 16 runs;
  - beside them, as fitted, on the megamat, validation2 and mixture inputs (`BESIDE_MODELS`): the lab's published
    hemibrain model, a uniform draw, the bouton baseline, the current reference model and the firing-ceiling fit: 15
    runs.
  A run with the response rate held is held first (E77's `hold_rates`, the threshold bisected with APL off, then the
  tied APL scales with APL on, breadth held at the lab input's), then run once more at the held parameters with every
  stage saved, which must reproduce the hold's outputs exactly. Each point is claimed by mkdir of
  results/E134_population_baseline/claims/<tag>/ (atomic on NFS; `--reclaim` takes over a dead process's), so any
  number of lanes on any hosts share the list (`--part commands`). Rows in runs_3_<host>_<lane>.csv, provenance in
  provenance_3_<host>_<lane>.json.
- score (no simulation): step 3's runs read and scored as step 2's: the stage table (each megamat stage's population
  error and views, the family end to end from the ORN input, as fitted and held, the beside models), the step table
  (the model's change of each view from stage to stage beside the flies' change and its interval, viki's step 1), the
  no-step yardstick beside, the held-out panels (validation2 and the mixtures at the KC somata against the APL16
  panel's value), the classes the entry fixes beforehand (`classify_stages`, `classify_steps`, `classify_panels`), and
  realism (`realism_rows`).
- commands: the lanes' command lines.

The classes (the entry's, fixed beforehand): a stage is "inside" if the model is no further from the flies' mean than
the most atypical fly is from the others on 6 of the 8 wirings, "far" if above 2 on 6 of 8, "not readable" if the class
differs between the stage's two readings; a step is "explained" if the model's change of every view lies inside the
interval of the flies' change, "no better than no step" if the model's downstream error is at or above the no-step
yardstick, otherwise "in part"; a held-out panel "holds" if the family's KC population error on it is within 0.2 of its
value on the APL16 panel, otherwise "panel-specific".

Realism (reported, not optimized): for every run on the APL16 panel the four bands from the saved tables; for the
megamat, validation2 and mixture runs the overall response rate beside the flies' (0.0915 in Remy's four megamat flies;
the control flies' 0.068), the largest spike counts per KC and odor, and the threshold after the re-hold beside the
fit's own.

Concentrations: the megamat recordings and the megamat ORN input are both at 10^-3; the validation2 ORN input and KC
recordings are one panel at the same concentrations; the mixtures are used as measured. On the APL16 panel the input of
paa, euc and -aPine is the ORN data within 1 log10 unit of the flies' concentrations, and every APL16 number involves
them.

Holds (CLAUDE.md): only Pratyush's hemibrain (the family's and the beside models' wirings); no olfsysm change; no new
connectome table; only Yang's processed products; Remy's processed data without the bouton fly 2022-10-21 fly02 (viki's
lock excludes it). Nothing is published.

Outputs:
- outputs/E134_population_baseline/ (local, large): step3/ (the runs: <tag>_<on|acute_off>_<pn|bouton|claw>_<win|pre|
  pre_settled>.parquet, <tag>_<on|acute_off>_apl.csv, _responses.parquet, _spike_counts.parquet, <tag>_units.json,
  <tag>_orn.parquet), populations/<part>/ (the populations read, one parquet per wiring, stage, panel, condition and
  reading); archived to data/model_outputs/E134_population_baseline/ (step3/ as the runs are made, populations/ when
  scored).
- <shared>/results/E134_population_baseline/ (compact; shared with viki's flies_*.csv and lock.json and hal's
  archive_*.csv): step2_scores_<reading>.csv, step2_wirings_<reading>.csv, step2_gate_wirings_<reading>.csv,
  runs_3_<host>_<lane>.csv, provenance_3_<host>_<lane>.json, claims/, score_stage_table.csv, score_step_table.csv,
  score_panels.csv, score_classes.json, score_realism.csv, score_table.md, curves_<part>.json, provenance_<part>.json.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_population_baseline.py --part step2
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E134_population_baseline.py --part step3 --lane <host>[-<n>] [--reverse]
    ./venv/bin/python -u experiments/E134_population_baseline.py --part score
    ./venv/bin/python -u experiments/E134_population_baseline.py --part commands
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import subprocess
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import population_fit as P
from kcq import stage_readings as R

EXPERIMENT = 'E134_population_baseline'
SCRIPT = 'experiments/E134_population_baseline.py'
# viki's committed copy of the lock (data_prep/stage_population_data.py LOCK_COPY)
LOCK_COPY = H.HERE / 'experiments' / 'E134_lock.json'

# the family (E83's S_WIRINGS, T_WIRINGS; phase 2 checks they are E83's)
FAMILY = 'C2'
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS

# step 2: the saved runs (source, experiment folder, part folder, tags, the input's panel, the stages saved)
Source = namedtuple('Source', 'source folder part tags input_panel stages')
E90_D1_MODELS: Tuple[str, ...] = ('fit1', 'P1_fc', 'P2_fc', 'M1_fc')
STEP2_SOURCES: Tuple[Source, ...] = (
    Source('E90:B', 'E90_model_stages', 'B', WIRINGS, 'apl16', R.STAGES),
    Source('E107:B', 'E107_validation2_panel', 'B',
           tuple(f'{p}_{r}_{w}' for p in ('megamat', 'validation2') for r in ('fitted', 'held') for w in WIRINGS),
           None, ('orn', 'kc_soma')),
    Source('E90:D1', 'E90_model_stages', 'D1', E90_D1_MODELS, 'apl16', R.STAGES),
    Source('E90:B5', 'E90_model_stages', 'B5', ('real_pulse2',), 'apl16', R.STAGES),
)
# the measure gate's second check: one wiring against the other seven as flies
GATE_WIRING_MAX = 1.2

# step 3: the points
Point = namedtuple('Point', 'panel reading model')
STEP3_FAMILY_PANELS: Tuple[str, ...] = ('megamat', 'mixtures')
STEP3_READINGS: Tuple[str, ...] = ('fitted', 'held')
BESIDE_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures')
# PROVISIONAL until phase 2 checks each against its saved outputs (where they exist) and viki's lock records it
BESIDE_MODELS: Dict[str, str] = {
    'published_hemibrain': "the lab's published hemibrain model: fit_mb_model with pn2kc_connections='hemibrain' and "
                           "pyolfsysm.util.hemibrain_paper_repro_kws, tuned to 0.0915 on the paper's megamat input "
                           "(paper_megamat_orn_deltas); its tuned thresholds and APL scale carried to the other inputs",
    'uniform_draw': "the published uniform model, one draw: pn2kc_connections='uniform', n_claws=7, seed 0, with "
                    "pyolfsysm.util.paper_repro_kws, tuned as the published hemibrain model",
    'bouton_baseline': "the bouton baseline: pyolfsysm.baselines' 'bouton_base' (E13's bouton base, tuned with APL on "
                       "to 5 % on the APL16 input), its parameters carried",
    'reference': "the current reference model: E44's operating point (E69's 'ref40': eight glomeruli at 1.5, 40 Hz, "
                 "6.8 %), its parameters carried",
    'firing_ceiling': "the firing-ceiling fit: E69's fit 1 (E82_protocol_runs' model_parts('fit1')), its parameters, "
                      "refractory period and adaptation carried",
}
# the held response rates: the megamat flies' (Remy's four megamat flies by her criterion; CLAUDE.md, "Realism"), the
# mixture flies' from viki's lock (None here)
HELD_RATE: Dict[str, Optional[float]] = {'megamat': 0.0915, 'mixtures': None}
# the held APL-off rate: the held APL-on rate times the wiring's own APL-off / APL-on ratio on its fit panel ('ratio'),
# or E107's rule, the wiring's own APL-off fraction ('e107')
HELD_OFF_RULES: Tuple[str, ...] = ('ratio', 'e107')
HELD_OFF_RULE = 'ratio'
# the classes (the entry's)
INSIDE_MIN_WIRINGS = 6
FAR_ERROR = 2.0
PANEL_HOLDS_WITHIN = 0.2
# the lanes: threads per lane by host, and the lanes of each host (a lane's name, reversed or not)
LANES: Dict[str, Tuple[int, Tuple[Tuple[str, bool], ...]]] = {
    'viki': (16, (('viki', False),)),
    'hal': (16, (('hal-1', False), ('hal-2', True))),
    'atlas': (8, (('atlas', True),)),
}


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E134_population_baseline/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    """outputs/E134_population_baseline/[<part>/] (created)."""
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E134_population_baseline/."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


# ---------------------------------------------------------------------------
# the guard

def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git and has no uncommitted change."""
    try:
        r = subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def require_lock() -> dict:
    """The committed lock, before any model is read: experiments/E134_lock.json exists, is tracked with no uncommitted
    change, and passes viki's check (`stage_population_data.require_lock`). Raises PermissionError otherwise.

    Phase 2 adds: the lock's readings of the model are `stage_readings.reading_definitions()`.
    """
    if not LOCK_COPY.exists():
        raise PermissionError(f'{LOCK_COPY.name}: no lock yet (viki writes it, data_prep/stage_population_data.py '
                              '--lock); no model is read before it is committed')
    if not _committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is not committed (or has uncommitted changes)')
    from data_prep import stage_population_data as D
    try:
        return D.require_lock()
    except NotImplementedError:
        raise PermissionError("viki's lock check (stage_population_data.require_lock) is not implemented yet") from None
    except RuntimeError as e:
        raise PermissionError(f'the lock check failed: {e}') from None


# ---------------------------------------------------------------------------
# step 2: the saved runs

def step2_runs(source: Source) -> List[R.SavedRun]:
    """A source's saved runs (`stage_readings.load_e90_run`, `load_e107_run`), each with its files' sha256; a run whose
    files are missing is left out and reported."""
    raise NotImplementedError


def score_family(readings: Mapping[str, R.ModelReading], reading_name: str) -> pd.DataFrame:
    """The family's wirings as flies, per (stage, panel, condition): each wiring's `P.Score`, their mean and the score
    of the wirings' mean views (hal's `wirings_as_flies`, against viki's flies of `stage_readings.cohort_of`), with the
    energy guard beside (hal's `energy_guard`); rows of model, wiring ('mean', 'common' for the summaries), stage,
    panel, condition, reading and the Score columns."""
    raise NotImplementedError


def gate_wirings(readings: Mapping[str, R.ModelReading], reading_name: str) -> pd.DataFrame:
    """The measure gate's second check: each wiring's population scored against the other seven wirings taken as the
    flies (hal's `fly_reference` on their populations), per (stage, panel, condition); `passes` where every wiring is at
    GATE_WIRING_MAX or less."""
    raise NotImplementedError


def part_step2(readings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 2: every source of STEP2_SOURCES read in each stage's default reading and its alternatives (or `readings`,
    by name), scored (`score_family`, `gate_wirings`); step2_scores_<reading>.csv, step2_wirings_<reading>.csv,
    step2_gate_wirings_<reading>.csv, curves_step2.json, provenance_step2.json; the populations under
    outputs/.../populations/step2/."""
    require_lock()
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 3: the points

def tag_of(pt: Point) -> str:
    """'<panel>_<reading>_<model>' (E107's `run_tag` form), e.g. 'megamat_fitted_real', 'validation2_fitted_reference'."""
    if pt.panel not in R.PANELS or pt.reading not in STEP3_READINGS:
        raise ValueError(f'{pt}')
    return f'{pt.panel}_{pt.reading}_{pt.model}'


def parse_tag(tag: str) -> Point:
    """The point of a tag (`tag_of`'s inverse)."""
    panel, reading, model = tag.split('_', 2)
    pt = Point(panel, reading, model)
    tag_of(pt)
    return pt


def planned_step3(models: Optional[Sequence[str]] = None, panels: Optional[Sequence[str]] = None,
                  readings: Optional[Sequence[str]] = None, reverse: bool = False) -> List[Point]:
    """Step 3's points, in the order a forward lane takes them: the family on the megamat input (as fitted, then held,
    wiring by wiring in WIRINGS order), then on the mixture input, then the beside models (BESIDE_MODELS order) as
    fitted on BESIDE_PANELS; restricted to `models` (wiring labels and beside models' names), `panels` and `readings`
    where given; `reverse` the same list from its other end."""
    pts = [Point(p, r, w) for p in STEP3_FAMILY_PANELS for r in STEP3_READINGS for w in WIRINGS]
    pts += [Point(p, 'fitted', m) for m in BESIDE_MODELS for p in BESIDE_PANELS]
    pts = [pt for pt in pts if (models is None or pt.model in models) and (panels is None or pt.panel in panels)
           and (readings is None or pt.reading in readings)]
    return pts[::-1] if reverse else pts


def claim(tag: str, lane: str, reclaim: bool = False, root: Optional[Path] = None) -> bool:
    """Claims a point by mkdir of results/E134_population_baseline/claims/<tag>/ (atomic on NFS); False if another
    process holds it (unless `reclaim`)."""
    d = (out_dir() if root is None else root) / 'claims' / tag
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(tag=tag, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(tag: str, root: Optional[Path] = None) -> None:
    """Removes a point's claim (after a failed run, so another lane may take it)."""
    shutil.rmtree((out_dir() if root is None else root) / 'claims' / tag, ignore_errors=True)


def panel_input(panel: str) -> pd.DataFrame:
    """The model's ORN input of a panel on the lab input's 36 glomeruli (H.apl16_orn_deltas' rows): 'megamat' and
    'validation2' E107's `panel_input`; 'mixtures' pyolfsysm's `natmix_orn_deltas()` (diagnostics dropped, the kiwi and
    control panels' odors and mixtures as measured); 'apl16' the lab's input. Raises if a glomerulus is missing or NaN."""
    raise NotImplementedError


def held_targets(pt: Point, size: int) -> Tuple[int, int, float, float]:
    """(target_off, target_on, fraction_off, fraction_on) of a held point as response counts on its KCs x odors
    (`size`): APL on at HELD_RATE (the mixtures' from viki's lock), APL off by HELD_OFF_RULE."""
    raise NotImplementedError


def build_point(pt: Point) -> Tuple:
    """A point's fixed-parameter simulator on its panel's input, APL's activity saved (`mp.kc.save_inh_sims`):
    the family's wiring through E83's `build_simulator(label, orn=panel_input(panel), breadth=<the lab input's>)`; a
    beside model through its own builder (BESIDE_MODELS). Returns (simulate, info) with info's 'mp', 'rv',
    'odor_order', 'W' (the claw x bouton wiring), 'params', 'params_off'."""
    raise NotImplementedError


def run_stages(pt: Point, sim, info: dict, params: Mapping[str, float], params_off: Mapping[str, float]) -> Dict[str, dict]:
    """Both conditions of a point at the given parameters, each read as E90's `read_stages` over MODEL_WINDOW with the
    pre-odor means over PRE_WINDOW and SETTLED_PRE_WINDOW ('<stage>_pre_settled'), with the responses and spike
    counts: {'on', 'off'}: the stage dict."""
    raise NotImplementedError


def check_identity(pt: Point, out: Mapping[str, dict], hold: Optional[dict] = None) -> dict:
    """A point's identity before its stages are saved: as fitted, the counts and responses equal E107's saved run of
    the same tag where it exists (megamat), and the saved fit on the odors the panels share (E107's `run_identity`);
    held, they equal the hold's final outputs. 0 differing counts and responses, or RuntimeError."""
    raise NotImplementedError


def save_run(tag: str, out: Mapping[str, dict], meta: Mapping[str, np.ndarray], orn: pd.DataFrame,
             archive: bool = True) -> None:
    """A point's stage matrices, responses, spike counts, APL means, units' metadata (E90's units.json plus each claw's
    bouton) and input, locally (outputs/.../step3/) and in the model-output archive, written atomically."""
    raise NotImplementedError


def load_run(tag: str) -> Optional[R.SavedRun]:
    """A saved step-3 point (locally, else the archive) as `stage_readings.SavedRun`; None if absent."""
    raise NotImplementedError


def run_point(pt: Point, lane: str, table: Path, prov_path: Path, reclaim: bool = False) -> Optional[dict]:
    """One point, claimed (None if done or claimed elsewhere): built, held if its reading is 'held', run with its
    stages, checked (`check_identity`), saved (`save_run`), its row appended to `table` (absolute parameters, factors
    and targets, response counts and rates, largest counts, seconds, host, lane, threads, the input's sha256, the
    commits) and its provenance to `prov_path`; the claim released on failure."""
    raise NotImplementedError


def part_step3(models: Optional[Sequence[str]] = None, panels: Optional[Sequence[str]] = None,
               readings: Optional[Sequence[str]] = None, lane: Optional[str] = None, reverse: bool = False,
               reclaim: bool = False) -> pd.DataFrame:
    """Step 3, resumable, one point claimed at a time (`planned_step3`, `run_point`); rows in
    runs_3_<host>_<lane>.csv."""
    require_lock()
    raise NotImplementedError


def runs_table() -> pd.DataFrame:
    """Every row of every host's and lane's runs_3_*.csv, the last row per tag."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# score

def stage_table(scores: pd.DataFrame) -> pd.DataFrame:
    """The stage table: per megamat stage, the flies' floor and one recording against the others (viki's step 1) and
    the model's population error with its views, each wiring, their mean and the common error, the family as fitted
    and held and the beside models (score_stage_table.csv)."""
    raise NotImplementedError


def step_table(readings: Mapping[str, R.ModelReading]) -> pd.DataFrame:
    """The step table: per step between neighbouring megamat stages and per part, the model's change (hal's
    `step_change` on each wiring's views and on their mean) beside the flies' change and its interval and the no-step
    yardstick (viki's flies_steps.csv, flies_no_step.csv) (score_step_table.csv)."""
    raise NotImplementedError


def classify_stages(stage_tab: pd.DataFrame, by_reading: Mapping[str, pd.DataFrame]) -> Dict[str, str]:
    """Each stage's class: 'inside' (no further from the flies' mean than the most atypical fly on INSIDE_MIN_WIRINGS
    of 8 wirings), 'far' (above FAR_ERROR on 6 of 8), 'not readable' (the class differs between the stage's two
    readings), otherwise 'between'."""
    raise NotImplementedError


def classify_steps(step_tab: pd.DataFrame) -> Dict[str, str]:
    """Each step's class: 'explained', 'no better than no step' or 'in part' (the module docstring)."""
    raise NotImplementedError


def classify_panels(scores: pd.DataFrame) -> Dict[str, str]:
    """Each held-out panel's class at the KC somata: 'holds' if the family's population error there is within
    PANEL_HOLDS_WITHIN of its APL16 value (APL on, against the control flies), else 'panel-specific'."""
    raise NotImplementedError


def realism_rows(run: R.SavedRun, fit_params: Optional[Mapping[str, float]] = None) -> dict:
    """A run's realism, reported only: overall response rates with APL on and off and their ratio, the largest spike
    count per KC and odor in each condition, the threshold beside the fit's own (held runs), and on the APL16 panel the
    four bands (kcq/apl_effect.py's)."""
    raise NotImplementedError


def part_score() -> dict:
    """Step 3's runs (and E107's validation2 runs of the family) read and scored as step 2's; the stage table, the step
    table, the held-out panels, the classes, realism, score_table.md and provenance_score.json."""
    require_lock()
    raise NotImplementedError


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the script, `H.model_provenance` (or `H.code_commits()` for parts without simulation),
    the lock's sha256, the readings (`stage_readings.reading_definitions`), host, time and `extra`."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# lanes and CLI

def print_commands() -> None:
    """The lanes' command lines (LANES): every lane runs the whole list and skips what another has claimed; the
    reversed lanes take it from the other end."""
    pre = 'OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    for host, (t, lanes) in LANES.items():
        for name, rev in lanes:
            print(f"{host}: {pre.format(t=t)} --part step3 --lane {name}{' --reverse' if rev else ''}")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step2', 'step3', 'score', 'commands'))
    ap.add_argument('--models', nargs='+', default=None,
                    help='step3: the wirings and beside models of the lane (default all)')
    ap.add_argument('--panels', nargs='+', default=None, choices=list(R.PANELS), help='step3: the panels (default all)')
    ap.add_argument('--readings', nargs='+', default=None,
                    help="step3: 'fitted' and / or 'held'; step2: the stage readings by name (default each stage's "
                         "default and alternatives)")
    ap.add_argument('--lane', default=None, help="the lane's name (default the host's)")
    ap.add_argument('--reverse', action='store_true', help='step3: take the points from the other end')
    ap.add_argument('--reclaim', action='store_true', help="step3: take over the points' claims (a process that died)")
    a = ap.parse_args(argv)
    if a.part == 'step2':
        part_step2(a.readings)
    elif a.part == 'step3':
        part_step3(a.models, a.panels, a.readings, a.lane or _host(), a.reverse, a.reclaim)
    elif a.part == 'score':
        part_score()
    elif a.part == 'commands':
        print_commands()


if __name__ == '__main__':
    main()
