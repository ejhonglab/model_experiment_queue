#!/usr/bin/env python3
"""E134, the models as populations at every recorded stage: the runner of step 2's saved runs and step 3's short runs
at fixed parameters, and their score.

The entry is item E134 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04), design steps 5 and 8. "The family" is C2's
breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `members`, `saved_fit`,
`build_simulator`); S is real, P1, P2, M1 and T is P3, P4, M2, M3. The code is split over three hosts, one owner per
file: viki's data_prep/stage_population_data.py (the flies' populations, the lock) and experiments/E134_flies.py (step
1); hal's kcq/population_fit.py (the views, one-fly units, the population error, the energy guard) and
experiments/E134_archive_score.py (step 2's archive at the KC somata, design steps 6 and 7); atlas's
kcq/stage_readings.py (the model read as each stage is recorded) and this runner. Nothing is fitted or searched, and no
model is changed.

The guard: no model is read against the flies until viki's lock (experiments/E134_lock.json) exists, is committed with
no uncommitted change, and passes viki's own check (`stage_readings.require_lock`); 'step2', 'step3' and 'score' call
it first.

Parts (CLI `--part`):
- step2 (no simulation; design step 5): the saved runs read as populations (kcq/stage_readings.py) and scored against
  the flies (kcq/population_fit.py), in each stage's reading and its alternatives (`STEP2_SOURCES`):
  - E90's part B: the family's eight wirings on the APL16 input, APL on and off, with every stage saved: the PN
    dendrites, the PN boutons and the KC claws on the panel's 12 megamat odors against the megamat flies, the somata
    on the APL16 panel against the control flies and APL output blocked (TNT) and on the 12 megamat odors against
    Remy's megamat somata; the ORN stage is the lab's input on the 12 odors;
  - E107's part B: the family's somata on the megamat (17 odors) and validation2 (22) inputs, as fitted and re-held at
    6.8 %, against Remy's megamat and validation2 somata;
  - E90's D1 (the four firing-ceiling fits) and B5 (the hemibrain wiring at the flies' 2 s pulse), at the stages they
    saved.
  Each wiring is read as one fly: every wiring's error, their mean and the error of the wirings' mean views (hal's
  `wirings_as_flies`), beside the flies' floor and their most atypical fly (hal's `floor`, `leave_one_out`), and the
  measure gate's second check, one wiring against the other wirings as flies (`gate_wirings`, GATE_WIRING_MAX).
- step3 (simulations, resumable, lanes; design step 8): 47 short runs at fixed parameters (`planned_step3`):
  - the family's eight wirings on the megamat input (17 odors) with every stage saved, as fitted (each must reproduce
    E107's saved megamat_fitted_<wiring> counts exactly) and with the response rate held at the megamat flies' (0.0915
    in Remy's four megamat flies by her criterion, the lock's 'remy_megamat_criterion_rate', the operating point
    CLAUDE.md's main goal names, by two spikes: HELD_RATE_SOURCE; the APL-off rate at the wiring's own APL-off /
    APL-on ratio on its fit panel, `HELD_OFF_RULE`): 16 runs;
  - the same on the kiwi / control mixture input ('mixtures': pyolfsysm's `natmix_orn_deltas` on the lab input's 36
    glomeruli and on the 16 stimuli of the two mixture panels, `mixture_input`; the measured mixtures used as
    measured), held at the mixture flies' response rate (the lock's KC somata of 'kmix' and 'cmix', `held_rate`), read
    on both panels: 16 runs. The mixtures share no odor at the same concentration with the fit's panel, so an as-fitted
    mixture run has no saved counts to reproduce (its identity row says so);
  - beside them, as fitted, on the megamat, validation2 and mixture inputs (`BESIDE_MODELS`, `build_beside`): the
    lab's published hemibrain model, a uniform draw, the bouton baseline, the current reference model and the
    firing-ceiling fit: 15 runs. Before a beside model's first point on a host, its build on its own fit panel must
    reproduce its saved outputs exactly (`replay_beside`, cached in results/.../replay_<model>_<host>.json).
  A run with the response rate held is held first (E77's `hold_rates`, the threshold bisected with APL off, then the
  tied APL scales with APL on, breadth held at the lab input's), then built again and run at the held parameters with
  every stage saved, which must reproduce the hold's outputs exactly (`check_identity`). Each point is claimed by mkdir
  of results/E134_population_baseline/claims/<tag>/ (atomic on NFS; `--reclaim` takes over a dead process's), so any
  number of lanes on any hosts share the list (`--part commands`); a point whose row is in a runs table is done. Rows in
  runs_3_<host>_<lane>.csv, provenance in provenance_3_<host>_<lane>.json.
- score (no simulation): step 3's runs read and scored as step 2's: the stage table (each megamat stage's population
  error and views, the family end to end from the ORN input, as fitted and held, the beside models), the step table
  (the model's change of each part from stage to stage beside the flies' change and its interval, and the no-step
  yardstick), the held-out panels (validation2, from E107's saved runs, and the two mixture panels at the KC somata
  against the APL16 panel's value from step 2), the classes the entry fixes beforehand (`classify_stages`, `classify_steps`,
  `classify_panels`), and realism (`realism_rows`).
- commands: the lanes' command lines.

The classes (the entry's, fixed beforehand): a stage is "inside" if the model is no further from the flies' mean than
the most atypical fly is from the others on 6 of the 8 wirings, "far" if above 2 on 6 of 8, "not readable" if the class
differs between the stage's two readings; a step is "explained" if the model's change of every part lies inside the
interval of the flies' change, "no better than no step" if the model's downstream error is at or above the no-step
yardstick, otherwise "in part"; a held-out panel "holds" if the family's KC population error on it is within 0.2 of its
value on the APL16 panel, otherwise "panel-specific".

The beside models (`BESIDE_MODELS`), as far as the code and saved outputs define them; where the entry's name was
ambiguous the closest saved configuration was taken, as each line says:
- 'published_hemibrain': the lab's published hemibrain model as pyolfsysm's paper-reproduction kwargs make it
  (pn2kc_connections='hemibrain', `hemibrain_paper_repro_kws`: weight divisor 20, no KCs dropped), at the tuned
  parameters it saved (pyolfsysm/data/sent_to_remy/2025-03-18/.../params.csv: fixed_thr 268.04, wAPLKC 4.623, wKCAPL
  0.00252, tuned to 0.0915 on the paper's megamat input). Its megamat and validation2 points run on the paper's inputs
  (`paper_megamat_orn_deltas`, `paper_validation2_orn_deltas`) at each panel's published parameters (validation2's own
  saved fit, 2025-03-19, threshold 283.18), and must reproduce the saved spike counts; the mixtures run on
  `mixture_input` (all its glomeruli) at the megamat parameters. It has no boutons or claws: its PN stage and
  somata are saved.
- 'uniform_draw': the paper's uniform model (pn2kc_connections='uniform', 7 claws, `paper_repro_kws`, glomeruli with a
  '+' kept as pyolfsysm's test_uniform_paper_repro keeps them to match the saved outputs), its first draw
  (seed UNIFORM_SEED, the first of the 100 in pyolfsysm/data/sent_to_anoop/2024-05-16), at that draw's saved tuned
  parameters, on the same inputs as the published hemibrain model (the megamat parameters throughout: no validation2
  draw is saved); its megamat responses must equal the saved ones of that seed. PN stage and somata.
- 'bouton_baseline': E3's bouton base at its tuned parameters (E12's `e3_bouton_params`, tuned with APL on to 5 % on
  the APL16 input; the hemibrain wiring, connectome APL weights without breadth weighting, `E10_joint_fit.base_kws`),
  as a fixed-parameter simulator; replayed on the APL16 input against hal's saved E3 run (bouton_pnapl_x1). Every
  stage.
- 'reference': "the current reference model" read as E69's 'ref40' (E44's operating point
  thr40_gain1.5_rate0.068_k2_cap3: eight glomeruli at gain 1.5, a 40 Hz per-claw threshold, 6.8 %), rebuilt as E69's
  `Model('ref40')` rebuilds it (E43's `fit_model`); replayed against its saved outputs. Its claws pass a per-claw
  threshold that E90's claw reading (the linear drive) does not model, so its claw stage is not saved; PN, boutons and
  somata are.
- 'firing_ceiling': E69's fit 1 (E82's `model_parts('fit1')`, its refractory period and adaptation), built as E82 and
  E90's D1 build it; replayed against its saved outputs. Every stage.

Realism (reported, not optimized): for the family on the APL16 panel the four bands from results/apl_effect/
by_variant.csv (`apl16_bands`); for every step-3 run the overall response rates with APL on and off and their ratio, the
largest spike counts per KC and odor, and the threshold after the re-hold beside the fit's own (`realism_rows`), beside
the flies' 0.0915 (Remy's four megamat flies) and the control flies' 0.068.

Concentrations: the megamat recordings and the megamat ORN input are both at 10^-3; the validation2 ORN input and KC
recordings are one panel at the same concentrations; the mixtures are used as measured. On the APL16 panel the input of
paa, euc and -aPine is the ORN data within 1 log10 unit of the flies' concentrations, and every APL16 number involves
them.

Holds (CLAUDE.md): only Pratyush's hemibrain (the family's and the beside models' wirings); no olfsysm change; no new
connectome table; only Yang's processed products; Remy's processed data without the bouton fly 2022-10-21 fly02 (viki's
lock excludes it). Nothing is published.

Outputs:
- outputs/E134_population_baseline/ (local, large): step3/ (the runs, `stage_readings.save_step3_run`'s layout),
  populations/<part>/ (the populations read under each stage's default reading, one parquet per model, wiring, stage,
  panel and condition), beside/ (the beside models' caches); step3/ is also written to the archive
  data/model_outputs/E134_population_baseline/step3/ as the runs are made.
- <shared>/results/E134_population_baseline/ (compact; shared with viki's flies_*.csv and lock.json and hal's
  archive_*.csv): step2_scores.csv, step2_gate_wirings.csv, curves_step2.json, provenance_step2.json;
  runs_3_<host>_<lane>.csv, provenance_3_<host>_<lane>.json, identity_3_<host>_<lane>.json,
  replay_<model>_<host>.json, claims/;
  score_scores.csv, score_stage_table.csv, score_step_table.csv, score_panels.csv, score_classes.json,
  score_realism.csv, score_table.md, provenance_score.json.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E134_population_baseline.py --part step2
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        experiments/E134_population_baseline.py --part step3 --lane <host>[-<n>] [--reverse]
    ./venv/bin/python -u experiments/E134_population_baseline.py --part score
    ./venv/bin/python -u experiments/E134_population_baseline.py --part commands
"""
from __future__ import annotations

import argparse
import ast
import gc
import json
import os
import shutil
import socket
import time
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from kcq import population_fit as P
from kcq import stage_readings as R

EXPERIMENT = 'E134_population_baseline'
SCRIPT = 'experiments/E134_population_baseline.py'
LOCK_COPY = R.LOCK_COPY

# the family (E83's S_WIRINGS, T_WIRINGS; `check_wirings` checks they are E83's)
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
# the measure gate's second check: one wiring against the other wirings as flies
GATE_WIRING_MAX = 1.2
N_FLOOR = 1000

# step 3: the points
Point = namedtuple('Point', 'panel reading model')
STEP3_FAMILY_PANELS: Tuple[str, ...] = ('megamat', 'mixtures')
STEP3_READINGS: Tuple[str, ...] = ('fitted', 'held')
BESIDE_PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures')
BESIDE_MODELS: Dict[str, str] = {
    'published_hemibrain': "the lab's published hemibrain model (pyolfsysm's paper-reproduction kwargs) at its saved "
                           "tuned parameters, on the paper's inputs",
    'uniform_draw': "the paper's uniform model, its first saved draw, at that draw's tuned parameters",
    'bouton_baseline': "E3's bouton base at its tuned parameters (5 % on the APL16 input)",
    'reference': "E69's 'ref40', E44's operating point thr40_gain1.5_rate0.068_k2_cap3, at its fitted parameters",
    'firing_ceiling': "E69's fit 1 (E82's model_parts('fit1')) at its fitted parameters",
}
# the panel the beside model was fitted (tuned) on, for its replay
BESIDE_FIT_PANEL: Dict[str, str] = {'published_hemibrain': 'megamat', 'uniform_draw': 'megamat',
                                    'bouton_baseline': 'apl16', 'reference': 'apl16', 'firing_ceiling': 'apl16'}
UNIFORM_SEED = 94895
REFERENCE_SPEC = 'E44_operating_point:thr40_gain1.5_rate0.068_k2_cap3'
E3_SAVED = ('data', 'hal_outputs', 'E3_drive_dependent_gain', 'bouton_pnapl_x1')
# the held response rates come from viki's lock (`held_rate`), never from constants here: 'remy_criterion' (megamat:
# Remy's four flies by her criterion, 'remy_megamat_criterion_rate', 0.0915, the operating point CLAUDE.md's main goal
# names, held by two spikes) or 'lock_criterion' (the panel's KC somata populations' rates at the lock's criteria: the
# two mixture panels, each population's mean over flies weighted by its number of stimuli; about 0.056 and 0.052 in
# the lock as the planner re-derived it, chat.md 2026-11-24 11:34, where Remy's flies' soma criterion is the level that
# gives 0.0915)
HELD_RATE_SOURCE: Dict[str, str] = {'megamat': 'remy_criterion', 'mixtures': 'lock_criterion'}
# the held APL-off rate: the held APL-on rate times the wiring's own APL-off / APL-on ratio on its fit panel ('ratio'),
# or E107's rule, the wiring's own APL-off fraction ('e107')
HELD_OFF_RULES: Tuple[str, ...] = ('ratio', 'e107')
HELD_OFF_RULE = 'ratio'
PARAM_KEYS = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale')
# the classes (the entry's)
INSIDE_MIN_WIRINGS = 6
FAR_ERROR = 2.0
PANEL_HOLDS_WITHIN = 0.2
MEGAMAT_STAGES = R.STAGES
# the stages of the megamat stage error (CLAUDE.md's main goal), the PN dendrites on four views
MEGAMAT_ERROR_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
# the step class's threshold: the 95th percentile of the resampled flies' distances; the consistency gate's minimum
STEP_NULL_Q = 0.95
GATE_CONSISTENCY_MIN = 0.8
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


# the columns a held point's row carries beyond a fitted point's, in `run_point`'s order (inserted before 'stages'); they
# also let `read_runs_csv` read a file whose header was written by a fitted row (viki's runs_3 file of 2026-11-24 12:26)
HELD_EXTRA_COLUMNS = ('thr_factor', 'apl_factor', 'target_on', 'target_off', 'frac_on_target', 'frac_off_target',
                      'at_bound_thr', 'at_bound_apl', 'hold_n_sims', 'held_off_rule', 'fit_fixed_thr')


def read_runs_csv(path: Path) -> pd.DataFrame:
    """A lane's runs_3 table. A file written before the union fix (`_append_csv`) can hold held rows under a fitted
    row's shorter header: those rows are read with `HELD_EXTRA_COLUMNS` inserted before 'stages'. ValueError for any
    other row length."""
    import csv
    import io
    try:
        return pd.read_csv(path)
    except pd.errors.ParserError:
        pass
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    header, body = rows[0], [r for r in rows[1:] if r]
    i = header.index('stages')
    held = header[:i] + list(HELD_EXTRA_COLUMNS) + header[i:]
    parts = []
    for names in (header, held):
        group = [r for r in body if len(r) == len(names)]
        if group:
            buf = io.StringIO()
            csv.writer(buf).writerows([names] + group)
            buf.seek(0)
            parts.append(pd.read_csv(buf).assign(_line=[body.index(r) for r in group]))
    if sum(len(p) for p in parts) != len(body):
        raise ValueError(f'{path}: rows of {sorted({len(r) for r in body})} fields under a header of {len(header)}')
    return pd.concat(parts, ignore_index=True, sort=False).sort_values('_line').drop(columns='_line').reset_index(drop=True)


def _append_csv(path: Path, row: Mapping) -> None:
    """Appends one row to a lane's table. When the row's columns differ from the file's, the file is rewritten with the
    union of the columns (the file's order, then the row's new ones) and replaced atomically; only the lane's own
    process writes its file."""
    new = pd.DataFrame([dict(row)])
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        new.to_csv(path, index=False)
        return
    old = read_runs_csv(path)
    if list(old.columns) == list(new.columns):
        new.to_csv(path, mode='a', header=False, index=False)
        return
    cols = list(old.columns) + [c for c in new.columns if c not in old.columns]
    df = pd.concat([old, new], ignore_index=True, sort=False)[cols]
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def _differ(a: pd.DataFrame, b: pd.DataFrame) -> int:
    """How many entries of two outputs differ (-1 if their shapes or indices differ)."""
    if a.shape != b.shape or not a.index.equals(b.index):
        return -1
    return int((a.to_numpy(dtype=float) != b.to_numpy(dtype=float)).sum())


def require_lock() -> dict:
    """`stage_readings.require_lock`: the committed lock, before any model is read against the flies; PermissionError
    otherwise."""
    return R.require_lock()


def check_wirings() -> None:
    """The family's wirings are E83's (S, T and their order); RuntimeError otherwise."""
    from experiments import E83_cross_wiring_fit as E83
    if (tuple(E83.S_WIRINGS), tuple(E83.T_WIRINGS)) != (S_WIRINGS, T_WIRINGS):
        raise RuntimeError(f"the family's wirings are not E83's: {E83.S_WIRINGS} / {E83.T_WIRINGS}")


# ---------------------------------------------------------------------------
# scoring (shared by step 2 and score)

def breadth_edges(n_odors: int):
    """The breadth bins of a panel read on `n_odors` odors: viki's `breadth_bins` (the lock's, scaled to the panel's
    size); hal's APL16 edges where viki's is not implemented and the panel has 16 odors."""
    from data_prep import stage_population_data as D
    try:
        return D.breadth_bins(int(n_odors))
    except NotImplementedError:
        if int(n_odors) == 16:
            return P.APL16_BREADTH_EDGES
        raise


def fly_side_read(flies) -> Tuple[List[P.Population], List[str]]:
    """The flies and odors a recorded population is read on, as viki's step 1 reads it (experiments/E134_flies.py
    `reading`, the lock's rule: validation2's KC somata on the largest odor set that at least 6 flies share, and those
    flies; every other population all its flies, on the odors they were all given), from these flies' tables (hal's
    `fly_populations`; so a half's population gives the half's)."""
    from experiments import E134_flies as F
    fp = P.fly_populations(flies)
    if flies.key == ('kc_soma', 'validation2', 'wild_type'):
        od, keep = F.v2_odor_set(fp)
        return [f for f in fp if f.name in keep], od
    return fp, P.common_odors(fp)


def fly_reference_for(flies, odors: Sequence[str]) -> Tuple[List[P.Population], P.FlyReference, List[str]]:
    """The flies' populations and odors as viki's step 1 reads them (`fly_side_read`), their reference on those of
    the odors the model has too, and those odors."""
    fp, have = fly_side_read(flies)
    od = [o for o in odors if o in set(have)]
    return fp, P.fly_reference(fp, od, breadth_edges(len(od))), od


def score_family(readings: Mapping[str, R.ModelReading], reading_name: str, flies: R.FlyCache,
                 n_floor: int = N_FLOOR) -> pd.DataFrame:
    """The family's wirings as flies, per (stage, panel, condition) with flies: each wiring's `P.Score`, their mean and
    the score of the wirings' mean views (hal's `wirings_as_flies`), the energy guard per wiring (hal's `energy_guard`,
    every population's pooled values mapped rank for rank onto the flies' pooled values), and beside them the flies'
    floor (`P.floor`, `n_floor` draws) and one recording against the others (the mean and largest of `leave_one_out`);
    rows of model, wiring ('mean', 'common' for the summaries), stage, panel, condition, reading, n_units and the
    Score columns."""
    keys = sorted({k for r in readings.values() for k in r.populations})
    rows = []
    for stage, panel, cond in keys:
        fl = flies.get((stage, panel, cond))
        if fl is None:
            continue
        pops = [r.populations[(stage, panel, cond)] for r in readings.values() if (stage, panel, cond) in r.populations]
        fp, ref, od = fly_reference_for(fl, list(pops[0].values.columns))
        tab = P.wirings_as_flies(pops, ref)
        refvals = R.pooled_values([f.values for f in fp], od)
        eg = {p.name: P.energy_guard(p, fp, od, reference=refvals) for p in pops}
        loo = P.leave_one_out(fp, ref=ref)
        fl_floor = P.floor_summary(P.floor(ref, n_boot=n_floor))
        for name, row in tab.iterrows():
            model, _, wiring = str(name).partition('/')
            rec = dict(model=(model if wiring else pops[0].name.split('/')[0]), wiring=wiring or str(name),
                       stage=stage, panel=panel, condition=cond, reading=reading_name, n_odors=len(od),
                       n_units=int(next((p.values.shape[0] for p in pops if p.name == name), -1)),
                       flies_floor=float(fl_floor['error']), flies_loo_mean=float(loo['error'].mean()),
                       flies_loo_max=float(loo['error'].max()), n_flies=len(fp))
            g = eg.get(name)
            if g is None and str(name) == 'mean':
                g = {k: float(np.mean([v[k] for v in eg.values()])) for k in ('to_flies', 'fly_to_fly', 'ratio')}
            rec.update({f'energy_{k}': (np.nan if g is None else g[k]) for k in ('to_flies', 'fly_to_fly', 'ratio')})
            rec.update(row.to_dict())
            rows.append(rec)
    return pd.DataFrame(rows)


def gate_wirings(readings: Mapping[str, R.ModelReading], reading_name: str) -> pd.DataFrame:
    """The measure gate's wiring check (the planner's reading, chat.md 2026-11-24 11:34): each wiring's population
    scored against the other wirings taken as the flies, in their own one-fly units (hal's `fly_reference` on the
    others, `score`), per (stage, panel, condition); `passes` where the mean over the wirings is at GATE_WIRING_MAX or
    less, the range over the wirings and the number above GATE_WIRING_MAX beside (the lock's rules -> gate ->
    A_wiring_check). No flies are read."""
    keys = sorted({k for r in readings.values() for k in r.populations})
    rows = []
    for stage, panel, cond in keys:
        pops = [r.populations[(stage, panel, cond)] for r in readings.values() if (stage, panel, cond) in r.populations]
        if len(pops) < 3:
            continue
        od = list(pops[0].values.columns)
        edges = breadth_edges(len(od))
        rk = []
        for i, pop in enumerate(pops):
            ref = P.fly_reference(pops[:i] + pops[i + 1:], od, edges)
            rk.append(dict(wiring=str(pop.name).partition('/')[2] or str(pop.name), stage=stage, panel=panel,
                           condition=cond, reading=reading_name, **P.score(pop, ref).as_series().to_dict()))
        errs = np.array([r['error'] for r in rk], dtype=float)
        for r in rk:
            r.update(mean_error=float(np.nanmean(errs)), min_error=float(np.nanmin(errs)),
                     max_error=float(np.nanmax(errs)), n_above=int((errs > GATE_WIRING_MAX).sum()),
                     passes=bool(np.nanmean(errs) <= GATE_WIRING_MAX))
        rows += rk
    return pd.DataFrame(rows)


def save_populations(readings: Mapping[str, R.ModelReading], part: str, source: str) -> None:
    """Every population read (values) as outputs/.../populations/<part>/<source>__<model>__<wiring>__<stage>__<panel>__
    <condition>.parquet."""
    d = local_dir(f'populations/{part}')
    for r in readings.values():
        for (stage, panel, cond), pop in r.populations.items():
            v = pop.values.copy()
            v.columns = [str(c) for c in v.columns]
            v.index = _str_index(v.index)
            v.to_parquet(d / f"{source.replace(':', '-')}__{r.model}__{r.wiring}__{stage}__{panel}__{cond}.parquet")


def _str_index(idx: pd.Index) -> pd.Index:
    """An index with its labels as strings for parquet: a MultiIndex (e.g. claws' (kc_id, claw_id) or pooled units)
    keeps its levels and names, each level's labels cast to str (a MultiIndex cannot be cast with astype)."""
    if isinstance(idx, pd.MultiIndex):
        return pd.MultiIndex.from_tuples([tuple(str(x) for x in t) for t in idx], names=idx.names)
    return idx.astype(str)


# ---------------------------------------------------------------------------
# step 2: the saved runs

def needs_claw_bouton(readings: Optional[Sequence[str]] = None) -> bool:
    """Whether the claws' readings to be read (each stage's default and alternatives, or those named) include a
    pool by bouton ('b<n>'), which needs each claw's bouton in a saved E90 run's metadata."""
    return any(str(rd.pool).startswith('b') for _, rd in _reading_sets(['kc_claw'], readings))


def step2_runs(source: Source, with_claw_bouton: bool = False) -> List[List[R.SavedRun]]:
    """A source's saved runs, grouped into families read together (E90's eight wirings; E107's eight per panel and
    reading; the D1 fits; B5's one run); a run whose files are missing is left out and reported. `with_claw_bouton`
    loads E90's runs with each claw's bouton (`R.load_e90_run`), for the claws' 'b<n>' pools (`needs_claw_bouton`)."""
    groups: Dict[tuple, List[R.SavedRun]] = {}
    for tag in source.tags:
        try:
            if source.source == 'E107:B':
                panel, reading, wiring = tag.split('_', 2)
                run, key = R.load_e107_run(panel, reading, wiring), (panel, reading)
            else:
                cb = with_claw_bouton and 'kc_claw' in source.stages
                run, key = R.load_e90_run(tag, source.part, with_claw_bouton=cb), (source.part,)
        except FileNotFoundError as e:
            print(f'{source.source} {tag}: not found ({e})', flush=True)
            continue
        groups.setdefault(key, []).append(run)
    return list(groups.values())


def _reading_sets(stages: Sequence[str], names: Optional[Sequence[str]]) -> List[Tuple[str, R.Reading]]:
    """(stage, reading) pairs: each stage's default and alternatives, or those named."""
    out = []
    for s in stages:
        for rd in R.readings_of(s):
            if names is None or rd.name in names:
                out.append((s, rd))
    return out


def part_step2(readings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 2: every source of STEP2_SOURCES read in each stage's default reading and its alternatives (or those named
    in `readings`), scored (`score_family`, `gate_wirings`); step2_scores.csv, step2_gate_wirings.csv,
    curves_step2.json, provenance_step2.json; the populations of the default readings under
    outputs/.../populations/step2/."""
    lock = require_lock()
    check_wirings()
    t0 = time.time()
    flies = R.FlyCache()
    scores, gates, curves = [], [], {}
    claw_bouton = needs_claw_bouton(readings)
    for src in STEP2_SOURCES:
        for runs in step2_runs(src, with_claw_bouton=claw_bouton):
            for stage, rd in _reading_sets([s for s in src.stages if s in runs[0].stages()], readings):
                fam = R.read_family(runs, stages=[stage], readings={stage: rd}, lock=lock, flies=flies)
                sc = score_family(fam, rd.name, flies)
                if len(sc):
                    sc.insert(0, 'source', src.source)
                    sc.insert(1, 'group', runs[0].tag if len(runs) == 1 else '_'.join(runs[0].tag.split('_')[:2])
                              if src.source == 'E107:B' else src.part)
                    scores.append(sc)
                g = gate_wirings(fam, rd.name)
                if len(g):
                    g.insert(0, 'source', src.source)
                    gates.append(g)
                for w, r in fam.items():
                    for k, cv in r.curves.items():
                        curves[f"{src.source}|{runs[0].tag}|{w}|{'|'.join(k)}|{rd.name}"] = cv.as_dict()
                if rd is R.DEFAULT_READINGS[stage]:
                    save_populations(fam, 'step2', src.source)
                print(f'step2 {src.source} {runs[0].tag} {stage} {rd.name}: {time.time() - t0:.0f} s', flush=True)
                gc.collect()
    tab = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    tab.to_csv(out_dir() / 'step2_scores.csv', index=False)
    (pd.concat(gates, ignore_index=True) if gates else pd.DataFrame()).to_csv(out_dir() / 'step2_gate_wirings.csv',
                                                                             index=False)
    _write_json(out_dir() / 'curves_step2.json', curves)
    provenance('step2', lock=lock, seconds=round(time.time() - t0, 1),
               sources=[s._asdict() for s in STEP2_SOURCES])
    return tab


# ---------------------------------------------------------------------------
# step 3: the points

def tag_of(pt: Point) -> str:
    """'<panel>_<reading>_<model>' (E107's `run_tag` form), e.g. 'megamat_fitted_real', 'validation2_fitted_reference'."""
    if pt.panel not in R.INPUT_PANELS or pt.reading not in STEP3_READINGS:
        raise ValueError(f'{pt}')
    return f'{pt.panel}_{pt.reading}_{pt.model}'


def parse_tag(tag: str) -> Point:
    """The point of a tag (`tag_of`'s inverse)."""
    panel, reading, model = tag.split('_', 2)
    pt = Point(panel, reading, model)
    tag_of(pt)
    return pt


def is_family(pt: Point) -> bool:
    return pt.model in WIRINGS


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
    if panel == 'apl16':
        return H.apl16_orn_deltas()
    if panel in ('megamat', 'validation2'):
        from experiments import E107_validation2_panel as E107
        return E107.panel_input(panel)
    if panel == 'mixtures':
        lab = H.apl16_orn_deltas()
        n = mixture_input()
        missing = [g for g in lab.index if g not in n.index]
        if missing:
            raise RuntimeError(f"the mixture input lacks the lab input's glomeruli {missing}")
        out = n.loc[list(lab.index)]
        if out.isna().any().any():
            raise RuntimeError('mixtures: NaN in the input')
        return out
    raise ValueError(f'{panel=}')


def mixture_input() -> pd.DataFrame:
    """pyolfsysm's `natmix_orn_deltas()` (diagnostics dropped) restricted to the two mixture panels' stimuli (viki's
    KMIX_ODORS and CMIX_ODORS, through `H.odor_abbrev`: the binary mixtures and the air mixes, which the flies' panels
    leave out, are not run), all its glomeruli."""
    from data_prep import stage_population_data as D
    from pyolfsysm.mb_model import natmix_orn_deltas
    n = natmix_orn_deltas()
    want = {o for p in R.MIXTURE_PANELS for o in D.panel_odors(p)}
    keep = [c for c in n.columns if H.odor_abbrev(c[-1] if isinstance(c, tuple) else c) in want]
    if len(keep) != len(want):
        raise RuntimeError(f'the mixture input has {len(keep)} of the panels\' {len(want)} stimuli')
    return n[keep]


def published_input(panel: str) -> pd.DataFrame:
    """The published models' input: the paper's megamat and validation2 inputs (`paper_megamat_orn_deltas`,
    `paper_validation2_orn_deltas`, all their glomeruli), `natmix_orn_deltas()` for the mixtures."""
    from pyolfsysm import mb_model as M
    if panel == 'megamat':
        return M.paper_megamat_orn_deltas()
    if panel == 'validation2':
        return M.paper_validation2_orn_deltas()
    if panel == 'mixtures':
        return mixture_input()
    raise ValueError(f'{panel=}')


def held_rate(panel: str, lock: Optional[Mapping] = None) -> float:
    """The APL-on response rate a held point is held at, from viki's lock, as HELD_RATE_SOURCE says: 'remy_criterion'
    Remy's four megamat flies' rate by her criterion ('remy_megamat_criterion_rate' -> 'mean'); 'lock_criterion' the
    panel's KC somata populations at the lock's criterion (each population's 'response_rate' -> 'mean_over_flies'; the
    two mixture panels weighted by their numbers of stimuli, as the model's mixture input holds both panels' stimuli).
    RuntimeError if the lock lacks them."""
    lk = lock if lock is not None else require_lock()
    src = HELD_RATE_SOURCE.get(panel)
    if panel == 'megamat' and src == 'remy_criterion':
        v = (lk.get('remy_megamat_criterion_rate') or {}).get('mean')
        if v is None:
            raise RuntimeError("the lock has no 'remy_megamat_criterion_rate'")
        return float(v)
    if panel in ('megamat', 'mixtures') and src == 'lock_criterion':
        num = den = 0.0
        for p in (R.MIXTURE_PANELS if panel == 'mixtures' else ('megamat',)):
            L = (lk.get('populations') or {}).get(f'kc_soma__{p}__wild_type')
            if L is None:
                raise RuntimeError(f'the lock has no KC somata of {p}')
            n = len(L.get('odors') or ())
            num += n * float(L['response_rate']['mean_over_flies'])
            den += n
        return num / den
    raise ValueError(f'{panel}: no held response rate')


def held_targets(pt: Point, size: int, lock: Optional[Mapping] = None) -> Tuple[int, int, float, float]:
    """(target_off, target_on, fraction_off, fraction_on) of a held point as response counts on its KCs x odors
    (`size`): APL on at `held_rate`, APL off by HELD_OFF_RULE ('ratio': times the wiring's own fit's APL-off /
    APL-on response counts; 'e107': the wiring's own fit's APL-off fraction)."""
    from experiments import E83_cross_wiring_fit as E83
    f_on = held_rate(pt.panel, lock)
    F = E83.saved_fit(pt.model)
    f_off = f_on * F['target_off'] / F['target_on'] if HELD_OFF_RULE == 'ratio' else F['target_off'] / F['size']
    return int(round(f_off * size)), int(round(f_on * size)), float(f_off), float(f_on)


@dataclass
class Built:
    """A point's model, built and ready to run at fixed parameters.

    Attributes:
        sim: params -> (responses, spike counts), the simulator's outputs as it gives them.
        params, params_off: the parameters of the APL-on and acutely APL-off runs (absolute).
        vars: () -> (olfsysm's mp, rv and the odor order: output column j is olfsysm's odor order[j]) after the last
            simulation.
        W: the claw x bouton wiring (None for models without boutons).
        glomeruli: the PN stage's glomeruli, in olfsysm's row order.
        stages: the stage matrices read (of R.STAGE_FILES).
        orn: the input.
        meta: the units' metadata (E90's units.json: bouton_glomerulus, claw_glomerulus, kc_of_claw, claw_bouton).
        info: the model's description, sources and kwargs.
    """
    sim: Callable
    params: Dict[str, float]
    params_off: Dict[str, float]
    vars: Callable
    W: Optional[pd.DataFrame]
    glomeruli: List[str]
    stages: Tuple[str, ...]
    orn: pd.DataFrame
    meta: Dict[str, np.ndarray] = field(default_factory=dict)
    info: Dict[str, object] = field(default_factory=dict)


def _units_meta(W: pd.DataFrame) -> Dict[str, np.ndarray]:
    from experiments import E90_model_stages as E90
    return dict(bouton_glomerulus=np.asarray([str(x) for x in W.columns.get_level_values('glomerulus')]),
                claw_glomerulus=np.asarray([str(x) for x in E90.glomerulus_of_claws(W)]),
                kc_of_claw=np.asarray([int(x) for x in W.index.get_level_values(H.KC_ID_LEVEL)]),
                claw_bouton=R.claw_boutons(W))


def _simulator_built(sim, cells: Mapping, W: pd.DataFrame, params: dict, orn: pd.DataFrame, stages: Tuple[str, ...],
                     info: dict) -> Built:
    """A `Built` of a fixed-parameter simulator (its closure's mp, rv and odor order; APL's activity saved)."""
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    mp.kc.save_inh_sims = True
    wp = np.asarray(rv.kc.wPNKC, dtype=float)
    if wp.shape != W.shape or not np.array_equal(wp, W.to_numpy(dtype=float)):
        raise RuntimeError(f"{info.get('model')}: olfsysm's wPNKC differs from the wiring's claw x bouton matrix")
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    gl = list(dict.fromkeys(str(g) for g in W.columns.get_level_values('glomerulus')))
    p = dict(params)
    return Built(sim=sim, params=p, params_off=dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                 vars=lambda: (mp, rv, order), W=W, glomeruli=gl, stages=stages, orn=orn, meta=_units_meta(W),
                 info=info)


def _override_record(orn: Optional[pd.DataFrame], extra_kws: Optional[Mapping]) -> dict:
    """What a run built with overrides records in its info: the input's sha256 when `orn` was given, and each of
    `extra_kws` as its value, or as its sha256 when it is a table (e.g. E135's pn_rates)."""
    from experiments.E77_input_screen import frame_sha256
    rec = {}
    if orn is not None:
        rec['orn_override_sha256'] = frame_sha256(orn)
    for k, v in dict(extra_kws or {}).items():
        rec[f'extra_kws.{k}'] = (f'sha256:{frame_sha256(v)}' if isinstance(v, (pd.DataFrame, pd.Series))
                                 else v if np.isscalar(v) or v is None else str(v))
    return rec


def build_family(label: str, panel: str, params: Optional[Mapping[str, float]] = None,
                 orn: Optional[pd.DataFrame] = None, extra_kws: Optional[Mapping] = None) -> Built:
    """One wiring of the family on a panel's input (E83's `build_simulator`, the breadth held at the lab input's, as
    E107's `run_fitted`), at its fit's parameters or `params`.

    Args:
        label: one of E83's wirings.
        panel: the panel whose input is used when `orn` is None (`panel_input`).
        params: the parameters (default the wiring's fit's).
        orn: an ORN input in place of the panel's (e.g. E135's solved input); the breadth stays the lab input's.
        extra_kws: further `make_fixed_param_kc_simulator` kwargs passed through E83's `build_simulator` (e.g.
            E135's `dict(pn_rates=<glomeruli x odors>)`, olfsysm's preset PN rates). Both overrides are recorded in
            the Built's info (`_override_record`).
    """
    from experiments import E83_cross_wiring_fit as E83
    Pp = E83.wiring_parts(label)
    given = orn
    orn = panel_input(panel) if orn is None else orn
    p = dict(Pp['p'] if params is None else params)
    sim, info = E83.build_simulator(label, orn=orn, breadth=Pp['breadth'], params=p,
                                    extra_kws=None if not extra_kws else dict(extra_kws))
    return _simulator_built(sim, info['cells'], info['W'], p, orn, R.STAGE_FILES,
                            dict(model=FAMILY, wiring=label, fit_folder=Pp['folder'], fit_tag=Pp['tag'],
                                 kws=info['kws'], breadth='held at the lab input\'s (E83 wiring_parts)',
                                 **_override_record(given, extra_kws)))


def _fit_mb_model_built(name: str, orn: pd.DataFrame, kws: Mapping, params: Mapping[str, float], info: dict) -> Built:
    """A `Built` that calls `fit_mb_model` at fixed parameters for each run (models the fixed-parameter simulator does
    not support, e.g. pyolfsysm's uniform wiring): one call per condition, the olfsysm variables kept; PN stage and
    somata only. APL off is wAPLKC 0."""
    from pyolfsysm import mb_model as M
    tuning = getattr(M, '_tuning_only_kws', ('target_sparsity', 'target_sparsity_factor_pre_APL', 'max_iters', 'sp_acc',
                                            'thr_sp_lr_coeff', 'sp_lr_coeff', 'linear_lr_falloff',
                                            'binary_search_on_overshoot', 'hardcode_initial_sp'))
    base = {k: v for k, v in kws.items() if k not in tuning}
    last: dict = {}

    def sim(p: Mapping[str, float]):
        r, s, wPNKC, pdict = M.fit_mb_model(orn, fixed_thr=float(p['fixed_thr']), wAPLKC=p['wAPLKC'],
                                            wKCAPL=p['wKCAPL'], return_olfsysm_vars=True, **base)
        last.update(mp=pdict['mp'], rv=pdict['rv'], wPNKC=wPNKC)
        return r, s

    p = dict(params)
    b = Built(sim=sim, params=p, params_off=dict(p, wAPLKC=0.0), vars=lambda: (last['mp'], last['rv'],
              np.arange(orn.shape[1])), W=None, glomeruli=[], stages=('pn',), orn=orn,
              info=dict(info, model=name, kws={k: str(v) for k, v in base.items()}))
    b.info['_last'] = last
    return b


def pn_glomeruli(built: Built, n_rows: int) -> List[str]:
    """The labels of the PN stage's rows: the wiring's glomeruli (in order of first appearance), or for a model built
    through `fit_mb_model` its wPNKC's glomeruli or its input's, whichever has the PN stage's number of rows."""
    if built.glomeruli:
        return list(built.glomeruli)
    cands = []
    w = (built.info.get('_last') or {}).get('wPNKC')
    if w is not None:
        cols = w.columns
        g = cols.get_level_values('glomerulus') if 'glomerulus' in (cols.names or []) else cols
        cands.append(list(dict.fromkeys(str(x) for x in g)))
    cands.append([str(x) for x in built.orn.index])
    for c in cands:
        if len(c) == n_rows:
            return c
    raise RuntimeError(f'{n_rows} PN rows match no glomerulus list ({[len(c) for c in cands]})')


def _paper_params(d: Path) -> Dict[str, float]:
    """fixed_thr, wAPLKC, wKCAPL of a published model's params.csv (two columns, key and value)."""
    s = pd.read_csv(d / 'params.csv', header=None, index_col=0).iloc[:, 0]
    return {k: float(s[k]) for k in ('fixed_thr', 'wAPLKC', 'wKCAPL')}


def _uniform_params(seed: int = UNIFORM_SEED) -> Dict[str, float]:
    """The paper's uniform draw `seed`'s tuned parameters (pyolfsysm/data/sent_to_anoop/2024-05-16: the params CSV's
    per-seed lists, in the order of the seeds of its responses CSV)."""
    from pyolfsysm.mb_model import paper_uniform_output_dir as d
    s = pd.read_csv(d / 'megamat_uniform_model_params.csv', header=None, index_col=0).iloc[:, 0]
    seeds = sorted(pd.read_csv(d / 'megamat_uniform_model_responses_n-seeds_100.csv', usecols=['seed'])['seed'].unique())
    i = seeds.index(int(seed))
    vals = {k: ast.literal_eval(str(s[k])) for k in ('fixed_thr', 'wAPLKC', 'wKCAPL')}
    if any(len(v) != len(seeds) for v in vals.values()):
        raise RuntimeError("the uniform model's parameter lists are not one per seed")
    return {k: float(v[i]) for k, v in vals.items()}


def build_beside(model: str, panel: str) -> Built:
    """A beside model on a panel's input (BESIDE_MODELS; the module docstring says how each is defined)."""
    from pyolfsysm import mb_model as M
    if model == 'published_hemibrain':
        from pyolfsysm.util import hemibrain_paper_repro_kws
        d = M.paper_validation2_hemibrain_output_dir if panel == 'validation2' else M.paper_hemibrain_output_dir
        return _fit_mb_model_built(model, published_input(panel), dict(hemibrain_paper_repro_kws,
                                   pn2kc_connections='hemibrain'), _paper_params(d),
                                   dict(params_from=str(d / 'params.csv')))
    if model == 'uniform_draw':
        from pyolfsysm.util import paper_repro_kws
        # (_drop_glom_with_plus=False: the setting under which pyolfsysm's test_uniform_paper_repro matches the saved
        # uniform outputs)
        return _fit_mb_model_built(model, published_input(panel), dict(paper_repro_kws, pn2kc_connections='uniform',
                                   n_claws=7, seed=UNIFORM_SEED, _drop_glom_with_plus=False), _uniform_params(),
                                   dict(params_from=str(M.paper_uniform_output_dir / 'megamat_uniform_model_params.csv'),
                                        seed=UNIFORM_SEED))
    orn = panel_input(panel)
    if model == 'bouton_baseline':
        from experiments import E10_joint_fit as J
        from experiments import E83_cross_wiring_fit as E83
        from experiments.E12_breadth_apl import e3_bouton_params
        p = e3_bouton_params()
        sim, _ = M.make_fixed_param_kc_simulator(orn, initial_params=p, **J.base_kws(None))
        cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
        return _simulator_built(sim, cells, E83.wiring_parts('real')['W'], p, orn, R.STAGE_FILES,
                                dict(model=model, params_from="E12_breadth_apl.e3_bouton_params (E3's bouton_pnapl_x1)",
                                     kws='E10_joint_fit.base_kws(None)'))
    if model == 'reference':
        from experiments import E69_kc_firing_ceiling as E69
        m = E69.Model('ref40')
        sim, _ = M.make_fixed_param_kc_simulator(orn, initial_params=m.p, **m.kws)
        cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
        return _simulator_built(sim, cells, m.W, m.p, orn, ('pn', 'bouton'),
                                dict(model=model, spec=REFERENCE_SPEC, base=m.base, info=m.info,
                                     note="its claws pass a 40 Hz per-claw threshold: no claw stage is saved"))
    if model == 'firing_ceiling':
        from experiments import E82_protocol_runs as E82
        m = E82.model_parts('fit1')
        W = None if m['wiring'] is None else E82.wiring_parts(m['wiring'])['W']
        sim, info = E82.build_simulator(m['kc'], wPNKC=W, breadth=E82.own_breadth(m), params=m['p'], orn=orn)
        return _simulator_built(sim, info['cells'], info['wPNKC'], m['p'], orn, R.STAGE_FILES,
                                dict(model=model, source=list(m['source']), kc=m['kc']))
    raise ValueError(f'{model=}')


def build_point(pt: Point, params: Optional[Mapping[str, float]] = None, orn: Optional[pd.DataFrame] = None,
                extra_kws: Optional[Mapping] = None) -> Built:
    """A point's model on its panel's input: the family's wiring (`build_family`, at `params` for a held point's second
    build, with `orn` / `extra_kws` overrides for other experiments' points, e.g. E135's step 2) or a beside model
    (`build_beside`, which takes no overrides)."""
    if is_family(pt):
        return build_family(pt.model, pt.panel, params, orn=orn, extra_kws=extra_kws)
    if params is not None or orn is not None or extra_kws:
        raise ValueError('a beside model runs at its own parameters and input')
    return build_beside(pt.model, pt.panel)


def read_model_stages(mp, rv, W: Optional[pd.DataFrame], order: np.ndarray, odors: Sequence[str],
                      glomeruli: Sequence[str], stages: Sequence[str]) -> Dict[str, pd.DataFrame]:
    """The stage matrices of the last simulation, E90's `read_stages` with the settled pre-odor window beside: per stage
    of `stages` ('pn' glomeruli, 'bouton' boutons, 'claw' claws; x odors) the means over MODEL_WINDOW ('<stage>_win'),
    E90's PRE_WINDOW ('<stage>_pre') and SETTLED_PRE_WINDOW ('<stage>_pre_settled'), the claws as E90's
    `claw_inputs`, and 'apl' (APL's activity over the three windows, per odor) where APL's activity was saved."""
    from experiments import E90_model_stages as E90
    from pyolfsysm.mb_model import get_sim_var
    t = E90.step_times(mp, rv)
    w = E90._windows(t, R.MODEL_WINDOW)
    win = dict(win=w['win'], pre=w['pre'],
               pre_settled=np.flatnonzero((t >= R.SETTLED_PRE_WINDOW[0] - 1e-9) & (t < R.SETTLED_PRE_WINDOW[1] - 1e-9)))
    pn = get_sim_var(rv.pn, 'pn_sims')
    if np.asarray(pn[0]).shape[0] != len(glomeruli):
        raise RuntimeError(f'{np.asarray(pn[0]).shape[0]} PN rows for {len(glomeruli)} glomeruli')
    want_b = 'bouton' in stages or 'claw' in stages
    if want_b:
        btn = get_sim_var(rv.pn, 'bouton_sims')
        inh = get_sim_var(rv.kc, 'inh_sims')
        wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
        Wm = W.to_numpy(dtype=float)
        if btn.shape[1] != W.shape[1] or len(wk) != W.shape[0]:
            raise RuntimeError(f'bouton_sims {btn.shape}, wAPLKC {len(wk)}, wiring {W.shape}')
    else:
        try:
            inh = get_sim_var(rv.kc, 'inh_sims')
        except Exception:
            inh = None
    out = {f'{k}_{x}': {} for k in stages for x in win}
    apl = {}
    for j, o in enumerate(odors):
        i = int(order[j])
        for x, idx in win.items():
            out[f'pn_{x}'][o] = np.asarray(pn[i])[:, idx].mean(axis=1)
            if want_b:
                b = np.asarray(btn[i])[:, idx]
                if np.isnan(b).any():
                    raise RuntimeError(f'NaN bouton rates in the {x} window')
                if 'bouton' in stages:
                    out[f'bouton_{x}'][o] = b.mean(axis=1)
                if 'claw' in stages:
                    out[f'claw_{x}'][o] = E90.claw_inputs(Wm, b, np.asarray(inh[i]).ravel()[idx - 1], wk).mean(axis=1)
        if inh is not None and len(inh):
            a = np.asarray(inh[i]).ravel()
            apl[o] = {x: float(a[idx].mean()) for x, idx in win.items()}
    res = {}
    for k, v in out.items():
        st = k.split('_')[0]
        idx = (pd.Index(list(glomeruli), name='glomerulus') if st == 'pn' else W.columns if st == 'bouton' else W.index)
        res[k] = pd.DataFrame(v, index=idx)
    if apl:
        res['apl'] = pd.DataFrame(apl).T
    return res


def run_stages(pt: Point, built: Built, params: Mapping[str, float], params_off: Mapping[str, float]) -> Dict[str, dict]:
    """Both conditions of a point at the given parameters, each read by `read_model_stages` (the stages its model has),
    with the responses and spike counts: {'on', 'off'}: the stage dict."""
    from pyolfsysm.mb_model import get_sim_var
    out = {}
    for c, p in (('on', params), ('off', params_off)):
        r, s = built.sim(dict(p))
        mp, rv, order = built.vars()
        odors = list(R._abbrev(r.iloc[:0]).columns)
        gl = pn_glomeruli(built, int(np.asarray(get_sim_var(rv.pn, 'pn_sims')[0]).shape[0]))
        st = read_model_stages(mp, rv, built.W, order, odors, gl, built.stages)
        st.update(responses=r, spike_counts=s)
        out[c] = st
    return out


def replay_beside(model: str) -> dict:
    """A beside model's replay on its fit panel on this host (cached in results/.../replay_<model>_<host>.json): built
    there (`build_beside`) and run at its parameters, its outputs must equal its saved ones exactly (spike counts, or
    responses where only those were saved). RuntimeError otherwise."""
    path = out_dir() / f'replay_{model}_{_host()}.json'
    rec = _read_json(path)
    if rec.get('exact'):
        return rec
    t0 = time.time()
    if BESIDE_FIT_PANEL[model] == 'megamat':
        b = build_beside(model, 'megamat')
        r, s = b.sim(dict(b.params))
        checks = published_checks(model, 'megamat', dict(responses=r, spike_counts=s))
        rec = dict(model=model, panel='megamat', checks=checks, exact=all(c['differ'] == 0 for c in checks),
                   seconds=round(time.time() - t0, 1), host=_host(), time=_now(), **H.code_commits())
        _write_json(path, rec)
        del b
        gc.collect()
        if not rec['exact']:
            raise RuntimeError(f'{model}: the build does not reproduce its published outputs on megamat: {rec}')
        return rec
    saved = saved_beside(model)
    b = build_beside(model, 'apl16')
    rec = dict(model=model, panel='apl16')
    for c, p in (('on', b.params), ('off', b.params_off)):
        r, s = b.sim(dict(p))
        rs, ss = saved[c]
        rec[c] = dict(counts_differ=_differ(R._abbrev(s), R._abbrev(ss)),
                      responses_differ=_differ(R._abbrev(r), R._abbrev(rs)))
    rec['exact'] = all(rec[c]['counts_differ'] == 0 and rec[c]['responses_differ'] == 0 for c in ('on', 'off'))
    rec.update(seconds=round(time.time() - t0, 1), host=_host(), time=_now(), **H.code_commits())
    _write_json(path, rec)
    del b
    gc.collect()
    if not rec['exact']:
        raise RuntimeError(f'{model}: the build does not reproduce its saved outputs on the APL16 input: {rec}')
    return rec


def saved_beside(model: str) -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    """A beside model's saved outputs on its fit panel ({'on', 'off'}: (responses, spike counts))."""
    if model == 'bouton_baseline':
        d = H.shared_dir().joinpath(*E3_SAVED[:-1])
        n = E3_SAVED[-1]
        return {c: tuple(pd.read_parquet(d / f'{n}_{cc}_{k}.parquet') for k in ('responses', 'spike_counts'))
                for c, cc in R.CONDITION_FILE.items()}
    if model == 'reference':
        from experiments import E43_model_trial_noise as E43
        from experiments import E69_kc_firing_ceiling as E69
        b = E69.bases()['ref40']
        fm = E43.fit_model(E43.model_from(f"{b['exp']}:{b['tag']}"), H.apl16_orn_deltas(), E69.base_dir())
        return {c: tuple(pd.read_parquet(p) for p in fm['saved'][cc]) for c, cc in R.CONDITION_FILE.items()}
    if model == 'firing_ceiling':
        from experiments import E82_protocol_runs as E82
        return dict(E82.model_parts('fit1')['saved'])
    raise ValueError(f'{model}: no saved outputs on the APL16 input')


def published_checks(model: str, panel: str, on: Optional[Mapping[str, pd.DataFrame]]) -> Optional[List[dict]]:
    """The published models' checks of an APL-on run against their saved outputs, by position (the saved tables are
    indexed by KC number): the published hemibrain model's spike counts on megamat and validation2, the uniform draw's
    responses of its seed on megamat. None where the model has no saved outputs on that panel (with `on` None, only
    whether it has)."""
    from pyolfsysm import mb_model as M
    if model == 'published_hemibrain' and panel in ('megamat', 'validation2'):
        if on is None:
            return []
        d = M.paper_validation2_hemibrain_output_dir if panel == 'validation2' else M.paper_hemibrain_output_dir
        sv, x, name = pd.read_csv(d / 'spike_counts.csv', index_col=0), on['spike_counts'], 'published spike counts'
    elif model == 'uniform_draw' and panel == 'megamat':
        if on is None:
            return []
        sv = M.paper_uniform_model_responses().xs(UNIFORM_SEED, level='seed')
        x, name = on['responses'], f'published responses of seed {UNIFORM_SEED}'
    else:
        return None
    return [dict(name=f'{name} (APL on)', differ=-1 if sv.shape != x.shape else
                 int((sv.to_numpy(dtype=float) != x.to_numpy(dtype=float)).sum()))]


def check_identity(pt: Point, out: Mapping[str, dict], hold: Optional[dict] = None) -> dict:
    """A point's identity before its stages are saved (0 differing counts and responses; `exact`):
    - a held point: both conditions equal the hold's final outputs;
    - the family as fitted: equal E107's saved run of the same tag where there is one (megamat), and the saved fit on
      the odors the panels share (E107's `run_identity`; the mixtures share none);
    - the published hemibrain model on megamat / validation2: its APL-on spike counts equal the saved ones;
      the uniform draw on megamat: its APL-on responses equal the saved ones of its seed;
    - the other beside models: their replay on the fit panel (`replay_beside`) stands for the build."""
    rec = dict(tag=tag_of(pt), checks=[])

    def add(name, a, b):
        rec['checks'].append(dict(name=name, differ=_differ(a, b)))

    if hold is not None:
        for c in ('on', 'off'):
            add(f'hold {c} counts', out[c]['spike_counts'], hold[c][1])
            add(f'hold {c} responses', out[c]['responses'], hold[c][0])
    elif is_family(pt):
        try:
            e = R.load_e107_run(pt.panel, 'fitted', pt.model)
        except (FileNotFoundError, ValueError):
            e = None
        if e is not None:
            for c in ('on', 'off'):
                add(f'E107 {c} counts', R._abbrev(out[c]['spike_counts']), R._abbrev(e.conditions[c]['spike_counts']))
                add(f'E107 {c} responses', R._abbrev(out[c]['responses']), R._abbrev(e.conditions[c]['responses']))
        if pt.panel in ('megamat', 'validation2'):
            from experiments import E107_validation2_panel as E107
            o = {c: (R._abbrev(out[c]['responses']), R._abbrev(out[c]['spike_counts'])) for c in ('on', 'off')}
            ri = E107.run_identity(pt.model, o, pt.panel)
            rec['saved_fit'] = ri
            rec['checks'].append(dict(name='saved fit, shared odors', differ=0 if ri['reproduces'] else 1))
        else:
            rec['note'] = 'the mixtures share no odor at the same concentration with the fit panel: no saved counts'
    elif published_checks(pt.model, pt.panel, None) is not None:
        rec['checks'] += published_checks(pt.model, pt.panel, out['on'])
    else:
        rec['note'] = 'no saved outputs on this panel; the replay on the fit panel stands for the build'
        if pt.model in BESIDE_MODELS:
            rec['replay'] = replay_beside(pt.model)
    rec['exact'] = all(c['differ'] == 0 for c in rec['checks'])
    return rec


def save_run(tag: str, out: Mapping[str, dict], built: Built, info: Mapping[str, object],
             archive: bool = True) -> List[Path]:
    """A point's stage matrices, responses, spike counts, APL means, units' metadata and input, locally
    (outputs/.../step3/) and in the model-output archive (`stage_readings.save_step3_run`)."""
    dirs = [local_dir(R.STEP3_PART)] + ([archive_dir() / R.STEP3_PART] if archive else [])
    return R.save_step3_run(tag, out, built.meta, built.orn, info, dirs)


def load_run(tag: str) -> Optional[R.SavedRun]:
    """A saved step-3 point (locally, else the archive) as `stage_readings.SavedRun`; None if absent."""
    try:
        return R.load_step3_run(tag)
    except FileNotFoundError:
        return None


def runs_table() -> pd.DataFrame:
    """Every row of every host's and lane's runs_3_*.csv, the last row per tag."""
    fs = [read_runs_csv(p) for p in sorted(out_dir().glob('runs_3_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def _counts_row(out: Mapping[str, dict]) -> dict:
    r = {}
    for c in ('on', 'off'):
        resp = out[c]['responses'].to_numpy(dtype=float)
        cnt = out[c]['spike_counts'].to_numpy(dtype=float)
        r.update({f'n_{c}': int((resp > 0.5).sum()), f'rate_{c}': float((resp > 0.5).mean()),
                  f'max_count_{c}': float(cnt.max())})
    r['rate_ratio'] = r['rate_off'] / r['rate_on'] if r['rate_on'] > 0 else np.nan
    return r


def run_point(pt: Point, lane: str, table: Path, prov_path: Path, reclaim: bool = False,
              lock: Optional[Mapping] = None, done: Optional[set] = None) -> Optional[dict]:
    """One point, claimed (None if done or claimed elsewhere): built, held if its reading is 'held' (then built again at
    the held parameters), run with its stages, checked (`check_identity`), saved (`save_run`), its row appended to
    `table` (absolute parameters, factors and targets, response counts and rates, largest counts, seconds, host, lane,
    threads, the input's sha256, the commits) and its provenance to `prov_path`; the claim is released on failure."""
    from experiments.E77_input_screen import frame_sha256, hold_rates
    tag = tag_of(pt)
    if done is not None and tag in done:
        return None
    if not claim(tag, lane, reclaim):
        print(f'E134 {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        hold, extra = None, {}
        if pt.reading == 'held':
            if not is_family(pt):
                raise ValueError(f'{tag}: only the family is held')
            from experiments import E83_cross_wiring_fit as E83
            Pp = E83.wiring_parts(pt.model)
            orn = panel_input(pt.panel)
            n_kc = E83.saved_fit(pt.model)['saved']['on'][0].shape[0]
            t_off, t_on, f_off, f_on = held_targets(pt, n_kc * orn.shape[1], lock)
            hold = hold_rates(lambda: E83.build_simulator(pt.model, orn=orn, breadth=Pp['breadth'], params=dict(Pp['p'])),
                              target_off=t_off, target_on=t_on, params=dict(Pp['p']), label=f'E134 {tag}')
            built = build_family(pt.model, pt.panel, hold['params'], orn=orn)
            params, params_off = dict(hold['params']), dict(hold['params_off'])
            extra = dict(thr_factor=hold['thr_factor'], apl_factor=hold['apl_factor'], target_on=t_on, target_off=t_off,
                         frac_on_target=f_on, frac_off_target=f_off, at_bound_thr=hold['at_bound_thr'],
                         at_bound_apl=hold['at_bound_apl'], hold_n_sims=hold['n_sims'], held_off_rule=HELD_OFF_RULE,
                         fit_fixed_thr=float(Pp['p']['fixed_thr']))
            assert tuple(extra) == HELD_EXTRA_COLUMNS
        else:
            built = build_point(pt)
            params, params_off = built.params, built.params_off
        out = run_stages(pt, built, params, params_off)
        ident = check_identity(pt, out, None if hold is None else {c: hold[c] for c in ('on', 'off')})
        id_path = out_dir() / f'identity_3_{_host()}_{lane}.json'
        ids = _read_json(id_path)
        ids[tag] = ident
        _write_json(id_path, ids)
        if not ident['exact']:
            raise RuntimeError(f'{tag}: the run is not exact: {ident}')
        info = dict(model=FAMILY if is_family(pt) else pt.model, wiring=pt.model if is_family(pt) else pt.model,
                    panel=pt.panel, reading=pt.reading, tag=tag, stages=list(built.stages), params=params,
                    params_off=params_off, host=_host(), lane=lane, time=_now())
        save_run(tag, out, built, info)
        row = dict(tag=tag, panel=pt.panel, reading=pt.reading, model=info['model'], wiring=info['wiring'],
                   **{k: float(v) for k, v in params.items() if np.isscalar(v)},
                   **{f'off_{k}': float(v) for k, v in params_off.items() if np.isscalar(v)},
                   n_kcs=int(out['on']['spike_counts'].shape[0]), n_odors=int(out['on']['spike_counts'].shape[1]),
                   **_counts_row(out), **extra, stages='+'.join(built.stages), seconds=round(time.time() - t0, 1),
                   host=_host(), lane=lane, time=_now(), threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=frame_sha256(built.orn), identity_exact=True,
                   **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
        _append_csv(table, row)
        prov = _read_json(prov_path) or dict(part='step3', lane=lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(info={k: v for k, v in built.info.items() if not str(k).startswith('_')},
                                   identity=ident, params=params, params_off=params_off,
                                   model=H.model_provenance(orn_deltas=built.orn))
        _write_json(prov_path, prov)
        print(f"E134 {tag}: {row['seconds']:.0f} s, rate on {row['rate_on']:.4f} / off {row['rate_off']:.4f}", flush=True)
        return row
    except BaseException:
        release(tag)
        raise
    finally:
        gc.collect()


def part_step3(models: Optional[Sequence[str]] = None, panels: Optional[Sequence[str]] = None,
               readings: Optional[Sequence[str]] = None, lane: Optional[str] = None, reverse: bool = False,
               reclaim: bool = False) -> pd.DataFrame:
    """Step 3, resumable, one point claimed at a time (`planned_step3`, `run_point`); rows in
    runs_3_<host>_<lane>.csv."""
    lock = require_lock()
    check_wirings()
    lane = lane or _host()
    table = out_dir() / f'runs_3_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_3_{_host()}_{lane}.json'
    for pt in planned_step3(models, panels, readings, reverse):
        done = set(runs_table()['tag'].astype(str))
        run_point(pt, lane, table, prov_path, reclaim, lock, done)
    return runs_table()


# ---------------------------------------------------------------------------
# score

def score_groups() -> Dict[str, List[R.SavedRun]]:
    """The runs the score reads, by group: the family's step-3 runs ('<panel>_<reading>'), E107's validation2 runs of
    the family ('validation2_fitted', 'validation2_held'), and each beside model's step-3 runs ('<model>')."""
    groups: Dict[str, List[R.SavedRun]] = {}
    for pt in planned_step3():
        run = load_run(tag_of(pt))
        if run is None:
            print(f'score: {tag_of(pt)} not run yet', flush=True)
            continue
        groups.setdefault(f'{pt.panel}_{pt.reading}' if is_family(pt) else f'{pt.model}_{pt.panel}', []).append(run)
    for reading in STEP3_READINGS:
        runs = []
        for w in WIRINGS:
            try:
                runs.append(R.load_e107_run('validation2', reading, w))
            except FileNotFoundError:
                pass
        if runs:
            groups[f'validation2_{reading}'] = runs
    return groups


def stage_table(scores: pd.DataFrame, panel: str = 'megamat') -> pd.DataFrame:
    """The stage table: the rows of `scores` on `panel` (every group, stage, condition and reading; each wiring, their
    mean and the common error beside the flies' floor and one recording against the others), in stage order."""
    if not len(scores):
        return scores
    t = scores[scores.panel == panel].copy()
    t['_s'] = t.stage.map({s: i for i, s in enumerate(R.STAGES)})
    return t.sort_values(['group', '_s', 'condition', 'reading'], kind='mergesort').drop(columns='_s')


def step_table(readings: Mapping[str, R.ModelReading], flies: R.FlyCache, panel: str = 'megamat',
               condition: str = 'on') -> Tuple[pd.DataFrame, Dict[str, tuple]]:
    """The step table: per step between neighbouring stages of `panel`, on the odors all have and the flies as viki's
    step 1 reads them (`fly_side_read`):
    - level 'entry': per part entry, the model's change (hal's `step_change` on the wirings' mean parts, and the
      wirings' range) beside the flies' change and its 90 % interval (hal's `step_interval`), the description;
    - level 'part_distance': per part and wiring, the RMS distance of the wiring's change from the flies' change over
      the part's entries, which `classify_steps` reads against `step_null_distance`;
    - level 'population': the no-step yardstick (hal's `no_step_score`) beside the model's downstream error (the
      wirings' mean, `wirings_as_flies`).

    Returns:
        (the table, {step: (upstream flies, downstream flies, odors, breadth bins)}, the inputs of
        `step_null_distance`).
    """
    rows, inputs = [], {}
    names = list(readings)
    for up, down in zip(MEGAMAT_STAGES[:-1], MEGAMAT_STAGES[1:]):
        ku, kd = (up, panel, condition), (down, panel, condition)
        fu, fd = flies.get(ku), flies.get(kd)
        ws = [w for w in names if ku in readings[w].populations and kd in readings[w].populations]
        if fu is None or fd is None or not ws:
            continue
        mu = [readings[w].populations[ku] for w in ws]
        md = [readings[w].populations[kd] for w in ws]
        (pu, ou), (pd_, odn) = fly_side_read(fu), fly_side_read(fd)
        od = [o for o in mu[0].values.columns if o in set(ou) and o in set(odn) and o in md[0].values.columns]
        edges = breadth_edges(len(od))
        step = f'{up}->{down}'
        inputs[step] = (pu, pd_, od, edges)
        iv = P.step_interval(pu, pd_, od, edges)
        views_u = [P.population_views(p, od, edges) for p in mu]
        views_d = [P.population_views(p, od, edges) for p in md]
        mean_ch = P.step_change(P.mean_views(views_u), P.mean_views(views_d))
        each = [P.step_change(a, b) for a, b in zip(views_u, views_d)]
        for part, df in iv.items():
            if part not in mean_ch:
                continue
            m = mean_ch[part].reindex(df.index)
            for e in df.index:
                rows.append(dict(step=step, level='entry', part=part, entry=str(e), model_change=float(m[e]),
                                 model_change_lo=float(np.nanmin([x[part].get(e, np.nan) for x in each])),
                                 model_change_hi=float(np.nanmax([x[part].get(e, np.nan) for x in each])),
                                 flies_change=float(df.loc[e, 'change']), flies_lo=float(df.loc[e, 'lo']),
                                 flies_hi=float(df.loc[e, 'hi']),
                                 inside=bool(df.loc[e, 'lo'] <= m[e] <= df.loc[e, 'hi'])))
            for w, x in zip(ws, each):
                d = x[part].reindex(df.index).to_numpy(dtype=float) - df['change'].to_numpy(dtype=float)
                rows.append(dict(step=step, level='part_distance', part=part, wiring=w,
                                 distance=float(np.sqrt(np.nanmean(d ** 2))) if np.isfinite(d).any() else np.nan))
        ref_d = P.fly_reference(pd_, od, edges)
        ns = P.no_step_score(pu, ref_d)
        model_err = float(P.wirings_as_flies(md, ref_d).loc['mean', 'error'])
        rows.append(dict(step=step, level='population', part='population', entry='error',
                         no_step_error=float(ns.error), model_downstream_error=model_err))
    return pd.DataFrame(rows), inputs


def step_nulls(inputs: Mapping[str, tuple]) -> Dict[str, Dict[str, float]]:
    """`step_null_distance` of every step of `step_table`'s inputs."""
    return {step: step_null_distance(*args) for step, args in inputs.items()}


def megamat_stage_error(scores: pd.DataFrame, group: str = 'megamat_fitted', condition: str = 'on') -> pd.DataFrame:
    """The megamat stage error (CLAUDE.md's main goal; the planner's reading, chat.md 2026-11-24 11:34) of a group:
    per wiring, the root mean square of the population errors of the megamat PN dendrites (on four views, hal's
    `error_with_matrix`), PN boutons, KC claws and KC somata (on three, `error`), each stage in its default reading;
    'stage_error_three' beside with the PN dendrites on three views too; the model's value, the mean over its wirings
    of each wiring's stage error (row 'model'); and beside it the stage error of the wirings' mean views (row 'common',
    the error common to the wirings; the lock's rules -> megamat_stage_error). A stage missing for a wiring makes its
    value NaN."""
    s = scores[(scores.group == group) & (scores.panel == 'megamat') & (scores.condition == condition)]
    wirings = [w for w in dict.fromkeys(s.wiring) if w not in ('mean', 'common')] + (
        ['common'] if (s.wiring == 'common').any() else [])
    rows = []
    for w in wirings:
        rec, e4, e3 = dict(group=group, wiring=w, condition=condition), [], []
        for st in MEGAMAT_ERROR_STAGES:
            x = s[(s.stage == st) & (s.wiring == w) & (s.reading == R.DEFAULT_READINGS[st].name)]
            v3 = float(x.error.iloc[0]) if len(x) else np.nan
            v4 = float(x.error_with_matrix.iloc[0]) if (len(x) and st == 'pn_dendrite') else v3
            rec[f'{st}_error'] = v4
            e4.append(v4)
            e3.append(v3)
        rec['stage_error'] = float(np.sqrt(np.mean(np.square(e4))))
        rec['stage_error_three'] = float(np.sqrt(np.mean(np.square(e3))))
        rows.append(rec)
    out = pd.DataFrame(rows)
    if len(out):
        each = out[out.wiring != 'common']
        m = each.drop(columns=['group', 'wiring', 'condition']).mean(skipna=False)
        out = pd.concat([each, pd.DataFrame([dict(group=group, wiring='model', condition=condition, **m.to_dict())]),
                         out[out.wiring == 'common']], ignore_index=True)
    return out


def consistency_gate(readings: Sequence[str] = ()) -> dict:
    """The measure's gate as hal's and viki's archive scoring writes it (experiments/E134_archive_score.py:
    results/.../archive_consistency_<reading>.csv, measure 'error' over the distinct variants inside the bands; the
    planner's reading, chat.md 2026-11-24 11:34): per reading the split-half ρ with APL on and with APL off, the gate
    passing on the default KC-soma reading at GATE_CONSISTENCY_MIN in each condition, 'raw' beside. Read, not
    recomputed; a reading whose file is missing is left out."""
    default = R.DEFAULT_READINGS['kc_soma'].name
    out = {}
    for rd in dict.fromkeys([default, 'raw', *readings]):
        p = out_dir() / f'archive_consistency_{rd}.csv'
        if not p.exists():
            continue
        c = pd.read_csv(p)
        g = c[(c.measure == 'error') & (c.variant_set == 'inside_bands')]
        rho = {str(r.condition): float(r.split_half_rho) for r in g.itertuples()}
        out[rd] = dict(split_half_rho=rho, passes=bool(rho and all(v >= GATE_CONSISTENCY_MIN for v in rho.values())
                                                         and len(rho) >= 2), file=str(p))
    if default in out:
        out['gate'] = out[default]['passes']
    return out


def classify_stages(scores: pd.DataFrame, panel: str = 'megamat', condition: str = 'on',
                    group: str = 'megamat_fitted') -> Dict[str, dict]:
    """Each stage's class on `panel` for the family group `group`, per reading and over its two readings (the stage's
    default and its first alternative): 'inside' (no further from the flies' mean than the most atypical fly on
    INSIDE_MIN_WIRINGS of 8 wirings), 'far' (above FAR_ERROR on as many), 'between' otherwise; 'not readable' when the
    two readings' classes differ."""
    out = {}
    s = scores[(scores.panel == panel) & (scores.condition == condition) & (scores.group == group)]
    for stage in MEGAMAT_STAGES:
        rs = R.readings_of(stage)
        cls = {}
        for rd in rs[:2]:
            x = s[(s.stage == stage) & (s.reading == rd.name) & s.wiring.isin(WIRINGS)]
            if len(x) < len(WIRINGS):
                continue
            inside = int((x.error <= x.flies_loo_max).sum())
            far = int((x.error > FAR_ERROR).sum())
            cls[rd.name] = ('inside' if inside >= INSIDE_MIN_WIRINGS else 'far' if far >= INSIDE_MIN_WIRINGS
                            else 'between')
        if cls:
            vals = list(cls.values())
            out[stage] = dict(by_reading=cls, cls=vals[0] if len(set(vals)) == 1 and len(vals) == 2
                              else ('not readable' if len(set(vals)) > 1 else f'{vals[0]} (one reading)'))
    return out


def step_null_distance(upstream: Sequence[P.Population], downstream: Sequence[P.Population], odors: Sequence[str],
                       edges, q: float = STEP_NULL_Q) -> Dict[str, float]:
    """Per part, the `q` quantile of the RMS distance of a resampled change (the flies of both stages resampled) from
    the flies' measured change: the threshold a model's change of that part must be within to be inside.

    Waits on hal's resampled step distance in kcq/population_fit.py (the planner's reading, chat.md 2026-11-24 11:34),
    not pushed when this was written; it raises until then."""
    raise NotImplementedError("the step class waits on hal's resampled step distance in kcq/population_fit.py "
                              "(chat.md 2026-11-24 11:34: per part, the RMS distance of the resampled flies' change from "
                              "the measured change), not pushed yet")


def classify_steps(step_tab: pd.DataFrame, nulls: Mapping[str, Mapping[str, float]]) -> Dict[str, dict]:
    """Each step's class, as the lock's rules fix it (rules -> step_class; the planner's reading 6, chat.md 2026-11-24
    11:34): a part is inside on a wiring when its RMS distance from the flies' change (`step_table`'s 'part_distance'
    rows) is no larger than `nulls[step][part]` (`step_null_distance`, the 95th percentile of the resampled flies'
    distances); a view is inside on a wiring when every part of it is; the step is 'explained' when every view is
    inside on INSIDE_MIN_WIRINGS of the 8 wirings, 'no better than no step' when the model's downstream error (the
    wirings' mean) is at or above the no-step yardstick, otherwise 'in part'. The per-entry intervals stay as the
    description (the entries inside, of all, per part)."""
    out = {}
    for step, df in step_tab.groupby('step', sort=False):
        pop = df[df.level == 'population']
        dist = df[df.level == 'part_distance'].copy()
        ent = df[df.level == 'entry']
        thr = nulls.get(step) or {}
        dist['inside'] = [bool(part in thr and np.isfinite(d) and d <= thr[part]) for part, d in zip(dist.part,
                                                                                                    dist.distance)]
        dist['view'] = [str(part).split('/')[0] for part in dist.part]
        views = {}
        for v, x in dist.groupby('view', sort=False):
            per_w = x.groupby('wiring').inside.all()
            n_in = int(per_w.sum())
            views[v] = dict(n_wirings_inside=n_in, n_wirings=int(len(per_w)),
                            inside=bool(n_in >= min(INSIDE_MIN_WIRINGS, len(per_w))))
        entries = {part: f'{int(x.inside.sum())} of {len(x)}' for part, x in ent.groupby('part', sort=False)}
        if len(pop) and float(pop.model_downstream_error.iloc[0]) >= float(pop.no_step_error.iloc[0]):
            cls = 'no better than no step'
        elif views and all(v['inside'] for v in views.values()):
            cls = 'explained'
        else:
            cls = 'in part'
        out[step] = dict(cls=cls, views=views, thresholds={k: float(v) for k, v in thr.items()},
                         entries_inside=entries)
    return out


def classify_panels(scores: pd.DataFrame, step2: Optional[pd.DataFrame] = None) -> Dict[str, dict]:
    """Each held-out panel's class at the KC somata (the family as fitted, APL on, the default reading, the wirings'
    mean): 'holds' if its population error is within PANEL_HOLDS_WITHIN of the APL16 panel's (step 2's E90 runs
    against the control flies), else 'panel-specific'."""
    rd = R.DEFAULT_READINGS['kc_soma'].name
    if step2 is None:
        p = out_dir() / 'step2_scores.csv'
        step2 = pd.read_csv(p) if p.exists() else pd.DataFrame()
    base = step2[(step2.get('source') == 'E90:B') & (step2.stage == 'kc_soma') & (step2.panel == 'apl16')
                 & (step2.condition == 'on') & (step2.reading == rd) & (step2.wiring == 'mean')] if len(step2) else []
    if not len(base):
        return {}
    ref = float(base.error.iloc[0])
    out = {}
    for panel in ('megamat', 'validation2') + R.MIXTURE_PANELS:
        inp = 'mixtures' if panel in R.MIXTURE_PANELS else panel
        x = scores[(scores.stage == 'kc_soma') & (scores.panel == panel) & (scores.condition == 'on')
                   & (scores.reading == rd) & (scores.wiring == 'mean') & (scores.group == f'{inp}_fitted')]
        if len(x):
            e = float(x.error.iloc[0])
            out[panel] = dict(error=e, apl16=ref, cls='holds' if abs(e - ref) <= PANEL_HOLDS_WITHIN else 'panel-specific')
    return out


def realism_rows(run: R.SavedRun, fit_params: Optional[Mapping[str, float]] = None) -> dict:
    """A run's realism, reported only: overall response rates with APL on and off and their ratio, the largest spike
    count per KC and odor in each condition, the threshold beside the fit's own (held runs)."""
    row = dict(tag=run.tag, model=run.model, wiring=run.wiring, panel=run.input_panel)
    row.update(_counts_row(run.conditions))
    thr = (run.info.get('params') or {}).get('fixed_thr') if isinstance(run.info, Mapping) else None
    row['fixed_thr'] = np.nan if thr is None else float(thr)
    if fit_params is not None:
        row['fit_fixed_thr'] = float(fit_params['fixed_thr'])
    return row


def apl16_bands() -> pd.DataFrame:
    """The family's four realism bands on the APL16 panel, from results/apl_effect/by_variant.csv (E72's fits of the
    eight wirings, E83's `members`)."""
    from experiments import E83_cross_wiring_fit as E83
    bv = pd.read_csv(H.shared_dir() / 'results' / 'apl_effect' / 'by_variant.csv')
    rows = []
    for w, (folder, tag) in E83.members().items():
        x = bv[(bv.folder == folder) & (bv.tag == tag)]
        if len(x):
            r = x.iloc[0]
            rows.append(dict(wiring=w, folder=folder, tag=tag, **{k: r.get(k) for k in (
                'rate_on', 'rate_off', 'spearman_rate_on_vs_data', 'opk_on_tv', 'opk_off_tv', 'realism_flag', 'realism')}))
    return pd.DataFrame(rows)


def part_score() -> dict:
    """Step 3's runs (and E107's validation2 runs of the family) read and scored as step 2's (default readings and
    alternatives); the stage table, the step table, the held-out panels, the classes, realism, score_table.md and
    provenance_score.json."""
    lock = require_lock()
    t0 = time.time()
    flies = R.FlyCache()
    scores, steps, step_classes = [], [], {}
    groups = score_groups()
    for g, runs in groups.items():
        for stage, rd in _reading_sets([s for s in R.STAGES if s in runs[0].stages()], None):
            fam = R.read_family(runs, stages=[stage], readings={stage: rd}, lock=lock, flies=flies)
            sc = score_family(fam, rd.name, flies)
            if len(sc):
                sc.insert(0, 'group', g)
                scores.append(sc)
            if rd is R.DEFAULT_READINGS[stage]:
                save_populations(fam, 'score', g)
        if g in ('megamat_fitted', 'megamat_held'):
            fam = R.read_family(runs, lock=lock, flies=flies)
            st, inputs = step_table(fam, flies)
            if len(st):
                st.insert(0, 'group', g)
                steps.append(st)
                try:
                    step_classes[g] = classify_steps(st, step_nulls(inputs))
                except NotImplementedError as e:
                    step_classes[g] = dict(error=str(e))
        print(f'score {g}: {time.time() - t0:.0f} s', flush=True)
    sc = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    stp = pd.concat(steps, ignore_index=True) if steps else pd.DataFrame()
    sc.to_csv(out_dir() / 'score_scores.csv', index=False)
    stp.to_csv(out_dir() / 'score_step_table.csv', index=False)
    stage_table(sc).to_csv(out_dir() / 'score_stage_table.csv', index=False)
    se = pd.concat([megamat_stage_error(sc, g) for g in dict.fromkeys(sc.group)], ignore_index=True) if len(sc) \
        else pd.DataFrame()
    se.to_csv(out_dir() / 'score_stage_error.csv', index=False)
    gate = consistency_gate()
    gw = out_dir() / 'step2_gate_wirings.csv'
    if gw.exists():
        g2 = pd.read_csv(gw)
        g2 = g2[g2.apply(lambda r: r.reading == R.DEFAULT_READINGS[r.stage].name, axis=1)] if len(g2) else g2
        gate['wirings'] = {f'{r.source}|{r.stage}|{r.panel}|{r.condition}': dict(mean_error=float(r.mean_error),
                           passes=bool(r.passes)) for r in g2.drop_duplicates(['source', 'stage', 'panel',
                                                                                'condition']).itertuples()}
    classes = dict(stages={g: classify_stages(sc, group=g) for g in ('megamat_fitted', 'megamat_held')} if len(sc) else {},
                   steps=step_classes, panels=classify_panels(sc) if len(sc) else {}, gate=gate,
                   megamat_stage_error={g: float(x[x.wiring == 'model'].stage_error.iloc[0])
                                        for g, x in se.groupby('group')} if len(se) else {})
    _write_json(out_dir() / 'score_classes.json', classes)
    rows = [realism_rows(r) for runs in groups.values() for r in runs if r.source == 'E134:step3']
    pd.DataFrame(rows).to_csv(out_dir() / 'score_realism.csv', index=False)
    try:
        apl16_bands().to_csv(out_dir() / 'score_apl16_bands.csv', index=False)
    except FileNotFoundError:
        pass
    write_table(sc, stp, classes)
    provenance('score', lock=lock, seconds=round(time.time() - t0, 1), groups={g: [r.tag for r in rs]
                                                                                for g, rs in groups.items()})
    return classes


def write_table(scores: pd.DataFrame, steps: pd.DataFrame, classes: Mapping) -> None:
    """score_table.md: the stage table (default readings, the family's wirings' mean and common error beside the flies'
    floor and most atypical fly), the step classes, the panel classes."""
    L = ['# E134: the stage table (default readings)', '']
    if len(scores):
        d = scores[scores.apply(lambda r: r.reading == R.DEFAULT_READINGS[r.stage].name, axis=1)]
        d = d[d.wiring.isin(['mean', 'common'])]
        L += ['| group | stage | panel | condition | wirings | error | strength | sharing | breadth | flies floor | '
              'atypical fly |', '|---|---|---|---|---|---:|---:|---:|---:|---:|---:|']
        for r in d.itertuples():
            L.append(f'| {r.group} | {r.stage} | {r.panel} | {r.condition} | {r.wiring} | {r.error:.3f} | '
                     f"{getattr(r, 'view_strength', np.nan):.3f} | {getattr(r, 'view_sharing', np.nan):.3f} | "
                     f"{getattr(r, 'view_breadth', np.nan):.3f} | {r.flies_floor:.3f} | {r.flies_loo_max:.3f} |")
    L += ['', '## Classes', '', '```', json.dumps(classes, indent=1, default=str), '```', '']
    (out_dir() / 'score_table.md').write_text('\n'.join(L) + '\n')


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, the lock's sha256, the readings
    (`stage_readings.reading_definitions`), host, time and `extra` (the lock itself is not copied)."""
    lock = extra.pop('lock', None)
    _write_json(out_dir() / f'provenance_{part}.json', dict(
        part=part, script=SCRIPT, **H.code_commits(), lock_sha256=(lock or {}).get('sha256') if isinstance(lock, Mapping)
        else None, lock_readings_agree=R.lock_readings_agree(lock) if isinstance(lock, Mapping) else None,
        readings=R.reading_definitions(), host=_host(), time=_now(), **extra))


# ---------------------------------------------------------------------------
# lanes and CLI

def print_commands() -> None:
    """The lanes' command lines (LANES): every lane runs the whole list and skips what another has claimed or done;
    the reversed lanes take it from the other end."""
    pre = 'OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u ' + SCRIPT
    for host, (t, lanes) in LANES.items():
        for name, rev in lanes:
            print(f"{host}: {pre.format(t=t)} --part step3 --lane {name}{' --reverse' if rev else ''}")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step2', 'step3', 'score', 'commands'))
    ap.add_argument('--models', nargs='+', default=None,
                    help='step3: the wirings and beside models of the lane (default all)')
    ap.add_argument('--panels', nargs='+', default=None, choices=list(R.INPUT_PANELS),
                    help='step3: the input panels (default all)')
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
