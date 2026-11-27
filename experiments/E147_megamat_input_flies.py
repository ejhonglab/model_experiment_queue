#!/usr/bin/env python3
"""E147, the megamat stage error read on model flies with each ORN fly's own input: whether the averaged input is a
fair stand-in at every megamat stage, what the model's own flies differ by, and what the classes' spreads leave out
(experiment_priorities.md, E147; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-27 09:01;
moved in by viki under the idle-host rule 2026-11-27 09:19; atlas: the lock, this module, step 0 and the score; hal: the
reads of steps 1 and 2, in experiments/E147_reads.py, hal's alone).

Nothing here simulates, fits, searches or weights. The reads (steps 1 and 2) are hal's runner's
(`./venv/bin/python -u experiments/E147_reads.py --part reads --step 1|2`); this module gives it the input flies, the
points, the tags, the folders and the lock, and scores what it saves.

Definitions (the entry's, fixed in experiments/E147_lock.json before any statistic is read):
- an input fly: one of the nine megamat ORN flies of olfdb dataset 1 as E134's copy has them (the lock's
  orn / megamat / wild_type flies, '2023-04-22/2' ... '2023-06-22/1'), labelled f0 ... f8 in that (sorted) order,
  which is E85's; its input is the lab's input plus the fly's reproducible deviation (E85's 'repro': the reliability
  r_g of each glomerulus, pooled over E85's 12 odors, times E85's common factor, times the fly's trial mean minus the
  mean over the flies that have the entry, on E85's scaling to spike deltas), for the odors and glomeruli measured in
  that fly, and the lab's value elsewhere (`input_fly`). On megamat (the lab's megamat input, E134's
  `panel_input('megamat')`, 17 odors) every odor the fly was given; on the APL16 input (`H.apl16_orn_deltas()`) the 12
  megamat odors, with paa, euc, -aPine and B-myr kept at the lab's input (their ORN recordings are other flies'). The
  mean over the nine flies of every entry is the lab's input exactly (aa and va were given to 7 of the 9);
- a model fly: one wiring with one input fly, at the setting's parameters, both response rates re-held at the panel's
  operating point as E136's `read_at` holds a read (E134's `held_targets` on megamat, the wiring's own fit's counts on
  APL16), the breadth the lab input's, the wiring's seed for pooling and noise kept; read exactly as a wiring is read
  (E134's default readings, E136's `Flies` and `StageLoss`);
- the population's error: the error of the mean view over the model flies; a model fly's error: its own; both against
  all the stage's flies, E136's held-out half beside;
- the terms: the SD over input flies of the mean over wirings (the input-fly term), the SD over wirings of the mean
  over input flies (the wiring term), and the remainder (`terms`);
- a model fly against the others: its view against the mean view of the model flies with another wiring and another
  input fly, in the stage's one-fly units (`fly_against_others`);
- the spread over model flies: the SD over all model flies of the quantity whose SD over the eight wirings the classes
  use (`model_fly_spreads`).

Points and tags (`Point`, `point_tag`, `parse_tag`, `planned_points`): '<setting>_<panel>_<inputfly>_<wiring>', e.g.
'family_megamat_f3_P4', with setting one of SETTINGS ('family': the unrefitted family at each wiring's held megamat
point, E136's FAMILY; 'A1': E136's A1 at its classed point, `E136.fitted_z('A1')`; 'familyfree' and 'A1free': the free
reading, each model fly at its wiring's held parameters with the input fly's input, nothing re-held), panel 'megamat' or
'apl16', input fly f0 ... f8, wiring one of E136's WIRINGS. Step 1: 'family' on megamat, the nine input flies on all
eight wirings, the T wirings first (72). Step 2, on T: 'A1' on megamat (36), then 'A1' and 'family' on APL16 (36 each).
The free reading ('familyfree', 'A1free', one run per model fly of steps 1 and 2) is beside, classing nothing, and is
the runner's too.

What hal's runner writes, which `--part score` reads:
- each read's run in E134's step-3 layout (`kcq.stage_readings.save_step3_run`, every stage saved, the units' metadata,
  the input as `<tag>_orn.parquet`) under `run_dir(step, setting, wiring, fly)` =
  outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/<fly>/ (local; `run_dirs` also gives the archive's
  data/model_outputs/E147_megamat_input_flies/... equivalent); a read whose `<tag>_run.json` is there is done;
- its row in results/E147_megamat_input_flies/reads_<lane>.csv (any name reads_*.csv; the last row per tag counts,
  `reads_table`), with the columns of READ_COLUMNS: tag, step, setting, panel, fly, wiring; the held parameters
  (PARAM_KEYS, absolute, and off_<key> of the APL-off run); thr_factor, apl_factor (against the wiring's held point the
  hold started from); n_on, n_off, target_on, target_off, frac_on_target, frac_off_target; at_bound_thr, at_bound_apl;
  held (True when n_on == target_on and n_off == target_off and neither factor is at a bound); rate_on, rate_off,
  max_count_on, max_count_off; input_sha256 (E77's `frame_sha256` of the input, which `--part score` checks against
  `input_fly`); hold_n_sims, seconds, host, lane, time, threads and the commits (`H.code_commits()`). A read that
  failed to hold writes its row with held False and its run as held so far, so that the score counts it unheld;
- E134's reading of the run beside its counts is optional: `--part score` reads every run itself, through E136's
  `StageLoss` and `Flies` (E134's default readings, seed index the wiring's).

Parts (CLI `--part`):
- lock (atlas): experiments/E147_lock.json and results/.../lock.json; commit the copy before any statistic is read;
- step0 (atlas, no simulation): the identities (E134's 0.918 and 0.894; E85's part A on the 12 shared odors; E140's aa
  vectors), the table of megamat's input flies on its 17 odors, E85's 72 saved model flies and its bounds on S scored
  at the APL16 somata as E134's archive scoring scores a run, every class of E136, E139 and E141 recomputed with the
  APL16 spreads over model flies, and prediction.json (a prediction per megamat stage, timestamped, before any run);
- step1, step2: hal's (experiments/E147_reads.py --part reads --step 1|2); here they only say so;
- score (atlas, no simulation; run once hal's reads are in): the classes per megamat stage, for the megamat stage error
  and at the APL16 somata; A1's megamat fall and its APL16 changes model fly by model fly against the family on the
  same wiring and input fly; the megamat spreads over model flies and every class of E136, E139 and E141 recomputed
  with them, as a description; the realism bands per model fly and per population; the beside items the reads allow.

Outputs: <shared>/results/E147_megamat_input_flies/ (compact): lock.json, step0_*.csv, classes_step0.json,
prediction.json, score_*.csv, classes.json, tables.md, provenance_<part>.json; hal's reads_<lane>.csv,
provenance_reads_<lane>.json, claims/. outputs/E147_megamat_input_flies/ (local, large): step1/, step2/ (the reads),
step0/ (intermediates).

Concentrations: megamat is at 10^-3 on both sides. Every APL16 number involves paa, euc and -aPine, whose model input
is the ORN data within 1 log10 unit of the flies' concentrations; those three and B-myr keep the lab's averaged input in
every model fly, so the input-fly term at the APL16 somata is a lower bound.

Holds (CLAUDE.md): Pratyush's hemibrain only (the family's wirings, 'first' claws, pinned); no olfsysm change; nothing
searched, fitted or weighted; validation2 is not read; E145's cells are not read; Remy's processed data without the
bouton fly 2022-10-21 fly02 and without her older megamat flies (E134's lock).

Usage (from the repository root):
    ./venv/bin/python -u experiments/E147_megamat_input_flies.py --part lock     # then commit experiments/E147_lock.json
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E147_megamat_input_flies.py --part step0
    OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2 nice -n 10 ./venv/bin/python -u experiments/E147_megamat_input_flies.py --part score
"""
from __future__ import annotations

import argparse
import json
import re
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E147_megamat_input_flies'
SCRIPT = 'experiments/E147_megamat_input_flies.py'
RUNNER = 'experiments/E147_reads.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E147_lock.json'
E134 = 'E134_population_baseline'
E85 = 'E85_orn_fly_population'
E136 = 'E136_stage_population_fit'

S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PANELS: Tuple[str, ...] = ('megamat', 'apl16')
SETTINGS: Tuple[str, ...] = ('family', 'A1', 'familyfree', 'A1free')
HELD_SETTINGS: Tuple[str, ...] = ('family', 'A1')
N_INPUT_FLIES = 9
FLY_LABELS: Tuple[str, ...] = tuple(f'f{i}' for i in range(N_INPUT_FLIES))
# the APL16 odors whose ORN recordings are other flies': they keep the lab's input in every model fly
APL16_KEPT_AT_LAB: Tuple[str, ...] = ('paa', 'euc', '-aPine', 'B-myr')
DEVIATION = 'repro'
PARAM_KEYS: Tuple[str, ...] = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale')
READ_COLUMNS: Tuple[str, ...] = (
    'tag', 'step', 'setting', 'panel', 'fly', 'wiring', *PARAM_KEYS, *(f'off_{k}' for k in PARAM_KEYS), 'thr_factor',
    'apl_factor', 'n_on', 'n_off', 'target_on', 'target_off', 'frac_on_target', 'frac_off_target', 'at_bound_thr',
    'at_bound_apl', 'held', 'rate_on', 'rate_off', 'max_count_on', 'max_count_off', 'input_sha256', 'hold_n_sims',
    'seconds', 'host', 'lane', 'time', 'threads', 'queue_commit', 'queue_dirty', 'pyolfsysm_commit', 'pyolfsysm_dirty')

Point = namedtuple('Point', 'step setting panel fly wiring')


# ---------------------------------------------------------------------------
# folders and the lock

def results_dir() -> Path:
    """<shared>/results/E147_megamat_input_flies/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    """outputs/E147_megamat_input_flies/[<part>/] in this repository (created)."""
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E147_megamat_input_flies/ (not created)."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def lock_path() -> Path:
    """The lock's shared copy, results/E147_megamat_input_flies/lock.json; the committed copy is LOCK_COPY."""
    return results_dir() / LOCK_NAME


def require_lock() -> dict:
    """The committed lock (LOCK_COPY tracked with no uncommitted change, and equal to the shared copy): PermissionError
    otherwise. Every part but 'lock' calls it before reading any statistic."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the input flies

def input_flies(panel: str = 'megamat') -> List[str]:
    """The input flies' labels, f0 ... f8 (E134's copy's orn / megamat / wild_type flies in its sorted order, which is
    E85's), the same for both panels."""
    if panel not in PANELS:
        raise ValueError(f'{panel=}')
    return list(FLY_LABELS)


def fly_key(fly: str) -> str:
    """An input fly's E134 key ('2023-04-22/2', ...) from its label ('f0', ...) or the key itself."""
    raise NotImplementedError


def lab_input(panel: str) -> pd.DataFrame:
    """The lab's input of a panel, glomeruli x odors, exactly as E134's held runs had it: E134's
    `panel_input('megamat')` (E107's, 17 odors, columns as loaded) or `H.apl16_orn_deltas()` (16 odors)."""
    raise NotImplementedError


def deviations(panel: str, variant: str = DEVIATION) -> np.ndarray:
    """flies x odors x glomeruli deviations to add to `lab_input(panel)` (its odor and glomerulus order; 0 where a fly
    did not measure an entry, and on the APL16 odors of APL16_KEPT_AT_LAB and every non-megamat odor)."""
    raise NotImplementedError


def input_fly(fly: str, panel: str = 'megamat', variant: str = DEVIATION) -> pd.DataFrame:
    """Input fly `fly`'s model input on `panel` ('megamat' or 'apl16'): `lab_input(panel)` (same index, columns and
    dtypes) plus the fly's reproducible deviation where measured, the lab's value elsewhere; on APL16 paa, euc, -aPine
    and B-myr at the lab's input. `fly` is a label (f0 ... f8) or an E134 key. Pass it to E134's `build_family(...,
    orn=...)` and E77's `hold_rates` in place of the panel's input."""
    raise NotImplementedError


def zero_deviation_input(panel: str = 'megamat') -> pd.DataFrame:
    """`input_fly`'s construction with every deviation 0: equal to `lab_input(panel)` (checked), for the runner's
    identity (an input fly with no deviation reproduces E134's held megamat run with 0 differing counts)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# points, tags and the reads' folders

def point_tag(setting: str, panel: str, fly: str, wiring: str) -> str:
    """'<setting>_<panel>_<inputfly>_<wiring>' (ValueError for an unknown part)."""
    if setting not in SETTINGS or panel not in PANELS or fly not in FLY_LABELS or wiring not in WIRINGS:
        raise ValueError(f'{setting=} {panel=} {fly=} {wiring=}')
    return f'{setting}_{panel}_{fly}_{wiring}'


def parse_tag(tag: str) -> Tuple[str, str, str, str]:
    """(setting, panel, fly, wiring) of a tag (`point_tag`'s inverse)."""
    m = re.fullmatch(r'([A-Za-z0-9]+)_(megamat|apl16)_(f\d)_([A-Za-z0-9]+)', tag)
    if not m:
        raise ValueError(tag)
    setting, panel, fly, wiring = m.groups()
    point_tag(setting, panel, fly, wiring)
    return setting, panel, fly, wiring


def planned_points(step: int, free: bool = False) -> List[Point]:
    """The reads of a step in the order a forward lane takes them: step 1 'family' on megamat, the T wirings first,
    then S, fly by fly within a wiring (72); step 2 on T 'A1' on megamat (36), then 'A1' on APL16 (36), then 'family'
    on APL16 (36). `free`: the same points under the free reading ('familyfree', 'A1free')."""
    fl = list(FLY_LABELS)
    if step == 1:
        pts = [Point(1, 'family', 'megamat', f, w) for w in T_WIRINGS + S_WIRINGS for f in fl]
    elif step == 2:
        pts = [Point(2, s, p, f, w) for s, p in (('A1', 'megamat'), ('A1', 'apl16'), ('family', 'apl16'))
               for w in T_WIRINGS for f in fl]
    else:
        raise ValueError(f'{step=}')
    if free:
        pts = [pt._replace(setting=f'{pt.setting}free') for pt in pts]
    return pts


def step_of(setting: str, panel: str) -> int:
    """The step a point belongs to: 1 for the family on megamat, 2 otherwise (the free reading as its held twin)."""
    base = setting[:-4] if setting.endswith('free') else setting
    return 1 if (base, panel) == ('family', 'megamat') else 2


def run_dir(step: int, setting: str, wiring: str, fly: str) -> Path:
    """Where a read's run is saved locally: outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/<fly>/
    (created)."""
    return local_dir(f'step{int(step)}/{setting}/{wiring}/{fly}')


def run_dirs(step: int, setting: str, wiring: str, fly: str) -> List[Path]:
    """`run_dir` and its archive twin (data/model_outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/
    <fly>/), in that order, for `kcq.stage_readings.save_step3_run` and `load_step3_run`."""
    return [run_dir(step, setting, wiring, fly), archive_dir() / f'step{int(step)}' / setting / wiring / fly]


def reads_table() -> pd.DataFrame:
    """Every lane's results/.../reads_*.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(results_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=list(READ_COLUMNS))
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def load_point(tag: str):
    """A saved read (`kcq.stage_readings.SavedRun`; locally, else the archive); None if absent."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

def part_lock() -> dict:
    """The lock: every file, fly set, odor set, seed, class threshold, the T-first order and the choices; written to
    results/.../lock.json and LOCK_COPY (commit the copy before any other part)."""
    raise NotImplementedError


def part_step0() -> dict:
    """Step 0 (the module docstring; no simulation). Stops (RuntimeError) if an identity fails."""
    raise NotImplementedError


def part_score() -> dict:
    """The score of hal's reads (the module docstring; no simulation)."""
    raise NotImplementedError


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'step0', 'step1', 'step2', 'score'))
    a = ap.parse_args(argv)
    if a.part in ('step1', 'step2'):
        raise SystemExit(f'{a.part} is the runner\'s: ./venv/bin/python -u {RUNNER} --part reads --step {a.part[-1]}')
    {'lock': part_lock, 'step0': part_step0, 'score': part_score}[a.part]()


if __name__ == '__main__':
    main()
