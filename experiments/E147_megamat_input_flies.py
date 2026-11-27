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
  held (the runner's own reading of its hold); rate_on, rate_off,
  max_count_on, max_count_off; input_sha256 (E77's `frame_sha256` of the input, which `--part score` checks against
  `input_fly`); hold_n_sims, seconds, host, lane, time, threads and the commits (`H.code_commits()`). The score decides
  "held" itself (`is_held`: both counts within HOLD_TOL, 0.5 %, of their targets and no factor at a bound), so the
  row's own 'held' is informative only. A read that failed to hold writes its row and its run as held so far, so that
  the score counts it unheld;
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

# the classes (the entry's, fixed beforehand)
TERM_SMALL, TERM_PART = 0.5, 1.0          # input-fly term / wiring term: "small" below, "a part" below, else "as large"
FLIES_LOW, FLIES_HIGH = 0.75, 1.25        # a model fly against the others: "too little" below, "as the flies", "too much"
IMPROVE_SE = 2.0                          # A1's fall: beyond this many standard errors
CONTROL_FALL = 0.040                      # and beyond the control fit's fall (E136 classes.json control_fall, 0.0403)
UNHELD_MAX_SHARE = 1.0 / 9.0              # no class for a population with more than this share of its model flies unheld
HOLD_TOL = 0.005                          # held: both counts within this share of their targets (E136's HOLD_TOL)
MEASURE = 'strength'                      # E136's adopted measure (MEASURES['strength'])
E136_MEASURE_DIR = 'measure_strength'
APL16_KEYS: Tuple[str, ...] = ('apl16/apl16/on', 'apl16/apl16/off')
# the class tables recomputed with the spreads over model flies (E143's CLASS_TABLES): (experiment, path under results/)
CLASS_TABLES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    (E136, (E136_MEASURE_DIR, 'classes.json')), ('E139_normalization_ladder', ('classes.json',)),
    ('E141_kckc_rise_forms', ('classes.json',)))
E134_IDENTITY = {'model': 0.918247443617637, 'common': 0.893708030399992}    # writeup_stage_error.csv, megamat_held
IDENTITY_TOL = 1e-9
# E85's saved model flies (runs.csv stages) and what each is in this experiment's step 0
E85_GROUPS: Dict[str, str] = {'B1': "'repro' held, all eight wirings (B1S + B1T): the 72 model flies",
                              'B2': "'measured' (s = 1) held, S: the upper bound on the deviations' size",
                              'B3': "'gainfree' held, S", 'B4': "'permuted' held, S: the deviations without their "
                              "glomerulus structure", 'B5': "'repro' free (the fits' parameters, nothing re-held)"}
CHOICES: Tuple[str, ...] = (
    "the deviation is E85's 'repro' exactly: r_g pooled over E85's 12 odors (not megamat's 17) times E85's common "
    "factor (part_a.json, 1.0), on E85's slope to spike deltas (99.77, the APL16 input on the flies' mean over the 12 "
    "shared odors, recomputed from E134's copy), so that a fly's deviation on a shared odor is the same on both panels "
    "and its APL16 input is E85's saved 'repro' input; the 17-odor slope (101.1, against the megamat input) and the "
    "17-odor r_g are reported beside, unused",
    "megamat's five odors outside APL16 (aa, va, B-cit, Lin, benz) get the deviation by the same rule; aa and va were "
    "not given to f0 and f1 (2023-04-22 flies 2 and 3), which keep the lab's value there",
    "a model fly's error and the population's are read against all the stage's flies (E134's reading: the curve, noise "
    "and criterion of all the stage's flies, seed index the wiring's position in E136's WIRINGS, the same for every input "
    "fly on that wiring); E136's held-out half beside (E136's Flies('held'): the half's curve and noise, all flies' "
    "units), as E136 reads a read",
    "the measure is strength by odor (E136's MEASURES['strength']); the megamat stage error is the root mean square of "
    "the four stages' strength errors (E134's writeup 'stage_error_strength'); every view's error beside",
    "the terms of a quantity E[wiring, input fly] (a model fly's error): input-fly term the SD over input flies (ddof 1) "
    "of the mean over wirings, wiring term the SD over wirings (ddof 1) of the mean over input flies, remainder the root "
    "mean square of the two-way residual (E - row mean - column mean + grand mean) scaled by sqrt(n_w n_f / ((n_w - 1) "
    "(n_f - 1))), on a complete wirings x input flies table only",
    "a model fly against the others: its views scored against the mean views of the model flies with another wiring "
    "and another input fly (P.score_views on the stage's reference with its mean replaced, so in the recorded flies' "
    "one-fly units); a stage's value is the mean over the model flies (min and max beside), the megamat stage error's "
    "the mean over model flies of the root mean square over the four stages",
    "'a fair stand-in': |the population's error - the error of the wirings' mean view at the lab's input| <= the spread "
    "as scored (E136's classes.json spreads: 'stage' for the megamat stage error, each stage's, 'apl16/apl16/on' and "
    "'/off' at the APL16 somata); 'nearer' / 'further' by its sign otherwise; the lab input's value is E134's held "
    "megamat runs on the same wirings (step 1: all eight; T only where the population is T's) and E90's part B on APL16",
    "the spread over model flies: the SD (ddof 1) over the model flies of a model fly's error, the quantity whose SD over "
    "the eight wirings (ddof 1) E136's spreads are; step 0's APL16 spreads from E85's 72 'repro' held model flies "
    "(breadth recomputed per input), the score's from step 2's family APL16 reads (breadth the lab input's) beside them",
    "the classes recomputed with other spreads: every comparison of E136's (measure_strength), E139's and E141's class "
    "tables keeps its per-wiring changes; its spread is replaced where a spread over model flies exists, its standard "
    "error recomputed by E136's compare (the larger of half the spread and the SD of the changes over sqrt(n)), and the "
    "class read by E143's `reclass` ladder (E141's comparisons with A1 and its added classes by E143's `against_a1` "
    "and `e141_entry`); with the spreads as scored it must give every scored class (step 0's check); the guard's "
    "comparisons (every view) take the all-view spreads over model flies",
    "A1's fall: per model fly (T wiring, input fly) A1's megamat stage error minus the family's on the same wiring and "
    "input fly; the mean over the 36; the standard error the larger of the SD over input flies of the per-fly means "
    "(over the four wirings) over sqrt(9) and the SD over wirings of the per-wiring means over sqrt(4); it stands if "
    "-mean > IMPROVE_SE x se and -mean > CONTROL_FALL; classed against E136's held-out half, where E136's class and "
    "the control's 0.040 are read, with all the flies beside",
    "A1 at the APL16 somata: per input fly the mean over the T wirings of A1's strength error minus the family's (step "
    "2's reads, both breadths the lab input's); the share of the nine on which it is above the spread as scored "
    "(0.051 with APL on, 0.090 off), and below minus it, beside",
    "the megamat spreads over model flies, for the classes recomputed: the SD over step 1's 72 model flies of each "
    "model fly's error against E136's held-out half (the quantity of E136's spreads: 'stage', 'downstream', each "
    "stage on strength, the guard's on every view), against all the flies beside; the APL16 spreads from step 0 (E85's "
    "72) and from step 2's family reads (36, on T); the classes recomputed with (a) the megamat spreads, the APL16 ones "
    "as scored, (b) both over model flies (E85's APL16), (c) both, step 2's APL16",
    "a fair stand-in at A1's population: against A1's own reads at the lab's input on T (E136's A1 megamat reads); the "
    "family's T-only population beside its eight-wiring one; at the APL16 somata against the family's own fits on the "
    "same wirings (E90's part B counts, E83's saved fits) and A1's E136 APL16 reads",
    "the beside items: each odor's strength entries by input fly at every megamat stage (the mean view over the "
    "wirings of each input fly's model flies minus the flies' mean, in one-fly units), as leads; the views' errors "
    "(sharing, breadth) of every population; the somata without aa's entry (E136's `stage_without_odor` rule: aa's "
    "term out of the squared error, nothing renormalized); each input fly's megamat input at the ORN stage in its "
    "one-fly units (step 0, from the inputs alone); step 2's family APL16 counts beside E85's saved ones on the same "
    "wiring and input fly",
    "unheld: a read whose APL-on or APL-off response count misses its target by more than HOLD_TOL (0.5 %, E136's "
    "hold tolerance) of the target, or whose threshold or APL factor is at a bound (`is_held`, from the row's n_on, "
    "n_off, target_on, target_off, at_bound_thr and at_bound_apl); unheld model flies are left out of every value and "
    "counted; a population with more than a ninth unheld has no class; the free reading is not held and has no unheld",
    "the prediction (prediction.json) is a rule written before any read: it is derived from step 0's APL16 numbers and "
    "the stage's place (the PN stages see the input before any KC parameter), and says per megamat stage the expected "
    "class of each of the three classed questions",
)


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


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (LOCK_COPY tracked with no uncommitted change, and the same bytes as the shared copy):
    PermissionError otherwise. Every part but 'lock' calls it before reading any statistic; its files are checked
    (`check_files`)."""
    from kcq import stage_readings as R
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it)')
    if not R._committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = lock_path() if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = sha256_file(copy)
    return lock


def sha256_file(path: Path) -> str:
    """The sha256 of a file's bytes."""
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def resolve(path: str) -> Path:
    """A path the lock records: '<shared>/...' under the shared directory, '<repo>/...' under this repository."""
    if path.startswith('<shared>/'):
        return H.shared_dir() / path[len('<shared>/'):]
    if path.startswith('<repo>/'):
        return H.HERE / path[len('<repo>/'):]
    return Path(path)


def file_record(path: Path) -> dict:
    """{'path' ('<shared>/...' or '<repo>/...' where it is under either), 'bytes', 'sha256'} of a file."""
    p = Path(path).resolve()
    rel = str(p)
    for root, tag in ((H.shared_dir().resolve(), '<shared>'), (H.HERE.resolve(), '<repo>')):
        try:
            rel = f'{tag}/{p.relative_to(root)}'
            break
        except ValueError:
            continue
    return dict(path=rel, bytes=p.stat().st_size, sha256=sha256_file(p))


def check_files(lock: Mapping) -> None:
    """Every file of the lock's 'files' has its recorded sha256; ValueError otherwise."""
    for name, e in lock['files'].items():
        got = sha256_file(resolve(e['path']))
        if got != e['sha256']:
            raise ValueError(f"{name} ({e['path']}): sha256 {got[:12]} differs from the lock's {e['sha256'][:12]}")


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
    keys = orn_flies()['keys']
    if fly in FLY_LABELS:
        return keys[FLY_LABELS.index(fly)]
    if fly in keys:
        return fly
    raise ValueError(f'{fly!r} is not an input fly')


def fly_label(fly: str) -> str:
    """An input fly's label (f0 ... f8) from its label or its E134 key."""
    return FLY_LABELS[orn_flies()['keys'].index(fly_key(fly))]


def lab_input(panel: str) -> pd.DataFrame:
    """The lab's input of a panel, glomeruli x odors, exactly as E134's held runs had it: E134's
    `panel_input('megamat')` (E107's, 17 odors, columns as loaded) or `H.apl16_orn_deltas()` (16 odors)."""
    from experiments import E134_population_baseline as B
    if panel not in PANELS:
        raise ValueError(f'{panel=}')
    return B.panel_input(panel)


def olfdb_copy() -> Path:
    """E134's copy of olfdb dataset 1 (data/stage_populations/sources/olfdb/ds1_long.csv in the shared directory)."""
    return H.shared_dir() / 'data' / 'stage_populations' / 'sources' / 'olfdb' / 'ds1_long.csv'


def e85_part_a() -> dict:
    """E85's part_a.json (its scaling, reliabilities and common factor)."""
    return json.loads((H.results_dir(E85) / 'part_a.json').read_text())


_ORN: Dict[str, object] = {}


def orn_flies(df: Optional[pd.DataFrame] = None) -> dict:
    """The nine megamat ORN flies on the lab input's 36 glomeruli and megamat's 17 odors, by E85's map, from E134's
    copy of olfdb dataset 1 (cached; `df` a long table in its place, for tests, not cached). Each fly's dF/F is scaled
    'to-avg-max' (its largest |trial mean| over all its panels and ROIs brought to pyolfsysm's avg_flymax) and by one
    slope to spike deltas: E85's, the least-squares slope through the origin of the lab's APL16 input on the flies' mean
    over the 12 shared odors (recomputed here and checked against E85's part_a.json by step 0's identity), the same
    for every odor of both panels. Returns 'keys' (E134's 'date/fly' keys, sorted, which is E85's order), 'labels',
    'odors' (megamat's 17 in the lab's megamat input's column order), 'odors12' (H.MEGAMAT12), 'gloms', 'x' (spike-delta
    trials, flies x odors x glomeruli x trials, NaN where not measured), 'slope', 'slope_r2', 'slope_r2_origin',
    'slope17' and 'slope17_r2_origin' (the same regression on megamat's 17 odors against the megamat input, beside),
    'fly_max', 'avg_flymax', 'repeats'."""
    planted = df is not None
    if not planted and _ORN:
        return _ORN
    from experiments import E85_orn_fly_population as E85m
    if not planted:
        df = pd.read_csv(olfdb_copy())
    lab16 = H.apl16_orn_deltas()
    meg = lab_input('megamat')
    gloms = [str(g) for g in lab16.index]
    odors = [H.odor_abbrev(c) for c in meg.columns]
    keys = sorted({(str(d), int(f)) for d, f in zip(df.date, df.fly_num)})
    kstr = [f'{d}/{f}' for d, f in keys]
    df = df.assign(fly=[keys.index((str(d), int(f))) for d, f in zip(df.date, df.fly_num)])
    tm = df.groupby(['fly', 'panel', 'odor1', 'roi']).value.mean()
    fly_max = tm.abs().groupby(level='fly').max().reindex(range(len(keys))).to_numpy(dtype=float)
    want = {f'{o} @ -3': o for o in odors}
    d = df[(df.panel == 'megamat') & df.odor1.isin(want) & df.roi.astype(str).isin(gloms)]
    if set(d.odor1) != set(want):
        raise RuntimeError(f'olfdb dataset 1: missing {sorted(set(want) - set(d.odor1))}')
    reps = sorted(d['repeat'].unique())
    x = np.full((len(keys), len(odors), len(gloms), len(reps)), np.nan)
    x[d.fly.to_numpy(), np.array([odors.index(want[o]) for o in d.odor1]),
      np.array([gloms.index(str(g)) for g in d.roi]), np.array([reps.index(r) for r in d['repeat']])] = \
        d.value.to_numpy(dtype=float)
    A = E85m.avg_flymax()
    x = x / fly_max[:, None, None, None] * A
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        M = np.nanmean(np.nanmean(x, axis=3), axis=0)           # odors x glomeruli

    def through_origin(Mx, lab):
        ok = ~np.isnan(Mx)
        a, b = Mx[ok], lab[ok]
        sl = float((a @ b) / (a @ a))
        res = b - sl * a
        return sl, float(1.0 - (res @ res) / ((b - b.mean()) @ (b - b.mean()))), float(1.0 - (res @ res) / (b @ b))

    i12 = [odors.index(o) for o in H.MEGAMAT12]
    lab12 = H.abbrev_columns(lab16)[list(H.MEGAMAT12)].T.to_numpy(dtype=float)
    slope, r2, r2o = through_origin(M[i12], lab12)
    s17, _, r2o17 = through_origin(M, meg.T.to_numpy(dtype=float))
    out = dict(keys=kstr, labels=list(FLY_LABELS[:len(keys)]), odors=odors, odors12=list(H.MEGAMAT12), gloms=gloms,
               x=x * slope, slope=slope, slope_r2=r2, slope_r2_origin=r2o, slope17=s17, slope17_r2_origin=r2o17,
               fly_max=fly_max, avg_flymax=A, repeats=[int(r) for r in reps])
    if len(keys) != N_INPUT_FLIES:
        raise RuntimeError(f'olfdb dataset 1 has {len(keys)} flies, not {N_INPUT_FLIES}')
    if planted:
        return out
    _ORN.update(out)
    return _ORN


def reliabilities(F: Optional[dict] = None) -> np.ndarray:
    """E85's r_g of each glomerulus (one minus the variance of a trial mean over the variance across flies of the trial
    means, pooled over E85's 12 odors, clipped to [0, 1]; `E85.trial_stats`), times E85's common factor (part_a.json,
    near_check.c), clipped to [0, 1]."""
    from experiments import E85_orn_fly_population as E85m
    F = orn_flies() if F is None else F
    i12 = [F['odors'].index(o) for o in F['odors12']]
    r_g = E85m.trial_stats(F['x'][:, i12])['r_g']
    c = float(e85_part_a()['near_check']['c'])
    return np.clip(c * r_g, 0.0, 1.0)


def megamat_deviations(variant: str = DEVIATION, F: Optional[dict] = None) -> np.ndarray:
    """flies x megamat's 17 odors (`orn_flies()['odors']` order) x 36 glomeruli deviations: 'repro' r_g (`reliabilities`)
    times the fly's trial mean minus the mean over the flies that have the entry, 0 where the fly did not measure it;
    'zero' all 0. The sum over the flies of every entry is 0."""
    from experiments import E85_orn_fly_population as E85m
    F = orn_flies() if F is None else F
    if variant == 'zero':
        return np.zeros(F['x'].shape[:3])
    if variant != 'repro':
        raise ValueError(f'{variant=}')
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        st = E85m.trial_stats(F['x'])
    return reliabilities(F)[None, None, :] * st['D']


def deviations(panel: str, variant: str = DEVIATION, F: Optional[dict] = None) -> np.ndarray:
    """flies x odors x glomeruli deviations to add to `lab_input(panel)` (its odor and glomerulus order; 0 where a fly
    did not measure an entry, and on the APL16 odors of APL16_KEPT_AT_LAB); `variant` 'repro' (DEVIATION) or 'zero'."""
    F = orn_flies() if F is None else F
    dev = megamat_deviations(variant, F)
    lab = lab_input(panel)
    if [str(g) for g in lab.index] != F['gloms']:
        raise RuntimeError(f"{panel}: the lab input's glomeruli differ from the flies'")
    out = np.zeros((dev.shape[0], lab.shape[1], lab.shape[0]))
    for j, c in enumerate(lab.columns):
        o = H.odor_abbrev(c)
        if o in F['odors'] and not (panel == 'apl16' and o in APL16_KEPT_AT_LAB):
            if panel == 'apl16' and o not in F['odors12']:
                continue
            out[:, j, :] = dev[:, F['odors'].index(o), :]
    return out


def input_fly(fly: str, panel: str = 'megamat', variant: str = DEVIATION, F: Optional[dict] = None) -> pd.DataFrame:
    """Input fly `fly`'s model input on `panel` ('megamat' or 'apl16'): `lab_input(panel)` (same index, columns and
    dtypes) plus the fly's reproducible deviation where measured, the lab's value elsewhere; on APL16 paa, euc, -aPine
    and B-myr at the lab's input. `fly` is a label (f0 ... f8) or an E134 key. Pass it to E134's `build_family(...,
    orn=...)` and E77's `hold_rates` in place of the panel's input."""
    f = FLY_LABELS.index(fly_label(fly)) if F is None else F['labels'].index(fly)
    lab = lab_input(panel)
    dev = deviations(panel, variant, F)[f]
    out = lab.astype(float).copy()
    out.iloc[:, :] = lab.to_numpy(dtype=float) + dev.T
    return out


def input_flies_inputs(panel: str = 'megamat', variant: str = DEVIATION) -> Dict[str, pd.DataFrame]:
    """{label: `input_fly(label, panel, variant)`} for the nine input flies."""
    return {f: input_fly(f, panel, variant) for f in FLY_LABELS}


def zero_deviation_input(panel: str = 'megamat') -> pd.DataFrame:
    """`input_fly`'s construction with every deviation 0: equal to `lab_input(panel)` (checked), for the runner's
    identity (an input fly with no deviation reproduces E134's held megamat run with 0 differing counts)."""
    out = input_fly('f0', panel, 'zero')
    lab = lab_input(panel)
    if not (out.to_numpy(dtype=float) == lab.to_numpy(dtype=float)).all() or not out.index.equals(lab.index) \
            or not out.columns.equals(lab.columns):
        raise RuntimeError(f'{panel}: the zero-deviation input differs from the lab input')
    return out


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


def load_point(tag: str, step: Optional[int] = None):
    """A saved read (`kcq.stage_readings.SavedRun`; locally, else the archive); None if absent."""
    from kcq import stage_readings as R
    setting, panel, fly, wiring = parse_tag(tag)
    step = step_of(setting, panel) if step is None else step
    try:
        return R.load_step3_run(tag, dirs=run_dirs(step, setting, wiring, fly))
    except FileNotFoundError:
        return None


# ---------------------------------------------------------------------------
# reading model flies (no simulation): the engine shared by step 0 and the score

class Side:
    """Model flies read at one (stage, panel, condition) against one reference: per part, a model flies x entries
    array of each model fly's parts (`P.reference_views`), with the model flies' tags, wirings and input flies.

    Attributes: key, ref (hal's `FlyReference`), tags, wiring, fly (arrays of str), arrays (part -> n x entries)."""

    def __init__(self, key, ref, tags, wiring, fly, arrays):
        self.key, self.ref = tuple(key), ref
        self.tags, self.wiring, self.fly = list(tags), np.asarray(wiring, dtype=object), np.asarray(fly, dtype=object)
        self.arrays = {k: np.asarray(v, dtype=float) for k, v in arrays.items()}

    def __len__(self):
        return len(self.tags)

    def subset(self, mask) -> 'Side':
        m = np.asarray(mask, dtype=bool)
        return Side(self.key, self.ref, [t for t, k in zip(self.tags, m) if k], self.wiring[m], self.fly[m],
                    {p: A[m] for p, A in self.arrays.items()})


def _nanmean0(A: np.ndarray) -> np.ndarray:
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(A, axis=0)


def side_scores(side: Side, means=None) -> pd.DataFrame:
    """Each model fly's `Score` columns against the side's reference (`means`: per part, a mean per row in place of
    the flies', e.g. `others_means`), index the tags."""
    from kcq import population_fit as P
    m = {p: side.ref.mean[p].to_numpy(dtype=float) for p in side.arrays} if means is None else means
    return P.score_rows(side.arrays, m, {p: side.ref.one_fly[p] for p in side.arrays}, index=side.tags)


def population_score(side: Side) -> pd.Series:
    """The population's error: the model flies' mean views (NaN-aware, entry by entry) scored against the reference."""
    from kcq import population_fit as P
    A = {p: _nanmean0(X)[None, :] for p, X in side.arrays.items()}
    return P.score_arrays(A, side.ref, index=['population']).iloc[0]


def others_means(side: Side) -> Dict[str, np.ndarray]:
    """Per part, for each model fly the mean views of the model flies with another wiring and another input fly."""
    out = {}
    other = (side.wiring[:, None] != side.wiring[None, :]) & (side.fly[:, None] != side.fly[None, :])
    for p, A in side.arrays.items():
        fin = np.isfinite(A)
        S = np.where(fin, A, 0.0)
        num = other.astype(float) @ S
        den = other.astype(float) @ fin.astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[p] = np.where(den > 0, num / den, np.nan)
    return out


def against_others(side: Side) -> pd.DataFrame:
    """Each model fly against the others (another wiring and another input fly), in the stage's one-fly units."""
    return side_scores(side, others_means(side))


def stage_column(stage: str, views: Optional[str] = MEASURE) -> str:
    """The error column a stage is read on: 'view_strength' under the measure, every view's error ('error', the PN
    dendrites 'error_with_matrix', E136's `error_column` for views None) for the guard."""
    if views == MEASURE:
        return 'view_strength'
    return 'error_with_matrix' if stage == 'pn_dendrite' else 'error'


def _truthy(x) -> bool:
    if isinstance(x, str):
        return x.strip().lower() in ('true', '1', 'yes')
    return bool(x) if x is not None and not (isinstance(x, float) and np.isnan(x)) else False


def is_held(row: Mapping) -> bool:
    """Whether a read held: both response counts within HOLD_TOL of their targets and neither factor at a bound (from
    the row's n_on, n_off, target_on, target_off, at_bound_thr, at_bound_apl)."""
    try:
        ok = all(abs(float(row[f'n_{c}']) - float(row[f'target_{c}'])) <= HOLD_TOL * float(row[f'target_{c}'])
                 for c in ('on', 'off'))
    except (KeyError, TypeError, ValueError):
        return False
    return bool(ok and not _truthy(row.get('at_bound_thr')) and not _truthy(row.get('at_bound_apl')))


def terms(E: pd.DataFrame) -> dict:
    """The terms of a wirings x input flies table of a model fly's quantity (complete tables only): 'input_fly' the
    SD over input flies of the mean over wirings, 'wiring' the SD over wirings of the mean over input flies (both ddof
    1), 'remainder' the two-way residual's root mean square times sqrt(n_w n_f / ((n_w - 1)(n_f - 1))), 'ratio' the
    input-fly term over the wiring term, and the table's size."""
    X = E.to_numpy(dtype=float)
    nw, nf = X.shape
    if nw < 2 or nf < 2 or not np.isfinite(X).all():
        return dict(input_fly=np.nan, wiring=np.nan, remainder=np.nan, ratio=np.nan, n_wirings=nw, n_flies=nf,
                    complete=False)
    rm, cm, g = X.mean(axis=1, keepdims=True), X.mean(axis=0, keepdims=True), X.mean()
    res = X - rm - cm + g
    fly_t, wir_t = float(np.std(cm, ddof=1)), float(np.std(rm, ddof=1))
    rem = float(np.sqrt((res ** 2).sum() / ((nw - 1) * (nf - 1))))
    return dict(input_fly=fly_t, wiring=wir_t, remainder=rem, ratio=fly_t / wir_t if wir_t > 0 else np.inf,
                n_wirings=nw, n_flies=nf, complete=True)


def term_class(ratio: float) -> str:
    """The input-fly term's class from its ratio to the wiring term."""
    if not np.isfinite(ratio) and ratio != np.inf:
        return 'no class'
    return 'small' if ratio < TERM_SMALL else 'a part' if ratio < TERM_PART else "as large as the wiring's"


def flies_class(v: float) -> str:
    """The model flies' class from a model fly against the others."""
    if not np.isfinite(v):
        return 'no class'
    return 'vary too little' if v < FLIES_LOW else 'vary too much' if v > FLIES_HIGH else 'as the flies'


def stand_in_class(pop: float, lab: float, spread: float) -> str:
    """'a fair stand-in' if |pop - lab| <= spread, else 'moves it nearer' / 'moves it further'."""
    if not (np.isfinite(pop) and np.isfinite(lab) and np.isfinite(spread)):
        return 'no class'
    if abs(pop - lab) <= spread:
        return 'a fair stand-in'
    return 'moves it nearer' if pop < lab else 'moves it further'


def table_of(values: pd.Series, side: Side) -> pd.DataFrame:
    """A per-model-fly Series (index the side's tags) as a wirings x input flies table (WIRINGS, FLY_LABELS order,
    those present)."""
    df = pd.DataFrame(dict(v=values.reindex(side.tags).to_numpy(dtype=float), w=side.wiring, f=side.fly))
    t = df.pivot_table(index='w', columns='f', values='v', aggfunc='first')
    return t.reindex(index=[w for w in WIRINGS if w in t.index], columns=[f for f in FLY_LABELS if f in t.columns])


def somata_run(tag: str, wiring: str, pair: Mapping[str, tuple], source: str, panel: str = 'apl16'):
    """A somata-only `SavedRun` from {'on', 'off'}: (responses, spike counts) (E85's and E83's saved layout)."""
    from kcq import stage_readings as R
    return R.SavedRun(source=source, tag=tag, model='C2', wiring=wiring, input_panel=panel,
                      conditions={c: dict(responses=pair[c][0], spike_counts=pair[c][1]) for c in ('on', 'off')})


_READ_CACHE: Dict[str, object] = {}


def _fly_cache(half: Optional[str], stages: Sequence[str]):
    """The flies a model fly is read against: all the stage's flies (a FlyCache, E134's) for None; E136's half
    ('held' or 'fit') otherwise (E136's `Flies`, the half's curve and noise, all flies' units)."""
    from kcq import stage_readings as R
    key = f'{half}:{",".join(stages)}'
    if key not in _READ_CACHE:
        if half is None:
            _READ_CACHE[key] = R.FlyCache()
        else:
            from experiments import E136_stage_population_fit as E136m
            _READ_CACHE[key] = E136m.Flies(half, list(stages), views=None)
    return _READ_CACHE[key]


def read_sides(runs: Mapping[str, Tuple[str, str, object]], stages: Sequence[str], panel: str,
               half: Optional[str] = None, conditions: Sequence[str] = ('on', 'off')) -> Dict[tuple, Side]:
    """The model flies `runs` ({tag: (wiring, input fly, SavedRun)}) read at `stages` on `panel` (E134's default
    readings, seed index the wiring's, `_fly_cache(half)`), one `Side` per (stage, panel, condition) with flies: the
    reference E134's `fly_reference_for` on all the stage's flies (every part), or for a half E136's `Flies.reference`
    (every part, the half's mean in all flies' units)."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as E136m
    from kcq import population_fit as P
    from kcq import stage_readings as R
    lock = R.require_lock()
    fc = _fly_cache(half, stages)
    cache = fc if half is None else fc.cache
    pops: Dict[tuple, list] = {}
    for tag, (w, f, run) in runs.items():
        rd = R.read_run(run, stages=list(stages), lock=lock, flies=cache, seed_index=E136m.seed_index(w))
        for k, pop in rd.populations.items():
            if k[1] == panel and k[2] in conditions:
                pops.setdefault(k, []).append((tag, w, f, pop))
        del rd
    out = {}
    for k, items in pops.items():
        odors = list(items[0][3].values.columns)
        if half is None:
            fl = cache.get(k)
            if fl is None:
                continue
            _, ref, _ = B.fly_reference_for(fl, odors)
        else:
            if k[2] != 'on':
                continue
            ref = fc.reference(k[0], odors)
        arrays = P.stack_views([P.reference_views(it[3], ref) for it in items], ref)
        out[k] = Side(k, ref, [it[0] for it in items], [it[1] for it in items], [it[2] for it in items], arrays)
    return out


def stage_error_frame(sides: Mapping[str, Side], views: Optional[str] = MEASURE, kind: str = 'each') -> pd.Series:
    """The megamat stage error over the sides' stages (root mean square of the stages' errors): per model fly
    ('each'), the population's ('population', one value), or each model fly against the others ('others'); stage
    errors on `stage_column`."""
    cols = []
    for st, sd in sides.items():
        c = stage_column(st, views)
        if kind == 'each':
            cols.append(side_scores(sd)[c].rename(st))
        elif kind == 'others':
            cols.append(against_others(sd)[c].rename(st))
        else:
            cols.append(pd.Series({'population': float(population_score(sd)[c])}, name=st))
    df = pd.concat(cols, axis=1)
    return np.sqrt((df ** 2).mean(axis=1, skipna=False))


# ---------------------------------------------------------------------------
# the classes of E136, E139 and E141 recomputed with other spreads (a description)

def class_tables() -> Dict[str, dict]:
    """{experiment: classes.json} of CLASS_TABLES."""
    return {e: json.loads(H.results_dir(e).joinpath(*path).read_text()) for e, path in CLASS_TABLES}


def recompare(scored: Mapping, spread: Optional[float]) -> dict:
    """E136's `compare` of a scored comparison (its per-wiring changes kept) with another spread: the standard error
    the larger of half the spread and the SD of the changes (ddof 1) over sqrt(n); `spread` None keeps the scored one."""
    out = dict(scored)
    if not scored or not scored.get('complete') or spread is None:
        return out
    ch = np.asarray(scored['changes'], dtype=float)
    sd = float(np.std(ch, ddof=1)) if len(ch) > 1 else np.nan
    se = float(max(spread / 2.0 if np.isfinite(spread) else 0.0, sd / np.sqrt(len(ch)) if np.isfinite(sd) else 0.0))
    out.update(spread=float(spread), se=se, improves=bool(float(scored['change']) < -IMPROVE_SE * se),
               worse_than_spread=bool(np.isfinite(spread) and float(scored['change']) > spread))
    return out


def respread(rec: Mapping, spreads: Mapping[str, float], guard_spreads: Optional[Mapping[str, float]] = None) -> dict:
    """A class record with the spreads of `spreads` (keys as E136's classes.json 'spreads': 'stage', 'downstream',
    each stage, each 'input/panel/condition') in its comparisons (megamat, stages, panels, panels beside, and E141's
    against_a1), and `guard_spreads` in its guard (every view), each recomputed by `recompare`; others as scored."""
    import copy
    r = copy.deepcopy(dict(rec))
    gs = guard_spreads or {}
    if 'megamat' in r:
        r['megamat'] = recompare(r['megamat'], spreads.get('stage'))
    for s_, c in list((r.get('stages') or {}).items()):
        r['stages'][s_] = recompare(c, spreads.get(s_))
    for grp in ('panels', 'panels_beside'):
        for k, c in list((r.get(grp) or {}).items()):
            r[grp][k] = recompare(c, spreads.get(k))
    for k, c in list((r.get('guard') or {}).items()):
        r['guard'][k] = recompare(c, gs.get(k))
    vs = r.get('against_a1')
    if vs:
        if vs.get('megamat'):
            vs['megamat'] = recompare(vs['megamat'], spreads.get('stage'))
        for k, c in list((vs.get('panels') or {}).items()):
            vs['panels'][k] = recompare(c, spreads.get(k))
    return r


def _changes(c: Optional[Mapping]) -> np.ndarray:
    """A comparison's per-wiring changes as one reference (1 x wirings); NaN where it is not complete (E143's `_cmp`
    then reads it as missing from its record, whatever the values)."""
    if c and c.get('complete') and 'changes' in c:
        return np.asarray(c['changes'], dtype=float)[None, :]
    return np.full((1, len(T_WIRINGS)), np.nan)


def reclass_record(rec: Mapping, table: str) -> dict:
    """One class record read by E143's ladder (`reclass`, and for E141 `against_a1` and `e141_entry`) on its own
    comparisons' changes, spreads and standard errors: {'cls', and E141's 'stands_in' / 'does_more_than_a1' where the
    record has them}."""
    from experiments import E143_megamat_fly_sampling as E143
    name, fam = 'x', 'family'

    def pair(c):
        a = _changes(c)
        return {name: a, fam: np.zeros_like(a)}

    mega = {'stage': pair(rec.get('megamat'))}
    mega.update({s_: pair(c) for s_, c in (rec.get('stages') or {}).items()})
    guard = {k: pair(c) for k, c in (rec.get('guard') or {}).items() if k in LOSS_STAGES}
    panels = {k: pair(c) for grp in ('panels', 'panels_beside') for k, c in (rec.get(grp) or {}).items()}
    gpanels = {k: pair(c) for k, c in (rec.get('guard') or {}).items() if k not in LOSS_STAGES}
    rv = E143.Resampled(mega, guard, panels, gpanels)
    cls, _ = E143.reclass(rec, name, rv, table)
    out = dict(cls=str(cls[0]))
    if table == 'E141_kckc_rise_forms' and 'against_a1' in rec:
        vs = rec.get('against_a1') or {}

        def pair_a1(c):
            a = _changes(c)
            return {name: a, 'E136:A1': np.zeros_like(a)}

        rv2 = E143.Resampled({'stage': pair_a1(vs.get('megamat'))}, {},
                             {k: pair_a1(c) for k, c in (vs.get('panels') or {}).items()}, {})
        vsc = E143.against_a1(rec, name, rv2)
        extra = E143.e141_entry(rec, cls, vsc)
        out.update({k: bool(v[0]) for k, v in extra.items()})
    return out


def classes_with_spreads(spreads: Mapping[str, float], guard_spreads: Optional[Mapping[str, float]] = None,
                         tables: Optional[Mapping[str, dict]] = None) -> pd.DataFrame:
    """Every class of E136's, E139's and E141's tables recomputed with `spreads` (and `guard_spreads`), as a
    description: per record its scored class, the class with the spreads as scored (which must be the scored one) and
    the class with the new spreads, and E141's added classes likewise."""
    tables = class_tables() if tables is None else tables
    rows = []
    for exp, tab in tables.items():
        for arm, rec in tab['classes'].items():
            same = reclass_record(rec, exp)
            new = reclass_record(respread(rec, spreads, guard_spreads), exp)
            row = dict(table=exp, arm=arm, cls_scored=rec.get('cls'), cls_as_scored=same['cls'], cls_new=new['cls'],
                       changed=new['cls'] != rec.get('cls'))
            for k in ('stands_in', 'does_more_than_a1'):
                if k in same:
                    row.update({f'{k}_scored': rec.get(k), f'{k}_as_scored': same[k], f'{k}_new': new[k]})
                    row['changed'] = row['changed'] or (new[k] != bool(rec.get(k)))
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the parts

def lock_files() -> Dict[str, Path]:
    """The files the lock records by sha256 (name -> path)."""
    r = H.shared_dir() / 'results'
    files = {
        'e134_lock': H.HERE / 'experiments' / 'E134_lock.json',
        'olfdb_dataset_1': olfdb_copy(),
        'e85_part_a': r / E85 / 'part_a.json',
        'e85_runs': r / E85 / 'runs.csv',
        'e85_flies': r / E85 / 'flies.csv',
        'e134_writeup_stage_error': r / E134 / 'writeup_stage_error.csv',
        'e134_score_classes': r / E134 / 'score_classes.json',
        'e136_classes': r / E136 / E136_MEASURE_DIR / 'classes.json',
        'e136_fit_a1': r / E136 / E136_MEASURE_DIR / 'fit_A1.json',
        'e136_score_panels': r / E136 / E136_MEASURE_DIR / 'score_panels.csv',
        'e136_score_megamat': r / E136 / E136_MEASURE_DIR / 'score_megamat.csv',
        'e139_classes': r / 'E139_normalization_ladder' / 'classes.json',
        'e141_classes': r / 'E141_kckc_rise_forms' / 'classes.json',
        'e140_input_aa': r / 'E140_aa_claw_soma_step' / 'input_aa.csv',
        'e140_provenance_input': r / 'E140_aa_claw_soma_step' / 'provenance_input.json',
        'e143_classes': r / 'E143_megamat_fly_sampling' / 'classes.json',
    }
    for f in FLY_LABELS:
        files[f'e85_input_repro_{f}'] = r / E85 / 'input_flies' / f'repro_{f}.csv'
    return files


def lock_contents() -> dict:
    """The lock's contents (`part_lock`)."""
    import socket
    from datetime import datetime
    F = orn_flies()
    e134 = json.loads((H.HERE / 'experiments' / 'E134_lock.json').read_text())
    orn = e134['populations']['orn__megamat__wild_type']
    if list(orn['flies']) != F['keys']:
        raise RuntimeError(f"E134's ORN flies {orn['flies']} are not the copy's {F['keys']}")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        m = np.nanmean(F['x'], axis=3)
    given = {F['labels'][i]: [o for j, o in enumerate(F['odors']) if np.isfinite(m[i, j]).any()]
             for i in range(len(F['labels']))}
    a1 = json.loads((H.results_dir(E136) / E136_MEASURE_DIR / 'fit_A1.json').read_text())
    lab16 = lab_input('apl16')
    return dict(
        experiment=EXPERIMENT, step='lock', script=SCRIPT, runner=RUNNER, host=socket.gethostname().split('.')[0],
        time=datetime.now().strftime('%Y-%m-%d %H:%M'), commits=H.code_commits(),
        files={k: file_record(v) for k, v in lock_files().items()},
        input_flies=dict(
            source="olfdb dataset 1 (diagnostics, megamat; 9 flies), E134's copy (data/stage_populations/sources/olfdb/"
                   "ds1_long.csv); E134's lock's orn / megamat / wild_type flies, sorted, labelled f0 ... f8 (E85's order)",
            flies=[dict(label=l, key=k, odors_given=given[l]) for l, k in zip(F['labels'], F['keys'])],
            deviation=DEVIATION, slope=F['slope'], slope_rule="E85's: least squares through the origin of the lab's "
            "APL16 input on the flies' mean over the 12 shared odors", slope17_beside=F['slope17'],
            common_factor=float(e85_part_a()['near_check']['c']), avg_flymax=F['avg_flymax'],
            r_g=dict(zip(F['gloms'], [float(x) for x in reliabilities(F)]))),
        odors=dict(megamat=F['odors'], apl16=[H.odor_abbrev(c) for c in lab16.columns], shared=F['odors12'],
                   apl16_kept_at_lab=list(APL16_KEPT_AT_LAB)),
        glomeruli=F['gloms'],
        wirings=dict(S=list(S_WIRINGS), T=list(T_WIRINGS), all=list(WIRINGS),
                     seed_index="E136's seed_index: the wiring's position in WIRINGS (pooling and noise), the same for "
                                "every input fly on that wiring", claws="'first', pinned (the family's fits)"),
        settings=dict(family="the unrefitted family at each wiring's held megamat point (E136's FAMILY, E134's held "
                             "targets on megamat, the wiring's own fit's counts on APL16)",
                      A1="E136's A1 at its classed point (measure_strength/fit_A1.json classed_z)",
                      A1_classed_z=a1.get('classed_z') or a1.get('best_z'),
                      free="'familyfree' / 'A1free': each model fly at its wiring's held parameters, nothing re-held"),
        points=dict(tag='<setting>_<panel>_<inputfly>_<wiring>',
                    step1=[point_tag(p.setting, p.panel, p.fly, p.wiring) for p in planned_points(1)],
                    step2=[point_tag(p.setting, p.panel, p.fly, p.wiring) for p in planned_points(2)],
                    order='T wirings first (step 1: T then S; step 2 on T only), fly by fly within a wiring'),
        reading=dict(readings="kcq.stage_readings DEFAULT_READINGS (E134's lock)", measure=MEASURE,
                     stages=list(LOSS_STAGES), flies='all the stage\'s flies (E134); E136\'s held-out half beside '
                     '(E136 Flies("held"), HALF_SEED 0)', half_seed=0),
        classes=dict(
            fair_stand_in="the population's error within the spread as scored of the error of the wirings' mean view at "
                          "the lab's input; otherwise 'the population reading moves it', nearer or further",
            input_fly_term=dict(small_below=TERM_SMALL, part_below=TERM_PART, rule='input-fly term / wiring term'),
            model_flies=dict(too_little_below=FLIES_LOW, too_much_above=FLIES_HIGH),
            a1_fall=dict(se_multiple=IMPROVE_SE, control_fall=CONTROL_FALL),
            a1_apl16='no class: the share of input flies on which A1 is worse than the spread as scored',
            unheld_max_share=UNHELD_MAX_SHARE, hold_tol=HOLD_TOL),
        revisions=[dict(time='2026-11-27 09:45', what="'unheld' read with E136's hold tolerance (HOLD_TOL, 0.5 % of each "
                        "target) in place of exact counts", why="step 0's first run (09:39) found that E77's hold_rates, "
                        "which E85's runs and E136's reads use, misses a target by one or two counts on 12 of E85's 72 "
                        "'repro' runs (at most 0.1 % of ~1900 APL-on and ~4100 APL-off responses), as a bisection on "
                        "discrete counts can; the entry's 'unheld' means a hold that failed, not a count step. No step 1 or "
                        "2 read existed. The first lock is commit 2a58895")],
        spreads_as_scored=json.loads((H.results_dir(E136) / E136_MEASURE_DIR / 'classes.json').read_text())['spreads'],
        e134_identity=E134_IDENTITY, identity_tol=IDENTITY_TOL, e85_groups=E85_GROUPS,
        choices=list(CHOICES),
        holds=['no simulation in the lock, step 0 or the score; no olfsysm change', 'nothing searched, fitted or weighted',
               'validation2 not read; E145\'s cells not read', "Pratyush's hemibrain only",
               "Remy's processed data without the bouton fly 2022-10-21 fly02 and her older megamat flies (E134's lock)"],
    )


def part_lock() -> dict:
    """The lock: every file, fly set, odor set, seed, class threshold, the T-first order and the choices; written to
    results/.../lock.json and LOCK_COPY (commit the copy before any other part)."""
    lock = lock_contents()
    txt = json.dumps(lock, indent=1, default=float) + '\n'
    LOCK_COPY.write_text(txt)
    lock_path().write_text(txt)
    print(f'wrote {LOCK_COPY} and {lock_path()}; commit {LOCK_COPY.relative_to(H.HERE)} before any other part',
          flush=True)
    return lock


def _host() -> str:
    import socket
    return socket.gethostname().split('.')[0]


def _now() -> str:
    from datetime import datetime
    return datetime.now().strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    tmp = path.with_name(f'.{path.name}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')
    tmp.replace(path)


def _json_default(x):
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    return str(x)


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    path = results_dir() / name
    tmp = path.with_name(f'.{path.name}.tmp')
    df.to_csv(tmp, index=False)
    tmp.replace(path)
    return path


def provenance(part: str, lock: Mapping, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time, the lock's sha256 and `extra`."""
    _write_json(results_dir() / f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
                     lock_sha256=lock.get('sha256'), **extra))


# ---------------------------------------------------------------------------
# step 0

def e134_held_runs(wirings: Sequence[str] = WIRINGS) -> Dict[str, Tuple[str, str, object]]:
    """E134's held megamat runs (the family at the lab's input, step 3's megamat_held_<wiring>) as model flies of the
    input fly 'lab'."""
    from experiments import E134_population_baseline as B
    out = {}
    for w in wirings:
        r = B.load_run(B.tag_of(B.Point('megamat', 'held', w)))
        if r is None:
            raise RuntimeError(f"E134's held megamat run on {w} is missing")
        out[f'lab_megamat_{w}'] = (w, 'lab', r)
    return out


def lab_apl16_runs(wirings: Sequence[str] = WIRINGS) -> Dict[str, Tuple[str, str, object]]:
    """The family at the lab's APL16 input (each wiring's own fit's saved outputs, E83's `saved_fit`) as somata runs."""
    from experiments import E83_cross_wiring_fit as E83
    return {f'lab_apl16_{w}': (w, 'lab', somata_run(f'lab_apl16_{w}', w, E83.saved_fit(w)['saved'], 'E83:saved_fit'))
            for w in wirings}


def identity_e134(lock: Mapping) -> dict:
    """E134's 0.918 and 0.894: the eight held megamat runs read by this module (all flies, strength by odor) give the
    mean over the wirings of each wiring's megamat stage error and the error of their mean views of E134's
    writeup_stage_error.csv (megamat_held), each stage's per-wiring and common values beside; and E136's held-out half
    (score_megamat.csv, the family's 'held' rows) read the same way."""
    runs = e134_held_runs()
    sides = {k[0]: v for k, v in read_sides(runs, LOSS_STAGES, 'megamat', conditions=('on',)).items()}
    each, pop = stage_error_frame(sides), float(stage_error_frame(sides, kind='population').iloc[0])
    w = pd.read_csv(resolve(lock['files']['e134_writeup_stage_error']['path']))
    w = w[w.group == 'megamat_held'].set_index('wiring')
    tags = [f'lab_megamat_{x}' for x in WIRINGS]
    diffs = {}
    for st, sd in sides.items():
        sc = side_scores(sd)['view_strength'].reindex(tags).to_numpy()
        diffs[st] = dict(each=float(np.abs(sc - w.loc[list(WIRINGS), f'{st}_strength'].to_numpy()).max()),
                         common=float(abs(float(population_score(sd)['view_strength']) - w.loc['common', f'{st}_strength'])))
    held = {k[0]: v for k, v in read_sides(runs, LOSS_STAGES, 'megamat', half='held', conditions=('on',)).items()}
    he, hp = stage_error_frame(held), float(stage_error_frame(held, kind='population').iloc[0])
    sm = pd.read_csv(resolve(lock['files']['e136_score_megamat']['path']))
    sm = sm[(sm.arm == 'family') & (sm.half == 'held')].set_index('wiring')
    hd = float(np.abs(he.reindex(tags).to_numpy() - sm.loc[list(WIRINGS), 'stage'].to_numpy()).max())
    hcd = float(abs(hp - sm.loc['common', 'stage']))
    res = dict(model=float(each.mean()), common=pop, expected=lock['e134_identity'],
               model_diff=float(abs(each.mean() - lock['e134_identity']['model'])),
               common_diff=float(abs(pop - lock['e134_identity']['common'])), stages=diffs,
               held_half=dict(model=float(he.mean()), common=hp, max_wiring_diff_from_e136=hd, common_diff_from_e136=hcd))
    tol = float(lock['identity_tol'])
    res['ok'] = bool(res['model_diff'] <= tol and res['common_diff'] <= tol
                     and all(v['each'] <= tol and v['common'] <= tol for v in diffs.values()) and hd <= tol and hcd <= tol)
    res['_sides'] = sides
    res['_held'] = held
    return res


def identity_e85(lock: Mapping) -> dict:
    """E85's part A on the 12 shared odors from E134's copy: the slope, every glomerulus's r_g (fly_spread.csv), the
    nine 'repro' APL16 inputs (input_flies/repro_<fly>.csv) and the mean over flies being the lab's input; beside, E134's
    copy against E85's own cache of olfdb dataset 1 where this host has it."""
    F = orn_flies()
    pa = json.loads(resolve(lock['files']['e85_part_a']['path']).read_text())
    sp = pd.read_csv(H.results_dir(E85) / 'fly_spread.csv')
    rg85 = sp.groupby('glomerulus').r_g.first().reindex(F['gloms']).to_numpy(dtype=float)
    from experiments import E85_orn_fly_population as E85m
    i12 = [F['odors'].index(o) for o in F['odors12']]
    rg = E85m.trial_stats(F['x'][:, i12])['r_g']
    inputs = {}
    for f in FLY_LABELS:
        mine = H.abbrev_columns(input_fly(f, 'apl16'))
        e85 = pd.read_csv(resolve(lock['files'][f'e85_input_repro_{f}']['path']), index_col=0)
        inputs[f] = float(np.abs(mine.loc[[str(g) for g in e85.index], list(e85.columns)].to_numpy() - e85.to_numpy()).max())
    mean_dev = {pn: float(np.abs(np.mean([input_fly(f, pn).to_numpy(dtype=float) for f in FLY_LABELS], axis=0)
                                 - lab_input(pn).to_numpy(dtype=float)).max()) for pn in PANELS}
    cache = H.HERE / 'outputs' / E85 / 'olfdb_dataset1.parquet'
    same_copy = None
    if cache.exists():
        a = pd.read_csv(olfdb_copy()).sort_values(['date', 'fly_num', 'roi', 'panel', 'odor1', 'repeat']).reset_index(drop=True)
        b = pd.read_parquet(cache)[list(a.columns)].sort_values(['date', 'fly_num', 'roi', 'panel', 'odor1', 'repeat']
                                                                ).reset_index(drop=True)
        b['date'] = b['date'].astype(str)
        same_copy = bool(len(a) == len(b) and (a.drop(columns='value').astype(str).to_numpy()
                                                == b.drop(columns='value').astype(str).to_numpy()).all()
                         and np.allclose(a.value.to_numpy(dtype=float), b.value.to_numpy(dtype=float), rtol=0, atol=1e-12))
    res = dict(slope=F['slope'], slope_e85=float(pa['scaling']['slope']),
               slope_rel_diff=float(abs(F['slope'] / float(pa['scaling']['slope']) - 1.0)),
               r_g_max_diff=float(np.nanmax(np.abs(rg - rg85))), inputs_max_abs_diff=inputs,
               mean_over_flies_minus_lab=mean_dev, e134_copy_equals_e85_cache=same_copy,
               common_factor=float(pa['near_check']['c']))
    res['ok'] = bool(res['slope_rel_diff'] <= 1e-12 and res['r_g_max_diff'] <= 1e-12
                     and max(inputs.values()) <= 1e-9 and max(mean_dev.values()) <= 1e-9 and same_copy is not False)
    return res


def identity_e140(lock: Mapping) -> dict:
    """E140's aa vectors: each fly's aa trial means on the 36 glomeruli by this module's map, at E140's slope (the
    17-odor regression, `orn_flies()['slope17']`, which must be E140's), equal E140's input_aa.csv orn_<fly> columns;
    the flies given aa are E140's."""
    F = orn_flies()
    prov = json.loads(resolve(lock['files']['e140_provenance_input']['path']).read_text())['orn_map']
    tab = pd.read_csv(resolve(lock['files']['e140_input_aa']['path'])).set_index('glomerulus')
    tab.index = tab.index.astype(str)
    j = F['odors'].index('aa')
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        m = np.nanmean(F['x'][:, j], axis=2) / F['slope'] * F['slope17']      # flies x glomeruli
    given = [F['keys'][i] for i in range(len(F['keys'])) if np.isfinite(m[i]).any()]
    diffs = {}
    for i, k in enumerate(F['keys']):
        col = f'orn_{k}'
        if col not in tab.columns:
            continue
        e = tab[col].reindex(F['gloms']).to_numpy(dtype=float)
        ok = np.isfinite(e) | np.isfinite(m[i])
        diffs[k] = float(np.nanmax(np.abs(np.where(ok, m[i], 0) - np.where(ok, e, 0)) / np.maximum(1.0, np.abs(e))))
    res = dict(slope17=F['slope17'], slope_e140=float(prov['slope']),
               slope_rel_diff=float(abs(F['slope17'] / float(prov['slope']) - 1.0)), flies_given_aa=given,
               flies_given_aa_e140=list(prov['flies_given_aa']), max_rel_diff=diffs)
    res['ok'] = bool(res['slope_rel_diff'] <= 1e-12 and given == list(prov['flies_given_aa'])
                     and diffs and max(diffs.values()) <= 1e-9)
    return res


def identity_apl16(lock: Mapping) -> dict:
    """The APL16 reading of this module: the family at the lab's APL16 input (E83's saved fits) read at the somata
    gives E136's score_panels.csv family rows (each wiring's strength error and the common error, APL on and off)."""
    sides = read_sides(lab_apl16_runs(), ['kc_soma'], 'apl16')
    sp = pd.read_csv(resolve(lock['files']['e136_score_panels']['path']))
    sp = sp[(sp.arm == 'family') & (sp.input == 'apl16') & (sp.panel == 'apl16')]
    out = {}
    tags = [f'lab_apl16_{w}' for w in WIRINGS]
    for k, sd in sides.items():
        ref = sp[sp.condition == k[2]].set_index('wiring')
        sc = side_scores(sd)['view_strength'].reindex(tags).to_numpy()
        out[k[2]] = dict(each=float(np.abs(sc - ref.loc[list(WIRINGS), 'error'].to_numpy()).max()),
                         common=float(abs(float(population_score(sd)['view_strength']) - ref.loc['common', 'error'])))
    res = dict(diffs=out, ok=bool(len(out) == 2 and all(v['each'] <= 1e-9 and v['common'] <= 1e-9 for v in out.values())))
    res['_sides'] = sides
    return res


def input_flies_table() -> pd.DataFrame:
    """Megamat's input flies on its 17 odors: per fly and odor whether it was given, its measured glomeruli, the root
    mean square of its deviation and of the lab's input over the 36 glomeruli, their ratio, the Pearson r of its input
    with the lab's, its input's 1 - Pearson distance from the lab's, and the entries below the lab input's minimum; a
    row per fly over all odors ('all')."""
    F = orn_flies()
    lab = lab_input('megamat')
    labA = H.abbrev_columns(lab)
    dev = megamat_deviations()
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        m = np.nanmean(F['x'], axis=3)
    rows = []
    lmin = float(lab.to_numpy(dtype=float).min())
    for i, f in enumerate(F['labels']):
        inp = H.abbrev_columns(input_fly(f, 'megamat'))
        for j, o in enumerate(F['odors']):
            d = dev[i, j]
            L = labA[o].reindex(F['gloms']).to_numpy(dtype=float)
            x = inp[o].reindex(F['gloms']).to_numpy(dtype=float)
            r = float(np.corrcoef(x, L)[0, 1])
            rows.append(dict(fly=f, key=F['keys'][i], odor=o, given=bool(np.isfinite(m[i, j]).any()),
                             n_glomeruli=int(np.isfinite(m[i, j]).sum()), rms_deviation=float(np.sqrt(np.mean(d ** 2))),
                             rms_lab=float(np.sqrt(np.mean(L ** 2))),
                             deviation_share=float(np.sqrt(np.mean(d ** 2)) / np.sqrt(np.mean(L ** 2))),
                             pearson_with_lab=r, distance_from_lab=1.0 - r, n_below_lab_min=int((x < lmin).sum()),
                             shared_with_apl16=o in F['odors12']))
        t = pd.DataFrame([q for q in rows if q['fly'] == f])
        rows.append(dict(fly=f, key=F['keys'][i], odor='all', given=bool(t.given.all()), n_glomeruli=int(t.n_glomeruli.sum()),
                         rms_deviation=float(np.sqrt(np.mean(dev[i] ** 2))),
                         rms_lab=float(np.sqrt(np.mean(labA.to_numpy(dtype=float) ** 2))),
                         deviation_share=float(np.sqrt(np.mean(dev[i] ** 2)) / np.sqrt(np.mean(labA.to_numpy(dtype=float) ** 2))),
                         pearson_with_lab=float(t.pearson_with_lab.mean()), distance_from_lab=float(t.distance_from_lab.mean()),
                         n_below_lab_min=int(t.n_below_lab_min.sum()), shared_with_apl16=False))
    return pd.DataFrame(rows)


def orn_stage_runs(inputs: Mapping[str, pd.DataFrame]) -> Dict[str, Tuple[str, str, object]]:
    """Inputs as runs read at the ORN stage alone (the model's ORN stage is its input; the wiring 'real' gives the seed
    index 0, the ORN stage takes no seed)."""
    from kcq import stage_readings as R
    out = {}
    for name, orn in inputs.items():
        o = H.abbrev_columns(orn.astype(float))
        dummy = pd.DataFrame(np.zeros((1, o.shape[1])), columns=list(o.columns))
        out[name] = ('real', name, R.SavedRun(source='E147:input', tag=name, model='input', wiring='real',
                                              input_panel='megamat', conditions={'on': dict(spike_counts=dummy)}, orn=o))
    return out


def orn_stage_table() -> pd.DataFrame:
    """Each input fly's megamat input, and the lab's, read at the ORN stage (E134's ORN reading) against the nine ORN
    flies, in the ORN stage's one-fly units: every view's error and the strength view's (beside, classing nothing; it
    needs no run)."""
    inputs = {'lab': lab_input('megamat'), **{f: input_fly(f, 'megamat') for f in FLY_LABELS}}
    sides = read_sides(orn_stage_runs(inputs), ['orn'], 'megamat', conditions=('on',))
    sd = sides[('orn', 'megamat', 'on')]
    sc = side_scores(sd)
    sub = sd.subset([t != 'lab' for t in sd.tags])
    pop = population_score(sub)
    rows = [dict(unit=t, **{k: float(v) for k, v in sc.loc[t].items()}) for t in sd.tags]
    rows.append(dict(unit='population_of_nine', **{k: float(v) for k, v in pop.items()}))
    return pd.DataFrame(rows)


def e85_runs(groups: Sequence[str] = tuple(E85_GROUPS)) -> Tuple[pd.DataFrame, Dict[str, Tuple[str, str, object]]]:
    """E85's saved model flies of `groups` (runs.csv, stage B1S and B1T as 'B1') as somata runs, and their rows with
    'group' and 'held' (counts at both targets and no factor at a bound; True for the free reading B5)."""
    from experiments import E85_orn_fly_population as E85m
    rt = pd.read_csv(H.results_dir(E85) / 'runs.csv')
    rt['group'] = [E85m.group_of(x) for x in rt.stage]
    rt = rt[rt.group.isin(groups)].reset_index(drop=True)
    rt['held'] = [True if r.reading == 'free' else is_held(r._asdict()) for r in rt.itertuples()]
    runs = {}
    for r in rt.itertuples():
        o = E85m.load_run('B', r.tag)
        if o is None:
            raise RuntimeError(f'E85 {r.tag}: outputs missing')
        runs[r.tag] = (str(r.wiring), str(r.fly), somata_run(r.tag, str(r.wiring), o, 'E85:B'))
    return rt, runs


def apl16_group_rows(group: str, sides: Mapping[tuple, Side], lab_sides: Mapping[tuple, Side], held: pd.Series,
                     spreads: Mapping[str, float]) -> Tuple[List[dict], pd.DataFrame]:
    """One group of model flies at the APL16 somata, per condition: the population's error against the wirings' mean
    view at the lab's input on the same wirings ('a fair stand-in'), the terms, a model fly against the others, the
    spreads over model flies, each on strength by odor and on every view; and the per-model-fly rows."""
    rows, per = [], []
    for cond in ('on', 'off'):
        key = ('kc_soma', 'apl16', cond)
        sd = sides[key]
        keep = np.array([bool(held.get(t, True)) for t in sd.tags])
        n_unheld = int((~keep).sum())
        sd = sd.subset(keep)
        ws = [w for w in WIRINGS if w in set(sd.wiring)]
        lab = lab_sides[key].subset([w in ws for w in lab_sides[key].wiring])
        each, oth = side_scores(sd), against_others(sd)
        lab_each = side_scores(lab)
        pop, labpop = population_score(sd), population_score(lab)
        k = f'apl16/apl16/{cond}'
        rec = dict(group=group, condition=cond, n_model_flies=len(sd), n_unheld=n_unheld,
                   no_class=bool(n_unheld > UNHELD_MAX_SHARE * (len(sd) + n_unheld)), wirings=' '.join(ws))
        for vname, col in (('strength', 'view_strength'), ('all_views', 'error')):
            T = table_of(each[col], sd)
            tm = terms(T)
            o = oth[col]
            rec.update({
                f'{vname}_population': float(pop[col]), f'{vname}_lab_common': float(labpop[col]),
                f'{vname}_population_minus_lab': float(pop[col] - labpop[col]),
                f'{vname}_mean_model_fly': float(each[col].mean()), f'{vname}_lab_mean_wiring': float(lab_each[col].mean()),
                f'{vname}_spread_model_flies': float(each[col].std(ddof=1)),
                f'{vname}_spread_wirings_lab': float(lab_each[col].std(ddof=1)) if len(lab_each) > 1 else np.nan,
                f'{vname}_term_input_fly': tm['input_fly'], f'{vname}_term_wiring': tm['wiring'],
                f'{vname}_term_remainder': tm['remainder'], f'{vname}_term_ratio': tm['ratio'],
                f'{vname}_terms_complete': tm['complete'],
                f'{vname}_against_others_mean': float(o.mean()), f'{vname}_against_others_min': float(o.min()),
                f'{vname}_against_others_max': float(o.max())})
        rec.update(spread_as_scored=float(spreads.get(k, np.nan)),
                   stand_in=('no class' if rec['no_class'] else
                             stand_in_class(rec['strength_population'], rec['strength_lab_common'], spreads.get(k, np.nan))),
                   input_fly_term=('no class' if rec['no_class'] or not rec['strength_terms_complete']
                                   else term_class(rec['strength_term_ratio'])),
                   model_flies=('no class' if rec['no_class'] else flies_class(rec['strength_against_others_mean'])))
        rows.append(rec)
        per.append(pd.DataFrame(dict(group=group, condition=cond, tag=sd.tags, wiring=sd.wiring, fly=sd.fly,
                                     strength=each['view_strength'].to_numpy(), all_views=each['error'].to_numpy(),
                                     sharing=each['view_sharing'].to_numpy(), breadth=each['view_breadth'].to_numpy(),
                                     against_others_strength=oth['view_strength'].to_numpy(),
                                     against_others_all_views=oth['error'].to_numpy())))
    return rows, pd.concat(per, ignore_index=True)


def apl16_realism_rows(rt: pd.DataFrame, runs: Mapping[str, Tuple[str, str, object]]) -> pd.DataFrame:
    """The four realism bands per model fly (E136's `apl16_realism`) and per population (E85's `population_realism` on
    the mean over its model flies), per E85 group; the wirings at the lab's input beside."""
    from experiments import E136_stage_population_fit as E136m
    from experiments import E85_orn_fly_population as E85m
    from experiments import E83_cross_wiring_fit as E83
    rows, parts = [], {}
    for r in rt.itertuples():
        run = runs[r.tag][2]
        on, off = run.conditions['on']['responses'], run.conditions['off']['responses']
        rr = E136m.apl16_realism(on, off)
        rows.append(dict(group=r.group, unit=r.tag, kind='model_fly', wiring=r.wiring, fly=r.fly, n_model_flies=1,
                         max_count_on=float(run.conditions['on']['spike_counts'].to_numpy(dtype=float).max()),
                         max_count_off=float(run.conditions['off']['spike_counts'].to_numpy(dtype=float).max()),
                         thr_factor=float(r.thr_factor), apl_factor=float(r.apl_factor), **rr))
        parts[r.tag] = E85m.run_realism_parts({'on': (on, run.conditions['on']['spike_counts']),
                                               'off': (off, run.conditions['off']['spike_counts'])})
    for w in WIRINGS:
        o = E83.saved_fit(w)['saved']
        rows.append(dict(group='lab', unit=w, kind='wiring', wiring=w, fly='lab', n_model_flies=1,
                         max_count_on=float(o['on'][1].to_numpy(dtype=float).max()),
                         max_count_off=float(o['off'][1].to_numpy(dtype=float).max()), **E136m.apl16_realism(o['on'][0], o['off'][0])))
    for g, gt in rt.groupby('group', sort=False):
        pr = E85m.population_realism([parts[t] for t in gt.tag])
        rows.append(dict(group=g, unit='population', kind='population', n_model_flies=len(gt),
                         rate_on=pr['rate_on'], rate_off=pr['rate_off'], rate_ratio=pr['ratio'],
                         rate_rho_on=pr['rate_rho_on'], rate_rho_off=pr['rate_rho_off'], opk_on_tv=pr['opk_on_tv'],
                         opk_off_tv=pr['opk_off_tv'], realism_flag=pr['realism_flag']))
    return pd.DataFrame(rows)


def probe_numbers(lock: Mapping) -> pd.DataFrame:
    """The planner's probe recomputed as the entry describes it (E134's archive scoring, archive_scores_curve_noise.csv,
    E85's B1 tags, view_strength; the archive reads every variant at noise seed 0), beside this module's reading."""
    a = pd.read_csv(H.results_dir(E134) / 'archive_scores_curve_noise.csv', dtype={'tag': str})
    a = a[(a.folder == E85) & a.tag.str.startswith('B1')]
    rows = []
    for cond, g in a.groupby('condition'):
        parts = g.tag.str.extract(r'B1[ST]_repro_(f\d)_(\w+)$')
        t = pd.DataFrame(dict(v=g.view_strength.to_numpy(), f=parts[0].to_numpy(), w=parts[1].to_numpy()))
        T = t.pivot_table(index='w', columns='f', values='v')
        rows.append(dict(condition=cond, n=len(g), mean=float(g.view_strength.mean()),
                         sd_over_input_flies=float(T.mean(axis=0).std(ddof=1)),
                         sd_over_wirings=float(T.mean(axis=1).std(ddof=1))))
    return pd.DataFrame(rows)


def model_fly_spreads_apl16(groups: pd.DataFrame, group: str = 'B1') -> Tuple[Dict[str, float], Dict[str, float]]:
    """(strength spreads, every-view spreads) over the model flies of `group` at the APL16 somata, keyed as E136's
    spreads ('apl16/apl16/on', '/off')."""
    g = groups[groups.group == group].set_index('condition')
    st = {f'apl16/apl16/{c}': float(g.loc[c, 'strength_spread_model_flies']) for c in ('on', 'off') if c in g.index}
    al = {f'apl16/apl16/{c}': float(g.loc[c, 'all_views_spread_model_flies']) for c in ('on', 'off') if c in g.index}
    return st, al


STAGE_SPREAD_KEYS: Tuple[str, ...] = ('stage', *LOSS_STAGES)


def prediction(lock: Mapping, groups: pd.DataFrame) -> dict:
    """The prediction per megamat stage, written by step 0 before any read (a rule on step 0's APL16 numbers and the
    stage's place; the lock's last choice). Its three classed questions:
    - 'a fair stand-in': the size of the population's change at the APL16 somata with APL on (E85's 72, |population -
      the wirings' mean view at the lab's input|) is taken as the size at every megamat stage, with its sign; a stage
      whose spread as scored is larger reads 'a fair stand-in', else 'moves it' in that direction;
    - the input-fly term: at the PN dendrites and boutons "as large as the wiring's" (the wirings differ there only by
      the seeds for pooling and noise, and at the boutons by APL's feedback alone); at the claws and somata the APL16
      somata's ratio with APL on, by the classes' thresholds;
    - the model's flies: the APL16 somata's model fly against the others with APL on, by the classes' thresholds, at
      every stage (the family's PN stage compresses the odors' differences: its spread of the odors' totals is a quarter
      of the flies', E90 and E120, so input differences are not expected to grow on the way);
    and A1's fall: "stands through the input flies" (its fall on T at the lab's input, 0.115, is eleven of its standard
    errors, and the input-fly term of a paired change is expected to be well under the APL16 term of an error)."""
    g = groups[(groups.group == 'B1') & (groups.condition == 'on')].iloc[0]
    delta = float(g['strength_population'] - g['strength_lab_common'])
    ratio = float(g['strength_term_ratio'])
    oth = float(g['strength_against_others_mean'])
    spreads = lock['spreads_as_scored']
    stages = {}
    for k in STAGE_SPREAD_KEYS:
        sp = float(spreads[k])
        stand = 'a fair stand-in' if abs(delta) <= sp else ('moves it nearer' if delta < 0 else 'moves it further')
        term = "as large as the wiring's" if k in ('pn_dendrite', 'pn_bouton') else term_class(ratio)
        stages['megamat_stage_error' if k == 'stage' else k] = dict(
            stand_in=stand, input_fly_term=term, model_flies=flies_class(oth), spread_as_scored=sp)
    return dict(experiment=EXPERIMENT, written=_now(), host=_host(), before_any_read=True,
                basis=dict(apl16_on_population_minus_lab=delta, apl16_on_term_ratio=ratio,
                           apl16_on_against_others=oth, source="step 0: E85's 72 'repro' held model flies at the APL16 "
                           "somata, read by this module"),
                rule=prediction.__doc__, stages=stages, a1_fall='stands through the input flies')


def part_step0() -> dict:
    """Step 0 (the module docstring; no simulation). Stops (RuntimeError) if an identity fails. Writes
    step0_identities.json, step0_input_flies.csv, step0_orn_stage.csv, step0_apl16_groups.csv,
    step0_apl16_model_flies.csv, step0_apl16_realism.csv, step0_probe.csv, step0_classes_recomputed.csv,
    classes_step0.json, prediction.json (once: an existing one is kept), tables_step0.md and provenance_step0.json."""
    import time
    t0 = time.time()
    lock = require_lock()
    check_files(lock)
    ids = {}
    ids['e134'] = identity_e134(lock)
    ids['e85_part_a'] = identity_e85(lock)
    ids['e140_aa'] = identity_e140(lock)
    ids['apl16_reading'] = identity_apl16(lock)
    tables = class_tables()
    same = classes_with_spreads({}, {}, tables)
    ids['classes_as_scored'] = dict(ok=bool((same.cls_as_scored == same.cls_scored).all()),
                                    differ=same[same.cls_as_scored != same.cls_scored][['table', 'arm']].to_dict('records'))
    lab_apl16 = ids['apl16_reading'].pop('_sides')
    for k in ('_sides', '_held'):
        ids['e134'].pop(k, None)
    _write_json(results_dir() / 'step0_identities.json', dict(time=_now(), **ids))
    failed = [k for k, v in ids.items() if not v.get('ok')]
    print(json.dumps({k: v.get('ok') for k, v in ids.items()}), flush=True)
    if failed:
        raise RuntimeError(f'step 0: identities failed: {failed} (step0_identities.json)')
    # the input flies
    _write_csv(input_flies_table(), 'step0_input_flies.csv')
    _write_csv(orn_stage_table(), 'step0_orn_stage.csv')
    # E85's saved model flies at the APL16 somata
    rt, runs = e85_runs()
    spreads = lock['spreads_as_scored']
    grows, per = [], []
    for g, gt in rt.groupby('group', sort=False):
        sub = {t: runs[t] for t in gt.tag}
        sides = read_sides(sub, ['kc_soma'], 'apl16')
        r, pf = apl16_group_rows(g, sides, lab_apl16, gt.set_index('tag').held, spreads)
        grows += r
        per.append(pf)
        print(f'step 0: {g} read ({len(gt)} model flies, {time.time() - t0:.0f} s)', flush=True)
    groups = pd.DataFrame(grows)
    _write_csv(groups, 'step0_apl16_groups.csv')
    _write_csv(pd.concat(per, ignore_index=True), 'step0_apl16_model_flies.csv')
    _write_csv(apl16_realism_rows(rt, runs), 'step0_apl16_realism.csv')
    _write_csv(probe_numbers(lock), 'step0_probe.csv')
    # the classes recomputed with the APL16 spreads over model flies
    st, al = model_fly_spreads_apl16(groups)
    rc = classes_with_spreads(st, al, tables)
    _write_csv(rc, 'step0_classes_recomputed.csv')
    b1 = groups[groups.group == 'B1'].set_index('condition')
    classes = dict(time=_now(), spreads_as_scored={k: spreads[k] for k in APL16_KEYS},
                   spreads_model_flies=st, spreads_model_flies_all_views=al,
                   apl16={c: dict(stand_in=b1.loc[c, 'stand_in'], input_fly_term=b1.loc[c, 'input_fly_term'],
                                  model_flies=b1.loc[c, 'model_flies'], n_unheld=int(b1.loc[c, 'n_unheld']))
                          for c in ('on', 'off')},
                   recomputed_changed=rc[rc.changed][['table', 'arm', 'cls_scored', 'cls_new']].to_dict('records'))
    _write_json(results_dir() / 'classes_step0.json', classes)
    pp = results_dir() / 'prediction.json'
    if not pp.exists():
        _write_json(pp, prediction(lock, groups))
    write_tables_step0(ids, groups, rc)
    provenance('step0', lock, seconds=round(time.time() - t0, 1), n_e85_runs=len(rt),
               inputs={k: v['sha256'] for k, v in lock['files'].items()})
    print(f'step 0 done in {time.time() - t0:.0f} s', flush=True)
    return dict(identities=ids, groups=groups, classes=classes)


def _md(df: pd.DataFrame, floatfmt: str = '.3f') -> str:
    cols = list(df.columns)
    out = ['| ' + ' | '.join(str(c) for c in cols) + ' |', '|' + '|'.join('---:' if pd.api.types.is_numeric_dtype(df[c])
                                                                         else '---' for c in cols) + '|']
    for r in df.itertuples(index=False):
        out.append('| ' + ' | '.join(format(v, floatfmt) if isinstance(v, (float, np.floating)) and np.isfinite(v)
                                     else ('–' if isinstance(v, (float, np.floating)) else str(v)) for v in r) + ' |')
    return '\n'.join(out)


def write_tables_step0(ids: Mapping, groups: pd.DataFrame, rc: pd.DataFrame) -> None:
    """tables_step0.md: the identities, the APL16 somata's model flies by group and the recomputed classes."""
    pred = json.loads((results_dir() / 'prediction.json').read_text())
    g = groups.copy()
    show = g[['group', 'condition', 'n_model_flies', 'n_unheld', 'strength_population', 'strength_lab_common',
              'spread_as_scored', 'stand_in', 'strength_mean_model_fly', 'strength_lab_mean_wiring',
              'strength_term_input_fly', 'strength_term_wiring', 'strength_term_remainder', 'input_fly_term',
              'strength_against_others_mean', 'model_flies', 'strength_spread_model_flies', 'strength_spread_wirings_lab']]
    lines = ['# E147 step 0: the megamat stage error read on model flies with each ORN fly\'s own input', '',
             f'Generated by {SCRIPT} at {H.code_commits()["queue_commit"][:7]} ({_now()}). No simulation.', '',
             '## Identities', '',
             f"- E134's megamat stage error, the family held, all flies: mean of the eight wirings {ids['e134']['model']:.6f}, "
             f"their mean view {ids['e134']['common']:.6f} (E134: {ids['e134']['expected']['model']:.6f} / "
             f"{ids['e134']['expected']['common']:.6f}); E136's held-out half {ids['e134']['held_half']['common']:.6f} "
             f"(the same as E136's score_megamat.csv to {ids['e134']['held_half']['common_diff_from_e136']:.1e}).",
             f"- E85's part A on the 12 shared odors from E134's copy: slope {ids['e85_part_a']['slope']:.4f} (E85 "
             f"{ids['e85_part_a']['slope_e85']:.4f}), r_g to {ids['e85_part_a']['r_g_max_diff']:.1e}, the nine APL16 inputs "
             f"to {max(ids['e85_part_a']['inputs_max_abs_diff'].values()):.1e}.",
             f"- E140's aa vectors: E140's slope {ids['e140_aa']['slope_e140']:.4f} (this map's 17-odor regression "
             f"{ids['e140_aa']['slope17']:.4f}), the seven flies given aa, every vector to "
             f"{max(ids['e140_aa']['max_rel_diff'].values()):.1e}.",
             f"- The APL16 somata reading: E136's family rows to {max(v['each'] for v in ids['apl16_reading']['diffs'].values()):.1e}; "
             f"every class of E136, E139 and E141 recomputed with its spreads as scored is the scored class: "
             f"{ids['classes_as_scored']['ok']}.", '',
             "## E85's saved model flies at the APL16 somata (strength by odor, one-fly units)", '',
             _md(show), '',
             '## Every class recomputed with the APL16 spreads over model flies (B1), as a description', '',
             _md(rc[['table', 'arm', 'cls_scored', 'cls_new', 'changed']]), '',
             '## The prediction (prediction.json, written before any read)', '',
             f"Written {pred['written']} on {pred['host']}.", '',
             _md(pd.DataFrame([dict(stage=k, **v) for k, v in pred['stages'].items()])), '',
             f"A1's fall: {pred['a1_fall']}.", '']
    (results_dir() / 'tables_step0.md').write_text('\n'.join(lines))


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
