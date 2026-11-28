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
  (E134's default readings, seed index the wiring's, E136's `Flies` for the held-out half; `read_sides_multi`);
- the population's error: the error of the mean view over the model flies; a model fly's error: its own; both against
  all the stage's flies, E136's held-out half beside;
- the terms: the SD over input flies of the mean over wirings (the input-fly term), the SD over wirings of the mean
  over input flies (the wiring term), and the remainder (`terms`);
- a model fly against the others: its view against the mean view of the model flies with another wiring and another
  input fly, in the stage's one-fly units (`against_others`);
- the spread over model flies: the SD over all model flies of the quantity whose SD over the eight wirings the classes
  use (`model_fly_spreads`).

Points and tags (`Point`, `point_tag`, `parse_tag`, `planned_points`): '<setting>_<panel>_<inputfly>_<wiring>', e.g.
'family_megamat_f3_P4', with setting one of SETTINGS ('family': the unrefitted family at each wiring's held megamat
point, E136's FAMILY; 'A1': E136's A1 at its classed point, `E136.fitted_z('A1')`; 'familyfree' and 'A1free': the free
reading, each model fly at its wiring's held parameters with the input fly's input, nothing re-held), panel 'megamat' or
'apl16', input fly f0 ... f8, wiring one of E136's WIRINGS. Step 1: 'family' on megamat, the nine input flies on all
eight wirings, the T wirings first (72). Step 2, on T: 'A1' on megamat (36), then 'A1' and 'family' on APL16 (36 each).
The free reading ('familyfree', 'A1free', one run per model fly of steps 1 and 2) is beside, classing nothing, and is
the runner's too. Step '6a' (the lock's revision STEP6A_REVISION; beside, classing nothing): the family on megamat
re-held as in step 1 with the input flies' deviations as measured ('familymeas', E85's 'measured') and their
reproducible deviations permuted over the glomeruli ('familyperm', E85's 'permuted'), T first and within a wiring
'familymeas' then 'familyperm' (144); `setting_parts` and SETTING_VARIANT give the runner each setting's re-held setting
and its input's variant (`input_fly(fly, 'megamat', SETTING_VARIANT[setting])`), and its runs go to
`run_dirs('6a', ...)` (folder step6a).

What hal's runner writes, which `--part score` reads:
- each read's run in E134's step-3 layout (`kcq.stage_readings.save_step3_run`, every stage saved, the units' metadata,
  the input as `<tag>_orn.parquet`) under `run_dir(step, setting, wiring, fly)` =
  outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/<fly>/ (local, step1, step2 or step6a; `run_dirs` also
  gives the archive's data/model_outputs/E147_megamat_input_flies/... equivalent); a read whose `<tag>_run.json` is
  there is done;
- its row in results/E147_megamat_input_flies/reads_<lane>.csv (any name reads_*.csv; the last row per tag counts,
  `reads_table`), with the columns of READ_COLUMNS: tag, step, setting, panel, fly, wiring; the held parameters
  (PARAM_KEYS, absolute, and off_<key> of the APL-off run); thr_factor, apl_factor (against the point the hold started
  from, which the row's 'start' column names); thr_vs_own, apl_vs_own (hal's runner: against the setting's own reading at
  the lab's input, which the realism tables give beside); n_on, n_off, target_on, target_off, frac_on_target, frac_off_target; at_bound_thr, at_bound_apl;
  held (the runner's own reading of its hold); rate_on, rate_off,
  max_count_on, max_count_off; input_sha256 (E77's `frame_sha256` of the input, which `--part score` checks against
  `input_fly`); hold_n_sims, seconds, host, lane, time, threads and the commits (`H.code_commits()`). The score decides
  "held" itself (`is_held`: both counts within HOLD_TOL, 0.5 %, of their targets and no factor at a bound), so the
  row's own 'held' is informative only; where a row lacks n_on or n_off the score counts them from the run's saved
  responses (`fill_counts`). A read that failed to hold writes its row and its run as held so far, so that the score
  counts it unheld;
- the free reading's points ('familyfree', 'A1free') in the same layout under the step of their held twin; the score
  reads them beside, with no class;
- E134's reading of the run beside its counts is optional: `--part score` reads every run itself
  (`read_sides_multi`: E134's default readings, seed index the wiring's, all the flies and E136's held-out half).

Parts (CLI `--part`):
- lock (atlas): experiments/E147_lock.json and results/.../lock.json; commit the copy before any statistic is read;
- step0 (atlas, no simulation): the identities (E134's 0.918 and 0.894; E85's part A on the 12 shared odors; E140's aa
  vectors), the table of megamat's input flies on its 17 odors, E85's 72 saved model flies and its bounds on S scored
  at the APL16 somata as E134 and E136 score a run (seed index the wiring's; the planner's probe, E134's archive
  scoring at seed 0, beside in step0_probe.csv), every class of E136, E139 and E141 recomputed with the
  APL16 spreads over model flies, and prediction.json (a prediction per megamat stage, timestamped, before any run);
- lock6a (atlas): the lock's step-6a revision appended to the committed lock before any step-6a read
  (`part_lock_step6a`); commit both copies;
- step1, step2, step6a: hal's (experiments/E147_reads.py --part reads --step 1|2|6a); here they only say so;
- score (atlas, no simulation; run once hal's reads are in): the classes per megamat stage, for the megamat stage error
  and at the APL16 somata; A1's megamat fall and its APL16 changes model fly by model fly against the family on the
  same wiring and input fly; the megamat spreads over model flies and every class of E136, E139 and E141 recomputed
  with them, as a description; the realism bands per model fly and per population; the beside items the reads allow;
  where step-6a reads exist, its two populations beside the reproducible deviations' on the same wirings and the
  lock's three readings per stage (`score_step6a`: score_step6a_*.csv, classes.json 'step6a', a section of tables.md),
  gating nothing; without them nothing of step 6a is written.

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
    ./venv/bin/python -u experiments/E147_megamat_input_flies.py --part lock6a   # step 6a's revision; commit it likewise
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
N_INPUT_FLIES = 9
FLY_LABELS: Tuple[str, ...] = tuple(f'f{i}' for i in range(N_INPUT_FLIES))
# the APL16 odors whose ORN recordings are other flies': they keep the lab's input in every model fly
APL16_KEPT_AT_LAB: Tuple[str, ...] = ('paa', 'euc', '-aPine', 'B-myr')
DEVIATION = 'repro'
# step 6a's two inputs (E85's 'measured' and 'permuted', on megamat's 17 odors; `megamat_deviations`)
VARIANT_MEASURED = 'measured'
VARIANT_PERMUTED = 'permuted'
DEVIATION_VARIANTS: Tuple[str, ...] = (DEVIATION, VARIANT_MEASURED, VARIANT_PERMUTED, 'zero')
STEP6A = '6a'
STEP6A_SETTINGS: Tuple[str, ...] = ('familymeas', 'familyperm')
SETTINGS: Tuple[str, ...] = ('family', 'A1', 'familyfree', 'A1free', *STEP6A_SETTINGS)
HELD_SETTINGS: Tuple[str, ...] = ('family', 'A1')
# per setting: the setting whose parameters and hold it reads, whether nothing is re-held, and its input flies' variant
SETTING_BASE: Dict[str, str] = {'family': 'family', 'A1': 'A1', 'familyfree': 'family', 'A1free': 'A1',
                                'familymeas': 'family', 'familyperm': 'family'}
SETTING_VARIANT: Dict[str, str] = {'family': DEVIATION, 'A1': DEVIATION, 'familyfree': DEVIATION, 'A1free': DEVIATION,
                                   'familymeas': VARIANT_MEASURED, 'familyperm': VARIANT_PERMUTED}
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

# folders in place of the defaults (tests): 'results', 'local', 'archive'
DIRS: Dict[str, Path] = {}


def results_dir() -> Path:
    """<shared>/results/E147_megamat_input_flies/ (created)."""
    d = Path(DIRS['results']) if 'results' in DIRS else H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    """outputs/E147_megamat_input_flies/[<part>/] in this repository (created)."""
    root = Path(DIRS['local']) if 'local' in DIRS else H.outputs_dir(EXPERIMENT)
    d = root / part if part else root
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E147_megamat_input_flies/ (not created)."""
    return Path(DIRS['archive']) if 'archive' in DIRS else H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


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
    """flies x megamat's 17 odors (`orn_flies()['odors']` order) x 36 glomeruli deviations, 0 where the fly did not
    measure an entry, with a sum over the flies of 0 for every entry:
    - 'repro' (DEVIATION): r_g (`reliabilities`) times the fly's trial mean minus the mean over the flies that have the
      entry;
    - 'measured' (VARIANT_MEASURED, E85's 'measured'): the trial mean minus that mean, with no reliability factor;
    - 'permuted' (VARIANT_PERMUTED, E85's 'permuted'): 'repro' permuted over each fly's measured glomeruli within each
      odor and re-centred per entry over the flies that have it (E85's `permute_deviations`, seeded by the fly's
      index, and `recenter`), the odors in this order, so that on the 12 shared odors the draws are not E85's;
    - 'zero': all 0."""
    from experiments import E85_orn_fly_population as E85m
    F = orn_flies() if F is None else F
    if variant == 'zero':
        return np.zeros(F['x'].shape[:3])
    if variant not in (DEVIATION, VARIANT_MEASURED, VARIANT_PERMUTED):
        raise ValueError(f'{variant=}')
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        st = E85m.trial_stats(F['x'])
    if variant == VARIANT_MEASURED:
        return st['D']
    d = reliabilities(F)[None, None, :] * st['D']
    if variant == VARIANT_PERMUTED:
        return E85m.recenter(E85m.permute_deviations(d, st['measured']), st['measured'])
    return d


def deviations(panel: str, variant: str = DEVIATION, F: Optional[dict] = None) -> np.ndarray:
    """flies x odors x glomeruli deviations to add to `lab_input(panel)` (its odor and glomerulus order; 0 where a fly
    did not measure an entry, and on the APL16 odors of APL16_KEPT_AT_LAB); `variant` one of DEVIATION_VARIANTS
    (`megamat_deviations`)."""
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
    dtypes) plus the fly's deviation of `variant` (the reproducible one by default; a setting's is
    SETTING_VARIANT[setting]) where measured, the lab's value elsewhere; on APL16 paa, euc, -aPine and B-myr at the
    lab's input. `fly` is a label (f0 ... f8) or an E134 key. Pass it to E134's `build_family(..., orn=...)` and E77's
    `hold_rates` in place of the panel's input."""
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


def step_key(step) -> object:
    """A step as the points carry it: STEP6A ('6a') for '6a' or 'step6a', else the int (1, 2; from 1, '1', ...)."""
    if str(step) in (STEP6A, f'step{STEP6A}'):
        return STEP6A
    return int(step)


def step_folder(step) -> str:
    """A step's folder name: 'step1', 'step2', 'step6a'."""
    return f'step{step_key(step)}'


def setting_parts(setting: str) -> Tuple[str, bool, str]:
    """(the re-held setting whose parameters and hold the reads use, whether the reading is free, the input flies'
    variant to pass to `input_fly`) of a setting of SETTINGS: 'familymeas' is ('family', False, 'measured'),
    'familyperm' ('family', False, 'permuted'), 'A1free' ('A1', True, 'repro')."""
    if setting not in SETTINGS:
        raise ValueError(f'{setting=} not in {SETTINGS}')
    return SETTING_BASE[setting], setting.endswith('free'), SETTING_VARIANT[setting]


def planned_points(step, free: bool = False) -> List[Point]:
    """The reads of a step in the order a forward lane takes them: step 1 'family' on megamat, the T wirings first,
    then S, fly by fly within a wiring (72); step 2 on T 'A1' on megamat (36), then 'A1' on APL16 (36), then 'family'
    on APL16 (36); step '6a' (STEP6A) on megamat the T wirings first, then S, within a wiring 'familymeas' (the
    deviations as measured) then 'familyperm' (permuted), fly by fly (144; Point.step '6a'). `free`: the same points
    under the free reading ('familyfree', 'A1free'; none for step 6a)."""
    fl = list(FLY_LABELS)
    step = step_key(step)
    if step == 1:
        pts = [Point(1, 'family', 'megamat', f, w) for w in T_WIRINGS + S_WIRINGS for f in fl]
    elif step == 2:
        pts = [Point(2, s, p, f, w) for s, p in (('A1', 'megamat'), ('A1', 'apl16'), ('family', 'apl16'))
               for w in T_WIRINGS for f in fl]
    elif step == STEP6A:
        if free:
            raise ValueError('step 6a has no free reading')
        pts = [Point(STEP6A, s, 'megamat', f, w) for w in T_WIRINGS + S_WIRINGS for s in STEP6A_SETTINGS for f in fl]
    else:
        raise ValueError(f'{step=}')
    if free:
        pts = [pt._replace(setting=f'{pt.setting}free') for pt in pts]
    return pts


def step_of(setting: str, panel: str) -> object:
    """The step a point belongs to: '6a' for STEP6A_SETTINGS, 1 for the family on megamat, 2 otherwise (the free
    reading as its held twin)."""
    if setting in STEP6A_SETTINGS:
        return STEP6A
    base = setting[:-4] if setting.endswith('free') else setting
    return 1 if (base, panel) == ('family', 'megamat') else 2


def run_dir(step, setting: str, wiring: str, fly: str) -> Path:
    """Where a read's run is saved locally: outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/<fly>/
    (`step_folder`: step1, step2, step6a; created)."""
    return local_dir(f'{step_folder(step)}/{setting}/{wiring}/{fly}')


def run_dirs(step, setting: str, wiring: str, fly: str) -> List[Path]:
    """`run_dir` and its archive twin (data/model_outputs/E147_megamat_input_flies/step<step>/<setting>/<wiring>/
    <fly>/), in that order, for `kcq.stage_readings.save_step3_run` and `load_step3_run`."""
    return [run_dir(step, setting, wiring, fly), archive_dir() / step_folder(step) / setting / wiring / fly]


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
    """The model flies `runs` ({tag: (wiring, input fly, SavedRun or a function returning it)}) read at `stages` on
    `panel` (E134's default readings, seed index the wiring's, `_fly_cache(half)`), one `Side` per (stage, panel,
    condition) with flies: the reference E134's `fly_reference_for` on all the stage's flies (every part), or for a half
    E136's `Flies.reference` (every part, the half's mean in all flies' units)."""
    return read_sides_multi(runs, stages, panel, (half,), conditions)[half]


def read_sides_multi(runs: Mapping[str, Tuple[str, str, object]], stages: Sequence[str], panel: str,
                     halves: Sequence[Optional[str]] = (None,), conditions: Sequence[str] = ('on', 'off')
                     ) -> Dict[Optional[str], Dict[tuple, Side]]:
    """`read_sides` for several readings of the flies at once ({half: {key: Side}}; None all the flies), each run
    loaded once (a function in `runs` is called and its run dropped after reading)."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as E136m
    from kcq import population_fit as P
    from kcq import stage_readings as R
    lock = R.require_lock()
    fcs = {h: _fly_cache(h, stages) for h in halves}
    pops: Dict[Optional[str], Dict[tuple, list]] = {h: {} for h in halves}
    for tag, (w, f, run) in runs.items():
        run = run() if callable(run) else run
        for h in halves:
            cache = fcs[h] if h is None else fcs[h].cache
            rd = R.read_run(run, stages=list(stages), lock=lock, flies=cache, seed_index=E136m.seed_index(w))
            for k, pop in rd.populations.items():
                if k[1] == panel and k[2] in conditions and (h is None or k[2] == 'on'):
                    pops[h].setdefault(k, []).append((tag, w, f, pop))
            del rd
        del run
    out: Dict[Optional[str], Dict[tuple, Side]] = {}
    for h in halves:
        out[h] = {}
        for k, items in pops[h].items():
            odors = list(items[0][3].values.columns)
            if h is None:
                fl = fcs[h].get(k)
                if fl is None:
                    continue
                _, ref, _ = B.fly_reference_for(fl, odors)
            else:
                ref = fcs[h].reference(k[0], odors)
            arrays = P.stack_views([P.reference_views(it[3], ref) for it in items], ref)
            out[h][k] = Side(k, ref, [it[0] for it in items], [it[1] for it in items], [it[2] for it in items], arrays)
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


REVISION_UNHELD = dict(
    time='2026-11-27 09:45', what="'unheld' read with E136's hold tolerance (HOLD_TOL, 0.5 % of each target) in place "
    "of exact counts", why="step 0's first run (09:39) found that E77's hold_rates, which E85's runs and E136's reads "
    "use, misses a target by one or two counts on 12 of E85's 72 'repro' runs (at most 0.1 % of ~1900 APL-on and ~4100 "
    "APL-off responses), as a bisection on discrete counts can; the entry's 'unheld' means a hold that failed, not a "
    "count step. No step 1 or 2 read existed. The first lock is commit 2a58895")
STEP6A_REVISION = dict(
    time='2026-11-27 20:37', what="step 6a added (the lock's 'step6a', points.step6a and settings 'familymeas' / "
    "'familyperm'): the family on megamat re-held as in step 1 with the nine input flies' deviations as measured (E85's "
    "'measured') and their reproducible deviations permuted over the glomeruli (E85's 'permuted'), 144 reads, T first; "
    "per stage and for the megamat stage error the readings \"the structure's\", \"the variance's\" and \"(e) stands at "
    "the upper bound\", beside, classing nothing", why="the planner's word on design item 6 (viki 2026-11-27 20:16; "
    "planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-27 19:55): the score's outcomes (c) and "
    "(e) (20:07) lack E85's two controls on megamat. No step-6a read existed (the runner's option was not written). "
    "Steps 0 to 2, their points, tags, settings, classes and rules are unchanged")
REVISIONS: Tuple[dict, ...] = (REVISION_UNHELD, STEP6A_REVISION)

STEP6A_LOCK = dict(
    step=STEP6A, folder=f'step{STEP6A}',
    planned_by="the planner's word on design item 6 (viki 2026-11-27 20:16; a Fable subagent, claude-fable-5-2, "
               "spawned from viki 2026-11-27 19:55)",
    settings={s: dict(base=SETTING_BASE[s], variant=SETTING_VARIANT[s], free=False) for s in STEP6A_SETTINGS},
    reads="the family on megamat (E136's FAMILY at each wiring's held megamat point), both response rates re-held as in "
          "step 1 (E134's held targets), the breadth the lab input's, the wiring's seed for pooling and noise; the "
          "input `input_fly(fly, 'megamat', SETTING_VARIANT[setting])`; every stage saved as in step 1 under "
          "run_dirs('6a', setting, wiring, fly) (folder step6a); the row's step '6a' and, where the runner writes one, "
          "its 'variant' column the setting's variant (the score checks it with the input's hash)",
    variants={
        VARIANT_MEASURED: "E85's 'measured' on megamat's 17 odors: each fly's trial mean minus the mean over the flies "
                          "that have the entry, 0 where unmeasured, no reliability factor (on the 12 shared odors E85's "
                          "exactly); the upper bound on the deviations' size",
        VARIANT_PERMUTED: "E85's 'permuted' on megamat's 17 odors: the reproducible deviation (`reliabilities` times "
                          "the above) permuted over each fly's measured glomeruli within each odor by E85's "
                          "`permute_deviations` (numpy default_rng seeded by the fly's index 0 ... 8, the odors in "
                          "megamat's order, `orn_flies()['odors']`, so that the draws on the 12 shared odors are not "
                          "E85's) and re-centred per entry over the flies that have it by E85's `recenter`; the "
                          "deviations without their glomerulus structure"},
    order="the T wirings first (P3, P4, M2, M3), then S (real, P1, P2, M1); within a wiring 'familymeas' then "
          "'familyperm', fly by fly f0 ... f8 (144 reads: 36 and 36 on T, then 36 and 36 on S)",
    n_points=144,
    readings=dict(
        population="per megamat stage (strength by odor) and for the megamat stage error, against all the stage's "
                   "flies (E136's held-out half beside), each of M's and P's populations by step 1's rule: its "
                   "population's error against the error of the wirings' mean view at the lab's input on the same "
                   "wirings (E134's held megamat runs), 'a fair stand-in' if the difference is within the spread as "
                   "scored (the lock's spreads_as_scored: 'stage' and each stage's), else 'moves it nearer' or 'moves "
                   "it further' by its sign; the reproducible deviations' population (R: step 1's held family reads) "
                   "on the same wirings beside",
        wiring_sets="read on T once a variant's 36 T reads are saved and on all eight once its 72 are (both "
                    "variants for the readings); a population on other wirings is reported with no reading",
        structure="\"the structure's\": at a stage where R moves it nearer, P is a fair stand-in",
        variance="\"the variance's\": at a stage where R moves it nearer, P moves it nearer too",
        neither="at a stage where R moves it nearer and P moves it further, neither; where R is a fair stand-in or moves "
                "it further, neither is read (M's and P's classes are given)",
        e_upper_bound="\"(e) stands at the upper bound\": at a stage where R's model flies vary too little (step 1's "
                      "reading: a model fly against the others, the mean over the model flies, under 0.75), it stands "
                      "if under M a model fly against the others (the same mean over M's model flies, on the same "
                      "wirings) stays under 0.75 (FLIES_LOW), and does not stand otherwise",
        unheld="step 1's: unheld model flies left out and counted; no reading for a population with more than a ninth "
               "of its model flies unheld"),
    gates="nothing: a description beside the classes; no class as scored, no reading of the stage table and no "
          "yardstick changes; which reading is a stage's under 'like for like' stays with Tom and Yang; no fit, search "
          "or weighting",
)


STEP6A_SETTING_TEXT: Dict[str, str] = {
    'familymeas': "step 6a: the family on megamat re-held as 'family' (step 1), the input flies' deviations as measured "
                  "(E85's 'measured')",
    'familyperm': "step 6a: the family on megamat re-held as 'family' (step 1), the input flies' reproducible "
                  "deviations permuted over the glomeruli (E85's 'permuted')"}


def step6a_lock_points() -> List[str]:
    """Step 6a's tags in `planned_points('6a')` order (the lock's points.step6a)."""
    return [point_tag(p.setting, p.panel, p.fly, p.wiring) for p in planned_points(STEP6A)]


def step6a_started() -> List[str]:
    """What of step 6a exists already: its tags among the reads' rows and the claims, and its run folders (local and
    the archive's) holding a saved run."""
    tags = set(step6a_lock_points())
    found = sorted(tags & set(reads_table().tag.astype(str)))
    claims = results_dir() / 'claims'
    if claims.exists():
        found += sorted(f'claim {d.name}' for d in claims.iterdir() if d.name in tags)
    for root in (local_dir(), archive_dir()):
        d = Path(root) / step_folder(STEP6A)
        if d.exists() and any(d.rglob('*_run.json')):
            found.append(str(d))
    return found


def part_lock_step6a() -> dict:
    """The lock's step-6a revision (STEP6A_REVISION), before any step-6a read: the committed lock (LOCK_COPY) with
    points.step6a, settings 'familymeas' and 'familyperm', 'step6a' (STEP6A_LOCK) and the revision added, everything
    else as committed; written to LOCK_COPY and results/.../lock.json (commit the copy before any step-6a read).
    PermissionError if a step-6a read, claim or run exists; a lock that has the revision already is left as it is."""
    lock = json.loads(LOCK_COPY.read_text())
    if any(r.get('time') == STEP6A_REVISION['time'] for r in lock.get('revisions', [])):
        print(f'{LOCK_COPY.name} has the step-6a revision already', flush=True)
        return lock
    if [dict(r) for r in lock.get('revisions', [])] != [dict(REVISION_UNHELD)]:
        raise RuntimeError(f"{LOCK_COPY.name}: its revisions are not this module's REVISIONS before step 6a")
    found = step6a_started()
    if found:
        raise PermissionError(f'step 6a has begun ({found[:5]}): its lock revision comes before any read')
    lock['points']['step6a'] = step6a_lock_points()
    lock['settings'].update(STEP6A_SETTING_TEXT)
    lock['step6a'] = json.loads(json.dumps(STEP6A_LOCK))
    lock['revisions'].append(dict(STEP6A_REVISION))
    txt = json.dumps(lock, indent=1, default=float) + '\n'
    LOCK_COPY.write_text(txt)
    lock_path().write_text(txt)
    print(f'wrote {LOCK_COPY} and {lock_path()} (revision {STEP6A_REVISION["time"]}); commit '
          f'{LOCK_COPY.relative_to(H.HERE)} before any step-6a read', flush=True)
    return lock


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
                      free="'familyfree' / 'A1free': each model fly at its wiring's held parameters, nothing re-held",
                      **STEP6A_SETTING_TEXT),
        points=dict(tag='<setting>_<panel>_<inputfly>_<wiring>',
                    step1=[point_tag(p.setting, p.panel, p.fly, p.wiring) for p in planned_points(1)],
                    step2=[point_tag(p.setting, p.panel, p.fly, p.wiring) for p in planned_points(2)],
                    order='T wirings first (step 1: T then S; step 2 on T only), fly by fly within a wiring',
                    step6a=step6a_lock_points()),
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
        revisions=[dict(r) for r in REVISIONS],
        step6a=STEP6A_LOCK,
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


# ---------------------------------------------------------------------------
# the score

READINGS: Tuple[Optional[str], ...] = (None, 'held')   # all the flies (the primary reading), E136's held-out half
DOWNSTREAM: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')


def reading_name(half: Optional[str]) -> str:
    return 'all' if half is None else half


def points_from_rows(rows: pd.DataFrame) -> pd.DataFrame:
    """The reads' rows with setting, panel, fly, wiring and step parsed from their tags (rows whose tag does not parse
    are dropped), 'held' by `is_held` (True for the free reading, which is not held) and 'group' '<setting>/<panel>';
    the last row per tag counts."""
    if rows.empty:
        return pd.DataFrame(columns=['tag', 'setting', 'panel', 'fly', 'wiring', 'step', 'held', 'group'])
    out = []
    for r in rows.drop_duplicates('tag', keep='last').to_dict('records'):
        try:
            st, pn, f, w = parse_tag(str(r['tag']))
        except ValueError:
            continue
        free = st.endswith('free')
        out.append(dict(r, setting=st, panel=pn, fly=f, wiring=w, step=step_of(st, pn),
                        held=True if free else is_held(r), free=free, group=f'{st}/{pn}'))
    return pd.DataFrame(out)


def unheld_rule(n_unheld: int, n_total: int) -> bool:
    """True when a population has no class: more than UNHELD_MAX_SHARE of its model flies unheld."""
    return bool(n_total > 0 and n_unheld > UNHELD_MAX_SHARE * n_total)


def paired(a: pd.Series, b: pd.Series) -> pd.DataFrame:
    """Per (wiring, input fly) the change a - b of two per-model-fly Series indexed by (wiring, fly), both present."""
    j = pd.concat([a.rename('a'), b.rename('b')], axis=1).dropna()
    j['change'] = j['a'] - j['b']
    return j


def a1_fall(changes: pd.DataFrame) -> dict:
    """A1's fall through the input flies from `paired` changes (A1 minus the family on the same wiring and input
    fly): the mean change, the standard error (the larger of the SD over input flies of the per-fly means over sqrt of
    their number and the SD over wirings of the per-wiring means over sqrt of theirs), the fall (minus the mean) and
    whether it stands (beyond IMPROVE_SE standard errors and beyond CONTROL_FALL)."""
    if changes.empty:
        return dict(n=0, stands=None)
    ch = changes['change']
    by_f = ch.groupby(level=1).mean()
    by_w = ch.groupby(level=0).mean()
    se_f = float(by_f.std(ddof=1) / np.sqrt(len(by_f))) if len(by_f) > 1 else np.nan
    se_w = float(by_w.std(ddof=1) / np.sqrt(len(by_w))) if len(by_w) > 1 else np.nan
    se = float(np.nanmax([se_f, se_w])) if np.isfinite([se_f, se_w]).any() else np.nan
    mean = float(ch.mean())
    fall = -mean
    return dict(n=int(len(ch)), n_wirings=int(len(by_w)), n_flies=int(len(by_f)), mean_change=mean, fall=fall,
                se=se, se_over_input_flies=se_f, se_over_wirings=se_w, se_multiples=fall / se if se > 0 else np.nan,
                control_fall=CONTROL_FALL, n_model_flies_falling=int((ch < 0).sum()),
                per_fly={str(k): float(v) for k, v in by_f.items()}, per_wiring={str(k): float(v) for k, v in by_w.items()},
                stands=bool(np.isfinite(se) and fall > IMPROVE_SE * se and fall > CONTROL_FALL))


def a1_apl16_shares(changes: pd.DataFrame, spread: float) -> dict:
    """A1 at the APL16 somata: per input fly the mean over wirings of A1's change, and the share of input flies on
    which it is worse than the spread as scored (above it), better than it (below minus it)."""
    if changes.empty:
        return dict(n_flies=0)
    by_f = changes['change'].groupby(level=1).mean()
    return dict(n_flies=int(len(by_f)), n_model_flies=int(len(changes)), mean_change=float(changes['change'].mean()),
                spread=float(spread), share_worse=float((by_f > spread).mean()), share_better=float((by_f < -spread).mean()),
                per_fly={str(k): float(v) for k, v in by_f.items()})


def per_fly_series(values: pd.Series, side: Side) -> pd.Series:
    """A per-model-fly Series (index the side's tags) re-indexed by (wiring, input fly)."""
    return pd.Series(values.reindex(side.tags).to_numpy(dtype=float),
                     index=pd.MultiIndex.from_arrays([side.wiring, side.fly], names=['wiring', 'fly']))


def megamat_values(sides: Mapping[str, Side]) -> Dict[str, pd.Series]:
    """Per model fly ((wiring, input fly) index): the megamat stage error ('stage'), the error after the PN dendrites
    ('downstream'), each stage's strength error, and each stage's error on every view ('guard/<stage>')."""
    out = {}
    each = {st: side_scores(sd) for st, sd in sides.items()}
    for st, sd in sides.items():
        out[st] = per_fly_series(each[st]['view_strength'], sd)
        out[f'guard/{st}'] = per_fly_series(each[st][stage_column(st, None)], sd)
    stg = [st for st in LOSS_STAGES if st in out]
    df = pd.concat({st: out[st] for st in stg}, axis=1)
    out['stage'] = np.sqrt((df ** 2).mean(axis=1, skipna=False))
    dn = [st for st in DOWNSTREAM if st in out]
    out['downstream'] = np.sqrt((df[dn] ** 2).mean(axis=1, skipna=False)) if dn else pd.Series(dtype=float)
    return out


def population_rows(name: str, sides: Mapping[str, Side], lab: Optional[Mapping[str, Side]], spreads: Mapping[str, float],
                    reading: str, n_unheld: int = 0, n_missing: int = 0, classed: bool = True) -> List[dict]:
    """One megamat population read one way: per stage and for the megamat stage error, the population's error against
    the wirings' mean view at the lab's input on the same wirings, each model fly's error (mean, min, max), the terms,
    a model fly against the others, the spread over model flies and the classes (none where `classed` is False, a
    reference is missing, or more than a ninth of the model flies are unheld); every view's error of the population
    beside."""
    rows = []
    total = len(next(iter(sides.values()))) + n_unheld if sides else 0
    no_class = unheld_rule(n_unheld, total) or not classed
    vals = megamat_values(sides)
    pops = {st: population_score(sd) for st, sd in sides.items()}
    oth = {st: against_others(sd) for st, sd in sides.items()}
    labpops = {st: population_score(sd) for st, sd in lab.items()} if lab else {}
    ws = sorted(set(next(iter(sides.values())).wiring), key=WIRINGS.index) if sides else []
    for st in [*[x for x in LOSS_STAGES if x in sides], 'megamat_stage_error']:
        if st == 'megamat_stage_error':
            pv = float(np.sqrt(np.mean([float(pops[x]['view_strength']) ** 2 for x in sides])))
            lv = float(np.sqrt(np.mean([float(labpops[x]['view_strength']) ** 2 for x in sides]))) if labpops else np.nan
            each = vals['stage']
            o = np.sqrt(pd.concat({x: per_fly_series(oth[x]['view_strength'], sides[x]) for x in sides}, axis=1).pow(2)
                        .mean(axis=1, skipna=False))
            spread = spreads.get('stage', np.nan)
            extra = {}
        else:
            pv = float(pops[st]['view_strength'])
            lv = float(labpops[st]['view_strength']) if st in labpops else np.nan
            each = vals[st]
            o = per_fly_series(oth[st]['view_strength'], sides[st])
            spread = spreads.get(st, np.nan)
            extra = {f'pop_{c}': float(pops[st][c]) for c in pops[st].index
                     if c.startswith('view_') or c in ('error', 'error_with_matrix')}
            extra.update({f'labpop_{c}': float(labpops[st][c]) for c in (labpops[st].index if st in labpops else [])
                          if c.startswith('view_') or c in ('error', 'error_with_matrix')})
        T = each.unstack('fly').reindex(index=[w for w in WIRINGS if w in ws], columns=[f for f in FLY_LABELS])
        T = T.dropna(axis=1, how='all')
        tm = terms(T)
        rec = dict(population=name, reading=reading, stage=st, n_model_flies=int(each.notna().sum()), n_unheld=n_unheld,
                   n_missing=n_missing, wirings=' '.join(ws), population_error=pv, lab_common=lv,
                   population_minus_lab=pv - lv, spread_as_scored=float(spread),
                   mean_model_fly=float(each.mean()), min_model_fly=float(each.min()), max_model_fly=float(each.max()),
                   spread_model_flies=float(each.std(ddof=1)) if each.notna().sum() > 1 else np.nan,
                   term_input_fly=tm['input_fly'], term_wiring=tm['wiring'], term_remainder=tm['remainder'],
                   term_ratio=tm['ratio'], terms_complete=tm['complete'], against_others_mean=float(o.mean()),
                   against_others_min=float(o.min()), against_others_max=float(o.max()), no_class=no_class, **extra)
        rec.update(stand_in='no class' if no_class or not np.isfinite(lv) else stand_in_class(pv, lv, spread),
                   input_fly_term='no class' if no_class or not tm['complete'] else term_class(tm['ratio']),
                   model_flies='no class' if no_class else flies_class(rec['against_others_mean']))
        rows.append(rec)
    return rows


def model_fly_spreads(vals: Mapping[str, pd.Series]) -> Tuple[Dict[str, float], Dict[str, float]]:
    """(the spreads keyed as E136's: 'stage', 'downstream', each stage; the guard's: each stage on every view), the SD
    (ddof 1) over the model flies of `megamat_values`."""
    st = {k: float(vals[k].std(ddof=1)) for k in ('stage', 'downstream', *LOSS_STAGES) if k in vals and len(vals[k]) > 1}
    gd = {k: float(vals[f'guard/{k}'].std(ddof=1)) for k in LOSS_STAGES if f'guard/{k}' in vals}
    return st, gd


def odor_entry_rows(name: str, sides: Mapping[str, Side]) -> pd.DataFrame:
    """Each odor's strength entries by input fly at every stage (leads): the mean over the population's wirings of each
    input fly's model flies' entry, the flies' mean, and their difference in the part's one-fly units; and the
    population's (all its model flies) beside as 'all'."""
    rows = []
    for st, sd in sides.items():
        for q in [q for q in ('strength/rate', 'strength/mean') if q in sd.arrays]:
            m = sd.ref.mean[q]
            u = float(sd.ref.one_fly[q])
            for f in [*FLY_LABELS, 'all']:
                mask = np.ones(len(sd), bool) if f == 'all' else (sd.fly == f)
                if not mask.any():
                    continue
                a = _nanmean0(sd.arrays[q][mask])
                for o, av, mv in zip(m.index, a, m.to_numpy(dtype=float)):
                    rows.append(dict(population=name, stage=st, part=q, fly=f, odor=o, model=float(av), flies=float(mv),
                                     diff_units=float((av - mv) / u)))
    return pd.DataFrame(rows)


def without_odor_rows(name: str, side: Side, odor: str = 'aa') -> pd.DataFrame:
    """The strength error of each model fly and of the population with and without one odor's entry (E138's
    `odor_entries` rule, as E136's `stage_without_odor` reads it: the odor's term out of each strength part's squared
    error, the number of odors kept, nothing renormalized)."""
    parts = [q for q in ('strength/rate', 'strength/mean') if q in side.arrays]
    labels = [*side.tags, 'population']
    with_, without = np.zeros(len(labels)), np.zeros(len(labels))
    for q in parts:
        m = side.ref.mean[q]
        if odor not in list(m.index):
            raise KeyError(f'{odor!r} not among the odors of {side.key}')
        A = np.vstack([side.arrays[q], _nanmean0(side.arrays[q])[None, :]])
        e = (A - m.to_numpy(dtype=float)[None, :]) / float(side.ref.one_fly[q])
        ok = np.isfinite(e)
        n = np.maximum(ok.sum(axis=1), 1)
        e2 = np.where(ok, e ** 2, 0.0)
        j = list(m.index).index(odor)
        with_ += e2.sum(axis=1) / n
        without += (e2.sum(axis=1) - e2[:, j]) / n
    k = len(parts)
    return pd.DataFrame(dict(population=name, stage=side.key[0], odor=odor, row=labels,
                             view_error=np.sqrt(with_ / k), view_error_without=np.sqrt(without / k)))


def megamat_realism_rows(pts: pd.DataFrame, runs_info: Mapping[str, dict]) -> pd.DataFrame:
    """Per megamat model fly: the response rates with APL on and off beside the flies' 0.0915 (E134's lock), the
    threshold and APL factors after the re-hold (thr_vs_own / apl_vs_own: against the setting's own reading at the lab's
    input; thr_factor / apl_factor: against the hold's start), the held parameters and the largest spike counts (the
    reads' rows, and the runs where the rows lack them)."""
    rows = []
    for r in pts.to_dict('records'):
        info = runs_info.get(r['tag'], {})
        rows.append(dict(tag=r['tag'], group=r['group'], wiring=r['wiring'], fly=r['fly'], held=r['held'],
                         rate_on=r.get('rate_on', info.get('rate_on')), rate_off=r.get('rate_off', info.get('rate_off')),
                         flies_rate=0.0915, thr_vs_own=r.get('thr_vs_own'), apl_vs_own=r.get('apl_vs_own'),
                         thr_factor=r.get('thr_factor'), apl_factor=r.get('apl_factor'), fixed_thr=r.get('fixed_thr'), wAPLKC_scale=r.get('wAPLKC_scale'),
                         max_count_on=r.get('max_count_on', info.get('max_count_on')),
                         max_count_off=r.get('max_count_off', info.get('max_count_off'))))
    return pd.DataFrame(rows)


def _loader(tag: str, step: int, info: dict):
    """A function loading a read (`load_point`) that records its rates and largest counts into `info`."""
    def load():
        run = load_point(tag, step)
        if run is None:
            raise FileNotFoundError(tag)
        rec = {}
        for c in ('on', 'off'):
            st = run.conditions[c]
            rec[f'rate_{c}'] = float((st['responses'].to_numpy(dtype=float) > 0.5).mean())
            rec[f'max_count_{c}'] = float(st['spike_counts'].to_numpy(dtype=float).max())
        info[tag] = rec
        return run
    return load


def fill_counts(pts: pd.DataFrame) -> pd.DataFrame:
    """The points with n_on and n_off (the runs' response counts, from their saved responses) where the rows lack them,
    and 'held' re-read by `is_held` for those."""
    pts = pts.copy()
    for c in ('n_on', 'n_off'):
        if c not in pts:
            pts[c] = np.nan
    for i, r in pts.iterrows():
        if np.isfinite(pd.to_numeric(r.get('n_on'), errors='coerce')) and \
                np.isfinite(pd.to_numeric(r.get('n_off'), errors='coerce')):
            continue
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            for d in run_dirs(r['step'], r['setting'], r['wiring'], r['fly']):
                f = Path(d) / f"{r['tag']}_{cc}_responses.parquet"
                if f.exists():
                    pts.at[i, f'n_{c}'] = float((pd.read_parquet(f).to_numpy(dtype=float) > 0.5).sum())
                    break
        if not r['free']:
            pts.at[i, 'held'] = is_held(pts.loc[i].to_dict())
    return pts


def available(pts: pd.DataFrame) -> pd.DataFrame:
    """The points whose run is saved (locally or in the archive)."""
    keep = []
    for r in pts.itertuples():
        d = run_dirs(r.step, r.setting, r.wiring, r.fly)
        keep.append(any((Path(x) / f'{r.tag}_run.json').exists() for x in d))
    return pts[np.asarray(keep, dtype=bool)] if len(pts) else pts


def a1_lab_runs(panel: str) -> Dict[str, Tuple[str, str, object]]:
    """E136's A1 reads at its classed point on the lab's input (megamat on T, read against the held-out half; APL16 on
    T), as model flies of the input fly 'lab'; missing reads are left out."""
    from experiments import E136_stage_population_fit as E136m
    before = E136m.MEASURE
    E136m.set_measure(MEASURE)
    try:
        out = {}
        for w in T_WIRINGS:
            r = E136m.load_read(E136m.ReadPoint('A1', panel, w, 'held' if panel == 'megamat' else None))
            if r is not None:
                out[f'lab_A1_{panel}_{w}'] = (w, 'lab', r)
    finally:
        E136m.set_measure(before)       # E136's measure is a module global: leave it as it was
    return out


INPUT_ATOL = 1e-9  # Hz: a saved input within this of `input_fly` is the same input (the hash is bitwise)


def saved_input(tag: str, step, setting: str, wiring: str, fly: str) -> Optional[pd.DataFrame]:
    """A read's saved ORN input (`<tag>_orn.parquet` in `run_dirs`), or None."""
    for d in run_dirs(step, setting, wiring, fly):
        f = Path(d) / f'{tag}_orn.parquet'
        if f.exists():
            return pd.read_parquet(f)
    return None


def check_inputs_detail(pts: pd.DataFrame) -> Tuple[List[str], Dict[str, float]]:
    """The reads whose input is not `input_fly` of their fly, panel and setting's variant (SETTING_VARIANT), and the
    largest difference of the others.

    A row whose input_sha256 equals the sha256 of `input_fly` passes. The hash is of the floats' bytes, so the same
    input built on another host can differ in its last bits (E85's map is a least-squares fit, whose rounding follows
    the host's BLAS): such a row is compared value by value with its saved `<tag>_orn.parquet` (by position: the saved
    frame keeps the glomeruli and odors in `input_fly`'s order, its columns flattened to the odor names). It passes
    when every entry is within INPUT_ATOL. A row with no hash match and no saved input, or a difference beyond
    INPUT_ATOL, is bad; so is a row whose 'variant' column (the runner's, where it writes one) is not its setting's.

    Returns:
        (bad tags, {'n_bitwise', 'n_within_atol', 'max_abs_diff'}).
    """
    from experiments.E77_input_screen import frame_sha256
    frames: Dict[tuple, pd.DataFrame] = {}
    want: Dict[tuple, str] = {}
    bad, n_bit, n_tol, worst = [], 0, 0, 0.0
    for r in pts.to_dict('records'):
        variant = SETTING_VARIANT.get(r['setting'], DEVIATION)
        if isinstance(r.get('variant'), str) and r['variant'] != variant:
            bad.append(r['tag'])
            continue
        key = (r['fly'], r['panel'], variant)
        if key not in frames:
            frames[key] = input_fly(r['fly'], r['panel'], variant)
            want[key] = frame_sha256(frames[key])
        if r.get('input_sha256') == want[key]:
            n_bit += 1
            continue
        x = saved_input(r['tag'], r['step'], r['setting'], r['wiring'], r['fly'])
        w = frames[key]
        if x is None or x.shape != w.shape or [str(i) for i in x.index] != [str(i) for i in w.index]:
            bad.append(r['tag'])
            continue
        diff = float(np.max(np.abs(x.to_numpy(dtype=float) - w.to_numpy(dtype=float))))
        if not np.isfinite(diff) or diff > INPUT_ATOL:
            bad.append(r['tag'])
        else:
            n_tol += 1
            worst = max(worst, diff)
    return bad, dict(n_bitwise=n_bit, n_within_atol=n_tol, max_abs_diff=worst)


def check_inputs(pts: pd.DataFrame) -> List[str]:
    """The tags of `check_inputs_detail`'s bad reads: those whose input is neither `input_fly`'s bitwise (by hash) nor
    within INPUT_ATOL of it (by their saved input)."""
    return check_inputs_detail(pts)[0]


def part_score(rows: Optional[pd.DataFrame] = None) -> dict:
    """The score of hal's reads (the module docstring; no simulation), on whatever reads are saved: writes
    score_populations.csv, score_model_flies.csv, score_a1.json, score_apl16_groups.csv, score_apl16_a1.csv,
    score_realism_megamat.csv, score_realism_apl16.csv, score_classes_recomputed.csv, score_odor_entries.csv,
    score_without_aa.csv, score_e85_beside.csv, classes.json, tables.md and provenance_score.json."""
    import time
    t0 = time.time()
    lock = require_lock()
    check_files(lock)
    spreads = lock['spreads_as_scored']
    pts = available(points_from_rows(reads_table() if rows is None else rows))
    if not pts.empty:
        pts = fill_counts(pts)
    if pts.empty:
        print('score: no reads yet', flush=True)
        return {}
    bad_inputs, input_check = check_inputs_detail(pts)
    info: Dict[str, dict] = {}
    res: Dict[str, object] = dict(time=_now(), n_reads=int(len(pts)), bad_inputs=bad_inputs, input_check=input_check,
                                  groups={g: int(n) for g, n in pts.group.value_counts().items()})
    poprows, flyrows, mega_vals = [], [], {}
    # megamat
    lab_cache: Dict[str, dict] = {}

    def lab_fam() -> dict:
        if 'fam' not in lab_cache:
            lab_cache['fam'] = read_sides_multi(e134_held_runs(), LOSS_STAGES, 'megamat', READINGS, ('on',))
        return lab_cache['fam']

    repro_sides: Dict[str, Dict[str, Side]] = {}
    for group in ('family/megamat', 'A1/megamat', 'familyfree/megamat', 'A1free/megamat'):
        g = pts[pts.group == group]
        if g.empty:
            continue
        setting = group.split('/')[0]
        unheld = int((~g.held.astype(bool)).sum())
        use = g[g.held.astype(bool)]
        runs = {r.tag: (r.wiring, r.fly, _loader(r.tag, r.step, info)) for r in use.itertuples()}
        sides = read_sides_multi(runs, LOSS_STAGES, 'megamat', READINGS, ('on',))
        base = setting[:-4] if setting.endswith('free') else setting
        ws = sorted(set(use.wiring), key=WIRINGS.index)
        if group == 'family/megamat':
            repro_sides = {reading_name(h): {k[0]: v for k, v in sides[h].items()} for h in READINGS}
        if base == 'family':
            lab = {h: {k: v.subset([w in ws for w in v.wiring]) for k, v in lab_fam()[h].items()} for h in READINGS}
        else:
            lab = read_sides_multi(a1_lab_runs('megamat'), LOSS_STAGES, 'megamat', READINGS, ('on',))
        for h in READINGS:
            sd = {k[0]: v for k, v in sides[h].items()}
            lb = {k[0]: v for k, v in lab[h].items()}
            name = group
            poprows += population_rows(name, sd, lb, spreads, reading_name(h), unheld, 0, classed=not setting.endswith('free'))
            if base == 'family' and set(ws) > set(T_WIRINGS):
                tmask = {k: v.subset([w in T_WIRINGS for w in v.wiring]) for k, v in sd.items()}
                lmask = {k: v.subset([w in T_WIRINGS for w in v.wiring]) for k, v in lb.items()}
                poprows += population_rows(f'{name} (T)', tmask, lmask, spreads, reading_name(h), unheld, 0,
                                           classed=not setting.endswith('free'))
            vals = megamat_values(sd)
            mega_vals[(group, reading_name(h))] = vals
            for k, v in vals.items():
                for (w, f), x in v.items():
                    flyrows.append(dict(population=group, reading=reading_name(h), quantity=k, wiring=w, fly=f,
                                        value=float(x)))
            if h is None:
                odor = odor_entry_rows(group, sd)
                _append_frame(res, 'odor_entries', odor)
                if 'kc_soma' in sd:
                    _append_frame(res, 'without_aa', without_odor_rows(group, sd['kc_soma']))
        print(f'score: {group} read ({len(use)} held, {unheld} unheld; {time.time() - t0:.0f} s)', flush=True)
    pop = pd.DataFrame(poprows)
    if len(pop):
        _write_csv(pop, 'score_populations.csv')
    if flyrows:
        _write_csv(pd.DataFrame(flyrows), 'score_model_flies.csv')
    for k, fname in (('odor_entries', 'score_odor_entries.csv'), ('without_aa', 'score_without_aa.csv')):
        if k in res:
            _write_csv(res.pop(k), fname)
    # step 6a: the deviations as measured and permuted, beside (nothing without its reads)
    s6 = score_step6a(pts, info, spreads, repro_sides, lab_fam, pts[pts.group == 'family/megamat'])
    if s6 is not None:
        _write_csv(s6['populations'], 'score_step6a_populations.csv')
        _write_csv(s6['readings_table'], 'score_step6a_readings.csv')
        if len(s6['model_flies']):
            _write_csv(s6['model_flies'], 'score_step6a_model_flies.csv')
        print(f"score: step 6a read ({s6['n_reads']}; {time.time() - t0:.0f} s)", flush=True)
    # A1's fall, model fly by model fly
    a1 = {}
    for rd in ('held', 'all'):
        A, F = mega_vals.get(('A1/megamat', rd)), mega_vals.get(('family/megamat', rd))
        if A is not None and F is not None:
            a1[rd] = dict(stage=a1_fall(paired(A['stage'], F['stage'])),
                          **{st: a1_fall(paired(A[st], F[st])) for st in LOSS_STAGES if st in A and st in F})
    if a1:
        _write_json(results_dir() / 'score_a1.json', a1)
    # the APL16 somata
    apl_rows, apl_per, apl_sides = [], [], {}
    for group in ('family/apl16', 'A1/apl16', 'familyfree/apl16', 'A1free/apl16'):
        g = pts[pts.group == group]
        if g.empty:
            continue
        setting = group.split('/')[0]
        base = setting[:-4] if setting.endswith('free') else setting
        runs = {r.tag: (r.wiring, r.fly, _loader(r.tag, r.step, info)) for r in g.itertuples()}
        sides = read_sides(runs, ['kc_soma'], 'apl16')
        apl_sides[group] = sides
        ws = sorted(set(g.wiring), key=WIRINGS.index)
        lab = read_sides(lab_apl16_runs(ws) if base == 'family' else a1_lab_runs('apl16'), ['kc_soma'], 'apl16')
        r_, pf = apl16_group_rows(group, sides, lab, g.set_index('tag').held, spreads)
        if setting.endswith('free'):
            for q in r_:
                q.update(stand_in='beside', input_fly_term='beside', model_flies='beside')
        apl_rows += r_
        apl_per.append(pf)
    groups = pd.DataFrame(apl_rows)
    if len(groups):
        _write_csv(groups, 'score_apl16_groups.csv')
        _write_csv(pd.concat(apl_per, ignore_index=True), 'score_apl16_model_flies.csv')
    a1_apl = {}
    if 'A1/apl16' in apl_sides and 'family/apl16' in apl_sides:
        for c in ('on', 'off'):
            k = ('kc_soma', 'apl16', c)
            sa, sf = apl_sides['A1/apl16'][k], apl_sides['family/apl16'][k]
            ch = paired(per_fly_series(side_scores(sa)['view_strength'], sa),
                        per_fly_series(side_scores(sf)['view_strength'], sf))
            a1_apl[c] = a1_apl16_shares(ch, float(spreads[f'apl16/apl16/{c}']))
        _write_json(results_dir() / 'score_apl16_a1.json', a1_apl)
    # realism
    _write_csv(megamat_realism_rows(pts[pts.panel == 'megamat'], info), 'score_realism_megamat.csv')
    if (pts.panel == 'apl16').any():
        _write_csv(apl16_realism_score(pts[pts.panel == 'apl16']), 'score_realism_apl16.csv')
    # E85's saved model flies beside step 2's family APL16 reads
    if 'family/apl16' in apl_sides:
        _write_csv(e85_beside(pts[pts.group == 'family/apl16'], apl_sides['family/apl16']), 'score_e85_beside.csv')
    # the classes recomputed with the spreads over model flies (a description)
    recomputed = {}
    fam = mega_vals.get(('family/megamat', 'held'))
    if fam is not None:
        mst, mgd = model_fly_spreads(fam)
        s0 = json.loads((results_dir() / 'classes_step0.json').read_text()) if (results_dir() / 'classes_step0.json').exists() else {}
        variants = {'a_megamat': ({}, {})}
        if s0:
            variants['b_megamat_apl16_e85'] = (s0.get('spreads_model_flies', {}), s0.get('spreads_model_flies_all_views', {}))
        if len(groups) and (groups.group == 'family/apl16').any():
            gg = groups[groups.group == 'family/apl16'].set_index('condition')
            variants['c_megamat_apl16_step2'] = (
                {f'apl16/apl16/{c}': float(gg.loc[c, 'strength_spread_model_flies']) for c in gg.index},
                {f'apl16/apl16/{c}': float(gg.loc[c, 'all_views_spread_model_flies']) for c in gg.index})
        tabs = class_tables()
        parts = []
        for vname, (ast, agd) in variants.items():
            t = classes_with_spreads(dict(mst, **ast), dict(mgd, **agd), tabs).assign(variant=vname)
            parts.append(t)
            recomputed[vname] = t[t.changed][['table', 'arm', 'cls_scored', 'cls_new']].to_dict('records')
        _write_csv(pd.concat(parts, ignore_index=True), 'score_classes_recomputed.csv')
        res['spreads_model_flies_megamat'] = dict(strength=mst, guard=mgd)
    # the classes
    classes = dict(time=_now(), n_reads=int(len(pts)), megamat={}, apl16={}, a1_fall={}, a1_apl16=a1_apl,
                   recomputed_changed=recomputed, bad_inputs=bad_inputs, input_check=input_check)
    if len(pop):
        for r in pop[pop.reading == 'all'].itertuples():
            classes['megamat'].setdefault(r.population, {})[r.stage] = dict(
                stand_in=r.stand_in, input_fly_term=r.input_fly_term, model_flies=r.model_flies, n_unheld=int(r.n_unheld))
    if len(groups):
        for r in groups.itertuples():
            classes['apl16'].setdefault(r.group, {})[r.condition] = dict(
                stand_in=r.stand_in, input_fly_term=r.input_fly_term, model_flies=r.model_flies, n_unheld=int(r.n_unheld))
    for rd, v in a1.items():
        classes['a1_fall'][rd] = dict(stands=v['stage'].get('stands'), fall=v['stage'].get('fall'), se=v['stage'].get('se'))
    classes['outcomes'] = outcomes(classes, pop)
    if s6 is not None:
        classes['step6a'] = s6['classes']
    _write_json(results_dir() / 'classes.json', dict(classes, **{k: v for k, v in res.items() if k not in classes}))
    write_tables_score(pop, groups, a1, a1_apl, recomputed, classes, s6)
    provenance('score', lock, seconds=round(time.time() - t0, 1), n_reads=int(len(pts)),
               reads=sorted(pts.tag), inputs={k: v['sha256'] for k, v in lock['files'].items()})
    print(f'score done in {time.time() - t0:.0f} s', flush=True)
    return classes


def _append_frame(res: dict, key: str, df: pd.DataFrame) -> None:
    res[key] = pd.concat([res[key], df], ignore_index=True) if key in res else df


def apl16_realism_score(pts: pd.DataFrame) -> pd.DataFrame:
    """The four realism bands per APL16 model fly (E136's `apl16_realism`) and per population (E85's
    `population_realism` over its model flies), per group, with the largest counts and the factors."""
    from experiments import E136_stage_population_fit as E136m
    from experiments import E85_orn_fly_population as E85m
    rows, parts = [], {}
    for r in pts.to_dict('records'):
        run = load_point(r['tag'], r['step'])
        if run is None:
            continue
        on, off = run.conditions['on'], run.conditions['off']
        rows.append(dict(group=r['group'], unit=r['tag'], kind='model_fly', wiring=r['wiring'], fly=r['fly'],
                         held=r['held'], n_model_flies=1, thr_vs_own=r.get('thr_vs_own'), apl_vs_own=r.get('apl_vs_own'),
                         thr_factor=r.get('thr_factor'), apl_factor=r.get('apl_factor'),
                         max_count_on=float(on['spike_counts'].to_numpy(dtype=float).max()),
                         max_count_off=float(off['spike_counts'].to_numpy(dtype=float).max()),
                         **E136m.apl16_realism(on['responses'], off['responses'])))
        parts.setdefault(r['group'], []).append(E85m.run_realism_parts({'on': (on['responses'], on['spike_counts']),
                                                                         'off': (off['responses'], off['spike_counts'])}))
    for g, ps in parts.items():
        pr = E85m.population_realism(ps)
        rows.append(dict(group=g, unit='population', kind='population', n_model_flies=len(ps), rate_on=pr['rate_on'],
                         rate_off=pr['rate_off'], rate_ratio=pr['ratio'], rate_rho_on=pr['rate_rho_on'],
                         rate_rho_off=pr['rate_rho_off'], opk_on_tv=pr['opk_on_tv'], opk_off_tv=pr['opk_off_tv'],
                         realism_flag=pr['realism_flag']))
    return pd.DataFrame(rows)


def e85_beside(pts: pd.DataFrame, sides: Mapping[tuple, Side]) -> pd.DataFrame:
    """Step 2's family APL16 reads (breadth the lab input's) beside E85's saved B1T model flies (breadth recomputed per
    input) on the same wiring and input fly: the entries whose spike counts differ, the response counts, and the
    strength error of each at the APL16 somata."""
    from experiments import E85_orn_fly_population as E85m
    rows = []
    e85_scores = pd.read_csv(results_dir() / 'step0_apl16_model_flies.csv') \
        if (results_dir() / 'step0_apl16_model_flies.csv').exists() else pd.DataFrame()
    for r in pts.to_dict('records'):
        tag85 = f"B1T_repro_{r['fly']}_{r['wiring']}" if r['wiring'] in T_WIRINGS else f"B1S_repro_{r['fly']}_{r['wiring']}"
        o = E85m.load_run('B', tag85)
        run = load_point(r['tag'], r['step'])
        if o is None or run is None:
            continue
        rec = dict(tag=r['tag'], e85_tag=tag85, wiring=r['wiring'], fly=r['fly'])
        for c, cc in (('on', 'on'), ('off', 'off')):
            a = H.abbrev_columns(run.conditions[c]['spike_counts'].astype(float))
            b = H.abbrev_columns(o[cc][1].astype(float))[list(a.columns)]
            rec[f'differing_counts_{c}'] = int((a.to_numpy() != b.to_numpy()).sum())
            rec[f'responses_{c}'] = int((run.conditions[c]['responses'].to_numpy(dtype=float) > 0.5).sum())
            rec[f'e85_responses_{c}'] = int((o[cc][0].to_numpy(dtype=float) > 0.5).sum())
            sd = sides[('kc_soma', 'apl16', c)]
            rec[f'strength_{c}'] = float(side_scores(sd)['view_strength'].get(r['tag'], np.nan))
            if len(e85_scores):
                m = e85_scores[(e85_scores.group == 'B1') & (e85_scores.condition == c) & (e85_scores.tag == tag85)]
                rec[f'e85_strength_{c}'] = float(m.strength.iloc[0]) if len(m) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def outcomes(classes: Mapping, pop: pd.DataFrame) -> dict:
    """The entry's outcomes read off the classes (a description; nothing follows by itself): (a) a fair stand-in at
    every megamat stage and the megamat stage error with a small input-fly term, (b) a fair stand-in with an input-fly
    term as large as the wiring's, (c) the population reading moves a megamat stage, (d) A1's fall does not stand,
    (e) the model's flies vary too much or too little at a stage, (f) a population without a class for its holds."""
    fam = (classes.get('megamat') or {}).get('family/megamat') or {}
    out: Dict[str, object] = {}
    if fam:
        stands = [v['stand_in'] for v in fam.values()]
        tterm = [v['input_fly_term'] for v in fam.values()]
        out['a'] = bool(all(x == 'a fair stand-in' for x in stands) and all(x == 'small' for x in tterm))
        out['b'] = bool(all(x == 'a fair stand-in' for x in stands) and any(x.startswith('as large') for x in tterm))
        out['c'] = sorted(k for k, v in fam.items() if v['stand_in'].startswith('moves'))
        out['e'] = {k: v['model_flies'] for k, v in fam.items() if v['model_flies'] in ('vary too little', 'vary too much')}
    a1 = (classes.get('a1_fall') or {}).get('held')
    if a1 and a1.get('stands') is not None:
        out['d'] = not a1['stands']
    out['f'] = sorted({p for p, st in (classes.get('megamat') or {}).items() for v in st.values()
                       if v['stand_in'] == 'no class' and v['n_unheld'] > 0}
                      | {g for g, st in (classes.get('apl16') or {}).items() for v in st.values()
                         if v['stand_in'] == 'no class' and v['n_unheld'] > 0})
    return out


# ---------------------------------------------------------------------------
# step 6a: the deviations as measured and permuted, beside (the lock's revision STEP6A_REVISION)

STEP6A_SETS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (('T', T_WIRINGS), ('all', WIRINGS))
STEP6A_POPULATIONS: Dict[str, str] = {DEVIATION: 'family/megamat', VARIANT_MEASURED: 'familymeas/megamat',
                                      VARIANT_PERMUTED: 'familyperm/megamat'}
STEP6A_DESCRIPTION = (
    f"step 6a, beside, classing nothing (the lock's revision of {STEP6A_REVISION['time']}): the family on megamat "
    "re-held as in step 1, the nine input flies with the deviations as measured (M, 'familymeas') and the reproducible ones "
    "permuted over the glomeruli (P, 'familyperm'); each population read as step 1's ('a fair stand-in' against the "
    "wirings' mean view at the lab's input on the same wirings, the spread as scored), beside the reproducible "
    "deviations' population (R, step 1's family reads) on the same wirings; read on T once both variants' 36 T reads "
    "are in and on all eight once their 72 are. The readings gate nothing: no class as scored, no reading of the stage "
    "table and no yardstick changes")


def step6a_sets(g: pd.DataFrame, setting: str) -> List[Tuple[str, Tuple[str, ...]]]:
    """The wiring sets a step-6a population is read on: ('T', T_WIRINGS) once every one of its 36 T points is saved,
    ('all', WIRINGS) once its 72 are; ('partial', the wirings with a saved read) when neither, read with no reading."""
    have = set(g.tag.astype(str))
    out = [(name, ws) for name, ws in STEP6A_SETS
           if all(point_tag(setting, 'megamat', f, w) in have for w in ws for f in FLY_LABELS)]
    if not out:
        ws = tuple(w for w in WIRINGS if w in set(g.wiring))
        if ws:
            out.append(('partial', ws))
    return out


def subset_sides(sides: Mapping[str, Side], wirings: Sequence[str]) -> Dict[str, Side]:
    """Each stage's Side restricted to the model flies on `wirings`."""
    return {st: sd.subset([w in wirings for w in sd.wiring]) for st, sd in sides.items()}


def _num(x) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def step6a_reading(r: Optional[Mapping], m: Optional[Mapping], p: Optional[Mapping]) -> dict:
    """The lock's three step-6a readings at one stage (or the megamat stage error), from the `population_rows` rows,
    all the flies, of the reproducible deviations' population (R, step 1's family reads), the deviations as measured
    (M) and permuted (P), all on the same wirings (None where a population is not read):
    - 'structure_or_variance': where R moves the stage nearer, "the structure's" if P is a fair stand-in, "the
      variance's" if P moves it nearer too, 'neither: P moves it further' otherwise; 'not read' where R is a fair
      stand-in or moves it further; 'no reading' where R or P has no class;
    - 'e_applies': whether R's model flies vary too little there (a model fly against the others under FLIES_LOW);
      'e_at_upper_bound': where it applies, whether under M a model fly against the others (the mean over M's model
      flies) stays under FLIES_LOW; None where it does not apply or M has no class.
    Each population's class, its model flies' class and its values beside."""
    def get(x, k, default='no class'):
        return (x or {}).get(k, default)

    rc, mc, pc = get(r, 'stand_in'), get(m, 'stand_in'), get(p, 'stand_in')
    if rc == 'no class' or pc == 'no class':
        sv = 'no reading'
    elif rc != 'moves it nearer':
        sv = 'not read'
    elif pc == 'a fair stand-in':
        sv = "the structure's"
    elif pc == 'moves it nearer':
        sv = "the variance's"
    else:
        sv = 'neither: P moves it further'
    e_applies = get(r, 'model_flies') == 'vary too little'
    m_oth = _num(get(m, 'against_others_mean', None))
    e = None if not e_applies or get(m, 'model_flies') == 'no class' or m_oth is None else bool(m_oth < FLIES_LOW)
    out = dict(structure_or_variance=sv, e_applies=bool(e_applies), e_at_upper_bound=e,
               spread_as_scored=_num(get(r or m or p, 'spread_as_scored', None)))
    for name, x in (('repro', r), ('measured', m), ('permuted', p)):
        out[f'{name}_stand_in'] = get(x, 'stand_in')
        out[f'{name}_model_flies'] = get(x, 'model_flies')
        for k in ('population_error', 'lab_common', 'population_minus_lab', 'against_others_mean', 'term_input_fly',
                  'term_wiring', 'n_model_flies', 'n_unheld'):
            out[f'{name}_{k}'] = _num(get(x, k, None))
    return out


def score_step6a(pts: pd.DataFrame, info: dict, spreads: Mapping[str, float], repro: Mapping[str, Mapping[str, Side]],
                 lab_fam, repro_pts: pd.DataFrame) -> Optional[dict]:
    """Step 6a's two populations beside the reproducible deviations' (no simulation): None without a step-6a read.

    Args:
        pts: the score's points (`available`, `fill_counts`).
        info: `_loader`'s record of the runs' rates and counts.
        spreads: the spreads as scored (the lock's).
        repro: step 1's held family reads on megamat, {reading name: {stage: Side}} (empty without them).
        lab_fam: a function giving E134's held megamat runs read as `read_sides_multi` gives them, {half: {key: Side}}.
        repro_pts: step 1's family points on megamat (for their unheld and missing counts per wiring set).

    Returns:
        dict(populations: `population_rows` of R, M and P per wiring set (`step6a_sets`) and reading, with 'setting',
        'variant' and 'wirings_set'; model_flies: each M and P model fly's `megamat_values`; readings_table: per
        complete wiring set ('T', 'all') and stage, `step6a_reading` (all the flies) and beside it on E136's held-out
        half ('reading' 'held'); classes: what classes.json keeps under 'step6a'; n_reads, n_unheld per setting).
    """
    six = pts[pts.setting.isin(STEP6A_SETTINGS)]
    if six.empty:
        return None
    lab = lab_fam()
    poprows, flyrows, wsets = [], [], []
    n_reads, n_unheld = {}, {}
    for setting in STEP6A_SETTINGS:
        group, variant = f'{setting}/megamat', SETTING_VARIANT[setting]
        g = six[six.setting == setting]
        n_reads[setting], n_unheld[setting] = int(len(g)), int((~g.held.astype(bool)).sum())
        if g.empty:
            continue
        use = g[g.held.astype(bool)]
        runs = {r.tag: (r.wiring, r.fly, _loader(r.tag, r.step, info)) for r in use.itertuples()}
        sides = read_sides_multi(runs, LOSS_STAGES, 'megamat', READINGS, ('on',)) if runs else {}
        for h in READINGS:
            sd = {k[0]: v for k, v in (sides.get(h) or {}).items()}
            if not sd:
                continue
            for k, v in megamat_values(sd).items():
                for (w, f), x in v.items():
                    flyrows.append(dict(population=group, variant=variant, reading=reading_name(h), quantity=k,
                                        wiring=w, fly=f, value=float(x)))
        for wset, ws in step6a_sets(g, setting):
            if (wset, ws) not in wsets:
                wsets.append((wset, ws))
            gs = g[g.wiring.isin(ws)]
            nu = int((~gs.held.astype(bool)).sum())
            for h in READINGS:
                sd = subset_sides({k[0]: v for k, v in (sides.get(h) or {}).items()}, ws)
                if not sd or any(len(v) == 0 for v in sd.values()):
                    continue
                lb = subset_sides({k[0]: v for k, v in lab[h].items()}, ws)
                rows = population_rows(group, sd, lb, spreads, reading_name(h), nu, len(ws) * N_INPUT_FLIES - len(gs),
                                       classed=wset != 'partial')
                poprows += [dict(r, setting=setting, variant=variant, wirings_set=wset) for r in rows]
    for wset, ws in wsets:
        rp = repro_pts[repro_pts.wiring.isin(ws)]
        nu = int((~rp.held.astype(bool)).sum()) if len(rp) else 0
        for h in READINGS:
            sd = subset_sides(repro.get(reading_name(h)) or {}, ws)
            if not sd or any(len(v) == 0 for v in sd.values()):
                continue
            lb = subset_sides({k[0]: v for k, v in lab[h].items()}, ws)
            rows = population_rows(STEP6A_POPULATIONS[DEVIATION], sd, lb, spreads, reading_name(h), nu,
                                   len(ws) * N_INPUT_FLIES - len(rp))
            poprows += [dict(r, setting='family', variant=DEVIATION, wirings_set=wset) for r in rows]
    pop = pd.DataFrame(poprows)
    readings, rrows = {}, []
    for wset, _ in wsets:
        if wset == 'partial' or wset in readings or pop.empty:
            continue
        for rd in ('all', 'held'):
            a = pop[(pop.wirings_set == wset) & (pop.reading == rd)]
            by = {v: a[a.variant == v].set_index('stage') for v in STEP6A_POPULATIONS}
            stages = [s for s in (*LOSS_STAGES, 'megamat_stage_error') if any(s in b.index for b in by.values())]
            for st in stages:
                x = step6a_reading(*(by[v].loc[st].to_dict() if st in by[v].index else None for v in STEP6A_POPULATIONS))
                rrows.append(dict(wirings_set=wset, reading=rd, stage=st, **x))
                if rd == 'all':
                    readings.setdefault(wset, {})[st] = x
    pops = {}
    if len(pop):
        for r in pop[pop.reading == 'all'].to_dict('records'):
            pops.setdefault(r['population'], {}).setdefault(r['wirings_set'], {})[r['stage']] = dict(
                stand_in=r['stand_in'], input_fly_term=r['input_fly_term'], model_flies=r['model_flies'],
                population_error=_num(r['population_error']), lab_common=_num(r['lab_common']),
                against_others_mean=_num(r['against_others_mean']), n_model_flies=int(r['n_model_flies']),
                n_unheld=int(r['n_unheld']), n_missing=int(r['n_missing']))
    classes = dict(description=STEP6A_DESCRIPTION, revision=STEP6A_REVISION['time'], gates=False,
                   settings={s: dict(base=SETTING_BASE[s], variant=SETTING_VARIANT[s]) for s in STEP6A_SETTINGS},
                   n_reads=n_reads, n_unheld=n_unheld, readings=readings, populations=pops)
    return dict(populations=pop, model_flies=pd.DataFrame(flyrows), readings_table=pd.DataFrame(rrows),
                classes=classes, n_reads=int(len(six)), n_unheld=n_unheld)


def step6a_tables(s6: Mapping) -> List[str]:
    """tables.md's step-6a section: the populations (all the flies) and the readings per wiring set and stage."""
    lines = ['## Step 6a: the deviations as measured (M) and permuted (P), beside the reproducible ones (R); classing '
             'nothing', '', STEP6A_DESCRIPTION + '.', '']
    pop = s6['populations']
    if len(pop):
        cols = ['wirings_set', 'population', 'stage', 'n_model_flies', 'n_unheld', 'n_missing', 'population_error',
                'lab_common', 'spread_as_scored', 'stand_in', 'mean_model_fly', 'term_input_fly', 'term_wiring',
                'input_fly_term', 'against_others_mean', 'model_flies']
        lines += ['### Populations (all the flies)', '', _md(pop[pop.reading == 'all'][cols]), '']
    rt = s6['readings_table']
    if len(rt):
        cols = ['wirings_set', 'reading', 'stage', 'repro_stand_in', 'measured_stand_in', 'permuted_stand_in',
                'structure_or_variance', 'repro_against_others_mean', 'measured_against_others_mean', 'e_applies',
                'e_at_upper_bound']
        lines += ['### The readings (the lock\'s revision; a description)', '',
                  "\"The structure's\": P a fair stand-in where R moves the stage nearer; \"the variance's\": P moves it "
                  "nearer too; \"(e) stands at the upper bound\": under M a model fly against the others stays under "
                  f"{FLIES_LOW}.", '', _md(rt[cols]), '']
    else:
        lines += ['No reading yet: neither variant has its 36 T reads (or all 72) saved.', '']
    return lines


def write_tables_score(pop: pd.DataFrame, groups: pd.DataFrame, a1: Mapping, a1_apl: Mapping, recomputed: Mapping,
                       classes: Mapping, step6a: Optional[Mapping] = None) -> None:
    """tables.md: the score's tables (populations, A1's fall, the APL16 somata, the recomputed classes, step 6a's
    populations and readings where its reads exist) above step 0's (tables_step0.md)."""
    lines = ["# E147: the megamat stage error read on model flies with each ORN fly's own input", '',
             f'Generated by {SCRIPT} at {H.code_commits()["queue_commit"][:7]} ({_now()}) from {classes["n_reads"]} '
             f'reads. Strength by odor in one-fly units, lower is better.', '']
    if len(pop):
        cols = ['population', 'reading', 'stage', 'n_model_flies', 'n_unheld', 'population_error', 'lab_common',
                'spread_as_scored', 'stand_in', 'mean_model_fly', 'term_input_fly', 'term_wiring', 'term_remainder',
                'input_fly_term', 'against_others_mean', 'model_flies', 'spread_model_flies']
        lines += ['## Megamat populations', '', _md(pop[cols]), '']
    if a1:
        lines += ["## A1's fall, model fly by model fly (A1 minus the family on the same wiring and input fly)", '']
        for rd, v in a1.items():
            st = v['stage']
            if st.get('n'):
                lines.append(f"- {rd} flies: fall {st['fall']:.4f} over {st['n']} model flies, se {st['se']:.4f} (over "
                             f"input flies {st['se_over_input_flies']:.4f}, over wirings {st['se_over_wirings']:.4f}); "
                             f"stands: {st['stands']}.")
        lines.append('')
    if len(groups):
        cols = ['group', 'condition', 'n_model_flies', 'n_unheld', 'strength_population', 'strength_lab_common',
                'spread_as_scored', 'stand_in', 'strength_term_input_fly', 'strength_term_wiring', 'input_fly_term',
                'strength_against_others_mean', 'model_flies', 'strength_spread_model_flies']
        lines += ['## The APL16 somata', '', _md(groups[cols]), '']
    if a1_apl:
        lines += ['## A1 at the APL16 somata, by input fly', '']
        for c, v in a1_apl.items():
            if v.get('n_flies'):
                lines.append(f"- APL {c}: mean change {v['mean_change']:+.4f}; worse than the spread ({v['spread']:.3f}) "
                             f"on {v['share_worse']:.0%} of the input flies, better on {v['share_better']:.0%}.")
        lines.append('')
    if recomputed:
        lines += ['## The classes recomputed with the spreads over model flies (changed only; a description)', '']
        for k, v in recomputed.items():
            lines.append(f"- {k}: " + ('; '.join(f"{x['table']} {x['arm']}: {x['cls_scored']} -> {x['cls_new']}" for x in v)
                                      if v else 'no class changes'))
        lines.append('')
    if step6a is not None:
        lines += step6a_tables(step6a)
    lines += ['## Outcomes (a description)', '', '```', json.dumps(classes.get('outcomes'), indent=1, default=_json_default),
              '```', '']
    st0 = results_dir() / 'tables_step0.md'
    if st0.exists():
        lines += ['---', '', st0.read_text()]
    (results_dir() / 'tables.md').write_text('\n'.join(lines))


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'lock6a', 'step0', 'step1', 'step2', 'step6a', 'score'))
    a = ap.parse_args(argv)
    if a.part in ('step1', 'step2', 'step6a'):
        raise SystemExit(f"{a.part} is the runner's: ./venv/bin/python -u {RUNNER} --part reads --step "
                         f"{a.part[len('step'):]}")
    {'lock': part_lock, 'lock6a': part_lock_step6a, 'step0': part_step0, 'score': part_score}[a.part]()


if __name__ == '__main__':
    main()
