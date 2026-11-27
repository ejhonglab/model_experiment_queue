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

# the classes (the entry's, fixed beforehand)
TERM_SMALL, TERM_PART = 0.5, 1.0          # input-fly term / wiring term: "small" below, "a part" below, else "as large"
FLIES_LOW, FLIES_HIGH = 0.75, 1.25        # a model fly against the others: "too little" below, "as the flies", "too much"
IMPROVE_SE = 2.0                          # A1's fall: beyond this many standard errors
CONTROL_FALL = 0.040                      # and beyond the control fit's fall (E136 classes.json control_fall, 0.0403)
UNHELD_MAX_SHARE = 1.0 / 9.0              # no class for a population with more than this share of its model flies unheld
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
    "-mean > IMPROVE_SE x se and -mean > CONTROL_FALL",
    "unheld: a read whose row has held False (its counts not at both targets, or a factor at a bound); unheld model "
    "flies are left out of every value and counted; a population with more than a ninth unheld has no class",
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
            unheld_max_share=UNHELD_MAX_SHARE),
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
