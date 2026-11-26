#!/usr/bin/env python3
"""E140, aa at the step from the claws to the somata: why the model's aa falls where the flies' rises, read in the
recordings, in every saved megamat run and in the lab's own measurements of aa's input (step 0, no simulation).

The entry is item E140 of experiment_priorities.md (proposed by viki 2026-11-25 19:00 from E138's first next step;
planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-25 18:39; moved in under the idle-host rule).
"The family" is C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings (E83's; S real, P1, P2, M1 and T P3,
P4, M2, M3), each read as one fly; the views, one-fly units, readings and flies are E134's (kcq/population_fit.py,
kcq/stage_readings.py, data_prep/stage_population_data.py and its lock experiments/E134_lock.json). Strength by odor
has two parts: each odor's response rate (part 'rate') and its mean response (part 'mean'), each as a multiple of the
mean over the panel's 17 odors. Nothing is run, fitted or searched: saved runs and archives are read again.

The guard: no statistic of aa (or of any odor) is computed before this experiment's lock exists, is committed with no
uncommitted change and equals its copy in the shared results (`require_lock`); E134's lock guards every reading of the
flies (`stage_readings.require_lock`). Every part but 'lock' calls both first. The lock (`lock_contents`) fixes:
- the readings of the flies' aa at the somata (`SOMA_READINGS`): the lock's (E134's), each of the three presentations
  alone, each half of the response window, and E134's alternative soma criteria, with which of them count for class
  (iv) and which go beside;
- "concentration" of an odor's responders at a stage (`concentration`): the share of the summed response carried by
  the top tenth of the responding units, with the responders' coefficient of variation beside;
- the recordings: the PN dendrites (7 flies), PN boutons (5; the bouton fly 2022-10-21 fly02 left out), the five claw
  recordings (Remy's four claw flies and the same-fly recording), the KC somata (Remy's four 2022 flies, and the same
  fly); the ORNs beside. Remy's older megamat flies (pyolfsysm's data/from_remy/2024-11-12) are never read;
- ranks among the 17 megamat odors (`ranks`): rank 1 the lowest; "the lowest third" ranks 1 to 6, "the highest third"
  ranks 12 to 17;
- the four classes' exact thresholds (`CLASSES`), the halves of the soma flies (E136's seeded split), and the saved
  model runs read (`MODEL_GROUPS`).

Parts (CLI `--part`):
- lock: writes results/E140_aa_claw_soma_step/lock.json and its committed copy experiments/E140_lock.json (commit it
  before any other part); refuses to overwrite a committed lock. Reads no statistic: file hashes, fly names and the
  seeded split only.
- reduce (no statistic of aa): the four megamat soma flies' frames (Remy's suite2p trials files on matrix, read only;
  sha256 checked against E134's files.csv) reduced to the two halves of the response window, each cell's per-trial
  window mean minus its pre-odor mean over its pre-odor frame SD, as data_prep/stage_population_data.py reduces the
  whole window (the whole window recomputed must equal its frames_reduced.parquet); soma_windows.parquet.
- flies (design item 1): aa's rate part, mean part and concentration at the dendrites, boutons, claws and somata of
  every recording, with ranks; the same-fly step; aa's somatic entry under every reading of the lock; who responds
  with aa (aa's row of the sharing view at every stage, and the number of odors aa's somatic responders respond to).
- models (design item 2): the same statistics for the saved megamat runs (`MODEL_GROUPS`): the family held on its
  eight wirings with APL on and off (and as fitted, beside), E136's five arms, E135's imposed dendrites and E139's
  settings as they arrive (`e139_runs`); beside them aa's somatic responders by the glomeruli of their claws (E138's
  IR-fed class against the rest) and by KC type, and each odor's entry at each stage under A1.
- follows (design item 3): over the 17 odors, the somata's strength against the input's sum, its largest entry and its
  breadth, and against the upstream stages, in the model and in the flies; aa's residual with aa left out of the fit.
- input (design item 4): aa's model input beside each of the seven ORN flies that were given aa (by E85's map of the
  lab's dF/F to spike deltas), beside Hallem's aa on the shared glomeruli, and beside the dendrite recordings on their
  28 glomeruli.
- classes: the classes (i) to (iv) as the lock fixes them, and the outcome they point to; classes.json, tables.md.
- figures: a few figures of the above.
- all: flies, models, follows, input, classes, figures.

Step 1 (28 runs, aa's input from each ORN fly's own vector on the four T wirings) runs only on class (i), after the
worker reads step 0; step 2 waits for Tom's or Yang's word. Neither is in this script yet.

Concentrations: megamat at 10^-3 in every recording and input used; Hallem's aa is at 10^-2 (within 1 log10 unit),
read beside only. No APL16 number is computed. Validation2 is not read.

Holds (CLAUDE.md): only Pratyush's hemibrain (the family's wirings) and, through E138's lock, nothing else; no olfsysm
change; no fit, no search, no simulation; only Yang's processed products (none are read here); Remy's processed data
without the bouton fly 2022-10-21 fly02 and without her older megamat flies. Nothing is published.

Outputs (<shared>/results/E140_aa_claw_soma_step/): lock.json; soma_windows.parquet, soma_windows.json;
flies_entries.csv, flies_sharing.csv, flies_sharing_summary.csv, flies_breadth.csv, flies_same_fly_step.csv,
flies_soma_readings.csv; models_entries.csv, models_sharing.csv, models_sharing_summary.csv, models_breadth.csv,
models_anatomy.csv, models_a1_stages.csv; follows.csv, follows_points.csv, input_features.csv; input_aa.csv,
input_summary.csv; classes.json, tables.md; fig_*.png; provenance_<part>.json.

Usage (from the repository root; one thread under nice):
    ./venv/bin/python -u experiments/E140_aa_claw_soma_step.py --part lock      # then commit experiments/E140_lock.json
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E140_aa_claw_soma_step.py \\
        --part reduce   # then --part all
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import socket
import subprocess
import time
from collections import namedtuple
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E140_aa_claw_soma_step'
SCRIPT = 'experiments/E140_aa_claw_soma_step.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E140_lock.json'
E134_LOCK_COPY = H.HERE / 'experiments' / 'E134_lock.json'
E138_LOCK_COPY = H.HERE / 'experiments' / 'E138_classes_lock.json'

ODOR = 'aa'
PANEL = 'megamat'
N_ODORS = 17
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
STRENGTH_PARTS: Tuple[str, ...] = ('strength/rate', 'strength/mean')
# the stages read (the ORNs beside: the input as the terminals are recorded); the entry's "boutons, claws and somata"
STAGES: Tuple[str, ...] = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
MAIN_STAGES: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')

# "concentration" of an odor's responders in one recording: the share of the responders' summed response carried by the
# top tenth of them (the Lorenz curve of the responders' values, largest first, at a tenth of their number, linear
# between units), with the responders' coefficient of variation (SD over mean, ddof 0) beside; undefined with fewer
# than MIN_RESPONDERS responders
TOP_SHARE = 0.1
MIN_RESPONDERS = 5
CV_DDOF = 0
# ranks among the 17 megamat odors: ascending (rank 1 the lowest), ties averaged; an odor undefined in a recording has
# no rank there and the others are ranked among the odors defined; "the lowest third" the THIRD lowest ranks (1 to 6),
# "the highest third" the THIRD highest (12 to 17 of 17)
THIRD = 6

# the soma flies' readings (design item 1; the lock's "readings of the flies' aa at the somata"): name -> what it is;
# 'counts' says whether the reading counts for class (iv) (the others go beside)
SOMA_RATE = 0.0915           # the megamat operating point (E134's lock: Remy's four flies by her criterion)
LOCK_SOMA_CRITERION = 5.9486  # E134's lock: Remy's megamat somata at the level that gives 0.0915
APL16_SOMA_CRITERION = 1.9684  # E134's lock's first soma criterion (the APL16 somata's)
FULL_WINDOW = (2.0, 8.0)       # E134's lock: Remy's megamat somata, response window (s after onset, closed)
PRE_WINDOW = (-10.0, -0.5)     # and pre-odor window (closed)
# the two halves of the response window (closed; Remy's frames are at 0.25 s: 2.00 to 4.75 s, 12 frames, and 5.00 to
# 8.00 s, 13 frames, every frame of the whole window in exactly one half)
HALF_WINDOWS: Dict[str, Tuple[float, float]] = {'early': (2.0, 4.8), 'late': (5.0, 8.0)}
SomaReading = namedtuple('SomaReading', 'name source criterion counts what')
SOMA_READINGS: Tuple[SomaReading, ...] = (
    SomaReading('lock', 'population', 'lock', True,
                "E134's lock: the population file's trial means (2 to 8 s minus -10 to -0.5 s over the pre-odor frame "
                'SD), criterion 5.9486'),
    *(SomaReading(f'presentation_{r + 1}', f'presentation:{r}', 'rate', True,
                  f'the {("first", "second", "third")[r]} presentation of each odor alone (the population file\'s '
                  'single trials), the criterion re-derived by the lock\'s rule (the level at which the mean over the '
                  'four flies of each fly\'s response rate over the 17 odors is 0.0915)') for r in range(3)),
    *(SomaReading(f'window_{w}', f'window:{w}', 'rate', True,
                  f'the {w} half of the response window ({HALF_WINDOWS[w][0]:g} to {HALF_WINDOWS[w][1]:g} s, closed), '
                  "trial means, each cell's pre-odor mean and SD as the lock's; the criterion re-derived by the lock's "
                  'rule') for w in HALF_WINDOWS),
    SomaReading('criterion_apl16', 'population', 'apl16', True,
                "E134's alternative soma criterion: the lock's first, the APL16 somata's 1.9684 (E134's 'at 1.9684' "
                'column)'),
    SomaReading('criterion_remy', 'remy_binary', 'remy', True,
                "E134's alternative soma criterion: Remy's own responder calls (her df_stim_responders_from_peak.csv "
                "per fly, the lock's remy_megamat_criterion_rate; the rate part only, the mean part the lock's)"),
    *(SomaReading(f'presentation_{r + 1}_at_lock', f'presentation:{r}', 'lock', False,
                  f'beside: presentation {r + 1} alone at the lock\'s criterion 5.9486') for r in range(3)),
    *(SomaReading(f'window_{w}_at_lock', f'window:{w}', 'lock', False,
                  f'beside: the {w} half window at the lock\'s criterion 5.9486') for w in HALF_WINDOWS),
    SomaReading('per_fly_level', 'population', 'per_fly', False,
                "beside: each fly at its own level giving 0.0915 (Remy's two levels, E134's next step)"),
)

# the classes, fixed before any statistic of aa (the entry's, made exact)
CLASSES = {
    'i': dict(name='few and strong in the flies, many and weak in the model',
              rule='in at least 4 of the 5 claw recordings aa is in the lowest third on the rate part and in the highest '
                   "third on concentration, and on 6 of 8 wirings the family's aa (held at 0.0915, APL on, E134's "
                   'default claw reading) is in the highest third on the rate part and the lowest third on concentration',
              flies_min=4, flies_n=5, model_min=6, model_n=8, model_group='family_held', condition='on',
              stage='kc_claw', reading='curve_noise'),
    'ii': dict(name='apart',
               rule="at the flies' somata (E134's lock reading, the four flies) aa's mean cosine distance to the 16 other "
                    "odors on who responds (the sharing view's part 'who', the half's mean over its flies) is among "
                    'the three largest of the 17 odors on both halves of the flies (E136\'s seeded split)',
               top=3, part='sharing/who', stage='kc_soma', reading='lock'),
    'iii': dict(name="inside the ORN flies' spread",
                rule="at least one ORN fly's aa vector brings the somata's aa entry (the rate part, the family held, APL on, "
                     "E134's default soma reading) into the soma flies' range (the min to max over the four flies of aa's "
                     "rate entry, the lock's reading) on 3 of the 4 T wirings",
                min_wirings=3, wirings=list(T_WIRINGS), read_in='step 1 (not run in step 0)'),
    'iv': dict(name="the reading's",
               rule="under a reading of the lock that counts (SOMA_READINGS 'counts'), the rank of the flies' mean rate "
                    "entry of aa among the 17 odors is 11 or lower (aa leaves the highest third)",
               part='strength/rate'),
}
OUTCOMES = {
    'a': '(i) holds and (iii) does not: the lab\'s aa input, as the model reads it, lacks a strong channel that few claws '
         'carry; Tom\'s to weigh (step 2 on his or Yang\'s word)',
    'b': "(iii) holds: the answer is inside what seven ORN flies showed; the choice of input is Tom's",
    'c': "(i) does not hold: the flies' claws see aa as the model's do, and the step itself differs; item 2's tables say "
         'which KCs the model loses',
    'd': "(iv) holds: aa's strength at the somata is the reading's",
}

# the saved megamat runs (design item 2): group -> (source, tags by wiring); E139's settings are a hook (`e139_runs`)
ModelGroup = namedtuple('ModelGroup', 'name source wirings what')
MODEL_GROUPS: Tuple[ModelGroup, ...] = (
    ModelGroup('family_held', 'E134:step3:megamat_held_{w}', WIRINGS,
               "the family held at the megamat flies' 0.0915 (E134's step 3), APL on and off"),
    ModelGroup('family_fitted', 'E134:step3:megamat_fitted_{w}', WIRINGS, "beside: the family as fitted (E134's step 3)"),
    *(ModelGroup(a, f'E136:read:{a}_megamat_{{w}}', T_WIRINGS,
                 f"E136's arm {a}, its megamat reads on the T wirings (the S wirings are its fit's evaluations)")
      for a in ('A0', 'A1', 'A2', 'A3', 'A4')),
    ModelGroup('e135_imposed', 'E135:step2:megamat_imposed_{w}', WIRINGS,
               "E135's recorded dendrites imposed (step 2, the table in noise units)"),
)
E139_EXPERIMENT = 'E139_normalization_ladder'
# E136's measure whose reads are read (the adopted measure's, measure_strength/ in its results and archive)
E136_MEASURE = 'strength'
# the KC-soma readings read beside the default on the family held (E134's alternatives)
MODEL_SOMA_BESIDE: Tuple[str, ...] = ('raw', 'affine_noise')
MODEL_CLAW_BESIDE: Tuple[str, ...] = ('curve_noise_b14', 'curve_noise_share25')

# the item-3 regressions: y = a + b x by least squares over the odors other than aa
FOLLOW_INPUT_FEATURES: Tuple[str, ...] = ('input_sum', 'input_max', 'input_breadth')
FOLLOW_STAGES: Tuple[str, ...] = ('kc_claw', 'pn_bouton', 'pn_dendrite', 'orn')
# item 4: Hallem's aa (drosolf's Hallem and Carlson 2006 table, deltas) and the threshold of a "responding" receptor
HALLEM_ODOR = 'acetic acid'
HALLEM_RESPONDING_HZ = 20.0


def out_dir() -> Path:
    """<shared>/results/E140_aa_claw_soma_step/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir() -> Path:
    """outputs/E140_aa_claw_soma_step/ (local, created)."""
    return H.outputs_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def file_entry(path: Path, sha: Optional[str] = None) -> dict:
    """{'path', 'bytes', 'mtime', 'sha256'} of a file (the sha256 computed unless given)."""
    p = Path(path)
    st = p.stat()
    return dict(path=str(p), bytes=int(st.st_size), mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime)),
                sha256=sha if sha is not None else _sha256(p))


def _json_default(x):
    if hasattr(x, 'item'):
        return x.item()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')
    os.replace(tmp, path)


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    p = out_dir() / name
    tmp = p.with_name(f'.{p.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, p)
    return p


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git with no uncommitted change."""
    try:
        if subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)],
                          capture_output=True).returncode != 0:
            return False
        return subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)],
                              capture_output=True).returncode == 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# the lock

def recordings() -> Dict[str, dict]:
    """The recordings read, by stage: their population (viki's (stage, panel, cohort) keys) and fly names, from the
    population files (names only); the claws' five recordings (Remy's four claw flies and the same-fly recording); the
    exclusions. No value is read."""
    from data_prep import stage_population_data as D
    out = {}
    for stage, cohorts in (('orn', ('wild_type',)), ('pn_dendrite', ('wild_type',)), ('pn_bouton', ('wild_type',)),
                           ('kc_claw', ('wild_type', 'same_fly')), ('kc_soma', ('wild_type', 'same_fly'))):
        for c in cohorts:
            p = D.load_population(stage, PANEL, c)
            out[f'{stage}__{c}'] = dict(stage=stage, panel=PANEL, cohort=c, flies=list(p.flies),
                                        criterion=float(p.criterion), units=p.units, odors=list(p.odors),
                                        role=('beside' if stage == 'orn' else 'read'))
    claws = [dict(cohort='wild_type', fly=f) for f in out['kc_claw__wild_type']['flies']]
    claws += [dict(cohort='same_fly', fly=f) for f in out['kc_claw__same_fly']['flies']]
    if len(claws) != 5:
        raise RuntimeError(f'{len(claws)} claw recordings, not 5')
    out['claw_recordings'] = dict(recordings=claws, note='the same-fly recording (2023-04-11 fly 3, megamat1_redo) is '
                                  "another acquisition of the fly of the claw recording 2023-04-11__fly03__megamat1: "
                                  'five recordings of four flies')
    out['excluded'] = dict(
        flies=[dict(stage='pn_bouton', fly=x['fly'], reason=x['reason']) for x in D.EXCLUDED_FLIES],
        older_megamat="Remy's older megamat flies (pyolfsysm's data/from_remy/2024-11-12, 2018-2019, partial panels): "
                      'never read (Tom, 2026-11-25 10:19, relayed by hal: "leave that data out still")')
    if any(f.startswith('2022-10-21__fly02') for f in out['pn_bouton__wild_type']['flies']):
        raise RuntimeError('the excluded bouton fly is in the bouton population')
    return out


def soma_halves(flies: Sequence[str]) -> Dict[str, List[str]]:
    """E136's seeded split of the soma flies (`population_fit.split_halves(n, 1, 0)[0]`, as E138's `half_names`): 'fit'
    the first part, 'held' the rest, in the population's order."""
    from kcq import population_fit as P
    a, b = P.split_halves(len(flies), 1, 0)[0]
    return dict(fit=[flies[i] for i in a], held=[flies[i] for i in b])


def soma_frame_files() -> Dict[str, dict]:
    """Per megamat soma fly, Remy's suite2p trials file on matrix (its bytes and sha256 from E134's sources/files.csv;
    checked when read, `part_reduce`) and the copy of her respvec file (sha256 computed)."""
    from data_prep import stage_population_data as D
    fl = pd.read_csv(D.SOURCES_DIR / 'files.csv')
    out = {}
    for acq in D.REMY_MEGAMAT_ACQS:
        src = D.megamat_src(acq)
        row = fl[fl.source == str(src['trials'])]
        if len(row) != 1:
            raise RuntimeError(f"{acq}: the trials file is not in E134's files.csv once")
        date, fly, _ = acq.split('/')
        rv = D.megamat_copy_dir() / D.megamat_acq_label(acq) / src['respvec'].name
        out[f'{date}/{int(fly)}'] = dict(acq=acq, trials=dict(path=str(src['trials']),
                                                              bytes=int(row.source_bytes.iloc[0]),
                                                              sha256=str(row.source_sha256.iloc[0]),
                                                              source='E134 sources/files.csv'),
                                         respvec=file_entry(rv))
    return out


def remy_binary_files() -> Dict[str, dict]:
    """Remy's per-fly responder calls (pyolfsysm's data/from_remy/response_rates/.../by_acq_csvs/<acq>/
    df_stim_responders_from_peak.csv) and the acquisitions' ledger df_acqs.csv, with sha256."""
    from pyolfsysm import util as U
    out = {p.parent.name: file_entry(p) for p in sorted(Path(U.remy_binary_response_dir).glob(
        f'*/{U.remy_fly_binary_response_fname}'))}
    out['df_acqs.csv'] = file_entry(Path(U.remy_sparsity_dir) / 'df_acqs.csv')
    return out


def hallem_file() -> dict:
    """drosolf's Hallem and Carlson 2006 table, with sha256."""
    from drosolf import orns
    return file_entry(Path(orns.hallem_csv))


def lock_contents() -> dict:
    """The lock: the inputs (files with sha256, as E137's lock records them), the recordings, the soma flies' halves,
    the definitions (concentration, ranks, thirds, the readings of the flies' somata), the classes and outcomes, the
    model groups and their readings, and item 3's and item 4's rules. No statistic is read."""
    from data_prep import stage_population_data as D
    from kcq import stage_readings as R
    from pyolfsysm import util as U
    e134 = json.loads(E134_LOCK_COPY.read_text())
    recs = recordings()
    used = [k for k in recs if '__' in k]
    pops = {k: dict(file=e134['populations'][f'{recs[k]["stage"]}__{PANEL}__{recs[k]["cohort"]}']['file'],
                    criterion=e134['populations'][f'{recs[k]["stage"]}__{PANEL}__{recs[k]["cohort"]}']['criterion'])
            for k in used}
    soma = recs['kc_soma__wild_type']['flies']
    e138 = json.loads(E138_LOCK_COPY.read_text())
    inputs = dict(
        e134_lock=file_entry(E134_LOCK_COPY), e138_classes_lock=file_entry(E138_LOCK_COPY),
        populations=pops,
        frames_reduced=file_entry(D.megamat_copy_dir() / 'frames_reduced.parquet'),
        soma_frames=soma_frame_files(),
        remy_binary=remy_binary_files(),
        olfdb_dataset_1=file_entry(D.SOURCES_DIR / 'olfdb' / 'ds1_long.csv'),
        hallem=hallem_file(),
        dff_to_spiking=file_entry(Path(U.data_root) / 'internal' / 'dff2spiking_model_choices.csv'),
        megamat_input=file_entry(Path(U.data_root) / 'internal' / 'megamat_mean_est_spike_deltas.csv'))
    definitions = dict(
        strength="kcq/population_fit.py's strength view on the 17 megamat odors of each recording (population_views): "
                 "part 'rate', each odor's response rate (the fraction of units above the stage's criterion) over the "
                 "mean over the 17 odors; part 'mean', each odor's mean response over all units likewise; the flies' "
                 'mean entry the mean over recordings, their range the min and max',
        concentration=dict(rule="the share of the responders' summed response carried by the top tenth of them: the "
                                "responders' values (units above the stage's criterion) sorted largest first, the "
                                'Lorenz curve at a tenth of their number, linear between units (a unit at the boundary '
                                'counts by its fraction)', top_share=TOP_SHARE, min_responders=MIN_RESPONDERS,
                           beside="the responders' coefficient of variation (SD over mean, ddof 0)",
                           model="on the model's read population (values after the stage's curve and noise, its "
                                 "responders above the flies' criterion), as the flies are read",
                           undefined=f'fewer than {MIN_RESPONDERS} responders: NaN, no rank'),
        ranks=dict(rule='ascending among the 17 megamat odors (rank 1 the lowest), ties averaged; an undefined odor has '
                        'no rank and the others are ranked among the odors defined',
                   lowest_third=f'the {THIRD} lowest ranks (1 to 6 of 17)',
                   highest_third=f'the {THIRD} highest ranks (12 to 17 of 17)', third=THIRD),
        sharing_row="aa's 16 cosine distances of the sharing view (parts 'who' and 'how_strongly'), each recording's and "
                    "the flies' mean; aa's mean over its 16 pairs, ranked among the 17 odors' means (rank 17 the most "
                    'apart)',
        responders_breadth="per soma recording, the number of the 17 odors each of an odor's responders responds to "
                           "(E134's breadth bins for 17 odors beside the mean)",
        same_fly_step="the same-fly recording's somata minus its claws, each odor's rate and mean entries (one fly, a "
                      'description); beside, the four soma flies\' mean minus the four claw flies\' mean',
        soma_readings=[r._asdict() for r in SOMA_READINGS],
        soma_rate=SOMA_RATE, lock_soma_criterion=LOCK_SOMA_CRITERION, apl16_soma_criterion=APL16_SOMA_CRITERION,
        criterion_search=dict(D.CRITERION_SEARCH, rule="data_prep/stage_population_data.py's criterion_at_rate (the "
                              "mean over flies of each fly's response rate over the 17 odors at 0.0915), rounded"),
        windows=dict(full=FULL_WINDOW, pre=PRE_WINDOW, halves=HALF_WINDOWS,
                     rule="each cell's per-trial mean over the half's frames minus its mean over the pre-odor frames, "
                          "over its pre-odor frame SD (the RMS over trials of each trial's SD), as "
                          "stage_population_data.reduce_megamat_trials reduces the whole window; Remy's good cells"),
        soma_halves=soma_halves(soma))
    model = dict(
        groups=[g._asdict() for g in MODEL_GROUPS],
        e139=f'a hook: {E139_EXPERIMENT} has no saved megamat runs yet; `e139_runs` reads any that arrive in '
             f'outputs/{E139_EXPERIMENT}/ or the model-output archive (stage_readings\' step-3 layout, megamat input)',
        readings=R.reading_definitions(),
        reading_rule="E134's default reading of each stage, against all the stage's recorded flies (stage_readings' "
                     'FlyCache); APL off read with the APL-on flies\' curve, noise and criterion (no megamat flies with '
                     "APL output blocked); each wiring's seed index its position in WIRINGS (T_WIRINGS for E136's arms)",
        soma_beside=list(MODEL_SOMA_BESIDE), claw_beside=list(MODEL_CLAW_BESIDE),
        rows="each wiring; 'common' the wirings' mean views (rate and mean parts) and the mean of the wirings' "
             "concentrations; 'mean' the mean over wirings of each wiring's entries",
        anatomy=dict(responders="the KCs with two or more spikes (the model's own criterion, at which the family is "
                                'held), each odor', ir_glomeruli=e138['classes']['ir']['glomeruli'],
                     ir_rule="E138's lock class 'ir'", claws="each responder's share of claws on an IR-fed glomerulus "
                                                             "(the run's units: claw_glomerulus, kc_of_claw)",
                     kc_types="the spike counts' kc_type level"))
    follows = dict(
        fit='y = a + b x by least squares over the 16 odors other than aa; aa\'s residual y - (a + b x), in the entry\'s '
            'units and over the residual SD of the 16 (ddof 2); Pearson r and Spearman rho over the 16 beside; every '
            "odor's leave-one-out residual and aa's rank by its absolute value",
        y=["the KC somata's rate entry (the flies' mean; each model group's 'common' row)",
           'beside: the mean entry'],
        x=dict(input_sum="the sum over the model's 36 glomeruli of the odor's megamat input (spike deltas, signed)",
               input_max='its largest entry', input_breadth="the participation ratio of its positive part, "
                                                            "(sum x+)^2 / sum x+^2",
               stages="the same side's rate entry at the claws, boutons, dendrites and ORNs (the flies' at theirs, the "
                      "model's at its own)"))
    input_ = dict(
        model="the megamat input's aa column on the 36 glomeruli (E134's panel_input('megamat'))",
        orn_flies="olfdb dataset 1 (E134's copy, sources/olfdb/ds1_long.csv), the flies given aa @ -3: each fly scaled "
                  "as E85 scales it ('to-avg-max': its largest |trial mean| over all its panels and ROIs brought to "
                  "the lab's avg_flymax), then the least-squares slope through the origin of the megamat input on the "
                  "flies' mean (36 glomeruli x 17 odors, each entry the mean over the flies that have it); E85's slope "
                  '(12 odors) beside',
        hallem=f"drosolf's Hallem and Carlson 2006 table (deltas, add_sfr False), odor {HALLEM_ODOR!r} (10^-2), on the "
               'glomeruli it shares with the model\'s 36 by name; a receptor "responds" at >= '
               f'{HALLEM_RESPONDING_HZ:g} spikes per s',
        dendrites="the PN-dendrite recordings' aa column (7 flies, trial-noise SD units) on E138's lock's 28 "
                  "glomeruli, beside the family held's PN stage (the PN rates before the curve, wiring 'real') and "
                  'the model input on them',
        dendrite28=e138['classes']['not_dendrite28']['dendrite28'])
    return dict(experiment=EXPERIMENT, script=SCRIPT, part='lock', written=_now(), host=_host(),
                commits=H.code_commits(), odor=ODOR, panel=PANEL, odors=list(D.MEGAMAT_ODORS), inputs=inputs,
                recordings=recs, definitions=definitions, classes=CLASSES, outcomes=OUTCOMES, model=model,
                follows=follows, input=input_,
                step1=dict(rule="only on class (i): aa's input taken from each of the seven ORN flies' own vectors by "
                                'the lab\'s map, the other odors unchanged, the family held, on the four T wirings (28 '
                                'runs); nothing searched', status='not in this script'),
                step2=dict(rule="held for Tom's or Yang's word: aa's input on E138's IR-fed class at 1.5, 2 and 3 times, "
                                'with seven other glomeruli of matched summed input as the control (24 runs)',
                           status='held'),
                holds=['no simulation, no model run, no olfsysm change, no input factor searched in step 0',
                       'validation2 not read', "only Pratyush's hemibrain (the family's wirings)",
                       "only Yang's processed products (none read here)",
                       "Remy's processed data without the bouton fly 2022-10-21 fly02 and her older megamat flies"],
                note='written before any statistic of aa (or of any odor) was computed; every part but lock refuses '
                     'until this file is committed (require_lock)')


def part_lock() -> dict:
    """results/E140_aa_claw_soma_step/lock.json and its committed copy experiments/E140_lock.json (commit it before any
    other part); refuses to overwrite a committed lock."""
    if LOCK_COPY.exists() and _committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is committed already; the lock is fixed')
    lock = lock_contents()
    text = json.dumps(lock, indent=1, default=_json_default) + '\n'
    (out_dir() / LOCK_NAME).write_text(text)
    LOCK_COPY.write_text(text)
    print(f'wrote {out_dir() / LOCK_NAME} and {LOCK_COPY} (commit the copy before any other part)', flush=True)
    return lock


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock: experiments/E140_lock.json exists, is tracked with no uncommitted change, and has the same
    bytes as results/.../lock.json; with E134's lock (`stage_readings.require_lock`) required too. PermissionError
    otherwise (no statistic is read before)."""
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it); no statistic of aa is read '
                              'before it is committed')
    if not _committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = (out_dir() / LOCK_NAME) if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = _sha256(copy)
    return lock


def require_all() -> Tuple[dict, dict]:
    """(this experiment's lock, E134's lock), both required before any statistic."""
    from kcq import stage_readings as R
    return require_lock(), R.require_lock()


# ---------------------------------------------------------------------------
# the statistics (pure functions)

RANKED: Tuple[str, ...] = ('rate', 'mean', 'concentration', 'cv')
ENTRY_COLUMNS: Tuple[str, ...] = ('rate', 'mean', 'response_rate', 'mean_response', 'concentration', 'cv',
                                  'n_responders', 'n_units')


def concentration(values: np.ndarray, top: float = TOP_SHARE, min_n: int = MIN_RESPONDERS) -> float:
    """The share of the summed values carried by the largest `top` fraction of them (the Lorenz curve, largest first,
    at top x n, linear between units: a unit at the boundary counts by its fraction); NaN with fewer than `min_n`
    values or a non-positive sum."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    n = len(v)
    if n < max(int(min_n), 1):
        return np.nan
    s = float(v.sum())
    if not s > 0:
        return np.nan
    v = np.sort(v)[::-1]
    k = float(top) * n
    m = int(np.floor(k))
    part = float(v[:m].sum()) + ((k - m) * float(v[m]) if m < n else 0.0)
    return part / s


def responders_cv(values: np.ndarray, min_n: int = MIN_RESPONDERS) -> float:
    """SD (ddof CV_DDOF) over mean of the responders' values; NaN with fewer than `min_n` (or two) values or a
    non-positive mean."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if len(v) < max(int(min_n), 2):
        return np.nan
    mu = float(v.mean())
    return float(v.std(ddof=CV_DDOF) / mu) if mu > 0 else np.nan


def ranks(x) -> pd.Series:
    """Ascending ranks (1 the lowest, ties averaged) among the defined entries; NaN where undefined."""
    return pd.Series(x, dtype=float).rank(method='average', ascending=True)


def third(rank: float, n_defined: int, k: int = THIRD) -> str:
    """'lowest' (rank <= k), 'highest' (rank > n_defined - k), else 'middle'; '' for an undefined rank."""
    if rank is None or not np.isfinite(rank):
        return ''
    if rank <= k:
        return 'lowest'
    if rank > n_defined - k:
        return 'highest'
    return 'middle'


def add_ranks(df: pd.DataFrame, cols: Sequence[str] = RANKED) -> pd.DataFrame:
    """`df` (one row per odor) with '<col>_rank', '<col>_third' and '<col>_n_ranked' for each of `cols`."""
    df = df.copy()
    for c in cols:
        r = ranks(df[c].to_numpy(dtype=float))
        n = int(r.notna().sum())
        df[f'{c}_rank'] = r.to_numpy()
        df[f'{c}_third'] = [third(v, n) for v in r]
        df[f'{c}_n_ranked'] = n
    return df


def odor_stats(pop, odors: Sequence[str]) -> pd.DataFrame:
    """Per odor of one population (hal's `Population`): the strength parts' entries ('rate', 'mean'; population_views
    on `odors`), the response rate and mean response, the responders' concentration and CV (graded populations only),
    the number of responders and of units, the mean entry over the rate entry, and each of RANKED's rank and third."""
    from kcq import population_fit as P
    od = list(odors)
    views = P.population_views(pop, od, parts=list(STRENGTH_PARTS))
    V = pop.values.reindex(columns=od).astype(float)
    Rr = pop.responds.reindex(columns=od).astype(float)
    rows = []
    for o in od:
        v, r = V[o].to_numpy(), Rr[o].to_numpy()
        ok = np.isfinite(v)
        resp = v[ok & (np.nan_to_num(r) > 0.5)]
        rows.append(dict(odor=o, rate=float(views['strength/rate'].get(o, np.nan)),
                         mean=float(views['strength/mean'].get(o, np.nan)) if 'strength/mean' in views else np.nan,
                         response_rate=float(np.mean(r[ok])) if ok.any() else np.nan,
                         mean_response=float(np.mean(v[ok])) if ok.any() else np.nan,
                         concentration=concentration(resp) if pop.graded else np.nan,
                         cv=responders_cv(resp) if pop.graded else np.nan,
                         n_responders=int(len(resp)), n_units=int(ok.sum())))
    df = pd.DataFrame(rows)
    with np.errstate(invalid='ignore', divide='ignore'):
        df['mean_over_rate'] = df['mean'] / df['rate']
    return add_ranks(df)


def mean_rows(per: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """The recordings' (or wirings') mean entry per odor: each of ENTRY_COLUMNS the NaN-aware mean over them (the
    flies' mean views of the strength parts, as hal's `fly_reference` and `mean_views` take them), the min and max over
    them of the rate and mean entries and of the concentration, the number of recordings defining the rate entry, the
    mean entry over the rate entry of the means; ranked again."""
    names = list(per)
    first = list(per[names[0]].odor)
    stack = pd.concat([per[k].set_index('odor')[list(ENTRY_COLUMNS)] for k in names], keys=names,
                      names=['recording', 'odor'])
    g = stack.groupby(level='odor', sort=False)
    m = g.mean()
    for c in ('rate', 'mean', 'concentration'):
        m[f'{c}_min'] = g[c].min()
        m[f'{c}_max'] = g[c].max()
    m['n_recordings'] = g['rate'].count()
    m = m.reindex(first).rename_axis('odor').reset_index()
    with np.errstate(invalid='ignore', divide='ignore'):
        m['mean_over_rate'] = m['mean'] / m['rate']
    return add_ranks(m)


def sharing_parts(pop, odors: Sequence[str]) -> Dict[str, pd.Series]:
    """The sharing view's two parts of one population on `odors` (hal's population_views; 'a|b' pairs)."""
    from kcq import population_fit as P
    return P.population_views(pop, list(odors), parts=['sharing/who', 'sharing/how_strongly'])


def mean_parts(parts: Sequence[Mapping[str, pd.Series]]) -> Dict[str, pd.Series]:
    """The entry-by-entry NaN-aware mean of several populations' parts (hal's `mean_views`)."""
    from kcq import population_fit as P
    return P.mean_views(list(parts))


def odor_row(part: pd.Series, odor: str = ODOR, sep: str = '|') -> pd.Series:
    """An odor's row of a pair part: its distance to every other odor (indexed by the other odor)."""
    out = {}
    for k, v in part.items():
        a, b = str(k).split(sep)
        if a == odor:
            out[b] = float(v)
        elif b == odor:
            out[a] = float(v)
    return pd.Series(out, dtype=float)


def odor_mean_distance(part: pd.Series, odors: Sequence[str]) -> pd.Series:
    """Each odor's mean distance to the other odors over its pairs (NaN pairs left out)."""
    return pd.Series({o: float(np.nanmean(odor_row(part, o).to_numpy())) if odor_row(part, o).notna().any()
                      else np.nan for o in odors})


def sharing_summary(parts: Mapping[str, pd.Series], odors: Sequence[str]) -> pd.DataFrame:
    """Per part and odor: its mean distance to the other odors, its rank among the odors (ascending: rank 17 the most
    apart) and whether it is among the CLASSES['ii'] top largest."""
    rows = []
    top = int(CLASSES['ii']['top'])
    for p, s in parts.items():
        md = odor_mean_distance(s, odors)
        r = ranks(md.to_numpy())
        n = int(r.notna().sum())
        for o, v, k in zip(md.index, md.to_numpy(), r.to_numpy()):
            rows.append(dict(part=p, odor=o, mean_distance=v, rank=k, n_ranked=n,
                             among_top=bool(np.isfinite(k) and k > n - top)))
    return pd.DataFrame(rows)


def responders_breadth(pop, odors: Sequence[str]) -> pd.DataFrame:
    """Per odor: the number of `odors` each of its responders responds to (mean, median) and the fraction of them in
    each of E134's breadth bins for that many odors; beside, all units' mean breadth."""
    from data_prep import stage_population_data as D
    od = list(odors)
    R = pop.responds.reindex(columns=od).astype(float)
    n = np.nansum(R.to_numpy(), axis=1)
    bins = D.breadth_bins(len(od))
    rows = []
    for o in od:
        m = np.nan_to_num(R[o].to_numpy()) > 0.5
        x = n[m]
        rec = dict(odor=o, n_responders=int(m.sum()), breadth_mean=float(x.mean()) if len(x) else np.nan,
                   breadth_median=float(np.median(x)) if len(x) else np.nan, all_units_breadth_mean=float(n.mean()))
        for lo, hi in bins:
            rec[f'bin_{lo}_{hi}'] = float(np.mean((x >= lo) & (x <= hi))) if len(x) else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    r = ranks(df.breadth_mean.to_numpy())
    df['breadth_mean_rank'] = r.to_numpy()
    df['breadth_mean_third'] = [third(v, int(r.notna().sum())) for v in r]
    return df


def loo_fit(x: pd.Series, y: pd.Series, leave: str = ODOR) -> dict:
    """y = a + b x by least squares over the odors other than `leave` (both defined); `leave`'s residual in y's units
    and over the residual SD of the others (ddof 2); Pearson r and Spearman rho over the others and over all; every
    odor's leave-one-out residual and `leave`'s rank by its absolute value (1 the largest)."""
    from scipy import stats
    d = pd.DataFrame({'x': pd.Series(x, dtype=float), 'y': pd.Series(y, dtype=float)}).dropna()
    out = dict(n=int(len(d)), x_aa=float(d.x.get(leave, np.nan)), y_aa=float(d.y.get(leave, np.nan)))
    o = d.drop(index=[leave], errors='ignore')
    if len(o) < 3 or leave not in d.index or o.x.std() == 0:
        out.update(slope=np.nan, intercept=np.nan, pred_aa=np.nan, resid_aa=np.nan, z_aa=np.nan, r=np.nan, rho=np.nan,
                   r_all=np.nan, rho_all=np.nan, loo_rank_aa=np.nan)
        return out
    b, a = np.polyfit(o.x.to_numpy(), o.y.to_numpy(), 1)
    res_o = o.y - (a + b * o.x)
    sd = float(np.sqrt((res_o ** 2).sum() / max(len(o) - 2, 1)))
    pred = float(a + b * d.x[leave])
    loo = {}
    for k in d.index:
        dk = d.drop(index=[k])
        bk, ak = np.polyfit(dk.x.to_numpy(), dk.y.to_numpy(), 1)
        loo[k] = float(d.y[k] - (ak + bk * d.x[k]))
    lr = pd.Series(loo).abs().rank(ascending=False, method='average')
    out.update(slope=float(b), intercept=float(a), pred_aa=pred, resid_aa=float(d.y[leave] - pred),
               z_aa=float((d.y[leave] - pred) / sd) if sd > 0 else np.nan, resid_sd_others=sd,
               r=float(np.corrcoef(o.x, o.y)[0, 1]), rho=float(stats.spearmanr(o.x, o.y).statistic),
               r_all=float(np.corrcoef(d.x, d.y)[0, 1]), rho_all=float(stats.spearmanr(d.x, d.y).statistic),
               loo_rank_aa=float(lr[leave]), loo_resid_aa=loo[leave])
    return out


def input_features(inp: pd.DataFrame) -> pd.DataFrame:
    """Per odor of an input (glomeruli x odors, spike deltas): the summed input (signed), the sum of its positive
    part, its largest entry and that glomerulus, its breadth (the participation ratio of the positive part), and the
    number of E138's feeders (the glomeruli carrying two thirds of the positive input)."""
    from experiments import E138_strength_by_odor_entries as E138
    rows = []
    for o in inp.columns:
        x = inp[o].astype(float)
        pos = x.clip(lower=0)
        rows.append(dict(odor=o, input_sum=float(x.sum()), input_pos_sum=float(pos.sum()), input_max=float(x.max()),
                         input_argmax=str(x.idxmax()),
                         input_breadth=float(pos.sum() ** 2 / (pos ** 2).sum()) if (pos ** 2).sum() > 0 else np.nan,
                         n_feeders=int(len(E138.feeders(x)))))
    return pd.DataFrame(rows)


def class_i(flies_claws: pd.DataFrame, model_claws: pd.DataFrame, rule: Mapping = CLASSES['i']) -> dict:
    """Class (i) from aa's rows at the claws: `flies_claws` one row per claw recording (columns recording, rate_third,
    concentration_third), `model_claws` one per wiring (wiring, rate_third, concentration_third)."""
    f = flies_claws.assign(ok=(flies_claws.rate_third == 'lowest') & (flies_claws.concentration_third == 'highest'))
    m = model_claws.assign(ok=(model_claws.rate_third == 'highest') & (model_claws.concentration_third == 'lowest'))
    nf, nm = int(f.ok.sum()), int(m.ok.sum())
    holds_f = len(f) == int(rule['flies_n']) and nf >= int(rule['flies_min'])
    holds_m = len(m) == int(rule['model_n']) and nm >= int(rule['model_min'])
    return dict(name=rule['name'], rule=rule['rule'], flies=f.to_dict('records'), model=m.to_dict('records'),
                flies_count=nf, flies_n=int(len(f)), model_count=nm, model_n=int(len(m)), flies_part=holds_f,
                model_part=holds_m, holds=bool(holds_f and holds_m))


def class_ii(halves: Mapping[str, pd.DataFrame], rule: Mapping = CLASSES['ii']) -> dict:
    """Class (ii) from each half's sharing summary (`sharing_summary` rows of part CLASSES['ii']['part'])."""
    per = {}
    for h, s in halves.items():
        x = s[(s.part == rule['part']) & (s.odor == ODOR)]
        per[h] = dict(mean_distance=float(x.mean_distance.iloc[0]), rank=float(x['rank'].iloc[0]),
                      n_ranked=int(x.n_ranked.iloc[0]), among_top=bool(x.among_top.iloc[0]),
                      top_odors=list(s[s.part == rule['part']].sort_values('mean_distance', ascending=False)
                                     .odor[:int(rule['top'])]))
    return dict(name=rule['name'], rule=rule['rule'], halves=per,
                holds=bool(len(per) == 2 and all(v['among_top'] for v in per.values())))


def class_iv(readings: pd.DataFrame, rule: Mapping = CLASSES['iv']) -> dict:
    """Class (iv) from aa's mean rows under each soma reading (`readings`: columns reading, counts, rate, rate_rank,
    rate_n_ranked, rate_min, rate_max, mean, mean_rank)."""
    per = {}
    leaves = []
    for r in readings.itertuples():
        top = r.rate_n_ranked - THIRD
        left = bool(np.isfinite(r.rate_rank) and r.rate_rank <= top)
        per[r.reading] = dict(counts=bool(r.counts), rate=float(r.rate), rate_rank=float(r.rate_rank),
                              rate_min=float(r.rate_min), rate_max=float(r.rate_max), mean=float(r.mean),
                              mean_rank=float(r.mean_rank), leaves_highest_third=left)
        if left and r.counts:
            leaves.append(r.reading)
    return dict(name=rule['name'], rule=rule['rule'], readings=per, leaving=leaves, holds=bool(leaves))


def outcomes_of(ci: bool, ciii: Optional[bool], civ: bool) -> dict:
    """The entry's outcomes the classes point to: (d) where (iv) holds; (c) where (i) does not; where (i) holds, (b) if
    (iii) holds, (a) if it does not, and 'step 1' while (iii) is unread (None)."""
    out = []
    if civ:
        out.append('d')
    if not ci:
        out.append('c')
    elif ciii is None:
        out.append('a or b (step 1 reads (iii))')
    else:
        out.append('b' if ciii else 'a')
    return dict(outcomes=out, step1_called=bool(ci), texts={k: OUTCOMES[k] for k in OUTCOMES if k in ''.join(out)})


def _tag(df: pd.DataFrame, **cols) -> pd.DataFrame:
    df = df.copy()
    for i, (k, v) in enumerate(cols.items()):
        df.insert(i, k, v)
    return df


def _check_file(entry: Mapping) -> Path:
    """The lock's file `entry` (path, sha256) after checking its sha256; ValueError if it changed."""
    p = Path(entry['path'])
    got = _sha256(p)
    if got != entry['sha256']:
        raise ValueError(f"{p}: sha256 {got[:12]} differs from the lock's {entry['sha256'][:12]}")
    return p


def provenance(part: str, lock: Mapping, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time, the lock's sha256 and `extra`."""
    _write_json(out_dir() / f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
                     lock_sha256=lock.get('sha256'), **extra))


# ---------------------------------------------------------------------------
# reduce: the soma flies' half windows

def reduce_windows(trials_path: Path, respvec_path: Path) -> Tuple[pd.DataFrame, dict]:
    """One megamat soma fly's frames reduced as data_prep/stage_population_data.py's `reduce_megamat_trials` reduces
    the whole window (her good cells, the pre-odor window and its frame SD), for the whole window and each of
    HALF_WINDOWS: long rows cell, trial, stim, sd, dfc_full, dfc_early, dfc_late."""
    import xarray as xr
    rv = xr.open_dataset(respvec_path)
    good = set(int(x) for x in np.atleast_1d(rv.attrs['good_xid']))
    keep_rv = np.isin(rv['xid0'].values, list(good)) & (rv['iscell'].values == 1)
    good_cells = set(int(c) for c in rv['cells'].values[keep_rv])
    rv.close()
    tr = xr.open_dataset(trials_path)
    t = tr['time'].values
    cells = tr['cells'].values
    sel = np.array([int(c) in good_cells for c in cells])
    Fc = tr['Fc'].values[:, sel]
    stims = [str(s) for s in tr['stim'].values]
    tr.close()

    def win(lo, hi):
        return (t >= lo) & (t <= hi)
    pre = win(*PRE_WINDOW)
    base = Fc[:, :, pre].mean(axis=2)
    sd = np.sqrt(np.mean(np.std(Fc[:, :, pre], axis=2) ** 2, axis=0))
    cs = cells[sel]
    cols = dict(cell=np.tile(cs, len(stims)), trial=np.repeat(np.arange(len(stims)), len(cs)),
                stim=np.repeat(stims, len(cs)), sd=np.tile(sd, len(stims)))
    masks = dict(full=win(*FULL_WINDOW), **{k: win(*v) for k, v in HALF_WINDOWS.items()})
    for k, m in masks.items():
        cols[f'dfc_{k}'] = (Fc[:, :, m].mean(axis=2) - base).reshape(-1)
    info = dict(n_cells=int(sel.sum()), n_frames={k: int(m.sum()) for k, m in masks.items()}, n_frames_pre=int(pre.sum()),
                frame_times={k: [float(t[m].min()), float(t[m].max())] for k, m in masks.items()},
                halves_partition_full=bool(np.array_equal(masks['full'], masks['early'] | masks['late'])
                                           and not (masks['early'] & masks['late']).any()))
    return pd.DataFrame(cols), info


def part_reduce() -> dict:
    """The four megamat soma flies' frames (Remy's suite2p trials files on matrix, sha256 checked against the lock's,
    read only) reduced to the whole window and its two halves (`reduce_windows`); the whole window must reproduce
    E134's frames_reduced.parquet (dfc and sd, to 1e-6 relative); soma_windows.parquet and soma_windows.json."""
    lock, _ = require_all()
    t0 = time.time()
    ref = pd.read_parquet(_check_file(lock['inputs']['frames_reduced']))
    parts, info = [], {}
    for fly, f in lock['inputs']['soma_frames'].items():
        tp = _check_file(f['trials'])
        rp = _check_file(f['respvec'])
        d, inf = reduce_windows(tp, rp)
        d.insert(0, 'fly', fly)
        r = ref[ref.fly == fly].sort_values(['trial', 'cell'])
        dd = d.sort_values(['trial', 'cell'])
        same = (len(r) == len(dd) and np.array_equal(r.cell.to_numpy(), dd.cell.to_numpy())
                and np.array_equal(r.trial.to_numpy(), dd.trial.to_numpy()))
        scale = max(1.0, float(np.nanmax(np.abs(r.dfc.to_numpy())))) if len(r) else 1.0
        diff = float(np.nanmax(np.abs(r.dfc.to_numpy() - dd.dfc_full.to_numpy()))) if same else np.inf
        dsd = float(np.nanmax(np.abs(r.sd.to_numpy() - dd.sd.to_numpy()))) if same else np.inf
        inf.update(same_rows=bool(same), max_abs_diff_dfc=diff, max_abs_diff_sd=dsd,
                   reproduces=bool(same and diff <= 1e-6 * scale and dsd <= 1e-6 * max(1.0, float(r.sd.max()))))
        info[fly] = inf
        print(f'{fly}: {inf["n_cells"]} cells, frames {inf["n_frames"]}, reproduces {inf["reproduces"]} '
              f'({diff:.2e}), {time.time() - t0:.0f} s', flush=True)
        if not inf['reproduces']:
            raise RuntimeError(f"{fly}: the whole window does not reproduce E134's frames_reduced.parquet")
        if not inf['halves_partition_full']:
            raise RuntimeError(f'{fly}: the two halves do not partition the whole window')
        parts.append(d)
    out = out_dir() / 'soma_windows.parquet'
    tmp = out.with_name(f'.{out.name}.{_host()}.{os.getpid()}.tmp')
    pd.concat(parts, ignore_index=True).to_parquet(tmp, index=False)
    os.replace(tmp, out)
    summ = dict(file=file_entry(out), flies=info, windows=dict(full=FULL_WINDOW, pre=PRE_WINDOW, halves=HALF_WINDOWS),
                seconds=round(time.time() - t0, 1))
    _write_json(out_dir() / 'soma_windows.json', summ)
    provenance('reduce', lock, soma_windows=summ['file'])
    return summ


# ---------------------------------------------------------------------------
# the flies (design item 1)

def recorded(stage: str, cohort: str, lock: Mapping, trials: bool = False):
    """viki's recorded population of a stage and cohort on megamat (its sha256 checked against E134's lock), with its
    flies checked against this lock's."""
    from data_prep import stage_population_data as D
    pop = D.load_population(stage, PANEL, cohort, trials=trials)
    want = lock['recordings'][f'{stage}__{cohort}']['flies']
    if list(pop.flies) != list(want):
        raise RuntimeError(f'{stage} {cohort}: flies {list(pop.flies)} differ from the lock\'s {want}')
    return pop


def remy_binary_pops(lock: Mapping) -> List:
    """Remy's own responder calls of the four soma flies (her df_stim_responders_from_peak.csv, files checked against
    the lock), as binary populations named by fly ('date/fly_num', through her df_acqs.csv), odors abbreviated."""
    from kcq import population_fit as P
    from pyolfsysm import util as U
    for k, e in lock['inputs']['remy_binary'].items():
        _check_file(e)
    b = U.load_remy_megamat_kc_binary_responses()
    acqs = pd.read_csv(Path(U.remy_sparsity_dir) / 'df_acqs.csv')
    name = {int(r.acq): f'{pd.Timestamp(r.date_imaged):%Y-%m-%d}/{int(r.fly_num)}' for r in acqs.itertuples()}
    out = {}
    for acq, x in b.groupby(level=U.remy_fly_id, sort=True):
        v = x.droplevel(U.remy_fly_id).astype(float)
        v.index = v.index.astype(str)
        out[name[int(acq)]] = P.Population(v, v, name=name[int(acq)], graded=False)
    return [out[f] for f in lock['recordings']['kc_soma__wild_type']['flies']]


def soma_reading_pops(rd: SomaReading, lock: Mapping) -> Tuple[List, dict]:
    """The four soma flies' populations under one of SOMA_READINGS, and its criteria per fly."""
    from data_prep import stage_population_data as D
    from kcq import population_fit as P
    odors = list(lock['odors'])
    if rd.source == 'remy_binary':
        pops = remy_binary_pops(lock)
        return pops, dict(criterion='her responder calls', n_units={p.name: int(len(p.values)) for p in pops})
    if rd.source.startswith('window:'):
        w = rd.source.split(':', 1)[1]
        p = out_dir() / 'soma_windows.parquet'
        summ = json.loads((out_dir() / 'soma_windows.json').read_text())
        if _sha256(p) != summ['file']['sha256']:
            raise ValueError(f'{p.name} differs from soma_windows.json')
        r = pd.read_parquet(p)
        r = r[np.isfinite(r.sd) & (r.sd > 0)]
        r = r.assign(odor=r.stim.str.split(' @ ').str[0], value=r[f'dfc_{w}'] / r.sd)
        r = r[r.odor.isin(odors)]
        tables = {f: x.groupby(['cell', 'odor']).value.mean().unstack('odor').reindex(columns=odors)
                  for f, x in r.groupby('fly', sort=False)}
        for t in tables.values():
            t.index = t.index.astype(str)
        tables = {f: tables[f] for f in lock['recordings']['kc_soma__wild_type']['flies']}
    else:
        pop = recorded('kc_soma', 'wild_type', lock, trials=rd.source.startswith('presentation:'))
        if rd.source.startswith('presentation:'):
            k = int(rd.source.split(':', 1)[1])
            tables = {f: pop.trials[f].xs(k, axis=1, level='repeat').reindex(columns=odors) for f in pop.flies}
        else:
            tables = {f: pop.tables[f].reindex(columns=odors) for f in pop.flies}
            if not np.isclose(float(pop.criterion), LOCK_SOMA_CRITERION):
                raise RuntimeError(f"the soma population's criterion {pop.criterion} is not the lock's")
    dec = int(D.CRITERION_SEARCH['decimals'])
    if rd.criterion == 'lock':
        crit = {f: LOCK_SOMA_CRITERION for f in tables}
    elif rd.criterion == 'apl16':
        crit = {f: APL16_SOMA_CRITERION for f in tables}
    elif rd.criterion == 'rate':
        c = round(D.criterion_at_rate(list(tables.values()), SOMA_RATE), dec)
        crit = {f: c for f in tables}
    elif rd.criterion == 'per_fly':
        crit = {f: round(D.criterion_at_rate([t], SOMA_RATE), dec) for f, t in tables.items()}
    else:
        raise ValueError(rd.criterion)
    pops = [P.Population.from_criterion(t, crit[f], name=f) for f, t in tables.items()]
    rates = {f: D.fly_rate(t, crit[f]) for f, t in tables.items()}
    return pops, dict(criterion=crit, response_rate=rates, response_rate_mean=float(np.mean(list(rates.values()))))


def part_flies() -> dict:
    """Design item 1 (the flies, no simulation): flies_entries.csv (every recording's per-odor entries, concentration
    and ranks at every stage, and the flies' mean rows, the five claw recordings' mean beside), flies_soma_readings.csv
    (the somata under every reading of the lock), flies_sharing.csv and flies_sharing_summary.csv (aa's row of the
    sharing view, and every odor's mean distance and rank, at every stage, each recording, all flies and the soma
    halves), flies_breadth.csv (the number of odors each odor's somatic responders respond to), flies_same_fly_step.csv
    (the same-fly step and the four flies' means beside)."""
    from kcq import population_fit as P
    lock, _ = require_all()
    t0 = time.time()
    odors = list(lock['odors'])
    rows, per, fps_of = [], {}, {}
    for key, rec in lock['recordings'].items():
        if '__' not in key:
            continue
        pop = recorded(rec['stage'], rec['cohort'], lock)
        fps = P.fly_populations(pop)
        d = {fp.name: odor_stats(fp, odors) for fp in fps}
        for n, s in d.items():
            rows.append(_tag(s, stage=rec['stage'], cohort=rec['cohort'], reading='lock', recording=n))
        if len(d) >= 2:
            rows.append(_tag(mean_rows(d), stage=rec['stage'], cohort=rec['cohort'], reading='lock', recording='mean'))
        per[(rec['stage'], rec['cohort'])] = d
        fps_of[(rec['stage'], rec['cohort'])] = fps
    five = {**per[('kc_claw', 'wild_type')], **per[('kc_claw', 'same_fly')]}
    rows.append(_tag(mean_rows(five), stage='kc_claw', cohort='five', reading='lock', recording='mean'))
    entries = pd.concat(rows, ignore_index=True)
    _write_csv(entries, 'flies_entries.csv')
    print(f'flies: entries ({time.time() - t0:.0f} s)', flush=True)

    # the somata under every reading of the lock
    srows, sinfo = [], {}
    lock_pops = fps_of[('kc_soma', 'wild_type')]
    for rd in SOMA_READINGS:
        pops, info = soma_reading_pops(rd, lock)
        d = {p.name: odor_stats(p, odors) for p in pops}
        if rd.source == 'remy_binary':
            # the rate part only: the mean part, the mean response and the concentration are the lock's reading's
            for p in lock_pops:
                s = odor_stats(p, odors).set_index('odor')
                x = d[p.name].set_index('odor')
                for c in ('mean', 'mean_response'):
                    x[c] = s[c]
                d[p.name] = add_ranks(x.reset_index())
        for n, s in d.items():
            srows.append(_tag(s, reading=rd.name, counts=rd.counts, recording=n,
                              criterion=json.dumps(info.get('criterion', {}).get(n) if isinstance(info.get('criterion'),
                                                                                                   dict)
                                                   else info.get('criterion'), default=_json_default)))
        srows.append(_tag(mean_rows(d), reading=rd.name, counts=rd.counts, recording='mean', criterion=''))
        sinfo[rd.name] = dict(info, what=rd.what, counts=rd.counts)
        print(f'flies: soma reading {rd.name} ({time.time() - t0:.0f} s)', flush=True)
    _write_csv(pd.concat(srows, ignore_index=True), 'flies_soma_readings.csv')

    # who responds with aa: the sharing view at every stage
    sh_rows, sh_sum = [], []
    for (stage, cohort), fps in fps_of.items():
        parts = {fp.name: sharing_parts(fp, odors) for fp in fps}
        sets = dict(parts)
        if len(parts) >= 2:
            sets['all'] = mean_parts(list(parts.values()))
        if (stage, cohort) == ('kc_soma', 'wild_type'):
            for h, names in lock['definitions']['soma_halves'].items():
                sets[f'half_{h}'] = mean_parts([parts[n] for n in names])
        for sname, pp in sets.items():
            for part, s in pp.items():
                r = odor_row(s)
                sh_rows += [dict(stage=stage, cohort=cohort, set=sname, part=part, odor=ODOR, other=o, distance=v)
                            for o, v in r.items()]
            sh_sum.append(_tag(sharing_summary(pp, odors), stage=stage, cohort=cohort, set=sname))
    _write_csv(pd.DataFrame(sh_rows), 'flies_sharing.csv')
    _write_csv(pd.concat(sh_sum, ignore_index=True), 'flies_sharing_summary.csv')

    # the number of odors aa's somatic responders respond to
    br = []
    for cohort in ('wild_type', 'same_fly'):
        d = {fp.name: responders_breadth(fp, odors) for fp in fps_of[('kc_soma', cohort)]}
        for n, s in d.items():
            br.append(_tag(s, cohort=cohort, recording=n))
        if len(d) >= 2:
            m = pd.concat(d.values()).groupby('odor', sort=False).mean(numeric_only=True).reindex(odors).reset_index()
            r = ranks(m.breadth_mean.to_numpy())
            m['breadth_mean_rank'] = r.to_numpy()
            m['breadth_mean_third'] = [third(v, int(r.notna().sum())) for v in r]
            br.append(_tag(m, cohort=cohort, recording='mean'))
    _write_csv(pd.concat(br, ignore_index=True), 'flies_breadth.csv')

    # the same-fly step, and the four flies' means beside
    def step(a: pd.DataFrame, b: pd.DataFrame, label: str) -> pd.DataFrame:
        a, b = a.set_index('odor'), b.set_index('odor')
        x = pd.DataFrame({'claws_rate': a.rate, 'somata_rate': b.rate, 'claws_mean': a['mean'], 'somata_mean': b['mean'],
                          'claws_concentration': a.concentration, 'somata_concentration': b.concentration})
        x['step_rate'] = x.somata_rate - x.claws_rate
        x['step_mean'] = x.somata_mean - x.claws_mean
        x = x.reindex(odors).rename_axis('odor').reset_index()
        for c in ('claws_rate', 'somata_rate', 'step_rate', 'step_mean'):
            r = ranks(x[c].to_numpy())
            x[f'{c}_rank'] = r.to_numpy()
        return _tag(x, set=label)
    sf_claw = per[('kc_claw', 'same_fly')][lock['recordings']['kc_claw__same_fly']['flies'][0]]
    sf_soma = per[('kc_soma', 'same_fly')][lock['recordings']['kc_soma__same_fly']['flies'][0]]
    mean_of = {k: e for k, e in ((k, entries[(entries.stage == k[0]) & (entries.cohort == k[1])
                                             & (entries.recording == 'mean')]) for k in
                                 (('kc_claw', 'wild_type'), ('kc_soma', 'wild_type'), ('kc_claw', 'five')))}
    sfs = pd.concat([step(sf_claw, sf_soma, 'same_fly'),
                     step(mean_of[('kc_claw', 'wild_type')], mean_of[('kc_soma', 'wild_type')], 'four_flies_means'),
                     step(mean_of[('kc_claw', 'five')], mean_of[('kc_soma', 'wild_type')], 'five_claws_four_somata')],
                    ignore_index=True)
    _write_csv(sfs, 'flies_same_fly_step.csv')
    provenance('flies', lock, soma_readings=sinfo, seconds=round(time.time() - t0, 1),
               soma_windows=json.loads((out_dir() / 'soma_windows.json').read_text())['file'])
    print(f'flies: done ({time.time() - t0:.0f} s)', flush=True)
    return dict(entries=entries)


# ---------------------------------------------------------------------------
# the saved megamat runs (design item 2)

def group_runs(g: ModelGroup) -> Tuple[Dict[str, object], List[str]]:
    """A model group's saved runs ({wiring: SavedRun}, WIRINGS order) and the tags not found."""
    from kcq import stage_readings as R
    exp, part, pattern = g.source.split(':', 2)
    if exp == 'E134':
        dirs = R.step3_dirs()
    elif exp == 'E136':
        # the arms' reads under the adopted measure (strength by odor), those E136's note reports: the archive's
        # measure_strength/read/ (E136's default measure 'all' points at a folder with no reads)
        from experiments import E136_stage_population_fit as E136
        E136.set_measure(E136_MEASURE)
        dirs = E136.read_dirs()
    elif exp == 'E135':
        from experiments import E135_antennal_lobe_dendrites as E135
        dirs = E135.step2_dirs()
    else:
        raise ValueError(g.source)
    runs, missing = {}, []
    for w in g.wirings:
        tag = pattern.format(w=w)
        try:
            runs[w] = R.load_step3_run(tag, dirs=dirs)
        except FileNotFoundError:
            missing.append(tag)
    return runs, missing


def e139_runs() -> Dict[str, Dict[str, object]]:
    """E139's saved megamat runs, when they arrive: {setting tag ('<setting>_inhadd<rung>'): {wiring: SavedRun}} from
    outputs/E139_normalization_ladder/read/ and the archive's read/ (stage_readings' step-3 layout, E139's
    `point_tag` '<setting>_inhadd<rung>_megamat_<wiring>'); {} until then."""
    from kcq import stage_readings as R
    out: Dict[str, Dict[str, object]] = {}
    dirs = R._archive_dirs(E139_EXPERIMENT, 'read')
    for d in dirs:
        if not d.is_dir():
            continue
        for p in sorted(d.glob('*_megamat_*_run.json')):
            tag = p.name[:-len('_run.json')]
            w = tag.rsplit('_', 1)[-1]
            setting = tag[:-len(f'_megamat_{w}')]
            if w not in WIRINGS or w in out.get(setting, {}):
                continue
            try:
                out.setdefault(setting, {})[w] = R.load_step3_run(tag, dirs=[d])
            except (FileNotFoundError, RuntimeError, KeyError) as e:
                print(f'E139 {tag}: not read ({e})', flush=True)
    return {s: {w: r[w] for w in WIRINGS if w in r} for s, r in out.items()}


def read_runs(runs: Mapping[str, object], e134: Mapping, flies, stages: Optional[Sequence[str]] = None,
              readings: Optional[Mapping] = None) -> Dict[str, object]:
    """Each run read as E134 and E136 read a wiring (atlas's `read_run`, the stages' default readings unless
    `readings`, seed index the wiring's position in WIRINGS, E136's `seed_index`)."""
    from kcq import stage_readings as R
    return {w: R.read_run(r, stages=stages, readings=readings, seed_index=WIRINGS.index(w), lock=dict(e134), flies=flies)
            for w, r in runs.items()}


def anatomy(run, condition: str, odors: Sequence[str], ir: Sequence[str]) -> pd.DataFrame:
    """Per odor of one run and condition: its responders (two or more spikes), their mean share of claws on an IR-fed
    glomerulus against all KCs', and per KC type the share of the responders of that type, the type's share of all KCs
    and the type's response rate."""
    from kcq import stage_readings as R
    X = R._abbrev(run.conditions[condition]['spike_counts'].astype(float))[list(odors)]
    kc = np.asarray(X.index.get_level_values(H.KC_ID_LEVEL)).astype(np.int64)
    types = np.asarray(X.index.get_level_values(H.KC_TYPE_LEVEL)).astype(str)
    cg = np.asarray(run.meta['claw_glomerulus']).astype(str)
    ck = np.asarray(run.meta['kc_of_claw']).astype(np.int64)
    claws = pd.DataFrame(dict(kc=ck, ir=np.isin(cg, list(ir)))).groupby('kc').ir.agg(['sum', 'count'])
    share = (claws['sum'] / claws['count']).reindex(kc).to_numpy(dtype=float)
    ut = sorted(set(types))
    rows = []
    for o in odors:
        m = X[o].to_numpy() >= 2.0
        rec = dict(odor=o, n_responders=int(m.sum()), response_rate=float(m.mean()),
                   ir_claw_share_responders=float(np.nanmean(share[m])) if m.any() else np.nan,
                   ir_claw_share_all=float(np.nanmean(share)),
                   n_claws_responders=float(claws['count'].reindex(kc).to_numpy()[m].mean()) if m.any() else np.nan)
        for t in ut:
            tm = types == t
            rec[f'type_{t}_share_responders'] = float(np.mean(tm[m])) if m.any() else np.nan
            rec[f'type_{t}_share_all'] = float(tm.mean())
            rec[f'type_{t}_response_rate'] = float(np.mean(m[tm])) if tm.any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def part_models() -> dict:
    """Design item 2 (the saved megamat runs, no simulation): models_entries.csv (each group's per-odor entries,
    concentration and ranks at every stage, APL on and off, each wiring and the wirings' mean 'common', in the
    default reading; the family held also in MODEL_SOMA_BESIDE and MODEL_CLAW_BESIDE), models_sharing.csv and
    models_sharing_summary.csv (aa's row of the sharing view), models_breadth.csv (the number of odors each odor's
    somatic responders respond to, in the default reading), models_anatomy.csv (each odor's somatic responders by
    the glomeruli of their claws and by KC type), models_a1_stages.csv (each odor's entry at each stage under A1 beside
    the flies'), and the groups not found."""
    from kcq import stage_readings as R
    lock, e134 = require_all()
    t0 = time.time()
    odors = list(lock['odors'])
    ir = lock['model']['anatomy']['ir_glomeruli']
    flies = R.FlyCache()
    groups = [(g.name, g.source) for g in MODEL_GROUPS]
    ent, sh, shs, an, brd, missing, read_info = [], [], [], [], [], {}, {}
    e139 = e139_runs()
    if not e139:
        missing['e139'] = [f'{E139_EXPERIMENT}: no saved megamat runs yet']

    def entries_of(name: str, readings: Mapping[str, object], reading_of: Callable[[str], str],
                   stages: Sequence[str], with_sharing: bool):
        for stage in stages:
            for cond in R.CONDITIONS:
                pops = {w: m.populations.get((stage, PANEL, cond)) for w, m in readings.items()}
                pops = {w: p for w, p in pops.items() if p is not None}
                if not pops:
                    continue
                d = {w: odor_stats(p, odors) for w, p in pops.items()}
                tags = dict(group=name, condition=cond, stage=stage, reading=reading_of(stage))
                for w, s in d.items():
                    ent.append(_tag(s, **tags, wiring=w))
                if len(d) >= 2:
                    ent.append(_tag(mean_rows(d), **tags, wiring='common'))
                if with_sharing and stage == 'kc_soma':
                    b = {w: responders_breadth(p, odors) for w, p in pops.items()}
                    for w, x in b.items():
                        brd.append(_tag(x, **tags, wiring=w))
                    if len(b) >= 2:
                        m = pd.concat(b.values()).groupby('odor', sort=False).mean(numeric_only=True).reindex(odors)
                        r = ranks(m.breadth_mean.to_numpy())
                        m['breadth_mean_rank'] = r.to_numpy()
                        m['breadth_mean_third'] = [third(v, int(r.notna().sum())) for v in r]
                        brd.append(_tag(m.reset_index(), **tags, wiring='mean'))
                if with_sharing:
                    parts = {w: sharing_parts(p, odors) for w, p in pops.items()}
                    sets = dict(parts, **({'common': mean_parts(list(parts.values()))} if len(parts) >= 2 else {}))
                    for sname, pp in sets.items():
                        for part, s in pp.items():
                            sh.extend(dict(**tags, wiring=sname, part=part, odor=ODOR, other=o, distance=v)
                                      for o, v in odor_row(s).items())
                        shs.append(_tag(sharing_summary(pp, odors), **tags, wiring=sname))

    for name, runs in [(g.name, group_runs(g)) for g in MODEL_GROUPS] + [(f'e139_{s}', (r, [])) for s, r in
                                                                       e139.items()]:
        runs, miss = runs
        if miss:
            missing[name] = miss
        if not runs:
            print(f'models: {name}: no runs', flush=True)
            continue
        rd = read_runs(runs, e134, flies)
        read_info[name] = dict(wirings=list(runs), files=[dict(path=f['path'], sha256=f.get('sha256'))
                                                          for r in runs.values() for f in r.files],
                               skipped={w: {f'{k[0]}|{k[1]}': v for k, v in m.skipped.items()} for w, m in rd.items()})
        entries_of(name, rd, lambda st: R.DEFAULT_READINGS[st].name, R.STAGES, True)
        for w, r in runs.items():
            for cond in R.CONDITIONS:
                an.append(_tag(anatomy(r, cond, odors, ir), group=name, condition=cond, wiring=w))
        if name == 'family_held':
            for st, names in (('kc_soma', MODEL_SOMA_BESIDE), ('kc_claw', MODEL_CLAW_BESIDE)):
                for nm in names:
                    rb = read_runs(runs, e134, flies, stages=[st], readings={st: R.reading_by_name(st, nm)})
                    entries_of(name, rb, lambda _s, nm=nm: nm, [st], False)
        print(f'models: {name} read ({len(runs)} wirings, {time.time() - t0:.0f} s)', flush=True)
        del rd
    entries = pd.concat(ent, ignore_index=True)
    _write_csv(entries, 'models_entries.csv')
    _write_csv(pd.DataFrame(sh), 'models_sharing.csv')
    _write_csv(pd.concat(shs, ignore_index=True), 'models_sharing_summary.csv')
    _write_csv(pd.concat(brd, ignore_index=True), 'models_breadth.csv')
    anat = pd.concat(an, ignore_index=True)
    am = anat.groupby(['group', 'condition', 'odor'], sort=False).mean(numeric_only=True).reset_index()
    _write_csv(pd.concat([anat, am.assign(wiring='mean')], ignore_index=True), 'models_anatomy.csv')
    # each odor's entry at each stage under A1, beside the flies' and the family held's
    fl = pd.read_csv(out_dir() / 'flies_entries.csv')
    fm = fl[(fl.recording == 'mean') & (fl.cohort == 'wild_type')].set_index(['stage', 'odor'])
    rows = []
    for st in R.STAGES:
        for o in odors:
            rec = dict(stage=st, odor=o)
            for c in ('rate', 'mean'):
                rec[f'flies_{c}'] = fm[c].get((st, o), np.nan)
                rec[f'flies_{c}_min'] = fm[f'{c}_min'].get((st, o), np.nan)
                rec[f'flies_{c}_max'] = fm[f'{c}_max'].get((st, o), np.nan)
            for g in ('A1', 'family_held'):
                x = entries[(entries.group == g) & (entries.wiring == 'common') & (entries.condition == 'on')
                            & (entries.stage == st) & (entries.odor == o)
                            & (entries.reading == R.DEFAULT_READINGS[st].name)]
                for c in ('rate', 'mean'):
                    rec[f'{g}_{c}'] = float(x[c].iloc[0]) if len(x) else np.nan
            rows.append(rec)
    _write_csv(pd.DataFrame(rows), 'models_a1_stages.csv')
    provenance('models', lock, e134_lock_sha256=e134.get('sha256'), groups=read_info, missing=missing,
               seconds=round(time.time() - t0, 1))
    print(f'models: done ({time.time() - t0:.0f} s); missing {missing}', flush=True)
    return dict(entries=entries, missing=missing)


# ---------------------------------------------------------------------------
# what each side's step follows (design item 3)

def megamat_input() -> pd.DataFrame:
    """The model's megamat input on its 36 glomeruli, odors abbreviated (E134's and E138's `panel_input`)."""
    from experiments import E138_strength_by_odor_entries as E138
    return E138.panel_input(PANEL)


def part_follows() -> dict:
    """Design item 3 (no simulation): over the 17 odors, each side's KC-soma rate entry (and mean entry) against the
    input's sum, largest entry and breadth and against the same side's upstream stages, the line fitted with aa left
    out (`loo_fit`); follows.csv, follows_points.csv."""
    from kcq import stage_readings as R
    lock, _ = require_all()
    odors = list(lock['odors'])
    feats = input_features(megamat_input()[odors]).set_index('odor')
    fl = pd.read_csv(out_dir() / 'flies_entries.csv')
    me = pd.read_csv(out_dir() / 'models_entries.csv')
    sides = {'flies': fl[(fl.recording == 'mean') & (fl.cohort == 'wild_type') & (fl.reading == 'lock')]}
    for g in me.group.unique():
        x = me[(me.group == g) & (me.wiring == 'common') & (me.condition == 'on')]
        sides[g] = x[[r == R.DEFAULT_READINGS[s].name for s, r in zip(x.stage, x.reading)]]
    rows, pts = [], []
    for side, e in sides.items():
        by = {st: e[e.stage == st].set_index('odor') for st in R.STAGES}
        if not len(by['kc_soma']):
            continue
        for part in ('rate', 'mean'):
            y = by['kc_soma'][part].reindex(odors)
            xs = {f: feats[f].reindex(odors) for f in FOLLOW_INPUT_FEATURES}
            xs.update({f'stage_{st}': by[st][part].reindex(odors) for st in FOLLOW_STAGES if len(by[st])})
            for xn, x in xs.items():
                rows.append(dict(side=side, y=f'kc_soma_{part}', x=xn, **loo_fit(x, y)))
                pts += [dict(side=side, y=f'kc_soma_{part}', x=xn, odor=o, x_value=float(x.get(o, np.nan)),
                             y_value=float(y.get(o, np.nan))) for o in odors]
    out = pd.DataFrame(rows)
    _write_csv(out, 'follows.csv')
    _write_csv(pd.DataFrame(pts), 'follows_points.csv')
    _write_csv(feats.reset_index(), 'input_features.csv')
    provenance('follows', lock, input_sha256=R.array_sha256(megamat_input()[odors].to_numpy()))
    return dict(follows=out)


# ---------------------------------------------------------------------------
# the input beside its measurements (design item 4)

def orn_flies(lock: Mapping) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The megamat ORN flies by E85's map (olfdb dataset 1, E134's copy, sha256 checked): each fly's trial-mean dF/F
    scaled 'to-avg-max' (its largest |trial mean| over all its panels and ROIs brought to the lab's avg_flymax), times
    the least-squares slope through the origin of the megamat input (36 glomeruli x 17 odors) on the flies' mean.
    Returns (aa's vectors on the 36 glomeruli, flies as columns, NaN where not measured; aa's vectors on every recorded
    ROI; info: the slope, its R^2, E85's slope, the flies given aa)."""
    from experiments import E85_orn_fly_population as E85
    df = pd.read_csv(_check_file(lock['inputs']['olfdb_dataset_1']))
    odors = list(lock['odors'])
    inp = megamat_input()[odors]
    gl = [str(g) for g in inp.index]
    df = df.assign(fly=[f'{d}/{int(n)}' for d, n in zip(df.date, df.fly_num)])
    tm = df.groupby(['fly', 'panel', 'odor1', 'roi']).value.mean()
    fly_max = tm.abs().groupby(level='fly').max()
    A = E85.avg_flymax()
    want = {f'{o} @ {-3:g}': o for o in odors}
    m = df[(df.panel == PANEL) & df.odor1.isin(want)]
    m = m.assign(odor=m.odor1.map(want))
    fm = m.groupby(['fly', 'odor', 'roi']).value.mean()
    fm = fm / fly_max.reindex(fm.index.get_level_values('fly')).to_numpy() * A
    M = fm.groupby(level=['odor', 'roi']).mean().unstack('odor').reindex(index=gl, columns=odors)
    ok = M.notna().to_numpy()
    a, b = M.to_numpy()[ok], inp.to_numpy(dtype=float)[ok]
    slope = float((a @ b) / (a @ a))
    res = b - slope * a
    r2o = float(1.0 - (res @ res) / (b @ b))
    aa = fm.xs(ODOR, level='odor').unstack('fly') * slope
    given = sorted(aa.columns)
    vec36 = aa.reindex(gl)
    e85 = json.loads((H.shared_dir() / 'results' / 'E85_orn_fly_population' / 'part_a.json').read_text())
    mean_aa = M[ODOR] * slope
    info = dict(slope=slope, r2_origin=r2o, e85_slope=float(e85['scaling']['slope']), avg_flymax=A,
                flies=sorted(fly_max.index), flies_given_aa=given, fly_max={k: float(v) for k, v in fly_max.items()},
                aa_mean_vs_input=dict(pearson=float(np.corrcoef(mean_aa.dropna(), inp[ODOR][mean_aa.notna()])[0, 1]),
                                      ratio_of_sums=float(inp[ODOR][mean_aa.notna()].sum() / mean_aa.dropna().sum())))
    return vec36, aa, info


def vector_summary(x: pd.Series, ref: Optional[pd.Series], ir: Sequence[str], top: int = 3) -> dict:
    """An input vector's sum, positive sum, largest entry and its glomerulus, top glomeruli, breadth (participation
    ratio of the positive part), share of the positive part on the IR-fed class, the number at or above
    HALLEM_RESPONDING_HZ, and its Pearson r and Spearman rho with `ref` on the glomeruli both have."""
    from scipy import stats
    x = pd.Series(x, dtype=float).dropna()
    pos = x.clip(lower=0)
    out = dict(n_glomeruli=int(len(x)), sum=float(x.sum()), pos_sum=float(pos.sum()),
               max=float(x.max()) if len(x) else np.nan, argmax=str(x.idxmax()) if len(x) else '',
               top=';'.join(x.sort_values(ascending=False).index[:top].astype(str)),
               breadth=float(pos.sum() ** 2 / (pos ** 2).sum()) if (pos ** 2).sum() > 0 else np.nan,
               ir_share=float(pos[[g for g in pos.index if g in set(ir)]].sum() / pos.sum()) if pos.sum() > 0 else np.nan,
               n_at_or_above_hz=int((x >= HALLEM_RESPONDING_HZ).sum()))
    if ref is not None:
        r = pd.Series(ref, dtype=float).reindex(x.index)
        ok = r.notna()
        if ok.sum() >= 3:
            out.update(pearson_with_ref=float(np.corrcoef(x[ok], r[ok])[0, 1]),
                       spearman_with_ref=float(stats.spearmanr(x[ok], r[ok]).statistic),
                       ref_sum_same_glomeruli=float(r[ok].sum()), n_shared_with_ref=int(ok.sum()))
    return out


def part_input() -> dict:
    """Design item 4 (no run): aa's model input on the 36 glomeruli beside each ORN fly given aa (E85's map), Hallem's
    aa on the glomeruli it shares with them (drosolf's table, deltas), the PN-dendrite recordings on E138's 28 glomeruli
    and the family held's PN stage on them; input_aa.csv (one row per glomerulus) and input_summary.csv (one row per
    vector)."""
    from drosolf import orns
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    lock, _ = require_all()
    odors = list(lock['odors'])
    ir = lock['model']['anatomy']['ir_glomeruli']
    d28 = lock['input']['dendrite28']
    inp = megamat_input()[odors][ODOR]
    vec36, vec_all, info = orn_flies(lock)
    _check_file(lock['inputs']['hallem'])
    hal = orns.orns(add_sfr=False, columns='glomerulus').loc[HALLEM_ODOR].astype(float)
    hal.index = hal.index.astype(str)
    hal36 = hal[[g for g in hal.index if g in set(inp.index)]]
    den = recorded('pn_dendrite', 'wild_type', lock)
    dvec = pd.DataFrame({f: den.tables[f][ODOR] for f in den.flies})
    dvec.index = dvec.index.astype(str)
    run = B.load_run(B.tag_of(B.Point(PANEL, 'held', 'real')))
    pn = R.stage_changes(run, 'pn_dendrite', 'on', odors)[ODOR]
    pn.index = pn.index.astype(str)
    gl = list(dict.fromkeys(list(inp.index.astype(str)) + list(vec_all.index.astype(str)) + list(dvec.index)))
    tab = pd.DataFrame(index=pd.Index(gl, name='glomerulus'))
    tab['in_model36'] = tab.index.isin(list(inp.index.astype(str)))
    tab['in_dendrite28'] = tab.index.isin(d28)
    tab['ir'] = tab.index.isin(ir)
    tab['model_input'] = inp.reindex(tab.index)
    for f in vec_all.columns:
        tab[f'orn_{f}'] = vec_all[f].reindex(tab.index)
    tab['orn_mean'] = vec_all.mean(axis=1).reindex(tab.index)
    tab['hallem'] = hal.reindex(tab.index)
    for f in dvec.columns:
        tab[f'dendrite_{f}'] = dvec[f].reindex(tab.index)
    tab['dendrite_mean'] = dvec.mean(axis=1).reindex(tab.index)
    tab['model_pn_family_held_real'] = pn.reindex(tab.index)
    _write_csv(tab.reset_index(), 'input_aa.csv')
    rows = [dict(source='model_input', **vector_summary(inp, None, ir))]
    for f in vec36.columns:
        rows.append(dict(source=f'orn_{f}', on='model36', **vector_summary(vec36[f], inp, ir)))
        rows.append(dict(source=f'orn_{f}', on='all_recorded', **vector_summary(vec_all[f], inp, ir)))
    rows.append(dict(source='orn_mean', on='model36', **vector_summary(vec36.mean(axis=1), inp, ir)))
    rows.append(dict(source='hallem', on='shared_with_model36', **vector_summary(hal36, inp, ir)))
    rows.append(dict(source='hallem', on='all_receptors', **vector_summary(hal, None, ir)))
    rows.append(dict(source='orn_mean', on='hallem_shared', **vector_summary(vec36.mean(axis=1).reindex(hal36.index),
                                                                             hal36, ir)))
    rows.append(dict(source='dendrite_mean', on='dendrite28',
                     **vector_summary(dvec.mean(axis=1).reindex([g for g in d28 if g in dvec.index]), inp, ir)))
    rows.append(dict(source='model_pn_family_held_real', on='dendrite28',
                     **vector_summary(pn.reindex([g for g in d28 if g in pn.index]),
                                      dvec.mean(axis=1), ir)))
    rows.append(dict(source='model_input', on='dendrite28', **vector_summary(inp.reindex(d28), None, ir)))
    summ = pd.DataFrame(rows)
    _write_csv(summ, 'input_summary.csv')
    provenance('input', lock, orn_map=info, model_pn_run=[dict(path=f['path'], sha256=f.get('sha256'))
                                                          for f in run.files],
               input_sha256=R.array_sha256(megamat_input()[odors].to_numpy()))
    return dict(table=tab, summary=summ, info=info)


# ---------------------------------------------------------------------------
# the classes

def aa_rows(df: pd.DataFrame, **sel) -> pd.DataFrame:
    """aa's rows of a table, selected by column values."""
    m = df.odor == ODOR
    for k, v in sel.items():
        m &= df[k].isin(v) if isinstance(v, (list, tuple, set)) else (df[k] == v)
    return df[m]


def part_classes() -> dict:
    """The classes (i), (ii) and (iv) read on step 0's tables as the lock fixes them, (iii) left to step 1 with its
    range recorded, and the outcomes they point to; classes.json and tables.md."""
    lock, _ = require_all()
    fl = pd.read_csv(out_dir() / 'flies_entries.csv')
    me = pd.read_csv(out_dir() / 'models_entries.csv')
    sr = pd.read_csv(out_dir() / 'flies_soma_readings.csv')
    ss = pd.read_csv(out_dir() / 'flies_sharing_summary.csv')
    ci = CLASSES['i']
    claws = [(r['cohort'], r['fly']) for r in lock['recordings']['claw_recordings']['recordings']]
    fc = pd.concat([aa_rows(fl, stage='kc_claw', cohort=c, recording=f, reading='lock') for c, f in claws])
    mc = aa_rows(me, group=ci['model_group'], condition=ci['condition'], stage=ci['stage'], reading=ci['reading'],
                 wiring=list(WIRINGS))
    keep = ['recording', 'cohort', 'rate', 'rate_rank', 'rate_third', 'concentration', 'concentration_rank',
            'concentration_third', 'cv', 'n_responders']
    c1 = class_i(fc[keep], mc[['wiring'] + keep[2:]])
    halves = {h: ss[(ss.stage == 'kc_soma') & (ss.cohort == 'wild_type') & (ss.set == f'half_{h}')]
              for h in lock['definitions']['soma_halves']}
    c2 = class_ii(halves)
    c2['all_flies'] = class_ii({'all': ss[(ss.stage == 'kc_soma') & (ss.cohort == 'wild_type') & (ss.set == 'all')]}
                               )['halves']['all']
    c2['beside_how_strongly'] = class_ii(halves, dict(CLASSES['ii'], part='sharing/how_strongly'))['halves']
    rd = aa_rows(sr, recording='mean')
    c4 = class_iv(rd)
    lk = aa_rows(fl, stage='kc_soma', cohort='wild_type', recording='mean', reading='lock')
    c3 = dict(name=CLASSES['iii']['name'], rule=CLASSES['iii']['rule'], holds=None, read_in='step 1',
              range=dict(lo=float(lk.rate_min.iloc[0]), hi=float(lk.rate_max.iloc[0])),
              model_now={w: float(v) for w, v in aa_rows(me, group='family_held', condition='on', stage='kc_soma',
                                                          reading='curve_noise', wiring=list(T_WIRINGS))
                         .set_index('wiring').rate.items()})
    out = dict(i=c1, ii=c2, iii=c3, iv=c4, outcome=outcomes_of(c1['holds'], None, c4['holds']),
               lock_sha256=lock['sha256'], time=_now(), host=_host())
    _write_json(out_dir() / 'classes.json', out)
    write_tables(out)
    provenance('classes', lock)
    print(json.dumps({k: (v.get('holds') if isinstance(v, dict) else v) for k, v in out.items()
                      if k in ('i', 'ii', 'iii', 'iv')}), out['outcome']['outcomes'], flush=True)
    return out


def _f(x, spec: str = '.2f') -> str:
    try:
        return format(float(x), spec) if np.isfinite(float(x)) else '–'
    except (TypeError, ValueError):
        return str(x)


def write_tables(classes: Mapping) -> None:
    """tables.md: aa along the stages (the flies, every recording, and each model group), the readings of the flies'
    somata, the same-fly step, the sharing rows, the anatomy, item 3 and item 4, and the classes."""
    from kcq import stage_readings as R
    d = out_dir()
    fl, me = pd.read_csv(d / 'flies_entries.csv'), pd.read_csv(d / 'models_entries.csv')
    sr, sfs = pd.read_csv(d / 'flies_soma_readings.csv'), pd.read_csv(d / 'flies_same_fly_step.csv')
    L = ['# E140: aa at the step from the claws to the somata (step 0, no simulation)', '',
         'Entries are multiples of the mean over the 17 megamat odors (rank 1 the lowest of 17); concentration is the '
         "share of the responders' summed response carried by the top tenth of them.", '',
         '## aa in the flies, every recording (the lock\'s reading)', '',
         '| stage | recording | rate (rank) | mean (rank) | concentration (rank) | CV | responders |',
         '|---|---|---:|---:|---:|---:|---:|']
    for r in aa_rows(fl, reading='lock').itertuples():
        L.append(f'| {r.stage} ({r.cohort}) | {r.recording} | {_f(r.rate)} ({_f(r.rate_rank, "g")}) | {_f(r.mean)} '
                 f'({_f(r.mean_rank, "g")}) | {_f(r.concentration, ".3f")} ({_f(r.concentration_rank, "g")}) | '
                 f'{_f(r.cv)} | {_f(r.n_responders, ".0f")} |')
    L += ['', '## aa in the models (the wirings\' mean views, E134\'s default readings)', '',
          '| group | condition | stage | reading | rate (rank) | mean (rank) | concentration (rank) |',
          '|---|---|---|---|---:|---:|---:|']
    for r in aa_rows(me, wiring='common').itertuples():
        L.append(f'| {r.group} | {r.condition} | {r.stage} | {r.reading} | {_f(r.rate)} ({_f(r.rate_rank, "g")}) | '
                 f'{_f(r.mean)} ({_f(r.mean_rank, "g")}) | {_f(r.concentration, ".3f")} '
                 f'({_f(r.concentration_rank, "g")}) |')
    L += ['', "## aa at the flies' somata under every reading (mean of the four flies)", '',
          '| reading | counts for (iv) | rate (rank) | range over flies | mean (rank) |', '|---|---|---:|---|---:|']
    for r in aa_rows(sr, recording='mean').itertuples():
        L.append(f'| {r.reading} | {"yes" if r.counts else "beside"} | {_f(r.rate)} ({_f(r.rate_rank, "g")}) | '
                 f'{_f(r.rate_min)} to {_f(r.rate_max)} | {_f(r.mean)} ({_f(r.mean_rank, "g")}) |')
    L += ['', '## The same-fly step and the four flies\' means', '',
          '| set | claws rate (rank) | somata rate (rank) | step, rate (rank) | step, mean |', '|---|---:|---:|---:|---:|']
    for r in aa_rows(sfs).itertuples():
        L.append(f'| {r.set} | {_f(r.claws_rate)} ({_f(r.claws_rate_rank, "g")}) | {_f(r.somata_rate)} '
                 f'({_f(r.somata_rate_rank, "g")}) | {_f(r.step_rate, "+.2f")} ({_f(r.step_rate_rank, "g")}) | '
                 f'{_f(r.step_mean, "+.2f")} |')
    if (d / 'models_breadth.csv').exists():
        fb, mb = pd.read_csv(d / 'flies_breadth.csv'), pd.read_csv(d / 'models_breadth.csv')
        L += ['', "## How many odors aa's somatic responders respond to (of 17)", '',
              '| side | responders | their mean breadth (rank) | all units\' mean breadth |', '|---|---:|---:|---:|']
        for r in aa_rows(fb, recording='mean').itertuples():
            L.append(f'| the flies ({r.cohort}) | {_f(r.n_responders, ".0f")} | {_f(r.breadth_mean)} '
                     f'({_f(r.breadth_mean_rank, "g")}) | {_f(r.all_units_breadth_mean)} |')
        for r in aa_rows(mb, wiring='mean').itertuples():
            L.append(f'| {r.group}, APL {r.condition} | {_f(r.n_responders, ".0f")} | {_f(r.breadth_mean)} '
                     f'({_f(r.breadth_mean_rank, "g")}) | {_f(r.all_units_breadth_mean)} |')
    if (d / 'follows.csv').exists():
        fo = pd.read_csv(d / 'follows.csv')
        L += ['', "## What each side's somata follow (aa left out of the fit; rate part)", '',
              '| side | x | r (16 odors) | aa: observed / predicted / residual (z) | aa\'s rank by leave-one-out residual |',
              '|---|---|---:|---|---:|']
        for r in fo[fo.y == 'kc_soma_rate'].itertuples():
            L.append(f'| {r.side} | {r.x} | {_f(r.r)} | {_f(r.y_aa)} / {_f(r.pred_aa)} / {_f(r.resid_aa, "+.2f")} '
                     f'({_f(r.z_aa, "+.1f")}) | {_f(r.loo_rank_aa, "g")} |')
    if (d / 'input_summary.csv').exists():
        su = pd.read_csv(d / 'input_summary.csv')
        L += ['', "## aa's input beside its measurements", '',
              '| source | on | glomeruli | sum | max (glomerulus) | top 3 | breadth | IR share | r with the model input |',
              '|---|---|---:|---:|---|---|---:|---:|---:|']
        for r in su.itertuples():
            L.append(f'| {r.source} | {getattr(r, "on", "") if isinstance(getattr(r, "on", ""), str) else ""} | '
                     f'{r.n_glomeruli} | {_f(r.sum, ".0f")} | {_f(r.max, ".0f")} ({r.argmax}) | {r.top} | '
                     f'{_f(r.breadth, ".1f")} | {_f(r.ir_share)} | {_f(getattr(r, "pearson_with_ref", np.nan))} |')
    L += ['', '## Classes', '', '```', json.dumps(classes, indent=1, default=_json_default), '```', '']
    (d / 'tables.md').write_text('\n'.join(L) + '\n')


# ---------------------------------------------------------------------------
# the figures

INK, MUTED = '#1f1f1e', '#a9a8a2'
COLORS = {'flies': INK, 'family_held': '#2a78d6', 'A1': '#eb6834', 'e135_imposed': '#1baf7a', 'aa': '#eb6834'}
STAGE_LABEL = {'orn': 'ORNs', 'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws',
               'kc_soma': 'KC somata'}


def part_figures() -> None:
    """fig_aa_stages.png (aa's rate entry along the stages, flies and models), fig_claws_rate_concentration.png (the
    rate entry against concentration at the claws, every odor, flies and the family held), fig_soma_follows.png (the
    somata against the input's sum in the model and against the claws in the flies), fig_aa_input.png (aa's input
    beside the ORN flies and Hallem)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from kcq import stage_readings as R
    lock, _ = require_all()
    d = out_dir()
    fl, me = pd.read_csv(d / 'flies_entries.csv'), pd.read_csv(d / 'models_entries.csv')
    st = ['pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma']
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})

    # aa along the stages
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, part in zip(axes, ('rate', 'mean')):
        for i, s in enumerate(st):
            x = aa_rows(fl, stage=s, cohort='wild_type', reading='lock')
            pts = x[x.recording != 'mean'][part]
            ax.plot(np.full(len(pts), i), pts, 'o', ms=4, color=MUTED, zorder=2)
        m = [float(aa_rows(fl, stage=s, cohort='wild_type', recording='mean', reading='lock')[part].iloc[0]) for s in st]
        ax.plot(range(len(st)), m, '-o', color=COLORS['flies'], lw=2, ms=6, label="the flies' mean", zorder=3)
        for g, lab in (('family_held', 'the family, held'), ('A1', "E136's A1"), ('e135_imposed', 'dendrites imposed')):
            y = []
            for s in st:
                x = aa_rows(me, group=g, wiring='common', condition='on', stage=s, reading=R.DEFAULT_READINGS[s].name)
                y.append(float(x[part].iloc[0]) if len(x) else np.nan)
            ax.plot(range(len(st)), y, '-s', color=COLORS[g], lw=2, ms=5, label=lab, zorder=3)
        ax.axhline(1.0, color=MUTED, lw=0.8, ls=':')
        ax.set_xticks(range(len(st)), [STAGE_LABEL[s] for s in st])
        ax.set_title(f"aa's {'response rate' if part == 'rate' else 'mean response'}, as a multiple of the panel's mean")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(d / 'fig_aa_stages.png', dpi=150)
    plt.close(fig)

    # the claws: rate entry against concentration, every odor
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=False)
    claws = [(r['cohort'], r['fly']) for r in lock['recordings']['claw_recordings']['recordings']]
    x = pd.concat([fl[(fl.stage == 'kc_claw') & (fl.cohort == c) & (fl.recording == f) & (fl.reading == 'lock')]
                   for c, f in claws])
    y = me[(me.group == 'family_held') & (me.stage == 'kc_claw') & (me.condition == 'on')
           & (me.reading == 'curve_noise') & me.wiring.isin(WIRINGS)]
    for ax, df, title in ((axes[0], x, 'KC claws, the five recordings'),
                          (axes[1], y, 'KC claws, the family held, eight wirings')):
        o = df[df.odor != ODOR]
        a = df[df.odor == ODOR]
        ax.plot(o.rate, o.concentration, 'o', ms=3, color=MUTED, label='the other 16 odors')
        ax.plot(a.rate, a.concentration, 'o', ms=6, color=COLORS['aa'], label='aa', zorder=3)
        ax.set_xlabel("response rate, as a multiple of the panel's mean")
        ax.set_ylabel("concentration (top tenth's share)")
        ax.set_title(title)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(d / 'fig_claws_rate_concentration.png', dpi=150)
    plt.close(fig)

    # what the somata follow
    if (d / 'follows_points.csv').exists():
        fp, fo = pd.read_csv(d / 'follows_points.csv'), pd.read_csv(d / 'follows.csv')
        panels = (('family_held', 'input_sum', 'the family held: KC somata against the summed input',
                   'summed input (spike deltas, Hz)'),
                  ('flies', 'stage_kc_claw', "the flies: KC somata against their claws", 'claws, response rate entry'))
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
        for ax, (side, xn, title, xl) in zip(axes, panels):
            p = fp[(fp.side == side) & (fp.x == xn) & (fp.y == 'kc_soma_rate')]
            f = fo[(fo.side == side) & (fo.x == xn) & (fo.y == 'kc_soma_rate')]
            o, a = p[p.odor != ODOR], p[p.odor == ODOR]
            ax.plot(o.x_value, o.y_value, 'o', ms=4, color=MUTED, label='the other 16 odors')
            ax.plot(a.x_value, a.y_value, 'o', ms=7, color=COLORS['aa'], label='aa', zorder=3)
            for r in o.itertuples():
                ax.annotate(r.odor, (r.x_value, r.y_value), fontsize=6, color='#6b6a64', xytext=(2, 2),
                            textcoords='offset points')
            if len(f) and np.isfinite(f.slope.iloc[0]):
                xs = np.linspace(np.nanmin(p.x_value), np.nanmax(p.x_value), 2)
                ax.plot(xs, f.intercept.iloc[0] + f.slope.iloc[0] * xs, '-', color=INK, lw=1,
                        label='line over the 16 (aa left out)')
            ax.set_title(title)
            ax.set_xlabel(xl)
            ax.set_ylabel("KC somata, response rate entry")
        axes[0].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(d / 'fig_soma_follows.png', dpi=150)
        plt.close(fig)

    # aa's input
    if (d / 'input_aa.csv').exists():
        t = pd.read_csv(d / 'input_aa.csv')
        t = t[t.in_model36].sort_values('model_input', ascending=False)
        fig, ax = plt.subplots(figsize=(10, 3.8))
        xs = np.arange(len(t))
        ax.bar(xs, t.model_input, color=COLORS['family_held'], width=0.7, label="the model's input")
        oc = [c for c in t.columns if c.startswith('orn_') and c != 'orn_mean']
        for c in oc:
            ax.plot(xs, t[c], '_', ms=8, mew=1.5, color=INK, alpha=0.6)
        ax.plot([], [], '_', ms=8, mew=1.5, color=INK, label='each ORN fly given aa (E85\'s map)')
        ax.plot(xs, t.hallem, 'D', ms=4, color=COLORS['aa'], label="Hallem's aa (10^-2), shared glomeruli")
        ax.set_xticks(xs, [f'{g}*' if ir else g for g, ir in zip(t.glomerulus, t.ir)], rotation=90, fontsize=7)
        ax.set_ylabel('spike deltas (Hz)')
        ax.set_title("aa's input on the model's 36 glomeruli (* IR-fed)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(d / 'fig_aa_input.png', dpi=150)
        plt.close(fig)
    provenance('figures', lock)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('lock', 'reduce', 'flies', 'models', 'follows', 'input', 'classes', 'figures', 'all'))
    a = ap.parse_args(argv)
    if a.part == 'all':
        for f in (part_flies, part_models, part_follows, part_input, part_classes, part_figures):
            f()
        return
    dict(lock=part_lock, reduce=part_reduce, flies=part_flies, models=part_models, follows=part_follows,
         input=part_input, classes=part_classes, figures=part_figures)[a.part]()


if __name__ == '__main__':
    main()
