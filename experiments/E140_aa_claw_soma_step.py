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
flies_entries.csv, flies_sharing.csv, flies_breadth.csv, flies_same_fly_step.csv, flies_soma_readings.csv;
models_entries.csv, models_sharing.csv, models_anatomy.csv, models_a1_stages.csv; follows.csv; input_aa.csv,
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

def concentration(values: np.ndarray, top: float = TOP_SHARE, min_n: int = MIN_RESPONDERS) -> float:
    """The share of the summed values carried by the largest `top` fraction of them (the Lorenz curve, largest first,
    at top x n, linear between units); NaN with fewer than `min_n` values or a non-positive sum."""
    raise NotImplementedError


def responders_cv(values: np.ndarray, min_n: int = MIN_RESPONDERS) -> float:
    """SD (ddof CV_DDOF) over mean of the responders' values; NaN with fewer than `min_n`."""
    raise NotImplementedError


def ranks(x: pd.Series) -> pd.Series:
    """Ascending ranks (1 the lowest, ties averaged) among the defined entries; NaN where undefined."""
    raise NotImplementedError


def third(rank: float, n_defined: int, k: int = THIRD) -> str:
    """'lowest' (rank <= k), 'highest' (rank > n_defined - k), else 'middle'; '' for an undefined rank."""
    raise NotImplementedError


def odor_stats(pop, odors: Sequence[str]) -> pd.DataFrame:
    """Per odor of one population (hal's `Population`): the strength parts' entries ('rate', 'mean'), the response rate
    and mean response, the responders' concentration and CV, the number of responders and units, and each column's
    rank and third."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the parts

def part_reduce() -> dict:
    """The four megamat soma flies' half-window reduction (soma_windows.parquet), after the lock."""
    raise NotImplementedError


def part_flies() -> dict:
    """Design item 1."""
    raise NotImplementedError


def e139_runs() -> Dict[str, list]:
    """E139's saved megamat runs, when they arrive ({setting: [SavedRun per wiring]}); {} until then."""
    raise NotImplementedError


def part_models() -> dict:
    """Design item 2."""
    raise NotImplementedError


def part_follows() -> dict:
    """Design item 3."""
    raise NotImplementedError


def part_input() -> dict:
    """Design item 4."""
    raise NotImplementedError


def part_classes() -> dict:
    """The classes (i) to (iv) and the outcome they point to; classes.json, tables.md."""
    raise NotImplementedError


def part_figures() -> None:
    """The figures."""
    raise NotImplementedError


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
