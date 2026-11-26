#!/usr/bin/env python3
"""E143, the megamat stage error through the flies' sampling: every split of the flies, an interval over flies on every
class since the stage fit, and how many flies a stage needs to tell models apart (no simulation; experiment_priorities.md,
E143; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 10:31; code claimed by hal 2026-11-26
10:55).

Nothing is run, fitted, searched or chosen: every saved megamat read of the archives with its four stages saved (the
settings, SETTINGS) is read once (`part_read`), and the flies are then split and resampled around those reads.

Definitions (the entry's, fixed in the lock before any statistic; experiments/E143_lock.json):
- the measure is E136's adopted one, strength by odor (MEASURE: kcq/population_fit.py's strength view, parts 'rate' and
  'mean', in one-fly units), with every view beside where a class's guard reads it (E136's rule 3 guard);
- a setting's value at a stage is the mean over its T wirings (P3, P4, M2, M3) of each wiring's error; its megamat stage
  error the mean over T of each wiring's root mean square of its four stages' errors;
- the model is read once: each run through E134's default reading (kcq.stage_readings' DEFAULT_READINGS) with the
  curve, noise and criterion of all the stage's flies, seed index the wiring's position in E136's WIRINGS
  (`read_setting`); a subset or a resample of the flies changes the reference (its mean view), never the model's
  reading or the one-fly unit (all the stage's flies'); the shortcut's check (`part_check`): on E136's own split it
  gives the family's and A1's held-half megamat values within CHECK_TOL of E136's score, else the splits are read
  E136's way on N_FALLBACK seeded joint splits (`fallback_splits`);
- a split: every way to divide a stage's flies as E136 did (`stage_splits`: the smaller half floor(n / 2) flies; 35 for
  the PN dendrites, 10 for the boutons, 3 each for the claws and the somata, a 2 / 2 split counted once); a joint split
  one split per stage: N_JOINT drawn without replacement from all of them (seed JOINT_SEED), and E136's own
  (`joint_splits`);
- consistency at a stage on a split: Spearman ρ over the settings (SETTINGS with role 'setting') between their values
  against the two halves (`consistency`);
- the interval over flies of a change (a setting minus the family, or minus another setting, on the same wirings): the
  5th to 95th percentile over N_BOOT resamples of the stage's or panel's flies with replacement (seed BOOT_SEED),
  megamat's four stages resampled together (`resamples`); the resampled flies' mean view is the reference, in all the
  flies' units;
- the calyx contrasts: CALYX_GROUPS, the nine settings at pn_inhadd 1000 and the five at the family's normalization;
- the fifth soma fly: the same-fly recording's somata (viki's population kc_soma / megamat / same_fly, 2023-04-11 fly 3,
  on E134's footing), in no fit, half or class; its claws (kc_claw / megamat / same_fly) a second acquisition of the
  claw fly 2023-04-11 fly03, read as that fly against itself.

Design (CLI `--part`):
- lock: the lock (files, flies, splits, seeds, settings, classes), written to results/ and experiments/E143_lock.json;
  commit the copy before any other part;
- read: every setting's runs read once (megamat on T; the held-out panels' KC somata on T where read), each wiring's
  parts against all flies saved to outputs/E143_megamat_fly_sampling/views/ (`part_read`; resumable, one file per
  setting and input);
- check: the read-once check (design item 1);
- splits: megamat on every split (item 2); calyx: megamat's calyx (item 3); intervals: the intervals over flies and the
  classes recomputed on each resample (item 4); count: how many flies (item 5); fifth: the fifth soma fly (item 6);
  no_aa: the megamat somata without aa's entry (item 7); classes: the entry's classes and the outcome; all: every part
  after read, in order.

Outputs: <shared>/results/E143_megamat_fly_sampling/: lock.json; check.json, check.csv, check_classes.csv (item 1);
splits.json, splits_consistency.csv, splits_stage_values.csv, splits_joint_consistency.csv, splits_beyond_bar.csv
(item 2); calyx.json, calyx_consistency.csv, calyx_contrasts.csv (item 3); intervals.json, intervals_changes.csv,
class_shares.csv (item 4); count.json, count_consistency.csv (item 5); fifth_fly.json, fifth_fly_settings.csv (item 6);
the files of items 2 to 4 with the suffix _no_aa (item 7); classes.json, results_tables.md, provenance.json and
provenance_<part>.json. outputs/E143_megamat_fly_sampling/: views/ (every setting's parts per input and the flies'),
read.log, read_summary_<host>.csv.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E143_megamat_fly_sampling --part lock    # then commit experiments/E143_lock.json
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E143_megamat_fly_sampling --part read [--jobs N]
    ./venv/bin/python -u -m experiments.E143_megamat_fly_sampling --part all
"""
from __future__ import annotations

import argparse
import contextlib
import itertools
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E143_megamat_fly_sampling'
SCRIPT = 'experiments/E143_megamat_fly_sampling.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E143_lock.json'
E134_LOCK_COPY = H.HERE / 'experiments' / 'E134_lock.json'
E134 = 'E134_population_baseline'
E135 = 'E135_antennal_lobe_dendrites'
E136 = 'E136_stage_population_fit'
E139 = 'E139_normalization_ladder'
E141 = 'E141_kckc_rise_forms'
PANEL = 'megamat'
MEASURE = 'strength'
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PN_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton')
CALYX_STAGES: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
# the held-out panels of the classes (E136's HOLDS_PANELS, as 'input/panel/condition') and the one beside (kmix)
HOLDS_PANELS: Tuple[str, ...] = ('validation2/validation2/on', 'mixtures/cmix/on', 'apl16/apl16/on', 'apl16/apl16/off')
BESIDE_PANELS: Tuple[str, ...] = ('mixtures/kmix/on',)
PANEL_INPUTS: Tuple[str, ...] = ('validation2', 'mixtures', 'apl16')


def side_key(class_key: str) -> str:
    """The side (flies) key 'kc_soma/<panel>/<condition>' of a class table's panel key 'input/panel/condition'."""
    _, pan, cond = class_key.split('/')
    return f'kc_soma/{pan}/{cond}'

# the read-once check (design item 1): E136's held-half megamat stage error on T, family and A1 (classes.json, measure
# strength), reproduced by the shortcut to within CHECK_TOL; otherwise E136's way on N_FALLBACK seeded joint splits
CHECK_ARMS: Tuple[str, ...] = ('family', 'E136:A1')
CHECK_TOL = 0.02
N_FALLBACK, FALLBACK_SEED = 20, 0
# the joint splits (design item 2), the resamples (item 4), the halves of item 5
N_JOINT, JOINT_SEED = 200, 0
N_BOOT, BOOT_SEED = 1000, 0
CI: Tuple[float, float] = (0.05, 0.95)
COUNT_SIZES: Tuple[int, ...] = (2, 3, 4, 5)
N_COUNT_DRAWS, COUNT_SEED = 200, 0
# the panels of item 5: (input/panel/condition, label). The entry names the control mixtures with 11 flies; E136 reads
# them on 4 (E134's fly_side_read of kc_soma / cmix / wild_type), 11 being the two mixture panels' KC flies together
# (the kiwi mixtures 7, the control mixtures 4). So the three panels with many flies are the APL16 somata with APL on
# (11) and off (8) and the kiwi mixtures (7, the mixture panel with the most flies; read as uploaded, a description of
# the panel, no verdict on a model); the count's and the panel's classes on the entry's three (the control mixtures,
# halves of 2 only) are given beside
COUNT_PANELS: Tuple[Tuple[str, str], ...] = (('apl16/apl16/on', 'APL16 somata, APL on'),
                                             ('apl16/apl16/off', 'APL16 somata, APL off'),
                                             ('mixtures/kmix/on', 'kiwi mixtures'),
                                             ('mixtures/cmix/on', 'control mixtures'),
                                             ('validation2/validation2/on', 'validation2'))
MANY_FLY_PANELS: Tuple[str, ...] = ('apl16/apl16/on', 'apl16/apl16/off', 'mixtures/kmix/on')
MANY_FLY_PANELS_AS_NAMED: Tuple[str, ...] = ('apl16/apl16/on', 'apl16/apl16/off', 'mixtures/cmix/on')
# the classes' thresholds (the entry's)
TELLS_MEDIAN, TELLS_MIN = 0.8, 0.5
STANDS_SHARE = 0.8
CLASS_SHARE = 0.8
NEAR_TIE_SHARE = 0.8      # "a near-tie class changes on more than a fifth of the resamples"
COUNT_LOW, COUNT_HIGH = 0.5, 0.8
# E136's classes' margin (IMPROVE_SE standard errors)
IMPROVE_SE = 2.0
# the megamat somata without this odor's entry (design item 7; E140's outcome (d))
WITHOUT_ODOR = 'aa'
# the realism flag that gets its interval (the entry's Realism): the shuffled matrix's response rate ρ with APL on, on
# P3, over the control recordings
FLAG_SETTING, FLAG_WIRING = 'E141:S_inhadd1000', 'P3'


@dataclass(frozen=True)
class Setting:
    """One saved megamat read (a model at fixed parameters on its wirings).

    Attributes:
        name: '<source>:<its name there>' ('family' for E134's held runs).
        source: the experiment whose runs these are ('E134', 'E136', 'E139', 'E141', 'E135').
        arm: its name in that experiment (E136's arm or registered variant, E135's arm).
        label: what it is, in words.
        role: 'setting' (in every consistency) or 'yardstick' (E135's permuted controls: scored, reported beside, in no
            consistency and no class).
        calyx: the calyx group it belongs to (CALYX_GROUPS' key) or ''.
        panels: the held-out inputs it was read on (PANEL_INPUTS), as saved.
    """
    name: str
    source: str
    arm: str
    label: str
    role: str = 'setting'
    calyx: str = ''
    panels: Tuple[str, ...] = PANEL_INPUTS


SETTINGS: Tuple[Setting, ...] = (
    Setting('family', 'E134', 'family', "the family held (E134's megamat_held runs; E136's reads and E90's part B on the "
            'held-out inputs)', calyx='family_norm'),
    Setting('E136:A0', 'E136', 'A0', 'A0, the family refitted', calyx='family_norm'),
    Setting('E136:A1', 'E136', 'A1', "A1, the antennal lobe's normalization searched", calyx='inhadd1000'),
    Setting('E136:A2', 'E136', 'A2', "A2, APL's inhibition of the claws as a gain", calyx='family_norm'),
    Setting('E136:A3', 'E136', 'A3', 'A3, the per-claw threshold searched', calyx='family_norm'),
    Setting('E136:A4', 'E136', 'A4', 'A4, without the inhibition among KCs', calyx='family_norm'),
    Setting('E139:F_inhadd100', 'E139', 'F_inhadd100', "F, the family's calyx factors, at pn_inhadd 100"),
    Setting('E139:F_inhadd314', 'E139', 'F_inhadd314', "F at pn_inhadd 314"),
    Setting('E139:F_inhadd1000', 'E139', 'F_inhadd1000', "F at pn_inhadd 1000", calyx='inhadd1000'),
    Setting('E139:F_inhadd3141', 'E139', 'F_inhadd3141', "F at pn_inhadd 3141"),
    Setting('E139:B_inhadd314', 'E139', 'B_inhadd314', 'B, APL->bouton x0.155, at pn_inhadd 314'),
    Setting('E139:B_inhadd1000', 'E139', 'B_inhadd1000', 'B at pn_inhadd 1000', calyx='inhadd1000'),
    Setting('E139:K_inhadd1000', 'E139', 'K_inhadd1000', 'K, KC-KC x3.24, at pn_inhadd 1000', calyx='inhadd1000'),
    Setting('E139:A0_inhadd314', 'E139', 'A0_inhadd314', "E139's refit, A0 at pn_inhadd 314"),
    Setting('E141:D1.8_inhadd1000', 'E141', 'D1.8_inhadd1000', 'D1.8, the matrix at x1.8', calyx='inhadd1000'),
    Setting('E141:D5.8_inhadd1000', 'E141', 'D5.8_inhadd1000', 'D5.8, the matrix at x5.83', calyx='inhadd1000'),
    Setting('E141:N_inhadd1000', 'E141', 'N_inhadd1000', 'N, the rise with APL on only', calyx='inhadd1000'),
    Setting('E141:G_inhadd1000', 'E141', 'G_inhadd1000', "G, APL's inhibition of the claws as a gain", calyx='inhadd1000'),
    Setting('E141:S_inhadd1000', 'E141', 'S_inhadd1000', 'the shuffled matrix (within type pairs)', calyx='inhadd1000'),
    Setting('E135:R1', 'E135', 'R1', "E135's R1, the PN stage's global parameters fitted", panels=('apl16',)),
    Setting('E135:rung', 'E135', 'rung', "E135's R3, the highest kept rung", panels=('apl16',)),
    Setting('E135:imposed', 'E135', 'imposed', "E135's recorded dendrites imposed", panels=('apl16',)),
    Setting('E135:imposed_dff', 'E135', 'imposed_dff', "E135's recorded dendrites imposed, table in dF/F", panels=()),
    Setting('E135:rung_glom_perm', 'E135', 'rung_glom_perm', "E135's R3, glomeruli permuted", role='yardstick',
            panels=()),
    Setting('E135:imposed_glom_perm', 'E135', 'imposed_glom_perm', "E135's imposed, glomeruli permuted",
            role='yardstick', panels=()),
    Setting('E135:imposed_odor_perm', 'E135', 'imposed_odor_perm', "E135's imposed, odors permuted", role='yardstick',
            panels=()),
    Setting('E135:imposed_glom_perm_dff', 'E135', 'imposed_glom_perm_dff', "E135's imposed in dF/F, glomeruli permuted",
            role='yardstick', panels=()),
    Setting('E135:imposed_odor_perm_dff', 'E135', 'imposed_odor_perm_dff', "E135's imposed in dF/F, odors permuted",
            role='yardstick', panels=()),
)
SETTING_BY_NAME: Dict[str, Setting] = {s.name: s for s in SETTINGS}
CALYX_GROUPS: Dict[str, Tuple[str, ...]] = {
    g: tuple(s.name for s in SETTINGS if s.calyx == g) for g in ('inhadd1000', 'family_norm')}
# the class tables whose changes get intervals and whose classes are recomputed on each resample (design item 4):
# (experiment, classes.json's path under results/, the settings' name prefix)
CLASS_TABLES: Tuple[Tuple[str, Tuple[str, ...], str], ...] = (
    (E136, ('measure_strength', 'classes.json'), 'E136:'),
    (E139, ('classes.json',), 'E139:'),
    (E141, ('classes.json',), 'E141:'),
)
# E141's arms whose reads stopped after their holds failed (no megamat reads on T): no class whatever the flies
STOPPED_ARMS: Tuple[str, ...] = ('E141:P_inhadd1000', 'E141:T_inhadd1000')
# the near-tie classes the entry names (outcome (d)), as (class table, its name there)
NEAR_TIES: Tuple[Tuple[str, str], ...] = ((E141, 'D1.8_inhadd1000'), (E139, 'K_inhadd1000'), (E141, 'K_inhadd1000'),
                                          (E141, 'N_inhadd1000'), (E141, 'G_inhadd1000'))
# beside (a description): the fit evaluations whose stages this host holds (E136's fits on the S wirings), read once
FIT_EVALS_BESIDE: Tuple[Tuple[str, str], ...] = ((E136, 'A3'),)
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')

# the classes and outcomes (the entry's, fixed beforehand)
CLASSES: Dict[str, str] = {
    'tells_apart': f'"A stage tells the saved models apart": consistency over the settings is {TELLS_MEDIAN:g} or more '
                   f'at the median split and {TELLS_MIN:g} or more on every split',
    'reads_calyx': '"Megamat reads the calyx": the same over the nine settings at pn_inhadd 1000, at the claws and at '
                   'the somata',
    'change_stands': f'"A change stands through the flies": its interval excludes zero; on megamat it is also beyond the '
                     f'bar against both halves on at least {STANDS_SHARE:.0%} of the joint splits',
    'class_stands': f'"A class stands through the flies": the scored class is read on at least {CLASS_SHARE:.0%} of the '
                    'resamples',
    'count': f'"The count\'s": on at least two of the three panels with many flies, halves of 2 give consistency under '
             f'{COUNT_LOW:g} and the largest halves read give {COUNT_HIGH:g} or more',
    'panel': f'"The panel\'s": halves of 2 already give {COUNT_HIGH:g} or more there (on at least two of the three)',
    'fifth_fly': 'the fifth fly is a description, never a class: it is one fly',
}
OUTCOMES: Dict[str, str] = {
    'a': 'megamat tells models apart at the PN stages only, and does not read the calyx',
    'b': 'all four stages tell models apart',
    'c': "A1's megamat fall or its \"holds\" does not stand through the flies",
    'd': 'near-tie classes change with the flies (at least one changes on more than a fifth of the resamples)',
    'e': "the count's or the panel's",
}
# the judgments the entry leaves open (the narrowest reading, closest to E136's code); each is in the lock
NARROW_READINGS: Tuple[str, ...] = (
    "the measure is E136's adopted one, strength by odor (E136's MEASURES['strength']); the PN dendrites read on it "
    "with E136's error column ('error', not the matrix view); E136's rule-3 guard reads every view (four at the PN "
    "dendrites) and is recomputed on each resample where a class reads it",
    'only the T wirings are read: every value is the mean over T, and every spread is the one scored (the family over '
    'its eight wirings, from the class tables); the fitted arms\' S wirings (their fits\' best evaluations) are not read',
    "E135's fitted (R1, R3) and imposed (from the trial-noise and the dF/F dendrites) arms are settings; their five "
    'permuted controls are yardsticks only (scored and reported, in no consistency and no class); consistency over the '
    "settings without E135's four is given beside, since their PN stages were set from the dendrite flies",
    "a split's halves: 'a' the smaller half (floor(n / 2) flies, E136's fit-half size), 'b' the rest; for a 2 / 2 split "
    "'a' the half with the stage's first fly, so each division is counted once; E136's own split is "
    "split_halves(n, 1, 0) as E136's split_names takes it",
    'the joint splits: JOINT_SEED draws N_JOINT distinct combinations (np.random.default_rng(JOINT_SEED).choice without '
    'replacement over the product of the four stages\' splits, in stage order), and E136\'s own split beside them',
    "consistency is Spearman ρ (scipy's, ties averaged) over the settings' values (the mean over T of each wiring's "
    'error) against the two halves, in all the stage\'s flies\' one-fly units; the median over a stage\'s splits '
    'decides "tells apart"',
    "beyond E136's bar against a half: the change (mean over T of the per-wiring changes) has the sign of the change as "
    "scored (the class table's, else against E136's held-out half) and exceeds IMPROVE_SE standard errors by E136's "
    "rule (compare: the larger of half the spread as scored and the SD of the per-wiring changes over T over 2)",
    'the resamples: each stage\'s and each panel\'s flies drawn with replacement, as many as the stage or panel has, '
    'every stage and panel drawn independently on the same resample index (one generator, BOOT_SEED, stages in '
    'LOSS_STAGES order then the panels in HOLDS_PANELS + BESIDE_PANELS order); the reference is the drawn flies\' mean '
    'view, in all the flies\' one-fly units',
    "a class recomputed on a resample: E136's classify ladder on the resampled values with every comparison's spread "
    "and bar (IMPROVE_SE standard errors) as scored in the class table, the permuted control's fall as scored, and the "
    "class table's reasons that do not depend on the flies (a bound, a realism band, a failed hold, reads stopped) "
    "kept; E141's added classes (stands in, does more than A1) recomputed the same way against A1",
    'the class against all the flies, unresampled, is given beside each share',
    "the changes minus another setting are E141's comparisons with A1 (against_a1: its megamat stage error and the "
    "held-out panels); the class tables' names map to the settings: 'A1' is E136's A1, E139's names E139's reads "
    "(E141's B_inhadd1000 and K_inhadd1000 among them), E141's names its reads",
    "fit evaluations beside (FIT_EVALS_BESIDE): every saved evaluation of E136's fits whose stages this host holds "
    "(hal: A3's), read once on S as the settings are on T, their consistency on every split given beside as a "
    'description (E136\'s own read 483 evaluations of all fits on one split); they enter no class',
    "item 5: the control mixtures have 4 flies as E136 reads them (the entry's 11 is the two mixture panels' KC flies "
    'together), so the third panel with many flies is the kiwi mixtures (7 flies; MANY_FLY_PANELS), and the classes on '
    "the entry's three (MANY_FLY_PANELS_AS_NAMED, the control mixtures at halves of 2 only) are given beside",
    "item 5 beside: megamat's four stages read the same way at the half sizes they allow (halves of 2 at the claws "
    'and somata are their three splits)',
    'item 5: disjoint halves of k flies drawn N_COUNT_DRAWS times (seed COUNT_SEED, each draw a permutation of the '
    "panel's flies, its first k one half and its next k the other); consistency over the settings read on that panel; "
    '"halves of 2" and "the largest halves" read at the median over draws',
    "the fifth fly: each setting's somata strength error against the same-fly recording's somata alone (its parts the "
    "reference's mean), in the four soma flies' one-fly units, beside each of the four against the other three and "
    "each setting against the four; the claw fly 2023-04-11 fly03 against its second acquisition likewise",
    "item 7: E136's stage_without_odor rule (E138's odor_entries: aa's term taken out of the squared error, the number "
    'of odors kept, nothing renormalized) applied to the somata parts read once; items 2 to 4 on the somata so read',
    "the realism flag's interval: the response rate ρ (Spearman over the 16 odors) of the shuffled matrix's APL-on read "
    "on P3 against the per-odor response fractions of resampled control recordings (each recording's KCs above the "
    "bands' noise criterion), N_BOOT draws with replacement, seed BOOT_SEED",
)
HOLDS: Tuple[str, ...] = (
    'no simulation, no model run, no olfsysm change, nothing fitted, searched or chosen',
    'no held-out panel enters a fit; nothing is chosen from E138\'s readings of validation2',
    "Remy's processed data only, without the bouton fly 2022-10-21 fly02 and without her older megamat flies",
    "Yang's processed products only", 'no connectome table',
)


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E143_megamat_fly_sampling/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    """outputs/E143_megamat_fly_sampling/[<part>/] (created)."""
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _e140():
    from experiments import E140_aa_claw_soma_step as E140
    return E140


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_e140()._json_default) + '\n')
    os.replace(tmp, path)


def results_file(experiment: str, *parts: str) -> Path:
    """<shared>/results/<experiment>/<parts...> (not created)."""
    return H.shared_dir() / 'results' / experiment / Path(*parts)


# ---------------------------------------------------------------------------
# the splits

def stage_splits(n: int) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """Every way to divide n flies as E136 did, as (half a, half b) index tuples: half a the smaller half (floor(n / 2)
    flies), for an even n the half with fly 0, so each division is counted once (n = 7: 35, 5: 10, 4: 3)."""
    if n < 2:
        raise ValueError(f'{n=}')
    k = n // 2
    out = []
    for a in itertools.combinations(range(n), k):
        if 2 * k == n and 0 not in a:
            continue
        out.append((tuple(a), tuple(i for i in range(n) if i not in a)))
    return out


def e136_split(n: int) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """E136's own split of n flies (hal's `split_halves(n, 1, 0)`: its fit half, then its held-out half)."""
    from kcq import population_fit as P
    a, b = P.split_halves(n, 1, 0)[0]
    return tuple(int(i) for i in a), tuple(int(i) for i in b)


def split_index(n: int, split: Tuple[Sequence[int], Sequence[int]]) -> int:
    """The index in `stage_splits(n)` of a division of n flies (either half may come first)."""
    a = tuple(sorted(int(i) for i in split[0]))
    b = tuple(sorted(int(i) for i in split[1]))
    for i, (x, y) in enumerate(stage_splits(n)):
        if (x, y) in ((a, b), (b, a)):
            return i
    raise ValueError(f'{split} is not a split of {n} flies')


def joint_splits(sizes: Mapping[str, int], n_joint: int = N_JOINT, seed: int = JOINT_SEED) -> pd.DataFrame:
    """The joint splits: one split per stage (an index into `stage_splits`), `n_joint` drawn without replacement from
    their product (stage order LOSS_STAGES), then E136's own split ('e136'). Columns: joint, <stage> (the index), and
    whether the row is E136's own."""
    stages = [s for s in LOSS_STAGES if s in sizes]
    counts = [len(stage_splits(int(sizes[s]))) for s in stages]
    total = int(np.prod(counts))
    pick = np.random.default_rng(seed).choice(total, size=min(int(n_joint), total), replace=False)
    idx = np.array(np.unravel_index(pick, counts)).T
    rows = [dict(joint=int(j), **{s: int(v) for s, v in zip(stages, r)}, e136=False) for j, r in enumerate(idx)]
    rows.append(dict(joint=len(rows), **{s: split_index(int(sizes[s]), e136_split(int(sizes[s]))) for s in stages},
                     e136=True))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the lock

def lock_inputs() -> dict:
    """The inputs the lock records with sha256: the locks and data files every statistic rests on ('data', checked when
    read) and the class tables and scores the parts reproduce or reread ('results', recorded)."""
    E140 = _e140()
    roi = H.shared_dir() / 'data' / 'remy_per_roi'
    data = dict(e134_lock=E140.file_entry(E134_LOCK_COPY),
                e135_lock=E140.file_entry(results_file(E135, 'lock.json')),
                same_fly_trials=E140.file_entry(roi / 'same_fly_trials.npz'),
                same_fly_checks=E140.file_entry(roi / 'same_fly_checks.json'))
    res = {}
    for k, (exp, parts) in dict(
            e136_classes=(E136, ('measure_strength', 'classes.json')),
            e136_score_megamat=(E136, ('measure_strength', 'score_megamat.csv')),
            e136_score_panels=(E136, ('measure_strength', 'score_panels.csv')),
            e136_consistency=(E136, ('measure_strength', 'consistency_strength.csv')),
            e139_classes=(E139, ('classes.json',)), e139_score_megamat=(E139, ('score_megamat.csv',)),
            e139_score_panels=(E139, ('score_panels.csv',)), e139_without_aa=(E139, ('score_without_aa.csv',)),
            e141_classes=(E141, ('classes.json',)), e141_score_megamat=(E141, ('score_megamat.csv',)),
            e141_score_panels=(E141, ('score_panels.csv',)), e141_score_apl16=(E141, ('score_apl16.csv',)),
            e141_without_aa=(E141, ('score_without_aa.csv',)),
            e135_score_scores=(E135, ('score_scores.csv',)), e135_verdict=(E135, ('score_verdict.json',)),
            e134_score_scores=(E134, ('score_scores.csv',))).items():
        p = results_file(exp, *parts)
        res[k] = E140.file_entry(p) if p.exists() else dict(path=str(p), missing=True)
    return dict(data=data, results=res)


def recordings() -> dict:
    """The recordings read (names only, no value): per megamat stage the flies, criterion and E136's halves; the
    held-out panels' flies as E134 reads them; the same-fly populations; the exclusions."""
    from data_prep import stage_population_data as D
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    out: Dict[str, dict] = {}
    for stage in LOSS_STAGES:
        p = D.load_population(stage, PANEL, 'wild_type')
        n = len(p.flies)
        out[f'{stage}/{PANEL}/on'] = dict(stage=stage, panel=PANEL, cohort='wild_type', flies=list(p.flies),
                                          criterion=float(p.criterion), units=p.units, odors=list(p.odors),
                                          e136_halves=F.split_names(p, F.HALF_SEED), n_splits=len(stage_splits(n)),
                                          e136_split_index=split_index(n, e136_split(n)))
    if any(f.startswith('2022-10-21__fly02') for f in out[f'pn_bouton/{PANEL}/on']['flies']):
        raise RuntimeError('the excluded bouton fly is in the bouton population')
    for key in HOLDS_PANELS + BESIDE_PANELS:
        _, panel, cond = key.split('/')
        fl = _fly_side('kc_soma', panel, cond)
        fp, have = B.fly_side_read(fl)
        out[key] = dict(stage='kc_soma', panel=panel, condition=cond, cohort=fl.key[2], flies=[f.name for f in fp],
                        criterion=float(fl.criterion), units=fl.units, odors=list(have))
    for stage in ('kc_claw', 'kc_soma'):
        p = D.load_population(stage, PANEL, 'same_fly')
        out[f'{stage}/{PANEL}/same_fly'] = dict(stage=stage, panel=PANEL, cohort='same_fly', flies=list(p.flies),
                                                criterion=float(p.criterion), units=p.units,
                                                n_units={f: int(p.tables[f].shape[0]) for f in p.flies})
    out['excluded'] = dict(
        flies=[dict(stage=x['stage'], fly=x['fly'], reason=x['reason']) for x in D.EXCLUDED_FLIES],
        older_megamat="Remy's older megamat flies (pyolfsysm's data/from_remy/2024-11-12): never read",
        same_fly='the same-fly claw and soma recording (2023-04-11 fly 3, megamat1_redo): read only in design item 6')
    return out


def _fly_side(stage: str, panel: str, cond: str):
    """The flies a model condition is read against (stage_readings' `fly_side`; RuntimeError if none)."""
    from kcq import stage_readings as R
    fl = R.fly_side(stage, panel, cond)
    if fl is None:
        raise RuntimeError(f'no recorded {stage} on {panel} ({cond})')
    return fl


def lock_contents() -> dict:
    """The lock: inputs, recordings, splits, seeds, settings, classes, outcomes, narrow readings and holds. No statistic
    is read."""
    from experiments import E136_stage_population_fit as F
    from kcq import stage_readings as R
    recs = recordings()
    sizes = {s: len(recs[f'{s}/{PANEL}/on']['flies']) for s in LOSS_STAGES}
    js = joint_splits(sizes)
    splits = dict(rule="every way to divide a stage's flies as E136 did (stage_splits)",
                  per_stage={s: [dict(a=list(a), b=list(b)) for a, b in stage_splits(n)] for s, n in sizes.items()},
                  e136={s: dict(zip(('fit', 'held'), map(list, e136_split(n)))) for s, n in sizes.items()},
                  joint=dict(n=N_JOINT, seed=JOINT_SEED, rows=js.to_dict('records')))
    seeds = dict(joint=JOINT_SEED, boot=BOOT_SEED, n_boot=N_BOOT, ci=list(CI), count=COUNT_SEED,
                 n_count_draws=N_COUNT_DRAWS, count_sizes=list(COUNT_SIZES), fallback=FALLBACK_SEED,
                 n_fallback=N_FALLBACK)
    settings = [dict(name=s.name, source=s.source, arm=s.arm, label=s.label, role=s.role, calyx=s.calyx,
                     panels=list(s.panels)) for s in SETTINGS]
    definitions = dict(
        measure=f"E136's MEASURES[{MEASURE!r}] (strength by odor: kcq/population_fit.py's parts strength/rate and "
                "strength/mean, each odor's entry relative to the panel's mean, the view's error the RMS of the two "
                "parts' errors in one-fly units); the guard every view (E136's error_column: the PN dendrites on four)",
        value='the mean over the T wirings of each wiring\'s error; the megamat stage error per wiring the RMS of its four '
              'stages\' errors',
        reading="each run once: stage_readings.read_run with E134's default readings against all the stage's (or "
                "panel's) flies (E136's Flies(None) for megamat, a FlyCache for the held-out panels), seed index "
                "E136's seed_index(wiring); the parts reference_views against all the flies' reference",
        reference="a subset or resample: the weighted mean of the flies' parts (NaN-aware), in all the flies' one-fly "
                  "units (kcq/population_fit.py's subset_reference(keep_units=True))",
        check=dict(arms=list(CHECK_ARMS), tol=CHECK_TOL, against="E136's classes.json (measure strength): the "
                   "megamat 'arm' value on T against the held-out half", fallback=dict(n=N_FALLBACK, seed=FALLBACK_SEED)),
        consistency="Spearman ρ over the settings (role 'setting') of their values against half a and half b",
        interval=f'{CI[0]:.0%} to {CI[1]:.0%} percentiles over {N_BOOT} resamples (seed {BOOT_SEED})',
        bar=f"IMPROVE_SE = {IMPROVE_SE:g} standard errors by E136's compare (se the larger of half the spread and the "
            "SD of the per-wiring changes over T over sqrt(4)); spreads as scored in the class tables",
        calyx_groups={k: list(v) for k, v in CALYX_GROUPS.items()},
        class_tables=[dict(experiment=e, path='/'.join(p), prefix=pre) for e, p, pre in CLASS_TABLES],
        stopped=list(STOPPED_ARMS), near_ties=list(NEAR_TIES),
        count=dict(panels=[k for k, _ in COUNT_PANELS], many=list(MANY_FLY_PANELS),
                   many_as_named=list(MANY_FLY_PANELS_AS_NAMED), sizes=list(COUNT_SIZES),
                   draws=N_COUNT_DRAWS, seed=COUNT_SEED, statistic='the median over draws'),
        without_odor=WITHOUT_ODOR, flag=dict(setting=FLAG_SETTING, wiring=FLAG_WIRING, statistic='response rate ρ, '
                                             'APL on, against the control recordings in baseline-noise units'),
        thresholds=dict(tells_median=TELLS_MEDIAN, tells_min=TELLS_MIN, stands_share=STANDS_SHARE,
                        class_share=CLASS_SHARE, near_tie_share=NEAR_TIE_SHARE, count_low=COUNT_LOW,
                        count_high=COUNT_HIGH))
    return dict(experiment=EXPERIMENT, script=SCRIPT, part='lock', written=_now(), host=_host(),
                commits=H.code_commits(), panel=PANEL, measure=MEASURE, stages=list(LOSS_STAGES),
                t_wirings=list(T_WIRINGS), e136_wirings=list(F.WIRINGS), readings=R.reading_definitions(),
                inputs=lock_inputs(), recordings=recs, splits=splits, seeds=seeds, settings=settings,
                definitions=definitions, classes=CLASSES, outcomes=OUTCOMES, narrow_readings=list(NARROW_READINGS),
                holds=list(HOLDS),
                note='written before any statistic of E143 was computed; every part but lock refuses until this file '
                     'is committed (require_lock)')


def part_lock() -> dict:
    """results/.../lock.json and its committed copy experiments/E143_lock.json (commit it before any other part);
    refuses to overwrite a committed lock."""
    E140 = _e140()
    if LOCK_COPY.exists() and E140._committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is committed already; the lock is fixed')
    lock = lock_contents()
    text = json.dumps(lock, indent=1, default=E140._json_default) + '\n'
    (out_dir() / LOCK_NAME).write_text(text)
    LOCK_COPY.write_text(text)
    print(f'wrote {out_dir() / LOCK_NAME} and {LOCK_COPY} (commit the copy before any other part)', flush=True)
    return lock


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (tracked, no uncommitted change, the same bytes as the shared copy); PermissionError
    otherwise."""
    E140 = _e140()
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it)')
    if not E140._committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = (out_dir() / LOCK_NAME) if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = E140._sha256(copy)
    return lock


def require_all() -> Tuple[dict, dict]:
    """(this experiment's lock, E134's lock), both required before any statistic; the lock's data files checked."""
    from kcq import stage_readings as R
    lock = require_lock()
    for e in lock['inputs']['data'].values():
        _e140()._check_file(e)
    return lock, R.require_lock()


# ---------------------------------------------------------------------------
# the read (each run once against all flies)

@contextlib.contextmanager
def e136():
    """E136's module on the measure MEASURE with E139's and E141's variants registered (E139's refit by its fit
    loader), for the block."""
    from experiments import E139_normalization_ladder as E139m
    from experiments import E141_kckc_rise_forms as E141m
    with E139m.e136() as F:
        for st in SETTINGS:
            if st.source == 'E139':
                if E139m.fit_rung(st.arm) is not None:
                    F.register_variant(st.arm, E139m.load_fit_run)
                else:
                    E139m.register_variants(F, [st.arm])
            elif st.source == 'E141':
                E141m.register_variants(F, [st.arm])
        yield F


def setting_runs(setting: Setting, panel: str = PANEL, wirings: Sequence[str] = T_WIRINGS) -> Dict[str, object]:
    """A setting's saved runs on one input, per wiring of `wirings` (missing wirings left out): E136's `load_read`
    (E134's held runs and E90's part B for the family, the reads of E136, E139 and E141), E135's step-2 runs."""
    from kcq import stage_readings as R
    out = {}
    if setting.source == 'E135':
        from experiments import E135_antennal_lobe_dendrites as E135m
        for w in wirings:
            try:
                out[w] = R.load_step3_run(E135m.tag_of(E135m.Arm(panel, setting.arm), w), dirs=E135m.step2_dirs())
            except FileNotFoundError:
                pass
        return out
    with e136() as F:
        for w in wirings:
            r = F.load_read(F.ReadPoint(setting.arm, panel, w, 'held' if panel == PANEL else None))
            if r is not None:
                out[w] = r
    return out


def fit_eval_runs(experiment: str, arm: str) -> Dict[str, Dict[str, object]]:
    """The saved evaluations of an E136 fit whose stages this host holds, per evaluation its S runs (complete only)."""
    from kcq import stage_readings as R
    out = {}
    with e136() as F:
        if experiment != E136:
            raise ValueError(experiment)
        for i in F.history_evals(arm):
            try:
                out[f'e{i}'] = {w: R.load_step3_run(f'{arm}_e{i}_{w}', dirs=F.fit_dirs(arm)) for w in S_WIRINGS}
            except FileNotFoundError:
                continue
    return out


def views_dir() -> Path:
    """outputs/E143_megamat_fly_sampling/views/ (created)."""
    return local_dir('views')


def _safe(name: str) -> str:
    return name.replace(':', '__').replace('/', '_')


def views_path(name: str, panel: str) -> Path:
    """The saved parts of one setting (or fit evaluation) on one input."""
    return views_dir() / f'{_safe(name)}__{panel}.parquet'


def flies_path(key: str) -> Path:
    """The saved parts and one-fly units of the flies of one stage or panel ('stage/panel/condition')."""
    return views_dir() / f'flies__{_safe(key)}.parquet'


_REFS: Dict[Tuple[str, Tuple[str, ...]], object] = {}


def _reference(key: str, odors: Sequence[str]):
    """All the flies' reference of a stage or panel on the model's odors (those the flies have), every part: megamat
    stages as E136's Flies(None) reads them, held-out panels as E136's score_panels (E134's fly_side_read); cached, and
    saved once to flies_path(key)."""
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    ck = (key, tuple(odors))
    if ck not in _REFS:
        stage, panel, cond = key.split('/')
        fl = _fly_side(stage, panel, cond)
        fp, have = B.fly_side_read(fl)
        od = [o for o in odors if o in set(have)]
        ref = P.fly_reference(fp, od, B.breadth_edges(len(od)))
        _REFS[ck] = ref
        path = flies_path(key)
        if not path.exists():
            rows = [pd.DataFrame(dict(key=key, part=q, fly=f, entry=list(df.columns),
                                      value=df.loc[f].to_numpy(dtype=float), unit=float(ref.one_fly[q])))
                    for q, df in ref.parts.items() for f in df.index]
            tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
            pd.concat(rows, ignore_index=True).to_parquet(tmp, index=False)
            os.replace(tmp, path)
    return _REFS[ck]


def views_rows(pops: Mapping[str, object], key: str, name: str, panel: str, wiring: str) -> List[pd.DataFrame]:
    """One wiring's populations (`read_run`'s) at a key read against all the flies' reference: long rows."""
    from kcq import population_fit as P
    pop = pops[key]
    ref = _reference(key, list(pop.values.columns))
    v = P.reference_views(pop, ref)
    return [pd.DataFrame(dict(setting=name, input=panel, key=key, wiring=wiring, part=q, entry=list(x.index),
                              value=x.to_numpy(dtype=float))) for q, x in v.items() if q in ref.parts]


def read_runs(name: str, panel: str, runs: Mapping[str, object]) -> pd.DataFrame:
    """Runs of one setting (or fit evaluation) on one input read once against all flies (E134's default readings, E136's
    seed index): megamat at the four stages, a held-out input at its KC somata on every panel and condition with
    flies."""
    from experiments import E136_stage_population_fit as F
    from kcq import stage_readings as R
    lock = R.require_lock()
    rows = []
    for w, run in runs.items():
        if panel == PANEL:
            flies = F.Flies(None, LOSS_STAGES, views=None)
            rd = R.read_run(run, stages=list(LOSS_STAGES), lock=lock, flies=flies.cache, seed_index=F.seed_index(w))
            keys = [f'{st}/{PANEL}/on' for st in LOSS_STAGES]
        else:
            rd = R.read_run(run, stages=['kc_soma'], lock=lock, flies=R.FlyCache(), seed_index=F.seed_index(w))
            want = {side_key(k) for k in HOLDS_PANELS + BESIDE_PANELS}
            keys = [k for k in ('/'.join(t) for t in rd.populations) if k in want]
        pops = {'/'.join(k): v for k, v in rd.populations.items()}
        for key in keys:
            rows += views_rows(pops, key, name, panel, w)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def read_setting(setting: Setting, panel: str = PANEL) -> Optional[pd.DataFrame]:
    """One setting's T runs on one input read once (`read_runs`), each wiring's parts (every view) against all the flies
    of each stage or panel read: long table (setting, input, key 'stage/panel/condition', wiring, part, entry, value),
    saved to views_path. None if no run is saved."""
    runs = setting_runs(setting, panel)
    if not runs:
        return None
    df = read_runs(setting.name, panel, runs)
    df['n_wirings'] = len(runs)
    path = views_path(setting.name, panel)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)
    return df


def read_fit_evals(experiment: str, arm: str) -> Optional[pd.DataFrame]:
    """The fit evaluations beside (FIT_EVALS_BESIDE) read once on S, setting '<experiment>:<arm>_<eval>'; saved."""
    path = views_path(f'fitevals_{experiment}_{arm}', PANEL)
    if path.exists():
        return pd.read_parquet(path)
    evs = fit_eval_runs(experiment, arm)
    if not evs:
        return None
    df = pd.concat([read_runs(f'fit:{arm}_{e}', PANEL, runs) for e, runs in evs.items()], ignore_index=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)
    return df


def read_jobs(settings: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """(setting, input) pairs to read, megamat first, then the held-out inputs."""
    want = [s for s in SETTINGS if settings is None or s.name in settings]
    return [(s.name, PANEL) for s in want] + [(s.name, p) for s in want for p in s.panels]


def _read_job(job: Tuple[str, str]) -> dict:
    name, panel = job
    path = views_path(name, panel)
    if path.exists():
        return dict(setting=name, input=panel, status='saved before', seconds=0.0)
    t0 = time.time()
    try:
        df = read_setting(SETTING_BY_NAME[name], panel)
    except Exception as e:      # reported in the summary, the other jobs go on
        return dict(setting=name, input=panel, status=f'failed: {type(e).__name__}: {e}', seconds=time.time() - t0)
    st = 'no runs' if df is None else f"read {int(df['n_wirings'].iloc[0])} wirings"
    print(f'E143 read {name} {panel}: {st}, {time.time() - t0:.0f} s', flush=True)
    return dict(setting=name, input=panel, status=st, seconds=round(time.time() - t0, 1))


def part_read(jobs: int = 1, settings: Optional[Sequence[str]] = None, fit_evals: bool = True) -> pd.DataFrame:
    """Every setting on every input it was read on (`read_setting`), skipping what is saved, in `jobs` processes, then
    the fit evaluations beside; a summary per file (read_summary_<host>.csv in the outputs)."""
    todo = read_jobs(settings)
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(jobs) as pool:
            rows = list(pool.imap_unordered(_read_job, todo))
    else:
        rows = [_read_job(j) for j in todo]
    if fit_evals and settings is None:
        for exp, arm in FIT_EVALS_BESIDE:
            t0 = time.time()
            df = read_fit_evals(exp, arm)
            rows.append(dict(setting=f'fitevals_{exp}_{arm}', input=PANEL,
                             status='no evaluations held' if df is None else f'{df.setting.nunique()} evaluations',
                             seconds=round(time.time() - t0, 1)))
    out = pd.DataFrame(rows)
    out.to_csv(local_dir() / f'read_summary_{_host()}.csv', index=False)
    from kcq import stage_readings as R
    provenance('read', require_lock(), e134_lock_sha256=R.require_lock().get('sha256'), jobs=int(jobs),
               threads=os.environ.get('OMP_NUM_THREADS'), views=str(views_dir()),
               files=sorted(f.name for f in views_dir().glob('*.parquet')),
               summary=out.to_dict('records'))
    return out


# ---------------------------------------------------------------------------
# the flies' side and the scoring of many references at once

@dataclass
class Side:
    """One stage or panel's flies and the settings' parts read against them.

    Attributes:
        key: 'stage/panel/condition' (megamat stages 'kc_soma/megamat/on', ...).
        flies: the flies' names.
        fly_parts: per part, flies x entries.
        units: per part, all the flies' one-fly unit.
        model: per part, (setting, wiring) rows x entries.
        rows: the model rows' (setting, wiring) index.
    """
    key: str
    flies: List[str]
    fly_parts: Dict[str, np.ndarray]
    units: Dict[str, float]
    model: Dict[str, np.ndarray]
    rows: pd.MultiIndex


def side_from_tables(key: str, flies: pd.DataFrame, model: pd.DataFrame) -> Side:
    """A side from long tables: `flies` (part, fly, entry, value, unit), `model` (setting, wiring, part, entry, value);
    each part's entries in the flies' order, the model's rows (setting, wiring) in the order first met."""
    from kcq import population_fit as P
    fparts, units, mparts = {}, {}, {}
    names = list(dict.fromkeys(flies['fly'].astype(str)))
    rows = pd.MultiIndex.from_frame(model[['setting', 'wiring']].drop_duplicates()) if len(model) else \
        pd.MultiIndex.from_tuples([], names=['setting', 'wiring'])
    for q in [q for q in P.ALL_PARTS if q in set(flies['part'])]:
        f = flies[flies['part'] == q]
        entries = list(dict.fromkeys(f['entry'].astype(str)))
        fparts[q] = f.pivot(index='fly', columns='entry', values='value').reindex(index=names, columns=entries
                                                                                 ).to_numpy(dtype=float)
        units[q] = float(f['unit'].iloc[0])
        m = model[model['part'] == q]
        M = m.pivot_table(index=['setting', 'wiring'], columns='entry', values='value', aggfunc='first', dropna=False)
        mparts[q] = M.reindex(index=rows, columns=entries).to_numpy(dtype=float)
    return Side(key=key, flies=names, fly_parts=fparts, units=units, model=mparts, rows=rows)


def load_side(key: str, names: Optional[Sequence[str]] = None, fit_evals: bool = False) -> Side:
    """A stage or panel's side from the saved reads (`read_setting`; the settings of `names`, default all saved, with
    the fit evaluations beside if `fit_evals`) and the flies' parts saved with them (`flies_path`)."""
    fl = pd.read_parquet(flies_path(key))
    files = sorted(views_dir().glob('*.parquet'))
    files = [f for f in files if not f.name.startswith('flies__') and (fit_evals or not f.name.startswith('fitevals_'))]
    dfs = []
    for f in files:
        d = pd.read_parquet(f, filters=[('key', '==', key)])
        if len(d) and (names is None or d['setting'].iloc[0] in set(names) or f.name.startswith('fitevals_')):
            dfs.append(d)
    model = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=['setting', 'wiring', 'part', 'entry',
                                                                                 'value'])
    if names is not None:
        order = {n: i for i, n in enumerate(names)}
        model = model.assign(_o=model['setting'].map(lambda x: order.get(x, len(order)))).sort_values(
            ['_o', 'wiring'], kind='stable').drop(columns='_o')
    return side_from_tables(key, fl, model)


def half_weights(n: int, halves: Sequence[Sequence[int]]) -> np.ndarray:
    """References x flies weights of halves (each a set of fly indices; 1 / its size on its flies)."""
    W = np.zeros((len(halves), n))
    for i, h in enumerate(halves):
        h = list(h)
        W[i, h] = 1.0 / len(h)
    return W


def boot_weights(n: int, draws: np.ndarray) -> np.ndarray:
    """References x flies weights of resamples (`draws`: resamples x n fly indices drawn with replacement)."""
    W = np.zeros((len(draws), n))
    for i, d in enumerate(draws):
        np.add.at(W[i], d, 1.0)
    return W / draws.shape[1]


def _part_errors(side: Side, q: str, W: np.ndarray, drop: Sequence[str] = (), chunk: int = 100) -> np.ndarray:
    """References x rows: each model row's error on part q (RMS over the entries both define, over the part's one-fly
    unit, kcq/population_fit.py's score_rows) against each weighted reference; entries of `drop` count as matching the
    reference (E138's 'without': the term out, the number of entries kept)."""
    X = side.fly_parts[q]
    fin = np.isfinite(X)
    X0 = np.where(fin, X, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = (W @ X0) / (W @ fin.astype(float))
    M = side.model[q]
    entries = _entries(side, q)
    di = [entries.index(e) for e in drop if e in entries]
    out = np.full((len(W), M.shape[0]), np.nan)
    for a in range(0, len(W), chunk):
        D = M[None, :, :] - mean[a:a + chunk, None, :]
        if di:
            D[:, :, di] = np.where(np.isfinite(D[:, :, di]), 0.0, D[:, :, di])
        ok = np.isfinite(D)
        with np.errstate(invalid='ignore', divide='ignore'):
            msd = np.where(ok, D ** 2, 0.0).sum(axis=2) / np.maximum(ok.sum(axis=2), 1)
        out[a:a + chunk] = np.where(ok.any(axis=2), np.sqrt(msd), np.nan) / side.units[q]
    return out


_ENTRIES: Dict[Tuple[str, str], List[str]] = {}


def _entries(side: Side, q: str) -> List[str]:
    """The entries of part q in the side's column order (from the saved flies' table)."""
    k = (side.key, q)
    if k not in _ENTRIES:
        fl = pd.read_parquet(flies_path(side.key)) if flies_path(side.key).exists() else None
        _ENTRIES[k] = list(dict.fromkeys(fl[fl['part'] == q]['entry'].astype(str))) if fl is not None else []
    return _ENTRIES[k]


def side_errors(side: Side, weights: np.ndarray, views: Optional[Sequence[str]] = ('strength',),
                error_col: str = 'error', drop: Sequence[str] = ()) -> np.ndarray:
    """Every model row's error against many references at once: `weights` (references x flies) gives each reference's
    mean as the weighted mean of the flies' parts (a half: 1 / k on its flies; a resample: counts / n), in all the
    flies' units; the views of `views` (None: every view), combined as kcq/population_fit.py's score_rows does
    ('error' over the three views read, 'error_with_matrix' with the matrix view too). Returns references x rows."""
    from kcq import population_fit as P
    view_err = {}
    for v in (*P.VIEWS, P.NAMED_VIEW):
        if views is not None and v not in views:
            continue
        qs = [f'{v}/{x}' for x in P.PARTS[v] if f'{v}/{x}' in side.fly_parts]
        if not qs:
            continue
        pe = [_part_errors(side, q, weights, drop) for q in qs]
        view_err[v] = np.sqrt(np.mean(np.square(pe), axis=0))
    three = [view_err[v] for v in P.VIEWS if v in view_err]
    if error_col == 'error_with_matrix':
        if P.NAMED_VIEW not in view_err:
            raise ValueError(f'{side.key}: no matrix view')
        return np.sqrt(np.mean(np.square(three + [view_err[P.NAMED_VIEW]]), axis=0))
    return np.sqrt(np.mean(np.square(three), axis=0))


def stage_error_col(key: str, views: Optional[Sequence[str]]) -> str:
    """E136's error column of a stage under a measure (`error_column`): the PN dendrites on four views where every
    view is read, else 'error'."""
    return 'error_with_matrix' if key.startswith('pn_dendrite/') and (views is None or 'matrix' in views) else 'error'


def setting_values(errors: np.ndarray, side: Side, wirings: Sequence[str] = T_WIRINGS) -> pd.DataFrame:
    """References x settings: each setting's value (the mean over `wirings` of each wiring's error; NaN unless all
    are read) from `side_errors`' output."""
    df = pd.DataFrame(errors, columns=side.rows)
    out = {}
    for name in dict.fromkeys(side.rows.get_level_values(0)):
        cols = [(name, w) for w in wirings if (name, w) in df.columns]
        out[name] = df[cols].mean(axis=1) if len(cols) == len(wirings) else pd.Series(np.nan, index=df.index)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# the parts

def _write_csv(df: pd.DataFrame, name: str) -> Path:
    path = out_dir() / name
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)
    return path


def provenance(part: str, lock: Mapping, **extra) -> None:
    """provenance_<part>.json: the script, `H.code_commits()`, host, time, the lock's sha256 and `extra`."""
    _write_json(out_dir() / f'provenance_{part}.json',
                dict(part=part, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
                     lock_sha256=lock.get('sha256'), **extra))


def setting_names(role: Optional[str] = 'setting', without_e135: bool = False) -> List[str]:
    """The settings' names (role 'setting' by default; None every one), in SETTINGS' order."""
    return [s.name for s in SETTINGS if (role is None or s.role == role) and not (without_e135 and s.source == 'E135')]


_SIDES: Dict[Tuple[str, bool], Side] = {}


def side(key: str, fit_evals: bool = False) -> Side:
    """`load_side` of every setting (SETTINGS' order), cached."""
    k = (key, fit_evals)
    if k not in _SIDES:
        _SIDES[k] = load_side(key, [s.name for s in SETTINGS], fit_evals=fit_evals)
    return _SIDES[k]


def megamat_key(stage: str) -> str:
    return f'{stage}/{PANEL}/on'


def lock_flies(lock: Mapping, key: str) -> List[str]:
    """The flies of a stage or panel in the lock (megamat stages by 'stage/megamat/on', panels by their class key)."""
    return list(lock['recordings'][key]['flies'])


def check_side_flies(sd: Side, flies: Sequence[str]) -> None:
    """The side's flies are the lock's, in its order."""
    if list(sd.flies) != [str(f) for f in flies]:
        raise RuntimeError(f'{sd.key}: the saved flies {sd.flies} differ from the lock\'s {list(flies)}')


@dataclass
class StageErrors:
    """Per megamat stage, the model rows' errors against references (references x rows), on one measure.

    Attributes:
        errors: per stage, references x rows.
        rows: the rows (setting, wiring), the same for every stage.
    """
    errors: Dict[str, np.ndarray]
    rows: pd.MultiIndex


def _aligned(sides: Mapping[str, Side]) -> pd.MultiIndex:
    rows = None
    for sd in sides.values():
        rows = sd.rows if rows is None else rows.intersection(sd.rows, sort=False)
    return rows


def stage_errors(weights: Mapping[str, np.ndarray], views: Optional[Sequence[str]] = ('strength',),
                 drop_soma: Sequence[str] = (), fit_evals: bool = False) -> StageErrors:
    """Per megamat stage, every row's error against the references of `weights` (per stage, references x flies), on
    the views of `views` with E136's error column (`stage_error_col`); the somata's entries of `drop_soma` left out
    (E138's rule)."""
    sides = {st: side(megamat_key(st), fit_evals) for st in weights}
    rows = _aligned(sides)
    out = {}
    for st, W in weights.items():
        sd = sides[st]
        e = side_errors(sd, W, views, stage_error_col(sd.key, views), drop=drop_soma if st == 'kc_soma' else ())
        idx = sd.rows.get_indexer(rows)
        out[st] = e[:, idx]
    return StageErrors(out, rows)


def rms_stages(se: StageErrors, stages: Sequence[str] = LOSS_STAGES) -> np.ndarray:
    """References x rows: each row's megamat stage error (the RMS of its stages' errors)."""
    return np.sqrt(np.mean(np.square([se.errors[s] for s in stages]), axis=0))


def values_of(errors: np.ndarray, rows: pd.MultiIndex, wirings: Sequence[str] = T_WIRINGS) -> pd.DataFrame:
    """References x settings: the mean over `wirings` of each wiring's value (NaN unless all are there)."""
    df = pd.DataFrame(errors, columns=rows)
    out = {}
    for name in dict.fromkeys(rows.get_level_values(0)):
        cols = [(name, w) for w in wirings if (name, w) in df.columns]
        out[name] = df[cols].mean(axis=1) if len(cols) == len(wirings) else pd.Series(np.nan, index=df.index)
    return pd.DataFrame(out)


def per_wiring(errors: np.ndarray, rows: pd.MultiIndex, name: str, wirings: Sequence[str] = T_WIRINGS) -> np.ndarray:
    """References x wirings: one setting's per-wiring values (NaN where a wiring is missing)."""
    df = pd.DataFrame(errors, columns=rows)
    return np.column_stack([df[(name, w)].to_numpy() if (name, w) in df.columns else np.full(len(df), np.nan)
                            for w in wirings])


def consistency(a: pd.Series, b: pd.Series) -> float:
    """Spearman ρ between two sets of settings' values (the settings both define)."""
    from scipy import stats
    ok = a.notna() & b.notna()
    if ok.sum() < 3:
        return np.nan
    return float(stats.spearmanr(a[ok].to_numpy(), b[ok].to_numpy())[0])


def split_weights(lock: Mapping, stage: str) -> Tuple[np.ndarray, np.ndarray]:
    """Halves a and b of every split of a megamat stage (`stage_splits`) as weights (splits x flies each)."""
    n = len(lock_flies(lock, megamat_key(stage)))
    sp = stage_splits(n)
    return half_weights(n, [a for a, _ in sp]), half_weights(n, [b for _, b in sp])


# ---------------------------------------------------------------------------
# design item 1: the check

def e136_classes() -> dict:
    return json.loads(results_file(E136, 'measure_strength', 'classes.json').read_text())


def part_check() -> dict:
    """Design item 1: on E136's own split, the shortcut's held-half megamat stage error on T for the family and A1
    against E136's score (classes.json, measure strength) and its per-stage values beside; passes within CHECK_TOL."""
    lock, _ = require_all()
    W = {}
    for st in LOSS_STAGES:
        n = len(lock_flies(lock, megamat_key(st)))
        check_side_flies(side(megamat_key(st)), lock_flies(lock, megamat_key(st)))
        fit, held = e136_split(n)
        W[st] = half_weights(n, [held, fit, tuple(range(n))])
    se = stage_errors(W)
    vals = {st: values_of(se.errors[st], se.rows) for st in LOSS_STAGES}
    vals['stage'] = values_of(rms_stages(se), se.rows)
    c = e136_classes()['classes']
    rows = []
    for arm in [a for a in ('A0', 'A1', 'A2', 'A3', 'A4') if a in c]:
        rec = c[arm]
        for what, scored_of in (('stage', lambda r: r['megamat']),
                                *((st, (lambda s: lambda r: r['stages'][s])(st)) for st in LOSS_STAGES)):
            for name, field in (('family', 'ref'), (f'E136:{arm}', 'arm')):
                scored = float(scored_of(rec)[field])
                rows.append(dict(setting=name, statistic=what, e136_held=scored,
                                 shortcut_held=float(vals[what][name].iloc[0]),
                                 shortcut_fit=float(vals[what][name].iloc[1]),
                                 shortcut_all=float(vals[what][name].iloc[2]), via=arm))
    df = pd.DataFrame(rows).drop_duplicates(['setting', 'statistic'])
    df['diff'] = df['shortcut_held'] - df['e136_held']
    _write_csv(df, 'check.csv')
    main = df[(df.statistic == 'stage') & df.setting.isin(CHECK_ARMS)]
    passed = bool(len(main) == len(CHECK_ARMS) and (main['diff'].abs() <= CHECK_TOL).all())
    out = dict(passed=passed, tol=CHECK_TOL, arms=list(CHECK_ARMS),
               values={r.setting: dict(e136_held=r.e136_held, shortcut_held=r.shortcut_held, diff=r.diff)
                       for r in main.itertuples()},
               largest_abs_diff_any=float(df['diff'].abs().max()),
               fallback='not needed' if passed else f'the splits are read E136\'s way on {N_FALLBACK} joint splits')
    rc = reclass_on_e136_split(lock)
    _write_csv(rc, 'check_classes.csv')
    out['reclass_on_e136_split'] = dict(n=int(len(rc)), reproduced=int(rc['same'].sum()),
                                        differ=rc[~rc['same']][['table', 'arm', 'cls_scored', 'cls_shortcut']
                                                               ].to_dict('records'))
    _write_json(out_dir() / 'check.json', out)
    provenance('check', lock)
    return out


def reclass_on_e136_split(lock: Mapping) -> pd.DataFrame:
    """Beside the check: every class of the class tables recomputed (`reclass`) on the shortcut's values against E136's
    held-out half (the held-out panels against all their flies, as scored), against the class as scored."""
    W = {}
    for st in LOSS_STAGES:
        n = len(lock_flies(lock, megamat_key(st)))
        W[megamat_key(st)] = half_weights(n, [e136_split(n)[1]])
    for k in HOLDS_PANELS + BESIDE_PANELS:
        n = len(lock_flies(lock, k))
        W[k] = half_weights(n, [tuple(range(n))])
    rv = weighted_values(lock, W)
    rows = []
    for exp, path, pre in CLASS_TABLES:
        tab = json.loads(results_file(exp, *path).read_text())
        for arm, rec in tab['classes'].items():
            cls, comps = reclass(rec, setting_of(exp, arm), rv, exp)
            m = comps.get('megamat')
            rows.append(dict(table=exp, arm=arm, cls_scored=rec.get('cls'), cls_shortcut=cls[0],
                             same=bool(cls[0] == rec.get('cls')),
                             change_scored=_scored_change(rec, 'megamat'),
                             change_shortcut=np.nan if m is None else float(m['change'][0])))
    return pd.DataFrame(rows)


def fallback_splits(n: int = N_FALLBACK, seed: int = FALLBACK_SEED) -> pd.DataFrame:
    """The check's fallback: N_FALLBACK seeded joint splits read E136's way (each half's own curve and noise); only if
    the check fails."""
    raise NotImplementedError('the check passed; the fallback is not needed (written only if it fails)')


def require_check() -> dict:
    """check.json, which must pass (the shortcut stands for E136's reading); RuntimeError otherwise."""
    p = out_dir() / 'check.json'
    if not p.exists():
        raise RuntimeError('run --part check first')
    c = json.loads(p.read_text())
    if not c['passed']:
        raise RuntimeError('the read-once check failed: the splits are to be read E136\'s way (fallback_splits)')
    return c


# ---------------------------------------------------------------------------
# design item 2: megamat on every split

def scored_sign(name: str, va: Optional[pd.Series] = None) -> float:
    """The sign of a setting's megamat change from the family as scored: E136's, E139's or E141's class table (the
    first that has it), else the shortcut's against E136's held-out half (`va`, the values of E136's split)."""
    for exp, path, pre in CLASS_TABLES:
        if not name.startswith(pre) and not (name == 'E136:A1'):
            continue
        p = results_file(exp, *path)
        if not p.exists():
            continue
        arm = name.split(':', 1)[1]
        rec = json.loads(p.read_text())['classes'].get(arm, {})
        ch = (rec.get('megamat') or {}).get('change')
        if ch is not None and np.isfinite(ch):
            return float(np.sign(ch))
    if va is not None and name in va and 'family' in va:
        return float(np.sign(va[name] - va['family']))
    return np.nan


def spreads_scored() -> Dict[str, float]:
    """The spreads as scored (E136's classes.json: the family's SD over its eight wirings, held-out half)."""
    return {k: float(v) for k, v in e136_classes()['spreads'].items()}


def beyond_bar(changes: np.ndarray, sign: float, spread: float) -> np.ndarray:
    """References x wirings of per-wiring changes -> whether the mean change is beyond E136's bar in the direction
    `sign` (E136's compare: se the larger of half the spread and the SD of the changes over sqrt(n))."""
    n = changes.shape[1]
    mean = changes.mean(axis=1)
    sd = changes.std(axis=1, ddof=1) if n > 1 else np.zeros(len(changes))
    se = np.maximum(spread / 2.0 if np.isfinite(spread) else 0.0, sd / np.sqrt(n))
    return np.isfinite(mean) & (np.sign(mean) == sign) & (np.abs(mean) > IMPROVE_SE * se)


def joint_weights(lock: Mapping) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], pd.DataFrame]:
    """Per stage, halves a and b of every joint split (the lock's rows) as weights (joint splits x flies)."""
    js = pd.DataFrame(lock['splits']['joint']['rows'])
    wa, wb = {}, {}
    for st in LOSS_STAGES:
        n = len(lock_flies(lock, megamat_key(st)))
        sp = stage_splits(n)
        wa[st] = half_weights(n, [sp[i][0] for i in js[st]])
        wb[st] = half_weights(n, [sp[i][1] for i in js[st]])
    return wa, wb, js


def part_splits(drop_soma: Sequence[str] = (), suffix: str = '') -> dict:
    """Design item 2: each setting's error at each stage against both halves of every split, consistency by stage and
    split, the megamat stage error on the joint splits, and per setting the share of joint splits on which its change
    from the family is beyond E136's bar against both halves (with `drop_soma`, the somata's entries left out: item 7)."""
    lock, _ = require_all()
    require_check()
    names, wo = setting_names(), setting_names(without_e135=True)
    stages = ('kc_soma',) if drop_soma else LOSS_STAGES
    vrows, crows = [], []
    for st in stages:
        Wa, Wb = split_weights(lock, st)
        se = stage_errors({st: np.vstack([Wa, Wb])}, drop_soma=drop_soma)
        v = values_of(se.errors[st], se.rows)
        n_sp = len(Wa)
        va, vb = v.iloc[:n_sp].reset_index(drop=True), v.iloc[n_sp:].reset_index(drop=True)
        n = len(lock_flies(lock, megamat_key(st)))
        e136_i = split_index(n, e136_split(n))
        for i in range(n_sp):
            crows.append(dict(stage=st, split=i, e136=i == e136_i, n_settings=int(va.loc[i, names].notna().sum()),
                              rho=consistency(va.loc[i, names], vb.loc[i, names]),
                              rho_without_e135=consistency(va.loc[i, wo], vb.loc[i, wo])))
            for nm in v.columns:
                vrows.append(dict(stage=st, split=i, setting=nm, half_a=float(va.loc[i, nm]), half_b=float(vb.loc[i, nm])))
    cons = pd.DataFrame(crows)
    _write_csv(cons, f'splits_consistency{suffix}.csv')
    _write_csv(pd.DataFrame(vrows), f'splits_stage_values{suffix}.csv')
    # the joint splits: the megamat stage error and each stage on them
    wa, wb, js = joint_weights(lock)
    sa = stage_errors(wa, drop_soma=drop_soma)
    sb = stage_errors(wb, drop_soma=drop_soma)
    ea, eb = rms_stages(sa), rms_stages(sb)
    va, vb = values_of(ea, sa.rows), values_of(eb, sb.rows)
    jrows = []
    for j in range(len(js)):
        r = dict(joint=int(js['joint'].iloc[j]), e136=bool(js['e136'].iloc[j]),
                 rho=consistency(va.loc[j, names], vb.loc[j, names]),
                 rho_without_e135=consistency(va.loc[j, wo], vb.loc[j, wo]))
        jrows.append(r)
    jdf = pd.DataFrame(jrows)
    _write_csv(jdf, f'splits_joint_consistency{suffix}.csv')
    # the share beyond the bar against both halves, per setting, on the megamat stage error and on each stage
    sp = spreads_scored()
    e136_row = int(np.flatnonzero(js['e136'].to_numpy())[0])
    srows = []
    stats_of = [('stage', ea, eb, sp['stage'])] + [(st, sa.errors[st], sb.errors[st], sp[st]) for st in LOSS_STAGES]
    va_e136 = values_of(ea, sa.rows).iloc[e136_row]
    drawn = ~js['e136'].to_numpy()
    for nm in setting_names(None):
        if nm == 'family':
            continue
        sign = scored_sign(nm, va_e136)
        rec = dict(setting=nm, sign=sign)
        for what, A, B, spread in stats_of:
            ca = per_wiring(A, sa.rows, nm) - per_wiring(A, sa.rows, 'family')
            cb = per_wiring(B, sb.rows, nm) - per_wiring(B, sb.rows, 'family')
            both = beyond_bar(ca, sign, spread) & beyond_bar(cb, sign, spread)
            rec[f'{what}_share_both'] = float(both[drawn].mean())
            rec[f'{what}_e136_both'] = bool(both[e136_row])
            rec[f'{what}_change_a_median'] = float(np.median(ca.mean(axis=1)[drawn]))
            rec[f'{what}_change_b_median'] = float(np.median(cb.mean(axis=1)[drawn]))
        srows.append(rec)
    shares = pd.DataFrame(srows)
    _write_csv(shares, f'splits_beyond_bar{suffix}.csv')
    summ = {}
    for st, g in cons.groupby('stage', sort=False):
        summ[st] = _tells(g['rho'].to_numpy(), g['rho_without_e135'].to_numpy(), float(g.loc[g.e136, 'rho'].iloc[0]))
    jd = jdf[~jdf.e136]
    summ['megamat_stage_error'] = _tells(jd['rho'].to_numpy(), jd['rho_without_e135'].to_numpy(),
                                         float(jdf.loc[jdf.e136, 'rho'].iloc[0]))
    out = dict(stages=summ, n_settings=len(names), suffix=suffix, drop_soma=list(drop_soma))
    if not drop_soma:
        out['fit_evals_beside'] = fit_evals_consistency(lock)
    _write_json(out_dir() / f'splits{suffix}.json', out)
    provenance(f'splits{suffix}', lock)
    return out


def fit_evals_consistency(lock: Mapping) -> dict:
    """Beside item 2 (a description): consistency over the fit evaluations this host holds (FIT_EVALS_BESIDE, read on
    S; each evaluation's value the mean over its S wirings) at each stage on every split and on the joint splits."""
    files = list(views_dir().glob('fitevals_*.parquet'))
    if not files:
        return dict(note='no fit evaluations held on this host')
    out = {}
    for st in LOSS_STAGES:
        Wa, Wb = split_weights(lock, st)
        sd = side(megamat_key(st), fit_evals=True)
        keep = [i for i, (nm, _) in enumerate(sd.rows) if str(nm).startswith('fit:')]
        rows = sd.rows[keep]
        e = side_errors(sd, np.vstack([Wa, Wb]), ('strength',), stage_error_col(sd.key, ('strength',)))[:, keep]
        v = values_of(e, rows, S_WIRINGS)
        va, vb = v.iloc[:len(Wa)].reset_index(drop=True), v.iloc[len(Wa):].reset_index(drop=True)
        rho = np.array([consistency(va.loc[i], vb.loc[i]) for i in range(len(Wa))])
        out[st] = dict(n_evals=int(v.shape[1]), median=float(np.nanmedian(rho)) if np.isfinite(rho).any() else None,
                       min=float(np.nanmin(rho)) if np.isfinite(rho).any() else None,
                       note='every evaluation has the same PN stage' if st == 'pn_dendrite' else '')
    return out


def tells_apart(rho: Sequence[float]) -> bool:
    """"A stage tells the saved models apart": consistency TELLS_MEDIAN or more at the median split and TELLS_MIN or
    more on every split."""
    r = np.asarray(rho, dtype=float)
    r = r[np.isfinite(r)]
    return bool(len(r) and np.median(r) >= TELLS_MEDIAN and np.min(r) >= TELLS_MIN)


def change_stands(lo: float, hi: float, share_both: Optional[float] = None) -> bool:
    """"A change stands through the flies": its interval excludes zero; on megamat (`share_both` given) it is also
    beyond the bar against both halves on at least STANDS_SHARE of the joint splits."""
    ex = bool(lo > 0 or hi < 0)
    if share_both is None:
        return ex
    return bool(ex and np.isfinite(share_both) and share_both >= STANDS_SHARE)


def class_stands(share: float) -> bool:
    """"A class stands through the flies": the scored class on at least CLASS_SHARE of the resamples."""
    return bool(np.isfinite(share) and share >= CLASS_SHARE)


def outcome_letters(tells: Mapping[str, bool], reads_calyx: bool, a1_fall: Sequence[bool], a1_holds: Sequence[bool],
                    near_shares: Sequence[float], count: bool, panel: bool) -> List[str]:
    """The entry's outcomes (a) to (e) from the classes: (a) the PN stages tell the saved models apart, the claws and
    somata do not, and megamat does not read the calyx; (b) all four stages tell them apart; (c) A1's megamat fall or
    its "holds" does not stand through the flies (in any class table); (d) a near-tie class is read on fewer than
    NEAR_TIE_SHARE of the resamples; (e) the count's or the panel's."""
    out = []
    if all(tells.get(st) for st in PN_STAGES) and not any(tells.get(st) for st in ('kc_claw', 'kc_soma')) \
            and not reads_calyx:
        out.append('a')
    if all(tells.get(st) for st in LOSS_STAGES):
        out.append('b')
    if not all(a1_fall) or not all(a1_holds):
        out.append('c')
    if any(sh < NEAR_TIE_SHARE for sh in near_shares):
        out.append('d')
    if count or panel:
        out.append('e')
    return out


def _tells(rho: np.ndarray, rho_wo: np.ndarray, e136: float) -> dict:
    """The consistency summary of one stage over its splits and its class "tells the saved models apart"."""
    r = rho[np.isfinite(rho)]
    med, mn = (float(np.median(r)), float(np.min(r))) if len(r) else (np.nan, np.nan)
    rw = rho_wo[np.isfinite(rho_wo)]
    return dict(n_splits=int(len(rho)), median=med, min=mn, max=float(np.max(r)) if len(r) else np.nan,
                share_ge_high=float(np.mean(r >= TELLS_MEDIAN)) if len(r) else np.nan,
                share_lt_low=float(np.mean(r < TELLS_MIN)) if len(r) else np.nan, e136_split=e136,
                median_without_e135=float(np.median(rw)) if len(rw) else np.nan,
                min_without_e135=float(np.min(rw)) if len(rw) else np.nan,
                tells_apart=tells_apart(rho), tells_apart_without_e135=tells_apart(rho_wo))


# ---------------------------------------------------------------------------
# design item 3: megamat's calyx

def part_calyx(drop_soma: Sequence[str] = (), suffix: str = '') -> dict:
    """Design item 3: consistency over the nine settings at pn_inhadd 1000 and the five at the family's normalization,
    at the boutons, claws and somata, on every split; the sign of each calyx contrast at the claws and somata, half
    against half (with `drop_soma`: the somata without those entries, item 7)."""
    lock, _ = require_all()
    require_check()
    stages = ('kc_soma',) if drop_soma else CALYX_STAGES
    crows, prows = [], []
    for st in stages:
        Wa, Wb = split_weights(lock, st)
        se = stage_errors({st: np.vstack([Wa, Wb])}, drop_soma=drop_soma)
        v = values_of(se.errors[st], se.rows)
        n_sp = len(Wa)
        va, vb = v.iloc[:n_sp].reset_index(drop=True), v.iloc[n_sp:].reset_index(drop=True)
        n = len(lock_flies(lock, megamat_key(st)))
        e136_i = split_index(n, e136_split(n))
        for g, members in CALYX_GROUPS.items():
            for i in range(n_sp):
                crows.append(dict(stage=st, group=g, split=i, e136=i == e136_i,
                                  rho=consistency(va.loc[i, list(members)], vb.loc[i, list(members)])))
            if st == 'pn_bouton':
                continue
            for x, y in itertools.combinations(members, 2):
                sa = np.sign(va[x] - va[y]).to_numpy()
                sb = np.sign(vb[x] - vb[y]).to_numpy()
                prows.append(dict(stage=st, group=g, a=x, b=y, share_splits_agree=float(np.mean(sa == sb)),
                                  e136_agree=bool(sa[e136_i] == sb[e136_i]),
                                  diff_a_median=float(np.median(va[x] - va[y])),
                                  diff_b_median=float(np.median(vb[x] - vb[y]))))
    cons, pairs = pd.DataFrame(crows), pd.DataFrame(prows)
    _write_csv(cons, f'calyx_consistency{suffix}.csv')
    _write_csv(pairs, f'calyx_contrasts{suffix}.csv')
    summ = {}
    for (st, g), x in cons.groupby(['stage', 'group'], sort=False):
        r = x['rho'].to_numpy()
        summ[f'{st}/{g}'] = dict(n_splits=int(len(r)), median=float(np.nanmedian(r)), min=float(np.nanmin(r)),
                                 e136_split=float(x.loc[x.e136, 'rho'].iloc[0]),
                                 tells_apart=tells_apart(r))
    for (st, g), x in pairs.groupby(['stage', 'group'], sort=False):
        summ[f'{st}/{g}'].update(n_pairs=int(len(x)), share_pairs_agree_mean=float(x['share_splits_agree'].mean()),
                                 share_pairs_agree_on_every_split=float((x['share_splits_agree'] == 1).mean()))
    reads = [summ.get(f'{st}/inhadd1000', {}).get('tells_apart') for st in ('kc_claw', 'kc_soma') if st in stages]
    out = dict(groups=summ, reads_calyx=bool(all(reads)) if reads else None, suffix=suffix)
    _write_json(out_dir() / f'calyx{suffix}.json', out)
    provenance(f'calyx{suffix}', lock)
    return out


# ---------------------------------------------------------------------------
# design item 4: intervals over flies and the classes recomputed

def resamples(lock: Mapping, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> Dict[str, np.ndarray]:
    """The resamples' draws per megamat stage ('stage/megamat/on') and held-out panel (class key): n_boot x flies
    indices, one generator, stages in LOSS_STAGES order then the panels in HOLDS_PANELS + BESIDE_PANELS order."""
    rng = np.random.default_rng(seed)
    out = {}
    for key in [megamat_key(s) for s in LOSS_STAGES] + list(HOLDS_PANELS + BESIDE_PANELS):
        n = len(lock_flies(lock, key))
        out[key] = rng.integers(0, n, size=(n_boot, n))
    return out


def _values_by_wiring(errors: np.ndarray, rows: pd.MultiIndex) -> Dict[str, np.ndarray]:
    """Per setting, references x T wirings."""
    return {nm: per_wiring(errors, rows, nm) for nm in dict.fromkeys(rows.get_level_values(0))}


@dataclass
class Resampled:
    """Every value a class reads, per reference (resample, or the unresampled flies at index -1 of each array).

    Attributes:
        mega: per statistic ('stage', 'downstream', each stage), per setting, references x T wirings (strength).
        guard: per stage, per setting, references x T wirings (every view; E136's error column).
        panels: per class key (HOLDS_PANELS + BESIDE_PANELS), per setting, references x T wirings (strength).
        guard_panels: the same on every view.
    """
    mega: Dict[str, Dict[str, np.ndarray]]
    guard: Dict[str, Dict[str, np.ndarray]]
    panels: Dict[str, Dict[str, np.ndarray]]
    guard_panels: Dict[str, Dict[str, np.ndarray]]


def resampled_values(lock: Mapping, draws: Mapping[str, np.ndarray], drop_soma: Sequence[str] = ()) -> Resampled:
    """Every value a class reads against each resample (`draws`), with the unresampled flies appended as the last
    reference."""
    W = {}
    for key, d in draws.items():
        n = d.shape[1]
        W[key] = np.vstack([boot_weights(n, d), np.full((1, n), 1.0 / n)])
    return weighted_values(lock, W, drop_soma)


def weighted_values(lock: Mapping, W: Mapping[str, np.ndarray], drop_soma: Sequence[str] = ()) -> Resampled:
    """Every value a class reads against the references of `W` (per megamat stage 'stage/megamat/on' and held-out
    panel class key, references x flies weights)."""
    ws = {st: W[megamat_key(st)] for st in LOSS_STAGES}
    se = stage_errors(ws, drop_soma=drop_soma)
    mega = {'stage': _values_by_wiring(rms_stages(se), se.rows),
            'downstream': _values_by_wiring(rms_stages(se, CALYX_STAGES), se.rows)}
    mega.update({st: _values_by_wiring(se.errors[st], se.rows) for st in LOSS_STAGES})
    ge = stage_errors(ws, views=None)
    guard = {st: _values_by_wiring(ge.errors[st], ge.rows) for st in LOSS_STAGES}
    panels, gpanels = {}, {}
    for key in HOLDS_PANELS + BESIDE_PANELS:
        sd = side(side_key(key))
        check_side_flies(sd, lock_flies(lock, key))
        panels[key] = _values_by_wiring(side_errors(sd, W[key], ('strength',)), sd.rows)
        gpanels[key] = _values_by_wiring(side_errors(sd, W[key], None), sd.rows)
    return Resampled(mega, guard, panels, gpanels)


def setting_of(table: str, arm: str) -> str:
    """The setting of a class table's name: 'A1' and E136's arms E136's, E139's names E139's reads (E141's B and K
    among them), E141's own variants E141's, the unrefitted family 'family'."""
    if arm in ('family', 'A1') or table == E136:
        return 'family' if arm == 'family' else f'E136:{arm}'
    if arm in ('B_inhadd1000', 'K_inhadd1000') or table == E139:
        return f'E139:{arm}'
    return f'E141:{arm}'


FLY_DEPENDENT_REASONS: Tuple[str, ...] = ('no better than the permuted control',)


def _cmp(a: Optional[np.ndarray], b: Optional[np.ndarray], scored: Mapping) -> Optional[dict]:
    """E136's compare on references x T wirings with the spread and the bar as scored: per reference the mean
    change, 'improves' (below -IMPROVE_SE scored standard errors) and 'worse' (above the scored spread). None if the
    scored comparison was not complete or a side is missing."""
    if a is None or b is None or not scored or not scored.get('complete'):
        return None
    ch = (a - b).mean(axis=1)
    se, spread = float(scored['se']), float(scored.get('spread', np.nan))
    return dict(change=ch, improves=ch < -IMPROVE_SE * se,
                worse=(ch > spread) if np.isfinite(spread) else np.zeros(len(ch), bool),
                better=(ch < -spread) if np.isfinite(spread) else np.zeros(len(ch), bool))


def reclass(rec: Mapping, name: str, rv: Resampled, table: str) -> Tuple[np.ndarray, Dict[str, dict]]:
    """One class recomputed on every reference of `rv`, as E136's classify ladder with the spreads, bars and
    fly-independent reasons as scored (`rec`, the class table's record), E139's and E141's failed holds kept. Returns
    (classes per reference, the comparisons by statistic)."""
    K = len(next(iter(rv.mega['stage'].values())))
    fam = 'family'
    if 'megamat' not in rec:          # reads stopped or missing: no class whatever the flies
        return np.array([rec.get('cls', 'no class')] * K, dtype=object), {}
    get = lambda d, nm: d.get(nm)
    comps: Dict[str, dict] = {}
    m = _cmp(get(rv.mega['stage'], name), get(rv.mega['stage'], fam), rec['megamat'])
    comps['megamat'] = m
    st = {s: _cmp(get(rv.mega[s], name), get(rv.mega[s], fam), rec['stages'].get(s, {})) for s in rec['stages']}
    comps.update({f'stage/{s}': c for s, c in st.items()})
    pn = {k: _cmp(get(rv.panels.get(k, {}), name), get(rv.panels.get(k, {}), fam), rec['panels'].get(k, {}))
          for k in rec.get('panels', {})}
    comps.update({f'panel/{k}': c for k, c in pn.items()})
    pb = {k: _cmp(get(rv.panels.get(k, {}), name), get(rv.panels.get(k, {}), fam), v)
          for k, v in (rec.get('panels_beside') or {}).items()}
    comps.update({f'beside/{k}': c for k, c in pb.items()})
    gd = rec.get('guard') or {}
    g_st = {s: _cmp(get(rv.guard[s], name), get(rv.guard[s], fam), gd[s]) for s in LOSS_STAGES if s in gd}
    g_pn = {k: _cmp(get(rv.guard_panels.get(k, {}), name), get(rv.guard_panels.get(k, {}), fam), gd[k])
            for k in HOLDS_PANELS + BESIDE_PANELS if k in gd}
    comps.update({f'guard/{s}': c for s, c in g_st.items()})
    comps.update({f'guard/{k}': c for k, c in g_pn.items()})
    fixed = [r for r in rec.get('reasons', []) if not r.startswith(FLY_DEPENDENT_REASONS)]
    if m is None or any(c is None for c in st.values()):
        return np.array(['no class'] * K, dtype=object), comps
    cf = rec.get('control_fall')
    no_better = (~(-m['change'] > cf)) if cf is not None and np.isfinite(cf) else np.zeros(K, bool)
    stage_worse = np.any([c['worse'] for c in st.values()], axis=0) if st else np.zeros(K, bool)
    guard_worse = np.any([c['worse'] for c in g_st.values() if c is not None], axis=0) if g_st else np.zeros(K, bool)
    complete = {k: c for k, c in pn.items() if c is not None}
    panel_worse = np.any([c['worse'] for c in complete.values()], axis=0) if complete else np.zeros(K, bool)
    for k, c in g_pn.items():
        if c is not None and k in HOLDS_PANELS:
            panel_worse = panel_worse | c['worse']
    panel_missing = len(complete) < len(rec.get('panels', {})) or bool(rec.get('panel_missing'))
    improves = m['improves'] & ~stage_worse & ~guard_worse
    cls = np.empty(K, dtype=object)
    for i in range(K):
        if fixed or no_better[i]:
            cls[i] = 'no class'
        elif improves[i] and not panel_worse[i] and not panel_missing:
            cls[i] = 'holds'
        elif improves[i] and panel_worse[i]:
            cls[i] = 'megamat only'
        elif improves[i]:
            cls[i] = 'improves (held-out panels not all read)'
        elif m['improves'][i] and stage_worse[i]:
            cls[i] = "one stage at another's cost"
        elif m['improves'][i] and guard_worse[i]:
            cls[i] = "strength at another view's cost"
        else:
            cls[i] = 'does not improve'
    return cls, comps


def against_a1(rec: Mapping, name: str, rv: Resampled) -> Dict[str, dict]:
    """E141's comparisons of a variant with A1 on every reference (megamat and the held-out panels; spreads and bars
    as scored in its against_a1 record)."""
    vs = rec.get('against_a1') or {}
    out = {}
    if vs.get('megamat'):
        out['vs_a1/megamat'] = _cmp(rv.mega['stage'].get(name), rv.mega['stage'].get('E136:A1'), vs['megamat'])
    for k, sc in (vs.get('panels') or {}).items():
        out[f'vs_a1/panel/{k}'] = _cmp(rv.panels.get(k, {}).get(name), rv.panels.get(k, {}).get('E136:A1'), sc)
    return out


def e141_entry(rec: Mapping, cls: np.ndarray, vs: Mapping[str, dict]) -> Dict[str, np.ndarray]:
    """E141's added classes on every reference: "stands in for the rise" (a form that holds), "does more than A1" (a
    dose that holds and improves on A1's megamat stage error, or is better on a held-out panel by more than the
    family's spread with none worse)."""
    out = {}
    holds = cls == 'holds'
    pn = [c for k, c in vs.items() if k.startswith('vs_a1/panel/') and c is not None]
    mm = vs.get('vs_a1/megamat')
    if rec.get('stands_in') is not None:
        out['stands_in'] = holds
    if rec.get('does_more_than_a1') is not None:
        better = np.any([c['better'] for c in pn], axis=0) if pn else np.zeros(len(cls), bool)
        worse = np.any([c['worse'] for c in pn], axis=0) if pn else np.zeros(len(cls), bool)
        imp = mm['improves'] if mm is not None else np.zeros(len(cls), bool)
        out['does_more_than_a1'] = holds & (imp | (better & ~worse))
    return out


def part_intervals(drop_soma: Sequence[str] = (), suffix: str = '') -> dict:
    """Design item 4: every megamat and held-out panel change of the class tables (minus the family; E141's minus A1)
    with its interval over flies, the classes recomputed on each resample and the share that gives the scored class;
    the realism flag's interval (with `drop_soma`: the somata without those entries, item 7)."""
    lock, _ = require_all()
    require_check()
    draws = resamples(lock)
    rv = resampled_values(lock, draws, drop_soma)
    bar = pd.read_csv(out_dir() / f'splits_beyond_bar{suffix}.csv').set_index('setting') \
        if (out_dir() / f'splits_beyond_bar{suffix}.csv').exists() else pd.DataFrame()
    crows, krows = [], []
    for exp, path, pre in CLASS_TABLES:
        tab = json.loads(results_file(exp, *path).read_text())
        for arm, rec in tab['classes'].items():
            name = setting_of(exp, arm)
            cls, comps = reclass(rec, name, rv, exp)
            extra = {}
            if exp == E141 and 'against_a1' in rec:
                vs = against_a1(rec, name, rv)
                comps.update(vs)
                extra = e141_entry(rec, cls, vs)
            boot, allf = cls[:-1], cls[-1]
            krows.append(dict(table=exp, arm=arm, setting=name, cls_scored=rec.get('cls'), cls_all_flies=allf,
                              share_scored=float(np.mean(boot == rec.get('cls'))),
                              stands=class_stands(float(np.mean(boot == rec.get('cls')))),
                              shares={c: float(np.mean(boot == c)) for c in sorted(set(boot))},
                              **{f'{k}_scored': rec.get(k) for k in extra},
                              **{f'{k}_share_scored': float(np.mean(v[:-1] == bool(rec.get(k)))) for k, v in extra.items()},
                              near_tie=(exp, arm) in NEAR_TIES))
            for what, c in comps.items():
                if c is None:
                    continue
                ch = c['change']
                lo, hi = np.quantile(ch[:-1], CI)
                scored = _scored_change(rec, what)
                row = dict(table=exp, arm=arm, setting=name, statistic=what, scored=scored, all_flies=float(ch[-1]),
                           lo=float(lo), hi=float(hi), excludes_zero=bool(lo > 0 or hi < 0),
                           share_improves=float(np.mean(c['improves'][:-1])), share_worse=float(np.mean(c['worse'][:-1])))
                if what == 'megamat' or what.startswith('stage/'):
                    col = 'stage' if what == 'megamat' else what.split('/', 1)[1]
                    sh = bar.loc[name, f'{col}_share_both'] if name in bar.index and f'{col}_share_both' in bar.columns \
                        else np.nan
                    row.update(share_beyond_bar_both=float(sh), stands=change_stands(lo, hi, float(sh)))
                else:
                    row['stands'] = change_stands(lo, hi)
                crows.append(row)
    ch_df = pd.DataFrame(crows)
    cl_df = pd.DataFrame(krows)
    cl_df['shares'] = cl_df['shares'].map(json.dumps)
    _write_csv(ch_df, f'intervals_changes{suffix}.csv')
    _write_csv(cl_df, f'class_shares{suffix}.csv')
    out = dict(n_boot=N_BOOT, ci=list(CI), suffix=suffix, n_changes=len(ch_df), n_classes=len(cl_df))
    if not drop_soma:
        out['realism_flag'] = flag_interval(lock, draws['apl16/apl16/on'])
    _write_json(out_dir() / f'intervals{suffix}.json', out)
    provenance(f'intervals{suffix}', lock)
    return out


def _scored_change(rec: Mapping, what: str) -> float:
    """The change of a comparison as scored in the class table (NaN if absent)."""
    try:
        if what == 'megamat':
            return float(rec['megamat']['change'])
        kind, _, key = what.partition('/')
        if kind == 'stage':
            return float(rec['stages'][key]['change'])
        if kind == 'panel':
            return float(rec['panels'][key]['change'])
        if kind == 'beside':
            return float(rec['panels_beside'][key]['change'])
        if kind == 'guard':
            return float(rec['guard'][key]['change'])
        if kind == 'vs_a1':
            sub, _, k = key.partition('/')
            a = rec['against_a1']
            return float(a['megamat']['change'] if sub == 'megamat' else a['panels'][k]['change'])
    except (KeyError, TypeError):
        pass
    return np.nan


def flag_interval(lock: Mapping, draws: np.ndarray) -> dict:
    """The realism flag's interval: the response rate ρ with APL on of FLAG_SETTING's APL16 read on FLAG_WIRING against
    the per-odor response fractions of resampled control recordings (each recording's KCs above its population's
    criterion, the bands' baseline-noise units), beside the band (kcq.apl_effect's RATE_RHO_MIN)."""
    from kcq import apl_effect as A
    fl = _fly_side('kc_soma', 'apl16', 'on')
    flies = [str(f) for f in fl.flies]
    if flies != lock_flies(lock, 'apl16/apl16/on'):
        raise RuntimeError('the control recordings differ from the lock\'s')
    fr = np.array([(fl.tables[f].astype(float) > float(fl.criterion)).mean(axis=0).reindex(H.APL16_ODORS).to_numpy()
                   for f in fl.flies])
    runs = setting_runs(SETTING_BY_NAME[FLAG_SETTING], 'apl16', (FLAG_WIRING,))
    r = runs[FLAG_WIRING].conditions['on']['responses']
    model = H.abbrev_columns(r.astype(float) > 0.5).mean().reindex(H.APL16_ODORS)
    rho_all = H.spearman(model, pd.Series(fr.mean(axis=0), index=H.APL16_ODORS))
    boots = np.array([H.spearman(model, pd.Series(fr[d].mean(axis=0), index=H.APL16_ODORS)) for d in draws])
    lo, hi = np.quantile(boots, CI)
    return dict(setting=FLAG_SETTING, wiring=FLAG_WIRING, rate_rho_on_all=float(rho_all),
                rate_rho_registry=float(H.spearman(model, A.frac_data('on', 'noise'))), lo=float(lo), hi=float(hi),
                band_min=float(A.RATE_RHO_MIN), share_below_band=float(np.mean(boots < A.RATE_RHO_MIN)))


# ---------------------------------------------------------------------------
# design item 5: how many flies

def disjoint_halves(n: int, k: int, n_draws: int = N_COUNT_DRAWS, seed: int = COUNT_SEED
                    ) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """`n_draws` pairs of disjoint halves of k of n flies: each draw a permutation of the flies (np.random.default_rng
    of `seed`), its first k one half and its next k the other."""
    if 2 * k > n:
        raise ValueError(f'{k=} {n=}')
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_draws):
        p = rng.permutation(n)
        out.append((tuple(sorted(int(i) for i in p[:k])), tuple(sorted(int(i) for i in p[k:2 * k]))))
    return out


def count_rows(sd: Side, label: str, key: str, names: Sequence[str], sizes: Sequence[int] = COUNT_SIZES) -> List[dict]:
    """Consistency over `names` between disjoint halves of each size the side's flies allow (`disjoint_halves`)."""
    n = len(sd.flies)
    rows = []
    for k in [k for k in sizes if 2 * k <= n]:
        hv = disjoint_halves(n, k)
        W = np.vstack([half_weights(n, [a for a, _ in hv]), half_weights(n, [b for _, b in hv])])
        col = stage_error_col(sd.key, ('strength',))
        v = setting_values(side_errors(sd, W, ('strength',), col), sd)
        va, vb = v.iloc[:len(hv)].reset_index(drop=True), v.iloc[len(hv):].reset_index(drop=True)
        have = [nm for nm in names if nm in v.columns and v[nm].notna().all()]
        rho = np.array([consistency(va.loc[i, have], vb.loc[i, have]) for i in range(len(hv))])
        r = rho[np.isfinite(rho)]
        rows.append(dict(key=key, label=label, n_flies=n, half=k, n_settings=len(have), n_draws=len(hv),
                         median=float(np.median(r)), lo=float(np.quantile(r, CI[0])), hi=float(np.quantile(r, CI[1])),
                         share_ge_high=float(np.mean(r >= COUNT_HIGH)), share_lt_low=float(np.mean(r < COUNT_LOW))))
    return rows


def _count_class(df: pd.DataFrame, panels: Sequence[str]) -> dict:
    """The count's and the panel's classes on `panels` (item 5's rows, by key)."""
    per = {}
    for key in panels:
        g = df[df.key == key].sort_values('half')
        if not len(g):
            continue
        two = g[g.half == 2]['median']
        big = g.iloc[-1]
        per[key] = dict(two=float(two.iloc[0]) if len(two) else np.nan, largest_half=int(big.half),
                        largest=float(big['median']),
                        count=bool(len(two) and two.iloc[0] < COUNT_LOW and big['median'] >= COUNT_HIGH
                                   and int(big.half) > 2),
                        panel=bool(len(two) and two.iloc[0] >= COUNT_HIGH),
                        first_half_at_high=next((int(r.half) for r in g.itertuples() if r.median >= COUNT_HIGH), None))
    return dict(per_panel=per, count=sum(v['count'] for v in per.values()) >= 2,
                panel=sum(v['panel'] for v in per.values()) >= 2)


def part_count() -> dict:
    """Design item 5: on the APL16 somata with APL on and off, the kiwi and control mixtures and validation2,
    consistency over the settings read there between disjoint halves of 2, 3, 4 and 5 flies, with megamat's four stages
    beside; the count's and the panel's classes on the three panels with many flies (and on the entry's three
    beside)."""
    lock, _ = require_all()
    require_check()
    names = setting_names()
    rows = []
    for key, label in COUNT_PANELS:
        sd = side(side_key(key))
        check_side_flies(sd, lock_flies(lock, key))
        rows += count_rows(sd, label, key, names)
    for st in LOSS_STAGES:
        sd = side(megamat_key(st))
        rows += count_rows(sd, f'megamat {st}', megamat_key(st), names)
    df = pd.DataFrame(rows)
    _write_csv(df, 'count_consistency.csv')
    out = dict(many=_count_class(df, MANY_FLY_PANELS), as_named=_count_class(df, MANY_FLY_PANELS_AS_NAMED),
               megamat_soma_two=float(df[(df.key == megamat_key('kc_soma')) & (df.half == 2)]['median'].iloc[0]))
    _write_json(out_dir() / 'count.json', out)
    provenance('count', lock)
    return out


# ---------------------------------------------------------------------------
# design item 6: the fifth soma fly

def _extra_fly_parts(stage: str, sd: Side) -> Dict[str, np.ndarray]:
    """The same-fly recording's parts at a megamat stage (viki's population kc_soma or kc_claw / megamat / same_fly),
    read as the stage's flies are (E134's fly_side_read rule, the reference's odors and bins), per part 1 x entries in
    the side's entry order."""
    from data_prep import stage_population_data as D
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    pop = D.load_population(stage, PANEL, 'same_fly')
    fp = P.fly_populations(pop)
    if len(fp) != 1:
        raise RuntimeError(f'{stage} same_fly: {len(fp)} flies')
    odors = _entries(sd, 'strength/rate')
    v = P.population_views(fp[0], odors, B.breadth_edges(len(odors)), list(sd.fly_parts))
    return {q: v[q].reindex(_entries(sd, q)).to_numpy(dtype=float)[None, :] for q in sd.fly_parts if q in v}


def _with_fly(sd: Side, parts: Mapping[str, np.ndarray], name: str) -> Side:
    """The side with one more fly appended (its units unchanged: the side's flies')."""
    return Side(key=sd.key, flies=list(sd.flies) + [name],
                fly_parts={q: np.vstack([sd.fly_parts[q], parts[q]]) for q in sd.fly_parts if q in parts},
                units=dict(sd.units), model=dict(sd.model), rows=sd.rows)


def _fly_vs(sd: Side, i: int, ref: Sequence[int], views=('strength',)) -> float:
    """Fly i's parts scored against the mean of the flies `ref`, in the side's units."""
    one = Side(key=sd.key, flies=sd.flies, fly_parts=sd.fly_parts, units=sd.units,
               model={q: sd.fly_parts[q][i:i + 1] for q in sd.fly_parts}, rows=pd.MultiIndex.from_tuples([('fly', i)]))
    W = half_weights(len(sd.flies), [list(ref)])
    return float(side_errors(one, W, views, 'error')[0, 0])


def part_fifth() -> dict:
    """Design item 6: each setting's somata error against the fifth soma fly (the same-fly recording) in the four flies'
    one-fly units, beside each of the four against the other three, the fifth against the four and each setting against
    the four; the claw fly 2023-04-11 fly03 against its second acquisition, beside each claw fly against the other three.
    A description (one fly), never a class."""
    lock, _ = require_all()
    require_check()
    out, rows = {}, []
    for stage, own in (('kc_soma', None), ('kc_claw', '2023-04-11__fly03__megamat1')):
        sd0 = side(megamat_key(stage))
        check_side_flies(sd0, lock_flies(lock, megamat_key(stage)))
        name = lock_flies(lock, f'{stage}/{PANEL}/same_fly')[0]
        sd = _with_fly(sd0, _extra_fly_parts(stage, sd0), name)
        n0 = len(sd0.flies)
        five = n0
        for views, lab in ((('strength',), 'strength'), (None, 'three views')):
            if stage == 'kc_soma':
                W = half_weights(n0 + 1, [[five], list(range(n0))] + [[i] for i in range(n0)])
                v = setting_values(side_errors(sd, W, views, 'error'), sd)
                for nm in v.columns:
                    single = {f'against_{f}': float(v.loc[2 + i, nm]) for i, f in enumerate(sd0.flies)}
                    rows.append(dict(stage=stage, views=lab, setting=nm, against_fifth=float(v.loc[0, nm]),
                                     against_four=float(v.loc[1, nm]),
                                     against_one_of_four_mean=float(np.mean(list(single.values()))), **single))
            loo = {f: _fly_vs(sd, i, [j for j in range(n0) if j != i], views) for i, f in enumerate(sd0.flies)}
            rec = dict(each_against_others=loo, fifth_against_all=_fly_vs(sd, five, range(n0), views))
            if own is not None:
                oi = sd0.flies.index(own)
                rec.update(second_against_first=_fly_vs(sd, five, [oi], views),
                           first_against_second=_fly_vs(_with_fly(sd0, {q: sd.fly_parts[q][five:five + 1]
                                                                        for q in sd.fly_parts}, name), oi, [n0], views),
                           second_against_others=_fly_vs(sd, five, [j for j in range(n0) if j != oi], views),
                           claw_fly=own)
            out[f'{stage}/{lab}'] = rec
    df = pd.DataFrame(rows)
    _write_csv(df, 'fifth_fly_settings.csv')
    for lab, g in df.groupby('views', sort=False):
        g = g[g.setting.isin(setting_names())].set_index('setting')
        out[f'kc_soma/{lab}'].update(
            order_fifth_vs_four=consistency(g['against_fifth'], g['against_four']),
            order_fifth_vs_each={f: consistency(g['against_fifth'], g[f'against_{f}'])
                                 for f in side(megamat_key('kc_soma')).flies})
    out['note'] = 'one fly: a description, never a class'
    _write_json(out_dir() / 'fifth_fly.json', out)
    provenance('fifth', lock)
    return out


# ---------------------------------------------------------------------------
# design item 7: the megamat somata without aa

def part_no_aa() -> dict:
    """Design item 7: items 2 to 4 for the megamat somata with aa's entry left out (E138's rule, E140's outcome (d)):
    the somata's consistency on its splits and the megamat stage error on the joint splits, the calyx at the somata,
    and the intervals and class shares, files with the suffix '_no_aa'."""
    drop = (WITHOUT_ODOR,)
    return dict(splits=part_splits(drop, '_no_aa'), calyx=part_calyx(drop, '_no_aa'),
                intervals=part_intervals(drop, '_no_aa'))


# ---------------------------------------------------------------------------
# the classes and the outcome

def _read(name: str):
    p = out_dir() / name
    if not p.exists():
        raise RuntimeError(f'{name} is missing: run its part first')
    return json.loads(p.read_text()) if name.endswith('.json') else pd.read_csv(p)


def part_classes() -> dict:
    """The entry's classes and the outcome from the parts' tables: classes.json and results_tables.md."""
    lock, _ = require_all()
    check = require_check()
    sp, cal, cnt = _read('splits.json'), _read('calyx.json'), _read('count.json')
    ch, cl = _read('intervals_changes.csv'), _read('class_shares.csv')
    tells = {st: v['tells_apart'] for st, v in sp['stages'].items()}
    a1 = ch[(ch.setting == 'E136:A1') & (ch.statistic == 'megamat')]
    a1_fall = {r.table: bool(r.stands) for r in a1.itertuples()}
    a1_cls = {r.table: bool(r.stands) for r in cl[cl.setting == 'E136:A1'].itertuples() if r.cls_scored == 'holds'}
    near = cl[cl.near_tie.astype(bool)]
    near_changes = {f'{r.table}/{r.arm}': dict(cls=r.cls_scored, share=float(r.share_scored)) for r in near.itertuples()}
    out = dict(
        check=dict(passed=check['passed'], values=check['values']),
        tells_apart=tells,
        reads_calyx=cal['reads_calyx'],
        a1=dict(fall_stands=a1_fall, holds_stands=a1_cls),
        near_ties=near_changes,
        count=cnt['many'], count_as_named=cnt['as_named'],
    )
    args = dict(reads_calyx=bool(cal['reads_calyx']), a1_fall=list(a1_fall.values()), a1_holds=list(a1_cls.values()),
                near_shares=[v['share'] for v in near_changes.values()], count=bool(cnt['many']['count']),
                panel=bool(cnt['many']['panel']))
    letters = outcome_letters(tells, **args)
    tells_wo = {st: v.get('tells_apart_without_e135') for st, v in sp['stages'].items()}
    out['tells_apart_without_e135'] = tells_wo
    out['outcomes_without_e135'] = outcome_letters(tells_wo, **args)
    out['outcomes_as_named_count'] = outcome_letters(tells, **dict(args, count=bool(cnt['as_named']['count']),
                                                                   panel=bool(cnt['as_named']['panel'])))
    out['outcomes'] = letters
    out['outcome_texts'] = {k: OUTCOMES[k] for k in letters}
    _write_json(out_dir() / 'classes.json', out)
    write_tables()
    provenance('classes', lock, outcomes=letters)
    _write_json(out_dir() / 'provenance.json', dict(
        experiment=EXPERIMENT, script=SCRIPT, **H.code_commits(), host=_host(), time=_now(),
        lock_sha256=lock.get('sha256'), e134_lock_sha256=_e140()._sha256(E134_LOCK_COPY),
        parts=sorted(p.name for p in out_dir().glob('provenance_*.json')),
        read_summary=str(local_dir() / f'read_summary_{_host()}.csv'), views=str(views_dir()),
        note='no simulation: every saved run read once (part read), then splits and resamples of the flies'))
    return out


def _md(df: pd.DataFrame, digits: int = 3) -> str:
    """A markdown table of `df` (numbers rounded to `digits`, right-aligned)."""
    def fmt(x):
        if isinstance(x, (bool, np.bool_)):
            return 'yes' if x else 'no'
        if isinstance(x, (float, np.floating)):
            return '–' if not np.isfinite(x) else f'{x:.{digits}f}'
        return str(x)
    cols = list(df.columns)
    num = [all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))
               for v in df[c]) for c in cols]
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |',
             '|' + '|'.join('---:' if n else '---' for n in num) + '|']
    lines += ['| ' + ' | '.join(fmt(v) for v in r) + ' |' for r in df.itertuples(index=False)]
    return '\n'.join(lines)


STAGE_WORDS = {'pn_dendrite': 'PN dendrites', 'pn_bouton': 'PN boutons', 'kc_claw': 'KC claws', 'kc_soma': 'KC somata',
               'megamat_stage_error': 'megamat stage error (joint splits)'}


def write_tables() -> Path:
    """results_tables.md: the parts' tables in short (every number from the parts' files)."""
    lock = require_lock()
    chk, sp, cal = _read('check.json'), _read('splits.json'), _read('calyx.json')
    cnt, ff, iv = _read('count.json'), _read('fifth_fly.json'), _read('intervals.json')
    sp0, cal0 = _read('splits_no_aa.json'), _read('calyx_no_aa.json')
    cl, ch = _read('class_shares.csv'), _read('intervals_changes.csv')
    cl0, ch0 = _read('class_shares_no_aa.csv'), _read('intervals_changes_no_aa.csv')
    cc, bar = _read('check.csv'), _read('splits_beyond_bar.csv')
    out = [f'# E143 results tables', '',
           f'From the parts of {_now()} ({_host()}, queue {H.code_commits().get("queue_commit", "")[:7]}); lock '
           f'{lock.get("sha256", "")[:12]}. Strength by odor in one-fly units (lower is better), the megamat panel, the '
           'T wirings\' mean unless stated. Every model run is read once against all the flies; a half or a resample '
           'changes only the flies\' mean view.', '']
    # the check
    c = cc[cc.statistic == 'stage'][['setting', 'e136_held', 'shortcut_held', 'diff', 'shortcut_fit', 'shortcut_all']]
    out += ['## 1. The read-once check', '',
            f"Megamat stage error on T against E136's held-out half: E136's score against the shortcut (tolerance "
            f"{CHECK_TOL}); passed: {'yes' if chk['passed'] else 'no'}. Beside: the shortcut against the fit half and "
            "against all flies (E136's own all-flies values are 0.919 and 0.817).", '',
            _md(c.rename(columns=dict(e136_held="E136, held-out half", shortcut_held='shortcut, held-out half',
                                      diff='difference', shortcut_fit='shortcut, fit half',
                                      shortcut_all='shortcut, all flies'))), '',
            f"E136's, E139's and E141's classes recomputed on the shortcut's values against E136's held-out half: "
            f"{chk['reclass_on_e136_split']['reproduced']} of {chk['reclass_on_e136_split']['n']} reproduced "
            f"(differ: {', '.join(f'{d['arm']} ({d['cls_scored']} -> {d['cls_shortcut']})' for d in chk['reclass_on_e136_split']['differ']) or 'none'}).", '']
    # item 2
    rows = []
    for st, v in sp['stages'].items():
        rows.append(dict(stage=STAGE_WORDS.get(st, st), splits=v['n_splits'], median=v['median'], min=v['min'],
                         e136_split=v['e136_split'], median_wo_e135=v['median_without_e135'],
                         min_wo_e135=v['min_without_e135'], tells_apart=v['tells_apart'],
                         tells_apart_wo_e135=v['tells_apart_without_e135']))
    out += ['## 2. Megamat on every split', '',
            f"Consistency: Spearman ρ over the {sp['n_settings']} settings between their values against the two halves "
            "(yardstick 0.8 at the median split, 0.5 on every split). \"Without E135\" leaves out E135's four arms, "
            'whose PN stage was set from the dendrite flies.', '',
            _md(pd.DataFrame(rows).rename(columns=dict(e136_split="E136's split", median_wo_e135='median without E135',
                                                       min_wo_e135='min without E135', tells_apart='tells apart',
                                                       tells_apart_wo_e135='tells apart without E135'))), '']
    fe = sp.get('fit_evals_beside') or {}
    if 'note' not in fe:
        out += ['Beside (a description): consistency over the fit evaluations of E136\'s A3 held on this host (S '
                'wirings): ' + '; '.join(f"{STAGE_WORDS[st]} median {v['median']:.3f}, min {v['min']:.3f} "
                                         f"({v['n_evals']} evaluations)" if v['median'] is not None else
                                         f"{STAGE_WORDS[st]}: not defined ({v['note'] or 'ties'})"
                                         for st, v in fe.items()) + '.', '']
    b = bar[['setting', 'sign', 'stage_share_both', 'stage_e136_both', 'stage_change_a_median', 'stage_change_b_median']]
    out += ["Per setting, its change of the megamat stage error from the family (minus, mean over T): the median over "
            "the 200 joint splits against each half, and the share of joint splits on which it is beyond E136's bar "
            "against both halves in the scored direction.", '',
            _md(b.rename(columns=dict(sign='scored sign', stage_share_both='share beyond the bar, both halves',
                                      stage_e136_both="on E136's split", stage_change_a_median='change, half a',
                                      stage_change_b_median='change, half b'))), '']
    # item 3
    rows = [dict(stage_group=k, splits=v['n_splits'], median=v['median'], min=v['min'], e136_split=v['e136_split'],
                 pairs_agreeing=v.get('share_pairs_agree_mean', np.nan),
                 pairs_agreeing_every_split=v.get('share_pairs_agree_on_every_split', np.nan)) for k, v in
            cal['groups'].items()]
    out += ["## 3. Megamat's calyx", '',
            "Consistency over the nine settings at pn_inhadd 1000 ('inhadd1000') and the five at the family's "
            "normalization ('family_norm'); beside, the calyx contrasts (pairs within a group): the mean share of "
            f"splits on which a pair's sign agrees between the halves, and the share of pairs agreeing on every split. "
            f"Megamat reads the calyx: {'yes' if cal['reads_calyx'] else 'no'}.", '',
            _md(pd.DataFrame(rows).rename(columns=dict(stage_group='stage / group', e136_split="E136's split",
                                                       pairs_agreeing='pairs: share of splits agreeing',
                                                       pairs_agreeing_every_split='pairs agreeing on every split'))), '']
    # item 4
    m = ch[ch.statistic == 'megamat'][['table', 'arm', 'scored', 'all_flies', 'lo', 'hi', 'share_beyond_bar_both',
                                       'stands']].copy()
    m['table'] = m['table'].str[:4]
    out += ['## 4. Intervals over flies and the classes through the flies', '',
            f"The megamat stage error's change from the family (mean over T; lower is better) as scored (the held-out "
            f"half), against all flies, and its {CI[0]:.0%} to {CI[1]:.0%} interval over {N_BOOT} resamples of the "
            "flies (the four stages together); \"stands\" also asks the bar on both halves on at least 80 % of the "
            'joint splits.', '',
            _md(m.rename(columns=dict(scored='scored (held-out half)', all_flies='all flies',
                                      share_beyond_bar_both='beyond the bar, both halves'))), '']
    k = cl[['table', 'arm', 'cls_scored', 'cls_all_flies', 'share_scored', 'stands', 'shares']].copy()
    k['table'] = k['table'].str[:4]
    k['shares'] = k['shares'].map(lambda x: ', '.join(f'{c} {v:.2f}' for c, v in json.loads(x).items()))
    out += ['Each class recomputed on every resample (spreads, bars and the fly-independent reasons as scored):', '',
            _md(k.rename(columns=dict(cls_scored='class as scored', cls_all_flies='class against all flies',
                                      share_scored='share giving the scored class', shares='every class\'s share'))),
            '']
    focus = [('E136_stage_population_fit', 'A1')] + list(NEAR_TIES) + [(E141, 'D5.8_inhadd1000')]
    pr = ch[np.array([(t, a) in focus for t, a in zip(ch.table, ch.arm)]) & ~ch.statistic.str.startswith('guard')
            & (ch.statistic != 'megamat')][['table', 'arm', 'statistic', 'scored', 'lo', 'hi', 'share_improves',
                                            'share_worse']].copy()
    pr['table'] = pr['table'].str[:4]
    out += ["The stages' and held-out panels' changes behind those classes (A1 and the near-ties; minus the family, "
            "E141's 'vs_a1' minus A1): the interval and the shares of resamples on which the change improves (below "
            "the scored bar) or is worse than the scored spread.", '', _md(pr), '']
    rf = iv['realism_flag']
    out += [f"Realism flag: the shuffled matrix's response rate ρ with APL on, wiring {rf['wiring']}: "
            f"{rf['rate_rho_on_all']:.3f} against all control recordings (the band's floor {rf['band_min']:.2f}), "
            f"interval {rf['lo']:.3f} to {rf['hi']:.3f} over the recordings; below the band on "
            f"{rf['share_below_band']:.0%} of the resamples.", '']
    # item 5
    cd = _read('count_consistency.csv')[['label', 'n_flies', 'half', 'n_settings', 'median', 'lo', 'hi',
                                         'share_ge_high', 'share_lt_low']]
    out += ['## 5. How many flies', '',
            "Consistency over the settings read on each panel between disjoint halves of each size (200 draws): the "
            "median, its 5 % to 95 % range and the shares at 0.8 or more and under 0.5. The kiwi mixtures stand in "
            "for the control mixtures among the panels with many flies (the control mixtures have 4 flies as E136 "
            "reads them).", '',
            _md(cd.rename(columns=dict(label='panel', n_flies='flies', half='flies per half', n_settings='settings',
                                       share_ge_high='share >= 0.8', share_lt_low='share < 0.5'))), '',
            f"The count's: {'yes' if cnt['many']['count'] else 'no'}; the panel's: "
            f"{'yes' if cnt['many']['panel'] else 'no'} (on the entry's three, with the control mixtures: "
            f"{'yes' if cnt['as_named']['count'] else 'no'} / {'yes' if cnt['as_named']['panel'] else 'no'}).", '']
    # item 6
    f6 = _read('fifth_fly_settings.csv')
    f6 = f6[(f6.views == 'strength') & f6.setting.isin(setting_names())]
    so, cw = ff['kc_soma/strength'], ff['kc_claw/strength']
    rows = [dict(statistic='each of the four soma flies against the other three', value=', '.join(
        f'{v:.3f}' for v in so['each_against_others'].values())),
            dict(statistic='the fifth soma fly against the four', value=f"{so['fifth_against_all']:.3f}"),
            dict(statistic='settings against the fifth fly, median (range)',
                 value=f"{f6.against_fifth.median():.3f} ({f6.against_fifth.min():.3f} to {f6.against_fifth.max():.3f})"),
            dict(statistic='settings against the four, median (range)',
                 value=f"{f6.against_four.median():.3f} ({f6.against_four.min():.3f} to {f6.against_four.max():.3f})"),
            dict(statistic='settings against one of the four, mean over the four, median',
                 value=f"{f6.against_one_of_four_mean.median():.3f}"),
            dict(statistic='order of the settings, against the fifth vs against the four (ρ)',
                 value=f"{so['order_fifth_vs_four']:.3f}"),
            dict(statistic='the claw fly 2023-04-11 fly03 against its second acquisition',
                 value=f"{cw['second_against_first']:.3f}"),
            dict(statistic='each of the four claw flies against the other three', value=', '.join(
                f'{v:.3f}' for v in cw['each_against_others'].values())),
            dict(statistic='the second acquisition against the other three claw flies',
                 value=f"{cw['second_against_others']:.3f}")]
    out += ['## 6. The fifth soma fly (a description)', '',
            "Strength by odor in the four flies' one-fly units (the somata: Remy's four wild-type flies; the claws: "
            "the four claw flies).", '', _md(pd.DataFrame(rows)), '']
    # item 7
    rows = []
    for (lab, a, b) in (('KC somata', sp['stages']['kc_soma'], sp0['stages']['kc_soma']),
                        ('megamat stage error (joint splits)', sp['stages']['megamat_stage_error'],
                         sp0['stages']['megamat_stage_error']),
                        ('KC somata, nine at pn_inhadd 1000', cal['groups']['kc_soma/inhadd1000'],
                         cal0['groups']['kc_soma/inhadd1000'])):
        rows.append(dict(statistic=lab, median_with_aa=a['median'], min_with_aa=a['min'], median_without=b['median'],
                         min_without=b['min']))
    out += ["## 7. The megamat somata without aa's entry", '', _md(pd.DataFrame(rows).rename(columns=dict(
        median_with_aa='median, with aa', min_with_aa='min, with aa', median_without='median, without aa',
        min_without='min, without aa'))), '']
    k0 = cl0[['table', 'arm', 'cls_scored', 'share_scored']].merge(
        cl[['table', 'arm', 'share_scored']], on=['table', 'arm'], suffixes=('_without_aa', '_with_aa'))
    k0['table'] = k0['table'].str[:4]
    out += ['The classes through the flies with and without aa\'s entry at the somata (the share of resamples giving '
            'the scored class):', '', _md(k0), '']
    # the classes
    cls = _read('classes.json')
    out += ['## The classes and the outcome', '',
            f"Tells the saved models apart (with E135's four / without): " + '; '.join(
                f"{STAGE_WORDS.get(st, st)} {'yes' if v else 'no'} / "
                f"{'yes' if cls['tells_apart_without_e135'].get(st) else 'no'}" for st, v in cls['tells_apart'].items())
            + '.',
            f"Megamat reads the calyx: {'yes' if cls['reads_calyx'] else 'no'}. A1's megamat fall stands: "
            f"{cls['a1']['fall_stands']}; A1's \"holds\" stands: {cls['a1']['holds_stands']}.",
            'Near-ties (share of resamples giving the scored class): ' + ', '.join(
                f"{k} {v['cls']} {v['share']:.2f}" for k, v in cls['near_ties'].items()) + '.',
            f"Outcomes: {', '.join(cls['outcomes']) or 'none'} (without E135's four: "
            f"{', '.join(cls['outcomes_without_e135']) or 'none'}).", '']
    path = out_dir() / 'results_tables.md'
    path.write_text('\n'.join(out) + '\n')
    return path


PARTS: Tuple[str, ...] = ('lock', 'read', 'check', 'splits', 'calyx', 'intervals', 'count', 'fifth', 'no_aa', 'classes',
                          'all')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--jobs', type=int, default=1, help='--part read: processes')
    ap.add_argument('--settings', nargs='+', help='--part read: only these settings')
    a = ap.parse_args(argv)
    if a.part == 'lock':
        part_lock()
        return
    require_lock()
    if a.part == 'read':
        print(part_read(a.jobs, a.settings).to_string(), flush=True)
        return
    order = ('check', 'splits', 'calyx', 'intervals', 'count', 'fifth', 'no_aa', 'classes')
    fns = dict(check=part_check, splits=part_splits, calyx=part_calyx, intervals=part_intervals, count=part_count,
               fifth=part_fifth, no_aa=part_no_aa, classes=part_classes)
    for p in (order if a.part == 'all' else (a.part,)):
        print(json.dumps(fns[p](), indent=1, default=str)[:4000], flush=True)


if __name__ == '__main__':
    main()
