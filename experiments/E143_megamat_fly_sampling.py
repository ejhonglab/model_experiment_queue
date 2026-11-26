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

Outputs: <shared>/results/E143_megamat_fly_sampling/ (lock.json, check.json, splits_*.csv, calyx_*.csv,
intervals_*.csv, class_shares.csv, count_*.csv, fifth_*.csv, no_aa_*.csv, classes.json, results_tables.md,
provenance_<part>.json); outputs/E143_megamat_fly_sampling/ (views/, logs).

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

def part_check() -> dict:
    """Design item 1: on E136's own split, the shortcut's held-half megamat stage error on T for the family and A1
    against E136's score (classes.json, measure strength) and its per-stage values beside; passes within CHECK_TOL."""
    raise NotImplementedError


def fallback_splits(n: int = N_FALLBACK, seed: int = FALLBACK_SEED) -> pd.DataFrame:
    """The check's fallback: N_FALLBACK seeded joint splits read E136's way (each half's own curve and noise); only if
    the check fails."""
    raise NotImplementedError


def consistency(a: pd.Series, b: pd.Series) -> float:
    """Spearman ρ between two sets of settings' values (the settings both define)."""
    raise NotImplementedError


def part_splits() -> dict:
    """Design item 2: each setting's error at each stage against both halves of every split, consistency by stage and
    split, the megamat stage error on the joint splits, and per setting the share of joint splits on which its change
    from the family is beyond E136's bar against both halves."""
    raise NotImplementedError


def part_calyx() -> dict:
    """Design item 3: consistency over the nine settings at pn_inhadd 1000 and the five at the family's normalization,
    at the boutons, claws and somata, on every split; the sign of each calyx contrast at the claws and somata, half
    against half."""
    raise NotImplementedError


def resamples(n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> Dict[str, np.ndarray]:
    """The resamples' weights per stage and panel (n_boot x flies, the counts of each fly over its number)."""
    raise NotImplementedError


def reclass(scored: Mapping, values: Mapping) -> str:
    """One class recomputed on resampled values with the spreads, bars and fly-independent reasons as scored."""
    raise NotImplementedError


def part_intervals() -> dict:
    """Design item 4: every megamat and held-out panel change of the class tables with its interval over flies, the
    classes recomputed on each resample and the share that gives the scored class; the realism flag's interval."""
    raise NotImplementedError


def part_count() -> dict:
    """Design item 5: on the APL16 somata with APL on and off, the control mixtures and validation2, consistency over
    the settings between disjoint halves of 2, 3, 4 and 5 flies, with megamat's somata at halves of 2 beside."""
    raise NotImplementedError


def part_fifth() -> dict:
    """Design item 6: each setting's somata error against the fifth soma fly, in the four flies' units, beside each of
    the four against the other three; the claw fly against its second acquisition."""
    raise NotImplementedError


def part_no_aa() -> dict:
    """Design item 7: items 2 to 4 for the megamat somata with aa's entry left out."""
    raise NotImplementedError


def part_classes() -> dict:
    """The entry's classes and the outcome from the parts' tables; classes.json and results_tables.md."""
    raise NotImplementedError


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
