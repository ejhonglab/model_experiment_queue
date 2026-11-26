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
# the panels of item 5: (input/panel/condition, label); the three with many flies first (the count's and the panel's
# classes read those), validation2 beside
COUNT_PANELS: Tuple[Tuple[str, str], ...] = (('apl16/apl16/on', 'APL16 somata, APL on'),
                                             ('apl16/apl16/off', 'APL16 somata, APL off'),
                                             ('mixtures/cmix/on', 'control mixtures'),
                                             ('validation2/validation2/on', 'validation2'))
MANY_FLY_PANELS: Tuple[str, ...] = ('apl16/apl16/on', 'apl16/apl16/off', 'mixtures/cmix/on')
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
# the near-tie classes the entry names (outcome (d))
NEAR_TIES: Tuple[str, ...] = ('E141:D1.8_inhadd1000', 'E139:K_inhadd1000', 'E141:K_inhadd1000', 'E141:N_inhadd1000',
                              'E141:G_inhadd1000')

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
    raise NotImplementedError


def e136_split(n: int) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """E136's own split of n flies (hal's `split_halves(n, 1, 0)`: its fit half, then its held-out half)."""
    raise NotImplementedError


def joint_splits(sizes: Mapping[str, int], n_joint: int = N_JOINT, seed: int = JOINT_SEED) -> pd.DataFrame:
    """The joint splits: one split per stage (an index into `stage_splits`), `n_joint` drawn without replacement from
    their product (stage order LOSS_STAGES), then E136's own split ('e136'). Columns: joint, <stage> (the index), and
    whether the row is E136's own."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the lock

def lock_inputs() -> dict:
    """The inputs the lock records with sha256: the locks and data files every statistic rests on ('data', checked when
    read) and the class tables and scores the parts reproduce or reread ('results', recorded)."""
    raise NotImplementedError


def recordings() -> dict:
    """The recordings read (names only, no value): per megamat stage the flies, criterion and E136's halves; the
    held-out panels' flies as E134 reads them; the same-fly populations; the exclusions."""
    raise NotImplementedError


def lock_contents() -> dict:
    """The lock: inputs, recordings, splits, seeds, settings, classes, outcomes, narrow readings and holds. No statistic
    is read."""
    raise NotImplementedError


def part_lock() -> dict:
    """results/.../lock.json and its committed copy experiments/E143_lock.json (commit it before any other part);
    refuses to overwrite a committed lock."""
    raise NotImplementedError


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock (tracked, no uncommitted change, the same bytes as the shared copy); PermissionError
    otherwise."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the read (each run once against all flies)

def setting_runs(setting: Setting, panel: str = PANEL) -> Dict[str, object]:
    """A setting's saved runs on one input, per T wiring (missing wirings left out)."""
    raise NotImplementedError


def read_setting(setting: Setting, panel: str = PANEL) -> Optional[pd.DataFrame]:
    """One setting's runs on one input read once (E134's default readings against all flies, E136's seed index), each
    wiring's parts (every view) against all the flies of each stage or panel read: long table (setting, input, key
    'stage/panel/condition', wiring, part, entry, value). Saved to outputs/.../views/<setting>__<input>.parquet."""
    raise NotImplementedError


def part_read(jobs: int = 1, settings: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Every setting on every input it was read on (`read_setting`), skipping what is saved; a summary per file."""
    raise NotImplementedError


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


def load_side(key: str, parts: Optional[Sequence[str]] = None) -> Side:
    """A stage or panel's side from the saved reads (`read_setting`) and the flies' parts (E134's fly_side_read)."""
    raise NotImplementedError


def side_errors(side: Side, weights: np.ndarray, views: Optional[Sequence[str]] = ('strength',),
                error_col: str = 'error') -> np.ndarray:
    """Every model row's error against many references at once: `weights` (references x flies) gives each reference's
    mean as the weighted mean of the flies' parts (a half: 1 / k on its flies; a resample: counts / n), in all the
    flies' units. Returns references x rows."""
    raise NotImplementedError


def setting_values(errors: np.ndarray, side: Side) -> pd.DataFrame:
    """Per reference, each setting's value (the mean over its T wirings) from `side_errors`' output."""
    raise NotImplementedError


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
