#!/usr/bin/env python3
"""E103: the flies' APL across odors in the held-out vial series. The readings of the flies' APL disagree about ms: in
the 12 recordings of July to September 2021 APL responds to ms at 10^-3 about as strongly as to paa, in the 2022 flies
paa evokes two to three times as much as ms; Yang's series of 2021-10/11 (30 processed recordings of 5 flies, ms, paa
and pa at three concentrations each, with two further vials of pa or of paa as her same-odor control) is read by rules
fixed before any of it is opened: which side those flies are on, the family's APL against them, what two vials of one
odor differ by in strength and in the map, and the plane. A diagnostic: it moves no loss; no simulation.

The entry is item E103 of experiment_priorities.md (proposed by viki 2026-11-15 22:10; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-15 21:48). The flies: viki's copy of Yang's processed files
(data/yang_apl_recordings/2021_oct_nov/ in the shared directory; only processed recordings, the rule of 2026-11-15 15:28;
SERIES lists the 30). The 12: E98's part A responses (results/E98_apl_response_strength/a1_responses.csv, from
data/yang_apl_calyx_maps/). The family: C2's (breadth exponent 2, 6.8 %, the 0.5 s pulse) on E72 part B's eight wirings,
read from E98's B1 runs (`E98.load_b1_all`) at the blank's level (carrying) and the settled level (partner (6)); no
variant is chosen or gated. No connectome table is read (rule of 2026-11-15 15:35).

Definitions (the entry's; E98's unless said):
- The reader: E98's part D reader (`E98.d_stim_table`, `E98.d_responses_from_means`, `E98.d_frame_means`), extended to
  blocks of three to five presentations (`n_first`) and runs with two paraffin-oil-only blocks (the order-matched
  baseline is the mean of the two at each position). Each stimulus's vial is read from the yaml's odor names
  (`vial_label`: 'pentyl acetate2' -> 'pa_vial2'; vial 1 carries as 'pa' / 'paa'; ms has one vial).
- A response (the carrying reading, READINGS['carrying']): the order-matched baseline (each presentation minus the
  paraffin-oil-only block's presentation at the same place, the mean of the two such blocks where a run has two), the
  mean of the block's first two odor presentations (N_FIRST), a stimulus in two blocks read at its first. A run without
  a paraffin-oil-only block has no carrying reading. The floor (MIN_RESPONSE 0.02 dF/F), the interpolation on log10
  concentration inside a recording's imaged range and the one-log rule are E98's (`E98.value_at`).
- A fly is a date and a fly number; a ratio is the geometric mean over flies (a fly's recordings averaged on the log
  scale) with a 90 % interval from 10,000 bootstrap draws of flies (`E98.fly_ratio`). The plane: 18 um carries
  (CARRYING_PLANE); vial 1 carries paa and pa in every ratio between odors.
- The fixed list (RATIOS): P1, P2, P3, Q1, Q2, F1, S1, S2; P1, Q1, Q2 and F1 are E98's V1 at 10^-4, its probe, its V5 at
  10^-4 and its V3 at 10^-4 (E98_IDS).
- The split rule (`split_rule`, E98's `pulse_permutation`): two sets of flies differ on a ratio if fewer than 1 in 10
  (SPLIT_SHARE) of the splits of all their flies into sets of the same sizes (1287 for 5 and 8) give a difference of
  mean log ratios at least as large (absolute).
- The family against an interval, per wiring (`class_ab`): inside if any ratio of its set over the indicator family
  (E98's INDICATOR_NS x FORMS) lies in the interval; above if all lie above its upper end; below if all below its lower
  end. A comparison takes a class when at least 6 of the 8 wirings do (`comparison_class`).
- Held (`held_class`): a class counts only if the carrying reading and every partner give it; otherwise it is
  reading-dependent. Partners: (1) E79's first-in-block baseline, which also admits the four runs without a
  paraffin-oil-only block; (2) every presentation of the block and both passes; (3) without the five marked runs
  (MARKED); (4) each of the five flies left out in turn; for the family's classes also (5) the interval from Student's
  t on the flies' log ratios (`t_interval`) and (6) the model's settled pre-odor level.
- X2 / X3 over the five estimates of ms's input (`x_over_estimates`): inside if inside is held under at least one; above
  (below) only if above (below) is held under all five; otherwise not decided.

Parts (CLI `--part`):
- 0, the identities (p0_identity.json; part A refuses to run unless all pass): the extended reader reproduces E98's part
  A responses of the 12 (a1_responses.csv) and its part D responses (d_responses.csv) exactly; the eligibility table
  rebuilt from the series' stimulus files (yaml only, no map opened) equals the entry's lists, and recordings.csv's
  frame count of every map equals its presentations (the maps themselves are checked again when part B opens them);
  each run's three-vial odor in its stimulus file agrees with the order of Yang's date notes ("change to Phe", "change
  to PA"); the family's ratio sets equal E98's c_wirings.csv, and the 12's P1 and Q1 equal E98's a3_ratios.csv.
- A, before any map of the series is opened: the 12 at the series' stimuli (the eight ratios under the carrying reading
  and partners (1) and (3), with their between-fly SDs; for the 12, partner (3) is E98's 'without_three' subset, the
  12's own marked recordings), the least difference that 5 flies against the 12's flies can show by the split rule per
  ratio (`least_detectable`), the 12's classes for X2 and X3 as descriptions (E98's verdict stands as written). Its
  tables are hashed into PART_A_LOCK (experiments/E103_part_a_lock.json), which is committed before part B runs:
  `require_part_a` refuses to open any map or response file of the series until the lock is committed and the tables
  match it.
- B, X1 (which side the held-out flies are on), C, the family against the held-out flies (X2, X3), D, the vials (Y1, M),
  E, beside (Y2, Y3, Z, Hallem's input): not implemented yet (stubs); their rules' pure functions are here and tested
  (`x1_side`, `y1_variance_split`, `vial_yardstick`, `m_share`, `m_class`, `z_rule`).

Every results table starts with Yang's caveat line ('# ' QUALITY_CAVEAT) and a line on what stays unmeasured
(UNMEASURED); `read_table` skips them. Outputs: results/E103_apl_across_odors_heldout/ (no model outputs).

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u -m experiments.E103_apl_across_odors_heldout --part 0
    nice -n 10 ./venv/bin/python -u -m experiments.E103_apl_across_odors_heldout --part A
"""
from __future__ import annotations

import argparse
import json
import socket
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E103_apl_across_odors_heldout'
Stim = Tuple[str, float]

# ---------------------------------------------------------------------------
# the series (the entry's lists; fixed before any map is opened)

SERIES_SET = 'B2021'
SERIES_FOLDER = ('data', 'yang_apl_recordings', '2021_oct_nov')
PLANES: Dict[int, Tuple[str, ...]] = {
    18: ('20211022_1FLfun3', '20211022_1FLfun5', '20211022_1FRfun3', '20211022_1FRfun5', '20211029_3FLfun2',
         '20211029_3FLfun6', '20211104_3FLfun2', '20211104_3FRfun2', '20211107_1FLfun2', '20211107_1FLfun3',
         '20211107_1FRfun2', '20211107_1FRfun4', '20211109_5FRfun1'),
    15: ('20211022_1FLfun2', '20211022_1FLfun4', '20211022_1FRfun2', '20211022_1FRfun4', '20211022_1FRfun6',
         '20211029_3FLfun1', '20211029_3FLfun3', '20211029_3FLfun5', '20211029_3FLfun7', '20211104_3FLfun1',
         '20211104_3FRfun1', '20211107_1FLfun1', '20211107_1FRfun1', '20211107_1FRfun3', '20211109_5FRfun3'),
    21: ('20211022_1FRfun7', '20211104_3FRfun3')}
SERIES: Tuple[str, ...] = tuple(sorted(x for v in PLANES.values() for x in v))
NO_PFO_BLOCK = ('20211029_3FLfun1', '20211029_3FLfun2', '20211029_3FLfun3', '20211104_3FRfun3')
THREE_VIALS_PAA = ('20211022_1FLfun4', '20211022_1FLfun5', '20211022_1FRfun5', '20211022_1FRfun6', '20211022_1FRfun7',
                   '20211029_3FLfun3', '20211029_3FLfun5', '20211029_3FLfun6', '20211104_3FRfun1', '20211104_3FRfun2',
                   '20211104_3FRfun3', '20211107_1FLfun3', '20211107_1FRfun3', '20211107_1FRfun4')   # pa in the other 16
TWO_PASSES = ('20211029_3FLfun5', '20211029_3FLfun6', '20211029_3FLfun7', '20211104_3FLfun1', '20211104_3FLfun2',
              '20211104_3FRfun1', '20211104_3FRfun2', '20211107_1FLfun1', '20211107_1FLfun2', '20211109_5FRfun1',
              '20211109_5FRfun3')
MARKED: Dict[str, str] = {   # by Yang's date notes
    '20211107_1FLfun1': 'a z drift of 6 um or more', '20211107_1FLfun2': 'a z drift of 6 um or more',
    '20211107_1FRfun2': 'frame averaging not at 3', '20211029_3FLfun3': 'the air tank found empty after the run',
    '20211022_1FRfun7': 'ended early; its name was used for two runs'}
EXTRA_OUT: Dict[str, Tuple[int, ...]] = {'20211022_1FLfun2': (15,), '20211022_1FLfun5': (15,)}   # blocks from 1: a typo
CARRYING_PLANE = 18
SERIES_FLIES = ('2021-10-22_fly1', '2021-10-29_fly3', '2021-11-04_fly3', '2021-11-07_fly1', '2021-11-09_fly5')
# the stimuli of every run (log10 concentrations): ms, paa vial 1, pa vial 1 (by date), vials 2 and 3 of pa or paa
MS_CONCS = (-4.386, -3.431, -2.0)
PAA_CONCS = (-4.386, -3.908, -3.431)
PA_CONCS_EARLY = (-5.34, -4.863, -4.386)      # 10-22 and 10-29
PA_CONCS_LATE = (-5.495, -4.796, -4.097)      # from 11-04
PRESENTATIONS_PER_BLOCK = {3: 2, 4: 3, 5: 4}   # nstim -> odor presentations after the block's paraffin oil

# ---------------------------------------------------------------------------
# the readings and the fixed list

N_FIRST = 2                                     # the carrying reading: the block's first two odor presentations
READINGS: Dict[str, dict] = {                   # the flies' readings of the series (part B on)
    'carrying': dict(baseline='order_matched', n_first=N_FIRST, which='first', without_marked=False),
    'first_in_block': dict(baseline='first_in_block', n_first=N_FIRST, which='first', without_marked=False),    # (1)
    'all_presentations': dict(baseline='order_matched', n_first=None, which='all', without_marked=False),      # (2)
    'without_marked': dict(baseline='order_matched', n_first=N_FIRST, which='first', without_marked=True)}     # (3)
# the 12's readings (part A; E98's a1_responses.csv): name -> ((baseline, blank subtracted), E98 subset)
READINGS_12: Dict[str, Tuple[Tuple[str, bool], str]] = {
    'carrying': (('order_matched', False), 'all'),
    'first_in_block': (('first_in_block', False), 'all'),                           # partner (1)
    'without_marked': (('order_matched', False), 'without_large_blank_and_3FRfun1')}  # partner (3): E98's 'without_three'
RATIOS: Tuple[Tuple[str, Stim, Stim, str], ...] = (
    ('P1', ('paa', -4.0), ('ms', -3.0), 'paa over ms'),
    ('P2', ('paa', -3.431), ('ms', -3.0), 'paa over ms'),
    ('P3', ('paa', -3.431), ('ms', -3.431), 'paa over ms'),
    ('Q1', ('pa', -4.386), ('ms', -3.0), 'pa over ms'),
    ('Q2', ('pa', -4.0), ('ms', -4.0), 'pa over ms'),
    ('F1', ('ms', -4.0), ('ms', -2.0), "ms's fall"),
    ('S1', ('paa', -4.386), ('paa', -3.431), "paa's fall"),
    ('S2', ('pa', -5.34), ('pa', -4.386), "pa's fall"))
RATIO_IDS = tuple(r[0] for r in RATIOS)
X1_RATIOS = ('P1', 'P2', 'P3', 'Q1', 'Q2', 'F1')          # the six with ms
E98_IDS = {'P1': 'V1_paa-4', 'Q1': 'probe_pa-4.386', 'Q2': 'V5_-4', 'F1': 'V3_ms-4'}
# set S (the 2022 single-colour flies; the entry's table, quotients of E100's ratios): ratio -> (lo, hi) of its value
SET_S: Dict[str, Tuple[float, float]] = {'P1': (2.2, 2.2), 'P2': (2.2, 3.2), 'P3': (3.2, 3.2), 'Q1': (2.6, 2.6),
                                         'Q2': (6.0, 6.0), 'F1': (0.17, 0.17), 'S1': (0.69, 0.69), 'S2': (0.68, 0.68)}
SPLIT_SHARE = 0.1                     # E98's PULSE_SHARE
N_SERIES_FLIES = 5                    # the held-out flies under the carrying reading at 18 um
MIN_WIRINGS = 6                       # of the 8
X_CLASSES = ('inside', 'above', 'below')
CI_LEVEL = 0.90
PARTNERS_12 = ('first_in_block', 'without_marked', 'leave_one_fly_out', 't_interval', 'settled')
X1_MAX_DIFFER_12 = 1                  # the 12's side: at most one of the six differs
X1_MIN_DIFFER_S = 4                   # the 2022 flies' side: at least four, each toward set S and nearer to it

# parts D and E (rules fixed by the entry)
Y1_YARDSTICK = 2.8                    # a difference between series beyond the vials: |log2| > 2.8 x the persistent SD
M_BLUR_PX = 5.0                       # carrying; 2, 10 and 20 px beside
M_BLURS_PX = (5.0, 2.0, 10.0, 20.0)
M_MAX_BLOCKS_APART = 4
M_MIN_R_REP = 0.3
M_ODOR_LO = 0.5                       # the odor's: the share's 90 % interval over flies above this ...
M_MIN_CALYCES = 7                     # ... and r_vial > r_odor in at least 7 of the 8 calyces
M_BLOCK_HI = 0.25                     # the block's: the interval below this
Z_MIN_HIGHER = 7                      # Z replicates: at least 7 of the 8 calyces higher at 15 um, the interval excludes 1
Z_MAX_HIGHER_NOT = 4                  # does not: 4 or fewer

QUALITY_CAVEAT = ("Yang has not checked these recordings' quality (her answer to question 19, 2026-11-15 18:08: all data "
                  "can be used, but she is unsure of quality questions); the 2021-10/11 series' quality is unchecked")
UNMEASURED = ("unmeasured: dF/F of one calyx plane against the model's single APL variable; nominal dilutions on another "
              "olfactometer than the ORN data's")
PART_A_TABLES = ('a_ratios_12.csv', 'a_per_fly_12.csv', 'a_least_detectable.csv', 'a_classes_12.csv',
                 'a_class_wirings_12.csv', 'a_summary.json')
PART_A_LOCK = Path(__file__).resolve().parent / 'E103_part_a_lock.json'


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def write_table(df: pd.DataFrame, path: Path) -> None:
    """A results table: its first line Yang's caveat ('# ' QUALITY_CAVEAT), its second what stays unmeasured, then the
    csv (`read_table` reads it; report_values skips leading '#' lines)."""
    raise NotImplementedError


def read_table(path: Path) -> pd.DataFrame:
    """A table written by `write_table` (its leading '#' lines skipped)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the guard: no map or response file of the series before part A's lock

def part_a_locked(results: Optional[Path] = None, lock: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether part A's tables exist in `results` (default `out_dir()`), equal the sha256 recorded in `lock` (default
    PART_A_LOCK), and the lock file is committed in this repository (tracked, no uncommitted change; `check_git`)."""
    raise NotImplementedError


def require_part_a(stem: str, **kw) -> None:
    """Raises PermissionError when `stem` is a recording of the series (set B2021) and `part_a_locked` is False."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the recordings and the reader

def series_rows(stems: Sequence[str]) -> pd.DataFrame:
    """viki's recordings.csv rows of `stems` (the series' 30, the 12 and E98's part D lists), indexed by stem, with
    'fly_id', 'side' and 'role' ('series', '12', 'D1', 'D2', 'D3'). Raises PermissionError on any other stem and on a
    row that is not processed (the rule of 2026-11-15 15:28)."""
    raise NotImplementedError


def vial_label(odor: str) -> Tuple[str, int]:
    """A yaml odor as abbreviated by E100 ('pentyl acetate2', 'pa', 'phenylacetaldehyde1', 'ms') -> (label, vial): vial
    1 as the odor ('pa', 1), vials 2 and 3 as '<odor>_vial<k>'; an odor without a vial digit is vial 1."""
    raise NotImplementedError


def stim_table(row) -> pd.DataFrame:
    """One recording's presentations from its yaml (`E98.d_stim_table`, blocks from 1, EXTRA_OUT left out), the odors
    relabelled by `vial_label` (column 'vial')."""
    raise NotImplementedError


def frame_means(row, source: str = 'maps') -> np.ndarray:
    """Per-presentation calyx means (`E98.d_frame_means`: 'maps' or Yang's 'yang' table); `require_part_a` first."""
    raise NotImplementedError


def responses_from_means(means, table: pd.DataFrame, baseline: str, which: str = 'first',
                         n_first: Optional[int] = None) -> pd.DataFrame:
    """`E98.d_responses_from_means` with `n_first` (the block's first `n_first` odor presentations; None: all)."""
    raise NotImplementedError


def read_responses(stems: Sequence[str], source: str = 'maps', baselines: Sequence[str] = ('order_matched', 'first_in_block'),
                   which: str = 'first', n_first: Optional[int] = None) -> pd.DataFrame:
    """Responses in E98's `stimulus_responses` format for `stems`: per recording, stimulus and baseline; a baseline a
    run cannot give (no paraffin-oil-only block for 'order_matched') is left out for that run."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part 0: the identities

def eligibility_table() -> pd.DataFrame:
    """The series' eligibility from its stimulus files (yaml) and recordings.csv alone (no map opened): per recording
    its fly, side, plane, nstim and presentations per block, its paraffin-oil-only blocks, the three-vial odor, two
    passes, the blocks left out, marked, the map's frame count (recordings.csv) against the presentations, the stimuli
    against the entry's, beside the entry's lists and whether they match."""
    raise NotImplementedError


def note_vial_order(path: Path) -> Dict[str, str]:
    """From one date note: each run header's three-vial odor ('pa' or 'paa') by the order of the note's "change to
    Phe" / "change to PA" lines (before the first change, the other odor of its target, or the note's own pin table)."""
    raise NotImplementedError


def vials_identity(el: pd.DataFrame) -> pd.DataFrame:
    """Per run of the series, the three-vial odor of its stimulus file against the date note's order."""
    raise NotImplementedError


def reader_identity() -> dict:
    """The extended reader against E98's a1_responses.csv (the 12) and d_responses.csv (part D's recordings)."""
    raise NotImplementedError


def family_identity() -> dict:
    """The family's ratio sets (`stimulus_classes` on B1's rows at the blank's and the settled levels) against E98's
    c_wirings.csv (V1's comparisons, carrying and settled)."""
    raise NotImplementedError


def twelve_identity(ratios_12: Optional[pd.DataFrame] = None) -> dict:
    """The 12's P1 and Q1 (and, beside, Q2 and F1) under the carrying reading against E98's a3_ratios.csv."""
    raise NotImplementedError


def part_0() -> dict:
    """Part 0: writes p0_eligibility.csv, p0_vials.csv, p0_reader_identity.csv, p0_family_identity.csv and
    p0_identity.json; raises if an identity fails."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# ratios, intervals and classes

def flies_ratios_12(reading: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The 12's eight ratios under one of READINGS_12 (`E98.flies_ratios` on a1_responses.csv): (rows, per recording)."""
    raise NotImplementedError


def t_interval(per_fly_log: Sequence[float], level: float = CI_LEVEL) -> Tuple[float, float, float]:
    """The geometric mean over flies and its `level` interval from Student's t on the flies' (natural) log ratios."""
    raise NotImplementedError


def split_rule(set_a: pd.Series, set_b: pd.Series) -> dict:
    """E98's `pulse_permutation` on two sets of flies' mean log ratios: 'differ' when fewer than SPLIT_SHARE of the
    splits of the pooled flies into sets of the same sizes give a difference at least as large."""
    raise NotImplementedError


def least_detectable(per_fly_log2: Sequence[float], n_new: int = N_SERIES_FLIES, sd_log2: Optional[float] = None) -> dict:
    """The least difference (log2, and as a factor, up and down) that `n_new` flies spread like the given flies (normal
    scores scaled to `sd_log2`, default their between-fly SD) can show against them by the split rule; beside, the same
    with both sets normal scores and the normal approximation 1.645 sd sqrt(1/n + 1/n_new)."""
    raise NotImplementedError


def class_ab(ratios: Sequence[float], lo: float, hi: float) -> str:
    """One wiring's set against [lo, hi]: 'inside', 'above', 'below' or 'mixed'."""
    raise NotImplementedError


def comparison_class(classes: Sequence[str], min_wirings: int = MIN_WIRINGS) -> str:
    """'inside', 'above' or 'below' if at least `min_wirings` wirings are; otherwise 'not decided'."""
    raise NotImplementedError


def held_class(categories: Dict[str, str], partners: Sequence[str]) -> str:
    """The carrying reading's class when it and every partner give it; otherwise 'reading-dependent'."""
    raise NotImplementedError


def x_over_estimates(held: Dict[str, str]) -> str:
    """X2 / X3 over the estimates of ms's input: 'inside' if held inside under any; 'above' / 'below' if held so under
    every one; otherwise 'not decided'."""
    raise NotImplementedError


def wiring_classes(on: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring the model's ratio set (`E98.stimulus_classes`' rows) with `class_ab` against [lo, hi]."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A

def part_a() -> dict:
    """Part A (no map of the series opened): the 12 at the series' stimuli, the least detectable differences, the 12's
    classes for X2 and X3 (descriptions); writes PART_A_TABLES and the lock."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the rules of parts B to E (pure functions)

def x1_side(rows: pd.DataFrame) -> str:
    """X1 from the six ratios with ms (columns ratio_id, differ, toward_s, nearer_s): "the 12's side" (at most one
    differs), "the 2022 flies' side" (at least four differ, each toward set S and nearer to set S's value on the log
    scale than to the 12's), else 'between'."""
    raise NotImplementedError


def y1_variance_split(rows: pd.DataFrame) -> dict:
    """Y1: from per-recording log2 ratios of two vials of one odor (columns day, recording, odor, conc, pair,
    log2_ratio): the SD over single recordings and the persistent SD (per day, odor, concentration and pair of vials the
    mean over that day's recordings; the SD of those means with the scatter of single recordings about them removed,
    floored at 0)."""
    raise NotImplementedError


def vial_yardstick(persistent_sd: float) -> float:
    """The log2 difference between two series beyond the vials: Y1_YARDSTICK x the persistent SD."""
    raise NotImplementedError


def m_share(r_vial: float, r_odor: float, r_rep: float) -> float:
    """M's share per calyx: (r_vial - r_odor) / (r_rep - r_odor)."""
    raise NotImplementedError


def m_class(share_lo: float, share_hi: float, n_vial_over_odor: int) -> str:
    """M: "the odor's" (the share's interval above M_ODOR_LO and r_vial > r_odor in at least M_MIN_CALYCES calyces),
    "the block's" (the interval below M_BLOCK_HI), else 'not decided'."""
    raise NotImplementedError


def z_rule(log_diff_per_calyx: Sequence[float], lo: float, hi: float) -> str:
    """Z: 'replicates' if at least Z_MIN_HIGHER calyces are higher at 15 um and the mean factor's interval [lo, hi]
    excludes 1; 'does not' if Z_MAX_HIGHER_NOT or fewer are higher; else 'not decided'."""
    raise NotImplementedError


def part_b() -> dict:
    raise NotImplementedError('part B opens the series\' maps: not implemented yet')


def part_c() -> dict:
    raise NotImplementedError


def part_d() -> dict:
    raise NotImplementedError


def part_e() -> dict:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'E'))
    a = ap.parse_args()
    {'0': part_0, 'A': part_a, 'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e}[a.part]()


if __name__ == '__main__':
    main()
