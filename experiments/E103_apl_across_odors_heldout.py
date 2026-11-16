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
variant is chosen or gated. No connectome table is read (rule of 2026-11-15 15:35): E79's data side only, through
E100's `_ds`.

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
- The split rule (`split_rule`, E98's `pulse_permutation`; `split_share` vectorized): two sets of flies differ on a
  ratio if fewer than 1 in 10 (SPLIT_SHARE) of the splits of all their flies into sets of the same sizes (1287 for 5 and
  8) give a difference of mean log ratios at least as large (absolute).
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
  frame count of every map equals its presentations (the maps themselves are checked again when a later part opens
  them, `read_responses`); each run's three-vial odor in its stimulus file agrees with the order of Yang's date notes
  ("change to Phe", "change to PA"); the family's ratio sets equal E98's c_wirings.csv, and the 12's P1 and Q1 equal
  E98's a3_ratios.csv.
- A, before any map of the series is opened: the 12 at the series' stimuli (the eight ratios under the carrying reading
  and partners (1) and (3), with their between-fly SDs; for the 12, partner (3) is E98's 'without_three' subset, the
  12's own marked recordings: the two with a large blank and 20210805_3FRfun1), the least difference that 5 flies
  against the 12's flies can show by the split rule per ratio (`least_detectable`), the 12's classes for X2 and X3 as
  descriptions (E98's verdict stands as written). Its tables are hashed into PART_A_LOCK
  (experiments/E103_part_a_lock.json), which is committed before any later part runs: `require_part_a` refuses to open
  any map or response file of the series until the lock is committed and the tables match it; a re-run of part A once
  the lock exists writes to a_recheck/ and never over the locked tables.
- B, X1 (which side the held-out flies are on), C, the family against the held-out flies (X2, X3), D, the vials (Y1, M),
  E, beside (Y2, Y3, Z, Hallem's input): not implemented yet (stubs); their rules' pure functions are here and tested
  (`x1_direction`, `x1_side`, `y1_variance_split`, `vial_yardstick`, `m_calyx`, `m_share`, `m_class`, `z_rule`).

Deviations and readings of the entry (each also in the tables or a_summary.json):
- Part 0's frame-count identity is read from recordings.csv (`peak_map_shape`, `frames_equal_presentations`, made when
  the copy was checked), so that no map of the series is opened before part A; `read_responses` checks every map's
  frame count again when it opens it.
- Partner (3) for the 12 (part A) is E98's own quality subset (`without_large_blank_and_3FRfun1`); partner (2) is the
  carrying reading itself for the 12 (two presentations per block, one pass), and is not listed again.
- Set S's values (SET_S) are the entry's table's estimates, fixed before the series is read; P2's is the range 2.2 to
  3.2, and "nearer to set S's value" is the log distance to that range.

Every results table starts with Yang's caveat line ('# ' QUALITY_CAVEAT) and a line on what stays unmeasured
(UNMEASURED); `read_table` skips them. Outputs: results/E103_apl_across_odors_heldout/ (no model outputs).

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u -m experiments.E103_apl_across_odors_heldout --part 0
    nice -n 10 ./venv/bin/python -u -m experiments.E103_apl_across_odors_heldout --part A
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import socket
import subprocess
import time
from itertools import combinations
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
SERIES_DATES = ('20211022', '20211029', '20211104', '20211105', '20211107', '20211109')   # every run of the set
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
N_CARRYING_RECORDINGS = 12                      # 18 um with a paraffin-oil-only block
# the stimuli of every run (log10 concentrations): ms, paa vial 1, pa vial 1 (by date), vials 2 and 3 of pa or paa
MS_CONCS = (-4.386, -3.431, -2.0)
PAA_CONCS = (-4.386, -3.908, -3.431)
PA_CONCS_EARLY = (-5.34, -4.863, -4.386)      # 10-22 and 10-29
PA_CONCS_LATE = (-5.495, -4.796, -4.097)      # from 11-04
PA_LATE_FROM = '20211104'
PRESENTATIONS_PER_BLOCK = {3: 2, 4: 3, 5: 4}   # nstim -> odor presentations after the block's paraffin oil
N_RUNS_BY_PRESENTATIONS = {2: 10, 3: 4, 4: 16}  # the entry: two times (10-22), three times (4 runs), four times (16)
PASS_DISTANCE = 16                              # two passes 16 blocks apart

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
IDENTITY_IDS = ('P1', 'Q1')                                # part 0's identity; Q2 and F1 beside
# the entry's table: the 12's values as the planner estimated them (ratio, lo, hi), and set S's (the 2022 single-colour
# flies, quotients of E100's ratios; lo, hi of its value, a point where lo == hi)
ENTRY_12: Dict[str, Tuple[float, float, float]] = {
    'P1': (0.81, 0.67, 0.96), 'P2': (0.89, 0.68, 1.12), 'P3': (1.05, 0.92, 1.19), 'Q1': (1.12, 1.01, 1.24),
    'Q2': (1.47, 1.32, 1.64), 'F1': (0.60, 0.49, 0.74), 'S1': (0.76, 0.64, 0.91), 'S2': (0.80, 0.71, 0.89)}
SET_S: Dict[str, Tuple[float, float]] = {'P1': (2.2, 2.2), 'P2': (2.2, 3.2), 'P3': (3.2, 3.2), 'Q1': (2.6, 2.6),
                                         'Q2': (6.0, 6.0), 'F1': (0.17, 0.17), 'S1': (0.69, 0.69), 'S2': (0.68, 0.68)}
LDD_ENTRY_P1 = (1.4, 0.46)            # the entry's estimate for P1: a factor of about 1.4 from a between-fly SD of 0.46 (log2)
SPLIT_SHARE = 0.1                     # E98's PULSE_SHARE
N_SERIES_FLIES = 5                    # the held-out flies under the carrying reading at 18 um
MIN_WIRINGS = 6                       # of the 8
X_CLASSES = ('inside', 'above', 'below')
CI_LEVEL = 0.90
NORMAL_Z90 = 1.6449                   # two-sided 90 %
LDD_TOL = 1e-7                        # the least detectable difference's bisection (log2)
X_IDS = {'X2': ('P1', 'P2'), 'X3': ('P1',)}   # X2: the hull of P1's and P2's intervals; X3: P1's alone
X1_MAX_DIFFER_12 = 1                  # the 12's side: at most one of the six differs
X1_MIN_DIFFER_S = 4                   # the 2022 flies' side: at least four, each toward set S and nearer to it
SIDES = ("the 12's side", "the 2022 flies' side", 'between')

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


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')


def _json_default(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return str(x)


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_table(df: pd.DataFrame, path: Path) -> None:
    """A results table: its first line Yang's caveat ('# ' QUALITY_CAVEAT), its second what stays unmeasured, then the
    csv (`read_table` reads it; report_values skips leading '#' lines)."""
    with open(path, 'w') as f:
        f.write(f'# {QUALITY_CAVEAT}\n# {UNMEASURED}\n')
        df.to_csv(f, index=False)


def read_table(path: Path) -> pd.DataFrame:
    """A table written by `write_table` (its leading '#' lines skipped)."""
    n = 0
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line.startswith('#'):
                break
            n += 1
    return pd.read_csv(path, skiprows=n, float_precision='round_trip')


def _e98():
    """E98's module, with E79's data paths repointed to this host's shared directory first (E100's `_ds`, which does
    not import E79_connectome_maps)."""
    from experiments import E100_apl_weak_odors as E100
    E100._ds()
    from experiments import E98_apl_response_strength as E98
    return E98


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, experiment=EXPERIMENT, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                                         **extra))
    _write_json(path, d)


# ---------------------------------------------------------------------------
# the guard: no map or response file of the series before part A's lock

def is_series(stem: str) -> bool:
    """Whether `stem` is a run of the 2021-10/11 series (any of its dates, processed or raw only)."""
    return str(stem)[:8] in SERIES_DATES


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by its git repository and has no uncommitted change."""
    repo = Path(path).resolve().parent
    try:
        r = subprocess.run(['git', '-C', str(repo), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(repo), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def part_a_locked(results: Optional[Path] = None, lock: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether part A's tables exist in `results` (default `out_dir()`), equal the sha256 recorded in `lock` (default
    PART_A_LOCK), and the lock file is committed in this repository (tracked, no uncommitted change; `check_git`)."""
    results = out_dir() if results is None else Path(results)
    lock = PART_A_LOCK if lock is None else Path(lock)
    if not lock.exists():
        return False
    try:
        tabs = json.loads(lock.read_text()).get('tables') or {}
    except ValueError:
        return False
    if set(tabs) != set(PART_A_TABLES):
        return False
    for name, sha in tabs.items():
        p = results / name
        if not p.exists() or _sha256(p) != sha:
            return False
    return _committed(lock) if check_git else True


def require_part_a(stem: str, **kw) -> None:
    """Raises PermissionError when `stem` is a run of the series and `part_a_locked(**kw)` is False: no map or response
    file of the series is opened before part A's tables are written and their lock committed."""
    if is_series(stem) and not part_a_locked(**kw):
        raise PermissionError(f"{stem}: no map or response file of the 2021-10/11 series is opened before part A's tables "
                              f"are written and their lock ({PART_A_LOCK.name}) committed")


# ---------------------------------------------------------------------------
# the recordings and the reader

def check_row(row) -> None:
    """Raises PermissionError on a recordings.csv row that is not processed (raw only: the rule of 2026-11-15 15:28)."""
    v = row.get('use_until_further_notice', '') if hasattr(row, 'get') else ''
    if not str(v).startswith('processed'):
        raise PermissionError(f"{row.get('stem')}: not processed ({v}); raw-only recordings are not read until further "
                              "notice (rule of 2026-11-15 15:28)")


def allowed_stems() -> Dict[str, str]:
    """The recordings this module reads, stem -> role: the series' 30 ('series'), the 12 ('12') and E98's part D lists
    ('D1', 'D2', 'D3'; part 0's reader identity only)."""
    E98 = _e98()
    out = {st: 'series' for st in SERIES}
    out.update({st: '12' for st in E98.twelve()})
    for st in E98.D_RECORDINGS:
        out.setdefault(st, 'D1' if st in E98.D1_SHUFFLED else 'D2' if st in E98.D2_ASCENDING else 'D3')
    return out


def series_rows(stems: Sequence[str]) -> pd.DataFrame:
    """viki's recordings.csv rows of `stems` (the series' 30, the 12 and E98's part D lists), indexed by stem, with
    'fly_id', 'side' and 'role' ('series', '12', 'D1', 'D2', 'D3'). Raises PermissionError on any other stem and on a
    row that is not processed (the rule of 2026-11-15 15:28)."""
    from experiments import E100_apl_weak_odors as E100
    E98 = _e98()
    roles = allowed_stems()
    bad = [x for x in stems if x not in roles]
    if bad:
        raise PermissionError(f'E103 reads only the series\' 30 processed runs, the 12 and E98\'s part D lists: {bad}')
    r = E100.read_recordings(sets=None)
    out = []
    for st in stems:
        g = r[(r['stem'] == st) & (r['date'].astype(str).str.replace('-', '') == st[:8])]
        if len(g) != 1:
            raise RuntimeError(f'{st}: {len(g)} rows in recordings.csv')
        row = g.iloc[0].to_dict()
        check_row(row)
        row['fly_id'] = E98.fly_id(row['date'], row['fly'])
        row['side'] = str(row['hemisphere'])
        row['role'] = roles[st]
        out.append(row)
    return pd.DataFrame(out).set_index('stem', drop=False)


def vial_label(odor: str) -> Tuple[str, int]:
    """A yaml odor as abbreviated by E100 ('pentyl acetate2', 'pa', 'phenylacetaldehyde1', 'ms') -> (label, vial): vial
    1 as the odor ('pa', 1), vials 2 and 3 as '<odor>_vial<k>'; an odor without a vial digit is vial 1."""
    from experiments import E100_apl_weak_odors as E100
    o = ' '.join(str(odor).split())
    m = re.fullmatch(r'(?P<b>.*[a-z])\s*(?P<v>[1-9])', o)
    if m is None or m.group('b').lower() not in E100.ODOR_ABBREV:    # only a full odor name carries a vial digit
        return E100.ODOR_ABBREV.get(o.lower(), o), 1
    ab, vial = E100.ODOR_ABBREV[m.group('b').lower()], int(m.group('v'))
    return (ab, 1) if vial == 1 else (f'{ab}_vial{vial}', vial)


def table_from_yaml(path: Path, stem: str, extra_out: Sequence[int] = ()) -> pd.DataFrame:
    """One recording's presentations from its yaml (`E98.d_stim_table`, blocks from 1, `extra_out` left out), the odors
    relabelled by `vial_label` ('odor', 'vial'; the solvent vial 0)."""
    E98 = _e98()
    t = E98.d_stim_table(path, stem, extra_out)
    lab = [(o, 0) if o == E98.SOLVENT else vial_label(o) for o in t['odor']]
    t['odor'] = [x for x, _ in lab]
    t['vial'] = [v for _, v in lab]
    return t


def stim_table(row) -> pd.DataFrame:
    """One recording's presentations from its yaml (`table_from_yaml`; the series' EXTRA_OUT, part D's D_EXTRA_OUT)."""
    E98 = _e98()
    st = row['stem']
    extra = EXTRA_OUT.get(st, ()) if is_series(st) else E98.D_EXTRA_OUT.get(st, ())
    return table_from_yaml(E98.d_yaml_path(row), st, extra)


def frame_means(row, source: str = 'maps') -> np.ndarray:
    """Per-presentation calyx means (`E98.d_frame_means`: 'maps', the peak maps inside Yang's outline, or 'yang', her
    roi_peak_dff table); `require_part_a` first, and `check_row`."""
    require_part_a(row['stem'])
    check_row(row)
    return _e98().d_frame_means(row, source)


def responses_from_means(means, table: pd.DataFrame, baseline: str, which: str = 'first',
                         n_first: Optional[int] = None) -> pd.DataFrame:
    """`E98.d_responses_from_means` with `n_first` (the block's first `n_first` odor presentations; None: all)."""
    return _e98().d_responses_from_means(means, table, baseline, which, n_first)


def read_responses(stems: Sequence[str], source: str = 'maps', baselines: Sequence[str] = ('order_matched', 'first_in_block'),
                   which: str = 'first', n_first: Optional[int] = None) -> pd.DataFrame:
    """Responses in E98's `stimulus_responses` / `d_resp` format for `stems` (`series_rows`; the series' runs only once
    part A is locked): per recording, stimulus (block) and baseline, with the vial, the plane and the role. Each map's
    frame count must equal its presentations. A run without a paraffin-oil-only block has no 'order_matched' rows."""
    rows = series_rows(stems)
    out = []
    for st, row in rows.iterrows():
        tab = stim_table(row)
        m = frame_means(row, source)
        if len(m) != len(tab):
            raise RuntimeError(f'{st}: {len(m)} frames for {len(tab)} presentations')
        vial = {(o, float(c)): int(v) for o, c, v in zip(tab['odor'], tab['log10_conc'], tab['vial'])}
        for bl in baselines:
            if bl == 'order_matched' and not tab.loc[~tab['out'], 'pfo_block'].any():
                continue
            r = responses_from_means(m, tab, bl, which, n_first)
            reps = sorted((c for c in r.columns if re.fullmatch(r'rep\d+', c)), key=lambda c: int(c[3:]))
            for x in r.to_dict('records'):
                out.append(dict(recording=st, fly=row['fly_id'], side=row['side'], pulse_s=float(row['pulse_s']),
                                plane_um=row['plane_depth_um'], role=row['role'], odor=x['odor'], conc=x['conc'],
                                vial=vial.get((x['odor'], float(x['conc'])), np.nan),
                                stimulus=f"{x['odor']} {x['conc']:g}", block=x['block1'] - 1, block1=x['block1'],
                                baseline=bl, blank_subtracted=False, blank=np.nan, response=x['response'],
                                **{k: x[k] for k in reps}, n_reps=x['n_reps'], source=source))
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# part 0: the identities

def _note_path(row) -> Path:
    return Path(row['folder']) / f"{row['stem'][:8]}.txt"


def _stimuli_expected(stem: str, three: str) -> set:
    pa = PA_CONCS_LATE if stem[:8] >= PA_LATE_FROM else PA_CONCS_EARLY
    concs = {'ms': MS_CONCS, 'paa': PAA_CONCS, 'pa': pa}
    s = {(o, round(c, 3)) for o in ('ms', 'paa', 'pa') for c in concs[o]}
    s |= {(f'{three}_vial{k}', round(c, 3)) for k in (2, 3) for c in concs[three]}
    return s


def eligibility_row(stem: str, table: pd.DataFrame, row) -> dict:
    """One run's eligibility from its stimulus table (`stim_table`) and recordings.csv row (no map opened), beside the
    entry's lists ('_expected') and whether each matches. 'stimuli_as_entry' compares the stimulus file's stimuli
    (blocks left out included) with the entry's; 'not_imaged_after_blocks_out' lists those only a left-out block held
    (20211022_1FLfun2's pa vial 2 and 1FLfun5's paa vial 2 at the middle concentration: block 15's typo)."""
    E98 = _e98()
    t = table
    kept = t[~t['out']]
    nstim = int(t['position'].max())
    blocks = kept.groupby('block1')
    pfo = sorted(int(b) for b, g in blocks if bool(g['pfo_block'].iloc[0]))
    # the design: every odor block is the solvent, then one stimulus nstim - 1 times
    design_ok = True
    for b, g in blocks:
        if bool(g['pfo_block'].iloc[0]):
            continue
        g = g.sort_values('position')
        o = g[g['odor'] != E98.SOLVENT]
        if not (g['odor'].iloc[0] == E98.SOLVENT and len(o) == nstim - 1 and o['position'].min() == 2
                and o[['odor', 'log10_conc']].drop_duplicates().shape[0] == 1):
            design_ok = False
    odors = kept[kept['odor'] != E98.SOLVENT]
    vials = sorted({o.split('_vial')[0] for o in odors['odor'] if '_vial' in o})
    three = vials[0] if len(vials) == 1 else ';'.join(vials) or 'none'
    stims = {(o, round(float(c), 3)) for o, c in zip(odors['odor'], odors['log10_conc'])}
    every = t[t['odor'] != E98.SOLVENT]           # the stimulus file's stimuli, blocks left out included
    stims_all = {(o, round(float(c), 3)) for o, c in zip(every['odor'], every['log10_conc'])}
    per_stim = odors.drop_duplicates(['odor', 'log10_conc', 'block1']).groupby(['odor', 'log10_conc'])['block1']
    nblocks = per_stim.nunique()
    two_passes = bool(len(nblocks) and (nblocks == 2).all())
    dists = sorted({int(g.max() - g.min()) for _, g in per_stim if g.nunique() == 2})
    shape = str(row.get('peak_map_shape', ''))
    map_frames = int(shape.split('x')[0]) if re.match(r'^\d+x', shape) else -1
    plane = float(row['plane_depth_um'])
    plane_exp = next((float(p) for p, v in PLANES.items() if stem in v), np.nan)
    rec = dict(stem=stem, fly=row['fly_id'], side=row['side'], plane_um=plane, pulse_s=float(row['pulse_s']),
               nblock=int(t['block1'].max()), nstim=nstim, odor_presentations_per_block=nstim - 1,
               n_presentations=int(len(t)), recordings_csv_n_presentations=float(row['n_presentations']),
               map_frames_recordings_csv=map_frames,
               frames_equal_presentations_recordings_csv=row.get('frames_equal_presentations'),
               pfo_only_blocks=';'.join(str(b) for b in pfo), n_pfo_only_blocks=len(pfo),
               three_vial_odor=three, two_passes=two_passes, pass_distances=';'.join(str(d) for d in dists),
               blocks_left_out=';'.join(str(b) for b in sorted(t.loc[t['out'], 'block1'].unique())),
               marked=MARKED.get(stem, ''), design_ok=design_ok,
               stimuli_as_entry=bool(stims_all == _stimuli_expected(stem, three)) if three in ('pa', 'paa') else False,
               not_imaged_after_blocks_out=';'.join(f'{o} {c:g}' for o, c in sorted(stims_all - stims)),
               plane_expected=plane_exp, no_pfo_block_expected=stem in NO_PFO_BLOCK,
               three_vial_odor_expected='paa' if stem in THREE_VIALS_PAA else 'pa', two_passes_expected=stem in TWO_PASSES,
               blocks_left_out_expected=';'.join(str(b) for b in EXTRA_OUT.get(stem, ())),
               carrying=bool(plane == CARRYING_PLANE and len(pfo) > 0))
    rec['frames_ok'] = bool(map_frames == len(t) and rec['recordings_csv_n_presentations'] == len(t)
                            and str(rec['frames_equal_presentations_recordings_csv']) == 'True')
    rec['matches'] = bool(plane == plane_exp and (len(pfo) == 0) == rec['no_pfo_block_expected']
                          and three == rec['three_vial_odor_expected'] and two_passes == rec['two_passes_expected']
                          and rec['blocks_left_out'] == rec['blocks_left_out_expected'] and design_ok
                          and rec['stimuli_as_entry'] and rec['frames_ok']
                          and (not two_passes or dists == [PASS_DISTANCE]))
    return rec


def eligibility_table() -> pd.DataFrame:
    """The series' eligibility from its stimulus files (yaml) and recordings.csv alone (no map opened): per recording
    its fly, side, plane, nstim and presentations per block, its paraffin-oil-only blocks, the three-vial odor, two
    passes, the blocks left out, marked, the map's frame count (recordings.csv) against the presentations, the stimuli
    against the entry's, beside the entry's lists and whether they match; then the set's raw-only runs from
    recordings.csv alone ('raw only: not read'; no file of theirs is opened)."""
    from experiments import E100_apl_weak_odors as E100
    rows = series_rows(list(SERIES))
    out = [eligibility_row(st, stim_table(row), row) for st, row in rows.iterrows()]
    r = E100.read_recordings(sets=[SERIES_SET])
    for _, x in r[~r['stem'].isin(SERIES)].iterrows():
        out.append(dict(stem=x['stem'], plane_um=x['plane_depth_um'], marked='',
                        raw_only=str(x.get('use_until_further_notice', '')), matches=np.nan))
    return pd.DataFrame(out)


_HEADER = re.compile(r'^\s*(\d+F[LR]fun\d+)\b')
_CHANGE = re.compile(r'\bchange\b.*?\bto\s+(phe|pa)\b', re.I)
_PIN_TABLE = re.compile(r'^\s*(3[89]|4[0-3])\s*:?\s+(pa|phe)\s*[23]\b', re.I)


def _note_odor(word: str) -> str:
    return 'paa' if word.lower().startswith('phe') else 'pa'


def note_vial_order(path: Path) -> Dict[str, str]:
    """From one date note: each run header's three-vial odor ('pa' or 'paa') by the order of the note's "change to
    Phe" / "change to PA" lines (before the first change, the note's pin table of pins 38 to 43 where it has one before
    the first run, else the other odor of the first change's target); a run whose headers disagree is 'conflict'."""
    lines = Path(path).read_text(errors='replace').splitlines()
    first_change = next((_note_odor(m.group(1)) for ln in lines for m in [_CHANGE.search(ln)] if m), None)
    state = None
    out: Dict[str, str] = {}
    seen_header = False
    for ln in lines:
        m = _PIN_TABLE.match(ln)
        if m and not seen_header:
            state = _note_odor(m.group(2))
            continue
        m = _CHANGE.search(ln)
        if m:
            state = _note_odor(m.group(1))
            continue
        m = _HEADER.match(ln)
        if m:
            if state is None:
                state = {'pa': 'paa', 'paa': 'pa'}.get(first_change, 'unknown')
            seen_header = True
            run = m.group(1)
            out[run] = state if out.get(run, state) == state else 'conflict'
    return out


def vials_identity(el: pd.DataFrame) -> pd.DataFrame:
    """Per run of the series, the three-vial odor of its stimulus file against the date note's order."""
    rows = series_rows(list(SERIES))
    notes: Dict[str, Dict[str, str]] = {}
    out = []
    for st, row in rows.iterrows():
        p = _note_path(row)
        if str(p) not in notes:
            notes[str(p)] = note_vial_order(p) if p.exists() else {}
        note = notes[str(p)].get(st[9:], 'not in the note')
        yml = el.set_index('stem').loc[st, 'three_vial_odor']
        out.append(dict(stem=st, note=p.name, three_vial_odor_yaml=yml, three_vial_odor_note=note, agrees=bool(yml == note)))
    return pd.DataFrame(out)


IDENTITY_TOL = 1e-12


def _compare(a: pd.DataFrame, b: pd.DataFrame, key: Sequence[str], cols: Sequence[str]) -> dict:
    a = a.set_index(list(key))[list(cols)].sort_index()
    b = b.set_index(list(key))[list(cols)].sort_index()
    same = a.index.equals(b.index)
    if not same:
        return dict(index_equal=False, n_a=int(len(a)), n_b=int(len(b)), max_abs_diff=np.nan, exact=False,
                    only_a=[list(x) for x in a.index.difference(b.index)[:10]],
                    only_b=[list(x) for x in b.index.difference(a.index)[:10]])
    d = (a.astype(float) - b.astype(float)).abs()
    both_nan = a.isna() & b.isna()
    d = d.where(~both_nan, 0.0)
    mx = float(np.nanmax(d.to_numpy())) if d.size else 0.0
    nan_mismatch = int((a.isna() ^ b.isna()).to_numpy().sum())
    return dict(index_equal=True, n=int(len(a)), max_abs_diff=mx, nan_mismatch=nan_mismatch,
                exact=bool(mx <= IDENTITY_TOL and nan_mismatch == 0), zero=bool(mx == 0.0 and nan_mismatch == 0))


def reader_identity() -> Tuple[dict, pd.DataFrame]:
    """The extended reader (E98's mode: every presentation of the block, `n_first` None) against E98's a1_responses.csv
    (the 12, both baselines, nothing subtracted: response, rep1, rep2) and d_responses.csv (part D's recordings and the
    12, as E98's part D read them: response, rep1 to rep3, n_reps), within IDENTITY_TOL; and, on the 12, the carrying reading (first two presentations)
    equal to E98's (their blocks have two). Returns (result, per-recording max differences)."""
    E98 = _e98()
    d98 = E98.out_dir()
    key = ['recording', 'stimulus', 'baseline']
    a1 = pd.read_csv(d98 / 'a1_responses.csv', float_precision='round_trip')
    a1 = a1[~a1['blank_subtracted'].astype(bool) & a1['baseline'].isin(E98.D_BASELINES)]
    mine12 = read_responses(E98.twelve(), 'maps', E98.D_BASELINES, 'first', None)
    mine12_2 = read_responses(E98.twelve(), 'maps', E98.D_BASELINES, 'first', N_FIRST)
    dres = pd.read_csv(d98 / 'd_responses.csv', float_precision='round_trip')
    dres = dres[dres['source'] == 'maps']
    mineD = read_responses(list(E98.D_RECORDINGS) + E98.twelve(), 'maps', E98.D_BASELINES, 'first', None)   # E98's part D read both
    res = dict(twelve=_compare(a1, mine12, key, ('response', 'rep1', 'rep2')),
               twelve_first_two=_compare(a1, mine12_2, key, ('response', 'rep1', 'rep2')),
               part_d=_compare(dres, mineD, key, ('response', 'rep1', 'rep2', 'rep3', 'n_reps')))
    per = []
    for name, (a, b, cols) in dict(twelve=(a1, mine12, ('response', 'rep1', 'rep2')),
                                   part_d=(dres, mineD, ('response', 'rep1', 'rep2', 'rep3'))).items():
        m = a.merge(b, on=key, suffixes=('_e98', '_e103'))
        for rec, g in m.groupby('recording'):
            dd = max(float((g[f'{c}_e98'] - g[f'{c}_e103']).abs().max()) for c in cols)
            per.append(dict(set=name, recording=rec, n_rows=int(len(g)), max_abs_diff=dd))
    res['passed'] = bool(res['twelve']['exact'] and res['part_d']['exact'] and res['twelve_first_two']['exact'])
    return res, pd.DataFrame(per)


def _on(b1: pd.DataFrame) -> pd.DataFrame:
    return b1[b1['cond'] == 'on'][['wiring', 'stimulus', 'win', 'pre']]


_FAMILY: Dict[str, pd.DataFrame] = {}


def family_levels() -> Dict[str, pd.DataFrame]:
    """The family's APL-on rows (wiring, stimulus, win, pre) at the blank's level ('blank', carrying) and the settled
    level ('settled', partner (6)), from E98's B1 runs (`load_b1_all`), read once; attrs 'sources' on each."""
    if not _FAMILY:
        E98 = _e98()
        b1 = E98.load_b1_all(E98.ALL_WIRINGS)
        missing = b1.attrs.get('missing') or []
        if missing:
            raise RuntimeError(f'B1 runs missing: {missing}')
        for lv in ('blank', 'settled'):
            _FAMILY[lv] = _on(E98.with_a0(b1, lv))
            _FAMILY[lv].attrs['sources'] = dict(b1.attrs.get('sources') or {})
    return _FAMILY


def family_identity() -> Tuple[dict, pd.DataFrame]:
    """The family's ratio sets (`E98.stimulus_classes` on B1's rows at the blank's and the settled levels) against E98's
    c_wirings.csv: every row of V1's comparisons (each estimate of ms's input) at those two levels, per wiring (a_num,
    a_den, a0, ratio_linear, ratio_min, ratio_max within IDENTITY_TOL, and the class)."""
    E98 = _e98()
    cw = pd.read_csv(E98.out_dir() / 'c_wirings.csv', float_precision='round_trip')
    cw = cw[cw['comparison'].str.startswith("V1, ms's input from ") & cw['model_level'].isin(('blank', 'settled'))]
    lv = family_levels()
    rows = []
    cols = ('a_num', 'a_den', 'a0', 'ratio_linear', 'ratio_min', 'ratio_max')
    for (comp, rs, level), g in cw.groupby(['comparison', 'reading_set', 'model_level'], sort=False):
        est = comp[len("V1, ms's input from "):]
        lab = E98.V1_ESTIMATES[est]
        lo, hi = float(g['flies_lo90'].iloc[0]), float(g['flies_hi90'].iloc[0])
        wc = E98.stimulus_classes(lv[level], E98.V1_NUM, lab, lo, hi).set_index('wiring')
        for _, x in g.iterrows():
            w = x['wiring']
            d = max(abs(float(wc.loc[w, c]) - float(x[c])) for c in cols) if w in wc.index else np.nan
            rows.append(dict(comparison=comp, reading_set=rs, model_level=level, wiring=w, max_abs_diff=d,
                             cls_e98=x['cls'], cls_e103=wc.loc[w, 'cls'] if w in wc.index else 'missing',
                             same=bool(w in wc.index and d <= IDENTITY_TOL and wc.loc[w, 'cls'] == x['cls'])))
    t = pd.DataFrame(rows)
    res = dict(n_rows=int(len(t)), n_comparisons=int(t[['comparison', 'reading_set']].drop_duplicates().shape[0]),
               max_abs_diff=float(t['max_abs_diff'].max()), passed=bool(len(t) and t['same'].all()),
               b1_sources=lv['blank'].attrs.get('sources'))
    return res, t


def twelve_identity(ratios_12: Optional[pd.DataFrame] = None) -> dict:
    """The 12's P1 and Q1 (and, beside, Q2 and F1) under the carrying reading against E98's a3_ratios.csv (reading
    order_matched, subset all): ratio, lo90, hi90, sd_log2_flies, n_flies, n_recordings within IDENTITY_TOL."""
    E98 = _e98()
    if ratios_12 is None:
        ratios_12 = flies_ratios_12('carrying')[0]
    a3 = pd.read_csv(E98.out_dir() / 'a3_ratios.csv', float_precision='round_trip')
    a3 = a3[(a3['reading'] == 'order_matched') & (a3['subset'] == 'all')].set_index('ratio_id')
    mine = ratios_12.set_index('ratio_id')
    out = {}
    for rid, eid in E98_IDS.items():
        d = max(abs(float(mine.loc[rid, c]) - float(a3.loc[eid, c]))
                for c in ('ratio', 'lo90', 'hi90', 'sd_log2_flies', 'n_flies', 'n_recordings'))
        out[rid] = dict(e98_id=eid, ratio=float(mine.loc[rid, 'ratio']), lo90=float(mine.loc[rid, 'lo90']),
                        hi90=float(mine.loc[rid, 'hi90']), e98_ratio=float(a3.loc[eid, 'ratio']),
                        e98_lo90=float(a3.loc[eid, 'lo90']), e98_hi90=float(a3.loc[eid, 'hi90']), max_abs_diff=d,
                        equal=bool(d <= IDENTITY_TOL), gates=rid in IDENTITY_IDS)
    out['passed'] = bool(all(out[r]['equal'] for r in IDENTITY_IDS))
    return out


def part_0() -> dict:
    """Part 0: writes p0_eligibility.csv, p0_vials.csv, p0_reader_identity.csv, p0_family_identity.csv and
    p0_identity.json (no map of the series opened: yaml, recordings.csv and the date notes only); raises if an identity
    fails."""
    out = out_dir()
    el = eligibility_table()
    write_table(el, out / 'p0_eligibility.csv')
    s = el[el['stem'].isin(SERIES)]
    by_pres = {int(k): int(v) for k, v in s['odor_presentations_per_block'].value_counts().items()}
    carrying = s[s['carrying'].astype(bool)]
    elig = dict(passed=bool(s['matches'].astype(bool).all() and len(s) == len(SERIES)
                            and by_pres == N_RUNS_BY_PRESENTATIONS
                            and int((s['n_pfo_only_blocks'] == 2).sum()) == 10 and int((s['n_pfo_only_blocks'] == 0).sum()) == 4
                            and len(carrying) == N_CARRYING_RECORDINGS
                            and sorted(carrying['fly'].unique()) == sorted(SERIES_FLIES)
                            and s['frames_ok'].astype(bool).all()),
                n_series=int(len(s)), n_matches=int(s['matches'].astype(bool).sum()),
                mismatches=s.loc[~s['matches'].astype(bool), 'stem'].tolist(), runs_by_odor_presentations=by_pres,
                n_two_pfo_blocks=int((s['n_pfo_only_blocks'] == 2).sum()), n_no_pfo_block=int((s['n_pfo_only_blocks'] == 0).sum()),
                n_frames_ok=int(s['frames_ok'].astype(bool).sum()), n_carrying_recordings=int(len(carrying)),
                carrying_flies=sorted(carrying['fly'].unique()), n_raw_only_listed=int(len(el) - len(s)),
                frames_source="recordings.csv's peak_map_shape and frames_equal_presentations (no map opened)")
    vi = vials_identity(el)
    write_table(vi, out / 'p0_vials.csv')
    vials = dict(passed=bool(vi['agrees'].all()), n=int(len(vi)), disagree=vi.loc[~vi['agrees'], 'stem'].tolist())
    rd, rd_per = reader_identity()
    write_table(rd_per, out / 'p0_reader_identity.csv')
    fam, fam_rows = family_identity()
    write_table(fam_rows, out / 'p0_family_identity.csv')
    tw = twelve_identity()
    res = dict(caveat=QUALITY_CAVEAT, eligibility=elig, vials=vials, reader=rd, family=fam, twelve=tw,
               passed=bool(elig['passed'] and vials['passed'] and rd['passed'] and fam['passed'] and tw['passed']),
               host=_host(), time=time.strftime('%Y-%m-%d %H:%M'))
    _write_json(out / 'p0_identity.json', res)
    provenance('0', passed=res['passed'], files={str(p): _sha256(p) for p in (
        _e98().out_dir() / 'a1_responses.csv', _e98().out_dir() / 'd_responses.csv', _e98().out_dir() / 'c_wirings.csv',
        _e98().out_dir() / 'a3_ratios.csv', H.shared_dir().joinpath('data', 'yang_apl_recordings', 'recordings.csv'))})
    if not res['passed']:
        raise RuntimeError(f"part 0: an identity fails (p0_identity.json): eligibility {elig['passed']}, vials "
                           f"{vials['passed']}, reader {rd['passed']}, family {fam['passed']}, the 12 {tw['passed']}")
    return res


def part_0_passed(path: Optional[Path] = None) -> bool:
    d = _read_json(out_dir() / 'p0_identity.json' if path is None else path)
    return bool(d and d.get('passed'))


# ---------------------------------------------------------------------------
# ratios, intervals and classes

_A1: Dict[str, pd.DataFrame] = {}


def a1_responses() -> pd.DataFrame:
    """E98's part A responses of the 12 (a1_responses.csv), read once (round-trip floats)."""
    if 'a1' not in _A1:
        _A1['a1'] = pd.read_csv(_e98().out_dir() / 'a1_responses.csv', float_precision='round_trip')
    return _A1['a1']


def flies_ratios_12(reading: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The 12's eight ratios under one of READINGS_12 (`E98.flies_ratios` on a1_responses.csv, the 12's design ranges):
    (one row per ratio, per recording rows)."""
    E98 = _e98()
    rd, subset = READINGS_12[reading]
    f, pr = E98.flies_ratios(a1_responses(), rd, subset, RATIOS)
    f.insert(0, 'reading_12', reading)
    pr.insert(0, 'reading_12', reading)
    return f, pr


def t_interval(per_fly_log: Sequence[float], level: float = CI_LEVEL) -> Tuple[float, float, float]:
    """The geometric mean over flies and its `level` interval from Student's t on the flies' (natural) log ratios
    (n - 1 degrees of freedom); NaN bounds with fewer than two flies."""
    from scipy.stats import t as student_t
    x = np.asarray(per_fly_log, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    m = float(x.mean())
    if len(x) < 2:
        return float(np.exp(m)), np.nan, np.nan
    h = float(student_t.ppf(0.5 + level / 2.0, len(x) - 1)) * float(x.std(ddof=1)) / np.sqrt(len(x))
    return float(np.exp(m)), float(np.exp(m - h)), float(np.exp(m + h))


def split_rule(set_a: pd.Series, set_b: pd.Series) -> dict:
    """E98's `pulse_permutation` on two sets of flies' mean (natural) log ratios (Series indexed by fly): over every
    split of the pooled flies into sets of the two sizes, the share whose difference of mean log ratios is at least as
    large (absolute) as the observed; 'differ' when that share is below SPLIT_SHARE. Adds the two sets' geometric means
    ('ratio_a', 'ratio_b')."""
    E98 = _e98()
    if abs(E98.PULSE_SHARE - SPLIT_SHARE) > 0:
        raise RuntimeError("E98's PULSE_SHARE is not the entry's split share")
    a, b = set_a.dropna(), set_b.dropna()
    s = pd.concat([pd.Series(a.to_numpy(float), index=[f'a:{f}' for f in a.index]),
                   pd.Series(b.to_numpy(float), index=[f'b:{f}' for f in b.index])])
    r = E98.pulse_permutation(s, [f'b:{f}' for f in b.index])
    r = {k: v for k, v in r.items() if k not in ('ratio_short', 'ratio_long')}
    r.update(n_a=int(len(a)), n_b=int(len(b)), ratio_a=float(np.exp(a.mean())) if len(a) else np.nan,
             ratio_b=float(np.exp(b.mean())) if len(b) else np.nan, differ=bool(r['disagree']))
    return r


_SPLITS: Dict[Tuple[int, int], np.ndarray] = {}


def split_share(a: Sequence[float], b: Sequence[float]) -> Tuple[float, int]:
    """The split rule's share, vectorized (the same count as `E98.pulse_permutation`): over every choice of len(b) of
    the pooled values as the second set, the share whose |difference of means| is at least the observed (minus 1e-12);
    and the number of splits."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n, k = len(a) + len(b), len(b)
    if (n, k) not in _SPLITS:
        _SPLITS[(n, k)] = np.array(list(combinations(range(n), k)), dtype=int)
    idx = _SPLITS[(n, k)]
    v = np.concatenate([a, b])
    tot = v.sum()
    sb = v[idx].sum(axis=1)
    d = np.abs(sb / k - (tot - sb) / (n - k))
    obs = abs(b.mean() - a.mean())
    return float((d >= obs - 1e-12).sum() / len(idx)), int(len(idx))


def normal_scores(n: int) -> np.ndarray:
    """Blom's normal scores of `n` values, standardized to mean 0 and sample SD 1."""
    from scipy.stats import norm
    z = norm.ppf((np.arange(1, n + 1) - 0.375) / (n + 0.25))
    return (z - z.mean()) / z.std(ddof=1)


def _least_shift(a: np.ndarray, b0: np.ndarray, sign: float, hi: float) -> float:
    """The least shift d >= 0 of `b0` (centred on a's mean) by sign * d at which the split rule says differ."""
    def differs(d):
        return split_share(a, b0 + sign * d)[0] < SPLIT_SHARE
    if not differs(hi):
        return np.nan
    lo = 0.0
    if differs(lo):
        return 0.0
    while hi - lo > LDD_TOL:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if differs(mid) else (mid, hi)
    return hi


def least_detectable(per_fly_log2: Sequence[float], n_new: int = N_SERIES_FLIES, sd_log2: Optional[float] = None) -> dict:
    """The least difference of mean log2 ratios that `n_new` new flies can show against the given flies by the split rule:
    the new flies spread like the given ones (normal scores scaled to `sd_log2`, default the given flies' between-fly
    SD), their mean the given flies' mean shifted up ('up') or down ('down'), bisected to LDD_TOL. Carrying: the given
    flies' own values; beside, both sets as normal scores ('_normal_scores', symmetric) and the normal approximation
    NORMAL_Z90 sd sqrt(1/n + 1/n_new) ('_normal_approx'). Factors are 2^difference."""
    a = np.asarray(per_fly_log2, dtype=float)
    a = a[np.isfinite(a)]
    n = len(a)
    sd = float(a.std(ddof=1)) if sd_log2 is None else float(sd_log2)
    out = dict(n_flies=n, n_new=int(n_new), sd_log2=sd, n_splits=int(math.comb(n + n_new, n_new)))
    if n < 2 or not np.isfinite(sd) or sd <= 0:
        return out
    hi = 20.0 * sd + 1.0
    b0 = a.mean() + normal_scores(n_new) * sd
    up, down = _least_shift(a, b0, 1.0, hi), _least_shift(a, b0, -1.0, hi)
    an = normal_scores(n) * sd
    ns = _least_shift(an, normal_scores(n_new) * sd, 1.0, hi)
    approx = NORMAL_Z90 * sd * np.sqrt(1.0 / n + 1.0 / n_new)
    out.update(ldd_up_log2=up, ldd_down_log2=down, factor_up=float(2 ** up), factor_down=float(2 ** -down),
               ldd_log2_normal_scores=ns, factor_normal_scores=float(2 ** ns), ldd_log2_normal_approx=float(approx),
               factor_normal_approx=float(2 ** approx))
    return out


def class_ab(ratios: Sequence[float], lo: float, hi: float) -> str:
    """One wiring's set against [lo, hi] (compared on the log scale, as E98's `classify`): 'inside' if any lies in it,
    'above' if all lie above `hi`, 'below' if all lie below `lo`, else 'mixed'."""
    r = np.asarray(ratios, dtype=float)
    if not np.all(np.isfinite(r) & (r > 0)):
        raise ValueError(f'ratios must be positive and finite: {ratios}')
    x = np.log(r)
    L, U = np.log(float(lo)), np.log(float(hi))
    if not L <= U:
        raise ValueError(f'empty interval [{lo}, {hi}]')
    if np.any((x >= L) & (x <= U)):
        return 'inside'
    if np.all(x > U):
        return 'above'
    if np.all(x < L):
        return 'below'
    return 'mixed'


def comparison_class(classes: Sequence[str], min_wirings: int = MIN_WIRINGS) -> str:
    """'inside', 'above' or 'below' if at least `min_wirings` wirings are; otherwise 'not decided'."""
    c = list(classes)
    for v in X_CLASSES:
        if sum(x == v for x in c) >= min_wirings:
            return v
    return 'not decided'


def held_class(categories: Dict[str, str], partners: Sequence[str]) -> str:
    """The carrying reading's class ('carrying') when it is one of X_CLASSES and every partner in `partners` gives it;
    otherwise 'reading-dependent'."""
    c = categories['carrying']
    if c in X_CLASSES and all(categories[k] == c for k in partners):
        return c
    return 'reading-dependent'


def x_over_estimates(held: Dict[str, str]) -> str:
    """X2 / X3 over the estimates of ms's input: 'inside' if held inside under any; 'above' / 'below' if held so under
    every one; otherwise 'not decided'."""
    v = list(held.values())
    if not v:
        raise ValueError("no estimate of ms's input")
    if 'inside' in v:
        return 'inside'
    for d in ('above', 'below'):
        if all(x == d for x in v):
            return d
    return 'not decided'


def wiring_classes(on: pd.DataFrame, num: str, den: str, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring the model's ratio set (`E98.stimulus_classes`' rows: a_num, a_den, a0, the linear ratio, the set's
    range, E98's class) with this entry's class of the set ('cls_ab') and of the linear reading ('cls_ab_linear')
    against [lo, hi]."""
    E98 = _e98()
    wc = E98.stimulus_classes(on, num, den, lo, hi)
    ab, ab_lin = [], []
    for _, x in wc.iterrows():
        if x['cls'] == 'no ratio':
            ab.append('no ratio')
            ab_lin.append('no ratio')
            continue
        rs = E98.model_ratio_set(x['a_num'], x['a_den'], x['a0'])['ratio']
        ab.append(class_ab(rs, lo, hi))
        ab_lin.append(class_ab([x['ratio_linear']], lo, hi))
    wc['cls_ab'] = ab
    wc['cls_ab_linear'] = ab_lin
    return wc


def _hull(ivs: Sequence[Tuple[float, float]]) -> Tuple[float, float]:
    return _e98().hull(ivs)


# ---------------------------------------------------------------------------
# part A

def ratios_12_tables() -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, pd.DataFrame]]:
    """The 12's eight ratios under READINGS_12 (rows with the between-fly SD, the Student's t interval of the carrying
    reading, the entry's estimate and set S's value beside), the per-fly mean log ratios, and the per-recording rows by
    reading."""
    rows, per_fly, per = [], [], {}
    for rk in READINGS_12:
        f, pr = flies_ratios_12(rk)
        per[rk] = pr
        for _, x in f.iterrows():
            rid = x['ratio_id']
            e = pr[(pr['ratio_id'] == rid) & pr['entered']]
            pf = e.groupby('fly')['log_ratio'].agg(['mean', 'size'])
            _, tlo, thi = t_interval(pf['mean'].to_numpy())
            d = {k: x[k] for k in ('reading_12', 'ratio_id', 'group', 'numerator', 'denominator', 'reading', 'subset',
                                   'ratio', 'lo90', 'hi90', 'sd_log2_flies', 'n_flies', 'n_recordings', 'n_flies_above_1',
                                   'n_recordings_in_subset', 'n_recordings_left_out', 'num_how', 'den_how')}
            d.update(t_lo90=tlo, t_hi90=thi, e98_id=E98_IDS.get(rid, ''), entry_ratio=ENTRY_12[rid][0],
                     entry_lo90=ENTRY_12[rid][1], entry_hi90=ENTRY_12[rid][2], set_s_lo=SET_S[rid][0],
                     set_s_hi=SET_S[rid][1])
            rows.append(d)
            for fly, y in pf.iterrows():
                per_fly.append(dict(reading_12=rk, ratio_id=rid, fly=fly, n_recordings=int(y['size']),
                                    mean_log_ratio=float(y['mean']), mean_log2_ratio=float(y['mean'] / np.log(2)),
                                    ratio=float(np.exp(y['mean']))))
    return pd.DataFrame(rows), pd.DataFrame(per_fly), per


def ldd_table(per_fly: pd.DataFrame, ratios: pd.DataFrame) -> pd.DataFrame:
    """Per ratio, the least difference 5 flies can show against the 12's flies by the split rule (`least_detectable` on
    the carrying reading's per-fly log2 ratios), with the 12's ratio and interval, the bounds a held-out ratio must pass
    (up and down) and where set S's value lies against them."""
    rows = []
    rt = ratios[ratios['reading_12'] == 'carrying'].set_index('ratio_id')
    for rid in RATIO_IDS:
        v = per_fly[(per_fly['reading_12'] == 'carrying') & (per_fly['ratio_id'] == rid)]['mean_log2_ratio'].to_numpy()
        ld = least_detectable(v)
        r12 = float(rt.loc[rid, 'ratio'])
        s_lo, s_hi = SET_S[rid]
        d = dict(ratio_id=rid, numerator=rt.loc[rid, 'numerator'], denominator=rt.loc[rid, 'denominator'],
                 ratio_12=r12, lo90_12=float(rt.loc[rid, 'lo90']), hi90_12=float(rt.loc[rid, 'hi90']), **ld)
        if 'ldd_up_log2' in ld:
            d.update(detectable_above=r12 * ld['factor_up'], detectable_below=r12 * ld['factor_down'],
                     set_s_lo=s_lo, set_s_hi=s_hi,
                     set_s_beyond=bool(s_lo >= r12 * ld['factor_up'] or s_hi <= r12 * ld['factor_down']),
                     entry_estimate=(f'a factor of about {LDD_ENTRY_P1[0]} from a between-fly SD of {LDD_ENTRY_P1[1]} '
                                     '(log2)') if rid == 'P1' else '')
        rows.append(d)
    return pd.DataFrame(rows)


def _flies_interval(pr: pd.DataFrame, rid: str, drop_fly: Optional[str] = None) -> Tuple[float, float]:
    e = pr[(pr['ratio_id'] == rid) & pr['entered']]
    if drop_fly is not None:
        e = e[e['fly'] != drop_fly]
    fr = _e98().fly_ratio(e)
    return fr['lo90'], fr['hi90']


def classes_12(ratios: pd.DataFrame, per: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The 12's classes for X2 (the hull of P1's and P2's intervals) and X3 (P1's), as descriptions: per estimate of
    ms's input (E98's V1_ESTIMATES) and reading (the carrying reading; partners (1) first-in-block and (3) E98's
    without_three on the flies' side, (4) each fly of the 12 left out, (5) Student's t, all with the model at the
    blank's level; (6) the model at the settled level with the flies' carrying reading), the class of the eight wirings'
    sets (`wiring_classes`, `comparison_class`) and of their linear readings, the held class per estimate and X2 / X3
    over the estimates (`x_over_estimates`). Returns (one row per X, estimate and reading; per wiring rows; summary)."""
    E98 = _e98()
    lv = family_levels()
    rt = {rk: ratios[ratios['reading_12'] == rk].set_index('ratio_id') for rk in READINGS_12}
    flies = sorted(per['carrying'].loc[per['carrying']['entered'], 'fly'].unique())
    pf = {rid: per['carrying'][(per['carrying']['ratio_id'] == rid) & per['carrying']['entered']]
          .groupby('fly')['log_ratio'].mean() for rid in RATIO_IDS}

    def intervals(rk: str) -> Dict[str, Tuple[float, float]]:
        if rk in READINGS_12:
            return {rid: (float(rt[rk].loc[rid, 'lo90']), float(rt[rk].loc[rid, 'hi90'])) for rid in RATIO_IDS}
        if rk == 't_interval':
            return {rid: t_interval(pf[rid].to_numpy())[1:] for rid in RATIO_IDS}
        if rk.startswith('without '):
            return {rid: _flies_interval(per['carrying'], rid, rk[len('without '):]) for rid in ('P1', 'P2')}
        raise ValueError(rk)

    readings = [('carrying', 'carrying', 'blank'), ('first_in_block', 'first_in_block', 'blank'),
                ('without_marked', 'without_marked', 'blank')]
    readings += [(f'without {f}', f'without {f}', 'blank') for f in flies]
    readings += [('t_interval', 't_interval', 'blank'), ('settled', 'carrying', 'settled')]
    ivs = {flies_key: intervals(flies_key) for _, flies_key, _ in readings}
    rows, wrows, summ = [], [], {}
    partners = [k for k, _, _ in readings if k != 'carrying']
    for x, ids in X_IDS.items():
        held = {}
        for est, lab in E98.V1_ESTIMATES.items():
            cats = {}
            for key, fk, level in readings:
                lo, hi = _hull([ivs[fk][i] for i in ids])
                wc = wiring_classes(lv[level], E98.V1_NUM, lab, lo, hi)
                cats[key] = comparison_class(wc['cls_ab'])
                rows.append(dict(x=x, estimate=est, reading=key, flies_reading=fk, model_level=level, flies_lo90=lo,
                                 flies_hi90=hi, n_wirings=int(len(wc)), category=cats[key],
                                 category_linear=comparison_class(wc['cls_ab_linear']),
                                 **{f'n_{k}': int((wc['cls_ab'] == k).sum()) for k in X_CLASSES + ('mixed',)},
                                 model_linear_min=float(wc['ratio_linear'].min()),
                                 model_linear_max=float(wc['ratio_linear'].max()),
                                 model_min=float(wc['ratio_min'].min()), model_max=float(wc['ratio_max'].max())))
                if key in ('carrying', 'settled'):
                    wrows.append(wc.assign(x=x, estimate=est, reading=key, model_level=level))
            held[est] = held_class(cats, partners)
            for r in rows:
                if r['x'] == x and r['estimate'] == est:
                    r['held_per_estimate'] = held[est]
        summ[x] = dict(over_estimates=x_over_estimates(held), held_per_estimate=held,
                       carrying_per_estimate={r['estimate']: r['category'] for r in rows
                                              if r['x'] == x and r['reading'] == 'carrying'})
    cl = pd.DataFrame(rows)
    for x in X_IDS:
        cl.loc[cl['x'] == x, 'x_over_estimates'] = summ[x]['over_estimates']
    wcols = ['x', 'estimate', 'reading', 'model_level', 'wiring', 'numerator', 'denominator', 'a_num', 'a_den', 'a0',
             'ratio_linear', 'ratio_min', 'ratio_max', 'flies_lo90', 'flies_hi90', 'cls_ab', 'cls_ab_linear', 'cls']
    w = pd.concat(wrows, ignore_index=True)[wcols].rename(columns={'cls': 'cls_e98'})
    summ['readings'] = [dict(key=k, flies=fk, model_level=lv_) for k, fk, lv_ in readings]
    summ['flies_12'] = flies
    return cl, w, summ


def part_a(results: Optional[Path] = None) -> dict:
    """Part A (no map of the series opened; refuses to run unless part 0 passed): the 12 at the series' stimuli
    (a_ratios_12.csv, a_per_fly_12.csv), the least detectable differences (a_least_detectable.csv), the 12's classes for
    X2 and X3 as descriptions (a_classes_12.csv, a_class_wirings_12.csv), a_summary.json; then the lock (a_lock.json in
    the results and PART_A_LOCK in this repository, the sha256 of every table). With the lock already in place it writes
    to a_recheck/ and reports whether the tables are the same, never over the locked ones."""
    E98 = _e98()
    if not part_0_passed():
        raise RuntimeError("part A runs after part 0's identities pass (p0_identity.json)")
    out = out_dir() if results is None else Path(results)
    recheck = PART_A_LOCK.exists()
    dest = out / 'a_recheck' if recheck else out
    dest.mkdir(parents=True, exist_ok=True)
    ratios, per_fly, per = ratios_12_tables()
    tw = twelve_identity(ratios[ratios['reading_12'] == 'carrying'])
    if not tw['passed']:
        raise RuntimeError(f"the 12's P1 and Q1 differ from E98's a3_ratios.csv: {tw}")
    ldd = ldd_table(per_fly, ratios)
    cl, wcl, csum = classes_12(ratios, per)
    write_table(ratios, dest / 'a_ratios_12.csv')
    write_table(per_fly, dest / 'a_per_fly_12.csv')
    write_table(ldd, dest / 'a_least_detectable.csv')
    write_table(cl, dest / 'a_classes_12.csv')
    write_table(wcl, dest / 'a_class_wirings_12.csv')
    v98 = _read_json(E98.out_dir() / 'verdict.json') or {}
    car = ratios[ratios['reading_12'] == 'carrying'].set_index('ratio_id')
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED,
                note=("Part A reads no map or response file of the 2021-10/11 series; the 12 are E98's part A responses. "
                      "The 12's classes are descriptions; E98's verdict stands as written."),
                e98_verdict=dict(V1=v98.get('V1'), V1_held_per_estimate=v98.get('V1_held_per_estimate')),
                ratios_12_carrying={rid: dict(ratio=float(car.loc[rid, 'ratio']), lo90=float(car.loc[rid, 'lo90']),
                                              hi90=float(car.loc[rid, 'hi90']), sd_log2_flies=float(car.loc[rid, 'sd_log2_flies']),
                                              n_flies=int(car.loc[rid, 'n_flies']), n_recordings=int(car.loc[rid, 'n_recordings']))
                                    for rid in RATIO_IDS},
                least_detectable={r['ratio_id']: dict(factor_up=r.get('factor_up'), factor_down=r.get('factor_down'),
                                                      factor_normal_approx=r.get('factor_normal_approx'),
                                                      sd_log2=r.get('sd_log2'), n_flies=r.get('n_flies'),
                                                      n_splits=r.get('n_splits'))
                                  for r in ldd.to_dict('records')},
                classes_12={x: dict(over_estimates=csum[x]['over_estimates'], held_per_estimate=csum[x]['held_per_estimate'],
                                    carrying_per_estimate=csum[x]['carrying_per_estimate']) for x in X_IDS},
                readings_12=dict(READINGS_12={k: dict(reading=list(v[0]), subset=v[1]) for k, v in READINGS_12.items()},
                                 class_readings=csum['readings'],
                                 partner_2='the carrying reading itself for the 12 (two presentations per block, one pass)'),
                twelve_identity={k: v for k, v in tw.items()},
                family_sources=family_levels()['blank'].attrs.get('sources'),
                settings=dict(N_BOOT=E98.N_BOOT, SEED=E98.SEED, MIN_RESPONSE=E98.MIN_RESPONSE, SPLIT_SHARE=SPLIT_SHARE,
                              N_SERIES_FLIES=N_SERIES_FLIES, MIN_WIRINGS=MIN_WIRINGS, CI_LEVEL=CI_LEVEL,
                              LDD_TOL=LDD_TOL, RATIOS=[list(r) for r in RATIOS], SET_S=SET_S))
    _write_json(dest / 'a_summary.json', summ)
    tables = {name: _sha256(dest / name) for name in PART_A_TABLES}
    res = dict(summary=summ, tables=tables, dest=str(dest))
    if recheck:
        locked = json.loads(PART_A_LOCK.read_text()).get('tables', {})
        res['recheck_same'] = {k: bool(locked.get(k) == v) for k, v in tables.items()}
        print('part A re-run into a_recheck/; same as the locked tables:', res['recheck_same'])
        return res
    lock = dict(experiment=EXPERIMENT, part='A', results=str(out), tables=tables, host=_host(),
                time=time.strftime('%Y-%m-%d %H:%M'), commits={k: v for k, v in H.code_commits().items()
                                                                if not isinstance(v, (dict, list))},
                note=("Part A's tables, written before any map or response file of the 2021-10/11 series was opened. "
                      "Commit this file before running any later part: `require_part_a` refuses to open the series' "
                      "maps until it is committed and the tables match it."))
    _write_json(out / 'a_lock.json', lock)
    _write_json(PART_A_LOCK, lock)
    provenance('A', tables=tables, lock=str(PART_A_LOCK))
    return res


# ---------------------------------------------------------------------------
# the rules of parts B to E (pure functions)

def x1_direction(ratio_id: str, series_ratio: float, ratio_12: float) -> dict:
    """X1 for one ratio: whether the held-out flies' ratio lies from the 12's toward set S's value ('toward_s'; set S's
    value is a point or, for P2, a range: the direction to its nearer end) and whether it is nearer to set S's value on
    the log scale than to the 12's ('nearer_s'; the distance to a range is 0 inside it)."""
    s_lo, s_hi = SET_S[ratio_id]
    x, r = np.log(float(series_ratio)), np.log(float(ratio_12))
    L, U = np.log(s_lo), np.log(s_hi)
    target = L if abs(L - r) <= abs(U - r) else U
    toward = bool(np.sign(x - r) == np.sign(target - r) and x != r)
    d_s = 0.0 if L <= x <= U else min(abs(x - L), abs(x - U))
    return dict(toward_s=toward, nearer_s=bool(d_s < abs(x - r)), log_dist_s=float(d_s), log_dist_12=float(abs(x - r)))


def x1_side(rows: pd.DataFrame) -> str:
    """X1 from the six ratios with ms (columns ratio_id, differ, toward_s, nearer_s): "the 12's side" (at most one
    differs), "the 2022 flies' side" (at least four differ, each toward set S and nearer to set S's value on the log
    scale than to the 12's), else 'between'."""
    r = rows.set_index('ratio_id').loc[list(X1_RATIOS)]
    d = r['differ'].astype(bool)
    if int(d.sum()) <= X1_MAX_DIFFER_12:
        return SIDES[0]
    dd = r[d]
    if int(d.sum()) >= X1_MIN_DIFFER_S and bool((dd['toward_s'].astype(bool) & dd['nearer_s'].astype(bool)).all()):
        return SIDES[1]
    return SIDES[2]


def y1_variance_split(rows: pd.DataFrame) -> dict:
    """Y1: from per-recording log2 ratios of two vials of one odor (columns day, recording, odor, conc, pair,
    log2_ratio): 'sd_single' the SD over single recordings; per day, odor, concentration and pair of vials the mean over
    that day's recordings; 'persistent_sd' the SD of those means with the scatter of single recordings about them removed
    (the pooled within-group variance over the mean group size's inverse), floored at 0."""
    d = rows.dropna(subset=['log2_ratio'])
    g = d.groupby(['day', 'odor', 'conc', 'pair'])['log2_ratio']
    means, n = g.mean(), g.size()
    dev = d['log2_ratio'] - g.transform('mean')
    df_w = int((n - 1).sum())
    s_w2 = float((dev ** 2).sum() / df_w) if df_w > 0 else np.nan
    var_means = float(means.var(ddof=1)) if len(means) > 1 else np.nan
    noise = s_w2 * float((1.0 / n).mean()) if np.isfinite(s_w2) else np.nan
    pv = var_means - noise if np.isfinite(var_means) and np.isfinite(noise) else np.nan
    return dict(n_values=int(len(d)), n_groups=int(len(means)), sd_single=float(d['log2_ratio'].std(ddof=1)),
                sd_within_groups=float(np.sqrt(s_w2)) if np.isfinite(s_w2) else np.nan,
                sd_group_means=float(np.sqrt(var_means)) if np.isfinite(var_means) else np.nan,
                persistent_var=pv, persistent_sd=float(np.sqrt(max(0.0, pv))) if np.isfinite(pv) else np.nan,
                mean_of_group_means=float(means.mean()))


def vial_yardstick(persistent_sd: float) -> float:
    """The log2 difference between two series beyond the vials: Y1_YARDSTICK x the persistent SD."""
    return float(Y1_YARDSTICK * persistent_sd)


def _blurred_in(m: np.ndarray, mask: np.ndarray, blur: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(np.asarray(m, dtype=float), blur)[mask]


def m_calyx(pres: Dict[object, List[np.ndarray]], common: Sequence[object], vial_pairs: Sequence[Tuple[object, object]],
            odor_pairs: Sequence[Tuple[object, object]], mask: np.ndarray, blur: float = M_BLUR_PX) -> dict:
    """M for one calyx: `pres` the presentation maps of each stimulus in its block; the common map the first singular
    vector of the blurred mean maps (means removed) of the `common` stimuli (ms and vial 1 of paa and pa; E79's
    `svd_shared`); a stimulus's remainder its blurred map, mean removed, without its part along the common map. r_rep:
    per stimulus of the pairs, r between the remainders of the first and second halves of its presentations; r_vial:
    per pair in `vial_pairs`, r between the two stimuli's remainders; r_odor: per pair in `odor_pairs` likewise. The
    calyx's values are the means over pairs (r_rep over the stimuli of `vial_pairs`), and its share `m_share`."""
    from experiments import E100_apl_weak_odors as E100
    DS = E100._ds()
    X = np.stack([_blurred_in(np.mean(pres[k], axis=0), mask, blur) for k in common])
    _, v = DS.svd_shared(X, True)

    def rem(m):
        return DS._residual(_blurred_in(m, mask, blur), v, True)

    def r(a, b):
        return float(np.corrcoef(a, b)[0, 1])
    remainder = {k: rem(np.mean(p, axis=0)) for k, p in pres.items()}
    reps = []
    for k in sorted({x for p in vial_pairs for x in p}, key=str):
        p = pres[k]
        h = len(p) // 2
        if h >= 1:
            reps.append(r(rem(np.mean(p[:h], axis=0)), rem(np.mean(p[h:], axis=0))))
    rv = [r(remainder[a], remainder[b]) for a, b in vial_pairs]
    ro = [r(remainder[a], remainder[b]) for a, b in odor_pairs]
    out = dict(r_rep=float(np.mean(reps)) if reps else np.nan, r_vial=float(np.mean(rv)) if rv else np.nan,
               r_odor=float(np.mean(ro)) if ro else np.nan, n_vial_pairs=len(rv), n_odor_pairs=len(ro))
    out['share'] = m_share(out['r_vial'], out['r_odor'], out['r_rep'])
    return out


def m_share(r_vial: float, r_odor: float, r_rep: float) -> float:
    """M's share per calyx: (r_vial - r_odor) / (r_rep - r_odor); 1 if all of what repeats within a block returns with
    another vial of the same odor in another block, 0 if none does."""
    den = float(r_rep) - float(r_odor)
    return float((float(r_vial) - float(r_odor)) / den) if den != 0 else np.nan


def m_class(share_lo: float, share_hi: float, n_vial_over_odor: int) -> str:
    """M: "the odor's" (the share's interval above M_ODOR_LO and r_vial > r_odor in at least M_MIN_CALYCES calyces),
    "the block's" (the interval below M_BLOCK_HI), else 'not decided'."""
    if share_lo > M_ODOR_LO and int(n_vial_over_odor) >= M_MIN_CALYCES:
        return "the odor's"
    if share_hi < M_BLOCK_HI:
        return "the block's"
    return 'not decided'


def z_rule(log_diff_per_calyx: Sequence[float], lo: float, hi: float) -> str:
    """Z: per calyx log P1 at 15 um minus at 18 um; 'replicates' if at least Z_MIN_HIGHER calyces are higher at 15 um
    and the mean factor's interval over flies [lo, hi] excludes 1; 'does not' if Z_MAX_HIGHER_NOT or fewer are higher;
    else 'not decided'."""
    x = np.asarray(log_diff_per_calyx, dtype=float)
    n_hi = int((x[np.isfinite(x)] > 0).sum())
    if n_hi >= Z_MIN_HIGHER and not (lo <= 1.0 <= hi):
        return 'replicates'
    if n_hi <= Z_MAX_HIGHER_NOT:
        return 'does not'
    return 'not decided'


# ---------------------------------------------------------------------------
# parts B to E: the series' responses, readings and ratios

READINGS_SERIES = ('carrying', 'first_in_block', 'all_presentations', 'without_marked')     # then each fly left out
TWELVE_FOR = {'carrying': 'carrying', 'first_in_block': 'first_in_block', 'all_presentations': 'carrying',
              'without_marked': 'without_marked'}   # the 12's reading paired with each of the series' (part A's (1), (3))
TEN22_PAA = tuple(x for x in THREE_VIALS_PAA if x.startswith('20211022'))     # vial 1 in the other half: vials 2, 3 carry
TEN22_PA = tuple(x for x in SERIES if x.startswith('20211022') and x not in THREE_VIALS_PAA)
VIAL_PAIRS = ((1, 2), (1, 3), (2, 3))
ODORS3 = ('ms', 'paa', 'pa')                  # vial 1 of paa and pa; ms has one vial
Y1_READINGS = ('carrying', 'first_in_block', 'without_marked')
M_VARIANTS = ('carrying', 'before_common_removed', 'ten22_pa_beside')
Y3_MAX_REP = 4

_RESP: Dict[tuple, pd.DataFrame] = {}
_TABLES: Dict[str, pd.DataFrame] = {}


def _require_lock(part: str) -> None:
    if not part_a_locked():
        raise PermissionError(f"part {part} opens the series' maps: it runs only while part A's lock is committed and its "
                              "tables match it")


def series_responses(source: str = 'maps', which: str = 'first', n_first: Optional[int] = N_FIRST) -> pd.DataFrame:
    """The series' 30 runs read once (`read_responses`, both baselines; the maps or Yang's tables; the stimulus's
    first block or every block; the first `n_first` presentations or all)."""
    key = (source, which, n_first)
    if key not in _RESP:
        _RESP[key] = read_responses(list(SERIES), source, ('order_matched', 'first_in_block'), which, n_first)
    return _RESP[key]


def series_table(stem: str) -> pd.DataFrame:
    """One run's stimulus table (`stim_table`), read once."""
    if stem not in _TABLES:
        _TABLES[stem] = stim_table(series_rows([stem]).iloc[0])
    return _TABLES[stem]


_KEYS = ['recording', 'fly', 'side', 'pulse_s', 'plane_um', 'role', 'odor', 'conc', 'vial', 'stimulus', 'baseline',
         'blank_subtracted', 'source']


def every_presentation(resp_all: pd.DataFrame) -> pd.DataFrame:
    """Partner (2): per recording, stimulus and baseline the mean over every presentation of every kept block of the
    stimulus (both passes; pa vial 1's two series in the runs of 10-22), from `read_responses(..., which='all',
    n_first=None)`'s rows (each a block's mean over its `n_reps` presentations)."""
    d = resp_all.assign(_w=resp_all['response'] * resp_all['n_reps'])
    out = d.groupby(_KEYS, dropna=False).agg(_w=('_w', 'sum'), n_reps=('n_reps', 'sum'), n_blocks=('block1', 'nunique'),
                                              block1=('block1', 'min')).reset_index()
    out['response'] = out['_w'] / out['n_reps']
    out['block'] = out['block1'] - 1
    out['blank'] = np.nan
    return out.drop(columns='_w')


def three_vial_odor(labels: Sequence[str]) -> Optional[str]:
    """The run's three-vial odor from its stimulus labels ('pa' if 'pa_vial2' is among them, ...)."""
    v = sorted({str(o).split('_vial')[0] for o in labels if '_vial' in str(o)})
    return v[0] if len(v) == 1 else None


def average_three_vials(resp: pd.DataFrame) -> pd.DataFrame:
    """Part E's reading: per run, baseline and concentration, vial 1's response of the run's three-vial odor replaced by
    the mean of the three vials' responses (those imaged); the other stimuli unchanged."""
    out = resp.copy()
    for rec, g in resp.groupby('recording'):
        o = three_vial_odor(g['odor'])
        if o is None:
            continue
        v = g[g['odor'].isin([o, f'{o}_vial2', f'{o}_vial3'])]
        m = v.groupby(['baseline', 'conc'])['response'].mean()
        idx = out.index[(out['recording'] == rec) & (out['odor'] == o)]
        out.loc[idx, 'response'] = [float(m[(b, c)]) for b, c in zip(out.loc[idx, 'baseline'], out.loc[idx, 'conc'])]
    return out


def _base_key(key: str) -> str:
    return 'carrying' if key.startswith('without 20') or key == 't_interval' else key


def reading_responses(key: str, source: str = 'maps', three_vials: bool = False) -> Tuple[pd.DataFrame, str]:
    """The series' responses and baseline of reading `key` (READINGS_SERIES, 'without <fly>')."""
    k = _base_key(key)
    if k == 'all_presentations':
        resp = every_presentation(series_responses(source, 'all', None))
    else:
        resp = series_responses(source)
    if three_vials:
        resp = average_three_vials(resp)
    return resp, READINGS[k]['baseline']


def reading_recordings(resp: pd.DataFrame, key: str, baseline: str, plane: Optional[float] = CARRYING_PLANE) -> List[str]:
    """The runs a reading reads: those with rows of its baseline (the order-matched one needs a paraffin-oil-only block),
    at `plane` (None: all planes), without the marked runs ('without_marked') or one fly ('without <fly>')."""
    r = resp[resp['baseline'] == baseline]
    if plane is not None:
        r = r[np.isclose(r['plane_um'].astype(float), float(plane))]
    recs = sorted(r['recording'].unique())
    if key == 'without_marked':
        recs = [x for x in recs if x not in MARKED]
    if key.startswith('without 20'):
        fly = key[len('without '):]
        recs = sorted(r.loc[r['fly'] != fly, 'recording'].unique())
    return recs


def series_ratios(key: str = 'carrying', plane: Optional[float] = CARRYING_PLANE, source: str = 'maps',
                  three_vials: bool = False, recordings: Optional[Sequence[str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The series' eight ratios under reading `key` (`E98.flies_ratios`: the series' design ranges, the floor, the
    interpolation and the one-log rule; vial 1 carries): (one row per ratio, per recording rows)."""
    E98 = _e98()
    resp, bl = reading_responses(key, source, three_vials)
    recs = list(recordings) if recordings is not None else reading_recordings(resp, key, bl, plane)
    f, pr = E98.flies_ratios(resp, (bl, False), 'all', RATIOS, None, recs)
    f.insert(0, 'series_reading', key)
    pr.insert(0, 'series_reading', key)
    return f, pr


def per_fly(pr: pd.DataFrame, rid: str) -> pd.Series:
    """The flies' mean (natural) log ratio of `rid` over their entered recordings."""
    e = pr[(pr['ratio_id'] == rid) & pr['entered'].astype(bool)]
    return e.groupby('fly')['log_ratio'].mean()


def twelve_locked() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Part A's locked tables of the 12 (a_ratios_12.csv, a_per_fly_12.csv); refuses unless the lock holds."""
    _require_lock('B')
    return read_table(out_dir() / 'a_ratios_12.csv'), read_table(out_dir() / 'a_per_fly_12.csv')


def twelve_side(reading_12: str) -> Tuple[pd.DataFrame, Dict[str, pd.Series]]:
    """The 12's ratio rows (indexed by ratio) and per-fly mean log ratios under one of READINGS_12, from part A's locked
    tables."""
    r, pf = twelve_locked()
    rows = r[r['reading_12'] == reading_12].set_index('ratio_id')
    p = pf[pf['reading_12'] == reading_12]
    return rows, {rid: p[p['ratio_id'] == rid].set_index('fly')['mean_log_ratio'] for rid in RATIO_IDS}


def x1_rows(f_series: pd.DataFrame, pr_series: pd.DataFrame, r12: pd.DataFrame, pf12: Dict[str, pd.Series]) -> pd.DataFrame:
    """Per ratio of the fixed list: the held-out flies' ratio and interval, the reference set's, the split rule
    (`split_rule`: its share is the p of the rule, 'differ' below SPLIT_SHARE) and X1's direction (`x1_direction`)."""
    fs = f_series.set_index('ratio_id')
    rows = []
    for rid in RATIO_IDS:
        s, a = per_fly(pr_series, rid), pf12[rid].dropna()
        d = dict(ratio_id=rid, numerator=fs.loc[rid, 'numerator'], denominator=fs.loc[rid, 'denominator'],
                 ratio=float(fs.loc[rid, 'ratio']), lo90=float(fs.loc[rid, 'lo90']), hi90=float(fs.loc[rid, 'hi90']),
                 sd_log2_flies=float(fs.loc[rid, 'sd_log2_flies']), n_flies=int(fs.loc[rid, 'n_flies']),
                 n_recordings=int(fs.loc[rid, 'n_recordings']), num_how=fs.loc[rid, 'num_how'], den_how=fs.loc[rid, 'den_how'],
                 ratio_ref=float(r12.loc[rid, 'ratio']), lo90_ref=float(r12.loc[rid, 'lo90']),
                 hi90_ref=float(r12.loc[rid, 'hi90']), n_flies_ref=int(len(a)), set_s_lo=SET_S[rid][0],
                 set_s_hi=SET_S[rid][1], in_x1=rid in X1_RATIOS)
        if len(s) and len(a):
            sr = split_rule(a, s)
            d.update(split_share=sr['share'], n_splits=sr['n_splits'], differ=sr['differ'],
                     **x1_direction(rid, d['ratio'], d['ratio_ref']))
        else:
            d.update(split_share=np.nan, n_splits=0, differ=False, toward_s=False, nearer_s=False)
        rows.append(d)
    return pd.DataFrame(rows)


def series_flies(key: str = 'carrying', plane: Optional[float] = CARRYING_PLANE) -> List[str]:
    resp, bl = reading_responses(key)
    recs = reading_recordings(resp, key, bl, plane)
    return sorted(resp.loc[resp['recording'].isin(recs), 'fly'].unique())


def x1_readings() -> List[str]:
    """The readings of X1: the carrying one, partners (1) to (3), and (4) each of the held-out flies left out."""
    return list(READINGS_SERIES) + [f'without {f}' for f in series_flies()]


def x1_all() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """X1 under every reading at 18 um against the 12's matching reading (TWELVE_FOR; part A's locked tables): (per
    reading and ratio, the series' ratio rows, the series' per-recording rows, the result: the side per reading and the
    held side)."""
    rows, rat, per, sides = [], [], [], {}
    for key in x1_readings():
        f, pr = series_ratios(key)
        r12, pf12 = twelve_side(TWELVE_FOR[_base_key(key)])
        x = x1_rows(f, pr, r12, pf12).assign(series_reading=key, twelve_reading=TWELVE_FOR[_base_key(key)])
        sides[key] = x1_side(x)
        rows.append(x.assign(side=sides[key]))
        rat.append(f)
        per.append(pr)
    c = sides['carrying']
    held = c if all(v == c for v in sides.values()) else 'reading-dependent'
    res = dict(sides=sides, held=held, carrying=c)
    return pd.concat(rows, ignore_index=True), pd.concat(rat, ignore_index=True), pd.concat(per, ignore_index=True), res


def pooled_intervals() -> pd.DataFrame:
    """The pooled intervals of P1 and P2 over the 13 flies (the held-out flies' entered recordings and the 12's per-fly
    means from part A's locked a_per_fly_12.csv; E98's `fly_ratio`, a fly's recordings averaged on the log scale), under
    each reading of X1 with the 12's matching reading, and Student's t on the carrying
    reading's 13 per-fly log ratios."""
    E98 = _e98()
    rows = []
    for key in x1_readings():
        _, pr = series_ratios(key)
        _, pf12 = twelve_side(TWELVE_FOR[_base_key(key)])
        for rid in ('P1', 'P2'):
            twelve = pf12[rid].rename('log_ratio').rename_axis('fly').reset_index()     # one row per fly: its mean
            e = pd.concat([pr[(pr['ratio_id'] == rid) & pr['entered'].astype(bool)][['fly', 'log_ratio']], twelve])
            fr = E98.fly_ratio(e)
            rows.append(dict(reading=key, ratio_id=rid, ratio=fr['ratio'], lo90=fr['lo90'], hi90=fr['hi90'],
                             n_flies=fr['n_flies'], n_recordings=fr['n_recordings'], sd_log2_flies=fr['sd_log2_flies']))
            if key == 'carrying':
                g, lo, hi = t_interval(e.groupby('fly')['log_ratio'].mean().to_numpy())
                rows.append(dict(reading='t_interval', ratio_id=rid, ratio=g, lo90=lo, hi90=hi, n_flies=fr['n_flies'],
                                 n_recordings=fr['n_recordings'], sd_log2_flies=fr['sd_log2_flies']))
    return pd.DataFrame(rows)


def d2_15um_responses() -> pd.DataFrame:
    """E98's part D2 recordings at 15 um read by the extended reader with the carrying reading (order-matched, the
    first two presentations)."""
    E98 = _e98()
    return read_responses(list(E98.D2_PLANES[15]), 'maps', ('order_matched',), 'first', N_FIRST)


def x1_beside() -> Tuple[pd.DataFrame, dict]:
    """X1's beside sets (carrying reading): the held-out flies at 18 um against the 12's 0.5 s flies alone (E98's
    subset 'pulse_0.5s', 3 flies); the held-out flies at 15 um against E98's D2 recordings at 15 um (the 12's design,
    E98's part D rule)."""
    E98 = _e98()
    rows, sides = [], {}
    f, pr = series_ratios('carrying')
    f12, pr12 = E98.flies_ratios(a1_responses(), ('order_matched', False), 'pulse_0.5s', RATIOS)
    x = x1_rows(f, pr, f12.set_index('ratio_id'), {rid: per_fly(pr12, rid) for rid in RATIO_IDS})
    sides['18 um against the 12\'s 0.5 s flies'] = x1_side(x)
    rows.append(x.assign(set="18 um against the 12's 0.5 s flies", side=x1_side(x)))
    f15, pr15 = series_ratios('carrying', plane=15)
    rd = d2_15um_responses()
    fd, prd = E98.flies_ratios(rd, ('order_matched', False), 'all', RATIOS, None, sorted(rd['recording'].unique()),
                               E98.twelve_design())
    x = x1_rows(f15, pr15, fd.set_index('ratio_id'), {rid: per_fly(prd, rid) for rid in RATIO_IDS})
    sides["15 um against E98's D2 at 15 um"] = x1_side(x)
    rows.append(x.assign(set="15 um against E98's D2 at 15 um", side=x1_side(x)))
    return pd.concat(rows, ignore_index=True), sides


def part_b() -> dict:
    """Part B, X1 (which side the held-out flies are on): the eight ratios of the held-out flies at 18 um under the
    carrying reading and partners (1) to (4) against the 12's matching reading (part A's locked tables), the split rule
    and the side per reading, the held side; S1 and S2 as the check; beside, the 12's 0.5 s flies and the 15 um sets;
    the pooled intervals of P1 and P2 over the 13 flies only if X1 is the 12's side, held. Writes b_x1.csv,
    b_ratios.csv, b_per_recording.csv, b_responses.csv, b_beside.csv, b_pooled.csv (if formed), b_summary.json."""
    _require_lock('B')
    out = out_dir()
    x1, rat, per, res = x1_all()
    write_table(x1, out / 'b_x1.csv')
    write_table(rat, out / 'b_ratios.csv')
    write_table(per, out / 'b_per_recording.csv')
    write_table(series_responses(), out / 'b_responses.csv')
    bes, bsides = x1_beside()
    write_table(bes, out / 'b_beside.csv')
    car = x1[x1['series_reading'] == 'carrying'].set_index('ratio_id')
    check = {rid: bool(car.loc[rid, 'differ']) for rid in ('S1', 'S2')}
    pooled = None
    if res['held'] == SIDES[0]:
        pooled = pooled_intervals()
        write_table(pooled, out / 'b_pooled.csv')
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED, x1=res['held'], sides_per_reading=res['sides'],
                carrying=car[['ratio', 'lo90', 'hi90', 'n_flies', 'ratio_ref', 'lo90_ref', 'hi90_ref', 'split_share',
                              'differ', 'toward_s', 'nearer_s']].to_dict('index'),
                s1_s2_differ=check, sets_differ_in_more_than_ms=bool(any(check.values())), beside=bsides,
                pooled_formed=pooled is not None,
                pooled=pooled.to_dict('records') if pooled is not None else 'not formed: X1 is not the 12\'s side, held',
                readings=dict(series={k: READINGS[_base_key(k)] for k in x1_readings()}, twelve=TWELVE_FOR),
                recordings_carrying=reading_recordings(series_responses(), 'carrying', 'order_matched'))
    _write_json(out / 'b_summary.json', summ)
    provenance('B', x1=res['held'])
    return summ


# ---------------------------------------------------------------------------
# part C: the family against the held-out flies

def x_classes(ivs: Dict[str, Dict[str, Tuple[float, float]]], readings: Sequence[Tuple[str, str, str]],
              levels: Optional[Dict[str, pd.DataFrame]] = None, num: Optional[str] = None,
              estimates: Optional[Dict[str, str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """X2 (the hull of P1's and P2's intervals) and X3 (P1's) per estimate of ms's input and reading: `ivs` the flies'
    intervals by flies' reading key and ratio; `readings` (key, flies' key, model level), the carrying one first and
    every other a partner; the class of the eight wirings' sets (`wiring_classes`, `comparison_class`) and of their
    linear readings; the held class per estimate (`held_class`) and X2 / X3 over the estimates (`x_over_estimates`).
    Returns (one row per X, estimate and reading; per wiring rows of the carrying and settled readings; summary)."""
    E98 = _e98()
    levels = family_levels() if levels is None else levels
    num = E98.V1_NUM if num is None else num
    estimates = E98.V1_ESTIMATES if estimates is None else estimates
    rows, wrows, summ = [], [], {}
    partners = [k for k, _, _ in readings if k != 'carrying']
    for x, ids in X_IDS.items():
        held = {}
        for est, lab in estimates.items():
            cats = {}
            for key, fk, level in readings:
                lo, hi = _hull([ivs[fk][i] for i in ids])
                wc = wiring_classes(levels[level], num, lab, lo, hi)
                cats[key] = comparison_class(wc['cls_ab'])
                rows.append(dict(x=x, estimate=est, reading=key, flies_reading=fk, model_level=level, flies_lo90=lo,
                                 flies_hi90=hi, n_wirings=int(len(wc)), category=cats[key],
                                 category_linear=comparison_class(wc['cls_ab_linear']),
                                 **{f'n_{k}': int((wc['cls_ab'] == k).sum()) for k in X_CLASSES + ('mixed',)},
                                 model_linear_min=float(wc['ratio_linear'].min()),
                                 model_linear_max=float(wc['ratio_linear'].max()),
                                 model_min=float(wc['ratio_min'].min()), model_max=float(wc['ratio_max'].max())))
                if key in ('carrying', 'settled'):
                    wrows.append(wc.assign(x=x, estimate=est, reading=key, model_level=level))
            held[est] = held_class(cats, partners)
            for r in rows:
                if r['x'] == x and r['estimate'] == est:
                    r['held_per_estimate'] = held[est]
        summ[x] = dict(over_estimates=x_over_estimates(held), held_per_estimate=held,
                       carrying_per_estimate={r['estimate']: r['category'] for r in rows
                                              if r['x'] == x and r['reading'] == 'carrying'})
    cl = pd.DataFrame(rows)
    for x in X_IDS:
        cl.loc[cl['x'] == x, 'x_over_estimates'] = summ[x]['over_estimates']
    wcols = ['x', 'estimate', 'reading', 'model_level', 'wiring', 'numerator', 'denominator', 'a_num', 'a_den', 'a0',
             'ratio_linear', 'ratio_min', 'ratio_max', 'flies_lo90', 'flies_hi90', 'cls_ab', 'cls_ab_linear', 'cls']
    w = pd.concat(wrows, ignore_index=True)[wcols].rename(columns={'cls': 'cls_e98'})
    return cl, w, summ


def series_intervals(keys: Sequence[str]) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """The held-out flies' intervals at 18 um per reading key (READINGS_SERIES, 'without <fly>', and 't_interval':
    Student's t on the carrying reading's per-fly log ratios)."""
    out = {}
    for k in keys:
        if k == 't_interval':
            _, pr = series_ratios('carrying')
            out[k] = {rid: t_interval(per_fly(pr, rid).to_numpy())[1:] for rid in RATIO_IDS}
        else:
            f, _ = series_ratios(k)
            f = f.set_index('ratio_id')
            out[k] = {rid: (float(f.loc[rid, 'lo90']), float(f.loc[rid, 'hi90'])) for rid in RATIO_IDS}
    return out


def class_readings(fly_keys: Sequence[str]) -> List[Tuple[str, str, str]]:
    """The readings of the family's classes: the carrying one and partners (1) to (4) with the model at the blank's
    level, (5) Student's t, (6) the settled level with the flies' carrying reading."""
    r = [(k, k, 'blank') for k in fly_keys if k != 't_interval']
    return r + [('t_interval', 't_interval', 'blank'), ('settled', 'carrying', 'settled')]


def widen(ivs: Dict[str, Dict[str, Tuple[float, float]]], sd_log2: float) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """Every interval widened by `sd_log2` (log2) at each end: [lo 2^-sd, hi 2^sd]."""
    f = 2.0 ** float(sd_log2)
    return {k: {rid: (lo / f, hi * f) for rid, (lo, hi) in v.items()} for k, v in ivs.items()}


def part_c() -> dict:
    """Part C, the family against the held-out flies: X2 and X3 per estimate of ms's input under the carrying reading
    and partners (1) to (6) (c_classes.csv, c_class_wirings.csv); again against the pooled intervals where part B forms
    them (c_pooled_classes.csv); described, each class with the intervals widened by Y1's persistent SD
    (c_widened_classes.csv). c_summary.json."""
    _require_lock('C')
    out = out_dir()
    keys = x1_readings() + ['t_interval']
    ivs = series_intervals(keys)
    readings = class_readings(keys)
    cl, w, summ = x_classes(ivs, readings)
    write_table(cl, out / 'c_classes.csv')
    write_table(w, out / 'c_class_wirings.csv')
    _, _, _, x1 = x1_all()
    pooled = None
    if x1['held'] == SIDES[0]:
        p = pooled_intervals()
        piv = {k: {rid: (float(g.set_index('ratio_id').loc[rid, 'lo90']), float(g.set_index('ratio_id').loc[rid, 'hi90']))
                   for rid in ('P1', 'P2')} for k, g in p.groupby('reading')}
        pc, _, pooled = x_classes(piv, class_readings(list(piv)))
        write_table(pc, out / 'c_pooled_classes.csv')
    y1 = y1_summary(y1_table())
    sd = y1['carrying']['all']['persistent_sd']
    wc, _, wsumm = x_classes(widen(ivs, sd), readings)
    write_table(wc.assign(widened_by_log2=sd), out / 'c_widened_classes.csv')
    changed = {x: dict(before=summ[x]['over_estimates'], after=wsumm[x]['over_estimates'],
                       held_before=summ[x]['held_per_estimate'], held_after=wsumm[x]['held_per_estimate'],
                       changes=bool(summ[x]['held_per_estimate'] != wsumm[x]['held_per_estimate']))
               for x in X_IDS}
    res = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED,
               classes={x: dict(over_estimates=summ[x]['over_estimates'], held_per_estimate=summ[x]['held_per_estimate'],
                                carrying_per_estimate=summ[x]['carrying_per_estimate']) for x in X_IDS},
               pooled=({x: dict(over_estimates=pooled[x]['over_estimates'], held_per_estimate=pooled[x]['held_per_estimate'])
                        for x in X_IDS} if pooled else "not formed: X1 is not the 12's side, held"),
               widened=dict(persistent_sd_log2=sd, changes=changed,
                            note=("each interval widened by Y1's persistent SD (log2) at each end, as if one pair of vials "
                                  "had served all five days; described")),
               readings=[dict(key=k, flies=fk, model_level=lv) for k, fk, lv in readings],
               family_sources=family_levels()['blank'].attrs.get('sources'), x1=x1['held'])
    _write_json(out / 'c_summary.json', res)
    provenance('C', classes={x: res['classes'][x]['over_estimates'] for x in X_IDS})
    return res


# ---------------------------------------------------------------------------
# part D: the vials

def y1_rows(resp: pd.DataFrame, baseline: str) -> pd.DataFrame:
    """Y1's rows: per run, concentration and pair of vials of the run's three-vial odor (VIAL_PAIRS), the log2 ratio
    of the lower-numbered vial's response over the other's (both above MIN_RESPONSE; each at its first kept block) and
    their block distance; 'carries' False for the pairs with vial 1 in the five runs of 10-22 with three vials of paa
    (vial 1 in the other half of the run: vials 2 and 3 carry there)."""
    E98 = _e98()
    r = resp[resp['baseline'] == baseline]
    rows = []
    for rec, g in r.groupby('recording'):
        o = three_vial_odor(g['odor'])
        if o is None:
            continue
        lab = {1: o, 2: f'{o}_vial2', 3: f'{o}_vial3'}
        gi = g.set_index(['odor', 'conc'])
        for c in sorted(g.loc[g['odor'] == o, 'conc'].unique()):
            for a, b in VIAL_PAIRS:
                ka, kb = (lab[a], c), (lab[b], c)
                if ka not in gi.index or kb not in gi.index:
                    continue
                va, vb = float(gi.loc[ka, 'response']), float(gi.loc[kb, 'response'])
                ok = va > E98.MIN_RESPONSE and vb > E98.MIN_RESPONSE
                rows.append(dict(recording=rec, day=rec[:8], fly=g['fly'].iloc[0], side=g['side'].iloc[0],
                                 plane_um=float(g['plane_um'].iloc[0]), odor=o, conc=float(c), pair=f'{a}/{b}',
                                 response_a=va, response_b=vb, block_a=int(gi.loc[ka, 'block1']),
                                 block_b=int(gi.loc[kb, 'block1']),
                                 block_distance=abs(int(gi.loc[ka, 'block1']) - int(gi.loc[kb, 'block1'])),
                                 log2_ratio=float(np.log2(va / vb)) if ok else np.nan, entered=bool(ok),
                                 carries=not (rec in TEN22_PAA and a == 1), marked=rec in MARKED, baseline=baseline))
    return pd.DataFrame(rows)


def y1_table() -> pd.DataFrame:
    """Y1's rows under the carrying reading (order-matched, the first two presentations; all planes), the first-in-block
    baseline (which admits the four runs without a paraffin-oil-only block) and without the marked runs."""
    resp = series_responses()
    out = []
    for key in Y1_READINGS:
        bl = READINGS[key]['baseline']
        t = y1_rows(resp, bl).assign(reading=key)
        if key == 'without_marked':
            t = t[~t['marked']]
        out.append(t)
    return pd.concat(out, ignore_index=True)


def y1_summary(t: pd.DataFrame) -> dict:
    """Y1 per reading over the carrying rows (entered): the variance split (`y1_variance_split`) over both odors and per
    odor, the yardstick (`vial_yardstick`), and the SD by block distance."""
    out = {}
    for key, g in t.groupby('reading'):
        g = g[g['carries'].astype(bool) & g['entered'].astype(bool)]
        d = {'all': y1_variance_split(g)}
        for o, go in g.groupby('odor'):
            d[o] = y1_variance_split(go)
        for k in list(d):
            d[k]['yardstick_log2'] = vial_yardstick(d[k]['persistent_sd']) if np.isfinite(d[k]['persistent_sd']) else np.nan
            d[k]['yardstick_factor'] = float(2 ** d[k]['yardstick_log2'])
        d['by_block_distance'] = {int(b): dict(n=int(len(x)), mean=float(x['log2_ratio'].mean()),
                                               sd=float(x['log2_ratio'].std(ddof=1)) if len(x) > 1 else np.nan,
                                               rms=float(np.sqrt((x['log2_ratio'] ** 2).mean())))
                                  for b, x in g.groupby('block_distance')}
        out[key] = d
    return out


def series_maps(row) -> Tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """One run's peak-map stack, stimulus table and calyx outline mask (`require_part_a` first; the frame count must
    equal the presentations)."""
    import tifffile
    from experiments import E100_apl_weak_odors as E100
    require_part_a(row['stem'])
    check_row(row)
    E100._ds()
    from data_prep import apl_calyx_maps as ACM
    p = E100.recording_paths(row)
    stack = tifffile.imread(p['peak_maps']).astype(float)
    mask = ACM.roi_mask(ACM.read_imagej_polygon(p['roi']), stack.shape[1:])
    table = stim_table(row)
    if stack.shape[0] != len(table):
        raise RuntimeError(f"{row['stem']}: {stack.shape[0]} maps for {len(table)} presentations")
    return stack, table, mask


def block_maps(stack: np.ndarray, table: pd.DataFrame) -> Optional[Dict[Tuple[str, float], dict]]:
    """Per stimulus at its first kept block, its presentation maps minus the order-matched baseline map (the mean of
    the paraffin-oil-only blocks' maps at the same position), every presentation; None without such a block."""
    E98 = _e98()
    t = table.reset_index(drop=True)
    kept = t[~t['out']]
    pb = kept[kept['pfo_block']]
    if pb.empty:
        return None
    base = {int(p): stack[g.index].mean(axis=0) for p, g in pb.groupby('position')}
    out = {}
    for _, sb in E98.d_stimulus_blocks(t, 'first').iterrows():
        g = t[(t['block1'] == sb['block1']) & (t['odor'] == sb['odor']) & (t['log10_conc'] == sb['log10_conc'])]
        out[(str(sb['odor']), round(float(sb['log10_conc']), 3))] = dict(
            block1=int(sb['block1']), maps=[stack[i] - base[int(p)] for i, p in zip(g.index, g['position'])])
    return out


def _r(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def m_recording(bm: Dict[Tuple[str, float], dict], mask: np.ndarray, blur: float = M_BLUR_PX,
                remove_common: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """M for one run: the common map from the nine stimuli of ms and vial 1 of paa and pa (E79's `svd_shared` on the
    blurred block-mean maps inside the outline, means removed); each stimulus's remainder (E79's `_residual`; with
    `remove_common` False the blurred map with its mean removed); r_rep per stimulus (the remainders of the first and
    second halves of its presentations in its block); the pairs: vials 2 and 3 of the three-vial odor at one
    concentration ('vial') and vial 2 or 3 against each other odor's vial 1 at the same concentration rank ('odor'),
    with r between the remainders and the block distance. Returns (per stimulus, per pair)."""
    from experiments import E100_apl_weak_odors as E100
    DS = E100._ds()
    concs = {}
    for o, c in bm:
        concs.setdefault(o, []).append(c)
    concs = {o: sorted(v) for o, v in concs.items()}
    common = [(o, c) for o in ODORS3 for c in concs.get(o, [])]
    v = None
    if remove_common:
        X = np.stack([_blurred_in(np.mean(bm[k]['maps'], axis=0), mask, blur) for k in common])
        v = DS.svd_shared(X, True)[1]

    def rem(m):
        x = _blurred_in(m, mask, blur)
        return DS._residual(x, v, True) if v is not None else x - x.mean()
    rems, rrep, stim = {}, {}, []
    for k, d in bm.items():
        p = d['maps']
        rems[k] = rem(np.mean(p, axis=0))
        h = len(p) // 2
        rrep[k] = _r(rem(np.mean(p[:h], axis=0)), rem(np.mean(p[h:], axis=0))) if h >= 1 else np.nan
        stim.append(dict(odor=k[0], conc=k[1], block1=d['block1'], n_presentations=len(p), r_rep=rrep[k]))
    o3 = three_vial_odor([o for o, _ in bm])
    pairs = []
    if o3 is not None:
        base = concs.get(o3, [])
        for rank, c in enumerate(base):
            a, b = (f'{o3}_vial2', c), (f'{o3}_vial3', c)
            if a in bm and b in bm:
                pairs.append(dict(kind='vial', a=f'{a[0]} {c:g}', b=f'{b[0]} {c:g}', rank=rank + 1, r=_r(rems[a], rems[b]),
                                  block_distance=abs(bm[a]['block1'] - bm[b]['block1']), r_rep_a=rrep[a], r_rep_b=rrep[b]))
            for k in (a, b):
                if k not in bm:
                    continue
                for other in (o for o in ODORS3 if o != o3):
                    if rank >= len(concs.get(other, [])):
                        continue
                    q = (other, concs[other][rank])
                    pairs.append(dict(kind='odor', a=f'{k[0]} {c:g}', b=f'{q[0]} {q[1]:g}', rank=rank + 1,
                                      r=_r(rems[k], rems[q]), block_distance=abs(bm[k]['block1'] - bm[q]['block1']),
                                      r_rep_a=rrep[k], r_rep_b=rrep[q]))
    return pd.DataFrame(stim), pd.DataFrame(pairs)


def m_calyces(pairs: pd.DataFrame, max_blocks: Optional[int] = M_MAX_BLOCKS_APART,
              max_blocks_odor: Optional[int] = M_MAX_BLOCKS_APART, min_r_rep: float = M_MIN_R_REP) -> pd.DataFrame:
    """Per calyx (fly, side): the pairs that enter (block distance at most `max_blocks` for vial pairs and
    `max_blocks_odor` for odor pairs, None for no limit; both stimuli's r_rep at least `min_r_rep`), r_vial and r_odor
    their means, r_rep the mean r_rep of the entering vial pairs' stimuli, and the share (`m_share`)."""
    p = pairs.copy()
    lim = np.where(p['kind'] == 'vial', np.inf if max_blocks is None else max_blocks,
                   np.inf if max_blocks_odor is None else max_blocks_odor)
    p['enters'] = (p['block_distance'] <= lim) & (p['r_rep_a'] >= min_r_rep) & (p['r_rep_b'] >= min_r_rep)
    rows = []
    for (fly, side), g in p.groupby(['fly', 'side']):
        e = g[g['enters']]
        v, o = e[e['kind'] == 'vial'], e[e['kind'] == 'odor']
        rv = float(v['r'].mean()) if len(v) else np.nan
        ro = float(o['r'].mean()) if len(o) else np.nan
        rr = float(np.mean(np.concatenate([v['r_rep_a'], v['r_rep_b']]))) if len(v) else np.nan
        sh = m_share(rv, ro, rr) if np.isfinite(rv) and np.isfinite(ro) and np.isfinite(rr) else np.nan
        rows.append(dict(fly=fly, side=side, calyx=f'{fly}_{side}', n_vial_pairs=int(len(v)), n_odor_pairs=int(len(o)),
                         n_recordings=int(e['recording'].nunique()), r_vial=rv, r_odor=ro, r_rep=rr, share=sh,
                         vial_over_odor=bool(np.isfinite(rv) and np.isfinite(ro) and rv > ro)))
    return pd.DataFrame(rows)


def m_result(cal: pd.DataFrame) -> dict:
    """M's class from the calyces: the share's mean over flies (a fly's calyces averaged) and its 90 % interval from
    E79's fly bootstrap (arithmetic, 10,000 draws), r_vial > r_odor counted over calyces, `m_class`."""
    E98 = _e98()
    from experiments import E79_apl_calyx_constraints as E79
    c = cal[np.isfinite(cal['share'])]
    pf = c.groupby('fly')['share'].mean()
    if len(pf) == 0:
        return dict(n_calyces=0, n_flies=0, share=np.nan, lo90=np.nan, hi90=np.nan, n_vial_over_odor=0, cls='not decided')
    m, lo, hi = E79.fly_bootstrap(pf.to_numpy(), n_boot=E98.N_BOOT, seed=E98.SEED)
    n = int(c['vial_over_odor'].sum())
    return dict(n_calyces=int(len(c)), n_flies=int(len(pf)), share=float(m[0]), lo90=float(lo[0]), hi90=float(hi[0]),
                n_vial_over_odor=n, cls=m_class(float(lo[0]), float(hi[0]), n))


def m_tables() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """M over every run with a paraffin-oil-only block (all planes), at each blur of M_BLURS_PX, with and without the
    common map removed: (per stimulus, per pair), each row with its run, fly, side, plane, blur and variant."""
    rows = series_rows(list(SERIES))
    stims, pairs = [], []
    for st, row in rows.iterrows():
        stack, table, mask = series_maps(row)
        bm = block_maps(stack, table)
        if bm is None:
            continue
        meta = dict(recording=st, fly=row['fly_id'], side=row['side'], plane_um=float(row['plane_depth_um']),
                    ten22_pa=st in TEN22_PA)
        for blur in M_BLURS_PX:
            for rc in (True, False):
                s, p = m_recording(bm, mask, blur, rc)
                kw = dict(meta, blur_px=blur, common_removed=rc)
                stims.append(s.assign(**kw))
                pairs.append(p.assign(**kw))
    return pd.concat(stims, ignore_index=True), pd.concat(pairs, ignore_index=True)


def m_summary(pairs: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per blur and variant (the carrying one: the common map removed, without the five 10-22 runs with three vials of
    pa; before the common map is removed; those five runs alone with no block limit on the odor pairs, beside): the
    calyces' table and M's result."""
    cal_rows, res_rows = [], []
    for blur, g in pairs.groupby('blur_px'):
        variants = {'carrying': (g[g['common_removed'] & ~g['ten22_pa']], M_MAX_BLOCKS_APART),
                    'before_common_removed': (g[~g['common_removed'] & ~g['ten22_pa']], M_MAX_BLOCKS_APART),
                    'ten22_pa_beside': (g[g['common_removed'] & g['ten22_pa']], None)}
        for name, (p, odor_lim) in variants.items():
            cal = m_calyces(p, M_MAX_BLOCKS_APART, odor_lim) if len(p) else pd.DataFrame()
            res = m_result(cal) if len(cal) else dict(n_calyces=0, cls='not decided')
            cal_rows.append(cal.assign(blur_px=blur, variant=name))
            res_rows.append(dict(blur_px=blur, variant=name, carries=bool(name == 'carrying' and blur == M_BLUR_PX), **res))
    return pd.concat(cal_rows, ignore_index=True), pd.DataFrame(res_rows)


def part_d() -> dict:
    """Part D, the vials: Y1 (d_y1_pairs.csv, the SDs, the persistent part and the yardstick in d_summary.json), M
    (d_m_stimuli.csv, d_m_pairs.csv, d_m_calyces.csv, d_m.csv: the share and class per blur and variant)."""
    _require_lock('D')
    out = out_dir()
    t = y1_table()
    write_table(t, out / 'd_y1_pairs.csv')
    y1 = y1_summary(t)
    stims, pairs = m_tables()
    write_table(stims, out / 'd_m_stimuli.csv')
    write_table(pairs, out / 'd_m_pairs.csv')
    cal, res = m_summary(pairs)
    write_table(cal, out / 'd_m_calyces.csv')
    write_table(res, out / 'd_m.csv')
    car = res[res['carries']].iloc[0].to_dict()
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED, y1=y1, m_carrying=car,
                m=res.to_dict('records'),
                rules=dict(y1_yardstick=Y1_YARDSTICK, m_blur_carrying=M_BLUR_PX, m_max_blocks=M_MAX_BLOCKS_APART,
                           m_min_r_rep=M_MIN_R_REP, m_odor_lo=M_ODOR_LO, m_min_calyces=M_MIN_CALYCES, m_block_hi=M_BLOCK_HI),
                notes=("Y1: every pair of the three vials at each concentration, each vial at its first kept block, in "
                       "every run with a reading (all planes); in the five 10-22 runs with three vials of paa only vials 2 "
                       "and 3 carry. M: every run with a paraffin-oil-only block, all planes, every presentation of the "
                       "stimulus's first block; the five 10-22 runs with three vials of pa are beside."))
    _write_json(out / 'd_summary.json', summ)
    provenance('D', y1_persistent_sd=y1['carrying']['all']['persistent_sd'], m=car.get('cls'))
    return summ


# ---------------------------------------------------------------------------
# part E: beside (Y2, Y3, Z, Hallem's input, all three vials, Yang's tables)

def y2_repeats() -> pd.DataFrame:
    """Y2: every stimulus read in two or more kept blocks (the two passes, 16 blocks apart; pa vial 1's two series in the
    runs of 10-22), each block under the carrying reading: the later block's response over the first's (natural log),
    with the block distance and the kind."""
    E98 = _e98()
    r = series_responses('maps', 'all', N_FIRST)
    r = r[r['baseline'] == 'order_matched']
    rows = []
    for (rec, o, c), g in r.groupby(['recording', 'odor', 'conc']):
        g = g.sort_values('block1')
        if len(g) < 2:
            continue
        first = g.iloc[0]
        for _, x in g.iloc[1:].iterrows():
            ok = first['response'] > E98.MIN_RESPONSE and x['response'] > E98.MIN_RESPONSE
            rows.append(dict(recording=rec, fly=first['fly'], side=first['side'], plane_um=float(first['plane_um']),
                             odor=o, conc=float(c), first_block=int(first['block1']), later_block=int(x['block1']),
                             block_distance=int(x['block1'] - first['block1']), first=float(first['response']),
                             later=float(x['response']),
                             log_later_over_first=float(np.log(x['response'] / first['response'])) if ok else np.nan,
                             kind='two passes' if rec in TWO_PASSES else '10-22 pa vial 1' if rec.startswith('20211022')
                             else 'other'))
    return pd.DataFrame(rows)


def drift_per_block(rep: pd.DataFrame) -> float:
    """The drift per block (natural log): the least-squares slope through the origin of the later-over-first log ratio
    on the block distance."""
    d = rep.dropna(subset=['log_later_over_first'])
    x, y = d['block_distance'].to_numpy(float), d['log_later_over_first'].to_numpy(float)
    return float((x * y).sum() / (x * x).sum()) if len(d) else np.nan


def y2_p1_corrected(drift: float) -> Tuple[pd.DataFrame, dict]:
    """P1 at 18 um (carrying reading) with its block lag corrected by the drift: per recording the lag (E98's
    `block_lag`, paa at 10^-4 over ms at 10^-3) and log P1 minus drift x lag; the flies' ratio before and after."""
    E98 = _e98()
    _, pr = series_ratios('carrying')
    e = pr[(pr['ratio_id'] == 'P1') & pr['entered'].astype(bool)].copy()
    e['lag'] = [E98.block_lag(series_table(rec), ('paa', -4.0), ('ms', -3.0)) for rec in e['recording']]
    e['log_ratio_corrected'] = e['log_ratio'] - drift * e['lag']
    before = E98.fly_ratio(e[['fly', 'log_ratio']])
    after = E98.fly_ratio(e[['fly']].assign(log_ratio=e['log_ratio_corrected']))
    return e, dict(drift_per_block=drift, before=before, after=after)


def y3_within_block() -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Y3: per run (all planes, carrying baseline) and stimulus at its first block, the k-th presentation over the
    first (k 2 to 4; both above MIN_RESPONSE; natural log); by odor, concentration and k the count, median and
    geometric mean; and whether paa and ms differ in the second over the first: per fly the mean log ratio of paa's
    (vial 1) minus ms's, its mean over flies and E79's fly bootstrap interval (differ when it excludes 0)."""
    E98 = _e98()
    from experiments import E79_apl_calyx_constraints as E79
    r = series_responses()
    r = r[r['baseline'] == 'order_matched']
    rows = []
    for _, x in r.iterrows():
        for k in range(2, Y3_MAX_REP + 1):
            col = f'rep{k}'
            if col not in r.columns or not np.isfinite(x.get(col, np.nan)):
                continue
            ok = x['rep1'] > E98.MIN_RESPONSE and x[col] > E98.MIN_RESPONSE
            rows.append(dict(recording=x['recording'], fly=x['fly'], plane_um=float(x['plane_um']), odor=x['odor'],
                             conc=float(x['conc']), k=k, first=float(x['rep1']), kth=float(x[col]),
                             log_kth_over_first=float(np.log(x[col] / x['rep1'])) if ok else np.nan))
    t = pd.DataFrame(rows)
    s = t.dropna(subset=['log_kth_over_first']).groupby(['odor', 'conc', 'k'])['log_kth_over_first'].agg(
        n='size', median='median', mean='mean').reset_index()
    s['geomean_ratio'] = np.exp(s['mean'])
    k2 = t[(t['k'] == 2)].dropna(subset=['log_kth_over_first'])
    pf = k2[k2['odor'].isin(['paa', 'ms'])].groupby(['fly', 'odor'])['log_kth_over_first'].mean().unstack()
    pf = pf.dropna()
    d = (pf['paa'] - pf['ms']) if len(pf) else pd.Series(dtype=float)
    res = dict(n_flies=int(len(d)))
    if len(d):
        m, lo, hi = E79.fly_bootstrap(d.to_numpy(), n_boot=E98.N_BOOT, seed=E98.SEED)
        res.update(mean_log_diff=float(m[0]), lo90=float(lo[0]), hi90=float(hi[0]),
                   factor=float(np.exp(m[0])), differ=bool(lo[0] > 0 or hi[0] < 0), per_fly=d.to_dict())
    return t, s, res


def calyx_runs() -> pd.DataFrame:
    """The series' processed runs per calyx in run order, with their planes."""
    rows = series_rows(list(SERIES))
    t = pd.DataFrame(dict(recording=rows.index, fly=rows['fly_id'].to_numpy(), side=rows['side'].to_numpy(),
                          plane_um=rows['plane_depth_um'].astype(float).to_numpy(),
                          run=[int(re.search(r'fun(\d+)', s).group(1)) for s in rows.index]))
    return t.sort_values(['fly', 'side', 'run']).reset_index(drop=True)


def alternating_calyces(runs: pd.DataFrame) -> List[str]:
    """The calyces whose planes alternate from run to run: at least two switches between 15 and 18 um among consecutive
    processed runs."""
    out = []
    for (fly, side), g in runs.groupby(['fly', 'side']):
        p = g['plane_um'].tolist()
        n = sum({a, b} == {15.0, 18.0} for a, b in zip(p[:-1], p[1:]))
        if n >= 2:
            out.append(f'{fly}_{side}')
    return out


def z_tables(rid: str = 'P1') -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Z for ratio `rid` (P1 carries, P3 beside), carrying reading, all planes: per calyx the mean log ratio at 15 um
    minus at 18 um (`z_rule` over the calyces, the mean factor's interval from E79's fly bootstrap of the flies' mean
    calyx differences); in the alternating calyces the difference between neighbouring runs of the two planes, by order."""
    E98 = _e98()
    from experiments import E79_apl_calyx_constraints as E79
    _, pr = series_ratios('carrying', plane=None)
    runs = calyx_runs()
    meta = runs.set_index('recording')
    e = pr[(pr['ratio_id'] == rid) & pr['entered'].astype(bool)].copy()
    e['side'] = e['recording'].map(meta['side'])
    e['plane_um'] = e['recording'].map(meta['plane_um'])
    rows = []
    for (fly, side), g in e.groupby(['fly', 'side']):
        a, b = g[g['plane_um'] == 15.0]['log_ratio'], g[g['plane_um'] == 18.0]['log_ratio']
        rows.append(dict(ratio_id=rid, fly=fly, side=side, calyx=f'{fly}_{side}', n_15=int(len(a)), n_18=int(len(b)),
                         log_15=float(a.mean()) if len(a) else np.nan, log_18=float(b.mean()) if len(b) else np.nan,
                         log_diff=float(a.mean() - b.mean()) if len(a) and len(b) else np.nan))
    cal = pd.DataFrame(rows)
    c = cal.dropna(subset=['log_diff'])
    pf = c.groupby('fly')['log_diff'].mean()
    m, lo, hi = E79.fly_bootstrap(pf.to_numpy(), n_boot=E98.N_BOOT, seed=E98.SEED)
    res = dict(ratio_id=rid, n_calyces=int(len(c)), n_higher_15=int((c['log_diff'] > 0).sum()), n_flies=int(len(pf)),
               factor=float(np.exp(m[0])), lo90=float(np.exp(lo[0])), hi90=float(np.exp(hi[0])),
               cls=z_rule(c['log_diff'].to_numpy(), float(np.exp(lo[0])), float(np.exp(hi[0]))))
    alt = alternating_calyces(runs)
    lr = e.set_index('recording')['log_ratio']
    nb = []
    for (fly, side), g in runs.groupby(['fly', 'side']):
        if f'{fly}_{side}' not in alt:
            continue
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            if {a['plane_um'], b['plane_um']} != {15.0, 18.0}:
                continue
            r15, r18 = (a, b) if a['plane_um'] == 15.0 else (b, a)
            d = (lr.get(r15['recording'], np.nan) - lr.get(r18['recording'], np.nan))
            nb.append(dict(ratio_id=rid, calyx=f'{fly}_{side}', earlier=a['recording'], later=b['recording'],
                           order='15 then 18' if a['plane_um'] == 15.0 else '18 then 15', log_diff_15_minus_18=float(d)))
    nbt = pd.DataFrame(nb)
    res['alternating_calyces'] = alt
    if len(nbt):
        res['neighbours_by_order'] = {o: dict(n=int(g['log_diff_15_minus_18'].notna().sum()),
                                              mean_factor=float(np.exp(g['log_diff_15_minus_18'].mean())),
                                              n_higher_15=int((g['log_diff_15_minus_18'] > 0).sum()))
                                      for o, g in nbt.groupby('order')}
    return cal, nbt, res


def hallem_rows_e(intervals: Dict[str, Dict[str, Tuple[float, float]]]) -> pd.DataFrame:
    """Q2 and F1 against the family's ratios on Hallem's rows (E98's B1 labels; as E98's V5 at 10^-4 and V3 at 10^-4),
    at the blank's level, per set of flies (`intervals`: set -> ratio -> interval). Described."""
    E98 = _e98()
    lv = family_levels()['blank']
    pairs = {'Q2': (E98.hallem_label('pa', -4), E98.hallem_label('ms', -4)),
             'F1': (E98.hallem_label('ms', -4), E98.hallem_label('ms', -2))}
    rows = []
    for name, iv in intervals.items():
        for rid, (num, den) in pairs.items():
            lo, hi = iv[rid]
            wc = wiring_classes(lv, num, den, lo, hi)
            rows.append(dict(set=name, ratio_id=rid, numerator=num, denominator=den, flies_lo90=lo, flies_hi90=hi,
                             n_wirings=int(len(wc)), category=comparison_class(wc['cls_ab']),
                             category_linear=comparison_class(wc['cls_ab_linear']),
                             model_linear_min=float(wc['ratio_linear'].min()) if len(wc) else np.nan,
                             model_linear_max=float(wc['ratio_linear'].max()) if len(wc) else np.nan,
                             **{f'n_{k}': int((wc['cls_ab'] == k).sum()) for k in X_CLASSES + ('mixed',)}))
    return pd.DataFrame(rows)


def alt_reading_x1(name: str, source: str = 'maps', three_vials: bool = False) -> Tuple[pd.DataFrame, str]:
    """X1's rows under the carrying reading with all three vials averaged in place of vial 1, or from Yang's tables in
    place of the maps (described), against the 12's locked carrying reading."""
    f, pr = series_ratios('carrying', source=source, three_vials=three_vials)
    r12, pf12 = twelve_side('carrying')
    x = x1_rows(f, pr, r12, pf12)
    return x.assign(reading=name, side=x1_side(x)), x1_side(x)


def maps_vs_yang() -> dict:
    """The series' calyx means from the maps against Yang's roi_peak_dff tables (every presentation of every run)."""
    rows = series_rows(list(SERIES))
    d = []
    for st, row in rows.iterrows():
        a, b = frame_means(row, 'maps'), frame_means(row, 'yang')
        d.append(dict(recording=st, max_abs_diff=float(np.max(np.abs(a - b))), r=float(np.corrcoef(a, b)[0, 1])))
    t = pd.DataFrame(d)
    return dict(max_abs_diff=float(t['max_abs_diff'].max()), min_r=float(t['r'].min()))


def part_e() -> dict:
    """Part E, beside (only Z has a class): Y2 (e_y2_repeats.csv, e_y2_p1.csv), Y3 (e_y3.csv, e_y3_summary.csv), Z
    (e_z_calyces.csv, e_z_neighbours.csv, P1 carrying and P3 beside), Hallem's input (e_hallem.csv), all three vials
    averaged and Yang's tables in place of the maps (e_alt_readings.csv); e_summary.json."""
    _require_lock('E')
    out = out_dir()
    rep = y2_repeats()
    write_table(rep, out / 'e_y2_repeats.csv')
    drift = drift_per_block(rep)
    drift_kind = {k: drift_per_block(g) for k, g in rep.groupby('kind')}
    p1c, p1res = y2_p1_corrected(drift)
    write_table(p1c, out / 'e_y2_p1.csv')
    y3, y3s, y3res = y3_within_block()
    write_table(y3, out / 'e_y3.csv')
    write_table(y3s, out / 'e_y3_summary.csv')
    zc, zn, zres = [], [], {}
    for rid in ('P1', 'P3'):
        c, n, r = z_tables(rid)
        zc.append(c)
        zn.append(n)
        zres[rid] = r
    write_table(pd.concat(zc, ignore_index=True), out / 'e_z_calyces.csv')
    write_table(pd.concat(zn, ignore_index=True), out / 'e_z_neighbours.csv')
    f, _ = series_ratios('carrying')
    f = f.set_index('ratio_id')
    r12, _ = twelve_side('carrying')
    ivs = {'held-out flies (18 um, carrying)': {rid: (float(f.loc[rid, 'lo90']), float(f.loc[rid, 'hi90'])) for rid in ('Q2', 'F1')},
           'the 12 (carrying)': {rid: (float(r12.loc[rid, 'lo90']), float(r12.loc[rid, 'hi90'])) for rid in ('Q2', 'F1')}}
    hal = hallem_rows_e(ivs)
    write_table(hal, out / 'e_hallem.csv')
    alt, alt_sides = [], {}
    for name, src, tv in (('all three vials averaged', 'maps', True), ("Yang's tables", 'yang', False)):
        x, s = alt_reading_x1(name, src, tv)
        alt.append(x)
        alt_sides[name] = s
    write_table(pd.concat(alt, ignore_index=True), out / 'e_alt_readings.csv')
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED,
                y2=dict(drift_per_block=drift, drift_factor_per_block=float(np.exp(drift)), drift_by_kind=drift_kind,
                        n_pairs=int(rep['log_later_over_first'].notna().sum()), p1=p1res),
                y3=y3res, z=zres, hallem=hal.to_dict('records'), alt_readings_x1=alt_sides,
                maps_vs_yang=maps_vs_yang())
    _write_json(out / 'e_summary.json', summ)
    provenance('E', z=zres['P1']['cls'])
    return summ


# ---------------------------------------------------------------------------
# the note's figure (no simulation; from the parts' saved tables)

SET_COLORS = {'the 12': '#2a78d6', 'held-out': '#eb6834', 'set S': '#1baf7a'}   # the reference palette's first three
RATIO_LABELS = {'P1': 'P1  paa 10^-4 / ms 10^-3', 'P2': 'P2  paa 10^-3.431 / ms 10^-3',
                'P3': 'P3  paa / ms, both 10^-3.431', 'Q1': 'Q1  pa 10^-4.386 / ms 10^-3',
                'Q2': 'Q2  pa / ms, both 10^-4', 'F1': 'F1  ms 10^-4 / ms 10^-2',
                'S1': 'S1  paa 10^-4.386 / paa 10^-3.431', 'S2': 'S2  pa 10^-5.34 / pa 10^-4.386'}


def figure_ratios(path: Optional[Path] = None) -> Path:
    """E103_ratios.png: per ratio of the fixed list, the 12 and the held-out flies (geometric mean over flies, 90 %
    interval; carrying reading, 18 um), set S's value (the entry's estimate; P2 a range), and the family's linear
    range (C2 on the eight wirings at the blank's level: over the five estimates of ms's input for P1 and P2, on
    Hallem's input for Q2 and F1). From b_x1.csv, c_classes.csv and e_hallem.csv."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = out_dir()
    path = out / 'E103_ratios.png' if path is None else Path(path)
    x = read_table(out / 'b_x1.csv')
    x = x[x['series_reading'] == 'carrying'].set_index('ratio_id')
    c = read_table(out / 'c_classes.csv')
    cc = c[(c['reading'] == 'carrying') & (c['x'] == 'X2')]
    fam = {rid: (float(cc['model_linear_min'].min()), float(cc['model_linear_max'].max())) for rid in ('P1', 'P2')}
    h = read_table(out / 'e_hallem.csv')
    h = h[h['set'].str.startswith('held-out')].set_index('ratio_id')
    fam.update({rid: (float(h.loc[rid, 'model_linear_min']), float(h.loc[rid, 'model_linear_max'])) for rid in ('Q2', 'F1')})
    fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=150)
    ids = list(RATIO_IDS)
    for i, rid in enumerate(ids):
        y = len(ids) - 1 - i
        if rid in fam:
            lo, hi = fam[rid]
            ax.fill_betweenx([y - 0.38, y + 0.38], lo, hi, color='#d7d7d7', lw=0, zorder=0,
                             label='the family, linear (C2, 8 wirings)' if i == 0 else None)
        for k, (name, r, lo, hi, dy) in enumerate((
                ('the 12', x.loc[rid, 'ratio_ref'], x.loc[rid, 'lo90_ref'], x.loc[rid, 'hi90_ref'], 0.2),
                ('held-out', x.loc[rid, 'ratio'], x.loc[rid, 'lo90'], x.loc[rid, 'hi90'], 0.0))):
            ax.plot([lo, hi], [y + dy, y + dy], color=SET_COLORS[name], lw=2, solid_capstyle='round', zorder=2)
            ax.plot([r], [y + dy], 'o', ms=6.5, color=SET_COLORS[name], mec='white', mew=1.2, zorder=3,
                    label={'the 12': 'the 12 recordings, Jul to Sep 2021 (8 flies)',
                           'held-out': 'the held-out series, Oct to Nov 2021 (5 flies)'}[name] if i == 0 else None)
        s_lo, s_hi = SET_S[rid]
        if s_hi > s_lo:
            ax.plot([s_lo, s_hi], [y - 0.2, y - 0.2], color=SET_COLORS['set S'], lw=2, zorder=2)
        ax.plot([np.sqrt(s_lo * s_hi)], [y - 0.2], 'D', ms=6, color=SET_COLORS['set S'], mec='white', mew=1.2, zorder=3,
                label='set S, the 2022 single-colour flies (estimate)' if i == 0 else None)
        if bool(x.loc[rid, 'differ']):
            ax.text(0.105, y, 'differs', fontsize=7.5, color='#333333', va='center', ha='left')
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels([RATIO_LABELS[r] for r in ids[::-1]], fontsize=8.5)
    ax.set_xscale('log')
    ax.set_xlim(0.1, 8.0)
    ticks = (0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:g}' for t in ticks])
    from matplotlib.ticker import NullLocator
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1.0, color='#999999', lw=0.8, zorder=0)
    ax.set_xlabel("ratio of APL's response (calyx-mean peak dF/F), 90 % interval over flies", fontsize=9)
    for sname in ('top', 'right'):
        ax.spines[sname].set_visible(False)
    ax.tick_params(colors='#333333', labelsize=8.5)
    ax.legend(loc='upper center', bbox_to_anchor=(0.45, 1.16), ncol=2, fontsize=7.8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, facecolor='white')
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'E', 'figures'))
    a = ap.parse_args()
    if a.part == 'figures':
        print(figure_ratios())
        return
    res = {'0': part_0, 'A': part_a, 'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e}[a.part]()
    print(json.dumps(res if a.part != 'A' else dict(tables=res['tables'], dest=res['dest']), indent=1,
                     default=_json_default)[:20000])


if __name__ == '__main__':
    main()
