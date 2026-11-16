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


def _later_part(name: str) -> dict:
    if not part_a_locked():
        raise PermissionError(f"part {name} opens the series' maps: it runs only after part A's lock is committed")
    raise NotImplementedError(f'part {name} is not implemented yet (code to come, after part A)')


def part_b() -> dict:
    """Part B, X1: which side the held-out flies are on (not implemented yet)."""
    return _later_part('B')


def part_c() -> dict:
    """Part C, the family against the held-out flies (X2, X3; not implemented yet)."""
    return _later_part('C')


def part_d() -> dict:
    """Part D, the vials (Y1, M; not implemented yet)."""
    return _later_part('D')


def part_e() -> dict:
    """Part E, beside (Y2, Y3, Z, Hallem's input; not implemented yet)."""
    return _later_part('E')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'E'))
    a = ap.parse_args()
    res = {'0': part_0, 'A': part_a, 'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e}[a.part]()
    print(json.dumps(res if a.part != 'A' else dict(tables=res['tables'], dest=res['dest']), indent=1,
                     default=_json_default)[:20000])


if __name__ == '__main__':
    main()
