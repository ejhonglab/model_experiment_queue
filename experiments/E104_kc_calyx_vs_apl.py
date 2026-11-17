#!/usr/bin/env python3
"""E104: the flies' KCs beside their APL along concentration. Whether the steep fall of the flies' APL from paa at 10^-4
to 10^-5 is in its claws' input or is APL's own (part A: Yang's KC-only calyx recordings, set K, against her APL flies,
by rules locked before any response of set K is read), and whether the KCs' recruitment with and without APL at the weak
odors agrees with it (part B: the family's APL scaled by a ladder of factors at the weak stimuli, against the KC data's
response rates with APL on and with APL output blocked). A diagnostic: it moves no loss; no variant is scored on the main
target.

The entry is item E104 of experiment_priorities.md (proposed by atlas 2026-11-16; planned by a Fable subagent,
claude-fable-5-2, spawned from atlas 2026-11-16 19:44; its text in outputs/_planning/atlas_2026-11-15/entry_E104.md).
The flies: viki's copy of Yang's processed recordings (data/yang_apl_recordings/ in the shared directory; only processed
recordings, the rule of 2026-11-15 15:28). The family: C2's (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s
pulse) on E72 part B's eight wirings at their fits, nothing re-tuned; m = 1 and the acute APL-off runs are E98's B1.

The recordings (the entry's lists; part 0 checks them against recordings.csv and the stimulus files, no response table
opened):
- K6F, set K's GCaMP6f flies (MB247dbLexA > LexAop-GCaMP6f, 330;365), carrying: 2022-09-18, 10-19, 11-04, 12-01 and
  12-02, 55 processed recordings of 5 flies (9 calyces); 45 hold a paa series, 46 pa, 37 ms.
- K1020, 2022-10-20 (10 processed recordings, one fly): read apart, as a KC fly (partner (7)) and as an APL fly, until
  Yang answers question 7.
- K7F, set K's MB247 SytGCaMP7f flies (2022-01-25 to 01-27), beside: 22 recordings of 4 flies once the left calyx of
  2022-01-25 fly 1 ("left side APL labeled") is left out.
- APL, carrying: set S (E100's 29 recordings, 3 flies, 4 calyces; read here by the same reader, which reproduces E100's
  saved values) and the APL channel of the dual-colour flies on Yang's list (E101's set L, its saved c_responses.csv,
  not read again). Beside: the 12 recordings of 2021 (E98), E101's V0, set H (E100's reading) for K7F.
The reader (`recording_rows`, `read_table_k`) refuses raw-only rows, every recording outside these lists, and any
response table of set K until part A's lock (PART_A_LOCK) is committed and matches the definitions (`require_part_a`).

Definitions (the entry's; E100's for set S unless said), all in `part_a_definitions` and hashed into the lock:
- A response: Yang's calyx table (`*_roi_peak_dff.csv`) per presentation; the mean of the block's presentations minus
  the solvent-only blocks' presentation at the same place (the order-matched baseline; the mean of the two such blocks
  where a run has two, E100's `presentation_readings`). E98's floor (MIN_RESPONSE, 0.02 dF/F) and its interpolation on
  log10 concentration (`E98.value_at`) for stimuli off the whole-log grid.
- A fly is a date and a fly number; its recordings, all planes, are averaged on the log scale; a set's value is the
  geometric mean over flies with a 90 % interval from 10,000 bootstrap draws of flies (`E98.fly_ratio`).
- The falls (lower over higher concentration, within a recording): G1 paa 10^-5 / 10^-4, G2 pa 10^-5 / 10^-4, G3 ms
  10^-4 / 10^-3; beside: paa 10^-4 / 10^-3, pa 10^-6 / 10^-4, ms 10^-4 / 10^-2, each odor 10^-7 / 10^-3 (arithmetic
  means over flies).
- Q of a fall: APL's over the KCs', its interval from joint bootstrap draws of both sets' flies (`q_of_fall`).
- The split rule (E98's): two sets differ on a fall if fewer than 1 in 10 of the splits of their flies into sets of the
  same sizes give a difference of mean log falls at least as large (`split_rule`).
- The class of a fall (`fall_class`): "APL steeper" if they differ and Q < 0.8; "KCs steeper" if they differ and Q >
  1.25; "one function" if they do not differ and Q's interval lies inside 0.67 to 1.5; otherwise not decided.
- Held (`held_class`): a class counts only if the carrying reading and every partner give it. Partners: (1) set S alone
  as the APL set; (2) the block's own solvent in place of the order-matched one; (3) the block's first presentation
  alone; (4) 18 um alone; (5) without the marked runs (MARKED: notes saying aborted, "significant" warping or dim); (6)
  each KC fly left out in turn; (7) 2022-10-20 added as a KC fly.

Parts (CLI `--part`):
- 0, the identities (p0_identity.json; a failure stops the experiment): the eligibility table (p0_eligibility.csv,
  from recordings.csv and the stimulus yamls, no response table of set K opened) gives the entry's counts and every
  eligible recording of set K has Yang's table (by name); the reader reproduces E100's saved set S values
  (b_per_stimulus.csv, Yang's table, every reading) exactly; the lock (PART_A_LOCK, `write_lock`); part B's runner on
  the real wiring at m = 1 and m = 0 against E98's saved B1 on and acute-off runs (0 count differences;
  `runner_identity`, `--runs`).
- A (no simulation; only once the lock is committed): A1 strengths, A2 X1 (the classes of G1 to G3), A3 X2 (the
  family's claws, E101's a3_family.csv) and Hallem's beside, A4 X3 (the weak end), A5 beside (K7F against set H; pa and
  paa over ms; E101's Q).
- B: the ladder (m in M_LADDER, both APL scales times m) on the eight wirings at the 11 stimuli (`b_stimuli`), d(m),
  B1 (the response rate with APL on as a function of d), B2 (d*, the strength the KCs imply), Y1 and Y2 with their
  held rule, B3 described.
- C: the verdicts X1, X2, X3, Y1 and Y2 in one table with every partner reading.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part 0
    env OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E104_kc_calyx_vs_apl --part 0 --runs --threads 11      # with the runner's identities
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part 0 --write-lock   # once, before any of set K is read
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part A
    env OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E104_kc_calyx_vs_apl --part B --threads 11

Outputs: results/E104_kc_calyx_vs_apl/ in the shared directory; part B's runs under outputs/E104_kc_calyx_vs_apl/B/,
archived to data/model_outputs/E104_kc_calyx_vs_apl/ (m and the absolute APL scales of every run in params.csv).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E104_kc_calyx_vs_apl'
Stim = Tuple[str, float]
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's / E98's

# ---------------------------------------------------------------------------
# the recordings (the entry's lists; fixed before any response table of set K is opened)

K_SET = 'K'
K6F_DATES = ('20220918', '20221019', '20221104', '20221201', '20221202')      # GCaMP6f in the KCs (330;365), carrying
K1020_DATE = '20221020'                                                        # read apart (question 7)
K7F_DATES = ('20220125', '20220126', '20220127')                               # MB247 SytGCaMP7f, beside
K7F_LEFT_OUT = (('20220125', 1, 'L'),)        # (date, fly, calyx): its note says "left side APL labeled"
ODORS = ('paa', 'pa', 'ms')
# the entry's counts (part 0 checks them; a difference stops the experiment)
ENTRY_COUNTS = dict(K6F=dict(n_recordings=55, n_flies=5, n_calyces=9, planes_um=(12.0, 21.0),
                             n_series=dict(paa=45, pa=46, ms=37)),
                    K1020=dict(n_recordings=10, n_flies=1),
                    K7F=dict(n_processed=27, n_flies=4, n_recordings=22, n_flies_kept=4),
                    S=dict(n_recordings=29, n_flies=3, n_calyces=4),
                    L=dict(n_recordings=20, n_flies=7))
SERIES_MIN_CONCS = 2          # a recording holds an odor's series if it presents that odor at this many concentrations

# ---------------------------------------------------------------------------
# the definitions of part A (hashed into the lock)

MIN_RESPONSE = 0.02           # E98's floor (dF/F)
N_BOOT = 10000
SEED = 0
CI_LEVEL = 0.90
FALLS: Tuple[Tuple[str, Stim, Stim], ...] = (
    ('G1', ('paa', -5.0), ('paa', -4.0)),
    ('G2', ('pa', -5.0), ('pa', -4.0)),
    ('G3', ('ms', -4.0), ('ms', -3.0)))
FALLS_BESIDE: Tuple[Tuple[str, Stim, Stim], ...] = (
    ('paa_-4_over_-3', ('paa', -4.0), ('paa', -3.0)),
    ('pa_-6_over_-4', ('pa', -6.0), ('pa', -4.0)),
    ('ms_-4_over_-2', ('ms', -4.0), ('ms', -2.0)),
    ('paa_-7_over_-3', ('paa', -7.0), ('paa', -3.0)),
    ('pa_-7_over_-3', ('pa', -7.0), ('pa', -3.0)),
    ('ms_-7_over_-3', ('ms', -7.0), ('ms', -3.0)))
SPLIT_SHARE = 0.1
Q_APL_STEEPER = 0.8           # Q below this (and the sets differ): APL steeper
Q_KC_STEEPER = 1.25           # Q above this (and the sets differ): KCs steeper
Q_ONE_FUNCTION = (0.67, 1.5)  # the sets do not differ and Q's interval lies inside this: one function
FALL_CLASSES = ('APL steeper', 'KCs steeper', 'one function', 'not decided')
X1_MIN = 2                    # X1: at least two falls of a class, none of the opposing ones
PARTNERS = ('set_S_alone', 'block_solvent', 'first_presentation', 'plane_18', 'without_marked', 'kc_fly_left_out',
            'k1020_as_kc')
X2_MIN_WIRINGS = 6            # of the 8
X3_MIN_FLIES = 4              # of the 5 KC flies: a response at 10^-7 positive
X3_MIN_ODORS = 2              # of the three
WEAK_CONC = -7.0
PANEL_CONC = -3.0

# ---------------------------------------------------------------------------
# part B (fixed by the entry; not in the lock, which is part A's)

M_LADDER = (0.15, 0.25, 0.35, 0.5, 0.75, 1.5, 2.0)
M_SAVED = (0.0, 1.0)          # from E98's B1: the acute APL-off run and the fit
WEAK_STIMULI: Tuple[Stim, ...] = (('2h', -6.0), ('ma', -7.0), ('farn', -2.0))   # the input's (2h -7, ma -7, farn -3)
FLIES_WEAK: Dict[str, Stim] = {'2h': ('2h', -7.0), 'ma': ('ma', -7.0), 'farn': ('farn', -3.0)}
REFERENCE = 'paa'             # the panel's paa (input at 10^-4 for the flies' 10^-3)
INDICATOR_NS = tuple(float(x) for x in np.arange(1.0, 2.5 + 1e-9, 0.25))
Y_CLASSES = ('too weak', 'inside', 'too strong')
Y_MIN_ODORS = 2
PANEL_KC_RATIO = 2.28         # the KC data's APL-off / APL-on response rate ratio at the panel (the registry's)

QUALITY_CAVEAT = ("Yang has not checked these recordings' quality (her answers to questions 13 and 19, 2026-11-15 18:08: "
                  "all data can be used, but she is unsure of quality questions)")
UNMEASURED = ("unmeasured: dF/F of one calyx plane against a model variable; GCaMP6f in two cell types at different "
              "signal levels (set S's APL peaks at 3 to 5 dF/F, near saturation); other flies for each cell type; a "
              "chronic block (TNT) against the model's acute APL off; three recordings with APL output blocked at the "
              "weak odors")
PART_A_LOCK = Path(__file__).resolve().parent / 'E104_part_a_lock.json'


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the lock: part A's definitions, committed before any response table of set K is opened

def part_a_definitions() -> dict:
    """Part A's definitions, lists and thresholds as one JSON-able dict (the recordings' lists, the response, the fly,
    the falls, Q and its classes, the split rule, X1's rule, the held rule's partners, MARKED, X2's and X3's rules)."""
    raise NotImplementedError


def definitions_sha(defs: dict) -> str:
    """sha256 of the definitions' canonical JSON (sorted keys, no spaces)."""
    raise NotImplementedError


def write_lock(path: Optional[Path] = None) -> dict:
    """Writes PART_A_LOCK: the definitions, their sha256, the eligibility table's sha256, host, time and commits.
    Refuses to overwrite a lock with other definitions."""
    raise NotImplementedError


def part_a_locked(lock: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether the lock exists, its sha256 is its definitions', those equal `part_a_definitions()` and (`check_git`) the
    lock file is committed (tracked, no uncommitted change)."""
    raise NotImplementedError


def require_part_a(stem: str, **kw) -> None:
    """Raises PermissionError when `stem` is a recording of set K and `part_a_locked(**kw)` is False."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the recordings and the reader

def is_set_k(stem: str) -> bool:
    """Whether `stem` is a recording of set K (any of its dates)."""
    raise NotImplementedError


def check_row(row) -> None:
    """Raises PermissionError on a recordings.csv row that is not processed (raw only, the rule of 2026-11-15 15:28)."""
    raise NotImplementedError


def k_group(row) -> str:
    """The role of a set K row: 'K6F', 'K1020', 'K7F', 'K7F_left_out' (its calyx is left out) or 'raw_only'."""
    raise NotImplementedError


def allowed_stems() -> Dict[str, str]:
    """The recordings this module reads, stem -> role: K6F, K1020 and K7F (set K's eligible recordings), 'S' (E100's set
    S) and 'H' (E100's set H, for A5)."""
    raise NotImplementedError


def recording_rows(stems: Sequence[str]) -> pd.DataFrame:
    """recordings.csv rows of `stems`, indexed by stem, with fly_id, calyx and role; PermissionError on a stem outside
    `allowed_stems` and on a raw-only row."""
    raise NotImplementedError


def yaml_path(row) -> Path:
    """The recording's stimulus yaml in viki's copy ('<stem>.yaml', '<stem>_copy.yaml', '<stem>yaml.yaml')."""
    raise NotImplementedError


def stim_table(row) -> pd.DataFrame:
    """One recording's presentations from its yaml (E100's `yaml_stimuli`): frame, block, position, pin, odor,
    log10_conc, is_solvent, stim."""
    raise NotImplementedError


def yang_table_path(row) -> Optional[Path]:
    """Yang's calyx table of the recording (`<stem>_<tag>_roi_peak_dff.csv`), found by name only (never opened here);
    the tag of the copied peak maps where several exist (E100's `recording_paths`), else the single table; None if
    absent."""
    raise NotImplementedError


def series_concs(table: pd.DataFrame) -> Dict[str, List[float]]:
    """Each odor's presented concentrations in a stimulus table."""
    raise NotImplementedError


def eligibility_row(row, table: Optional[pd.DataFrame]) -> dict:
    """One recording's eligibility from recordings.csv and its yaml: group, fly, calyx, plane, pulse, design (nstim,
    nblock, the solvent-only blocks), each odor's concentrations and whether it holds a series, Yang's table by name,
    marked."""
    raise NotImplementedError


def eligibility_table() -> pd.DataFrame:
    """Every row of set K (raw-only rows from recordings.csv alone), plus set S's rows."""
    raise NotImplementedError


def eligibility_counts(el: pd.DataFrame) -> dict:
    """The counts the entry gives (ENTRY_COUNTS) from the eligibility table, each beside the entry's and whether equal."""
    raise NotImplementedError


def marked_runs(rows: pd.DataFrame) -> Dict[str, str]:
    """The marked runs of set K and set S: recordings.csv's notes saying aborted, "significant" warping or dim."""
    raise NotImplementedError


def read_table_k(row) -> pd.DataFrame:
    """One recording's presentations with Yang's calyx response ('yang_peak'): `require_part_a` and `check_row` first;
    the table's rows and pins must match the yaml's presentations."""
    raise NotImplementedError


def recording_readings(row) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per presentation, per stimulus) readings of one recording (E100's `presentation_readings` on Yang's table:
    raw, minus_recording, minus_block, order_matched; the first presentation of the block beside)."""
    raise NotImplementedError


def s_identity() -> Tuple[dict, pd.DataFrame]:
    """Part 0: the reader reproduces E100's saved set S values (b_per_stimulus.csv, the yang_* columns) exactly."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B's runner

def b_stimuli() -> List[str]:
    """The 11 B1 labels of a ladder run: the three weak stimuli under each session's estimate, paa and the blank."""
    raise NotImplementedError


def ladder_params(p: Dict[str, float], m: float) -> Dict[str, float]:
    """The fit's parameters with APL->KC and APL->bouton scales both times `m`."""
    raise NotImplementedError


def ladder_input() -> pd.DataFrame:
    """E98's B1 input restricted to `b_stimuli` (the panel's paa and ten added stimuli)."""
    raise NotImplementedError


def run_ladder_point(wiring: str, m: float, er=None) -> dict:
    """One wiring at `m` on `ladder_input` (E98's `ExtendedRun`): responses, spike counts, APL's means."""
    raise NotImplementedError


def runner_identity(wiring: str = 'real') -> dict:
    """Part 0: the runner at m = 1 and m = 0 against E98's saved B1 'on' and 'acute_off' runs on the 11 stimuli: 0
    differences in spike counts (and responses), APL's means beside."""
    raise NotImplementedError


def part_0(runs: bool = False) -> dict:
    """Part 0: p0_eligibility.csv, p0_counts.json, p0_s_identity.csv, p0_identity.json (and with `runs` the runner's
    identity on the real wiring); raises if an identity fails."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part A: the flies' KCs' calyx signal against their APL

def fly_falls(per_stim: pd.DataFrame, fall: Tuple[str, Stim, Stim], value: str) -> pd.Series:
    """Per fly, the mean (natural) log fall over its recordings (`E98.value_at` per recording; both above
    MIN_RESPONSE)."""
    raise NotImplementedError


def split_rule(a: pd.Series, b: pd.Series) -> dict:
    """E98's split rule on two sets of flies' mean log falls (E103's `split_rule`)."""
    raise NotImplementedError


def q_of_fall(apl: pd.Series, kc: pd.Series, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Q = APL's fall over the KCs' (geometric means over flies) with its 90 % interval from joint bootstrap draws of
    both sets' flies."""
    raise NotImplementedError


def fall_class(differ: bool, q: float, q_lo: float, q_hi: float) -> str:
    """The class of a fall (FALL_CLASSES) from the split rule's verdict and Q with its interval."""
    raise NotImplementedError


def held_class(by_reading: Dict[str, str]) -> str:
    """The carrying reading's class if every partner gives the same, else 'not decided'."""
    raise NotImplementedError


def x1_verdict(classes: Dict[str, str]) -> str:
    """X1 from the held classes of G1, G2 and G3."""
    raise NotImplementedError


def x3_verdict(per_odor_positive: Dict[str, int]) -> dict:
    """X3 from the number of KC flies whose response at 10^-7 is positive, per odor."""
    raise NotImplementedError


def part_a() -> dict:
    """Part A (A1 to A5): refuses to run until `part_a_locked()`."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part B: the KCs' somata at the weak odors

def d_of_m(apl: pd.DataFrame, m: float, stim: str) -> float:
    """APL's strength at a stimulus: m times APL's evoked activity at m (window mean minus the blank's at the same m)
    over paa's evoked activity at m = 1."""
    raise NotImplementedError


def interp_monotone(x: Sequence[float], y: Sequence[float], x0: float) -> float:
    """y at x0 by linear interpolation between the ladder's points, the points first made monotone (isotonic) in x."""
    raise NotImplementedError


def d_star(d: Sequence[float], ratio: Sequence[float], target: float) -> float:
    """The strength at which the family's response-rate ratio (APL off over on) equals `target` (interpolated)."""
    raise NotImplementedError


def flies_ratio_interval(n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """The flies' APL-off / APL-on response-rate ratio per weak odor with its interval: bootstrap of the control flies,
    each recording with APL output blocked left out in turn, the hull."""
    raise NotImplementedError


def y_class(strength: Tuple[float, float], dstar: Tuple[float, float]) -> str:
    """'too weak' if the strength's interval lies below d*'s, 'too strong' if above, otherwise 'inside'."""
    raise NotImplementedError


def y_verdict(per_odor: Dict[str, str]) -> str:
    """Y1 / Y2 over the three odors."""
    raise NotImplementedError


def part_b(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """Part B: the ladder's runs (resumable) and the tables B1 to B3."""
    raise NotImplementedError


def part_c() -> dict:
    """Part C: the verdict table."""
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C'))
    ap.add_argument('--runs', action='store_true', help="part 0: also the runner's identities on the real wiring")
    ap.add_argument('--write-lock', action='store_true', help="part 0: write part A's lock (commit it before part A)")
    ap.add_argument('--wirings', default=','.join(ALL_WIRINGS), help='part B: comma-separated wirings')
    ap.add_argument('--threads', type=int, default=None,
                    help='runs: OMP threads (OMP_NUM_THREADS and OMP_THREAD_LIMIT must match)')
    a = ap.parse_args()
    raise NotImplementedError


if __name__ == '__main__':
    main()
