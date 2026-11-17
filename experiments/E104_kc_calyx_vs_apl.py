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

Readings of the entry (fixed in part 0, before any response of set K; each also in p0_identity.json or the lock):
- "45 hold a paa series, 46 pa and 37 ms" are the K6F recordings whose stimulus file presents the odor ('n_holding').
  Ten of the 46 present pa at 10^-2 only (2022-12-01 fun5 to 7 and 2022-12-02 fun6 to 9 and FR1 to 3: Yang's "PA -2!!!"
  pin), so 36 hold a pa series of two or more concentrations ('n_series', beside); they enter no fall.
- The plane: recordings.csv's (Yang's yaml / date note); 20220918_1FRfun2 has none there and takes its Experiment.xml's
  piezo position (18 um, as its date note's '1FRfun1_001'). Three recordings' yaml and Experiment.xml disagree by 3 um
  (PLANE_UNCERTAIN); partner (4), 18 um alone, leaves them out.
- Yang's table where a recording has several (one per registration): TABLE_TAG_ORDER (rigid motion-corrected first, as
  the copy chose the peak maps; for set S this is the table E100 read). 20221019_1FLfun8 has the plain and the
  piecewise-rigid ('els_moco') tables only; the motion-corrected one is read.
- The marked runs (MARKED) are those MARK_RULES find in recordings.csv's notes; 20221104_1FLfun3's note says the run
  was aborted in block 10 and re-run under the same name, and it is marked.
- Set L (E101's c_responses.csv, not read again) is E101's carrying reading (the map-based calyx mean minus the
  first-in-block solvent presentations); partner (2) is its 'block' reading and partner (3) its 'first_two', the
  nearest E101 saved.
- X2's model ratio set: per wiring, the three claw versions of E101's a3_family.csv through E98's indicator family (n 1
  to 2.5, both forms), "inside" if any lies in set K's interval (E98's `classify`).

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
K_DATES = K7F_DATES + K6F_DATES[:2] + (K1020_DATE,) + K6F_DATES[2:]
ODORS = ('paa', 'pa', 'ms')
# the lists, fixed 2026-11-16 from recordings.csv (processed rows) and checked by part 0 (`eligibility_identity`)
K6F: Tuple[str, ...] = (
    '20220918_1FLfun1', '20220918_1FLfun2', '20220918_1FLfun3', '20220918_1FLfun4', '20220918_1FLfun5', '20220918_1FRfun1',
    '20220918_1FRfun2', '20220918_1FRfun3', '20221019_1FLfun1', '20221019_1FLfun2', '20221019_1FLfun4', '20221019_1FLfun5',
    '20221019_1FLfun6', '20221019_1FLfun7', '20221019_1FLfun8', '20221019_1FLfun9', '20221019_1FRfun1', '20221019_1FRfun2',
    '20221019_1FRfun3', '20221019_1FRfun4', '20221104_1FLfun1', '20221104_1FLfun2', '20221104_1FLfun3', '20221104_1FLfun4',
    '20221104_1FLfun5', '20221104_1FLfun6', '20221104_1FLfun7', '20221104_1FLfun8', '20221104_1FRfun1', '20221104_1FRfun2',
    '20221104_1FRfun3', '20221104_1FRfun4', '20221201_1FLfun1', '20221201_1FLfun2', '20221201_1FLfun3', '20221201_1FLfun4',
    '20221201_1FLfun5', '20221201_1FLfun6', '20221201_1FLfun7', '20221201_1FLfun8', '20221201_1FLfun9', '20221202_1FLfun1',
    '20221202_1FLfun2', '20221202_1FLfun3', '20221202_1FLfun4', '20221202_1FLfun5', '20221202_1FLfun6', '20221202_1FLfun7',
    '20221202_1FLfun8', '20221202_1FLfun9', '20221202_1FRfun1', '20221202_1FRfun2', '20221202_1FRfun3', '20221202_1FRfun4',
    '20221202_1FRfun5')
K1020: Tuple[str, ...] = (
    '20221020_1FLfun1', '20221020_1FLfun2', '20221020_1FLfun3', '20221020_1FLfun4', '20221020_1FLfun5', '20221020_1FLfun6',
    '20221020_1FLfun7', '20221020_1FRfun2', '20221020_1FRfun3', '20221020_1FRfun4')
K7F: Tuple[str, ...] = (
    '20220125_1FRfun1', '20220125_1FRfun3', '20220125_1FRfun4', '20220126_1FLfun1', '20220126_1FLfun2', '20220126_1FLfun3',
    '20220126_1FRfun1', '20220126_1FRfun2', '20220126_1FRfun3', '20220127_1MLfun1', '20220127_1MLfun2', '20220127_1MLfun3',
    '20220127_1MLfun4', '20220127_1MRfun1', '20220127_1MRfun2', '20220127_1MRfun3', '20220127_2FLfun1', '20220127_2FLfun2',
    '20220127_2FLfun3', '20220127_2FRfun1', '20220127_2FRfun2', '20220127_2FRfun3')
K7F_LEFT_OUT_STEMS: Tuple[str, ...] = (
    '20220125_1FLfun1', '20220125_1FLfun2', '20220125_1FLfun3', '20220125_1FLfun4', '20220125_1FLfun5')
S_STEMS: Tuple[str, ...] = (     # E100's set S (b_recordings.csv), read here by the same reader
    '20221031_1FLfun1', '20221031_1FLfun2', '20221031_1FLfun3', '20221031_1FLfun4', '20221031_1FLfun5', '20221031_1FLfun6',
    '20221103_1FRfun2', '20221103_1FRfun3', '20221103_1FRfun4', '20221103_1FRfun5', '20221103_1FRfun6', '20221103_1FRfun7',
    '20221103_1FRfun8', '20221103_2FLfun1', '20221103_2FLfun2', '20221103_2FLfun3', '20221103_2FLfun4', '20221103_2FLfun5',
    '20221103_2FLfun6', '20221103_2FLfun7', '20221103_2FRfun1', '20221103_2FRfun2', '20221103_2FRfun3', '20221103_2FRfun4',
    '20221103_2FRfun5', '20221103_2FRfun6', '20221103_2FRfun7', '20221103_2FRfun8', '20221103_2FRfun9')
# the marked runs (the entry: notes saying aborted, "significant" warping or dim), by MARK_RULES on recordings.csv's
# quality_and_notes of the lists above; part 0 checks that the rules give exactly these
MARK_RULES: Tuple[Tuple[str, str], ...] = (('aborted', r'abort'), ('significant warping', r'sign\w*cant\s+warp'),
                                            ('dim', r'\bdim\b'))
MARKED: Dict[str, str] = {
    '20220918_1FLfun1': 'significant warping', '20221019_1FLfun1': 'significant warping',
    '20221019_1FLfun2': 'significant warping', '20221104_1FLfun1': 'aborted', '20221104_1FLfun3': 'aborted',
    '20221201_1FLfun1': 'significant warping', '20221202_1FLfun1': 'significant warping', '20221202_1FLfun3': 'dim'}
# Yang's calyx table of a recording with several (one per registration of its peak maps): the copy's rule for the peak
# maps (the rigid motion-corrected one where it exists), then the other motion-corrected ones, then the plain one
TABLE_TAG_ORDER = ('rig_moco', 'moco', 'els_moco', '')
# the plane: recordings.csv's (Yang's yaml z_position / date note); the Experiment.xml's piezo position where it has none
# (20220918_1FRfun2: 18 um, as the date note's '1FRfun1_001, z at 18 um'). These three have a yaml and an
# Experiment.xml that disagree by 3 um (yaml 18 / 15 / 18, xml 21 / 18 / 21); partner (4) leaves them out
PLANE_UNCERTAIN = ('20221019_1FRfun3', '20221104_1FLfun2', '20221104_1FLfun6')
# the entry's counts (part 0 checks them; a difference stops the experiment). 'n_holding': recordings whose stimulus
# file presents the odor (the entry's "hold a ... series"); 'n_series', at two or more concentrations, is given beside
ENTRY_COUNTS = dict(K6F=dict(n_recordings=55, n_flies=5, n_calyces=9, plane_min_um=12.0, plane_max_um=21.0,
                             n_holding=dict(paa=45, pa=46, ms=37), n_with_table=55),
                    K1020=dict(n_recordings=10, n_flies=1, n_with_table=10),
                    K7F=dict(n_processed=27, n_flies_processed=4, n_recordings=22, n_flies=4, n_with_table=22),
                    S=dict(n_recordings=29, n_flies=3, n_calyces=4),
                    L=dict(n_recordings=20, n_flies=7))
SERIES_MIN_CONCS = 2          # 'n_series': a recording presents the odor at this many concentrations or more

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
WEAK_FALL_IDS = ('paa_-7_over_-3', 'pa_-7_over_-3', 'ms_-7_over_-3')   # arithmetic means over flies
SPLIT_SHARE = 0.1
Q_APL_STEEPER = 0.8           # Q below this (and the sets differ): APL steeper
Q_KC_STEEPER = 1.25           # Q above this (and the sets differ): KCs steeper
Q_ONE_FUNCTION = (0.67, 1.5)  # the sets do not differ and Q's interval lies inside this: one function
FALL_CLASSES = ('APL steeper', 'KCs steeper', 'one function', 'not decided')
X1_MIN = 2                    # X1: at least two falls of a class, none of the opposing ones
PARTNERS = ('set_S_alone', 'block_solvent', 'first_presentation', 'plane_18', 'without_marked', 'kc_fly_left_out',
            'k1020_as_kc')
CARRYING_PLANE_PARTNER = 18.0
X2_MIN_WIRINGS = 6            # of the 8
X2_CLAW_COLUMNS = ('claw_all_00', 'claw_all_25', 'claw_all_50')   # E101's a3_family.csv, all claws, its three versions
X3_MIN_FLIES = 4              # of the 5 KC flies: a response at 10^-7 positive
X3_MIN_ODORS = 2              # of the three
WEAK_CONC = -7.0
PANEL_CONC = -3.0
A5_ODORS = ('2h', 'ma', 'pa', 'paa')
A5_ACROSS: Tuple[Tuple[str, Stim, Stim], ...] = (('pa_over_ms', ('pa', -3.0), ('ms', -3.0)),
                                                 ('paa_over_ms', ('paa', -3.0), ('ms', -3.0)))

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
RUNNER_CONDS = ((1.0, 'on'), (0.0, 'acute_off'))   # part 0: the runner at m against E98's saved B1 run of that tag

QUALITY_CAVEAT = ("Yang has not checked these recordings' quality (her answers to questions 13 and 19, 2026-11-15 18:08: "
                  "all data can be used, but she is unsure of quality questions)")
UNMEASURED = ("unmeasured: dF/F of one calyx plane against a model variable; GCaMP6f in two cell types at different "
              "signal levels (set S's APL peaks at 3 to 5 dF/F, near saturation); other flies for each cell type; a "
              "chronic block (TNT) against the model's acute APL off; three recordings with APL output blocked at the "
              "weak odors")
PART_A_LOCK = Path(__file__).resolve().parent / 'E104_part_a_lock.json'
IDENTITY_TOL = 0.0            # the reader against E100's saved set S values: exact


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _json_default(x):
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    Path(path).write_text(json.dumps(obj, indent=1, default=_json_default) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def write_table(df: pd.DataFrame, path: Path) -> None:
    """A results table: its first line Yang's caveat ('# ' QUALITY_CAVEAT), its second what stays unmeasured, then the
    csv (`read_table` reads it)."""
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


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, experiment=EXPERIMENT, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits={k: v for k, v in H.code_commits().items()
                                                  if not isinstance(v, (dict, list))},
                                         host=_host(), time=time.strftime('%Y-%m-%d %H:%M'), **extra))
    _write_json(path, d)


def _e98():
    """E98's module, E79's data paths repointed to this host's shared directory first (E100's `_ds`)."""
    from experiments import E100_apl_weak_odors as E100
    E100._ds()
    from experiments import E98_apl_response_strength as E98
    return E98


def _e100():
    from experiments import E100_apl_weak_odors as E100
    return E100


# ---------------------------------------------------------------------------
# the lock: part A's definitions, committed before any response table of set K is opened

def part_a_definitions() -> dict:
    """Part A's definitions, lists and thresholds as one JSON-able dict (constants of this module only, so that the lock
    can be checked without the shared directory): the recordings' lists and their roles, the marked runs, Yang's table
    and the response, the fly and the set's value, the falls, Q and its classes, the split rule, X1's rule, the held
    rule's partners, X2's, X3's and A5's rules."""
    return dict(
        experiment=EXPERIMENT,
        recordings=dict(
            K6F=list(K6F), K1020=list(K1020), K7F=list(K7F), K7F_left_out=list(K7F_LEFT_OUT_STEMS), S=list(S_STEMS),
            roles=dict(K6F="the KCs, carrying (set K's GCaMP6f flies)", K1020="2022-10-20: partner (7) as a KC fly, "
                       "and described as an APL fly; in no carrying set", K7F="beside (A5), against set H",
                       S="APL, carrying, with L", L="APL, carrying, with S: E101's c_responses.csv, group 'L', not read "
                       "again", H="beside (A5): E100's b_per_stimulus.csv, set H, label_unsure False, not read again"),
            L=dict(source='results/E101_apl_vs_kc_claws/c_responses.csv', group='L', b_mult=0.0, column='apl',
                   carrying='carrying', block_solvent='block', first_presentation='first_two (E101 saved no single '
                   'first presentation)'),
            H=dict(source='results/E100_apl_weak_odors/b_per_stimulus.csv', set='H', label_unsure=False,
                   carrying='yang_order_matched', block_solvent='yang_minus_block'),
            beside=dict(twelve_2021="E98's a1_responses.csv", V0="E101's c_responses.csv, group 'V0'")),
        marked=dict(rules=[list(r) for r in MARK_RULES], runs=dict(sorted(MARKED.items()))),
        response=dict(
            table="Yang's `<stem>_<tag>_roi_peak_dff.csv`, column roi_peak_dff, one row per presentation in the yaml's "
                  "odor_stim_seq order (rows and pins must match the yaml's, else the recording has no value)",
            table_tag_order=list(TABLE_TAG_ORDER),
            carrying="order_matched: each presentation minus the mean of the solvent-only blocks' presentations at "
                     "the same position in the block (E100's `presentation_readings`); a stimulus's response the mean "
                     "over its presentations (both passes where a run has two); no value in a run without a "
                     "solvent-only block",
            block_solvent="minus_block: each presentation minus its own block's solvent presentation(s)",
            first_presentation="the order-matched value of the stimulus's first presentation (its first block)",
            min_response=MIN_RESPONSE,
            interpolation="E98's `value_at` on log10 concentration inside a recording's imaged range; outside it the "
                          "nearest end within one log unit, inclusive"),
        fly=dict(id="date and fly number (E98's `fly_id`)", recordings="all planes; a fly's falls are the mean of its "
                 "recordings' natural-log falls (both responses above min_response)", value="the geometric mean over "
                 "flies (the arithmetic mean for the 10^-7 / 10^-3 falls, WEAK_FALL_IDS, and for X3's responses)",
                 n_boot=N_BOOT, seed=SEED, ci_level=CI_LEVEL, bootstrap="E98's `fly_ratio` (E79's `fly_bootstrap`)"),
        falls=[[f, list(a), list(b)] for f, a, b in FALLS],
        falls_beside=[[f, list(a), list(b)] for f, a, b in FALLS_BESIDE], weak_fall_ids=list(WEAK_FALL_IDS),
        q=dict(definition="APL's fall over the KCs' (geometric means over flies)", interval="joint bootstrap: "
               "n_boot draws, each resampling both sets' flies with replacement; the 5th and 95th percentiles"),
        split_rule=dict(share=SPLIT_SHARE, definition="two sets differ on a fall if fewer than share of the splits of "
                        "all their flies into sets of the same sizes give a difference of mean log falls at least as "
                        "large (absolute; E98's `pulse_permutation`)"),
        classes=dict(apl_steeper=Q_APL_STEEPER, kcs_steeper=Q_KC_STEEPER, one_function=list(Q_ONE_FUNCTION),
                     labels=list(FALL_CLASSES)),
        x1=dict(min_falls=X1_MIN, rule="'APL steeper' if at least min_falls of G1 to G3 are held APL steeper and none "
                "KCs steeper; 'one function' if at least min_falls are held one function and none steeper either way; "
                "otherwise not decided"),
        held=dict(partners=list(PARTNERS), plane_um=CARRYING_PLANE_PARTNER,
                  plane="recordings.csv's plane_depth_um (Yang's yaml / date note), the Experiment.xml's piezo "
                        "position where it has none; partner (4) leaves out PLANE_UNCERTAIN (yaml and xml disagree)",
                  plane_uncertain=list(PLANE_UNCERTAIN),
                  rule="a class counts only if the carrying reading and every partner give it; partner (6) is one "
                       "reading per KC fly left out, each of which must give it"),
        x2=dict(table='results/E101_apl_vs_kc_claws/a3_family.csv', cond='on', columns=list(X2_CLAW_COLUMNS),
                level="the window mean minus the blank's (same wiring and condition)", fall='G1',
                numerator="paa at 10^-5 under each of the three sessions' estimates", denominator="the panel's paa",
                indicator_ns=list(INDICATOR_NS), forms=['baseline', 'evoked'], min_wirings=X2_MIN_WIRINGS,
                rule="per wiring E98's `classify` of the model's ratio set (the three claw versions through E98's "
                     "indicator family) against set K's G1 interval (inside, flatter, steeper); a class takes "
                     "min_wirings of 8; over the estimates E100's rule (inside if inside under any; flatter or "
                     "steeper only if so under all)"),
        x3=dict(conc=WEAK_CONC, reference_conc=PANEL_CONC, min_flies=X3_MIN_FLIES, min_odors=X3_MIN_ODORS,
                rule="an odor's KC response at 10^-7 (carrying reading) is above the solvent's if it is positive in "
                     "at least min_flies of the 5 KC flies; X3 holds if so for at least min_odors of paa, pa and ms"),
        a5=dict(odors=list(A5_ODORS), across=[[f, list(a), list(b)] for f, a, b in A5_ACROSS],
                rule="K7F against set H by the same rules (described; changes no class)"))


def _canonical(defs: dict) -> str:
    return json.dumps(defs, sort_keys=True, separators=(',', ':'), default=_json_default)


def definitions_sha(defs: dict) -> str:
    """sha256 of the definitions' canonical JSON (sorted keys, no spaces)."""
    return hashlib.sha256(_canonical(defs).encode()).hexdigest()


def write_lock(path: Optional[Path] = None, eligibility: Optional[Path] = None) -> dict:
    """Writes the lock (default PART_A_LOCK): the definitions, their sha256, the sha256 of part 0's eligibility table
    and identity file (recorded, not checked), host, time and commits. Refuses to overwrite a lock with other
    definitions."""
    path = PART_A_LOCK if path is None else Path(path)
    defs = part_a_definitions()
    sha = definitions_sha(defs)
    old = _read_json(path)
    if old is not None and old.get('sha256') != sha:
        raise RuntimeError(f'{path.name} exists with other definitions ({old.get("sha256")}); not overwritten')
    files = {}
    for p in ([Path(eligibility)] if eligibility is not None else [out_dir() / 'p0_eligibility.csv',
                                                                    out_dir() / 'p0_identity.json']):
        if p.exists():
            files[p.name] = _sha256(p)
    lock = dict(experiment=EXPERIMENT, part='A', sha256=sha, definitions=defs, part_0_files=files,
                host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                commits={k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))},
                note="Part A's definitions, written by part 0 before any response table of set K was opened. "
                     "`require_part_a` refuses to open set K's tables until this file is committed, its sha256 is its "
                     "definitions' and they equal the module's `part_a_definitions()`.")
    _write_json(path, lock)
    return lock


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


def part_a_locked(lock: Optional[Path] = None, check_git: bool = True) -> bool:
    """Whether the lock exists, its sha256 is its definitions', those equal `part_a_definitions()` and (`check_git`) the
    lock file is committed (tracked, no uncommitted change)."""
    lock = PART_A_LOCK if lock is None else Path(lock)
    if not lock.exists():
        return False
    try:
        d = json.loads(lock.read_text())
    except ValueError:
        return False
    defs = d.get('definitions')
    if not isinstance(defs, dict) or d.get('sha256') != definitions_sha(defs):
        return False
    if _canonical(defs) != _canonical(part_a_definitions()):
        return False
    return _committed(lock) if check_git else True


def require_part_a(stem: str, **kw) -> None:
    """Raises PermissionError when `stem` is a recording of set K and `part_a_locked(**kw)` is False: no response table
    of set K is opened before part A's lock is committed."""
    if is_set_k(stem) and not part_a_locked(**kw):
        raise PermissionError(f"{stem}: no response table of set K is opened before part A's lock "
                              f"({PART_A_LOCK.name}) is committed and matches the definitions")


# ---------------------------------------------------------------------------
# the recordings and the reader

def is_set_k(stem: str) -> bool:
    """Whether `stem` is a recording of set K (any of its dates)."""
    return str(stem)[:8] in K_DATES


def check_row(row) -> None:
    """Raises PermissionError on a recordings.csv row that is not processed (raw only, the rule of 2026-11-15 15:28)."""
    v = row.get('use_until_further_notice', '') if hasattr(row, 'get') else ''
    if not str(v).startswith('processed'):
        raise PermissionError(f"{row.get('stem')}: not processed ({v}); raw-only recordings are not read until further "
                              "notice (rule of 2026-11-15 15:28)")


def k_group(row) -> str:
    """The role of a set K row of recordings.csv: 'raw_only' (not processed), 'K6F', 'K1020', 'K7F', 'K7F_left_out'
    (its calyx is left out, K7F_LEFT_OUT) or 'other' (a date outside the entry's)."""
    if not str(row.get('use_until_further_notice', '')).startswith('processed'):
        return 'raw_only'
    d = str(row['date']).replace('-', '')
    if d in K6F_DATES:
        return 'K6F'
    if d == K1020_DATE:
        return 'K1020'
    if d in K7F_DATES:
        key = (d, int(float(row['fly'])), str(row['hemisphere']))
        return 'K7F_left_out' if key in set(K7F_LEFT_OUT) else 'K7F'
    return 'other'


def allowed_stems() -> Dict[str, str]:
    """The recordings this module's reader opens, stem -> role: K6F, K1020 and K7F (set K's eligible recordings) and
    'S' (E100's set S). Set H, set L and the 2021 recordings are read from E100's, E101's and E98's saved tables."""
    out = {st: 'K6F' for st in K6F}
    out.update({st: 'K1020' for st in K1020})
    out.update({st: 'K7F' for st in K7F})
    out.update({st: 'S' for st in S_STEMS})
    return out


def read_recordings() -> pd.DataFrame:
    """viki's recordings.csv, every set (E100's `read_recordings`: with 'folder' and 'usable')."""
    return _e100().read_recordings(sets=None)


def recording_rows(stems: Sequence[str], recs: Optional[pd.DataFrame] = None,
                   allowed: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """recordings.csv rows of `stems` (`recs`, default `read_recordings()`), indexed by stem, with fly_id, calyx and
    role; PermissionError on a stem outside `allowed` (default `allowed_stems()`) and on a raw-only row."""
    allowed = allowed_stems() if allowed is None else allowed
    bad = [x for x in stems if x not in allowed]
    if bad:
        raise PermissionError(f"E104 reads only set K's eligible recordings and set S: {bad}")
    recs = read_recordings() if recs is None else recs
    out = []
    for st in stems:
        g = recs[(recs['stem'] == st) & (recs['date'].astype(str).str.replace('-', '') == st[:8])]
        if len(g) != 1:
            raise RuntimeError(f'{st}: {len(g)} rows in recordings.csv')
        row = g.iloc[0].to_dict()
        check_row(row)
        row['fly_id'] = f"{str(row['date']).replace('-', '')[:4]}-{str(row['date']).replace('-', '')[4:6]}-" \
                        f"{str(row['date']).replace('-', '')[6:8]}_fly{int(float(row['fly']))}"
        row['calyx'] = f"{row['fly_id']}_{row['hemisphere']}"
        row['role'] = allowed[st]
        out.append(row)
    return pd.DataFrame(out).set_index('stem', drop=False)


def yaml_path(row) -> Optional[Path]:
    """The recording's stimulus yaml in viki's copy ('<stem>.yaml', '<stem>_copy.yaml', '<stem>yaml.yaml'), or None."""
    f, st = Path(row['folder']), row['stem']
    return next((x for x in (f / f'{st}.yaml', f / f'{st}_copy.yaml', f / f'{st}yaml.yaml') if x.exists()), None)


def stim_table(row) -> pd.DataFrame:
    """One recording's presentations from its yaml (E100's `yaml_stimuli`: frame, block, position, pin, odor_name,
    odor, log10_conc), with is_solvent and stim as E100's `read_recording` sets them."""
    E100 = _e100()
    p = yaml_path(row)
    if p is None:
        raise FileNotFoundError(f"{row['stem']}: no stimulus yaml")
    t = E100.yaml_stimuli(p, row['stem'])
    t['is_solvent'] = t['odor'] == E100.SOLVENT
    t['stim'] = [E100.SOLVENT if sv else E100.stim_label((o, c)) for sv, o, c in
                 zip(t['is_solvent'], t['odor'], t['log10_conc'])]
    return t


def yang_tables(row) -> Dict[str, Path]:
    """Yang's calyx tables of the recording by name only (never opened here): tag ('' for '<stem>_roi_peak_dff.csv')
    -> path."""
    f, st = Path(row['folder']), row['stem']
    out = {}
    for p in sorted(f.glob(f'{st}_*roi_peak_dff.csv')) + sorted(f.glob(f'{st}_roi_peak_dff.csv')):
        tag = p.name[len(st) + 1:-len('roi_peak_dff.csv')].rstrip('_')
        out[tag] = p
    return out


def yang_table_path(row) -> Optional[Path]:
    """The table the reader uses (TABLE_TAG_ORDER over `yang_tables`), or None."""
    tabs = yang_tables(row)
    return next((tabs[t] for t in TABLE_TAG_ORDER if t in tabs), None)


def series_concs(table: pd.DataFrame) -> Dict[str, List[float]]:
    """Each odor's presented concentrations in a stimulus table (solvent left out)."""
    o = table[~table['is_solvent'].astype(bool)]
    return {od: sorted(float(round(c, 3)) for c in set(g['log10_conc'])) for od, g in o.groupby('odor')}


def xml_plane_um(row) -> float:
    """The plane from the recording's Experiment.xml: ThorImage's piezo ('ThorZPiezo') setupPositionMM in um; NaN if
    absent."""
    import xml.etree.ElementTree as ET
    f, st = Path(row['folder']), row['stem']
    x = next(iter(sorted(f.glob(f'{st}*Experiment.xml'))), None)
    if x is None:
        return np.nan
    for z in ET.parse(x).getroot().iter('ZStage'):
        if z.get('name') == 'ThorZPiezo' and z.get('setupPositionMM') is not None:
            return round(float(z.get('setupPositionMM')) * 1000.0, 3)
    return np.nan


def mark_reason(notes) -> str:
    """The first of MARK_RULES whose pattern is in `notes` (case-insensitive), or ''."""
    if not isinstance(notes, str):
        return ''
    return next((name for name, pat in MARK_RULES if re.search(pat, notes, re.I)), '')


def eligibility_row(row, table: Optional[pd.DataFrame], group: str) -> dict:
    """One recording's eligibility from recordings.csv and its yaml (`table`; None for a raw-only row, whose files are
    not read): group, fly, calyx, plane (recordings.csv's; the Experiment.xml's piezo position where recordings.csv has
    none, and beside), pulse, design (nstim, nblock, the solvent-only blocks, two passes), each odor's concentrations,
    whether it presents the odor ('holds_<odor>') and a series of it ('series_<odor>'), Yang's tables by name, the one
    read, and the marked reason."""
    d = str(row['date']).replace('-', '')
    fly_id = f'{d[:4]}-{d[4:6]}-{d[6:8]}_fly{int(float(row["fly"]))}'
    rec = dict(stem=row['stem'], set=row['set'], group=group, date=d, fly=fly_id,
               calyx=f"{fly_id}_{row['hemisphere']}" if isinstance(row['hemisphere'], str) else '',
               use=str(row.get('use_until_further_notice', '')), notes=row.get('quality_and_notes'))
    if table is None:
        return rec
    p_csv = float(row['plane_depth_um']) if pd.notna(row.get('plane_depth_um')) else np.nan
    p_xml = xml_plane_um(row)
    concs = series_concs(table)
    blocks = table.groupby('block')['is_solvent']
    solv_only = sorted(int(b) for b, g in blocks if g.all())
    per_stim = table[~table['is_solvent']].groupby('stim')['block'].nunique()
    tabs = yang_tables(row)
    used = yang_table_path(row)
    rec.update(plane_um=p_csv if np.isfinite(p_csv) else p_xml,
               plane_source='recordings.csv' if np.isfinite(p_csv) else ('Experiment.xml' if np.isfinite(p_xml) else ''),
               plane_xml_um=p_xml, pulse_s=float(row['pulse_s']), nstim=int(table['position'].max()),
               nblock=int(table['block'].max()) + 1, n_presentations=int(len(table)),
               n_presentations_recordings_csv=float(row['n_presentations']) if pd.notna(row.get('n_presentations'))
               else np.nan, design=_e100().block_design(table), solvent_only_blocks=';'.join(map(str, solv_only)),
               n_solvent_only_blocks=len(solv_only), two_passes=bool(len(per_stim) and (per_stim > 1).any()),
               odors=';'.join(f"{o} {' '.join(f'{c:g}' for c in v)}" for o, v in sorted(concs.items())),
               yang_tables=';'.join(p.name for p in tabs.values()), yang_table=used.name if used is not None else '',
               has_yang_table=used is not None, marked=mark_reason(row.get('quality_and_notes')))
    for o in ODORS + ('2h', 'ma'):
        rec[f'n_concs_{o}'] = len(concs.get(o, []))
        rec[f'holds_{o}'] = bool(concs.get(o))
        rec[f'series_{o}'] = bool(len(concs.get(o, [])) >= SERIES_MIN_CONCS)
    return rec


def eligibility_table(recs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Every row of set K (raw-only rows and other dates from recordings.csv alone, no file of theirs read) and set S's
    29, from recordings.csv and the stimulus yamls; no response table is opened."""
    recs = read_recordings() if recs is None else recs
    out = []
    k = recs[recs['set'] == K_SET]
    for _, r in k.iterrows():
        row = r.to_dict()
        g = k_group(row)
        out.append(eligibility_row(row, stim_table(row) if g != 'raw_only' else None, g))
    s = recs[recs['stem'].isin(S_STEMS) & (recs['set'] == 'S')]
    for _, r in s.iterrows():
        row = r.to_dict()
        out.append(eligibility_row(row, stim_table(row) if str(row.get('use_until_further_notice', '')).startswith(
            'processed') else None, 'S'))
    return pd.DataFrame(out)


def eligibility_counts(el: pd.DataFrame, l_counts: Optional[dict] = None) -> dict:
    """The counts the entry gives (ENTRY_COUNTS) from the eligibility table (and set L's from E101's table, `l_counts`),
    each beside the entry's; 'equal' per group and overall ('passed'). Also the counts the entry does not give but the
    reading depends on: each odor's series (two or more concentrations) and the recordings without a solvent-only
    block."""
    def grp(name):
        return el[el['group'] == name]

    def n_flies(x):
        return int(x['fly'].nunique())

    k6 = grp('K6F')
    got = dict(K6F=dict(n_recordings=int(len(k6)), n_flies=n_flies(k6), n_calyces=int(k6['calyx'].nunique()),
                        plane_min_um=float(k6['plane_um'].min()) if len(k6) else np.nan,
                        plane_max_um=float(k6['plane_um'].max()) if len(k6) else np.nan,
                        n_holding={o: int(k6[f'holds_{o}'].astype(bool).sum()) for o in ODORS} if len(k6) else {},
                        n_with_table=int(k6['has_yang_table'].astype(bool).sum()) if len(k6) else 0),
               K1020=dict(n_recordings=int(len(grp('K1020'))), n_flies=n_flies(grp('K1020')),
                          n_with_table=int(grp('K1020')['has_yang_table'].astype(bool).sum()) if len(grp('K1020')) else 0),
               K7F=dict(n_processed=int(len(grp('K7F')) + len(grp('K7F_left_out'))),
                        n_flies_processed=n_flies(pd.concat([grp('K7F'), grp('K7F_left_out')])),
                        n_recordings=int(len(grp('K7F'))), n_flies=n_flies(grp('K7F')),
                        n_with_table=int(grp('K7F')['has_yang_table'].astype(bool).sum()) if len(grp('K7F')) else 0),
               S=dict(n_recordings=int(len(grp('S'))), n_flies=n_flies(grp('S')), n_calyces=int(grp('S')['calyx'].nunique())))
    if l_counts is not None:
        got['L'] = dict(l_counts)
    res = {}
    for g, want in ENTRY_COUNTS.items():
        if g not in got:
            res[g] = dict(entry=want, got=None, equal=False)
            continue
        res[g] = dict(entry=want, got=got[g], equal=bool(all(got[g].get(k) == v for k, v in want.items())))
    beside = dict(K6F_n_series={o: int(k6[f'series_{o}'].astype(bool).sum()) for o in ODORS} if len(k6) else {},
                  K6F_holding_one_concentration_only={o: sorted(k6.loc[k6[f'holds_{o}'].astype(bool)
                                                                       & ~k6[f'series_{o}'].astype(bool), 'stem'])
                                                      for o in ODORS} if len(k6) else {},
                  K6F_without_solvent_only_block=sorted(k6.loc[k6['n_solvent_only_blocks'] == 0, 'stem']) if len(k6) else [],
                  K6F_planes_from_xml=sorted(k6.loc[k6['plane_source'] == 'Experiment.xml', 'stem']) if len(k6) else [],
                  n_raw_only=int((el['group'] == 'raw_only').sum()), n_other=int((el['group'] == 'other').sum()))
    return dict(passed=bool(all(r['equal'] for r in res.values())), groups=res, beside=beside)


def eligibility_identity(el: pd.DataFrame) -> dict:
    """The module's lists against the eligibility table: K6F, K1020, K7F, K7F_LEFT_OUT_STEMS and S_STEMS equal its
    groups; every eligible recording has Yang's table (by name) and a yaml whose presentations equal recordings.csv's;
    the Experiment.xml's plane differs from recordings.csv's exactly at PLANE_UNCERTAIN; MARK_RULES give MARKED on the
    lists' notes."""
    lists = dict(K6F=K6F, K1020=K1020, K7F=K7F, K7F_left_out=K7F_LEFT_OUT_STEMS, S=S_STEMS)
    res = {g: dict(equal=bool(sorted(el.loc[el['group'] == g, 'stem']) == sorted(v)),
                   missing=sorted(set(v) - set(el.loc[el['group'] == g, 'stem'])),
                   extra=sorted(set(el.loc[el['group'] == g, 'stem']) - set(v))) for g, v in lists.items()}
    elig = el[el['group'].isin(['K6F', 'K1020', 'K7F', 'S'])]
    no_tab = sorted(elig.loc[~elig['has_yang_table'].astype(bool), 'stem'])
    pres = elig[elig['n_presentations_recordings_csv'].notna()]
    pres_bad = sorted(pres.loc[pres['n_presentations'] != pres['n_presentations_recordings_csv'], 'stem'])
    both = elig[elig['plane_source'].eq('recordings.csv') & elig['plane_xml_um'].notna()]
    plane_bad = sorted(both.loc[(both['plane_um'] - both['plane_xml_um']).abs() > 1e-6, 'stem'])
    marked = {st: r for st, r in zip(elig['stem'], elig['marked']) if r}
    out = dict(lists=res, without_table=no_tab, presentations_differ=pres_bad,
               n_planes_checked=int(len(both)), planes_differ=plane_bad,
               planes_differ_as_listed=bool(plane_bad == sorted(PLANE_UNCERTAIN)), marked=marked,
               marked_equal=bool(marked == MARKED))
    out['passed'] = bool(all(r['equal'] for r in res.values()) and not no_tab and not pres_bad
                         and out['planes_differ_as_listed'] and out['marked_equal'])
    return out


def set_l_counts() -> dict:
    """Set L's recordings and flies in E101's c_responses.csv (group 'L', its carrying reading at b = 0)."""
    c = pd.read_csv(H.shared_dir() / 'results' / 'E101_apl_vs_kc_claws' / 'c_responses.csv', comment='#')
    c = c[(c['group'] == 'L') & (c['reading'] == 'carrying') & (c['b_mult'] == 0.0)]
    return dict(n_recordings=int(c['recording'].nunique()), n_flies=int(c['fly'].nunique()))


def read_table_k(row, table: Optional[pd.DataFrame] = None) -> Tuple[pd.DataFrame, str]:
    """One recording's presentations (`stim_table`, or `table`) with Yang's calyx response ('yang_peak') and the
    table's status ('ok', or why it has no value: rows or pins differ from the yaml's). `require_part_a`, the allowed
    lists and `check_row` first: no table of set K before the lock, nothing outside the lists, no raw-only row."""
    st = row['stem']
    require_part_a(st)
    if st not in allowed_stems():
        raise PermissionError(f"{st}: outside E104's lists")
    check_row(row)
    fr = stim_table(row) if table is None else table.copy()
    p = yang_table_path(row)
    if p is None:
        return fr, 'missing'
    t = pd.read_csv(p)
    if len(t) != len(fr):
        return fr, f'{len(t)} rows for {len(fr)} presentations'
    if 'pin_n' in t and not (t['pin_n'].to_numpy(dtype=int) == fr['pin'].to_numpy(dtype=int)).all():
        return fr, "pins differ from the yaml's sequence"
    fr['yang_peak'] = t['roi_peak_dff'].to_numpy(dtype=float)
    return fr, 'ok'


def recording_readings(frames: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per presentation, per stimulus) readings of one recording's Yang values (E100's `presentation_readings` and
    `recording_b`'s means, so that set S's values are E100's to the bit): raw, minus_recording, minus_block,
    order_matched, and 'first', the order-matched value of the stimulus's first presentation; per stimulus also the SD
    over presentations of the raw response, the number of presentations and the first block."""
    E100 = _e100()
    o = E100.presentation_readings(frames, 'yang_peak')
    g = o.groupby('stim', sort=False)
    st = pd.DataFrame(dict(odor=g['odor'].first(), conc=g['log10_conc'].first(), n_presentations=g.size(),
                           first_block=g['block'].min()))
    for rd in E100.READINGS_B:
        st[f'yang_{rd}'] = g[rd].mean()
    st['yang_sd_presentations'] = g['raw'].std(ddof=1)
    st['yang_first'] = o[o['rep'] == 0].set_index('stim')['order_matched'].reindex(st.index)
    return o, st.reset_index()


S_ID_COLS = ('n_presentations', 'first_block', 'yang_raw', 'yang_minus_recording', 'yang_minus_block',
             'yang_order_matched', 'yang_sd_presentations')


def s_identity() -> Tuple[dict, pd.DataFrame]:
    """Part 0: the reader (`read_table_k`, `recording_readings`) reproduces E100's saved set S values
    (b_per_stimulus.csv, S_ID_COLS) exactly on every recording and stimulus, from the same Yang table (the file E100's
    b_inputs.csv records)."""
    E100 = _e100()
    res_dir = H.results_dir(E100.EXPERIMENT)
    ref = pd.read_csv(res_dir / 'b_per_stimulus.csv', float_precision='round_trip')
    ref = ref[ref['set'] == 'S']
    inputs = pd.read_csv(res_dir / 'b_inputs.csv')
    e100_tab = {r: Path(p).name for r, k, p in zip(inputs['recording'], inputs['kind'], inputs['path'])
                if k == 'roi_table'}
    rows = recording_rows(list(S_STEMS))
    per = []
    for st, row in rows.iterrows():
        fr, status = read_table_k(row)
        used = yang_table_path(row)
        d = dict(recording=st, table=used.name if used is not None else '', e100_table=e100_tab.get(st, ''),
                 table_status=status)
        d['same_table'] = bool(d['table'] == d['e100_table'])
        if status != 'ok':
            per.append(dict(d, n_stimuli=0, max_abs_diff=np.inf))
            continue
        _, mine = recording_readings(fr)
        r = ref[ref['recording'] == st].set_index('stim')
        m = mine.set_index('stim')
        d['same_stimuli'] = bool(sorted(m.index) == sorted(r.index))
        diffs = {}
        for c in S_ID_COLS:
            a = m[c].reindex(r.index).astype(float).to_numpy()
            b = r[c].astype(float).to_numpy()
            both_nan = np.isnan(a) & np.isnan(b)
            dd = np.where(both_nan, 0.0, np.abs(a - b))
            diffs[c] = float(np.nanmax(np.where(np.isnan(dd), np.inf, dd))) if len(dd) else 0.0
        per.append(dict(d, n_stimuli=int(len(r)), max_abs_diff=max(diffs.values()), **{f'diff_{k}': v for k, v in
                                                                                       diffs.items()}))
    df = pd.DataFrame(per)
    worst = float(df['max_abs_diff'].max()) if len(df) else np.inf
    res = dict(passed=bool(len(df) == len(S_STEMS) and worst <= IDENTITY_TOL and df['same_table'].all()
                           and df.get('same_stimuli', pd.Series([False])).all()),
               n_recordings=int(len(df)), n_stimuli=int(df['n_stimuli'].sum()), max_abs_diff=worst,
               columns=list(S_ID_COLS), tables_differ=df.loc[~df['same_table'], 'recording'].tolist())
    return res, df


# ---------------------------------------------------------------------------
# part B's runner

def b_stimuli() -> List[str]:
    """The 11 B1 labels of a ladder run, in B1's order: the panel's paa, the three weak stimuli under each session's
    estimate, the blank."""
    E98 = _e98()
    want = {REFERENCE, E98.BLANK_LABEL}
    want |= {E98.session_label(o, c, s, E98.DIAG) for o, c in WEAK_STIMULI for s in E98.SESSIONS}
    order = [sp['label'] for sp in E98.b1_specs()]
    out = [x for x in order if x in want]
    if len(out) != 11:
        raise RuntimeError(f"the ladder's stimuli are not all in B1: {sorted(want - set(out))}")
    return out


def ladder_params(p: Dict[str, float], m: float) -> Dict[str, float]:
    """The fit's parameters with APL->KC and APL->bouton scales both times `m` (m = 0: E98's acute APL off, the scales
    0; m = 1: the fit's)."""
    q = dict(p)
    for k in ('wAPLKC_scale', 'wAPLPN_scale'):
        if k not in q:
            raise KeyError(f"{k} not in the fit's parameters")
        if m == 0:
            q[k] = 0.0
        elif m != 1:
            q[k] = float(p[k]) * float(m)
    return q


def ladder_input() -> pd.DataFrame:
    """E98's B1 input restricted to `b_stimuli` (the panel's paa column and ten added stimuli, B1's order)."""
    E98 = _e98()
    orn, _ = E98.b1_input()
    labels = E98.input_labels(orn)
    want = set(b_stimuli())
    keep = [i for i, lab in enumerate(labels) if lab in want]
    return orn.iloc[:, keep]


def ladder_run(wiring: str, orn: Optional[pd.DataFrame] = None):
    """E98's `ExtendedRun` of `wiring` at its fit on `ladder_input` (APL's activity saved)."""
    return _e98().ExtendedRun(wiring, ladder_input() if orn is None else orn)


def run_ladder_point(er, m: float) -> dict:
    """One simulation of a ladder run (`ladder_run`) at `m`: responses, spike counts (columns the B1 labels), APL's
    means per stimulus ('apl': win, pre, pre_settled), its trace, the PN stage's means, and the parameters."""
    E98 = _e98()
    params = ladder_params(er.p, m)
    r, s = er.sim(params)
    gl = pd.Index(list(dict.fromkeys(er.W.columns.get_level_values('glomerulus'))), name='glomerulus')
    o = dict(E98.read_apl_pn(er.mp, er.rv, er.order, er.labels, gl), **E98.relabel(r, s, er.labels))
    o['params'] = params
    return o


def save_run(d: Path, wiring: str, m: float, o: dict) -> None:
    """A ladder run's outputs under `d`: <wiring>_m<m>_{spike_counts,responses,apl_trace,pn_win,pn_pre}.parquet and
    _apl.csv."""
    tag = f'{wiring}_m{m:g}'
    for k in ('responses', 'spike_counts', 'pn_win', 'pn_pre', 'apl_trace'):
        x = o[k].copy()
        x.columns = [str(c) for c in x.columns]
        x.to_parquet(d / f'{tag}_{k}.parquet')
    o['apl'].to_csv(d / f'{tag}_apl.csv')


def runner_identity(wiring: str = 'real', save: bool = True) -> dict:
    """Part 0: the runner (`ladder_run`, `run_ladder_point`) at m = 1 and m = 0 against E98's saved B1 'on' and
    'acute_off' runs of `wiring` on the 11 stimuli: the number of KC x stimulus spike counts (and responses) that
    differ (0 to pass), APL's window, pre-odor and settled means' largest difference beside; the parameters at m = 0
    and 1 against E98's. Saves the two runs under outputs/E104_kc_calyx_vs_apl/B/identity/."""
    E98 = _e98()
    saved = E98.load_b1_run(wiring)
    if saved is None:
        raise RuntimeError(f"{wiring}: no saved B1 run of E98's (b1_roots)")
    t0 = time.time()
    er = ladder_run(wiring)
    labs = list(er.labels)
    if labs != b_stimuli():
        raise RuntimeError(f"the run's stimuli are not the ladder's: {labs}")
    d = local_dir('B/identity') if save else None
    res = dict(wiring=wiring, saved_dir=saved['on']['dir'], stimuli=labs, threads=os.environ.get('OMP_NUM_THREADS'),
               thread_limit=os.environ.get('OMP_THREAD_LIMIT'))
    for m, ct in RUNNER_CONDS:
        t1 = time.time()
        o = run_ladder_point(er, m)
        sv = saved[ct]
        s_new, s_old = o['spike_counts'][labs].astype(float), sv['spike_counts'][labs].astype(float)
        r_new, r_old = o['responses'][labs].astype(float), sv['responses'][labs].astype(float)
        idx_ok = bool(len(s_new) == len(s_old) and np.array_equal(np.asarray(s_new.index), np.asarray(s_old.index)))
        nd = int((s_new.to_numpy() != s_old.to_numpy()).sum()) if idx_ok else -1
        nr = int((r_new.to_numpy() != r_old.to_numpy()).sum()) if idx_ok else -1
        a_new, a_old = o['apl'].loc[labs], sv['apl'].loc[labs]
        apl_diff = {c: float((a_new[c].astype(float) - a_old[c].astype(float)).abs().max())
                    for c in ('win', 'pre', 'pre_settled') if c in a_old.columns}
        e98_params = er.p if ct == 'on' else dict(er.p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
        res[f'm{m:g}'] = dict(saved_tag=ct, index_equal=idx_ok, counts_differ=nd, responses_differ=nr,
                              n_cells=int(s_new.size), apl_max_abs_diff=apl_diff,
                              params_equal_e98=bool(o['params'] == e98_params),
                              wAPLKC_scale=float(o['params']['wAPLKC_scale']),
                              wAPLPN_scale=float(o['params']['wAPLPN_scale']),
                              max_spike_count=float(s_new.to_numpy().max()), seconds=round(time.time() - t1, 1))
        if save:
            save_run(d, wiring, m, o)
        print(f"{time.strftime('%H:%M:%S')} {wiring} m={m:g} against E98's {ct}: {nd} counts differ, {nr} responses "
              f"differ ({time.time() - t1:.0f} s)", flush=True)
    res['passed'] = bool(all(res[f'm{m:g}']['counts_differ'] == 0 and res[f'm{m:g}']['responses_differ'] == 0
                             and res[f'm{m:g}']['params_equal_e98'] for m, _ in RUNNER_CONDS))
    res['seconds'] = round(time.time() - t0, 1)
    return res


def part_0(runs: bool = False, lock: bool = False) -> dict:
    """Part 0 (no response table of set K opened): p0_eligibility.csv, p0_s_identity.csv and p0_identity.json (the
    entry's counts, the lists, set S's reader identity; with `runs` the runner's identity on the real wiring, kept from
    an earlier run otherwise); with `lock`, then the lock (`write_lock`), only if every identity passed. Raises if one
    fails."""
    out = out_dir()
    el = eligibility_table()
    write_table(el, out / 'p0_eligibility.csv')
    counts = eligibility_counts(el, set_l_counts())
    ident = eligibility_identity(el)
    s_res, s_rows = s_identity()
    write_table(s_rows, out / 'p0_s_identity.csv')
    old = _read_json(out / 'p0_identity.json') or {}
    runner = runner_identity('real') if runs else old.get('runner')
    res = dict(caveat=QUALITY_CAVEAT, counts=counts, lists=ident, s_identity=s_res, runner=runner,
               no_k_table_opened=True, definitions_sha256=definitions_sha(part_a_definitions()),
               host=_host(), time=time.strftime('%Y-%m-%d %H:%M'))
    res['passed'] = bool(counts['passed'] and ident['passed'] and s_res['passed'] and runner is not None
                         and runner.get('passed'))
    _write_json(out / 'p0_identity.json', res)
    provenance('0', runs=runs, passed=res['passed'], files={str(p): _sha256(p) for p in (
        H.shared_dir() / 'data' / 'yang_apl_recordings' / 'recordings.csv',
        H.results_dir('E100_apl_weak_odors') / 'b_per_stimulus.csv',
        H.shared_dir() / 'results' / 'E101_apl_vs_kc_claws' / 'c_responses.csv')})
    if not res['passed']:
        raise RuntimeError(f"part 0: an identity fails or has not run (p0_identity.json): counts {counts['passed']}, "
                           f"lists {ident['passed']}, set S {s_res['passed']}, runner "
                           f"{None if runner is None else runner.get('passed')}")
    if lock:
        res['lock'] = dict(path=str(PART_A_LOCK), sha256=write_lock()['sha256'])
    return res


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
    ap.add_argument('--write-lock', action='store_true',
                    help="part 0: write part A's lock once every identity passes (commit it before part A)")
    ap.add_argument('--wirings', default=','.join(ALL_WIRINGS), help='part B: comma-separated wirings')
    ap.add_argument('--threads', type=int, default=None,
                    help='runs (part 0 with --runs, part B): OMP threads; OMP_NUM_THREADS and OMP_THREAD_LIMIT must match')
    a = ap.parse_args()
    simulates = (a.part == '0' and a.runs) or a.part == 'B'
    if simulates:
        if a.threads is None:
            ap.error('runs need --threads')
        for v in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT'):
            if os.environ.get(v) != str(a.threads):
                ap.error(f'run with OMP_NUM_THREADS={a.threads} OMP_THREAD_LIMIT={a.threads} in the environment')
    if a.part == '0':
        print(json.dumps(part_0(runs=a.runs, lock=a.write_lock), indent=1, default=_json_default), flush=True)
    elif a.part == 'A':
        print(json.dumps(part_a(), indent=1, default=_json_default), flush=True)
    elif a.part == 'B':
        print(json.dumps(part_b(tuple(x for x in a.wirings.split(',') if x)), indent=1, default=_json_default),
              flush=True)
    else:
        print(json.dumps(part_c(), indent=1, default=_json_default), flush=True)


if __name__ == '__main__':
    main()
