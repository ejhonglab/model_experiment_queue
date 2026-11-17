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
- D (added with the write-up; no simulation, from parts A's and B's saved tables): the counts outcomes (e) and (f) ask
  for. X1: per fall and reading, the new flies (KC flies alone, or the same number in each set) at which a class would
  be reached by the typical draw and with 80 % power if the present Q and between-fly SDs were the truth
  (`flies_to_decide`, the normal approximation of the split rule and of Q's interval, as E98's `flies_needed`; checked
  against the split rule itself, E103's `least_detectable`). Y1 (and Y2): the flies' ratio above which an odor is too
  weak (`y_thresholds`), and the new blocked recordings (the control flies as they are) or new flies of each genotype at
  which it would be (`y_counts`). d_x1_counts.csv, d_y_counts.csv, d_summary.json.

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
- X3's count is literal: at least 4 KC flies with a positive response at 10^-7. 2022-09-18 imaged paa and ms only down
  to 10^-6, so for those two odors 4 flies have data, and leaving out any other fly (partner (6)) leaves 3: that
  partner cannot give "holds".
- A5's falls (described): each odor's lower concentrations of the SytGCaMP7f design over its highest, against set H
  (E100's reading); E101's Q pairs set L's APL channel with its own KC channel at b = 0 and b_G.
- Part B's flies: the KC rates are recomputed from the weak-odor recordings' trials (diag_mix_long.csv, trial mean above
  1 z; `kc_identity` checks them against the saved fractions and E49's means); the ratio is the recording means' (as
  E49's), its interval from bootstrap draws of the 12 control flies (a drawn fly brings its recordings) with each of the
  three blocked recordings left out in turn, the hull. The flies' APL over paa at 10^-3 is E100's set W (its
  score_weak_table.csv, arithmetic means over 8 flies). "Relative to the panel": the family's ratio over its own ratio
  on the 16 panel odors at its fit (mean rates, acute off over on).
- d*: the family's ratio (APL off over on) made non-decreasing in d (pool-adjacent-violators) and inverted linearly;
  an interval end beyond the ladder's top (m = 2) is +inf (censored). Y1's strength is the flies' APL interval, Y2's the
  family's own d at m = 1 under the same reading of d (the blank's level, or the settled level).
- B3's widening by E103's vial SD (beside): each side of the flies' APL interval combined in quadrature with 1.645 x
  0.461 on the log2 scale.

Usage (from the repository root):
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part 0
    env OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E104_kc_calyx_vs_apl --part 0 --runs --threads 11      # with the runner's identities
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part 0 --write-lock   # once, before any of set K is read
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part A
    env OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E104_kc_calyx_vs_apl --part B --threads 11
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part C
    ./venv/bin/python -u -m experiments.E104_kc_calyx_vs_apl --part D

Outputs: results/E104_kc_calyx_vs_apl/ in the shared directory; part B's runs under outputs/E104_kc_calyx_vs_apl/B/,
archived to data/model_outputs/E104_kc_calyx_vs_apl/ (m and the absolute APL scales of every run in params.csv).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
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

READINGS_A = ('carrying', 'block_solvent', 'first_presentation')     # the value columns of a_per_stimulus.csv
KC_SETS_CARRYING = ('K6F',)
APL_SETS_CARRYING = ('S', 'L')
E101_C = ('results', 'E101_apl_vs_kc_claws', 'c_responses.csv')
E101_A3_FAMILY = ('results', 'E101_apl_vs_kc_claws', 'a3_family.csv')
E100_B = ('results', 'E100_apl_weak_odors')
E101_READINGS = {'carrying': 'carrying', 'block_solvent': 'block', 'first_presentation': 'first_two'}
X2_ESTIMATE_SESSIONS = ('kiwi-control', 'megamat', 'validation2')
HALLEM_BESIDE = (('pa_-6_over_-4', 'pa [hallem] @ -6', 'pa [hallem] @ -4'),
                 ('ms_-4_over_-2', 'ms [hallem] @ -4', 'ms [hallem] @ -2'))


def _ci(x: np.ndarray) -> Tuple[float, float]:
    lo, hi = np.percentile(x, [100 * (1 - CI_LEVEL) / 2, 100 * (1 + CI_LEVEL) / 2])
    return float(lo), float(hi)


def read_flies_k_s(stems: Optional[Sequence[str]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Part A's reading of sets K and S from Yang's tables (`read_table_k`; the lock is required for set K): per
    recording and stimulus the readings (READINGS_A as columns 'carrying' (order-matched), 'block_solvent' (minus its
    block's solvent), 'first_presentation' and 'raw'), with set, fly, calyx, plane and marked; and every file read with
    its sha256 and the table's status."""
    stems = list(K6F + K1020 + K7F + S_STEMS) if stems is None else list(stems)
    rows = recording_rows(stems)
    out, files = [], []
    for st, row in rows.iterrows():
        fr, status = read_table_k(row)
        p = yang_table_path(row)
        files.append(dict(recording=st, set=row['role'], kind='yang_table', path=str(p) if p else '',
                          sha256=_sha256(p) if p else '', status=status, n_presentations=int(len(fr))))
        y = yaml_path(row)
        files.append(dict(recording=st, set=row['role'], kind='yaml', path=str(y), sha256=_sha256(y), status='',
                          n_presentations=int(len(fr))))
        if status != 'ok':
            continue
        plane = float(row['plane_depth_um']) if pd.notna(row.get('plane_depth_um')) else xml_plane_um(row)
        _, st_tab = recording_readings(fr)
        for x in st_tab.to_dict('records'):
            out.append(dict(set=row['role'], recording=st, fly=row['fly_id'], calyx=row['calyx'], plane_um=plane,
                            plane_uncertain=st in PLANE_UNCERTAIN, marked=MARKED.get(st, ''), odor=x['odor'],
                            conc=float(x['conc']), stim=x['stim'], carrying=x['yang_order_matched'],
                            block_solvent=x['yang_minus_block'], first_presentation=x['yang_first'],
                            raw=x['yang_raw'], n_presentations=int(x['n_presentations'])))
        print(f"{time.strftime('%H:%M:%S')} {st} ({row['role']}): Yang's table {status}", flush=True)
    return pd.DataFrame(out), pd.DataFrame(files)


def read_set_l(path: Optional[Path] = None) -> pd.DataFrame:
    """Set L from E101's c_responses.csv (group 'L', b = 0): the APL channel's 'apl' under E101's readings mapped to
    READINGS_A (E101_READINGS), and the KC channel's 'kc' at b = 0 and b_G (for E101's Q, A5); fly ids E101's."""
    p = H.shared_dir().joinpath(*E101_C) if path is None else Path(path)
    c = pd.read_csv(p, comment='#')
    c = c[c['group'] == 'L']
    out = []
    for (rec, fly, plane, stim, odor, conc), g in c.groupby(['recording', 'fly', 'plane_um', 'stim', 'odor', 'log10_conc'],
                                                            sort=False):
        d = dict(set='L', recording=rec, fly=str(fly), calyx='', plane_um=float(plane), plane_uncertain=False, marked='',
                 odor=odor, conc=float(conc), stim=stim, raw=np.nan, n_presentations=np.nan)
        for ours, theirs in E101_READINGS.items():
            x = g[(g['reading'] == theirs) & (g['b_mult'] == 0.0)]
            d[ours] = float(x['apl'].iloc[0]) if len(x) else np.nan
            for bm in (0.0, 1.0):
                y = g[(g['reading'] == theirs) & (g['b_mult'] == bm)]
                d[f'kc_{ours}_b{bm:g}'] = float(y['kc'].iloc[0]) if len(y) else np.nan
        out.append(d)
    return pd.DataFrame(out)


def read_set_h() -> pd.DataFrame:
    """Set H from E100's b_per_stimulus.csv (label_unsure False): carrying = yang_order_matched, block_solvent =
    yang_minus_block, raw = yang_raw (A5, beside; first_presentation not saved)."""
    b = pd.read_csv(H.shared_dir().joinpath(*E100_B, 'b_per_stimulus.csv'), float_precision='round_trip')
    b = b[(b['set'] == 'H') & ~b['label_unsure'].astype(bool)]
    return pd.DataFrame(dict(set='H', recording=b['recording'], fly=[_fly_from_e100(f) for f in b['fly']], calyx='',
                             plane_um=b['plane_um'].astype(float), plane_uncertain=False, marked='', odor=b['odor'],
                             conc=b['conc'].astype(float), stim=b['stim'], carrying=b['yang_order_matched'],
                             block_solvent=b['yang_minus_block'], first_presentation=np.nan, raw=b['yang_raw'],
                             n_presentations=b['n_presentations']))


def _fly_from_e100(f: str) -> str:
    d, n = str(f).split('_')
    return f'{d[:4]}-{d[4:6]}-{d[6:8]}_fly{int(n)}'


def recording_fall(g: pd.DataFrame, num: Stim, den: Stim, value: str) -> dict:
    """One recording's fall num / den on `value` (E98's `value_at` on each odor's imaged concentrations of the
    recording): the two values and how each was read; 'log_fall' when both exceed MIN_RESPONSE, 'ratio' (arithmetic)
    when the denominator does."""
    E98 = _e98()
    out = {}
    for k, (o, c) in (('num', num), ('den', den)):
        x = g[(g['odor'] == o) & g[value].notna()]
        v, how = E98.value_at(x['conc'].to_numpy(float), x[value].to_numpy(float), c) if len(x) else (np.nan, 'missing')
        out[k], out[f'{k}_how'] = v, how
    ok = bool(np.isfinite(out['num']) and np.isfinite(out['den']))
    out['ratio'] = out['num'] / out['den'] if ok and out['den'] > MIN_RESPONSE else np.nan
    out['log_fall'] = float(np.log(out['num'] / out['den'])) if ok and out['num'] > MIN_RESPONSE and \
        out['den'] > MIN_RESPONSE else np.nan
    return out


def recording_falls(vals: pd.DataFrame, falls: Sequence[Tuple[str, Stim, Stim]], value: str) -> pd.DataFrame:
    """`recording_fall` for every recording of `vals` and every fall (rows: set, recording, fly, plane, fall, ...)."""
    rows = []
    for (set_, rec), g in vals.groupby(['set', 'recording'], sort=False):
        for fid, num, den in falls:
            r = recording_fall(g, num, den, value)
            rows.append(dict(set=set_, recording=rec, fly=g['fly'].iloc[0], plane_um=g['plane_um'].iloc[0],
                             plane_uncertain=bool(g['plane_uncertain'].iloc[0]), marked=g['marked'].iloc[0],
                             fall=fid, value=value, **r))
    return pd.DataFrame(rows)


def fly_falls(rec_falls: pd.DataFrame, arithmetic: bool = False) -> pd.Series:
    """Per fly the mean of its recordings' natural-log falls (`arithmetic`: of their ratios), recordings without one
    left out; a Series indexed by fly."""
    col = 'ratio' if arithmetic else 'log_fall'
    d = rec_falls.dropna(subset=[col])
    return d.groupby('fly')[col].mean()


def set_value(per_fly: pd.Series, arithmetic: bool = False, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """A set's value from per-fly falls: the geometric mean over flies (per-fly log falls) or the arithmetic mean
    (`arithmetic`), its 90 % interval from `n_boot` bootstrap draws of flies (E79's `fly_bootstrap`), the counts."""
    x = per_fly.dropna().to_numpy(float)
    if len(x) == 0:
        return dict(value=np.nan, lo90=np.nan, hi90=np.nan, n_flies=0)
    from experiments import E79_apl_calyx_constraints as E79
    m, lo, hi = E79.fly_bootstrap(x, n_boot=n_boot, seed=seed)
    f = (lambda v: float(v)) if arithmetic else (lambda v: float(np.exp(v)))
    return dict(value=f(m[0]), lo90=f(lo[0]), hi90=f(hi[0]), n_flies=int(len(x)))


def split_rule(a: pd.Series, b: pd.Series) -> dict:
    """E98's split rule on two sets of flies' mean log falls (E103's `split_share`, the same count as E98's
    `pulse_permutation`): the share of the splits of all flies into sets of the two sizes whose difference of means is at
    least the observed (absolute); 'differ' when below SPLIT_SHARE."""
    from experiments import E103_apl_across_odors_heldout as E103
    a, b = a.dropna(), b.dropna()
    if len(a) == 0 or len(b) == 0:
        return dict(share=np.nan, n_splits=0, differ=False, n_a=int(len(a)), n_b=int(len(b)))
    share, n = E103.split_share(a.to_numpy(float), b.to_numpy(float))
    return dict(share=share, n_splits=n, differ=bool(share < SPLIT_SHARE), n_a=int(len(a)), n_b=int(len(b)))


def q_of_fall(apl: pd.Series, kc: pd.Series, arithmetic: bool = False, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Q = APL's fall over the KCs' (geometric means over flies of per-fly log falls; `arithmetic`: the arithmetic
    means' quotient) with its 90 % interval from joint bootstrap draws (each draw resamples both sets' flies with
    replacement)."""
    a, k = apl.dropna().to_numpy(float), kc.dropna().to_numpy(float)
    if len(a) == 0 or len(k) == 0:
        return dict(q=np.nan, q_lo90=np.nan, q_hi90=np.nan)
    rng = np.random.default_rng(seed)
    ia = rng.integers(0, len(a), size=(n_boot, len(a)))
    ik = rng.integers(0, len(k), size=(n_boot, len(k)))
    if arithmetic:
        q = a.mean() / k.mean()
        qb = a[ia].mean(axis=1) / k[ik].mean(axis=1)
    else:
        q = float(np.exp(a.mean() - k.mean()))
        qb = np.exp(a[ia].mean(axis=1) - k[ik].mean(axis=1))
    lo, hi = _ci(qb)
    return dict(q=float(q), q_lo90=lo, q_hi90=hi)


def fall_class(differ: bool, q: float, q_lo: float, q_hi: float) -> str:
    """The class of a fall (FALL_CLASSES) from the split rule's verdict and Q with its interval."""
    if not np.isfinite(q):
        return 'not decided'
    if differ and q < Q_APL_STEEPER:
        return 'APL steeper'
    if differ and q > Q_KC_STEEPER:
        return 'KCs steeper'
    if not differ and np.isfinite(q_lo) and np.isfinite(q_hi) and Q_ONE_FUNCTION[0] <= q_lo and q_hi <= Q_ONE_FUNCTION[1]:
        return 'one function'
    return 'not decided'


def held_class(by_reading: Dict[str, str], carrying: str = 'carrying') -> str:
    """The carrying reading's class if every other reading in `by_reading` gives the same; else 'not decided'."""
    c = by_reading.get(carrying, 'not decided')
    return c if all(v == c for v in by_reading.values()) else 'not decided'


def x1_verdict(classes: Dict[str, str]) -> str:
    """X1 from the (held) classes of G1, G2 and G3: 'APL steeper' if at least X1_MIN are APL steeper and none KCs
    steeper; 'one function' if at least X1_MIN are one function and none steeper either way; otherwise 'not decided'."""
    v = list(classes.values())
    n = {k: sum(x == k for x in v) for k in FALL_CLASSES}
    if n['APL steeper'] >= X1_MIN and n['KCs steeper'] == 0:
        return 'APL steeper'
    if n['one function'] >= X1_MIN and n['APL steeper'] == 0 and n['KCs steeper'] == 0:
        return 'one function'
    return 'not decided'


def x3_verdict(per_odor_positive: Dict[str, int], min_flies: int = X3_MIN_FLIES, min_odors: int = X3_MIN_ODORS) -> dict:
    """X3 from the number of KC flies whose response at 10^-7 is positive, per odor: an odor is above the solvent with
    at least `min_flies`; X3 'holds' with at least `min_odors` such odors, else 'fails'."""
    above = {o: bool(n >= min_flies) for o, n in per_odor_positive.items()}
    return dict(above=above, n_above=int(sum(above.values())),
                verdict='holds' if sum(above.values()) >= min_odors else 'fails')


def reading_specs(kc_flies: Sequence[str]) -> Dict[str, dict]:
    """The carrying reading and the partners (PARTNERS; partner (6) one reading per KC fly left out): the value
    column, the APL and KC sets, the plane, the marked runs, the KC fly left out."""
    base = dict(value='carrying', apl_sets=APL_SETS_CARRYING, kc_sets=KC_SETS_CARRYING, plane=None,
                without_marked=False, drop_kc_fly=None)
    out = {'carrying': base,
           'set_S_alone': dict(base, apl_sets=('S',)),
           'block_solvent': dict(base, value='block_solvent'),
           'first_presentation': dict(base, value='first_presentation'),
           'plane_18': dict(base, plane=CARRYING_PLANE_PARTNER),
           'without_marked': dict(base, without_marked=True)}
    for f in sorted(kc_flies):
        out[f'kc_fly_left_out:{f}'] = dict(base, drop_kc_fly=f)
    out['k1020_as_kc'] = dict(base, kc_sets=('K6F', 'K1020'))
    return out


def select(vals: pd.DataFrame, sets: Sequence[str], spec: dict, side: str) -> pd.DataFrame:
    """The rows of `vals` of `sets` under a reading's spec (plane, marked runs, the KC fly left out)."""
    v = vals[vals['set'].isin(sets)]
    if spec.get('plane') is not None:
        v = v[np.isclose(v['plane_um'].astype(float), float(spec['plane'])) & ~v['plane_uncertain'].astype(bool)]
    if spec.get('without_marked'):
        v = v[v['marked'].fillna('').astype(str) == '']
    if side == 'kc' and spec.get('drop_kc_fly'):
        v = v[v['fly'] != spec['drop_kc_fly']]
    return v


def fall_reading(vals: pd.DataFrame, fall: Tuple[str, Stim, Stim], spec: dict, arithmetic: bool = False,
                 apl_sets: Optional[Sequence[str]] = None, kc_sets: Optional[Sequence[str]] = None) -> dict:
    """One fall under one reading: each set's value with its interval (`set_value`), the split rule, Q and the class."""
    a = select(vals, apl_sets or spec['apl_sets'], spec, 'apl')
    k = select(vals, kc_sets or spec['kc_sets'], spec, 'kc')
    fa = fly_falls(recording_falls(a, [fall], spec['value']), arithmetic)
    fk = fly_falls(recording_falls(k, [fall], spec['value']), arithmetic)
    va, vk = set_value(fa, arithmetic), set_value(fk, arithmetic)
    sr = split_rule(fa, fk)
    q = q_of_fall(fa, fk, arithmetic)
    cls = fall_class(sr['differ'], q['q'], q['q_lo90'], q['q_hi90']) if not arithmetic else 'described'
    return dict(fall=fall[0], apl=va['value'], apl_lo90=va['lo90'], apl_hi90=va['hi90'], apl_n_flies=va['n_flies'],
                kc=vk['value'], kc_lo90=vk['lo90'], kc_hi90=vk['hi90'], kc_n_flies=vk['n_flies'],
                split_share=sr['share'], n_splits=sr['n_splits'], differ=sr['differ'], **q, cls=cls,
                mean='arithmetic' if arithmetic else 'geometric')


def falls_table(vals: pd.DataFrame) -> pd.DataFrame:
    """A2 (and the beside falls): every fall of FALLS and FALLS_BESIDE under the carrying reading and every partner
    (`reading_specs`); the 10^-7 / 10^-3 falls with arithmetic means (described); 2022-10-20 as an APL fly beside."""
    kc_flies = sorted(vals.loc[vals['set'].isin(KC_SETS_CARRYING), 'fly'].unique())
    specs = reading_specs(kc_flies)
    rows = []
    for fall in FALLS + FALLS_BESIDE:
        arith = fall[0] in WEAK_FALL_IDS
        for name, spec in specs.items():
            rows.append(dict(reading=name, beside=fall[0] not in dict((f, 1) for f, _, _ in FALLS),
                             **fall_reading(vals, fall, spec, arith)))
        rows.append(dict(reading='k1020_as_apl (described)', beside=True,
                         **fall_reading(vals, fall, specs['carrying'], arith, apl_sets=('S', 'L', 'K1020'))))
    return pd.DataFrame(rows)


def x1_table(falls: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Per reading the classes of G1 to G3 and X1 under it; the held class of each fall and the held X1 (`x1_verdict`
    of the held classes)."""
    g = falls[falls['fall'].isin([f for f, _, _ in FALLS]) & ~falls['reading'].str.contains('described')]
    rows = []
    for rd, x in g.groupby('reading', sort=False):
        cls = dict(zip(x['fall'], x['cls']))
        rows.append(dict(reading=rd, **cls, x1=x1_verdict(cls)))
    t = pd.DataFrame(rows)
    held = {f: held_class(dict(zip(t['reading'], t[f]))) for f, _, _ in FALLS}
    return t, dict(held_falls=held, x1=x1_verdict(held), x1_carrying=t.loc[t['reading'] == 'carrying', 'x1'].iloc[0])


def x2_family(path: Optional[Path] = None) -> pd.DataFrame:
    """E101's a3_family.csv, APL on: per wiring and stimulus the claw columns' window means (X2_CLAW_COLUMNS)."""
    p = H.shared_dir().joinpath(*E101_A3_FAMILY) if path is None else Path(path)
    a = pd.read_csv(p, comment='#', float_precision='round_trip')
    return a[a['cond'] == 'on']


def x2_wiring_ratios(fam: pd.DataFrame, num: str, den: str) -> pd.DataFrame:
    """Per wiring the model ratio set of `num` over `den`: each claw version's window means with the blank's window
    mean as the level, through E98's indicator family (`model_ratio_set`); the positive finite ratios kept."""
    E98 = _e98()
    rows = []
    for w, g in fam.groupby('wiring', sort=False):
        g = g.set_index('stimulus')
        for col in X2_CLAW_COLUMNS:
            a1, a2, a0 = (float(g.loc[s, f'{col}_win']) for s in (num, den, _e98().BLANK_LABEL))
            rs = E98.model_ratio_set(a1, a2, a0, INDICATOR_NS, ('baseline', 'evoked'))
            for r in rs.to_dict('records'):
                rows.append(dict(wiring=w, claw=col, a_num=a1, a_den=a2, a_blank=a0, **r))
    return pd.DataFrame(rows)


def x2_classes(ratios: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring E98's `classify` of its ratio set (positive finite ratios) against [lo, hi]; 'no ratio' if none."""
    E98 = _e98()
    rows = []
    for w, g in ratios.groupby('wiring', sort=False):
        r = g['ratio'].to_numpy(float)
        r = r[np.isfinite(r) & (r > 0)]
        lin = g[(g['n'] == 1.0) & (g['form'] == 'evoked')]['ratio']
        rows.append(dict(wiring=w, ratio_min=float(r.min()) if len(r) else np.nan,
                         ratio_max=float(r.max()) if len(r) else np.nan,
                         ratio_linear_min=float(lin.min()), ratio_linear_max=float(lin.max()),
                         n_ratios=int(len(r)), cls=E98.classify(r, lo, hi) if len(r) and np.isfinite(lo) and
                         np.isfinite(hi) and lo > 0 else 'no ratio'))
    return pd.DataFrame(rows)


def x2_table(vals: pd.DataFrame, fam: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    """X2: per reading of set K (`reading_specs`; partner (1) changes no KC value), set K's G1 interval, per session's
    estimate of paa at 10^-5 the classes of the eight wirings and E98's comparison verdict (X2_MIN_WIRINGS), and the
    verdict over the estimates (E100's `ratio_over_estimates`); the held X2. Beside, on Hallem's input: pa 10^-6 /
    10^-4 and ms 10^-4 / 10^-2 against set K's beside falls (carrying reading)."""
    E98, E100 = _e98(), _e100()
    kc_flies = sorted(vals.loc[vals['set'].isin(KC_SETS_CARRYING), 'fly'].unique())
    g1 = FALLS[0]
    per_w, rows, by_reading = [], [], {}
    for rd, spec in reading_specs(kc_flies).items():
        k = select(vals, spec['kc_sets'], spec, 'kc')
        v = set_value(fly_falls(recording_falls(k, [g1], spec['value'])))
        per_est = {}
        for s in X2_ESTIMATE_SESSIONS:
            num = E98.session_label('paa', -5, s, E98.DIAG)
            c = x2_classes(x2_wiring_ratios(fam, num, 'paa'), v['lo90'], v['hi90'])
            per_w.append(c.assign(reading=rd, estimate=s, kc_g1=v['value'], kc_lo90=v['lo90'], kc_hi90=v['hi90']))
            per_est[s] = E98.comparison_verdict(c['cls'], X2_MIN_WIRINGS)
        verdict = E100.ratio_over_estimates(per_est)
        by_reading[rd] = verdict
        rows.append(dict(reading=rd, kc_g1=v['value'], kc_lo90=v['lo90'], kc_hi90=v['hi90'], kc_n_flies=v['n_flies'],
                         **{f'verdict_{s}': x for s, x in per_est.items()}, x2=verdict))
    t = pd.DataFrame(rows)
    hall = []
    spec = reading_specs(kc_flies)['carrying']
    k = select(vals, spec['kc_sets'], spec, 'kc')
    for fid, num, den in HALLEM_BESIDE:
        fall = next(f for f in FALLS_BESIDE if f[0] == fid)
        v = set_value(fly_falls(recording_falls(k, [fall], 'carrying')))
        c = x2_classes(x2_wiring_ratios(fam, num, den), v['lo90'], v['hi90'])
        hall.append(c.assign(fall=fid, source="Hallem's receptor data (another laboratory)", kc=v['value'],
                             kc_lo90=v['lo90'], kc_hi90=v['hi90'],
                             verdict=E98.comparison_verdict(c['cls'], X2_MIN_WIRINGS)))
    return t, pd.concat(per_w + hall, ignore_index=True), dict(x2=held_class(by_reading), x2_carrying=by_reading['carrying'],
                                                               by_reading=by_reading)


def x3_table(vals: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """X3: per reading (partners (2) to (7); (1) changes no KC value) and odor, each KC fly's response at 10^-7 (the
    mean over its recordings that imaged it), the flies with a positive one; X3 per reading and held."""
    kc_flies = sorted(vals.loc[vals['set'].isin(KC_SETS_CARRYING), 'fly'].unique())
    rows, by_reading = [], {}
    for rd, spec in reading_specs(kc_flies).items():
        k = select(vals, spec['kc_sets'], spec, 'kc')
        pos = {}
        for o in ODORS:
            x = k[(k['odor'] == o) & np.isclose(k['conc'], WEAK_CONC)].dropna(subset=[spec['value']])
            per_fly = x.groupby('fly')[spec['value']].mean()
            pos[o] = int((per_fly > 0).sum())
            rows.append(dict(reading=rd, odor=o, n_flies=int(len(per_fly)), n_positive=pos[o],
                             fly_values=';'.join(f'{f}:{v:.4f}' for f, v in per_fly.items())))
        x3 = x3_verdict(pos)
        by_reading[rd] = x3['verdict']
        for r in rows[-len(ODORS):]:
            r.update(above=x3['above'][r['odor']], x3=x3['verdict'])
    return pd.DataFrame(rows), dict(x3=held_class(by_reading), x3_carrying=by_reading['carrying'], by_reading=by_reading)


def a1_strengths(vals: pd.DataFrame) -> pd.DataFrame:
    """A1: per set, odor and concentration the carrying response (arithmetic mean over flies of each fly's mean over
    its recordings) with its 90 % interval from bootstrap draws of flies, and the raw dF/F's mean beside."""
    rows = []
    for (set_, o, c), g in vals.groupby(['set', 'odor', 'conc']):
        pf = g.groupby('fly')['carrying'].mean().dropna()
        v = set_value(pf, arithmetic=True)
        raw = g.groupby('fly')['raw'].mean().dropna()
        rows.append(dict(set=set_, odor=o, conc=float(c), response=v['value'], lo90=v['lo90'], hi90=v['hi90'],
                         n_flies=v['n_flies'], n_recordings=int(g['recording'].nunique()),
                         raw_dff=float(raw.mean()) if len(raw) else np.nan,
                         raw_dff_max_recording=float(g['raw'].max()) if g['raw'].notna().any() else np.nan))
    return pd.DataFrame(rows)


def a5_falls() -> Tuple[Tuple[str, Stim, Stim], ...]:
    """A5's falls (K7F against set H): each odor's lower concentrations of the SytGCaMP7f design over its highest."""
    design = {'2h': (-6.771, -5.34, -4.386, -3.908), 'ma': (-7.248, -5.817, -4.863, -3.908),
              'pa': (-6.194, -5.495, -4.796, -4.097), 'paa': (-4.863, -4.386, -3.908, -3.431)}
    out = []
    for o in A5_ODORS:
        cs = design[o]
        for c in cs[:-1]:
            out.append((f'{o}_{c:g}_over_{cs[-1]:g}', (o, c), (o, cs[-1])))
    return tuple(out)


def a5_tables(vals: pd.DataFrame) -> pd.DataFrame:
    """A5 (described): K7F against set H on A5's falls; pa / ms and paa / ms at 10^-3 per set (K6F, S, L); E101's Q in
    the dual-colour flies (set L's APL channel against its KC channel at b = 0 and b_G) on G1 to G3."""
    spec = reading_specs([])['carrying']
    rows = []
    for fall in a5_falls():
        rows.append(dict(part='K7F vs H', **fall_reading(vals, fall, spec, apl_sets=('H',), kc_sets=('K7F',))))
    for fid, num, den in A5_ACROSS:
        for set_ in ('K6F', 'S', 'L'):
            v = set_value(fly_falls(recording_falls(vals[vals['set'] == set_], [(fid, num, den)], 'carrying')))
            rows.append(dict(part=f'across odors, {set_}', fall=fid, apl=v['value'], apl_lo90=v['lo90'],
                             apl_hi90=v['hi90'], apl_n_flies=v['n_flies'], cls='described'))
    lv = vals[vals['set'] == 'L']
    for bm in (0.0, 1.0):
        kcv = lv.assign(carrying=lv[f'kc_carrying_b{bm:g}'], set='L_kc')
        both = pd.concat([lv, kcv], ignore_index=True)
        for fall in FALLS:
            rows.append(dict(part=f"E101's Q, dual-colour flies (KC at b = {bm:g} b_G)",
                             **fall_reading(both, fall, spec, apl_sets=('L',), kc_sets=('L_kc',))))
    return pd.DataFrame(rows)


def part_a() -> dict:
    """Part A (A1 to A5): refuses to run until `part_a_locked()`; reads sets K and S from Yang's tables, sets L and H
    from E101's and E100's saved tables; writes a_per_stimulus.csv, a_inputs.csv, a1_strengths.csv, a2_falls.csv,
    a2_x1.csv, a3_x2.csv, a3_x2_wirings.csv, a4_x3.csv, a5_beside.csv and a_summary.json."""
    if not part_a_locked():
        raise PermissionError(f"part A's lock ({PART_A_LOCK.name}) is not committed or does not match the definitions")
    lock = json.loads(PART_A_LOCK.read_text())
    out = out_dir()
    t0 = time.time()
    ks, files = read_flies_k_s()
    vals = pd.concat([ks, read_set_l(), read_set_h()], ignore_index=True)
    for p, k in ((H.shared_dir().joinpath(*E101_C), 'E101 c_responses'),
                 (H.shared_dir().joinpath(*E100_B, 'b_per_stimulus.csv'), 'E100 b_per_stimulus'),
                 (H.shared_dir().joinpath(*E101_A3_FAMILY), 'E101 a3_family')):
        files = pd.concat([files, pd.DataFrame([dict(recording='', set='', kind=k, path=str(p), sha256=_sha256(p),
                                                     status='', n_presentations=np.nan)])], ignore_index=True)
    write_table(vals, out / 'a_per_stimulus.csv')
    write_table(files, out / 'a_inputs.csv')
    a1 = a1_strengths(vals)
    write_table(a1, out / 'a1_strengths.csv')
    falls = falls_table(vals)
    write_table(falls, out / 'a2_falls.csv')
    x1t, x1 = x1_table(falls)
    write_table(x1t, out / 'a2_x1.csv')
    x2t, x2w, x2 = x2_table(vals, x2_family())
    write_table(x2t, out / 'a3_x2.csv')
    write_table(x2w, out / 'a3_x2_wirings.csv')
    x3t, x3 = x3_table(vals)
    write_table(x3t, out / 'a4_x3.csv')
    a5 = a5_tables(vals)
    write_table(a5, out / 'a5_beside.csv')
    st = files[files['kind'] == 'yang_table']
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED, lock_sha256=lock['sha256'],
                lock_file_sha256=_sha256(PART_A_LOCK),
                tables=dict(n=int(len(st)), status=st['status'].value_counts().to_dict(),
                            not_ok=st.loc[st['status'] != 'ok', 'recording'].tolist()),
                n_flies={s: int(vals.loc[vals['set'] == s, 'fly'].nunique()) for s in vals['set'].unique()},
                x1=x1, x2=x2, x3=x3, seconds=round(time.time() - t0, 1))
    _write_json(out / 'a_summary.json', summ)
    provenance('A', lock_sha256=lock['sha256'], lock_file_sha256=_sha256(PART_A_LOCK),
               inputs_sha256={r['kind'] if not r['recording'] else f"{r['recording']}:{r['kind']}": r['sha256']
                              for r in files.to_dict('records')})
    return summ


# ---------------------------------------------------------------------------
# part B: the KCs' somata at the weak odors

LADDER_SUB = 'B/ladder'
E101_A3_ARCHIVE = ('data', 'model_outputs', 'E101_apl_vs_kc_claws', 'A3')
WEAK_KC_PANEL = ('data', 'yang_weak_odor_kc', 'diag_mix')
E49_FRACTIONS = ('results', 'E72_wiring_population', 'low_conc', 'data_fractions.csv')
E98_A3_RATIOS = ('results', 'E98_apl_response_strength', 'a3_ratios.csv')
KC_RESPONSE_Z = 1.0           # the KC data's criterion: trial mean above 1 z (E49's)
FIRST_TRIALS = 3              # partner: the first three trials
SOLVENT_ODOR = 'pfo @ 0'
Y_READINGS = ('carrying', 'first_three_trials', 'solvent_subtracted', 'relative_to_panel', 'one_spike', 'settled')
VIAL_SD_LOG2 = 0.461          # E103's vial SD (beside, waits on question 22)
NORMAL_Z90 = 1.6449


def flies_stim_label(odor: str) -> str:
    """The KC data's and set W's label of a weak odor ('2h @ -7')."""
    o, c = FLIES_WEAK[odor]
    return f'{o} @ {c:g}'


def model_label(odor: str, session: str) -> str:
    """The B1 label of a weak odor's input under a session's estimate."""
    E98 = _e98()
    o, c = next(s for s in WEAK_STIMULI if s[0] == odor)
    return E98.session_label(o, c, session, E98.DIAG)


def saved_run(wiring: str) -> Dict[str, dict]:
    """A wiring's saved B1 run (m = 1 'on', m = 0 'acute_off'): E98's for the real wiring, E101's archived A3 runs for
    the other seven (as E98's `load_b1_all` reads them)."""
    E98 = _e98()
    roots = None if wiring == 'real' else [H.shared_dir().joinpath(*E101_A3_ARCHIVE)]
    run = E98.load_b1_run(wiring, roots)
    if run is None:
        raise RuntimeError(f'{wiring}: no saved B1 run')
    return run


def compare_runs(o: dict, sv: dict, labs: Sequence[str]) -> dict:
    """Spike counts and responses of a ladder run against a saved run on `labs`: the cells that differ."""
    s_new, s_old = o['spike_counts'][labs].astype(float), sv['spike_counts'][labs].astype(float)
    r_new, r_old = o['responses'][labs].astype(float), sv['responses'][labs].astype(float)
    ok = bool(len(s_new) == len(s_old) and np.array_equal(np.asarray(s_new.index), np.asarray(s_old.index)))
    return dict(index_equal=ok, counts_differ=int((s_new.to_numpy() != s_old.to_numpy()).sum()) if ok else -1,
                responses_differ=int((r_new.to_numpy() != r_old.to_numpy()).sum()) if ok else -1,
                apl_win_max_abs_diff=float((o['apl'].loc[labs, 'win'] - sv['apl'].loc[labs, 'win']).abs().max()))


def run_ladder(wirings: Sequence[str] = ALL_WIRINGS, force: bool = False) -> pd.DataFrame:
    """The ladder's runs (resumable): per wiring, m in M_LADDER at its fit (`ladder_params`) on the 11 stimuli, saved
    under outputs/E104_kc_calyx_vs_apl/B/ladder/ (`save_run`) with a runs.csv row each (m, the absolute APL scales, the
    fit's, every parameter, commits, host, threads). Each wiring other than the real one (checked in part 0) is first
    run at m = 1 against its saved run (<wiring>_m1_identity.json); a difference stops the ladder."""
    import gc
    d = local_dir(LADDER_SUB)
    commits = {k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}
    if commits.get('queue_dirty') or commits.get('pyolfsysm_dirty'):
        raise RuntimeError(f'part B runs from clean checkouts only: {commits}')
    for w in wirings:
        todo = [m for m in M_LADDER if force or not (d / f'{w}_m{m:g}_apl.csv').exists()]
        ident_path = d / f'{w}_m1_identity.json'
        need_ident = w != 'real' and not ident_path.exists()
        if not todo and not need_ident:
            continue
        t0 = time.time()
        er = ladder_run(w)
        labs = list(er.labels)
        if labs != b_stimuli():
            raise RuntimeError(f"{w}: the run's stimuli are not the ladder's")
        if need_ident:
            o = run_ladder_point(er, 1.0)
            ident = dict(wiring=w, **compare_runs(o, saved_run(w)['on'], labs), host=_host(),
                         time=time.strftime('%Y-%m-%d %H:%M'))
            ident['passed'] = bool(ident['counts_differ'] == 0 and ident['responses_differ'] == 0)
            _write_json(ident_path, ident)
            print(f"{time.strftime('%H:%M:%S')} {w} m=1 against its saved run: {ident['counts_differ']} counts differ",
                  flush=True)
            if not ident['passed']:
                raise RuntimeError(f'{w}: the runner at m = 1 does not reproduce the saved run ({ident})')
        for m in todo:
            t1 = time.time()
            o = run_ladder_point(er, m)
            save_run(d, w, m, o)
            p = o['params']
            row = dict(wiring=w, m=m, wAPLKC_scale=float(p['wAPLKC_scale']), wAPLPN_scale=float(p['wAPLPN_scale']),
                       fit_wAPLKC_scale=float(er.p['wAPLKC_scale']), fit_wAPLPN_scale=float(er.p['wAPLPN_scale']),
                       params_json=json.dumps({k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v))
                                               for k, v in p.items()}, sort_keys=True),
                       n_stimuli=len(labs), stimuli=';'.join(labs), max_spike_count=float(o['spike_counts'].to_numpy().max()),
                       runtime_s=round(time.time() - t1, 1), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                       threads=os.environ.get('OMP_NUM_THREADS'), thread_limit=os.environ.get('OMP_THREAD_LIMIT'),
                       **commits)
            pd.DataFrame([row]).to_csv(d / 'runs.csv', mode='a', header=not (d / 'runs.csv').exists(), index=False)
            print(f"{time.strftime('%H:%M:%S')} {w} m={m:g}: {time.time() - t1:.0f} s", flush=True)
            del o
        del er
        gc.collect()
        print(f"{time.strftime('%H:%M:%S')} {w}: ladder done in {time.time() - t0:.0f} s", flush=True)
    return pd.read_csv(d / 'runs.csv')


def _run_rows(w: str, m: float, src: str, apl: pd.DataFrame, resp: pd.DataFrame, counts: pd.DataFrame,
              labs: Sequence[str]) -> List[dict]:
    rows = []
    for s in labs:
        c = counts[s].astype(float)
        rows.append(dict(wiring=w, m=m, source=src, stimulus=s, win=float(apl.loc[s, 'win']), pre=float(apl.loc[s, 'pre']),
                         pre_settled=float(apl.loc[s, 'pre_settled']) if 'pre_settled' in apl.columns else np.nan, rate=float(resp[s].astype(float).mean()),
                         rate_1spike=float((c >= 1).mean()), max_spike_count=float(c.max())))
    return rows


def ladder_model_table(wirings: Sequence[str] = ALL_WIRINGS) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The family's runs at every m: the ladder's (outputs/.../B/ladder/) and m = 1 / m = 0 from the saved runs
    (`saved_run`); per wiring, m and stimulus APL's window, pre-odor and settled means, the KCs' response rate (two or
    more spikes), the rate at one spike or more, and the largest spike count. And per wiring its panel ratio (the 16
    panel odors' mean response rate, acute off over on, at the fit)."""
    d = local_dir(LADDER_SUB)
    labs = b_stimuli()
    rows, panel = [], []
    for w in wirings:
        sv = saved_run(w)
        for m, ct in ((1.0, 'on'), (0.0, 'acute_off')):
            rows += _run_rows(w, m, f"saved:{sv[ct]['dir']}", sv[ct]['apl'], sv[ct]['responses'], sv[ct]['spike_counts'],
                              labs)
        p16 = list(H.APL16_ODORS)
        on, off = (sv[ct]['responses'][p16].astype(float).to_numpy().mean() for ct in ('on', 'acute_off'))
        on1, off1 = ((sv[ct]['spike_counts'][p16].astype(float) >= 1).to_numpy().mean() for ct in ('on', 'acute_off'))
        panel.append(dict(wiring=w, panel_rate_on=on, panel_rate_off=off, panel_ratio=off / on,
                          panel_rate_on_1spike=on1, panel_rate_off_1spike=off1, panel_ratio_1spike=off1 / on1))
        for m in M_LADDER:
            tag = f'{w}_m{m:g}'
            if not (d / f'{tag}_apl.csv').exists():
                raise RuntimeError(f'{tag}: the ladder run is missing (run part B)')
            apl = pd.read_csv(d / f'{tag}_apl.csv', index_col=0)
            rows += _run_rows(w, m, 'ladder', apl, pd.read_parquet(d / f'{tag}_responses.parquet'),
                              pd.read_parquet(d / f'{tag}_spike_counts.parquet'), labs)
    return pd.DataFrame(rows).sort_values(['wiring', 'stimulus', 'm'], kind='stable').reset_index(drop=True), \
        pd.DataFrame(panel)


def add_d(model: pd.DataFrame) -> pd.DataFrame:
    """d(m) at every run: m times APL's evoked activity at m (window mean minus the blank's at the same m: 'd'; minus
    the stimulus's settled pre-odor mean: 'd_settled'), over paa's at m = 1; and the response-rate ratios, APL off (m =
    0) over APL on at m ('ratio', 'ratio_1spike')."""
    E98 = _e98()
    x = model.copy()
    key = ['wiring', 'm']
    blank = x.loc[x['stimulus'] == E98.BLANK_LABEL, key + ['win']].rename(columns={'win': 'blank_win'})
    x = x.merge(blank, on=key, how='left', validate='many_to_one')
    x['evoked'] = x['win'] - x['blank_win']
    x['evoked_settled'] = x['win'] - x['pre_settled']
    ref = x[(x['stimulus'] == REFERENCE) & (x['m'] == 1.0)].set_index('wiring')
    x['d'] = x['m'] * x['evoked'] / x['wiring'].map(ref['evoked']).astype(float)
    x['d_settled'] = x['m'] * x['evoked_settled'] / x['wiring'].map(ref['evoked_settled']).astype(float)
    off = x[x['m'] == 0.0].set_index(['wiring', 'stimulus'])
    idx = pd.MultiIndex.from_arrays([x['wiring'], x['stimulus']])
    x['rate_off'] = off['rate'].reindex(idx).to_numpy(float)
    x['rate_off_1spike'] = off['rate_1spike'].reindex(idx).to_numpy(float)
    x['ratio'] = x['rate_off'] / x['rate']
    x['ratio_1spike'] = x['rate_off_1spike'] / x['rate_1spike']
    return x


def isotonic(y: Sequence[float], increasing: bool = True) -> np.ndarray:
    """Pool-adjacent-violators: the least-squares monotone (non-decreasing, or non-increasing) fit of `y` in order."""
    v = np.asarray(y, dtype=float)
    if not increasing:
        return -isotonic(-v, True)
    blocks = []                          # [sum, count]
    for val in v:
        blocks.append([val, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s, n = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += n
    return np.concatenate([np.full(n, s / n) for s, n in blocks]) if blocks else v


def interp_monotone(x: Sequence[float], y: Sequence[float], x0: float, increasing: bool = True) -> float:
    """y at x0 by linear interpolation between the points sorted by x, y first made monotone (`isotonic`); NaN outside
    the points' range."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) == 0 or not (x.min() <= x0 <= x.max()):
        return np.nan
    o = np.argsort(x, kind='stable')
    return float(np.interp(x0, x[o], isotonic(y[o], increasing)))


def _inverse(x: np.ndarray, f: np.ndarray, y: float, side: str) -> float:
    """For f non-decreasing in x (sorted), the x at which the piecewise-linear f reaches y: the first ('left') or the last
    ('right') such x; -inf below f's range, +inf above it."""
    if y < f[0]:
        return -np.inf
    if y > f[-1]:
        return np.inf
    if side == 'left':
        i = int(np.argmax(f >= y))
        if i == 0 or f[i] == f[i - 1]:
            return float(x[i])
        return float(x[i - 1] + (y - f[i - 1]) * (x[i] - x[i - 1]) / (f[i] - f[i - 1]))
    j = len(f) - 1 - int(np.argmax(f[::-1] <= y))
    if j == len(f) - 1 or f[j + 1] == f[j]:
        return float(x[j])
    return float(x[j] + (y - f[j]) * (x[j + 1] - x[j]) / (f[j + 1] - f[j]))


def d_star(d: Sequence[float], ratio: Sequence[float], target: float, lo: Optional[float] = None,
           hi: Optional[float] = None) -> dict:
    """The strength at which the family's response-rate ratio (APL off over on; made non-decreasing in d) equals
    `target`, with its interval from [lo, hi] (the lower end the first d reaching lo, the upper the last d not above
    hi); +inf where the ratio stays below the value over the ladder (d* beyond m = 2), 0 below the ratio at m = 0 (1)."""
    x, y = np.asarray(d, float), np.asarray(ratio, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    o = np.argsort(x, kind='stable')
    x, f = x[o], isotonic(y[o], True)
    if len(x) == 0:
        return dict(d_star=np.nan, d_star_lo=np.nan, d_star_hi=np.nan, censored='')
    pt = 0.5 * (_inverse(x, f, target, 'left') + _inverse(x, f, target, 'right'))
    a = _inverse(x, f, target if lo is None else lo, 'left')
    b = _inverse(x, f, target if hi is None else hi, 'right')
    def clip(v):
        return v if v == np.inf else max(v, 0.0)
    cens = 'above the ladder' if np.isinf(b) and b > 0 else ''
    return dict(d_star=clip(pt), d_star_lo=clip(a), d_star_hi=clip(b), censored=cens)


def y_class(strength: Tuple[float, float], dstar: Tuple[float, float]) -> str:
    """'too weak' if the strength's interval lies wholly below d*'s interval, 'too strong' if wholly above, else
    'inside' (a point strength has lo == hi); 'no d*' when d* is not defined."""
    s_lo, s_hi = strength
    d_lo, d_hi = dstar
    if np.isnan(d_lo) or np.isnan(d_hi):
        return 'no d*'
    if s_hi < d_lo:
        return 'too weak'
    if s_lo > d_hi:
        return 'too strong'
    return 'inside'


def y_over_wirings(classes: Sequence[str], min_wirings: int = X2_MIN_WIRINGS) -> str:
    """A class under an estimate: the class at least `min_wirings` of the 8 wirings give, else 'not decided'."""
    c = list(classes)
    return next((k for k in Y_CLASSES if sum(x == k for x in c) >= min_wirings), 'not decided')


def y_over_estimates(per_estimate: Dict[str, str]) -> str:
    """Over the estimates (E100's rule): inside if inside under any; too weak or too strong only if so under all."""
    v = list(per_estimate.values())
    if 'inside' in v:
        return 'inside'
    for k in ('too weak', 'too strong'):
        if v and all(x == k for x in v):
            return k
    return 'not decided'


def y_verdict(per_odor: Dict[str, str], min_odors: int = Y_MIN_ODORS) -> str:
    """Y1 / Y2 over the three odors: 'too weak' if at least two are too weak and none too strong; 'too strong' the
    reverse; 'inside' if at least two are inside; otherwise 'not decided'."""
    v = list(per_odor.values())
    n = {k: sum(x == k for x in v) for k in Y_CLASSES}
    if n['too weak'] >= min_odors and n['too strong'] == 0:
        return 'too weak'
    if n['too strong'] >= min_odors and n['too weak'] == 0:
        return 'too strong'
    if n['inside'] >= min_odors:
        return 'inside'
    return 'not decided'


def weak_kc_rates(first_n: Optional[int] = None, path: Optional[Path] = None) -> pd.DataFrame:
    """The weak-odor KC recordings' (diag_mix) response rates per recording and odor (the three weak odors and the
    solvent): the fraction of KCs whose trial mean (over all trials, or the first `first_n`) exceeds KC_RESPONSE_Z,
    with the genotype group and the fly (date, fly number)."""
    p = H.shared_dir().joinpath(*WEAK_KC_PANEL, 'diag_mix_long.csv') if path is None else Path(path)
    want = [flies_stim_label(o) for o in FLIES_WEAK] + [SOLVENT_ODOR]
    t = pd.read_csv(p, usecols=['date', 'fly_num', 'roi', 'odor', 'repeat', 'value', 'recording', 'genotype_group'])
    t = t[t['odor'].isin(want)]
    if first_n is not None:
        t = t[t['repeat'] < first_n]
    m = t.groupby(['recording', 'genotype_group', 'date', 'fly_num', 'odor', 'roi'], sort=False)['value'].mean()
    r = (m > KC_RESPONSE_Z).groupby(level=['recording', 'genotype_group', 'date', 'fly_num', 'odor']).mean()
    r = r.rename('rate').reset_index()
    r['fly'] = r['date'].astype(str) + '_fly' + r['fly_num'].astype(str)
    return r


def flies_rate_ratio(rates: pd.DataFrame, odor_label: str, solvent: bool = False, n_boot: int = N_BOOT,
                     seed: int = SEED) -> dict:
    """The flies' response-rate ratio, APL output blocked (TNT) over control flies, at one odor (recording means, as
    E49's), with its interval: `n_boot` bootstrap draws of the control flies (a drawn fly brings its recordings), each of
    the blocked recordings left out in turn, the 90 % interval of each, their hull. `solvent`: each recording's solvent
    rate subtracted first."""
    x = rates[rates['odor'] == odor_label].set_index('recording')
    if solvent:
        s = rates[rates['odor'] == SOLVENT_ODOR].set_index('recording')['rate']
        x = x.assign(rate=x['rate'] - s.reindex(x.index).to_numpy(float))
    c = x[x['genotype_group'] == 'control']
    b = x[x['genotype_group'] == 'apl_output_blocked']
    ctrl, tnt = float(c['rate'].mean()), float(b['rate'].mean())
    flies = sorted(c['fly'].unique())
    by_fly = [c.loc[c['fly'] == f, 'rate'].to_numpy(float) for f in flies]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(flies), size=(n_boot, len(flies)))
    sums = np.array([v.sum() for v in by_fly])
    ns = np.array([len(v) for v in by_fly])
    ctrl_b = sums[idx].sum(axis=1) / ns[idx].sum(axis=1)
    ivs = []
    for r in b.index:
        loo = float(b.drop(index=r)['rate'].mean())
        ivs.append(_ci(loo / ctrl_b))
    lo, hi = min(a for a, _ in ivs), max(z for _, z in ivs)
    return dict(stimulus=odor_label, ctrl=ctrl, tnt=tnt, ratio=tnt / ctrl, lo90=lo, hi90=hi, n_ctrl_recordings=int(len(c)),
                n_ctrl_flies=len(flies), n_tnt_recordings=int(len(b)), n_tnt_flies=int(b['fly'].nunique()),
                solvent_subtracted=solvent)


def flies_apl_ratios() -> pd.DataFrame:
    """The flies' APL at the weak odors over paa at 10^-3 (set W, E100's score_weak_table.csv: arithmetic means over 8
    flies, 90 % intervals)."""
    w = pd.read_csv(H.shared_dir().joinpath(*E100_B, 'score_weak_table.csv'))
    rows = []
    for o in FLIES_WEAK:
        x = w[(w['set'] == 'W') & (w['stimulus'] == flies_stim_label(o))].iloc[0]
        rows.append(dict(odor=o, stimulus=flies_stim_label(o), apl=float(x['flies_ratio']), apl_lo90=float(x['flies_lo90']),
                         apl_hi90=float(x['flies_hi90']), n_flies=int(x['n_flies']), mean=x['flies_mean']))
    return pd.DataFrame(rows)


def widened(lo: float, mid: float, hi: float, sd_log2: float = VIAL_SD_LOG2) -> Tuple[float, float]:
    """An interval widened by a vial SD on the log2 scale (each side's half-width combined in quadrature with
    NORMAL_Z90 x sd_log2)."""
    a, b = np.log2(mid) - np.log2(lo), np.log2(hi) - np.log2(mid)
    w = NORMAL_Z90 * sd_log2
    return float(2 ** (np.log2(mid) - np.hypot(a, w))), float(2 ** (np.log2(mid) + np.hypot(b, w)))


def flies_targets() -> pd.DataFrame:
    """The flies' KC response-rate ratios per weak odor under every Y reading that changes them (carrying, the first
    three trials, the solvent subtracted, relative to the panel's PANEL_KC_RATIO), with the control flies' rate."""
    rows = []
    all_t, first = weak_kc_rates(), weak_kc_rates(FIRST_TRIALS)
    for o in FLIES_WEAK:
        lab = flies_stim_label(o)
        for rd, rates, solv in (('carrying', all_t, False), ('first_three_trials', first, False),
                                ('solvent_subtracted', all_t, True)):
            rows.append(dict(odor=o, reading=rd, **flies_rate_ratio(rates, lab, solv)))
        c = dict(rows[-3])
        rows.append(dict(c, reading='relative_to_panel', ratio=c['ratio'] / PANEL_KC_RATIO,
                         lo90=c['lo90'] / PANEL_KC_RATIO, hi90=c['hi90'] / PANEL_KC_RATIO))
        for rd in ('one_spike', 'settled'):
            rows.append(dict(c, reading=rd))
    return pd.DataFrame(rows)


def kc_identity(rates: pd.DataFrame) -> dict:
    """The recomputed rates (all trials) against data/yang_weak_odor_kc's response_fraction_per_recording.csv
    (frac_trialmean_gt1) and E49's recording means (data_fractions.csv)."""
    ref = pd.read_csv(H.shared_dir().joinpath(*WEAK_KC_PANEL, 'computed', 'response_fraction_per_recording.csv'))
    m = rates.merge(ref[['recording', 'odor', 'frac_trialmean_gt1']], on=['recording', 'odor'], how='left')
    d1 = float((m['rate'] - m['frac_trialmean_gt1']).abs().max())
    e49 = pd.read_csv(H.shared_dir().joinpath(*E49_FRACTIONS))
    e49 = e49[e49['panel'] == 'diag_mix']
    means = rates.groupby(['genotype_group', 'odor'])['rate'].mean()
    d2 = max(abs(float(means.get((g, o), np.nan)) - float(f)) for g, o, f in zip(e49['group'], e49['odor'], e49['frac'])
             if (g, o) in means.index)
    return dict(passed=bool(d1 <= 1e-12 and d2 <= 1e-12), max_abs_diff_per_recording=d1, max_abs_diff_e49_means=d2,
                n=int(len(m)))


def reading_curve(x: pd.DataFrame, reading: str, panel: pd.DataFrame) -> pd.DataFrame:
    """The model columns a Y reading uses: 'dd' (d, or d_settled), 'rr' (ratio, ratio_1spike, or ratio over the
    wiring's panel ratio), 'rate_on' (2+ or 1+ spikes)."""
    y = x.copy()
    y['dd'] = y['d_settled'] if reading == 'settled' else y['d']
    y['rr'] = y['ratio_1spike'] if reading == 'one_spike' else y['ratio']
    y['rate_on'] = y['rate_1spike'] if reading == 'one_spike' else y['rate']
    if reading == 'relative_to_panel':
        y['rr'] = y['rr'] / y['wiring'].map(panel.set_index('wiring')['panel_ratio']).astype(float)
    return y


def y_tables(x: pd.DataFrame, panel: pd.DataFrame, targets: pd.DataFrame, apl: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """B2 and the classes: per reading, odor, estimate and wiring d* with its interval, the flies' APL strength (Y1)
    and the family's own d at m = 1 (Y2) and their classes; the verdicts over wirings, estimates and odors, and held."""
    rows = []
    for rd in Y_READINGS:
        y = reading_curve(x, rd, panel)
        for o in FLIES_WEAK:
            t = targets[(targets['odor'] == o) & (targets['reading'] == rd)].iloc[0]
            a = apl[apl['odor'] == o].iloc[0]
            for s in X2_ESTIMATE_SESSIONS:
                lab = model_label(o, s)
                for w, g in y[y['stimulus'] == lab].groupby('wiring', sort=False):
                    ds = d_star(g['dd'], g['rr'], t['ratio'], t['lo90'], t['hi90'])
                    d1 = float(g.loc[g['m'] == 1.0, 'dd'].iloc[0])
                    rows.append(dict(reading=rd, odor=o, estimate=s, stimulus=lab, wiring=w, flies_ratio=t['ratio'],
                                     flies_lo90=t['lo90'], flies_hi90=t['hi90'], **ds, apl=a['apl'], apl_lo90=a['apl_lo90'],
                                     apl_hi90=a['apl_hi90'], d_m1=d1,
                                     family_ratio_m1=float(g.loc[g['m'] == 1.0, 'rr'].iloc[0]),
                                     y1=y_class((a['apl_lo90'], a['apl_hi90']), (ds['d_star_lo'], ds['d_star_hi'])),
                                     y2=y_class((d1, d1), (ds['d_star_lo'], ds['d_star_hi'])),
                                     max_spike_count=float(g['max_spike_count'].max())))
    t = pd.DataFrame(rows)
    res = {}
    for yk in ('y1', 'y2'):
        by_reading = {}
        for rd, g in t.groupby('reading', sort=False):
            per_odor = {}
            for o, go in g.groupby('odor', sort=False):
                per_est = {s: y_over_wirings(ge[yk]) for s, ge in go.groupby('estimate', sort=False)}
                per_odor[o] = dict(per_estimate=per_est, over_estimates=y_over_estimates(per_est))
            by_reading[rd] = dict(per_odor=per_odor, verdict=y_verdict({o: v['over_estimates'] for o, v in per_odor.items()}))
        res[yk] = dict(by_reading=by_reading, held=held_class({k: v['verdict'] for k, v in by_reading.items()}),
                       carrying=by_reading['carrying']['verdict'])
    return t, res


def b1_table(x: pd.DataFrame, apl: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    """B1 (described): per odor, estimate and wiring the family's response rate with APL on at the flies' APL ratio,
    read linearly and through the indicator family (the ratio to the power 1/n, INDICATOR_NS), interpolated in d (made
    non-increasing), against the control flies' rate."""
    rows = []
    for o in FLIES_WEAK:
        a = float(apl.loc[apl['odor'] == o, 'apl'].iloc[0])
        ctrl = float(targets[(targets['odor'] == o) & (targets['reading'] == 'carrying')]['ctrl'].iloc[0])
        for s in X2_ESTIMATE_SESSIONS:
            for w, g in x[x['stimulus'] == model_label(o, s)].groupby('wiring', sort=False):
                for n in INDICATOR_NS:
                    dd = a ** (1.0 / n)
                    rows.append(dict(odor=o, estimate=s, wiring=w, n=n, d=dd,
                                     rate_on=interp_monotone(g['d'], g['rate'], dd, increasing=False),
                                     rate_on_m1=float(g.loc[g['m'] == 1.0, 'rate'].iloc[0]),
                                     d_max=float(g['d'].max()), ctrl_flies=ctrl))
    return pd.DataFrame(rows)


def b3_table(yt: pd.DataFrame) -> pd.DataFrame:
    """B3 (described): per odor, estimate and wiring the indicator exponent at which the flies' APL ratio meets d*
    (n = ln ratio / ln d*; none where d* is at or above 1, paa's own strength, which a ratio below 1 cannot reach), E98's
    V2 read with it (V2^(1/n)), and Y1 with the flies' APL interval widened by E103's vial SD (beside; waits on
    question 22)."""
    v2 = pd.read_csv(H.shared_dir().joinpath(*E98_A3_RATIOS), comment='#')
    v2 = v2[(v2['ratio_id'] == 'V2') & (v2['reading'] == 'order_matched') & (v2['subset'] == 'all')].iloc[0]
    t = yt[yt['reading'] == 'carrying'].copy()
    ok = np.isfinite(t['d_star']) & (t['d_star'] > 0) & (t['d_star'] < 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        t['n_star'] = np.where(ok, np.log(t['apl']) / np.log(t['d_star']), np.nan)
        t['v2_read'] = float(v2['ratio']) ** (1.0 / t['n_star'])
    t['n_star_note'] = np.where(ok, '', 'no exponent: d* at or above paa\'s (1), which a ratio below 1 reaches with none')
    t['v2'] = float(v2['ratio'])
    wid = [widened(a, m, b) for a, m, b in zip(t['apl_lo90'], t['apl'], t['apl_hi90'])]
    t['apl_lo90_widened'] = [a for a, _ in wid]
    t['apl_hi90_widened'] = [b for _, b in wid]
    t['y1_widened'] = [y_class((a, b), (lo, hi)) for a, b, lo, hi in
                       zip(t['apl_lo90_widened'], t['apl_hi90_widened'], t['d_star_lo'], t['d_star_hi'])]
    return t[['odor', 'estimate', 'wiring', 'apl', 'd_star', 'n_star', 'n_star_note', 'v2', 'v2_read', 'apl_lo90_widened',
              'apl_hi90_widened', 'd_star_lo', 'd_star_hi', 'y1_widened']]


def part_b(wirings: Sequence[str] = ALL_WIRINGS, runs: bool = True) -> dict:
    """Part B: the ladder's runs (`run_ladder`, resumable; `runs` False: tables from the saved runs only), then
    b_model.csv (every run's APL and response rates, d, the ratios), b_flies.csv (the flies' targets), b2_dstar.csv (d*
    and the classes), b1_rates.csv, b3_exponents.csv and b_summary.json."""
    if runs:
        run_ladder(wirings)
    out = out_dir()
    model, panel = ladder_model_table(wirings)
    x = add_d(model)
    rates = weak_kc_rates()
    kid = kc_identity(rates)
    if not kid['passed']:
        raise RuntimeError(f"the KC rates do not reproduce the saved fractions: {kid}")
    targets = flies_targets()
    apl = flies_apl_ratios()
    yt, ys = y_tables(x, panel, targets, apl)
    write_table(x, out / 'b_model.csv')
    write_table(panel, out / 'b_panel.csv')
    write_table(targets, out / 'b_flies.csv')
    write_table(apl, out / 'b_flies_apl.csv')
    write_table(yt, out / 'b2_dstar.csv')
    write_table(b1_table(x, apl, targets), out / 'b1_rates.csv')
    b3 = b3_table(yt)
    write_table(b3, out / 'b3_exponents.csv')
    ident = {w: _read_json(local_dir(LADDER_SUB) / f'{w}_m1_identity.json') for w in wirings if w != 'real'}
    summ = dict(caveat=QUALITY_CAVEAT, unmeasured=UNMEASURED, kc_identity=kid, m1_identities=ident,
                dstar_median=yt[yt['reading'] == 'carrying'].groupby('odor')[['d_star', 'd_star_lo', 'd_star_hi']]
                .median().to_dict('index'),
                y1=ys['y1'], y2=ys['y2'],
                y1_widened=y_verdict({o: y_over_estimates({s: y_over_wirings(g2['y1_widened'])
                                                           for s, g2 in g.groupby('estimate')})
                                      for o, g in b3.groupby('odor')}),
                max_spike_count=float(x['max_spike_count'].max()))
    _write_json(out / 'b_summary.json', summ)
    provenance('B', wirings=list(wirings), runs=runs, threads=os.environ.get('OMP_NUM_THREADS'))
    return summ


# ---------------------------------------------------------------------------
# part D: the counts that outcomes (e) and (f) ask for (no simulation; added with the write-up, 2026-11-17). How many
# flies would decide X1 on each fall, and how many recordings with APL output blocked would bound Y1, if the present
# means and between-fly spreads were the truth (normal approximations, as E98's `flies_needed`)

POWER_Z_B = 0.8416            # one-sided 80 % (E98's POWER_Z)
DECIDE_P = (0.5, 0.8)         # a class reached by the typical draw (P 0.5: the next flies like the present), or 80 % power
COUNT_MAX = 5000              # the most new flies (or recordings) searched
SPLIT_CHECK_MAX_SPLITS = 300_000   # the exact split rule's check (E103's `least_detectable`) while the splits number this
X_ROUTES = ('APL steeper', 'one function', 'KCs steeper')
Y_SCENARIOS = ('tnt_own', 'tnt_ctrl_spread', 'both_own', 'both_ctrl_spread')
LN2 = float(np.log(2.0))


def fall_class_probabilities(delta: float, se, z: float = NORMAL_Z90) -> Dict[str, np.ndarray]:
    """The chance of each class of a fall (`fall_class`) when the estimate of ln Q is normal about `delta` (the true ln
    Q) with SD `se` (scalar or array): the split rule says differ when |estimate| > z se (two-sided 10 %, the normal
    approximation E98's `flies_needed` makes), and Q's 90 % interval is the estimate +- z se. 'APL steeper': the
    estimate below both ln Q_APL_STEEPER and -z se; 'KCs steeper': above both ln Q_KC_STEEPER and z se; 'one function':
    within +- z se with the interval inside Q_ONE_FUNCTION; 'not decided': the rest."""
    from scipy.stats import norm
    d = float(delta)
    s = np.asarray(se, dtype=float)
    if not np.isfinite(d) or np.any(~np.isfinite(s) | (s <= 0)):
        raise ValueError(f'delta must be finite and se positive: {delta}, {se}')
    zs = z * s
    p_apl = norm.cdf((np.minimum(np.log(Q_APL_STEEPER), -zs) - d) / s)
    p_kc = norm.sf((np.maximum(np.log(Q_KC_STEEPER), zs) - d) / s)
    lo = np.maximum(-zs, np.log(Q_ONE_FUNCTION[0]) + zs)
    hi = np.minimum(zs, np.log(Q_ONE_FUNCTION[1]) - zs)
    p_one = np.where(hi > lo, norm.cdf((hi - d) / s) - norm.cdf((lo - d) / s), 0.0)
    return {'APL steeper': p_apl, 'KCs steeper': p_kc, 'one function': p_one,
            'not decided': np.clip(1.0 - p_apl - p_kc - p_one, 0.0, 1.0)}


def limit_class(q: float) -> str:
    """The class of a fall once the split rule detects every difference (flies without bound): 'APL steeper' below
    Q_APL_STEEPER, 'KCs steeper' above Q_KC_STEEPER, 'one function' at Q = 1 exactly, otherwise 'not decided' (Q between
    the cuts: the sets differ, and Q is too near 1 for either steeper class)."""
    if q < Q_APL_STEEPER:
        return 'APL steeper'
    if q > Q_KC_STEEPER:
        return 'KCs steeper'
    return 'one function' if q == 1.0 else 'not decided'


def _reach(m: np.ndarray, p: np.ndarray, level: float) -> Tuple[float, float]:
    """The first and the last count in `m` at which `p` is at least `level` (inf, inf if never; the last is inf while p
    is still at the level at the search's end)."""
    i = np.flatnonzero(p >= level)
    if not len(i):
        return np.inf, np.inf
    return float(m[i[0]]), (np.inf if i[-1] == len(m) - 1 else float(m[i[-1]]))


def flies_to_decide(apl_log: Sequence[float], kc_log: Sequence[float], n_max: int = COUNT_MAX) -> pd.DataFrame:
    """How many flies would decide a fall if the present flies were typical: ln Q's true value the present difference
    of mean (natural) log falls, each set's between-fly SD its present one (`fall_class_probabilities`). Scenarios:
    'kc', more KC flies with the APL set as it is; 'both', the same number of new flies in each set. Per scenario and
    route (X_ROUTES): the chance at the present flies; per level of DECIDE_P the least number of new flies at which the
    route reaches it ('new_flies_p50', with each set's flies in all) and the last ('_last', inf while it still holds at
    `n_max`); the largest chance and where. With the present Q, the SDs (log2) and the class at many flies ('limit')."""
    a = np.asarray(apl_log, dtype=float)
    k = np.asarray(kc_log, dtype=float)
    a, k = a[np.isfinite(a)], k[np.isfinite(k)]
    na, nk = len(a), len(k)
    if na < 2 or nk < 2:
        raise ValueError(f'need two flies or more in each set: {na} APL, {nk} KC')
    delta = float(a.mean() - k.mean())
    sa, sk = float(a.std(ddof=1)), float(k.std(ddof=1))
    if sa == 0 and sk == 0:
        raise ValueError('both sets have no spread')
    meta = dict(n_apl=na, n_kc=nk, q=float(np.exp(delta)), sd_apl_log2=sa / LN2, sd_kc_log2=sk / LN2,
                limit=limit_class(float(np.exp(delta))))
    m = np.arange(0, int(n_max) + 1)
    p_now = fall_class_probabilities(delta, np.sqrt(sa ** 2 / na + sk ** 2 / nk))
    rows = []
    for sc, n_a, n_k in (('kc', np.full(m.shape, na), nk + m), ('both', na + m, nk + m)):
        p = fall_class_probabilities(delta, np.sqrt(sa ** 2 / n_a + sk ** 2 / n_k))
        for c in X_ROUTES:
            row = dict(scenario=sc, cls=c, **meta, p_present=float(p_now[c]))
            for lev in DECIDE_P:
                tag = f'p{int(round(100 * lev))}'
                first, last = _reach(m, p[c], lev)
                row[f'new_flies_{tag}'] = first
                row[f'apl_flies_{tag}'] = float(na + (first if sc == 'both' else 0))
                row[f'kc_flies_{tag}'] = float(nk + first)
                row[f'new_flies_{tag}_last'] = last
            j = int(np.argmax(p[c]))
            row.update(p_max=float(p[c][j]), new_flies_at_max=float(m[j]))
            rows.append(row)
    return pd.DataFrame(rows)


def split_rule_check(apl_log: Sequence[float], kc_log: Sequence[float],
                     max_splits: int = SPLIT_CHECK_MAX_SPLITS) -> dict:
    """The normal approximation against the split rule itself (E103's `least_detectable`): the APL flies as they are,
    n KC flies spread like the present ones (normal scores at their SD) about APL's mean shifted toward the KCs' side.
    The least detectable Q at the present KC count, and the least KC count (searched while the splits number at most
    `max_splits`) at which the present Q is detected: the typical draw's count for 'differ', which both steeper classes
    need (inf if not found in the search)."""
    from experiments import E103_apl_across_odors_heldout as E103
    a = np.asarray(apl_log, dtype=float)
    k = np.asarray(kc_log, dtype=float)
    a, k = a[np.isfinite(a)] / LN2, k[np.isfinite(k)] / LN2
    obs = float(k.mean() - a.mean())             # log2, positive when APL's fall is the steeper
    key = 'ldd_up_log2' if obs >= 0 else 'ldd_down_log2'
    sd = float(k.std(ddof=1))

    def ldd(n):
        return float(E103.least_detectable(a, n_new=n, sd_log2=sd).get(key, np.nan))

    def q_of(v):
        return float(2 ** (-v)) if obs >= 0 else float(2 ** v)
    now = ldd(len(k))
    out = dict(ldd_q_present=q_of(now), n_kc_split_rule=np.inf, n_kc_searched=len(k))
    n = len(k)
    while math.comb(len(a) + n, n) <= max_splits:
        v = ldd(n) if n != len(k) else now
        out['n_kc_searched'] = n
        if np.isfinite(v) and v <= abs(obs):
            out['n_kc_split_rule'] = float(n)
            break
        n += 1
    return out


def x1_counts(vals: pd.DataFrame, n_max: int = COUNT_MAX, split_check: Sequence[str] = ('carrying',)
              ) -> Tuple[pd.DataFrame, dict]:
    """Per reading of `reading_specs` and fall of FALLS, `flies_to_decide` on the per-fly log falls (with
    `split_rule_check` for the readings in `split_check`). The summary: per reading, scenario and level the new flies at
    which X1 would reach each route on X1_MIN falls (the X1_MIN-th smallest over the falls), and held, the most over the
    readings (new flies count alike in every reading: a reading that leaves a KC fly out still gains the same new
    ones)."""
    kc_flies = sorted(vals.loc[vals['set'].isin(KC_SETS_CARRYING), 'fly'].unique())
    specs = reading_specs(kc_flies)
    rows = []
    for name, spec in specs.items():
        for fall in FALLS:
            a = select(vals, spec['apl_sets'], spec, 'apl')
            k = select(vals, spec['kc_sets'], spec, 'kc')
            fa = fly_falls(recording_falls(a, [fall], spec['value']))
            fk = fly_falls(recording_falls(k, [fall], spec['value']))
            t = flies_to_decide(fa.to_numpy(float), fk.to_numpy(float), n_max)
            t.insert(0, 'fall', fall[0])
            t.insert(0, 'reading', name)
            if name in split_check:
                sc = split_rule_check(fa.to_numpy(float), fk.to_numpy(float))
                for c, v in sc.items():
                    t[f'split_rule_{c}'] = v
            rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    summ = {}
    for (rd, sc), g in t.groupby(['reading', 'scenario'], sort=False):
        d = {}
        for c in X_ROUTES:
            x = g[g['cls'] == c]
            d[c] = {f'p{int(round(100 * lev))}': float(np.sort(x[f'new_flies_p{int(round(100 * lev))}'])[X1_MIN - 1])
                    for lev in DECIDE_P}
        summ.setdefault(rd, {})[sc] = d
    held = {sc: {c: {lv: max(summ[rd][sc][c][lv] for rd in summ) for lv in summ['carrying'][sc][c]}
                 for c in X_ROUTES} for sc in ('kc', 'both')}
    return t, dict(by_reading=summ, held=held)


def ratio_spread(rates: pd.DataFrame, odor_label: str, solvent: bool = False, n_boot: int = N_BOOT,
                 seed: int = SEED) -> dict:
    """The parts of the flies' response-rate ratio's spread (`flies_rate_ratio`'s reading, its control draws): the
    ratio; 'se_ctrl', the SD of the log of the control flies' mean rate over `n_boot` bootstrap draws of flies;
    'cv_tnt_own', the blocked recordings' SD over their mean; 'cv_ctrl_flies', the control flies' (SD of fly means over
    their mean), the spread new blocked flies would have if they varied as the control flies do; the counts."""
    x = rates[rates['odor'] == odor_label].set_index('recording')
    if solvent:
        s = rates[rates['odor'] == SOLVENT_ODOR].set_index('recording')['rate']
        x = x.assign(rate=x['rate'] - s.reindex(x.index).to_numpy(float))
    c = x[x['genotype_group'] == 'control']
    b = x[x['genotype_group'] == 'apl_output_blocked']
    flies = sorted(c['fly'].unique())
    by_fly = [c.loc[c['fly'] == f, 'rate'].to_numpy(float) for f in flies]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(flies), size=(n_boot, len(flies)))
    sums = np.array([v.sum() for v in by_fly])
    ns = np.array([len(v) for v in by_fly])
    ctrl_b = sums[idx].sum(axis=1) / ns[idx].sum(axis=1)
    fm = np.array([v.mean() for v in by_fly])
    ctrl, tnt = float(c['rate'].mean()), float(b['rate'].mean())
    return dict(stimulus=odor_label, ratio=tnt / ctrl, se_ctrl=float(np.std(np.log(ctrl_b))),
                cv_tnt_own=float(b['rate'].std(ddof=1) / tnt), cv_ctrl_flies=float(fm.std(ddof=1) / fm.mean()),
                n_ctrl_flies=len(flies), n_tnt_recordings=int(len(b)), n_tnt_flies=int(b['fly'].nunique()))


def y_wiring_bounds(y: pd.DataFrame, odor: str, strength, labels: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Per estimate and wiring the family's response-rate ratio (`reading_curve`'s 'rr' against 'dd', made
    non-decreasing) at a strength's two ends: 'a' at its top, 'b' at its bottom, each clipped to the ladder's range.
    `strength(wiring, g)` gives (lo, hi); `labels` maps each estimate to the odor's stimulus label (default
    `model_label`). A flies' ratio interval [lo, hi] makes a wiring 'too weak' when lo > a (then d*'s lower end lies
    above the strength), 'too strong' when hi < b, as `y_class` on `d_star` does."""
    labels = {s: model_label(odor, s) for s in X2_ESTIMATE_SESSIONS} if labels is None else labels
    rows = []
    for s, lab in labels.items():
        for w, g in y[y['stimulus'] == lab].groupby('wiring', sort=False):
            lo, hi = strength(w, g)
            dd = g['dd'].to_numpy(float)
            rr = g['rr'].to_numpy(float)
            ok = np.isfinite(dd) & np.isfinite(rr)

            def f(v):
                return interp_monotone(dd, rr, float(np.clip(v, dd[ok].min(), dd[ok].max())), increasing=True)
            rows.append(dict(estimate=s, wiring=w, s_lo=float(lo), s_hi=float(hi), a=f(hi), b=f(lo)))
    return pd.DataFrame(rows)


def y_odor_class(lo: float, hi: float, bounds: pd.DataFrame) -> str:
    """An odor's class from a flies' ratio interval [lo, hi] and `y_wiring_bounds`: per wiring too weak, too strong or
    inside; over the wirings (`y_over_wirings`) and the estimates (`y_over_estimates`)."""
    per_est = {}
    for s, g in bounds.groupby('estimate', sort=False):
        cls = np.where(lo > g['a'], 'too weak', np.where(hi < g['b'], 'too strong', 'inside'))
        per_est[s] = y_over_wirings(list(cls))
    return y_over_estimates(per_est)


def y_thresholds(bounds: pd.DataFrame, min_wirings: int = X2_MIN_WIRINGS) -> Tuple[float, float]:
    """The flies' ratio above which an odor is too weak under every estimate (A: the most over the estimates of the
    `min_wirings`-th smallest 'a') and below which it is too strong under every estimate (B: the least of the
    (n - `min_wirings` + 1)-th smallest 'b'); `y_odor_class` gives 'too weak' iff lo > A and 'too strong' iff hi < B."""
    A, B = -np.inf, np.inf
    for _, g in bounds.groupby('estimate', sort=False):
        a, b = np.sort(g['a'].to_numpy(float)), np.sort(g['b'].to_numpy(float))
        A = max(A, float(a[min_wirings - 1]))
        B = min(B, float(b[len(b) - min_wirings]))
    return A, B


def _se_y(sp: dict, scenario: str, m: np.ndarray) -> np.ndarray:
    """The SD of the log of the flies' ratio with `m` new blocked recordings ('tnt_*': the control flies as they are)
    or `m` new flies of each genotype ('both_*'), new blocked flies spread as the present three ('_own') or as the
    control flies ('_ctrl_spread')."""
    cv = sp['cv_tnt_own'] if scenario.endswith('_own') else sp['cv_ctrl_flies']
    nb = sp['n_tnt_recordings'] + m
    if scenario.startswith('tnt'):
        sc2 = np.full(m.shape, sp['se_ctrl'] ** 2)
    else:
        sc2 = sp['se_ctrl'] ** 2 * sp['n_ctrl_flies'] / (sp['n_ctrl_flies'] + m)
    return np.sqrt(sc2 + cv ** 2 / nb)


def y_counts(x: pd.DataFrame, panel: pd.DataFrame, targets: pd.DataFrame, apl: pd.DataFrame,
             rates_all: Optional[pd.DataFrame] = None, rates_first: Optional[pd.DataFrame] = None,
             n_max: int = COUNT_MAX) -> Tuple[pd.DataFrame, dict]:
    """How many blocked recordings would bound Y1 (and Y2) if the flies' present ratio were the truth: per Y reading,
    verdict and odor, the thresholds A and B (`y_thresholds`), the class now (the saved interval, `y_odor_class`; the
    identity with part B is checked) and at the limit (lo = hi = the ratio), and per scenario (Y_SCENARIOS, `_se_y`) and
    class (too weak, too strong) the chance with nothing added and the least m at which it reaches each of DECIDE_P (the
    ratio's estimate lognormal about the present one with SD `_se_y`, its 90 % interval +- NORMAL_Z90 SD). The summary:
    per reading, verdict, scenario and level the m at which Y_MIN_ODORS odors are too weak (the Y_MIN_ODORS-th smallest
    over the odors), and held, the most over the readings."""
    from scipy.stats import norm
    rates_all = weak_kc_rates() if rates_all is None else rates_all
    rates_first = weak_kc_rates(FIRST_TRIALS) if rates_first is None else rates_first
    saved = _read_json(out_dir() / 'b_summary.json')
    m = np.arange(0, int(n_max) + 1)
    rows = []
    for rd in Y_READINGS:
        y = reading_curve(x, rd, panel)
        for o in FLIES_WEAK:
            lab = flies_stim_label(o)
            t = targets[(targets['odor'] == o) & (targets['reading'] == rd)].iloc[0]
            sp = ratio_spread(rates_first if rd == 'first_three_trials' else rates_all, lab,
                              solvent=rd == 'solvent_subtracted')
            r = float(t['ratio'])
            ap = apl[apl['odor'] == o].iloc[0]
            for yk, strength in (('y1', lambda w, g, ap=ap: (ap['apl_lo90'], ap['apl_hi90'])),
                                 ('y2', lambda w, g: (float(g.loc[g['m'] == 1.0, 'dd'].iloc[0]),) * 2)):
                bd = y_wiring_bounds(y, o, strength)
                A, B = y_thresholds(bd)
                now = y_odor_class(float(t['lo90']), float(t['hi90']), bd)
                part_b = (saved or {}).get(yk, {}).get('by_reading', {}).get(rd, {}).get('per_odor', {}).get(o, {})
                base = dict(reading=rd, verdict=yk.upper(), odor=o, flies_ratio=r, flies_lo90=float(t['lo90']),
                            flies_hi90=float(t['hi90']), A_too_weak_above=A, B_too_strong_below=B, class_now=now,
                            class_part_b=part_b.get('over_estimates', ''), class_limit=y_odor_class(r, r, bd),
                            se_ctrl=sp['se_ctrl'], cv_tnt_own=sp['cv_tnt_own'], cv_ctrl_flies=sp['cv_ctrl_flies'],
                            n_ctrl_flies=sp['n_ctrl_flies'], n_tnt_recordings=sp['n_tnt_recordings'])
                for scn in Y_SCENARIOS:
                    se = _se_y(sp, scn, m)
                    for cls, gap in (('too weak', np.log(r) - np.log(A)), ('too strong', np.log(B) - np.log(r))):
                        p = norm.cdf((gap - NORMAL_Z90 * se) / se) if np.isfinite(gap) else np.zeros(m.shape)
                        row = dict(base, scenario=scn, cls=cls, log_gap=float(gap), se_present=float(se[0]),
                                   p_present=float(p[0]))
                        for lev in DECIDE_P:
                            row[f'new_p{int(round(100 * lev))}'] = _reach(m, p, lev)[0]
                        row['se_floor'] = sp['se_ctrl'] if scn.startswith('tnt') else 0.0
                        rows.append(row)
    t = pd.DataFrame(rows)
    ident = bool((t['class_now'] == t['class_part_b']).all()) if saved else None
    summ = {}
    for (rd, yk, scn), g in t[t['cls'] == 'too weak'].groupby(['reading', 'verdict', 'scenario'], sort=False):
        summ.setdefault(yk, {}).setdefault(rd, {})[scn] = {
            f'p{int(round(100 * lev))}': float(np.sort(g[f'new_p{int(round(100 * lev))}'])[Y_MIN_ODORS - 1])
            for lev in DECIDE_P}
    held = {yk: {scn: {lv: max(summ[yk][rd][scn][lv] for rd in summ[yk]) for lv in summ[yk]['carrying'][scn]}
                 for scn in Y_SCENARIOS} for yk in summ}
    limit = {yk: {rd: y_verdict(dict(zip(g['odor'], g['class_limit'])))
                  for rd, g in t[(t['verdict'] == yk) & (t['scenario'] == Y_SCENARIOS[0]) & (t['cls'] == 'too weak')]
                  .groupby('reading', sort=False)} for yk in ('Y1', 'Y2')}
    return t, dict(identity_with_part_b=ident, too_weak_by_reading=summ, too_weak_held=held, limit=limit)


def part_d() -> dict:
    """Part D (no simulation): the counts that outcomes (e) and (f) ask for, from part A's and part B's saved tables:
    d_x1_counts.csv (`x1_counts`), d_y_counts.csv (`y_counts`) and d_summary.json."""
    out = out_dir()
    vals = read_table(out / 'a_per_stimulus.csv')
    x1t, x1s = x1_counts(vals)
    write_table(x1t, out / 'd_x1_counts.csv')
    x = read_table(out / 'b_model.csv')
    panel = read_table(out / 'b_panel.csv')
    targets = read_table(out / 'b_flies.csv')
    apl = read_table(out / 'b_flies_apl.csv')
    yt, ys = y_counts(x, panel, targets, apl)
    write_table(yt, out / 'd_y_counts.csv')
    summ = dict(caveat=QUALITY_CAVEAT, decide_p=list(DECIDE_P), count_max=COUNT_MAX,
                note=('normal approximations with the present means and between-fly spreads as the truth; new flies '
                      '(X1: per set; Y: blocked recordings, or flies of each genotype) beyond the present ones; inf: '
                      'not reached within count_max'),
                x1=x1s, y=ys)
    _write_json(out / 'd_summary.json', summ)
    provenance('D')
    return summ


VARIANT_PARAM_COLS = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'kc_kc_inh',
                      'kc_kc_tau', 'pn_inhadd', 'pn_inhsc', 'k')


def write_variants() -> pd.DataFrame:
    """results/E104_kc_calyx_vs_apl/variants.csv, which the archive builds params.csv from: one row per saved ladder
    run (tag = the files' variant, '<wiring>_m<m>'), with m, the absolute APL scales of the run and the fit's, every
    scalar parameter of the run (base_params, JSON) and its columns of VARIANT_PARAM_COLS, the stimuli and the largest
    spike count; m = 1 and m = 0 (acute off) are E98's and E101's saved runs, archived with those experiments."""
    runs = pd.read_csv(local_dir(LADDER_SUB) / 'runs.csv')
    runs = runs.drop_duplicates(['wiring', 'm'], keep='last')
    rows = []
    for r in runs.to_dict('records'):
        p = json.loads(r['params_json'])
        tag = f"{r['wiring']}_m{float(r['m']):g}"
        row = dict(tag=tag, variant=tag, part='B', wiring=r['wiring'], m=float(r['m']),
                   fit_wAPLKC_scale=r['fit_wAPLKC_scale'], fit_wAPLPN_scale=r['fit_wAPLPN_scale'],
                   base_params=r['params_json'], stimuli=r['stimuli'], max_count_on=r['max_spike_count'],
                   ref_experiment='E72_wiring_population (part B fits; E83 wiring labels)', host=r['host'],
                   queue_commit=r.get('queue_commit'), pyolfsysm_commit=r.get('pyolfsysm_commit'))
        for c in VARIANT_PARAM_COLS:
            if c in p:
                try:
                    row[c] = float(p[c])
                except (TypeError, ValueError):
                    row[c] = p[c]
        rows.append(row)
    v = pd.DataFrame(rows)
    v.to_csv(out_dir() / 'variants.csv', index=False)
    return v


def part_c() -> dict:
    """Part C: the verdicts X1 (with G1 to G3), X2, X3, Y1 and Y2 in one table (c_verdicts.csv: per verdict and
    reading, its class; the held class), and verdict.json."""
    out = out_dir()
    a = _read_json(out / 'a_summary.json')
    b = _read_json(out / 'b_summary.json')
    if a is None or b is None:
        raise RuntimeError('part C needs parts A and B')
    x1 = read_table(out / 'a2_x1.csv')
    rows = []
    for _, r in x1.iterrows():
        for k in [f for f, _, _ in FALLS] + ['x1']:
            rows.append(dict(verdict='X1' if k == 'x1' else f'X1: {k}', reading=r['reading'], cls=r[k]))
    for k in ('x2', 'x3'):
        for rd, v in a[k]['by_reading'].items():
            rows.append(dict(verdict=k.upper(), reading=rd, cls=v))
    for k in ('y1', 'y2'):
        for rd, v in b[k]['by_reading'].items():
            rows.append(dict(verdict=k.upper(), reading=rd, cls=v['verdict'],
                             per_odor=';'.join(f"{o}: {p['over_estimates']}" for o, p in v['per_odor'].items())))
    held = {'X1': a['x1']['x1'], **{f'X1: {f}': a['x1']['held_falls'][f] for f, _, _ in FALLS}, 'X2': a['x2']['x2'],
            'X3': a['x3']['x3'], 'Y1': b['y1']['held'], 'Y2': b['y2']['held']}
    for k, v in held.items():
        rows.append(dict(verdict=k, reading='held', cls=v))
    t = pd.DataFrame(rows)
    write_table(t, out / 'c_verdicts.csv')
    v = dict(caveat=QUALITY_CAVEAT, held=held, y1_widened_beside=b.get('y1_widened'),
             lock_sha256=a.get('lock_sha256'), time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(out / 'verdict.json', v)
    provenance('C')
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'variants'))
    ap.add_argument('--runs', action='store_true', help="part 0: also the runner's identities on the real wiring")
    ap.add_argument('--write-lock', action='store_true',
                    help="part 0: write part A's lock once every identity passes (commit it before part A)")
    ap.add_argument('--wirings', default=','.join(ALL_WIRINGS), help='part B: comma-separated wirings')
    ap.add_argument('--no-runs', action='store_true', help="part B: the tables from the saved runs only")
    ap.add_argument('--threads', type=int, default=None,
                    help='runs (part 0 with --runs, part B): OMP threads; OMP_NUM_THREADS and OMP_THREAD_LIMIT must match')
    a = ap.parse_args()
    simulates = (a.part == '0' and a.runs) or (a.part == 'B' and not a.no_runs)
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
        print(json.dumps(part_b(tuple(x for x in a.wirings.split(',') if x), runs=not a.no_runs), indent=1,
                         default=_json_default), flush=True)
    elif a.part == 'D':
        print(json.dumps(part_d(), indent=1, default=_json_default), flush=True)
    elif a.part == 'variants':
        v = write_variants()
        print(f'variants.csv: {len(v)} rows', flush=True)
    else:
        print(json.dumps(part_c(), indent=1, default=_json_default), flush=True)


if __name__ == '__main__':
    main()
