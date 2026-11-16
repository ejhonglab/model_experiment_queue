#!/usr/bin/env python3
"""E101: the flies' APL against their KC claws in the same fly and plane. Yang's dual-colour recordings (GCaMP6f in APL,
jRGECO1a in the KCs; pa, paa and ms as concentration series), both channels, read for whether APL's response is one function
of the claws' summed response (X1), how the claws and APL change across the three odors and along concentration against the
best family's claw stage and APL (X2), and whether APL's calcium follows the claws' activity place by place (X3). A
diagnostic: it moves no loss.

The entry is item E101 of experiment_priorities.md (proposed by hal 2026-11-15 16:01; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-15 15:32, at Tom's instruction in hal's terminal, 15:30, that both channels of
the dual-colour recordings be used). C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each of
E72 part B's eight wirings at its own fit (E83's `saved_fit` / `build_simulator`), nothing re-tuned. Two rules of
2026-11-15 from viki's terminal hold throughout: only recordings with Yang's processed products are read (15:28;
recordings.csv's `use_until_further_notice`), and no connectome table made by Yang's Claude is read (15:35; none is needed).

Reused, not duplicated:
- E100 (experiments/E100_apl_weak_odors.py): the readers of viki's copy (data/yang_apl_recordings/ in the shared directory;
  `read_recordings`, `yaml_stimuli`, `experiment_xml`, `usable_row`, `released_2021`), E79's outline reader, and the two
  refusals (a recording Yang has not released; a raw-only row).
- E98 (experiments/E98_apl_response_strength.py): the indicator family (`indicator_ratio`, `model_ratio_set`), `classify`,
  `comparison_verdict`, `hull`, `v1_verdict`, the extended input of B1 (`b1_input`, `b1_specs`) and its runner
  (`ExtendedRun`), B1's saved runs (`load_b1_run`) for the identity of part A3.
- E90 (experiments/E90_model_stages.py): the stage reading (`read_stages`, MODEL_WINDOW, PRE_WINDOW), the claw versions
  (`claw_signal`, CLAW_SPIKE_SHARES) and the claws' glomeruli (`glomerulus_of_claws`).

The groups of two-channel recordings (the entry's table, fixed here; GROUPS):
- L, Yang's own list (notes/yang_apl_calyx_metadata.md, answer 7): 20 of the 21 she named (20221016_4FLfun5 is not among
  the copied recordings), 7 flies, 7 calyces; the verdicts.
- L+ (9 recordings of 5 of the same flies, not in her list), V0 (the development flies of 2022-07-15 to 09-08, 23
  recordings, 5 flies, 8 calyces), VT (the VTUG flies, 33 recordings, 5 flies): beside, decide nothing.
- G (green only, both channels recorded: 2022-10-05 fly 1, 2022-12-17 fly 6; 6 recordings): the leak of green into red.
- R (red only: 2022-07-15 fly 4's right calyx, 1 recording; 2022-10-05 fly 2, 12 recordings): the reverse leak.
- M (MB-DsRed with GCaMP6f in APL: 2022-07-05 fly 1; 1 recording): movement.
- Left out: 20220719_2FRfun1 (a plane at 87 um, one channel processed). Out under the raw-only rule, listed in A1 and not
  read: the fall 2024 sets (F1, F2), R2023 and the unprocessed runs of sets D, V and DEV.

Definitions (the entry's, fixed before any reading):
- The channels: `_C1` is APL (GCaMP6f, green PMT) and `_C2` the KCs (jRGECO1a, red PMT), as viki's copy infers; nothing is
  read for results before Yang confirms it (question 2; CHANNELS_CONFIRMED).
- A recording's maps: per presentation and channel, Yang's peak dF/F map in the copied variant (APL's rigid
  motion-corrected; the KC channel moved by the same shifts), inside her calyx outline, in cells of CELL_UM (4 um) made with
  the recording's own pixel size (CELL_UM_BESIDE, 2 and 8 um, beside); a cell enters when at least CELL_MIN_INSIDE of its
  pixels are inside the outline.
- A calyx response: the mean of a map inside the outline. A stimulus's response: the mean over its presentations minus the
  mean of the recording's first-in-block solvent presentations (the carrying reading, 'carrying'); partners: minus its own
  block's solvent ('block'), and the first two presentations only ('first_two'). READINGS.
- The KC channel with the leak removed: K_b = C2 - b C1, presentation by presentation, for maps and for calyx responses.
  b_G is the median, over cells and over group G's recordings, of the least-squares slope of the KC channel on the APL
  channel across stimuli. Every KC quantity is read at b = 0 and b_G, with 2 b_G beside (B_MULTIPLES).
- Kept stimuli of a recording: those whose APL response and KC response (at b_G) each exceed KEPT_SD_MULT (2) times the SD
  of that channel's solvent responses in the recording (every solvent presentation's calyx response).
- A fly: its recordings averaged (ratios on the log scale). Flies are the units; intervals are 90 %, from N_BOOT bootstrap
  draws of flies.
- The model's claws: E90's claw stage (the wiring times the bouton rate minus APL's inhibition, clipped at 0; window mean
  minus pre-odor mean), summed over the claws of the 36 glomeruli with input ('input'; all claws beside, 'all'), with APL on,
  in E90's three versions of what a claw's calcium reports (CLAW_SHARES: the drive alone, and with 0.25 and 0.5 of the
  KC's spike signal; scales from the wiring's APL-on run on the 16 panel odors, as E90). The model's APL: as in E98 (window
  mean over 0-0.5 s and pre-odor mean).
- The indicators: APL through E98's family (APL_NS, n 1 to 2.5, both forms); the claws through CLAW_NS (n 1 to 2 in steps of
  0.25), both forms, applied to the claw total's window and pre-odor means as E98 applies APL's (jRGECO1a's Hill
  coefficient is about 1.9: Dana et al. 2016, eLife 5:e12727).
- Held. A class is held if b = 0 and b = b_G, the carrying reading and both partners give it, and the flies the planner had
  not probed give it by the same rule: for X1 and X2 the five flies of 2022-12 (UNPROBED_X12; at least HELD_X12_MIN, 4, of
  the five must remain), for X3 the three flies whose maps were not opened (UNPROBED_X3; HELD_X3_MIN, 2, of the three must
  agree). A class that is not held is reading-dependent and counts as not decided (`held_class`).

The verdicts (fixed here, on set L with flies as units; `x1_verdict`, `x2_quantity_verdict`, `x3_verdict`):
- X1. For each of the three pairs of odors (X1_PAIRS), the flies' mean d (the difference of the two odors' mean residuals
  from the reduced-major-axis line of log2 APL on log2 K_b over the kept stimuli) with its interval. One function if all
  three intervals lie within +-X1_BOUND (0.585, a factor of 1.5); depends on the odor if some pair's interval lies wholly
  beyond it; otherwise not decided; decided only with X1_MIN_FLIES (5) or more flies in every pair. A recording enters a
  pair only with X1_MIN_KEPT_PER_ODOR (2) or more kept stimuli of each odor.
- X2. Per ratio (X2_RATIOS: r1 pa over ms near 10^-3, r2 paa over ms near 10^-3, r3 paa at 10^-5 over 10^-4) and quantity
  (the claws, APL, Q = APL's ratio over the claws' ratio within the fly), E98's class over the eight wirings (6 of 8). A
  quantity is flatter, inside or steeper if at least X2_MIN_RATIOS (2) of the three ratios are and the third is not outside
  on the other side (for 'inside': the third is anything; for 'flatter': not 'steeper'; for 'steeper': not 'flatter');
  otherwise not decided. r1 and r2 under each estimate of ms's input: inside if inside under any, outside only if outside
  under all (E98's V1). The model's Q is read over every combination of the two indicator families, their forms and the
  three claw versions; a wiring is inside if any combination lies in the flies' interval.
- X3. A fly passes at b if its A(b) (mean over its recordings) is at least X3_PASS_A (0.3) and above the mean of its
  recordings' null percentiles; it fails if A(b) is below X3_FAIL_A (0.15) or below that percentile while both channels'
  repeat r is at least X3_FAIL_MIN_REPEAT_R (0.5). "Follows" if B1's cross-check passed and at least X3_MIN_FLIES (5) of
  the 7 flies pass at b = 0, b_G and 2 b_G, X3_UNPROBED_MIN (2) of the three unopened flies among them ("at b = 0 and b_G"
  alone when A4's rule says so, a4_rule.json); "follows at the controls' leak" if that holds at 0 and b_G but not at 2 b_G
  (counts as not decided until baseline images settle b); "does not follow" if at least 5 of the 7 flies fail at b = 0;
  otherwise not decided; with fewer than 5 flies holding a recording with X3_MIN_KEPT (4) kept stimuli from X3_MIN_ODORS (2)
  odors, not decided.

Parts (CLI `--part`):
- 0, the identities; a failure stops the experiment and is reported (identity.json). Those that need no new recording:
  the reader gives E100's calyx means on one released 2021 recording; the stage run on the extended input reproduces E98's
  saved B1 counts with 0 differences and its APL means (where E98 has them; checked in A3), and E90's saved claw-stage sums
  on the 16 panel odors (`e90_claw_identity`; reads E90's saved stage matrices); a blank evokes nothing at the model's claws
  or in its APL. After Yang's release: each channel's calyx means against Yang's own table of that channel; each
  recording's frames equal its presentations; both channels' maps have one shape and one outline; pixel size, channels and
  gains from each Experiment.xml.
- A (before any new recording is read): A1 the table of the two-channel recordings by group (a1_recordings.csv,
  a1_groups.csv; recordings.csv, the stimulus files and the date notes only, no map values); A2 the family's premise from
  E90's saved stage matrices (no simulation): evoked APL against the summed claw stage over the 16 panel odors on the eight
  wirings, the line of log APL on log claws with each odor's residual, and pa over ms, paa over ms and paa over pa for the
  claws, APL, the PN totals and the KCs' response rates (a2_family.csv, a2_lines.csv, a2_ratios.csv); A3 the fixed-parameter
  runs of E90's stage reading on E98's extended input (B1's stimuli), eight wirings, APL on and acutely off, with the claw
  stage saved (outputs/E101_apl_vs_kc_claws/A3/; a3_runs.csv, a3_stimuli.csv, a3_family.csv, a3_ratios.csv,
  a3_identity.json); A4 what the place test can detect on made-up two-channel maps (a4_power.csv, a4_rule.json).
- B (after Yang's release): B1 group G's slope per cell, b_G per fly and pooled, and the cross-check between the two
  green-only flies; B2 group R's APL channel where APL carries no indicator; B3 group M's red marker; B4 set L's solvent
  maps; B5 set L's own slope against b_G (below twice b_G: outcome (h)).
- C (after Yang's release): C1 per recording and stimulus APL's response and K_b; C2 X1's statistic; C3 X2's three ratios
  for the claws, APL and Q; C4 beside (V0, L+, VT, Hallem's input, APL off, Yang's single-colour group tables).
- D (after Yang's release): X3's statistic per recording of set L, with its null, ceiling, cross-correlogram and E79's
  common-map share.
- E (after Yang's release; described): the claws' own odor-odor distances from K_b's maps.
- score: the verdicts (verdict.json).

Deviations from the entry (each also in the note):
- A3's stimuli are B1's as `b1_input` gives them (46 of `b1_specs`' 52: Hallem's six rows below 10^-2 have no source and are
  left out, as in E98's first B1 runs, and listed in a3_stimuli.csv).
- A3's check against E90's saved claw-stage sums and A2 read E90's saved stage matrices; the code subagent that wrote this
  could not read them (the permission system refused the read, as it had the planning subagent's), so `--part A2` and
  `--part 0 --e90` are left for hal's session. A3 reproduces the saved fit's counts on the 16 panel odors (E83's
  `saved_fit`) and E98's B1 where saved.
- The blank's evoked APL is read against the settled pre-odor level (E98's SETTLED_PRE_WINDOW, -0.25 to 0 s) as well as
  E90's pre-odor window, which holds the KC layer's start-up (E98's part 0 found 0.51 with E90's window on the real wiring).
- Maps in part D are the stimulus maps minus the recording's solvent map (the mean of its first-in-block solvent maps), as
  the carrying reading subtracts the solvent's calyx response; the partners subtract the block's solvent map, or use the
  first two presentations (whose halves are then the 1st and the 2nd).
- A4's made-up flies: 3 recordings each, 10 kept stimuli of 3 odors (4, 3, 3) on a 15 x 15 grid of cells (A4_*); the
  planted leak is A4_B_G (0.034, the planner's green-only slopes 0.028 and 0.039, the entry's Probes) with the KC channel's
  odor-evoked amplitude A4_KC_OVER_APL (0.2) of APL's; "the made-up flies" of A4's rule are pooled over the whole grid of
  reliabilities (0.3 to 0.9 per channel, every pair equally weighted), the literal reading of the entry.

CLI (from the repository root):
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A1
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A2              # reads E90's saved stage matrices
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E101_apl_vs_kc_claws --part A3
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A4
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part B           # after Yang's release (YANG_RELEASED)
    --dry-run lists what a part would read and write, and exits.

Outputs: results/E101_apl_vs_kc_claws/ in the shared directory (the tables above, identity.json, verdict.json,
provenance_<part>.json); A3's runs under outputs/E101_apl_vs_kc_claws/A3/ (to be archived to
data/model_outputs/E101_apl_vs_kc_claws/).
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

EXPERIMENT = 'E101_apl_vs_kc_claws'
ALL_WIRINGS = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E83's / E98's

# ---------------------------------------------------------------------------
# the groups (the entry's table; fixed before any reading)

GROUP_L: Tuple[str, ...] = (
    '20221016_3FLfun2', '20221016_4FRfun2', '20221016_4FRfun3', '20221016_4FRfun4',
    '20221207_1FLfun2', '20221207_1FLfun4', '20221207_1FLfun5', '20221207_2FRfun1', '20221207_2FRfun2',
    '20221208_1FRfun1', '20221208_1FRfun2',
    '20221213_2FLfun2', '20221213_2FLfun4', '20221213_2FLfun5',
    '20221217_2FLfun3', '20221217_2FLfun4', '20221217_2FLfun5', '20221217_2FLfun6', '20221217_2FLfun7', '20221217_2FLfun8')
YANG_LIST_NOT_COPIED = ('20221016_4FLfun5',)       # named by Yang (answer 7), not among the copied recordings
GROUP_L_PLUS: Tuple[str, ...] = (
    '20221016_3FLfun3', '20221016_3FLfun4', '20221016_3FLfun5', '20221016_4FRfun1', '20221016_4FRfun5',
    '20221207_1FLfun3', '20221207_2FRfun3', '20221213_2FLfun1', '20221213_2FLfun6')
GROUP_V0: Tuple[str, ...] = (
    '20220715_4FLfun1', '20220715_4FLfun2', '20220715_4FLfun4',
    '20220719_4FLfun1', '20220719_4FLfun2',
    '20220818_4FLfun1', '20220818_4FLfun2', '20220818_4FLfun3', '20220818_4FRfun1', '20220818_4FRfun2', '20220818_4FRfun3',
    '20220906_1FLfun1', '20220906_1FLfun2', '20220906_1FLfun3', '20220906_1FRfun1',
    '20220908_3FLfun1', '20220908_3FLfun2', '20220908_3FLfun3', '20220908_3FLfun4',
    '20220908_3FRfun1', '20220908_3FRfun2', '20220908_3FRfun3', '20220908_3FRfun4')
GROUP_VT: Tuple[str, ...] = (
    '20220826_2FLfun1', '20220826_2FLfun2', '20220826_2FLfun3', '20220826_2FLfun4',
    '20220826_2FRfun1', '20220826_2FRfun2', '20220826_2FRfun3',
    '20220915_1FLfun1', '20220915_1FLfun2', '20220915_1FLfun3', '20220915_1FLfun4',
    '20220915_1FRfun1', '20220915_1FRfun2', '20220915_1FRfun3',
    '20221006_1FLfun2', '20221006_1FLfun3', '20221006_1FLfun4', '20221006_1FLfun5', '20221006_1FRfun1', '20221006_1FRfun2',
    '20230126_1FLfun2', '20230126_1FLfun3', '20230126_1FLfun4', '20230126_1FLfun5',
    '20230210_1FLfun2', '20230210_1FLfun3', '20230210_1FLfun4', '20230210_1FLfun5', '20230210_1FLfun6', '20230210_1FLfun7',
    '20230210_1FRfun1', '20230210_1FRfun2', '20230210_1FRfun3')
GROUP_G: Tuple[str, ...] = ('20221005_1MRfun2', '20221005_1MRfun3', '20221005_1MRfun4',
                            '20221217_6FRfun2', '20221217_6FRfun3', '20221217_6FRfun5')
GROUP_R: Tuple[str, ...] = ('20220715_4FRfun2',
                            '20221005_2FLfun1', '20221005_2FLfun2', '20221005_2FLfun3', '20221005_2FLfun4',
                            '20221005_2FRfun1', '20221005_2FRfun2', '20221005_2FRfun3', '20221005_2FRfun4',
                            '20221005_2FRfun5', '20221005_2FRfun6', '20221005_2FRfun7', '20221005_2FRfun8')
GROUP_M: Tuple[str, ...] = ('20220705_1FRfun1',)
LEFT_OUT: Tuple[str, ...] = ('20220719_2FRfun1',)  # a plane at 87 um, one channel processed
GROUPS: Dict[str, Tuple[str, ...]] = dict(L=GROUP_L, **{'L+': GROUP_L_PLUS}, V0=GROUP_V0, VT=GROUP_VT, G=GROUP_G, R=GROUP_R,
                                          M=GROUP_M)
GROUP_USE = {'L': 'the verdicts', 'L+': 'beside', 'V0': 'the repetition; decides nothing',
             'VT': 'described apart, once Yang has said what VTUG labels', 'G': "the leak's size", 'R': 'the reverse leak',
             'M': 'movement'}
# the sets whose two-channel recordings are raw only (listed in A1, not read)
RAW_ONLY_SETS = ('F1', 'F2', 'R2023')
# the flies the planner did not probe: for X1 and X2 the five flies of 2022-12 (it read both channels' tables only for the
# two flies of 2022-10-16), for X3 the three flies whose maps were not opened
UNPROBED_X12 = ('20221207_1', '20221207_2', '20221208_1', '20221213_2', '20221217_2')
UNPROBED_X3 = ('20221016_4', '20221207_2', '20221208_1')
HELD_X12_MIN = 4                         # of the five
HELD_X3_MIN = 2                          # of the three

# ---------------------------------------------------------------------------
# the definitions

CHANNELS = ('C1', 'C2')                  # C1 APL (GCaMP6f, green PMT), C2 the KCs (jRGECO1a, red PMT), as viki's copy infers
CHANNEL_ROLE = {'C1': 'APL', 'C2': 'KC'}
CELL_UM = 4.0
CELL_UM_BESIDE = (2.0, 8.0)
CELL_MIN_INSIDE = 0.5                    # a cell enters when at least this share of its pixels is inside the outline
READINGS = ('carrying', 'block', 'first_two')
CARRYING = 'carrying'
B_MULTIPLES = (0.0, 1.0, 2.0)            # b = 0, b_G, 2 b_G
KEPT_SD_MULT = 2.0
N_BOOT = 10000
SEED = 0
CI = (5.0, 95.0)                         # 90 % intervals
MIN_RESPONSE = 0.02                      # E98's: a ratio's denominator must exceed this (dF/F)
MATCH_LOG_UNITS = 1.0                    # Tom's rule, inclusive
CLAW_SHARES = (0.0, 0.25, 0.5)           # E90's CLAW_SPIKE_SHARES
CLAW_SUBSETS = ('input', 'all')          # the claws of the 36 glomeruli with input; all claws (beside)
APL_NS = tuple(float(x) for x in np.arange(1.0, 2.5 + 1e-9, 0.25))     # E98's INDICATOR_NS
CLAW_NS = tuple(float(x) for x in np.arange(1.0, 2.0 + 1e-9, 0.25))
FORMS = ('baseline', 'evoked')
VERDICT_MIN_WIRINGS = 6                  # of the 8 wirings (E98's)

# X1
X1_PAIRS: Tuple[Tuple[str, str], ...] = (('pa', 'ms'), ('paa', 'ms'), ('paa', 'pa'))
X1_BOUND = float(np.log2(1.5))           # 0.585
X1_MIN_FLIES = 5
X1_MIN_KEPT_PER_ODOR = 2

# X2: (id, numerator readings, denominator, note); a stimulus is (odor, log10 concentration); the flies' interval is the hull
# of the readings' intervals; ms at 10^-3 is read between its imaged neighbours (E98's `value_at`)
Stim = Tuple[str, float]
X2_RATIOS: Tuple[Tuple[str, Tuple[Stim, ...], Stim], ...] = (
    ('r1', (('pa', -3.806), ('pa', -2.0)), ('ms', -3.0)),      # pa over ms near 10^-3 (pa at 10^-3.8 and 10^-2)
    ('r2', (('paa', -4.0), ('paa', -3.0)), ('ms', -3.0)),      # paa over ms near 10^-3 (paa at 10^-4 and 10^-3)
    ('r3', (('paa', -5.0),), ('paa', -4.0)),                   # paa at 10^-5 over 10^-4, both imaged
)
X2_QUANTITIES = ('claws', 'APL', 'Q')
X2_MIN_RATIOS = 2
# the model's stimuli for X2 (B1's labels): the numerator and the denominator estimates
X2_MODEL: Dict[str, Tuple[Tuple[str, ...], Tuple[str, ...]]] = {
    'r1': (('pa',), ('ms', 'ms [kiwi-control/diag] @ -3', 'ms [megamat/diag] @ -3', 'ms [validation2/diag] @ -3',
                     'ms [kiwi-control/control] @ -3')),
    'r2': (('paa',), ('ms', 'ms [kiwi-control/diag] @ -3', 'ms [megamat/diag] @ -3', 'ms [validation2/diag] @ -3',
                      'ms [kiwi-control/control] @ -3')),
    'r3': (('paa [validation2/diag] @ -5',), ('paa',)),
}
X2_MODEL_BESIDE = {'r3': ('paa [kiwi-control/diag] @ -5', 'paa [megamat/diag] @ -5')}

# X3
X3_PASS_A = 0.3
X3_FAIL_A = 0.15
X3_FAIL_MIN_REPEAT_R = 0.5
X3_MIN_FLIES = 5                         # of the 7
X3_UNPROBED_MIN = 2                      # of the three unopened flies
X3_MIN_KEPT = 4
X3_MIN_ODORS = 2
N_PERM = 1000
NULL_PCT = 95.0
B5_SLOPE_MULT = 2.0                      # set L's median slope below this times b_G: outcome (h)

# A4: the place test's power on made-up two-channel maps
A4_RELIABILITIES = tuple(round(x, 2) for x in np.arange(0.3, 0.9 + 1e-9, 0.1))
A4_TRUE_SHARES = (0.3, 0.5, 0.7)         # the true agreement as a share of the ceiling
A4_B_G = 0.034                           # the planted leak (the entry's Probes: green-only slopes 0.028 and 0.039)
A4_KC_OVER_APL = 0.2                     # the KC channel's odor-evoked amplitude over APL's (Probes: slopes 0.13 to 0.30)
A4_MODULATION = 0.3                      # the odor-specific part's SD over cells as a share of the map's mean
A4_GRID = 15                             # cells per side (60 um at 4 um)
A4_ODORS = (4, 3, 3)                     # kept stimuli per odor in a made-up recording
A4_RECORDINGS = 3                        # per made-up fly
A4_N_FLIES = 100                         # made-up flies per condition
A4_N_PERM = 1000
A4_POWER_MIN = 0.8
A4_SMOOTH_CELLS = 1.5                    # the spatial scale of the made-up patterns (Gaussian SD, cells)

# Yang's release and her answers (questions 2 and 13 to 18): parts B to E read nothing until these are set
YANG_RELEASED = False                    # question 13: her agreement to the use of the dual-colour recordings, both channels
CHANNELS_CONFIRMED = False               # question 2: C1 = APL (green), C2 = KCs (red)
RELEASED_STEMS: Tuple[str, ...] = ()     # if she releases only some recordings: their stems (empty: all, once released)


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def copy_dir() -> Path:
    """viki's copy of Yang's processed APL recordings (data/yang_apl_recordings/ in the shared directory)."""
    return H.shared_dir() / 'data' / 'yang_apl_recordings'


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(path.read_text()) if path.exists() else None


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(),
                                         time=time.strftime('%Y-%m-%d %H:%M'), **extra))
    _write_json(path, d)


def group_of(stem: str) -> Optional[str]:
    """The group (GROUPS) of a recording, 'left out' for LEFT_OUT, None otherwise."""
    raise NotImplementedError


def fly_of(stem: str) -> str:
    """A recording's fly, '<date>_<fly number>' ('20221207_1FLfun2' -> '20221207_1')."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A1: the table of the two-channel recordings (recordings.csv, the stimulus files and the date notes; no map values)

def block_design(frames: pd.DataFrame) -> dict:
    """A recording's design from its presentations (E100's `yaml_stimuli`): the stimuli, the presentations per stimulus,
    whether each block is one solvent then repeats of one stimulus ('solvent_then_repeats'), the order of the stimulus
    blocks, and whether each odor's blocks ascend in concentration ('ascending') or not ('shuffled')."""
    raise NotImplementedError


def note_remarks(row) -> dict:
    """A recording's remarks from recordings.csv's quality_and_notes and its date note's lines naming the run: the warp
    in um where a remark gives one, the laser wavelength (nm) where a remark gives one, and the remarks themselves."""
    raise NotImplementedError


def part_a1() -> dict:
    """A1: a1_recordings.csv (one row per two-channel recording of GROUPS, LEFT_OUT and the raw-only sets: group, fly,
    calyx, plane, pulse, stimuli, block order, gains, wavelength, Yang's remarks, whether in her list, usable) and
    a1_groups.csv (per group: recordings, flies, calyces, use)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the model: claw totals, APL, ratios (A2, A3)

def claw_totals(st_on: dict, st_off: Optional[dict], glomeruli_of_claws: np.ndarray, kc_of_claw: np.ndarray,
                input_glomeruli: Sequence[str], panel: Sequence[str]) -> pd.DataFrame:
    """Per stimulus and condition, the claw stage's totals: for each subset (CLAW_SUBSETS) and claw version
    (CLAW_SHARES), the sum over claws of the window mean ('win'), of the pre-odor mean ('pre') and of the change ('change',
    win minus pre); E90's scales (the root mean square of the claw change and of the spike counts over all claws and the
    `panel` stimuli of the APL-on run)."""
    raise NotImplementedError


def rma_line(x: np.ndarray, y: np.ndarray) -> dict:
    """The reduced-major-axis line of y on x (slope sign(r) sd(y) / sd(x), through the means): slope, intercept,
    residuals, Pearson r."""
    raise NotImplementedError


def family_line(fam: pd.DataFrame, x: str = 'claw_change', y: str = 'apl_evoked') -> pd.DataFrame:
    """A2's line per wiring: log2 of `y` on log2 of `x` over the odors (least squares and reduced major axis), each odor's
    residual (RMA); one row per wiring and odor."""
    raise NotImplementedError


def family_ratios(fam: pd.DataFrame, pairs: Sequence[Tuple[str, str]] = X1_PAIRS) -> pd.DataFrame:
    """A2's ratios per wiring: pa over ms, paa over ms, paa over pa for the claws (each subset and version, linear and the
    claw indicator's range), APL (linear, evoked, and E98's indicator range), the PN totals and the KCs' response rates
    (APL on and off)."""
    raise NotImplementedError


def part_a2() -> dict:
    """A2 (no simulation): from E90's saved stage matrices (data/model_outputs/E90_model_stages/B/), per wiring and panel
    odor: APL's window and pre-odor means and evoked activity, the claw totals (`claw_totals`), the PN stage's total change
    and the KCs' response rates (a2_family.csv); the line of log APL on log claws (a2_lines.csv); the ratios
    (a2_ratios.csv)."""
    raise NotImplementedError


class StageExtendedRun:
    """E98's `ExtendedRun` (one wiring at its fit on B1's extended input, APL's activity saved) with E90's stage reading
    after each run: `run(cond)` returns E98's outputs plus E90's 'claw_win', 'claw_pre', 'bouton_win', 'bouton_pre' (units x
    stimuli)."""

    def __init__(self, label: str, orn: Optional[pd.DataFrame] = None):
        raise NotImplementedError

    def run(self, cond: str) -> dict:
        raise NotImplementedError


def run_a3(label: str, force: bool = False) -> dict:
    """A3 on one wiring: APL on and acutely off on B1's input with the stages saved (outputs/E101_apl_vs_kc_claws/A3/),
    with the run's identities (`a3_identity`). Resumable."""
    raise NotImplementedError


def a3_identity(label: str, out: Dict[str, dict]) -> dict:
    """A3's identities on one wiring: the 16 panel odors' counts and responses equal the saved fit's (0 differences), B1's
    saved counts, responses and APL means equal (where E98 has the wiring's B1 run), and the blank's evoked APL and claw
    change (E90's pre-odor window and the settled one)."""
    raise NotImplementedError


def part_a3(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A3: the runs (`run_a3`) and the tables a3_runs.csv, a3_stimuli.csv, a3_family.csv, a3_ratios.csv,
    a3_identity.json."""
    raise NotImplementedError


def model_x2(fam: pd.DataFrame) -> pd.DataFrame:
    """The model's side of X2 per wiring, ratio (X2_MODEL) and denominator estimate: the claws' ratio (each version,
    linear and over CLAW_NS and FORMS), APL's (linear and over APL_NS and FORMS) and Q's range over every combination."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# A4: the place test's power on made-up two-channel maps

def smooth_field(rng: np.random.Generator, n: int, sigma: float) -> np.ndarray:
    """A random smooth field on an n x n grid (white noise blurred with a Gaussian of SD `sigma` cells), unit SD."""
    raise NotImplementedError


def made_up_recording(rng: np.random.Generator, r_apl: float, r_kc: float, true_corr: float, leak: float,
                      odors: Sequence[int] = A4_ODORS, n: int = A4_GRID) -> dict:
    """One made-up two-channel recording: per odor a true APL pattern and a KC pattern correlated `true_corr`, per kept
    stimulus four presentations of each channel (calyx-mean response times one plus A4_MODULATION times the pattern, plus
    per-cell noise set so that each channel's odor-specific parts have repeat r `r_apl` / `r_kc` between halves), the KC
    channel with `leak` times the APL channel added; returns the maps (presentations x cells per channel), the stimulus
    and odor of each presentation."""
    raise NotImplementedError


def part_a4(n_flies: int = A4_N_FLIES, n_perm: int = A4_N_PERM, seed: int = SEED) -> dict:
    """A4: the share of made-up flies that pass X3's rule at b = 0, b_G and 2 b_G, per pair of reliabilities and true
    agreement (a4_power.csv), and the rule that decides whether X3's "follows" is read at b_G alone (a4_rule.json)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the readers (parts B to E; viki's copy, both channels)

def released(stem: str) -> bool:
    """Whether Yang has released a recording for E101 (YANG_RELEASED and RELEASED_STEMS), or it is one of the 12 2021
    recordings she released before (data/yang_apl_calyx_maps/)."""
    raise NotImplementedError


def channel_paths(row) -> dict:
    """One recording's files in viki's copy, per channel: the peak-map stacks (C1's copied variant and C2's moved by C1's
    shifts; a single-channel recording's one stack as C1), Yang's ROI table of each channel, the outline, the yaml and
    the Experiment.xml."""
    raise NotImplementedError


def read_two_channel(row, maps: bool = True) -> dict:
    """One recording read for parts B to E: its presentations (E100's `yaml_stimuli`), both channels' stacks (one shape,
    one outline), each presentation's calyx response per channel with Yang's own value beside, pixel size and channels from
    the Experiment.xml. Refuses a raw-only row (rule of 2026-11-15 15:28) and a recording Yang has not released
    (`released`)."""
    raise NotImplementedError


def cell_grid(mask: np.ndarray, pixel_um: float, cell_um: float = CELL_UM, min_inside: float = CELL_MIN_INSIDE) -> dict:
    """The cells of `cell_um` over the outline's bounding box: each pixel's cell, the cells kept (at least `min_inside` of
    their pixels inside the outline) and their shape."""
    raise NotImplementedError


def cell_maps(stack: np.ndarray, grid: dict) -> np.ndarray:
    """Per presentation, the mean of each kept cell's inside pixels (presentations x kept cells)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the statistics (parts B to E)

def stimulus_responses(frames: pd.DataFrame, values: np.ndarray, reading: str = CARRYING) -> pd.DataFrame:
    """Per stimulus, the response under `reading` (READINGS) from per-presentation values (calyx means, or maps:
    presentations x cells): 'carrying' the mean over presentations minus the mean of the recording's first-in-block solvent
    presentations; 'block' each presentation minus its block's solvent; 'first_two' the first two presentations only, minus
    the recording's first-in-block solvent."""
    raise NotImplementedError


def leak_slopes(apl: np.ndarray, kc: np.ndarray) -> np.ndarray:
    """Per cell, the least-squares slope of the KC channel on the APL channel across stimuli (stimuli x cells each)."""
    raise NotImplementedError


def remove_leak(kc, apl, b: float):
    """K_b = C2 - b C1 (presentation by presentation, maps or calyx means)."""
    raise NotImplementedError


def kept_stimuli(resp_apl: pd.Series, resp_kc: pd.Series, solvent_apl: np.ndarray, solvent_kc: np.ndarray,
                 mult: float = KEPT_SD_MULT) -> List[str]:
    """The stimuli whose APL response and KC response (at b_G) each exceed `mult` times the SD of that channel's solvent
    responses."""
    raise NotImplementedError


def x1_recording(apl: pd.Series, kc: pd.Series, pairs: Sequence[Tuple[str, str]] = X1_PAIRS) -> dict:
    """X1 in one recording from its kept stimuli (stimulus label -> response, both positive): the RMA line of log2 APL on
    log2 K_b, each odor's mean residual, and d for each pair with X1_MIN_KEPT_PER_ODOR or more kept stimuli of each odor."""
    raise NotImplementedError


def fly_mean_interval(per_fly: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED) -> Tuple[float, float, float]:
    """The mean over flies with its 90 % interval from bootstrap draws of flies."""
    raise NotImplementedError


def x1_verdict(intervals: Dict[str, Tuple[float, float]], n_flies: Dict[str, int], bound: float = X1_BOUND,
               min_flies: int = X1_MIN_FLIES) -> str:
    """X1 from each pair's interval of d: 'one function', 'depends on the odor' or 'not decided'."""
    raise NotImplementedError


def flies_ratio(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """A ratio over flies from per-recording numerator and denominator (columns fly, num, den): recordings whose
    denominator is at most MIN_RESPONSE left out, each fly's recordings averaged on the log scale, the geometric mean over
    flies with its 90 % interval."""
    raise NotImplementedError


def flies_q(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Q, APL's ratio over the claws' ratio within the fly (columns fly, apl_num, apl_den, kc_num, kc_den), over flies."""
    raise NotImplementedError


def x2_quantity_verdict(ratio_classes: Dict[str, str], min_ratios: int = X2_MIN_RATIOS) -> str:
    """A quantity's X2 class from its three ratios' classes."""
    raise NotImplementedError


def q_classes(apl: pd.DataFrame, claws: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    """Per wiring, Q's class against the flies' interval over every combination of the two indicator families, their forms
    and the claw versions (`apl` and `claws`: per wiring the numerator's and denominator's window and pre-odor means)."""
    raise NotImplementedError


def odor_specific_parts(maps: np.ndarray, stims: Sequence[str]) -> np.ndarray:
    """The odor-specific parts of kept stimuli's maps (stimuli x cells): each map divided by its mean over cells, minus
    the mean of those maps over the stimuli."""
    raise NotImplementedError


def x3_recording(apl_halves: Tuple[np.ndarray, np.ndarray], kc_halves: Tuple[np.ndarray, np.ndarray],
                 odors: Sequence[str], n_perm: int = N_PERM, seed: int = SEED) -> dict:
    """X3 in one recording: A(b) (the mean over kept stimuli and the two crossings of the halves of the Pearson r over
    cells between APL's odor-specific part from one half and K_b's from the other), the null (N_PERM assignments giving each
    stimulus the K_b part of a stimulus of another odor; its NULL_PCT percentile), each channel's repeat r and the
    ceiling."""
    raise NotImplementedError


def fly_x3(per_recording: pd.DataFrame) -> pd.DataFrame:
    """Per fly and b: A(b) and the null percentile (means over its recordings), both channels' repeat r, and 'pass' /
    'fail' / 'neither' by X3's rule."""
    raise NotImplementedError


def x3_verdict(fly_classes: pd.DataFrame, cross_check_passed: bool, b_g_alone: bool = False) -> dict:
    """X3 from one row per fly and b multiple (columns fly, b_mult, cls): 'follows', 'follows at the controls\\' leak',
    'does not follow' or 'not decided'."""
    raise NotImplementedError


def cross_correlogram(a: np.ndarray, k: np.ndarray, cell_um: float) -> dict:
    """The normalized cross-correlogram of two odor-specific parts on the cell grid (2-D, NaN outside the outline set to 0):
    the peak's lag (cells) and the length (um) at which the radial mean falls to 1/e of the peak."""
    raise NotImplementedError


def common_map_share(parts: np.ndarray, kc_parts: Optional[np.ndarray] = None) -> dict:
    """E79's measure on normalized maps (stimuli x cells, each map over its mean): the share of APL's sum of squares (about
    each map's mean) held by one map common to all stimuli, and the share of the rest that the claws' rest accounts for."""
    raise NotImplementedError


def held_class(classes: Dict[Tuple[float, str], str], unprobed: str, n_unprobed: int, min_unprobed: int) -> str:
    """A class held over b (0, b_G) and the readings (carrying and partners), with the unprobed flies' class by the same
    rule and at least `min_unprobed` of them remaining; otherwise 'not decided'."""
    raise NotImplementedError


def _waiting(part: str) -> None:
    raise NotImplementedError(
        f"part {part} reads Yang's two-channel recordings from viki's copy only after her release (question 13) and her "
        "confirmation of the channels (question 2): YANG_RELEASED and CHANNELS_CONFIRMED; its readers and statistics are "
        "written and tested on made-up maps")


def part_0(e90: bool = False) -> dict:
    """Part 0's identities that need no new recording (identity.json)."""
    raise NotImplementedError


def part_b() -> dict:
    _waiting('B')


def part_c() -> dict:
    _waiting('C')


def part_d() -> dict:
    _waiting('D')


def part_e() -> dict:
    _waiting('E')


def score() -> dict:
    _waiting('score')


PARTS = ('0', 'A1', 'A2', 'A3', 'A4', 'B', 'C', 'D', 'E', 'score')
READS = {
    '0': "one released 2021 recording (E100's reader); with --e90, E90's saved stage matrices",
    'A1': "viki's recordings.csv, the stimulus yaml files and Yang's date notes (no map values)",
    'A2': "E90's saved stage matrices (data/model_outputs/E90_model_stages/B/)",
    'A3': "C2's family on E98's B1 input (E83's saved fits); E98's saved B1 runs for the identity",
    'A4': 'nothing (made-up maps)',
}
WRITES = {
    '0': 'identity.json', 'A1': 'a1_recordings.csv, a1_groups.csv',
    'A2': 'a2_family.csv, a2_lines.csv, a2_ratios.csv',
    'A3': 'a3_runs.csv, a3_stimuli.csv, a3_family.csv, a3_ratios.csv, a3_lines.csv, a3_identity.json; '
          'outputs/E101_apl_vs_kc_claws/A3/',
    'A4': 'a4_power.csv, a4_rule.json',
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--wirings', default=','.join(ALL_WIRINGS), help='A3: comma-separated wirings')
    ap.add_argument('--e90', action='store_true', help="part 0: also check against E90's saved stage matrices")
    ap.add_argument('--n-flies', type=int, default=A4_N_FLIES, help='A4: made-up flies per condition')
    ap.add_argument('--n-perm', type=int, default=A4_N_PERM, help='A4: null assignments per made-up recording')
    ap.add_argument('--dry-run', action='store_true', help='say what the part reads and writes, and exit')
    a = ap.parse_args()
    if a.dry_run:
        print(f'part {a.part}: reads {READS.get(a.part, "viki copy of Yang recordings (after her release)")}; writes '
              f'{WRITES.get(a.part, "not written yet")} under {out_dir()}')
        return
    if a.part == '0':
        provenance('0', e90=a.e90)
        print(json.dumps(part_0(e90=a.e90), indent=1, default=float), flush=True)
    elif a.part == 'A1':
        provenance('A1')
        print(json.dumps(part_a1(), indent=1, default=float), flush=True)
    elif a.part == 'A2':
        provenance('A2')
        print(json.dumps(part_a2(), indent=1, default=float), flush=True)
    elif a.part == 'A3':
        w = tuple(x for x in a.wirings.split(',') if x)
        provenance('A3', wirings=list(w))
        print(json.dumps(part_a3(w), indent=1, default=float), flush=True)
    elif a.part == 'A4':
        provenance('A4', n_flies=a.n_flies, n_perm=a.n_perm)
        print(json.dumps(part_a4(a.n_flies, a.n_perm), indent=1, default=float), flush=True)
    else:
        {'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e, 'score': score}[a.part]()


if __name__ == '__main__':
    main()
