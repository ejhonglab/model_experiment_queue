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
- The channels: `_C1` is APL (GCaMP6f, green PMT) and `_C2` the KCs (jRGECO1a, red PMT), as viki's copy infers. Yang's
  question 2 is open (YANG_CHANNELS_CONFIRMED); by Tom's decision of 2026-11-15 17:33 (CHANNEL_DECISION) E101 reads them
  as inferred, checked by the controls first (CHANNEL_MODE; `channel_check`, part B): parts C to E and the score run only
  if the check passed.
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
- B1's cross-check: if it fails, the slope does not carry between flies, and every KC quantity at b_G (X1; X2's claws and
  Q) is not decided (`bg_undecided`), so none of them is held; X3's "follows" needs the cross-check passed. The classes as
  computed at b_G are reported beside (if_slope_carried), deciding nothing.
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
  on the 16 panel odors (`e90_claw_identity`; reads E90's saved stage matrices); in place of "a blank evokes nothing"
  (failed as written, the amendment of 2026-11-15 18:08), `a3_blank_checks` on every A3 run. After Yang's release: each channel's calyx means against Yang's own table of that channel; each
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
- B (after Yang's release): first the channel check (`channel_check`; channel_check.csv, channel_check.json), then B1
  group G's slope per cell, b_G per fly and pooled, and the cross-check between the two green-only flies; B2 group R's APL
  channel where APL carries no indicator; B3 group M's red marker; B4 set L's solvent maps; B5 set L's own slope against
  b_G (below twice b_G: outcome (h)). Part B runs whole whatever the check says.
- C (after part B's channel check passed): C1 per recording and stimulus APL's response and K_b; C2 X1's statistic; C3
  X2's three ratios for the claws, APL and Q; C4 beside (V0, L+, VT, Hallem's input, APL off, Yang's single-colour group
  tables).
- D (after the check passed): X3's statistic per recording of set L, with its null, ceiling, cross-correlogram and E79's
  common-map share.
- E (after the check passed; described): the claws' own odor-odor distances from K_b's maps (e_claw_distances.csv), with
  the family's claw stage and Remy's claws beside (e_model_distances.csv).
- score (after the check passed): the verdicts (verdict.json, score_x2.csv, score_x3_flies.csv); X2 under the model's
  carrying level (A0_CARRYING) with its partner (A0_PARTNER) and E90's level shown (A0_SHOWN), with a verdict line per
  quantity; the readings beside.
- A2-blank (no simulation; from A3's family table): A2's line and ratios again at the model's carrying level, with E90's
  level beside, APL on and acutely off (a2_blank_lines.csv, a2_blank_ratios.csv, a2_blank.json).
- C4-model (no simulation; from A3's family table): C4's model side, each condition at its own run's blank: Hallem's
  input for pa and ms along concentration (c4_hallem_levels.csv, c4_hallem_ratios.csv) and X2's ratios with APL acutely
  off beside APL on (c4_apl_off.csv); c4_model.json.
- C-strengths (after the check passed): the flies' calyx responses per stimulus in sets L and V0, APL's and K_b's, and
  each over its odor's 10^-2 (c_strengths.csv; described).
- variants: variants.csv, the archive's params.csv source (one row per wiring's A3 runs, the absolute parameters).

The channel check (Tom, 2026-11-15 17:33; thresholds fixed in code before any map of the controls was read, `channel_check`):
per two-channel recording of the controls, under the carrying reading of calyx responses, the channel that must respond
(CHANNEL_EXPECTED: C1 in group G, the green-only flies; C2 in group R, the red-only flies) responds if its response to its
strongest stimulus exceeds CHANNEL_RESPOND_SD (2) times the SD of its solvent responses and CHANNEL_RESPOND_MIN (0.02 dF/F);
the other channel stays at the leak if its response to that stimulus is at most CHANNEL_OTHER_MAX (0.15) of it in absolute
value. A recording is 'as inferred' if both hold, 'swapped' if the same holds with the channels' roles exchanged, and
'neither' otherwise (also when both hold). A control fly is 'as inferred' ('swapped') if at least CHANNEL_FLY_MIN_SHARE
(2/3) of its two-channel recordings are; otherwise 'neither'. The check passes if every control fly of groups G and R with
a two-channel recording is 'as inferred', with at least one such fly in each group; it reads 'swapped' if every one is
'swapped'; otherwise it fails. Why these numbers: the entry's Probes give the green-only flies' KC channel at 0.028 and
0.039 of the APL channel (slopes) and the red-only fly's APL channel at 0.015 dF/F at most, against slopes of 0.13 to 0.30
of the KC channel on the APL channel in the dual-colour flies; 0.15 lies about four times above the leak and below the
dual-colour flies' coupling, and a swap puts the ratio near 1 / 0.035.

Deviations from the entry (each also in the note):
- A3's stimuli are B1's as `b1_input` gives them: all 52 of `b1_specs` since E98's 11884b2 (Hallem's rows from olfsysm's
  hc_data.csv). A first A3 run on the 46 stimuli of E98's first B1 runs (Hallem's six rows below 10^-2 without a source)
  was made before that commit and is kept apart (outputs/E101_apl_vs_kc_claws/A3_b1_46/); the tables read the 52.
- A3's check against E90's saved claw-stage sums and A2 read E90's saved stage matrices; the code subagent that wrote this
  could not read them (the permission system refused the read, as it had the planning subagent's), so `--part A2` and
  `--part 0 --e90` are left for hal's session. A3 reproduces the saved fit's counts on the 16 panel odors (E83's
  `saved_fit`) and E98's B1 where saved.
- The blank's evoked APL is read against the settled pre-odor level (E98's SETTLED_PRE_WINDOW, -0.25 to 0 s) as well as
  E90's pre-odor window, which holds the KC layer's start-up (E98's part 0 found 0.51 with E90's window on the real wiring).
- Maps in part D are the stimulus maps minus the recording's solvent map (the mean of its first-in-block solvent maps), as
  the carrying reading subtracts the solvent's calyx response; the partners subtract the block's solvent map, or use the
  first two presentations (whose halves are then the 1st and the 2nd).
- X3's odor-specific parts divide each map by its mean over cells; in the controls whose KC channel carries no indicator
  (B1's cross-check between the green-only flies, B3's marker), the corrected KC maps have no positive mean, so their parts
  are the maps minus their mean over the stimuli, without the division (Pearson's r does not depend on a map's scale).
- The kept stimuli of a recording are found once, under the carrying reading at b_G, and used for every reading and b.
- A4's made-up flies: 3 recordings each, 10 kept stimuli of 3 odors (4, 3, 3) on a 15 x 15 grid of cells (A4_*); the
  planted leak is A4_B_G (0.034, the planner's green-only slopes 0.028 and 0.039, the entry's Probes) with the KC channel's
  odor-evoked amplitude A4_KC_OVER_APL (0.2) of APL's; "the made-up flies" of A4's rule are pooled over the whole grid of
  reliabilities (0.3 to 0.9 per channel, every pair equally weighted), the literal reading of the entry.
- The channels are read as inferred while Yang's question 2 is open, by Tom's decision of 2026-11-15 17:33 (the entry
  asked for her confirmation first); part B's channel check takes its place as a gate (above). The check is per control
  fly: every green-only and red-only fly with a two-channel recording must agree, the literal reading of "in the green-only
  flies ... in the red-only fly" (group R has two flies, 2022-07-15 fly 4's right calyx and 2022-10-05 fly 2).
- The model's pre-odor level for X2 (the amendment of 2026-11-15 18:08, planned by a Fable subagent from hal 18:01 before
  the flies' X2 values were seen; open issue 3e15a91): the carrying level is the blank's window mean in the same run and
  condition, for the claws (each subset and version) and APL ('blank'); the partner is APL's settled pre-odor level with
  the claws at the blank's window ('settled'; the claws have no second level, A3 having saved none); E90's pre-odor mean
  ('e90', the entry's literal reading) is shown on every verdict line and decides nothing. APL's and Q's classes are held
  only if the partner gives the same class under the flies' carrying reading (at b = 0 and b_G; `partner_held`). The
  earlier 'settled' reading of the claws (at E90's level) is gone. Part 0's item "a blank evokes nothing" is recorded as
  failed as written, its stop lifted, and three checks on every A3 run take its place (`a3_blank_checks`).
- Yang's answer of 2026-11-15 18:08 (questions 1 and 13): the data can be used, but she has not checked the dual-colour
  recordings' quality; QUALITY_CAVEAT heads every results table of parts B to E and the score (`write_table`) and is in
  their json files.
- Part E adds the family's claw stage (A3: the correlation distance over the claws of the input's glomeruli, window minus
  E90's pre-odor mean, in E90's three claw versions) and Remy's claws for pa and ms (E90's readers, as measured and with the
  trial noise removed) as e_model_distances.csv, beside the flies'.
- Parts C to E skip a recording that lacks one of the two channels' maps (listed in their outputs).
- The green-only fly 2022-10-05 fly 1 (group G) is written '1MR' in Yang's processed files and '1FR' in its raw run folders,
  whose Experiment.xml files viki copied; its pixel size is read from the run with the same date, fly, side and run number
  (`experiment_xml_path`), with the frame size checked against the maps.
- Yang answered question 18 (2026-11-15 17:53): no run is left out, so the runs at 925 nm and those marked as warped by 5 um
  or more stay in every carrying reading, as the entry has them; the score reports the sensitivity readings without them
  beside (and the planes apart, V0 as the repetition, X3 at 2 and 8 um), deciding nothing (`beside`). Whether C2 moved by
  C1's shifts is aligned with C1 pixel for pixel is not answered; the score says where X3 depends on it (ALIGNMENT_NOTE).

CLI (from the repository root):
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A1
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A2              # reads E90's saved stage matrices
    env OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E101_apl_vs_kc_claws --part A3
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A4
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part B           # after Yang's release (YANG_RELEASED)
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part C           # C, D, E and score: after the check passed
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part score
    ./venv/bin/python -u -m experiments.E101_apl_vs_kc_claws --part A2-blank    # and C4-model, C-strengths, variants
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
SENSITIVITY_WARP_UM = 5.0                # the sensitivity readings leave out runs Yang's notes mark as warped this much or more
SENSITIVITY_NM = 925.0                   # or as taken at this wavelength
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
YANG_RELEASED = True                     # question 13: "yes, all" (2026-11-15 17:13, questions_for_yang/answers_2026-11-15_1614__2026-11-15_171352.txt)
YANG_CHANNELS_CONFIRMED = False          # question 2: C1 = APL (green), C2 = KCs (red); open on 2026-11-15
RELEASED_STEMS: Tuple[str, ...] = ()     # if she releases only some recordings: their stems (empty: all, once released)

# how the channels are read (replaces the gate on Yang's confirmation, 2026-11-15): 'waiting' (nothing read), 'inferred,
# checked by the controls' (part B's channel check gates parts C to E and the score), 'confirmed by Yang'
CHANNEL_MODES = ('waiting', 'inferred, checked by the controls', 'confirmed by Yang')
INFERRED_CHECKED = CHANNEL_MODES[1]
CHANNEL_MODE = INFERRED_CHECKED
CHANNEL_DECISION = ("Tom, 2026-11-15 17:33 (hal's terminal), while Yang's question 2 is open: E101 goes ahead with C1 = APL "
                    "(GCaMP6f) and C2 = KCs (jRGECO1a), as inferred, and part B's controls check it first: in the green-only "
                    "flies (group G) C1 must respond and C2 carry only the leak, in the red-only fly (group R) C2 must respond "
                    "and C1 stay near zero; otherwise E101 stops after part B and reports.")
# the channel check's thresholds, fixed 2026-11-15 before any map of the controls was read (module docstring)
CHANNEL_EXPECTED = {'G': ('C1', 'C2'), 'R': ('C2', 'C1')}   # (the channel that must respond, the one that must not)
CHANNEL_RESPOND_SD = 2.0                 # responds: above this many SDs of its solvent responses ...
CHANNEL_RESPOND_MIN = 0.02               # ... and above this (dF/F; E98's MIN_RESPONSE)
CHANNEL_OTHER_MAX = 0.15                 # the other channel at that stimulus: |response| at most this share of it
CHANNEL_FLY_MIN_SHARE = 2.0 / 3.0        # a control fly agrees when at least this share of its recordings do

# the model's pre-odor level for X2, fixed by the entry's amendment of 2026-11-15 18:08 (planned by a Fable subagent,
# claude-fable-5-2, from hal 18:01, before the flies' X2 values were seen; open issue 3e15a91 in unresolved_issues.md):
# carrying, the blank's window mean in the same run and condition for the claws and APL; partner, APL's settled pre-odor
# level with the claws at the blank's window (held: APL's and Q's classes only if the partner gives the same); E90's
# pre-odor mean shown on every verdict line, deciding nothing
A0_CARRYING = 'blank'
A0_PARTNER = 'settled'
A0_SHOWN = 'e90'
A0_ROLE = {A0_CARRYING: 'carrying', A0_PARTNER: 'partner', A0_SHOWN: "E90's level, decides nothing"}
PARTNER_QUANTITIES = ('APL', 'Q')        # the claws have no second level (A3 saved no settled level for them)
# Yang's answer of 2026-11-15 18:08 (questions 1 and 13): the data can be used, but she would have to check their quality
QUALITY_CAVEAT = "Yang has not checked these recordings' quality (2026-11-15 18:08)"


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


def write_table(df: pd.DataFrame, path: Path) -> None:
    """A results table of parts B to E or the score: its first line the quality caveat (`# ` QUALITY_CAVEAT), then the
    csv (`read_table` reads it)."""
    with open(path, 'w') as f:
        f.write(f'# {QUALITY_CAVEAT}\n')
        df.to_csv(f, index=False)


def read_table(path: Path) -> pd.DataFrame:
    """A results table written by `write_table` (its caveat line skipped), or one without it."""
    with open(path) as f:
        first = f.readline()
    return pd.read_csv(path, skiprows=1 if first.startswith('#') else 0)


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: the part's run record (commits, host, time, settings), earlier runs kept."""
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script=H.script_path(__file__))
    d.setdefault('runs', []).append(dict(commits=dict(H.code_commits()), host=_host(),
                                         time=time.strftime('%Y-%m-%d %H:%M'), **extra))
    _write_json(path, d)


def data_provenance() -> dict:
    """For parts B to E and the score: the sha256 of viki's recordings.csv, the channel mode and Tom's decision."""
    import hashlib
    f = copy_dir() / 'recordings.csv'
    return dict(recordings_csv_sha256=hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None,
                channel_mode=CHANNEL_MODE, channel_decision=CHANNEL_DECISION,
                yang_channels_confirmed=YANG_CHANNELS_CONFIRMED)


def group_of(stem: str) -> Optional[str]:
    """The group (GROUPS) of a recording, 'left out' for LEFT_OUT, None otherwise."""
    for g, stems in GROUPS.items():
        if stem in stems:
            return g
    return 'left out' if stem in LEFT_OUT else None


def _stem_parts(stem: str) -> dict:
    """'20221207_1FLfun2' -> date '20221207', fly '1', sex 'F', side 'L', run '2', run_id '1FLfun2' (sex and side may be
    missing, '20220908_1Ffun1'; a suffix is kept apart, '20230126_1FLfun1_000')."""
    import re
    m = re.fullmatch(r'(?P<date>\d{8})_(?P<fly>\d+)(?P<sex>[FM]?)(?P<side>[LR]?)fun(?P<run>\d+)(?P<suffix>.*)', str(stem))
    if m is None:
        raise ValueError(f'not a recording stem: {stem!r}')
    return dict(m.groupdict(), run_id=f"{m.group('fly')}{m.group('sex')}{m.group('side')}fun{m.group('run')}")


def fly_of(stem: str) -> str:
    """A recording's fly, '<date>_<fly number>' ('20221207_1FLfun2' -> '20221207_1')."""
    p = _stem_parts(stem)
    return f"{p['date']}_{p['fly']}"


def calyx_of(stem: str) -> str:
    """A recording's calyx, '<date>_<fly number><side as Yang writes it>' ('20221207_1FLfun2' -> '20221207_1L')."""
    p = _stem_parts(stem)
    return f"{p['date']}_{p['fly']}{p['side']}"


# ---------------------------------------------------------------------------
# A1: the table of the two-channel recordings (recordings.csv, the stimulus files and the date notes; no map values)

def block_design(frames: pd.DataFrame) -> dict:
    """A recording's design from its presentations (E100's `yaml_stimuli`: frame, block, position, odor, log10_conc; the
    solvent as E100's SOLVENT): the stimuli ('<odor> @ <conc>'), the presentations per stimulus, whether every block is
    either all solvent or one solvent then repeats of one stimulus ('solvent_then_repeats'; otherwise 'other'), the order
    of the stimuli's first presentations, and whether each odor's stimuli are first presented in ascending concentration
    ('ascending'; otherwise 'shuffled')."""
    from experiments import E100_apl_weak_odors as E100
    f = frames.sort_values('frame').reset_index(drop=True)
    sol = f['odor'].astype(str) == E100.SOLVENT
    lab = [E100.SOLVENT if s else E100.stim_label((o, float(c))) for s, o, c in zip(sol, f['odor'], f['log10_conc'])]
    f = f.assign(stim=lab, is_solvent=sol)
    odor_rows = f[~f['is_solvent']]
    order = list(dict.fromkeys(odor_rows['stim']))
    design = 'solvent_then_repeats'
    for _, g in f.groupby('block', sort=True):
        g = g.sort_values('position')
        if g['is_solvent'].all():
            continue
        if not (g['is_solvent'].iloc[0] and not g['is_solvent'].iloc[1:].any() and g['stim'].iloc[1:].nunique() == 1):
            design = 'other'
            break
    per_odor = {}
    for s in order:
        o, c = E100.parse_stim(s)
        per_odor.setdefault(o, []).append(c)
    ascending = all(list(cs) == sorted(cs) for cs in per_odor.values())
    n_pres = odor_rows.groupby('stim').size()
    return dict(stimuli=order, n_stimuli=len(order), presentations=dict(n_pres.reindex(order).astype(int)),
                n_solvent=int(sol.sum()), design=design, order='ascending' if ascending else 'shuffled',
                odors={o: sorted(cs) for o, cs in per_odor.items()})


WAVELENGTH_RE = r'(?<![\d.])(9\d{2}|10\d{2})(?=\s*(?:nm|@|!|\s+for\b))'
WARP_RE = r'warp\w*[^0-9\n|]{0,30}?(<)?\s*~?\s*(\d+(?:\.\d+)?)\s*um'
RUN_HEADER_RE = r'^\s*(\d+[FM][LR]fun\d+)\b'


def date_note_runs(text: str) -> Dict[str, dict]:
    """Yang's date note read run by run: the laser wavelength in effect for each run and the run's remarks. The note's
    header sets the wavelength (its 'MaiTai <nm>@' line); a run section starts at a line '<fly><sex><side>fun<run> ...'. A
    wavelength named in the run's header line or in its settings (before its channelA sequence) is the run's; one named
    after its sequence holds from the next run on (Yang writes "Change to 925 for testing" before the run it applies to).
    The run's remarks are its section's lines that are not settings or the sequence. A line 're-run as <run>' after a run's
    sequence starts that run with the same settings. Parsed text: check against the note."""
    import re
    lines = text.splitlines()
    state = None
    for ln in lines:
        if re.match(RUN_HEADER_RE, ln):
            break
        m = re.search(r'(?i)mai\s*tai\s*(\d{3,4})\s*@', ln)
        if m:
            state = int(m.group(1))
    runs: Dict[str, dict] = {}
    cur, phase = None, None
    setting = re.compile(r'(?i)^\s*(pockels?|gain|z at|frame rate|average_num|zoom|two-way|const |seq\d*|\{|\}|[\d,;\s{}]+$)')
    for ln in lines:
        m = re.match(RUN_HEADER_RE, ln)
        if m:
            if cur is not None and runs[cur]['wavelength_nm'] is None:
                runs[cur]['wavelength_nm'] = state
            cur, phase = m.group(1), 'settings'
            runs[cur] = dict(wavelength_nm=None, remarks=[], header=ln.strip())
            w = re.search(WAVELENGTH_RE, ln[m.end():])
            if w:
                state = int(w.group(1))
            rest = ln[m.end():].strip()
            rest = re.sub(r'^\d+\s*x\s*\d+\s*', '', rest)
            if rest:
                runs[cur]['remarks'].append(rest)
            continue
        if cur is None:
            continue
        rr = re.search(r'(?i)re-?\s*run\s+as\s+(\d+[FM][LR]fun\d+)', ln)
        if rr and phase == 'after':
            # 're-run as 2FLfun3': the new run repeats the current one's settings
            prev = cur
            if runs[prev]['wavelength_nm'] is None:
                runs[prev]['wavelength_nm'] = state
            cur = rr.group(1)
            runs[cur] = dict(wavelength_nm=runs[prev]['wavelength_nm'], remarks=[f're-run of {prev}'], header=ln.strip())
            continue
        w = re.search(WAVELENGTH_RE, ln)
        if phase == 'settings':
            if w:
                state = int(w.group(1))
            if 'channelA' in ln:
                runs[cur]['wavelength_nm'] = state
                phase = 'sequence'
                continue
        elif phase == 'sequence':
            if '}' in ln:
                phase = 'after'
            continue
        elif w:
            state = int(w.group(1))
        if ln.strip() and not setting.match(ln):
            runs[cur]['remarks'].append(ln.strip())
    if cur is not None and runs[cur]['wavelength_nm'] is None:
        runs[cur]['wavelength_nm'] = state
    return runs


def warp_um(text: str) -> Optional[float]:
    """The largest warp in um that a remark gives ('warpped a lot maybe ~9um' -> 9; 'warpped for about <3um' -> 3), NaN
    when a remark names a warp without a size, None when none is named."""
    import re
    if not isinstance(text, str) or 'warp' not in text.lower():
        return None
    v = [float(m.group(2)) for m in re.finditer(WARP_RE, text, flags=re.IGNORECASE)]
    return max(v) if v else np.nan


def note_remarks(row) -> dict:
    """A recording's remarks: recordings.csv's quality_and_notes (viki's parse of Yang's yaml Note and the date note's
    remarks for the run) and the run's section of its date note (`date_note_runs`; the date folder's '<date>.txt', else the
    first other note there that has the run): the laser wavelength in effect (nm), the largest warp named (um), and the
    remarks."""
    stem = str(row['stem'])
    p = _stem_parts(stem)
    folder = H.shared_dir() / str(row['shared_path']) if isinstance(row.get('shared_path'), str) else None
    run = None
    note = ''
    if folder is not None and folder.is_dir():
        cands = [folder / f"{p['date']}.txt"] + sorted(x for x in folder.glob('*.txt')
                                                       if x.name != f"{p['date']}.txt" and 'powermeasurement' not in x.name)
        for c in cands:
            if c.exists():
                runs = date_note_runs(c.read_text(errors='replace'))
                if p['run_id'] in runs:
                    run, note = runs[p['run_id']], c.name
                    break
    q = row.get('quality_and_notes')
    q = q if isinstance(q, str) else ''
    parts = [x.strip() for x in q.split(' | ')] + (run['remarks'] if run else [])
    parts = [x[len('date note:'):].strip() if x.startswith('date note:') else x for x in parts]
    remarks = ' | '.join(dict.fromkeys(x for x in parts if x))
    return dict(wavelength_nm=run['wavelength_nm'] if run else None, warp_um=warp_um(remarks), remarks=remarks,
                date_note=note)


def a1_rows(recs: pd.DataFrame) -> pd.DataFrame:
    """A1's rows: the recordings of GROUPS and LEFT_OUT, and the two-channel recordings that are raw only (the raw-only
    sets RAW_ONLY_SETS and the unprocessed runs of D, V and DEV with a dual-colour label), with their group ('raw only' for
    the last) and fly / calyx."""
    from experiments import E100_apl_weak_odors as E100
    named = set(s for v in GROUPS.values() for s in v) | set(LEFT_OUT)
    dual = recs['content'].astype(str).str.contains('dual-colour')
    raw = ~recs.apply(E100.usable_row, axis=1)
    pick = recs['stem'].isin(named) | (dual & raw & recs['set'].isin(RAW_ONLY_SETS + ('D', 'V', 'DEV')))
    out = recs[pick].copy()
    out['group'] = [group_of(s) or 'raw only' for s in out['stem']]
    out['fly_id'] = [fly_of(s) for s in out['stem']]
    out['calyx'] = [calyx_of(s) for s in out['stem']]
    out['in_yang_list'] = out['stem'].isin(GROUP_L + YANG_LIST_NOT_COPIED)
    out['usable'] = out.apply(E100.usable_row, axis=1)
    return out.reset_index(drop=True)


def part_a1() -> dict:
    """A1: a1_recordings.csv (one row per recording of `a1_rows`: group, fly, calyx, plane, pulse, blocks, the stimuli and
    their presentations, the design and block order from the stimulus yaml (usable recordings only; the raw-only ones from
    recordings.csv's `odors`), PMT gains, the laser wavelength and the warp from Yang's notes, her remarks, whether in her
    list, usable) and a1_groups.csv (per group: recordings, flies, calyces, usable recordings, planes, the odors'
    concentrations, use). Reads recordings.csv, the stimulus yaml files and the date notes; no map values."""
    from experiments import E100_apl_weak_odors as E100
    recs = pd.read_csv(copy_dir() / 'recordings.csv')
    rows = []
    for r in a1_rows(recs).to_dict('records'):
        d = dict(group=r['group'], stem=r['stem'], set=r['set'], fly=r['fly_id'], calyx=r['calyx'],
                 in_yang_list=bool(r['in_yang_list']), usable=bool(r['usable']), content=r['content'],
                 labelling_inferred=r['labelling_inferred'], genotype_as_written=r['genotype_as_written'],
                 plane_depth_um=r['plane_depth_um'], pulse_s=r['pulse_s'], nblock=r['nblock'], nstim=r['nstim'],
                 n_presentations=r['n_presentations'], channels_xml=r['channels_xml'], has_C1=r['has_C1'],
                 has_C2=r['has_C2'], peak_map_copied=r['peak_map_copied'],
                 frames_equal_presentations=r['frames_equal_presentations'], pixel_size_um=r['pixel_size_um'],
                 pmt_gain_A_B=r['pmt_gain_A_B'], note_gain=r['note_gain'], channel_identity=r['channel_identity'])
        yml = None
        if d['usable'] and isinstance(r.get('shared_path'), str):
            folder = H.shared_dir() / r['shared_path']
            yml = next((f for f in (folder / f"{r['stem']}.yaml", folder / f"{r['stem']}_copy.yaml") if f.exists()), None)
        if yml is not None:
            bd = block_design(E100.yaml_stimuli(yml, r['stem']))
            d.update(stimulus_source='yaml', stimuli='; '.join(f"{s} x{bd['presentations'][s]}" for s in bd['stimuli']),
                     n_stimuli=bd['n_stimuli'], n_solvent=bd['n_solvent'], design=bd['design'], block_order=bd['order'],
                     odors='; '.join(f"{o} {','.join(f'{c:g}' for c in cs)}" for o, cs in bd['odors'].items()))
        else:
            d.update(stimulus_source="recordings.csv's odors", stimuli=r['odors'], n_stimuli=np.nan, n_solvent=np.nan,
                     design='', block_order='', odors='')
        d.update(note_remarks(r))
        w_, nm = pd.to_numeric(d['warp_um'], errors='coerce'), pd.to_numeric(d['wavelength_nm'], errors='coerce')
        d['out_of_sensitivity_readings'] = '; '.join(x for x, ok in (
            (f"warp {w_:g} um", np.isfinite(w_) and w_ >= SENSITIVITY_WARP_UM),
            (f"{nm:g} nm", np.isfinite(nm) and abs(nm - SENSITIVITY_NM) < 1)) if ok)
        d['why_not_read'] = ('' if d['usable'] and d['group'] != 'left out' else
                             'a plane at 87 um, one channel processed' if d['group'] == 'left out' else
                             'raw only (rule of 2026-11-15 15:28)')
        rows.append(d)
    t = pd.DataFrame(rows)
    order = {g: i for i, g in enumerate(list(GROUPS) + ['left out', 'raw only'])}
    t = t.sort_values(['group', 'stem'], key=lambda s: s.map(order) if s.name == 'group' else s).reset_index(drop=True)
    t.to_csv(out_dir() / 'a1_recordings.csv', index=False)
    g = []
    for name, x in t.groupby('group', sort=False):
        u = x[x['usable']]
        g.append(dict(group=name, use=GROUP_USE.get(name, 'not read'), recordings=len(x), usable=len(u),
                      flies=x['fly'].nunique(), calyces=x['calyx'].nunique(),
                      in_yang_list=int(x['in_yang_list'].sum()),
                      planes_um=', '.join(f'{v:g}' for v in sorted(pd.to_numeric(x['plane_depth_um'], errors='coerce')
                                                                   .dropna().unique())),
                      designs=', '.join(sorted(set(u['design'].dropna().astype(str)) - {''})),
                      block_orders=', '.join(f'{k} {v}' for k, v in u['block_order'].value_counts().items()),
                      wavelengths_nm=', '.join(f'{k:g} x{v}' for k, v in
                                               pd.to_numeric(x['wavelength_nm'], errors='coerce').value_counts().items()),
                      warp_5um_or_more=int((pd.to_numeric(x['warp_um'], errors='coerce') >= SENSITIVITY_WARP_UM).sum()),
                      out_of_sensitivity_readings=int((x['out_of_sensitivity_readings'] != '').sum()),
                      stimulus_sets=' || '.join(sorted(set(u['odors'].dropna().astype(str)) - {''}))))
    gt = pd.DataFrame(g)
    gt.to_csv(out_dir() / 'a1_groups.csv', index=False)
    return dict(n_recordings=len(t), groups=gt[['group', 'recordings', 'usable', 'flies', 'calyces']].to_dict('records'))


# ---------------------------------------------------------------------------
# the model: claw totals, APL, ratios (A2, A3)

def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def claw_col(subset: str, share: float, k: str) -> str:
    """A claw total's column: 'claw_<subset>_<share x 100, 2 digits>_<win|pre|change>' ('claw_input_00_change')."""
    return f'claw_{subset}_{int(round(share * 100)):02d}_{k}'


def claw_totals(st_on: dict, st_off: Optional[dict], glomeruli_of_claws: np.ndarray, kc_of_claw: np.ndarray,
                input_glomeruli: Sequence[str], panel: Sequence[str]) -> pd.DataFrame:
    """Per stimulus and condition, the claw stage's totals: for each subset (CLAW_SUBSETS: the claws of
    `input_glomeruli`, all claws) and claw version (CLAW_SHARES), the sum over claws of the window mean ('win'), the
    pre-odor mean ('pre') and the change ('change', win minus pre), columns `claw_col`. A version with a share s of the
    KC's spike signal is E90's `claw_signal`: (1 - s) x the claw's drive / scale0 + s x its KC's spike count / scale1, the
    scales the root mean square of the claw change and of the spike counts over all claws (KCs) and the `panel` stimuli of
    the APL-on run (E90's, with the 16 panel odors); the spike signal has no pre-odor part. `st_*`: dicts with 'claw_win',
    'claw_pre' (claws x stimuli) and 'spike_counts' (KCs x stimuli, index level H.KC_ID_LEVEL), stimulus columns the same
    labels. One row per (cond, stimulus)."""
    panel = list(panel)
    ch_on = st_on['claw_win'].astype(float) - st_on['claw_pre'].astype(float)
    scales = (_rms(ch_on[panel].to_numpy()), _rms(st_on['spike_counts'][panel].to_numpy(dtype=float)))
    masks = {'input': np.isin(np.asarray(glomeruli_of_claws).astype(str), [str(g) for g in input_glomeruli]),
             'all': np.ones(len(glomeruli_of_claws), dtype=bool)}
    rows = []
    for cond, st in (('on', st_on), ('off', st_off)):
        if st is None:
            continue
        win, pre = st['claw_win'].astype(float), st['claw_pre'].astype(float)
        if len(win) != len(glomeruli_of_claws):
            raise RuntimeError(f'{len(win)} claws for {len(glomeruli_of_claws)} claw glomeruli')
        counts = st['spike_counts'].astype(float)
        ix = pd.Index(counts.index.get_level_values(H.KC_ID_LEVEL)).get_indexer(np.asarray(kc_of_claw))
        if (ix < 0).any():
            raise RuntimeError('claws whose KC is not in the counts')
        stims = list(win.columns)
        K = counts[stims].to_numpy()[ix]                       # claws x stimuli: each claw's KC's spike count
        W, P = win[stims].to_numpy(), pre[stims].to_numpy()
        tot = {}
        for subset, m in masks.items():
            for share in CLAW_SHARES:
                w_ = (1.0 - share) * W[m] / (scales[0] if share else 1.0) + (share * K[m] / scales[1] if share else 0.0)
                p_ = (1.0 - share) * P[m] / (scales[0] if share else 1.0)
                tot[claw_col(subset, share, 'win')] = w_.sum(axis=0)
                tot[claw_col(subset, share, 'pre')] = p_.sum(axis=0)
                tot[claw_col(subset, share, 'change')] = (w_ - p_).sum(axis=0)
        df = pd.DataFrame(tot, index=pd.Index(stims, name='stimulus')).reset_index()
        df.insert(0, 'cond', cond)
        df['claw_scale_drive'], df['claw_scale_spikes'] = scales
        df['n_claws_input'], df['n_claws_all'] = int(masks['input'].sum()), int(len(m))
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def family_table(label: str, run: Dict[str, dict], meta: dict, panel: Sequence[str]) -> pd.DataFrame:
    """One wiring's family table (A2 from E90's runs, A3 from its own): per condition and stimulus APL's window, pre-odor
    (E90's window) and settled pre-odor means (where saved) and evoked activity (window minus E90's pre-odor mean), the claw
    totals (`claw_totals`), the PN stage's total change (sum over glomeruli of window minus pre-odor mean) and the KCs'
    response rate (two or more spikes). `run`: {'on': ..., 'acute_off': ...} with 'apl', 'claw_win', 'claw_pre',
    'spike_counts', 'responses', 'pn_win', 'pn_pre' (stimulus columns the same labels); `meta`: 'claw_glomerulus',
    'kc_of_claw'."""
    on, off = run['on'], run['acute_off']
    ct = claw_totals(on, off, meta['claw_glomerulus'], meta['kc_of_claw'], list(H.apl16_orn_deltas().index), panel)
    rows = []
    for cond, o in (('on', on), ('off', off)):
        a = o['apl']
        pn = (o['pn_win'].astype(float) - o['pn_pre'].astype(float)).sum(axis=0)
        rate = o['responses'].astype(float).mean(axis=0)
        for s in o['claw_win'].columns:
            rows.append(dict(wiring=label, cond=cond, stimulus=str(s), apl_win=float(a.loc[s, 'win']),
                             apl_pre=float(a.loc[s, 'pre']),
                             apl_pre_settled=float(a.loc[s, 'pre_settled']) if 'pre_settled' in a.columns else np.nan,
                             apl_evoked=float(a.loc[s, 'win'] - a.loc[s, 'pre']), pn_total=float(pn[s]),
                             kc_rate=float(rate[s])))
    df = pd.DataFrame(rows).merge(ct, on=['cond', 'stimulus'], how='left', validate='one_to_one')
    return df


def rma_line(x: np.ndarray, y: np.ndarray) -> dict:
    """The reduced-major-axis line of y on x (slope sign(r) sd(y) / sd(x), through the means): slope, intercept,
    residuals (y minus the line), Pearson r, and the least-squares slope beside."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return dict(slope=np.nan, intercept=np.nan, residuals=np.full(len(x), np.nan), r=np.nan, slope_ols=np.nan)
    r = float(np.corrcoef(x, y)[0, 1])
    slope = float(np.sign(r) * np.std(y, ddof=1) / np.std(x, ddof=1))
    intercept = float(y.mean() - slope * x.mean())
    return dict(slope=slope, intercept=intercept, residuals=y - (intercept + slope * x), r=r,
                slope_ols=float(r * np.std(y, ddof=1) / np.std(x, ddof=1)))


def family_line(fam: pd.DataFrame, x: str = 'claw_input_00_change', y: str = 'apl_evoked',
                odors: Sequence[str] = tuple(H.APL16_ODORS)) -> pd.DataFrame:
    """The family's line per wiring (APL on): log2 of `y` on log2 of `x` over `odors` (the 16 panel odors) with both
    positive (reduced major axis, the least-squares slope beside), each odor's residual (RMA; NaN for an odor left out);
    one row per wiring and odor, with the line's slope, intercept, r and the number of odors used."""
    rows = []
    for w, g in fam[fam['cond'] == 'on'].groupby('wiring', sort=False):
        g = g.set_index('stimulus').reindex(list(odors))
        ok = (g[x] > 0) & (g[y] > 0)
        ln = rma_line(np.log2(g.loc[ok, x]), np.log2(g.loc[ok, y]))
        res = pd.Series(np.nan, index=g.index)
        res[ok] = ln['residuals']
        for o in g.index:
            rows.append(dict(wiring=w, odor=o, x=x, y=y, x_value=float(g.loc[o, x]), y_value=float(g.loc[o, y]),
                             residual_log2=float(res[o]), slope_rma=ln['slope'], slope_ols=ln['slope_ols'],
                             intercept_rma=ln['intercept'], r=ln['r'], n_odors=int(ok.sum())))
    return pd.DataFrame(rows)


def _claw_levels(g: pd.Series, subset: str, share: float, baseline: str, blank: Optional[pd.Series]) -> Tuple[float, float]:
    """A claw total's window value and its baseline: E90's pre-odor mean ('e90') or the blank's window ('blank')."""
    win = float(g[claw_col(subset, share, 'win')])
    if baseline == 'e90':
        return win, float(g[claw_col(subset, share, 'pre')])
    return win, float(blank[claw_col(subset, share, 'win')])


def _apl_levels(g: pd.Series, baseline: str, blank: Optional[pd.Series]) -> Tuple[float, float]:
    """APL's window mean and its baseline: E90's pre-odor mean ('e90'), the settled pre-odor mean ('settled') or the
    blank's window ('blank'); E98's three pre-odor levels."""
    if baseline == 'e90':
        return float(g['apl_win']), float(g['apl_pre'])
    if baseline == 'settled':
        return float(g['apl_win']), float(g['apl_pre_settled'])
    return float(g['apl_win']), float(blank['apl_win'])


def ratio_sets(fam_w: pd.DataFrame, num: str, den: str, baselines: Sequence[str] = ('e90',),
               cond: str = 'on') -> pd.DataFrame:
    """One wiring's ratios of stimulus `num` over `den` (in condition `cond`, APL on by default, 'off' the acute-off
    run; `fam_w` one wiring's `family_table`): APL through E98's indicator family (APL_NS, FORMS) and the claws through
    CLAW_NS and FORMS for each subset and version, each under the `baselines` (A0_CHOICES: 'blank', the carrying level,
    the blank's window for APL and the claws in the same run and condition; 'settled', the partner, APL's settled
    pre-odor level with the claws at the blank's window, A3 having saved no settled level for the claws; 'e90', E90's
    pre-odor mean for both, shown and deciding nothing; the amendment of 2026-11-15 18:08); n 1 with the 'evoked' form
    is the linear ratio. A run without a blank is refused for 'blank' and 'settled'. Long rows: quantity ('APL',
    'claws'), baseline, subset, share, n, form, ratio."""
    from experiments import E98_apl_response_strength as E98
    on = fam_w[fam_w['cond'] == cond].set_index('stimulus')
    blank = on.loc[BLANK_LABEL] if BLANK_LABEL in on.index else None
    a, b = on.loc[num], on.loc[den]
    rows = []
    for bl in baselines:
        if bl not in A0_CHOICES:
            raise ValueError(f'unknown pre-odor level {bl!r}')
        if bl in ('blank', 'settled') and blank is None:
            raise ValueError(f"the level {bl!r} needs the blank ({BLANK_LABEL!r}) in the run, which has none")
        if bl == 'settled' and not np.isfinite(a['apl_pre_settled']):
            raise ValueError("the level 'settled' needs APL's settled pre-odor level, which the run has not saved")
        a1, a0 = _apl_levels(a, bl, blank)
        a2, _ = _apl_levels(b, bl, blank)
        for form in FORMS:
            for n in APL_NS:
                rows.append(dict(quantity='APL', baseline=bl, subset='', share=np.nan, n=n, form=form,
                                 ratio=float(E98.indicator_ratio(a1, a2, a0, n, form))))
        cb = 'blank' if bl == 'settled' else bl
        for subset in CLAW_SUBSETS:
            for share in CLAW_SHARES:
                c1, c0 = _claw_levels(a, subset, share, cb, blank)
                c2, c0b = _claw_levels(b, subset, share, cb, blank)
                if not np.isclose(c0, c0b, rtol=1e-9, atol=1e-9):
                    raise RuntimeError(f'{num} / {den}: the claws\' baselines differ ({c0} / {c0b})')
                for form in FORMS:
                    for n in CLAW_NS:
                        rows.append(dict(quantity='claws', baseline=bl, subset=subset, share=share, n=n, form=form,
                                         ratio=float(E98.indicator_ratio(c1, c2, c0, n, form))))
    return pd.DataFrame(rows)


def ratio_summary(rs: pd.DataFrame) -> pd.DataFrame:
    """Per quantity, baseline, subset and version: the linear ratio (n 1, 'evoked') and the range over the indicator
    family."""
    rows = []
    for key, g in rs.groupby(['quantity', 'baseline', 'subset', 'share'], dropna=False, sort=False):
        lin = g[(g['n'] == 1.0) & (g['form'] == 'evoked')]['ratio']
        v = g['ratio'][np.isfinite(g['ratio'])]
        rows.append(dict(zip(['quantity', 'baseline', 'subset', 'share'], key), linear=float(lin.iloc[0]),
                         indicator_min=float(v.min()) if len(v) else np.nan,
                         indicator_max=float(v.max()) if len(v) else np.nan))
    return pd.DataFrame(rows)


def family_ratios(fam: pd.DataFrame, pairs: Sequence[Tuple[str, str]] = X1_PAIRS) -> pd.DataFrame:
    """A2's ratios per wiring (APL on; the panel's stimuli): pa over ms, paa over ms, paa over pa for APL and the claws
    (`ratio_sets` summarized: linear and the indicator family's range, E90's pre-odor level), the PN totals and the KCs'
    response rates with APL on and off (linear)."""
    rows = []
    for w, g in fam.groupby('wiring', sort=False):
        on = g[g['cond'] == 'on'].set_index('stimulus')
        off = g[g['cond'] == 'off'].set_index('stimulus')
        for a, b in pairs:
            key = dict(wiring=w, pair=f'{a}/{b}', numerator=a, denominator=b)
            for r in ratio_summary(ratio_sets(g, a, b)).to_dict('records'):
                rows.append(dict(key, **r))
            for q, src, col in (('PN total', on, 'pn_total'), ('KC response rate, APL on', on, 'kc_rate'),
                                ('KC response rate, APL off', off, 'kc_rate')):
                v = float(src.loc[a, col] / src.loc[b, col])
                rows.append(dict(key, quantity=q, baseline='', subset='', share=np.nan, linear=v, indicator_min=np.nan,
                                 indicator_max=np.nan))
    return pd.DataFrame(rows)


E90_EXPERIMENT = 'E90_model_stages'


def e90_run(label: str) -> Tuple[Dict[str, dict], dict]:
    """E90's saved stage run of one wiring (APL on and acutely off; `E90_model_stages.load_stage`, the local outputs then
    the archive), its outputs' odor columns abbreviated, and its claws' metadata (E90's <label>_units.json beside the
    runs, else from the wiring, `unit_meta`)."""
    from experiments import E90_model_stages as E90
    run = {}
    for _, ct in A3_CONDS:
        st = E90.load_stage('B', f'{label}_{ct}')
        if st is None:
            raise FileNotFoundError(f"E90's stage run {label}_{ct} not found")
        for k in ('responses', 'spike_counts'):
            st[k] = H.abbrev_columns(st[k])
        run[ct] = st
    for d in (H.HERE / 'outputs' / E90_EXPERIMENT / 'B', H.shared_dir() / 'data' / 'model_outputs' / E90_EXPERIMENT / 'B'):
        f = d / f'{label}_units.json'
        if f.exists():
            return run, {k: np.asarray(v) for k, v in json.loads(f.read_text()).items()}
    from experiments import E83_cross_wiring_fit as E83
    return run, {k: np.asarray(v) for k, v in unit_meta(E83.wiring_parts(label)['W']).items()}


def part_a2(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A2 (no simulation): from E90's saved stage matrices (`e90_run`), per wiring, condition and panel odor: APL's window
    and pre-odor means and evoked activity, the claw totals, the PN stage's total change and the KCs' response rates
    (`family_table`; a2_family.csv); the line of log2 APL on log2 claws over the 16 odors for each claw subset and version
    (`family_line`; a2_lines.csv); pa over ms, paa over ms and paa over pa (`family_ratios`; a2_ratios.csv)."""
    fams = []
    for w in wirings:
        run, meta = e90_run(w)
        fams.append(family_table(w, run, meta, list(H.APL16_ODORS)))
    fam = pd.concat(fams, ignore_index=True)
    out = out_dir()
    fam.to_csv(out / 'a2_family.csv', index=False)
    lines = pd.concat([family_line(fam, x=claw_col(sub, sh, 'change')) for sub in CLAW_SUBSETS for sh in CLAW_SHARES]
                      + [family_line(fam, x='pn_total')], ignore_index=True)
    lines.to_csv(out / 'a2_lines.csv', index=False)
    ratios = family_ratios(fam)
    ratios.to_csv(out / 'a2_ratios.csv', index=False)
    main = lines[lines['x'] == claw_col('input', 0.0, 'change')].drop_duplicates('wiring')
    return dict(n_rows=len(fam), slope_rma=dict(zip(main['wiring'], main['slope_rma'])),
                r=dict(zip(main['wiring'], main['r'])))


A3_CONDS = (('on', 'on'), ('off', 'acute_off'))          # (condition, its file tag; E90's and E98's)
A3_STAGE_KEYS = ('claw_win', 'claw_pre', 'bouton_win', 'bouton_pre')
A3_E98_KEYS = ('responses', 'spike_counts', 'pn_win', 'pn_pre', 'apl_trace')
APL_TOL = 1e-9                           # A3: APL's means against E98's B1 (absolute)
BLANK_LABEL = 'blank @ 0'                # E98's


class StageExtendedRun:
    """E98's `ExtendedRun` (one wiring at its fit on B1's extended input, APL's activity saved) with E90's stage reading
    after each run: `run(cond)` returns E98's outputs plus E90's 'claw_win', 'claw_pre', 'bouton_win', 'bouton_pre' (units x
    stimuli, columns the stimulus labels) and 'apl_stage' (APL's means as `read_stages` reads them, which must equal E98's
    'apl')."""

    def __init__(self, label: str, orn: Optional[pd.DataFrame] = None):
        from experiments import E98_apl_response_strength as E98
        self.er = E98.ExtendedRun(label, orn)
        self.label, self.W, self.labels = label, self.er.W, self.er.labels
        wp = np.asarray(self.er.rv.kc.wPNKC, dtype=float)
        if not (wp.shape == self.W.shape and np.array_equal(wp, self.W.to_numpy(dtype=float))):
            raise RuntimeError(f"{label}: olfsysm's wPNKC differs from the wiring's claw x bouton matrix")

    def params(self, cond: str) -> Dict[str, float]:
        return self.er.params(cond)

    def run(self, cond: str) -> dict:
        from experiments import E90_model_stages as E90
        o = self.er.run(cond)
        st = E90.read_stages(self.er.mp, self.er.rv, self.W, self.er.order, self.labels, E90.MODEL_WINDOW)
        for k in A3_STAGE_KEYS:
            o[k] = st[k]
        o['apl_stage'] = st['apl']
        d = (st['apl'][['win', 'pre']].astype(float) - o['apl'].loc[st['apl'].index, ['win', 'pre']].astype(float)).abs()
        if float(d.to_numpy().max()) > APL_TOL:
            raise RuntimeError(f"{self.label} {cond}: E90's and E98's readings of APL differ by {float(d.to_numpy().max())}")
        return o


def a3_dir() -> Path:
    return local_dir('A3')


def unit_meta(W: pd.DataFrame) -> dict:
    """A wiring's claws: each claw's glomerulus (E90's `glomerulus_of_claws`) and KC, and each bouton's glomerulus."""
    from experiments import E90_model_stages as E90
    return dict(claw_glomerulus=[str(x) for x in E90.glomerulus_of_claws(W)],
                kc_of_claw=[int(x) for x in W.index.get_level_values(H.KC_ID_LEVEL)],
                bouton_glomerulus=[str(x) for x in W.columns.get_level_values('glomerulus')])


def save_a3(d: Path, label: str, cond_tag: str, o: dict) -> None:
    """One run's outputs under `d`: E98's B1 layout (`save_b1`: responses, spike counts, PN stage, APL's means and trace)
    plus the claw and bouton stages (`<label>_<cond>_<stage>_<win|pre>.parquet`, E90's layout)."""
    from experiments import E98_apl_response_strength as E98
    E98.save_b1(d, label, cond_tag, o)
    for k in A3_STAGE_KEYS:
        x = o[k].copy()
        x.columns = [str(c) for c in x.columns]
        if isinstance(x.index, pd.MultiIndex):
            x.index = x.index.set_names([str(n) for n in x.index.names])
        x.to_parquet(d / f'{label}_{cond_tag}_{k}.parquet')


def load_a3_run(label: str, d: Optional[Path] = None) -> Optional[Dict[str, dict]]:
    """One wiring's saved A3 run ({cond_tag: {'apl', 'responses', 'spike_counts', 'pn_win', 'pn_pre', 'claw_win',
    'claw_pre', 'bouton_win', 'bouton_pre'}}, and 'meta' the units' metadata), or None."""
    d = a3_dir() if d is None else Path(d)
    if not all((d / f'{label}_{ct}_claw_win.parquet').exists() for _, ct in A3_CONDS):
        return None
    out = {}
    for _, ct in A3_CONDS:
        o = {k: pd.read_parquet(d / f'{label}_{ct}_{k}.parquet') for k in
             ('responses', 'spike_counts', 'pn_win', 'pn_pre') + A3_STAGE_KEYS}
        o['apl'] = pd.read_csv(d / f'{label}_{ct}_apl.csv', index_col=0)
        out[ct] = o
    out['meta'] = {k: np.asarray(v) for k, v in json.loads((d / f'{label}_units.json').read_text()).items()}
    return out


def a3_identity(label: str, out: Dict[str, dict]) -> dict:
    """A3's identities on one wiring (`out`: {cond_tag: outputs}): the 16 panel odors' counts and responses equal the
    saved fit's (E83's `saved_fit`; 0 differences); B1's saved counts and responses equal on every stimulus both carry (0
    differences) and its APL means (within APL_TOL), where E98 has the wiring's B1 run (`load_b1_run`); the blank's evoked
    APL (window minus E90's pre-odor mean, and minus the settled pre-odor mean) and its claw change summed over the claws
    of the input's glomeruli and over all claws (window minus E90's pre-odor mean)."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments import E98_apl_response_strength as E98
    saved = E83.saved_fit(label)['saved']
    b1 = E98.load_b1_run(label)
    res = dict(wiring=label, b1_dir=None if b1 is None else b1['on']['dir'])
    for c, ct in A3_CONDS:
        o = out[ct]
        rs, ss = (H.abbrev_columns(x.astype(float)) for x in saved[c])
        s16 = o['spike_counts'][list(H.APL16_ODORS)].astype(float)
        r16 = o['responses'][list(H.APL16_ODORS)].astype(float)
        idx_ok = bool(len(s16) == len(ss) and np.array_equal(np.asarray(s16.index), np.asarray(ss.index)))
        row = dict(index_equal=idx_ok,
                   counts_differ=int((s16.to_numpy() != ss[list(H.APL16_ODORS)].to_numpy()).sum()) if idx_ok else -1,
                   responses_differ=int((r16.to_numpy() != rs[list(H.APL16_ODORS)].to_numpy()).sum()) if idx_ok else -1)
        if b1 is not None:
            bo = b1[ct]
            common = [s for s in o['spike_counts'].columns if s in bo['spike_counts'].columns]
            same_idx = bool(np.array_equal(np.asarray(o['spike_counts'].index), np.asarray(bo['spike_counts'].index)))
            row.update(b1_stimuli_compared=len(common), b1_index_equal=same_idx,
                       b1_counts_differ=int((o['spike_counts'][common].astype(float).to_numpy()
                                             != bo['spike_counts'][common].astype(float).to_numpy()).sum()) if same_idx else -1,
                       b1_responses_differ=int((o['responses'][common].astype(float).to_numpy()
                                                != bo['responses'][common].astype(float).to_numpy()).sum()) if same_idx else -1,
                       b1_apl_win_max_abs_diff=float((o['apl'].loc[common, 'win'] - bo['apl'].loc[common, 'win']).abs().max()),
                       b1_apl_pre_max_abs_diff=float((o['apl'].loc[common, 'pre'] - bo['apl'].loc[common, 'pre']).abs().max()))
        a = o['apl']
        if BLANK_LABEL in a.index:
            b = a.loc[BLANK_LABEL]
            row.update(blank_apl_evoked_e90_pre=float(b['win'] - b['pre']),
                       blank_apl_evoked_settled_pre=float(b['win'] - b['pre_settled']) if 'pre_settled' in b.index else None,
                       blank_response_rate=float(o['responses'][BLANK_LABEL].astype(float).mean()))
            ch = o['claw_win'][BLANK_LABEL].astype(float) - o['claw_pre'][BLANK_LABEL].astype(float)
            m = np.isin(np.asarray(out['meta']['claw_glomerulus']), list(H.apl16_orn_deltas().index))
            row.update(blank_claw_change_input=float(ch.to_numpy()[m].sum()), blank_claw_change_all=float(ch.sum()),
                       blank_claw_win_input=float(o['claw_win'][BLANK_LABEL].astype(float).to_numpy()[m].sum()))
        res[c] = row
    res['reproduces_fit'] = bool(all(res[c]['counts_differ'] == 0 and res[c]['responses_differ'] == 0 for c, _ in A3_CONDS))
    if b1 is not None:
        res['reproduces_b1'] = bool(all(res[c]['b1_counts_differ'] == 0 and res[c]['b1_responses_differ'] == 0
                                        and res[c]['b1_apl_win_max_abs_diff'] <= APL_TOL
                                        and res[c]['b1_apl_pre_max_abs_diff'] <= APL_TOL for c, _ in A3_CONDS))
    return res


def _process_commits() -> Dict[str, object]:
    from experiments import E98_apl_response_strength as E98
    return E98.process_commits()


def run_a3(label: str, force: bool = False) -> dict:
    """A3 on one wiring: APL on and acutely off on B1's input (E98's `b1_input`) with the stages saved
    (outputs/E101_apl_vs_kc_claws/A3/, `save_a3`), the units' metadata (<label>_units.json), the stimulus table with every
    input vector (a3_stimuli.csv), one row per condition in runs.csv (the absolute parameters, the input's hash, the
    commits read at the process's start, host, time) and the run's identities (`a3_identity`, <label>_identity.json).
    Resumable: a wiring whose runs exist is skipped unless `force`."""
    import gc
    import os
    from experiments import E98_apl_response_strength as E98
    d = a3_dir()
    if not force and load_a3_run(label, d) is not None:
        return _read_json(d / f'{label}_identity.json') or a3_identity(label, load_a3_run(label, d))
    t0 = time.time()
    orn, table = E98.b1_input()
    table.to_csv(d / 'a3_stimuli.csv', index=False)
    sr = StageExtendedRun(label, orn)
    (d / f'{label}_units.json').write_text(json.dumps(unit_meta(sr.W)) + '\n')
    commits = _process_commits()
    for c, ct in A3_CONDS:
        t1 = time.time()
        o = sr.run(c)
        save_a3(d, label, ct, o)
        row = dict(wiring=label, cond=ct, n_stimuli=len(sr.labels), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                   runtime_s=time.time() - t1, threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=E98._frame_sha(orn), **{f'param_{k}': v for k, v in sr.params(c).items()}, **commits)
        E98._append_csv(d / 'runs.csv', row)
        print(f'{label} {ct}: {time.time() - t1:.0f} s', flush=True)
        del o
    del sr
    gc.collect()
    ident = a3_identity(label, load_a3_run(label, d))
    _write_json(d / f'{label}_identity.json', ident)
    print(f"{label}: A3 done in {time.time() - t0:.0f} s; reproduces the fit: {ident['reproduces_fit']}; B1: "
          f"{ident.get('reproduces_b1', 'no B1 run')}", flush=True)
    return ident


def part_a3_runs(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A3's runs (`run_a3`, the real wiring first), each checked; a wiring that does not reproduce its fit or B1 stops the
    part. The provenance (`H.model_provenance` with B1's input) goes to outputs/.../A3/model_provenance.json."""
    from experiments import E98_apl_response_strength as E98
    d = a3_dir()
    if not (d / 'model_provenance.json').exists():
        _write_json(d / 'model_provenance.json', dict(script=H.script_path(__file__),
                                                      model=H.model_provenance(orn_deltas=E98.b1_input()[0])))
    out = {}
    for w in wirings:
        ident = run_a3(w)
        if not ident.get('reproduces_fit') or ident.get('reproduces_b1') is False:
            raise RuntimeError(f'{w}: the A3 run does not reproduce the saved fit or B1: {ident}')
        out[w] = ident
    return out


def part_a3(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A3: the runs (`part_a3_runs`) and the tables (`a3_tables`): a3_runs.csv, a3_stimuli.csv, a3_family.csv,
    a3_ratios.csv, a3_lines.csv, a3_identity.json."""
    ident = part_a3_runs(wirings)
    _write_json(out_dir() / 'a3_identity.json', ident)
    return dict(apl=a3_apl_tables(wirings), claws=a3_tables(wirings))


def model_x2(fam: pd.DataFrame, baselines: Sequence[str] = ('e90', 'settled', 'blank')) -> pd.DataFrame:
    """The model's side of X2 (APL on; `fam` the A3 family table): per wiring, ratio (X2_MODEL; X2_MODEL_BESIDE with
    role 'beside') and estimate of the input (the denominator's for r1 and r2, the numerator's for r3), `ratio_sets`
    under each baseline. Long rows; Q's set is the quotient of every APL row by every claw row of the same wiring, ratio,
    estimate and baseline (`q_classes`)."""
    rows = []
    for w, g in fam.groupby('wiring', sort=False):
        have = set(g['stimulus'])
        for rid, (nums, dens) in X2_MODEL.items():
            pairs = [(nu, de, 'main') for nu in nums for de in dens]
            pairs += [(nu, dens[0], 'beside') for nu in X2_MODEL_BESIDE.get(rid, ())]
            for nu, de, role in pairs:
                if nu not in have or de not in have:
                    continue
                rs = ratio_sets(g, nu, de, baselines)
                est = de if rid in ('r1', 'r2') else nu
                rows.append(rs.assign(wiring=w, ratio_id=rid, numerator=nu, denominator=de, estimate=est, role=role))
    return pd.concat(rows, ignore_index=True)


def x2_model_summary(sets: pd.DataFrame) -> pd.DataFrame:
    """Per wiring, ratio, estimate and baseline: APL's and the claws' linear ratios (claws: input subset, the drive alone)
    and their ranges over the indicator families and versions, and Q's range over every combination (each claw subset)."""
    rows = []
    keys = ['wiring', 'ratio_id', 'numerator', 'denominator', 'estimate', 'role', 'baseline']
    for key, g in sets.groupby(keys, sort=False):
        apl = g[g['quantity'] == 'APL']
        cl = g[g['quantity'] == 'claws']
        lin = lambda x: float(x[(x['n'] == 1.0) & (x['form'] == 'evoked')]['ratio'].iloc[0])
        row = dict(zip(keys, key), apl_linear=lin(apl), apl_min=float(apl['ratio'].min()), apl_max=float(apl['ratio'].max()))
        for sub in CLAW_SUBSETS:
            c = cl[cl['subset'] == sub]
            row.update({f'claws_{sub}_linear': lin(c[c['share'] == 0.0]), f'claws_{sub}_min': float(c['ratio'].min()),
                        f'claws_{sub}_max': float(c['ratio'].max())})
            q = np.divide.outer(apl['ratio'].to_numpy(), c['ratio'].to_numpy())
            row.update({f'q_{sub}_linear': row['apl_linear'] / row[f'claws_{sub}_linear'],
                        f'q_{sub}_min': float(np.nanmin(q)), f'q_{sub}_max': float(np.nanmax(q))})
        rows.append(row)
    return pd.DataFrame(rows)


def a3_tables(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A3's tables from the saved runs (`load_a3_run`): a3_family.csv (`family_table` on every stimulus, both
    conditions), a3_lines.csv (the line of log2 APL on log2 claws over the 16 panel odors, as A2), a3_x2_sets.csv
    (`model_x2`), a3_x2_summary.csv (`x2_model_summary`), and copies of runs.csv, the stimulus table (every input vector)
    and the identities (a3_runs.csv, a3_stimuli.csv, a3_identity.json)."""
    import shutil
    d = a3_dir()
    fams = []
    for w in wirings:
        run = load_a3_run(w, d)
        if run is None:
            raise FileNotFoundError(f'{w}: no A3 run under {d}')
        fams.append(family_table(w, run, run['meta'], list(H.APL16_ODORS)))
    fam = pd.concat(fams, ignore_index=True)
    out = out_dir()
    fam.to_csv(out / 'a3_family.csv', index=False)
    lines = pd.concat([family_line(fam, x=claw_col(sub, sh, 'change')) for sub in CLAW_SUBSETS for sh in CLAW_SHARES],
                      ignore_index=True)
    lines.to_csv(out / 'a3_lines.csv', index=False)
    sets = model_x2(fam)
    sets.to_csv(out / 'a3_x2_sets.csv.gz', index=False)
    summ = x2_model_summary(sets)
    summ.to_csv(out / 'a3_x2_summary.csv', index=False)
    shutil.copyfile(d / 'runs.csv', out / 'a3_runs.csv')
    shutil.copyfile(d / 'a3_stimuli.csv', out / 'a3_stimuli.csv')
    _write_json(out / 'a3_identity.json', {w: _read_json(d / f'{w}_identity.json') for w in wirings})
    return dict(n_rows=len(fam), n_x2_rows=len(summ))


def a3_apl_tables(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A3's APL side, from the saved runs' APL means and responses only (E98's B1 layout, read by E98's `load_b1`; no
    claw stage): a3_apl.csv (per wiring, condition and stimulus APL's window, pre-odor and settled pre-odor means, its
    evoked activity and the KCs' response rate) and a3_apl_x2.csv (per wiring, X2 ratio and estimate of the input, under
    the three pre-odor levels A0_CHOICES, APL's linear ratio and its range over E98's indicator family)."""
    from experiments import E98_apl_response_strength as E98
    b = E98.load_b1(wirings, roots=[a3_dir()])
    out = out_dir()
    b.to_csv(out / 'a3_apl.csv', index=False)
    rows = []
    for w, g in b[b['cond'] == 'on'].groupby('wiring', sort=False):
        g = g.set_index('stimulus')
        blank = g.loc[BLANK_LABEL] if BLANK_LABEL in g.index else None
        for rid, (nums, dens) in X2_MODEL.items():
            pairs = [(nu, de, 'main') for nu in nums for de in dens]
            pairs += [(nu, dens[0], 'beside') for nu in X2_MODEL_BESIDE.get(rid, ())]
            for nu, de, role in pairs:
                if nu not in g.index or de not in g.index:
                    continue
                for a0 in A0_CHOICES:
                    base = dict(e90=g.loc[nu, 'pre'], settled=g.loc[nu, 'pre_settled'],
                                blank=np.nan if blank is None else blank['win'])[a0]
                    if not np.isfinite(base):
                        continue
                    rs = E98.model_ratio_set(g.loc[nu, 'win'], g.loc[de, 'win'], base, APL_NS)
                    rows.append(dict(wiring=w, ratio_id=rid, numerator=nu, denominator=de, role=role,
                                     estimate=de if rid in ('r1', 'r2') else nu, a0=a0,
                                     apl_linear=float(E98.indicator_ratio(g.loc[nu, 'win'], g.loc[de, 'win'], base, 1.0,
                                                                          'evoked')),
                                     apl_min=float(rs['ratio'].min()), apl_max=float(rs['ratio'].max()),
                                     kc_rate_num=float(g.loc[nu, 'response_rate']),
                                     kc_rate_den=float(g.loc[de, 'response_rate'])))
    x2 = pd.DataFrame(rows)
    x2.to_csv(out / 'a3_apl_x2.csv', index=False)
    return dict(n_rows=len(b), n_x2_rows=len(x2), missing=b.attrs.get('missing'))


# ---------------------------------------------------------------------------
# A2 again at the model's carrying level, C4's model side and the archive's variants: from A3's tables (a3_family.csv,
# a3_runs.csv, a3_identity.json), no simulation and no stage matrix read

HALLEM_ODORS = ('pa', 'ms')
HALLEM_CONCS = (-2.0, -4.0, -6.0, -8.0)  # B1's rows of Hallem's input (olfsysm's hc_data.csv)
LEVELS_SHOWN = (A0_CARRYING, A0_SHOWN)   # A2 again: the carrying level ('blank'), E90's ('e90') beside


def hallem_label(odor: str, conc: float) -> str:
    """B1's label of a stimulus on Hallem's input ('pa [hallem] @ -4')."""
    return f'{odor} [hallem] @ {int(round(conc))}'


def blank_col(subset: str, share: float) -> str:
    """A claw total's column at the blank's level: 'claw_<subset>_<share>_blank'."""
    return claw_col(subset, share, 'blank')


def at_blank(fam: pd.DataFrame, cond: str = 'on') -> pd.DataFrame:
    """One condition's rows of a family table (`family_table`; A3's holds the blank) at the model's carrying level
    (A0_CARRYING, the amendment of 2026-11-15 18:08): per wiring and stimulus, APL's response ('apl_blank': its window
    mean minus the blank's window mean in the same run and condition), each claw total's (`blank_col`: its window minus
    the blank's window, for each subset and claw version), the PN stage's total change minus the blank's ('pn_blank'),
    and the KCs' response rate. A wiring whose run has no blank is refused."""
    rows = []
    for w, g in fam[fam['cond'] == cond].groupby('wiring', sort=False):
        g = g.set_index('stimulus')
        if BLANK_LABEL not in g.index:
            raise ValueError(f'{w} {cond}: no blank ({BLANK_LABEL!r}) in the run')
        b = g.loc[BLANK_LABEL]
        d = pd.DataFrame(dict(apl_win=g['apl_win'], apl_blank=g['apl_win'] - float(b['apl_win'])), index=g.index)
        for sub in CLAW_SUBSETS:
            for sh in CLAW_SHARES:
                c = claw_col(sub, sh, 'win')
                d[blank_col(sub, sh)] = g[c] - float(b[c])
        d['pn_blank'] = g['pn_total'] - float(b['pn_total'])
        d['kc_rate'] = g['kc_rate']
        rows.append(d.reset_index().assign(wiring=w, cond=cond))
    out = pd.concat(rows, ignore_index=True)
    return out[['wiring', 'cond', 'stimulus'] + [c for c in out.columns if c not in ('wiring', 'cond', 'stimulus')]]


def ratio_table(fam: pd.DataFrame, pairs: Sequence[dict], baseline: str = 'blank', cond: str = 'on',
                subset_q: str = 'input') -> pd.DataFrame:
    """Per wiring and pair, in condition `cond` under the model's level `baseline` ('blank' or 'e90'): APL's and each claw
    total's (subset, version) linear ratio and range over the indicator families (`ratio_sets`, `ratio_summary`); Q's,
    APL's ratio over the claws' of `subset_q` (linear: with the drive alone; the range over every combination of the two
    indicator families, their forms and the three claw versions, as `q_classes` reads it); the PN stage's total change
    (minus the blank's under 'blank') and the KCs' response rates (linear). `pairs`: dicts with 'numerator',
    'denominator' and any labels to carry (pair, ratio_id, estimate, role, series, ...); pairs a wiring's run lacks are
    skipped. Long rows: the labels, wiring, cond, baseline, quantity, subset, share, linear, indicator_min, indicator_max,
    kc_rate_num, kc_rate_den."""
    if baseline not in ('blank', 'e90'):
        raise ValueError(f'ratio_table reads the blank\'s level or E90\'s, not {baseline!r}')
    lin = lambda x: float(x[(x['n'] == 1.0) & (x['form'] == 'evoked')]['ratio'].iloc[0])
    rows = []
    for w, g in fam.groupby('wiring', sort=False):
        c = g[g['cond'] == cond].set_index('stimulus')
        for p in pairs:
            nu, de = p['numerator'], p['denominator']
            if nu not in c.index or de not in c.index:
                continue
            key = dict(p, wiring=w, cond=cond, baseline=baseline, kc_rate_num=float(c.loc[nu, 'kc_rate']),
                       kc_rate_den=float(c.loc[de, 'kc_rate']))
            rs = ratio_sets(g, nu, de, (baseline,), cond=cond)
            for r in ratio_summary(rs).to_dict('records'):
                rows.append(dict(key, **{k: v for k, v in r.items() if k != 'baseline'}))
            a = rs[rs['quantity'] == 'APL']
            cl = rs[(rs['quantity'] == 'claws') & (rs['subset'] == subset_q)]
            q = np.divide.outer(a['ratio'].to_numpy(float), cl['ratio'].to_numpy(float)).ravel()
            q = q[np.isfinite(q)]
            rows.append(dict(key, quantity='Q', subset=subset_q, share=np.nan, linear=lin(a) / lin(cl[cl['share'] == 0.0]),
                             indicator_min=float(q.min()) if len(q) else np.nan,
                             indicator_max=float(q.max()) if len(q) else np.nan))
            pn = c['pn_total'] - (float(c.loc[BLANK_LABEL, 'pn_total']) if baseline == 'blank' else 0.0)
            rows.append(dict(key, quantity='PN total', subset='', share=np.nan, linear=float(pn[nu] / pn[de]),
                             indicator_min=np.nan, indicator_max=np.nan))
            rows.append(dict(key, quantity='KC response rate', subset='', share=np.nan,
                             linear=float(c.loc[nu, 'kc_rate'] / c.loc[de, 'kc_rate']) if c.loc[de, 'kc_rate'] > 0
                             else np.nan, indicator_min=np.nan, indicator_max=np.nan))
    return pd.DataFrame(rows)


def ranges_over_wirings(t: pd.DataFrame, by: Sequence[str]) -> List[dict]:
    """Per group of `by` (a `ratio_table`'s labels), the range over the wirings of the linear ratio and of the indicator
    family's ends (the lowest minimum and the highest maximum), with the number of wirings."""
    out = []
    for key, g in t.groupby(list(by), dropna=False, sort=False):
        out.append(dict(zip(by, key), n_wirings=int(g['wiring'].nunique()), linear_min=float(g['linear'].min()),
                        linear_max=float(g['linear'].max()), indicator_min=float(g['indicator_min'].min()),
                        indicator_max=float(g['indicator_max'].max())))
    return out


def part_a2_blank(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """A2 again at the model's carrying level (the amendment of 2026-11-15 18:08: "the note gives A2's line and the three
    ratios again under the carrying level from A3's runs, which hold the 16 panel odors and the blank"), from A3's family
    table (a3_family.csv; no simulation). A3's claw stage equals E90's on the 16 panel odors (part 0), so E90's level
    here reproduces A2. Writes a2_blank_lines.csv (per level in LEVELS_SHOWN, the line of log2 APL on log2 of each claw
    total and of the PN total over the 16 panel odors, `family_line`: slope, r and each odor's residual), a2_blank_ratios.csv
    (pa over ms, paa over ms and paa over pa, `ratio_table`, APL on and acutely off, each condition at its own run's
    blank, with E90's level beside) and a2_blank.json (the ranges over the wirings)."""
    out = out_dir()
    fam = pd.read_csv(out / 'a3_family.csv')
    panel = list(H.APL16_ODORS)
    fam = fam[fam['wiring'].isin(wirings) & fam['stimulus'].isin(panel + [BLANK_LABEL])]
    lines = []
    for level in LEVELS_SHOWN:
        if level == 'blank':
            f = at_blank(fam, 'on')
            xs = [(blank_col(sub, sh), sub, sh) for sub in CLAW_SUBSETS for sh in CLAW_SHARES] + [('pn_blank', 'PN', np.nan)]
            y = 'apl_blank'
        else:
            f = fam
            xs = [(claw_col(sub, sh, 'change'), sub, sh) for sub in CLAW_SUBSETS for sh in CLAW_SHARES] + \
                 [('pn_total', 'PN', np.nan)]
            y = 'apl_evoked'
        for x, sub, sh in xs:
            lines.append(family_line(f, x=x, y=y, odors=panel).assign(level=level, subset=sub, share=sh))
    lines = pd.concat(lines, ignore_index=True)
    lines.to_csv(out / 'a2_blank_lines.csv', index=False)
    pairs = [dict(pair=f'{a}/{b}', numerator=a, denominator=b) for a, b in X1_PAIRS]
    ratios = pd.concat([ratio_table(fam, pairs, level, cond) for level in LEVELS_SHOWN for cond in ('on', 'off')],
                       ignore_index=True)
    ratios.to_csv(out / 'a2_blank_ratios.csv', index=False)
    ln = lines.drop_duplicates(['level', 'subset', 'share', 'wiring'])
    summ = dict(
        note=("A2 at the model's carrying level (the blank's window mean in the same run and condition), from A3's runs "
              "(16 panel odors and the blank); E90's level beside reproduces A2 (a2_lines.csv, a2_ratios.csv)"),
        lines=[dict(level=lv, subset=sub, share=None if pd.isna(sh) else float(sh), n_wirings=int(g['wiring'].nunique()),
                    slope_rma_min=float(g['slope_rma'].min()), slope_rma_max=float(g['slope_rma'].max()),
                    r_min=float(g['r'].min()), r_max=float(g['r'].max()), n_odors_min=int(g['n_odors'].min()))
               for (lv, sub, sh), g in ln.groupby(['level', 'subset', 'share'], dropna=False, sort=False)],
        residual_extremes=[dict(level=lv, odor=o, residual_log2_min=float(g['residual_log2'].min()),
                                residual_log2_max=float(g['residual_log2'].max()))
                           for (lv, o), g in lines[(lines['subset'] == 'input') & (lines['share'] == 0.0)]
                           .groupby(['level', 'odor'], sort=False)],
        ratios=ranges_over_wirings(ratios, ['baseline', 'cond', 'pair', 'quantity', 'subset', 'share']))
    _write_json(out / 'a2_blank.json', summ)
    main = [r for r in summ['lines'] if r['subset'] == 'input' and r['share'] == 0.0]
    return dict(n_lines=len(lines), n_ratios=len(ratios), lines_input_drive=main)


def hallem_pairs() -> List[dict]:
    """C4's pairs on Hallem's input: pa over ms at each concentration (HALLEM_CONCS), and each odor at each lower
    concentration over the same odor at 10^-2."""
    pairs = [dict(series='pa over ms', log10_conc=c, numerator=hallem_label('pa', c), denominator=hallem_label('ms', c))
             for c in HALLEM_CONCS]
    pairs += [dict(series=f'{o} over {o} at 10^-2', log10_conc=c, numerator=hallem_label(o, c),
                   denominator=hallem_label(o, -2.0)) for o in HALLEM_ODORS for c in HALLEM_CONCS[1:]]
    return pairs


def x2_pairs() -> List[dict]:
    """X2's model pairs (X2_MODEL, with X2_MODEL_BESIDE as role 'beside') as `ratio_table` pairs: ratio_id, role and the
    estimate of the input (the denominator's for r1 and r2, the numerator's for r3)."""
    pairs = []
    for rid, (nums, dens) in X2_MODEL.items():
        pairs += [dict(ratio_id=rid, role='main', estimate=de if rid in ('r1', 'r2') else nu, numerator=nu, denominator=de)
                  for nu in nums for de in dens]
        pairs += [dict(ratio_id=rid, role='beside', estimate=nu, numerator=nu, denominator=dens[0])
                  for nu in X2_MODEL_BESIDE.get(rid, ())]
    return pairs


def part_c4_model(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """C4's model side (beside the verdicts, deciding nothing), from A3's family table (no simulation); each condition at
    its own run's blank (the amendment of 2026-11-15 18:08: "C4's readings with APL off take the acute-off run's own
    blank"), E90's level not read:
    - Hallem's input for pa and ms along concentration (HALLEM_CONCS), APL on and acutely off: per stimulus APL's and the
      claws' responses and the KCs' response rate (`at_blank`; c4_hallem_levels.csv, with the panel's pa and ms beside),
      and pa over ms at each concentration and each odor over its 10^-2 for APL, the claws, Q, the PN total and the KCs'
      response rate (`ratio_table`; c4_hallem_ratios.csv);
    - X2's ratios (`x2_pairs`) with APL acutely off, beside the same with APL on (c4_apl_off.csv);
    - c4_model.json, the ranges over the wirings."""
    out = out_dir()
    fam = pd.read_csv(out / 'a3_family.csv')
    fam = fam[fam['wiring'].isin(wirings)]
    hs = [hallem_label(o, c) for o in HALLEM_ODORS for c in HALLEM_CONCS]
    lv = pd.concat([at_blank(fam, cond) for cond in ('on', 'off')], ignore_index=True)
    lv = lv[lv['stimulus'].isin(hs + ['pa', 'ms'])]
    lv.to_csv(out / 'c4_hallem_levels.csv', index=False)
    hr = pd.concat([ratio_table(fam, hallem_pairs(), 'blank', cond) for cond in ('on', 'off')], ignore_index=True)
    hr.to_csv(out / 'c4_hallem_ratios.csv', index=False)
    off = pd.concat([ratio_table(fam, x2_pairs(), 'blank', cond) for cond in ('on', 'off')], ignore_index=True)
    off.to_csv(out / 'c4_apl_off.csv', index=False)
    neg = lv[(lv['apl_blank'] <= 0) | (lv[blank_col('input', 0.0)] <= 0)]
    summ = dict(
        note=("C4's model side, beside and deciding nothing: each condition at its own run's blank (APL acutely off: the "
              "acute-off run's blank); Hallem's input is olfsysm's hc_data.csv (B1's rows)"),
        hallem_levels=[dict(cond=c, stimulus=s, apl_blank_min=float(g['apl_blank'].min()),
                            apl_blank_max=float(g['apl_blank'].max()),
                            claws_input_drive_min=float(g[blank_col('input', 0.0)].min()),
                            claws_input_drive_max=float(g[blank_col('input', 0.0)].max()),
                            kc_rate_min=float(g['kc_rate'].min()), kc_rate_max=float(g['kc_rate'].max()))
                       for (c, s), g in lv.groupby(['cond', 'stimulus'], sort=False)],
        not_positive=sorted({f"{r['wiring']} {r['cond']} {r['stimulus']}" for _, r in neg.iterrows()}),
        hallem_ratios=ranges_over_wirings(hr, ['cond', 'series', 'log10_conc', 'quantity', 'subset', 'share']),
        x2_ratios=ranges_over_wirings(off, ['cond', 'ratio_id', 'role', 'estimate', 'quantity', 'subset', 'share']))
    _write_json(out / 'c4_model.json', summ)
    return dict(n_levels=len(lv), n_hallem_ratios=len(hr), n_x2_ratios=len(off), not_positive=summ['not_positive'])


def flies_strengths(c: pd.DataFrame, groups: Sequence[str] = ('L', 'V0'), multiples: Sequence[float] = (0.0, 1.0),
                    n_boot: int = N_BOOT) -> pd.DataFrame:
    """The flies' calyx responses per stimulus (part C1's rows `c`, the carrying reading), described: per group, b and
    stimulus the numbers of flies and recordings, the mean over flies (each fly's recordings averaged) of APL's response
    and of K_b's (dF/F) with 90 % intervals, and each stimulus over its odor's 10^-2 for APL and for K_b (`flies_ratio`:
    the geometric mean over flies, recordings whose responses are at most MIN_RESPONSE left out)."""
    rows = []
    for grp in groups:
        for bm in multiples:
            x = c[(c['group'] == grp) & (c['reading'] == CARRYING) & (c['b_mult'] == bm)]
            fly = x.drop_duplicates('recording').set_index('recording')['fly']
            for s, h in x.groupby('stim', sort=False):
                odor, conc = str(h['odor'].iloc[0]), float(h['log10_conc'].iloc[0])
                row = dict(group=grp, b_mult=float(bm), stim=s, odor=odor, log10_conc=conc,
                           n_flies=int(h['fly'].nunique()), n_recordings=int(h['recording'].nunique()))
                ref = x[(x['odor'] == odor) & np.isclose(x['log10_conc'], -2.0)].set_index('recording')
                for col in ('apl', 'kc'):
                    m, lo, hi = fly_mean_interval(h.groupby('fly')[col].mean().to_numpy(), n_boot)
                    row.update({f'{col}_mean': m, f'{col}_lo': lo, f'{col}_hi': hi})
                    pr = pd.DataFrame(dict(num=h.set_index('recording')[col], den=ref[col])).dropna()
                    pr['fly'] = fly.reindex(pr.index).to_numpy()
                    r = flies_ratio(pr, n_boot) if len(pr) else dict(ratio=np.nan, lo=np.nan, hi=np.nan, n_flies=0)
                    row.update({f'{col}_over_2': r['ratio'], f'{col}_over_2_lo': r['lo'], f'{col}_over_2_hi': r['hi'],
                                f'{col}_over_2_n_flies': r['n_flies']})
                rows.append(row)
    return pd.DataFrame(rows)


def part_c_strengths(n_boot: int = N_BOOT) -> dict:
    """The flies' strengths per stimulus (`flies_strengths` on c_responses.csv, sets L and V0; c_strengths.csv, with
    Yang's quality caveat first). Described: the APL channel's responses are what outcome (h) reports; the KC channel's
    (K_b) are not shown to be told from the leak."""
    _gate('C')
    out = out_dir()
    t = flies_strengths(read_table(out / 'c_responses.csv'), n_boot=n_boot)
    write_table(t, out / 'c_strengths.csv')
    return dict(n_rows=len(t))


def write_variants(wirings: Sequence[str] = ALL_WIRINGS) -> pd.DataFrame:
    """variants.csv (the archive's params.csv source; one row per wiring's A3 runs, tag = the wiring as in the archived
    file names under A3/): the absolute parameters with APL on (a3_runs.csv; the acute-off run zeroes wAPLKC_scale and
    wAPLPN_scale), the fit it runs at with nothing re-tuned (E83's `members`: E67's or E72's folder and tag), C2's
    breadth exponent, KC-KC time constant and pulse (as E90's variants.csv gives them), the input's sha256, the number of
    stimuli, the commits read at the run, the KCs' response rates over the 16 panel odors with APL on and acutely off,
    and whether the runs reproduce the fit and E98's B1 (a3_identity.json)."""
    from experiments import E83_cross_wiring_fit as E83
    out = out_dir()
    runs = pd.read_csv(out / 'a3_runs.csv')
    fam = pd.read_csv(out / 'a3_family.csv')
    ident = _read_json(out / 'a3_identity.json') or {}
    mem = E83.members()
    rows = []
    for w in wirings:
        r = runs[(runs['wiring'] == w) & (runs['cond'] == 'on')]
        if r.empty:
            continue
        r = r.iloc[-1]
        rate = fam[(fam['wiring'] == w) & fam['stimulus'].isin(list(H.APL16_ODORS))].groupby('cond')['kc_rate'].mean()
        folder, tag = mem[w]
        iw = ident.get(w) or {}
        rows.append(dict(tag=w, variant=w, part='A3', wiring=w, ref_experiment=folder, ref_fit=tag, k=2.0, kc_kc_tau=0.05,
                         time_stim_end=0.5, **{k[len('param_'):]: float(r[k]) for k in r.index if k.startswith('param_')},
                         n_stimuli=int(r['n_stimuli']), input_sha256=r['input_sha256'], queue_commit=r['queue_commit'],
                         pyolfsysm_commit=r['pyolfsysm_commit'], rate_on=float(rate['on']),
                         rate_acute_off=float(rate['off']), reproduces_fit=iw.get('reproduces_fit'),
                         reproduces_b1=iw.get('reproduces_b1')))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    return df


def e90_claw_identity(wirings: Sequence[str] = ALL_WIRINGS, tol: float = 1e-9) -> dict:
    """Part 0 with --e90: A3's claw stage on the 16 panel odors against E90's saved one (`e90_run`), per wiring and
    condition: the largest absolute difference of a claw's window and pre-odor means (relative to the largest value), and
    of the claw totals (`claw_totals`); passed when every difference is at most `tol`."""
    rows = []
    for w in wirings:
        a3 = load_a3_run(w)
        if a3 is None:
            rows.append(dict(wiring=w, ok=False, why='no A3 run'))
            continue
        e90, _ = e90_run(w)
        for _, ct in A3_CONDS:
            d = {}
            for k in ('claw_win', 'claw_pre'):
                x = a3[ct][k][list(H.APL16_ODORS)].to_numpy(dtype=float)
                y = e90[ct][k][list(H.APL16_ODORS)].to_numpy(dtype=float)
                d[k] = float(np.max(np.abs(x - y)) / max(np.max(np.abs(y)), 1e-300)) if x.shape == y.shape else np.inf
            rows.append(dict(wiring=w, cond=ct, claw_win_rel_diff=d['claw_win'], claw_pre_rel_diff=d['claw_pre'],
                             ok=bool(max(d.values()) <= tol)))
    return dict(passed=bool(rows and all(r['ok'] for r in rows)), per_wiring=rows)


# ---------------------------------------------------------------------------
# X3's statistic (part D; A4 on made-up maps)

def _row_z(X: np.ndarray) -> np.ndarray:
    """Rows centred and scaled to unit norm (NaN rows stay NaN), so that Z1 @ Z2.T is Pearson's r between rows."""
    X = np.asarray(X, dtype=float)
    Z = X - X.mean(axis=-1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return Z / np.linalg.norm(Z, axis=-1, keepdims=True)


def odor_specific_parts(maps: np.ndarray, normalize: bool = True) -> np.ndarray:
    """The odor-specific parts of kept stimuli's maps (..., stimuli, cells): each map divided by its mean over cells,
    minus the mean of those maps over the stimuli. A map whose mean over cells is not positive gives NaN. Without
    `normalize` (a channel without an indicator: the controls' KC channel once the leak is removed, whose maps have no
    positive mean), the maps minus their mean over the stimuli."""
    M = np.asarray(maps, dtype=float)
    if not normalize:
        return M - M.mean(axis=-2, keepdims=True)
    m = M.mean(axis=-1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        P = np.where(m > 0, M / np.where(m > 0, m, 1.0), np.nan)
    return P - P.mean(axis=-2, keepdims=True)


def other_odor_assignments(odors: Sequence[str], n_perm: int = N_PERM, seed: int = SEED) -> np.ndarray:
    """X3's null: `n_perm` assignments (n_perm x stimuli) giving each stimulus a stimulus of another odor, drawn uniformly
    and independently per stimulus."""
    od = np.asarray([str(o) for o in odors])
    rng = np.random.default_rng(seed)
    out = np.empty((n_perm, len(od)), dtype=int)
    for s, o in enumerate(od):
        others = np.flatnonzero(od != o)
        if not len(others):
            raise ValueError('X3 needs stimuli of at least two odors')
        out[:, s] = rng.choice(others, size=n_perm)
    return out


def x3_from_parts(oa: Tuple[np.ndarray, np.ndarray], ok: Tuple[np.ndarray, np.ndarray], perms: np.ndarray) -> dict:
    """X3 from the odor-specific parts of the two halves (each (..., stimuli, cells); leading axes are recordings):
    A, the mean over stimuli and the two crossings of the Pearson r between APL's part from one half and K_b's from the
    other; the null (each stimulus given K_b's part of the stimulus `perms` assigns) with its NULL_PCT percentile and mean;
    each channel's repeat r (the mean over stimuli of the r between its two halves) and the ceiling (the square root of
    their product, 0 when either is negative)."""
    za1, za2, zk1, zk2 = (_row_z(x) for x in (oa[0], oa[1], ok[0], ok[1]))
    c12 = np.einsum('...sc,...tc->...st', za1, zk2)
    c21 = np.einsum('...sc,...tc->...st', za2, zk1)
    c = (c12 + c21) / 2.0
    S = c.shape[-1]
    A = np.nanmean(np.diagonal(c, axis1=-2, axis2=-1), axis=-1)
    null = np.nanmean(c[..., np.arange(S)[None, :], perms], axis=-1)            # (..., n_perm)
    ra = np.nanmean(np.einsum('...sc,...sc->...s', za1, za2), axis=-1)
    rk = np.nanmean(np.einsum('...sc,...sc->...s', zk1, zk2), axis=-1)
    ceiling = np.sqrt(np.clip(ra, 0, None) * np.clip(rk, 0, None))
    return dict(A=A, null_pct=np.percentile(null, NULL_PCT, axis=-1), null_mean=null.mean(axis=-1), repeat_r_apl=ra,
                repeat_r_kc=rk, ceiling=ceiling)


def x3_recording(apl_halves: Tuple[np.ndarray, np.ndarray], kc_halves: Tuple[np.ndarray, np.ndarray],
                 odors: Sequence[str], n_perm: int = N_PERM, seed: int = SEED, normalize_kc: bool = True) -> dict:
    """X3 in one recording: the halves' maps of the kept stimuli (stimuli x cells; K_b already formed) made odor-specific
    (`odor_specific_parts`, each half apart; the KC side without the division by its mean when not `normalize_kc`),
    then `x3_from_parts` with `other_odor_assignments`; with the counts."""
    perms = other_odor_assignments(odors, n_perm, seed)
    oa = tuple(odor_specific_parts(h) for h in apl_halves)
    ok = tuple(odor_specific_parts(h, normalize_kc) for h in kc_halves)
    r = {k: float(v) for k, v in x3_from_parts(oa, ok, perms).items()}
    return dict(r, n_stimuli=len(odors), n_odors=len(set(odors)))


def fly_x3(per_recording: pd.DataFrame) -> pd.DataFrame:
    """Per fly and b (`per_recording`: columns fly, b_mult, A, null_pct, repeat_r_apl, repeat_r_kc): the means over its
    recordings and X3's class: 'pass' if A is at least X3_PASS_A and above the null percentile; 'fail' if A is below
    X3_FAIL_A or below the null percentile while both channels' repeat r is at least X3_FAIL_MIN_REPEAT_R; else
    'neither'."""
    g = per_recording.groupby(['fly', 'b_mult'], sort=False)[['A', 'null_pct', 'repeat_r_apl', 'repeat_r_kc']].mean()
    g['n_recordings'] = per_recording.groupby(['fly', 'b_mult'], sort=False).size()
    ok_r = (g['repeat_r_apl'] >= X3_FAIL_MIN_REPEAT_R) & (g['repeat_r_kc'] >= X3_FAIL_MIN_REPEAT_R)
    passed = (g['A'] >= X3_PASS_A) & (g['A'] > g['null_pct'])
    failed = ((g['A'] < X3_FAIL_A) | (g['A'] < g['null_pct'])) & ok_r
    g['cls'] = np.where(passed, 'pass', np.where(failed, 'fail', 'neither'))
    return g.reset_index()


def x3_verdict(fly_classes: pd.DataFrame, cross_check_passed: bool, b_g_alone: bool = False,
               unprobed: Sequence[str] = UNPROBED_X3) -> dict:
    """X3 from one row per fly and b multiple (columns fly, b_mult, cls; the flies holding a recording with X3_MIN_KEPT
    kept stimuli of X3_MIN_ODORS odors): 'follows' if B1's cross-check passed and at least X3_MIN_FLIES flies pass at b =
    0, b_G and 2 b_G (0 and b_G with `b_g_alone`, A4's rule), X3_UNPROBED_MIN of the `unprobed` flies among them; 'follows
    at the controls' leak' if that holds at 0 and b_G but not at 2 b_G (counts as not decided until baseline images
    settle b); 'does not follow' if at least X3_MIN_FLIES flies fail at b = 0, X3_UNPROBED_MIN of the unprobed among them;
    otherwise 'not decided' (also with fewer than X3_MIN_FLIES flies)."""
    flies = sorted(fly_classes['fly'].unique())
    res = dict(n_flies=len(flies), cross_check_passed=bool(cross_check_passed), b_g_alone=bool(b_g_alone))

    def at(b, cls):
        return set(fly_classes[(fly_classes['b_mult'] == b) & (fly_classes['cls'] == cls)]['fly'])

    def holds(s):
        return len(s) >= X3_MIN_FLIES and len(s & set(unprobed)) >= X3_UNPROBED_MIN
    p01 = at(0.0, 'pass') & at(1.0, 'pass')
    p012 = p01 & at(2.0, 'pass')
    res.update(n_pass_0_bg=len(p01), n_pass_0_bg_2bg=len(p012), n_fail_0=len(at(0.0, 'fail')))
    if len(flies) < X3_MIN_FLIES:
        res['verdict'] = 'not decided'
        res['why'] = f'{len(flies)} flies with enough kept stimuli'
    elif cross_check_passed and holds(p01 if b_g_alone else p012):
        res['verdict'] = 'follows'
    elif cross_check_passed and not b_g_alone and holds(p01):
        res['verdict'] = "follows at the controls' leak"
    elif holds(at(0.0, 'fail')):
        res['verdict'] = 'does not follow'
    else:
        res['verdict'] = 'not decided'
    return res


def cross_correlogram(a: np.ndarray, k: np.ndarray, cell_um: float) -> dict:
    """The normalized cross-correlogram of two odor-specific parts on the cell grid (2-D arrays; NaN, outside the outline,
    set to 0 after each part is centred on its inside cells): the peak's lag (cells, rows and columns) and whether it lies
    within one cell of zero lag, the correlogram at zero lag, and the length (um) at which its radial mean falls to 1/e of
    the zero-lag value (NaN if it never does within the grid)."""
    from scipy.signal import correlate
    a, k = np.asarray(a, dtype=float), np.asarray(k, dtype=float)
    m = np.isfinite(a) & np.isfinite(k)
    A = np.where(m, a - np.nanmean(a[m]), 0.0)
    K = np.where(m, k - np.nanmean(k[m]), 0.0)
    den = np.sqrt((A ** 2).sum() * (K ** 2).sum())
    cc = correlate(A, K, mode='full', method='direct') / (den if den > 0 else np.nan)
    c0 = np.array(cc.shape) // 2
    pk = np.unravel_index(np.nanargmax(cc), cc.shape)
    lag = (int(pk[0] - c0[0]), int(pk[1] - c0[1]))
    yy, xx = np.indices(cc.shape)
    rr = np.hypot(yy - c0[0], xx - c0[1])
    zero = float(cc[tuple(c0)])
    length = np.nan
    if zero > 0:
        for r in range(1, int(rr.max()) + 1):
            ring = cc[(rr > r - 0.5) & (rr <= r + 0.5)]
            if ring.size and np.nanmean(ring) <= zero / np.e:
                length = float(r * cell_um)
                break
    return dict(peak_lag_cells=lag, peak_within_one_cell=bool(max(abs(lag[0]), abs(lag[1])) <= 1), zero_lag=zero,
                length_um=length)


def common_map_share(parts: np.ndarray, kc_parts: Optional[np.ndarray] = None) -> dict:
    """E79's measure on normalized maps (stimuli x cells, each map over its mean over cells): the share of APL's sum of
    squares about each map's mean held by one map common to all stimuli (the mean map); with `kc_parts` (the claws'
    normalized maps, other presentations), the share of APL's rest (each map minus the common map) that the claws' rest
    accounts for (one least-squares coefficient over all stimuli and cells)."""
    P = np.asarray(parts, dtype=float)
    P = P - P.mean(axis=1, keepdims=True)
    common = P.mean(axis=0, keepdims=True)
    tot = float((P ** 2).sum())
    out = dict(common_share=float((np.broadcast_to(common, P.shape) ** 2).sum() / tot) if tot > 0 else np.nan)
    if kc_parts is not None:
        Kp = np.asarray(kc_parts, dtype=float)
        Kp = Kp - Kp.mean(axis=1, keepdims=True)
        ra, rk = P - common, Kp - Kp.mean(axis=0, keepdims=True)
        beta = float((ra * rk).sum() / (rk ** 2).sum()) if (rk ** 2).sum() > 0 else np.nan
        rest = float((ra ** 2).sum())
        out.update(claws_share_of_rest=float(1.0 - ((ra - beta * rk) ** 2).sum() / rest) if rest > 0 else np.nan,
                   claws_coefficient=beta)
    return out


# ---------------------------------------------------------------------------
# A4: the place test's power on made-up two-channel maps

def smooth_field(rng: np.random.Generator, n: int, sigma: float, size: Tuple[int, ...] = ()) -> np.ndarray:
    """Random smooth fields on an n x n grid (white noise blurred with a Gaussian of SD `sigma` cells, wrapping at the
    edges), each with mean 0 and SD 1 over cells; `size` leading axes, cells flattened last."""
    from scipy.ndimage import gaussian_filter
    x = rng.standard_normal(tuple(size) + (n, n))
    x = gaussian_filter(x, sigma=(0,) * len(size) + (sigma, sigma), mode='wrap')
    x = x.reshape(tuple(size) + (n * n,))
    x = x - x.mean(axis=-1, keepdims=True)
    return x / x.std(axis=-1, keepdims=True)


def made_up_recordings(rng: np.random.Generator, n_rec: int, r_apl: float, r_kc: float, true_corr: float, leak: float,
                       odors: Sequence[int] = A4_ODORS, n: int = A4_GRID) -> dict:
    """`n_rec` made-up two-channel recordings with the same design: per recording and odor a true APL pattern U and a KC
    pattern V = c U + sqrt(1 - c^2) W (`smooth_field`, c `true_corr`); per kept stimulus (`odors`: stimuli per odor, at
    calyx-mean amplitudes from 0.4 to 1 along each odor's series) four presentations of each channel: APL a (1 +
    A4_MODULATION U + e), the KCs A4_KC_OVER_APL a (1 + A4_MODULATION V + e') plus `leak` times the APL presentation; e and
    e' independent per cell and presentation, their SD set so that each channel's odor-specific parts have repeat r
    `r_apl` / `r_kc` between halves (the mean of two presentations). Returns 'apl' and 'kc' (recordings x presentations x
    stimuli x cells; presentations 1 to 4), 'odors' (per stimulus)."""
    od = np.repeat(np.arange(len(odors)), odors)
    S, C = len(od), n * n
    amp = np.concatenate([np.linspace(0.4, 1.0, k) for k in odors])
    U = smooth_field(rng, n, A4_SMOOTH_CELLS, (n_rec, len(odors)))
    W = smooth_field(rng, n, A4_SMOOTH_CELLS, (n_rec, len(odors)))
    V = true_corr * U + np.sqrt(max(1.0 - true_corr ** 2, 0.0)) * W
    U, V = U[:, od], V[:, od]                                                   # recordings x stimuli x cells

    def sigma(pattern, r):
        clean = odor_specific_parts(1.0 + A4_MODULATION * pattern)
        sig = np.nanmean(clean.var(axis=-1))
        return float(np.sqrt(2.0 * sig * (1.0 - r) / r))                        # a half averages two presentations
    sa, sk = sigma(U, r_apl), sigma(V, r_kc)
    e = rng.standard_normal((n_rec, 4, S, C))
    ek = rng.standard_normal((n_rec, 4, S, C))
    apl = amp[None, None, :, None] * (1.0 + A4_MODULATION * U[:, None] + sa * e)
    kc = A4_KC_OVER_APL * amp[None, None, :, None] * (1.0 + A4_MODULATION * V[:, None] + sk * ek) + leak * apl
    return dict(apl=apl, kc=kc, odors=[str(x) for x in od])


def halves(stack: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """X3's halves from presentations 1 to 4 (axis -3): the mean of the 1st and 3rd, and of the 2nd and 4th."""
    return (stack[..., [0, 2], :, :].mean(axis=-3), stack[..., [1, 3], :, :].mean(axis=-3))


def a4_condition(rng: np.random.Generator, r_apl: float, r_kc: float, share: float, perms: np.ndarray,
                 n_flies: int = A4_N_FLIES, leak: float = A4_B_G) -> List[dict]:
    """One A4 condition: `n_flies` made-up flies of A4_RECORDINGS recordings, the planted leak `leak`; at each b of
    B_MULTIPLES x `leak`, each fly's class by X3's rule (`fly_x3`); rows of the share of flies passing and failing, the
    mean A, ceiling and null percentile, and the noise-free A of the same maps."""
    rec = made_up_recordings(rng, n_flies * A4_RECORDINGS, r_apl, r_kc, share, leak)
    fly = np.repeat(np.arange(n_flies), A4_RECORDINGS)
    ha, hk = halves(rec['apl']), halves(rec['kc'])
    rows = []
    for bm in B_MULTIPLES:
        b = bm * leak
        oa = tuple(odor_specific_parts(h) for h in ha)
        ok = tuple(odor_specific_parts(k - b * a) for k, a in zip(hk, ha))
        r = x3_from_parts(oa, ok, perms)
        pr = pd.DataFrame(dict(fly=fly, b_mult=bm, A=r['A'], null_pct=r['null_pct'], repeat_r_apl=r['repeat_r_apl'],
                               repeat_r_kc=r['repeat_r_kc'], ceiling=r['ceiling']))
        fx = fly_x3(pr)
        rows.append(dict(r_apl=r_apl, r_kc=r_kc, true_share=share, b_mult=bm, power_pass=float((fx['cls'] == 'pass').mean()),
                         share_fail=float((fx['cls'] == 'fail').mean()), mean_A=float(pr['A'].mean()),
                         mean_ceiling=float(pr['ceiling'].mean()), mean_null_pct=float(pr['null_pct'].mean()),
                         mean_repeat_r_apl=float(pr['repeat_r_apl'].mean()), mean_repeat_r_kc=float(pr['repeat_r_kc'].mean())))
    return rows


def a4_rule(power: pd.DataFrame, min_power: float = A4_POWER_MIN) -> dict:
    """A4's rule, decided before any map is read: the share of made-up flies (pooled over the whole grid of reliabilities,
    every pair equally weighted) passing at 2 b_G with a true agreement of half the ceiling; X3's "follows" is read at b_G
    alone (b = 0 and b_G) if it is below `min_power`. The same pooled shares at b = 0 and b_G and the other true agreements
    beside."""
    pooled = power.groupby(['true_share', 'b_mult'])['power_pass'].mean()
    p = float(pooled.loc[(0.5, 2.0)])
    return dict(rule=("If fewer than 0.8 of the made-up flies pass at 2 b_G with a true agreement of half the ceiling, X3's "
                      "'follows' is read at b_G alone (b = 0 and b_G); 'the made-up flies' are pooled over the whole grid of "
                      "repeat reliabilities (0.3 to 0.9 per channel, every pair equally weighted)."),
                power_2bg_half_ceiling=p, min_power=float(min_power), b_g_alone=bool(p < min_power),
                pooled_power={f'share {s:g}, b {b:g} b_G': float(v) for (s, b), v in pooled.items()},
                written=time.strftime('%Y-%m-%d %H:%M'), host=_host(),
                note="decided on made-up maps only, before any of the dual-colour recordings' maps is read")


def part_a4(n_flies: int = A4_N_FLIES, n_perm: int = A4_N_PERM, seed: int = SEED) -> dict:
    """A4: the share of made-up flies that pass (and fail) X3's rule at b = 0, b_G and 2 b_G, per pair of channel
    reliabilities (A4_RELIABILITIES) and true agreement (A4_TRUE_SHARES of the ceiling), with the planted leak A4_B_G
    (a4_power.csv), and the rule on X3's "follows" (`a4_rule`; a4_rule.json, not overwritten once written unless the
    settings changed)."""
    rng = np.random.default_rng(seed)
    od = [str(x) for x in np.repeat(np.arange(len(A4_ODORS)), A4_ODORS)]
    perms = other_odor_assignments(od, n_perm, seed)
    rows = []
    for ra in A4_RELIABILITIES:
        for rk in A4_RELIABILITIES:
            for sh in A4_TRUE_SHARES:
                rows += a4_condition(rng, ra, rk, sh, perms, n_flies)
    power = pd.DataFrame(rows)
    out = out_dir()
    power.to_csv(out / 'a4_power.csv', index=False)
    rule = a4_rule(power)
    rule['settings'] = dict(A4_RELIABILITIES=A4_RELIABILITIES, A4_TRUE_SHARES=A4_TRUE_SHARES, A4_B_G=A4_B_G,
                            A4_KC_OVER_APL=A4_KC_OVER_APL, A4_MODULATION=A4_MODULATION, A4_GRID=A4_GRID,
                            A4_ODORS=A4_ODORS, A4_RECORDINGS=A4_RECORDINGS, n_flies=n_flies, n_perm=n_perm, seed=seed,
                            A4_SMOOTH_CELLS=A4_SMOOTH_CELLS)
    path = out / 'a4_rule.json'
    old = _read_json(path)
    if old is not None and old.get('settings') == json.loads(json.dumps(rule['settings'], default=float)):
        rule = old
    else:
        _write_json(path, rule)
    return dict(rule=rule, n_conditions=len(power))


# ---------------------------------------------------------------------------
# the readers (parts B to E; viki's copy, both channels)

def released(stem: str) -> bool:
    """Whether Yang has released a recording for E101 (YANG_RELEASED, and RELEASED_STEMS when she names some), or it is
    one of the 12 2021 recordings she released before (data/yang_apl_calyx_maps/)."""
    from experiments import E100_apl_weak_odors as E100
    if stem in E100.released_2021():
        return True
    return bool(YANG_RELEASED and (not RELEASED_STEMS or stem in RELEASED_STEMS))


def check_readable(row) -> None:
    """The two refusals: a raw-only row (rule of 2026-11-15 15:28, recordings.csv's use_until_further_notice) and a
    recording Yang has not released (`released`). How the channels are read is gated per part (`_gate`, CHANNEL_MODE)."""
    from experiments import E100_apl_weak_odors as E100
    if not E100.usable_row(row):
        raise PermissionError(f"{row['stem']}: raw only, not used until further notice (rule of 2026-11-15 15:28, "
                              "recordings.csv's use_until_further_notice)")
    if not released(str(row['stem'])):
        raise PermissionError(f"{row['stem']}: Yang has not released this recording for E101 (question 13); set "
                              "YANG_RELEASED (and RELEASED_STEMS) after her answer")


def copied_variants(row) -> List[str]:
    """The peak-map variants viki copied for a recording (recordings.csv's peak_map_copied, 'peak5_ave_<tag>; ...', a
    trailing remark in parentheses dropped), as tags ('C1_rig_moco', 'C2_rig_moco_C1shift', 'C1', 'moco')."""
    import re
    v = row.get('peak_map_copied')
    if not isinstance(v, str):
        return []
    v = re.sub(r'\s*\(.*\)\s*$', '', v)
    return [x.strip()[len('peak5_ave_'):] for x in v.split(';') if x.strip().startswith('peak5_ave_')]


def matching_c1_tag(c2_tag: str) -> str:
    """The C1 registration that a C2 variant was made with: 'C2_rig_moco_C1shift' -> 'C1_rig_moco', 'C2' -> 'C1',
    'C2_moco' -> 'C1_moco'."""
    t = 'C1' + c2_tag[2:]
    return t[:-len('_C1shift')] if t.endswith('_C1shift') else t


def channel_paths(row) -> dict:
    """One recording's files in viki's copy, per channel, from the variants copied (`copied_variants`): C2's peak-map
    stack (the one moved by C1's shifts where copied) and C1's with the matching registration (`matching_c1_tag`; else
    the copied C1 variant, with 'registration_matched' False), Yang's ROI table of each, the outline matching C1's tag
    ('<stem>_<tag>.roi', else '<stem>.roi'), the yaml and the Experiment.xml. A single-channel recording's one stack is C1.
    The 12 released 2021 recordings are read from data/yang_apl_calyx_maps/ (C1 only, E100's `recording_paths`). Missing
    files are None."""
    from experiments import E100_apl_weak_odors as E100
    stem = str(row['stem'])
    if stem in E100.released_2021():
        f = E100.recording_paths(row)
        return dict(yaml=f['yaml'], roi=f['roi'], experiment_xml=f['experiment_xml'], C1=f['peak_maps'], C1_tag=f['tag'],
                    C1_table=f['roi_table'], C2=None, C2_tag='', C2_table=None, registration_matched=True)
    folder = Path(row['folder']) if isinstance(row.get('folder'), str) and row['folder'] else Path('/nonexistent')
    tags = copied_variants(row)
    c2 = sorted((t for t in tags if t.startswith('C2')), key=lambda t: ('C1shift' not in t, t))
    c1 = [t for t in tags if not t.startswith('C2')]
    c2_tag = c2[0] if c2 else ''
    c1_tag = ''
    if c1:
        want = matching_c1_tag(c2_tag) if c2_tag else None
        c1_tag = want if want in c1 else c1[0]

    def stack(tag):
        f = folder / f'{stem}_{tag}_dff_peak5frames_ave.tif'
        return f if tag and f.exists() else None

    def table(tag):
        f = folder / f'{stem}_{tag}_roi_peak_dff.csv'
        return f if tag and f.exists() else None
    roi = next((f for f in ([folder / f'{stem}_{c1_tag}.roi'] if c1_tag else []) + [folder / f'{stem}.roi'] if f.exists()),
               None)
    yml = next((f for f in (folder / f'{stem}.yaml', folder / f'{stem}_copy.yaml') if f.exists()), None)
    xml, xml_from = experiment_xml_path(folder, stem)
    return dict(yaml=yml, roi=roi, experiment_xml=xml, experiment_xml_from=xml_from, C1=stack(c1_tag), C1_tag=c1_tag,
                C1_table=table(c1_tag),
                C2=stack(c2_tag), C2_tag=c2_tag, C2_table=table(c2_tag),
                registration_matched=bool(not c2_tag or not c1_tag or matching_c1_tag(c2_tag) == c1_tag))


def experiment_xml_path(folder: Path, stem: str) -> Tuple[Optional[Path], str]:
    """A recording's Experiment.xml in its folder ('<stem>*Experiment.xml'), and the stem it was found under. Where there is
    none, the same date, fly, side and run with the other sex letter (Yang's processed files of 2022-10-05 fly 1, a male,
    are named '1MR...' and its raw run folders '1FR...'; `read_two_channel` checks the frame size against the maps)."""
    xml = next(iter(sorted(Path(folder).glob(f'{stem}*Experiment.xml'))), None)
    if xml is not None:
        return xml, stem
    p = _stem_parts(stem)
    for sx in ('F', 'M', ''):
        if sx == p['sex']:
            continue
        alt = f"{p['date']}_{p['fly']}{sx}{p['side']}fun{p['run']}"
        cand = sorted(Path(folder).glob(f'{alt}_Experiment.xml'))
        if cand:
            return cand[0], alt
    return None, ''


def read_two_channel(row, maps: bool = True) -> dict:
    """One recording read for parts B to E (after `check_readable`): its presentations (E100's `yaml_stimuli`; frames past
    the last complete block dropped where the stack is shorter, as 20221213_2FLfun2), each channel's stack (one shape),
    the outline mask, each presentation's calyx response per channel ('C1', 'C2') with Yang's own value beside where her
    table has one row per presentation ('yang_C1', 'yang_C2'), the pixel size and channels from the Experiment.xml."""
    import tifffile
    from experiments import E100_apl_weak_odors as E100
    check_readable(row)
    _, DS = E100._e79()
    f = channel_paths(row)
    if (f['C1'] is None and f['C2'] is None) or f['yaml'] is None or f['roi'] is None:
        raise FileNotFoundError(f"{row['stem']}: peak maps, yaml or outline missing ({f})")
    stacks = {c: tifffile.imread(f[c]).astype(float) for c in CHANNELS if f[c] is not None}
    shapes = {c: s.shape for c, s in stacks.items()}
    if len({s[1:] for s in shapes.values()}) != 1:
        raise ValueError(f"{row['stem']}: the channels' maps differ in shape: {shapes}")
    if len({s[0] for s in shapes.values()}) != 1:
        raise ValueError(f"{row['stem']}: the channels have different numbers of maps: {shapes}")
    frames = E100.yaml_stimuli(f['yaml'], row['stem'])
    n = next(iter(stacks.values())).shape[0]
    nstim = int(frames['position'].max())
    truncated = False
    if n != len(frames):
        if n > len(frames) or n % nstim:
            raise ValueError(f"{row['stem']}: {n} peak maps for {len(frames)} presentations")
        frames = frames.iloc[:n].copy()                 # read to its last complete block
        truncated = True
    mask = DS.ACM.roi_mask(DS.ACM.read_imagej_polygon(f['roi']), next(iter(stacks.values())).shape[1:])
    frames['is_solvent'] = frames['odor'] == E100.SOLVENT
    frames['stim'] = [E100.SOLVENT if sv else E100.stim_label((o, c)) for sv, o, c in
                      zip(frames['is_solvent'], frames['odor'], frames['log10_conc'])]
    for c, st in stacks.items():
        frames[c] = E100.calyx_responses(st, mask)
        tab = f.get(f'{c}_table')
        if tab is not None:
            t = pd.read_csv(tab)
            if len(t) >= len(frames):
                frames[f'yang_{c}'] = t['roi_peak_dff'].to_numpy(dtype=float)[:len(frames)]
    xml = E100.experiment_xml(f['experiment_xml']) if f['experiment_xml'] is not None else None
    if xml is not None and f.get('experiment_xml_from', row['stem']) != row['stem']:
        shape = next(iter(stacks.values())).shape[1:]
        if (xml['pixel_y'], xml['pixel_x']) != tuple(shape):
            raise ValueError(f"{row['stem']}: the Experiment.xml of {f['experiment_xml_from']} has frames "
                             f"{xml['pixel_y']} x {xml['pixel_x']}, the maps {shape}")
    out = dict(stem=row['stem'], fly=fly_of(row['stem']), calyx=calyx_of(row['stem']), group=group_of(row['stem']),
               plane_um=row.get('plane_depth_um'), frames=frames, mask=mask, files=f, xml=xml, truncated=truncated,
               channels=sorted(stacks), registration_matched=f['registration_matched'],
               experiment_xml_from=f.get('experiment_xml_from', ''),
               pixel_um=float(xml['pixel_um']) if xml is not None else float(row.get('pixel_size_um', np.nan)))
    if not np.isfinite(out['pixel_um']):
        raise ValueError(f"{row['stem']}: no pixel size (no Experiment.xml, none in recordings.csv)")
    if maps:
        out['stacks'] = stacks
    return out


def cell_grid(mask: np.ndarray, pixel_um: float, cell_um: float = CELL_UM, min_inside: float = CELL_MIN_INSIDE) -> dict:
    """The cells of `cell_um` (rounded to whole pixels) over the outline's bounding box: each pixel's cell index (-1 outside
    the box or in a cell left out), the cells kept (at least `min_inside` of their pixels inside the outline) as
    (row, column) on the cell grid, the grid's shape, the cell size in pixels and um."""
    mask = np.asarray(mask, dtype=bool)
    k = max(int(round(float(cell_um) / float(pixel_um))), 1)
    rows, cols = np.flatnonzero(mask.any(axis=1)), np.flatnonzero(mask.any(axis=0))
    r0, c0 = int(rows[0]), int(cols[0])
    nr, nc = (int(rows[-1]) - r0) // k + 1, (int(cols[-1]) - c0) // k + 1
    rr, cc = np.indices(mask.shape)
    gr, gc = (rr - r0) // k, (cc - c0) // k
    in_box = (rr >= r0) & (cc >= c0) & (gr < nr) & (gc < nc)
    cell = np.where(in_box, gr * nc + gc, -1)
    inside = np.bincount(cell[in_box & mask], minlength=nr * nc)
    total = np.bincount(cell[in_box], minlength=nr * nc)
    keep = (total > 0) & (inside >= min_inside * total)
    idx = np.full(nr * nc, -1)
    idx[keep] = np.arange(int(keep.sum()))
    pix = np.where(in_box & mask, idx[np.clip(cell, 0, None)], -1)
    pix[~(in_box & mask)] = -1
    kept = np.flatnonzero(keep)
    return dict(pixel_cell=pix, kept_rc=np.stack([kept // nc, kept % nc], axis=1), shape=(nr, nc), cell_px=k,
                cell_um=k * float(pixel_um), n_cells=int(keep.sum()))


def cell_maps(stack: np.ndarray, grid: dict) -> np.ndarray:
    """Per presentation, the mean of each kept cell's inside pixels (presentations x kept cells)."""
    st = np.asarray(stack, dtype=float)
    pc = grid['pixel_cell']
    m = pc >= 0
    lab = pc[m]
    cnt = np.bincount(lab, minlength=grid['n_cells']).astype(float)
    out = np.stack([np.bincount(lab, weights=fr[m], minlength=grid['n_cells']) for fr in st]) / cnt
    return out


def cells_to_grid(values: np.ndarray, grid: dict) -> np.ndarray:
    """Kept cells' values back on the cell grid (NaN elsewhere), for the cross-correlogram."""
    g = np.full(grid['shape'], np.nan)
    g[grid['kept_rc'][:, 0], grid['kept_rc'][:, 1]] = values
    return g


# ---------------------------------------------------------------------------
# the statistics (parts B to E)

def first_in_block_solvent(frames: pd.DataFrame) -> np.ndarray:
    """The presentations (positions in `frames`) that are a solvent first in its block."""
    f = frames.reset_index(drop=True)
    return np.flatnonzero(f['is_solvent'].to_numpy() & (f['position'].to_numpy() == 1))


def stimulus_responses(frames: pd.DataFrame, values: np.ndarray, reading: str = CARRYING) -> pd.DataFrame:
    """Per stimulus, the response under `reading` (READINGS) from per-presentation values (`values`: presentations, or
    presentations x cells, in `frames`' order): 'carrying' the mean over its presentations minus the mean of the
    recording's first-in-block solvent presentations; 'block' the mean over its presentations of each minus its block's
    first solvent; 'first_two' the mean of its first two presentations minus the recording's first-in-block solvent. Index
    the stimulus labels (no solvent); columns 'value' (or the cells, for maps) and 'n'."""
    f = frames.reset_index(drop=True)
    v = np.asarray(values, dtype=float)
    v2 = v.reshape(len(f), -1)
    fib = first_in_block_solvent(f)
    if not len(fib):
        raise ValueError('no first-in-block solvent presentation')
    base = v2[fib].mean(axis=0)
    blk = {int(f.loc[i, 'block']): v2[i] for i in fib}
    rows, idx, ns = [], [], []
    for s, g in f[~f['is_solvent']].groupby('stim', sort=False):
        ii = g.index.to_numpy()
        if reading == 'carrying':
            x = v2[ii].mean(axis=0) - base
        elif reading == 'block':
            x = np.mean([v2[i] - blk.get(int(f.loc[i, 'block']), base) for i in ii], axis=0)
        elif reading == 'first_two':
            ii = ii[:2]
            x = v2[ii].mean(axis=0) - base
        else:
            raise ValueError(reading)
        rows.append(x)
        idx.append(s)
        ns.append(len(ii))
    out = pd.DataFrame(np.vstack(rows), index=pd.Index(idx, name='stim'),
                       columns=['value'] if v.ndim == 1 else list(range(v2.shape[1])))
    out['n'] = ns
    return out


def solvent_responses(frames: pd.DataFrame, values: np.ndarray) -> np.ndarray:
    """Every solvent presentation's value minus the mean of the first-in-block solvent presentations (the spread that the
    kept-stimulus rule compares with)."""
    f = frames.reset_index(drop=True)
    v = np.asarray(values, dtype=float)
    return v[f['is_solvent'].to_numpy()] - v[first_in_block_solvent(f)].mean()


def leak_slopes(apl: np.ndarray, kc: np.ndarray) -> np.ndarray:
    """Per cell, the least-squares slope (with an intercept) of the KC channel on the APL channel across stimuli
    (stimuli x cells each)."""
    a, k = np.asarray(apl, dtype=float), np.asarray(kc, dtype=float)
    ac, kc_ = a - a.mean(axis=0), k - k.mean(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (ac * kc_).sum(axis=0) / (ac ** 2).sum(axis=0)


def remove_leak(kc, apl, b: float):
    """K_b = C2 - b C1 (presentation by presentation, maps or calyx means)."""
    return np.asarray(kc, dtype=float) - float(b) * np.asarray(apl, dtype=float)


def kept_stimuli(resp_apl: pd.Series, resp_kc: pd.Series, solvent_apl: np.ndarray, solvent_kc: np.ndarray,
                 mult: float = KEPT_SD_MULT) -> List[str]:
    """The stimuli whose APL response and KC response (at b_G) each exceed `mult` times the SD of that channel's solvent
    responses (`solvent_responses`)."""
    ta = mult * float(np.std(np.asarray(solvent_apl, dtype=float), ddof=1))
    tk = mult * float(np.std(np.asarray(solvent_kc, dtype=float), ddof=1))
    return [s for s in resp_apl.index if resp_apl[s] > ta and s in resp_kc.index and resp_kc[s] > tk]


def _odor(stim: str) -> str:
    return str(stim).split(' @ ')[0]


def x1_recording(apl: pd.Series, kc: pd.Series, pairs: Sequence[Tuple[str, str]] = X1_PAIRS,
                 min_kept: int = X1_MIN_KEPT_PER_ODOR) -> dict:
    """X1 in one recording from its kept stimuli (stimulus label -> response, both positive): the RMA line of log2 APL on
    log2 K_b (`rma_line`), each odor's mean residual, and d for each pair (the first odor's mean residual minus the
    second's) where the recording has `min_kept` or more kept stimuli of each odor (else NaN)."""
    s = [x for x in apl.index if x in kc.index and apl[x] > 0 and kc[x] > 0]
    x, y = np.log2(kc[s].to_numpy(dtype=float)), np.log2(apl[s].to_numpy(dtype=float))
    ln = rma_line(x, y)
    res = pd.Series(ln['residuals'], index=s)
    od = pd.Series([_odor(v) for v in s], index=s)
    mean_res = res.groupby(od).mean()
    count = od.value_counts()
    d = {}
    for a, b in pairs:
        ok = count.get(a, 0) >= min_kept and count.get(b, 0) >= min_kept
        d[f'{a}/{b}'] = float(mean_res[a] - mean_res[b]) if ok else np.nan
    return dict(slope=ln['slope'], intercept=ln['intercept'], r=ln['r'], n_kept=len(s), mean_residual=mean_res.to_dict(),
                d=d)


def fly_mean_interval(per_fly: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED) -> Tuple[float, float, float]:
    """The mean over flies with its 90 % interval from `n_boot` bootstrap draws of flies (NaN for no fly)."""
    v = np.asarray(per_fly, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = v[rng.integers(0, len(v), size=(n_boot, len(v)))].mean(axis=1)
    lo, hi = np.percentile(boot, CI)
    return float(v.mean()), float(lo), float(hi)


def x1_verdict(intervals: Dict[str, Tuple[float, float]], n_flies: Dict[str, int], bound: float = X1_BOUND,
               min_flies: int = X1_MIN_FLIES) -> str:
    """X1 from each pair's interval of d (all three pairs of X1_PAIRS): 'one function' if every interval lies within
    +-`bound`; 'depends on the odor' if some pair's interval lies wholly beyond it; otherwise 'not decided'; 'not decided'
    also when a pair has fewer than `min_flies` flies."""
    pairs = [f'{a}/{b}' for a, b in X1_PAIRS]
    if any(int(n_flies.get(p, 0)) < min_flies or p not in intervals for p in pairs):
        return 'not decided'
    iv = [intervals[p] for p in pairs]
    if all(-bound <= lo and hi <= bound for lo, hi in iv):
        return 'one function'
    if any(lo > bound or hi < -bound for lo, hi in iv):
        return 'depends on the odor'
    return 'not decided'


def flies_ratio(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """A ratio over flies from per-recording numerator and denominator (columns fly, num, den): recordings whose
    numerator or denominator is at most MIN_RESPONSE left out, each fly's recordings averaged on the log scale, the
    geometric mean over flies with its 90 % interval."""
    d = per_recording[(per_recording['den'] > MIN_RESPONSE) & (per_recording['num'] > MIN_RESPONSE)]
    per_fly = np.log(d['num'] / d['den']).groupby(d['fly']).mean()
    m, lo, hi = fly_mean_interval(per_fly.to_numpy(), n_boot, seed)
    return dict(ratio=float(np.exp(m)), lo=float(np.exp(lo)), hi=float(np.exp(hi)), n_flies=int(len(per_fly)),
                n_recordings=int(len(d)), n_left_out=int(len(per_recording) - len(d)))


def flies_q(per_recording: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Q, APL's ratio over the claws' ratio within the fly (columns fly, apl_num, apl_den, kc_num, kc_den; a recording
    enters only if all four exceed MIN_RESPONSE), each fly's recordings averaged on the log scale, the geometric mean over
    flies with its 90 % interval."""
    d = per_recording[(per_recording[['apl_num', 'apl_den', 'kc_num', 'kc_den']] > MIN_RESPONSE).all(axis=1)]
    lq = np.log(d['apl_num'] / d['apl_den']) - np.log(d['kc_num'] / d['kc_den'])
    per_fly = lq.groupby(d['fly']).mean()
    m, lo, hi = fly_mean_interval(per_fly.to_numpy(), n_boot, seed)
    return dict(ratio=float(np.exp(m)), lo=float(np.exp(lo)), hi=float(np.exp(hi)), n_flies=int(len(per_fly)),
                n_recordings=int(len(d)), n_left_out=int(len(per_recording) - len(d)))


def flies_at(per_stim: pd.DataFrame, stim: Stim, col: str) -> pd.Series:
    """Per recording, the response at `stim` (`per_stim`: columns recording, odor, log10_conc and `col`): imaged, or ms
    read between its imaged neighbours on log concentration (E98's `value_at`, the recording's own range)."""
    from experiments import E98_apl_response_strength as E98
    out = {}
    for rec, g in per_stim[per_stim['odor'] == stim[0]].groupby('recording'):
        v, _ = E98.value_at(g['log10_conc'].to_numpy(float), g[col].to_numpy(float), stim[1])
        out[rec] = v
    return pd.Series(out, dtype=float)


def x2_quantity_verdict(ratio_classes: Dict[str, str], min_ratios: int = X2_MIN_RATIOS) -> str:
    """A quantity's X2 class from its three ratios' classes: 'flatter', 'inside' or 'steeper' if at least `min_ratios` of
    the three are and the third is not outside on the other side ('flatter' with a 'steeper' third, or 'steeper' with a
    'flatter' third, is not decided; 'inside' takes any third); otherwise 'not decided'."""
    v = list(ratio_classes.values())
    for k, other in (('flatter', 'steeper'), ('inside', None), ('steeper', 'flatter')):
        if sum(x == k for x in v) >= min_ratios and (other is None or other not in v):
            return k
    return 'not decided'


def q_classes(sets: pd.DataFrame, lo: float, hi: float, baseline: str = 'e90', subset: str = 'input') -> pd.DataFrame:
    """Per wiring and estimate, Q's class against the flies' interval [lo, hi] over every combination of the two indicator
    families, their forms and the claw versions (`sets`: `model_x2`'s rows of one ratio; Q = an APL row's ratio over a
    claw row's), E98's `classify`; with Q's range and its linear value."""
    from experiments import E98_apl_response_strength as E98
    rows = []
    s = sets[(sets['baseline'] == baseline)]
    for (w, est), g in s.groupby(['wiring', 'estimate'], sort=False):
        a = g[g['quantity'] == 'APL']
        c = g[(g['quantity'] == 'claws') & (g['subset'] == subset)]
        q = np.divide.outer(a['ratio'].to_numpy(float), c['ratio'].to_numpy(float)).ravel()
        q = q[np.isfinite(q) & (q > 0)]
        lin = float(a[(a['n'] == 1.0) & (a['form'] == 'evoked')]['ratio'].iloc[0]
                    / c[(c['n'] == 1.0) & (c['form'] == 'evoked') & (c['share'] == 0.0)]['ratio'].iloc[0])
        rows.append(dict(wiring=w, estimate=est, q_linear=lin, q_min=float(q.min()), q_max=float(q.max()),
                         cls=E98.classify(q, lo, hi)))
    return pd.DataFrame(rows)


def held_class(classes: Dict[object, str], unprobed: Optional[str], n_unprobed: int, min_unprobed: int) -> str:
    """A class held over the readings (`classes`: one class per b in (0, b_G) and reading in READINGS, keyed as the
    caller likes): the class if every reading gives it, it is decided, and the unprobed flies give it by the same rule
    (`unprobed`, None to skip) with at least `min_unprobed` of them remaining; otherwise 'not decided'."""
    v = set(classes.values())
    if len(v) != 1:
        return 'not decided'
    c = v.pop()
    if c == 'not decided':
        return c
    if unprobed is not None and (unprobed != c or int(n_unprobed) < min_unprobed):
        return 'not decided'
    return c


def _waiting(part: str) -> None:
    raise NotImplementedError(
        f"part {part} reads Yang's two-channel recordings from viki's copy only after her release (question 13) and with "
        "the channels either confirmed by her (question 2) or read as inferred and checked by the controls (Tom, "
        "2026-11-15 17:33): YANG_RELEASED and CHANNEL_MODE; its readers and statistics are written and tested on made-up "
        "maps")


class ChannelCheckNotPassed(RuntimeError):
    """Parts C to E and the score refuse to run while the channels are read as inferred and part B's channel check has not
    passed (channel_check.json)."""


def reader_identity_released(stem: Optional[str] = None) -> dict:
    """Part 0: on one released 2021 recording (default the first of E100's `released_2021`), the two-channel reader's
    calyx means of C1 (`read_two_channel`) equal E100's reader's (`read_recording`) and its presentations are E100's;
    the difference from Yang's own table is reported beside (E100 found 2 to 8 %; read as her answer to question 10
    says), not required."""
    from experiments import E100_apl_weak_odors as E100
    recs = E100.read_recordings(sets=None)
    stem = E100.released_2021()[0] if stem is None else stem
    row = recs[recs['stem'] == stem].iloc[0]
    a = read_two_channel(row, maps=False)['frames']
    b = E100.read_recording(row)['frames']
    d = float(np.max(np.abs(a['C1'].to_numpy() - b['response'].to_numpy())))
    dy = float(np.max(np.abs(a['C1'].to_numpy() - a['yang_C1'].to_numpy()))) if 'yang_C1' in a else np.nan
    ry = float(np.median(a['C1'].to_numpy() / a['yang_C1'].to_numpy())) if 'yang_C1' in a else np.nan
    same = bool((a['stim'].to_numpy() == b['stim'].to_numpy()).all())
    return dict(passed=bool(d <= 1e-12 and same), recording=stem, max_abs_diff_e100=d, presentations_equal=same,
                max_abs_diff_yang_table=dy, median_ratio_to_yang_table=ry)


BLANK_SETTLED_TOL = 0.005               # the blank's APL window mean within this share of APL's settled level
PRE_SAME_TOL = 1e-9                     # a pre-odor level the same for every stimulus, relative to its size (at least 1)


def a3_blank_checks(label: str, run: Dict[str, dict]) -> List[dict]:
    """Part 0's checks that replace its blank item (the amendment of 2026-11-15 18:08), on one A3 run (`load_a3_run`), APL
    on and acutely off: the blank's APL window mean lies within BLANK_SETTLED_TOL of APL's settled level; APL's settled
    level and E90's pre-odor mean of each claw total (the claws of the input's glomeruli and all claws; E90's versions
    scale the drive's pre-odor total by a constant, so the drive's total stands for them) are the same for every stimulus
    of the run, to PRE_SAME_TOL; the KCs' response rate to the blank is 0. One row per condition."""
    gl = np.asarray(run['meta']['claw_glomerulus']).astype(str)
    masks = {'input': np.isin(gl, [str(g) for g in H.apl16_orn_deltas().index]), 'all': np.ones(len(gl), dtype=bool)}
    rows = []
    for _, ct in A3_CONDS:
        o = run[ct]
        a = o['apl']
        if BLANK_LABEL not in a.index:
            rows.append(dict(wiring=label, cond=ct, passed=False, why='no blank in the run'))
            continue
        st = a['pre_settled'].astype(float)
        b = a.loc[BLANK_LABEL]
        rel = float(abs(b['win'] - b['pre_settled']) / abs(b['pre_settled']))
        rng = float(st.max() - st.min())
        row = dict(wiring=label, cond=ct, blank_apl_win=float(b['win']), blank_apl_settled=float(b['pre_settled']),
                   blank_apl_rel_diff=rel, blank_apl_ok=bool(rel <= BLANK_SETTLED_TOL), apl_settled_range=rng,
                   apl_settled_same=bool(rng <= PRE_SAME_TOL * max(1.0, abs(float(st.mean())))))
        P = o['claw_pre'].astype(float).to_numpy()
        for sub, m in masks.items():
            tot = P[m].sum(axis=0)
            r_ = float(tot.max() - tot.min())
            row[f'claw_pre_{sub}_range'] = r_
            row[f'claw_pre_{sub}_same'] = bool(r_ <= PRE_SAME_TOL * max(1.0, abs(float(tot.mean()))))
        rate = float(o['responses'][BLANK_LABEL].astype(float).mean())
        row.update(blank_response_rate=rate, blank_rate_zero=bool(rate == 0.0))
        row['passed'] = bool(all(row[k] for k in ('blank_apl_ok', 'apl_settled_same', 'claw_pre_input_same',
                                                  'claw_pre_all_same', 'blank_rate_zero')))
        rows.append(row)
    return rows


def part_0(e90: bool = False) -> dict:
    """Part 0's identities that need no new recording (identity.json): E100's reader on one 2021 recording (E100's
    `reader_identity_2021`, E79's calyx means), the two-channel reader against E100's on a released recording
    (`reader_identity_released`), A3's run identities (each wiring's <label>_identity.json: the saved fit's counts and
    responses, E98's B1 where saved, the blank) and, with `e90`, A3's claw stage against E90's saved one
    (`e90_claw_identity`)."""
    from experiments import E100_apl_weak_odors as E100
    res = dict(reader_2021=E100.reader_identity_2021(), reader_released=reader_identity_released())
    a3 = {w: _read_json(a3_dir() / f'{w}_identity.json') for w in ALL_WIRINGS}
    a3 = {w: v for w, v in a3.items() if v is not None}
    blank = {w: dict(apl_evoked_e90_pre=v['on'].get('blank_apl_evoked_e90_pre'),
                     apl_evoked_settled_pre=v['on'].get('blank_apl_evoked_settled_pre'),
                     claw_change_input_e90_pre=v['on'].get('blank_claw_change_input'),
                     response_rate_on=v['on'].get('blank_response_rate')) for w, v in a3.items()}
    res['a3_runs'] = dict(passed=bool(a3 and all(v.get('reproduces_fit') and v.get('reproduces_b1', True) is not False
                                                 for v in a3.values())),
                          wirings=sorted(a3), b1_compared=sorted(w for w, v in a3.items() if 'reproduces_b1' in v),
                          blank=blank,
                          note=("the blank's evoked APL and claw change against E90's pre-odor window hold the KC layer's "
                                "start-up (E98's part 0); against the settled pre-odor level APL's is within 0.02 of 0"))
    blank_rows = []
    for w in sorted(a3):
        run = load_a3_run(w)
        blank_rows += a3_blank_checks(w, run) if run is not None else [dict(wiring=w, passed=False, why='no A3 run')]
    res['a3_blank'] = dict(
        passed=bool(blank_rows and all(r['passed'] for r in blank_rows)), rows=blank_rows,
        item_as_written=("'A blank evokes nothing at the model's claws or in its APL': failed as written (against E90's "
                         "pre-odor window, which holds the KC layer's start-up), its stop lifted by the amendment of "
                         "2026-11-15 18:08; these checks take its place"),
        not_checked="whether the blank's claw window lies at the claws' settled level (A3 saved none)")
    if e90:
        res['e90_claws'] = e90_claw_identity(sorted(a3))
    res['passed'] = all(v['passed'] for v in res.values() if isinstance(v, dict) and v.get('passed') is not None)
    path = out_dir() / 'identity.json'
    old = _read_json(path) or {}
    old.update(res)
    _write_json(path, old)
    return res


# ---------------------------------------------------------------------------
# parts B to E and the score (gated: they read Yang's two-channel recordings only after her release)

A0_CHOICES = ('e90', 'settled', 'blank')     # the model's pre-odor levels for X2 (A0_CARRYING, A0_PARTNER, A0_SHOWN)
B_SCAN = tuple(np.round(np.arange(0.0, 4.0 + 1e-9, 0.25), 2))   # D: multiples of b_G scanned for each fly's largest pass


def _gate(part: str) -> None:
    """Parts B to E and the score: nothing is read before Yang's release or while CHANNEL_MODE is 'waiting'; with the
    channels read as inferred (Tom, 2026-11-15 17:33), parts C to E and the score need part B's channel check passed."""
    if not YANG_RELEASED or CHANNEL_MODE not in CHANNEL_MODES[1:]:
        _waiting(part)
    if part != 'B' and CHANNEL_MODE == INFERRED_CHECKED:
        chk = _read_json(out_dir() / 'channel_check.json')
        if chk is None or chk.get('verdict') != 'passed':
            raise ChannelCheckNotPassed(
                f"part {part}: the channels are read as inferred (C1 = APL, C2 = KCs; {CHANNEL_DECISION}) and part B's "
                f"channel check has {'not run' if chk is None else 'given ' + repr(chk.get('verdict'))} "
                "(channel_check.json); E101 stops after part B")


# ---- the channel check (part B, first)

def channel_responses(frames: pd.DataFrame) -> dict:
    """A two-channel recording's calyx responses per stimulus under the carrying reading (`stimulus_responses`) and the SD
    of each channel's solvent responses (`solvent_responses`), per channel."""
    out = {}
    for c in CHANNELS:
        v = frames[c].to_numpy(dtype=float)
        out[c] = stimulus_responses(frames, v)['value']
        out[f'{c}_solvent_sd'] = float(np.std(solvent_responses(frames, v), ddof=1))
    return out


def _responds(resp: pd.Series, sd: float) -> Tuple[str, float, float, bool]:
    """The strongest stimulus of a channel, its response, the threshold, and whether it responds (above
    CHANNEL_RESPOND_SD SDs of its solvent responses and CHANNEL_RESPOND_MIN)."""
    s = str(resp.idxmax())
    v = float(resp[s])
    thr = max(CHANNEL_RESPOND_SD * float(sd), CHANNEL_RESPOND_MIN)
    return s, v, thr, bool(np.isfinite(v) and v > thr)


def channel_class(resp: dict, expected: str, other: str) -> dict:
    """One control recording's class (`resp` from `channel_responses`): 'as inferred' if `expected` responds to its
    strongest stimulus (`_responds`) and `other`'s response to that stimulus is at most CHANNEL_OTHER_MAX of it in absolute
    value; 'swapped' if the same holds with the two channels' roles exchanged; otherwise 'neither' (also if both hold).
    With the numbers that decide it."""
    row = {}
    ok = {}
    for role, a, b in (('expected', expected, other), ('swapped', other, expected)):
        s, v, thr, r = _responds(resp[a], resp[f'{a}_solvent_sd'])
        frac = float(abs(resp[b][s]) / v) if v > 0 else np.nan
        ok[role] = bool(r and np.isfinite(frac) and frac <= CHANNEL_OTHER_MAX)
        row.update({f'{role}_channel': a, f'{role}_strongest': s, f'{role}_response': v, f'{role}_threshold': thr,
                    f'{role}_responds': r, f'{role}_solvent_sd': float(resp[f'{a}_solvent_sd']),
                    f'{role}_other_response': float(resp[b][s]), f'{role}_other_fraction': frac})
    row['cls'] = ('neither' if ok['expected'] and ok['swapped'] else 'as inferred' if ok['expected'] else
                  'swapped' if ok['swapped'] else 'neither')
    return row


def channel_verdict(rows: pd.DataFrame) -> dict:
    """The check from its per-recording rows (columns group, fly, cls; groups G and R; recordings without both channels
    absent): each fly 'as inferred' or 'swapped' if at least CHANNEL_FLY_MIN_SHARE of its recordings are, else 'neither';
    'passed' if every fly is 'as inferred' with at least one fly in each of G and R, 'swapped' if every fly is 'swapped'
    (the same condition), otherwise 'failed'."""
    flies = []
    for (g, f), x in rows.groupby(['group', 'fly'], sort=True):
        n = len(x)
        si, ss = float((x['cls'] == 'as inferred').mean()), float((x['cls'] == 'swapped').mean())
        flies.append(dict(group=g, fly=f, n_recordings=n, share_as_inferred=si, share_swapped=ss,
                          cls='as inferred' if si >= CHANNEL_FLY_MIN_SHARE - 1e-12 else
                          'swapped' if ss >= CHANNEL_FLY_MIN_SHARE - 1e-12 else 'neither'))
    fl = pd.DataFrame(flies, columns=['group', 'fly', 'n_recordings', 'share_as_inferred', 'share_swapped', 'cls'])
    both = all((fl['group'] == g).any() for g in CHANNEL_EXPECTED)
    if both and (fl['cls'] == 'as inferred').all():
        v = 'passed'
    elif both and (fl['cls'] == 'swapped').all():
        v = 'swapped'
    else:
        v = 'failed'
    why = '' if v == 'passed' else (
        'a control group has no fly with a two-channel recording' if not both else
        '; '.join(f"{r['fly']} ({r['group']}): {r['cls']}" for r in flies if r['cls'] != 'as inferred'))
    return dict(verdict=v, why=why, flies=flies)


def channel_check(recs: Sequence[dict]) -> Tuple[pd.DataFrame, dict]:
    """The channel check over the control recordings `recs` (read by `read_two_channel`, groups G and R): per recording
    `channel_class` (a recording without both channels, or whose responses cannot be formed, is listed with
    cls 'not evaluable' and left out of the verdict), then `channel_verdict`, with the thresholds and Tom's decision."""
    rows = []
    for r in recs:
        g = r.get('group')
        if g not in CHANNEL_EXPECTED:
            continue
        base = dict(group=g, recording=r['stem'], fly=r.get('fly', fly_of(r['stem'])))
        if 'skipped' in r or not set(CHANNELS) <= set(r['frames'].columns):
            rows.append(dict(base, cls='not evaluable', why=r.get('skipped', 'one channel only')))
            continue
        try:
            resp = channel_responses(r['frames'])
        except ValueError as e:
            rows.append(dict(base, cls='not evaluable', why=str(e)))
            continue
        rows.append(dict(base, **channel_class(resp, *CHANNEL_EXPECTED[g]), why=''))
    t = pd.DataFrame(rows)
    ev = t[t['cls'] != 'not evaluable'] if len(t) else t
    res = channel_verdict(ev if len(ev) else pd.DataFrame(columns=['group', 'fly', 'cls']))
    res.update(mode=CHANNEL_MODE, decision=CHANNEL_DECISION, yang_question_2_confirmed=YANG_CHANNELS_CONFIRMED,
               thresholds=dict(CHANNEL_EXPECTED={k: list(v) for k, v in CHANNEL_EXPECTED.items()},
                               CHANNEL_RESPOND_SD=CHANNEL_RESPOND_SD, CHANNEL_RESPOND_MIN=CHANNEL_RESPOND_MIN,
                               CHANNEL_OTHER_MAX=CHANNEL_OTHER_MAX, CHANNEL_FLY_MIN_SHARE=CHANNEL_FLY_MIN_SHARE,
                               reading=CARRYING),
               n_recordings=int(len(ev)), n_not_evaluable=int(len(t) - len(ev)),
               written=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    return t, res


def group_rows(groups: Sequence[str]) -> pd.DataFrame:
    """viki's recordings.csv rows of the recordings of `groups` (E100's `read_recordings`, every set), with 'group'."""
    from experiments import E100_apl_weak_odors as E100
    recs = E100.read_recordings(sets=None)
    stems = {s: g for g in groups for s in GROUPS[g]}
    out = recs[recs['stem'].isin(stems)].copy()
    out['group'] = out['stem'].map(stems)
    return out.reset_index(drop=True)


def read_groups(groups: Sequence[str], maps: bool = True) -> List[dict]:
    """Every readable recording of `groups` (`read_two_channel`); a recording that cannot be read is listed in the
    returned dicts' 'skipped' entries (with the reason) and left out."""
    out = []
    for row in group_rows(groups).to_dict('records'):
        try:
            rec = read_two_channel(row, maps=maps)
        except (FileNotFoundError, ValueError) as e:
            out.append(dict(stem=row['stem'], group=row['group'], skipped=str(e)))
            continue
        rec['group'] = row['group']
        out.append(rec)
    return out


def two_channel(recs: Sequence[dict]) -> Tuple[List[dict], List[dict]]:
    """The recordings read with both channels, and the others ({stem, group, why}): skipped by the reader, or one channel
    only."""
    ok, left = [], []
    for r in recs:
        if 'skipped' in r:
            left.append(dict(stem=r['stem'], group=r.get('group'), why=r['skipped']))
        elif not set(CHANNELS) <= set(r['frames'].columns):
            left.append(dict(stem=r['stem'], group=r.get('group'), why='one channel only'))
        else:
            ok.append(r)
    return ok, left


def recording_cells(rec: dict, cell_um: float = CELL_UM) -> dict:
    """A recording's cell grid (`cell_grid` with its own pixel size) and each channel's cell maps (presentations x
    cells)."""
    grid = cell_grid(rec['mask'], rec['pixel_um'], cell_um)
    return dict(grid=grid, **{c: cell_maps(rec['stacks'][c], grid) for c in rec['stacks']})


def _baselines(frames: pd.DataFrame, values: np.ndarray, reading: str) -> np.ndarray:
    """Each presentation's solvent baseline under `reading` (presentations x ...): the mean of the recording's
    first-in-block solvent presentations, or for 'block' its block's own first solvent (else the recording's)."""
    f = frames.reset_index(drop=True)
    v = np.asarray(values, dtype=float)
    fib = first_in_block_solvent(f)
    base = v[fib].mean(axis=0)
    if reading != 'block':
        return np.broadcast_to(base, v.shape)
    blk = {int(f.loc[i, 'block']): v[i] for i in fib}
    return np.stack([blk.get(int(b), base) for b in f['block']])


def stim_halves(frames: pd.DataFrame, values: np.ndarray, reading: str = CARRYING) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Per stimulus, X3's two halves of its presentations, each presentation minus its baseline (`_baselines`): the mean of
    its 1st and 3rd and of its 2nd and 4th presentations ('first_two': its 1st and its 2nd); stimuli with fewer
    presentations are left out."""
    f = frames.reset_index(drop=True)
    v = np.asarray(values, dtype=float) - _baselines(f, values, reading)
    out = {}
    for s, g in f[~f['is_solvent']].groupby('stim', sort=False):
        ii = g.sort_values('frame').index.to_numpy()
        if reading == 'first_two':
            if len(ii) >= 2:
                out[s] = (v[ii[0]], v[ii[1]])
        elif len(ii) >= 4:
            out[s] = (v[ii[[0, 2]]].mean(axis=0), v[ii[[1, 3]]].mean(axis=0))
    return out


def recording_kept(rec: dict, b_g: float, apl_only: bool = False) -> List[str]:
    """A recording's kept stimuli under the carrying reading: APL's calyx response and K_b's at `b_g` each above
    KEPT_SD_MULT SDs of that channel's solvent responses (`kept_stimuli`); `apl_only` (the green-only and marker controls,
    whose KC channel carries no indicator) applies the APL side only."""
    f = rec['frames']
    a = stimulus_responses(f, f['C1'].to_numpy())['value']
    k_v = remove_leak(f['C2'].to_numpy(), f['C1'].to_numpy(), b_g)
    k = stimulus_responses(f, k_v)['value']
    sa, sk = solvent_responses(f, f['C1'].to_numpy()), solvent_responses(f, k_v)
    if apl_only:
        sk = np.zeros_like(sk)
        k = pd.Series(np.inf, index=k.index)
    return kept_stimuli(a, k, sa, sk)


def x3_rows(rec: dict, cells: dict, b_g: float, kept: Sequence[str], multiples: Sequence[float] = B_MULTIPLES,
            readings: Sequence[str] = READINGS, n_perm: int = N_PERM, normalize_kc: bool = True) -> List[dict]:
    """X3 in one recording (`x3_recording`) for each reading and multiple of `b_g` on its `kept` stimuli; rows with the
    recording's fly, group and plane. A reading without X3_MIN_KEPT kept stimuli of X3_MIN_ODORS odors gives no row."""
    f = rec['frames']
    rows = []
    for rd in readings:
        ha = stim_halves(f, cells['C1'], rd)
        hk = stim_halves(f, cells['C2'], rd)
        st = [s for s in kept if s in ha and s in hk]
        od = [_odor(s) for s in st]
        if len(st) < X3_MIN_KEPT or len(set(od)) < X3_MIN_ODORS:
            continue
        A = tuple(np.stack([ha[s][h] for s in st]) for h in (0, 1))
        K = tuple(np.stack([hk[s][h] for s in st]) for h in (0, 1))
        for bm in multiples:
            x = x3_recording(A, tuple(remove_leak(k, a, bm * b_g) for k, a in zip(K, A)), od, n_perm=n_perm,
                             seed=SEED + zlib_seed(rec['stem']), normalize_kc=normalize_kc)
            rows.append(dict(recording=rec['stem'], fly=rec['fly'], group=rec.get('group'), plane_um=rec.get('plane_um'),
                             registration_matched=rec.get('registration_matched'), reading=rd, b_mult=float(bm),
                             b=float(bm * b_g), **x))
    return rows


def zlib_seed(text: str) -> int:
    import zlib
    return int(zlib.crc32(str(text).encode()) % 100000)


# ---- part B

def b1_leak(recs: Sequence[dict], cells: Dict[str, dict], n_perm: int = N_PERM) -> dict:
    """B1 on group G: per recording the least-squares slope of the KC channel on the APL channel across stimuli per cell
    (stimulus maps under the carrying reading; `leak_slopes`); b_G per fly and pooled (the median over cells and
    recordings); the cross-check: each green-only fly's KC channel corrected with the other fly's b, its X3 agreement A
    against its own null's percentile (APL-kept stimuli, `x3_rows`; the KC channel without an indicator is not divided by
    its mean, `odor_specific_parts`); passed when every fly's A is below it."""
    slopes, per_rec = {}, []
    for r in recs:
        c = cells[r['stem']]
        a = stimulus_responses(r['frames'], c['C1']).drop(columns='n').to_numpy()
        k = stimulus_responses(r['frames'], c['C2']).drop(columns='n').to_numpy()
        sl = leak_slopes(a, k)
        slopes.setdefault(r['fly'], []).append(sl)
        per_rec.append(dict(recording=r['stem'], fly=r['fly'], median_slope=float(np.nanmedian(sl)),
                            n_cells=int(np.isfinite(sl).sum())))
    per_fly = {f: float(np.nanmedian(np.concatenate(v))) for f, v in slopes.items()}
    b_g = float(np.nanmedian(np.concatenate([x for v in slopes.values() for x in v])))
    cross = []
    flies = sorted(per_fly)
    for f in flies:
        others = [per_fly[o] for o in flies if o != f]
        if not others:
            continue
        b_other = float(np.median(others))
        for r in (x for x in recs if x['fly'] == f):
            kept = recording_kept(r, 0.0, apl_only=True)
            for row in x3_rows(r, cells[r['stem']], 1.0, kept, multiples=(0.0, b_other), readings=(CARRYING,),
                               n_perm=n_perm, normalize_kc=False):
                cross.append(dict(row, b_other_fly=b_other))
    cr = pd.DataFrame(cross)
    passed = None
    if len(cr):
        corr = cr[cr['b_mult'] > 0].groupby('fly')[['A', 'null_pct']].mean()
        passed = bool(len(corr) == len(flies) and (corr['A'] < corr['null_pct']).all())
    return dict(b_G=b_g, per_fly=per_fly, per_recording=per_rec, cross_check=cr, cross_check_passed=passed)


def part_b(n_perm: int = N_PERM) -> dict:
    """Part B (the controls, first of all recordings): first the channel check (`channel_check` on groups G and R;
    channel_check.csv, channel_check.json; its verdict gates parts C to E while the channels are read as inferred), then
    B1 (`b1_leak`; b_leak.json, b1_slopes.csv, b1_cross_check.csv);
    B2, group R's APL channel where APL carries no indicator, per stimulus and as a share of set L's median APL response
    at the same stimulus, and the per-cell slope of that channel on the KC channel (b2_reverse.csv); B3, group M's red
    marker's response and its agreement with APL's map (b3_marker.csv); B4, set L's solvent maps: the r over cells between
    one channel's mean solvent map from the odd solvent presentations and the other's from the even ones (b4_solvent.csv);
    B5, set L's own slope against b_G (b5_slopes.csv; outcome (h) if its median is below B5_SLOPE_MULT b_G)."""
    _gate('B')
    out = out_dir()
    G_all = read_groups(('G',))
    R_all = read_groups(('R',))
    chk_rows, chk = channel_check(G_all + R_all)
    chk['caveat'] = QUALITY_CAVEAT
    write_table(chk_rows, out / 'channel_check.csv')
    _write_json(out / 'channel_check.json', chk)
    print(f"channel check: {chk['verdict']}{' (' + chk['why'] + ')' if chk['why'] else ''}", flush=True)
    G = [r for r in G_all if 'skipped' not in r and set(CHANNELS) <= set(r['stacks'])]
    cells = {r['stem']: recording_cells(r) for r in G}
    b1 = b1_leak(G, cells, n_perm)
    write_table(pd.DataFrame(b1['per_recording']), out / 'b1_slopes.csv')
    write_table(b1['cross_check'], out / 'b1_cross_check.csv')
    b_g = b1['b_G']
    L, l_left = two_channel(read_groups(('L',)))
    lc = {r['stem']: recording_cells(r) for r in L}
    l_apl = pd.concat([stimulus_responses(r['frames'], r['frames']['C1'].to_numpy())['value'].rename(r['stem'])
                       for r in L], axis=1).median(axis=1)
    rows = []
    for r in R_all:
        if 'skipped' in r or 'C1' not in r.get('stacks', {}) or 'C2' not in r.get('stacks', {}):
            rows.append(dict(recording=r['stem'], skipped=r.get('skipped', 'one channel only')))
            continue
        c = recording_cells(r)
        a = stimulus_responses(r['frames'], r['frames']['C1'].to_numpy())['value']
        k = stimulus_responses(r['frames'], r['frames']['C2'].to_numpy())['value']
        sl = leak_slopes(stimulus_responses(r['frames'], c['C2']).drop(columns='n').to_numpy(),
                         stimulus_responses(r['frames'], c['C1']).drop(columns='n').to_numpy())
        for s in a.index:
            o, cc = s.split(' @ ')
            lm = l_apl[[x for x in l_apl.index if x.split(' @ ')[0] == o]]
            near = min(lm.index, key=lambda x: abs(float(x.split(' @ ')[1]) - float(cc))) if len(lm) else None
            ok = near is not None and abs(float(near.split(' @ ')[1]) - float(cc)) <= MATCH_LOG_UNITS + 1e-9
            rows.append(dict(recording=r['stem'], fly=r['fly'], stim=s, apl_channel=float(a[s]), kc_channel=float(k[s]),
                             set_l_apl=float(lm[near]) if ok else np.nan, set_l_stim=near if ok else '',
                             share_of_set_l=float(a[s] / lm[near]) if ok else np.nan,
                             median_slope_apl_on_kc=float(np.nanmedian(sl))))
    write_table(pd.DataFrame(rows), out / 'b2_reverse.csv')
    rows = []
    for r in read_groups(('M',)):
        if 'skipped' in r:
            rows.append(dict(recording=r['stem'], skipped=r['skipped']))
            continue
        c = recording_cells(r)
        k = stimulus_responses(r['frames'], r['frames']['C2'].to_numpy())['value']
        kept = recording_kept(r, 0.0, apl_only=True)
        x = x3_rows(r, c, 1.0, kept, multiples=(0.0,), readings=(CARRYING,), n_perm=n_perm, normalize_kc=False)
        rows.append(dict(recording=r['stem'], registration_matched=r['registration_matched'],
                         marker_response_median=float(k.median()), marker_response_max_abs=float(k.abs().max()),
                         **({k_: v for k_, v in x[0].items() if k_ in ('A', 'null_pct', 'ceiling', 'repeat_r_apl',
                                                                         'repeat_r_kc')} if x else {})))
    write_table(pd.DataFrame(rows), out / 'b3_marker.csv')
    rows, b5 = [], []
    for r in L:
        f = r['frames'].reset_index(drop=True)
        sv = np.flatnonzero(f['is_solvent'].to_numpy())
        c = lc[r['stem']]
        for x, y in (('C1', 'C2'), ('C2', 'C1')):
            mx, my = c[x][sv[0::2]].mean(axis=0), c[y][sv[1::2]].mean(axis=0)
            rows.append(dict(recording=r['stem'], fly=r['fly'], first=x, second=y,
                             r_cells=float(np.corrcoef(mx, my)[0, 1])))
        a = stimulus_responses(f, c['C1']).drop(columns='n').to_numpy()
        k = stimulus_responses(f, c['C2']).drop(columns='n').to_numpy()
        b5.append(dict(recording=r['stem'], fly=r['fly'], median_slope=float(np.nanmedian(leak_slopes(a, k)))))
    write_table(pd.DataFrame(rows), out / 'b4_solvent.csv')
    b5 = pd.DataFrame(b5)
    write_table(b5, out / 'b5_slopes.csv')
    set_l = float(b5.groupby('fly')['median_slope'].median().median()) if len(b5) else np.nan
    res = dict(caveat=QUALITY_CAVEAT, channel_check=chk['verdict'], b_G=b_g, b_G_per_fly=b1['per_fly'],
               cross_check_passed=b1['cross_check_passed'],
               set_l_median_slope=set_l,
               # outcome (h), "the controls fail (B1's cross-check, or B5)": either clause; B5's alone beside
               outcome_h_b5=bool(np.isfinite(set_l) and set_l < B5_SLOPE_MULT * b_g),
               outcome_h=bool((not b1['cross_check_passed'])
                              or (np.isfinite(set_l) and set_l < B5_SLOPE_MULT * b_g)),
               skipped=[dict(stem=r['stem'], why=r['skipped']) for r in G_all + R_all if 'skipped' in r] + l_left)
    _write_json(out / 'b_leak.json', res)
    return res


# ---- part C

def c_rows(rec: dict, b_g: float, multiples: Sequence[float] = B_MULTIPLES,
           readings: Sequence[str] = READINGS) -> pd.DataFrame:
    """Part C1 for one recording: per reading, multiple of b_G and stimulus, APL's calyx response and K_b's, the odor and
    concentration, and whether the stimulus is kept (`recording_kept`, the carrying reading at b_G)."""
    f = rec['frames']
    kept = set(recording_kept(rec, b_g))
    rows = []
    for rd in readings:
        a = stimulus_responses(f, f['C1'].to_numpy(), rd)['value']
        for bm in multiples:
            k = stimulus_responses(f, remove_leak(f['C2'].to_numpy(), f['C1'].to_numpy(), bm * b_g), rd)['value']
            for s in a.index:
                o, c = s.split(' @ ')
                rows.append(dict(recording=rec['stem'], fly=rec['fly'], group=rec.get('group'), plane_um=rec.get('plane_um'),
                                 reading=rd, b_mult=float(bm), stim=s, odor=o, log10_conc=float(c), apl=float(a[s]),
                                 kc=float(k[s]), kept=s in kept))
    return pd.DataFrame(rows)


def x1_table(c: pd.DataFrame, flies: Optional[Sequence[str]] = None, min_flies: int = X1_MIN_FLIES,
             n_boot: int = N_BOOT) -> Tuple[pd.DataFrame, Dict[Tuple[str, float], dict]]:
    """X1 per reading and b (`c`: `c_rows` of set L): per recording `x1_recording` on its kept stimuli, each fly's d the
    mean over its recordings, the flies' mean with its interval per pair and `x1_verdict` (`flies` restricts the flies;
    `min_flies` the verdict's minimum). Returns the per-recording rows and {(reading, b_mult): result}."""
    c = c if flies is None else c[c['fly'].isin(flies)]
    rows, res = [], {}
    for (rd, bm), g in c.groupby(['reading', 'b_mult'], sort=False):
        for rec, h in g[g['kept']].groupby('recording', sort=False):
            x = x1_recording(h.set_index('stim')['apl'], h.set_index('stim')['kc'])
            rows.append(dict(reading=rd, b_mult=bm, recording=rec, fly=h['fly'].iloc[0], slope=x['slope'], r=x['r'],
                             n_kept=x['n_kept'], **{f'd_{p}': v for p, v in x['d'].items()}))
        per = pd.DataFrame([r for r in rows if r['reading'] == rd and r['b_mult'] == bm])
        iv, nf, summ = {}, {}, {}
        for a, b in X1_PAIRS:
            p = f'{a}/{b}'
            v = per.groupby('fly')[f'd_{p}'].mean().dropna() if len(per) else pd.Series(dtype=float)
            m, lo, hi = fly_mean_interval(v.to_numpy(), n_boot)
            iv[p], nf[p] = (lo, hi), int(len(v))
            summ[p] = dict(mean=m, lo=lo, hi=hi, n_flies=int(len(v)))
        res[(rd, float(bm))] = dict(pairs=summ, verdict=x1_verdict(iv, nf, min_flies=min_flies))
    return pd.DataFrame(rows), res


def x2_flies(c: pd.DataFrame, flies: Optional[Sequence[str]] = None, n_boot: int = N_BOOT) -> pd.DataFrame:
    """X2's flies' side per reading, b and ratio (X2_RATIOS): for the claws (K_b), APL and Q, per numerator reading the
    flies' ratio with its interval (`flies_ratio`, `flies_q`; ms at 10^-3 read between its neighbours, `flies_at`), and the
    hull of the numerator readings' intervals (lo, hi)."""
    from experiments import E98_apl_response_strength as E98
    c = c if flies is None else c[c['fly'].isin(flies)]
    rows = []
    for (rd, bm), g in c.groupby(['reading', 'b_mult'], sort=False):
        fly = g.drop_duplicates('recording').set_index('recording')['fly']
        for rid, nums, den in X2_RATIOS:
            per_q = {}
            for num in nums:
                vals = {col: (flies_at(g, num, col), flies_at(g, den, col)) for col in ('apl', 'kc')}
                for q, col in (('APL', 'apl'), ('claws', 'kc')):
                    n_, d_ = vals[col]
                    pr = pd.DataFrame(dict(num=n_, den=d_)).dropna()
                    pr['fly'] = fly.reindex(pr.index).to_numpy()
                    per_q.setdefault(q, []).append(dict(num=num, **flies_ratio(pr, n_boot)))
                pq = pd.DataFrame(dict(apl_num=vals['apl'][0], apl_den=vals['apl'][1], kc_num=vals['kc'][0],
                                       kc_den=vals['kc'][1])).dropna()
                pq['fly'] = fly.reindex(pq.index).to_numpy()
                per_q.setdefault('Q', []).append(dict(num=num, **flies_q(pq, n_boot)))
            for q, lst in per_q.items():
                ok = [x for x in lst if np.isfinite(x['lo']) and np.isfinite(x['hi'])]
                lo, hi = E98.hull([(x['lo'], x['hi']) for x in ok]) if ok else (np.nan, np.nan)
                rows.append(dict(reading=rd, b_mult=bm, ratio_id=rid, quantity=q, lo=lo, hi=hi,
                                 n_flies=min((x['n_flies'] for x in lst), default=0),
                                 readings='; '.join(f"{stim_label_(x['num'])}: {x['ratio']:.3g} [{x['lo']:.3g}, "
                                                    f"{x['hi']:.3g}] ({x['n_flies']} flies)" for x in lst)))
    return pd.DataFrame(rows)


def stim_label_(s: Stim) -> str:
    return f'{s[0]} @ {float(s[1]):g}'


def part_c(n_boot: int = N_BOOT) -> dict:
    """Part C: c_responses.csv (`c_rows` for sets L, L+, V0 and VT at b = 0, b_G and 2 b_G under the three readings; b_G
    from part B's b_leak.json), c_x1_recordings.csv and c_x1.json (`x1_table` on set L, and on its unprobed flies), and
    c_x2_flies.csv (`x2_flies` on set L and its unprobed flies; V0, L+ and VT beside). The model's side and the classes are
    formed in `score`."""
    _gate('C')
    out = out_dir()
    bl = _read_json(out / 'b_leak.json')
    if bl is None:
        raise RuntimeError('part C needs part B first (b_leak.json)')
    b_g = float(bl['b_G'])
    recs, left = two_channel(read_groups(('L', 'L+', 'V0', 'VT'), maps=False))
    c = pd.concat([c_rows(r, b_g) for r in recs], ignore_index=True)
    write_table(c, out / 'c_responses.csv')
    write_table(pd.DataFrame(left, columns=['stem', 'group', 'why']), out / 'c_left_out.csv')
    L = c[c['group'] == 'L']
    x1r, x1 = x1_table(L, n_boot=n_boot)
    _, x1u = x1_table(L, flies=UNPROBED_X12, min_flies=HELD_X12_MIN, n_boot=n_boot)
    write_table(x1r, out / 'c_x1_recordings.csv')
    _write_json(out / 'c_x1.json', dict(caveat=QUALITY_CAVEAT, all={f'{k[0]} b{k[1]:g}': v for k, v in x1.items()},
                                        unprobed={f'{k[0]} b{k[1]:g}': v for k, v in x1u.items()}))
    x2 = pd.concat([x2_flies(L, n_boot=n_boot).assign(flies='L'),
                    x2_flies(L, flies=UNPROBED_X12, n_boot=n_boot).assign(flies='L unprobed')]
                   + [x2_flies(c[c['group'] == g], n_boot=n_boot).assign(flies=g) for g in ('L+', 'V0', 'VT')
                      if (c['group'] == g).any()], ignore_index=True)
    write_table(x2, out / 'c_x2_flies.csv')
    return dict(n_rows=len(c), x1={f'{k[0]} b{k[1]:g}': v['verdict'] for k, v in x1.items()})


# ---- part D

def largest_b_passing(rows: pd.DataFrame) -> float:
    """The largest multiple of b_G at which a fly still passes (its recordings' rows over B_SCAN), NaN if none."""
    f = fly_x3(rows)
    ok = f[f['cls'] == 'pass']['b_mult']
    return float(ok.max()) if len(ok) else np.nan


def d_described(rec: dict, cells: dict, b_g: float, kept: Sequence[str]) -> dict:
    """Part D's descriptions for one recording (carrying reading, b_G): the cross-correlogram of the two channels'
    odor-specific parts (each the mean of its halves; `cross_correlogram` per stimulus, then the mean length and the share
    of stimuli whose peak lies within one cell of zero lag), and E79's common-map share of APL with the claws' share of the
    rest (APL from the first halves, K_b from the second; `common_map_share`)."""
    f = rec['frames']
    ha, hk = stim_halves(f, cells['C1']), stim_halves(f, cells['C2'])
    st = [s for s in kept if s in ha and s in hk]
    if len(st) < 2:
        return {}
    A = [np.stack([ha[s][h] for s in st]) for h in (0, 1)]
    K = [remove_leak(np.stack([hk[s][h] for s in st]), A[h], b_g) for h in (0, 1)]
    oa = odor_specific_parts((A[0] + A[1]) / 2)
    ok = odor_specific_parts((K[0] + K[1]) / 2)
    cc = [cross_correlogram(cells_to_grid(oa[i], cells['grid']), cells_to_grid(ok[i], cells['grid']),
                            cells['grid']['cell_um']) for i in range(len(st))]
    norm = lambda M_: M_ / M_.mean(axis=1, keepdims=True)
    cm = common_map_share(norm(A[0]), norm(K[1]))
    return dict(cc_length_um_mean=float(np.nanmean([x['length_um'] for x in cc])),
                cc_peak_within_one_cell_share=float(np.mean([x['peak_within_one_cell'] for x in cc])),
                cc_zero_lag_mean=float(np.nanmean([x['zero_lag'] for x in cc])), **cm)


def part_d(n_perm: int = N_PERM, cell_ums: Sequence[float] = (CELL_UM,) + CELL_UM_BESIDE) -> dict:
    """Part D (X3's statistic): for every recording of sets L, L+ and V0 with enough kept stimuli, `x3_rows` at b = 0, b_G
    and 2 b_G under the three readings and each cell size (d_x3.csv), the scan over B_SCAN on the carrying reading at
    4 um (d_x3_scan.csv), and the descriptions (`d_described`; d_described.csv). The flies' classes and the verdict are
    formed in `score`."""
    _gate('D')
    out = out_dir()
    bl = _read_json(out / 'b_leak.json')
    if bl is None:
        raise RuntimeError('part D needs part B first (b_leak.json)')
    b_g = float(bl['b_G'])
    rows, scan, desc = [], [], []
    recs, left = two_channel(read_groups(('L', 'L+', 'V0')))
    write_table(pd.DataFrame(left, columns=['stem', 'group', 'why']), out / 'd_left_out.csv')
    for r in recs:
        kept = recording_kept(r, b_g)
        for cu in cell_ums:
            cells = recording_cells(r, cu)
            rows += [dict(x, cell_um=cu) for x in x3_rows(r, cells, b_g, kept, n_perm=n_perm)]
            if cu == CELL_UM:
                scan += x3_rows(r, cells, b_g, kept, multiples=B_SCAN, readings=(CARRYING,), n_perm=n_perm)
                desc.append(dict(recording=r['stem'], fly=r['fly'], group=r['group'], n_kept=len(kept),
                                 registration_matched=r.get('registration_matched'),
                                 C1_tag=r.get('files', {}).get('C1_tag'), C2_tag=r.get('files', {}).get('C2_tag'),
                                 **d_described(r, cells, b_g, kept)))
    write_table(pd.DataFrame(rows), out / 'd_x3.csv')
    write_table(pd.DataFrame(scan), out / 'd_x3_scan.csv')
    write_table(pd.DataFrame(desc), out / 'd_described.csv')
    return dict(n_rows=len(rows))


# ---- part E

def e_rows(rec: dict, cells: dict, b_g: float, multiples: Sequence[float] = (0.0, 1.0)) -> List[dict]:
    """Part E for one recording: K_b's maps (4 um cells) of pa, paa and ms at their two highest concentrations and of
    every concentration of pa, the correlation distance between stimuli across units (cells) as measured (the means of
    all presentations) and with the trial noise removed through the halves (the cross-half r between two stimuli over the
    geometric mean of their repeat r, E90's `disattenuate`)."""
    from experiments import E90_model_stages as E90
    f = rec['frames']
    ha, hk = stim_halves(f, cells['C1']), stim_halves(f, cells['C2'])
    by_odor: Dict[str, List[Tuple[float, str]]] = {}
    for s in hk:
        o, c = s.split(' @ ')
        by_odor.setdefault(o, []).append((float(c), s))
    st = []
    for o in ('pa', 'paa', 'ms'):
        st += [s for _, s in sorted(by_odor.get(o, []))[-2:]]
    st += [s for _, s in sorted(by_odor.get('pa', [])) if s not in st]
    rows = []
    for bm in multiples:
        K = [np.stack([remove_leak(hk[s][h], ha[s][h], bm * b_g) for s in st]) for h in (0, 1)]
        z1, z2 = _row_z(K[0]), _row_z(K[1])
        r = (z1 @ z2.T + z2 @ z1.T) / 2.0
        nr = E90.disattenuate(r)
        meas = 1.0 - np.corrcoef((K[0] + K[1]) / 2)
        for i in range(len(st)):
            for j in range(i + 1, len(st)):
                rows.append(dict(recording=rec['stem'], fly=rec['fly'], b_mult=bm, stim_a=st[i], stim_b=st[j],
                                 distance_measured=float(meas[i, j]), distance_noise_removed=float(nr[i, j])))
    return rows


E_MODEL_STIMULI = ('pa', 'ms', 'paa', 'paa [validation2/diag] @ -5', 'paa [kiwi-control/diag] @ -5',
                   'paa [megamat/diag] @ -5')    # B1's labels: pa, ms and paa at the panel's input, paa at 10^-5 (three sessions)


def e_model_rows(wirings: Sequence[str] = ALL_WIRINGS) -> List[dict]:
    """Part E's model side (described): per wiring (A3's saved runs, APL on) and claw version (CLAW_SHARES, E90's
    `claw_signal` scales as `claw_totals` forms them), the correlation distance (1 - Pearson r over the claws of the 36
    input glomeruli) between E_MODEL_STIMULI, each claw's window mean minus E90's pre-odor mean (the model has no trial
    noise). Each claw's vector at the model's carrying level (its window mean minus its window mean at the blank, the
    amendment of 2026-11-15 18:08), with E90's level beside. Wirings without an A3 run are left out."""
    rows = []
    input_gl = [str(g) for g in H.apl16_orn_deltas().index]
    for w in wirings:
        run = load_a3_run(w)
        if run is None:
            continue
        on, meta = run['on'], run['meta']
        ch = on['claw_win'].astype(float) - on['claw_pre'].astype(float)
        counts = on['spike_counts'].astype(float)
        panel = list(H.APL16_ODORS)
        s0, s1 = _rms(ch[panel].to_numpy()), _rms(counts[panel].to_numpy())
        ix = pd.Index(counts.index.get_level_values(H.KC_ID_LEVEL)).get_indexer(np.asarray(meta['kc_of_claw']))
        m = np.isin(np.asarray(meta['claw_glomerulus']).astype(str), input_gl)
        st = [s for s in E_MODEL_STIMULI if s in ch.columns]
        win = on['claw_win'].astype(float)
        if BLANK_LABEL not in win.columns:
            raise ValueError(f'{w}: no blank in the A3 run')
        K = counts[st].to_numpy()[ix]
        Kb = counts[[BLANK_LABEL]].to_numpy()[ix]
        for level in (A0_CARRYING, A0_SHOWN):
            # the carrying level: each claw's window mean minus its window mean at the blank (the spike signal likewise);
            # E90's level beside: window minus E90's pre-odor mean (the spike signal has no pre-odor part)
            Xd = (win[st].to_numpy() - win[[BLANK_LABEL]].to_numpy()) if level == A0_CARRYING else ch[st].to_numpy()
            Ks = (K - Kb) if level == A0_CARRYING else K
            for share in CLAW_SHARES:
                X = (1.0 - share) * Xd / s0 + share * Ks / s1 if share else Xd
                D = 1.0 - np.corrcoef(X[m].T)
                for i in range(len(st)):
                    for j in range(i + 1, len(st)):
                        rows.append(dict(source='model', wiring=w, level=level, level_role=A0_ROLE[level],
                                         claw_share=share, stim_a=st[i], stim_b=st[j], distance=float(D[i, j]),
                                         n_claws=int(m.sum())))
    return rows


def e_remy_rows(pairs: Sequence[Tuple[str, str]] = (('pa', 'ms'),)) -> List[dict]:
    """Remy's claws (E90's readers; megamat panel, her window): per fly the correlation distance between two odors as
    measured (`remy_window_means`) and with the trial noise removed (`remy_noise_removed`). Empty when her matrices are not
    there."""
    from experiments import E90_model_stages as E90
    try:
        _, _, odors = E90.remy_matrices('claw', 'stim')
        meas, nr = E90.remy_window_means('claw'), E90.remy_noise_removed('claw')
    except FileNotFoundError:
        return []
    rows = []
    for a, b in pairs:
        if a not in odors or b not in odors:
            continue
        i, j = odors.index(a), odors.index(b)
        for f in range(meas.shape[0]):
            rows.append(dict(source='Remy claws', wiring='', fly=f'fly{f + 1}', stim_a=a, stim_b=b,
                             distance=float(meas[f, i, j]), distance_noise_removed=float(nr[f, i, j])))
    return rows


def part_e() -> dict:
    """Part E (described, decides nothing): `e_rows` for every recording of set L (e_claw_distances.csv); beside, the
    family's claw stage (`e_model_rows`) and Remy's claws (`e_remy_rows`) (e_model_distances.csv)."""
    _gate('E')
    out = out_dir()
    bl = _read_json(out / 'b_leak.json')
    if bl is None:
        raise RuntimeError('part E needs part B first (b_leak.json)')
    rows = []
    recs, left = two_channel(read_groups(('L',)))
    for r in recs:
        rows += e_rows(r, recording_cells(r), float(bl['b_G']))
    write_table(pd.DataFrame(rows), out / 'e_claw_distances.csv')
    model = e_model_rows() + e_remy_rows()
    write_table(pd.DataFrame(model), out / 'e_model_distances.csv')
    return dict(n_rows=len(rows), n_model_rows=len(model), left_out=left)


# ---- the score

def load_model_x2() -> pd.DataFrame:
    """A3's model side of X2 (a3_x2_sets.csv.gz)."""
    return pd.read_csv(out_dir() / 'a3_x2_sets.csv.gz')


def x2_model_classes(sets: pd.DataFrame, rid: str, quantity: str, lo: float, hi: float, a0: str,
                     subset: str = 'input') -> dict:
    """The model's class for one ratio and quantity against the flies' interval: per estimate and wiring E98's `classify`
    over the quantity's set (APL over APL_NS and FORMS; the claws over CLAW_NS, FORMS and the versions; Q over every
    combination, `q_classes`), the comparison verdict over the eight wirings per estimate, and over the estimates E98's
    V1 rule (`v1_verdict`; r3 has one)."""
    from experiments import E98_apl_response_strength as E98
    s = sets[(sets['ratio_id'] == rid) & (sets['role'] == 'main') & (sets['baseline'] == a0)]
    per_est, per_wiring = {}, []
    if quantity == 'Q':
        q = q_classes(s, lo, hi, a0, subset)
        for est, g in q.groupby('estimate', sort=False):
            per_est[est] = E98.comparison_verdict(g['cls'])
            per_wiring += g.assign(quantity='Q').to_dict('records')
    else:
        x = s[(s['quantity'] == quantity) & ((s['subset'] == subset) if quantity == 'claws' else True)]
        for (est, w), g in x.groupby(['estimate', 'wiring'], sort=False):
            r = g['ratio'].to_numpy(float)
            r = r[np.isfinite(r) & (r > 0)]
            per_wiring.append(dict(wiring=w, estimate=est, quantity=quantity, cls=E98.classify(r, lo, hi)))
        pw = pd.DataFrame(per_wiring)
        for est, g in pw.groupby('estimate', sort=False):
            per_est[est] = E98.comparison_verdict(g['cls'])
    cls = E98.v1_verdict(per_est) if per_est else 'not decided'
    return dict(cls=cls, per_estimate=per_est, per_wiring=per_wiring)


KC_QUANTITIES = ('claws', 'Q')           # X2's quantities that read the KC channel (K_b); APL's does not


def bg_undecided(cls: str, b_mult: float, reads_kc: bool, cross_check_passed: bool) -> str:
    """B1's rule: if the cross-check failed, the slope does not carry between flies and every KC quantity at b_G is not
    decided (the entry's B1); other classes pass unchanged."""
    return 'not decided' if (not cross_check_passed and reads_kc and float(b_mult) == 1.0) else cls


def x2_score(sets: pd.DataFrame, flies: pd.DataFrame, a0: str, keys: Sequence[Tuple[str, float]],
             cross_check_passed: bool = True) -> Tuple[dict, List[dict]]:
    """X2 under one pre-odor level of the model `a0`: per quantity, reading and b (`keys`), each ratio's class
    (`x2_model_classes` against the flies' hull of c_x2_flies.csv) for set L and for its unprobed flies, the quantity's
    class (`x2_quantity_verdict`), and the class held over the readings, b and the unprobed flies (`held_class`; at least
    HELD_X12_MIN of them remaining). With B1's cross-check failed the claws' and Q's classes at b_G are not decided
    (`bg_undecided`; the class as computed is kept in the rows as cls_if_slope_carried). Returns {quantity: {held,
    by_reading}} and the long rows of score_x2.csv."""
    x2, rows = {}, []
    for q in X2_QUANTITIES:
        by, by_un = {}, {}
        for (rd, bm) in keys:
            for who, store in (('L', by), ('L unprobed', by_un)):
                f = flies[(flies['flies'] == who) & (flies['reading'] == rd) & (flies['b_mult'] == bm)
                          & (flies['quantity'] == q)]
                rc = {}
                for rid, _, _ in X2_RATIOS:
                    row = f[f['ratio_id'] == rid]
                    if row.empty or not np.isfinite(row['lo'].iloc[0]):
                        rc[rid] = 'not decided'
                        rows.append(dict(a0=a0, quantity=q, flies=who, reading=rd, b_mult=bm, ratio_id=rid, lo=np.nan,
                                         hi=np.nan, n_flies=int(row['n_flies'].iloc[0]) if len(row) else 0,
                                         cls='not decided', per_estimate=''))
                        continue
                    lo, hi = float(row['lo'].iloc[0]), float(row['hi'].iloc[0])
                    mc = x2_model_classes(sets, rid, q, lo, hi, a0)
                    rc[rid] = bg_undecided(mc['cls'], bm, q in KC_QUANTITIES, cross_check_passed)
                    rows.append(dict(a0=a0, quantity=q, flies=who, reading=rd, b_mult=bm, ratio_id=rid, lo=lo, hi=hi,
                                     n_flies=int(row['n_flies'].iloc[0]), cls=rc[rid], cls_if_slope_carried=mc['cls'],
                                     per_estimate='; '.join(f'{e}: {c}' for e, c in mc['per_estimate'].items()),
                                     per_wiring='; '.join(f"{w['wiring']}/{w['estimate']}: {w['cls']}"
                                                          for w in mc['per_wiring'])))
                store[(rd, bm)] = dict(ratios=rc, cls=x2_quantity_verdict(rc))
        un = {v['cls'] for v in by_un.values()}
        fu = flies[(flies['flies'] == 'L unprobed') & (flies['quantity'] == q)]
        n_un = int(fu['n_flies'].min()) if len(fu) else 0
        x2[q] = dict(held=held_class({kk: v['cls'] for kk, v in by.items()}, un.pop() if len(un) == 1 else 'not decided',
                                     n_un, HELD_X12_MIN),
                     unprobed_classes=sorted({v['cls'] for v in by_un.values()}), n_unprobed_flies_min=n_un,
                     by_reading={f'{a} b{b:g}': v for (a, b), v in by.items()},
                     by_reading_unprobed={f'{a} b{b:g}': v for (a, b), v in by_un.items()})
    return x2, rows


CARRYING_KEYS = ((CARRYING, 0.0), (CARRYING, 1.0))   # the flies' carrying reading at b = 0 and b_G


def partner_held(held: str, partner_classes: Sequence[str], quantity: str) -> str:
    """The held clause of the amendment of 2026-11-15 18:08: APL's and Q's class is held only if the model's partner level
    (A0_PARTNER) gives the same class (`partner_classes`: the quantity's classes under the partner with the flies'
    carrying reading at b = 0 and b_G); otherwise the quantity is reading-dependent and not decided. The claws' class
    passes unchanged (they have no second level)."""
    if quantity not in PARTNER_QUANTITIES or held == 'not decided':
        return held
    return held if all(c == held for c in partner_classes) else 'not decided'


def _ratio_details(sets: pd.DataFrame, flies: pd.DataFrame, a0: str, q: str, rd: str, bm: float, cross: bool) -> dict:
    """Per ratio of X2_RATIOS, one quantity's class against set L's flies (reading `rd`, b `bm`) under the model's level
    `a0`: the flies' interval, the class (`bg_undecided` applied), and the classes per estimate and per wiring."""
    f = flies[(flies['flies'] == 'L') & (flies['reading'] == rd) & (flies['b_mult'] == bm) & (flies['quantity'] == q)]
    out = {}
    for rid, _, _ in X2_RATIOS:
        row = f[f['ratio_id'] == rid]
        if row.empty or not np.isfinite(row['lo'].iloc[0]):
            out[rid] = dict(cls='not decided', lo=np.nan, hi=np.nan, per_estimate={}, per_wiring={})
            continue
        lo, hi = float(row['lo'].iloc[0]), float(row['hi'].iloc[0])
        mc = x2_model_classes(sets, rid, q, lo, hi, a0)
        out[rid] = dict(cls=bg_undecided(mc['cls'], bm, q in KC_QUANTITIES, cross), lo=lo, hi=hi,
                        n_flies=int(row['n_flies'].iloc[0]), per_estimate=mc['per_estimate'],
                        per_wiring={f"{w['estimate']} | {w['wiring']}": w['cls'] for w in mc['per_wiring']})
    return out


def _counts(per_wiring: Dict[str, str]) -> str:
    """'<estimate>: inside 7, flatter 1; ...' from {'<estimate> | <wiring>': class}."""
    from collections import Counter
    by: Dict[str, Counter] = {}
    for k, c in per_wiring.items():
        by.setdefault(k.split(' | ')[0], Counter())[c] += 1
    return '; '.join(f"{e}: {', '.join(f'{c} {n}' for c, n in sorted(cnt.items(), key=lambda x: -x[1]))}"
                     for e, cnt in by.items())


def x2_levels(sets: pd.DataFrame, flies: pd.DataFrame, keys, cross: bool, unprobed: bool = True,
              lines: bool = True) -> Tuple[dict, List[dict]]:
    """X2 under the amendment of 2026-11-15 18:08, per quantity: the class held over the flies' readings and b and (with
    `unprobed`) the unprobed flies under the model's carrying level (A0_CARRYING), then the partner clause
    (`partner_held`: APL and Q, the partner A0_PARTNER under the flies' carrying reading at b = 0 and b_G), and E90's level
    (A0_SHOWN) shown, deciding nothing; with `lines`, the verdict line (`x2_line`). Returns the per-quantity dict and the
    long rows of every level (score_x2.csv)."""
    main, rows_m = x2_score(sets, flies, A0_CARRYING, keys, cross)
    part, rows_p = x2_score(sets, flies, A0_PARTNER, CARRYING_KEYS, cross)
    shown, rows_s = x2_score(sets, flies, A0_SHOWN, keys, cross)
    out = {}
    for q in X2_QUANTITIES:
        def held(x):
            return x[q]['held'] if unprobed else held_class(
                {kk: r['cls'] for kk, r in x[q]['by_reading'].items()}, None, 0, 0)
        h0 = held(main)
        pc = [part[q]['by_reading'][f'{CARRYING} b{bm:g}']['cls'] for _, bm in CARRYING_KEYS]
        out[q] = dict(held=partner_held(h0, pc, q), held_before_partner=h0,
                      set_l_classes=sorted({r['cls'] for r in main[q]['by_reading'].values()}),
                      unprobed_classes=main[q]['unprobed_classes'],
                      n_unprobed_flies_min=main[q]['n_unprobed_flies_min'],
                      partner_classes={f'{CARRYING} b{bm:g}': c for (_, bm), c in zip(CARRYING_KEYS, pc)},
                      e90_held=held(shown), carrying=main[q], partner=part[q], e90=shown[q])
        if lines:
            out[q]['line'] = x2_line(sets, flies, q, out[q], cross)
    rows = ([dict(r, level_role=A0_ROLE[A0_CARRYING]) for r in rows_m] + [dict(r, level_role=A0_ROLE[A0_PARTNER])
            for r in rows_p] + [dict(r, level_role=A0_ROLE[A0_SHOWN]) for r in rows_s])
    return out, rows


def x2_line(sets: pd.DataFrame, flies: pd.DataFrame, q: str, res: dict, cross: bool) -> dict:
    """One quantity's verdict line (the amendment's order): the class under the flies' carrying reading and the model's
    carrying level, with each ratio's class and count of wirings, at b = 0 and b_G; the flies' partners; the model's
    partner level; E90's level, deciding nothing, and where it gives another class the wirings that change."""
    det = {bm: _ratio_details(sets, flies, A0_CARRYING, q, CARRYING, bm, cross) for bm in (0.0, 1.0)}
    e90 = {bm: _ratio_details(sets, flies, A0_SHOWN, q, CARRYING, bm, cross) for bm in (0.0, 1.0)}
    c = res['carrying']['by_reading']
    txt = [f"{q}: {res['held']} (held{'; before the partner clause ' + res['held_before_partner'] if res['held'] != res['held_before_partner'] else ''})."]
    for bm in (0.0, 1.0):
        rs = '; '.join(f"{rid} {d['cls']} [{d['lo']:.3g}, {d['hi']:.3g}] ({_counts(d['per_wiring'])})"
                       for rid, d in det[bm].items())
        txt.append(f"Carrying (flies' carrying reading, model at the blank), b {bm:g} b_G: {c[f'{CARRYING} b{bm:g}']['cls']}: {rs}.")
    txt.append("Flies' partners (model at the blank): " + ', '.join(
        f"{rd} b {bm:g} b_G {c[f'{rd} b{bm:g}']['cls']}" for rd in READINGS if rd != CARRYING for bm in (0.0, 1.0)) + '.')
    if q in PARTNER_QUANTITIES:
        txt.append("Model's partner level (APL settled, the claws at the blank), flies' carrying reading: " + ', '.join(
            f'{k} {v}' for k, v in res['partner_classes'].items()) + '.')
    else:
        txt.append("Model's partner level: none for the claws (A3 saved no settled level for them).")
    ec = res['e90']['by_reading']
    changes = []
    for bm in (0.0, 1.0):
        for rid in det[bm]:
            pw0, pw1 = det[bm][rid]['per_wiring'], e90[bm][rid]['per_wiring']
            ch = [f"{k} {pw0[k]} -> {pw1[k]}" for k in pw0 if k in pw1 and pw0[k] != pw1[k]]
            if ch:
                changes.append(f"b {bm:g} b_G {rid}: " + '; '.join(ch))
    differs = any(ec[f'{CARRYING} b{bm:g}']['cls'] != c[f'{CARRYING} b{bm:g}']['cls'] for bm in (0.0, 1.0)) or \
        res['e90_held'] != res['held_before_partner']
    txt.append("E90's level (decides nothing), flies' carrying reading: " + ', '.join(
        f"b {bm:g} b_G {ec[f'{CARRYING} b{bm:g}']['cls']}" for bm in (0.0, 1.0)) + f"; held {res['e90_held']}"
        + ('; another class than the carrying level' if differs else '; the same classes') + '.')
    return dict(text=' '.join(txt), e90_differs=bool(differs), wirings_that_change_under_e90=changes,
                carrying_ratios={f'b{bm:g}': {rid: dict(cls=d['cls'], lo=d['lo'], hi=d['hi'], counts=_counts(d['per_wiring']),
                                                        per_estimate=d['per_estimate'])
                                              for rid, d in det[bm].items()} for bm in (0.0, 1.0)},
                e90_ratios={f'b{bm:g}': {rid: dict(cls=d['cls'], counts=_counts(d['per_wiring']))
                                         for rid, d in e90[bm].items()} for bm in (0.0, 1.0)})


def score(model_sets: Optional[pd.DataFrame] = None) -> dict:
    """The verdicts (verdict.json, score_x2.csv, score_x3_flies.csv), from parts A3, A4, B, C and D, on set L with flies as
    units, each held (`held_class`) over b in (0, b_G), the carrying reading and both partners, and the unprobed flies:
    X1 (`x1_table`'s verdicts), X2 under the model's carrying level with the partner clause and E90's level shown
    (`x2_levels`, the amendment of 2026-11-15 18:08), X3 (`fly_x3` and `x3_verdict` per reading, with B1's cross-check and
    A4's rule); B1's cross-check rule (`bg_undecided`); the channel check's verdict and Yang's quality caveat recorded;
    the readings beside (`beside`)."""
    _gate('score')
    out = out_dir()
    bl, x1, rule = _read_json(out / 'b_leak.json'), _read_json(out / 'c_x1.json'), _read_json(out / 'a4_rule.json')
    if bl is None or x1 is None or rule is None:
        raise RuntimeError('score needs parts A4, B and C first')
    chk = _read_json(out / 'channel_check.json')
    keys = [(rd, bm) for bm in (0.0, 1.0) for rd in READINGS]
    res = dict(caveat=QUALITY_CAVEAT, channel_mode=CHANNEL_MODE,
               channel_check=None if chk is None else chk.get('verdict'),
               model_levels=dict(carrying=A0_CARRYING, partner=A0_PARTNER, shown_deciding_nothing=A0_SHOWN),
               b_G=bl['b_G'], outcome_h=bl['outcome_h'], outcome_h_b5=bl.get('outcome_h_b5'),
               cross_check_passed=bl['cross_check_passed'],
               b_g_alone=rule['b_g_alone'])
    k = lambda rd, bm: f'{rd} b{bm:g}'
    cross = bool(bl['cross_check_passed'])
    res['cross_check_rule'] = (None if cross else
                               "B1's cross-check failed: the slope does not carry between flies, so every KC quantity at "
                               "b_G (X1; X2's claws and Q) is not decided, and X3's 'follows' cannot be given (the entry's "
                               "B1 and X3); the classes as computed at b_G are kept beside (if_slope_carried)")
    x1_raw = {kk: x1['all'][k(*kk)]['verdict'] for kk in keys if k(*kk) in x1['all']}
    x1_cls = {kk: bg_undecided(v, kk[1], True, cross) for kk, v in x1_raw.items()}
    x1_un = {kk: bg_undecided(x1['unprobed'][k(*kk)]['verdict'], kk[1], True, cross) for kk in keys
             if k(*kk) in x1['unprobed']}
    un = set(x1_un.values())
    n_un = min((int(v['n_flies']) for kk in keys if k(*kk) in x1['unprobed']
                for v in x1['unprobed'][k(*kk)]['pairs'].values()), default=0)
    res['X1'] = dict(held=held_class(x1_cls, un.pop() if len(un) == 1 else 'not decided', n_un, HELD_X12_MIN),
                     unprobed_classes=sorted(set(x1_un.values())), n_unprobed_flies_min=n_un,
                     by_reading={f'{a} b{b:g}': v for (a, b), v in x1_cls.items()},
                     by_reading_if_slope_carried={f'{a} b{b:g}': v for (a, b), v in x1_raw.items()},
                     by_reading_unprobed={f'{a} b{b:g}': v for (a, b), v in x1_un.items()},
                     pairs={kk_: x1['all'][kk_]['pairs'] for kk_ in x1['all']})
    sets = load_model_x2() if model_sets is None else model_sets
    flies = read_table(out / 'c_x2_flies.csv')
    x2, x2_rows = x2_levels(sets, flies, keys, cross)
    write_table(pd.DataFrame(x2_rows), out / 'score_x2.csv')
    res['X2'] = dict(levels=dict(carrying=A0_CARRYING, partner=A0_PARTNER, shown_deciding_nothing=A0_SHOWN),
                     held={q: v['held'] for q, v in x2.items()}, lines={q: v['line']['text'] for q, v in x2.items()},
                     quantities={q: {kk_: vv for kk_, vv in v.items() if kk_ not in ('carrying', 'partner', 'e90')}
                                 for q, v in x2.items()},
                     by_reading={q: dict(carrying=v['carrying']['by_reading'], partner=v['partner']['by_reading'],
                                         e90=v['e90']['by_reading'],
                                         carrying_unprobed=v['carrying']['by_reading_unprobed'])
                                 for q, v in x2.items()})
    d = read_table(out / 'd_x3.csv')
    d = d[(d['group'] == 'L') & (d['cell_um'] == CELL_UM)]
    x3 = {}
    for rd, g in d.groupby('reading', sort=False):
        fx = fly_x3(g)
        x3[rd] = dict(x3_verdict(fx, bool(bl['cross_check_passed']), bool(rule['b_g_alone'])),
                      flies=fx.to_dict('records'))
    res['X3'] = dict(held=held_class({rd: v['verdict'] for rd, v in x3.items()} if len(x3) == len(READINGS)
                                     else {'missing': 'not decided'}, None, 0, 0),
                     by_reading={rd: {kk: v for kk, v in r.items() if kk != 'flies'} for rd, r in x3.items()})
    fl = [dict(f, reading=rd) for rd, r in x3.items() for f in r['flies']]
    write_table(pd.DataFrame(fl), out / 'score_x3_flies.csv')
    if not cross:
        # the X1 and X2 classes as they would be held were the slope to carry (described, deciding nothing)
        x1u_raw = {x1['unprobed'][k(*kk)]['verdict'] for kk in keys if k(*kk) in x1['unprobed']}
        res['X1']['held_if_slope_carried'] = held_class(x1_raw, x1u_raw.pop() if len(x1u_raw) == 1 else 'not decided',
                                                        n_un, HELD_X12_MIN)
        x2c, _ = x2_levels(sets, flies, keys, True, lines=False)
        res['X2']['held_if_slope_carried'] = {q: x2c[q]['held'] for q in KC_QUANTITIES}
    res.update(beside(sets, keys, cross, bool(rule['b_g_alone']), res))
    _write_json(out / 'verdict.json', res)
    return res


ALIGNMENT_NOTE = (
    "Yang's question 18, second half (whether C2 moved by C1's shifts is aligned with C1 pixel for pixel), is not answered "
    "(2026-11-15 17:53). X3 depends on that alignment throughout: its agreement A is a Pearson r over cells between APL's "
    "odor-specific map and the KC channel's from other presentations, and K_b = C2 - b C1 is formed cell by cell, so a "
    "shift of the KC map against APL's lowers A and leaves part of the leak in K_b where APL's map changes steeply (the "
    "cross-check of B1 and the scan over b read the same cells). X1 and X2 read calyx means inside one outline and do not "
    "depend on it. The data's own checks, described: the cross-correlogram's peak lag (d_described.csv: the share of "
    "stimuli whose peak lies within one cell of zero), the 8 um cells beside (less sensitive to a shift of a fraction of a "
    "cell), and the recordings whose two channels' maps were not registered with the same shifts (registration_matched "
    "False).")


def _x1_classes(c: pd.DataFrame, flies: Optional[Sequence[str]] = None, min_flies: int = X1_MIN_FLIES,
                n_boot: int = N_BOOT, cross: bool = True) -> Tuple[dict, dict]:
    """X1's class per reading and b in (0, b_G) on the response rows `c` (`x1_table`; `bg_undecided` with B1's
    cross-check `cross`), and the pairs' summaries."""
    _, x1 = x1_table(c, flies=flies, min_flies=min_flies, n_boot=n_boot)
    keep = {k: v for k, v in x1.items() if k[1] in (0.0, 1.0)}
    return ({f'{rd} b{bm:g}': bg_undecided(v['verdict'], bm, True, cross) for (rd, bm), v in keep.items()},
            {f'{rd} b{bm:g}': v['pairs'] for (rd, bm), v in keep.items()})


def _x3_classes(d: pd.DataFrame, cross: bool, b_g_alone: bool, unprobed: Sequence[str]) -> dict:
    """X3 per reading on part D's rows `d` (`fly_x3`, `x3_verdict`), without the per-fly rows."""
    return {rd: x3_verdict(fly_x3(g), cross, b_g_alone, unprobed) for rd, g in d.groupby('reading', sort=False)}


def beside(sets: pd.DataFrame, keys, cross: bool, b_g_alone: bool, main: dict) -> dict:
    """What the score reports beside the verdicts, deciding nothing: the sensitivity readings (set L without the runs
    that Yang's notes mark as warped by SENSITIVITY_WARP_UM or more or as taken at SENSITIVITY_NM, a1_recordings.csv's
    out_of_sensitivity_readings; Yang, question 18, 2026-11-15 17:53: no run is left out of the carrying readings; and
    the planes apart), V0 as the repetition ("not repeated in the development flies" where its class differs from a held
    one), X3 at the other cell sizes, and the alignment of the two channels' maps (ALIGNMENT_NOTE)."""
    out = out_dir()
    c = read_table(out / 'c_responses.csv')
    d_all = read_table(out / 'd_x3.csv')
    a1 = pd.read_csv(out / 'a1_recordings.csv')
    outs = set(a1[a1['out_of_sensitivity_readings'].fillna('').astype(str) != '']['stem'])
    cL, dL = c[c['group'] == 'L'], d_all[(d_all['group'] == 'L') & (d_all['cell_um'] == CELL_UM)]
    # sensitivity: without the warped and 925 nm runs
    cs, ds = cL[~cL['recording'].isin(outs)], dL[~dL['recording'].isin(outs)]
    x1s, _ = _x1_classes(cs, cross=cross)
    x1su, _ = _x1_classes(cs, flies=UNPROBED_X12, min_flies=HELD_X12_MIN, cross=cross)
    un = set(x1su.values())
    fs = pd.concat([x2_flies(cs).assign(flies='L'), x2_flies(cs, flies=UNPROBED_X12).assign(flies='L unprobed')],
                   ignore_index=True)
    x3s = _x3_classes(ds, cross, b_g_alone, UNPROBED_X3)
    sens = dict(left_out=sorted(outs & set(cL['recording'])),
                left_out_flies=sorted({fly_of(x) for x in outs & set(cL['recording'])}),
                X1=dict(by_reading=x1s, held=held_class(x1s, un.pop() if len(un) == 1 else 'not decided', HELD_X12_MIN,
                                                         HELD_X12_MIN)),
                X2={q: v['held'] for q, v in x2_levels(sets, fs, keys, cross, lines=False)[0].items()},
                X3=dict(by_reading={rd: v['verdict'] for rd, v in x3s.items()},
                        held=held_class({rd: v['verdict'] for rd, v in x3s.items()} if len(x3s) == len(READINGS)
                                        else {'missing': 'not decided'}, None, 0, 0)))
    # the planes apart (carrying reading)
    planes = {}
    for pl in sorted(pd.to_numeric(cL['plane_um'], errors='coerce').dropna().unique()):
        cp = cL[pd.to_numeric(cL['plane_um'], errors='coerce') == pl]
        dp = dL[(pd.to_numeric(dL['plane_um'], errors='coerce') == pl) & (dL['reading'] == CARRYING)]
        x1p, pairs = _x1_classes(cp[cp['reading'] == CARRYING], cross=cross)
        fx = fly_x3(dp) if len(dp) else pd.DataFrame(columns=['fly', 'b_mult', 'cls'])
        n_pass = (len(set(fx[(fx['b_mult'] == 0.0) & (fx['cls'] == 'pass')]['fly'])
                      & set(fx[(fx['b_mult'] == 1.0) & (fx['cls'] == 'pass')]['fly'])) if len(fx) else 0)
        planes[f'{pl:g} um'] = dict(n_recordings=int(cp['recording'].nunique()), n_flies=int(cp['fly'].nunique()),
                                    X1=x1p, X1_pairs=pairs,
                                    X3=x3_verdict(fx, cross, b_g_alone)['verdict'] if len(fx) else 'not decided',
                                    X3_n_flies=int(fx['fly'].nunique()) if len(fx) else 0, X3_n_pass_0_bG=n_pass)
    sens['planes'] = planes
    # V0, the repetition (its own flies; no unprobed rule)
    cV = c[c['group'] == 'V0']
    dV = d_all[(d_all['group'] == 'V0') & (d_all['cell_um'] == CELL_UM)]
    v0_flies = sorted(cV['fly'].unique())
    x1v, x1v_pairs = _x1_classes(cV, cross=cross)
    fv = read_table(out / 'c_x2_flies.csv')
    fv = fv[fv['flies'] == 'V0'].assign(flies='L')
    x3v = _x3_classes(dV, cross, b_g_alone, v0_flies)
    v0 = dict(n_flies=len(v0_flies), n_recordings=int(cV['recording'].nunique()),
              X1=dict(by_reading=x1v, held=held_class(x1v, None, 0, 0), pairs=x1v_pairs),
              X2={q: v['held'] for q, v in x2_levels(sets, fv, keys, cross, unprobed=False, lines=False)[0].items()}
              if len(fv) else {},
              X3=dict(by_reading={rd: v['verdict'] for rd, v in x3v.items()},
                      held=held_class({rd: v['verdict'] for rd, v in x3v.items()} if len(x3v) == len(READINGS)
                                      else {'missing': 'not decided'}, None, 0, 0)))
    flags, undecided = [], []

    def compare(name, h, hv):
        if h == 'not decided':
            return
        if hv == 'not decided':
            undecided.append(f'{name} {h}: the development flies do not decide it')
        elif hv != h:
            flags.append(f'{name} {h}: not repeated in the development flies ({hv})')
    compare('X1', main['X1']['held'], v0['X1']['held'])
    for q in X2_QUANTITIES:
        compare(f'X2 {q}', main['X2']['held'][q], v0['X2'].get(q, 'not decided'))
    compare('X3', main['X3']['held'], v0['X3']['held'])
    v0['not_repeated'] = flags
    v0['not_decided_in_v0'] = undecided
    # X3 at the other cell sizes (carrying reading's classes per reading)
    cells = {}
    for cu in CELL_UM_BESIDE:
        dc = d_all[(d_all['group'] == 'L') & (d_all['cell_um'] == cu)]
        if len(dc):
            cells[f'{cu:g} um'] = {rd: v['verdict'] for rd, v in _x3_classes(dc, cross, b_g_alone, UNPROBED_X3).items()}
    desc = read_table(out / 'd_described.csv') if (out / 'd_described.csv').exists() else pd.DataFrame()
    dl = desc[desc['group'] == 'L'] if len(desc) else desc
    align = dict(note=ALIGNMENT_NOTE,
                 set_l_not_registered_together=sorted(dl[dl['registration_matched'] == False]['recording'])   # noqa: E712
                 if len(dl) and 'registration_matched' in dl else [],
                 cc_peak_within_one_cell_share_median=float(dl['cc_peak_within_one_cell_share'].median())
                 if len(dl) and 'cc_peak_within_one_cell_share' in dl else np.nan,
                 X3_other_cell_sizes=cells)
    return dict(sensitivity=sens, V0=v0, alignment=align)


PARTS = ('0', 'A1', 'A2', 'A3', 'A4', 'B', 'C', 'D', 'E', 'score', 'A2-blank', 'C4-model', 'C-strengths', 'variants')
READS = {
    '0': "one released 2021 recording (E100's reader); with --e90, E90's saved stage matrices",
    'A1': "viki's recordings.csv, the stimulus yaml files and Yang's date notes (no map values)",
    'A2': "E90's saved stage matrices (data/model_outputs/E90_model_stages/B/)",
    'A3': "C2's family on E98's B1 input (E83's saved fits); E98's saved B1 runs for the identity",
    'A4': 'nothing (made-up maps)',
    'B': "viki's copy: groups G, R, M and L, both channels (after Yang's release)",
    'C': "viki's copy: sets L, L+, V0 and VT, calyx means of both channels (after the channel check passed)",
    'D': "viki's copy: sets L, L+ and V0, both channels' maps (after the channel check passed)",
    'E': "viki's copy: set L's maps; A3's claw stage; Remy's claw matrices (after the channel check passed)",
    'score': 'the tables of parts A3, A4, B, C and D (after the channel check passed)',
    'A2-blank': "A3's family table (a3_family.csv; no simulation)",
    'C4-model': "A3's family table (a3_family.csv; no simulation)",
    'C-strengths': 'c_responses.csv (after the channel check passed)',
    'variants': "a3_runs.csv, a3_family.csv, a3_identity.json and E83's members (the fits' folders and tags)",
}
WRITES = {
    '0': 'identity.json', 'A1': 'a1_recordings.csv, a1_groups.csv',
    'A2': 'a2_family.csv, a2_lines.csv, a2_ratios.csv',
    'A3': 'a3_runs.csv, a3_stimuli.csv, a3_family.csv, a3_ratios.csv, a3_lines.csv, a3_identity.json; '
          'outputs/E101_apl_vs_kc_claws/A3/',
    'A4': 'a4_power.csv, a4_rule.json',
    'B': 'channel_check.csv, channel_check.json, b_leak.json, b1_slopes.csv, b1_cross_check.csv, b2_reverse.csv, '
         'b3_marker.csv, b4_solvent.csv, b5_slopes.csv',
    'C': 'c_responses.csv, c_left_out.csv, c_x1_recordings.csv, c_x1.json, c_x2_flies.csv',
    'D': 'd_x3.csv, d_x3_scan.csv, d_described.csv, d_left_out.csv',
    'E': 'e_claw_distances.csv, e_model_distances.csv',
    'score': 'verdict.json, score_x2.csv, score_x3_flies.csv',
    'A2-blank': 'a2_blank_lines.csv, a2_blank_ratios.csv, a2_blank.json',
    'C4-model': 'c4_hallem_levels.csv, c4_hallem_ratios.csv, c4_apl_off.csv, c4_model.json',
    'C-strengths': 'c_strengths.csv',
    'variants': 'variants.csv (the archive\'s params.csv source)',
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=PARTS)
    ap.add_argument('--wirings', default=','.join(ALL_WIRINGS), help='A3: comma-separated wirings')
    ap.add_argument('--e90', action='store_true', help="part 0: also check against E90's saved stage matrices")
    ap.add_argument('--runs-only', action='store_true', help="A3: the runs and their identities only, no tables")
    ap.add_argument('--apl-only', action='store_true',
                    help="A3: the runs, their identities and the APL side's tables (no claw stage read)")
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
        provenance('A3', wirings=list(w), runs_only=a.runs_only)
        _process_commits()
        if a.runs_only or a.apl_only:
            res = part_a3_runs(w)
            _write_json(out_dir() / 'a3_identity.json', res)
            if a.apl_only:
                res = dict(identity=res, apl=a3_apl_tables(w))
        else:
            res = part_a3(w)
        print(json.dumps(res, indent=1, default=float), flush=True)
    elif a.part == 'A4':
        provenance('A4', n_flies=a.n_flies, n_perm=a.n_perm)
        print(json.dumps(part_a4(a.n_flies, a.n_perm), indent=1, default=float), flush=True)
    elif a.part in ('A2-blank', 'C4-model', 'variants'):
        provenance(a.part)
        fn = {'A2-blank': part_a2_blank, 'C4-model': part_c4_model,
              'variants': lambda: dict(n_rows=len(write_variants()))}[a.part]
        print(json.dumps(fn(), indent=1, default=float), flush=True)
    elif a.part == 'C-strengths':
        provenance(a.part, **data_provenance())
        print(json.dumps(part_c_strengths(), indent=1, default=float), flush=True)
    elif a.part == 'score':
        provenance('score', levels=dict(carrying=A0_CARRYING, partner=A0_PARTNER, shown=A0_SHOWN), **data_provenance())
        print(json.dumps(score(), indent=1, default=float), flush=True)
    else:
        provenance(a.part, **data_provenance())
        print(json.dumps({'B': part_b, 'C': part_c, 'D': part_d, 'E': part_e}[a.part](), indent=1, default=float),
              flush=True)


if __name__ == '__main__':
    main()
