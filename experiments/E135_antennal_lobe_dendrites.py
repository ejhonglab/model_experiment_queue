#!/usr/bin/env python3
"""E135, the antennal lobe held to the recorded PN dendrites: the step from ORNs to PNs read glomerulus by glomerulus,
the measured step and then global parameters and terms per glomerulus fitted to the dendrites and read on held-out
flies and odors, and the recorded dendrites put in the model's PN stage.

The entry is item E135 of experiment_priorities.md (proposed by viki 2026-11-24 09:36 at Yang's request; planned by a
Fable subagent, claude-fable-5-2; moved in by viki 2026-11-24 10:04, revised on moving after Tom's words of 09:38, about
09:41 and 09:55, Yang's OK of 09:43 and hal's notes on Zocchi et al. 2022). It builds on E134's shared code and does not
copy it: viki's data_prep/stage_population_data.py and its committed lock (experiments/E134_lock.json), hal's
kcq/population_fit.py (the views and one-fly units), atlas's kcq/stage_readings.py (the model read as each stage is
recorded) and experiments/E134_population_baseline.py (lanes, claims, the lock guard), viki's experiments/E134_flies.py
(the flies' recorded steps and the no-step yardstick).

The readings, fixed before any run (from the entry's Definitions):
- The shared glomeruli are the 28 that the ORN recordings, the PN-dendrite recordings and the model's input all hold
  (E90's `shared_glomeruli`, checked against the lock's populations).
- The ORN input of every arm is the model's own megamat input: the lab's recordings at the ORN axon terminals (olfdb
  dataset 1, 9 flies), converted to rates (`megamat_orn_deltas`). It already holds the lateral signal among ORNs. No arm
  uses Hallem's spiking data. On the APL16 input, paa, euc, -aPine and B-myr come from the validation2 terminal
  recordings.
- The stage's curve is one monotone curve from PN rate to dendrite signal, for all glomeruli and odors (stage_readings'
  `fit_stage_curve`, set on pooled values only).
- Imposed: on the shared glomeruli each PN's odor-evoked rate is set from the recorded dendrites, rank for rank on the
  model's own pooled distribution of PN rates. The arrangement over glomeruli and odors is the flies'; the distribution
  of rates and the time course are the model's. The dendrites are the fly means in the lock's trial-noise units, the
  rates read as olfsysm's option reads them (`pn_evoked_rates`), and a rate outside what its PN can reach is set at
  that PN's bound (silent or saturated). The other glomeruli keep the model's PN stage.
- The global parameters are pn_inhsc, pn_inhadd, ln_inhsc, ln_inhadd and the scale of the PN's saturating curve
  (olfsysm's pn.tanhsc).
- A PN's odor-evoked rate is its mean over the stage window (MODEL_WINDOW, closed) minus its mean over the settled
  pre-odor window (SETTLED_PRE_WINDOW, closed on the left), as stage_readings reads the model's PN dendrites.

The planner's rulings on the ten open definitions (claude-fable-5-2, spawned from hal 2026-11-24 12:56; applied to the
entry by hal 13:16) are coded as follows; the lock (`--part lock`) records each.
1. Units. Fits, held-out reads and the imposed table read the dendrites in the lock's trial-noise units (FIT_UNITS). The
   measured step the rungs use reads the ORN terminals in dF/F (the model's input is the terminals' dF/F times one
   factor) and the dendrites in the fit's units (RUNG_STEP_READING). Step 0 gives the measured step in three readings
   (STEP_READINGS: dF/F on both sides, noise units on both sides, and the rungs' reading); a gain whose side of 1
   differs between readings is "unit-dependent" (`unit_dependence`). Before outcome (e) or (f) is declared, the three
   imposed megamat arms are repeated with the table set from the dendrites in dF/F (DFF_REPEAT_ARMS, 24 runs); if the
   class differs, the outcome is "unit-dependent" (`outcome`).
2. A glomerulus's dendrites do not follow its terminals where the 90 % interval of its slope (dendrite on terminal,
   over resampled flies) does not lie above 0 (`follows_rule`), applied in each fit on the flies and odors that fit
   sees; its measured gain is then 1 in R2m (R2f and R3 may still fit it), and each fold's record lists it.
3. Folds and directions: 17 folds of one odor (N_FOLDS); each half of the PN flies (the seeded split of step 0, three
   and four flies) fitted in turn and read on the other; a rung's error is the mean of the two directions.
4. The fit's loss is the affine curve (`matrix_loss`: the matrix view's error in one-fly units, no noise, the 28 shared
   glomeruli, the curve set on the training odors). Every kept decision and reported error reads the lock's default
   reading (KEEP_READING, the quantile curve and the flies' noise); its 'affine_noise' alternative (BESIDE_READING) is
   read beside, and a rung kept in one reading and not the other is "reading-dependent" (the default decides).
5. R1 and a kept rung enter step 2 as tables of PN rates, the folds mapped before assembly. Each fit's table is held at
   the model's own rate scale by one positive factor (`rung_scale_factors`: the pooled SD of its evoked rates over the
   28 shared glomeruli and its training odors equals the unchanged stage's there), set once per fit on the megamat
   input and applied to that fit's APL16 table too; a factor outside SCALE_FACTOR_RANGE is named. The tables come from
   fits to all seven PN flies (`rung_fits_all_flies`): each megamat odor from the fold fit that left it out; paa, euc,
   -aPine and B-myr from the fit to all 17 odors.
6. A rung is kept if it is lower than the highest kept rung below it (R0 counts as kept) by more than KEEP_SE_MULT
   standard errors on the mean of the two directions, and no higher in either direction (`keep_rung`); rungs are still
   fitted in order, each on the one below. A lateral term must gain more than at least LATERAL_PERM_MIN_EXCEEDED of
   N_PERM_LATERAL fits of the same term to odor-permuted dendrites; R3 stops at MAX_LATERAL_TERMS added terms, and
   records whether a fifth would have passed. The permuted refits of R1 and R2 (N_PERM_CONTROL) are reported and decide
   nothing; they run as separate processes by permutation index (`--part step1_control --perm i`), and step 2 does not
   wait for them. In step 2 every permuted arm takes one seeded permutation per wiring (`arm_seed`).
7. Rates outside what a PN can reach are clipped to `pn_rate_bounds` on both sides, BOUNDS_MARGIN_HZ inside them, for
   every table entering by route 1 (imposed, its permutations, R1, the rung; `clip_table`); the same clipped table is
   route 2's target. Each run records the count and the largest excursion; an arm with more than CLIP_NAME_SHARE of
   its set entries clipped is named.
8. The imposed table is set in `pn_evoked_rates`' reading ('preset'); scoring stays E134's reading. Each imposed run
   reports the largest difference between its table and its PN rates in E134's reading and its PN-dendrite matrix error
   (`imposed_run_check`); a difference above IMPOSED_CHECK_HZ is named.
9. Route 1 defines the arms; the route-2 repeat is a check of the option (one megamat imposed run on the real wiring,
   the same clipped table). The glomeruli without input are not solved; their largest drift and the LN's change are
   reported. The routes "agree" where each megamat stage's population error (boutons, claws, somata) differs between
   them by less than that stage's reset threshold (`routes_agree`). No class rests on route 2; a non-converged solve is
   reported, not forced.
10. The thresholds in one-fly units, the somata's 0.15 and +-0.05, are reset to THRESHOLD_SE_MULT standard errors over
   the eight wirings (SD with ddof 1 over sqrt(8)) of the megamat somata's strength-by-odor error (view_strength) in
   E134's runs held at the megamat flies' response rate, each stage's default reading (`class_thresholds` with its
   column; 0.034 in E134's stage table of 15:02). They are frozen with the table's sha256 before step 2's first arm
   (`freeze_thresholds`, THRESHOLDS_FILE) and never re-read. The halving, the half shares, six of eight wirings and the
   300 Hz exclusion stay. The carry classes are read only when "upstream" holds (the rulings of 15:41 below; until
   then this ruling read the somata's population error).
Budget (the rulings'): step 1 up to BUDGET_LANE_HOURS lane-hours. Before the lock is written, one R1 fold fit (and one
of each other fit) is timed on planted data of the real size (`time_step1`, `--part timing`), and if step 1 projects
above the budget the controls drop from 5 to 3, then the folds from 17 to 6 (seeded) (`choose_budget`); the lock
records the choice, and it is never changed after (`write_lock` refuses an existing lock).

The measured step (step 0, recordings only): for each shared glomerulus, the slope through the origin of its dendrite
response on its own ORN-terminal response across odors (fly means), as a multiple of the median glomerulus's (its
measured gain); how closely the dendrites follow the terminals (Pearson r over odors); the regression on its own
response and on the odor's summed terminal response (over the model input's glomeruli); the residual by odor; and the
other input glomeruli whose terminal responses that residual follows best (its candidate partners, ranked by the share
of the residual's variance they explain, with the chance level of the best of them under permuted odor labels).
Intervals from resampled flies (both stages' flies, independently). Each is given on all flies and on both halves of
the PN flies, with their agreement, in each of the three readings. DM1, DL5 and VA2 are shown first. paa and B-myr have
no dendrite recording (the dendrites hold the megamat panel only), so of the three weak-input odors only 2-but is read.

The glomerulus terms: a gain is one factor per shared glomerulus on the drive its PN reads from its own ORN terminals,
between 0.25 and 4; a lateral term adds the terminal response (the ORN delta) of one other glomerulus to that drive,
with a weight between -1 and 1, at most two per glomerulus. olfsysm's LN layer reads only the mean ORN rate over
glomeruli and its PN layer reads its own glomerulus's ORN delta, so the terms act on the PN's drive only and the
normalization keeps the recorded ORN input (`PNStage`, a port of olfsysm's ORN, LN and PN layers checked against them).
What a term stands for: something between the recorded ORN terminals and the PN dendrites; no term is an identified
interaction.

Held to the dendrites: every parameter of this step is fitted to the recorded PN dendrites and to nothing downstream.
No parameter per odor is fitted, and no other stage gains one. The halves, the folds and everything above are fixed in
this experiment's lock (experiments/E135_lock.json, `--part lock`) before any fit.

The rungs: R0, the model's stage as it is. R1, the global parameters. R2, R1 with the gains: first set from the
measured gains with no gain fitted (R2m), then fitted from there (R2f); where the fitted gains are not kept against the
measured ones, the measured gains are the rung. R3, R2 with lateral terms, tried only among each glomerulus's two best
measured candidates and added one at a time for as long as the held-out error falls.

The gates (the entry's Design 1 and Status): step 0 reads recordings only and runs on E134's committed lock. No fit
(step 1) and no model run (step 2) before E134's measure gate reads "usable" and its megamat runs are scored
(`require_e134_gate`); a planner then re-reads the thresholds. The timing helper fits planted data only. Exception for
step 1 (Tom, 2026-11-24 about 15:21, in hal's terminal, after E134's gate read "not usable" at 15:19): "Run step 1 now:
fit the PN stage to the dendrites now on hal's and disco's lanes. Step 2's model runs still wait for the adopted
measure." Step 1 and its controls (and the family's build they read the PN stage from) run when the gate is usable or
when results/.../step1_authorized.json records Tom's authorization (`require_step1_allowed`; hal's session writes it);
step 2 and the score keep `require_e134_gate`. The lock never needed the gate.

Step 2 (up to 73 runs at fixed scales, both response rates re-held, on the family's eight wirings): on the megamat
input R1, the highest kept rung (each odor predicted by the fit that left it out), that rung with its terms' glomerulus
labels permuted, the dendrites imposed, imposed with glomerulus labels permuted, imposed with odor labels permuted; on
the APL16 input R1, the highest kept rung, and the dendrites imposed for its 12 megamat odors; one run repeats an
imposed run by the second route. R1 runs whether kept or not; the rung's 24 runs only if R2 or R3 is the highest kept.
The dF/F repeat's 24 runs only before a closure (when the score's verdict asks for them). Nothing of the family is
refitted. Each point (`run_point2`): the family's PN stage on its input and wiring (`family_stage`; a lane's first point
on a wiring first checks that the unchanged family reproduces its saved run exactly, `lane_identity`), the arm's table
(`arm_table`, clipped), entered by route 1 or 2 (`arm_input`), both response rates re-held (E77's `hold_rates`; on
megamat at E134's held targets, on APL16 at the wiring's own fit's, HOLD_RULES), the held point built again and run with
every stage (E134's `build_point` and `run_stages`), reproducing the hold exactly, its PN rates checked against its table
(`run_pn_check`), an imposed run's check in E134's reading (ruling 8), and the run saved in stage_readings' layout
(`save_point`, locally and in the archive) with its row and provenance.

The score (`part_score`): every lane's runs read and scored as E134 scores its own (atlas's `read_run`, hal's measures
through E134's `score_family`), the class rows against E134's runs of the same wiring held at the megamat flies'
response rate (`class_rows`, `e134_base`), the exclusions (`exclusions`), the classes at the frozen thresholds, the
dF/F repeat's and the second reading's classes, rule 3's guard, the outcome, step 1 read on strength, the routes'
agreement, the megamat stage error per arm, the imposed checks, realism and the APL16 runs' APL effect beside the
family's own fits (`score_tables`).

The two routes into the model's PN stage: route 1, olfsysm's option for the PN stage's rates given directly (viki's
Work claim, Yang's OK of 09:43; olfsysm b1e89d7's fit_mb_model `pn_rates`, `route1_pn_rates`; the odor's LN stays the
model's, and rates outside a PN's reach are clipped first, `clip_table`); route 2, the ORN input solved glomerulus by
glomerulus so that the unchanged PN stage returns the target rates to within SOLVE_TOL_HZ (the clipped table on the set
glomeruli and the model's own rates on the other input glomeruli; Tom, 09:38: an inversion fixed by the recordings, not
a search; `solve_input`, with its identity check `solve_identity`).

The classes (on the eight wirings, against E134's own runs held at the megamat flies' response rate), on strength by
odor, the measure adopted after E134's gate (the rulings of 15:41, from a Fable planner, claude-fable-5-2, spawned from
hal 2026-11-24 15:32): every class reads that view's error in one-fly units in each stage's default reading
(`classify`). "Upstream": the imposed dendrites lower the somata's strength-by-odor error by the threshold or more on
six of eight wirings, by more than both permuted controls, and halve the boutons' strength-by-odor error. "Upstream to
the boutons only": the bouton condition holds and the somata's stays within ± the threshold. "Not upstream" otherwise,
and (f) where the permuted controls gain as much on the somata's strength-by-odor error. "The global parameters carry"
/ "the glomerulus terms carry": R1 / the kept rung give at least half of the imposed fall of the strength-by-odor error
at the boutons and at the somata (read only when "upstream" holds). No class for a run outside a band on the APL16
input or with more than 1 % of bouton responses above 300 Hz (`exclusions`). Beside every arm, wiring and megamat
stage, entering no class: sharing, breadth, the matrix view at the dendrites, the three-view population error and the
two-view error (`class_rows`); the megamat stage error on strength over the four stages with the three-view error
beside, and for the imposed arms both over boutons, claws and somata (`stage_error_table`). Rule 3's guard
(`guard`), for R1 and the kept rung against the family's run of the same wiring: raised where, on six of eight wirings,
the strength-by-odor error or the population error rises by more than 3 standard errors (over the eight wirings of the
same quantity in E134's runs) at the megamat boutons, claws or somata or at the APL16 somata with APL on or off (E134's
step 2, E90:B), or where step 1's held-out strength-by-odor error at the dendrites is above R0's in both directions; a
raised guard is named in the outcome ("at a cost beside") and that stage is not called an improvement. A second reading
(SECOND_READINGS: the somata on raw counts, the boutons on 'curve_noise_g3'): an outcome that differs is
"reading-dependent", the default decides, and (e) or (f) is not declared on it. The routes' check reads strength by odor
with the population error beside. Step 1 is read on strength too: a kept rung above the highest kept rung below it on
held-out strength in both directions is named "kept on the matrix, higher on strength" (`step1_strength_names`; no
refit). The outcomes (a) to (f), "unit-dependent" and "reading-dependent" (`outcome`).

Concentrations: megamat throughout, at 10^-3 on both sides. The APL16 input's paa, euc and -aPine are the ORN data
within 1 log10 unit of the flies' concentrations; the imposed arm does not touch them, a kept rung's terms are applied
to them, and every APL16 number involves them.

Outputs: <shared>/results/E135_antennal_lobe_dendrites/ (compact: step0_*.csv, step0_summary.json, timing_step1.json,
lock.json, step1_*.csv, step1_control_<i>.json, runs_2_<host>_<lane>.csv, provenance_2_<host>_<lane>.json,
identity_2_<host>_<lane>.json, claims/, score_scores.csv, score_class_rows.csv, score_drops.csv, score_stage_error.csv,
score_imposed_checks.csv, score_realism.csv, score_apl_effect.csv, score_verdict.json, score_table.md,
provenance_*.json); outputs/E135_antennal_lobe_dendrites/step2/ (local, large: step 2's runs, also in the archive's
data/model_outputs/E135_antennal_lobe_dendrites/step2/).

Usage (from the repository root, as modules, so that the checkout's kcq is the one imported):
    ./venv/bin/python -u -m experiments.E135_antennal_lobe_dendrites --part step0
    nice -n 10 ./venv/bin/python -u -m experiments.E135_antennal_lobe_dendrites --part timing
    ./venv/bin/python -u -m experiments.E135_antennal_lobe_dendrites --part lock
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E135_antennal_lobe_dendrites --part step1 | --part step1_control --perm i | --part step2 \\
        [--lane <host>[-<n>]] [--reverse]
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E135_antennal_lobe_dendrites --part score [--jobs <n>]
    ./venv/bin/python -u -m experiments.E135_antennal_lobe_dendrites --part commands|variants
"""
from __future__ import annotations

import argparse
import dataclasses
import gc
import hashlib
import itertools
import json
import os
import socket
import subprocess
import time
import warnings
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, stats

from kcq import harness as H
from kcq import population_fit as P

EXPERIMENT = 'E135_antennal_lobe_dendrites'
SCRIPT = 'experiments/E135_antennal_lobe_dendrites.py'
LOCK_COPY = H.HERE / 'experiments' / 'E135_lock.json'
E134 = 'E134_population_baseline'
# E134's verdict files, read before any fit or model run: the measure gate's verdict (written by E134's owner when the
# gate is read; {"verdict": "usable", ...}) and the score of its megamat runs (atlas's `--part score`)
E134_GATE_FILE = 'measure_gate.json'
E134_SCORED_FILE = 'score_classes.json'
E134_STAGE_TABLE = 'score_stage_table.csv'

ORN_KEY = ('orn', 'megamat', 'wild_type')
PN_KEY = ('pn_dendrite', 'megamat', 'wild_type')
MEGAMAT_ODORS: Tuple[str, ...] = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms',
                                  'B-cit', 'Lin', 'aa', 'benz', 'va')
N_SHARED = 28
NAMED_GLOMERULI: Tuple[str, ...] = ('DM1', 'DL5', 'VA2')
WEAK_INPUT_ODORS: Tuple[str, ...] = ('paa', 'B-myr', '2-but')    # only 2-but has a dendrite recording
# the units a stage's recordings are read in: dF/F as recorded, or the lock's trial-noise units
UNITS: Tuple[str, ...] = ('dff', 'noise')
# the readings of the measured step, as (ORN terminals' units, PN dendrites' units) (ruling 1): dF/F on both sides
# (the first look's), noise units on both sides, and the rungs' reading, the terminals in dF/F (the model's input is
# their dF/F times one factor) and the dendrites in the fit's units
STEP_READINGS: Dict[str, Tuple[str, str]] = {'dff': ('dff', 'dff'), 'noise': ('noise', 'noise'),
                                             'orn_dff_pn_noise': ('dff', 'noise')}
STEP_UNITS: Tuple[str, ...] = tuple(STEP_READINGS)
RUNG_STEP_READING = 'orn_dff_pn_noise'
# the units the fits, the held-out reads and the imposed table read the dendrites in (the lock's trial-noise units),
# and the units of the ORN side of every fit's measured step
FIT_UNITS = 'noise'
FIT_ORN_UNITS = 'dff'

GLOBAL_PARAMS: Tuple[str, ...] = ('pn_inhsc', 'pn_inhadd', 'ln_inhsc', 'ln_inhadd', 'pn_tanhsc')
OLFSYSM_DEFAULTS: Dict[str, float] = {'pn_inhsc': 368.6631, 'pn_inhadd': 31.4088, 'ln_inhsc': 500.0, 'ln_inhadd': 200.0,
                                      'pn_tanhsc': 5.3395}
# R1's bounds, as multiples of the defaults (log-spaced search); pn_inhadd up to the ladder's top (10^6, the gain fixed)
GLOBAL_BOUNDS: Dict[str, Tuple[float, float]] = {'pn_inhsc': (0.05, 20.0), 'pn_inhadd': (0.03, 3.2e4),
                                                 'ln_inhsc': (0.05, 20.0), 'ln_inhadd': (0.05, 20.0),
                                                 'pn_tanhsc': (0.05, 20.0)}
GAIN_BOUNDS: Tuple[float, float] = (0.25, 4.0)
LATERAL_BOUNDS: Tuple[float, float] = (-1.0, 1.0)
MAX_LATERAL = 2
N_CANDIDATES = 2
KEEP_SE_MULT = 2.0
RUNGS: Tuple[str, ...] = ('R0', 'R1', 'R2', 'R3')
SOLVE_TOL_HZ = 0.5
N_BOOT = 1000
N_BOOT_KEEP = 200
N_PERM_PARTNERS = 1000
SEED = 0
CI = 0.9
# folds and directions (ruling 3): one odor per fold, each half fitted and read on the other
N_FOLDS = 17
N_FOLDS_FALLBACK = 6
DIRECTIONS: Tuple[Tuple[str, str], ...] = (('half_a', 'half_b'), ('half_b', 'half_a'))
# the reading every kept decision and reported error uses (the lock's default for the PN dendrites) and the one read
# beside it (ruling 4)
KEEP_READING = 'curve_noise'
BESIDE_READING = 'affine_noise'
# permutations (ruling 6): a lateral term must gain more than at least 19 of 20 fits of the same term to odor-permuted
# dendrites; R1's and R2's refits to odor-permuted dendrites (5, or 3 on cost) are reported and decide nothing
N_PERM_LATERAL = 20
LATERAL_PERM_MIN_EXCEEDED = 19
N_PERM_CONTROL = 5
N_PERM_CONTROL_FALLBACK = 3
CONTROL_SEED_TAG = 1
MAX_LATERAL_TERMS = 4
# the tables' scale (ruling 5): a factor outside this range is named
SCALE_FACTOR_RANGE: Tuple[float, float] = (0.5, 2.0)
SCALE_READING = 'preset'
# rates outside a PN's reach (ruling 7): clipped this far inside `pn_rate_bounds` (well inside SOLVE_TOL_HZ, and
# enough that the silenced floor is reached by a finite ORN delta inside route 2's bracket); an arm with more than
# CLIP_NAME_SHARE of its set entries clipped is named
CLIP_POLICY = 'clip'
BOUNDS_MARGIN_HZ = 0.1
CLIP_NAME_SHARE = 0.10
# the imposed runs' check against E134's reading (ruling 8)
IMPOSED_CHECK_HZ = 0.5
PN_MATRIX_FLOOR_NOTE = 0.31
# step 1's budget (the rulings'): lane-hours, and the fallback on cost only, (controls, folds) in order
BUDGET_LANE_HOURS = 8.0
BUDGET_STEPS: Tuple[Tuple[int, int], ...] = ((N_PERM_CONTROL, N_FOLDS), (N_PERM_CONTROL_FALLBACK, N_FOLDS),
                                            (N_PERM_CONTROL_FALLBACK, N_FOLDS_FALLBACK))
# the real size of step 1's fits, for the timing helper: the family's PN stage has the wiring's 54 glomeruli, the 36 of
# the lab input with input and 18 without (ruling 9's "18 glomeruli without input"); 28 shared; 17 odors; 7 PN flies
# of 3 trials (the halves 3 and 4)
REAL_SIZE: Dict[str, int] = dict(n_glomeruli=54, n_input=36, n_shared=N_SHARED, n_odors=17, n_pn_flies=7, n_trials=3)

# the stage windows (s), as kcq/stage_readings.py reads the model's PN dendrites
MODEL_WINDOW: Tuple[float, float] = (0.0, 0.5)
SETTLED_PRE_WINDOW: Tuple[float, float] = (-0.25, 0.0)

# step 2: the arms (input panel, arm); 8 wirings each, plus one repeat of an imposed run by the second route
Arm = namedtuple('Arm', 'panel arm')
MEGAMAT_ARMS: Tuple[str, ...] = ('R1', 'rung', 'rung_glom_perm', 'imposed', 'imposed_glom_perm', 'imposed_odor_perm')
APL16_ARMS: Tuple[str, ...] = ('R1', 'rung', 'imposed')
RUNG_ARMS: Tuple[str, ...] = ('rung', 'rung_glom_perm')
REPEAT_ARM = Arm('megamat', 'imposed_route2_repeat')
# the dF/F repeat of the three imposed megamat arms, before outcome (e) or (f) is declared (ruling 1's guard)
DFF_SUFFIX = '_dff'
IMPOSED_ARMS: Tuple[str, ...] = ('imposed', 'imposed_glom_perm', 'imposed_odor_perm')
DFF_REPEAT_ARMS: Tuple[str, ...] = tuple(a + DFF_SUFFIX for a in IMPOSED_ARMS)
WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1', 'P3', 'P4', 'M2', 'M3')    # E134's S then T
# the stages whose population errors the classes and the routes' agreement read
CLASS_STAGES: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')
E134_HELD_GROUP = 'megamat_held'

# the classes (the entry's; the thresholds are reset from E134's table before any arm runs)
MIN_WIRINGS = 6
SOMA_DROP_UPSTREAM = 0.15
SOMA_WITHIN_BOUTONS_ONLY = 0.05
BOUTON_STRENGTH_HALVED = 0.5
CARRY_SHARE = 0.5
THRESHOLD_SE_MULT = 3.0
# the forms of a stage's error the megamat stage error is given in (E134's gate of 2026-11-24 15:19 read the three-view
# population error not usable; Tom, about 15:21, and Yang, 15:28, adopted strength by odor alone, with sharing, breadth
# and the matrix view beside and the two-view error tested beside; measure_gate.json 'measure'), as E136's MEASURES:
# 'all' the population error (three views; the PN dendrites on four), else the root mean square of the named views'
# errors
MEASURES: Dict[str, Optional[Tuple[str, ...]]] = {
    'all': None,
    'strength': ('strength',),
    'strength_sharing': ('strength', 'sharing'),
}
# the classes' column (the rulings of 15:41, from a Fable planner, claude-fable-5-2, spawned from hal 2026-11-24 15:32):
# strength by odor, the adopted measure, at the somata as at the boutons, in each stage's default reading; the
# thresholds are frozen before step 2's first arm (THRESHOLDS_FILE, `freeze_thresholds`) from E134's stage table and its
# step-2 scores, and never re-read
CLASS_COLUMN = 'view_strength'
THRESHOLDS_FILE = 'step2_thresholds.json'
E134_STEP2_SCORES = 'step2_scores.csv'
E134_APL16_SOURCE = 'E90:B'
THRESHOLD_STAGES: Tuple[str, ...] = ('pn_dendrite',) + CLASS_STAGES
# the classes' second reading (the rulings): the somata on raw counts, the boutons in their first alternative; a class
# that differs is "reading-dependent"
SECOND_READINGS: Dict[str, str] = {'kc_soma': 'raw', 'pn_bouton': 'curve_noise_g3'}
# rule 3's guard (the rulings): for R1 and the kept rung against the family's run of the same wiring, a rise of the
# strength-by-odor error or of the population error by more than 3 standard errors on six of eight wirings
GUARD_ARMS: Tuple[str, ...] = ('R1', 'rung')
GUARD_STAGES: Tuple[str, ...] = CLASS_STAGES
GUARD_COLUMNS: Tuple[str, ...] = ('view_strength', 'error')
BOUTON_RATE_CAP_HZ = 300.0
BOUTON_RATE_CAP_SHARE = 0.01
# the lanes: threads per lane by host, and the lanes of each host (a lane's name, reversed or not); disco (a lab
# member's machine, 16 cores) and atlas take half their cores
LANES: Dict[str, Tuple[int, Tuple[Tuple[str, bool], ...]]] = {
    'viki': (16, (('viki', False),)),
    'hal': (16, (('hal-1', False), ('hal-2', True))),
    'atlas': (8, (('atlas', True),)),
    'disco': (8, (('disco', False),)),
}
MODULE = 'experiments.E135_antennal_lobe_dendrites'
# step 1's lanes (Tom, 2026-11-24 about 15:21: "fit the PN stage to the dendrites now on hal's and disco's lanes"): the
# real rungs on the first, the controls on the others in turn
STEP1_LANES: Tuple[str, ...] = ('hal-1', 'hal-2', 'disco')


# ---------------------------------------------------------------------------
# paths, small helpers and guards

def out_dir() -> Path:
    """<shared>/results/E135_antennal_lobe_dendrites/ (created)."""
    return H.results_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _jsonable(x):
    if hasattr(x, 'item'):
        return x.item()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=_jsonable) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text()) if Path(path).exists() else {}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def _nanmean(A: np.ndarray, axis: int = 0) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(A, axis=axis)


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by this repository's git with no uncommitted change."""
    try:
        r = subprocess.run(['git', '-C', str(H.HERE), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(H.HERE), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


def require_e134_lock() -> dict:
    """E134's committed lock (atlas's `E134_population_baseline.require_lock`, which runs viki's check); every part reads
    the recordings through it."""
    from experiments import E134_population_baseline as E134R
    return E134R.require_lock()


def e134_results_dir() -> Path:
    """<shared>/results/E134_population_baseline/: E134's gate, its stage table and its stage errors."""
    return H.shared_dir() / 'results' / E134


def e134_gate_status(results: Optional[Path] = None) -> dict:
    """What E134 has recorded of its measure gate: the verdict file's contents (E134_GATE_FILE) and whether its megamat
    runs are scored (E134_SCORED_FILE), in <shared>/results/E134_population_baseline/ (or `results`)."""
    d = Path(results) if results is not None else e134_results_dir()
    gate = _read_json(d / E134_GATE_FILE)
    return dict(dir=str(d), gate_file=E134_GATE_FILE, gate=gate, verdict=gate.get('verdict'),
                scored=(d / E134_SCORED_FILE).exists(), scored_file=E134_SCORED_FILE)


def require_e134_gate(results: Optional[Path] = None) -> dict:
    """E134's measure gate, read before any fit or model run: refuses (PermissionError, naming what it waits for) until
    E134 has recorded that its measure is "usable" and its megamat runs are scored."""
    st = e134_gate_status(results)
    if st['verdict'] != 'usable' or not st['scored']:
        raise PermissionError(
            "E135's fits and model runs wait on E134's measure gate (the entry's Design 1 and Status): they need "
            f"{st['dir']}/{E134_GATE_FILE} with \"verdict\": \"usable\" (now: {st['verdict']!r}), written by E134's "
            f"owner once the gate is read, and {st['dir']}/{E134_SCORED_FILE} (E134's megamat runs scored, its runner's "
            f"--part score; now: {'present' if st['scored'] else 'absent'}). A planner then re-reads E135's thresholds. "
            'Step 0 (recordings only) does not wait.')
    return st


# Tom's authorization of step 1 before E134's gate (2026-11-24 about 15:21, in hal's terminal), as hal's session
# records it in results/E135_antennal_lobe_dendrites/: {"authorized_by": "Tom", "when": ..., "words": ...}
STEP1_AUTH_FILE = 'step1_authorized.json'
STEP1_AUTH_KEYS: Tuple[str, ...] = ('authorized_by', 'when', 'words')


def step1_authorization(root: Optional[Path] = None) -> dict:
    """Tom's authorization of step 1 before E134's gate (results/.../STEP1_AUTH_FILE, or in `root`): its contents
    where it is Tom's and has every key of STEP1_AUTH_KEYS non-empty, else {}."""
    a = _read_json((out_dir() if root is None else Path(root)) / STEP1_AUTH_FILE)
    if a.get('authorized_by') != 'Tom' or not all(a.get(k) for k in STEP1_AUTH_KEYS):
        return {}
    return a


def require_step1_allowed(results: Optional[Path] = None, root: Optional[Path] = None) -> dict:
    """Step 1 may run (its fits, its controls and the family's build they read the PN stage from): E134's measure gate
    reads "usable" with its megamat runs scored (`require_e134_gate`, `results` its folder), or Tom has authorized
    step 1 before the gate (`step1_authorization`, `root` its folder; Tom, 2026-11-24 about 15:21: "Run step 1 now: fit
    the PN stage to the dendrites now on hal's and disco's lanes. Step 2's model runs still wait for the adopted
    measure."). PermissionError naming both otherwise. Step 2 and the score do not use this: they keep
    `require_e134_gate`.

    Returns:
        dict with 'how' ('gate' or 'authorized'), the gate's status and, where authorized, the authorization.
    """
    try:
        return dict(how='gate', gate=require_e134_gate(results))
    except PermissionError as e:
        auth = step1_authorization(root)
        if auth:
            return dict(how='authorized', authorization=auth, gate=e134_gate_status(results))
        where = (out_dir() if root is None else Path(root)) / STEP1_AUTH_FILE
        raise PermissionError(f'{e} Step 1 may also run on Tom\'s authorization before the gate: {where} with '
                              '{"authorized_by": "Tom", "when": ..., "words": ...} (written by hal\'s session); it is '
                              'absent or incomplete.') from e


# ---------------------------------------------------------------------------
# step 0: the recorded step (recordings only)

def shared_glomeruli(check: bool = True) -> List[str]:
    """The 28 glomeruli the ORN recordings, the PN dendrites and the model's input share (E90's `shared_glomeruli`),
    checked against the lock's populations (a glomerulus any fly of each stage has) when `check`."""
    from experiments import E90_model_stages as E90
    gl = list(E90.shared_glomeruli())
    if check:
        from data_prep import stage_population_data as D
        require_e134_lock()
        have = []
        for key in (ORN_KEY, PN_KEY):
            pop = D.load_population(*key)
            have.append(set().union(*[set(map(str, t.index)) for t in pop.tables.values()]))
        model = set(H.apl16_orn_deltas().index)
        alt = sorted(have[0] & have[1] & model)
        if alt != sorted(gl):
            raise RuntimeError(f"E90's shared glomeruli {sorted(gl)} differ from the lock's populations' {alt}")
    if len(gl) != N_SHARED:
        raise RuntimeError(f'{len(gl)} shared glomeruli, not {N_SHARED}')
    return sorted(gl)


def recorded_stage(key: Tuple[str, str, str], units: str = 'dff', tol: float = 1e-9) -> Dict[str, pd.DataFrame]:
    """Per fly, units x the 17 megamat odors of trial means, in dF/F (`units='dff'`) or in the lock's trial-noise
    units ('noise'), read through the data module against E134's lock.

    The population file (sha256 checked against the lock) holds the trial-noise values only; dF/F is the same rows
    times their noise SD, from the data module's own reader of the population's source (`orn_long` / `pn_dendrite_long`),
    whose trial means in noise units must equal the population's to `tol`."""
    from data_prep import stage_population_data as D
    if units not in UNITS:
        raise ValueError(f'{units=}')
    require_e134_lock()
    pop = D.load_population(*key)
    if key == ORN_KEY:
        long = D.orn_long(D.OLFDB_DATASETS['orn_megamat'], 'megamat')
    elif key == PN_KEY:
        long = D.pn_dendrite_long()
    else:
        raise KeyError(key)
    long = long.assign(dff=long.value * long.noise_sd)
    out = {}
    for fly in pop.flies:
        sub = long[long.fly == fly]
        t = pop.tables[fly]
        g = sub.groupby(['unit', 'odor'], sort=False)
        noise = g.value.mean().unstack('odor').reindex(index=t.index, columns=list(pop.odors))
        a, b = noise.to_numpy(dtype=float), t.to_numpy(dtype=float)
        if not (np.array_equal(np.isnan(a), np.isnan(b)) and np.nanmax(np.abs(a - b)) <= tol):
            raise RuntimeError(f'{key} {fly}: the reader differs from the locked population')
        x = noise if units == 'noise' else g.dff.mean().unstack('odor').reindex(index=t.index, columns=list(pop.odors))
        x.index = x.index.astype(str)
        out[str(fly)] = x.astype(float)
    return out


def fly_stack(tables: Mapping[str, pd.DataFrame], glomeruli: Sequence[str], odors: Sequence[str]) -> np.ndarray:
    """flies x glomeruli x odors (NaN where a fly lacks a glomerulus or an odor)."""
    return np.stack([t.reindex(index=list(glomeruli), columns=list(odors)).to_numpy(dtype=float)
                     for t in tables.values()])


def _slopes(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per row: the slope through the origin of y on x, Pearson r and the number of jointly finite entries."""
    ok = np.isfinite(x) & np.isfinite(y)
    n = ok.sum(axis=1)
    x0, y0 = np.where(ok, x, 0.0), np.where(ok, y, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        slope = (x0 * y0).sum(axis=1) / (x0 * x0).sum(axis=1)
        mx, my = x0.sum(axis=1) / n, y0.sum(axis=1) / n
        dx, dy = np.where(ok, x - mx[:, None], 0.0), np.where(ok, y - my[:, None], 0.0)
        r = (dx * dy).sum(axis=1) / np.sqrt((dx * dx).sum(axis=1) * (dy * dy).sum(axis=1))
    slope = np.where(n >= 3, slope, np.nan)
    r = np.where(n >= 3, r, np.nan)
    return slope, r, n


def _relative(s: np.ndarray) -> np.ndarray:
    m = np.nanmedian(s[np.isfinite(s) & (s > 0)]) if np.any(np.isfinite(s) & (s > 0)) else np.nan
    return s / m


def measured_step(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                  total: Optional[np.ndarray] = None) -> pd.DataFrame:
    """The measured step on fly means.

    Args:
        orn, pn: flies x glomeruli x odors of the ORN terminals and the PN dendrites (the same glomeruli and odors).
        glomeruli, odors: their labels.
        total: the odor's summed terminal response (odors,), for the regression on own and summed response.

    Returns:
        Per glomerulus: slope (through the origin, PN on own ORN), gain (the slope over the median glomerulus's), r
        (Pearson over odors), n_odors; with `total` also own_coef, total_coef, intercept (least squares of PN on own
        ORN, the summed ORN and 1), r2_own (r squared) and r2_own_total.
    """
    x, y = _nanmean(orn), _nanmean(pn)
    slope, r, n = _slopes(x, y)
    df = pd.DataFrame(dict(slope=slope, gain=_relative(slope), r=r, n_odors=n), index=pd.Index(list(glomeruli),
                      name='glomerulus'))
    if total is not None:
        rows = []
        for i in range(len(glomeruli)):
            ok = np.isfinite(x[i]) & np.isfinite(y[i]) & np.isfinite(total)
            if ok.sum() < 4:
                rows.append((np.nan,) * 4)
                continue
            A = np.column_stack([x[i, ok], total[ok], np.ones(ok.sum())])
            coef, *_ = np.linalg.lstsq(A, y[i, ok], rcond=None)
            res = y[i, ok] - A @ coef
            ss = ((y[i, ok] - y[i, ok].mean()) ** 2).sum()
            rows.append((coef[0], coef[1], coef[2], 1.0 - (res ** 2).sum() / ss if ss > 0 else np.nan))
        R = np.array(rows)
        df['own_coef'], df['total_coef'], df['intercept'], df['r2_own_total'] = R.T
        df['r2_own'] = df.r ** 2
    return df


def residuals(orn: np.ndarray, pn: np.ndarray, step: pd.DataFrame, glomeruli: Sequence[str],
              odors: Sequence[str]) -> pd.DataFrame:
    """The residual by odor of each glomerulus's dendrites (fly mean) about its own slope (glomeruli x odors)."""
    x, y = _nanmean(orn), _nanmean(pn)
    s = step.loc[list(glomeruli), 'slope'].to_numpy(dtype=float)
    return pd.DataFrame(y - s[:, None] * x, index=pd.Index(list(glomeruli), name='glomerulus'), columns=list(odors))


def _corr_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Pearson r between every row of A and every row of B (no NaN), A.shape[0] x B.shape[0]."""
    a = A - A.mean(axis=1, keepdims=True)
    b = B - B.mean(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (a @ b.T) / np.sqrt((a * a).sum(axis=1)[:, None] * (b * b).sum(axis=1)[None, :])


def candidate_partners(resid: pd.DataFrame, orn_mean: pd.DataFrame, own: Optional[pd.DataFrame] = None,
                       dendrites: Optional[pd.DataFrame] = None, n: int = N_CANDIDATES,
                       n_perm: int = N_PERM_PARTNERS, seed: int = SEED) -> pd.DataFrame:
    """Per glomerulus, the other glomeruli whose terminal responses its residual follows best.

    Args:
        resid: glomeruli x odors, the residuals (`residuals`).
        orn_mean: partners x odors, the terminal responses (fly means) of the candidate partners (the model input's
            glomeruli the ORN recordings have); a glomerulus is never its own partner.
        own, dendrites: glomeruli x odors of the own terminal responses and the dendrites (fly means), for each pair's
            joint least-squares weights (dendrites on own and partner, through the origin).
        n: partners kept per glomerulus.
        n_perm: permutations of the residual's odor labels for the chance level of the best partner's r squared.
        seed: of the permutations.

    Returns:
        Rows glomerulus, rank (1 = best), partner, r, r2, chance_r2_95 (the 95th percentile of the best partner's r
        squared under permuted odor labels), and with `own` and `dendrites` own_weight, partner_weight and
        weight_rel (partner over own weight).
    """
    odors = list(resid.columns)
    X = orn_mean.reindex(columns=odors).to_numpy(dtype=float)
    partners = list(orn_mean.index.astype(str))
    rng = np.random.default_rng(seed)
    rows = []
    for g in resid.index:
        e = resid.loc[g].to_numpy(dtype=float)
        ok = np.isfinite(e) & np.isfinite(X).all(axis=0)
        cand = [i for i, p in enumerate(partners) if p != str(g)]
        if ok.sum() < 4 or not cand:
            continue
        Xc = X[np.ix_(cand, np.flatnonzero(ok))]
        r = _corr_rows(e[ok][None, :], Xc)[0]
        if n_perm > 0:
            perm = np.array([rng.permutation(e[ok]) for _ in range(n_perm)])
            chance = float(np.nanpercentile(np.nanmax(_corr_rows(perm, Xc) ** 2, axis=1), 95))
        else:
            chance = np.nan
        order = np.argsort(-(r ** 2))
        for k, j in enumerate(order[:n]):
            row = dict(glomerulus=str(g), rank=k + 1, partner=partners[cand[j]], r=float(r[j]), r2=float(r[j] ** 2),
                       chance_r2_95=chance)
            if own is not None and dendrites is not None:
                xo = own.loc[g, odors].to_numpy(dtype=float)
                y = dendrites.loc[g, odors].to_numpy(dtype=float)
                ok2 = ok & np.isfinite(xo) & np.isfinite(y)
                A = np.column_stack([xo[ok2], X[cand[j], ok2]])
                w, *_ = np.linalg.lstsq(A, y[ok2], rcond=None)
                row.update(own_weight=float(w[0]), partner_weight=float(w[1]),
                           weight_rel=float(w[1] / w[0]) if w[0] != 0 else np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def _total(orn_all: np.ndarray) -> np.ndarray:
    """The odor's summed terminal response over glomeruli (fly means; glomeruli no fly has are left out)."""
    return np.nansum(_nanmean(orn_all), axis=0)


def step_intervals(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                   n_boot: int = N_BOOT, seed: int = SEED, ci: float = CI) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Intervals from resampling the flies of both stages (each with replacement, independently), centered on the
    point values as the analysis of 2026-11-17 centers them (`centered`: the point minus the draws' mean less their
    lower quantile, the point plus their upper quantile less their mean): a resample with repeated flies has noisier fly
    means, so its r sits below the measured one and raw percentiles would not hold the point.

    Returns:
        (per glomerulus gain_lo, gain_hi, r_lo, r_hi, slope_lo, slope_hi, gain_boot_sd, log2_gain_boot_sd;
        per glomerulus x odor resid_lo, resid_hi), at coverage `ci`.
    """
    rng = np.random.default_rng(seed)
    G, O = len(glomeruli), len(odors)
    gains, rs, sl = np.full((n_boot, G), np.nan), np.full((n_boot, G), np.nan), np.full((n_boot, G), np.nan)
    res = np.full((n_boot, G, O), np.nan)
    for b in range(n_boot):
        io = rng.integers(0, orn.shape[0], orn.shape[0])
        ip = rng.integers(0, pn.shape[0], pn.shape[0])
        x, y = _nanmean(orn[io]), _nanmean(pn[ip])
        s, r, _ = _slopes(x, y)
        gains[b], rs[b], sl[b] = _relative(s), r, s
        res[b] = y - s[:, None] * x
    q = ((1 - ci) / 2, 1 - (1 - ci) / 2)
    x0, y0 = _nanmean(orn), _nanmean(pn)
    s0, r0, _ = _slopes(x0, y0)
    point = dict(gain=_relative(s0), r=r0, slope=s0)
    res0 = y0 - s0[:, None] * x0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        cols = {}
        for k, B in (('gain', gains), ('r', rs), ('slope', sl)):
            cols[f'{k}_lo'], cols[f'{k}_hi'] = centered(point[k], B, q)
        cols['r_lo'], cols['r_hi'] = np.clip(cols['r_lo'], -1.0, 1.0), np.clip(cols['r_hi'], -1.0, 1.0)
        g = pd.DataFrame(dict(cols, gain_boot_sd=np.nanstd(gains, axis=0),
                              log2_gain_boot_sd=np.nanstd(np.log2(np.where(gains > 0, gains, np.nan)), axis=0)),
                         index=pd.Index(list(glomeruli), name='glomerulus'))
        lo, hi = centered(res0, res, q)
    e = pd.DataFrame(dict(glomerulus=np.repeat(list(glomeruli), O), odor=np.tile(list(odors), G),
                          resid_lo=lo.ravel(), resid_hi=hi.ravel()))
    return g, e


def centered(point: np.ndarray, draws: np.ndarray, q: Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]:
    """The interval of the draws (axis 0) centered on the point: point - (mean - lower quantile), point + (upper
    quantile - mean)."""
    m = np.nanmean(draws, axis=0)
    return point - (m - np.nanquantile(draws, q[0], axis=0)), point + (np.nanquantile(draws, q[1], axis=0) - m)


def halves(pn_flies: Sequence[str], seed: int = SEED) -> Tuple[List[str], List[str]]:
    """The two halves of the PN flies (sizes floor(n/2) and the rest; hal's `split_halves` with one split), each in the
    flies' order."""
    a, b = P.split_halves(len(pn_flies), n_splits=1, seed=seed)[0]
    return [pn_flies[i] for i in a], [pn_flies[i] for i in b]


def halves_agreement(sa: pd.DataFrame, sb: pd.DataFrame) -> Dict[str, float]:
    """The agreement of two halves' measured steps: Pearson r and Spearman ρ of log2 gains over glomeruli, the RMS
    difference of log2 gains, and Pearson r of the halves' r over glomeruli."""
    ga, gb = np.log2(sa.gain.where(sa.gain > 0)), np.log2(sb.gain.where(sb.gain > 0))
    ok = ga.notna() & gb.notna()
    ra, rb = sa.r, sb.r
    okr = ra.notna() & rb.notna()
    return dict(log2_gain_pearson=float(np.corrcoef(ga[ok], gb[ok])[0, 1]) if ok.sum() > 2 else np.nan,
                log2_gain_spearman=float(stats.spearmanr(ga[ok], gb[ok]).correlation) if ok.sum() > 2 else np.nan,
                log2_gain_rms_diff=float(np.sqrt(np.mean((ga[ok] - gb[ok]) ** 2))) if ok.any() else np.nan,
                r_pearson=float(np.corrcoef(ra[okr], rb[okr])[0, 1]) if okr.sum() > 2 else np.nan,
                n_glomeruli=int(ok.sum()))


def all_splits_agreement(orn: np.ndarray, pn: np.ndarray, glomeruli: Sequence[str], odors: Sequence[str],
                         total: Optional[np.ndarray] = None) -> pd.DataFrame:
    """`halves_agreement` over every split of the PN flies into floor(n/2) and the rest (35 for 7 flies)."""
    n = pn.shape[0]
    rows = []
    for a in itertools.combinations(range(n), n // 2):
        b = [i for i in range(n) if i not in a]
        sa = measured_step(orn, pn[list(a)], glomeruli, odors)
        sb = measured_step(orn, pn[b], glomeruli, odors)
        rows.append(dict(half_a=list(a), **halves_agreement(sa, sb)))
    return pd.DataFrame(rows)


def follows_rule(step: pd.DataFrame) -> pd.Series:
    """A glomerulus's dendrites follow its terminals where the 90 % interval of its slope (dendrite on terminal, over
    resampled flies; `step_intervals`' slope_lo) lies above 0 (ruling 2). Where they do not, its measured gain is left
    at 1 (R2m; R2f and R3 may still fit it). Applied in each fit on the flies and odors that fit sees."""
    return step['slope_lo'] > 0


def unit_dependence(step: pd.DataFrame, readings: Sequence[str] = STEP_UNITS) -> pd.DataFrame:
    """Each glomerulus's measured gain (all flies) in each reading of the measured step, the side of 1 it is on, and
    whether that side differs between the readings ("unit-dependent", ruling 1).

    Args:
        step: step 0's measured step (rows with columns units, flies, glomerulus, gain, gain_lo, gain_hi, rank_gain).
        readings: the readings compared (STEP_READINGS' names).

    Returns:
        One row per glomerulus: gain_<reading>, gain_lo_<reading>, gain_hi_<reading>, rank_<reading>, side_<reading>
        ('above', 'below' or 'at' 1), unit_dependent (the sides differ) and interval_dependent (some reading's interval
        lies wholly above 1 and another's wholly below).
    """
    s = step[step.flies == 'all']
    out = None
    for r in readings:
        x = s[s.units == r].set_index('glomerulus')
        t = pd.DataFrame({f'gain_{r}': x.gain, f'gain_lo_{r}': x.gain_lo, f'gain_hi_{r}': x.gain_hi,
                          f'rank_{r}': x.rank_gain,
                          f'side_{r}': np.where(x.gain > 1, 'above', np.where(x.gain < 1, 'below', 'at'))},
                         index=x.index)
        out = t if out is None else out.join(t, how='outer')
    sides = out[[f'side_{r}' for r in readings]]
    out['unit_dependent'] = sides.nunique(axis=1) > 1
    above = pd.concat([out[f'gain_lo_{r}'] > 1 for r in readings], axis=1).any(axis=1)
    below = pd.concat([out[f'gain_hi_{r}'] < 1 for r in readings], axis=1).any(axis=1)
    out['interval_dependent'] = above & below
    return out.reset_index()


def no_step_rows() -> pd.DataFrame:
    """E134's no-step yardstick and recorded-step summary of the ORN -> PN dendrite step (viki's step 1, the flies
    alone): flies_no_step.csv and flies_steps_summary.csv rows of that step."""
    d = H.shared_dir() / 'results' / E134
    step = f'{"__".join(ORN_KEY)} -> {"__".join(PN_KEY)}'
    out = []
    for name in ('flies_no_step.csv', 'flies_steps_summary.csv'):
        p = d / name
        if not p.exists():
            continue
        t = pd.read_csv(p)
        t = t[t.step == step].copy()
        t.insert(0, 'source', name)
        out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def read_units(reading: str, glomeruli: Sequence[str], odors: Sequence[str] = MEGAMAT_ODORS,
               cache: Optional[Dict[tuple, Dict[str, pd.DataFrame]]] = None) -> dict:
    """The recordings of both stages in a reading of the measured step (STEP_READINGS: the ORN terminals' and the PN
    dendrites' units): per-fly tables, the stacks on the shared glomeruli, the ORN stack on the model input's glomeruli
    the ORN recordings have (the partners and the summed response, in the ORN side's units), and the flies. `cache`
    holds `recorded_stage`'s tables by (stage key, units) across readings."""
    orn_units, pn_units = STEP_READINGS[reading]
    cache = {} if cache is None else cache
    for key, u in ((ORN_KEY, orn_units), (PN_KEY, pn_units)):
        if (key, u) not in cache:
            cache[(key, u)] = recorded_stage(key, u)
    orn_t, pn_t = cache[(ORN_KEY, orn_units)], cache[(PN_KEY, pn_units)]
    model_gl = list(H.apl16_orn_deltas().index)
    orn_have = sorted(set().union(*[set(t.index) for t in orn_t.values()]))
    partners = [g for g in model_gl if g in orn_have]
    return dict(orn_tables=orn_t, pn_tables=pn_t, orn=fly_stack(orn_t, glomeruli, odors),
                pn=fly_stack(pn_t, glomeruli, odors), orn_all=fly_stack(orn_t, partners, odors), partners=partners,
                orn_flies=list(orn_t), pn_flies=list(pn_t), orn_units=orn_units, pn_units=pn_units)


def step0_units(units: str, glomeruli: Sequence[str], odors: Sequence[str] = MEGAMAT_ODORS, n_boot: int = N_BOOT,
                seed: int = SEED, data: Optional[dict] = None, cache: Optional[dict] = None) -> dict:
    """Step 0 in one reading of the measured step (`units`, one of STEP_READINGS; `part_step0`'s tables for it)."""
    R = read_units(units, glomeruli, odors, cache) if data is None else data
    orn, pn, orn_all = R['orn'], R['pn'], R['orn_all']
    total = _total(orn_all)
    A, B = halves(R['pn_flies'], seed)
    ia = [R['pn_flies'].index(f) for f in A]
    ib = [R['pn_flies'].index(f) for f in B]
    steps, partners_rows = [], []
    allstep = measured_step(orn, pn, glomeruli, odors, total)
    iv, ev = step_intervals(orn, pn, glomeruli, odors, n_boot, seed)
    allstep = allstep.join(iv)
    allstep['follows'] = follows_rule(allstep)
    allstep['gain_for_r2m'] = allstep.gain.where(allstep.follows, 1.0)
    allstep['gain_interval_excludes_1'] = (allstep.gain_lo > 1) | (allstep.gain_hi < 1)
    allstep['rank_gain'] = allstep.gain.rank(ascending=False, method='min')
    allstep['rank_r'] = allstep.r.rank(ascending=False, method='min')
    sets = {'all': (orn, pn, allstep), 'half_a': (orn, pn[ia], None), 'half_b': (orn, pn[ib], None)}
    orn_mean = pd.DataFrame(_nanmean(orn_all), index=R['partners'], columns=list(odors))
    for name, (o, p, st) in sets.items():
        if st is None:
            st = measured_step(o, p, glomeruli, odors, total)
            hv, _ = step_intervals(o, p, glomeruli, odors, n_boot, seed)
            st = st.join(hv)
            st['follows'] = follows_rule(st)
        res = residuals(o, p, st, glomeruli, odors)
        own = pd.DataFrame(_nanmean(o), index=list(glomeruli), columns=list(odors))
        den = pd.DataFrame(_nanmean(p), index=list(glomeruli), columns=list(odors))
        cp = candidate_partners(res, orn_mean, own, den, N_CANDIDATES, N_PERM_PARTNERS, seed)
        cp.insert(0, 'flies', name)
        partners_rows.append(cp)
        st = st.reset_index()
        st.insert(0, 'flies', name)
        steps.append(st)
        if name == 'all':
            resid_all, own_all, den_all = res, own, den
    step_df = pd.concat(steps, ignore_index=True)
    step_df.insert(0, 'units', units)
    step_df.insert(1, 'orn_units', R['orn_units'])
    step_df.insert(2, 'pn_units', R['pn_units'])
    part_df = pd.concat(partners_rows, ignore_index=True)
    part_df.insert(0, 'units', units)
    rl = resid_all.stack().rename('resid').reset_index()
    rl.columns = ['glomerulus', 'odor', 'resid']
    rl['orn'] = own_all.stack().to_numpy()
    rl['pn'] = den_all.stack().to_numpy()
    rl['fit'] = rl.pn - rl.resid
    rl = rl.merge(ev, on=['glomerulus', 'odor'], how='left')
    rl['departs'] = (rl.resid_lo > 0) | (rl.resid_hi < 0)
    rl.insert(0, 'units', units)
    sa = step_df[step_df.flies == 'half_a'].set_index('glomerulus')
    sb = step_df[step_df.flies == 'half_b'].set_index('glomerulus')
    agree = halves_agreement(sa, sb)
    spl = all_splits_agreement(orn, pn, glomeruli, odors)
    agree_rows = []
    for k in ('log2_gain_pearson', 'log2_gain_spearman', 'log2_gain_rms_diff', 'r_pearson'):
        agree_rows.append(dict(units=units, measure=k, locked_split=agree[k], mean_all_splits=float(spl[k].mean()),
                               min_all_splits=float(spl[k].min()), max_all_splits=float(spl[k].max()),
                               n_splits=int(len(spl)), half_a=';'.join(A), half_b=';'.join(B)))
    return dict(step=step_df, partners=part_df, residuals=rl, agreement=pd.DataFrame(agree_rows),
                halves=dict(half_a=A, half_b=B), n_orn_flies=int(orn.shape[0]), n_pn_flies=int(pn.shape[0]),
                partners_pool=R['partners'], orn_units=R['orn_units'], pn_units=R['pn_units'])


def named_rows(step: pd.DataFrame, resid: pd.DataFrame, partners: pd.DataFrame, units: str,
               names: Sequence[str] = NAMED_GLOMERULI, top: int = 4) -> List[dict]:
    """The named glomeruli's readings (all flies): where each stands among the 28, its largest residuals with
    intervals, and its candidate partners."""
    s = step[(step.units == units) & (step.flies == 'all')].set_index('glomerulus')
    out = []
    for g in names:
        if g not in s.index:
            out.append(dict(glomerulus=g, present=False))
            continue
        r = s.loc[g]
        rr = resid[(resid.units == units) & (resid.glomerulus == g)].copy()
        rr['abs'] = rr.resid.abs()
        rr = rr.sort_values('abs', ascending=False).head(top)
        cp = partners[(partners.units == units) & (partners.flies == 'all') & (partners.glomerulus == g)]
        out.append(dict(glomerulus=g, present=True, units=units, gain=float(r.gain), gain_lo=float(r.gain_lo),
                        gain_hi=float(r.gain_hi), r=float(r.r), r_lo=float(r.r_lo), r_hi=float(r.r_hi),
                        rank_gain=int(r.rank_gain), rank_r=int(r.rank_r), n_glomeruli=int(s.gain.notna().sum()),
                        total_coef=float(r.total_coef), r2_own=float(r.r2_own), r2_own_total=float(r.r2_own_total),
                        follows=bool(r.follows), slope_lo=float(r.slope_lo), slope_hi=float(r.slope_hi),
                        largest_residuals=[dict(odor=o, resid=float(v), lo=float(lo), hi=float(hi), departs=bool(d))
                                           for o, v, lo, hi, d in zip(rr.odor, rr.resid, rr.resid_lo, rr.resid_hi,
                                                                      rr.departs)],
                        partners=[dict(partner=p, r=float(x), r2=float(x2), chance_r2_95=float(c))
                                  for p, x, x2, c in zip(cp.partner, cp.r, cp.r2, cp.chance_r2_95)]))
    return out


def part_step0(units: Sequence[str] = STEP_UNITS, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Step 0: the measured step on all flies and on both halves with their agreement, in each reading of STEP_READINGS
    (dF/F on both sides, noise units on both sides, and the rungs' reading, RUNG_STEP_READING), the named glomeruli
    first, the glomeruli that do not follow their terminals (`follows_rule`), the gains whose side of 1 differs between
    the readings (`unit_dependence`), the entries that depart beyond the flies' spread, E134's no-step yardstick, and
    what cannot be read. Writes step0_measured_step.csv, step0_partners.csv, step0_residuals.csv, step0_halves.csv,
    step0_unit_dependence.csv, step0_no_step.csv, step0_summary.json and provenance_step0.json.

    The saved PN stages of the earlier experiments against the dendrites (the entry's step 0, item 2) are
    `saved_stage_rows`, a separate part (`--part step0_saved`), since it reads saved model outputs."""
    t0 = time.time()
    lock = require_e134_lock()
    glom = shared_glomeruli()
    cache: Dict[tuple, Dict[str, pd.DataFrame]] = {}
    res = {u: step0_units(u, glom, MEGAMAT_ODORS, n_boot, seed, cache=cache) for u in units}
    out = out_dir()
    step = pd.concat([r['step'] for r in res.values()], ignore_index=True)
    partners = pd.concat([r['partners'] for r in res.values()], ignore_index=True)
    resid = pd.concat([r['residuals'] for r in res.values()], ignore_index=True)
    agree = pd.concat([r['agreement'] for r in res.values()], ignore_index=True)
    dep = unit_dependence(step, units)
    nos = no_step_rows()
    step.to_csv(out / 'step0_measured_step.csv', index=False)
    partners.to_csv(out / 'step0_partners.csv', index=False)
    resid.to_csv(out / 'step0_residuals.csv', index=False)
    agree.to_csv(out / 'step0_halves.csv', index=False)
    dep.to_csv(out / 'step0_unit_dependence.csv', index=False)
    nos.to_csv(out / 'step0_no_step.csv', index=False)
    named_dep = dep.set_index('glomerulus').reindex(list(NAMED_GLOMERULI))
    summ = dict(experiment=EXPERIMENT, part='step0', host=_host(), time=_now(), glomeruli=glom,
                odors=list(MEGAMAT_ODORS), units=list(units),
                readings={u: dict(orn=STEP_READINGS[u][0], pn=STEP_READINGS[u][1]) for u in units},
                rung_reading=RUNG_STEP_READING, n_boot=n_boot, seed=seed,
                halves=res[units[0]]['halves'], n_orn_flies=res[units[0]]['n_orn_flies'],
                n_pn_flies=res[units[0]]['n_pn_flies'], partners_pool=res[units[0]]['partners_pool'],
                not_readable=dict(odors=[o for o in WEAK_INPUT_ODORS if o not in MEGAMAT_ODORS],
                                  read=[o for o in WEAK_INPUT_ODORS if o in MEGAMAT_ODORS],
                                  why='the dendrites hold the megamat panel only'),
                rules=dict(follows='the 90 % interval of the slope (dendrite on terminal) over resampled flies lies '
                                   'above 0 (`follows_rule`, ruling 2); otherwise the measured gain is 1 in R2m',
                           unit_dependent='a gain whose side of 1 (its point) differs between the readings (ruling 1); '
                                          'interval_dependent: one reading\'s interval wholly above 1 and another\'s '
                                          'wholly below',
                           total='the summed terminal response over the model input glomeruli the ORN recordings have',
                           flies_noise='fly means (the noise between flies averaged out)'),
                unit_dependent=sorted(dep.glomerulus[dep.unit_dependent.astype(bool)]),
                interval_dependent=sorted(dep.glomerulus[dep.interval_dependent.astype(bool)]),
                named_by_reading={g: {u: dict(gain=float(named_dep.loc[g, f'gain_{u}']),
                                              gain_lo=float(named_dep.loc[g, f'gain_lo_{u}']),
                                              gain_hi=float(named_dep.loc[g, f'gain_hi_{u}']),
                                              rank=int(named_dep.loc[g, f'rank_{u}']),
                                              side=str(named_dep.loc[g, f'side_{u}'])) for u in units}
                                  | dict(unit_dependent=bool(named_dep.loc[g, 'unit_dependent']))
                                  for g in NAMED_GLOMERULI if g in named_dep.index
                                  and np.isfinite(named_dep.loc[g, f'gain_{units[0]}'])})
    for u in units:
        s = step[(step.units == u) & (step.flies == 'all')].set_index('glomerulus')
        r = resid[resid.units == u]
        summ[u] = dict(
            orn_units=STEP_READINGS[u][0], pn_units=STEP_READINGS[u][1],
            named=named_rows(step, resid, partners, u),
            not_following=sorted(s.index[~s.follows.astype(bool)]),
            not_following_halves={h: sorted(step[(step.units == u) & (step.flies == h) & ~step.follows.astype(bool)]
                                            .glomerulus) for h in ('half_a', 'half_b')},
            gain_departs_from_median=sorted(s.index[s.gain_interval_excludes_1.astype(bool)]),
            entries_departing=int(r.departs.sum()), n_entries=int(r.resid.notna().sum()),
            departing_by_glomerulus=r[r.departs].groupby('glomerulus').odor.apply(list).to_dict(),
            gain_range=[float(s.gain.min()), float(s.gain.max())],
            partners_above_chance=sorted(partners[(partners.units == u) & (partners.flies == 'all')
                                                  & (partners['rank'] == 1)
                                                  & (partners.r2 > partners.chance_r2_95)].glomerulus),
            median_r=float(s.r.median()),
            agreement=agree[agree.units == u].set_index('measure').locked_split.to_dict(),
            two_but=r[r.odor == '2-but'][['glomerulus', 'resid', 'resid_lo', 'resid_hi', 'departs']].pipe(
                lambda d: d[d.departs]).to_dict('records'))
    _write_json(out / 'step0_summary.json', summ)
    from data_prep import stage_population_data as D
    pops = lock['populations']
    prov = dict(script=H.script_path(__file__), part='step0', host=_host(), time=_now(), seconds=time.time() - t0,
                **H.code_commits(),
                inputs=dict(lock=dict(path=str(LOCK_E134()), sha256=_sha256(LOCK_E134())),
                            populations={D.pop_name(k): pops[D.pop_name(k)]['file'] for k in (ORN_KEY, PN_KEY)},
                            sources={str(p): _sha256(p) for p in (D.olfdb_path(D.OLFDB_DATASETS['orn_megamat']),
                                                                  D.gh146_copy('GH146'))},
                            e134_flies={n: _sha256(H.shared_dir() / 'results' / E134 / n)
                                        for n in ('flies_no_step.csv', 'flies_steps_summary.csv')
                                        if (H.shared_dir() / 'results' / E134 / n).exists()}),
                readings={u: STEP_READINGS[u] for u in units}, n_boot=n_boot, seed=seed)
    _write_json(out / 'provenance_step0.json', prov)
    print(f'step 0: wrote {out} ({time.time() - t0:.0f} s)', flush=True)
    return summ


def LOCK_E134() -> Path:
    return H.HERE / 'experiments' / 'E134_lock.json'


def saved_stage_rows(sources: Optional[Mapping[str, pd.DataFrame]] = None, units: str = FIT_UNITS) -> pd.DataFrame:
    """The entry's step 0, item 2: saved model PN stages (glomeruli x odors of evoked rates, on the odors both have)
    against the dendrites, per source: the model's own gain per glomerulus over the median glomerulus's (its slope on the
    model's ORN input), the measured gain over it (`data_over_model`; the measured step in the rungs' reading, the
    terminals in dF/F and the dendrites in `units`), and the matrix view's error in one-fly units
    after the stage's curve (the lock's default pn_dendrite reading, no noise). `sources` maps a label to its table
    (default the analysis of 2026-11-17's versions.csv and E90's saved PN stage of the family's real wiring, APL on);
    E102's and E114's saved PN summaries (b_pn.csv, pn_stage.csv) hold no glomerulus table and are quoted as they are
    by `saved_stage_summaries`. Reads saved model outputs: run as its own part (`--part step0_saved`).

    Returns one row per source and glomerulus (model_gain, measured_gain, data_over_model on the source's odors) and,
    per source, a row with glomerulus '_all' holding the matrix view's error and the three views' (`read_prediction`
    against all PN flies on the source's odors and the shared glomeruli)."""
    from data_prep import stage_population_data as D
    glom = shared_glomeruli()
    if sources is None:
        sources = default_saved_sources()
    pop = pn_population(units)
    orn_t = recorded_stage(ORN_KEY, FIT_ORN_UNITS)
    rows = []
    for label, (pn, orn_in) in sources.items():
        od = [o for o in MEGAMAT_ODORS if o in pn.columns and o in orn_in.columns]
        gl = [g for g in glom if g in pn.index and g in orn_in.index]
        x = orn_in.loc[gl, od].to_numpy(dtype=float)
        y = pn.loc[gl, od].to_numpy(dtype=float)
        sm, _, _ = _slopes(x, y)
        meas = measured_step(fly_stack(orn_t, gl, od), fly_stack(pop.tables, gl, od), gl, od)
        mrel = _relative(sm)
        for i, g in enumerate(gl):
            rows.append(dict(source=label, glomerulus=g, n_odors=len(od), model_gain=float(mrel[i]),
                             measured_gain=float(meas.loc[g, 'gain']),
                             data_over_model=float(meas.loc[g, 'gain'] / mrel[i]) if mrel[i] > 0 else np.nan))
        tgt = Target.make(pop, od, gl)
        sc = read_prediction(pn.loc[gl, od], tgt)
        rows.append(dict(source=label, glomerulus='_all', n_odors=len(od),
                         **{k: v for k, v in sc.as_series().items()}))
    return pd.DataFrame(rows)


def default_saved_sources() -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    """The saved PN stages with glomerulus tables: each version of the analysis of 2026-11-17 (results/
    yang_orn_pn_gh146/versions.csv, in sample, 12 odors; its ORN input `apl16_orn_deltas`) and E90's family on the
    real wiring, APL on (stage_readings' `load_e90_run` and `stage_changes`, the settled pre-odor level; its own input),
    as (PN table, ORN input), odors abbreviated."""
    from kcq import stage_readings as R
    orn = H.apl16_orn_deltas()
    orn.columns = [H.odor_abbrev(c) for c in orn.columns]
    out = {}
    v = pd.read_csv(H.shared_dir() / 'results' / 'yang_orn_pn_gh146' / 'versions.csv')
    v = v[(v['sample'] == 'in_sample') & ~v.version.isin(['orn_input', 'orn_dff', 'gh146'])]
    for name, sub in v.groupby('version'):
        out[f'yang_orn_pn_gh146:{name}'] = (sub.pivot(index='glomerulus', columns='odor', values='value'), orn)
    run = R.load_e90_run('real', 'B')
    od = [o for o in R.run_odors(run) if o in MEGAMAT_ODORS]
    pn = R.stage_changes(run, 'pn_dendrite', 'on', od)
    pn.index = pn.index.astype(str)
    out['E90:B:real:on'] = (pn, R.stage_changes(run, 'orn', 'on', od))
    return out


def saved_stage_summaries() -> pd.DataFrame:
    """E102's and E114's saved PN-stage summaries against the dendrites (results/E102_antennal_lobe_normalization/
    b_pn.csv, d_pn_*.csv, e_pn_*.csv; results/E114_input_compression/pn_stage.csv), quoted with their file."""
    d = H.shared_dir() / 'results'
    out = []
    for f in sorted((d / 'E102_antennal_lobe_normalization').glob('*_pn*.csv')) + [d / 'E114_input_compression' /
                                                                                  'pn_stage.csv']:
        if f.exists():
            t = pd.read_csv(f, comment='#')
            t.insert(0, 'file', str(f.relative_to(d)))
            out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def part_step0_saved(units: str = FIT_UNITS) -> pd.DataFrame:
    """The saved PN stages against the dendrites (`saved_stage_rows`, `saved_stage_summaries`):
    step0_saved_stages.csv, step0_saved_summaries.csv, provenance_step0_saved.json. Reads saved model outputs, no
    simulation."""
    require_e134_lock()
    t0 = time.time()
    rows = saved_stage_rows(units=units)
    rows.to_csv(out_dir() / 'step0_saved_stages.csv', index=False)
    saved_stage_summaries().to_csv(out_dir() / 'step0_saved_summaries.csv', index=False)
    _write_json(out_dir() / 'provenance_step0_saved.json', dict(script=H.script_path(__file__), part='step0_saved',
                host=_host(), time=_now(), seconds=time.time() - t0, units=units, **H.code_commits()))
    return rows


def pn_population(units: str = FIT_UNITS):
    """viki's PN-dendrite population (lock-checked) in `units`: 'noise' as locked; 'dff' the same rows (and single
    trials) times their noise SD (`recorded_stage`'s reader), the criterion then not in the trials' units."""
    from data_prep import stage_population_data as D
    require_e134_lock()
    pop = D.load_population(*PN_KEY, trials=True)
    if units == 'noise':
        return pop
    if units != 'dff':
        raise ValueError(f'{units=}')
    long = D.pn_dendrite_long()
    recorded_stage(PN_KEY, 'noise')                          # the reader against the locked tables
    long = long.assign(value=long.value * long.noise_sd).drop(columns='noise_sd')
    spec = dict(units='dF/F', criterion=pop.criterion, window=pop.window, notes=pop.notes,
                recording_site=pop.recording_site)
    return D.population_from_long(PN_KEY, long, trials=True, spec=spec)


# ---------------------------------------------------------------------------
# this experiment's lock: the halves and the folds, before any fit

def odor_folds(n_folds: int, odors: Sequence[str] = MEGAMAT_ODORS, seed: int = SEED) -> List[List[str]]:
    """The odors split into `n_folds` folds of nearly equal size (a seeded permutation), each odor in one fold."""
    if not 2 <= int(n_folds) <= len(odors):
        raise ValueError(f'{n_folds=} for {len(odors)} odors')
    perm = np.random.default_rng(seed).permutation(len(odors))
    return [[odors[i] for i in sorted(f)] for f in np.array_split(perm, int(n_folds))]


def lock_rules(n_folds: int, n_perm_control: int) -> dict:
    """What the lock fixes beside the halves and the folds (ruling 3's list): the directions, the units of both sides,
    the loss and the readings, the keep rule, the permutation counts, R3's limit, the tables' scale, the clip policy
    and the imposed runs' check."""
    return dict(
        directions=[list(d) for d in DIRECTIONS], n_folds=int(n_folds),
        units=dict(fit_pn=FIT_UNITS, fit_orn=FIT_ORN_UNITS, measured_step_reading=RUNG_STEP_READING,
                   measured_step=dict(orn=STEP_READINGS[RUNG_STEP_READING][0],
                                      pn=STEP_READINGS[RUNG_STEP_READING][1]),
                   imposed=FIT_UNITS, imposed_dff_repeat='dff', step0_readings={k: list(v) for k, v in
                                                                                STEP_READINGS.items()}),
        loss=dict(fit='the matrix view\'s error in one-fly units after an affine curve set on the pooled values of the '
                      'training odors, no noise, the 28 shared glomeruli (`matrix_loss`)',
                  keep_reading=KEEP_READING, beside_reading=BESIDE_READING,
                  reading_dependent='a rung kept in one reading and not the other; the keep reading decides'),
        keep=dict(rule='lower than the highest kept rung below it (R0 counts as kept) by more than se_mult standard '
                       'errors of the paired difference over resampled flies, on the mean of the two directions, and '
                       'no higher in either direction', se_mult=KEEP_SE_MULT, n_boot=N_BOOT_KEEP),
        permutations=dict(lateral=N_PERM_LATERAL, lateral_min_exceeded=LATERAL_PERM_MIN_EXCEEDED,
                          control=int(n_perm_control), control_seed='np.random.default_rng([SEED, CONTROL_SEED_TAG, i])',
                          step2='one seeded permutation per wiring, seed SEED + the wiring\'s index in WIRINGS'),
        r3=dict(max_terms=MAX_LATERAL_TERMS, candidates=N_CANDIDATES, max_per_glomerulus=MAX_LATERAL,
                fifth='checked and reported, never added'),
        follows='slope_lo > 0 (the 90 % interval of the slope over resampled flies above 0), in each fit',
        scale=dict(rule='one positive factor per fit: the pooled SD of its evoked rates over the 28 shared glomeruli '
                        'and its training odors equals the unchanged stage\'s there, set on the megamat input and '
                        'applied to the APL16 table too', reading=SCALE_READING, named_outside=list(SCALE_FACTOR_RANGE)),
        clip=dict(policy=CLIP_POLICY, margin_hz=BOUNDS_MARGIN_HZ, name_share=CLIP_NAME_SHARE,
                  applies_to='every table entering by route 1 (imposed, its permutations, R1, the rung); the same '
                             'clipped table is route 2\'s target'),
        imposed_check=dict(max_diff_hz=IMPOSED_CHECK_HZ, pn_matrix_floor_note=PN_MATRIX_FLOOR_NOTE))


def lock_contents(n_folds: int = N_FOLDS, seed: int = SEED, n_perm_control: int = N_PERM_CONTROL,
                  budget: Optional[Mapping] = None) -> dict:
    """The halves of the PN flies and the odor folds (each odor left out once), the rules (`lock_rules`), the budget's
    choice (`choose_budget`), with the commits and E134's lock."""
    from data_prep import stage_population_data as D
    require_e134_lock()
    flies = [str(f) for f in D.load_population(*PN_KEY).flies]
    A, B = halves(flies, seed)
    return dict(experiment=EXPERIMENT, written=_now(), host=_host(), seed=seed, pn_flies=flies,
                halves=dict(half_a=A, half_b=B), folds=odor_folds(n_folds, MEGAMAT_ODORS, seed),
                odors=list(MEGAMAT_ODORS), glomeruli=shared_glomeruli(), fit_units=FIT_UNITS,
                **lock_rules(n_folds, n_perm_control), budget=dict(budget or {}),
                e134_lock_sha256=_sha256(LOCK_E134()), commits=H.code_commits())


def write_lock(n_folds: Optional[int] = None, seed: int = SEED, timing_path: Optional[Path] = None,
               n_perm_control: Optional[int] = None) -> dict:
    """Writes results/.../lock.json and its copy experiments/E135_lock.json (then commit it), once: refuses if the copy
    exists (the folds and the controls are decided before the lock and never after).

    The number of folds and of controls come from the budget's rule (`choose_budget` on the timing of one fit of each
    kind on planted data of the real size, `time_step1`: read from `timing_path`, default results/.../
    timing_step1.json, or timed now when absent): 17 folds and 5 controls unless step 1 projects above
    BUDGET_LANE_HOURS. `n_folds` / `n_perm_control` given override the rule, and the lock says so."""
    if LOCK_COPY.exists():
        raise FileExistsError(f'{LOCK_COPY} exists: the lock is written once, before any fit, and never changed')
    tp = Path(timing_path) if timing_path is not None else out_dir() / 'timing_step1.json'
    timing = _read_json(tp) if tp.exists() else part_timing()
    budget = choose_budget(timing['unit_seconds'])
    budget['timing_file'] = str(tp)
    nc, nf = budget['n_perm_control'], budget['n_folds']
    if n_folds is not None or n_perm_control is not None:
        budget['overridden'] = dict(n_folds=n_folds, n_perm_control=n_perm_control)
        nf = nf if n_folds is None else int(n_folds)
        nc = nc if n_perm_control is None else int(n_perm_control)
    c = lock_contents(nf, seed, nc, budget)
    _write_json(out_dir() / 'lock.json', c)
    _write_json(LOCK_COPY, c)
    print(f'wrote {LOCK_COPY} and {out_dir() / "lock.json"} ({nf} folds, {nc} controls, step 1 projected at '
          f'{budget["projected_lane_hours"]:.1f} lane-hours); commit {LOCK_COPY.name} before step 1', flush=True)
    return c


def read_lock(path: Path = LOCK_COPY) -> dict:
    """The committed lock (PermissionError if absent or uncommitted)."""
    if not path.exists():
        raise PermissionError(f'{path.name}: no lock yet (--part lock, then commit it); no fit before it')
    if path == LOCK_COPY and not _committed(path):
        raise PermissionError(f'{path.name} is not committed (or has uncommitted changes)')
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# step 1: the PN stage alone

@dataclass
class StageParams:
    """The olfsysm parameters of the ORN, LN and PN layers (`from_mp`); the global parameters are R1's."""
    pn_inhsc: float = 368.6631
    pn_inhadd: float = 31.4088
    pn_tanhsc: float = 5.3395
    pn_offset: float = 2.941
    pn_taum: float = 0.01
    ln_inhsc: float = 500.0
    ln_inhadd: float = 200.0
    ln_taum: float = 0.01
    ln_tauGA: float = 0.1
    ln_tauGB: float = 0.4
    ln_thr: float = 1.0
    orn_taum: float = 0.01
    dt: float = 0.0005
    pre_start: float = -2.0
    start: float = -0.5
    end: float = 0.75
    stim_start: float = 0.0
    stim_end: float = 0.5
    n_physical_gloms: int = 51
    noise_mean: float = 0.0
    noise_sd: float = 0.0

    @classmethod
    def from_mp(cls, mp) -> 'StageParams':
        return cls(pn_inhsc=mp.pn.inhsc, pn_inhadd=mp.pn.inhadd, pn_tanhsc=mp.pn.tanhsc, pn_offset=mp.pn.offset,
                   pn_taum=mp.pn.taum, ln_inhsc=mp.ln.inhsc, ln_inhadd=mp.ln.inhadd, ln_taum=mp.ln.taum,
                   ln_tauGA=mp.ln.tauGA, ln_tauGB=mp.ln.tauGB, ln_thr=mp.ln.thr, orn_taum=mp.orn.taum, dt=mp.time_dt,
                   pre_start=mp.time_pre_start, start=mp.time_start, end=mp.time_end, stim_start=mp.time_stim_start,
                   stim_end=mp.time_stim_end, n_physical_gloms=int(mp.orn.n_physical_gloms),
                   noise_mean=mp.pn.noise.mean, noise_sd=mp.pn.noise.sd)

    def with_global(self, params: Optional[Mapping[str, float]]) -> 'StageParams':
        if not params:
            return self
        bad = [k for k in params if k not in GLOBAL_PARAMS]
        if bad:
            raise KeyError(f'not global parameters: {bad}')
        return dataclasses.replace(self, **{k: float(v) for k, v in params.items()})

    def global_values(self) -> Dict[str, float]:
        return {k: float(getattr(self, k)) for k in GLOBAL_PARAMS}


def n_steps(sp: StageParams) -> int:
    """olfsysm's steps_all: unsigned((end - pre_start) / dt)."""
    return int((sp.end - sp.pre_start) / sp.dt)


def step_times(sp: StageParams) -> np.ndarray:
    """The time of every step, pre_start + k dt, rounded to 1e-9 s (pyolfsysm's `get_time_index`)."""
    return np.round(sp.pre_start + np.arange(n_steps(sp)) * sp.dt, 9)


def stim_row(sp: StageParams) -> np.ndarray:
    T = n_steps(sp)
    r = np.zeros(T)
    r[int((sp.stim_start - sp.pre_start) / sp.dt):int((sp.stim_end - sp.pre_start) / sp.dt)] = 1.0
    return r


def windows(sp: StageParams, reading: str = 'stage') -> Tuple[np.ndarray, np.ndarray]:
    """Boolean masks over the steps of a reading's window and pre-odor window.

    Args:
        sp: the layers' parameters.
        reading: 'stage', as kcq/stage_readings.py reads the model's PN dendrites (MODEL_WINDOW closed, the settled
            pre-odor window SETTLED_PRE_WINDOW closed on the left); 'preset', olfsysm's `pn_evoked_rate`, the reading of
            fit_mb_model's `pn_rates` (the steps the odor is on, [stim start, stim end), minus the steps from time_start
            to odor onset, by olfsysm's truncating step arithmetic).
    """
    T = n_steps(sp)
    if reading == 'stage':
        t = step_times(sp)
        win = (t >= MODEL_WINDOW[0] - 1e-9) & (t <= MODEL_WINDOW[1] + 1e-9)
        pre = (t >= SETTLED_PRE_WINDOW[0] - 1e-9) & (t < SETTLED_PRE_WINDOW[1] - 1e-9)
        return win, pre
    if reading == 'preset':
        start = int((sp.start - sp.pre_start) / sp.dt)
        on = int((sp.stim_start - sp.pre_start) / sp.dt)
        off = int((sp.stim_end - sp.pre_start) / sp.dt)
        k = np.arange(T)
        return (k >= on) & (k < off), (k >= start) & (k < on)
    raise ValueError(f'{reading=}')


def orn_layer(sp: StageParams, spont: np.ndarray, delta: np.ndarray) -> np.ndarray:
    """olfsysm's sim_ORN_layer for every odor at once: odors x glomeruli x steps (delta: glomeruli x odors)."""
    T = n_steps(sp)
    s = stim_row(sp)
    odor = spont[None, :, None] + delta.T[:, :, None] * s[None, None, :]
    wsize = 0.02 / sp.dt
    ext = 2.0 / (wsize + 1.0) if wsize > 1.0 else wsize
    for i in range(1, T):
        odor[:, :, i] = ext * odor[:, :, i] + (1 - ext) * odor[:, :, i - 1]
    orn = np.empty_like(odor)
    orn[:, :, 0] = spont[None, :]
    mul = sp.dt / sp.orn_taum
    for t in range(1, T):
        orn[:, :, t] = orn[:, :, t - 1] * (1.0 - mul) + odor[:, :, t] * mul
    return orn


def ln_layer(sp: StageParams, orn: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """olfsysm's sim_LN_layer for every odor at once: (inhA, inhB), odors x steps."""
    return ln_layer_mean(sp, orn.mean(axis=1), orn.shape[1])


def ln_layer_mean(sp: StageParams, m: np.ndarray, n_gloms: int) -> Tuple[np.ndarray, np.ndarray]:
    """olfsysm's sim_LN_layer from the mean ORN rate over the `n_gloms` glomeruli (odors x steps), which is all it
    reads: (inhA, inhB), odors x steps."""
    O, T = m.shape
    scaling = float(n_gloms) / float(sp.n_physical_gloms)
    pot = np.full(O, 300.0)
    resp = np.ones(O)
    inhA = np.full((O, T), 50.0)
    inhB = np.full((O, T), 50.0)
    inh_ln = np.zeros(O)
    for t in range(1, T):
        dA = -inhA[:, t - 1] + resp
        dB = -inhB[:, t - 1] + resp
        dLN = -pot + np.power(m[:, t - 1] * scaling, 3.0) / scaling / 2.0 * inh_ln
        inhA[:, t] = inhA[:, t - 1] + dA * sp.dt / sp.ln_tauGA
        inhB[:, t] = inhB[:, t - 1] + dB * sp.dt / sp.ln_tauGB
        inh_ln = sp.ln_inhsc / (sp.ln_inhadd + inhA[:, t])
        pot = pot + dLN * sp.dt / sp.ln_taum
        resp = (pot - sp.ln_thr) * (pot > sp.ln_thr)
    return inhA, inhB


def pn_layer(sp: StageParams, spont: np.ndarray, orn: np.ndarray, inhA: np.ndarray, inhB: np.ndarray,
             M: Optional[np.ndarray] = None, home: Optional[np.ndarray] = None, full: bool = False,
             reading: str = 'stage'):
    """olfsysm's sim_PN_layer for every odor at once, on virtual PNs: PN n reads the drive orn_delta @ M[:, n] (orn_delta
    the ORN rates minus their spontaneous rates, odors x glomeruli) in place of its own glomerulus's ORN delta, with the
    spontaneous rates of its home glomerulus `home[n]`. M the identity (default) is olfsysm's PN layer.

    Returns the evoked rates (virtual PNs x odors) in `reading` (`windows`), or with `full` the rates over time (odors x
    PNs x steps)."""
    if sp.noise_sd != 0 or sp.noise_mean != 0:
        raise ValueError('the port has no PN noise (olfsysm adds it per step); the model runs with none')
    O, G, T = orn.shape
    if M is None:
        M = np.eye(G)
        home = np.arange(G)
    home = np.asarray(home, dtype=int)
    N = M.shape[1]
    spont_h = spont[home]
    spont_pn = spont_h * sp.pn_inhsc / (spont.sum() + sp.pn_inhadd)
    inh_pn_t = sp.pn_inhsc / (sp.pn_inhadd + 0.25 * inhA + 0.75 * inhB)     # odors x steps
    win, pre = windows(sp, reading)
    pn = np.broadcast_to(spont_h[None, :], (O, N)).copy()
    acc_w = pn * win[0]
    acc_p = pn * pre[0]
    out = np.empty((O, N, T)) if full else None
    if full:
        out[:, :, 0] = pn
    inh = np.zeros((O, 1))
    identity = M is not None and N == G and np.array_equal(M, np.eye(G)) and np.array_equal(home, np.arange(G))
    for t in range(1, T):
        od = orn[:, :, t - 1] - spont[None, :]
        drive = od if identity else od @ M
        d = -pn + spont_pn[None, :] + 200.0 * np.tanh((drive + sp.pn_offset) * sp.pn_tanhsc / 200.0 * inh)
        inh = inh_pn_t[:, t:t + 1]
        pn = pn + d * sp.dt / sp.pn_taum
        pn = np.where(0.0 < pn, pn, 0.0)
        if win[t]:
            acc_w += pn
        if pre[t]:
            acc_p += pn
        if full:
            out[:, :, t] = pn
    if full:
        return out
    return (acc_w / win.sum() - acc_p / pre.sum()).T


@dataclass
class Terms:
    """The glomerulus terms of a rung: gains per shared glomerulus (1 where absent) and lateral terms
    (glomerulus, partner, weight)."""
    gains: Dict[str, float] = field(default_factory=dict)
    laterals: List[Tuple[str, str, float]] = field(default_factory=list)

    def matrix(self, glomeruli: Sequence[str]) -> np.ndarray:
        """The drive matrix M (glomeruli x glomeruli): PN g reads sum_h orn_delta_h M[h, g]."""
        idx = {g: i for i, g in enumerate(glomeruli)}
        M = np.eye(len(glomeruli))
        for g, v in self.gains.items():
            M[idx[g], idx[g]] = float(v)
        n_lat: Dict[str, int] = {}
        for g, p, w in self.laterals:
            if g == p:
                raise ValueError(f'{g}: a lateral term from itself')
            n_lat[g] = n_lat.get(g, 0) + 1
            if n_lat[g] > MAX_LATERAL:
                raise ValueError(f'{g}: more than {MAX_LATERAL} lateral terms')
            M[idx[p], idx[g]] += float(w)
        return M

    def as_dict(self) -> dict:
        return dict(gains={k: float(v) for k, v in self.gains.items()},
                    laterals=[[g, p, float(w)] for g, p, w in self.laterals])

    @classmethod
    def from_dict(cls, d: Mapping) -> 'Terms':
        return cls(gains=dict(d.get('gains', {})), laterals=[tuple(x) for x in d.get('laterals', [])])


class PNStage:
    """olfsysm's ORN, LN and PN layers ported to numpy (`orn_layer`, `ln_layer`, `pn_layer`), on one ORN input, with
    the glomerulus terms on each PN's drive: the LN layer reads the recorded input (the mean ORN rate over glomeruli),
    so the normalization keeps it. Identical to olfsysm's layers with no terms (`check_identity`).

    Args:
        delta: glomeruli x odors, the ORN input's deltas in the model's glomerulus order (all of mp.orn.data.delta's
            rows, those without ORN data at 0).
        spont: the glomeruli's spontaneous ORN rates (mp.orn.data.spont).
        sp: the layers' parameters.
    """

    def __init__(self, delta: pd.DataFrame, spont: pd.Series, sp: Optional[StageParams] = None):
        self.delta = delta.astype(float)
        self.glomeruli = [str(g) for g in delta.index]
        self.odors = [str(o) for o in delta.columns]
        self.spont = spont.reindex(delta.index).to_numpy(dtype=float)
        if np.isnan(self.spont).any():
            raise ValueError('spontaneous rates missing for some glomeruli')
        self.sp = StageParams() if sp is None else sp
        self._orn: Dict[str, np.ndarray] = {}
        self._ln: Dict[tuple, Tuple[np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_mp(cls, mp, glomeruli: Sequence[str], odors: Sequence[str]) -> 'PNStage':
        """From olfsysm's ModelParams (its delta's rows are `glomeruli`, its columns `odors` in olfsysm's order)."""
        d = np.asarray(mp.orn.data.delta, dtype=float)
        s = np.asarray(mp.orn.data.spont, dtype=float).ravel()
        return cls(pd.DataFrame(d, index=list(glomeruli), columns=list(odors)), pd.Series(s, index=list(glomeruli)),
                   StageParams.from_mp(mp))

    def _key(self, delta: np.ndarray) -> str:
        return hashlib.sha1(np.ascontiguousarray(delta).tobytes()).hexdigest()

    def orn(self, delta: Optional[np.ndarray] = None) -> np.ndarray:
        d = self.delta.to_numpy() if delta is None else np.asarray(delta, dtype=float)
        k = self._key(d)
        if k not in self._orn:
            if len(self._orn) > 4:
                self._orn.clear()
            self._orn[k] = orn_layer(self.sp, self.spont, d)
        return self._orn[k]

    def ln(self, sp: StageParams, delta: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        d = self.delta.to_numpy() if delta is None else np.asarray(delta, dtype=float)
        k = (self._key(d), sp.ln_inhsc, sp.ln_inhadd, sp.ln_taum, sp.ln_tauGA, sp.ln_tauGB, sp.ln_thr,
             sp.n_physical_gloms)
        if k not in self._ln:
            if len(self._ln) > 64:
                self._ln.clear()
            self._ln[k] = ln_layer(sp, self.orn(d))
        return self._ln[k]

    def run(self, params: Optional[Mapping[str, float]] = None, terms: Optional[Terms] = None,
            delta: Optional[pd.DataFrame] = None, reading: str = 'stage') -> pd.DataFrame:
        """The PN's odor-evoked rates in `reading` (`windows`), glomeruli x odors (Hz), at the global parameters
        `params` (default the stage's) with the glomerulus terms `terms`, on the stage's input (or `delta`, which the
        LN then reads too)."""
        sp = self.sp.with_global(params)
        d = None if delta is None else delta.reindex(index=self.glomeruli, columns=self.odors).to_numpy(dtype=float)
        orn = self.orn(d)
        inhA, inhB = self.ln(sp, d)
        M = None if terms is None else terms.matrix(self.glomeruli)
        home = None if terms is None else np.arange(len(self.glomeruli))
        x = pn_layer(sp, self.spont, orn, inhA, inhB, M, home, reading=reading)
        return pd.DataFrame(x, index=self.glomeruli, columns=self.odors)

    def run_virtual(self, M: np.ndarray, home: np.ndarray, params: Optional[Mapping[str, float]] = None,
                    reading: str = 'stage') -> np.ndarray:
        """Virtual PNs (`pn_layer`) on the stage's input: evoked rates, virtual PNs x odors."""
        sp = self.sp.with_global(params)
        inhA, inhB = self.ln(sp)
        return pn_layer(sp, self.spont, self.orn(), inhA, inhB, M, home, reading=reading)

    def rate_bounds(self, params: Optional[Mapping[str, float]] = None, reading: str = 'preset',
                    delta: Optional[pd.DataFrame] = None, shift: float = 1e12) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """The evoked rates a PN can reach at the stage's LN (lo: its glomerulus's ORN delta far down, the PN silenced
        from odor onset; hi: far up, its tanh saturated), as olfsysm's `pn_rate_bounds` reads them with the LN kept:
        (lo, hi), glomeruli x odors."""
        sp = self.sp.with_global(params)
        d = None if delta is None else delta.reindex(index=self.glomeruli, columns=self.odors).to_numpy(dtype=float)
        inhA, inhB = self.ln(sp, d)
        k = orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
        rows, cols = np.arange(len(self.glomeruli)), np.arange(len(self.odors))
        one = np.ones((len(rows), len(cols)))
        lo = _rates_at(sp, self.spont, rows, cols, -shift * one, k, inhA, inhB, reading)
        hi = _rates_at(sp, self.spont, rows, cols, shift * one, k, inhA, inhB, reading)
        mk = lambda a: pd.DataFrame(a, index=self.glomeruli, columns=self.odors)
        return mk(lo), mk(hi)

    def check_identity(self, mp, rv, tol: float = 1e-9, params: Optional[Mapping[str, float]] = None) -> dict:
        """The port against olfsysm's layers (`run_ORN_LN_sims`, `run_PN_sims` on `mp`, `rv`, whose input is this
        stage's): the ORN rates, the LN's inhA and inhB and the PN rates over time, at the stage's parameters or
        `params`. Raises RuntimeError beyond `tol`."""
        import olfsysm as osm
        sp = self.sp.with_global(params)
        keep = {k: getattr(mp.pn, k) for k in ('inhsc', 'inhadd', 'tanhsc')}
        keep_ln = {k: getattr(mp.ln, k) for k in ('inhsc', 'inhadd')}
        try:
            mp.pn.inhsc, mp.pn.inhadd, mp.pn.tanhsc = sp.pn_inhsc, sp.pn_inhadd, sp.pn_tanhsc
            mp.ln.inhsc, mp.ln.inhadd = sp.ln_inhsc, sp.ln_inhadd
            osm.run_ORN_LN_sims(mp, rv)
            osm.run_PN_sims(mp, rv)
            o_ref = np.stack([np.asarray(x, dtype=float) for x in rv.orn.sims])
            a_ref = np.stack([np.asarray(x, dtype=float).ravel() for x in rv.ln.ln_inhA.sims])
            b_ref = np.stack([np.asarray(x, dtype=float).ravel() for x in rv.ln.ln_inhB.sims])
            p_ref = np.stack([np.asarray(x, dtype=float) for x in rv.pn.pn_sims])
        finally:
            for k, v in keep.items():
                setattr(mp.pn, k, v)
            for k, v in keep_ln.items():
                setattr(mp.ln, k, v)
        orn = self.orn()
        inhA, inhB = self.ln(sp)
        pn = pn_layer(sp, self.spont, orn, inhA, inhB, full=True)
        if o_ref.shape[-1] != orn.shape[-1]:
            raise RuntimeError(f'olfsysm has {o_ref.shape[-1]} steps, the port {orn.shape[-1]} (pretime removed?)')
        diffs = dict(orn=float(np.abs(o_ref - orn).max()), inhA=float(np.abs(a_ref - inhA).max()),
                     inhB=float(np.abs(b_ref - inhB).max()), pn=float(np.abs(p_ref - pn).max()))
        if hasattr(osm, 'pn_evoked_rates'):
            mp.pn.inhsc, mp.pn.inhadd, mp.pn.tanhsc = sp.pn_inhsc, sp.pn_inhadd, sp.pn_tanhsc
            try:
                e_ref = np.asarray(osm.pn_evoked_rates(mp, rv), dtype=float)
            finally:
                for k, v in keep.items():
                    setattr(mp.pn, k, v)
            e = pn_layer(sp, self.spont, orn, inhA, inhB, reading='preset')
            diffs['evoked_preset'] = float(np.abs(e_ref - e).max())
        res = dict(max_abs_diff=diffs, tol=tol, ok=all(v <= tol for v in diffs.values()), params=sp.global_values())
        if not res['ok']:
            raise RuntimeError(f'the port differs from olfsysm: {res}')
        return res


def flies_subset(pop, flies: Sequence[str]):
    """viki's D.Population restricted to `flies` (tables and trials)."""
    flies = [f for f in pop.flies if f in set(flies)]
    return dataclasses.replace(pop, flies=tuple(flies), tables={f: pop.tables[f] for f in flies},
                               trials=None if pop.trials is None else {f: pop.trials[f] for f in flies
                                                                       if f in pop.trials})


def fit_reading():
    """The reading inside a fit (ruling 4): the lock's pn_dendrite reading with an affine curve (smooth in the
    parameters), no noise, on the 28 shared glomeruli; `matrix_loss` computes it in closed form, the curve set on the
    fit's training odors."""
    from kcq import stage_readings as R
    return R.Reading('fit_affine_nonoise', 'pn_dendrite', curve='affine', noise='none', glomeruli='shared28')


def read_reading(name: str = KEEP_READING):
    """The reading of a held-out prediction: the lock's pn_dendrite reading called `name` (ruling 4: KEEP_READING, the
    default reading, the quantile curve and the flies' noise, decides; BESIDE_READING, its 'affine_noise' alternative,
    is read beside)."""
    from kcq import stage_readings as R
    if name == R.DEFAULT_READINGS['pn_dendrite'].name:
        return R.DEFAULT_READINGS['pn_dendrite']
    return R.reading_by_name('pn_dendrite', name)


@dataclass
class Target:
    """The dendrites a fit or a reading is held to: viki's population of the flies of one half, the odors, the shared
    glomeruli, and hal's reference of those flies (in their own one-fly units unless `units_from` is given)."""
    flies: object
    odors: List[str]
    glomeruli: List[str]
    ref: P.FlyReference
    _loss: Optional[tuple] = field(default=None, repr=False)

    def loss_terms(self) -> tuple:
        """(the flies' pooled mean and variance, the mean matrix part in glomeruli x odors order, its one-fly unit),
        computed once, for `matrix_loss`."""
        if self._loss is None:
            from kcq import stage_readings as R
            fv, fw = R.pooled_fly_values(self.flies, self.odors)
            w = fw / fw.sum()
            mf = float(np.sum(w * fv))
            vf = float(np.sum(w * (fv - mf) ** 2))
            m = self.ref.mean['matrix/values'].reindex([f'{g}{P.PAIR_SEP}{o}' for g in self.glomeruli
                                                        for o in self.odors]).to_numpy(dtype=float)
            self._loss = (mf, vf, m, float(self.ref.one_fly['matrix/values']))
        return self._loss

    @classmethod
    def make(cls, flies, odors: Sequence[str], glomeruli: Sequence[str], parts: Optional[Sequence[str]] = None,
             units_from: Optional[P.FlyReference] = None) -> 'Target':
        from data_prep import stage_population_data as D
        fl = P.fly_populations(flies)
        fl = [dataclasses.replace(f, values=f.values.reindex(index=list(glomeruli)),
                                  responds=f.responds.reindex(index=list(glomeruli))) for f in fl]
        ref = P.fly_reference(fl, list(odors), D.breadth_bins(len(odors)), parts)
        if units_from is not None:
            ref = dataclasses.replace(ref, one_fly={p: units_from.one_fly[p] for p in ref.one_fly})
        return cls(flies, list(odors), list(glomeruli), ref)


def read_prediction(pred: pd.DataFrame, target: Target, reading=None, seed_index: int = 0) -> P.Score:
    """A predicted table of PN rates (glomeruli x odors) read as the dendrites are recorded (the reading's curve set on
    the prediction's pooled values against the target flies' pooled values, then its noise and criterion) and scored
    against the target in one-fly units (the matrix view and the three views)."""
    from kcq import stage_readings as R
    reading = read_reading() if reading is None else reading
    U = pred.reindex(index=target.glomeruli, columns=target.odors).astype(float)
    curve = R.set_curve(R.pooled_values([U]), target.flies, target.odors, reading, seed_index)
    pop, _, _ = R.finish_population(U, 'pn_dendrite', 'on', reading, curve, target.flies, target.odors, seed_index,
                                    'model')
    return P.score(pop, target.ref)


def matrix_loss(pred: pd.DataFrame, target: Target) -> float:
    """The fit's loss: the matrix view's error in one-fly units after an affine curve set on the pooled values (no
    noise; `fit_reading`; stage_readings' 'affine' curve: the model's pooled mean and SD onto the flies')."""
    return matrix_loss_array(pred.reindex(index=target.glomeruli, columns=target.odors).to_numpy(dtype=float), target)


def matrix_loss_array(U: np.ndarray, target: Target) -> float:
    """`matrix_loss` of an array in the target's glomeruli x odors order."""
    mf, vf, m, unit = target.loss_terms()
    x = U[np.isfinite(U)]
    sm = float(x.std())
    scale = float(np.sqrt(vf) / sm) if sm > 0 else 0.0
    d = (scale * U + (mf - scale * float(x.mean()))).ravel() - m
    ok = np.isfinite(d)
    return float(np.sqrt(np.mean(d[ok] ** 2)) / unit)


def _log_bounds(name: str) -> Tuple[float, float]:
    lo, hi = GLOBAL_BOUNDS[name]
    d = OLFSYSM_DEFAULTS[name]
    return float(np.log(lo * d)), float(np.log(hi * d))


def fit_global(stage: PNStage, target: Target, start: Optional[Mapping[str, float]] = None,
               names: Sequence[str] = GLOBAL_PARAMS, maxiter: int = 400) -> dict:
    """R1: the global parameters fitted to the target (`matrix_loss`), Nelder-Mead on their logs inside GLOBAL_BOUNDS
    (a point outside scores as 10), from `start` (default the stage's)."""
    x0 = stage.sp.global_values() if start is None else dict(start)
    lo = np.array([_log_bounds(n)[0] for n in names])
    hi = np.array([_log_bounds(n)[1] for n in names])

    def f(z):
        if np.any(z < lo) or np.any(z > hi):
            return 10.0
        p = dict(zip(names, np.exp(z)))
        return matrix_loss(stage.run(p).reindex(index=target.glomeruli), target)
    z0 = np.log([x0[n] for n in names])
    r = optimize.minimize(f, z0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-6, maxiter=maxiter))
    p = dict(zip(names, map(float, np.exp(r.x))))
    at = [n for n, z, a, b in zip(names, r.x, lo, hi) if min(z - a, b - z) < 0.02]
    return dict(params=p, loss=float(r.fun), nfev=int(r.nfev), at_bounds=at)


def gain_grid(stage: PNStage, glomeruli: Sequence[str], params: Mapping[str, float], base: Terms,
              grid: np.ndarray) -> np.ndarray:
    """Each shared glomerulus's PN rates with its gain at every grid value and its laterals as in `base` (virtual PNs,
    one batched run): glomeruli x grid x odors. PNs are independent given the LN, so one run gives every glomerulus."""
    G = stage.glomeruli
    idx = {g: i for i, g in enumerate(G)}
    M0 = base.matrix(G)
    cols, home = [], []
    for g in glomeruli:
        for v in grid:
            c = M0[:, idx[g]].copy()
            c[idx[g]] = v
            cols.append(c)
            home.append(idx[g])
    X = stage.run_virtual(np.column_stack(cols), np.array(home), params)
    return X.reshape(len(glomeruli), len(grid), -1)


def fit_gains(stage: PNStage, target: Target, params: Mapping[str, float], start: Mapping[str, float],
              base: Optional[Terms] = None, n_grid: int = 41, maxiter: int = 200) -> dict:
    """R2's fit: the gains of the shared glomeruli (log, inside GAIN_BOUNDS) fitted to the target by L-BFGS-B, with the
    PN rates interpolated from one batched run over a log-spaced grid of gains (`gain_grid`) and re-run exactly at the
    end; from `start` (the measured gains)."""
    base = Terms() if base is None else base
    gl = target.glomeruli
    grid = np.exp(np.linspace(np.log(GAIN_BOUNDS[0]), np.log(GAIN_BOUNDS[1]), n_grid))
    oi = [stage.odors.index(o) for o in target.odors]
    Y = gain_grid(stage, gl, params, base, grid)[:, :, oi]
    lg = np.log(grid)
    step = lg[1] - lg[0]
    ar = np.arange(len(gl))

    def f(z):
        pos = np.clip((z - lg[0]) / step, 0.0, len(lg) - 1.0)
        i0 = np.minimum(np.floor(pos).astype(int), len(lg) - 2)
        w = (pos - i0)[:, None]
        return matrix_loss_array(Y[ar, i0] * (1 - w) + Y[ar, i0 + 1] * w, target)
    z0 = np.log(np.clip([start.get(g, 1.0) for g in gl], *GAIN_BOUNDS))
    r = optimize.minimize(f, z0, method='L-BFGS-B', bounds=[(lg[0], lg[-1])] * len(gl),
                          options=dict(maxiter=maxiter, eps=1e-4))
    gains = {g: float(np.exp(v)) for g, v in zip(gl, r.x)}
    terms = Terms(gains=gains, laterals=list(base.laterals))
    exact = matrix_loss(stage.run(params, terms), target)
    at = [g for g, v in zip(gl, r.x) if min(v - lg[0], lg[-1] - v) < 1e-3]
    return dict(terms=terms, loss=exact, loss_interp=float(r.fun), nit=int(r.nit), at_bounds=at)


def fit_lateral(stage: PNStage, target: Target, params: Mapping[str, float], terms: Terms, g: str, partner: str,
                n_grid: int = 41) -> dict:
    """One lateral term's weight (inside LATERAL_BOUNDS) with everything else held, by a grid over the weight and a
    bounded scalar search on the interpolated rates, re-run exactly at the end."""
    G = stage.glomeruli
    idx = {x: i for i, x in enumerate(G)}
    M0 = terms.matrix(G)
    grid = np.linspace(LATERAL_BOUNDS[0], LATERAL_BOUNDS[1], n_grid)
    cols = []
    for w in grid:
        c = M0[:, idx[g]].copy()
        c[idx[partner]] += w
        cols.append(c)
    Y = stage.run_virtual(np.column_stack(cols), np.full(n_grid, idx[g]), params)
    rest = stage.run(params, terms)

    def f(w):
        X = rest.copy()
        X.loc[g] = [np.interp(w, grid, Y[:, j]) for j in range(Y.shape[1])]
        return matrix_loss(X, target)
    r = optimize.minimize_scalar(f, bounds=LATERAL_BOUNDS, method='bounded', options=dict(xatol=1e-3))
    t = Terms(gains=dict(terms.gains), laterals=list(terms.laterals) + [(g, partner, float(r.x))])
    return dict(terms=t, weight=float(r.x), loss=matrix_loss(stage.run(params, t), target))


def model_gains(stage: PNStage, params: Mapping[str, float], glomeruli: Sequence[str], odors: Sequence[str]) -> pd.Series:
    """The model stage's own gain per glomerulus: the slope through the origin of its evoked PN rate on its ORN input
    over `odors`, over the median glomerulus's."""
    x = stage.delta.loc[list(glomeruli), list(odors)].to_numpy(dtype=float)
    y = stage.run(params).loc[list(glomeruli), list(odors)].to_numpy(dtype=float)
    s, _, _ = _slopes(x, y)
    return pd.Series(_relative(s), index=list(glomeruli))


def measured_gain_start(measured: pd.DataFrame, model_rel: pd.Series) -> Dict[str, float]:
    """R2's starting gains: the measured gain over the model stage's own gain (both over their median glomerulus),
    1 where the dendrites do not follow their terminals (`follows_rule`, the 'follows' column) or either is undefined,
    clipped to GAIN_BOUNDS."""
    out = {}
    for g in model_rel.index:
        ok = g in measured.index and bool(measured.loc[g, 'follows'])
        v = float(measured.loc[g, 'gain'] / model_rel[g]) if ok else 1.0
        out[g] = float(np.clip(v if np.isfinite(v) and v > 0 else 1.0, *GAIN_BOUNDS))
    return out


def keep_rung(err_rung: np.ndarray, err_below: np.ndarray, point_rung: float, point_below: float,
              mult: float = KEEP_SE_MULT, dir_rung: Optional[Sequence[float]] = None,
              dir_below: Optional[Sequence[float]] = None) -> Tuple[bool, float]:
    """Kept (ruling 6) if the rung's held-out error (the mean of the two directions) is lower than the comparator's (the
    highest kept rung below it) by more than `mult` standard errors of the paired difference over resamplings of the
    flies, and, with the per-direction errors `dir_rung` / `dir_below` given, no higher in either direction. Returns
    (kept, the standard error)."""
    d = np.asarray(err_below, dtype=float) - np.asarray(err_rung, dtype=float)
    se = float(np.nanstd(d, ddof=1))
    kept = (point_below - point_rung) > mult * se
    if dir_rung is not None and dir_below is not None:
        kept = kept and all(float(a) <= float(b) for a, b in zip(dir_rung, dir_below))
    return bool(kept), se


def _directions_of(o: Mapping) -> List[str]:
    return [k for k in o if '->' in k]


def keep_against(o: Mapping, below: Mapping, reading: str = KEEP_READING, mult: float = KEEP_SE_MULT) -> dict:
    """`keep_rung` of two cross-fits (`cross_fit`'s outputs) in one reading: kept, the SE, the margin, and whether the
    rung is no higher in either direction."""
    a, b = o['by_reading'][reading], below['by_reading'][reading]
    ds = _directions_of(o)
    da, db = [a['directions'][d] for d in ds], [b['directions'][d] for d in ds]
    kept, se = keep_rung(a['draws'], b['draws'], a['point'], b['point'], mult, da, db)
    return dict(kept=kept, se=se, diff=float(b['point'] - a['point']),
                no_higher_each_direction=bool(all(x <= y for x, y in zip(da, db))))


def permute_odor_labels(flies, rng: np.random.Generator, odors: Optional[Sequence[str]] = None):
    """viki's population with its odor labels permuted among `odors` (default all the panel's; one permutation of the
    columns, for every fly and its trials), the other columns kept."""
    odors = list(flies.odors if odors is None else odors)
    perm = [odors[i] for i in rng.permutation(len(odors))]
    ren = dict(zip(perm, odors))
    tables = {}
    for f, t in flies.tables.items():
        x = t.copy()
        x[odors] = t[perm].to_numpy()
        tables[f] = x
    trials = None
    if flies.trials:
        trials = {}
        for f, t in flies.trials.items():
            x = t.copy()
            x.columns = pd.MultiIndex.from_tuples([(ren.get(o, o), r) for o, r in t.columns], names=t.columns.names)
            trials[f] = x.sort_index(axis=1)
    return dataclasses.replace(flies, tables=tables, trials=trials)


def permute_glomerulus_labels(x: pd.DataFrame, glomeruli: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    """A glomeruli x odors table with its shared glomeruli's rows permuted among themselves."""
    g = list(glomeruli)
    perm = [g[i] for i in rng.permutation(len(g))]
    out = x.copy()
    out.loc[g] = x.loc[perm].to_numpy()
    return out


def permute_terms(terms: Terms, glomeruli: Sequence[str], rng: np.random.Generator) -> Terms:
    """A rung's terms with their glomerulus labels permuted among the shared glomeruli (gains and laterals' targets
    move together; a lateral's partner keeps its label)."""
    g = list(glomeruli)
    m = dict(zip(g, [g[i] for i in rng.permutation(len(g))]))
    lat = [(m[a], b, w) for a, b, w in terms.laterals if m[a] != b]
    return Terms(gains={m[k]: v for k, v in terms.gains.items()}, laterals=lat)


def held_out_draws(pred: pd.DataFrame, target: Target, n_boot: int = N_BOOT_KEEP, seed: int = SEED,
                   reading=None) -> np.ndarray:
    """The matrix view's error of a held-out prediction against resampled sets of the target's flies (with replacement)
    in the measured set's one-fly units, one per draw (the same draws for every prediction at the same seed)."""
    from kcq import stage_readings as R
    reading = read_reading() if reading is None else reading
    U = pred.reindex(index=target.glomeruli, columns=target.odors).astype(float)
    curve = R.set_curve(R.pooled_values([U]), target.flies, target.odors, reading, 0)
    pop, _, _ = R.finish_population(U, 'pn_dendrite', 'on', reading, curve, target.flies, target.odors, 0, 'model')
    v = P.population_views(pop, target.odors, target.ref.breadth_edges, ['matrix/values'])['matrix/values']
    A = target.ref.parts['matrix/values'].to_numpy(dtype=float)
    m = v.reindex(target.ref.parts['matrix/values'].columns).to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    out = np.empty(n_boot)
    for b in range(n_boot):
        mb = _nanmean(A[rng.integers(0, n, n)])
        d = m - mb
        out[b] = np.sqrt(np.nanmean(d ** 2)) / target.ref.one_fly['matrix/values']
    return out


def fold_map(table: pd.DataFrame, target: Target) -> pd.DataFrame:
    """A fold's table of PN rates in the fit flies' units: through the fit's own affine curve (`matrix_loss`'s), set on
    the fold's training odors only, applied to every odor. A curve absorbs the rates' scale, so the global parameters
    can change it from fold to fold; the folds' held-out columns are comparable only through these maps, as in the
    analysis of 2026-11-17 (its `loo`)."""
    mf, vf, _, _ = target.loss_terms()
    x = table.reindex(index=target.glomeruli, columns=target.odors).to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    sm = float(x.std())
    scale = float(np.sqrt(vf) / sm) if sm > 0 else 0.0
    return table * scale + (mf - scale * float(x.mean()))


def cross_fit(fit_one: Callable[[Target, List[str]], Tuple[pd.DataFrame, dict]], folds: Sequence[Sequence[str]],
              fit_targets: Mapping[str, Callable[[List[str]], Target]], read_targets: Mapping[str, Target],
              directions: Sequence[Tuple[str, str]] = DIRECTIONS, n_boot: int = N_BOOT_KEEP, seed: int = SEED,
              map_folds: bool = True, readings: Sequence[str] = (KEEP_READING, BESIDE_READING)) -> dict:
    """A rung fitted fold by fold on one half and read on the other: per direction (fit half, read half), for each
    fold the fit on the fit half's flies and the other folds' odors (`fit_one(target, train_odors)` -> (the table of PN
    rates on every odor, the fit's record)), the held-out columns assembled into one table (each odor by the fit that
    left it out; with `map_folds` each fold's table first in the fit flies' units, `fold_map`), read on the read half in
    each of `readings` (`read_prediction`) with bootstrap draws of its matrix error (`held_out_draws`).

    Returns:
        Per direction ('<fit>-><read>') the table (and 'table_hz', the folds' own rates assembled unmapped), the folds'
        records, 'by_reading' ({reading: score, draws, point}) and, of the first reading (the keep reading), 'score',
        'draws' and 'point'; at the top 'point' and 'draws' (the mean over directions, first reading) and 'by_reading'
        ({reading: point, draws, directions: {direction: point}}).
    """
    out = {}
    for fit_half, read_half in directions:
        cols, cols_hz, recs = {}, {}, []
        for k, fold in enumerate(folds):
            train = [o for f in folds for o in f if o not in set(fold)]
            t = fit_targets[fit_half](train)
            tab, rec = fit_one(t, train)
            mapped = fold_map(tab, t) if map_folds else tab
            for o in fold:
                cols[o] = mapped[o]
                cols_hz[o] = tab[o]
            recs.append(dict(fold=k, held_out=list(fold), **rec))
        tgt = read_targets[read_half]
        table = pd.DataFrame(cols)[tgt.odors]
        by = {}
        for name in readings:
            rd = read_reading(name)
            sc = read_prediction(table, tgt, rd)
            by[name] = dict(score=sc, draws=held_out_draws(table, tgt, n_boot, seed, rd),
                            point=float(sc.parts['matrix/values']))
        first = by[readings[0]]
        out[f'{fit_half}->{read_half}'] = dict(table=table, table_hz=pd.DataFrame(cols_hz)[tgt.odors], folds=recs,
                                                by_reading=by, score=first['score'], draws=first['draws'],
                                                point=first['point'])
    keys = list(out)
    out['by_reading'] = {name: dict(point=float(np.mean([out[k]['by_reading'][name]['point'] for k in keys])),
                                    draws=np.mean([out[k]['by_reading'][name]['draws'] for k in keys], axis=0),
                                    directions={k: out[k]['by_reading'][name]['point'] for k in keys})
                         for name in readings}
    out['point'] = out['by_reading'][readings[0]]['point']
    out['draws'] = out['by_reading'][readings[0]]['draws']
    return out


N_BOOT_MEASURED = 200


def measured_for(target: Target, orn: np.ndarray, orn_partners: np.ndarray, partner_names: Sequence[str],
                 n_boot: int = N_BOOT_MEASURED, seed: int = SEED, odors_all: Sequence[str] = MEGAMAT_ODORS) -> dict:
    """The measured step on a fit target's flies and odors (ruling 1's reading: the ORN terminals in dF/F, the
    dendrites in the fit's units, the target's): its gains with the follows rule (`follows_rule`, ruling 2, on the
    flies and odors this fit sees) and its candidate partners. `orn` and `orn_partners` are the ORN recordings in dF/F
    (flies x the target's glomeruli / the partner pool x `odors_all`)."""
    oi = [list(odors_all).index(o) for o in target.odors]
    pn = fly_stack(target.flies.tables, target.glomeruli, target.odors)
    o = orn[:, :, oi]
    st = measured_step(o, pn, target.glomeruli, target.odors)
    iv, _ = step_intervals(o, pn, target.glomeruli, target.odors, n_boot, seed)
    st = st.join(iv)
    st['follows'] = follows_rule(st)
    res = residuals(o, pn, st, target.glomeruli, target.odors)
    om = pd.DataFrame(_nanmean(orn_partners[:, :, oi]), index=list(partner_names), columns=target.odors)
    cp = candidate_partners(res, om, n=N_CANDIDATES, n_perm=0)
    return dict(step=st, partners=cp)


def _key(t: Target) -> tuple:
    return (tuple(t.flies.flies), tuple(t.odors))


def control_seed(i: int, seed: int = SEED) -> int:
    """The seed of R1's and R2's odor-permuted refit number `i` (a function of the index only, so that each control
    runs as a process of its own)."""
    return int(np.random.default_rng([int(seed), CONTROL_SEED_TAG, int(i)]).integers(1 << 30))


def lateral_min_exceeded(n_perm: int) -> int:
    """How many of `n_perm` permuted fits a lateral term must out-gain: LATERAL_PERM_MIN_EXCEEDED of N_PERM_LATERAL
    (19 of 20), the same share for other counts (rounded up)."""
    return int(np.ceil(LATERAL_PERM_MIN_EXCEEDED / N_PERM_LATERAL * int(n_perm) - 1e-9)) if n_perm else 0


class Step1:
    """Step 1's shared pieces: the PN stage, the PN population in the fit's units (its odors named as the stage's), the
    ORN recordings in dF/F (the measured step's ORN side, ruling 1), the lock's glomeruli, folds, halves and directions,
    the targets of each half, and the cross-fit on them (`cf`) with the fold fitters of R1, R2m and R2f (`fitters`).

    Args:
        stage: the PN stage (`family_stage`, or planted).
        pop: viki's PN population in the fit's units.
        orn, orn_partners, partner_names: the ORN recordings (flies x the lock's glomeruli x the megamat odors, flies x
            the partner pool x the odors) and the pool's names.
        lock: gives 'glomeruli', 'folds', 'halves' and (optionally) 'directions'.
    """

    def __init__(self, stage: PNStage, pop, orn: np.ndarray, orn_partners: np.ndarray, partner_names: Sequence[str],
                 lock: Mapping, n_boot: int = N_BOOT_KEEP, seed: int = SEED, global_maxiter: int = 400,
                 readings: Sequence[str] = (KEEP_READING, BESIDE_READING), verbose: bool = True):
        self.stage, self.pop, self.orn, self.orn_partners = stage, pop, orn, orn_partners
        self.partner_names = list(partner_names)
        self.gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
        self.folds = [list(f) for f in lock['folds']]
        self.halves = dict(lock['halves'])
        self.directions = [tuple(d) for d in lock.get('directions', DIRECTIONS)]
        self.odors = [o for f in self.folds for o in f]
        self.n_boot, self.seed, self.global_maxiter = n_boot, seed, global_maxiter
        self.readings = tuple(readings)
        self.verbose = verbose
        gl = self.gl
        self.fit_t = {h: (lambda train, fl=fl: Target.make(flies_subset(pop, fl), train, gl))
                      for h, fl in self.halves.items()}
        self.read_t = {h: Target.make(flies_subset(pop, fl), [o for o in pop.odors if o in self.odors], gl)
                       for h, fl in self.halves.items()}

    def say(self, msg: str) -> None:
        if self.verbose:
            print(f'{_now()} step 1: {msg}', flush=True)

    def cf(self, fit_one, targets: Optional[Mapping] = None) -> dict:
        return cross_fit(fit_one, self.folds, targets or self.fit_t, self.read_t, self.directions, self.n_boot,
                         self.seed, readings=self.readings)

    def permuted_targets(self, pseed: int) -> Dict[str, Callable[[List[str]], Target]]:
        """Fit targets whose dendrites have their odor labels permuted among each fold's training odors (seed
        `pseed`, the same permutation of each fold's odors for both halves)."""
        gl, pop = self.gl, self.pop
        return {h: (lambda train, fl=fl: Target.make(
            permute_odor_labels(flies_subset(pop, fl), np.random.default_rng(pseed), train), train, gl))
            for h, fl in self.halves.items()}

    def fitters(self, ns: Dict[tuple, dict]):
        """R1, R2m and R2f's fold fits, keeping each fold's state in `ns` (the real fits' or a control's own)."""
        stage, gl = self.stage, self.gl

        def f1(t, train):
            r = fit_global(stage, t, maxiter=self.global_maxiter)
            ns.setdefault(_key(t), {})['R1'] = r['params']
            return stage.run(r['params']), dict(params=r['params'], loss=r['loss'], at_bounds=r['at_bounds'])

        def f2m(t, train):
            s_ = ns[_key(t)]
            meas = measured_for(t, self.orn, self.orn_partners, self.partner_names, seed=self.seed)
            start = measured_gain_start(meas['step'], model_gains(stage, s_['R1'], gl, train))
            s_.update(meas=meas, R2m=Terms(gains=start))
            nf = sorted(str(g) for g in meas['step'].index[~meas['step'].follows.astype(bool)])
            return stage.run(s_['R1'], s_['R2m']), dict(gains=start, not_following=nf)

        def f2f(t, train):
            s_ = ns[_key(t)]
            r = fit_gains(stage, t, s_['R1'], s_['R2m'].gains)
            s_['R2f'] = r['terms']
            return stage.run(s_['R1'], r['terms']), dict(gains=r['terms'].gains, loss=r['loss'],
                                                         at_bounds=r['at_bounds'])
        return f1, f2m, f2f


def run_control(ctx: Step1, i: int) -> dict:
    """R1's and R2's refits to odor-permuted dendrites, number `i` (seed `control_seed(i)`; the permutation among each
    fold's training odors), cross-fitted and read on the real held-out dendrites as the rungs are: the points of R1, R2m
    and R2f in each reading. Reported; they decide nothing (ruling 6)."""
    pseed = control_seed(i, ctx.seed)
    pt = ctx.permuted_targets(pseed)
    c1, c2m, c2f = ctx.fitters({})
    row = dict(perm=int(i), seed=pseed)
    for name, f in (('R1', c1), ('R2m', c2m), ('R2f', c2f)):
        ctx.say(f'control {i}: {name}')
        o = ctx.cf(f, pt)
        row[name] = o['point']
        for rd in ctx.readings[1:]:
            row[f'{name}_{rd}'] = o['by_reading'][rd]['point']
        for d in _directions_of(o):
            row[f'{name}_{d}'] = o[d]['point']
    return row


def run_rungs(stage: PNStage, pop, orn: np.ndarray, orn_partners: np.ndarray, partner_names: Sequence[str],
              lock: Mapping, directions: Optional[Sequence[Tuple[str, str]]] = None, n_boot: int = N_BOOT_KEEP,
              seed: int = SEED, n_perm_lateral: int = N_PERM_LATERAL, n_perm_control: int = 0,
              max_lateral_terms: int = MAX_LATERAL_TERMS, check_fifth: bool = True, global_maxiter: int = 400,
              readings: Sequence[str] = (KEEP_READING, BESIDE_READING), verbose: bool = True) -> dict:
    """Step 1's rungs on the PN stage alone, in order, each fitted on the one below and cross-fitted (`cross_fit`: fold
    by fold on one half of the PN flies, read on the other half on the odors it did not see, both directions).

    R0 the stage as it is; R1 the global parameters (`fit_global`); R2 with the gains, measured (R2m: the measured gain
    over the R1 stage's own, `measured_gain_start`, no fit) and fitted from there (R2f, `fit_gains`), R2 being R2f if
    it is kept against R2m (`keep_rung`), else R2m; R3 with lateral terms, forward-selected among each glomerulus's
    N_CANDIDATES best measured partners in each fold (a term is a glomerulus and a partner rank, so each fold uses its
    own partners; its weight fitted fold by fold, `fit_lateral`), added one at a time while the held-out error falls
    and the term gains more than at least `lateral_min_exceeded(n_perm_lateral)` of `n_perm_lateral` fits of the same
    term to dendrites with odor labels permuted among the fold's training odors, up to `max_lateral_terms` terms (with
    `check_fifth`, whether one more would have passed is recorded, never added).

    A rung is compared with the highest kept rung below it (R0 counts as kept) by `keep_against` in each reading of
    `readings`; the first (KEEP_READING) decides, and a rung kept in one reading and not the other is
    "reading-dependent". `n_perm_control` odor-permuted refits of R1 and R2 (`run_control`, indices 0 to n - 1) run
    here too when asked; step 1 runs them as processes of their own (`part_step1_control`). `pop` is viki's PN
    population in the fit's units, its odors named as the stage's; `orn` and `orn_partners` the ORN recordings in dF/F;
    `lock` gives 'halves', 'folds', 'glomeruli' (and 'directions').

    Returns:
        dict with 'rungs' (one row per rung: the point and per-direction held-out matrix errors, the view errors, the
        comparator 'below', its SE, 'kept' and the beside reading's, 'reading_dependent'), 'controls', 'selected' (R3's
        terms), 'excluded', 'lateral_steps' (each step's best term, its gain, the permuted gains, how many it
        exceeded, passed), 'fifth', 'r2', 'highest_kept', 'folds' (per rung the folds' records), 'tables' (per rung the
        assembled held-out tables per direction), 'state'.
    """
    lk = dict(lock)
    if directions is not None:
        lk['directions'] = [list(d) for d in directions]
    ctx = Step1(stage, pop, orn, orn_partners, partner_names, lk, n_boot, seed, global_maxiter, readings, verbose)
    gl = ctx.gl
    st: Dict[tuple, dict] = {}
    res: Dict[str, dict] = {}

    def f0(t, train):
        return stage.run(), {}

    f1, f2m, f2f = ctx.fitters(st)
    for name, f in (('R0', f0), ('R1', f1), ('R2m', f2m), ('R2f', f2f)):
        ctx.say(name)
        res[name] = ctx.cf(f)
    r2f_vs_r2m = keep_against(res['R2f'], res['R2m'], readings[0])
    r2 = 'R2f' if r2f_vs_r2m['kept'] else 'R2m'
    for k in st:
        st[k]['R2'] = st[k][r2]
        st[k]['R3'] = Terms(gains=dict(st[k]['R2'].gains), laterals=[])

    def partner(k, g, rank):
        cp = st[k]['meas']['partners']
        x = cp[(cp.glomerulus == g) & (cp['rank'] == rank)]
        return None if not len(x) else str(x.partner.iloc[0])

    def lateral_fit(term, store=None):
        g, rank = term

        def f(t, train):
            s_ = st[_key(t)]
            p_ = partner(_key(t), g, rank)
            if p_ is None or p_ not in stage.glomeruli:
                return stage.run(s_['R1'], s_['R3']), dict(partner=None)
            r = fit_lateral(stage, t, s_['R1'], s_['R3'], g, p_)
            if store is not None:
                store[_key(t)] = r['terms']
            return stage.run(s_['R1'], r['terms']), dict(partner=p_, weight=r['weight'], loss=r['loss'])
        return f

    current = res[r2]
    selected: List[tuple] = []
    excluded: List[tuple] = []
    steps_rec: List[dict] = []
    fifth: Optional[dict] = None
    cands = [(g, r) for g in gl for r in range(1, N_CANDIDATES + 1)]
    need = lateral_min_exceeded(n_perm_lateral)
    prng = np.random.default_rng([int(seed), 2])
    cache: Optional[Dict[tuple, tuple]] = None              # the candidates' cross-fits while the state is unchanged
    while True:
        dry = len(selected) >= max_lateral_terms
        if dry and not check_fifth:
            break
        n_of = {g: sum(1 for x in selected if x[0] == g) for g in gl}
        trial = [c for c in cands if c not in selected and c not in excluded and n_of[c[0]] < MAX_LATERAL]
        if not trial:
            if dry:
                fifth = dict(would_pass=False, why='no candidate term left')
            break
        if cache is None:
            ctx.say(f'R3 step {len(selected) + 1}{" (the fifth, checked only)" if dry else ""}: {len(trial)} '
                    'candidate terms')
            cache = {}
            for c in trial:
                store: Dict[tuple, Terms] = {}
                cache[c] = (ctx.cf(lateral_fit(c, store=store)), store)
        best = min(trial, key=lambda c: cache[c][0]['point'])
        best_out, best_store = cache[best]
        gain = current['point'] - best_out['point']
        if gain <= 0:
            if dry:
                fifth = dict(would_pass=False, why='no candidate lowers the held-out error')
            break
        perm_gains, perm_seeds = [], []
        for _ in range(n_perm_lateral):
            pseed = int(prng.integers(1 << 30))
            out_p = ctx.cf(lateral_fit(best), ctx.permuted_targets(pseed))
            perm_gains.append(current['point'] - out_p['point'])
            perm_seeds.append(pseed)
        n_exc = int(sum(gain > pg for pg in perm_gains))
        passes = n_exc >= need
        steps_rec.append(dict(step=len(selected) + 1, glomerulus=best[0], rank=best[1], gain=float(gain),
                              point=float(best_out['point']), perm_gains=[float(x) for x in perm_gains],
                              perm_seeds=perm_seeds, n_exceeded=n_exc, needed=need, passes=bool(passes),
                              checked_only=bool(dry)))
        if not passes:
            excluded.append(best)
            ctx.say(f'R3: {best} gains {gain:.4f}, more than {n_exc} of {n_perm_lateral} odor-permuted fits '
                    f'(needs {need}); left out')
            continue
        if dry:
            fifth = dict(would_pass=True, glomerulus=best[0], rank=best[1], gain=float(gain), n_exceeded=n_exc)
            ctx.say(f'R3: a fifth term {best} would have passed (gain {gain:.4f}); not added')
            break
        selected.append(best)
        for k, terms in best_store.items():
            st[k]['R3'] = terms
        current = best_out
        cache = None
        ctx.say(f'R3: {best} added, held-out error {current["point"]:.4f} (gain {gain:.4f}; more than {n_exc} of '
                f'{n_perm_lateral} permuted)')
    res['R3'] = current if selected else res[r2]

    controls = [run_control(ctx, i) for i in range(int(n_perm_control))]

    # each rung against the highest kept rung below it (R0 counts as kept), in each reading; the first decides
    rows, highest = [], {rd: 'R0' for rd in readings}
    below_of: Dict[str, Dict[str, str]] = {}

    def compare(r, against=None):
        cmp = {}
        for rd in readings:
            b = highest[rd] if against is None else against
            cmp[rd] = dict(below=b, **keep_against(res[r], res[b], rd))
        return cmp

    cmps: Dict[str, dict] = {}
    cmps['R1'] = compare('R1')
    for rd in readings:
        if cmps['R1'][rd]['kept']:
            highest[rd] = 'R1'
    cmps['R2m'] = compare('R2m')
    cmps['R2f'] = compare('R2f')
    vs_r2m = {rd: keep_against(res['R2f'], res['R2m'], rd) for rd in readings}
    r2_by = {rd: ('R2f' if vs_r2m[rd]['kept'] else 'R2m') for rd in readings}
    for rd in readings:
        if cmps[r2_by[rd]][rd]['kept']:
            highest[rd] = r2_by[rd]
    if selected:
        cmps['R3'] = compare('R3')
        for rd in readings:
            if cmps['R3'][rd]['kept']:
                highest[rd] = 'R3'
    for r in ['R0', 'R1', 'R2m', 'R2f', 'R3']:
        o = res[r]
        row = dict(rung=r, point=o['point'])
        for rd in readings[1:]:
            row[f'point_{rd}'] = o['by_reading'][rd]['point']
        for d in _directions_of(o):
            row[f'point_{d}'] = o[d]['point']
            sc = o[d]['score']
            row.update({f'{d}_{k}': v for k, v in sc.as_series().items() if k.startswith(('error', 'view_'))})
        if r in cmps:
            c0 = cmps[r][readings[0]]
            row.update(below=c0['below'], se_vs_below=c0['se'], diff_vs_below=c0['diff'], kept=c0['kept'],
                       no_higher_each_direction=c0['no_higher_each_direction'])
            for rd in readings[1:]:
                row.update({f'below_{rd}': cmps[r][rd]['below'], f'kept_{rd}': cmps[r][rd]['kept']})
            row['reading_dependent'] = bool(len({cmps[r][rd]['kept'] for rd in readings}) > 1)
        elif r == 'R0':
            row.update(kept=True)
        if r == 'R2f':
            row.update(kept_vs_R2m=vs_r2m[readings[0]]['kept'], se_vs_R2m=vs_r2m[readings[0]]['se'],
                       **{f'kept_vs_R2m_{rd}': vs_r2m[rd]['kept'] for rd in readings[1:]})
        if r == 'R3' and not selected:
            row.update(kept=False, note='no lateral term passed; R3 is R2')
        row['is_r2'] = r == r2
        rows.append(row)
    rungs = pd.DataFrame(rows)
    return dict(rungs=rungs, controls=pd.DataFrame(controls), selected=[list(x) for x in selected],
                excluded=[list(x) for x in excluded], lateral_steps=steps_rec, fifth=fifth, r2=r2,
                r2_by_reading=r2_by, highest_kept=highest[readings[0]], highest_kept_by_reading=highest,
                folds={r: {d: res[r][d]['folds'] for d in _directions_of(res[r])} for r in res},
                tables={r: {d: res[r][d]['table'] for d in _directions_of(res[r])} for r in res},
                state={str(k): dict(R1=v.get('R1'), R2=v['R2'].as_dict() if 'R2' in v else None,
                                    R3=v['R3'].as_dict() if 'R3' in v else None) for k, v in st.items()})


def step2_rung(highest_kept: str) -> Optional[str]:
    """The rung step 2 runs from step 1's highest kept rung (ruling 6): 'R2' (the fits' measured or fitted gains, as
    step 1 chose) or 'R3', or None when R0 or R1 is the highest kept (the rung's 24 runs are dropped; R1 runs anyway)."""
    return {'R2m': 'R2', 'R2f': 'R2', 'R3': 'R3'}.get(highest_kept)


def step1_inputs(lock: Mapping, panel: str = 'megamat') -> Tuple[Step1, PNStage, dict]:
    """The family's PN stage on `panel` (`family_stage` for step 1: after E134's gate or Tom's authorization), viki's PN population in the lock's fit
    units, and the ORN recordings in dF/F (ruling 1) on the lock's glomeruli and on the partner pool: (`Step1` context,
    the stage, `family_stage`'s info)."""
    stage, info = family_stage(panel, purpose='step1')
    units = lock.get('units', {})
    pop = pn_population(units.get('fit_pn', lock.get('fit_units', FIT_UNITS)))
    orn_t = recorded_stage(ORN_KEY, units.get('fit_orn', FIT_ORN_UNITS))
    gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
    model_gl = list(H.apl16_orn_deltas().index)
    pool = [g for g in model_gl if g in set().union(*[set(t.index) for t in orn_t.values()]) and g in stage.glomeruli]
    ctx = Step1(stage, pop, fly_stack(orn_t, gl, MEGAMAT_ODORS), fly_stack(orn_t, pool, MEGAMAT_ODORS), pool, lock)
    return ctx, stage, info


def part_step1() -> dict:
    """Step 1, the real rungs (one lane; they alone decide step 2): fitted in order to the dendrites and read on the
    held-out flies and odors, R3's lateral terms with their permuted fits (the lock's count), and the rungs refitted on
    all seven PN flies for step 2 with their scale factors (`rung_fits_all_flies`). The permuted refits of R1 and R2 run
    as processes of their own (`part_step1_control`). Refuses before E134's gate and this experiment's committed lock.
    Writes step1_rungs.csv, step1_lateral.csv, step1_tables.csv, step1_fits.json and provenance_step1.json. Runs when
    E134's gate is usable or on Tom's authorization (`require_step1_allowed`), which the provenance records."""
    allowed = require_step1_allowed()
    lock = read_lock()
    t0 = time.time()
    ctx, stage, info = step1_inputs(lock)
    perms = lock.get('permutations', {})
    out = run_rungs(stage, ctx.pop, ctx.orn, ctx.orn_partners, ctx.partner_names, lock,
                    n_perm_lateral=int(perms.get('lateral', N_PERM_LATERAL)), n_perm_control=0,
                    max_lateral_terms=int(lock.get('r3', {}).get('max_terms', MAX_LATERAL_TERMS)))
    d = out_dir()
    out['rungs'].to_csv(d / 'step1_rungs.csv', index=False)
    pd.DataFrame(out['lateral_steps']).to_csv(d / 'step1_lateral.csv', index=False)
    tabs = []
    for r, by in out['tables'].items():
        for direction, t in by.items():
            x = t.stack().rename('value').reset_index()
            x.columns = ['glomerulus', 'odor', 'value']
            tabs.append(x.assign(rung=r, direction=direction))
    pd.concat(tabs, ignore_index=True).to_csv(d / 'step1_tables.csv', index=False)
    fits = rung_fits_all_flies(stage, ctx.pop, ctx.orn, ctx.orn_partners, ctx.partner_names, lock, out)
    _write_json(d / 'step1_fits.json', dict(folds=out['folds'], state=out['state'], selected=out['selected'],
                                            excluded=out['excluded'], fifth=out['fifth'], r2=out['r2'],
                                            highest_kept=out['highest_kept'],
                                            highest_kept_by_reading=out['highest_kept_by_reading'],
                                            step2_rung=step2_rung(out['highest_kept']), all_flies=fits,
                                            scale=scale_record(fits)))
    _write_json(d / 'provenance_step1.json', dict(script=H.script_path(__file__), part='step1', host=_host(),
                time=_now(), seconds=time.time() - t0, lock_sha256=_sha256(LOCK_COPY), identity=info['identity'],
                allowed=allowed, model_provenance=H.model_provenance(orn_deltas=stage.delta), **H.code_commits()))
    return out


def part_step1_control(perms: Sequence[int], lane: Optional[str] = None, reclaim: bool = False) -> pd.DataFrame:
    """R1's and R2's odor-permuted refits by index (`run_control`), each a process of its own on any lane (claimed by
    mkdir of claims/step1_control_<i>/); writes step1_control_<i>.json and gathers every one written so far into
    step1_controls.csv (`collect_controls`). Refuses before E134's gate (or Tom's authorization of step 1,
    `require_step1_allowed`) and the committed lock, and for an index the lock does not count."""
    allowed = require_step1_allowed()
    lock = read_lock()
    n = int(lock.get('permutations', {}).get('control', N_PERM_CONTROL))
    bad = [i for i in perms if not 0 <= int(i) < n]
    if bad:
        raise ValueError(f'the lock counts {n} controls (0 to {n - 1}); not {bad}')
    ctx = None
    lane = lane or _host()
    for i in perms:
        tag = f'step1_control_{int(i)}'
        if (out_dir() / f'{tag}.json').exists() or not claim(tag, lane, reclaim):
            print(f'{tag}: done or claimed elsewhere, skipped', flush=True)
            continue
        try:
            t0 = time.time()
            if ctx is None:
                ctx, stage, info = step1_inputs(lock)
            row = run_control(ctx, int(i))
            _write_json(out_dir() / f'{tag}.json', dict(row=row, host=_host(), lane=lane, time=_now(),
                                                        seconds=time.time() - t0, lock_sha256=_sha256(LOCK_COPY),
                                                        allowed=allowed, **H.code_commits()))
        except BaseException:
            release(tag)
            raise
    return collect_controls()


def collect_controls(root: Optional[Path] = None) -> pd.DataFrame:
    """Every step1_control_<i>.json written so far, as step1_controls.csv (one row per permutation)."""
    d = out_dir() if root is None else Path(root)
    rows = [_read_json(p)['row'] for p in sorted(d.glob('step1_control_*.json'))]
    t = pd.DataFrame(rows)
    if len(t):
        t = t.sort_values('perm')
        t.to_csv(d / 'step1_controls.csv', index=False)
    return t


def family_stage(panel: str = 'megamat', wiring: str = 'real', purpose: str = 'step2') -> Tuple[PNStage, dict]:
    """The family's PN stage on a panel's input (`panel_orn`: megamat, apl16, validation2, mixtures): E83's
    `build_simulator(wiring, orn=..., breadth=<the lab input's>)`, its olfsysm variables read into a `PNStage` (rows
    the wiring's glomeruli, columns the input's odors by `odor_labels`, in olfsysm's order through the simulator's
    'odor_order'), checked against olfsysm's own layers (`check_identity`) on a copy of the build's variables. A model
    build: for `purpose` 'step1' after E134's gate or Tom's authorization of step 1 (`require_step1_allowed`), for any
    other purpose (step 2, E136) after E134's gate only (`require_e134_gate`)."""
    if purpose == 'step1':
        require_step1_allowed()
    else:
        require_e134_gate()
    from experiments import E83_cross_wiring_fit as X
    orn = panel_orn(panel)
    sim, info = X.build_simulator(wiring, orn=orn, breadth=X.wiring_parts(wiring)['breadth'])
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    names = [None] * len(order)
    labels = odor_labels(orn.columns)
    for j, i in enumerate(order):
        names[int(i)] = labels[j]
    gl = [str(g) for g in dict.fromkeys(info['W'].columns.get_level_values('glomerulus'))]
    stage = PNStage.from_mp(mp, gl, names)
    want = orn.copy()
    want.columns = labels
    want = want.reindex(index=gl, columns=names).fillna(0.0)
    if not np.allclose(stage.delta.to_numpy(), want.to_numpy(), atol=1e-12):
        raise RuntimeError(f"{wiring} {panel}: olfsysm's delta is not the panel's input on the wiring's glomeruli")
    ident = stage.check_identity(mp, rv)
    return stage, dict(sim=sim, info=info, mp=mp, rv=rv, identity=ident, orn=orn)


def panel_orn(panel: str) -> pd.DataFrame:
    """The model's ORN input of a panel: 'megamat' E107's `panel_input('megamat')` (17 odors on the lab input's 36
    glomeruli), 'apl16' the lab's input (`apl16_orn_deltas`), 'validation2' and 'mixtures' E134's `panel_input`."""
    if panel == 'megamat':
        from experiments import E107_validation2_panel as E107
        return E107.panel_input('megamat')
    if panel == 'apl16':
        return H.apl16_orn_deltas()
    if panel in ('validation2', 'mixtures'):
        from experiments import E134_population_baseline as E134R
        return E134R.panel_input(panel)
    raise ValueError(f'{panel=}')


def rung_fits_all_flies(stage: PNStage, pop, orn: np.ndarray, orn_partners: np.ndarray,
                        partner_names: Sequence[str], lock: Mapping, rungs_out: Mapping) -> dict:
    """The rungs refitted on all seven PN flies for step 2 (ruling 5): fold by fold (each megamat odor's table from the
    fit that left it out) and once on all 17 odors (for the APL16 input's paa, euc, -aPine and B-myr); R1's
    parameters, R2's gains (measured or fitted, as step 1 chose), R3's selected terms (each fold's own partner of each
    selected glomerulus and rank). `stage` is the megamat stage: each fit's scale factors are set on it
    (`rung_scale_factor`) and recorded under 'scale' with the fit's training odors."""
    gl = [g for g in lock['glomeruli'] if g in stage.glomeruli]
    folds = [list(f) for f in lock['folds']]
    odors = [o for f in folds for o in f]
    sets = [(k, [o for o in odors if o not in set(f)]) for k, f in enumerate(folds)] + [('all', odors)]
    base = stage.run(reading=SCALE_READING)
    out = {}
    for k, train in sets:
        t = Target.make(pop, train, gl)
        r1 = fit_global(stage, t)['params']
        meas = measured_for(t, orn, orn_partners, partner_names)
        start = measured_gain_start(meas['step'], model_gains(stage, r1, gl, train))
        terms = Terms(gains=start)
        if rungs_out['r2'] == 'R2f':
            terms = fit_gains(stage, t, r1, start)['terms']
        r2 = terms
        for g, rank in rungs_out['selected']:
            cp = meas['partners']
            x = cp[(cp.glomerulus == g) & (cp['rank'] == rank)]
            if len(x) and str(x.partner.iloc[0]) in stage.glomeruli:
                terms = fit_lateral(stage, t, r1, terms, g, str(x.partner.iloc[0]))['terms']
        f = dict(held_out=[o for o in odors if o not in set(train)], train=list(train), R1=r1, R2=r2.as_dict(),
                 R3=terms.as_dict(),
                 not_following=sorted(str(g) for g in meas['step'].index[~meas['step'].follows.astype(bool)]))
        f['scale'] = {r: rung_scale_factor(stage, f, r, train, gl, base=base) for r in ('R1', 'R2', 'R3')}
        out[str(k)] = f
    return out


def rung_scale_factor(stage: PNStage, fit: Mapping, rung: str, train: Sequence[str], glomeruli: Sequence[str],
                      reading: str = SCALE_READING, base: Optional[pd.DataFrame] = None) -> float:
    """The one positive factor that holds a fit's table at the model's own rate scale (ruling 5): the pooled SD of the
    unchanged stage's evoked rates over `glomeruli` (the 28 shared) and the fit's training odors over the fit's own
    there. Set on the megamat input and applied to the fit's APL16 table too."""
    b = (stage.run(reading=reading) if base is None else base).loc[list(glomeruli), list(train)].to_numpy(dtype=float)
    x = stage.run(fit['R1'], None if rung == 'R1' else Terms.from_dict(fit[rung]), reading=reading)
    x = x.loc[list(glomeruli), list(train)].to_numpy(dtype=float)
    sx = float(np.nanstd(x))
    if not sx > 0:
        raise ValueError(f'{rung}: the fit\'s rates do not vary over the shared glomeruli and training odors')
    return float(np.nanstd(b) / sx)


def scale_record(fits: Mapping[str, dict], rungs: Sequence[str] = ('R1', 'R2', 'R3')) -> dict:
    """The scale factors of every fit and rung, and those outside SCALE_FACTOR_RANGE, named (ruling 5)."""
    lo, hi = SCALE_FACTOR_RANGE
    fac = {r: {k: float(f['scale'][r]) for k, f in fits.items() if 'scale' in f and r in f['scale']} for r in rungs}
    named = {r: sorted(k for k, v in d.items() if not lo <= v <= hi) for r, d in fac.items()}
    return dict(factors=fac, range=[lo, hi], named=named)


def rung_table(stage: PNStage, fits: Mapping[str, dict], rung: str, reading: str = 'preset',
               scaled: bool = True, fold_odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """A rung's table of PN rates on a stage's input for step 2 (`rung_fits_all_flies`' fits): each odor a fold held
    out from that fold's fit, every other odor (the APL16 input's paa, euc, -aPine, B-myr) from the fit on all odors;
    rates in `reading`, each fit's table times its scale factor for `rung` (set on the megamat input, ruling 5) unless
    not `scaled`. `fold_odors` limits the odors taken from the fold fits (default every fold's held-out odors; () takes
    every odor from the fit on all odors)."""
    keep = None if fold_odors is None else set(fold_odors)
    held = {o: k for k, f in fits.items() if k != 'all' for o in f['held_out'] if keep is None or o in keep}
    cache: Dict[str, pd.DataFrame] = {}
    cols = {}
    for o in stage.odors:
        k = held.get(o, 'all')
        if k not in fits:
            continue
        if k not in cache:
            f = fits[k]
            x = stage.run(f['R1'], None if rung == 'R1' else Terms.from_dict(f[rung]), reading=reading)
            if scaled:
                if rung not in f.get('scale', {}):
                    raise KeyError(f'fit {k}: no scale factor for {rung} (rung_fits_all_flies sets them)')
                x = x * float(f['scale'][rung])
            cache[k] = x
        cols[o] = cache[k][o]
    return pd.DataFrame(cols)[[o for o in stage.odors if o in cols]]


def odor_labels(columns: Sequence) -> List[str]:
    """An input's odor labels as the stage names them: abbreviated (`H.odor_abbrev`) where the abbreviations are
    unique (the megamat, APL16 and validation2 inputs), else the full names (the mixtures, whose stimuli share odors at
    several concentrations)."""
    ab = [H.odor_abbrev(c) for c in columns]
    return ab if len(set(ab)) == len(ab) else [str(c) for c in columns]


def stage_on_input(stage: PNStage, orn: pd.DataFrame) -> PNStage:
    """The family's PN stage with another ORN input: `orn` (glomeruli x odors of deltas, e.g. E134's `panel_input`)
    on the stage's glomeruli (glomeruli the input lacks at 0, as olfsysm's delta holds them), columns named by
    `odor_labels`, with the stage's spontaneous rates and layer parameters. Those do not depend on the panel: the
    model's spontaneous rates are Hallem's per glomerulus on the wiring's glomeruli (fit_mb_model's sfr), and the
    family's ORN, LN and PN parameters are olfsysm's defaults on every panel."""
    d = orn.copy()
    d.columns = odor_labels(orn.columns)
    d.index = d.index.astype(str)
    d = d.reindex(index=stage.glomeruli).fillna(0.0)
    return PNStage(d, pd.Series(stage.spont, index=stage.glomeruli), stage.sp)


def kept_stage_table(fits: Mapping[str, dict], rung: str, stage: PNStage, orn: Optional[pd.DataFrame] = None,
                     panel: Optional[str] = None, fold_odors: Optional[Sequence[str]] = None) -> Tuple[pd.DataFrame, dict]:
    """E135's PN stage as a table of PN rates on any ORN input, as step 2's route 1 takes it (for E136's arm A1 on
    validation2 and the mixtures as much as on megamat and APL16).

    The table is in `pn_evoked_rates`' reading ('preset'), each fit's rates times its scale factor (ruling 5: set once
    per fit on the megamat input), and clipped to the stage's `pn_rate_bounds` on that input, BOUNDS_MARGIN_HZ inside
    (ruling 7). Its columns are the stage's odor labels (`odor_labels`: abbreviations where unique), which
    fit_mb_model's `pn_rates` does not take: pass it as `extra_kws=dict(pn_rates=input_columns(route1_pn_rates(table,
    list(table.index)), orn))` (the input's own columns) to E83's `build_simulator` (or E134's `build_point`) on the
    same input `orn`; the odor's LN stays the model's.

    Args:
        fits: step 1's fits on all seven PN flies (step1_fits.json's 'all_flies', `rung_fits_all_flies`).
        rung: 'R1', or the kept rung: 'R2' / 'R3' (step1_fits.json's 'step2_rung'; 'R2m' and 'R2f' read as 'R2').
        stage: the family's PN stage on its wiring (`family_stage(panel)` after E134's gate, or any panel's build of
            the same wiring when `orn` is given).
        orn: the ORN input (glomeruli x odors of deltas), put on `stage` (`stage_on_input`); default `stage`'s own.
        panel: the input's panel: on 'megamat' and 'apl16' each megamat odor comes from the fold fit that left it out
            (the megamat odors at the fits' 10^-3), every other odor from the fit to all 17 odors; on any other panel
            (validation2, the mixtures) every odor from the fit to all 17 odors.
        fold_odors: overrides which odors come from their fold fits.

    Returns:
        (the table: the stage's glomeruli x the input's odors, Hz; the record: rung, panel, the odors from fold fits,
        the scale factors (`scale_record`, those outside 0.5 to 2 named) and the clip's (`clip_table`)).
    """
    r = {'R2m': 'R2', 'R2f': 'R2'}.get(rung, rung)
    if r not in ('R1', 'R2', 'R3'):
        raise ValueError(f'{rung=}: R1, R2 or R3')
    st = stage if orn is None else stage_on_input(stage, orn)
    if fold_odors is None:
        fold_odors = list(MEGAMAT_ODORS) if panel in ('megamat', 'apl16') else []
    tab = rung_table(st, fits, r, fold_odors=fold_odors)
    tab, crec = clip_table(tab, st.rate_bounds(), list(st.glomeruli), None, CLIP_POLICY)
    held = {o for k, f in fits.items() if k != 'all' for o in f['held_out']}
    used = sorted(o for o in tab.columns if o in set(fold_odors) and o in held)
    return tab, dict(rung=r, panel=panel, fold_odors=used, scale=scale_record(fits, (r,)), clip=crec,
                     reading='preset', margin_hz=BOUNDS_MARGIN_HZ)


# ---------------------------------------------------------------------------
# the timing helper and the budget (before the lock)

def planted_step1(size: Optional[Mapping[str, int]] = None, seed: int = SEED, n_orn_flies: int = 9) -> dict:
    """Planted data of the real size for timing step 1 (`time_step1`): a PN stage (olfsysm's default layers, as the
    family's) on a random input of size['n_glomeruli'] glomeruli, the first size['n_input'] with input and the first
    size['n_shared'] of those shared; size['n_pn_flies'] PN flies of size['n_trials'] trials made from its rates at
    planted gains (an affine map plus trial noise, in viki's layout and the trial-noise units' footing), the ORN
    recordings of `n_orn_flies` flies (the input's deltas times a factor plus noise), the halves (seeded, as the lock's)
    and the folds of one odor. Returns dict(stage, pop, orn, orn_partners, partner_names, lock)."""
    from data_prep import stage_population_data as D
    z = dict(REAL_SIZE, **dict(size or {}))
    rng = np.random.default_rng(seed)
    G, Gi, Gs, O = z['n_glomeruli'], z['n_input'], z['n_shared'], z['n_odors']
    gl = [f'g{i:02d}' for i in range(G)]
    odors = list(MEGAMAT_ODORS)[:O]
    delta = np.zeros((G, O))
    delta[:Gi] = rng.uniform(-20, 150, (Gi, O))
    stage = PNStage(pd.DataFrame(delta, index=gl, columns=odors), pd.Series(rng.uniform(5, 30, G), index=gl))
    shared = gl[:Gs]
    true = Terms(gains={g: float(v) for g, v in zip(shared, np.exp(rng.normal(0.0, 0.5, Gs)))})
    rates = stage.run(terms=true).loc[shared, odors].to_numpy()
    rows = []
    for f in range(z['n_pn_flies']):
        X = (0.05 * rates + 0.2)[None] + rng.normal(scale=0.15, size=(z['n_trials'], Gs, O))
        rows.append(D.long_from_trials(f'f{f}', shared, odors, X))
    spec = dict(units=D.UNITS_TRIAL, criterion=2.0, window={}, notes=[], recording_site='dendrites (planted)')
    pop = D.population_from_long(PN_KEY, pd.concat(rows, ignore_index=True), spec=spec, trials=True)
    orn_all = (0.01 * delta[:Gi])[None] + rng.normal(scale=0.05, size=(n_orn_flies, Gi, O))
    flies = [str(f) for f in pop.flies]
    A, B = halves(flies, seed)
    lock = dict(glomeruli=shared, folds=odor_folds(O, odors, seed), halves=dict(half_a=A, half_b=B),
                directions=[list(d) for d in DIRECTIONS])
    return dict(stage=stage, pop=pop, orn=orn_all[:, :Gs], orn_partners=orn_all, partner_names=gl[:Gi], lock=lock,
                size=z)


def time_step1(size: Optional[Mapping[str, int]] = None, seed: int = SEED, global_maxiter: int = 400,
               case: Optional[dict] = None) -> dict:
    """One fit of each kind of step 1, timed on planted data of the real size (`planted_step1`; the rulings' "time one
    R1 fold fit on planted data of the real size before writing the lock"): one R1 fold fit (the larger half's flies,
    the training odors of one fold, `fit_global` at `global_maxiter`), the measured step of a fold (`measured_for` and
    the model's gains), one gain fit (`fit_gains`), one lateral fit (`fit_lateral`), one target (`Target.make`), one
    stage run, and one read of an assembled table in both readings (`read_prediction`, `held_out_draws`). Returns
    dict(unit_seconds, nfev_R1, size, threads)."""
    c = planted_step1(size, seed) if case is None else case
    ctx = Step1(c['stage'], c['pop'], c['orn'], c['orn_partners'], c['partner_names'], c['lock'], verbose=False,
                global_maxiter=global_maxiter)
    stage, gl = ctx.stage, ctx.gl
    train = [o for f in ctx.folds[1:] for o in f]
    big = max(ctx.halves, key=lambda h: len(ctx.halves[h]))
    u: Dict[str, float] = {}
    t0 = time.time()
    t = ctx.fit_t[big](train)
    u['target'] = time.time() - t0
    t0 = time.time()
    stage.run()
    u['run'] = time.time() - t0
    t0 = time.time()
    r1 = fit_global(stage, t, maxiter=global_maxiter)
    u['R1'] = time.time() - t0
    t0 = time.time()
    meas = measured_for(t, ctx.orn, ctx.orn_partners, ctx.partner_names)
    start = measured_gain_start(meas['step'], model_gains(stage, r1['params'], gl, train))
    u['measured'] = time.time() - t0
    t0 = time.time()
    g = fit_gains(stage, t, r1['params'], start)
    u['gains'] = time.time() - t0
    t0 = time.time()
    fit_lateral(stage, t, r1['params'], g['terms'], gl[0], gl[1])
    u['lateral'] = time.time() - t0
    read_half = [h for h in ctx.halves if h != big][0]
    tab = stage.run(r1['params'], g['terms'])
    t0 = time.time()
    for name in ctx.readings:
        rd = read_reading(name)
        read_prediction(tab, ctx.read_t[read_half], rd)
        held_out_draws(tab, ctx.read_t[read_half], ctx.n_boot, ctx.seed, rd)
    u['read'] = time.time() - t0
    return dict(unit_seconds={k: float(v) for k, v in u.items()}, nfev_R1=int(r1['nfev']), size=dict(c['size']),
                global_maxiter=int(global_maxiter), threads=os.environ.get('OMP_NUM_THREADS'),
                thread_limit=os.environ.get('OMP_THREAD_LIMIT'))


def project_step1(unit: Mapping[str, float], n_folds: int = N_FOLDS, n_perm_control: int = N_PERM_CONTROL,
                  n_perm_lateral: int = N_PERM_LATERAL, max_terms: int = MAX_LATERAL_TERMS, check_fifth: bool = True,
                  n_candidates: int = N_CANDIDATES * N_SHARED, n_directions: int = len(DIRECTIONS)) -> dict:
    """Step 1's projected cost in lane-hours from the unit times of `time_step1`, at the worst case of R3 (every step
    to `max_terms` terms and the fifth's check, each step's best term passing its permutations at once; an excluded
    term adds its permutations again).

    Per cross-fit: n_directions x n_folds fold fits (each a target, its fit and a stage run) and n_directions reads (in
    both readings). The real lane: R0, R1, R2m, R2f; R3's steps (every candidate cross-fitted, then n_perm_lateral
    permuted cross-fits of the best, each fold making two targets); the fits on all seven flies (n_folds + 1 sets of
    R1, the measured step, the gains and max_terms lateral fits). Each control: R1, R2m and R2f cross-fitted on
    permuted targets. Returns the lane-hours of the real lane, of each control, the total and the breakdown (s)."""
    u = dict(unit)
    F = n_directions * int(n_folds)
    rd = n_directions * u['read']
    fold = lambda fit: F * (u['target'] + u[fit] + u['run']) + rd if fit else F * (u['target'] + u['run']) + rd
    br = dict(R0=fold(None), R1=fold('R1'), R2m=fold('measured'), R2f=fold('gains'))
    steps = int(max_terms) + (1 if check_fifth else 0)
    r3 = 0.0
    for k in range(steps):
        nc = max(int(n_candidates) - k, 0)
        r3 += nc * (F * (u['target'] + u['lateral'] + u['run']) + rd)
        r3 += int(n_perm_lateral) * (F * (2 * u['target'] + u['lateral'] + u['run']) + rd)
    br['R3'] = r3
    br['all_flies'] = (int(n_folds) + 1) * (u['target'] + u['R1'] + u['measured'] + u['gains']
                                            + int(max_terms) * u['lateral'])
    real = sum(br.values())
    control = fold('R1') + fold('measured') + fold('gains')
    total = real + int(n_perm_control) * control
    return dict(n_folds=int(n_folds), n_perm_control=int(n_perm_control), n_perm_lateral=int(n_perm_lateral),
                max_terms=int(max_terms), check_fifth=bool(check_fifth), real_lane_hours=real / 3600,
                control_lane_hours=control / 3600, total_lane_hours=total / 3600,
                breakdown_hours={k: v / 3600 for k, v in br.items()}, fold_fits_per_cross_fit=F)


def choose_budget(unit: Mapping[str, float], budget: float = BUDGET_LANE_HOURS,
                  steps: Sequence[Tuple[int, int]] = BUDGET_STEPS) -> dict:
    """The rulings' fallback on cost only, decided before the lock: the first of (controls, folds) in `steps` (5 and
    17, then 3 and 17, then 3 and 6) whose projected step 1 (`project_step1`) is within `budget` lane-hours; the last if
    none is (flagged). Returns the choice, its projection and every projection tried."""
    tried = []
    for nc, nf in steps:
        p = project_step1(unit, nf, nc)
        tried.append(p)
        if p['total_lane_hours'] <= budget:
            break
    last = tried[-1]
    return dict(n_perm_control=last['n_perm_control'], n_folds=last['n_folds'],
                projected_lane_hours=last['total_lane_hours'], within_budget=bool(last['total_lane_hours'] <= budget),
                budget_lane_hours=float(budget), projection=last, tried=tried, unit_seconds=dict(unit))


def part_timing(size: Optional[Mapping[str, int]] = None, seed: int = SEED) -> dict:
    """The timing helper (`time_step1`, planted data only: no recording and no model of the experiment), its
    projections and the budget's choice (`choose_budget`), written to results/.../timing_step1.json for `write_lock`."""
    t0 = time.time()
    tm = time_step1(size, seed)
    ch = choose_budget(tm['unit_seconds'])
    rec = dict(experiment=EXPERIMENT, part='timing', host=_host(), time=_now(), seconds=time.time() - t0, **tm,
               choice=ch, **H.code_commits())
    _write_json(out_dir() / 'timing_step1.json', rec)
    p = ch['projection']
    print(f"timing: R1 fold fit {tm['unit_seconds']['R1']:.1f} s ({tm['nfev_R1']} evaluations); step 1 projected "
          f"{p['total_lane_hours']:.2f} lane-hours at {ch['n_folds']} folds and {ch['n_perm_control']} controls "
          f"(real lane {p['real_lane_hours']:.2f}, each control {p['control_lane_hours']:.2f})", flush=True)
    return rec


# ---------------------------------------------------------------------------
# the imposed dendrites, both routes

def rank_assign(values: pd.DataFrame, reference: np.ndarray) -> pd.DataFrame:
    """Rank for rank: where `values` and `reference` have as many finite entries, the entry of rank r gets the
    reference's r-th smallest value (ties in `values` share their values' mean), so the result holds exactly the
    reference's values in `values`' arrangement; otherwise hal's `rank_map` (quantile matching)."""
    ref = np.sort(np.asarray(reference, dtype=float)[np.isfinite(reference)])
    v = values.to_numpy(dtype=float)
    ok = np.isfinite(v)
    if ok.sum() != len(ref):
        return P.rank_map(values, ref)
    flat = v[ok]
    order = np.argsort(flat, kind='stable')
    out_flat = np.empty_like(flat)
    out_flat[order] = ref
    ranks = stats.rankdata(flat)
    for r in np.unique(ranks[ranks != np.round(ranks)]):
        tie = ranks == r
        out_flat[tie] = out_flat[tie].mean()
    out = np.full(v.shape, np.nan)
    out[ok] = out_flat
    return pd.DataFrame(out, index=values.index, columns=values.columns)


OUT_OF_BOUNDS_POLICIES: Tuple[Optional[str], ...] = (None, 'clip')


def clip_table(table: pd.DataFrame, bounds: Tuple[pd.DataFrame, pd.DataFrame], glomeruli: Sequence[str],
               odors: Optional[Sequence[str]] = None, policy: Optional[str] = CLIP_POLICY,
               margin: float = BOUNDS_MARGIN_HZ) -> Tuple[pd.DataFrame, dict]:
    """A table of PN rates held inside what each PN can reach (ruling 7): its set entries (`glomeruli` x `odors`,
    default every column; finite entries only) below a PN's floor or above its ceiling (`bounds`, (lo, hi) as olfsysm's
    `pn_rate_bounds` or `PNStage.rate_bounds` give them at the run's LN) are counted and, with `policy` 'clip', moved
    `margin` Hz inside the bounds, on both sides; with `policy` None any such entry raises ValueError.

    Returns:
        (the table, the record: n_set, n_below, n_above, n_clipped, share_clipped, the largest excursion beyond a
        bound in Hz, named (share_clipped above CLIP_NAME_SHARE), the policy and margin, and the entries outside).
    """
    if policy not in OUT_OF_BOUNDS_POLICIES:
        raise ValueError(f'{policy=}')
    gl = [g for g in glomeruli if g in table.index]
    od = list(table.columns if odors is None else odors)
    m = table.loc[gl, od].to_numpy(dtype=float)
    lo = bounds[0].loc[gl, od].to_numpy(dtype=float)
    hi = bounds[1].loc[gl, od].to_numpy(dtype=float)
    ok = np.isfinite(m)
    below, above = ok & (m < lo + margin), ok & (m > hi - margin)
    out_lo, out_hi = ok & (m <= lo), ok & (m >= hi)
    n_set = int(ok.sum())
    rec = dict(n_set=n_set, n_below=int(below.sum()), n_above=int(above.sum()), n_clipped=int((below | above).sum()),
               n_outside=int((out_lo | out_hi).sum()),
               largest_excursion_hz=float(max(np.max(np.where(out_lo, lo - m, 0.0), initial=0.0),
                                              np.max(np.where(out_hi, m - hi, 0.0), initial=0.0))),
               policy=policy, margin_hz=float(margin),
               outside=[dict(glomerulus=gl[i], odor=od[j], rate=float(m[i, j]), lo=float(lo[i, j]), hi=float(hi[i, j]))
                        for i, j in zip(*np.nonzero(out_lo | out_hi))])
    rec['share_clipped'] = rec['n_clipped'] / n_set if n_set else 0.0
    rec['named'] = bool(rec['share_clipped'] > CLIP_NAME_SHARE)
    if rec['n_outside'] and policy is None:
        raise ValueError(f"{rec['n_outside']} rates outside the PN stage's reach (largest by "
                         f"{rec['largest_excursion_hz']:.2f} Hz); choose a policy (policy='clip')")
    out = table.copy().astype(float)
    if policy == 'clip' and rec['n_clipped']:
        out.loc[gl, od] = np.where(ok, np.clip(m, lo + margin, hi - margin), m)
    return out, rec


def imposed_rates(model_rates: pd.DataFrame, dendrites: pd.DataFrame, glomeruli: Sequence[str],
                  odors: Optional[Sequence[str]] = None,
                  bounds: Optional[Tuple[pd.DataFrame, pd.DataFrame]] = None,
                  out_of_bounds: Optional[str] = None, margin: float = BOUNDS_MARGIN_HZ) -> Tuple[pd.DataFrame, dict]:
    """The imposed table: on the shared glomeruli and `odors` (default the dendrites' columns) the dendrites'
    arrangement mapped rank for rank onto the model's own pooled PN rates there (`rank_assign`: exactly the model's values,
    in the dendrites' order); every other entry keeps the model's.

    A mapped rate can lie outside what that glomerulus's PN can reach for that odor (below its floor, about minus its
    spontaneous rate, since the pooled rates reach the most inhibited glomerulus's floor): with `bounds` (lo, hi) such
    entries are counted and `out_of_bounds` says what is done (`clip_table`): None raises ValueError, 'clip' (ruling 7,
    what step 2 passes) moves them `margin` Hz inside the bounds.

    Returns (the table, a record: entries mapped and `clip_table`'s record)."""
    odors = list(dendrites.columns) if odors is None else list(odors)
    gl = list(glomeruli)
    ref = model_rates.loc[gl, odors].to_numpy(dtype=float).ravel()
    mapped = rank_assign(dendrites.loc[gl, odors].astype(float), ref)
    out = model_rates.copy().astype(float)
    out.loc[gl, odors] = mapped.to_numpy()
    rec = dict(n_mapped=int(np.isfinite(mapped.to_numpy()).sum()), n_outside=0, largest_excursion_hz=0.0,
               policy=out_of_bounds)
    if bounds is not None:
        out, crec = clip_table(out, bounds, gl, odors, out_of_bounds, margin)
        rec.update(crec)
    return out, rec


def solve_input(stage: PNStage, target: pd.DataFrame, glomeruli: Optional[Sequence[str]] = None,
                params: Optional[Mapping[str, float]] = None, tol: float = SOLVE_TOL_HZ, max_outer: int = 30,
                n_iter: int = 60, bracket: Tuple[float, float] = (-400.0, 2000.0),
                reading: str = 'preset', delta_tol: float = 1e-3) -> Tuple[pd.DataFrame, dict]:
    """Route 2: the ORN input solved glomerulus by glomerulus so that the unchanged PN stage returns `target`.

    The deltas of `glomeruli` (default the target's rows: the input's glomeruli) are solved, every other row of the
    stage's input kept. The ORN layer is linear in the delta, a PN's evoked rate rises monotonically with its own
    glomerulus's delta at a given LN, and the LN reads only the mean ORN rate over glomeruli: so for each odor the LN
    is fixed by one number, the mean delta over glomeruli. The solve finds, odor by odor, the mean at which the deltas
    solved at that mean's LN (every entry at once, by the Illinois method inside `bracket`) have that mean (a secant
    iteration on the mean), and then checks the stage, re-run exactly with its LN on the solved input, against the
    target: to within `tol` Hz on every solved entry. Rates are in `reading` (default 'preset', the reading of route 1's
    `pn_rates`, so that both routes aim at the same table).

    Returns:
        (the solved input, glomeruli x odors as the stage's delta; the solve's record: outer (secant) iterations, the
        largest difference from the target, entries at the bracket's edge (a target outside the stage's reach), the last
        change of the mean delta, and the largest change of the other glomeruli's rates).
    """
    gl = list(target.index if glomeruli is None else glomeruli)
    od = [o for o in stage.odors if o in target.columns]
    if not od:
        raise ValueError("the target has none of the stage's odors")
    sp = stage.sp.with_global(params)
    G = stage.glomeruli
    rows = np.array([G.index(g) for g in gl])
    cols = np.array([stage.odors.index(o) for o in od])
    T = target.loc[gl, od].to_numpy(dtype=float)
    d = stage.delta.to_numpy(dtype=float).copy()
    base = stage.run(params, reading=reading)
    others = [i for i in range(len(G)) if i not in set(rows)]
    k = orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
    fixed = d[others][:, cols].sum(axis=0) if others else np.zeros(len(cols))
    nG = len(G)
    inhA_full, inhB_full = stage.ln(sp)
    inhA_full, inhB_full = inhA_full.copy(), inhB_full.copy()
    edge = dict(n=0)

    def solved_at(m):
        mean_t = stage.spont.mean() + m[:, None] * k[None, :]
        a, b = ln_layer_mean(sp, mean_t, nG)
        inhA_full[cols], inhB_full[cols] = a, b
        x, n_edge = _solve_entries(sp, stage.spont, rows, cols, T, k, inhA_full, inhB_full, reading, bracket, n_iter,
                                   tol / 100)
        edge['n'] = n_edge
        return (x.sum(axis=0) + fixed) / nG - m, x

    # the mean delta of each odor at which the solved deltas have that mean: F(m) = mean(solved(m)) - m, bracketed
    # outward from the input's own mean, then the Illinois method, every odor at once
    thr = delta_tol / nG
    m0 = d[:, cols].mean(axis=0)
    F0, x0 = solved_at(m0)
    it = 1
    done = np.abs(F0) <= thr
    lo_m, flo, hi_m, fhi = m0.copy(), F0.copy(), m0.copy(), F0.copy()
    step = np.where(np.abs(F0) > 1e-3, F0, np.sign(F0) * 1e-3)
    found = done.copy()
    # (outward on both sides: F falls with m where the solved mean rises more slowly than m, and rises where faster)
    for _ in range(max_outer):
        if found.all():
            break
        for sgn in (1.0, -1.0):
            m1 = np.where(found, m0, m0 + sgn * step)
            F1, _ = solved_at(m1)
            it += 1
            cross = ~found & (np.sign(F1) != np.sign(F0))
            lo_m, flo = np.where(cross, m0, lo_m), np.where(cross, F0, flo)
            hi_m, fhi = np.where(cross, m1, hi_m), np.where(cross, F1, fhi)
            found = found | cross
            if found.all():
                break
        step = np.where(found, step, 2.0 * step)
    x, moved = x0, float(np.abs(F0).max())
    if not done.all():
        side = np.zeros(len(cols))
        for _ in range(max_outer):
            den = fhi - flo
            m = np.where(done, m0, np.where(den != 0, hi_m - fhi * (hi_m - lo_m) / np.where(den != 0, den, 1.0),
                                            0.5 * (lo_m + hi_m)))
            Fm, x = solved_at(m)
            it += 1
            moved = float(np.abs(Fm[~done]).max())
            if moved <= thr:
                break
            same = (np.sign(Fm) == np.sign(flo)) & ~done
            other = ~same & ~done
            lo_m, flo = np.where(same, m, lo_m), np.where(same, Fm, flo)
            hi_m, fhi = np.where(other, m, hi_m), np.where(other, Fm, fhi)
            fhi = np.where(same & (side == -1), fhi / 2, fhi)
            flo = np.where(other & (side == 1), flo / 2, flo)
            side = np.where(same, -1.0, np.where(other, 1.0, side))
    def check(xs):
        dd = d.copy()
        dd[np.ix_(rows, cols)] = xs
        got = PNStage(pd.DataFrame(dd, index=G, columns=stage.odors), pd.Series(stage.spont, index=G), sp).run(
            reading=reading)
        return dd, got, float(np.nanmax(np.abs(got.loc[gl, od].to_numpy() - T)))

    d_best, got, diff = check(x)
    edges, how = int(edge['n']), 'bracketed'
    if diff > tol or edges:
        # no self-consistent mean found cleanly (strong normalization can make the rates nearly blind to the mean):
        # plain alternation from the input's own mean, keeping the deltas whose exact rates are nearest the target
        m = m0.copy()
        for _ in range(max_outer):
            Fm, xm = solved_at(m)
            it += 1
            dd, gg, df = check(xm)
            if edge['n'] == 0 and (df < diff or edges):
                d_best, got, diff, edges, how = dd, gg, df, 0, 'alternation'
            if diff <= tol:
                break
            m = m + Fm
    rec = dict(outer=it, max_abs_diff=diff, at_edge=edges, tol=tol, last_mean_move=moved, method=how,
               others_max_change=float(np.abs(got.iloc[others][od].to_numpy() - base.iloc[others][od].to_numpy()
                                              ).max()) if others else 0.0,
               ln_change=ln_change(stage, pd.DataFrame(d_best, index=G, columns=stage.odors), params, od))
    rec['converged'] = bool(rec['max_abs_diff'] <= tol and rec['at_edge'] == 0)
    rec['fixed_point'] = bool(moved <= thr)
    rec['bracketed'] = bool(found.all())
    return pd.DataFrame(d_best, index=G, columns=stage.odors), rec


def ln_change(stage: PNStage, delta: pd.DataFrame, params: Optional[Mapping[str, float]] = None,
              odors: Optional[Sequence[str]] = None) -> dict:
    """How far the LN moves when the stage's input is replaced by `delta` (route 2's solved input, which the LN reads;
    ruling 9): the largest absolute change of inhA and inhB over the steps and `odors`, and the largest relative change
    of the PN's inhibition factor (pn_inhsc / (pn_inhadd + 0.25 inhA + 0.75 inhB)) over the odor's window."""
    sp = stage.sp.with_global(params)
    cols = [stage.odors.index(o) for o in (stage.odors if odors is None else odors)]
    a0, b0 = stage.ln(sp)
    a1, b1 = stage.ln(sp, delta.reindex(index=stage.glomeruli, columns=stage.odors).to_numpy(dtype=float))
    f = lambda a, b: sp.pn_inhsc / (sp.pn_inhadd + 0.25 * a + 0.75 * b)
    win, _ = windows(sp, 'preset')
    f0, f1 = f(a0[cols], b0[cols])[:, win], f(a1[cols], b1[cols])[:, win]
    return dict(inhA_max_abs=float(np.abs(a1[cols] - a0[cols]).max()),
                inhB_max_abs=float(np.abs(b1[cols] - b0[cols]).max()),
                pn_inh_factor_max_rel=float(np.max(np.abs(f1 / f0 - 1.0))))


def _solve_entries(sp: StageParams, spont: np.ndarray, rows: np.ndarray, cols: np.ndarray, T: np.ndarray,
                   kernel: np.ndarray, inhA: np.ndarray, inhB: np.ndarray, reading: str,
                   bracket: Tuple[float, float], n_iter: int, ftol: float) -> Tuple[np.ndarray, int]:
    """Every entry's own delta at which its PN's evoked rate is the target, at a given LN, by the Illinois method
    (regula falsi with the retained end halved) inside `bracket`, all entries at once; entries whose target lies outside
    the rates the bracket reaches are put at its edge and counted."""
    f = lambda x: _rates_at(sp, spont, rows, cols, x, kernel, inhA, inhB, reading) - T
    lo = np.full(T.shape, float(bracket[0]))
    hi = np.full(T.shape, float(bracket[1]))
    flo, fhi = f(lo), f(hi)
    below, above = flo >= 0, fhi <= 0
    inside = ~(below | above)
    side = np.zeros(T.shape)
    x = 0.5 * (lo + hi)
    for _ in range(n_iter):
        den = fhi - flo
        x = np.where(den > 0, hi - fhi * (hi - lo) / np.where(den > 0, den, 1.0), 0.5 * (lo + hi))
        x = np.clip(x, lo, hi)
        fx = f(x)
        neg = fx < 0
        lo, flo = np.where(neg, x, lo), np.where(neg, fx, flo)
        hi, fhi = np.where(neg, hi, x), np.where(neg, fhi, fx)
        fhi = np.where(neg & (side == -1), fhi / 2, fhi)
        flo = np.where(~neg & (side == 1), flo / 2, flo)
        side = np.where(neg, -1.0, 1.0)
        if not inside.any() or float(np.nanmax(np.abs(fx[inside]))) <= ftol:
            break
    x = np.where(below, bracket[0], np.where(above, bracket[1], x))
    return x, int(below.sum() + above.sum())


def _rates_at(sp: StageParams, spont: np.ndarray, rows: np.ndarray, cols: np.ndarray, deltas: np.ndarray,
              kernel: np.ndarray, inhA: np.ndarray, inhB: np.ndarray, reading: str = 'stage') -> np.ndarray:
    """Evoked PN rates of glomeruli `rows` at odors `cols` with their own ORN deltas set to `deltas` (rows x cols), at
    the given LN (odors x steps): each entry a PN reading delta x kernel(t) as its ORN delta over time, its rate in
    `reading`."""
    win, pre = windows(sp, reading)
    R_, C_ = deltas.shape
    s = spont[rows][:, None] * np.ones((1, C_))
    spont_pn = s * sp.pn_inhsc / (spont.sum() + sp.pn_inhadd)
    inh_t = sp.pn_inhsc / (sp.pn_inhadd + 0.25 * inhA[cols] + 0.75 * inhB[cols])        # cols x steps
    pn = s.copy()
    acc_w = pn * win[0]
    acc_p = pn * pre[0]
    inh = np.zeros((1, C_))
    for t in range(1, len(kernel)):
        drive = deltas * kernel[t - 1]
        dd = -pn + spont_pn + 200.0 * np.tanh((drive + sp.pn_offset) * sp.pn_tanhsc / 200.0 * inh)
        inh = inh_t[:, t][None, :]
        pn = pn + dd * sp.dt / sp.pn_taum
        pn = np.where(0.0 < pn, pn, 0.0)
        if win[t]:
            acc_w += pn
        if pre[t]:
            acc_p += pn
    return acc_w / win.sum() - acc_p / pre.sum()


def solve_identity(stage: PNStage, glomeruli: Optional[Sequence[str]] = None, tol: float = SOLVE_TOL_HZ,
                   reading: str = 'preset') -> dict:
    """The solve's identity check: the stage's own rates as the target give back the stage's input (where the rate
    depends on the delta) and the same rates."""
    gl = list(stage.glomeruli if glomeruli is None else glomeruli)
    own = stage.run(reading=reading)
    solved, rec = solve_input(stage, own.loc[gl], gl, tol=tol / 10, reading=reading)
    d0 = stage.delta.loc[gl].to_numpy(dtype=float)
    d1 = solved.loc[gl].to_numpy(dtype=float)
    back = PNStage(solved, pd.Series(stage.spont, index=stage.glomeruli), stage.sp).run(reading=reading)
    rate_diff = float(np.abs(back.loc[gl].to_numpy() - own.loc[gl].to_numpy()).max())
    ok_rows = own.loc[gl].to_numpy() > 0
    res = dict(solve=rec, max_rate_diff=rate_diff, max_delta_diff=float(np.abs(d1 - d0)[ok_rows].max())
               if ok_rows.any() else 0.0, tol=tol)
    res['ok'] = bool(rate_diff <= tol and rec['converged'])
    return res


def route1_pn_rates(table: pd.DataFrame, set_glomeruli: Sequence[str],
                    set_odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Route 1, olfsysm's option for the PN stage's rates given directly (viki's, olfsysm 7528f6d / b1e89d7, Yang's OK
    of 2026-11-24 09:43): fit_mb_model's `pn_rates`, a glomeruli x odors table in the reading of `pn_evoked_rates`
    ('preset' here), NaN keeping the model's PN. Returns `table` on `set_glomeruli` and `set_odors` (default all its
    columns), NaN elsewhere, to pass as `extra_kws=dict(pn_rates=...)` to E83's `build_simulator` (or to
    fit_mb_model). The odor's LN inhibition stays the model's (route 2's solved input changes it)."""
    odors = list(table.columns if set_odors is None else set_odors)
    out = pd.DataFrame(np.nan, index=list(table.index), columns=list(table.columns))
    gl = [g for g in table.index if g in set(set_glomeruli)]
    out.loc[gl, odors] = table.loc[gl, odors].to_numpy(dtype=float)
    if not np.isfinite(out.to_numpy()).any():
        raise ValueError('no rate set')
    if np.isinf(out.to_numpy()).any():
        raise ValueError('infinite rates')
    return out


def route1_check(mp, rv, wPNKC: pd.DataFrame, odors: Sequence, table: pd.DataFrame, tol: float = 1e-6) -> dict:
    """A route-1 run's identity: its PN rates (pyolfsysm's `pn_evoked_rates` of the run's `mp`, `rv`) equal the set
    entries of `table` to `tol` Hz, and the unset entries are the model's own (checked by the caller against the run
    without `pn_rates`). Raises RuntimeError beyond `tol`."""
    from pyolfsysm.mb_model import pn_evoked_rates
    got = pn_evoked_rates(mp, rv, wPNKC=wPNKC, odors=odors)
    t = table.reindex(index=got.index, columns=got.columns)
    m = np.isfinite(t.to_numpy())
    d = float(np.abs(got.to_numpy()[m] - t.to_numpy()[m]).max()) if m.any() else 0.0
    if d > tol:
        raise RuntimeError(f'route 1: the run\'s PN rates differ from the set table by {d} Hz')
    return dict(n_set=int(m.sum()), max_abs_diff=d, tol=tol)


# ---------------------------------------------------------------------------
# step 2: the runs, and the score

def planned_step2(reverse: bool = False, with_rung: bool = True, with_dff_repeat: bool = False) -> List[Tuple[Arm, str]]:
    """Step 2's points (arm, wiring): 48 on the megamat input, 24 on the APL16 input, one repeat by the second route
    (73); without the rung's 24 if neither R2 nor R3 is the highest kept (`step2_rung`; R1 runs anyway); with
    `with_dff_repeat` the three imposed megamat arms again with the table set from the dendrites in dF/F
    (DFF_REPEAT_ARMS, 24 more), run only before outcome (e) or (f) is declared."""
    pts = [(Arm('megamat', a), w) for a in MEGAMAT_ARMS for w in WIRINGS]
    pts += [(Arm('apl16', a), w) for a in APL16_ARMS for w in WIRINGS]
    pts += [(REPEAT_ARM, WIRINGS[0])]
    if not with_rung:
        pts = [(a, w) for a, w in pts if a.arm not in RUNG_ARMS]
    if with_dff_repeat:
        pts += [(Arm('megamat', a), w) for a in DFF_REPEAT_ARMS for w in WIRINGS]
    return pts[::-1] if reverse else pts


def tag_of(arm: Arm, wiring: str) -> str:
    """'<panel>_<arm>_<wiring>', e.g. 'megamat_imposed_glom_perm_P1'."""
    return f'{arm.panel}_{arm.arm}_{wiring}'


def arm_seed(wiring: str, seed: int = SEED) -> int:
    """A permuted arm's seed on a wiring (ruling 6): one seeded permutation per wiring, `seed` plus the wiring's index
    in WIRINGS."""
    return int(seed) + WIRINGS.index(wiring)


def arm_units(arm: Arm) -> str:
    """The units an arm's table is set from: the dendrites in dF/F for the repeat's arms (DFF_REPEAT_ARMS), otherwise
    the lock's trial-noise units (FIT_UNITS)."""
    return 'dff' if arm.arm.endswith(DFF_SUFFIX) else FIT_UNITS


def claim(tag: str, lane: str, reclaim: bool = False, root: Optional[Path] = None) -> bool:
    """Claims a point by mkdir of results/E135_antennal_lobe_dendrites/claims/<tag>/ (atomic on NFS), as E134's runner
    does (its `claim`, rooted here)."""
    from experiments import E134_population_baseline as E134R
    return E134R.claim(tag, lane, reclaim=reclaim, root=out_dir() if root is None else root)


def release(tag: str, root: Optional[Path] = None) -> None:
    from experiments import E134_population_baseline as E134R
    E134R.release(tag, root=out_dir() if root is None else root)


def dendrite_means(units: str = FIT_UNITS, glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The dendrites the imposed arms are set from: the fly mean of all PN flies (glomeruli x the megamat odors), in
    `units` (the lock's trial-noise units, or dF/F for the repeat), on `glomeruli` (default the shared 28)."""
    gl = shared_glomeruli() if glomeruli is None else list(glomeruli)
    return pd.DataFrame(_nanmean(fly_stack(recorded_stage(PN_KEY, units), gl, MEGAMAT_ODORS)), index=gl,
                        columns=list(MEGAMAT_ODORS))


def arm_table(arm: Arm, stage: PNStage, dendrites: Optional[pd.DataFrame], glomeruli: Sequence[str],
              fits: Optional[Mapping[str, dict]] = None, rung: Optional[str] = None, seed: int = SEED,
              out_of_bounds: Optional[str] = CLIP_POLICY,
              bounds: Optional[Tuple[pd.DataFrame, pd.DataFrame]] = None) -> Tuple[pd.DataFrame, List[str], dict]:
    """An arm's table of PN rates on a stage's input ('preset' reading, route 1's), the glomeruli it sets, and its
    record.

    'R1': the R1 fits' stage (`rung_table`, at each fit's scale factor); 'rung' / 'rung_glom_perm': the highest kept
    rung (its terms' glomerulus labels permuted, `permute_terms`, for the control); 'imposed' and its permutations (and
    their dF/F repeats, the same with `dendrites` in dF/F): the dendrites (fly mean of all PN flies, glomeruli x odors)
    on the shared glomeruli, glomerulus or odor labels permuted for the controls, rank for rank on the stage's own rates
    (`imposed_rates`). `seed` is the permutation's (`arm_seed` of the wiring). Every table is then held inside the
    stage's reachable rates (`clip_table` with `out_of_bounds`, ruling 7; `bounds` default `stage.rate_bounds()`), and
    the record carries the clip's counts, the share clipped and whether the arm is named, for the rungs the scale
    factors (`scale_record`), and the odors the table sets ('set_odors': every odor for R1 and the rungs, the
    dendrites' odors the stage has for the imposed arms; on the APL16 input paa, euc, -aPine and B-myr keep the
    model's PNs there)."""
    rng = np.random.default_rng(seed)
    gl = [g for g in glomeruli if g in stage.glomeruli]
    rec: dict = dict(arm=arm.arm, panel=arm.panel, seed=int(seed), units=arm_units(arm))
    bounds = stage.rate_bounds() if bounds is None else bounds
    base = arm.arm[:-len(DFF_SUFFIX)] if arm.arm.endswith(DFF_SUFFIX) else arm.arm
    if base in ('R1', 'rung', 'rung_glom_perm'):
        r = 'R1' if base == 'R1' else rung
        if base == 'rung_glom_perm':
            fits = {k: dict(f, **{rung: permute_terms(Terms.from_dict(f[rung]), gl, np.random.default_rng(seed)
                                                      ).as_dict()}) for k, f in fits.items()}
        tab = rung_table(stage, fits, r)
        rec['rung'] = r
        rec['scale'] = scale_record(fits, (r,))
        tab, crec = clip_table(tab, bounds, list(stage.glomeruli), None, out_of_bounds)
        rec['clip'] = crec
        rec['set_odors'] = list(tab.columns)
        return tab, list(stage.glomeruli), rec
    own = stage.run(reading='preset')
    den = dendrites.loc[gl, [o for o in dendrites.columns if o in stage.odors]].astype(float)
    if base == 'imposed_glom_perm':
        den = permute_glomerulus_labels(den, gl, rng)
    elif base == 'imposed_odor_perm':
        cols = list(den.columns)
        perm = [cols[i] for i in rng.permutation(len(cols))]
        den = pd.DataFrame(den[perm].to_numpy(), index=den.index, columns=cols)
    elif base not in ('imposed', 'imposed_route2_repeat'):
        raise ValueError(arm)
    tab, r = imposed_rates(own, den, gl, bounds=bounds, out_of_bounds=out_of_bounds)
    rec['n_mapped'] = r['n_mapped']
    rec['clip'] = {k: v for k, v in r.items() if k != 'n_mapped'}
    rec['set_odors'] = list(den.columns)
    return tab, gl, rec


def arm_input(stage: PNStage, table: pd.DataFrame, set_glomeruli: Sequence[str], route: int,
              set_odors: Optional[Sequence[str]] = None) -> Tuple[dict, dict]:
    """How an arm's (clipped) table enters the family's simulator: route 1, extra_kws=dict(pn_rates=<the set entries>)
    on the panel's input (`route1_pn_rates`; the table's columns are the stage's odor labels, `input_columns` names
    them as the input does); route 2, the check of the option (ruling 9), the input solved on the stage
    (`solve_input`) for the glomeruli with input only: the set glomeruli at the table's rates, the other input glomeruli
    at their own rates there, the LN on the solved input. The glomeruli without input are not solved; their largest
    drift and the LN's change are in the record. A solve that does not converge is reported in the record
    ('converged' False), not forced. `set_odors` (default every column) are the odors the table sets (`arm_table`'s
    record); the others keep the model's PNs. Returns (build_simulator's keyword arguments, the record)."""
    if route == 1:
        return dict(extra_kws=dict(pn_rates=route1_pn_rates(table, set_glomeruli, set_odors))), dict(route=1)
    if route != 2:
        raise ValueError(f'{route=}')
    own = stage.run(reading='preset')
    has_input = [g for g in stage.glomeruli if np.any(stage.delta.loc[g].to_numpy() != 0)]
    setg = [g for g in set_glomeruli if g in set(has_input)]
    od = [o for o in (own.columns if set_odors is None else set_odors) if o in set(own.columns)]
    target = own.loc[has_input].copy()
    target.loc[setg, od] = table.loc[setg, od].to_numpy()
    solved, rec = solve_input(stage, target, has_input)
    rec = dict(rec, route=2, n_solved=len(has_input), n_without_input=len(stage.glomeruli) - len(has_input),
               no_input_max_drift_hz=rec['others_max_change'],
               set_without_input=[g for g in set_glomeruli if g not in set(has_input)])
    return dict(orn=solved), rec


def predicted_stage_rates(stage: PNStage, table: pd.DataFrame, set_glomeruli: Sequence[str],
                          bracket: Tuple[float, float] = (-1e4, 1e4)) -> pd.DataFrame:
    """What a route-1 run's PN rates should read in E134's reading ('stage': MODEL_WINDOW minus the settled pre-odor
    level), from the port: olfsysm's preset rates keep the LN and shift a PN's drive by a constant times the unit ORN
    kernel, so each set entry is the PN at the own delta whose 'preset' rate is the table's (`_solve_entries` at the
    stage's LN), read in 'stage'; the other entries are the stage's own. For `imposed_run_check` before or beside a
    run."""
    sp = stage.sp
    gl = [g for g in set_glomeruli if g in stage.glomeruli]
    rows = np.array([stage.glomeruli.index(g) for g in gl])
    cols = np.arange(len(stage.odors))
    k = orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
    inhA, inhB = stage.ln(sp)
    T = table.loc[gl, stage.odors].to_numpy(dtype=float)
    x, _ = _solve_entries(sp, stage.spont, rows, cols, T, k, inhA, inhB, 'preset', bracket, 200, 1e-8)
    out = stage.run(reading='stage')
    out.loc[gl] = _rates_at(sp, stage.spont, rows, cols, x, k, inhA, inhB, 'stage')
    return out


def imposed_run_check(table: pd.DataFrame, pn_rates_e134: pd.DataFrame, set_glomeruli: Sequence[str],
                      odors: Optional[Sequence[str]] = None, pn_matrix_error: Optional[float] = None,
                      flies_floor: Optional[float] = None) -> dict:
    """An imposed run's check (ruling 8): the largest difference between its table (set in `pn_evoked_rates`' reading)
    and its PN rates in E134's reading (`pn_rates_e134`, glomeruli x odors, from the run's saved PN stage or
    `predicted_stage_rates`) over the set entries, named above IMPOSED_CHECK_HZ, beside its PN-dendrite matrix error
    (which should sit near the flies' floor, about PN_MATRIX_FLOOR_NOTE)."""
    gl = [g for g in set_glomeruli if g in table.index and g in pn_rates_e134.index]
    od = [o for o in (table.columns if odors is None else odors) if o in pn_rates_e134.columns]
    d = np.abs(table.loc[gl, od].to_numpy(dtype=float) - pn_rates_e134.loc[gl, od].to_numpy(dtype=float))
    i, j = np.unravel_index(int(np.nanargmax(d)), d.shape) if np.isfinite(d).any() else (0, 0)
    mx = float(np.nanmax(d)) if np.isfinite(d).any() else np.nan
    return dict(max_abs_diff_hz=mx, at=dict(glomerulus=gl[i], odor=od[j]) if np.isfinite(mx) else None,
                mean_abs_diff_hz=float(np.nanmean(d)), named=bool(mx > IMPOSED_CHECK_HZ), threshold_hz=IMPOSED_CHECK_HZ,
                pn_matrix_error=None if pn_matrix_error is None else float(pn_matrix_error),
                flies_floor=PN_MATRIX_FLOOR_NOTE if flies_floor is None else float(flies_floor))


def routes_agree(errors_route1: Mapping[str, float], errors_route2: Mapping[str, float],
                 thresholds: Mapping[str, Mapping[str, float]], stages: Sequence[str] = CLASS_STAGES) -> dict:
    """Whether the route-2 repeat agrees with its route-1 run (ruling 9): each megamat stage's population error
    (boutons, claws, somata) differs between the routes by less than that stage's reset threshold (`class_thresholds`).
    No class rests on route 2; a disagreement is reported."""
    rows = {}
    for st in stages:
        a, b = errors_route1.get(st), errors_route2.get(st)
        thr = float(thresholds[st]['threshold']) if st in thresholds else np.nan
        diff = abs(float(a) - float(b)) if a is not None and b is not None else np.nan
        rows[st] = dict(route1=a, route2=b, diff=diff, threshold=thr, agree=bool(diff < thr))
    return dict(stages=rows, agree=bool(all(r['agree'] for r in rows.values())))


STEP2_PART = 'step2'
STEP1_FITS_FILE = 'step1_fits.json'
SCORE_VERDICT_FILE = 'score_verdict.json'
FAMILY = 'C2'
# the response rates a step-2 point is re-held at (the entry's "both response rates re-held", E77's `hold_rates`: the
# threshold with APL off, then the tied APL scales with APL on; the APL->bouton and KC-KC scales stay the fit's)
HOLD_RULES: Dict[str, str] = {
    'megamat': "E134's megamat held runs' targets (`held_targets`: the megamat flies' rate by Remy's criterion, APL off "
               "at the wiring's own fit's APL-off / APL-on ratio), so that the classes compare like with like",
    'apl16': "the wiring's own fit's response counts (E83 saved_fit target_on / target_off), the operating point its "
             "bands and APL effect are read at",
}
# a route-1 run's PN rates against its table (olfsysm sets them to within 1e-9 Hz) and its unset entries against the
# unchanged stage (the port, identical to olfsysm to about 1e-11 Hz)
ROUTE1_TOL_HZ = 1e-6


def step2_dirs() -> List[Path]:
    """Where step 2's runs are written and read: outputs/E135_antennal_lobe_dendrites/step2/ (local), then the
    model-output archive's data/model_outputs/E135_antennal_lobe_dendrites/step2/ (stage_readings' step-3 layout,
    `save_step3_run`, so that the score on any host reads every lane's runs)."""
    from kcq import stage_readings as R
    return R._archive_dirs(EXPERIMENT, STEP2_PART)


def read_step1_fits(root: Optional[Path] = None) -> dict:
    """Step 1's fits as step 2 takes them (results/.../step1_fits.json: 'all_flies', `rung_fits_all_flies`' fits on all
    seven PN flies with their scale factors; 'highest_kept'; 'step2_rung'). PermissionError until step 1 has written
    them."""
    p = (out_dir() if root is None else Path(root)) / STEP1_FITS_FILE
    d = _read_json(p)
    if not d.get('all_flies') or 'highest_kept' not in d:
        raise PermissionError(f"{p}: step 1's results are not written yet (--part step1); step 2 runs R1 and the kept "
                              'rung from its fits on all seven PN flies')
    return d


def require_step2() -> Tuple[dict, dict]:
    """Step 2's guards, in order: E134's measure gate reads "usable" and its megamat runs are scored
    (`require_e134_gate`), this experiment's committed lock (`read_lock`), step 1's fits (`read_step1_fits`). Each
    refuses with PermissionError. Returns (the lock, step 1's fits record)."""
    require_e134_gate()
    lock = read_lock()
    return lock, read_step1_fits()


def dff_repeat_asked(root: Optional[Path] = None) -> bool:
    """Whether the score has asked for the dF/F repeat (ruling 1's guard): its verdict's outcome is 'needs the dF/F
    repeat' (outcome (e) or (f) pending, results/.../score_verdict.json)."""
    v = _read_json((out_dir() if root is None else Path(root)) / SCORE_VERDICT_FILE)
    return (v.get('outcome') or {}).get('outcome') == 'needs the dF/F repeat'


def step2_rung_of(fits_rec: Mapping) -> Optional[str]:
    """The rung step 2 runs ('R2', 'R3' or None), as step 1 recorded it."""
    return fits_rec.get('step2_rung', step2_rung(fits_rec.get('highest_kept')))


def step2_points(fits_rec: Mapping, reverse: bool = False, dff: Optional[bool] = None) -> List[Tuple[Arm, str]]:
    """Step 2's points for this lane: `planned_step2` with the rung's 24 only if step 1 kept a rung above R1, and the
    dF/F repeat's 24 only when the score asks for them (`dff_repeat_asked`, unless `dff` is given)."""
    return planned_step2(reverse, with_rung=step2_rung_of(fits_rec) is not None,
                         with_dff_repeat=dff_repeat_asked() if dff is None else bool(dff))


def input_columns(table: pd.DataFrame, orn: pd.DataFrame) -> pd.DataFrame:
    """A table on the stage's odor labels (`odor_labels`) with the input's own columns, as fit_mb_model's `pn_rates`
    takes them (the odors of orn_deltas, with or without its 'panel' level, not their abbreviations); the input's odors
    the table lacks are NaN (the model's PNs), the table's odors the input lacks are dropped."""
    t = table.reindex(columns=odor_labels(orn.columns))
    return pd.DataFrame(t.to_numpy(dtype=float), index=t.index, columns=orn.columns)


def input_layout(delta: pd.DataFrame, orn: pd.DataFrame) -> pd.DataFrame:
    """A solved input on the stage's rows and odor labels (route 2, `arm_input`) in the panel input's layout (its
    glomeruli and columns), for build_simulator's `orn`. RuntimeError if a glomerulus outside the input was given a
    nonzero delta (the glomeruli without input are not solved)."""
    labels = odor_labels(orn.columns)
    rows = [str(g) for g in orn.index]
    extra = [g for g in delta.index if g not in set(rows)]
    if extra and np.any(delta.loc[extra, labels].to_numpy(dtype=float) != 0):
        raise RuntimeError(f'the solved input gives the glomeruli outside the input a delta: {extra}')
    x = delta.reindex(index=rows, columns=labels)
    if x.isna().any().any():
        raise RuntimeError("the solved input lacks some of the input's glomeruli or odors")
    return pd.DataFrame(x.to_numpy(dtype=float), index=orn.index, columns=orn.columns)


def family_parts(wiring: str) -> dict:
    """A wiring's fit as step 2 holds it (E83's `wiring_parts`, one simulation of the bouton base per process): the
    fit's parameters 'p', the lab input's 'breadth', the saved outputs, the number of KCs, the fit's response counts
    ('target_on', 'target_off', on 'size' KCs x odors of the APL16 input), its folder and tag."""
    from experiments import E83_cross_wiring_fit as E83
    Pp = E83.wiring_parts(wiring)
    return dict(p=dict(Pp['p']), breadth=Pp['breadth'], saved=Pp['saved'], n_kc=int(Pp['saved']['on'][0].shape[0]),
                target_on=int(Pp['target_on']), target_off=int(Pp['target_off']), size=int(Pp['size']),
                folder=Pp['folder'], tag=Pp['tag'])


def family_saved(wiring: str) -> Dict[str, tuple]:
    """A wiring's saved fit on the APL16 input ({'on', 'off'}: (responses, spike counts); E83's `saved_fit`, no
    simulation)."""
    from experiments import E83_cross_wiring_fit as E83
    return E83.saved_fit(wiring)['saved']


def step2_targets(panel: str, wiring: str, n_odors: int, fam: Mapping, lock134: Optional[Mapping] = None) -> dict:
    """A point's response counts to re-hold (HOLD_RULES): on the megamat input E134's `held_targets` of the same
    wiring's held point, on its KCs x `n_odors`; on the APL16 input the wiring's own fit's counts."""
    size = int(fam['n_kc']) * int(n_odors)
    if panel == 'megamat':
        from experiments import E134_population_baseline as E134R
        t_off, t_on, f_off, f_on = E134R.held_targets(E134R.Point('megamat', 'held', wiring), size, lock134)
    elif panel == 'apl16':
        if size != int(fam['size']):
            raise RuntimeError(f'{wiring}: {size} KCs x odors on the APL16 input, the fit has {fam["size"]}')
        t_off, t_on = int(fam['target_off']), int(fam['target_on'])
        f_off, f_on = t_off / size, t_on / size
    else:
        raise ValueError(f'{panel=}')
    return dict(target_off=int(t_off), target_on=int(t_on), frac_off_target=float(f_off), frac_on_target=float(f_on),
                size=size, hold_rule=HOLD_RULES[panel])


def family_identity(panel: str, wiring: str, sim: Callable) -> dict:
    """The unchanged family on a panel's input at its fit's parameters, APL on and acutely off (APL->KC and
    APL->bouton scales 0), on `family_stage`'s simulator, reproduced exactly: on megamat E134's check of an as-fitted
    point (E107's saved run and the saved fit on the odors the panels share, `E134_population_baseline.check_identity`),
    on APL16 the saved fit itself (every count and response). The record's 'exact' says whether it holds."""
    from experiments import E134_population_baseline as E134R
    from kcq import stage_readings as R
    fam = family_parts(wiring)
    p = dict(fam['p'])
    out = {}
    for c, q in (('on', p), ('off', dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))):
        r, s = sim(dict(q))
        out[c] = dict(responses=r, spike_counts=s)
    if panel == 'megamat':
        rec = E134R.check_identity(E134R.Point('megamat', 'fitted', wiring), out)
    elif panel == 'apl16':
        saved = fam['saved']
        rec = dict(tag=f'apl16_fitted_{wiring}', checks=[])
        for c in ('on', 'off'):
            rec['checks'].append(dict(name=f'saved fit {c} counts', differ=E134R._differ(
                R._abbrev(out[c]['spike_counts']), R._abbrev(saved[c][1]))))
            rec['checks'].append(dict(name=f'saved fit {c} responses', differ=E134R._differ(
                R._abbrev(out[c]['responses']), R._abbrev(saved[c][0]))))
        rec['exact'] = all(x['differ'] == 0 for x in rec['checks'])
    else:
        raise ValueError(f'{panel=}')
    rec.update(panel=panel, wiring=wiring)
    return rec


def identity_path(lane: str) -> Path:
    """results/.../identity_2_<host>_<lane>.json: the lane's identity of the unchanged family per wiring ('family:
    <wiring>') and each point's hold identity (by tag)."""
    return out_dir() / f'identity_2_{_host()}_{lane}.json'


def _identity_put(lane: str, key: str, rec: Mapping) -> None:
    path = identity_path(lane)
    ids = _read_json(path)
    ids[key] = dict(rec)
    _write_json(path, ids)


def lane_identity(panel: str, wiring: str, lane: str, finfo: Mapping) -> dict:
    """The identity check of a lane's first point on each wiring: the unchanged family reproduced exactly
    (`family_identity` on `family_stage`'s simulator), once per host, lane and wiring (cached in `identity_path`; the
    panel is that of the lane's first point on the wiring). RuntimeError if it does not hold."""
    key = f'family:{wiring}'
    old = _read_json(identity_path(lane)).get(key)
    if old and old.get('exact'):
        return dict(old, cached=True)
    rec = family_identity(panel, wiring, finfo['sim'])
    rec.update(host=_host(), lane=lane, time=_now(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _identity_put(lane, key, rec)
    if not rec.get('exact'):
        raise RuntimeError(f'{wiring} on {panel}: the unchanged family does not reproduce its saved run on lane {lane}: '
                           f'{rec}')
    return rec


class Step2:
    """What step 2's points share in one process: the lock, step 1's fits on all seven PN flies and its rung, E134's
    lock, the family's PN stage per (panel, wiring) (`family_stage`, built once each; the lane's identity on the
    wiring is checked on that build, `lane_identity`, and the build is then freed), and the dendrite means per units
    (`dendrite_means`, on the lock's glomeruli)."""

    def __init__(self, lock: Mapping, fits_rec: Mapping, lane: str, lock134: Optional[Mapping] = None):
        self.lock = lock
        self.fits_rec = fits_rec
        self.fits = fits_rec['all_flies']
        self.rung = step2_rung_of(fits_rec)
        self.glomeruli = [str(g) for g in lock['glomeruli']]
        self.lane = lane
        self.lock134 = lock134
        self._stages: Dict[Tuple[str, str], Tuple[PNStage, dict]] = {}
        self._den: Dict[str, pd.DataFrame] = {}

    def stage(self, panel: str, wiring: str) -> Tuple[PNStage, dict]:
        key = (panel, wiring)
        if key not in self._stages:
            stage, finfo = family_stage(panel, wiring)
            ident = lane_identity(panel, wiring, self.lane, finfo)
            port = finfo.get('identity')
            finfo = None
            gc.collect()
            self._stages[key] = (stage, dict(port=port, lane_identity=ident))
        return self._stages[key]

    def dendrites(self, units: str) -> pd.DataFrame:
        if units not in self._den:
            self._den[units] = dendrite_means(units, self.glomeruli)
        return self._den[units]


def _hold(wiring: str, orn: pd.DataFrame, fam: Mapping, extra: Optional[Mapping], targets: Mapping,
          label: str) -> dict:
    """Both response rates re-held (E77's `hold_rates`) on the family's simulator of `wiring` on `orn` with `extra`
    (route 1's pn_rates), the breadth held at the lab input's, from the fit's parameters."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import hold_rates
    p = dict(fam['p'])
    return hold_rates(lambda: E83.build_simulator(wiring, orn=orn, breadth=fam['breadth'], params=dict(p),
                                                  extra_kws=None if not extra else dict(extra)),
                      target_off=targets['target_off'], target_on=targets['target_on'], params=p, label=label)


def _build(panel: str, wiring: str, params: Mapping[str, float], orn: Optional[pd.DataFrame],
           extra: Optional[Mapping]):
    """The held point built again at the held parameters (E134's `build_point` with its overrides, d2c7c9d)."""
    from experiments import E134_population_baseline as E134R
    return E134R.build_point(E134R.Point(panel, 'held', wiring), dict(params), orn=orn, extra_kws=extra)


def _run(panel: str, wiring: str, built, params: Mapping[str, float], params_off: Mapping[str, float]) -> dict:
    """Both conditions with every stage saved (E134's `run_stages`)."""
    from experiments import E134_population_baseline as E134R
    return E134R.run_stages(E134R.Point(panel, 'held', wiring), built, params, params_off)


def _hold_identity(panel: str, wiring: str, out: Mapping, hold: Mapping) -> dict:
    """The run at the held parameters reproduces the hold's final outputs exactly (E134's `check_identity`)."""
    from experiments import E134_population_baseline as E134R
    return E134R.check_identity(E134R.Point(panel, 'held', wiring), out, {c: hold[c] for c in ('on', 'off')})


def _provenance(orn: pd.DataFrame) -> dict:
    return H.model_provenance(orn_deltas=orn)


def olfsysm_odor_names(built) -> List[str]:
    """A built point's odor labels (`odor_labels` of its input) in olfsysm's column order (its simulator's 'odor_order':
    output column j is olfsysm's odor order[j])."""
    _, _, order = built.vars()
    labels = odor_labels(built.orn.columns)
    names: List[Optional[str]] = [None] * len(order)
    for j, i in enumerate(order):
        names[int(i)] = labels[j]
    return [str(x) for x in names]


def run_pn_check(built, stage: PNStage, table: pd.DataFrame, set_glomeruli: Sequence[str],
                 set_odors: Optional[Sequence[str]] = None, tol: float = ROUTE1_TOL_HZ, strict: bool = True) -> dict:
    """A run's PN rates as olfsysm's option reads them (pyolfsysm's `pn_evoked_rates` of the run's variables, 'preset')
    against its table on the set entries, and against the unchanged stage's own rates (the port, `PNStage.run`) on the
    others. Route 1: both to `tol`, RuntimeError beyond it with `strict`; route 2 (the solved input): reported only
    (the set entries to the solve's tolerance, the others as the LN on the solved input leaves them)."""
    from pyolfsysm.mb_model import pn_evoked_rates
    mp, rv, _ = built.vars()
    got = pn_evoked_rates(mp, rv, wPNKC=built.W, odors=olfsysm_odor_names(built))
    got.index = [str(g) for g in got.index]
    gl = [g for g in set_glomeruli if g in set(got.index)]
    od = [o for o in (table.columns if set_odors is None else set_odors) if o in set(got.columns)]
    d_set = np.abs(got.loc[gl, od].to_numpy(dtype=float) - table.loc[gl, od].to_numpy(dtype=float))
    own = stage.run(reading='preset').reindex(index=got.index, columns=got.columns)
    unset = np.ones(got.shape, dtype=bool)
    unset[np.ix_([list(got.index).index(g) for g in gl], [list(got.columns).index(o) for o in od])] = False
    d_un = np.abs(got.to_numpy(dtype=float) - own.to_numpy(dtype=float))[unset]
    rec = dict(n_set=int(np.isfinite(d_set).sum()), max_abs_diff_set=float(np.nanmax(d_set)) if d_set.size else 0.0,
               n_unset=int(unset.sum()), max_abs_diff_unset=float(np.nanmax(d_un)) if d_un.size else 0.0, tol=float(tol),
               strict=bool(strict))
    rec['ok'] = bool(rec['max_abs_diff_set'] <= tol and rec['max_abs_diff_unset'] <= tol)
    if strict and not rec['ok']:
        raise RuntimeError(f"the run's PN rates are not its table's (set) or the model's own (unset): {rec}")
    return rec


def e134_pn_rates(cond: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """A condition's PN rates in E134's reading (stage_readings': the MODEL_WINDOW mean minus the settled pre-odor
    mean), from `run_stages`' output ('pn_pre_settled') or a loaded run (whose 'pn_pre' is the settled level)."""
    pre = cond['pn_pre_settled'] if 'pn_pre_settled' in cond else cond['pn_pre']
    x = cond['pn_win'] - pre
    x.index = [str(g) for g in x.index]
    return x


def stage_rates(conds: Mapping[str, Mapping[str, pd.DataFrame]]) -> dict:
    """A run's largest PN and bouton rates (window means, Hz) per condition and the share of bouton x odor window
    means above BOUTON_RATE_CAP_HZ (E127's measure), the larger of the two conditions' beside; 'bouton_cap_excluded'
    where it is above BOUTON_RATE_CAP_SHARE (no class)."""
    rec = {}
    for c in ('on', 'off'):
        st = conds.get(c, {})
        if 'pn_win' in st:
            rec[f'max_pn_rate_{c}'] = float(np.nanmax(st['pn_win'].to_numpy(dtype=float)))
        if 'bouton_win' in st:
            b = st['bouton_win'].to_numpy(dtype=float)
            rec[f'max_bouton_rate_{c}'] = float(np.nanmax(b))
            rec[f'bouton_share_above_300_{c}'] = float(np.mean(b > BOUTON_RATE_CAP_HZ))
    sh = [v for k, v in rec.items() if k.startswith('bouton_share_above_300_')]
    rec['bouton_share_above_300'] = float(max(sh)) if sh else np.nan
    rec['bouton_cap_excluded'] = bool(sh and max(sh) > BOUTON_RATE_CAP_SHARE)
    return rec


def save_point(tag: str, out: Mapping[str, dict], built, info: Mapping[str, object], table: pd.DataFrame,
               dirs: Optional[Sequence[Path]] = None) -> List[Path]:
    """A point's table of PN rates (<tag>_pn_table.parquet, the stage's labels, 'preset' reading) and its run in
    stage_readings' step-3 layout (`save_step3_run`: every stage matrix, responses, spike counts, APL means, the units'
    metadata, the input, and <tag>_run.json last), locally and in the archive (`step2_dirs`)."""
    from kcq import stage_readings as R
    dirs = [Path(d) for d in (step2_dirs() if dirs is None else dirs)]
    t = table.copy()
    t.columns = [str(c) for c in t.columns]
    t.index = [str(g) for g in t.index]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / f'.{tag}_pn_table.parquet.{os.getpid()}.tmp'
        t.to_parquet(tmp)
        os.replace(tmp, d / f'{tag}_pn_table.parquet')
    return R.save_step3_run(tag, out, built.meta, built.orn, info, dirs)


def runs2_table(root: Optional[Path] = None) -> pd.DataFrame:
    """Every row of every host's and lane's runs_2_*.csv, the last row per tag."""
    d = out_dir() if root is None else Path(root)
    fs = [pd.read_csv(p) for p in sorted(d.glob('runs_2_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


# what each step-2 arm is, for the archive's params.csv (`write_variants`)
ARM_KINDS: Dict[str, str] = {
    'R1': "R1, the global parameters fitted to the PN dendrites (step 1), its table of PN rates set by route 1",
    'rung': "the highest kept rung of step 1, its table of PN rates set by route 1",
    'rung_glom_perm': "the highest kept rung with its terms' glomerulus labels permuted, set by route 1",
    'imposed': "the recorded PN dendrites (fly mean, the lock's noise units) imposed rank for rank on the model's rates, "
               "route 1",
    'imposed_glom_perm': "the recorded PN dendrites imposed, glomerulus labels permuted, route 1",
    'imposed_odor_perm': "the recorded PN dendrites imposed, odor labels permuted, route 1",
    'imposed_route2_repeat': "the recorded PN dendrites imposed by route 2 (an ORN input solved for the table)",
}
VARIANT_COLUMNS: Tuple[str, ...] = (
    'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale', 'off_fixed_thr',
    'off_wAPLKC_scale', 'off_wKCAPL_scale', 'off_wAPLPN_scale', 'off_wPNAPL_scale', 'off_wKCKC_scale', 'thr_factor',
    'apl_factor', 'target_on', 'target_off', 'frac_on_target', 'frac_off_target', 'at_bound_thr', 'at_bound_apl',
    'fit_fixed_thr', 'n_kcs', 'n_odors', 'rate_on', 'rate_off', 'rate_ratio', 'max_count_on', 'max_count_off',
    'n_set_glomeruli', 'n_set_odors', 'clip_n_clipped', 'table_sha256', 'input_sha256', 'stages', 'host', 'lane', 'time',
    'queue_commit', 'pyolfsysm_commit')


def _point_built(tag: str, root: Optional[Path] = None) -> dict:
    """A step-2 run's model build as its lane recorded it (provenance_2_*.json, points[tag]['built_info']: model,
    wiring, fit_folder, fit_tag, kws, breadth); {} if no lane recorded it."""
    d = out_dir() if root is None else Path(root)
    for f in sorted(d.glob('provenance_2_*.json')):
        b = _read_json(f).get('points', {}).get(tag, {}).get('built_info')
        if b:
            return b
    return {}


def write_variants(root: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv: one row per saved step-2 run (`runs2_table`; the archive's step2/<tag>_<on|acute_off>_...), for the
    model-output archive's params.csv: the tag (the files' variant name), panel, arm, wiring, route, units, rung and
    seed, what the arm is (ARM_KINDS; a dF/F repeat's table set from the dendrites in dF/F), the absolute parameters of
    both conditions (VARIANT_COLUMNS), and the family's build from the lane's provenance (kws, fit_folder, fit_tag,
    breadth)."""
    d = out_dir() if root is None else Path(root)
    runs = runs2_table(d)
    rows = []
    for r in runs.to_dict('records'):
        tag = str(r['tag'])
        arm = str(r.get('arm'))
        base = arm[:-len(DFF_SUFFIX)] if arm.endswith(DFF_SUFFIX) else arm
        kind = ARM_KINDS.get(base, base)
        if arm != base:
            kind += '; the table set from the dendrites in dF/F (ruling 1\'s repeat)'
        b = _point_built(tag, d)
        rows.append(dict(tag=tag, variant=tag, part=STEP2_PART, subdir=STEP2_PART, panel=r.get('panel'), arm=arm,
                         model=b.get('model', 'C2'), wiring=r.get('wiring'), route=r.get('route'),
                         units=r.get('units'), rung=r.get('rung'), seed=r.get('seed'),
                         kind=f"C2's family on wiring {r.get('wiring')}, on the {r.get('panel')} input, both response "
                              f"rates re-held: {kind}",
                         **{c: r.get(c) for c in VARIANT_COLUMNS if c in r},
                         kws=json.dumps(b.get('kws'), sort_keys=True) if b.get('kws') is not None else None,
                         fit_folder=b.get('fit_folder'), fit_tag=b.get('fit_tag'), breadth=b.get('breadth')))
    df = pd.DataFrame(rows)
    df.to_csv(d / 'variants.csv', index=False)
    print(f'variants.csv: {len(df)} runs -> {d / "variants.csv"}', flush=True)
    return df


def _base_arm(arm: Arm) -> str:
    return arm.arm[:-len(DFF_SUFFIX)] if arm.arm.endswith(DFF_SUFFIX) else arm.arm


def run_point2(arm: Arm, wiring: str, ctx: Step2, table_path: Path, prov_path: Path, reclaim: bool = False,
               done: Optional[set] = None) -> Optional[dict]:
    """One point of step 2, claimed (None if done or claimed elsewhere; the claim is released on failure).

    1. The family's PN stage on the point's input and wiring (`Step2.stage`; the lane's first point on the wiring checks
       that the unchanged family reproduces its saved run exactly, `lane_identity`).
    2. The arm's table of PN rates (`arm_table`: R1 or step 1's kept rung at their scale factors, each megamat odor
       from the fold fit that left it out; or the imposed dendrites, permuted for the controls with the wiring's
       seed, `arm_seed`; every table clipped inside the stage's reach, ruling 7).
    3. Route 1 (olfsysm's pn_rates on the set entries, `input_columns`) or, for the repeat, route 2 (the input solved,
       `arm_input`, put back in the input's layout, `input_layout`; ruling 9: reported, not forced).
    4. Both response rates re-held (`_hold`, `step2_targets`), the point built again at the held parameters and run
       with every stage (`_build`, `_run`), which must reproduce the hold exactly (`_hold_identity`).
    5. The run's PN rates against its table (`run_pn_check`: exact for route 1), an imposed run's check in E134's
       reading (ruling 8, `imposed_run_check`), the largest PN and bouton rates (`stage_rates`).
    6. Saved (`save_point`), its row appended to `table_path` and its provenance (`_provenance`, the table's, the
       dendrites', the locks' and step 1's sha256) to `prov_path`.
    """
    from experiments import E134_population_baseline as E134R
    from experiments.E77_input_screen import frame_sha256
    tag = tag_of(arm, wiring)
    if done is not None and tag in done:
        return None
    if not claim(tag, ctx.lane, reclaim):
        print(f'E135 {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        base = _base_arm(arm)
        if base in RUNG_ARMS and ctx.rung is None:
            raise RuntimeError(f'{tag}: step 1 kept no rung above R1, so the rung\'s points are not run')
        stage, srec = ctx.stage(arm.panel, wiring)
        den = None if base in ('R1',) + RUNG_ARMS else ctx.dendrites(arm_units(arm))
        table, set_gl, trec = arm_table(arm, stage, den, ctx.glomeruli, fits=ctx.fits, rung=ctx.rung,
                                        seed=arm_seed(wiring))
        set_od = list(trec['set_odors'])
        route = 2 if arm == REPEAT_ARM else 1
        kw, rrec = arm_input(stage, table, set_gl, route, set_odors=set_od)
        orn = panel_orn(arm.panel)
        if route == 1:
            extra, orn_in = dict(pn_rates=input_columns(kw['extra_kws']['pn_rates'], orn)), None
        else:
            extra, orn_in = None, input_layout(kw['orn'], orn)
        fam = family_parts(wiring)
        tg = step2_targets(arm.panel, wiring, orn.shape[1], fam, ctx.lock134)
        hold = _hold(wiring, orn if orn_in is None else orn_in, fam, extra, tg, f'E135 {tag}')
        params, params_off = dict(hold['params']), dict(hold['params_off'])
        built = _build(arm.panel, wiring, params, orn_in, extra)
        out = _run(arm.panel, wiring, built, params, params_off)
        ident = _hold_identity(arm.panel, wiring, out, hold)
        _identity_put(ctx.lane, tag, ident)
        if not ident['exact']:
            raise RuntimeError(f'{tag}: the run at the held parameters is not the hold\'s: {ident}')
        pn = run_pn_check(built, stage, table, set_gl, set_od, tol=ROUTE1_TOL_HZ if route == 1 else SOLVE_TOL_HZ,
                          strict=route == 1)
        imp = (imposed_run_check(table, e134_pn_rates(out['on']), set_gl, set_od)
               if base in IMPOSED_ARMS + (REPEAT_ARM.arm,) else None)
        sr = stage_rates(out)
        table_sha = frame_sha256(table)
        info = dict(model=FAMILY, wiring=wiring, panel=arm.panel, reading='held', tag=tag, experiment=EXPERIMENT,
                    arm=arm.arm, route=route, units=arm_units(arm), rung=trec.get('rung'), seed=trec['seed'],
                    set_glomeruli=list(set_gl), set_odors=set_od, table_sha256=table_sha, stages=list(built.stages),
                    params=params, params_off=params_off, targets=tg, host=_host(), lane=ctx.lane, time=_now())
        save_point(tag, out, built, info, table)
        clip = trec.get('clip', {})
        scale = (trec.get('scale') or {}).get('factors', {}).get(trec.get('rung'), {})
        row = dict(tag=tag, panel=arm.panel, arm=arm.arm, wiring=wiring, route=route, units=arm_units(arm),
                   rung=trec.get('rung'), seed=trec['seed'],
                   **{k: float(v) for k, v in params.items() if np.isscalar(v)},
                   **{f'off_{k}': float(v) for k, v in params_off.items() if np.isscalar(v)},
                   n_kcs=int(out['on']['spike_counts'].shape[0]), n_odors=int(out['on']['spike_counts'].shape[1]),
                   **E134R._counts_row(out),
                   thr_factor=hold['thr_factor'], apl_factor=hold['apl_factor'], target_on=tg['target_on'],
                   target_off=tg['target_off'], frac_on_target=tg['frac_on_target'],
                   frac_off_target=tg['frac_off_target'], at_bound_thr=hold['at_bound_thr'],
                   at_bound_apl=hold['at_bound_apl'], hold_n_sims=hold['n_sims'],
                   fit_fixed_thr=float(fam['p']['fixed_thr']),
                   n_set_glomeruli=len(set_gl), n_set_odors=len(set_od), n_mapped=trec.get('n_mapped'),
                   clip_n_set=clip.get('n_set'), clip_n_clipped=clip.get('n_clipped'),
                   clip_share=clip.get('share_clipped'), clip_largest_excursion_hz=clip.get('largest_excursion_hz'),
                   clip_named=clip.get('named'),
                   scale_min=min(scale.values()) if scale else np.nan, scale_max=max(scale.values()) if scale else np.nan,
                   scale_named=','.join((trec.get('scale') or {}).get('named', {}).get(trec.get('rung'), [])) or None,
                   pn_check_max_set_hz=pn['max_abs_diff_set'], pn_check_max_unset_hz=pn['max_abs_diff_unset'],
                   imposed_check_max_hz=None if imp is None else imp['max_abs_diff_hz'],
                   imposed_check_named=None if imp is None else imp['named'],
                   **({} if route == 1 else dict(
                       solve_converged=rrec['converged'], solve_max_abs_diff_hz=rrec['max_abs_diff'],
                       solve_at_edge=rrec['at_edge'], solve_method=rrec['method'],
                       solve_no_input_max_drift_hz=rrec['no_input_max_drift_hz'],
                       ln_inhA_max_abs=rrec['ln_change']['inhA_max_abs'],
                       ln_pn_inh_factor_max_rel=rrec['ln_change']['pn_inh_factor_max_rel'])),
                   **sr, table_sha256=table_sha, input_sha256=frame_sha256(built.orn), stages='+'.join(built.stages),
                   seconds=round(time.time() - t0, 1), host=_host(), lane=ctx.lane, time=_now(),
                   threads=os.environ.get('OMP_NUM_THREADS'), identity_exact=True,
                   lane_identity_exact=bool(srec['lane_identity'].get('exact')),
                   **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
        E134R._append_csv(table_path, row)
        prov = _read_json(prov_path) or dict(part='step2', lane=ctx.lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(
            info=info, built_info={k: v for k, v in built.info.items() if not str(k).startswith('_')}, arm=trec,
            route=rrec, pn_check=pn, imposed_check=imp, identity=ident, lane_identity=srec['lane_identity'],
            port_identity=srec['port'], targets=tg, hold=dict(thr_factor=hold['thr_factor'],
                                                              apl_factor=hold['apl_factor'], n_sims=hold['n_sims']),
            table_sha256=table_sha, dendrites_sha256=None if den is None else frame_sha256(den),
            lock_sha256=_sha256(LOCK_COPY) if LOCK_COPY.exists() else None,
            e134_lock_sha256=_sha256(LOCK_E134()) if LOCK_E134().exists() else None,
            step1_fits_sha256=_sha256(out_dir() / STEP1_FITS_FILE) if (out_dir() / STEP1_FITS_FILE).exists() else None,
            thresholds_sha256=_sha256(out_dir() / THRESHOLDS_FILE) if (out_dir() / THRESHOLDS_FILE).exists() else None,
            model=_provenance(built.orn))
        _write_json(prov_path, prov)
        print(f"E135 {tag}: {row['seconds']:.0f} s, rate on {row['rate_on']:.4f} / off {row['rate_off']:.4f}, "
              f"clipped {clip.get('n_clipped')} of {clip.get('n_set')}", flush=True)
        del built, out, hold
        return row
    except BaseException:
        release(tag)
        raise
    finally:
        gc.collect()


def part_step2(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 2, resumable, one point claimed at a time (`step2_points`, `run_point2`); every lane on any host runs the
    whole list and skips what another has claimed or done. Refuses before E134's gate, this experiment's lock and step
    1's fits (`require_step2`). The rung's 24 points run only if step 1 kept a rung above R1, the dF/F repeat's 24
    only once the score asks for them. Before the first arm the classes' thresholds are frozen (`freeze_thresholds`:
    written once; the file present at step 2's start decides). Rows in runs_2_<host>_<lane>.csv, provenance in
    provenance_2_<host>_<lane>.json, identities in identity_2_<host>_<lane>.json, the runs in `step2_dirs`."""
    lock, fits_rec = require_step2()
    from experiments import E134_population_baseline as E134R
    lock134 = require_e134_lock()
    E134R.check_wirings()
    if tuple(WIRINGS) != tuple(E134R.WIRINGS):
        raise RuntimeError(f"E135's wirings {WIRINGS} are not E134's {E134R.WIRINGS}")
    lane = lane or _host()
    thr = freeze_thresholds()
    ctx = Step2(lock, fits_rec, lane, lock134)
    pts = step2_points(fits_rec, reverse)
    table = out_dir() / f'runs_2_{_host()}_{lane}.csv'
    prov = out_dir() / f'provenance_2_{_host()}_{lane}.json'
    print(f'E135 step 2, lane {lane}: {len(pts)} points (rung {ctx.rung}, dF/F repeat '
          f'{any(a.arm in DFF_REPEAT_ARMS for a, _ in pts)}; thresholds frozen {thr.get("written")}, somata '
          f'{thr["classes"]["kc_soma"]["threshold"]:.4f} on strength by odor)', flush=True)
    for arm, w in pts:
        done = set(runs2_table()['tag'].astype(str))
        run_point2(arm, w, ctx, table, prov, reclaim, done)
    return runs2_table()


def read_e134_table(name: str) -> pd.DataFrame:
    """One of E134's score tables (`e134_results_dir`: score_stage_table.csv, step2_scores.csv, ...);
    FileNotFoundError until E134 has written it."""
    p = e134_results_dir() / name
    if not p.exists():
        raise FileNotFoundError(f"{p}: E134 has not written it yet")
    return pd.read_csv(p)


def e134_held_rows(stage_table: pd.DataFrame, group: str = E134_HELD_GROUP,
                   readings: Optional[Mapping[str, str]] = None) -> pd.DataFrame:
    """E134's runs held at the megamat flies' response rate, as the classes read them: the rows of its stage table of
    group `group`, the megamat panel, APL on, each stage's default reading (or the reading `readings` names for it) and
    the eight wirings (one row per stage and wiring)."""
    from kcq import stage_readings as R
    t = stage_table
    need = {'stage', 'wiring', 'error'}
    if not need <= set(t.columns):
        raise KeyError(f'the stage table needs columns {sorted(need)}')
    if 'group' in t.columns:
        t = t[t.group == group]
    if 'panel' in t.columns:
        t = t[t.panel == 'megamat']
    if 'condition' in t.columns:
        t = t[t.condition == 'on']
    if 'reading' in t.columns:
        want = {st: dict(readings or {}).get(st, R.DEFAULT_READINGS[st].name) for st in R.DEFAULT_READINGS}
        t = t[np.array([rd == want.get(st) for st, rd in zip(t.stage, t.reading)], dtype=bool)]
    return t[t.wiring.isin(WIRINGS)].drop_duplicates(['stage', 'wiring'])


def measure_values(t: pd.DataFrame, measure: str = 'all', matrix_for_dendrites: bool = False) -> np.ndarray:
    """Each row's error under a measure (MEASURES): 'all' its population error ('error'; with `matrix_for_dendrites`
    the PN dendrites' 'error_with_matrix', as E134's megamat stage error reads them), else the root mean square of the
    named views' errors ('view_<view>'; 'strength' is the adopted measure's view, 'strength_sharing' E134_archive_score's
    two-view error)."""
    if measure not in MEASURES:
        raise ValueError(f'measure {measure!r} not in {sorted(MEASURES)}')
    views = MEASURES[measure]
    if views is None:
        e = t['error'].to_numpy(dtype=float)
        if matrix_for_dendrites and 'error_with_matrix' in t.columns:
            e = np.where(t['stage'].to_numpy() == 'pn_dendrite', t['error_with_matrix'].to_numpy(dtype=float), e)
        return e
    missing = [f'view_{v}' for v in views if f'view_{v}' not in t.columns]
    if missing:
        raise KeyError(f'measure {measure!r} needs the columns {missing}')
    return np.sqrt(np.mean([t[f'view_{v}'].to_numpy(dtype=float) ** 2 for v in views], axis=0))


def _column_values(t: pd.DataFrame, column: str) -> np.ndarray:
    """A column of a scores table, 'error_strength_sharing' computed where absent (`measure_values`)."""
    if column in t.columns:
        return t[column].to_numpy(dtype=float)
    if column == 'error_strength_sharing':
        return measure_values(t, 'strength_sharing')
    raise KeyError(f'no column {column!r}')


# the columns reported beside the classes' measure for every arm, wiring and megamat stage (the rulings of 15:41): the
# three-view population error, the PN dendrites' with the matrix view, the two-view error and the views themselves
BESIDE_COLUMNS: Tuple[str, ...] = ('error', 'error_with_matrix', 'error_strength_sharing', 'view_strength',
                                   'view_sharing', 'view_breadth', 'view_matrix')


def e134_base(stage_table: Optional[pd.DataFrame] = None, group: str = E134_HELD_GROUP,
              readings: Optional[Mapping[str, str]] = None) -> pd.DataFrame:
    """The base of the classes, E134's own run of each wiring held at the megamat flies' response rate
    (`e134_held_rows`, in each stage's default reading or `readings`): per (stage, wiring) its strength-by-odor view's
    error ('strength_base', the classes' measure), its population error ('error_base', the guard's) and beside them
    '<column>_base' of BESIDE_COLUMNS."""
    t = e134_held_rows(read_e134_table(E134_STAGE_TABLE) if stage_table is None else stage_table, group, readings)
    out = t[['stage', 'wiring']].copy()
    for c in BESIDE_COLUMNS:
        name = dict(error='error_base', view_strength='strength_base').get(c, f'{c}_base')
        try:
            out[name] = _column_values(t, c)
        except KeyError:
            out[name] = np.nan
    return out.reset_index(drop=True)


def class_thresholds(stage_table: Optional[pd.DataFrame] = None, group: str = E134_HELD_GROUP,
                     stages: Sequence[str] = CLASS_STAGES, column: str = 'error',
                     readings: Optional[Mapping[str, str]] = None) -> dict:
    """Thresholds in one-fly units (ruling 10, and the rulings of 15:41): per megamat stage, THRESHOLD_SE_MULT
    standard errors over the eight wirings (the SD, ddof 1, of the wirings' values over sqrt(8)) of `column` in E134's
    runs held at the megamat flies' response rate (its stage table's group 'megamat_held', the megamat panel, APL on,
    each stage's default reading or `readings`). The classes read 'view_strength' (the somata's replaces the entry's
    0.15 and ±0.05 in `classify`, each stage's is the routes' agreement threshold); the guard reads 'view_strength' and
    'error' (`threshold_values`). Reads <shared>/results/E134_population_baseline/score_stage_table.csv by default."""
    t = e134_held_rows(read_e134_table(E134_STAGE_TABLE) if stage_table is None else stage_table, group, readings)
    out = {}
    for st, sub in t[t.stage.isin(stages)].groupby('stage'):
        e = _column_values(sub.drop_duplicates('wiring'), column)
        n = int(np.isfinite(e).sum())
        se = float(np.nanstd(e, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
        out[st] = dict(se=se, threshold=THRESHOLD_SE_MULT * se, n_wirings=n, group=group, column=column,
                       reading=str(sub.reading.iloc[0]) if 'reading' in sub.columns else None)
    if 'kc_soma' in stages and 'kc_soma' not in out:
        raise KeyError(f'no megamat kc_soma rows of group {group!r} in the stage table')
    return out


def e134_apl16_rows(step2_scores: pd.DataFrame, source: str = E134_APL16_SOURCE) -> pd.DataFrame:
    """The family's own runs on the APL16 input in E134's step 2 (source 'E90:B'): its somata on the APL16 panel, each
    condition, the default reading, the eight wirings (the guard's base at the APL16 somata)."""
    from kcq import stage_readings as R
    t = step2_scores
    t = t[(t.source == source) & (t.stage == 'kc_soma') & (t.panel == 'apl16')
          & (t.reading == R.DEFAULT_READINGS['kc_soma'].name) & t.wiring.isin(WIRINGS)]
    return t.drop_duplicates(['condition', 'wiring'])


def threshold_values(stage_table: pd.DataFrame, step2_scores: Optional[pd.DataFrame] = None) -> dict:
    """Every threshold step 2's score reads (the rulings of 15:41), from E134's tables: 'classes' (view_strength of
    each megamat stage in its default reading; the somata's is the classes' threshold), 'error' (the population error,
    for the guard), 'second_reading' (view_strength in the classes' second reading, SECOND_READINGS) and
    'apl16_somata' (per condition, view_strength and the population error of the family's APL16 somata in E134's step
    2, for the guard)."""
    rec = dict(rule=f'{THRESHOLD_SE_MULT:g} x SD (ddof 1) / sqrt(8) over the eight wirings in E134\'s runs', mult=THRESHOLD_SE_MULT,
               classes=class_thresholds(stage_table, stages=THRESHOLD_STAGES, column='view_strength'),
               error=class_thresholds(stage_table, stages=THRESHOLD_STAGES, column='error'),
               second_reading=class_thresholds(stage_table, stages=CLASS_STAGES, column='view_strength',
                                               readings=SECOND_READINGS),
               second_readings=dict(SECOND_READINGS), apl16_somata={})
    if step2_scores is not None and len(step2_scores):
        a = e134_apl16_rows(step2_scores)
        for cond, sub in a.groupby('condition'):
            rec['apl16_somata'][str(cond)] = {}
            for col in GUARD_COLUMNS:
                e = sub[col].to_numpy(dtype=float)
                n = int(np.isfinite(e).sum())
                se = float(np.nanstd(e, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
                rec['apl16_somata'][str(cond)][col] = dict(se=se, threshold=THRESHOLD_SE_MULT * se, n_wirings=n,
                                                           source=E134_APL16_SOURCE)
    return rec


def freeze_thresholds(root: Optional[Path] = None) -> dict:
    """The thresholds frozen before step 2's first arm (the rulings of 15:41): `threshold_values` of E134's
    score_stage_table.csv and step2_scores.csv, written with their sha256 to results/.../THRESHOLDS_FILE once. An
    existing file is returned as it is and never recomputed: the file present at step 2's start decides."""
    p = (out_dir() if root is None else Path(root)) / THRESHOLDS_FILE
    if p.exists():
        return _read_json(p)
    d = e134_results_dir()
    rec = dict(threshold_values(read_e134_table(E134_STAGE_TABLE), read_e134_table(E134_STEP2_SCORES)),
               e134_stage_table_sha256=_sha256(d / E134_STAGE_TABLE),
               e134_step2_scores_sha256=_sha256(d / E134_STEP2_SCORES), written=_now(), host=_host(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _write_json(p, rec)
    return rec


def read_thresholds(root: Optional[Path] = None) -> dict:
    """The thresholds step 2 froze at its start (`freeze_thresholds`); PermissionError if there are none (step 2 has
    not started)."""
    p = (out_dir() if root is None else Path(root)) / THRESHOLDS_FILE
    if not p.exists():
        raise PermissionError(f"{p}: no thresholds frozen yet; step 2 freezes them before its first arm")
    return _read_json(p)


def classify(scores: pd.DataFrame, thresholds: Optional[Mapping[str, Mapping[str, float]]] = None,
             soma_col: str = 'strength') -> dict:
    """The entry's classes from step 2's scores against E134's runs, on strength by odor (the adopted measure; the
    rulings of 15:41).

    "Upstream": the imposed dendrites lower the somata's strength-by-odor error by the threshold or more on six of
    eight wirings, by more than both permuted controls, and halve the boutons' strength-by-odor error. "Upstream to the
    boutons only": the bouton condition holds and the somata's strength-by-odor error stays within ± the threshold.
    (f): the imposed arm alone would be upstream, but the permuted controls gain as much. The carry classes, read only
    when "upstream" holds: R1 (or the kept rung, by more than its glomerulus-permuted control) gives at least half of
    the imposed fall of the strength-by-odor error at the boutons and at the somata.

    Args:
        scores: one row per (arm, wiring, stage) of the megamat input with columns arm, wiring, stage
            ('pn_bouton' / 'kc_soma'), strength (the strength-by-odor view's error) and strength_base (E134's run of
            the same wiring), error and error_base (the population error, read with `soma_col` 'error'), and excluded
            (no class: outside a band on the APL16 input, or more than 1 % of bouton responses above 300 Hz).
        thresholds: per stage {'threshold': ...} (`class_thresholds`); default the entry's fixed numbers.
        soma_col: the somata's column, 'strength' (the rulings) or 'error' (the entry's first text).

    Returns:
        dict with 'upstream', 'boutons_only', 'permuted_as_much', 'not_upstream' (bools), 'global_carry' and
        'terms_carry', the per-arm counts behind them and the thresholds.
    """
    s = scores[~scores.get('excluded', pd.Series(False, index=scores.index)).astype(bool)]
    thr_soma = SOMA_DROP_UPSTREAM if thresholds is None else float(thresholds['kc_soma']['threshold'])
    thr_within = SOMA_WITHIN_BOUTONS_ONLY if thresholds is None else float(thresholds['kc_soma']['threshold'])

    def drop(arm, stage, col):
        x = s[(s.arm == arm) & (s.stage == stage)].set_index('wiring')
        return x[f'{col}_base'] - x[col]

    def n_ok(mask):
        return int(np.sum(np.asarray(mask, dtype=bool)))

    out = {}
    soma = drop('imposed', 'kc_soma', soma_col)
    soma_g = drop('imposed_glom_perm', 'kc_soma', soma_col).reindex(soma.index)
    soma_o = drop('imposed_odor_perm', 'kc_soma', soma_col).reindex(soma.index)
    b = s[(s.arm == 'imposed') & (s.stage == 'pn_bouton')].set_index('wiring')
    halved = (b.strength <= BOUTON_STRENGTH_HALVED * b.strength_base) if len(b) else pd.Series(dtype=bool)
    n_soma_alone = n_ok(soma >= thr_soma)
    n_soma = n_ok((soma >= thr_soma) & (soma > soma_g) & (soma > soma_o))
    n_halved = n_ok(halved)
    out['upstream'] = bool(n_soma >= MIN_WIRINGS and n_halved >= MIN_WIRINGS)
    n_within = n_ok(soma.abs() <= thr_within)
    out['boutons_only'] = bool(not out['upstream'] and n_halved >= MIN_WIRINGS and n_within >= MIN_WIRINGS)
    out['permuted_as_much'] = bool(not out['upstream'] and not out['boutons_only'] and n_soma_alone >= MIN_WIRINGS
                                   and n_halved >= MIN_WIRINGS)
    out['not_upstream'] = bool(not out['upstream'] and not out['boutons_only'])
    gain_b = drop('imposed', 'pn_bouton', 'strength')
    r1_b, r1_s = drop('R1', 'pn_bouton', 'strength'), drop('R1', 'kc_soma', soma_col)
    rg_b, rg_s = drop('rung', 'pn_bouton', 'strength'), drop('rung', 'kc_soma', soma_col)
    rp_b, rp_s = drop('rung_glom_perm', 'pn_bouton', 'strength'), drop('rung_glom_perm', 'kc_soma', soma_col)
    half = lambda x, ref: (x >= CARRY_SHARE * ref.reindex(x.index))
    n_r1 = n_ok(half(r1_b, gain_b) & half(r1_s, soma)) if len(r1_b) else 0
    out['global_carry'] = bool(out['upstream'] and n_r1 >= MIN_WIRINGS)
    n_rung = (n_ok(half(rg_b, gain_b) & half(rg_s, soma) & (rg_b > rp_b.reindex(rg_b.index))
                   & (rg_s > rp_s.reindex(rg_s.index))) if len(rg_b) else 0)
    out['terms_carry'] = bool(out['upstream'] and not out['global_carry'] and n_rung >= MIN_WIRINGS)
    out['counts'] = dict(soma_upstream=n_soma, soma_alone=n_soma_alone, boutons_halved=n_halved, soma_within=n_within,
                         r1_carry=n_r1, rung_carry=n_rung)
    out['thresholds'] = dict(soma_drop=thr_soma, soma_within=thr_within, bouton_halved=BOUTON_STRENGTH_HALVED,
                             reset=thresholds is not None, soma_col=soma_col)
    return out


def class_of(c: Mapping) -> str:
    """The upstream class of `classify`'s result: 'upstream', 'boutons_only', 'permuted_as_much' or 'not_upstream'."""
    for k in ('upstream', 'boutons_only', 'permuted_as_much'):
        if c[k]:
            return k
    return 'not_upstream'


def outcome_letter(classes: Mapping) -> str:
    """The outcome's letter of `classify`'s classes, before the guards: (a) to (f)."""
    c = class_of(classes)
    if c == 'upstream':
        return 'a' if classes['global_carry'] else 'b' if classes['terms_carry'] else 'c'
    return 'd' if c == 'boutons_only' else 'f' if c == 'permuted_as_much' else 'e'


def dff_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Step 2's scores with the dF/F repeat's arms (DFF_REPEAT_ARMS) in place of the three imposed arms, for
    `classify` (the other arms kept)."""
    rest = scores[~scores.arm.isin(IMPOSED_ARMS) & ~scores.arm.isin(DFF_REPEAT_ARMS)]
    dff = scores[scores.arm.isin(DFF_REPEAT_ARMS)].copy()
    dff['arm'] = dff.arm.str[:-len(DFF_SUFFIX)]
    return pd.concat([rest, dff], ignore_index=True)


def outcome(classes: Mapping, classes_dff: Optional[Mapping] = None, classes_second: Optional[Mapping] = None,
            guard_rec: Optional[Mapping] = None) -> dict:
    """The entry's outcome from `classify`'s classes: (a) upstream and the global parameters carry, (b) upstream and the
    glomerulus terms carry, (c) upstream and no fitted stage carries, (d) upstream to the boutons only, (e) not upstream,
    (f) the permuted controls gain as much as the imposed dendrites.

    Before (e) or (f) is declared, the class read with the dF/F repeat (`classes_dff`, `classify` of `dff_scores`) is
    compared (ruling 1's guard): without it the outcome waits on the repeat ('needs the dF/F repeat'), and if it
    differs the outcome is 'unit-dependent' (nothing is closed; it goes to Tom and Yang). The classes read in the second
    reading (`classes_second`, SECOND_READINGS) are compared too: an outcome that differs is 'reading_dependent', the
    default reading decides, and (e) or (f) is not declared on it ('reading-dependent'). Rule 3's guard
    (`guard_rec`, `guard`) is carried: a raised guard of the arm the outcome hands on (R1 for (a), the rung for (b)) is
    named ('at_a_cost_beside', with the views and stages) and that stage is not called an improvement."""
    c = class_of(classes)
    o = outcome_letter(classes)
    res = dict(outcome=o, cls=c)
    if classes_second is not None:
        o2 = outcome_letter(classes_second)
        res.update(cls_second=class_of(classes_second), outcome_second=o2, reading_dependent=bool(o2 != o))
    if o in ('e', 'f'):
        if classes_dff is None:
            res.update(outcome='needs the dF/F repeat', pending=o)
        else:
            cd = class_of(classes_dff)
            res['cls_dff'] = cd
            if cd != c:
                res.update(outcome='unit-dependent', pending=o)
            elif res.get('reading_dependent'):
                res.update(outcome='reading-dependent', pending=o)
    if guard_rec is not None:
        res['guard_raised'] = {a: bool(g.get('raised')) for a, g in guard_rec.items()}
        arm = {'a': 'R1', 'b': 'rung'}.get(o)
        if arm is not None:
            g = guard_rec.get(arm, {})
            res['improvement'] = not bool(g.get('raised'))
            if g.get('raised'):
                res['at_a_cost_beside'] = list(g.get('raised_items', []))
    return res


# ---------------------------------------------------------------------------
# step 2's score

def step2_groups(dirs: Optional[Sequence[Path]] = None) -> Dict[Tuple[str, str], list]:
    """Step 2's saved runs (`step2_dirs`, stage_readings' `load_step3_run`) by (input panel, arm), each group's runs in
    WIRINGS order; points not run yet are left out."""
    from kcq import stage_readings as R
    groups: Dict[Tuple[str, str], list] = {}
    for arm, w in planned_step2(with_rung=True, with_dff_repeat=True):
        try:
            run = R.load_step3_run(tag_of(arm, w), dirs=step2_dirs() if dirs is None else dirs)
        except FileNotFoundError:
            continue
        groups.setdefault((arm.panel, arm.arm), []).append(run)
    return groups


def _flies():
    from kcq import stage_readings as R
    return R.FlyCache()


def score_group(runs: Sequence, flies, lock134: Optional[Mapping] = None, readings: Optional[Sequence[str]] = None,
                n_floor: Optional[int] = None) -> pd.DataFrame:
    """One group's runs read as E134's score reads its runs (atlas's `read_run`, each wiring as one fly with its seed
    index, E83's wiring order) in each stage's default reading and its alternatives (or those named), and scored with
    hal's measures as E134 scores them (`score_family`: each wiring's population error and views, their mean and the
    common error, beside the flies' floor and one fly against the others)."""
    from experiments import E134_population_baseline as E134R
    from kcq import stage_readings as R
    stages = [s for s in R.STAGES if s in runs[0].stages()]
    out = []
    for stage, rd in E134R._reading_sets(stages, readings):
        if rd.curve_scope == 'family':
            fam = R.read_family(list(runs), stages=[stage], readings={stage: rd}, lock=lock134, flies=flies)
        else:
            fam = {r.wiring: R.read_run(r, stages=[stage], readings={stage: rd}, seed_index=WIRINGS.index(r.wiring),
                                        lock=lock134, flies=flies) for r in runs}
        sc = E134R.score_family(fam, rd.name, flies, E134R.N_FLOOR if n_floor is None else int(n_floor))
        if len(sc):
            out.append(sc)
        gc.collect()
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def run_bands(o: Mapping[str, tuple]) -> dict:
    """An APL16 run's realism (E83's `realism`, E77's `realism_of`): the response rates and their ratio, the response
    rate ρ, the odors-per-KC distances and 'realism_flag' (empty inside the four bands)."""
    from experiments import E83_cross_wiring_fit as E83
    return E83.realism(dict(o))


def apl16_effect(o: Mapping[str, tuple]) -> dict:
    """The APL-effect statistics of an APL16 run (a held-out test, reported and not gated): `apl_effect_stats_readout`
    on raw counts ('raw_') and through the calibrated readout, K set on the run's APL-on counts ('cal_', E83's
    `readout_K`): the loss, its mean and pattern parts, the amplitude and shape parts, ρ and Pearson r."""
    from experiments import E83_cross_wiring_fit as E83
    on, off = o['on'][1], o['off'][1]
    row = {f'raw_{k}': v for k, v in H.apl_effect_stats_readout(on, off, 'raw').items()}
    K = E83.readout_K(dict(o))
    row.update({f'cal_{k}': v for k, v in H.apl_effect_stats_readout(on, off, 'saturating', K).items()}, cal_K=K)
    return row


def _outputs(run) -> Dict[str, tuple]:
    return {c: (run.conditions[c]['responses'], run.conditions[c]['spike_counts']) for c in ('on', 'off')}


def totals_cv(run, glomeruli: Sequence[str], odors: Sequence[str] = tuple(H.MEGAMAT12)) -> dict:
    """A description beside (E119's and E120's statistic, here on the model's own units, neither sampled nor pooled):
    the coefficient of variation over `odors` (default the 12 megamat odors both inputs have) of each odor's mean
    change (window minus settled pre-odor level, signed), APL on, at the PN dendrites of `glomeruli` (the shared 28)
    and at the PN boutons of those glomeruli. The entry's expectations are given in it (the family 0.074 against the
    flies' 0.312 at the dendrites, 0.103 against 0.375 at the boutons)."""
    c = run.conditions['on']
    rec = {}
    gl = set(str(g) for g in glomeruli)
    for st in ('pn', 'bouton'):
        if f'{st}_win' not in c:
            continue
        x = c[f'{st}_win'] - c[f'{st}_pre']
        od = [o for o in odors if o in x.columns]
        if st == 'pn':
            keep = np.array([str(g) in gl for g in x.index])
        else:
            keep = np.array([str(g) in gl for g in np.asarray(run.meta['bouton_glomerulus'])])
        tot = x.loc[keep, od].mean(axis=0).to_numpy(dtype=float)
        rec[f'totals_cv_{st}'] = float(np.std(tot) / np.mean(tot)) if len(tot) and np.mean(tot) != 0 else np.nan
    return rec


def run_realism(run, glomeruli: Sequence[str], lock134: Optional[Mapping] = None) -> dict:
    """A step-2 run's realism, reported only: the overall response rates with APL on and off and their ratio, the
    largest spike counts, the parameters after the re-hold (threshold and APL scales) beside the fit's threshold, the
    largest PN and bouton rates and the share of bouton responses above 300 Hz (`stage_rates`), the totals' CV
    (`totals_cv`); on the APL16 input the four bands (`run_bands`), on megamat the flies' held rate beside."""
    from experiments import E134_population_baseline as E134R
    info = run.info if isinstance(run.info, Mapping) else {}
    row = dict(tag=run.tag, input_panel=run.input_panel, arm=info.get('arm'), wiring=run.wiring,
               route=info.get('route'))
    row.update(E134R._counts_row(run.conditions))
    p = info.get('params') or {}
    row.update({k: p.get(k, np.nan) for k in E134R.PARAM_KEYS})
    row.update(stage_rates(run.conditions))
    row.update(totals_cv(run, glomeruli))
    if run.input_panel == 'apl16':
        b = run_bands(_outputs(run))
        row.update({k: v for k, v in b.items() if k not in row})
        row['outside_band'] = bool(b.get('realism_flag'))
    else:
        row['realism_flag'] = ''
        row['outside_band'] = False
        if run.input_panel == 'megamat' and lock134 is not None:
            row['rate_on_flies'] = E134R.held_rate('megamat', lock134)
    return row


def exclusions(realism: pd.DataFrame) -> Dict[Tuple[str, str], List[str]]:
    """(arm, wiring) of the megamat input -> why it gets no class (the entry's): its APL16 run (the arms run on both
    inputs: R1, the rung and the imposed dendrites) is outside a band, or its own megamat run has more than 1 % of its
    bouton responses above 300 Hz (`stage_rates`)."""
    out: Dict[Tuple[str, str], List[str]] = {}
    for r in realism.to_dict('records'):
        key = (str(r['arm']), str(r['wiring']))
        if r['input_panel'] == 'apl16' and bool(r.get('outside_band')):
            out.setdefault(key, []).append(f"its APL16 run is outside a band ({r.get('realism_flag')})")
        if r['input_panel'] == 'megamat' and bool(r.get('bouton_cap_excluded')):
            out.setdefault(key, []).append(f"{r['bouton_share_above_300']:.4f} of its bouton responses above "
                                           f'{BOUTON_RATE_CAP_HZ:.0f} Hz')
    return out


def class_rows(scores: pd.DataFrame, base: pd.DataFrame,
               excluded: Optional[Mapping[Tuple[str, str], Sequence[str]]] = None,
               readings: Optional[Mapping[str, str]] = None) -> pd.DataFrame:
    """The rows `classify` and `guard` read, one per (arm, wiring, stage) of the megamat input (the classes' stages and
    the PN dendrites), the megamat panel, APL on, each stage's default reading (or `readings`, the second reading):
    'strength' (the strength-by-odor view's error, the classes' measure), 'error' (the population error, three
    views), the columns reported beside (BESIDE_COLUMNS: the two-view error, sharing, breadth, the matrix view), E134's
    same wiring beside (`e134_base` in the same readings) with the changes ('<column>_change'), the falls ('drop' of
    the population error, 'strength_drop'), 'strength_ratio', and 'excluded' with its reasons (`exclusions`)."""
    from kcq import stage_readings as R
    excluded = excluded or {}
    want = {st: dict(readings or {}).get(st, R.DEFAULT_READINGS[st].name) for st in R.DEFAULT_READINGS}
    s = scores[(scores.input_panel == 'megamat') & (scores.panel == 'megamat') & (scores.condition == 'on')
               & scores.wiring.isin(WIRINGS) & scores.stage.isin(CLASS_STAGES + ('pn_dendrite',))]
    s = s[np.array([rd == want.get(st) for st, rd in zip(s.stage, s.reading)], dtype=bool)]
    keep = [c for c in ('error', 'error_with_matrix') if c in s.columns] + [c for c in s.columns if c.startswith('view_')]
    out = s[['arm', 'wiring', 'stage', 'reading'] + keep].copy()
    out['strength'] = out['view_strength'] if 'view_strength' in out.columns else np.nan
    try:
        out['error_strength_sharing'] = _column_values(out, 'error_strength_sharing')
    except KeyError:
        out['error_strength_sharing'] = np.nan
    out = out.merge(base, on=['stage', 'wiring'], how='left')
    out['drop'] = out.error_base - out.error
    out['strength_drop'] = out.strength_base - out.strength
    out['strength_ratio'] = out.strength / out.strength_base
    for c in BESIDE_COLUMNS:
        b = dict(error='error_base', view_strength='strength_base').get(c, f'{c}_base')
        if c in out.columns and b in out.columns:
            out[f'{c}_change'] = out[c] - out[b]
    why = [list(excluded.get((str(a), str(w)), ())) for a, w in zip(out.arm, out.wiring)]
    out['excluded'] = [bool(x) for x in why]
    out['excluded_why'] = ['; '.join(x) for x in why]
    return out.reset_index(drop=True)


def stage_error_measure(t: pd.DataFrame, measure: str, stages: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The megamat stage error under `measure` (`measure_values`; 'all' with the PN dendrites on four views, E134's
    definition) per (group, wiring) of a scores table or E134's stage table: the root mean square over `stages` (default
    E134's four megamat stages, MEGAMAT_ERROR_STAGES; NaN where one is missing), each in its default reading, APL on,
    the eight wirings; and per group the mean over its wirings (wiring 'model')."""
    from experiments import E134_population_baseline as E134R
    from kcq import stage_readings as R
    st = tuple(E134R.MEGAMAT_ERROR_STAGES if stages is None else stages)
    s = t[(t.panel == 'megamat') & (t.condition == 'on') & t.wiring.isin(WIRINGS) & t.stage.isin(st)]
    s = s[np.array([rd == R.DEFAULT_READINGS[x].name for x, rd in zip(s.stage, s.reading)], dtype=bool)]
    if not len(s):
        return pd.DataFrame(columns=['group', 'wiring', 'value'])
    s = s.drop_duplicates(['group', 'wiring', 'stage'])
    s = s.assign(_v=measure_values(s, measure, matrix_for_dendrites=True))
    rows = []
    for (g, w), x in s.groupby(['group', 'wiring'], sort=False):
        v = x._v.to_numpy(dtype=float)
        rows.append(dict(group=g, wiring=w, value=float(np.sqrt(np.mean(v ** 2))) if len(v) == len(st) else np.nan))
    out = pd.DataFrame(rows)
    model = out.groupby('group', sort=False).value.mean().reset_index().assign(wiring='model')
    return pd.concat([out, model[['group', 'wiring', 'value']]], ignore_index=True)


# the megamat stage error's forms (the rulings of 15:41): on strength over the four stages (the adopted measure), the
# three-view error beside (E134's definition), the two-view error beside; for the imposed arms both also over the
# boutons, claws and somata, since their dendrites are the flies' by construction
STAGE_ERROR_FORMS: Dict[str, Tuple[str, Optional[Tuple[str, ...]]]] = {
    'stage_error_strength': ('strength', None),
    'stage_error_three_view': ('all', None),
    'stage_error_strength_sharing': ('strength_sharing', None),
    'stage_error_strength_3': ('strength', CLASS_STAGES),
    'stage_error_three_view_3': ('all', CLASS_STAGES),
}


def stage_error_table(scores: pd.DataFrame, stage_table: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """The megamat stage error of every group of step 2 in each form of STAGE_ERROR_FORMS (`stage_error_measure`), per
    wiring and their mean ('model'), beside E134's held runs of the same wiring ('<form>_base', from `stage_table`,
    its score_stage_table.csv, group 'megamat_held') and the change (lower is better). The three-stage forms are given
    for the imposed arms only. On the APL16 input the megamat stages are read on the 12 megamat odors it has."""
    if not len(scores):
        return pd.DataFrame()
    keys = scores[['input_panel', 'arm', 'group']].drop_duplicates()
    out = None
    for form, (m, stages) in STAGE_ERROR_FORMS.items():
        try:
            v = stage_error_measure(scores, m, stages)
        except KeyError:
            continue
        v = v.rename(columns={'value': form})
        if stage_table is not None and len(stage_table):
            try:
                b = stage_error_measure(stage_table[stage_table.group == E134_HELD_GROUP], m, stages)
                v[f'{form}_base'] = v.wiring.map(b.set_index('wiring').value)
                v[f'{form}_change'] = v[form] - v[f'{form}_base']
            except KeyError:
                pass
        out = v if out is None else out.merge(v, on=['group', 'wiring'], how='outer')
    if out is None:
        return pd.DataFrame()
    out = keys.merge(out, on='group', how='right')
    imposed = out.arm.map(lambda a: _base_arm(Arm('megamat', str(a))) in IMPOSED_ARMS + (REPEAT_ARM.arm,))
    for c in [c for c in out.columns if c.startswith(('stage_error_strength_3', 'stage_error_three_view_3'))]:
        out.loc[~imposed, c] = np.nan
    order = {w: i for i, w in enumerate(list(WIRINGS) + ['model'])}
    return out.sort_values(['input_panel', 'arm', 'wiring'], key=lambda x: x.map(order) if x.name == 'wiring' else x,
                           kind='mergesort').reset_index(drop=True)


def routes_check(rows: pd.DataFrame, thresholds: Mapping, runs: Optional[pd.DataFrame] = None) -> dict:
    """Ruling 9 on the adopted measure: the route-2 repeat against its route-1 run ('imposed' on the first wiring),
    each megamat stage's strength-by-odor error (`routes_agree`, at the frozen view_strength thresholds), with the
    population error's difference beside, and the solve's record from the runs table. No class rests on it."""
    w = WIRINGS[0]
    a = rows[(rows.arm == 'imposed') & (rows.wiring == w)].set_index('stage')
    b = rows[(rows.arm == REPEAT_ARM.arm) & (rows.wiring == w)].set_index('stage')
    if not len(a) or not len(b):
        return dict(available=False)
    res = routes_agree({k: float(v) for k, v in a.strength.items()}, {k: float(v) for k, v in b.strength.items()},
                       thresholds)
    res['available'] = True
    res['column'] = 'view_strength'
    res['population_error'] = {st: dict(route1=float(a.error.get(st, np.nan)), route2=float(b.error.get(st, np.nan)),
                                        diff=float(abs(a.error.get(st, np.nan) - b.error.get(st, np.nan))))
                               for st in CLASS_STAGES}
    if runs is not None and len(runs) and 'tag' in runs.columns:
        r = runs[runs.tag == tag_of(REPEAT_ARM, w)]
        if len(r):
            res['solve'] = {k: r.iloc[0][k] for k in r.columns if str(k).startswith(('solve_', 'ln_'))}
    return res


def apl16_guard_rows(scores: pd.DataFrame, step2_scores: Optional[pd.DataFrame]) -> pd.DataFrame:
    """The guard's rows at the APL16 somata: for R1 and the kept rung on the APL16 input, per wiring and condition, the
    somata's strength-by-odor error and population error (the APL16 panel, the default reading) beside the family's own
    run of the same wiring in E134's step 2 (`e134_apl16_rows`)."""
    from kcq import stage_readings as R
    if not len(scores) or step2_scores is None or not len(step2_scores):
        return pd.DataFrame()
    s = scores[(scores.input_panel == 'apl16') & scores.arm.isin(GUARD_ARMS) & (scores.stage == 'kc_soma')
               & (scores.panel == 'apl16') & (scores.reading == R.DEFAULT_READINGS['kc_soma'].name)
               & scores.wiring.isin(WIRINGS)]
    b = e134_apl16_rows(step2_scores)[['condition', 'wiring'] + list(GUARD_COLUMNS)]
    b = b.rename(columns={c: f'{c}_base' for c in GUARD_COLUMNS})
    return s[['arm', 'wiring', 'condition'] + list(GUARD_COLUMNS)].merge(b, on=['condition', 'wiring'], how='left')


def step1_rung_of_arm(arm: str, highest_kept: Optional[str]) -> Optional[str]:
    """The step-1 rung behind a step-2 arm: R1's is 'R1', the rung's is step 1's highest kept rung ('R2m', 'R2f' or
    'R3')."""
    return 'R1' if arm == 'R1' else highest_kept if arm == 'rung' else None


def _strength_by_direction(step1_rungs: pd.DataFrame, rung: str) -> Dict[str, float]:
    r = step1_rungs[step1_rungs.rung == rung]
    if not len(r):
        return {}
    return {c[:-len('_view_strength')]: float(r.iloc[0][c]) for c in step1_rungs.columns
            if c.endswith('_view_strength') and '->' in c}


def step1_dendrite_guard(step1_rungs: Optional[pd.DataFrame], rung: Optional[str]) -> dict:
    """The guard at the dendrites (the rulings): raised where the rung's held-out strength-by-odor error in step 1
    (step1_rungs.csv, '<direction>_view_strength') is above R0's in both directions."""
    if step1_rungs is None or not len(step1_rungs) or rung is None:
        return dict(available=False, stage='pn_dendrite', column='view_strength')
    x, r0 = _strength_by_direction(step1_rungs, rung), _strength_by_direction(step1_rungs, 'R0')
    if not x or not r0:
        return dict(available=False, stage='pn_dendrite', column='view_strength', rung=rung)
    above = {d: bool(x[d] > r0[d]) for d in x if d in r0}
    return dict(available=True, stage='pn_dendrite', column='view_strength', rung=rung, values=x, r0=r0, above=above,
                raised=bool(above and all(above.values())))


def step1_strength_names(step1_rungs: Optional[pd.DataFrame]) -> List[dict]:
    """Step 1 read on the adopted measure (the rulings): each kept rung above R0 whose held-out strength-by-odor error
    is above that of the highest kept rung below it ('below') in both directions, named "kept on the matrix, higher on
    strength". No refit: it is still step 2's rung."""
    out = []
    if step1_rungs is None or not len(step1_rungs) or 'kept' not in step1_rungs.columns:
        return out
    for r in step1_rungs.to_dict('records'):
        kept = r.get('kept') is True or str(r.get('kept')) == 'True'
        below = r.get('below')
        if r['rung'] == 'R0' or not kept or not isinstance(below, str):
            continue
        x, b = _strength_by_direction(step1_rungs, r['rung']), _strength_by_direction(step1_rungs, below)
        if x and b and all(x[d] > b[d] for d in x if d in b):
            out.append(dict(rung=r['rung'], below=below, values=x, below_values=b,
                            name='kept on the matrix, higher on strength'))
    return out


def guard(rows: pd.DataFrame, thresholds: Mapping, apl16: Optional[pd.DataFrame] = None,
          step1_rungs: Optional[pd.DataFrame] = None, highest_kept: Optional[str] = None,
          arms: Sequence[str] = GUARD_ARMS) -> dict:
    """Rule 3's guard beside the classes (the rulings of 15:41), for R1 and the kept rung against the family's run of
    the same wiring: raised where, on six of eight wirings, the strength-by-odor error or the population error rises
    by more than three standard errors over the eight wirings of the same quantity in E134's runs (the frozen
    thresholds: 'classes' and 'error' at the megamat boutons, claws and somata, from `rows`, `class_rows` in the
    default readings; 'apl16_somata' at the APL16 somata with APL on and off, from `apl16`, `apl16_guard_rows`); at the
    dendrites where step 1's held-out strength-by-odor error is above R0's in both directions (`step1_dendrite_guard`).
    Every wiring with a run counts (the classes' exclusions do not apply to a cost).

    Returns:
        {arm: dict(raised, items (one per stage, panel, condition and column: the number of wirings risen beyond the
        threshold, the mean rise), dendrites, raised_items (names of the raised ones))}.
    """
    thr = {'view_strength': thresholds.get('classes', {}), 'error': thresholds.get('error', {})}
    out = {}
    for arm in arms:
        items = []
        x = rows[rows.arm == arm]
        for st in GUARD_STAGES:
            y = x[x.stage == st]
            for col, c in (('view_strength', 'strength'), ('error', 'error')):
                t = thr[col].get(st, {}).get('threshold', np.nan)
                rise = (y[c] - y[f'{c}_base']).to_numpy(dtype=float) if len(y) else np.array([])
                n = int(np.sum(rise > t)) if np.isfinite(t) else 0
                items.append(dict(stage=st, panel='megamat', condition='on', column=col, threshold=float(t), n_risen=n,
                                  n_wirings=int(np.isfinite(rise).sum()),
                                  mean_rise=float(np.nanmean(rise)) if np.isfinite(rise).any() else np.nan,
                                  raised=bool(n >= MIN_WIRINGS)))
        if apl16 is not None and len(apl16):
            z = apl16[apl16.arm == arm]
            for cond in ('on', 'off'):
                zz = z[z.condition == cond]
                for col in GUARD_COLUMNS:
                    t = thresholds.get('apl16_somata', {}).get(cond, {}).get(col, {}).get('threshold', np.nan)
                    rise = (zz[col] - zz[f'{col}_base']).to_numpy(dtype=float) if len(zz) else np.array([])
                    n = int(np.sum(rise > t)) if np.isfinite(t) else 0
                    items.append(dict(stage='kc_soma', panel='apl16', condition=cond, column=col, threshold=float(t),
                                      n_risen=n, n_wirings=int(np.isfinite(rise).sum()),
                                      mean_rise=float(np.nanmean(rise)) if np.isfinite(rise).any() else np.nan,
                                      raised=bool(n >= MIN_WIRINGS)))
        dend = step1_dendrite_guard(step1_rungs, step1_rung_of_arm(arm, highest_kept))
        names = [f"{i['column']} at {i['stage']} ({i['panel']}, APL {i['condition']})" for i in items if i['raised']]
        if dend.get('raised'):
            names.append('view_strength at pn_dendrite (step 1, held out, above R0 in both directions)')
        out[arm] = dict(raised=bool(names), items=items, dendrites=dend, raised_items=names)
    return out


def imposed_check_table(runs: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Ruling 8 per imposed run: the largest difference between its table and its PN rates in E134's reading
    (recorded by the runner), named above IMPOSED_CHECK_HZ, beside its PN-dendrite matrix error (the matrix view in
    one-fly units, the default reading, APL on, megamat panel) and the flies' floor note."""
    from kcq import stage_readings as R
    if not len(runs) or 'imposed_check_max_hz' not in runs.columns:
        return pd.DataFrame()
    r = runs[runs.imposed_check_max_hz.notna()].copy()
    if not len(scores):
        return r
    s = scores[(scores.stage == 'pn_dendrite') & (scores.panel == 'megamat') & (scores.condition == 'on')
               & (scores.reading == R.DEFAULT_READINGS['pn_dendrite'].name) & scores.wiring.isin(WIRINGS)]
    vm = s['view_matrix'] if 'view_matrix' in s.columns else pd.Series(np.nan, index=s.index)
    m = {(p, a, w): v for p, a, w, v in zip(s.input_panel, s.arm, s.wiring, vm)}
    r['pn_matrix_error'] = [m.get((p, a, w), np.nan) for p, a, w in zip(r.panel, r.arm, r.wiring)]
    r['flies_floor_note'] = PN_MATRIX_FLOOR_NOTE
    keep = ['tag', 'panel', 'arm', 'wiring', 'route', 'imposed_check_max_hz', 'imposed_check_named', 'pn_matrix_error',
            'flies_floor_note', 'pn_check_max_set_hz', 'clip_n_clipped', 'clip_share', 'clip_named']
    return r[[c for c in keep if c in r.columns]].reset_index(drop=True)


def apl16_effect_table(groups: Mapping[Tuple[str, str], Sequence], with_family: bool = True) -> pd.DataFrame:
    """The APL-effect statistics (`apl16_effect`) of every APL16 run of step 2, and beside them the family's own fit
    on each wiring (arm 'family', the saved fit, `family_saved`): a held-out test, reported and not gated."""
    rows = []
    for (panel, arm), runs in groups.items():
        if panel != 'apl16':
            continue
        for r in runs:
            rows.append(dict(arm=arm, wiring=r.wiring, tag=r.tag, **apl16_effect(_outputs(r))))
    if with_family and rows:
        for w in WIRINGS:
            rows.append(dict(arm='family', wiring=w, tag=f'fit_{w}', **apl16_effect(family_saved(w))))
    return pd.DataFrame(rows)


def drops_summary(rows: pd.DataFrame) -> pd.DataFrame:
    """Per arm and stage: the mean over wirings (not excluded) of the falls of the population error and of the
    strength-by-odor view's error from E134's held runs, the number of wirings and of excluded wirings."""
    if not len(rows):
        return pd.DataFrame()
    out = []
    for (arm, st), x in rows.groupby(['arm', 'stage'], sort=False):
        k = x[~x.excluded]
        out.append(dict(arm=arm, stage=st, n_wirings=int(len(k)), n_excluded=int(x.excluded.sum()),
                        drop_mean=float(k['drop'].mean()) if len(k) else np.nan,
                        strength_drop_mean=float(k.strength_drop.mean()) if len(k) else np.nan,
                        strength_ratio_mean=float(k.strength_ratio.mean()) if len(k) else np.nan))
    return pd.DataFrame(out)


def score_tables(scores: pd.DataFrame, stage_table: pd.DataFrame, realism: pd.DataFrame,
                 runs: Optional[pd.DataFrame] = None, with_rung: bool = True, thresholds: Optional[Mapping] = None,
                 step2_scores: Optional[pd.DataFrame] = None, step1_rungs: Optional[pd.DataFrame] = None,
                 highest_kept: Optional[str] = None) -> dict:
    """Step 2's classes from its scores (`score_group` per group, with columns input_panel and arm) against E134's held
    runs (`stage_table`, E134's score_stage_table.csv), on strength by odor (the rulings of 15:41):
    - the frozen thresholds (`thresholds`, `freeze_thresholds`; computed from the tables, `threshold_values`, when not
      given), the exclusions (`exclusions` of `realism`, `run_realism` per run), the class rows (`class_rows`, every
      view and the population and two-view errors beside), the classes (`classify` on view_strength), the dF/F
      repeat's classes once its 24 runs exist, the classes in the second reading (SECOND_READINGS);
    - rule 3's guard for R1 and the kept rung (`guard`: the megamat boutons, claws and somata, the APL16 somata from
      E134's step 2 `step2_scores`, and the dendrites from step 1's `step1_rungs`);
    - the outcome (`outcome`, carrying the dF/F and second-reading checks and the guard; 'incomplete' until every
      planned point has a run), step 1 read on strength (`step1_strength_names`), the routes (`routes_check`), the
      megamat stage error in its forms (`stage_error_table`) and the imposed checks (`imposed_check_table`).

    Returns:
        dict of 'verdict' (JSON-able), 'class_rows', 'class_rows_second', 'drops', 'stage_error', 'imposed_checks',
        'guard_items'.
    """
    thr = dict(thresholds) if thresholds is not None else threshold_values(stage_table, step2_scores)
    cls_thr = thr['classes']
    excl = exclusions(realism) if len(realism) else {}
    have = set(realism.tag) if len(realism) else set()
    planned = planned_step2(with_rung=with_rung)
    missing = [tag_of(a, w) for a, w in planned if tag_of(a, w) not in have]
    dff_missing = [tag_of(Arm('megamat', a), w) for a in DFF_REPEAT_ARMS for w in WIRINGS
                   if tag_of(Arm('megamat', a), w) not in have]
    empty = pd.DataFrame(columns=['arm', 'wiring', 'stage', 'error', 'error_base', 'strength', 'strength_base',
                                  'excluded'])
    rows = class_rows(scores, e134_base(stage_table), excl) if len(scores) else empty
    classes = classify(rows[rows.arm.isin(MEGAMAT_ARMS)], cls_thr)
    classes_dff = None
    if not dff_missing:
        classes_dff = classify(dff_scores(rows[rows.arm.isin(MEGAMAT_ARMS + DFF_REPEAT_ARMS)]), cls_thr)
    rows2 = (class_rows(scores, e134_base(stage_table, readings=SECOND_READINGS), excl, SECOND_READINGS)
             if len(scores) else empty)
    second = rows2[rows2.arm.isin(MEGAMAT_ARMS)]
    classes_second = None
    if len(second) and set(SECOND_READINGS) <= set(second.stage) and 'kc_soma' in thr.get('second_reading', {}):
        classes_second = classify(second, thr['second_reading'])
    g = guard(rows, thr, apl16_guard_rows(scores, step2_scores), step1_rungs, highest_kept)
    oc = outcome(classes, classes_dff, classes_second, g)
    if missing:
        oc = dict(outcome='incomplete', n_missing=len(missing), provisional=oc)
    stage_err = stage_error_table(scores, stage_table)
    imposed = imposed_check_table(runs if runs is not None else pd.DataFrame(), scores)
    drops = drops_summary(rows)
    se_mean = {}
    if len(stage_err):
        m = stage_err[stage_err.wiring == 'model']
        cols = [c for c in m.columns if c.startswith('stage_error_')]
        se_mean = {f'{p}_{a}': {c: (None if not np.isfinite(v) else float(v)) for c, v in zip(cols, r)}
                   for p, a, r in zip(m.input_panel, m.arm, m[cols].to_numpy(dtype=float))}
    named = {}
    if runs is not None and len(runs):
        for col in ('clip_named', 'imposed_check_named'):
            if col in runs.columns:
                named[col] = sorted(str(t) for t, v in zip(runs.tag, runs[col]) if v is True or str(v) == 'True')
        if 'scale_named' in runs.columns:
            named['scale_named'] = {str(t): str(v) for t, v in zip(runs.tag, runs.scale_named)
                                    if isinstance(v, str) and v}
    items = pd.DataFrame([dict(arm=a, **i) for a, x in g.items() for i in x['items']])
    verdict = dict(complete=not missing, n_runs=len(have), n_planned=len(planned), missing=missing,
                   dff_repeat=dict(complete=not dff_missing, n_missing=len(dff_missing)),
                   measure=dict(column=CLASS_COLUMN, reading="each stage's default",
                                rule='strength by odor, adopted after E134\'s gate (the rulings of 15:41)'),
                   with_rung=bool(with_rung), thresholds=cls_thr,
                   thresholds_frozen={k: thr.get(k) for k in ('e134_stage_table_sha256', 'e134_step2_scores_sha256',
                                                               'written', 'host')},
                   classes=classes, cls=class_of(classes), classes_dff=classes_dff, classes_second=classes_second,
                   second_readings=dict(SECOND_READINGS), outcome=oc,
                   guard={a: dict(raised=x['raised'], raised_items=x['raised_items'], dendrites=x['dendrites'])
                          for a, x in g.items()},
                   step1_strength=step1_strength_names(step1_rungs), routes=routes_check(rows, cls_thr, runs),
                   excluded={f'{a}|{w}': v for (a, w), v in excl.items()},
                   megamat_stage_error=se_mean, named=named,
                   realism_flags={str(t): str(f) for t, f in zip(realism.get('tag', []), realism.get('realism_flag', []))
                                  if isinstance(f, str) and f})
    return dict(verdict=verdict, class_rows=rows, class_rows_second=rows2, drops=drops, stage_error=stage_err,
                imposed_checks=imposed, guard_items=items)


def write_score_table(res: Mapping, path: Path) -> None:
    """score_table.md: the outcome and classes (strength by odor), the thresholds, the guard, the second reading, the
    falls per arm and stage with the population error beside, the megamat stage error per arm against E134's held runs,
    the routes' agreement and the exclusions."""
    v = res['verdict']
    oc = v['outcome']
    L = ['# E135 step 2: the classes on strength by odor (against E134\'s runs held at the megamat flies\' response '
         'rate)', '',
         f"- outcome: {oc}", f"- class: {v['cls']}; complete: {v['complete']} ({v['n_runs']} runs of "
                             f"{v['n_planned']} planned, dF/F repeat complete: {v['dff_repeat']['complete']})",
         "- thresholds (strength by odor, one-fly units, frozen): " + ', '.join(
             f"{k} {d['threshold']:.4f}" for k, d in v['thresholds'].items()),
         f"- counts: {v['classes']['counts']}",
         f"- second reading ({v['second_readings']}): "
         f"{class_of(v['classes_second']) if v['classes_second'] else 'not read'}",
         '- guard: ' + '; '.join(f"{a}: {'raised (' + ', '.join(g['raised_items']) + ')' if g['raised'] else 'not raised'}"
                                 for a, g in v['guard'].items()),
         f"- step 1 on strength: {v['step1_strength'] or 'none named'}", '']
    d = res['drops']
    if len(d):
        L += ['| arm | stage | wirings | excluded | fall of strength by odor | strength by odor / E134\'s | '
              'fall of the population error (beside) |', '|---|---|---:|---:|---:|---:|---:|']
        for r in d.itertuples():
            L.append(f'| {r.arm} | {r.stage} | {r.n_wirings} | {r.n_excluded} | {r.strength_drop_mean:+.3f} | '
                     f'{r.strength_ratio_mean:.3f} | {r.drop_mean:+.3f} |')
        L.append('')
    if v['megamat_stage_error']:
        L += ['| group | stage error, strength (four stages) | change | three-view (beside) | change |',
              '|---|---:|---:|---:|---:|']
        f = lambda x: '—' if x is None else f'{x:+.3f}'
        g = lambda x: '—' if x is None else f'{x:.3f}'
        for k, x in v['megamat_stage_error'].items():
            L.append(f"| {k} | {g(x.get('stage_error_strength'))} | {f(x.get('stage_error_strength_change'))} | "
                     f"{g(x.get('stage_error_three_view'))} | {f(x.get('stage_error_three_view_change'))} |")
        L.append('')
    L += ['## Routes, exclusions, named', '', '```', json.dumps(dict(routes=v['routes'], excluded=v['excluded'],
                                                                     named=v['named']), indent=1, default=str), '```', '']
    path.write_text('\n'.join(L) + '\n')


def _score_tags(tags: Sequence[str], lock134: Optional[Mapping], readings: Optional[Sequence[str]],
                n_floor: Optional[int]) -> pd.DataFrame:
    """One group's `score_group` in a worker process of `score_groups`: the group's runs loaded again by tag from
    `step2_dirs` (so that no run is pickled), with the worker's own fly cache."""
    from kcq import stage_readings as R
    runs = [R.load_step3_run(t, dirs=step2_dirs()) for t in tags]
    return score_group(runs, _flies(), lock134, readings, n_floor)


def _score_pool(jobs: int):
    """The process pool of `score_groups` (forkserver: the parent has loaded olfsysm's OpenMP runtime, which a forked
    child must not inherit)."""
    import concurrent.futures as cf
    import multiprocessing as mp
    return cf.ProcessPoolExecutor(max_workers=jobs, mp_context=mp.get_context('forkserver'))


def score_groups(groups: Mapping[Tuple[str, str], list], lock134: Optional[Mapping],
                 readings: Optional[Sequence[str]] = None, n_floor: Optional[int] = None,
                 jobs: int = 1) -> List[pd.DataFrame]:
    """`score_group` of every group of `step2_groups`, with columns input_panel, arm and group put first, in the
    groups' order. With `jobs` > 1 the groups are scored in that many worker processes (`_score_tags`); the results are
    the serial ones, since every draw inside (the flies' floor, the halves) is seeded per call. Prints each group's
    end with the seconds since the start."""
    t0 = time.time()
    keys = list(groups)
    out: Dict[Tuple[str, str], pd.DataFrame] = {}

    def done(k, sc):
        if len(sc):
            sc.insert(0, 'input_panel', k[0])
            sc.insert(1, 'arm', k[1])
            sc.insert(2, 'group', f'{k[0]}_{k[1]}')
        out[k] = sc
        print(f'E135 score {k[0]} {k[1]}: {len(groups[k])} runs, {time.time() - t0:.0f} s', flush=True)

    if jobs <= 1 or len(keys) <= 1:
        flies = _flies()
        for k in keys:
            done(k, score_group(groups[k], flies, lock134, readings, n_floor))
    else:
        import concurrent.futures as cf
        with _score_pool(min(int(jobs), len(keys))) as ex:
            futs = {ex.submit(_score_tags, [r.tag for r in groups[k]], lock134, readings, n_floor): k for k in keys}
            for f in cf.as_completed(futs):
                done(futs[f], f.result())
    return [out[k] for k in keys if len(out[k])]


def part_score(readings: Optional[Sequence[str]] = None, n_floor: Optional[int] = None, jobs: int = 1) -> dict:
    """Step 2's runs (every lane's, `step2_groups`) read and scored as E134 scores its runs (`score_group`), their
    realism (`run_realism`), the APL16 runs' APL effect (`apl16_effect_table`), and the classes, the guard and the
    outcome on strength by odor (`score_tables`, at the thresholds frozen at step 2's start, `read_thresholds`; E134's
    step 2 for the APL16 somata's guard; step 1's step1_rungs.csv for the dendrites' guard). Refuses before E134's
    gate, the lock, step 1's fits and the frozen thresholds. Writes to results/.../: score_scores.csv,
    score_class_rows.csv, score_class_rows_second.csv, score_guard.csv, score_drops.csv, score_stage_error.csv,
    score_imposed_checks.csv, score_realism.csv, score_apl_effect.csv, score_verdict.json (which `dff_repeat_asked`
    reads), score_table.md, provenance_score.json. `jobs` > 1 scores the groups in that many processes
    (`score_groups`; the same results, about one group's time instead of the sum: 6 h for the 73 runs on one core)."""
    lock, fits_rec = require_step2()
    thresholds = read_thresholds()
    lock134 = require_e134_lock()
    t0 = time.time()
    groups = step2_groups()
    scores = score_groups(groups, lock134, readings, n_floor, jobs)
    sc = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    gl = [str(g) for g in lock['glomeruli']]
    realism = pd.DataFrame([run_realism(r, gl, lock134) for runs in groups.values() for r in runs])
    effect = apl16_effect_table(groups)
    stage_table = read_e134_table(E134_STAGE_TABLE)
    try:
        step2_scores = read_e134_table(E134_STEP2_SCORES)
    except FileNotFoundError:
        step2_scores = None
    d = out_dir()
    s1 = pd.read_csv(d / 'step1_rungs.csv') if (d / 'step1_rungs.csv').exists() else None
    res = score_tables(sc, stage_table, realism, runs2_table(), with_rung=step2_rung_of(fits_rec) is not None,
                       thresholds=thresholds, step2_scores=step2_scores, step1_rungs=s1,
                       highest_kept=fits_rec.get('highest_kept'))
    sc.to_csv(d / 'score_scores.csv', index=False)
    res['class_rows'].to_csv(d / 'score_class_rows.csv', index=False)
    res['class_rows_second'].to_csv(d / 'score_class_rows_second.csv', index=False)
    res['guard_items'].to_csv(d / 'score_guard.csv', index=False)
    res['drops'].to_csv(d / 'score_drops.csv', index=False)
    res['stage_error'].to_csv(d / 'score_stage_error.csv', index=False)
    res['imposed_checks'].to_csv(d / 'score_imposed_checks.csv', index=False)
    realism.to_csv(d / 'score_realism.csv', index=False)
    effect.to_csv(d / 'score_apl_effect.csv', index=False)
    v = dict(res['verdict'], time=_now(), host=_host())
    _write_json(d / SCORE_VERDICT_FILE, v)
    write_score_table(res, d / 'score_table.md')
    from kcq import stage_readings as R
    _write_json(d / 'provenance_score.json', dict(
        part='score', script=SCRIPT, host=_host(), time=_now(), seconds=round(time.time() - t0, 1), jobs=int(jobs),
        **H.code_commits(), lock_sha256=_sha256(LOCK_COPY) if LOCK_COPY.exists() else None,
        e134_lock_sha256=_sha256(LOCK_E134()) if LOCK_E134().exists() else None,
        e134_stage_table_sha256=_sha256(e134_results_dir() / E134_STAGE_TABLE),
        thresholds_sha256=_sha256(d / THRESHOLDS_FILE),
        step1_fits_sha256=_sha256(d / STEP1_FITS_FILE) if (d / STEP1_FITS_FILE).exists() else None,
        step1_rungs_sha256=_sha256(d / 'step1_rungs.csv') if s1 is not None else None,
        readings=R.reading_definitions(), readings_scored=list(readings) if readings else 'default and alternatives',
        runs={r.tag: [f for f in r.files if str(f.get('path', '')).endswith('_run.json')]
              for rs in groups.values() for r in rs}))
    print(f"E135 score: outcome {v['outcome']}, class {v['cls']} ({time.time() - t0:.0f} s)", flush=True)
    return v


def print_commands(n_controls: Optional[int] = None) -> None:
    """The lanes' command lines (LANES), as module invocations from the repository's root: step 1's real rungs on the
    first of STEP1_LANES (hal's and disco's, Tom's words), its controls (`n_controls`, default the lock's count) one
    process each on the others in turn (step 2 does not wait for them), and step 2 on every lane: hal's two 16-thread
    lanes, viki's, atlas's and disco's at 8.
    Each step-2 lane runs the whole list and skips what another has claimed or done; any lane may take `--reverse` (or
    drop it) to start from the other end."""
    pre = ('OMP_NUM_THREADS={t} OMP_THREAD_LIMIT={t} PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u -m ' + MODULE)
    lanes = [(host, t, name, rev) for host, (t, ls) in LANES.items() for name, rev in ls]
    if n_controls is None:
        try:
            n_controls = int(read_lock().get('permutations', {}).get('control', N_PERM_CONTROL))
        except PermissionError:
            n_controls = N_PERM_CONTROL
    first = next(x for x in lanes if x[2] == STEP1_LANES[0])
    others = [x for x in lanes if x[2] in STEP1_LANES[1:]]
    print("# step 1 (on hal's and disco's lanes: the real rungs on one, the controls one process each)")
    print(f'{first[0]}: {pre.format(t=first[1])} --part step1 --lane {first[2]}')
    for i in range(int(n_controls)):
        host, t, name, _ = others[i % len(others)]
        print(f'{host}: {pre.format(t=t)} --part step1_control --perm {i} --lane {name}')
    print("# step 2 (after step 1's fits; every lane runs the whole list, skipping what is claimed or done; "
          '--reverse takes it from the other end)')
    for host, t, name, rev in lanes:
        print(f"{host}: {pre.format(t=t)} --part step2 --lane {name}{' --reverse' if rev else ''}")
    print('# the score (no simulation; after the lanes)')
    print(f'{first[0]}: ./venv/bin/python -u -m {MODULE} --part score')


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('step0', 'step0_saved', 'timing', 'lock', 'step1', 'step1_control', 'step2', 'score',
                             'commands', 'gate', 'variants'))
    ap.add_argument('--n-folds', type=int, default=None,
                    help='lock: override the budget rule\'s number of odor folds (default: 17, or 6 on cost)')
    ap.add_argument('--n-controls', type=int, default=None,
                    help='lock: override the budget rule\'s number of R1 / R2 controls (default: 5, or 3 on cost)')
    ap.add_argument('--timing', default=None, help='lock: the timing file (default results/.../timing_step1.json)')
    ap.add_argument('--perm', default=None,
                    help="step1_control: the control indices, comma-separated, or 'all' (the lock's count)")
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--lane', default=None)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--readings', nargs='+', default=None,
                    help="score: the stage readings by name (default each stage's default and its alternatives)")
    ap.add_argument('--n-floor', type=int, default=None, help="score: the flies' floor's draws (default E134's)")
    ap.add_argument('--jobs', type=int, default=1,
                    help='score: the groups scored in this many processes at once (default 1; the same results)')
    a = ap.parse_args(argv)
    if a.part == 'step0':
        part_step0(n_boot=a.n_boot, seed=a.seed)
    elif a.part == 'step0_saved':
        part_step0_saved()
    elif a.part == 'timing':
        part_timing(seed=a.seed)
    elif a.part == 'lock':
        write_lock(a.n_folds, a.seed, a.timing, a.n_controls)
    elif a.part == 'step1':
        part_step1()
    elif a.part == 'step1_control':
        if a.perm is None:
            ap.error('step1_control needs --perm')
        if a.perm == 'all':
            n = int(read_lock().get('permutations', {}).get('control', N_PERM_CONTROL))
            perms = list(range(n))
        else:
            perms = [int(x) for x in a.perm.split(',')]
        part_step1_control(perms, a.lane, a.reclaim)
    elif a.part == 'step2':
        part_step2(a.lane, a.reverse, a.reclaim)
    elif a.part == 'score':
        part_score(a.readings, a.n_floor, a.jobs)
    elif a.part == 'commands':
        print_commands()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'gate':
        print(json.dumps(dict(e134_gate_status(), step1_authorization=step1_authorization()), indent=1))


if __name__ == '__main__':
    main()
