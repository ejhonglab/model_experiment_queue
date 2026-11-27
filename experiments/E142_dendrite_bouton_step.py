#!/usr/bin/env python3
"""E142, the step from the PN dendrites to the PN boutons by itself, with the recorded dendrites under the stage fit's
calyx: which odors carry the boutons' error, whether the connectomes' bouton counts make the flies' step, and (step 1)
the recorded dendrites imposed under E136's A1 calyx factors.

The entry is item E142 of experiment_priorities.md (proposed by viki 2026-11-26 05:10 from E139's next steps, E134's
step table, E135's imposed arm and a read of E138's and E140's saved entries; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-26 04:48; moved in under the idle-host rule). "The family", S, T, the halves
and the readings are those of the best family fitted to all megamat stages at once (E136): C2's breadth-2 APL-effect fit
(E67) on E72 part B's eight wirings (E83's; S real, P1, P2, M1 and T P3, P4, M2, M3), each read as one fly; the views,
one-fly units and readings are E134's (kcq/population_fit.py, kcq/stage_readings.py, data_prep/stage_population_data.py
and its lock experiments/E134_lock.json). Strength by odor has two parts, each odor's response rate (part 'rate') and
its mean response (part 'mean'), each as a multiple of the mean over the panel's 17 odors; a part's error is the root
mean square over odors of the model's (or a yardstick's) entry minus the flies' mean entry in the flies' one-fly unit,
and the strength view's error the root mean square of its two parts' errors (E136's adopted measure, 'strength').

The guard: no statistic is computed before this experiment's lock exists, is committed with no uncommitted change and
equals its copy in the shared results (`require_lock`); E134's lock guards every reading of the flies
(`stage_readings.require_lock`). Every part but 'lock' calls both first.

Step 0 (no simulation; `--part step0`, its pieces in experiments/E142_step0.py, so that step 1's code and step 0's do not
touch the same lines):
- the flies' step by odor from the dendrites to the boutons and from the boutons to the claws, on all flies and on
  E136's halves, with intervals from resampling the flies of both stages; the odors whose step is outside zero on both
  halves; the no-step yardstick's error by odor;
- every saved run's entries by odor at the dendrites, boutons and claws: each odor's share of the squared strength
  error on all flies and on halves (E138's "concentrated"), and the model's step by odor against the flies' (r over
  the 17 odors);
- the count-weighted yardsticks (the recorded dendrites with each glomerulus counted by its boutons, or its claws, in
  the hemibrain and in the male CNS's two calyces) with their null of permuted counts; beside, the same with the
  glomeruli the ORN recordings hold outside the dendrite recordings carried by the median measured gain;
- what the model's bouton reading samples: its units per glomerulus beside the connectomes' counts, and the share of
  units on the 18 glomeruli without input;
- the equal-weight reading of the boutons (each of E128's tuning groups one unit in the flies, each glomerulus's
  boutons one unit in the model), the flies and every saved run, and E136's and E139's classes read again with it and
  with the megamat stage error over the three other stages;
- the planner's probe numbers reproduced (`PROBE`), and step 0's classes ("the counts make the step", "the reading's").

Step 1 (56 runs at fixed scales on E135's step-2 runner and table) is hal's module, experiments/E142_imposed_calyx_runs.py
(`--part step1` calls its `part_step1`); its arms, A1's factors and its classes are fixed in the lock (the module's
LOCK_ARMS, `CLASSES`); its score (`--part score`) is reserved for hal. The lock's second amendment (2026-11-26 14:22,
after the first holds and before any score; AMENDMENTS) adds the main arm on A1's own PN stage and the arm with
APL->bouton at none and KC-KC at the family's (40 runs, the module's AMENDED_ARMS) and the classes' amended reading.

Judgments the entry leaves open (the narrowest reading; each is in the lock's 'narrow_readings'):
see NARROW_READINGS.

Concentrations: megamat at 10^-3 in every recording and input used. Step 1's APL16 arm involves paa, euc and -aPine,
whose model input is the ORN data within 1 log10 unit of the flies' concentrations (step 0 reads no APL16 number).

Holds (CLAUDE.md, the entry's Limits): no bouton region is given a glomerulus or an identity from position or tuning
(counts are applied to named glomeruli of the dendrite recordings only; tuning groups stay unnamed); only Pratyush's
hemibrain and male CNS (olfsysm's loaders of his tables, through E138's and E130's code); only processed data; Remy's
processed data without the bouton fly 2022-10-21 fly02 and without her older megamat flies (pyolfsysm's
data/from_remy/2024-11-12); validation2 is not read; nothing is fitted or searched; no olfsysm change and no new default
(step 1's dendrites enter by `pn_rates`, olfsysm b1e89d7, the option E135 used with Yang's OK). Nothing is published.

Outputs (<shared>/results/E142_dendrite_bouton_step/): lock.json; step 0's tables (E142_step0.py's docstring lists
them), classes_step0.json, tables.md, provenance_<part>.json.

Usage (from the repository root; step 0 on one thread under nice):
    ./venv/bin/python -u experiments/E142_dendrite_bouton_step.py --part lock      # then commit experiments/E142_lock.json
    OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 nice -n 10 ./venv/bin/python -u experiments/E142_dendrite_bouton_step.py \\
        --part step0
"""
from __future__ import annotations

import argparse
import json
import socket
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from kcq import harness as H

EXPERIMENT = 'E142_dendrite_bouton_step'
SCRIPT = 'experiments/E142_dendrite_bouton_step.py'
LOCK_NAME = 'lock.json'
LOCK_COPY = H.HERE / 'experiments' / 'E142_lock.json'
E134_LOCK_COPY = H.HERE / 'experiments' / 'E134_lock.json'
E138_LOCK_COPY = H.HERE / 'experiments' / 'E138_classes_lock.json'
E140_LOCK_COPY = H.HERE / 'experiments' / 'E140_lock.json'

PANEL = 'megamat'
N_ODORS = 17
S_WIRINGS: Tuple[str, ...] = ('real', 'P1', 'P2', 'M1')
T_WIRINGS: Tuple[str, ...] = ('P3', 'P4', 'M2', 'M3')
WIRINGS: Tuple[str, ...] = S_WIRINGS + T_WIRINGS
STRENGTH_PARTS: Tuple[str, ...] = ('strength/rate', 'strength/mean')
# the stages step 0 reads by odor (the dendrites for the model's step; the boutons and the claws the entry's), the
# stages of the megamat stage error (E136's readable stages), and the steps
ENTRY_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw')
LOSS_STAGES: Tuple[str, ...] = ('pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
STEPS: Tuple[Tuple[str, str], ...] = (('pn_dendrite', 'pn_bouton'), ('pn_bouton', 'kc_claw'))
# the fly sets: all flies, and E136's halves (hal's split_halves(n, 1, 0), the first part 'fit'; each stage split
# alone, a step's halves paired by name)
HALVES: Tuple[str, ...] = ('fit', 'held')
FLY_SETS: Tuple[str, ...] = ('all',) + HALVES
HALF_SEED = 0

# the flies' step by odor: hal's `step_interval` (each stage's flies resampled with replacement, independently)
STEP_N_BOOT = 1000
STEP_CI = 0.9
STEP_SEED = 0

# the count-weighted yardsticks: the calyces (E138's: the family's hemibrain wiring and Pratyush's male CNS, left and
# right), the counts, the null and the class's margins
CALYCES: Tuple[str, ...] = ('hemibrain-first', 'malecns-left', 'malecns-right')
MALE_CALYCES: Tuple[str, ...] = CALYCES[1:]
CONNECTOMES: Dict[str, Tuple[str, ...]] = {'hemibrain': ('hemibrain-first',), 'male_cns': MALE_CALYCES}
COUNT_OF: Dict[str, str] = {'pn_bouton': 'n_boutons', 'kc_claw': 'n_claws'}
N_PERM = 1000
PERM_SEED = 0
COUNTS_MARGIN = 0.10           # below the no-step yardstick by 0.10 or more on both halves
NULL_PERCENTILE = 5.0          # and below the 5th percentile of the permuted counts
# beside: the glomeruli the ORN recordings hold outside the dendrite recordings, carried by the median measured gain
# (E135's step 0: the slope through the origin of a glomerulus's dendrite response on its own ORN-terminal response
# across odors, fly means, in trial-noise SD units, all flies; the median over the 28 shared glomeruli)
E135_GAIN_UNITS = 'noise'
E135_GAIN_FLIES = 'all'

# the equal-weight reading (E128's tuning groups, E119's A4): responsive ROIs (any odor above the criterion) grouped by
# their tuning on E119's 12 odors (average linkage on 1 - Pearson r, cut at 0.3) from the mean of two repeats, read on
# the third; the remainder one more unit
GROUP_CRITERION = 2.0
GROUP_CUT = 0.3

# E138's "concentrated": the two largest odors on all flies carry more than half of the squared strength error on both
# halves (E138's CONCENTRATED_TOP, CONCENTRATED_SHARE; read on the row 'common', the wirings' mean views)
CONCENTRATED_TOP = 2
CONCENTRATED_SHARE = 0.5

# the planner's probe numbers (Section 3 of its answer; E142's entry), reproduced in step 0 before anything rests on
# them: the lead table (E140's models_entries.csv / flies_entries.csv, reading curve_noise, APL on, column 'rate'; the
# family and the imposed dendrites the mean over eight wirings, A1 over its four T wirings; the flies' 'mean' row of the
# wild-type cohort), E139's stage table, the room's arithmetic, E134's and E135's bouton numbers
PROBE_TOL = 0.005      # the numbers are given to two decimals (three for the stage table)
PROBE = dict(
    lead_rate={
        # odor: (flies' dendrites, flies' boutons, the family held, A1, the dendrites imposed), boutons for the models
        'pa': (1.07, 1.47, 0.84, 0.94, 1.00), 'IaA': (1.09, 1.48, 1.07, 1.18, 0.98),
        'ep': (1.22, 1.55, 1.10, 1.11, 1.21), 'eb': (1.19, 1.51, 1.07, 1.10, 1.28),
        '1-6ol': (1.62, 1.10, 1.12, 1.24, 1.73), 'benz': (1.31, 0.90, 1.00, 1.00, 1.35),
        'ms': (0.67, 0.48, 1.05, 0.89, 0.58), 'va': (0.82, 0.39, 0.82, 0.79, 0.74)},
    lead_sd=(0.33, 0.37, 0.13, 0.18, 0.38),
    lead_r=(None, None, 0.44, 0.63, 0.62),
    lead_columns=('flies_dendrites', 'flies_boutons', 'family_held', 'A1', 'e135_imposed'),
    stage_table={
        # stage: (half to half, the family, A1), strength on T against the held-out half (E139's table, E136's scores)
        'pn_dendrite': (0.798, 1.012, 0.858), 'pn_bouton': (0.817, 1.194, 1.005), 'kc_claw': (0.830, 0.972, 0.905),
        'kc_soma': (0.970, 1.003, 0.961), 'stage': (0.856, 1.049, 0.934)},
    room=dict(boutons_at_half=0.887, boutons_and_claws_at_half=0.868),
    e134_no_step=dict(family_boutons=1.06, no_step=1.05),
    e135_imposed=dict(imposed_boutons=1.14, family_boutons=1.07),
)

# the saved runs (the entry's "Definitions"): group -> source, wirings, what; sources 'E134:<tag pattern>' (E134's
# step-3 runs), 'E136:<pattern>' (E136's reads under the adopted measure), 'E139:<pattern>' (E139's reads, a registered
# variant of E136's), 'E135:<pattern>' (E135's step-2 runs)
ModelGroup = namedtuple('ModelGroup', 'name source wirings what')
MODEL_GROUPS: Tuple[ModelGroup, ...] = (
    ModelGroup('family_held', 'E134:megamat_held_{w}', WIRINGS,
               "the family held at the megamat flies' 0.0915 (E134's step 3), the unrefitted family of every class"),
    *(ModelGroup(a, f'E136:{a}_megamat_{{w}}', T_WIRINGS,
                 f"E136's arm {a}, its megamat reads on the T wirings (its S wirings are its fit's evaluations)")
      for a in ('A0', 'A1', 'A4')),
    *(ModelGroup(v, f'E139:{v}_megamat_{{w}}', WIRINGS, f"E139's {v} ({d})")
      for v, d in (('F_inhadd314', "the family's calyx factors, pn_inhadd 314"),
                   ('F_inhadd1000', "the family's calyx factors, pn_inhadd 1000"),
                   ('B_inhadd314', 'APL->bouton x0.155, KC-KC x1, pn_inhadd 314'),
                   ('B_inhadd1000', 'APL->bouton x0.155, KC-KC x1, pn_inhadd 1000'),
                   ('K_inhadd1000', 'KC-KC x3.24, APL->bouton x1, pn_inhadd 1000'))),
    ModelGroup('e135_imposed', 'E135:megamat_imposed_{w}', WIRINGS,
               "E135's recorded dendrites imposed at the family's calyx (step 2, the table in noise units)"),
    ModelGroup('e135_imposed_glom_perm', 'E135:megamat_imposed_glom_perm_{w}', WIRINGS,
               "E135's control: the imposed table with its glomeruli permuted"),
    ModelGroup('e135_imposed_odor_perm', 'E135:megamat_imposed_odor_perm_{w}', WIRINGS,
               "E135's control: the imposed table with its odors permuted"),
)
# the arms whose classes are read again with the equal-weight reading (E136's classes.json under measure_strength, and
# E139's classes.json), and the arm of "the reading's"
E136_EXPERIMENT = 'E136_stage_population_fit'
E136_MEASURE = 'strength'
E139_EXPERIMENT = 'E139_normalization_ladder'
E135_EXPERIMENT = 'E135_antennal_lobe_dendrites'
E134_EXPERIMENT = 'E134_population_baseline'
E138_EXPERIMENT = 'E138_strength_by_odor_entries'
E140_EXPERIMENT = 'E140_aa_claw_soma_step'
E128_EXPERIMENT = 'E128_bouton_groups'
E119_EXPERIMENT = 'E119_stages_per_roi'
READING_ARM = 'A1'
FAMILY = 'family'

# step 1 (hal's module experiments/E142_imposed_calyx_runs.py, queue d6b58ff): its arms (`STEP1_ARMS` there: panel,
# name, E135's table, and each calyx factor 'A1' (A1's), 'family' (x1, the wiring's own fit's scale) or 'none' (0)), tags
# '<panel>_<arm>_<wiring>', 7 arms x 8 wirings; A1's factors E136's classed point of A1 under strength (fit_A1.json's
# 'classed_z': APL->bouton x0.155, KC-KC x3.238), the lock records them as read when it was written
STEP1_MODULE = 'experiments/E142_imposed_calyx_runs.py'
STEP1_APL_ARMS: Tuple[str, ...] = ('megamat_imposed_a1', 'megamat_imposed_aplpn', 'megamat_imposed_kckc_aplpn_none')
STEP1_CARRY_ARM = 'megamat_imposed_a1'
STEP1_CONTROLS: Tuple[str, ...] = ('megamat_imposed_a1_glom_perm', 'megamat_imposed_a1_odor_perm')
N_STEP1_RUNS = 56
# the lock's amendments (each before any statistic of E142 was computed)
AMENDMENTS: Tuple[str, ...] = (
    "2026-11-26 05:23, viki (before any statistic): step 1's arms, names and readings taken from hal's runner "
    "(E142_imposed_calyx_runs.py, d6b58ff) in place of the first lock's (059753e), which had named the arm with "
    "APL->bouton at none 'imposed_kckc_aplpn0' and the APL16 arm 'imposed_a1_apl16' and read the factors as the entry's "
    "rounded x0.155 and x3.24; step 0's definitions unchanged",
    "2026-11-26 14:22, hal (after the first holds were read and before any step-1 run was scored; the entry's *Step 1 "
    "amended* of 13:59, planned by a Fable subagent, claude-fable-5-2, spawned from hal 2026-11-26 13:48): on the "
    "lock's table (the dendrites ranked onto the family's own PN rates, pn_inhadd 31.4) KC-KC x3.24 cannot be held at "
    "megamat's operating point (megamat_imposed_a1's threshold at the hold's bound, x0.1 of the wiring's fit, on 8 of "
    "8 wirings; every megamat arm with KC-KC x3.24 on that table has the same APL-off run), so the lock's arms on the "
    "family's PN stage are kept as the record: megamat_imposed_a1 'no class: a hold fails'; megamat_imposed_kckc and "
    "megamat_imposed_kckc_aplpn_none not run (placeholder claims, lane 'stopped'), 'no class: a hold fails'; the two "
    "permuted tables read once (M3), then stopped; megamat_imposed_aplpn and apl16_imposed_a1 run as they are. Added, "
    "40 runs (step 1 is 66: 26 of the lock's 56 and 40 new): the main arm on A1's own PN stage, A1's classed point "
    "whole (E136's pn_norm_kws(fitted_z('A1')): pn_inhadd 1017.3 and pn_inhsc 598.85; APL->bouton x0.155 and KC-KC "
    "x3.24), its table E135's imposed table built on A1's own PN stage (the dendrites' arrangement ranked onto A1's "
    "PN rates on the 28 glomeruli, every other entry A1's), as megamat_a1stage_imposed, its two permuted controls "
    "megamat_a1stage_imposed_glom_perm and megamat_a1stage_imposed_odor_perm (E135's seeds) and apl16_a1stage_imposed "
    "(its 12 megamat odors set), eight wirings each, by pn_rates, each held as A1's own reads were (E136's hold_at "
    "from the wiring's held megamat point, read_targets' counts, the same brackets); and 'APL->bouton at none' "
    "replaced by megamat_imposed_aplpn_none (the family's PN stage, KC-KC at the family's, APL->bouton 0, the lock's "
    "hold), eight runs. Not run: each factor alone on A1's stage. Refused: a wider threshold bracket, and the APL-off "
    "condition at KC-KC x1. The classes as amended (thresholds, spreads and the bar as locked): 'the dendrites carry "
    "under the calyx that holds' is read on megamat_a1stage_imposed against A1's own megamat reads on T, with the "
    "A1-stage permuted controls and the A1-stage APL16 arm's bands; 'APL at the boutons hides the step' on the "
    "family's PN stage, in megamat_imposed_aplpn or megamat_imposed_aplpn_none against E135's imposed arm, an arm "
    "whose hold fails left out and named; beside, classing nothing: the amended main arm's boutons against the "
    "no-step yardstick, and E143's interval over resampled flies on its change from A1 with the share of resamples on "
    "which the class returns. Outcomes (d) and (e) are read on the amended main arm; (g), new, 'not readable': the "
    "amended main arm's hold fails or it leaves a realism band, and nothing further is tried in this entry. Code: "
    "experiments/E142_imposed_calyx_runs.py (the arms' pn_stage, LOCK_ARMS and AMENDED_ARMS, classify_step1, "
    "step1_outcome)",
)
# step 1's points as amended, recorded in the lock beside the amendments (2026-11-26, after the score, changing no run
# or class): the amendment of 14:22 added five arms (AMENDED_ARMS) to the lock's seven; step1's n_runs stays the lock's
N_STEP1_RUNS_AS_AMENDED = 96
AMENDMENT_COUNTS = dict(
    n_runs_as_amended=N_STEP1_RUNS_AS_AMENDED, n_written=59, n_stopped_by_claim=37,
    note="the amendment of 2026-11-26 14:22 added 5 arms to the lock's 7, so step 1 planned 12 arms x 8 wirings = 96 "
         "points: 59 written, 37 stopped by claim (step1_holds.csv); step1.n_runs (56) is the lock's own count and its "
         "'step 1 is 66' the runs the amendment meant to run (7 of them later stopped too); recorded 2026-11-26 after "
         "the score, changing no run or class")

# the classes, fixed before any statistic (the entry's, made exact; NARROW_READINGS say where a reading was chosen)
CLASSES = {
    'counts_make_step': dict(
        name='the counts make the step', step=0,
        rule='for a connectome, the count-weighted yardstick (the recorded dendrites, all seven flies, each glomerulus '
             'counted by its boutons in the connectome, as the prediction for the bouton flies) is below the no-step '
             "yardstick (the same with every glomerulus once) by COUNTS_MARGIN or more in the strength view's error, "
             "against each half of the bouton flies (E136's halves, each half's mean in all five flies' one-fly units), "
             'and below the NULL_PERCENTILE-th percentile of the same error with the counts permuted among the '
             'glomeruli (N_PERM permutations, PERM_SEED), against all five bouton flies',
        margin=COUNTS_MARGIN, null_percentile=NULL_PERCENTILE, n_perm=N_PERM, perm_seed=PERM_SEED,
        calyces=list(CALYCES), connectomes={k: list(v) for k, v in CONNECTOMES.items()},
        connectome_rule="the hemibrain holds when 'hemibrain-first' does; the male CNS when both of its calyces do "
                        '(each calyx given)',
        beside=["the claws with claw counts (against the claw flies, both no-step yardsticks)",
                'the glomeruli outside the dendrite recordings carried by the median measured gain']),
    'apl_hides_step': dict(
        name='APL at the boutons hides the step', step=1,
        rule="the imposed dendrites' boutons read below the no-step yardstick on 6 of 8 wirings in an arm with "
             "APL->bouton at x0.155 or none (megamat_imposed_a1, megamat_imposed_aplpn, "
             "megamat_imposed_kckc_aplpn_none: any of them), and not with it at the family's (E135's saved imposed arm, "
             're-scored the same way)',
        min_wirings=6, n_wirings=8, arms=list(STEP1_APL_ARMS), reference_arm='e135_imposed'),
    'dendrites_carry': dict(
        name='the dendrites carry under the calyx that holds', step=1,
        rule="with A1's factors (megamat_imposed_a1), the root mean square of the boutons', claws' and somata's "
             "strength errors is below A1's own by more than twice its standard error over T, no stage is worse than "
             "A1's by more than its spread, and every stage is below both permuted controls "
             '(megamat_imposed_a1_glom_perm, megamat_imposed_a1_odor_perm)',
        se_multiple=2.0, arm=STEP1_CARRY_ARM, reference=READING_ARM, controls=list(STEP1_CONTROLS)),
    'readings': dict(
        name="the reading's", step=0,
        rule="under the equal-weight reading of the boutons, A1's change at the boutons from the family (the mean over "
             'T of A1 minus the unrefitted family on the same wiring, against the held-out half) is inside the stage\'s '
             "spread (the SD over the family's eight wirings of its equal-weight bouton error against the held-out "
             'half) or of the other sign (above 0): change >= -spread',
        arm=READING_ARM),
    'no_class': dict(name='no class', rule='a run outside a realism band or whose holds fail gets no class'),
}
OUTCOMES = {
    'a': "the counts make the step, and the model's reading does not sample them: the boutons' error is in part how "
         "the model's boutons are read; the reading is corrected like for like and the saved reads scored again",
    'b': "the counts do not, and no arm of step 1 moves the boutons: the flies' step is not in the dendrites, the counts "
         "or APL's bouton weights; whether the boutons keep an equal weight in the megamat stage error goes to Tom and "
         'Yang',
    'c': "APL at the boutons hides the step: APL's bouton weights by glomerulus are the lead",
    'd': "the dendrites carry under A1's calyx: the ORN->PN route is open",
    'e': 'they do not carry: no term at the ORN->PN step held to the 28 recorded glomeruli can lower the stages after '
         'the dendrites',
    'f': "the reading's: A1's gain at the boutons depends on how regions are counted; to Tom and Yang with the megamat "
         'stage error under both readings',
}

NARROW_READINGS: Tuple[str, ...] = (
    "the saved runs: 'E139's F and B at 314 and 1000 and K at 1000' read as F and B each at pn_inhadd 314 and 1000, and "
    "K at 1000 (five variants; F at 100 and 3141 and E139's step-3 fit A0_inhadd314 are not read)",
    "E136's arms A0, A1, A4 are read on their four T wirings (their megamat reads; the S wirings are fit evaluations), "
    'the family held, E139\'s variants and E135\'s three imposed arms on all eight; "common" is the wirings\' mean views',
    "E135's imposed arm is its noise-unit table (megamat_imposed_<w>), not the dF/F repeat; its two permuted controls "
    'are megamat_imposed_glom_perm_<w> and megamat_imposed_odor_perm_<w>',
    "every run is read in E134's default readings (the lock's), APL on only (no megamat flies with APL output blocked), "
    "against each fly set as E136's Flies reads a half: the curve, noise and criterion from that set's flies, the "
    "reference its mean views in all the stage's flies' one-fly units; the seed index the wiring's position in WIRINGS",
    "the claws' flies are E136's: the four wild-type claw flies (the same-fly recording is not read), halves 2 / 2",
    "the flies' step on halves pairs each stage's half by name (the dendrites' fit half with the boutons' fit half, "
    "and so on), each stage split alone by E136's split_names",
    "a step is outside zero on a set when its 90 % interval (hal's step_interval, 1000 resamples, seed 0, both stages' "
    'flies resampled independently, uncorrected) excludes 0; "outside zero on both halves" when it is on each half '
    'with the same sign',
    "the no-step yardstick and the count-weighted yardsticks predict from all seven dendrite flies (their mean views) "
    "against all bouton flies and against each bouton half (in all five flies' units), as E134's no_step_score does",
    'the count-weighted yardstick weighs each glomerulus of each dendrite fly by its count (equivalent to repeating its '
    'row that many times): the rate part the count-weighted fraction of glomeruli responding, the mean part the '
    "count-weighted mean response, each relative to its mean over the 17 odors; the dendrite recordings' 29 glomeruli "
    "(E135's 28 and VA4) as recorded, every one of which each calyx has",
    'the null permutes the count vector among the 29 dendrite glomeruli (one set of N_PERM permutations from seed 0, '
    'the same for every calyx and stage); the percentile is of the error with all bouton (or claw) flies',
    "the male CNS's two calyces are read each as a connectome's reading; the male CNS 'makes the step' only when both "
    'calyces do',
    "beside: the glomeruli the ORN recordings hold outside the 29 dendrite glomeruli (DC4, DL2d, DL2v, DP1l, VC3, VC4, "
    "VL1, VM5d, VM5v; VA4 keeps its recorded dendrites) are carried as the median over E135's 28 of its measured slope "
    "(units 'noise', all flies) times the ORN flies' mean response (trial-noise SD units, over the ORN flies that have "
    "it), the same rows added to every dendrite fly, responding above the dendrites' criterion 2.0, each counted by its "
    'count; the null then permutes among the 38',
    "the claws' yardsticks: the dendrites weighted by claw counts against the claw flies, beside the dendrites "
    "unweighted (no-step from the dendrites) and the bouton flies (E134's no-step for the step to the claws); no class",
    "entries by odor are E138's odor_entries (each odor's share of the squared strength view error, 'view_share'); "
    "E138's concentrated on the row 'common' (for a group of T wirings, their mean views)",
    "the model's step by odor: per part, the model's entry downstream minus upstream on the row 'common', against the "
    "flies' mean entries' difference; Pearson r over the 17 odors (Spearman beside), on all flies",
    "the equal-weight reading of the flies: per bouton fly and held-out repeat, E128's fly_groups on the mean of the "
    "other two repeats restricted to E119's 12 odors (E90's ODORS; criterion 2, cut 0.3), each group one unit valued "
    "at its ROIs' mean on the held-out repeat on the 17 odors, the remainder (ROIs responsive to no odor) one more unit, "
    'a unit responding above 2.0; the fly\'s strength parts the mean over its three held-out repeats\' parts; the '
    "groups' sizes on the trial means checked against E119's a4_groups.csv",
    "the equal-weight reading of the model: the default bouton reading (pool 1, the 28 shared glomeruli, silent "
    "boutons dropped, curve and noise from the fly set), then each glomerulus's read boutons averaged into one unit, "
    'responding above 2.0; the glomerulus labels are the run\'s own (meta bouton_glomerulus), never the flies\'',
    "the classes read again with the equal-weight reading keep from the original records what does not involve the "
    "bouton strength error (the held-out panels, the all-views guard, realism, bounds, the control fit's fall); the "
    "bouton stage's spread and the megamat stage error's are recomputed from the family's eight wirings under the "
    "reading; 'over the three other stages' drops the boutons from the stage error and the stage comparisons",
    "the model's bouton sampling is read on each wiring's saved runs of the family held (the units per glomerulus, APL "
    "on, after the silent drop); the share of units on the 18 glomeruli without input counts every bouton of the "
    "wiring, and the default reading's share (subset 'shared': 0 by construction) is given beside",
    "step 1 (hal's runner, agreed): A1's calyx factors are E136's classed point of A1 under strength (fit_A1.json's "
    "'classed_z', x0.155 and x3.238; the entry's x3.24 is its rounding), because A1's own reads, against which 'the "
    "dendrites carry under the calyx that holds' is classed, are at that point; on each wiring's own fit's wAPLPN_scale and "
    "wKCKC_scale; A1's weaker antennal-lobe normalization is not taken (the glomeruli and odors not imposed, and the "
    "distribution the dendrites are ranked onto, are the family's, as in E135), so the change from E135's runs is the "
    'calyx alone; "APL->bouton at none" is wAPLPN_scale 0; the hold is E135\'s (E134\'s held megamat targets, the '
    "APL16 arm at each wiring's own fit's counts); the permuted controls reuse E135's permutations (per-wiring seeds); "
    "no run before results/E142_dendrite_bouton_step/step1_authorized.json records Yang's word on pn_rates",
    "step 1's 'below the no-step yardstick on 6 of 8 wirings': each wiring's bouton strength error against the "
    "held-out half of the bouton flies, below the no-step yardstick against the same half",
    "step 1's 'twice its standard error over T': E136's compare (the larger of half the spread and the SD of the "
    "per-wiring changes over sqrt(4)), the spread the SD over the family's eight wirings of the downstream error (the "
    "boutons', claws' and somata's root mean square); 'no stage worse than A1's by more than its spread' E136's "
    "spreads per stage; 'every stage below both permuted controls' the mean over T of each stage's error",
)


def out_dir() -> Path:
    """<shared>/results/E142_dendrite_bouton_step/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir() -> Path:
    """outputs/E142_dendrite_bouton_step/ (local, created)."""
    return H.outputs_dir(EXPERIMENT)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _e140():
    from experiments import E140_aa_claw_soma_step as E140
    return E140


def results_file(experiment: str, *parts: str) -> Path:
    """A file of another experiment's results in the shared directory."""
    return H.shared_dir().joinpath('results', experiment, *parts)


# ---------------------------------------------------------------------------
# the lock

def lock_inputs() -> dict:
    """The inputs the lock records with sha256 (E137's and E140's form): the locks and data files every statistic rests
    on ('data', checked when read: a change stops the run) and the saved results step 0 reads or reproduces ('results',
    recorded: a later change is reported, not refused)."""
    E140 = _e140()
    data = dict(e134_lock=E140.file_entry(E134_LOCK_COPY), e138_classes_lock=E140.file_entry(E138_LOCK_COPY),
                e135_lock=E140.file_entry(results_file(E135_EXPERIMENT, 'lock.json')),
                e135_measured_step=E140.file_entry(results_file(E135_EXPERIMENT, 'step0_measured_step.csv')))
    res = {}
    for k, (exp, parts) in dict(
            e138_feeders=(E138_EXPERIMENT, ('feeders.csv',)), e138_stages=(E138_EXPERIMENT, ('stages.csv',)),
            e140_models_entries=(E140_EXPERIMENT, ('models_entries.csv',)),
            e140_flies_entries=(E140_EXPERIMENT, ('flies_entries.csv',)),
            e119_a4_groups=(E119_EXPERIMENT, ('a4_groups.csv',)), e119_a4_summary=(E119_EXPERIMENT, ('a4_summary.csv',)),
            e128_a_repeat=(E128_EXPERIMENT, ('a_repeat.csv',)),
            e134_flies_no_step=(E134_EXPERIMENT, ('flies_no_step.csv',)),
            e134_score_scores=(E134_EXPERIMENT, ('score_scores.csv',)),
            e135_score_scores=(E135_EXPERIMENT, ('score_scores.csv',)),
            e135_score_stage_error=(E135_EXPERIMENT, ('score_stage_error.csv',)),
            e136_classes=(E136_EXPERIMENT, (f'measure_{E136_MEASURE}', 'classes.json')),
            e136_score_megamat=(E136_EXPERIMENT, (f'measure_{E136_MEASURE}', 'score_megamat.csv')),
            e136_score_panels=(E136_EXPERIMENT, (f'measure_{E136_MEASURE}', 'score_panels.csv')),
            e136_score_apl16=(E136_EXPERIMENT, (f'measure_{E136_MEASURE}', 'score_apl16.csv')),
            e139_classes=(E139_EXPERIMENT, ('classes.json',)),
            e139_score_megamat=(E139_EXPERIMENT, ('score_megamat.csv',))).items():
        p = results_file(exp, *parts)
        res[k] = E140.file_entry(p) if p.exists() else dict(path=str(p), missing=True)
    return dict(data=data, results=res)


def recordings() -> Dict[str, dict]:
    """The recordings read, by stage (viki's populations; names only, no value is read): the ORNs (for the carried
    glomeruli), the PN dendrites, the PN boutons, the KC claws (the four wild-type flies, E136's) and the KC somata (for
    the classes read again, through E136's saved scores); E136's halves; the exclusions."""
    from data_prep import stage_population_data as D
    from experiments import E136_stage_population_fit as E136
    out = {}
    for stage in ('orn',) + LOSS_STAGES:
        p = D.load_population(stage, PANEL, 'wild_type')
        out[stage] = dict(stage=stage, panel=PANEL, cohort='wild_type', flies=list(p.flies), criterion=float(p.criterion),
                          units=p.units, named_units=bool(p.named_units), odors=list(p.odors),
                          glomeruli=(sorted(set().union(*[set(map(str, p.tables[f].index)) for f in p.flies]))
                                     if p.named_units else None),
                          halves=(E136.split_names(p, HALF_SEED) if stage != 'orn' else None))
    if any(f.startswith('2022-10-21__fly02') for f in out['pn_bouton']['flies']):
        raise RuntimeError('the excluded bouton fly is in the bouton population')
    out['excluded'] = dict(
        flies=[dict(stage='pn_bouton', fly=x['fly'], reason=x['reason']) for x in D.EXCLUDED_FLIES],
        older_megamat="Remy's older megamat flies (pyolfsysm's data/from_remy/2024-11-12): never read",
        same_fly="the same-fly claw and soma recording (2023-04-11 fly 3, megamat1_redo): not read (E136's flies)",
        validation2='not read')
    return out


def glomerulus_sets() -> dict:
    """The named glomeruli the yardsticks use: E135's 28 (E138's lock), the dendrite recordings' union, the ORN
    recordings' union, the carried glomeruli (ORN-recorded, not dendrite-recorded), the model's 36 input glomeruli and
    E90's shared 28 (the bouton reading's subset)."""
    from experiments import E90_model_stages as E90
    e138 = json.loads(E138_LOCK_COPY.read_text())
    recs = recordings()
    d28 = list(e138['classes']['not_dendrite28']['dendrite28'])
    den = list(recs['pn_dendrite']['glomeruli'])
    orn = list(recs['orn']['glomeruli'])
    shared = [str(g) for g in E90.shared_glomeruli()]
    if sorted(shared) != sorted(d28):
        raise RuntimeError("E90's shared glomeruli differ from E138's lock's dendrite28")
    return dict(dendrite28=d28, dendrites=den, dendrites_outside28=[g for g in den if g not in d28], orn=orn,
                carried=[g for g in orn if g not in den], model36=list(e138['universe']['glomeruli']),
                bouton_reading_subset=shared)


def lock_contents() -> dict:
    """The lock: the inputs (with sha256), the recordings and halves, the glomerulus sets, the definitions, the saved
    runs and their readings, the classes (step 0's and step 1's) and outcomes, step 1's arms, the narrow readings, the
    probe numbers, and the holds. No statistic is read."""
    from data_prep import stage_population_data as D
    from kcq import stage_readings as R
    definitions = dict(
        strength="kcq/population_fit.py's strength view on the 17 megamat odors (population_views): part 'rate', each "
                 "odor's response rate over the mean over the 17 odors; part 'mean', each odor's mean response likewise; "
                 "a part's error the RMS over odors of the entry minus the flies' mean entry over the part's one-fly "
                 "unit; the view's error the RMS of the two parts' errors (E136's measure 'strength')",
        step_by_odor=dict(rule="an odor's entry in a strength part at the downstream stage minus its entry upstream "
                               "(hal's step_change of the flies' mean parts); interval from hal's step_interval",
                          n_boot=STEP_N_BOOT, ci=STEP_CI, seed=STEP_SEED, draw='sets', corrected=False,
                          steps=[list(s) for s in STEPS], sets=list(FLY_SETS),
                          outside_zero='lo > 0 or hi < 0; on both halves: on each half, with the same sign'),
        no_step="hal's no_step_score: the upstream flies' (all seven dendrite flies') mean views as the prediction for "
                "the downstream flies (against all of them, or a half's mean in all the stage's flies' units); its "
                "entries by odor E138's odor_entries of the upstream flies against the downstream reference (row "
                "'common')",
        count_weighted=dict(rule="each dendrite fly's glomeruli weighted by their counts in a calyx (COUNT_OF: boutons "
                                 "for the boutons, claws for the claws): rate the weighted fraction responding, mean the "
                                 'weighted mean response (over the glomeruli the fly recorded for the odor), each '
                                 "relative to its mean over the 17 odors; the flies' mean views scored as the no-step",
                            calyces=list(CALYCES), counts=COUNT_OF, source="E138's all_counts (E130's "
                                                                           "step0_wiring_only, olfsysm's loaders of "
                                                                           "Pratyush's tables), checked against "
                                                                           "E138's feeders.csv",
                            n_perm=N_PERM, perm_seed=PERM_SEED, null_percentile=NULL_PERCENTILE,
                            carried=dict(gain_file='E135 step0_measured_step.csv', units=E135_GAIN_UNITS,
                                         flies=E135_GAIN_FLIES, statistic='median slope over the 28')),
        entries=dict(rule="E138's odor_entries (per row: each wiring, and 'common' the wirings' mean views)",
                     concentrated=dict(top=CONCENTRATED_TOP, share=CONCENTRATED_SHARE, row='common',
                                       rule="E138's concentrated: the two largest odors by view_share on all flies "
                                            'carry more than half of the squared strength error on both halves')),
        model_step="per part, the model's 'common' entry downstream minus upstream, against the flies' mean entries' "
                   'difference; Pearson r over the 17 odors, Spearman beside; on all flies',
        bouton_sampling="per wiring (the family held's runs, APL on): the wiring's boutons per glomerulus (meta "
                        "bouton_glomerulus), the default reading's units per glomerulus (subset 'shared', pool 1, after "
                        "the silent drop), beside each calyx's counts; the share of the wiring's boutons on the "
                        "glomeruli without model input, and of the read units",
        equal_weight=dict(criterion=GROUP_CRITERION, cut=GROUP_CUT, group_odors="E119's ODORS (E90's 12)",
                          flies="E128's fly_groups on two repeats, read on the third, every group one unit and the "
                                "remainder one more; the fly's parts the mean over its three held-out repeats",
                          model="each glomerulus's read boutons one unit (their mean), responding above "
                                f'{GROUP_CRITERION:g}'),
        stage_error="the root mean square over LOSS_STAGES of each stage's strength error (E136's 'stage'); "
                    "'downstream' over the boutons, claws and somata; 'without boutons' over the dendrites, claws and "
                    'somata')
    model = dict(groups=[g._asdict() for g in MODEL_GROUPS], readings=R.reading_definitions(),
                 reading_rule="E134's default reading of each stage (stage_readings' DEFAULT_READINGS), APL on, "
                              "against each fly set as E136's Flies reads it; seed index the wiring's position in "
                              'WIRINGS', fly_sets=list(FLY_SETS), entry_stages=list(ENTRY_STAGES))
    from experiments import E142_imposed_calyx_runs as S1
    a1 = S1.a1_factors()
    # the lock's seven arms as written (the amendment of 2026-11-26 14:22 added S1.AMENDED_ARMS, in AMENDMENTS)
    step1 = dict(module=STEP1_MODULE, arms=[{k: v for k, v in a._asdict().items() if k != 'pn_stage'}
                                            for a in S1.LOCK_ARMS],
                 tags=[S1.tag_of(a, w) for a, w in S1.planned_step1(arms=S1.LOCK_ARMS)], n_runs=N_STEP1_RUNS,
                 a1_factors=a1,
                 runner="E135's step-2 runner and table (pn_rates, olfsysm b1e89d7) with the two calyx scales changed "
                        "before the hold, both response rates re-held at E134's held megamat targets (the APL16 arm at "
                        "the wiring's own fit's counts)",
                 score="E136's reading of the boutons, claws and somata on T against the held-out half; E135's saved "
                       'runs at the family\'s calyx re-scored the same way, with A1 beside',
                 realism="the megamat response rate beside the flies' 0.0915, the four bands on the APL16 arm, the "
                         'largest spike counts and bouton rates, the entries clipped, and where APL\'s inhibition is '
                         'made',
                 guard="results/E142_dendrite_bouton_step/step1_authorized.json (Yang's word on pn_rates)",
                 status="hal's code (E142_imposed_calyx_runs.py); not in step 0")
    if len(step1['tags']) != N_STEP1_RUNS:
        raise RuntimeError(f"step 1 plans {len(step1['tags'])} runs, not {N_STEP1_RUNS}")
    holds = ['no simulation, no model run, no olfsysm change, no new default, nothing fitted or searched in step 0',
             'no bouton region given a glomerulus or an identity from position or tuning; tuning groups unnamed',
             "only Pratyush's hemibrain and male CNS", 'only processed data',
             "Remy's processed data without the bouton fly 2022-10-21 fly02 and without her older megamat flies",
             'validation2 not read']
    return dict(experiment=EXPERIMENT, script=SCRIPT, part='lock', written=_now(), host=_host(),
                commits=H.code_commits(), panel=PANEL, odors=list(D.panel_odors(PANEL)),
                inputs=lock_inputs(), recordings=recordings(), glomeruli=glomerulus_sets(), definitions=definitions,
                model=model, classes=CLASSES, outcomes=OUTCOMES, step1=step1, narrow_readings=list(NARROW_READINGS),
                probe=PROBE, probe_tol=PROBE_TOL, holds=holds, amendments=list(AMENDMENTS),
                amendment_counts=dict(AMENDMENT_COUNTS),
                note='written before any statistic of E142 was computed; every part but lock refuses until this file '
                     'is committed (require_lock)')


def part_lock() -> dict:
    """results/E142_dendrite_bouton_step/lock.json and its committed copy experiments/E142_lock.json (commit it before
    any other part); refuses to overwrite a committed lock."""
    E140 = _e140()
    if LOCK_COPY.exists() and E140._committed(LOCK_COPY):
        raise PermissionError(f'{LOCK_COPY.name} is committed already; the lock is fixed')
    lock = lock_contents()
    text = json.dumps(lock, indent=1, default=E140._json_default) + '\n'
    (out_dir() / LOCK_NAME).write_text(text)
    LOCK_COPY.write_text(text)
    print(f'wrote {out_dir() / LOCK_NAME} and {LOCK_COPY} (commit the copy before any other part)', flush=True)
    return lock


def require_lock(copy: Optional[Path] = None, shared: Optional[Path] = None) -> dict:
    """The committed lock: experiments/E142_lock.json exists, is tracked with no uncommitted change, and has the same
    bytes as results/.../lock.json. PermissionError otherwise (no statistic is read before)."""
    E140 = _e140()
    copy = LOCK_COPY if copy is None else Path(copy)
    if not copy.exists():
        raise PermissionError(f'{copy.name}: no lock yet (--part lock, then commit it)')
    if not E140._committed(copy):
        raise PermissionError(f'{copy.name} is not committed (or has uncommitted changes)')
    shared = (out_dir() / LOCK_NAME) if shared is None else Path(shared)
    if not shared.exists() or shared.read_bytes() != copy.read_bytes():
        raise PermissionError(f'{shared} differs from the committed {copy.name}')
    lock = json.loads(copy.read_text())
    lock['sha256'] = E140._sha256(copy)
    return lock


def require_all() -> Tuple[dict, dict]:
    """(this experiment's lock, E134's lock), both required before any statistic; the lock's data files checked."""
    from kcq import stage_readings as R
    lock = require_lock()
    E140 = _e140()
    for k, e in lock['inputs']['data'].items():
        E140._check_file(e)
    return lock, R.require_lock()


# ---------------------------------------------------------------------------
# step 0 (experiments/E142_step0.py)

def part_step0(pieces: Optional[Sequence[str]] = None) -> dict:
    """Step 0, no simulation: E142_step0.run (every piece in order, or `pieces`)."""
    from experiments import E142_step0 as S0
    return S0.run(pieces)


# ---------------------------------------------------------------------------
# step 1: hal's module (E142_imposed_calyx_runs.py), its runs and its score

def run_step1(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False):
    """Step 1's runs on one lane: hal's `part_step1` (E135's step-2 runner with A1's calyx factors; guarded by
    step1_authorized.json)."""
    from experiments import E142_imposed_calyx_runs as S1
    return S1.part_step1(lane, reverse, reclaim)


def score_step1() -> dict:
    """Step 1's score and classes ('apl_hides_step', 'dendrites_carry', 'no_class'): E136's reading of the boutons,
    claws and somata on T against the held-out half, E135's saved runs at the family's calyx re-scored the same way,
    A1 beside; realism: hal's `part_score` (E142_imposed_calyx_runs.py; step1_classes.json and the score_step1_*.csv
    tables in the results folder)."""
    from experiments import E142_imposed_calyx_runs as S1
    return S1.part_score()


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'step0', 'step1', 'score'))
    ap.add_argument('--pieces', nargs='*', default=None,
                    help="step 0's pieces to run (default all, in order; E142_step0.PIECES)")
    ap.add_argument('--lane', default=None, help="step 1's lane (hal's)")
    ap.add_argument('--reverse', action='store_true', help="step 1: run the planned points in reverse")
    ap.add_argument('--reclaim', action='store_true', help="step 1: reclaim stale claims")
    a = ap.parse_args(argv)
    if a.part == 'lock':
        part_lock()
    elif a.part == 'step0':
        part_step0(a.pieces)
    elif a.part == 'step1':
        run_step1(a.lane, a.reverse, a.reclaim)
    else:
        score_step1()


if __name__ == '__main__':
    main()
