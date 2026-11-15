#!/usr/bin/env python3
"""E99: the evenness of APL's weights among KCs and over a KC's claws: what loses the mean APL effect when APL's counted
synapses take the bouton path's place at the KCs (the weaker growth of a KC's amount with its number of claws, the
inequality among KCs with the same number of claws, or the uneven placement over a KC's claws, each changed alone at fixed
scales), and whether each KC's count of all of APL's synapses, spread evenly over its claws, holds the mean.

The entry is item E99 of experiment_priorities.md (proposed by viki 2026-11-15 12:48; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-15 12:27). Built on E96 (`experiments/E96_offclaw_apl.py`: the counts, the
components of its modes, `e96_weights`, `reading96`, the per-KC tables of its part 0, the readings, the choice, the gate,
the equal rule and part D's plan, starts and scoring), E95 (`StaticModel`, `kc_losses`, `d_unit`, `d_realism`, its part D
runs), E94 (the stage readers, D0) and E92 (the bracket).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B.
S: real, P1, P2, M1. T: P3, P4, M2, M3 (only part C's gate and part E). Halves A / B E77's, C / D E86's. The input is the
lab's in every run. Everything as E96's part B: the bouton path at dose 0.1; a claw's APL->KC weight is the family's
plus d x unit x the component (the component per claw has a mean of 1 over the claws, times the connectome's mean
normalized weight, E10's convention, so its size d is in E95's unit); the tied APL scale bisected to the wiring's APL-on
count in E92's bracket; the KC-KC scale and the somatic threshold at the fit's; the APL-off run of every point the fit's.

The settings (the entry's Definitions; `shape99`, from the counts of E96's `wiring_counts` with the loader's filters):
for KC k with n_k claws, a_k its APL->KC synapses on no claw in the calyx (a^c_k; with the peduncle's a^c_k + a^p_k),
c_kj its synapses on kept claw j and c_k their sum (the model's weights). A claw-count group G is the KCs with the same
number of claws, a group with fewer than 20 KCs pooled with the next lower count (`claw_groups`; on the real wiring the
10 KCs with 12 claws join the 35 with 11); mu_G = sum of a_k over G / sum of n_k over G. Each component is scaled to a
mean of 1 over the claws.
- K1 (in proportion): (a_k / n_k) / mu_G on each claw of KC k (every group the same amount per claw).
- K2 (the group's mean): mu_G on each claw of a KC of G. K1 x K2 is E96's off-claw component, claw by claw.
- K3 (half the inequality): mu_G + (a_k / n_k - mu_G) / 2 (all of it: E96's off-claw component; none: K2).
- K4 (presence): the same on every claw of a KC with a_k > 0, 0 on the others.
- W1: c_k / n_k on each claw (the on-claw counts per KC, evenly). W2: c_kj / (c_k / n_k) (equal amounts per KC, placed
  as the on-claw counts lie; a KC with no synapse on a kept claw 1 on each claw). W1 x W2 is the connectome's per-claw
  component (E95's D2), claw by claw.
- T1: (c_k + a^c_k + a^p_k) / n_k on each claw (E10's --aplkc-flat allsyn_even); T2: (c_k + a^c_k) / n_k
  (allsyn_even_calyx).
K1 to K4 on the calyx's counts ('offclaw') or with the peduncle's ('offclawped' in names).

Variants: K<1-4>_<offclaw|offclawped>_d<d>, W<1|2>_d<d>, T<1|2>_d<d>; the saved corners keep their names (E95's
D1_uniform_d<d>, D2_connectome_d1; E96's B1_offclaw_d<d>, B2_offclawped_d1, B6_allsyn_d<d>, the two B4 permutations;
E94's D0; 'ref' the unchanged fits).

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B, C and D refuse to run until it passes): through this
  experiment's path (`point_weights99`, E96's `reading96`) on the real wiring, with 0 count differences and the APL
  factor within 1e-12: the uniform component at d = 1 reproduces E95's D1_uniform_d1, the connectome's E95's
  D2_connectome_d1, the off-claw E96's B1_offclaw_d1, all counted synapses E96's B6_allsyn_d1 (4 runs, about 5 min);
  without simulations on each of S's wirings (`shape_identities`, from E96's cached per-KC and per-claw tables and the
  live counts): every component has a mean of 1 over the claws; K1 x K2 is the off-claw component and W1 x W2 the
  connectome's, claw by claw; K1's mean per claw is 1 in every group; K2 has one value within a group and K2 and K4 one
  within a KC (K4 one value on the claws it does not leave at 0); K3 with all of the inequality is the off-claw component
  and with none K2; W1's, T1's and T2's sums per KC are the counts' (in the components' unit); W2's sum per KC is n_k;
  the cached tables are the live counts; T1 and T2 are E10's new modes exactly.
- A (no simulations; one thread): A1 the shapes (a1_shapes.csv: for E95's and E96's five saved components, E96's two
  permutations (their a_k as the runs saved them) and the ten new settings, on each of S's wirings, the variance of the
  weight per claw split between claw-count groups, among the KCs of a group and within KCs, the share of claws and of KCs
  at zero, the log-log slope of a KC's total on its claw count, the CV among KCs with the same number of claws, r with
  E95's b_k); A2 the predictions (a2_predictions.csv / a2_predictions.json: the two descriptions, (i) the mean APL effect
  linear in the variance per claw as a whole and (ii) in its three parts, fitted to the seven saved means at size 1 of
  S's population under the calibrated readout (E96's readings.csv), their predictions for every new setting at size 1
  and the classes of part C they predict; written before part B's first point, with its time, and not rewritten once
  part B has a row); A3 who responds (a3_fifths.csv, a3_claws.csv: in the saved runs at size 1 on S's wirings, E95's D1
  and D2 and E96's B1, B2 and B6, the odors per KC with APL on and off and the share of the responses by fifth of the KC's
  calyx off-claw count per claw relative to its claw-count group (the entry's probe) and of the run's own component per
  claw relative to its group, by claw class and by KC type; under D2 and B6 the claw stage on each KC's most and least
  weighted claw); A4 the bound's yardstick (a4_yardstick.csv: the CV among KCs with the same number of claws of the
  off-claw, on-claw and peduncle counts in E96's a2_connectomes_kc.csv beside the model's); a_summary.json.
- B (`--block I|II`, default both; `--wirings`, default S; resumable, runs_B.csv; runs under
  outputs/E99_apl_weight_evenness/runs/ with their stage matrices, each component per KC and per claw and the APL->KC
  weights per claw): block I (24 points) K1 and K2 on the calyx's counts, W1, W2, T1 and T2, all at d = 1; block II (40)
  K1 and K2 on the calyx's counts at d = 2 and on the counts with the peduncle's at d = 1, K3 and K4 on the calyx's counts
  at d = 1, T1 and T2 at d = 0.5 and 2. Refuses until part 0 has passed and A2's predictions are written.
- score (`--no-co`, `--no-zcell`): readings.csv per variant and unit (each wiring; S and T where all four have a run) and
  readout ('cal', 'calK', 'raw', 'zcell', 'zcell_lin', as E96) against every target of `targets89` (E92's `score_row92`
  with E94's group effects), this experiment's points and the saved corners (E95's D1 at the three sizes and D2 at d = 1,
  E96's B1 at d = 1 and 2, B2 at d = 1, B4 and B6 at the three sizes, E94's D0, the unchanged fits) read the same way;
  pairs.csv; line.csv (E92's line); coresponse.csv (the alcohol pairs' co-response change); b_fifths.csv and b_claws.csv
  (A3's reading for every point); part C's verdicts (c_shares.csv, c_verdicts.json: the shares and classes between KCs at
  d = 1 on the calyx's counts with the rule by wiring and the two replicates, within a KC, T1 and T2 against B6 at each
  size, the real wiring against P1, P2 and M1 for W2, D2 and B6, how even is even enough, the KCs with none); the choice,
  gate and equal rule for W1, T1 and T2 (choice.csv, gate.csv, equal.csv, vs_uniform.csv); part D's rule (d_rule.json);
  verdict.json, the figure and comparison_table.md.
- C: the choice per half among W1, T1, T2 and the unchanged fits; a chosen setting with something to confirm is run on
  T (runs_C.csv) with E95's uniform component at the same size as rule (v)'s control (run here on T through this path);
  then `score` without the co-response. The gate is E83's as E91 applied it (E92's `gate_direction92` / `gate92`).
- Dcheck (once, before part D): E10_joint_fit.py --aplkc-flat allsyn_even / allsyn_even_calyx at size 0 against E94's D1
  fit's first evaluation on the real wiring (d_check.json; E96's `d_check` on these modes).
- D (`--pop S|T`, `--start own|e96`, `--wirings` for a lane, `--dry-run`, `--force`): the family fitted again at the
  connectome's bouton share with all counted synapses per KC spread evenly (`E10_joint_fit.py -k 2 --aplpn-hold F
  --aplkc-flat allsyn_even|allsyn_even_calyx --aplkc-flat-unit U --aplkc-flat-max 8 --rate-target ctrl_frac --loss apl
  --apl-readout both`, 200 evaluations), only by its rule (`d_rule99`): a setting of T1 or T2 on S against all recordings
  with a mean APL effect of +0.168 or more, a pattern part within 0.0015 of E95's uniform component's at the same size or
  lower, its four runs inside the bands and not at a bound; T2 if it meets the rule and its best pattern part (among its
  settings that meet it) is within 0.0005 of T1's, otherwise the one that meets it with the lower pattern part. Two
  starts per S wiring: 'own' (the wiring's best point of that component in part B; tag suffix _g1) and 'e96' (E96's kept
  refit of that wiring, its parameters and fitted size, the component exchanged; _g1e). T's fits (`--pop T`, start 'own'
  only) only on an equal of S's fitted family: each wiring's own fit with the bouton path at a tenth and the component at
  the median of S's fitted sizes (the entry's part E's refits).
- Dscore: per S wiring the fit with the lower final value of its own loss (d_fits rows in d_scores.csv); S's population
  by E94's D1 rule (an equal as a fitted family), as good as the uniform component (within 0.0015 of E95's refit, or
  lower), beside E96's refit with the calyx's off-claw synapses and E94's; the eight by E91's part D rule.
- E: documented stub (`part_e`): runs only on an equal as a fitted family.
- steps (`--variants`): the step from the claws to the somata for chosen or equal settings (E96's `part_steps`).
- variants: variants.csv and model_pair_distances_long.csv for the archive.
- table: comparison_table.md and the figure from what exists.

Deviations from the entry, forced by the code or the data:
- The counts are E96's, with the loader's filters (E96's first deviation): a^c_k is the counts' no_claw (7,234 synapses on
  the model's KCs), c_kj the loader's kept counts (including the 238 synapses with claw ID -1 on the wiring's claws -1).
- The pooling of small claw-count groups runs from the highest count down (a group still short after pooling joins the
  next lower count again); a lowest group still short would join the next higher one (it does not occur: every wiring has
  the real wiring's claws). A group with no off-claw synapse would get K1 = 1 on its claws (it does not occur).
- "K2 and K4 have one value within a group and within a KC": K2 has one value within a group; K4 has one value within a
  KC and one value on every claw it does not leave at 0 (its value within a group is that or 0).
- A3's "by fifth of the KC's component per claw relative to its claw-count group" is given two ways: by the calyx's
  off-claw count per claw relative to the group (K1's value: the same fifths for every run, the entry's probe table) and
  by the run's own component per claw relative to the group (none for the uniform component). "The same for a KC's most
  and least weighted claws under D2 and B6" is read on the claw stage: per KC with at least two claws and an uneven
  component, its most and least weighted claw's window state with APL on (the run) and off (the fit), the share of
  (claw, odor) above 0 and the mean, by claw class (a3_claws.csv).
- A1's permutations are read from the a_k that E96's B4 runs saved per KC (outputs/E96_offclaw_apl/runs/ or the archive),
  so that they are the components those runs had.
- How even is even enough: the loss against the CV among KCs with the same number of claws (the median over claw counts
  2 to 10 of the KC totals, E96's `_within`) is a power law through O and K3 (K2 at 0), and the CV at which 0.01 of the
  mean is lost is read from it (extrapolated when the gap is below 0.01; said in c_verdicts.json).
- T's refits of the entry's part E run as `--part D --pop T` (E96's practice); part E's other readings are a stub.
- Part D's 'e96' start reads E96's kept refit per wiring from E96's d_scores.csv ('picked'), its parameters from E96's
  fits.csv, and refuses if E96's unit differs from this experiment's.

Outputs: results/E99_apl_weight_evenness/; spike counts, responses, stage matrices, components and weights under
outputs/E99_apl_weight_evenness/ (runs/, D/ logs); part D's fits where E10_joint_fit.py writes them for this experiment.

Usage (from the repository root; one 16-thread process at a time on viki; check `uptime` first; long runs with -u; the
simulation parts refuse to start while a checkout has uncommitted changes, since their rows record the commits):
    nice -n 19 env OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part A                              # A1 to A4; no simulations; A2's predictions before part B's first point
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part 0 > outputs/E99_apl_weight_evenness/0.log 2>&1
  part B (64 points, about 1 h; resumable; block I first, then II, or both):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part B --block I > outputs/E99_apl_weight_evenness/B_I.log 2>&1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part B --block II > outputs/E99_apl_weight_evenness/B_II.log 2>&1
  then the readings (co-response readings in 4 processes):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part score > outputs/E99_apl_weight_evenness/score.log 2>&1
  (C) only on a choice (it reads the readings; nothing runs when nothing is to confirm):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part C
  part D, only by its rule (`--dry-run` prints the rule, the starts and the commands without running):
    ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part D --dry-run
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part Dcheck > outputs/E99_apl_weight_evenness/Dcheck.log 2>&1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part D --pop S --start own > outputs/E99_apl_weight_evenness/D_own.log 2>&1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part D --pop S --start e96 > outputs/E99_apl_weight_evenness/D_e96.log 2>&1
  (`--wirings real,P1` runs a lane.) If the second start runs on hal as a lane, copy its fits' outputs afterwards:
  `mkdir -p <shared>/data/model_outputs/E99_apl_weight_evenness/from_hal && cp outputs/E99_apl_weight_evenness/*_g1e_*.parquet
  <shared>/data/model_outputs/E99_apl_weight_evenness/from_hal/` (Dscore and variants read them there).
    ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part Dscore
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E99_apl_weight_evenness \
        --part D --pop T                      # only on an equal as a fitted family
    ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part steps
    ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part variants
    ./venv/bin/python -u -m experiments.E99_apl_weight_evenness --part table
    ./venv/bin/python archive_tools/archive_model_outputs.py E99_apl_weight_evenness
    ./venv/bin/python archive_tools/migrate_model_output_params.py --apply
"""
from __future__ import annotations

import argparse
import json
import re
import socket
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E89_apl_breadth_weighting as E89
from experiments import E91_summed_drive_inhibition as E91
from experiments import E92_gain as E92G
from experiments import E94_bouton_share as E94
from experiments import E95_static_apl as E95
from experiments import E96_offclaw_apl as E96

EXPERIMENT = 'E99_apl_weight_evenness'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
KC = H.KC_ID_LEVEL

# ---------------------------------------------------------------------------
# the settings (the entry's Design and Definitions; fixed before any run)

APL_BRACKET = E92G.APL_BRACKET                      # (0.01, 10): the tied APL scale, E92's bracket
RAW_CRIT = E92G.RAW_CRIT                            # a model response is two or more spikes
CLAW_CLASSES = E94.CLAW_GROUPS                      # 1-4, 5-8, 9+ claws
D_DOSE = E96.D_DOSE                                 # the bouton path at a tenth (E94's anatomical dose)
GROUP_MIN = 20                                      # a claw-count group with fewer KCs is pooled with the next lower count
K3_FRAC = 0.5                                       # K3: half the inequality
SHAPES = ('K1', 'K2', 'K3', 'K4', 'W1', 'W2', 'T1', 'T2')
K_SHAPES, W_SHAPES, T_SHAPES = SHAPES[:4], ('W1', 'W2'), ('T1', 'T2')
COUNTS = ('offclaw', 'offclaw_ped')                 # K1 to K4: the calyx's off-claw counts, or with the peduncle's
COUNTS_CODE = {'offclaw': 'offclaw', 'offclaw_ped': 'offclawped'}
CODE_COUNTS = {v: k for k, v in COUNTS_CODE.items()}
T_MODE = {'T1': 'allsyn_even', 'T2': 'allsyn_even_calyx'}   # E10_joint_fit.py's --aplkc-flat modes
MODE_T = {v: k for k, v in T_MODE.items()}
CONSTRUCTIONS = ('K1', 'K2', 'K3', 'K4', 'W2')      # change one property of the counts: never chosen, equal or fitted
CHOICE_ARMS = ('W1', 'T1', 'T2')                    # counts as counted, spread evenly: the choice, the equal rule
SIZES = (0.5, 1.0, 2.0)


def vname(shape: str, d: float, counts: Optional[str] = None) -> str:
    """A setting's name: K<1-4>_<offclaw|offclawped>_d<d> (K needs `counts`), W<1|2>_d<d>, T<1|2>_d<d>."""
    if shape in K_SHAPES:
        if counts not in COUNTS:
            raise ValueError(f'{shape} needs the counts ({COUNTS})')
        return f'{shape}_{COUNTS_CODE[counts]}_d{d:g}'
    if shape not in SHAPES or counts is not None:
        raise ValueError(f'{shape!r} with counts {counts!r}')
    return f'{shape}_d{d:g}'


BLOCKS: Dict[str, Tuple[str, ...]] = {
    'I': (vname('K1', 1, 'offclaw'), vname('K2', 1, 'offclaw'), vname('W1', 1), vname('W2', 1), vname('T1', 1),
          vname('T2', 1)),
    'II': (vname('K1', 2, 'offclaw'), vname('K2', 2, 'offclaw'), vname('K1', 1, 'offclaw_ped'), vname('K2', 1, 'offclaw_ped'),
           vname('K3', 1, 'offclaw'), vname('K4', 1, 'offclaw'), vname('T1', 0.5), vname('T1', 2), vname('T2', 0.5),
           vname('T2', 2)),
}
REF = 'ref'                                         # the unchanged fits (source 'saved')
D0 = E95.D0                                         # dose 0.1 without a component (E94's point on S)
UNIFORM = {d: f'D1_uniform_d{d:g}' for d in SIZES}  # E95's uniform component (U at each size)
CONNECTOME = 'D2_connectome_d1'                     # E95's per-claw weights at d = 1 (D)
OFFCLAW = {1.0: 'B1_offclaw_d1', 2.0: 'B1_offclaw_d2'}   # E96's calyx off-claw component (O)
OFFCLAW_PED = 'B2_offclawped_d1'                    # E96's off-claw with the peduncle's at d = 1 (the replicate's O)
ALLSYN = {d: f'B6_allsyn_d{d:g}' for d in SIZES}    # E96's all counted synapses, the on-claw part where it lies
PERMS = ('B4_offclawped_d1_permclaw1', 'B4_offclawped_d1_permall1')
CORNERS: Dict[str, Tuple[str, ...]] = {
    'E95': tuple(UNIFORM.values()) + (CONNECTOME,),
    'E96': (OFFCLAW[1.0], OFFCLAW[2.0], OFFCLAW_PED) + tuple(ALLSYN.values()) + PERMS,
}

# part C (the entry's rules)
CLASS_HI, CLASS_LO = 0.6, 0.25                      # a share of 0.6 or more carries it; below 0.25 does not
BY_WIRING_MIN = 3                                   # the class holds by wiring if 3 of the 4 wirings give it
BOUND_LOSS = 0.01                                   # the CV among KCs of a claw count at which 0.01 of the mean is lost
BOUND_MIN_SHARE = 0.25                              # read only if K2 recovers a quarter of the gap or more
NONE_HI, NONE_LO = 0.5, 0.25                        # U - K4: half of the gap or more / a quarter or less
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
HALVES = E91.HALVES
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
FLIES_MEAN = E89.FLIES_MEAN                         # +0.178

# part D
D_MEAN_MIN = FLIES_MEAN - 0.01                      # +0.168: within 0.01 of the flies'
D_VS_UNIFORM = 0.0015                               # the pattern part within this of the uniform component's, or lower
D_T2_TIE = 0.0005                                   # T2 if within this of T1's
D_FLAT_BOUNDS = (0.0, 8.0)
D_MAXFEV = 200
D_TAGS = {'own': '_g1', 'e96': '_g1e'}
D_START_FILE = 'd_start_params.json'
D_AS_GOOD = 0.0015
E94_REFIT_S = E95.E94_REFIT_S                       # 0.0593
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = E95.UNCHANGED_8                       # 0.0482

_CACHE: Dict[tuple, object] = {}


def _write_json(path: Path, obj) -> None:
    E89._write_json(path, obj)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def out_dir() -> Path:
    return H.results_dir(EXPERIMENT)


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _commits(strict: bool = True) -> dict:
    return E89._commits_clean() if strict else dict(H.code_commits())


def _wirings(arg: str) -> Tuple[str, ...]:
    return E91._wirings(arg)


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


# ---------------------------------------------------------------------------
# names

_RE_K = re.compile(r'^(K[1-4])_(offclaw|offclawped)_d([0-9.]+)$')
_RE_WT = re.compile(r'^([WT][12])_d([0-9.]+)$')


def _size(v: str, d: str) -> float:
    if f'{float(d):g}' != d or not float(d) > 0:
        raise ValueError(f'variant {v!r}: the size is not written canonically or not positive')
    return float(d)


def parse_variant(v: str) -> dict:
    """A variant's fields: shape (SHAPES, or the saved component's mode), counts (K1 to K4: 'offclaw' / 'offclaw_ped'),
    part ('B' this experiment's settings; 'e96' / 'e95' the saved corners; 'ref'; 'control' D0), arm and screen (the
    shape for this experiment's settings; 'e96', 'e95', 'ref', 'control' otherwise), mode (the shape, or the saved
    component's mode), size, pool (0), dose, permute and seed (E96's B4), construction (K1 to K4 and W2: a construction
    that is never chosen, called an equal or fitted), equal (W1, T1, T2: the choice and the equal rule apply)."""
    out = dict(variant=v, arm=None, part=None, shape=None, counts=None, mode=None, size=0.0, pool=0.0, dose=D_DOSE,
               permute=None, seed=None, control=False, control_of='', screen=None, construction=False, equal=False)
    m = _RE_K.match(v)
    if m:
        shape, code, d = m.groups()
        out.update(arm=shape, part='B', shape=shape, counts=CODE_COUNTS[code], mode=shape, size=_size(v, d), screen=shape,
                   construction=True)
        return out
    m = _RE_WT.match(v)
    if m:
        shape, d = m.groups()
        out.update(arm=shape, part='B', shape=shape, mode=shape, size=_size(v, d), screen=shape,
                   construction=shape in CONSTRUCTIONS, equal=shape in CHOICE_ARMS)
        return out
    p = E96.parse_variant(v)                         # E96's names, E95's part D settings, ref, D0 (raises otherwise)
    if p['part'] == 'ref':
        out.update(arm=REF, part=REF, dose=1.0, screen=REF)
    elif v == D0:
        out.update(arm=D0, part='control', control=True, control_of='C', screen='control')
    elif p['part'] == 'e95':
        out.update(arm=p['arm'], part='e95', shape=p['mode'], mode=p['mode'], size=float(p['size']), dose=float(p['dose']),
                   screen='e95')
    elif p['part'] == 'B':
        out.update(arm=f"e96_{p['arm']}", part='e96', shape=p['mode'], mode=p['mode'], size=float(p['size']),
                   pool=float(p['pool']), dose=float(p['dose']), permute=p['permute'], seed=p['seed'], screen='e96')
    else:
        raise ValueError(f"variant {v!r}: E96's {p['part']} settings are not read here")
    return out


def block_variants(block: str) -> List[str]:
    return list(BLOCKS[block])


def control_of99(variant: str) -> str:
    """The rule (v) control of a chosen setting: E95's uniform component at the same size (D1_uniform_d<d>)."""
    pv = parse_variant(variant)
    if pv['part'] != 'B' or pv['screen'] not in CHOICE_ARMS:
        raise ValueError(f'{variant}: not a setting of the choice')
    return f'D1_uniform_d{pv["size"]:g}'


# ---------------------------------------------------------------------------
# the counts, the groups and the shapes

def count_tables(label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per KC, per claw) of a wiring's live counts (E96's `wiring_counts`, the loader's filters): per KC kc_id, n_claws,
    kept, no_claw, peduncle (and the counts' other columns), per claw kc_id, claw_id, kept (the wiring's claws' order)."""
    from experiments import E10_joint_fit as J
    wn = E96.wiring_wnorm(label)
    per_kc, per_claw = J.aplkc_flat_counts(wn.index)
    return tables_from_counts(per_kc, per_claw)


def tables_from_counts(per_kc: pd.DataFrame, per_claw: pd.Series) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The (per KC, per claw) frames of `shape99` from pyolfsysm's counts (per KC indexed by KC id; per claw a Series
    indexed by the claws with kc_id and claw_id levels)."""
    k = per_kc.copy()
    k.index = k.index.rename(KC)
    k = k.reset_index()
    ix = per_claw.index
    claws = pd.DataFrame({KC: np.asarray(ix.get_level_values(KC)), 'claw_id': np.asarray(ix.get_level_values('claw_id')),
                          'kept': per_claw.to_numpy(dtype=float)})
    return k, claws


def cached_tables(label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """E96's part 0 cache of a wiring (outputs/E96_offclaw_apl/W/): per KC (the outputs' order; with b and c) and per
    claw (with wnorm and E96's components comp_<mode>)."""
    key = ('tables', label)
    if key not in _CACHE:
        _CACHE[key] = E96.load_kc_tables96(label)
    kc, claws = _CACHE[key]
    return kc.copy(), claws.copy()


def claw_groups(n_claws: np.ndarray, min_kcs: int = GROUP_MIN) -> np.ndarray:
    """Each KC's claw-count group G, labelled by the group's lowest claw count: the KCs with the same number of claws, a
    group with fewer than `min_kcs` KCs pooled with the next lower count (from the highest count down, so that a run of
    small groups pools until it has `min_kcs`); a lowest group still short joins the next higher one."""
    n = np.asarray(n_claws).astype(int)
    vc = pd.Series(n).value_counts().sort_index()
    groups: List[List[int]] = []
    acc: List[int] = []
    for k in reversed(vc.index.tolist()):
        acc.append(int(k))
        if int(vc[acc].sum()) >= min_kcs:
            groups.append(acc)
            acc = []
    if acc:
        if groups:
            groups[-1].extend(acc)
        else:
            groups.append(acc)
    label = {k: min(g) for g in groups for k in g}
    return np.array([label[int(k)] for k in n])


def _check_tables(kc: pd.DataFrame, claws: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray]:
    k = kc.set_index(KC)
    if k.index.has_duplicates:
        raise ValueError('the per-KC table has a KC twice')
    kc_of = claws[KC].to_numpy()
    n = pd.Series(kc_of).value_counts()
    if set(n.index) != set(k.index):
        raise ValueError("the per-KC table's KCs are not the claws' KCs")
    if not np.array_equal(n.reindex(k.index).to_numpy(dtype=float), k['n_claws'].to_numpy(dtype=float)):
        raise ValueError("the per-KC table's claw counts are not the claws'")
    return k, kc_of


def group_of_kc(kc: pd.DataFrame) -> pd.Series:
    """Each KC's claw-count group (`claw_groups`), indexed by KC id."""
    return pd.Series(claw_groups(kc['n_claws'].to_numpy()), index=kc[KC].to_numpy(), name='group')


def offclaw_amounts(kc: pd.DataFrame, counts: str = 'offclaw') -> pd.Series:
    """a_k per KC id: the calyx's off-claw counts a^c_k ('offclaw') or with the peduncle's ('offclaw_ped')."""
    from experiments import E10_joint_fit as J
    if counts not in COUNTS:
        raise ValueError(f'counts {counts!r}')
    return J.aplkc_flat_offclaw(counts, kc.set_index(KC))


def group_means(a: pd.Series, n: pd.Series, grp: pd.Series) -> pd.Series:
    """mu_G per KC: the sum of a over the KC's group over the sum of n (all indexed by KC id)."""
    sa = a.groupby(grp).sum()
    sn = n.astype(float).groupby(grp).sum()
    return (sa / sn).reindex(grp.to_numpy()).set_axis(a.index)


def shape99(name: str, kc: pd.DataFrame, claws: pd.DataFrame, counts: str = 'offclaw',
            frac: Optional[float] = None) -> np.ndarray:
    """The component per claw (the claws' order, mean 1 over the claws) of a setting `name` (SHAPES) from a wiring's
    per-KC table (kc_id, n_claws, no_claw, peduncle, kept) and per-claw table (kc_id, kept): K1 to K4 on `counts`
    ('offclaw' a^c_k, 'offclaw_ped' a^c_k + a^p_k), K3 with `frac` of the inequality kept (default K3_FRAC; 1 is E96's
    off-claw component, 0 is K2); W1, W2 from the on-claw counts; T1, T2 E10's allsyn_even / allsyn_even_calyx
    (`aplkc_flat_shape` of the KC's whole count, `aplkc_flat_per_kc`)."""
    from experiments import E10_joint_fit as J
    k, kc_of = _check_tables(kc, claws)
    n = k['n_claws'].astype(float)
    if name in K_SHAPES:
        a = offclaw_amounts(kc, counts)
        grp = group_of_kc(kc)
        mu = group_means(a, n, grp)
        r = a / n
        if name == 'K1':
            v = (r / mu).where(mu > 0, 1.0)
        elif name == 'K2':
            v = mu
        elif name == 'K3':
            f = K3_FRAC if frac is None else float(frac)
            v = mu + f * (r - mu)
        else:
            v = (a > 0).astype(float)
        x = v.reindex(kc_of).to_numpy(dtype=float)
    elif name in W_SHAPES:
        kept = claws['kept'].to_numpy(dtype=float)
        c = pd.Series(kept).groupby(kc_of).sum().reindex(k.index)
        per = (c / n).reindex(kc_of).to_numpy(dtype=float)
        if name == 'W1':
            x = per
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                x = np.where(per > 0, kept / np.where(per > 0, per, 1.0), 1.0)
    elif name in T_SHAPES:
        mode = T_MODE[name]
        return J.aplkc_flat_shape(mode, kc_of, J.aplkc_flat_per_kc(mode, k))
    else:
        raise ValueError(f'unknown setting {name!r}')
    if (x < 0).any() or not np.isfinite(x).all():
        raise ValueError(f'{name}: the component must be finite and not negative')
    m = float(x.mean())
    if not m > 0:
        raise ValueError(f'{name}: the component is zero on every claw')
    return x / m


def setting_shape(v: str, kc: pd.DataFrame, claws: pd.DataFrame) -> np.ndarray:
    """The per-claw shape (mean 1) of a variant of this experiment (`parse_variant`'s shape and counts)."""
    pv = parse_variant(v)
    if pv['part'] != 'B':
        raise ValueError(f'{v}: not a setting of this experiment')
    return shape99(pv['shape'], kc, claws, pv['counts'] or 'offclaw')


def kc_totals(x: np.ndarray, kc_of_claw: np.ndarray, ids: Sequence) -> np.ndarray:
    """A per-claw component summed over each KC's claws, in the order `ids`."""
    return pd.Series(np.asarray(x, dtype=float)).groupby(np.asarray(kc_of_claw)).sum().reindex(list(ids)).to_numpy(dtype=float)


def variance_split(x: np.ndarray, kc_of_claw: np.ndarray, group_of_claw: np.ndarray) -> dict:
    """The variance over the claws of a per-claw component, relative to its squared mean (so that, at a mean of 1, the
    parts add up to the squared CV per claw): 'between' the claw-count groups (the claws' groups' means about the mean),
    'among' the KCs of a group (the KCs' means about their group's), 'within' KCs (the claws about their KC's mean),
    'total' their sum, all claw-weighted (population variances)."""
    x = np.asarray(x, dtype=float).ravel()
    m = float(x.mean())
    s = pd.DataFrame(dict(x=x, k=np.asarray(kc_of_claw), g=np.asarray(group_of_claw)))
    km = s.groupby('k').x.transform('mean').to_numpy()
    gm = s.groupby('g').x.transform('mean').to_numpy()
    m2 = m * m
    out = dict(between=float(((gm - m) ** 2).mean() / m2), among=float(((km - gm) ** 2).mean() / m2),
               within=float(((x - km) ** 2).mean() / m2))
    out['total'] = float(((x - m) ** 2).mean() / m2)
    out['parts_sum'] = out['between'] + out['among'] + out['within']
    return out


# ---------------------------------------------------------------------------
# the saved components (A1's comparisons) from E96's cache and E96's saved runs

SAVED_SHAPES = {'uniform': 'D1_uniform_d1', 'connectome': 'D2_connectome_d1', 'offclaw': 'B1_offclaw_d1',
                'offclaw_ped': 'B2_offclawped_d1', 'allsyn': 'B6_allsyn_d1', 'perm_claw': PERMS[0], 'perm_all': PERMS[1]}
NEW_SHAPES = {'K1': ('K1', 'offclaw'), 'K2': ('K2', 'offclaw'), 'K3': ('K3', 'offclaw'), 'K4': ('K4', 'offclaw'),
              'K1p': ('K1', 'offclaw_ped'), 'K2p': ('K2', 'offclaw_ped'), 'W1': ('W1', None), 'W2': ('W2', None),
              'T1': ('T1', None), 'T2': ('T2', None)}


def new_variant(key: str, d: float = 1.0) -> str:
    shape, counts = NEW_SHAPES[key]
    return vname(shape, d, counts)


def _e96_saved_ak(tag: str) -> Optional[pd.DataFrame]:
    for d in E96._dirs('runs'):
        p = d / f'{tag}_component.parquet'
        if p.exists():
            return pd.read_parquet(p)
    return None


def saved_shape(name: str, kc: pd.DataFrame, claws: pd.DataFrame, wiring: str) -> np.ndarray:
    """A saved component's per-claw shape (mean 1) on a wiring: 'uniform', 'connectome' (wnorm over its mean), E96's
    modes (the cache's comp_<mode>), E96's permutations ('perm_claw', 'perm_all': a_k as E96's B4 run on that wiring saved
    it, spread evenly, or recomputed from the counts in pyolfsysm's order when that run is not on this host)."""
    from experiments import E10_joint_fit as J
    if name == 'uniform':
        return np.ones(len(claws))
    if name == 'connectome':
        w = claws['wnorm'].to_numpy(dtype=float)
        return w / w.mean()
    if name in E96.MODES:
        return claws[f'comp_{name}'].to_numpy(dtype=float)
    if name in ('perm_claw', 'perm_all'):
        tag = run_tag(SAVED_SHAPES[name], wiring)
        f = _e96_saved_ak(tag)
        kc_of = claws[KC].to_numpy()
        if f is not None:
            a = pd.Series(f['a_k'].to_numpy(dtype=float), index=f[KC].to_numpy())
        else:
            per_kc, _ = J.aplkc_flat_counts(pd.MultiIndex.from_frame(claws[[KC, 'claw_id']]))
            a = E96.offclaw_counts('offclaw_ped', per_kc, 'claw' if name == 'perm_claw' else 'all', E96.PERM_SEED)
        return J.aplkc_flat_shape('offclaw_ped', kc_of, a)
    raise ValueError(f'unknown saved component {name!r}')


# ---------------------------------------------------------------------------
# A1: the shapes

def shape_row(x: np.ndarray, kc: pd.DataFrame, claws: pd.DataFrame) -> dict:
    """A1's row of a per-claw component: the variance split (`variance_split` on the claw-count groups), the CV per
    claw, the share of claws and of KCs at zero, a KC total's log-log slope on its claw count, its CV among KCs with the
    same number of claws (E96's `_within`, claw counts 2 to 10) and its r with E95's b_k (where `kc` has b)."""
    kc_of = claws[KC].to_numpy()
    grp = group_of_kc(kc)
    split = variance_split(x, kc_of, grp.reindex(kc_of).to_numpy())
    tot = kc_totals(x, kc_of, kc[KC])
    n = kc['n_claws'].to_numpy(dtype=float)
    b = kc['b'].to_numpy(dtype=float) if 'b' in kc else None
    row = dict(mean=float(np.mean(x)), **split, per_claw_cv=float(np.std(x) / np.mean(x)),
               claws_zero_share=float((np.asarray(x) == 0).mean()), kcs_zero_share=float((tot == 0).mean()),
               loglog_slope_n=E96._loglog(n, tot), cv_kcs=float(tot.std() / tot.mean()),
               r_n_claws=E95._corr(n, tot)[0])
    w = E96._within(n, tot, b)
    row.update(w)
    if b is not None:
        row['r_b'] = E95._corr(b, tot)[0] if np.ptp(tot) > 0 else np.nan
    return row


def a1_shapes(out: Path, wirings: Sequence[str] = S_WIRINGS) -> pd.DataFrame:
    """A1: a1_shapes.csv, per wiring of `wirings` and component (the saved five, E96's two permutations, the ten new
    settings) `shape_row`, with the variant that runs it at size 1."""
    rows = []
    for w in wirings:
        kc, claws = cached_tables(w)
        for name, v in SAVED_SHAPES.items():
            rows.append(dict(wiring=w, setting=name, variant=v, kind='saved',
                             **shape_row(saved_shape(name, kc, claws, w), kc, claws)))
        for key, (shape, counts) in NEW_SHAPES.items():
            rows.append(dict(wiring=w, setting=key, variant=new_variant(key), kind='construction' if shape in CONSTRUCTIONS
                             else 'counted', **shape_row(shape99(shape, kc, claws, counts or 'offclaw'), kc, claws)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a1_shapes.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# A2: the two descriptions and their predictions (written before part B's first point)

SAVED7 = ('uniform', 'offclaw_ped', 'perm_claw', 'perm_all', 'allsyn', 'offclaw', 'connectome')
PARTS3 = ('between', 'among', 'within')
PRED_KEYS = ('K1', 'K2', 'K3', 'K4', 'K1p', 'K2p', 'W1', 'W2', 'T1', 'T2')


def fit_descriptions(var: pd.DataFrame, means: pd.Series) -> dict:
    """The two descriptions of the mean APL effect at size 1 (`var` per setting: between, among, within, total; `means`
    the saved means, same index): (i) a + b x total, (ii) a + b1 between + b2 among + b3 within, by least squares; their
    coefficients, residuals and root mean square residuals."""
    idx = [s for s in means.index if s in var.index and np.isfinite(means[s])]
    y = means.reindex(idx).to_numpy(dtype=float)
    X1 = np.column_stack([np.ones(len(idx)), var.loc[idx, 'total'].to_numpy(dtype=float)])
    X2 = np.column_stack([np.ones(len(idx))] + [var.loc[idx, c].to_numpy(dtype=float) for c in PARTS3])
    out = {}
    for name, X, labels in (('i', X1, ('intercept', 'total')), ('ii', X2, ('intercept',) + PARTS3)):
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        res = y - X @ beta
        out[name] = dict(coef={l: float(b) for l, b in zip(labels, beta)}, residuals=dict(zip(idx, map(float, res))),
                         rms=float(np.sqrt(np.mean(res ** 2))), n_points=len(idx), n_coef=len(labels))
    return out


def predict(desc: dict, row: dict) -> Tuple[float, float]:
    """The two descriptions' predictions of the mean APL effect for a component's variance parts."""
    c1, c2 = desc['i']['coef'], desc['ii']['coef']
    p1 = c1['intercept'] + c1['total'] * row['total']
    p2 = c2['intercept'] + sum(c2[c] * row[c] for c in PARTS3)
    return float(p1), float(p2)


def class_between(slope_share: float, ineq_share: float) -> str:
    """Part C's class between KCs: 'inequality' (its share 0.6 or more, the slope's below 0.25), 'slope' (the other way
    round), 'neither' (both below 0.25), 'both' (any other pair); 'not read' when a share is missing."""
    if not (np.isfinite(slope_share) and np.isfinite(ineq_share)):
        return 'not read'
    if ineq_share >= CLASS_HI and slope_share < CLASS_LO:
        return 'inequality'
    if slope_share >= CLASS_HI and ineq_share < CLASS_LO:
        return 'slope'
    if slope_share < CLASS_LO and ineq_share < CLASS_LO:
        return 'neither'
    return 'both'


def class_within(placement_share: float, amounts_share: float) -> str:
    """Part C's class within a KC: 'placement', 'amounts', 'neither', 'both' (the same thresholds as `class_between`)."""
    c = class_between(amounts_share, placement_share)
    return {'inequality': 'placement', 'slope': 'amounts'}.get(c, c)


def shares_between(U: float, O: float, K1: float, K2: float) -> dict:
    """The gap U - O and the shares of it that K1 (the slope set to proportional) and K2 (the inequality removed)
    recover, with the class."""
    gap = U - O
    ok = np.isfinite(gap) and gap > 0
    s = (K1 - O) / gap if ok else np.nan
    i = (K2 - O) / gap if ok else np.nan
    return dict(U=U, O=O, K1=K1, K2=K2, gap=gap, slope_share=s, ineq_share=i,
                cls=class_between(s, i) if ok else ('no gap' if np.isfinite(gap) else 'not read'))


def shares_within(U: float, D: float, W1: float, W2: float) -> dict:
    """The gap U - D and the shares of it that W1 (each KC's on-claw count spread evenly: the placement) and W2 (every
    KC the same amount per claw: the amounts) recover, with the class."""
    gap = U - D
    ok = np.isfinite(gap) and gap > 0
    p = (W1 - D) / gap if ok else np.nan
    a = (W2 - D) / gap if ok else np.nan
    return dict(U=U, D=D, W1=W1, W2=W2, gap=gap, placement_share=p, amounts_share=a,
                cls=class_within(p, a) if ok else ('no gap' if np.isfinite(gap) else 'not read'))


def none_verdict(U: float, O: float, K4: float) -> dict:
    """The KCs with none: (U - K4) / (U - O); half or more: the KCs without an off-claw synapse carry the loss; a
    quarter or less: presence alone holds the mean."""
    gap = U - O
    sh = (U - K4) / gap if np.isfinite(gap) and gap > 0 else np.nan
    if not np.isfinite(sh):
        v = 'not read'
    elif sh >= NONE_HI:
        v = 'the KCs without an off-claw synapse carry the loss'
    elif sh <= NONE_LO:
        v = 'presence alone holds the mean'
    else:
        v = 'between'
    return dict(U=U, O=O, K4=K4, lost_share=sh, verdict=v)


def evenness_bound(U: float, O: float, K2: float, K3: float, cv_O: float, cv_K3: float, loss: float = BOUND_LOSS) -> dict:
    """How even is even enough: read only if K2 recovers a quarter of the gap U - O or more. r3 = (K3 - O) / (K2 - O),
    the share that K3 (half the inequality) recovers of what K2 (none) recovers (0.75 if the loss goes with the square of
    the inequality, 0.5 if in proportion to it); the loss L(cv) = K2 - mean as a power law of the CV among KCs with the
    same number of claws through (cv_K3, K2 - K3) and (cv_O, K2 - O), so its exponent is log((K2 - O) / (K2 - K3)) /
    log(cv_O / cv_K3); 'cv_at_loss' the CV at which `loss` of the mean is lost ('extrapolated' when K2 - O < `loss`)."""
    gap = U - O
    k2_share = (K2 - O) / gap if np.isfinite(gap) and gap > 0 else np.nan
    out = dict(U=U, O=O, K2=K2, K3=K3, cv_O=cv_O, cv_K3=cv_K3, k2_share=k2_share,
               applies=bool(np.isfinite(k2_share) and k2_share >= BOUND_MIN_SHARE), loss=loss)
    LO, L3 = K2 - O, K2 - K3
    out['r3'] = float((K3 - O) / LO) if np.isfinite(LO) and LO > 0 else np.nan
    ok = (np.isfinite(LO) and np.isfinite(L3) and L3 > 0 and LO > L3 and np.isfinite(cv_O) and np.isfinite(cv_K3)
          and cv_O > cv_K3 > 0)
    if ok:
        g = float(np.log(LO / L3) / np.log(cv_O / cv_K3))
        out.update(exponent=g, cv_at_loss=float(cv_O * (loss / LO) ** (1.0 / g)), extrapolated=bool(LO < loss))
    else:
        out.update(exponent=np.nan, cv_at_loss=np.nan, extrapolated=False,
                   note='no power law: K3 not between O and K2, or the CVs not ordered')
    return out


def by_wiring(classes: Dict[str, str], s_class: str, min_same: int = BY_WIRING_MIN) -> dict:
    """The class holds by wiring if at least `min_same` of S's wirings give S's class from their own values."""
    same = [w for w, c in classes.items() if c == s_class]
    return dict(per_wiring=dict(classes), n_same=len(same), holds_by_wiring=bool(len(same) >= min_same),
                reported_as='by wiring' if len(same) >= min_same else "S's alone")


def _e96_readings() -> pd.DataFrame:
    p = H.results_dir(E96.EXPERIMENT) / 'readings.csv'
    if not p.exists():
        raise RuntimeError(f"E96's readings.csv is not in {p.parent}")
    return pd.read_csv(p)


def saved_means(rd96: pd.DataFrame, unit: str = 'S', readout: str = 'cal') -> pd.Series:
    """The saved corners' mean APL effect (E96's readings.csv: E96's settings and E95's part D settings read beside),
    by A1's setting name (SAVED7) and by variant."""
    g = rd96[(rd96.unit == unit) & (rd96.readout == readout)].drop_duplicates('variant').set_index('variant').mean_effect
    return pd.Series({s: float(g.get(SAVED_SHAPES[s], np.nan)) for s in SAVED7})


def a2_predictions(out: Path, a1: Optional[pd.DataFrame] = None, force: bool = False) -> dict:
    """A2 (module docstring): a2_predictions.csv and a2_predictions.json (the descriptions fitted to the seven saved
    means of S's population under the calibrated readout at size 1, on the real wiring's variance parts; their
    predictions for the ten new settings at size 1 and the classes of part C they predict; written_at, the commits, the
    sha256 of E96's readings.csv). Not rewritten once part B has a row (unless `force`)."""
    import hashlib
    path = out / 'a2_predictions.json'
    if path.exists() and (out / RUNS_NAMES['B']).exists() and not force:
        print(f'E99 A2: part B has rows; the predictions of {_read_json(path).get("written_at")} are kept', flush=True)
        return _read_json(path)
    a1 = pd.read_csv(out / 'a1_shapes.csv') if a1 is None else a1
    var = a1[a1.wiring == 'real'].set_index('setting')
    rd96 = _e96_readings()
    means = saved_means(rd96)
    desc = fit_descriptions(var, means)
    rows = []
    for s in SAVED7:
        r = var.loc[s]
        p1, p2 = predict(desc, r)
        rows.append(dict(setting=s, variant=SAVED_SHAPES[s], kind='saved', **{c: float(r[c]) for c in PARTS3 + ('total',)},
                         loglog_slope_n=float(r['loglog_slope_n']), observed=float(means[s]), pred_i=p1, pred_ii=p2))
    pred = {}
    for s in PRED_KEYS:
        r = var.loc[s]
        p1, p2 = predict(desc, r)
        pred[s] = (p1, p2)
        rows.append(dict(setting=s, variant=new_variant(s), kind='new', **{c: float(r[c]) for c in PARTS3 + ('total',)},
                         loglog_slope_n=float(r['loglog_slope_n']), observed=np.nan, pred_i=p1, pred_ii=p2))
    df = pd.DataFrame(rows)
    U, O, Op, D = means['uniform'], means['offclaw'], means['offclaw_ped'], means['connectome']
    classes = {}
    for j, name in ((0, 'i'), (1, 'ii')):
        p = {s: pred[s][j] for s in PRED_KEYS}
        cv_O = float(var.loc['offclaw', 'within_cv_median'])
        cv_K3 = float(var.loc['K3', 'within_cv_median'])
        classes[name] = dict(
            between=shares_between(U, O, p['K1'], p['K2']),
            replicate_peduncle=shares_between(U, Op, p['K1p'], p['K2p']),
            within=shares_within(U, D, p['W1'], p['W2']),
            none=none_verdict(U, O, p['K4']),
            bound=evenness_bound(U, O, p['K2'], p['K3'], cv_O, cv_K3),
            d_rule_mean={t: dict(mean=p[t], meets_mean=bool(p[t] >= D_MEAN_MIN)) for t in T_SHAPES})
    df.to_csv(out / 'a2_predictions.csv', index=False)
    res = dict(written_at=time.strftime('%Y-%m-%d %H:%M:%S'), host=_host(), commits=dict(H.code_commits()),
               source=f"{H.results_dir(E96.EXPERIMENT) / 'readings.csv'} (unit S, readout cal, mean_effect)",
               source_sha256=hashlib.sha256((H.results_dir(E96.EXPERIMENT) / 'readings.csv').read_bytes()).hexdigest(),
               variances='a1_shapes.csv, the real wiring (the shapes are the same on every wiring: the same claws and counts)',
               saved_means={s: float(means[s]) for s in SAVED7}, descriptions=desc,
               predictions={s: dict(variant=new_variant(s), i=pred[s][0], ii=pred[s][1]) for s in PRED_KEYS},
               classes=classes, flies_mean=FLIES_MEAN, d_mean_min=D_MEAN_MIN)
    _write_json(path, res)
    return res


def predictions_written() -> bool:
    return bool(_read_json(out_dir() / 'a2_predictions.json').get('written_at'))


# ---------------------------------------------------------------------------
# A3: who responds

def _rank_fifths(v: np.ndarray) -> np.ndarray:
    """Fifths 1 (lowest) to 5 of `v` by rank (ties in the rows' order)."""
    r = pd.Series(np.asarray(v, dtype=float)).rank(method='first').to_numpy()
    return np.ceil(r / len(r) * 5).astype(int).clip(1, 5)


def relative_per_claw(x: np.ndarray, kc: pd.DataFrame, claws: pd.DataFrame) -> np.ndarray:
    """Each KC's component per claw (its total over its claws) relative to its claw-count group's (the group's total over
    its claws), in `kc`'s order."""
    kc_of = claws[KC].to_numpy()
    ids = kc[KC].to_numpy()
    tot = pd.Series(kc_totals(x, kc_of, ids), index=ids)
    n = kc.set_index(KC)['n_claws'].astype(float)
    grp = group_of_kc(kc)
    mu = group_means(tot, n, grp)
    with np.errstate(divide='ignore', invalid='ignore'):
        return ((tot / n) / mu).to_numpy(dtype=float)


def _claw_class(n: np.ndarray) -> np.ndarray:
    return E95._claw_group(np.asarray(n, dtype=float))


def kc_response_frame(on_counts: pd.DataFrame, off_counts: pd.DataFrame, kc: pd.DataFrame, claws: pd.DataFrame,
                      x: Optional[np.ndarray]) -> pd.DataFrame:
    """Per KC (`kc`'s order): its type, claws, claw class, the odors it responds to with APL on and off (more than
    RAW_CRIT spikes), the fifth of its calyx off-claw count per claw relative to its group (K1's value) and, with `x`
    (the run's component per claw), the fifth of the component per claw relative to its group (none when that is the
    same for every KC)."""
    ids = kc[KC].to_numpy()

    def odors(c: pd.DataFrame) -> np.ndarray:
        got = np.asarray(c.index.get_level_values(KC))
        s = pd.Series((c.to_numpy(dtype=float) > RAW_CRIT).sum(axis=1), index=got)
        if s.index.has_duplicates or not set(ids) <= set(s.index):
            raise RuntimeError("the run's KCs are not the table's")
        return s.reindex(ids).to_numpy(dtype=float)
    f = pd.DataFrame({KC: ids, 'kc_type': kc['kc_type'].astype(str).to_numpy() if 'kc_type' in kc else '',
                      'n_claws': kc['n_claws'].to_numpy(dtype=float)})
    f['claw_class'] = _claw_class(f['n_claws'].to_numpy())
    f['odors_on'] = odors(on_counts)
    f['odors_off'] = odors(off_counts)
    k1 = shape99('K1', kc, claws, 'offclaw')                 # (a_k / n_k) / mu_G on each claw: the probe's ranking
    f['rel_offclaw'] = kc_totals(k1, claws[KC].to_numpy(), ids) / f['n_claws'].to_numpy()
    f['fifth_offclaw'] = _rank_fifths(f['rel_offclaw'].to_numpy())
    if x is not None:
        rel = relative_per_claw(x, kc, claws)
        f['rel_own'] = rel
        f['fifth_own'] = _rank_fifths(rel) if np.nanmax(rel) - np.nanmin(rel) > 1e-9 else 0
    return f


def summarize_responses(f: pd.DataFrame, labels: dict) -> List[dict]:
    """A3's rows of a per-KC frame (`kc_response_frame`, one wiring or several pooled): by fifth (of the off-claw count
    and of the run's own component), by claw class and by KC type, the KCs, the odors per KC with APL on and off and the
    share of the responses with APL on and off."""
    rows = []
    t_on, t_off = float(f.odors_on.sum()), float(f.odors_off.sum())
    keys = [('fifth_offclaw', 'fifth'), ('fifth_own', 'fifth'), ('claw_class', 'claw class'), ('kc_type', 'KC type')]
    for col, kind in keys:
        if col not in f or (col == 'fifth_own' and not (f[col] > 0).any()):
            continue
        for g, q in f.groupby(col, sort=True):
            rows.append(dict(labels, by=col, kind=kind, group=str(g), n_kcs=int(len(q)), odors_on=float(q.odors_on.mean()),
                             odors_off=float(q.odors_off.mean()), share_on=float(q.odors_on.sum() / t_on) if t_on else np.nan,
                             share_off=float(q.odors_off.sum() / t_off) if t_off else np.nan))
    rows.append(dict(labels, by='all', kind='all', group='all', n_kcs=int(len(f)), odors_on=float(f.odors_on.mean()),
                     odors_off=float(f.odors_off.mean()), share_on=1.0, share_off=1.0))
    return rows


def claw_extremes(stage_on: pd.DataFrame, stage_off: pd.DataFrame, x: np.ndarray, claws: pd.DataFrame) -> pd.DataFrame:
    """Per KC with two or more claws and an uneven component `x` (per claw, the claws' order): its most and least weighted
    claw (the first of ties) and their window state with APL on (`stage_on`, claws x odors) and off: the mean over odors
    and the share of odors above 0."""
    ix = pd.MultiIndex.from_frame(claws[[KC, 'claw_id']])

    def align(st: pd.DataFrame) -> np.ndarray:
        got = pd.MultiIndex.from_arrays([st.index.get_level_values(KC), st.index.get_level_values('claw_id')])
        pos = got.get_indexer(ix)
        if (pos < 0).any():
            raise RuntimeError("the stage's claws are not the wiring's")
        return st.to_numpy(dtype=float)[pos]
    on, off = align(stage_on), align(stage_off)
    s = pd.DataFrame(dict(k=claws[KC].to_numpy(), x=np.asarray(x, dtype=float), i=np.arange(len(claws))))
    g = s.groupby('k')
    rows = []
    for k, q in g:
        if len(q) < 2 or np.ptp(q.x.to_numpy()) <= 1e-12:
            continue
        hi = int(q.i.to_numpy()[np.argmax(q.x.to_numpy())])
        lo = int(q.i.to_numpy()[np.argmin(q.x.to_numpy())])
        for which, i in (('most', hi), ('least', lo)):
            rows.append(dict(kc_id=k, which=which, n_claws=len(q), weight=float(s.x.iloc[i]), mean_on=float(on[i].mean()),
                             mean_off=float(off[i].mean()), active_on=float((on[i] > 0).mean()),
                             active_off=float((off[i] > 0).mean())))
    return pd.DataFrame(rows)


def summarize_claws(c: pd.DataFrame, labels: dict) -> List[dict]:
    rows = []
    if c.empty:
        return rows
    c = c.assign(claw_class=_claw_class(c.n_claws.to_numpy()))
    for (cls, which), q in list(c.groupby(['claw_class', 'which'])) + [(('all', w), q) for w, q in c.groupby('which')]:
        rows.append(dict(labels, claw_class=str(cls), which=which, n_kcs=int(len(q)), weight=float(q.weight.mean()),
                         mean_on=float(q.mean_on.mean()), mean_off=float(q.mean_off.mean()),
                         kept_share=float(q.mean_on.sum() / q.mean_off.sum()) if q.mean_off.sum() > 0 else np.nan,
                         active_on=float(q.active_on.mean()), active_off=float(q.active_off.mean())))
    return rows


A3_RUNS = {'D1_uniform_d1': ('E95', 'uniform'), 'D2_connectome_d1': ('E95', 'connectome'),
           'B1_offclaw_d1': ('E96', 'offclaw'), 'B2_offclawped_d1': ('E96', 'offclaw_ped'), 'B6_allsyn_d1': ('E96', 'allsyn')}
A3_CLAWS = ('D2_connectome_d1', 'B6_allsyn_d1')


def who_responds(specs: Sequence[Tuple[str, str, str]], wirings: Sequence[str] = S_WIRINGS,
                 claw_variants: Sequence[str] = ()) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3's reading of runs: `specs` (variant, source, shape: SAVED_SHAPES' name or a variant of this experiment) on each
    wiring and S's population (the four wirings' KCs pooled); the claw stage of `claw_variants`."""
    rows, crows = [], []
    for v, source, shape in specs:
        frames, cl = [], []
        for w in wirings:
            kc, claws = cached_tables(w)
            tag = run_tag(v, w)
            o = load_any99(source, tag, w)
            off = E83.saved_fit(w)['saved']['off'][1]
            x = saved_shape(shape, kc, claws, w) if shape in SAVED_SHAPES else setting_shape(shape, kc, claws)
            f = kc_response_frame(o['on'][1], off, kc, claws, None if shape == 'uniform' else x)
            frames.append(f)
            rows += summarize_responses(f, dict(variant=v, source=source, unit=w))
            if v in claw_variants:
                st = {'E99': load_stage99, 'E96': E96.load_stage96, 'E95': E95.load_stage95}[source](tag)
                if st is None:
                    raise RuntimeError(f'{tag}: no stage matrices')
                odors = E94._odors(st['claw_win'])
                st_off = E94.fit_stage(w, 'off')['claw_win']
                c = claw_extremes(st['claw_win'][odors], st_off[odors], x, claws)
                cl.append(c)
                crows += summarize_claws(c, dict(variant=v, source=source, unit=w))
        if len(frames) == len(wirings) and tuple(wirings) == tuple(S_WIRINGS):
            pooled = pd.concat(frames, ignore_index=True)
            rows += summarize_responses(pooled, dict(variant=v, source=source, unit='S'))
            if cl:
                crows += summarize_claws(pd.concat(cl, ignore_index=True), dict(variant=v, source=source, unit='S'))
    return pd.DataFrame(rows), pd.DataFrame(crows)


def a3_who_responds(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3: a3_fifths.csv and a3_claws.csv for the saved runs at size 1 (A3_RUNS; the claw stage of D2 and B6)."""
    specs = [(v, src, shape) for v, (src, shape) in A3_RUNS.items()]
    f, c = who_responds(specs, S_WIRINGS, A3_CLAWS)
    f.to_csv(out / 'a3_fifths.csv', index=False)
    c.to_csv(out / 'a3_claws.csv', index=False)
    return f, c


# ---------------------------------------------------------------------------
# A4: the bound's yardstick

A4_KINDS = {'calyx_off_claws': 'no_claw', 'calyx_on_claws': 'kept', 'peduncle': 'peduncle',
            'off_claws_and_peduncle': 'no_claw+peduncle'}


def a4_yardstick(out: Path) -> pd.DataFrame:
    """A4: a4_yardstick.csv, the CV among KCs with the same number of claws (E96's `_within`) of the off-claw, on-claw and
    peduncle counts in the four connectomes' derived tables (E96's a2_connectomes_kc.csv, their own claw counts) beside the
    model's (the real wiring's counts, Pratyush's claws)."""
    p = H.results_dir(E96.EXPERIMENT) / 'a2_connectomes_kc.csv'
    rows = []
    if p.exists():
        d = pd.read_csv(p)
        for r in d.itertuples():
            rows.append(dict(source='connectome', dataset=r.dataset, apl=r.apl, kind=r.kind, n_kcs=int(r.n_kcs),
                             within_cv_min=r.within_cv_min, within_cv_median=r.within_cv_median,
                             within_cv_max=r.within_cv_max, loglog_slope_n=r.loglog_slope_n, kcs_none_share=r.kcs_none_share))
    kc, _ = cached_tables('real')
    n = kc['n_claws'].to_numpy(dtype=float)
    for kind, col in A4_KINDS.items():
        v = E96._kind_values(kc, col)
        w = E96._within(n, v)
        rows.append(dict(source='model', dataset='hemibrain (the model\'s claws and counts)', apl='APL_R', kind=kind,
                         n_kcs=int(len(kc)), within_cv_min=w['within_cv_min'], within_cv_median=w['within_cv_median'],
                         within_cv_max=w['within_cv_max'], loglog_slope_n=E96._loglog(n, v),
                         kcs_none_share=float((v == 0).mean())))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a4_yardstick.csv', index=False)
    return df


def part_a(out: Path, force_a2: bool = False) -> dict:
    """A1 to A4 and a_summary.json."""
    t0 = time.time()
    a1 = a1_shapes(out)
    a2 = a2_predictions(out, a1, force=force_a2)
    f, c = a3_who_responds(out)
    a4 = a4_yardstick(out)
    real = a1[a1.wiring == 'real'].set_index('setting')
    s = f[(f.unit == 'S')]
    summ = dict(a1_real={k: {c_: float(real.loc[k, c_]) for c_ in ('between', 'among', 'within', 'total', 'loglog_slope_n',
                                                                     'within_cv_median', 'claws_zero_share', 'kcs_zero_share')}
                         for k in real.index},
                a1_same_on_s=bool(all(np.allclose(a1[a1.wiring == w].set_index('setting')[list(PARTS3)].to_numpy(),
                                                  real[list(PARTS3)].to_numpy(), atol=1e-12) for w in S_WIRINGS)),
                a2=dict(written_at=a2.get('written_at'), rms={k: a2['descriptions'][k]['rms'] for k in ('i', 'ii')},
                        coef={k: a2['descriptions'][k]['coef'] for k in ('i', 'ii')},
                        predictions=a2['predictions'],
                        classes={k: {kk: (vv.get('cls') or vv.get('verdict') or vv) for kk, vv in a2['classes'][k].items()}
                                 for k in ('i', 'ii')}),
                a3_S_fifths={v: q[q.by == 'fifth_offclaw'][['group', 'odors_on', 'share_on']].to_dict('records')
                             for v, q in s.groupby('variant')},
                a3_S_claw_class={v: q[q.by == 'claw_class'][['group', 'odors_on', 'odors_off']].to_dict('records')
                                 for v, q in s.groupby('variant')},
                a4=a4[['source', 'dataset', 'apl', 'kind', 'within_cv_median']].to_dict('records'),
                seconds=time.time() - t0)
    _write_json(out / 'a_summary.json', summ)
    print(f'E99 A in {time.time() - t0:.0f} s', flush=True)
    return summ


# ---------------------------------------------------------------------------
# part 0: the identities

def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('passes'))


def shape_identities(kc: pd.DataFrame, claws: pd.DataFrame, rtol: float = 1e-12) -> dict:
    """The shapes' identities on one wiring's tables (module docstring, part 0; no simulation)."""
    from experiments import E10_joint_fit as J
    kc_of = claws[KC].to_numpy()
    ids = kc[KC].to_numpy()
    grp = group_of_kc(kc)
    gclaw = grp.reindex(kc_of).to_numpy()
    n = kc.set_index(KC)['n_claws'].astype(float)
    sh = {f'{s}_{c}': shape99(s, kc, claws, c) for s in K_SHAPES for c in COUNTS}
    sh.update({s: shape99(s, kc, claws) for s in W_SHAPES + T_SHAPES})
    res = dict(n_kcs=int(len(kc)), n_claws=int(len(claws)), groups={int(g): int((grp == g).sum()) for g in sorted(set(grp))})
    res['means_1'] = bool(all(abs(x.mean() - 1) < 1e-12 for x in sh.values()))
    k = kc.set_index(KC)
    off = {c: J.aplkc_flat_shape(c, kc_of, J.aplkc_flat_offclaw(c, k)) for c in COUNTS}
    kept = claws['kept'].to_numpy(dtype=float)
    conn = kept / kept.mean()
    res['k1_k2_is_offclaw'] = {c: bool(np.allclose(sh[f'K1_{c}'] * sh[f'K2_{c}'], off[c], rtol=rtol, atol=0)) for c in COUNTS}
    res['w1_w2_is_connectome'] = bool(np.allclose(sh['W1'] * sh['W2'], conn, rtol=rtol, atol=0))
    if 'wnorm' in claws:
        w = claws['wnorm'].to_numpy(dtype=float)
        res['w1_w2_is_wnorm_shape'] = bool(np.allclose(sh['W1'] * sh['W2'], w / w.mean(), rtol=1e-10, atol=0))
    if 'comp_offclaw' in claws:
        res['offclaw_is_e96_cache'] = bool(np.allclose(off['offclaw'], claws['comp_offclaw'].to_numpy(), rtol=rtol, atol=0)
                                           and np.allclose(off['offclaw_ped'], claws['comp_offclaw_ped'].to_numpy(),
                                                           rtol=rtol, atol=0))
    gm = {c: pd.Series(sh[f'K1_{c}']).groupby(gclaw).mean() for c in COUNTS}
    res['k1_group_means_1'] = {c: bool(np.allclose(gm[c].to_numpy(), 1.0, rtol=0, atol=1e-12)) for c in COUNTS}

    def one_value(x, by):
        s = pd.Series(np.asarray(x, dtype=float)).groupby(np.asarray(by))
        return bool(((s.max() - s.min()) <= 1e-12).all())
    res['k2_one_value_in_group'] = {c: one_value(sh[f'K2_{c}'], gclaw) for c in COUNTS}
    res['k2_k4_one_value_in_kc'] = {c: bool(one_value(sh[f'K2_{c}'], kc_of) and one_value(sh[f'K4_{c}'], kc_of)) for c in COUNTS}
    res['k4_one_nonzero_value'] = {c: bool(np.ptp(sh[f'K4_{c}'][sh[f'K4_{c}'] > 0]) <= 1e-12) for c in COUNTS}
    res['k4_zero_kcs_share'] = {c: float((offclaw_amounts(kc, c) == 0).mean()) for c in COUNTS}
    res['k3_all_is_offclaw'] = {c: bool(np.allclose(shape99('K3', kc, claws, c, frac=1.0), off[c], rtol=rtol, atol=0))
                                for c in COUNTS}
    res['k3_none_is_k2'] = {c: bool(np.allclose(shape99('K3', kc, claws, c, frac=0.0), sh[f'K2_{c}'], rtol=rtol, atol=0))
                            for c in COUNTS}
    # the sums per KC: W1 c_k, T1 c_k + a^c_k + a^p_k, T2 c_k + a^c_k (in the components' unit: the count per claw's mean)
    c_k = pd.Series(kept).groupby(kc_of).sum().reindex(ids).to_numpy()
    tot = {'W1': c_k, 'T1': c_k + k['no_claw'].reindex(ids).to_numpy() + k['peduncle'].reindex(ids).to_numpy(),
           'T2': c_k + k['no_claw'].reindex(ids).to_numpy()}
    res['sums_are_counts'] = {}
    for s, t in tot.items():
        got = kc_totals(sh[s], kc_of, ids)
        res['sums_are_counts'][s] = bool(np.allclose(got * (t.sum() / len(claws)), t, rtol=rtol, atol=1e-9))
    res['w2_sum_is_n'] = bool(np.allclose(kc_totals(sh['W2'], kc_of, ids), n.reindex(ids).to_numpy(), rtol=rtol, atol=0))
    res['w2_kcs_without_kept'] = int((c_k == 0).sum())
    # T1 and T2 are E10's new modes (the fit's components) exactly
    per_claw = pd.Series(kept, index=pd.MultiIndex.from_frame(claws[[KC, 'claw_id']]))
    wn = pd.Series(conn, index=per_claw.index)
    for s in T_SHAPES:
        comp = J.aplkc_flat_component(wn, T_MODE[s], (k, per_claw))
        res[f'{s}_is_e10'] = bool(np.array_equal(comp, float(wn.to_numpy().mean()) * sh[s]))
    flat = [res['means_1'], res['w1_w2_is_connectome'], res['w2_sum_is_n'], res['T1_is_e10'], res['T2_is_e10']]
    flat += [res.get('w1_w2_is_wnorm_shape', True), res.get('offclaw_is_e96_cache', True)]
    for key in ('k1_k2_is_offclaw', 'k1_group_means_1', 'k2_one_value_in_group', 'k2_k4_one_value_in_kc',
                'k4_one_nonzero_value', 'k3_all_is_offclaw', 'k3_none_is_k2', 'sums_are_counts'):
        flat += list(res[key].values())
    res['passes'] = bool(all(flat))
    return res


def tables_identity(label: str) -> dict:
    """E96's cached per-KC and per-claw tables against the wiring's live counts (the same KCs, claws and counts)."""
    kc0, cl0 = cached_tables(label)
    kc1, cl1 = count_tables(label)
    k0, k1 = kc0.set_index(KC), kc1.set_index(KC)
    cols = ['n_claws', 'kept', 'no_claw', 'peduncle']
    same_kcs = set(k0.index) == set(k1.index)
    res = dict(same_kcs=bool(same_kcs),
               same_counts=bool(same_kcs and np.array_equal(k0[cols].to_numpy(dtype=float),
                                                            k1.reindex(k0.index)[cols].to_numpy(dtype=float))),
               same_claws=bool(np.array_equal(cl0[[KC, 'claw_id']].to_numpy(), cl1[[KC, 'claw_id']].to_numpy())),
               same_kept=bool(np.array_equal(cl0['kept'].to_numpy(dtype=float), cl1['kept'].to_numpy(dtype=float))))
    res['passes'] = bool(all(res.values()))
    return res


IDENTITY_RUNS = {'uniform_real': ('D1_uniform_d1', 'E95'), 'connectome_real': ('D2_connectome_d1', 'E95'),
                 'offclaw_real': ('B1_offclaw_d1', 'E96'), 'allsyn_real': ('B6_allsyn_d1', 'E96')}


def _saved_row(source: str, tag: str) -> pd.Series:
    p = {'E95': H.results_dir(E95.EXPERIMENT) / 'runs_D.csv', 'E96': H.results_dir(E96.EXPERIMENT) / 'runs_B.csv'}[source]
    return pd.read_csv(p).drop_duplicates('tag', keep='last').set_index('tag').loc[tag]


def _saved_weights(source: str, tag: str) -> Optional[np.ndarray]:
    dirs = {'E95': E95._dirs('runs'), 'E96': E96._dirs('runs')}[source]
    for d in dirs:
        p = d / f'{tag}_wAPLKC.parquet'
        if p.exists():
            return pd.read_parquet(p)['wAPLKC'].to_numpy(dtype=float)
    return None


def part_0(out: Path, wirings: Sequence[str] = S_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per key); raises if one fails."""
    import gc
    cc = _commits()
    path = out / 'identity.json'
    res = _read_json(path)
    res['commits'] = cc
    for w in wirings:
        key = f'shapes_{w}'
        if not res.get(key, {}).get('passes'):
            t0 = time.time()
            kc, claws = cached_tables(w)
            res[key] = dict(shape_identities(kc, claws), tables=tables_identity(w), seconds=time.time() - t0)
            res[key]['passes'] = bool(res[key]['passes'] and res[key]['tables']['passes'])
            print(f"E99 0 {key}: {json.dumps(res[key], default=str)[:1500]}", flush=True)
            _write_json(path, res)
    if 'real' in wirings and not all(res.get(k, {}).get('passes') for k in IDENTITY_RUNS):
        m = E95.StaticModel('real', check=False)
        L = E95.kc_losses('real', m.saved['on'][1].index)
        unit = E95.d_unit(m, L)
        for key, (v, source) in IDENTITY_RUNS.items():
            if res.get(key, {}).get('passes'):
                continue
            t0 = time.time()
            tag = run_tag(v, 'real')
            pv = parse_variant(v)
            wk, _, _ = point_weights99(m, pv, unit)
            rd = E96.reading96(m, pv, wk, label=f'identity {key}')
            o = load_any99(source, tag, 'real')
            dd = E96._counts_diff(rd['on'], o['on'])
            row = _saved_row(source, tag)
            f0 = float(row['apl_factor'])
            ws = _saved_weights(source, tag)
            same_w = None if ws is None else bool(np.array_equal(wk.to_numpy(dtype=float), ws))
            ok = dict(variant=v, source=source, counts=dd, apl_factor=rd['factor'], apl_factor_saved=f0,
                      factor_rel_diff=abs(rd['factor'] - f0) / f0, unit=unit, unit_saved=float(row['component_unit']),
                      same_weights_as_saved=same_w, n_sims=rd['n_sims'], seconds=time.time() - t0)
            ok['passes'] = bool(dd['exact'] and ok['factor_rel_diff'] <= 1e-12 and same_w is not False
                                and abs(unit - ok['unit_saved']) <= 1e-12 * unit)
            res[key] = ok
            print(f'E99 0 {key}: {ok}', flush=True)
            _write_json(path, res)
        del m
        gc.collect()
    keys = [f'shapes_{w}' for w in S_WIRINGS] + list(IDENTITY_RUNS)
    have_all = all(k in res for k in keys)
    failed = [k for k in keys if k in res and not res[k].get('passes')]
    res.update(passes=bool(have_all and not failed), complete=bool(have_all), failed=failed)
    _write_json(path, res)
    if failed:
        raise RuntimeError(f'E99 part 0: the identity fails for {failed} (identity.json); the experiment stops here')
    print(f"E99 0: {'passes' if res['passes'] else 'so far exact (not all keys yet)'}", flush=True)
    return res


# ---------------------------------------------------------------------------
# part B (and C's runs): the points

def point_weights99(model: E95.StaticModel, pv: dict, unit: float
                    ) -> Tuple[Optional[pd.Series], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """A point's APL->KC weights (None: the fit's, size 0), its component per KC (the outputs' KCs: n_claws, the amount
    the component spreads or places, the component's sum over the KC's claws in unscaled weight units, the family's sum)
    and per claw (kc_id, claw_id, the shape with mean 1). This experiment's settings from `shape99` of the wiring's live
    counts; the saved corners (E95's uniform and connectome, E96's modes) from E96's `component_values`."""
    w = model.label
    if not pv['size'] > 0:
        return None, None, None
    wn = E96.wiring_wnorm(w)
    kc, claws = count_tables(w)
    if not np.array_equal(claws[[KC, 'claw_id']].to_numpy(), wn.index.to_frame(index=False)[[KC, 'claw_id']].to_numpy()):
        raise RuntimeError(f"{w}: the counts' claws are not the weights' claws")
    if pv['part'] == 'B':
        shape = shape99(pv['shape'], kc, claws, pv['counts'] or 'offclaw')
        comp = float(wn.to_numpy(dtype=float).mean()) * shape
        if pv['shape'] in T_SHAPES:                       # the fit's component (E10's mode) exactly
            from experiments import E10_joint_fit as J
            per_kc, per_claw = J.aplkc_flat_counts(wn.index)
            c10 = J.aplkc_flat_component(wn, T_MODE[pv['shape']], (per_kc, per_claw))
            if not np.array_equal(c10, comp):
                raise RuntimeError(f"{pv['variant']}: the component is not E10's {T_MODE[pv['shape']]}")
    elif pv['part'] in ('e95', 'e96') and pv['permute'] is None:
        counts = E96.wiring_counts(w) if pv['mode'] in E96.MODES else None
        comp = E96.component_values(pv['mode'], wn, counts)
        shape = comp / float(wn.to_numpy(dtype=float).mean())
    else:
        raise ValueError(f"{pv['variant']}: not run through this experiment's path")
    wk = E96.e96_weights(model.wAPLKC, wn, comp, pv['size'], unit)
    kc_of = np.asarray(wn.index.get_level_values(KC))
    ids = np.asarray(model.saved['on'][1].index.get_level_values(KC))
    add = pd.Series(wk.to_numpy(dtype=float) - model.wAPLKC.to_numpy(dtype=float)).groupby(kc_of).sum().reindex(ids)
    k = kc.set_index(KC).reindex(ids)
    frame = pd.DataFrame({KC: ids, 'n_claws': pd.Series(kc_of).value_counts().reindex(ids).to_numpy(),
                          'component_sum': add.to_numpy(dtype=float),
                          'w_fam_sum': pd.Series(model.wAPLKC.to_numpy(dtype=float)).groupby(kc_of).sum().reindex(ids).to_numpy(),
                          'shape_sum': kc_totals(shape, kc_of, ids)})
    for c in ('kept', 'no_claw', 'peduncle'):
        frame[c] = k[c].to_numpy(dtype=float)
    per_claw = pd.DataFrame({KC: kc_of, 'claw_id': np.asarray(wn.index.get_level_values('claw_id')), 'shape': shape})
    return wk, frame, per_claw


def _dirs(sub: str = 'runs') -> List[Path]:
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [H.outputs_dir(EXPERIMENT) / sub, arc / sub, arc / 'from_hal' / sub]


STAGE_KEYS = E94.STAGE_KEYS


def save_run99(tag: str, o: Dict[str, tuple], sub: str = 'runs', stages: Optional[dict] = None,
               drive: Optional[pd.DataFrame] = None, loss: Optional[pd.DataFrame] = None,
               weights: Optional[pd.Series] = None, component: Optional[pd.DataFrame] = None,
               component_claw: Optional[pd.DataFrame] = None) -> None:
    """A run's outputs (APL on, and the saved APL-off run beside), its stage matrices with APL on, each KC's summed claw
    drive, E92's loss of drive at the claws, the APL->KC weights per claw, the component per KC and per claw, under
    outputs/E99_apl_weight_evenness/<sub>/ (E96's `save_run96` layout)."""
    d = local_dir(sub)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')
    if stages is not None:
        for k in STAGE_KEYS:
            x = stages[k].copy()
            x.columns = [str(c) for c in x.columns]
            if isinstance(x.index, pd.MultiIndex):
                x.index = x.index.set_names([str(n) for n in x.index.names])
            x.to_parquet(d / f'{tag}_on_{k}.parquet')
        stages['apl'].to_csv(d / f'{tag}_on_apl.csv')
    for name, x in (('drive', drive), ('loss', loss)):
        if x is not None:
            x = x.copy()
            x.columns = [str(c) for c in x.columns]
            x.to_parquet(d / f'{tag}_{name}.parquet')
    if weights is not None:
        weights.rename('wAPLKC').reset_index(drop=True).to_frame().to_parquet(d / f'{tag}_wAPLKC.parquet')
    if component is not None:
        component.to_parquet(d / f'{tag}_component.parquet')
    if component_claw is not None:
        component_claw.to_parquet(d / f'{tag}_component_claw.parquet')


def load_run99(tag: str, sub: str = 'runs') -> Optional[Dict[str, tuple]]:
    for d in _dirs(sub):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def load_stage99(tag: str, sub: str = 'runs') -> Optional[dict]:
    for d in _dirs(sub):
        if (d / f'{tag}_on_claw_win.parquet').exists():
            st = {k: pd.read_parquet(d / f'{tag}_on_{k}.parquet') for k in STAGE_KEYS}
            st['apl'] = pd.read_csv(d / f'{tag}_on_apl.csv', index_col=0)
            return st
    return None


def _append(path: Path, row: dict) -> None:
    E83._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def run_point99(model: E95.StaticModel, v: str, part: str = '', cc: Optional[dict] = None, sub: str = 'runs'
                ) -> Tuple[dict, dict]:
    """One variant on a wiring at its fit (E96's `reading96` with `point_weights99`), saved with its stage matrices, drive,
    loss, APL->KC weights and component per KC and per claw; the row as E96's `run_point` writes it, with the shape and
    counts. Returns (the run's row, the reading)."""
    w = model.label
    tag = run_tag(v, w)
    pv = parse_variant(v)
    if pv['part'] == REF:
        raise ValueError('the reference is the saved fit; it is not run')
    cc = H.code_commits() if cc is None else cc
    kc_index = model.saved['on'][1].index
    L = E95.kc_losses(w, kc_index)
    unit = E95.d_unit(model, L)
    wk, comp_frame, comp_claw = point_weights99(model, pv, unit)
    res = E96.reading96(model, pv, wk, label=tag)
    o = {'on': res['on'], 'off': res['off']}
    st = res['stages']
    odors = E94._odors(st['bouton_win'])
    drive = E94.kc_drive(model.W, st['bouton_win'][odors].to_numpy(dtype=float), kc_index, odors)
    save_run99(tag, o, sub, stages=st, drive=drive, loss=E96._frame(res['gain']['L'], res['on'][1]), weights=wk,
               component=comp_frame, component_claw=comp_claw)
    resp_off = E89._responders(model.saved['off'][1])
    p = res['params']
    per_kc = E96.wiring_counts(w)[0]
    syn_mode = pv['counts'] if pv['part'] == 'B' and pv['counts'] else (pv['mode'] if pv['mode'] in E96.MODES else None)
    row = dict(part=part, arm=pv['arm'], tag=tag, variant=v, source='E99', screen=pv['screen'], shape=pv['shape'],
               counts=pv['counts'], mode=pv['mode'], size=pv['size'], pool=pv['pool'], dose=pv['dose'],
               construction=pv['construction'], wiring=w, host=_host(), fit_folder=model.folder, fit_tag=model.tag,
               factor=res['factor'], apl_factor=res['factor'], at_bound=res['at_bound'], bracket_lo=res['bracket'][0],
               bracket_hi=res['bracket'][1], n_on=E83._n(res['on'][0]), target_on=model.target_on,
               n_off=E83._n(res['off'][0]), target_off=model.target_off, n_sims=res['n_sims'], seconds=res['seconds'],
               **{k: p[k] for k in E83.PARAM_KEYS}, aplpn_fit_multiple=E94.aplpn_fit_multiple(w),
               aplpn_multiple=E94.aplpn_multiple(w, pv['dose']), wAPLKC_sha256=E89.E80_sha(model.wAPLKC if wk is None else wk),
               K_fit=model.K, K_run=E83.readout_K(o), component_unit=unit,
               synapse_share=E96.synapse_share(syn_mode, per_kc) if syn_mode else np.nan, **res['gain']['row'], **cc)
    row.update(E94.stage_extras(model.W, w, st))
    w_all = (model.wAPLKC if wk is None else wk).to_numpy(dtype=float) * p['wAPLKC_scale']
    w_fam = model.wAPLKC.to_numpy(dtype=float) * p['wAPLKC_scale']
    off_st = E94.fit_stage(w, 'off')
    row.update(E95.d_realism(model.W.to_numpy(dtype=float), off_st['bouton_win'][odors].to_numpy(dtype=float),
                             st['bouton_win'][odors].to_numpy(dtype=float), st['claw_win'][odors].to_numpy(dtype=float),
                             w_fam, w_all - w_fam, E94.kc_of_claw(model.W), resp_off))
    row.update(E83.realism(o))
    print(f"E99 {tag}: APL scale x {res['factor']:.4g}{' (at a bound)' if res['at_bound'] else ''}, {row['n_on']} / "
          f"{model.target_on} responses, {res['n_sims']} sims, {res['seconds']:.0f}s; shares at the boutons / component / "
          f"weighted {row.get('bouton_share_d', np.nan):.3f} / {row.get('comp_share', np.nan):.3f} / "
          f"{row.get('weighted_share', np.nan):.3f}; flag {row.get('realism_flag')!r}", flush=True)
    return row, res


RUNS_NAMES = {'B': 'runs_B.csv', 'C': 'runs_C.csv'}


def _provenance_add(out: Path, part: str, tag: str, prov: dict) -> None:
    path = out / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E99_apl_weight_evenness.py')
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def points(out: Path, variants: Sequence[str], wirings: Sequence[str], runs_name: str, part: str) -> pd.DataFrame:
    """Variants on wirings, each by `run_point99` (resumable: rows in `runs_name`). Needs part 0 to have passed and A2's
    predictions to be written."""
    import gc
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identities must hold before any scored point)')
    if not predictions_written():
        raise RuntimeError("run --part A first (A2's predictions are written before part B's first point)")
    for v in variants:
        parse_variant(v)
    cc = _commits()
    path = out / runs_name
    for w in wirings:
        todo = [v for v in variants if run_tag(v, w) not in _done(path)]
        if not todo:
            continue
        model = E95.StaticModel(w, check=False)
        for v in todo:
            row, _ = run_point99(model, v, part, cc)
            _append(path, row)
            _provenance_add(out, part, run_tag(v, w), dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **cc))
        del model
        gc.collect()
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def part_b(out: Path, blocks: Optional[Sequence[str]] = None, wirings: Sequence[str] = S_WIRINGS) -> pd.DataFrame:
    """The blocks' points (default I then II) on `wirings` (resumable; runs_B.csv)."""
    blocks = list(BLOCKS) if not blocks else list(blocks)
    bad = [b for b in blocks if b not in BLOCKS]
    if bad:
        raise ValueError(f'{bad} are not blocks of part B ({list(BLOCKS)})')
    return points(out, [v for b in blocks for v in BLOCKS[b]], wirings, RUNS_NAMES['B'], 'B')


# ---------------------------------------------------------------------------
# the runs table and the readings

def corner_rows() -> pd.DataFrame:
    """The saved corners on S: E96's runs_B.csv rows of CORNERS['E96'] (source 'E96') and E95's runs_D.csv rows of
    CORNERS['E95'] (source 'E95')."""
    parts = []
    for src, p in (('E96', H.results_dir(E96.EXPERIMENT) / 'runs_B.csv'), ('E95', H.results_dir(E95.EXPERIMENT) / 'runs_D.csv')):
        if not p.exists():
            continue
        d = pd.read_csv(p).drop_duplicates('tag', keep='last')
        d = d[d.variant.isin(CORNERS[src]) & d.wiring.isin(S_WIRINGS)].copy()
        d['source'] = src
        parts.append(d)
    return pd.concat(parts, ignore_index=True, sort=False) if parts else pd.DataFrame()


def runs_table(out: Path) -> pd.DataFrame:
    """Every run of the readings: this experiment's (runs_B / runs_C.csv, the latest row per tag) first, then the saved
    corners and the references (E95's `reference_rows`: the unchanged fits on every wiring, E94's D0 on S); with
    'flagged', 'at_bound' and the parsed setting."""
    parts = []
    for n in RUNS_NAMES.values():
        if (out / n).exists():
            try:
                parts.append(pd.read_csv(out / n).drop_duplicates('tag', keep='last'))
            except pd.errors.EmptyDataError:
                pass
    parts.append(corner_rows())
    parts.append(E95.reference_rows())
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates(['variant', 'wiring'], keep='first').reset_index(drop=True)
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('arm', 'part', 'screen', 'shape', 'counts', 'mode', 'size', 'pool', 'dose', 'permute', 'control', 'construction',
              'equal'):
        df[c] = pv[c].to_numpy()
    df['flagged'] = [E94._flagged(x) for x in df.realism_flag]
    df['at_bound'] = [E94._as_bool(x) for x in df.at_bound]
    return df


def load_any99(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E99' this experiment's, otherwise E96's `load_any` ('E96', 'E95', 'E94', 'saved')."""
    if source == 'E99':
        o = load_run99(tag)
        if o is None:
            raise RuntimeError(f'no saved run E99:{tag}')
        return o
    return E96.load_any(source, tag, wiring)


READOUTS = E95.READOUTS


def _distances(source: str, tag: str, wiring: str, readouts: Sequence[str] = READOUTS) -> Dict[str, Dict[str, pd.Series]]:
    key = ('dist', source, tag, wiring)
    have = _CACHE.setdefault(key, {})
    need = [r for r in readouts if r not in have]
    if need:
        o = load_any99(source, tag, wiring)
        K = E83.saved_fit(wiring)['K']
        for ro in need:
            if ro == 'calK':
                have[ro] = E83.pair_distances(o, 'calK', K)
            elif ro == 'zcell':
                have[ro] = E95.zcell_distances_map(o, 'saturating')
            elif ro == 'zcell_lin':
                have[ro] = E95.zcell_distances_map(o, 'linear')
            else:
                have[ro] = E83.pair_distances(o, ro)
    return have


def readings99(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = READOUTS) -> pd.DataFrame:
    """readings.csv: per variant and unit (each wiring; S and T where all four have a run) and readout E92's
    `score_row92` against every target of E89's `targets89` with E94's group effects, the runs outside the bands and at a
    bracket's bound; pairs.csv: each population's per-pair distances and effect under 'cal' and 'raw' (E96's
    `readings96` on this experiment's runs table)."""
    runs = runs_table(out) if runs is None else runs
    rows, prs = [], []
    for variant, g in runs.groupby('variant', sort=False):
        ds = {r.wiring: _distances(r.source, r.tag, r.wiring, readouts) for r in g.itertuples()}
        meta = {r.wiring: r for r in g.itertuples()}
        pv = parse_variant(variant)
        units = [(w, (w,)) for w in ds] + [(p, ws) for p, ws in (('S', S_WIRINGS), ('T', T_WIRINGS)) if all(w in ds for w in ws)]
        for unit, ws in units:
            for ro in readouts:
                d = ds[unit][ro] if len(ws) == 1 else E83.population([ds[w][ro] for w in ws])
                rows.append(dict(variant=variant, arm=pv['arm'], part=pv['part'], screen=pv['screen'], shape=pv['shape'],
                                 counts=pv['counts'], mode=pv['mode'], size=pv['size'], pool=pv['pool'], dose=pv['dose'],
                                 permute=pv['permute'], control=pv['control'], construction=pv['construction'],
                                 equal=pv['equal'], unit=unit, n_wirings=len(ws), readout=ro,
                                 n_flagged=int(sum(meta[w].flagged for w in ws)),
                                 n_at_bound=int(sum(meta[w].at_bound for w in ws)),
                                 **E92G.score_row92(d), **E94.group_effects(d)))
                if len(ws) > 1 and ro in ('cal', 'raw'):
                    me = E83._cs(H.apl_effect(d['on'], d['off']))
                    prs.append(pd.DataFrame(dict(variant=variant, unit=unit, readout=ro, odor_a=me.index.get_level_values(0),
                                                 odor_b=me.index.get_level_values(1),
                                                 on=E83._cs(d['on']).reindex(me.index).to_numpy(),
                                                 off=E83._cs(d['off']).reindex(me.index).to_numpy(), effect=me.to_numpy())))
    rd = pd.DataFrame(rows)
    rd.to_csv(out / 'readings.csv', index=False)
    if prs:
        pd.concat(prs, ignore_index=True).to_csv(out / 'pairs.csv', index=False)
    return rd


def _co_job99(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and the geometric mean of E88's two estimates at L1 and L3 (E96's
    `_co_job96` on this experiment's loader)."""
    import os
    import warnings
    source, tag, wiring, K, cand, rel, gap = job
    warnings.simplefilter('ignore')
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    o = load_any99(source, tag, wiring)
    mr = E89.ModelReading89(o, label=tag, K=K, rel_target=rel, gap_target=gap)
    res = dict(tag=tag, wiring=wiring, readings={})
    for rdn in ('r1', 'r1b'):
        F = mr.mean_factors(rdn)
        res['readings'][rdn] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    for level in ('L1', 'L3'):
        v = {c: mr.level_values(level, c, None)[0] for c in E89.CONDS}
        F = {c: E89.cand_factors(v[c], cand) for c in E89.CONDS}
        res['readings'][f'cand_{level}'] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    return res


def coresponse99(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse.csv: per S setting (all four S runs) the changes of L (co-response) and R (recruitment) under APL of the
    four nearest ester pairs, the ester and the alcohol groups on S's population, against the flies' (E96's
    `coresponse96` on this experiment's runs; jobs cached under outputs/.../coresponse_cache.pkl)."""
    import multiprocessing
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    from experiments import E88_responder_sets as E88
    runs = runs_table(out)
    runs = runs[runs.wiring.isin(S_WIRINGS)]
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    cache = local_dir() / 'coresponse_cache.pkl'
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    cand = E92G.CO_CANDIDATE
    full = {v for v, g in runs.groupby('variant') if set(g.wiring) == set(S_WIRINGS)}
    jobs = [(r.source, r.tag, r.wiring, E83.saved_fit(r.wiring)['K'], cand, rel, gap) for r in runs.itertuples()
            if (r.source, r.tag) not in done and r.variant in full]
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            for res, j in zip(ex.map(_co_job99, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E99 coresponse {j[1]}', flush=True)
    fl = {}
    fv = E89.flies_values({})
    for subset in E88.SUBSETS:
        tnt = E88.flies_subset(subset)
        rf = E88.flies_factors('z', 'r1', tnt)
        F = {c: rf[c].mean() for c in E89.CONDS}
        fl[('r1', subset)] = (E88.log_changes(F, ('L', 'R')), [E88.log_changes(b, ('L', 'R')) for b in E88.flies_boot(rf)])
        fl[('cand', subset)] = (E88.log_changes(E89.flies_cand(fv, cand, None, tnt), ('L', 'R')),
                                [E88.log_changes(b, ('L', 'R')) for b in E89.flies_cand_boot(fv, cand, None, tnt)])
    known = E92G.gm_known_error()
    rows = []
    for variant, g in runs.groupby('variant', sort=False):
        if variant not in full:
            continue
        for rdn in ('r1', 'r1b', 'cand_L1', 'cand_L3'):
            Fs = [done[(r.source, r.tag)]['readings'][rdn] for r in g.itertuples()]
            Fp = {c: E88.Factors(*(E88._mean_rows(np.stack([f[c][k] for f in Fs])) for k in ('s', 'O', 'R'))) for c in E89.CONDS}
            ch = E88.log_changes(Fp, ('L', 'R'))
            for subset in E88.SUBSETS:
                fch, fb = fl[('r1' if rdn.startswith('r1') else 'cand', subset)]
                for f in ('L', 'R'):
                    for grp in ('nearest4', 'ester', 'alcohol'):
                        k = f'chg_{f}_{grp}'
                        lo, hi = np.nanpercentile([b[k] for b in fb], E88.INTERVAL)
                        x = ch[k]
                        differs = bool(not (lo <= x <= hi) and abs(x - fch[k]) >= np.log(E88.B1_FACTOR))
                        rows.append(dict(variant=variant, reading=rdn, candidate=cand if rdn.startswith('cand') else '',
                                         subset=subset, factor=f, group=grp, ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])),
                                         flies_lo=float(np.exp(lo)), flies_hi=float(np.exp(hi)), differs=differs,
                                         known_error=known.get(rdn.split('_')[1], np.nan) if rdn.startswith('cand') else np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'coresponse.csv', index=False)
    return df


def fifths_of_points(out: Path, runs: Optional[pd.DataFrame] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3's reading for every point of this experiment on S (b_fifths.csv; b_claws.csv for W2, whose placement is
    uneven) with the saved corners at size 1 beside (A3_RUNS)."""
    runs = runs_table(out) if runs is None else runs
    mine = runs[(runs.source == 'E99') & runs.wiring.isin(S_WIRINGS)]
    full = [v for v, g in mine.groupby('variant', sort=False) if set(g.wiring) == set(S_WIRINGS)]
    specs = [(v, 'E99', v) for v in full] + [(v, src, shape) for v, (src, shape) in A3_RUNS.items()]
    f, c = who_responds(specs, S_WIRINGS, [v for v in full if v.startswith('W2_')])
    f.to_csv(out / 'b_fifths.csv', index=False)
    c.to_csv(out / 'b_claws.csv', index=False)
    return f, c


# ---------------------------------------------------------------------------
# part C: the verdicts

def mean_effect(rd: pd.DataFrame, variant: str, unit: str = 'S', readout: str = 'cal', bound_ok: bool = False) -> float:
    """A setting's mean APL effect on a unit under a readout; NaN when not read or (unless `bound_ok`) when a run of it
    on the unit is at a bracket's bound (such a point is neither a choice nor part of a share)."""
    q = rd[(rd.variant == variant) & (rd.unit == unit) & (rd.readout == readout)]
    if not len(q):
        return np.nan
    r = q.iloc[0]
    if not bound_ok and int(r.n_at_bound) > 0:
        return np.nan
    return float(r.mean_effect)


BETWEEN = {  # (name, U, O, K1, K2)
    'calyx_d1': (UNIFORM[1.0], OFFCLAW[1.0], vname('K1', 1, 'offclaw'), vname('K2', 1, 'offclaw')),
    'calyx_d2': (UNIFORM[2.0], OFFCLAW[2.0], vname('K1', 2, 'offclaw'), vname('K2', 2, 'offclaw')),
    'peduncle_d1': (UNIFORM[1.0], OFFCLAW_PED, vname('K1', 1, 'offclaw_ped'), vname('K2', 1, 'offclaw_ped')),
}
WITHIN = (UNIFORM[1.0], CONNECTOME, vname('W1', 1), vname('W2', 1))


def c_verdicts(rd: pd.DataFrame, a1: Optional[pd.DataFrame] = None, out: Optional[Path] = None) -> Tuple[dict, pd.DataFrame]:
    """Part C's verdicts on the mean APL effect under the calibrated readout (raw counts beside): the shares and classes
    between KCs at d = 1 on the calyx's counts on S and on each of S's wirings (with the rule by wiring), the two
    replicates (agreeing with the class or not), within a KC, T1 and T2 against B6 at each size, the real wiring against
    P1, P2 and M1 for W2, D2 and B6, how even is even enough and the KCs with none. Returns (c_verdicts, c_shares rows)."""
    rows = []
    res: dict = {}
    for key, (u, o, k1, k2) in BETWEEN.items():
        per = {}
        for unit in ('S',) + S_WIRINGS:
            for ro in ('cal', 'raw'):
                s = shares_between(*(mean_effect(rd, v, unit, ro) for v in (u, o, k1, k2)))
                rows.append(dict(kind='between', case=key, unit=unit, readout=ro, **s))
                if ro == 'cal':
                    per[unit] = s
        res[f'between_{key}'] = dict(S=per['S'], **by_wiring({w: per[w]['cls'] for w in S_WIRINGS}, per['S']['cls']))
    s_cls = res['between_calyx_d1']['S']['cls']
    for key in ('calyx_d2', 'peduncle_d1'):
        res[f'between_{key}']['agrees_with_class'] = bool(res[f'between_{key}']['S']['cls'] == s_cls)
    per = {}
    for unit in ('S',) + S_WIRINGS:
        for ro in ('cal', 'raw'):
            s = shares_within(*(mean_effect(rd, v, unit, ro) for v in WITHIN))
            rows.append(dict(kind='within', case='d1', unit=unit, readout=ro, **s))
            if ro == 'cal':
                per[unit] = s
    res['within_d1'] = dict(S=per['S'], **by_wiring({w: per[w]['cls'] for w in S_WIRINGS}, per['S']['cls']))
    # T1 and T2 against B6 at each size (what spreading the on-claw part evenly adds to all counted synapses)
    tb = []
    for t in T_SHAPES:
        for d in SIZES:
            v, b = vname(t, d), ALLSYN[d]
            r = dict(setting=v, b6=b, size=d)
            for ro in ('cal', 'raw'):
                r[f'mean_{ro}'] = mean_effect(rd, v, 'S', ro, bound_ok=True)
                r[f'b6_mean_{ro}'] = mean_effect(rd, b, 'S', ro, bound_ok=True)
                r[f'd_mean_{ro}'] = r[f'mean_{ro}'] - r[f'b6_mean_{ro}']
            r['d_pattern'] = E96._pat(rd, v, 'S', 'z8') - E96._pat(rd, b, 'S', 'z8')
            r['d_raw_loss'] = E96._pat(rd, v, 'S', 'z8', 'loss', 'raw') - E96._pat(rd, b, 'S', 'z8', 'loss', 'raw')
            tb.append(r)
    res['t_vs_b6'] = tb
    # the real wiring against P1, P2 and M1 for W2, D2 and B6 (a description of which claws)
    res['which_claws'] = {v: {w: dict(mean_cal=mean_effect(rd, v, w, 'cal', bound_ok=True),
                                      mean_raw=mean_effect(rd, v, w, 'raw', bound_ok=True),
                                      pattern=E96._pat(rd, v, w, 'z8')) for w in S_WIRINGS}
                          for v in (vname('W2', 1), CONNECTOME, ALLSYN[1.0])}
    # how even is even enough
    out = out_dir() if out is None else out
    a1 = a1 if a1 is not None else (pd.read_csv(out / 'a1_shapes.csv') if (out / 'a1_shapes.csv').exists() else None)
    cv_O = cv_K3 = np.nan
    if a1 is not None:
        real = a1[a1.wiring == 'real'].set_index('setting')
        cv_O, cv_K3 = float(real.loc['offclaw', 'within_cv_median']), float(real.loc['K3', 'within_cv_median'])
    U, O = mean_effect(rd, UNIFORM[1.0]), mean_effect(rd, OFFCLAW[1.0])
    K2, K3, K4 = (mean_effect(rd, vname(s, 1, 'offclaw')) for s in ('K2', 'K3', 'K4'))
    res['bound'] = evenness_bound(U, O, K2, K3, cv_O, cv_K3)
    a4p = out / 'a4_yardstick.csv'
    if a4p.exists():
        a4 = pd.read_csv(a4p)
        res['bound']['yardstick'] = a4[a4.kind.isin(['calyx_off_claws', 'calyx_on_claws', 'peduncle', 'off_claws_and_peduncle'])][
            ['source', 'dataset', 'apl', 'kind', 'within_cv_median']].to_dict('records')
    res['none'] = none_verdict(U, O, K4)
    res['flies_mean'] = FLIES_MEAN
    return res, pd.DataFrame(rows)


def compare_predictions(res: dict, a2: dict) -> dict:
    """The observed classes beside A2's predictions under each description."""
    out = {}
    if not a2:
        return out
    for name in ('i', 'ii'):
        c = a2['classes'][name]
        out[name] = dict(between=dict(predicted=c['between']['cls'], observed=res['between_calyx_d1']['S']['cls']),
                         peduncle=dict(predicted=c['replicate_peduncle']['cls'],
                                       observed=res['between_peduncle_d1']['S']['cls']),
                         within=dict(predicted=c['within']['cls'], observed=res['within_d1']['S']['cls']),
                         none=dict(predicted=c['none']['verdict'], observed=res['none']['verdict']))
    return out


# the choice, the gate, the equal rule (E96's for W1, T1 and T2)

PARTS_COLS = E92G.PARTS_COLS


def _settings99(rd: pd.DataFrame, screens: Sequence[str] = (REF,) + CHOICE_ARMS) -> pd.DataFrame:
    return rd[(rd.unit == 'S') & (rd.readout == 'cal') & rd.screen.isin(list(screens)) & ~rd.control.astype(bool)]


def choose99(rd: pd.DataFrame, half: str) -> dict:
    """The choice on S against half `half` (E96's `choose96`): the setting with the lowest pattern part of S's population
    under the calibrated readout among the unchanged fits and CHOICE_ARMS' settings whose four runs are inside the bands
    and not at a bracket's bound; something to confirm if it is a setting and gains CHOICE_MIN_GAIN or more."""
    g = _settings99(rd)
    ok = g[(g.variant == REF) | ((g.n_flagged == 0) & (g.n_at_bound == 0))]
    if ok.empty or REF not in set(g.variant):
        return dict(half=half, chosen=None, confirm=False, note='runs missing')
    best = ok.loc[ok[f'half{half}_pattern'].idxmin()]
    rr = g[g.variant == REF].iloc[0]
    gain = float(rr[f'half{half}_pattern']) - float(best[f'half{half}_pattern'])
    confirm = bool(best.screen in CHOICE_ARMS and gain >= CHOICE_MIN_GAIN)
    row = dict(half=half, chosen=best.variant, screen=best.screen, pattern=float(best[f'half{half}_pattern']),
               ref_pattern=float(rr[f'half{half}_pattern']), gain=gain, confirm=confirm, n_eligible=int(len(ok)),
               n_settings=int(len(g)), excluded=json.dumps(sorted(set(g.variant) - set(ok.variant))),
               control=control_of99(best.variant) if best.variant != REF else '')
    for c in PARTS_COLS[1:]:
        row[f'chosen_{c}'] = float(best[f'half{half}_{c}'])
        row[f'ref_{c}'] = float(rr[f'half{half}_{c}'])
    return row


def plan_c(rd: pd.DataFrame) -> Dict[str, dict]:
    return {h: choose99(rd, h) for h in HALVES}


def equal_rule99(rd: pd.DataFrame) -> List[dict]:
    """E94's equal rule at fixed scales for W1, T1 and T2 (E96's `equal_rule96`)."""
    g = _settings99(rd)
    out = []
    if REF not in set(g.variant):
        return out
    me0 = float(g[g.variant == REF].mean_effect.iloc[0])
    for r in g.itertuples():
        if r.screen not in CHOICE_ARMS:
            continue
        d = {t: E96._pat(rd, r.variant, 'S', t) - E96._pat(rd, REF, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = E96._pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - E96._pat(rd, REF, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - me0
        rest = bool(int(r.n_flagged) == 0 and int(r.n_at_bound) == 0 and abs(raw) <= EQUAL_RAW and abs(me) <= EQUAL_MEAN)
        out.append(dict(variant=r.variant, screen=r.screen, size=float(r.size), **{f'd_pattern_{t}': x for t, x in d.items()},
                        d_raw_loss=raw, d_mean_effect=me, n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound),
                        equal=bool(all(abs(x) <= EQUAL_PATTERN for x in d.values()) and rest),
                        better=bool(all(x < 0 for x in d.values()) and rest)))
    return out


def vs_uniform(rd: pd.DataFrame) -> List[dict]:
    """Each setting of this experiment against E95's uniform component at the same size (on S under the readout against
    all recordings and halves A and B, the mean APL effect)."""
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    have = set(g.variant)
    out = []
    for r in g[g.part == 'B'].itertuples():
        u = f'D1_uniform_d{float(r.size):g}'
        row = dict(variant=r.variant, uniform=u, uniform_read=u in have)
        if u in have:
            for t in ('z8', 'halfA', 'halfB'):
                row[f'd_pattern_{t}'] = E96._pat(rd, r.variant, 'S', t) - E96._pat(rd, u, 'S', t)
            row['d_mean_effect'] = float(r.mean_effect) - float(g[g.variant == u].mean_effect.iloc[0])
        out.append(row)
    return out


def gate_rows(rd: pd.DataFrame, runs: pd.DataFrame, plan: Dict[str, dict]) -> List[dict]:
    """The gate's four directions (E92's `gate_direction92` with the rule (v) control, E95's uniform component at the same
    size)."""
    dirs = []
    for X in HALVES:
        c = plan[X]
        chosen = c['chosen'] if c['confirm'] else None
        r = E92G.gate_direction92(rd, runs, X, chosen, REF, [control_of99(chosen)] if chosen else [])
        if chosen:
            r['rule_v_control'] = control_of99(chosen)
        dirs.append(r)
    return dirs


def part_c(out: Path) -> dict:
    """C: the choice per half (choice.csv); each chosen setting with something to confirm run on T (runs_C.csv) with its
    rule (v) control on T (E95's uniform component at the same size, through this experiment's path); then `score`
    without the co-response."""
    runs = runs_table(out)
    rd = readings99(out, runs, readouts=('cal', 'raw'))
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    chosen = sorted({c['chosen'] for c in plan.values() if c['confirm']})
    if not chosen:
        print('E99 C: nothing to confirm', flush=True)
    for v in chosen:
        points(out, [v], T_WIRINGS, RUNS_NAMES['C'], 'C')
        points(out, [control_of99(v)], T_WIRINGS, RUNS_NAMES['C'], 'C')
    return score(out, coresponse=False)


# ---------------------------------------------------------------------------
# part D: the rule, the starts, the commands, the scores

def d_rule99(rd: pd.DataFrame) -> dict:
    """Part D's rule (the entry's): a setting of T1 or T2 on S under the readout against all recordings meets it with a
    mean APL effect of D_MEAN_MIN or more, a pattern part within D_VS_UNIFORM of E95's uniform component's at the same
    size or lower, and its four runs inside the bands and not at a bound. The component fitted: T2 if it meets the rule
    and its best pattern part (among its settings that meet it) is within D_T2_TIE of T1's (or T1 has none that meets
    it); otherwise the one that meets it with the lower pattern part."""
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    rows, best = [], {}
    for r in g[g.screen.isin(list(T_SHAPES))].itertuples():
        u = f'D1_uniform_d{float(r.size):g}'
        up = E96._pat(rd, u, 'S', 'z8')
        meets = bool(float(r.mean_effect) >= D_MEAN_MIN and np.isfinite(up) and float(r.z8_pattern) <= up + D_VS_UNIFORM
                     and int(r.n_flagged) == 0 and int(r.n_at_bound) == 0)
        rows.append(dict(variant=r.variant, shape=r.screen, size=float(r.size), mean_effect=float(r.mean_effect),
                         pattern=float(r.z8_pattern), uniform=u, uniform_pattern=up, n_flagged=int(r.n_flagged),
                         n_at_bound=int(r.n_at_bound), meets=meets))
        if meets and (r.screen not in best or float(r.z8_pattern) < best[r.screen]['pattern']):
            best[r.screen] = rows[-1]
    out = dict(rows=rows, best={k: v['variant'] for k, v in best.items()}, mean_min=D_MEAN_MIN, vs_uniform=D_VS_UNIFORM,
               runs=bool(best))
    if 'T2' in best and ('T1' not in best or best['T2']['pattern'] <= best['T1']['pattern'] + D_T2_TIE):
        out['component'] = 'T2'
    elif best:
        out['component'] = min(best, key=lambda k: best[k]['pattern'])
    else:
        out['component'] = None
    out['mode'] = T_MODE.get(out['component']) if out['component'] else None
    return out


def d_tag(label: str, comp: str, start: str) -> str:
    """The tag E10_joint_fit.py gives a part D fit (its prefix rules: the wiring, rate0.068_, aplboth_, aplpnhold_,
    flatt1_ / flatt2_, then k2 and D_TAGS[start])."""
    from experiments import E10_joint_fit as J
    a = E91.d_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return (f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.claw_inh_pool_tag(None)}{J.aplpn_hold_tag(1.0)}"
            f"{J.aplkc_flat_tag(dict(kind=T_MODE[comp]))}k2{D_TAGS[start]}")


def d_start_key(label: str, comp: str, start: str) -> str:
    return f'{label}_{comp}_{start}'


def own_start(rd: pd.DataFrame, runs: pd.DataFrame, label: str, comp: str) -> dict:
    """The 'own' start of a wiring: the absolute parameters of its best point of `comp` (T1 / T2) in part B (the lowest
    pattern part of that wiring under the readout against all recordings, inside the bands first, none at a bracket's
    bound), with that point's size and unit."""
    q = runs[(runs.screen == comp) & (runs.wiring == label) & ~runs.at_bound]
    if q.empty:
        raise RuntimeError(f'{label}: no {comp} point inside the bracket (run --part B first)')
    sc = rd[(rd.unit == label) & (rd.readout == 'cal')].drop_duplicates('variant').set_index('variant').z8_pattern
    q = q.assign(pattern=q.variant.map(sc), flag_order=q.flagged.astype(int)).sort_values(['flag_order', 'pattern'])
    r = q.iloc[0]
    return dict({k: float(r[k]) for k in E83.PARAM_KEYS}, aplkc_flat=float(r['size']), from_tag=r.tag,
                pattern=float(r.pattern), flagged=bool(r.flagged), component_unit=float(r['component_unit']), start='own')


def e96_kept_refit(label: str) -> dict:
    """E96's kept refit of a wiring (its d_scores.csv row picked for mode offclaw without the pooled share; the fit with
    the lower final value of its own loss if 'picked' is missing) with its fits.csv row."""
    out96 = H.results_dir(E96.EXPERIMENT)
    s = pd.read_csv(out96 / 'd_scores.csv')
    q = s[(s['mode'] == 'offclaw') & (s['pool'].fillna(0) == 0) & (s.unit == label) & s.tag.notna()]
    if q.empty:
        raise RuntimeError(f"{label}: E96 has no refit in d_scores.csv")
    if 'picked' in q and q.picked.astype(str).isin(['True', 'true', '1']).any():
        r = q[q.picked.astype(str).isin(['True', 'true', '1'])].iloc[0]
    else:
        r = q.loc[q.fit_loss.astype(float).idxmin()]
    f = pd.read_csv(out96 / 'fits.csv').drop_duplicates('tag', keep='last').set_index('tag')
    if r.tag not in f.index:
        raise RuntimeError(f"{label}: E96's refit {r.tag} is not in its fits.csv")
    return dict(tag=r.tag, start=r.start, row=f.loc[r.tag])


def e96_start(label: str, unit: float) -> dict:
    """The 'e96' start of a wiring: E96's kept refit of that wiring (`e96_kept_refit`), its absolute parameters and fitted
    size, the component exchanged; the unit is this experiment's for the wiring (`unit`, checked against E96's)."""
    k = e96_kept_refit(label)
    r = k['row']
    u96 = float(r['aplkc_flat_unit'])
    if abs(u96 - unit) > 1e-9 * unit:
        raise RuntimeError(f"{label}: E96's unit {u96} is not this experiment's {unit}")
    p = {c: float(r[c]) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')}
    p['wKCKC_scale'] = float(r['kc_kc_inh'])
    return dict(p, aplkc_flat=float(r['aplkc_flat_size']), from_tag=f"{E96.EXPERIMENT}:{k['tag']}", component_unit=float(unit),
                unit_e96=u96, e96_start=k['start'], start='e96')


def t_start(label: str, size: float) -> dict:
    """T's start: the wiring's own fit with the APL->bouton scale at dose D_DOSE and the component at `size` (the median
    of S's fitted sizes); the unit the wiring's."""
    return E96.t_start(label, size)


def write_start(out: Path, key: str, start: dict) -> None:
    path = out / D_START_FILE
    d = _read_json(path)
    d[key] = start
    _write_json(path, d)


def d_command(label: str, comp: str, start: str, unit: float, maxfev: int = D_MAXFEV) -> List[str]:
    """One part D fit's E10_joint_fit.py command (E96's `d_command` with --aplkc-flat allsyn_even / allsyn_even_calyx, no
    pooled share, into this experiment's folders, from d_start_params.json's `d_start_key`)."""
    m = E91.d_member(label)
    F = E94.aplpn_multiple(label, D_DOSE)
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
           repr(float(F)), '--aplkc-flat', T_MODE[comp], '--aplkc-flat-unit', repr(float(unit)), '--aplkc-flat-max',
           f'{D_FLAT_BOUNDS[1]:g}', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag',
           D_TAGS[start], '--maxfev', str(int(maxfev)), '--experiment', EXPERIMENT,
           '--init-params', f'results/{EXPERIMENT}/{D_START_FILE}:{d_start_key(label, comp, start)}']
    return cmd + m['wiring_args']


def _fits(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def score_key(unit: str, comp: str) -> str:
    return f'{unit}_{comp}'


def d_allowed(out: Path, pop: str, comp: Optional[str] = None) -> Tuple[bool, dict]:
    """S's fits: by part D's rule (`d_rule99`); T's: only if S's fitted population is an equal or better (`d_scores`)."""
    if not (out / 'readings.csv').exists():
        return False, dict(note='not scored (run --part score)')
    rule = d_rule99(pd.read_csv(out / 'readings.csv'))
    comp = comp or rule.get('component')
    if pop == 'S':
        return bool(rule['runs'] and comp is not None), dict(rule, component=comp)
    s = d_scores(out, write=False).get(score_key('S', comp)) if comp else None
    if not s:
        return False, dict(rule, note="S's fits not scored")
    return bool(s.get('equal') or s.get('better')), dict(rule, component=comp, s_scores=s)


def part_d(out: Path, pop: str = 'S', start: str = 'own', wirings: Optional[Sequence[str]] = None,
           comp: Optional[str] = None, force: bool = False, dry_run: bool = False) -> List[str]:
    """Part D's fits of population `pop` from `start` (or `wirings` of it: a lane), one process after another
    (`d_command`; logs under outputs/.../D/), skipping those in fits.csv; refuses what the entry's rule does not allow
    unless `force`. `comp` None: the rule's component."""
    import subprocess
    ok, plan = d_allowed(out, pop, comp)
    comp = comp or plan.get('component')
    print(f"E99 D {pop} {start}: {'allowed' if ok else 'not allowed'}; component {comp}; "
          f"{json.dumps({k: plan.get(k) for k in ('best', 'runs', 'mean_min', 'vs_uniform', 'note')}, default=str)}",
          flush=True)
    for r in plan.get('rows', []):
        print(f"  {r['variant']}: mean {r['mean_effect']:+.3f}, pattern {r['pattern']:.4f} (uniform {r['uniform_pattern']:.4f}), "
              f"flagged {r['n_flagged']}, at bound {r['n_at_bound']}: {'meets' if r['meets'] else 'does not meet'}", flush=True)
    if (not ok and not force) or comp is None:
        return []
    if pop == 'T' and start != 'own':
        raise ValueError("T's fits start from each wiring's own fit ('own')")
    ws = [w for w in _wirings(pop) if wirings is None or w in wirings]
    if not dry_run:
        _commits()
    rd = pd.read_csv(out / 'readings.csv')
    runs = runs_table(out)
    logs = local_dir('D')
    tags = []
    sizes = None
    if pop == 'T':
        sizes = (d_scores(out, write=False).get(score_key('S', comp)) or {}).get('sizes')
        if not sizes:
            raise RuntimeError("T's fits need S's fitted sizes (run --part Dscore)")
    for w in ws:
        tag = d_tag(w, comp, start)
        tags.append(tag)
        if tag in set(_fits(out).tag):
            print(f'E99 D {tag}: done', flush=True)
            continue
        if pop == 'T':
            st = t_start(w, float(np.median(sizes)))
        elif start == 'own':
            st = own_start(rd, runs, w, comp)
        else:
            unit = float(runs[(runs.wiring == w) & (runs.source == 'E99')].component_unit.dropna().iloc[0])
            st = e96_start(w, unit)
        key = d_start_key(w, comp, start)
        cmd = d_command(w, comp, start, st['component_unit'])
        print(f'E99 D {key}: start {json.dumps(st, default=str)}', flush=True)
        print('E99 D:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        write_start(out, key, st)
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E99 D {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits(out).tag):
            raise RuntimeError(f'E99 D: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def load_fit99(tag: str) -> Optional[Dict[str, tuple]]:
    """A part D fit's outputs: this host's outputs, the archive, the archive's from_hal/."""
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    for d in (H.outputs_dir(EXPERIMENT), arc, arc / 'from_hal'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def e96_refit_s() -> dict:
    """E96's refit with the calyx's off-claw synapses, S's population (d_scores.json 'S_offclaw'; the entry's 0.0541,
    +0.165, 0.0838)."""
    s = _read_json(H.results_dir(E96.EXPERIMENT) / 'd_scores.json').get('S_offclaw', {})
    return dict(pattern_all=float(s.get('pattern_all', 0.0541)), mean_effect=float(s.get('mean_effect', 0.165)),
                raw_loss=float(s.get('raw_loss', 0.0838)))


def d_scores(out: Path, write: bool = True) -> dict:
    """Part D's scores (E96's `d_scores` for T1 and T2): per wiring and start the readings under the readout and on raw
    counts against the unchanged fits of the same wirings, with realism; per wiring the fit with the lower final value of
    its own loss (`E96.pick_start`); S by E94's D1 rule (an equal / better as a fitted family), beside E95's refit (as
    good as the uniform component: within D_AS_GOOD of it, or lower), E96's refit and E94's; the eight against
    UNCHANGED_8 by the same rule and E91's part D pass rule. d_scores.csv, d_scores.json."""
    fits = _fits(out)
    if fits.empty:
        return {}
    fits = fits.drop_duplicates('tag', keep='last')
    rows, res = [], {}
    u95 = E96.e95_refit_s()
    e96r = e96_refit_s()
    for comp in T_SHAPES:
        st = {}
        for w in ALL_WIRINGS:
            cand = []
            for start in D_TAGS:
                try:
                    tag = d_tag(w, comp, start)
                except Exception:
                    continue
                if tag not in set(fits.tag):
                    continue
                o = load_fit99(tag)
                if o is None:
                    continue
                d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
                d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
                sc, sc0 = E92G.score_row92(d['cal']), E92G.score_row92(d0['cal'])
                raw, raw0 = E83.effect_scores(d['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
                fr = fits[fits.tag == tag].iloc[-1]
                row = dict(component=comp, mode=T_MODE[comp], unit=w, start=start, tag=tag, **E94._d_row(sc, sc0, raw, raw0),
                           lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                           **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale',
                                                     'wPNAPL_scale', 'kc_kc_inh', 'aplpn_hold', 'aplpn_factor_vs_e3',
                                                     'aplkc_flat_size', 'aplkc_flat_unit', 'init_params', 'n_evals',
                                                     'initial_loss')},
                           fit_loss=float(fr.get('loss')), **E83.realism(o))
                cand.append((row, d, d0))
            if not cand:
                continue
            pick = E96.pick_start([c[0] for c in cand])
            for c in cand:
                c[0]['picked'] = c[0] is pick
                rows.append(c[0])
            st[w] = next(c for c in cand if c[0] is pick)
        for unit, ws in (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS)):
            if not all(w in st for w in ws):
                continue
            dp = {ro: E83.population([st[w][1][ro] for w in ws]) for ro in ('cal', 'raw')}
            d0 = {ro: E83.population([st[w][2][ro] for w in ws]) for ro in ('cal', 'raw')}
            sc, sc0 = E92G.score_row92(dp['cal']), E92G.score_row92(d0['cal'])
            raw, raw0 = E83.effect_scores(dp['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
            r = dict(E94._d_row(sc, sc0, raw, raw0), n_lower=int(sum(st[w][0]['lower'] for w in ws)),
                     n_flagged=int(sum(E94._flagged(st[w][0].get('realism_flag')) for w in ws)),
                     sizes=[float(st[w][0].get('aplkc_flat_size', np.nan)) for w in ws],
                     starts=[st[w][0]['start'] for w in ws])
            if unit == 'S':
                r['equal'], r['better'] = E94.d1_rule(r)
                r['vs_e95_refit'] = float(r['pattern_all'] - u95)
                r['as_good_as_uniform'] = bool(r['pattern_all'] <= u95 + D_AS_GOOD)
                r['vs_e96_refit'] = float(r['pattern_all'] - e96r['pattern_all'])
                r['vs_e94_refit'] = float(r['pattern_all'] - E94_REFIT_S)
            if unit == 'all8':
                rest = bool(r['n_flagged'] == 0 and abs(r['raw_loss'] - r['raw_loss_unchanged']) <= E94.D_RAW_SLACK
                            and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
                r['equal'] = bool(abs(r['pattern_all'] - UNCHANGED_8) <= E94.D_EQUAL and rest)
                r['passes'] = bool(r['pattern_all'] <= UNCHANGED_8 - E94.D_PASS_GAIN and r['n_lower'] >= E94.D_WIRINGS_LOWER
                                   and r['n_flagged'] == 0 and r['raw_loss'] <= E94.D_RAW_REF + E94.D_RAW_SLACK
                                   and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
            res[score_key(unit, comp)] = r
            rows.append(dict(component=comp, mode=T_MODE[comp], unit=unit,
                             **{k: x for k, x in r.items() if k not in ('sizes', 'starts')}))
    res['beside'] = dict(e95_refit_S=u95, e96_refit_S=e96r, e94_refit_S=E94_REFIT_S)
    if write:
        pd.DataFrame(rows).to_csv(out / 'd_scores.csv', index=False)
        _write_json(out / 'd_scores.json', res)
    return res


def d_check(out: Path, modes: Sequence[str] = tuple(T_MODE.values())) -> dict:
    """E10_joint_fit.py with --aplkc-flat allsyn_even / allsyn_even_calyx at size 0 against E94's D1 fit's first
    evaluation on the real wiring (E96's `d_check` on these modes, into this experiment's results folder, tag suffix
    _e10check); each first evaluation's loss terms must equal E94's history row 0 exactly. d_check.json."""
    import subprocess
    from experiments import E10_joint_fit as J
    hist94 = pd.read_csv(H.results_dir(E94.EXPERIMENT) / f"history_{E94.d1_tag('real')}.csv")
    m = E91.d_member('real')
    base = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
            repr(float(E94.aplpn_multiple('real', float(_read_json(H.results_dir(E94.EXPERIMENT) / E94.D_START_FILE)[
                E94.d_start_key('real')]['dose'])))), '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both',
            '--maxfev', '5', '--smoke', '--experiment', EXPERIMENT,
            '--init-params', f'results/{E94.EXPERIMENT}/{E94.D_START_FILE}:{E94.d_start_key("real")}'] + m['wiring_args']
    L = E95.kc_losses('real')
    wn = E96.wiring_wnorm('real')
    unit = E95.component_unit(L['bclaw_mean'], E83.saved_fit('real')['p']['wAPLKC_scale'], L['apl_mean'], float(wn.mean()))
    res = dict(reference=f"{E94.EXPERIMENT}: history_{E94.d1_tag('real')}.csv row 0", unit=unit)
    cols = [c for c in hist94.columns if c not in ('time_s', 'eval', 'Unnamed: 0')]
    for mode in modes:
        cmd = base + ['--tag', '_e10check', '--aplkc-flat', mode, '--aplkc-flat-unit', repr(unit), '--aplkc-flat-init', '0']
        log = local_dir('D') / f'd_check_{mode}.log'
        with open(log, 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E99 d_check {mode}: E10_joint_fit.py exited {r.returncode} (log {log})')
        tag = E94.d1_tag('real').replace(E94.D1_TAG, '_e10check').replace(
            'aplpnhold_', 'aplpnhold_' + J.aplkc_flat_tag(dict(kind=mode)))
        h = pd.read_csv(out / f'history_{tag}.csv')
        diff = {c: [float(h[c].iloc[0]), float(hist94[c].iloc[0])] for c in cols if c in h.columns
                and not E95.E10_same(h[c].iloc[0], hist94[c].iloc[0])}
        res[mode] = dict(tag=tag, n_cols=int(sum(c in h.columns for c in cols)), differs=diff, exact=not diff,
                         flat_first=float(h['aplkc_flat'].iloc[0]) if 'aplkc_flat' in h else None,
                         flat_later=[float(x) for x in h['aplkc_flat'].iloc[1:]] if 'aplkc_flat' in h else None)
        print(f'E99 d_check {mode}: {json.dumps(res[mode], default=str)[:600]}', flush=True)
    res['passes'] = bool(all(res[m_]['exact'] and res[m_]['flat_first'] == 0.0 and res[m_]['n_cols'] > 0 for m_ in modes))
    _write_json(out / 'd_check.json', res)
    return res


def part_e(out: Path) -> None:
    """Part E (only on an equal as a fitted family): T's four refits run as `--part D --pop T`; the eight are read
    against 0.0482 by E91's part D rule in `d_scores`; the stage readings of the fits as in E96's part E are not written
    (they run only by that rule)."""
    ds = _read_json(out / 'd_scores.json')
    due = any(bool(x.get('equal') or x.get('better')) for k, x in ds.items() if k.startswith('S_') and isinstance(x, dict))
    print(f'E99 E: due {due}', flush=True)
    raise NotImplementedError("E99 part E's stage readings run only on an equal as a fitted family; not written")


# ---------------------------------------------------------------------------
# the verdict and score

def verdict99(out: Path, rd: Optional[pd.DataFrame] = None, runs: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: part C's verdicts (c_verdicts.json, c_shares.csv) with A2's predictions beside, the choices, the
    gate in four directions, the equals, each setting against E95's uniform component, S's settings with their parts,
    part D's rule (d_rule.json) and scores, the flagged runs and those at a bracket's bound, the alcohol pairs'
    co-response change."""
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table(out) if runs is None else runs
    cv, cs = c_verdicts(rd, out=out)
    cv['vs_predictions'] = compare_predictions(cv, _read_json(out / 'a2_predictions.json'))
    _write_json(out / 'c_verdicts.json', cv)
    cs.to_csv(out / 'c_shares.csv', index=False)
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    dirs = gate_rows(rd, runs, plan)
    pd.DataFrame(dirs).to_csv(out / 'gate.csv', index=False)
    gate = E92G.gate92(dirs)
    eqs = equal_rule99(rd)
    pd.DataFrame(eqs).to_csv(out / 'equal.csv', index=False)
    vu = vs_uniform(rd)
    pd.DataFrame(vu).to_csv(out / 'vs_uniform.csv', index=False)
    rule = d_rule99(rd)
    _write_json(out / 'd_rule.json', rule)
    gains_s = any(c['confirm'] for c in plan.values())
    v = dict(thresholds=dict(CLASS_HI=CLASS_HI, CLASS_LO=CLASS_LO, BY_WIRING_MIN=BY_WIRING_MIN, BOUND_LOSS=BOUND_LOSS,
                             BOUND_MIN_SHARE=BOUND_MIN_SHARE, NONE_HI=NONE_HI, NONE_LO=NONE_LO,
                             CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=E91.GATE_GAIN, EQUAL_PATTERN=EQUAL_PATTERN,
                             EQUAL_RAW=EQUAL_RAW, EQUAL_MEAN=EQUAL_MEAN, APL_BRACKET=APL_BRACKET, D_MEAN_MIN=D_MEAN_MIN,
                             D_VS_UNIFORM=D_VS_UNIFORM, D_T2_TIE=D_T2_TIE, D_AS_GOOD=D_AS_GOOD),
             readout="cal (the calibrated readout, K set on each run's APL-on counts): the verdicts, the choice and the gate",
             reference=REF, classes=dict(between=cv['between_calyx_d1']['S']['cls'],
                                         between_by_wiring=cv['between_calyx_d1']['holds_by_wiring'],
                                         replicate_d2=cv['between_calyx_d2']['S']['cls'],
                                         replicate_peduncle=cv['between_peduncle_d1']['S']['cls'],
                                         within=cv['within_d1']['S']['cls'], within_by_wiring=cv['within_d1']['holds_by_wiring'],
                                         none=cv['none']['verdict']),
             choice=plan, gate=gate,
             verdict='passes' if gate['passes'] else ('gains on S, not through the gate' if gains_s else 'level'),
             equals=eqs, equal_settings=[e['variant'] for e in eqs if e['equal']], vs_uniform=vu,
             D_rule=dict(runs=rule['runs'], component=rule['component'], best=rule['best']))
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    v['S'] = {r.variant: dict(screen=r.screen, halfA=r.halfA_pattern, halfB=r.halfB_pattern, halfC=r.halfC_pattern,
                              halfD=r.halfD_pattern, all=r.z8_pattern, six=r.z6_pattern, loss=r.z8_loss, mean_effect=r.mean_effect,
                              amplitude=r.z8_amplitude, shape=r.z8_shape, sd_ratio=r.z8_sd_ratio, pearson=r.z8_pearson,
                              rho=r.z8_rho, n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound)) for r in g.itertuples()}
    if (out / 'd_scores.json').exists():
        v['D'] = _read_json(out / 'd_scores.json')
    mine = runs[runs.source == 'E99']
    v['flagged_runs'] = mine[mine.flagged][['tag', 'realism_flag']].values.tolist()
    v['at_bound_runs'] = mine[mine.at_bound].tag.tolist()
    if (out / 'coresponse.csv').exists():
        try:
            co = pd.read_csv(out / 'coresponse.csv')
            reach = co[(co.factor == 'L') & (co.group == 'alcohol')]
            v['alcohol_coresponse'] = {f'{r.variant}:{r.reading}:{r.subset}': dict(ratio=r.ratio, differs=bool(r.differs))
                                       for r in reach.itertuples()}
        except pd.errors.EmptyDataError:
            pass
    _write_json(out / 'verdict.json', v)
    print(f"E99 verdict: between {v['classes']['between']} (by wiring {v['classes']['between_by_wiring']}), within "
          f"{v['classes']['within']}, none: {v['classes']['none']}; choice {v['verdict']}; equals {v['equal_settings']}; "
          f"D {'runs' if rule['runs'] else 'does not run'} ({rule['component']})", flush=True)
    return v


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """readings.csv, pairs.csv, line.csv, coresponse.csv (with `coresponse`), b_fifths.csv, b_claws.csv, c_shares.csv,
    c_verdicts.json, choice.csv, gate.csv, equal.csv, vs_uniform.csv, d_rule.json, verdict.json, the figure and
    comparison_table.md."""
    runs = runs_table(out)
    rd = readings99(out, runs, readouts=READOUTS if zcell else tuple(r for r in READOUTS if not r.startswith('zcell')))
    try:
        E95.line95(out)
    except Exception as e:                       # (E92 part A's line needs its line.json; reported, not fatal)
        print(f'E99 score: line.csv left out ({e!r})', flush=True)
    if coresponse:
        coresponse99(out, n_workers)
    try:
        fifths_of_points(out, runs)
    except Exception as e:                       # (a run missing its saved outputs on this host; reported, not fatal)
        print(f'E99 score: b_fifths.csv left out ({e!r})', flush=True)
    v = verdict99(out, rd, runs)
    figure(out)
    comparison_table(out)
    return v


# ---------------------------------------------------------------------------
# the steps from the claws to the somata, for a chosen or equal setting

def steps_variants(out: Path) -> List[str]:
    v = _read_json(out / 'verdict.json')
    vs = [c['chosen'] for c in v.get('choice', {}).values() if c.get('confirm')]
    vs += list(v.get('equal_settings', []))
    return sorted(dict.fromkeys(vs))


def part_steps(out: Path, variants: Optional[Sequence[str]] = None) -> Optional[pd.DataFrame]:
    """The step from the claws to the somata in wild-type flies for `variants` (default every chosen or equal setting)
    beside C2's family on S (E96's `part_steps` on this experiment's runs): somata_steps.csv, claws_steps.csv, steps.csv."""
    from experiments import E90_model_stages as E90
    from experiments import E92_steps_like_for_like as A
    variants = list(variants) if variants else steps_variants(out)
    if not variants:
        print('E99 steps: no chosen or equal setting', flush=True)
        return None
    t0 = time.time()
    runs = runs_table(out)
    soma, claw = [], []
    loaders = {'E99': load_stage99, 'E96': E96.load_stage96, 'E95': E95.load_stage95, 'E94': E94.load_stage94}
    i = 0
    for v in variants:
        for w in S_WIRINGS:
            q = runs[(runs.variant == v) & (runs.wiring == w)]
            if not len(q):
                raise RuntimeError(f'{v} has no run on {w}')
            r = q.iloc[0]
            o = load_any99(r.source, r.tag, w)
            soma.append(A.somata_pairs(o['on'][1], o['off'][1]).assign(family=v, wiring=w, source=r.source, tag=r.tag))
            st_on = loaders[r.source](r.tag) if r.source in loaders else dict(E94.fit_stage(w, 'on'))
            if st_on is None:
                raise RuntimeError(f'{r.tag}: no stage matrices')
            if 'spike_counts' not in st_on:
                st_on['responses'], st_on['spike_counts'] = o['on']
            W, _ = E94.wiring_W(w)
            meta = {k: np.asarray(x) for k, x in E92G.units_meta(W).items()}
            rr = E90.read_wiring(w, st_on, dict(E94.fit_stage(w, 'off')), meta,
                                 np.random.default_rng(E90.NOISE_SEED + E94.STEPS_NOISE_SEED + 300 + i))
            claw.append(pd.DataFrame(rr).assign(family=v, tag=r.tag, wiring=w))
            i += 1
    for w in S_WIRINGS:
        o = E83.saved_fit(w)['saved']
        soma.append(A.somata_pairs(o['on'][1], o['off'][1]).assign(family='C2', wiring=w, source='saved', tag=w))
    sm = E92G.f_soma_means(pd.concat(soma, ignore_index=True))
    sm.to_csv(out / 'somata_steps.csv', index=False)
    cl = E92G.claw_groups(pd.concat(claw, ignore_index=True))
    c2 = A.model_claws()
    c2 = c2[c2.wiring.isin(S_WIRINGS)].assign(family='C2')
    claws = pd.concat([c2, cl], ignore_index=True)
    claws.to_csv(out / 'claws_steps.csv', index=False)
    e90v = _read_json(H.results_dir('E90_model_stages') / 'verdict.json')
    rows = E92G.f_step_rows(sm, claws, A.check_passes(), bool(e90v.get('noise_check_passed', {}).get('claw', False)))
    rows.to_csv(out / 'steps.csv', index=False)
    print(f'E99 steps in {time.time() - t0:.0f} s', flush=True)
    return rows


# ---------------------------------------------------------------------------
# the archive, the figure and the table

SHAPE_NAMES = {
    'K1': "the off-claw counts in proportion: each claw of KC k gets (a_k / n_k) / mu_G (every claw-count group the same "
          'amount per claw, each KC its count relative to its group)',
    'K2': "the group's mean: each claw of a KC of claw-count group G gets mu_G (the counts' slope, no inequality within a "
          'group)',
    'K3': 'half the inequality: mu_G + (a_k / n_k - mu_G) / 2 on each claw',
    'K4': 'presence: the same on every claw of a KC with an off-claw synapse, nothing on the others',
    'W1': "each KC's on-claw count c_k spread evenly over its claws",
    'W2': 'equal amounts per KC (n_k), placed as the on-claw counts lie (c_kj / (c_k / n_k))',
    'T1': "each KC's count of all of APL's synapses (on its kept claws, on no claw in the calyx, in the peduncle) spread "
          'evenly over its claws',
    'T2': "each KC's count of APL's synapses on its kept claws and on no claw in the calyx spread evenly over its claws",
}
COUNTS_NAMES = {'offclaw': "the calyx's off-claw counts (a^c_k)", 'offclaw_ped': "the off-claw counts with the peduncle's "
                                                                                  '(a^c_k + a^p_k)'}


def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike counts,
    the harness's primary model RDM): one row per run of this experiment with its absolute parameters, the component,
    its size and unit, the dose, realism; part D's fits' distances, which E10_joint_fit.py wrote, are kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E99']
    rows, parts = [], []
    for _, r in runs.iterrows():
        pv = parse_variant(r.variant)
        if pv['part'] == 'B':
            desc = SHAPE_NAMES[pv['shape']] + (f", on {COUNTS_NAMES[pv['counts']]}" if pv['counts'] else '')
        else:
            desc = E96.MODE_NAMES.get(pv['mode'], f"E95's {pv['mode']} component")
        kind = ("C2's family, the tied APL scale re-held; the APL->KC weights the family's"
                + f" plus a component of size {pv['size']:g} (unit {r.get('component_unit'):.6g}, mean 1 over the claws "
                  f"times the connectome's mean weight): {desc}"
                + f"; APL->bouton scale x {pv['dose']:g} of the fit's; the APL-off run the fit's")
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=r.arm, kind=kind, config=r.variant, shape=pv['shape'],
                         counts=pv['counts'], construction=pv['construction'], size=r['size'], dose=r.dose, wiring=r.wiring,
                         ref_experiment=r.fit_folder, ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, factor=r.factor, aplpn_multiple=r.aplpn_multiple,
                         component_unit=r.get('component_unit'), synapse_share=r.get('synapse_share'), rate_on=r.rate_on,
                         rate_off=r.rate_off, max_count_on=r.max_count_on, max_count_off=r.max_count_off,
                         wAPLKC_sha256=r.wAPLKC_sha256, bouton_share=r.get('bouton_share'), comp_share=r.get('comp_share'),
                         realism_flag=r.realism_flag, host=r.get('host'), subdir='runs'))
        o = load_run99(r.tag)
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    long_path = out / 'model_pair_distances_long.csv'
    keep = []
    if long_path.exists():
        old = pd.read_csv(long_path)
        keep = [old[old.tag.isin(set(_fits(out).tag))]]
    if keep or parts:
        pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def figure(out: Path) -> Optional[Path]:
    """E99_settings.png: S's mean APL effect under the calibrated readout per setting (the flies' dotted; the uniform
    component's and the calyx off-claw component's at d = 1 dashed) and the pattern part against all recordings and
    halves A and B (the unchanged fits dashed); hollow: a setting with a run outside the bands or at a bracket's bound."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not (out / 'readings.csv').exists():
        return None
    rd = pd.read_csv(out / 'readings.csv')
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal') & (rd.variant != REF)]
    ref = rd[(rd.unit == 'S') & (rd.readout == 'cal') & (rd.variant == REF)]
    if g.empty:
        return None
    order = [v for b in BLOCKS.values() for v in b] + [v for vs in CORNERS.values() for v in vs] + [D0]
    g = g.set_index('variant').reindex([v for v in order if v in set(g.variant)]).reset_index()
    x = np.arange(len(g))
    bad = ((g.n_flagged > 0) | (g.n_at_bound > 0)).to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    ax = axes[0]
    ax.scatter(x, g.mean_effect, s=26, facecolors=['none' if b else '0.2' for b in bad], edgecolors='0.2')
    ax.axhline(FLIES_MEAN, color='k', ls=':', lw=0.8, label="the flies' mean APL effect")
    for v, c in ((UNIFORM[1.0], '#2e86ab'), (OFFCLAW[1.0], '#d1495b')):
        if v in set(g.variant):
            ax.axhline(float(g[g.variant == v].mean_effect.iloc[0]), color=c, ls='--', lw=0.7, label=v)
    ax.set_ylabel('mean APL effect (calibrated readout)')
    ax.legend(fontsize=7)
    ax = axes[1]
    cols = {'halfA': '#d1495b', 'halfB': '#2e86ab', 'z8': '0.2'}
    for t, c in cols.items():
        ax.scatter(x, g[f'{t}_pattern'], s=26, facecolors=['none' if b else c for b in bad], edgecolors=c,
                   label='all recordings' if t == 'z8' else f'half {t[-1]}')
        if len(ref):
            ax.axhline(float(ref[f'{t}_pattern'].iloc[0]), color=c, ls='--', lw=0.6)
    ax.set_ylabel('pattern part of the APL-effect loss')
    ax.legend(fontsize=7)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(g.variant, rotation=70, ha='right', fontsize=7)
    fig.suptitle("S's population, the bouton path at a tenth (hollow: outside the bands or at a bound)", fontsize=9)
    fig.tight_layout()
    path = out / 'E99_settings.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _md(df: pd.DataFrame, fmt: str = '.4f') -> str:
    return E94._md(df, fmt)


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    L = ["# E99: the evenness of APL's weights among KCs and over a KC's claws", '',
         f'Generated by experiments/E99_apl_weight_evenness.py at {H.code_commits().get("queue_commit", "")[:7]}.', '',
         "What is the connectome's and what is not: K1 to K4 and W2 are constructions that change one property of the "
         'counts (not candidates; none is chosen, called an equal or fitted). W1, T1 and T2 are counts as counted, spread '
         'evenly over a KC\'s claws (where on the dendrite a synapse lies is not modelled); T1 counts the peduncle\'s '
         'synapses as acting on the KC\'s drive, T2 does not. The weighted part of the weights keeps the family\'s '
         'breadth weighting. The counts are one hemibrain\'s.', '']
    if (out / 'identity.json').exists():
        s = _read_json(out / 'identity.json')
        L += ['## Part 0: the identities', '', f"passes: {s.get('passes')}", '']
        for k, x in s.items():
            if not isinstance(x, dict) or k == 'commits':
                continue
            keep = {kk: x.get(kk) for kk in ('passes', 'means_1', 'k1_k2_is_offclaw', 'w1_w2_is_connectome', 'k1_group_means_1',
                                              'k3_all_is_offclaw', 'k3_none_is_k2', 'sums_are_counts', 'w2_sum_is_n',
                                              'T1_is_e10', 'T2_is_e10', 'tables', 'apl_factor', 'apl_factor_saved',
                                              'factor_rel_diff', 'same_weights_as_saved') if kk in x}
            if isinstance(x.get('counts'), dict):
                keep['counts_differ'] = x['counts'].get('counts_differ')
            L.append(f'- {k}: {json.dumps(keep, default=str)}')
        L.append('')
    if (out / 'a1_shapes.csv').exists():
        d = pd.read_csv(out / 'a1_shapes.csv')
        d = d[d.wiring == 'real'][['setting', 'variant', 'kind', 'between', 'among', 'within', 'total', 'per_claw_cv',
                                   'claws_zero_share', 'kcs_zero_share', 'loglog_slope_n', 'within_cv_median', 'r_b']]
        L += ['## A1: the shapes (the real wiring; the same on every wiring but r with b_k)', '', _md(d, '.3f'), '']
    if (out / 'a2_predictions.json').exists():
        a2 = _read_json(out / 'a2_predictions.json')
        L += [f"## A2: the predictions (written {a2.get('written_at')})", '',
              f"(i) mean = {a2['descriptions']['i']['coef']['intercept']:+.4f} {a2['descriptions']['i']['coef']['total']:+.4f}"
              f" x total, rms residual {a2['descriptions']['i']['rms']:.4f}; (ii) mean = "
              f"{a2['descriptions']['ii']['coef']['intercept']:+.4f} "
              + ' '.join(f"{a2['descriptions']['ii']['coef'][c]:+.4f} x {c}" for c in PARTS3)
              + f", rms residual {a2['descriptions']['ii']['rms']:.4f}", '']
        if (out / 'a2_predictions.csv').exists():
            L += [_md(pd.read_csv(out / 'a2_predictions.csv'), '.4f'), '']
        for k in ('i', 'ii'):
            c = a2['classes'][k]
            L.append(f"- predicted by ({k}): between {c['between']['cls']} (slope {c['between']['slope_share']:.2f}, "
                     f"inequality {c['between']['ineq_share']:.2f}); with the peduncle's {c['replicate_peduncle']['cls']}; "
                     f"within {c['within']['cls']} (placement {c['within']['placement_share']:.2f}, amounts "
                     f"{c['within']['amounts_share']:.2f}); none: {c['none']['verdict']} ({c['none']['lost_share']:.2f}); "
                     f"K3/K2 {c['bound']['r3']:.2f}")
        L.append('')
    for name, title in (('a3_fifths.csv', 'A3: who responds (saved runs at size 1)'),
                        ('a3_claws.csv', "A3: the most and least weighted claws (D2, B6)"),
                        ('a4_yardstick.csv', "A4: the CV among KCs with the same number of claws")):
        if (out / name).exists():
            d = pd.read_csv(out / name)
            if name == 'a3_fifths.csv':
                d = d[(d.unit == 'S') & d.by.isin(['fifth_offclaw', 'claw_class', 'all'])]
            if name == 'a3_claws.csv':
                d = d[d.unit == 'S']
            L += [f'## {title}', '', _md(d, '.3f'), '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        runs = runs_table(out)
        g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
        graw = rd[(rd.unit == 'S') & (rd.readout == 'raw')].drop_duplicates('variant').set_index('variant')
        order = [REF] + [v for b in BLOCKS.values() for v in b] + [v for vs in CORNERS.values() for v in vs] + [D0]
        g = g.drop_duplicates('variant').set_index('variant').reindex([v for v in order if v in set(g.variant)]).reset_index()
        rr = runs.groupby('variant')
        L += ["## S's population (calibrated readout; the parts of the pattern part against all recordings)", '',
              '| setting | mean APL effect | mean, raw | pattern, half A | half B | all | six | amplitude | shape | spread ratio | '
              'Pearson r | ρ | loss | raw-count loss | factor | share at the boutons / component / weighted | flagged / '
              'at bound |', '|---|' + '---:|' * 16]

        def rng(q, c, spec='.3g'):
            if not len(q) or c not in q or not pd.to_numeric(q[c], errors='coerce').notna().any():
                return '–'
            v = pd.to_numeric(q[c], errors='coerce')
            return f"{v.min():{spec}}-{v.max():{spec}}"
        for r in g.itertuples():
            q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()
            sh = ' / '.join(rng(q, c, '.3f') for c in ('bouton_share_d', 'comp_share', 'weighted_share'))
            L.append(f'| {r.variant} | {r.mean_effect:+.3f} | {graw.mean_effect.get(r.variant, np.nan):+.3f} | '
                     f'{r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | {r.z6_pattern:.4f} | '
                     f'{r.z8_amplitude:.4f} | {r.z8_shape:.4f} | {r.z8_sd_ratio:.2f} | {r.z8_pearson:.2f} | {r.z8_rho:.2f} | '
                     f'{r.z8_loss:.4f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | {rng(q, "factor")} | {sh} | '
                     f'{r.n_flagged} / {r.n_at_bound} |')
        L.append('')
        mine = runs[runs.source == 'E99']
        if len(mine):
            L += ['## Realism of the runs (rule of 2026-11-10 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | '
                  'K | share at the boutons / component / weighted | APL->KC weight p90/p10 over responders | flag | flag '
                  '(earlier rule) |', '|---|' + '---:|' * 9 + '---|---|']
            for r in mine.itertuples():
                sh = ' / '.join(E94._f(getattr(r, c, np.nan), '.3f') for c in ('bouton_share_d', 'comp_share', 'weighted_share'))
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | '
                         f'{r.opk_on_tv:.3f} / {r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | '
                         f'{E94._f(getattr(r, "K_run", np.nan), ".3g")} | {sh} | '
                         f'{E94._f(getattr(r, "wkc_all_resp_p90_over_p10", np.nan), ".3g")} | '
                         f'{r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
    if (out / 'c_verdicts.json').exists():
        cv = _read_json(out / 'c_verdicts.json')
        L += ['## Part C: the verdicts (mean APL effect of S\'s population, calibrated readout)', '']
        for k in ('between_calyx_d1', 'between_calyx_d2', 'between_peduncle_d1', 'within_d1'):
            x = cv.get(k, {})
            s = x.get('S', {})
            L.append(f"- {k}: {s.get('cls')} ({json.dumps({kk: s.get(kk) for kk in s if kk != 'cls'}, default=str)}); "
                     f"by wiring {x.get('per_wiring')} -> {x.get('reported_as')}"
                     + (f"; agrees with the class: {x.get('agrees_with_class')}" if 'agrees_with_class' in x else ''))
        L.append(f"- none: {json.dumps(cv.get('none'), default=str)}")
        L.append(f"- bound: {json.dumps({k: v for k, v in cv.get('bound', {}).items() if k != 'yardstick'}, default=str)}")
        L.append(f"- vs the predictions: {json.dumps(cv.get('vs_predictions'), default=str)}")
        L.append('')
        L += ['### T1 and T2 against B6', '', _md(pd.DataFrame(cv.get('t_vs_b6', [])), '.4f'), '']
    for name in ('b_fifths.csv', 'line.csv', 'choice.csv', 'gate.csv', 'equal.csv', 'vs_uniform.csv', 'd_scores.csv',
                 'steps.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                L += [f'## {name}', '', '(none)', '']
                continue
            if name == 'b_fifths.csv':
                d = d[(d.unit == 'S') & (d.by == 'fifth_offclaw')]
            if name == 'line.csv':
                d = d[[c for c in d.columns if c in ('variant', 'unit', 'readout', 'ester_rel', 'alcohol_rel', 'resid',
                                                       'off_toward_flies', 'off_toward_flies_nr')]]
            if name == 'steps.csv' and len(d) > 200 and 'read' in d:
                d = d[d.read.astype(bool)]
            L += [f'## {name}', '', _md(d), '']
    if (out / 'd_rule.json').exists():
        r = _read_json(out / 'd_rule.json')
        L += ["## Part D's rule", '', _md(pd.DataFrame(r.get('rows', [])), '.4f'), '',
              f"runs: {r.get('runs')}; component {r.get('component')}; best {r.get('best')}", '']
    if (out / 'verdict.json').exists():
        vc = _read_json(out / 'verdict.json')
        short = {k: vc.get(k) for k in ('classes', 'verdict', 'gate', 'equal_settings', 'D_rule', 'D', 'flagged_runs',
                                        'at_bound_runs')}
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json's dated records: the commits, olfsysm's build, the model's provenance with the lab's input
    for the simulation parts, the settings and `extra`."""
    import olfsysm
    path = out / f'provenance_{part}.json'
    old = _read_json(path)
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('0', 'B', 'C', 'D', 'Dcheck'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    old.update(script='experiments/E99_apl_weight_evenness.py', part=part, runs=runs,
               settings=dict(BLOCKS={k: list(v) for k, v in BLOCKS.items()}, CORNERS={k: list(v) for k, v in CORNERS.items()},
                             APL_BRACKET=APL_BRACKET, D_DOSE=D_DOSE, GROUP_MIN=GROUP_MIN, K3_FRAC=K3_FRAC, T_MODE=T_MODE,
                             CLASS_HI=CLASS_HI, CLASS_LO=CLASS_LO, BY_WIRING_MIN=BY_WIRING_MIN, BOUND_LOSS=BOUND_LOSS,
                             NONE=dict(HI=NONE_HI, LO=NONE_LO), CHOICE_MIN_GAIN=CHOICE_MIN_GAIN,
                             EQUAL=dict(PATTERN=EQUAL_PATTERN, RAW=EQUAL_RAW, MEAN=EQUAL_MEAN), D_MEAN_MIN=D_MEAN_MIN,
                             D_VS_UNIFORM=D_VS_UNIFORM, D_T2_TIE=D_T2_TIE, D_FLAT_BOUNDS=D_FLAT_BOUNDS, D_MAXFEV=D_MAXFEV,
                             D_TAGS=D_TAGS, D_AS_GOOD=D_AS_GOOD, SHAPE_NAMES=SHAPE_NAMES,
                             counts_source='pyolfsysm connectome_APL_KC_synapse_counts (APL>KC table of 2025-12-14, '
                                           'APL_R, both sites in the ROI), through E96\'s wiring_counts'))
    _write_json(path, old)


def main() -> None:
    import warnings
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'A1', 'A2', 'A3', 'A4', 'B', 'score', 'C', 'Dcheck', 'D',
                                                      'Dscore', 'E', 'steps', 'variants', 'table'))
    ap.add_argument('--block', default=None, help='B: I, II or I,II (default both)')
    ap.add_argument('--wirings', default=None, help="'S', 'T', 'all8' or a comma-separated list (0 and B: S; D: a lane)")
    ap.add_argument('--pop', default='S', choices=('S', 'T'), help='D: which fits')
    ap.add_argument('--start', default='own', choices=tuple(D_TAGS), help="D: the start ('own': the wiring's best point of "
                    "part B; 'e96': E96's kept refit with the component exchanged)")
    ap.add_argument('--component', default=None, choices=T_SHAPES, help="D: the component (default the rule's)")
    ap.add_argument('--force', action='store_true', help="D: run fits the entry's rule does not allow")
    ap.add_argument('--force-a2', action='store_true', help="A: rewrite A2's predictions although part B has rows")
    ap.add_argument('--dry-run', action='store_true', help='D: print the rule, the starts and the commands only')
    ap.add_argument('--workers', type=int, default=4, help='score: co-response processes')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    ap.add_argument('--no-zcell', action='store_true', help="score: skip E73's per-KC readings")
    ap.add_argument('--variants', default=None, help='steps: comma-separated variants (default: chosen and equal settings)')
    a = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = out_dir()
    t0 = time.time()
    if a.part == '0':
        ws = _wirings(a.wirings or 'S')
        provenance(out, '0', wirings=list(ws))
        part_0(out, ws)
    elif a.part == 'A':
        provenance(out, 'A')
        part_a(out, force_a2=a.force_a2)
    elif a.part == 'A1':
        a1_shapes(out)
    elif a.part == 'A2':
        a2_predictions(out, force=a.force_a2)
    elif a.part == 'A3':
        a3_who_responds(out)
    elif a.part == 'A4':
        a4_yardstick(out)
    elif a.part == 'B':
        ws = _wirings(a.wirings or 'S')
        blocks = [x for x in a.block.split(',') if x] if a.block else None
        provenance(out, 'B', wirings=list(ws), blocks=blocks)
        part_b(out, blocks, ws)
    elif a.part == 'score':
        score(out, a.workers, coresponse=not a.no_co, zcell=not a.no_zcell)
    elif a.part == 'C':
        provenance(out, 'C')
        part_c(out)
    elif a.part == 'Dcheck':
        provenance(out, 'Dcheck')
        d_check(out)
    elif a.part == 'D':
        if not a.dry_run:
            provenance(out, 'D', pop=a.pop, start=a.start, wirings=a.wirings, component=a.component)
        part_d(out, a.pop, a.start, None if a.wirings is None else _wirings(a.wirings), a.component, force=a.force,
               dry_run=a.dry_run)
    elif a.part == 'Dscore':
        d_scores(out)
        if (out / 'readings.csv').exists():
            verdict99(out)
        comparison_table(out)
    elif a.part == 'E':
        part_e(out)
    elif a.part == 'steps':
        part_steps(out, None if a.variants is None else [x for x in a.variants.split(',') if x])
        comparison_table(out)
    elif a.part == 'variants':
        write_variants(out)
    elif a.part == 'table':
        figure(out)
        comparison_table(out)
    print(f'E99 {a.part} done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
