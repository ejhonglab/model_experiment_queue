#!/usr/bin/env python3
"""E132: the KC-KC synapses with a cholinergic sign: recurrent excitation among KCs through the hemibrain's claw-to-claw
synapses, with a pooled inhibition and alone, read as a dose through zero. Steps 0 and 1 (no olfsysm change) and step 2
(on olfsysm's option for negative KC-KC weights, 45e5cf9).

The entry is item E132 of experiment_priorities.md (moved into the list 2026-11-23; claimed by hal 21:00). "The family" is
C2's breadth-2 APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S is real,
P1, P2, M1 and T is P3, P4, M2, M3. The family's bouton base pins multibouton_claws='first' (E10_joint_fit.base_kws).
"KC-KC at zero" is E110's part A (KC-KC 0 in both conditions, both response rates held at fixed scales).

The dose and the pooled neuron. A dose d is the KC-KC scale (`kc_kc_inh`, the simulator's 'wKCKC_scale') on the
loader's normalized hemibrain matrix HB (post x pre, mean row sum 1, `connectome_KCKC_weights`), in absolute units
(the family's fitted scales s are +23 to +46); a negative dose is excitation. The pooled neuron is E122's uniform neuron
U (every entry HB's total / N^2, the diagonal included), filtered with 50 ms as the matrix is; its scale b ('wU_scale')
is absolute, and b / s is E129's factor (the pooled neuron's weight as a fraction of the matrix's at the fitted scale).

How the pooled neuron goes in (the entry leaves the route to the code): as a rank-1 part of the KC-KC matrix, E122's
route. The simulator's KC-KC matrix is set before every call to d HB + b U (summed dense, then CSC, as fit_mb_model
builds it) with 'wKCKC_scale' 1 (`Runner`). At b 0 the matrix olfsysm receives is bit for bit the family's (s h + 0 is
s h, with HB's sparsity); at d 0 it is bit for bit E129's pooled arm's (b u + 0). The feedback-cell path of olfsysm
fa3543d (`mbc1_inh` with uniform cell weights) is the same in exact arithmetic but filters the pooled spikes once
instead of each KC's activity, so it is not E129's run in floating point; the check part measures it on 'real'
(`route_diagnostic`). The checks below must pass with 0 differing spike counts before any run.

The route's cost and the alternative (2026-11-23, hal). The feedback-cell route also gave 0 differing spike counts
against E129's S1_P_real at E129's factor, but only empirically, at one point; the matrix route is kept because its
input to olfsysm is identical by construction, so U0's whole bisection follows E129's. Its price is speed: U is dense
(N^2 entries), and a simulation with the pooled neuron present takes about 2.5 times the family's (35 s against 14 s
on a quiet hal at 16 threads, 85 s against 33 s at a load of about 48), memory-bound, so two lanes at once do not
halve the total.

Readings the entry leaves open, fixed here before any run:
1. The pooled hold (step 1's arms; E129's, extended to both conditions; `pooled_hold`): the threshold stays the fit's.
   With APL off (APL->KC and APL->bouton 0) the pooled scale b = s x f is bisected (E69's `bisect_factor` through
   pyolfsysm.tuning.hold_rates with thr_key 'wU_scale', from f 1, E129's SCALE_BRACKET 1/50 to 50 and SCALE_STEP 1.25)
   to the family's APL-off count; then the tied APL scales (APL->KC, KC->APL, bouton->APL; E83.TIED_ON) are bisected
   with APL on (E77's APL_BRACKET 0.1 to 10, APL_STEP 1.25) to the family's APL-on count, the pooled neuron at that b.
   U0's first stage is E129's step-1a bisection of arm P exactly, so its APL-off run must equal E129's S1_P_<w> run
   (0 differing spike counts, checked in the run; RuntimeError otherwise).
2. The threshold hold (step 2's X arms; E110's part A): no pooled neuron; E77_input_screen.hold_rates (the threshold
   bisected with APL off, then the tied APL scales with APL on).
3. A corner: a bisected scale at a bound of its bracket, or a response count more than RATE_TOL (0.2 %) from its target
   (E129's `held_ok`). An arm with corners on more than CORNER_WIRINGS_MAX (two) wirings gets no class.
4. The classes (E123's, E126's and E131's, unchanged; calibrated pattern part, S's and T's populations apart, against
   all recordings and each half A to D; `classify`): "below the family" (lower by WITHIN 0.003 or more on both
   populations and on 6 of 8 wirings on every target, the mean APL effect under the readout within MEAN_TOL 0.02 of the
   flies' +0.178, the raw-count loss not higher by RAW_LOSS_TOL 0.005 on both populations, every run in the bands);
   "equal" (E123's, one-sided); "costs" (higher by more than WITHIN on both populations against all recordings);
   "neither" otherwise; "corner (no class)" by reading 3. "Keeps the gain" (E123's flag): at or below the half-way point
   between the family and KC-KC at zero (0.0554 on S, 0.0577 on T against all recordings) on every target.
5. Two runs "differ" only by E126's rule (E129's `differ`): one lower on 7 of 8 wirings against all recordings and on
   both populations against all five targets. Step 1 reads each mirror dose against U0 (its dose 0) and against the
   family; the sign's own effect (a negative dose against its mirror and its dose 0) is step 2's.
6. The stop rule (`stop_rule`, read on S): a dose passes if (i) S's population against all recordings is not above the
   family's by more than STOP_SLACK (0.005), (ii) the largest spike count per KC and odor with APL off is at most
   MAX_COUNT_FACTOR (2) times the family's on each of the eight wirings, (iii) every run is in the four bands with
   corners on at most two wirings, and (iv, added on moving) it is lower than its own dose 0 (U0 for the pooled arms,
   E110's part A for the arm alone) on S's population against all recordings and on 3 of S's 4 wirings. The mirror's
   doses are read by it and reported; only step 2's negative doses are candidates for step 3.
7. Where the excitation would land (step 0, a description; `landing`): each KC's KC-KC input per odor is E129's step-0
   reader (s HB @ n / 0.5 s, Hz; n the family's saved spike counts of the condition), with APL on and off; responders
   are the saved responses (two or more spikes) of the same condition. For an ordered pair (a, b): the mean input under
   odor a onto the KCs that respond to b and not to a, over the mean input under a onto the KCs that respond to
   neither; a pair's ratio is the mean of its two directed ratios. Group means over E90's groups of the megamat pairs
   (ester 15, alcohol 6, ester with alcohol 24, other 21) and over all 120 pairs, per wiring and for S4, T4 and all8 as
   the mean of their wirings' group means.
8. The matrix's side of stability (step 0; `stability`): the largest real eigenvalue (the Perron root; HB is
   nonnegative) of HB and of each KC type's diagonal block (γ, α′β′, αβ), and how much of the leading right and left
   eigenvectors is uniform over the KCs (cos² of their angle with the all-ones vector), the part a pooled neuron can
   cancel; beside it each wiring's s times the root and the doses 2, 5 and 10 times it.
9. The slope toward zero (step 0; `slopes`): under the threshold hold, the family (s), the matrix at s / 2 (E131's
   H50: fixed scales, both rates re-held) and KC-KC at zero (E110's part A). Per wiring, each reading's change per unit
   of scale on the segments s -> s / 2 and s / 2 -> 0; per population the same over the mean of its wirings' s. Carried
   through zero (the expectation on record for the excitation alone) the half->zero segment's slope: r(0) - slope x d at
   d 2, 5, 10 (doses -2, -5, -10).
10. Late firing (`late_share`): the share of spikes at or after the pulse's end (0.5 s), with APL on and off. The
    family's from E121's step-0 rebuild (`E121_onset_race.run_wiring`, committed code; part late); each step-1 run's
    from its held point simulated again with spike times, which must reproduce the held counts.
11. Sizes (per run): b / s; the threshold; the mean over KCs and odors of each part's input (d HB @ n / 0.5 s and
    b U @ n / 0.5 s, Hz) with APL on and off, over the fit's threshold; the largest count per KC and odor and the KCs
    with 20 or more spikes to some odor (E110's run_extras).

Step 2 (the excitation; Yang's OK for the option 2026-11-23 21:25; code 2026-11-24). It follows readings 1 to 6 as they
stand; what is written below is how they are carried to negative doses, nothing re-planned:
12. The option. olfsysm 45e5cf9 adds `mp.kc.allow_negative_kc_kc` (C++, default false) and fit_mb_model's
    `allow_negative_kc_kc`. A step-2 runner is built with it on (`Runner(w, allow_negative=True)`: E83's
    `build_simulator(w, extra_kws=dict(allow_negative_kc_kc=True))`, its mp's flag checked true); a step-1 runner without
    it (the flag checked false where the build has it). The runner still feeds d HB + b U with 'wKCKC_scale' 1, so only
    olfsysm's own check (`check_kc_kc_inh`) sees the negative entries; `combined_matrix` lets them through only for a
    runner with the flag, b stays >= 0, and a negative dose on a build without the option raises SignedDoseError
    (`negative_kckc_available`). Nothing negative is bisected: the pooled hold bisects b, the threshold hold the
    threshold, both the tied APL scales, each of which lowers the response count as it rises (E69's `bisect_factor`).
13. The flag check (part flagcheck, before any step-2 run): on 'real', a runner with the flag on and non-negative
    matrices reproduces the family's saved run (matrix s, pooled 0) and step 1's S1_UI2_real at its held parameters,
    APL on and off, with 0 differing spike counts and responses (flagcheck_<host>.json). Step 2 refuses until it has
    passed with a build that has the option, and until step 1's checks have passed on every wiring.
14. UX2, UX5, UX10 (the matrix at -2, -5, -10 beside the pooled neuron) by the pooled hold (reading 1, from b = s);
    X2, X5, X10 (the matrix alone) by the threshold hold (reading 2, from the fit's threshold, d HB and no pooled
    neuron). Order inside a lane (`planned_step2`): -2 on 'real' first in both arms (UX2 then X2, the stability smoke)
    when the lane has 'real', then the UX arms wiring by wiring in the order the lane gives, then the X arms likewise;
    `--reverse` reverses the whole list (a lane joining from the other end meets the others last). Every run saves its
    spike times, as step 1's do (locally and in the archive, step2/), unless a condition has more than MAX_SPIKE_TIMES
    spikes (a runaway): then its counts only, and its late share is NaN, recorded in the row ('spike_times_saved').
15. Step 2's scoring (score; readings 3 to 6 unchanged): each arm's class and the stop rule (its dose 0 U0 for the UX
    arms, E110's part A for the X arms); the candidates are the negative doses that pass, and if several pass the one
    with the lowest pattern part on S's population against all recordings goes on (the entry's rule); if none passes
    the experiment stops after step 2. "Differ" (reading 5): each UX against its mirror UI at the same size and against
    U0 (the sign's own effect), each X against E110's part A and against the saved half matrix H50 (the entry's step
    2), and UX against X at the same dose (outcome (d)'s "not told apart"). Beside them, no gate: the corners and the
    thresholds under 60 Hz, the amplitude and shape parts, the spread ratio, Pearson r and ρ on the per-pair APL effect,
    the mean APL effect under the readout and on raw counts (the latter also minus the family's and minus the dose 0's,
    for outcome (b)), the raw-count loss, the largest counts, the KCs with 20 or more spikes, the late shares and γ's
    off / on ratio.
Steps 3 and 4 run only for a dose that passes (verdict.json 'step2'); their code is written then (`part_step3`,
`part_step4` say which).

Parts (CLI `--part`):
- check (simulations, `--wirings`): (i) the runner at the fit's parameters (matrix s, pooled 0) reproduces the
  wiring's saved family run with 0 differing spike counts and responses, APL on and off; (ii) at matrix 0 and the
  pooled neuron at E129's bisected factor (its runs table) the APL-off run equals E129's saved S1_P_<w> run; on 'real'
  (iii) the route diagnostic (the feedback-cell path at the same b, its differing counts against E129's run, recorded,
  not a gate). check_<host>_<lane>.json; the runs refuse until every wiring has passed (i) and (ii). A simulation with
  the pooled neuron present (U is dense, N^2 entries) takes about 35 s on hal at 16 threads, against 14 s for the
  family's sparse matrix.
- step0 (no simulation): step0_references.csv (readings of the references, REFS: the family, H50, E110's parts A to
  D, the uniform neuron at the whole weight, E129's pooled arm, E125's PO50; `references`), step0_ref_runs.csv,
  step0_slopes.csv, step0_landing_pairs.csv, step0_landing.csv, step0_stability.csv, step0.json,
  provenance_step0.json.
- late (simulations): the family's late firing from E121's rebuild, late_family.csv (spike times under
  outputs/E132_kckc_excitation/late/).
- runs (simulations, resumable, per lane `--wirings`; a point claimed by mkdir of results/E132_kckc_excitation/claims/
  <tag>/, `--reclaim` to take over a dead process's): step 1's U0, UI2, UI5, UI10. Rows in runs_<host>_<lane>.csv,
  provenance in provenance_<host>_<lane>.json, runs and spike times saved locally (outputs/E132_kckc_excitation/step1/)
  and in the model-output archive (data/model_outputs/E132_kckc_excitation/step1/).
- flagcheck (simulations, `--wirings`, default real): reading 13, flagcheck_<host>.json.
- step2 (simulations, resumable, per lane `--wirings`, `--arms`, `--reverse`; claims as for runs): step 2's UX2, UX5,
  UX10, X2, X5, X10. Rows in runs_<host>_step2-<lane>.csv, provenance in provenance_<host>_step2-<lane>.json, runs and
  spike times under outputs/E132_kckc_excitation/step2/ and data/model_outputs/E132_kckc_excitation/step2/.
- score (no simulation): readings.csv, runs.csv, verdict.json, comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts, tagged by the run tag, the archive's file
  name, so kcq/apl_effect.py scores the realism bands).
- plan: the points a lane would run (`--step 2` for step 2's).

The holds (CLAUDE.md, "Connectome data"): only the model's hemibrain wiring and E72's seven re-draws of it, and the
loader's hemibrain KC->KC counts and KC types; no olfsysm change but the option of reading 12; nothing is published.

Usage (from the repository root; 16 threads, under nice):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part check
    ./venv/bin/python -u experiments/E132_kckc_excitation.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part late
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part runs --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part runs --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E132_kckc_excitation.py --part score && ./venv/bin/python -u experiments/E132_kckc_excitation.py --part variants
  step 2 (on a build with the option):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part flagcheck
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part step2 --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E132_kckc_excitation.py --part step2 --wirings P3 P4 M2 M3
    (another host, from the other end:) ... --part step2 --wirings real P1 P2 M1 P3 P4 M2 M3 --reverse
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import shutil
import socket
import time
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E132_kckc_excitation'
SCRIPT = 'experiments/E132_kckc_excitation.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
C_WIRINGS = ('real', 'P1', 'P3', 'M2')
POPS = ('S4', 'T4')
POP_WIRINGS = {'S4': S_WIRINGS, 'T4': T_WIRINGS, 'all8': WIRINGS, 'C4': C_WIRINGS}
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
CONDS = ('on', 'off')
# arm -> step, matrix dose d (absolute KC-KC scale), pooled neuron, hold
ARMS: Dict[str, dict] = {
    'U0': dict(step=1, dose=0.0, pooled=True, hold='pooled'),
    'UI2': dict(step=1, dose=2.0, pooled=True, hold='pooled'),
    'UI5': dict(step=1, dose=5.0, pooled=True, hold='pooled'),
    'UI10': dict(step=1, dose=10.0, pooled=True, hold='pooled'),
    'UX2': dict(step=2, dose=-2.0, pooled=True, hold='pooled'),
    'UX5': dict(step=2, dose=-5.0, pooled=True, hold='pooled'),
    'UX10': dict(step=2, dose=-10.0, pooled=True, hold='pooled'),
    'X2': dict(step=2, dose=-2.0, pooled=False, hold='threshold'),
    'X5': dict(step=2, dose=-5.0, pooled=False, hold='threshold'),
    'X10': dict(step=2, dose=-10.0, pooled=False, hold='threshold'),
}
STEP1_ARMS = ('U0', 'UI2', 'UI5', 'UI10')
STEP2_ARMS = ('UX2', 'UX5', 'UX10', 'X2', 'X5', 'X10')
MIRROR = {'UX2': 'UI2', 'UX5': 'UI5', 'UX10': 'UI10'}
ARM_LABELS = {
    'U0': "the pooled neuron alone in the matrix's place (matrix 0), pooled hold",
    'UI2': 'the matrix at +2 beside the pooled neuron, pooled hold',
    'UI5': 'the matrix at +5 beside the pooled neuron, pooled hold',
    'UI10': 'the matrix at +10 beside the pooled neuron, pooled hold',
    'UX2': 'the matrix at -2 (excitation) beside the pooled neuron, pooled hold',
    'UX5': 'the matrix at -5 (excitation) beside the pooled neuron, pooled hold',
    'UX10': 'the matrix at -10 (excitation) beside the pooled neuron, pooled hold',
    'X2': 'the matrix at -2 (excitation) alone, threshold hold',
    'X5': 'the matrix at -5 (excitation) alone, threshold hold',
    'X10': 'the matrix at -10 (excitation) alone, threshold hold',
}
# references read in step 0 and score: kind -> what it is
REFS = {
    'family': "the family's saved fits (matrix at s)",
    'H50': "the matrix at s / 2, both rates re-held (E131's H50)",
    'E110_D': "the matrix at s / 2 held, the rest refitted (E110's part D; S's wirings)",
    'E110_A': 'KC-KC 0, both rates re-held (E110 part A)',
    'E110_B': "KC-KC 0, the rest refitted (E110's part B, counted start)",
    'E110_C': "the return fits from KC-KC 2 (E110's part C; real, P1, P3, M2)",
    'U': "the uniform neuron at the matrix's whole weight, both rates re-held (E122 step 1 on S, E123 on T)",
    'E129_P': "E129's pooled arm: the family's APL-on run, the pooled neuron in the matrix's place with APL off",
    'PO50': "E125's pooled survivor at 0.5 (PO50)",
}
# the "differ" comparisons of step 1 (a, b); descriptions, no gate
COMPARISONS = (('UI2', 'U0'), ('UI5', 'U0'), ('UI10', 'U0'), ('UI10', 'UI2'), ('U0', 'E129_P'), ('U0', 'U'),
               ('U0', 'E110_A'), ('UI10', 'H50'))
# step 2's (reading 15): each UX against its mirror and against U0, each X against E110's part A and H50, UX against X
STEP2_COMPARISONS = tuple((x, MIRROR[x]) for x in ('UX2', 'UX5', 'UX10')) + tuple(
    (x, 'U0') for x in ('UX2', 'UX5', 'UX10')) + tuple((x, 'E110_A') for x in ('X2', 'X5', 'X10')) + tuple(
    (x, 'H50') for x in ('X2', 'X5', 'X10')) + (('UX2', 'X2'), ('UX5', 'X5'), ('UX10', 'X10'))
STEP2_DOSES = {2: ('UX2', 'X2'), 5: ('UX5', 'X5'), 10: ('UX10', 'X10')}
SMOKE = (('UX2', 'real'), ('X2', 'real'))
FLAG_CHECK_ARM = 'UI2'
MAX_SPIKE_TIMES = 20_000_000
PULSE_S = 0.5
STIM_END = 0.5
LANDING_DOSES = (2.0, 5.0, 10.0)
SCALE_BRACKET, SCALE_STEP = (1 / 50, 50.0), 1.25
RATE_TOL = 0.002
WITHIN = 0.003
N_WIRINGS_NEEDED = 6
MEAN_TOL = 0.02
RAW_LOSS_TOL = 0.005
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7
STOP_SLACK = 0.005
MAX_COUNT_FACTOR = 2.0
ENTRY_FAMILY_MAX_COUNT_OFF = 72
DOSE0_WIRINGS = 3
REALISM_THR_HZ = 60.0
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS + ('wU_scale',)
ROUTE = ("rank-1 part of the KC-KC matrix: the simulator's KC-KC matrix d HB + b U (dense sum, CSC) with "
         "'wKCKC_scale' 1 (E122's route)")
# readings of the slope (name, readings.csv column, readout)
SLOPE_READINGS = (('pattern', 'all_pattern', 'cal'), ('amplitude', 'all_amplitude', 'cal'),
                  ('shape', 'all_shape', 'cal'), ('sd_ratio', 'all_sd_ratio', 'cal'),
                  ('mean_effect_cal', 'mean_effect', 'cal'), ('mean_effect_raw', 'mean_effect', 'raw'),
                  ('loss_raw', 'all_loss', 'raw'), ('rho_on', 'rho_on', 'raw'), ('rho_off', 'rho_off', 'raw'),
                  ('max_count_off', 'run_max_count_off', 'cal'), ('n_kc_ge20_off', 'run_n_kc_ge20_off', 'cal'),
                  ('ratio', 'run_ratio', 'cal'))
_CACHE: dict = {}


class SignedDoseError(NotImplementedError):
    """A negative KC-KC dose where it cannot run: on a build without olfsysm's option for negative KC-KC weights, on a
    runner built without the flag, or in a part for step 1's arms."""


class GateClosedError(NotImplementedError):
    """A step that runs only for a dose that passes the stop rule (steps 3 and 4)."""


def negative_kckc_available() -> bool:
    """Whether the installed olfsysm has the option for negative KC-KC weights (olfsysm 45e5cf9): the C++ flag
    `mp.kc.allow_negative_kc_kc` and fit_mb_model's `allow_negative_kc_kc`."""
    import inspect
    import olfsysm
    from pyolfsysm.mb_model import fit_mb_model
    return bool(hasattr(olfsysm.ModelParams().kc, 'allow_negative_kc_kc')
                and 'allow_negative_kc_kc' in inspect.signature(fit_mb_model).parameters)


def olfsysm_build() -> Dict[str, object]:
    """The installed olfsysm's commit, its dirty flags and whether it has the option."""
    import olfsysm
    bi = olfsysm.build_info()
    return dict(olfsysm_git_hash=bi.get('git_hash'), olfsysm_git_dirty=bi.get('git_dirty'),
                olfsysm_cpp_dirty=bi.get('cpp_dirty'), olfsysm_source_sha256=bi.get('source_sha256'),
                negative_kckc_option=negative_kckc_available())


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E132_kckc_excitation/ (created)."""
    return H.results_dir(EXPERIMENT)


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=str) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def lane_name(wirings: Sequence[str]) -> str:
    return '-'.join(wirings)


def tag_of(arm: str, w: str) -> str:
    return f"S{ARMS[arm]['step']}_{arm}_{w}"


def n_resp(responses: pd.DataFrame) -> int:
    return int((responses.to_numpy(dtype=float) > 0.5).sum())


def kc_types_of(kc_index: pd.Index) -> np.ndarray:
    return np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)


def odor_name(col) -> str:
    """'2h @ -3' -> '2h'."""
    return str(col).split(' @ ')[0]


def sha256_array(a: np.ndarray) -> str:
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=float), dtype='<f8').tobytes()).hexdigest()


def _load_from(experiment: str, part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """Another experiment's saved run ({'on', 'off'}): this host's outputs/<experiment>/<part>/, else the archive's
    data/model_outputs/<experiment>/<part>/ (no folder is created)."""
    for d in (H.HERE / 'outputs' / experiment / part, H.shared_dir() / 'data' / 'model_outputs' / experiment / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def _ndiff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> Dict[str, int]:
    return dict(counts=int((a[1].to_numpy(dtype=float) != b[1].to_numpy(dtype=float)).sum()),
                responses=int((a[0].to_numpy(dtype=float) != b[0].to_numpy(dtype=float)).sum()))


def _off(p: Mapping[str, float]) -> Dict[str, float]:
    return dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)


def _preload() -> None:
    """Imports every module the simulation parts use before E108's commit check."""
    from experiments import E10_joint_fit  # noqa: F401
    from experiments import E69_kc_firing_ceiling  # noqa: F401
    from experiments import E77_input_screen  # noqa: F401
    from experiments import E108_paa_extrapolation  # noqa: F401
    from experiments import E110_no_kckc  # noqa: F401
    from experiments import E123_kckc_form  # noqa: F401
    from experiments import E129_block_survivor_static  # noqa: F401


# ---------------------------------------------------------------------------
# the matrices and the runner

def require_dose(d: float, b: float = 0.0, allow_negative: bool = False) -> None:
    """Refuses a non-finite dose or pooled scale and a negative pooled scale (ValueError: a pooled neuron is an
    inhibition), and a negative dose unless `allow_negative` (SignedDoseError: it needs a runner built with olfsysm's
    option for negative KC-KC weights, reading 12)."""
    d, b = float(d), float(b)
    if not (np.isfinite(d) and np.isfinite(b)):
        raise ValueError(f'dose {d}, pooled scale {b}')
    if d < 0 and not allow_negative:
        raise SignedDoseError(f'a negative KC-KC dose ({d:g}) needs a runner built with the olfsysm option for negative '
                              'KC-KC weights (Runner(w, allow_negative=True); olfsysm 45e5cf9 or later)')
    if b < 0:
        raise ValueError(f'the pooled scale is {b:g}; a pooled neuron is an inhibition')


def require_nonnegative(d: float, b: float = 0.0) -> None:
    """`require_dose` without negative doses (step 1)."""
    require_dose(d, b, allow_negative=False)


def uniform_neuron(ref: np.ndarray) -> np.ndarray:
    """E122's U (E129's `uniform_neuron`, the same operations): every entry ref's total / N^2, the diagonal included."""
    n = ref.shape[0]
    return np.full((n, n), float(ref.sum()) / n / n)


def combined_matrix(HB: np.ndarray, U: np.ndarray, d: float, b: float, allow_negative: bool = False):
    """d HB + b U as a CSC matrix (post x pre), the KC-KC matrix of the route; refuses a negative pooled scale, and a
    negative dose or negative entries unless `allow_negative` (`require_dose`).

    At b 0 its entries and sparsity are those of the family's matrix (d HB + 0 = d HB exactly, also for d < 0: HB's
    zeros give -0.0 + 0.0 = 0.0, dropped); at d 0 they are b U's.
    """
    from scipy.sparse import csc_matrix
    require_dose(d, b, allow_negative)
    HB, U = np.asarray(HB, dtype=float), np.asarray(U, dtype=float)
    if HB.shape != U.shape or HB.shape[0] != HB.shape[1]:
        raise ValueError(f'matrices {HB.shape}, {U.shape}')
    if (HB < 0).any() or (U < 0).any():
        raise ValueError('HB and U must be nonnegative (the sign is the dose\'s)')
    M = float(d) * HB + float(b) * U
    if not np.isfinite(M).all() or (not allow_negative and (M < 0).any()):
        raise ValueError('the KC-KC matrix must be finite (and nonnegative without the option)')
    return csc_matrix(M)


class Runner:
    """One wiring's family fit as a fixed-parameter simulator (E83's `build_simulator` at the fit) whose KC-KC matrix is
    d HB + b U: params 'wKCKC_scale' is the dose d and 'wU_scale' the pooled scale b (0 when absent); every other
    parameter as E83's simulator takes it. Before each call the simulator's unscaled KC-KC matrix (its closure cell
    'wKCKC_unscaled') is replaced by `combined_matrix(HB, U, d, b)` and 'wKCKC_scale' passed as 1, then restored.

    With `allow_negative` (step 2, reading 12) the simulator is built with olfsysm's option for negative KC-KC weights
    on (E83's `build_simulator(w, extra_kws=dict(allow_negative_kc_kc=True))`), and negative doses are let through;
    without it, the flag is checked off where the build has it.

    Attributes: w, allow_negative, p (the fit's parameters), s (its KC-KC scale), saved, target_on, target_off, size,
    kc_index, types, HB (the loader's matrix, checked equal to the simulator's), U, hb_sha256, u_sha256, mp, rv, sim,
    info."""

    def __init__(self, w: str, allow_negative: bool = False):
        from experiments import E10_joint_fit as J
        from experiments import E129_block_survivor_static as E129
        _preload()
        if allow_negative and not negative_kckc_available():
            raise SignedDoseError('the installed olfsysm has no option for negative KC-KC weights (olfsysm 45e5cf9: '
                                  "mp.kc.allow_negative_kc_kc); pull olfsysm and run ./venv/bin/pip install -e ../olfsysm")
        P = E83.wiring_parts(w)
        self.w = w
        self.allow_negative = bool(allow_negative)
        self.p = {k: float(P['p'][k]) for k in E83.PARAM_KEYS}
        self.s = float(self.p['wKCKC_scale'])
        self.saved = P['saved']
        self.target_on, self.target_off, self.size = int(P['target_on']), int(P['target_off']), int(P['size'])
        self.kc_index = self.saved['on'][1].index
        self.types = kc_types_of(self.kc_index)
        if self.allow_negative:
            self.sim, self.info = E83.build_simulator(w, extra_kws=dict(allow_negative_kc_kc=True))
        else:
            self.sim, self.info = E83.build_simulator(w)
        cells = self.info['cells']
        self._cell = cells['wKCKC_unscaled']
        self._W0 = self._cell.cell_contents
        self.mp, self.rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        flag = getattr(self.mp.kc, 'allow_negative_kc_kc', None)
        if (flag is not None or self.allow_negative) and bool(flag) != self.allow_negative:
            raise RuntimeError(f"{w}: the simulator's mp.kc.allow_negative_kc_kc is {flag}, not {self.allow_negative}")
        HB = E129.hemibrain_matrix(self.kc_index).to_numpy(dtype=float)
        W0 = self._W0.toarray() if hasattr(self._W0, 'toarray') else np.asarray(self._W0, dtype=float)
        if not np.array_equal(W0, HB):
            raise RuntimeError(f"{w}: the simulator's KC-KC matrix is not the loader's hemibrain matrix")
        self.HB, self.U = HB, uniform_neuron(HB)
        self.hb_sha256, self.u_sha256 = sha256_array(HB), sha256_array(self.U)
        if float(self.mp.kc.kc_kc_tau) != float(J.TAU):
            raise RuntimeError(f'{w}: kc_kc_tau is {self.mp.kc.kc_kc_tau}, not {J.TAU}')
        self._key, self._M = None, None

    def matrix(self, d: float, b: float):
        key = (float(d), float(b))
        if key != self._key:
            self._M = None
            self._M = combined_matrix(self.HB, self.U, d, b, self.allow_negative)
            self._key = key
        return self._M

    def params(self, d: float, b: float, **kw) -> Dict[str, float]:
        return dict(self.p, wKCKC_scale=float(d), wU_scale=float(b), **kw)

    def __call__(self, params: Mapping[str, float], spike_times: bool = False):
        d, b = float(params['wKCKC_scale']), float(params.get('wU_scale', 0.0))
        M = self.matrix(d, b)
        q = {k: float(v) for k, v in params.items() if k != 'wU_scale'}
        q['wKCKC_scale'] = 1.0
        self._cell.cell_contents = M
        try:
            return self.sim(q, spike_times=spike_times)
        finally:
            self._cell.cell_contents = self._W0

    def close(self) -> None:
        self._M = None
        self._key = None
        self.sim = None
        self.info = None
        gc.collect()


def get_runner(w: str, allow_negative: bool = False) -> Runner:
    """The Runner of a wiring (with or without the option's flag), cached (one at a time)."""
    c = _CACHE.get('runner')
    if c is not None and c.w == w and c.allow_negative == bool(allow_negative):
        return c
    if c is not None:
        c.close()
    _CACHE.pop('runner', None)
    gc.collect()
    _CACHE['runner'] = Runner(w, allow_negative=allow_negative)
    return _CACHE['runner']


def held_ok(n: int, target: int) -> bool:
    """Reading 3: the response count within RATE_TOL of the target."""
    return bool(abs(int(n) - int(target)) <= RATE_TOL * int(target))


def _finish_hold(res: dict, first: str) -> dict:
    res['held_on'] = held_ok(res['n_on'], res['target_on'])
    res['held_off'] = held_ok(res['n_off'], res['target_off'])
    res['corner'] = bool(res['at_bound_thr'] or res['at_bound_apl'] or not (res['held_on'] and res['held_off']))
    res['first_stage'] = first
    return res


def pooled_hold(runner: Runner, d: float, label: str) -> dict:
    """Reading 1: the pooled hold of one point (pyolfsysm.tuning.hold_rates with thr_key 'wU_scale' from b = s; the
    threshold the fit's). Adds 'pooled_factor' (b / s), 'pooled_scale', 'held_on', 'held_off', 'corner', 'seconds'."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E77_input_screen as E77
    from pyolfsysm import tuning
    require_dose(d, allow_negative=runner.allow_negative)
    t0 = time.time()
    p = runner.params(d, runner.s)
    res = tuning.hold_rates(runner, p, target_off=runner.target_off, target_on=runner.target_on, size=runner.size,
                            thr_key='wU_scale', off_params=dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                            tied_on=E83.TIED_ON, thr_bracket=SCALE_BRACKET, thr_step=SCALE_STEP,
                            apl_bracket=E77.APL_BRACKET, apl_step=E77.APL_STEP, label=label,
                            max_iters=E69.BISECT_ITERS)
    res['pooled_factor'] = float(res['thr_factor'])
    res['pooled_scale'] = float(res['params']['wU_scale'])
    res['seconds'] = time.time() - t0
    return _finish_hold(res, 'wU_scale')


def threshold_hold(runner: Runner, d: float, label: str) -> dict:
    """Reading 2: the threshold hold of one point (no pooled neuron; E77's brackets and steps, as E77_input_screen's
    `hold_rates` passes them to pyolfsysm.tuning.hold_rates), for step 2's arms alone; a negative dose needs a runner
    with the option (reading 12)."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E77_input_screen as E77
    from pyolfsysm import tuning
    require_dose(d, allow_negative=runner.allow_negative)
    t0 = time.time()
    p = runner.params(d, 0.0)
    res = tuning.hold_rates(runner, p, target_off=runner.target_off, target_on=runner.target_on, size=runner.size,
                            thr_key='fixed_thr', off_params=dict(wAPLKC_scale=0.0, wAPLPN_scale=0.0),
                            tied_on=E83.TIED_ON, thr_bracket=E77.THR_BRACKET, thr_step=E77.THR_STEP,
                            apl_bracket=E77.APL_BRACKET, apl_step=E77.APL_STEP, label=label,
                            max_iters=E69.BISECT_ITERS)
    res['pooled_factor'], res['pooled_scale'] = 0.0, 0.0
    res['seconds'] = time.time() - t0
    return _finish_hold(res, 'fixed_thr')


def late_share(spikes: pd.DataFrame, stim_end: float = STIM_END) -> float:
    """Reading 10: the share of spikes at or after the pulse's end (NaN without spikes)."""
    t = spikes['time_s'].to_numpy(dtype=float)
    return float((t >= stim_end - 1e-12).mean()) if t.size else float('nan')


def part_inputs(counts: pd.DataFrame, HB: np.ndarray, d: float, b: float, fixed_thr: float) -> dict:
    """Reading 11: the mean over KCs and odors of the matrix's input (d HB @ n / 0.5 s) and the pooled neuron's
    (b U @ n / 0.5 s), in Hz and over the fit's threshold, from spike counts n (KCs x odors)."""
    n = counts.to_numpy(dtype=float)
    N = n.shape[0]
    mat = float(d) * (np.asarray(HB) @ n) / PULSE_S
    pooled_per_odor = float(b) * float(np.asarray(HB).sum()) / N / N * n.sum(axis=0) / PULSE_S
    return dict(kckc_in_mean=float(mat.mean()), kckc_in_p90=float(np.quantile(mat, 0.9)),
                kckc_in_over_thr=float(mat.mean() / fixed_thr), pooled_in_mean=float(pooled_per_odor.mean()),
                pooled_in_over_thr=float(pooled_per_odor.mean() / fixed_thr))


def run_realism(runner: Runner, res: Mapping, d: float, max_spike_times: int = MAX_SPIKE_TIMES
                ) -> Tuple[dict, Dict[str, pd.DataFrame]]:
    """Readings 10 and 11 of a held point: simulated again with spike times (which must reproduce the held counts) for
    the late shares; the parts' inputs; returns (row, {'on', 'off': spike times}). A condition whose held run has more
    than `max_spike_times` spikes (a runaway) is not simulated again: no spike times, late share NaN,
    'spike_times_saved_<c>' False."""
    row, spikes = {}, {}
    for c, q in (('on', res['params']), ('off', res['params_off'])):
        n_held = int(res[c][1].to_numpy(dtype=float).sum())
        row[f'spike_times_saved_{c}'] = n_held <= max_spike_times
        if n_held > max_spike_times:
            row[f'late_share_{c}'] = float('nan')
            row[f'n_spikes_{c}'] = n_held
            s = res[c][1]
        else:
            r, s, spk = runner(dict(q), spike_times=True)
            dd = _ndiff((r, s), res[c])
            if dd['counts'] or dd['responses']:
                raise RuntimeError(f'the held point simulated again differs from the held run ({c}: {dd})')
            spikes[c] = spk
            row[f'late_share_{c}'] = late_share(spk)
            row[f'n_spikes_{c}'] = int(len(spk))
        for k, v in part_inputs(s, runner.HB, d, float(q.get('wU_scale', 0.0)), float(q['fixed_thr'])).items():
            row[f'{k}_{c}'] = v
    return row, spikes


# ---------------------------------------------------------------------------
# step 0 (no simulation)

def time_mean_input(counts: np.ndarray, M: np.ndarray, scale: float) -> np.ndarray:
    """E129's step-0 reader (`E129_block_survivor_static.time_mean_input`): scale x M @ counts / 0.5 s (Hz)."""
    from experiments import E129_block_survivor_static as E129
    return E129.time_mean_input(counts, M, scale, PULSE_S)


def landing(counts: pd.DataFrame, responses: pd.DataFrame, M: np.ndarray, scale: float) -> pd.DataFrame:
    """Reading 7 for one wiring and condition: per pair of odors (canonical order of the odor names), the directed
    ratios 'ratio_ab' (input under a onto the KCs responding to b and not a, over the input under a onto the KCs
    responding to neither) and 'ratio_ba', their mean 'ratio', and the KC counts and mean inputs beside."""
    from itertools import combinations
    I = time_mean_input(counts.to_numpy(dtype=float), M, scale)
    R = responses.to_numpy(dtype=float) > 0.5
    names = [odor_name(c) for c in counts.columns]
    rows = []
    for i, j in combinations(range(len(names)), 2):
        a, b = sorted(((names[i], i), (names[j], j)))
        (na, ia), (nb, ib) = a, b
        r = dict(odor_a=na, odor_b=nb)
        for key, x, y in (('ab', ia, ib), ('ba', ib, ia)):
            only = R[:, y] & ~R[:, x]
            neither = ~R[:, x] & ~R[:, y]
            num = float(I[only, x].mean()) if only.any() else np.nan
            den = float(I[neither, x].mean()) if neither.any() else np.nan
            r.update({f'n_only_{key}': int(only.sum()), f'input_only_{key}': num, f'input_neither_{key}': den,
                      f'ratio_{key}': num / den if den > 0 else np.nan})
        r['n_neither'] = int((~R[:, ia] & ~R[:, ib]).sum())
        r['ratio'] = float(np.nanmean([r['ratio_ab'], r['ratio_ba']]))
        rows.append(r)
    return pd.DataFrame(rows)


def pair_groups(pairs: pd.DataFrame) -> Dict[str, np.ndarray]:
    """E90's groups of the megamat pairs (ester, alcohol, ester_alcohol, other) and 'all120', over a table's
    (odor_a, odor_b) rows."""
    from experiments import E90_model_stages as E90
    idx = pd.MultiIndex.from_arrays([pairs['odor_a'], pairs['odor_b']], names=['odor_a', 'odor_b'])
    m = E90.group_masks(idx)
    out = {g: np.asarray(m[g]) for g in ('ester', 'alcohol', 'ester_alcohol', 'other')}
    out['all120'] = np.ones(len(pairs), dtype=bool)
    return out


LANDING_GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'all120')


def landing_groups(pairs: pd.DataFrame) -> Dict[str, float]:
    """The mean pair ratio over each group (`pair_groups`)."""
    g = pair_groups(pairs)
    return {k: float(np.nanmean(pairs.loc[g[k], 'ratio'])) for k in LANDING_GROUPS}


def perron(A: np.ndarray) -> dict:
    """The largest real eigenvalue of a nonnegative matrix (its Perron root), its multiplicity within 1e-9, the
    spectral radius, and the uniform share (cos² with the all-ones vector) of the leading right and left
    eigenvectors."""
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    out = {}
    for side, X in (('right', A), ('left', A.T)):
        w, V = np.linalg.eig(X)
        k = int(np.argmax(w.real))
        v = np.real(V[:, k])
        out[f'uniform_share_{side}'] = float(v.sum() ** 2 / (n * (v ** 2).sum())) if (v ** 2).sum() > 0 else np.nan
        if side == 'right':
            root = float(w.real[k])
            out.update(perron_root=root, n_at_root=int((np.abs(w - root) < 1e-9).sum()),
                       spectral_radius=float(np.abs(w).max()), imag_at_root=float(abs(w.imag[k])))
    return out


def stability(HB: np.ndarray, types: np.ndarray) -> pd.DataFrame:
    """Reading 8: `perron` of HB and of each KC type's diagonal block, with the block's size and mean row sum."""
    HB = np.asarray(HB, dtype=float)
    types = np.asarray(types, dtype=object)
    rows = []
    for t in ('all',) + TYPES:
        m = np.ones(len(types), dtype=bool) if t == 'all' else types == t
        B = HB[np.ix_(m, m)]
        rows.append(dict(block='all' if t == 'all' else TYPE_COL[t], n_kcs=int(m.sum()),
                         row_sum_mean=float(B.sum(axis=1).mean()), **perron(B)))
    return pd.DataFrame(rows)


def ref_runs_table(experiment: str, pattern: str) -> pd.DataFrame:
    """Another experiment's run rows (results/<experiment>/<pattern>), the last per tag."""
    ps = sorted((H.shared_dir() / 'results' / experiment).glob(pattern))
    if not ps:
        return pd.DataFrame(columns=['tag'])
    return pd.concat([pd.read_csv(p) for p in ps], ignore_index=True, sort=False).drop_duplicates(
        'tag', keep='last').set_index('tag')


def run_meta_row(kind: str, w: str, o: Mapping[str, tuple], **meta) -> dict:
    """A run's row: kind, wiring, meta, realism (E83's), the types' response rates (E129's), E110's run extras and
    the four bands (E123's)."""
    from experiments import E110_no_kckc as E110
    from experiments import E123_kckc_form as E123
    from experiments import E129_block_survivor_static as E129
    r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E129.type_rates(o), **E110.run_extras(o))
    r.update(E123.bands(r))
    return r


def references() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """The references' outputs and runs rows (REFS), each on the wirings it has, read from the saved runs."""
    from experiments import E110_no_kckc as E110
    outs, rows = {}, []

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        rows.append(run_meta_row(kind, w, o, **meta))

    def pk(r):
        return {k: float(r[k]) for k in E83.PARAM_KEYS if k in r}

    ra = E110.runs_table_of('A').set_index('tag')
    r131 = ref_runs_table('E131_mbc1_cells', 'runs_*.csv')
    r125 = ref_runs_table('E125_kckc_survivor', 'runs_1_*.csv')
    r122 = ref_runs_table('E122_feedback_neuron', 'runs_1.csv')
    r123 = ref_runs_table('E123_kckc_form', 'runs_1_*.csv')
    r129 = ref_runs_table('E129_block_survivor_static', 'runs_step1a_*.csv')
    for w in WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], corner=False, **F['p'])
        o = _load_from('E131_mbc1_cells', 'step1', f'S1_H50_{w}')
        if o is not None and f'S1_H50_{w}' in r131.index:
            r = r131.loc[f'S1_H50_{w}']
            add('H50', w, o, tag=f'S1_H50_{w}', folder='E131_mbc1_cells/step1', corner=bool(r['corner']), **pk(r))
        a = E110.load_a(w)
        if a is not None and f'A_{w}' in ra.index:
            add('E110_A', w, a, tag=f'A_{w}', folder='E110_no_kckc/A', corner=False, **pk(ra.loc[f'A_{w}']))
        b = E110.counted_b(w)
        if b is not None:
            o = E110.load_fit_outputs(b['name'], b['tag'])
            if o is not None:
                add('E110_B', w, o, tag=b['tag'], folder=f"{E110.FIT_EXPERIMENT_PREFIX}/{b['name']}", corner=False,
                    **E110.params_of_row(b['row']))
        if w in S_WIRINGS:
            u = ('E122_feedback_neuron', 'step1', r122)
        else:
            u = ('E123_kckc_form', 'step1', r123)
        o = _load_from(u[0], u[1], f'S1_U_{w}')
        if o is not None and f'S1_U_{w}' in u[2].index:
            add('U', w, o, tag=f'S1_U_{w}', folder=f'{u[0]}/step1', corner=False, **pk(u[2].loc[f'S1_U_{w}']))
        o = _load_from('E129_block_survivor_static', 'step1a', f'S1_P_{w}')
        if o is not None and f'S1_P_{w}' in r129.index:
            r = r129.loc[f'S1_P_{w}']
            add('E129_P', w, o, tag=f'S1_P_{w}', folder='E129_block_survivor_static/step1a', corner=bool(r['corner']),
                pooled_factor=float(r['factor']), kckc_scale_off=float(r['kckc_scale_off']), **pk(r))
        o = _load_from('E125_kckc_survivor', 'step1', f'S1_PO50_{w}')
        if o is not None and f'S1_PO50_{w}' in r125.index:
            r = r125.loc[f'S1_PO50_{w}']
            add('PO50', w, o, tag=f'S1_PO50_{w}', folder='E125_kckc_survivor/step1',
                corner=bool(r['at_bound_thr'] or r['at_bound_apl']), **pk(r))
    for part in ('C', 'D'):
        for spec in E110.fit_specs(part):
            tag = E110.fit_tag(spec)
            r = E110.fit_row(spec.name, tag)
            o = None if r is None else E110.load_fit_outputs(spec.name, tag)
            if o is not None:
                add(f'E110_{part}', spec.wiring, o, tag=tag, folder=f'{E110.FIT_EXPERIMENT_PREFIX}/{spec.name}',
                    corner=False, **E110.params_of_row(r))
    runs = pd.DataFrame(rows)
    runs['corner'] = runs['corner'].fillna(False).astype(bool)
    return outs, runs


def _dist(o: Mapping[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (E83.pair_distances), cached by the sha256 of its spike counts and the readout."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    path = local_dir() / 'distances_cache.pkl'
    if 'dist' not in _CACHE:
        try:
            _CACHE['dist'] = pickle.loads(path.read_bytes()) if path.exists() else {}
        except Exception:
            _CACHE['dist'] = {}
    key = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if key not in _CACHE['dist']:
        _CACHE['dist'][key] = E83.pair_distances(dict(o), readout)
        _CACHE['dist_dirty'] = True
    return _CACHE['dist'][key]


def _save_dist_cache() -> None:
    import pickle
    if _CACHE.get('dist_dirty'):
        try:
            (local_dir() / 'distances_cache.pkl').write_bytes(pickle.dumps(_CACHE['dist']))
        except OSError:
            pass
        _CACHE['dist_dirty'] = False


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'n_kc_ge20_on', 'n_kc_ge20_off', 'silent_on', 'silent_off', 'fixed_thr', 'wAPLKC_scale',
            'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale', 'wU_scale', 'pooled_factor',
            'late_share_on', 'late_share_off', 'kckc_in_over_thr_on', 'kckc_in_over_thr_off', 'pooled_in_over_thr_on',
            'pooled_in_over_thr_off') + tuple(
    f'{k}_{TYPE_COL[t]}' for t in TYPES for k in ('rate_on', 'rate_off', 'ratio', 'count_on', 'count_off'))


def readings(outs: Mapping[Tuple[str, str], Mapping[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8', 'C4' when all of the set have one) and readout:
    E110's `scores_of`, the runs' means of RUN_COLS ('run_<col>'), the number of runs outside the bands."""
    from experiments import E110_no_kckc as E110
    rows = []
    for kind in dict.fromkeys(runs['kind']):
        have = [w for w in WIRINGS if (kind, w) in outs]
        sub = runs[runs['kind'] == kind].drop_duplicates('wiring', keep='last').set_index('wiring')
        units = [(w, (w,)) for w in have] + [(u, ws) for u, ws in POP_WIRINGS.items() if all(w in have for w in ws)]
        for ro in READOUTS:
            dd = {w: _dist(outs[(kind, w)], ro) for w in have}
            for unit, ws in units:
                d = dd[unit] if len(ws) == 1 else E83.population([dd[w] for w in ws])
                s = sub.loc[list(ws)]
                rows.append(dict(kind=kind, unit=unit, n_wirings=len(ws), readout=ro,
                                 n_outside_bands=int((~s['in_bands'].astype(bool)).sum()),
                                 **{f'run_{c}': float(pd.to_numeric(s[c], errors='coerce').mean())
                                    for c in RUN_COLS if c in s.columns},
                                 **E110.scores_of(d)))
    _save_dist_cache()
    return pd.DataFrame(rows)


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def slopes(rd: pd.DataFrame, scales: Mapping[str, float]) -> pd.DataFrame:
    """Reading 9: per unit (each wiring, S4, T4, all8) and reading (SLOPE_READINGS), the family's ('at_s'), H50's
    ('at_half') and E110 part A's ('at_zero') values, the change per unit of scale on the segments s -> s / 2
    ('slope_hi') and s / 2 -> 0 ('slope_lo') and over s -> 0 ('slope_all'), and carried through zero at doses -2, -5,
    -10 by the half->zero segment ('at_m<d>') and by the whole one ('at_m<d>_all'). `scales`: each wiring's s; a
    population's is the mean of its wirings'."""
    rows = []
    units = list(WIRINGS) + ['S4', 'T4', 'all8']
    for u in units:
        ws = POP_WIRINGS.get(u, (u,))
        s = float(np.mean([scales[w] for w in ws]))
        for name, col, ro in SLOPE_READINGS:
            f, h, z = (_v(rd, k, u, col, ro) for k in ('family', 'H50', 'E110_A'))
            r = dict(unit=u, reading=name, column=col, readout=ro, s=s, at_s=f, at_half=h, at_zero=z,
                     slope_hi=(h - f) / (s / 2 - s), slope_lo=(z - h) / (0 - s / 2), slope_all=(z - f) / (0 - s))
            for d in LANDING_DOSES:
                r[f'at_m{d:g}'] = z - r['slope_lo'] * d
                r[f'at_m{d:g}_all'] = z - r['slope_all'] * d
            rows.append(r)
    return pd.DataFrame(rows)


def part_step0() -> dict:
    """Step 0 (no simulation): step0_references.csv, step0_ref_runs.csv, step0_slopes.csv, step0_landing_pairs.csv,
    step0_landing.csv, step0_stability.csv, step0.json, provenance_step0.json."""
    from experiments import E129_block_survivor_static as E129
    out = out_dir()
    outs, runs = references()
    rd = readings(outs, runs)
    rd.to_csv(out / 'step0_references.csv', index=False)
    runs.to_csv(out / 'step0_ref_runs.csv', index=False)
    scales = {w: float(E83.saved_fit(w)['p']['wKCKC_scale']) for w in WIRINGS}
    sl = slopes(rd, scales)
    sl.to_csv(out / 'step0_slopes.csv', index=False)
    idx = E83.saved_fit('real')['saved']['on'][1].index
    for w in WIRINGS:
        if not E83.saved_fit(w)['saved']['on'][1].index.equals(idx):
            raise RuntimeError(f"{w}'s KCs are not the hemibrain wiring's, in its order")
    HB = E129.hemibrain_matrix(idx).to_numpy(dtype=float)
    pair_rows, grp_rows, inputs = [], [], {}
    for w in WIRINGS:
        F = E83.saved_fit(w)
        s = scales[w]
        thr = float(F['p']['fixed_thr'])
        for c in CONDS:
            r, n = F['saved'][c]
            lp = landing(n, r, HB, s)
            grp_rows.append(dict(wiring=w, condition=c, unit=w, **landing_groups(lp)))
            lp.insert(0, 'condition', c)
            lp.insert(0, 'wiring', w)
            pair_rows.append(lp)
            I = time_mean_input(n.to_numpy(dtype=float), HB, s)
            inputs[f'{w}|{c}'] = dict(mean_hz=float(I.mean()), over_thr=float(I.mean() / thr), fixed_thr=thr, s=s)
    pairs = pd.concat(pair_rows, ignore_index=True)
    pairs.to_csv(out / 'step0_landing_pairs.csv', index=False)
    grp = pd.DataFrame(grp_rows)
    for u in ('S4', 'T4', 'all8'):
        for c in CONDS:
            sub = grp[(grp['condition'] == c) & grp['wiring'].isin(POP_WIRINGS[u])]
            grp_rows.append(dict(wiring='', condition=c, unit=u, **{g: float(sub[g].mean()) for g in LANDING_GROUPS}))
    grp = pd.DataFrame(grp_rows)
    grp['alcohol_above_ester'] = grp['alcohol'] > grp['ester']
    grp.to_csv(out / 'step0_landing.csv', index=False)
    st = stability(HB, kc_types_of(idx))
    root = float(st.set_index('block').loc['all', 'perron_root'])
    for d in LANDING_DOSES:
        st[f'dose{d:g}_x_root'] = d * st['perron_root']
    st.to_csv(out / 'step0_stability.csv', index=False)
    refs = {}
    for k in REFS:
        refs[k] = {}
        for u in ('S4', 'T4', 'all8', 'C4'):
            if not ((rd['kind'] == k) & (rd['unit'] == u)).any():
                continue
            refs[k][u] = dict(pattern=_v(rd, k, u, 'all_pattern'), amplitude=_v(rd, k, u, 'all_amplitude'),
                              shape=_v(rd, k, u, 'all_shape'), sd_ratio=_v(rd, k, u, 'all_sd_ratio'),
                              mean_effect_cal=_v(rd, k, u, 'mean_effect'),
                              mean_effect_raw=_v(rd, k, u, 'mean_effect', 'raw'), loss_raw=_v(rd, k, u, 'all_loss', 'raw'),
                              max_count_off=_v(rd, k, u, 'run_max_count_off'),
                              n_kc_ge20_off=_v(rd, k, u, 'run_n_kc_ge20_off'),
                              wKCKC_scale=_v(rd, k, u, 'run_wKCKC_scale'), fixed_thr=_v(rd, k, u, 'run_fixed_thr'))
    lg = grp[grp['unit'].isin(['S4', 'T4', 'all8'])]
    v = dict(generated=_now(), host=_host(), route=ROUTE, references=refs,
             ref_wirings={k: sorted(runs.loc[runs['kind'] == k, 'wiring']) for k in REFS},
             slope={u: {r['reading']: dict(slope_lo=r['slope_lo'], slope_all=r['slope_all'],
                                           **{f'at_m{d:g}': r[f'at_m{d:g}'] for d in LANDING_DOSES})
                        for r in sl[sl['unit'] == u].to_dict('records')} for u in ('S4', 'T4', 'all8')},
             landing={f"{r['unit']}|{r['condition']}": {g: r[g] for g in LANDING_GROUPS} for r in lg.to_dict('records')},
             landing_alcohol_above_ester={f"{r['unit']}|{r['condition']}": bool(r['alcohol_above_ester'])
                                          for r in lg.to_dict('records')},
             family_kckc_input=inputs,
             stability=dict(perron_root=root, blocks=st.to_dict('records'),
                            s_x_root={w: scales[w] * root for w in WIRINGS},
                            dose_x_root={f'{d:g}': d * root for d in LANDING_DOSES}),
             scales=scales, hb_sha256=sha256_array(HB))
    _write_json(out / 'step0.json', v)
    _write_json(out / 'provenance_step0.json', dict(
        generated=_now(), host=_host(), script=SCRIPT, part='step0', commits=H.code_commits(),
        inputs=dict(family={w: list(E83.members()[w]) for w in WIRINGS},
                    references={k: sorted(set(runs.loc[runs['kind'] == k, 'folder'])) for k in REFS},
                    kckc='pyolfsysm.mb_model.connectome_KCKC_weights(model KC IDs) (E129.hemibrain_matrix)')))
    print(json.dumps(dict(references={k: {u: round(x['pattern'], 4) for u, x in d.items()} for k, d in refs.items()},
                          slope_S4=v['slope'].get('S4', {}).get('pattern'), landing=v['landing'],
                          perron_root=root), indent=1, default=str), flush=True)
    return v


# ---------------------------------------------------------------------------
# the checks and late firing (simulations)

def e129_factor(w: str) -> float:
    """E129's bisected factor of its pooled arm P on a wiring, from its provenance (JSON, exact)."""
    vals = set()
    for p in sorted((H.shared_dir() / 'results' / 'E129_block_survivor_static').glob('provenance_step1a_*.json')):
        pt = _read_json(p).get('points', {}).get(f'S1_P_{w}')
        if pt is not None:
            vals.add(float(pt['held']['factor']))
    if len(vals) != 1:
        raise RuntimeError(f"E129's S1_P_{w} factor: {sorted(vals)} in its provenance")
    return vals.pop()


def route_diagnostic(w: str, b: float, target: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    """The feedback-cell route (olfsysm fa3543d): one uniform cell (input weights 1 / N, output weights HB's total / N,
    so that output x input is U) at scale b with the matrix 0, APL off, against `target` (E129's run). Not a gate."""
    from experiments import E129_block_survivor_static as E129
    P = E83.wiring_parts(w)
    idx = P['saved']['on'][1].index
    ids = pd.Index(idx.get_level_values(H.KC_ID_LEVEL))
    HB = E129.hemibrain_matrix(idx).to_numpy(dtype=float)
    N = len(ids)
    w_in = pd.DataFrame(np.full((N, 1), 1.0 / N), index=ids, columns=['U'])
    w_out = pd.DataFrame(np.full((N, 1), float(HB.sum()) / N), index=ids, columns=['U'])
    p = {k: float(P['p'][k]) for k in E83.PARAM_KEYS}
    sim, _ = E83.build_simulator(w, params=dict(p, wMBC1_scale=0.0),
                                 extra_kws=dict(mbc1_inh=0.0, _wMBC1KC=w_out, _wKCMBC1=w_in))
    out = sim(_off(dict(p, wKCKC_scale=0.0, wMBC1_scale=float(b))))
    d = _ndiff(out, target)
    del sim
    gc.collect()
    return dict(route='feedback cell (mbc1_inh, one uniform cell)', b=float(b), n_differing=d['counts'],
                n_responses_differing=d['responses'], n_responses=n_resp(out[0]), target_n=n_resp(target[0]))


def part_check(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, diagnostic: bool = True) -> dict:
    """The checks (the module docstring): check_<host>_<lane>.json."""
    from experiments import E108_paa_extrapolation as E108
    _preload()
    lane = lane or lane_name(wirings)
    path = out_dir() / f'check_{_host()}_{lane}.json'
    res = dict(_read_json(path), host=_host(), lane=lane, script=SCRIPT, route=ROUTE, commits=E108._commits_clean())
    res.setdefault('wirings', {})
    for w in wirings:
        r = res['wirings'].get(w, {})
        if r.get('passes') and (w != 'real' or not diagnostic or 'route_diagnostic' in r):
            print(f'E132 check {w}: done', flush=True)
            continue
        t0 = time.time()
        R = get_runner(w)
        t1 = time.time()
        fam = {}
        for c, q in (('on', R.params(R.s, 0.0)), ('off', _off(R.params(R.s, 0.0)))):
            ts = time.time()
            fam[c] = dict(_ndiff(R(q), R.saved[c]), seconds=time.time() - ts)
        f = e129_factor(w)
        b = R.s * f
        e129 = _load_from('E129_block_survivor_static', 'step1a', f'S1_P_{w}')
        if e129 is None:
            raise RuntimeError(f"E129's S1_P_{w} run is neither local nor in the archive")
        ts = time.time()
        pooled = R(_off(R.params(0.0, b)))
        pooled_d = dict(_ndiff(pooled, e129['off']), seconds=time.time() - ts)
        r = dict(family=fam, e129_pooled=dict(pooled_d, factor=f, b=b, n_responses=n_resp(pooled[0]),
                                              target_off=R.target_off),
                 hb_sha256=R.hb_sha256, u_sha256=R.u_sha256, build_seconds=t1 - t0)
        r['passes'] = bool(all(fam[c][k] == 0 for c in CONDS for k in ('counts', 'responses'))
                           and pooled_d['counts'] == 0 and pooled_d['responses'] == 0)
        if w == 'real' and diagnostic:
            _CACHE.pop('runner', None)
            R.close()
            del R
            gc.collect()
            try:
                r['route_diagnostic'] = route_diagnostic(w, b, e129['off'])
            except Exception as e:      # a diagnostic, not a gate
                r['route_diagnostic'] = dict(error=f'{type(e).__name__}: {e}')
        r.update(seconds=time.time() - t0, checked=_now())
        res['wirings'][w] = r
        _write_json(path, res)
        print(f'E132 check {w}: {json.dumps(r, default=str)}', flush=True)
    res['passes'] = all(res['wirings'].get(w, {}).get('passes') for w in wirings)
    _write_json(path, res)
    return res


def check_wirings() -> set:
    """The wirings that passed checks (i) and (ii) in any host's check_*.json."""
    ok = set()
    for p in out_dir().glob('check_*.json'):
        for w, r in _read_json(p).get('wirings', {}).items():
            if r.get('passes'):
                ok.add(w)
    return ok


def check_passed() -> bool:
    return check_wirings() >= set(WIRINGS)


def part_late(wirings: Sequence[str] = WIRINGS) -> pd.DataFrame:
    """The family's late firing (reading 10) from E121's step-0 rebuild (`E121_onset_race.run_wiring`, which checks that
    it reproduces the archive): late_family.csv, its spike times under outputs/E132_kckc_excitation/late/."""
    from experiments import E121_onset_race as E121
    out = out_dir() / 'late_family.csv'
    have = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=['wiring'])
    rows = [r for r in have.to_dict('records') if r['wiring'] not in wirings]
    root = H.shared_dir() / 'data' / 'model_outputs'
    for w in [x for x in WIRINGS if x in wirings]:
        lab = 'S' if w == 'real' else w
        exp, tag = E121.WIRINGS[lab]
        t0 = time.time()
        res = E121.run_wiring(lab, root / exp / 'rebuild' / tag, local_dir('late'))
        F = E83.saved_fit(w)
        for c in CONDS:
            spk = res[f'spikes_{c}']
            same = bool(np.array_equal(res[f'counts_{c}'].to_numpy(dtype=float), F['saved'][c][1].to_numpy(dtype=float)))
            if not same:
                raise RuntimeError(f"{w} {c}: E121's rebuild is not the family's saved run")
            spk.to_parquet(local_dir('late') / f'{w}_{c}_spike_times.parquet')
            rows.append(dict(wiring=w, condition=c, n_spikes=int(len(spk)), late_share=late_share(spk),
                             n_late=int((spk['time_s'].to_numpy(dtype=float) >= STIM_END - 1e-12).sum()),
                             same_as_saved=same, rebuild=f'{exp}/rebuild/{tag}', seconds=time.time() - t0))
        del res
        gc.collect()
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f'E132 late {w}: {[(r["condition"], round(r["late_share"], 4)) for r in rows if r["wiring"] == w]}',
              flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    _write_json(out_dir() / 'provenance_late.json', dict(generated=_now(), host=_host(), script=SCRIPT, part='late',
                                                         commits=H.code_commits(),
                                                         source='experiments/E121_onset_race.py run_wiring'))
    return df


# ---------------------------------------------------------------------------
# the runs (simulations)

def planned(wirings: Sequence[str] = WIRINGS, arms: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """The (arm, wiring) points of a lane: wiring by wiring (one runner), step 1's arms in STEP1_ARMS order."""
    arms = STEP1_ARMS if arms is None else tuple(arms)
    bad = [a for a in arms if a not in STEP1_ARMS]
    if bad:
        raise SignedDoseError(f"{bad}: --part runs runs step 1's arms; step 2's run with --part step2")
    return [(a, w) for w in WIRINGS if w in wirings for a in STEP1_ARMS if a in arms]


def planned_step2(wirings: Sequence[str] = WIRINGS, arms: Optional[Sequence[str]] = None,
                  reverse: bool = False) -> List[Tuple[str, str]]:
    """Step 2's (arm, wiring) points of a lane (reading 14): the smoke (UX2, then X2, on 'real') first when the lane has
    'real', then the UX arms wiring by wiring in the order `wirings` gives (UX2, UX5, UX10 on each), then the X arms
    likewise; `reverse` reverses the whole list."""
    arms = STEP2_ARMS if arms is None else tuple(arms)
    bad = [a for a in arms if a not in STEP2_ARMS]
    if bad:
        raise ValueError(f"{bad}: step 2's arms are {STEP2_ARMS}")
    ws = list(dict.fromkeys(wirings))
    unknown = [w for w in ws if w not in WIRINGS]
    if unknown:
        raise ValueError(f'wirings {unknown}')
    pts = [p for p in SMOKE if p[0] in arms and p[1] in ws]
    for group in (('UX2', 'UX5', 'UX10'), ('X2', 'X5', 'X10')):
        pts += [(a, w) for w in ws for a in group if a in arms and (a, w) not in pts]
    return pts[::-1] if reverse else pts


def step2_lane(wirings: Sequence[str], reverse: bool = False) -> str:
    return 'step2-' + lane_name(wirings) + ('-rev' if reverse else '')


def runs_table() -> pd.DataFrame:
    """Every row of every host and lane (runs_*.csv), the last row per tag."""
    fs = [pd.read_csv(p) for p in sorted(out_dir().glob('runs_*.csv')) if p.name != 'runs.csv']
    if not fs:
        return pd.DataFrame(columns=['tag', 'arm', 'wiring'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def claim(tag: str, lane: str, reclaim: bool = False) -> bool:
    """Claims a point by mkdir of results/<EXPERIMENT>/claims/<tag>/ (atomic on NFS); False if another process holds
    it (unless `reclaim`)."""
    d = out_dir() / 'claims' / tag
    try:
        d.mkdir(parents=True)
    except FileExistsError:
        if not reclaim:
            return False
    _write_json(d / 'claim.json', dict(tag=tag, host=_host(), lane=lane, pid=os.getpid(), time=_now()))
    return True


def release(tag: str) -> None:
    shutil.rmtree(out_dir() / 'claims' / tag, ignore_errors=True)


def save_run(tag: str, o: Mapping[str, tuple], spikes: Optional[Mapping[str, pd.DataFrame]] = None,
             part: str = 'step1') -> None:
    """The run's conditions (and spike times), locally and in the model-output archive."""
    for d in (local_dir(part), archive_dir() / part):
        d.mkdir(parents=True, exist_ok=True)
        items = [(f'{tag}_{cc}_{k}.parquet', x) for c, cc in (('on', 'on'), ('off', 'acute_off'))
                 for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1]))]
        items += [(f'{tag}_{c}_spike_times.parquet', x) for c, x in (spikes or {}).items()]
        for name, x in items:
            p = d / name
            tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
            x.to_parquet(tmp)
            os.replace(tmp, p)


def load_run(tag: str, part: str = 'step1') -> Optional[Dict[str, tuple]]:
    for d in (H.HERE / 'outputs' / EXPERIMENT / part, archive_dir() / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_row(tag: str, arm: str, w: str, res: Mapping, o: Mapping[str, tuple], real: Mapping, cc: dict, lane: str,
            runner: Runner) -> dict:
    a = ARMS[arm]
    po = res['params_off']
    row = dict(tag=tag, arm=arm, step=a['step'], hold=a['hold'], dose=a['dose'], pooled=a['pooled'], wiring=w,
               allow_negative_kc_kc=runner.allow_negative, fit_kckc_scale=runner.s,
               fit_folder=E83.members()[w][0], fit_tag=E83.members()[w][1], pooled_factor=res['pooled_factor'],
               pooled_scale=res['pooled_scale'], thr_factor=res['thr_factor'], apl_factor=res['apl_factor'],
               n_on=res['n_on'], n_off=res['n_off'], target_on=res['target_on'], target_off=res['target_off'],
               at_bound_thr=bool(res['at_bound_thr']), at_bound_apl=bool(res['at_bound_apl']), held_on=res['held_on'],
               held_off=res['held_off'], corner=res['corner'], first_stage=res['first_stage'], n_sims=res['n_sims'],
               seconds=res.get('seconds'), **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in po.items()}), route=ROUTE, hb_sha256=runner.hb_sha256,
               u_sha256=runner.u_sha256, host=_host(), lane=lane)
    row.update(real)
    r = run_meta_row(arm, w, o)
    r.pop('kind', None)
    r.pop('wiring', None)
    row.update(r)
    row['thr_below_60'] = bool(float(row['fixed_thr']) < REALISM_THR_HZ)
    row.update(cc)
    return row


def hold_change(arm: str) -> str:
    """The provenance's description of an arm's hold."""
    a = ARMS[arm]
    if a['hold'] == 'pooled':
        return (f"arm {arm} ({ARM_LABELS[arm]}): the KC-KC matrix d HB + b U, d {a['dose']:g}; the pooled hold (b "
                "bisected from the fit's KC-KC scale with APL off, then the tied APL scales with APL on; threshold the "
                "fit's)")
    return (f"arm {arm} ({ARM_LABELS[arm]}): the KC-KC matrix d HB, d {a['dose']:g}, no pooled neuron; the threshold "
            "hold (the threshold bisected from the fit's with APL off, then the tied APL scales with APL on; E110's "
            "part A, E77_input_screen.hold_rates' brackets)")


def run_point(arm: str, w: str, lane: str, table: Path, prov_path: Path, cc: dict, reclaim: bool = False,
              part_name: str = 'runs') -> Optional[dict]:
    """One (arm, wiring) point of step 1 or 2, claimed (None if done or claimed elsewhere): the arm's hold on the
    wiring's runner (with the option's flag for step 2), simulated again with spike times, saved (step<n>/), its row
    appended to `table` and its provenance to `prov_path`. Returns the row."""
    from experiments import E108_paa_extrapolation as E108
    a = ARMS[arm]
    tag = tag_of(arm, w)
    if tag in set(runs_table()['tag']):
        print(f'E132 {tag}: done', flush=True)
        return None
    if not claim(tag, lane, reclaim):
        print(f'E132 {tag}: claimed by another process, skipped', flush=True)
        return None
    try:
        E108._still_clean(cc)
        P = E83.wiring_parts(w)
        runner = get_runner(w, allow_negative=a['step'] >= 2)
        d = float(a['dose'])
        hold = pooled_hold if a['hold'] == 'pooled' else threshold_hold
        res = hold(runner, d, f'E132 {tag}')
        o = {'on': res['on'], 'off': res['off']}
        extra = {}
        if arm == 'U0':
            e129 = _load_from('E129_block_survivor_static', 'step1a', f'S1_P_{w}')
            dd = _ndiff(res['off'], e129['off'])
            if dd['counts'] or dd['responses'] or res['pooled_factor'] != e129_factor(w):
                raise RuntimeError(f"{tag}: the APL-off run is not E129's S1_P_{w} ({dd}, factor "
                                   f"{res['pooled_factor']!r} against {e129_factor(w)!r})")
            extra['off_is_e129'] = True
        real, spikes = run_realism(runner, res, d)
        real.update(extra)
        if a['step'] >= 2:
            real.update(olfsysm_build())
        save_run(tag, o, spikes, part=f"step{a['step']}")
        row = run_row(tag, arm, w, res, o, real, cc, lane, runner)
        E83._append_row(table, row)
        prov = _read_json(prov_path) or dict(part=part_name, lane=lane, script=SCRIPT, points={})
        held = dict(change=hold_change(arm), route=ROUTE, fit=f"{P['folder']}:{P['tag']}", dose=d,
                    pooled_scale=res['pooled_scale'], pooled_factor=res['pooled_factor'], hb_sha256=runner.hb_sha256,
                    u_sha256=runner.u_sha256, allow_negative_kc_kc=runner.allow_negative, history=res['history'])
        if a['hold'] == 'pooled':
            held.update(bracket_pooled=list(SCALE_BRACKET), step_pooled=SCALE_STEP)
        prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=held)
        _write_json(prov_path, prov)
        late = '/'.join('nan' if not np.isfinite(row[f'late_share_{c}']) else f"{row[f'late_share_{c}']:.3f}"
                        for c in CONDS)
        print(f"E132 {tag}: b/s {res['pooled_factor']:.4f}, threshold {row['fixed_thr']:.1f} ({a['hold']} hold, "
              f"first stage x{res['thr_factor']:.4f}), APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} "
              f"on, {res['n_off']}/"
              f"{res['target_off']} off); corner {res['corner']}; bands {row['in_bands']}; max off "
              f"{row['max_count_off']}; KCs >= 20 off {row['n_kc_ge20_off']}; late on/off {late}; "
              f"{res.get('seconds', 0):.0f}s", flush=True)
        del res, o, spikes
        gc.collect()
        return row
    except BaseException:
        release(tag)
        raise


def _close_runner() -> None:
    c = _CACHE.pop('runner', None)
    if c is not None:
        c.close()


def part_runs(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
              arms: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 1's arms (the module docstring), resumable, one point claimed at a time."""
    from experiments import E108_paa_extrapolation as E108
    pts = planned(wirings, arms)
    if not check_passed():
        raise SystemExit('E132 runs refused: checks (i) and (ii) have not passed on every wiring (run --part check)')
    _preload()
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in pts:
        run_point(arm, w, lane, table, prov_path, cc, reclaim, 'runs')
    _close_runner()
    return runs_table()


def flag_check_passed() -> bool:
    """Reading 13: a flagcheck_*.json (any host) that passed on 'real' with a build that has the option and the
    installed build's C++ sources (olfsysm's source_sha256)."""
    sha = olfsysm_build()['olfsysm_source_sha256']
    for p in out_dir().glob('flagcheck_*.json'):
        r = _read_json(p)
        b = r.get('olfsysm', {})
        if (r.get('wirings', {}).get('real', {}).get('passes') and b.get('negative_kckc_option')
                and b.get('olfsysm_source_sha256') == sha):
            return True
    return False


def part_flag_check(wirings: Sequence[str] = ('real',)) -> dict:
    """Reading 13: on each wiring, a runner with the option's flag on reproduces the family's saved run (matrix s, pooled
    0) and step 1's FLAG_CHECK_ARM point at its held parameters (its runs row's params_on / params_off), APL on and off,
    with 0 differing spike counts and responses: flagcheck_<host>.json."""
    from experiments import E108_paa_extrapolation as E108
    if not negative_kckc_available():
        raise SignedDoseError('the installed olfsysm has no option for negative KC-KC weights (olfsysm 45e5cf9)')
    _preload()
    path = out_dir() / f'flagcheck_{_host()}.json'
    res = dict(_read_json(path), host=_host(), script=SCRIPT, route=ROUTE, arm=FLAG_CHECK_ARM,
               commits=E108._commits_clean(), olfsysm=olfsysm_build())
    res.setdefault('wirings', {})
    rt = runs_table().drop_duplicates('tag', keep='last').set_index('tag')
    for w in wirings:
        t0 = time.time()
        R = get_runner(w, allow_negative=True)
        fam = {c: _ndiff(R(q), R.saved[c]) for c, q in (('on', R.params(R.s, 0.0)), ('off', _off(R.params(R.s, 0.0))))}
        tag = tag_of(FLAG_CHECK_ARM, w)
        saved = load_run(tag, 'step1')
        if tag not in rt.index or saved is None:
            raise RuntimeError(f"{tag}: step 1's row or run is missing (local outputs, archive)")
        q = {'on': json.loads(rt.loc[tag, 'params_on']), 'off': json.loads(rt.loc[tag, 'params_off'])}
        step1 = {c: _ndiff(R(q[c]), saved[c]) for c in CONDS}
        r = dict(family=fam, step1=dict(tag=tag, params=q, **step1), flag=bool(R.mp.kc.allow_negative_kc_kc),
                 seconds=time.time() - t0, checked=_now())
        r['passes'] = bool(r['flag'] and all(x[k] == 0 for x in (*fam.values(), *step1.values())
                                             for k in ('counts', 'responses')))
        res['wirings'][w] = r
        _write_json(path, res)
        print(f'E132 flag check {w}: {json.dumps(dict(family=fam, step1=step1, passes=r["passes"]))}', flush=True)
    _close_runner()
    res['passes'] = all(res['wirings'].get(w, {}).get('passes') for w in wirings)
    _write_json(path, res)
    return res


def part_step2(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
               arms: Optional[Sequence[str]] = None, reverse: bool = False) -> pd.DataFrame:
    """Step 2 (reading 14): UX2, UX5, UX10 with the pooled neuron (pooled hold) and X2, X5, X10 alone (threshold
    hold), resumable, one point claimed at a time; needs the option, step 1's checks and the flag check."""
    from experiments import E108_paa_extrapolation as E108
    pts = planned_step2(wirings, arms, reverse)
    if not negative_kckc_available():
        raise SignedDoseError('step 2 needs the olfsysm option for negative KC-KC weights (olfsysm 45e5cf9); pull '
                              'olfsysm and run ./venv/bin/pip install -e ../olfsysm')
    if not check_passed():
        raise SystemExit('E132 step 2 refused: checks (i) and (ii) have not passed on every wiring (run --part check)')
    if not flag_check_passed():
        raise SystemExit('E132 step 2 refused: no flag check has passed with this build (run --part flagcheck)')
    _preload()
    lane = lane or step2_lane(wirings, reverse)
    table = out_dir() / f'runs_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    print(f'E132 step 2 lane {lane}: {[tag_of(a, w) for a, w in pts]}', flush=True)
    for arm, w in pts:
        run_point(arm, w, lane, table, prov_path, cc, reclaim, 'step2')
    _close_runner()
    return runs_table()


def step2_outcome() -> dict:
    """verdict.json's 'step2' (empty before step 2 is scored)."""
    v = _read_json(out_dir() / 'verdict.json').get('step2')
    return v if isinstance(v, dict) else {}


def part_step3(*args, **kwargs):
    """Step 3, only for a dose that passes the stop rule (the entry's design item 6): the matrix shuffled within each
    pair of KC types (E122's shuffle, seed 0) at that dose and hold on the eight wirings, and the male CNS's two
    calyces beside, each at its own E130 fit with its own claw-to-claw matrix, at the dose and at its dose 0. Its code
    is written once a dose passes; until then this says why it does not run."""
    s = step2_outcome()
    if not s.get('complete'):
        raise GateClosedError('step 3 waits on step 2 (--part step2 on all eight wirings, then --part score)')
    if s.get('goes_on') is None:
        raise GateClosedError('no dose of step 2 passes the stop rule: the experiment stops after step 2 (no step 3)')
    raise GateClosedError(f"{s['goes_on']} passes the stop rule: step 3's code (the shuffle control and the male "
                          'calyces) is to be written')


def part_step4(*args, **kwargs):
    """Step 4, only if step 3's control keeps less than half of the passing dose's gain over its dose 0 (design item
    7): four refits on S with a signed KC-KC scale and the pooled scale free (E10_joint_fit.py needs the signed scale),
    read on T, then the sign control. Written only if that gate opens."""
    s = step2_outcome()
    if s.get('goes_on') is None:
        raise GateClosedError('step 4 runs only after a dose passes the stop rule and step 3 (none has)')
    raise GateClosedError('step 4 runs only if step 3\'s control keeps less than half of the dose\'s gain; its code '
                          '(a signed KC-KC scale in E10_joint_fit.py) is written then')


# ---------------------------------------------------------------------------
# score

def flies_mean_effect() -> float:
    """The flies' mean APL effect against all recordings (E108's 'all' target)."""
    from experiments import E108_paa_extrapolation as E108
    return float(E108.effect_targets()['all'].mean())


def targets() -> Tuple[str, ...]:
    from experiments import E108_paa_extrapolation as E108
    return ('all',) + tuple(E108.HALVES)


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run: the references (`references`) and step 1's and step 2's runs."""
    outs, runs = references()
    rows = runs.to_dict('records')
    for r in runs_table().to_dict('records'):
        o = load_run(r['tag'], f"step{int(r['step'])}")
        if o is None:
            continue
        outs[(r['arm'], r['wiring'])] = o
        rows.append(dict(r, kind=r['arm']))
    out = pd.DataFrame(rows)
    out['corner'] = out['corner'].fillna(False).astype(bool)
    out['thr_below_60'] = pd.to_numeric(out['fixed_thr'], errors='coerce') < REALISM_THR_HZ
    return outs, out


def complete(runs: pd.DataFrame, kind: str) -> bool:
    return all(((runs['kind'] == kind) & (runs['wiring'] == w)).any() for w in WIRINGS)


def corners(runs: pd.DataFrame, kind: str) -> dict:
    """Reading 3: the wirings on which a kind's run is a corner, and whether there are more than CORNER_WIRINGS_MAX."""
    sub = runs[runs['kind'] == kind]
    c = sub['corner'].astype(bool) if 'corner' in sub else pd.Series(False, index=sub.index)
    ws = sorted(sub.loc[c, 'wiring'])
    return dict(corner_wirings=ws, corner=len(ws) > CORNER_WIRINGS_MAX)


def below_family(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, flies_mean: float) -> dict:
    """Reading 4's "below the family" (E131's)."""
    pat, per_w = True, {}
    for t in targets():
        col = f'{t}_pattern'
        pat &= all(_v(rd, arm, u, col) - _v(rd, 'family', u, col) <= -WITHIN for u in POPS)
        d = {w: _v(rd, arm, w, col) - _v(rd, 'family', w, col) for w in WIRINGS}
        per_w[t] = int(sum(v <= -WITHIN for v in d.values()))
        pat &= per_w[t] >= N_WIRINGS_NEEDED
    mean = {u: _v(rd, arm, u, 'mean_effect') for u in POPS}
    mean_ok = all(abs(v - flies_mean) <= MEAN_TOL for v in mean.values())
    raw = {u: _v(rd, arm, u, 'all_loss', 'raw') - _v(rd, 'family', u, 'all_loss', 'raw') for u in POPS}
    raw_ok = all(v <= RAW_LOSS_TOL for v in raw.values())
    bands = bool(runs.loc[runs['kind'] == arm, 'in_bands'].astype(bool).all())
    return dict(below=bool(pat and mean_ok and raw_ok and bands), pattern_lower=bool(pat), wirings_lower=per_w,
                mean_effect_cal=mean, mean_effect_ok=bool(mean_ok), raw_loss_minus_family=raw, raw_loss_ok=bool(raw_ok),
                all_in_bands=bands)


def classify(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, flies_mean: float) -> dict:
    """Reading 4: the class ('below the family', 'equal', 'costs', 'neither', 'corner (no class)') with E123's flags
    ('keeps_gain', 'equal'), the parts and the realism beside; incomplete without all eight wirings."""
    from experiments import E123_kckc_form as E123
    if not complete(runs, arm):
        return dict(arm=arm, complete=False)
    c = E123.classify(rd, runs, arm)
    b = below_family(rd, runs, arm, flies_mean)
    costs = all(_v(rd, arm, u, 'all_pattern') - _v(rd, 'family', u, 'all_pattern') > WITHIN for u in POPS)
    cls = 'below the family' if b['below'] else 'equal' if c['equal'] else 'costs' if costs else 'neither'
    c.update(below=b, costs=bool(costs), **corners(runs, arm))
    if c['corner']:
        c['cls_before_corner'] = cls
        cls = 'corner (no class)'
    c['cls'] = cls
    c['pattern'] = {u: _v(rd, arm, u, 'all_pattern') for u in POPS + ('all8',)}
    c['d_family'] = {u: _v(rd, arm, u, 'all_pattern') - _v(rd, 'family', u, 'all_pattern') for u in POPS + ('all8',)}
    c['parts'] = {u: {k: _v(rd, arm, u, f'all_{k}') for k in ('amplitude', 'shape', 'sd_ratio', 'pearson_effect',
                                                               'rho_effect')} for u in POPS}
    c['mean_effect'] = {ro: _v(rd, arm, 'all8', 'mean_effect', ro) for ro in READOUTS}
    sub = runs[runs['kind'] == arm]
    c['runs_outside_bands'] = sorted(sub.loc[~sub['in_bands'].astype(bool), 'wiring'])
    c['thr_below_60'] = sorted(sub.loc[sub['thr_below_60'].astype(bool), 'wiring'])
    for col in ('fixed_thr', 'wAPLKC_scale', 'wKCKC_scale', 'wU_scale', 'pooled_factor', 'max_count_off',
                'n_kc_ge20_off', 'late_share_on', 'late_share_off', 'kckc_in_over_thr_off', 'pooled_in_over_thr_off',
                'ratio_g', 'rate_on_g'):
        if col in sub:
            c[f'{col}_mean'] = float(pd.to_numeric(sub[col], errors='coerce').mean())
    return c


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    """Reading 5: E126's "differ" rule with this experiment's corners: 'verdict' '<a> lower', '<b> lower', 'not told
    apart', 'no class (corner)' or 'incomplete'."""
    if not (complete(runs, a) and complete(runs, b)):
        return dict(a=a, b=b, verdict='incomplete')
    per_w = {w: _v(rd, a, w, 'all_pattern') - _v(rd, b, w, 'all_pattern') for w in WIRINGS}
    per_t = {t: {u: _v(rd, a, u, f'{t}_pattern') - _v(rd, b, u, f'{t}_pattern') for u in POPS} for t in targets()}
    n_a, n_b = sum(v < 0 for v in per_w.values()), sum(v > 0 for v in per_w.values())
    a_all = all(v < 0 for t in per_t.values() for v in t.values())
    b_all = all(v > 0 for t in per_t.values() for v in t.values())
    corner = [x for x in (a, b) if corners(runs, x)['corner']]
    if corner:
        verdict = 'no class (corner)'
    elif n_a >= DIFFER_WIRINGS and a_all:
        verdict = f'{a} lower'
    elif n_b >= DIFFER_WIRINGS and b_all:
        verdict = f'{b} lower'
    else:
        verdict = 'not told apart'
    return dict(a=a, b=b, verdict=verdict, n_wirings_a_lower=int(n_a), n_wirings_b_lower=int(n_b),
                a_lower_on_all_targets=bool(a_all), b_lower_on_all_targets=bool(b_all), corner_arms=corner,
                d_per_wiring=per_w, d_per_target=per_t)


def dose0_of(arm: str) -> Optional[str]:
    """A dose's own dose 0: U0 for the pooled arms, KC-KC at zero (E110_A) for the arm alone; None for U0."""
    if arm == 'U0':
        return None
    return 'U0' if ARMS[arm]['pooled'] else 'E110_A'


def stop_rule(rd: pd.DataFrame, runs: pd.DataFrame, arm: str) -> dict:
    """Reading 6 (read on S; the module docstring): conditions (i) to (iv) and 'passes' (None for U0, which is no
    dose; 'candidate' whether a pass would go to step 3, only the negative doses)."""
    if not complete(runs, arm):
        return dict(arm=arm, complete=False)
    fam, a = _v(rd, 'family', 'S4', 'all_pattern'), _v(rd, arm, 'S4', 'all_pattern')
    i = bool(a - fam <= STOP_SLACK)
    mx = {}
    for w in WIRINGS:
        f = runs[(runs['kind'] == 'family') & (runs['wiring'] == w)]['max_count_off']
        x = runs[(runs['kind'] == arm) & (runs['wiring'] == w)]['max_count_off']
        mx[w] = dict(arm=int(x.iloc[-1]), family=int(f.iloc[-1]))
    ii = all(v['arm'] <= MAX_COUNT_FACTOR * v['family'] for v in mx.values())
    sub = runs[runs['kind'] == arm]
    cw = corners(runs, arm)
    iii = bool(sub['in_bands'].astype(bool).all() and len(cw['corner_wirings']) <= CORNER_WIRINGS_MAX)
    z = dose0_of(arm)
    iv, d0 = None, None
    if z is not None and complete(runs, z):
        zs = _v(rd, z, 'S4', 'all_pattern')
        per_w = {w: _v(rd, arm, w, 'all_pattern') - _v(rd, z, w, 'all_pattern') for w in S_WIRINGS}
        n = int(sum(v < 0 for v in per_w.values()))
        iv = bool(a < zs and n >= DOSE0_WIRINGS)
        d0 = dict(kind=z, pattern_S4=zs, d_S4=a - zs, n_s_wirings_lower=n, d_per_wiring=per_w)
    passes = None if z is None else bool(i and ii and iii and bool(iv))
    return dict(arm=arm, complete=True, dose=ARMS[arm]['dose'], pattern_S4=a, family_S4=fam, d_family_S4=a - fam,
                i_pattern=i, ii_max_count=bool(ii), max_count_off=mx,
                max_count_off_vs_entry=bool(all(v['arm'] <= MAX_COUNT_FACTOR * ENTRY_FAMILY_MAX_COUNT_OFF
                                                for v in mx.values())),
                iii_bands_corners=iii, corner_wirings=cw['corner_wirings'], iv_dose0=iv, dose0=d0, passes=passes,
                candidate=bool(passes and ARMS[arm]['dose'] < 0))


KIND_COLS = ('kind', 'dose', 'pattern S4', 'pattern T4', 'd family S4', 'd family T4', 'amplitude S4', 'shape S4',
             'spread ratio S4', 'Pearson r S4', 'ρ effect S4', 'mean effect (readout)', 'mean effect (raw)', 'loss (raw)',
             'b / s', 'threshold', 'thr < 60 Hz', 'max count off', 'KCs ≥ 20 off', 'γ off / on', 'late on / off',
             'corners', 'outside bands', 'class', 'stop rule')


def comparison_md(rd: pd.DataFrame, runs: pd.DataFrame, kinds: Sequence[str], classes: Mapping[str, dict],
                  stops: Mapping[str, dict], comps: Optional[Mapping[str, dict]] = None) -> str:
    from experiments import E123_kckc_form as E123
    lines = ['# E132 comparison tables (generated)', '',
             'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
             '## The kinds (calibrated readout, all recordings)', '',
             '| ' + ' | '.join(KIND_COLS) + ' |',
             '|---|' + '---:|' * (len(KIND_COLS) - 4) + '---|---|---|']
    for k in kinds:
        if not ((rd['kind'] == k) & (rd['unit'] == 'S4')).any():
            continue
        c = classes.get(k, {})
        sub = runs[runs['kind'] == k]

        def mean(col):
            return float(pd.to_numeric(sub[col], errors='coerce').mean()) if col in sub else np.nan
        sr = stops.get(k, {}).get('passes')
        n60 = int(sub['thr_below_60'].astype(bool).sum()) if 'thr_below_60' in sub else 0
        lines.append(
            f"| {k} | {ARMS[k]['dose'] if k in ARMS else ''} | {_v(rd, k, 'S4', 'all_pattern'):.4f} | "
            f"{_v(rd, k, 'T4', 'all_pattern'):.4f} | {_v(rd, k, 'S4', 'all_pattern') - _v(rd, 'family', 'S4', 'all_pattern'):+.4f} | "
            f"{_v(rd, k, 'T4', 'all_pattern') - _v(rd, 'family', 'T4', 'all_pattern'):+.4f} | "
            f"{_v(rd, k, 'S4', 'all_amplitude'):.4f} | {_v(rd, k, 'S4', 'all_shape'):.4f} | "
            f"{_v(rd, k, 'S4', 'all_sd_ratio'):.3f} | {_v(rd, k, 'S4', 'all_pearson_effect'):.3f} | "
            f"{_v(rd, k, 'S4', 'all_rho_effect'):.3f} | {_v(rd, k, 'S4', 'mean_effect'):.3f} | "
            f"{_v(rd, k, 'S4', 'mean_effect', 'raw'):.3f} | {_v(rd, k, 'S4', 'all_loss', 'raw'):.4f} | "
            f"{mean('pooled_factor'):.3f} | {mean('fixed_thr'):.1f} | {n60} | {mean('max_count_off'):.1f} | "
            f"{mean('n_kc_ge20_off'):.0f} | {mean('ratio_g'):.2f} | "
            f"{mean('late_share_on'):.3f} / {mean('late_share_off'):.3f} | "
            f"{len(corners(runs, k)['corner_wirings'])} | {int((~sub['in_bands'].astype(bool)).sum())} | "
            f"{c.get('cls', '')} | {'' if sr is None else sr} |")
    lines.append('')
    st = [(a, s) for a, s in stops.items() if s.get('complete')]
    if st:
        lines += ['## The stop rule (read on S; the negative doses are the candidates)', '',
                  '| arm | dose | pattern S4 | d family S4 | (i) within 0.005 | (ii) max count ≤ 2× | (iii) bands, corners '
                  '| dose 0 | d dose 0 S4 | S wirings lower | (iv) below dose 0 | passes | candidate |',
                  '|---|---:|---:|---:|---|---|---|---|---:|---:|---|---|---|']
        for a, s in st:
            d0 = s.get('dose0') or {}
            lines.append(f"| {a} | {s['dose']:g} | {s['pattern_S4']:.4f} | {s['d_family_S4']:+.4f} | {s['i_pattern']} | "
                         f"{s['ii_max_count']} | {s['iii_bands_corners']} | {d0.get('kind', '')} | "
                         f"{d0.get('d_S4', np.nan):+.4f} | {d0.get('n_s_wirings_lower', '')} | {s['iv_dose0']} | "
                         f"{s['passes']} | {s['candidate']} |")
        lines.append('')
    if comps:
        lines += ['## "Differ" (E126\'s rule: lower on 7 of 8 wirings and on both populations against all five '
                  'targets)', '', '| a | b | verdict | wirings a lower | wirings b lower |', '|---|---|---|---:|---:|']
        for c in comps.values():
            lines.append(f"| {c['a']} | {c['b']} | {c['verdict']} | {c.get('n_wirings_a_lower', '')} | "
                         f"{c.get('n_wirings_b_lower', '')} |")
        lines.append('')
    for u in ('S4', 'T4', 'all8'):
        ks = [k for k in kinds if ((rd['kind'] == k) & (rd['unit'] == u)).any()]
        lines += [f'## {u}', '', E123.comparison_table(rd, ks, u), '']
    lines += ['## Runs and bands', '', E123.runs_bands_table(runs), '']
    return '\n'.join(lines)


def step2_summary(rd: pd.DataFrame, runs: pd.DataFrame, classes: Mapping[str, dict], comps: Mapping[str, dict],
                  stops: Mapping[str, dict]) -> dict:
    """Reading 15's summary of step 2 for verdict.json: the candidates and the dose that goes on ('goes_on', the
    lowest pattern part on S's population among the passing negative doses, once all six arms are complete; None
    otherwise), each arm's class, stop-rule conditions, parts and realism beside, the sign's own effect (UX against its
    mirror and U0), the arm alone against E110's part A and H50, and the pooled neuron's part per dose (UX against X,
    and each arm's corners, bands and largest counts). Descriptions; the only decision is the stop rule's."""
    have = [a for a in STEP2_ARMS if complete(runs, a)]
    done = len(have) == len(STEP2_ARMS)
    cands = [a for a in have if stops.get(a, {}).get('candidate')]
    goes_on = min(cands, key=lambda a: stops[a]['pattern_S4']) if (done and cands) else None
    fam_raw = {u: _v(rd, 'family', u, 'mean_effect', 'raw') for u in POPS + ('all8',)}
    arms = {}
    for a in have:
        z = dose0_of(a)
        s, c = stops[a], classes.get(a, {})
        raw = {u: _v(rd, a, u, 'mean_effect', 'raw') for u in POPS + ('all8',)}
        arms[a] = dict(
            dose=ARMS[a]['dose'], hold=ARMS[a]['hold'], pooled=ARMS[a]['pooled'], dose0=z, cls=c.get('cls'),
            keeps_gain=c.get('keeps_gain'), passes=s.get('passes'), candidate=s.get('candidate'),
            conditions=dict(i_pattern=s.get('i_pattern'), ii_max_count=s.get('ii_max_count'),
                            iii_bands_corners=s.get('iii_bands_corners'), iv_dose0=s.get('iv_dose0')),
            pattern={u: _v(rd, a, u, 'all_pattern') for u in POPS + ('all8',)},
            d_family={u: _v(rd, a, u, 'all_pattern') - _v(rd, 'family', u, 'all_pattern') for u in POPS + ('all8',)},
            d_dose0={u: _v(rd, a, u, 'all_pattern') - _v(rd, z, u, 'all_pattern') for u in POPS + ('all8',)},
            parts=c.get('parts'), mean_effect_cal={u: _v(rd, a, u, 'mean_effect') for u in POPS + ('all8',)},
            mean_effect_raw=raw, mean_effect_raw_minus_family={u: raw[u] - fam_raw[u] for u in raw},
            mean_effect_raw_minus_dose0={u: raw[u] - _v(rd, z, u, 'mean_effect', 'raw') for u in raw},
            loss_raw={u: _v(rd, a, u, 'all_loss', 'raw') for u in POPS + ('all8',)},
            corner_wirings=c.get('corner_wirings'), runs_outside_bands=c.get('runs_outside_bands'),
            thr_below_60=c.get('thr_below_60'), max_count_off=s.get('max_count_off'),
            **{k: c.get(k) for k in ('fixed_thr_mean', 'pooled_factor_mean', 'max_count_off_mean', 'n_kc_ge20_off_mean',
                                     'late_share_on_mean', 'late_share_off_mean', 'ratio_g_mean',
                                     'kckc_in_over_thr_off_mean', 'pooled_in_over_thr_off_mean')})

    def verdict(a, b):
        return comps.get(f'{a}|{b}', {}).get('verdict')
    sign = {a: dict(vs_mirror=verdict(a, MIRROR[a]), mirror=MIRROR[a], vs_U0=verdict(a, 'U0'))
            for a in ('UX2', 'UX5', 'UX10') if a in have}
    alone = {a: dict(vs_E110_A=verdict(a, 'E110_A'), vs_H50=verdict(a, 'H50')) for a in ('X2', 'X5', 'X10') if a in have}
    pooled_part = {}
    for dose, (ux, x) in STEP2_DOSES.items():
        if ux in have and x in have:
            pooled_part[f'-{dose}'] = dict(differ=verdict(ux, x), **{
                k: dict(corner_wirings=arms[k]['corner_wirings'], runs_outside_bands=arms[k]['runs_outside_bands'],
                        ii_max_count=arms[k]['conditions']['ii_max_count'], passes=arms[k]['passes'])
                for k in (ux, x)})
    return dict(complete=done, arms_complete=have, candidates=cands, goes_on=goes_on,
                stops_after_step2=bool(done and not cands), arms=arms, sign_effect=sign, alone=alone,
                pooled_part=pooled_part)


def score() -> dict:
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    kinds = list(dict.fromkeys(runs['kind']))
    mine = [k for k in kinds if k in ARMS]
    fm = flies_mean_effect()
    classes = {a: classify(rd, runs, a, fm) for a in kinds if a != 'family'}
    comps = {f'{a}|{b}': differ(rd, runs, a, b) for a, b in COMPARISONS + STEP2_COMPARISONS
             if a in kinds and b in kinds}
    stops = {a: stop_rule(rd, runs, a) for a in mine}
    late = pd.read_csv(out / 'late_family.csv') if (out / 'late_family.csv').exists() else None
    mirror_lower = [a for a in ('UI2', 'UI5', 'UI10') if comps.get(f'{a}|U0', {}).get('verdict') == f'{a} lower']
    step2 = (step2_summary(rd, runs, classes, comps, stops) if any(a in kinds for a in STEP2_ARMS)
             else 'not run yet')
    v = dict(generated=_now(), host=_host(), route=ROUTE, checks_passed=check_passed(), flies_mean_effect=fm,
             classes=classes, comparisons=comps, stop_rule=stops,
             step1=dict(mirror_lower_than_U0=mirror_lower,
                        pattern_S4_by_dose={a: _v(rd, a, 'S4', 'all_pattern') for a in ('U0',) + STEP1_ARMS[1:]
                                            if a in kinds},
                        corners={a: classes[a].get('corner_wirings', []) for a in mine if classes[a].get('complete')},
                        outside_bands={a: classes[a].get('runs_outside_bands', []) for a in mine
                                       if classes[a].get('complete')}),
             late_family=None if late is None else late.to_dict('records'),
             step2=step2, flag_check_passed=flag_check_passed() if any(a in kinds for a in STEP2_ARMS) else None,
             rules=dict(within=WITHIN, wirings_needed=N_WIRINGS_NEEDED, mean_tol=MEAN_TOL, raw_loss_tol=RAW_LOSS_TOL,
                        corner_wirings_max=CORNER_WIRINGS_MAX, differ_wirings=DIFFER_WIRINGS, rate_tol=RATE_TOL,
                        stop_slack=STOP_SLACK, max_count_factor=MAX_COUNT_FACTOR, dose0_wirings=DOSE0_WIRINGS,
                        targets=list(targets()), populations=list(POPS)),
             n_runs=int(len(runs)))
    _write_json(out / 'verdict.json', v)
    (out / 'comparison_table.md').write_text(comparison_md(rd, runs, kinds, classes, stops, comps))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT,
                                                    commits=H.code_commits()))
    print(json.dumps(dict(classes={a: c.get('cls') for a, c in classes.items()},
                          stop_rule={a: s.get('passes') for a, s in stops.items()}, step1=v['step1'],
                          step2=({k: step2[k] for k in ('complete', 'candidates', 'goes_on', 'stops_after_step2',
                                                         'sign_effect', 'alone')}
                                 if isinstance(step2, dict) else step2)), default=str), flush=True)
    return v


def distances_tag(tag: str) -> str:
    """The tag of a run in model_pair_distances_long.csv: its run tag (S1_U0_real), the name its responses have in the
    model-output archive (`save_run`), so that kcq/apl_effect.py finds them for the realism bands."""
    return tag


def write_variants() -> pd.DataFrame:
    """variants.csv (every E132 run, absolute parameters) and model_pair_distances_long.csv (raw spike counts, keyed by
    the run tag) for the archive."""
    out = out_dir()
    rt = runs_table()
    rows, parts = [], []
    for r in rt.to_dict('records'):
        sub = f"step{int(r['step'])}"
        o = load_run(r['tag'], sub)
        if o is None:
            continue
        arm = r['arm']
        neg = bool(r.get('allow_negative_kc_kc')) if pd.notna(r.get('allow_negative_kc_kc')) else False
        rows.append(dict(
            tag=r['tag'], variant=r['tag'], part=f"S{r['step']}", arm=arm, arm_label=ARM_LABELS[arm], dose=r['dose'],
            hold=r['hold'], base='bouton', k=E83.BREADTH_K, wiring=r['wiring'],
            ref_experiment=E83.members()[r['wiring']][0], ref_fit=E83.members()[r['wiring']][1],
            kind=f"C2's family's {r['wiring']}: {ARM_LABELS[arm]}",
            fixed_thr=r['fixed_thr'], wAPLKC_scale=r['wAPLKC_scale'], wKCAPL_scale=r['wKCAPL_scale'],
            wAPLPN_scale=r['wAPLPN_scale'], wPNAPL_scale=r['wPNAPL_scale'], wKCKC_scale=r['wKCKC_scale'],
            kc_kc_inh=r['wKCKC_scale'], kc_kc_tau=0.05, wU_scale=r['wU_scale'], pooled_factor=r['pooled_factor'],
            kckc_matrix=("wKCKC_scale x hemibrain (connectome_KCKC_weights) + wU_scale x E122's uniform neuron U "
                         "(HB's total / N^2 everywhere), both conditions, 50 ms"
                         + ('; wKCKC_scale < 0 is excitation (allow_negative_kc_kc)' if neg else '')),
            allow_negative_kc_kc=neg, olfsysm_git_hash=r.get('olfsysm_git_hash'),
            multibouton_claws='first', thr_factor=r['thr_factor'], apl_factor=r['apl_factor'], rate_on=r['rate_on'],
            rate_off=r['rate_off'], max_count_on=r['max_count_on'], max_count_off=r['max_count_off'],
            late_share_on=r.get('late_share_on'), late_share_off=r.get('late_share_off'),
            realism_flag=r.get('realism_flag'), corner=bool(r['corner']), subdir=sub))
        key = distances_tag(r['tag'])
        parts.append(pd.concat({(key, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (key, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', 'step0', 'late', 'runs', 'score', 'variants', 'plan',
                                                      'flagcheck', 'step2', 'step3', 'step4'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS),
                    help='the lane\'s wirings (step 2: in this order; flagcheck: default real)')
    ap.add_argument('--arms', nargs='+', default=None, choices=list(ARMS))
    ap.add_argument('--lane', default=None, help='the lane name in the tables (default the wirings joined by -; '
                                                 'step 2: step2-<wirings>[-rev])')
    ap.add_argument('--reclaim', action='store_true', help='take over claims left by a process that died')
    ap.add_argument('--reverse', action='store_true', help="step2 / plan --step 2: run the lane's points in reverse")
    ap.add_argument('--step', type=int, default=1, choices=(1, 2), help='plan: which step\'s points')
    ap.add_argument('--no-diagnostic', action='store_true', help="check: skip the route diagnostic on 'real'")
    return ap


def main() -> None:
    a = build_parser().parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'check':
        part_check(ws, a.lane, diagnostic=not a.no_diagnostic)
    elif a.part == 'step0':
        part_step0()
    elif a.part == 'late':
        part_late(ws)
    elif a.part == 'runs':
        part_runs(ws, a.lane, a.reclaim, a.arms)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'plan':
        pts = planned_step2(ws, a.arms, a.reverse) if a.step == 2 else planned(ws, a.arms)
        for arm, w in pts:
            print(tag_of(arm, w))
    elif a.part == 'flagcheck':
        part_flag_check(list(a.wirings or ('real',)))
    elif a.part == 'step2':
        part_step2(ws, a.lane, a.reclaim, a.arms, a.reverse)
    elif a.part == 'step3':
        part_step3()
    elif a.part == 'step4':
        part_step4()


if __name__ == '__main__':
    main()
