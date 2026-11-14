#!/usr/bin/env python3
"""E96: APL's synapses on the KCs outside the claws, in the calyx and in the peduncle, which the model's weights leave out,
as a component of the APL->KC weights per KC in the feedback APL with the bouton path at a tenth; whether weights that the
connectome has do what E95's uniform component did; and the family fitted again that way.

The entry is item E96 of experiment_priorities.md (proposed by viki 2026-11-14 10:21; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 10:02). Built on E95 (`experiments/E95_static_apl.py`: its part D runner,
readings, choice, gate, equal rule and part F), E94 (`WiringModel94`, the stage readers, D1), E91 (the pooled share),
E83 (the saved fits, the gate), E92 (the bracket) and E69 (`bisect_factor`).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B.
S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's, C / D the second split's (E89's
`targets89`). The input is the lab's in every run.

The counts (pyolfsysm's `connectome_APL_KC_synapse_counts`, olfsysm 9881846, Python only): Pratyush's APL->KC table of
2025-12-14 (the one `connectome_APL_weights` reads), with the loader's filters (APL_R; both sites in the ROI), per KC of
the model's claw wiring (1,732 KCs, 9,867 claw rows, the same on all eight wirings): 'kept' (in the calyx on a (KC, claw)
of the wiring: the loader's weights exactly), 'dropped_claw' (in the calyx on a claw the wiring does not have), 'no_claw'
(in the calyx on no claw: claw ID -1 where the wiring has no claw -1 for the KC, or no claw ID; a^c_k), 'peduncle' (both
sites in PED(R); a^p_k). The component (E10_joint_fit.py's `aplkc_flat_shape` / `aplkc_flat_component`): a_k / n_k on
each of KC k's n_k claws, scaled to mean 1 over the claws and then to the connectome's mean (normalized) weight, so its
size d is in E95's unit (`E95.component_unit`: at d = 1 the mean claw's component takes, at the fit's scale and APL
activity, the fit's mean loss a claw at the boutons). A claw's APL->KC weight is the family's plus d x unit x the
component. Modes: 'offclaw' a_k = a^c_k (B1), 'offclaw_ped' a_k = a^c_k + a^p_k (B2 to B5), 'allsyn' (B6) each kept claw's
own count plus B2's a_k / n_k, in synapses (the proportion in which they are counted), scaled the same way.

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B, C and D refuse to run until it passes):
  - zero_real: the weights with a component of size 0 at dose 0.1 are the family's and reproduce E94's saved point
    B_sub_k2_d0.1 on the real wiring (its factor and its counts, 0 differences);
  - uniform_real: the uniform mode at d = 1 through this script's weights reproduces E95's saved D1_uniform_d1 on the real
    wiring (factor and counts);
  - pool_real: the pooled share s = 1 with no component at the fit's dose (1) reproduces E91's saved B1_F_k2_s1 on the
    real wiring (factor and counts);
  - counts_<wiring> (all eight): the per-claw kept counts are the loader's weights exactly and proportional to the
    wiring's normalized APL->KC weights; the loader's own dropped counts (`connectome_APL_weights(dropped_counts=)`)
    are the counts' no_claw and dropped_claw per KC; the totals beside the entry's probe (7,479 on no claw, 14,576 on
    kept claws); the KCs, claws and weights the same on every wiring and as E95's a4_weights.csv (claws without a weight,
    KCs without one); the components of the three modes (mean 1, a KC's sum proportional to a_k). Each wiring's per-KC
    and per-claw tables (counts, b_k, c_k, the components) are cached under outputs/E96_offclaw_apl/W/ for part A.
  - `--part Dcheck` (run once before part D): E10_joint_fit.py with --aplkc-flat offclaw (and offclaw_ped, allsyn
    beside) at size 0 against E94's D1 fit on the real wiring's saved first evaluation (d_check.json).
- A (no simulations): A1 the counts (a1_counts.csv: the entry's table, exact, for the model's KCs, on the real wiring's
  b_k; a1_by_type_class.csv by KC type and claw class; a1_wirings.csv a^c_k, a^p_k and their sum against the claw count
  and E95's b_k and c_k on each of the eight wirings: r, the log-log slope, the spread within a claw count;
  a1_within.csv per claw count); A2 where the synapses on no claw lie (a2_distances.csv: in Pratyush's assignments, each
  APL->KC calyx synapse's distance to its KC's nearest PN->KC synapse on a claw, by kind; in the derived tables of four
  connectomes, the APL->KC calyx synapses' distance to the KC's nearest derived claw, the share within 2 um, per KC the
  synapses off the claws against the claw count, a2_connectomes.csv; read only); A3 the KC->APL synapses split the same
  way (a3_kc2apl.csv); A4 the component's shape among E95's static forms and its uniform and connectome components
  (a4_shapes.csv: the per-KC sum's log-log slope on the claw count, the spread within a claw count, the per-claw spread);
  a_summary.json.
- B (`--arm` B1..B6, default all; `--wirings`, default S; resumable, runs_B.csv; runs under outputs/E96_offclaw_apl/runs/
  with their stage matrices, each component per KC and the APL->KC weights per claw): at fixed scales on S with the
  bouton path at dose 0.1 (E94's anatomical dose) and the tied APL scale (APL->KC, KC->APL, bouton->APL) bisected to the
  wiring's APL-on count in E92's bracket (0.01 to 10); B1 offclaw and B2 offclaw_ped at d = 0.5, 1, 2; B3 B2 with E91's
  pooled share s = 1 (`claw_inh_pool`) at the three sizes, with B3_pool1 (s = 1, no component) its reference; B4 B2 at
  d = 1 with the counts permuted among the KCs with the same number of claws and over all KCs (seed 1); B5 B2 at d = 1
  without the bouton path (dose 0); B6 allsyn at the three sizes. 64 points. The KC-KC scale and the somatic threshold
  stay at the fit's; the APL-off run of every point is the saved fit's.
- score (`--no-co`, `--no-zcell`): readings.csv per variant and unit (each wiring; S and T where all four have a run) and
  readout ('cal' the calibrated readout, K set on each run's APL-on counts, the reading of the choice and the gate;
  'calK' the fit's K; 'raw'; 'zcell' and 'zcell_lin' E73's per-KC reading under the saturating and the linear map) against
  every target of `targets89` (E92's `score_row92`: loss, pattern part with its amplitude and shape parts, spread ratio,
  Pearson r and ρ on the per-pair APL effect, the two groups' errors, the mean APL effect; E94's group effects); E95's
  saved part D runs on S (D1_uniform_d*, D2_connectome_d*, D3_uniform_d1_dose0) and E94's D0 are read beside as the
  comparison; pairs.csv; line.csv (E92's part A line); coresponse.csv (E89's readings of the alcohol pairs' co-response
  change); then the tables of part C (choice.csv, gate.csv, equal.csv), part D's condition and verdict.json, the figure
  and comparison_table.md.
- C: the choice per half (`choose96`; among B1, B2, B3, B5, B6 and the unchanged fits), each chosen setting with
  something to confirm run on T (runs_C.csv), with its controls on T: the rule (v) control E95's uniform component at the
  same size, pooled share and dose (`control_of96`: D1_uniform_d<d>, D3_uniform_d1_dose0, or U_uniform_d<d>_pool1, which
  is run on S as well), and beside it the same setting with its counts permuted among KCs with the same number of claws
  (<v>_permclaw1) and the bouton path at a tenth without a component (D0); then `score` without the co-response. The gate
  is E83's as E91 applied it (E92's `gate_direction92` / `gate92`), with the beside-controls' gains in gate.csv.
- Dcheck: see part 0.
- D (`--pop S|T`, `--start own|e95`, `--wirings` for a lane, `--dry-run`, `--force`): the family fitted again at the
  connectome's bouton share with the component (`E10_joint_fit.py -k 2 --aplpn-hold F --aplkc-flat <mode>
  [--claw-inh-pool 1] --aplkc-flat-unit U --aplkc-flat-max 8 --rate-target ctrl_frac --loss apl --apl-readout both`,
  200 evaluations; F = 0.1 x the wiring's fitted APL->bouton multiple of E3's base; U the wiring's `component_unit`),
  only on its condition (`d_plan`): S's best setting of B1 or B2 within 0.008 of the unchanged fits' 0.0500; the mode
  the one of B1 and B2 with the lower pattern part at its best size (B1 if within 0.0005); the pooled share 1 only if
  B3 is lower than B2 at the mode's best size by 0.0015 or more; B6 (`--mode allsyn`) only on an equal of B6 at fixed
  scales. Two starts per S wiring: 'own' (the wiring's best point of the fitted arm of part B; tag suffix _g1) and 'e95'
  (E95's refit of that wiring, its parameters and fitted size with the component exchanged; _g1e). T's fits (`--pop T`,
  start 'own' only) only on an equal of S's fitted population: from each wiring's own fit with the APL->bouton scale at
  dose 0.1 and the component at the median of S's fitted sizes.
- Dscore: per S wiring the fit with the lower final value of its own loss (both reported, d_fits.csv); S's population
  by E94's D1 rule (equal / better), beside E95's refit (as good as the uniform component: within 0.0015 of its 0.0530,
  or lower) and E94's 0.0593; the eight against 0.0482 (E91's part D pass rule). d_scores.csv, d_scores.json.
- E: documented stubs (`part_e`): runs only by the entry's rules.
- steps (`--variants`): the step from the claws to the somata in wild-type flies for chosen or equal settings (E94's
  practice, E95's `part_steps`).
- variants: variants.csv and model_pair_distances_long.csv for the archive (`archive_tools/archive_model_outputs.py
  E96_offclaw_apl`, then `archive_tools/migrate_model_output_params.py --apply`).
- table: comparison_table.md and the figure from what exists.

Deviations from the entry, forced by the code or the data:
- The counts are read with the loader's filters, both synaptic sites in the ROI (`filter_synapses_to_roi`). The entry's
  probe (7,479 on no claw, 14,576 on kept claws, 2,820 on dropped claws, 8,713 in the peduncle) appears to have used the
  presynaptic site alone and to have counted the synapses with claw ID -1 that the loader keeps on the wiring's 62 claws
  with ID -1 (Pratyush's claw of spatially unclustered synapses from one PN) as on no claw. With the loader's filters the
  loader keeps 14,808 (14,570 on claws with an ID of 0 or more and 238 with ID -1 on the wiring's -1 claws), drops 7,234
  on no claw (7,233 with ID -1, 1 without an ID) and 2,816 on claws the wiring drops; 8,701 lie in the peduncle. Part 0
  reports both; the components use the loader's split (the 238 are in the model's weights already, so they are not
  added again).
- The pooled share's rule ("B3 lower than B2 at the same size by 0.0015 or more") is read at the size of the fitted
  mode's best point (B3 has the offclaw_ped counts; with mode offclaw the same comparison is used for the pool).
- B3_pool1 (s = 1, no component) is B3's reference, not a setting of the choice or of the equal rule.
- Part D's second start's tag suffix is _g1e (the first's _g1, as the entry): E10_joint_fit.py's tag must differ between
  two fits of one wiring. T's fits start from each wiring's own fit with the APL->bouton scale at dose 0.1 ('own').
- The uniform control with the pooled share (U_uniform_d<d>_pool1) is not one of E95's runs; it is run here on S and T
  when B3 is chosen.
- In B3's runs the shares of the claws' loss (at the boutons, by the component, by the weighted part) are of the claws'
  clipped states, which the pooled share does not change; what the pooled share takes beyond them is not in the shares.

Outputs: results/E96_offclaw_apl/; spike counts, responses, stage matrices, components and weights under
outputs/E96_offclaw_apl/ (runs/, W/, D/ logs); part D's fits where E10_joint_fit.py writes them for this experiment.

Usage (from the repository root; one 16-thread process at a time on viki; check `uptime` first; long runs with -u):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part 0
    ./venv/bin/python -u -m experiments.E96_offclaw_apl --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part Dcheck
  part B (64 points, about 1 h; resumable, one arm at a time or all):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part B \
        > outputs/E96_offclaw_apl/B.log 2>&1
  then the readings (co-response readings in 4 processes):
    OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part score
  (C) only on a choice (it reads choice.csv; nothing runs when nothing is to confirm):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part C
  part D, only on its condition (`--dry-run` prints the plan and the commands without running):
    ./venv/bin/python -u -m experiments.E96_offclaw_apl --part D --dry-run
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part D \
        --pop S --start own > outputs/E96_offclaw_apl/D_own.log 2>&1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part D \
        --pop S --start e95 > outputs/E96_offclaw_apl/D_e95.log 2>&1
  (`--wirings real,P1` runs a lane; add `--mode allsyn` for B6's fits on their rule.) On hal, the second starts as a
  lane: pull, `./venv/bin/pip install -e .` if needed, then
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E96_offclaw_apl --part D \
        --pop S --start e95 > outputs/E96_offclaw_apl/D_e95_hal.log 2>&1
  (its fits.csv rows and histories go to the shared results folder as viki's do; it needs viki's runs_B.csv and
  readings.csv there, and its own part 0 cache: run --part 0 --wirings real first, about 3 min), then copy its fits'
  outputs: `mkdir -p <shared>/data/model_outputs/E96_offclaw_apl/from_hal && cp outputs/E96_offclaw_apl/*_g1e_*.parquet
  <shared>/data/model_outputs/E96_offclaw_apl/from_hal/` (Dscore and variants read them there).
    ./venv/bin/python -u -m experiments.E96_offclaw_apl --part Dscore
    ./venv/bin/python -u -m experiments.E96_offclaw_apl --part steps
    ./venv/bin/python -u -m experiments.E96_offclaw_apl --part variants
    ./venv/bin/python archive_tools/archive_model_outputs.py E96_offclaw_apl
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

EXPERIMENT = 'E96_offclaw_apl'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

APL_BRACKET = E92G.APL_BRACKET                      # (0.01, 10): the tied APL scale, E92's bracket
APL_STEP = E83.APL_STEP
RAW_CRIT = E92G.RAW_CRIT                            # a model response is two or more spikes
CLAW_GROUPS = E94.CLAW_GROUPS                       # 1-4, 5-8, 9+ claws
D_DOSE = 0.1                                        # the bouton path at a tenth (E94's anatomical dose)
B_SIZES = (0.5, 1.0, 2.0)                           # the component's size d in E95's unit (as E95's D1 / D2)
B4_SIZE = B5_SIZE = 1.0
POOL = 1.0                                          # B3: E91's pooled share
PERM_SEED = 1
MODES = ('offclaw', 'offclaw_ped', 'allsyn')
MODE_CODE = {'offclaw': 'offclaw', 'offclaw_ped': 'offclawped', 'allsyn': 'allsyn'}
CODE_MODE = {v: k for k, v in MODE_CODE.items()}
ARM_MODE = {'B1': 'offclaw', 'B2': 'offclaw_ped', 'B3': 'offclaw_ped', 'B4': 'offclaw_ped', 'B5': 'offclaw_ped',
            'B6': 'allsyn'}
PROBE = dict(no_claw=7479, kept=14576, dropped_claw=2820, peduncle=8713)   # the entry's probe (estimates)
E91_POOL_REF = 'B1_F_k2_s1'                         # E91's saved point at s = 1 (part 0's identity)
E95_UNIFORM_REF = 'D1_uniform_d1'                   # E95's saved point (part 0's identity)


def _vname(arm: str, mode: str, d: float, pool: bool = False, dose0: bool = False, perm: Optional[str] = None) -> str:
    return (f'{arm}_{MODE_CODE[mode]}_d{d:g}' + ('_pool1' if pool else '') + ('_dose0' if dose0 else '')
            + (f'_perm{perm}{PERM_SEED}' if perm else ''))


ARMS: Dict[str, Tuple[str, ...]] = {
    'B1': tuple(_vname('B1', 'offclaw', d) for d in B_SIZES),
    'B2': tuple(_vname('B2', 'offclaw_ped', d) for d in B_SIZES),
    'B3': ('B3_pool1',) + tuple(_vname('B3', 'offclaw_ped', d, pool=True) for d in B_SIZES),
    'B4': tuple(_vname('B4', 'offclaw_ped', B4_SIZE, perm=p) for p in ('claw', 'all')),
    'B5': (_vname('B5', 'offclaw_ped', B5_SIZE, dose0=True),),
    'B6': tuple(_vname('B6', 'allsyn', d) for d in B_SIZES),
}
REF = 'ref'                                         # the unchanged fits (source 'saved')
D0 = E95.D0                                         # dose 0.1 without a component (E94's point on S)
B3_REF = 'B3_pool1'

# part C (the entry's rules: E91's, as in E94 and E95)
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
CONTROL_RATIO = E91.CONTROL_RATIO
HALVES = E91.HALVES
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
CHOICE_ARMS = ('B1', 'B2', 'B3', 'B5', 'B6')        # B4 (the permutations) and B3's reference excepted

# part D
D_COND = 0.008                                      # S's best of B1 or B2 within this of the unchanged fits' 0.0500
D_MODE_TIE = 0.0005                                 # B1 and B2 within this: B1
D_POOL_GAIN = 0.0015                                # B3 this far below B2 at the same size: the pooled share 1
D_FLAT_BOUNDS = (0.0, 8.0)                          # the component's size, a fifth free parameter
D_MAXFEV = 200
D_TAGS = {'own': '_g1', 'e95': '_g1e'}
D_START_FILE = 'd_start_params.json'
D_AS_GOOD = 0.0015                                  # as good as the uniform component's refit: within this, or lower
E94_REFIT_S = E95.E94_REFIT_S                       # 0.0593
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = E95.UNCHANGED_8                       # 0.0482
FLIES_MEAN = E89.FLIES_MEAN

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


def e95_refit_s() -> float:
    """E95's refit with the uniform component, S's pattern part against all recordings (f_scores.json; 0.0530)."""
    s = _read_json(H.results_dir(E95.EXPERIMENT) / 'f_scores.json').get('S_uniform', {})
    return float(s.get('pattern_all', 0.0530))


# ---------------------------------------------------------------------------
# names

_RE_B = re.compile(r'^(B[1-6])_(offclaw|offclawped|allsyn)_d([0-9.]+)(_pool1)?(_dose0)?(?:_perm(claw|all)(\d+))?$')
_RE_U = re.compile(r'^U_uniform_d([0-9.]+)_pool1$')
E95_REF_ARMS = ('D1', 'D2', 'D3')


def parse_variant(v: str) -> dict:
    """A variant's fields: arm, part ('B', 'ref', 'e95' E95's saved part D settings, 'control'), mode (the component:
    MODES, 'uniform', 'connectome' or None), size, pool (0 or 1), dose, permute ('claw' / 'all' / None) and seed, control
    and control_of, screen (the arm; 'B3ref'; 'ref'; 'e95'; 'control'), equal (whether the equal rule applies).
    Names: B<n>_<offclaw|offclawped|allsyn>_d<d>[_pool1][_dose0][_perm<claw|all><seed>] (B4 the permutations at d = 1;
    a permutation of another arm's setting is a control of part C), B3_pool1, ref, D0 (E94's point), E95's part D
    settings (D1_uniform_d<d>, D2_connectome_d<d>, D3_uniform_d1_dose0), U_uniform_d<d>_pool1 (the uniform control with
    the pooled share)."""
    out = dict(variant=v, arm=None, part=None, mode=None, size=0.0, pool=0.0, dose=D_DOSE, permute=None, seed=None,
               control=False, control_of='', screen=None, equal=False)
    if v == REF:
        out.update(arm=REF, part=REF, dose=1.0, screen=REF)
        return out
    if v == D0:
        out.update(arm=D0, part='control', control=True, control_of='C', screen='control')
        return out
    if v == B3_REF:
        out.update(arm='B3', part='B', pool=POOL, screen='B3ref')
        return out
    m = _RE_U.match(v)
    if m:
        d = m.group(1)
        if f'{float(d):g}' != d:
            raise ValueError(f'variant {v!r}: the size is not written canonically')
        out.update(arm='U', part='control', mode='uniform', size=float(d), pool=POOL, control=True, control_of='C',
                   screen='control')
        return out
    if v.startswith(E95_REF_ARMS):
        p = E95.parse_variant(v)
        if p['part'] != 'D':
            raise ValueError(f'variant {v!r}: not one of E95\'s part D settings')
        out.update(arm=f"e95_{p['arm']}", part='e95', mode=p['component'], size=float(p['size']), dose=float(p['dose']),
                   screen='e95')
        return out
    m = _RE_B.match(v)
    if not m:
        raise ValueError(f'variant {v!r}: not one of this experiment\'s')
    arm, code, d, pool, dose0, perm, seed = m.groups()
    if f'{float(d):g}' != d or float(d) <= 0:
        raise ValueError(f'variant {v!r}: the size is not written canonically or not positive')
    mode = CODE_MODE[code]
    out.update(arm=arm, part='B', mode=mode, size=float(d), pool=POOL if pool else 0.0, dose=0.0 if dose0 else D_DOSE,
               permute=perm, seed=None if seed is None else int(seed), screen=arm)
    core = v[:v.index('_perm')] if perm else v
    if arm == 'B4':
        if not (perm and mode == 'offclaw_ped' and float(d) == B4_SIZE and not pool and not dose0):
            raise ValueError(f'variant {v!r}: B4 is B2 at d = {B4_SIZE:g} with its counts permuted')
    elif perm:
        # a permutation of a setting of another arm: part C's control (its counts permuted)
        parse_variant(core)
        out.update(part='control', control=True, control_of=core, screen='control')
    if arm != 'B4' and not perm:
        want = {'B1': ('offclaw', False, False), 'B2': ('offclaw_ped', False, False), 'B3': ('offclaw_ped', True, False),
                'B5': ('offclaw_ped', False, True), 'B6': ('allsyn', False, False)}[arm]
        if (mode, bool(pool), bool(dose0)) != want or (arm == 'B5' and float(d) != B5_SIZE):
            raise ValueError(f'variant {v!r}: not a setting of {arm}')
    if arm != 'B4' and perm and arm not in ('B1', 'B2', 'B3', 'B6'):
        raise ValueError(f'variant {v!r}: no permuted control of {arm}')
    out['equal'] = bool(out['part'] == 'B' and arm in CHOICE_ARMS)
    return out


def control_of96(variant: str) -> str:
    """The rule (v) control of a chosen setting: E95's uniform component at the same size, pooled share and dose
    (D1_uniform_d<d> at dose 0.1, D3_uniform_d1_dose0 without the bouton path, U_uniform_d<d>_pool1 with the pooled
    share)."""
    pv = parse_variant(variant)
    if pv['part'] != 'B' or pv['screen'] not in CHOICE_ARMS:
        raise ValueError(f'{variant}: not a setting of the choice')
    if pv['dose'] == 0.0:
        return f'D3_uniform_d{pv["size"]:g}_dose0'
    if pv['pool'] > 0:
        return f'U_uniform_d{pv["size"]:g}_pool1'
    return f'D1_uniform_d{pv["size"]:g}'


def beside_controls(variant: str) -> List[str]:
    """The controls run on T beside rule (v)'s: the setting with its counts permuted among KCs with the same number of
    claws, and the bouton path at a tenth without a component (D0)."""
    return [f'{variant}_permclaw{PERM_SEED}', D0]


def arm_variants(arm: str) -> List[str]:
    return list(ARMS[arm])


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


# ---------------------------------------------------------------------------
# the counts and the component

def wiring_counts(label: str = 'real') -> Tuple[pd.DataFrame, pd.Series]:
    """(per KC, per claw) of a wiring's claws (E10's `aplkc_flat_counts` of its normalized APL->KC weights' index)."""
    from experiments import E10_joint_fit as J
    return J.aplkc_flat_counts(E83.wiring_parts(label)['wnorm']['wAPLKC'].index)


def offclaw_counts(mode: str, per_kc: pd.DataFrame, permute: Optional[str] = None, seed: int = PERM_SEED) -> pd.Series:
    """a_k of a mode per KC (E10's `aplkc_flat_offclaw`), permuted among the KCs with the same number of claws ('claw')
    or over all KCs ('all'; E95's `permute_values`, one generator from `seed`)."""
    from experiments import E10_joint_fit as J
    a = J.aplkc_flat_offclaw(mode, per_kc)
    if permute is None:
        return a
    if permute not in ('claw', 'all'):
        raise ValueError(f'permute {permute!r}')
    v = E95.permute_values(a.to_numpy(dtype=float), seed, per_kc['n_claws'].to_numpy() if permute == 'claw' else None)
    return pd.Series(v, index=a.index, name=a.name)


def component_shape(mode: str, wnorm: pd.Series, counts: Optional[Tuple[pd.DataFrame, pd.Series]] = None,
                    permute: Optional[str] = None, seed: int = PERM_SEED) -> np.ndarray:
    """The component per claw in units of the connectome's mean weight (mean 1 over the claws): 'uniform' 1 on every claw;
    'connectome' each claw's connectome weight over its mean; MODES E10's `aplkc_flat_shape` of the counts (`counts`,
    `wiring_counts` of `wnorm`'s claws when None), the off-claw counts permuted with `permute`."""
    from experiments import E10_joint_fit as J
    w = wnorm.to_numpy(dtype=float)
    if mode == 'uniform':
        return np.ones(w.size)
    if mode == 'connectome':
        return w / w.mean()
    if mode not in MODES:
        raise ValueError(f'unknown component {mode!r}')
    per_kc, per_claw = J.aplkc_flat_counts(wnorm.index) if counts is None else counts
    if not per_claw.index.equals(wnorm.index):
        raise ValueError("the counts' claws are not the weights' claws")
    kc_of = np.asarray(wnorm.index.get_level_values(H.KC_ID_LEVEL))
    a = offclaw_counts(mode, per_kc, permute, seed)
    return J.aplkc_flat_shape(mode, kc_of, a, per_claw.to_numpy(dtype=float))


def component_values(mode: str, wnorm: pd.Series, counts: Optional[Tuple[pd.DataFrame, pd.Series]] = None,
                     permute: Optional[str] = None, seed: int = PERM_SEED) -> np.ndarray:
    """The unscaled component per claw per unit of size and of the unit: E10's `aplkc_flat_component` (the connectome's
    mean weight times the shape; E95's uniform and connectome components exactly), or with the counts permuted the same
    built from `component_shape`."""
    from experiments import E10_joint_fit as J
    if mode in ('uniform', 'connectome') or (mode in MODES and permute is None):
        cnt = None
        if mode in MODES:
            cnt = J.aplkc_flat_counts(wnorm.index) if counts is None else counts
        return J.aplkc_flat_component(wnorm, mode, cnt)
    return float(wnorm.to_numpy(dtype=float).mean()) * component_shape(mode, wnorm, counts, permute, seed)


def e96_weights(w_fam: pd.Series, wnorm: pd.Series, comp: np.ndarray, d: float, unit: float) -> pd.Series:
    """The APL->KC weights per claw (unscaled): the family's (`w_fam`) plus d x unit x `comp` (`component_values`), as
    E95's `d_weights`; d = 0 gives `w_fam` itself."""
    if not wnorm.index.equals(w_fam.index):
        raise ValueError('the connectome weights are not indexed as the family\'s')
    comp = np.asarray(comp, dtype=float).ravel()
    if comp.size != len(w_fam) or not np.isfinite(comp).all() or (comp < 0).any():
        raise ValueError(f'component: {comp.size} values for {len(w_fam)} claws, finite and not negative')
    if not (np.isfinite(d) and d >= 0):
        raise ValueError(f'component size {d}')
    if d == 0:
        return w_fam.copy()
    return pd.Series(w_fam.to_numpy(dtype=float) + float(d) * float(unit) * comp, index=w_fam.index, name=w_fam.name)


def synapse_share(mode: str, per_kc: pd.DataFrame) -> float:
    """The synapses' share of the component's mode: the off-claw synapses it adds over those plus the loader's kept ones
    (allsyn: 1 minus nothing; reported as the share of the off-claw part of its counts)."""
    a = float((per_kc['no_claw'] + (per_kc['peduncle'] if mode != 'offclaw' else 0)).sum())
    k = float(per_kc['kept'].sum())
    return a / (a + k) if a + k > 0 else np.nan


# ---------------------------------------------------------------------------
# the simulator: the weights, the dose and the pooled share

def reading96(model: E95.StaticModel, pv: dict, wk: Optional[pd.Series], label: str = '') -> dict:
    """One setting on a wiring at its fit (E95's `reading95` in its part D mode): the APL->KC weights `wk` (None: the
    fit's), the dose on the APL->bouton scale, the pooled share (olfsysm's `claw_inh_pool`, E91's `set_pool`); the tied
    APL scale bisected to the fit's APL-on count in APL_BRACKET, the run again with APL's activity saved, its stage
    matrices and E92's `gain_stats`. The pooled share is reset on return, also on an error."""
    try:
        E91.set_pool(model.mp, float(pv['pool']))
        return E95.reading95(model, dict(mode='kcs', dose=float(pv['dose'])), None, wk, label=label)
    finally:
        E91.set_pool(model.mp, 0.0)


# ---------------------------------------------------------------------------
# saving and loading runs

STAGE_KEYS = E94.STAGE_KEYS


def _dirs(sub: str = 'runs') -> List[Path]:
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [H.outputs_dir(EXPERIMENT) / sub, arc / sub, arc / 'from_hal' / sub]


def save_run96(tag: str, o: Dict[str, tuple], sub: str = 'runs', stages: Optional[dict] = None,
               drive: Optional[pd.DataFrame] = None, loss: Optional[pd.DataFrame] = None,
               weights: Optional[pd.Series] = None, component: Optional[pd.DataFrame] = None) -> None:
    """A run's outputs (APL on, and the saved APL-off run beside), its stage matrices with APL on, each KC's summed claw
    drive from the boutons, E92's loss of drive at the claws, the APL->KC weights per claw where they are not the fit's
    and the component per KC, under outputs/E96_offclaw_apl/<sub>/."""
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


def load_run96(tag: str, sub: str = 'runs') -> Optional[Dict[str, tuple]]:
    for d in _dirs(sub):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def load_stage96(tag: str, sub: str = 'runs') -> Optional[dict]:
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


def _frame(a: np.ndarray, like: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(np.asarray(a, dtype=float), index=like.index, columns=[str(c) for c in like.columns])


def _counts_diff(a, b) -> dict:
    return E94._counts_diff(a, b)


# ---------------------------------------------------------------------------
# the per-KC and per-claw tables (part 0's cache, part A's input)

def _w_dir() -> Path:
    return local_dir('W')


def kc_tables96(label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per KC, per claw) of a wiring at its fit: per KC (the outputs' rows) its type, claws, the counts (kept,
    kept_unassigned, dropped_claw, no_claw, no_claw_nan, peduncle), E95's b_k and c_k (`kc_losses`), each mode's
    component summed over the KC's claws (units of the connectome's mean weight) and the connectome's and the family's
    summed weights; per claw its KC, the kept count, the connectome's and the family's weights and each mode's
    component."""
    P = E83.wiring_parts(label)
    wn = P['wnorm']['wAPLKC']
    per_kc, per_claw = wiring_counts(label)
    saved_on = E83.saved_fit(label)['saved']['on'][1]
    kc_index = saved_on.index
    ids = np.asarray(kc_index.get_level_values(H.KC_ID_LEVEL))
    L = E95.kc_losses(label, kc_index)
    from experiments import E10_joint_fit as J
    w_fam = J.breadth_weights(wn, P['breadth'], E83.BREADTH_K)
    kc_of = np.asarray(wn.index.get_level_values(H.KC_ID_LEVEL))
    claws = pd.DataFrame({H.KC_ID_LEVEL: kc_of, 'claw_id': np.asarray(wn.index.get_level_values('claw_id')),
                          'kept': per_claw.to_numpy(dtype=float), 'wnorm': wn.to_numpy(dtype=float),
                          'w_fam': w_fam.to_numpy(dtype=float)})
    for mode in ('uniform', 'connectome') + MODES:
        claws[f'comp_{mode}'] = component_shape(mode, wn, (per_kc, per_claw))
    g = claws.groupby(H.KC_ID_LEVEL)
    kc = per_kc.reindex(ids).reset_index()
    if kc.isna().any().any():
        raise RuntimeError(f'{label}: the counts do not cover the outputs\' KCs')
    kc.insert(1, H.KC_TYPE_LEVEL, np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL)))
    kc['b'] = L['b'].to_numpy(dtype=float)
    kc['c'] = L['c'].to_numpy(dtype=float)
    for mode in ('uniform', 'connectome') + MODES:
        kc[f'comp_{mode}_sum'] = g[f'comp_{mode}'].sum().reindex(ids).to_numpy(dtype=float)
    kc['wnorm_sum'] = g.wnorm.sum().reindex(ids).to_numpy(dtype=float)
    kc['w_fam_sum'] = g.w_fam.sum().reindex(ids).to_numpy(dtype=float)
    return kc, claws


def save_kc_tables96(label: str, kc: pd.DataFrame, claws: pd.DataFrame) -> None:
    kc.to_parquet(_w_dir() / f'{label}_kc.parquet')
    claws.to_parquet(_w_dir() / f'{label}_claw.parquet')


def load_kc_tables96(label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ps = [_w_dir() / f'{label}_{k}.parquet' for k in ('kc', 'claw')]
    if not all(p.exists() for p in ps):
        raise RuntimeError(f'{label}: no per-KC table (run --part 0 first: it caches them under {_w_dir()})')
    return tuple(pd.read_parquet(p) for p in ps)


# ---------------------------------------------------------------------------
# part 0: the identities

def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('passes'))


def counts_identity(label: str) -> dict:
    """The counts' identities on one wiring (module docstring, part 0): the per-claw kept counts against the loader's
    weights and the wiring's normalized weights, the loader's dropped counts against the counts, the totals beside the
    entry's probe, the claws and weights against the real wiring's and E95's a4_weights.csv, the components."""
    from pyolfsysm.mb_model import connectome_APL_weights, connectome_wPNKC, connectome_APL_KC_synapse_counts
    P = E83.wiring_parts(label)
    wn = P['wnorm']['wAPLKC']
    per_kc, per_claw = wiring_counts(label)
    res = dict(wiring=label, n_kcs=int(len(per_kc)), n_claws=int(len(per_claw)))
    totals = {c: int(per_kc[c].sum()) for c in per_kc.columns if c != 'n_claws'}
    res['totals'] = totals
    res['kept_on_claws_id_ge0'] = int(totals['kept'] - totals['kept_unassigned'])
    res['probe'] = dict(PROBE)
    # the loader's weights and drops, on the wiring the model builds (the base's claw wiring; the claws are the same on
    # every wiring)
    if 'loader' not in _CACHE:
        w0 = connectome_wPNKC(prat_claws=True, prat_boutons=True)
        dropped: dict = {}
        wa = connectome_APL_weights(wPNKC=w0, prat_claws=True, prat_boutons=True, dropped_counts=dropped)[0]
        _CACHE['loader'] = (w0, wa, dropped)
    w0, wa, dropped = _CACHE['loader']
    def claw_keys(ix):
        return ix.to_frame(index=False)[[H.KC_ID_LEVEL, 'claw_id']].to_numpy()
    same_index = bool(np.array_equal(claw_keys(w0.index), claw_keys(wn.index)))
    res['wiring_index_is_loaders'] = same_index
    res['wiring_index_equals_loaders'] = bool(w0.index.equals(wn.index))
    res['kept_is_loaders'] = bool(same_index and np.array_equal(per_claw.to_numpy(), wa.to_numpy(dtype=float)))
    k = per_claw.to_numpy(dtype=float)
    w = wn.to_numpy(dtype=float)
    scale = float(w.sum() / k.sum())
    res['wnorm_over_kept'] = scale
    res['wnorm_proportional_to_kept'] = bool(np.allclose(w, k * scale, rtol=1e-12, atol=1e-12))
    d = dropped['wAPLKC'].reindex(per_kc.index)
    res['loader_dropped'] = {c: int(d[c].sum()) for c in d.columns}
    res['loader_drops_are_counts'] = bool((d.nan_claw + d.unassigned_claw == per_kc.no_claw).all()
                                          and (d.other_claw == per_kc.dropped_claw).all()
                                          and (d.nan_claw == per_kc.no_claw_nan).all())
    # the presynaptic site alone (the probe's apparent rule), for the record
    if 'pre_only' not in _CACHE:
        from pyolfsysm.mb_model import prat_hemibrain_seg_dir, APL_KC_SYNAPSE_TABLES
        t = pd.read_parquet(prat_hemibrain_seg_dir / APL_KC_SYNAPSE_TABLES['APL>KC'])
        pre = t.copy()
        pre['roi_post'] = pre['roi_pre']
        a, _ = connectome_APL_KC_synapse_counts(wn.index, 'APL>KC', synapses=pre)
        _CACHE['pre_only'] = {c: int(a[c].sum()) for c in a.columns if c != 'n_claws'}
    pre = _CACHE['pre_only']
    res['presynaptic_site_only'] = dict(pre, no_claw_with_kept_unassigned=int(pre['no_claw'] + pre['kept_unassigned']),
                                        kept_on_claws_id_ge0=int(pre['kept'] - pre['kept_unassigned']))
    # the same claws and weights on every wiring; E95's a4_weights.csv
    if label != 'real':
        r = E83.wiring_parts('real')['wnorm']['wAPLKC']
        res['same_as_real'] = bool(r.index.equals(wn.index) and np.array_equal(r.to_numpy(), wn.to_numpy()))
    else:
        res['same_as_real'] = True
    try:
        a4 = pd.read_csv(H.results_dir(E95.EXPERIMENT) / 'a4_weights.csv').set_index('wiring').loc[label]
        zero = float((w == 0).mean())
        kc_zero = float((per_kc['kept'] == 0).mean())
        res['e95_a4'] = dict(n_claws=int(a4['n_claws']), zero_share=float(a4['zero_share']), zero_share_here=zero,
                             kc_zero_share=float(a4['kc_zero_share']), kc_zero_share_here=kc_zero)
        res['same_as_e95_a4'] = bool(int(a4['n_claws']) == len(w) and abs(zero - float(a4['zero_share'])) < 1e-12
                                     and abs(kc_zero - float(a4['kc_zero_share'])) < 1e-12)
    except (OSError, KeyError) as e:
        res['same_as_e95_a4'] = None
        res['e95_a4_note'] = repr(e)
    # the components: mean 1 over the claws; a KC's sum in proportion to a_k (and, for allsyn, a_k + its kept count)
    from experiments import E10_joint_fit as J
    kc_of = np.asarray(wn.index.get_level_values(H.KC_ID_LEVEL))
    comps = {}
    for mode in MODES:
        sh = component_shape(mode, wn, (per_kc, per_claw))
        a = J.aplkc_flat_offclaw(mode, per_kc)
        tot = a + (per_kc['kept'] if mode == 'allsyn' else 0.0)
        s = pd.Series(sh).groupby(kc_of).sum().reindex(per_kc.index)
        ratio = (s / tot)[tot > 0]
        comp = J.aplkc_flat_component(wn, mode, (per_kc, per_claw))
        comps[mode] = dict(mean=float(sh.mean()), kc_sum_proportional=bool(np.allclose(ratio, ratio.iloc[0], rtol=1e-9)),
                           kc_zero=int((s == 0).sum()), claws_zero=int((sh == 0).sum()), cv_claws=float(sh.std() / sh.mean()),
                           e10_component_mean_over_wnorm_mean=float(comp.mean() / w.mean()),
                           e10_is_shape=bool(np.allclose(comp, w.mean() * sh, rtol=1e-12)))
    res['components'] = comps
    res['passes'] = bool(res['kept_is_loaders'] and res['wnorm_proportional_to_kept'] and res['loader_drops_are_counts']
                         and res['same_as_real'] and res['same_as_e95_a4'] is not False
                         and all(abs(c['mean'] - 1) < 1e-12 and c['kc_sum_proportional'] and c['e10_is_shape']
                                 for c in comps.values()))
    return res


def part_0(out: Path, wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per key); raises if one fails."""
    import gc
    cc = _commits()
    path = out / 'identity.json'
    res = _read_json(path)
    res['commits'] = cc
    for w in wirings:
        key = f'counts_{w}'
        if not res.get(key, {}).get('passes'):
            t0 = time.time()
            res[key] = dict(counts_identity(w), seconds=time.time() - t0)
            print(f"E96 0 {key}: {json.dumps(res[key], default=str)[:1500]}", flush=True)
            save_kc_tables96(w, *kc_tables96(w))
            _write_json(path, res)
    if 'real' in wirings and not all(res.get(k, {}).get('passes') for k in ('zero_real', 'uniform_real', 'pool_real')):
        m = E95.StaticModel('real', check=False)
        L = E95.kc_losses('real', m.saved['on'][1].index)
        wn = E83.wiring_parts('real')['wnorm']['wAPLKC']
        unit = E95.d_unit(m, L)
        # (1) size 0 at dose 0.1: E94's B_sub_k2_d0.1
        if not res.get('zero_real', {}).get('passes'):
            t0 = time.time()
            comp = component_values('offclaw_ped', wn)
            wk = e96_weights(m.wAPLKC, wn, comp, 0.0, unit)
            r94 = pd.read_csv(H.results_dir(E94.EXPERIMENT) / 'runs_B1.csv').drop_duplicates('tag', keep='last')
            r94 = r94.set_index('tag').loc[f'{E95.E94_D0_TAG}_real']
            o94 = E94.load_run(f'{E95.E94_D0_TAG}_real')
            if o94 is None:
                raise RuntimeError(f"E94's point {E95.E94_D0_TAG} on the real wiring is not saved")
            rd = reading96(m, dict(pool=0.0, dose=D_DOSE), wk, label='identity zero real')
            dd = _counts_diff(rd['on'], o94['on'])
            f94 = float(r94['apl_factor'])
            ok = dict(same_weights=bool(np.array_equal(wk.to_numpy(), m.wAPLKC.to_numpy())), counts=dd,
                      apl_factor=rd['factor'], apl_factor_e94=f94, factor_rel_diff=abs(rd['factor'] - f94) / f94,
                      n_sims=rd['n_sims'], unit=unit, weights_reset=bool(np.array_equal(m._kc, m._kc0)),
                      pool_after=float(m.mp.kc.claw_inh_pool), seconds=time.time() - t0)
            ok['passes'] = bool(ok['same_weights'] and dd['exact'] and ok['factor_rel_diff'] <= 1e-12 and ok['weights_reset']
                                and ok['pool_after'] == 0.0)
            res['zero_real'] = ok
            print(f"E96 0 zero_real: {ok}", flush=True)
            _write_json(path, res)
        # (2) the uniform mode at d = 1: E95's D1_uniform_d1
        if not res.get('uniform_real', {}).get('passes'):
            t0 = time.time()
            comp = component_values('uniform', wn)
            wk = e96_weights(m.wAPLKC, wn, comp, 1.0, unit)
            w95 = E95.d_weights(m.wAPLKC, wn, 'uniform', 1.0, unit)
            r95 = pd.read_csv(H.results_dir(E95.EXPERIMENT) / 'runs_D.csv').drop_duplicates('tag', keep='last')
            r95 = r95.set_index('tag').loc[f'{E95_UNIFORM_REF}_real']
            o95 = E95.load_run95(f'{E95_UNIFORM_REF}_real')
            if o95 is None:
                raise RuntimeError(f"E95's point {E95_UNIFORM_REF} on the real wiring is not saved")
            rd = reading96(m, dict(pool=0.0, dose=D_DOSE), wk, label='identity uniform real')
            dd = _counts_diff(rd['on'], o95['on'])
            f95 = float(r95['apl_factor'])
            ok = dict(same_weights_as_e95=bool(np.array_equal(wk.to_numpy(), w95.to_numpy())), counts=dd,
                      apl_factor=rd['factor'], apl_factor_e95=f95, factor_rel_diff=abs(rd['factor'] - f95) / f95,
                      unit=unit, unit_e95=float(r95['component_unit']), n_sims=rd['n_sims'], seconds=time.time() - t0)
            ok['passes'] = bool(ok['same_weights_as_e95'] and dd['exact'] and ok['factor_rel_diff'] <= 1e-12
                                and abs(unit - ok['unit_e95']) <= 1e-12 * unit)
            res['uniform_real'] = ok
            print(f"E96 0 uniform_real: {ok}", flush=True)
            _write_json(path, res)
        # (3) the pooled share s = 1, no component, at the fit's dose: E91's B1_F_k2_s1
        if not res.get('pool_real', {}).get('passes'):
            t0 = time.time()
            r91 = pd.read_csv(H.results_dir(E91.EXPERIMENT) / 'runs_B1.csv').drop_duplicates('tag', keep='last')
            r91 = r91.set_index('tag').loc[f'{E91_POOL_REF}_real']
            o91 = E91.load_run(f'{E91_POOL_REF}_real')
            if o91 is None:
                raise RuntimeError(f"E91's point {E91_POOL_REF} on the real wiring is not saved")
            rd = reading96(m, dict(pool=POOL, dose=1.0), None, label='identity pool real')
            dd = _counts_diff(rd['on'], o91['on'])
            f91 = float(r91['apl_factor'])
            ok = dict(counts=dd, apl_factor=rd['factor'], apl_factor_e91=f91, factor_rel_diff=abs(rd['factor'] - f91) / f91,
                      n_sims=rd['n_sims'], pool_after=float(m.mp.kc.claw_inh_pool), seconds=time.time() - t0,
                      note="E91 bisected in E83's bracket (0.1, 10), this run in E92's (0.01, 10); the path from 1 is the "
                           'same unless a bound is reached')
            ok['passes'] = bool(dd['exact'] and ok['factor_rel_diff'] <= 1e-12 and ok['pool_after'] == 0.0)
            res['pool_real'] = ok
            print(f"E96 0 pool_real: {ok}", flush=True)
            _write_json(path, res)
        del m
        gc.collect()
    keys = [f'counts_{w}' for w in ALL_WIRINGS] + ['zero_real', 'uniform_real', 'pool_real']
    have_all = all(k in res for k in keys)
    failed = [k for k in keys if k in res and not res[k].get('passes')]
    res['passes'] = bool(have_all and not failed)
    res['complete'] = bool(have_all)
    res['failed'] = failed
    _write_json(path, res)
    if failed:
        raise RuntimeError(f'E96 part 0: the identity fails for {failed} (identity.json); the experiment stops here')
    print(f"E96 0: {'passes' if res['passes'] else 'so far exact (not all keys yet)'}", flush=True)
    return res


# ---------------------------------------------------------------------------
# part A (no simulations)

A1_KINDS = {'kept': ('kept',), 'kept_id_ge0': ('kept', '-kept_unassigned'), 'kept_unassigned': ('kept_unassigned',),
            'dropped_claw': ('dropped_claw',), 'no_claw': ('no_claw',), 'peduncle': ('peduncle',),
            'no_claw+peduncle': ('no_claw', 'peduncle')}
A1_NAMES = {'kept': 'calyx, on a claw that the wiring keeps (the model\'s weights)',
            'kept_id_ge0': 'of which on claws with an ID (0 or more)',
            'kept_unassigned': 'of which with claw ID -1, on the wiring\'s claw -1',
            'dropped_claw': 'calyx, on a claw that the wiring drops', 'no_claw': 'calyx, on no claw (a^c_k)',
            'peduncle': 'peduncle (a^p_k)', 'no_claw+peduncle': 'on no claw, calyx and peduncle together',
            'b': 'for comparison: what the bouton path takes from a KC in the fit (b_k, drive units)',
            'c': 'for comparison: what the claws take (c_k, drive units)'}


def _kind_values(kc: pd.DataFrame, kind: str) -> np.ndarray:
    if kind in ('b', 'c'):
        return kc[kind].to_numpy(dtype=float)
    v = np.zeros(len(kc))
    for c in A1_KINDS[kind]:
        v = v - kc[c[1:]].to_numpy(dtype=float) if c.startswith('-') else v + kc[c].to_numpy(dtype=float)
    return v


def _loglog(n: np.ndarray, v: np.ndarray) -> float:
    ok = (v > 0) & (n > 0)
    return float(np.polyfit(np.log(n[ok]), np.log(v[ok]), 1)[0]) if ok.sum() > 2 else np.nan


def _within(n: np.ndarray, v: np.ndarray, b: Optional[np.ndarray] = None, counts=range(2, 11)) -> dict:
    """The spread of v within a claw count (CV over the KCs with that count, counts with 10 KCs or more) and, with b,
    the correlation of v with b within a claw count: the range and the median over the counts."""
    cvs, rs = [], []
    for k in counts:
        m = n == k
        if m.sum() < 10:
            continue
        x = v[m]
        if x.mean() > 0:
            cvs.append(x.std() / x.mean())
        if b is not None:
            rs.append(E95._corr(x, b[m])[0])
    rs = [r for r in rs if np.isfinite(r)]
    return dict(within_cv_min=float(min(cvs)) if cvs else np.nan, within_cv_max=float(max(cvs)) if cvs else np.nan,
                within_cv_median=float(np.median(cvs)) if cvs else np.nan,
                within_r_b_min=float(min(rs)) if rs else np.nan, within_r_b_max=float(max(rs)) if rs else np.nan,
                within_r_b_median=float(np.median(rs)) if rs else np.nan)


def a1_row(kc: pd.DataFrame, kind: str) -> dict:
    """One row of A1's table: the synapses, the mean per KC, the CV over KCs, the share of KCs with none, r with the
    claw count, the log-log slope on it, the means by claw class, r with b_k (the wiring's), the spread within a claw
    count."""
    v = _kind_values(kc, kind)
    n = kc['n_claws'].to_numpy(dtype=float)
    b = kc['b'].to_numpy(dtype=float)
    row = dict(kind=kind, name=A1_NAMES[kind], synapses=float(v.sum()) if kind not in ('b', 'c') else np.nan,
               per_kc_mean=float(v.mean()), cv_kcs=float(v.std() / v.mean()) if v.mean() > 0 else np.nan,
               kcs_none_share=float((v == 0).mean()), r_n_claws=E95._corr(n, v)[0], loglog_slope_n=_loglog(n, v),
               r_b=E95._corr(b, v)[0], r_c=E95._corr(kc['c'].to_numpy(dtype=float), v)[0])
    for g, lo, hi in CLAW_GROUPS:
        m = (n >= lo) & (n <= hi)
        row[f'mean_{g}'] = float(v[m].mean()) if m.any() else np.nan
    row.update(_within(n, v, b))
    return row


def a1_counts(out: Path) -> Dict[str, pd.DataFrame]:
    """A1 (module docstring): a1_counts.csv (the entry's table, exact, on the real wiring's b_k), a1_by_type_class.csv,
    a1_wirings.csv (every wiring: each count against the claw count and that wiring's b_k and c_k), a1_within.csv (per
    claw count), a1_per_claw.csv (the per-claw spread of each component and of b_k / n_k)."""
    kc, claws = load_kc_tables96('real')
    rows = [a1_row(kc, k) for k in list(A1_KINDS) + ['b', 'c']]
    t = pd.DataFrame(rows)
    calyx = float(kc[['kept', 'dropped_claw', 'no_claw']].to_numpy().sum())
    t['share_of_calyx'] = t.synapses / calyx
    t['share_of_calyx_and_peduncle'] = t.synapses / (calyx + float(kc.peduncle.sum()))
    t.to_csv(out / 'a1_counts.csv', index=False)
    # by KC type and claw class
    trows = []
    grp = E95._claw_group(kc['n_claws'].to_numpy(dtype=float))
    for by, keys in ((H.KC_TYPE_LEVEL, kc[H.KC_TYPE_LEVEL].astype(str).to_numpy()), ('claw_class', grp)):
        for g in pd.unique(keys):
            m = keys == g
            r = dict(by=by, group=str(g), n_kcs=int(m.sum()), n_claws_mean=float(kc.n_claws[m].mean()))
            for k in list(A1_KINDS) + ['b', 'c']:
                r[f'{k}_mean'] = float(_kind_values(kc, k)[m].mean())
            trows.append(r)
    pd.DataFrame(trows).to_csv(out / 'a1_by_type_class.csv', index=False)
    # every wiring
    wrows = []
    for w in ALL_WIRINGS:
        kw, _ = load_kc_tables96(w)
        n = kw['n_claws'].to_numpy(dtype=float)
        b, c = kw['b'].to_numpy(dtype=float), kw['c'].to_numpy(dtype=float)
        for k in ('no_claw', 'peduncle', 'no_claw+peduncle', 'kept', 'b'):
            v = _kind_values(kw, k)
            wrows.append(dict(wiring=w, kind=k, r_n_claws=E95._corr(n, v)[0], loglog_slope_n=_loglog(n, v),
                              r_b=E95._corr(b, v)[0], r_c=E95._corr(c, v)[0], spearman_b=E95._corr(b, v)[1],
                              **_within(n, v, b)))
    pd.DataFrame(wrows).to_csv(out / 'a1_wirings.csv', index=False)
    # per claw count
    n = kc['n_claws'].to_numpy(dtype=float)
    b = kc['b'].to_numpy(dtype=float)
    crows = []
    for k in sorted(set(n.astype(int))):
        m = n == k
        r = dict(n_claws=int(k), n_kcs=int(m.sum()))
        for kind in ('no_claw', 'peduncle', 'no_claw+peduncle', 'kept', 'b'):
            v = _kind_values(kc, kind)[m]
            r[f'{kind}_mean'] = float(v.mean())
            r[f'{kind}_cv'] = float(v.std() / v.mean()) if v.mean() > 0 else np.nan
            if kind != 'b':
                r[f'{kind}_r_b'] = E95._corr(v, b[m])[0] if m.sum() > 2 else np.nan
        crows.append(r)
    pd.DataFrame(crows).to_csv(out / 'a1_within.csv', index=False)
    # the per-claw spread
    prow = []
    for col in [c for c in claws.columns if c.startswith('comp_')] + ['b_per_claw']:
        if col == 'b_per_claw':
            v = (kc.set_index(H.KC_ID_LEVEL)['b'] / kc.set_index(H.KC_ID_LEVEL)['n_claws']).reindex(
                claws[H.KC_ID_LEVEL]).to_numpy(dtype=float)
        else:
            v = claws[col].to_numpy(dtype=float)
        prow.append(dict(quantity=col, cv_claws=float(v.std() / v.mean()), zero_share=float((v == 0).mean())))
    pd.DataFrame(prow).to_csv(out / 'a1_per_claw.csv', index=False)
    return dict(counts=t, wirings=pd.DataFrame(wrows))


def a2_pratyush(kc_index: pd.MultiIndex) -> pd.DataFrame:
    """A2 in Pratyush's assignments: each APL->KC calyx synapse of the wiring's KCs (the loader's filters), its kind
    (kept / dropped_claw / no_claw, as the counts) and its distance (um) to the nearest postsynaptic site of a PN->KC
    synapse of the same KC on a claw with an ID of 0 or more (any claw, and a claw the wiring keeps)."""
    from scipy.spatial import cKDTree
    from pyolfsysm.mb_model import prat_hemibrain_seg_dir, APL_KC_SYNAPSE_TABLES
    t = pd.read_parquet(prat_hemibrain_seg_dir / APL_KC_SYNAPSE_TABLES['APL>KC'])
    claws = kc_index.to_frame(index=False)[[H.KC_ID_LEVEL, 'claw_id']]
    kcs = set(claws[H.KC_ID_LEVEL])
    keep = pd.MultiIndex.from_frame(claws.drop_duplicates())
    t = t[(t.instance_pre == 'APL_R') & (t.roi_pre == 'CA(R)') & (t.roi_post == 'CA(R)') & t.bodyId_post.isin(kcs)].copy()
    cl = t.anatomical_claw_corrected.to_numpy(dtype=float)
    has = ~np.isnan(cl)
    ci = np.where(has, cl, -2).astype(int)
    on = has & np.asarray(pd.MultiIndex.from_arrays([t.bodyId_post.to_numpy(), ci]).isin(keep))
    kind = np.where(on, 'kept', np.where(has & (ci != -1), 'dropped_claw', 'no_claw'))
    p = pd.read_parquet(prat_hemibrain_seg_dir / 'PNbouton-to-KCclaw_Connectivity_2025-12-09-19-47-36.parquet')
    p = p[(p.roi_post == 'CA(R)') & p.bodyId_post.isin(kcs) & (p.anatomical_claw_corrected >= 0)].copy()
    p['claw'] = p.anatomical_claw_corrected.astype(int)
    p['kept_claw'] = np.asarray(pd.MultiIndex.from_arrays([p.bodyId_post.to_numpy(), p.claw.to_numpy()]).isin(keep))
    um = 8e-3                                                       # hemibrain voxels of 8 nm
    d_any = np.full(len(t), np.nan)
    d_kept = np.full(len(t), np.nan)
    xyz = t[['x_post', 'y_post', 'z_post']].to_numpy(dtype=float) * um
    kc_of = t.bodyId_post.to_numpy()
    pg = {k: g for k, g in p.groupby('bodyId_post')}
    for k in np.unique(kc_of):
        idx = np.flatnonzero(kc_of == k)
        g = pg.get(k)
        if g is None:
            continue
        pts = g[['x_post', 'y_post', 'z_post']].to_numpy(dtype=float) * um
        d_any[idx] = cKDTree(pts).query(xyz[idx])[0]
        gk = g[g.kept_claw]
        if len(gk):
            d_kept[idx] = cKDTree(gk[['x_post', 'y_post', 'z_post']].to_numpy(dtype=float) * um).query(xyz[idx])[0]
    rows = []
    for k in ('kept', 'dropped_claw', 'no_claw', 'all'):
        m = np.ones(len(t), dtype=bool) if k == 'all' else kind == k
        for lab, dd in (('any_claw', d_any), ('kept_claw', d_kept)):
            v = dd[m]
            v = v[np.isfinite(v)]
            q = np.quantile(v, [0.1, 0.25, 0.5, 0.75, 0.9]) if v.size else [np.nan] * 5
            rows.append(dict(source='pratyush_hemibrain', kind=k, to=lab, n=int(m.sum()), n_with_distance=int(v.size),
                             d_p10=q[0], d_p25=q[1], d_p50=q[2], d_p75=q[3], d_p90=q[4],
                             within_1um=float((v <= 1).mean()) if v.size else np.nan,
                             within_2um=float((v <= 2).mean()) if v.size else np.nan,
                             within_3um=float((v <= 3).mean()) if v.size else np.nan))
    return pd.DataFrame(rows)


A2_DATASETS = ('hemibrain', 'flywire', 'malecns', 'banc')
A2_DIR = Path('/local/tensor/Tom/claude_APL_modeling_YZ4TC/Claude_APL_analysis/apl_synapses')


def a2_connectomes(base: Path = A2_DIR) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A2 in the derived tables of four connectomes (read only): per dataset and APL, the APL->KC synapses in the calyx
    (mb_region CA): the share within 2 um of a derived claw of the KC (claw_id given), the distance quantiles; per KC
    (the KCs with a derived claw), the synapses off the claws in the calyx and in the peduncle against the number of
    derived claws (r, log-log slope), the mean, the CV, the share of KCs with none, the share off the claws."""
    srows, krows = [], []
    for ds in A2_DATASETS:
        try:
            t = pd.read_csv(base / f'{ds}_APL_MB_synapses.csv.gz', low_memory=False,
                            usecols=lambda c: c in ('direction', 'apl_id', 'apl_instance', 'apl_side', 'partner_id',
                                                    'partner_class', 'partner_type', 'mb_region', 'claw_id',
                                                    'claw_dist_nm'))
            cl = pd.read_csv(base / f'{ds}_KC_claws_derived.csv.gz', usecols=['kc_id', 'claw_id'])
        except OSError as e:
            srows.append(dict(dataset=ds, note=repr(e)))
            continue
        t = t[t.direction.astype(str).str.startswith('APL_output') & (t.partner_class.astype(str).str.lower() == 'kenyon_cell')]
        n_claws = cl.groupby('kc_id').size()
        for apl, g in t.groupby('apl_id'):
            ca, ped = g[g.mb_region == 'CA'], g[g.mb_region == 'PED']
            dist = pd.to_numeric(ca.claw_dist_nm, errors='coerce').to_numpy(dtype=float) / 1000.0
            dist = dist[np.isfinite(dist)]
            q = np.quantile(dist, [0.25, 0.5, 0.75]) if dist.size else [np.nan] * 3
            on = ca.claw_id.notna()
            inst = g.apl_instance.iloc[0] if 'apl_instance' in g else np.nan
            side = g.apl_side.iloc[0] if 'apl_side' in g else np.nan
            srows.append(dict(dataset=ds, apl_id=apl, apl=str(inst) if isinstance(inst, str) else f'APL ({side})',
                              calyx_synapses=int(len(ca)), peduncle_synapses=int(len(ped)), share_within_2um=float(on.mean()),
                              d_p25=q[0], d_p50=q[1], d_p75=q[2]))
            kcs = pd.Index(sorted(set(ca.partner_id) & set(n_claws.index)))
            if not len(kcs):
                continue
            off = ca[~on].groupby('partner_id').size().reindex(kcs).fillna(0).to_numpy(dtype=float)
            onc = ca[on].groupby('partner_id').size().reindex(kcs).fillna(0).to_numpy(dtype=float)
            pc = ped.groupby('partner_id').size().reindex(kcs).fillna(0).to_numpy(dtype=float)
            n = n_claws.reindex(kcs).to_numpy(dtype=float)
            for lab, v in (('calyx_off_claws', off), ('calyx_on_claws', onc), ('peduncle', pc), ('off_claws_and_peduncle', off + pc)):
                krows.append(dict(dataset=ds, apl_id=apl, apl=srows[-1]['apl'], kind=lab, n_kcs=int(len(kcs)),
                                  synapses=float(v.sum()),
                                  per_kc_mean=float(v.mean()), cv_kcs=float(v.std() / v.mean()) if v.mean() > 0 else np.nan,
                                  kcs_none_share=float((v == 0).mean()), r_n_claws=E95._corr(n, v)[0],
                                  loglog_slope_n=_loglog(n, v), n_claws_mean=float(n.mean()),
                                  share_off_claws=float(off.sum() / (off.sum() + onc.sum())) if lab == 'calyx_off_claws' else np.nan,
                                  **{k: x for k, x in _within(n, v).items() if k.startswith('within_cv')}))
    return pd.DataFrame(srows), pd.DataFrame(krows)


def a2_distances(out: Path) -> Dict[str, pd.DataFrame]:
    """A2 (module docstring): a2_distances.csv (Pratyush's assignments) and a2_connectomes.csv / a2_connectomes_kc.csv
    (the four connectomes' derived tables)."""
    wn = E83.wiring_parts('real')['wnorm']['wAPLKC']
    a = a2_pratyush(wn.index)
    a.to_csv(out / 'a2_distances.csv', index=False)
    s, k = a2_connectomes()
    s.to_csv(out / 'a2_connectomes.csv', index=False)
    k.to_csv(out / 'a2_connectomes_kc.csv', index=False)
    return dict(pratyush=a, connectomes=s, connectomes_kc=k)


def a3_kc2apl(out: Path) -> pd.DataFrame:
    """A3 (module docstring): a3_kc2apl.csv, the KC->APL synapses of the same table split the same way (pyolfsysm's
    counts 'KC>APL' on the real wiring's claws), A1's columns, for the record."""
    from pyolfsysm.mb_model import connectome_APL_KC_synapse_counts
    wn = E83.wiring_parts('real')['wnorm']['wAPLKC']
    per_kc, _ = connectome_APL_KC_synapse_counts(wn.index, 'KC>APL')
    kc, _ = load_kc_tables96('real')
    k2 = per_kc.reindex(kc[H.KC_ID_LEVEL]).reset_index()
    for c in ('b', 'c', H.KC_TYPE_LEVEL):
        k2[c] = kc[c].to_numpy()
    t = pd.DataFrame([a1_row(k2, k) for k in A1_KINDS])
    calyx = float(k2[['kept', 'dropped_claw', 'no_claw']].to_numpy().sum())
    t['share_of_calyx'] = t.synapses / calyx
    t.insert(0, 'direction', 'KC>APL')
    t.to_csv(out / 'a3_kc2apl.csv', index=False)
    return t


def a4_shapes(out: Path) -> pd.DataFrame:
    """A4 (module docstring): a4_shapes.csv, per wiring and form, the per-KC value's log-log slope on the claw count, its
    CV over the KCs, the spread within a claw count, r with b_k, and the per-claw spread: the components (uniform,
    connectome, offclaw, offclaw_ped, allsyn; the per-KC sum of the per-claw shape) and E95's static forms (B1 b_k, B2
    beta n_k, B3 n_k^e, B4 b_k permuted within a claw count and over all KCs; per claw the per-KC value over n_k)."""
    rows = []
    for w in ALL_WIRINGS:
        kc, claws = load_kc_tables96(w)
        n = kc['n_claws'].to_numpy(dtype=float)
        b = kc['b'].to_numpy(dtype=float)
        forms = {f'comp_{m}': kc[f'comp_{m}_sum'].to_numpy(dtype=float) for m in ('uniform', 'connectome') + MODES}
        forms['B1_bk'] = b
        forms['B2_nk'] = E95.slope_through_origin(n, b) * n
        for e in E95.B3_EXPONENTS:
            x = n ** e
            forms[f'B3_nk{e:g}'] = x * (b.mean() / x.mean())
        forms['B4_bk_permclaw'] = E95.permute_values(b, E95.PERM_SEED, n)
        forms['B4_bk_permall'] = E95.permute_values(b, E95.PERM_SEED)
        for name, v in forms.items():
            pcl = (claws[name].to_numpy(dtype=float) if name.startswith('comp_')
                   else pd.Series(v / n, index=kc[H.KC_ID_LEVEL]).reindex(claws[H.KC_ID_LEVEL]).to_numpy(dtype=float))
            rows.append(dict(wiring=w, form=name, loglog_slope_n=_loglog(n, v), cv_kcs=float(v.std() / v.mean()),
                             r_n_claws=E95._corr(n, v)[0], r_b=E95._corr(b, v)[0], per_claw_cv=float(pcl.std() / pcl.mean()),
                             per_claw_zero_share=float((pcl == 0).mean()), **_within(n, v, b)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a4_shapes.csv', index=False)
    return df


def part_a(out: Path) -> dict:
    """A1 to A4 and a_summary.json."""
    t0 = time.time()
    a1 = a1_counts(out)
    a2 = a2_distances(out)
    a3 = a3_kc2apl(out)
    a4 = a4_shapes(out)
    t = a1['counts'].set_index('kind')
    summ = dict(a1={k: {c: (None if pd.isna(t.loc[k, c]) else float(t.loc[k, c])) for c in t.columns if c != 'name'}
                    for k in t.index},
                a1_wirings={k: {c: dict(min=float(g[c].min()), max=float(g[c].max()))
                                for c in ('r_n_claws', 'loglog_slope_n', 'r_b', 'r_c', 'within_r_b_median')}
                            for k, g in a1['wirings'].groupby('kind')},
                a2_pratyush=a2['pratyush'].to_dict('records'), a2_connectomes=a2['connectomes'].to_dict('records'),
                a3=a3.set_index('kind')[['synapses', 'per_kc_mean', 'r_n_claws']].to_dict('index'),
                a4_real=a4[a4.wiring == 'real'].set_index('form')[['loglog_slope_n', 'within_cv_median', 'per_claw_cv',
                                                                    'r_b']].to_dict('index'),
                seconds=time.time() - t0)
    _write_json(out / 'a_summary.json', summ)
    print(f'E96 A in {time.time() - t0:.0f} s', flush=True)
    return summ


# ---------------------------------------------------------------------------
# part B (and C's runs): the points

def wiring_wnorm(label: str) -> pd.Series:
    return E83.wiring_parts(label)['wnorm']['wAPLKC']


def point_weights(model: E95.StaticModel, pv: dict, unit: float) -> Tuple[Optional[pd.Series], Optional[pd.DataFrame]]:
    """A point's APL->KC weights (None: the fit's, size 0) and its component per KC (the outputs' KCs: n_claws, a_k as
    used, the component's sum over the KC's claws in unscaled weight units, the family's sum)."""
    w = model.label
    if not pv['size'] > 0:
        return None, None
    wn = wiring_wnorm(w)
    counts = wiring_counts(w) if pv['mode'] in MODES else None
    comp = component_values(pv['mode'], wn, counts, pv['permute'], pv['seed'] if pv['seed'] is not None else PERM_SEED)
    wk = e96_weights(model.wAPLKC, wn, comp, pv['size'], unit)
    kc_of = np.asarray(wn.index.get_level_values(H.KC_ID_LEVEL))
    ids = np.asarray(model.saved['on'][1].index.get_level_values(H.KC_ID_LEVEL))
    add = pd.Series(wk.to_numpy(dtype=float) - model.wAPLKC.to_numpy(dtype=float)).groupby(kc_of).sum().reindex(ids)
    frame = pd.DataFrame({H.KC_ID_LEVEL: ids, 'n_claws': pd.Series(kc_of).value_counts().reindex(ids).to_numpy(),
                          'component_sum': add.to_numpy(dtype=float),
                          'w_fam_sum': pd.Series(model.wAPLKC.to_numpy(dtype=float)).groupby(kc_of).sum().reindex(ids).to_numpy()})
    if pv['mode'] in MODES:
        a = offclaw_counts(pv['mode'], counts[0], pv['permute'], pv['seed'] if pv['seed'] is not None else PERM_SEED)
        frame['a_k'] = a.reindex(ids).to_numpy(dtype=float)
    return wk, frame


def run_point(model: E95.StaticModel, v: str, part: str = '', cc: Optional[dict] = None, sub: str = 'runs'
              ) -> Tuple[dict, dict]:
    """One variant on a wiring at its fit (`reading96`), saved with its stage matrices, drive, loss, APL->KC weights and
    component per KC under outputs/.../<sub>/. Returns (the run's row, the reading)."""
    w = model.label
    tag = run_tag(v, w)
    pv = parse_variant(v)
    if pv['part'] in (REF, 'e95'):
        raise ValueError(f'{v}: the saved fits and E95\'s points are not run here')
    cc = H.code_commits() if cc is None else cc
    kc_index = model.saved['on'][1].index
    L = E95.kc_losses(w, kc_index)
    unit = E95.d_unit(model, L)
    wk, comp_frame = point_weights(model, pv, unit)
    res = reading96(model, pv, wk, label=tag)
    o = {'on': res['on'], 'off': res['off']}
    st = res['stages']
    odors = E94._odors(st['bouton_win'])
    drive = E94.kc_drive(model.W, st['bouton_win'][odors].to_numpy(dtype=float), kc_index, odors)
    save_run96(tag, o, sub, stages=st, drive=drive, loss=_frame(res['gain']['L'], res['on'][1]), weights=wk,
               component=comp_frame)
    resp_off = E89._responders(model.saved['off'][1])
    p = res['params']
    per_kc = wiring_counts(w)[0]
    row = dict(part=part, arm=pv['arm'], tag=tag, variant=v, source='E96', screen=pv['screen'], mode=pv['mode'],
               size=pv['size'], pool=pv['pool'], dose=pv['dose'], permute=pv['permute'], seed=pv['seed'],
               control_of=pv['control_of'], wiring=w, host=_host(), fit_folder=model.folder, fit_tag=model.tag,
               factor=res['factor'], apl_factor=res['factor'], at_bound=res['at_bound'], bracket_lo=res['bracket'][0],
               bracket_hi=res['bracket'][1], n_on=E83._n(res['on'][0]), target_on=model.target_on,
               n_off=E83._n(res['off'][0]), target_off=model.target_off, n_sims=res['n_sims'], seconds=res['seconds'],
               **{k: p[k] for k in E83.PARAM_KEYS}, aplpn_fit_multiple=E94.aplpn_fit_multiple(w),
               aplpn_multiple=E94.aplpn_multiple(w, pv['dose']), wAPLKC_sha256=E89.E80_sha(model.wAPLKC if wk is None else wk),
               K_fit=model.K, K_run=E83.readout_K(o), component_unit=unit,
               synapse_share=synapse_share(pv['mode'], per_kc) if pv['mode'] in MODES else np.nan,
               **res['gain']['row'], **cc)
    row.update(E94.stage_extras(model.W, w, st))
    w_all = (model.wAPLKC if wk is None else wk).to_numpy(dtype=float) * p['wAPLKC_scale']
    w_fam = model.wAPLKC.to_numpy(dtype=float) * p['wAPLKC_scale']
    off_st = E94.fit_stage(w, 'off')
    row.update(E95.d_realism(model.W.to_numpy(dtype=float), off_st['bouton_win'][odors].to_numpy(dtype=float),
                             st['bouton_win'][odors].to_numpy(dtype=float), st['claw_win'][odors].to_numpy(dtype=float),
                             w_fam, w_all - w_fam, E94.kc_of_claw(model.W), resp_off))
    row.update(E83.realism(o))
    print(f"E96 {tag}: APL scale x {res['factor']:.4g}{' (at a bound)' if res['at_bound'] else ''}, {row['n_on']} / "
          f"{model.target_on} responses, {res['n_sims']} sims, {res['seconds']:.0f}s; shares at the boutons / component / "
          f"weighted {row.get('bouton_share_d', np.nan):.3f} / {row.get('comp_share', np.nan):.3f} / "
          f"{row.get('weighted_share', np.nan):.3f}; flag {row.get('realism_flag')!r}", flush=True)
    return row, res


RUNS_NAMES = {'B': 'runs_B.csv', 'C': 'runs_C.csv'}


def _provenance_add(out: Path, part: str, tag: str, prov: dict) -> None:
    path = out / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E96_offclaw_apl.py')
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def points(out: Path, variants: Sequence[str], wirings: Sequence[str], runs_name: str, part: str) -> pd.DataFrame:
    """Variants on wirings, each by `run_point` (resumable: rows in `runs_name`). Needs part 0 to have passed."""
    import gc
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identities must hold before any scored point)')
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
            row, _ = run_point(model, v, part, cc)
            _append(path, row)
            _provenance_add(out, part, run_tag(v, w), dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **cc))
        del model
        gc.collect()
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def part_b(out: Path, arms: Optional[Sequence[str]] = None, wirings: Sequence[str] = S_WIRINGS) -> pd.DataFrame:
    """The arms' points (default B1 to B6) on `wirings` (resumable; runs_B.csv)."""
    arms = list(ARMS) if not arms else list(arms)
    bad = [a for a in arms if a not in ARMS]
    if bad:
        raise ValueError(f'{bad} are not arms of part B ({list(ARMS)})')
    return points(out, [v for a in arms for v in ARMS[a]], wirings, RUNS_NAMES['B'], 'B')


# ---------------------------------------------------------------------------
# the references and the runs table

def e95_rows() -> pd.DataFrame:
    """E95's part D runs on S (runs_D.csv; source 'E95'): its uniform component (D1), the connectome's per-claw weights
    (D2) and D1 at d = 1 without the bouton path (D3), the comparison of part B."""
    p = H.results_dir(E95.EXPERIMENT) / 'runs_D.csv'
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p).drop_duplicates('tag', keep='last')
    d = d[d.variant.astype(str).str.startswith(E95_REF_ARMS)].copy()
    d['source'] = 'E95'
    return d


def runs_table(out: Path) -> pd.DataFrame:
    """Every run of the readings: this experiment's (runs_B / runs_C.csv, the latest row per tag) first, then E95's part D
    runs on S and the references (E95's `reference_rows`: the unchanged fits on every wiring, E94's D0 on S); with
    'flagged', 'at_bound' and the parsed setting."""
    parts = []
    for n in RUNS_NAMES.values():
        if (out / n).exists():
            try:
                parts.append(pd.read_csv(out / n).drop_duplicates('tag', keep='last'))
            except pd.errors.EmptyDataError:
                pass
    parts.append(e95_rows())
    parts.append(E95.reference_rows())
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates(['variant', 'wiring'], keep='first').reset_index(drop=True)
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('arm', 'part', 'screen', 'mode', 'size', 'pool', 'dose', 'permute', 'control', 'equal'):
        df[c] = pv[c].to_numpy()
    df['flagged'] = [E94._flagged(x) for x in df.realism_flag]
    df['at_bound'] = [E94._as_bool(x) for x in df.at_bound]
    return df


def load_any(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E96' this experiment's, 'E95' E95's, 'E94' E94's, 'saved' the wiring's fit."""
    if source == 'saved':
        return E83.saved_fit(wiring)['saved']
    o = {'E96': load_run96, 'E95': E95.load_run95, 'E94': E94.load_run}[source](tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


# ---------------------------------------------------------------------------
# the readings

READOUTS = E95.READOUTS


def _distances(source: str, tag: str, wiring: str, readouts: Sequence[str] = READOUTS) -> Dict[str, Dict[str, pd.Series]]:
    key = ('dist', source, tag, wiring)
    have = _CACHE.setdefault(key, {})
    need = [r for r in readouts if r not in have]
    if need:
        o = load_any(source, tag, wiring)
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


def readings96(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = READOUTS) -> pd.DataFrame:
    """readings.csv: per variant and unit (each wiring; S and T where all four have a run) and readout E92's
    `score_row92` against every target of E89's `targets89` with E94's group effects, the runs outside the bands and at a
    bracket's bound; pairs.csv: each population's per-pair distances and effect under 'cal' and 'raw'."""
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
                rows.append(dict(variant=variant, arm=pv['arm'], part=pv['part'], screen=pv['screen'], mode=pv['mode'],
                                 size=pv['size'], pool=pv['pool'], dose=pv['dose'], permute=pv['permute'],
                                 control=pv['control'], equal=pv['equal'], unit=unit, n_wirings=len(ws), readout=ro,
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


def _co_job96(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and the geometric mean of E88's two estimates at L1 and L3 (E95's
    `_co_job95` on this experiment's loader)."""
    import os
    import warnings
    source, tag, wiring, K, cand, rel, gap = job
    warnings.simplefilter('ignore')
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    o = load_any(source, tag, wiring)
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


def coresponse96(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse.csv: per S setting (all four S runs) the changes of L (co-response) and R (recruitment) under APL of the
    four nearest ester pairs, the ester and the alcohol groups on S's population, against the flies' (E95's
    `coresponse95` on this experiment's runs; jobs cached under outputs/.../coresponse_cache.pkl)."""
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
            for res, j in zip(ex.map(_co_job96, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E96 coresponse {j[1]}', flush=True)
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


# ---------------------------------------------------------------------------
# part C: the choice, the gate, the controls, the equal rule

PARTS_COLS = E92G.PARTS_COLS


def _settings96(rd: pd.DataFrame, screens: Sequence[str] = (REF,) + CHOICE_ARMS) -> pd.DataFrame:
    return rd[(rd.unit == 'S') & (rd.readout == 'cal') & rd.screen.isin(list(screens)) & ~rd.control.astype(bool)]


def _pat(rd: pd.DataFrame, variant: str, unit: str, target: str, col: str = 'pattern', readout: str = 'cal') -> float:
    return E89._pat(rd, variant, unit, target, col, readout)


def choose96(rd: pd.DataFrame, half: str) -> dict:
    """The choice on S against half `half`: the setting with the lowest pattern part of S's population under the
    calibrated readout among the unchanged fits and CHOICE_ARMS' settings whose four runs are inside the bands and not at
    a bracket's bound; something to confirm if it is a setting and gains CHOICE_MIN_GAIN or more over the unchanged fits.
    The parts of the pattern part beside (the thresholds stay on the pattern part)."""
    g = _settings96(rd)
    ok = g[(g.variant == REF) | ((g.n_flagged == 0) & (g.n_at_bound == 0))]
    if ok.empty or REF not in set(g.variant):
        return dict(half=half, chosen=None, confirm=False, note='runs missing')
    best = ok.loc[ok[f'half{half}_pattern'].idxmin()]
    rr = g[g.variant == REF].iloc[0]
    gain = float(rr[f'half{half}_pattern']) - float(best[f'half{half}_pattern'])
    confirm = bool(best.screen in CHOICE_ARMS and gain >= CHOICE_MIN_GAIN)
    row = dict(half=half, chosen=best.variant, screen=best.screen, mode=best['mode'], pattern=float(best[f'half{half}_pattern']),
               ref_pattern=float(rr[f'half{half}_pattern']), gain=gain, confirm=confirm, n_eligible=int(len(ok)),
               n_settings=int(len(g)), excluded=json.dumps(sorted(set(g.variant) - set(ok.variant))),
               control=control_of96(best.variant) if best.variant != REF else '')
    for c in PARTS_COLS[1:]:
        row[f'chosen_{c}'] = float(best[f'half{half}_{c}'])
        row[f'ref_{c}'] = float(rr[f'half{half}_{c}'])
    return row


def plan_c(rd: pd.DataFrame) -> Dict[str, dict]:
    return {h: choose96(rd, h) for h in HALVES}


def equal_rule96(rd: pd.DataFrame) -> List[dict]:
    """E94's equal rule at fixed scales: a setting of CHOICE_ARMS is an equal of the unchanged fits if, on S under the
    readout, its pattern part is within EQUAL_PATTERN of theirs against all recordings and halves A and B, its four runs
    are inside the bands and none at a bracket's bound, its raw-count loss within EQUAL_RAW and its mean APL effect within
    EQUAL_MEAN of theirs; 'better' if its pattern part is lower on all three (and the rest holds)."""
    g = _settings96(rd)
    out = []
    if REF not in set(g.variant):
        return out
    me0 = float(g[g.variant == REF].mean_effect.iloc[0])
    for r in g.itertuples():
        if r.screen not in CHOICE_ARMS:
            continue
        d = {t: _pat(rd, r.variant, 'S', t) - _pat(rd, REF, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = _pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - _pat(rd, REF, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - me0
        rest = bool(int(r.n_flagged) == 0 and int(r.n_at_bound) == 0 and abs(raw) <= EQUAL_RAW and abs(me) <= EQUAL_MEAN)
        out.append(dict(variant=r.variant, arm=r.arm, screen=r.screen, mode=r.mode, size=float(r.size),
                        **{f'd_pattern_{t}': x for t, x in d.items()}, d_raw_loss=raw, d_mean_effect=me,
                        n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound),
                        equal=bool(all(abs(x) <= EQUAL_PATTERN for x in d.values()) and rest),
                        better=bool(all(x < 0 for x in d.values()) and rest)))
    return out


def vs_uniform(rd: pd.DataFrame) -> List[dict]:
    """Each setting of B1, B2, B3, B5, B6 against E95's uniform component at the same size, pooled share and dose (its
    rule (v) control; on S under the readout against all recordings and halves A and B, the mean APL effect), and
    against the same setting's counts permuted (B4, for B2 at d = 1)."""
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    have = set(g.variant)
    out = []
    for r in g[g.screen.isin(CHOICE_ARMS) & ~g.control.astype(bool)].itertuples():
        u = control_of96(r.variant)
        row = dict(variant=r.variant, uniform=u, uniform_read=u in have)
        if u in have:
            for t in ('z8', 'halfA', 'halfB'):
                row[f'd_pattern_{t}'] = _pat(rd, r.variant, 'S', t) - _pat(rd, u, 'S', t)
            row['d_mean_effect'] = float(r.mean_effect) - float(g[g.variant == u].mean_effect.iloc[0])
        out.append(row)
    return out


def part_c(out: Path) -> dict:
    """C: the choice per half (choice.csv); each chosen setting with something to confirm run on T (runs_C.csv), then its
    controls on T (`control_of96`, run on S too where it is not one of E95's; `beside_controls`); then `score` without the
    co-response."""
    runs = runs_table(out)
    rd = readings96(out, runs, readouts=('cal', 'raw'))
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    chosen = sorted({c['chosen'] for c in plan.values() if c['confirm']})
    if not chosen:
        print('E96 C: nothing to confirm', flush=True)
    for v in chosen:
        points(out, [v], T_WIRINGS, RUNS_NAMES['C'], 'C')
    for v in chosen:
        u = control_of96(v)
        points(out, [u], (S_WIRINGS + T_WIRINGS) if u.startswith('U_') else T_WIRINGS, RUNS_NAMES['C'], 'C')
        points(out, beside_controls(v), T_WIRINGS, RUNS_NAMES['C'], 'C')
    return score(out, coresponse=False)


def gate_rows(rd: pd.DataFrame, runs: pd.DataFrame, plan: Dict[str, dict]) -> List[dict]:
    """The gate's four directions (E92's `gate_direction92` with the rule (v) control, the uniform component), each with
    the beside controls' gains on T against the other half ('beside_gains', JSON)."""
    dirs = []
    for X in HALVES:
        c = plan[X]
        chosen = c['chosen'] if c['confirm'] else None
        r = E92G.gate_direction92(rd, runs, X, chosen, REF, [control_of96(chosen)] if chosen else [])
        if chosen:
            Y = E89.HALF_OTHER[X]
            r0 = _pat(rd, REF, 'T', f'half{Y}')
            bg = {}
            for b in beside_controls(chosen):
                if ((rd.variant == b) & (rd.unit == 'T')).any():
                    bg[b] = float(r0 - _pat(rd, b, 'T', f'half{Y}'))
            r['beside_gains'] = json.dumps(bg)
            r['rule_v_control'] = control_of96(chosen)
        dirs.append(r)
    return dirs


# ---------------------------------------------------------------------------
# part D's plan

def _best(rd: pd.DataFrame, screens: Sequence[str], variants: Optional[Sequence[str]] = None) -> Optional[pd.Series]:
    g = _settings96(rd, screens)
    if variants is not None:
        g = g[g.variant.isin(list(variants))]
    if g.empty:
        return None
    return g.loc[g.z8_pattern.idxmin()]


def d_plan(rd: pd.DataFrame, eqs: Optional[Sequence[dict]] = None) -> dict:
    """Part D's plan (the entry's rules): it runs if S's best setting of B1 or B2 (the lowest pattern part of S's
    population under the readout against all recordings) is within D_COND of the unchanged fits' UNCHANGED_S['all'];
    the mode is the one of B1 and B2 with the lower pattern part at its best size (B1 when within D_MODE_TIE); the pooled
    share is 1 only if B3 is lower than B2 at that size by D_POOL_GAIN or more; B6 (allsyn) is fitted only if a setting of
    B6 is an equal at fixed scales (`equal_rule96`)."""
    out = dict(cond=D_COND, unchanged=UNCHANGED_S['all'],
               unchanged_here=_pat(rd, REF, 'S', 'z8') if REF in set(rd.variant) else np.nan)
    b1, b2 = _best(rd, ('B1',)), _best(rd, ('B2',))
    if b1 is None or b2 is None:
        return dict(out, runs=False, note='B1 / B2 not run on S')
    best = b1 if float(b1.z8_pattern) <= float(b2.z8_pattern) else b2
    out.update(best=best.variant, best_pattern=float(best.z8_pattern), above=float(best.z8_pattern) - UNCHANGED_S['all'],
               b1_best=b1.variant, b1_pattern=float(b1.z8_pattern), b2_best=b2.variant, b2_pattern=float(b2.z8_pattern))
    out['runs'] = bool(out['above'] <= D_COND)
    tie = abs(float(b1.z8_pattern) - float(b2.z8_pattern)) <= D_MODE_TIE
    sel = b1 if (tie or float(b1.z8_pattern) < float(b2.z8_pattern)) else b2
    out.update(mode=sel['mode'], mode_arm=sel.arm, mode_size=float(sel['size']), mode_tie=bool(tie))
    d = float(sel['size'])
    v2, v3 = _vname('B2', 'offclaw_ped', d), _vname('B3', 'offclaw_ped', d, pool=True)
    have = set(_settings96(rd, ('B2', 'B3')).variant)
    if v2 in have and v3 in have:
        diff = _pat(rd, v3, 'S', 'z8') - _pat(rd, v2, 'S', 'z8')
        out.update(pool_compare=[v3, v2], pool_diff=float(diff), pool=POOL if diff <= -D_POOL_GAIN else None)
    else:
        out.update(pool=None, pool_note=f'{v3} or {v2} not read')
    eqs = equal_rule96(rd) if eqs is None else eqs
    b6 = [e['variant'] for e in eqs if e['arm'] == 'B6' and e['equal']]
    out.update(allsyn=bool(b6), b6_equals=b6)
    return out


def d_arm(mode: str, pool: Optional[float]) -> str:
    """The arm of part B whose points start a fit of `mode` (with the pooled share)."""
    if mode == 'allsyn':
        return 'B6'
    if pool:
        return 'B3'
    return 'B1' if mode == 'offclaw' else 'B2'


def d_tag(label: str, mode: str, pool: Optional[float], start: str) -> str:
    """The tag E10_joint_fit.py gives a part D fit (its prefix rules: the wiring, rate0.068_, aplboth_, pool1_,
    aplpnhold_, flato_ / flatp_ / flata_, then k2 and D_TAGS[start])."""
    from experiments import E10_joint_fit as J
    a = E91.d_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return (f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.claw_inh_pool_tag(pool)}{J.aplpn_hold_tag(1.0)}"
            f"{J.aplkc_flat_tag(dict(kind=mode))}k2{D_TAGS[start]}")


def d_start_key(label: str, mode: str, pool: Optional[float], start: str) -> str:
    return f"{label}_{MODE_CODE[mode]}{'_pool1' if pool else ''}_{start}"


def own_start(rd: pd.DataFrame, runs: pd.DataFrame, label: str, mode: str, pool: Optional[float]) -> dict:
    """The 'own' start of a wiring: the absolute parameters of its best point of the fitted arm of part B (the lowest
    pattern part of that wiring under the readout against all recordings, inside the bands first, none at a bracket's
    bound), with that point's size and unit."""
    arm = d_arm(mode, pool)
    q = runs[(runs.screen == arm) & (runs.wiring == label) & ~runs.at_bound]
    if q.empty:
        raise RuntimeError(f'{label}: no {arm} point inside the bracket (run --part B --arm {arm} first)')
    sc = rd[(rd.unit == label) & (rd.readout == 'cal')].set_index('variant').z8_pattern
    q = q.assign(pattern=q.variant.map(sc), flag_order=q.flagged.astype(int)).sort_values(['flag_order', 'pattern'])
    r = q.iloc[0]
    return dict({k: float(r[k]) for k in E83.PARAM_KEYS}, aplkc_flat=float(r['size']), from_tag=r.tag,
                pattern=float(r.pattern), flagged=bool(r.flagged), component_unit=float(r['component_unit']), start='own')


def e95_start(label: str, unit: float) -> dict:
    """The 'e95' start of a wiring: E95's refit of that wiring with the uniform component (its fits.csv row): its
    absolute parameters and fitted size, the component exchanged; the unit is this experiment's for the wiring (`unit`,
    checked against E95's)."""
    tag = E95.f_tag(label, 'uniform')
    f = pd.read_csv(H.results_dir(E95.EXPERIMENT) / 'fits.csv').drop_duplicates('tag', keep='last').set_index('tag')
    if tag not in f.index:
        raise RuntimeError(f'{label}: E95 has no refit {tag}')
    r = f.loc[tag]
    u95 = float(r['aplkc_flat_unit'])
    if abs(u95 - unit) > 1e-9 * unit:
        raise RuntimeError(f"{label}: E95's unit {u95} is not this experiment's {unit}")
    p = {k: float(r[k]) for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')}
    p['wKCKC_scale'] = float(r['kc_kc_inh'])
    return dict(p, aplkc_flat=float(r['aplkc_flat_size']), from_tag=f'{E95.EXPERIMENT}:{tag}', component_unit=float(unit),
                unit_e95=u95, start='e95')


def t_start(label: str, size: float) -> dict:
    """T's start: the wiring's own fit with the APL->bouton scale at dose D_DOSE and the component at `size` (the median
    of S's fitted sizes); the unit the wiring's."""
    m = E95.StaticModel(label, check=False)
    p = m.params_on(1.0, D_DOSE)
    unit = E95.d_unit(m)
    return dict({k: float(p[k]) for k in E83.PARAM_KEYS}, aplkc_flat=float(size), from_tag=f'{m.folder}:{m.tag}',
                component_unit=float(unit), start='own')


def write_start(out: Path, key: str, start: dict) -> None:
    path = out / D_START_FILE
    d = _read_json(path)
    d[key] = start
    _write_json(path, d)


def d_command(label: str, mode: str, pool: Optional[float], start: str, unit: float, maxfev: int = D_MAXFEV) -> List[str]:
    """One part D fit's E10_joint_fit.py command: C2's family (`-k 2 --rate-target ctrl_frac --loss apl --apl-readout
    both`) with the APL->bouton scale held at D_DOSE times the wiring's fitted multiple of E3's base (`--aplpn-hold`, as
    E94's D1), the component's size fitted in D_FLAT_BOUNDS (`--aplkc-flat <mode>`, its unit `unit`), the pooled share
    with `pool`, into this experiment's folders, from d_start_params.json's `d_start_key` (`--init-params`, which also
    gives the size's start)."""
    m = E91.d_member(label)
    F = E94.aplpn_multiple(label, D_DOSE)
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
           repr(float(F)), '--aplkc-flat', mode, '--aplkc-flat-unit', repr(float(unit)), '--aplkc-flat-max',
           f'{D_FLAT_BOUNDS[1]:g}', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag',
           D_TAGS[start], '--maxfev', str(int(maxfev)), '--experiment', EXPERIMENT,
           '--init-params', f'results/{EXPERIMENT}/{D_START_FILE}:{d_start_key(label, mode, pool, start)}']
    if pool:
        cmd += ['--claw-inh-pool', f'{float(pool):g}']
    return cmd + m['wiring_args']


def _fits(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def d_allowed(out: Path, pop: str, mode: Optional[str] = None) -> Tuple[bool, dict]:
    """S's fits: on part D's condition (`d_plan`), the mode its rule (or B6's allsyn on its own rule); T's: only if S's
    fitted population is an equal or better (`d_scores`)."""
    rd = pd.read_csv(out / 'readings.csv') if (out / 'readings.csv').exists() else None
    if rd is None:
        return False, dict(note='not scored (run --part score)')
    plan = d_plan(rd)
    if pop == 'S':
        if mode == 'allsyn':
            return bool(plan.get('runs')) and bool(plan.get('allsyn')), plan
        return bool(plan.get('runs')), plan
    m_ = mode or plan.get('mode')
    s = d_scores(out, write=False).get(score_key('S', m_, None if m_ == 'allsyn' else plan.get('pool')))
    if not s:
        return False, dict(plan, note="S's fits not scored")
    return bool(s.get('equal') or s.get('better')), dict(plan, s_scores=s)


def part_d(out: Path, pop: str = 'S', start: str = 'own', wirings: Optional[Sequence[str]] = None,
           mode: Optional[str] = None, force: bool = False, dry_run: bool = False) -> List[str]:
    """Part D's fits of population `pop` from `start` (or `wirings` of it: a lane), one process after another
    (`d_command`; logs under outputs/.../D/), skipping those in fits.csv; refuses what the entry's rules do not allow
    unless `force`. `mode` None: the plan's; 'allsyn': B6's fits on their rule."""
    import subprocess
    ok, plan = d_allowed(out, pop, mode)
    mode = mode or plan.get('mode')
    pool = None if mode == 'allsyn' else plan.get('pool')
    print(f"E96 D {pop} {start}: {'allowed' if ok else 'not allowed'}; mode {mode}, pool {pool}; "
          f"{json.dumps({k: plan.get(k) for k in ('best', 'best_pattern', 'above', 'runs', 'mode_arm', 'mode_size', 'pool_diff', 'allsyn')}, default=str)}",
          flush=True)
    if (not ok and not force) or mode is None:
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
        s = d_scores(out, write=False).get(score_key('S', mode, pool), {})
        sizes = s.get('sizes')
        if not sizes:
            raise RuntimeError("T's fits need S's fitted sizes (run --part Dscore)")
    for w in ws:
        tag = d_tag(w, mode, pool, start)
        tags.append(tag)
        if tag in set(_fits(out).tag):
            print(f'E96 D {tag}: done', flush=True)
            continue
        if pop == 'T':
            st = t_start(w, float(np.median(sizes)))
        elif start == 'own':
            st = own_start(rd, runs, w, mode, pool)
        else:
            unit = float(runs[(runs.wiring == w) & (runs.source == 'E96')].component_unit.dropna().iloc[0])
            st = e95_start(w, unit)
        key = d_start_key(w, mode, pool, start)
        cmd = d_command(w, mode, pool, start, st['component_unit'])
        print(f'E96 D {key}: start {json.dumps(st, default=str)}', flush=True)
        print('E96 D:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        write_start(out, key, st)
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E96 D {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits(out).tag):
            raise RuntimeError(f'E96 D: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def score_key(unit: str, mode: str, pool: Optional[float]) -> str:
    """d_scores.json's key of a population (S, T, all8) of fits of `mode` with the pooled share `pool`."""
    return f"{unit}_{mode}{'_pool1' if pool else ''}"


def load_fit96(tag: str) -> Optional[Dict[str, tuple]]:
    """A part D fit's outputs: this host's outputs, the archive, the archive's from_hal/."""
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    for d in (H.outputs_dir(EXPERIMENT), arc, arc / 'from_hal'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def pick_start(rows: Sequence[dict]) -> Optional[dict]:
    """Of a wiring's fits from the two starts, the one with the lower final value of its own loss ('fit_loss'); None
    when none was fitted."""
    rows = [r for r in rows if r is not None and np.isfinite(float(r.get('fit_loss', np.nan)))]
    if not rows:
        return None
    return min(rows, key=lambda r: float(r['fit_loss']))


def d_scores(out: Path, write: bool = True) -> dict:
    """Part D's scores: per wiring, start and mode the readings under the readout and on raw counts against the unchanged
    fits of the same wirings, with realism (d_fits.csv); per wiring the fit with the lower final value of its own loss
    (`pick_start`); S by E94's D1 rule (`E94.d1_rule`), beside E95's refit (as good as the uniform component: within
    D_AS_GOOD of it, or lower) and E94's; the eight against UNCHANGED_8 by the same rule and E91's part D pass rule.
    d_scores.csv, d_scores.json."""
    fits = _fits(out)
    if fits.empty:
        return {}
    fits = fits.drop_duplicates('tag', keep='last')
    rows, res = [], {}
    u95 = e95_refit_s()
    for mode in MODES:
        for pool in (None, POOL):
            st = {}
            for w in ALL_WIRINGS:
                cand = []
                for start in D_TAGS:
                    try:
                        tag = d_tag(w, mode, pool, start)
                    except Exception:
                        continue
                    if tag not in set(fits.tag):
                        continue
                    o = load_fit96(tag)
                    if o is None:
                        continue
                    d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
                    d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
                    sc, sc0 = E92G.score_row92(d['cal']), E92G.score_row92(d0['cal'])
                    raw, raw0 = E83.effect_scores(d['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
                    fr = fits[fits.tag == tag].iloc[-1]
                    row = dict(mode=mode, pool=pool or 0.0, unit=w, start=start, tag=tag, **E94._d_row(sc, sc0, raw, raw0),
                               lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                               **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale',
                                                         'wPNAPL_scale', 'kc_kc_inh', 'aplpn_hold', 'aplpn_factor_vs_e3',
                                                         'aplkc_flat_size', 'aplkc_flat_unit', 'init_params', 'n_evals',
                                                         'initial_loss')},
                               fit_loss=float(fr.get('loss')), **E83.realism(o))
                    cand.append((row, d, d0))
                if not cand:
                    continue
                pick = pick_start([c[0] for c in cand])
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
                    r['vs_e94_refit'] = float(r['pattern_all'] - E94_REFIT_S)
                if unit == 'all8':
                    rest = bool(r['n_flagged'] == 0 and abs(r['raw_loss'] - r['raw_loss_unchanged']) <= E94.D_RAW_SLACK
                                and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
                    r['equal'] = bool(abs(r['pattern_all'] - UNCHANGED_8) <= E94.D_EQUAL and rest)
                    r['passes'] = bool(r['pattern_all'] <= UNCHANGED_8 - E94.D_PASS_GAIN and r['n_lower'] >= E94.D_WIRINGS_LOWER
                                       and r['n_flagged'] == 0 and r['raw_loss'] <= E94.D_RAW_REF + E94.D_RAW_SLACK
                                       and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
                res[score_key(unit, mode, pool)] = r
                rows.append(dict(mode=mode, pool=pool or 0.0, unit=unit, **{k: x for k, x in r.items() if k not in ('sizes', 'starts')}))
    if write:
        pd.DataFrame(rows).to_csv(out / 'd_scores.csv', index=False)
        _write_json(out / 'd_scores.json', res)
    return res


def d_check(out: Path, modes: Sequence[str] = MODES) -> dict:
    """E10_joint_fit.py with --aplkc-flat <mode> at size 0 against a saved fit's first evaluation: E94's D1 fit on the
    real wiring (its command, from its own start) run with --smoke (5 evaluations, into this experiment's results folder,
    tag suffix _e10check), each mode at size 0 with this wiring's unit; each first evaluation's loss terms must equal
    E94's history row 0 exactly. d_check.json."""
    import subprocess
    hist94 = pd.read_csv(H.results_dir(E94.EXPERIMENT) / f"history_{E94.d1_tag('real')}.csv")
    start = _read_json(H.results_dir(E94.EXPERIMENT) / E94.D_START_FILE)[E94.d_start_key('real')]
    m = E91.d_member('real')
    base = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
            repr(float(E94.aplpn_multiple('real', float(start['dose'])))), '--rate-target', 'ctrl_frac', '--loss', 'apl',
            '--apl-readout', 'both', '--maxfev', '5', '--smoke', '--experiment', EXPERIMENT,
            '--init-params', f'results/{E94.EXPERIMENT}/{E94.D_START_FILE}:{E94.d_start_key("real")}'] + m['wiring_args']
    L = E95.kc_losses('real')
    wn = wiring_wnorm('real')
    unit = E95.component_unit(L['bclaw_mean'], E83.saved_fit('real')['p']['wAPLKC_scale'], L['apl_mean'], float(wn.mean()))
    from experiments import E10_joint_fit as J
    res = dict(reference=f"{E94.EXPERIMENT}: history_{E94.d1_tag('real')}.csv row 0", unit=unit)
    cols = [c for c in hist94.columns if c not in ('time_s', 'eval', 'Unnamed: 0')]
    for mode in modes:
        cmd = base + ['--tag', '_e10check', '--aplkc-flat', mode, '--aplkc-flat-unit', repr(unit), '--aplkc-flat-init', '0']
        log = local_dir('D') / f'd_check_{mode}.log'
        with open(log, 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E96 d_check {mode}: E10_joint_fit.py exited {r.returncode} (log {log})')
        tag = E94.d1_tag('real').replace(E94.D1_TAG, '_e10check').replace(
            'aplpnhold_', 'aplpnhold_' + J.aplkc_flat_tag(dict(kind=mode)))
        h = pd.read_csv(out / f'history_{tag}.csv')
        diff = {c: [float(h[c].iloc[0]), float(hist94[c].iloc[0])] for c in cols if c in h.columns
                and not E95.E10_same(h[c].iloc[0], hist94[c].iloc[0])}
        res[mode] = dict(tag=tag, n_cols=int(sum(c in h.columns for c in cols)), differs=diff, exact=not diff,
                         flat_first=float(h['aplkc_flat'].iloc[0]) if 'aplkc_flat' in h else None,
                         flat_later=[float(x) for x in h['aplkc_flat'].iloc[1:]] if 'aplkc_flat' in h else None)
        print(f'E96 d_check {mode}: {json.dumps(res[mode], default=str)[:600]}', flush=True)
    res['passes'] = bool(all(res[m_]['exact'] and res[m_]['flat_first'] == 0.0 and res[m_]['n_cols'] > 0 for m_ in modes))
    _write_json(out / 'd_check.json', res)
    return res


def part_e(out: Path) -> None:
    """Part E (only by the entry's rules): on an equal as a fitted family, the stage readings of the fits (the shares of
    the claws' loss made at the boutons, by the component and by the weighted part, beside the synapses' shares), the
    step from the claws to the somata and the alcohol pairs' co-response change; on a pass by (C) or as a fitted family
    of eight, the firing-ceiling family with the response rate with APL on re-held, S at the flies' 2 s pulse with both
    response rates re-held, and S through the trial noise as E72's part E ran it (each holds if the population keeps at
    least half of its gain). Not written: it runs only by those rules."""
    v = _read_json(out / 'verdict.json')
    ds = _read_json(out / 'd_scores.json')
    due = dict(equal_fitted=any(bool(x.get('equal') or x.get('better')) for k, x in ds.items() if k.startswith('S_')),
               passes=bool(v.get('verdict') == 'passes' or any(x.get('passes') for k, x in ds.items() if k.startswith('all8_'))))
    print(f'E96 E: {due}', flush=True)
    raise NotImplementedError("E96 part E runs only by the entry's rules; it is not written")


# ---------------------------------------------------------------------------
# the verdict and score

def verdict96(out: Path, rd: Optional[pd.DataFrame] = None, runs: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: the choices, the gate in four directions, the equals, each setting against E95's uniform component,
    S's settings with their parts, part D's plan and scores, the flagged runs and those at a bracket's bound, the alcohol
    pairs' co-response change."""
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table(out) if runs is None else runs
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    dirs = gate_rows(rd, runs, plan)
    pd.DataFrame(dirs).to_csv(out / 'gate.csv', index=False)
    gate = E92G.gate92(dirs)
    eqs = equal_rule96(rd)
    pd.DataFrame(eqs).to_csv(out / 'equal.csv', index=False)
    vu = vs_uniform(rd)
    pd.DataFrame(vu).to_csv(out / 'vs_uniform.csv', index=False)
    gains_s = any(c['confirm'] for c in plan.values())
    v = dict(thresholds=dict(CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=E91.GATE_GAIN, GATE_SINGLES_SHARE=E91.GATE_SINGLES_SHARE,
                             GATE_RAW_SLACK=E91.GATE_RAW_SLACK, GATE_MEAN_SLACK=E91.GATE_MEAN_SLACK, GATE_SQ_ALL=E91.GATE_SQ_ALL,
                             GATE_SQ_SIX=E91.GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO, EQUAL_PATTERN=EQUAL_PATTERN,
                             EQUAL_RAW=EQUAL_RAW, EQUAL_MEAN=EQUAL_MEAN, APL_BRACKET=APL_BRACKET, D_COND=D_COND,
                             D_MODE_TIE=D_MODE_TIE, D_POOL_GAIN=D_POOL_GAIN, D_AS_GOOD=D_AS_GOOD),
             readout="cal (the calibrated readout, K set on each run's APL-on counts; E83's readings and gate)",
             reference=REF, choice=plan, gate=gate,
             verdict='passes' if gate['passes'] else ('gains on S, not through the gate' if gains_s else 'level'),
             equals=eqs, equal_settings=[e['variant'] for e in eqs if e['equal']], vs_uniform=vu)
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    v['S'] = {r.variant: dict(arm=r.arm, screen=r.screen, halfA=r.halfA_pattern, halfB=r.halfB_pattern, halfC=r.halfC_pattern,
                              halfD=r.halfD_pattern, all=r.z8_pattern, six=r.z6_pattern, loss=r.z8_loss, mean_effect=r.mean_effect,
                              amplitude=r.z8_amplitude, shape=r.z8_shape, sd_ratio=r.z8_sd_ratio, pearson=r.z8_pearson,
                              rho=r.z8_rho, err_ester=r.z8_err_ester, err_alcohol=r.z8_err_alcohol, n_flagged=int(r.n_flagged),
                              n_at_bound=int(r.n_at_bound)) for r in g.itertuples()}
    v['D_plan'] = d_plan(rd, eqs)
    if (out / 'd_scores.json').exists():
        v['D'] = _read_json(out / 'd_scores.json')
    mine = runs[runs.source == 'E96']
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
    print(f"E96 verdict: {v['verdict']}; " + '; '.join(f"{h} {plan[h].get('chosen')} (gain {plan[h].get('gain', np.nan):.4f})"
                                                     for h in HALVES)
          + f"; equals {v['equal_settings']}; D {'runs' if v['D_plan'].get('runs') else 'does not run'} "
            f"(mode {v['D_plan'].get('mode')}, pool {v['D_plan'].get('pool')}, allsyn {v['D_plan'].get('allsyn')})", flush=True)
    return v


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """readings.csv, pairs.csv, line.csv, coresponse.csv (with `coresponse`), choice.csv, gate.csv, equal.csv,
    vs_uniform.csv, verdict.json, the figure and comparison_table.md."""
    runs = runs_table(out)
    rd = readings96(out, runs, readouts=READOUTS if zcell else tuple(r for r in READOUTS if not r.startswith('zcell')))
    try:
        E95.line95(out)
    except Exception as e:                       # (E92 part A's line needs its line.json; reported, not fatal)
        print(f'E96 score: line.csv left out ({e!r})', flush=True)
    if coresponse:
        coresponse96(out, n_workers)
    v = verdict96(out, rd, runs)
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
    beside C2's family on S (E95's `part_steps` on this experiment's runs): somata_steps.csv, claws_steps.csv, steps.csv."""
    from experiments import E90_model_stages as E90
    from experiments import E92_steps_like_for_like as A
    variants = list(variants) if variants else steps_variants(out)
    if not variants:
        print('E96 steps: no chosen or equal setting', flush=True)
        return None
    t0 = time.time()
    runs = runs_table(out)
    soma, claw = [], []
    i = 0
    loaders = {'E96': load_stage96, 'E95': E95.load_stage95, 'E94': E94.load_stage94}
    for v in variants:
        for w in S_WIRINGS:
            q = runs[(runs.variant == v) & (runs.wiring == w)]
            if not len(q):
                raise RuntimeError(f'{v} has no run on {w}')
            r = q.iloc[0]
            o = load_any(r.source, r.tag, w)
            soma.append(A.somata_pairs(o['on'][1], o['off'][1]).assign(family=v, wiring=w, source=r.source, tag=r.tag))
            st_on = loaders[r.source](r.tag) if r.source in loaders else dict(E94.fit_stage(w, 'on'))
            if st_on is None:
                raise RuntimeError(f'{r.tag}: no stage matrices')
            if 'spike_counts' not in st_on:
                st_on['responses'], st_on['spike_counts'] = o['on']
            W, _ = E94.wiring_W(w)
            meta = {k: np.asarray(x) for k, x in E92G.units_meta(W).items()}
            rr = E90.read_wiring(w, st_on, dict(E94.fit_stage(w, 'off')), meta,
                                 np.random.default_rng(E90.NOISE_SEED + E94.STEPS_NOISE_SEED + 200 + i))
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
    print(f'E96 steps in {time.time() - t0:.0f} s', flush=True)
    return rows


# ---------------------------------------------------------------------------
# the archive, the figure and the table

MODE_NAMES = {'offclaw': "each KC's APL->KC synapses on no claw in the calyx (a^c_k) spread evenly over its claws",
              'offclaw_ped': "each KC's APL->KC synapses on no claw in the calyx and in the peduncle (a^c_k + a^p_k) spread "
                             'evenly over its claws',
              'allsyn': "the synapses on each kept claw plus the KC's synapses on no claw (calyx and peduncle) spread evenly "
                        'over its claws, in the proportion in which they are counted',
              'uniform': "E95's uniform component (the connectome's mean weight on every claw)"}


def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike counts,
    the harness's primary model RDM): one row per run of this experiment with its absolute parameters, the component,
    its size and unit, the pooled share, the dose, realism; part D's fits' distances, which E10_joint_fit.py wrote, are
    kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E96']
    rows, parts = [], []
    for _, r in runs.iterrows():
        pv = parse_variant(r.variant)
        kind = ("C2's family, the tied APL scale re-held; the APL->KC weights the family's"
                + (f" plus a component of size {pv['size']:g} (unit {r.get('component_unit'):.6g}): {MODE_NAMES[pv['mode']]}"
                   if pv['size'] > 0 else '')
                + (f"; its counts permuted {'within a claw count' if pv['permute'] == 'claw' else 'over all KCs'} (seed {pv['seed']})"
                   if pv['permute'] else '')
                + f"; APL->bouton scale x {pv['dose']:g} of the fit's"
                + (f"; pooled share {pv['pool']:g} (olfsysm claw_inh_pool)" if pv['pool'] else '')
                + "; the APL-off run the fit's")
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=r.arm, kind=kind, config=r.variant, mode=r['mode'],
                         size=r['size'], pool=r.pool, dose=r.dose, permute=r.permute, seed=r.seed, wiring=r.wiring,
                         ref_experiment=r.fit_folder, ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, factor=r.factor, aplpn_multiple=r.aplpn_multiple,
                         component_unit=r.get('component_unit'), synapse_share=r.get('synapse_share'), rate_on=r.rate_on,
                         rate_off=r.rate_off, max_count_on=r.max_count_on, max_count_off=r.max_count_off,
                         wAPLKC_sha256=r.wAPLKC_sha256, control_of=r.control_of, bouton_share=r.get('bouton_share'),
                         comp_share=r.get('comp_share'), realism_flag=r.realism_flag, host=r.get('host'), subdir='runs'))
        o = load_run96(r.tag)
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
    """E96_settings.png: S's pattern part under the readout against all recordings and halves A and B per setting of part
    B and E95's part D settings (the unchanged fits dashed; hollow: a setting with a run outside the bands or at a
    bracket's bound)."""
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
    groups = [('part B', g[g.part == 'B']), ("E95's part D and D0", g[g.part.isin(['e95', 'control'])])]
    groups = [(t, q) for t, q in groups if len(q)]
    fig, axes = plt.subplots(1, len(groups), figsize=(6.5 * len(groups), 4.8), squeeze=False)
    cols = {'halfA': '#d1495b', 'halfB': '#2e86ab', 'z8': '0.2'}
    for ax, (title, q) in zip(axes[0], groups):
        q = q.sort_values('variant')
        x = np.arange(len(q))
        bad = ((q.n_flagged > 0) | (q.n_at_bound > 0)).to_numpy()
        for t, c in cols.items():
            ax.scatter(x, q[f'{t}_pattern'], s=26, facecolors=['none' if b else c for b in bad], edgecolors=c,
                       label='all recordings' if t == 'z8' else f'half {t[-1]}')
            if len(ref):
                ax.axhline(float(ref[f'{t}_pattern'].iloc[0]), color=c, ls='--', lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(q.variant, rotation=70, ha='right', fontsize=7)
        ax.set_ylabel('pattern part of the APL-effect loss')
        ax.set_title(f"{title}, S's population, calibrated readout\n(dashed: the unchanged fits; hollow: outside the bands "
                     'or at a bound)', fontsize=8)
        ax.legend(fontsize=7)
    fig.tight_layout()
    path = out / 'E96_settings.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _md(df: pd.DataFrame, fmt: str = '.4f') -> str:
    return E94._md(df, fmt)


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    L = ["# E96: APL's synapses on the KCs outside the claws", '',
         f'Generated by experiments/E96_offclaw_apl.py at {H.code_commits().get("queue_commit", "")[:7]}.', '',
         'What stays non-anatomical in every setting: the part of the weights on the claws keeps the family\'s breadth '
         'weighting; the component is spread evenly over a KC\'s claws (where on the dendrite a synapse lies is not '
         'modelled); the counts are one hemibrain\'s. The peduncle\'s synapses are counted as acting on the KC\'s drive '
         '(B2 to B6); B1 does not rest on that.', '']
    if (out / 'identity.json').exists():
        s = _read_json(out / 'identity.json')
        L += ['## Part 0: the identities', '', f"passes: {s.get('passes')}", '']
        for k, x in s.items():
            if not isinstance(x, dict) or k == 'commits':
                continue
            keep = {kk: x.get(kk) for kk in ('passes', 'totals', 'presynaptic_site_only', 'loader_dropped', 'kept_is_loaders',
                                              'wnorm_proportional_to_kept', 'loader_drops_are_counts', 'same_as_real',
                                              'same_as_e95_a4', 'apl_factor', 'apl_factor_e94', 'apl_factor_e95',
                                              'apl_factor_e91', 'factor_rel_diff', 'same_weights', 'same_weights_as_e95')
                    if kk in x}
            if isinstance(x.get('counts'), dict):
                keep['counts_differ'] = x['counts'].get('counts_differ')
            L.append(f'- {k}: {json.dumps(keep, default=str)}')
        L.append('')
    if (out / 'd_check.json').exists():
        s = _read_json(out / 'd_check.json')
        L += ['## Dcheck: E10_joint_fit.py --aplkc-flat at size 0 against E94\'s D1 fit\'s first evaluation', '',
              '```', json.dumps(s, indent=1, default=str)[:4000], '```', '']
    for name, title, fmt in (('a1_counts.csv', 'A1: APL->KC synapses on the model\'s KCs (the real wiring\'s b_k)', '.3f'),
                             ('a1_by_type_class.csv', 'A1: by KC type and claw class', '.3f'),
                             ('a1_wirings.csv', 'A1: against the claw count, b_k and c_k on each wiring', '.3f'),
                             ('a1_per_claw.csv', 'A1: the per-claw spread', '.3f'),
                             ('a2_distances.csv', 'A2: distance to the KC\'s nearest claw (Pratyush\'s assignments, um)', '.3f'),
                             ('a2_connectomes.csv', 'A2: the four connectomes\' derived tables, calyx', '.3f'),
                             ('a2_connectomes_kc.csv', 'A2: per KC in the four connectomes', '.3f'),
                             ('a3_kc2apl.csv', 'A3: KC->APL synapses split the same way', '.3f'),
                             ('a4_shapes.csv', 'A4: the component\'s shape among E95\'s static forms', '.3f')):
        if (out / name).exists():
            d = pd.read_csv(out / name)
            if name == 'a4_shapes.csv':
                d = d[d.wiring == 'real']
            L += [f'## {title}', '', _md(d, fmt), '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        runs = runs_table(out)
        g = rd[(rd.unit == 'S') & (rd.readout == 'cal')].sort_values(['part', 'arm', 'variant'])
        graw = rd[(rd.unit == 'S') & (rd.readout == 'raw')].set_index('variant')
        rr = runs.groupby('variant')
        L += ["## S's population (calibrated readout; the parts of the pattern part against all recordings)", '',
              '| setting | pattern, half A | half B | all | six | amplitude | shape | spread ratio | Pearson r | ρ | loss | '
              'mean APL effect | raw-count loss | ester / alcohol error | factor | share at the boutons / component / '
              'weighted | flagged / at bound |', '|---|' + '---:|' * 16]

        def rng(q, c, spec='.3g'):
            if not len(q) or c not in q or not pd.to_numeric(q[c], errors='coerce').notna().any():
                return '–'
            v = pd.to_numeric(q[c], errors='coerce')
            return f"{v.min():{spec}}-{v.max():{spec}}"
        for r in g.itertuples():
            q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()
            sh = ' / '.join(rng(q, c, '.3f') for c in ('bouton_share_d', 'comp_share', 'weighted_share'))
            L.append(f'| {r.variant} | {r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | {r.z6_pattern:.4f} | '
                     f'{r.z8_amplitude:.4f} | {r.z8_shape:.4f} | {r.z8_sd_ratio:.2f} | {r.z8_pearson:.2f} | {r.z8_rho:.2f} | '
                     f'{r.z8_loss:.4f} | {r.mean_effect:+.3f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | '
                     f'{r.z8_err_ester:+.3f} / {r.z8_err_alcohol:+.3f} | {rng(q, "factor")} | {sh} | '
                     f'{r.n_flagged} / {r.n_at_bound} |')
        L.append('')
        mine = runs[runs.source == 'E96']
        if len(mine):
            L += ['## Realism of the runs (rule of 2026-11-10 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | '
                  'K | share at the boutons / component / weighted | synapses\' share | APL->KC weight p90/p10 over responders '
                  '| APL->bouton x E3 | flag | flag (earlier rule) |', '|---|' + '---:|' * 11 + '---|---|']
            for r in mine.itertuples():
                sh = ' / '.join(E94._f(getattr(r, c, np.nan), '.3f') for c in ('bouton_share_d', 'comp_share', 'weighted_share'))
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | '
                         f'{r.opk_on_tv:.3f} / {r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | '
                         f'{E94._f(getattr(r, "K_run", np.nan), ".3g")} | {sh} | {E94._f(getattr(r, "synapse_share", np.nan), ".3f")} | '
                         f'{E94._f(getattr(r, "wkc_all_resp_p90_over_p10", np.nan), ".3g")} | '
                         f'{E94._f(getattr(r, "aplpn_multiple", np.nan), ".3g")} | '
                         f'{r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
    for name in ('line.csv', 'choice.csv', 'gate.csv', 'equal.csv', 'vs_uniform.csv', 'd_scores.csv', 'steps.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                L += [f'## {name}', '', '(none)', '']
                continue
            if name == 'line.csv':
                d = d[[c for c in d.columns if c in ('variant', 'unit', 'readout', 'ester_rel', 'alcohol_rel', 'resid',
                                                       'off_toward_flies', 'off_toward_flies_nr')]]
            if name == 'steps.csv' and len(d) > 200 and 'read' in d:
                d = d[d.read.astype(bool)]
            L += [f'## {name}', '', _md(d), '']
    if (out / 'verdict.json').exists():
        vc = _read_json(out / 'verdict.json')
        short = {k: vc.get(k) for k in ('verdict', 'gate', 'equal_settings', 'D_plan', 'D', 'flagged_runs', 'at_bound_runs')}
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
    old.update(script='experiments/E96_offclaw_apl.py', part=part, runs=runs,
               settings=dict(ARMS={k: list(v) for k, v in ARMS.items()}, APL_BRACKET=APL_BRACKET, APL_STEP=APL_STEP,
                             D_DOSE=D_DOSE, B_SIZES=B_SIZES, POOL=POOL, PERM_SEED=PERM_SEED, MODES=MODES,
                             CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, CONTROL_RATIO=CONTROL_RATIO,
                             EQUAL=dict(PATTERN=EQUAL_PATTERN, RAW=EQUAL_RAW, MEAN=EQUAL_MEAN), D_COND=D_COND,
                             D_MODE_TIE=D_MODE_TIE, D_POOL_GAIN=D_POOL_GAIN, D_FLAT_BOUNDS=D_FLAT_BOUNDS, D_MAXFEV=D_MAXFEV,
                             D_TAGS=D_TAGS, D_AS_GOOD=D_AS_GOOD, MODE_NAMES=MODE_NAMES,
                             counts_source='pyolfsysm connectome_APL_KC_synapse_counts (APL>KC table of 2025-12-14, '
                                           'APL_R, both sites in the ROI)'))
    _write_json(path, old)


def main() -> None:
    import warnings
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'A1', 'A2', 'A3', 'A4', 'B', 'score', 'C', 'Dcheck', 'D',
                                                      'Dscore', 'E', 'steps', 'variants', 'table'))
    ap.add_argument('--arm', default=None, help='B: B1..B6 (comma-separated; default all)')
    ap.add_argument('--wirings', default=None, help="'S', 'T', 'all8' or a comma-separated list (0: all8; B: S; D: a lane)")
    ap.add_argument('--pop', default='S', choices=('S', 'T'), help='D: which fits')
    ap.add_argument('--start', default='own', choices=tuple(D_TAGS), help="D: the start ('own': the wiring's best point of "
                    "part B; 'e95': E95's refit with the component exchanged)")
    ap.add_argument('--mode', default=None, choices=MODES, help="D: the component (default the plan's; 'allsyn' for B6's "
                    'fits on their rule)')
    ap.add_argument('--force', action='store_true', help="D: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='D: print the plan and the commands only')
    ap.add_argument('--workers', type=int, default=4, help='score: co-response processes')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    ap.add_argument('--no-zcell', action='store_true', help="score: skip E73's per-KC readings")
    ap.add_argument('--variants', default=None, help='steps: comma-separated variants (default: chosen and equal settings)')
    a = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = out_dir()
    t0 = time.time()
    if a.part == '0':
        ws = _wirings(a.wirings or 'all8')
        provenance(out, '0', wirings=list(ws))
        part_0(out, ws)
    elif a.part == 'A':
        provenance(out, 'A')
        part_a(out)
    elif a.part in ('A1', 'A2', 'A3', 'A4'):
        {'A1': a1_counts, 'A2': a2_distances, 'A3': a3_kc2apl, 'A4': a4_shapes}[a.part](out)
    elif a.part == 'B':
        ws = _wirings(a.wirings or 'S')
        arms = [x for x in a.arm.split(',') if x] if a.arm else None
        provenance(out, 'B', wirings=list(ws), arms=arms)
        part_b(out, arms, ws)
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
            provenance(out, 'D', pop=a.pop, start=a.start, wirings=a.wirings, mode=a.mode)
        part_d(out, a.pop, a.start, None if a.wirings is None else _wirings(a.wirings), a.mode, force=a.force,
               dry_run=a.dry_run)
    elif a.part == 'Dscore':
        d_scores(out)
        if (out / 'readings.csv').exists():
            verdict96(out)
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
    print(f'E96 {a.part} done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
