#!/usr/bin/env python3
"""E95: what the best family's APL is for a KC: its inhibition at the PN boutons, and then all of it, replaced by a static
inhibition per KC that is present only with APL on (each KC's own loss of drive, or its number of claws alone); the same
offset made at the KCs in the feedback APL with the bouton path at a tenth; and the family fitted again that way.

The entry is item E95 of experiment_priorities.md (proposed by viki 2026-11-14 06:01; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 05:18).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E94's `WiringModel94`: E83's `WiringModel` through E91's `WiringModel91`, with a dose on the APL->bouton scale). S: real,
P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's, C / D the second split's (E89's `targets89`).

Definitions, per wiring at its fit, from E90's archived stage runs (window means over the 0.5 s pulse) and the simulator's
claw x bouton weights (E94's `wiring_W`, `kc_drive`): L^b(k, o) the loss of KC k's summed claw drive made at the boutons
for odor o (the drive from the boutons without APL minus that with APL), L^c(k, o) the loss made at the claws (the drive
from the boutons with APL minus the claws' states), b_k and c_k their means over the 16 odors, n_k the KC's number of claws
(`kc_losses`). A static inhibition is olfsysm's `ext_inh` (subtracted from a KC's input while the odor is on; KCs by
odors), set for the APL-on run only (`set_ext_inh_ko`, which maps the outputs' odor columns to olfsysm's odor order and
refuses a negative or non-finite value in Python: olfsysm's own check is inside its parallel region and aborts the
process). The APL-off run of every point is the saved fit's. The KC-KC scale and the somatic threshold stay at the fit's.

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; parts B, C and D refuse to run until it passes): on every
  wiring, a static inhibition of zero (a KCs x odors matrix of zeros, set) at the fit's dose and weights reproduces the
  saved APL-on run, and APL removed (every APL scale 0) with a static inhibition of zero the saved APL-off run, with 0
  count differences; on the real wiring the KC-by-odor setter with the same value for every odor is E91's per-KC setter
  (the same olfsysm matrix and the same counts), a matrix non-zero for one odor changes only that odor's column, and part
  D's weights with no unweighted component at dose 0.1 reproduce E94's saved point B_sub_k2_d0.1 (its factor and its
  counts). Each wiring's per-KC thresholds (fixed_thr + 2 x spont_in), claws, connectome and family APL->KC weights are
  cached under outputs/E95_static_apl/W/ for part A.
- A (no simulations; part 0's cache and E94's wiring cache): A1 the split of L^b and L^c into a value per KC, per odor and
  the rest, by claw class and over the KCs that respond without APL (a1_split.csv); A2 b_k against n_k (the slope a claw,
  the scatter, the log-log slope, what carries the rest: the APL->bouton weights of the KC's boutons, its glomeruli
  against E75's map, its type) and c_k against the breadth factor and n_k (a2_per_kc.csv, a2_glomeruli.csv); A3 b_k, c_k
  and their sum against each KC's threshold, and the responses without APL a threshold raise of b_k, or b_k + c_k, would
  remove by the window-mean drive, by odor group and the co-responses by group of pairs (a3_sizes.csv, a3_removed.csv; a
  count-level description only); A4 APL's activity by odor and before the odor in the fits and in E94's saved runs, the
  connectome's per-claw APL->KC weights on every wiring (a4_apl.csv, a4_weights.csv); A5 where E87's and E88's static
  forms, E91's static control and E80's threshold difference sit (a5_archive.csv). a_summary.json.
- B, C, D: the point lists (`--arm`, `--wirings`, resumable; runs_B.csv, runs_C.csv, runs_D.csv; runs under
  outputs/E95_static_apl/runs/ with their stage matrices and each static inhibition saved per KC):
  - B, the bouton path as a static inhibition (dose 0; the claw path's feedback kept with the family's weights; the tied
    APL scale bisected in APL_BRACKET): B1 b_k; B2 beta n_k (beta the wiring's least-squares slope through the origin);
    B3 n_k^e for e in B3_EXPONENTS scaled to b's mean; B4 b_k permuted among the KCs with the same number of claws and
    over all KCs (seed PERM_SEED); B5 L^b(k, o).
  - C, the whole APL as a static inhibition (every APL scale 0; the static inhibition times a factor bisected in
    STATIC_BRACKET): C1 b_k + c_k; C2 L^b + L^c; C3 beta' n_k (beta' the slope of b + c); C4 beta n_k + c_k; C5 C1's
    values permuted over all KCs; C6 the claw path alone static (dose 1, the APL->KC weights 0, c_k, the tied scale
    bisected in APL_BRACKET).
  - D, the offset at the KCs in the feedback APL (dose D_DOSE; the APL->KC weight of a claw the family's plus an
    unweighted component, d in units of what the bouton path took, `component_unit`; the tied scale bisected): D1 the
    same weight on every claw, D2 the connectome's per-claw weights unweighted by breadth, both at D_SIZES; D3 D1 at
    d = 1 with no bouton path. D0 (dose 0.1, no component: E94's B_sub_k2_d0.1) is the control of part D's settings: on S
    E94's saved runs, on T run here (part E).
- score (`--no-co`, `--no-zcell`): readings.csv per variant and unit (each wiring, S and T where all four have a run) and
  readout ('cal' the calibrated readout, K set on each run's APL-on counts, the reading of the choice and the gate;
  'calK' the fit's K; 'raw'; 'zcell' and 'zcell_lin' E73's per-KC reading under the saturating and the linear map) against
  every target of `targets89` (E92's `score_row92`: loss, pattern part with its amplitude and shape parts, spread ratio,
  Pearson r and ρ on the per-pair APL effect, the two groups' errors, the mean APL effect); pairs.csv, line.csv (E92's
  part A line), coresponse.csv (E89's readings of the alcohol pairs' co-response change); then part E's tables.
- E: the choice per half (`choose95`), T's runs of a chosen setting and its control (runs_E.csv), the gate (E92's
  `gate_direction92` / `gate92`, E91's rules), the equal rule in three kinds (`equal_rule95`), verdict.json.
- F: the family fitted again at the connectome's bouton share with the unweighted component (`E10_joint_fit.py -k 2
  --aplpn-hold F --aplkc-flat uniform ...`, `f_command`), only on its condition (`f_condition`); `--part Fscore` scores
  them by E94's D1 rule (`f_scores`).
- G: documented stubs (`part_g`): run only by the entry's rules.
- steps (`--variants`): the step from the claws to the somata in wild-type flies for chosen or equal settings (E94's
  practice).
- variants: variants.csv and model_pair_distances_long.csv for the archive.
- table: comparison_table.md and the figure.

- Fcheck: E10_joint_fit.py with and without --aplkc-flat (at size 0) against E94's D1 fit on the real wiring's saved first
  evaluation (`f_check`: its command with --smoke, 5 evaluations, history into this experiment's results folder);
  f_check.json. Run once before part F.

Deviations from the entry, forced by the code:
- Part F's fifth parameter is outside olfsysm's simulator, which optimize_mb_model_params could not fit: pyolfsysm's
  optimize_mb_model_params gained `extra_params` (olfsysm 1bb4b90, Python only, no C++ change; None, the default, changes
  nothing), each set through the simulator's closure before every evaluation; E10_joint_fit.py's --aplkc-flat uses it.
  The size's unit (d = 1: a claw loses the fit's mean loss a claw at the boutons) needs the fit's stage run, which E10 does
  not read, so E10 takes it as --aplkc-flat-unit (`component_unit`, from part D's runs) with --aplkc-flat-init (default
  --init-params' 'aplkc_flat') and --aplkc-flat-max (8).
- "Every APL scale 0" (part C) sets all four (APL->KC, KC->APL, APL->bouton, bouton->APL); part 0 checks that this with a
  static inhibition of zero is the saved APL-off run (whose two APL output scales are 0) exactly.
- "The static inhibition against the somatic threshold" is against each KC's threshold as the simulator sets it,
  fixed_thr + 2 x its spontaneous input (median about 2.5 x fixed_thr); against fixed_thr alone beside (the premise's "about
  half the threshold" is against fixed_thr).
- D0 (dose 0.1 without the component, part D's control) on S is E94's saved runs (part 0 reproduces the real wiring's
  exactly), on T run here.
- The gate's control of a static setting ('<v>_ctrl', its values permuted over all KCs with seed PERM_SEED) is run as its
  own point even where it is the same permutation as B4's or C5's.
- Part F starts each wiring from its own best point of D1 (the lowest pattern part of that wiring under the readout
  against all recordings; points inside the bands first, none at a bracket's bound).
- A3's "responses a threshold raise would remove" are read on the window-mean summed drive against each KC's threshold
  (a count-level description, as the entry says); A2's "against E75's map" is the glomeruli's mean loss a claw at the
  boutons against E75's ridge profile (E83's `e75_profile('ridge_x1')`); the per-KC reading "under both maps" is E93's
  (saturating, K calibrated) and the linear map (noise and noise ratio calibrated on E73's targets).
- olfsysm's odor order is the outputs' column order on these simulators, so the one-odor identity checks the middle column
  (it would check the column whose olfsysm index differs most otherwise).

Outputs: results/E95_static_apl/; spike counts, responses, stage matrices and static inhibitions under
outputs/E95_static_apl/ (runs/, W/, F/).

Usage (from the repository root; one 16-thread process at a time on viki; check `uptime` first):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E95_static_apl --part 0
    ./venv/bin/python -u -m experiments.E95_static_apl --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E95_static_apl --part B --arm B1
    ./venv/bin/python -u -m experiments.E95_static_apl --part score
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

EXPERIMENT = 'E95_static_apl'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

APL_BRACKET = E92G.APL_BRACKET                      # (0.01, 10): the tied APL scale (B, C6, D), E92's bracket
STATIC_BRACKET = (0.05, 20.0)                       # C: the static inhibition's factor (the probes: 1.24 to 1.48)
APL_STEP = E83.APL_STEP
RAW_CRIT = E92G.RAW_CRIT                            # a model response is two or more spikes
CLAW_GROUPS = E94.CLAW_GROUPS                       # 1-4, 5-8, 9+ claws
B3_EXPONENTS = (0.0, 0.5, 1.5, 2.0)                 # B3: n_k^e scaled to b's mean (1 is B2)
PERM_SEED = 1                                       # B4, C5 and the gate's controls: one seed
D_DOSE = 0.1                                        # D: E94's anatomical dose
D_SIZES = (0.5, 1.0, 2.0)                           # D1, D2: the component's size in units of the bouton path
D3_SIZE = 1.0
COMPONENTS = ('uniform', 'connectome')
CONNECTOME_SHARE_RANGE = (0.096, 0.118)             # the connectome's share of the claws' loss made at the boutons
E94_D0_TAG = 'B_sub_k2_d0.1'                        # E94's point at dose 0.1 (D's control; part 0's identity)

# the variants of each arm (88 points on S: B 36, C 24, D 28)
ARMS: Dict[str, Tuple[str, ...]] = {
    'B1': ('B1_bk',),
    'B2': ('B2_nk',),
    'B3': tuple(f'B3_nk{e:g}' for e in B3_EXPONENTS),
    'B4': (f'B4_bk_permclaw{PERM_SEED}', f'B4_bk_permall{PERM_SEED}'),
    'B5': ('B5_bko',),
    'C1': ('C1_bck',),
    'C2': ('C2_bcko',),
    'C3': ('C3_nk',),
    'C4': ('C4_nkck',),
    'C5': (f'C5_bck_permall{PERM_SEED}',),
    'C6': ('C6_ck',),
    'D1': tuple(f'D1_uniform_d{d:g}' for d in D_SIZES),
    'D2': tuple(f'D2_connectome_d{d:g}' for d in D_SIZES),
    'D3': (f'D3_uniform_d{D3_SIZE:g}_dose0',),
}
PART_OF_ARM = {a: a[0] for a in ARMS}
REF = 'ref'                                         # the unchanged fits (source 'saved')
D0 = 'D0'                                           # dose 0.1 without the component (E94's point): the control of D

# part E (the entry's rules: E91's, as in E94)
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
CONTROL_RATIO = E91.CONTROL_RATIO
HALVES = E91.HALVES
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
EQUAL_KINDS = {'bouton_static': ('B1', 'B2', 'B5'), 'whole_static': ('C1', 'C2', 'C4', 'C6'), 'at_kcs': ('D1', 'D2', 'D3')}
CHOICE_ARMS = ('B1', 'B2', 'B3', 'B5', 'C1', 'C2', 'C3', 'C4', 'C6', 'D1', 'D2', 'D3')   # B4 and C5 excepted

# part F
F_COND = 0.008                                      # S's best of D1 or D3 within this of the unchanged fits' 0.0500
F_FLAT_BOUNDS = (0.0, 8.0)                          # the component's size, a free parameter
F_MAXFEV = 200
F_TAG = '_f1'
F_START_FILE = 'f_start_params.json'
E94_REFIT_S = 0.0593                                # E94's refit at the same bouton share without the component (S)
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
FLIES_MEAN = E89.FLIES_MEAN

_CACHE: Dict[tuple, object] = {}
APL_SCALES = ('wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')   # C: every APL scale 0


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


def _flagged(x) -> bool:
    return E94._flagged(x)


def _as_bool(x) -> bool:
    return E94._as_bool(x)


def _commits(strict: bool = True) -> dict:
    """H.code_commits() at the start of a simulation part, which must be clean (E89's `_commits_clean`)."""
    return E89._commits_clean() if strict else dict(H.code_commits())


def _wirings(arg: str) -> Tuple[str, ...]:
    return E91._wirings(arg)


# ---------------------------------------------------------------------------
# names

_FIXED = {'B1_bk': ('B1', 'b'), 'B2_nk': ('B2', 'nk_b'), 'B5_bko': ('B5', 'b_o'), 'C1_bck': ('C1', 'bc'),
          'C2_bcko': ('C2', 'bc_o'), 'C3_nk': ('C3', 'nk_bc'), 'C4_nkck': ('C4', 'nk_b+c'), 'C6_ck': ('C6', 'c')}
_RE_B3 = re.compile(r'^B3_nk([0-9.]+)$')
_RE_PERM = re.compile(r'^(B4_bk|C5_bck)_perm(claw|all)(\d+)$')
_RE_D = re.compile(r'^D([12])_(uniform|connectome)_d([0-9.]+)$')
_RE_D3 = re.compile(r'^D3_uniform_d([0-9.]+)_dose0$')
STATIC_NAMES = {'b': 'b_k (each KC\'s mean loss at the boutons)', 'nk_b': 'beta n_k (b\'s slope through the origin)',
                'nke_b': 'n_k^e scaled to b\'s mean', 'b_o': 'L^b(k, o)', 'bc': 'b_k + c_k', 'bc_o': 'L^b + L^c',
                'nk_bc': 'beta\' n_k (b + c\'s slope)', 'nk_b+c': 'beta n_k + c_k', 'c': 'c_k (each KC\'s mean loss at the claws)'}
MODE_NAMES = {'bouton': 'the bouton path removed (dose 0), a static inhibition, the claw path\'s feedback kept, the tied APL '
                        'scale re-held',
              'whole': 'every APL scale 0, a static inhibition times a factor re-held',
              'claw': 'the bouton path\'s feedback as fitted, the APL->KC weights 0, a static inhibition c_k, the tied APL '
                      'scale re-held',
              'kcs': 'the APL->KC weights the family\'s plus an unweighted component, the bouton path at a dose, the tied APL '
                     'scale re-held', 'ref': 'the unchanged fit'}


def _kind_of_arm(arm: str) -> Optional[str]:
    return next((k for k, arms in EQUAL_KINDS.items() if arm in arms), None)


def parse_variant(v: str) -> dict:
    """A variant's fields: arm, part ('B' / 'C' / 'D' / 'ref'), mode ('bouton' the bouton path static, 'whole' every APL
    scale 0, 'claw' C6, 'kcs' part D, 'ref'), static (the static inhibition's form, STATIC_NAMES, or None), exponent,
    permute ('claw' / 'all' / None) and seed, component, size and dose (part D), control and control_of, screen (the arm,
    'ref' or 'control') and the kind of equal it can be (EQUAL_KINDS; None for B3, C3, the permutations and the
    controls). '<v>_ctrl' is a static setting's control (its values permuted over all KCs, seed PERM_SEED); D0 the control
    of part D's settings (dose 0.1 without the component)."""
    out = dict(variant=v, arm=None, part=None, mode=None, static=None, exponent=np.nan, permute=None, seed=None,
               component=None, size=np.nan, dose=np.nan, control=False, control_of='', screen=None, equal_kind=None)
    if v == REF:
        out.update(arm=REF, part=REF, mode=REF, dose=1.0, screen=REF)
        return out
    if v == D0:
        out.update(arm=D0, part='D', mode='kcs', dose=D_DOSE, size=0.0, control=True, control_of='D', screen='control')
        return out
    ctl = v.endswith('_ctrl')
    core = v[:-len('_ctrl')] if ctl else v
    if core in _FIXED:
        arm, st = _FIXED[core]
        out.update(arm=arm, static=st)
    elif _RE_B3.match(core):
        e = _RE_B3.match(core).group(1)
        if f'{float(e):g}' != e or float(e) not in B3_EXPONENTS:
            raise ValueError(f'variant {v!r}: B3\'s exponents are {B3_EXPONENTS}')
        out.update(arm='B3', static='nke_b', exponent=float(e))
    elif _RE_PERM.match(core):
        m = _RE_PERM.match(core)
        arm = 'B4' if m.group(1) == 'B4_bk' else 'C5'
        if arm == 'C5' and m.group(2) != 'all':
            raise ValueError(f'variant {v!r}: C5 permutes over all KCs')
        out.update(arm=arm, static='b' if arm == 'B4' else 'bc', permute=m.group(2), seed=int(m.group(3)))
    elif _RE_D.match(core) or _RE_D3.match(core):
        m = _RE_D.match(core)
        if m:
            arm, comp, d = f'D{m.group(1)}', m.group(2), m.group(3)
            if (arm == 'D1') != (comp == 'uniform'):
                raise ValueError(f'variant {v!r}: D1 is the uniform component, D2 the connectome\'s')
            dose = D_DOSE
        else:
            arm, comp, d, dose = 'D3', 'uniform', _RE_D3.match(core).group(1), 0.0
        if f'{float(d):g}' != d:
            raise ValueError(f'variant {v!r}: the size is not written canonically')
        out.update(arm=arm, component=comp, size=float(d), dose=dose)
    else:
        raise ValueError(f'variant {v!r}: not one of this experiment\'s')
    arm = out['arm']
    out['part'] = arm[0]
    out['mode'] = 'kcs' if out['part'] == 'D' else ('claw' if arm == 'C6' else ('bouton' if out['part'] == 'B' else 'whole'))
    if out['part'] != 'D':
        out['dose'] = 1.0 if out['mode'] == 'claw' else 0.0
    if ctl:
        if out['part'] == 'D' or out['permute'] is not None:
            raise ValueError(f'variant {v!r}: part D\'s control is D0; the permutations have none')
        out.update(control=True, control_of=core, permute='all', seed=PERM_SEED, screen='control')
    else:
        out['screen'] = arm
        out['equal_kind'] = None if out['permute'] else _kind_of_arm(arm)
    return out


def control_of(variant: str) -> str:
    """The gate's control of a setting: a static setting's values permuted over the KCs ('<v>_ctrl', seed PERM_SEED);
    part D's settings: D0 (the same setting without the component)."""
    pv = parse_variant(variant)
    if pv['control'] or pv['permute'] or pv['part'] == REF:
        raise ValueError(f'{variant}: no control (a control, a permutation or the reference)')
    return D0 if pv['part'] == 'D' else f'{variant}_ctrl'


def arm_variants(arm: str) -> List[str]:
    return list(ARMS[arm])


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


# ---------------------------------------------------------------------------
# the losses per KC, the static matrices, the component

def kc_losses(wiring: str, kc_index: Optional[pd.Index] = None) -> dict:
    """A wiring's fit, per KC (the outputs' row order `kc_index`, default the saved fit's) and odor (the outputs' column
    order, abbreviated), from E90's stage runs (window means over the pulse) and E94's wiring cache: 'Lb' the loss of the
    summed claw drive made at the boutons (the drive from the boutons without APL minus that with APL), 'Lc' the loss made
    at the claws (the drive from the boutons with APL minus the claws' states), 'Doff' the summed claw drive without APL,
    'b' and 'c' their means over the odors, 'n' the claws per KC; 'bclaw_mean' the mean loss a claw at the boutons (over
    claws and odors), 'apl_mean' APL's mean window activity over the odors, 'odors', 'kc_index', 'W'."""
    key = ('losses', wiring)
    if key in _CACHE:
        return _CACHE[key]
    W, _ = E94.wiring_W(wiring)
    on, off = E94.fit_stage(wiring, 'on'), E94.fit_stage(wiring, 'off')
    for st, c in ((on, 'on'), (off, 'off')):
        E94.check_stage_shapes(W, st, f'{wiring} {c}')
    saved = E83.saved_fit(wiring)['saved']['on'][1]
    kc_index = saved.index if kc_index is None else kc_index
    odors = E94._odors(on['bouton_win'])
    if odors != [H.odor_abbrev(c) for c in saved.columns]:
        raise RuntimeError(f"{wiring}: the stage runs' odors {odors} are not the outputs' columns")
    Wm = W.to_numpy(dtype=float)
    x_off = Wm @ off['bouton_win'][odors].to_numpy(dtype=float)
    x_on = Wm @ on['bouton_win'][odors].to_numpy(dtype=float)
    c_on = on['claw_win'][odors].to_numpy(dtype=float)
    c_off = off['claw_win'][odors].to_numpy(dtype=float)
    kc_of = E94.kc_of_claw(W)
    ids = np.asarray(kc_index.get_level_values(H.KC_ID_LEVEL) if isinstance(kc_index, pd.MultiIndex) else kc_index)

    def per_kc(M: np.ndarray) -> pd.DataFrame:
        d = pd.DataFrame(M, columns=odors).groupby(kc_of).sum().reindex(ids)
        if d.isna().any().any():
            raise RuntimeError(f'{wiring}: {int(d.isna().any(axis=1).sum())} KCs of the outputs have no claws in the wiring')
        d.index = kc_index
        return d
    Lb, Lc, Doff = per_kc(x_off - x_on), per_kc(x_on - c_on), per_kc(c_off)
    n = pd.Series(pd.Series(kc_of).value_counts().reindex(ids).to_numpy(dtype=float), index=kc_index, name='n_claws')
    apl = on['apl'].reindex(odors)['win'].astype(float)
    res = dict(Lb=Lb, Lc=Lc, Doff=Doff, b=Lb.mean(axis=1), c=Lc.mean(axis=1), n=n, odors=odors, kc_index=kc_index, W=W,
               bclaw_mean=float((x_off - x_on).mean()), apl_mean=float(apl.mean()), apl_win=apl, kc_of_claw=kc_of)
    _CACHE[key] = res
    return res


def slope_through_origin(x: np.ndarray, y: np.ndarray) -> float:
    """The least-squares slope of y on x through the origin."""
    x, y = np.asarray(x, dtype=float).ravel(), np.asarray(y, dtype=float).ravel()
    xx = float(x @ x)
    if not xx > 0:
        raise ValueError('slope through the origin of an all-zero x')
    return float(x @ y) / xx


def permute_values(v: np.ndarray, seed: int, groups: Optional[np.ndarray] = None) -> np.ndarray:
    """`v`'s rows (a vector or a KCs x odors matrix) permuted over all rows, or within each group of `groups` (e.g. the KCs'
    claw counts; groups in sorted order, one generator `np.random.default_rng(seed)`)."""
    v = np.asarray(v)
    rng = np.random.default_rng(int(seed))
    if groups is None:
        return v[rng.permutation(len(v))]
    g = np.asarray(groups)
    if len(g) != len(v):
        raise ValueError(f'{len(g)} groups for {len(v)} rows')
    out = v.copy()
    for key in np.unique(g):
        idx = np.flatnonzero(g == key)
        out[idx] = v[idx[rng.permutation(idx.size)]]
    return out


def static_matrix(pv: dict, losses: dict, record: Optional[dict] = None) -> pd.DataFrame:
    """The static inhibition (KCs x odors, the outputs' order, clipped at 0) of a parsed variant (`parse_variant`) from
    `kc_losses`: 'b' b_k; 'nk_b' beta n_k (beta b's slope through the origin on n_k); 'nke_b' n_k^e scaled to b's mean over
    the KCs; 'b_o' L^b; 'bc' b_k + c_k; 'bc_o' L^b + L^c; 'nk_bc' beta' n_k (the slope of b + c); 'nk_b+c' beta n_k + c_k;
    'c' c_k; the per-KC forms the same for every odor; then the variant's permutation of the KCs (within a claw count or
    over all). `record` gets the clipping ('n_clipped', 'min_before_clip') and the slope where one is used."""
    st = pv['static']
    n = losses['n'].to_numpy(dtype=float)
    b, c = losses['b'].to_numpy(dtype=float), losses['c'].to_numpy(dtype=float)
    n_od = len(losses['odors'])
    rec = {} if record is None else record
    v = None
    if st == 'b':
        v = b
    elif st == 'nk_b':
        rec['slope'] = slope_through_origin(n, b)
        v = rec['slope'] * n
    elif st == 'nke_b':
        x = n ** float(pv['exponent'])
        v = x * (b.mean() / x.mean())
    elif st == 'bc':
        v = b + c
    elif st == 'nk_bc':
        rec['slope'] = slope_through_origin(n, b + c)
        v = rec['slope'] * n
    elif st == 'nk_b+c':
        rec['slope'] = slope_through_origin(n, b)
        v = rec['slope'] * n + c
    elif st == 'c':
        v = c
    if v is not None:
        A = np.repeat(np.asarray(v, dtype=float)[:, None], n_od, axis=1)
    elif st == 'b_o':
        A = losses['Lb'].to_numpy(dtype=float).copy()
    elif st == 'bc_o':
        A = (losses['Lb'] + losses['Lc']).to_numpy(dtype=float)
    else:
        raise ValueError(f'no static inhibition of form {st!r}')
    if pv.get('permute'):
        A = permute_values(A, pv['seed'], n if pv['permute'] == 'claw' else None)
    rec['n_clipped'] = int((A < 0).sum())
    rec['min_before_clip'] = float(A.min())
    return pd.DataFrame(np.clip(A, 0.0, None), index=losses['kc_index'], columns=list(losses['odors']))


def component_unit(bclaw_mean: float, apl_scale: float, apl_mean: float, w_mean: float = 1.0) -> float:
    """The unscaled APL->KC weight, per unit of a component's mean weight `w_mean`, that at the fit's APL->KC scale and
    APL's mean activity takes from a claw the fit's mean loss a claw at the boutons (`bclaw_mean`), wherever it is not
    clipped: the component's size d = 1 (inhibition = scale x weight x APL's activity)."""
    den = float(apl_scale) * float(apl_mean) * float(w_mean)
    if not (np.isfinite(den) and den > 0 and np.isfinite(bclaw_mean)):
        raise ValueError(f'component unit: scale {apl_scale}, APL {apl_mean}, mean weight {w_mean}')
    return float(bclaw_mean) / den


def d_weights(w_fam: pd.Series, wnorm: pd.Series, kind: str, d: float, unit: float) -> pd.Series:
    """Part D's APL->KC weights per claw (unscaled): the family's (`w_fam`) plus d x unit x the component ('uniform': the
    connectome's mean weight on every claw; 'connectome': each claw's connectome weight, unweighted by breadth), so both
    components have the same mean; d = 0 gives `w_fam` itself."""
    if not wnorm.index.equals(w_fam.index):
        raise ValueError('the connectome weights are not indexed as the family\'s')
    if not (np.isfinite(d) and d >= 0):
        raise ValueError(f'component size {d}')
    if d == 0:
        return w_fam.copy()
    if kind == 'uniform':
        comp = np.full(len(wnorm), float(wnorm.mean()))
    elif kind == 'connectome':
        comp = wnorm.to_numpy(dtype=float)
    else:
        raise ValueError(f'unknown component {kind!r}')
    return pd.Series(w_fam.to_numpy(dtype=float) + float(d) * float(unit) * comp, index=w_fam.index, name=w_fam.name)


# ---------------------------------------------------------------------------
# the simulator: the static inhibition KCs x odors

def set_ext_inh_ko(mp, rv, M: Optional[np.ndarray], order: Optional[np.ndarray] = None, tol: float = 1e-9) -> dict:
    """Set olfsysm's static inhibition (`rv.kc.ext_inh`, subtracted from each KC's input while the odor is on, times
    `mp.kc.ext_inh_scale`) from a KCs x odors matrix `M` in the outputs' column order (`order[j]`: olfsysm's odor index
    of output column j, the simulator's `odor_order`; default the identity). Checked in Python, since olfsysm's own check
    runs inside its parallel region, where a refusal aborts the process: the KC count, finite values, values below
    -tol x max refused and the rest clipped at 0, ext_inh_scale 1, n_claws_active_to_spike not positive. None removes it.
    Returns {'n_clipped'}."""
    if M is None:
        rv.kc.ext_inh = np.zeros((0, 0))
        return {}
    A = np.asarray(M, dtype=float)
    n_kc = int(mp.kc.N)
    if A.ndim != 2 or A.shape[0] != n_kc:
        raise ValueError(f'ext_inh: {A.shape} for {n_kc} KCs (KCs x odors)')
    n_od = A.shape[1]
    order = np.arange(n_od) if order is None else np.asarray(order, dtype=int)
    if order.shape != (n_od,) or not np.array_equal(np.sort(order), np.arange(n_od)):
        raise ValueError(f'ext_inh: the odor order {order} is not a permutation of the {n_od} columns')
    if not np.isfinite(A).all():
        raise ValueError('ext_inh: not finite')
    scale = max(1.0, float(np.abs(A).max()))
    if (A < -tol * scale).any():
        raise ValueError(f'ext_inh: negative values (min {A.min():.4g}); olfsysm refuses them')
    if float(mp.kc.ext_inh_scale) != 1.0:
        raise RuntimeError(f'ext_inh_scale is {mp.kc.ext_inh_scale}, not 1')
    if int(getattr(mp.kc, 'n_claws_active_to_spike', -1)) > 0:
        raise RuntimeError('ext_inh with n_claws_active_to_spike > 0: olfsysm refuses it')
    n_clip = int((A < 0).sum())
    out = np.zeros((n_kc, n_od))
    out[:, order] = np.clip(A, 0.0, None)
    rv.kc.ext_inh = out
    return dict(n_clipped=n_clip)


class StaticModel(E94.WiringModel94):
    """E94's `WiringModel94` (one wiring at its fit with the lab's input; a dose on the APL->bouton scale; olfsysm's `mp`,
    `rv`, `order`) with this experiment's runs, each with the static inhibition `ext` (KCs x odors, the outputs' order;
    None for none) set for that run only and removed after: `run_on(x)` APL on with the tied APL scales times x at the
    dose; `run_static(f)` every APL scale 0 with the static inhibition times f; `run_off()` the fit's APL-off parameters,
    never with a static inhibition."""

    def __init__(self, label: str, check: bool = True, **model_kws):
        super().__init__(label, check=check, **model_kws)
        self.ext: Optional[np.ndarray] = None

    def params_static(self) -> Dict[str, float]:
        return dict(self.p, **{k: 0.0 for k in APL_SCALES})

    def _run(self, params: Dict[str, float], M: Optional[np.ndarray]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        set_ext_inh_ko(self.mp, self.rv, M, None if M is None else self.order)
        try:
            return self.sim(params)
        finally:
            set_ext_inh_ko(self.mp, self.rv, None)

    def run_on(self, x: float = 1.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self._run(self.params_on(x), self.ext)

    def run_static(self, f: float = 1.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if self.ext is None:
            raise RuntimeError('run_static without a static inhibition')
        return self._run(self.params_static(), float(f) * np.asarray(self.ext, dtype=float))

    def run_off(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self._run(self.p_off(), None)


def _as_kc_odor(M, model) -> np.ndarray:
    """A static matrix as a KCs x odors array in the model's outputs' order (a DataFrame is checked against them)."""
    if isinstance(M, pd.DataFrame):
        on = model.saved['on'][1]
        if not (M.index.equals(on.index) and list(M.columns) == [H.odor_abbrev(c) for c in on.columns]):
            raise ValueError('the static inhibition is not in the outputs\' KC and odor order')
    A = np.asarray(M, dtype=float)
    if A.shape != tuple(model.saved['on'][1].shape):
        raise ValueError(f'static inhibition {A.shape}, outputs {model.saved["on"][1].shape}')
    return A


def reading95(model: StaticModel, pv: dict, M: Optional[pd.DataFrame], wk: Optional[pd.Series], label: str = '') -> dict:
    """One setting (`parse_variant`) on a wiring at its fit: the APL->KC weights `wk` (None: the fit's), the variant's
    dose, the static inhibition `M` (KCs x odors, APL-on runs only); the bisection to the fit's APL-on count (E69's
    `bisect_factor` from 1: the tied APL scale in APL_BRACKET, or with every APL scale 0 ('whole') the static factor in
    STATIC_BRACKET), then the run again with APL's activity saved (its counts must be equal), read for E90's stage
    matrices and E92's `gain_stats`. The weights, the dose, the static inhibition and `save_inh_sims` are reset on
    return, also on an error."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E90_model_stages as E90
    t0 = time.time()
    whole = pv['mode'] == 'whole'
    if whole and M is None:
        raise ValueError('every APL scale 0 needs a static inhibition')
    ext = None if M is None else _as_kc_odor(M, model)
    try:
        E89.set_weights(model, wk)
        model.dose = 1.0 if whole else float(pv['dose'])
        model.ext = ext
        E94.set_form94(model.mp, 'sub')
        run, bracket = (model.run_static, STATIC_BRACKET) if whole else (model.run_on, APL_BRACKET)
        x, on, hist, at_bound = E69.bisect_factor(run, model.target_on, model.size, bracket, APL_STEP,
                                                  f"{label} {'static factor' if whole else 'APL scale'}")
        params = model.params_static() if whole else model.params_on(x)
        model.mp.kc.save_inh_sims = True
        try:
            on2 = run(x)
        finally:
            model.mp.kc.save_inh_sims = False
        ns = int((on2[1].to_numpy(dtype=float) != on[1].to_numpy(dtype=float)).sum())
        if ns:
            raise RuntimeError(f"{model.label}: the run again with APL's activity saved differs in {ns} counts")
        odors = [H.odor_abbrev(c) for c in on[1].columns]
        st = E90.read_stages(model.mp, model.rv, model.W, model.order, odors, E90.MODEL_WINDOW)
        gs = E92G.gain_stats(model.rv, model.mp, model.order, on[1].to_numpy(dtype=float) > RAW_CRIT, 'sub')
        ext = None if model.ext is None else (float(x) if whole else 1.0) * model.ext
    finally:
        model.ext = None
        set_ext_inh_ko(model.mp, model.rv, None)
        model.dose = 1.0
        E89.set_weights(model, None)
    return dict(on=on, off=model.saved['off'], factor=float(x), at_bound=bool(at_bound), n_sims=len(hist) + 1,
                params=params, seconds=time.time() - t0, history=hist, stages=st, gain=gs, bracket=tuple(bracket), ext=ext,
                whole=whole)


# ---------------------------------------------------------------------------
# realism beside the four bands

def _q(v: np.ndarray, prefix: str) -> dict:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if not v.size:
        return {f'{prefix}_{n}': np.nan for n in ('p10', 'p50', 'p90', 'max')}
    q = np.quantile(v, [0.1, 0.5, 0.9, 1.0])
    return {f'{prefix}_{n}': float(x) for n, x in zip(('p10', 'p50', 'p90', 'max'), q)}


def static_realism(M: Optional[np.ndarray], thr: np.ndarray, n_claws: np.ndarray, d_off: np.ndarray,
                   on_counts: Optional[np.ndarray] = None, off_counts: Optional[np.ndarray] = None,
                   fixed_thr: Optional[float] = None) -> dict:
    """The static inhibition (KCs x odors, as run; None: none) per KC (its mean over the odors) as a multiple of the KC's
    somatic threshold `thr` (fixed_thr + 2 x its spontaneous input): 10th percentile, median, 90th, largest, and the median
    by claw class; 'static_over_max_drive' the share of KCs whose static inhibition exceeds their largest summed claw drive
    without APL (`d_off`, window means; any odor); 'silent_by_mean_drive' the share whose largest drive without APL minus
    the static inhibition stays below the threshold for every odor; with the runs' counts, the shares of KCs that respond
    to no odor with APL on and off (two or more spikes); with `fixed_thr`, the median as a multiple of it beside."""
    row = {}
    thr = np.asarray(thr, dtype=float).ravel()
    n = np.asarray(n_claws, dtype=float).ravel()
    D = np.asarray(d_off, dtype=float)
    if M is None:
        row.update(_q([], 'static_over_thr'))
        row.update({f'static_over_thr_p50_{g}': np.nan for g, _, _ in CLAW_GROUPS})
        row.update(static_over_max_drive=np.nan, silent_by_mean_drive=np.nan, static_mean=np.nan, static_max=np.nan,
                   static_over_fixed_thr_p50=np.nan)
    else:
        A = np.asarray(M, dtype=float)
        s = A.mean(axis=1)
        r = s / thr
        row.update(_q(r, 'static_over_thr'))
        for g, lo, hi in CLAW_GROUPS:
            m = (n >= lo) & (n <= hi)
            row[f'static_over_thr_p50_{g}'] = float(np.median(r[m])) if m.any() else np.nan
        row['static_over_max_drive'] = float((s > D.max(axis=1)).mean())
        row['silent_by_mean_drive'] = float(((D - A).max(axis=1) < thr).mean())
        row['static_mean'] = float(A.mean())
        row['static_max'] = float(A.max())
        row['static_over_fixed_thr_p50'] = float(np.median(s) / fixed_thr) if fixed_thr else np.nan
    for lab, C in (('on', on_counts), ('off', off_counts)):
        if C is not None:
            row[f'silent_share_{lab}'] = float((~(np.asarray(C, dtype=float) > RAW_CRIT).any(axis=1)).mean())
    return row


def d_realism(W, b_off: np.ndarray, b_on: np.ndarray, c_on: np.ndarray, w_fam: np.ndarray, w_comp: np.ndarray,
              kc_of: np.ndarray, responders: Optional[np.ndarray] = None) -> dict:
    """Part D's realism from a run's stage matrices (window means; claws x odors, boutons x odors) and the scaled APL->KC
    weights per claw split into the family's `w_fam` and the unweighted component `w_comp`: the share of the claws' loss
    made at the boutons ('bouton_share', E94's `claw_loss_shares`), by the component ('comp_share') and by the weighted
    part ('weighted_share'; each claw's loss at the claws split in proportion to its two weights; the three add to 1);
    the spread of the per-KC summed APL->KC weight (the family's alone and with the component) over the KCs `responders`
    (KC ids that respond without APL; all KCs if None): 10th percentile, median, 90th and the 90th over the 10th."""
    Wm = np.asarray(W, dtype=float)
    x_off, x_on = Wm @ np.asarray(b_off, dtype=float), Wm @ np.asarray(b_on, dtype=float)
    c_on = np.asarray(c_on, dtype=float)
    wf, wc = np.asarray(w_fam, dtype=float).ravel(), np.asarray(w_comp, dtype=float).ravel()
    if not (wf.size == wc.size == x_on.shape[0]):
        raise ValueError(f'{wf.size} / {wc.size} weights for {x_on.shape[0]} claws')
    tot = float((x_off - c_on).sum())
    at_claw = x_on - c_on
    wsum = wf + wc
    frac = np.divide(wc, wsum, out=np.zeros_like(wsum), where=wsum > 0)
    comp = float((at_claw * frac[:, None]).sum())
    row = dict(bouton_share_d=float((x_off - x_on).sum()) / tot if tot > 0 else np.nan,
               comp_share=comp / tot if tot > 0 else np.nan,
               weighted_share=float(at_claw.sum() - comp) / tot if tot > 0 else np.nan)
    kc_of = np.asarray(kc_of)
    for lab, w in (('fam', wf), ('all', wsum)):
        per = pd.Series(w).groupby(kc_of).sum()
        if responders is not None:
            per = per.reindex(np.asarray(responders)).dropna()
        q = np.quantile(per.to_numpy(dtype=float), [0.1, 0.5, 0.9]) if len(per) else [np.nan] * 3
        row.update({f'wkc_{lab}_resp_p10': float(q[0]), f'wkc_{lab}_resp_p50': float(q[1]), f'wkc_{lab}_resp_p90': float(q[2]),
                    f'wkc_{lab}_resp_p90_over_p10': float(q[2] / q[0]) if q[0] > 0 else np.nan})
    return row


# ---------------------------------------------------------------------------
# the per-KC and per-claw tables (part 0's cache, part A's input)

def _w_dir() -> Path:
    return local_dir('W')


def kc_tables(model: StaticModel) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per KC, per claw) of a wiring at its fit: per KC (the outputs' rows) its type, claws, spontaneous input,
    threshold (fixed_thr + 2 x spont_in, as the simulator sets it), the family's breadth factor (exponent 2) and its
    summed unscaled APL->KC weights (the family's and the connectome's); per claw its KC, the connectome's unscaled weight
    and the family's."""
    P = E83.wiring_parts(model.label)
    wn = P['wnorm']['wAPLKC']
    if not wn.index.equals(model.wAPLKC.index):
        raise RuntimeError(f'{model.label}: the connectome APL->KC weights are not indexed as the family\'s')
    kc_index = model.saved['on'][1].index
    ids = np.asarray(kc_index.get_level_values(H.KC_ID_LEVEL))
    spont = np.asarray(model.rv.kc.spont_in, dtype=float).ravel()
    if spont.size != len(ids):
        raise RuntimeError(f'{spont.size} spontaneous inputs for {len(ids)} KCs')
    kc_of = np.asarray(wn.index.get_level_values(H.KC_ID_LEVEL))
    W = model.W
    if not np.array_equal(np.asarray(W.index.get_level_values(H.KC_ID_LEVEL)), kc_of):
        raise RuntimeError(f'{model.label}: the APL->KC weights are not in the wiring\'s claw order')
    claws = pd.DataFrame({H.KC_ID_LEVEL: kc_of, 'wnorm': wn.to_numpy(dtype=float),
                          'w_fam': model.wAPLKC.to_numpy(dtype=float)})
    g = claws.groupby(H.KC_ID_LEVEL)
    fac = E89.weighting_factors(P['breadth'], 2.0).reindex(ids)
    kc = pd.DataFrame({H.KC_ID_LEVEL: ids, H.KC_TYPE_LEVEL: np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL)),
                       'n_claws': g.size().reindex(ids).to_numpy(dtype=float), 'spont_in': spont,
                       'thr': float(model.p['fixed_thr']) + 2.0 * spont, 'factor': fac.to_numpy(dtype=float),
                       'w_fam_sum': g.w_fam.sum().reindex(ids).to_numpy(dtype=float),
                       'wnorm_sum': g.wnorm.sum().reindex(ids).to_numpy(dtype=float)})
    return kc, claws


def save_kc_tables(label: str, kc: pd.DataFrame, claws: pd.DataFrame) -> None:
    kc.to_parquet(_w_dir() / f'{label}_kc.parquet')
    claws.to_parquet(_w_dir() / f'{label}_claw.parquet')


def load_kc_tables(label: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ps = [_w_dir() / f'{label}_{k}.parquet' for k in ('kc', 'claw')]
    if not all(p.exists() for p in ps):
        raise RuntimeError(f'{label}: no per-KC table (run --part 0 first: it caches them under {_w_dir()})')
    return tuple(pd.read_parquet(p) for p in ps)


# ---------------------------------------------------------------------------
# saving and loading runs

STAGE_KEYS = E94.STAGE_KEYS


def _dirs(sub: str = 'runs') -> List[Path]:
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    return [H.outputs_dir(EXPERIMENT) / sub, arc / sub, arc / 'from_hal' / sub]


def save_run95(tag: str, o: Dict[str, tuple], sub: str = 'runs', stages: Optional[dict] = None,
               drive: Optional[pd.DataFrame] = None, loss: Optional[pd.DataFrame] = None,
               static: Optional[pd.DataFrame] = None, weights: Optional[pd.Series] = None) -> None:
    """A run's outputs (APL on, and the saved APL-off run beside), its stage matrices with APL on, each KC's summed claw
    drive from the boutons, E92's loss of drive at the claws, the static inhibition per KC and odor as run, and the
    APL->KC weights per claw where they are not the fit's, under outputs/E95_static_apl/<sub>/."""
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
    for name, x in (('drive', drive), ('loss', loss), ('static', static)):
        if x is not None:
            x = x.copy()
            x.columns = [str(c) for c in x.columns]
            x.to_parquet(d / f'{tag}_{name}.parquet')
    if weights is not None:
        weights.rename('wAPLKC').reset_index(drop=True).to_frame().to_parquet(d / f'{tag}_wAPLKC.parquet')


def load_run95(tag: str, sub: str = 'runs') -> Optional[Dict[str, tuple]]:
    for d in _dirs(sub):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def load_stage95(tag: str, sub: str = 'runs') -> Optional[dict]:
    for d in _dirs(sub):
        if (d / f'{tag}_on_claw_win.parquet').exists():
            st = {k: pd.read_parquet(d / f'{tag}_on_{k}.parquet') for k in STAGE_KEYS}
            st['apl'] = pd.read_csv(d / f'{tag}_on_apl.csv', index_col=0)
            return st
    return None


def load_frame95(tag: str, name: str, sub: str = 'runs') -> pd.DataFrame:
    for d in _dirs(sub):
        p = d / f'{tag}_{name}.parquet'
        if p.exists():
            return pd.read_parquet(p)
    raise RuntimeError(f'no {name} saved for {tag}')


def _append(path: Path, row: dict) -> None:
    E83._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def _frame(a: np.ndarray, like: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(np.asarray(a, dtype=float), index=like.index, columns=[str(c) for c in like.columns])


# ---------------------------------------------------------------------------
# part 0: the identities

def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('passes'))


def _counts_diff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    return E94._counts_diff(a, b)


def _column_changes(a: pd.DataFrame, b: pd.DataFrame) -> List[int]:
    """The output columns in which two spike-count frames differ."""
    return [int(j) for j in np.flatnonzero((a.to_numpy(dtype=float) != b.to_numpy(dtype=float)).any(axis=0))]


def d_unit(model: StaticModel, losses: Optional[dict] = None) -> float:
    """`component_unit` of a wiring at its fit: its mean loss a claw at the boutons, the fit's APL->KC scale, APL's mean
    window activity (E90's stage run) and the connectome's mean APL->KC weight."""
    L = kc_losses(model.label, model.saved['on'][1].index) if losses is None else losses
    wn = E83.wiring_parts(model.label)['wnorm']['wAPLKC']
    return component_unit(L['bclaw_mean'], model.p['wAPLKC_scale'], L['apl_mean'], float(wn.mean()))


def part_0(out: Path, wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per wiring); raises if one fails."""
    import gc
    cc = _commits()
    path = out / 'identity.json'
    res = _read_json(path)
    res.setdefault('commits', cc)
    real_keys = ('setter_real', 'one_odor_real', 'd0_real')
    for w in wirings:
        key = f'zero_{w}'
        if res.get(key, {}).get('exact') and (w != 'real' or all(res.get(k, {}).get('passes') for k in real_keys)):
            continue
        t0 = time.time()
        m = StaticModel(w, check=False)
        on_fit = m.run_on(1.0)                                      # no static inhibition (sets the thresholds)
        save_kc_tables(w, *kc_tables(m))
        Z = np.zeros(m.saved['on'][1].shape)
        m.ext = Z
        seen = {}
        sim0 = m.sim

        def spy(params, _sim=sim0, _seen=seen):
            _seen.setdefault('ext_shapes', []).append(list(np.asarray(m.rv.kc.ext_inh).shape))
            return _sim(params)
        m.sim = spy
        try:
            on = m.run_on(1.0)
            off = m.run_static(1.0)
        finally:
            m.sim = sim0
            m.ext = None
        d_on, d_off = _counts_diff(on, m.saved['on']), _counts_diff(off, m.saved['off'])
        reset = int(np.asarray(m.rv.kc.ext_inh).size) == 0
        res[key] = dict(wiring=w, on=d_on, off=d_off, fit_on=_counts_diff(on_fit, m.saved['on']),
                        ext_set_in_runs=seen.get('ext_shapes'), ext_reset=reset,
                        exact=bool(d_on['exact'] and d_off['exact'] and reset
                                   and seen.get('ext_shapes') == [list(Z.shape)] * 2), seconds=time.time() - t0)
        print(f"E95 0 {key}: {res[key]}", flush=True)
        _write_json(path, res)
        if w == 'real':
            L = kc_losses(w, m.saved['on'][1].index)
            n_od = len(m.order)
            v = L['b'].to_numpy(dtype=float)
            # the KC-by-odor setter with one value per KC against E91's per-KC setter
            E91.set_ext_inh(m.mp, m.rv, v, n_od)
            A91 = np.array(m.rv.kc.ext_inh, dtype=float, copy=True)
            set_ext_inh_ko(m.mp, m.rv, np.repeat(v[:, None], n_od, axis=1), m.order)
            Ako = np.array(m.rv.kc.ext_inh, dtype=float, copy=True)
            set_ext_inh_ko(m.mp, m.rv, None)
            m.dose = 0.0
            try:
                E91.set_ext_inh(m.mp, m.rv, v, n_od)
                try:
                    a = m.sim(m.params_on(1.0))
                finally:
                    E91.set_ext_inh(m.mp, m.rv, None)
                m.ext = np.repeat(v[:, None], n_od, axis=1)
                b = m.run_on(1.0)
            finally:
                m.ext, m.dose = None, 1.0
            ds = _counts_diff(a, b)
            res['setter_real'] = dict(same_matrix=bool(np.array_equal(A91, Ako)), counts=ds,
                                      changed_from_fit=_counts_diff(b, m.saved['on'])['counts_differ'],
                                      passes=bool(np.array_equal(A91, Ako) and ds['exact']))
            print(f"E95 0 setter_real: {res['setter_real']}", flush=True)
            # one odor: the column of the odor whose olfsysm index differs most from its output column (the middle
            # column when the two orders are the same)
            order = np.asarray(m.order, dtype=int)
            j = int(np.argmax(np.abs(order - np.arange(n_od)))) if (order != np.arange(n_od)).any() else n_od // 2
            Mj = np.zeros((len(v), n_od))
            Mj[:, j] = (L['Lb'] + L['Lc']).to_numpy(dtype=float)[:, j]
            m.ext = Mj
            try:
                r = m.run_on(1.0)
            finally:
                m.ext = None
            ch = _column_changes(r[1], m.saved['on'][1])
            res['one_odor_real'] = dict(column=j, odor=L['odors'][j], olfsysm_index=int(order[j]),
                                        order_is_identity=bool((order == np.arange(n_od)).all()), changed_columns=ch,
                                        n_changed=int((r[1].to_numpy() != m.saved['on'][1].to_numpy()).sum()),
                                        ext_reset=int(np.asarray(m.rv.kc.ext_inh).size) == 0,
                                        passes=bool(ch == [j] and int(np.asarray(m.rv.kc.ext_inh).size) == 0))
            print(f"E95 0 one_odor_real: {res['one_odor_real']}", flush=True)
            _write_json(path, res)
            # part D's weights with no component at dose 0.1 against E94's saved point
            t2 = time.time()
            wn = E83.wiring_parts(w)['wnorm']['wAPLKC']
            unit = d_unit(m, L)
            wk = d_weights(m.wAPLKC, wn, 'uniform', 0.0, unit)
            r94 = pd.read_csv(H.results_dir(E94.EXPERIMENT) / 'runs_B1.csv').drop_duplicates('tag', keep='last')
            r94 = r94.set_index('tag').loc[f'{E94_D0_TAG}_real']
            o94 = E94.load_run(f'{E94_D0_TAG}_real')
            if o94 is None:
                raise RuntimeError(f"E94's point {E94_D0_TAG} on the real wiring is not saved")
            rd = reading95(m, parse_variant(D0), None, wk, label='identity D0 real')
            dd = _counts_diff(rd['on'], o94['on'])
            f94 = float(r94['apl_factor'])
            res['d0_real'] = dict(same_weights=bool(np.array_equal(wk.to_numpy(), m.wAPLKC.to_numpy())), counts=dd,
                                  apl_factor=rd['factor'], apl_factor_e94=f94, factor_rel_diff=abs(rd['factor'] - f94) / f94,
                                  n_sims=rd['n_sims'], unit=unit, weights_reset=bool(np.array_equal(m._kc, m._kc0)),
                                  dose_after=float(m.dose), seconds=time.time() - t2)
            res['d0_real']['passes'] = bool(res['d0_real']['same_weights'] and dd['exact']
                                            and res['d0_real']['factor_rel_diff'] <= 1e-12
                                            and res['d0_real']['weights_reset'] and m.dose == 1.0)
            print(f"E95 0 d0_real: {res['d0_real']}", flush=True)
            _write_json(path, res)
        del m
        gc.collect()
    keys = [f'zero_{w}' for w in ALL_WIRINGS]
    have_all = all(k in res for k in keys) and all(k in res for k in real_keys)
    failed = [k for k in keys if k in res and not res[k]['exact']] + [k for k in real_keys if k in res and not res[k]['passes']]
    res['passes'] = bool(have_all and not failed)
    res['complete'] = bool(have_all)
    _write_json(path, res)
    if failed:
        raise RuntimeError(f'E95 part 0: the identity fails for {failed} (identity.json); the experiment stops here')
    print(f"E95 0: {'passes' if res['passes'] else 'so far exact (not all wirings yet)'}", flush=True)
    return res


# ---------------------------------------------------------------------------
# part A (no simulations)

def two_way_split(A: np.ndarray) -> dict:
    """The variance of a KCs x odors matrix split into a value per KC, a value per odor and the rest (shares adding to 1),
    and the variation over the odors at the median KC (the coefficient of variation per KC, over KCs with a positive
    mean)."""
    A = np.asarray(A, dtype=float)
    gm = A.mean()
    kc = A.mean(axis=1, keepdims=True) - gm
    od = A.mean(axis=0, keepdims=True) - gm
    rest = A - gm - kc - od
    tot = float(((A - gm) ** 2).sum())
    mu = A.mean(axis=1)
    ok = mu > 0
    cv = A[ok].std(axis=1) / mu[ok]
    if not tot > 0:
        return dict(share_kc=np.nan, share_odor=np.nan, share_rest=np.nan, cv_odors_median_kc=np.nan)
    return dict(share_kc=float((kc ** 2).sum() * A.shape[1] / tot), share_odor=float((od ** 2).sum() * A.shape[0] / tot),
                share_rest=float((rest ** 2).sum() / tot), cv_odors_median_kc=float(np.median(cv)) if cv.size else np.nan)


def _claw_group(n: np.ndarray) -> np.ndarray:
    out = np.full(len(n), '', dtype=object)
    for g, lo, hi in CLAW_GROUPS:
        out[(n >= lo) & (n <= hi)] = g
    return out


def _responders_mask(wiring: str) -> np.ndarray:
    off = E83.saved_fit(wiring)['saved']['off'][1].to_numpy(dtype=float)
    return (off > RAW_CRIT).any(axis=1)


def _corr(x, y) -> Tuple[float, float]:
    from scipy import stats
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan, np.nan
    return float(np.corrcoef(x[ok], y[ok])[0, 1]), float(stats.spearmanr(x[ok], y[ok])[0])


def a1_split(out: Path) -> pd.DataFrame:
    """A1 (module docstring): a1_split.csv, per wiring, site ('boutons' L^b, 'claws' L^c) and KC set ('all', 'responders'
    the KCs that respond without APL): the two-way split, the mean per KC and by claw class."""
    rows = []
    for w in ALL_WIRINGS:
        L = kc_losses(w)
        resp = _responders_mask(w)
        n = L['n'].to_numpy()
        thr = float(E83.saved_fit(w)['p']['fixed_thr'])
        for site, X in (('boutons', L['Lb']), ('claws', L['Lc'])):
            A = X.to_numpy(dtype=float)
            for kset, m in (('all', np.ones(len(A), dtype=bool)), ('responders', resp)):
                mu = A[m].mean(axis=1)
                row = dict(wiring=w, site=site, kcs=kset, n_kcs=int(m.sum()), fixed_thr=thr, mean_per_kc=float(mu.mean()),
                           **two_way_split(A[m]))
                for g, lo, hi in CLAW_GROUPS:
                    mg = (n[m] >= lo) & (n[m] <= hi)
                    row[f'mean_{g}'] = float(mu[mg].mean()) if mg.any() else np.nan
                    row[f'n_kcs_{g}'] = int(mg.sum())
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a1_split.csv', index=False)
    return df


def _kc_bouton_apl(W: pd.DataFrame, wapl: np.ndarray, kc_index: pd.Index) -> pd.Series:
    """Per KC the mean unscaled APL->bouton weight of the boutons its claws sample, weighted by the claws' wiring
    weights."""
    Wm = W.to_numpy(dtype=float)
    num = pd.Series(Wm @ np.asarray(wapl, dtype=float)).groupby(E94.kc_of_claw(W)).sum()
    den = pd.Series(Wm.sum(axis=1)).groupby(E94.kc_of_claw(W)).sum()
    ids = np.asarray(kc_index.get_level_values(H.KC_ID_LEVEL))
    return pd.Series((num / den).reindex(ids).to_numpy(dtype=float), index=kc_index)


def a2_per_kc(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A2 (module docstring): a2_per_kc.csv (per wiring: b_k against n_k, what carries the rest, c_k against the breadth
    factor, n_k and the summed weights) and a2_glomeruli.csv (per wiring and glomerulus: the mean loss a claw at the
    boutons over its claws, relative to the wiring's mean, the glomerulus' mean APL->bouton weight, E75's map)."""
    from scipy import stats
    rows, grows = [], []
    try:
        e75 = E83.e75_profile('ridge_x1')
    except Exception as e:                       # (a description; reported without the map when E75's file is absent)
        print(f'E95 A2: E75\'s map not read ({e!r})', flush=True)
        e75 = {}
    for w in ALL_WIRINGS:
        L = kc_losses(w)
        kc, claws = load_kc_tables(w)
        W, wapl = E94.wiring_W(w)
        ids = np.asarray(L['kc_index'].get_level_values(H.KC_ID_LEVEL))
        if not np.array_equal(kc[H.KC_ID_LEVEL].to_numpy(), ids):
            raise RuntimeError(f'{w}: the per-KC table is not in the outputs\' order')
        n, b, c = L['n'].to_numpy(), L['b'].to_numpy(), L['c'].to_numpy()
        beta = slope_through_origin(n, b)
        res = b - beta * n
        ok = b > 0
        pc = b / n
        p = E83.saved_fit(w)['p']
        wb = _kc_bouton_apl(W, wapl, L['kc_index']).to_numpy()
        pot = p['wAPLPN_scale'] * L['apl_mean'] * pd.Series(W.to_numpy(dtype=float) @ wapl).groupby(
            E94.kc_of_claw(W)).sum().reindex(ids).to_numpy(dtype=float)
        drive_claw = L['Doff'].to_numpy().mean(axis=1) / n
        types = kc[H.KC_TYPE_LEVEL].astype(str).to_numpy()
        tm = pd.Series(res).groupby(types).transform('mean').to_numpy()
        row = dict(wiring=w, n_kcs=len(n), slope_per_claw=beta,
                   r2_origin=float(1 - ((b - beta * n) ** 2).sum() / ((b - b.mean()) ** 2).sum()),
                   pearson_b_n=_corr(n, b)[0], spearman_b_n=_corr(n, b)[1],
                   loglog_slope=float(np.polyfit(np.log(n[ok]), np.log(b[ok]), 1)[0]),
                   per_claw_mean=float(pc.mean()), per_claw_cv=float(pc.std() / pc.mean()),
                   resid_sd=float(res.std()), resid_over_b_sd=float(res.std() / b.std()),
                   pearson_resid_bouton_apl=_corr(res, wb)[0], pearson_perclaw_bouton_apl=_corr(pc, wb)[0],
                   pearson_b_potential=_corr(b, pot)[0], pearson_perclaw_drive=_corr(pc, drive_claw)[0],
                   pearson_resid_drive=_corr(res, drive_claw)[0],
                   resid_type_eta2=float(((tm - res.mean()) ** 2).sum() / ((res - res.mean()) ** 2).sum()),
                   c_mean=float(c.mean()), pearson_c_factor=_corr(kc.factor, c)[0], spearman_c_factor=_corr(kc.factor, c)[1],
                   pearson_c_n=_corr(n, c)[0], spearman_c_n=_corr(n, c)[1],
                   pearson_c_wfam=_corr(kc.w_fam_sum, c)[0], pearson_c_wnorm=_corr(kc.wnorm_sum, c)[0],
                   slope_c_per_claw=slope_through_origin(n, c), slope_bc_per_claw=slope_through_origin(n, b + c))
        rows.append(row)
        # per glomerulus: the claws' loss at the boutons (each claw's bouton of largest weight)
        Wm = W.to_numpy(dtype=float)
        on, off = E94.fit_stage(w, 'on'), E94.fit_stage(w, 'off')
        od = L['odors']
        lc = (Wm @ (off['bouton_win'][od].to_numpy(dtype=float) - on['bouton_win'][od].to_numpy(dtype=float))).mean(axis=1)
        gl = np.asarray(W.columns.get_level_values('glomerulus')).astype(str)[Wm.argmax(axis=1)]
        bw = pd.Series(np.asarray(wapl, dtype=float), index=np.asarray(W.columns.get_level_values('glomerulus')).astype(str))
        gd = pd.DataFrame(dict(glomerulus=gl, loss=lc)).groupby('glomerulus').loss.agg(['mean', 'size'])
        for g, r in gd.iterrows():
            grows.append(dict(wiring=w, glomerulus=g, n_claws=int(r['size']), claw_loss_rel=float(r['mean'] / lc.mean()),
                              apl_bouton_weight_rel=float(bw[g].mean() / bw.mean()), e75_score=e75.get(g, np.nan)))
    df, gdf = pd.DataFrame(rows), pd.DataFrame(grows)
    if len(gdf) and gdf.e75_score.notna().any():
        for w, g in gdf.groupby('wiring'):
            q = g[g.e75_score.notna()]
            i = df.index[df.wiring == w][0]
            df.loc[i, 'spearman_glom_loss_e75'] = float(stats.spearmanr(q.claw_loss_rel, q.e75_score)[0]) if len(q) > 2 else np.nan
            df.loc[i, 'spearman_glom_loss_aplw'] = float(stats.spearmanr(g.claw_loss_rel, g.apl_bouton_weight_rel)[0])
    df.to_csv(out / 'a2_per_kc.csv', index=False)
    gdf.to_csv(out / 'a2_glomeruli.csv', index=False)
    return df, gdf


def a3_sizes(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3 (module docstring): a3_sizes.csv (per wiring and quantity b, c, b + c: per KC as a multiple of its threshold,
    fixed_thr + 2 x its spontaneous input, 10th percentile, median, 90th, largest, and the median by claw class; beside, as
    a multiple of fixed_thr alone, the fits' thresholds the entry's premise quotes) and a3_removed.csv (per wiring, raise (b,
    b + c) and group: the responses without APL, those whose window-mean drive is at or above the threshold, and those of
    them a threshold raise would bring below it; per odor group (E90's ester and alcohol odors, the others) and, for the
    co-responses (KCs responding to both odors of a pair), per group of pairs (E92's part A groups). A count-level
    description: the window-mean drive does not decide a response, transients do)."""
    from experiments import E90_model_stages as E90
    from experiments import E92_steps_like_for_like as A
    srows, rrows = [], []
    for w in ALL_WIRINGS:
        L = kc_losses(w)
        kc, _ = load_kc_tables(w)
        thr = kc.thr.to_numpy(dtype=float)
        n = L['n'].to_numpy()
        grp = _claw_group(n)
        b, c = L['b'].to_numpy(), L['c'].to_numpy()
        fthr = float(E83.saved_fit(w)['p']['fixed_thr'])
        for q, v in (('b', b), ('c', c), ('b+c', b + c)):
            r = v / thr
            row = dict(wiring=w, quantity=q, mean_abs=float(v.mean()), thr_median=float(np.median(thr)), fixed_thr=fthr,
                       **_q(r, 'over_thr'), **_q(v / fthr, 'over_fixed_thr'))
            for g, _, _ in CLAW_GROUPS:
                row[f'over_thr_p50_{g}'] = float(np.median(r[grp == g])) if (grp == g).any() else np.nan
            srows.append(row)
        off = E83.saved_fit(w)['saved']['off'][1]
        R = off.to_numpy(dtype=float) > RAW_CRIT
        D = L['Doff'].to_numpy(dtype=float)
        odors = list(L['odors'])
        margin = D - thr[:, None]
        og = np.array(['ester' if o in E90.ESTER else ('alcohol' if o in E90.ALCOHOL else 'other') for o in odors])
        for raise_name, rv in (('b', b), ('b+c', b + c)):
            above = R & (margin >= 0)
            removed = above & (margin < rv[:, None])
            for g in ('ester', 'alcohol', 'other', 'all'):
                cols = np.ones(len(odors), dtype=bool) if g == 'all' else (og == g)
                rrows.append(dict(wiring=w, raise_=raise_name, kind='responses', group=g, n=int(R[:, cols].sum()),
                                  n_above_thr=int(above[:, cols].sum()), n_removed=int(removed[:, cols].sum()),
                                  share_removed_of_above=float(removed[:, cols].sum() / max(above[:, cols].sum(), 1))))
            pairs = pd.MultiIndex.from_tuples([tuple(sorted((odors[i], odors[j]))) for i in range(len(odors))
                                               for j in range(i + 1, len(odors))])
            ij = [(odors.index(a), odors.index(bb)) for a, bb in pairs]
            co = np.array([(above[:, i] & above[:, j]).sum() for i, j in ij])
            co_kept = np.array([((above[:, i] & ~removed[:, i]) & (above[:, j] & ~removed[:, j])).sum() for i, j in ij])
            gm = A.group_masks(pairs)
            for g, mk in list(gm.items()) + [('all120', np.ones(len(ij), dtype=bool))]:
                if not mk.any():
                    continue
                rrows.append(dict(wiring=w, raise_=raise_name, kind='co_responses', group=g, n_above_thr=int(co[mk].sum()),
                                  n_removed=int((co[mk] - co_kept[mk]).sum()),
                                  share_removed_of_above=float((co[mk] - co_kept[mk]).sum() / max(co[mk].sum(), 1))))
    s, r = pd.DataFrame(srows), pd.DataFrame(rrows).rename(columns={'raise_': 'raise'})
    s.to_csv(out / 'a3_sizes.csv', index=False)
    r.to_csv(out / 'a3_removed.csv', index=False)
    return s, r


E94_STAGE_DOSES = (0.0, 0.1, 0.25, 0.5, 2.0, 4.0)


def _apl_stats(apl: pd.DataFrame) -> dict:
    win, pre = apl['win'].astype(float), apl['pre'].astype(float)
    return dict(apl_win_mean=float(win.mean()), apl_win_cv=float(win.std() / win.mean()) if win.mean() > 0 else np.nan,
                apl_win_rel_max=float((win / win.mean() - 1).abs().max()) if win.mean() > 0 else np.nan,
                apl_pre_over_win=float(pre.mean() / win.mean()) if win.mean() > 0 else np.nan)


def a4_apl_weights(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A4 (module docstring): a4_apl.csv (APL's activity by odor and before the odor: the fits' stage runs on every wiring
    and E94's saved runs of B1 on S at its doses) and a4_weights.csv (per wiring the connectome's per-claw APL->KC weights:
    claws without a weight, the spread over the claws, the per-KC sum against the claw count and against b_k, the family's
    factor's quantiles, the inhibition a claw of the connectome's mean weight takes at the fit against the mean loss a claw
    at the boutons, the component's unit)."""
    arows, wrows = [], []
    for w in ALL_WIRINGS:
        arows.append(dict(wiring=w, source='fit', dose=1.0, **_apl_stats(E94.fit_stage(w, 'on')['apl'])))
        for d in E94_STAGE_DOSES:
            st = E94.load_stage94(f"{E94.variant_name('sub', 2.0, d)}_{w}")
            if st is not None:
                arows.append(dict(wiring=w, source='E94_B1', dose=d, **_apl_stats(st['apl'])))
        L = kc_losses(w)
        kc, claws = load_kc_tables(w)
        v = claws.wnorm.to_numpy(dtype=float)
        p = E83.saved_fit(w)['p']
        f = kc.factor.to_numpy(dtype=float)
        unw = p['wAPLKC_scale'] * v.mean() * L['apl_mean']
        wrows.append(dict(wiring=w, n_claws=int(v.size), zero_share=float((v == 0).mean()), cv_claws=float(v.std() / v.mean()),
                          kc_zero_share=float((kc.wnorm_sum == 0).mean()),
                          pearson_kcsum_n=_corr(kc.wnorm_sum, kc.n_claws)[0], pearson_kcsum_b=_corr(kc.wnorm_sum, L['b'])[0],
                          factor_p50=float(np.median(f)), factor_p90=float(np.quantile(f, 0.9)),
                          factor_p99=float(np.quantile(f, 0.99)), unweighted_inh_per_claw=float(unw),
                          bclaw_mean=L['bclaw_mean'], bclaw_over_unweighted=float(L['bclaw_mean'] / unw),
                          component_unit=component_unit(L['bclaw_mean'], p['wAPLKC_scale'], L['apl_mean'], float(v.mean())),
                          wAPLKC_scale=p['wAPLKC_scale'], apl_mean=L['apl_mean']))
    a, b = pd.DataFrame(arows), pd.DataFrame(wrows)
    a.to_csv(out / 'a4_apl.csv', index=False)
    b.to_csv(out / 'a4_weights.csv', index=False)
    return a, b


A5_COLS = ('z8_pattern', 'halfA_pattern', 'halfB_pattern', 'z6_pattern', 'z8_loss', 'mean_effect', 'z8_rho')


def a5_archive(out: Path) -> pd.DataFrame:
    """A5 (module docstring): a5_archive.csv, the static forms of the archive beside the unchanged fits: E87's part D
    (readings_D.csv: the whole APL replaced by a threshold raise by the odor's strength and the KC's breadth; S, T and the
    eight), E88's part D1 (readings_D1.csv, S), E91's static control (readings.csv, the controls), E80's threshold lower
    with APL on (screen.csv, B1 and B1_fit0, the real wiring), each under the calibrated readout and on raw counts where
    the table has them."""
    rows = []
    rdir = H.shared_dir() / 'results'

    def take(src, df, label_col, keep_units, extra=None):
        for r in df.to_dict('records'):
            if r.get('unit') not in keep_units:
                continue
            row = dict(source=src, form=str(r.get(label_col)), unit=r.get('unit'), readout=r.get('readout'))
            row.update({c: r.get(c) for c in A5_COLS})
            row.update({f'unchanged_{c}': r.get(f'unchanged_{c}') for c in A5_COLS if f'unchanged_{c}' in r})
            row.update(extra or {})
            rows.append(row)
    try:
        take('E87_D', pd.read_csv(rdir / 'E87_inhibition_form' / 'readings_D.csv'), 'form', ('S', 'T', 'all8', 'real'))
    except OSError as e:
        rows.append(dict(source='E87_D', note=repr(e)))
    try:
        take('E88_D1', pd.read_csv(rdir / 'E88_responder_sets' / 'readings_D1.csv'), 'form', ('S', 'real'))
    except OSError as e:
        rows.append(dict(source='E88_D1', note=repr(e)))
    try:
        d = pd.read_csv(rdir / 'E91_summed_drive_inhibition' / 'readings.csv')
        d = d[d.control.astype(str).isin(['True', '1']) | d.variant.isin(['B1_F_k2', 'B2_F_k2', 'B2_F_k1_s1'])]
        take('E91', d, 'variant', ('S', 'T', 'real'))
    except OSError as e:
        rows.append(dict(source='E91', note=repr(e)))
    try:
        s = pd.read_csv(rdir / 'E80_remaining_candidates' / 'screen.csv')
        s = s[s.handle.isin(['B1', 'B1_fit0']) | s.point_id.isin(['base', 'base_fit0'])]
        for r in s.to_dict('records'):
            for ro in ('cal', 'raw'):
                rows.append(dict(source='E80_B1', form=r['point_id'], unit='real', readout=ro, z8_pattern=r.get(f'{ro}_pattern'),
                                 z8_loss=r.get(f'{ro}_loss'), mean_effect=r.get(f'{ro}_mean'), z8_rho=r.get(f'{ro}_rho'),
                                 value=r.get('value')))
    except OSError as e:
        rows.append(dict(source='E80_B1', note=repr(e)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a5_archive.csv', index=False)
    return df


def part_a(out: Path) -> dict:
    """A1 to A5 and a_summary.json (the ranges over the eight wirings of what the entry's premise quotes)."""
    t0 = time.time()
    a1 = a1_split(out)
    a2, _ = a2_per_kc(out)
    a3, a3r = a3_sizes(out)
    a4, a4w = a4_apl_weights(out)
    a5 = a5_archive(out)

    def rng(df, col):
        v = pd.to_numeric(df[col], errors='coerce')
        return dict(min=float(v.min()), max=float(v.max()), mean=float(v.mean()))
    summ = dict(seconds=None)
    for site in ('boutons', 'claws'):
        q = a1[(a1.site == site) & (a1.kcs == 'all')]
        summ[f'a1_{site}'] = {c: rng(q, c) for c in ('share_kc', 'share_odor', 'share_rest', 'cv_odors_median_kc',
                                                      'mean_per_kc', 'mean_1-4', 'mean_5-8', 'mean_9+', 'fixed_thr')}
        q = a1[(a1.site == site) & (a1.kcs == 'responders')]
        summ[f'a1_{site}_responders'] = {c: rng(q, c) for c in ('share_kc', 'share_odor', 'share_rest', 'mean_per_kc')}
    summ['a2'] = {c: rng(a2, c) for c in a2.columns if c != 'wiring' and pd.api.types.is_numeric_dtype(a2[c])}
    summ['a3'] = {q: {c: rng(a3[a3.quantity == q], c) for c in ('over_thr_p10', 'over_thr_p50', 'over_thr_p90', 'over_thr_max',
                                                                 'over_thr_p50_1-4', 'over_thr_p50_5-8', 'over_thr_p50_9+')}
                  for q in ('b', 'c', 'b+c')}
    rr = a3r[(a3r.kind == 'responses') & (a3r.group == 'all')]
    summ['a3_removed'] = {r_: rng(rr[rr['raise'] == r_], 'share_removed_of_above') for r_ in ('b', 'b+c')}
    summ['a4_apl_fits'] = {c: rng(a4[a4.source == 'fit'], c) for c in ('apl_win_cv', 'apl_win_rel_max', 'apl_pre_over_win')}
    summ['a4_weights'] = {c: rng(a4w, c) for c in a4w.columns if c != 'wiring'}
    summ['a5_rows'] = int(len(a5))
    summ['seconds'] = time.time() - t0
    _write_json(out / 'a_summary.json', summ)
    print(f'E95 A in {time.time() - t0:.0f} s', flush=True)
    return summ


# ---------------------------------------------------------------------------
# parts B, C, D

def wiring_wnorm(label: str) -> pd.Series:
    return E83.wiring_parts(label)['wnorm']['wAPLKC']


def run_point(model: StaticModel, v: str, part: str = '', cc: Optional[dict] = None, sub: str = 'runs'
              ) -> Tuple[dict, dict]:
    """One variant on a wiring at its fit (`reading95`), saved with its stage matrices, drive, loss, static inhibition and
    (part D, C6) APL->KC weights under outputs/.../<sub>/. A control '<v>_ctrl' is the setting with its static inhibition
    permuted over the KCs; D0 the bouton path at dose 0.1 without the component. Returns (the run's row, the reading)."""
    w = model.label
    tag = run_tag(v, w)
    pv = parse_variant(v)
    if pv['part'] == REF:
        raise ValueError('the reference is the saved fit; it is not run')
    cc = H.code_commits() if cc is None else cc
    kc_index = model.saved['on'][1].index
    L = kc_losses(w, kc_index)
    rec: dict = {}
    M = static_matrix(pv, L, rec) if pv['static'] is not None else None
    wk, unit = None, np.nan
    if pv['mode'] == 'claw':
        wk = pd.Series(0.0, index=model.wAPLKC.index, name=model.wAPLKC.name)
    elif pv['mode'] == 'kcs':
        unit = d_unit(model, L)
        if pv['size'] > 0:
            wk = d_weights(model.wAPLKC, wiring_wnorm(w), pv['component'], pv['size'], unit)
    res = reading95(model, pv, M, wk, label=tag)
    o = {'on': res['on'], 'off': res['off']}
    st = res['stages']
    odors = E94._odors(st['bouton_win'])
    drive = E94.kc_drive(model.W, st['bouton_win'][odors].to_numpy(dtype=float), kc_index, odors)
    static = None if res['ext'] is None else _frame(res['ext'], res['on'][1])
    save_run95(tag, o, sub, stages=st, drive=drive, loss=_frame(res['gain']['L'], res['on'][1]), static=static,
               weights=wk)
    thr = np.asarray(model.rv.kc.thr, dtype=float).ravel()
    resp_off = E89._responders(model.saved['off'][1])
    p = res['params']
    dose = 0.0 if res['whole'] else float(pv['dose'])
    row = dict(part=part, arm=pv['arm'], tag=tag, variant=v, source='E95', screen=pv['screen'], mode=pv['mode'],
               static=pv['static'], exponent=pv['exponent'], permute=pv['permute'], seed=pv['seed'],
               component=pv['component'], size=pv['size'], dose=dose, control_of=pv['control_of'], wiring=w, host=_host(),
               fit_folder=model.folder, fit_tag=model.tag, factor=res['factor'],
               factor_kind='static' if res['whole'] else 'tied_apl',
               apl_factor=np.nan if res['whole'] else res['factor'], static_factor=res['factor'] if res['whole'] else np.nan,
               at_bound=res['at_bound'], bracket_lo=res['bracket'][0], bracket_hi=res['bracket'][1],
               n_on=E83._n(res['on'][0]), target_on=model.target_on, n_off=E83._n(res['off'][0]), target_off=model.target_off,
               n_sims=res['n_sims'], seconds=res['seconds'], **{k: p[k] for k in E83.PARAM_KEYS},
               aplpn_fit_multiple=E94.aplpn_fit_multiple(w), aplpn_multiple=E94.aplpn_multiple(w, dose),
               wAPLKC_sha256=E89.E80_sha(model.wAPLKC if wk is None else wk), K_fit=model.K, K_run=E83.readout_K(o),
               static_slope=rec.get('slope', np.nan), static_n_clipped=rec.get('n_clipped', np.nan),
               static_min_before_clip=rec.get('min_before_clip', np.nan), component_unit=unit,
               **res['gain']['row'], **cc)
    row.update(E94.stage_extras(model.W, w, st))
    row.update(static_realism(res['ext'], thr, L['n'].to_numpy(), L['Doff'].to_numpy(), res['on'][1].to_numpy(),
                              res['off'][1].to_numpy(), float(p['fixed_thr'])))
    if pv['part'] == 'D':
        w_all = (model.wAPLKC if wk is None else wk).to_numpy(dtype=float) * p['wAPLKC_scale']
        w_fam = model.wAPLKC.to_numpy(dtype=float) * p['wAPLKC_scale']
        off_st = E94.fit_stage(w, 'off')
        row.update(d_realism(model.W.to_numpy(dtype=float), off_st['bouton_win'][odors].to_numpy(dtype=float),
                             st['bouton_win'][odors].to_numpy(dtype=float), st['claw_win'][odors].to_numpy(dtype=float),
                             w_fam, w_all - w_fam, E94.kc_of_claw(model.W), resp_off))
    row.update(E83.realism(o))
    print(f"E95 {tag}: {row['factor_kind']} factor x {res['factor']:.4g}"
          f"{' (at a bound)' if res['at_bound'] else ''}, {row['n_on']} / {model.target_on} responses, {res['n_sims']} sims, "
          f"{res['seconds']:.0f}s; static / threshold p50 {row['static_over_thr_p50']:.3g}; share of the claws' loss at the "
          f"boutons {row.get('bouton_share', np.nan):.3f}; flag {row.get('realism_flag')!r}", flush=True)
    return row, res


RUNS_NAMES = {'B': 'runs_B.csv', 'C': 'runs_C.csv', 'D': 'runs_D.csv', 'E': 'runs_E.csv'}


def _provenance_add(out: Path, part: str, tag: str, prov: dict) -> None:
    path = out / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E95_static_apl.py')
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def points(out: Path, variants: Sequence[str], wirings: Sequence[str], runs_name: str, part: str) -> pd.DataFrame:
    """Variants on wirings, each by `run_point` (resumable: rows in `runs_name`). Needs part 0 to have passed."""
    import gc
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identities must hold before any scored point)')
    cc = _commits()
    path = out / runs_name
    for w in wirings:
        todo = [v for v in variants if run_tag(v, w) not in _done(path)]
        if not todo:
            continue
        model = StaticModel(w, check=False)
        for v in todo:
            row, _ = run_point(model, v, part, cc)
            _append(path, row)
            _provenance_add(out, part, run_tag(v, w), dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **cc))
        del model
        gc.collect()
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def part_points(out: Path, part: str, arms: Optional[Sequence[str]] = None, wirings: Sequence[str] = S_WIRINGS
                ) -> pd.DataFrame:
    """The arms' points (default every arm of the part) on `wirings` (resumable; runs_<part>.csv)."""
    arms = [a for a in ARMS if PART_OF_ARM[a] == part] if not arms else list(arms)
    bad = [a for a in arms if PART_OF_ARM.get(a) != part]
    if bad:
        raise ValueError(f'{bad} are not arms of part {part}')
    vs = [v for a in arms for v in ARMS[a]]
    return points(out, vs, wirings, RUNS_NAMES[part], part)


# ---------------------------------------------------------------------------
# the references and the runs table

def reference_rows() -> pd.DataFrame:
    """The references: the unchanged fits ('ref', source 'saved', every wiring, with E90's stage runs' bouton and claw
    readings) and D0 on S (E94's saved B_sub_k2_d0.1 runs, source 'E94', with E94's rows); realism recomputed from the
    runs."""
    if 'refs' in _CACHE:
        return _CACHE['refs'].copy()
    rows = []
    for w in ALL_WIRINGS:
        f = E83.saved_fit(w)
        r = dict(part='saved', arm=REF, tag=f"{f['folder']}:{f['tag']}", variant=REF, source='saved', screen=REF, mode=REF,
                 dose=1.0, wiring=w, fit_folder=f['folder'], fit_tag=f['tag'], factor=1.0, apl_factor=1.0, at_bound=False,
                 **{k: f['p'][k] for k in E83.PARAM_KEYS}, aplpn_fit_multiple=E94.aplpn_fit_multiple(w),
                 aplpn_multiple=E94.aplpn_fit_multiple(w), K_run=f['K'], **E83.realism(f['saved']))
        try:
            W, _ = E94.wiring_W(w)
            r.update(E94.stage_extras(W, w, E94.fit_stage(w, 'on')))
        except Exception as e:                    # (the bouton readings beside; not needed for the scores)
            r['stage_note'] = repr(e)
        rows.append(r)
    p94 = H.results_dir(E94.EXPERIMENT) / 'runs_B1.csv'
    if p94.exists():
        d = pd.read_csv(p94).drop_duplicates('tag', keep='last')
        for q in d[d.variant == E94_D0_TAG].to_dict('records'):
            o = E94.load_run(q['tag'])
            if o is None or q['wiring'] not in S_WIRINGS:
                continue
            keep = {k: q.get(k) for k in E83.PARAM_KEYS + ('bouton_share', 'bouton_share_removed', 'bouton_zero_share',
                                                           'aplpn_multiple', 'apl_factor', 'n_sims')}
            rows.append(dict(part='E94_B1', arm=D0, tag=q['tag'], variant=D0, source='E94', screen='control', mode='kcs',
                             dose=D_DOSE, size=0.0, wiring=q['wiring'], factor=q.get('apl_factor'), at_bound=_as_bool(q.get('at_bound')),
                             aplpn_fit_multiple=E94.aplpn_fit_multiple(q['wiring']), K_run=E83.readout_K(o), **keep,
                             **E83.realism(o)))
    df = pd.DataFrame(rows)
    _CACHE['refs'] = df
    return df.copy()


def runs_table(out: Path) -> pd.DataFrame:
    """Every run of the readings: this experiment's (runs_B / C / D / E.csv, the latest row per tag) first, then the
    references; with 'flagged', 'at_bound' and the parsed setting."""
    parts = []
    for n in RUNS_NAMES.values():
        if (out / n).exists():
            try:
                parts.append(pd.read_csv(out / n).drop_duplicates('tag', keep='last'))
            except pd.errors.EmptyDataError:
                pass
    parts.append(reference_rows())
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates(['variant', 'wiring'], keep='first').reset_index(drop=True)
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('arm', 'screen', 'mode', 'static', 'exponent', 'permute', 'component', 'size', 'control', 'equal_kind'):
        df[c] = pv[c].to_numpy()
    df['flagged'] = [_flagged(x) for x in df.realism_flag]
    df['at_bound'] = [_as_bool(x) for x in df.at_bound]
    return df


def load_any(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E95' this experiment's, 'saved' the wiring's fit, 'E94' E94's runs."""
    if source == 'saved':
        return E83.saved_fit(wiring)['saved']
    o = {'E95': load_run95, 'E94': E94.load_run}[source](tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


# ---------------------------------------------------------------------------
# the readings

READOUTS = ('cal', 'calK', 'raw', 'zcell', 'zcell_lin')


def zcell_distances_map(o: Dict[str, tuple], amp: str = 'saturating') -> Dict[str, pd.Series]:
    """E73's per-KC reading (the session z-score, graded, noise in expectation) under the saturating map (E93's
    `zcell_distances`, K calibrated on the run's APL-on counts) or the linear map (only the noise and the noise ratio
    calibrated, on the same targets): {'on', 'off'} pair distances (canonical order)."""
    from experiments import E93_missing_glomeruli_input as E93
    if amp == 'saturating':
        return E93.zcell_distances(o)
    from experiments import E73_per_kc_normalization as E73
    on = H.abbrev_columns(o['on'][1].astype(float))[list(H.APL16_ODORS)].to_numpy()
    off = H.abbrev_columns(o['off'][1].astype(float))[list(H.APL16_ODORS)].to_numpy()
    t, q, drift = E73.refs()['targets'], E73.q_values()[0], E73.refs()['drift']
    cal = H.calibrate_session_zscore(on, q=q, amp='linear', target_rate=t['rate'], target_silent_sd=t['silent_sd'], drift=drift)
    kw = dict(q=float(q), noise=float(cal['noise']), amp='linear', noise_ratio=float(cal['noise_ratio']), drift=float(drift))
    return {c: E83._cs(H.zcell_pair_distances(H.session_zscore(x, **kw), 'graded')) for c, x in (('on', on), ('off', off))}


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
                have[ro] = zcell_distances_map(o, 'saturating')
            elif ro == 'zcell_lin':
                have[ro] = zcell_distances_map(o, 'linear')
            else:
                have[ro] = E83.pair_distances(o, ro)
    return have


def readings95(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = READOUTS) -> pd.DataFrame:
    """readings.csv: per variant and unit (each wiring; S and T where all four have a run) and readout (READOUTS) E92's
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
                rows.append(dict(variant=variant, arm=pv['arm'], screen=pv['screen'], mode=pv['mode'], static=pv['static'],
                                 exponent=pv['exponent'], permute=pv['permute'], component=pv['component'], size=pv['size'],
                                 dose=pv['dose'], control=pv['control'], equal_kind=pv['equal_kind'], unit=unit,
                                 n_wirings=len(ws), readout=ro, n_flagged=int(sum(meta[w].flagged for w in ws)),
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


def line95(out: Path, pr: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """line.csv: each population's ester and alcohol groups' APL-on distances relative to the 120 pairs' mean, placed
    against E92's part A line of S (E94's `line94` with this experiment's reference)."""
    from experiments import E92_steps_like_for_like as A
    pr = pd.read_csv(out / 'pairs.csv') if pr is None else pr
    pts = {}
    for (v, unit, ro), g in pr.groupby(['variant', 'unit', 'readout'], sort=False):
        on = pd.Series(g.on.to_numpy(), index=pd.MultiIndex.from_arrays([g.odor_a, g.odor_b]))
        pts[(v, unit, ro)] = A.relative_distances(on)
    rows = []
    for (v, unit, ro), (est, alc) in pts.items():
        L = A.line(ro)['S']
        ref = 0.0
        if unit != 'S' and (REF, unit, ro) in pts:
            e0, a0 = pts[(REF, unit, ro)]
            ref = float(e0 - (L['intercept'] + L['slope'] * a0))
        rows.append(dict(variant=v, unit=unit, readout=ro, ester_rel=est, alcohol_rel=alc, slope=L['slope'],
                         intercept=L['intercept'], scatter=L['scatter'], max_abs_resid=L['max_abs_resid'], ref_resid=ref,
                         **A.place_on_line(est, alc, L, ref)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'line.csv', index=False)
    return df


def _co_job95(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3 (E94's
    `_co_job94` on this experiment's runs)."""
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


def coresponse95(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse.csv: per S setting (all four S runs) the changes of L (co-response) and R (recruitment) under APL of the
    four nearest ester pairs, the ester and the alcohol groups on S's population, against the flies' (E94's
    `coresponse94` on this experiment's runs; jobs cached under outputs/.../coresponse_cache.pkl)."""
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
            for res, j in zip(ex.map(_co_job95, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E95 coresponse {j[1]}', flush=True)
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


def c1_groups(out: Path, pr: Optional[pd.DataFrame] = None, variant: str = 'C1_bck') -> pd.DataFrame:
    """c1_groups.csv: what the feedback adds, read on C1 against the unchanged fits: per readout ('cal', 'raw') and group of
    pairs (E92's part A groups and the 120), S's mean APL-on distance, mean APL effect and their differences (C1 minus the
    fits)."""
    from experiments import E92_steps_like_for_like as A
    pr = pd.read_csv(out / 'pairs.csv') if pr is None else pr
    rows = []
    for ro in ('cal', 'raw'):
        g = {v: pr[(pr.variant == v) & (pr.unit == 'S') & (pr.readout == ro)] for v in (variant, REF)}
        if any(x.empty for x in g.values()):
            continue
        s = {v: x.set_index(['odor_a', 'odor_b']) for v, x in g.items()}
        idx = s[REF].index
        m = A.group_masks(pd.MultiIndex.from_tuples([tuple(sorted(i)) for i in idx]))
        m['all120'] = np.ones(len(idx), dtype=bool)
        for grp, mk in m.items():
            if not mk.any():
                continue
            row = dict(readout=ro, group=grp, n_pairs=int(mk.sum()))
            for c in ('on', 'effect'):
                a, b = s[variant][c].reindex(idx).to_numpy()[mk], s[REF][c].to_numpy()[mk]
                row.update({f'{c}_c1': float(np.mean(a)), f'{c}_fits': float(np.mean(b)), f'{c}_diff': float(np.mean(a - b))})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'c1_groups.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# part E: the choice, the gate, the controls, the equal rule

PARTS_COLS = E92G.PARTS_COLS


def _settings95(rd: pd.DataFrame, screens: Sequence[str] = (REF,) + CHOICE_ARMS) -> pd.DataFrame:
    return rd[(rd.unit == 'S') & (rd.readout == 'cal') & rd.screen.isin(list(screens)) & ~rd.control.astype(bool)]


def _pat(rd: pd.DataFrame, variant: str, unit: str, target: str, col: str = 'pattern', readout: str = 'cal') -> float:
    return E89._pat(rd, variant, unit, target, col, readout)


def choose95(rd: pd.DataFrame, half: str, runs: Optional[pd.DataFrame] = None) -> dict:
    """The choice on S against half `half`: the setting with the lowest pattern part of S's population under the
    calibrated readout among the unchanged fits and CHOICE_ARMS' settings whose four runs are inside the bands and not at
    a bracket's bound; something to confirm if it is a setting and gains CHOICE_MIN_GAIN or more over the unchanged fits.
    The parts of the pattern part beside (the thresholds stay on the pattern part)."""
    g = _settings95(rd)
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
               control=control_of(best.variant) if best.variant != REF else '')
    for c in PARTS_COLS[1:]:
        row[f'chosen_{c}'] = float(best[f'half{half}_{c}'])
        row[f'ref_{c}'] = float(rr[f'half{half}_{c}'])
    return row


def plan_e(rd: pd.DataFrame, runs: Optional[pd.DataFrame] = None) -> Dict[str, dict]:
    return {h: choose95(rd, h, runs) for h in HALVES}


def equal_rule95(rd: pd.DataFrame) -> List[dict]:
    """E94's equal rule in the entry's three kinds (EQUAL_KINDS): a setting of one of those arms is an equal of the
    unchanged fits if, on S under the readout, its pattern part is within EQUAL_PATTERN of theirs against all recordings
    and halves A and B, its four runs are inside the bands and none at a bracket's bound, its raw-count loss within
    EQUAL_RAW and its mean APL effect within EQUAL_MEAN of theirs; 'better' if its pattern part is lower on all three (and
    the rest holds). The other settings of CHOICE_ARMS are listed with kind None."""
    g = _settings95(rd)
    out = []
    if REF not in set(g.variant):
        return out
    me0 = float(g[g.variant == REF].mean_effect.iloc[0])
    for r in g.itertuples():
        if r.screen not in CHOICE_ARMS:
            continue
        kind = r.equal_kind if isinstance(r.equal_kind, str) and r.equal_kind else None
        d = {t: _pat(rd, r.variant, 'S', t) - _pat(rd, REF, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = _pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - _pat(rd, REF, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - me0
        rest = bool(int(r.n_flagged) == 0 and int(r.n_at_bound) == 0 and abs(raw) <= EQUAL_RAW and abs(me) <= EQUAL_MEAN)
        out.append(dict(variant=r.variant, arm=r.arm, screen=r.screen, kind=kind, **{f'd_pattern_{t}': x for t, x in d.items()},
                        d_raw_loss=raw, d_mean_effect=me, n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound),
                        equal=bool(kind is not None and all(abs(x) <= EQUAL_PATTERN for x in d.values()) and rest),
                        better=bool(kind is not None and all(x < 0 for x in d.values()) and rest)))
    return out


def part_e(out: Path) -> dict:
    """E: the choice per half (choice.csv); each chosen setting with something to confirm run on T (runs_E.csv), then its
    control (a static setting's values permuted over the KCs on S and T; part D's D0 on T, S's being E94's saved runs);
    then `score` without the co-response."""
    runs = runs_table(out)
    rd = readings95(out, runs, readouts=('cal', 'raw'))
    plan = plan_e(rd, runs)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    chosen = sorted({c['chosen'] for c in plan.values() if c['confirm']})
    if not chosen:
        print('E95 E: nothing to confirm', flush=True)
    for v in chosen:
        points(out, [v], T_WIRINGS, RUNS_NAMES['E'], 'E')
    for v in chosen:
        c = control_of(v)
        for ws in ((T_WIRINGS,) if c == D0 else (S_WIRINGS, T_WIRINGS)):
            points(out, [c], ws, RUNS_NAMES['E'], 'E')
    return score(out, coresponse=False)


def f_condition(rd: pd.DataFrame, eqs: Optional[Sequence[dict]] = None) -> dict:
    """Whether part F runs (the entry): S's best setting of D1 or D3 (the lowest pattern part of S's population under the
    readout against all recordings) within F_COND of the unchanged fits' UNCHANGED_S['all']; the connectome's own weights
    (`--aplkc-flat connectome`) only if a setting of D2 is an equal at fixed scales (`equal_rule95`)."""
    g = _settings95(rd, ('D1', 'D3'))
    out = dict(cond=F_COND, unchanged=UNCHANGED_S['all'],
               unchanged_here=_pat(rd, REF, 'S', 'z8') if REF in set(rd.variant) else np.nan)
    if g.empty:
        return dict(out, runs=False, note='D1 / D3 not run on S')
    best = g.loc[g.z8_pattern.idxmin()]
    eqs = equal_rule95(rd) if eqs is None else eqs
    d2 = [e['variant'] for e in eqs if e['arm'] == 'D2' and e['equal']]
    return dict(out, best=best.variant, best_pattern=float(best.z8_pattern), above=float(best.z8_pattern) - UNCHANGED_S['all'],
                runs=bool(float(best.z8_pattern) - UNCHANGED_S['all'] <= F_COND), connectome=bool(d2), d2_equals=d2)


def verdict95(out: Path, rd: Optional[pd.DataFrame] = None, runs: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: the choices, the gate in four directions, the equals by kind, S's settings with their parts, part F's
    condition and scores, the flagged runs and those at a bracket's bound, the alcohol pairs' co-response change."""
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table(out) if runs is None else runs
    plan = plan_e(rd, runs)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    dirs = []
    for X in HALVES:
        c = plan[X]
        chosen = c['chosen'] if c['confirm'] else None
        dirs.append(E92G.gate_direction92(rd, runs, X, chosen, REF, [control_of(chosen)] if chosen else []))
    pd.DataFrame(dirs).to_csv(out / 'gate.csv', index=False)
    gate = E92G.gate92(dirs)
    eqs = equal_rule95(rd)
    pd.DataFrame(eqs).to_csv(out / 'equal.csv', index=False)
    gains_s = any(c['confirm'] for c in plan.values())
    v = dict(thresholds=dict(CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=E91.GATE_GAIN, GATE_SINGLES_SHARE=E91.GATE_SINGLES_SHARE,
                             GATE_RAW_SLACK=E91.GATE_RAW_SLACK, GATE_MEAN_SLACK=E91.GATE_MEAN_SLACK, GATE_SQ_ALL=E91.GATE_SQ_ALL,
                             GATE_SQ_SIX=E91.GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO, EQUAL_PATTERN=EQUAL_PATTERN,
                             EQUAL_RAW=EQUAL_RAW, EQUAL_MEAN=EQUAL_MEAN, APL_BRACKET=APL_BRACKET, STATIC_BRACKET=STATIC_BRACKET,
                             F_COND=F_COND),
             readout="cal (the calibrated readout, K set on each run's APL-on counts; E83's readings and gate)",
             reference=REF, choice=plan, gate=gate,
             verdict='passes' if gate['passes'] else ('gains on S, not through the gate' if gains_s else 'level'),
             equals=eqs, equals_by_kind={k: [e['variant'] for e in eqs if e['equal'] and e['kind'] == k] for k in EQUAL_KINDS})
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    v['S'] = {r.variant: dict(arm=r.arm, screen=r.screen, halfA=r.halfA_pattern, halfB=r.halfB_pattern, halfC=r.halfC_pattern,
                              halfD=r.halfD_pattern, all=r.z8_pattern, six=r.z6_pattern, loss=r.z8_loss, mean_effect=r.mean_effect,
                              amplitude=r.z8_amplitude, shape=r.z8_shape, sd_ratio=r.z8_sd_ratio, pearson=r.z8_pearson,
                              rho=r.z8_rho, err_ester=r.z8_err_ester, err_alcohol=r.z8_err_alcohol, n_flagged=int(r.n_flagged),
                              n_at_bound=int(r.n_at_bound)) for r in g.itertuples()}
    v['F_condition'] = f_condition(rd, eqs)
    if (out / 'f_scores.json').exists():
        v['F'] = _read_json(out / 'f_scores.json')
    mine = runs[runs.source == 'E95']
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
    print(f"E95 E: {v['verdict']}; " + '; '.join(f"{h} {plan[h].get('chosen')} (gain {plan[h].get('gain', np.nan):.4f})"
                                                for h in HALVES)
          + f"; equals {v['equals_by_kind']}; F {'runs' if v['F_condition'].get('runs') else 'does not run'}", flush=True)
    return v


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """readings.csv, pairs.csv, line.csv, c1_groups.csv, coresponse.csv (with `coresponse`), choice.csv, gate.csv,
    equal.csv, verdict.json, the figure and comparison_table.md."""
    runs = runs_table(out)
    rd = readings95(out, runs, readouts=READOUTS if zcell else tuple(r for r in READOUTS if not r.startswith('zcell')))
    try:
        line95(out)
    except Exception as e:                       # (E92 part A's line needs its line.json; reported, not fatal)
        print(f'E95 score: line.csv left out ({e!r})', flush=True)
    c1_groups(out)
    if coresponse:
        coresponse95(out, n_workers)
    v = verdict95(out, rd, runs)
    figure(out)
    comparison_table(out)
    return v


# ---------------------------------------------------------------------------
# part F: the family fitted again at the connectome's bouton share with the component

def f_tag(label: str, component: str = 'uniform') -> str:
    """The tag E10_joint_fit.py gives a part F fit (its prefix rules: the wiring, rate0.068_, aplboth_, aplpnhold_, flatu_ /
    flatc_, then k2_f1)."""
    from experiments import E10_joint_fit as J
    a = E91.d_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return (f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.aplpn_hold_tag(1.0)}"
            f"{J.aplkc_flat_tag(dict(kind=component))}k2{F_TAG}")


def f_start_key(label: str, component: str = 'uniform') -> str:
    return f'{label}_f1_{component}'


def f_start_params(out: Path, label: str, component: str = 'uniform') -> dict:
    """Part F's start on a wiring: the absolute parameters of its best point of D1 (D2 for the connectome's component;
    the lowest pattern part of the wiring under the readout against all recordings, inside the bands where one is, not at
    a bracket's bound) and that point's size as 'aplkc_flat', written into f_start_params.json under `f_start_key`."""
    arm = 'D1' if component == 'uniform' else 'D2'
    rd = pd.read_csv(out / 'readings.csv')
    runs = runs_table(out)
    q = runs[(runs.arm == arm) & (runs.wiring == label) & ~runs.at_bound]
    if q.empty:
        raise RuntimeError(f'{label}: no {arm} point inside the bracket (run --part D --arm {arm} on this wiring first)')
    sc = rd[(rd.unit == label) & (rd.readout == 'cal')].set_index('variant').z8_pattern
    q = q.assign(pattern=q.variant.map(sc), flag_order=q.flagged.astype(int)).sort_values(['flag_order', 'pattern'])
    r = q.iloc[0]
    start = {k: float(r[k]) for k in E83.PARAM_KEYS}
    path = out / F_START_FILE
    d = _read_json(path)
    d[f_start_key(label, component)] = dict(start, aplkc_flat=float(r['size']), from_tag=r.tag, pattern=float(r.pattern),
                                            flagged=bool(r.flagged), component_unit=float(r['component_unit']))
    _write_json(path, d)
    return d[f_start_key(label, component)]


def f_command(label: str, component: str = 'uniform', maxfev: int = F_MAXFEV, unit: Optional[float] = None) -> List[str]:
    """The E10_joint_fit.py command of one part F fit: C2's family (`-k 2 --rate-target ctrl_frac --loss apl --apl-readout
    both`) with the APL->bouton scale held at D_DOSE times the wiring's fitted multiple of E3's base (`--aplpn-hold`, as
    E94's D1), the unweighted component's size fitted in F_FLAT_BOUNDS (`--aplkc-flat`, its unit the wiring's
    `component_unit`, from f_start_params.json unless given), into this experiment's folders, from the wiring's best point
    of D1 (`--init-params`, which also gives the size's start)."""
    m = E91.d_member(label)
    F = E94.aplpn_multiple(label, D_DOSE)
    if unit is None:
        unit = float(_read_json(out_dir() / F_START_FILE)[f_start_key(label, component)]['component_unit'])
    return [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
            repr(float(F)), '--aplkc-flat', component, '--aplkc-flat-unit', repr(float(unit)), '--aplkc-flat-max',
            f'{F_FLAT_BOUNDS[1]:g}', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', F_TAG,
            '--maxfev', str(int(maxfev)), '--experiment', EXPERIMENT,
            '--init-params', f'results/{EXPERIMENT}/{F_START_FILE}:{f_start_key(label, component)}'] + m['wiring_args']


def _fits(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def f_allowed(out: Path, pop: str, component: str = 'uniform') -> Tuple[bool, str]:
    """S's fits: on part F's condition (`f_condition`; the connectome's component only on a D2 equal). T's: only if S's
    fitted population is an equal or better by E94's D1 rule (`f_scores`)."""
    if pop == 'S':
        rd = pd.read_csv(out / 'readings.csv') if (out / 'readings.csv').exists() else None
        if rd is None:
            return False, 'not scored (run --part score)'
        c = f_condition(rd)
        ok = bool(c.get('runs')) and (component == 'uniform' or bool(c.get('connectome')))
        return ok, json.dumps({k: c.get(k) for k in ('best', 'best_pattern', 'above', 'runs', 'connectome')})
    s = f_scores(out, write=False).get(f'S_{component}')
    if not s:
        return False, "S's fits not scored"
    return bool(s.get('equal') or s.get('better')), f"S's population: equal {s.get('equal')}, better {s.get('better')}"


def part_f(out: Path, pop: str = 'S', wirings: Optional[Sequence[str]] = None, component: str = 'uniform',
           force: bool = False, dry_run: bool = False) -> List[str]:
    """Part F's fits of population `pop` (or `wirings` of it: a lane), one process after another (`f_command`; logs under
    outputs/.../F/), skipping those in fits.csv; refuses what the entry's rules do not allow unless `force`."""
    import subprocess
    ok, why = f_allowed(out, pop, component)
    print(f'E95 F {pop} {component}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
    if not ok and not force:
        return []
    ws = [w for w in _wirings(pop) if wirings is None or w in wirings]
    if not dry_run:
        _commits()
    logs = local_dir('F')
    tags = []
    for w in ws:
        tag = f_tag(w, component)
        tags.append(tag)
        if tag in set(_fits(out).tag):
            print(f'E95 F {tag}: done', flush=True)
            continue
        if not dry_run:
            f_start_params(out, w, component)
        try:
            cmd = f_command(w, component)
        except KeyError:
            cmd = f_command(w, component, unit=float('nan'))
        print('E95 F:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E95 F {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits(out).tag):
            raise RuntimeError(f'E95 F: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def f_scores(out: Path, write: bool = True) -> dict:
    """Part F's scores: per wiring and population (S, T, the eight) and component, the readings under the readout and on
    raw counts against the unchanged fits of the same wirings, with realism; S by E94's D1 rule (`E94.d1_rule`, beside
    E94's refit at the same bouton share, E94_REFIT_S); the eight against UNCHANGED_8 by the same rule and by E91's part D
    pass rule. f_scores.csv, f_scores.json."""
    from experiments import E72_wiring_population as E72
    fits = _fits(out)
    rows, res = [], {}
    for component in COMPONENTS:
        st = {}
        for w in ALL_WIRINGS:
            try:
                tag = f_tag(w, component)
            except Exception:
                continue
            if tag not in set(fits.tag):
                continue
            o = E72.load_fit(EXPERIMENT, tag)
            if o is None:
                continue
            d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
            d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
            sc, sc0 = E92G.score_row92(d['cal']), E92G.score_row92(d0['cal'])
            raw, raw0 = E83.effect_scores(d['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
            fr = fits[fits.tag == tag].iloc[-1]
            row = dict(component=component, unit=w, tag=tag, **E94._d_row(sc, sc0, raw, raw0),
                       lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                       **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                                 'kc_kc_inh', 'aplpn_hold', 'aplpn_factor_vs_e3', 'aplkc_flat_size',
                                                 'aplkc_flat_unit', 'init_params', 'n_evals')},
                       fit_loss=fr.get('loss'), **E83.realism(o))
            rows.append(row)
            st[w] = (d, d0, row)
        for unit, ws in (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS)):
            if not all(w in st for w in ws):
                continue
            dp = {ro: E83.population([st[w][0][ro] for w in ws]) for ro in ('cal', 'raw')}
            d0 = {ro: E83.population([st[w][1][ro] for w in ws]) for ro in ('cal', 'raw')}
            sc, sc0 = E92G.score_row92(dp['cal']), E92G.score_row92(d0['cal'])
            raw, raw0 = E83.effect_scores(dp['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
            r = dict(E94._d_row(sc, sc0, raw, raw0), n_lower=int(sum(st[w][2]['lower'] for w in ws)),
                     n_flagged=int(sum(_flagged(st[w][2].get('realism_flag')) for w in ws)),
                     sizes=[float(st[w][2].get('aplkc_flat_size', np.nan)) for w in ws])
            if unit == 'S':
                r['equal'], r['better'] = E94.d1_rule(r)
                r['vs_e94_refit'] = float(r['pattern_all'] - E94_REFIT_S)
            if unit == 'all8':
                rest = bool(r['n_flagged'] == 0 and abs(r['raw_loss'] - r['raw_loss_unchanged']) <= E94.D_RAW_SLACK
                            and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
                r['equal'] = bool(abs(r['pattern_all'] - UNCHANGED_8) <= E94.D_EQUAL and rest)
                r['passes'] = bool(r['pattern_all'] <= UNCHANGED_8 - E94.D_PASS_GAIN and r['n_lower'] >= E94.D_WIRINGS_LOWER
                                   and r['n_flagged'] == 0 and r['raw_loss'] <= E94.D_RAW_REF + E94.D_RAW_SLACK
                                   and abs(r['mean_effect'] - FLIES_MEAN) <= E94.D_MEAN_SLACK)
            res[f'{unit}_{component}'] = r
            rows.append(dict(component=component, unit=unit, **{k: x for k, x in r.items() if k != 'sizes'}))
    if write:
        pd.DataFrame(rows).to_csv(out / 'f_scores.csv', index=False)
        _write_json(out / 'f_scores.json', res)
    return res


def f_check(out: Path) -> dict:
    """E10_joint_fit.py with and without --aplkc-flat against a saved fit's first evaluation: E94's D1 fit on the real
    wiring (its command, from its own start) run with --smoke (5 evaluations, into this experiment's results folder, tag
    suffix _e10check), without the option and with it at size 0 (uniform, this wiring's unit); each first evaluation's
    loss terms must equal E94's history row 0 exactly. f_check.json."""
    import subprocess
    hist94 = pd.read_csv(H.results_dir(E94.EXPERIMENT) / f"history_{E94.d1_tag('real')}.csv")
    start = _read_json(H.results_dir(E94.EXPERIMENT) / E94.D_START_FILE)[E94.d_start_key('real')]
    m = E91.d_member('real')
    base = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
            repr(float(E94.aplpn_multiple('real', float(start['dose'])))), '--rate-target', 'ctrl_frac', '--loss', 'apl',
            '--apl-readout', 'both', '--maxfev', '5', '--smoke', '--experiment', EXPERIMENT,
            '--init-params', f'results/{E94.EXPERIMENT}/{E94.D_START_FILE}:{E94.d_start_key("real")}'] + m['wiring_args']
    L = kc_losses('real')
    wn = wiring_wnorm('real')
    unit = component_unit(L['bclaw_mean'], E83.saved_fit('real')['p']['wAPLKC_scale'], L['apl_mean'], float(wn.mean()))
    runs = {'off': base + ['--tag', '_e10check'],
            'flat0': base + ['--tag', '_e10check', '--aplkc-flat', 'uniform', '--aplkc-flat-unit', repr(unit),
                             '--aplkc-flat-init', '0']}
    res = dict(reference=f"{E94.EXPERIMENT}: history_{E94.d1_tag('real')}.csv row 0", unit=unit)
    cols = [c for c in hist94.columns if c not in ('time_s', 'eval', 'Unnamed: 0')]
    for name, cmd in runs.items():
        log = local_dir('F') / f'f_check_{name}.log'
        with open(log, 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E95 f_check {name}: E10_joint_fit.py exited {r.returncode} (log {log})')
        tag = E94.d1_tag('real').replace(E94.D1_TAG, '_e10check')
        if name == 'flat0':
            tag = tag.replace('aplpnhold_', 'aplpnhold_flatu_')
        h = pd.read_csv(out / f'history_{tag}.csv')
        diff = {c: [float(h[c].iloc[0]), float(hist94[c].iloc[0])] for c in cols if c in h.columns
                and not E10_same(h[c].iloc[0], hist94[c].iloc[0])}
        res[name] = dict(tag=tag, n_cols=int(sum(c in h.columns for c in cols)), differs=diff, exact=not diff,
                         flat_first=float(h['aplkc_flat'].iloc[0]) if 'aplkc_flat' in h else None)
    h_off = pd.read_csv(out / f"history_{res['off']['tag']}.csv")
    h_flat = pd.read_csv(out / f"history_{res['flat0']['tag']}.csv")
    both = [c for c in h_off.columns if c in h_flat.columns and c not in ('time_s', 'eval', 'Unnamed: 0')]
    res['flat0_vs_off'] = {c: [float(h_flat[c].iloc[0]), float(h_off[c].iloc[0])] for c in both
                           if not E10_same(h_flat[c].iloc[0], h_off[c].iloc[0])}
    res['passes'] = bool(res['off']['exact'] and res['flat0']['exact'] and res['flat0']['flat_first'] == 0.0
                         and not res['flat0_vs_off'])
    _write_json(out / 'f_check.json', res)
    print(f'E95 f_check: {json.dumps(res, default=str)[:2000]}', flush=True)
    return res


def E10_same(a, b) -> bool:
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    return (np.isnan(fa) and np.isnan(fb)) or fa == fb


def part_g(out: Path) -> None:
    """Part G (only by the entry's rules): on a pass by (E) or as a fitted family of eight, the firing-ceiling family with
    the response rate with APL on re-held, S at the flies' 2 s pulse with both response rates re-held, and S through the
    trial noise as E72's part E ran it (each holds if the population keeps at least half of its gain); on an equal of C1,
    C1 on S at the flies' 2 s pulse with both response rates re-held (E80's part A), on raw counts and under the readout,
    beside E80's readings of the feedback family. Not written: it runs only by those rules."""
    v = _read_json(out / 'verdict.json')
    due = dict(passes=bool(v.get('verdict') == 'passes' or v.get('F', {}).get('all8_uniform', {}).get('passes')),
               c1_equal='C1_bck' in v.get('equals_by_kind', {}).get('whole_static', []))
    print(f'E95 G: {due}', flush=True)
    raise NotImplementedError('E95 part G runs only by the entry\'s rules (a pass, or an equal of C1); it is not written')


# ---------------------------------------------------------------------------
# the steps from the claws to the somata, for a chosen or equal setting

def steps_variants(out: Path) -> List[str]:
    v = _read_json(out / 'verdict.json')
    vs = [c['chosen'] for c in v.get('choice', {}).values() if c.get('confirm')]
    vs += [e['variant'] for e in v.get('equals', []) if e.get('equal')]
    return sorted(dict.fromkeys(vs))


def part_steps(out: Path, variants: Optional[Sequence[str]] = None) -> Optional[pd.DataFrame]:
    """The step from the claws to the somata in wild-type flies for `variants` (default every chosen or equal setting)
    beside C2's family on S (E94's `part_steps` on this experiment's runs): somata_steps.csv, claws_steps.csv, steps.csv."""
    from experiments import E90_model_stages as E90
    from experiments import E92_steps_like_for_like as A
    variants = list(variants) if variants else steps_variants(out)
    if not variants:
        print('E95 steps: no chosen or equal setting', flush=True)
        return None
    t0 = time.time()
    runs = runs_table(out)
    soma, claw = [], []
    i = 0
    for v in variants:
        for w in S_WIRINGS:
            q = runs[(runs.variant == v) & (runs.wiring == w)]
            if not len(q):
                raise RuntimeError(f'{v} has no run on {w}')
            r = q.iloc[0]
            o = load_any(r.source, r.tag, w)
            soma.append(A.somata_pairs(o['on'][1], o['off'][1]).assign(family=v, wiring=w, source=r.source, tag=r.tag))
            st_on = {'E95': load_stage95, 'E94': E94.load_stage94}.get(r.source, lambda t: dict(E94.fit_stage(w, 'on')))(r.tag)
            if st_on is None:
                raise RuntimeError(f'{r.tag}: no stage matrices')
            if 'spike_counts' not in st_on:
                # load_stage95 saves the stages only; E90's reading also needs the run's APL-on responses and counts
                st_on['responses'], st_on['spike_counts'] = o['on']
            W, _ = E94.wiring_W(w)
            meta = {k: np.asarray(x) for k, x in E92G.units_meta(W).items()}
            rr = E90.read_wiring(w, st_on, dict(E94.fit_stage(w, 'off')), meta,
                                 np.random.default_rng(E90.NOISE_SEED + E94.STEPS_NOISE_SEED + 100 + i))
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
    print(f'E95 steps in {time.time() - t0:.0f} s', flush=True)
    return rows


# ---------------------------------------------------------------------------
# the archive, the figure and the table

def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike counts,
    the harness's primary model RDM): one row per run of this experiment with its absolute parameters, the static
    inhibition's form and size, the component, realism; part F's fits' distances, which E10_joint_fit.py wrote, are kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E95']
    rows, parts = [], []
    for _, r in runs.iterrows():
        pv = parse_variant(r.variant)
        kind = MODE_NAMES[pv['mode']]
        if pv['static']:
            kind += f"; static inhibition {STATIC_NAMES[pv['static']]}" + (f" (exponent {pv['exponent']:g})" if pv['static'] == 'nke_b'
                                                                           else '')
            if pv['permute']:
                kind += f"; permuted {'within a claw count' if pv['permute'] == 'claw' else 'over all KCs'} (seed {pv['seed']})"
        if pv['part'] == 'D':
            kind += (f"; component {pv['component']}, size {pv['size']:g} (unit {r.get('component_unit'):.6g}), "
                     f"APL->bouton scale x {pv['dose']:g} of the fit's")
        kind += "; the APL-off run the fit's"
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=r.arm, kind=kind, config=r.variant, mode=r['mode'],
                         static=r.static, exponent=r.exponent, permute=r.permute, seed=r.seed, component=r.component,
                         size=r['size'], dose=r.dose, wiring=r.wiring, ref_experiment=r.fit_folder, ref_fit=r.fit_tag,
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         factor=r.factor, factor_kind=r.factor_kind, aplpn_multiple=r.aplpn_multiple,
                         component_unit=r.get('component_unit'), static_mean=r.get('static_mean'),
                         static_over_thr_p50=r.get('static_over_thr_p50'), rate_on=r.rate_on, rate_off=r.rate_off,
                         max_count_on=r.max_count_on, max_count_off=r.max_count_off, wAPLKC_sha256=r.wAPLKC_sha256,
                         control_of=r.control_of, bouton_share=r.get('bouton_share'), realism_flag=r.realism_flag,
                         host=r.get('host'), subdir='runs'))
        o = load_run95(r.tag)
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
    """E95_settings.png: S's pattern part under the readout against all recordings and halves A and B per setting (the
    unchanged fits dashed; hollow: a setting with a run outside the bands or at a bracket's bound), parts B, C and D side by
    side."""
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
    parts = [p for p in ('B', 'C', 'D') if g.arm.astype(str).str.startswith(p).any() or (p == 'D' and (g.variant == D0).any())]
    fig, axes = plt.subplots(1, max(len(parts), 1), figsize=(5.5 * max(len(parts), 1), 4.8), squeeze=False)
    cols = {'halfA': '#d1495b', 'halfB': '#2e86ab', 'z8': '0.2'}
    for ax, part in zip(axes[0], parts):
        q = g[g.arm.astype(str).str.startswith(part) | ((g.variant == D0) & (part == 'D'))]
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
        ax.set_title(f"part {part}, S's population, calibrated readout\n(dashed: the unchanged fits; hollow: outside the bands "
                     'or at a bound)', fontsize=8)
        ax.legend(fontsize=7)
    fig.tight_layout()
    path = out / 'E95_settings.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _md(df: pd.DataFrame, fmt: str = '.4f') -> str:
    return E94._md(df, fmt)


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    L = ["# E95: APL's inhibition as a static threshold raise per KC", '',
         f'Generated by experiments/E95_static_apl.py at {H.code_commits().get("queue_commit", "")[:7]}.', '']
    if (out / 'identity.json').exists():
        s = _read_json(out / 'identity.json')
        L += ['## Part 0: the identities', '', f"passes: {s.get('passes')}", '']
        for k, x in s.items():
            if not isinstance(x, dict) or k == 'commits':
                continue
            keep = {kk: x.get(kk) for kk in ('exact', 'passes', 'same_matrix', 'same_weights', 'changed_columns', 'column',
                                              'odor', 'order_is_identity', 'apl_factor', 'apl_factor_e94', 'factor_rel_diff',
                                              'ext_reset') if kk in x}
            for c in ('on', 'off', 'counts'):
                if isinstance(x.get(c), dict):
                    keep[f'{c}_counts_differ'] = x[c].get('counts_differ')
            L.append(f'- {k}: {json.dumps(keep, default=str)}')
        L.append('')
    for name, title, fmt in (('a1_split.csv', 'A1: the split of the loss per KC and odor (boutons L^b, claws L^c)', '.3f'),
                             ('a2_per_kc.csv', 'A2: b_k against the claw count; c_k against the breadth factor', '.3f'),
                             ('a3_sizes.csv', 'A3: the sizes against each KC\'s threshold', '.3f'),
                             ('a4_apl.csv', 'A4: APL\'s activity by odor and before the odor', '.3f'),
                             ('a4_weights.csv', 'A4: the connectome\'s per-claw APL->KC weights', '.3g'),
                             ('a5_archive.csv', 'A5: the archive\'s static forms', '.4f')):
        if (out / name).exists():
            d = pd.read_csv(out / name)
            if name == 'a3_sizes.csv' and (out / 'a3_removed.csv').exists():
                rr = pd.read_csv(out / 'a3_removed.csv')
                rr = rr[rr.group.isin(['all', 'ester', 'alcohol', 'all120'])]
                L += [f'## {title}', '', _md(d, fmt), '', '### A3: what a threshold raise would remove (window-mean drive; a '
                      'count-level description)', '', _md(rr, fmt), '']
                continue
            L += [f'## {title}', '', _md(d, fmt), '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        runs = runs_table(out)
        g = rd[(rd.unit == 'S') & (rd.readout == 'cal')].sort_values(['arm', 'variant'])
        graw = rd[(rd.unit == 'S') & (rd.readout == 'raw')].set_index('variant')
        rr = runs.groupby('variant')
        L += ["## S's population (calibrated readout; the parts of the pattern part against all recordings)", '',
              '| setting | pattern, half A | half B | all | six | amplitude | shape | spread ratio | Pearson r | ρ | loss | '
              'mean APL effect | raw-count loss | ester / alcohol error | factor | static / threshold p50 | share at the '
              'boutons | flagged / at bound |', '|---|' + '---:|' * 17]

        def rng(q, c, spec='.3g'):
            if not len(q) or c not in q or not pd.to_numeric(q[c], errors='coerce').notna().any():
                return '–'
            v = pd.to_numeric(q[c], errors='coerce')
            return f"{v.min():{spec}}-{v.max():{spec}}"
        for r in g.itertuples():
            q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()
            L.append(f'| {r.variant} | {r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | {r.z6_pattern:.4f} | '
                     f'{r.z8_amplitude:.4f} | {r.z8_shape:.4f} | {r.z8_sd_ratio:.2f} | {r.z8_pearson:.2f} | {r.z8_rho:.2f} | '
                     f'{r.z8_loss:.4f} | {r.mean_effect:+.3f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | '
                     f'{r.z8_err_ester:+.3f} / {r.z8_err_alcohol:+.3f} | {rng(q, "factor")} | '
                     f'{rng(q, "static_over_thr_p50", ".3f")} | {rng(q, "bouton_share", ".3f")} | {r.n_flagged} / {r.n_at_bound} |')
        L.append('')
        mine = runs[runs.source == 'E95']
        if len(mine):
            L += ['## Realism of the runs (rule of 2026-11-10 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | '
                  'K | static / threshold p10 / p50 / p90 / max | silent (static > largest drive) | silent on / off | '
                  'share at the boutons / component / weighted | APL->KC weight p90/p10 over responders | flag | flag '
                  '(earlier rule) |', '|---|' + '---:|' * 12 + '---|---|']
            for r in mine.itertuples():
                sp = ' / '.join(E94._f(getattr(r, f'static_over_thr_{n}', np.nan), '.3f') for n in ('p10', 'p50', 'p90', 'max'))
                sh = ' / '.join(E94._f(getattr(r, c, np.nan), '.3f') for c in ('bouton_share', 'comp_share', 'weighted_share'))
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | '
                         f'{r.opk_on_tv:.3f} / {r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | '
                         f'{E94._f(getattr(r, "K_run", np.nan), ".3g")} | {sp} | {E94._f(getattr(r, "static_over_max_drive", np.nan), ".3f")} | '
                         f'{E94._f(getattr(r, "silent_share_on", np.nan), ".3f")} / {E94._f(getattr(r, "silent_share_off", np.nan), ".3f")} | '
                         f'{sh} | {E94._f(getattr(r, "wkc_all_resp_p90_over_p10", np.nan), ".3g")} | '
                         f'{r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
            if (mine.arm.astype(str).str.startswith('D')).any():
                L += ['Part D: the unweighted component puts an APL->KC weight on claws for which the tables hold no APL synapse '
                      '(30 % of the claws on the real wiring; a4_weights.csv, zero_share).', '']
    for name in ('c1_groups.csv', 'line.csv', 'choice.csv', 'gate.csv', 'equal.csv', 'f_scores.csv', 'steps.csv'):
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
        short = {k: vc.get(k) for k in ('verdict', 'gate', 'equals_by_kind', 'F_condition', 'F', 'flagged_runs', 'at_bound_runs')}
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json's dated records: the commits, olfsysm's build, the model's provenance with the lab's input
    for the simulation parts, the settings and `extra` (the point lists' per-point records are in its 'points')."""
    import olfsysm
    path = out / f'provenance_{part}.json'
    old = _read_json(path)
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('0', 'B', 'C', 'D', 'E', 'F'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    old.update(script='experiments/E95_static_apl.py', part=part, runs=runs,
               settings=dict(ARMS={k: list(v) for k, v in ARMS.items()}, APL_BRACKET=APL_BRACKET, STATIC_BRACKET=STATIC_BRACKET,
                             APL_STEP=APL_STEP, B3_EXPONENTS=B3_EXPONENTS, PERM_SEED=PERM_SEED, D_DOSE=D_DOSE, D_SIZES=D_SIZES,
                             D3_SIZE=D3_SIZE, CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, CONTROL_RATIO=CONTROL_RATIO,
                             EQUAL=dict(PATTERN=EQUAL_PATTERN, RAW=EQUAL_RAW, MEAN=EQUAL_MEAN), EQUAL_KINDS=EQUAL_KINDS,
                             F_COND=F_COND, F_FLAT_BOUNDS=F_FLAT_BOUNDS, F_MAXFEV=F_MAXFEV, MODE_NAMES=MODE_NAMES,
                             STATIC_NAMES=STATIC_NAMES))
    _write_json(path, old)


def main() -> None:
    import warnings
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'A1', 'A2', 'A3', 'A4', 'A5', 'B', 'C', 'D', 'score', 'E',
                                                      'Fcheck', 'F', 'Fscore', 'G', 'steps', 'variants', 'table'))
    ap.add_argument('--arm', default=None, help='B: B1..B5, C: C1..C6, D: D1..D3 (comma-separated; default all of the part)')
    ap.add_argument('--wirings', default=None, help="'S', 'T', 'all8' or a comma-separated list (0: all8; B / C / D: S)")
    ap.add_argument('--pop', default='S', choices=('S', 'T'), help='F: which fits')
    ap.add_argument('--component', default='uniform', choices=COMPONENTS, help="F: the component ('connectome' only on "
                    "its rule)")
    ap.add_argument('--force', action='store_true', help="F: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='F: print the commands only')
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
    elif a.part in ('A1', 'A2', 'A3', 'A4', 'A5'):
        {'A1': a1_split, 'A2': a2_per_kc, 'A3': a3_sizes, 'A4': a4_apl_weights, 'A5': a5_archive}[a.part](out)
    elif a.part in ('B', 'C', 'D'):
        ws = _wirings(a.wirings or 'S')
        arms = [x for x in a.arm.split(',') if x] if a.arm else None
        provenance(out, a.part, wirings=list(ws), arms=arms)
        part_points(out, a.part, arms, ws)
    elif a.part == 'score':
        score(out, a.workers, coresponse=not a.no_co, zcell=not a.no_zcell)
    elif a.part == 'E':
        provenance(out, 'E')
        part_e(out)
    elif a.part == 'Fcheck':
        provenance(out, 'Fcheck')
        f_check(out)
    elif a.part == 'F':
        if not a.dry_run:
            provenance(out, 'F', pop=a.pop, wirings=a.wirings, component=a.component)
        part_f(out, a.pop, None if a.wirings is None else _wirings(a.wirings), a.component, force=a.force,
               dry_run=a.dry_run)
    elif a.part == 'Fscore':
        f_scores(out)
        comparison_table(out)
    elif a.part == 'G':
        part_g(out)
    elif a.part == 'steps':
        part_steps(out, None if a.variants is None else [x for x in a.variants.split(',') if x])
        comparison_table(out)
    elif a.part == 'variants':
        write_variants(out)
    elif a.part == 'table':
        figure(out)
        comparison_table(out)
    print(f'E95 {a.part} done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
