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

Deviations from the entry, forced by the code: (filled in as the parts are written)

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
    E89.set_weights(model, wk)
    model.dose = 1.0 if whole else float(pv['dose'])
    model.ext = None if M is None else _as_kc_odor(M, model)
    E94.set_form94(model.mp, 'sub')
    try:
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
                   on_counts: Optional[np.ndarray] = None, off_counts: Optional[np.ndarray] = None) -> dict:
    """The static inhibition (KCs x odors, as run; None: none) per KC (its mean over the odors) as a multiple of the KC's
    somatic threshold `thr` (fixed_thr + 2 x its spontaneous input): 10th percentile, median, 90th, largest, and the median
    by claw class; 'static_over_max_drive' the share of KCs whose static inhibition exceeds their largest summed claw drive
    without APL (`d_off`, window means; any odor); 'silent_by_mean_drive' the share whose largest drive without APL minus
    the static inhibition stays below the threshold for every odor; with the runs' counts, the shares of KCs that respond
    to no odor with APL on and off (two or more spikes)."""
    row = {}
    thr = np.asarray(thr, dtype=float).ravel()
    n = np.asarray(n_claws, dtype=float).ravel()
    D = np.asarray(d_off, dtype=float)
    if M is None:
        row.update(_q([], 'static_over_thr'))
        row.update({f'static_over_thr_p50_{g}': np.nan for g, _, _ in CLAW_GROUPS})
        row.update(static_over_max_drive=np.nan, silent_by_mean_drive=np.nan, static_mean=np.nan, static_max=np.nan)
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
            # one odor: the column of the odor whose olfsysm index differs most from its output column
            order = np.asarray(m.order, dtype=int)
            j = int(np.argmax(np.abs(order - np.arange(n_od)))) if (order != np.arange(n_od)).any() else 0
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
    10th percentile, median, 90th, largest, and the median by claw class) and a3_removed.csv (per wiring, raise (b,
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
        for q, v in (('b', b), ('c', c), ('b+c', b + c)):
            r = v / thr
            row = dict(wiring=w, quantity=q, mean_abs=float(v.mean()), thr_median=float(np.median(thr)), **_q(r, 'over_thr'))
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
                              res['off'][1].to_numpy()))
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
# the readings, part E

def runs_table(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def readings95(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = ()) -> pd.DataFrame:
    raise NotImplementedError


def choose95(rd: pd.DataFrame, half: str, runs: Optional[pd.DataFrame] = None) -> dict:
    """The choice on S against half `half`: the lowest pattern part of S's population under the readout among CHOICE_ARMS'
    settings (inside the bands on all four wirings, not at a bracket's bound) and the unchanged fits; something to confirm
    if a setting gains CHOICE_MIN_GAIN or more over the unchanged fits."""
    raise NotImplementedError


def equal_rule95(rd: pd.DataFrame) -> List[dict]:
    """E94's rule in the three kinds of EQUAL_KINDS."""
    raise NotImplementedError


def part_e(out: Path) -> dict:
    raise NotImplementedError


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    raise NotImplementedError


# ---------------------------------------------------------------------------
# part F, part G

def f_condition(rd: pd.DataFrame) -> dict:
    """Whether part F runs: S's best setting of D1 or D3 within F_COND of UNCHANGED_S['all'] against all recordings;
    D2's (`--aplkc-flat connectome`) only if a D2 setting is an equal at fixed scales."""
    raise NotImplementedError


def f_command(label: str, component: str = 'uniform', maxfev: int = F_MAXFEV) -> List[str]:
    raise NotImplementedError


def part_f(out: Path, pop: str = 'S', wirings: Optional[Sequence[str]] = None, component: str = 'uniform',
           force: bool = False, dry_run: bool = False) -> List[str]:
    raise NotImplementedError


def f_scores(out: Path, write: bool = True) -> dict:
    raise NotImplementedError


def part_g(out: Path) -> None:
    """Part G (only by the entry's rules): on a pass by (E) or as a fitted family of eight, the firing-ceiling family with
    the response rate with APL on re-held, S at the flies' 2 s pulse with both response rates re-held, and S through the
    trial noise as E72's part E ran it; on an equal of C1, C1 on S at the flies' 2 s pulse with both response rates
    re-held (E80's part A). Not written: it runs only by those rules."""
    raise NotImplementedError('E95 part G runs only by the entry\'s rules (a pass, or an equal of C1); it is not written')


# ---------------------------------------------------------------------------
# the archive, the table

def write_variants(out: Path) -> pd.DataFrame:
    raise NotImplementedError


def comparison_table(out: Path) -> str:
    raise NotImplementedError


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'B', 'C', 'D', 'score', 'E', 'F', 'Fscore', 'G', 'steps',
                                                      'variants', 'table'))
    ap.add_argument('--arm', default=None, help='B: B1..B5, C: C1..C6, D: D1..D3 (comma-separated; default all of the part)')
    ap.add_argument('--wirings', default=None, help="'S', 'T', 'all8' or a comma-separated list (0: all8; B / C / D: S)")
    a = ap.parse_args()
    raise NotImplementedError(f'E95 part {a.part}: not written yet')


if __name__ == '__main__':
    main()
