#!/usr/bin/env python3
"""E131: MB-C1's two cells beside the inhibition among KCs: a γ-directed feedback that the block spares, at its own
scale next to the KC-KC matrix, in half of the matrix's place, and with a time constant and a threshold of its own.

The entry is item E131 of experiment_priorities.md (proposed by viki 2026-11-23 16:10, moved into the list under the
idle-host rule; planned by a Fable subagent, claude-fable-5-2, spawned from viki 15:53). "The family" is C2's breadth-2
APL-effect fit (E67) on E72 part B's eight wirings, each at its own fit (E83's `saved_fit`); S is real, P1, P2, M1 and T
is P3, P4, M2, M3. "KC-KC at zero" is E110's part A. MB-C1 is olfsysm fa3543d's option: the hemibrain's two MB-C1_R
cells as feedback cells of their own (`fit_mb_model(mbc1_inh=, mbc1_tau=, mbc1_thr=, mbc1_pooled=)`, the simulator's
'wMBC1_scale'), alongside the KC-KC matrix, each with its own scale, present with APL on and with APL off.

The protocol (design item 1; nothing fitted). Each run is one of E83's simulators of the wiring at its fit, built with
the arm's MB-C1 options (`build_runner`), both response rates re-held to the family's by E77_input_screen.hold_rates
(the threshold bisected with APL off, then the tied APL scales with APL on). The KC-KC scale ('wKCKC_scale') and
MB-C1's ('wMBC1_scale') are simulator parameters, so one simulator serves every arm with the plain cells.

The arms (tags S1_<arm>_<wiring>; doses m are multiples of the wiring's fitted KC-KC scale s):
- B (24 runs), beside the matrix: the matrix at s, MB-C1 at m 0.25, 1 and 4 (B25, B100, B400).
- R (16 runs), in half of the matrix's place: the matrix at s / 2 with MB-C1 at m 0.5 (R50), and the matrix at s / 2
  alone as the reference (H50).
- X (8 runs), E125's split with a named survivor: the matrix at s / 2 with APL on and at 0 with APL off, MB-C1 at m 0.5
  in both (X50); read against E125's pooled neuron (PO50) and neuron per type (TY50) at the same d.
- The cell's own properties (16 runs, 8 more if gated), beside the matrix at m 1: late (L100: MB-C1's time constant
  1 s, its scale step 0's), strong drive only (D100: MB-C1's threshold step 0's), pooled (P100: the two cells as one,
  only if step 0's gate opens).
- References, read and not run: the family; KC-KC at zero (E110's part A, 'E110_A'); E125's PO50 and TY50.
- Not written: step 2 (the fit) and the permuted control; they run only if a gate opens (`--part score` reports it),
  and their code is written then.

Readings the entry leaves open, fixed here before any run:
1. The dose. m is MB-C1's scale as a multiple of the wiring's fitted KC-KC scale: wMBC1_scale = m x s (`mbc1_scale`),
   also in R and X (R50's MB-C1 at 0.5 s beside the matrix at 0.5 s). The loader's weights (pyolfsysm's
   `connectome_MBC1_cell_weights`, normalize=True) give each cell input weights summing to 1 and an equivalent KC x KC
   matrix (cells' output x input, summed over the cells) with mean row sum 1, as the hemibrain matrix has; so for
   uniform KC activity a the cells' total output is m s a N, the matrix's s a N, and m 1 delivers the same total. Step 0
   checks this on the model's KCs (`dose_check`: both totals N, each cell's input weights summing to 1, the equivalent
   matrix equal to `connectome_MBC1_KCKC_weights(pooled=False)`).
2. H50 is the plain cells' simulator at wMBC1_scale 0, which check (ii) shows equals the option off exactly.
3. X50's APL-off condition has wKCKC_scale 0 (hold_rates' off_extra, as E123's HBab), MB-C1 at 0.5 s in both.
4. A cell's activity (step 0, `cell_activity`): olfsysm's recursion m(t) = (1 - dt / tau) m(t - 1) + w_in @ spikes(t - 1)
   / tau, applied to the family's saved spike steps (E121 step 0's spike times, checked to add up to the saved counts),
   in Hz. Per cell and odor: its mean over the pulse (the steps with the odor on, 0 to 0.5 s) and its peak over the
   trial.
5. D100's threshold, per wiring: the largest activity of either cell over the 16 odors and every time step with APL on
   in the family's run (tau 50 ms), so that the cells are silent wherever the family's APL is intact. Check (iv) runs
   APL on at the family's parameters with MB-C1 at m 1 and this threshold and must reproduce the family's APL-on run
   exactly, else D100 is refused.
6. L100's scale, per wiring: s times the ratio of the cells' summed output over the pulse with APL off in the family's
   run at 50 ms to the same at 1 s (each cell's pulse mean weighted by its total output weight, summed over the cells
   and odors; `late_ratio`), so that the late cells' mean output over the pulse with APL off equals the plain cells'.
7. The pooled arm's gate: the Pearson r over the 16 odors of the two cells' pulse means, per wiring and condition (APL
   on and off); P100 runs if any of the 16 values is below POOLED_GATE_R (0.95).
8. A corner: the threshold or the APL scale at a bound of its bisection, or either response count more than RATE_TOL
   (0.2 %) from its target (E129's `held_ok`). An arm with corners on more than two wirings gets no class.
9. The classes (calibrated pattern part, S's and T's populations apart, against all recordings and each half A to D):
   - "below the family": lower than the family by WITHIN (0.003) or more on both populations and on 6 of 8 wirings, on
     every target; the mean APL effect under the readout within MEAN_TOL (0.02) of the flies' (E108's 'all' target,
     +0.178) on both populations; the raw-count loss against all recordings not higher than the family's by more than
     RAW_LOSS_TOL (0.005) on both populations; every run in the four bands;
   - "equal": E123's (`E123.classify`, one-sided);
   - "costs": the pattern part against all recordings higher than the family's by more than WITHIN on both
     populations;
   - "neither" otherwise; in that order of precedence; a corner arm "corner (no class)".
   - Beside the class: "keeps the gain" (E123's flag, read for R50 and X50), "stands in" (R50: E123's equal flag and
     R50 lower than H50 by E126's "differ" rule, E129's `differ`).
10. The outcomes are flags for the write-up: (a) every B dose costs, R50 does not stand in, X50 does not keep the gain;
    (b) B25 equal and every B dose above the largest tolerated one (the largest dose with every dose up to it equal)
    costs; (c) R50 stands in or X50 keeps the gain; (d) a B dose below the family (the fit's half, on T, is step 2's);
    (e) L100, D100 or P100 against B100 by the "differ" rule; (f) a corner, or a run outside the bands. (f) is a flag:
    bands enter the classes only where the classes name them ("equal", "below"), so a dose that leaves a band can
    still cost.
11. Step 2's gate (reported, not run): a B dose "below the family", R50 "stands in" or X50 "keeps the gain", none a
    corner.
12. Realism of MB-C1 per run (`cell_realism`, from the held point simulated again with its spike times and APL's
    trace, which must reproduce the held counts): each cell's activity in Hz (pulse mean over the odors, largest peak)
    with APL on and off; MB-C1's weight per output synapse against a KC->KC synapse's; and with APL on, the share of
    the inhibition delivered to γ KCs over the pulse that is MB-C1's against MB-C1's and APL's together (APL's: each
    KC's APL->claw weights times APL's trace), beside MB-C1's share of the inhibitory synapses onto γ against APL's
    (step 0's tables).
13. MB-C1's output synapses are those of Pratyush's MBC1>KC table (MB-C1 presynaptic; 12,135 onto the model's KCs) and
    its input those of the KC>MBC1 table (4,716). The entry's "4,716 output synapses", its dose of about 0.16 and its
    γ share of about 0.3 against APL were computed from the input table; step 0 recomputes all three from the output
    table (`dose_arithmetic`: 0.40, and 0.54 against APL's synapses on the claws, 0.42 against all of APL's in the
    calyx), and reading 12 compares with step 0's numbers. The doses m are multiples of the fitted KC-KC scale and do
    not depend on this.

Parts (CLI `--part`):
- check (simulations, `--wirings`): on every wiring, at the family's saved parameters, (i) the option off and (ii) the
  option on at scale 0 reproduce the saved family run (APL on and off) with 0 differing spike counts; on 'real', (iii)
  the option at m 1 changes the counts; once step 0 is done, (iv) D100's threshold check. check_<host>_<lane>.json;
  the runs refuse until every wiring has passed (i) and (ii).
- step0 (no arm simulated; the hemibrain wiring's claw index, for APL's synapse counts, comes from E83's
  `wiring_parts`, which runs the bouton base once): step0_synapses.csv (each cell's synapses in and out by KC type,
  with APL's and the KC->KC synapses beside), step0_activity.csv (each cell's activity per odor, wiring and
  condition), step0_gamma_input.csv (each γ KC's MB-C1 input at m 1 beside its KC-KC input), step0.json (the two fixed
  numbers per wiring, the cells' agreement and the pooled gate, the dose arithmetic, the dose check),
  provenance_step0.json.
- runs (simulations, resumable, per lane `--wirings`, `--arms`; a point claimed by mkdir of
  results/E131_mbc1_cells/claims/<tag>/, skipped when done or claimed, `--reclaim` to take over a dead process's): rows
  in runs_<host>_<lane>.csv, provenance (`H.model_provenance` per point) in provenance_<host>_<lane>.json, runs saved
  locally (outputs/E131_mbc1_cells/step1/) and in the model-output archive (data/model_outputs/E131_mbc1_cells/step1/).
- score (no simulation): readings.csv, runs.csv, verdict.json (classes, comparisons, outcomes, step 2's gate),
  comparison_table.md, provenance_score.json.
- variants: variants.csv and model_pair_distances_long.csv (raw spike counts) for the archive.
- plan: the points a lane would run.

The holds (CLAUDE.md, "Connectome data"): only the model's hemibrain wiring and E72's seven re-draws of it, the
loader's hemibrain tables (KC->KC, APL<->KC, KC types) and Pratyush's hemibrain MB-C1 tables as pyolfsysm reads them;
the eight wirings keep the hemibrain's 1732 KCs, so MB-C1's weights are the same on each. No male CNS table, no olfsysm
change; nothing is published.

Usage (from the repository root; 16 threads, under nice; step 0 needs E121 step 0's spike times, on viki; the check
runs again after step 0 for (iv); one lane for all eight wirings, or S's and T's wirings apart on two hosts):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part check
    ./venv/bin/python -u experiments/E131_mbc1_cells.py --part step0
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part check
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part runs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part runs --wirings real P1 P2 M1
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u experiments/E131_mbc1_cells.py --part runs --wirings P3 P4 M2 M3
    ./venv/bin/python -u experiments/E131_mbc1_cells.py --part score && ./venv/bin/python -u experiments/E131_mbc1_cells.py --part variants
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

EXPERIMENT = 'E131_mbc1_cells'
SCRIPT = 'experiments/E131_mbc1_cells.py'
S_WIRINGS, T_WIRINGS, WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS
POPS = ('S4', 'T4')
POP_WIRINGS = {'S4': S_WIRINGS, 'T4': T_WIRINGS, 'all8': WIRINGS}
TYPES = ('g', "a'b'", 'ab')
TYPE_COL = {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}
CONDS = ('on', 'off')
# E121 step 0 names the hemibrain wiring 'S'
E121_LABEL = {'real': 'S'}
# arm -> group, KC-KC scale with APL on / off as multiples of s, MB-C1 dose m, MB-C1 form
ARMS: Dict[str, dict] = {
    'B25': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=0.25, form='plain'),
    'B100': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=1.0, form='plain'),
    'B400': dict(group='B', kckc_on=1.0, kckc_off=1.0, m=4.0, form='plain'),
    'R50': dict(group='R', kckc_on=0.5, kckc_off=0.5, m=0.5, form='plain'),
    'H50': dict(group='R', kckc_on=0.5, kckc_off=0.5, m=0.0, form='plain'),
    'X50': dict(group='X', kckc_on=0.5, kckc_off=0.0, m=0.5, form='plain'),
    'L100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='late'),
    'D100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='strong'),
    'P100': dict(group='forms', kckc_on=1.0, kckc_off=1.0, m=1.0, form='pooled'),
}
GROUPS = ('B', 'R', 'X', 'forms')
B_ARMS = ('B25', 'B100', 'B400')
FORM_ARMS = ('L100', 'D100', 'P100')
STEP0_ARMS = ('L100', 'D100', 'P100')         # arms that need step 0's numbers or gate
ARM_LABELS = {
    'B25': "the matrix at the fit's KC-KC scale s, MB-C1's two cells at 0.25 s, both conditions",
    'B100': "the matrix at s, MB-C1's two cells at s, both conditions",
    'B400': "the matrix at s, MB-C1's two cells at 4 s, both conditions",
    'R50': 'the matrix at s / 2, MB-C1 at s / 2, both conditions',
    'H50': 'the matrix at s / 2 alone, both conditions (MB-C1 at scale 0)',
    'X50': 'the matrix at s / 2 with APL on and 0 with APL off, MB-C1 at s / 2 in both',
    'L100': "the matrix at s, MB-C1 late (time constant 1 s) at step 0's scale, both conditions",
    'D100': "the matrix at s, MB-C1 at s with step 0's threshold (strong drive only), both conditions",
    'P100': 'the matrix at s, the two MB-C1 cells pooled as one at s, both conditions',
}
E125_READ = ('PO50', 'TY50')
REFS = ('family', 'E110_A') + E125_READ
LATE_TAU_S = 1.0
PLAIN_TAU_S = 0.05
# olfsysm's time grid of the family's runs (checked against the simulator in --part check and against E121's spikes)
TIME_PRE_START, TIME_DT, TIME_END = -2.0, 0.0005, 0.75
STIM_START, STIM_END = 0.0, 0.5
POOLED_GATE_R = 0.95
RATE_TOL = 0.002
WITHIN = 0.003
N_WIRINGS_NEEDED = 6
MEAN_TOL = 0.02
RAW_LOSS_TOL = 0.005
CORNER_WIRINGS_MAX = 2
DIFFER_WIRINGS = 7
REALISM_THR_HZ = 60.0
READOUTS = ('raw', 'cal')
PARAM_KEYS = E83.PARAM_KEYS + ('wMBC1_scale',)
# the "differ" comparisons, (a, b)
COMPARISONS = (('R50', 'H50'), ('X50', 'PO50'), ('X50', 'TY50'), ('L100', 'B100'), ('D100', 'B100'),
               ('P100', 'B100'))
_CACHE: dict = {}


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E131_mbc1_cells/ (created)."""
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
    return f'S1_{arm}_{w}'


def n_resp(responses: pd.DataFrame) -> int:
    return int((responses.to_numpy(dtype=float) > 0.5).sum())


def kc_types_of(kc_index: pd.Index) -> np.ndarray:
    return np.asarray(kc_index.get_level_values(H.KC_TYPE_LEVEL), dtype=object)


def kc_ids_of(kc_index: pd.Index) -> pd.Index:
    return pd.Index(kc_index.get_level_values(H.KC_ID_LEVEL))


def model_kc_index(w: str = 'real') -> pd.Index:
    """The model's KC index (kc_id, kc_type) of a wiring's saved fit."""
    return E83.saved_fit(w)['saved']['on'][1].index


def n_steps() -> int:
    return int(round((TIME_END - TIME_PRE_START) / TIME_DT))


def pulse_steps() -> Tuple[int, int]:
    """[first, last + 1) of the steps with the odor on."""
    return (int(round((STIM_START - TIME_PRE_START) / TIME_DT)), int(round((STIM_END - TIME_PRE_START) / TIME_DT)))


def mbc1_scale(arm: str, s: float, late_ratio: Optional[float] = None) -> float:
    """Reading 1: MB-C1's 'wMBC1_scale' of an arm on a wiring with fitted KC-KC scale `s` (L100: times `late_ratio`)."""
    a = ARMS[arm]
    v = float(a['m']) * float(s)
    if a['form'] == 'late':
        if late_ratio is None or not (np.isfinite(late_ratio) and late_ratio > 0):
            raise ValueError(f'{arm} needs step 0\'s late ratio (got {late_ratio})')
        v *= float(late_ratio)
    return v


def arm_params(arm: str, p: Mapping[str, float], late_ratio: Optional[float] = None) -> Tuple[Dict[str, float], dict]:
    """(the arm's starting parameters: the fit's with wKCKC_scale kckc_on x s and wMBC1_scale; hold_rates' off_extra:
    wKCKC_scale kckc_off x s where it differs from APL on, else {})."""
    a = ARMS[arm]
    s = float(p['wKCKC_scale'])
    q = {k: float(p[k]) for k in E83.PARAM_KEYS}
    q['wKCKC_scale'] = float(a['kckc_on']) * s
    q['wMBC1_scale'] = mbc1_scale(arm, s, late_ratio)
    off = {} if a['kckc_off'] == a['kckc_on'] else {'wKCKC_scale': float(a['kckc_off']) * s}
    return q, off


def form_kws(form: str, w: str, step0: Optional[Mapping] = None) -> dict:
    """The MB-C1 model kwargs of a form, besides mbc1_inh: 'plain' none, 'late' mbc1_tau 1 s, 'strong' mbc1_thr (step 0's,
    per wiring), 'pooled' mbc1_pooled; 'off' is the option off (no MB-C1 kwarg at all)."""
    if form in ('plain', 'off'):
        return {}
    if form == 'late':
        return dict(mbc1_tau=LATE_TAU_S)
    if form == 'pooled':
        return dict(mbc1_pooled=True)
    if form == 'strong':
        if step0 is None:
            raise ValueError('the strong-drive form needs step 0')
        return dict(mbc1_thr=float(step0['threshold'][w]))
    raise ValueError(f'unknown form {form!r}')


def build_kws(arm: str, w: str, step0: Optional[Mapping] = None) -> dict:
    """The MB-C1 kwargs of an arm's simulator (`make_fixed_param_kc_simulator` model kwargs): mbc1_inh (the arm's scale,
    which 'wMBC1_scale' overrides), and mbc1_tau (L100), mbc1_thr (D100, step 0's per wiring) or mbc1_pooled (P100)."""
    a = ARMS[arm]
    lr = None
    if a['form'] == 'late':
        if step0 is None:
            raise ValueError(f'{arm} needs step 0')
        lr = float(step0['late_ratio'][w])
    s = float(E83.saved_fit(w)['p']['wKCKC_scale'])
    return dict(mbc1_inh=mbc1_scale(arm, s, lr), **form_kws(a['form'], w, step0))


def held_ok(n: int, target: int) -> bool:
    """Reading 8: the response count within RATE_TOL of the target."""
    return bool(abs(int(n) - int(target)) <= RATE_TOL * int(target))


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
# step 0 (no simulation)

def cell_weights(kc_index: pd.Index, pooled: bool = False, normalize: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """pyolfsysm's `connectome_MBC1_cell_weights` for the model's KCs, in their order (normalize=True: as fit_mb_model
    sets them; False: synapse counts): (out, KCs x cells; in, KCs x cells)."""
    from pyolfsysm.mb_model import connectome_MBC1_cell_weights
    ids = kc_ids_of(kc_index)
    w_out, w_in = connectome_MBC1_cell_weights(ids, pooled=pooled, normalize=normalize)
    if not (np.array_equal(w_out.index.to_numpy(), ids.to_numpy()) and np.array_equal(w_in.index.to_numpy(), ids.to_numpy())):
        raise RuntimeError("MB-C1's cell weights are not in the model KCs' order")
    return w_out.astype(float), w_in.astype(float)


def dose_check(kc_index: pd.Index) -> dict:
    """Reading 1's check on the model's KCs, for the two cells and for the pooled cell: the equivalent matrix's total
    and the hemibrain matrix's (both N), each cell's input weights' sum (1), the equivalent matrix against
    `connectome_MBC1_KCKC_weights` (same pooling), and the total output for uniform activity against the matrix's."""
    from pyolfsysm.mb_model import connectome_KCKC_weights, connectome_MBC1_KCKC_weights
    ids = kc_ids_of(kc_index)
    n = len(ids)
    HB = connectome_KCKC_weights(ids).to_numpy(dtype=float)
    row = dict(n_kcs=n, hb_total=float(HB.sum()), hb_uniform_output=float((HB @ np.ones(n)).sum()))
    ok = bool(np.isclose(HB.sum(), n, rtol=1e-9))
    for name, pooled in (('cells', False), ('pooled', True)):
        o, i = cell_weights(kc_index, pooled=pooled)
        Wo, Wi = o.to_numpy(), i.to_numpy()
        E = Wo @ Wi.T
        M = connectome_MBC1_KCKC_weights(ids, pooled=pooled).to_numpy(dtype=float)
        uni = float((Wo @ (Wi.T @ np.ones(n))).sum())
        r = dict(n_cells=int(Wo.shape[1]), equiv_total=float(E.sum()), in_sums=[float(x) for x in Wi.sum(axis=0)],
                 out_totals=[float(x) for x in Wo.sum(axis=0)], equiv_vs_kckc_weights_max_diff=float(np.abs(E - M).max()),
                 uniform_output=uni, uniform_output_over_matrix=uni / row['hb_uniform_output'])
        r['ok'] = bool(np.isclose(E.sum(), n, rtol=1e-9) and np.allclose(Wi.sum(axis=0), 1.0, rtol=1e-12)
                       and np.allclose(E, M, rtol=1e-9, atol=1e-12) and np.isclose(uni, row['hb_uniform_output'], rtol=1e-9))
        ok &= r['ok']
        row[name] = r
    row['ok'] = bool(ok)
    return row


def synapse_table(kc_index: pd.Index, claw_index: pd.MultiIndex) -> pd.DataFrame:
    """Each MB-C1 cell's synapses in and out on the model's KCs by KC type ('all' first), with APL's synapses onto them
    (on the wiring's claws, as the model's APL->KC weights count them, and all in the calyx) and the KC->KC synapses
    among them beside: totals ('<x>') and per KC ('<x>_per_kc')."""
    from pyolfsysm.mb_model import connectome_APL_KC_synapse_counts, connectome_KCKC_weights
    ids = kc_ids_of(kc_index)
    types = kc_types_of(kc_index)
    o, i = cell_weights(kc_index, normalize=False)
    cells = [str(c) for c in o.columns]
    KK = connectome_KCKC_weights(ids, normalize=False).to_numpy(dtype=float)       # post x pre
    apl, _ = connectome_APL_KC_synapse_counts(claw_index, 'APL>KC')
    apl = apl.reindex(ids)
    if apl.isna().any().any():
        raise RuntimeError("APL's synapse counts lack some model KCs")
    cols = {'apl_out_claws': apl['kept'].to_numpy(dtype=float),
            'apl_out_calyx': (apl['kept'] + apl['dropped_claw'] + apl['no_claw']).to_numpy(dtype=float),
            'kckc_in': KK.sum(axis=1), 'kckc_out': KK.sum(axis=0)}
    for j, c in enumerate(cells):
        cols[f'mbc1_out_{c}'] = o.iloc[:, j].to_numpy()
        cols[f'mbc1_in_{c}'] = i.iloc[:, j].to_numpy()
    cols['mbc1_out'] = o.to_numpy().sum(axis=1)
    cols['mbc1_in'] = i.to_numpy().sum(axis=1)
    rows = []
    for t in ('all',) + TYPES:
        m = np.ones(len(ids), dtype=bool) if t == 'all' else types == t
        r = dict(kc_type=t, n_kcs=int(m.sum()))
        for k, v in cols.items():
            r[k] = float(v[m].sum())
            r[f'{k}_per_kc'] = float(v[m].mean()) if m.any() else np.nan
        rows.append(r)
    df = pd.DataFrame(rows)
    tot = df.set_index('kc_type').loc['all']
    for k in ('mbc1_out', 'mbc1_in'):
        df[f'{k}_share'] = df[k] / float(tot[k])
    df.attrs['cells'] = cells
    return df


def dose_arithmetic(syn: pd.DataFrame) -> dict:
    """The dose at which an MB-C1 output synapse weighs what a KC->KC synapse does, for uniform KC activity: on average
    over MB-C1's output synapses (their total over the KC->KC synapses), and per cell c (the equivalent matrix's total
    over cell c's input total times the KC->KC synapses; reading 1's normalization); and MB-C1's share of the inhibitory
    synapses onto γ against APL's (on the claws and in the whole calyx)."""
    t = syn.set_index('kc_type')
    cells = [k[len('mbc1_out_'):] for k in t.columns if k.startswith('mbc1_out_') and not k.endswith(('_per_kc', '_share'))]
    a = t.loc['all']
    kk = float(a['kckc_in'])
    pairs = sum(float(a[f'mbc1_out_{c}']) * float(a[f'mbc1_in_{c}']) for c in cells)
    g = t.loc['g']
    return dict(kckc_synapses=kk, mbc1_out_synapses=float(a['mbc1_out']), mbc1_in_synapses=float(a['mbc1_in']),
                m_equal_synapse=float(a['mbc1_out']) / kk,
                m_equal_synapse_per_cell={c: pairs / (float(a[f'mbc1_in_{c}']) * kk) for c in cells},
                mbc1_out_per_g_kc=float(g['mbc1_out_per_kc']), apl_out_claws_per_g_kc=float(g['apl_out_claws_per_kc']),
                apl_out_calyx_per_g_kc=float(g['apl_out_calyx_per_kc']),
                g_share_mbc1_vs_apl_claws=float(g['mbc1_out'] / (g['mbc1_out'] + g['apl_out_claws'])),
                g_share_mbc1_vs_apl_calyx=float(g['mbc1_out'] / (g['mbc1_out'] + g['apl_out_calyx'])),
                mbc1_out_share_g=float(g['mbc1_out_share']), mbc1_in_share_g=float(g['mbc1_in_share']),
                mbc1_over_apl_claws_out=float(a['mbc1_out'] / a['apl_out_claws']))


def spike_times(w: str, cond: str, counts: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """E121 step 0's spike times of a wiring and condition (outputs/E121_onset_race/step0/: kc and odor positions in the
    outputs' order, step, time_s), checked to add up to `counts` (default the family's saved spike counts) and to be on
    the time grid."""
    path = H.outputs_dir('E121_onset_race') / 'step0' / f'{E121_LABEL.get(w, w)}_{cond}_spike_times.parquet'
    if not path.exists():
        raise RuntimeError(f'{path} is missing (E121 step 0 on viki)')
    spk = pd.read_parquet(path)
    C = E83.saved_fit(w)['saved'][cond][1] if counts is None else counts
    check_spikes(spk, C)
    return spk


def check_spikes(spk: pd.DataFrame, counts: pd.DataFrame) -> None:
    """Spikes (kc, odor, step, time_s) add up to `counts` (KCs x odors) and lie on the time grid (RuntimeError)."""
    n = np.zeros(counts.shape)
    np.add.at(n, (spk['kc'].to_numpy(), spk['odor'].to_numpy()), 1.0)
    if not np.array_equal(n, counts.to_numpy(dtype=float)):
        raise RuntimeError(f'the spikes do not add up to the spike counts ({int((n != counts.to_numpy()).sum())} differ)')
    t = TIME_PRE_START + spk['step'].to_numpy(dtype=float) * TIME_DT
    if len(spk) and not np.allclose(t, spk['time_s'].to_numpy(dtype=float), atol=1e-9):
        raise RuntimeError('the spikes are not on the time grid')
    if len(spk) and not (0 <= spk['step'].min() and spk['step'].max() < n_steps()):
        raise RuntimeError('spikes outside the simulated steps')


def activity_trace(spikes: pd.DataFrame, w_in: np.ndarray, n_odors: int, tau: float, dt: float = TIME_DT,
                   steps: Optional[int] = None) -> np.ndarray:
    """Reading 4: each cell's activity (cells x odors x steps, Hz) by olfsysm's recursion m(t) = (1 - dt / tau) m(t - 1)
    + w_in @ spikes(t - 1) / tau, m(0) = 0; `w_in` cells x KCs (rv.kc.wKCMBC1's layout)."""
    w_in = np.atleast_2d(np.asarray(w_in, dtype=float))
    steps = n_steps() if steps is None else int(steps)
    C = w_in.shape[0]
    X = np.zeros((C, n_odors, steps))
    kc, od, st = (spikes[k].to_numpy(dtype=int) for k in ('kc', 'odor', 'step'))
    for c in range(C):
        np.add.at(X[c], (od, st), w_in[c, kc])
    decay = 1.0 - dt / tau
    tr = np.zeros((C, n_odors, steps))
    m = np.zeros((C, n_odors))
    for t in range(1, steps):
        m = decay * m + X[:, :, t - 1] / tau
        tr[:, :, t] = m
    return tr


def cell_activity(spikes: pd.DataFrame, w_in: np.ndarray, n_odors: int, tau: float, thr: float = 0.0,
                  dt: float = TIME_DT, steps: Optional[int] = None,
                  pulse: Optional[Tuple[int, int]] = None) -> Dict[str, np.ndarray]:
    """Reading 4: each cell's activity (cells x odors): 'pulse_mean' (over the steps with the odor on) and 'peak' (over
    the trial), in Hz, and 'pulse_output' (the integral over the pulse of max(activity - thr, 0), in Hz s)."""
    tr = activity_trace(spikes, w_in, n_odors, tau, dt, steps)
    p0, p1 = pulse_steps() if pulse is None else pulse
    seg = tr[:, :, p0:p1]
    return dict(pulse_mean=seg.mean(axis=2), peak=tr.max(axis=2), pulse_output=np.maximum(seg - thr, 0.0).sum(axis=2) * dt)


def late_ratio(act50: Mapping[str, np.ndarray], act_late: Mapping[str, np.ndarray], out_w: np.ndarray) -> float:
    """Reading 6: the cells' summed output over the pulse at 50 ms over the same at 1 s (each cell's pulse means summed
    over the odors, weighted by its total output weight `out_w`)."""
    out_w = np.asarray(out_w, dtype=float).ravel()
    num = float(out_w @ np.asarray(act50['pulse_mean']).sum(axis=1))
    den = float(out_w @ np.asarray(act_late['pulse_mean']).sum(axis=1))
    if not den > 0:
        raise ValueError(f'the late cells have no output over the pulse ({den})')
    return num / den


def cells_agreement(pulse_mean: np.ndarray) -> float:
    """Reading 7: the Pearson r over the odors of the two cells' pulse means (nan for one cell)."""
    a = np.atleast_2d(np.asarray(pulse_mean, dtype=float))
    if a.shape[0] < 2:
        return float('nan')
    if a.shape[0] != 2:
        raise ValueError(f'{a.shape[0]} cells')
    return float(np.corrcoef(a[0], a[1])[0, 1])


def pooled_gate(agreement: Mapping[str, float]) -> dict:
    """Reading 7: open if any agreement ({'<wiring>|<cond>': r}) is below POOLED_GATE_R (or not finite)."""
    low = {k: float(v) for k, v in agreement.items() if not (np.isfinite(v) and v >= POOLED_GATE_R)}
    return dict(open=bool(low), below=low, min_r=float(np.nanmin(list(agreement.values()))) if agreement else np.nan,
                rule=f'P100 runs if the two cells agree at r < {POOLED_GATE_R} on any wiring (APL on or off)')


def gamma_input(counts: pd.DataFrame, s: float, equiv: np.ndarray, hb: np.ndarray) -> pd.DataFrame:
    """Each γ KC's MB-C1 input at m 1 beside its KC-KC input: trial integrals summed over the odors, s x equiv @ counts
    and s x hb @ counts (KCs x KCs, post x pre), with their ratio."""
    C = counts.to_numpy(dtype=float)
    types = kc_types_of(counts.index)
    g = types == 'g'
    mb = float(s) * (np.asarray(equiv)[g] @ C).sum(axis=1)
    kk = float(s) * (np.asarray(hb)[g] @ C).sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(kk > 0, mb / kk, np.nan)
    return pd.DataFrame(dict(kc_id=kc_ids_of(counts.index)[g], mbc1_input=mb, kckc_input=kk, mbc1_over_kckc=ratio))


def part_step0() -> dict:
    """Step 0 (the module docstring's step0 part): step0_synapses.csv, step0_activity.csv, step0_gamma_input.csv,
    step0.json, provenance_step0.json."""
    from pyolfsysm.mb_model import connectome_KCKC_weights
    out = out_dir()
    idx = model_kc_index('real')
    ids = kc_ids_of(idx)
    for w in WIRINGS:
        if not kc_ids_of(model_kc_index(w)).equals(ids):
            raise RuntimeError(f"{w}'s KCs are not the hemibrain wiring's, in its order")
    dc = dose_check(idx)
    if not dc['ok']:
        raise RuntimeError(f'reading 1\'s dose check failed: {json.dumps(dc, default=str)}')
    claw_index = E83.wiring_parts('real')['W'].index
    syn = synapse_table(idx, claw_index)
    syn.to_csv(out / 'step0_synapses.csv', index=False)
    arith = dose_arithmetic(syn)
    w_out, w_in = cell_weights(idx)
    cells = [str(c) for c in w_out.columns]
    out_w = w_out.to_numpy().sum(axis=0)
    Wi = w_in.to_numpy().T                                           # cells x KCs, as rv.kc.wKCMBC1
    equiv = w_out.to_numpy() @ w_in.to_numpy().T
    hb = connectome_KCKC_weights(ids).to_numpy(dtype=float)
    act_rows, g_rows = [], []
    thr, lr, agree, summ = {}, {}, {}, {}
    for w in WIRINGS:
        F = E83.saved_fit(w)
        s = float(F['p']['wKCKC_scale'])
        n_od = F['saved']['on'][1].shape[1]
        odors = [str(c) for c in F['saved']['on'][1].columns]
        acts = {}
        for c in CONDS:
            spk = spike_times(w, c)
            for tau in (PLAIN_TAU_S, LATE_TAU_S):
                a = cell_activity(spk, Wi, n_od, tau)
                acts[(c, tau)] = a
                for j, cell in enumerate(cells):
                    for o in range(n_od):
                        act_rows.append(dict(wiring=w, condition=c, tau=tau, cell=cell, odor=odors[o],
                                             pulse_mean=float(a['pulse_mean'][j, o]), peak=float(a['peak'][j, o])))
            agree[f'{w}|{c}'] = cells_agreement(acts[(c, PLAIN_TAU_S)]['pulse_mean'])
            gi = gamma_input(F['saved'][c][1], s, equiv, hb)
            gi.insert(0, 'condition', c)
            gi.insert(0, 'wiring', w)
            g_rows.append(gi)
            summ[f'{w}|{c}'] = dict(mbc1_over_kckc_total=float(gi['mbc1_input'].sum() / gi['kckc_input'].sum()),
                                    mbc1_over_kckc_median=float(np.nanmedian(gi['mbc1_over_kckc'])),
                                    pulse_mean_hz={cell: float(acts[(c, PLAIN_TAU_S)]['pulse_mean'][j].mean())
                                                   for j, cell in enumerate(cells)},
                                    peak_hz={cell: float(acts[(c, PLAIN_TAU_S)]['peak'][j].max())
                                             for j, cell in enumerate(cells)})
        thr[w] = float(acts[('on', PLAIN_TAU_S)]['peak'].max())
        lr[w] = late_ratio(acts[('off', PLAIN_TAU_S)], acts[('off', LATE_TAU_S)], out_w)
    pd.DataFrame(act_rows).to_csv(out / 'step0_activity.csv', index=False)
    pd.concat(g_rows, ignore_index=True).to_csv(out / 'step0_gamma_input.csv', index=False)
    pg = pooled_gate(agree)
    v = dict(generated=_now(), host=_host(), threshold=thr, late_ratio=lr, agreement=agree, pooled_gate=pg,
             dose_check=dc, dose_arithmetic=arith, cells=cells, gamma_input=summ,
             time_grid=dict(pre_start=TIME_PRE_START, dt=TIME_DT, end=TIME_END, stim=[STIM_START, STIM_END],
                            pulse_steps=list(pulse_steps())),
             readings=dict(threshold='reading 5: the largest activity of either cell over the odors and the trial, APL '
                                     'on, the family\'s run, tau 50 ms',
                           late_ratio='reading 6: the cells\' summed pulse output at 50 ms over that at 1 s, APL off',
                           agreement='reading 7: Pearson r over the 16 odors of the two cells\' pulse means'),
             kc_ids_sha256=sha256_array(ids.to_numpy(dtype=float)))
    _write_json(out / 'step0.json', v)
    _write_json(out / 'provenance_step0.json', dict(
        generated=_now(), host=_host(), script=SCRIPT, part='step0', commits=H.code_commits(),
        inputs=dict(family={w: list(E83.members()[w]) for w in WIRINGS},
                    spike_times='outputs/E121_onset_race/step0/<w>_<cond>_spike_times.parquet (viki)',
                    mbc1=mbc1_provenance())))
    print(json.dumps(dict(threshold=thr, late_ratio=lr, agreement=agree, pooled_gate=pg['open'],
                          dose=arith['m_equal_synapse'], g_share=arith['g_share_mbc1_vs_apl_claws']), default=str),
          flush=True)
    return v


def sha256_array(a: np.ndarray) -> str:
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=float), dtype='<f8').tobytes()).hexdigest()


def step0_fixed() -> dict:
    """step0.json (raises if step 0 has not run)."""
    v = _read_json(out_dir() / 'step0.json')
    if not v:
        raise RuntimeError('step 0 has not run (--part step0)')
    return v


def mbc1_provenance() -> dict:
    """The MB-C1 tables pyolfsysm reads (`mbc1_synapse_paths`, file names and sha256) and the function used."""
    import hashlib
    from pyolfsysm.mb_model import mbc1_synapse_paths
    return dict(source='pyolfsysm.mb_model.connectome_MBC1_cell_weights(kc_ids, pooled=..., normalize=True)',
                tables={k: dict(file=Path(mbc1_synapse_paths[k]).name,
                                sha256=hashlib.sha256(Path(mbc1_synapse_paths[k]).read_bytes()).hexdigest())
                        for k in ('KC>MBC1', 'MBC1>KC')})


# ---------------------------------------------------------------------------
# the runner and the checks (simulations)

def simulator(w: str, form: str, params: Mapping[str, float], step0: Optional[Mapping] = None) -> Tuple[Callable, dict]:
    """(simulate, info): E83's simulator of wiring `w` at `params`, with MB-C1 of a form ('off': the option off; else
    mbc1_inh params['wMBC1_scale'] and `form_kws`)."""
    if form == 'off':
        return E83.build_simulator(w, params={k: float(params[k]) for k in E83.PARAM_KEYS})
    extra = dict(mbc1_inh=float(params['wMBC1_scale']), **form_kws(form, w, step0))
    return E83.build_simulator(w, params=dict(params), extra_kws=extra)


def build_runner(arm: str, w: str, params: Mapping[str, float], step0: Optional[Mapping] = None) -> Tuple[Callable, dict]:
    """(simulate, info) for an arm on a wiring, cached for the arms of the same form on the same wiring (one simulator
    at a time): its scales are parameters, so the cached simulator serves any of them."""
    form = ARMS[arm]['form']
    key = (w, form, json.dumps(form_kws(form, w, step0), sort_keys=True))
    c = _CACHE.get('runner')
    if c is not None and c[0] == key:
        return c[1], c[2]
    _CACHE.pop('runner', None)
    gc.collect()
    sim, info = simulator(w, form, params, step0)
    _CACHE['runner'] = (key, sim, info)
    return sim, info


def _ndiff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> Dict[str, int]:
    return dict(counts=int((a[1].to_numpy(dtype=float) != b[1].to_numpy(dtype=float)).sum()),
                responses=int((a[0].to_numpy(dtype=float) != b[0].to_numpy(dtype=float)).sum()))


def _off(p: Mapping[str, float]) -> Dict[str, float]:
    return dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)


def _check_grid(info: dict) -> dict:
    mp = info['cells']['mp'].cell_contents
    got = dict(pre_start=float(mp.time_pre_start), dt=float(mp.time_dt), end=float(mp.time_end),
               stim_start=float(mp.time_stim_start), stim_end=float(mp.time_stim_end))
    want = dict(pre_start=TIME_PRE_START, dt=TIME_DT, end=TIME_END, stim_start=STIM_START, stim_end=STIM_END)
    if not all(np.isclose(got[k], want[k]) for k in want):
        raise RuntimeError(f'the time grid is {got}, not {want}')
    return got


def weights_match(info: dict, kc_index: pd.Index, pooled: bool = False) -> dict:
    """The simulator's MB-C1 weights (rv.kc.wKCMBC1, the unscaled output weights) against `cell_weights` exactly."""
    cells = info['cells']
    rv = cells['rv'].cell_contents
    w_out, w_in = cell_weights(kc_index, pooled=pooled)
    a_in = np.asarray(rv.kc.wKCMBC1, dtype=float)
    a_out = np.asarray(cells['wMBC1KC_unscaled'].cell_contents, dtype=float)
    return dict(w_in=bool(np.array_equal(a_in, w_in.to_numpy().T)), w_out=bool(np.array_equal(a_out, w_out.to_numpy())))


def part_check(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None) -> dict:
    """The checks (the module docstring): check_<host>_<lane>.json."""
    from experiments import E108_paa_extrapolation as E108
    _preload()
    lane = lane or lane_name(wirings)
    path = out_dir() / f'check_{_host()}_{lane}.json'
    res = dict(_read_json(path), host=_host(), lane=lane, script=SCRIPT, commits=E108._commits_clean())
    res.setdefault('wirings', {})
    st = _read_json(out_dir() / 'step0.json') or None
    for w in wirings:
        r = res['wirings'].get(w, {})
        t0 = time.time()
        P = E83.wiring_parts(w)
        p, saved = dict(P['p']), P['saved']
        s = float(p['wKCKC_scale'])
        if not r.get('passes'):
            sim, info = simulator(w, 'off', p)
            grid = _check_grid(info)
            off_ = {c: _ndiff(sim(p if c == 'on' else _off(p)), saved[c]) for c in CONDS}
            del sim, info
            gc.collect()
            q = dict(p, wMBC1_scale=0.0)
            sim, info = simulator(w, 'plain', q)
            zero = {c: _ndiff(sim(q if c == 'on' else _off(q)), saved[c]) for c in CONDS}
            wm = weights_match(info, saved['on'][1].index)
            r = dict(option_off=off_, option_scale0=zero, weights_match=wm, time_grid=grid)
            r['passes'] = bool(all(v == 0 for d in (off_, zero) for x in d.values() for v in x.values()) and all(wm.values()))
            if w == 'real':
                live = sim(_off(dict(p, wMBC1_scale=s)))
                r['option_m1_off'] = dict(n_differing=_ndiff(live, saved['off'])['counts'], n_responses=n_resp(live[0]),
                                          target_off=int(P['target_off']))
                r['passes'] = bool(r['passes'] and r['option_m1_off']['n_differing'] > 0)
            del sim, info
            gc.collect()
        if st and not r.get('threshold_checked'):
            q = dict(p, wMBC1_scale=s)
            sim, info = simulator(w, 'strong', q, st)
            on = _ndiff(sim(q), saved['on'])
            off_ = _ndiff(sim(_off(q)), saved['off'])
            r.update(threshold=float(st['threshold'][w]), threshold_on=on, threshold_off_n_differing=off_['counts'],
                     threshold_passes=bool(on['counts'] == 0 and on['responses'] == 0), threshold_checked=_now())
            del sim, info
            gc.collect()
        r.update(seconds=time.time() - t0, checked=_now())
        res['wirings'][w] = r
        _write_json(path, res)
        print(f'E131 check {w}: {json.dumps(r, default=str)}', flush=True)
    res['passes'] = all(res['wirings'].get(w, {}).get('passes') for w in wirings)
    _write_json(path, res)
    return res


def check_wirings(key: str = 'passes') -> set:
    """The wirings that passed a check (`key` 'passes' for (i) to (iii), 'threshold_passes' for (iv)) in any host's
    check_*.json."""
    ok = set()
    for p in out_dir().glob('check_*.json'):
        for w, r in _read_json(p).get('wirings', {}).items():
            if r.get(key):
                ok.add(w)
    return ok


def check_passed(need_threshold: bool = False) -> bool:
    ok = check_wirings() >= set(WIRINGS)
    return bool(ok and (not need_threshold or check_wirings('threshold_passes') >= set(WIRINGS)))


def held(arm: str, w: str, label: str, step0: Optional[Mapping] = None) -> Tuple[dict, Callable, dict]:
    """One held point of an arm on a wiring: E77_input_screen.hold_rates on `build_runner`'s simulator, from the arm's
    parameters (`arm_params`; X50's APL-off KC-KC scale through off_extra). Returns (hold_rates' result with 'held_on',
    'held_off', 'corner', 'off_extra' added; the simulator; its info)."""
    from experiments.E77_input_screen import hold_rates
    P = E83.wiring_parts(w)
    lr = float(step0['late_ratio'][w]) if ARMS[arm]['form'] == 'late' else None
    q, off_extra = arm_params(arm, P['p'], lr)
    box = {}

    def build():
        box['sim'], box['info'] = build_runner(arm, w, q, step0)
        return box['sim'], box['info']

    res = hold_rates(build, target_off=P['target_off'], target_on=P['target_on'], params=dict(q), label=label,
                     off_extra=off_extra or None)
    res['held_on'] = held_ok(res['n_on'], res['target_on'])
    res['held_off'] = held_ok(res['n_off'], res['target_off'])
    res['corner'] = bool(res['at_bound_thr'] or res['at_bound_apl'] or not (res['held_on'] and res['held_off']))
    res['off_extra'] = dict(off_extra)
    return res, box['sim'], box['info']


def kc_apl_weights(info: dict, kc_index: pd.Index) -> np.ndarray:
    """Each KC's unscaled APL->claw weights summed over its claws, in the model KCs' order."""
    s = info['wAPLKC']
    w = pd.Series(s.to_numpy(dtype=float), index=s.index.get_level_values(H.KC_ID_LEVEL)).groupby(level=0).sum()
    return w.reindex(kc_ids_of(kc_index)).fillna(0.0).to_numpy(dtype=float)


def pulse_integral(spikes: pd.DataFrame, n_kcs: int, n_odors: int, tau: float, dt: float = TIME_DT,
                   pulse: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """Each KC's filtered activity (tau) integrated over the pulse (KCs x odors): a spike at step s adds
    decay^(a - s - 1) (1 - decay^(p1 - a)), a = max(s + 1, p0), where a < p1 (the recursion of `activity_trace`)."""
    p0, p1 = pulse_steps() if pulse is None else pulse
    decay = 1.0 - dt / tau
    s = spikes['step'].to_numpy(dtype=int)
    a = np.maximum(s + 1, p0)
    v = np.where(a < p1, decay ** (a - s - 1) * (1.0 - decay ** np.maximum(p1 - a, 0)), 0.0)
    out = np.zeros((n_kcs, n_odors))
    np.add.at(out, (spikes['kc'].to_numpy(dtype=int), spikes['odor'].to_numpy(dtype=int)), v)
    return out


def inhibition_shares(mbc1: np.ndarray, apl: np.ndarray, kckc: np.ndarray, types: np.ndarray) -> dict:
    """Reading 12: the inhibition delivered over the pulse (KCs x odors each), per KC type and over all KCs: MB-C1's
    share against MB-C1's and APL's ('share_mbc1_vs_apl_<t>') and each source's share of all three."""
    row = {}
    for t in ('all',) + TYPES:
        m = np.ones(len(types), dtype=bool) if t == 'all' else np.asarray(types, dtype=object) == t
        k = 'all' if t == 'all' else TYPE_COL[t]
        mb, ap, kk = (float(np.asarray(x)[m].sum()) for x in (mbc1, apl, kckc))
        tot = mb + ap + kk
        row[f'share_mbc1_vs_apl_{k}'] = mb / (mb + ap) if mb + ap > 0 else np.nan
        for name, x in (('mbc1', mb), ('apl', ap), ('kckc', kk)):
            row[f'share3_{name}_{k}'] = x / tot if tot > 0 else np.nan
    return row


def cell_realism(sim: Callable, info: dict, res: Mapping, arm: str) -> dict:
    """Reading 12: the held point simulated again with its spike times (and APL's trace with APL on), which must
    reproduce the held counts; each cell's activity with APL on and off, and the inhibition shares with APL on."""
    from pyolfsysm.mb_model import get_sim_var
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    kc_index = res['on'][1].index
    types = kc_types_of(kc_index)
    n_kc, n_od = res['on'][1].shape
    w_in = np.asarray(rv.kc.wKCMBC1, dtype=float)
    w_out = np.asarray(cells['wMBC1KC_unscaled'].cell_contents, dtype=float)
    tau, thr = float(mp.kc.mbc1_tau), float(mp.kc.mbc1_thr)
    p0, p1 = pulse_steps()
    row = dict(mbc1_tau=tau, mbc1_thr=thr, mbc1_pooled=bool(w_in.shape[0] == 1), mbc1_n_cells=int(w_in.shape[0]))
    for c, q in (('on', res['params']), ('off', res['params_off'])):
        mp.kc.save_inh_sims = c == 'on'
        try:
            r, s, spk = sim(dict(q), spike_times=True)
            inh = np.asarray(get_sim_var(rv.kc, 'inh_sims'), dtype=float) if c == 'on' else None
        finally:
            mp.kc.save_inh_sims = False
        d = _ndiff((r, s), res[c])
        if d['counts'] or d['responses']:
            raise RuntimeError(f'{arm}: the held point simulated again differs from the held run ({c}: {d})')
        a = cell_activity(spk, w_in, n_od, tau, thr)
        for j in range(w_in.shape[0]):
            row[f'mbc1_act_{c}_c{j}'] = float(a['pulse_mean'][j].mean())
            row[f'mbc1_peak_{c}_c{j}'] = float(a['peak'][j].max())
        row[f'mbc1_act_{c}'] = float(a['pulse_mean'].mean())
        row[f'mbc1_peak_{c}'] = float(a['peak'].max())
        if c == 'on':
            mb = float(q['wMBC1_scale']) * (w_out @ a['pulse_output'])                       # KCs x odors
            inh = inh.reshape(inh.shape[0], -1)[order]                                       # odors x steps
            apl_int = inh[:, p0 - 1:p1 - 1].sum(axis=1) * TIME_DT                           # APL(t - 1) over the pulse
            ap = np.outer(float(q['wAPLKC_scale']) * kc_apl_weights(info, kc_index), apl_int)
            Wkk = cells['wKCKC_unscaled'].cell_contents
            I = pulse_integral(spk, n_kc, n_od, float(mp.kc.kc_kc_tau))
            kk = float(q['wKCKC_scale']) * np.asarray(Wkk @ I)
            row.update(inhibition_shares(mb, ap, kk, types))
    return row


def planned(wirings: Sequence[str] = WIRINGS, arms: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """The (arm, wiring) points of a lane: group by group (B, R, X, the forms), each wiring's arms of the group
    together (one simulator per form and wiring)."""
    arms = tuple(ARMS) if arms is None else tuple(arms)
    ws = [w for w in WIRINGS if w in wirings]
    out = []
    for g in GROUPS:
        ga = [a for a in ARMS if ARMS[a]['group'] == g and a in arms]
        if g == 'forms':
            out += [(a, w) for a in ga for w in ws]
        else:
            out += [(a, w) for w in ws for a in ga]
    return out


def runnable_arms(step0: Optional[Mapping] = None, threshold_wirings: Optional[set] = None) -> Tuple[str, ...]:
    """The arms that may run: without step 0 none of L100, D100, P100; P100 only if step 0's pooled gate is open; D100
    only once check (iv) has passed on every wiring."""
    out = []
    for a in ARMS:
        if a in STEP0_ARMS and not step0:
            continue
        if a == 'P100' and not step0.get('pooled_gate', {}).get('open'):
            continue
        if a == 'D100' and not (threshold_wirings or set()) >= set(WIRINGS):
            continue
        out.append(a)
    return tuple(out)


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


def save_run(tag: str, o: Mapping[str, tuple]) -> None:
    """The run's conditions, locally and in the model-output archive (which another host's score part reads)."""
    for d in (local_dir('step1'), archive_dir() / 'step1'):
        d.mkdir(parents=True, exist_ok=True)
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            for k, x in (('responses', o[c][0]), ('spike_counts', o[c][1])):
                p = d / f'{tag}_{cc}_{k}.parquet'
                tmp = p.with_name(f'.{p.name}.{os.getpid()}.tmp')
                x.to_parquet(tmp)
                os.replace(tmp, p)


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    for d in (local_dir('step1'), archive_dir() / 'step1'):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_row(tag: str, arm: str, w: str, res: Mapping, o: Mapping[str, tuple], real: Mapping, cc: dict, lane: str,
            step0: Optional[Mapping] = None) -> dict:
    from experiments import E123_kckc_form as E123
    from experiments import E129_block_survivor_static as E129
    P = E83.wiring_parts(w)
    a = ARMS[arm]
    s = float(P['p']['wKCKC_scale'])
    po = res['params_off']
    row = dict(tag=tag, arm=arm, group=a['group'], form=a['form'], m=a['m'], wiring=w, kckc_on_factor=a['kckc_on'],
               kckc_off_factor=a['kckc_off'], kckc_scale=float(res['params']['wKCKC_scale']),
               kckc_scale_off=float(po['wKCKC_scale']), fit_kckc_scale=s, fit_folder=P['folder'], fit_tag=P['tag'],
               thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], n_on=res['n_on'], n_off=res['n_off'],
               target_on=res['target_on'], target_off=res['target_off'], at_bound_thr=bool(res['at_bound_thr']),
               at_bound_apl=bool(res['at_bound_apl']), held_on=res['held_on'], held_off=res['held_off'],
               corner=res['corner'], n_sims=res['n_sims'], seconds=res.get('seconds'),
               **{k: float(res['params'][k]) for k in PARAM_KEYS},
               params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
               params_off=json.dumps({k: float(v) for k, v in po.items()}), host=_host(), lane=lane)
    row['mbc1_inh'] = row['wMBC1_scale']
    row['mbc1_claws_only'] = False
    arith = (step0 or {}).get('dose_arithmetic', {})
    if arith:
        # MB-C1's weight per output synapse against a KC->KC synapse's at the fit's KC-KC scale (reading 12)
        row['mbc1_per_synapse_vs_kckc'] = row['wMBC1_scale'] / s / float(arith['m_equal_synapse'])
    row.update(real)
    row.update(E83.realism(o))
    row.update(E129.type_rates(o))
    row.update(E123.bands(row))
    row['thr_below_60'] = bool(float(row['fixed_thr']) < REALISM_THR_HZ)
    row.update(cc)
    return row


def part_runs(wirings: Sequence[str] = WIRINGS, lane: Optional[str] = None, reclaim: bool = False,
              arms: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The fixed-scale arms (the module docstring), resumable, one point claimed at a time."""
    from experiments import E108_paa_extrapolation as E108
    if not check_passed():
        raise SystemExit('E131 runs refused: checks (i) to (iii) have not passed on every wiring (run --part check)')
    _preload()
    st = _read_json(out_dir() / 'step0.json') or None
    ok = runnable_arms(st, check_wirings('threshold_passes'))
    arms = tuple(a for a in (ARMS if arms is None else arms) if a in ARMS)
    skipped = [a for a in arms if a not in ok]
    if skipped:
        print(f'E131: {skipped} not run now (step 0, its pooled gate or check (iv) not there yet)', flush=True)
    arms = tuple(a for a in arms if a in ok)
    lane = lane or lane_name(wirings)
    table = out_dir() / f'runs_{_host()}_{lane}.csv'
    prov_path = out_dir() / f'provenance_{_host()}_{lane}.json'
    cc = E108._commits_clean()
    for arm, w in planned(wirings, arms):
        tag = tag_of(arm, w)
        if tag in set(runs_table()['tag']):
            print(f'E131 {tag}: done', flush=True)
            continue
        if not claim(tag, lane, reclaim):
            print(f'E131 {tag}: claimed by another process, skipped', flush=True)
            continue
        try:
            E108._still_clean(cc)
            P = E83.wiring_parts(w)
            res, sim, info = held(arm, w, f'E131 {tag}', st)
            o = {'on': res['on'], 'off': res['off']}
            real = cell_realism(sim, info, res, arm)
            save_run(tag, o)
            row = run_row(tag, arm, w, res, o, real, cc, lane, st)
            E83._append_row(table, row)
            prov = _read_json(prov_path) or dict(part='runs', lane=lane, script=SCRIPT, points={})
            prov['points'][tag] = dict(H.model_provenance(orn_deltas=P['orn']), held=dict(
                change=(f"arm {arm} ({ARM_LABELS[arm]}); both response rates re-held (E77_input_screen.hold_rates: "
                        "fixed_thr bisected with APL off, then the tied APL scales with APL on)"),
                fit=f"{P['folder']}:{P['tag']}", kws=res.get('kws'), off_extra=res['off_extra'],
                mbc1=dict(build_kws(arm, w, st), **mbc1_provenance()), history=res['history']))
            _write_json(prov_path, prov)
            print(f"E131 {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/"
                  f"{res['target_on']} on, {res['n_off']}/{res['target_off']} off); corner {res['corner']}; thr "
                  f"{row['fixed_thr']:.1f}; bands {row['in_bands']}; γ share {row.get('share_mbc1_vs_apl_g', np.nan):.3f}; "
                  f"{res.get('seconds', 0):.0f}s", flush=True)
            del res, o
            gc.collect()
        except BaseException:
            release(tag)
            raise
    _CACHE.pop('runner', None)
    return runs_table()


# ---------------------------------------------------------------------------
# score

def flies_mean_effect() -> float:
    """The flies' mean APL effect against all recordings (E108's 'all' target)."""
    from experiments import E108_paa_extrapolation as E108
    return float(E108.effect_targets()['all'].mean())


def _v(rd: pd.DataFrame, kind: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd['kind'] == kind) & (rd['unit'] == unit) & (rd['readout'] == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def targets() -> Tuple[str, ...]:
    from experiments import E108_paa_extrapolation as E108
    return ('all',) + tuple(E108.HALVES)


def collect() -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run: ({(kind, wiring): outputs}, runs table). Kinds: 'family', 'E110_A' (read), E125's PO50 and TY50
    (read), and this experiment's arms."""
    from experiments import E110_no_kckc as E110
    from experiments import E123_kckc_form as E123
    from experiments import E125_kckc_survivor as E125
    from experiments import E129_block_survivor_static as E129
    outs, rows = {}, []

    def add(kind, w, o, **meta):
        outs[(kind, w)] = o
        r = dict(kind=kind, wiring=w, **meta, **E83.realism(o), **E129.type_rates(o))
        r.update(E123.bands(r))
        rows.append(r)

    ra = E110.runs_table_of('A').set_index('tag')
    r125 = E125.runs_table().drop_duplicates('tag', keep='last').set_index('tag')
    for w in WIRINGS:
        F = E83.saved_fit(w)
        add('family', w, F['saved'], tag=F['tag'], folder=F['folder'], corner=False, **F['p'])
        a = E110.load_a(w)
        if a is not None:
            add('E110_A', w, a, tag=f'A_{w}', folder='E110_no_kckc/A', corner=False,
                **{k: float(ra.loc[f'A_{w}', k]) for k in E83.PARAM_KEYS})
        for arm in E125_READ:
            t = f'S1_{arm}_{w}'
            o = E125.load_run(t)
            if o is not None and t in r125.index:
                r = r125.loc[t]
                add(arm, w, o, tag=t, folder=f'{E125.EXPERIMENT}/step1', corner=bool(r['at_bound_thr'] or r['at_bound_apl']),
                    **{k: float(r[k]) for k in E83.PARAM_KEYS})
    for r in runs_table().to_dict('records'):
        o = load_run(r['tag'])
        if o is None:
            continue
        meta = {k: r[k] for k in ('tag', 'group', 'form', 'm', 'corner', 'at_bound_thr', 'at_bound_apl', 'held_on',
                                  'held_off', 'thr_factor', 'apl_factor', 'kckc_scale_off', 'host') if k in r}
        meta.update({k: r[k] for k in r if k.startswith(('mbc1_', 'share_', 'share3_'))})
        add(r['arm'], r['wiring'], o, folder=f'{EXPERIMENT}/step1', **{k: float(r[k]) for k in PARAM_KEYS}, **meta)
    runs = pd.DataFrame(rows)
    runs['corner'] = runs['corner'].fillna(False).astype(bool)
    runs['thr_below_60'] = runs['fixed_thr'].astype(float) < REALISM_THR_HZ
    return outs, runs


RUN_COLS = ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'rate_rho_off', 'opk_on_tv', 'opk_off_tv', 'max_count_on',
            'max_count_off', 'fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale', 'wKCKC_scale',
            'wMBC1_scale', 'mbc1_act_on', 'mbc1_act_off', 'mbc1_peak_on', 'mbc1_peak_off', 'share_mbc1_vs_apl_g',
            'share_mbc1_vs_apl_all', 'mbc1_per_synapse_vs_kckc') + tuple(
    f'{k}_{TYPE_COL[t]}' for t in TYPES for k in ('rate_on', 'rate_off', 'ratio', 'count_on', 'count_off'))


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


def readings(outs: Mapping[Tuple[str, str], Mapping[str, tuple]], runs: pd.DataFrame) -> pd.DataFrame:
    """Per kind, unit (each wiring with a run; 'S4', 'T4', 'all8' when all of the set have one) and readout: E110's
    `scores_of`, the runs' means of RUN_COLS ('run_<col>'), the number of runs outside the bands."""
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


def corners(runs: pd.DataFrame, kind: str) -> dict:
    """Reading 8: the wirings on which a kind's run is a corner, and whether there are more than CORNER_WIRINGS_MAX."""
    from experiments import E129_block_survivor_static as E129
    return E129.corners(runs, kind)


def below_family(rd: pd.DataFrame, runs: pd.DataFrame, arm: str, flies_mean: float) -> dict:
    """Reading 9's "below the family", each condition apart."""
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
    """Reading 9: the class of an arm ('below the family', 'equal', 'costs', 'neither', 'corner (no class)') with E123's
    flags, the parts and the realism beside."""
    from experiments import E123_kckc_form as E123
    c = E123.classify(rd, runs, arm)
    if not c.get('complete'):
        return c
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
    for col in ('fixed_thr', 'wAPLKC_scale', 'wKCKC_scale', 'wMBC1_scale', 'mbc1_act_on', 'mbc1_act_off',
                'share_mbc1_vs_apl_g', 'mbc1_per_synapse_vs_kckc', 'ratio_g', 'rate_on_g'):
        if col in sub:
            c[f'{col}_mean'] = float(pd.to_numeric(sub[col], errors='coerce').mean())
    return c


def differ(rd: pd.DataFrame, runs: pd.DataFrame, a: str, b: str) -> dict:
    """E126's "differ" rule (E129's `differ`, with this experiment's corners)."""
    from experiments import E129_block_survivor_static as E129
    return E129.differ(rd, runs, a, b)


def flags(classes: Mapping[str, dict], comps: Mapping[str, dict]) -> Dict[str, dict]:
    """"keeps the gain" (E123's flag, not a corner) and "stands in" (R50: E123's equal flag, not a corner, and R50
    lower than H50 by the "differ" rule)."""
    out = {}
    for a, c in classes.items():
        if not c.get('complete'):
            continue
        ok = not c.get('corner')
        f = dict(keeps_gain=bool(ok and c.get('keeps_gain')))
        if a == 'R50':
            f['stands_in'] = bool(ok and c.get('equal') and comps.get('R50|H50', {}).get('verdict') == 'R50 lower')
        out[a] = f
    return out


def outcomes(classes: Mapping[str, dict], comps: Mapping[str, dict], fl: Mapping[str, dict],
             runs: pd.DataFrame) -> dict:
    """Reading 10's flags; None where an arm they need is not complete."""
    def done(a):
        return bool(classes.get(a, {}).get('complete'))

    def cls(a):
        return classes.get(a, {}).get('cls')

    b_done = all(done(a) for a in B_ARMS)
    tolerated = None
    if b_done:
        for a in B_ARMS:
            if cls(a) != 'equal':
                break
            tolerated = a
    above = list(B_ARMS[B_ARMS.index(tolerated) + 1:]) if tolerated else []
    stands = fl.get('R50', {}).get('stands_in')
    keeps_x = fl.get('X50', {}).get('keeps_gain')
    forms = {a: comps.get(f'{a}|B100', {}).get('verdict') for a in FORM_ARMS if f'{a}|B100' in comps}
    mine = [a for a in classes if a in ARMS and done(a)]
    return dict(
        a=(bool(all(cls(a) == 'costs' for a in B_ARMS) and not stands and not keeps_x)
           if b_done and done('R50') and done('X50') else None),
        b=(bool(tolerated is not None and all(cls(a) == 'costs' for a in above)) if b_done else None),
        largest_tolerated_dose=None if tolerated is None else ARMS[tolerated]['m'],
        c=(bool(stands or keeps_x) if done('R50') or done('X50') else None),
        d_b_below=[a for a in B_ARMS if cls(a) == 'below the family'] if b_done else None,
        d=None,
        d_note='(d) needs step 2 (the fit) to hold on T; d_b_below lists the B doses below the family',
        e=forms, e_lower=[a for a, v in forms.items() if v == f'{a} lower'],
        f=bool(any(classes[a].get('corner_wirings') for a in mine)
               or bool((~runs.loc[runs['kind'].isin(mine), 'in_bands'].astype(bool)).any())),
        corner_wirings={a: classes[a].get('corner_wirings', []) for a in mine},
        outside_bands={a: classes[a].get('runs_outside_bands', []) for a in mine})


def step2_gate(classes: Mapping[str, dict], fl: Mapping[str, dict]) -> dict:
    """Reading 11: open if a B dose is "below the family", R50 "stands in" or X50 "keeps the gain" (none a corner).
    Step 2 and the permuted control are not written; they are written only when this opens."""
    why = [a for a in B_ARMS if classes.get(a, {}).get('cls') == 'below the family']
    if fl.get('R50', {}).get('stands_in'):
        why.append('R50 stands in')
    if fl.get('X50', {}).get('keeps_gain'):
        why.append('X50 keeps the gain')
    return dict(open=bool(why), opened_by=why, note='step 2 (the fit, E10_joint_fit.py with a free MB-C1 scale) and '
                                                    'the permuted control are not written yet')


def comparison_md(rd: pd.DataFrame, runs: pd.DataFrame, kinds: Sequence[str], classes: Mapping[str, dict],
                  fl: Mapping[str, dict]) -> str:
    from experiments import E123_kckc_form as E123
    lines = ['# E131 comparison tables (generated)', '',
             'Classes: ' + ', '.join(f"{a} {c.get('cls', 'incomplete')}" for a, c in classes.items()), '',
             '## The arms (calibrated readout, all recordings)', '',
             '| kind | m | pattern S4 | pattern T4 | d family S4 | d family T4 | amplitude S4 | shape S4 | spread ratio S4 | '
             'Pearson r S4 | ρ effect S4 | mean effect (readout) | mean effect (raw) | loss (raw) | MB-C1 Hz on / off | '
             'γ share vs APL | per synapse vs KC->KC | threshold | corners | outside bands | class | keeps gain | '
             'stands in |', '|---|' + '---:|' * 19 + '---|---|---|']
    for k in kinds:
        if not ((rd['kind'] == k) & (rd['unit'] == 'all8')).any():
            continue
        c = classes.get(k, {})
        sub = runs[runs['kind'] == k]

        def mean(col):
            return float(pd.to_numeric(sub[col], errors='coerce').mean()) if col in sub else np.nan
        f = fl.get(k, {})
        lines.append(
            f"| {k} | {ARMS[k]['m'] if k in ARMS else ''} | {_v(rd, k, 'S4', 'all_pattern'):.4f} | "
            f"{_v(rd, k, 'T4', 'all_pattern'):.4f} | {_v(rd, k, 'S4', 'all_pattern') - _v(rd, 'family', 'S4', 'all_pattern'):+.4f} | "
            f"{_v(rd, k, 'T4', 'all_pattern') - _v(rd, 'family', 'T4', 'all_pattern'):+.4f} | "
            f"{_v(rd, k, 'S4', 'all_amplitude'):.4f} | {_v(rd, k, 'S4', 'all_shape'):.4f} | "
            f"{_v(rd, k, 'S4', 'all_sd_ratio'):.3f} | {_v(rd, k, 'S4', 'all_pearson_effect'):.3f} | "
            f"{_v(rd, k, 'S4', 'all_rho_effect'):.3f} | {_v(rd, k, 'all8', 'mean_effect'):.3f} | "
            f"{_v(rd, k, 'all8', 'mean_effect', 'raw'):.3f} | {_v(rd, k, 'all8', 'all_loss', 'raw'):.4f} | "
            f"{mean('mbc1_act_on'):.2f} / {mean('mbc1_act_off'):.2f} | {mean('share_mbc1_vs_apl_g'):.3f} | "
            f"{mean('mbc1_per_synapse_vs_kckc'):.2f} | {mean('fixed_thr'):.1f} | "
            f"{len(c.get('corner_wirings', []))} | {len(c.get('runs_outside_bands', []))} | {c.get('cls', '')} | "
            f"{f.get('keeps_gain', '')} | {f.get('stands_in', '')} |")
    lines.append('')
    for u in ('S4', 'T4', 'all8'):
        lines += [f'## {u}', '', E123.comparison_table(rd, list(kinds), u), '']
    lines += ['## Runs and bands', '', E123.runs_bands_table(runs), '']
    return '\n'.join(lines)


def score() -> dict:
    from experiments import E129_block_survivor_static as E129
    out = out_dir()
    outs, runs = collect()
    rd = readings(outs, runs)
    rd.to_csv(out / 'readings.csv', index=False)
    runs.to_csv(out / 'runs.csv', index=False)
    kinds = list(dict.fromkeys(runs['kind']))
    fm = flies_mean_effect()
    classes = {a: classify(rd, runs, a, fm) for a in kinds if a not in ('family', 'E110_A')}
    comps = {f'{a}|{b}': differ(rd, runs, a, b) for a, b in COMPARISONS if a in kinds and b in kinds}
    fl = flags(classes, comps)
    st = _read_json(out / 'step0.json')
    v = dict(generated=_now(), host=_host(), checks_passed=check_passed(),
             threshold_check_wirings=sorted(check_wirings('threshold_passes')), half_way=E129.half_way_check(rd),
             flies_mean_effect=fm, classes=classes, comparisons=comps, flags=fl,
             outcomes=outcomes(classes, comps, fl, runs), step2_gate=step2_gate(classes, fl),
             step0=dict(threshold=st.get('threshold'), late_ratio=st.get('late_ratio'),
                        pooled_gate=st.get('pooled_gate'), dose_arithmetic=st.get('dose_arithmetic')),
             rules=dict(within=WITHIN, wirings_needed=N_WIRINGS_NEEDED, mean_tol=MEAN_TOL, raw_loss_tol=RAW_LOSS_TOL,
                        corner_wirings_max=CORNER_WIRINGS_MAX, differ_wirings=DIFFER_WIRINGS, rate_tol=RATE_TOL,
                        targets=list(targets()), populations=list(POPS), pooled_gate_r=POOLED_GATE_R),
             n_runs=int(len(runs)))
    _write_json(out / 'verdict.json', v)
    (out / 'comparison_table.md').write_text(comparison_md(rd, runs, [k for k in kinds], classes, fl))
    _write_json(out / 'provenance_score.json', dict(generated=_now(), host=_host(), script=SCRIPT,
                                                    commits=H.code_commits()))
    print(json.dumps(dict(classes={a: c.get('cls') for a, c in classes.items()}, flags=fl, outcomes=v['outcomes'],
                          step2_gate=v['step2_gate']), default=str), flush=True)
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (every E131 run, absolute parameters, the MB-C1 options included) and
    model_pair_distances_long.csv (raw spike counts) for the archive."""
    out = out_dir()
    rt = runs_table()
    rows, parts = [], []
    for r in rt.to_dict('records'):
        o = load_run(r['tag'])
        if o is None:
            continue
        arm = r['arm']
        rows.append(dict(
            tag=r['tag'], variant=r['tag'], part='S1', arm=arm, arm_label=ARM_LABELS[arm], form=r['form'], m=r['m'],
            base='bouton', k=E83.BREADTH_K, wiring=r['wiring'], ref_experiment=E83.members()[r['wiring']][0],
            ref_fit=E83.members()[r['wiring']][1],
            kind=f"C2's family's {r['wiring']}: {ARM_LABELS[arm]}, both response rates held",
            fixed_thr=r['fixed_thr'], wAPLKC_scale=r['wAPLKC_scale'], wKCAPL_scale=r['wKCAPL_scale'],
            wAPLPN_scale=r['wAPLPN_scale'], wPNAPL_scale=r['wPNAPL_scale'], wKCKC_scale=r['wKCKC_scale'],
            kc_kc_inh=r['wKCKC_scale'], kc_kc_tau=PLAIN_TAU_S, off_kckc_scale=r['kckc_scale_off'],
            kckc_matrix='hemibrain (connectome_KCKC_weights)', mbc1_inh=r['wMBC1_scale'], wMBC1_scale=r['wMBC1_scale'],
            mbc1_tau=r['mbc1_tau'], mbc1_thr=r['mbc1_thr'], mbc1_pooled=bool(r['mbc1_pooled']), mbc1_claws_only=False,
            multibouton_claws='first', thr_factor=r['thr_factor'], apl_factor=r['apl_factor'], rate_on=r['rate_on'],
            rate_off=r['rate_off'], max_count_on=r['max_count_on'], max_count_off=r['max_count_off'],
            realism_flag=r.get('realism_flag'), corner=bool(r['corner']), subdir='step1'))
        parts.append(pd.concat({(r['tag'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['tag'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out / 'model_pair_distances_long.csv', index=False)
    return df


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('check', 'step0', 'runs', 'score', 'variants', 'plan'))
    ap.add_argument('--wirings', nargs='+', default=None, choices=list(WIRINGS))
    ap.add_argument('--arms', nargs='+', default=None, choices=list(ARMS))
    ap.add_argument('--lane', default=None, help='the lane name in the tables (default the wirings joined by -)')
    ap.add_argument('--reclaim', action='store_true', help='take over claims left by a process that died')
    return ap


def main() -> None:
    a = build_parser().parse_args()
    ws = list(a.wirings or WIRINGS)
    if a.part == 'check':
        part_check(ws, a.lane)
    elif a.part == 'step0':
        part_step0()
    elif a.part == 'runs':
        part_runs(ws, a.lane, a.reclaim, a.arms)
    elif a.part == 'score':
        score()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'plan':
        st = _read_json(out_dir() / 'step0.json') or None
        ok = runnable_arms(st, check_wirings('threshold_passes'))
        for arm, w in planned(ws, a.arms):
            print(tag_of(arm, w), '' if arm in ok else '(not runnable yet)')


if __name__ == '__main__':
    main()
