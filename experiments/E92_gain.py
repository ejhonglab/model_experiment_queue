#!/usr/bin/env python3
"""E92, parts B to F: APL's inhibition of the KCs as a gain on each claw's drive.

The entry is item E92 of experiment_priorities.md (proposed by hal 2026-11-13 19:04; planned by a Fable subagent,
claude-fable-5-2, spawned from hal 2026-11-13 18:38). Part A (no simulations) is experiments.E92_steps_like_for_like;
this module holds the parts that simulate. olfsysm's `apl_inh_divisive` (`mp.kc.apl_inh_divisive`) divides each claw's
PN drive x by 1 + max(i, 0), i its APL inhibition (APL->claw weight times APL's activity), in place of the clipped
subtraction max(x - i, 0); `bouton_inh_divisive` (`mp.pn.bouton_inh_divisive`) does the same at the PN boutons. With APL
off both are the subtractive runs (x / (1 + 0) = x - 0), so every APL-off run is the saved fit's (checked per wiring and
form).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`, through E91's `WiringModel91`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the
recordings E77's (`recording_halves()`, seed 0), C / D the second split's (seed 1; E89's `targets89`). Every point is on
the fits as they are (one arm, 'B'); its APL->KC weights are E89's at the family's measure F and the point's exponent
(`apl_kc_weights`; F at exponent 2 gives the fit's weights exactly), and the tied APL scale (APL->KC, KC->APL,
bouton->APL; E83's TIED_ON) is bisected to the wiring's APL-on count (E89's `weighting_reading`, E69's `bisect_factor`)
in APL_BRACKET, wider than E83's because the scale changes units under a divisor; a point at a bound is reported
(`at_bound`) and left out of the choice.

The reference (`--reference`, default 'C2'): C2's clipped family, E89's saved C1 points at the same weightings
(variants B_F_k2 = the unchanged fits, B_F_k1, B_F_k0.5, B_none) and on T the unchanged fits. If E91's refitted family
passes as a population of eight, the reference is that family ('E91_refits'): a stub (`reference_rows_e91_refits`),
since olfsysm refuses the divisive form together with `claw_inh_pool`, so that branch compares against the pooled form's
points at the refits' parameters, and E91's runs are not read before its note is posted.

Parts (CLI `--part`):
- B1: `apl_inh_divisive=True` at F with exponents 2, 1, 0.5 and no weighting on S (16 points); runs_B1.csv.
- B2, for the record (E91's outcome (e)): `bouton_inh_divisive=True` at exponent 2, and both divisive forms at exponents
  2 and 1 (12 points); runs_B2.csv. Read and reported beside B1, not chosen from.
- B3: only if part A's line (A4) names a variant off the line that no gate has read (at most 8 points on S); a stub.
- smoke: B1's point at exponent 2 on the real wiring alone, into smoke/ (runs_smoke.csv and smoke.json; not runs_B1.csv):
  the simulator's identity at the fit, the APL-off run with the divisive form equal to the saved fit's, the bisection
  inside the bracket, the parameter set in place equal to a simulator built with it, the flags reset after.
- co: the co-response changes under APL of the four nearest ester pairs, the ester and the alcohol groups beside every S
  setting (E89's `coresponse_c` readings 1 and 1b, and the geometric mean of E88's two estimates at L1 and L3 with its
  known error from E89's estimator check); coresponse_B.csv.
- C: the choice per half of both splits (`choose92`: the lowest pattern part of S's population under the calibrated
  readout among B1's settings and the references, inside the bands and not at a bracket's bound); a divisive choice that
  gains CHOICE_MIN_GAIN or more over the reference is run on T and read with E83's gate against the other half in four
  directions (`gate92`, E91's `gate91`); the control on S and T (`control_matrix92`: the reference's subtractive form at
  the chosen weighting with a static inhibition at the soma, olfsysm's `ext_inh_kc`, present only with APL on, equal for
  every odor to each KC's mean loss of drive under the gain at the chosen setting, the APL scale re-held); then `score`.
- score: readings_B.csv (every setting per unit and readout, with the amplitude part, the shape part, the spread ratio
  and Pearson r beside every pattern part), beside_B.csv (part A's readings beside every S setting: the two groups'
  distances with APL on in A2's readings and the place against A4's line, when part A has them), coresponse_B.csv
  (unless --no-co), choice.csv, gate.csv, equal.csv, verdict_B.json, comparison_table_B.md, E92_gain_screen.png.
- D: the family fitted again, only on a choice by C or an equal of B1's (`d_plan`): `E10_joint_fit.py -k <exponent>
  --aplpn-max 3 --rate-target ctrl_frac --loss apl --apl-readout both --tag _cap3 --apl-divisive`, 150 evaluations, from
  each wiring's own fit of C2's family with its tied APL scale at B1's re-held factor for that exponent (`--init-params`
  on d_start_params.json; the point is run first on a wiring that lacks it; `--d-start exact` starts at the fit's own
  scales instead, as E91's D did), one process at a time; `--pop S` first, `--pop T` only if S's fitted population is
  D_S_GAIN or more below the reference against all recordings. `--part Dscore` scores them (d_scores.csv, .json).
- E: only on a pass by C or D (a stub).
- F: the summed form (E91, its note posted 2026-11-13 19:38) read as part A read C2's family, on E91's archived runs:
  E91's refitted family (part D: the four fits on S at s = 1, exponent 2), the summed form on the fits as they are
  (B1_F_k2_s1, S), its choice on the base with the flies' recruitment (B2_F_k1_s1, S and T) and that base's own family
  (B2_F_k2: E89's C2 runs on S, E91's on T), with C2's clipped family (E90's runs) as the reference (F_FAMILIES). F1:
  their somata in part A's A2 readings (`somata_pairs`), compared with the flies by part A's rule (6 of 8 wirings on the
  side for a population of eight, 3 of 4 for S's population), somata_F.csv and somata_F_means.csv. F2: every one of
  E91's archived variants (and E89's base runs on S) placed against A4's line (part A's `variant_rows`: S and T
  populations, the pattern part with its amplitude and shape parts, the spread ratio, Pearson r and ρ), line_F.csv. F3:
  the refits and the B2 choice on S run again with the stage matrices saved (E90's reading; `claw_inh_pool` 1; the B2
  choice odor by odor at its per-odor thresholds and APL scales, `FPool`), each reproducing its archived counts and
  responses exactly (identity_F.json), the matrices under outputs/.../F/ (the counts stay the archived run's). F4: their
  claws in E90's readings (`read_wiring`) and the two steps like for like (part A's `step_verdict_rows` on S),
  claws_F.csv and steps_F.csv. Then verdict_F.json, comparison_table_F.md, E92_F_line.png. `--part Ftable` rebuilds the
  last three.
- variants: variants.csv and model_pair_distances_long.csv for the archive: `./venv/bin/python
  archive_tools/archive_model_outputs.py E92_steps_like_for_like`, then `./venv/bin/python
  archive_tools/migrate_model_output_params.py --apply`.
- table: comparison_table_B.md and the figure.

Every run row carries the absolute parameters (the tied APL scale's factor and the scales as run), the form, the
weighting and the spread of its weights over the responding KCs (with APL off, E89's `factor_summary`, and with APL on),
the realism (E83's `realism`: the four bands under the rule of 2026-11-10 16:10, `realism_flag`, with the earlier rule
beside, `realism_flag_zcrit`), the largest spike counts, the divisor 1 + i at the responders' claws during the pulse
(10th percentile, median, 90th percentile, largest, and the share above DIVISOR_HIGH; at the boutons feeding them for
the bouton form), the loss of drive over the responders, and `H.code_commits()` of the process that ran it (which must
be clean).

Outputs: results/E92_steps_like_for_like/; spike counts and responses under outputs/E92_steps_like_for_like/B/ (each
run's loss of drive per KC and odor beside its counts, <tag>_loss.parquet); D's fits where E10_joint_fit.py writes them
for this experiment name (fits.csv in the results folder, counts at the outputs folder's root, logs under D/).

Usage (from the repository root; one 16-thread process per lane, two lanes at most on hal):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E92_gain --part B1 --wirings S
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E87_inhibition_form as E87
from experiments import E88_responder_sets as E88
from experiments import E89_apl_breadth_weighting as E89
from experiments import E91_summed_drive_inhibition as E91

EXPERIMENT = 'E92_steps_like_for_like'
ARM = 'B'                                           # every point is on the fits as they are

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

FORMS = {'sub': (False, False), 'div': (True, False), 'bdiv': (False, True), 'both': (True, True)}
FORM_FLAGS = ('apl_inh_divisive', 'bouton_inh_divisive')
B1_SETTINGS = ((('F', 2.0), 'div'), (('F', 1.0), 'div'), (('F', 0.5), 'div'), (('F', 0.0), 'div'))
B2_SETTINGS = ((('F', 2.0), 'bdiv'), (('F', 2.0), 'both'), (('F', 1.0), 'both'))
REF_WEIGHTINGS = tuple(w for w, _ in B1_SETTINGS)   # E89's C1 points at these are the references
REFERENCE = E89.REFERENCE                           # ('F', 2.0): the family's weighting
REFERENCE_FAMILIES = ('C2', 'E91_refits')
RAW_CRIT = E87.RAW_CRIT                             # a model response is two or more spikes
PULSE = E91.PULSE                                   # the odor pulse (s; [start, end))
# The tied APL scale's bracket. olfsysm's tests tune divisive inhibition to ~5-15x smaller APL scales than subtractive
# (each scale tuned; test_mb_model's _DIVISIVE_TUNING_KWS). Here one factor multiplies the three tied scales (APL->KC,
# KC->APL, bouton->APL), so the inhibition at a claw goes about with its square (15x less inhibition: a factor near
# 0.26); 0.01 covers 15x even if the inhibition went linearly with the factor (0.067), with a margin of more than 6.
# The top is E83's. Points at a bound are reported (at_bound) and are not choices.
APL_BRACKET = (0.01, 10.0)
APL_STEP = E83.APL_STEP
DIVISOR_HIGH = 10.0                                 # the entry's realism note: a shunt above about tenfold
CO_CANDIDATE = E91.CO_CANDIDATE                     # the geometric mean of E88's two estimates of the co-response
CHOICE_SCREENS = ('ref', 'B1')                      # the choice is among the references and B1's settings

# part C (the entry's rules: E91's, whose gate91 reads its own module's thresholds)
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
CONTROL_RATIO = E91.CONTROL_RATIO
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN, EQUAL_MAX_EXPONENT = (E91.EQUAL_PATTERN, E91.EQUAL_RAW, E91.EQUAL_MEAN,
                                                            E91.EQUAL_MAX_EXPONENT)
HALVES = E91.HALVES
SPLIT_OF = E91.SPLIT_OF

# part D
D_MAXFEV = 150
D_STARTS = ('reheld', 'exact')
D_S_GAIN = 0.0015                                   # S's fitted population this far below the reference -> T's fits
D_PASS_GAIN = 0.003                                 # the population of eight this far below UNCHANGED_8
D_WIRINGS_LOWER = 6                                 # of the eight
D_RAW_REF, D_RAW_SLACK, D_MEAN_SLACK = 0.073, 0.005, 0.01
D_EXPONENTS = (2.0, 1.0, 0.5, 0.0)
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
FLIES_MEAN = E89.FLIES_MEAN

# the names part A's module is expected to provide (experiments.E92_steps_like_for_like), read beside every point
PART_A_NAMES = ('somata_readings', 'line', 'place_on_line', 'relative_distances')

_CACHE: Dict[tuple, object] = {}


def _write_json(path: Path, obj) -> None:
    E89._write_json(path, obj)


# ---------------------------------------------------------------------------
# names

def variant_name(w: Tuple[str, float], form: str = 'sub', control: bool = False) -> str:
    """A setting's variant name: 'B_<weighting>' (E89's `wtag`) for the subtractive form, so that E89's references read
    as 'B_F_k2' etc.; '_<form>' for a divisive form ('div', 'bdiv', 'both'); '_ctrl' for the control of the setting
    named before it (which itself runs in the subtractive form)."""
    if form not in FORMS:
        raise ValueError(f'unknown form {form!r}')
    if control and form == 'sub':
        raise ValueError('a control belongs to a divisive setting')
    v = f'{ARM}_{E89.wtag((w[0], float(w[1])))}'
    if form != 'sub':
        v += f'_{form}'
    if control:
        v += '_ctrl'
    return v


def parse_variant(v: str) -> dict:
    """`variant_name`'s fields back: arm, measure ('none' for no weighting), exponent, weighting ((measure, exponent),
    E89's NO_WEIGHTING for none), wtag, form (as run: 'sub' for a control), form_of (the form in the name), the two
    flags as run (apl_inh_divisive, bouton_inh_divisive), control, screen ('ref', 'B1', 'B2' or 'control')."""
    arm, _, rest = v.partition('_')
    if arm != ARM or not rest:
        raise ValueError(f'variant {v!r}: not one of this experiment\'s')
    parts = rest.split('_')
    if parts[0] == 'none':
        w, i = E89.NO_WEIGHTING, 1
    else:
        if len(parts) < 2 or not re.fullmatch(r'k[0-9.]+', parts[1]):
            raise ValueError(f'variant {v!r}: no weighting')
        w, i = (parts[0], float(parts[1][1:])), 2
    form, control = 'sub', False
    for p in parts[i:]:
        if p in FORMS and p != 'sub' and form == 'sub' and not control:
            form = p
        elif p == 'ctrl' and form != 'sub' and not control:
            control = True
        else:
            raise ValueError(f'variant {v!r}: unknown part {p!r}')
    run_form = 'sub' if control else form
    screen = 'control' if control else ('ref' if form == 'sub' else ('B1' if form == 'div' else 'B2'))
    return dict(arm=arm, measure='none' if w[1] == 0 else w[0], exponent=float(w[1]), weighting=w, wtag=E89.wtag(w),
                form=run_form, form_of=form, apl_inh_divisive=FORMS[run_form][0], bouton_inh_divisive=FORMS[run_form][1],
                control=control, screen=screen)


def control_of(variant: str) -> str:
    """The control's name of a divisive setting."""
    pv = parse_variant(variant)
    return variant_name(pv['weighting'], pv['form_of'], control=True)


def ref_variant() -> str:
    """The family's reference: F at exponent 2, subtractive (B_F_k2, the unchanged fits)."""
    return variant_name(REFERENCE)


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


def b1_variants() -> List[str]:
    return [variant_name(w, f) for w, f in B1_SETTINGS]


def b2_variants() -> List[str]:
    return [variant_name(w, f) for w, f in B2_SETTINGS]


def _wirings(arg: str) -> Tuple[str, ...]:
    return E91._wirings(arg)


# ---------------------------------------------------------------------------
# the divisive forms on a simulator

def form_kws(form: str) -> dict:
    """`fit_mb_model` kwargs of a form (only the True flags, as fit_mb_model records them)."""
    return {k: True for k, f in zip(FORM_FLAGS, FORMS[form]) if f}


def set_form(mp, form: str = 'sub') -> None:
    """Set the divisive flags of `form` (FORMS) on a simulator's olfsysm parameters (they stay set for its later
    simulations; `set_form(mp, 'sub')` resets them). Refuses the divisive KC form with a pooled share set (olfsysm
    refuses the two together, inside its parallel region) and the bouton form without preset APL->bouton weights."""
    if form not in FORMS:
        raise ValueError(f'unknown form {form!r}')
    kc, bt = FORMS[form]
    if kc and float(getattr(mp.kc, 'claw_inh_pool', 0.0)) != 0.0:
        raise RuntimeError('apl_inh_divisive with claw_inh_pool set: olfsysm refuses the two together')
    if bt and not bool(mp.pn.preset_wAPLPN):
        raise RuntimeError('bouton_inh_divisive needs the boutons\' APL weights (prat_boutons)')
    mp.kc.apl_inh_divisive = bool(kc)
    mp.pn.bouton_inh_divisive = bool(bt)


def get_form(mp) -> str:
    """The form a simulator's olfsysm parameters hold."""
    flags = (bool(mp.kc.apl_inh_divisive), bool(mp.pn.bouton_inh_divisive))
    return next(k for k, v in FORMS.items() if v == flags)


# ---------------------------------------------------------------------------
# the claws' loss of drive and the divisor, from a run's saved bouton rates and APL activity

def claw_loss(x: np.ndarray, inh: np.ndarray, kc_div: bool) -> Tuple[np.ndarray, np.ndarray]:
    """(claw state, loss of drive) from the claws' PN drive `x` and APL inhibition `inh` (APL->claw weight times APL's
    activity at the step before), as olfsysm forms the claw state: divisive max(x / (1 + max(i, 0)), 0), whose loss
    x - x / (1 + i) = state * i; subtractive max(x - i, 0), loss max(x, 0) - state."""
    x = np.asarray(x, dtype=float)
    inh = np.asarray(inh, dtype=float)
    if kc_div:
        i = np.maximum(inh, 0.0)
        state = np.maximum(x / (1.0 + i), 0.0)
        return state, state * i
    state = np.maximum(x - inh, 0.0)
    return state, np.maximum(x, 0.0) - state


QUANTILES = (('p10', 0.1), ('p50', 0.5), ('p90', 0.9), ('max', 1.0))


def _quantile_row(prefix: str, pool: List[np.ndarray], high: Optional[float] = None) -> dict:
    v = np.concatenate(pool) if pool else np.zeros(0)
    if not v.size:
        row = {f'{prefix}_{n}': np.nan for n, _ in QUANTILES}
        if high is not None:
            row[f'{prefix}_share_above{high:g}'] = np.nan
        return row
    q = np.quantile(v, [x for _, x in QUANTILES])
    row = {f'{prefix}_{n}': float(x) for (n, _), x in zip(QUANTILES, q)}
    if high is not None:
        row[f'{prefix}_share_above{high:g}'] = float((v > high).mean())
    return row


def gain_stats(rv, mp, odors: Sequence[int], responders: np.ndarray, form: str) -> dict:
    """From the last APL-on run of a simulator (APL's activity saved: `mp.kc.save_inh_sims` was on), per KC and odor
    (`odors`: olfsysm odor indices, one per output column) the window means over PULSE of the loss of drive ('L',
    `claw_loss` summed over the KC's claws) and of the summed drive before the KCs' inhibition ('D', sum of max(x, 0));
    and, over the pulse's time steps at the claws of the `responders` (KCs x odors), the summary 'row': i (inh_kc_*),
    the divisor 1 + i when the KC form is divisive (divisor_kc_*), the divisor at the boutons feeding those claws when
    the bouton form is (divisor_bouton_*), the mean loss, drive and loss over drive over the responders. x = wPNKC @
    bouton rates(t) (the saved rates, after the boutons' own inhibition), i = wAPLKC * APL(t - 1), as olfsysm forms them."""
    from pyolfsysm.mb_model import get_time_index
    kc_div, bt_div = FORMS[form]
    t = np.asarray(get_time_index(mp, rv), dtype=float)
    Wp = np.asarray(rv.kc.wPNKC, dtype=float)
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    c2k = np.asarray(rv.kc.claw_to_kc, dtype=int).ravel()
    n_kc = int(mp.kc.N)
    if Wp.shape[0] != c2k.size or wk.size != c2k.size:
        raise RuntimeError(f'{Wp.shape[0]} claws in wPNKC, {wk.size} APL>claw weights, {c2k.size} claw KCs')
    responders = np.asarray(responders, dtype=bool).reshape(n_kc, len(odors))
    ih = rv.kc.inh_sims
    if len(ih) == 0 or np.asarray(ih[int(odors[0])]).size == 0:
        raise RuntimeError("APL's activity was not saved (mp.kc.save_inh_sims)")
    wpn = None
    if bt_div:
        wpn = np.asarray(rv.pn.wAPLPN, dtype=float).ravel()
        if wpn.size != Wp.shape[1]:
            raise RuntimeError(f'{wpn.size} APL>bouton weights for {Wp.shape[1]} boutons')
    cols = np.flatnonzero((t >= PULSE[0] - 1e-12) & (t < PULSE[1] - 1e-12))
    prev = np.maximum(cols - 1, 0)
    A = E91._kc_sum(c2k, n_kc)
    shape = (n_kc, len(odors))
    out = {'L': np.full(shape, np.nan), 'D': np.full(shape, np.nan)}
    pools: Dict[str, List[np.ndarray]] = {'inh': [], 'bouton': []}
    for j, oi in enumerate(odors):
        bt = np.asarray(rv.pn.bouton_sims[int(oi)], dtype=float)
        a = np.asarray(ih[int(oi)], dtype=float).reshape(-1, t.size)[0][prev]
        x = Wp @ bt[:, cols]
        inh = wk[:, None] * a[None, :]
        _, loss = claw_loss(x, inh, kc_div)
        out['L'][:, j] = np.asarray(A @ loss).mean(axis=1)
        out['D'][:, j] = np.asarray(A @ np.maximum(x, 0.0)).mean(axis=1)
        rk = np.flatnonzero(responders[:, j])
        if rk.size:
            rc = np.flatnonzero(np.isin(c2k, rk))
            pools['inh'].append(np.maximum(inh[rc], 0.0).ravel())
            if bt_div:
                bts = np.flatnonzero((Wp[rc] != 0).any(axis=0))
                pools['bouton'].append((1.0 + np.maximum(wpn[bts][:, None] * a[None, :], 0.0)).ravel())
        del x, inh, loss
    row = _quantile_row('inh_kc', pools['inh'])
    if kc_div:
        row.update(_quantile_row('divisor_kc', [1.0 + p for p in pools['inh']], DIVISOR_HIGH))
    else:
        row.update(_quantile_row('divisor_kc', [], DIVISOR_HIGH))
    row.update(_quantile_row('divisor_bouton', pools['bouton'] if bt_div else [], DIVISOR_HIGH))
    with np.errstate(divide='ignore', invalid='ignore'):
        ld = np.where(out['D'] > 0, out['L'] / out['D'], np.nan)
    resp = responders
    row.update(loss_resp_mean=float(np.nanmean(out['L'][resp])) if resp.any() else np.nan,
               drive_resp_mean=float(np.nanmean(out['D'][resp])) if resp.any() else np.nan,
               loss_over_drive_resp=float(np.nanmean(ld[resp])) if resp.any() else np.nan,
               loss_all_mean=float(np.nanmean(out['L'])))
    out['row'] = row
    return out


def control_matrix92(loss: pd.DataFrame) -> pd.DataFrame:
    """The control's static inhibition at the soma (olfsysm `ext_inh_kc`, a KCs x odors matrix): per KC its mean over the
    odors of the loss of drive under the gain (`loss`: KCs x odors, window means over the pulse of the KC's claws' loss,
    the chosen setting's run), the same for every odor (E91's `control_matrix` at s = 1). It keeps how much each KC is
    inhibited and removes the dependence on what the odor drives through it."""
    return E91.control_matrix(loss, 1.0)


# ---------------------------------------------------------------------------
# saving and loading runs

def _b_dir(sub: str = 'B') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_run(tag: str, o: Dict[str, tuple], loss: Optional[pd.DataFrame] = None, sub: str = 'B') -> None:
    d = _b_dir(sub)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')
    if loss is not None:
        loss.to_parquet(d / f'{tag}_loss.parquet')


def _dirs(sub: str) -> Tuple[Path, Path]:
    return _b_dir(sub), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / sub


def load_run(tag: str, sub: str = 'B') -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment ({'on', 'off'}: (responses, spike counts)): the local outputs, else the archive."""
    for d in _dirs(sub):
        out = {}
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            ps = [d / f'{tag}_{cc}_{k}.parquet' for k in ('responses', 'spike_counts')]
            if not all(p.exists() for p in ps):
                out = None
                break
            out[c] = tuple(pd.read_parquet(p) for p in ps)
        if out:
            return out
    return None


def load_loss(tag: str, sub: str = 'B') -> pd.DataFrame:
    for d in _dirs(sub):
        p = d / f'{tag}_loss.parquet'
        if p.exists():
            return pd.read_parquet(p)
    raise RuntimeError(f'no loss of drive saved for {tag}')


def load_any(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E92' this experiment's, 'E89' E89's C runs (the references), 'saved' the wiring's fit."""
    if source == 'saved':
        return E83.saved_fit(wiring)['saved']
    o = E89.load_run(tag) if source == 'E89' else load_run(tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


def _append(path: Path, row: dict) -> None:
    E83._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def _frame(a: np.ndarray, like: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(a, index=like.index, columns=like.columns)


def _commits(strict: bool = True) -> dict:
    """H.code_commits() at the start of a simulation part, which must be clean (E89's `_commits_clean`); with `strict`
    False (the smoke run only) a dirty checkout is recorded, with its changed files, instead of refused."""
    if strict:
        return E89._commits_clean()
    cc = dict(H.code_commits())
    if cc.get('queue_dirty'):
        st = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=H.HERE, capture_output=True,
                            text=True).stdout
        cc['queue_dirty_files'] = ' '.join(ln[3:] for ln in st.splitlines())
    return cc


# ---------------------------------------------------------------------------
# B: the settings on the fits as they are

def off_check(model, form: str, sub: str = 'B') -> dict:
    """With `form` set, the APL-off run (all APL scales 0) equals the saved fit's (x / (1 + 0) = x - 0); cached per
    wiring and form (outputs/.../<sub>/off_check_<wiring>_<form>.json); raises otherwise."""
    path = _b_dir(sub) / f'off_check_{model.label}_{form}.json'
    if path.exists():
        return json.loads(path.read_text())
    set_form(model.mp, form)
    try:
        off = model.sim(model.p_off())
    finally:
        set_form(model.mp, 'sub')
    ns, nr = E83._same({'on': off, 'off': off}, {'on': model.saved['off'], 'off': model.saved['off']})
    d = dict(wiring=model.label, form=form, count_diffs=ns // 2 if ns >= 0 else ns, response_diffs=nr // 2 if nr >= 0 else nr)
    if ns or nr:
        raise RuntimeError(f'{model.label}: the APL-off run with form {form} differs from the saved one ({d})')
    path.write_text(json.dumps(d) + '\n')
    return d


def b_reading92(model, wk: Optional[pd.Series], form: str, ext: Optional[np.ndarray] = None, label: str = '',
                bracket: Tuple[float, float] = APL_BRACKET) -> dict:
    """One setting on a wiring at its fit: the APL->KC weights `wk` (None: the fit's), the form's divisive flags, the
    control's static inhibition `ext` (per KC, the outputs' row order; None for none); the tied APL scale (E83's TIED_ON)
    bisected to the fit's APL-on count (E69's `bisect_factor` in `bracket`, E83's APL_STEP, from 1); then the run again
    with APL's activity saved for `gain_stats` (its counts must equal). The APL-off run is the saved fit's. Everything
    is reset on return, also on an error."""
    from experiments import E69_kc_firing_ceiling as E69
    t0 = time.time()
    E89.set_weights(model, wk)
    set_form(model.mp, form)
    E91.set_ext_inh(model.mp, model.rv, ext, len(model.order))
    try:
        x, on, hist, at_bound = E69.bisect_factor(model.run_on, model.target_on, model.size, bracket, APL_STEP,
                                                  f'{label} APL scale')
        model.mp.kc.save_inh_sims = True
        try:
            on2 = model.run_on(x)
        finally:
            model.mp.kc.save_inh_sims = False
        ns = int((on2[1].to_numpy(dtype=float) != on[1].to_numpy(dtype=float)).sum())
        if ns:
            raise RuntimeError(f"{model.label}: the run again with APL's activity saved differs in {ns} counts")
        st = gain_stats(model.rv, model.mp, model.order, on[1].to_numpy(dtype=float) > RAW_CRIT, form)
    finally:
        E91.set_ext_inh(model.mp, model.rv, None)
        set_form(model.mp, 'sub')
        E89.set_weights(model, None)
    p = model.p
    return dict(on=on, off=model.saved['off'], apl_factor=float(x), at_bound_apl=bool(at_bound), n_sims=len(hist) + 1,
                params=dict(p, **{k: p[k] * x for k in E83.TIED_ON}), seconds=time.time() - t0, history=hist, stats=st,
                bracket=tuple(bracket))


def run_point(model, v: str, sub: str = 'B', part: str = '', cc: Optional[dict] = None) -> Tuple[dict, dict]:
    """One setting (variant name `v`) on a wiring's model at its fit (`b_reading92`), its counts and loss of drive saved
    under outputs/.../<sub>/; a divisive form's `off_check` first; a control (`..._ctrl`) runs in the subtractive form at
    its setting's weighting with the static inhibition `control_matrix92` of its setting's run on the same wiring (which
    must exist). Returns (the run's row, `b_reading92`'s result)."""
    wiring = model.label
    tag = run_tag(v, wiring)
    pv = parse_variant(v)
    cc = H.code_commits() if cc is None else cc
    resp_off = E89._responders(model.saved['off'][1])
    offd = off_check(model, pv['form'], sub) if pv['form'] != 'sub' else {}
    wk, fac = E91._setting_weights(model, pv['weighting'])
    ext, src = None, ''
    if pv['control']:
        src = variant_name(pv['weighting'], pv['form_of'])
        cm = control_matrix92(load_loss(run_tag(src, wiring), sub)).reindex(model.saved['on'][1].index)
        if cm.isna().any().any():
            raise RuntimeError(f"{tag}: the control matrix does not cover the model's KCs")
        ext = cm.iloc[:, 0].to_numpy(dtype=float)
    res = b_reading92(model, wk, pv['form'], ext=ext, label=tag)
    o = {'on': res['on'], 'off': res['off']}
    save_run(tag, o, _frame(res['stats']['L'], res['on'][1]), sub)
    resp_on = E89._responders(res['on'][1])
    fac_on = {k.replace('fac_resp_', 'fac_respon_'): x for k, x in E89.factor_summary(fac, resp_on).items()
              if k.startswith('fac_resp_')}
    row = dict(part=part, arm=ARM, tag=tag, variant=v, source='E92', screen=pv['screen'], weighting=pv['wtag'],
               measure=pv['measure'], exponent=pv['exponent'], form=pv['form'], form_of=pv['form_of'],
               apl_inh_divisive=pv['apl_inh_divisive'], bouton_inh_divisive=pv['bouton_inh_divisive'],
               control_of=src, wiring=wiring, permute_seed=np.nan, fit_folder=model.folder, fit_tag=model.tag,
               apl_factor=res['apl_factor'], apl_min=res['apl_factor'], apl_max=res['apl_factor'],
               at_bound=res['at_bound_apl'], bracket_lo=res['bracket'][0], bracket_hi=res['bracket'][1],
               n_on=E83._n(res['on'][0]), target_on=model.target_on, n_off=E83._n(res['off'][0]),
               target_off=model.target_off, off_check_count_diffs=offd.get('count_diffs', np.nan),
               n_sims=res['n_sims'], seconds=res['seconds'], **{k: res['params'][k] for k in E83.PARAM_KEYS},
               wAPLKC_sha256=E89.E80_sha(wk), K_fit=model.K,
               ext_inh_kc_max=float(ext.max()) if ext is not None else np.nan,
               ext_inh_kc_mean=float(ext.mean()) if ext is not None else np.nan,
               **E89.factor_summary(fac, resp_off), **fac_on, **res['stats']['row'], **cc)
    row.update(E89._realism_row(o))
    print(f'E92 {tag}: APL scale x {res["apl_factor"]:.4g}{" (at a bound)" if res["at_bound_apl"] else ""}, '
          f'{row["n_on"]} / {model.target_on} responses, {res["n_sims"]} sims, {res["seconds"]:.0f}s; divisor at '
          f'the responders\' claws p50 {row["divisor_kc_p50"]:.3g}, max {row["divisor_kc_max"]:.3g}; '
          f'flag {row.get("realism_flag")!r}', flush=True)
    return row, res


def b_points(out: Path, wirings: Sequence[str], variants: Sequence[str], runs_name: str, part: str, sub: str = 'B',
             check: bool = False) -> pd.DataFrame:
    """Settings (variant names) on the fits as they are (B1, B2, C's runs, D's starts), each by `run_point`. With
    `check`, the simulator's identity at the fit (E83's `_check_saved`). Resumable (rows in `runs_name`), runs under
    outputs/.../<sub>/."""
    import gc
    cc = _commits()
    path = out / runs_name
    for wiring in wirings:
        todo = [v for v in variants if run_tag(v, wiring) not in _done(path)]
        if not todo:
            continue
        model = E91.WiringModel91(wiring, check=check)
        for v in todo:
            row, _ = run_point(model, v, sub, part, cc)
            _append(path, row)
        del model
        gc.collect()
    return pd.read_csv(path)


def part_b1(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B1: apl_inh_divisive at F with exponents 2, 1, 0.5 and no weighting (16 points on S)."""
    return b_points(out, wirings, b1_variants(), 'runs_B1.csv', 'B1')


def part_b2(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B2, for the record: bouton_inh_divisive at exponent 2, both forms at exponents 2 and 1 (12 points on S)."""
    return b_points(out, wirings, b2_variants(), 'runs_B2.csv', 'B2')


def part_b3(out: Path, wirings: Sequence[str]) -> pd.DataFrame:
    """B3: only if part A's line (A4, line_variants.csv) names a variant off the line toward the flies that no gate has
    read, at most 8 points of it on S, at fixed scales with the tied APL scale re-held, read as B1's (a stub until A4
    names one)."""
    raise NotImplementedError('E92 part B3 runs only if part A4 names a variant off the line that no gate has read; '
                              'not written')


def part_smoke(out: Path) -> dict:
    """The smoke run: B1's point at exponent 2 on the real wiring alone (outputs/.../smoke/, results' runs_smoke.csv and
    smoke.json): the simulator's identity at the fit (E83's `_check_saved`), the APL-off run with the divisive form equal
    to the saved fit's (`off_check`), the bisection inside the bracket, the model's flags, weights and static inhibition
    reset after the point, and the point again on a simulator built with `apl_inh_divisive=True` (its counts equal to
    the in-place run's, its APL-off run the saved fit's). A dirty checkout is recorded, not refused."""
    import gc
    t0 = time.time()
    v = variant_name(REFERENCE, 'div')
    tag = run_tag(v, 'real')
    cc = _commits(strict=False)
    m = E91.WiringModel91('real', check=True)
    row, res_p = run_point(m, v, 'smoke', 'smoke', cc)
    path = out / 'runs_smoke.csv'
    if path.exists():
        path.unlink()
    _append(path, row)
    lo, hi = APL_BRACKET
    res = dict(tag=tag, apl_factor=row['apl_factor'], bracket=[lo, hi], at_bound=bool(row['at_bound']),
               inside_bracket=bool(not row['at_bound'] and lo < row['apl_factor'] < hi), n_on=int(row['n_on']),
               target_on=int(row['target_on']), off_check_count_diffs=row['off_check_count_diffs'],
               seconds_point=float(row['seconds']), bisection=res_p['history'],
               form_after=get_form(m.mp), weights_reset=bool(np.array_equal(m._kc, m._kc0)),
               ext_inh_after=int(np.asarray(m.rv.kc.ext_inh).size), save_inh_sims_after=bool(m.mp.kc.save_inh_sims))
    del m
    gc.collect()
    mk = E91.WiringModel91('real', check=False, apl_inh_divisive=True)
    res['built_flag'] = bool(mk.mp.kc.apl_inh_divisive)
    on_b, off_b = mk.run_on(float(row['apl_factor'])), mk.sim(mk.p_off())
    o = load_run(tag, 'smoke')
    res['built_on_count_diffs'] = int((on_b[1].to_numpy(dtype=float) != o['on'][1].to_numpy(dtype=float)).sum())
    res['built_off_count_diffs'] = int((off_b[1].to_numpy(dtype=float) != mk.saved['off'][1].to_numpy(dtype=float)).sum())
    del mk
    gc.collect()
    res['seconds_total'] = time.time() - t0
    res['passes'] = bool(res['inside_bracket'] and res['built_on_count_diffs'] == 0 and res['built_off_count_diffs'] == 0
                         and row['off_check_count_diffs'] == 0 and res['built_flag'] and res['form_after'] == 'sub'
                         and res['weights_reset'] and res['ext_inh_after'] == 0 and not res['save_inh_sims_after'])
    res['row'] = {k: (None if isinstance(x, float) and not np.isfinite(x) else x) for k, x in row.items()}
    _write_json(out / 'smoke.json', res)
    print(f"E92 smoke: {json.dumps({k: x for k, x in res.items() if k not in ('row', 'bisection')}, default=str)}", flush=True)
    return res


# ---------------------------------------------------------------------------
# the references and the runs table

def reference_rows_e91_refits() -> pd.DataFrame:
    """The reference if E91's refitted family passes as a population of eight (its outcome (b)): that family at its fits,
    and the number to beat its population's. olfsysm refuses the divisive form together with `claw_inh_pool`, so B1 then
    compares the gain, run from the refits' parameters with the pooled share 0, against the summed subtraction at the
    refits; parts A2 and A4 are repeated on the refits first (part F). Not written: E91's runs are not read before its
    note is posted."""
    raise NotImplementedError("the E91-refits reference is written only after E91's note is posted and its refits pass")


def reference_rows(family: str = 'C2') -> pd.DataFrame:
    """The references (source 'E89' or 'saved'): E89's C1 rows at REF_WEIGHTINGS as B_<weighting> (S's wirings), and the
    unchanged fits as B_F_k2 on every wiring E89 lacks (T's; their realism computed here)."""
    if family != 'C2':
        if family == 'E91_refits':
            return reference_rows_e91_refits()
        raise ValueError(f'unknown reference family {family!r}')
    key = ('reference_rows', family)
    if key in _CACHE:
        return _CACHE[key].copy()
    rows = []
    p = H.results_dir(E89.EXPERIMENT) / 'runs_C1.csv'
    if p.exists():
        d = pd.read_csv(p)
        d = d[d.permute_seed.isna() & d.weighting.isin([E89.wtag(w) for w in REF_WEIGHTINGS])].drop_duplicates('tag', keep='last')
        for r in d.to_dict('records'):
            r.update(part='E89_C1', arm=ARM, variant=variant_name(E89._wt_of(f"C1_{r['weighting']}")), source='E89',
                     screen='ref', form='sub', form_of='sub', apl_inh_divisive=False, bouton_inh_divisive=False, control_of='')
            rows.append(r)
    have = {(r['variant'], r['wiring']) for r in rows}
    for w in E83.ALL_WIRINGS:
        if (ref_variant(), w) in have:
            continue
        f = E83.saved_fit(w)
        rows.append(dict(part='saved', arm=ARM, tag=f'{f["folder"]}:{f["tag"]}', variant=ref_variant(), source='saved',
                         screen='ref', weighting=E89.wtag(REFERENCE), measure='F', exponent=2.0, form='sub', form_of='sub',
                         apl_inh_divisive=False, bouton_inh_divisive=False, control_of='', wiring=w, permute_seed=np.nan,
                         fit_folder=f['folder'], fit_tag=f['tag'], apl_factor=1.0, at_bound=False,
                         **{k: f['p'][k] for k in E83.PARAM_KEYS}, **E83.realism(f['saved'])))
    df = pd.DataFrame(rows)
    _CACHE[key] = df
    return df.copy()


RUNS_FILES = ('runs_B1.csv', 'runs_B2.csv', 'runs_B3.csv', 'runs_C.csv', 'runs_Dstart.csv')


def runs_table(out: Path, family: str = 'C2') -> pd.DataFrame:
    """Every run of the readings: this experiment's (RUNS_FILES) and the references, with 'variant', the parsed setting
    (screen, form, form_of, control, exponent), 'flagged' and 'at_bound'."""
    parts = [pd.read_csv(out / n) for n in RUNS_FILES if (out / n).exists()]
    parts.append(reference_rows(family))
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    if not len(df):
        return df
    df = df.drop_duplicates(['variant', 'wiring'], keep='first').reset_index(drop=True)
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('screen', 'form', 'form_of', 'control', 'exponent'):
        df[c] = pv[c].to_numpy()
    df['permute_seed'] = np.nan
    df['flagged'] = df.realism_flag.fillna('').astype(str).str.strip().replace('nan', '') != ''
    df['at_bound'] = df.at_bound.fillna(False).astype(bool)
    return df


# ---------------------------------------------------------------------------
# the readings

def _canonical_target(t: str, de: pd.Series) -> pd.Series:
    key = ('target', t, tuple(de.index), de.to_numpy(dtype=float).tobytes())
    if key not in _CACHE:
        _CACHE[key] = H.canonical_pairs(de)
    return _CACHE[key]


def spread_row(d: Dict[str, pd.Series], tg: Dict[str, pd.Series]) -> dict:
    """Per target t, the pattern part's two parts and their companions (`H.apl_effect_spread_shape` and
    `H.apl_effect_stats` of the distances `d` against t's per-pair effect): '<t>_amplitude', '<t>_shape', '<t>_sd_ratio'
    (the model's spread over the pairs as a multiple of the flies'), '<t>_pearson' (Pearson r on the per-pair APL
    effect); pattern^2 = amplitude^2 + shape^2. The harness pairs the data in its own canonical orientation
    (`H.canonical_pairs`), E83's `_cs` in another, and pyolfsysm's align_pairs drops pairs whose orientations differ
    without a word, so the target and the model's distances go in through `H.canonical_pairs` and every pair of the
    target must be aligned (raises otherwise)."""
    row = {}
    on, off = H.canonical_pairs(d['on']), H.canonical_pairs(d['off'])
    for t, de in tg.items():
        data = _canonical_target(t, de)
        ss = H.apl_effect_spread_shape(on, off, data)
        st = H.apl_effect_stats(on, off, data)
        if int(st['apl_effect_n_pairs']) != len(de) or len(de) != len(d['on']):
            raise RuntimeError(f'target {t}: {st["apl_effect_n_pairs"]} pairs aligned of {len(de)} (model {len(d["on"])})')
        row.update({f'{t}_amplitude': ss['apl_effect_amplitude'], f'{t}_shape': ss['apl_effect_shape'],
                    f'{t}_sd_ratio': ss['apl_effect_sd_ratio'], f'{t}_pearson': st['apl_effect_pearson'],
                    f'{t}_pattern_h': st['apl_effect_rmse_centered']})
    return row


def score_row92(d: Dict[str, pd.Series]) -> dict:
    """E89's `_score_row` (E83's `effect_scores` against every target of `targets89`, the block of the per-pair effect)
    with `spread_row` beside; raises unless, for every target, the harness's pattern part is E83's and pattern^2 =
    amplitude^2 + shape^2 (a check that no pair was misaligned)."""
    row = E89._score_row(d)
    sp = spread_row(d, E89.targets89())
    for t in E89.targets89():
        pat = row[f'{t}_pattern']
        if not (np.isclose(sp.pop(f'{t}_pattern_h'), pat, rtol=1e-9, atol=1e-12)
                and np.isclose(pat ** 2, sp[f'{t}_amplitude'] ** 2 + sp[f'{t}_shape'] ** 2, rtol=1e-9, atol=1e-12)):
            raise RuntimeError(f'{t}: the amplitude and shape parts do not add up to the pattern part (pairs misaligned?)')
    row.update(sp)
    return row


def _distances(source: str, tag: str, wiring: str) -> Dict[str, Dict[str, pd.Series]]:
    key = ('dist', source, tag, wiring)
    if key not in _CACHE:
        o = load_any(source, tag, wiring)
        K = E83.saved_fit(wiring)['K']
        _CACHE[key] = {'cal': E83.pair_distances(o, 'cal'), 'calK': E83.pair_distances(o, 'calK', K),
                       'raw': E83.pair_distances(o, 'raw')}
    return _CACHE[key]


def readings92(out: Path, runs: Optional[pd.DataFrame] = None, family: str = 'C2') -> pd.DataFrame:
    """readings_B.csv: per variant and unit (each wiring; S and T where all four wirings have a run) and readout ('cal'
    the calibrated readout with K set on each run's APL-on counts, the readings and the gate; 'calK' the fit's K; 'raw'
    raw counts) `score_row92` against every target of E89's `targets89`; per unit the number of runs outside the bands
    and at a bracket's bound; the setting's fields."""
    runs = runs_table(out, family) if runs is None else runs
    rows = []
    for variant, g in runs.groupby('variant', sort=False):
        ds = {r.wiring: _distances(r.source, r.tag, r.wiring) for r in g.itertuples()}
        meta = {r.wiring: dict(flagged=bool(r.flagged), at_bound=bool(r.at_bound)) for r in g.itertuples()}
        pv = parse_variant(variant)
        units = [(w, (w,)) for w in ds] + [(p, ws) for p, ws in (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS))
                                           if all(w in ds for w in ws)]
        for unit, ws in units:
            for ro in ('cal', 'calK', 'raw'):
                d = ds[unit][ro] if len(ws) == 1 else E83.population([ds[w][ro] for w in ws])
                rows.append(dict(arm=ARM, variant=variant, screen=pv['screen'], weighting=pv['wtag'], exponent=pv['exponent'],
                                 form=pv['form'], form_of=pv['form_of'], control=pv['control'], permute_seed=np.nan,
                                 unit=unit, n_wirings=len(ws), readout=ro,
                                 n_flagged=int(sum(meta[w]['flagged'] for w in ws)),
                                 n_at_bound=int(sum(meta[w]['at_bound'] for w in ws)), **score_row92(d)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'readings_B.csv', index=False)
    return df


def _pat(rd: pd.DataFrame, variant: str, unit: str, target: str, col: str = 'pattern', readout: str = 'cal') -> float:
    return E89._pat(rd, variant, unit, target, col, readout)


def _part_a():
    """Part A's module if it has the names this module reads (PART_A_NAMES), else None (logged)."""
    try:
        from experiments import E92_steps_like_for_like as A
    except Exception as e:                       # (part A is written by another session; not fatal here)
        print(f'E92 B: part A\'s module does not import ({e!r}); its readings are left out', flush=True)
        return None
    missing = [n for n in PART_A_NAMES if not hasattr(A, n)]
    if missing:
        print(f'E92 B: part A\'s module lacks {missing}; its readings are left out', flush=True)
        return None
    return A


def _s_line(ln: dict) -> dict:
    """The line of S's population from `line(readout)` ({'S': ..., 'real': ..., ...}, or one line)."""
    return ln['S'] if isinstance(ln.get('S'), dict) else ln


def beside_b(out: Path, runs: Optional[pd.DataFrame] = None, rd: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """beside_B.csv: per S setting (all four S runs) the two groups' errors (readings_B.csv, all recordings and both
    halves, calibrated readout and raw counts), and, where part A has them, the ester and the alcohol groups' distances
    with APL on in A2's readings (`somata_readings`; the mean over the four wirings) and the place against A4's line of
    S's population on raw counts and under the readout (`relative_distances` of the population's APL-on distances,
    `place_on_line` against `line(readout)['S']`). Part A is another module: any of its functions that fails (not
    implemented yet, or changed) leaves its columns out, with a note in beside_B_notes.json; it never stops the
    scoring."""
    runs = runs_table(out) if runs is None else runs
    rd = pd.read_csv(out / 'readings_B.csv') if rd is None else rd
    A = _part_a()
    notes: Dict[str, str] = {} if A is not None else {'part_a': 'part A\'s module or its names missing'}
    lines = {}
    if A is not None:
        for ro in ('raw', 'cal'):
            try:
                lines[ro] = _s_line(A.line(ro))
            except Exception as e:
                notes[f'line_{ro}'] = f'part A4 line({ro!r}): {e!r}'
    rows = []
    S = list(E83.S_WIRINGS)
    for variant, g in runs[runs.wiring.isin(S)].groupby('variant', sort=False):
        if set(g.wiring) != set(S):
            continue
        pv = parse_variant(variant)
        row = dict(variant=variant, screen=pv['screen'], form=pv['form'], form_of=pv['form_of'], exponent=pv['exponent'])
        for ro in ('cal', 'raw'):
            for t in ('z8', 'halfA', 'halfB'):
                for grp in ('ester', 'alcohol'):
                    row[f'err_{grp}_{t}_{ro}'] = _pat(rd, variant, 'S', t, f'err_{grp}', ro)
        if A is not None and 'somata' not in notes:
            try:
                acc: Dict[str, List[float]] = {}
                for r in g.itertuples():
                    o = load_any(r.source, r.tag, r.wiring)
                    sr = A.somata_readings(o['on'][1], o['off'][1])
                    for q in sr[sr['cond'] == 'on'].to_dict('records'):
                        for grp in ('ester', 'alcohol'):
                            acc.setdefault(f'a2_{q["reading"]}_{grp}_on', []).append(float(q[grp]))
                row.update({k: float(np.mean(x)) for k, x in acc.items() if len(x) == len(S)})
            except Exception as e:
                notes['somata'] = f'part A2 somata_readings: {e!r}'
        for ro, ln in lines.items():
            if f'place_{ro}' in notes:
                continue
            try:
                d_on = E83.population([_distances(r.source, r.tag, r.wiring)[ro] for r in g.itertuples()])['on']
                est, alc = A.relative_distances(d_on)
                pl = A.place_on_line(est, alc, ln)
                row.update({f'line_{ro}_ester_rel': est, f'line_{ro}_alcohol_rel': alc})
                row.update({f'line_{ro}_{k}': (json.dumps(x) if isinstance(x, (dict, list)) else x) for k, x in pl.items()})
            except Exception as e:
                notes[f'place_{ro}'] = f'part A4 relative_distances / place_on_line ({ro}): {e!r}'
        rows.append(row)
    for k, n in notes.items():
        print(f'E92 B: {n}; left out of beside_B.csv', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'beside_B.csv', index=False)
    _write_json(out / 'beside_B_notes.json', dict(notes=notes, part_a_names=list(PART_A_NAMES)))
    return df


# ---------------------------------------------------------------------------
# the co-response beside every S setting (E89's coresponse_c, as E91's coresponse91)

def _co_job92(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3."""
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


def gm_known_error() -> Dict[str, float]:
    """The geometric mean's known error on the co-response (E89's estimator check: the largest error of its L factor on
    the breadth-2 population, per level), {level: error}; empty if E89's estimator_choice.json lacks it."""
    p = H.results_dir(E89.EXPERIMENT) / 'estimator_choice.json'
    try:
        c = json.loads(p.read_text())['populations']['b2']['candidates'][CO_CANDIDATE]['worst']
        return {lv: float(c[lv]['L']) for lv in ('L1', 'L3') if lv in c}
    except (OSError, KeyError, ValueError):
        return {}


def coresponse92(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse_B.csv: per S setting (every variant with all four S runs) the changes of L (co-response) and R
    (recruitment) under APL of the four nearest ester pairs, the ester and the alcohol groups on S's population (E88's
    `log_changes`) in readings 1 and 1b and with the geometric mean (CO_CANDIDATE) at L1 and L3 (its known error beside,
    `gm_known_error`), against the flies' (90 % intervals over resampled recordings; E88's B1 rule: 'differs' outside
    the flies' interval and by a factor of B1_FACTOR or more); 'reaches_flies' if it differs in none of the readings.
    Jobs cached in coresponse_B_cache.pkl."""
    import multiprocessing
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    runs = runs_table(out)
    runs = runs[runs.wiring.isin(E83.S_WIRINGS)]
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    cache = out / 'coresponse_B_cache.pkl'
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    jobs = [(r.source, r.tag, r.wiring, E83.saved_fit(r.wiring)['K'], CO_CANDIDATE, rel, gap) for r in runs.itertuples()
            if (r.source, r.tag) not in done]
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            for res, j in zip(ex.map(_co_job92, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'coresponse {j[1]}', flush=True)
    fl = {}
    fv = E89.flies_values({})
    for subset in E88.SUBSETS:
        tnt = E88.flies_subset(subset)
        rf = E88.flies_factors('z', 'r1', tnt)
        F = {c: rf[c].mean() for c in E89.CONDS}
        fl[('r1', subset)] = (E88.log_changes(F, ('L', 'R')), [E88.log_changes(b, ('L', 'R')) for b in E88.flies_boot(rf)])
        fl[('cand', subset)] = (E88.log_changes(E89.flies_cand(fv, CO_CANDIDATE, None, tnt), ('L', 'R')),
                                [E88.log_changes(b, ('L', 'R')) for b in E89.flies_cand_boot(fv, CO_CANDIDATE, None, tnt)])
    known = gm_known_error()
    rows = []
    for variant, g in runs.groupby('variant', sort=False):
        if set(g.wiring) != set(E83.S_WIRINGS):
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
                        rows.append(dict(arm=ARM, variant=variant, weighting=g.weighting.iloc[0], reading=rdn,
                                         candidate=CO_CANDIDATE if rdn.startswith('cand') else '', subset=subset, factor=f,
                                         group=grp, ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])),
                                         flies_lo=float(np.exp(lo)), flies_hi=float(np.exp(hi)), differs=differs,
                                         known_error=known.get(rdn.split('_')[1], np.nan) if rdn.startswith('cand') else np.nan))
    df = pd.DataFrame(rows)
    if len(df):
        keys = ['arm', 'variant', 'subset', 'factor', 'group']
        reach = df.groupby(keys).differs.any().rename('differs_somewhere').reset_index()
        reach['reaches_flies'] = ~reach.differs_somewhere
        df = df.merge(reach[keys + ['reaches_flies']], on=keys)
    df.to_csv(out / 'coresponse_B.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# part C: the choice, the gate, the control, the equal rule

def _settings92(rd: pd.DataFrame, screens: Sequence[str] = CHOICE_SCREENS) -> pd.DataFrame:
    """The settings of the choice on S under the readout (the references and B1's, by default; no controls)."""
    return rd[(rd.unit == 'S') & (rd.readout == 'cal') & rd.screen.isin(list(screens)) & ~rd.control.astype(bool)]


PARTS_COLS = ('pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson', 'rho')


def choose92(rd: pd.DataFrame, half: str, screens: Sequence[str] = CHOICE_SCREENS, ref: Optional[str] = None) -> dict:
    """The choice on S against half `half`: the setting with the lowest pattern part of S's population under the
    calibrated readout among the reference and the settings whose four runs are inside the bands and not at a bracket's
    bound; it goes to T if its form is divisive and it gains CHOICE_MIN_GAIN or more over the reference. The chosen and
    the reference's parts of the pattern part against the half are beside (thresholds stay on the pattern part)."""
    g = _settings92(rd, screens)
    ref = ref or ref_variant()
    ok = g[(g.variant == ref) | ((g.n_flagged == 0) & (g.n_at_bound == 0))]
    if ok.empty or ref not in set(g.variant):
        return dict(arm=ARM, half=half, chosen=None, confirm=False, note='runs missing')
    best = ok.loc[ok[f'half{half}_pattern'].idxmin()]
    rr = g[g.variant == ref].iloc[0]
    gain = float(rr[f'half{half}_pattern']) - float(best[f'half{half}_pattern'])
    confirm = bool(best.form != 'sub' and gain >= CHOICE_MIN_GAIN)
    row = dict(arm=ARM, half=half, chosen=best.variant, form=best.form, weighting=best.weighting, exponent=float(best.exponent),
               pattern=float(best[f'half{half}_pattern']), ref_pattern=float(rr[f'half{half}_pattern']), gain=gain,
               confirm=confirm, n_eligible=int(len(ok)), n_settings=int(len(g)),
               excluded=json.dumps(sorted(set(g.variant) - set(ok.variant))))
    for c in PARTS_COLS[1:]:
        row[f'chosen_{c}'] = float(best[f'half{half}_{c}'])
        row[f'ref_{c}'] = float(rr[f'half{half}_{c}'])
    return row


def gate_direction92(rd: pd.DataFrame, runs: pd.DataFrame, X: str, chosen: Optional[str], ref: str,
                     controls: Sequence[str]) -> dict:
    """One direction of the gate on T (E91's `gate_direction`: chosen on half X, read against the other half), with the
    parts of the pattern part beside: for the chosen setting, the reference and each control on T under the readout,
    against the other half ('other'), all recordings ('z8') and the six typical ('z6'), the pattern part, the amplitude
    part, the shape part, the spread ratio, Pearson r and ρ on the per-pair APL effect (a description)."""
    r = E91.gate_direction(rd, runs, ARM, X, chosen, ref, controls)
    Y = E89.HALF_OTHER[X]
    for lab, v in [('chosen', chosen), ('ref', ref)] + [(f'control{i}', c) for i, c in enumerate(controls)]:
        if v is None:
            continue
        q = rd[(rd.variant == v) & (rd.unit == 'T') & (rd.readout == 'cal')]
        if not len(q):
            continue
        q = q.iloc[0]
        for key, t in (('other', f'half{Y}'), ('z8', 'z8'), ('z6', 'z6')):
            for c in PARTS_COLS:
                r[f'{lab}_{key}_{c}'] = float(q[f'{t}_{c}'])
    return r


def gate92(rows: Sequence[dict]) -> dict:
    """E83's gate on four directions with E91's rules (`E91.gate91`): (i) the mean gain 0.003 or more; (ii) 12 of the
    16 single-wiring cells lower; (iii) every run inside the bands, the raw-count loss not above the reference's by more
    than 0.005, the mean APL effect within 0.01 of the flies'; (iv) the squared pattern error 15 % lower against all
    recordings and 7.5 % against the six typical; the control gaining less than half as much."""
    return E91.gate91(rows)


def equal_rule92(rd: pd.DataFrame, ref: Optional[str] = None, screens: Sequence[str] = ('B1', 'B2')) -> List[dict]:
    """The equal rule (E91's): a divisive setting at an exponent of EQUAL_MAX_EXPONENT or less (no weighting included) is
    an equal of the family's (the reference, F at exponent 2, subtractive) if, on S under the readout, its pattern part is
    within EQUAL_PATTERN of the reference's against all recordings and both halves, its four runs are inside the bands,
    its raw-count loss within EQUAL_RAW and its mean APL effect within EQUAL_MEAN of the reference's; 'better' if its
    pattern part is lower on all three. B2's settings are read for the record (`screen`)."""
    ref = ref or ref_variant()
    g = _settings92(rd, ('ref',) + tuple(screens))
    out = []
    if ref not in set(g.variant):
        return out
    me0 = float(g[g.variant == ref].mean_effect.iloc[0])
    for r in g.itertuples():
        if r.variant == ref or r.form == 'sub' or float(r.exponent) > EQUAL_MAX_EXPONENT:
            continue
        d = {t: _pat(rd, r.variant, 'S', t) - _pat(rd, ref, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = _pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - _pat(rd, ref, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - me0
        eq = bool(all(abs(x) <= EQUAL_PATTERN for x in d.values()) and int(r.n_flagged) == 0 and int(r.n_at_bound) == 0
                  and abs(raw) <= EQUAL_RAW and abs(me) <= EQUAL_MEAN)
        out.append(dict(arm=ARM, variant=r.variant, screen=r.screen, form=r.form, exponent=float(r.exponent),
                        **{f'd_pattern_{t}': x for t, x in d.items()}, d_raw_loss=raw, d_mean_effect=me,
                        n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound), equal=eq,
                        better=bool(all(x < 0 for x in d.values()))))
    return out


def plan_c(rd: pd.DataFrame) -> Dict[str, dict]:
    return {h: choose92(rd, h) for h in HALVES}


def part_c(out: Path, family: str = 'C2') -> dict:
    """C: the choice per half (choice.csv); every divisive choice with something to confirm run on T at each wiring's fit
    with the APL scale re-held, then its control on S and T (runs_C.csv; resumable); then `score`."""
    if family != 'C2':
        reference_rows(family)
    rd = readings92(out, family=family)
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    chosen = sorted({c['chosen'] for c in plan.values() if c['confirm']})
    if not chosen:
        print('E92 C: nothing to confirm', flush=True)
    for v in chosen:
        b_points(out, E83.T_WIRINGS, [v], 'runs_C.csv', 'C')
    for v in chosen:
        for ws in (E83.S_WIRINGS, E83.T_WIRINGS):
            b_points(out, ws, [control_of(v)], 'runs_C.csv', 'C')
    return score(out, coresponse=False, family=family)


# ---------------------------------------------------------------------------
# scoring and the verdict

def score(out: Path, n_workers: int = 4, coresponse: bool = True, family: str = 'C2') -> dict:
    """readings_B.csv, beside_B.csv, coresponse_B.csv (with `coresponse`), choice.csv, gate.csv, equal.csv,
    verdict_B.json, the figure, comparison_table_B.md."""
    runs = runs_table(out, family)
    rd = readings92(out, runs)
    beside_b(out, runs, rd)
    if coresponse:
        coresponse92(out, n_workers)
    v = verdict92(out, rd, runs)
    figure(out)
    comparison_table(out)
    return v


def _divisor_summary(runs: pd.DataFrame) -> Dict[str, dict]:
    rr = runs[runs.source == 'E92']
    out = {}
    for variant, g in rr.groupby('variant', sort=False):
        d = {}
        for c in ('divisor_kc_p50', 'divisor_kc_p90', 'divisor_kc_max', 'divisor_bouton_p50', 'divisor_bouton_max'):
            if c in g and g[c].notna().any():
                d[c] = float(g[c].max())
        if 'divisor_kc_max' in d:
            d['above_divisor_high'] = bool(d['divisor_kc_max'] > DIVISOR_HIGH)
        out[variant] = d
    return out


def verdict92(out: Path, rd: Optional[pd.DataFrame] = None, runs: Optional[pd.DataFrame] = None) -> dict:
    """verdict_B.json: the choices (four halves), the gate against the reference ('passes' / 'gains on S, not through
    the gate' / 'level'), the equals, B2 against the reference and against the choice (for the record), the flagged runs,
    the divisors, D's plan and scores (if fitted); every S setting's numbers."""
    rd = pd.read_csv(out / 'readings_B.csv') if rd is None else rd
    runs = runs_table(out) if runs is None else runs
    plan = plan_c(rd)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    ref = ref_variant()
    dirs = []
    for X in HALVES:
        c = plan[X]
        chosen = c['chosen'] if c['confirm'] else None
        dirs.append(gate_direction92(rd, runs, X, chosen, ref, [control_of(chosen)] if chosen else []))
    pd.DataFrame(dirs).to_csv(out / 'gate.csv', index=False)
    gate = gate92(dirs)
    eqs = equal_rule92(rd)
    pd.DataFrame(eqs).to_csv(out / 'equal.csv', index=False)
    gains_s = any(c['confirm'] for c in plan.values())
    v = dict(thresholds=dict(CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=E91.GATE_GAIN, GATE_SINGLES_SHARE=E91.GATE_SINGLES_SHARE,
                             GATE_RAW_SLACK=E91.GATE_RAW_SLACK, GATE_MEAN_SLACK=E91.GATE_MEAN_SLACK, GATE_SQ_ALL=E91.GATE_SQ_ALL,
                             GATE_SQ_SIX=E91.GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO, EQUAL_PATTERN=EQUAL_PATTERN,
                             EQUAL_RAW=EQUAL_RAW, EQUAL_MEAN=EQUAL_MEAN, APL_BRACKET=APL_BRACKET, DIVISOR_HIGH=DIVISOR_HIGH),
             readout="cal (the calibrated readout, K set on each run's APL-on counts; E83's readings and gate)",
             reference=ref, choice=plan, gate=gate,
             verdict='passes' if gate['passes'] else ('gains on S, not through the gate' if gains_s else 'level'),
             equals=eqs)
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    v['S'] = {r.variant: dict(screen=r.screen, form=r.form, exponent=r.exponent, halfA=r.halfA_pattern, halfB=r.halfB_pattern,
                              halfC=r.halfC_pattern, halfD=r.halfD_pattern, all=r.z8_pattern, loss=r.z8_loss,
                              mean_effect=r.mean_effect, amplitude=r.z8_amplitude, shape=r.z8_shape, sd_ratio=r.z8_sd_ratio,
                              pearson=r.z8_pearson, rho=r.z8_rho, err_ester=r.z8_err_ester, err_alcohol=r.z8_err_alcohol,
                              n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound)) for r in g.itertuples()}
    # B2 for the record: against the reference and against the best of the choice
    b2 = {}
    for vv in b2_variants():
        if vv in v['S']:
            b2[vv] = {f'd_pattern_{t}': _pat(rd, vv, 'S', t) - _pat(rd, ref, 'S', t) for t in ('halfA', 'halfB', 'z8')}
            b2[vv]['lower_than_every_choice'] = bool(all(
                plan[h]['chosen'] is None or _pat(rd, vv, 'S', f'half{h}') < plan[h]['pattern'] for h in HALVES))
    v['B2'] = b2
    rr = runs[runs.flagged & (runs.source == 'E92')]
    v['flagged_runs'] = rr[['tag', 'realism_flag']].values.tolist()
    v['divisors'] = _divisor_summary(runs)
    v['at_bound_runs'] = runs[runs.at_bound & (runs.source == 'E92')].tag.tolist()
    try:
        v['D_plan'] = d_plan(out, rd)
    except Exception as e:                      # (D's plan needs the shared fits; reported, not fatal)
        v['D_plan'] = dict(error=repr(e))
    if (out / 'd_scores.json').exists():
        v['D'] = json.loads((out / 'd_scores.json').read_text())
    if (out / 'coresponse_B.csv').exists():
        co = pd.read_csv(out / 'coresponse_B.csv')
        if len(co):
            reach = co[(co.factor == 'L') & (co.group == 'alcohol')].drop_duplicates(['variant', 'subset'])
            v['alcohol_coresponse_reaches_flies'] = {f'{r.variant}:{r.subset}': bool(r.reaches_flies) for r in reach.itertuples()}
    _write_json(out / 'verdict_B.json', v)
    print(f"E92 C: {v['verdict']}; " + '; '.join(f"{h} {plan[h].get('chosen')} (gain {plan[h].get('gain', np.nan):.4f})"
                                                for h in HALVES), flush=True)
    return v


# ---------------------------------------------------------------------------
# part D: the family fitted again

def d_plan(out: Path, rd: Optional[pd.DataFrame] = None) -> dict:
    """Whether D runs, and at which exponents: B1's choices that go to T (`choose92`'s confirm) and B1's equals (equal or
    better, `equal_rule92`), each at its setting's exponent."""
    rd = pd.read_csv(out / 'readings_B.csv') if rd is None else rd
    why: Dict[float, List[str]] = {}
    for h, c in plan_c(rd).items():
        if c['confirm'] and parse_variant(c['chosen'])['screen'] == 'B1':
            why.setdefault(float(c['exponent']), []).append(f"choice on half {h}: {c['chosen']}")
    for e in equal_rule92(rd, screens=('B1',)):
        if e['equal'] or e['better']:
            why.setdefault(float(e['exponent']), []).append(f"{'equal' if e['equal'] else 'better'}: {e['variant']}")
    return dict(allowed=bool(why), exponents=sorted(why, reverse=True), why={f'{k:g}': x for k, x in why.items()})


def d_tag(label: str, k: float) -> str:
    """The tag E10_joint_fit.py gives a D fit (its prefix rules: the wiring, rate0.068_, aplboth_, div_, then
    k<k>_cap3)."""
    from experiments import E10_joint_fit as J
    a = E91.d_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.apl_divisive_tag(True)}k{k:g}_cap3"


def d_start_key(label: str, k: float) -> str:
    return f'{label}_k{k:g}'


D_START_FILE = 'd_start_params.json'


def d_start_params(out: Path, label: str, k: float) -> dict:
    """D's start on a wiring at exponent k: its C2 fit's absolute parameters with the tied APL scales (E83's TIED_ON)
    times the factor of B1's re-held point of the divisive form at that exponent on that wiring (runs_B1 / runs_C /
    runs_Dstart; run into runs_Dstart.csv first if no file has it), written into results/.../d_start_params.json under
    `d_start_key` (E10_joint_fit.py --init-params)."""
    v = variant_name(('F', float(k)), 'div')
    tag = run_tag(v, label)
    rows = [pd.read_csv(out / n) for n in RUNS_FILES if (out / n).exists()]
    have = pd.concat(rows, ignore_index=True, sort=False) if rows else pd.DataFrame(columns=['tag'])
    have = have[have.tag == tag]
    if not len(have):
        have = b_points(out, [label], [v], 'runs_Dstart.csv', 'Dstart')
        have = have[have.tag == tag]
    r = have.iloc[0]
    if bool(r.at_bound):
        raise RuntimeError(f'{tag}: the re-held point is at a bound of the APL bracket; no start for D')
    p = E83.saved_fit(label)['p']
    x = float(r.apl_factor)
    start = dict(p, **{kk: p[kk] * x for kk in E83.TIED_ON})
    path = out / D_START_FILE
    d = json.loads(path.read_text()) if path.exists() else {}
    d[d_start_key(label, k)] = dict(start, apl_factor=x, from_tag=tag)
    _write_json(path, d)
    return d[d_start_key(label, k)]


def d_command(label: str, k: float, start: str = 'reheld', maxfev: int = D_MAXFEV) -> List[str]:
    """The E10_joint_fit.py command of one D fit: C2's family (`-k <k> --aplpn-max 3 --rate-target ctrl_frac --loss apl
    --apl-readout both --tag _cap3`) with `--apl-divisive`, into this experiment's folders, from the wiring's own fit:
    'reheld' (the default) at its parameters with the tied APL scales at B1's re-held factor (`--init-params` on
    d_start_params.json, which `d_start_params` writes), 'exact' at its own scales (`--init-from`, exactly at the same
    exponent, as E91's D)."""
    if start not in D_STARTS:
        raise ValueError(f'start {start!r}: one of {D_STARTS}')
    m = E91.d_member(label)
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', f'{k:g}', '--aplpn-max', '3',
           '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', '_cap3', '--maxfev', str(int(maxfev)),
           '--apl-divisive', '--experiment', EXPERIMENT] + m['wiring_args']
    if start == 'reheld':
        cmd += ['--init-params', f'results/{EXPERIMENT}/{D_START_FILE}:{d_start_key(label, k)}']
    else:
        cmd += ['--init-from', m['tag'], '--init-from-experiment', m['folder']]
        if float(k) == m['k']:
            cmd.append('--init-from-exact')
    return cmd


def _fits(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def d_allowed(out: Path, pop: str, k: float) -> Tuple[bool, str]:
    """The entry's rules for D: S's fits at exponent k only if `d_plan` names k; T's only if S's fitted population at k
    is D_S_GAIN or more below the reference's pattern part against all recordings under the readout."""
    if pop == 'S':
        if not (out / 'readings_B.csv').exists():
            return False, 'no readings (run --part score after B1)'
        pl = d_plan(out)
        ok = float(k) in pl['exponents']
        return ok, f"exponent {k:g}: {pl['why'].get(f'{k:g}', 'not named by a choice or an equal')}"
    sd = score_d(out, write=False)
    p = sd.get(f'k{k:g}', {}).get('S', {}).get('pattern_all')
    if p is None:
        return False, f"S's fits at exponent {k:g} not scored"
    ok = bool(p <= UNCHANGED_S['all'] - D_S_GAIN)
    return ok, f"S's population pattern part {p:.4f} against {UNCHANGED_S['all']} - {D_S_GAIN}"


def part_d(out: Path, pop: str = 'S', k: Optional[float] = None, start: str = 'reheld', force: bool = False,
           dry_run: bool = False) -> List[str]:
    """Runs D's fits of population `pop` at exponent `k` (default: every exponent `d_plan` names) one after another
    (`d_command`; a process each, its log under outputs/.../D/), skipping those in this experiment's fits.csv; refuses
    what the entry's rules do not allow (`d_allowed`) unless `force`. Returns the tags."""
    ks = [float(k)] if k is not None else d_plan(out)['exponents']
    if not ks:
        print('E92 D: nothing to fit (no choice by C and no equal)', flush=True)
        if not force:
            return []
    if not dry_run:
        _commits()
    logs = H.outputs_dir(EXPERIMENT) / 'D'
    logs.mkdir(parents=True, exist_ok=True)
    tags = []
    for kk in ks:
        ok, why = d_allowed(out, pop, kk)
        print(f'E92 D {pop} k{kk:g}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
        if not ok and not force:
            continue
        for w in _wirings(pop):
            tag = d_tag(w, kk)
            tags.append(tag)
            if tag in set(_fits(out).tag):
                print(f'E92 D {tag}: done', flush=True)
                continue
            if start == 'reheld' and not dry_run:
                d_start_params(out, w, kk)
            cmd = d_command(w, kk, start)
            print('E92 D:', ' '.join(cmd), flush=True)
            if dry_run:
                continue
            with open(logs / f'{tag}.log', 'w') as fh:
                r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
            if r.returncode:
                raise RuntimeError(f'E92 D {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
            if tag not in set(_fits(out).tag):
                raise RuntimeError(f'E92 D: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def _d_row(sc: dict, sc0: dict, raw: dict, raw0: dict) -> dict:
    return dict(pattern_all=sc['z8_pattern'], pattern_all_unchanged=sc0['z8_pattern'], halfA=sc['halfA_pattern'],
                halfB=sc['halfB_pattern'], loss=sc['z8_loss'], mean_effect=sc['mean_effect'], amplitude=sc['z8_amplitude'],
                shape=sc['z8_shape'], sd_ratio=sc['z8_sd_ratio'], pearson=sc['z8_pearson'], rho=sc['z8_rho'],
                raw_loss=raw['z8_loss'], raw_loss_unchanged=raw0['z8_loss'], err_ester=sc['z8_err_ester'],
                err_alcohol=sc['z8_err_alcohol'])


def score_d(out: Path, write: bool = True) -> dict:
    """D's scores: per exponent fitted, wiring and population (S, T and the eight) the readings under the calibrated
    readout and raw counts (`score_row92` against `targets89`), against the unchanged fits of the same wirings, with the
    realism; the family as a population of eight passes by E91's part D rule: its pattern part D_PASS_GAIN or more below
    UNCHANGED_8, at least D_WIRINGS_LOWER of the eight lower, every fit inside the bands, the raw-count loss not above
    D_RAW_REF by more than D_RAW_SLACK, the mean APL effect within D_MEAN_SLACK of the flies'. d_scores.csv,
    d_scores.json."""
    from experiments import E72_wiring_population as E72
    fits = _fits(out)
    res, rows = {}, []
    for k in D_EXPONENTS:
        st = {}
        for w in E83.ALL_WIRINGS:
            try:
                tag = d_tag(w, k)
            except Exception:
                continue
            if tag not in set(fits.tag):
                continue
            o = E72.load_fit(EXPERIMENT, tag)
            if o is None:
                continue
            d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
            d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
            sc, sc0 = score_row92(d['cal']), score_row92(d0['cal'])
            raw, raw0 = E83.effect_scores(d['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
            fr = fits[fits.tag == tag].iloc[-1]
            row = dict(k=k, unit=w, tag=tag, **_d_row(sc, sc0, raw, raw0), lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                       **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                                 'kc_kc_inh', 'apl_inh_divisive', 'init_params', 'n_evals')},
                       fit_loss=fr.get('loss'),
                       **E83.realism(o))
            rows.append(row)
            st[w] = (d, d0, row)
        if not st:
            continue
        rs = {}
        for unit, ws in (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS), ('all8', E83.ALL_WIRINGS)):
            if not all(w in st for w in ws):
                continue
            dp = {ro: E83.population([st[w][0][ro] for w in ws]) for ro in ('cal', 'raw')}
            d0 = {ro: E83.population([st[w][1][ro] for w in ws]) for ro in ('cal', 'raw')}
            sc, sc0 = score_row92(dp['cal']), score_row92(d0['cal'])
            raw, raw0 = E83.effect_scores(dp['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
            rs[unit] = dict(_d_row(sc, sc0, raw, raw0), n_lower=int(sum(st[w][2]['lower'] for w in ws)),
                            n_flagged=int(sum(str(st[w][2].get('realism_flag') or '').strip() not in ('', 'nan') for w in ws)))
            rows.append(dict(k=k, unit=unit, **rs[unit]))
        if 'all8' in rs:
            e = rs['all8']
            e['passes'] = bool(e['pattern_all'] <= UNCHANGED_8 - D_PASS_GAIN and e['n_lower'] >= D_WIRINGS_LOWER
                               and e['n_flagged'] == 0 and e['raw_loss'] <= D_RAW_REF + D_RAW_SLACK
                               and abs(e['mean_effect'] - FLIES_MEAN) <= D_MEAN_SLACK)
        res[f'k{k:g}'] = rs
    if write:
        pd.DataFrame(rows).to_csv(out / 'd_scores.csv', index=False)
        _write_json(out / 'd_scores.json', res)
    return res


# ---------------------------------------------------------------------------
# parts E and F (stubs)

def part_e(out: Path) -> None:
    """Part E (only on a pass by C or D), for the passing setting: the eight wirings again with the stage matrices saved
    (E90's part B) for the two steps like for like and E90's step target; the firing-ceiling family with the response
    rate with APL on re-held; S at the flies' 2 s pulse with both response rates re-held; S through the trial noise as
    E72's part E ran it; each holds if the population keeps at least half of its gain. Not written: it runs only on a
    pass."""
    raise NotImplementedError('E92 part E runs only on a pass by part C or part D; not written')


# ---------------------------------------------------------------------------
# part F: the summed form (E91, posted 2026-11-13 19:38) read as part A read C2's family

F_SOURCE = E91.EXPERIMENT
F_K, F_POOL = 2.0, E91.D_POOL                       # E91's part D (stage S2): C2's family fitted again at s = 1, exponent 2
F_B1 = E91.variant_name('B1', REFERENCE, 1.0)       # 'B1_F_k2_s1': s = 1 on the fits as they are, the APL scale re-held
F_B2_REF = E91.ref_variant('B2')                    # 'B2_F_k2': the base with the flies' recruitment, its own family (s = 0)
F_B2_CHOICE = E91.variant_name('B2', ('F', 1.0), 1.0)   # 'B2_F_k1_s1': E91's choice on that base (exponent 1, s = 1)
F_FAMILIES = {
    'C2': "C2's clipped family, the unchanged fits (E90's runs of them)",
    'E91_B1_k2_s1': 'the summed form (s = 1) on the fits as they are at exponent 2, the APL scale re-held (E91 B1)',
    'E91_refit': "C2's family fitted again with the summed form (s = 1, exponent 2; E91's part D)",
    'E91_B2_k2': "the base with the flies' recruitment at exponent 2, clipped (E89's C2 runs on S, E91's on T)",
    'E91_B2_k1_s1': "E91's choice on the base with the flies' recruitment: exponent 1, s = 1 (E91 B2 on S, C on T)",
}
F_STAGE_FAMILIES = ('E91_refit', 'E91_B2_k1_s1')    # run again on S with the stage matrices saved (E90's reading)
F_SAME_SIDE_4 = 3                                   # a population of four: 3 of its 4 wirings on the side (part A's D1)
F_NOISE_SEED = 300                                  # E90's NOISE_SEED + this + the run's index, for the claws' noisy reading
F_DIR = 'F'
STAGE_KEYS = ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')
# the line readings of part F beyond E91's own folder: E89's C runs of the base with the flies' recruitment on S
F_LINE_E89 = ('C2_F_k2', 'C2_F_k1')
# E91's variants as they appear in line_F.csv (folder, family), per F family
F_LINE_OF = {'C2': ('E90_model_stages', 'C2_unchanged'), 'E91_B1_k2_s1': (F_SOURCE, F_B1),
             'E91_refit': (F_SOURCE, 'E91_refit'), 'E91_B2_k2': ((E89.EXPERIMENT, 'C2_F_k2'), (F_SOURCE, F_B2_REF)),
             'E91_B2_k1_s1': (F_SOURCE, F_B2_CHOICE)}


def _A():
    from experiments import E92_steps_like_for_like as A
    return A


def _f_dir() -> Path:
    d = H.outputs_dir(EXPERIMENT) / F_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def f_refit_tags() -> Dict[str, str]:
    """{E91's part D fit tag: S wiring} (E91's `d_tag` at exponent F_K and s = F_POOL)."""
    return {E91.d_tag(w, F_K, F_POOL): w for w in E83.S_WIRINGS}


def f_runs() -> pd.DataFrame:
    """The runs part F reads: family (F_FAMILIES), wiring, source ('E90' E90's part B run of C2's fit, 'E91' E91's B / C
    runs, 'E91fit' E91's part D fits, 'E89' E89's C runs), tag."""
    rows = [dict(family='C2', wiring=w, source='E90', tag=w) for w in E83.ALL_WIRINGS]
    for w in E83.S_WIRINGS:
        rows.append(dict(family='E91_B1_k2_s1', wiring=w, source='E91', tag=E91.run_tag(F_B1, w)))
    rows += [dict(family='E91_refit', wiring=w, source='E91fit', tag=t) for t, w in f_refit_tags().items()]
    rows += [dict(family='E91_B2_k2', wiring=w, source='E89', tag=E89.run_tag('C2', REFERENCE, w)) for w in E83.S_WIRINGS]
    rows += [dict(family='E91_B2_k2', wiring=w, source='E91', tag=E91.run_tag(F_B2_REF, w)) for w in E83.T_WIRINGS]
    rows += [dict(family='E91_B2_k1_s1', wiring=w, source='E91', tag=E91.run_tag(F_B2_CHOICE, w)) for w in E83.ALL_WIRINGS]
    return pd.DataFrame(rows)


def f_load(source: str, tag: str) -> Dict[str, tuple]:
    """One of part F's runs ({'on', 'off'}: (responses, spike counts)), from the local outputs or the archive."""
    from experiments import E72_wiring_population as E72
    if source == 'E90':
        return _A().e90_run('B', tag)
    o = {'E91': E91.load_run, 'E89': E89.load_run, 'E91fit': lambda t: E72.load_fit(F_SOURCE, t)}[source](tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


# --- F1: the somata in part A's A2 readings

def _f_soma_job(job: tuple) -> pd.DataFrame:
    """One run's somata in every reading of A2 (part A's `somata_pairs`: the 66 megamat pairs)."""
    family, wiring, source, tag = job
    warnings.simplefilter('ignore')
    o = f_load(source, tag)
    return _A().somata_pairs(o['on'][1], o['off'][1]).assign(family=family, wiring=wiring, source=source, tag=tag)


F_POPULATIONS = (('S', E83.S_WIRINGS), ('T', E83.T_WIRINGS), ('all8', E83.ALL_WIRINGS))


def f_soma_means(pairs: pd.DataFrame) -> pd.DataFrame:
    """Group means (part A's GROUPS) of A2's pair rows per run ('run', wiring) and per population of each family
    (F_POPULATIONS where all its wirings are present: the mean over its runs of the pair distances, as part A's
    `somata_model_table`)."""
    A = _A()
    runs = A._groups_of_pairs(pairs, ['family', 'wiring', 'cond', 'reading'])
    runs.insert(0, 'kind', 'run')
    out = [runs]
    for fam, g in pairs.groupby('family', sort=False):
        for name, ws in F_POPULATIONS:
            sub = g[g.wiring.isin(ws)]
            if set(sub.wiring) != set(ws):
                continue
            p = sub.groupby(['cond', 'reading', 'odor_a', 'odor_b'], sort=False)['value'].mean().reset_index()
            gm = A._groups_of_pairs(p, ['cond', 'reading'])
            gm.insert(0, 'kind', 'population')
            gm.insert(1, 'family', fam)
            gm.insert(2, 'wiring', name)
            out.append(gm)
    return pd.concat(out, ignore_index=True)


def _n_same(n_wirings: int) -> int:
    """The wirings required on the population's side: E90's 6 of 8; 3 of 4 for a population of four."""
    if n_wirings == len(E83.ALL_WIRINGS):
        return _A().VERDICT_SAME_SIDE
    if n_wirings == len(E83.S_WIRINGS):
        return F_SAME_SIDE_4
    raise ValueError(f'a population of {n_wirings} wirings')


def f_soma_rows(sm: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool], flies=None) -> pd.DataFrame:
    """Every (family, population, group, condition, reading) comparison of the somata with the flies' (part A's
    `compare`, `_n_same` wirings on the side): the like-for-like readings (set 'like_for_like'; a noise-removed
    comparison of a group that failed A1's check is kept with read False) and, beside, raw counts and the calibrated
    readout against the flies in baseline-noise units (set 'beside'). Populations: all8 where a family has eight
    wirings, and S."""
    A = _A()
    flies = A.flies_group if flies is None else flies
    comps = [(p, u, fr, ms, 'like_for_like') for p, u, fr, ms in A.like_for_like()]
    comps += [(ro, 'noise', fr, ms, 'beside') for ro in ('raw', 'cal') for fr in A.FLY_READINGS for ms in A.MEASURES]
    rows = []
    pops = sm[sm.kind == 'population']
    for fam in dict.fromkeys(sm.family):
        for name, ws in F_POPULATIONS:
            if name == 'T' or not len(pops[(pops.family == fam) & (pops.wiring == name)]):
                continue
            for prefix, u, fr, ms, cset in comps:
                reading = f'{prefix}_{ms}'
                for cond in A.CONDS:
                    x = sm[(sm.family == fam) & (sm.reading == reading) & (sm.cond == cond)]
                    p = x[(x.kind == 'population') & (x.wiring == name)]
                    w = x[(x.kind == 'run') & x.wiring.isin(ws)]
                    if len(p) != 1 or len(w) != len(ws):
                        raise RuntimeError(f'{fam} {name} {reading} {cond}: {len(p)} population rows, {len(w)} runs')
                    for r in A.compare(p.iloc[0], w, flies(u, A.GENO[cond], fr, ms), n_same=_n_same(len(ws))):
                        read = bool(fr != 'noise_removed' or passes.get((r['group'], ms, u), False))
                        rows.append(dict(family=fam, population=name, element='soma', set=cset, cond=cond, reading=reading,
                                         flies_units=u, flies_reading=fr, measure=ms, read=read, **r))
    return pd.DataFrame(rows)


def part_f1(out: Path, workers: int = 4) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """F1: somata_F.csv (the comparisons, `f_soma_rows`) and somata_F_means.csv (the group means, `f_soma_means`); the
    pair rows under outputs/.../F/somata_F_pairs.parquet."""
    A = _A()
    t0 = time.time()
    runs = f_runs()
    res = A._map(_f_soma_job, [tuple(r) for r in runs[['family', 'wiring', 'source', 'tag']].itertuples(index=False)],
                 workers)
    pairs = pd.concat(res, ignore_index=True)
    pairs.to_parquet(_f_dir() / 'somata_F_pairs.parquet')
    sm = f_soma_means(pairs)
    sm.to_csv(out / 'somata_F_means.csv', index=False)
    rows = f_soma_rows(sm, A.check_passes())
    rows.to_csv(out / 'somata_F.csv', index=False)
    print(f'F1 in {time.time() - t0:.0f} s', flush=True)
    return sm, rows


# --- F2: the line

def _f_line_row(folder: str, tag: str, family: str, wiring: str, files: Dict[str, List[Path]]) -> dict:
    A = _A()
    paths = files.get(tag, [])
    path = off = None
    for f in sorted(paths, key=lambda x: len(x.parts)):
        o = f.with_name(f'{tag}_acute_off_spike_counts.parquet')
        if f.with_name(f'{tag}_on_responses.parquet').exists() and o.exists() and A._resp(o).exists():
            path, off = f, o
            break
    reason = '' if path is not None else 'no saved APL-on counts with an APL-off run'
    return dict(folder=folder, tag=tag, family=family, wiring=wiring, kind='oracle' if 'oracle' in tag else 'variant',
                path=str(path) if path else '', off_path=str(off) if path else '', off_from='own' if path else '',
                n_paths=len(paths), included=not reason, reason=reason)


def f_line_variants() -> pd.DataFrame:
    """The archived runs F2 places against A4's line, in part A's `a4_variants` layout: every tag of E91's params.csv
    (its part D fits as family 'E91_refit' on their S wirings, `f_refit_tags`; the others by part A's `_variant_key`),
    E89's C runs of F_LINE_E89 on S, and E90's eight runs of C2 (the references for T and single runs)."""
    A = _A()
    rows = []
    p = pd.read_csv(A.ARCHIVE / F_SOURCE / 'params.csv', low_memory=False).drop_duplicates('tag')
    files = A._find_files(F_SOURCE)
    refit = f_refit_tags()
    for tag in p.tag:
        fam, w = ('E91_refit', refit[tag]) if tag in refit else A._variant_key(F_SOURCE, tag, None)
        rows.append(_f_line_row(F_SOURCE, tag, fam, w, files))
    f89 = A._find_files(E89.EXPERIMENT)
    for fam in F_LINE_E89:
        for w in E83.S_WIRINGS:
            rows.append(_f_line_row(E89.EXPERIMENT, f'{fam}_{w}', fam, w, f89))
    for w in E83.ALL_WIRINGS:
        d = next(x for x in A._e90_dir('B') if (x / f'{w}_on_spike_counts.parquet').exists())
        rows.append(dict(folder='E90_model_stages', tag=w, family='C2_unchanged', wiring=w, kind='reference',
                         path=str(d / f'{w}_on_spike_counts.parquet'), off_path=str(d / f'{w}_acute_off_spike_counts.parquet'),
                         off_from='own', n_paths=1, included=True, reason=''))
    return pd.DataFrame(rows)


def part_f2(out: Path, workers: int = 4) -> pd.DataFrame:
    """F2: line_F.csv, part A's `variant_rows` of `f_line_variants` against A4's lines (line.json): per family the S and
    T populations (single runs otherwise), their relative distances, place against the line, pattern part with its
    amplitude and shape parts, spread ratio, Pearson r and ρ, the loss, the mean APL effect and the realism."""
    A = _A()
    t0 = time.time()
    v = f_line_variants()
    inc = v[v.included]
    res = {r['idx']: r for r in A._map(A._a4_job, list(zip(inc.index, inc.path, inc.off_path)), workers)}
    errs = {int(k): r['error'] for k, r in res.items() if r['error']}
    if errs:
        raise RuntimeError(f'F2: {errs}')
    lines = {ro: A.line(ro) for ro in A.READOUTS}
    vr = A.variant_rows(v, res, lines)
    vr.to_csv(out / 'line_F.csv', index=False)
    v.to_csv(_f_dir() / 'line_F_variants.csv', index=False)
    print(f'F2: {len(inc)} runs of {len(v)} in {time.time() - t0:.0f} s', flush=True)
    return vr


# --- F3: the stage runs

def read_stages_odor(mp, rv, W: pd.DataFrame, i: int, window: Tuple[float, float]) -> Dict[str, object]:
    """E90's `read_stages` of one odor (olfsysm's odor index i) read from olfsysm's per-odor lists directly, so that it
    holds after a run of that odor alone (`mp.sim_only`), where the other odors' saved variables can be empty: the window
    and pre-odor means of the PNs (glomeruli), the boutons and the claws (`E90.claw_inputs`: wiring x bouton rate minus
    the claw's APL weight times APL's activity at the step before, clipped at 0), and APL's ('apl_win', 'apl_pre')."""
    from experiments import E90_model_stages as E90
    t = E90.step_times(mp, rv)
    w = E90._windows(t, window)
    gl = list(dict.fromkeys(W.columns.get_level_values('glomerulus')))
    pn = np.asarray(rv.pn.pn_sims[i], dtype=float)
    bt = np.asarray(rv.pn.bouton_sims[i], dtype=float)
    a = np.asarray(rv.kc.inh_sims[i], dtype=float).ravel()
    wk = np.asarray(rv.kc.wAPLKC, dtype=float).ravel()
    if pn.shape[0] != len(gl) or bt.shape[0] != W.shape[1] or a.size != t.size or wk.size != W.shape[0]:
        raise RuntimeError(f'odor {i}: PNs {pn.shape}, boutons {bt.shape}, APL {a.shape}, {wk.size} APL>claw weights for '
                           f'the wiring {W.shape} and {len(gl)} glomeruli, {t.size} steps')
    Wm = W.to_numpy(dtype=float)
    out: Dict[str, object] = {}
    for k in ('win', 'pre'):
        idx = w[k]
        out[f'pn_{k}'] = pn[:, idx].mean(axis=1)
        b = bt[:, idx]
        if np.isnan(b).any():
            raise RuntimeError(f'NaN bouton rates in the {k} window')
        out[f'bouton_{k}'] = b.mean(axis=1)
        out[f'claw_{k}'] = E90.claw_inputs(Wm, b, a[idx - 1], wk).mean(axis=1)
    out['apl_win'], out['apl_pre'] = float(a[w['win']].mean()), float(a[w['pre']].mean())
    return out


def stage_frames(per_odor: Sequence[Dict[str, object]], odors: Sequence[str], W: pd.DataFrame) -> dict:
    """E90's stage dict (`read_stages`' layout: '<stage>_<win|pre>' units x odors, 'apl' odors x {win, pre}) from
    `read_stages_odor`'s records, one per odor in the order of `odors` (abbreviated)."""
    gl = list(dict.fromkeys(W.columns.get_level_values('glomerulus')))
    idx = {'pn': pd.Index(gl, name='glomerulus'), 'bouton': W.columns, 'claw': W.index}
    st = {k: pd.DataFrame(np.column_stack([np.asarray(r[k], dtype=float) for r in per_odor]), index=idx[k.split('_')[0]],
                          columns=list(odors)) for k in STAGE_KEYS}
    st['apl'] = pd.DataFrame({o: dict(win=float(r['apl_win']), pre=float(r['apl_pre'])) for o, r in zip(odors, per_odor)}).T
    return st


def _f_b2_worker(conn, label: str, parts: dict, fit: dict, js: List[int]) -> None:
    """A worker of `FPool` (E87's WiringPool protocol): wiring `label`'s simulator on the ORN columns `js` (as E91's
    `_c91_worker`); for each task ({'thr', 'apl': per odor of js; 'wAPLKC': per claw; 'claw_inh_pool'}) each odor alone
    (`mp.sim_only`) at its threshold factor, APL on at its tied APL scale (E80's `_tied`) and APL off, APL's activity
    saved: its responses, spike counts and stage matrices (`read_stages_odor`); the weights, the pooled share and the
    saving reset after each task; until None."""
    try:
        from experiments import E80_screen as S80
        from experiments import E90_model_stages as E90
        E83._PARTS[label] = parts
        E83._FITS[label] = fit
        sim, info = E83.build_simulator(label, orn=parts['orn'].iloc[:, js], breadth=parts['breadth'])
        cells = info['cells']
        mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        names = [str(c) for c in cells['odor_index'].cell_contents]
        kc = cells['wAPLKC_unscaled'].cell_contents
        kc0 = kc.copy()
        n_kc = int(fit['saved']['on'][0].shape[0])
        W, p = info['W'], dict(info['params'])
        conn.send(('ready', names))
    except Exception as e:
        conn.send(('error', repr(e)))
        return
    while True:
        msg = conn.recv()
        if msg is None:
            break
        try:
            wk = msg.get('wAPLKC')
            kc[...] = kc0 if wk is None else np.asarray(wk, dtype=float).reshape(kc0.shape)
            recs = []
            try:
                E91.set_pool(mp, float(msg.get('claw_inh_pool', 0.0)))
                mp.kc.save_inh_sims = True
                for jj in range(len(order)):
                    rec = dict(odor=names[jj])
                    for cond, q in (('on', S80._tied(p, float(msg['apl'][jj]))),
                                    ('off', dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))):
                        mp.sim_only = [int(order[jj])]
                        try:
                            r, s = sim(q, kc_thr_factor=np.full(n_kc, float(msg['thr'][jj])))
                        finally:
                            mp.sim_only = []
                        rec[cond] = dict(read_stages_odor(mp, rv, W, int(order[jj]), E90.MODEL_WINDOW),
                                         responses=r.iloc[:, jj].to_numpy(dtype=float),
                                         counts=s.iloc[:, jj].to_numpy(dtype=float))
                    recs.append(rec)
            finally:
                kc[...] = kc0
                E91.set_pool(mp, 0.0)
                mp.kc.save_inh_sims = False
                mp.sim_only = []
            conn.send(('ok', recs))
        except Exception as e:
            conn.send(('error', repr(e)))


class FPool(E87.WiringPool):
    """E87's WiringPool with `_f_b2_worker`: `run_f(thr, apl, wk, s)` runs every odor alone at its threshold factor
    and tied APL scale (lists in the saved outputs' column order), with the APL->KC weights `wk` (a per-claw Series, or
    None for the fit's) and the pooled share s, and returns the per-odor records."""

    _WORKER = staticmethod(_f_b2_worker)

    def run_f(self, thr: Sequence[float], apl: Sequence[float], wk: Optional[pd.Series], s: float) -> List[dict]:
        arr = None if wk is None else wk.to_numpy(dtype=float)
        for conn, js in zip(self.conns, self.groups):
            conn.send(dict(thr=[float(thr[j]) for j in js], apl=[float(apl[j]) for j in js], wAPLKC=arr,
                           claw_inh_pool=float(s)))
        recs: List[Optional[dict]] = [None] * len(self.cols)
        for conn, js in zip(self.conns, self.groups):
            status, v = conn.recv()
            if status != 'ok':
                raise RuntimeError(f'{self.label}: a worker failed: {v}')
            for j, rec in zip(js, v):
                if rec['odor'] != self.cols[j]:
                    raise RuntimeError(f"odor {j}: {rec['odor']} where {self.cols[j]} was expected")
                recs[j] = rec
        return recs


def b2_odor_params(row: pd.Series, cols: Sequence[str]) -> Tuple[List[float], List[float]]:
    """The per-odor threshold factors and tied APL scales of one of E91's runs on the base with the flies' recruitment,
    as run (its row's thr_per_odor and apl_per_odor), in the order of the odor columns `cols`."""
    thr, apl = json.loads(row['thr_per_odor']), json.loads(row['apl_per_odor'])
    missing = [c for c in cols if c not in thr or c not in apl]
    if missing or len(thr) != len(cols) or len(apl) != len(cols):
        raise RuntimeError(f"{row['tag']}: per-odor parameters for {sorted(thr)}, the run's odors are {list(cols)}")
    return [float(thr[c]) for c in cols], [float(apl[c]) for c in cols]


def e91_run_row(tag: str) -> pd.Series:
    """E91's runs_*.csv row of a run (the last if repeated)."""
    for name in ('runs_B2.csv', 'runs_C.csv', 'runs_B1.csv', 'runs_B3.csv'):
        p = H.results_dir(F_SOURCE) / name
        if p.exists():
            d = pd.read_csv(p)
            d = d[d.tag == tag]
            if len(d):
                return d.iloc[-1]
    raise RuntimeError(f"no row of E91's runs for {tag}")


def f_stage_tag(family: str, wiring: str) -> str:
    return f'F_{family}_{wiring}'


def units_meta(W: pd.DataFrame) -> dict:
    """E90's `unit_meta` of a wiring (not cached): the boutons' and the claws' glomeruli, the claws' KCs."""
    from experiments import E90_model_stages as E90
    return dict(bouton_glomerulus=[str(x) for x in W.columns.get_level_values('glomerulus')],
                claw_glomerulus=[str(x) for x in E90.glomerulus_of_claws(W)],
                kc_of_claw=[int(x) for x in W.index.get_level_values(H.KC_ID_LEVEL)])


def save_stage_f(tag: str, cond: str, st: dict) -> None:
    """A stage run's matrices under outputs/.../F/ as E90 saves them (<tag>_<on|acute_off>_<stage>_<win|pre>.parquet,
    _apl.csv), without the responses and counts: they are the archived run's of the source experiment exactly (checked
    when made; `f_stage_identity`), which `load_stage_f` adds back."""
    d = _f_dir()
    cc = 'on' if cond == 'on' else 'acute_off'
    for k in STAGE_KEYS:
        x = st[k].copy()
        x.columns = [str(c) for c in x.columns]
        if isinstance(x.index, pd.MultiIndex):
            x.index = x.index.set_names([str(n) for n in x.index.names])
        x.to_parquet(d / f'{tag}_{cc}_{k}.parquet')
    st['apl'].to_csv(d / f'{tag}_{cc}_apl.csv')


def load_stage_f(tag: str, cond: str, source: Optional[Tuple[str, str]] = None) -> Optional[dict]:
    """A stage run's matrices (`save_stage_f`; the local outputs, else the archive), with the responses and counts of
    its source run (`f_load(*source)`) when given; None if absent."""
    cc = 'on' if cond == 'on' else 'acute_off'
    for d in (_f_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / F_DIR):
        if (d / f'{tag}_{cc}_claw_win.parquet').exists():
            st = {k: pd.read_parquet(d / f'{tag}_{cc}_{k}.parquet') for k in STAGE_KEYS}
            st['apl'] = pd.read_csv(d / f'{tag}_{cc}_apl.csv', index_col=0)
            if source is not None:
                o = f_load(*source)[cond]
                st['responses'], st['spike_counts'] = o[0], o[1]
            return st
    return None


def _identity(got: Tuple[pd.DataFrame, pd.DataFrame], ref: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    nd = int((got[1].to_numpy(dtype=float) != ref[1].to_numpy(dtype=float)).sum())
    nr = int((got[0].to_numpy(dtype=float) != ref[0].to_numpy(dtype=float)).sum())
    same_index = bool(got[1].index.equals(ref[1].index) and list(map(str, got[1].columns)) == list(map(str, ref[1].columns)))
    return dict(counts_differ=nd, responses_differ=nr, index_equal=same_index, exact=bool(nd == 0 and nr == 0 and same_index))


def f_stage_refit(w: str) -> Tuple[Dict[str, dict], dict]:
    """E91's part D fit on wiring w (s = F_POOL, exponent F_K) run again at its parameters (E72's `fitted_params` of
    E91's fits.csv) on E91's simulator with the pooled share (`build_simulator91`), APL on and off, APL's activity saved;
    each run's counts and responses must equal the archived fit's exactly; E90's stage matrices (`read_stages`)."""
    from experiments import E72_wiring_population as E72
    from experiments import E90_model_stages as E90
    tag = E91.d_tag(w, F_K, F_POOL)
    saved = f_load('E91fit', tag)
    p = E72.fitted_params(F_SOURCE, tag)
    sim, info = E91.build_simulator91(w, claw_inh_pool=F_POOL)
    cells = info['cells']
    mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
    order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
    if float(mp.kc.claw_inh_pool) != F_POOL:
        raise RuntimeError(f'{w}: the simulator\'s claw_inh_pool is {mp.kc.claw_inh_pool}')
    W = info['W']
    if not np.array_equal(np.asarray(rv.kc.wPNKC, dtype=float), W.to_numpy(dtype=float)):
        raise RuntimeError(f"{w}: olfsysm's wPNKC differs from the wiring's claw x bouton matrix")
    mp.kc.save_inh_sims = True
    out, ident = {}, dict(source=f'E91fit:{tag}', params=p, claw_inh_pool=F_POOL)
    for c in ('on', 'off'):
        r, s = sim(p if c == 'on' else dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0))
        ident[c] = _identity((r, s), saved[c])
        if not ident[c]['exact']:
            raise RuntimeError(f'{w} refit {c}: the run differs from the archived fit: {ident[c]}')
        st = E90.read_stages(mp, rv, W, order, [H.odor_abbrev(x) for x in r.columns], E90.MODEL_WINDOW)
        st.update(responses=r, spike_counts=s)
        out[c] = st
    ident['n_steps'] = out['on']['n_steps']
    return out, ident


def f_stage_b2(w: str, n_workers: int = 16) -> Tuple[Dict[str, dict], dict]:
    """E91's choice on the base with the flies' recruitment (F_B2_CHOICE: exponent 1, s = 1) on wiring w run again odor by
    odor (`FPool`) at its per-odor threshold factors and APL scales as run (`b2_odor_params` of its runs row), with the
    weights at exponent 1 (E89's `apl_kc_weights`), APL on and off; the counts and responses must equal the archived
    run's exactly; E90's stage matrices, odor by odor (`read_stages_odor`)."""
    tag = E91.run_tag(F_B2_CHOICE, w)
    saved = f_load('E91', tag)
    row = e91_run_row(tag)
    pv = E91.parse_variant(F_B2_CHOICE)
    parts = E83.wiring_parts(w)
    wk = E89.apl_kc_weights(parts['wnorm']['wAPLKC'], parts['breadth'], pv['exponent'])
    if str(row.get('wAPLKC_sha256')) != E89.E80_sha(wk):
        raise RuntimeError(f'{tag}: the APL->KC weights are not the run\'s')
    pool = FPool(w, n_workers)
    try:
        thr, apl = b2_odor_params(row, pool.cols)
        recs = pool.run_f(thr, apl, wk, pv['s'])
    finally:
        pool.close()
    odors = [H.odor_abbrev(c) for c in pool.cols]
    out, ident = {}, dict(source=f'E91:{tag}', thr=dict(zip(pool.cols, thr)), apl=dict(zip(pool.cols, apl)),
                          claw_inh_pool=pv['s'], exponent=pv['exponent'])
    for c in ('on', 'off'):
        like = saved[c]
        R = pd.DataFrame(np.column_stack([r[c]['responses'] for r in recs]), index=like[0].index, columns=like[0].columns)
        S = pd.DataFrame(np.column_stack([r[c]['counts'] for r in recs]), index=like[1].index, columns=like[1].columns)
        ident[c] = _identity((R, S), like)
        if not ident[c]['exact']:
            raise RuntimeError(f'{w} {tag} {c}: the run differs from the archived one: {ident[c]}')
        st = stage_frames([r[c] for r in recs], odors, parts['W'])
        st.update(responses=R, spike_counts=S)
        out[c] = st
    return out, ident


F_IDENTITY = 'identity_F.json'


def part_f3(out: Path, families: Sequence[str] = F_STAGE_FAMILIES, n_workers: int = 16) -> dict:
    """F3: the stage runs (F_STAGE_FAMILIES on S, resumable: a run whose matrices exist is skipped), each reproducing its
    archived run exactly (identity_F.json) before its matrices are saved (`save_stage_f`) with its units' metadata."""
    import gc
    path = out / F_IDENTITY
    ids = json.loads(path.read_text()) if path.exists() else {}
    cc = None
    for fam in families:
        for w in E83.S_WIRINGS:
            tag = f_stage_tag(fam, w)
            if load_stage_f(tag, 'on') is not None and load_stage_f(tag, 'off') is not None and tag in ids:
                continue
            cc = _commits() if cc is None else cc
            t0 = time.time()
            st, ident = f_stage_refit(w) if fam == 'E91_refit' else f_stage_b2(w, n_workers)
            for c in ('on', 'off'):
                save_stage_f(tag, c, st[c])
            (_f_dir() / f'{tag}_units.json').write_text(json.dumps(units_meta(E83.wiring_parts(w)['W'])) + '\n')
            ids[tag] = dict(ident, family=fam, wiring=w, seconds=round(time.time() - t0, 1), **cc)
            _write_json(path, ids)
            print(f'F3 {tag}: exact ({ident["on"]}, {ident["off"]}), {time.time() - t0:.0f} s', flush=True)
            del st
            gc.collect()
    return ids


# --- F4: the claws and the two steps

def f_claw_pairs(families: Sequence[str] = F_STAGE_FAMILIES) -> pd.DataFrame:
    """E90's readings (`read_wiring`) of the stage runs' boutons and claws, every subset and share, APL on and off
    (outputs/.../F/stage_pairs.parquet); one generator per run (E90's NOISE_SEED + F_NOISE_SEED + its index)."""
    from experiments import E90_model_stages as E90
    path = _f_dir() / 'stage_pairs.parquet'
    rows = []
    i = 0
    for fam in families:
        src = 'E91fit' if fam == 'E91_refit' else 'E91'
        for w in E83.S_WIRINGS:
            tag = f_stage_tag(fam, w)
            stag = E91.d_tag(w, F_K, F_POOL) if fam == 'E91_refit' else E91.run_tag(F_B2_CHOICE, w)
            st = {c: load_stage_f(tag, c, (src, stag)) for c in ('on', 'off')}
            if any(v is None for v in st.values()):
                raise RuntimeError(f'{tag}: no stage matrices (run --part F first)')
            meta = None
            for d in (_f_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / F_DIR):
                if (d / f'{tag}_units.json').exists():
                    meta = {k: np.asarray(v) for k, v in json.loads((d / f'{tag}_units.json').read_text()).items()}
                    break
            if meta is None:
                raise RuntimeError(f'{tag}: no units metadata')
            rr = E90.read_wiring(w, st['on'], st['off'], meta, np.random.default_rng(E90.NOISE_SEED + F_NOISE_SEED + i))
            rows.append(pd.DataFrame(rr).assign(family=fam, tag=tag))
            i += 1
    df = pd.concat(rows, ignore_index=True)
    df.to_parquet(path)
    return df


def claw_groups(pairs: pd.DataFrame) -> pd.DataFrame:
    """The claws as part A's steps read them (all claws, the drive alone, APL on; 'session_z' and 'noisy_measured'),
    group means per family and wiring."""
    c = pairs[(pairs.stage == 'claw') & (pairs.subset == 'all') & (pairs.share == 0.0) & (pairs.cond == 'on')
              & pairs.reading.isin(['session_z', 'noisy_measured'])]
    return _A()._groups_of_pairs(c, ['family', 'wiring', 'reading'])


def f_step_rows(sm: pd.DataFrame, claws: pd.DataFrame, passes: Dict[Tuple[str, str, str], bool], claw_ok: bool,
                flies_step=None) -> pd.DataFrame:
    """The two steps like for like (part A's `step_verdict_rows` on S's populations): per family with claws (`claws`:
    family, wiring, reading, groups) and like-for-like reading, the somata (`sm`, S's population and its runs) minus the
    claws (part A's CLAW_READING: the session z without noise, the noisy reading as measured through it), against the
    flies' steps (`flies_step(units, genotype, reading, measure, claw reading)`, part A's `flies_step_groups`), 3 of 4
    wirings on the side; cond 'on' the step with APL on, 'off' the step without APL, both from the claws with APL on. A
    noise-removed comparison is read if the group passed A1's check and E90's check of the claws' noise removal
    passed (`claw_ok`)."""
    A = _A()
    if flies_step is None:
        cache: Dict[tuple, pd.DataFrame] = {}
        fc = A.flies_claws()

        def flies_step(u, g, fr, ms, fcr):
            k = (u, g, fr, ms, fcr)
            if k not in cache:
                cache[k] = A.flies_step_groups(A.flies_samples(u)[(g, fr, ms)], fc[fcr])
            return cache[k]
    ws = list(E83.S_WIRINGS)
    G = list(A.GROUPS)
    rows = []
    for fam in dict.fromkeys(claws.family):
        for prefix, u, fr, ms in A.like_for_like():
            mcr, fcr = A.CLAW_READING[prefix.split('_')[-1]]
            reading = f'{prefix}_{ms}'
            cl = claws[(claws.family == fam) & (claws.reading == mcr)].set_index('wiring')
            if set(cl.index) != set(ws):
                raise RuntimeError(f'{fam}: claws on {sorted(cl.index)}')
            cl = cl.loc[ws, G]
            for cond in A.CONDS:
                x = sm[(sm.family == fam) & (sm.reading == reading) & (sm.cond == cond)]
                p = x[(x.kind == 'population') & (x.wiring == 'S')]
                wr = x[(x.kind == 'run') & x.wiring.isin(ws)].set_index('wiring')
                if len(p) != 1 or len(wr) != len(ws):
                    raise RuntimeError(f'{fam} {reading} {cond}: {len(p)} population rows, {len(wr)} runs')
                w_step = wr.loc[ws, G] - cl
                p_step = {g: float(p.iloc[0][g]) - float(cl[g].mean()) for g in G}
                fl = flies_step(u, A.GENO[cond], fr, ms, fcr)
                for r in A.compare(p_step, w_step, fl, n_same=F_SAME_SIDE_4):
                    read = bool(fr != 'noise_removed' or (passes.get((r['group'], ms, u), False) and claw_ok))
                    rows.append(dict(family=fam, population='S', element='step', set='like_for_like', cond=cond,
                                     reading=reading, flies_units=u, flies_reading=fr, measure=ms, read=read,
                                     claws_model=mcr, claws_flies=fcr, claw_model=float(cl[r['group']].mean()), **r))
    return pd.DataFrame(rows)


def part_f4(out: Path, sm: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """F4: the claws of the stage runs (`f_claw_pairs`) and C2's on S (E90's, part A's `model_claws`), claws_F.csv;
    steps_F.csv (`f_step_rows`) for C2's family on S and each stage family."""
    A = _A()
    sm = pd.read_csv(out / 'somata_F_means.csv') if sm is None else sm
    t0 = time.time()
    path = _f_dir() / 'stage_pairs.parquet'
    pairs = pd.read_parquet(path) if path.exists() else f_claw_pairs()
    cl = claw_groups(pairs)
    c2 = A.model_claws()
    c2 = c2[c2.wiring.isin(E83.S_WIRINGS)].assign(family='C2')
    claws = pd.concat([c2, cl], ignore_index=True)
    claws.to_csv(out / 'claws_F.csv', index=False)
    e90v = json.loads((H.results_dir('E90_model_stages') / 'verdict.json').read_text())
    rows = f_step_rows(sm, claws, A.check_passes(), bool(e90v['noise_check_passed']['claw']))
    rows.to_csv(out / 'steps_F.csv', index=False)
    print(f'F4 in {time.time() - t0:.0f} s', flush=True)
    return rows


# --- the verdict, the table and the figure

def _cell(e: dict) -> str:
    sg = ''.join('+' if x > 0 else '-' for x in e['signs'])
    return f"{e['n_differ']}/{e['n_readings']}" + (f' ({sg})' if sg else '')


def verdict_f(out: Path) -> dict:
    """verdict_F.json: per family and population the somata's verdicts per group (part A's `verdicts` over the
    like-for-like readings read), per stage family the steps', each family's place against the line and its scores
    (line_F.csv), the stage runs' identities."""
    A = _A()
    so = pd.read_csv(out / 'somata_F.csv')
    res: dict = dict(families=F_FAMILIES, runs=f_runs().to_dict(orient='records'), somata={}, steps={}, line={})
    for (fam, pop), x in so[so.set == 'like_for_like'].groupby(['family', 'population'], sort=False):
        res['somata'][f'{fam}|{pop}'] = {g: dict(verdict=v['verdict'], on=_cell(v['on']), off=_cell(v['off']),
                                                 mean_diff_on=v['mean_diff_on'], mean_diff_off=v['mean_diff_off'])
                                         for g, v in A.verdicts(x).items()}
    if (out / 'steps_F.csv').exists():
        st = pd.read_csv(out / 'steps_F.csv')
        for fam, x in st.groupby('family', sort=False):
            res['steps'][fam] = {g: dict(verdict=v['verdict'], on=_cell(v['on']), off=_cell(v['off']),
                                         mean_diff_on=v['mean_diff_on'], mean_diff_off=v['mean_diff_off'])
                                 for g, v in A.verdicts(x).items()}
    if (out / 'line_F.csv').exists():
        lc = pd.read_csv(out / 'line_F.csv')
        # C2's own residual on S (the line's 13 weightings include it): a place relative to C2 on S beside the rule's
        c2s = {r.readout: float(r.resid) for r in lc[(lc.folder == 'E90_model_stages') & (lc.population == 'S')].itertuples()}
        for fam, keys in F_LINE_OF.items():
            keys = keys if isinstance(keys[0], tuple) else (keys,)
            for folder, lf in keys:
                for r in lc[(lc.folder == folder) & (lc.family == lf) & lc.population.isin(['S', 'T'])].itertuples():
                    res['line'][f'{fam}|{r.population}|{r.readout}'] = dict(
                        folder=folder, family=lf, ester_rel=r.ester_rel, alcohol_rel=r.alcohol_rel, resid=r.resid,
                        resid_vs_c2=float(r.resid - (c2s.get(r.readout, np.nan) if r.population == 'S' else 0.0)),
                        off_toward_flies=bool(r.off_toward_flies), pattern=r.pattern, amplitude=r.amplitude, shape=r.shape,
                        sd_ratio=r.sd_ratio, pearson_r=r.pearson_r, rho=r.rho, loss=r.loss, mean_effect=r.mean_effect,
                        err_ester=r.err_ester, err_alcohol=r.err_alcohol, in_bands=bool(r.in_bands))
        off = lc[lc.in_bands & lc.off_toward_flies & (lc.kind != 'line') & (lc.folder == F_SOURCE)]
        res['e91_off_line'] = off[['readout', 'family', 'population', 'resid', 'pattern', 'loss', 'mean_effect']].to_dict(
            orient='records')
    if (out / F_IDENTITY).exists():
        res['identity'] = {k: {c: v.get(c) for c in ('on', 'off', 'source', 'seconds')}
                           for k, v in json.loads((out / F_IDENTITY).read_text()).items()}
    _write_json(out / 'verdict_F.json', res)
    return res


F_TABLE_READINGS = ('saturating_z_none_cosine', 'linear_z_none_cosine', 'saturating_z_none_pearson',
                    'linear_z_none_pearson')


def comparison_table_f(out: Path) -> str:
    """comparison_table_F.md from part F's tables."""
    A = _A()
    v = verdict_f(out)
    G = A.TABLE_GROUPS
    L = ["# E92 part F: the summed form (E91) read as part A read C2's family", '',
         f'Generated by experiments/E92_gain.py at {H.code_commits().get("queue_commit", "")[:7]}.', '',
         '## Families', '', '| family | what | runs |', '|---|---|---|']
    runs = f_runs()
    for fam, what in F_FAMILIES.items():
        r = runs[runs.family == fam]
        L.append(f"| {fam} | {what} | {', '.join(r.source + ':' + r.tag)} |")
    if v['line']:
        L += ['', '## The main target and the line (S and T populations; calibrated readout and raw counts)', '',
              '| family | pop. | readout | pattern | amplitude | shape | spread ratio | Pearson r | ρ | mean APL effect | loss | '
              'ester / alcohol error | ester rel. | alcohol rel. | residual | residual relative to C2 | off toward the flies |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
        for k, x in v['line'].items():
            fam, pop, ro = k.split('|')
            L.append(f"| {fam} | {pop} | {ro} | {x['pattern']:.4f} | {x['amplitude']:.4f} | {x['shape']:.4f} | "
                     f"{x['sd_ratio']:.2f} | {x['pearson_r']:.2f} | {x['rho']:.2f} | {x['mean_effect']:+.3f} | {x['loss']:.4f} | "
                     f"{x['err_ester']:+.3f} / {x['err_alcohol']:+.3f} | {x['ester_rel']:+.3f} | {x['alcohol_rel']:+.3f} | "
                     f"{x['resid']:+.3f} | {x['resid_vs_c2']:+.3f} | {x['off_toward_flies']} |")
        L += ['', f"E91's variants off the line toward the flies (inside the bands): {json.dumps(v.get('e91_off_line', []))}"]
    L += ['', '## The somata: like-for-like verdicts (differ / read, with the sign)', '',
          '| family | pop. | ' + ' | '.join(f'{A.TABLE_HEAD[g]}: on | off' for g in G) + ' |',
          '|---|---|' + '---:|---:|' * len(G)]
    for k, x in v['somata'].items():
        fam, pop = k.split('|')
        L.append(f'| {fam} | {pop} | ' + ' | '.join(f"{x[g]['on']} | {x[g]['off']}" for g in G) + ' |')
    sm = pd.read_csv(out / 'somata_F_means.csv')
    so = pd.read_csv(out / 'somata_F.csv')
    L += ['', "## The somata, S's populations, in four like-for-like readings (the flies' noise-removed values in Yang's z "
          "beside: the control flies with APL on, APL output blocked (TNT) with APL off)", '',
          '| APL | reading | family | ' + ' | '.join(A.TABLE_HEAD[g] for g in G) + ' |', '|---|---|---|' + '---:|' * len(G)]
    for cond in A.CONDS:
        for rd in F_TABLE_READINGS:
            fl = so[(so.reading == rd) & (so.cond == cond) & (so.population == 'S') & (so.family == 'C2')].set_index('group')
            L.append(f'| {cond} | {rd} | flies | ' + ' | '.join(_f(fl.loc[g, 'flies'], '.3f') for g in G) + ' |')
            for fam in F_FAMILIES:
                x = sm[(sm.kind == 'population') & (sm.wiring == 'S') & (sm.family == fam) & (sm.reading == rd)
                       & (sm.cond == cond)]
                if len(x):
                    L.append(f'| {cond} | {rd} | {fam} | ' + ' | '.join(_f(x.iloc[0][g], '.3f') for g in G) + ' |')
    if v['steps']:
        L += ['', "## The two steps like for like (S's populations, 3 of 4 wirings on the side)", '',
              '| family | ' + ' | '.join(f'{A.TABLE_HEAD[g]}: with APL on | without APL' for g in G) + ' |',
              '|---|' + '---:|---:|' * len(G)]
        for fam, x in v['steps'].items():
            L.append(f'| {fam} | ' + ' | '.join(f"{x[g]['on']} | {x[g]['off']}" for g in G) + ' |')
        cl = pd.read_csv(out / 'claws_F.csv')
        L += ['', "## The claws with APL on (all claws, the drive alone), S's wirings' mean", '',
              '| family | reading | ' + ' | '.join(A.TABLE_HEAD[g] for g in G) + ' |', '|---|---|' + '---:|' * len(G)]
        for (fam, rd), x in cl.groupby(['family', 'reading'], sort=False):
            L.append(f'| {fam} | {rd} | ' + ' | '.join(_f(x[g].mean(), '.3f') for g in G) + ' |')
    if 'identity' in v:
        L += ['', '## The stage runs: identities against the archived runs', '', '```',
              json.dumps(v['identity'], indent=1, default=str), '```']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table_F.md').write_text(text)
    return text


def figure_f(out: Path) -> Optional[Path]:
    """E92_F_line.png: A4's line (S's line through E89's 13 weightings, with the weightings' largest residual shaded)
    with C2's family and E91's S populations (the families of F named), and the control flies' points, raw counts and
    the calibrated readout."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not (out / 'line_F.csv').exists():
        return None
    A = _A()
    lc = pd.read_csv(out / 'line_F.csv')
    ink, ink2, grid = '#0b0b0b', '#52514e', '#e6e5e1'
    cols = {'C2': '#0b0b0b', 'E91_B1_k2_s1': '#2a78d6', 'E91_refit': '#c2449f', 'E91_B2_k2': '#5b4cc4',
            'E91_B2_k1_s1': '#138a60'}
    labels = {'C2': "C2's family (clipped)", 'E91_B1_k2_s1': 'summed form, the fits as they are',
              'E91_refit': 'summed form, the family fitted again', 'E91_B2_k2': "the flies' recruitment base, clipped",
              'E91_B2_k1_s1': "that base's choice: exponent 1, summed"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, ro in zip(axes, A.READOUTS):
        ln = A.line(ro)['S']
        x = lc[(lc.readout == ro) & (lc.population == 'S') & (lc.folder == F_SOURCE)]
        ax.scatter(x.alcohol_rel, x.ester_rel, s=14, c='#9aa7b4', label="E91's other S populations", zorder=2)
        pts = np.array(list(ln['points'].values()))
        ax.scatter(pts[:, 0], pts[:, 1], s=22, c='#7fb3e8', label="E89's 13 weightings", zorder=3)
        xs = np.linspace(min(pts[:, 0].min(), x.alcohol_rel.min()) - 0.02, max(pts[:, 0].max(), x.alcohol_rel.max()) + 0.02, 10)
        ys = ln['intercept'] + ln['slope'] * xs
        ax.plot(xs, ys, '-', color='#7fb3e8', lw=1)
        ax.fill_between(xs, ys - ln['max_abs_resid'], ys + ln['max_abs_resid'], color='#7fb3e8', alpha=0.12, lw=0,
                        label="the weightings' largest residual")
        for fam, keys in F_LINE_OF.items():
            keys = keys if isinstance(keys[0], tuple) else (keys,)
            folder, lf = keys[0]
            r = lc[(lc.readout == ro) & (lc.population == 'S') & (lc.folder == folder) & (lc.family == lf)]
            if len(r):
                ax.scatter(r.alcohol_rel, r.ester_rel, s=60, marker='D', facecolors='white', edgecolors=cols[fam], lw=1.6,
                           label=labels[fam], zorder=5)
        for k, mk in (('trial_mean', '*'), ('noise_removed', 'P')):
            f_ = ln['flies'][f'{k}_all']
            ax.plot(f_['alcohol'], f_['ester'], mk, color=ink if k == 'trial_mean' else ink2, ms=11,
                    label=f"control flies, {k.replace('_', ' ')}")
        ax.set_xlabel("alcohol pairs' APL-on distance minus the mean pair's (cosine)", fontsize=9)
        ax.set_ylabel("ester pairs' APL-on distance minus the mean pair's (cosine)", fontsize=9)
        ax.set_title(f"{'raw counts' if ro == 'raw' else 'calibrated readout'}: S's populations", fontsize=10)
        ax.grid(True, color=grid, lw=0.6)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    path = out / 'E92_F_line.png'
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def part_f(out: Path, workers: int = 4, n_workers: int = 16, stages: bool = True) -> dict:
    """Part F (E91's note posted 2026-11-13 19:38): F1 the somata of E91's families in A2's readings (somata_F.csv), F2
    their place against A4's line (line_F.csv), F3 the stage runs on S (E91's refits and its choice on the base with the
    flies' recruitment, each reproducing its archived run exactly), F4 the two steps like for like (steps_F.csv); then
    verdict_F.json, comparison_table_F.md and E92_F_line.png. With `stages` False, F3 and F4 are skipped."""
    sm, _ = part_f1(out, workers)
    part_f2(out, workers)
    if stages:
        part_f3(out, n_workers=n_workers)
        part_f4(out, sm)
    comparison_table_f(out)
    figure_f(out)
    return verdict_f(out)


# ---------------------------------------------------------------------------
# the archive, the figure and the table

def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): one row per run of this experiment (B, C and D's starts) with absolute
    parameters (the APL scales as run), the form, the weighting, the control's static inhibition, the factors' range,
    the divisor and realism; D's fits' distances, which E10_joint_fit.py wrote, are kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E92']
    rows, parts = [], []
    for _, r in runs.iterrows():
        form = {'sub': 'subtractive (max(x - i, 0))', 'div': 'divisive at the claws (x / (1 + i); olfsysm apl_inh_divisive)',
                'bdiv': 'divisive at the boutons (r / (1 + i); olfsysm bouton_inh_divisive)',
                'both': 'divisive at the claws and at the boutons'}[r.form]
        ctl = str(r.control_of) not in ('', 'nan')
        kind = (f"the fit as it is; APL->KC weights x ((F + 1) / mean(F + 1))^{r.exponent:g}; APL's inhibition {form}"
                + (f"; the control of {r.control_of}: a static inhibition at the soma (ext_inh_kc) per KC, its mean loss of "
                   "drive under the gain" if ctl else '') + "; the APL scale re-held; APL-off run the fit's")
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=ARM, kind=kind, config=r.variant, form=r.form,
                         breadth_def='F', breadth_k=r.exponent, k=r.exponent, apl_inh_divisive=bool(r.apl_inh_divisive),
                         bouton_inh_divisive=bool(r.bouton_inh_divisive), wiring=r.wiring, ref_experiment=r.fit_folder,
                         ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale,
                         wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale,
                         apl_factor=r.apl_factor, scale_min=r.apl_min, scale_max=r.apl_max, rate_on=r.rate_on,
                         rate_off=r.rate_off, max_count_on=r.max_count_on, max_count_off=r.max_count_off,
                         wAPLKC_sha256=r.wAPLKC_sha256, ext_inh_kc_max=r.get('ext_inh_kc_max'), control_of=r.control_of,
                         divisor_kc_max=r.get('divisor_kc_max'), realism_flag=r.realism_flag))
        o = load_run(r.tag)
        parts.append(pd.concat({(r.tag, 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r.tag, 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
    df = pd.DataFrame(rows)
    df.to_csv(out / 'variants.csv', index=False)
    long_path = out / 'model_pair_distances_long.csv'
    keep = []
    if long_path.exists():
        old = pd.read_csv(long_path)
        keep = [old[old.tag.isin(set(_fits(out).tag))]]           # D's fits, written by E10_joint_fit.py
    pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def figure(out: Path) -> Optional[Path]:
    """E92_gain_screen.png: S's pattern part under the readout against half A, half B and all recordings by setting
    (the references, B1 and B2; the family's reference dashed; hollow: a setting with a run outside the bands or at a
    bracket's bound)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not (out / 'readings_B.csv').exists():
        return None
    rd = pd.read_csv(out / 'readings_B.csv')
    g = _settings92(rd, ('ref', 'B1', 'B2'))
    if g.empty:
        return None
    order = {'ref': 0, 'B1': 1, 'B2': 2}
    g = g.assign(_o=g.screen.map(order)).sort_values(['_o', 'form', 'exponent'], ascending=[True, True, False])
    fig, ax = plt.subplots(figsize=(10, 4.8))
    x = np.arange(len(g))
    bad = ((g.n_flagged > 0) | (g.n_at_bound > 0)).to_numpy()
    for t, col, lab in (('halfA', '#d1495b', 'half A'), ('halfB', '#2e86ab', 'half B'), ('z8', '0.2', 'all recordings')):
        y = g[f'{t}_pattern'].to_numpy()
        ax.plot(x, y, '-', color=col, lw=0.8, label=lab)
        ax.scatter(x, y, s=28, facecolors=['none' if b else col for b in bad], edgecolors=col)
        ref = g[g.variant == ref_variant()][f'{t}_pattern']
        if len(ref):
            ax.axhline(float(ref.iloc[0]), color=col, ls='--', lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([v.split('_', 1)[1] for v in g.variant], rotation=60, fontsize=8)
    ax.set_ylabel('pattern part of the APL-effect loss')
    ax.set_title("APL's inhibition as a gain, the fits as they are: S's population\n"
                 "(calibrated readout; dashed: C2's family; hollow: a run outside the bands or at a bound)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / 'E92_gain_screen.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _f(x, spec: str = '.4f') -> str:
    return E87._f(x, spec)


def comparison_table(out: Path) -> str:
    """comparison_table_B.md from the parts' tables."""
    L = ["# E92 parts B to D: APL's inhibition of the KCs as a gain on the claw's drive", '',
         f'Generated by experiments/E92_gain.py at {H.code_commits().get("queue_commit", "")[:7]}.', '']
    if (out / 'smoke.json').exists():
        s = json.loads((out / 'smoke.json').read_text())
        L += ['## Smoke run', '', '```', json.dumps({k: x for k, x in s.items() if k != 'row'}, indent=1, default=str), '```', '']
    if (out / 'readings_B.csv').exists():
        rd = pd.read_csv(out / 'readings_B.csv')
        runs = runs_table(out)
        g = _settings92(rd, ('ref', 'B1', 'B2', 'control'))
        graw = rd[(rd.unit == 'S') & (rd.readout == 'raw')].set_index('variant')
        rr = runs.groupby('variant')
        L += ["## B: S's population (calibrated readout; the parts of the pattern part against all recordings)", '',
              '| setting | pattern, half A | half B | all | amplitude | shape | spread ratio | Pearson r | ρ | loss | '
              'mean APL effect | raw-count loss | ester / alcohol error | APL scale factor | divisor p50 / max | flagged / at bound |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for r in g.itertuples():
            q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()
            fac = f"{q.apl_factor.min():.3g}-{q.apl_factor.max():.3g}" if len(q) and 'apl_factor' in q else '–'
            dv = (f"{_f(q.divisor_kc_p50.median(), '.3g')} / {_f(q.divisor_kc_max.max(), '.3g')}"
                  if len(q) and 'divisor_kc_p50' in q and q.divisor_kc_p50.notna().any() else '–')
            L.append(f'| {r.variant} | {r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | '
                     f'{r.z8_amplitude:.4f} | {r.z8_shape:.4f} | {r.z8_sd_ratio:.2f} | {r.z8_pearson:.2f} | {r.z8_rho:.2f} | '
                     f'{r.z8_loss:.4f} | {r.mean_effect:+.3f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | '
                     f'{r.z8_err_ester:+.3f} / {r.z8_err_alcohol:+.3f} | {fac} | {dv} | {r.n_flagged} / {r.n_at_bound} |')
        L.append('')
        mine = runs[runs.source == 'E92']
        if len(mine):
            L += ['## Realism of the runs (rule of 2026-11-10 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | '
                  'divisor p10 / p50 / p90 / max | flag | flag (earlier rule) |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|---|---|']
            for r in mine.itertuples():
                dv = ' / '.join(_f(getattr(r, f'divisor_kc_{n}', np.nan), '.3g') for n, _ in QUANTILES)
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | '
                         f'{r.opk_on_tv:.3f} / {r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | {dv} | '
                         f'{r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
    for name in ('beside_B.csv', 'choice.csv', 'gate.csv', 'equal.csv', 'd_scores.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                L += [f'## {name}', '', '(none)', '']
                continue
            if name == 'beside_B.csv':
                # (the main columns; the file has every reading)
                keep = [c for c in d.columns if c == 'variant' or re.fullmatch(
                    r'err_(ester|alcohol)_z8_(cal|raw)|line_(raw|cal)_(ester_rel|alcohol_rel|resid|off_toward_flies(_nr)?)'
                    r'|a2_saturating_z_(none|expected)_(pearson|cosine)_(ester|alcohol)_on', c)]
                d = d[keep]
            try:
                txt = d.to_markdown(index=False, floatfmt='.4f')
            except ImportError:
                txt = '```\n' + d.to_string() + '\n```'
            L += [f'## {name}', '', txt, '']
    if (out / 'verdict_B.json').exists():
        vc = json.loads((out / 'verdict_B.json').read_text())
        short = {k: vc.get(k) for k in ('verdict', 'gate', 'B2', 'D_plan', 'D', 'flagged_runs', 'at_bound_runs')}
        short['equals'] = [e for e in vc.get('equals', []) if e.get('equal')]
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table_B.md').write_text(text)
    return text


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json: the script, the commits (H.code_commits()), olfsysm's build, the model's provenance with
    the lab's input (H.model_provenance) for the simulation parts, the settings and `extra`; a part run again adds a
    dated record."""
    import olfsysm
    path = out / f'provenance_{part}.json'
    old = json.loads(path.read_text()) if path.exists() else {}
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('B1', 'B2', 'B3', 'smoke', 'C', 'D', 'F'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    _write_json(path, dict(script='experiments/E92_gain.py', part=part,
                           settings=dict(FORMS=FORMS, B1_SETTINGS=B1_SETTINGS, B2_SETTINGS=B2_SETTINGS, PULSE=PULSE,
                                         APL_BRACKET=APL_BRACKET, APL_STEP=APL_STEP, DIVISOR_HIGH=DIVISOR_HIGH,
                                         CO_CANDIDATE=CO_CANDIDATE, CHOICE_SCREENS=CHOICE_SCREENS, D_MAXFEV=D_MAXFEV,
                                         **(dict(F_SOURCE=F_SOURCE, F_K=F_K, F_POOL=F_POOL, F_FAMILIES=F_FAMILIES,
                                                 F_STAGE_FAMILIES=F_STAGE_FAMILIES, F_SAME_SIDE_4=F_SAME_SIDE_4,
                                                 F_NOISE_SEED=F_NOISE_SEED, F_LINE_E89=F_LINE_E89)
                                            if part == 'F' else {})),
                           runs=runs))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('B1', 'B2', 'B3', 'smoke', 'co', 'C', 'score', 'D', 'Dscore', 'E', 'F', 'Ftable', 'variants',
                             'table'))
    ap.add_argument('--wirings', default='S', help="B1 / B2 / B3: 'S', 'T' or a comma-separated list")
    ap.add_argument('--workers', type=int, default=4, help='co / score: processes (4)')
    ap.add_argument('--reference', default='C2', choices=REFERENCE_FAMILIES,
                    help="C / score: the reference family (C2's clipped family; 'E91_refits' is a stub)")
    ap.add_argument('--pop', default='S', choices=('S', 'T'), help='D: which fits')
    ap.add_argument('--k', type=float, default=None, help="D: the exponent (default: every one D's plan names)")
    ap.add_argument('--d-start', default='reheld', choices=D_STARTS, help="D: the fits' start")
    ap.add_argument('--force', action='store_true', help="D: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='D: print the commands only')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    ap.add_argument('--f-workers', type=int, default=16, help="F: the B2 stage runs' one-odor worker processes (16)")
    ap.add_argument('--no-stages', action='store_true', help='F: skip the stage runs and the steps (F3, F4)')
    args = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    out = H.results_dir(EXPERIMENT)
    ws = _wirings(args.wirings)
    if args.part in ('B1', 'B2', 'B3', 'smoke', 'C', 'D', 'F') and not (args.part == 'D' and args.dry_run):
        provenance(out, args.part, wirings=list(ws), reference=args.reference,
                   **(dict(pop=args.pop, k=args.k, d_start=args.d_start) if args.part == 'D' else {}))
    if args.part == 'B1':
        part_b1(out, ws)
    elif args.part == 'B2':
        part_b2(out, ws)
    elif args.part == 'B3':
        part_b3(out, ws)
    elif args.part == 'smoke':
        part_smoke(out)
    elif args.part == 'co':
        coresponse92(out, args.workers)
    elif args.part == 'C':
        part_c(out, args.reference)
    elif args.part == 'score':
        score(out, args.workers, coresponse=not args.no_co, family=args.reference)
    elif args.part == 'D':
        part_d(out, args.pop, args.k, args.d_start, force=args.force, dry_run=args.dry_run)
    elif args.part == 'Dscore':
        score_d(out)
        verdict92(out)
        comparison_table(out)
    elif args.part == 'E':
        part_e(out)
    elif args.part == 'F':
        part_f(out, args.workers, args.f_workers, stages=not args.no_stages)
    elif args.part == 'Ftable':
        comparison_table_f(out)
        figure_f(out)
    elif args.part == 'variants':
        write_variants(out)
    elif args.part == 'table':
        figure(out)
        comparison_table(out)


if __name__ == '__main__':
    main()
