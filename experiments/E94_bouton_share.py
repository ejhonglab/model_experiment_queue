#!/usr/bin/env python3
"""E94: the PN boutons' share of APL's inhibition in C2's family, from none to four times the fitted dose, as a
subtraction and as a gain, read with the gate for fitted changes; and the family fitted again with that share at the
connectome's.

The entry is item E94 of experiment_priorities.md (proposed by viki 2026-11-14 01:20; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-14 00:51).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `WiringModel`, through E91's `WiringModel91`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the
recordings E77's, C / D the second split's (E89's `targets89`). A dose m multiplies the wiring's fitted APL->bouton scale
(`wAPLPN_scale`) and nothing else (`WiringModel94.run_on`); the tied APL scale (APL->KC, KC->APL, bouton->APL; E83's
TIED_ON) is then bisected to the wiring's APL-on count in E92's bracket (APL_BRACKET, 0.01 to 10; E69's `bisect_factor`
from 1, E83's APL_STEP); a point at a bound of the bracket is reported (`at_bound`) and is not a choice. The APL->KC
weights are E89's `apl_kc_weights` at the family's measure F (exponent 2, or 1 in B2), the somatic threshold and the
KC-KC scale the fit's. The APL-off run of every point is the saved fit's (the path is absent without APL). The gain is
olfsysm's `bouton_inh_divisive` (a bouton's rate r divided by 1 + max(i, 0), i its APL->bouton weight times APL's
activity), set and reset by E92's `set_form` ('bdiv'). "E3's base" is the APL->bouton scale every *_factor_vs_e3 column
refers to (E10_joint_fit's `bouton_base` `initial`, 0.1343); a wiring's fitted multiple of it is its fit's
`aplpn_factor_vs_e3` (2.43 to 2.80 over the eight), and a point's multiple is its dose times that (`aplpn_multiple`).
Variants are named 'B_<form>_k<exponent>_d<dose>' ('B_sub_k2_d0.25', 'B_gain_k2_d0.05'), a control '<setting>_ctrl';
the references are 'B_sub_k2_d1' (the unchanged fits, every wiring), 'B_sub_k1_d1' (E89's saved C1 point at exponent 1,
S) and, for the record and not a choice, 'B_gain_k2_d1' (E92's B2 point, the gain at the fitted weight, S).

Parts (CLI `--part`):
- 0, the identities, before any scored point (identity.json; the parts B and C refuse to run until it passes): on each
  wiring (default all eight) dose 1 as a subtraction, through the dose's path, reproduces the saved APL-on and APL-off
  runs with 0 count differences (and olfsysm's wPNKC is the wiring's claw x bouton matrix); on the real wiring dose 0
  reproduces a simulator built without the path (E83's `build_simulator` at the fit's parameters with wAPLPN_scale 0)
  exactly, and the gain at dose 1 at exponent 2 reproduces E92's B2 point (B_F_k2_bdiv_real: its re-held factor and
  its counts and responses), with the form, the dose, the weights, APL's saved activity and the static inhibition reset
  after. Each wiring's claw x bouton matrix and unscaled APL->bouton weights are cached under outputs/E94_bouton_share/W/
  for part A.
- A (no scored simulations; the wirings' matrices from part 0's cache, else the bouton base's one setup simulation):
  A1 the path in the fits, per wiring, from E90's archived stage runs (APL on and off, window means over the 0.5 s pulse)
  and the simulator's own claw x bouton weights (a1_path.csv, a1_by_odor.csv): APL's subtraction at a bouton (10th
  percentile, median, 90th), the share of the boutons' summed rate it removes, the bouton-by-odor entries at zero, the
  boutons' spontaneous rate APL off / on, the subtraction's variation over the odors at the median bouton, APL's activity
  per odor against its mean, the share of the claws' summed drive APL removes, and the share of that loss made at the
  boutons (`claw_loss_shares`), overall, by odor and by a KC's number of claws (1 to 4, 5 to 8, 9 or more). A2 the
  connectome's share (a2_connectome.csv): APL's synapses onto PN boutons over those onto boutons and KCs in the calyx
  (roi_pre and roi_post CA(R)), from the tables the model reads (pyolfsysm's from_pratyush/2025-12-05 APLunit-to-PNbouton
  and APLunit-to-KCclaw parquets): all, those assigned to a bouton / a claw, and KCs with a type, beside Yang's own
  Claude's 10 to 11 %. A3 the bouton stage under other doses (a3_bouton_stage.csv, a3_summary.csv): the saved
  subtraction (APL off minus APL on, window means) times the dose, applied to the saved APL-off bouton rates and
  clipped at 0, and the gain (rate / (1 + dose x scaled APL->bouton weight x APL's window activity)), at APL's activity
  as saved (which the re-held runs do not keep); the cosine and Pearson distances between the odors' bouton patterns, on
  rates and on evoked rates (minus the pre-odor window's, transformed alike), by group (E92 part A's `group_masks`: the
  ester group, the alcohol group, an ester with an alcohol, the other megamat pairs, the six near pairs, the near ester
  pairs; and the 120). A4 the archive (a4_archive.csv): E3's sweep of the APL->bouton scale on the plain bouton base
  (x0, 0.5, 1, 2), E25's unbounded 70 Hz structure fit against E31's bounded and fixed ones, and C2's real fit, on raw
  counts and under the calibrated readout (K re-calibrated per run): the pattern part with its parts, the loss, the mean
  APL effect, ρ and Pearson r on the per-pair APL effect, the two groups' errors and effects, realism.
- B, `--arm B1|B2|B3 --wirings S|T|all8|<list>` (resumable; runs_<arm>.csv; runs under outputs/E94_bouton_share/B/):
  B1 the subtraction at exponent 2, doses B1_DOSES; B2 the subtraction at exponent 1, doses B2_DOSES; B3 the gain at
  exponent 2, doses B3_DOSES. Every run is saved with its stage matrices with APL on (E90's `read_stages`, E90's window;
  a superset of the entry's "doses 0.1, 0.25 and 0.5 and any setting the rules use"), each KC's summed claw drive from
  the boutons (<tag>_drive.parquet, the control's input) and E92's loss of drive at the claws (<tag>_loss.parquet).
  Each row: the absolute parameters (fixed_thr, wAPLKC_scale, wKCAPL_scale, wAPLPN_scale, wPNAPL_scale, wKCKC_scale), the
  dose, the form, the exponent, the tied factor, the APL->bouton scale as a multiple of E3's base and of E34's bound (and
  'above_bound'), the realism (E83's `realism`: the four bands under the rule of 2026-11-10 16:10, `realism_flag`, the
  earlier rule beside, `realism_flag_zcrit`, the largest spike counts), the share of the boutons' summed rate APL removes,
  the bouton-by-odor entries at zero, the share of the claws' loss made at the boutons (overall, by a KC's number of claws,
  by odor as JSON) beside the connectome's, E92's `gain_stats` (for the gain, the divisor 1 + i at the boutons feeding
  the responders' claws: 10th percentile, median, 90th, largest, the share above DIVISOR_HIGH), the APL->KC factors'
  spread (E89's `factor_summary`), `H.code_commits()` of the process (which must be clean). Provenance per part
  (provenance_<part>.json, `H.model_provenance` with the lab's input).
- hostcheck (`--point`, `--wiring`, default B_sub_k2_d0.5 on the real wiring): another host runs one point of viki's into
  outputs/.../hostcheck_<host>/ (runs_hostcheck_<host>.csv) and compares it with every other copy found (local B/, the
  archive's B/, its from_hal/B/, other hosts' hostcheck folders): hostcheck_<host>.json, the count differences and the
  factor's. Run it again to compare once the other copy exists.
- score (`--no-co`, `--no-zcell`): readings.csv (per variant and unit, each wiring and S / T where all four have a run;
  readouts 'cal' the calibrated readout with K set on each run's APL-on counts, the reading the choice and the gate are
  on; 'calK' the fit's K; 'raw' raw counts; 'zcell' E73's per-KC reading, E93's `zcell_distances`): E92's `score_row92`
  against every target of E89's `targets89` (halves A to D, all recordings, the six typical, the noise-free effect,
  baseline-noise units): loss, pattern part with its amplitude and shape parts (`H.apl_effect_spread_shape` with the
  target through `H.canonical_pairs`, checked to add up), spread ratio, Pearson r and ρ on the per-pair APL effect, the
  two groups' errors, the mean APL effect; pairs.csv (populations' per-pair distances and effect); line.csv (each S
  population's ester and alcohol APL-on distances against E92's part A line, `relative_distances` / `place_on_line`);
  coresponse.csv (E89's co-response readings, as E92's `coresponse92`); then part C's tables (choice.csv, gate.csv,
  equal.csv, verdict.json), comparison_table.md and the figure.
- C: the choice per half (`choose94`: among the references and B1 to B3, inside the bands on all four of S's wirings and
  not at a bracket's bound, the lowest pattern part of S's population under the readout against that half; nothing to
  confirm unless it is a B setting and gains CHOICE_MIN_GAIN over the unchanged fits); each chosen setting with something
  to confirm run on T (runs_C.csv), then its control on S and T (`control_matrix94`); then `score`, whose gate is E83's
  as E91 applied it (E92's `gate_direction92` and `gate92`), and the equal rule (`equal_rule94`: three kinds).
- D: the family fitted again (`E10_joint_fit.py -k 2 --aplpn-hold F --rate-target ctrl_frac --loss apl --apl-readout
  both --tag _d1`, D_MAXFEV evaluations), D1 at the anatomical dose (`anatomical_dose`: of B1's 0.1, 0.25 and 0.5 the one
  whose share of the claws' loss made at the boutons, mean over S's four runs, is nearest ANATOMICAL_SHARE), F the dose
  times the wiring's fitted multiple, from the wiring's B1 point at that dose (`--init-params` on d_start_params.json);
  `--pop S` (allowed whatever (C) says), `--pop T` only on an equal or better of S's (`score_d`); `--wirings` runs a
  subset (a lane); `--part Dscore` scores them (d_scores.csv, d_scores.json). D2 (`--part D2`) is a documented stub:
  `d2_plan` says whether it is due, `d2_command` the fit for B1 / B2 (`--aplpn-max 12`); B3's needs a new E10 option
  (`--bouton-divisive`), not written.
- E: only on a pass by (C) or D2 (a stub).
- steps (`--variants`, default every chosen or equal setting of verdict.json): the step from the claws to the somata in
  wild-type flies on S (E90's practice, as E92's part F4: E90's `read_wiring` of the saved stage matrices, E92 part A's
  somata and steps; steps.csv, somata_steps.csv, claws_steps.csv).
- variants: variants.csv and model_pair_distances_long.csv for the archive (`archive_tools/archive_model_outputs.py
  E94_bouton_share`, then `archive_tools/migrate_model_output_params.py --apply`).
- table: comparison_table.md and the figure from what exists.

Deviations from the entry, forced by the code:
- D1's command uses `--aplpn-hold F`, a new E10_joint_fit.py option (off by default), not `--aplpn-fixed F`: the latter
  (E31's) also fixes the bouton->APL scale at F times E3's base and ties only KC->APL to APL->KC, whereas C2's family
  ties bouton->APL to APL->KC, which the entry keeps ("a dose m multiplies the wiring's fitted APL->bouton scale and
  nothing else"). `--aplpn-hold` holds the APL->bouton scale only, with the family's tie.
- The control: olfsysm's `ext_inh_kc` must not be negative (checked inside its parallel region, which would abort the
  process). The entry's control (the reference's bouton path, dose 1, with a static inhibition equal to each KC's mean
  change of drive under the chosen setting: the reference's summed claw drive from the boutons minus the chosen
  setting's, mean over the odors) is run whenever that change is not negative for any KC (a chosen dose above 1). Where
  it is (a dose below 1, or the gain, which takes more from some KCs and less from others), the control is the bouton
  path removed (dose 0) with a static inhibition equal to each KC's mean loss of drive at the boutons under the chosen
  setting (APL-off drive minus the chosen's, never negative): the same question (the boutons' odor-specific action
  against a threshold raise by claw count with the same mean). Each control row says which (`control_form`).
- The choice's references (the unchanged fits, E89's point at exponent 1) are never "something to confirm": E89's point
  was read by E89's own gate. E92's gain point is read beside and is not a choice.
- The stage matrices are saved for every run, not only for the doses the entry names.

Outputs: results/E94_bouton_share/; spike counts, responses and stage matrices under outputs/E94_bouton_share/ (B/,
C's runs in B/ too, hostcheck_<host>/, W/), D's fits where E10_joint_fit.py writes them for this experiment.

Usage (from the repository root; one 16-thread process at a time on viki; check `uptime` first):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E94_bouton_share --part 0
    ./venv/bin/python -u -m experiments.E94_bouton_share --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E94_bouton_share --part B --arm B1
    ./venv/bin/python -u -m experiments.E94_bouton_share --part score
"""
from __future__ import annotations

import argparse
import json
import os
import re
import socket
import subprocess
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83
from experiments import E89_apl_breadth_weighting as E89
from experiments import E91_summed_drive_inhibition as E91
from experiments import E92_gain as E92G

EXPERIMENT = 'E94_bouton_share'
ARM = 'B'
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# the settings (the entry's Design; fixed before any run)

FORMS = {'sub': False, 'gain': True}                # olfsysm's bouton_inh_divisive
FORM92 = {'sub': 'sub', 'gain': 'bdiv'}             # E92's form names
B1_DOSES = (0.0, 0.1, 0.25, 0.5, 2.0, 4.0)          # the subtraction at exponent 2 (dose 1: the saved fits)
B2_DOSES = (2.0, 4.0)                               # the subtraction at exponent 1 (dose 1: E89's saved point)
B3_DOSES = (0.02, 0.05, 0.15)                       # the gain at exponent 2 (dose 1: E92's B2 point)
ARMS = {'B1': ('sub', 2.0, B1_DOSES), 'B2': ('sub', 1.0, B2_DOSES), 'B3': ('gain', 2.0, B3_DOSES)}
APL_BRACKET = E92G.APL_BRACKET                      # (0.01, 10.0)
APL_STEP = E83.APL_STEP
RAW_CRIT = E92G.RAW_CRIT                            # a model response is two or more spikes
BOUND_E34 = 3.0                                     # the APL->bouton bound of E34, in multiples of E3's base
CONNECTOME_SHARE = (0.10, 0.11)                     # Yang's own Claude's tables (hemibrain, male CNS)
ANATOMICAL_SHARE = 0.105                            # D1's target share of the claws' loss made at the boutons
ANATOMICAL_DOSES = (0.1, 0.25, 0.5)
ZERO_SHARE_NOTE = 0.2                               # a chosen setting with more bouton-by-odor entries at zero is named
ZERO_EPS = 1e-9                                     # a bouton's window-mean rate at or below this is at zero
DIVISOR_HIGH = E92G.DIVISOR_HIGH                    # 10: E92's rule for a shunt
CLAW_GROUPS = (('1-4', 1, 4), ('5-8', 5, 8), ('9+', 9, 10 ** 9))
A3_SUB_DOSES = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
A3_GAIN_DOSES = (0.02, 0.05, 0.15, 1.0)

# part C (the entry's rules: E91's)
CHOICE_MIN_GAIN = E91.CHOICE_MIN_GAIN               # 0.0015
CONTROL_RATIO = E91.CONTROL_RATIO
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
EQUAL_KINDS = ('lower_share', 'weaker_weighting', 'gain')
LOWER_SHARE_MAX_DOSE = 0.5                          # an equal "at a lower bouton share": B1 at dose 0.5 or less
HALVES = E91.HALVES
CHOICE_SCREENS = ('ref', 'B1', 'B2', 'B3')

# part D
D_MAXFEV = 150
D1_TAG = '_d1'
D2_TAG = '_d2'
D2_APLPN_MAX = 12.0
D_START_FILE = 'd_start_params.json'
UNCHANGED_S = E89.UNCHANGED_S                       # {'A': 0.0569, 'B': 0.0552, 'all': 0.0500}
UNCHANGED_8 = 0.0482
D_EQUAL = 0.003                                     # D1: S's fitted population within this of the unchanged fits'
D_S_GAIN = 0.0015                                   # D2: T's fits only if S's fitted population is this far below
D_PASS_GAIN, D_WIRINGS_LOWER, D_RAW_REF, D_RAW_SLACK, D_MEAN_SLACK = 0.003, 6, 0.073, 0.005, 0.01
FLIES_MEAN = E89.FLIES_MEAN
STEPS_NOISE_SEED = 400                              # E90's NOISE_SEED + this + the run's index (the claws' noisy reading)

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


def _flagged(x) -> bool:
    return isinstance(x, str) and x.strip() not in ('', 'nan')


def _as_bool(x) -> bool:
    """A CSV cell as a bool: True only for True / 'True' / 1 (NaN and None are False)."""
    if isinstance(x, str):
        return x.strip().lower() in ('true', '1')
    try:
        return bool(x) and not (isinstance(x, float) and np.isnan(x))
    except (TypeError, ValueError):
        return False


def _commits(strict: bool = True) -> dict:
    """H.code_commits() at the start of a simulation part, which must be clean (E89's `_commits_clean`)."""
    return E89._commits_clean() if strict else dict(H.code_commits())


# ---------------------------------------------------------------------------
# names

def variant_name(form: str, k: float, dose: float, control: bool = False) -> str:
    """'B_<form>_k<k>_d<dose>' ('B_sub_k2_d0.25'), '_ctrl' for the control of that setting."""
    if form not in FORMS:
        raise ValueError(f'unknown form {form!r}')
    if not (np.isfinite(dose) and dose >= 0 and np.isfinite(k) and k >= 0):
        raise ValueError(f'dose {dose}, exponent {k}')
    return f'{ARM}_{form}_k{float(k):g}_d{float(dose):g}' + ('_ctrl' if control else '')


_VARIANT_RE = re.compile(r'^B_(sub|gain)_k([0-9.]+)_d([0-9.]+)(_ctrl)?$')


def parse_variant(v: str) -> dict:
    """`variant_name`'s fields back: form, exponent, dose, control, control_of (the setting a control belongs to),
    screen ('ref' the unchanged fits or E89's point at exponent 1; 'E92' E92's gain point; 'B1' / 'B2' / 'B3'; 'control'),
    weighting (('F', exponent)), the APL->bouton divisive flag, and the kind of equal it can be (`equal_kind`)."""
    m = _VARIANT_RE.match(v)
    if not m:
        raise ValueError(f'variant {v!r}: not one of this experiment\'s')
    form, k, dose, ctl = m.group(1), float(m.group(2)), float(m.group(3)), bool(m.group(4))
    if ctl:
        screen = 'control'
    elif form == 'sub' and dose == 1.0 and k in (2.0, 1.0):
        screen = 'ref'
    elif form == 'gain' and dose == 1.0 and k == 2.0:
        screen = 'E92'
    elif form == 'sub' and k == 2.0:
        screen = 'B1'
    elif form == 'sub' and k == 1.0:
        screen = 'B2'
    elif form == 'gain' and k == 2.0:
        screen = 'B3'
    else:
        screen = 'other'
    out = dict(variant=v, form=form, exponent=k, dose=dose, control=ctl, screen=screen, weighting=('F', k),
               wtag=E89.wtag(('F', k)), bouton_inh_divisive=FORMS[form],
               control_of=variant_name(form, k, dose) if ctl else '')
    out['equal_kind'] = equal_kind(out)
    return out


def equal_kind(pv: dict) -> Optional[str]:
    """The kind of equal a setting can be (the entry's three): 'lower_share' B1 at a dose of LOWER_SHARE_MAX_DOSE or less,
    'weaker_weighting' B2, 'gain' B3; None otherwise (an equal above dose 1 at exponent 2 is not one: less defensible and
    gains nothing)."""
    if pv.get('control'):
        return None
    if pv['screen'] == 'B1':
        return 'lower_share' if pv['dose'] <= LOWER_SHARE_MAX_DOSE else None
    return {'B2': 'weaker_weighting', 'B3': 'gain'}.get(pv['screen'])


def ref_variant() -> str:
    """The unchanged fits: the subtraction at exponent 2 at dose 1."""
    return variant_name('sub', 2.0, 1.0)


def ref_k1_variant() -> str:
    """E89's saved point at exponent 1 (dose 1)."""
    return variant_name('sub', 1.0, 1.0)


def e92_variant() -> str:
    """E92's B2 point: the gain at the fitted weight (dose 1), exponent 2."""
    return variant_name('gain', 2.0, 1.0)


REF_VARIANTS = (ref_variant(), ref_k1_variant())


def arm_variants(arm: str) -> List[str]:
    form, k, doses = ARMS[arm]
    return [variant_name(form, k, d) for d in doses]


def run_tag(variant: str, wiring: str) -> str:
    return f'{variant}_{wiring}'


def _wirings(arg: str) -> Tuple[str, ...]:
    return E91._wirings(arg)


def aplpn_fit_multiple(wiring: str) -> float:
    """The wiring's fitted APL->bouton scale as a multiple of E3's base (its fits.csv row's aplpn_factor_vs_e3)."""
    key = ('multiple', wiring)
    if key not in _CACHE:
        folder, tag = E83.members()[wiring]
        r = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[tag]
        _CACHE[key] = float(r['aplpn_factor_vs_e3'])
    return _CACHE[key]


def aplpn_multiple(wiring: str, dose: float) -> float:
    """A point's APL->bouton scale as a multiple of E3's base: its dose times the wiring's fitted multiple."""
    return float(dose) * aplpn_fit_multiple(wiring)


# ---------------------------------------------------------------------------
# the simulator: a dose on the APL->bouton scale, the form at the boutons

class WiringModel94(E91.WiringModel91):
    """E91's `WiringModel91` (one wiring at its fit with the lab's input; olfsysm's `mp`, `rv` and `order` at hand) with a
    dose on the APL->bouton scale: `run_on(x)` runs APL on with the tied scales (E83's TIED_ON) times x and
    wAPLPN_scale times `dose` (1, the fit's, unless set), every other parameter the fit's. The APL-off run (`p_off`) is
    the saved fit's."""

    def __init__(self, label: str, check: bool = True, **model_kws):
        super().__init__(label, check=check, **model_kws)
        self.dose = 1.0

    def params_on(self, x: float = 1.0, dose: Optional[float] = None) -> Dict[str, float]:
        """The APL-on parameters at tied factor x and the dose (default the model's)."""
        p = dict(self.p)
        for k in E83.TIED_ON:
            p[k] = self.p[k] * x
        p['wAPLPN_scale'] = self.p['wAPLPN_scale'] * float(self.dose if dose is None else dose)
        return p

    def run_on(self, x: float = 1.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self.sim(self.params_on(x))


def set_form94(mp, form: str = 'sub') -> None:
    """The form of APL's inhibition at the boutons on a simulator's olfsysm parameters ('sub' the subtraction, 'gain'
    olfsysm's `bouton_inh_divisive`; E92's `set_form`, the claws subtractive in both)."""
    if form not in FORMS:
        raise ValueError(f'unknown form {form!r}')
    E92G.set_form(mp, FORM92[form])


def get_form94(mp) -> str:
    f = E92G.get_form(mp)
    inv = {v: k for k, v in FORM92.items()}
    if f not in inv:
        raise RuntimeError(f'the simulator holds form {f!r}, not one of this experiment\'s')
    return inv[f]


def b_reading94(model: WiringModel94, wk: Optional[pd.Series], dose: float, form: str, ext: Optional[np.ndarray] = None,
                label: str = '', bracket: Tuple[float, float] = APL_BRACKET) -> dict:
    """One setting on a wiring at its fit: the APL->KC weights `wk` (None: the fit's), the dose on the APL->bouton scale,
    the form at the boutons, a static inhibition at the soma `ext` (per KC in the outputs' row order, while the odor is on;
    None for none); the tied APL scale bisected to the fit's APL-on count (E69's `bisect_factor` in `bracket`, from 1);
    then the run again with APL's activity saved (its counts must be equal), read for E90's stage matrices
    (`read_stages`, E90's window) and E92's `gain_stats`. The APL-off run is the saved fit's. The weights, dose, form,
    static inhibition and `save_inh_sims` are reset on return, also on an error."""
    from experiments import E69_kc_firing_ceiling as E69
    from experiments import E90_model_stages as E90
    t0 = time.time()
    E89.set_weights(model, wk)
    model.dose = float(dose)
    set_form94(model.mp, form)
    E91.set_ext_inh(model.mp, model.rv, ext, len(model.order))
    try:
        x, on, hist, at_bound = E69.bisect_factor(model.run_on, model.target_on, model.size, bracket, APL_STEP,
                                                  f'{label} APL scale')
        params = model.params_on(x)
        model.mp.kc.save_inh_sims = True
        try:
            on2 = model.run_on(x)
        finally:
            model.mp.kc.save_inh_sims = False
        ns = int((on2[1].to_numpy(dtype=float) != on[1].to_numpy(dtype=float)).sum())
        if ns:
            raise RuntimeError(f"{model.label}: the run again with APL's activity saved differs in {ns} counts")
        odors = [H.odor_abbrev(c) for c in on[1].columns]
        st = E90.read_stages(model.mp, model.rv, model.W, model.order, odors, E90.MODEL_WINDOW)
        gs = E92G.gain_stats(model.rv, model.mp, model.order, on[1].to_numpy(dtype=float) > RAW_CRIT, FORM92[form])
    finally:
        E91.set_ext_inh(model.mp, model.rv, None)
        set_form94(model.mp, 'sub')
        model.dose = 1.0
        E89.set_weights(model, None)
    return dict(on=on, off=model.saved['off'], apl_factor=float(x), at_bound_apl=bool(at_bound), n_sims=len(hist) + 1,
                params=params, seconds=time.time() - t0, history=hist, stages=st, gain=gs, bracket=tuple(bracket),
                dose=float(dose), form=form)


# ---------------------------------------------------------------------------
# the wirings' matrices and E90's stage runs of the fits

def _w_dir() -> Path:
    return local_dir('W')


def save_wiring(label: str, W: pd.DataFrame, wAPLPN: pd.Series) -> None:
    """A wiring's claw x bouton matrix (sparse, long: claw, bouton, weight; the index and columns beside) and its unscaled
    APL->bouton weights under outputs/E94_bouton_share/W/."""
    d = _w_dir()
    Wm = W.to_numpy(dtype=float)
    r, c = np.nonzero(Wm)
    pd.DataFrame(dict(claw=r, bouton=c, weight=Wm[r, c])).to_parquet(d / f'{label}_W_long.parquet')
    W.index.to_frame(index=False).to_parquet(d / f'{label}_W_index.parquet')
    W.columns.to_frame(index=False).to_parquet(d / f'{label}_W_columns.parquet')
    pd.DataFrame(dict(wAPLPN=np.asarray(wAPLPN, dtype=float).ravel())).to_parquet(d / f'{label}_wAPLPN.parquet')


def load_wiring(label: str) -> Optional[Tuple[pd.DataFrame, np.ndarray]]:
    d = _w_dir()
    ps = [d / f'{label}_{k}.parquet' for k in ('W_long', 'W_index', 'W_columns', 'wAPLPN')]
    if not all(p.exists() for p in ps):
        return None
    lg, ix, cx, wp = (pd.read_parquet(p) for p in ps)
    idx = pd.MultiIndex.from_frame(ix) if ix.shape[1] > 1 else pd.Index(ix.iloc[:, 0])
    cols = pd.MultiIndex.from_frame(cx) if cx.shape[1] > 1 else pd.Index(cx.iloc[:, 0])
    Wm = np.zeros((len(idx), len(cols)))
    Wm[lg.claw.to_numpy(), lg.bouton.to_numpy()] = lg.weight.to_numpy()
    return pd.DataFrame(Wm, index=idx, columns=cols), wp.wAPLPN.to_numpy(dtype=float)


def wiring_W(label: str) -> Tuple[pd.DataFrame, np.ndarray]:
    """(claw x bouton matrix, unscaled APL->bouton weights) of a wiring: part 0's cache, else E83's `wiring_parts` (the
    bouton base's one setup simulation), cached then."""
    got = load_wiring(label)
    if got is None:
        P = E83.wiring_parts(label)
        save_wiring(label, P['W'], P['wnorm']['wAPLPN'])
        got = load_wiring(label)
    return got


def fit_stage(wiring: str, cond: str) -> dict:
    """E90's archived stage run of a wiring's fit ('on' / 'off'), cached; raises if absent."""
    from experiments import E90_model_stages as E90
    key = ('stage', wiring, cond)
    if key not in _CACHE:
        st = E90.load_stage('B', f"{wiring}_{'on' if cond == 'on' else 'acute_off'}")
        if st is None:
            raise RuntimeError(f"E90's stage run of {wiring} ({cond}) not found")
        _CACHE[key] = st
    return _CACHE[key]


def _odors(df: pd.DataFrame) -> List[str]:
    return [str(c) for c in df.columns]


def kc_of_claw(W: pd.DataFrame) -> np.ndarray:
    return np.asarray(W.index.get_level_values(H.KC_ID_LEVEL))


def claw_counts(kc_of: np.ndarray) -> np.ndarray:
    """Per claw, its KC's number of claws."""
    s = pd.Series(np.asarray(kc_of))
    return s.map(s.value_counts()).to_numpy()


def check_stage_shapes(W: pd.DataFrame, st: dict, label: str = '') -> None:
    """The stage matrices are the wiring's: boutons (and their glomeruli) and claws (and their KCs) in its order."""
    b, c = st['bouton_win'], st['claw_win']
    if b.shape[0] != W.shape[1] or c.shape[0] != W.shape[0]:
        raise RuntimeError(f'{label}: stage {b.shape[0]} boutons / {c.shape[0]} claws, wiring {W.shape[1]} / {W.shape[0]}')
    gw = [str(x) for x in W.columns.get_level_values('glomerulus')]
    gs = [str(x) for x in b.index.get_level_values('glomerulus')]
    kw = np.asarray(W.index.get_level_values(H.KC_ID_LEVEL)).astype(np.int64)
    ks = np.asarray(c.index.get_level_values(H.KC_ID_LEVEL)).astype(np.int64)
    if gw != gs or not np.array_equal(kw, ks):
        raise RuntimeError(f'{label}: the stage matrices are not in the wiring\'s order')


def claw_loss_shares(W, b_off: np.ndarray, b_on: np.ndarray, claw_on: np.ndarray,
                     kc_of: Optional[np.ndarray] = None) -> dict:
    """The claws' loss of drive under APL and the share of it made at the boutons, from window means (claws x odors /
    boutons x odors): x = W @ b is a claw's PN drive (linear in the bouton rates, so the window mean of the drive is W
    times the window-mean rates); with APL off it is the claw's whole drive (no claw inhibition), with APL on the claws'
    states `claw_on` are what is left after APL's inhibition at the claws (E90's `claw_inputs`). The loss is
    x_off - claw_on, of which x_off - x_on is made at the boutons and x_on - claw_on at the claws. Returns
    'drive_removed' (the share of the claws' summed drive APL removes in all), 'bouton_share' (the share of that loss
    made at the boutons), 'claw_part' (the rest), 'by_odor' (bouton_share per odor), and, with `kc_of` (each claw's KC),
    'bouton_share_<group>' by a KC's number of claws (CLAW_GROUPS) and 'drive_removed_<group>'."""
    Wm = np.asarray(W, dtype=float)
    x_off = Wm @ np.asarray(b_off, dtype=float)
    x_on = Wm @ np.asarray(b_on, dtype=float)
    c_on = np.asarray(claw_on, dtype=float)
    if c_on.shape != x_off.shape:
        raise ValueError(f'claw states {c_on.shape}, drives {x_off.shape}')
    tot, at_b = x_off - c_on, x_off - x_on

    def share(a, b):
        b = float(b)
        return float(a) / b if b > 0 else np.nan
    out = dict(drive_removed=share(tot.sum(), x_off.sum()), bouton_share=share(at_b.sum(), tot.sum()),
               claw_part=share((x_on - c_on).sum(), tot.sum()),
               by_odor=[share(a, b) for a, b in zip(at_b.sum(axis=0), tot.sum(axis=0))])
    if kc_of is not None:
        n = claw_counts(kc_of)
        for name, lo, hi in CLAW_GROUPS:
            m = (n >= lo) & (n <= hi)
            out[f'bouton_share_{name}'] = share(at_b[m].sum(), tot[m].sum())
            out[f'drive_removed_{name}'] = share(tot[m].sum(), x_off[m].sum())
            out[f'n_claws_{name}'] = int(m.sum())
    return out


def bouton_realism(b_off: np.ndarray, b_on: np.ndarray) -> dict:
    """APL's action on the boutons from window means (boutons x odors): the share of the boutons' summed rate it removes,
    the share of bouton-by-odor entries at zero with APL on (and off), the removed rate per entry (10th percentile, median,
    90th; APL off minus APL on)."""
    off, on = np.asarray(b_off, dtype=float), np.asarray(b_on, dtype=float)
    sub = off - on
    q = np.quantile(sub, [0.1, 0.5, 0.9])
    return dict(bouton_share_removed=float(1.0 - on.sum() / off.sum()) if off.sum() > 0 else np.nan,
                bouton_zero_share=float((on <= ZERO_EPS).mean()), bouton_zero_share_off=float((off <= ZERO_EPS).mean()),
                bouton_sub_p10=float(q[0]), bouton_sub_p50=float(q[1]), bouton_sub_p90=float(q[2]))


def kc_drive(W: pd.DataFrame, b: np.ndarray, kc_index: pd.Index, odors: Sequence[str]) -> pd.DataFrame:
    """Each KC's summed claw drive from the boutons (sum over its claws of W @ b; window means), KCs in `kc_index`'s order
    (the outputs' rows, matched on the KC id level), odors as columns."""
    x = np.asarray(W, dtype=float) @ np.asarray(b, dtype=float)
    d = pd.DataFrame(x, columns=list(odors)).groupby(kc_of_claw(W)).sum()
    ids = kc_index.get_level_values(H.KC_ID_LEVEL) if isinstance(kc_index, pd.MultiIndex) else kc_index
    out = d.reindex(np.asarray(ids))
    if out.isna().any().any():
        raise RuntimeError(f'{int(out.isna().any(axis=1).sum())} KCs of the outputs have no claws in the wiring')
    out.index = kc_index
    return out


def fit_drive(wiring: str, cond: str, kc_index: pd.Index) -> pd.DataFrame:
    """`kc_drive` of the wiring's fit ('on' / 'off'; E90's stage runs)."""
    W, _ = wiring_W(wiring)
    st = fit_stage(wiring, cond)
    check_stage_shapes(W, st, f'{wiring} {cond}')
    return kc_drive(W, st['bouton_win'].to_numpy(dtype=float), kc_index, _odors(st['bouton_win']))


def control_matrix94(d_ref: pd.DataFrame, d_chosen: pd.DataFrame, d_off: pd.DataFrame
                     ) -> Tuple[pd.DataFrame, str, float]:
    """The control's static inhibition at the soma (olfsysm `ext_inh_kc`, KCs x odors, the same for every odor) and the
    dose of the bouton path it runs with. The entry's ('ref_plus_change', dose 1): the reference's bouton path with each
    KC's mean over the odors of its change of drive under the chosen setting (the reference's summed claw drive from the
    boutons minus the chosen setting's), used when no KC's change is negative (beyond rounding); otherwise
    ('none_plus_removal', dose 0): the bouton path removed with each KC's mean loss of drive at the boutons under the
    chosen setting (APL-off drive minus the chosen setting's, never negative), olfsysm refusing a negative `ext_inh_kc`.
    The three frames are KCs x odors in the same order."""
    for d in (d_chosen, d_off):
        if not (d.index.equals(d_ref.index) and list(d.columns) == list(d_ref.columns)):
            raise ValueError('the drives are not in the same KC and odor order')
    ch = (d_ref - d_chosen).astype(float).mean(axis=1)
    tol = 1e-9 * max(1.0, float(np.abs(ch).max()) if len(ch) else 1.0)
    if (ch >= -tol).all():
        v, form, dose = ch.clip(lower=0.0), 'ref_plus_change', 1.0
    else:
        v, form, dose = (d_off - d_chosen).astype(float).mean(axis=1).clip(lower=0.0), 'none_plus_removal', 0.0
    m = pd.DataFrame(np.repeat(v.to_numpy()[:, None], d_ref.shape[1], axis=1), index=d_ref.index, columns=d_ref.columns)
    return m, form, dose


# ---------------------------------------------------------------------------
# saving and loading runs

STAGE_KEYS = ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')


def _dirs(sub: str = 'B') -> List[Path]:
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    ds = [H.outputs_dir(EXPERIMENT) / sub, arc / sub, arc / 'from_hal' / sub]
    return ds


def save_run(tag: str, o: Dict[str, tuple], sub: str = 'B', stages: Optional[dict] = None,
             drive: Optional[pd.DataFrame] = None, loss: Optional[pd.DataFrame] = None) -> None:
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
    if drive is not None:
        x = drive.copy()
        x.columns = [str(c) for c in x.columns]
        x.to_parquet(d / f'{tag}_drive.parquet')
    if loss is not None:
        loss.to_parquet(d / f'{tag}_loss.parquet')


def load_run(tag: str, sub: str = 'B') -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}: (responses, spike counts)): the local outputs, else the archive (and its from_hal/)."""
    for d in _dirs(sub):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def load_stage94(tag: str, sub: str = 'B') -> Optional[dict]:
    """A run's stage matrices with APL on (and its responses and counts), or None."""
    for d in _dirs(sub):
        if (d / f'{tag}_on_claw_win.parquet').exists():
            st = {k: pd.read_parquet(d / f'{tag}_on_{k}.parquet') for k in STAGE_KEYS}
            st['apl'] = pd.read_csv(d / f'{tag}_on_apl.csv', index_col=0)
            o = E83._load_pair(d, tag)
            if o is not None:
                st['responses'], st['spike_counts'] = o['on']
            return st
    return None


def load_drive(tag: str, sub: str = 'B') -> pd.DataFrame:
    for d in _dirs(sub):
        p = d / f'{tag}_drive.parquet'
        if p.exists():
            return pd.read_parquet(p)
    raise RuntimeError(f'no drive saved for {tag}')


def load_any(source: str, tag: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: 'E94' this experiment's, 'saved' the wiring's fit, 'E89' E89's C runs, 'E92' E92's B runs."""
    if source == 'saved':
        return E83.saved_fit(wiring)['saved']
    o = {'E94': load_run, 'E89': E89.load_run, 'E92': E92G.load_run}[source](tag)
    if o is None:
        raise RuntimeError(f'no saved run {source}:{tag}')
    return o


def _append(path: Path, row: dict) -> None:
    E83._append_row(path, row)


def _done(path: Path) -> set:
    return set(pd.read_csv(path).tag) if path.exists() else set()


def _frame(a: np.ndarray, like: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(a, index=like.index, columns=[str(c) for c in like.columns])


# ---------------------------------------------------------------------------
# part 0: the identities

def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('passes'))


def _counts_diff(a: Tuple[pd.DataFrame, pd.DataFrame], b: Tuple[pd.DataFrame, pd.DataFrame]) -> dict:
    same = bool(a[1].index.equals(b[1].index) and list(map(str, a[1].columns)) == list(map(str, b[1].columns)))
    nd = int((a[1].to_numpy(dtype=float) != b[1].to_numpy(dtype=float)).sum()) if same else -1
    nr = int((a[0].to_numpy(dtype=float) != b[0].to_numpy(dtype=float)).sum()) if same else -1
    return dict(counts_differ=nd, responses_differ=nr, index_equal=same, exact=bool(same and nd == 0 and nr == 0))


def part_0(out: Path, wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """The identities (module docstring, part 0); identity.json (resumable per wiring); raises if one fails."""
    import gc
    cc = _commits()
    path = out / 'identity.json'
    res = _read_json(path)
    res.setdefault('commits', cc)
    for w in wirings:
        key = f'sub_d1_{w}'
        if res.get(key, {}).get('exact') and (w != 'real' or (res.get('sub_d0_real', {}).get('exact')
                                                                and res.get('gain_d1_real', {}).get('passes'))):
            continue
        t0 = time.time()
        m = WiringModel94(w, check=False)
        save_wiring(w, m.W, m.wAPLPN)
        wp = np.asarray(m.rv.kc.wPNKC, dtype=float)
        wiring_equal = bool(wp.shape == m.W.shape and np.array_equal(wp, m.W.to_numpy(dtype=float)))
        on, off = m.run_on(1.0), m.sim(m.p_off())
        d_on, d_off = _counts_diff(on, m.saved['on']), _counts_diff(off, m.saved['off'])
        res[key] = dict(wiring=w, on=d_on, off=d_off, wiring_equal=wiring_equal, params_on_equal_fit=bool(
            m.params_on(1.0) == dict(m.p)), exact=bool(d_on['exact'] and d_off['exact'] and wiring_equal),
            seconds=time.time() - t0)
        print(f"E94 0 {key}: {res[key]}", flush=True)
        _write_json(path, res)
        if w == 'real':
            t1 = time.time()
            p0 = dict(m.p, wAPLPN_scale=0.0)
            sim0, _ = E83.build_simulator('real', params=p0)
            a = sim0(p0)
            m.dose = 0.0
            try:
                b = m.run_on(1.0)
            finally:
                m.dose = 1.0
            d0 = _counts_diff(b, a)
            res['sub_d0_real'] = dict(d0, changed_from_fit=_counts_diff(b, m.saved['on'])['counts_differ'],
                                      seconds=time.time() - t1)
            print(f"E94 0 sub_d0_real: {res['sub_d0_real']}", flush=True)
            del sim0
            gc.collect()
            _write_json(path, res)
            t2 = time.time()
            r92 = pd.read_csv(H.results_dir(E92G.EXPERIMENT) / 'runs_B2.csv').set_index('tag').loc['B_F_k2_bdiv_real']
            o92 = E92G.load_run('B_F_k2_bdiv_real')
            if o92 is None:
                raise RuntimeError("E92's B2 point on the real wiring (B_F_k2_bdiv_real) is not saved")
            wk, _ = E91._setting_weights(m, ('F', 2.0))
            rg = b_reading94(m, wk, 1.0, 'gain', label='identity gain d1 real')
            dg = _counts_diff(rg['on'], o92['on'])
            f92 = float(r92['apl_factor'])
            reset = dict(form_after=get_form94(m.mp), dose_after=float(m.dose), weights_reset=bool(np.array_equal(m._kc, m._kc0)),
                         ext_inh_after=int(np.asarray(m.rv.kc.ext_inh).size), save_inh_sims_after=bool(m.mp.kc.save_inh_sims))
            res['gain_d1_real'] = dict(dg, apl_factor=rg['apl_factor'], apl_factor_e92=f92,
                                       factor_rel_diff=abs(rg['apl_factor'] - f92) / f92, n_sims=rg['n_sims'], **reset,
                                       divisor_bouton_p50=rg['gain']['row'].get('divisor_bouton_p50'),
                                       divisor_bouton_p50_e92=float(r92.get('divisor_bouton_p50', np.nan)),
                                       seconds=time.time() - t2)
            res['gain_d1_real']['passes'] = bool(dg['exact'] and res['gain_d1_real']['factor_rel_diff'] <= 1e-12
                                                 and reset['form_after'] == 'sub' and reset['dose_after'] == 1.0
                                                 and reset['weights_reset'] and reset['ext_inh_after'] == 0
                                                 and not reset['save_inh_sims_after'])
            print(f"E94 0 gain_d1_real: {res['gain_d1_real']}", flush=True)
            _write_json(path, res)
        del m
        gc.collect()
    keys = [f'sub_d1_{w}' for w in ALL_WIRINGS]
    have_all = all(k in res for k in keys) and 'sub_d0_real' in res and 'gain_d1_real' in res
    res['passes'] = bool(have_all and all(res[k]['exact'] for k in keys) and res['sub_d0_real']['exact']
                         and res['gain_d1_real']['passes'])
    res['complete'] = bool(have_all)
    _write_json(path, res)
    failed = [k for k in keys + ['sub_d0_real'] if k in res and not res[k]['exact']]
    if 'gain_d1_real' in res and not res['gain_d1_real']['passes']:
        failed.append('gain_d1_real')
    if failed:
        raise RuntimeError(f'E94 part 0: the identity fails for {failed} (identity.json); the experiment stops here')
    print(f"E94 0: {'passes' if res['passes'] else 'so far exact (not all wirings yet)'}", flush=True)
    return res


# ---------------------------------------------------------------------------
# part A (no scored simulations)

def a1_path(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A1 (module docstring): a1_path.csv (per wiring) and a1_by_odor.csv."""
    rows, odor_rows = [], []
    for w in ALL_WIRINGS:
        W, _ = wiring_W(w)
        on, off = fit_stage(w, 'on'), fit_stage(w, 'off')
        for st, c in ((on, 'on'), (off, 'off')):
            check_stage_shapes(W, st, f'{w} {c}')
        odors = _odors(on['bouton_win'])
        b_on, b_off = on['bouton_win'][odors].to_numpy(dtype=float), off['bouton_win'][odors].to_numpy(dtype=float)
        c_on, c_off = on['claw_win'][odors].to_numpy(dtype=float), off['claw_win'][odors].to_numpy(dtype=float)
        x_off = W.to_numpy(dtype=float) @ b_off
        sub = b_off - b_on
        mu = sub.mean(axis=1)
        cv = np.where(mu > 0, sub.std(axis=1) / np.where(mu > 0, mu, 1.0), np.nan)
        apl = on['apl'].reindex(odors)['win'].astype(float)
        sh = claw_loss_shares(W, b_off, b_on, c_on, kc_of_claw(W))
        row = dict(wiring=w, n_boutons=int(W.shape[1]), n_claws=int(W.shape[0]), n_kcs=int(len(np.unique(kc_of_claw(W)))),
                   aplpn_fit_multiple=aplpn_fit_multiple(w), **bouton_realism(b_off, b_on),
                   spont_off=float(off['bouton_pre'].to_numpy(dtype=float).mean()),
                   spont_on=float(on['bouton_pre'].to_numpy(dtype=float).mean()),
                   sub_cv_over_odors_median_bouton=float(np.nanmedian(cv)),
                   apl_rel_max=float((apl / apl.mean() - 1.0).abs().max()), apl_win_mean=float(apl.mean()),
                   claw_off_check_max_abs=float(np.abs(x_off - c_off).max()),
                   **{k: v for k, v in sh.items() if k != 'by_odor'})
        rows.append(row)
        odor_rows += [dict(wiring=w, odor=o, bouton_share=s, apl_win=float(apl[o])) for o, s in zip(odors, sh['by_odor'])]
    a1 = pd.DataFrame(rows)
    a1.to_csv(out / 'a1_path.csv', index=False)
    bo = pd.DataFrame(odor_rows)
    bo.to_csv(out / 'a1_by_odor.csv', index=False)
    return a1, bo


A2_TABLES = {'pn': 'APLunit-to-PNbouton_Connectivity_2025-12-05-13-39-22.parquet',
             'kc': 'APLunit-to-KCclaw_Connectivity_2025-12-14-21-12-48.parquet'}   # as pyolfsysm's connectome_APL_weights


def a2_connectome(out: Path) -> pd.DataFrame:
    """A2 (module docstring): a2_connectome.csv."""
    from pyolfsysm import mb_model as M
    d = Path(M.prat_hemibrain_seg_dir)
    pn = pd.read_parquet(d / A2_TABLES['pn'])
    kc = pd.read_parquet(d / A2_TABLES['kc'])
    pn = pn[(pn.roi_pre == 'CA(R)') & (pn.roi_post == 'CA(R)')]
    kc = kc[(kc.roi_pre == 'CA(R)') & (kc.roi_post == 'CA(R)')]
    sets = {'calyx, all': (pn, kc),
            'calyx, assigned to a bouton / a claw': (pn[pn.anatomical_bouton >= 0], kc[kc.anatomical_claw_corrected >= 0]),
            'calyx, KCs with a type': (pn, kc[kc.type_post.notna()]),
            'calyx, assigned, KCs with a type': (pn[pn.anatomical_bouton >= 0],
                                                 kc[(kc.anatomical_claw_corrected >= 0) & kc.type_post.notna()])}
    rows = [dict(set=k, n_apl_to_pn=int(len(a)), n_apl_to_kc=int(len(b)), share_pn=float(len(a) / (len(a) + len(b))),
                 tables=f"{A2_TABLES['pn']}; {A2_TABLES['kc']}") for k, (a, b) in sets.items()]
    rows.append(dict(set="Yang's own Claude's tables (hemibrain, male CNS)", share_pn_lo=CONNECTOME_SHARE[0],
                     share_pn_hi=CONNECTOME_SHARE[1], share_pn=float(np.mean(CONNECTOME_SHARE))))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a2_connectome.csv', index=False)
    return df


def _pairs120(D: np.ndarray, odors: Sequence[str]) -> pd.Series:
    from experiments import E92_steps_like_for_like as A
    return A._square_pairs(D, list(odors), A.PAIRS120)


A3_GROUPS = ('ester', 'alcohol', 'ester_alcohol', 'other', 'near', 'near_ester', 'all66', 'all120')


def a3_group_means(pairs: pd.Series) -> Dict[str, float]:
    from experiments import E92_steps_like_for_like as A
    m = A.group_masks(pairs.index)
    v = pairs.to_numpy(dtype=float)
    out = {g: float(np.nanmean(v[m[g]])) for g in A3_GROUPS if g in m and m[g].any()}
    out['all120'] = float(np.nanmean(v))
    return out


def a3_transform(O: np.ndarray, On: np.ndarray, form: str, dose: float, wa: Optional[np.ndarray] = None) -> np.ndarray:
    """A3's bouton rates at a dose: the subtraction (O - dose x (O - On), clipped at 0: the saved subtraction, as realized
    with APL on, times the dose) or the gain (O / (1 + dose x wa), wa the scaled APL->bouton weight times APL's activity,
    boutons x odors)."""
    O, On = np.asarray(O, dtype=float), np.asarray(On, dtype=float)
    if form == 'sub':
        return np.clip(O - float(dose) * (O - On), 0.0, None)
    if form == 'gain':
        return O / (1.0 + float(dose) * np.clip(np.asarray(wa, dtype=float), 0.0, None))
    raise ValueError(form)


def a3_bouton_stage(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """A3 (module docstring): a3_bouton_stage.csv and a3_summary.csv (the mean over the eight wirings)."""
    from experiments import E92_steps_like_for_like as A
    rows = []
    for w in ALL_WIRINGS:
        W, wapl = wiring_W(w)
        on, off = fit_stage(w, 'on'), fit_stage(w, 'off')
        odors = _odors(on['bouton_win'])
        O, On = off['bouton_win'][odors].to_numpy(dtype=float), on['bouton_win'][odors].to_numpy(dtype=float)
        PO, Pon = off['bouton_pre'][odors].to_numpy(dtype=float), on['bouton_pre'][odors].to_numpy(dtype=float)
        apl = on['apl'].reindex(odors)
        wsc = wapl * E83.saved_fit(w)['p']['wAPLPN_scale']
        wa, wa_pre = np.outer(wsc, apl['win'].to_numpy(dtype=float)), np.outer(wsc, apl['pre'].to_numpy(dtype=float))
        base = {}
        pts = [('sub', d) for d in A3_SUB_DOSES] + [('gain', d) for d in A3_GAIN_DOSES]
        for form, dose in pts:
            Y = a3_transform(O, On, form, dose, wa)
            PY = a3_transform(PO, Pon, form, dose, wa_pre)
            br = bouton_realism(O, Y)
            for ms in ('cosine', 'pearson'):
                for sig, X in (('rate', Y), ('evoked', Y - PY)):
                    gm = a3_group_means(_pairs120(A.distance_matrix(X, ms), odors))
                    for g, v in gm.items():
                        if form == 'sub' and dose == 0.0:
                            base[(ms, sig, g)] = v
                        rows.append(dict(wiring=w, form=form, dose=dose, measure=ms, signal=sig, group=g, distance=v,
                                         change=v - base[(ms, sig, g)], share_removed=br['bouton_share_removed'],
                                         zero_share=br['bouton_zero_share']))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a3_bouton_stage.csv', index=False)
    sm = df.groupby(['form', 'dose', 'measure', 'signal', 'group'], sort=False)[
        ['distance', 'change', 'share_removed', 'zero_share']].mean().reset_index()
    sm.to_csv(out / 'a3_summary.csv', index=False)
    return df, sm


A4_RUNS = (('E3_drive_dependent_gain', 'bouton_pnapl_x0', 0.0), ('E3_drive_dependent_gain', 'bouton_pnapl_x0.5', 0.5),
           ('E3_drive_dependent_gain', 'bouton_pnapl_x1', 1.0), ('E3_drive_dependent_gain', 'bouton_pnapl_x2', 2.0),
           ('E25_joint_fit_structure', 'thr70_k2_ext', None), ('E31_aplpn_bound', 'thr70_k2_aplpnmax3', None),
           ('E31_aplpn_bound', 'thr70_k2_aplpnfixed1', None), ('E31_aplpn_bound', 'thr70_k2_aplpntied', None),
           ('E31_aplpn_bound', 'thr70_k2_aplkcmin0.8', None))
A4_COLS = ('loss', 'pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson', 'rho', 'err_ester', 'err_alcohol')


def _a4_multiple(folder: str, tag: str, given: Optional[float]) -> float:
    if given is not None:
        return float(given)
    try:
        r = pd.read_csv(H.results_dir(folder) / 'fits.csv').set_index('tag').loc[tag]
        return float(r['aplpn_factor_vs_e3'])
    except (OSError, KeyError, ValueError):
        return np.nan


def group_effects(d: Dict[str, pd.Series]) -> Dict[str, float]:
    """The model's mean per-pair APL effect over E77's ester and alcohol groups."""
    from experiments import E77_input_side_common_error as E77
    me = E83._cs(H.apl_effect(d['on'], d['off']))
    m = E77.kind_masks(me.index)
    return dict(eff_ester=float(me[m['ester']].mean()), eff_alcohol=float(me[m['alcohol']].mean()))


def a4_archive(out: Path) -> pd.DataFrame:
    """A4 (module docstring): a4_archive.csv."""
    from experiments import E77_input_side_common_error as E77
    arc = H.shared_dir() / 'data' / 'model_outputs'
    tg = E89.targets89()
    runs = [(f, t, x) for f, t, x in A4_RUNS] + [('C2', 'real', 1.0)]
    rows = []
    for folder, tag, mult in runs:
        if folder == 'C2':
            o = E83.saved_fit('real')['saved']
            m = aplpn_fit_multiple('real')
        else:
            o = E83._load_pair(arc / folder, tag)
            m = _a4_multiple(folder, tag, mult)
        if o is None:
            rows.append(dict(folder=folder, tag=tag, note='not in the archive'))
            continue
        real = E83.realism(o)
        for ro in ('raw', 'cal'):
            d = E83.pair_distances(o, ro)
            sc = E92G.score_row92(d)
            row = dict(folder=folder, tag=tag, aplpn_multiple=m, readout=ro, mean_effect=sc['mean_effect'],
                       **group_effects(d))
            for t in ('z8', 'halfA', 'halfB'):
                row.update({f'{t}_{c}': sc[f'{t}_{c}'] for c in A4_COLS})
            de = E83._cs(tg['z8'])
            mk = E77.kind_masks(de.index)
            row.update(data_eff_ester=float(de[mk['ester']].mean()), data_eff_alcohol=float(de[mk['alcohol']].mean()))
            row.update({k: real[k] for k in ('rate_on', 'rate_off', 'ratio', 'rate_rho_on', 'opk_on_tv', 'opk_off_tv',
                                             'max_count_on', 'realism_flag')})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'a4_archive.csv', index=False)
    return df


def part_a(out: Path) -> dict:
    """Part A: A1 to A4 and a_summary.json (the ranges over the eight wirings)."""
    t0 = time.time()
    a1, bo = a1_path(out)
    a2 = a2_connectome(out)
    _, a3 = a3_bouton_stage(out)
    a4 = a4_archive(out)
    cols = ['bouton_sub_p10', 'bouton_sub_p50', 'bouton_sub_p90', 'bouton_share_removed', 'bouton_zero_share', 'spont_off',
            'spont_on', 'sub_cv_over_odors_median_bouton', 'apl_rel_max', 'drive_removed', 'bouton_share',
            'bouton_share_1-4', 'bouton_share_5-8', 'bouton_share_9+', 'aplpn_fit_multiple']
    summ = dict(a1={c: dict(min=float(a1[c].min()), max=float(a1[c].max()), mean=float(a1[c].mean())) for c in cols},
                a1_by_odor=dict(min=float(bo.bouton_share.min()), max=float(bo.bouton_share.max())),
                a2={r['set']: r['share_pn'] for r in a2.to_dict('records')}, seconds=time.time() - t0)
    _write_json(out / 'a_summary.json', summ)
    print(f'E94 A in {time.time() - t0:.0f} s', flush=True)
    return summ


# ---------------------------------------------------------------------------
# part B: the doses and forms at fixed scales

def stage_extras(W: pd.DataFrame, wiring: str, st_on: dict) -> dict:
    """A run's bouton and claw realism (`bouton_realism`, `claw_loss_shares`) against the fit's APL-off stage run."""
    off = fit_stage(wiring, 'off')
    check_stage_shapes(W, st_on, wiring)
    check_stage_shapes(W, off, f'{wiring} off')
    odors = _odors(st_on['bouton_win'])
    b_off = off['bouton_win'][odors].to_numpy(dtype=float)
    b_on = st_on['bouton_win'][odors].to_numpy(dtype=float)
    sh = claw_loss_shares(W, b_off, b_on, st_on['claw_win'][odors].to_numpy(dtype=float), kc_of_claw(W))
    row = dict(bouton_realism(b_off, b_on))
    row.update({(f'claw_{k}' if k.startswith('drive_removed') else k): v for k, v in sh.items() if k != 'by_odor'})
    row['bouton_share_by_odor'] = json.dumps({o: (None if not np.isfinite(s) else round(float(s), 6))
                                              for o, s in zip(odors, sh['by_odor'])})
    row['bouton_share_odor_min'] = float(np.nanmin(sh['by_odor']))
    row['bouton_share_odor_max'] = float(np.nanmax(sh['by_odor']))
    row['connectome_share'] = ANATOMICAL_SHARE
    apl = st_on['apl']['win'].astype(float)
    row['apl_win_mean'] = float(apl.mean())
    row['apl_rel_max'] = float((apl / apl.mean() - 1.0).abs().max())
    return row


def run_point(model: WiringModel94, v: str, sub: str = 'B', part: str = '', cc: Optional[dict] = None) -> Tuple[dict, dict]:
    """One setting (variant `v`) on a wiring at its fit (`b_reading94`), saved under outputs/.../<sub>/ with its stage
    matrices, drive and loss; a control ('..._ctrl') runs in the subtraction at the chosen setting's weighting with
    `control_matrix94`'s static inhibition and dose, from the chosen setting's run on the same wiring (which must exist).
    Returns (the run's row, the reading)."""
    w = model.label
    tag = run_tag(v, w)
    pv = parse_variant(v)
    cc = H.code_commits() if cc is None else cc
    wk, fac = E91._setting_weights(model, pv['weighting'])
    kc_index = model.saved['on'][1].index
    dose, form, ext, cform, cm = pv['dose'], pv['form'], None, '', None
    if pv['control']:
        d_ch = load_drive(run_tag(pv['control_of'], w), sub)
        d_ref, d_off = fit_drive(w, 'on', kc_index), fit_drive(w, 'off', kc_index)
        if len(d_ch) != len(kc_index):
            raise RuntimeError(f"{tag}: the chosen setting's drive has {len(d_ch)} KCs, the model {len(kc_index)}")
        d_ch.index = kc_index
        cm, cform, dose = control_matrix94(d_ref, d_ch[list(d_ref.columns)], d_off)
        ext = cm.iloc[:, 0].to_numpy(dtype=float)
        form = 'sub'
    res = b_reading94(model, wk, dose, form, ext=ext, label=tag)
    o = {'on': res['on'], 'off': res['off']}
    st = res['stages']
    odors = _odors(st['bouton_win'])
    drive = kc_drive(model.W, st['bouton_win'][odors].to_numpy(dtype=float), kc_index, odors)
    save_run(tag, o, sub, stages=st, drive=drive, loss=_frame(res['gain']['L'], res['on'][1]))
    resp_off = E89._responders(model.saved['off'][1])
    resp_on = E89._responders(res['on'][1])
    fac_on = {k.replace('fac_resp_', 'fac_respon_'): x for k, x in E89.factor_summary(fac, resp_on).items()
              if k.startswith('fac_resp_')}
    mult = aplpn_multiple(w, dose)
    row = dict(part=part, arm=ARM, tag=tag, variant=v, source='E94', screen=pv['screen'], form=pv['form'], run_form=form,
               exponent=pv['exponent'], dose=pv['dose'], run_dose=dose, control_of=pv['control_of'], control_form=cform,
               wiring=w, host=_host(), fit_folder=model.folder, fit_tag=model.tag, apl_factor=res['apl_factor'],
               at_bound=res['at_bound_apl'], bracket_lo=res['bracket'][0], bracket_hi=res['bracket'][1],
               n_on=E83._n(res['on'][0]), target_on=model.target_on, n_off=E83._n(res['off'][0]),
               target_off=model.target_off, n_sims=res['n_sims'], seconds=res['seconds'],
               **{k: res['params'][k] for k in E83.PARAM_KEYS}, aplpn_fit_multiple=aplpn_fit_multiple(w),
               aplpn_multiple=mult, aplpn_over_bound=mult / BOUND_E34, above_bound=bool(mult > BOUND_E34 + 1e-9),
               wAPLKC_sha256=E89.E80_sha(wk), K_fit=model.K,
               ext_inh_kc_max=float(ext.max()) if ext is not None else np.nan,
               ext_inh_kc_mean=float(ext.mean()) if ext is not None else np.nan,
               **E89.factor_summary(fac, resp_off), **fac_on, **res['gain']['row'], **cc)
    row.update(stage_extras(model.W, w, st))
    row.update(E83.realism(o))
    print(f"E94 {tag}: dose {dose:g} ({mult:.3g} x E3's base{', above the bound' if row['above_bound'] else ''}), "
          f"APL scale x {res['apl_factor']:.4g}{' (at a bound)' if res['at_bound_apl'] else ''}, "
          f"{row['n_on']} / {model.target_on} responses, {res['n_sims']} sims, {res['seconds']:.0f}s; boutons: "
          f"{row['bouton_share_removed']:.3f} removed, {row['bouton_zero_share']:.3f} at zero; share of the claws' loss at "
          f"the boutons {row['bouton_share']:.3f}; flag {row.get('realism_flag')!r}", flush=True)
    return row, res


def b_points(out: Path, wirings: Sequence[str], variants: Sequence[str], runs_name: str, part: str, sub: str = 'B'
             ) -> pd.DataFrame:
    """Settings (variant names) on the fits as they are, each by `run_point`; resumable (rows in `runs_name`). Needs part
    0 to have passed."""
    import gc
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identities must hold before any scored point)')
    cc = _commits()
    path = out / runs_name
    for w in wirings:
        todo = [v for v in variants if run_tag(v, w) not in _done(path)]
        if not todo:
            continue
        model = WiringModel94(w, check=False)
        for v in todo:
            row, _ = run_point(model, v, sub, part, cc)
            _append(path, row)
            _provenance_add(out, part, run_tag(v, w), dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **cc))
        del model
        gc.collect()
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def part_b(out: Path, arm: str, wirings: Sequence[str], points: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """B1, B2 or B3 (ARMS) on `wirings` (runs_<arm>.csv); `points` restricts to those variants of the arm."""
    vs = arm_variants(arm)
    if points:
        bad = set(points) - set(vs)
        if bad:
            raise ValueError(f'{sorted(bad)} are not points of {arm} ({vs})')
        vs = [v for v in vs if v in points]
    return b_points(out, wirings, vs, f'runs_{arm}.csv', arm)


def _provenance_add(out: Path, part: str, tag: str, prov: dict) -> None:
    path = out / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E94_bouton_share.py')
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def part_hostcheck(out: Path, v: str = 'B_sub_k2_d0.5', wiring: str = 'real') -> dict:
    """The host check (module docstring): this host's run of point `v` on `wiring` (run once) against every other copy."""
    host = _host()
    sub = f'hostcheck_{host}'
    tag = run_tag(v, wiring)
    path = out / f'runs_hostcheck_{host}.csv'
    if tag not in _done(path):
        b_points(out, [wiring], [v], path.name, 'hostcheck', sub)
    mine = load_run(tag, sub)
    r0 = pd.read_csv(path).set_index('tag').loc[tag]
    res = dict(host=host, tag=tag, apl_factor=float(r0['apl_factor']), others={})
    arc = H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT
    cands = _dirs('B') + sorted(p for p in list(H.outputs_dir(EXPERIMENT).glob('hostcheck_*')) + list(arc.glob('hostcheck_*'))
                                + list((arc / 'from_hal').glob('hostcheck_*')) if p.name != sub)
    for d in cands:
        o = E83._load_pair(d, tag)
        if o is None:
            continue
        res['others'][str(d)] = dict(on=_counts_diff(mine['on'], o['on']), off=_counts_diff(mine['off'], o['off']))
    for n in [p for p in out.glob('runs_*.csv') if p.name != path.name]:
        try:
            r = pd.read_csv(n)
        except pd.errors.EmptyDataError:
            continue
        r = r[r.tag == tag] if 'tag' in r else r.iloc[0:0]
        if len(r):
            res.setdefault('factors', {})[n.name] = float(r.apl_factor.iloc[-1])
    res['exact'] = bool(res['others'] and all(x['on']['exact'] and x['off']['exact'] for x in res['others'].values()))
    _write_json(out / f'hostcheck_{host}.json', res)
    print(f"E94 hostcheck {host} {tag}: {json.dumps(res, default=str)}", flush=True)
    return res


# ---------------------------------------------------------------------------
# the references and the runs table

def reference_rows() -> pd.DataFrame:
    """The references: the unchanged fits (B_sub_k2_d1, source 'saved', every wiring, with E90's stage runs' bouton and
    claw readings), E89's C1 point at exponent 1 (B_sub_k1_d1, source 'E89', S) and E92's B2 point (B_gain_k2_d1, source
    'E92', S); realism recomputed from the runs."""
    if 'refs' in _CACHE:
        return _CACHE['refs'].copy()
    rows = []
    for w in ALL_WIRINGS:
        f = E83.saved_fit(w)
        r = dict(part='saved', arm=ARM, tag=f"{f['folder']}:{f['tag']}", variant=ref_variant(), source='saved', screen='ref',
                 form='sub', run_form='sub', exponent=2.0, dose=1.0, run_dose=1.0, control_of='', wiring=w,
                 fit_folder=f['folder'], fit_tag=f['tag'], apl_factor=1.0, at_bound=False,
                 **{k: f['p'][k] for k in E83.PARAM_KEYS}, aplpn_fit_multiple=aplpn_fit_multiple(w),
                 aplpn_multiple=aplpn_fit_multiple(w), above_bound=False, **E83.realism(f['saved']))
        try:
            W, _ = wiring_W(w)
            r.update(stage_extras(W, w, fit_stage(w, 'on')))
        except Exception as e:                    # (the bouton readings beside; not needed for the scores)
            r['stage_note'] = repr(e)
        rows.append(r)
    p89 = H.results_dir(E89.EXPERIMENT) / 'runs_C1.csv'
    if p89.exists():
        d = pd.read_csv(p89)
        d = d[d.permute_seed.isna() & (d.weighting == 'F_k1')].drop_duplicates('tag', keep='last')
        for q in d.to_dict('records'):
            o = E89.load_run(q['tag'])
            if o is None or q['wiring'] not in S_WIRINGS:
                continue
            rows.append(dict(part='E89_C1', arm=ARM, tag=q['tag'], variant=ref_k1_variant(), source='E89', screen='ref',
                             form='sub', run_form='sub', exponent=1.0, dose=1.0, run_dose=1.0, control_of='',
                             wiring=q['wiring'], apl_factor=q.get('apl_factor'), at_bound=_as_bool(q.get('at_bound')),
                             **{k: q.get(k) for k in E83.PARAM_KEYS}, aplpn_fit_multiple=aplpn_fit_multiple(q['wiring']),
                             aplpn_multiple=aplpn_fit_multiple(q['wiring']), above_bound=False, **E83.realism(o)))
    p92 = H.results_dir(E92G.EXPERIMENT) / 'runs_B2.csv'
    if p92.exists():
        d = pd.read_csv(p92)
        d = d[d.variant == 'B_F_k2_bdiv'].drop_duplicates('tag', keep='last') if 'variant' in d else d.iloc[0:0]
        for q in d.to_dict('records'):
            o = E92G.load_run(q['tag'])
            if o is None:
                continue
            rows.append(dict(part='E92_B2', arm=ARM, tag=q['tag'], variant=e92_variant(), source='E92', screen='E92',
                             form='gain', run_form='gain', exponent=2.0, dose=1.0, run_dose=1.0, control_of='',
                             wiring=q['wiring'], apl_factor=q.get('apl_factor'), at_bound=_as_bool(q.get('at_bound')),
                             **{k: q.get(k) for k in E83.PARAM_KEYS}, aplpn_fit_multiple=aplpn_fit_multiple(q['wiring']),
                             aplpn_multiple=aplpn_fit_multiple(q['wiring']), above_bound=False,
                             **{c: q.get(c) for c in q if str(c).startswith('divisor_bouton_')}, **E83.realism(o)))
    df = pd.DataFrame(rows)
    _CACHE['refs'] = df
    return df.copy()


RUNS_FILES = ('runs_B1.csv', 'runs_B2.csv', 'runs_B3.csv', 'runs_C.csv')


def runs_table(out: Path) -> pd.DataFrame:
    """Every run of the readings: this experiment's (RUNS_FILES, the latest row per tag) and the references, with
    'flagged', 'at_bound' and the parsed setting."""
    parts = []
    for n in RUNS_FILES:
        if (out / n).exists():
            try:
                parts.append(pd.read_csv(out / n).drop_duplicates('tag', keep='last'))
            except pd.errors.EmptyDataError:
                pass
    parts.append(reference_rows())
    df = pd.concat([p for p in parts if len(p)], ignore_index=True, sort=False)
    df = df.drop_duplicates(['variant', 'wiring'], keep='first').reset_index(drop=True)
    pv = pd.DataFrame([parse_variant(v) for v in df.variant])
    for c in ('screen', 'form', 'exponent', 'dose', 'control', 'equal_kind'):
        df[c] = pv[c].to_numpy()
    df['flagged'] = [_flagged(x) for x in df.realism_flag]
    df['at_bound'] = [_as_bool(x) for x in df.at_bound]
    df['above_bound'] = [_as_bool(x) for x in df.above_bound]
    return df


# ---------------------------------------------------------------------------
# the readings

READOUTS = ('cal', 'calK', 'raw', 'zcell')


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
                from experiments import E93_missing_glomeruli_input as E93
                have[ro] = E93.zcell_distances(o)
            else:
                have[ro] = E83.pair_distances(o, ro)
    return have


def readings94(out: Path, runs: Optional[pd.DataFrame] = None, readouts: Sequence[str] = READOUTS) -> pd.DataFrame:
    """readings.csv (module docstring, score) and pairs.csv (each population's per-pair distances and effect under 'cal'
    and 'raw')."""
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
                rows.append(dict(arm=ARM, variant=variant, screen=pv['screen'], form=pv['form'], exponent=pv['exponent'],
                                 dose=pv['dose'], control=pv['control'], equal_kind=pv['equal_kind'], unit=unit,
                                 n_wirings=len(ws), readout=ro, n_flagged=int(sum(meta[w].flagged for w in ws)),
                                 n_at_bound=int(sum(meta[w].at_bound for w in ws)),
                                 n_above_bound=int(sum(meta[w].above_bound for w in ws)),
                                 **E92G.score_row92(d), **group_effects(d)))
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


def line94(out: Path, pr: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """line.csv: each population's ester and alcohol groups' APL-on distances relative to the 120 pairs' mean (E92 part
    A's `relative_distances`), placed against E92's part A line of S (`line(readout)['S']`, `place_on_line`); on T,
    relative to the unchanged fits' residual on T (E93's `line_reading`)."""
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
        if unit != 'S' and (ref_variant(), unit, ro) in pts:
            e0, a0 = pts[(ref_variant(), unit, ro)]
            ref = float(e0 - (L['intercept'] + L['slope'] * a0))
        rows.append(dict(variant=v, unit=unit, readout=ro, ester_rel=est, alcohol_rel=alc, slope=L['slope'],
                         intercept=L['intercept'], scatter=L['scatter'], max_abs_resid=L['max_abs_resid'], ref_resid=ref,
                         **A.place_on_line(est, alc, L, ref)))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'line.csv', index=False)
    return df


def _co_job94(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3 (E92's
    `_co_job92`)."""
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


def coresponse94(out: Path, n_workers: int = 4) -> pd.DataFrame:
    """coresponse.csv: per S setting (all four S runs) the changes of L (co-response) and R (recruitment) under APL of the
    four nearest ester pairs, the ester and the alcohol groups on S's population, against the flies' (E92's
    `coresponse92` on this experiment's runs; jobs cached under outputs/.../coresponse_cache.pkl)."""
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
            for res, j in zip(ex.map(_co_job94, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E94 coresponse {j[1]}', flush=True)
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
                        rows.append(dict(arm=ARM, variant=variant, reading=rdn, candidate=cand if rdn.startswith('cand') else '',
                                         subset=subset, factor=f, group=grp, ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])),
                                         flies_lo=float(np.exp(lo)), flies_hi=float(np.exp(hi)), differs=differs,
                                         known_error=known.get(rdn.split('_')[1], np.nan) if rdn.startswith('cand') else np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'coresponse.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# part C: the choice, the gate, the control, the equal rule

def _settings94(rd: pd.DataFrame, screens: Sequence[str] = CHOICE_SCREENS) -> pd.DataFrame:
    return rd[(rd.unit == 'S') & (rd.readout == 'cal') & rd.screen.isin(list(screens)) & ~rd.control.astype(bool)]


def _pat(rd: pd.DataFrame, variant: str, unit: str, target: str, col: str = 'pattern', readout: str = 'cal') -> float:
    return E89._pat(rd, variant, unit, target, col, readout)


PARTS_COLS = E92G.PARTS_COLS


def choose94(rd: pd.DataFrame, half: str, runs: Optional[pd.DataFrame] = None) -> dict:
    """The choice on S against half `half`: the setting with the lowest pattern part of S's population under the
    calibrated readout among the references (REF_VARIANTS) and B1 to B3's settings whose four runs are inside the bands and
    not at a bracket's bound; something to confirm if it is a B setting and gains CHOICE_MIN_GAIN or more over the
    unchanged fits. The parts of the pattern part, the APL->bouton multiple and the bouton-by-odor entries at zero of the
    chosen setting beside (`runs`, its S runs' mean; named above ZERO_SHARE_NOTE)."""
    g = _settings94(rd)
    ref = ref_variant()
    ok = g[g.variant.isin(REF_VARIANTS) | ((g.n_flagged == 0) & (g.n_at_bound == 0))]
    if ok.empty or ref not in set(g.variant):
        return dict(arm=ARM, half=half, chosen=None, confirm=False, note='runs missing')
    best = ok.loc[ok[f'half{half}_pattern'].idxmin()]
    rr = g[g.variant == ref].iloc[0]
    gain = float(rr[f'half{half}_pattern']) - float(best[f'half{half}_pattern'])
    confirm = bool(best.screen in ('B1', 'B2', 'B3') and gain >= CHOICE_MIN_GAIN)
    row = dict(arm=ARM, half=half, chosen=best.variant, screen=best.screen, form=best.form, exponent=float(best.exponent),
               dose=float(best.dose), pattern=float(best[f'half{half}_pattern']), ref_pattern=float(rr[f'half{half}_pattern']),
               gain=gain, confirm=confirm, n_eligible=int(len(ok)), n_settings=int(len(g)),
               n_above_bound=int(best.n_above_bound), excluded=json.dumps(sorted(set(g.variant) - set(ok.variant))))
    for c in PARTS_COLS[1:]:
        row[f'chosen_{c}'] = float(best[f'half{half}_{c}'])
        row[f'ref_{c}'] = float(rr[f'half{half}_{c}'])
    if runs is not None:
        q = runs[(runs.variant == best.variant) & runs.wiring.isin(S_WIRINGS)]
        if len(q) and 'bouton_zero_share' in q:
            z = float(q.bouton_zero_share.mean())
            row.update(zero_share=z, zero_share_named=bool(z > ZERO_SHARE_NOTE),
                       aplpn_multiple=float(q.aplpn_multiple.mean()) if 'aplpn_multiple' in q else np.nan)
    return row


def plan_c(rd: pd.DataFrame, runs: Optional[pd.DataFrame] = None) -> Dict[str, dict]:
    return {h: choose94(rd, h, runs) for h in HALVES}


def control_of(variant: str) -> str:
    pv = parse_variant(variant)
    return variant_name(pv['form'], pv['exponent'], pv['dose'], control=True)


def equal_rule94(rd: pd.DataFrame) -> List[dict]:
    """The equal rule: a B setting of one of the three kinds (`equal_kind`) is an equal of the unchanged fits if, on S
    under the readout, its pattern part is within EQUAL_PATTERN of theirs against all recordings and halves A and B, its
    four runs are inside the bands and none at a bracket's bound, its raw-count loss within EQUAL_RAW and its mean APL
    effect within EQUAL_MEAN of theirs; 'better' if its pattern part is lower on all three (and the rest holds). B1's
    settings above dose 1 are listed with kind None (not an equal: less defensible and nothing gained)."""
    ref = ref_variant()
    g = _settings94(rd)
    out = []
    if ref not in set(g.variant):
        return out
    me0 = float(g[g.variant == ref].mean_effect.iloc[0])
    for r in g.itertuples():
        if r.screen not in ('B1', 'B2', 'B3'):
            continue
        kind = r.equal_kind if isinstance(r.equal_kind, str) and r.equal_kind else None
        d = {t: _pat(rd, r.variant, 'S', t) - _pat(rd, ref, 'S', t) for t in ('halfA', 'halfB', 'z8')}
        raw = _pat(rd, r.variant, 'S', 'z8', 'loss', 'raw') - _pat(rd, ref, 'S', 'z8', 'loss', 'raw')
        me = float(r.mean_effect) - me0
        rest = bool(int(r.n_flagged) == 0 and int(r.n_at_bound) == 0 and abs(raw) <= EQUAL_RAW and abs(me) <= EQUAL_MEAN)
        eq = bool(kind is not None and all(abs(x) <= EQUAL_PATTERN for x in d.values()) and rest)
        out.append(dict(arm=ARM, variant=r.variant, screen=r.screen, form=r.form, exponent=float(r.exponent), dose=float(r.dose),
                        kind=kind, **{f'd_pattern_{t}': x for t, x in d.items()}, d_raw_loss=raw, d_mean_effect=me,
                        n_flagged=int(r.n_flagged), n_at_bound=int(r.n_at_bound), n_above_bound=int(r.n_above_bound),
                        equal=eq, better=bool(kind is not None and all(x < 0 for x in d.values()) and rest)))
    return out


def part_c(out: Path) -> dict:
    """C: the choice per half (choice.csv); every chosen setting with something to confirm run on T at each wiring's fit
    (runs_C.csv), then its control on S and T (resumable); then `score` (without the co-response)."""
    runs = runs_table(out)
    rd = readings94(out, runs, readouts=('cal', 'raw'))
    plan = plan_c(rd, runs)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    chosen = sorted({c['chosen'] for c in plan.values() if c['confirm']})
    if not chosen:
        print('E94 C: nothing to confirm', flush=True)
    for v in chosen:
        b_points(out, T_WIRINGS, [v], 'runs_C.csv', 'C')
    for v in chosen:
        for ws in (S_WIRINGS, T_WIRINGS):
            b_points(out, ws, [control_of(v)], 'runs_C.csv', 'C')
    return score(out, coresponse=False)


def verdict94(out: Path, rd: Optional[pd.DataFrame] = None, runs: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: the choices, the gate in four directions, the equals by kind, the dose response on S, the anatomical
    dose and D's plans, the flagged runs, the runs at a bracket's bound and above E34's bound, the gain's divisors."""
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table(out) if runs is None else runs
    plan = plan_c(rd, runs)
    pd.DataFrame(plan.values()).to_csv(out / 'choice.csv', index=False)
    ref = ref_variant()
    dirs = []
    for X in HALVES:
        c = plan[X]
        chosen = c['chosen'] if c['confirm'] else None
        dirs.append(E92G.gate_direction92(rd, runs, X, chosen, ref, [control_of(chosen)] if chosen else []))
    pd.DataFrame(dirs).to_csv(out / 'gate.csv', index=False)
    gate = E92G.gate92(dirs)
    eqs = equal_rule94(rd)
    pd.DataFrame(eqs).to_csv(out / 'equal.csv', index=False)
    gains_s = any(c['confirm'] for c in plan.values())
    v = dict(thresholds=dict(CHOICE_MIN_GAIN=CHOICE_MIN_GAIN, GATE_GAIN=E91.GATE_GAIN, GATE_SINGLES_SHARE=E91.GATE_SINGLES_SHARE,
                             GATE_RAW_SLACK=E91.GATE_RAW_SLACK, GATE_MEAN_SLACK=E91.GATE_MEAN_SLACK, GATE_SQ_ALL=E91.GATE_SQ_ALL,
                             GATE_SQ_SIX=E91.GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO, EQUAL_PATTERN=EQUAL_PATTERN,
                             EQUAL_RAW=EQUAL_RAW, EQUAL_MEAN=EQUAL_MEAN, APL_BRACKET=APL_BRACKET, BOUND_E34=BOUND_E34,
                             ZERO_SHARE_NOTE=ZERO_SHARE_NOTE, ANATOMICAL_SHARE=ANATOMICAL_SHARE),
             readout="cal (the calibrated readout, K set on each run's APL-on counts; E83's readings and gate)",
             reference=ref, choice=plan, gate=gate,
             verdict='passes' if gate['passes'] else ('gains on S, not through the gate' if gains_s else 'level'),
             equals=eqs, equals_by_kind={k: [e['variant'] for e in eqs if e['equal'] and e['kind'] == k] for k in EQUAL_KINDS})
    g = rd[(rd.unit == 'S') & (rd.readout == 'cal')]
    v['S'] = {r.variant: dict(screen=r.screen, form=r.form, exponent=r.exponent, dose=r.dose, halfA=r.halfA_pattern,
                              halfB=r.halfB_pattern, halfC=r.halfC_pattern, halfD=r.halfD_pattern, all=r.z8_pattern,
                              six=r.z6_pattern, loss=r.z8_loss, mean_effect=r.mean_effect, amplitude=r.z8_amplitude,
                              shape=r.z8_shape, sd_ratio=r.z8_sd_ratio, pearson=r.z8_pearson, rho=r.z8_rho,
                              err_ester=r.z8_err_ester, err_alcohol=r.z8_err_alcohol, n_flagged=int(r.n_flagged),
                              n_at_bound=int(r.n_at_bound), n_above_bound=int(r.n_above_bound)) for r in g.itertuples()}
    p0 = v['S'].get(ref, {}).get('all')
    b1 = [x for x in v['S'] if parse_variant(x)['screen'] == 'B1' and 0.1 <= parse_variant(x)['dose'] <= 4.0]
    if p0 is not None and b1:
        v['dose_within_equal'] = {x: bool(abs(v['S'][x]['all'] - p0) <= EQUAL_PATTERN) for x in b1}
    v['flagged_runs'] = runs[runs.flagged & (runs.source == 'E94')][['tag', 'realism_flag']].values.tolist()
    v['at_bound_runs'] = runs[runs.at_bound & (runs.source == 'E94')].tag.tolist()
    v['above_bound_runs'] = runs[runs.above_bound & (runs.source == 'E94')].tag.tolist()
    rr = runs[(runs.source == 'E94') & (runs.form == 'gain')]
    v['divisors'] = {vv: {c: float(q[c].max()) for c in ('divisor_bouton_p50', 'divisor_bouton_p90', 'divisor_bouton_max',
                                                           'divisor_bouton_share_above10') if c in q and q[c].notna().any()}
                     for vv, q in rr.groupby('variant')}
    try:
        v['anatomical_dose'] = anatomical_dose(runs)
    except Exception as e:                       # (reported, not fatal: B1's low doses may not be run yet)
        v['anatomical_dose'] = dict(error=repr(e))
    v['D2_plan'] = d2_plan(plan, eqs)
    if (out / 'd_scores.json').exists():
        v['D'] = _read_json(out / 'd_scores.json')
    if (out / 'coresponse.csv').exists():
        try:
            co = pd.read_csv(out / 'coresponse.csv')
            reach = co[(co.factor == 'L') & (co.group == 'alcohol')]
            v['alcohol_coresponse'] = {f'{r.variant}:{r.reading}:{r.subset}': dict(ratio=r.ratio, differs=bool(r.differs))
                                       for r in reach.itertuples()}
        except pd.errors.EmptyDataError:
            pass
    _write_json(out / 'verdict.json', v)
    print(f"E94 C: {v['verdict']}; " + '; '.join(f"{h} {plan[h].get('chosen')} (gain {plan[h].get('gain', np.nan):.4f})"
                                                for h in HALVES), flush=True)
    return v


def score(out: Path, n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """readings.csv, pairs.csv, line.csv, coresponse.csv (with `coresponse`), choice.csv, gate.csv, equal.csv,
    verdict.json, the figure and comparison_table.md."""
    runs = runs_table(out)
    rd = readings94(out, runs, readouts=READOUTS if zcell else tuple(r for r in READOUTS if r != 'zcell'))
    try:
        line94(out)
    except Exception as e:                       # (E92 part A's line needs its line.json; reported, not fatal)
        print(f'E94 score: line.csv left out ({e!r})', flush=True)
    if coresponse:
        coresponse94(out, n_workers)
    v = verdict94(out, rd, runs)
    figure(out)
    comparison_table(out)
    return v


# ---------------------------------------------------------------------------
# part D: the family fitted again

def anatomical_dose(runs: pd.DataFrame) -> dict:
    """D1's dose: of B1's ANATOMICAL_DOSES the one whose share of the claws' loss made at the boutons, mean over S's four
    runs, is nearest ANATOMICAL_SHARE (all four runs needed for each dose)."""
    shares = {}
    for d in ANATOMICAL_DOSES:
        q = runs[(runs.variant == variant_name('sub', 2.0, d)) & runs.wiring.isin(S_WIRINGS)]
        if set(q.wiring) != set(S_WIRINGS):
            raise RuntimeError(f'B1 at dose {d:g} is not run on all of S ({sorted(q.wiring)})')
        shares[d] = float(q.bouton_share.astype(float).mean())
    best = min(shares, key=lambda d: abs(shares[d] - ANATOMICAL_SHARE))
    return dict(dose=float(best), share=shares[best], shares={f'{k:g}': x for k, x in shares.items()}, target=ANATOMICAL_SHARE,
                holds={w: aplpn_multiple(w, best) for w in ALL_WIRINGS})


def d1_tag(label: str) -> str:
    """The tag E10_joint_fit.py gives a D1 fit (its prefix rules: the wiring, rate0.068_, aplboth_, aplpnhold_, then
    k2_d1)."""
    from experiments import E10_joint_fit as J
    a = E91.d_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return f"{wt}rate{H.data_target('ctrl_frac'):.3f}_aplboth_{J.aplpn_hold_tag(1.0)}k2{D1_TAG}"


def d_start_key(label: str, kind: str = 'd1') -> str:
    return f'{label}_{kind}'


def d1_start_params(out: Path, label: str, dose: float) -> dict:
    """D1's start on a wiring: the absolute parameters of its B1 point at `dose` (runs_B1.csv; the APL->bouton scale at
    the dose, the tied APL scales re-held), written into d_start_params.json under `d_start_key`."""
    rb = pd.read_csv(out / 'runs_B1.csv').drop_duplicates('tag', keep='last').set_index('tag')
    tag = run_tag(variant_name('sub', 2.0, dose), label)
    if tag not in rb.index:
        raise RuntimeError(f'{tag}: no B1 point (run --part B --arm B1 on this wiring first)')
    r = rb.loc[tag]
    if bool(r['at_bound']):
        raise RuntimeError(f'{tag}: the re-held point is at a bound of the APL bracket; no start for D1')
    start = {k: float(r[k]) for k in E83.PARAM_KEYS}
    path = out / D_START_FILE
    d = _read_json(path)
    d[d_start_key(label)] = dict(start, apl_factor=float(r['apl_factor']), dose=float(dose), from_tag=tag)
    _write_json(path, d)
    return d[d_start_key(label)]


def d1_command(label: str, dose: float, maxfev: int = D_MAXFEV) -> List[str]:
    """The E10_joint_fit.py command of one D1 fit: C2's family (`-k 2 --rate-target ctrl_frac --loss apl --apl-readout
    both`) with the APL->bouton scale held at the dose times the wiring's fitted multiple of E3's base (`--aplpn-hold`;
    the bouton->APL scale tied as in the family), into this experiment's folders, from the wiring's B1 point at the dose
    (`--init-params` on d_start_params.json)."""
    m = E91.d_member(label)
    F = aplpn_multiple(label, dose)
    return [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-hold',
            repr(float(F)), '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', D1_TAG,
            '--maxfev', str(int(maxfev)), '--experiment', EXPERIMENT,
            '--init-params', f'results/{EXPERIMENT}/{D_START_FILE}:{d_start_key(label)}'] + m['wiring_args']


def _fits(out: Path) -> pd.DataFrame:
    p = out / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


def d_allowed(out: Path, pop: str) -> Tuple[bool, str]:
    """D1 on S: whatever (C) says, once B1's low doses are run on S (the anatomical dose). On T: only if S's D1 population
    is an equal or better (`score_d`)."""
    if pop == 'S':
        try:
            ad = anatomical_dose(runs_table(out))
        except Exception as e:
            return False, f'no anatomical dose ({e})'
        return True, f"anatomical dose {ad['dose']:g} (share {ad['share']:.3f}, target {ANATOMICAL_SHARE})"
    sd = score_d(out, write=False)
    s = sd.get('S')
    if not s:
        return False, "S's D1 fits not scored"
    ok = bool(s.get('equal') or s.get('better'))
    return ok, f"S's D1 population: equal {s.get('equal')}, better {s.get('better')}"


def part_d(out: Path, pop: str = 'S', wirings: Optional[Sequence[str]] = None, force: bool = False,
           dry_run: bool = False) -> List[str]:
    """D1's fits of population `pop` (or `wirings` of it: a lane), one process after another (`d1_command`; logs under
    outputs/.../D/), skipping those in fits.csv; refuses what the entry's rules do not allow unless `force`."""
    ok, why = d_allowed(out, pop)
    print(f'E94 D1 {pop}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
    if not ok and not force:
        return []
    dose = anatomical_dose(runs_table(out))['dose']
    ws = [w for w in _wirings(pop) if wirings is None or w in wirings]
    if not dry_run:
        _commits()
    logs = local_dir('D')
    tags = []
    for w in ws:
        tag = d1_tag(w)
        tags.append(tag)
        if tag in set(_fits(out).tag):
            print(f'E94 D1 {tag}: done', flush=True)
            continue
        if not dry_run:
            d1_start_params(out, w, dose)
        cmd = d1_command(w, dose)
        print('E94 D1:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E94 D1 {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits(out).tag):
            raise RuntimeError(f'E94 D1: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def _d_row(sc: dict, sc0: dict, raw: dict, raw0: dict) -> dict:
    return dict(pattern_all=sc['z8_pattern'], pattern_all_unchanged=sc0['z8_pattern'], halfA=sc['halfA_pattern'],
                halfB=sc['halfB_pattern'], halfA_unchanged=sc0['halfA_pattern'], halfB_unchanged=sc0['halfB_pattern'],
                loss=sc['z8_loss'], mean_effect=sc['mean_effect'], amplitude=sc['z8_amplitude'], shape=sc['z8_shape'],
                sd_ratio=sc['z8_sd_ratio'], pearson=sc['z8_pearson'], rho=sc['z8_rho'], raw_loss=raw['z8_loss'],
                raw_loss_unchanged=raw0['z8_loss'], err_ester=sc['z8_err_ester'], err_alcohol=sc['z8_err_alcohol'])


def d1_rule(r: dict) -> Tuple[bool, bool]:
    """(equal, better) of S's D1 population (`score_d`'s row): within D_EQUAL of UNCHANGED_S against all recordings and
    halves A and B ('better': below all three), every fit inside the bands, the raw-count loss within D_RAW_SLACK of the
    unchanged fits', the mean APL effect within D_MEAN_SLACK of the flies'."""
    d = {'all': r['pattern_all'] - UNCHANGED_S['all'], 'A': r['halfA'] - UNCHANGED_S['A'], 'B': r['halfB'] - UNCHANGED_S['B']}
    rest = bool(r['n_flagged'] == 0 and abs(r['raw_loss'] - r['raw_loss_unchanged']) <= D_RAW_SLACK
                and abs(r['mean_effect'] - FLIES_MEAN) <= D_MEAN_SLACK)
    return bool(all(abs(x) <= D_EQUAL for x in d.values()) and rest), bool(all(x < 0 for x in d.values()) and rest)


def score_d(out: Path, write: bool = True) -> dict:
    """D1's scores: per wiring and population (S, T, the eight) the readings under the readout and raw counts against the
    unchanged fits of the same wirings, with realism; S by `d1_rule`; the eight against UNCHANGED_8 by the same rule
    (within D_EQUAL) and by E91's part D pass rule. d_scores.csv, d_scores.json."""
    from experiments import E72_wiring_population as E72
    fits = _fits(out)
    rows, st = [], {}
    for w in ALL_WIRINGS:
        try:
            tag = d1_tag(w)
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
        row = dict(unit=w, tag=tag, **_d_row(sc, sc0, raw, raw0), lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                   **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                             'kc_kc_inh', 'aplpn_hold', 'aplpn_factor_vs_e3', 'init_params', 'n_evals')},
                   fit_loss=fr.get('loss'), **E83.realism(o))
        rows.append(row)
        st[w] = (d, d0, row)
    res = {}
    for unit, ws in (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS)):
        if not all(w in st for w in ws):
            continue
        dp = {ro: E83.population([st[w][0][ro] for w in ws]) for ro in ('cal', 'raw')}
        d0 = {ro: E83.population([st[w][1][ro] for w in ws]) for ro in ('cal', 'raw')}
        sc, sc0 = E92G.score_row92(dp['cal']), E92G.score_row92(d0['cal'])
        raw, raw0 = E83.effect_scores(dp['raw'], E89.targets89()), E83.effect_scores(d0['raw'], E89.targets89())
        r = dict(_d_row(sc, sc0, raw, raw0), n_lower=int(sum(st[w][2]['lower'] for w in ws)),
                 n_flagged=int(sum(_flagged(st[w][2].get('realism_flag')) for w in ws)))
        if unit == 'S':
            r['equal'], r['better'] = d1_rule(r)
        if unit == 'all8':
            rest = bool(r['n_flagged'] == 0 and abs(r['raw_loss'] - r['raw_loss_unchanged']) <= D_RAW_SLACK
                        and abs(r['mean_effect'] - FLIES_MEAN) <= D_MEAN_SLACK)
            r['equal'] = bool(abs(r['pattern_all'] - UNCHANGED_8) <= D_EQUAL and rest)
            r['passes'] = bool(r['pattern_all'] <= UNCHANGED_8 - D_PASS_GAIN and r['n_lower'] >= D_WIRINGS_LOWER
                               and r['n_flagged'] == 0 and r['raw_loss'] <= D_RAW_REF + D_RAW_SLACK
                               and abs(r['mean_effect'] - FLIES_MEAN) <= D_MEAN_SLACK)
        res[unit] = r
        rows.append(dict(unit=unit, **r))
    if write:
        pd.DataFrame(rows).to_csv(out / 'd_scores.csv', index=False)
        _write_json(out / 'd_scores.json', res)
    return res


def d2_plan(plan: Dict[str, dict], eqs: Sequence[dict]) -> dict:
    """Whether D2 is due (the entry): on a choice by (C) at a dose above 1 (B1 or B2) or of the gain (B3), or on an equal
    of B2's or B3's; per arm the reasons."""
    why: Dict[str, List[str]] = {}
    for h, c in plan.items():
        if not c.get('confirm'):
            continue
        sc = c.get('screen')
        if (sc in ('B1', 'B2') and float(c.get('dose', 0)) > 1.0) or sc == 'B3':
            why.setdefault(sc, []).append(f"choice on half {h}: {c['chosen']}")
    for e in eqs:
        if e['equal'] and e['screen'] in ('B2', 'B3'):
            why.setdefault(e['screen'], []).append(f"equal ({e['kind']}): {e['variant']}")
    return dict(due=bool(why), arms=sorted(why), why=why)


def d2_command(label: str, arm: str, maxfev: int = D_MAXFEV) -> List[str]:
    """D2's fit for B1 (`-k 2 --aplpn-max 12`) or B2 (`-k 1 --aplpn-max 12`), from the wiring's own fit (E91's `d_member`
    `--init-from`, exactly at the same exponent). B3's needs `--bouton-divisive` (a new E10 option, not written)."""
    if arm == 'B3':
        raise NotImplementedError("D2 for the gain needs E10_joint_fit.py --bouton-divisive, written only if D2 is due for it")
    k = {'B1': 2.0, 'B2': 1.0}[arm]
    m = E91.d_member(label)
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', f'{k:g}', '--aplpn-max',
           f'{D2_APLPN_MAX:g}', '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', D2_TAG,
           '--maxfev', str(int(maxfev)), '--experiment', EXPERIMENT, '--init-from', m['tag'],
           '--init-from-experiment', m['folder']] + m['wiring_args']
    if k == m['k']:
        cmd.append('--init-from-exact')
    return cmd


def part_d2(out: Path) -> None:
    """D2 (a stub): prints whether it is due (verdict.json's D2_plan) and S's commands; S first, T's four only if S's
    fitted population is D_S_GAIN or more below UNCHANGED_S['all']; it passes as a population of eight by E91's part D rule.
    Written when it is due."""
    v = _read_json(out / 'verdict.json')
    p = v.get('D2_plan', {})
    print(f'E94 D2: {json.dumps(p)}', flush=True)
    for arm in p.get('arms', []):
        for w in S_WIRINGS:
            try:
                print('E94 D2:', ' '.join(d2_command(w, arm)), flush=True)
            except NotImplementedError as e:
                print(f'E94 D2 {arm}: {e}', flush=True)
    raise NotImplementedError('E94 part D2 is written when it is due (a choice above dose 1 or of the gain, or an equal of '
                              'B2 or B3)')


def part_e(out: Path) -> None:
    """Part E (only on a pass by (C) or D2): the firing-ceiling family with the response rate with APL on re-held, S at the
    flies' 2 s pulse with both response rates re-held, S through the trial noise as E72's part E ran it, and the summed
    form at exponent 1 (E91's next step 2) on the new family; each holds if the population keeps at least half of its
    gain. Not written: it runs only on a pass."""
    raise NotImplementedError('E94 part E runs only on a pass by part C or D2; it is not written')


# ---------------------------------------------------------------------------
# the steps from the claws to the somata, for a chosen or equal setting

def steps_variants(out: Path) -> List[str]:
    """The settings the steps are read for: every chosen setting with something to confirm and every equal (verdict.json)."""
    v = _read_json(out / 'verdict.json')
    vs = [c['chosen'] for c in v.get('choice', {}).values() if c.get('confirm')]
    vs += [e['variant'] for e in v.get('equals', []) if e.get('equal')]
    return sorted(dict.fromkeys(vs))


def part_steps(out: Path, variants: Optional[Sequence[str]] = None) -> Optional[pd.DataFrame]:
    """The steps (module docstring) for `variants` (default `steps_variants`) beside C2's family on S: somata_steps.csv,
    claws_steps.csv, steps.csv."""
    from experiments import E90_model_stages as E90
    from experiments import E92_steps_like_for_like as A
    variants = list(variants) if variants else steps_variants(out)
    if not variants:
        print('E94 steps: no chosen or equal setting', flush=True)
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
            st_on = load_stage94(r.tag) if r.source == 'E94' else (dict(fit_stage(w, 'on')) if r.source == 'saved' else None)
            if st_on is None:
                raise RuntimeError(f'{r.tag}: no stage matrices')
            st_off = dict(fit_stage(w, 'off'))
            W, _ = wiring_W(w)
            meta = {k: np.asarray(x) for k, x in E92G.units_meta(W).items()}
            rr = E90.read_wiring(w, st_on, st_off, meta, np.random.default_rng(E90.NOISE_SEED + STEPS_NOISE_SEED + i))
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
    print(f'E94 steps in {time.time() - t0:.0f} s', flush=True)
    return rows


# ---------------------------------------------------------------------------
# the archive, the figure and the table

def write_variants(out: Path) -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): one row per run of this experiment with its absolute parameters, the dose,
    the form, the weighting, the control's static inhibition, the bouton readings and realism; D's fits' distances, which
    E10_joint_fit.py wrote, are kept."""
    runs = runs_table(out)
    runs = runs[runs.source == 'E94']
    rows, parts = [], []
    for _, r in runs.iterrows():
        form = {'sub': 'a subtraction at the boutons (max(r - i, 0))',
                'gain': 'a gain at the boutons (r / (1 + i); olfsysm bouton_inh_divisive)'}[r.run_form]
        ctl = str(r.control_of) not in ('', 'nan')
        kind = (f"the fit as it is; APL->KC weights x ((F + 1) / mean(F + 1))^{r.exponent:g}; APL->bouton scale x {r.run_dose:g} "
                f"of the fit's ({r.aplpn_multiple:.3g} x E3's base); APL's inhibition of the boutons {form}"
                + (f"; the control of {r.control_of} ({r.control_form}): a static inhibition at the soma (ext_inh_kc) per KC"
                   if ctl else '') + "; the tied APL scale re-held; APL-off run the fit's")
        rows.append(dict(tag=r.tag, variant=r.tag, part=r.part, arm=ARM, kind=kind, config=r.variant, form=r.form,
                         run_form=r.run_form, breadth_def='F', breadth_k=r.exponent, k=r.exponent, dose=r.dose,
                         run_dose=r.run_dose, bouton_inh_divisive=bool(FORMS[r.run_form]), wiring=r.wiring,
                         ref_experiment=r.fit_folder, ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, apl_factor=r.apl_factor, aplpn_multiple=r.aplpn_multiple,
                         above_bound=bool(r.above_bound), rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, wAPLKC_sha256=r.wAPLKC_sha256, ext_inh_kc_max=r.get('ext_inh_kc_max'),
                         control_of=r.control_of, control_form=r.get('control_form'),
                         bouton_share_removed=r.get('bouton_share_removed'), bouton_zero_share=r.get('bouton_zero_share'),
                         bouton_share=r.get('bouton_share'), divisor_bouton_max=r.get('divisor_bouton_max'),
                         realism_flag=r.realism_flag, host=r.get('host'), subdir='B'))
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
        keep = [old[old.tag.isin(set(_fits(out).tag))]]
    if keep or parts:
        pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def figure(out: Path) -> Optional[Path]:
    """E94_dose.png: left, S's pattern part under the readout against the dose (B1, B2, B3; against half A, half B and
    all recordings; the unchanged fits dashed; hollow: a setting with a run outside the bands or at a bracket's bound);
    right, the share of the claws' loss made at the boutons against the dose (S's runs' mean), the connectome's share
    shaded."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not (out / 'readings.csv').exists():
        return None
    rd = pd.read_csv(out / 'readings.csv')
    g = _settings94(rd, ('ref', 'B1', 'B2', 'B3'))
    if g.empty:
        return None
    runs = runs_table(out)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    cols = {'halfA': '#d1495b', 'halfB': '#2e86ab', 'z8': '0.2'}
    styles = {'B1': ('-', 'o', 'subtraction, exponent 2'), 'B2': ('--', 's', 'subtraction, exponent 1'),
              'B3': (':', '^', 'gain, exponent 2')}
    ref = g[g.variant == ref_variant()]
    for scr, (ls, mk, lab) in styles.items():
        q = g[g.screen == scr].sort_values('dose')
        if scr in ('B1', 'B2'):
            r1 = g[g.variant == (ref_variant() if scr == 'B1' else ref_k1_variant())]
            q = pd.concat([q, r1]).sort_values('dose')
        if q.empty:
            continue
        bad = ((q.n_flagged > 0) | (q.n_at_bound > 0)).to_numpy()
        for t, c in cols.items():
            y = q[f'{t}_pattern'].to_numpy()
            ax.plot(q.dose, y, ls, color=c, lw=0.9, label=f'{lab}, {"all recordings" if t == "z8" else "half " + t[-1]}')
            ax.scatter(q.dose, y, s=26, marker=mk, facecolors=['none' if b else c for b in bad], edgecolors=c)
    for t, c in cols.items():
        if len(ref):
            ax.axhline(float(ref[f'{t}_pattern'].iloc[0]), color=c, ls='--', lw=0.6)
    ax.set_xscale('symlog', linthresh=0.05)
    ax.set_xlabel("dose (x the fitted APL->bouton scale)")
    ax.set_ylabel('pattern part of the APL-effect loss')
    ax.set_title("S's population, calibrated readout\n(dashed: the unchanged fits; hollow: outside the bands or at a bound)",
                 fontsize=9)
    ax.legend(fontsize=6, ncol=1)
    ax = axes[1]
    e = runs[runs.wiring.isin(S_WIRINGS) & ~runs.control.astype(bool)]
    if 'bouton_share' in e:
        for scr, (ls, mk, lab) in styles.items():
            refv = {'B1': ref_variant(), 'B2': ref_k1_variant()}.get(scr, '')
            q = e[(e.screen == scr) | (e.variant == refv)]
            q = q.groupby('dose').bouton_share.mean().dropna()
            if len(q):
                ax.plot(q.index, q.to_numpy(), ls, marker=mk, color='0.25', label=lab)
    ax.axhspan(CONNECTOME_SHARE[0], CONNECTOME_SHARE[1], color='#8fbc8f', alpha=0.4, label="connectome's share")
    ax.set_xscale('symlog', linthresh=0.05)
    ax.set_xlabel("dose (x the fitted APL->bouton scale)")
    ax.set_ylabel("share of the claws' loss made at the boutons")
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = out / 'E94_dose.png'
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _f(x, spec: str = '.4f') -> str:
    return E92G._f(x, spec)


def _md(df: pd.DataFrame, fmt: str = '.4f') -> str:
    try:
        return df.to_markdown(index=False, floatfmt=fmt)
    except ImportError:
        return '```\n' + df.to_string() + '\n```'


def comparison_table(out: Path) -> str:
    """comparison_table.md from the parts' tables."""
    L = ["# E94: the PN boutons' share of APL's inhibition", '',
         f'Generated by experiments/E94_bouton_share.py at {H.code_commits().get("queue_commit", "")[:7]}.', '']
    if (out / 'identity.json').exists():
        s = _read_json(out / 'identity.json')
        L += ['## Part 0: the identities', '', f"passes: {s.get('passes')}", '']
        for k, x in s.items():
            if not isinstance(x, dict) or k == 'commits':
                continue
            keep = {kk: x.get(kk) for kk in ('exact', 'passes', 'counts_differ', 'responses_differ', 'wiring_equal',
                                              'apl_factor', 'apl_factor_e92', 'factor_rel_diff', 'changed_from_fit')
                    if kk in x}
            if 'on' in x:
                keep.update(on_counts_differ=x['on']['counts_differ'], off_counts_differ=x['off']['counts_differ'])
            L.append(f'- {k}: {json.dumps(keep, default=str)}')
        L.append('')
    if (out / 'a1_path.csv').exists():
        a1 = pd.read_csv(out / 'a1_path.csv')
        keep = ['wiring', 'aplpn_fit_multiple', 'bouton_sub_p10', 'bouton_sub_p50', 'bouton_sub_p90', 'bouton_share_removed',
                'bouton_zero_share', 'spont_off', 'spont_on', 'sub_cv_over_odors_median_bouton', 'apl_rel_max',
                'drive_removed', 'bouton_share', 'bouton_share_1-4', 'bouton_share_5-8', 'bouton_share_9+']
        L += ['## A1: the bouton path in the fits (E90\'s stage runs, window means over the pulse)', '',
              _md(a1[[c for c in keep if c in a1]], '.3f'), '']
    if (out / 'a2_connectome.csv').exists():
        L += ["## A2: the connectome's share of APL's calyx synapses onto PN boutons", '',
              _md(pd.read_csv(out / 'a2_connectome.csv').drop(columns=['tables'], errors='ignore'), '.4f'), '']
    if (out / 'a3_summary.csv').exists():
        a3 = pd.read_csv(out / 'a3_summary.csv')
        q = a3[a3.group.isin(['ester', 'alcohol', 'ester_alcohol', 'near', 'other'])]
        pv = q.pivot_table(index=['form', 'dose'], columns=['measure', 'signal', 'group'], values='change')
        pv.columns = ['/'.join(map(str, c)) for c in pv.columns]
        sh = a3.groupby(['form', 'dose'])[['share_removed', 'zero_share']].first()
        L += ["## A3: the bouton stage under other doses (mean of the eight wirings; change against APL off)", '',
              _md(sh.join(pv).reset_index(), '.3f'), '']
    if (out / 'a4_archive.csv').exists():
        a4 = pd.read_csv(out / 'a4_archive.csv')
        keep = ['folder', 'tag', 'aplpn_multiple', 'readout', 'z8_pattern', 'z8_amplitude', 'z8_shape', 'z8_sd_ratio',
                'z8_pearson', 'z8_rho', 'mean_effect', 'z8_err_ester', 'z8_err_alcohol', 'eff_ester', 'eff_alcohol',
                'rate_on', 'realism_flag']
        L += ['## A4: the archive (all recordings)', '', _md(a4[[c for c in keep if c in a4]], '.3f'), '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        runs = runs_table(out)
        g = _settings94(rd, ('ref', 'E92', 'B1', 'B2', 'B3', 'control')).sort_values(['screen', 'dose'])
        graw = rd[(rd.unit == 'S') & (rd.readout == 'raw')].set_index('variant')
        rr = runs.groupby('variant')
        L += ["## B: S's population (calibrated readout; the parts of the pattern part against all recordings)", '',
              '| setting | pattern, half A | half B | all | amplitude | shape | spread ratio | Pearson r | ρ | loss | '
              'mean APL effect | raw-count loss | ester / alcohol error | APL->bouton x E3 | tied factor | share at the '
              'boutons | boutons at zero | flagged / at bound / above bound |',
              '|---|' + '---:|' * 17]
        for r in g.itertuples():
            q = rr.get_group(r.variant) if r.variant in rr.groups else pd.DataFrame()

            def rng(c, spec='.3g'):
                if not len(q) or c not in q or not q[c].notna().any():
                    return '–'
                return f"{q[c].min():{spec}}-{q[c].max():{spec}}"
            L.append(f'| {r.variant} | {r.halfA_pattern:.4f} | {r.halfB_pattern:.4f} | {r.z8_pattern:.4f} | '
                     f'{r.z8_amplitude:.4f} | {r.z8_shape:.4f} | {r.z8_sd_ratio:.2f} | {r.z8_pearson:.2f} | {r.z8_rho:.2f} | '
                     f'{r.z8_loss:.4f} | {r.mean_effect:+.3f} | {graw.z8_loss.get(r.variant, np.nan):.4f} | '
                     f'{r.z8_err_ester:+.3f} / {r.z8_err_alcohol:+.3f} | {rng("aplpn_multiple")} | {rng("apl_factor")} | '
                     f'{rng("bouton_share", ".3f")} | {rng("bouton_zero_share", ".3f")} | '
                     f'{r.n_flagged} / {r.n_at_bound} / {r.n_above_bound} |')
        L.append('')
        mine = runs[runs.source == 'E94']
        if len(mine):
            L += ['## Realism of the runs (rule of 2026-11-10 16:10; the earlier rule beside)', '',
                  '| run | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | largest counts on / off | '
                  'APL->bouton x E3 | boutons removed / at zero | share at the boutons (1-4 / 5-8 / 9+ claws) | '
                  'bouton divisor p10 / p50 / p90 / max | flag | flag (earlier rule) |',
                  '|---|' + '---:|' * 10 + '---|---|']
            for r in mine.itertuples():
                dv = ' / '.join(_f(getattr(r, f'divisor_bouton_{n}', np.nan), '.3g') for n in ('p10', 'p50', 'p90', 'max'))
                L.append(f'| {r.tag} | {r.rate_on:.4f} | {r.rate_off:.4f} | {r.ratio:.2f} | {r.rate_rho_on:.2f} | '
                         f'{r.opk_on_tv:.3f} / {r.opk_off_tv:.3f} | {int(r.max_count_on)} / {int(r.max_count_off)} | '
                         f'{r.aplpn_multiple:.3g}{" (above the bound (E34))" if r.above_bound else ""} | '
                         f'{_f(r.bouton_share_removed, ".3f")} / {_f(r.bouton_zero_share, ".3f")} | {_f(r.bouton_share, ".3f")} '
                         f"({_f(getattr(r, 'bouton_share_1-4', np.nan), '.3f')} / {_f(getattr(r, 'bouton_share_5-8', np.nan), '.3f')}"
                         f" / {_f(getattr(r, 'bouton_share_9+', np.nan), '.3f')}) | {dv} | "
                         f'{r.realism_flag if isinstance(r.realism_flag, str) else ""} | '
                         f'{r.realism_flag_zcrit if isinstance(r.realism_flag_zcrit, str) else ""} |')
            L.append('')
    for name in ('line.csv', 'choice.csv', 'gate.csv', 'equal.csv', 'd_scores.csv', 'steps.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                L += [f'## {name}', '', '(none)', '']
                continue
            if name == 'line.csv':
                d = d[[c for c in d.columns if c in ('variant', 'unit', 'readout', 'ester_rel', 'alcohol_rel', 'resid',
                                                       'off_toward_flies', 'off_toward_flies_nr')]]
            if name == 'steps.csv' and len(d) > 200:
                d = d[d.read.astype(bool)] if 'read' in d else d
            L += [f'## {name}', '', _md(d), '']
    if (out / 'verdict.json').exists():
        vc = _read_json(out / 'verdict.json')
        short = {k: vc.get(k) for k in ('verdict', 'gate', 'equals_by_kind', 'dose_within_equal', 'anatomical_dose', 'D2_plan',
                                         'D', 'flagged_runs', 'at_bound_runs', 'above_bound_runs', 'divisors')}
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(out: Path, part: str, **extra) -> None:
    """provenance_<part>.json's dated records: the commits, olfsysm's build, the model's provenance with the lab's input
    for the simulation parts, the settings and `extra` (part B's per-point records are in its 'points')."""
    import olfsysm
    path = out / f'provenance_{part}.json'
    old = _read_json(path)
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), host=_host(), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('0', 'B1', 'B2', 'B3', 'hostcheck', 'C', 'D'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    old.update(script='experiments/E94_bouton_share.py', part=part, runs=runs,
               settings=dict(FORMS=FORMS, ARMS={k: [v[0], v[1], list(v[2])] for k, v in ARMS.items()}, APL_BRACKET=APL_BRACKET,
                             APL_STEP=APL_STEP, BOUND_E34=BOUND_E34, CONNECTOME_SHARE=CONNECTOME_SHARE,
                             ANATOMICAL_SHARE=ANATOMICAL_SHARE, ANATOMICAL_DOSES=ANATOMICAL_DOSES, ZERO_SHARE_NOTE=ZERO_SHARE_NOTE,
                             DIVISOR_HIGH=DIVISOR_HIGH, CHOICE_MIN_GAIN=CHOICE_MIN_GAIN,
                             EQUAL=dict(PATTERN=EQUAL_PATTERN, RAW=EQUAL_RAW, MEAN=EQUAL_MEAN), D_MAXFEV=D_MAXFEV,
                             D_EQUAL=D_EQUAL, D2_APLPN_MAX=D2_APLPN_MAX))
    _write_json(path, old)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'A', 'A1', 'A2', 'A3', 'A4', 'B', 'hostcheck', 'score', 'C', 'D',
                                                      'Dscore', 'D2', 'E', 'steps', 'variants', 'table'))
    ap.add_argument('--arm', choices=tuple(ARMS), default=None, help='B: B1, B2 or B3')
    ap.add_argument('--wirings', default=None, help="0 / B / D: 'S', 'T', 'all8' or a comma-separated list (0: all8; B: S)")
    ap.add_argument('--points', default=None, help='B: comma-separated variants of the arm (default: all of them)')
    ap.add_argument('--point', default='B_sub_k2_d0.5', help='hostcheck: the variant')
    ap.add_argument('--wiring', default='real', help='hostcheck: the wiring')
    ap.add_argument('--pop', default='S', choices=('S', 'T'), help='D: which fits')
    ap.add_argument('--force', action='store_true', help="D: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='D: print the commands only')
    ap.add_argument('--workers', type=int, default=4, help='score: co-response processes')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    ap.add_argument('--no-zcell', action='store_true', help="score: skip E73's per-KC reading")
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
        {'A1': a1_path, 'A2': a2_connectome, 'A3': a3_bouton_stage, 'A4': a4_archive}[a.part](out)
    elif a.part == 'B':
        if a.arm is None:
            ap.error('--part B needs --arm')
        ws = _wirings(a.wirings or 'S')
        pts = [x for x in a.points.split(',') if x] if a.points else None
        provenance(out, a.arm, wirings=list(ws), points=pts)
        part_b(out, a.arm, ws, pts)
    elif a.part == 'hostcheck':
        provenance(out, 'hostcheck', point=a.point, wiring=a.wiring)
        part_hostcheck(out, a.point, a.wiring)
    elif a.part == 'score':
        score(out, a.workers, coresponse=not a.no_co, zcell=not a.no_zcell)
    elif a.part == 'C':
        provenance(out, 'C')
        part_c(out)
    elif a.part == 'D':
        if not a.dry_run:
            provenance(out, 'D', pop=a.pop, wirings=a.wirings)
        part_d(out, a.pop, None if a.wirings is None else _wirings(a.wirings), force=a.force, dry_run=a.dry_run)
    elif a.part == 'Dscore':
        score_d(out)
        verdict94(out)
        comparison_table(out)
    elif a.part == 'D2':
        part_d2(out)
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
    print(f'E94 {a.part} done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
