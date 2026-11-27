#!/usr/bin/env python3
"""E144, the weaker antennal-lobe normalization and the rise of the inhibition among KCs on the male CNS's two calyces
(experiment_priorities.md, E144; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 10:31).

A1's three changes (the antennal lobe's normalization weakened to pn_inhadd 1000 with pn_inhsc by E77's rule,
APL->bouton at x0.155, the KC-KC scale at x3.24) read on a second fly's calyces, where nothing of A1 was chosen: the male
CNS's left and right calyx, each with its own KC-KC matrix (E130's), and the hemibrain with 'argmax' claws as the
control. Nothing is fitted: E130's fits of the family on each calyx are used as they are.

Definitions (the entry's, fixed in the lock before any read, `lock_dict`; results/E144_a1_malecns_calyces/
E144_lock.json and its committed copy experiments/E144_lock.json):
- the calyces (CALYCES): H the hemibrain with 'argmax' claws on its real wiring (E130's refit); L and R the male CNS's
  left and right calyx, each on its real wiring and three within-type claw permutations (E133's seeds 1 to 3;
  WIRINGS_OF). A permutation runs at its calyx's fit's scales;
- the settings (SETTINGS), each from the calyx's own family (E130's counted fit of the calyx, `calyx_fit`): 'family'
  (olfsysm's pn_inhadd, 31.4), F (pn_inhadd 1000, pn_inhsc by E77's rule, `pn_terms`), B (F with APL->bouton at
  x APLPN), D1.8, A1 and D5.8 (B with the calyx's own matrix at x RISE**0.5, x RISE and x RISE**1.5 of its fitted
  scale, in both conditions) and Dabs on L only (B with L's matrix at the absolute scale A1 reaches on the hemibrain's
  real wiring, `dabs_scale`);
- the operating points (E136's, `read_targets`): megamat at the lock's 0.0915, validation2 at the lock's 0.050, the
  mixtures at E134's held mixture rate (0.054), each with the calyx's fit's APL-off / APL-on ratio for APL off, and the
  APL16 input at the calyx's fit's two response counts; both rates re-held at every read (E77's `hold_rates`: the
  threshold with APL off, then the APL->KC scale with KC->APL and bouton->APL tied);
- a point is a setting read on one calyx's wiring and one panel's input (`Point`, `point_tag`); a read is that point
  held, built again at the held parameters and run with every stage saved (E134's `run_stages`), which must give the
  hold's outputs exactly, saved in E139's layout (kcq.stage_readings' `save_step3_run`) so that E141's score reads it.

Judgments the entry leaves open (the narrowest reading; each is in the lock's 'judgments'):
- the start of every hold is the calyx's fit's absolute parameters (its fits.csv row) with the setting's factors, the
  threshold and the tied APL scales then held (E136's reads start from the family's held megamat point, which the male
  calyces do not have; the hold is exact, so the start moves only its path);
- the breadth behind the APL->KC weights is the calyx wiring's APL-off breadth on the lab's APL16 input (E133's
  `calyx_tables`: recomputed on a permutation, as E133's step 1 and the family's P wirings), on every panel and
  setting, as E136's and E139's reads hold it at the lab input's;
- RISE and APLPN are E141's rounded x3.24 and x0.155; Dabs is RISE times the family's KC-KC scale on the hemibrain's real
  wiring (E83's `saved_fit('real')`, 26.43: 85.6);
- the readings' seed index of a wiring is its position in the calyx's WIRINGS_OF (real 0, P1 1, P2 2, P3 3), so that
  H's real wiring is pooled and noised as the family's real wiring is.

Design:
- step 0 (`--part step0`): the identity (each calyx's builder at its counted E130 fit reproduces the fit's saved APL-on
  and APL-off runs with 0 differing counts and responses at the present olfsysm; every glomerulus of the megamat,
  validation2 and mixture inputs has boutons in each calyx); the KC-KC matrices beside the hemibrain's shuffle within
  type pairs (E122's SH, E141's arm S); the absolute scale each dose reaches on each calyx beside the hemibrain's eight
  wirings'; the prediction per calyx, written before any read (`write_prediction`);
- step 1 (`--part read --step 1`): the 148 reads (`planned_points`), the families' megamat reads first;
- step 2 (`--part read --step 2`, only for a male calyx on which A1 does not hold): D1.8 and D5.8 on its three
  permutations (24 reads per calyx);
- the score (`--part score`): E141's (E136's classes, E139's price), per calyx, each setting against the calyx's own
  family on the same wiring; H wiring for wiring against the saved reads of the hemibrain's real wiring with 'first'
  claws (E139's and E141's);
- variants (`--part variants`): variants.csv for the model-output archive's params.csv.

Outputs: results/E144_a1_malecns_calyces/ (E144_lock.json, step0_*.csv, step0.json, prediction.json,
reads_<host>_<lane>.csv, provenance_read_<host>_<lane>.json, claims/, score_*.csv, classes.json, variants.csv); runs in
outputs/E144_a1_malecns_calyces/read/ and data/model_outputs/E144_a1_malecns_calyces/read/ (E134's step-3 layout).

Holds (CLAUDE.md): Pratyush's hemibrain and male CNS only (olfsysm's model tables and his KCclaw-to-KCclaw table, through
E130's matrices); no olfsysm change (pn_inhadd and pn_inhsc are options at olfsysm b1e89d7); nothing fitted or
searched; nothing chosen from E138's readings of validation2. Nothing is published.

Concentrations: megamat is at 10^-3 on both sides; validation2 is one panel at matched concentrations; the mixtures are
used as measured; every APL16 number involves paa, euc and -aPine, whose model input is the ORN data within 1 log10
unit of the flies' concentrations.

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E144_a1_malecns_calyces --part read --step 1 --lane <host>[-<n>] [--reverse]
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u -m experiments.E144_a1_malecns_calyces \\
        --part step0
    ./venv/bin/python -u -m experiments.E144_a1_malecns_calyces --part lock | --part score | --part variants
"""
from __future__ import annotations

import argparse
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E144_a1_malecns_calyces'
SCRIPT = 'experiments/E144_a1_malecns_calyces.py'
LOCK_NAME = 'E144_lock.json'
PLANNER = 'a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 10:31'

# the calyces: code -> E130's connectome name
CALYCES: Dict[str, str] = {'H': 'hemibrain-argmax', 'L': 'malecns-left', 'R': 'malecns-right'}
MALE: Tuple[str, ...] = ('L', 'R')
# each calyx's wirings: label -> E133's seed (0 the real wiring; n the within-type claw permutation of seed n)
WIRING_SEEDS: Dict[str, int] = {'real': 0, 'P1': 1, 'P2': 2, 'P3': 3}
WIRINGS_OF: Dict[str, Tuple[str, ...]] = {'H': ('real',), 'L': tuple(WIRING_SEEDS), 'R': tuple(WIRING_SEEDS)}
PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
# E141's rounded factors (E139's K and B): the rise of the KC-KC scale and APL->bouton's factor
RISE = 3.24
APLPN = 0.155
RUNG = 1000.0


@dataclass(frozen=True)
class Setting:
    """One of the entry's settings on a calyx, from the calyx's own family.

    Attributes:
        name: 'family', 'F', 'B', 'D1.8', 'A1', 'D5.8' or 'Dabs'.
        label: what it is, in words (ASCII).
        rung: pn_inhadd of the antennal lobe's normalization (pn_inhsc by E77's rule), or None for olfsysm's defaults.
        aplpn: the factor on the fit's APL->bouton scale.
        kckc: the factor on the fit's KC-KC scale (both conditions), or None where `kckc_abs` sets it.
        kckc_abs: an absolute KC-KC scale ('dabs': `dabs_scale`), or None.
        calyces: the calyces the setting is read on.
    """
    name: str
    label: str
    rung: Optional[float]
    aplpn: float
    kckc: Optional[float]
    kckc_abs: Optional[str] = None
    calyces: Tuple[str, ...] = tuple(CALYCES)


SETTINGS: Dict[str, Setting] = {s.name: s for s in (
    Setting('family', "the calyx's family (E130's fit; olfsysm's pn_inhadd 31.4)", None, 1.0, 1.0),
    Setting('F', 'pn_inhadd 1000 with pn_inhsc by E77\'s rule', RUNG, 1.0, 1.0),
    Setting('B', f'F with APL->bouton at x{APLPN:g}', RUNG, APLPN, 1.0),
    Setting('D1.8', f"B with the calyx's matrix at x{RISE ** 0.5:.3g} of its fitted scale", RUNG, APLPN, RISE ** 0.5),
    Setting('A1', f"B with the calyx's matrix at x{RISE:g} of its fitted scale", RUNG, APLPN, RISE),
    Setting('D5.8', f"B with the calyx's matrix at x{RISE ** 1.5:.3g} of its fitted scale", RUNG, APLPN, RISE ** 1.5),
    Setting('Dabs', "B with L's matrix at the absolute scale A1 reaches on the hemibrain's real wiring", RUNG, APLPN,
            None, kckc_abs='dabs', calyces=('L',)),
)}
DOSES: Tuple[str, ...] = ('B', 'D1.8', 'A1', 'D5.8')


@dataclass(frozen=True)
class Point:
    """A setting read on one calyx's wiring and one panel's input."""
    calyx: str
    setting: str
    wiring: str
    panel: str


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E144_a1_malecns_calyces/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E144_a1_malecns_calyces/."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def local_dir(part: str = '') -> Path:
    """outputs/E144_a1_malecns_calyces/[<part>/] (created)."""
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_dirs() -> List[Path]:
    """Where the reads' runs are: outputs/.../read/, then the archive's read/."""
    return [local_dir('read'), archive_dir() / 'read']


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


# ---------------------------------------------------------------------------
# names

def point_tag(pt: Point) -> str:
    """'<calyx>_<setting>_<panel>_<wiring>', e.g. 'L_A1_megamat_P2' (ValueError for a point outside the design)."""
    check_point(pt)
    return f'{pt.calyx}_{pt.setting}_{pt.panel}_{pt.wiring}'


def parse_tag(tag: str) -> Point:
    """The point of a tag (`point_tag`)."""
    c, s, p, w = str(tag).split('_')
    pt = Point(c, s, w, p)
    check_point(pt)
    return pt


def check_point(pt: Point) -> None:
    """ValueError unless the calyx, setting, wiring and panel are the design's and the setting is read on the calyx."""
    if (pt.calyx not in CALYCES or pt.setting not in SETTINGS or pt.panel not in PANELS
            or pt.wiring not in WIRINGS_OF[pt.calyx] or pt.calyx not in SETTINGS[pt.setting].calyces):
        raise ValueError(f'{pt} is not a point of {EXPERIMENT}')


def variant(calyx: str, setting: str) -> str:
    """A setting on a calyx as a variant name, '<calyx>_<setting>' (e.g. 'L_A1')."""
    if calyx not in CALYCES or setting not in SETTINGS or calyx not in SETTINGS[setting].calyces:
        raise ValueError(f'{calyx=} {setting=}')
    return f'{calyx}_{setting}'


def seed_index(calyx: str, wiring: str) -> int:
    """A wiring's seed index for the stage readings' pooling and noise: its position in WIRINGS_OF (the lock's)."""
    return list(WIRING_SEEDS).index(wiring) if wiring in WIRINGS_OF[calyx] else _bad(calyx, wiring)


def _bad(calyx, wiring):
    raise ValueError(f'{wiring!r} is not a wiring of calyx {calyx!r}')


# ---------------------------------------------------------------------------
# the plan

def planned_points(steps: Sequence[str] = ('1',), calyces: Optional[Sequence[str]] = None) -> List[Point]:
    """The reads of the numbered steps in a forward lane's order. Step 1's 148 (the entry's design item 2) in four
    groups: (1) the families' megamat reads (L's and R's real wirings, then their permutations, then H's), the first
    stage table on a male calyx; (2) B's and A1's megamat reads, in the same order; (3) the family, B and A1 on the
    validation2, mixture and APL16 inputs (L and R on four wirings each, H on its real wiring); (4) F, D1.8 and D5.8 on
    the real wirings of L, R and H on all four panels, and Dabs on L's. Step 2's (only for `calyces`, default
    `step2_calyces`: the male calyces on which A1 does not hold): D1.8 and D5.8 on the calyx's three permutations on all
    four panels (24 per calyx)."""
    out: List[Point] = []
    for s in steps:
        if s == '1':
            out += _step1_points()
        elif s == '2':
            cs = tuple(calyces) if calyces is not None else step2_calyces()
            for c in cs:
                if c not in MALE:
                    raise ValueError(f'step 2 reads the male calyces only, not {c!r}')
                out += [Point(c, st, w, p) for st in ('D1.8', 'D5.8') for p in PANELS
                        for w in WIRINGS_OF[c] if w != 'real']
        else:
            raise ValueError(f'step {s!r}')
    return out


def _wiring_order() -> List[Tuple[str, str]]:
    """(calyx, wiring) in a lane's order: L's and R's real wirings, then their permutations by seed, then H's."""
    perms = [(c, w) for w in WIRINGS_OF['L'] if w != 'real' for c in MALE]
    return [(c, 'real') for c in MALE] + perms + [('H', 'real')]


def _step1_points() -> List[Point]:
    order = _wiring_order()
    pts = [Point(c, 'family', w, 'megamat') for c, w in order]
    pts += [Point(c, s, w, 'megamat') for s in ('B', 'A1') for c, w in order]
    pts += [Point(c, s, w, p) for s in ('family', 'B', 'A1') for p in PANELS if p != 'megamat' for c, w in order]
    pts += [Point(c, s, 'real', p) for s in ('F', 'D1.8', 'D5.8') for c in ('L', 'R', 'H') for p in PANELS]
    pts += [Point('L', 'Dabs', 'real', p) for p in PANELS]
    for pt in pts:
        check_point(pt)
    return pts


def step2_calyces() -> Tuple[str, ...]:
    """Step 2's calyces: the male calyces on which A1 does not hold (classes.json's 'calyces' -> <calyx> -> 'a1_holds'
    False); RuntimeError before step 1's score."""
    p = out_dir() / 'classes.json'
    if not p.exists():
        raise RuntimeError("step 2's calyces wait on step 1's score (classes.json)")
    cl = json.loads(p.read_text()).get('calyces') or {}
    return tuple(c for c in MALE if cl.get(c, {}).get('a1_holds') is False)


# ---------------------------------------------------------------------------
# the calyces' fits, settings and targets

PARAM5: Tuple[str, ...] = ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale')
PARAM_KEYS: Tuple[str, ...] = PARAM5 + ('wKCKC_scale',)
_CACHE: dict = {}


def calyx_fit(calyx: str) -> dict:
    """The calyx's family: E130's counted fit of the calyx (`E130.counted`: of the done starts the lower calibrated
    APL-effect loss). Returns dict(connectome, name (E130's fit name), tag (E10_joint_fit.py's), folder (results of
    the fit, relative to the shared directory), params (absolute: its fits.csv row's five scales and kc_kc_inh as
    wKCKC_scale), saved ({'on', 'off'}: (responses, spike counts)), target_on / target_off (response counts, two or
    more spikes), size (KCs x odors), n_kcs). RuntimeError while the fit is not counted or its outputs are missing."""
    if calyx in _CACHE.setdefault('fit', {}):
        return _CACHE['fit'][calyx]
    from experiments import E130_malecns_family as E130
    con = CALYCES[calyx]
    cnt = E130.counted(con)
    if cnt is None:
        raise RuntimeError(f"E130's fit of {con} is not counted")
    spec, row = cnt['spec'], cnt['row']
    o = E130.load_fit_outputs(spec)
    if o is None:
        raise RuntimeError(f"E130's fit {spec.name} has no saved outputs here or in the archive")
    p = {k: float(row[k]) for k in PARAM5}
    p['wKCKC_scale'] = float(row['kc_kc_inh'])
    n = lambda r: int((r.to_numpy(dtype=float) > 0.5).sum())
    rec = dict(connectome=con, name=spec.name, tag=E130.e10_tag(spec), folder=f'results/{E130.FIT_PREFIX}/{spec.name}',
               params=p, saved=o, target_on=n(o['on'][0]), target_off=n(o['off'][0]), size=int(o['on'][0].size),
               n_kcs=int(o['on'][0].shape[0]))
    _CACHE['fit'][calyx] = rec
    return rec


def pn_terms() -> Dict[str, float]:
    """F's antennal-lobe normalization (E141's `pn_terms`: pn_inhadd RUNG, pn_inhsc by E77's `inhsc_for`)."""
    from experiments import E141_kckc_rise_forms as E141
    return E141.pn_terms()


def hemibrain_real_kckc() -> float:
    """The family's KC-KC scale on the hemibrain's real wiring ('first' claws; E83's `saved_fit('real')`)."""
    from experiments import E83_cross_wiring_fit as E83
    return float(E83.saved_fit('real')['p']['wKCKC_scale'])


def dabs_scale() -> float:
    """Dabs's absolute KC-KC scale: RISE times the family's KC-KC scale on the hemibrain's real wiring (85.6)."""
    return float(RISE * hemibrain_real_kckc())


def setting_params(calyx: str, setting: str) -> Tuple[Dict[str, float], dict]:
    """A setting's start on a calyx and its model kwargs beside the calyx's: the calyx's fit's absolute parameters
    (`calyx_fit`) with APL->bouton times the setting's factor and the KC-KC scale times its factor (or at its absolute
    scale, Dabs's `dabs_scale`); the kwargs `pn_terms()` where the setting weakens the normalization, else {}. The
    threshold and the tied APL scales are then held (`hold_point`)."""
    st = SETTINGS[setting]
    if calyx not in st.calyces:
        raise ValueError(f'{setting} is not read on {calyx}')
    p = dict(calyx_fit(calyx)['params'])
    p['wAPLPN_scale'] = p['wAPLPN_scale'] * float(st.aplpn)
    if st.kckc_abs == 'dabs':
        p['wKCKC_scale'] = dabs_scale()
    elif st.kckc_abs is not None:
        raise ValueError(f'{st.kckc_abs=}')
    else:
        p['wKCKC_scale'] = p['wKCKC_scale'] * float(st.kckc)
    return p, (dict(pn_terms()) if st.rung is not None else {})


def read_targets(calyx: str, panel: str, n_kcs: int, n_odors: int, lock: Optional[Mapping] = None
                 ) -> Tuple[int, int, float, float]:
    """(target_off, target_on, fraction_off, fraction_on) of a read: APL on at E136's operating point of the panel
    (megamat E134's held rate, the lock's Remy criterion 0.0915; validation2 the lock's rate of its flies, 0.050; the
    mixtures E134's held mixture rate, 0.054), APL off at that times the calyx's fit's APL-off / APL-on response counts;
    the APL16 input at the calyx's fit's two counts as fractions of its size (a permutation keeps every KC, so the
    counts themselves). Counts are rounded fractions of `n_kcs` x `n_odors`."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    fit = calyx_fit(calyx)
    size = int(n_kcs) * int(n_odors)
    if panel == 'apl16':
        f_on, f_off = fit['target_on'] / fit['size'], fit['target_off'] / fit['size']
        return int(round(f_off * size)), int(round(f_on * size)), float(f_off), float(f_on)
    lk = lock if lock is not None else R.require_lock()
    if panel == 'megamat':
        f_on = float(B.held_rate('megamat', lk))
    elif panel == 'validation2':
        f_on = float(lk['populations']['kc_soma__validation2__wild_type']['response_rate']['mean_over_flies'])
    elif panel == 'mixtures':
        f_on = float(B.held_rate('mixtures', lk))
    else:
        raise ValueError(f'{panel=}')
    f_off = f_on * fit['target_off'] / fit['target_on']
    return int(round(f_off * size)), int(round(f_on * size)), float(f_off), float(f_on)


# ---------------------------------------------------------------------------
# the builder and the read

def tables(calyx: str, wiring: str) -> dict:
    """A calyx's tables on one of its wirings: E133's `calyx_tables` at the wiring's seed (the lab's APL16 input, the
    calyx's base, normalized APL weights, the claw x bouton wiring, its APL-off breadth, the male calyx's KC-KC matrix
    checked against E130's matrices.json; cached in E133's process cache)."""
    from experiments import E133_calyx_transplants as E133
    if wiring not in WIRINGS_OF[calyx] and not (calyx == 'H' and wiring in WIRING_SEEDS):
        raise ValueError(f'{wiring!r} is not a wiring of {calyx!r}')
    return E133.calyx_tables(CALYCES[calyx], WIRING_SEEDS[wiring])


def drop_tables(calyx: str, wiring: str) -> None:
    """Frees a wiring's cached tables (E133's `drop_tables`)."""
    from experiments import E133_calyx_transplants as E133
    E133.drop_tables(CALYCES[calyx], WIRING_SEEDS[wiring])


def build_simulator(calyx: str, wiring: str, params: Mapping[str, float], orn: Optional[pd.DataFrame] = None,
                    extra_kws: Optional[Mapping] = None) -> Tuple:
    """The calyx's family on a wiring as a fixed-parameter simulator at `params`: E133's `build_simulator` (E130's on
    given tables: the connectome's kwargs, the wiring, the breadth-weighted APL->KC weights at exponent 2, the other
    normalized APL weights, the male calyx's KC-KC matrix as `_wKCKC`) with `orn` in place of the tables' APL16 input
    and `extra_kws` (the setting's normalization terms) beside the connectome's kwargs; the breadth stays the tables'
    (the APL16 input's). With neither, the call is E133's exactly (the test checks the arguments). Returns (simulate,
    info): 'params', 'kws' (scalar model kwargs), 'connectome', 'multibouton_claws', 'breadth', 'wAPLKC', 'W', 'orn',
    'cells' (the simulator's closure cells: 'mp', 'rv', 'odor_order', for E134's `Built`)."""
    from experiments import E10_joint_fit as J
    from experiments import E133_calyx_transplants as E133
    from experiments.E77_input_screen import _scalar_kws
    from pyolfsysm.mb_model import make_fixed_param_kc_simulator
    T = tables(calyx, wiring)
    con = CALYCES[calyx]
    p = dict(params)
    orn = T['orn'] if orn is None else orn
    wnorm = T['wnorm']
    wAPLKC = J.breadth_weights(wnorm['wAPLKC'], T['breadth'], E133.K)
    extra = {} if T.get('kckc') is None else {'_wKCKC': T['kckc']}
    kws = dict(J.base_kws(None, connectome=con), **dict(extra_kws or {}))
    sim, _ = make_fixed_param_kc_simulator(orn, initial_params=p, **kws, **extra, _wPNKC=T['W'], _wAPLKC=wAPLKC,
                                           _wKCAPL=wnorm['wKCAPL'], _wAPLPN=wnorm['wAPLPN'], _wPNAPL=wnorm['wPNAPL'],
                                           kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)
    cells = dict(zip(sim.__code__.co_freevars, sim.__closure__))
    return sim, dict(params=p, kws=_scalar_kws(dict(kws, kc_kc_inh=p['wKCKC_scale'], kc_kc_tau=J.TAU)), connectome=con,
                     multibouton_claws=kws['multibouton_claws'], breadth=T['breadth'], wAPLKC=wAPLKC, W=T['W'], orn=orn,
                     cells=cells)


def run_wiring(calyx: str, wiring: str) -> str:
    """The wiring as a run records it ('<calyx>-<wiring>', e.g. 'L-P2'), so that no male run is named as a hemibrain
    wiring of the family's."""
    return f'{calyx}-{wiring}'


def build(calyx: str, wiring: str, panel: str, params: Mapping[str, float], extra_kws: Optional[Mapping] = None,
          orn: Optional[pd.DataFrame] = None):
    """E134's `Built` of the calyx's wiring on a panel's input (E134's `panel_input`, or `orn`) at `params`: the
    fixed-parameter simulator of `build_simulator` with its olfsysm variables, the wiring, the PN stage's glomeruli and
    the units' metadata, so that E134's `run_stages` saves every stage (PN, bouton, claw, soma, APL)."""
    from experiments import E134_population_baseline as B
    from kcq import stage_readings as R
    orn = B.panel_input(panel) if orn is None else orn
    sim, info = build_simulator(calyx, wiring, params, orn, extra_kws)
    T = tables(calyx, wiring)
    return B._simulator_built(sim, info['cells'], info['W'], dict(params), orn, R.STAGE_FILES, dict(
        model=B.FAMILY, wiring=run_wiring(calyx, wiring), calyx=calyx, connectome=CALYCES[calyx], calyx_wiring=wiring,
        wiring_seed=WIRING_SEEDS[wiring], kws=info['kws'], W_sha256=T['meta'].get('W_sha256'),
        breadth_sha256=T['meta'].get('breadth_sha256'), kckc_sha256=T['meta'].get('kckc_sha256'),
        breadth="the calyx wiring's APL-off breadth on the lab's APL16 input (E133's calyx_tables)",
        **B._override_record(orn, extra_kws)))


def hold_point(pt: Point, lock: Optional[Mapping] = None) -> dict:
    """Both response rates of a point held exactly at its operating point (`read_targets`), from the setting's start
    (`setting_params`): E77's `hold_rates` on `build_simulator` with the panel's input and the setting's kwargs, the
    threshold bisected with APL off (APL->KC and APL->bouton 0), then the APL->KC scale with KC->APL and bouton->APL
    tied with APL on, within E136's APL bracket (HOLD_RANGE); KC-KC and APL->bouton as set. Returns dict(hold:
    hold_rates' dict, start, kws, orn, targets: (target_off, target_on, fraction_off, fraction_on), n_kcs)."""
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    from experiments.E77_input_screen import hold_rates
    from kcq import stage_readings as R
    check_point(pt)
    lk = lock if lock is not None else R.require_lock()
    p0, kws = setting_params(pt.calyx, pt.setting)
    orn = B.panel_input(pt.panel)
    T = tables(pt.calyx, pt.wiring)
    n_kcs = int(T['breadth'].shape[0])
    if n_kcs != calyx_fit(pt.calyx)['n_kcs']:
        raise RuntimeError(f"{point_tag(pt)}: {n_kcs} KCs on the wiring, {calyx_fit(pt.calyx)['n_kcs']} in the fit")
    targets = read_targets(pt.calyx, pt.panel, n_kcs, orn.shape[1], lk)
    hold = hold_rates(lambda: build_simulator(pt.calyx, pt.wiring, dict(p0), orn, kws or None),
                      target_off=targets[0], target_on=targets[1], params=dict(p0), label=f'E144 {point_tag(pt)}',
                      apl_bracket=F.HOLD_RANGE)
    return dict(hold=hold, start=p0, kws=kws, orn=orn, targets=tuple(targets), n_kcs=n_kcs)


def read_extra(pt: Point) -> dict:
    """The fields a read records beside the hold's (in its row and its run's info)."""
    st = SETTINGS[pt.setting]
    return dict(variant=variant(pt.calyx, pt.setting), calyx=pt.calyx, connectome=CALYCES[pt.calyx],
                setting=pt.setting, calyx_wiring=pt.wiring, wiring_seed=WIRING_SEEDS[pt.wiring],
                seed_index=seed_index(pt.calyx, pt.wiring), rung=st.rung, aplpn_factor=st.aplpn,
                kckc_factor=st.kckc, kckc_abs=st.kckc_abs, fit=calyx_fit(pt.calyx)['name'])


def _append_row(path: Path, row: Mapping) -> None:
    from experiments import E136_stage_population_fit as F
    F._append_row(path, row)


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None,
               dirs: Optional[Sequence[Path]] = None, root: Optional[Path] = None) -> Optional[dict]:
    """One read, claimed by mkdir of <root>/claims/<tag>/ (E134's `claim`; root default out_dir(); None if done or
    claimed elsewhere): held (`hold_point`), built again at the held parameters (`build`) and run with every stage
    saved (E134's `run_stages`, E90's subtractive claw reader), which must give the hold's outputs exactly; saved in
    each of `dirs` (default `read_dirs`; E134's step-3 layout, `save_step3_run`), its row in
    <root>/reads_<host>_<lane>.csv and its provenance (`H.model_provenance` with the input) in
    <root>/provenance_read_<host>_<lane>.json. The claim is released on failure."""
    import gc
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    from experiments.E77_input_screen import frame_sha256
    from kcq import stage_readings as R
    tag = point_tag(pt)
    root = out_dir() if root is None else Path(root)
    dirs = read_dirs() if dirs is None else list(dirs)
    if done is not None and tag in done:
        return None
    if not B.claim(tag, lane, reclaim, root=root):
        print(f'E144 {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        lk = R.require_lock()
        h = hold_point(pt, lk)
        hold, kws, orn = h['hold'], h['kws'], h['orn']
        t_off, t_on, f_off, f_on = h['targets']
        p, poff = dict(hold['params']), dict(hold['params_off'])
        built = build(pt.calyx, pt.wiring, pt.panel, p, kws or None, orn=orn)
        out = B.run_stages(B.Point(pt.panel, 'read', run_wiring(pt.calyx, pt.wiring)), built, p, poff)
        for c in ('on', 'off'):
            if B._differ(out[c]['spike_counts'], hold[c][1]) or B._differ(out[c]['responses'], hold[c][0]):
                raise RuntimeError(f"{tag} {c}: the saved run differs from the hold's outputs")
        extra = read_extra(pt)
        info = dict(model=B.FAMILY, wiring=run_wiring(pt.calyx, pt.wiring), panel=pt.panel, reading='read',
                    arm=extra['variant'], params=p, params_off=poff, start=h['start'], pn_kws=kws, tag=tag,
                    targets=dict(on=t_on, off=t_off, frac_on=f_on, frac_off=f_off), host=_host(), lane=lane,
                    time=_now(), **extra)
        R.save_step3_run(tag, out, built.meta, built.orn, info, dirs)
        rr = {}
        for c in ('on', 'off'):
            r = out[c]['responses'].to_numpy(dtype=float) > 0.5
            rr.update({f'n_{c}': int(r.sum()), f'rate_{c}': float(r.mean()),
                       f'max_count_{c}': float(out[c]['spike_counts'].to_numpy(dtype=float).max())})
        row = dict(tag=tag, panel=pt.panel, wiring=pt.wiring, half='held' if pt.panel == 'megamat' else None, **extra,
                   **{k: float(v) for k, v in p.items()}, **{f'off_{k}': float(v) for k, v in poff.items()},
                   **{f'start_{k}': float(v) for k, v in h['start'].items()},
                   **{k: float(v) for k, v in kws.items()}, **rr,
                   target_on=t_on, target_off=t_off, frac_on_target=f_on, frac_off_target=f_off,
                   thr_factor=hold['thr_factor'], apl_factor=hold['apl_factor'], at_bound_thr=hold['at_bound_thr'],
                   at_bound_apl=hold['at_bound_apl'], hold_n_sims=hold['n_sims'], seconds=round(time.time() - t0, 1),
                   host=_host(), lane=lane, time=_now(), threads=os.environ.get('OMP_NUM_THREADS'),
                   input_sha256=frame_sha256(built.orn), **F._commits())
        _append_row(root / f'reads_{_host()}_{lane}.csv', row)
        prov_path = root / f'provenance_read_{_host()}_{lane}.json'
        prov = _read_json(prov_path) or dict(part='read', lane=lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(info={k: v for k, v in built.info.items() if not str(k).startswith('_')},
                                   params=p, params_off=poff, start=h['start'], kws=kws, hold_kws=hold.get('kws'),
                                   history=hold.get('history'), model=H.model_provenance(orn_deltas=built.orn))
        _write_json(prov_path, prov)
        print(f"E144 {tag}: {row['seconds']:.0f} s, rate on {rr['rate_on']:.4f} / off {rr['rate_off']:.4f} "
              f"(targets {f_on:.4f} / {f_off:.4f}), thr x{hold['thr_factor']:.4f}, APL x{hold['apl_factor']:.4f}",
              flush=True)
        del out, built, h, hold
        gc.collect()
        return row
    except BaseException:
        B.release(tag, root=root)
        raise


def reads_table() -> pd.DataFrame:
    """Every lane's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def clean_commits() -> dict:
    """`H.code_commits()` at the start of a lane: RuntimeError unless every checkout is committed and clean
    (queue_dirty included: start lanes from a clean checkout or worktree), so that every read's provenance names the
    code that ran."""
    cc = H.code_commits()
    dirty = [k for k, v in cc.items() if k.endswith('_dirty') and v]
    if dirty:
        raise RuntimeError(f'commit first or run from a clean checkout: {dirty}')
    return cc


def still_clean(cc: Mapping) -> None:
    """Before each read: no file this process has loaded has changed since the lane began (E108's `_still_clean`)."""
    from experiments import E108_paa_extrapolation as E108
    E108._still_clean(dict(cc))


def part_read(lane: str, steps: Sequence[str] = ('1',), reverse: bool = False, reclaim: bool = False,
              calyces: Optional[Sequence[str]] = None) -> None:
    """A lane of reads: every planned point of `steps` (`planned_points`; step 2's `calyces`), reversed with `reverse`,
    skipping what is done or claimed elsewhere; refuses unless the lock is in (`require_lock`), step 0's identity has
    passed (step0.json) and the checkouts are clean (`clean_commits`)."""
    import gc
    require_lock()
    s0 = _read_json(out_dir() / 'step0.json')
    if not (s0.get('identity_passes') and s0.get('glomeruli_pass')):
        raise RuntimeError("step 0's identity and glomerulus checks have not passed (step0.json): run --part step0")
    cc = clean_commits()
    pts = planned_points(steps, calyces)
    for pt in (pts[::-1] if reverse else pts):
        still_clean(cc)
        read_point(pt, lane, reclaim, set(reads_table().tag))
        gc.collect()


# ---------------------------------------------------------------------------
# step 0

def step0_identity(calyces: Sequence[str] = tuple(CALYCES)) -> pd.DataFrame:
    """Identity: each calyx's counted E130 fit reproduced by this module's read path at the present olfsysm: `build` on
    the real wiring and the APL16 input at the fit's absolute parameters, both conditions run by E134's `run_stages`
    (APL off: APL->KC and APL->bouton 0, as E10_joint_fit.py runs it), against the fit's saved APL-on and APL-off
    responses and spike counts (0 differing entries, the same KCs and odors). One row per calyx and condition, with the
    olfsysm build; a stage reading of the APL16 run (atlas's `read_run` on the megamat odors it holds) checks that the
    stage readings work on the calyx's wiring (its populations' sizes)."""
    import gc
    from experiments import E133_calyx_transplants as E133
    from experiments import E134_population_baseline as B
    from experiments import E136_stage_population_fit as F
    from kcq import stage_readings as R
    rows = []
    for c in calyces:
        t0 = time.time()
        fit = calyx_fit(c)
        p = dict(fit['params'])
        poff = dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
        built = build(c, 'real', 'apl16', p)
        out = B.run_stages(B.Point('apl16', 'identity', run_wiring(c, 'real')), built, p, poff)
        run = F.in_memory_run(f'{c}_identity', out, built, dict(model=B.FAMILY, wiring=run_wiring(c, 'real'),
                                                                 panel='apl16'))
        try:
            rd = R.read_run(run, seed_index=seed_index(c, 'real'))
            pops = {f'{k[0]}/{k[1]}/{k[2]}': int(v.values.shape[0]) for k, v in rd.populations.items()}
            read_ok, read_err = bool(pops), None
        except Exception as e:  # reported, not hidden: the reads need the stage readings on this wiring
            pops, read_ok, read_err = {}, False, f'{type(e).__name__}: {e}'
        for cond, k in (('on', 'on'), ('off', 'off')):
            r_ref, s_ref = fit['saved'][k]
            rows.append(dict(calyx=c, connectome=CALYCES[c], fit=fit['name'], tag=fit['tag'], condition=cond,
                             counts_differ=B._differ(out[cond]['spike_counts'], s_ref),
                             responses_differ=B._differ(out[cond]['responses'], r_ref),
                             same_columns=bool(list(out[cond]['spike_counts'].columns) == list(s_ref.columns)),
                             n_kcs=int(s_ref.shape[0]), n_odors=int(s_ref.shape[1]),
                             stages_saved=','.join(sorted(k2 for k2 in out[cond] if k2.endswith('_win'))),
                             stage_reading_ok=read_ok, stage_reading_error=read_err,
                             stage_populations=json.dumps(pops, sort_keys=True), seconds=round(time.time() - t0, 1),
                             host=_host(), **E133.olfsysm_build()))
        print(f"E144 step 0 identity {c}: " + ', '.join(f"{r['condition']} {r['counts_differ']} counts / "
                                                        f"{r['responses_differ']} responses differ"
                                                        for r in rows[-2:]) + f'; stage reading {read_ok}', flush=True)
        del out, built, run
        gc.collect()
    df = pd.DataFrame(rows)
    df['passes'] = (df['counts_differ'] == 0) & (df['responses_differ'] == 0) & df['same_columns']
    return df


def step0_glomeruli() -> pd.DataFrame:
    """Every glomerulus of the megamat, validation2, mixture and APL16 inputs (E134's `panel_input`: the lab input's
    rows) against each calyx's wiring (E130's `step0_wiring_only`; a claw permutation keeps the boutons): its boutons and
    claws per calyx, and whether it is in the input of each panel; the wiring's glomeruli without input beside (rows
    with in_input False). 'has_boutons' must be True for every input glomerulus on every calyx."""
    from experiments import E130_malecns_family as E130
    from experiments import E134_population_baseline as B
    inputs = {p: [str(g) for g in B.panel_input(p).index] for p in PANELS}
    allg = list(dict.fromkeys(g for p in PANELS for g in inputs[p]))
    rows = []
    for c, con in CALYCES.items():
        W = E130.step0_wiring_only(con)
        bg = pd.Series([str(x) for x in W.columns.get_level_values('glomerulus')])
        A = W.to_numpy(dtype=float) > 0
        claws = pd.Series(A.sum(axis=0), index=bg.to_numpy()).groupby(level=0).sum()
        nb = bg.value_counts()
        for g in list(dict.fromkeys(allg + sorted(set(bg) - set(allg)))):
            rows.append(dict(calyx=c, connectome=con, glomerulus=g, in_input=g in allg,
                             **{f'in_{p}': g in inputs[p] for p in PANELS}, n_boutons=int(nb.get(g, 0)),
                             n_claws=int(claws.get(g, 0)), has_boutons=bool(nb.get(g, 0) > 0)))
    return pd.DataFrame(rows)


def _kc_glomeruli(W: pd.DataFrame) -> Tuple[pd.Index, np.ndarray]:
    """The wiring's KC IDs (first appearance) and a KC x glomerulus 0/1 matrix of the glomeruli each KC's claws take
    input from (any positive weight)."""
    kc = W.index.get_level_values(H.KC_ID_LEVEL)
    ids = pd.Index(pd.unique(kc))
    g = pd.Index(list(dict.fromkeys(str(x) for x in W.columns.get_level_values('glomerulus'))))
    G = pd.DataFrame((W.to_numpy(dtype=float) > 0).astype(float), index=kc,
                     columns=[str(x) for x in W.columns.get_level_values('glomerulus')])
    G = G.T.groupby(level=0).max().T.groupby(level=0).max()
    return ids, (G.reindex(index=ids, columns=g).fillna(0.0).to_numpy() > 0)


def matrix_table(name: str, M: pd.DataFrame, W: pd.DataFrame, types: np.ndarray) -> dict:
    """One KC x KC matrix (post x pre, normalized to a mean row sum of 1) on its wiring `W` and its KCs' `types`:
    E122's `matrix_stats` (entries, inputs per KC, KCs without input, the CV of the KCs' summed input, the share of
    reciprocal pairs among nonzero entries) against itself, and beside it: the Gini coefficient and the top-10 % share of
    the KCs' summed input (unevenness); the weighted reciprocity (sum of min(M, M^T) over the total); the weight between
    KCs that share an input glomerulus (the share of pairs that share one, the share of the total weight on them, and
    the ratio of the two: 1 is no preference); and the weight by KC type pair (E122's `block_table`: each (post, pre)
    block's share of the total, and the share within type beside the share of pairs within type)."""
    from experiments import E122_feedback_neuron as E122
    A = M.to_numpy(dtype=float)
    n = len(A)
    ids, G = _kc_glomeruli(W)
    order = ids.get_indexer(M.index)
    if (order < 0).any():
        raise RuntimeError(f"{name}: the matrix's KCs are not the wiring's")
    G = G[order]
    S = (G.astype(float) @ G.T.astype(float)) > 0
    off = ~np.eye(n, dtype=bool)
    st = E122.matrix_stats(A, A)
    st.pop('row_sum_r_hemibrain', None)
    rs = np.sort(A.sum(axis=1))
    cum = np.cumsum(rs)
    gini = float((n + 1 - 2 * (cum.sum() / cum[-1])) / n) if cum[-1] > 0 else np.nan
    top = float(rs[-max(1, n // 10):].sum() / rs.sum())
    share_pairs = float(S[off].mean())
    share_w = float(A[S & off].sum() / A[off].sum())
    blk = E122.block_table(A, types)
    same = types[:, None] == types[None, :]
    row = dict(matrix=name, n_kcs=n, **st, row_sum_gini=gini, row_sum_top10_share=top,
               reciprocity_weighted=float(np.minimum(A, A.T)[off].sum() / A[off].sum()),
               shared_glomerulus_pair_share=share_pairs, shared_glomerulus_weight_share=share_w,
               shared_glomerulus_ratio=share_w / share_pairs,
               within_type_weight_share=float(A[same & off].sum() / A[off].sum()),
               within_type_pair_share=float((same & off).sum() / off.sum()))
    for r in blk.itertuples():
        row[f'block_{_tcol(r.post_type)}_from_{_tcol(r.pre_type)}'] = float(r.weight_share)
    return row


def _tcol(t: str) -> str:
    return {'g': 'g', "a'b'": 'apbp', 'ab': 'ab'}[t]


def step0_matrices() -> pd.DataFrame:
    """The KC-KC matrices (normalized, mean row sum 1), one row each (`matrix_table`): the hemibrain's on the family's
    KCs and wiring ('first' claws; 'hemibrain'), its shuffle within type pairs (E122's SH, E141's arm S; 'hemibrain_SH'),
    and the three calyces' own (H: olfsysm's hemibrain matrix among H's KCs; L and R: E130's male matrices)."""
    from experiments import E122_feedback_neuron as E122
    from experiments import E130_malecns_family as E130
    from pyolfsysm.mb_model import connectome_KCKC_weights
    rows = []
    Wf = E130.step0_wiring_only(E130.FIRST)
    idx = E122.model_kc_index()
    ref = E122.hemibrain_kckc(idx.get_level_values(H.KC_ID_LEVEL))
    tf = E122.kc_types_of(idx)
    rows.append(dict(matrix_table('hemibrain', ref, Wf, tf), calyx='hemibrain (first)'))
    rows.append(dict(matrix_table('hemibrain_SH', E122.arm_matrix('SH', ref, tf), Wf, tf),
                     calyx='hemibrain (first), shuffled within type pairs'))
    for c, con in CALYCES.items():
        W = E130.step0_wiring_only(con)
        ids = E130.model_kc_ids(W)
        M = connectome_KCKC_weights(pd.Index(ids)) if c == 'H' else E130.malecns_kckc(con, ids)
        types = E130.kc_types_of_wiring(W).reindex(M.index).to_numpy(dtype=object)
        rows.append(dict(matrix_table(c, M, W, types), calyx=con))
    return pd.DataFrame(rows)


def step0_scales() -> pd.DataFrame:
    """The absolute scales each setting reaches before the hold (`setting_params`; the threshold and the tied APL scales
    are then held): per calyx and setting, the KC-KC scale and APL->bouton scale with their factors on the fit's; the
    hemibrain's eight wirings beside (E83's `saved_fit`, 'first' claws: the family's own scales times the doses' factors,
    as E136's A1 and E141's doses apply them)."""
    from experiments import E83_cross_wiring_fit as E83
    rows = []
    for c in CALYCES:
        fp = calyx_fit(c)['params']
        for s, st in SETTINGS.items():
            if c not in st.calyces:
                continue
            p, kws = setting_params(c, s)
            rows.append(dict(calyx=c, wiring='real (and its permutations)' if c in MALE else 'real', setting=s,
                             kckc_fit=fp['wKCKC_scale'], kckc=p['wKCKC_scale'], kckc_factor=p['wKCKC_scale'] / fp['wKCKC_scale'],
                             kckc_over_hemibrain_real=p['wKCKC_scale'] / hemibrain_real_kckc(),
                             aplpn_fit=fp['wAPLPN_scale'], aplpn=p['wAPLPN_scale'],
                             pn_inhadd=kws.get('pn_inhadd'), pn_inhsc=kws.get('pn_inhsc'), source=calyx_fit(c)['name']))
    from experiments import E136_stage_population_fit as F
    for w in F.WIRINGS:
        fp = E83.saved_fit(w)['p']
        for s, st in SETTINGS.items():
            if st.kckc is None:
                continue
            rows.append(dict(calyx='hemibrain family (first)', wiring=w, setting=s, kckc_fit=fp['wKCKC_scale'],
                             kckc=fp['wKCKC_scale'] * st.kckc, kckc_factor=st.kckc,
                             kckc_over_hemibrain_real=fp['wKCKC_scale'] * st.kckc / hemibrain_real_kckc(),
                             aplpn_fit=fp['wAPLPN_scale'], aplpn=fp['wAPLPN_scale'] * st.aplpn,
                             pn_inhadd=pn_terms()['pn_inhadd'] if st.rung else None,
                             pn_inhsc=pn_terms()['pn_inhsc'] if st.rung else None, source=f"E83 saved_fit('{w}')"))
    return pd.DataFrame(rows)


# the entry's numbers the prediction rests on (its "How it moves" and "Why" paragraphs, the planner's arithmetic)
PREDICTION_CONTEXT: Dict[str, object] = {
    'hemibrain_heldout': dict(measure='strength by odor, one-fly units, megamat stage error, T wirings against the '
                                      'held-out half', family=1.049, A1=0.934, bar=-0.020),
    'dendrites_fall': dict(family=1.012, A1=0.858, contribution=0.035,
                           note='the same on every calyx; above the bar by itself (the planner\'s arithmetic)'),
    'test': 'after the dendrites: the boutons, claws and somata, and the held-out panels\' somata',
    'apl_effect_pattern_parts_e130': dict(H=0.0586, L=0.0715, R=0.0608, family=0.056, floor=0.029,
                                          a1_price_hemibrain=0.0164),
    'bouton_share_of_apl_synapses': dict(H=0.128, L=0.126, R=0.119),
    'e142_step0_bouton_count_weighting': dict(L=1.088, R=1.148, unweighted=1.052, hemibrain=0.997,
                                              reading='male counts read worse than unweighted: a prediction for '
                                                      'the boutons'),
    'kckc_scales': dict(L=7.7, R=25.8, hemibrain=26.4, dabs=85.6, L_at_a1=25.0),
    'family_outside_hemibrain_wirings': 'about 1 in 2 at the boutons (outcome (e))',
}


def write_prediction() -> dict:
    """The planner's prediction per calyx (PREDICTION, the entry's table) with the numbers it rests on
    (PREDICTION_CONTEXT) and the outcomes' odds (OUTCOMES), into prediction.json before any read: kept as first
    written; 'reads_before' records how many reads existed when it was written (0 expected)."""
    path = out_dir() / 'prediction.json'
    cur = _read_json(path)
    if cur:
        return cur
    rec = dict(experiment=EXPERIMENT, planner=PLANNER, prediction=PREDICTION, context=PREDICTION_CONTEXT,
               outcomes=OUTCOMES, reads_before=int(len(reads_table())), written=_now(), host=_host(),
               commits=H.code_commits())
    _write_json(path, rec)
    return rec


def _f(x, sp: str = '.3f') -> str:
    try:
        return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(float(x), sp)
    except (TypeError, ValueError):
        return str(x)


MATRIX_ROWS: Tuple[Tuple[str, str, str], ...] = (
    ('n_kcs', 'KCs', '.0f'), ('inputs_per_kc', 'KC-KC inputs per KC (nonzero entries)', '.2f'),
    ('kcs_without_input', 'KCs with no KC-KC input', '.0f'),
    ('row_sum_cv', "unevenness: CV of the KCs' summed input", '.3f'),
    ('row_sum_gini', "unevenness: Gini of the KCs' summed input", '.3f'),
    ('row_sum_top10_share', "unevenness: the top 10 % of KCs' share of the summed input", '.3f'),
    ('shared_glomerulus_pair_share', 'KC pairs that share an input glomerulus (share of pairs)', '.3f'),
    ('shared_glomerulus_weight_share', 'weight between KCs that share an input glomerulus (share of the total)', '.3f'),
    ('shared_glomerulus_ratio', 'the same, over the share of pairs (1: no preference)', '.3f'),
    ('reciprocal_share', 'reciprocity: nonzero entries whose transpose is nonzero', '.3f'),
    ('reciprocity_weighted', 'reciprocity: sum of min(M, transpose) over the total', '.3f'),
    ('within_type_weight_share', 'weight within KC type (share of the total)', '.3f'),
    ('within_type_pair_share', 'KC pairs within type (share of pairs)', '.3f'),
) + tuple((f'block_{a}_from_{b}', f'weight from {b} to {a} KCs (share of the total)', '.3f')
          for a in ('g', 'apbp', 'ab') for b in ('g', 'apbp', 'ab'))
MATRIX_NAMES = {'hemibrain': "hemibrain ('first' KCs)", 'hemibrain_SH': 'hemibrain shuffled within type pairs',
                'H': "H (hemibrain, 'argmax')", 'L': 'L (male CNS left)', 'R': 'R (male CNS right)'}


def step0_table(res: Mapping[str, pd.DataFrame]) -> str:
    """Step 0's tables as markdown: the identity, the input glomeruli's boutons, the matrices (statistic down the first
    column, the matrices across) and the absolute KC-KC scales per setting (the calyces beside the hemibrain's eight
    wirings' range)."""
    L: List[str] = ['# E144 step 0', '']
    if 'identity' in res:
        d = res['identity']
        L += ['## Identity against E130\'s counted fits (APL16 input, real wiring)', '',
              '| calyx | fit | condition | counts differing | responses differing | stage reading |',
              '|---|---|---|---:|---:|---|']
        L += [f'| {r.calyx} | {r.fit} | {r.condition} | {r.counts_differ} | {r.responses_differ} | '
              f'{"ok" if r.stage_reading_ok else r.stage_reading_error} |' for r in d.itertuples()]
        L.append('')
    if 'glomeruli' in res:
        g = res['glomeruli']
        g = g[g['in_input']]
        L += ['## The input glomeruli\'s boutons per calyx', '',
              '| calyx | input glomeruli | without boutons | boutons per input glomerulus, min / median / max |',
              '|---|---:|---:|---:|']
        for c, x in g.groupby('calyx'):
            L.append(f"| {c} | {len(x)} | {int((~x['has_boutons']).sum())} | {int(x.n_boutons.min())} / "
                     f"{x.n_boutons.median():g} / {int(x.n_boutons.max())} |")
        L.append('')
    if 'matrices' in res:
        m = res['matrices'].set_index('matrix')
        cols = [k for k in MATRIX_NAMES if k in m.index]
        L += ['## The KC-KC matrices (normalized to a mean row sum of 1)', '',
              '| statistic | ' + ' | '.join(MATRIX_NAMES[k] for k in cols) + ' |', '|---|' + '---:|' * len(cols)]
        for col, lab, sp in MATRIX_ROWS:
            if col in m.columns:
                L.append(f'| {lab} | ' + ' | '.join(_f(m.loc[k, col], sp) for k in cols) + ' |')
        L.append('')
    if 'scales' in res:
        sc = res['scales']
        own = sc[sc['calyx'].isin(list(CALYCES))]
        hb = sc[~sc['calyx'].isin(list(CALYCES))]
        L += ["## The absolute KC-KC scale each setting reaches (before the hold; APL->bouton beside)", '',
              "| setting | H | L | R | the hemibrain's eight wirings (min to max) | APL->bouton H / L / R |",
              '|---|---:|---:|---:|---:|---:|']
        for s in SETTINGS:
            v = {c: own[(own.calyx == c) & (own.setting == s)] for c in CALYCES}
            hh = hb[hb.setting == s]
            rng = f'{hh.kckc.min():.1f} to {hh.kckc.max():.1f}' if len(hh) else ''
            L.append(f'| {s} | ' + ' | '.join(_f(v[c].kckc.iloc[0], '.1f') if len(v[c]) else '' for c in CALYCES)
                     + f' | {rng} | ' + ' / '.join(_f(v[c].aplpn.iloc[0], '.3f') if len(v[c]) else '-'
                                                    for c in CALYCES) + ' |')
        L.append('')
    return '\n'.join(L) + '\n'


def write_step0_table() -> str:
    """step0_table.md from the step0_<part>.csv files in the results folder."""
    res = {k: pd.read_csv(out_dir() / f'step0_{k}.csv') for k in ('identity', 'glomeruli', 'matrices', 'scales')
           if (out_dir() / f'step0_{k}.csv').exists()}
    t = step0_table(res)
    (out_dir() / 'step0_table.md').write_text(t)
    return t


def part_step0(parts: Sequence[str] = ('glomeruli', 'matrices', 'scales', 'identity')) -> dict:
    """Step 0 (design item 1): the prediction first (`write_prediction`), then step0_<part>.csv for each of `parts`
    (identity, glomeruli, matrices, scales), and step0.json (identity_passes, glomeruli_pass and the tables' files),
    which the read lanes require."""
    require_lock()
    pred = write_prediction()
    fns = dict(glomeruli=step0_glomeruli, matrices=step0_matrices, scales=step0_scales, identity=step0_identity)
    rec = _read_json(out_dir() / 'step0.json')
    for k in parts:
        df = fns[k]()
        df.to_csv(out_dir() / f'step0_{k}.csv', index=False)
        print(f'E144 step 0: step0_{k}.csv ({len(df)} rows)', flush=True)
        if k == 'identity':
            rec['identity_passes'] = bool(df['passes'].all())
            rec['stage_reading_ok'] = bool(df['stage_reading_ok'].all())
            rec['identity'] = json.loads(df[['calyx', 'condition', 'counts_differ', 'responses_differ',
                                             'passes']].to_json(orient='records'))
        if k == 'glomeruli':
            g = df[df['in_input']]
            rec['glomeruli_pass'] = bool(g['has_boutons'].all())
            rec['glomeruli_without_boutons'] = json.loads(g.loc[~g['has_boutons'], ['calyx', 'glomerulus']]
                                                          .to_json(orient='records'))
    rec.update(prediction_written=pred.get('written'), time=_now(), host=_host(), script=SCRIPT,
               commits=H.code_commits())
    _write_json(out_dir() / 'step0.json', rec)
    print(write_step0_table(), flush=True)
    return {k: v for k, v in rec.items() if k not in ('commits',)}


def smoke_point(pt: Point, scratch: Path) -> dict:
    """One read into a scratch folder (claims, row, provenance and run there; nothing in results or the archive), then
    the saved run loaded back and read at every stage (atlas's `read_run`): the read path and the stage readings on the
    point's wiring. Returns the row and the populations' sizes."""
    from kcq import stage_readings as R
    scratch = Path(scratch)
    row = read_point(pt, 'smoke', root=scratch, dirs=[scratch / 'read'], reclaim=True)
    run = R.load_step3_run(point_tag(pt), dirs=[scratch / 'read'])
    rd = R.read_run(run, seed_index=seed_index(pt.calyx, pt.wiring))
    return dict(row=row, populations={'/'.join(k): int(v.values.shape[0]) for k, v in rd.populations.items()},
                skipped={'/'.join(k): v for k, v in rd.skipped.items()})


# ---------------------------------------------------------------------------
# the score and the archive's table

def part_score() -> dict:
    """The score (design item 3, with item 4's readings beside): experiments/E144_score.py's `part_score`, which only
    this part imports, so that the score can be committed while read lanes run (E108's `_still_clean` guard)."""
    from experiments import E144_score as S
    return S.part_score()


def variant_kind(pt: Point) -> str:
    """A read's description for the archive's params.csv (ASCII)."""
    st = SETTINGS[pt.setting]
    return (f"C2's family on {CALYCES[pt.calyx]} (E130's fit of the calyx), setting {pt.setting} ({st.label}), wiring "
            f"{pt.wiring}" + (f" (E133's within-type claw permutation, seed {WIRING_SEEDS[pt.wiring]})"
                              if pt.wiring != 'real' else '')
            + f', read on the {pt.panel} input, both response rates held')


VARIANT_PARAMS: Tuple[str, ...] = PARAM_KEYS


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv in out_dir(), for the model-output archive's params.csv: one row per read the archive folder holds
    (`arch`, default archive_dir(); <tag>_run.json in read/): the tag, the folder, the variant, calyx, connectome,
    setting, panel and wiring (with its seed), a description (`variant_kind`), both conditions' absolute parameters
    ('off_<name>' with APL off), the start (the calyx's fit's parameters with the setting's factors, 'start_<name>'),
    the normalization's terms (pn_inhadd, pn_inhsc; empty for the family), the factors, the fit, and where and when it
    ran. The tag is unique. The family's runs record no normalization (their PN stage is built with olfsysm's defaults,
    no pn_inhadd / pn_inhsc kwargs), so those are filled with E135's OLFSYSM_DEFAULTS, as E142's table does."""
    from experiments.E135_antennal_lobe_dendrites import OLFSYSM_DEFAULTS
    arch = Path(arch) if arch is not None else archive_dir()
    rows = []
    for p in sorted((arch / 'read').glob('*_run.json')):
        d = json.loads(p.read_text())
        pt = parse_tag(d['tag'])
        pn = d.get('pn_kws') if isinstance(d.get('pn_kws'), dict) else {}
        if pt.setting == 'family' and not pn:
            pn = {k: OLFSYSM_DEFAULTS[k] for k in ('pn_inhadd', 'pn_inhsc')}
        r = dict(tag=d['tag'], variant=d.get('variant'), part='read', subdir='read', calyx=pt.calyx,
                 connectome=CALYCES[pt.calyx], setting=pt.setting, panel=pt.panel, wiring=pt.wiring,
                 wiring_seed=WIRING_SEEDS[pt.wiring], model=d.get('model'), reading=d.get('reading'),
                 kind=variant_kind(pt), fit=d.get('fit'), aplpn_factor=d.get('aplpn_factor'),
                 kckc_factor=d.get('kckc_factor'), kckc_abs=d.get('kckc_abs'),
                 pn_inhadd=pn.get('pn_inhadd'), pn_inhsc=pn.get('pn_inhsc'))
        r.update({k: (d.get('params') or {}).get(k) for k in VARIANT_PARAMS})
        r.update({f'off_{k}': (d.get('params_off') or {}).get(k) for k in VARIANT_PARAMS})
        r.update({f'start_{k}': (d.get('start') or {}).get(k) for k in VARIANT_PARAMS})
        r.update(host=d.get('host'), lane=d.get('lane'), time=d.get('time'))
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        raise FileNotFoundError(f'no <tag>_run.json under {arch / "read"}')
    if df.tag.duplicated().any():
        raise ValueError(f'runs named twice: {sorted(df.tag[df.tag.duplicated()])}')
    df.to_csv(out_dir() / 'variants.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# the lock

# the planner's expectation per calyx (the entry's table): megamat falls with A1's changes, A1 "holds", what would fail it
PREDICTION: Dict[str, dict] = {
    'H': dict(megamat_falls='9 in 10', a1_holds='3 in 4', would_fail_it='the claw default'),
    'L': dict(megamat_falls='9 in 10', a1_holds='1 in 3',
              would_fail_it='the APL16 somata with APL on; its low KC-KC scale (7.7 against 26.4)'),
    'R': dict(megamat_falls='9 in 10', a1_holds='1 in 2', would_fail_it='validation2'),
}
OUTCOMES: Dict[str, str] = {
    'a': "A1 replicates on both male calyces (about 1 in 5): the hold is not one fly's; A1, with E141's doses, goes to "
         "Tom and Yang as a candidate base read on three calyces.",
    'b': "Megamat falls on both and A1 holds on neither (about 1 in 4): the rise keeps the held-out panels only with the "
         "hemibrain's matrix; step 0's table says what that matrix has that the male ones lack, and a planner weighs "
         "moving it between calyces (E133's arms) under the stage target.",
    'c': "It replicates on one calyx (about 1 in 3): the two male calyces differ in APL's synapses at the boutons and in "
         "their fitted KC-KC scale (7.7 and 25.8); step 2 reads the doses on the other's permutations.",
    'd': "B holds on a male calyx without the rise (about 1 in 8): the need for the rise is the hemibrain's, and so is "
         "E141's question of what carries it.",
    'e': "A male calyx's family is outside the hemibrain's wirings at a megamat stage (about 1 in 2 at the boutons; "
         "beside any of the above): the stage table on the male CNS goes to Tom, whose shared models are built on it.",
    'f': "The claw default changes a class on the hemibrain (about 1 in 5): the male calyces are then read against H, "
         "not against the eight wirings, and the note says so.",
}
CLASSES: Dict[str, str] = {
    'e136': "E136's classes with E139's price, on L's and R's four wirings (each setting against the calyx's own family "
            "on the same wiring; megamat against the held-out half, all flies beside).",
    'a1_replicates': 'A1 replicates on a calyx: A1 holds there.',
    'rise_needed': 'The rise is needed there: B does not hold, and D1.8, A1 or D5.8 does.',
    'family_inside': "The calyx's family is inside the hemibrain's wirings: the family's strength error on the calyx's "
                     "real wiring is inside the range of the hemibrain family's eight wirings, at the boutons, the claws "
                     "and the somata.",
    'claw_default': "The claw default changes it: on H, a megamat stage or a held-out panel of a setting differs from the "
                    "'first' real wiring's by more than its spread; the classes are then read against H.",
    'dose_ladder': "A description: validation2's change along B, D1.8, A1 and D5.8 on the real wiring, beside the "
                   "hemibrain's (+0.080, +0.021, -0.016, -0.036).",
    'no_class': "No class for a wiring outside a realism band or whose holds fail. The number of calyces and settings "
                "tried is reported with any that holds. No calyx, dose or setting becomes a base without Tom's or "
                "Yang's word.",
}
# amendments to the lock after its first commit (each a dict: id, when, by whom, what, why, what it gates), in order;
# each written by `write_lock` before any number it defines is read
AMENDMENTS: Tuple[dict, ...] = (
    dict(id='e143_resampled', when='2026-11-26 16:41',
         by="viki's Opus code subagent (claude-opus-5-5, xhigh), at the brief of viki's queue session",
         what="Design item 4's \"E143's intervals over flies on every change, if they are in\" made concrete (they are "
              "in: E143 completed 2026-11-26 11:34). Every read is read once against all the flies of each stage and "
              "panel, as E143 reads saved runs (E143's read_runs: E134's default readings with the curve, noise and "
              "criterion of all the stage's flies, this experiment's seed index), into outputs/"
              "E144_a1_malecns_calyces/views/. E143's 1,000 resamples of the flies (its resamples: seed 0, N_BOOT "
              "1000, each megamat stage and each held-out panel drawn with replacement from one generator, so that "
              "megamat's four stages are resampled together; the flies of E143's lock) are the references, in all the "
              "flies' one-fly units, the unresampled flies last. For every change in classes.json's class tables (each "
              "setting against its calyx's family on the wirings it was read on): its 5th to 95th percentile over the "
              "resamples (E143's CI), the share of resamples better than the bar and worse than the spread; each class "
              "recomputed on every resample (E143's reclass: E136's ladder with the spreads, bars and the reasons that "
              "do not depend on the flies as scored) and the share of resamples that give the scored class; the "
              "derived classes (A1 replicates, the rise is needed, B holds, megamat falls with A1) and the outcomes "
              "(a) to (d) recomputed the same way. E143's joint splits are not read, so its \"a change stands\" is "
              "given for the held-out panels only (the interval excluding zero); 'stands' is None on megamat's rows. "
              "Tables: e143_intervals_changes.csv, e143_class_shares.csv, e143_intervals.json "
              "(experiments/E144_score.py's part_resampled).",
         why="The entry's design item 4 names E143's intervals beside the score; E143's note found that A1's \"holds\" "
             "on the hemibrain stands on 54 % of its resamples, so a class on a male calyx is read with its share.",
         gates='nothing: a description beside the classes, which stay as scored against the held-out half'),
)
JUDGMENTS: Tuple[str, ...] = (
    "The start of every hold is the calyx's fit's absolute parameters (E130's counted fit, its fits.csv row) with the "
    "setting's factors; the threshold (APL off) and the APL->KC scale with KC->APL and bouton->APL tied (APL on) are then "
    "held exactly (E77's hold_rates, E136's APL bracket 0.1 to 10). E136's reads start from the family's held megamat "
    "point, which the male calyces do not have.",
    "A permutation (E133's within-type claw permutation of seed 1, 2 or 3) runs at its calyx's fit's scales and "
    "operating points (the same KCs, so the same counts on the APL16 input).",
    "The breadth behind the APL->KC weights is the calyx wiring's APL-off breadth on the lab's APL16 input (E133's "
    "calyx_tables; recomputed on a permutation), on every panel and setting, as E136's and E139's reads hold it.",
    "RISE and APLPN are E141's rounded x3.24 and x0.155 (E139's K and B). Dabs is RISE times the family's KC-KC scale on "
    "the hemibrain's real wiring (E83's saved_fit('real'), 26.43): 85.6.",
    "F's and B's normalization is E141's pn_terms (pn_inhadd 1000, pn_inhsc by E77's inhsc_for), E139's rung 1000.",
    "A read's tag is '<calyx>_<setting>_<panel>_<wiring>'; a megamat read's row has half 'held' (the score reads it "
    "against the held-out half, all flies beside), the others none.",
    "The stage readings' seed index of a wiring is its position in (real, P1, P2, P3), so that H's real wiring is "
    "pooled and noised as the family's real wiring (E136's seed_index 0).",
    "The KC-KC matrix of a male calyx is E130's (Pratyush's KCclaw-to-KCclaw table, matrices.json's sha256); H's is "
    "olfsysm's hemibrain matrix among H's KCs ('argmax' claws), as E130's refit has it.",
)


def lock_dict() -> dict:
    """The lock (the entry's definitions, settings, wirings, operating points, classes and outcomes, the judgments
    above, the planner's provenance and the amendments, none at first), with each calyx's family as read
    (`calyx_fit`: fit, tag, absolute parameters, counts) and Dabs's scale."""
    from kcq import stage_readings as R
    from experiments import E77_input_screen as E77
    lk = R.require_lock()
    fits = {}
    for c in CALYCES:
        f = calyx_fit(c)
        fits[c] = dict(connectome=f['connectome'], fit=f['name'], tag=f['tag'], folder=f['folder'], params=f['params'],
                       target_on=f['target_on'], target_off=f['target_off'], size=f['size'], n_kcs=f['n_kcs'],
                       off_on_ratio=f['target_off'] / f['target_on'],
                       outputs_sha256={c2: E77.frame_sha256(f['saved'][c2][1]) for c2 in ('on', 'off')})
    from experiments import E134_population_baseline as B
    return dict(
        experiment=EXPERIMENT, script=SCRIPT, part='lock', planner=PLANNER,
        entry='experiment_priorities.md item E144 (proposed and moved by viki 2026-11-26 10:53)',
        definitions=dict(
            calyces={c: dict(connectome=con, wirings={w: WIRING_SEEDS[w] for w in WIRINGS_OF[c]},
                             wiring_rule=('the real wiring' if c == 'H' else
                                          "the real wiring and E133's within-type claw permutations of seeds 1 to 3"))
                     for c, con in CALYCES.items()},
            settings={k: dict(label=v.label, rung=v.rung, aplpn=v.aplpn, kckc=v.kckc, kckc_abs=v.kckc_abs,
                              calyces=list(v.calyces)) for k, v in SETTINGS.items()},
            rise=RISE, aplpn=APLPN, rung=RUNG, pn_terms=pn_terms(), dabs_scale=dabs_scale(),
            hemibrain_real_kckc=hemibrain_real_kckc(),
            operating_points=dict(megamat=float(B.held_rate('megamat', lk)),
                                  validation2=float(lk['populations']['kc_soma__validation2__wild_type']
                                                    ['response_rate']['mean_over_flies']),
                                  mixtures=float(B.held_rate('mixtures', lk)),
                                  apl16="the calyx's fit's two response counts",
                                  apl_off="APL on times the calyx's fit's APL-off / APL-on response counts"),
            panels=list(PANELS),
            bar="A calyx's bar is twice the standard error of a change over its four wirings.",
            spreads=("The hemibrain family's over its eight wirings (E136's), the yardstick of every class so far; the "
                     "calyx's own SD over four wirings beside."),
            reads=dict(step1=148, step2_per_calyx=24)),
        fits=fits, classes=CLASSES, outcomes=OUTCOMES, prediction=PREDICTION, judgments=list(JUDGMENTS),
        e134_lock_sha256=lk.get('sha256'), amendments=list(AMENDMENTS))


VOLATILE = ('written', 'host', 'commits', 'amended', 'amended_host', 'amended_commits')


def lock_path() -> Path:
    return out_dir() / LOCK_NAME


def committed_lock_path() -> Path:
    return H.HERE / 'experiments' / LOCK_NAME


def _same_lock(a: Mapping, b: Mapping) -> bool:
    strip = lambda d: json.loads(json.dumps({k: v for k, v in d.items() if k not in VOLATILE}, default=str))
    return strip(a) == strip(b)


def _amended_only(cur: Mapping, new: Mapping) -> bool:
    """Whether `new` is the lock `cur` with amendments appended to its own (AMENDMENTS) and nothing else changed."""
    strip = lambda d: json.loads(json.dumps({k: v for k, v in d.items() if k not in VOLATILE + ('amendments',)},
                                            default=str))
    a = json.loads(json.dumps(list(cur.get('amendments') or []), default=str))
    b = json.loads(json.dumps(list(new.get('amendments') or []), default=str))
    return strip(cur) == strip(new) and len(b) > len(a) and b[:len(a)] == a


def write_lock(force: bool = False) -> dict:
    """results/E144_a1_malecns_calyces/E144_lock.json and its committed copy experiments/E144_lock.json (`lock_dict`
    with the time, host and commits). A lock already written is kept when it is the same; a different one is refused
    once any read is claimed or done, unless `force` or the only change is amendments appended (AMENDMENTS): the lock is
    then written with its first writing's time, host and commits kept and the amendment's beside ('amended',
    'amended_host', 'amended_commits')."""
    cur = _read_json(lock_path())
    new = lock_dict()
    if cur and _same_lock(cur, new):
        return cur
    if cur and _amended_only(cur, new):
        new = dict(new, written=cur.get('written'), host=cur.get('host'), commits=cur.get('commits'), amended=_now(),
                   amended_host=_host(), amended_commits=H.code_commits())
        _write_json(lock_path(), new)
        committed_lock_path().write_text(json.dumps(new, indent=1, default=str) + '\n')
        return new
    if cur and not force and (len(reads_table()) or any((out_dir() / 'claims').glob('*'))):
        raise RuntimeError(f'{lock_path()} differs from the code\'s lock and reads exist: add an amendment instead')
    new = dict(new, written=_now(), host=_host(), commits=H.code_commits())
    _write_json(lock_path(), new)
    committed_lock_path().write_text(json.dumps(new, indent=1, default=str) + '\n')
    return new


def require_lock() -> dict:
    """The lock as written (results), checked against the code's (`lock_dict`, AMENDMENTS included) and against its
    committed copy; RuntimeError otherwise."""
    cur = _read_json(lock_path())
    if not cur:
        raise RuntimeError(f'no lock yet ({lock_path()}): run --part lock first')
    com = _read_json(committed_lock_path())
    if not com or not _same_lock(cur, com):
        raise RuntimeError(f'the lock differs from its committed copy {committed_lock_path()}')
    if not _same_lock(cur, lock_dict()):
        raise RuntimeError("the lock differs from the code's definitions (lock_dict)")
    return cur


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('lock', 'step0', 'step0_table', 'read', 'smoke', 'score',
                                                       'resampled', 'variants'))
    ap.add_argument('--step', nargs='+', default=['1'], choices=('1', '2'))
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--point', help="--part smoke: the point's tag (e.g. L_family_megamat_real)")
    ap.add_argument('--scratch', help='--part smoke: the scratch folder')
    a = ap.parse_args(argv)
    if a.part == 'lock':
        print(json.dumps(write_lock(), indent=1, default=str)[:2000])
    elif a.part == 'step0':
        print(json.dumps(part_step0(), indent=1, default=str))
    elif a.part == 'step0_table':
        print(write_step0_table())
    elif a.part == 'read':
        part_read(a.lane, a.step, a.reverse, a.reclaim)
    elif a.part == 'smoke':
        if not a.point or not a.scratch:
            raise SystemExit('--part smoke needs --point and --scratch')
        print(json.dumps(smoke_point(parse_tag(a.point), Path(a.scratch)), indent=1, default=str))
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'resampled':
        from experiments import E144_score as S
        print(json.dumps(S.part_resampled(), indent=1, default=str))
    elif a.part == 'variants':
        print(f'variants.csv: {len(write_variants())} runs', flush=True)


if __name__ == '__main__':
    main()
