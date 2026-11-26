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
    """A calyx's tables on one of its wirings (E133's `calyx_tables` at the wiring's seed)."""
    raise NotImplementedError


def build_simulator(calyx: str, wiring: str, params: Mapping[str, float], orn: Optional[pd.DataFrame] = None,
                    extra_kws: Optional[Mapping] = None) -> Tuple:
    """The calyx's family on a wiring as a fixed-parameter simulator (E133's `build_simulator` with a panel's input in
    place of the lab's and the setting's kwargs beside the calyx's); (simulate, info) with the closure's cells."""
    raise NotImplementedError


def build(calyx: str, wiring: str, panel: str, params: Mapping[str, float], extra_kws: Optional[Mapping] = None):
    """E134's `Built` of the calyx's wiring on a panel's input at `params` (every stage saved by `run_stages`)."""
    raise NotImplementedError


def hold_point(pt: Point, lock: Optional[Mapping] = None) -> dict:
    """Both response rates of a point held exactly at its operating point (E77's `hold_rates`)."""
    raise NotImplementedError


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None,
               dirs: Optional[Sequence[Path]] = None, root: Optional[Path] = None) -> Optional[dict]:
    """One read, claimed under <root>/claims/<tag>/ (E134's `claim`): held (`hold_point`), built again at the held
    parameters and run with every stage saved, saved in `dirs` (default `read_dirs`), its row in
    reads_<host>_<lane>.csv and its provenance in provenance_read_<host>_<lane>.json. None if done or claimed."""
    raise NotImplementedError


def reads_table() -> pd.DataFrame:
    """Every lane's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def part_read(lane: str, steps: Sequence[str] = ('1',), reverse: bool = False, reclaim: bool = False) -> None:
    """A lane of reads: every planned point of `steps`, reversed with `reverse`, skipping what is done or claimed."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# step 0

def step0_identity() -> pd.DataFrame:
    """Each calyx's builder at its counted fit's parameters against the fit's saved runs (0 differing counts and
    responses, APL on and off), at the present olfsysm."""
    raise NotImplementedError


def step0_glomeruli() -> pd.DataFrame:
    """Every glomerulus of the megamat, validation2 and mixture inputs (and the APL16 input) against each calyx's
    boutons: its boutons and claws per calyx; the check that none has no bouton."""
    raise NotImplementedError


def step0_matrices() -> pd.DataFrame:
    """The KC-KC matrices of the three calyces and the hemibrain's ('first' KCs) beside the hemibrain's shuffle within
    type pairs (E122's SH): the unevenness of the KCs' summed input, the weight between KCs that share an input
    glomerulus, reciprocity, and the weight by KC type pair."""
    raise NotImplementedError


def step0_scales() -> pd.DataFrame:
    """The absolute KC-KC and APL->bouton scales each setting reaches on each calyx, beside the hemibrain's eight
    wirings'."""
    raise NotImplementedError


def write_prediction() -> dict:
    """The planner's prediction per calyx (the entry's table), into prediction.json before any read."""
    raise NotImplementedError


def part_step0() -> dict:
    """Step 0: the tables above as step0_<name>.csv, step0.json, prediction.json."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# the score and the archive's table

def part_score() -> dict:
    """The score (E141's, per calyx; the module docstring). Not written yet."""
    raise NotImplementedError


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv for the model-output archive's params.csv: one row per read in the archive's read/."""
    raise NotImplementedError


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
# amendments to the lock after its first commit (each a dict: when, by whom, what, why); none yet
AMENDMENTS: Tuple[dict, ...] = ()
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


VOLATILE = ('written', 'host', 'commits')


def lock_path() -> Path:
    return out_dir() / LOCK_NAME


def committed_lock_path() -> Path:
    return H.HERE / 'experiments' / LOCK_NAME


def _same_lock(a: Mapping, b: Mapping) -> bool:
    strip = lambda d: json.loads(json.dumps({k: v for k, v in d.items() if k not in VOLATILE}, default=str))
    return strip(a) == strip(b)


def write_lock(force: bool = False) -> dict:
    """results/E144_a1_malecns_calyces/E144_lock.json and its committed copy experiments/E144_lock.json (`lock_dict`
    with the time, host and commits). A lock already written is kept when it is the same; a different one is refused
    once any read is claimed or done (a change after that goes in as an amendment, AMENDMENTS), unless `force`."""
    cur = _read_json(lock_path())
    new = lock_dict()
    if cur and _same_lock(cur, new):
        return cur
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
    ap.add_argument('--part', required=True, choices=('lock', 'step0', 'read', 'score', 'variants'))
    ap.add_argument('--step', nargs='+', default=['1'], choices=('1', '2'))
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    a = ap.parse_args(argv)
    if a.part == 'lock':
        print(json.dumps(write_lock(), indent=1, default=str)[:2000])
    elif a.part == 'step0':
        print(json.dumps(part_step0(), indent=1, default=str))
    elif a.part == 'read':
        part_read(a.lane, a.step, a.reverse, a.reclaim)
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'variants':
        print(f'variants.csv: {len(write_variants())} runs', flush=True)


if __name__ == '__main__':
    main()
