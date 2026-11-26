#!/usr/bin/env python3
"""E142's step 1 (design item 2): the recorded PN dendrites imposed on the model's PN stage under the calyx of A1, the
one stage-fit arm that holds, at fixed scales with both response rates re-held: the lock's 56 runs on the family's
eight wirings and, since the amendment of 2026-11-26 13:59 (below), 40 more, the main arm on A1's own PN stage.

The entry is item E142 of experiment_priorities.md (proposed and moved in by viki 2026-11-26 05:10 under the idle-host
rule; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 04:48). Step 0 (no simulation) is
viki's, in experiments/E142_dendrite_bouton_step.py, whose `--part step1` calls `part_step1` here; this module holds
step 1's runner (hal's claim, 2026-11-26 05:12). It is E135's step-2 runner (the antennal lobe held to the recorded PN
dendrites, experiments/E135_antennal_lobe_dendrites.py, `run_point2`) with the calyx's two scales changed before the
hold, and it calls E135's functions rather than copying them.

The lock's arms (LOCK_ARMS, the first seven of STEP1_ARMS; the entry's design item 2), each on the eight wirings:
- on the megamat input, the dendrites imposed (E135's 'imposed' table) with A1's calyx factors (APL->bouton and KC-KC);
  with each factor alone; with KC-KC at A1's and APL->bouton at none (its scale 0, in both conditions); with A1's
  factors and the table's glomerulus labels permuted, and with its odor labels permuted (E135's 'imposed_glom_perm'
  and 'imposed_odor_perm', each wiring's one seeded permutation, `E135.arm_seed`, so that each pairs with E135's own
  permuted run of the wiring);
- on the APL16 input, the dendrites imposed for its 12 megamat odors with A1's factors (E135's APL16 'imposed' table:
  paa, euc, -aPine and B-myr keep the model's PNs), for the bands and the APL effect.

Readings taken where the entry leaves them open for the lock's arms (the narrowest; each is also in the docstrings; the
amendment below replaces the first for the main arm):
- A1's calyx factors are E136's classed point of A1 under the adopted measure (`a1_factors`: fit_A1.json's
  'classed_z' in results/E136_stage_population_fit/measure_strength/, APL->bouton x0.155 and KC-KC x3.24), applied as
  E136 applies its factors: each wiring's own fit's wAPLPN_scale and wKCKC_scale times the factor, in both conditions
  (APL->bouton is 0 with APL off anyway). A1's antennal lobe (its weaker normalization) is not taken: the PN stage
  that is not imposed (the glomeruli outside the 28 and the odors the table does not set) and the distribution the
  dendrites are ranked onto are the family's, as in E135's imposed arm, so that the change from E135's runs is the
  calyx alone.
- "APL->bouton at none" is wAPLPN_scale 0 (APL's output to the boutons removed; the boutons' input to APL, tied to the
  held APL scale, stays).
- The hold is E135's: E77's `hold_rates` from the wiring's fit's parameters with the factors applied (the somatic
  threshold with APL off, then the tied APL->KC, KC->APL and bouton->APL scales with APL on), on megamat at E134's held
  targets (E136's reads' operating point), on APL16 at the wiring's own fit's counts (`E135.step2_targets`).

Each point (`run_point1`), as E135's `run_point2` with the factors: the family's PN stage on the point's input and
wiring (`E135.family_stage`; the lane's first point on a wiring first checks that the unchanged family reproduces its
saved run exactly, `lane_identity`, with E135's `family_identity`); the arm's table (`E135.arm_table`, clipped, ruling
7 of E135); route 1, olfsysm's `pn_rates` on the set entries (`E135.arm_input`, `E135.input_columns`); the calyx
factors (`arm_params`); both response rates re-held (`hold`); the point built again at the held parameters and run with
every stage (`E135._build`, `E135._run`), which must reproduce the hold exactly (`E135._hold_identity`); its PN rates
against its table, exact (`E135.run_pn_check`); the imposed run's check in E134's reading (`E135.imposed_run_check`);
the largest PN and bouton rates (`E135.stage_rates`); saved in stage_readings' layout, locally and in the archive
(`E135.save_point` into `step1_dirs`), with its row in runs_1_<host>_<lane>.csv and its provenance in
provenance_1_<host>_<lane>.json.

Step 1 amended (2026-11-26 13:59, after the first holds were read and before any step-1 run was scored; planned by a
Fable subagent, claude-fable-5-2, spawned from hal 2026-11-26 13:48; the entry's *Step 1 amended* and the lock's
amendment). On the lock's table (the dendrites ranked onto the family's own PN rates, pn_inhadd 31.4) KC-KC x3.24
cannot be held at megamat's operating point: megamat_imposed_a1's threshold ended at the hold's bound, x0.1 of the
wiring's fit, on all eight wirings, and every megamat arm with KC-KC x3.24 on that table has the same APL-off run. So:
- the lock's seven arms (LOCK_ARMS) are kept as the record: imposed_a1 "no class: a hold fails"; imposed_kckc and
  imposed_kckc_aplpn_none not run (placeholder claims, lane 'stopped'; their APL-off run is imposed_a1's), "no class: a
  hold fails"; the two permuted tables read once (M3), then stopped; imposed_aplpn and apl16 imposed_a1 run as they are;
- the amended main arm replaces the reading above in which the change from E135's runs was the calyx alone: the
  dendrites imposed on A1's own PN stage (AMENDED_ARMS with pn_stage 'A1'), A1's classed point whole (`a1_point`:
  pn_inhadd 1017.3 and pn_inhsc 598.85 from E136's `pn_norm_kws(fitted_z('A1'))`, APL->bouton x0.155, KC-KC x3.24).
  Its table is E135's imposed table built on A1's own PN stage (`a1_stage`: E135's `family_stage` with A1's
  normalization): the dendrites' arrangement ranked onto A1's PN rates on the 28 glomeruli, every other entry A1's, so
  the arm differs from A1 in the arrangement of those entries alone. Megamat imposed, its two permuted controls (E135's
  seeds) and the APL16 input, eight wirings each (32 runs), each held as A1's own reads were (`hold_a1`: E136's
  `hold_at` of A1 at its classed point from the wiring's held megamat point, `start_params`, to `read_targets`' counts,
  the same brackets, with route 1's pn_rates beside A1's kwargs). A lane's first A1-stage point on a wiring and panel
  checks that the unchanged A1 (no table) reproduces its saved run there (`a1_identity`: megamat on T and APL16 on every
  wiring its reads, megamat on S its fit's best evaluation when that is its classed point);
- "APL->bouton at none" replaced (imposed_aplpn_none): the family's PN stage with KC-KC at the family's and APL->bouton
  0, the lock's hold, eight runs;
- not run: each factor alone on A1's stage. Refused: a wider threshold bracket, and the APL-off condition at KC-KC x1.
Step 1 is then 66 runs (26 of the lock's 56, 40 new). The score's classes as amended (`classify_step1`): "the dendrites
carry under the calyx that holds" on the amended main arm against A1's own megamat reads on T, with the A1-stage
permuted controls and the A1-stage APL16 arm's bands; "APL at the boutons hides the step" on imposed_aplpn and
imposed_aplpn_none against E135's imposed arm, an arm whose hold fails left out and named; outcome (g) "not readable"
when the amended main arm's hold fails or it leaves a band (`step1_outcome`); beside, classing nothing, the amended
main arm's boutons against the no-step yardstick and E143's interval over resampled flies on its change from A1
(`e143_beside`); the lock's own reading and a table of the holds (`holds_table`) as the record.

The guards: E134's measure gate (E135's `family_stage` checks it); E136's fit of A1 (`a1_factors`); and the
authorization of step 1's runs (`require_step1_allowed`): E135's runs entered the dendrites by olfsysm's `pn_rates`
with Yang's OK of 2026-11-24 09:43, and viki asks Yang whether that OK covers these runs (the entry's Status), so no
run starts until results/E142_dendrite_bouton_step/step1_authorized.json records the answer (written by the session
that relays it; `write_step1_authorization`). A dry run (`--part step1_dry`) builds every table on one wiring without
any simulation of the mushroom body.

Outputs: <shared>/results/E142_dendrite_bouton_step/ (runs_1_<host>_<lane>.csv, provenance_1_<host>_<lane>.json,
identity_1_<host>_<lane>.json, claims/, step1_dry_<wiring>.json); outputs/E142_dendrite_bouton_step/step1/ (local,
large: the runs, also in the archive's data/model_outputs/E142_dendrite_bouton_step/step1/).

Usage (from the repository root, as modules; `--tags` runs only the points named, e.g. a first read):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E142_imposed_calyx_runs --part step1 --lane <host>[-<n>] [--reverse] [--tags <tag> ...]
    ./venv/bin/python -u -m experiments.E142_imposed_calyx_runs --part step1_dry [--wiring real] [--tags <tag> ...]
    ./venv/bin/python -u -m experiments.E142_imposed_calyx_runs --part routes_check [--wirings real P1 ...] --lane <l>
    ./venv/bin/python -u -m experiments.E142_imposed_calyx_runs --part routes_summary

The routes' check (`part_routes_check`, `part_routes_summary`; the entry's "Step 1's route" of 2026-11-26 10:53):
whether route 2, E135's solved ORN input, returns every table step 1 uses to within 0.5 Hz with no mushroom-body run
(check 1), with where each odor's fixed point lies (`route2_mean_map`); checks 2 and 3 (8 route-2 runs of E135's
imposed arm against its route-1 runs) only once check 1 passes. Outputs: routes_check_1_<host>_<lane>.csv,
routes_check_1_odors_<host>_<lane>.csv and routes_check.json in the results folder.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import socket
import time
from collections import namedtuple
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E142_dendrite_bouton_step'
SCRIPT = 'experiments/E142_imposed_calyx_runs.py'
PART = 'step1'
AUTH_FILE = 'step1_authorized.json'
# E136's measure under which A1's classed point is read (E136's adopted measure, strength by odor; E139's MEASURE)
E136_MEASURE = 'strength'
A1 = 'A1'

# an arm of step 1: its input panel, its name, the E135 table it is built as, the two calyx factors, each 'A1' (A1's
# factor), 'family' (x1, the wiring's own fit's scale) or 'none' (the scale 0), and its PN stage (PN_STAGES): 'family'
# (the lock's: the family's PN stage, the hold from the wiring's fit's parameters) or 'A1' (the amendment's: A1's own
# PN stage, A1's classed point whole, held as A1's reads were; both factors A1's)
Step1Arm = namedtuple('Step1Arm', 'panel name table aplpn kckc pn_stage', defaults=('family',))
PN_STAGES: Tuple[str, ...] = ('family', 'A1')
LOCK_ARMS: Tuple[Step1Arm, ...] = (
    Step1Arm('megamat', 'imposed_a1', 'imposed', 'A1', 'A1'),
    Step1Arm('megamat', 'imposed_aplpn', 'imposed', 'A1', 'family'),
    Step1Arm('megamat', 'imposed_kckc', 'imposed', 'family', 'A1'),
    Step1Arm('megamat', 'imposed_kckc_aplpn_none', 'imposed', 'none', 'A1'),
    Step1Arm('megamat', 'imposed_a1_glom_perm', 'imposed_glom_perm', 'A1', 'A1'),
    Step1Arm('megamat', 'imposed_a1_odor_perm', 'imposed_odor_perm', 'A1', 'A1'),
    Step1Arm('apl16', 'imposed_a1', 'imposed', 'A1', 'A1'),
)
# the amendment of 2026-11-26 13:59 (module docstring), after the lock's seven so that their tags and the lanes' order
# stay as they were
AMENDED_ARMS: Tuple[Step1Arm, ...] = (
    Step1Arm('megamat', 'a1stage_imposed', 'imposed', 'A1', 'A1', 'A1'),
    Step1Arm('megamat', 'a1stage_imposed_glom_perm', 'imposed_glom_perm', 'A1', 'A1', 'A1'),
    Step1Arm('megamat', 'a1stage_imposed_odor_perm', 'imposed_odor_perm', 'A1', 'A1', 'A1'),
    Step1Arm('apl16', 'a1stage_imposed', 'imposed', 'A1', 'A1', 'A1'),
    Step1Arm('megamat', 'imposed_aplpn_none', 'imposed', 'none', 'family'),
)
STEP1_ARMS: Tuple[Step1Arm, ...] = LOCK_ARMS + AMENDED_ARMS
N_LOCK_RUNS = 56
N_AMENDED_RUNS = 40
FACTOR_PARAMS: Dict[str, str] = {'aplpn': 'wAPLPN_scale', 'kckc': 'wKCKC_scale'}
FACTOR_CHOICES: Tuple[str, ...] = ('A1', 'family', 'none')
# how an A1-stage point is held (`hold_a1`)
A1_HOLD_RULE = ("E136's hold_at of A1 at its classed point (the calyx factors on the wiring's held megamat point, "
                "start_params; A1's normalization), to read_targets' counts (megamat: E134's held targets; APL16: the "
                "wiring's own fit's counts), E77's hold_rates with E136's brackets, route 1's pn_rates beside A1's "
                "kwargs, as A1's own reads were held")
# what each arm is, for the archive's params.csv
ARM_KINDS: Dict[str, str] = {
    'imposed_a1': "the recorded PN dendrites imposed (E135's table, route 1) with A1's APL->bouton and KC-KC factors",
    'imposed_aplpn': "the recorded PN dendrites imposed with A1's APL->bouton factor alone (KC-KC the family's)",
    'imposed_kckc': "the recorded PN dendrites imposed with A1's KC-KC factor alone (APL->bouton the family's)",
    'imposed_kckc_aplpn_none': "the recorded PN dendrites imposed with A1's KC-KC factor and APL->bouton at none (0)",
    'imposed_a1_glom_perm': "the recorded PN dendrites imposed, glomerulus labels permuted, with A1's two factors",
    'imposed_a1_odor_perm': "the recorded PN dendrites imposed, odor labels permuted, with A1's two factors",
    'a1stage_imposed': "the recorded PN dendrites imposed on A1's own PN stage (ranked onto A1's PN rates, route 1), "
                       "A1's classed point whole (its normalization and both calyx factors), held as A1's reads were",
    'a1stage_imposed_glom_perm': "the recorded PN dendrites imposed on A1's own PN stage, glomerulus labels permuted, "
                                 "A1's classed point whole",
    'a1stage_imposed_odor_perm': "the recorded PN dendrites imposed on A1's own PN stage, odor labels permuted, A1's "
                                 "classed point whole",
    'imposed_aplpn_none': "the recorded PN dendrites imposed (E135's table, route 1) with APL->bouton at none (0) and "
                          "KC-KC the family's",
}


def _e135():
    from experiments import E135_antennal_lobe_dendrites as E135
    return E135


def out_dir() -> Path:
    """<shared>/results/E142_dendrite_bouton_step/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def step1_dirs() -> List[Path]:
    """Where step 1's runs are written and read: outputs/E142_dendrite_bouton_step/step1/ (local), then the archive's
    data/model_outputs/E142_dendrite_bouton_step/step1/ (stage_readings' step-3 layout)."""
    from kcq import stage_readings as R
    return R._archive_dirs(EXPERIMENT, PART)


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def _read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text()) if Path(path).exists() else {}


def _write_json(path: Path, obj) -> None:
    path = Path(path)
    tmp = path.with_name(f'.{path.name}.{_host()}.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')
    os.replace(tmp, path)


def tag_of(arm: Step1Arm, wiring: str) -> str:
    """'<panel>_<arm>_<wiring>', e.g. 'megamat_imposed_a1_glom_perm_P1'."""
    return f'{arm.panel}_{arm.name}_{wiring}'


def planned_step1(reverse: bool = False, arms: Optional[Sequence[Step1Arm]] = None) -> List[Tuple[Step1Arm, str]]:
    """Step 1's points (arm, wiring), the arms of `arms` (default STEP1_ARMS: the lock's 56 points, then the
    amendment's 40) in their order, each on E135's eight wirings (E134's S then T); reversed with `reverse`. The lock's
    points that are not run are blocked by placeholder claims (lane 'stopped'), which every lane skips."""
    pts = [(a, w) for a in (STEP1_ARMS if arms is None else arms) for w in _e135().WIRINGS]
    return pts[::-1] if reverse else pts


def select_tags(pts: Sequence[Tuple[Step1Arm, str]], tags: Optional[Sequence[str]]) -> List[Tuple[Step1Arm, str]]:
    """The points of `pts` whose tags are in `tags` (all of them for None or empty), in `pts`' order. ValueError for a
    tag that names no planned point."""
    if not tags:
        return list(pts)
    want = set(tags)
    unknown = sorted(want - {tag_of(a, w) for a, w in pts})
    if unknown:
        raise ValueError(f'not planned step-1 points: {unknown}')
    return [(a, w) for a, w in pts if tag_of(a, w) in want]


def a1_factors() -> dict:
    """A1's two calyx factors: E136's classed point of A1 under E136_MEASURE (`fitted_z`, fit_A1.json's 'classed_z'),
    as factors exp(z) of 'aplpn' and 'kckc', with the file's path and sha256. RuntimeError if A1 is not fitted."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        z = F.fitted_z(A1)
        path = F.out_dir() / f'fit_{A1}.json'
    finally:
        F.set_measure(old)
    missing = [k for k in FACTOR_PARAMS if k not in z]
    if missing:
        raise RuntimeError(f"A1's classed point lacks {missing}: {z}")
    E135 = _e135()
    return dict(aplpn=float(np.exp(z['aplpn'])), kckc=float(np.exp(z['kckc'])), z={k: float(z[k]) for k in FACTOR_PARAMS},
                source=str(path), sha256=E135._sha256(path) if path.exists() else None, measure=E136_MEASURE)


def a1_point() -> dict:
    """A1's classed point whole, for the amendment's A1-stage arms: E136's `fitted_z('A1')` under E136_MEASURE (every
    coordinate: aplpn, kckc and the normalization's inhadd, inhsc), the normalization's kwargs (E136's `pn_norm_kws`:
    pn_inhadd, pn_inhsc), the calyx factors exp(z) of 'aplpn' and 'kckc', the fit record's path and sha256, and whether
    its classed point is its best evaluation ('classed_is_best'; then its S wirings' megamat runs are at the point).
    RuntimeError unless E135's verdict leaves A1 in its 'search' mode (A1's PN stage its normalization) and the point
    has all four coordinates."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        mode = F.require_e135()['mode']
        if mode != 'search':
            raise RuntimeError(f"E135's verdict puts A1 in mode {mode!r}: its PN stage is not its normalization")
        z = F.fitted_z(A1)
        path = F.out_dir() / f'fit_{A1}.json'
        rec = _read_json(path)
    finally:
        F.set_measure(old)
    missing = [k for k in (*FACTOR_PARAMS, *F.PN_COORDS) if k not in z]
    if missing:
        raise RuntimeError(f"A1's classed point lacks {missing}: {z}")
    kws = F.pn_norm_kws(z)
    if set(kws) != {'pn_inhadd', 'pn_inhsc'}:
        raise RuntimeError(f"A1's normalization kwargs are {kws}")
    best, classed = rec.get('best_z'), rec.get('classed_z')
    return dict(z={k: float(v) for k, v in z.items()}, pn_kws={k: float(v) for k, v in kws.items()},
                aplpn=float(np.exp(z['aplpn'])), kckc=float(np.exp(z['kckc'])), source=str(path),
                sha256=_e135()._sha256(path) if path.exists() else None, measure=E136_MEASURE, e135_mode=mode,
                best_eval=rec.get('best_eval'), classed_is_best=bool(best is not None and (classed is None or
                                                                                          classed == best)))


def arm_factors(arm: Step1Arm, a1: Mapping[str, float]) -> Dict[str, float]:
    """The arm's factor on each of the two calyx scales: A1's ('A1'), 1 ('family') or 0 ('none'). An A1-stage arm
    (pn_stage 'A1') has A1's two factors and no other."""
    if arm.pn_stage not in PN_STAGES:
        raise ValueError(f'{arm}: pn_stage={arm.pn_stage!r}')
    if arm.pn_stage == 'A1' and (arm.aplpn, arm.kckc) != ('A1', 'A1'):
        raise ValueError(f"{arm}: an arm on A1's PN stage is A1's classed point whole (both factors 'A1')")
    out = {}
    for k in FACTOR_PARAMS:
        c = getattr(arm, k)
        if c not in FACTOR_CHOICES:
            raise ValueError(f'{arm}: {k}={c!r}')
        out[k] = float(a1[k]) if c == 'A1' else (1.0 if c == 'family' else 0.0)
    return out


def arm_params(p: Mapping[str, float], factors: Mapping[str, float]) -> Dict[str, float]:
    """The wiring's fit's parameters `p` with each calyx scale (FACTOR_PARAMS) times its factor; the others as they
    are."""
    q = dict(p)
    for k, name in FACTOR_PARAMS.items():
        if name not in q:
            raise KeyError(f"the fit's parameters have no {name}")
        q[name] = float(q[name]) * float(factors[k])
    return q


def require_step1_allowed(root: Optional[Path] = None) -> dict:
    """The authorization of step 1's runs (results/.../step1_authorized.json: who allowed them, when, in whose words,
    recorded by which host): E135 entered the dendrites by olfsysm's `pn_rates` with Yang's OK of 2026-11-24 09:43, and
    the entry's Status has viki ask Yang whether it covers these runs. PermissionError until it is written
    (`write_step1_authorization`)."""
    p = (out_dir() if root is None else Path(root)) / AUTH_FILE
    d = _read_json(p)
    if not d.get('allowed'):
        raise PermissionError(f"{p}: step 1's runs wait for Yang's (or Tom's) word that the OK for pn_rates in E135 "
                              f"covers them (the entry's Status); `write_step1_authorization` records it")
    return d


def write_step1_authorization(who: str, when: str, words: str, relayed_by: str, root: Optional[Path] = None) -> dict:
    """Records the word that allows step 1's runs (`require_step1_allowed`): `who` (Yang or Tom), `when`, their `words`
    as quoted, and the host whose session relayed them."""
    rec = dict(allowed=True, who=who, when=when, words=words, relayed_by=relayed_by, recorded=_now(), host=_host())
    _write_json((out_dir() if root is None else Path(root)) / AUTH_FILE, rec)
    return rec


def claim(tag: str, lane: str, reclaim: bool = False) -> bool:
    """Claims a point by mkdir of results/E142_dendrite_bouton_step/claims/<tag>/ (E134's `claim`, rooted here)."""
    from experiments import E134_population_baseline as E134R
    return E134R.claim(tag, lane, reclaim=reclaim, root=out_dir())


def release(tag: str) -> None:
    from experiments import E134_population_baseline as E134R
    E134R.release(tag, root=out_dir())


def identity_path(lane: str) -> Path:
    """results/.../identity_1_<host>_<lane>.json: the lane's identity of the unchanged family per wiring ('family:
    <wiring>') and each point's hold identity (by tag)."""
    return out_dir() / f'identity_1_{_host()}_{lane}.json'


def _identity_put(lane: str, key: str, rec: Mapping) -> None:
    path = identity_path(lane)
    ids = _read_json(path)
    ids[key] = dict(rec)
    _write_json(path, ids)


def lane_identity(panel: str, wiring: str, lane: str, finfo: Mapping) -> dict:
    """The identity check of a lane's first point on a wiring (E135's `family_identity` on `family_stage`'s simulator:
    the unchanged family reproduces its saved run exactly), cached per host, lane and wiring in `identity_path`.
    RuntimeError if it does not hold."""
    key = f'family:{wiring}'
    old = _read_json(identity_path(lane)).get(key)
    if old and old.get('exact'):
        return dict(old, cached=True)
    rec = _e135().family_identity(panel, wiring, finfo['sim'])
    rec.update(host=_host(), lane=lane, time=_now(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _identity_put(lane, key, rec)
    if not rec.get('exact'):
        raise RuntimeError(f'{wiring} on {panel}: the unchanged family does not reproduce its saved run on lane {lane}: '
                           f'{rec}')
    return rec


def a1_stage(panel: str, wiring: str, point: Mapping) -> Tuple[object, dict]:
    """A1's own PN stage on a panel's input and a wiring: E135's `family_stage` built with A1's normalization
    (`point['pn_kws']`, `a1_point`), read into E135's `PNStage` and checked there against olfsysm's own layers at that
    normalization. RuntimeError unless the build's layers carry exactly A1's pn_inhadd and pn_inhsc. Returns (the stage,
    `family_stage`'s record)."""
    stage, finfo = _e135().family_stage(panel, wiring, extra_kws=dict(point['pn_kws']))
    got = dict(pn_inhadd=float(stage.sp.pn_inhadd), pn_inhsc=float(stage.sp.pn_inhsc))
    if any(got[k] != float(v) for k, v in point['pn_kws'].items()):
        raise RuntimeError(f"{wiring} {panel}: the A1 build's PN layer has {got}, not A1's {point['pn_kws']}")
    return stage, finfo


def a1_saved(panel: str, wiring: str, point: Mapping):
    """A1's saved run at its classed point on a panel and wiring, for `a1_identity`: its read (E136's `load_read`, under
    E136_MEASURE: megamat on T, APL16 on every wiring); on megamat's S wirings its fit's best evaluation when that is
    its classed point (`point['classed_is_best']`; E136's `megamat_runs`); None otherwise or if absent."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        if panel == F.PANEL and wiring in F.S_WIRINGS:
            if not point.get('classed_is_best'):
                return None
            runs = F.megamat_runs(A1, [wiring], 'held')
            return None if runs is None else runs[wiring]
        return F.load_read(F.ReadPoint(A1, panel, wiring, 'held' if panel == F.PANEL else None))
    finally:
        F.set_measure(old)


def _n_differ(a: pd.DataFrame, b: pd.DataFrame) -> int:
    """How many entries of two outputs differ, odor columns abbreviated and the second put in the first's column order
    (-1 if their shapes or odors differ)."""
    from kcq import stage_readings as R
    a, b = R._abbrev(a), R._abbrev(b)
    if a.shape != b.shape or set(a.columns) != set(b.columns):
        return -1
    return int((a.to_numpy(dtype=float) != b[list(a.columns)].to_numpy(dtype=float)).sum())


def a1_identity(panel: str, wiring: str, sim, point: Mapping) -> dict:
    """The unchanged A1 (no imposed table) on the A1 build's simulator `sim` (`a1_stage`'s) reproduces A1's saved run
    of the panel and wiring (`a1_saved`) exactly: both conditions at that run's held parameters give its spike counts and
    responses, and the run is at A1's classed point (its coordinates and normalization are `point`'s). The hold itself
    is E136's `hold_at`, whose start and targets are the read's. 'exact' is None where no saved run is at the point."""
    run = a1_saved(panel, wiring, point)
    rec = dict(panel=panel, wiring=wiring, pn_stage='A1', checks=[])
    if run is None:
        return dict(rec, exact=None, note="no saved run of A1 at its classed point on this panel and wiring")
    info = run.info or {}
    z = {k: float(v) for k, v in (info.get('z') or {}).items()}
    kws = {k: float(v) for k, v in dict(info.get('pn_kws') or {k: info[k] for k in ('pn_inhadd', 'pn_inhsc')
                                                                if k in info}).items()}

    def same(a: Mapping[str, float], b: Mapping[str, float]) -> bool:
        return set(a) == set(b) and all(np.isclose(a[k], b[k], rtol=1e-12, atol=0.0) for k in a)
    rec.update(tag=run.tag, z=z, pn_kws=kws, z_equal=same(z, point['z']), pn_kws_equal=same(kws, point['pn_kws']))
    for c, key in (('on', 'params'), ('off', 'params_off')):
        r, s = sim(dict(info[key]))
        rec['checks'].append(dict(name=f'{c} counts', differ=_n_differ(s, run.conditions[c]['spike_counts'])))
        rec['checks'].append(dict(name=f'{c} responses', differ=_n_differ(r, run.conditions[c]['responses'])))
    rec['exact'] = bool(rec['z_equal'] and rec['pn_kws_equal'] and all(x['differ'] == 0 for x in rec['checks']))
    return rec


def a1_lane_identity(panel: str, wiring: str, lane: str, finfo: Mapping, point: Mapping) -> dict:
    """The identity check of a lane's first A1-stage point on a wiring and panel (`a1_identity` on the A1 build's
    simulator), cached per host, lane, panel and wiring in `identity_path` ('A1:<panel>:<wiring>'). RuntimeError if a
    saved run exists and is not reproduced."""
    key = f'A1:{panel}:{wiring}'
    old = _read_json(identity_path(lane)).get(key)
    if old and old.get('exact') is not False and 'exact' in old:
        return dict(old, cached=True)
    rec = a1_identity(panel, wiring, finfo['sim'], point)
    rec.update(host=_host(), lane=lane, time=_now(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
    _identity_put(lane, key, rec)
    if rec.get('exact') is False:
        raise RuntimeError(f'{wiring} on {panel}: the unchanged A1 does not reproduce its saved run on lane {lane}: {rec}')
    return rec


class Step1:
    """What step 1's points share in one process: E135's lock (its 28 shared glomeruli), E134's lock, A1's factors, A1's
    classed point (`a1_point`, read when an A1-stage arm first needs it), the PN stage per (panel, wiring, PN stage)
    (the family's, or A1's own: `a1_stage`; built once each, the lane's identity on the wiring checked on that build,
    the build then freed unless `keep_vars`, for the dry run's PN-layer check) and the dendrite means in E135's fit
    units."""

    def __init__(self, lane: str, lock134: Optional[Mapping] = None, a1: Optional[Mapping] = None,
                 identity: bool = True, a1pt: Optional[Mapping] = None, keep_vars: bool = False):
        E135 = _e135()
        self.lane = lane
        self.lock = E135.read_lock()
        self.glomeruli = [str(g) for g in self.lock['glomeruli']]
        self.lock134 = lock134
        self.a1 = dict(a1) if a1 is not None else a1_factors()
        self.a1pt = None if a1pt is None else dict(a1pt)
        self.identity = identity
        self.keep_vars = keep_vars
        self._stages: Dict[Tuple[str, str, str], Tuple[object, dict]] = {}
        self._vars: Dict[Tuple[str, str, str], dict] = {}
        self._den: Optional[pd.DataFrame] = None

    def point(self) -> dict:
        """A1's classed point whole (`a1_point`), read once."""
        if self.a1pt is None:
            self.a1pt = a1_point()
        return self.a1pt

    def stage(self, panel: str, wiring: str, pn_stage: str = 'family'):
        key = (panel, wiring, pn_stage)
        if key not in self._stages:
            if pn_stage == 'family':
                stage, finfo = _e135().family_stage(panel, wiring)
                ident = lane_identity(panel, wiring, self.lane, finfo) if self.identity else dict(skipped=True)
            elif pn_stage == 'A1':
                pt = self.point()
                stage, finfo = a1_stage(panel, wiring, pt)
                ident = a1_lane_identity(panel, wiring, self.lane, finfo, pt) if self.identity else dict(skipped=True)
            else:
                raise ValueError(f'{pn_stage=}')
            port = finfo.get('identity')
            if self.keep_vars:
                self._vars[key] = finfo
            finfo = None
            gc.collect()
            self._stages[key] = (stage, dict(port=port, lane_identity=ident, pn_stage=pn_stage,
                                             pn_inhadd=float(stage.sp.pn_inhadd), pn_inhsc=float(stage.sp.pn_inhsc)))
        return self._stages[key]

    def free_vars(self) -> None:
        """Drops the builds kept for the dry run."""
        self._vars.clear()
        gc.collect()

    def dendrites(self) -> pd.DataFrame:
        if self._den is None:
            E135 = _e135()
            self._den = E135.dendrite_means(E135.FIT_UNITS, self.glomeruli)
        return self._den


def point_table(arm: Step1Arm, wiring: str, ctx: Step1) -> Tuple[pd.DataFrame, List[str], dict, object, dict]:
    """A point's table of PN rates as E135 builds it for the arm's E135 table on the wiring (`E135.arm_table` with
    E135's seed for the wiring, clipped), on the arm's PN stage (the family's, or A1's own: the dendrites ranked onto
    that stage's own rates, its other entries its own), the glomeruli it sets, its record, the PN stage and the stage's
    record."""
    E135 = _e135()
    stage, srec = ctx.stage(arm.panel, wiring, arm.pn_stage)
    table, set_gl, trec = E135.arm_table(E135.Arm(arm.panel, arm.table), stage, ctx.dendrites(), ctx.glomeruli,
                                         seed=E135.arm_seed(wiring))
    return table, set_gl, trec, stage, srec


def hold(wiring: str, orn: pd.DataFrame, fam: Mapping, params: Mapping[str, float], extra: Mapping,
         targets: Mapping, label: str) -> dict:
    """Both response rates re-held (E77's `hold_rates`) on the family's simulator of `wiring` on `orn` with `extra`
    (route 1's pn_rates), the breadth at the lab input's, from `params` (the fit's with the arm's factors): the
    APL->bouton and KC-KC scales stay `params`'."""
    from experiments import E83_cross_wiring_fit as E83
    from experiments.E77_input_screen import hold_rates
    p = dict(params)
    return hold_rates(lambda: E83.build_simulator(wiring, orn=orn, breadth=fam['breadth'], params=dict(p),
                                                  extra_kws=dict(extra)),
                      target_off=targets['target_off'], target_on=targets['target_on'], params=p, label=label)


def hold_a1(wiring: str, panel: str, extra: Mapping, point: Mapping, label: str) -> dict:
    """Both response rates of an A1-stage point re-held as A1's own reads were (A1_HOLD_RULE): E136's `hold_at` of A1
    at its classed point (`point['z']`) on the wiring and the panel's input, with `extra` (route 1's pn_rates) beside A1's
    kwargs. Returns `hold_at`'s dict (hold, kws, orn, targets, p0, ...)."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        return F.hold_at(A1, dict(point['z']), wiring, panel, label=label, extra_kws=dict(extra))
    finally:
        F.set_measure(old)


def pn_only_route1(finfo: Mapping, stage, table: pd.DataFrame, set_gl: Sequence[str], set_od: Sequence[str],
                   tol: Optional[float] = None) -> dict:
    """Route 1 on the PN layer alone, with no simulation of the mushroom body: a PN-stage build's olfsysm variables
    (`family_stage`'s mp, rv and wiring) with the table's set entries as rv.pn.preset_rates (glomeruli and odors in
    olfsysm's order, as fit_mb_model sets them), olfsysm's ORN, LN and PN layers run again (`run_ORN_LN_sims`,
    `run_PN_sims`), and pyolfsysm's `pn_evoked_rates` against the table on the set entries and against the stage's own
    rates (the port, 'preset') on the others, both to `tol` (E135's ROUTE1_TOL_HZ); the preset rates are cleared after.
    On A1's build this says whether `pn_rates` composes with A1's pn_inhadd and pn_inhsc."""
    import olfsysm as osm
    from pyolfsysm.mb_model import pn_evoked_rates
    E135 = _e135()
    tol = E135.ROUTE1_TOL_HZ if tol is None else float(tol)
    mp, rv, W = finfo['mp'], finfo['rv'], finfo['info']['W']
    want = E135.route1_pn_rates(table, set_gl, set_od).reindex(index=stage.glomeruli, columns=stage.odors)
    rv.pn.preset_rates = want.to_numpy(dtype=float)
    try:
        osm.run_ORN_LN_sims(mp, rv)
        osm.run_PN_sims(mp, rv)
        got = pn_evoked_rates(mp, rv, wPNKC=W, odors=list(stage.odors))
    finally:
        rv.pn.preset_rates = np.zeros((0, 0))
    got.index = [str(g) for g in got.index]
    if list(got.index) != list(stage.glomeruli):
        raise RuntimeError("olfsysm's PN rows are not the stage's glomeruli")
    m = np.isfinite(want.to_numpy(dtype=float))
    G, T = got.to_numpy(dtype=float), want.to_numpy(dtype=float)
    own = stage.run(reading='preset').reindex(index=got.index, columns=got.columns).to_numpy(dtype=float)
    d_set = np.abs(G[m] - T[m])
    d_un = np.abs(G[~m] - own[~m])
    rec = dict(n_set=int(m.sum()), max_abs_diff_set=float(d_set.max()) if d_set.size else 0.0,
               n_unset=int((~m).sum()), max_abs_diff_unset=float(d_un.max()) if d_un.size else 0.0, tol=tol,
               pn_inhadd=float(mp.pn.inhadd), pn_inhsc=float(mp.pn.inhsc))
    rec['ok'] = bool(rec['max_abs_diff_set'] <= tol and rec['max_abs_diff_unset'] <= tol)
    return rec


def runs1_table(root: Optional[Path] = None) -> pd.DataFrame:
    """Every row of every host's and lane's runs_1_*.csv, the last row per tag."""
    d = out_dir() if root is None else Path(root)
    fs = [pd.read_csv(p) for p in sorted(d.glob('runs_1_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def run_point1(arm: Step1Arm, wiring: str, ctx: Step1, table_path: Path, prov_path: Path, reclaim: bool = False,
               done: Optional[set] = None) -> Optional[dict]:
    """One point of step 1, claimed (None if done or claimed elsewhere; the claim is released on failure): E135's
    `run_point2` for the arm's E135 table, route 1, on the arm's PN stage. On the family's PN stage (the lock's arms)
    the calyx factors are applied to the wiring's fit's parameters before the lock's hold (`hold`, E135's targets); on
    A1's own PN stage (the amendment's) the point is held as A1's reads were (`hold_a1`: A1's classed point on the
    wiring's held megamat point, `read_targets`), and built and run with A1's normalization beside the table (module
    docstring)."""
    from experiments import E134_population_baseline as E134R
    from experiments.E77_input_screen import frame_sha256
    E135 = _e135()
    tag = tag_of(arm, wiring)
    if done is not None and tag in done:
        return None
    if not claim(tag, ctx.lane, reclaim):
        print(f'E142 {tag}: claimed by another process, skipped', flush=True)
        return None
    t0 = time.time()
    try:
        factors = arm_factors(arm, ctx.a1)
        table, set_gl, trec, stage, srec = point_table(arm, wiring, ctx)
        set_od = list(trec['set_odors'])
        kw, rrec = E135.arm_input(stage, table, set_gl, 1, set_odors=set_od)
        orn = E135.panel_orn(arm.panel)
        extra = dict(pn_rates=E135.input_columns(kw['extra_kws']['pn_rates'], orn))
        fam = E135.family_parts(wiring)
        point = None
        if arm.pn_stage == 'A1':
            point = ctx.point()
            factors = dict(aplpn=float(point['aplpn']), kckc=float(point['kckc']))
            ha = hold_a1(wiring, arm.panel, extra, point, f'E142 {tag}')
            h, p0, kws = ha['hold'], dict(ha['p0']), dict(ha['kws'])
            off = {k: kws.get(k) for k in point['pn_kws'] if float(kws.get(k, np.nan)) != float(point['pn_kws'][k])}
            if off or 'pn_rates' not in kws:
                raise RuntimeError(f"{tag}: the hold's kwargs are not A1's normalization with the table: {off}")
            t_off, t_on, f_off, f_on = ha['targets']
            tg = dict(target_off=int(t_off), target_on=int(t_on), frac_off_target=float(f_off),
                      frac_on_target=float(f_on), size=int(fam['n_kc']) * int(orn.shape[1]), hold_rule=A1_HOLD_RULE)
            start = "the wiring's held megamat point (E136's start_params)"
        else:
            p0 = arm_params(fam['p'], factors)
            tg = E135.step2_targets(arm.panel, wiring, orn.shape[1], fam, ctx.lock134)
            h = hold(wiring, orn, fam, p0, extra, tg, f'E142 {tag}')
            kws = dict(extra)
            start = "the wiring's fit (E83's wiring_parts)"
        params, params_off = dict(h['params']), dict(h['params_off'])
        for k, name in FACTOR_PARAMS.items():
            if float(params[name]) != float(p0[name]):
                raise RuntimeError(f'{tag}: the hold moved {name} ({p0[name]} -> {params[name]})')
        built = E135._build(arm.panel, wiring, params, None, kws)
        out = E135._run(arm.panel, wiring, built, params, params_off)
        ident = E135._hold_identity(arm.panel, wiring, out, h)
        _identity_put(ctx.lane, tag, ident)
        if not ident['exact']:
            raise RuntimeError(f"{tag}: the run at the held parameters is not the hold's: {ident}")
        pn = E135.run_pn_check(built, stage, table, set_gl, set_od, tol=E135.ROUTE1_TOL_HZ, strict=True)
        imp = E135.imposed_run_check(table, E135.e134_pn_rates(out['on']), set_gl, set_od)
        sr = E135.stage_rates(out)
        table_sha = frame_sha256(table)
        info = dict(model=E135.FAMILY, wiring=wiring, panel=arm.panel, reading='held', tag=tag, experiment=EXPERIMENT,
                    arm=arm.name, table_arm=arm.table, route=1, units=E135.FIT_UNITS, seed=trec['seed'],
                    factors=factors, factor_choice=dict(aplpn=arm.aplpn, kckc=arm.kckc), a1=ctx.a1,
                    pn_stage=arm.pn_stage, pn_inhadd=srec.get('pn_inhadd'), pn_inhsc=srec.get('pn_inhsc'),
                    a1_point=point, start=start, start_params=p0,
                    set_glomeruli=list(set_gl), set_odors=set_od, table_sha256=table_sha, stages=list(built.stages),
                    params=params, params_off=params_off, targets=tg, host=_host(), lane=ctx.lane, time=_now())
        E135.save_point(tag, out, built, info, table, dirs=step1_dirs())
        clip = trec.get('clip', {})
        row = dict(tag=tag, panel=arm.panel, arm=arm.name, table_arm=arm.table, wiring=wiring, route=1,
                   units=E135.FIT_UNITS, seed=trec['seed'], aplpn=arm.aplpn, kckc=arm.kckc,
                   aplpn_factor=factors['aplpn'], kckc_factor=factors['kckc'], pn_stage=arm.pn_stage,
                   pn_inhadd=srec.get('pn_inhadd'), pn_inhsc=srec.get('pn_inhsc'), start=start,
                   start_fixed_thr=float(p0['fixed_thr']),
                   **{k: float(v) for k, v in params.items() if np.isscalar(v)},
                   **{f'off_{k}': float(v) for k, v in params_off.items() if np.isscalar(v)},
                   n_kcs=int(out['on']['spike_counts'].shape[0]), n_odors=int(out['on']['spike_counts'].shape[1]),
                   **E134R._counts_row(out),
                   thr_factor=h['thr_factor'], apl_factor=h['apl_factor'], target_on=tg['target_on'],
                   target_off=tg['target_off'], frac_on_target=tg['frac_on_target'],
                   frac_off_target=tg['frac_off_target'], at_bound_thr=h['at_bound_thr'],
                   at_bound_apl=h['at_bound_apl'], hold_n_sims=h['n_sims'],
                   fit_fixed_thr=float(fam['p']['fixed_thr']),
                   n_set_glomeruli=len(set_gl), n_set_odors=len(set_od), n_mapped=trec.get('n_mapped'),
                   clip_n_set=clip.get('n_set'), clip_n_clipped=clip.get('n_clipped'),
                   clip_share=clip.get('share_clipped'), clip_largest_excursion_hz=clip.get('largest_excursion_hz'),
                   clip_named=clip.get('named'),
                   pn_check_max_set_hz=pn['max_abs_diff_set'], pn_check_max_unset_hz=pn['max_abs_diff_unset'],
                   imposed_check_max_hz=imp['max_abs_diff_hz'], imposed_check_named=imp['named'],
                   **sr, table_sha256=table_sha, input_sha256=frame_sha256(built.orn), stages='+'.join(built.stages),
                   seconds=round(time.time() - t0, 1), host=_host(), lane=ctx.lane, time=_now(),
                   threads=os.environ.get('OMP_NUM_THREADS'), identity_exact=True,
                   lane_identity_exact=bool(srec['lane_identity'].get('exact')),
                   **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
        E134R._append_csv(table_path, row)
        prov = _read_json(prov_path) or dict(part=PART, lane=ctx.lane, script=SCRIPT, points={})
        prov['points'][tag] = dict(
            info=info, built_info={k: v for k, v in built.info.items() if not str(k).startswith('_')}, arm=trec,
            route=rrec, pn_check=pn, imposed_check=imp, identity=ident, lane_identity=srec['lane_identity'],
            port_identity=srec['port'], targets=tg, hold=dict(thr_factor=h['thr_factor'], apl_factor=h['apl_factor'],
                                                              n_sims=h['n_sims'], start=start,
                                                              at_bound_thr=h['at_bound_thr'],
                                                              at_bound_apl=h['at_bound_apl']),
            table_sha256=table_sha, dendrites_sha256=frame_sha256(ctx.dendrites()),
            e135_lock_sha256=E135._sha256(E135.LOCK_COPY) if E135.LOCK_COPY.exists() else None,
            e134_lock_sha256=E135._sha256(E135.LOCK_E134()) if E135.LOCK_E134().exists() else None,
            authorization=_read_json(out_dir() / AUTH_FILE), model=H.model_provenance(orn_deltas=built.orn))
        _write_json(prov_path, prov)
        print(f"E142 {tag}: {row['seconds']:.0f} s, rate on {row['rate_on']:.4f} / off {row['rate_off']:.4f}, "
              f"APL factor {h['apl_factor']:.3f}, clipped {clip.get('n_clipped')} of {clip.get('n_set')}", flush=True)
        del built, out, h
        return row
    except BaseException:
        release(tag)
        raise
    finally:
        gc.collect()


def part_step1(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False,
               tags: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Step 1, resumable, one point claimed at a time (`planned_step1`, `run_point1`); every lane on any host runs the
    whole list (or only the points of `tags`, `select_tags`, e.g. a first read) and skips what another has claimed or
    done, the placeholder claims of the points not run included. Refuses before the authorization
    (`require_step1_allowed`), E136's fit of A1 (`a1_factors`) and E134's gate (`E135.family_stage`); an A1-stage point
    also needs E135's verdict to leave A1 in its 'search' mode (`a1_point`). Rows in runs_1_<host>_<lane>.csv,
    provenance in provenance_1_<host>_<lane>.json, identities in identity_1_<host>_<lane>.json, the runs in
    `step1_dirs`."""
    from experiments import E134_population_baseline as E134R
    E135 = _e135()
    require_step1_allowed()
    E135.require_e134_gate()
    E134R.check_wirings()
    if tuple(E135.WIRINGS) != tuple(E134R.WIRINGS):
        raise RuntimeError(f"E135's wirings {E135.WIRINGS} are not E134's {E134R.WIRINGS}")
    lock134 = E135.require_e134_lock()
    lane = lane or _host()
    pts = select_tags(planned_step1(reverse), tags)
    ctx = Step1(lane, lock134)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    msg = (f"E142 step 1, lane {lane}: {len(pts)} points; A1's factors APL->bouton x{ctx.a1['aplpn']:.4f}, KC-KC "
           f"x{ctx.a1['kckc']:.4f}")
    if any(a.pn_stage == 'A1' for a, _ in pts):
        msg += f"; A1's PN stage {ctx.point()['pn_kws']}"
    print(msg, flush=True)
    for arm, w in pts:
        done = set(runs1_table()['tag'].astype(str))
        run_point1(arm, w, ctx, table, prov, reclaim, done)
    return runs1_table()


def part_step1_dry(wiring: str = 'real', tags: Optional[Sequence[str]] = None) -> dict:
    """Every arm's table on one wiring (or the arms of `tags` there), without any simulation of the mushroom body: per
    PN stage (the family's; A1's own, `a1_stage`) the stage is built once (one fit_mb_model call, as every build), its
    olfsysm layers checked against E135's port at the stage's normalization, and no lane identity is run. Per arm: the
    set glomeruli and odors, the clip's counts, the table's range on the set glomeruli beside the stage's own rates
    there (the values it is ranked onto), the arm's absolute calyx scales (from the wiring's fit, or for an A1-stage arm
    from the wiring's held megamat point with A1's factors), and route 1 on the PN layer alone (`pn_only_route1`: whether
    olfsysm's `pn_rates` sets exactly the table's rates on the set entries at the stage's normalization, and leaves the
    others at the stage's own); step1_dry_<wiring>.json. Needs E134's gate and A1's fit, not the authorization."""
    from experiments import E136_stage_population_fit as F
    E135 = _e135()
    ctx = Step1('dry', identity=False, keep_vars=True)
    fam = E135.family_parts(wiring)
    arms = [a for a, _ in select_tags([(a, wiring) for a in STEP1_ARMS], tags)]
    out = dict(wiring=wiring, a1=ctx.a1, time=_now(), host=_host(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}, arms={}, stages={})
    if any(a.pn_stage == 'A1' for a in arms):
        out['a1_point'] = ctx.point()
    for key in dict.fromkeys((a.panel, a.pn_stage) for a in arms):
        for arm in [a for a in arms if (a.panel, a.pn_stage) == key]:
            table, set_gl, trec, stage, srec = point_table(arm, wiring, ctx)
            f = arm_factors(arm, ctx.a1)
            if arm.pn_stage == 'A1':
                p0, _ = F.start_params(wiring)
                p = dict(p0, **{name: float(p0[name]) * float(np.exp(ctx.point()['z'][k]))
                                for k, name in FACTOR_PARAMS.items()})
            else:
                p = arm_params(fam['p'], f)
            clip = trec.get('clip', {})
            od = list(trec['set_odors'])
            own = stage.run(reading='preset').loc[set_gl, od].to_numpy(dtype=float)
            tab = table.loc[set_gl, od].to_numpy(dtype=float)
            chk = pn_only_route1(ctx._vars[(arm.panel, wiring, arm.pn_stage)], stage, table, set_gl, od)
            out['stages'][f'{arm.panel}/{arm.pn_stage}'] = dict(pn_inhadd=srec['pn_inhadd'], pn_inhsc=srec['pn_inhsc'],
                                                                port=srec['port'])
            out['arms'][tag_of(arm, wiring)] = dict(
                pn_stage=arm.pn_stage, pn_inhadd=srec['pn_inhadd'], pn_inhsc=srec['pn_inhsc'],
                n_set_glomeruli=len(set_gl), set_odors=od, n_mapped=trec.get('n_mapped'),
                clip_n_set=clip.get('n_set'), clip_n_clipped=clip.get('n_clipped'), clip_named=clip.get('named'),
                clip_largest_excursion_hz=clip.get('largest_excursion_hz'),
                clipped=[(x['glomerulus'], x['odor']) for x in clip.get('outside', [])],
                factors=f, **{name: p[name] for name in FACTOR_PARAMS.values()},
                table_min_hz=float(np.nanmin(tab)), table_max_hz=float(np.nanmax(tab)),
                own_min_hz=float(np.nanmin(own)), own_max_hz=float(np.nanmax(own)),
                n_entries_changed=int(np.sum(np.abs(tab - own) > 1e-9)),
                pn_only_route1=chk)
            print(f"E142 dry {tag_of(arm, wiring)}: {len(set_gl)} glomeruli x {len(od)} odors set, clipped "
                  f"{clip.get('n_clipped')} of {clip.get('n_set')}, table {np.nanmin(tab):.2f} to {np.nanmax(tab):.2f} Hz "
                  f"(own {np.nanmin(own):.2f} to {np.nanmax(own):.2f}); route 1 on the PN layer: set "
                  f"{chk['max_abs_diff_set']:.2e} Hz, unset {chk['max_abs_diff_unset']:.2e} Hz, ok {chk['ok']}",
                  flush=True)
        ctx.free_vars()
    _write_json(out_dir() / f'step1_dry_{wiring}.json', out)
    print(json.dumps({k: dict(ok=v['pn_only_route1']['ok'], clipped=v['clip_n_clipped'], n_set=v['clip_n_set'])
                      for k, v in out['arms'].items()}, indent=1, default=str), flush=True)
    return out


# ---------------------------------------------------------------------------
# the routes' check (the entry's "Step 1's route", added 2026-11-26 10:53; planned by a Fable subagent,
# claude-fable-5-2, spawned from viki 2026-11-26 10:31): route 2, E135's ORN input solved so that the family's PN stage
# returns the table, may carry step 1's runs in place of olfsysm's `pn_rates` only if (1) with no mushroom-body
# simulation the solve returns every table step 1 uses to within 0.5 Hz on every solved entry and `solve_identity`
# passes, (2) E135's imposed arm at the family's calyx, run by route 2 on the eight wirings, has its PN rates within
# 0.5 Hz of the table, and (3) its stage errors agree with E135's saved route-1 runs by E135's ruling 9. Checks 2 and 3
# run only once check 1 passes (`part_routes_summary` says which), so this part runs check 1.

ROUTES_PART = 'routes_check'
ROUTES_FILE = 'routes_check.json'
# the input's own mean delta of each odor times these factors: where route 2's fixed point is looked for, from no
# input at all to nine times the input's (beyond that the solved entries reach the bracket's edge), closer near the
# input's own mean, where a target near the stage's own rates puts it
ROUTE2_MEAN_GRID: Tuple[float, ...] = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 1.5,
                                       2.0, 3.0, 5.0, 9.0)
# |F| (in the delta's units, averaged over the stage's glomeruli) at or under which a grid point counts as a fixed
# point: a target near the stage's own rates can put two roots closer together than the grid's points
ROUTE2_ROOT_TOL = 0.05
# E135's `solve_input` bracket of an entry's delta
ROUTE2_BRACKET: Tuple[float, float] = (-400.0, 2000.0)


def route2_tables() -> List[Step1Arm]:
    """The tables the lock's step 1 uses (LOCK_ARMS; the routes' check was read on these, before the amendment), one arm
    per distinct (panel, E135 table): megamat's imposed table, its two permuted tables (each wiring's own seeded
    permutation) and APL16's imposed table."""
    seen, out = set(), []
    for arm in LOCK_ARMS:
        if (arm.panel, arm.table) not in seen:
            seen.add((arm.panel, arm.table))
            out.append(arm)
    return out


def route2_target(stage, table: pd.DataFrame, set_gl: Sequence[str],
                  set_od: Sequence[str]) -> Tuple[List[str], List[str], pd.DataFrame]:
    """Route 2's target as E135's `arm_input` builds it: on the glomeruli with input, the table's rates on the set
    glomeruli and odors and the stage's own rates elsewhere ('preset' reading). Returns (the glomeruli with input, the
    set odors, the target)."""
    own = stage.run(reading='preset')
    has_input = [g for g in stage.glomeruli if np.any(stage.delta.loc[g].to_numpy() != 0)]
    setg = [g for g in set_gl if g in set(has_input)]
    od = [o for o in set_od if o in set(own.columns)]
    target = own.loc[has_input].copy()
    target.loc[setg, od] = table.loc[setg, od].to_numpy()
    return has_input, od, target


def route2_mean_map(stage, target: pd.DataFrame, glomeruli: Sequence[str], odors: Sequence[str],
                    grid: Sequence[float] = ROUTE2_MEAN_GRID, bracket: Tuple[float, float] = ROUTE2_BRACKET,
                    n_iter: int = 60, root_tol: float = ROUTE2_ROOT_TOL) -> pd.DataFrame:
    """Where route 2's fixed point lies, odor by odor. E135's `solve_input` looks, for each odor, for the mean delta m
    over the stage's glomeruli at which the deltas solved at m's LN (each entry at its target, E135's `_solve_entries`)
    have the mean m: a root of F(m) = mean(solved deltas) - m. F is read at m = factor x the input's own mean delta for
    each factor of `grid`, at the points where no entry reaches the bracket's edge.

    An odor whose F stays above 0 at every point has no input, under the stage's LN, that returns its target: its PNs
    ask for more input than the LN that input drives lets through (the stage's divisive normalization).

    Returns one row per odor: the stage's own and the target's PN rates summed over `glomeruli`, F at the input's own
    mean, the least and the largest F, whether F changes sign between neighbouring points or comes within `root_tol` of
    0 at one (a root inside the grid), and the number of points read."""
    E135 = _e135()
    sp = stage.sp
    G = stage.glomeruli
    nG = len(G)
    rows = np.array([G.index(g) for g in glomeruli])
    cols = np.array([stage.odors.index(o) for o in odors])
    T = target.loc[list(glomeruli), list(odors)].to_numpy(dtype=float)
    d = stage.delta.to_numpy(dtype=float)
    others = [i for i in range(nG) if i not in set(rows.tolist())]
    fixed = d[others][:, cols].sum(axis=0) if others else np.zeros(len(cols))
    k = E135.orn_layer(sp, np.zeros(1), np.ones((1, 1)))[0, 0]
    inhA, inhB = stage.ln(sp)
    inhA, inhB = inhA.copy(), inhB.copy()
    m0 = d[:, cols].mean(axis=0)
    Fs, ok = [], []
    for f in grid:
        m = m0 * float(f)
        a, b = E135.ln_layer_mean(sp, stage.spont.mean() + m[:, None] * k[None, :], nG)
        inhA[cols], inhB[cols] = a, b
        x, n_edge = E135._solve_entries(sp, stage.spont, rows, cols, T, k, inhA, inhB, 'preset', bracket, n_iter,
                                        E135.SOLVE_TOL_HZ / 100)
        Fs.append((x.sum(axis=0) + fixed) / nG - m)
        ok.append(n_edge == 0)
    Fs, ok = np.asarray(Fs), np.asarray(ok)
    own = stage.run(reading='preset')
    at_own = [i for i, f in enumerate(grid) if float(f) == 1.0]
    out = []
    for j, o in enumerate(odors):
        fj = Fs[ok, j]
        out.append(dict(odor=o, own_total_hz=float(own.loc[list(glomeruli), o].sum()),
                        target_total_hz=float(target.loc[list(glomeruli), o].sum()),
                        F_at_own_mean=float(Fs[at_own[0], j]) if at_own and ok[at_own[0]] else np.nan,
                        min_F=float(fj.min()) if len(fj) else np.nan, max_F=float(fj.max()) if len(fj) else np.nan,
                        root_in_grid=bool(len(fj) and ((np.sign(fj[:-1]) != np.sign(fj[1:])).any()
                                                       or (np.abs(fj) <= root_tol).any())),
                        n_grid=int(ok.sum())))
    return pd.DataFrame(out)


def route2_odor_diffs(stage, solved: pd.DataFrame, target: pd.DataFrame, glomeruli: Sequence[str],
                      odors: Sequence[str]) -> pd.Series:
    """The largest difference of each odor's PN rates from the target (Hz, 'preset' reading, over `glomeruli`), with
    the stage re-run exactly on route 2's solved input (its LN on that input)."""
    E135 = _e135()
    got = E135.PNStage(solved, pd.Series(stage.spont, index=stage.glomeruli), stage.sp).run(reading='preset')
    gl, od = list(glomeruli), list(odors)
    return (got.loc[gl, od] - target.loc[gl, od]).abs().max(axis=0)


def part_routes_check(wirings: Optional[Sequence[str]] = None, lane: Optional[str] = None) -> pd.DataFrame:
    """The routes' check 1, with no mushroom-body simulation: on each wiring, for each table step 1 uses
    (`route2_tables`), E135's route 2 (`arm_input` with route 2: `solve_input` on the glomeruli with input), its record,
    each odor's largest difference from the target and where its fixed point lies (`route2_mean_map`); and on each
    (panel, wiring) the solve's identity on the stage's own rates (`E135.solve_identity`). The family's PN stage is
    built as step 1 builds it (`Step1.stage`, without the lane identity, which runs the mushroom body). Rows in
    routes_check_1_<host>_<lane>.csv and routes_check_1_odors_<host>_<lane>.csv; `part_routes_summary` reads them."""
    from experiments import E136_stage_population_fit as F
    E135 = _e135()
    E135.require_e134_gate()
    lane = lane or _host()
    ctx = Step1(lane, identity=False)
    wirings = list(F.WIRINGS if not wirings else wirings)
    tpath = out_dir() / f'routes_check_1_{_host()}_{lane}.csv'
    opath = out_dir() / f'routes_check_1_odors_{_host()}_{lane}.csv'
    prov = out_dir() / f'provenance_routes_check_1_{_host()}_{lane}.json'
    _write_json(prov, dict(script=SCRIPT, part=ROUTES_PART, wirings=wirings, lane=lane, host=_host(), time=_now(),
                           grid=list(ROUTE2_MEAN_GRID), bracket=list(ROUTE2_BRACKET), tol_hz=E135.SOLVE_TOL_HZ,
                           **H.code_commits()))
    trows, orows, idents = [], [], {}
    for w in wirings:
        for arm in route2_tables():
            t0 = time.time()
            table, set_gl, trec, stage, _ = point_table(arm, w, ctx)
            has_input, od, target = route2_target(stage, table, set_gl, trec['set_odors'])
            key = (arm.panel, w)
            if key not in idents:
                ident = E135.solve_identity(stage, glomeruli=has_input)
                idents[key] = dict(ok=ident['ok'], max_rate_diff_hz=ident['max_rate_diff'],
                                   max_delta_diff=ident['max_delta_diff'])
            kw, rec = E135.arm_input(stage, table, set_gl, 2, set_odors=trec['set_odors'])
            diffs = route2_odor_diffs(stage, kw['orn'], target, has_input, od)
            mp_ = route2_mean_map(stage, target, has_input, od)
            mp_['solve_max_abs_diff_hz'] = mp_['odor'].map(diffs).astype(float)
            mp_.insert(0, 'table', arm.table)
            mp_.insert(0, 'panel', arm.panel)
            mp_.insert(0, 'wiring', w)
            orows.append(mp_)
            ln = rec.get('ln_change') or {}
            trows.append(dict(
                wiring=w, panel=arm.panel, table=arm.table, n_set_glomeruli=len(set_gl), n_set_odors=len(od),
                n_solved=rec.get('n_solved'), converged=bool(rec['converged']),
                max_abs_diff_hz=float(rec['max_abs_diff']), at_edge=int(rec['at_edge']), method=rec['method'],
                fixed_point=bool(rec['fixed_point']), bracketed=bool(rec['bracketed']), outer=int(rec['outer']),
                no_input_max_drift_hz=float(rec['no_input_max_drift_hz']),
                ln_inhA_max_abs=ln.get('inhA_max_abs'), ln_pn_inh_factor_max_rel=ln.get('pn_inh_factor_max_rel'),
                n_odors_without_root=int((~mp_['root_in_grid']).sum()),
                odors_without_root=';'.join(mp_.loc[~mp_['root_in_grid'], 'odor'].astype(str)),
                identity_ok=idents[key]['ok'], identity_max_rate_diff_hz=idents[key]['max_rate_diff_hz'],
                seconds=round(time.time() - t0, 1), host=_host(), lane=lane, time=_now(),
                queue_commit=H.code_commits().get('queue_commit')))
            pd.DataFrame(trows).to_csv(tpath, index=False)
            pd.concat(orows, ignore_index=True).to_csv(opath, index=False)
            r = trows[-1]
            print(f"E142 routes' check 1 {w} {arm.panel} {arm.table}: converged {r['converged']}, largest difference "
                  f"{r['max_abs_diff_hz']:.2f} Hz ({r['method']}), {r['n_odors_without_root']} of {len(od)} odors "
                  f"without a fixed point, identity {r['identity_ok']} ({r['seconds']:.0f} s)", flush=True)
    return pd.DataFrame(trows)


def routes_verdict(tables: pd.DataFrame, n_expected: int, tol: float) -> dict:
    """Check 1's verdict from its rows: it passes when every table step 1 uses was read on every wiring
    (`n_expected` rows) and each converged to within `tol` Hz, at no bracket edge, with the solve's identity passing.
    Checks 2 and 3 run only once it passes."""
    n = len(tables)
    bad = tables[~(tables['converged'].astype(bool) & (tables['max_abs_diff_hz'] <= tol) & (tables['at_edge'] == 0)
                   & tables['identity_ok'].astype(bool))] if n else tables
    passes = bool(n == n_expected and n and bad.empty)
    return dict(check1=dict(passes=passes, n_tables=n, n_expected=n_expected, tol_hz=tol, n_failing=int(len(bad)),
                            failing=[f"{r.panel} {r.table} {r.wiring}" for r in bad.itertuples()]),
                check2='to run' if passes else 'not run: check 1 failed',
                check3='to run' if passes else 'not run: check 1 failed',
                route='2 after checks 2 and 3' if passes else "1 (olfsysm's pn_rates), on Yang's word")


def part_routes_summary() -> dict:
    """Check 1's rows from every lane, its verdict (`routes_verdict`) and, per table, the solve's largest difference and
    the odors without a fixed point, into routes_check.json; every lane's rows together, in a fixed order, into
    routes_check_tables.csv and routes_check_odors.csv (names that do not depend on the lanes, for the Summary's keys)."""
    from experiments import E136_stage_population_fit as F
    E135 = _e135()
    d = out_dir()
    tabs = [pd.read_csv(p) for p in sorted(d.glob('routes_check_1_*.csv')) if '_odors_' not in p.name]
    ods = [pd.read_csv(p) for p in sorted(d.glob('routes_check_1_odors_*.csv'))]
    tables = pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame()
    odors = pd.concat(ods, ignore_index=True) if ods else pd.DataFrame()
    order = {w: i for i, w in enumerate(F.WIRINGS)}
    tab_order = {(a.panel, a.table): i for i, a in enumerate(route2_tables())}
    for df, name in ((tables, 'routes_check_tables.csv'), (odors, 'routes_check_odors.csv')):
        if len(df):
            key = df.apply(lambda r: (tab_order.get((r['panel'], r['table']), 99), order.get(r['wiring'], 99)), axis=1)
            df.iloc[np.argsort(key.to_numpy(), kind='stable')].to_csv(d / name, index=False)
    v = routes_verdict(tables, len(F.WIRINGS) * len(route2_tables()), E135.SOLVE_TOL_HZ)
    if len(tables):
        v['tables'] = tables[['wiring', 'panel', 'table', 'converged', 'max_abs_diff_hz', 'method',
                              'n_odors_without_root', 'odors_without_root', 'identity_ok']].to_dict('records')
    if len(odors):
        o = odors.assign(no_root=~odors['root_in_grid'].astype(bool),
                         target_over_own=odors['target_total_hz'] / odors['own_total_hz'])
        g = o.groupby(['panel', 'table', 'odor'], sort=False)
        v['odors'] = g.agg(n_wirings=('wiring', 'size'), n_without_root=('no_root', 'sum'),
                           target_over_own=('target_over_own', 'mean'), least_min_F=('min_F', 'min'),
                           solve_max_abs_diff_hz=('solve_max_abs_diff_hz', 'max')).reset_index().to_dict('records')
    v.update(time=_now(), host=_host(), script=SCRIPT, **{k: x for k, x in H.code_commits().items()
                                                         if not isinstance(x, (dict, list))})
    _write_json(d / ROUTES_FILE, v)
    print(json.dumps({k: v[k] for k in ('check1', 'check2', 'check3', 'route')}, indent=1, default=str), flush=True)
    return v


def load_point(tag: str):
    """A saved step-1 run (stage_readings' `load_step3_run` from `step1_dirs`); FileNotFoundError if absent."""
    from kcq import stage_readings as R
    return R.load_step3_run(tag, dirs=step1_dirs())


def _pn_stage_of(row: Mapping) -> str:
    """A runs row's PN stage ('family' for the rows written before the amendment, which have none)."""
    v = row.get('pn_stage')
    return 'family' if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)


def write_variants(root: Optional[Path] = None) -> pd.DataFrame:
    """step1_variants.csv: one row per saved step-1 run (`runs1_table`), for the model-output archive's params.csv:
    the tag, panel, arm, the E135 table, wiring, seed, the factors, the PN stage and its normalization, what the arm is
    (ARM_KINDS) and the absolute parameters of both conditions (E135's VARIANT_COLUMNS)."""
    E135 = _e135()
    d = out_dir() if root is None else Path(root)
    rows = []
    for r in runs1_table(d).to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=r['tag'], part=PART, subdir=PART, panel=r.get('panel'), arm=r.get('arm'),
                         table_arm=r.get('table_arm'), model=E135.FAMILY, wiring=r.get('wiring'), route=1,
                         seed=r.get('seed'), aplpn_factor=r.get('aplpn_factor'), kckc_factor=r.get('kckc_factor'),
                         pn_stage=_pn_stage_of(r), pn_inhadd=r.get('pn_inhadd'), pn_inhsc=r.get('pn_inhsc'),
                         kind=f"C2's family on wiring {r.get('wiring')}, on the {r.get('panel')} input, both response "
                              f"rates re-held: {ARM_KINDS.get(str(r.get('arm')), r.get('arm'))}",
                         **{c: r.get(c) for c in E135.VARIANT_COLUMNS if c in r}))
    df = pd.DataFrame(rows)
    df.to_csv(d / 'step1_variants.csv', index=False)
    return df


# ---------------------------------------------------------------------------
# the score (design item 3) and step 1's classes (the lock's 'apl_hides_step', 'dendrites_carry', 'no_class')

# the variants the score registers with E136's score functions: step 1's arms as 'E142_<arm>' (the APL16 arm shares
# its name with the megamat arm of the same factors), E135's saved imposed arms at the family's calyx as 'E135_<arm>'
E142_PREFIX = 'E142_'
E135_VARIANTS: Dict[str, str] = {'E135_imposed': 'imposed', 'E135_imposed_glom_perm': 'imposed_glom_perm',
                                 'E135_imposed_odor_perm': 'imposed_odor_perm'}
DOWNSTREAM: Tuple[str, ...] = ('pn_bouton', 'kc_claw', 'kc_soma')
# step 1's classes as amended (2026-11-26 13:59; module docstring): "APL at the boutons hides the step" on the family's
# PN stage, in imposed_aplpn or imposed_aplpn_none against E135's imposed arm; "the dendrites carry under the calyx that
# holds" on the main arm on A1's own PN stage against A1's own megamat reads on T, with the A1-stage permuted controls
APL_ARMS: Tuple[str, ...] = ('imposed_aplpn', 'imposed_aplpn_none')
CARRY_ARM = 'a1stage_imposed'
CONTROL_ARMS: Tuple[str, ...] = ('a1stage_imposed_glom_perm', 'a1stage_imposed_odor_perm')
CARRY_REFERENCE = 'A1'
# the lock's step-1 classes as written (experiments/E142_dendrite_bouton_step.py CLASSES), kept as the record
LOCK_APL_ARMS: Tuple[str, ...] = ('imposed_a1', 'imposed_aplpn', 'imposed_kckc_aplpn_none')
LOCK_CARRY_ARM = 'imposed_a1'
LOCK_CONTROL_ARMS: Tuple[str, ...] = ('imposed_a1_glom_perm', 'imposed_a1_odor_perm')
APL_REFERENCE = 'E135_imposed'
MIN_WIRINGS = 6
SE_MULTIPLE = 2.0
# the lock's arms whose APL-off run is another arm's exactly (APL->bouton is 0 with APL off; the same KC-KC and table),
# so that its hold's threshold is that arm's
SAME_OFF_RUN: Dict[str, str] = {'imposed_kckc': 'imposed_a1', 'imposed_kckc_aplpn_none': 'imposed_a1'}
# E143's name of A1's saved megamat views (its setting), for `e143_beside`
E143_A1 = 'E136:A1'


def variant_of(arm_name: str) -> str:
    """A step-1 arm's variant name for E136's score functions, 'E142_<arm>'."""
    return E142_PREFIX + arm_name


def step1_variants() -> List[str]:
    """Every step-1 arm's variant name, once (STEP1_ARMS' order)."""
    return list(dict.fromkeys(variant_of(a.name) for a in STEP1_ARMS))


def load_variant_run(pt, claw_after: bool = False):
    """The run of a registered variant's point (E136's `ReadPoint` with the variant as its arm): a step-1 arm's from
    `step1_dirs` ('<panel>_<arm>_<wiring>'), E135's imposed arms from E135's step-2 folders; None if absent (also for
    the claws read after a threshold, which none of these has)."""
    from kcq import stage_readings as R
    if claw_after:
        return None
    if pt.arm in E135_VARIANTS:
        tag, dirs = f'{pt.panel}_{E135_VARIANTS[pt.arm]}_{pt.wiring}', _e135().step2_dirs()
    elif str(pt.arm).startswith(E142_PREFIX):
        tag, dirs = f'{pt.panel}_{pt.arm[len(E142_PREFIX):]}_{pt.wiring}', step1_dirs()
    else:
        raise ValueError(f'{pt.arm!r} is not a variant of E142 step 1')
    try:
        return R.load_step3_run(tag, dirs=dirs)
    except FileNotFoundError:
        return None


def holds_failed(rows: pd.DataFrame, tol: float) -> List[str]:
    """The runs of a runs table (E142's runs_1 or E135's runs_2) whose hold failed: at a bracket's bound, or a response
    rate more than `tol` (relative) from its target; as tags."""
    bad = []
    for r in rows.to_dict('records'):
        miss = [abs(float(r[f'rate_{c}']) - float(r[f'frac_{c}_target'])) / float(r[f'frac_{c}_target'])
                for c in ('on', 'off')]
        if str(r.get('at_bound_thr')).lower() in ('true', '1', '1.0') or \
                str(r.get('at_bound_apl')).lower() in ('true', '1', '1.0') or not max(miss) <= tol:
            bad.append(str(r['tag']))
    return bad


def no_step_yardstick(F, stage: str = 'pn_bouton', half: Optional[str] = 'held') -> float:
    """The no-step yardstick of a stage against a half of its flies: all seven dendrite flies' mean strength views as
    the prediction (hal's `no_step_score`), scored against the half's mean views in all the stage's flies' one-fly
    units (E136's `Flies`, strength views), as step 0 reads it (yardsticks.csv, glomeruli 'dendrites')."""
    from data_prep import stage_population_data as D
    from experiments import E134_population_baseline as B
    from kcq import population_fit as P
    stages = ('pn_dendrite', stage)
    fl = F.Flies(half, stages, views=('strength',))
    pops = B.fly_side_read(fl.full['pn_dendrite'])[0]
    ref = fl.reference(stage, list(D.panel_odors(F.PANEL)))
    return float(P.no_step_score(pops, ref).views['strength'])


def mean_over(vals: Mapping[str, Mapping[str, float]], key: str, wirings: Sequence[str]) -> float:
    return float(np.mean([vals[w][key] for w in wirings])) if all(w in vals for w in wirings) else np.nan


def classify_step1(mega: Mapping[str, Mapping[str, Dict[str, float]]], spreads: Mapping[str, float], no_step: float,
                   bad: Mapping[str, Sequence[str]], flagged: Mapping[str, Sequence[str]], F,
                   apl_arms: Sequence[str] = APL_ARMS, carry_arm: str = CARRY_ARM,
                   control_arms: Sequence[str] = CONTROL_ARMS, carry_reference: str = CARRY_REFERENCE,
                   leave_out_failed: bool = True) -> dict:
    """Step 1's classes from the megamat values against the held-out half (`mega`: per variant, per wiring, E136's
    `megamat_values`), the family's spreads over its eight wirings, the no-step yardstick of the boutons against the
    held-out half, the runs whose holds failed (`bad`: per variant, tags) and the wirings whose APL16 run leaves a realism
    band where the family's does not (`flagged`: per variant). The defaults are the amendment's reading (2026-11-26
    13:59); the lock's own, kept as the record, is LOCK_APL_ARMS, LOCK_CARRY_ARM and LOCK_CONTROL_ARMS with
    `leave_out_failed` False.

    - 'apl_hides_step': in at least one of `apl_arms`, the boutons' strength error is below the no-step yardstick on
      MIN_WIRINGS of the eight wirings, and in E135's imposed arm (APL->bouton at the family's) it is not. An arm whose
      hold fails is left out and named ('left_out', with `leave_out_failed`; otherwise its failed hold is a reason); an
      arm with a flagged wiring or runs missing gives no class, and so does no arm left.
    - 'dendrites_carry': `carry_arm`'s error over the boutons, claws and somata ('downstream') is below
      `carry_reference`'s on T by more than SE_MULTIPLE standard errors (E136's `compare`: the larger of half the spread
      and the SD of the per-wiring changes over sqrt(4)), none of the three stages is worse than the reference's by more
      than its spread, and each of the three is below both `control_arms` (the mean over T); no class with a failed hold
      or a flagged wiring among the arm's or the controls' runs on T.
    The PN dendrites are given beside (they are imposed, so the controls' rise there is by construction)."""
    out: Dict[str, dict] = {}
    # APL at the boutons hides the step
    arms = {}
    for a in tuple(apl_arms) + (None,):
        v = APL_REFERENCE if a is None else variant_of(a)
        vals = mega.get(v, {})
        below = {w: bool(vals[w]['pn_bouton'] < no_step) for w in F.WIRINGS if w in vals}
        reasons = []
        if len(below) < len(F.WIRINGS):
            reasons.append('runs missing')
        if bad.get(v):
            reasons.append('holds fail: ' + ', '.join(bad[v]))
        if flagged.get(v):
            reasons.append('leaves a realism band on ' + ', '.join(flagged[v]))
        arms[v] = dict(n_below=int(sum(below.values())), below=below, reasons=reasons,
                       bouton_errors={w: float(vals[w]['pn_bouton']) for w in vals},
                       left_out=bool(leave_out_failed and a is not None and bad.get(v)),
                       meets=bool(not reasons and sum(below.values()) >= MIN_WIRINGS))
    ref = arms[APL_REFERENCE]
    left_out = [v for v in arms if arms[v]['left_out']]
    kept = [v for v in arms if v != APL_REFERENCE and not arms[v]['left_out']]
    met = [v for v in kept if arms[v]['meets']]
    complete = bool(kept) and all(not arms[v]['reasons'] for v in kept) and not ref['reasons']
    holds = bool(met) and not ref['meets'] and not ref['reasons']
    if holds:
        cls = 'holds'
    elif ref['reasons'] or not kept:
        cls = 'no class'
    elif ref['meets'] or complete:
        cls = 'does not hold'
    else:
        cls = 'no class'
    out['apl_hides_step'] = dict(holds=holds, arms_meeting=met, reference_meets=ref['meets'], arms=arms,
                                 left_out=left_out, kept=kept, no_step=no_step, min_wirings=MIN_WIRINGS,
                                 complete=complete, cls=cls)
    # the dendrites carry under the calyx that holds
    a, r = mega.get(variant_of(carry_arm), {}), mega.get(carry_reference, {})
    dn = F.compare({w: v['downstream'] for w, v in a.items()}, {w: v['downstream'] for w, v in r.items()},
                   spreads.get('downstream', np.nan))
    if dn.get('complete'):
        dn['improves_2se'] = bool(dn['change'] < -SE_MULTIPLE * dn['se'])
    st = {s: F.compare({w: v[s] for w, v in a.items()}, {w: v[s] for w, v in r.items()}, spreads.get(s, np.nan))
          for s in DOWNSTREAM}
    worse = [s for s, v in st.items() if v.get('worse_than_spread')]
    ctl = {}
    for c in control_arms:
        cv = mega.get(variant_of(c), {})
        ctl[variant_of(c)] = {s: dict(arm=mean_over(a, s, F.T_WIRINGS), control=mean_over(cv, s, F.T_WIRINGS))
                              for s in ('pn_dendrite',) + DOWNSTREAM}
        for s, x in ctl[variant_of(c)].items():
            x['below'] = bool(x['arm'] < x['control'])
    below_controls = all(ctl[c][s]['below'] for c in ctl for s in DOWNSTREAM)
    reasons = []
    names = [variant_of(carry_arm)] + [variant_of(c) for c in control_arms]
    if not dn.get('complete') or not all(v.get('complete') for v in st.values()) or \
            any(not np.isfinite(ctl[c][s]['control']) for c in ctl for s in DOWNSTREAM):
        reasons.append('runs missing')
    for v in names:
        tb = [t for t in bad.get(v, []) if t.rsplit('_', 1)[-1] in F.T_WIRINGS]
        if tb:
            reasons.append(f'{v}: holds fail: ' + ', '.join(tb))
        tf = [w for w in flagged.get(v, []) if w in F.T_WIRINGS]
        if tf:
            reasons.append(f'{v}: leaves a realism band on ' + ', '.join(tf))
    carry = bool(not reasons and dn.get('improves_2se') and not worse and below_controls)
    out['dendrites_carry'] = dict(holds=carry, arm=variant_of(carry_arm), reference=carry_reference, downstream=dn,
                                  stages=st, stage_worse=worse, controls=ctl, below_controls=below_controls,
                                  reasons=reasons, se_multiple=SE_MULTIPLE,
                                  cls='no class' if reasons else ('holds' if carry else 'does not hold'))
    return out


def step1_outcome(classes: Mapping[str, dict], bad: Mapping[str, Sequence[str]],
                  flagged: Mapping[str, Sequence[str]], carry_arm: str = CARRY_ARM) -> dict:
    """Step 1's outcomes as amended: (g) "not readable" when the amended main arm's hold fails (any of its runs, megamat
    or APL16, on any wiring) or it leaves a realism band (a flagged wiring of its APL16 runs); otherwise (d) when "the
    dendrites carry under the calyx that holds" holds, (e) when it does not, and none while it has no class for another
    reason (runs missing, a control's hold); (c) beside, when "APL at the boutons hides the step" holds. Step 0's (a),
    (b) and (f) are not read here."""
    v = variant_of(carry_arm)
    failed, bands = list(bad.get(v, [])), list(flagged.get(v, []))
    carry = classes['dendrites_carry']['cls']
    letters = []
    if failed or bands:
        letters.append('g')
    elif carry == 'holds':
        letters.append('d')
    elif carry == 'does not hold':
        letters.append('e')
    if classes['apl_hides_step']['cls'] == 'holds':
        letters.append('c')
    return dict(letters=letters, not_readable=bool(failed or bands), holds_failed=failed, bands_left=bands,
                dendrites_carry=carry, apl_hides_step=classes['apl_hides_step']['cls'],
                pending=bool(not (failed or bands) and carry == 'no class'))


def stopped_claims(root: Optional[Path] = None) -> Dict[str, str]:
    """The placeholder claims of the points not run (claims/<tag>/claim.json with lane 'stopped'): tag -> reason."""
    d = (out_dir() if root is None else Path(root)) / 'claims'
    out = {}
    for c in (sorted(d.glob('*/claim.json')) if d.exists() else []):
        r = _read_json(c)
        if r.get('lane') == 'stopped':
            out[str(r.get('tag') or c.parent.name)] = str(r.get('reason', ''))
    return out


def holds_table(r1: pd.DataFrame, tol: float, stopped: Mapping[str, str], wirings: Sequence[str]) -> pd.DataFrame:
    """One row per step-1 arm (STEP1_ARMS): its PN stage, whether it is the lock's, the runs read, those whose holds fail
    (`holds_failed`), those not run (`stopped_claims`), the bounds reached, the response rates against their targets and
    the hold's factors (least and largest), and its record: 'no class: a hold fails' (also for an arm whose APL-off run
    is a failed arm's, SAME_OFF_RUN), 'held on every wiring', 'stopped after <n> read', or '<n> of 8 read'."""
    def num(rr: pd.DataFrame, col: str, f) -> float:
        return float(f(pd.to_numeric(rr[col], errors='coerce'))) if col in rr and len(rr) else np.nan

    def truth(x) -> bool:
        return str(x).lower() in ('true', '1', '1.0')
    base, failed_of = [], {}
    for arm in STEP1_ARMS:
        tags = [tag_of(arm, w) for w in wirings]
        rr = r1[r1['tag'].astype(str).isin(tags)] if len(r1) else r1
        failed = holds_failed(rr, tol) if len(rr) else []
        failed_of[(arm.panel, arm.name)] = failed
        base.append((arm, tags, rr, failed))
    rows = []
    for arm, tags, rr, failed in base:
        st = [x for x in tags if x in stopped]
        same = SAME_OFF_RUN.get(arm.name) if arm in LOCK_ARMS else None
        if failed:
            rec = 'no class: a hold fails'
        elif same and failed_of.get((arm.panel, same)):
            rec = f"no class: a hold fails (its APL-off run is {arm.panel}_{same}'s)"
        elif len(rr) == len(tags):
            rec = 'held on every wiring'
        elif st and len(rr) + len(st) == len(tags):
            rec = f'stopped after {len(rr)} read'
        else:
            rec = f'{len(rr)} of {len(tags)} read'
        rows.append(dict(panel=arm.panel, arm=arm.name, variant=variant_of(arm.name), pn_stage=arm.pn_stage,
                         lock=arm in LOCK_ARMS, aplpn=arm.aplpn, kckc=arm.kckc, n_wirings=len(tags), n_read=len(rr),
                         n_failed=len(failed), failed=';'.join(failed), n_stopped=len(st), stopped=';'.join(st),
                         n_at_bound_thr=int(sum(truth(x) for x in rr['at_bound_thr'])) if 'at_bound_thr' in rr else 0,
                         n_at_bound_apl=int(sum(truth(x) for x in rr['at_bound_apl'])) if 'at_bound_apl' in rr else 0,
                         rate_off_min=num(rr, 'rate_off', np.min), rate_off_max=num(rr, 'rate_off', np.max),
                         target_off_min=num(rr, 'frac_off_target', np.min),
                         target_off_max=num(rr, 'frac_off_target', np.max),
                         rate_on_min=num(rr, 'rate_on', np.min), rate_on_max=num(rr, 'rate_on', np.max),
                         target_on_min=num(rr, 'frac_on_target', np.min), target_on_max=num(rr, 'frac_on_target', np.max),
                         thr_factor_min=num(rr, 'thr_factor', np.min), thr_factor_max=num(rr, 'thr_factor', np.max),
                         apl_factor_min=num(rr, 'apl_factor', np.min), apl_factor_max=num(rr, 'apl_factor', np.max),
                         record=rec))
    return pd.DataFrame(rows)


def bouton_beside(mega: Mapping[str, Mapping[str, Dict[str, float]]], no_step: float, F,
                  carry_arm: str = CARRY_ARM, reference: str = CARRY_REFERENCE) -> dict:
    """Beside, classing nothing: the amended main arm's boutons against the no-step yardstick (each wiring's bouton
    strength error against the held-out half; the wirings below it, on T and on all eight; the mean over T), with
    A1's own beside (already below it on T: it does not separate the dendrites from the normalization there)."""
    out = {}
    for v in (variant_of(carry_arm), reference):
        vals = mega.get(v, {})
        errs = {w: float(vals[w]['pn_bouton']) for w in F.WIRINGS if w in vals}
        out[v] = dict(bouton_errors=errs, n_below_T=int(sum(errs[w] < no_step for w in F.T_WIRINGS if w in errs)),
                      n_below_all=int(sum(e < no_step for e in errs.values())), n_wirings=len(errs),
                      mean_T=mean_over(vals, 'pn_bouton', F.T_WIRINGS) if vals else np.nan)
    return dict(no_step=no_step, variants=out)


def e143_beside(spreads: Mapping[str, float], scored_cls: Optional[str], carry_arm: str = CARRY_ARM,
                control_arms: Sequence[str] = CONTROL_ARMS) -> dict:
    """Beside, classing nothing (the amendment's): E143's interval over resampled flies on the amended main arm's change
    from A1, at each megamat stage, over the boutons, claws and somata ('downstream') and over all four ('stage'), with
    the share of resamples on which "the dendrites carry under the calyx that holds" returns its scored class.

    E143's reading: the arm's and its controls' T runs read once against all the flies (E143's `read_runs`: E134's
    default readings, E136's seed index); A1 its saved views (E143's setting E143_A1); each of E143's resamples of a
    stage's flies (`resamples`: N_BOOT, its seed, the stages drawn from one generator) a reference, in all the flies'
    one-fly units, the unresampled flies the last; the change the mean over T of the arm minus A1 on the same wiring,
    its interval E143's CI over the resamples; the class on each reference E136's `compare` with the scored spreads
    (the family's over its eight wirings) and the controls' means over T. Needs E143's committed lock and its saved
    views and flies (outputs/E143_megamat_fly_sampling/views/, on the host that read them); otherwise 'status' says why
    it is not given."""
    from experiments import E143_megamat_fly_sampling as S
    lock = S.require_lock()
    a1_path = S.views_path(E143_A1, S.PANEL)
    if not a1_path.exists():
        return dict(status=f'not given: {a1_path} is not on this host (E143 read it on hal)')
    names = [variant_of(carry_arm)] + [variant_of(c) for c in control_arms]
    reads = {}
    for v in names:
        runs = {}
        for w in S.T_WIRINGS:
            try:
                runs[w] = load_point(f'{S.PANEL}_{v[len(E142_PREFIX):]}_{w}')
            except FileNotFoundError:
                pass
        if len(runs) < len(S.T_WIRINGS):
            return dict(status=f'not given: {v} has {len(runs)} of {len(S.T_WIRINGS)} T runs saved')
        reads[v] = S.read_runs(v, S.PANEL, runs)
    a1 = pd.read_parquet(a1_path)
    draws = S.resamples(lock)
    vals: Dict[str, Dict[str, np.ndarray]] = {}
    for st in S.LOSS_STAGES:
        key = S.megamat_key(st)
        if not S.flies_path(key).exists():
            return dict(status=f'not given: {S.flies_path(key)} is not on this host')
        fl = pd.read_parquet(S.flies_path(key))
        model = pd.concat([a1[a1['key'] == key]] + [reads[v][reads[v]['key'] == key] for v in names],
                          ignore_index=True)
        sd = S.side_from_tables(key, fl, model)
        S.check_side_flies(sd, S.lock_flies(lock, key))
        n = len(sd.flies)
        W = np.vstack([S.boot_weights(n, draws[key]), np.full((1, n), 1.0 / n)])
        e = S.side_errors(sd, W, ('strength',), S.stage_error_col(key, ('strength',)))
        vals[st] = {nm: S.per_wiring(e, sd.rows, nm) for nm in [E143_A1, *names]}

    def rms(stages: Sequence[str], nm: str) -> np.ndarray:
        return np.sqrt(np.mean(np.square([vals[s][nm] for s in stages]), axis=0))
    stat = dict(stage={nm: rms(S.LOSS_STAGES, nm) for nm in [E143_A1, *names]},
                downstream={nm: rms(DOWNSTREAM, nm) for nm in [E143_A1, *names]}, **vals)
    main = names[0]
    changes = {}
    for k, d in stat.items():
        ch = (d[main] - d[E143_A1]).mean(axis=1)
        lo, hi = np.quantile(ch[:-1], S.CI)
        changes[k] = dict(all_flies=float(ch[-1]), lo=float(lo), hi=float(hi), excludes_zero=bool(lo > 0 or hi < 0))
    nT = len(S.T_WIRINGS)

    def change(k: str) -> Tuple[np.ndarray, np.ndarray]:
        diff = stat[k][main] - stat[k][E143_A1]
        return diff.mean(axis=1), diff.std(axis=1, ddof=1)
    m, s = change('downstream')
    sp = float(spreads.get('downstream', np.nan))
    se = np.maximum(sp / 2.0 if np.isfinite(sp) else 0.0, s / np.sqrt(nT))
    improves = m < -SE_MULTIPLE * se
    worse = np.any([change(st)[0] > float(spreads.get(st, np.nan)) for st in DOWNSTREAM
                    if np.isfinite(float(spreads.get(st, np.nan)))] or [np.zeros(len(m), bool)], axis=0)
    below = np.all([stat[st][main].mean(axis=1) < stat[st][c].mean(axis=1) for st in DOWNSTREAM for c in names[1:]],
                   axis=0)
    holds = improves & ~worse & below
    cls = np.where(holds, 'holds', 'does not hold')
    share = float(np.mean(cls[:-1] == scored_cls)) if scored_cls in ('holds', 'does not hold') else np.nan
    return dict(status='given', n_boot=int(len(cls) - 1), ci=list(S.CI), boot_seed=S.BOOT_SEED,
                reference="E143's resamples of all the stage's flies (the last reference all the flies, unresampled)",
                arm=main, a1=E143_A1, controls=names[1:], changes=changes, class_all_flies=str(cls[-1]),
                share_holds=float(np.mean(holds[:-1])), class_scored=scored_cls, share_scored=share,
                e143_lock_sha256=lock.get('sha256'))


def part_score() -> dict:
    """Step 1's score (the entry's design item 3; the step-1 classes as amended 2026-11-26 13:59, with the lock's own
    reading kept as the record): every step-1 arm and E135's three saved
    imposed arms at the family's calyx registered with E136's score functions (`load_variant_run`), read as E136 reads
    its arms under the adopted measure (strength by odor), with the unrefitted family and A1 beside: the megamat stages
    on the eight wirings against the held-out half and all flies (score_step1_megamat.csv; per wiring, 'mean_S',
    'mean_T' and 'common' the wirings' mean views); the APL16 somata's KC population errors (score_step1_panels.csv,
    E136's `score_panels`); the APL16 rows with the realism bands and the APL effect, its calibrated parts beside (E136's
    `apl16_rows`, E139's `apl_effect_parts`; score_step1_apl16.csv); the no-step yardstick of the boutons against the
    held-out half (`no_step_yardstick`, checked against step 0's yardsticks.csv where it exists); the classes as amended
    (`classify_step1`) and the lock's own reading of them ('lock_reading'), with the holds (`holds_failed`, E136's
    HOLD_TOL) and the bands; the outcomes (`step1_outcome`, (g) among them); the holds table (`holds_table`,
    step1_holds.csv); beside, classing nothing, the amended main arm's boutons against the no-step yardstick
    (`bouton_beside`) and E143's interval over resampled flies on its change from A1 (`e143_beside`).
    step1_classes.json and provenance_score_step1.json."""
    from experiments import E136_stage_population_fit as F
    from experiments import E139_normalization_ladder as L
    t0 = time.time()
    E135 = _e135()
    old = F.MEASURE
    F.set_measure(E136_MEASURE)
    try:
        F.require_gate()
        stages = tuple(F.LOSS_STAGES)
        names = [*E135_VARIANTS, *step1_variants()]
        for v in names:
            F.register_variant(v, load_variant_run)
        allv = [F.FAMILY, A1, *names]
        mrows, mega = [], {}
        for v in allv:
            for half in ('held', None):
                sc = F.score_megamat(v, stages, half, F.WIRINGS)
                if sc is None:
                    continue
                vals = F.megamat_values(sc, stages, F.WIRINGS)
                if half == 'held':
                    mega[v] = vals
                base = dict(variant=v, half=half or 'all')
                for w, x in vals.items():
                    mrows.append(dict(base, wiring=w, **x))
                for lab, ws in (('mean_S', F.S_WIRINGS), ('mean_T', F.T_WIRINGS)):
                    mrows.append(dict(base, wiring=lab, **pd.DataFrame([vals[w] for w in ws]).mean().to_dict()))
                mrows.append(dict(base, wiring='common', stage=sc['stage_error'],
                                  downstream=sc['stage_error_downstream'], **{s: sc[f'{s}_error'] for s in stages}))
            print(f'E142 score: {v} ({time.time() - t0:.0f} s)', flush=True)
        megadf = pd.DataFrame(mrows)
        megadf.to_csv(out_dir() / 'score_step1_megamat.csv', index=False)
        fam = mega.get(F.FAMILY, {})
        spreads = dict(stage=F._sd(x['stage'] for x in fam.values()),
                       downstream=F._sd(x['downstream'] for x in fam.values()),
                       **{s: F._sd(x[s] for x in fam.values()) for s in stages})
        prow = [F.score_panels(v) for v in allv]
        panels = pd.concat([p for p in prow if len(p)], ignore_index=True) if any(len(p) for p in prow) else \
            pd.DataFrame()
        panels.to_csv(out_dir() / 'score_step1_panels.csv', index=False)
        apl_rows, flagged = [], {}
        fam_flags = {}
        for v in allv:
            r = F.apl16_rows(v)
            if not len(r):
                continue
            parts = {}
            for w in r.wiring:
                run = F.load_read(F.ReadPoint(v, 'apl16', w, None))
                parts[w] = L.apl_effect_parts(L.run_outputs(run)) if run is not None else {}
            r = r.assign(variant=v)
            for k in sorted({k for d in parts.values() for k in d}):
                r[f'cal_{k}'] = r.wiring.map({w: d.get(k, np.nan) for w, d in parts.items()})
            apl_rows.append(r)
            fl = F.realism_of(r)
            if v == F.FAMILY:
                fam_flags = fl
            flagged[v] = [w for w, f in fl.items() if f and not fam_flags.get(w)]
        apl = pd.concat(apl_rows, ignore_index=True) if apl_rows else pd.DataFrame()
        apl.to_csv(out_dir() / 'score_step1_apl16.csv', index=False)
        # the APL16 arm's flags stand for the model of its factors: imposed_a1 on megamat is the same model
        no_step = no_step_yardstick(F)
        check = None
        ys = out_dir() / 'yardsticks.csv'
        if ys.exists():
            y = pd.read_csv(ys)
            y = y[(y.stage == 'pn_bouton') & (y.glomeruli == 'dendrites') & (y.fly_set == 'held')]
            if len(y):
                check = float(y.no_step.iloc[0])
                if not np.isclose(check, no_step, rtol=0, atol=1e-9):
                    raise RuntimeError(f"the no-step yardstick {no_step} is not step 0's {check}")
        tol = float(F.HOLD_TOL)
        r1, r2 = runs1_table(), E135.runs2_table()
        bad = {}
        for v in names:
            if v in E135_VARIANTS:
                rows = r2[r2.tag.astype(str).isin([f'megamat_{E135_VARIANTS[v]}_{w}' for w in F.WIRINGS])]
            else:
                arm = v[len(E142_PREFIX):]
                rows = r1[r1.tag.astype(str).isin([f'{p}_{arm}_{w}' for p in ('megamat', 'apl16') for w in F.WIRINGS])]
            bad[v] = holds_failed(rows, tol) if len(rows) else []
        classes = classify_step1(mega, spreads, no_step, bad, flagged, F)
        lock_reading = classify_step1(mega, spreads, no_step, bad, flagged, F, apl_arms=LOCK_APL_ARMS,
                                      carry_arm=LOCK_CARRY_ARM, control_arms=LOCK_CONTROL_ARMS, leave_out_failed=False)
        outcome = step1_outcome(classes, bad, flagged)
        stopped = stopped_claims()
        holds = holds_table(r1, tol, stopped, F.WIRINGS)
        holds.to_csv(out_dir() / 'step1_holds.csv', index=False)
        beside = dict(boutons=bouton_beside(mega, no_step, F))
        try:
            beside['e143'] = e143_beside(spreads, classes['dendrites_carry']['cls'])
        except Exception as e:      # beside only: the score stands without it
            beside['e143'] = dict(status=f'not given: {type(e).__name__}: {e}')
        rec = dict(classes=classes, outcome=outcome, lock_reading=lock_reading, beside=beside,
                   holds_table=holds.to_dict('records'), stopped=stopped, spreads=spreads, no_step=no_step,
                   no_step_step0=check, holds_tol=tol, holds_failed=bad, flagged=flagged, variants=allv,
                   stages=list(stages), measure=E136_MEASURE, a1=a1_factors(), a1_point=a1_point(),
                   time=_now(), host=_host(), seconds=round(time.time() - t0, 1),
                   **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
        _write_json(out_dir() / 'step1_classes.json', rec)
        _write_json(out_dir() / 'provenance_score_step1.json',
                    dict(part='score', script=SCRIPT, time=_now(), host=_host(), code=H.code_commits(),
                         e135_lock_sha256=E135._sha256(E135.LOCK_COPY) if E135.LOCK_COPY.exists() else None))
        print(json.dumps(dict({k: v['cls'] for k, v in classes.items()}, outcome=outcome['letters'],
                              lock_reading={k: v['cls'] for k, v in lock_reading.items()},
                              e143=beside['e143'].get('status')), indent=1), flush=True)
        return rec
    finally:
        F.set_measure(old)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True,
                    choices=('step1', 'step1_dry', 'score', 'variants', 'routes_check', 'routes_summary'))
    ap.add_argument('--lane')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--wiring', default='real')
    ap.add_argument('--wirings', nargs='*', help="routes_check: the wirings to read (default all eight)")
    ap.add_argument('--tags', nargs='*', default=None,
                    help="step1: run only these points (e.g. a first read, megamat_a1stage_imposed_P3); step1_dry: "
                         "only these arms' tables (on --wiring)")
    a = ap.parse_args(argv)
    if a.part == 'step1':
        part_step1(a.lane, a.reverse, a.reclaim, a.tags)
    elif a.part == 'step1_dry':
        part_step1_dry(a.wiring, a.tags)
    elif a.part == 'routes_check':
        part_routes_check(a.wirings, a.lane)
    elif a.part == 'routes_summary':
        part_routes_summary()
    elif a.part == 'score':
        part_score()
    else:
        write_variants()


if __name__ == '__main__':
    main()
