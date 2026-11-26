#!/usr/bin/env python3
"""E142's step 1 (design item 2): the recorded PN dendrites imposed on the model's PN stage under the calyx of A1, the
one stage-fit arm that holds, at fixed scales with both response rates re-held, 56 runs on the family's eight wirings.

The entry is item E142 of experiment_priorities.md (proposed and moved in by viki 2026-11-26 05:10 under the idle-host
rule; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26 04:48). Step 0 (no simulation) is
viki's, in experiments/E142_dendrite_bouton_step.py, whose `--part step1` calls `part_step1` here; this module holds
step 1's runner (hal's claim, 2026-11-26 05:12). It is E135's step-2 runner (the antennal lobe held to the recorded PN
dendrites, experiments/E135_antennal_lobe_dendrites.py, `run_point2`) with the calyx's two scales changed before the
hold, and it calls E135's functions rather than copying them.

The arms (STEP1_ARMS; the entry's design item 2), each on the eight wirings:
- on the megamat input, the dendrites imposed (E135's 'imposed' table) with A1's calyx factors (APL->bouton and KC-KC);
  with each factor alone; with KC-KC at A1's and APL->bouton at none (its scale 0, in both conditions); with A1's
  factors and the table's glomerulus labels permuted, and with its odor labels permuted (E135's 'imposed_glom_perm'
  and 'imposed_odor_perm', each wiring's one seeded permutation, `E135.arm_seed`, so that each pairs with E135's own
  permuted run of the wiring);
- on the APL16 input, the dendrites imposed for its 12 megamat odors with A1's factors (E135's APL16 'imposed' table:
  paa, euc, -aPine and B-myr keep the model's PNs), for the bands and the APL effect.

Readings taken where the entry leaves them open (the narrowest; each is also in the docstrings):
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

The guards: E134's measure gate (E135's `family_stage` checks it); E136's fit of A1 (`a1_factors`); and the
authorization of step 1's runs (`require_step1_allowed`): E135's runs entered the dendrites by olfsysm's `pn_rates`
with Yang's OK of 2026-11-24 09:43, and viki asks Yang whether that OK covers these runs (the entry's Status), so no
run starts until results/E142_dendrite_bouton_step/step1_authorized.json records the answer (written by the session
that relays it; `write_step1_authorization`). A dry run (`--part step1_dry`) builds every table on one wiring without
any simulation of the mushroom body.

Outputs: <shared>/results/E142_dendrite_bouton_step/ (runs_1_<host>_<lane>.csv, provenance_1_<host>_<lane>.json,
identity_1_<host>_<lane>.json, claims/, step1_dry_<wiring>.json); outputs/E142_dendrite_bouton_step/step1/ (local,
large: the runs, also in the archive's data/model_outputs/E142_dendrite_bouton_step/step1/).

Usage (from the repository root, as modules):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E142_imposed_calyx_runs --part step1 --lane <host>[-<n>] [--reverse]
    ./venv/bin/python -u -m experiments.E142_imposed_calyx_runs --part step1_dry [--wiring real]
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

# an arm of step 1: its input panel, its name, the E135 table it is built as, and the two calyx factors, each 'A1' (A1's
# factor), 'family' (x1, the wiring's own fit's scale) or 'none' (the scale 0)
Step1Arm = namedtuple('Step1Arm', 'panel name table aplpn kckc')
STEP1_ARMS: Tuple[Step1Arm, ...] = (
    Step1Arm('megamat', 'imposed_a1', 'imposed', 'A1', 'A1'),
    Step1Arm('megamat', 'imposed_aplpn', 'imposed', 'A1', 'family'),
    Step1Arm('megamat', 'imposed_kckc', 'imposed', 'family', 'A1'),
    Step1Arm('megamat', 'imposed_kckc_aplpn_none', 'imposed', 'none', 'A1'),
    Step1Arm('megamat', 'imposed_a1_glom_perm', 'imposed_glom_perm', 'A1', 'A1'),
    Step1Arm('megamat', 'imposed_a1_odor_perm', 'imposed_odor_perm', 'A1', 'A1'),
    Step1Arm('apl16', 'imposed_a1', 'imposed', 'A1', 'A1'),
)
FACTOR_PARAMS: Dict[str, str] = {'aplpn': 'wAPLPN_scale', 'kckc': 'wKCKC_scale'}
FACTOR_CHOICES: Tuple[str, ...] = ('A1', 'family', 'none')
# what each arm is, for the archive's params.csv
ARM_KINDS: Dict[str, str] = {
    'imposed_a1': "the recorded PN dendrites imposed (E135's table, route 1) with A1's APL->bouton and KC-KC factors",
    'imposed_aplpn': "the recorded PN dendrites imposed with A1's APL->bouton factor alone (KC-KC the family's)",
    'imposed_kckc': "the recorded PN dendrites imposed with A1's KC-KC factor alone (APL->bouton the family's)",
    'imposed_kckc_aplpn_none': "the recorded PN dendrites imposed with A1's KC-KC factor and APL->bouton at none (0)",
    'imposed_a1_glom_perm': "the recorded PN dendrites imposed, glomerulus labels permuted, with A1's two factors",
    'imposed_a1_odor_perm': "the recorded PN dendrites imposed, odor labels permuted, with A1's two factors",
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


def planned_step1(reverse: bool = False) -> List[Tuple[Step1Arm, str]]:
    """Step 1's 56 points (arm, wiring), the arms in STEP1_ARMS' order, each on E135's eight wirings (E134's S then
    T); reversed with `reverse`."""
    pts = [(a, w) for a in STEP1_ARMS for w in _e135().WIRINGS]
    return pts[::-1] if reverse else pts


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


def arm_factors(arm: Step1Arm, a1: Mapping[str, float]) -> Dict[str, float]:
    """The arm's factor on each of the two calyx scales: A1's ('A1'), 1 ('family') or 0 ('none')."""
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


class Step1:
    """What step 1's points share in one process: E135's lock (its 28 shared glomeruli), E134's lock, A1's factors, the
    family's PN stage per (panel, wiring) (built once each, the lane's identity on the wiring checked on that build,
    the build then freed) and the dendrite means in E135's fit units."""

    def __init__(self, lane: str, lock134: Optional[Mapping] = None, a1: Optional[Mapping] = None,
                 identity: bool = True):
        E135 = _e135()
        self.lane = lane
        self.lock = E135.read_lock()
        self.glomeruli = [str(g) for g in self.lock['glomeruli']]
        self.lock134 = lock134
        self.a1 = dict(a1) if a1 is not None else a1_factors()
        self.identity = identity
        self._stages: Dict[Tuple[str, str], Tuple[object, dict]] = {}
        self._den: Optional[pd.DataFrame] = None

    def stage(self, panel: str, wiring: str):
        key = (panel, wiring)
        if key not in self._stages:
            stage, finfo = _e135().family_stage(panel, wiring)
            ident = lane_identity(panel, wiring, self.lane, finfo) if self.identity else dict(skipped=True)
            port = finfo.get('identity')
            finfo = None
            gc.collect()
            self._stages[key] = (stage, dict(port=port, lane_identity=ident))
        return self._stages[key]

    def dendrites(self) -> pd.DataFrame:
        if self._den is None:
            E135 = _e135()
            self._den = E135.dendrite_means(E135.FIT_UNITS, self.glomeruli)
        return self._den


def point_table(arm: Step1Arm, wiring: str, ctx: Step1) -> Tuple[pd.DataFrame, List[str], dict, object, dict]:
    """A point's table of PN rates as E135 builds it for the arm's E135 table on the wiring (`E135.arm_table` with
    E135's seed for the wiring, clipped), the glomeruli it sets, its record, the family's PN stage and the stage's
    record."""
    E135 = _e135()
    stage, srec = ctx.stage(arm.panel, wiring)
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
    `run_point2` for the arm's E135 table, route 1, with the calyx factors applied before the hold (module docstring)."""
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
        table, set_gl, trec, stage, srec = point_table(arm, wiring, ctx)
        set_od = list(trec['set_odors'])
        kw, rrec = E135.arm_input(stage, table, set_gl, 1, set_odors=set_od)
        orn = E135.panel_orn(arm.panel)
        extra = dict(pn_rates=E135.input_columns(kw['extra_kws']['pn_rates'], orn))
        fam = E135.family_parts(wiring)
        factors = arm_factors(arm, ctx.a1)
        p0 = arm_params(fam['p'], factors)
        tg = E135.step2_targets(arm.panel, wiring, orn.shape[1], fam, ctx.lock134)
        h = hold(wiring, orn, fam, p0, extra, tg, f'E142 {tag}')
        params, params_off = dict(h['params']), dict(h['params_off'])
        for k, name in FACTOR_PARAMS.items():
            if float(params[name]) != float(p0[name]):
                raise RuntimeError(f'{tag}: the hold moved {name} ({p0[name]} -> {params[name]})')
        built = E135._build(arm.panel, wiring, params, None, extra)
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
                    set_glomeruli=list(set_gl), set_odors=set_od, table_sha256=table_sha, stages=list(built.stages),
                    params=params, params_off=params_off, targets=tg, host=_host(), lane=ctx.lane, time=_now())
        E135.save_point(tag, out, built, info, table, dirs=step1_dirs())
        clip = trec.get('clip', {})
        row = dict(tag=tag, panel=arm.panel, arm=arm.name, table_arm=arm.table, wiring=wiring, route=1,
                   units=E135.FIT_UNITS, seed=trec['seed'], aplpn=arm.aplpn, kckc=arm.kckc,
                   aplpn_factor=factors['aplpn'], kckc_factor=factors['kckc'],
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
                                                              n_sims=h['n_sims']),
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


def part_step1(lane: Optional[str] = None, reverse: bool = False, reclaim: bool = False) -> pd.DataFrame:
    """Step 1, resumable, one point claimed at a time (`planned_step1`, `run_point1`); every lane on any host runs the
    whole list and skips what another has claimed or done. Refuses before the authorization (`require_step1_allowed`),
    E136's fit of A1 (`a1_factors`) and E134's gate (`E135.family_stage`). Rows in runs_1_<host>_<lane>.csv, provenance
    in provenance_1_<host>_<lane>.json, identities in identity_1_<host>_<lane>.json, the runs in `step1_dirs`."""
    from experiments import E134_population_baseline as E134R
    E135 = _e135()
    require_step1_allowed()
    E135.require_e134_gate()
    E134R.check_wirings()
    if tuple(E135.WIRINGS) != tuple(E134R.WIRINGS):
        raise RuntimeError(f"E135's wirings {E135.WIRINGS} are not E134's {E134R.WIRINGS}")
    lock134 = E135.require_e134_lock()
    lane = lane or _host()
    ctx = Step1(lane, lock134)
    pts = planned_step1(reverse)
    table = out_dir() / f'runs_1_{_host()}_{lane}.csv'
    prov = out_dir() / f'provenance_1_{_host()}_{lane}.json'
    print(f"E142 step 1, lane {lane}: {len(pts)} points; A1's factors APL->bouton x{ctx.a1['aplpn']:.4f}, KC-KC "
          f"x{ctx.a1['kckc']:.4f}", flush=True)
    for arm, w in pts:
        done = set(runs1_table()['tag'].astype(str))
        run_point1(arm, w, ctx, table, prov, reclaim, done)
    return runs1_table()


def part_step1_dry(wiring: str = 'real') -> dict:
    """Every arm's table on one wiring, without any simulation of the mushroom body (the family's PN stage is built,
    its olfsysm layers checked against E135's port, and no lane identity is run): per arm the set glomeruli and odors,
    the clip's counts and the arm's absolute calyx scales from the wiring's fit; step1_dry_<wiring>.json. Needs E134's
    gate and A1's fit, not the authorization."""
    E135 = _e135()
    ctx = Step1('dry', identity=False)
    fam = E135.family_parts(wiring)
    out = dict(wiring=wiring, a1=ctx.a1, time=_now(), host=_host(),
               **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))}, arms={})
    for arm in STEP1_ARMS:
        table, set_gl, trec, _, _ = point_table(arm, wiring, ctx)
        f = arm_factors(arm, ctx.a1)
        p = arm_params(fam['p'], f)
        clip = trec.get('clip', {})
        out['arms'][tag_of(arm, wiring)] = dict(
            n_set_glomeruli=len(set_gl), set_odors=list(trec['set_odors']), n_mapped=trec.get('n_mapped'),
            clip_n_set=clip.get('n_set'), clip_n_clipped=clip.get('n_clipped'), clip_named=clip.get('named'),
            factors=f, **{name: p[name] for name in FACTOR_PARAMS.values()},
            table_min_hz=float(np.nanmin(table.loc[set_gl].to_numpy(dtype=float))),
            table_max_hz=float(np.nanmax(table.loc[set_gl].to_numpy(dtype=float))))
    _write_json(out_dir() / f'step1_dry_{wiring}.json', out)
    print(json.dumps(out['arms'], indent=1, default=str), flush=True)
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
    """The tables step 1 uses, one arm per distinct (panel, E135 table): megamat's imposed table, its two permuted
    tables (each wiring's own seeded permutation) and APL16's imposed table."""
    seen, out = set(), []
    for arm in STEP1_ARMS:
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
    the odors without a fixed point, into routes_check.json."""
    from experiments import E136_stage_population_fit as F
    E135 = _e135()
    d = out_dir()
    tabs = [pd.read_csv(p) for p in sorted(d.glob('routes_check_1_*.csv')) if '_odors_' not in p.name]
    ods = [pd.read_csv(p) for p in sorted(d.glob('routes_check_1_odors_*.csv'))]
    tables = pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame()
    odors = pd.concat(ods, ignore_index=True) if ods else pd.DataFrame()
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


def write_variants(root: Optional[Path] = None) -> pd.DataFrame:
    """step1_variants.csv: one row per saved step-1 run (`runs1_table`), for the model-output archive's params.csv:
    the tag, panel, arm, the E135 table, wiring, seed, the factors, what the arm is (ARM_KINDS) and the absolute
    parameters of both conditions (E135's VARIANT_COLUMNS)."""
    E135 = _e135()
    d = out_dir() if root is None else Path(root)
    rows = []
    for r in runs1_table(d).to_dict('records'):
        rows.append(dict(tag=r['tag'], variant=r['tag'], part=PART, subdir=PART, panel=r.get('panel'), arm=r.get('arm'),
                         table_arm=r.get('table_arm'), model=E135.FAMILY, wiring=r.get('wiring'), route=1,
                         seed=r.get('seed'), aplpn_factor=r.get('aplpn_factor'), kckc_factor=r.get('kckc_factor'),
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
# the lock's step-1 classes (experiments/E142_dendrite_bouton_step.py CLASSES, read from the lock at score time)
APL_ARMS: Tuple[str, ...] = ('imposed_a1', 'imposed_aplpn', 'imposed_kckc_aplpn_none')
CARRY_ARM = 'imposed_a1'
CONTROL_ARMS: Tuple[str, ...] = ('imposed_a1_glom_perm', 'imposed_a1_odor_perm')
APL_REFERENCE = 'E135_imposed'
MIN_WIRINGS = 6
SE_MULTIPLE = 2.0


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
                   bad: Mapping[str, Sequence[str]], flagged: Mapping[str, Sequence[str]], F) -> dict:
    """Step 1's classes as the lock fixes them, from the megamat values against the held-out half (`mega`: per variant,
    per wiring, E136's `megamat_values`), the family's spreads over its eight wirings, the no-step yardstick of the
    boutons against the held-out half, the runs whose holds failed (`bad`: per variant, wirings) and the wirings whose
    APL16 run leaves a realism band where the family's does not (`flagged`: per variant).

    - 'apl_hides_step': in at least one of APL_ARMS, the boutons' strength error is below the no-step yardstick on
      MIN_WIRINGS of the eight wirings, and in E135's imposed arm (APL->bouton at the family's) it is not; an arm with a
      failed hold or a flagged wiring gives no class.
    - 'dendrites_carry': with A1's factors (CARRY_ARM), the error over the boutons, claws and somata ('downstream') is
      below A1's on T by more than SE_MULTIPLE standard errors (E136's `compare`: the larger of half the spread and the
      SD of the per-wiring changes over sqrt(4)), none of the three stages is worse than A1's by more than its spread,
      and each of the three is below both permuted controls (the mean over T); no class with a failed hold or a flagged
      wiring among the arm's or the controls' runs on T.
    The PN dendrites are given beside (they are imposed, so the controls' rise there is by construction)."""
    out: Dict[str, dict] = {}
    # APL at the boutons hides the step
    arms = {}
    for a in APL_ARMS + (None,):
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
                       meets=bool(not reasons and sum(below.values()) >= MIN_WIRINGS))
    ref = arms[APL_REFERENCE]
    met = [v for v in arms if v != APL_REFERENCE and arms[v]['meets']]
    complete = all(not arms[v]['reasons'] for v in arms if v != APL_REFERENCE) and not ref['reasons']
    holds = bool(met) and not ref['meets'] and not ref['reasons']
    out['apl_hides_step'] = dict(holds=holds, arms_meeting=met, reference_meets=ref['meets'], arms=arms,
                                 no_step=no_step, min_wirings=MIN_WIRINGS, complete=complete,
                                 cls='no class' if (not holds and not complete and not met) else
                                 ('holds' if holds else 'does not hold'))
    # the dendrites carry under the calyx that holds
    a, r = mega.get(variant_of(CARRY_ARM), {}), mega.get('A1', {})
    dn = F.compare({w: v['downstream'] for w, v in a.items()}, {w: v['downstream'] for w, v in r.items()},
                   spreads.get('downstream', np.nan))
    if dn.get('complete'):
        dn['improves_2se'] = bool(dn['change'] < -SE_MULTIPLE * dn['se'])
    st = {s: F.compare({w: v[s] for w, v in a.items()}, {w: v[s] for w, v in r.items()}, spreads.get(s, np.nan))
          for s in DOWNSTREAM}
    worse = [s for s, v in st.items() if v.get('worse_than_spread')]
    ctl = {}
    for c in CONTROL_ARMS:
        cv = mega.get(variant_of(c), {})
        ctl[variant_of(c)] = {s: dict(arm=mean_over(a, s, F.T_WIRINGS), control=mean_over(cv, s, F.T_WIRINGS))
                              for s in ('pn_dendrite',) + DOWNSTREAM}
        for s, x in ctl[variant_of(c)].items():
            x['below'] = bool(x['arm'] < x['control'])
    below_controls = all(ctl[c][s]['below'] for c in ctl for s in DOWNSTREAM)
    reasons = []
    names = [variant_of(CARRY_ARM)] + [variant_of(c) for c in CONTROL_ARMS]
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
    out['dendrites_carry'] = dict(holds=carry, downstream=dn, stages=st, stage_worse=worse, controls=ctl,
                                  below_controls=below_controls, reasons=reasons, se_multiple=SE_MULTIPLE,
                                  cls='no class' if reasons else ('holds' if carry else 'does not hold'))
    return out


def part_score() -> dict:
    """Step 1's score (the entry's design item 3; the lock's step-1 classes): every step-1 arm and E135's three saved
    imposed arms at the family's calyx registered with E136's score functions (`load_variant_run`), read as E136 reads
    its arms under the adopted measure (strength by odor), with the unrefitted family and A1 beside: the megamat stages
    on the eight wirings against the held-out half and all flies (score_step1_megamat.csv; per wiring, 'mean_S',
    'mean_T' and 'common' the wirings' mean views); the APL16 somata's KC population errors (score_step1_panels.csv,
    E136's `score_panels`); the APL16 rows with the realism bands and the APL effect, its calibrated parts beside (E136's
    `apl16_rows`, E139's `apl_effect_parts`; score_step1_apl16.csv); the no-step yardstick of the boutons against the
    held-out half (`no_step_yardstick`, checked against step 0's yardsticks.csv where it exists); the classes
    (`classify_step1`) with the holds (`holds_failed`, E136's HOLD_TOL) and the bands. step1_classes.json and
    provenance_score_step1.json."""
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
        rec = dict(classes=classes, spreads=spreads, no_step=no_step, no_step_step0=check, holds_tol=tol, holds_failed=bad,
                   flagged=flagged, variants=allv, stages=list(stages), measure=E136_MEASURE, a1=a1_factors(),
                   time=_now(), host=_host(), seconds=round(time.time() - t0, 1),
                   **{k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))})
        _write_json(out_dir() / 'step1_classes.json', rec)
        _write_json(out_dir() / 'provenance_score_step1.json',
                    dict(part='score', script=SCRIPT, time=_now(), host=_host(), code=H.code_commits(),
                         e135_lock_sha256=E135._sha256(E135.LOCK_COPY) if E135.LOCK_COPY.exists() else None))
        print(json.dumps({k: v['cls'] for k, v in classes.items()}, indent=1), flush=True)
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
    a = ap.parse_args(argv)
    if a.part == 'step1':
        part_step1(a.lane, a.reverse, a.reclaim)
    elif a.part == 'step1_dry':
        part_step1_dry(a.wiring)
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
