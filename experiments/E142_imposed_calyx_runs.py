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


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--part', required=True, choices=('step1', 'step1_dry', 'variants'))
    ap.add_argument('--lane')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--wiring', default='real')
    a = ap.parse_args(argv)
    if a.part == 'step1':
        part_step1(a.lane, a.reverse, a.reclaim)
    elif a.part == 'step1_dry':
        part_step1_dry(a.wiring)
    else:
        write_variants()


if __name__ == '__main__':
    main()
