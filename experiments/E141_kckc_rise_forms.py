#!/usr/bin/env python3
"""E141, the rise of the inhibition among KCs beside the weaker antennal-lobe normalization, along its dose and in other
forms (experiment_priorities.md, E141; planned by a Fable subagent, claude-fable-5-2, spawned from viki 2026-11-26
04:48).

E136's A1 holds through both of its calyx factors, and E139 found that KC-KC at x3.24 is what keeps validation2 and the
control mixtures (E139's outcome (c)). This asks what that rise does at the held-out panels and what could carry it:
its dose, and in its place a pooled neuron, an inhibition present with APL on only, and APL's inhibition of the claws as
a gain. Nothing is fitted first.

Definitions (fixed before any run; the entry's):
- the base point: pn_inhadd RUNG with pn_inhsc by E77's rule (`inhsc_for`) and APL->bouton at x APLPN, E139's B at 1000;
  A1 is taken as this point with the rise, as E139 took it for its two-by-two;
- the rise: the hemibrain's KC-KC matrix at x RISE of each wiring's own scale, in both conditions;
- the arms (ARMS): step 1's D1.8 and D5.8 (the matrix at x RISE**0.5 and x RISE**1.5), P (E122's uniform neuron, every
  KC inhibiting every KC alike, in the matrix's place at the summed weight of the matrix at x RISE), N (the matrix at x
  RISE with APL on and at x1 with APL off) and G (APL's inhibition of the claws as a gain, E136's A2 form
  `apl_inh_divisive`, with KC-KC at x1); step 2's M (MB-C1's two cells with their counted synapses beside the matrix at
  x1, at 2.24 times the matrix's summed weight, E131's unit), T and S (E122's one neuron per KC type, and the matrix
  shuffled within type pairs, each at the summed weight of x RISE);
- a point is an arm read on one wiring and one panel's input exactly as E136 reads an arm and E139 a setting
  (`read_point`: E136's `read_at`, both rates held at the panel's operating point: megamat 0.0915, validation2 0.050,
  the mixtures 0.054, each with the wiring's own APL-off ratio, the APL16 input at the wiring's own fit's two counts),
  on all eight wirings; a variant is an arm at the base point on every wiring and panel ('<arm>_inhadd<rung>',
  `variant`).

Judgments the entry leaves open (the narrowest reading; each is also in the docstring that applies it):
- the doses are read through E136's A1 at coordinates (`arm_z`), as E139's B and K were, so that D1.8 and D5.8 differ
  from E139's B and K at 1000 in the KC-KC factor alone;
- the forms are E136 arms of their own (`arm_def`, inside E136's `extra_arms`) with the base point's normalization in
  their kwargs: a matrix in the hemibrain's place as E122 passes it (`_wKCKC`, E122's `load_matrix`, the same total
  weight as the hemibrain matrix as fit_mb_model normalizes it, so that the KC-KC scale x RISE gives the summed weight
  of the matrix at x RISE); N's APL-off condition through E136's 'kckc_off_factor' (E77's `hold_rates` `off_extra`:
  the threshold is held with APL off at the APL-off KC-KC scale); G through E136's A2 kwargs and claw reader;
- RISE and APLPN are the entry's rounded x3.24 and x0.155 (E139's K and B), not A1's fitted 3.2378 and 0.15504.

Design:
- step 0 (no simulation, `--part step0`): E139's saved reads of F, B, K and A1 at each panel's somata, the strength error
  split into level, spread and order on all flies and on halves with each odor's entry (E136's `strength_split`); what
  the holds reach (the somatic thresholds, the KC-KC and APL->KC scales, the response rates by KC type, the odors per
  KC); A1's saved evaluations binned by the KC-KC factor (a lead);
- step 1 (`--part read --step 1`): D1.8, D5.8, P, N and G on the megamat, validation2, mixture and APL16 inputs, eight
  wirings each (160 reads), N and P first;
- step 2 (`--part read --step 2`): M if P stands in, else T and S (32 or 64 reads);
- beside (`--part read --step beside`): D1.8, D5.8 and every arm that holds on E98's extended input, read by E102's rules;
- the score (`--part score`): E136's, with the unrefitted family as the reference and A1, B and K beside (E139's reads),
  the megamat somata with and without aa's entry, the kiwi mixtures and both mixture panels' second-block readings
  beside, and the entry's added classes ("stands in for the rise", "does more than A1", "a lower price");
- step 3 (`--part fit`, only on a stand-in): A1's fit repeated with that form.

Outputs: results/E141_kckc_rise_forms/ (reads_<host>_<lane>.csv, provenance_read_<host>_<lane>.json, claims/, step0_*,
score_*.csv, classes.json, provenance_score.json, variants.csv); runs in outputs/E141_kckc_rise_forms/read/ and
data/model_outputs/E141_kckc_rise_forms/read/ (E134's step-3 layout, written in place, as E136's and E139's).

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E141_kckc_rise_forms --part read --step 1 --lane <host>[-<n>] [--reverse]
    ./venv/bin/python -u -m experiments.E141_kckc_rise_forms --part step0 | --part score | --part variants
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

EXPERIMENT = 'E141_kckc_rise_forms'
SCRIPT = 'experiments/E141_kckc_rise_forms.py'
# the base point (E139's B at 1000) and the rise (E139's K)
RUNG = 1000.0
APLPN = 0.155
RISE = 3.24
A1 = 'A1'
PANELS: Tuple[str, ...] = ('megamat', 'validation2', 'mixtures', 'apl16')
# "a lower price": the calibrated pattern part below A1's by more than this against all recordings and each half
LOWER_PRICE = 0.003


@dataclass(frozen=True)
class ArmDef:
    """One of the entry's arms at the base point.

    Attributes:
        name: the arm's short name ('N', 'P', 'D1.8', ...).
        label: what it is, in words.
        kckc: the KC-KC factor with APL on, of each wiring's own scale (of the form's matrix where `matrix` is set).
        form: 'dose' (the hemibrain matrix, read through E136's A1), 'matrix' (E122's matrix `matrix` in its place),
            'on_only' (the hemibrain matrix at `kckc` with APL on and `kckc_off` with APL off), 'gain' (E136's A2 form)
            or 'mbc1' (E131's MB-C1 cells beside the matrix).
        step: '1' or '2'.
        kckc_off: the KC-KC factor with APL off ('on_only'; None: the same as with APL on).
        matrix: E122's arm whose matrix is used ('U', 'TY', 'SH'; 'matrix' only).
    """
    name: str
    label: str
    kckc: float
    form: str
    step: str
    kckc_off: Optional[float] = None
    matrix: Optional[str] = None


ARMS: Dict[str, ArmDef] = {a.name: a for a in (
    ArmDef('N', f'the rise with APL on only: the matrix at x{RISE:g} with APL on and x1 with APL off', RISE, 'on_only',
           '1', kckc_off=1.0),
    ArmDef('P', f"E122's uniform neuron in the matrix's place, at the summed weight of the matrix at x{RISE:g}", RISE,
           'matrix', '1', matrix='U'),
    ArmDef('D1.8', f'the matrix at x{RISE ** 0.5:.3g} (the square root of {RISE:g})', RISE ** 0.5, 'dose', '1'),
    ArmDef('D5.8', f'the matrix at x{RISE ** 1.5:.3g} (the 1.5th power of {RISE:g})', RISE ** 1.5, 'dose', '1'),
    ArmDef('G', "APL's inhibition of the claws as a gain (E136's A2 form), KC-KC at x1", 1.0, 'gain', '1'),
    ArmDef('M', "MB-C1's two cells with their counted synapses beside the matrix at x1, at 2.24 times the matrix's "
                "summed weight (E131's unit)", 1.0, 'mbc1', '2'),
    ArmDef('T', f"E122's one neuron per KC type, at the summed weight of the matrix at x{RISE:g}", RISE, 'matrix', '2',
           matrix='TY'),
    ArmDef('S', f"the matrix shuffled within type pairs (E122's control), at the summed weight of the matrix at "
                f"x{RISE:g}", RISE, 'matrix', '2', matrix='SH'),
)}
# step 1's arms in a forward lane's order (the entry: N and P first); step 2's by step 1's outcome
STEP1: Tuple[str, ...] = ('N', 'P', 'D1.8', 'D5.8', 'G')
STEP2_IF_P: Tuple[str, ...] = ('M',)
STEP2_ELSE: Tuple[str, ...] = ('T', 'S')


@dataclass(frozen=True)
class Point:
    """An arm at the base point, read on one wiring and one panel's input."""
    arm: str
    panel: str
    wiring: str


# ---------------------------------------------------------------------------
# paths and small helpers

def out_dir() -> Path:
    """<shared>/results/E141_kckc_rise_forms/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E141_kckc_rise_forms/."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def local_dir(part: str = '') -> Path:
    """outputs/E141_kckc_rise_forms/[<part>/] (created)."""
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


# ---------------------------------------------------------------------------
# the arms, their variants and their models

def variant(arm: str) -> str:
    """An arm at the base point as a variant name, '<arm>_inhadd<RUNG>' (e.g. 'N_inhadd1000', 'D1.8_inhadd1000')."""
    if arm not in ARMS:
        raise ValueError(f'{arm=}')
    return f'{arm}_inhadd{RUNG:g}'


def parse_variant(name: str) -> str:
    """The arm of a variant name (`variant`)."""
    arm, sep, rung = str(name).partition('_inhadd')
    if not sep or arm not in ARMS or float(rung) != RUNG:
        raise ValueError(f'{name!r} is not a variant of {EXPERIMENT}')
    return arm


def point_tag(pt: Point) -> str:
    """'<variant>_<panel>_<wiring>', e.g. 'N_inhadd1000_megamat_P3'."""
    return f'{variant(pt.arm)}_{pt.panel}_{pt.wiring}'


def pn_terms() -> Dict[str, float]:
    """The base point's PN normalization terms (E136's `pn_norm_kws` at inhadd = ln(RUNG / olfsysm's default), inhsc 0:
    pn_inhsc by `inhsc_for`), as E139's rungs read them."""
    from experiments import E136_stage_population_fit as F
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    return F.pn_norm_kws(dict(inhadd=float(np.log(RUNG / PN_INHADD_DEFAULT)), inhsc=0.0))


def arm_z(arm: str) -> Dict[str, float]:
    """The arm's coordinates (log factors of each wiring's own scales): a dose's in E136's A1 (aplpn, kckc, and the
    normalization's inhadd = ln(RUNG / olfsysm's default), inhsc 0, as E139's `point_z`); a form's in its own E136 arm
    (aplpn and kckc; the normalization is in the arm's kwargs, `arm_def`)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    a = ARMS[arm]
    z = dict(aplpn=float(np.log(APLPN)), kckc=float(np.log(a.kckc)))
    if a.form == 'dose':
        z.update(inhadd=float(np.log(RUNG / PN_INHADD_DEFAULT)), inhsc=0.0)
    return z


def form_kws(arm: str) -> dict:
    """The model kwargs of a form beside the base point's normalization: 'matrix', dict(_wKCKC=<E122's matrix>);
    'on_only', dict(kckc_off_factor=kckc_off / kckc) (E136's APL-off condition of its own, no model kwarg); 'gain', E136's
    A2 kwargs; 'mbc1', not written yet (step 2, only if P stands in); a dose has none (it is read through A1)."""
    from experiments import E136_stage_population_fit as F
    a = ARMS[arm]
    if a.form == 'dose':
        return {}
    if a.form == 'matrix':
        from experiments import E122_feedback_neuron as E122
        return dict(_wKCKC=E122.load_matrix(a.matrix))
    if a.form == 'on_only':
        return dict(kckc_off_factor=float(a.kckc_off) / float(a.kckc))
    if a.form == 'gain':
        return dict(F.ARMS['A2'].kws)
    raise NotImplementedError(f"{arm}: the form {a.form!r} is not written yet (step 2's M runs only if P stands in)")


def arm_def(arm: str):
    """A form's E136 `Arm` (for E136's `extra_arms`), named as its variant: A0's two searched coordinates, the base
    point's normalization (`pn_terms`) and the form's kwargs (`form_kws`). ValueError for a dose (read through A1)."""
    from experiments import E136_stage_population_fit as F
    a = ARMS[arm]
    if a.form == 'dose':
        raise ValueError(f'{arm} is a dose, read through E136 A1')
    return F.Arm(variant(arm), a.label, tuple(F.ARMS['A0'].searched), dict(pn_terms(), **form_kws(arm)))


def planned_points(steps: Sequence[str] = ('1',), arms: Optional[Sequence[str]] = None) -> List[Point]:
    """The reads of the numbered steps (or of `arms`) in a forward lane's order: per arm (STEP1's order; step 2's by
    `step2_arms`), per panel, the eight wirings."""
    from experiments import E136_stage_population_fit as F
    names = list(arms) if arms is not None else [a for s in steps for a in (STEP1 if s == '1' else step2_arms())]
    return [Point(a, p, w) for a in names for p in PANELS for w in F.WIRINGS]


def step2_arms() -> Tuple[str, ...]:
    """Step 2's arms by step 1's score: M if P stands in (classes.json's P class 'holds'), else T and S; RuntimeError
    before the score of step 1."""
    p = out_dir() / 'classes.json'
    if not p.exists():
        raise RuntimeError("step 2's arms wait on step 1's score (classes.json)")
    cls = (json.loads(p.read_text()).get('classes') or {}).get(variant('P'), {}).get('cls')
    return STEP2_IF_P if cls == 'holds' else STEP2_ELSE


# ---------------------------------------------------------------------------
# the reads

def reads_table() -> pd.DataFrame:
    """Every lane's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def read_extra(arm: str) -> dict:
    """The fields a read records beside E136's (in its row and its run's info): the variant, the arm (as 'form_arm';
    E136's row has its own 'arm'), the form, the KC-KC factors, the rung, and a form's normalization terms (a dose is
    read through E136's A1, whose `read_at` records them itself: given twice, its row could not be built)."""
    a = ARMS[arm]
    return dict(variant=variant(arm), form_arm=arm, form=a.form, kckc_factor=float(a.kckc),
                kckc_off_factor=float(a.kckc if a.kckc_off is None else a.kckc_off), rung=RUNG,
                **({} if a.form == 'dose' else pn_terms()))


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read, claimed under out_dir()/claims/ as E136's and E139's reads are: E136's `read_at` of the arm's model at
    `arm_z` on the point's wiring and panel, both rates held at the panel's operating point (E136's `read_targets`),
    saved in read_dirs(); a dose through E136's A1 (RuntimeError unless E135's verdict leaves A1 in its 'search' mode),
    a form through its own arm (`arm_def`, inside E136's `extra_arms`). Its row (with `read_extra`'s fields) in
    reads_<host>_<lane>.csv. None if done or claimed elsewhere."""
    from experiments import E136_stage_population_fit as F
    if pt.arm not in ARMS or pt.panel not in F.READ_PANELS or pt.wiring not in F.WIRINGS:
        raise ValueError(f'{pt}')
    tag = point_tag(pt)
    if done is not None and tag in done:
        return None
    a = ARMS[pt.arm]
    kw = dict(reclaim=reclaim, done=done, extra=read_extra(pt.arm), label='E141', script=SCRIPT)
    if a.form == 'dose':
        mode = F.require_e135()['mode']
        if mode != 'search':
            raise RuntimeError(f"E135's verdict puts A1 in mode {mode!r}: its normalization is not the searched one")
        return F.read_at(tag, A1, arm_z(pt.arm), pt.wiring, pt.panel, None, lane, out_dir(), read_dirs(), **kw)
    with F.extra_arms(arm_def(pt.arm)):
        return F.read_at(tag, variant(pt.arm), arm_z(pt.arm), pt.wiring, pt.panel, None, lane, out_dir(), read_dirs(),
                         **kw)


def part_read(lane: str, steps: Sequence[str], arms: Optional[Sequence[str]] = None, reverse: bool = False,
              reclaim: bool = False) -> None:
    """A lane of reads: every planned point of the numbered steps of `steps` (or of `arms`; `planned_points`), then with
    'beside' the beside step's reads (`beside_points`: the doses and every arm that holds, or `arms`' with only 'beside'
    in `steps`), the whole list reversed with `reverse`, skipping what is done or claimed elsewhere."""
    nums = [s for s in steps if s in ('1', '2')]
    jobs = [('point', pt) for pt in planned_points(nums, arms)] if nums or (arms and 'beside' not in steps) else []
    if 'beside' in steps:
        jobs += [('beside', bp) for bp in beside_points(arms if not nums else None)]
    for kind, x in (jobs[::-1] if reverse else jobs):
        done = set(reads_table().tag)
        if kind == 'point':
            read_point(x, lane, reclaim, done)
        else:
            read_beside(x[0], x[1], lane, reclaim, done)


# ---------------------------------------------------------------------------
# step 0 (experiments/E141_step0.py), the beside step, the score, step 3 (a stub) and the archive's table

def step0_soma_splits() -> pd.DataFrame:
    """Step 0: E139's saved reads of F, B, K and A1 (and the family) at each panel's somata, the strength error split into
    level, spread and order (E136's `strength_split`), on all flies and on halves, with each odor's entry
    (experiments/E141_step0.py's `soma_splits`)."""
    from experiments import E141_step0 as Z
    return Z.soma_splits()


def step0_holds() -> pd.DataFrame:
    """Step 0: what E139's holds reach per variant and wiring: the somatic thresholds, the KC-KC and APL->KC scales,
    the response rates by KC type and the odors per KC (experiments/E141_step0.py's `holds`)."""
    from experiments import E141_step0 as Z
    return Z.holds()


def step0_a1_by_kckc() -> pd.DataFrame:
    """Step 0: A1's saved evaluations binned by the KC-KC factor, each stage's error along it (a lead: the other factors
    differ), as E139's `a1_by_rung` bins them by the additive term (experiments/E141_step0.py's `a1_by_kckc`)."""
    from experiments import E141_step0 as Z
    return Z.a1_by_kckc()


EXTENDED = 'extended'


def extended_root(name: str, where: str = 'local') -> Path:
    """The root of a variant's extended runs: outputs/E141_kckc_rise_forms/extended/<variant>/ ('local') or the
    archive's extended/<variant>/ ('archive')."""
    if where not in ('local', 'archive'):
        raise ValueError(f'{where=}')
    return (local_dir(EXTENDED) if where == 'local' else archive_dir() / EXTENDED) / name


def extended_dir(name: str, pn_inhadd: float, where: str = 'local') -> Path:
    """A variant's folder of extended runs, <root>/B/<E102's step tag of its pn_inhadd>/ (E102's part-B layout, as
    E139's `extended_dir`)."""
    from experiments import E102_antennal_lobe_normalization as E
    return extended_root(name, where) / 'B' / E.step_tag(float(pn_inhadd))


def same_as_read(name: str, wiring: str, hold: Mapping) -> Optional[dict]:
    """A beside hold against the variant's saved APL16 read (the same hold on the same input), if saved: whether the
    held parameters and the 16-odor counts are the same (as E139's `_same_as_read`); None if no read is saved."""
    from kcq import stage_readings as R
    try:
        run = R.load_step3_run(point_tag(Point(parse_variant(name), 'apl16', wiring)), dirs=read_dirs())
    except FileNotFoundError:
        return None
    if run is None:
        return None
    pa = {k: float(v) for k, v in (run.info.get('params') or {}).items()}
    pb = {k: float(v) for k, v in hold['params'].items()}
    same = all(np.array_equal(H.abbrev_columns(run.conditions[c]['spike_counts']).to_numpy(dtype=float),
                              H.abbrev_columns(hold[c][1]).to_numpy(dtype=float)) for c in ('on', 'off'))
    return dict(tag=run.tag, params_equal=pa == pb, counts_equal=bool(same))


def beside_context():
    """E139's `BesideContext` for this experiment: its folders, its arms as the rows' setting, its saved APL16 reads."""
    from experiments import E139_normalization_ladder as E139
    return E139.BesideContext(root=out_dir(), ext_dir=lambda n, a, where: extended_dir(n, a, where),
                              label_of=lambda n: (parse_variant(n), RUNG),
                              same_as_read=lambda n, w, h: same_as_read(n, w, h), script=SCRIPT, label='E141')


def beside_arms() -> Tuple[str, ...]:
    """The beside step's arms: the doses, and every arm that holds once step 1 or 2 is scored (classes.json)."""
    p = out_dir() / 'classes.json'
    cl = (json.loads(p.read_text()).get('classes') or {}) if p.exists() else {}
    return DOSES + tuple(a for a in ARMS if a not in DOSES and cl.get(variant(a), {}).get('cls') == 'holds')


def beside_points(arms: Optional[Sequence[str]] = None) -> List[Tuple[str, str]]:
    """The beside step's reads in a forward lane's order: (arm, wiring), each arm on the eight wirings."""
    from experiments import E136_stage_population_fit as F
    return [(a, w) for a in (beside_arms() if arms is None else arms) for w in F.WIRINGS]


def read_beside(arm: str, wiring: str, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read of the beside step, E139's `read_beside` in this experiment's folders (`beside_context`): the arm's model
    (a dose through E136's A1, a form through its own arm inside `extra_arms`) at `arm_z`, held on the lab's APL16 input
    as its read holds it, then run on E98's extended input by E102's rules; tag '<variant>_extended_<wiring>'."""
    from experiments import E136_stage_population_fit as F
    from experiments import E139_normalization_ladder as E139
    name, ctx = variant(arm), beside_context()
    if ARMS[arm].form == 'dose':
        return E139.read_beside(name, wiring, lane, reclaim, done, arm=A1, z=arm_z(arm), ctx=ctx)
    with F.extra_arms(arm_def(arm)):
        return E139.read_beside(name, wiring, lane, reclaim, done, arm=name, z=arm_z(arm), ctx=ctx)


# E139's variants scored beside (the base point without the rise, and the rise without APL->bouton's factor)
E139_BESIDE: Tuple[str, ...] = ('B_inhadd1000', 'K_inhadd1000')
# the forms whose holding is "stands in for the rise", and the doses whose holding is read against A1
FORMS: Tuple[str, ...] = tuple(a for a, d in ARMS.items() if d.form != 'dose')
DOSES: Tuple[str, ...] = tuple(a for a, d in ARMS.items() if d.form == 'dose')


def load_variant_run(pt, claw_after: bool = False):
    """The run of a variant's point for E136's score functions (E136's `ReadPoint` with the variant as its arm): E141's
    read (`read_dirs`); None if absent."""
    from kcq import stage_readings as R
    try:
        return R.load_step3_run(point_tag(Point(parse_variant(pt.arm), pt.panel, pt.wiring)), dirs=read_dirs())
    except FileNotFoundError:
        return None


def register_variants(F, names: Sequence[str]) -> None:
    """Registers each variant with E136's score functions (`register_variant`, loader `load_variant_run`)."""
    for v in names:
        parse_variant(v)
        F.register_variant(v, load_variant_run)


def label_of(v: str) -> dict:
    """The label columns of a name's score rows: the arm (or E139's setting, 'family', 'A1') and the rung."""
    try:
        return dict(setting=parse_variant(v), rung=RUNG)
    except ValueError:
        pass
    from experiments import E139_normalization_ladder as E139
    return dict(setting=E139._setting_of(v), rung=E139._rung_of(v))


def read_variants() -> List[str]:
    """The variants with at least one read (reads_table), in ARMS' order."""
    have = set(reads_table().get('variant', pd.Series(dtype=object)).dropna())
    return [variant(a) for a in ARMS if variant(a) in have]


def stages_read(F) -> Tuple[str, ...]:
    """The stages the score reads, as E139's `stages_read` (E136's `stage_classes_on_measure`, written to this
    experiment's stage_classes.json); RuntimeError if a stage is not readable."""
    cl = F.stage_classes_on_measure(F.MEASURES[F.MEASURE])
    _write_json(out_dir() / 'stage_classes.json', dict(measure=F.MEASURE, classes=cl, time=_now(), host=_host(),
                                                       **F._commits()))
    bad = [s for s in F.LOSS_STAGES if cl[s]['cls'] in ('not readable', 'not scored')]
    if bad:
        raise RuntimeError(f'on the measure {F.MEASURE!r} the stages {bad} are not readable: {cl}')
    return tuple(F.LOSS_STAGES)


def against_a1(F, v: str, mega, pv, spreads, parts) -> dict:
    """The entry's comparisons of a variant with A1, on T against the held-out half (E136's `compare`, A1 in the
    family's place): the megamat stage error ('megamat': below A1's by more than the bar is its 'improves'); each held-out
    panel of the class ('panels': 'better' where it is below A1's by more than the family's spread, 'worse' where above
    by more); and the calibrated pattern part against each target ('price': per target the mean over T of the
    per-wiring differences; 'lower' where every target's is below -LOWER_PRICE)."""
    from experiments import E139_normalization_ladder as E139
    rec: dict = {}
    if v in mega and A1 in mega:
        rec['megamat'] = F.compare({w: x['stage'] for w, x in mega[v].items()},
                                   {w: x['stage'] for w, x in mega[A1].items()}, spreads['stage'])
    pn = {}
    for key in F.HOLDS_PANELS:
        k = '/'.join(key)
        c = F.compare(pv.get(v, {}).get(key, {}), pv.get(A1, {}).get(key, {}), spreads.get(k, np.nan))
        if c.get('complete'):
            c['better'] = bool(np.isfinite(c['spread']) and c['change'] < -c['spread'])
            c['worse'] = bool(c['worse_than_spread'])
        pn[k] = c
    rec['panels'] = pn
    a, b = parts.get(v, {}), parts.get(A1, {})
    if all(w in a and w in b for w in F.T_WIRINGS):
        by = {t: float(np.mean([a[w][f'{t}_pattern'] - b[w][f'{t}_pattern'] for w in F.T_WIRINGS]))
              for t in E139.PRICE_TARGETS}
        rec['price'] = dict(by_target=by, threshold=LOWER_PRICE, lower=bool(max(by.values()) < -LOWER_PRICE))
    else:
        rec['price'] = dict(by_target=None, threshold=LOWER_PRICE, lower=None)
    return rec


def entry_classes(v: str, rec: Mapping, vs: Mapping) -> dict:
    """The entry's added classes of a variant from its E136 class (`rec`) and its comparison with A1 (`vs`): "stands in
    for the rise" (a form that holds), "does more than A1" (a dose that holds and is below A1's megamat stage error by
    more than the bar, or better than A1 on a held-out panel by more than its spread with none worse), "a lower price"
    (the calibrated pattern part below A1's by more than LOWER_PRICE against all recordings and each half; any arm)."""
    holds = rec.get('cls') == 'holds'
    arm = parse_variant(v)
    pn = [c for c in vs.get('panels', {}).values() if c.get('complete')]
    more = holds and (bool(vs.get('megamat', {}).get('improves'))
                      or (any(c['better'] for c in pn) and not any(c['worse'] for c in pn)))
    return dict(stands_in=bool(holds) if arm in FORMS else None, does_more_than_a1=bool(more) if arm in DOSES else None,
                lower_price=(vs.get('price') or {}).get('lower'))


def part_score(names: Optional[Sequence[str]] = None) -> dict:
    """The score, as E139's (`E139.score_sets`, the same values for every name): E141's variants (default every one with
    a read) through E136's score functions as registered variants, the unrefitted family as the reference and A1 (its S
    wirings its fit's best evaluation, as E136 reads it), E139's B and K at 1000 beside (E139_BESIDE, its reads):
    - the files of E139's score in this experiment's folder: score_megamat.csv, score_panels.csv,
      score_panels_all_views.csv, score_panels_e137.csv, score_apl16.csv, score_without_aa.csv;
    - E136's classes (`classify`, on T against the held-out half, the spreads over the family's eight wirings, rule 3's
      guard on every view; no permuted control and no bound: nothing is fitted; A1 with its own fit and control), no
      class for a variant whose holds fail (E139's `holds_failed`) or that leaves a realism band;
    - the price on the APL effect against the family (E139's `price_of`) and the comparisons with A1 (`against_a1`);
    - the entry's added classes (`entry_classes`), with the number of arms read beside any arm that holds.
    classes.json and provenance_score.json."""
    from kcq import stage_readings as R
    from experiments import E139_normalization_ladder as E139
    t0 = time.time()
    with E139.e136() as F:
        gate = F.require_gate()
        stages = stages_read(F)
        names = list(names) if names is not None else read_variants()
        register_variants(F, names)
        E139.register_variants(F, E139_BESIDE)
        allv = [F.FAMILY, *names, *E139_BESIDE, A1]
        S = E139.score_sets(F, allv, stages, label_of=label_of)
        mega, pv, spreads, guard = S['mega'], S['pv'], S['spreads'], S['guard']
        parts, apl16_runs, realism = S['parts'], S['apl16_runs'], S['realism']
        reads = pd.concat([reads_table(), E139.reads_table()], ignore_index=True, sort=False)
        e136_classes = json.loads((F.out_dir() / 'classes.json').read_text()) if (F.out_dir() / 'classes.json').exists() \
            else {}
        classes = {}
        for v in allv[1:]:
            if v not in mega or F.FAMILY not in mega:
                classes[v] = dict(cls='reads missing')
                continue
            if v == A1:
                fit = json.loads((F.out_dir() / f'fit_{A1}.json').read_text())
                control_fall = e136_classes.get('control_fall')
            else:
                fit, control_fall = {}, None
            rec = F.classify(v, stages, mega, pv, spreads, realism, fit, control_fall, guard)
            failed = E139.holds_failed(v, reads, apl16_runs.get(v, {}), F.HOLD_TOL)
            if failed:
                rec['reasons'].append('a hold fails: ' + ', '.join(failed))
                rec['cls'] = 'no class'
            rec['holds_failed'] = failed
            rec['flagged_wirings'] = {w: f for w, f in realism.get(v, {}).items() if f}
            rec['price_on_apl_effect'] = E139.price_of(parts, v, F.FAMILY)
            rec.update(label_of(v))
            if v in names:
                rec['against_a1'] = against_a1(F, v, mega, pv, spreads, parts)
                rec.update(entry_classes(v, rec, rec['against_a1']))
            classes[v] = rec
        holding = [v for v in names if classes.get(v, {}).get('cls') == 'holds']
        out = dict(stages=list(stages), measure=F.MEASURE, variants=names, beside_e139=list(E139_BESIDE),
                   reference=F.FAMILY, beside_arm=A1, base_point=dict(rung=RUNG, aplpn=APLPN, rise=RISE, **pn_terms()),
                   spreads=spreads, price_rule=dict(threshold=E139.PRICE, targets=list(E139.PRICE_TARGETS),
                                                    decided_on='T', lower_than_a1=LOWER_PRICE),
                   arms_read=len(names), arms_holding=holding, beside=S['beside'], classes=classes)
        pd.DataFrame(S['mrows']).to_csv(out_dir() / 'score_megamat.csv', index=False)
        S['panels'].to_csv(out_dir() / 'score_panels.csv', index=False)
        S['g_panels'].to_csv(out_dir() / 'score_panels_all_views.csv', index=False)
        S['e137'].to_csv(out_dir() / 'score_panels_e137.csv', index=False)
        S['apl'].to_csv(out_dir() / 'score_apl16.csv', index=False)
        if len(S['wo']):
            S['wo'].to_csv(out_dir() / f'score_without_{E139.WITHOUT_ODOR}.csv', index=False)
        # beside: the extended runs read by E102's rules (they gate nothing), E141's arms' beside E139's B at 1000 and
        # A1 and E102's own rows of its identity point and of the family's calyx at 1000
        own, own_w = E139.beside_readings([variant(a) for a in ARMS if extended_root(variant(a), 'archive').is_dir()
                                           or extended_root(variant(a), 'local').is_dir()],
                                          root_of=extended_root, source='E141:beside', e102_rows=False)
        ref, ref_w = E139.beside_readings(['B_inhadd1000', A1])
        ext = pd.concat([own, ref], ignore_index=True, sort=False)
        ext_w = pd.concat([own_w, ref_w], ignore_index=True, sort=False)
        ext.to_csv(out_dir() / 'score_beside_extended.csv', index=False)
        ext_w.to_csv(out_dir() / 'score_beside_extended_wirings.csv', index=False)
        out['beside']['extended'] = json.loads(ext.to_json(orient='records'))
        _write_json(out_dir() / 'classes.json', out)
        _write_json(out_dir() / 'provenance_score.json', dict(
            part='score', script=SCRIPT, lock=R.require_lock().get('sha256'), gate_measure=gate['gate'].get('measure'),
            seconds=round(time.time() - t0, 1), host=_host(), time=_now(), variants=names, **F._commits()))
    return {v: c.get('cls') for v, c in classes.items()}


def part_fit(arm: str, lane: str = 'fit', reclaim: bool = False) -> dict:
    """Step 3, only on a stand-in: A1's fit repeated with the arm's form (the additive term with the scale by E77's
    rule, APL->bouton and the form's own scale), on S against the fit half, with E136's budget, returns, reads and
    control rule."""
    raise NotImplementedError


def variant_kind(arm: str, reading: str, panel, wiring) -> str:
    """A run's description for the archive's params.csv (ASCII): the arm at the base point, then what the run is."""
    a = ARMS[arm]
    words = a.label.replace('→', '->').replace('×', 'x').replace('–', '-')
    what = f"arm {arm} ({words}) at pn_inhadd {RUNG:g} with APL->bouton x{APLPN:g}"
    if reading == EXTENDED:
        return (f"{what}, on E98's extended input (the beside step), wiring {wiring}, held at the wiring's own fit's "
                "APL16 counts")
    return f'{what}, read on the {panel} input, wiring {wiring}, both response rates held'


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv in out_dir(), for the model-output archive's params.csv (as E139's `write_variants`): one row per
    run the archive folder holds (`arch`, default archive_dir()): every read in read/ (<tag>_run.json) and the beside
    step's runs in extended/<variant>/B/<step>/ (<wiring>_point.json). A row has the tag (the prefix of the run's file
    names: the wiring, for a beside run), the folder (subdir, relative to `arch`; part, its top folder), the variant,
    the arm and its form, the panel ('extended' for a beside run) and wiring, the reading ('read' or 'extended'), a
    description (`variant_kind`), the coordinates (z_<name>), both conditions' absolute parameters (E136's
    VARIANT_PARAMS; 'off_<name>' with APL off: N's KC-KC scale differs there), the PN normalization's terms
    (`pn_terms`, which every arm has) and where and when it ran. (subdir, tag) is unique."""
    from experiments import E136_stage_population_fit as F
    arch = Path(arch) if arch is not None else archive_dir()
    runs = [(p, 'read') for p in sorted((arch / 'read').glob('*_run.json'))]
    runs += [(p, EXTENDED) for p in sorted(arch.glob(f'{EXTENDED}/*/B/*/*_point.json'))]
    rows = []
    for p, reading in runs:
        d = json.loads(p.read_text())
        name = d.get('variant')
        arm = d.get('form_arm') or parse_variant(name)
        panel = EXTENDED if reading == EXTENDED else d.get('panel')
        sub = p.parent.relative_to(arch).as_posix()
        r = dict(tag=p.name[:-len('_point.json')] if reading == EXTENDED else d['tag'], variant=name, arm=arm,
                 form=ARMS[arm].form, panel=panel, wiring=d.get('wiring'), reading=reading, model=d.get('model'),
                 subdir=sub, part=sub.split('/')[0], kind=variant_kind(arm, reading, panel, d.get('wiring')))
        r.update({f'z_{k}': v for k, v in (d.get('z') or {}).items()})
        r.update({k: (d.get('params') or {}).get(k) for k in F.VARIANT_PARAMS})
        r.update({f'off_{k}': (d.get('params_off') or {}).get(k) for k in F.VARIANT_PARAMS})
        r.update(pn_terms())
        r.update(host=d.get('host'), lane=d.get('lane'), time=d.get('time'))
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        raise FileNotFoundError(f'no <tag>_run.json or <wiring>_point.json under {arch}')
    if df.duplicated(['subdir', 'tag']).any():
        raise ValueError(f"runs named twice: {sorted(df.loc[df.duplicated(['subdir', 'tag']), 'tag'])}")
    df.to_csv(out_dir() / 'variants.csv', index=False)
    return df


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'read', 'score', 'fit', 'variants'))
    ap.add_argument('--step', nargs='+', default=['1'], choices=('1', '2', 'beside'))
    ap.add_argument('--arms', nargs='+', choices=tuple(ARMS), help='--part read: these arms in place of the steps\'')
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--arm', choices=tuple(ARMS), help='--part fit: the stand-in to fit')
    a = ap.parse_args(argv)
    if a.part == 'step0':
        for fn in (step0_soma_splits, step0_holds, step0_a1_by_kckc):
            df = fn()
            name = fn.__name__
            df.to_csv(out_dir() / f'{name}.csv', index=False)
            print(f'E141 step 0: {name}.csv ({len(df)} rows)', flush=True)
    elif a.part == 'read':
        from experiments import E136_stage_population_fit as F
        old = F.MEASURE
        F.set_measure('strength')
        try:
            F.require_gate()
        finally:
            F.set_measure(old)
        part_read(a.lane, a.step, a.arms, a.reverse, a.reclaim)
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'fit':
        if a.arm is None:
            raise SystemExit('--part fit needs --arm')
        print(json.dumps(part_fit(a.arm, lane=a.lane if a.lane != '1' else 'fit', reclaim=a.reclaim), indent=1,
                         default=str))
    elif a.part == 'variants':
        print(f'variants.csv: {len(write_variants())} runs', flush=True)


if __name__ == '__main__':
    main()
