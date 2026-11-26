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
    E136's row has its own 'arm'), the form, the KC-KC factors, the rung and the normalization's terms."""
    a = ARMS[arm]
    return dict(variant=variant(arm), form_arm=arm, form=a.form, kckc_factor=float(a.kckc),
                kckc_off_factor=float(a.kckc if a.kckc_off is None else a.kckc_off), rung=RUNG, **pn_terms())


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
    """A lane of reads: every planned point of the numbered steps of `steps` (or of `arms`; `planned_points`), the list
    reversed with `reverse`, skipping what is done or claimed elsewhere. The beside step is not written yet."""
    if 'beside' in steps:
        raise NotImplementedError("the beside step's reads are not written yet")
    pts = planned_points([s for s in steps if s in ('1', '2')], arms)
    for pt in (pts[::-1] if reverse else pts):
        read_point(pt, lane, reclaim, set(reads_table().tag))


# ---------------------------------------------------------------------------
# step 0, the beside step, the score, step 3 and the archive's table (stubs)

def step0_soma_splits() -> pd.DataFrame:
    """Step 0: E139's saved reads of F, B, K and A1 (and the family) at each panel's somata, the strength error split into
    level, spread and order (E136's `strength_split`), on all flies and on halves, with each odor's entry."""
    raise NotImplementedError


def step0_holds() -> pd.DataFrame:
    """Step 0: what E139's holds reach per variant and wiring: the somatic thresholds, the KC-KC and APL->KC scales,
    the response rates by KC type and the odors per KC."""
    raise NotImplementedError


def step0_a1_by_kckc() -> pd.DataFrame:
    """Step 0: A1's saved evaluations binned by the KC-KC factor, each stage's error along it (a lead: the other factors
    differ), as E139's `a1_by_rung` bins them by the additive term."""
    raise NotImplementedError


def read_beside(arm: str, wiring: str, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read of the beside step: the arm held on the APL16 input as its read holds it, then run on E98's extended
    input by E102's rules (as E139's `read_beside`), in this experiment's folders."""
    raise NotImplementedError


def part_score() -> dict:
    """The score: E136's, on T against the held-out half, with the unrefitted family as the reference and A1, E139's B
    and K at 1000 beside; the megamat somata with and without aa's entry; the kiwi mixtures and both mixture panels'
    second-block readings beside; the price on the APL effect (E139's rule) and the entry's added classes ("stands in
    for the rise", "does more than A1", "a lower price"); classes.json and provenance_score.json."""
    raise NotImplementedError


def part_fit(arm: str, lane: str = 'fit', reclaim: bool = False) -> dict:
    """Step 3, only on a stand-in: A1's fit repeated with the arm's form (the additive term with the scale by E77's
    rule, APL->bouton and the form's own scale), on S against the fit half, with E136's budget, returns, reads and
    control rule."""
    raise NotImplementedError


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv in out_dir(), for the model-output archive's params.csv (as E139's `write_variants`)."""
    raise NotImplementedError


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
