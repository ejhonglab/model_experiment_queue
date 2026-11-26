#!/usr/bin/env python3
"""E139, the weaker antennal-lobe normalization without the calyx refit, along its strength: the megamat stages, the
held-out panels and the APL effect at four rungs of the PN normalization's additive term, and the two calyx factors of
E136's A1 one at a time (experiment_priorities.md, E139; planned by a Fable subagent, claude-fable-5-2, spawned from
viki 2026-11-25 18:39).

E136's A1 holds (the megamat stage error on the T wirings against the held-out half falls from the unrefitted family's
1.049 to 0.934) at a point that joins three changes: the additive term at 32 times the family's (1017 against 31.4),
APL->bouton at x0.155 and KC-KC at x3.24. This reads them apart at fixed scales; nothing is fitted first.

Definitions (fixed before any run):
- a rung is a value of pn_inhadd with pn_inhsc following by E77's rule (`inhsc_for`), as in E102's ladder and A1's
  coordinate (E136's `pn_norm_kws` at z_inhsc = 0): RUNGS;
- the settings (SETTINGS): F, the family's own calyx factors; B, APL->bouton at x0.155 with KC-KC at x1; K, KC-KC at
  x3.24 with APL->bouton at x1 (E136's A1 has both);
- a point is a setting at a rung, read on one wiring and one panel's input exactly as E136 reads an arm (`read_point`:
  E136's `read_at` of A1's model at `point_z`, both rates held at the panel's operating point: megamat 0.0915,
  validation2 0.050, the mixtures 0.054, each with the wiring's own APL-off ratio, the APL16 input at the wiring's own
  fit's two counts), on all eight wirings; a variant is a setting at a rung on every wiring and panel ('<setting>_inhadd
  <rung>', `variant`), scored through E136's score functions as a registered variant (`register_variants`,
  `load_variant_run`);
- a setting is read on megamat against the fit half, the held-out half and all flies on all eight wirings; a rung is
  named on S against the fit half; every class is E136's, on T against the held-out half, against the unrefitted
  family; "a price on the APL effect": the calibrated pattern part above the family's by more than PRICE (E102's rule)
  against all recordings or either of halves A and B; "serves both": holds with no price.

Judgments the entry leaves open (the narrowest reading; each is also in the docstring that applies it):
- the measure is E136's adopted one, strength by odor (MEASURE), with E136's guard on every view beside;
- F's APL16 rows are E102's saved runs, as the entry says. E102 recomputed the breadth behind the APL->KC weights at each
  rung (its `point_breadth`), where E136's reads (and so E139's own reads, B and K on the APL16 input among them) hold
  it at the lab input's; F's APL16 rows and B's, K's and A1's are therefore not the same model in that one respect;
- the price is read per wiring (E67's calibrated readout, K set on each run's APL-on counts, as E136's APL16 rows and
  E102's readings read it) against E83's targets z8 (all recordings), halfA and halfB, and decided on the mean over the T
  wirings of each wiring's change from the family's (the classes' wirings); the mean over all eight wirings (the entry's
  quoted 0.056 and 0.072) is given beside;
- "whose holds fail": a read whose hold ended at a bracket's bound or whose rate is more than E136's HOLD_TOL from its
  target, or an E102 run at a bound or off its counts, on any wiring of the variant;
- step 3's candidate: a rung below CANDIDATE_BELOW at F or B whose class improves and whose price (the largest of the
  three targets' mean changes on T) is under CANDIDATE_PRICE_SHARE of A1's, read the same way.

Design:
- step 0 (no simulation, `--part step0`): A1's saved evaluations and their rescore binned by rung (`a1_by_rung`); E102's
  saved ladder scored as E136 scores a read (`e102_ladder`); the PN stage alone along the ladder on E135's port
  (`pn_stage_ladder`); the recorded step's signature, each odor's residual off a proportional step against its summed
  ORN input, beside each rung's (`recorded_step_signature`);
- step 1 (`--part read --step 1`): F at the four rungs on the megamat, validation2 and mixture inputs, eight wirings
  each (96 reads); F's APL16 rows are E102's saved runs;
- step 2 (`--part read --step 2`): B at 314 and 1000 and K at 1000 on the megamat, validation2, mixture and APL16
  inputs, eight wirings each (96 reads);
- beside (`--step beside`): B at 1000 and A1 on E98's extended input, read by E102's rules for the flies' APL along
  concentration and the KCs at low concentration (16 reads): not written yet (`part_read` refuses it);
- the score (`--part score`): E136's, with the unrefitted family as the reference, E136's A1 beside, the kiwi mixtures
  and both mixture panels' second-block readings (E137) beside every setting (`part_score`);
- step 3 (`--part fit`, only on a candidate): if a rung below 1000 "improves" at F or B with a price under half of
  A1's, A0's two-factor fit with the normalization fixed at that rung, on S against the fit half, then E136's reads
  (`part_fit`: not written yet).

Outputs: results/E139_normalization_ladder/ (reads_<host>_<lane>.csv, provenance_read_<host>_<lane>.json, claims/,
step0_*.csv and step0_*.json, stage_classes.json, score_*.csv, classes.json, provenance_score.json, variants.csv from
`--part variants` for the archive's params.csv); runs in
outputs/E139_normalization_ladder/read/ and data/model_outputs/E139_normalization_ladder/read/ (E134's step-3 layout,
written in place, as E136's).

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 PYTHONUNBUFFERED=1 nice -n 10 ./venv/bin/python -u \\
        -m experiments.E139_normalization_ladder --part read --step 1 --lane <host>[-<n>] [--reverse]
    ./venv/bin/python -u -m experiments.E139_normalization_ladder --part step0 | --part score
"""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E139_normalization_ladder'
SCRIPT = 'experiments/E139_normalization_ladder.py'
E102 = 'E102_antennal_lobe_normalization'
# the measure E139 scores on: E136's adopted one (E134's gate of 2026-11-24 15:28: strength by odor alone, with E136's
# guard on every view beside)
MEASURE = 'strength'
# the rungs of pn_inhadd (olfsysm's default 31.4088 is the family's)
RUNGS = (100.0, 314.0, 1000.0, 3141.0)
# the calyx factors of the settings, as multiples of the family's (E136's A1: APL->bouton x0.155, KC-KC x3.24)
SETTINGS: Dict[str, Dict[str, float]] = {'F': {}, 'B': {'aplpn': 0.155}, 'K': {'kckc': 3.24}}
# per step: (setting, rungs) pairs and the panels read (F's APL16 rows are E102's saved runs)
STEPS: Dict[str, tuple] = {
    '1': ((('F', RUNGS),), ('megamat', 'validation2', 'mixtures')),
    '2': ((('B', (314.0, 1000.0)), ('K', (1000.0,))), ('megamat', 'validation2', 'mixtures', 'apl16')),
}
# "a price on the APL effect": the calibrated pattern part above the family's by more than this (E102's rule), against
# any of these targets (E83's: all recordings, halves A and B of E77's split)
PRICE = 0.003
PRICE_TARGETS: Tuple[str, ...] = ('z8', 'halfA', 'halfB')
# step 3's candidate: a rung below this "improves" at F or B with a price under this share of A1's
CANDIDATE_BELOW, CANDIDATE_PRICE_SHARE = 1000.0, 0.5
# beside every setting, the megamat somata without this odor's entry (E140's outcome (d), E138's definition)
WITHOUT_ODOR = 'aa'
# E136's arm whose model the points are read with (its normalization searched: E135's outcome (e), mode 'search'), and
# the arm read beside as the two-by-two's fourth cell
A1 = 'A1'
# step 0's PN-stage ladder: the family's point, A1's (E136's fit, its gain on pn_inhsc included) and the rungs
PN_POINTS: Tuple[str, ...] = ('family', A1) + tuple(f'inhadd{r:g}' for r in RUNGS)
# the PN stage's statistics read on the megamat odors E90's totals use (the 12 both panels have) and the bootstrap of
# the recorded step's signature
N_BOOT_SIGNATURE, SIGNATURE_SEED, SIGNATURE_CI = 1000, 0, 0.9
# the flies' spread of the odors' totals at the PN dendrites, as the entry quotes it (E90, E120; 12 megamat odors)
FLIES_TOTALS_CV = 0.312


@dataclass(frozen=True)
class Point:
    """A setting at a rung, read on one wiring and one panel's input."""
    setting: str
    rung: float
    panel: str
    wiring: str


def out_dir() -> Path:
    """<shared>/results/E139_normalization_ladder/ (created)."""
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def archive_dir() -> Path:
    """<shared>/data/model_outputs/E139_normalization_ladder/."""
    return H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT


def local_dir(part: str = '') -> Path:
    """outputs/E139_normalization_ladder/[<part>/] (created)."""
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


@contextlib.contextmanager
def e136() -> Iterator[object]:
    """E136's module with its measure set to MEASURE for the block (restored after)."""
    from experiments import E136_stage_population_fit as F
    old = F.MEASURE
    F.set_measure(MEASURE)
    try:
        yield F
    finally:
        F.set_measure(old)


def variant(setting: str, rung: float) -> str:
    """A setting at a rung as a variant name, '<setting>_inhadd<rung>' (e.g. 'F_inhadd314')."""
    if setting not in SETTINGS:
        raise ValueError(f'{setting=}')
    return f'{setting}_inhadd{float(rung):g}'


def parse_variant(name: str) -> Tuple[str, float]:
    """(setting, rung) of a variant name (`variant`)."""
    setting, sep, rung = str(name).partition('_inhadd')
    if not sep or setting not in SETTINGS:
        raise ValueError(f'{name!r} is not a variant')
    return setting, float(rung)


def variants(steps: Sequence[str] = ('1', '2')) -> List[str]:
    """The variants of `steps`, in their planned order."""
    return list(dict.fromkeys(variant(s, r) for st in steps for s, rungs in STEPS[st][0] for r in rungs))


def point_tag(pt: Point) -> str:
    """'<setting>_inhadd<rung>_<panel>_<wiring>', e.g. 'F_inhadd314_megamat_P3'."""
    if pt.setting not in SETTINGS:
        raise ValueError(f'{pt}')
    return f'{variant(pt.setting, pt.rung)}_{pt.panel}_{pt.wiring}'


def point_z(setting: str, rung: float) -> Dict[str, float]:
    """The point's coordinates in E136's A1 (log factors of each wiring's own scales; `pn_norm_kws`): aplpn and kckc
    from SETTINGS, inhadd = ln(rung / olfsysm's default), inhsc = 0 (pn_inhsc follows by `inhsc_for`)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    z = {k: float(np.log(v)) for k, v in SETTINGS[setting].items()}
    z.update(aplpn=z.get('aplpn', 0.0), kckc=z.get('kckc', 0.0), inhadd=float(np.log(rung / PN_INHADD_DEFAULT)),
             inhsc=0.0)
    return z


def planned_points(steps: Sequence[str] = ('1', '2')) -> List[Point]:
    """The reads of steps 1 and 2 in the order a forward lane takes them: per step, per setting and rung, per panel,
    the eight wirings."""
    from experiments import E136_stage_population_fit as F
    pts = []
    for s in steps:
        groups, panels = STEPS[s]
        for setting, rungs in groups:
            for r in rungs:
                for p in panels:
                    pts += [Point(setting, float(r), p, w) for w in F.WIRINGS]
    return pts


# ---------------------------------------------------------------------------
# the reads

def reads_table() -> pd.DataFrame:
    """Every lane's reads_<host>_<lane>.csv, the last row per tag."""
    fs = [pd.read_csv(f) for f in sorted(out_dir().glob('reads_*.csv'))]
    if not fs:
        return pd.DataFrame(columns=['tag'])
    return pd.concat(fs, ignore_index=True, sort=False).drop_duplicates('tag', keep='last').reset_index(drop=True)


def read_point(pt: Point, lane: str, reclaim: bool = False, done: Optional[set] = None) -> Optional[dict]:
    """One read, claimed under out_dir()/claims/ as E136's reads are: E136's read (`read_at`) of A1's model at
    `point_z(pt.setting, pt.rung)` on the point's wiring and panel, both rates held at the panel's operating point (E136's
    `read_targets`), saved in read_dirs() (the local outputs and the archive's read/); its row (with the variant, the
    setting and the rung) in reads_<host>_<lane>.csv. None if done or claimed elsewhere. RuntimeError unless E135's
    verdict leaves A1 in its 'search' mode (the normalization's two terms, `pn_norm_kws`)."""
    from experiments import E136_stage_population_fit as F
    if pt.panel not in F.READ_PANELS or pt.wiring not in F.WIRINGS or not pt.rung > 0:
        raise ValueError(f'{pt}')
    tag = point_tag(pt)
    if done is not None and tag in done:
        return None
    mode = F.require_e135()['mode']
    if mode != 'search':
        raise RuntimeError(f"E135's verdict puts A1 in mode {mode!r}: its normalization is not the searched one")
    return F.read_at(tag, A1, point_z(pt.setting, pt.rung), pt.wiring, pt.panel, None, lane, out_dir(), read_dirs(),
                     reclaim=reclaim, done=done,
                     extra=dict(variant=variant(pt.setting, pt.rung), setting=pt.setting, rung=float(pt.rung)),
                     label='E139', script=SCRIPT)


def part_read(lane: str, steps: Sequence[str], reverse: bool = False, reclaim: bool = False) -> None:
    """A lane of reads: every planned point of `steps` (`planned_points`, reversed with `reverse`), skipping what is
    done or claimed elsewhere. The beside step (E98's extended input) is not written yet: NotImplementedError."""
    if 'beside' in steps:
        raise NotImplementedError(
            "E139's beside step (B at 1000 and A1 on E98's extended input, 16 reads, read by E102's rules for the flies' "
            "APL along concentration and the KCs at low concentration) is not written: it needs E102's `run_extended` "
            "at a read's held parameters with E136's breadth, and E102's step readings on those runs")
    pts = planned_points(list(steps))
    for pt in (pts[::-1] if reverse else pts):
        read_point(pt, lane, reclaim, set(reads_table().tag))


# ---------------------------------------------------------------------------
# the variants' runs for E136's score: E139's reads, and for F on the APL16 input E102's saved ladder

def e102_step(rung: Optional[float]) -> str:
    """E102's part-B folder of a rung ('identity' for the family's own normalization; E102's `step_tag`)."""
    from experiments import E102_antennal_lobe_normalization as E
    return E.step_tag(E.STEP_IDENTITY if rung is None else float(rung))


def load_e102_run(rung: Optional[float], wiring: str) -> Optional['object']:
    """E102's saved part-B point at a rung (None: its identity point, the family's) on one wiring, on the lab's APL16
    input, as a somata-only run of kcq.stage_readings (`SavedRun`, source 'E102:B'): its 16-odor held outputs (the rates
    re-held at the wiring's fit's counts; '<wiring>_held_<on|acute_off>_<responses|spike_counts>.parquet' in
    outputs/E102_antennal_lobe_normalization/B/<step>/ on the host that ran it, 'B_<step>_<wiring>_held_...' in the
    archive), its point record (<wiring>_point.json) as the run's info. None if not saved."""
    from kcq import stage_readings as R
    step = e102_step(rung)
    cands = [(H.HERE / 'outputs' / E102 / 'B' / step, ''),
             (H.shared_dir() / 'data' / 'model_outputs' / E102 / 'B' / step, f'B_{step}_')]
    for d, prefix in cands:
        paths = {c: {k: d / f'{prefix}{wiring}_held_{R.CONDITION_FILE[c]}_{k}.parquet'
                     for k in ('responses', 'spike_counts')} for c in R.CONDITIONS}
        if not all(p.exists() for ps in paths.values() for p in ps.values()):
            continue
        conds = {c: {k: pd.read_parquet(p) for k, p in ps.items()} for c, ps in paths.items()}
        pj = d / f'{prefix}{wiring}_point.json'
        info = json.loads(pj.read_text()) if pj.exists() else {}
        files = [R.file_record(p) for ps in paths.values() for p in ps.values()] + ([R.file_record(pj)]
                                                                                   if pj.exists() else [])
        return R.SavedRun(source='E102:B', tag=f'B_{step}_{wiring}_held', model='C2', wiring=wiring,
                          input_panel='apl16', conditions=conds, orn=R._abbrev(H.apl16_orn_deltas()), files=files,
                          settled=True, info=dict(info, e102_step=step))
    return None


def load_variant_run(pt, claw_after: bool = False):
    """The run of a variant's point for E136's score functions (E136's `ReadPoint` with the variant as its arm): F on
    the APL16 input E102's saved ladder (`load_e102_run`), every other point E139's read (`read_dirs`); None if
    absent."""
    from kcq import stage_readings as R
    setting, rung = parse_variant(pt.arm)
    if setting == 'F' and pt.panel == 'apl16':
        return load_e102_run(rung, pt.wiring)
    try:
        return R.load_step3_run(point_tag(Point(setting, rung, pt.panel, pt.wiring)), dirs=read_dirs())
    except FileNotFoundError:
        return None


def register_variants(F, names: Sequence[str]) -> None:
    """Registers each variant with E136's score functions (`register_variant`, loader `load_variant_run`)."""
    for v in names:
        parse_variant(v)
        F.register_variant(v, load_variant_run)


def run_outputs(run) -> Dict[str, tuple]:
    """{'on', 'off'}: (responses, spike counts) of a saved run."""
    return {c: (run.conditions[c]['responses'], run.conditions[c]['spike_counts']) for c in ('on', 'off')}


# ---------------------------------------------------------------------------
# the price on the APL effect

def _canonical(s: pd.Series) -> pd.Series:
    from experiments import E83_cross_wiring_fit as E83
    return E83._cs(s)


def apl_effect_parts(o: Mapping[str, tuple], targets: Optional[Mapping[str, pd.Series]] = None) -> Dict[str, float]:
    """One run's APL effect under the calibrated readout (E83's `pair_distances(o, 'cal')`: the saturating readout, K
    set on the run's APL-on counts), against each target of PRICE_TARGETS (E83's `targets`): the loss, the pattern part,
    its amplitude and shape parts, the spread ratio, ρ and Pearson r on the per-pair APL effect ('<target>_<name>'),
    the mean APL effect and the number of pairs. A pair the run leaves undefined (a silent odor) is left out on both
    sides, each side's mean taken over the pairs kept (as E102's `mask_pairs`); with every pair defined the values are
    E83's `effect_scores` and H's `apl_effect_spread_shape` / `apl_effect_stats`."""
    from experiments import E83_cross_wiring_fit as E83
    tg = E83.targets() if targets is None else targets
    d = E83.pair_distances(dict(o), 'cal')
    me = _canonical(H.apl_effect(d['on'], d['off']))
    ok = np.isfinite(me.to_numpy(dtype=float))
    row: Dict[str, float] = {'mean_effect': float(me[ok].mean()), 'n_pairs': int(ok.sum())}
    keep = me.index[ok]
    on, off = _canonical(d['on']).reindex(keep), _canonical(d['off']).reindex(keep)
    for t in PRICE_TARGETS:
        de = _canonical(tg[t]).reindex(keep)
        m = me.reindex(keep)
        err = (m - m.mean()) - (de - de.mean())
        row[f'{t}_loss'] = float(np.sqrt(np.mean((m - de).to_numpy(dtype=float) ** 2)))
        row[f'{t}_pattern'] = float(np.sqrt(np.mean(err.to_numpy(dtype=float) ** 2)))
        row[f'{t}_rho'] = float(H.spearman(m, de))
        de_c = H.canonical_pairs(tg[t])
        sp = H.apl_effect_spread_shape(on, off, de_c)
        st = H.apl_effect_stats(on, off, de_c)
        row.update({f'{t}_amplitude': float(sp['apl_effect_amplitude']), f'{t}_shape': float(sp['apl_effect_shape']),
                    f'{t}_sd_ratio': float(sp['apl_effect_sd_ratio']), f'{t}_pearson': float(st['apl_effect_pearson'])})
    return row


def price_of(parts: Mapping[str, Mapping[str, Mapping[str, float]]], name: str, family: str) -> dict:
    """The price on the APL effect of one variant against the family (`parts`: per name, per wiring, `apl_effect_parts`):
    per target of PRICE_TARGETS the per-wiring changes of the calibrated pattern part, their mean over T (deciding) and
    over all eight wirings (beside); 'price' where a T mean is above PRICE; 'price_value' the largest T mean (step 3's
    candidate rule reads it). Incomplete where a wiring is missing."""
    from experiments import E136_stage_population_fit as F
    a, f = parts.get(name, {}), parts.get(family, {})
    rec = dict(complete=all(w in a and w in f for w in F.WIRINGS), threshold=PRICE, by_target={})
    if not rec['complete']:
        rec.update(price=None, price_value=None, missing=[w for w in F.WIRINGS if w not in a or w not in f])
        return rec
    for t in PRICE_TARGETS:
        ch = {w: a[w][f'{t}_pattern'] - f[w][f'{t}_pattern'] for w in F.WIRINGS}
        rec['by_target'][t] = dict(changes=ch, change_T=float(np.mean([ch[w] for w in F.T_WIRINGS])),
                                   change_8=float(np.mean(list(ch.values()))),
                                   value_T=float(np.mean([a[w][f'{t}_pattern'] for w in F.T_WIRINGS])),
                                   value_8=float(np.mean([a[w][f'{t}_pattern'] for w in F.WIRINGS])),
                                   family_T=float(np.mean([f[w][f'{t}_pattern'] for w in F.T_WIRINGS])),
                                   family_8=float(np.mean([f[w][f'{t}_pattern'] for w in F.WIRINGS])))
    vals = [v['change_T'] for v in rec['by_target'].values()]
    rec.update(price=bool(max(vals) > PRICE), price_value=float(max(vals)),
               price_8=bool(max(v['change_8'] for v in rec['by_target'].values()) > PRICE))
    return rec


# ---------------------------------------------------------------------------
# step 0

def _e136_results() -> Path:
    """E136's results on the adopted measure (read only)."""
    with e136() as F:
        return F.out_dir()


def rung_bins(inhadd: np.ndarray) -> np.ndarray:
    """Each pn_inhadd's bin: the nearest in log of the family's (olfsysm's default) and RUNGS (the edges at the
    geometric midpoints)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    centers = np.array((PN_INHADD_DEFAULT,) + RUNGS)
    x = np.log(np.asarray(inhadd, dtype=float))
    return centers[np.argmin(np.abs(x[:, None] - np.log(centers)[None, :]), axis=1)]


def a1_by_rung(results: Optional[Path] = None) -> pd.DataFrame:
    """Step 0: A1's saved evaluations (fit_history_A1.csv) and their rescore against both halves (rescore_A1.csv,
    strength by odor; the fit half's equals the history's), binned by the additive term (`rung_bins`): per bin and half,
    'best' (the bin's evaluation with the lowest megamat stage error against the fit half, the S wirings' mean views,
    read against that half) and 'median' (the median over the bin's evaluations of each column), with the number of
    evaluations, pn_inhadd, the other coordinates as factors, the megamat stage error, the error over the boutons, claws
    and somata alone (stage_error_downstream) and each stage's error. Single evaluations whose other factors differ: a
    lead."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT
    d = Path(results) if results is not None else _e136_results()
    hist = pd.read_csv(d / f'fit_history_{A1}.csv', float_precision='round_trip')
    resc = pd.read_csv(d / f'rescore_{A1}.csv', float_precision='round_trip')
    resc = resc[resc.measure == MEASURE]
    z = hist.set_index('eval')[[f'z_{k}' for k in ('aplpn', 'kckc', 'inhadd', 'inhsc')]]
    cols = ['stage_error', 'stage_error_downstream', 'pn_dendrite_error', 'pn_bouton_error', 'kc_claw_error',
            'kc_soma_error']
    by_half = {h: resc[resc.half == h].set_index('eval')[cols] for h in ('fit', 'held')}
    inhadd = PN_INHADD_DEFAULT * np.exp(z['z_inhadd'])
    bins = pd.Series(rung_bins(inhadd.to_numpy()), index=z.index)
    rows = []
    for b in sorted(set(bins)):
        ev = [i for i in bins.index[bins == b] if i in by_half['fit'].index]
        if not ev:
            continue
        best = int(by_half['fit'].loc[ev, 'stage_error'].idxmin())
        for half, tab in by_half.items():
            have = [i for i in ev if i in tab.index]
            base = dict(bin_rung=float(b), half=half, n_evals=len(have))
            for stat in ('best', 'median'):
                if stat == 'best':
                    if best not in tab.index:
                        continue
                    vals = tab.loc[best]
                    zz = z.loc[best]
                    row = dict(base, statistic='best', eval=best, pn_inhadd=float(inhadd.loc[best]))
                else:
                    vals = tab.loc[have].median()
                    zz = z.loc[have].median()
                    row = dict(base, statistic='median', eval=np.nan, pn_inhadd=float(inhadd.loc[have].median()))
                row.update(aplpn_factor=float(np.exp(zz['z_aplpn'])), kckc_factor=float(np.exp(zz['z_kckc'])),
                           inhsc_gain=float(np.exp(zz['z_inhsc'])))
                row.update({c.replace('_error', '') if c.endswith('_error') and c != 'stage_error' else c: float(vals[c])
                            for c in cols})
                rows.append(row)
    return pd.DataFrame(rows)


def _e102_readings_by_step() -> pd.DataFrame:
    """E102's saved per-step readings (b_steps.csv): the flies' APL along concentration (E98's V2 in its readings, V1,
    E100's W1), the KCs at ma 10^-7 (the median over wirings, the flies', the class), silent odors and flags."""
    p = H.results_dir(E102) / 'b_steps.csv'
    b = pd.read_csv(p, comment='#')
    keep = ['step', 'value', 'v2', 'v2_held', 'v2_linear', 'v2_linear_min', 'v2_linear_max', 'v1', 'w1', 'n_flagged',
            'silent', 'silent_odors', 'low_on', 'low_off', 'low_off_over_on', 'low_flies_on', 'low_flies_off',
            'low_class_on', 'low_class_off']
    return b[[c for c in keep if c in b.columns]].rename(columns={c: f'e102_{c}' for c in keep if c != 'step'})


def _e102_low_conc_by_wiring() -> pd.DataFrame:
    """E102's saved per-wiring KC response rates at its carrying low-concentration stimulus (ma 10^-7 of the
    kiwi-control session, b_low_conc.csv)."""
    from experiments import E102_antennal_lobe_normalization as E
    lc = pd.read_csv(H.results_dir(E102) / 'b_low_conc.csv', comment='#')
    lc = lc[lc.stimulus == E.LOW_CARRYING]
    return lc[['setting', 'wiring', 'rate_on', 'rate_off', 'flies_on', 'flies_off']].rename(
        columns={'setting': 'step', 'rate_on': 'low_rate_on', 'rate_off': 'low_rate_off', 'flies_on': 'low_flies_on',
                 'flies_off': 'low_flies_off'})


def _panel_rows(F, name: str, views) -> pd.DataFrame:
    p = F.score_panels(name, views=views)
    return p[(p.input == 'apl16') & (p.panel == 'apl16')] if len(p) else p


def e102_ladder() -> pd.DataFrame:
    """Step 0: E102's saved ladder (its identity point, the family, and F at RUNGS on the APL16 input; `load_e102_run`)
    scored as E136 scores a read, wiring by wiring with the T and eight-wiring means: the APL16 somata's population
    error on strength by odor (E136's measure) and on every view, with APL on and off (`score_panels`; the model's
    APL-off population against APL output blocked (TNT)); E136's APL16 row (E67's `full_row`: the APL-effect statistics
    under the raw and the calibrated readouts, the realism bands, the largest counts; `apl16_rows`); the price's
    statistics against E83's z8, halfA and halfB (`apl_effect_parts`) with each wiring's change from the family; and,
    from E102's saved score, the flies' APL readings along concentration and the KCs at ma 10^-7, rung by rung
    (`_e102_readings_by_step`, b_steps.csv) and wiring by wiring (b_low_conc.csv). These are F's APL16 rows of the
    score. E102 recomputed the breadth behind the APL->KC weights at each rung (its `point_breadth`), where E136's and
    E139's reads hold it at the lab input's."""
    with e136() as F:
        names = [F.FAMILY] + [variant('F', r) for r in RUNGS]
        register_variants(F, names[1:])
        rows = []
        for name in names:
            rung = None if name == F.FAMILY else parse_variant(name)[1]
            apl = F.apl16_rows(name)
            apl = apl.set_index('wiring') if len(apl) else pd.DataFrame()
            pans = {v: _panel_rows(F, name, v) for v in ('measure', None)}
            for w in F.WIRINGS:
                run = F.load_read(F.ReadPoint(name, 'apl16', w, None))
                if run is None:
                    continue
                row = dict(variant=name, step=e102_step(rung), rung=rung if rung is not None else np.nan, wiring=w)
                for views, p in pans.items():
                    lab = MEASURE if views == 'measure' else 'all'
                    for c in ('on', 'off'):
                        x = p[(p.condition == c) & (p.wiring == w)] if len(p) else p
                        if not len(x):
                            continue
                        r = x.iloc[0]
                        row[f'{lab}_error_{c}'] = float(r['error'])
                        row.update({f'{lab}_{k}_{c}': float(r[k]) for k in r.index if str(k).startswith('view_')})
                if w in apl.index:
                    row.update({f'apl16_{k}': v for k, v in apl.loc[w].items() if k != 'arm'})
                info = run.info or {}
                row.update(hold_at_bound=bool(info.get('at_bound_thr') or info.get('at_bound_apl')),
                           hold_off_counts=bool(info.get('n_on') != info.get('target_on')
                                                or info.get('n_off') != info.get('target_off')) if info else np.nan)
                row.update({f'price_{k}': v for k, v in apl_effect_parts(run_outputs(run)).items()})
                rows.append(row)
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    fam = df[df.variant == names[0]].set_index('wiring')
    for t in PRICE_TARGETS:
        col = f'price_{t}_pattern'
        df[f'd_{col}'] = df[col] - df.wiring.map(fam[col])
    for col in [c for c in df.columns if c.startswith((f'{MEASURE}_error', 'all_error'))]:
        df[f'd_{col}'] = df[col] - df.wiring.map(fam[col])
    df = df.merge(_e102_readings_by_step(), on='step', how='left').merge(_e102_low_conc_by_wiring(),
                                                                        on=['step', 'wiring'], how='left')
    num = [c for c in df.columns if c not in ('variant', 'step', 'wiring') and not c.startswith('e102_')
           and pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c])]
    means = []
    with e136() as F:
        for (name, step), g in df.groupby(['variant', 'step'], sort=False):
            for lab, ws in (('mean_T', F.T_WIRINGS), ('mean_8', F.WIRINGS)):
                x = g[g.wiring.isin(ws)]
                if len(x) == len(ws):
                    means.append(dict(variant=name, step=step, wiring=lab, **x[num].mean().to_dict(),
                                      **{c: g[c].iloc[0] for c in df.columns if c.startswith('e102_')}))
    return pd.concat([df, pd.DataFrame(means)], ignore_index=True, sort=False)


_STAGE: dict = {}


def family_pn_stage():
    """The family's PN stage on the megamat input (E135's `family_stage` on the real wiring, checked against olfsysm
    there; built once per process)."""
    from experiments import E135_antennal_lobe_dendrites as AL
    if 'stage' not in _STAGE:
        _STAGE['stage'], _ = AL.family_stage('megamat', 'real', purpose='E139 step 0')
    return _STAGE['stage']


def pn_stage_points() -> Dict[str, Dict[str, float]]:
    """The normalization of each point of the PN-stage ladder (PN_POINTS): the family's (olfsysm's defaults), A1's (E136's
    fitted pn_inhadd and pn_inhsc, its gain on pn_inhsc included) and each rung's (pn_inhsc by `inhsc_for`)."""
    from experiments.E3_drive_dependent_gain import PN_INHADD_DEFAULT, PN_INHSC_DEFAULT
    from experiments.E77_input_screen import inhsc_for
    fit = json.loads((_e136_results() / f'fit_{A1}.json').read_text())
    a1 = fit.get('pn_kws') or {}
    if not a1:
        raise RuntimeError(f"E136's fit_{A1}.json has no pn_kws")
    out = {'family': dict(pn_inhadd=float(PN_INHADD_DEFAULT), pn_inhsc=float(PN_INHSC_DEFAULT)),
           A1: dict(pn_inhadd=float(a1['pn_inhadd']), pn_inhsc=float(a1['pn_inhsc']))}
    for r in RUNGS:
        out[f'inhadd{r:g}'] = dict(pn_inhadd=float(r), pn_inhsc=float(inhsc_for(r)))
    return out


def _pn_port_run(rates: pd.DataFrame, wiring: str, tag: str):
    """The port's evoked PN rates (glomeruli x odors) as a PN-dendrite-only saved run of kcq.stage_readings on the
    megamat input (window = the evoked rates, pre-odor = 0, so the reading's change is the rates; an empty spike-count
    table names the odors), read as E136 reads a run's PN dendrites."""
    from kcq import stage_readings as R
    x = rates.astype(float).copy()
    x.index = x.index.astype(str)
    st = {'pn_win': x, 'pn_pre': x * 0.0, 'spike_counts': pd.DataFrame(np.zeros((1, x.shape[1])), columns=x.columns)}
    return R.SavedRun(source='E139:pn_port', tag=tag, model='C2', wiring=wiring, input_panel='megamat',
                      conditions={'on': st}, settled=True)


def totals_cv(rates: pd.DataFrame, glomeruli: Sequence[str], odors: Sequence[str] = tuple(H.MEGAMAT12)) -> float:
    """The spread of the odors' totals (E135's `totals_cv`, E119's and E120's statistic): the coefficient of variation
    over `odors` of each odor's mean evoked rate over `glomeruli` (SD with ddof 0 over the mean)."""
    gl = [g for g in glomeruli if g in set(rates.index.astype(str))]
    x = rates.copy()
    x.index = x.index.astype(str)
    tot = x.loc[gl, [o for o in odors if o in x.columns]].mean(axis=0).to_numpy(dtype=float)
    return float(np.std(tot) / np.mean(tot)) if len(tot) and np.mean(tot) != 0 else np.nan


def pn_stage_ladder(stage=None) -> pd.DataFrame:
    """Step 0: the PN stage alone along the ladder (E135's numpy port of olfsysm's ORN, LN and PN layers, `PNStage`, on
    the family's megamat input; `family_stage` on the real wiring, checked against olfsysm there; the PN stage does not
    depend on the wiring), at each point of `pn_stage_points`: its evoked rates read as E136 reads the PN dendrites
    (`_pn_port_run`, E136's StageLoss reading on each wiring's seed index, the curve and noise of each half's flies,
    every view, all flies' one-fly units), against the fit half, the held-out half and all flies: the mean over the
    eight wirings of each wiring's population error on every view (error_with_matrix) and on strength by odor, and of
    its strength, matrix, sharing and breadth views; the spread of the odors' totals over the 12 megamat odors on the 28
    shared glomeruli (`totals_cv`; the flies' 0.312, FLIES_TOTALS_CV). The family's point against the held-out half is
    checked against E136's score_megamat.csv (the family's PN-dendrite error on each wiring, strength by odor):
    'identity_max_abs_diff' in the table's attrs."""
    from experiments import E135_antennal_lobe_dendrites as AL
    from kcq import population_fit as P
    from kcq import stage_readings as R
    with e136() as F:
        stage = family_pn_stage() if stage is None else stage
        pts = pn_stage_points()
        gl28 = AL.shared_glomeruli()
        lock = R.require_lock()
        rows, ident = [], {}
        fam_rates = {}
        for label, kws in pts.items():
            rates = stage.run(kws)
            fam_rates[label] = rates
            runs = {w: _pn_port_run(rates, w, f'pn_{label}_{w}') for w in F.WIRINGS}
            for half in ('fit', 'held', None):
                for views in (None, F.MEASURES[MEASURE]):
                    flies = F.Flies(half, ('pn_dendrite',), views=views)
                    reads = F.StageLoss(('pn_dendrite',), flies, lock).read(runs)
                    pops = [reads[w].populations[('pn_dendrite', F.PANEL, 'on')] for w in F.WIRINGS]
                    odors = list(pops[0].values.columns)
                    tab = P.wirings_as_flies(pops, flies.reference('pn_dendrite', odors))
                    lab = 'all' if views is None else MEASURE
                    for which in ('mean', 'common'):
                        r = tab.loc[which]
                        row = dict(point=label, half=half or 'all', views=lab, row=which, **kws,
                                   error=float(r['error']), totals_cv=totals_cv(rates, gl28))
                        if views is None:
                            row['error_with_matrix'] = float(r['error_with_matrix'])
                        row.update({k: float(v) for k, v in r.items() if str(k).startswith('view_')})
                        rows.append(row)
                    if label == 'family' and half == 'held' and views is not None:
                        ident = {w: float(tab.loc[p.name, 'error']) for w, p in zip(F.WIRINGS, pops)}
    df = pd.DataFrame(rows)
    df['flies_totals_cv'] = FLIES_TOTALS_CV
    sc = pd.read_csv(_e136_results() / 'score_megamat.csv')
    fam = sc[(sc.arm == 'family') & (sc.half == 'held')].set_index('wiring')['pn_dendrite']
    diffs = {w: abs(v - float(fam[w])) for w, v in ident.items() if w in fam.index}
    df.attrs['summary'] = dict(identity_against='E136 score_megamat.csv, family, held-out half, pn_dendrite',
                               identity_max_abs_diff=max(diffs.values()) if diffs else None, identity_by_wiring=diffs,
                               points=pts, glomeruli_totals=list(gl28))
    return df


def _signature(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """A step's signature by odor: with each glomerulus's proportional step (y on x through the origin, E135's
    `_slopes`, over its odors), each odor's mean residual over the glomeruli as a fraction of the step's mean prediction
    over all glomeruli and odors (scale-free, so the recorded dendrites' units and the model's Hz compare); x, y
    glomeruli x odors."""
    from experiments import E135_antennal_lobe_dendrites as AL
    s, _, _ = AL._slopes(x, y)
    pred = s[:, None] * x
    res = y - pred
    ok = np.isfinite(res)
    with np.errstate(invalid='ignore'):
        scale = float(np.nanmean(np.abs(np.where(ok, pred, np.nan))))
        return np.array([np.nanmean(res[ok[:, j], j]) if ok[:, j].any() else np.nan for j in range(res.shape[1])]) / scale


def recorded_step_signature(stage=None, n_boot: int = N_BOOT_SIGNATURE, seed: int = SIGNATURE_SEED) -> pd.DataFrame:
    """Step 0: the recorded step's own signature beside each point's (`pn_stage_points`). On the 28 shared glomeruli
    and the 17 megamat odors, in E135's rungs' reading (the ORN terminals in dF/F, the dendrites in the lock's trial-noise
    units; RUNG_STEP_READING): each odor's residual off a proportional step from the recorded terminals to the recorded
    dendrites (fly means; `_signature`), with a 90 % interval from resampling the flies of both stages (E135's
    `centered`), against the odor's summed terminal response (over the model input's glomeruli the recordings hold,
    E135's `_total`); beside it each point's signature, the step from the model's own megamat input (its deltas, the
    terminals' dF/F times one factor) to its evoked PN rates (E135's port). The check first: the family's and A1's
    signatures can be told apart through the flies' noise where the rms over odors of their difference exceeds the rms
    of the flies' interval's half-widths ('told_apart' in the attrs' summary, with the odors where the difference exceeds
    the half-width). One row per odor; the summary gives, per point, ρ of its signature with the summed input and the rms
    difference from the flies', with the flies' ρ and its interval."""
    from scipy import stats
    from experiments import E135_antennal_lobe_dendrites as AL
    gl = AL.shared_glomeruli()
    odors = list(AL.MEGAMAT_ODORS)
    u_orn, u_pn = AL.STEP_READINGS[AL.RUNG_STEP_READING]
    orn_t = AL.recorded_stage(AL.ORN_KEY, u_orn)
    pn_t = AL.recorded_stage(AL.PN_KEY, u_pn)
    orn, pn = AL.fly_stack(orn_t, gl, odors), AL.fly_stack(pn_t, gl, odors)
    model_gl = list(H.apl16_orn_deltas().index)
    pool = [g for g in model_gl if g in set().union(*[set(t.index) for t in orn_t.values()])]
    total = AL._total(AL.fly_stack(orn_t, pool, odors))
    sig = _signature(AL._nanmean(orn), AL._nanmean(pn))
    rng = np.random.default_rng(seed)
    draws = np.full((n_boot, len(odors)), np.nan)
    rho_draws = np.full(n_boot, np.nan)
    for b in range(n_boot):
        io = rng.integers(0, orn.shape[0], orn.shape[0])
        ip = rng.integers(0, pn.shape[0], pn.shape[0])
        draws[b] = _signature(AL._nanmean(orn[io]), AL._nanmean(pn[ip]))
        ok = np.isfinite(draws[b])
        rho_draws[b] = stats.spearmanr(draws[b][ok], total[ok])[0] if ok.sum() > 2 else np.nan
    q = ((1 - SIGNATURE_CI) / 2, 1 - (1 - SIGNATURE_CI) / 2)
    lo, hi = AL.centered(sig, draws, q)
    half = (hi - lo) / 2.0
    stage = family_pn_stage() if stage is None else stage
    x_model = stage.delta.reindex(index=gl, columns=odors).to_numpy(dtype=float)
    df = pd.DataFrame(dict(odor=odors, total_orn=total, flies=sig, flies_lo=lo, flies_hi=hi))
    pts = pn_stage_points()
    summary = {}
    ok_f = np.isfinite(sig)
    rho_f = float(stats.spearmanr(sig[ok_f], total[ok_f])[0])
    rlo, rhi = AL.centered(np.array([rho_f]), rho_draws[:, None], q)
    summary['flies'] = dict(rho_with_total=rho_f, rho_lo=float(rlo[0]), rho_hi=float(rhi[0]))
    for label, kws in pts.items():
        y = stage.run(kws).reindex(index=gl, columns=odors).to_numpy(dtype=float)
        s = _signature(x_model, y)
        df[label] = s
        ok = np.isfinite(s) & ok_f
        summary[label] = dict(rho_with_total=float(stats.spearmanr(s[ok], total[ok])[0]),
                              rms_from_flies=float(np.sqrt(np.mean((s[ok] - sig[ok]) ** 2))),
                              n_odors_inside=int(((s >= lo) & (s <= hi))[ok].sum()), **kws)
    d = df[A1].to_numpy() - df['family'].to_numpy()
    ok = np.isfinite(d) & np.isfinite(half)
    summary['check'] = dict(rms_difference_family_A1=float(np.sqrt(np.mean(d[ok] ** 2))),
                            rms_flies_half_width=float(np.sqrt(np.mean(half[ok] ** 2))),
                            odors_apart=[o for o, dd, hh in zip(odors, d, half) if abs(dd) > hh],
                            reading=AL.RUNG_STEP_READING, n_boot=n_boot, seed=seed, ci=SIGNATURE_CI)
    summary['check']['told_apart'] = bool(summary['check']['rms_difference_family_A1']
                                          > summary['check']['rms_flies_half_width'])
    df.attrs['summary'] = summary
    return df


# ---------------------------------------------------------------------------
# the score

def stages_read(F) -> Tuple[str, ...]:
    """The stages the score reads: E136's `readable_stages` on the measure in force, without its write into E136's
    folder (E136's `stage_classes_on_measure`, written to this experiment's stage_classes.json); RuntimeError if a stage
    is not readable."""
    cl = F.stage_classes_on_measure(F.MEASURES[F.MEASURE])
    _write_json(out_dir() / 'stage_classes.json', dict(measure=F.MEASURE, classes=cl, time=_now(), host=_host(),
                                                       **F._commits()))
    bad = [s for s in F.LOSS_STAGES if cl[s]['cls'] in ('not readable', 'not scored')]
    if bad:
        raise RuntimeError(f'on the measure {F.MEASURE!r} the stages {bad} are not readable: {cl}')
    return tuple(F.LOSS_STAGES)


def _truthy(x) -> bool:
    """A flag as read back from a CSV (a bool, or its text where the column holds missing values too)."""
    return str(x).strip().lower() in ('true', '1', '1.0')


def holds_failed(name: str, reads: pd.DataFrame, apl16_runs: Mapping[str, object], tol: float) -> List[str]:
    """The reads of a variant whose holds fail: E139's reads at a bracket's bound or whose rate is more than `tol`
    (relative) from its target, and for F's APL16 rows E102's points at a bound or off their counts; as tags."""
    bad = []
    if len(reads) and 'variant' in reads.columns:
        for r in reads[reads.variant == name].itertuples():
            miss = [abs(float(getattr(r, f'rate_{c}')) - float(getattr(r, f'frac_{c}_target')))
                    / float(getattr(r, f'frac_{c}_target')) for c in ('on', 'off')]
            if _truthy(r.at_bound_thr) or _truthy(r.at_bound_apl) or not max(miss) <= tol:
                bad.append(str(r.tag))
    for w, run in apl16_runs.items():
        info = getattr(run, 'info', None) or {}
        if getattr(run, 'source', '') == 'E102:B' and (info.get('at_bound_thr') or info.get('at_bound_apl')
                                                         or info.get('n_on') != info.get('target_on')
                                                         or info.get('n_off') != info.get('target_off')):
            bad.append(run.tag)
    return bad


def part_score(names: Optional[Sequence[str]] = None) -> dict:
    """The score: every variant (default those of steps 1 and 2) through E136's score functions as a registered variant
    (`register_variants`), the unrefitted family as the reference and E136's A1 beside (its S wirings its fit's best
    evaluation, as E136 reads it): the megamat stages on all eight wirings against the fit half, the held-out half and
    all flies (score_megamat.csv; per wiring, 'common' the eight wirings' mean views, 'mean_S' / 'mean_T' the means of
    each wiring's errors); the held-out panels' KC population errors on the measure and on every view
    (score_panels.csv, score_panels_all_views.csv); the APL16 rows (E136's `apl16_rows` with the price's statistics,
    `apl_effect_parts`; F's are E102's saved runs; score_apl16.csv); the mixture panels under E137's readings
    (score_panels_e137.csv; `score_e137_beside`, `e137_beside`); the megamat somata with and without aa's entry
    (score_without_aa.csv, classes.json's beside.kc_soma_without_aa; E136's `stage_without_odor`); E136's classes (`classify`, on T against the held-out
    half, the spreads over the family's eight wirings, rule 3's guard on every view; no permuted control, no bound:
    nothing is fitted), with no class for a variant whose holds fail (`holds_failed`); the price on the APL effect
    (`price_of`), "serves both" (holds, no price); the rung named per setting on S against the fit half; step 3's
    candidates. classes.json and provenance_score.json."""
    from kcq import stage_readings as R
    t0 = time.time()
    with e136() as F:
        gate = F.require_gate()
        stages = stages_read(F)
        names = list(names) if names is not None else variants()
        register_variants(F, names)
        allv = [F.FAMILY, *names, A1]
        mrows, mega, mega_fit = [], {}, {}
        for v in allv:
            for half in ('held', None, 'fit'):
                sc = F.score_megamat(v, stages, half, F.WIRINGS)
                if sc is None:
                    continue
                vals = F.megamat_values(sc, stages, F.WIRINGS)
                if half == 'held':
                    mega[v] = vals
                if half == 'fit':
                    mega_fit[v] = vals
                base = dict(variant=v, setting=_setting_of(v), rung=_rung_of(v), half=half or 'all')
                for w, x in vals.items():
                    mrows.append(dict(base, wiring=w, **x))
                for lab, ws in (('mean_S', F.S_WIRINGS), ('mean_T', F.T_WIRINGS)):
                    mrows.append(dict(base, wiring=lab, **pd.DataFrame([vals[w] for w in ws]).mean().to_dict()))
                mrows.append(dict(base, wiring='common', stage=sc['stage_error'], downstream=sc['stage_error_downstream'],
                                  **{s: sc[f'{s}_error'] for s in stages},
                                  **{k: x for k, x in sc.items() if '_split_' in k
                                     and not any(k.endswith(f'_{w}') for w in F.WIRINGS)}))
        prow = [F.score_panels(v) for v in allv]
        panels = pd.concat([p for p in prow if len(p)], ignore_index=True) if any(len(p) for p in prow) else pd.DataFrame()
        pv = {v: F.panel_values(panels, v) for v in allv}
        fam = mega.get(F.FAMILY, {})
        spreads = dict(stage=F._sd(x['stage'] for x in fam.values()), downstream=F._sd(x['downstream'] for x in fam.values()),
                       **{s: F._sd(x[s] for x in fam.values()) for s in stages})
        for key in F.SCORED_PANELS:
            spreads['/'.join(key)] = F._sd(pv[F.FAMILY].get(key, {}).values())
        # rule 3's guard: the same reads under every view, with their own spreads
        g_mega, g_rows = {}, []
        for v in allv:
            sc = F.score_megamat(v, stages, 'held', F.WIRINGS, views=None)
            if sc is not None:
                g_mega[v] = F.megamat_values(sc, stages, F.WIRINGS)
            p = F.score_panels(v, views=None)
            if len(p):
                g_rows.append(p)
        g_panels = pd.concat(g_rows, ignore_index=True) if g_rows else pd.DataFrame()
        g_pv = {v: F.panel_values(g_panels, v) for v in allv}
        gf = g_mega.get(F.FAMILY, {})
        g_spreads = dict(stage=F._sd(x['stage'] for x in gf.values()), **{s: F._sd(x[s] for x in gf.values())
                                                                          for s in stages})
        for key in F.SCORED_PANELS:
            g_spreads['/'.join(key)] = F._sd(g_pv[F.FAMILY].get(key, {}).values())
        guard = dict(mega=g_mega, panels=g_pv, spreads=g_spreads)
        # the APL16 rows and the price
        apl_rows, parts, apl16_runs = [], {}, {}
        for v in allv:
            r = F.apl16_rows(v)
            runs = {w: F.load_read(F.ReadPoint(v, 'apl16', w, None)) for w in F.WIRINGS}
            apl16_runs[v] = {w: x for w, x in runs.items() if x is not None}
            parts[v] = {w: apl_effect_parts(run_outputs(x)) for w, x in apl16_runs[v].items()}
            if len(r):
                r = r.assign(variant=v, setting=_setting_of(v), rung=_rung_of(v))
                for k in sorted({k for d in parts[v].values() for k in d}):
                    r[f'price_{k}'] = r.wiring.map({w: d.get(k, np.nan) for w, d in parts[v].items()})
                apl_rows.append(r)
        apl = pd.concat(apl_rows, ignore_index=True) if apl_rows else pd.DataFrame()
        realism = {v: F.realism_of(apl[apl.variant == v]) if len(apl) else {} for v in allv}
        # beside: the two-view error, the mixture panels under E137's readings, E136's yardsticks
        tv = {}
        for v in allv:
            sc = F.score_megamat(v, stages, 'held', F.T_WIRINGS, views=F.MEASURES['strength_sharing'])
            if sc is not None:
                tv[v] = {w: sc[f'stage_error_{w}'] for w in F.T_WIRINGS}
        e_rows = [F.score_e137_beside(v) for v in allv]
        e137 = pd.concat([r for r in e_rows if len(r)], ignore_index=True) if any(len(r) for r in e_rows) \
            else pd.DataFrame()
        beside = dict(two_view={v: F.compare(tv[v], tv[F.FAMILY], np.nan) for v in allv[1:] if v in tv and F.FAMILY in tv},
                      e137_readings=F.e137_beside(e137, allv[1:]), panels_beside=['/'.join(k) for k in F.BESIDE_PANELS],
                      spreads_all_views=g_spreads, yardsticks=F.yardsticks(stages))
        # beside: the megamat somata with and without aa's entry (E140's outcome (d); E136's `stage_without_odor`, E138's
        # definition), on T against the held-out half, with the change from the family on the same wirings
        wo_rows = [F.stage_without_odor(v, odor=WITHOUT_ODOR) for v in allv]
        wo = pd.concat([r for r in wo_rows if r is not None], ignore_index=True) if any(r is not None for r in wo_rows) \
            else pd.DataFrame()
        if len(wo):
            wo.to_csv(out_dir() / f'score_without_{WITHOUT_ODOR}.csv', index=False)
            ft = wo[(wo.arm == F.FAMILY) & wo.row.isin(F.T_WIRINGS)].set_index('row')
            rec_wo = {}
            for v, g in wo.groupby('arm', sort=False):
                tw = g[g.row.isin(F.T_WIRINGS)].set_index('row')
                rec_wo[v] = dict(with_mean=float(tw.view_error.mean()), without_mean=float(tw.view_error_without.mean()),
                                 share_mean=float(tw.view_share.mean()))
                if v != F.FAMILY and len(ft):
                    rec_wo[v].update(change_with=float((tw.view_error - ft.view_error.reindex(tw.index)).mean()),
                                     change_without=float((tw.view_error_without
                                                           - ft.view_error_without.reindex(tw.index)).mean()))
            beside[f'kc_soma_without_{WITHOUT_ODOR}'] = rec_wo
        # the classes
        reads = reads_table()
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
            failed = holds_failed(v, reads, apl16_runs.get(v, {}), F.HOLD_TOL)
            if failed:
                rec['reasons'].append('a hold fails: ' + ', '.join(failed))
                rec['cls'] = 'no class'
            rec['holds_failed'] = failed
            rec['flagged_wirings'] = {w: f for w, f in realism.get(v, {}).items() if f}
            rec['price_on_apl_effect'] = price_of(parts, v, F.FAMILY)
            pr = rec['price_on_apl_effect'].get('price')
            rec['serves_both'] = bool(rec['cls'] == 'holds' and pr is False)
            rec['setting'], rec['rung'] = _setting_of(v), _rung_of(v)
            classes[v] = rec
        named = name_rungs(mega_fit, names, F.S_WIRINGS)
        a1_price = (classes.get(A1, {}).get('price_on_apl_effect') or {}).get('price_value')
        cands = []
        for v in names:
            s, r = parse_variant(v)
            c = classes.get(v, {})
            pval = (c.get('price_on_apl_effect') or {}).get('price_value')
            if s in ('F', 'B') and r < CANDIDATE_BELOW and c.get('improves') and pval is not None and a1_price is not None \
                    and pval < CANDIDATE_PRICE_SHARE * a1_price:
                cands.append(dict(variant=v, setting=s, rung=r, price_value=pval, a1_price_value=a1_price))
        e136_spreads = e136_classes.get('spreads') or {}
        sp_diff = {k: abs(spreads[k] - e136_spreads[k]) for k in spreads if k in e136_spreads
                   and np.isfinite(spreads[k]) and np.isfinite(e136_spreads[k])}
        mdf = pd.DataFrame(mrows)
        mdf.to_csv(out_dir() / 'score_megamat.csv', index=False)
        panels.to_csv(out_dir() / 'score_panels.csv', index=False)
        g_panels.to_csv(out_dir() / 'score_panels_all_views.csv', index=False)
        e137.to_csv(out_dir() / 'score_panels_e137.csv', index=False)
        apl.to_csv(out_dir() / 'score_apl16.csv', index=False)
        out = dict(stages=list(stages), measure=F.MEASURE, variants=names, reference=F.FAMILY, beside_arm=A1,
                   spreads=spreads, spreads_max_abs_diff_from_e136=max(sp_diff.values()) if sp_diff else None,
                   price_rule=dict(threshold=PRICE, targets=list(PRICE_TARGETS), decided_on='T'),
                   named_rungs=named, step3_candidates=cands, beside=beside, classes=classes)
        _write_json(out_dir() / 'classes.json', out)
        _write_json(out_dir() / 'provenance_score.json', dict(
            part='score', script=SCRIPT, lock=R.require_lock().get('sha256'), gate_measure=gate['gate'].get('measure'),
            seconds=round(time.time() - t0, 1), host=_host(), time=_now(), variants=names, **F._commits()))
    return {v: c.get('cls') for v, c in classes.items()}


def _setting_of(v: str) -> str:
    try:
        return parse_variant(v)[0]
    except ValueError:
        return v


def _rung_of(v: str) -> float:
    try:
        return parse_variant(v)[1]
    except ValueError:
        return np.nan


def name_rungs(mega_fit: Mapping[str, Mapping[str, Mapping[str, float]]], names: Sequence[str],
               wirings: Sequence[str]) -> Dict[str, dict]:
    """Per setting, its rung named on the S wirings against the fit half: the rung whose mean over `wirings` of each
    wiring's megamat stage error is lowest (with every rung's value; settings with one rung read it alone)."""
    out: Dict[str, dict] = {}
    for v in names:
        if v not in mega_fit or not all(w in mega_fit[v] for w in wirings):
            continue
        s, r = parse_variant(v)
        out.setdefault(s, dict(values={}))['values'][f'{r:g}'] = float(np.mean([mega_fit[v][w]['stage']
                                                                              for w in wirings]))
    for s, rec in out.items():
        rec['named'] = float(min(rec['values'], key=rec['values'].get))
    return out


def write_variants(arch: Optional[Path] = None) -> pd.DataFrame:
    """variants.csv in out_dir(), for the model-output archive's params.csv (as E136's `write_variants`): one row per
    read the archive holds (`arch`, default archive_dir()/read/), from its <tag>_run.json: the tag, the variant, setting
    and rung, the panel and wiring, a description, the coordinates (z_<name>), both conditions' absolute parameters
    (E136's VARIANT_PARAMS; 'off_<name>' with APL off), the PN normalization's terms, and where and when it ran."""
    from experiments import E136_stage_population_fit as F
    arch = Path(arch) if arch is not None else archive_dir() / 'read'
    rows = []
    for p in sorted(arch.glob('*_run.json')):
        d = json.loads(p.read_text())
        r = dict(tag=d['tag'], variant=d.get('variant'), setting=d.get('setting'), rung=d.get('rung'), arm=d.get('arm'),
                 panel=d.get('panel'), wiring=d.get('wiring'), reading=d.get('reading'), model=d.get('model'),
                 kind=f"setting {d.get('setting')} (the family's calyx factors x {SETTINGS.get(d.get('setting'), {})}) at "
                      f"pn_inhadd {d.get('rung')}, read on the {d.get('panel')} input, wiring {d.get('wiring')}, both "
                      f"response rates held")
        r.update({f'z_{k}': v for k, v in (d.get('z') or {}).items()})
        r.update({k: (d.get('params') or {}).get(k) for k in F.VARIANT_PARAMS})
        r.update({f'off_{k}': (d.get('params_off') or {}).get(k) for k in F.VARIANT_PARAMS})
        r.update(pn_inhadd=d.get('pn_inhadd'), pn_inhsc=d.get('pn_inhsc'), host=d.get('host'), lane=d.get('lane'),
                 time=d.get('time'))
        rows.append(r)
    df = pd.DataFrame(rows)
    if df.empty:
        raise FileNotFoundError(f'no <tag>_run.json under {arch}')
    df.to_csv(out_dir() / 'variants.csv', index=False)
    return df


def part_fit(rung: float) -> dict:
    """Step 3, only on a candidate: A0's two-factor fit (E136's `fit_arm`) with the normalization fixed at `rung`, on S
    against the fit half. Not written yet: NotImplementedError."""
    raise NotImplementedError(
        f"E139's step 3 (A0's two-factor fit with pn_inhadd fixed at {rung:g} and pn_inhsc by inhsc_for) is not written: "
        "it needs an arm of E136 that is A0 with the normalization's kwargs fixed (ArmModel takes the arm's kwargs from "
        "E136's ARMS), its fit's records in this experiment's folder rather than E136's (fit_arm writes under E136's "
        "out_dir), and E136's reads after it; it runs only on a candidate (classes.json 'step3_candidates')")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('step0', 'read', 'score', 'fit', 'variants'))
    ap.add_argument('--step', nargs='+', default=['1', '2'], choices=('1', '2', 'beside'))
    ap.add_argument('--lane', default='1')
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--reclaim', action='store_true')
    ap.add_argument('--rung', type=float)
    ap.add_argument('--tables', nargs='+', choices=('a1_by_rung', 'e102_ladder', 'pn_stage_ladder',
                                                    'recorded_step_signature'),
                    help='--part step0: the tables to write (default all)')
    a = ap.parse_args(argv)
    if a.part == 'step0':
        fns = {f.__name__: f for f in (a1_by_rung, e102_ladder, pn_stage_ladder, recorded_step_signature)}
        for name in (a.tables or list(fns)):
            df = fns[name]()
            df.to_csv(out_dir() / f'step0_{name}.csv', index=False)
            if df.attrs.get('summary') is not None:
                _write_json(out_dir() / f'step0_{name}.json', dict(df.attrs['summary'], time=_now(), host=_host()))
            print(f'E139 step 0: step0_{name}.csv ({len(df)} rows)', flush=True)
    elif a.part == 'read':
        with e136() as F:
            F.require_gate()
        part_read(a.lane, a.step, a.reverse, a.reclaim)
    elif a.part == 'score':
        print(json.dumps(part_score(), indent=1, default=str))
    elif a.part == 'fit':
        if a.rung is None:
            raise SystemExit('--part fit needs --rung')
        print(json.dumps(part_fit(a.rung), indent=1, default=str))
    elif a.part == 'variants':
        print(f'variants.csv: {len(write_variants())} runs', flush=True)


if __name__ == '__main__':
    main()
