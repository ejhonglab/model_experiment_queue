#!/usr/bin/env python3
"""E86: the reduced searched input as a gated candidate. E77's input searched on the real wiring against one half of the
recordings, at 0.6 of its log size, read on the four held-out wirings (T: P3, P4, M2, M3) against the half its search
never saw: the label-permuted and the inverted controls, the search repeated with a second draw and on a second split of
the recordings, and a verdict fixed before any run.

The entry is item E86 of experiment_priorities.md (planned by a Fable subagent, claude-fable-5-2, spawned from hal
2026-11-12 11:08). C2's family throughout (the breadth-2 APL-effect fit of E67), each wiring at its own fit of E72's part
B, the 0.5 s pulse; S: real, P1, P2, M1; T: P3, P4, M2, M3 (E83's division). Halves A and B are E77's
(`recording_halves()`, seed 0); halves C and D are `recording_halves(1)`'s 'A' and 'B'. Every input is read at 0.6 of its
log factors (E83's `input_at`) with both response rates re-held (E83's `input_reading`: the breadth recomputed for the
input, the threshold bisected to the wiring's APL-off response count, then the APL scale, KC->APL and bouton->APL tied,
to its APL-on count), its breadth run cached under outputs/E86_reduced_input_candidate/ (`held_reading`). For an input, g
is T's pattern part under the calibrated readout, unchanged minus changed, against the half of its split that its search
never saw (positive is better). The candidate's readings on T are E83's part A2 runs (archived), read, not rerun.

Parts (CLI `--part`):
- A (no simulations): the four fitted inputs (E77's two searched inputs, E83 arm I's two) at 0.6 and at full size, by
  the pairs they move (size, projection on what the unseen half asks, the ester / alcohol / other groups' shares of the
  squared pattern error), wiring by wiring as paired differences, the agreement of their factors and of their per-pair
  changes, from E83's pairs.csv; the searches' step sizes and factors at a bound; E85's part A4 table cited, not
  recomputed. Writes partA_pairs.csv, partA_agreement.csv, partA_searches.csv, part_a.json.
- identity: T's simulators reproduce their saved fits; E83's part A2 point (P3, half A's input at 0.6) is reproduced
  exactly (counts, responses, both factors, the breadth). Writes identity.json.
- B1: E83's `permuted_factors` of the candidate's factors, seeds 1 and 2, at 0.6 on T, both directions (16 points).
- B2: the inverted input, the candidate's factors at -0.6 of their log size (8 points).
- B3: seeds 3 and 4, only if a B1 input has g above 0 (16 points).
- B: B1 and B2, scored, then B3 if B1 asks for it.
- C-check: E77's part E search (`E77_input_side_common_error.input_oracle`), first split, generator seed 0, the first
  generation of eb only: must reproduce E77's oracle_history_A.csv to the last digit (check_search.json).
- C (`--search C1_A|C1_B|C2_C|C2_D`): E77's part E search as it was, on the real wiring: C1 on halves A and B with
  generator seeds 10 and 11, C2 on halves C and D with generator seeds 20 and 21 (input_oracle seeds its generator with
  `seed` for the split's first half and `seed + 1` for its second, so C1 passes seed=10 and C2 seed=20 for both halves;
  SEARCHES, GENERATOR_SEEDS); a checkpoint per generation, resumed by default. Writes oracle_history_<name>.csv,
  oracle_factors_<name>.csv, inputs/oracle_<name>.csv, provenance_C_<name>.json.
- C-read: each of the four new inputs at 0.6 on T, both response rates re-held (16 points).
- score: readings.csv, pairs.csv, two_by_two.csv and inputs.csv (one row per input: g and the gate's statistics), the
  archive tables (model_pair_distances_long.csv, variants.csv), runs.csv.
- D: score, then the verdict (`verdict`): gate.csv, replication.csv, null.csv, verdict.json, comparison_table.md.
- E (only on "passes" or "carries"): what the input is. The six inputs searched on the real wiring (SIX: E77's two, the
  candidate, whose readings at 0.6 on all eight wirings are E83's part A2 runs, read, not rerun; part C's four, on T
  from C-read) are scored on every population and wiring (partE_readings.csv, partE_pairs.csv, partE_two_by_two.csv),
  the seven re-drawn wirings (P1, P2, M1, P3, P4, M2, M3) as paired differences against the unseen half of each input's
  split (partE_wirings.csv, partE_inputs.csv); the agreement of the six (log2 factors over all entries and per odor,
  per-pair changes on T, the seven and all eight; partE_agreement.csv, partE_agreement_by_odor.csv); the entries on
  which five of the six agree in sign with a mean log2 factor beyond +-0.5 (the list a measurement could test, named
  only with ten or more; partE_entries.csv, partE_five_of_six.csv, with the rule's label-permuted counts beside it in
  partE_null.csv); the pairs the six move (partE_pairs_moved.csv); the six against the ORN flies' spread and Hallem's
  data through E85's part A4 code (partE_vs_flies*.csv); part_e.json, the "Part E" section of comparison_table.md.
  These tables need no simulation and are written first; then the four new inputs at 0.6 on S (16 held readings,
  runs_E.csv, tags E_<input>_<wiring>), and the tables again on all eight wirings. E-tables: the tables only.
- F: only on "passes" (not written).

Outputs: results/E86_reduced_input_candidate/ (shared directory; also runs_<part>.csv with absolute parameters and
realism per point, provenance_<part>.json with `H.model_provenance` per point, inputs/<change>.csv in E38's layout) and
outputs/E86_reduced_input_candidate/ (local: spike counts and responses per run under <part>/ (tag <change>_<wiring>,
conditions on / acute_off), the searches' checkpoints under C/, the breadth runs).

Usage (from the repository root):
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E86_reduced_input_candidate.py --part B
    OMP_NUM_THREADS=8 OMP_THREAD_LIMIT=8 nice -n 10 ./venv/bin/python -u experiments/E86_reduced_input_candidate.py --part C --search C1_A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E86_reduced_input_candidate.py --part E
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E77_input_side_common_error as E77
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E86_reduced_input_candidate'
SCRIPT = 'experiments/E86_reduced_input_candidate.py'
SIZE = 0.6
INVERT_SIZE = -0.6
T_WIRINGS = ('P3', 'P4', 'M2', 'M3')
B1_SEEDS = (1, 2)
B3_SEEDS = (3, 4)
# the halves of each split of the recordings: split seed -> {half name here: the half of recording_halves(seed)}
SPLITS = {0: {'A': 'A', 'B': 'B'}, 1: {'C': 'A', 'D': 'B'}}
UNSEEN = {'A': 'B', 'B': 'A', 'C': 'D', 'D': 'C'}
# part C's searches: name -> the split, the half searched (its name here) and input_oracle's `seed` (input_oracle seeds
# its generator with seed for the split's first half and seed + 1 for its second: GENERATOR_SEEDS)
SEARCHES = {'C1_A': dict(split_seed=0, half='A', seed=10), 'C1_B': dict(split_seed=0, half='B', seed=10),
            'C2_C': dict(split_seed=1, half='C', seed=20), 'C2_D': dict(split_seed=1, half='D', seed=20)}
GENERATOR_SEEDS = {'C1_A': 10, 'C1_B': 11, 'C2_C': 20, 'C2_D': 21}
# part C's check: E77's search of half A (first split, generator seed 0), the first generation of its first odor
CHECK = dict(half='A', seed=0, split_seed=0, odors=('eb',), gens=1)
CHECK_COLUMNS = ('half', 'odor', 'generation', 'best_obj', 'parent_obj', 'success', 'sigma', 'pattern_rms_half',
                 'pattern_rms_other', 'frac_on', 'frac_off', 'c2_frac_on', 'c2_frac_off')
# the columns computed from the pattern error (floating-point sums over the 120 pairs), and the relative difference they
# may have from E77's when everything else (the generator's path, every count-derived column) is exact: the readout's
# distances changed in their last digits with the move of the harness's functions to hong2p and pyolfsysm (2026-11-12;
# C2's own pattern part on the real wiring against half A is 0.0609269329706763 in E83's readings.csv of 2026-11-11 and
# 0.06092693297067637 now, on hal, with the same numpy and pandas); on hal the check's differ by 4 to 8 ulps (1e-15)
CHECK_FLOAT_COLUMNS = ('best_obj', 'parent_obj', 'pattern_rms_half', 'pattern_rms_other')
CHECK_REL_TOL = 1e-12
# the four fitted inputs of part A: name -> (results folder, factor table, E83's change prefix, half searched, label)
FITTED = {'E77_A': ('E77_input_side_common_error', 'oracle_factors_A.csv', 'A2_halfA', 'A', "E77's, half A"),
          'E77_B': ('E77_input_side_common_error', 'oracle_factors_B.csv', 'A2_halfB', 'B', "E77's, half B"),
          'I_A': ('E83_cross_wiring_fit', 'input_factors_A.csv', 'DI_halfA', 'A', "E83's arm I, half A"),
          'I_B': ('E83_cross_wiring_fit', 'input_factors_B.csv', 'DI_halfB', 'B', "E83's arm I, half B")}
FITTED_SEARCHED_ON = {'E77_A': 'the real wiring', 'E77_B': 'the real wiring', 'I_A': 'real, P1, P2, M1',
                      'I_B': 'real, P1, P2, M1'}
SEARCH_HISTORIES = {'E77_A': ('E77_input_side_common_error', 'oracle_history_A.csv'),
                    'E77_B': ('E77_input_side_common_error', 'oracle_history_B.csv'),
                    'I_A': ('E83_cross_wiring_fit', 'search_history_I_A.csv'),
                    'I_B': ('E83_cross_wiring_fit', 'search_history_I_B.csv')}
E85_A4 = ('E85_orn_fly_population', 'searched_inputs_summary.csv')
E85_A4_NAMES = {'E77_A': 'E77_halfA', 'E77_B': 'E77_halfB', 'I_A': 'E83_armI_halfA', 'I_B': 'E83_armI_halfB'}
N_NULL = 1000
# the verdict's thresholds (the entry's part D, fixed before any run)
CONTROL_SHARE, B3_MEAN_SHARE = 0.5, 0.25
GATE_GAIN, GATE_SINGLES = 0.003, 12
GATE_RAW_SLACK, GATE_MEAN_SLACK, GATE_SQ_ALL, GATE_SQ_SIX = 0.005, 0.01, 0.15, 0.075
PASS_MEAN, CARRY_MEAN, MIN_POSITIVE = 0.002, 0.001, 3
# the groups of inputs: the candidate (E77's inputs at 0.6, E83's archived part A2 runs), the controls, the new searches
GROUPS = ('candidate', 'B1', 'B2', 'B3', 'C1', 'C2')
PARTS = ('B1', 'B2', 'B3', 'Cread')


# ---------------------------------------------------------------------------
# halves, targets, inputs

_HALVES: Dict[str, Dict[str, list]] = {}
_TARGETS: Dict[str, pd.Series] = {}


def halves() -> Dict[str, Dict[str, list]]:
    """The four halves of the recordings: {'A', 'B', 'C', 'D'}: {'ctrl' / 'tnt': [recording keys]} (A, B:
    `recording_halves(0)`; C, D: `recording_halves(1)`'s 'A' and 'B')."""
    if not _HALVES:
        for seed, names in SPLITS.items():
            hv = E77.recording_halves(seed)
            for name, h in names.items():
                _HALVES[name] = hv[h]
    return dict(_HALVES)


def split_of(half: str) -> int:
    """The split seed a half belongs to."""
    for seed, names in SPLITS.items():
        if half in names:
            return seed
    raise ValueError(half)


def targets() -> Dict[str, pd.Series]:
    """E83's targets (halfA, halfB, z8, z6, noise_free, noise8) and halfC, halfD (canonical pair order)."""
    if not _TARGETS:
        t = dict(E83.targets())
        hv = halves()
        for h in ('C', 'D'):
            t[f'half{h}'] = E83._cs(E77.half_effect(hv[h]))
        for k, v in t.items():
            if len(v) != 120 or v.isna().any():
                raise RuntimeError(f'target {k}: {len(v)} pairs, {int(v.isna().sum())} missing')
        _TARGETS.update(t)
    return dict(_TARGETS)


def candidate_factors(half: str) -> pd.DataFrame:
    """E77's searched factors on half `half` ('A' or 'B'), glomeruli x odors (E83's `read_factors`)."""
    return E83.read_factors(H.results_dir('E77_input_side_common_error') / f'oracle_factors_{half}.csv')


def control_factors(kind: str, half: str, seed: Optional[int] = None) -> Tuple[pd.DataFrame, float]:
    """A control input's factors and size: kind 'perm' E83's `permuted_factors` of the candidate's factors (seed) at
    SIZE; 'inverted' the candidate's factors at INVERT_SIZE."""
    F = candidate_factors(half)
    if kind == 'perm':
        if seed is None:
            raise ValueError("a permuted control needs a seed")
        return E83.permuted_factors(F, int(seed)), SIZE
    if kind == 'inverted':
        return F, INVERT_SIZE
    raise ValueError(kind)


def search_factors(out: Path, name: str) -> pd.DataFrame:
    """A part-C search's factors (oracle_factors_<name>.csv), glomeruli x odors; raises unless all ten odors are done."""
    p = out / f'oracle_factors_{name}.csv'
    if not p.exists():
        raise FileNotFoundError(p)
    done = set(pd.read_csv(p).odor)
    if done != set(E77.E_ODORS):
        raise RuntimeError(f'{p}: odors {sorted(done)} of {sorted(E77.E_ODORS)} searched')
    return E83.read_factors(p)


def _local(part: Optional[str] = None) -> Path:
    d = H.outputs_dir(EXPERIMENT)
    if part is not None:
        d = d / part
        d.mkdir(parents=True, exist_ok=True)
    return d


def build(label: str, orn: pd.DataFrame, input_tag: str):
    """E83's `build_simulator` with the breadth recomputed for `orn`, its APL-off breadth run cached under this
    experiment's outputs (by `input_tag`); the wiring and the normalized APL weights must be the lab input's."""
    from experiments import E10_joint_fit as J
    P = E83.wiring_parts(label)
    _, wnorm, W, breadth = J.bouton_base(orn, _local(), None, input_tag=input_tag, **E83.wiring_kws(label))
    if not W.equals(P['W']):
        raise RuntimeError(f"{label}: the wiring built with input {input_tag} differs from the lab input's")
    for k in ('wAPLKC', 'wKCAPL', 'wAPLPN', 'wPNAPL'):
        if not P['wnorm'][k].equals(wnorm[k]):
            raise RuntimeError(f"{label}: {k} built with input {input_tag} differs from the lab input's")
    return E83.build_simulator(label, orn, breadth=breadth)


def held_reading(label: str, orn: pd.DataFrame, input_tag: str) -> dict:
    """E83's `input_reading` (breadth recomputed, threshold then APL scale re-held; `E77_input_screen.hold_rates`) with
    this experiment's breadth cache."""
    from experiments import E77_input_screen as S
    P = E83.saved_fit(label)
    res = S.hold_rates(lambda: build(label, orn, input_tag), target_off=P['target_off'], target_on=P['target_on'],
                       params=P['p'], label=f'{label} {input_tag}')
    res['input_tag'] = input_tag
    return res


def save_run(part: str, tag: str, o: Dict[str, tuple]) -> None:
    """Spike counts and responses of a run under outputs/E86_reduced_input_candidate/<part>/."""
    d = _local(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(part: str, tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run of this experiment ({'on', 'off'}): the local outputs, else the model-output archive; or None."""
    for d in (H.outputs_dir(EXPERIMENT) / part, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def _rel(path: Path) -> str:
    """A path relative to the shared directory when it is inside it."""
    try:
        return str(Path(path).relative_to(H.shared_dir()))
    except ValueError:
        return str(path)


def _provenance_add(path: Path, part: str, tag: str, prov: dict) -> None:
    d = json.loads(path.read_text()) if path.exists() else dict(part=part, script=SCRIPT, points={})
    d['points'][tag] = prov
    path.write_text(json.dumps(d, indent=1, default=str) + '\n')


# ---------------------------------------------------------------------------
# part A (no simulations)

def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def _e83_effects(readout: str = 'cal', path: Optional[Path] = None) -> Dict[Tuple[str, str], pd.Series]:
    """{(change, unit): per-pair APL effect (canonical order)} from E83's pairs.csv (or `path`, a table of the same
    layout, e.g. part E's partE_pairs.csv)."""
    pr = pd.read_csv(H.results_dir(E83.EXPERIMENT) / 'pairs.csv' if path is None else path, float_precision='round_trip')
    pr = pr[pr.readout == readout]
    out = {}
    for (ch, unit), g in pr.groupby(['change', 'unit'], sort=False):
        out[(ch, unit)] = E83._cs(g.set_index(['odor_a', 'odor_b'])['effect'])
    return out


def change_stats(m0: pd.Series, m1: pd.Series, target: pd.Series) -> dict:
    """A change of the per-pair APL effect (m0 unchanged, m1 changed; canonical order) against a target effect: the
    pattern part unchanged / changed and the gain; the change's rms (raw and with its mean removed: the size); its
    correlation with what the target asks (the unchanged pattern error with its sign turned, data minus model, both
    means removed) and the share of that request it delivers (regression coefficient of the centered change on it);
    the change of the squared pattern error in the ester group (15 pairs), the alcohol group (6) and the other 99, each
    as a share of the unchanged total (they sum to the total's relative change); the change of the mean effect."""
    m0, m1 = E83._cs(m0), E83._cs(m1).reindex(E83._cs(m0).index)
    de = E83._cs(target).reindex(m0.index)
    e0 = (m0 - m0.mean()) - (de - de.mean())
    e1 = (m1 - m1.mean()) - (de - de.mean())
    d = (m1 - m0).to_numpy()
    dc = d - d.mean()
    a = -e0.to_numpy()
    masks = E77.kind_masks(m0.index)
    other = ~(masks['ester'] | masks['alcohol'])
    tot = float((e0.to_numpy() ** 2).sum())
    row = dict(unchanged=_rms(e0), changed=_rms(e1), gain=_rms(e0) - _rms(e1), rms_change=_rms(d), size_change=_rms(dc),
               corr_asked=float(np.corrcoef(d, a)[0, 1]), share_asked=float((dc * a).sum() / (a * a).sum()),
               d_mean_effect=float(np.mean(d)), sq_change_total=float((e1.to_numpy() ** 2).sum()) / tot - 1.0)
    for name, mk in (('ester', masks['ester']), ('alcohol', masks['alcohol']), ('other', other)):
        row[f'sq_change_{name}'] = float((e1.to_numpy()[mk] ** 2).sum() - (e0.to_numpy()[mk] ** 2).sum()) / tot
        row[f'n_{name}'] = int(mk.sum())
    return row


def variance_shares(L: pd.DataFrame) -> Dict[str, float]:
    """The share of a table of log factors' variance (glomeruli x odors, every entry used) explained by a per-glomerulus
    term, a per-odor term, and both (additive two-way fit), as E83's `agreement`."""
    y = L.to_numpy(dtype=float)
    tot = float(((y - y.mean()) ** 2).sum())
    g = y.mean(axis=1, keepdims=True)
    o = y.mean(axis=0, keepdims=True)
    both = g + o - y.mean()
    return dict(glomerulus=1.0 - float(((y - g) ** 2).sum()) / tot, odor=1.0 - float(((y - o) ** 2).sum()) / tot,
                both=1.0 - float(((y - both) ** 2).sum()) / tot)


def fitted_factors(name: str) -> pd.DataFrame:
    """A fitted input's factors (FITTED), glomeruli x the ten searched odors."""
    folder, fname, _, _, _ = FITTED[name]
    return E83.read_factors(H.results_dir(folder) / fname)[list(E77.E_ODORS)]


def part_a(out: Path) -> dict:
    """Part A (no simulations): the four fitted inputs by the pairs they move (partA_pairs.csv: per input, size, unit
    and target, `change_stats`), the agreement of their factors and per-pair changes (partA_agreement.csv), the searches'
    step sizes and factors at a bound with the variance shares and their permutation null (partA_searches.csv), E85's
    part A4 table cited (not recomputed), part_a.json and the part-A section of comparison_table.md."""
    eff = _e83_effects('cal')
    tg = targets()
    rows = []
    for name, (_, _, prefix, X, label) in FITTED.items():
        Y = UNSEEN[X]
        for size in (SIZE, 1.0):
            ch = f'{prefix}_s{size:g}'
            for unit in list(E83.POPULATIONS) + list(E83.ALL_WIRINGS):
                if (ch, unit) not in eff:
                    continue
                kind = 'population' if unit in E83.POPULATIONS else 'wiring'
                for against, t in (('unseen', f'half{Y}'), ('searched', f'half{X}'), ('all', 'z8'), ('six', 'z6')):
                    rows.append(dict(input=name, label=label, searched_on=X, unseen=Y, size=size, change=ch, unit=unit,
                                     kind=kind, against=against, target=t,
                                     **change_stats(eff[('none', unit)], eff[(ch, unit)], tg[t])))
    pa = pd.DataFrame(rows)
    pa.to_csv(out / 'partA_pairs.csv', index=False)

    # agreement: factors and per-pair changes, pairwise
    names = list(FITTED)
    L = {n: np.log2(fitted_factors(n)) for n in names}
    arows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            per_odor = [float(np.corrcoef(L[a][o], L[b][o])[0, 1]) for o in E77.E_ODORS]
            row = dict(input_a=a, input_b=b, factors_corr=float(np.corrcoef(L[a].to_numpy().ravel(), L[b].to_numpy().ravel())[0, 1]),
                       factors_corr_per_odor_median=float(np.median(per_odor)), factors_corr_per_odor_min=float(np.min(per_odor)),
                       factors_corr_per_odor_max=float(np.max(per_odor)))
            for size in (SIZE, 1.0):
                for unit in ('T', 'all8'):
                    ka, kb = (f'{FITTED[n][2]}_s{size:g}' for n in (a, b))
                    if (ka, unit) in eff and (kb, unit) in eff:
                        base = eff[('none', unit)]
                        da = (eff[(ka, unit)] - base).reindex(base.index)
                        db = (eff[(kb, unit)] - base).reindex(base.index)
                        row[f'changes_corr_{unit}_s{size:g}'] = float(np.corrcoef(da, db)[0, 1])
            arows.append(row)
    ag = pd.DataFrame(arows)
    ag.to_csv(out / 'partA_agreement.csv', index=False)

    # the searches: final step sizes, factors at a bound, variance shares against permuted factors
    srows = []
    for n in names:
        folder, fname = SEARCH_HISTORIES[n]
        hist = pd.read_csv(H.results_dir(folder) / fname)
        last = hist.sort_values('generation').groupby('odor').tail(1)
        F = fitted_factors(n)
        orn = H.abbrev_columns(H.apl16_orn_deltas())[list(E77.E_ODORS)]
        resp = orn.to_numpy(dtype=float) != 0
        f = F.to_numpy()[resp]
        at_bound = np.isclose(f, E77.E_BOUNDS[0], rtol=1e-9) | np.isclose(f, E77.E_BOUNDS[1], rtol=1e-9)
        vs = variance_shares(L[n])
        Ffull = E83.read_factors(H.results_dir(FITTED[n][0]) / FITTED[n][1])
        null = np.array([variance_shares(np.log2(E83.permuted_factors(Ffull, s)[list(E77.E_ODORS)]))['glomerulus']
                         for s in range(1000, 1000 + N_NULL)])
        srows.append(dict(input=n, n_entries=int(resp.sum()), n_responding_all=bool(resp.all()),
                          sigma_final_median=float(last.sigma.median()), sigma_final_min=float(last.sigma.min()),
                          sigma_final_max=float(last.sigma.max()), log_bound=float(np.log(E77.E_BOUNDS[1])),
                          share_at_bound=float(at_bound.mean()), var_glomerulus=vs['glomerulus'], var_odor=vs['odor'],
                          var_both=vs['both'], var_glomerulus_perm_mean=float(null.mean()),
                          var_glomerulus_perm_p95=float(np.quantile(null, 0.95)), n_perm=N_NULL,
                          var_glomerulus_perm_frac_ge=float((null >= vs['glomerulus']).mean())))
    se = pd.DataFrame(srows)
    se.to_csv(out / 'partA_searches.csv', index=False)

    e85 = pd.read_csv(H.results_dir(E85_A4[0]) / E85_A4[1])
    res = dict(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
               sources=dict(pairs=f'results/{E83.EXPERIMENT}/pairs.csv', e85_a4=f'results/{E85_A4[0]}/{E85_A4[1]}',
                            factors={n: f'results/{v[0]}/{v[1]}' for n, v in FITTED.items()},
                            histories={n: f'results/{v[0]}/{v[1]}' for n, v in SEARCH_HISTORIES.items()}))
    tu = pa[(pa.unit == 'T') & (pa.against == 'unseen')]
    res['T_unseen'] = {f'{r.input}_s{r.size:g}': dict(gain=r.gain, size=r.size, corr_asked=r.corr_asked)
                       for r in tu.itertuples()}
    res['mean_gain_T_unseen_0.6'] = float(tu[np.isclose(tu['size'], SIZE)].gain.mean())
    (out / 'part_a.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    _comparison_table_a(out, pa, ag, se, e85)
    print(tu[['input', 'size', 'gain', 'size_change', 'corr_asked', 'share_asked', 'sq_change_ester', 'sq_change_alcohol',
              'sq_change_other']].round(4).to_string(index=False), flush=True)
    return res


def _fmt(x, spec='.4f') -> str:
    return '–' if x is None or (isinstance(x, (float, np.floating)) and not np.isfinite(x)) else format(x, spec)


def _section(out: Path, heading: str, lines: List[str]) -> None:
    """Replace (or append) the comparison_table.md section starting with `heading` ('## ...')."""
    path = out / 'comparison_table.md'
    text = path.read_text() if path.exists() else '# E86, the reduced searched input as a gated candidate\n\n'
    sec = '\n'.join([heading] + lines)
    start = text.find(heading)
    if start >= 0:
        nxt = text.find('\n## ', start + 1)
        text = text[:start] + sec + ('\n' if nxt >= 0 else '') + ('' if nxt < 0 else text[nxt + 1:])
    else:
        text = text.rstrip('\n') + '\n\n' + sec
    path.write_text(text.rstrip('\n') + '\n')


def _comparison_table_a(out: Path, pa: pd.DataFrame, ag: pd.DataFrame, se: pd.DataFrame, e85: pd.DataFrame) -> None:
    L = ['', 'Pattern part of the APL-effect loss under the calibrated readout (lower is better), from E83\'s readings '
         '(pairs.csv; each wiring at its own fit, both response rates re-held). g: unchanged minus changed (positive is '
         'better). Size: the rms over the 120 pairs of the change of the per-pair APL effect, its mean removed. What the '
         'unseen half asks: its per-pair effect minus the unchanged model\'s, both means removed. The change of the '
         'squared pattern error per group is a share of the unchanged total (the three sum to its relative change).', '',
         '### On T (P3, P4, M2, M3) against the half each search never saw', '',
         '| input | searched on | share of its log size | g | size of the change | correlation with what the unseen half asks | '
         'share of it delivered | squared pattern error: ester (15) / alcohol (6) / other (99) |',
         '|---|---|---:|---:|---:|---:|---:|---:|']
    tu = pa[(pa.unit == 'T') & (pa.against == 'unseen')]
    for r in tu.to_dict('records'):
        L.append(f"| {r['label']} | {FITTED_SEARCHED_ON[r['input']]} | {r['size']:g} | {_fmt(r['gain'], '+.4f')} | "
                 f"{_fmt(r['size_change'])} | {_fmt(r['corr_asked'], '.2f')} | "
                 f"{_fmt(r['share_asked'], '.2f')} | {_fmt(100 * r['sq_change_ester'], '+.1f')} % / "
                 f"{_fmt(100 * r['sq_change_alcohol'], '+.1f')} % / {_fmt(100 * r['sq_change_other'], '+.1f')} % |")
    L += ['', '### Wiring by wiring, at 0.6, against the unseen half (changed minus unchanged; negative is better)', '',
          '| input | ' + ' | '.join(E83.ALL_WIRINGS) + ' | re-drawn wirings lower | mean over the re-drawn |',
          '|---|' + '---:|' * (len(E83.ALL_WIRINGS) + 2)]
    for name in FITTED:
        q = pa[(pa.input == name) & np.isclose(pa['size'], SIZE) & (pa.against == 'unseen') & (pa.kind == 'wiring')].set_index('unit')
        cells = [_fmt(-float(q.loc[w, 'gain']), '+.4f') if w in q.index else '–' for w in E83.ALL_WIRINGS]
        dr = q.reindex(list(E83.ALL_WIRINGS[1:])).gain.dropna()
        L.append(f"| {FITTED[name][4]} | " + ' | '.join(cells) + f" | {int((dr > 0).sum())} of {len(dr)} | {_fmt(float(dr.mean()), '+.4f')} |")
    L += ['', '### Agreement', '',
          '| inputs | log2 factors, all 360 entries | per odor, median (min to max) | per-pair change on T at 0.6 | '
          'on all eight at 0.6 | on T at full size |', '|---|---:|---:|---:|---:|---:|']
    for r in ag.to_dict('records'):
        L.append(f"| {FITTED[r['input_a']][4]} / {FITTED[r['input_b']][4]} | {_fmt(r['factors_corr'], '.2f')} | "
                 f"{_fmt(r['factors_corr_per_odor_median'], '.2f')} ({_fmt(r['factors_corr_per_odor_min'], '.2f')} to "
                 f"{_fmt(r['factors_corr_per_odor_max'], '.2f')}) | {_fmt(r.get('changes_corr_T_s0.6', np.nan), '.2f')} | "
                 f"{_fmt(r.get('changes_corr_all8_s0.6', np.nan), '.2f')} | {_fmt(r.get('changes_corr_T_s1', np.nan), '.2f')} |")
    L += ['', '### The searches', '',
          '| input | final step size, median over odors (natural log; bound ±0.69) | factors at a bound | variance of the '
          'log factors explained by a per-glomerulus term | the same for permuted factors, mean / 95th percentile |',
          '|---|---:|---:|---:|---:|']
    for r in se.to_dict('records'):
        L.append(f"| {FITTED[r['input']][4]} | {_fmt(r['sigma_final_median'], '.2f')} | {_fmt(100 * r['share_at_bound'], '.0f')} % | "
                 f"{_fmt(r['var_glomerulus'], '.2f')} | {_fmt(r['var_glomerulus_perm_mean'], '.2f')} / {_fmt(r['var_glomerulus_perm_p95'], '.2f')} |")
    L += ['', f"### Against the ORN flies (E85's part A4, cited from results/{E85_A4[0]}/{E85_A4[1]}, not recomputed)", '',
          '| input | entries within one SD of the flies\' spread (log2) | within two | rms log2 change | share of the squared '
          'size inside the nine flies\' span (permuted mean) | odors above the permuted 95th percentile |',
          '|---|---:|---:|---:|---:|---:|']
    e85 = e85.set_index('input')
    for name in FITTED:
        k = E85_A4_NAMES[name]
        if k not in e85.index:
            continue
        r = e85.loc[k]
        L.append(f"| {FITTED[name][4]} | {_fmt(100 * r.within1_log2, '.0f')} % | {_fmt(100 * r.within2_log2, '.1f')} % | "
                 f"{_fmt(r.rms_log2_change, '.2f')} | {_fmt(r.span_share, '.2f')} ({_fmt(r.span_share_perm_mean, '.2f')}) | "
                 f"{int(r.n_odors_above_perm_p95)} of {int(r.n_odors)} |")
    _section(out, '## Part A: the four fitted inputs, by the pairs they move (no simulations)', L + [''])


# ---------------------------------------------------------------------------
# identity

def identity(out: Path, wirings: Sequence[str] = T_WIRINGS) -> dict:
    """The identity checks: every wiring of `wirings` reproduces its saved fit with the lab's input (E83's
    `build_simulator` and `_check_saved`); E83's part A2 point A2_halfA_s0.6_P3 (half A's input at 0.6 on P3) is
    reproduced exactly by this experiment's reading (`held_reading`, its own breadth cache): the input to the last digit
    (E83's inputs/A2_halfA_s0.6.csv, read with round-trip float parsing), the breadth (sha256), both factors, every spike
    count and response. Writes
    identity.json; raises on a mismatch."""
    res = dict(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
               wirings={}, a2_point={})
    for w in wirings:
        t0 = time.time()
        sim, info = E83.build_simulator(w)
        E83._check_saved(w, sim, info['params'], E83.saved_fit(w)['saved'])
        del sim, info
        gc.collect()
        res['wirings'][w] = dict(fit=f"{E83.saved_fit(w)['folder']}:{E83.saved_fit(w)['tag']}", reproduced=True,
                                 seconds=time.time() - t0)
        print(f'{w}: the lab input reproduces {res["wirings"][w]["fit"]} exactly', flush=True)
    e83 = H.results_dir(E83.EXPERIMENT)
    orn_in = E83.input_at(candidate_factors('A'), SIZE)
    lay = H.abbrev_columns(orn_in)
    saved_in = pd.read_csv(e83 / 'inputs' / 'A2_halfA_s0.6.csv', index_col=0, float_precision='round_trip')
    if not (list(saved_in.index) == list(lay.index) and list(saved_in.columns) == list(lay.columns)
            and np.array_equal(saved_in.to_numpy(dtype=float), lay.to_numpy(dtype=float))):
        raise RuntimeError("half A's input at 0.6 differs from E83's inputs/A2_halfA_s0.6.csv")
    res['a2_input_equal'] = True
    q = pd.read_csv(e83 / 'runs_A2.csv', float_precision='round_trip').set_index('tag').loc['A2_halfA_s0.6_P3']
    t0 = time.time()
    r = held_reading('P3', orn_in, E83.input_tag_for(orn_in, 'E86_identity_A2_halfA_s0.6'))
    o = {'on': r['on'], 'off': r['off']}
    ns, nr = E83._same(o, E83.load_run('A2', 'A2_halfA_s0.6_P3'))
    a2 = dict(tag='A2_halfA_s0.6_P3', count_diffs=ns, response_diffs=nr, thr_factor=r['thr_factor'],
              thr_factor_e83=float(q.thr_factor), apl_factor=r['apl_factor'], apl_factor_e83=float(q.apl_factor),
              breadth_sha256=r['breadth_sha256'], breadth_sha256_e83=q.breadth_sha256, n_sims=r['n_sims'],
              seconds=time.time() - t0)
    a2['reproduced'] = bool(ns == 0 and nr == 0 and a2['thr_factor'] == a2['thr_factor_e83']
                            and a2['apl_factor'] == a2['apl_factor_e83'] and a2['breadth_sha256'] == a2['breadth_sha256_e83'])
    res['a2_point'] = a2
    (out / 'identity.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(json.dumps(res, indent=1, default=float), flush=True)
    if not a2['reproduced']:
        raise RuntimeError(f"E83's point A2_halfA_s0.6_P3 not reproduced: {a2}")
    return res


# ---------------------------------------------------------------------------
# parts B and C-read: held readings of inputs on T

def plan_b(part: str) -> List[dict]:
    """The points of part B1, B2 or B3: [dict(part, change, group, kind, seed, half, unseen, size, wiring, tag)], one
    input per direction (and seed), each on the four wirings of T."""
    if part not in ('B1', 'B2', 'B3'):
        raise ValueError(part)
    plan = []
    for X in ('A', 'B'):
        specs = ([('perm', s) for s in (B1_SEEDS if part == 'B1' else B3_SEEDS)] if part != 'B2' else [('inverted', None)])
        for kind, seed in specs:
            ch = f'{part}_perm{seed}_{X}' if kind == 'perm' else f'B2_inv_{X}'
            size = SIZE if kind == 'perm' else INVERT_SIZE
            plan += [dict(part=part, change=ch, group=part, kind=kind, seed=seed, half=X, unseen=UNSEEN[X], size=size,
                          wiring=w, tag=f'{ch}_{w}') for w in T_WIRINGS]
    return plan


def plan_c_read(out: Path) -> List[dict]:
    """The points of C-read: each finished search's input at 0.6 on T (searches not finished are left out)."""
    plan = []
    for name, cfg in SEARCHES.items():
        try:
            search_factors(out, name)
        except (FileNotFoundError, RuntimeError) as e:
            print(f'{name}: not read ({e})', flush=True)
            continue
        plan += [dict(part='Cread', change=name, group=name.split('_')[0], kind='search', seed=GENERATOR_SEEDS[name],
                      half=cfg['half'], unseen=UNSEEN[cfg['half']], size=SIZE, wiring=w, tag=f'{name}_{w}')
                 for w in T_WIRINGS]
    return plan


def point_input(out: Path, pt: Mapping) -> pd.DataFrame:
    """The ORN input of a planned point."""
    if pt['kind'] == 'search':
        return E83.input_at(search_factors(out, pt['change']), pt['size'])
    F, size = control_factors(pt['kind'], pt['half'], pt['seed'])
    if size != pt['size']:
        raise RuntimeError(f"{pt['change']}: size {size} planned as {pt['size']}")
    return E83.input_at(F, size)


def run_points(out: Path, part: str, plan: Sequence[dict]) -> pd.DataFrame:
    """Held readings of a plan's points (resumable: points in runs_<part>.csv are skipped; outputs under
    outputs/E86_reduced_input_candidate/<part>/, provenance per point in provenance_<part>.json, inputs in
    results/.../inputs/<change>.csv). Each wiring's simulator is first checked against its saved fit."""
    rp, pp = out / f'runs_{part}.csv', out / f'provenance_{part}.json'
    done = set(E83._read_rows(rp).get('tag', pd.Series(dtype=str)))
    todo = [p for p in plan if p['tag'] not in done]
    for w in [w for w in E83.ALL_WIRINGS if any(p['wiring'] == w for p in todo)]:
        sim, info = E83.build_simulator(w)
        E83._check_saved(w, sim, info['params'], E83.saved_fit(w)['saved'])
        del sim, info
        gc.collect()
        print(f"{w}: the simulator reproduces {E83.saved_fit(w)['folder']}:{E83.saved_fit(w)['tag']}", flush=True)
        for pt in [p for p in todo if p['wiring'] == w]:
            orn_in = point_input(out, pt)
            path = E83._write_input(out, pt['change'], orn_in)
            res = held_reading(w, orn_in, E83.input_tag_for(orn_in, f"E86_{pt['change']}"))
            o = {'on': res['on'], 'off': res['off']}
            save_run(part, pt['tag'], o)
            E83._append_row(rp, E83._run_row(part, pt['change'], w, pt['tag'], res, o, group=pt['group'], kind=pt['kind'],
                                             control_seed=pt['seed'], half=pt['half'], unseen=pt['unseen'],
                                             size=pt['size'], input_tag=res['input_tag'], orn_input=_rel(path)))
            _provenance_add(pp, part, pt['tag'], H.model_provenance(orn_deltas=orn_in))
            t = targets()[f"half{pt['unseen']}"]
            d0, d1 = E83.pair_distances(E83.saved_fit(w)['saved'], 'cal'), E83.pair_distances(o, 'cal')
            print(f"{pt['tag']}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f}; pattern under the readout "
                  f"against the unseen half {pt['unseen']} {_rms(E83.pattern_error(d0, t)):.4f} -> "
                  f"{_rms(E83.pattern_error(d1, t)):.4f}; {res['seconds']:.0f}s", flush=True)
            del res, o
            gc.collect()
    return E83._read_rows(rp)


# ---------------------------------------------------------------------------
# part C: the searches

def _search_dir() -> Path:
    return _local('C')


def _ulps(a: float, b: float) -> float:
    """|a - b| in units in the last place of the larger."""
    return float(abs(a - b) / np.spacing(max(abs(a), abs(b)))) if a != b else 0.0


def check_search(out: Path, simulate: bool = True) -> dict:
    """Part C's check: E77's `input_oracle` with its defaults' split and generator seed (first split, seed 0, half A),
    one generation of eb, written under outputs/E86_reduced_input_candidate/C/check/ (tag check_A), against the first
    row of E77's oracle_history_A.csv in every column but the time (both read with round-trip float parsing). 'exact':
    every column equal to the last digit; 'reproduced': the generator's path and the count-derived columns (success,
    sigma, the response fractions) equal to the last digit and the pattern-error columns (CHECK_FLOAT_COLUMNS) within a
    relative CHECK_REL_TOL (their units in the last place are reported). Writes check_search.json; raises unless
    reproduced. `simulate=False` compares the history the last run wrote, without running it again."""
    d = _search_dir() / 'check'
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if simulate:
        E77.input_oracle(d, CHECK['half'], odors=CHECK['odors'], gens=CHECK['gens'], seed=CHECK['seed'],
                         split_seed=CHECK['split_seed'], tag='check_A')
    hpath = E77.oracle_paths(d, 'A', 'check_A')['history']
    got = pd.read_csv(hpath, float_precision='round_trip').iloc[0]
    ref = pd.read_csv(H.results_dir(E77.EXPERIMENT) / 'oracle_history_A.csv', float_precision='round_trip').iloc[0]
    cols = {c: dict(e86=got[c], e77=ref[c], equal=bool(got[c] == ref[c])) for c in CHECK_COLUMNS}
    for c in CHECK_FLOAT_COLUMNS:
        a, b = float(got[c]), float(ref[c])
        cols[c].update(ulps=_ulps(a, b), rel_diff=abs(a - b) / max(abs(a), abs(b)))
    exact = all(v['equal'] for v in cols.values())
    within = all(v['equal'] for c, v in cols.items() if c not in CHECK_FLOAT_COLUMNS) and \
        all(cols[c]['rel_diff'] <= CHECK_REL_TOL for c in CHECK_FLOAT_COLUMNS)
    res = dict(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
               check=dict(CHECK, odors=list(CHECK['odors'])), history=str(hpath), simulated=simulate,
               seconds=time.time() - t0, columns=cols, exact=exact,
               max_ulps=max(cols[c]['ulps'] for c in CHECK_FLOAT_COLUMNS),
               max_rel_diff=max(cols[c]['rel_diff'] for c in CHECK_FLOAT_COLUMNS), rel_tol=CHECK_REL_TOL,
               reproduced=bool(exact or within))
    (out / 'check_search.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
    print(json.dumps(res, indent=1, default=str), flush=True)
    if not res['reproduced']:
        raise RuntimeError(f"the first generation of eb differs from E77's: {[c for c, v in cols.items() if not v['equal']]}")
    print("the first generation of eb reproduces E77's " + ('exactly' if exact else
          f"every column but the pattern-error ones exactly, those to {res['max_ulps']:g} ulps (relative "
          f"{res['max_rel_diff']:.1g})"), flush=True)
    return res


def search_call(name: str) -> dict:
    """The `input_oracle` arguments of a search: half (the half of `recording_halves(split_seed)`), seed, split_seed, tag;
    checks that the generator seed is GENERATOR_SEEDS[name]."""
    cfg = SEARCHES[name]
    half = SPLITS[cfg['split_seed']][cfg['half']]
    if E77.oracle_generator_seed(cfg['seed'], half) != GENERATOR_SEEDS[name]:
        raise RuntimeError(f'{name}: generator seed {E77.oracle_generator_seed(cfg["seed"], half)}, planned {GENERATOR_SEEDS[name]}')
    return dict(half_name=half, seed=cfg['seed'], split_seed=cfg['split_seed'], tag=name)


def search(out: Path, name: str, resume: bool = True) -> pd.DataFrame:
    """One of part C's searches (SEARCHES): E77's part E search as it was (`input_oracle`: C2 at its fitted parameters
    and breadth on the real wiring, the ten odors in E77's order, factors in [0.5, 2], (1 + 32), eight generations per
    odor, K held at C2's, the penalty outside +-15 % of C2's response fractions) against the named half, its files tagged
    with `name`, checkpointed every generation under outputs/E86_reduced_input_candidate/C/checkpoint_<name>.pkl and
    resumed from it unless `resume` is False. Writes provenance_C_<name>.json first."""
    call = search_call(name)
    ck = _search_dir() / f'checkpoint_{name}.pkl'
    if not resume:
        for p in [ck] + list(E77.oracle_paths(out, call['half_name'], name).values()):
            if p.exists():
                p.unlink()
    hv = E77.recording_halves(call['split_seed'])
    prov = dict(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), search=name,
                half=SEARCHES[name]['half'], half_in_split=call['half_name'], split_seed=call['split_seed'],
                seed_argument=call['seed'], generator_seed=GENERATOR_SEEDS[name], checkpoint=str(ck),
                halves={h: {g: [list(map(str, r)) for r in v] for g, v in hv[h].items()} for h in hv},
                settings=dict(odors=list(E77.E_ODORS), pop=E77.E_POP, generations=E77.E_GENS, bounds=E77.E_BOUNDS,
                              rate_tol=E77.E_RATE_TOL, penalty=E77.E_PENALTY, sigma0=E77.E_SIGMA0,
                              sigma_min=E77.E_SIGMA_MIN, sigma_max=E77.E_SIGMA_MAX),
                fit=dict(folder=E77.C2[0], tag=E77.C2[1]), **H.code_commits(),
                model=H.model_provenance(orn_deltas=H.apl16_orn_deltas()), approx_conc_note=H.APPROX_CONC_NOTE)
    (out / f'provenance_C_{name}.json').write_text(json.dumps(prov, indent=1, default=str) + '\n')
    return E77.input_oracle(out, call['half_name'], seed=call['seed'], split_seed=call['split_seed'], tag=name,
                            checkpoint=ck)


# ---------------------------------------------------------------------------
# scoring

def _change_meta() -> Dict[str, dict]:
    meta = {f'cand_{X}': dict(group='candidate', kind='candidate', seed=E77.oracle_generator_seed(0, X), half=X,
                              unseen=UNSEEN[X], size=SIZE, part='A2') for X in ('A', 'B')}
    for part in ('B1', 'B2', 'B3'):
        for p in plan_b(part):
            meta[p['change']] = {k: p[k] for k in ('group', 'kind', 'seed', 'half', 'unseen', 'size', 'part')}
    for name, cfg in SEARCHES.items():
        meta[name] = dict(group=name.split('_')[0], kind='search', seed=GENERATOR_SEEDS[name], half=cfg['half'],
                          unseen=UNSEEN[cfg['half']], size=SIZE, part='Cread')
    return meta


def collect_runs(out: Path) -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """Every run so far: {(change, wiring): outputs} (the candidate's: E83's archived part A2 runs on T, changes cand_A
    and cand_B) and the runs table (the candidate's rows from E83's runs_A2.csv, then this experiment's runs_<part>.csv)."""
    runs, rows = {}, []
    a2 = pd.read_csv(H.results_dir(E83.EXPERIMENT) / 'runs_A2.csv', float_precision='round_trip')
    for X in ('A', 'B'):
        for w in T_WIRINGS:
            tag = f'A2_half{X}_s{SIZE:g}_{w}'
            o = E83.load_run('A2', tag)
            if o is None:
                raise RuntimeError(f"E83's run {tag} not found")
            runs[(f'cand_{X}', w)] = o
            r = a2[a2.tag == tag].iloc[0].to_dict()
            rows.append(dict(r, part='A2', change=f'cand_{X}', group='candidate', kind='candidate', unseen=UNSEEN[X],
                             source=f"E83 part A2 ({tag}, archived, not rerun)"))
    for part in PARTS:
        rt = E83._read_rows(out / f'runs_{part}.csv')
        for r in rt.to_dict('records'):
            o = load_run(part, r['tag'])
            if o is None:
                print(f"{r['tag']}: outputs missing", flush=True)
                continue
            runs[(r['change'], r['wiring'])] = o
            rows.append(r)
    return runs, pd.DataFrame(rows)


def two_by_two(readings: pd.DataFrame, meta: Mapping[str, dict], units: Sequence[str] = ('T',) + T_WIRINGS
               ) -> pd.DataFrame:
    """Per input (change) on `units` (T and its single wirings by default; part E: every population and wiring), under
    the calibrated readout: the pattern part unchanged and changed against the half searched on, the half never seen,
    all recordings and the six typical ones, the gain and the reduction of the squared pattern error (E83's `two_by_two`
    columns, with this experiment's halves)."""
    r = readings[readings.readout == 'cal']
    base = r[r.change == 'none'].set_index('unit')
    rows = []
    for ch, m in meta.items():
        g = r[r.change == ch].set_index('unit')
        if g.empty:
            continue
        for unit in units:
            if unit not in g.index:
                continue
            for which, t in (('searched', f"half{m['half']}"), ('unseen', f"half{m['unseen']}"), ('all', 'z8'), ('six', 'z6')):
                u0, u1 = float(base.loc[unit, f'{t}_pattern']), float(g.loc[unit, f'{t}_pattern'])
                rows.append(dict(change=ch, group=m['group'], searched_on=m['half'], size=m['size'], unit=unit,
                                 kind=g.loc[unit, 'kind'], against=which, target=t, unchanged=u0, changed=u1, gain=u0 - u1,
                                 sq_reduction=1.0 - (u1 ** 2) / (u0 ** 2)))
    return pd.DataFrame(rows, columns=['change', 'group', 'searched_on', 'size', 'unit', 'kind', 'against', 'target',
                                       'unchanged', 'changed', 'gain', 'sq_reduction'])


def inputs_table(readings: pd.DataFrame, tb: pd.DataFrame, rt: pd.DataFrame, meta: Mapping[str, dict]) -> pd.DataFrame:
    """One row per input read on all four wirings of T: g (T against the unseen half) and the gate's statistics: the
    single wirings lower (of 4), T's raw-count loss change against all recordings, T's mean APL effect under the
    readout, the reductions of T's squared pattern error against all recordings and the six typical ones, the points
    flagged by the four bands (the rule of 16:10)."""
    rows = []
    for ch, m in meta.items():
        t = tb[(tb.change == ch) & (tb.unit == 'T')].set_index('against')
        if t.empty:
            continue
        singles = tb[(tb.change == ch) & tb.unit.isin(T_WIRINGS) & (tb.against == 'unseen')]
        rawT = readings[(readings.change == ch) & (readings.unit == 'T') & (readings.readout == 'raw')]
        raw0 = readings[(readings.change == 'none') & (readings.unit == 'T') & (readings.readout == 'raw')]
        calT = readings[(readings.change == ch) & (readings.unit == 'T') & (readings.readout == 'cal')]
        rr = rt[rt.change == ch]
        flags = (rr.realism_flag.fillna('').astype(str).str.strip() != '') if 'realism_flag' in rr.columns else pd.Series(dtype=bool)
        rows.append(dict(change=ch, group=m['group'], kind=m['kind'], seed=m['seed'], searched_on=m['half'],
                         unseen=m['unseen'], split=split_of(m['half']), size=m['size'],
                         g=float(t.loc['unseen', 'gain']), gain_searched=float(t.loc['searched', 'gain']),
                         pattern_unseen_unchanged=float(t.loc['unseen', 'unchanged']),
                         pattern_unseen_changed=float(t.loc['unseen', 'changed']),
                         n_singles_lower=int((singles.gain > 0).sum()), n_singles=int(len(singles)),
                         raw_loss_change=float(rawT.z8_loss.iloc[0] - raw0.z8_loss.iloc[0]),
                         mean_effect_cal=float(calT.mean_effect.iloc[0]),
                         sq_reduction_all=float(t.loc['all', 'sq_reduction']), sq_reduction_six=float(t.loc['six', 'sq_reduction']),
                         n_flagged=int(flags.sum()), n_points=int(len(rr))))
    return pd.DataFrame(rows)


def _archive_tables(out: Path, rt: pd.DataFrame, of_parts: Sequence[str] = PARTS) -> None:
    """model_pair_distances_long.csv rows (raw spike counts, the harness's primary model RDM) and variants.csv rows
    (absolute parameters) of this experiment's own runs (those of `of_parts`), for kcq/apl_effect.py and the model-output
    archive."""
    own = rt[rt.part.isin(of_parts)] if len(rt) else rt
    parts, vrows = [], []
    for r in own.to_dict('records'):
        o = load_run(r['part'], r['tag'])
        if o is None:
            continue
        parts.append(pd.concat({(r['tag'], 'on'): H.model_rdms(*o['on'])[H.PRIMARY_MODEL_RDM],
                                (r['tag'], 'acute_off'): H.model_rdms(*o['off'])[H.PRIMARY_MODEL_RDM]},
                               names=['tag', 'condition', 'odor_a', 'odor_b']).rename('distance').reset_index())
        vrows.append(dict(tag=r['tag'], variant=r['tag'], part=r['part'], change=r['change'], wiring=r['wiring'],
                          base='bouton', k=E83.BREADTH_K, fit_folder=r['fit_folder'], fit_tag=r['fit_tag'],
                          kc_kc_inh=r['wKCKC_scale'], **{k: r[k] for k in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale',
                                                                          'wAPLPN_scale', 'wPNAPL_scale', 'thr_factor',
                                                                          'apl_factor', 'rate_on', 'rate_off')},
                          orn_input=r.get('orn_input'), subdir=r['part']))
    if parts:
        new = pd.concat(parts, ignore_index=True)
        E83._replace_rows(out / 'model_pair_distances_long.csv', new, 'tag', set(new.tag))
    if vrows:
        v = pd.DataFrame(vrows)
        E83._replace_rows(out / 'variants.csv', v, 'tag', set(v.tag))


def score(out: Path) -> pd.DataFrame:
    """Scores every run so far (`collect_runs`) with E83's `score_runs` against this experiment's targets: readings.csv,
    pairs.csv, two_by_two.csv, inputs.csv, runs.csv and the archive tables. The candidate's g must equal E83's
    two_by_two_A2.csv (raises otherwise)."""
    runs, rt = collect_runs(out)
    meta = _change_meta()
    # (score_runs' rows carry the metadata; its own 'kind' is the unit's)
    sm = {k: dict({a: b for a, b in v.items() if a != 'kind'}, input_kind=v['kind']) for k, v in meta.items()}
    readings, pairs = E83.score_runs(runs, sm, tg=targets())
    readings.to_csv(out / 'readings.csv', index=False)
    pairs.to_csv(out / 'pairs.csv', index=False)
    tb = two_by_two(readings, meta)
    tb.to_csv(out / 'two_by_two.csv', index=False)
    it = inputs_table(readings, tb, rt, meta)
    e83 = pd.read_csv(H.results_dir(E83.EXPERIMENT) / 'two_by_two_A2.csv', float_precision='round_trip')
    for X in ('A', 'B'):
        ref = e83[(e83.searched_on == X) & np.isclose(e83['size'], SIZE) & (e83.unit == 'T') & (e83.against == 'unseen')]
        got = float(it[it.change == f'cand_{X}'].g.iloc[0])
        if abs(got - float(ref.gain.iloc[0])) > 1e-12:
            raise RuntimeError(f"the candidate's g (half {X}) is {got}, E83's {float(ref.gain.iloc[0])}")
    it.to_csv(out / 'inputs.csv', index=False)
    rt.to_csv(out / 'runs.csv', index=False)
    _archive_tables(out, rt)
    print(it[['change', 'group', 'searched_on', 'unseen', 'g', 'n_singles_lower', 'raw_loss_change', 'mean_effect_cal',
              'sq_reduction_all', 'sq_reduction_six', 'n_flagged']].round(4).to_string(index=False), flush=True)
    return it


# ---------------------------------------------------------------------------
# the verdict (part D)

def b3_needed(inputs: pd.DataFrame) -> bool:
    """Whether part B3 runs: one of B1's permuted inputs has g above 0."""
    b1 = inputs[inputs.group == 'B1']
    return bool((b1.g > 0).any())


def verdict(inputs: pd.DataFrame, flies_mean: float) -> dict:
    """The entry's verdict (part D), fixed before any run, from a table with one row per input (inputs.csv's columns:
    change, group, searched_on, g, n_singles_lower, n_singles, raw_loss_change, mean_effect_cal, sq_reduction_all,
    sq_reduction_six, n_flagged; groups 'candidate' (two rows, searched on A and B), 'B1' (four), 'B2' (two), 'B3' (four,
    when run), 'C1' (two, A and B), 'C2' (two, C and D)).

    - The controls hold if (v) every B1 input has g below CONTROL_SHARE of the candidate's in its direction; if B3 ran,
      in each direction the candidate's g is above all four permuted inputs' (B1 and B3) and their mean is below
      B3_MEAN_SHARE of it; (vi) both inverted inputs (B2) have g of 0 or less.
    - The replication: the mean g of the four new inputs (C1, C2) and how many are positive.
    - E83's gate on the mean of the two splits (the candidate's two inputs and C2's two): (i) mean g of GATE_GAIN or
      more; (ii) at least GATE_SINGLES of the 16 single-wiring cells of T lower; (iii) every point inside the four bands,
      T's raw-count loss (the mean over the four inputs of its change) not above its unchanged value by more than
      GATE_RAW_SLACK, T's mean APL effect under the readout (the mean over the four) within GATE_MEAN_SLACK of the
      flies' (`flies_mean`); (iv) T's squared pattern error against all recordings GATE_SQ_ALL lower and GATE_SQ_SIX
      against the six typical recordings (means over the four).
    - 'passes': the controls hold, the gate holds, and the four new inputs have a mean g of PASS_MEAN or more with at
      least MIN_POSITIVE positive. 'carries': the controls hold and the four new inputs have a mean g of CARRY_MEAN or
      more with at least MIN_POSITIVE positive, without a pass. 'falls': anything else. 'incomplete' while a required
      input is missing (B3 is required when `b3_needed`).

    Returns a dict with 'outcome', 'missing', 'controls', 'replication', 'gate' and 'replicate_kinds' (C1's and C2's
    means and counts, for outcome (e))."""
    grp = {g: inputs[inputs.group == g] for g in GROUPS}
    missing = []
    need = {'candidate': 2, 'B1': 4, 'B2': 2, 'C1': 2, 'C2': 2}
    for g, n in need.items():
        if len(grp[g]) < n:
            missing.append(f'{g}: {len(grp[g])} of {n} inputs')
    cand = {r.searched_on: float(r.g) for r in grp['candidate'].itertuples()}
    res = dict(outcome=None, missing=missing, flies_mean=float(flies_mean), candidate_g=cand)

    # the controls
    ctl = dict(v=None, b3_ran=bool(len(grp['B3'])), b3_needed=b3_needed(inputs), b3=None, vi=None, hold=None)
    if len(grp['B1']) and set(cand) >= {'A', 'B'}:
        rows = [dict(change=r.change, searched_on=r.searched_on, g=float(r.g), threshold=CONTROL_SHARE * cand[r.searched_on],
                     holds=bool(r.g < CONTROL_SHARE * cand[r.searched_on])) for r in grp['B1'].itertuples()]
        ctl['v_rows'] = rows
        ctl['v'] = bool(len(rows) >= need['B1'] and all(r['holds'] for r in rows))
    if ctl['b3_needed'] and len(grp['B3']) < 4:
        missing.append(f"B3: {len(grp['B3'])} of 4 inputs (a B1 input has g above 0)")
    if len(grp['B3']) and set(cand) >= {'A', 'B'}:
        b3 = {}
        for X in ('A', 'B'):
            perm = pd.concat([grp['B1'], grp['B3']])
            perm = perm[perm.searched_on == X]
            b3[X] = dict(n=int(len(perm)), max_g=float(perm.g.max()), mean_g=float(perm.g.mean()), candidate_g=cand[X],
                         above_all=bool(cand[X] > perm.g.max()),
                         mean_below=bool(perm.g.mean() < B3_MEAN_SHARE * cand[X]))
            b3[X]['holds'] = bool(b3[X]['n'] >= 4 and b3[X]['above_all'] and b3[X]['mean_below'])
        ctl['b3_rows'] = b3
        ctl['b3'] = bool(all(v['holds'] for v in b3.values()))
    if len(grp['B2']):
        ctl['vi_rows'] = [dict(change=r.change, searched_on=r.searched_on, g=float(r.g), holds=bool(r.g <= 0))
                          for r in grp['B2'].itertuples()]
        ctl['vi'] = bool(len(grp['B2']) >= need['B2'] and (grp['B2'].g <= 0).all())
    if ctl['v'] is not None and ctl['vi'] is not None and (ctl['b3'] is not None or not ctl['b3_needed']):
        ctl['hold'] = bool(ctl['v'] and ctl['vi'] and (ctl['b3'] if ctl['b3_ran'] else True))
    res['controls'] = ctl

    # the replication
    rep_rows = pd.concat([grp['C1'], grp['C2']])
    rep = dict(n=int(len(rep_rows)), mean_g=float(rep_rows.g.mean()) if len(rep_rows) else None,
               n_positive=int((rep_rows.g > 0).sum()), rows=[dict(change=r.change, searched_on=r.searched_on, g=float(r.g))
                                                             for r in rep_rows.itertuples()])
    res['replication'] = rep
    res['replicate_kinds'] = {g: dict(n=int(len(grp[g])), mean_g=float(grp[g].g.mean()) if len(grp[g]) else None,
                                      n_positive=int((grp[g].g > 0).sum()),
                                      carries=bool(len(grp[g]) == 2 and grp[g].g.mean() >= CARRY_MEAN and (grp[g].g > 0).all()))
                              for g in ('C1', 'C2')}

    # E83's gate on the mean of the two splits
    gate_rows = pd.concat([grp['candidate'], grp['C2']])
    gate = dict(inputs=list(gate_rows.change), n=int(len(gate_rows)))
    if len(gate_rows):
        gate.update(mean_g=float(gate_rows.g.mean()), mean_g_split0=float(grp['candidate'].g.mean()) if len(grp['candidate']) else None,
                    mean_g_split1=float(grp['C2'].g.mean()) if len(grp['C2']) else None,
                    n_singles_lower=int(gate_rows.n_singles_lower.sum()), n_singles=int(gate_rows.n_singles.sum()),
                    n_flagged=int(gate_rows.n_flagged.sum()), raw_loss_change_mean=float(gate_rows.raw_loss_change.mean()),
                    raw_loss_change_max=float(gate_rows.raw_loss_change.max()),
                    mean_effect_cal_mean=float(gate_rows.mean_effect_cal.mean()),
                    mean_effect_cal_farthest=float(gate_rows.mean_effect_cal.iloc[int(np.argmax(np.abs(gate_rows.mean_effect_cal - flies_mean)))]),
                    sq_reduction_all_mean=float(gate_rows.sq_reduction_all.mean()),
                    sq_reduction_six_mean=float(gate_rows.sq_reduction_six.mean()))
        gate['i'] = bool(gate['mean_g'] >= GATE_GAIN)
        gate['ii'] = bool(gate['n_singles_lower'] >= GATE_SINGLES)
        gate['iii'] = bool(gate['n_flagged'] == 0 and gate['raw_loss_change_mean'] <= GATE_RAW_SLACK
                           and abs(gate['mean_effect_cal_mean'] - flies_mean) <= GATE_MEAN_SLACK)
        gate['iv'] = bool(gate['sq_reduction_all_mean'] >= GATE_SQ_ALL and gate['sq_reduction_six_mean'] >= GATE_SQ_SIX)
        gate['holds'] = bool(gate['n'] >= 4 and gate['i'] and gate['ii'] and gate['iii'] and gate['iv'])
    res['gate'] = gate

    if missing or ctl['hold'] is None:
        res['outcome'] = 'incomplete'
        return res
    rep_ok = lambda m: bool(rep['n'] >= 4 and rep['mean_g'] >= m and rep['n_positive'] >= MIN_POSITIVE)
    if ctl['hold'] and gate['holds'] and rep_ok(PASS_MEAN):
        res['outcome'] = 'passes'
    elif ctl['hold'] and rep_ok(CARRY_MEAN):
        res['outcome'] = 'carries'
    else:
        res['outcome'] = 'falls'
    return res


def write_verdict(out: Path) -> dict:
    """The verdict from inputs.csv: gate.csv (the gate's four inputs, then a row 'gate' with its conditions), replication.csv
    (the four new inputs, then a row 'mean'), null.csv (the controls with their thresholds), verdict.json and the
    verdict's section of comparison_table.md."""
    it = pd.read_csv(out / 'inputs.csv', float_precision='round_trip')
    flies_mean = float(targets()['z8'].mean())
    v = verdict(it, flies_mean)
    v.update(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
             thresholds=dict(control_share=CONTROL_SHARE, b3_mean_share=B3_MEAN_SHARE, gate_gain=GATE_GAIN,
                             gate_singles=GATE_SINGLES, gate_raw_slack=GATE_RAW_SLACK, gate_mean_slack=GATE_MEAN_SLACK,
                             gate_sq_all=GATE_SQ_ALL, gate_sq_six=GATE_SQ_SIX, pass_mean=PASS_MEAN, carry_mean=CARRY_MEAN,
                             min_positive=MIN_POSITIVE))
    cols = ['change', 'group', 'searched_on', 'unseen', 'split', 'g', 'n_singles_lower', 'n_singles', 'raw_loss_change',
            'mean_effect_cal', 'sq_reduction_all', 'sq_reduction_six', 'n_flagged', 'n_points']
    gr = it[it.group.isin(['candidate', 'C2'])][cols]
    g = v['gate']
    summ = dict(change='gate', group='mean of the two splits', g=g.get('mean_g'), n_singles_lower=g.get('n_singles_lower'),
                n_singles=g.get('n_singles'), raw_loss_change=g.get('raw_loss_change_mean'),
                mean_effect_cal=g.get('mean_effect_cal_mean'), sq_reduction_all=g.get('sq_reduction_all_mean'),
                sq_reduction_six=g.get('sq_reduction_six_mean'), n_flagged=g.get('n_flagged'),
                i=g.get('i'), ii=g.get('ii'), iii=g.get('iii'), iv=g.get('iv'), holds=g.get('holds'))
    pd.concat([gr, pd.DataFrame([summ])], ignore_index=True).to_csv(out / 'gate.csv', index=False)
    rp = it[it.group.isin(['C1', 'C2'])][cols]
    rs = dict(change='mean', group='C1 and C2', g=v['replication']['mean_g'], n_positive=v['replication']['n_positive'],
              n=v['replication']['n'])
    pd.concat([rp, pd.DataFrame([rs])], ignore_index=True).to_csv(out / 'replication.csv', index=False)
    nl = it[it.group.isin(['B1', 'B2', 'B3'])][cols + ['seed']].copy()
    cand = v['candidate_g']
    nl['candidate_g'] = nl.searched_on.map(cand)
    # B1: below half the candidate's g ((v)); B3: below the candidate's (its mean condition is in verdict.json);
    # B2: 0 or less ((vi))
    share = {'B1': CONTROL_SHARE, 'B3': 1.0, 'B2': 0.0}
    nl['threshold'] = [share[gr_] * cand.get(x, np.nan) for gr_, x in zip(nl.group, nl.searched_on)]
    nl['holds'] = [bool(gv <= 0) if gr_ == 'B2' else bool(gv < th) for gr_, gv, th in zip(nl.group, nl.g, nl.threshold)]
    nl.to_csv(out / 'null.csv', index=False)
    (out / 'verdict.json').write_text(json.dumps(v, indent=1, default=float) + '\n')
    _comparison_table_d(out, it, v)
    print(json.dumps({k: v[k] for k in ('outcome', 'missing', 'candidate_g', 'replication', 'replicate_kinds')}, indent=1,
                     default=float), flush=True)
    return v


GROUP_LABELS = {'candidate': "the candidate (E77's input at 0.6; E83's part A2 runs)", 'B1': 'labels permuted (B1)',
                'B2': 'inverted (B2)', 'B3': 'labels permuted (B3)', 'C1': 'second draw on the same halves (C1)',
                'C2': 'second split of the flies (C2)'}


def _comparison_table_d(out: Path, it: pd.DataFrame, v: dict) -> None:
    L = ['', f"Verdict: **{v['outcome']}**" + (f" (missing: {'; '.join(v['missing'])})" if v['missing'] else '') + '.', '',
         'On T (P3, P4, M2, M3) at 0.6, each wiring at its own fit, both response rates re-held; the pattern part of the '
         'APL-effect loss under the calibrated readout against the half each input\'s search never saw (lower is better); '
         'g = unchanged minus changed (positive is better).', '',
         '| input | searched on | scored on | g | T\'s wirings lower | raw-count loss change | mean APL effect | '
         'squared pattern error lower: all recordings / six typical | points flagged |',
         '|---|---|---|---:|---:|---:|---:|---:|---:|']
    for r in it.to_dict('records'):
        L.append(f"| {GROUP_LABELS.get(r['group'], r['group'])}: {r['change']} | {r['searched_on']} | {r['unseen']} | "
                 f"{_fmt(r['g'], '+.4f')} | {r['n_singles_lower']} of {r['n_singles']} | {_fmt(r['raw_loss_change'], '+.4f')} | "
                 f"{_fmt(r['mean_effect_cal'], '+.3f')} | {_fmt(100 * r['sq_reduction_all'], '.1f')} % / "
                 f"{_fmt(100 * r['sq_reduction_six'], '.1f')} % | {r['n_flagged']} of {r['n_points']} |")
    L.append(f"| KC data | | | | | | {_fmt(v['flies_mean'], '+.3f')} | | |")
    c, rep, g = v['controls'], v['replication'], v['gate']
    yn = lambda x: '–' if x is None else ('holds' if x else 'fails')
    L += ['', f"- Controls: (v) {yn(c.get('v'))}; B3 {'ran' if c.get('b3_ran') else ('needed' if c.get('b3_needed') else 'not needed')}"
          f"{(', ' + yn(c.get('b3'))) if c.get('b3_ran') else ''}; (vi) {yn(c.get('vi'))}; together {yn(c.get('hold'))}.",
          f"- Replication: mean g of the four new inputs {_fmt(rep.get('mean_g'), '+.4f')}, {rep.get('n_positive')} of "
          f"{rep.get('n')} positive (passes needs {PASS_MEAN} and {MIN_POSITIVE}; carries {CARRY_MEAN} and {MIN_POSITIVE}); "
          f"C1 {_fmt(v['replicate_kinds']['C1']['mean_g'], '+.4f')}, C2 {_fmt(v['replicate_kinds']['C2']['mean_g'], '+.4f')}."]
    if g.get('n'):
        L.append(f"- E83's gate on the mean of the two splits ({', '.join(g['inputs'])}): (i) mean g {_fmt(g['mean_g'], '+.4f')} "
                 f"(needs {GATE_GAIN}): {yn(g['i'])}; (ii) {g['n_singles_lower']} of {g['n_singles']} single wirings lower "
                 f"(needs {GATE_SINGLES}): {yn(g['ii'])}; (iii) {g['n_flagged']} points flagged, raw-count loss change "
                 f"{_fmt(g['raw_loss_change_mean'], '+.4f')} (max {_fmt(g['raw_loss_change_max'], '+.4f')}; needs at most "
                 f"{GATE_RAW_SLACK}), mean APL effect {_fmt(g['mean_effect_cal_mean'], '+.3f')} (farthest "
                 f"{_fmt(g['mean_effect_cal_farthest'], '+.3f')}; needs within {GATE_MEAN_SLACK} of {_fmt(v['flies_mean'], '+.3f')}): "
                 f"{yn(g['iii'])}; (iv) squared pattern error {_fmt(100 * g['sq_reduction_all_mean'], '.1f')} % lower against all "
                 f"recordings (needs 15), {_fmt(100 * g['sq_reduction_six_mean'], '.1f')} % against the six typical (needs 7.5): "
                 f"{yn(g['iv'])}; the gate {yn(g.get('holds'))}.")
    _section(out, '## Part D: the controls, the replication and the verdict', L + [''])


# ---------------------------------------------------------------------------
# part E (on "passes" or "carries"): what the input is

PART_E = 'E'
# the six inputs searched on the real wiring: E77's two (the candidate; E83's part A2 read them at SIZE on all eight
# wirings, archived, read, not rerun) and part C's four (on T from C-read, on S from this part)
SIX = ('cand_A', 'cand_B', 'C1_A', 'C1_B', 'C2_C', 'C2_D')
SIX_LABELS = {'cand_A': "the candidate, half A (E77's)", 'cand_B': "the candidate, half B (E77's)",
              'C1_A': 'second draw, half A (C1_A)', 'C1_B': 'second draw, half B (C1_B)',
              'C2_C': 'second split, half C (C2_C)', 'C2_D': 'second split, half D (C2_D)'}
# the six under E85's part A4 code: input -> (E85's name, (results folder, factor table)); E77's two keep E85's names,
# so that E85's rows of them are reproduced
E85_SIX = {'cand_A': ('E77_halfA', ('E77_input_side_common_error', 'oracle_factors_A.csv')),
           'cand_B': ('E77_halfB', ('E77_input_side_common_error', 'oracle_factors_B.csv')),
           **{n: (f'E86_{n}', (EXPERIMENT, f'oracle_factors_{n}.csv')) for n in SEARCHES}}
DRAWS7 = tuple(E83.POPULATIONS['draws7'])
E_UNITS = tuple(E83.POPULATIONS) + tuple(E83.ALL_WIRINGS)
# the list a measurement could test (fixed before any run): the entries (glomerulus x odor) on which at least E_AGREE of
# the six agree in sign with a mean log2 factor beyond +-E_LOG2 (the searched factors, at their full log size; the same
# rule on the change as read, SIZE times the log2 factor, is given beside); with fewer than E_MIN_ENTRIES the input is
# described by the pairs it moves and no measurement is named
E_AGREE, E_LOG2, E_MIN_ENTRIES = 5, 0.5, 10
# beside the rule, not part of it: its counts with each of the six tables label-permuted (E83's `permuted_factors`;
# table i of draw d with seed E_NULL_SEED0 + len(SIX) * d + i)
N_NULL_E, E_NULL_SEED0 = 1000, 100000
E_TOP_PAIRS = 12


def require_e(out: Path) -> dict:
    """verdict.json; raises unless its outcome is 'passes' or 'carries' (part E runs only then)."""
    v = json.loads((out / 'verdict.json').read_text())
    if v['outcome'] not in ('passes', 'carries'):
        raise RuntimeError(f"the verdict is {v['outcome']!r}: part E runs only on 'passes' or 'carries'")
    return v


def six_meta() -> Dict[str, dict]:
    """The six inputs' metadata (`_change_meta`'s)."""
    meta = _change_meta()
    return {n: meta[n] for n in SIX}


def six_factors(out: Path) -> Dict[str, pd.DataFrame]:
    """The six inputs' searched factors, glomeruli x the 16 odors (1 where not searched)."""
    F = {f'cand_{X}': candidate_factors(X) for X in ('A', 'B')}
    F.update({n: search_factors(out, n) for n in SEARCHES})
    return {n: F[n] for n in SIX}


def log2_factors(F: Mapping[str, pd.DataFrame], odors: Sequence[str] = E77.E_ODORS) -> Dict[str, pd.DataFrame]:
    """log2 of factor tables, restricted to `odors` (the ten searched)."""
    return {n: np.log2(f[list(odors)].astype(float)) for n, f in F.items()}


def plan_e(out: Path) -> List[dict]:
    """Part E's points: part C's four new inputs at SIZE on S (16 points, tags E_<input>_<wiring>; the change keeps the
    input's name, so that its readings on S and on T (C-read) make one input on all eight wirings)."""
    plan = []
    for name, cfg in SEARCHES.items():
        search_factors(out, name)
        plan += [dict(part=PART_E, change=name, group=name.split('_')[0], kind='search', seed=GENERATOR_SEEDS[name],
                      half=cfg['half'], unseen=UNSEEN[cfg['half']], size=SIZE, wiring=w, tag=f'E_{name}_{w}')
                 for w in E83.S_WIRINGS]
    return plan


def relation(a: Mapping, b: Mapping) -> str:
    """How two inputs' searches relate, from their 'half': 'same half' (the same recordings searched), 'same split' (the
    two halves of one split) or 'other split'."""
    if a['half'] == b['half']:
        return 'same half'
    return 'same split' if split_of(a['half']) == split_of(b['half']) else 'other split'


def factor_agreement(L: Mapping[str, pd.DataFrame], meta: Mapping[str, Mapping]) -> pd.DataFrame:
    """Pairwise agreement of log2 factor tables (glomeruli x odors, one layout): the correlation over all entries and,
    per odor, over its glomeruli (each odor's in factors_corr_<odor>; their median, min and max), with the pair's
    `relation`."""
    names = list(L)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            A, B = L[a], L[b].reindex(index=L[a].index, columns=L[a].columns)
            per_odor = [float(np.corrcoef(A[o].to_numpy(), B[o].to_numpy())[0, 1]) for o in A.columns]
            rows.append(dict(input_a=a, input_b=b, relation=relation(meta[a], meta[b]),
                             factors_corr=float(np.corrcoef(A.to_numpy().ravel(), B.to_numpy().ravel())[0, 1]),
                             factors_corr_per_odor_median=float(np.median(per_odor)),
                             factors_corr_per_odor_min=float(np.min(per_odor)),
                             factors_corr_per_odor_max=float(np.max(per_odor)),
                             **{f'factors_corr_{o}': c for o, c in zip(A.columns, per_odor)}))
    return pd.DataFrame(rows)


def change_effect(effects: Mapping[Tuple[str, str], pd.Series], name: str, unit: str) -> Optional[pd.Series]:
    """Input `name`'s change of the per-pair APL effect on `unit` (changed minus unchanged, the unchanged one's pair
    order), or None without both readings."""
    if (name, unit) not in effects or ('none', unit) not in effects:
        return None
    base = effects[('none', unit)]
    return (effects[(name, unit)].reindex(base.index) - base)


def changes_agreement(effects: Mapping[Tuple[str, str], pd.Series], names: Sequence[str],
                      units: Sequence[str]) -> pd.DataFrame:
    """Pairwise correlation of the inputs' per-pair changes (`change_effect`) on each unit (changes_corr_<unit>; NaN
    where either lacks a reading on it)."""
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            row = dict(input_a=a, input_b=b)
            for u in units:
                da, db = change_effect(effects, a, u), change_effect(effects, b, u)
                row[f'changes_corr_{u}'] = (np.nan if da is None or db is None else
                                            float(np.corrcoef(da.to_numpy(), db.to_numpy())[0, 1]))
            rows.append(row)
    return pd.DataFrame(rows)


def _rule(st: np.ndarray, min_agree: int, beyond: float, scale: float = 1.0) -> Tuple[np.ndarray, ...]:
    """The rule on a stack of log2 factors (tables first): (number above 0, below 0, the majority sign, how many share
    it, the mean, meets). An entry meets it when at least `min_agree` share the sign and `scale` times the mean is
    beyond +-`beyond` in that sign (0 counts for neither sign)."""
    npos, nneg = (st > 0).sum(axis=0), (st < 0).sum(axis=0)
    sign = np.sign(npos - nneg)
    n_agree = np.maximum(npos, nneg)
    mean = st.mean(axis=0)
    meets = (n_agree >= min_agree) & (sign != 0) & (sign * mean * scale > beyond)
    return npos, nneg, sign, n_agree, mean, meets


def five_of_six(L: Mapping[str, pd.DataFrame], min_agree: int = E_AGREE, beyond: float = E_LOG2,
                scale_read: float = SIZE) -> pd.DataFrame:
    """Per entry (odor, glomerulus) of a set of log2 factor tables (glomeruli x odors, one layout): each table's log2
    factor (a column per table), how many are above and below 0, the majority sign and how many share it (n_agree), the
    mean (mean_log2; mean_log2_read at `scale_read`), and the rule (`_rule`): `meets` on the log2 factors, `meets_read`
    on the change as read (`scale_read` times them)."""
    names = list(L)
    ref = L[names[0]]
    st = np.stack([L[n].reindex(index=ref.index, columns=ref.columns).to_numpy(dtype=float) for n in names])
    npos, nneg, sign, n_agree, mean, meets = _rule(st, min_agree, beyond)
    meets_read = _rule(st, min_agree, beyond, scale_read)[-1]
    rows = []
    for j, o in enumerate(ref.columns):
        for i, g in enumerate(ref.index):
            rows.append(dict(odor=o, glomerulus=g, **{n: float(st[k, i, j]) for k, n in enumerate(names)},
                             n_above=int(npos[i, j]), n_below=int(nneg[i, j]), sign=int(sign[i, j]),
                             n_agree=int(n_agree[i, j]), mean_log2=float(mean[i, j]),
                             mean_log2_read=float(scale_read * mean[i, j]), meets=bool(meets[i, j]),
                             meets_read=bool(meets_read[i, j])))
    return pd.DataFrame(rows)


def _factors_sha256(F: Mapping[str, pd.DataFrame]) -> str:
    import hashlib
    h = hashlib.sha256()
    for n, f in F.items():
        h.update(n.encode())
        h.update(f.to_csv(float_format='%.17g').encode())
    return h.hexdigest()


def five_of_six_null(F: Mapping[str, pd.DataFrame], n: int = N_NULL_E, seed0: int = E_NULL_SEED0,
                     odors: Sequence[str] = E77.E_ODORS, min_agree: int = E_AGREE, beyond: float = E_LOG2,
                     scale_read: float = SIZE) -> pd.DataFrame:
    """The rule's counts with every table's factors label-permuted (E83's `permuted_factors`: each odor's factors over
    its responding glomeruli), table i of draw d with seed seed0 + len(F) * d + i: one row per draw (n_meets,
    n_meets_read, n_agree: the entries with at least `min_agree` sharing the sign)."""
    names = list(F)
    rows = []
    for d in range(n):
        st = np.stack([np.log2(E83.permuted_factors(F[k], seed0 + len(names) * d + i)[list(odors)].to_numpy(dtype=float))
                       for i, k in enumerate(names)])
        r = _rule(st, min_agree, beyond)
        rows.append(dict(draw=d, n_meets=int(r[-1].sum()), n_meets_read=int(_rule(st, min_agree, beyond, scale_read)[-1].sum()),
                         n_agree=int((r[3] >= min_agree).sum())))
    return pd.DataFrame(rows)


def null_summary(observed: float, null: pd.Series) -> dict:
    """An observed count against its null: the null's mean, 95th percentile and the share of draws at or above it."""
    null = pd.Series(null, dtype=float)
    return dict(observed=float(observed), null_mean=float(null.mean()), null_p95=float(np.quantile(null, 0.95)),
                null_frac_ge=float((null >= observed).mean()), n_draws=int(len(null)))


def paired_differences(tb: pd.DataFrame, wirings: Sequence[str] = DRAWS7, against: str = 'unseen') -> pd.DataFrame:
    """Per input (change) of a `two_by_two` table, its single wirings `wirings` as paired differences against
    `against`: changed minus unchanged on each wiring (a column per wiring, NaN without a reading; negative is better),
    how many have a reading (n) and are lower (n_lower), their mean, and whether every wiring has one (complete)."""
    t = tb[(tb.against == against) & tb.unit.isin(list(wirings))]
    rows = []
    for ch in dict.fromkeys(t.change):
        g = t[t.change == ch].set_index('unit')
        d = (g.changed - g.unchanged).reindex(list(wirings))
        have = d.dropna()
        rows.append(dict(change=ch, against=against, **{w: float(d[w]) for w in wirings}, n=int(len(have)),
                         n_lower=int((have < 0).sum()), mean=float(have.mean()) if len(have) else np.nan,
                         complete=bool(len(have) == len(wirings))))
    return pd.DataFrame(rows)


def pairs_moved(effects: Mapping[Tuple[str, str], pd.Series], names: Sequence[str], unit: str,
                target: pd.Series) -> pd.DataFrame:
    """Per odor pair, the inputs' changes of the per-pair APL effect on `unit` (`change_effect`; inputs without a
    reading there left out), their mean, how many inputs share the mean's sign, what `target` asks (its per-pair effect
    minus the unchanged model's, both means removed) and the pair's group (ester, alcohol or other); sorted by the size
    of the mean change."""
    ch = {n: change_effect(effects, n, unit) for n in names}
    ch = {n: v for n, v in ch.items() if v is not None}
    if not ch:
        return pd.DataFrame()
    base = effects[('none', unit)]
    idx = base.index
    D = pd.DataFrame({n: v.reindex(idx).to_numpy() for n, v in ch.items()}, index=idx)
    mean = D.mean(axis=1)
    same = (np.sign(D.to_numpy()) == np.sign(mean.to_numpy())[:, None]).sum(axis=1)
    t = E83._cs(target).reindex(idx)
    asked = (t - t.mean()) - (base - base.mean())
    masks = E77.kind_masks(idx)
    group = np.where(masks['ester'], 'ester', np.where(masks['alcohol'], 'alcohol', 'other'))
    res = pd.DataFrame(dict(odor_a=idx.get_level_values(0), odor_b=idx.get_level_values(1), unit=unit, group=group,
                            **{n: D[n].to_numpy() for n in D.columns}, mean_change=mean.to_numpy(),
                            n_inputs=D.shape[1], n_same_sign=same, asked=asked.to_numpy()))
    return res.iloc[np.argsort(-np.abs(res.mean_change.to_numpy()), kind='stable')].reset_index(drop=True)


def agreement_by_odor(ag: pd.DataFrame, ent: pd.DataFrame, odors: Sequence[str] = E77.E_ODORS) -> pd.DataFrame:
    """Per odor: the inputs' pairwise correlations of its log2 factors (`factor_agreement`'s per-odor columns: mean,
    min, max over the pairs, and the mean over pairs searched on the same half), and of its entries (`five_of_six`) how
    many have at least E_AGREE sharing the sign and how many meet the rule."""
    rows = []
    for o in odors:
        c = ag[f'factors_corr_{o}']
        e = ent[ent.odor == o]
        rows.append(dict(odor=o, corr_mean=float(c.mean()), corr_min=float(c.min()), corr_max=float(c.max()),
                         corr_mean_same_half=float(c[ag.relation == 'same half'].mean()), n_entries=int(len(e)),
                         n_agree=int((e.n_agree >= E_AGREE).sum()), n_meets=int(e.meets.sum()),
                         n_meets_read=int(e.meets_read.sum())))
    return pd.DataFrame(rows)


def collect_six(out: Path) -> Tuple[Dict[Tuple[str, str], Dict[str, tuple]], pd.DataFrame]:
    """The six inputs' runs so far, {(change, wiring): outputs}, and their rows: the candidate's E83 part A2 runs at
    SIZE on all eight wirings (archived, not rerun; rows from E83's runs_A2.csv), part C's four on T (runs_Cread.csv)
    and on S (runs_E.csv)."""
    runs, rows = {}, []
    a2 = pd.read_csv(H.results_dir(E83.EXPERIMENT) / 'runs_A2.csv', float_precision='round_trip')
    for X in ('A', 'B'):
        for w in E83.ALL_WIRINGS:
            tag = f'A2_half{X}_s{SIZE:g}_{w}'
            o = E83.load_run('A2', tag)
            if o is None:
                raise RuntimeError(f"E83's run {tag} not found")
            runs[(f'cand_{X}', w)] = o
            rows.append(dict(a2[a2.tag == tag].iloc[0].to_dict(), part='A2', change=f'cand_{X}', group='candidate',
                             kind='candidate', half=X, unseen=UNSEEN[X],
                             source=f"E83 part A2 ({tag}, archived, not rerun)"))
    for part in ('Cread', PART_E):
        for r in E83._read_rows(out / f'runs_{part}.csv').to_dict('records'):
            if r.get('change') not in SIX:
                continue
            o = load_run(part, r['tag'])
            if o is None:
                print(f"{r['tag']}: outputs missing", flush=True)
                continue
            runs[(r['change'], r['wiring'])] = o
            rows.append(r)
    return runs, pd.DataFrame(rows)


def score_e(out: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The six inputs' runs so far (`collect_six`) scored with E83's `score_runs` against this experiment's targets
    (populations only where all their wirings have a reading): partE_readings.csv, partE_pairs.csv and
    partE_two_by_two.csv (`two_by_two` on every population and wiring). T's g of each input must equal part D's
    inputs.csv (raises otherwise). Returns (readings, pairs, two_by_two, runs rows)."""
    runs, rt = collect_six(out)
    meta = six_meta()
    sm = {k: dict({a: b for a, b in v.items() if a != 'kind'}, input_kind=v['kind']) for k, v in meta.items()}
    readings, pairs = E83.score_runs(runs, sm, tg=targets())
    readings.to_csv(out / 'partE_readings.csv', index=False)
    pairs.to_csv(out / 'partE_pairs.csv', index=False)
    tb = two_by_two(readings, meta, units=E_UNITS)
    tb.to_csv(out / 'partE_two_by_two.csv', index=False)
    it = pd.read_csv(out / 'inputs.csv', float_precision='round_trip').set_index('change')
    for n in SIX:
        got = tb[(tb.change == n) & (tb.unit == 'T') & (tb.against == 'unseen')].gain
        if len(got) != 1 or abs(float(got.iloc[0]) - float(it.loc[n, 'g'])) > 1e-12:
            raise RuntimeError(f"{n}: T's g {list(got)} differs from part D's {float(it.loc[n, 'g'])}")
    return readings, pairs, tb, rt


def inputs_e(readings: pd.DataFrame, tb: pd.DataFrame, rt: pd.DataFrame) -> pd.DataFrame:
    """One row per input of the six: g (the pattern part against the unseen half of its split, unchanged minus changed)
    on S, T, the seven re-drawn wirings as a population (draws7) and all eight, and against all recordings; the seven
    as paired differences (`paired_differences`: number lower, mean, against the unseen half and against all
    recordings); the raw-count loss change against all recordings and the mean APL effect under the readout on the
    seven and on all eight; the readings flagged by the four bands (the rule of 16:10) of all its readings."""
    meta = six_meta()
    pu = paired_differences(tb, DRAWS7, 'unseen').set_index('change')
    pa = paired_differences(tb, DRAWS7, 'all').set_index('change')
    rows = []
    for n in SIX:
        m = meta[n]
        row = dict(input=n, label=SIX_LABELS[n], group=m['group'], searched_on=m['half'], unseen=m['unseen'],
                   split=split_of(m['half']))
        for u in ('S', 'T', 'draws7', 'all8'):
            q = tb[(tb.change == n) & (tb.unit == u)].set_index('against')
            row[f'g_{u}'] = float(q.loc['unseen', 'gain']) if 'unseen' in q.index else np.nan
            row[f'g_all_{u}'] = float(q.loc['all', 'gain']) if 'all' in q.index else np.nan
        for p, k in ((pu, 'draws7'), (pa, 'draws7_all')):
            row[f'{k}_n'] = int(p.loc[n, 'n']) if n in p.index else 0
            row[f'{k}_lower'] = int(p.loc[n, 'n_lower']) if n in p.index else 0
            row[f'{k}_mean_diff'] = float(p.loc[n, 'mean']) if n in p.index else np.nan
        for u in ('draws7', 'all8'):
            raw = readings[(readings.change == n) & (readings.unit == u) & (readings.readout == 'raw')]
            raw0 = readings[(readings.change == 'none') & (readings.unit == u) & (readings.readout == 'raw')]
            cal = readings[(readings.change == n) & (readings.unit == u) & (readings.readout == 'cal')]
            row[f'raw_loss_change_{u}'] = float(raw.z8_loss.iloc[0] - raw0.z8_loss.iloc[0]) if len(raw) else np.nan
            row[f'mean_effect_cal_{u}'] = float(cal.mean_effect.iloc[0]) if len(cal) else np.nan
        rr = rt[rt.change == n] if len(rt) else rt
        flags = rr.realism_flag.fillna('').astype(str).str.strip() != '' if len(rr) else pd.Series(dtype=bool)
        row.update(n_points=int(len(rr)), n_flagged=int(flags.sum()),
                   wirings=','.join(w for w in E83.ALL_WIRINGS if w in set(rr.get('wiring', []))))
        rows.append(row)
    return pd.DataFrame(rows)


def vs_flies() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """The six inputs against the ORN flies' spread and Hallem's data through E85's part A4 code
    (`E85_orn_fly_population.searched_tables` and `searched_summary`, at E85's SEARCH_SIZE, which is SIZE): per entry,
    per input and odor, and per input, with the correlation of the change (log2, as read) with Hallem's ratio to the
    lab's input over the entries E38's table has. E77's two keep E85's names, and their summary rows are compared with
    E85's searched_inputs_summary.csv (the largest absolute difference is returned with the tables)."""
    from experiments import E85_orn_fly_population as E85
    if float(E85.SEARCH_SIZE) != SIZE:
        raise RuntimeError(f"E85's part A4 reads the inputs at {E85.SEARCH_SIZE}, not {SIZE}")
    searched = {E85_SIX[n][0]: E85_SIX[n][1] for n in SIX}
    entries, span = E85.searched_tables(searched)
    summ = E85.searched_summary(entries, span)
    hc = []
    for name, e in entries.groupby('input', sort=False):
        h = e.dropna(subset=['hallem_log2_vs_lab'])
        hc.append(dict(input=name, n_hallem=int(len(h)), hallem_odors=','.join(sorted(set(h.odor))),
                       hallem_corr=float(np.corrcoef(h.log2_change, h.hallem_log2_vs_lab)[0, 1]) if len(h) >= 3 else np.nan))
    summ = summ.merge(pd.DataFrame(hc), on='input', how='left')
    inv = {v[0]: k for k, v in E85_SIX.items()}
    for df in (entries, span, summ):
        df.insert(0, 'e86_input', df.input.map(inv))
    ref = pd.read_csv(H.results_dir(E85.EXPERIMENT) / 'searched_inputs_summary.csv', float_precision='round_trip')
    ref = ref.set_index('input')
    s = summ.set_index('input')
    cols = [c for c in ref.columns if c in s.columns and pd.api.types.is_numeric_dtype(ref[c])]
    diffs = {k: float(np.nanmax(np.abs(s.loc[k, cols].to_numpy(dtype=float) - ref.loc[k, cols].to_numpy(dtype=float))))
             for k in (E85_SIX['cand_A'][0], E85_SIX['cand_B'][0]) if k in ref.index}
    check = dict(e85_rows=list(diffs), columns=cols, max_abs_diff=max(diffs.values()) if diffs else None,
                 reproduced=bool(diffs) and max(diffs.values()) <= 1e-12)
    return entries, span, summ, check


def annotate_entries(ent: pd.DataFrame) -> pd.DataFrame:
    """`five_of_six`'s entries with the lab's input, the ORN flies' mean and spread (E85's `spread_table`: SD across
    flies of log2(fly / mean)), the mean change as read in units of that spread, the entry's share of the flies'
    largest mean response, and Hallem's log2 ratio to the lab's input where E38's table has it (E85's `hallem_log2`)."""
    from experiments import E85_orn_fly_population as E85
    sp = E85.spread_table().set_index(['odor', 'glomerulus'])
    hal = E85.hallem_log2()
    lab = H.abbrev_columns(H.apl16_orn_deltas())
    keys = list(zip(ent.odor, ent.glomerulus))
    get = lambda k, c: float(sp.loc[k, c]) if k in sp.index else np.nan
    res = ent.copy()
    res['lab_input'] = [float(lab.loc[g, o]) for o, g in keys]
    res['flies_mean'] = [get(k, 'mean') for k in keys]
    res['flies_share_of_largest'] = [get(k, 'share_of_largest') for k in keys]
    res['flies_sd_log2'] = [get(k, 'sd_log2_ratio') for k in keys]
    res['read_in_flies_sd'] = res.mean_log2_read / res.flies_sd_log2
    res['hallem_log2_vs_lab'] = [float(hal.get(k, np.nan)) for k in keys]
    return res


def part_e_tables(out: Path, n_null: int = N_NULL_E) -> dict:
    """Part E's tables from whatever runs exist (the six on T exist before any of this part's runs): partE_readings.csv,
    partE_pairs.csv, partE_two_by_two.csv (`score_e`); partE_inputs.csv (`inputs_e`); partE_wirings.csv (the seven
    re-drawn wirings as paired differences, against the unseen half and all recordings); partE_agreement.csv (pairwise:
    log2 factors over all entries and per odor, per-pair changes on T, the seven and all eight) and
    partE_agreement_by_odor.csv; partE_entries.csv (every entry under the rule, annotated) and partE_five_of_six.csv
    (those that meet it); partE_null.csv (the rule's label-permuted counts; reused while the six factor tables are
    unchanged); partE_pairs_moved.csv (on T and the seven, against all recordings); partE_vs_flies.csv,
    partE_vs_flies_span.csv, partE_vs_flies_summary.csv (E85's part A4 code); the archive tables' rows of this part's
    runs; part_e.json and the "Part E" section of comparison_table.md."""
    v = require_e(out)
    t0 = time.time()
    readings, pairs, tb, rt = score_e(out)
    ie = inputs_e(readings, tb, rt)
    ie.to_csv(out / 'partE_inputs.csv', index=False)
    pdiff = pd.concat([paired_differences(tb, DRAWS7, a) for a in ('unseen', 'all')], ignore_index=True)
    pdiff.to_csv(out / 'partE_wirings.csv', index=False)
    meta = six_meta()
    F = six_factors(out)
    L = log2_factors(F)
    effects = _e83_effects('cal', out / 'partE_pairs.csv')
    change_units = ('T', 'draws7', 'all8')
    ag = factor_agreement(L, meta).merge(changes_agreement(effects, SIX, change_units), on=['input_a', 'input_b'])
    ag.to_csv(out / 'partE_agreement.csv', index=False)
    ent = annotate_entries(five_of_six(L))
    ent.to_csv(out / 'partE_entries.csv', index=False)
    listed = ent[ent.meets].copy()
    listed = listed.iloc[np.argsort(-np.abs(listed.mean_log2.to_numpy()), kind='stable')]
    listed.to_csv(out / 'partE_five_of_six.csv', index=False)
    bo = agreement_by_odor(ag, ent)
    bo.to_csv(out / 'partE_agreement_by_odor.csv', index=False)
    sha = _factors_sha256(F)
    np_path = out / 'partE_null.csv'
    null = pd.read_csv(np_path) if np_path.exists() else pd.DataFrame()
    if not (len(null) == n_null and 'factors_sha256' in null.columns and (null.factors_sha256 == sha).all()):
        print(f'the rule under {n_null} label permutations of the six', flush=True)
        null = five_of_six_null(F, n_null).assign(factors_sha256=sha)
        null.to_csv(np_path, index=False)
    pm = pd.concat([pairs_moved(effects, SIX, u, targets()['z8']) for u in ('T', 'draws7')], ignore_index=True)
    pm.to_csv(out / 'partE_pairs_moved.csv', index=False)
    fe, fs, fsum, fcheck = vs_flies()
    fe.to_csv(out / 'partE_vs_flies.csv', index=False)
    fs.to_csv(out / 'partE_vs_flies_span.csv', index=False)
    fsum.to_csv(out / 'partE_vs_flies_summary.csv', index=False)
    rt_e = rt[rt.part == PART_E] if len(rt) else rt
    if len(rt_e):
        _archive_tables(out, rt_e, of_parts=(PART_E,))

    rule = dict(min_agree=E_AGREE, beyond_log2=E_LOG2, min_entries=E_MIN_ENTRIES, n_entries=int(len(ent)),
                n_agree=int((ent.n_agree >= E_AGREE).sum()), n_meets=int(ent.meets.sum()),
                n_meets_read=int(ent.meets_read.sum()),
                null_meets=null_summary(ent.meets.sum(), null.n_meets),
                null_meets_read=null_summary(ent.meets_read.sum(), null.n_meets_read),
                null_agree=null_summary((ent.n_agree >= E_AGREE).sum(), null.n_agree),
                null_seed0=E_NULL_SEED0, factors_sha256=sha)
    rule['description'] = 'list' if rule['n_meets'] >= E_MIN_ENTRIES else 'pairs'
    by_rel = lambda col: {r: float(g[col].mean()) for r, g in ag.groupby('relation') if g[col].notna().any()}
    agreement = dict(factors_corr_mean=float(ag.factors_corr.mean()), factors_corr_min=float(ag.factors_corr.min()),
                     factors_corr_max=float(ag.factors_corr.max()), factors_corr_by_relation=by_rel('factors_corr'),
                     factors_corr_per_odor_median_mean=float(ag.factors_corr_per_odor_median.mean()))
    for u in change_units:
        c = ag[f'changes_corr_{u}']
        agreement[f'changes_corr_{u}_mean'] = float(c.mean()) if c.notna().any() else None
        agreement[f'changes_corr_{u}_by_relation'] = by_rel(f'changes_corr_{u}') if c.notna().any() else None
    moved = {}
    for u, g in pm.groupby('unit', sort=False):
        moved[u] = dict(n_inputs=int(g.n_inputs.iloc[0]), corr_mean_change_asked=float(np.corrcoef(g.mean_change, g.asked)[0, 1]),
                        n_pairs_5of6_same_sign=int((g.n_same_sign >= E_AGREE).sum()),
                        groups={k: dict(n=int(len(q)), mean_change=float(q.mean_change.mean()), asked=float(q.asked.mean()))
                                for k, q in g.groupby('group')})
    new = ie[ie.input.isin(list(SEARCHES))]
    res = dict(script=SCRIPT, generated=datetime.now(timezone.utc).isoformat(timespec='seconds'), **H.code_commits(),
               verdict=v['outcome'], size=SIZE, six=list(SIX), seconds=time.time() - t0,
               runs=dict(n_E=int(len(rt_e)), planned_E=len(SEARCHES) * len(E83.S_WIRINGS),
                         n_flagged_E=int((rt_e.realism_flag.fillna('').astype(str).str.strip() != '').sum()) if len(rt_e) else 0),
               inputs={r['input']: {k: r[k] for k in ('g_S', 'g_T', 'g_draws7', 'g_all8', 'g_all_draws7', 'draws7_n',
                                                      'draws7_lower', 'draws7_mean_diff', 'raw_loss_change_draws7',
                                                      'mean_effect_cal_draws7', 'n_points', 'n_flagged')}
                       for r in ie.to_dict('records')},
               means=dict(six_g_draws7=float(ie.g_draws7.mean()), new_g_draws7=float(new.g_draws7.mean()),
                          six_draws7_mean_diff=float(ie.draws7_mean_diff.mean()),
                          six_draws7_lower=int(ie.draws7_lower.sum()), six_draws7_n=int(ie.draws7_n.sum())),
               agreement=agreement, rule=rule, pairs_moved=moved,
               flies=dict(e85_check=fcheck, summary=fsum.drop(columns=['input']).set_index('e86_input').to_dict('index')),
               flies_mean_effect=float(targets()['z8'].mean()),
               sources=dict(candidate_runs=f'data/model_outputs/{E83.EXPERIMENT}/A2 (E83 part A2, archived)',
                            candidate_rows=f'results/{E83.EXPERIMENT}/runs_A2.csv',
                            factors={n: f'results/{E85_SIX[n][1][0]}/{E85_SIX[n][1][1]}' for n in SIX},
                            e85_code='experiments/E85_orn_fly_population.py (searched_tables, searched_summary, '
                                     'spread_table, hallem_log2)'))
    (out / 'part_e.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    _comparison_table_e(out, ie, pdiff, ag, bo, listed, pm, fsum, fcheck, res)
    print(json.dumps(dict(runs=res['runs'], means=res['means'], agreement=agreement,
                          rule={k: rule[k] for k in ('n_agree', 'n_meets', 'n_meets_read', 'description')},
                          null=rule['null_meets'], e85_check=fcheck), indent=1, default=float), flush=True)
    print(ie[['input', 'g_S', 'g_T', 'g_draws7', 'g_all8', 'draws7_lower', 'draws7_n', 'draws7_mean_diff',
              'n_flagged']].round(4).to_string(index=False), flush=True)
    return res


def _comparison_table_e(out: Path, ie: pd.DataFrame, pdiff: pd.DataFrame, ag: pd.DataFrame, bo: pd.DataFrame,
                        listed: pd.DataFrame, pm: pd.DataFrame, fsum: pd.DataFrame, fcheck: dict, res: dict) -> None:
    rule = res['rule']
    n_e, planned = res['runs']['n_E'], res['runs']['planned_E']
    L = ['', f"The verdict is **{res['verdict']}**, so part E runs. The six inputs searched on the real wiring: E77's two "
         "(the candidate; E83's part A2 read them at 0.6 on all eight wirings, archived, not rerun) and part C's four "
         f"(on T from C-read, on S from this part: {n_e} of {planned} readings so far). Each at 0.6 of its log size, both "
         'response rates re-held. Pattern part of the APL-effect loss under the calibrated readout, against the half of '
         "each input's split that its search never saw (lower is better).", '',
         "### The six on the seven re-drawn wirings (paired differences against the unseen half)", '',
         'Changed minus unchanged on each wiring (negative is better); – where a wiring has no reading yet.', '',
         '| input | searched on / scored on | ' + ' | '.join(DRAWS7) + ' | lower | mean |',
         '|---|---|' + '---:|' * (len(DRAWS7) + 2)]
    pu = pdiff[pdiff.against == 'unseen'].set_index('change')
    for r in ie.to_dict('records'):
        q = pu.loc[r['input']] if r['input'] in pu.index else None
        cells = [_fmt(q[w], '+.4f') if q is not None else '–' for w in DRAWS7]
        L.append(f"| {r['label']} | {r['searched_on']} / {r['unseen']} | " + ' | '.join(cells) +
                 f" | {r['draws7_lower']} of {r['draws7_n']} | {_fmt(r['draws7_mean_diff'], '+.4f')} |")
    L += ['', 'As populations (the mean over their wirings of the pair distances per condition); g = unchanged minus '
          'changed (positive is better):', '',
          '| input | g: S | g: T | g: the seven re-drawn | g: all eight | g against all recordings: the seven | '
          'raw-count loss change: the seven | mean APL effect: the seven | readings flagged |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in ie.to_dict('records'):
        L.append(f"| {r['label']} | {_fmt(r['g_S'], '+.4f')} | {_fmt(r['g_T'], '+.4f')} | {_fmt(r['g_draws7'], '+.4f')} | "
                 f"{_fmt(r['g_all8'], '+.4f')} | {_fmt(r['g_all_draws7'], '+.4f')} | "
                 f"{_fmt(r['raw_loss_change_draws7'], '+.4f')} | {_fmt(r['mean_effect_cal_draws7'], '+.3f')} | "
                 f"{r['n_flagged']} of {r['n_points']} |")
    L.append(f"| KC data | | | | | | | {_fmt(res['flies_mean_effect'], '+.3f')} | |")
    m = res['means']
    L += ['', f"- Over the six: mean g on the seven {_fmt(m['six_g_draws7'], '+.4f')} (the four new "
          f"{_fmt(m['new_g_draws7'], '+.4f')}); {m['six_draws7_lower']} of {m['six_draws7_n']} single re-drawn wirings "
          f"lower, mean change {_fmt(m['six_draws7_mean_diff'], '+.4f')}."]
    a = res['agreement']
    L += ['', '### The agreement of the six', '',
          '| inputs | searched | log2 factors, all 360 entries | per odor, median (min to max) | per-pair change on T | '
          'on the seven re-drawn | on all eight |', '|---|---|---:|---:|---:|---:|---:|']
    for r in ag.to_dict('records'):
        L.append(f"| {SIX_LABELS[r['input_a']]} / {SIX_LABELS[r['input_b']]} | {r['relation']} | "
                 f"{_fmt(r['factors_corr'], '.2f')} | {_fmt(r['factors_corr_per_odor_median'], '.2f')} "
                 f"({_fmt(r['factors_corr_per_odor_min'], '.2f')} to {_fmt(r['factors_corr_per_odor_max'], '.2f')}) | "
                 f"{_fmt(r['changes_corr_T'], '.2f')} | {_fmt(r['changes_corr_draws7'], '.2f')} | "
                 f"{_fmt(r['changes_corr_all8'], '.2f')} |")
    rel = lambda d: '; '.join(f'{k} {_fmt(x, ".2f")}' for k, x in (d or {}).items())
    L += ['', f"- Mean over the 15 pairs: log2 factors {_fmt(a['factors_corr_mean'], '.2f')} ({rel(a['factors_corr_by_relation'])}); "
          f"per-pair changes on T {_fmt(a['changes_corr_T_mean'], '.2f')} ({rel(a['changes_corr_T_by_relation'])})"
          + (f"; on the seven {_fmt(a['changes_corr_draws7_mean'], '.2f')} ({rel(a['changes_corr_draws7_by_relation'])})"
             if a.get('changes_corr_draws7_mean') is not None else '') + '.', '',
          '| odor | log2 factors: mean correlation over the 15 pairs (min to max) | same half | entries with five of six '
          'agreeing in sign (of 36) | meeting the rule | the rule at 0.6 |', '|---|---:|---:|---:|---:|---:|']
    for r in bo.to_dict('records'):
        L.append(f"| {r['odor']} | {_fmt(r['corr_mean'], '.2f')} ({_fmt(r['corr_min'], '.2f')} to {_fmt(r['corr_max'], '.2f')}) | "
                 f"{_fmt(r['corr_mean_same_half'], '.2f')} | {r['n_agree']} | {r['n_meets']} | {r['n_meets_read']} |")
    nm, nr = rule['null_meets'], rule['null_meets_read']
    L += ['', '### The entries on which five of the six agree', '',
          f"- The rule (fixed before any run): at least {E_AGREE} of the six agree in sign and the mean log2 factor is "
          f"beyond ±{E_LOG2}. {rule['n_meets']} of {rule['n_entries']} entries meet it ({rule['n_agree']} have five of six "
          f"agreeing in sign). On the change as read (0.6 of the log2 factor): {rule['n_meets_read']}.",
          f"- Beside the rule, not part of it: with each input's factors label-permuted over each odor's glomeruli (E83's "
          f"`permuted_factors`, {nm['n_draws']} draws) {_fmt(nm['null_mean'], '.1f')} entries meet it on average (95th "
          f"percentile {_fmt(nm['null_p95'], '.0f')}; {_fmt(100 * nm['null_frac_ge'], '.1f')} % of draws as many or more); "
          f"at 0.6, {_fmt(nr['null_mean'], '.1f')} (95th percentile {_fmt(nr['null_p95'], '.0f')}; "
          f"{_fmt(100 * nr['null_frac_ge'], '.1f')} %).",
          ('- Ten or more entries: the list below is the one a measurement could test.' if rule['description'] == 'list' else
           '- Fewer than ten entries: the input is described by the pairs it moves (below), and no measurement is named.')]
    if rule['description'] == 'list':
        L += ['', '| odor | glomerulus | log2 factors of the six (' + ' / '.join(SIX) + ') | mean | as read (0.6) | lab input '
              '(spike deltas) | the flies\' SD (log2) | as read, in the flies\' SDs | Hallem, log2 ratio to the lab input |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|']
        for r in listed.to_dict('records'):
            L.append(f"| {r['odor']} | {r['glomerulus']} | " + ' / '.join(_fmt(r[n], '+.2f') for n in SIX) +
                     f" | {_fmt(r['mean_log2'], '+.2f')} | {_fmt(r['mean_log2_read'], '+.2f')} | {_fmt(r['lab_input'], '.1f')} | "
                     f"{_fmt(r['flies_sd_log2'], '.2f')} | {_fmt(r['read_in_flies_sd'], '+.1f')} | "
                     f"{_fmt(r['hallem_log2_vs_lab'], '+.2f')} |")
    L += ['', '### The pairs the six move (against all recordings)', '']
    for u, name in (('T', 'T'), ('draws7', 'the seven re-drawn wirings')):
        g = pm[pm.unit == u] if len(pm) else pm
        if not len(g):
            continue
        mv = res['pairs_moved'][u]
        L += [f"- On {name} ({mv['n_inputs']} inputs): the correlation of the mean change of the per-pair APL effect with "
              f"what all recordings ask (their effect minus the unchanged model's, both means removed) "
              f"{_fmt(mv['corr_mean_change_asked'], '.2f')}; {mv['n_pairs_5of6_same_sign']} of 120 pairs moved the same way "
              f"by at least {E_AGREE} inputs; by group (mean change / asked): " +
              '; '.join(f"{k} ({q['n']}) {_fmt(q['mean_change'], '+.4f')} / {_fmt(q['asked'], '+.4f')}"
                        for k, q in mv['groups'].items()) + '.']
    u = 'draws7' if 'draws7' in res['pairs_moved'] and res['pairs_moved']['draws7']['n_inputs'] == len(SIX) else 'T'
    g = pm[pm.unit == u].head(E_TOP_PAIRS) if len(pm) else pm
    if len(g):
        L += ['', f"The {E_TOP_PAIRS} pairs with the largest mean change on {'the seven re-drawn wirings' if u == 'draws7' else 'T'}:", '',
              '| pair | group | mean change | inputs sharing its sign | what all recordings ask |', '|---|---|---:|---:|---:|']
        for r in g.to_dict('records'):
            L.append(f"| {r['odor_a']}, {r['odor_b']} | {r['group']} | {_fmt(r['mean_change'], '+.4f')} | "
                     f"{r['n_same_sign']} of {r['n_inputs']} | {_fmt(r['asked'], '+.4f')} |")
    L += ['', "### Against the ORN flies and Hallem's data (E85's part A4 code)", '',
          '| input | entries within one SD of the flies\' spread (log2) | within two | rms log2 change | share of the squared '
          'size inside the nine flies\' span (permuted mean) | odors above the permuted 95th percentile | correlation with '
          'Hallem\'s ratio to the lab input (entries) |', '|---|---:|---:|---:|---:|---:|---:|']
    for r in fsum.to_dict('records'):
        L.append(f"| {SIX_LABELS[r['e86_input']]} | {_fmt(100 * r['within1_log2'], '.0f')} % | {_fmt(100 * r['within2_log2'], '.1f')} % | "
                 f"{_fmt(r['rms_log2_change'], '.2f')} | {_fmt(r['span_share'], '.2f')} ({_fmt(r['span_share_perm_mean'], '.2f')}) | "
                 f"{int(r['n_odors_above_perm_p95'])} of {int(r['n_odors'])} | {_fmt(r['hallem_corr'], '.2f')} ({int(r['n_hallem'])}) |")
    L += ['', f"- E85's rows of E77's two inputs are {'reproduced' if fcheck['reproduced'] else 'not reproduced'} by the same "
          f"code here (largest difference {_fmt(fcheck['max_abs_diff'], '.1g')})."]
    _section(out, '## Part E: what the input is (no simulations, then the four new inputs on S)', L + [''])


def part_e(out: Path, simulate: bool = True) -> dict:
    """Part E, only on "passes" or "carries" (`require_e`): the tables that need no simulation first
    (`part_e_tables` on the runs that exist), then (with `simulate`) the four new inputs on S (`plan_e`, 16 held
    readings; resumable: runs_E.csv, provenance_E.json, outputs under outputs/E86_reduced_input_candidate/E/), then the
    tables again on all eight wirings."""
    require_e(out)
    print('part E: the tables that need no simulation', flush=True)
    res = part_e_tables(out)
    if simulate:
        run_points(out, PART_E, plan_e(out))
        print('part E: the tables on all eight wirings', flush=True)
        res = part_e_tables(out)
    return res


# ---------------------------------------------------------------------------
# part F (only on "passes")

def part_f(out: Path) -> None:
    """Part F, only on "passes" (not written; needs E83's build_simulator with the family and the pulse as arguments)."""
    raise NotImplementedError('part F runs only on "passes"; not written yet')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=['A', 'identity', 'B1', 'B2', 'B3', 'B', 'C-check', 'C', 'C-read',
                                                      'score', 'D', 'E', 'E-tables', 'F'])
    ap.add_argument('--search', choices=list(SEARCHES), help='part C: the search to run')
    ap.add_argument('--no-resume', action='store_true', help='part C: start the search over (ignore its checkpoint)')
    ap.add_argument('--no-simulate', action='store_true', help="part C-check: compare the last run's history only")
    args = ap.parse_args()
    out = H.results_dir(EXPERIMENT)
    if args.part == 'A':
        part_a(out)
    elif args.part == 'identity':
        identity(out)
    elif args.part in ('B1', 'B2', 'B3'):
        run_points(out, args.part, plan_b(args.part))
    elif args.part == 'B':
        for part in ('B1', 'B2'):
            run_points(out, part, plan_b(part))
        it = score(out)
        if b3_needed(it):
            print('a B1 input has g above 0: B3 runs', flush=True)
            run_points(out, 'B3', plan_b('B3'))
            score(out)
        else:
            print('no B1 input has g above 0: B3 does not run', flush=True)
    elif args.part == 'C-check':
        check_search(out, simulate=not args.no_simulate)
    elif args.part == 'C':
        if args.search is None:
            ap.error('--part C needs --search ' + '|'.join(SEARCHES))
        search(out, args.search, resume=not args.no_resume)
    elif args.part == 'C-read':
        run_points(out, 'Cread', plan_c_read(out))
    elif args.part == 'score':
        score(out)
    elif args.part == 'D':
        score(out)
        write_verdict(out)
    elif args.part == 'E':
        part_e(out)
    elif args.part == 'E-tables':
        part_e(out, simulate=False)
    else:
        part_f(out)


if __name__ == '__main__':
    main()
