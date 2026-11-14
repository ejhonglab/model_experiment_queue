#!/usr/bin/env python3
"""E93: the measured input for the glomeruli that the model's input leaves at zero.

The entry is item E93 of experiment_priorities.md (proposed by viki 2026-11-13 20:37; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-13 20:12). The model's ORN input (`H.apl16_orn_deltas()`) has the 36
glomeruli both of the lab's panels have; the hemibrain wiring has 54 and `fit_mb_model` zero-fills the other 18. For ten
of them responses to the panel's odors have been measured: VA4 and VL1 in the lab's megamat ORN recordings
(`megamat_orn_deltas()`, the 12 megamat odors), DL4 and VM4 in its validation recordings (`validation2_orn_deltas()`,
paa, -aPine, euc, B-myr), and DA4m, DA3, DA4l, VA1v, DL3, VA1d (15 odors, all but euc) and DL4 (the 12 megamat odors) in
Hallem's receptor data (drosolf, as E38's `hallem_deltas`; 10^-2), each odor's values multiplied by the slope through the
origin of the lab's input on Hallem's over the 16 glomeruli both have for that odor. Unmeasured entries stay 0 (26 of the
160).

C2's family (the breadth-2 APL-effect fit of E67, 6.8 %, the 0.5 s pulse), each wiring at its own fit of E72's part B
(E83's `saved_fit`). S: real, P1, P2, M1. T: P3, P4, M2, M3. Halves A / B of the recordings E77's (`recording_halves()`,
seed 0), C / D the second split's (seed 1; E89's `targets89`). Every input is fixed before any run and none is fitted to
the KC recordings.

The inputs (INPUTS; glomeruli x the 16 odors, the lab's 36 rows unchanged in all of them; inputs/<name>.csv with sha256):
- I0, the lab's input (the reference: the saved runs); I0z, I0 with the ten rows written out as zeros (part 0).
- F, the fuller input (the primary): I0 plus the ten rows above.
- Bounds, on S only: L (the lab's four rows, no Hallem), U (F with Hallem's values as recorded at 10^-2), V (VA4's row
  alone).
- Controls, on T: PL<seed> (F with the ten rows' labels permuted among the ten glomeruli, no row left in place) and
  PO<seed> (F with each added row's values permuted over the odors, megamat and validation odors separately). Seed 1
  each; seed 2 each too if F moves T's population (pattern part under the readout against all recordings) by
  CONTROL_DRAW2 or more either way.
- DD, only if part A2 says "in the list's direction": F with the PN dendrites' departures put on the input at the
  entries whose halves agree in sign, bounded to a factor of DD_BOUND (on all eight wirings).

Parts (CLI `--part`):
- inputs: the inputs, added_rows.csv, hallem_scaling.csv, inputs/ (no simulations).
- 0, the identity: on each wiring, I0z's simulator (the breadth recomputed for it, `E85.build`, which also checks that
  `J.bouton_base` returns the same wiring and APL weights for a 46-row input) reproduces the saved APL-on and APL-off runs
  with 0 count differences, and its breadth equals the lab input's exactly (identity.json). If not, the experiment stops
  there (outcome (h)).
- A (no simulations): A1 what is added (a1_reliability.csv: halves of the flies for VA4, VL1 (olfdb dataset 1), DL4,
  VM4 (dataset 4), VA4 in the ORNs against the PN dendrites; hallem_scaling.csv; input_shares.csv: each odor's added
  share of claw-weighted positive input; input_distances.csv: the input distances by group with and without the rows;
  who_responds.csv: in the saved runs of the eight wirings the response rate and the odors per KC by a KC's share of
  zero-filled claws, of claws on the ten, on VA4 and on DA4l or DA3). A2 the PN dendrites at the searched input's
  entries (E86's partE_entries.csv on the glomeruli the GH146, pebbled and model tables share): the departure, log2 of
  the dendrites' relative response over the ORNs' (E90's `al_table`), normalized within a glomerulus, within an odor and
  with the model's own PN stage (E90's saved PN rates) in place of the ORNs; intervals from DEP_BOOT resamples of both
  tables' flies; the sign agreement and ρ with the six searches' mean log2 factor against a null of the departures
  permuted within each odor (NULL_DRAWS); the halves' agreement; the rule (A2_RULE, fixed in the entry).
  dendrite_departures.csv, dendrite_reading.json.
- B: the inputs at fixed scales with both response rates re-held (E85's `held_reading`: the breadth recomputed for the
  input, the somatic threshold bisected to the fit's APL-off count, the tied APL scale to its APL-on count); F on all
  eight, L, U, V on S, the controls on T (and DD if A2 says so). Resumable; runs_B.csv (absolute parameters, the factors,
  realism under both rules, the largest spike counts, the named odors' response fractions, the share of KCs responding to
  no odor, `H.code_commits()` of the process, which must be clean); runs under outputs/E93_missing_glomeruli_input/B/.
- score (C): readings.csv (per input and unit, every target of E89's `targets89`, the 120 and the 66 megamat pairs, raw
  counts, the calibrated readout and E73's per-KC reading), spread.csv (the pattern part's amplitude and shape parts
  beside every pattern part, `H.apl_effect_spread_shape` with the target through `H.canonical_pairs`), apl_free.csv (E87's
  target without APL), coresponse.csv (E89's co-response change), line.csv (the groups' APL-on distances against E92's
  part A line), who_responds_runs.csv, gate.csv (C's gate on T), equal.csv, bounds.csv, verdict.json,
  comparison_table.md, the figure.
- D: the eight wirings at F's held points again with the stage matrices saved (E90's `read_stages`; each run's counts
  checked against part B's); `--part Dscore`: E90's claw readings and verdict by group on all claws, on the claws of the
  46 glomeruli with input and on those of the 36 (stages.csv, the verdict's 'D').
- E: the family fitted again with F (`E10_joint_fit.py --orn-input <F> -k 2 --aplpn-max 3 --rate-target ctrl_frac --loss
  apl --apl-readout both --tag _cap3`, E_MAXFEV evaluations, from each wiring's own fit), on S (stage S), then on T (stage
  T) only if S's fitted population is E_S_GAIN or more below the reference against all recordings; refused unless the
  gate passes or F is an equal (`--force` overrides); `--part Escore`: fits.csv's rows scored, the population of eight by
  E91's part D rule.
- F: only on a pass by (C) or (E) (not written).
- variants: variants.csv and model_pair_distances_long.csv for the archive:
  `./venv/bin/python archive_tools/archive_model_outputs.py E93_missing_glomeruli_input`, then
  `./venv/bin/python archive_tools/migrate_model_output_params.py --apply`.
- table: comparison_table.md and the figure from what exists.

Outputs: results/E93_missing_glomeruli_input/; spike counts and responses under outputs/E93_missing_glomeruli_input/
(B/ for part B, D/ for the stage runs, E's fits at the folder's root as E10_joint_fit.py writes them).

Usage (from the repository root; one 16-thread process at a time on viki):
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part inputs
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part 0
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part A
    OMP_NUM_THREADS=16 OMP_THREAD_LIMIT=16 nice -n 10 ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part B
    ./venv/bin/python -u experiments/E93_missing_glomeruli_input.py --part score
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H
from experiments import E83_cross_wiring_fit as E83

EXPERIMENT = 'E93_missing_glomeruli_input'

# ---------------------------------------------------------------------------
# the inputs (the entry's Design; fixed before any run)

ADDED = ('VA4', 'VL1', 'DL4', 'VM4', 'DA4m', 'DA3', 'DA4l', 'VA1v', 'DL3', 'VA1d')
LAB_MEGAMAT = ('VA4', 'VL1')            # megamat_orn_deltas(), the 12 megamat odors
LAB_VALIDATION = ('DL4', 'VM4')         # validation2_orn_deltas(), paa, -aPine, euc, B-myr
HALLEM_ONLY = ('DA4m', 'DA3', 'DA4l', 'VA1v', 'DL3', 'VA1d')   # Hallem's 15 odors
HALLEM_MEGAMAT = ('DL4',)               # Hallem for the 12 megamat odors (the validation odors are the lab's)
HALLEM_NAMES = {'2h': '2-heptanone', 'IaA': 'isopentyl acetate', 'pa': 'pentyl acetate', '2-but': '2-butanone',
                'eb': 'ethyl butyrate', 'ep': 'ethyl propionate', '6al': 'hexanal', 't2h': 'E2-hexenal', '1-8ol': '1-octanol',
                '1-5ol': '1-pentanol', '1-6ol': '1-hexanol', 'ms': 'methyl salicylate', 'paa': 'phenylacetaldehyde',
                '-aPine': 'a-pinene', 'B-myr': 'b-myrcene'}
HALLEM_CONC = -2.0
INPUTS = ('I0', 'I0z', 'F', 'L', 'U', 'V')
CONTROL_KINDS = ('PL', 'PO')            # the ten rows' labels permuted; each row's values permuted over the odors
CONTROL_SEEDS = (1, 2)
NAMED_ODORS = ('2h', 'IaA', 'pa', '1-5ol', '1-6ol')
S_WIRINGS, T_WIRINGS, ALL_WIRINGS = E83.S_WIRINGS, E83.T_WIRINGS, E83.ALL_WIRINGS

# ---------------------------------------------------------------------------
# part A

SHARE_BINS = (-0.01, 0.0, 0.2, 0.4, 1.0)
SHARE_LABELS = ('none', 'up to 0.2', '0.2 to 0.4', 'above 0.4')
OLFDB_MEGAMAT, OLFDB_VALIDATION = 1, 4  # the lab's ORN datasets (olfdb), megamat and validation2
N_HALF_SPLITS = 200
DEP_FLOOR = 0.05                        # the departure's floor on the relative responses
DEP_NORMS = ('glom', 'odor', 'model_pn')
DEP_BOOT = 1000
NULL_DRAWS = 2000
A2_SUBSETS = ('meets', 'meets_read')    # E86's rule's entries; at 0.6
A2_MIN_NORMS = 2                        # the rule: in at least two of the three normalizations
DD_BOUND = 2.0                          # the extra arm's factors bounded to [1 / DD_BOUND, DD_BOUND]
SEED = 0

# ---------------------------------------------------------------------------
# part C, the gate (fixed in the entry)

REF_ALL8 = 0.0482                       # the reference's eight wirings (E83), pattern part, readout, all recordings
GATE_GAIN = 0.003                       # (i) T's population against all recordings and the mean over the four halves
GATE_CELLS = 12                         # (ii) of the 16 cells (T's four wirings by the four halves)
GATE_RAW_SLACK, GATE_MEAN_SLACK = 0.005, 0.01   # (iii)
GATE_SQ_ALL, GATE_SQ_SIX = 0.15, 0.075  # (iv)
CONTROL_RATIO = 0.5                     # (v)
CONTROL_DRAW2 = 0.0015                  # a second draw of each control if F moves T this much either way
BOUND_GAIN = 0.0015                     # a bound is preferred to F only this far below on S against one half
EQUAL_PATTERN, EQUAL_RAW, EQUAL_MEAN = 0.003, 0.005, 0.01
HALVES = ('A', 'B', 'C', 'D')
HALF_OTHER = {'A': 'B', 'B': 'A', 'C': 'D', 'D': 'C'}
FLIES_MEAN = 0.178

# ---------------------------------------------------------------------------
# part E

E_MAXFEV = 150
E_STAGES = {'S': S_WIRINGS, 'T': T_WIRINGS}
E_S_GAIN = 0.0015
E_PASS_GAIN = 0.003
E_WIRINGS_LOWER = 6
E_RAW_REF, E_RAW_SLACK, E_MEAN_SLACK = 0.073, 0.005, 0.01


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(part: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / part if part else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


_CACHE: Dict[object, object] = {}


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _commits_clean() -> dict:
    """`H.code_commits()` at the start of a part, which must be clean (E89's `_commits_clean`)."""
    from experiments import E89_apl_breadth_weighting as E89
    return E89._commits_clean()


def _wirings(arg: str) -> Tuple[str, ...]:
    return {'S': S_WIRINGS, 'T': T_WIRINGS, 'all8': ALL_WIRINGS}.get(arg) or tuple(x for x in arg.split(',') if x)


def _rms(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x ** 2)))


def _flagged(x) -> bool:
    return isinstance(x, str) and x.strip() not in ('', 'nan')


# ---------------------------------------------------------------------------
# the inputs

def lab_input() -> pd.DataFrame:
    """The lab's input (`H.apl16_orn_deltas()`: 36 glomeruli x (panel, odor))."""
    if 'lab' not in _CACHE:
        _CACHE['lab'] = H.apl16_orn_deltas().astype(float)
    return _CACHE['lab'].copy()


def _panel_cols(lab: pd.DataFrame, panel: str) -> list:
    return [c for c in lab.columns if c[0] == panel]


def lab_rows() -> pd.DataFrame:
    """VA4 and VL1 (megamat odors) and DL4 and VM4 (validation odors) from the lab's tables, in the lab input's columns;
    NaN where the panel did not measure the glomerulus."""
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas
    lab = lab_input()
    out = pd.DataFrame(np.nan, index=pd.Index(LAB_MEGAMAT + LAB_VALIDATION, name=lab.index.name), columns=lab.columns)
    for rows, panel, src in ((LAB_MEGAMAT, 'megamat', megamat_orn_deltas()), (LAB_VALIDATION, 'validation2', validation2_orn_deltas())):
        cols = _panel_cols(lab, panel)
        out.loc[list(rows), cols] = src.loc[list(rows), cols].to_numpy(dtype=float)
    return out


def hallem_table() -> pd.DataFrame:
    """Hallem's spike deltas at 10^-2 (drosolf, `add_sfr=False`), the uniglomerular receptors' glomeruli x the 15 odors
    of HALLEM_NAMES (abbreviated)."""
    if 'hallem' not in _CACHE:
        from drosolf import orns
        d = orns.orns(add_sfr=False, drop_sfr=False, columns='glomerulus')
        d = d.loc[:, [g for g in d.columns if '+' not in g]]
        if d.columns.duplicated().any():
            raise RuntimeError(f'Hallem: glomeruli with several receptors {sorted(d.columns[d.columns.duplicated()])}')
        out = pd.DataFrame({a: d.loc[name].astype(float) for a, name in HALLEM_NAMES.items()})
        out.index.name = 'glomerulus'
        _CACHE['hallem'] = out
    return _CACHE['hallem'].copy()


def slope_through_origin(x, y) -> float:
    """The least-squares slope of y on x through the origin (sum x y / sum x x), over the entries where both are
    finite."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    den = float(x[ok] @ x[ok])
    return float(x[ok] @ y[ok]) / den if den > 0 else np.nan


def hallem_scaling(lab: Optional[pd.DataFrame] = None, hal: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per odor of Hallem's 15 (and 'pooled'): the slope through the origin of the lab's input on Hallem's over the
    glomeruli both have, Pearson r, the number of glomeruli and their largest values (hallem_scaling.csv). `lab`:
    glomeruli x abbreviated odors (default the lab's input); `hal`: `hallem_table()`'s layout."""
    lab = H.abbrev_columns(lab_input()) if lab is None else lab
    hal = hallem_table() if hal is None else hal
    shared = [g for g in hal.index if g in lab.index]
    rows, xs, ys = [], [], []
    for o in hal.columns:
        x, y = hal.loc[shared, o].to_numpy(dtype=float), lab.loc[shared, o].to_numpy(dtype=float)
        xs.append(x)
        ys.append(y)
        rows.append(dict(odor=o, slope=slope_through_origin(x, y), r=float(np.corrcoef(x, y)[0, 1]), n_glomeruli=len(shared),
                         hallem_max=float(x.max()), lab_max=float(y.max()), glomeruli=' '.join(shared)))
    x, y = np.concatenate(xs), np.concatenate(ys)
    rows.append(dict(odor='pooled', slope=slope_through_origin(x, y), r=float(np.corrcoef(x, y)[0, 1]),
                     n_glomeruli=len(shared), hallem_max=float(x.max()), lab_max=float(y.max()), glomeruli=' '.join(shared)))
    return pd.DataFrame(rows).set_index('odor')


def added_rows(kind: str = 'F') -> pd.DataFrame:
    """The ten added rows (ADDED x the lab input's columns): 'F' (the lab's rows and Hallem's scaled per odor), 'U'
    (Hallem's as recorded), 'L' (the lab's four rows only), 'V' (VA4's row only), 'zero' (all 0); unmeasured entries
    0."""
    if kind not in ('F', 'U', 'L', 'V', 'zero'):
        raise ValueError(kind)
    lab = lab_input()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    rows = pd.DataFrame(0.0, index=pd.Index(ADDED, name=lab.index.name), columns=lab.columns)
    if kind == 'zero':
        return rows
    lr = lab_rows().fillna(0.0)
    for g in (('VA4',) if kind == 'V' else LAB_MEGAMAT + LAB_VALIDATION):
        rows.loc[g] = lr.loc[g].to_numpy(dtype=float)
    if kind in ('V', 'L'):
        return rows
    hal = hallem_table()
    slope = hallem_scaling(hal=hal)['slope']
    for g in HALLEM_ONLY + HALLEM_MEGAMAT:
        for j, o in enumerate(ab):
            if o not in hal.columns or (g in HALLEM_MEGAMAT and o not in H.MEGAMAT12):
                continue
            rows.iloc[ADDED.index(g), j] = float(hal.loc[g, o]) * (float(slope[o]) if kind == 'F' else 1.0)
    return rows


def measured_mask() -> pd.DataFrame:
    """Which of F's 160 added entries are measured (bool, ADDED x the lab input's columns; 134 True)."""
    lab = lab_input()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    m = pd.DataFrame(False, index=pd.Index(ADDED, name=lab.index.name), columns=lab.columns)
    lr = lab_rows()
    for g in LAB_MEGAMAT + LAB_VALIDATION:
        m.loc[g] = lr.loc[g].notna().to_numpy()
    hal = set(hallem_table().columns)
    for j, o in enumerate(ab):
        for g in HALLEM_ONLY:
            m.iloc[ADDED.index(g), j] = o in hal
        if o in H.MEGAMAT12:
            m.iloc[ADDED.index('DL4'), j] = m.iloc[ADDED.index('DL4'), j] or o in hal
    return m


def entry_sources(kind: str = 'F') -> pd.DataFrame:
    """The source of every added entry of `kind` (ADDED x abbreviated odors, strings)."""
    lab = lab_input()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    mm = measured_mask()
    out = pd.DataFrame('unmeasured (0)', index=list(ADDED), columns=ab)
    hal = set(hallem_table().columns)
    for i, g in enumerate(ADDED):
        for j, o in enumerate(ab):
            if g in LAB_MEGAMAT and o in H.MEGAMAT12 and mm.iloc[i, j]:
                src = 'the lab\'s megamat ORN recordings'
            elif g in LAB_VALIDATION and o in H.VALIDATION2_4 and mm.iloc[i, j]:
                src = 'the lab\'s validation ORN recordings'
            elif (g in HALLEM_ONLY or (g in HALLEM_MEGAMAT and o in H.MEGAMAT12)) and o in hal:
                src = {'F': 'Hallem x the odor\'s slope', 'U': 'Hallem as recorded (10^-2)'}.get(kind, 'not in this input (0)')
            else:
                src = 'unmeasured (0)'
            if kind == 'V' and g != 'VA4' and src != 'unmeasured (0)':
                src = 'not in this input (0)'
            out.iloc[i, j] = src
    return out


def full_input(rows: pd.DataFrame, lab: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """The lab's 36 rows followed by `rows` (the lab input's columns)."""
    lab = lab_input() if lab is None else lab
    if not rows.columns.equals(lab.columns):
        raise ValueError("the added rows' columns are not the lab input's")
    if set(rows.index) & set(lab.index):
        raise ValueError(f'rows already in the lab input: {sorted(set(rows.index) & set(lab.index))}')
    out = pd.concat([lab.astype(float), rows.astype(float)])
    out.index.name = lab.index.name
    return out


def permute_labels(rows: pd.DataFrame, seed: int) -> pd.DataFrame:
    """The rows' labels permuted among themselves with no row left in place (the first derangement drawn with `seed`):
    row i gets the values of row p[i]."""
    rng = np.random.default_rng(seed)
    n = len(rows)
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            break
    return pd.DataFrame(rows.to_numpy()[p], index=rows.index, columns=rows.columns)


def permute_odors(rows: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Each row's values permuted over the odors (one generator `seed`, rows in order, the megamat odors first), the
    megamat and the validation odors separately."""
    rng = np.random.default_rng(seed)
    ab = [H.odor_abbrev(c) for c in rows.columns]
    blocks = [np.array([j for j, o in enumerate(ab) if o in H.MEGAMAT12]),
              np.array([j for j, o in enumerate(ab) if o in H.VALIDATION2_4])]
    x = rows.to_numpy(dtype=float).copy()
    for i in range(x.shape[0]):
        for b in blocks:
            x[i, b] = x[i, b[rng.permutation(len(b))]]
    return pd.DataFrame(x, index=rows.index, columns=rows.columns)


def control_names(draws: Sequence[int] = CONTROL_SEEDS) -> List[str]:
    return [f'{k}{s}' for s in draws for k in CONTROL_KINDS]


def input_of(name: str) -> pd.DataFrame:
    """An input by name: INPUTS, 'PL<seed>' / 'PO<seed>' (the controls) or 'DD' (part A2's arm)."""
    if name in _CACHE.get('inputs', {}):
        return _CACHE['inputs'][name].copy()
    if name == 'I0':
        orn = lab_input()
    elif name == 'I0z':
        orn = full_input(added_rows('zero'))
    elif name in ('F', 'L', 'U', 'V'):
        orn = full_input(added_rows(name))
    elif name[:2] in CONTROL_KINDS and name[2:].isdigit():
        f = permute_labels if name[:2] == 'PL' else permute_odors
        orn = full_input(f(added_rows('F'), int(name[2:])))
    elif name == 'DD':
        orn = dendrite_input()
    else:
        raise ValueError(f'unknown input {name!r}')
    _CACHE.setdefault('inputs', {})[name] = orn
    return orn.copy()


def input_path(name: str) -> Path:
    return out_dir() / 'inputs' / f'{name}.csv'


def write_input(name: str) -> dict:
    """inputs/<name>.csv (glomeruli x abbreviated odors) and its record (sha256 of the file and of the frame)."""
    import hashlib
    from experiments.E77_input_screen import frame_sha256
    orn = input_of(name)
    lay = H.abbrev_columns(orn)
    lay.index.name = 'glomerulus'
    p = input_path(name)
    p.parent.mkdir(exist_ok=True)
    lay.to_csv(p)
    rec = dict(file=str(p.relative_to(H.shared_dir())), file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
               frame_sha256=frame_sha256(orn), n_glomeruli=int(orn.shape[0]))
    shas = _read_json(p.parent / 'sha256.json')
    shas[name] = rec
    _write_json(p.parent / 'sha256.json', shas)
    return rec


def write_inputs(names: Optional[Sequence[str]] = None) -> dict:
    """inputs/<name>.csv with their sha256 (inputs/sha256.json), added_rows.csv (every added entry of F, U, L and V with
    its source and whether measured), hallem_scaling.csv."""
    out = out_dir()
    names = list(INPUTS) + control_names() if names is None else list(names)
    recs = {n: write_input(n) for n in names}
    lab = lab_input()
    ab = [H.odor_abbrev(c) for c in lab.columns]
    mm = measured_mask()
    sl = hallem_scaling()
    rows = []
    for kind in ('F', 'U', 'L', 'V'):
        r = added_rows(kind)
        src = entry_sources(kind)
        for i, g in enumerate(ADDED):
            for j, o in enumerate(ab):
                rows.append(dict(input=kind, glomerulus=g, odor=o, value=float(r.iloc[i, j]), measured=bool(mm.iloc[i, j]),
                                 source=src.iloc[i, j], hallem_slope=float(sl.loc[o, 'slope']) if o in sl.index else np.nan))
    pd.DataFrame(rows).to_csv(out / 'added_rows.csv', index=False)
    sl.to_csv(out / 'hallem_scaling.csv')
    print(f'E93 inputs: {", ".join(names)}; {int(mm.to_numpy().sum())} of {mm.size} added entries measured', flush=True)
    return recs


# ---------------------------------------------------------------------------
# part 0

def build93(label: str, orn: pd.DataFrame, name: str):
    """E85's `build` with this experiment's cache: ((simulate, info), input tag). Raises if `J.bouton_base` gives the
    input another wiring or other APL weights than the lab input's."""
    from experiments import E85_orn_fly_population as E85
    tag = E83.input_tag_for(orn, f'E93_{name}')
    return E85.build(label, orn, tag, cache_dir=local_dir()), tag


def identity(wirings: Sequence[str] = ALL_WIRINGS) -> dict:
    """Part 0: I0z (the lab's input with the ten rows written out as zeros) on each wiring, its breadth recomputed (the
    wiring and APL weights checked against the lab input's by `build93`), at the fit's parameters APL on and off: 0 count
    and response differences from the saved fit, and the breadth equal to the lab input's. Resumable per wiring;
    identity.json; raises after all wirings if anything differs (outcome (h))."""
    import gc
    path = out_dir() / 'identity.json'
    res = _read_json(path)
    res.update(_commits_clean())
    orn = input_of('I0z')
    res.update(input='I0z', n_input_rows=int(orn.shape[0]), input_record=write_input('I0z'))
    res.setdefault('wirings', {})
    for w in wirings:
        if w in res['wirings'] and 'count_diffs' in res['wirings'][w]:
            continue
        t0 = time.time()
        P = E83.wiring_parts(w)
        try:
            (sim, info), tag = build93(w, orn, 'I0z')
        except RuntimeError as e:
            res['wirings'][w] = dict(wiring_equal=False, error=str(e), count_diffs=-1, response_diffs=-1, breadth_equal=False)
            _write_json(path, res)
            print(f'E93 0 {w}: {e}', flush=True)
            continue
        b = info['breadth']
        beq = bool(b.index.equals(P['breadth'].index) and np.array_equal(b.to_numpy(dtype=float), P['breadth'].to_numpy(dtype=float)))
        on = sim(dict(P['p']))
        off = sim(dict(P['p'], wAPLKC_scale=0.0, wAPLPN_scale=0.0))
        ns, nr = E83._same({'on': on, 'off': off}, P['saved'])
        res['wirings'][w] = dict(wiring_equal=True, breadth_equal=beq, breadth_max_abs_diff=float(np.max(np.abs(
            b.to_numpy(dtype=float) - P['breadth'].reindex(b.index).to_numpy(dtype=float)))), count_diffs=ns,
            response_diffs=nr, input_tag=tag, seconds=round(time.time() - t0, 1))
        _write_json(path, res)
        print(f'E93 0 {w}: {res["wirings"][w]}', flush=True)
        del sim, info, on, off
        gc.collect()
    done = [res['wirings'][w] for w in ALL_WIRINGS if w in res['wirings']]
    res['all_exact'] = bool(len(done) == len(ALL_WIRINGS) and all(
        d['wiring_equal'] and d['breadth_equal'] and d['count_diffs'] == 0 and d['response_diffs'] == 0 for d in done))
    _write_json(path, res)
    bad = {w: d for w, d in res['wirings'].items() if not (d['wiring_equal'] and d['count_diffs'] == 0 and d['response_diffs'] == 0)}
    if bad:
        raise RuntimeError(f'E93 part 0: I0z does not reproduce the saved fits: {bad} (outcome (h); the experiment stops here)')
    return res


def identity_passed() -> bool:
    return bool(_read_json(out_dir() / 'identity.json').get('all_exact'))


# ---------------------------------------------------------------------------
# part A1

def unit_meta(label: str) -> dict:
    """A wiring's claws' glomeruli and KCs and its boutons' glomeruli: E90's units (this host's outputs or the
    archive), else from E83's `wiring_parts` (one simulation)."""
    from experiments import E90_model_stages as E90
    for d in (H.outputs_dir(E90.EXPERIMENT) / 'B', H.shared_dir() / 'data' / 'model_outputs' / E90.EXPERIMENT / 'B'):
        p = d / f'{label}_units.json'
        if p.exists():
            return {k: np.asarray(v) for k, v in json.loads(p.read_text()).items()}
    W = E83.wiring_parts(label)['W']
    return dict(bouton_glomerulus=np.asarray([str(x) for x in W.columns.get_level_values('glomerulus')]),
                claw_glomerulus=np.asarray([str(x) for x in E90.glomerulus_of_claws(W)]),
                kc_of_claw=np.asarray(W.index.get_level_values(H.KC_ID_LEVEL)).astype(int))


def claw_counts(label: str) -> pd.DataFrame:
    """KC x glomerulus claw counts of a wiring (all of the wiring's glomeruli; `unit_meta`)."""
    m = unit_meta(label)
    kc = pd.Index(np.asarray(m['kc_of_claw']).astype(int), name=H.KC_ID_LEVEL)
    gl = pd.Index(np.asarray(m['claw_glomerulus']).astype(str), name='glomerulus')
    M = pd.crosstab(kc, gl)
    return M.reindex(columns=sorted(set(np.asarray(m['bouton_glomerulus']).astype(str)) | set(M.columns)), fill_value=0)


def kc_classes(M: pd.DataFrame, lab_glomeruli: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Per KC: its claws, its share of claws on the zero-filled glomeruli (the wiring's glomeruli outside the lab
    input's 36) and on the ten, whether it has a claw on VA4 and on DA4l or DA3, and the share bins (SHARE_BINS)."""
    lab = set(lab_input().index if lab_glomeruli is None else lab_glomeruli)
    zero = [g for g in M.columns if g not in lab]
    ten = [g for g in ADDED if g in M.columns]
    n = M.sum(axis=1).astype(float)
    df = pd.DataFrame(dict(n_claws=n, share_zero=M[zero].sum(axis=1) / n, share_ten=M[ten].sum(axis=1) / n,
                           va4=(M['VA4'] > 0) if 'VA4' in M else False,
                           da4l_or_da3=(M.reindex(columns=['DA4l', 'DA3'], fill_value=0).sum(axis=1) > 0)), index=M.index)
    for k in ('zero', 'ten'):
        df[f'bin_{k}'] = pd.cut(df[f'share_{k}'], list(SHARE_BINS), labels=list(SHARE_LABELS)).astype(str)
    return df


def _kc_ids(df: pd.DataFrame) -> pd.Index:
    return pd.Index(df.index.get_level_values(H.KC_ID_LEVEL) if isinstance(df.index, pd.MultiIndex) else df.index).astype(int)


def kc_table(o: Dict[str, tuple], M: pd.DataFrame) -> pd.DataFrame:
    """Per KC of a run: `kc_classes`, and per condition the odors it responds to and its responses to NAMED_ODORS."""
    cls = kc_classes(M)
    parts = [cls]
    for c in ('on', 'off'):
        r = H.abbrev_columns(o[c][0].astype(float) > 0.5)
        r.index = _kc_ids(r)
        missing = r.index.difference(cls.index)
        if len(missing):
            raise RuntimeError(f'{len(missing)} KCs of the run have no claws in the wiring')
        r = r.reindex(cls.index)
        parts.append(pd.DataFrame({f'n_odors_{c}': r[list(H.APL16_ODORS)].sum(axis=1).astype(float),
                                   **{f'{o_}_{c}': r[o_].astype(float) for o_ in NAMED_ODORS}}, index=cls.index))
    return pd.concat(parts, axis=1)


def who_responds_rows(t: pd.DataFrame, **meta) -> List[dict]:
    """Rows of the response rate (over KCs and the 16 odors) and the odors per KC by class of a `kc_table` (the bins of
    the zero-filled share and of the ten's share; a claw on VA4 or not; on DA4l or DA3 or not), APL on and off, with the
    named odors' response rates."""
    rows = []
    n_odors = len(H.APL16_ODORS)
    for kind, col in (('zero_filled_share', 'bin_zero'), ('ten_share', 'bin_ten'), ('VA4_claw', 'va4'), ('DA4l_or_DA3_claw', 'da4l_or_da3')):
        for lab_, g in t.groupby(col, sort=True):
            row = dict(meta, kind=kind, cls=str(lab_), n_kcs=int(len(g)), mean_claws=float(g.n_claws.mean()))
            for c in ('on', 'off'):
                row[f'rate_{c}'] = float(g[f'n_odors_{c}'].mean() / n_odors)
                row[f'odors_per_kc_{c}'] = float(g[f'n_odors_{c}'].mean())
                row[f'silent_share_{c}'] = float((g[f'n_odors_{c}'] == 0).mean())
                for o_ in NAMED_ODORS:
                    row[f'rate_{o_}_{c}'] = float(g[f'{o_}_{c}'].mean())
            rows.append(row)
    return rows


def who_responds(runs: Dict[str, Dict[str, tuple]], **meta) -> List[dict]:
    """`who_responds_rows` per wiring of {wiring: run} and on the KCs of all of them pooled ('pooled')."""
    rows, tabs = [], []
    for w, o in runs.items():
        t = kc_table(o, claw_counts(w))
        tabs.append(t)
        rows += who_responds_rows(t, unit=w, **meta)
    if len(tabs) > 1:
        rows += who_responds_rows(pd.concat(tabs, ignore_index=True), unit='pooled', **meta)
    return rows


def input_shares(orn: pd.DataFrame, base: pd.DataFrame, claws: pd.Series) -> pd.Series:
    """Per odor the added share of claw-weighted positive input: (sum over glomeruli of claws x max(orn, 0) minus the
    same for `base`) over the former (columns abbreviated)."""
    def tot(x):
        x = H.abbrev_columns(x).clip(lower=0.0).reindex(claws.index).fillna(0.0)
        return x.mul(claws.astype(float), axis=0).sum()
    a, b = tot(orn), tot(base)
    return (a - b) / a


def weighted_distances(orn: pd.DataFrame, claws: Optional[pd.Series] = None, odors: Sequence[str] = tuple(H.MEGAMAT12)
                       ) -> pd.Series:
    """1 - the claw-weighted (`claws`; None: unweighted) Pearson correlation over the glomeruli with input for `odors`
    (a nonzero entry among them; default the 12 megamat odors), per pair of `odors` (canonical order)."""
    from itertools import combinations
    X = H.abbrev_columns(orn).astype(float)[list(odors)]
    X = X.loc[(X != 0).any(axis=1)]
    w = np.ones(len(X)) if claws is None else claws.reindex(X.index).fillna(0.0).to_numpy(dtype=float)
    w = w / w.sum()
    out = {}
    for a, b in combinations(X.columns, 2):
        x, y = X[a].to_numpy(), X[b].to_numpy()
        mx, my = w @ x, w @ y
        r = (w * (x - mx) * (y - my)).sum() / np.sqrt((w * (x - mx) ** 2).sum() * (w * (y - my) ** 2).sum())
        out[(a, b)] = 1.0 - float(r)
    return E83._cs(pd.Series(out))


def half_agreement(values: np.ndarray, n_splits: int = N_HALF_SPLITS, seed: int = SEED) -> float:
    """The mean over random splits of the flies in two halves of Pearson r between the halves' mean profiles (flies x
    odors, NaN where not measured; odors measured in both halves)."""
    v = np.asarray(values, dtype=float)
    n = v.shape[0]
    if n < 2:
        return np.nan
    rng = np.random.default_rng(seed)
    rs = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for _ in range(n_splits):
            p = rng.permutation(n)
            a, b = np.nanmean(v[p[:n // 2]], axis=0), np.nanmean(v[p[n // 2:]], axis=0)
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() >= 3 and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
                rs.append(float(np.corrcoef(a[ok], b[ok])[0, 1]))
    return float(np.mean(rs)) if rs else np.nan


def olfdb_long(ds: int) -> pd.DataFrame:
    """An olfdb dataset's long table, cached under outputs/E93_missing_glomeruli_input/ (fetched from http://hal:8080)."""
    from data_prep import orn16_distances_from_olfdb as O
    p = local_dir() / f'olfdb_dataset{ds}.parquet'
    if not p.exists():
        df, meta = O.fetch(ds)
        df.to_parquet(p)
        _write_json(p.with_suffix('.json'), {k: meta.get(k) for k in ('id', 'title', 'data_sha256', 'n_flies', 'created_at')})
    return pd.read_parquet(p)


def fly_profiles(df: pd.DataFrame, panel: str, roi: str, odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Flies x odors trial-mean responses (dF/F) of one ROI in a panel; `odors` abbreviations at the model input's
    concentrations (None: every odor of the panel at its concentration, solvent mixtures excluded)."""
    d = df[(df.panel == panel) & (df.roi == roi)]
    if odors is not None:
        want = {f'{o} @ {H.ORN_INPUT_CONC[o]:g}': o for o in odors}
        d = d[d.odor1.isin(want)]
    d = d.assign(fly=[f'{a}/{b}' for a, b in zip(d.date, d.fly_num)])
    if 'odor2' in d.columns:
        d = d[d.odor2.isna() | (d.odor2.astype(str).isin(['solvent', 'nan', '']))]
    return d.groupby(['fly', 'odor1']).value.mean().unstack('odor1')


def a1_reliability() -> pd.DataFrame:
    """a1_reliability.csv: per glomerulus of the lab's four rows, the agreement between halves of the ORN flies (olfdb
    dataset 1 for VA4 and VL1, the 12 megamat odors; dataset 4 for DL4 and VM4, the validation panel's odors and the
    model's four), the al_analysis tables' (pebbled ORNs, GH146 PN dendrites; megamat), and VA4 and VL1 in the ORNs
    against the PN dendrites (Pearson r of the fly means over the 12 odors)."""
    from experiments import E90_model_stages as E90
    rows = []
    for ds, panel, gloms, odor_sets in ((OLFDB_MEGAMAT, 'megamat', LAB_MEGAMAT, (('model odors', list(H.MEGAMAT12)),)),
                                        (OLFDB_VALIDATION, 'validation2', LAB_VALIDATION,
                                         (('model odors', list(H.VALIDATION2_4)), ('panel odors', None)))):
        try:
            df = olfdb_long(ds)
        except Exception as e:      # noqa: BLE001 (olfdb unreachable: recorded, not fatal)
            rows += [dict(glomerulus=g, source=f'olfdb dataset {ds}', odors=k, error=repr(e)) for g in gloms for k, _ in odor_sets]
            continue
        for g in gloms:
            for k, ods in odor_sets:
                x = fly_profiles(df, panel, g, ods)
                rows.append(dict(glomerulus=g, source=f'olfdb dataset {ds} ({panel})', odors=k, n_flies=int(len(x)),
                                 n_odors=int(x.shape[1]), half_r=half_agreement(x.to_numpy())))
    al = {k: E90.al_table(k) for k in ('pebbled', 'GH146')}
    for g in LAB_MEGAMAT:
        prof = {}
        for k, t in al.items():
            if g not in t.columns:
                rows.append(dict(glomerulus=g, source=f'al_analysis {k}', odors='model odors', n_flies=0))
                continue
            x = t[g].unstack('odor').reindex(columns=list(H.MEGAMAT12)).dropna(how='all')
            prof[k] = x.mean()
            rows.append(dict(glomerulus=g, source=f'al_analysis {k}', odors='model odors', n_flies=int(len(x)),
                             n_odors=int(x.shape[1]), half_r=half_agreement(x.to_numpy())))
        if len(prof) == 2:
            a, b = prof['pebbled'], prof['GH146']
            ok = a.notna() & b.notna()
            rows.append(dict(glomerulus=g, source='ORNs (pebbled) against PN dendrites (GH146), fly means', odors='model odors',
                             n_odors=int(ok.sum()), r_orn_pn=float(np.corrcoef(a[ok], b[ok])[0, 1])))
    return pd.DataFrame(rows)


def part_a1() -> dict:
    """A1: a1_reliability.csv, hallem_scaling.csv, added_rows.csv (via `write_inputs`), input_shares.csv,
    input_distances.csv, who_responds.csv (the saved runs of the eight wirings)."""
    from experiments import E90_model_stages as E90
    out = out_dir()
    write_inputs()
    rel = a1_reliability()
    rel.to_csv(out / 'a1_reliability.csv', index=False)
    M = claw_counts('real')
    claws = M.sum(axis=0)
    lab = lab_input()
    sh = pd.DataFrame({n: input_shares(input_of(n), lab, claws) for n in ('F', 'U', 'L', 'V')})
    sh.index.name = 'odor'
    sh.to_csv(out / 'input_shares.csv')
    drows = []
    for n in ('I0', 'F', 'U', 'L', 'V'):
        for wn, wt in (('claw-weighted', claws), ('unweighted', None)):
            gm = E90.group_means(weighted_distances(input_of(n), wt))
            drows.append(dict(input=n, weighting=wn, **gm))
    dist = pd.DataFrame(drows)
    dist.to_csv(out / 'input_distances.csv', index=False)
    wr = pd.DataFrame(who_responds({w: E83.saved_fit(w)['saved'] for w in ALL_WIRINGS}, input='I0'))
    wr.to_csv(out / 'who_responds.csv', index=False)
    cls = kc_classes(M)
    res = dict(**H.code_commits(), claws_real=dict(total=int(M.to_numpy().sum()), on_ten=int(M[[g for g in ADDED if g in M]].to_numpy().sum()),
               zero_filled=int(M[[g for g in M.columns if g not in set(lab.index)]].to_numpy().sum())),
               kcs_real=dict(n=int(len(cls)), with_zero_filled_claw=float((cls.share_zero > 0).mean()),
                             with_claw_on_ten=float((cls.share_ten > 0).mean()), with_VA4=float(cls.va4.mean()),
                             with_DA4l_or_DA3=float(cls.da4l_or_da3.mean())),
               added_share_max=float(sh['F'].max()), n_measured=int(measured_mask().to_numpy().sum()))
    _write_json(out / 'part_a1.json', res)
    print(json.dumps({k: res[k] for k in ('claws_real', 'kcs_real')}, indent=1), flush=True)
    return res


# ---------------------------------------------------------------------------
# part A2

def departure(P: pd.DataFrame, O: pd.DataFrame, norm: str, floor: float = DEP_FLOOR) -> pd.DataFrame:
    """log2 of the dendrites' relative response over the ORNs' (odors x glomeruli, both clipped at 0): relative within a
    glomerulus over the odors ('glom'), within an odor over the glomeruli ('odor'), or over the whole table ('global'),
    with `floor` added to both. The model's PN stage in place of the ORNs is 'glom' with `O` the model's PN rate
    changes."""
    P, O = P.clip(lower=0.0), O.clip(lower=0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        if norm == 'glom':
            Pn, On = P / P.mean(axis=0), O / O.mean(axis=0)
        elif norm == 'odor':
            Pn, On = P.div(P.mean(axis=1), axis=0), O.div(O.mean(axis=1), axis=0)
        elif norm == 'global':
            Pn, On = P / float(np.nanmean(P.to_numpy())), O / float(np.nanmean(O.to_numpy()))
        else:
            raise ValueError(norm)
        return np.log2((Pn + floor) / (On + floor))


def sign_agreement(dep: np.ndarray, factor: np.ndarray) -> float:
    """The share of entries whose departure and searched log2 factor have the same sign (zeros and NaN count as
    disagreement); NaN without entries."""
    d, f = np.asarray(dep, dtype=float), np.asarray(factor, dtype=float)
    if d.size == 0:
        return np.nan
    return float(np.mean((np.sign(d) == np.sign(f)) & (np.sign(d) != 0)))


def null_agreement(dep: pd.Series, odors: Sequence[str], factor: np.ndarray, mask: np.ndarray, n: int = NULL_DRAWS,
                   seed: int = SEED) -> np.ndarray:
    """The sign agreement over `mask`'s entries with the departures (`dep`, one per entry of the list, aligned with
    `odors` and `factor`) permuted within each odor over the list's entries, `n` draws."""
    rng = np.random.default_rng(seed)
    d = np.asarray(dep, dtype=float)
    od = np.asarray(odors)
    m = np.asarray(mask, dtype=bool)
    f = np.asarray(factor, dtype=float)
    groups = [np.flatnonzero(od == o) for o in pd.unique(od)]
    out = np.empty(n)
    for k in range(n):
        dp = d.copy()
        for g in groups:
            dp[g] = d[g[rng.permutation(len(g))]]
        out[k] = sign_agreement(dp[m], f[m])
    return out


def a2_rule(results: Sequence[dict]) -> str:
    """The entry's rule over the normalizations' results (dicts with 'agree', 'null_p5', 'null_p95'): "in the list's
    direction" if the agreement is above the null's 95th percentile in at least A2_MIN_NORMS, "opposite" if below its
    5th in at least A2_MIN_NORMS, "neither" otherwise."""
    up = sum(1 for r in results if np.isfinite(r['agree']) and r['agree'] > r['null_p95'])
    down = sum(1 for r in results if np.isfinite(r['agree']) and r['agree'] < r['null_p5'])
    if up >= A2_MIN_NORMS:
        return "in the list's direction"
    if down >= A2_MIN_NORMS:
        return 'opposite'
    return 'neither'


def _fly_mean(t: pd.DataFrame, gloms: Sequence[str], flies: Optional[Sequence[str]] = None) -> pd.DataFrame:
    x = t if flies is None else t[t.index.get_level_values('fly').isin(list(flies))]
    return x.groupby(level='odor').mean().reindex(list(H.MEGAMAT12))[list(gloms)]


def model_pn_changes() -> pd.DataFrame:
    """The model's PN stage with APL on (E90's saved PN rates, window minus pre-odor mean), the mean over the wirings
    E90 saved, odors x glomeruli (megamat odors)."""
    from experiments import E90_model_stages as E90
    xs = []
    for w in ALL_WIRINGS:
        st = E90.load_stage('B', f'{w}_on')
        if st is not None:
            xs.append(E90.stage_units(st, 'pn')[list(H.MEGAMAT12)])
    if not xs:
        raise RuntimeError("E90's PN rates are neither in this host's outputs nor in the archive")
    return (sum(xs) / len(xs)).T


def a2_entries() -> pd.DataFrame:
    """E86's part E list on the glomeruli the GH146, pebbled and model tables share (E90's `shared_glomeruli`)."""
    from experiments import E90_model_stages as E90
    e = pd.read_csv(H.results_dir('E86_reduced_input_candidate') / 'partE_entries.csv')
    sh = set(E90.shared_glomeruli())
    e = e[e.glomerulus.isin(sh) & e.odor.isin(H.MEGAMAT12)].reset_index(drop=True)
    for c in A2_SUBSETS:
        e[c] = e[c].astype(str).str.lower().isin(['true', '1'])
    return e


def _dep_at(D: pd.DataFrame, e: pd.DataFrame) -> np.ndarray:
    return np.array([D.loc[o, g] for o, g in zip(e.odor, e.glomerulus)], dtype=float)


def part_a2(n_boot: int = DEP_BOOT, n_null: int = NULL_DRAWS, n_splits: int = N_HALF_SPLITS) -> dict:
    """A2 (the module docstring): dendrite_departures.csv (per entry and normalization: the departure, its 90 % interval
    over resampled flies, the share of splits whose halves agree in sign with it), dendrite_reading.json (per
    normalization and subset: entries, sign agreement, ρ, the null's mean, 5th and 95th percentiles and p; the halves'
    r; the rule)."""
    from experiments import E90_model_stages as E90
    out = out_dir()
    gh, pb = E90.al_table('GH146'), E90.al_table('pebbled')
    gl = sorted(set(E90.shared_glomeruli()))
    e = a2_entries()
    pn = model_pn_changes().reindex(columns=gl)
    P, O = _fly_mean(gh, gl), _fly_mean(pb, gl)
    fl_gh, fl_pb = list(gh.index.get_level_values('fly').unique()), list(pb.index.get_level_values('fly').unique())

    def dep_of(Pm, Om, norm):
        return departure(Pm, pn if norm == 'model_pn' else Om, 'glom' if norm == 'model_pn' else norm)

    rng = np.random.default_rng(SEED)
    rows, reading = [], dict(rule_inputs={}, subsets={}, n_entries=int(len(e)), glomeruli=gl)
    f = e.mean_log2.to_numpy(dtype=float)
    results_meets = []
    for norm in DEP_NORMS:
        D = dep_of(P, O, norm)
        d = _dep_at(D, e)
        boot = np.empty((n_boot, len(e)))
        for k in range(n_boot):
            bg = rng.choice(fl_gh, len(fl_gh))
            bp = rng.choice(fl_pb, len(fl_pb))
            Pb = pd.concat([gh.xs(x, level='fly') for x in bg]).groupby(level=0).mean().reindex(list(H.MEGAMAT12))[gl]
            Ob = pd.concat([pb.xs(x, level='fly') for x in bp]).groupby(level=0).mean().reindex(list(H.MEGAMAT12))[gl]
            boot[k] = _dep_at(dep_of(Pb, Ob, norm), e)
        lo, hi = np.nanpercentile(boot, [5, 95], axis=0)
        same, rs = np.zeros(len(e)), []
        for k in range(n_splits):
            ag, ap = rng.permutation(fl_gh), rng.permutation(fl_pb)
            d1 = _dep_at(dep_of(_fly_mean(gh, gl, ag[:len(ag) // 2]), _fly_mean(pb, gl, ap[:len(ap) // 2]), norm), e)
            d2 = _dep_at(dep_of(_fly_mean(gh, gl, ag[len(ag) // 2:]), _fly_mean(pb, gl, ap[len(ap) // 2:]), norm), e)
            same += (np.sign(d1) == np.sign(d)) & (np.sign(d2) == np.sign(d))
            ok = np.isfinite(d1) & np.isfinite(d2) & e.meets.to_numpy()
            if ok.sum() >= 3:
                rs.append(float(np.corrcoef(d1[ok], d2[ok])[0, 1]))
        share = same / n_splits
        for i, r in e.iterrows():
            rows.append(dict(odor=r.odor, glomerulus=r.glomerulus, norm=norm, mean_log2=float(r.mean_log2), meets=bool(r.meets),
                             meets_read=bool(r.meets_read), departure=float(d[i]), lo=float(lo[i]), hi=float(hi[i]),
                             halves_same_sign_share=float(share[i]), halves_agree=bool(share[i] >= 0.5)))
        reading['subsets'][norm] = {}
        for sub in A2_SUBSETS:
            m = e[sub].to_numpy()
            nul = null_agreement(pd.Series(d), e.odor.to_numpy(), f, m, n_null, SEED)
            agree = sign_agreement(d[m], f[m])
            ok = m & np.isfinite(d)
            rho = float(H.spearman(pd.Series(f[ok]), pd.Series(d[ok]))) if ok.sum() >= 3 else np.nan
            res = dict(n=int(m.sum()), agree=agree, rho=rho, null_mean=float(np.nanmean(nul)),
                       null_p5=float(np.nanpercentile(nul, 5)), null_p95=float(np.nanpercentile(nul, 95)),
                       p_ge=float(np.mean(nul >= agree)), p_le=float(np.mean(nul <= agree)))
            reading['subsets'][norm][sub] = res
            if sub == 'meets':
                results_meets.append(res)
        reading['subsets'][norm]['halves_r'] = float(np.mean(rs)) if rs else np.nan
    reading['rule'] = a2_rule(results_meets)
    reading['rule_text'] = ("in the list's direction: the sign agreement over the entries that meet E86's rule above the null's "
                            "95th percentile in at least two of the three normalizations; opposite: below its 5th in at least "
                            'two; neither otherwise')
    reading['dd_arm'] = reading['rule'] == "in the list's direction"
    reading.update(**H.code_commits(), settings=dict(DEP_FLOOR=DEP_FLOOR, DEP_NORMS=DEP_NORMS, DEP_BOOT=n_boot,
                   NULL_DRAWS=n_null, N_HALF_SPLITS=n_splits, SEED=SEED, model_pn="E90's saved PN rates, mean over wirings, "
                   'APL on, window minus pre-odor mean; the dendrites and the model normalized within a glomerulus'))
    pd.DataFrame(rows).to_csv(out / 'dendrite_departures.csv', index=False)
    _write_json(out / 'dendrite_reading.json', reading)
    print(f"E93 A2: {reading['rule']}; " + '; '.join(f"{n}: agree {reading['subsets'][n]['meets']['agree']:.2f} "
          f"(null {reading['subsets'][n]['meets']['null_p5']:.2f} to {reading['subsets'][n]['meets']['null_p95']:.2f})"
          for n in DEP_NORMS), flush=True)
    return reading


def dendrite_input() -> pd.DataFrame:
    """DD: F with the PN dendrites' departures (the 'glom' normalization) put on the input at the entries that meet E86's
    rule and whose halves agree in sign: the entry times 2 ** departure, bounded to [1 / DD_BOUND, DD_BOUND]. Only if
    part A2's rule is "in the list's direction"."""
    rd = _read_json(out_dir() / 'dendrite_reading.json')
    if not rd.get('dd_arm'):
        raise RuntimeError("the dendrite arm runs only if part A2's rule is \"in the list's direction\"")
    dep = pd.read_csv(out_dir() / 'dendrite_departures.csv')
    dep = dep[(dep.norm == 'glom') & dep.meets & dep.halves_agree]
    orn = input_of('F')
    ab = {H.odor_abbrev(c): c for c in orn.columns}
    lim = np.log2(DD_BOUND)
    for r in dep.itertuples():
        orn.loc[r.glomerulus, ab[r.odor]] *= float(2.0 ** np.clip(r.departure, -lim, lim))
    return orn


def part_a() -> dict:
    return dict(A1=part_a1(), A2=part_a2())


# ---------------------------------------------------------------------------
# part B

def held93(label: str, orn: pd.DataFrame, name: str) -> dict:
    """E85's `held_reading` (E83's reading with both response rates held) with this experiment's cache (the breadth's
    APL-off run under outputs/E93.../, by the input's tag)."""
    from experiments import E85_orn_fly_population as E85
    tag = E83.input_tag_for(orn, f'E93_{name}')
    return E85.held_reading(label, orn, tag, cache_dir=local_dir())


def run_tag(name: str, wiring: str) -> str:
    return f'B_{name}_{wiring}'


def save_run(tag: str, o: Dict[str, tuple], part: str = 'B') -> None:
    d = local_dir(part)
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(tag: str, part: str = 'B') -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}): this host's outputs, else the model-output archive; or None."""
    for d in (H.outputs_dir(EXPERIMENT) / part, H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / part):
        o = E83._load_pair(d, tag)
        if o is not None:
            return o
    return None


def run_outputs(name: str, wiring: str) -> Dict[str, tuple]:
    """A run of the readings: I0 is the wiring's saved fit; the others part B's."""
    if name == 'I0':
        return E83.saved_fit(wiring)['saved']
    o = load_run(run_tag(name, wiring))
    if o is None:
        raise RuntimeError(f'no saved run {run_tag(name, wiring)}')
    return o


def runs_b() -> pd.DataFrame:
    p = out_dir() / 'runs_B.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag', 'input', 'wiring'])


def _have_in(runs: pd.DataFrame, name: str, wirings: Sequence[str]) -> bool:
    """Whether `runs` (a table with 'input' and 'wiring') has the input on every one of `wirings` (I0 always)."""
    if name == 'I0':
        return True
    got = set(zip(runs['input'], runs['wiring']))
    return all((name, w) in got for w in wirings)


def _have(name: str, wirings: Sequence[str], rb: Optional[pd.DataFrame] = None) -> bool:
    return _have_in(runs_b() if rb is None else rb, name, wirings)


def _dist(name: str, wiring: str, readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances (`distances`), cached in this process (by input, wiring and readout) and in
    outputs/.../distances_cache.pkl (by the sha256 of the run's APL-on and APL-off spike counts and the readout, so a run
    made again is read again)."""
    import pickle
    from experiments.E77_input_screen import frame_sha256
    mem = _CACHE.setdefault('dist_mem', {})
    key = (name, wiring, readout)
    if key in mem:
        return mem[key]
    path = local_dir() / 'distances_cache.pkl'
    if 'dist_store' not in _CACHE:
        _CACHE['dist_store'] = pickle.loads(path.read_bytes()) if path.exists() else {}
    store = _CACHE['dist_store']
    o = run_outputs(name, wiring)
    skey = (frame_sha256(o['on'][1]), frame_sha256(o['off'][1]), readout)
    if skey not in store:
        store[skey] = distances(o, readout)
        try:
            path.write_bytes(pickle.dumps(store))
        except OSError:
            pass
    mem[key] = store[skey]
    return mem[key]


def pop_pattern(name: str, wirings: Sequence[str], target: str = 'z8', readout: str = 'cal') -> float:
    """The pattern part of a population's per-pair APL effect against a target (E83's `pattern_error`)."""
    d = E83.population([_dist(name, w, readout) for w in wirings])
    return _rms(E83.pattern_error(d, targets93()[target]).to_numpy())


def b_plan(rb: Optional[pd.DataFrame] = None) -> List[Tuple[str, str]]:
    """[(input, wiring)] of part B as far as it is decided now: F on all eight; the controls' first draws on T; their
    second draws if F moves T's population by CONTROL_DRAW2 or more either way (once F's T runs exist); L, U, V on S;
    a bound on T if the rule for fitted changes selects it on S (`bound_choices`, once the bounds' S runs exist); DD on
    all eight if part A2 says so."""
    rb = runs_b() if rb is None else rb
    plan = [('F', w) for w in ALL_WIRINGS]
    plan += [(c, w) for c in control_names(CONTROL_SEEDS[:1]) for w in T_WIRINGS]
    if _have('F', T_WIRINGS, rb):
        g = pop_pattern('I0', T_WIRINGS) - pop_pattern('F', T_WIRINGS)
        if abs(g) >= CONTROL_DRAW2:
            plan += [(c, w) for c in control_names(CONTROL_SEEDS[1:]) for w in T_WIRINGS]
    plan += [(b, w) for b in ('L', 'U', 'V') for w in S_WIRINGS]
    if all(_have(b, S_WIRINGS, rb) for b in ('L', 'U', 'V')) and _have('F', S_WIRINGS, rb):
        for b in sorted({c['bound'] for c in bound_choices() if c['to_T']}):
            plan += [(b, w) for w in T_WIRINGS]
    if _read_json(out_dir() / 'dendrite_reading.json').get('dd_arm'):
        plan += [('DD', w) for w in ALL_WIRINGS]
    return list(dict.fromkeys(plan))


def run_extras(o: Dict[str, tuple]) -> dict:
    """The named odors' response fractions and the share of KCs responding to no odor, APL on and off."""
    row = {}
    for c in ('on', 'off'):
        r = H.abbrev_columns(o[c][0].astype(float) > 0.5)
        for od in NAMED_ODORS:
            row[f'frac_{od}_{c}'] = float(r[od].mean())
        row[f'silent_share_{c}'] = float((~r[list(H.APL16_ODORS)].any(axis=1)).mean())
    return row


def _provenance_add(part: str, tag: str, prov: dict) -> None:
    path = out_dir() / f'provenance_{part}.json'
    d = _read_json(path) or dict(part=part, script='experiments/E93_missing_glomeruli_input.py')
    # (the file may already hold `provenance`'s dated records, written at the part's start, without 'points')
    d.setdefault('points', {})[tag] = prov
    _write_json(path, d)


def part_b(inputs: Optional[Sequence[str]] = None, wirings: Optional[Sequence[str]] = None, max_points: Optional[int] = None
           ) -> pd.DataFrame:
    """Part B (the module docstring): the plan's points in order (re-planned after each point, so the conditional ones
    follow from the runs before them), restricted to `inputs` / `wirings` when given; resumable (runs_B.csv, runs under
    outputs/.../B/, provenance per point in provenance_B.json). Needs part 0 to have passed."""
    import gc
    if not identity_passed():
        raise RuntimeError('run --part 0 first (the identity must hold before any reading)')
    cc = _commits_clean()
    path = out_dir() / 'runs_B.csv'
    n = 0
    while True:
        rb = runs_b()
        done = set(rb['tag'])
        todo = [(nm, w) for nm, w in b_plan(rb) if run_tag(nm, w) not in done
                and (inputs is None or nm in inputs) and (wirings is None or w in wirings)]
        if not todo or (max_points is not None and n >= max_points):
            break
        name, w = todo[0]
        tag = run_tag(name, w)
        orn = input_of(name)
        rec = write_input(name)
        res = held93(w, orn, name)
        o = {'on': res['on'], 'off': res['off']}
        save_run(tag, o)
        row = E83._run_row('B', name, w, tag, res, o, input=name, input_file_sha256=rec['file_sha256'],
                           input_frame_sha256=rec['frame_sha256'], input_tag=res['input_tag'], n_input_rows=int(orn.shape[0]),
                           orn_input=rec['file'], params_on=json.dumps({k: float(v) for k, v in res['params'].items()}),
                           params_off=json.dumps({k: float(v) for k, v in res['params_off'].items()}))
        row.update(run_extras(o))
        row.update(cc)
        E83._append_row(path, row)
        _provenance_add('B', tag, H.model_provenance(orn_deltas=orn))
        p0 = _rms(E83.pattern_error(E83.pair_distances(E83.saved_fit(w)['saved'], 'cal'), targets93()['z8']).to_numpy())
        p1 = _rms(E83.pattern_error(E83.pair_distances(o, 'cal'), targets93()['z8']).to_numpy())
        print(f"E93 B {tag}: thr x{res['thr_factor']:.4f}, APL x{res['apl_factor']:.4f} ({res['n_on']}/{res['target_on']} on, "
              f"{res['n_off']}/{res['target_off']} off); pattern part (readout, all recordings) {p0:.4f} -> {p1:.4f}; "
              f"flag '{row.get('realism_flag') or ''}'; {res['seconds']:.0f}s", flush=True)
        n += 1
        del res, o
        gc.collect()
    return runs_b()


# ---------------------------------------------------------------------------
# scoring

def targets93() -> Dict[str, pd.Series]:
    """E89's `targets89`: halves A to D, all recordings, the six typical, the noise-free effect, baseline-noise units."""
    from experiments import E89_apl_breadth_weighting as E89
    return E89.targets89()


def zcell_distances(o: Dict[str, tuple]) -> Dict[str, pd.Series]:
    """E73's per-KC reading (the session z-score, saturating map, calibrated on the run's APL-on counts by E74's
    `_calibrate_zcell`), graded, noise in expectation: {'on', 'off'} pair distances (canonical order)."""
    from experiments import E74_apl_effect_by_odors_per_kc as E74
    on = H.abbrev_columns(o['on'][1].astype(float))[list(H.APL16_ODORS)].to_numpy()
    off = H.abbrev_columns(o['off'][1].astype(float))[list(H.APL16_ODORS)].to_numpy()
    cal = E74._calibrate_zcell(on)
    kw = dict(q=float(cal['q']), noise=float(cal['noise']), amp='saturating', K=float(cal['K']),
              noise_ratio=float(cal['noise_ratio']), drift=float(cal['drift']))
    return {c: E83._cs(H.zcell_pair_distances(H.session_zscore(x, **kw), 'graded')) for c, x in (('on', on), ('off', off))}


def distances(o: Dict[str, tuple], readout: str) -> Dict[str, pd.Series]:
    """A run's pair distances: 'raw', 'cal' (E83's `pair_distances`), 'zcell' (`zcell_distances`)."""
    if readout == 'zcell':
        return zcell_distances(o)
    return E83.pair_distances(o, readout)


def _megamat_mask(index: pd.MultiIndex) -> np.ndarray:
    a, b = index.get_level_values(0), index.get_level_values(1)
    return np.asarray(a.isin(H.MEGAMAT12) & b.isin(H.MEGAMAT12))


def unit_scores(d: Dict[str, pd.Series], tg: Dict[str, pd.Series]) -> dict:
    """E83's `effect_scores` plus, per target, the pattern part, the loss and ρ over the 66 megamat pairs; the mean
    distances, ρ on odor-odor distances with APL on (control flies) and off (APL output blocked (TNT)), and the APL effect
    of 1-5ol with 1-6ol (the named pair)."""
    from experiments import E85_orn_fly_population as E85
    row = E83.effect_scores(d, tg)
    on, off = E83._cs(d['on']), E83._cs(d['off'])
    me = E83._cs(H.apl_effect(on, off))
    mm = _megamat_mask(me.index)
    for t, de0 in tg.items():
        de = E83._cs(de0).reindex(me.index)
        a, b = me[mm], de[mm]
        row[f'{t}_pattern66'] = _rms(((a - a.mean()) - (b - b.mean())).to_numpy())
        row[f'{t}_loss66'] = _rms((a - b).to_numpy())
        row[f'{t}_rho66'] = H.spearman(a, b)
    dd = E85.data_distances()
    row.update(mean_on=float(on.mean()), mean_off=float(off.mean()), rho_dist_on=H.spearman(on, dd['on'].reindex(on.index)),
               rho_dist_off=H.spearman(off, dd['off'].reindex(off.index)),
               effect_15ol_16ol=float(me.get(tuple(sorted(('1-5ol', '1-6ol'))), np.nan)))
    return row


def spread_rows(d: Dict[str, pd.Series], tg: Dict[str, pd.Series], pairs: str = 'all', **meta) -> List[dict]:
    """Per target, the pattern part with its amplitude and shape parts, the spread ratio, Pearson r and ρ (on the 120
    pairs, or the 66 megamat pairs with `pairs='megamat'`; `H.apl_effect_spread_shape` with the target through
    `H.canonical_pairs`, as E91's `_spread_rows`; the parts checked to add up to the pattern part)."""
    rows = []
    on, off = E83._cs(d['on']), E83._cs(d['off'])
    if pairs == 'megamat':
        m = _megamat_mask(on.index)
        on, off = on[m], off[m]
    elif pairs != 'all':
        raise ValueError(pairs)
    me = E83._cs(H.apl_effect(on, off))
    for t, de0 in tg.items():
        de = E83._cs(de0).reindex(me.index)
        sp = H.apl_effect_spread_shape(on, off, H.canonical_pairs(de))
        err = (me - me.mean()) - (de - de.mean())
        pattern = _rms(err.to_numpy())
        if not np.isclose(pattern ** 2, sp['apl_effect_amplitude'] ** 2 + sp['apl_effect_shape'] ** 2, rtol=1e-9, atol=1e-12):
            raise RuntimeError(f'{t}: the amplitude and shape parts do not add up to the pattern part (pairs misaligned?)')
        rows.append(dict(meta, pairs=pairs, n_pairs=int(len(me)), target=t, pattern=pattern, amplitude=sp['apl_effect_amplitude'],
                         shape=sp['apl_effect_shape'], sd_model=sp['apl_effect_sd_model'], sd_data=sp['apl_effect_sd_data'],
                         sd_ratio=sp['apl_effect_sd_ratio'], pearson=float(np.corrcoef(me.to_numpy(), de.to_numpy())[0, 1]),
                         rho=H.spearman(me, de)))
    return rows


def runs_table() -> pd.DataFrame:
    """Every run of the readings: part B's (runs_B.csv) and I0's (the saved fits, with their realism), with 'flagged'."""
    rb = runs_b()
    rows = []
    for w in ALL_WIRINGS:
        P = E83.saved_fit(w)
        rows.append(dict(tag=f'I0_{w}', input='I0', wiring=w, at_bound_thr=False, at_bound_apl=False, thr_factor=1.0,
                         apl_factor=1.0, fit_folder=P['folder'], fit_tag=P['tag'], **E83.realism(P['saved']),
                         **run_extras(P['saved'])))
    df = pd.concat([pd.DataFrame(rows), rb], ignore_index=True, sort=False) if len(rb) else pd.DataFrame(rows)
    df = df.drop_duplicates(['input', 'wiring'], keep='last').reset_index(drop=True)
    df['flagged'] = [_flagged(x) for x in df.realism_flag]
    df['at_bound'] = df.at_bound_thr.fillna(False).astype(bool) | df.at_bound_apl.fillna(False).astype(bool)
    return df


UNITS = (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS))


def readings(readouts: Sequence[str] = ('cal', 'raw', 'zcell')) -> pd.DataFrame:
    """readings.csv (per input, unit and readout `unit_scores` against every target, the runs outside the bands and at a
    bracket's bound, and the change against I0's same unit and readout, 'd_<column>'), spread.csv (`spread_rows` under
    the readout and on raw counts, the 120 and the 66 megamat pairs) and pairs.csv (the per-pair APL effect per input,
    population and readout)."""
    out = out_dir()
    tg = targets93()
    runs = runs_table()
    rows, srows, prs = [], [], []
    for name, g in runs.groupby('input', sort=False):
        ws = list(g.wiring)
        meta = {r.wiring: r for r in g.itertuples()}
        units = [(w, (w,)) for w in ws] + [(u, wu) for u, wu in UNITS if all(w in ws for w in wu)]
        for ro in readouts:
            dd = {w: _dist(name, w, ro) for w in ws}
            for unit, wu in units:
                d = dd[unit] if len(wu) == 1 and unit in dd else E83.population([dd[w] for w in wu])
                kind = 'wiring' if len(wu) == 1 and unit in dd else 'population'
                rows.append(dict(input=name, unit=unit, kind=kind, n_wirings=len(wu), readout=ro,
                                 n_flagged=int(sum(meta[w].flagged for w in wu)), n_at_bound=int(sum(meta[w].at_bound for w in wu)),
                                 **unit_scores(d, tg)))
                if ro in ('cal', 'raw'):
                    for pairs in ('all', 'megamat'):
                        srows += spread_rows(d, tg, pairs, input=name, unit=unit, kind=kind, readout=ro)
                if kind == 'population':
                    me = E83._cs(H.apl_effect(d['on'], d['off']))
                    prs.append(pd.DataFrame(dict(input=name, unit=unit, readout=ro, odor_a=me.index.get_level_values(0),
                                                 odor_b=me.index.get_level_values(1), on=E83._cs(d['on']).reindex(me.index).to_numpy(),
                                                 off=E83._cs(d['off']).reindex(me.index).to_numpy(), effect=me.to_numpy())))
    rd = pd.DataFrame(rows)
    base = rd[rd.input == 'I0'].set_index(['unit', 'readout'])
    cols = [c for c in rd.columns if c.endswith(('_loss', '_pattern', '_rho', '_pattern66', '_loss66', '_rho66'))
            or c in ('mean_effect', 'mean_on', 'mean_off', 'rho_dist_on', 'rho_dist_off', 'effect_15ol_16ol')]
    for c in cols:
        rd[f'd_{c}'] = [v - base[c].get((u, ro), np.nan) for u, ro, v in zip(rd.unit, rd.readout, rd[c])]
    rd.to_csv(out / 'readings.csv', index=False)
    pd.DataFrame(srows).to_csv(out / 'spread.csv', index=False)
    if prs:
        pd.concat(prs, ignore_index=True).to_csv(out / 'pairs.csv', index=False)
    return rd


def _val(rd: pd.DataFrame, name: str, unit: str, col: str, readout: str = 'cal') -> float:
    r = rd[(rd.input == name) & (rd.unit == unit) & (rd.readout == readout)]
    return float(r[col].iloc[0]) if len(r) and col in r.columns else np.nan


def _flags(runs: pd.DataFrame, name: str, wirings: Sequence[str]) -> int:
    r = runs[(runs.input == name) & runs.wiring.isin(wirings)]
    return int(r.flagged.sum())


def gate93(rd: pd.DataFrame, runs: pd.DataFrame, name: str = 'F', ref: str = 'I0', flies_mean: Optional[float] = None) -> dict:
    """The gate (C), fixed in the entry, on T for input `name` against the reference (readout 'cal' unless named): (i)
    T's population gains GATE_GAIN or more against all recordings and as the mean over the four halves; (ii) at least
    GATE_CELLS of T's 16 wiring-by-half cells lower; (iii) every point of the input inside the four bands (its runs on
    all eight wirings), T's raw-count loss not above the reference's by more than GATE_RAW_SLACK, T's mean APL effect
    within GATE_MEAN_SLACK of the flies'; (iv) T's squared pattern error GATE_SQ_ALL lower against all recordings and
    GATE_SQ_SIX against the six typical recordings; (v) each control kind's mean gain over its draws (T, all recordings)
    less than CONTROL_RATIO times the input's; (vi) on the 66 megamat pairs alone T's population lower too. Beside: S's
    gains, the number of the eight wirings lower, the parts of the pattern part. `runs`: `runs_table`'s layout ('input',
    'wiring', 'flagged'); `flies_mean` the flies' mean APL effect (default the main target's)."""
    fm = float(targets93()['z8'].mean()) if flies_mean is None else float(flies_mean)
    out = dict(input=name, reference=ref, complete=bool(_have_in(runs, name, T_WIRINGS)))
    if not out['complete']:
        out.update(passes=False, note='T runs missing')
        return out
    hs = [f'half{h}' for h in HALVES]
    g = {t: _val(rd, ref, 'T', f'{t}_pattern') - _val(rd, name, 'T', f'{t}_pattern') for t in ('z8', 'z6') + tuple(hs)}
    out.update(pattern_T=_val(rd, name, 'T', 'z8_pattern'), pattern_T_ref=_val(rd, ref, 'T', 'z8_pattern'), gain_all=g['z8'],
               gain_halves_mean=float(np.mean([g[t] for t in hs])), **{f'gain_{t}': g[t] for t in hs})
    out['i'] = bool(out['gain_all'] >= GATE_GAIN and out['gain_halves_mean'] >= GATE_GAIN)
    cells = [(w, t) for w in T_WIRINGS for t in hs]
    out['n_cells_lower'] = int(sum(_val(rd, name, w, f'{t}_pattern') < _val(rd, ref, w, f'{t}_pattern') for w, t in cells))
    out['n_cells'] = len(cells)
    out['ii'] = bool(out['n_cells_lower'] >= GATE_CELLS)
    out['n_flagged_T'] = _flags(runs, name, T_WIRINGS)
    out['n_flagged_S'] = _flags(runs, name, S_WIRINGS)
    out['n_flagged'] = _flags(runs, name, ALL_WIRINGS)
    out['raw_loss_change'] = _val(rd, name, 'T', 'z8_loss', 'raw') - _val(rd, ref, 'T', 'z8_loss', 'raw')
    out['mean_effect_cal'] = _val(rd, name, 'T', 'mean_effect')
    out['flies_mean'] = fm
    out['iii'] = bool(out['n_flagged'] == 0 and out['raw_loss_change'] <= GATE_RAW_SLACK
                      and abs(out['mean_effect_cal'] - fm) <= GATE_MEAN_SLACK)
    out['sq_reduction_all'] = 1.0 - _val(rd, name, 'T', 'z8_pattern') ** 2 / _val(rd, ref, 'T', 'z8_pattern') ** 2
    out['sq_reduction_six'] = 1.0 - _val(rd, name, 'T', 'z6_pattern') ** 2 / _val(rd, ref, 'T', 'z6_pattern') ** 2
    out['iv'] = bool(out['sq_reduction_all'] >= GATE_SQ_ALL and out['sq_reduction_six'] >= GATE_SQ_SIX)
    ctl = {}
    for k in CONTROL_KINDS:
        gs = [_val(rd, ref, 'T', 'z8_pattern') - _val(rd, f'{k}{s}', 'T', 'z8_pattern') for s in CONTROL_SEEDS
              if _have_in(runs, f'{k}{s}', T_WIRINGS)]
        ctl[k] = dict(n_draws=len(gs), gains=gs, mean_gain=float(np.mean(gs)) if gs else np.nan)
    out['controls'] = ctl
    out['v'] = bool(out['gain_all'] > 0 and all(c['n_draws'] > 0 and c['mean_gain'] < CONTROL_RATIO * out['gain_all']
                                               for c in ctl.values()))
    out['gain_all_66'] = _val(rd, ref, 'T', 'z8_pattern66') - _val(rd, name, 'T', 'z8_pattern66')
    out['vi'] = bool(out['gain_all_66'] > 0)
    out['passes_i_to_v'] = bool(out['i'] and out['ii'] and out['iii'] and out['iv'] and out['v'])
    out['passes'] = bool(out['passes_i_to_v'] and out['vi'])
    out['passes_120_not_66'] = bool(out['passes_i_to_v'] and not out['vi'])
    # beside
    out['gain_S_all'] = _val(rd, ref, 'S', 'z8_pattern') - _val(rd, name, 'S', 'z8_pattern') if _have_in(runs, name, S_WIRINGS) else np.nan
    if _have_in(runs, name, ALL_WIRINGS):
        out['n_wirings_lower'] = int(sum(_val(rd, name, w, 'z8_pattern') < _val(rd, ref, w, 'z8_pattern') for w in ALL_WIRINGS))
        out['gain_all8'] = _val(rd, ref, 'all8', 'z8_pattern') - _val(rd, name, 'all8', 'z8_pattern')
    sp = _spread_of(name, 'T')
    sp0 = _spread_of(ref, 'T')
    out.update(amplitude_T=sp.get('amplitude'), shape_T=sp.get('shape'), sd_ratio_T=sp.get('sd_ratio'), pearson_T=sp.get('pearson'),
               amplitude_T_ref=sp0.get('amplitude'), shape_T_ref=sp0.get('shape'), sd_ratio_T_ref=sp0.get('sd_ratio'),
               pearson_T_ref=sp0.get('pearson'))
    return out


def _spread_of(name: str, unit: str, target: str = 'z8', readout: str = 'cal', pairs: str = 'all') -> dict:
    p = out_dir() / 'spread.csv'
    if not p.exists():
        return {}
    sp = pd.read_csv(p)
    r = sp[(sp.input == name) & (sp.unit == unit) & (sp.target == target) & (sp.readout == readout) & (sp.pairs == pairs)]
    return r.iloc[0].to_dict() if len(r) else {}


def equal93(rd: pd.DataFrame, runs: pd.DataFrame, name: str = 'F', ref: str = 'I0') -> dict:
    """The equal rule: the input within EQUAL_PATTERN of the reference on S and on T against all recordings and all four
    halves (readout), all its runs inside the bands, within EQUAL_RAW in raw-count loss and EQUAL_MEAN in mean APL effect
    (on S and on T)."""
    if not (_have_in(runs, name, S_WIRINGS) and _have_in(runs, name, T_WIRINGS)):
        return dict(input=name, equal=False, note='runs missing')
    d = {f'{u}_{t}': _val(rd, name, u, f'{t}_pattern') - _val(rd, ref, u, f'{t}_pattern')
         for u in ('S', 'T') for t in ('z8',) + tuple(f'half{h}' for h in HALVES)}
    raw = {u: _val(rd, name, u, 'z8_loss', 'raw') - _val(rd, ref, u, 'z8_loss', 'raw') for u in ('S', 'T')}
    me = {u: _val(rd, name, u, 'mean_effect') - _val(rd, ref, u, 'mean_effect') for u in ('S', 'T')}
    nf = _flags(runs, name, ALL_WIRINGS)
    eq = bool(all(abs(x) <= EQUAL_PATTERN for x in d.values()) and nf == 0 and all(abs(x) <= EQUAL_RAW for x in raw.values())
              and all(abs(x) <= EQUAL_MEAN for x in me.values()))
    return dict(input=name, reference=ref, **{f'd_pattern_{k}': v for k, v in d.items()},
                **{f'd_raw_loss_{u}': v for u, v in raw.items()}, **{f'd_mean_effect_{u}': v for u, v in me.items()},
                n_flagged=nf, max_abs_d_pattern=float(max(abs(x) for x in d.values())), equal=eq,
                better_everywhere=bool(all(x < 0 for x in d.values())))


def bound_choices(rd: Optional[pd.DataFrame] = None) -> List[dict]:
    """The rule for fitted changes on the bounds, per half X of the four: the bound (L, U, V) with the lowest pattern
    part on S against half X; it goes to T ('to_T') if it is BOUND_GAIN or more below F's on S against half X."""
    out = []
    for X in HALVES:
        t = f'half{X}'
        if rd is not None:
            p = {b: _val(rd, b, 'S', f'{t}_pattern') for b in ('L', 'U', 'V', 'F')}
        else:
            p = {b: pop_pattern(b, S_WIRINGS, t) for b in ('L', 'U', 'V', 'F')}
        b = min(('L', 'U', 'V'), key=lambda k: p[k])
        out.append(dict(half=X, other=HALF_OTHER[X], bound=b, pattern_S=p[b], gain_over_F=p['F'] - p[b],
                        to_T=bool(p['F'] - p[b] >= BOUND_GAIN), **{f'pattern_S_{k}': v for k, v in p.items()}))
    return out


def bounds93(rd: pd.DataFrame, runs: pd.DataFrame, flies_mean: Optional[float] = None) -> List[dict]:
    """bounds.csv: `bound_choices`, and for each choice that goes to T, the gate read on T against the other half
    (against I0: (i) the gain against half Y at least GATE_GAIN; (ii) at least GATE_CELLS / 16 of T's four wirings lower
    against Y; (iii) as the gate's; (iv) the squared pattern error GATE_SQ_ALL lower against Y and GATE_SQ_SIX against
    the six typical; (vi) lower on the 66 megamat pairs against Y), and whether T's pattern part against Y is below F's
    ('preferred_to_F', with the gate passed)."""
    fm = float(targets93()['z8'].mean()) if flies_mean is None else float(flies_mean)
    out = []
    for c in bound_choices(rd):
        c = dict(c)
        b, Y = c['bound'], f"half{c['other']}"
        if c['to_T'] and _have_in(runs, b, T_WIRINGS):
            g = _val(rd, 'I0', 'T', f'{Y}_pattern') - _val(rd, b, 'T', f'{Y}_pattern')
            lower = sum(_val(rd, b, w, f'{Y}_pattern') < _val(rd, 'I0', w, f'{Y}_pattern') for w in T_WIRINGS)
            raw = _val(rd, b, 'T', 'z8_loss', 'raw') - _val(rd, 'I0', 'T', 'z8_loss', 'raw')
            mean = _val(rd, b, 'T', 'mean_effect')
            sqY = 1.0 - _val(rd, b, 'T', f'{Y}_pattern') ** 2 / _val(rd, 'I0', 'T', f'{Y}_pattern') ** 2
            sq6 = 1.0 - _val(rd, b, 'T', 'z6_pattern') ** 2 / _val(rd, 'I0', 'T', 'z6_pattern') ** 2
            g66 = _val(rd, 'I0', 'T', f'{Y}_pattern66') - _val(rd, b, 'T', f'{Y}_pattern66')
            nf = _flags(runs, b, list(S_WIRINGS) + list(T_WIRINGS))
            c.update(read_on_T=True, gain_T_other=g, i=bool(g >= GATE_GAIN), n_T_lower=int(lower),
                     ii=bool(lower >= GATE_CELLS / 16 * len(T_WIRINGS)), raw_loss_change=raw, mean_effect_cal=mean, n_flagged=nf,
                     iii=bool(nf == 0 and raw <= GATE_RAW_SLACK and abs(mean - fm) <= GATE_MEAN_SLACK), sq_reduction_other=sqY,
                     sq_reduction_six=sq6, iv=bool(sqY >= GATE_SQ_ALL and sq6 >= GATE_SQ_SIX), gain_T_other_66=g66, vi=bool(g66 > 0))
            c['passes'] = bool(c['i'] and c['ii'] and c['iii'] and c['iv'] and c['vi'])
            c['preferred_to_F'] = bool(c['passes'] and _val(rd, b, 'T', f'{Y}_pattern') < _val(rd, 'F', 'T', f'{Y}_pattern'))
        else:
            c.update(read_on_T=False, passes=False, preferred_to_F=False)
        out.append(c)
    return out


def apl_free(units: Sequence[Tuple[str, Sequence[str]]] = UNITS) -> pd.DataFrame:
    """apl_free.csv: E87's target without APL (the effect of E87's chosen form within the recordings with APL output
    blocked (TNT), z, with the trial noise removed and as trial means; apl_free_target_pairs.csv) against each input's
    populations: E87's `ModelResponses` (raw counts and the readout) through the same form, the population's APL-off
    distances after the form minus its own (E87's `scores`)."""
    from experiments import E87_inhibition_form as E87
    f = E87.chosen_form(H.results_dir(E87.EXPERIMENT))
    t = pd.read_csv(H.results_dir(E87.EXPERIMENT) / 'apl_free_target_pairs.csv')
    tg = {}
    for rd_ in ('noise_free', 'trial_means'):
        s = t[(t.unit == 'z') & (t.reading == rd_) & (t.form_label == 'chosen')]
        tg[rd_] = E87._pairs_to_triu(pd.Series(s.effect.to_numpy(), index=pd.MultiIndex.from_arrays([s.odor_a, s.odor_b])))
    runs = runs_table()
    rows = []
    for name, g in runs.groupby('input', sort=False):
        ws = set(g.wiring)
        for ro in ('raw', 'cal'):
            mrs = {w: E87.ModelResponses(run_outputs(name, w), ro, label=w) for w in ws}
            for unit, wu in UNITS:
                if not all(w in ws for w in wu):
                    continue
                pred = E87._pop([mrs[w].form(f)[0] for w in wu]) - E87._pop([mrs[w].d_off for w in wu])
                for rd_, tv in tg.items():
                    rows.append(dict(input=name, unit=unit, readout=ro, reading=rd_, form=f.name, **E87.scores(pred, tv)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'apl_free.csv', index=False)
    return df


def line_reading(rd_pairs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """line.csv: each input's populations' ester and alcohol groups' mean APL-on distances relative to the mean over the
    120 pairs (E92's `relative_distances`), placed against E92's part A line (`E92.line(readout)['S']`, fitted through
    E89's 13 weightings as populations of S's four wirings, under each readout) by E92's `place_on_line`: the residual,
    signed toward each flies' point, and whether it is off the line toward the flies against halves A and B (beyond the
    weightings' largest residual; beyond their SD beside). On T and the eight, where the line was not fitted, the residual
    is taken relative to I0's on the same unit (E92's `ref_resid`)."""
    from experiments import E92_steps_like_for_like as E92
    pr = pd.read_csv(out_dir() / 'pairs.csv') if rd_pairs is None else rd_pairs
    pts = {}
    for (name, unit, ro), g in pr.groupby(['input', 'unit', 'readout'], sort=False):
        if ro not in ('raw', 'cal'):
            continue
        on = pd.Series(g.on.to_numpy(), index=pd.MultiIndex.from_arrays([g.odor_a, g.odor_b]))
        pts[(name, unit, ro)] = E92.relative_distances(on)
    rows = []
    for (name, unit, ro), (est, alc) in pts.items():
        L = E92.line(ro)['S']
        ref = 0.0
        if unit != 'S' and ('I0', unit, ro) in pts:
            e0, a0 = pts[('I0', unit, ro)]
            ref = float(e0 - (L['intercept'] + L['slope'] * a0))
        rows.append(dict(input=name, unit=unit, readout=ro, ester_rel=est, alcohol_rel=alc, slope=L['slope'],
                         intercept=L['intercept'], scatter=L['scatter'], max_abs_resid=L['max_abs_resid'], ref_resid=ref,
                         **E92.place_on_line(est, alc, L, ref)))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'line.csv', index=False)
    return df


def _co_job93(job: tuple) -> dict:
    """One run read by E88's readings 1 and 1b and by the geometric mean of E88's two estimates at L1 and L3 (E91's
    `_co_job91`)."""
    import os
    from experiments import E89_apl_breadth_weighting as E89
    name, wiring, K, cand, rel, gap = job
    warnings.simplefilter('ignore')
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(k, '1')
    o = run_outputs(name, wiring)
    mr = E89.ModelReading89(o, label=f'{name}_{wiring}', K=K, rel_target=rel, gap_target=gap)
    res = dict(input=name, wiring=wiring, readings={})
    for rdn in ('r1', 'r1b'):
        F = mr.mean_factors(rdn)
        res['readings'][rdn] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    for level in ('L1', 'L3'):
        v = {c: mr.level_values(level, c, None)[0] for c in E89.CONDS}
        F = {c: E89.cand_factors(v[c], cand) for c in E89.CONDS}
        res['readings'][f'cand_{level}'] = {c: dict(s=F[c].s, O=F[c].O, R=F[c].R) for c in E89.CONDS}
    return res


CO_CANDIDATE = 'gm'


def coresponse93(n_workers: int = 4, inputs: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """coresponse.csv: per input and population (S, T) the changes of L (co-response) and R (recruitment) under APL of the
    four nearest ester pairs, the ester and the alcohol groups (E88's `log_changes`) in readings 1 and 1b and with the
    geometric mean of E88's two estimates at L1 and L3, against the flies' (90 % intervals over resampled recordings;
    E88's B1 rule), as E91's `coresponse91` reads E89's `coresponse_c`. Jobs cached in coresponse_cache.pkl."""
    import multiprocessing
    import pickle
    from concurrent.futures import ProcessPoolExecutor
    from experiments import E88_responder_sets as E88
    from experiments import E89_apl_breadth_weighting as E89
    out = out_dir()
    runs = runs_table()
    if inputs is not None:
        runs = runs[runs.input.isin(list(inputs))]
    rel = E88.flies_reliability('z')
    gap = {c: v['thresholded'] for c, v in E89.flies_gap_target('z').items()}
    cache = local_dir() / 'coresponse_cache.pkl'
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    units = {n: [(u, wu) for u, wu in UNITS[:2] if set(wu) <= set(g.wiring)] for n, g in runs.groupby('input')}
    jobs = [(r.input, r.wiring, E83.saved_fit(r.wiring)['K'], CO_CANDIDATE, rel, gap) for r in runs.itertuples()
            if (r.input, r.wiring) not in done and any(r.wiring in wu for _, wu in units[r.input])]
    if jobs:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=multiprocessing.get_context('spawn')) as ex:
            for res, j in zip(ex.map(_co_job93, jobs), jobs):
                done[(j[0], j[1])] = res
                cache.write_bytes(pickle.dumps(done))
                print(f'E93 coresponse {j[0]} {j[1]}', flush=True)
    fl = {}
    fv = E89.flies_values({})
    for subset in E88.SUBSETS:
        tnt = E88.flies_subset(subset)
        rf = E88.flies_factors('z', 'r1', tnt)
        F = {c: rf[c].mean() for c in E89.CONDS}
        fl[('r1', subset)] = (E88.log_changes(F, ('L', 'R')), [E88.log_changes(b, ('L', 'R')) for b in E88.flies_boot(rf)])
        fl[('cand', subset)] = (E88.log_changes(E89.flies_cand(fv, CO_CANDIDATE, None, tnt), ('L', 'R')),
                                [E88.log_changes(b, ('L', 'R')) for b in E89.flies_cand_boot(fv, CO_CANDIDATE, None, tnt)])
    rows = []
    for name, us in units.items():
        for unit, wu in us:
            for rdn in ('r1', 'r1b', 'cand_L1', 'cand_L3'):
                Fs = [done[(name, w)]['readings'][rdn] for w in wu]
                Fp = {c: E88.Factors(*(E88._mean_rows(np.stack([f[c][k] for f in Fs])) for k in ('s', 'O', 'R'))) for c in E89.CONDS}
                ch = E88.log_changes(Fp, ('L', 'R'))
                for subset in E88.SUBSETS:
                    fch, fb = fl[('r1' if rdn.startswith('r1') else 'cand', subset)]
                    for f in ('L', 'R'):
                        for grp in ('nearest4', 'ester', 'alcohol'):
                            k = f'chg_{f}_{grp}'
                            lo, hi = np.nanpercentile([b[k] for b in fb], E88.INTERVAL)
                            x = ch[k]
                            differs = bool(not (lo <= x <= hi) and abs(x - fch[k]) >= np.log(E89.B1_FACTOR))
                            rows.append(dict(input=name, unit=unit, reading=rdn, candidate=CO_CANDIDATE if rdn.startswith('cand') else '',
                                             subset=subset, factor=f, group=grp, ratio=float(np.exp(x)), flies_ratio=float(np.exp(fch[k])),
                                             flies_lo=float(np.exp(lo)), flies_hi=float(np.exp(hi)), differs=differs))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'coresponse.csv', index=False)
    return df


def who_responds_runs(inputs: Sequence[str] = ('I0', 'F')) -> pd.DataFrame:
    """who_responds_runs.csv: A1's table of who responds on each input's runs (all wirings it has, and pooled)."""
    runs = runs_table()
    rows = []
    for name in inputs:
        ws = list(runs[runs.input == name].wiring)
        if ws:
            rows += who_responds({w: run_outputs(name, w) for w in ws}, input=name)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'who_responds_runs.csv', index=False)
    return df


def verdict(rd: Optional[pd.DataFrame] = None) -> dict:
    """verdict.json: the identity; part A2's rule; F's gate (C), the equal rule, 'worse' (T's population above the
    reference by more than EQUAL_PATTERN at fixed scales) with the controls' gains, outcome (d) (passes on the 120 pairs
    and not on the 66); the bounds; E87's target without APL and the alcohol pairs' co-response beside; the named pair;
    the runs outside the bands; D's and E's results when they exist; and the outcome."""
    out = out_dir()
    rd = pd.read_csv(out / 'readings.csv') if rd is None else rd
    runs = runs_table()
    v = dict(**H.code_commits(), identity=dict(all_exact=identity_passed()),
             A2=_read_json(out / 'dendrite_reading.json').get('rule'))
    g = gate93(rd, runs, 'F')
    eq = equal93(rd, runs, 'F')
    worse = bool(g.get('complete') and g['gain_all'] <= -EQUAL_PATTERN)
    v['F'] = dict(gate=g, equal=eq, worse=worse)
    if worse:
        # (c): the permuted rows cost as much (each control kind loses at least half as much as F) -> part E's four fits
        # on S are run once to see whether a refit absorbs it
        ctl = g.get('controls', {})
        v['F']['controls_cost_as_much'] = bool(ctl and all(c['n_draws'] > 0 and c['mean_gain'] <= CONTROL_RATIO * g['gain_all']
                                                          for c in ctl.values()))
    bd = bounds93(rd, runs)
    pd.DataFrame(bd).to_csv(out / 'bounds.csv', index=False)
    v['bounds'] = bd
    pd.DataFrame([{k: (json.dumps(x) if isinstance(x, (dict, list)) else x) for k, x in g.items()}]).to_csv(out / 'gate.csv', index=False)
    pd.DataFrame([eq]).to_csv(out / 'equal.csv', index=False)
    if g.get('passes'):
        outcome = 'pass'
    elif g.get('passes_120_not_66'):
        outcome = 'passes on the 120 pairs, not on the 66 megamat pairs (not a pass)'
    elif eq.get('equal'):
        outcome = 'equal'
    elif worse:
        outcome = 'worse'
    elif not g.get('complete'):
        outcome = 'incomplete'
    else:
        outcome = 'neither (not a pass, not an equal)'
    v['outcome'] = outcome
    v['bound_preferred'] = [b for b in bd if b.get('preferred_to_F')]
    p = out / 'apl_free.csv'
    if p.exists():
        af = pd.read_csv(p)
        a = af[(af.reading == 'noise_free') & (af.readout == 'raw')].set_index(['input', 'unit'])
        v['apl_free'] = {f'{n}_{u}': float(a.loc[(n, u), 'pattern']) for n, u in a.index if n in ('I0', 'F')}
    p = out / 'coresponse.csv'
    if p.exists():
        co = pd.read_csv(p)
        c = co[(co.factor == 'L') & (co.group == 'alcohol') & (co.reading == 'r1') & (co.subset == co.subset.iloc[0])]
        v['coresponse_alcohol'] = {f'{r.input}_{r.unit}': dict(ratio=float(r.ratio), flies=float(r.flies_ratio)) for r in c.itertuples()}
    v['pair_15ol_16ol'] = {f'{n}_{u}': _val(rd, n, u, 'effect_15ol_16ol') for n in ('I0', 'F') for u in ('S', 'T', 'all8')
                           if len(rd[(rd.input == n) & (rd.unit == u)])}
    v['pair_15ol_16ol']['data'] = float(E83._cs(targets93()['z8']).get(tuple(sorted(('1-5ol', '1-6ol'))), np.nan))
    v['flagged_runs'] = [dict(tag=r.tag, flag=r.realism_flag) for r in runs.itertuples() if r.flagged]
    for k, f in (('D', 'd_verdict.json'), ('E', 'e_scores.json')):
        if (out / f).exists():
            v[k] = _read_json(out / f)
    v['E_allowed'] = {s: e_allowed(s, v) for s in E_STAGES}
    _write_json(out / 'verdict.json', v)
    print(json.dumps(dict(outcome=outcome, gate={k: g.get(k) for k in ('gain_all', 'gain_halves_mean', 'n_cells_lower', 'i', 'ii',
                                                                       'iii', 'iv', 'v', 'vi', 'passes')},
                          equal=eq.get('equal'), worse=worse), indent=1, default=str), flush=True)
    return v


def score(n_workers: int = 4, coresponse: bool = True, zcell: bool = True) -> dict:
    """Part C's readings and verdict on what part B has run: readings.csv, spread.csv, pairs.csv, apl_free.csv,
    line.csv, coresponse.csv (with `coresponse`), who_responds_runs.csv, gate.csv, equal.csv, bounds.csv, verdict.json,
    the figure, comparison_table.md."""
    rd = readings(('cal', 'raw', 'zcell') if zcell else ('cal', 'raw'))
    apl_free()
    line_reading()
    if coresponse:
        coresponse93(n_workers)
    who_responds_runs()
    v = verdict(rd)
    figure()
    comparison_table()
    return v


# ---------------------------------------------------------------------------
# part D

D_SUBSETS = {'all': {'bouton': ('shared',), 'claw': ('all',)},
             'input46': {'bouton': ('shared',), 'claw': ('input46',)},
             'input36': {'bouton': ('shared',), 'claw': ('input',)}}


def part_d(wirings: Sequence[str] = ALL_WIRINGS, name: str = 'F') -> None:
    """D: each wiring at the input's held point (runs_B.csv's absolute parameters) run again, APL on and off, with APL's
    activity saved and the stage matrices read (E90's `read_stages` over E90's MODEL_WINDOW); each run's counts and
    responses must equal part B's exactly; saved under outputs/.../D/ (E90's layout). Resumable per wiring."""
    import gc
    from experiments import E90_model_stages as E90
    cc = _commits_clean()
    rb = runs_b().set_index('tag')
    d = local_dir('D')
    path = out_dir() / 'identity_D.json'
    ident = _read_json(path)
    orn = input_of(name)
    for w in wirings:
        tags = [f'{name}_{w}_{cc_}' for _, cc_ in E90.CONDS]
        if all((d / f'{t}_claw_win.parquet').exists() for t in tags):
            continue
        tb = run_tag(name, w)
        if tb not in rb.index:
            raise RuntimeError(f'{tb}: part B has not run')
        row = rb.loc[tb]
        saved = load_run(tb)
        t0 = time.time()
        (sim, info), _ = build93(w, orn, name)
        cells = info['cells']
        mp, rv = cells['mp'].cell_contents, cells['rv'].cell_contents
        order = np.asarray(cells['odor_order'].cell_contents, dtype=int)
        mp.kc.save_inh_sims = True
        W = info['W']
        params = {'on': json.loads(row.params_on), 'off': json.loads(row.params_off)}
        for c, cc_ in E90.CONDS:
            r, s = sim(params[c])
            rs, ss = saved[c]
            nd = int((s.to_numpy(dtype=float) != ss.to_numpy(dtype=float)).sum())
            nr = int((r.to_numpy(dtype=float) != rs.to_numpy(dtype=float)).sum())
            ident[f'{name}_{w}_{c}'] = dict(counts_differ=nd, responses_differ=nr, exact=bool(nd == 0 and nr == 0), **cc)
            _write_json(path, ident)
            if nd or nr:
                raise RuntimeError(f'{name} {w} {c}: the stage run differs from part B\'s in {nd} counts, {nr} responses')
            st = E90.read_stages(mp, rv, W, order, [H.odor_abbrev(x) for x in r.columns], E90.MODEL_WINDOW)
            E90._save_stage(d, f'{name}_{w}_{cc_}', st, (r, s))
        print(f'E93 D {w}: exact, {time.time() - t0:.0f}s', flush=True)
        del sim, info, cells, mp, rv
        gc.collect()


def load_stage93(tag: str) -> Optional[dict]:
    """A D stage run (E90's `load_stage` layout) from this host's outputs or the archive."""
    for d in (H.outputs_dir(EXPERIMENT) / 'D', H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'D'):
        if (d / f'{tag}_claw_win.parquet').exists():
            st = {k: pd.read_parquet(d / f'{tag}_{k}.parquet') for k in
                  ('pn_win', 'pn_pre', 'bouton_win', 'bouton_pre', 'claw_win', 'claw_pre')}
            st['apl'] = pd.read_csv(d / f'{tag}_apl.csv', index_col=0)
            st['responses'] = pd.read_parquet(d / f'{tag}_responses.parquet')
            st['spike_counts'] = pd.read_parquet(d / f'{tag}_spike_counts.parquet')
            return st
    return None


def score_d(name: str = 'F') -> dict:
    """Dscore: E90's readings of the stage runs (`read_wiring` on E90's subsets and 'input46', the claws of the 46
    glomeruli with input; `pn_rows`; `soma_rows`), its noise check, the population's group means (stages.csv), the steps
    (d_steps.csv), and E90's verdict by group read with the claws on all claws, on the 46 glomeruli's and on the 36's
    (d_verdict.json), each group's error in the effect being the input's population's (readout, all recordings); E90's
    own verdict beside."""
    from experiments import E90_model_stages as E90
    out = out_dir()
    E90.EXTRA_SUBSETS['input46'] = list(input_of(name).index)
    subsets = tuple(E90.SUBSETS) + ('input46',)
    rows = []
    for i, w in enumerate(ALL_WIRINGS):
        st = {c: load_stage93(f'{name}_{w}_{cc_}') for c, cc_ in E90.CONDS}
        if any(x is None for x in st.values()):
            print(f'E93 Dscore: {w} has no stage runs', flush=True)
            continue
        rows += E90.read_wiring(w, st['on'], st['off'], unit_meta(w), np.random.default_rng(E90.NOISE_SEED + i), subsets=subsets)
        for c in ('on', 'off'):
            rows += E90.pn_rows(w, st[c], c, subsets=subsets)
        K = E83.readout_K({'on': (st['on']['responses'], st['on']['spike_counts'])})
        rows += E90.soma_rows(w, {c: st[c]['spike_counts'] for c in ('on', 'off')}, K)
    mp = pd.DataFrame(rows)
    if mp.empty:
        raise RuntimeError('no stage runs (run --part D)')
    mp.to_parquet(local_dir() / 'd_model_pairs.parquet')
    nc = E90.noise_check(mp, path=out / 'd_noise_check.csv')
    passed = E90.check_passes(nc)
    steps = E90.model_steps_of(mp)
    gw = E90._group_table(mp, E90.KEY + ['wiring'])
    gp = E90._group_table(E90.populations(mp), E90.KEY)
    gp['wiring'] = 'population'
    sm = pd.concat([gp, gw], ignore_index=True)
    sm.to_csv(out / 'stages.csv', index=False)
    E90._group_table(steps, ['step', 'soma', 'subset', 'share', 'reading', 'wiring']).to_csv(out / 'd_steps.csv', index=False)
    ws = sorted(set(mp.wiring), key=list(ALL_WIRINGS).index)
    err = E83.pattern_error(E83.population([_dist(name, w, 'cal') for w in ws]), targets93()['z8'])
    m = E90.group_masks(err.index)
    ee = {g: float(err[m[g]].mean()) for g in E90.GROUPS}
    res = dict(noise_check_passed=passed, effect_error=ee, wirings=ws, readings={})
    vrs = []
    for key, subs in D_SUBSETS.items():
        vr = E90.verdict_rows(sm, steps, passed, subs)
        vrs.append(vr.assign(claws=key))
        res['readings'][key] = {}
        for g in E90.GROUPS:
            er = {el: E90.element_result(vr[(vr.group == g) & (vr.element == el)]) for el in ('bouton', 'claw', 'step_on', 'step_off')}
            res['readings'][key][g] = dict(verdict=E90.group_verdict(er, ee[g]), elements=er)
    pd.concat(vrs, ignore_index=True).to_csv(out / 'd_verdict_rows.csv', index=False)
    e90 = _read_json(H.results_dir(E90.EXPERIMENT) / 'verdict.json')
    res['E90'] = {g: x.get('verdict') for g, x in e90.get('groups', {}).items()}
    res['ester_alcohol_upstream_all_claws'] = res['readings']['all']['ester_alcohol']['verdict'] == 'upstream'
    res['E90_ester_alcohol'] = res['E90'].get('ester_alcohol')
    res.update(**H.code_commits())
    _write_json(out / 'd_verdict.json', res)
    print(json.dumps({k: {g: x['verdict'] for g, x in v.items()} for k, v in res['readings'].items()}, indent=1), flush=True)
    return res


# ---------------------------------------------------------------------------
# part E

def e_member(label: str) -> dict:
    """A wiring's fit of C2's family with its E10_joint_fit.py wiring arguments (E91's `d_member`)."""
    from experiments import E91_summed_drive_inhibition as E91
    return E91.d_member(label)


def e_init_params(name: str = 'F') -> Path:
    """E/init_params_<name>.json in this experiment's results: per wiring, the absolute parameters of part B's held point
    with input `name` (runs_B.csv's params_on of B_<name>_<wiring>: the wiring's own fit with the somatic threshold and the
    tied APL scale re-held for the input), the start of each refit. An exact start at the fit's own parameters
    (`--init-from-exact`) is refused by E10_joint_fit.py for another input, its base being another (2026-11-13)."""
    rb = pd.read_csv(out_dir() / 'runs_B.csv')
    d = {}
    for r in rb[rb.change == name].itertuples():
        d[r.wiring] = {k: float(v) for k, v in json.loads(r.params_on).items()}
    path = out_dir() / 'E' / f'init_params_{name}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(path, d)
    return path


def e_command(label: str, name: str = 'F', maxfev: int = E_MAXFEV) -> List[str]:
    """The E10_joint_fit.py command of one refit: C2's family (`-k 2 --aplpn-max 3 --rate-target ctrl_frac --loss apl
    --apl-readout both --tag _cap3`) with `--orn-input` the input's CSV, into this experiment's folders, from part B's
    held point of the input on the wiring (`--init-params`, `e_init_params`)."""
    m = e_member(label)
    init = e_init_params(name)
    if label not in json.loads(init.read_text()):
        raise RuntimeError(f'no held point of {name} on {label} in runs_B.csv (run part B first)')
    cmd = [str(H.HERE / 'venv' / 'bin' / 'python'), '-u', 'experiments/E10_joint_fit.py', '-k', '2', '--aplpn-max', '3',
           '--rate-target', 'ctrl_frac', '--loss', 'apl', '--apl-readout', 'both', '--tag', '_cap3', '--maxfev', str(int(maxfev)),
           '--orn-input', str(input_path(name).relative_to(H.shared_dir())), '--experiment', EXPERIMENT,
           '--init-params', f'{init}:{label}'] + m['wiring_args']
    return cmd


def e_tag(label: str, name: str = 'F') -> str:
    """The tag E10_joint_fit.py gives a refit (its prefix rules: the wiring, input-<stem>_, rate0.068_, aplboth_, then
    k2_cap3)."""
    a = e_member(label)['wiring_args']
    wt = ''
    if '--wiring' in a:
        wt = f"wiring-{a[a.index('--wiring') + 1]}-seed{a[a.index('--wiring-seed') + 1]}_"
    if '--wiring-tag' in a:
        wt += f"wiring-{a[a.index('--wiring-tag') + 1]}_"
    return f"{wt}input-{name}_rate{H.data_target('ctrl_frac'):.3f}_aplboth_k2_cap3"


def _fits() -> pd.DataFrame:
    p = out_dir() / 'fits.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=['tag'])


REF_S_ALL = 0.0500      # S's unchanged fits, pattern part under the readout against all recordings (E83)


def e_allowed(stage: str, v: Optional[dict] = None) -> Tuple[bool, str]:
    """The entry's rules for the refits: S only on a pass by the gate or an equal (or, outcome (c), when F costs more
    than EQUAL_PATTERN and the permuted controls cost as much, once); T only if S's fitted population is E_S_GAIN or
    more below the reference against all recordings under the readout."""
    v = _read_json(out_dir() / 'verdict.json') if v is None else v
    if not v:
        return False, 'no verdict (run --part score)'
    if stage == 'S':
        f = v.get('F', {})
        if f.get('gate', {}).get('passes'):
            return True, 'the gate passes'
        if f.get('equal', {}).get('equal'):
            return True, 'F is an equal'
        if f.get('worse') and f.get('controls_cost_as_much'):
            return True, 'F costs and the permuted controls cost as much (outcome (c)): the four fits on S once'
        return False, f"the gate does not pass and F is not an equal (outcome: {v.get('outcome')})"
    se = score_e(write=False)
    p = se.get('S', {}).get('pattern_all')
    if p is None:
        return False, 'S not fitted / scored'
    ok = bool(p <= REF_S_ALL - E_S_GAIN)
    return ok, f'S fitted population pattern part {p:.4f} against {REF_S_ALL} - {E_S_GAIN}'


def part_e(stage: str = 'S', force: bool = False, dry_run: bool = False, name: str = 'F') -> List[str]:
    """Runs the stage's refits one after another (`e_command`; a process each, its log under outputs/.../E/), skipping
    those in this experiment's fits.csv; refuses a stage the entry's rules do not allow (`e_allowed`) unless `force`."""
    import subprocess
    ok, why = e_allowed(stage)
    print(f'E93 E {stage}: {"allowed" if ok else "not allowed"} ({why})', flush=True)
    if not ok and not force:
        return []
    _commits_clean()
    write_input(name)
    logs = local_dir('E')
    tags = []
    for w in E_STAGES[stage]:
        tag = e_tag(w, name)
        tags.append(tag)
        if tag in set(_fits().tag):
            print(f'E93 E {tag}: done', flush=True)
            continue
        cmd = e_command(w, name)
        print('E93 E:', ' '.join(cmd), flush=True)
        if dry_run:
            continue
        with open(logs / f'{tag}.log', 'w') as fh:
            r = subprocess.run(cmd, cwd=H.HERE, stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f'E93 E {tag}: E10_joint_fit.py exited {r.returncode} (log {logs / (tag + ".log")})')
        if tag not in set(_fits().tag):
            raise RuntimeError(f'E93 E: the fit wrote no row {tag} (log {logs / (tag + ".log")})')
    return tags


def score_e(write: bool = True, name: str = 'F') -> dict:
    """E's scores: per stage, wiring and population (S, T) E83's `effect_scores` under the readout and on raw counts,
    against the wirings' saved fits (I0), with realism; the family as a population of eight by E91's part D rule
    (E_PASS_GAIN below REF_ALL8, at least E_WIRINGS_LOWER of the eight lower, every fit inside the bands, the raw-count
    loss not above E_RAW_REF by more than E_RAW_SLACK, the mean APL effect within E_MEAN_SLACK of the flies').
    e_scores.csv, e_scores.json, and e_readings.csv: each population, refitted (fit 'refit') and at the saved fits (fit
    'I0'), read like part B's populations (`unit_scores` and `spread_rows` against all recordings, on the 120 and the 66
    megamat pairs, under the readout and on raw counts)."""
    from experiments import E72_wiring_population as E72
    fits = _fits()
    tg = targets93()
    res, rows, got = {}, [], {}
    for stage, ws in E_STAGES.items():
        for w in ws:
            tag = e_tag(w, name)
            if tag in set(fits.tag):
                o = E72.load_fit(EXPERIMENT, tag)
                if o is not None:
                    got[w] = (tag, o)
    for w, (tag, o) in got.items():
        d = {ro: E83.pair_distances(o, ro) for ro in ('cal', 'raw')}
        d0 = {ro: E83.pair_distances(E83.saved_fit(w)['saved'], ro) for ro in ('cal', 'raw')}
        sc, sc0 = E83.effect_scores(d['cal'], tg), E83.effect_scores(d0['cal'], tg)
        fr = fits[fits.tag == tag].iloc[-1]
        rows.append(dict(unit=w, tag=tag, pattern_all=sc['z8_pattern'], pattern_all_I0=sc0['z8_pattern'], loss=sc['z8_loss'],
                         mean_effect=sc['mean_effect'], raw_loss=E83.effect_scores(d['raw'], tg)['z8_loss'],
                         raw_loss_I0=E83.effect_scores(d0['raw'], tg)['z8_loss'], lower=bool(sc['z8_pattern'] < sc0['z8_pattern']),
                         **{c: fr.get(c) for c in ('fixed_thr', 'wAPLKC_scale', 'wKCAPL_scale', 'wAPLPN_scale', 'wPNAPL_scale',
                                                   'kc_kc_inh', 'n_evals')}, fit_loss=fr.get('loss'), **E83.realism(o)))
    by = {r['unit']: r for r in rows}
    readings = []
    for unit, ws in (('S', S_WIRINGS), ('T', T_WIRINGS), ('all8', ALL_WIRINGS)):
        if not all(w in got for w in ws):
            continue
        dp = {ro: E83.population([E83.pair_distances(got[w][1], ro) for w in ws]) for ro in ('cal', 'raw')}
        d0 = {ro: E83.population([E83.pair_distances(E83.saved_fit(w)['saved'], ro) for w in ws]) for ro in ('cal', 'raw')}
        for fit, dd in (('refit', dp), ('I0', d0)):
            for ro in ('cal', 'raw'):
                us = unit_scores(dd[ro], tg)
                scal = {k: us[k] for k in ('mean_effect', 'z8_loss', 'z8_pattern', 'z8_rho', 'z8_pattern66', 'z8_err_ester',
                                           'z8_err_alcohol', 'halfA_pattern', 'halfB_pattern', 'halfC_pattern', 'halfD_pattern',
                                           'mean_on', 'mean_off', 'effect_15ol_16ol') if k in us}
                for pairs in ('all', 'megamat'):
                    for sp in spread_rows(dd[ro], {'z8': tg['z8']}, pairs):
                        readings.append(dict(input=name if fit == 'refit' else 'I0', fit=fit, unit=unit, readout=ro, **{k: v for k, v in sp.items()
                                             if k != 'target'}, **scal))
        sc, sc0 = E83.effect_scores(dp['cal'], tg), E83.effect_scores(d0['cal'], tg)
        r = dict(unit=unit, pattern_all=sc['z8_pattern'], pattern_all_I0=sc0['z8_pattern'], mean_effect=sc['mean_effect'],
                 loss=sc['z8_loss'], raw_loss=E83.effect_scores(dp['raw'], tg)['z8_loss'],
                 raw_loss_I0=E83.effect_scores(d0['raw'], tg)['z8_loss'], n_lower=int(sum(by[w]['lower'] for w in ws)),
                 n_flagged=int(sum(_flagged(by[w].get('realism_flag')) for w in ws)))
        rows.append(r)
        res[unit] = r
    if 'all8' in res:
        e = res['all8']
        res['eight_passes'] = bool(e['pattern_all'] <= REF_ALL8 - E_PASS_GAIN and e['n_lower'] >= E_WIRINGS_LOWER and e['n_flagged'] == 0
                                   and e['raw_loss'] <= E_RAW_REF + E_RAW_SLACK and abs(e['mean_effect'] - FLIES_MEAN) <= E_MEAN_SLACK)
    if write:
        pd.DataFrame(rows).to_csv(out_dir() / 'e_scores.csv', index=False)
        if readings:
            pd.DataFrame(readings).to_csv(out_dir() / 'e_readings.csv', index=False)
        _write_json(out_dir() / 'e_scores.json', res)
    return res


def part_f() -> None:
    """Part F (only on a pass by (C) or (E)): the firing-ceiling family with the response rate with APL on re-held, S at
    the flies' 2 s pulse with both response rates re-held, and S through the trial noise as E72's part E ran it; each
    holds if the population keeps at least half of its gain. Not written: it runs only on a pass."""
    raise NotImplementedError('E93 part F runs only on a pass by part C or part E; it is not written')


# ---------------------------------------------------------------------------
# the archive, the figure and the table

def write_variants() -> pd.DataFrame:
    """variants.csv for the model-output archive and model_pair_distances_long.csv for kcq/apl_effect.py (raw spike
    counts, the harness's primary model RDM): one row per part B run with its absolute parameters, the input, its file and
    sha256, the factors and realism; E's fits (fits.csv rows), whose distances E10_joint_fit.py wrote, are kept."""
    out = out_dir()
    rb = runs_b()
    rows, parts = [], []
    for _, r in rb.iterrows():
        rows.append(dict(tag=r.tag, variant=r.tag, part='B', arm=r.input, base='bouton', k=E83.BREADTH_K, breadth_k=E83.BREADTH_K,
                         kind=(f"C2's family on {r.wiring} with the ORN input {r.input} (results/{EXPERIMENT}/inputs/{r.input}.csv, "
                               f"{int(r.n_input_rows)} glomeruli); breadth recomputed for the input; the somatic threshold and the "
                               "tied APL scale re-held at the fit's APL-off and APL-on counts (E85's held reading)"),
                         config=r.input, wiring=r.wiring, ref_experiment=r.fit_folder, ref_fit=r.fit_tag, orn_input=r.orn_input,
                         fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale, wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale,
                         wPNAPL_scale=r.wPNAPL_scale, wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale, thr_factor=r.thr_factor,
                         apl_factor=r.apl_factor, rate_on=r.rate_on, rate_off=r.rate_off, max_count_on=r.max_count_on,
                         max_count_off=r.max_count_off, realism_flag=r.realism_flag, subdir='B'))
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
        keep = [old[old.tag.isin(set(_fits().tag))]]
    if keep or parts:
        pd.concat(keep + parts, ignore_index=True).to_csv(long_path, index=False)
    return df


def figure() -> Optional[Path]:
    """E93_readings.png: per wiring the pattern part (readout, all recordings) of I0 and F (and the other inputs where
    run), the populations S, T and all eight; and the per-pair APL effect of the eight wirings' populations against the
    KC data."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p = out_dir() / 'readings.csv'
    if not p.exists():
        return None
    rd = pd.read_csv(p)
    r = rd[rd.readout == 'cal']
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    ax = axs[0]
    units = list(ALL_WIRINGS) + ['S', 'T', 'all8']
    names = [n for n in pd.unique(r.input)]
    for k, n in enumerate(names):
        y = [_val(rd, n, u, 'z8_pattern') for u in units]
        ax.plot(np.arange(len(units)) + 0.06 * (k - len(names) / 2), y, 'o', ms=4 if n not in ('I0', 'F') else 6, label=n,
                color={'I0': 'k', 'F': '#d1495b'}.get(n))
    ax.set(xticks=range(len(units)), xticklabels=units, ylabel='pattern part of the APL-effect loss\n(calibrated readout, all recordings)',
           title='each wiring and the populations (lower is better)')
    ax.legend(fontsize=7, ncol=2)
    ax = axs[1]
    pp = out_dir() / 'pairs.csv'
    if pp.exists():
        pr = pd.read_csv(pp)
        de = E83._cs(targets93()['z8'])
        order = np.argsort(de.to_numpy())
        x = np.arange(len(de))
        ax.plot(x, de.to_numpy()[order], color='k', label='KC data')
        for n, col in (('I0', '0.5'), ('F', '#d1495b')):
            q = pr[(pr.input == n) & (pr.unit == 'all8') & (pr.readout == 'cal')]
            if len(q):
                e = E83._cs(pd.Series(q.effect.to_numpy(), index=pd.MultiIndex.from_arrays([q.odor_a, q.odor_b]))).reindex(de.index)
                ax.plot(x, e.to_numpy()[order], lw=0.9, color=col, label=f'eight wirings, {n}')
        ax.set(xlabel='odor pair (ordered by the KC data)', ylabel='APL effect (APL on minus off)\ncosine distance',
               title='per-pair APL effect, calibrated readout')
        ax.legend(fontsize=8)
    path = out_dir() / 'E93_readings.png'
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def _f(x, spec: str = '.4f') -> str:
    try:
        return '–' if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def _md(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False, floatfmt='.4f')
    except ImportError:
        return '```\n' + df.to_string() + '\n```'


def comparison_table() -> str:
    """comparison_table.md from the parts' tables."""
    out = out_dir()
    L = ['# E93: the measured input for the glomeruli that the model\'s input leaves at zero', '',
         f'Generated by experiments/E93_missing_glomeruli_input.py at {H.code_commits().get("queue_commit", "")[:7]}.', '']
    if (out / 'identity.json').exists():
        idj = _read_json(out / 'identity.json')
        L += ['## 0: the identity (I0 with the ten rows as zeros against the saved fits)', '',
              '| wiring | wiring equal | breadth equal | count differences | response differences |', '|---|---|---|---:|---:|']
        for w, d in idj.get('wirings', {}).items():
            L.append(f"| {w} | {d.get('wiring_equal')} | {d.get('breadth_equal')} | {d.get('count_diffs')} | {d.get('response_diffs')} |")
        L += ['', f"All exact: {idj.get('all_exact')}.", '']
    if (out / 'hallem_scaling.csv').exists():
        L += ["## A1: Hallem's data against the lab's input on the shared glomeruli (slope through the origin)", '',
              _md(pd.read_csv(out / 'hallem_scaling.csv').drop(columns=['glomeruli'], errors='ignore')), '']
    for name, title in (('a1_reliability.csv', 'A1: the added rows\' reliability'),
                        ('input_shares.csv', "A1: each odor's added share of claw-weighted positive input (real wiring)"),
                        ('input_distances.csv', 'A1: the input distances by group (66 megamat pairs)')):
        if (out / name).exists():
            L += [f'## {title}', '', _md(pd.read_csv(out / name)), '']
    if (out / 'who_responds.csv').exists():
        w = pd.read_csv(out / 'who_responds.csv')
        w = w[w.unit.isin(['real', 'pooled'])]
        L += ["## A1: who responds in the saved runs (the real wiring and the eight pooled)", '',
              _md(w[['unit', 'kind', 'cls', 'n_kcs', 'rate_off', 'rate_on', 'odors_per_kc_off', 'odors_per_kc_on',
                     'silent_share_off']]), '']
    if (out / 'dendrite_reading.json').exists():
        rdj = _read_json(out / 'dendrite_reading.json')
        L += ['## A2: the PN dendrites at the searched input\'s entries', '',
              '| normalization | entries | sign agreement | ρ | null 5th to 95th | p (null >=) | halves r |', '|---|---:|---:|---:|---:|---:|---:|']
        for n, s in rdj.get('subsets', {}).items():
            for sub in A2_SUBSETS:
                x = s.get(sub, {})
                L.append(f"| {n}, {sub} | {x.get('n')} | {_f(x.get('agree'), '.3f')} | {_f(x.get('rho'), '.2f')} | "
                         f"{_f(x.get('null_p5'), '.3f')} to {_f(x.get('null_p95'), '.3f')} | {_f(x.get('p_ge'), '.3f')} | "
                         f"{_f(s.get('halves_r'), '.2f')} |")
        L += ['', f"Rule: **{rdj.get('rule')}**.", '']
    if (out / 'readings.csv').exists():
        rd = pd.read_csv(out / 'readings.csv')
        cols = [(n, u) for n in pd.unique(rd.input) for u in ('S', 'T', 'all8') if len(rd[(rd.input == n) & (rd.unit == u)])]
        from experiments import E85_orn_fly_population as E85
        ds = E85.data_scores()
        stats = [('pattern part, readout, all recordings (lower is better)', 'cal', 'z8_pattern', '.4f', None),
                 ('pattern part, readout, half A / B / C / D', 'cal', None, None, None),
                 ('APL-effect loss, readout', 'cal', 'z8_loss', '.4f', None),
                 ('mean APL effect, readout', 'cal', 'mean_effect', '+.3f', float(targets93()['z8'].mean())),
                 ('ρ on the per-pair APL effect, readout', 'cal', 'z8_rho', '.2f', None),
                 ('pattern part, readout, the 66 megamat pairs', 'cal', 'z8_pattern66', '.4f', None),
                 ('pattern part, readout, six typical TNT recordings', 'cal', 'z6_pattern', '.4f', None),
                 ('pattern part, readout, noise-free effect', 'cal', 'noise_free_pattern', '.4f', None),
                 ('pattern part, per-KC reading', 'zcell', 'z8_pattern', '.4f', None),
                 ('APL-effect loss, raw counts', 'raw', 'z8_loss', '.4f', None),
                 ('pattern part, raw counts', 'raw', 'z8_pattern', '.4f', None),
                 ('mean distance, APL on / control flies, readout', 'cal', 'mean_on', '.3f', ds['mean_on']),
                 ('mean distance, APL off / APL output blocked (TNT), readout', 'cal', 'mean_off', '.3f', ds['mean_off']),
                 ('ρ on odor-odor distances, APL on, readout', 'cal', 'rho_dist_on', '.2f', None),
                 ('ρ on odor-odor distances, APL off, readout', 'cal', 'rho_dist_off', '.2f', None),
                 ('APL effect of 1-5ol with 1-6ol, readout', 'cal', 'effect_15ol_16ol', '+.3f',
                  float(E83._cs(targets93()['z8']).get(('1-5ol', '1-6ol'), np.nan)))]
        L += ['## B: the populations', '', '| statistic | KC data | ' + ' | '.join(f'{n}<br>{u}' for n, u in cols) + ' |',
              '|---|---:|' + '---:|' * len(cols)]
        for lab_, ro, c, spec, dv in stats:
            if c is None:
                cells = [' / '.join(_f(_val(rd, n, u, f'half{h}_pattern', 'cal'), '.4f') for h in HALVES) for n, u in cols]
            else:
                cells = [_f(_val(rd, n, u, c, ro), spec) for n, u in cols]
            L.append(f'| {lab_} | {_f(dv, spec) if dv is not None else "–"} | ' + ' | '.join(cells) + ' |')
        L.append('')
        sp = out / 'spread.csv'
        if sp.exists():
            s = pd.read_csv(sp)
            s = s[(s.kind == 'population') & (s.target == 'z8') & (s.readout == 'cal')]
            L += ['## The pattern part\'s amplitude and shape parts (readout, all recordings)', '',
                  _md(s[['input', 'unit', 'pairs', 'pattern', 'amplitude', 'shape', 'sd_ratio', 'pearson', 'rho']]), '']
    rt = runs_table()
    if rt is not None and len(rt):
        L += ['## Realism and the factors of every run (rule of 16:10; the earlier rule beside)', '',
              '| run | thr factor | APL factor | rate on | rate off | ratio | response rate ρ on | odors per KC TV on / off | '
              'largest counts on / off | 2h / IaA / pa / 1-5ol / 1-6ol on | silent on / off | flag | flag (earlier rule) |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|']
        for _, r in rt.iterrows():
            fr = ' / '.join(_f(r.get(f'frac_{o}_on'), '.3f') for o in NAMED_ODORS)
            fl, flz = r.get('realism_flag'), r.get('realism_flag_zcrit')
            L.append(f'| {r.tag} | {_f(r.thr_factor, ".4f")} | {_f(r.apl_factor, ".4f")} | {_f(r.rate_on, ".4f")} | {_f(r.rate_off, ".4f")} | '
                     f'{_f(r.ratio, ".2f")} | {_f(r.rate_rho_on, ".2f")} | {_f(r.opk_on_tv, ".3f")} / {_f(r.opk_off_tv, ".3f")} | '
                     f'{_f(r.max_count_on, ".0f")} / {_f(r.max_count_off, ".0f")} | {fr} | {_f(r.silent_share_on, ".3f")} / '
                     f'{_f(r.silent_share_off, ".3f")} | {fl if _flagged(fl) else ""} | {flz if _flagged(flz) else ""} |')
        L.append('')
    for name in ('gate.csv', 'equal.csv', 'bounds.csv', 'apl_free.csv', 'line.csv', 'e_scores.csv'):
        if (out / name).exists():
            try:
                d = pd.read_csv(out / name)
            except pd.errors.EmptyDataError:
                continue
            if name == 'apl_free.csv':
                d = d[d.reading == 'noise_free'][['input', 'unit', 'readout', 'pattern', 'loss', 'mean', 'rho']]
            L += [f'## {name}', '', _md(d), '']
    if (out / 'coresponse.csv').exists():
        co = pd.read_csv(out / 'coresponse.csv')
        co = co[(co.factor == 'L') & (co.group == 'alcohol')]
        L += ["## The alcohol pairs' co-response change under APL", '',
              _md(co[['input', 'unit', 'reading', 'subset', 'ratio', 'flies_ratio', 'flies_lo', 'flies_hi', 'differs']]), '']
    if (out / 'd_verdict.json').exists():
        dv = _read_json(out / 'd_verdict.json')
        L += ['## D: the stages with F (E90\'s verdict by group)', '', '| group | all claws | the 46 glomeruli\'s claws | the 36\'s | E90 |',
              '|---|---|---|---|---|']
        for g in dv['readings']['all']:
            L.append(f"| {g} | {dv['readings']['all'][g]['verdict']} | {dv['readings']['input46'][g]['verdict']} | "
                     f"{dv['readings']['input36'][g]['verdict']} | {dv.get('E90', {}).get(g, '')} |")
        L.append('')
    if (out / 'verdict.json').exists():
        v = _read_json(out / 'verdict.json')
        short = {k: v.get(k) for k in ('outcome', 'A2', 'identity', 'E_allowed', 'pair_15ol_16ol', 'apl_free', 'coresponse_alcohol')}
        short['gate'] = {k: v['F']['gate'].get(k) for k in ('gain_all', 'gain_halves_mean', 'n_cells_lower', 'i', 'ii', 'iii', 'iv', 'v',
                                                            'vi', 'passes', 'controls')} if 'F' in v else None
        short['equal'] = v.get('F', {}).get('equal', {}).get('equal')
        L += ['## Verdict', '', '```', json.dumps(short, indent=1, default=str)[:12000], '```', '']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json's dated records: the commits, olfsysm's build, the model's provenance with the lab's input
    for the simulation parts, the settings and `extra` (part B's per-point provenance is in its 'points')."""
    import olfsysm
    path = out_dir() / f'provenance_{part}.json'
    old = _read_json(path)
    rec = dict(time=time.strftime('%Y-%m-%d %H:%M'), **H.code_commits(), **extra,
               olfsysm_build={k: v for k, v in olfsysm.build_info().items() if k not in ('git_diff', 'cpp_diff')})
    if part in ('0', 'B', 'D', 'E'):
        rec['model_provenance'] = H.model_provenance(orn_deltas=H.apl16_orn_deltas())
    runs = old.get('runs', [])
    runs.append(rec)
    old.update(script='experiments/E93_missing_glomeruli_input.py', part=part, runs=runs,
               settings=dict(ADDED=ADDED, HALLEM_CONC=HALLEM_CONC, CONTROL_SEEDS=CONTROL_SEEDS, CONTROL_DRAW2=CONTROL_DRAW2,
                             GATE=dict(GAIN=GATE_GAIN, CELLS=GATE_CELLS, RAW_SLACK=GATE_RAW_SLACK, MEAN_SLACK=GATE_MEAN_SLACK,
                                       SQ_ALL=GATE_SQ_ALL, SQ_SIX=GATE_SQ_SIX, CONTROL_RATIO=CONTROL_RATIO),
                             EQUAL=dict(PATTERN=EQUAL_PATTERN, RAW=EQUAL_RAW, MEAN=EQUAL_MEAN), BOUND_GAIN=BOUND_GAIN,
                             A2=dict(DEP_FLOOR=DEP_FLOOR, DEP_NORMS=DEP_NORMS, DEP_BOOT=DEP_BOOT, NULL_DRAWS=NULL_DRAWS,
                                     A2_MIN_NORMS=A2_MIN_NORMS, DD_BOUND=DD_BOUND), E_MAXFEV=E_MAXFEV))
    _write_json(path, old)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('inputs', '0', 'A', 'A1', 'A2', 'B', 'score', 'C', 'D', 'Dscore', 'E',
                                                      'Escore', 'F', 'variants', 'table'))
    ap.add_argument('--wirings', default='all8', help="0 / B / D: 'S', 'T', 'all8' or a comma-separated list")
    ap.add_argument('--inputs', default=None, help='B: comma-separated inputs (default: the plan)')
    ap.add_argument('--max-points', type=int, default=None, help='B: stop after this many points')
    ap.add_argument('--stage', default='S', choices=list(E_STAGES), help='E: which fits')
    ap.add_argument('--force', action='store_true', help="E: run fits the entry's rules do not allow")
    ap.add_argument('--dry-run', action='store_true', help='E: print the commands only')
    ap.add_argument('--workers', type=int, default=4, help='score: co-response processes')
    ap.add_argument('--no-co', action='store_true', help='score: skip the co-response readings')
    ap.add_argument('--no-zcell', action='store_true', help="score: skip E73's per-KC reading")
    a = ap.parse_args()
    warnings.simplefilter('ignore', RuntimeWarning)
    t0 = time.time()
    ws = _wirings(a.wirings)
    if a.part in ('0', 'A', 'A1', 'A2', 'B', 'D', 'E'):
        provenance(a.part, wirings=list(ws), inputs=a.inputs, stage=a.stage if a.part == 'E' else None)
    if a.part == 'inputs':
        write_inputs()
    elif a.part == '0':
        identity(ws)
    elif a.part == 'A':
        part_a()
    elif a.part == 'A1':
        part_a1()
    elif a.part == 'A2':
        part_a2()
    elif a.part == 'B':
        part_b(None if a.inputs is None else [x for x in a.inputs.split(',') if x], None if a.wirings == 'all8' else ws, a.max_points)
    elif a.part in ('score', 'C'):
        score(a.workers, coresponse=not a.no_co, zcell=not a.no_zcell)
    elif a.part == 'D':
        part_d(ws)
    elif a.part == 'Dscore':
        score_d()
        verdict()
        comparison_table()
    elif a.part == 'E':
        part_e(a.stage, force=a.force, dry_run=a.dry_run)
    elif a.part == 'Escore':
        score_e()
        verdict()
        comparison_table()
    elif a.part == 'F':
        part_f()
    elif a.part == 'variants':
        write_variants()
    elif a.part == 'table':
        figure()
        comparison_table()
    print(f'E93 {a.part} done in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
