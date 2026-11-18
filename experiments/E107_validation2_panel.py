#!/usr/bin/env python3
"""E107: the four validation2 odors in the panel their ORN input was recorded in: the best family against Remy's
validation2 KC recordings. A diagnostic: it changes no model and moves no loss.

The entry is item E107 of experiment_priorities.md (proposed by viki 2026-11-17 15:40; planned by a Fable subagent,
claude-fable-5-2, spawned from viki 2026-11-17 15:08; moved into the list 2026-11-17 16:34). In Yang's panel the family
holds paa and B-myr 0.50 apart with APL on where the control flies have 0.93; this reads Remy's wild-type KC somata of
the validation2 panel (22 odors, every one recorded in the same ORN session as the model's input for paa, -aPine, euc
and B-myr) to say whether the gap is paa's concentration (a), the family's KC stage (b), the ORN recordings (c) or the
joining of two ORN sessions (d). "The family" is the breadth-2 APL-effect fit (C2, E67) on its eight wirings (E72's
part B: the hemibrain wiring and seven re-drawn), each at its own fit (E83's `saved_fit`, the absolute parameters of
E72's / E67's fits.csv), read as the mean matrix over the wirings (E72's population).

Parts (CLI `--part`):
- 0, the lock (no distance of Remy's opened): the files with their sha256, the map of Remy's odor names onto the ORN
  input's (olfdb's odors table), the odors of each set of pairs, the flies behind each pair by concentration, the least
  number of flies for a pair to enter a per-fly statistic, and every choice the entry leaves open (`lock_definitions`);
  written to experiments/E107_part0_lock.json (LOCK_PATH), which is committed before any reading code runs
  (`require_lock`).
- sets: the odors of Remy's sets of 136 and 153 pairs, which her manifests do not name, identified by the lock's rule
  from her mean table's non-neural columns only (`identify_sets`, sets.json; no KC or ORN distance parsed).
- data: the files read, copied to data/remy_validation2_kc/ in the shared directory with a README (sources, sizes,
  sha256); the 41 MB per-trial file is not copied.
- A (no simulations): Remy's validation2 KC matrix (her mean table and the per-fly matrices), the validation2 ORN input
  alone against it (ρ on odor-odor distances over each set of pairs and over the pairs of paa, B-myr, euc and -aPine),
  the distance of paa from B-myr, euc's and -aPine's pairs by concentration, the ceilings from halves of the flies, her
  megamat somata as the yardstick and the level check against Yang's control flies.
- B (simulations at fixed parameters): the family's eight wirings on the validation2 ORN input and on the megamat ORN
  input (`panel_input`), APL on and acutely off, as fitted (`run_fitted`) and with the response rates re-held
  (`run_held`); each as-fitted run is checked against the saved fit on the odors both panels share (0 count
  differences); resumable; runs under outputs/E107_validation2_panel/B/.
- C (no simulations): the runs read against Remy's matrices (ρ over each set of pairs, E61's rank residuals and shares
  by odor, the pairs of paa and B-myr), the outcome class (`classify`), the realism numbers with no band applied,
  verdict.json, comparison_table.md, the figures.
- variants: variants.csv for the model-output archive.

Outputs: <shared>/results/E107_validation2_panel/; spike counts and responses under outputs/E107_validation2_panel/.

Usage (from the repository root):
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part 0
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part sets
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part data
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part A
    OMP_NUM_THREADS=11 OMP_THREAD_LIMIT=11 nice -n 10 ./venv/bin/python -u experiments/E107_validation2_panel.py --part B
    ./venv/bin/python -u experiments/E107_validation2_panel.py --part C
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
import socket
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from kcq import harness as H

EXPERIMENT = 'E107_validation2_panel'
LOCK_PATH = H.HERE / 'experiments' / 'E107_part0_lock.json'
SETS_PATH = H.HERE / 'experiments' / 'E107_sets.json'
VERDICT_PAIR = ('paa', 'B-myr')
FOUR = ('paa', 'B-myr', 'euc', '-aPine')
CLOSE, APART = 0.70, 0.85          # the entry's bands for the distance of paa from B-myr
ORDER_TOL = 0.05                   # "within 0.05 of its megamat level"
LEVEL_TOL = 0.05                   # her megamat somata against Yang's control flies, difference of the means

# Remy's files (relative to `remy_dir`): what part 0 may read (manifests, stimulus table, quality-control list, panel
# table, column groups), and the distance files the reading parts open once the lock is committed
MANIFESTS = ('fly_counts_stim_rdm_concat.csv', 'stims_per_fly_long.csv', 'Odor-Panels.tex',
             'validation2/kc_soma/stim_to_drop_by_qcthresh.csv', 'figure-05/resampled/column_groups.yaml')
DISTANCE_FILES = ('validation2/kc_soma/xrda_stim_rdm_concat__rmap_nls_cluster__correlation__0.3.nc',
                  'validation2/kc_soma/xrds_stim_rdm_concat.nc', 'figure-05/resampled/tidy_abbrev_dists.tsv',
                  'megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc')
NOT_READ = ('validation2/kc_soma/xrds_trial_rdm_concat.nc',)    # the 41 MB per-trial file: not read, not copied
MEAN_TABLE = 'figure-05/resampled/tidy_abbrev_dists.tsv'

# the lock's choices (`lock_definitions`)
EXCLUDED_STIMULI = ('MethOH @ -2.0', 'pfo @ 0.0', 'geos @ -5.0')
POOLED = {'euc': (-2.0, -1.75), '-aPine': (-2.0, -1.5), 'ger': (-2.0, -3.0)}
CELLS = 'rmap_nls_cluster'
QC_LIST = dict(cell_mask_coord='rmap_nls_cluster', metric='correlation', qc_thresh=0.3)
MIN_FLIES = 4
N_SPLITS, SEED = 1000, 0
MIN_PER_HALF = 1
MEGAMAT_CELL_PREFERENCE = ('rmap_nls_cluster', 'iscell_nls', 'iscell_good_xid0', 'none')
MEGAMAT_METRIC_PREFERENCE = ('cosine', 'correlation')
SET_SIZES = {'136': 17, '153': 18, '231': 22}
NON_NEURAL_GROUPS = ('chem_cols', 'metacyc_cols', 'vcf_cols_remy', 'vcf_cols_anoop', 'vcf_cols_remy_hmds_329',
                     'vcf_cols_remy_hmds_452', 'lotus_cols')
PANELS = ('validation2', 'megamat')
READINGS = ('fitted', 'held')
PARTNERS = ('held', 'other_set', 'megamat66', 'fly_mean', 'correlation', 'orn_cosine')
OLFDB = 'http://hal:8080'


def _host() -> str:
    return socket.gethostname().split('.')[0]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def _canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), default=lambda x: x.item() if hasattr(x, 'item') else str(x))


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(path.read_text()) if path.exists() else None


def out_dir() -> Path:
    d = H.results_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def local_dir(sub: str = '') -> Path:
    d = H.outputs_dir(EXPERIMENT) / sub if sub else H.outputs_dir(EXPERIMENT)
    d.mkdir(parents=True, exist_ok=True)
    return d


def remy_dir() -> Path:
    """Remy's folder on tensor (`Remy/MB ms FA2023/data_for_autocheck/`; the shared directory is
    <tensor>/Tom/claude_kc_remapping), read only."""
    return H.shared_dir().parents[1] / 'Remy' / 'MB ms FA2023' / 'data_for_autocheck'


# ---------------------------------------------------------------------------
# part 0, the lock (manifests and stimulus table only)

def file_table(paths: Sequence[str], root: Optional[Path] = None) -> List[dict]:
    """[{'path', 'size', 'sha256'}] of Remy's files (paths relative to `root`, default `remy_dir`)."""
    root = remy_dir() if root is None else Path(root)
    return [dict(path=p, size=int((root / p).stat().st_size), sha256=_sha256(root / p)) for p in paths]


def olfdb_odors(url: str = OLFDB) -> pd.DataFrame:
    """olfdb's odors table (`<url>/api/odors`): id, name, abbrev, aliases, cid."""
    import urllib.request
    with urllib.request.urlopen(f'{url}/api/odors', timeout=30) as r:
        rows = json.loads(r.read().decode())['odors']
    return pd.DataFrame([dict(id=o['id'], name=o['name'], abbrev=o['abbrev'], aliases=list(o.get('aliases') or []),
                              cid=o.get('cid')) for o in rows])


def _resolve(name: str, olfdb: pd.DataFrame) -> int:
    """The olfdb row (position) whose abbreviation or an alias is `name`: exact first, then ignoring case."""
    for fold in (False, True):
        f = (lambda s: str(s).lower()) if fold else str
        hits = [i for i, r in enumerate(olfdb.itertuples()) if f(name) in {f(r.abbrev), *map(f, r.aliases)}]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError(f'{name!r}: several olfdb odors ({[olfdb.iloc[i].abbrev for i in hits]})')
    raise ValueError(f'{name!r}: no olfdb odor by abbreviation or alias')


def odor_map(names: Sequence[str], olfdb: pd.DataFrame, input_names: Sequence[str]) -> Dict[str, dict]:
    """{Remy's name: {'input': the ORN input's name, 'olfdb_id', 'olfdb_abbrev', 'olfdb_name', 'cid'}}, both names
    resolved to one olfdb odor by abbreviation or alias (case-sensitive first, then case-insensitive); raises when a name
    resolves to no odor or to several, or two input names to one odor. Names with no input odor map to None."""
    by_id: Dict[int, str] = {}
    for n in input_names:
        i = _resolve(n, olfdb)
        oid = int(olfdb.iloc[i].id)
        if oid in by_id and by_id[oid] != n:
            raise ValueError(f'input names {by_id[oid]!r} and {n!r} are one olfdb odor')
        by_id[oid] = n
    out = {}
    for n in names:
        r = olfdb.iloc[_resolve(n, olfdb)]
        out[n] = dict(input=by_id.get(int(r.id)), olfdb_id=int(r.id), olfdb_abbrev=str(r.abbrev), olfdb_name=str(r['name']),
                      cid=None if pd.isna(r.cid) else int(r.cid))
    return out


def split_stim(stim: str) -> Tuple[str, float]:
    """'euc @ -1.75' -> ('euc', -1.75)."""
    name, conc = str(stim).split(' @ ')
    return name, float(conc)


def stimulus_table(panel: str = 'validation2', root: Optional[Path] = None) -> pd.DataFrame:
    """Remy's stimulus table (stims_per_fly_long.csv) for the KC somata of `panel` ('validation2' / 'megamat17'), the
    excluded stimuli (EXCLUDED_STIMULI) out: fly, stim, name (hers), conc, presented, in_rdm (in the manuscript's
    per-fly matrix, after her quality control), qc_dropped, fly_in_manuscript."""
    root = remy_dir() if root is None else Path(root)
    s = pd.read_csv(root / 'stims_per_fly_long.csv')
    s = s[(s.panel == panel) & (s.region == 'kc_soma') & ~s.stim.isin(EXCLUDED_STIMULI)].copy()
    s['name'] = [split_stim(x)[0] for x in s.stim]
    s['conc'] = [split_stim(x)[1] for x in s.stim]
    for c in ('presented', 'in_rdm', 'qc_dropped', 'fly_in_manuscript'):
        s[c] = s[c].astype(bool)
    return s[['fly', 'stim', 'name', 'conc', 'presented', 'in_rdm', 'qc_dropped', 'fly_in_manuscript']].reset_index(drop=True)


def qc_drops(root: Optional[Path] = None, **sel) -> pd.DataFrame:
    """Remy's quality-control drop list (stim_to_drop_by_qcthresh.csv) at `sel` (default QC_LIST): fly
    ('date/fly_num'), stim."""
    root = remy_dir() if root is None else Path(root)
    q = pd.read_csv(root / 'validation2/kc_soma/stim_to_drop_by_qcthresh.csv', sep='\t')
    sel = dict(QC_LIST, **sel)
    q = q[(q.cell_mask_coord == sel['cell_mask_coord']) & (q.metric == sel['metric']) &
          np.isclose(q.qc_thresh.astype(float), float(sel['qc_thresh']))]
    fly = [f'{d}/{int(n)}' for d, n in zip(q.date_imaged, q.fly_num)]
    return pd.DataFrame(dict(fly=fly, stim=q.stim_to_drop.astype(str).to_numpy())).drop_duplicates().reset_index(drop=True)


def qc_identity(stims: pd.DataFrame, drops: pd.DataFrame) -> dict:
    """The stimulus table's qc_dropped flags against the drop list: the (fly, stim) dropped by one and not the other,
    over the stimuli the table holds (presented)."""
    a = set(map(tuple, stims.loc[stims.qc_dropped, ['fly', 'stim']].to_numpy()))
    pres = set(map(tuple, stims.loc[stims.presented, ['fly', 'stim']].to_numpy()))
    b = {x for x in map(tuple, drops[['fly', 'stim']].to_numpy()) if x in pres}
    in_rdm_ok = bool(((stims.presented & ~stims.qc_dropped) == stims.in_rdm).all())
    return dict(table_only=sorted(a - b), list_only=sorted(b - a), equal=bool(a == b), in_rdm_is_presented_not_dropped=in_rdm_ok,
                n_dropped=len(a))


def pooled_name(name: str, conc: float, by_conc=False) -> str:
    """A stimulus's odor name in a matrix: the name, or for the pooled odors read by concentration '<name> @ <conc>'
    (`by_conc` True for all of POOLED, or the odors to read so)."""
    split = set(POOLED) if by_conc is True else set(by_conc or ())
    return f'{name} @ {conc:g}' if name in split and name in POOLED else name


def fly_odors(stims: pd.DataFrame, names: Dict[str, dict], by_conc: bool = False) -> Dict[str, set]:
    """{fly: the odors (the input's names) its matrix holds}."""
    out: Dict[str, set] = {}
    for r in stims[stims.in_rdm].itertuples():
        out.setdefault(r.fly, set()).add(pooled_name(names[r.name]['input'], r.conc, by_conc))
    return out


def flies_by_pair(stims: pd.DataFrame, names: Dict[str, dict]) -> pd.DataFrame:
    """Per pair of odors (the input's names, canonical order) and per combination of the two odors' concentrations:
    the flies whose matrix holds the pair (both stimuli in_rdm) and their number; a row 'all' per pair with every fly."""
    held = {}
    for r in stims[stims.in_rdm].itertuples():
        held.setdefault(r.fly, {})[names[r.name]['input']] = r.conc
    odors = H.odor_order(sorted({names[n]['input'] for n in stims.name}))
    rows = []
    for a, b in H.pair_list(odors):
        flies = sorted(f for f, h in held.items() if a in h and b in h)
        rows.append(dict(odor_a=a, odor_b=b, conc_a='all', conc_b='all', n_flies=len(flies), flies=';'.join(flies)))
        combos = {}
        for f in flies:
            combos.setdefault((held[f][a], held[f][b]), []).append(f)
        if len(combos) > 1:
            for (ca, cb), fs in sorted(combos.items()):
                rows.append(dict(odor_a=a, odor_b=b, conc_a=f'{ca:g}', conc_b=f'{cb:g}', n_flies=len(fs), flies=';'.join(fs)))
    return pd.DataFrame(rows)


def odor_concentrations(stims: pd.DataFrame, names: Dict[str, dict]) -> pd.DataFrame:
    """Per odor (the input's name) and concentration: the flies presented it and the flies whose matrix holds it."""
    rows = []
    for (n, c), g in stims.groupby(['name', 'conc']):
        rows.append(dict(odor=names[n]['input'], remy_name=n, conc=c, n_presented=int(g.presented.sum()),
                         n_in_rdm=int(g.in_rdm.sum()), n_qc_dropped=int(g.qc_dropped.sum())))
    return pd.DataFrame(rows).sort_values(['odor', 'conc']).reset_index(drop=True)


def panel_table_cids(root: Optional[Path] = None) -> Dict[str, Optional[int]]:
    """{abbreviation: CID} from Remy's Odor-Panels.tex (both tables; the CID is the column of integers after the
    abbreviation or after its concentration)."""
    root = remy_dir() if root is None else Path(root)
    out = {}
    for line in (root / 'Odor-Panels.tex').read_text().splitlines():
        if line.count('&') < 4 or line.strip().startswith('name'):
            continue
        cells = [c.strip() for c in line.rstrip('\\ ').split('&')]
        cid = next((int(c) for c in cells[2:4] if re.fullmatch(r'\d+', c)), None)
        out[cells[1]] = cid
    return out


def input_names(panel: str) -> List[str]:
    """The odor names of the model's ORN input of `panel` ('validation2' / 'megamat'), its column order."""
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas
    df = validation2_orn_deltas() if panel == 'validation2' else megamat_orn_deltas()
    return [H.odor_abbrev(c) for c in df.columns]


def input_concentrations(panel: str) -> Dict[str, float]:
    """{odor: log10 concentration} of the model's ORN input of `panel`."""
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas
    df = validation2_orn_deltas() if panel == 'validation2' else megamat_orn_deltas()
    return {H.odor_abbrev(c): float(str(c[-1]).split(' @ ')[1]) for c in df.columns}


def pair_sets_known(odors: Sequence[str]) -> Dict[str, List[str]]:
    """The sets of pairs the manifests name: '231', the 22 odors. The 136 and the 153 are identified after the lock
    (`identify_sets`)."""
    return {'231': H.odor_order(odors)}


def lock_definitions() -> dict:
    """Every choice of the lock as one JSON-able dict (constants of this module only)."""
    return dict(
        experiment=EXPERIMENT,
        files=dict(manifests=list(MANIFESTS), distance_files=list(DISTANCE_FILES), not_read=list(NOT_READ),
                   note='Part 0 reads the manifests only. The distance files are hashed (bytes) and opened by no code '
                        'before this lock is committed.'),
        stimuli=dict(
            excluded=list(EXCLUDED_STIMULI),
            excluded_why="Remy's conventions (fly_counts_stim_rdm_concat.csv: 'ms grid excludes MethOH/pfo/geos@-5'; "
                         "stims_per_fly_tables.tex: 'MethOH @ -2.0 is excluded throughout'): the solvent, methanol and "
                         "geos at 10^-5 are not odors of the manuscript's matrix",
            pooled={k: list(v) for k, v in POOLED.items()},
            pooled_why="her mean table pools each odor's two concentrations under one name (no fly received both)",
            names="Remy's names mapped onto the ORN input's (olfdb's odors table, abbreviation or alias; `odor_map`): "
                  "PAA -> paa, 1-prop -> 1-3ol, every other name unchanged; every table names odors by the input's names"),
        measures=dict(
            verdict="the cosine distance of the manuscript's cells: the mean table's column kc-cosine-rmap_nls_cluster "
                    "(figure-05/resampled/tidy_abbrev_dists.tsv), the measure of the control flies' 0.93 and the "
                    "family's 0.50",
            beside="the correlation distance of the same cells (kc-correlation-rmap_nls_cluster)",
            band_rule=f"close: {CLOSE} or less; apart: {APART} or more; middle between. The band of the verdict is the "
                      "cosine distance's when the correlation distance falls in the same band, else 'middle' (the "
                      "outcome is then read on the pairs' order alone)",
            remy_matrix_for_rho="her mean table's kc-cosine-rmap_nls_cluster (pooled concentrations, as the "
                                "manuscript's); kc-correlation-rmap_nls_cluster beside (partner 'correlation'); the "
                                "mean over flies of her per-fly matrices (pairs held by MIN_FLIES or more) beside "
                                "(partner 'fly_mean')",
            per_fly=f"xrds_stim_rdm_concat.nc (before her quality control) at cell_mask_coord {CELLS}, each metric, "
                    f"with the manuscript's drop list ({QC_LIST}) applied to both metrics, so that both are read on the "
                    "same flies and stimuli; checked against the manuscript's per-fly correlation file "
                    "(xrda_stim_rdm_concat__rmap_nls_cluster__correlation__0.3.nc), which it must equal",
            orn_input="the ORN input alone as the program reads the model's input: 1 - Pearson across the 36 glomeruli "
                      "of the lab's input (`H.apl16_orn_deltas`' rows; validation2's DL4 and VM4 left out, as in the "
                      "input), of validation2_orn_deltas() (22 odors) or megamat_orn_deltas() (17 odors); cosine beside "
                      "(partner 'orn_cosine'); Remy's own ORN columns (orn-correlation-01234, orn-cosine-01234, her "
                      "5 manuscript flies) a cross-check beside, deciding nothing",
            model="cosine distance on spike counts with APL on (raw counts, the program's primary model matrix); the "
                  "family is the mean over its eight wirings of the pair distances (E72's population), each wiring "
                  "given beside"),
        sets=dict(
            known=dict(s231='the 22 odors of the validation2 panel (her mean table, the manuscript grid)'),
            s136_s153=("Remy's sets of 136 and 153 pairs (17 and 18 odors; figure-05 file names). The manifests do not "
                       "name their odors. `identify_sets` identifies them after this lock is committed and before any "
                       "KC or ORN distance is read, from her mean table's odor columns and its non-neural columns only "
                       f"(the column groups {list(NON_NEURAL_GROUPS)} of column_groups.yaml; no kc-, orn- or model "
                       "column parsed): per column, the largest set of odors all of whose pairs are non-missing (odors "
                       "removed one at a time, the one with most missing pairs first, ties by name); per group, the "
                       "intersection over its columns; and the intersection over all of them. A set of n odors (17, "
                       "18) is identified when exactly one distinct n-odor set arises; otherwise it is 'not identified'"),
            carrying="ρ over 'her pairs' is read on the 136 (identified), with the 153 and the 231 beside; if the 136 is "
                     "not identified, on the 231, and the note says so",
            pb="the pairs of paa and B-myr: the pairs of the carrying set involving paa or B-myr (if the carrying set "
               "lacks either odor, those of the 231)",
            four="the pairs of the four odors: those of the set involving paa, B-myr, euc or -aPine; the 6 pairs among "
                 "the four given as values",
            min_flies=MIN_FLIES,
            min_flies_why="a pair enters a per-fly statistic (the mean over flies, the halves, the per-concentration "
                          "readings) when MIN_FLIES or more flies' matrices hold it; 12 of the 231 pairs are held by 2 or "
                          "3 flies (the stimulus table, after her quality control)",
            by_concentration="euc's and -aPine's pairs by concentration from the per-fly matrices, each concentration's "
                             "flies (MIN_FLIES or more), and pooled as in her mean table beside; every table with them "
                             "says that 7 of her 8 flies with euc and 4 of her 9 with -aPine had the odor at 10^-2, "
                             "against an ORN input recorded at 10^-1.75 and 10^-1.5; the verdict does not rest on them. "
                             "ger, pooled too (5 flies at 10^-3, 5 at 10^-2, the input at 10^-2), is read pooled only, "
                             "and its tables say so"),
        halves=dict(n_splits=N_SPLITS, seed=SEED, min_per_half=MIN_PER_HALF,
                    rule="random halves of the flies (floor(n/2) and the rest); each split's two fly means read on "
                         "the pairs both halves hold (MIN_PER_HALF or more flies each) of the set, the pairs held by "
                         "MIN_FLIES or more flies in all; ρ between the halves, its mean over splits; the Spearman-Brown "
                         "value for all flies and its square root, the ceiling (as hong2p.stats.split_half_reliability)"),
        megamat=dict(
            file='megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc (the manifest: 18 flies, the manuscript\'s)',
            reading=f"the mean over her 18 flies of the per-fly matrices (pairs held by MIN_FLIES or more), at the "
                    f"metric (preference {list(MEGAMAT_METRIC_PREFERENCE)}) and cell selection (preference "
                    f"{list(MEGAMAT_CELL_PREFERENCE)}, else the first the file lists) the file holds; no quality-control "
                    "list for megamat in her folder, so none applied",
            level='ρ of a reader against her megamat somata over the 136 pairs of the 17 megamat odors (the megamat '
                  'ORN input\'s panel); over the 66 pairs of the 12 odors of Yang\'s panel beside (partner '
                  '\'megamat66\')',
            like_for_like='the validation2 ρ\'s that are compared with a megamat level use Remy\'s matrix in the '
                          'metric the megamat level is read in'),
        level_check=dict(
            rule=f"her megamat somata against Yang's control flies (H.yang_mean_rdm('TNT3c', <metric>)), the same "
                 f"metric, on the 66 pairs of the 12 odors both have: the two means; flagged when they differ by more "
                 f"than {LEVEL_TOL}",
            consequence="when flagged, the note says so beside the verdict, and the class on the pairs' order alone (as "
                        "for a distance in the middle band) is given beside; the outcome is then held only if both agree"),
        family=dict(
            members="C2 (E67's rate0.068_aplboth_k2_cap3) on the hemibrain wiring and E72 part B's seven re-drawn "
                    "wirings (P1 to P4, M1 to M3), each at its own fit: E83's `saved_fit` (E72's `fitted_params`, the "
                    "absolute fixed_thr, wAPLKC_scale, wKCAPL_scale, wAPLPN_scale, wPNAPL_scale and kc_kc_inh of "
                    "E67's and E72's fits.csv)",
            runner="E83's `build_simulator` (the runner of E93 and E98), breadth held at the lab input's "
                   "(`wiring_parts(label)['breadth']`), so the APL->KC weights are the fit's (E98's choice); APL acutely "
                   "off: the APL->KC and APL->bouton scales 0",
            identity="each as-fitted run reproduces the saved fit's spike counts and responses on the odors its panel "
                     "shares with the fit's panel (validation2: paa, -aPine, euc, B-myr; megamat: the 12) with 0 "
                     "differences, APL on and off; a failure stops part B",
            fitted="as fitted: the fit's parameters, nothing re-tuned (carries the verdict)",
            held="re-held: E77's `hold_rates` on the panel (the somatic threshold bisected with APL off to the "
                 "wiring's own APL-off fraction on its fit panel, target_off / size, then the tied APL scales to the "
                 "control flies' response rate, the registry's ctrl_frac, 6.8 %, with APL on), the counts rounded on "
                 "the panel's KCs x odors; partner 'held'",
            per_wiring="each wiring's own class given beside (descriptive)"),
        outcome_rules=dict(
            quantities="d = the distance of paa from B-myr (verdict measure); pb(r) = ρ of reader r against Remy's "
                       "validation2 matrix over the pairs of paa and B-myr; all(r) = over the carrying set; mm(r) = r's "
                       "megamat level; readers: the ORN input alone ('orn') and the family ('fam')",
            orders="orn_orders: pb(orn) >= mm(orn) - ORDER_TOL; fam_orders_pb: pb(fam) >= mm(fam) - ORDER_TOL; "
                   "fam_orders_all: all(fam) >= mm(fam) - ORDER_TOL",
            classes="band close -> (a). Otherwise: fam_orders_pb and fam_orders_all -> (d); else orn_orders and not "
                    "fam_orders_pb -> (b); else not orn_orders and not fam_orders_pb -> (c); else 'none' (the family "
                    "orders the pairs of paa and B-myr but not the panel). In the middle band the class is read on the "
                    "pairs' order alone (the same rules without the distance)",
            order_tol=ORDER_TOL, close=CLOSE, apart=APART,
            held=f"the class holds when the carrying reading and every partner {list(PARTNERS)} give it (and, if the "
                 "level check is flagged, the class on order alone too); otherwise 'reading-dependent', naming the "
                 "partners that differ. Partners: held (the re-held family); other_set (the 231 when the 136 carries, "
                 "the 136 when the 231 carries and the 136 is identified); megamat66 (the megamat levels on the 66 "
                 "pairs); fly_mean (Remy's mean over flies in place of her mean table, for d and the ρ's); correlation "
                 "(Remy's correlation distance for the ρ's, the megamat level read likewise); orn_cosine (the ORN input "
                 "as cosine distance)"),
        realism="the family's response rates with APL on and acutely off (two or more spikes), their ratio, the "
                "per-odor response fractions, the distribution of odors per KC (OPK_EDGES' classes over the panel's "
                "odors) and the largest spike counts, on each panel and reading, per wiring and for the population; "
                "no band applied; Remy's files hold distances only, so her response rates are not beside them",
    )


def definitions_sha(defs: dict) -> str:
    return hashlib.sha256(_canonical(defs).encode()).hexdigest()


def lock_contents(root: Optional[Path] = None, olfdb: Optional[pd.DataFrame] = None) -> dict:
    """What part 0 reads and derives (no distance file opened): the files' sha256, the odor maps, the stimulus tables'
    summaries, the known set, the flies by pair and by concentration, the quality-control identity."""
    root = remy_dir() if root is None else Path(root)
    olfdb = olfdb_odors() if olfdb is None else olfdb
    v = stimulus_table('validation2', root)
    m = stimulus_table('megamat17', root)
    cids = panel_table_cids(root)
    maps = {}
    for panel, st, key in (('validation2', v, 'validation2'), ('megamat', m, 'megamat17')):
        mp = odor_map(sorted(set(st.name)), olfdb, input_names(panel))
        missing = [n for n, x in mp.items() if x['input'] is None]
        if missing:
            raise RuntimeError(f'{panel}: odors without input {missing}')
        for n, x in mp.items():
            x['cid_remy_table'] = cids.get(n)
            x['cid_agrees'] = (x['cid_remy_table'] == x['cid']) if x['cid_remy_table'] is not None else None
        maps[key] = mp
    conc_in = {p: input_concentrations(p) for p in PANELS}
    odors_v = H.odor_order(sorted({maps['validation2'][n]['input'] for n in v.name}))
    odors_m = H.odor_order(sorted({maps['megamat17'][n]['input'] for n in m.name}))
    fp_v = flies_by_pair(v, maps['validation2'])
    fp_m = flies_by_pair(m, maps['megamat17'])
    conc_v = odor_concentrations(v, maps['validation2'])
    conc_v['input_conc'] = [conc_in['validation2'][o] for o in conc_v.odor]
    conc_v['within_one_log'] = (conc_v.conc - conc_v.input_conc).abs() <= 1.0
    conc_m = odor_concentrations(m, maps['megamat17'])
    conc_m['input_conc'] = [conc_in['megamat'][o] for o in conc_m.odor]
    qc = qc_identity(v, qc_drops(root))
    allp = fp_v[fp_v.conc_a == 'all']
    allm = fp_m[fp_m.conc_a == 'all']
    return dict(
        files=dict(manifests=file_table(MANIFESTS, root), distance_files=file_table(DISTANCE_FILES, root),
                   root=str(root)),
        odor_map=maps,
        olfdb=dict(url=f'{OLFDB}/api/odors', n_odors=int(len(olfdb))),
        validation2=dict(
            n_flies=int(v.fly.nunique()), flies=sorted(v.fly.unique()), odors=odors_v,
            sets=pair_sets_known(odors_v), concentrations=conc_v.to_dict(orient='records'),
            flies_by_pair=fp_v.to_dict(orient='records'),
            pairs_by_n_flies={str(k): int(c) for k, c in allp.n_flies.value_counts().sort_index().items()},
            n_pairs_min_flies=int((allp.n_flies >= MIN_FLIES).sum()), qc_identity=qc,
            verdict_pair_flies=fp_v[(fp_v.odor_a == 'paa') & (fp_v.odor_b == 'B-myr')].to_dict(orient='records')),
        megamat17=dict(
            n_flies=int(m.fly.nunique()), flies=sorted(m.fly.unique()), odors=odors_m,
            concentrations=conc_m.to_dict(orient='records'), flies_by_pair=fp_m.to_dict(orient='records'),
            pairs_by_n_flies={str(k): int(c) for k, c in allm.n_flies.value_counts().sort_index().items()},
            n_pairs_min_flies=int((allm.n_flies >= MIN_FLIES).sum()),
            yang_shared=[o for o in odors_m if o in H.MEGAMAT12]),
        input_concentrations=conc_in,
    )


def write_lock(force: bool = False, root: Optional[Path] = None) -> dict:
    """Part 0: writes LOCK_PATH (the definitions with their sha256, `lock_contents`, host, time, commits); refuses to
    overwrite an existing lock unless `force`."""
    if LOCK_PATH.exists() and not force:
        raise RuntimeError(f'{LOCK_PATH.name} exists; not overwritten (pass force)')
    defs = lock_definitions()
    lock = dict(experiment=EXPERIMENT, part='0', definitions_sha256=definitions_sha(defs), definitions=defs,
                **lock_contents(root), host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
                commits={k: v for k, v in H.code_commits().items() if not isinstance(v, (dict, list))},
                note="Written by part 0 from Remy's manifests, stimulus table, quality-control list and panel table "
                     "alone, before any distance file was opened (the distance files' bytes are hashed only). The "
                     "reading code refuses to run until this file is committed and unchanged (`require_lock`).")
    lock['sha256'] = definitions_sha({k: v for k, v in lock.items() if k != 'sha256'})
    _write_json(LOCK_PATH, lock)
    return lock


def _committed(path: Path) -> bool:
    """Whether `path` is tracked by its git repository and has no uncommitted change."""
    repo = Path(path).resolve().parent
    try:
        r = subprocess.run(['git', '-C', str(repo), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
        if r.returncode != 0:
            return False
        r = subprocess.run(['git', '-C', str(repo), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
        return r.returncode == 0
    except OSError:
        return False


_LOCK_CACHE: Dict[tuple, dict] = {}


def require_lock(check_git: bool = True, check_files: bool = True) -> dict:
    """The committed lock, after checking that its sha256 is its contents', its definitions are the module's
    (`lock_definitions`), it is committed and unchanged in git (`check_git`) and Remy's files still have the locked
    sha256 (`check_files`); raises PermissionError otherwise (the reading code calls this before opening any
    distance)."""
    key = (check_git, check_files)
    if key in _LOCK_CACHE:
        return _LOCK_CACHE[key]
    lock = _read_json(LOCK_PATH)
    if lock is None:
        raise PermissionError(f'{LOCK_PATH.name}: no lock; run part 0 and commit it first')
    if lock.get('sha256') != definitions_sha({k: v for k, v in lock.items() if k != 'sha256'}):
        raise PermissionError(f'{LOCK_PATH.name}: its sha256 is not its contents\'')
    if _canonical(lock['definitions']) != _canonical(lock_definitions()):
        raise PermissionError(f"{LOCK_PATH.name}: its definitions are not the module's (lock_definitions)")
    if check_git and not _committed(LOCK_PATH):
        raise PermissionError(f'{LOCK_PATH.name} is not committed (or has uncommitted changes)')
    if check_files:
        root = Path(lock['files']['root'])
        for kind in ('manifests', 'distance_files'):
            for f in lock['files'][kind]:
                if _sha256(root / f['path']) != f['sha256']:
                    raise PermissionError(f"{f['path']}: sha256 differs from the lock's")
    _LOCK_CACHE[key] = lock
    return lock


# ---------------------------------------------------------------------------
# the sets of 136 and 153 pairs (the lock's rule, after the lock)

def column_groups(root: Optional[Path] = None) -> Dict[str, List[str]]:
    """Remy's column_groups.yaml: {group: [columns]}."""
    import yaml
    root = remy_dir() if root is None else Path(root)
    return yaml.safe_load((root / 'figure-05/resampled/column_groups.yaml').read_text())


def complete_odors(table: pd.DataFrame, col: str, odors: Sequence[str]) -> List[str]:
    """The largest set of `odors` all of whose pairs have a non-missing `col` in `table` (columns 'a', 'b', both orders
    of a pair allowed), by the lock's rule: odors removed one at a time, the one with most missing pairs among those
    kept first, ties by name."""
    ok = {}
    for a, b, v in zip(table['a'], table['b'], table[col]):
        if a != b:
            ok[frozenset((a, b))] = ok.get(frozenset((a, b)), False) or bool(pd.notna(v))
    keep = sorted(set(odors))
    while True:
        miss = {o: sum(not ok.get(frozenset((o, x)), False) for x in keep if x != o) for o in keep}
        worst = max(miss.values()) if miss else 0
        if worst == 0:
            return keep
        keep.remove(sorted(o for o, n in miss.items() if n == worst)[0])


def identify_sets(root: Optional[Path] = None, write: bool = True) -> dict:
    """The lock's rule for the 136 and the 153 (`lock_definitions()['sets']['s136_s153']`): reads the mean table's odor
    columns and its non-neural columns only (`usecols`; no kc-, orn- or model column parsed), per column the largest
    complete set of odors (`complete_odors`), per group the intersection over its columns, and over all of them; a size
    is identified when exactly one distinct set of that size arises. Writes SETS_PATH (`write`)."""
    lock = require_lock()
    root = remy_dir() if root is None else Path(root)
    names = lock['odor_map']['validation2']
    groups = column_groups(root)
    header = pd.read_csv(root / MEAN_TABLE, sep='\t', nrows=0).columns.tolist()
    odor_cols = header[:2]
    cols = {g: [c for c in groups.get(g, []) if c in header] for g in NON_NEURAL_GROUPS}
    use = odor_cols + sorted({c for cs in cols.values() for c in cs})
    if any(c.startswith(('kc-', 'orn-')) or c in groups.get('model_cols', []) for c in use):
        raise RuntimeError('a KC, ORN or model column among the columns to read')
    t = pd.read_csv(root / MEAN_TABLE, sep='\t', usecols=use)
    unknown = sorted((set(t[odor_cols[0]]) | set(t[odor_cols[1]])) - set(names))
    if unknown:
        raise RuntimeError(f'odor names not in the lock\'s map: {unknown}')
    t['a'] = [names[x]['input'] for x in t[odor_cols[0]]]
    t['b'] = [names[x]['input'] for x in t[odor_cols[1]]]
    odors = lock['validation2']['odors']
    per_col = {c: complete_odors(t, c, odors) for cs in cols.values() for c in cs}
    cands = {}
    for g, cs in cols.items():
        if cs:
            cands[g] = sorted(set.intersection(*[set(per_col[c]) for c in cs]))
    cands['all'] = sorted(set.intersection(*[set(v) for v in per_col.values()]))
    out = dict(rule=lock['definitions']['sets']['s136_s153'], odor_columns=odor_cols, columns=cols,
               per_column={c: len(v) for c, v in per_col.items()}, candidates=cands, sets={},
               outside_rule={c: H.odor_order(v) for c, v in per_col.items() if len(v) in (17, 18)},
               lock_sha256=lock['sha256'],
               mean_table_sha256=_sha256(root / MEAN_TABLE), time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    for key, n in SET_SIZES.items():
        if key == '231':
            out['sets'][key] = dict(identified=True, odors=H.odor_order(odors), source='the lock (all 22 odors)')
            continue
        distinct = sorted({tuple(v) for v in cands.values() if len(v) == n})
        if len(distinct) == 1:
            src = [g for g, v in cands.items() if tuple(v) == distinct[0]]
            out['sets'][key] = dict(identified=True, odors=H.odor_order(distinct[0]), source=src,
                                    left_out=H.odor_order(sorted(set(odors) - set(distinct[0]))))
        else:
            out['sets'][key] = dict(identified=False, n_candidates=len(distinct),
                                    candidates=[H.odor_order(d) for d in distinct])
    if write:
        _write_json(SETS_PATH, out)
    return out


def require_sets() -> dict:
    """The committed sets file (SETS_PATH) after `require_lock`; raises PermissionError when it is missing, uncommitted
    or written under another lock."""
    lock = require_lock()
    s = _read_json(SETS_PATH)
    if s is None or s.get('lock_sha256') != lock['sha256']:
        raise PermissionError(f'{SETS_PATH.name}: missing or written under another lock; run --part sets')
    if not _committed(SETS_PATH):
        raise PermissionError(f'{SETS_PATH.name} is not committed')
    return s


# ---------------------------------------------------------------------------
# reading Remy's matrices (only after the lock)

def pkey(a: str, b: str) -> Tuple[str, str]:
    """A pair's canonical orientation (`hong2p.rdm.pair_list` with the APL16 order first, then by name)."""
    return tuple(H.odor_order([a, b]))


def pair_series(values: Dict[Tuple[str, str], float], name: Optional[str] = None) -> pd.Series:
    """A pair Series from {(a, b): value} (either orientation), keys canonical, in `H.pair_list` order of its odors
    (pairs absent are left out, not NaN)."""
    d = {pkey(a, b): float(v) for (a, b), v in values.items() if a != b}
    odors = H.odor_order(sorted({o for p in d for o in p}))
    keys = [p for p in H.pair_list(odors) if p in d]
    return pd.Series([d[p] for p in keys], index=pd.MultiIndex.from_tuples(keys, names=['odor_a', 'odor_b']), name=name)


def _open(rel: str):
    import xarray as xr
    require_lock()
    return xr.open_dataset(remy_dir() / rel)


def _drops(qc: str) -> set:
    """The (fly, stim) dropped: qc 'lock' (the manuscript's list, QC_LIST), 'cosine' (her list for the cosine
    metric, the same cells and threshold; beside only) or 'none'."""
    if qc == 'none':
        return set()
    sel = {} if qc == 'lock' else dict(metric='cosine')
    return set(map(tuple, qc_drops(**sel)[['fly', 'stim']].to_numpy()))


def _fly_rows(X: np.ndarray, flies: Sequence[str], stims: Sequence[str], names: Dict[str, dict], drops: set,
              by_conc=False, exclude: Sequence[str] = EXCLUDED_STIMULI) -> pd.DataFrame:
    """Flies x pairs from per-fly stimulus x stimulus matrices X (flies x stims x stims): the excluded stimuli and the
    dropped (fly, stim) out, a stimulus a fly lacks (NaN diagonal) out, the odors named by the input's names (the
    pooled ones by concentration with `by_conc`); raises if a fly holds one odor at two concentrations."""
    rows = {}
    for k, f in enumerate(flies):
        keep = [i for i, s in enumerate(stims) if s not in exclude and (f, s) not in drops and not np.isnan(X[k, i, i])]
        nm = {}
        for i in keep:
            n, c = split_stim(stims[i])
            o = pooled_name(names[n]['input'], c, by_conc)
            if o in nm.values():
                raise RuntimeError(f'{f}: {o} twice')
            nm[i] = o
        vals = {}
        for i, j in itertools.combinations(keep, 2):
            v = X[k, i, j]
            if not np.isnan(v):
                vals[(nm[i], nm[j])] = v
        rows[f] = pair_series(vals)
    out = pd.DataFrame(rows).T
    out.index.name = 'fly'
    odors = H.odor_order(sorted({o for p in out.columns for o in p}))
    return out[[p for p in H.pair_list(odors) if p in out.columns]]


def per_fly(measure: str = 'cosine', cells: str = CELLS, by_conc=False, qc: str = 'lock') -> pd.DataFrame:
    """Remy's per-fly validation2 matrices (xrds_stim_rdm_concat.nc, Fc_zscore, at `cells` and `measure`, with the
    lock's drop list, `qc`) as flies x pairs (NaN where a fly lacks a pair); with `by_conc` (True: euc, -aPine and ger;
    or the odors to split), those odors keep their concentration in the odor name."""
    lock = require_lock()
    d = _open('validation2/kc_soma/xrds_stim_rdm_concat.nc')
    try:
        r = [i for i, (c, m) in enumerate(zip(d.cell_mask_coord.values, d.metric.values))
             if str(c) == cells and str(m) == measure]
        if len(r) != 1:
            raise RuntimeError(f'{cells} / {measure}: {len(r)} matrices')
        X = np.asarray(d['Fc_zscore'].values[:, r[0]], dtype=float)
        stims = [str(s) for s in d.stim_row.values]
        if stims != [str(s) for s in d.stim_col.values]:
            raise RuntimeError('rows and columns differ')
        flies = [f'{a}/{int(b)}' for a, b in zip(d.date_imaged.values, d.fly_num.values)]
    finally:
        d.close()
    return _fly_rows(X, flies, stims, lock['odor_map']['validation2'], _drops(qc), by_conc)


def manuscript_per_fly() -> pd.DataFrame:
    """The manuscript's per-fly correlation matrices (xrda_stim_rdm_concat__rmap_nls_cluster__correlation__0.3.nc) as
    flies x pairs."""
    lock = require_lock()
    d = _open('validation2/kc_soma/xrda_stim_rdm_concat__rmap_nls_cluster__correlation__0.3.nc')
    try:
        X = np.asarray(d['Fc_zscore'].values, dtype=float)
        stims = [str(s) for s in d.stim_row.values]
        flies = [f'{a}/{int(b)}' for a, b in zip(d.date_imaged.values, d.fly_num.values)]
    finally:
        d.close()
    return _fly_rows(X, flies, stims, lock['odor_map']['validation2'], set())


def mean_table(col: str = f'kc-cosine-{CELLS}') -> pd.Series:
    """One column of Remy's mean table (figure-05/resampled/tidy_abbrev_dists.tsv; every ordered pair, both orders
    equal) as a canonical pair Series over the 22 odors (the input's names; concentrations pooled as in her table)."""
    lock = require_lock()
    names = lock['odor_map']['validation2']
    t = pd.read_csv(remy_dir() / MEAN_TABLE, sep='\t', usecols=['Unnamed: 0', 'Unnamed: 1', col])
    t = t[t['Unnamed: 0'] != t['Unnamed: 1']]
    vals: Dict[Tuple[str, str], List[float]] = {}
    for a, b, v in zip(t['Unnamed: 0'], t['Unnamed: 1'], t[col]):
        vals.setdefault(pkey(names[a]['input'], names[b]['input']), []).append(float(v))
    bad = [p for p, v in vals.items() if len(v) != 2 or not (np.isclose(v[0], v[1]) or np.isnan(v).all())]
    if bad:
        raise RuntimeError(f'{col}: pairs without two equal orders: {bad[:5]}')
    s = pair_series({p: v[0] for p, v in vals.items()}, name=col)
    return H.canonical_pairs(s.dropna()) if not s.isna().any() else s


def fly_mean(rows: pd.DataFrame, min_flies: int = MIN_FLIES) -> pd.Series:
    """The mean over flies (rows) of each pair held by at least `min_flies` flies."""
    n = rows.notna().sum()
    return rows.mean()[n >= min_flies]


def split_half_ceiling(rows: pd.DataFrame, pairs: Optional[Sequence[Tuple[str, str]]] = None, n_splits: int = N_SPLITS,
                       seed: int = SEED, min_per_half: int = MIN_PER_HALF, min_flies: int = MIN_FLIES,
                       method: str = 'spearman') -> Dict[str, float]:
    """The ceiling of ρ on odor-odor distances from random halves of the flies (rows): the pairs of `pairs` (all if
    None) held by `min_flies` or more flies; each split's two halves' means read on the pairs both halves hold
    (`min_per_half` or more flies each); the mean split-half correlation, its 5th and 95th percentiles, the mean number
    of pairs, the Spearman-Brown value for all flies and its square root (the ceiling), as
    `hong2p.stats.split_half_reliability`."""
    from scipy import stats
    cols = list(rows.columns if pairs is None else [p for p in rows.columns if p in set(map(tuple, pairs))])
    R = rows[cols]
    R = R.loc[:, R.notna().sum() >= min_flies]
    X = R.to_numpy(dtype=float)
    ok = ~np.isnan(X)
    rng = np.random.default_rng(seed)
    n = len(R)
    corr = stats.spearmanr if method == 'spearman' else stats.pearsonr
    vals, npairs = [], []
    for _ in range(n_splits):
        idx = rng.permutation(n)
        h = [idx[:n // 2], idx[n // 2:]]
        cnt = [ok[i].sum(axis=0) for i in h]
        keep = (cnt[0] >= min_per_half) & (cnt[1] >= min_per_half)
        if keep.sum() < 3:
            vals.append(np.nan)
            continue
        with np.errstate(invalid='ignore'):
            m = [np.nanmean(X[i][:, keep], axis=0) for i in h]
        vals.append(float(corr(m[0], m[1])[0]))
        npairs.append(int(keep.sum()))
    vals = np.asarray(vals, dtype=float)
    r = float(np.nanmean(vals))
    full = 2 * r / (1 + r)
    return dict(n_flies=int(n), n_pairs=int(R.shape[1]), mean_pairs_per_split=float(np.mean(npairs)) if npairs else np.nan,
                split_half=r, split_half_lo=float(np.nanpercentile(vals, 5)), split_half_hi=float(np.nanpercentile(vals, 95)),
                spearman_brown_full=float(full), ceiling=float(np.sqrt(full)) if full > 0 else np.nan)


def megamat_somata(metric: str = 'cosine') -> Tuple[pd.DataFrame, dict]:
    """Remy's megamat KC somata per fly (megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc, Fc_zscore) as flies x pairs, at
    `metric` if the file holds it (else the next of MEGAMAT_METRIC_PREFERENCE) and the first cell selection of
    MEGAMAT_CELL_PREFERENCE it holds; and what was read."""
    lock = require_lock()
    d = _open('megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc')
    try:
        pairs = [(str(c), str(m)) for c, m in zip(d.cell_mask_coord.values, d.metric.values)]
        cells = next((c for c in MEGAMAT_CELL_PREFERENCE if any(p[0] == c for p in pairs)), pairs[0][0])
        metrics = [metric] + [m for m in MEGAMAT_METRIC_PREFERENCE if m != metric]
        met = next(m for m in metrics if (cells, m) in pairs)
        r = pairs.index((cells, met))
        X = np.asarray(d['Fc_zscore'].values[:, r], dtype=float)
        stims = [str(s) for s in d.stim_row.values]
        flies = [f'{a}/{int(b)}' for a, b in zip(d.date_imaged.values, d.fly_num.values)]
        info = dict(file='megamat17/kc_soma/xrds_stim_rdm_concat_mm.nc', cells=cells, metric=met, asked=metric,
                    held=pairs, n_flies=len(flies), panel_names=[str(x) for x in d.panel_name.values])
    finally:
        d.close()
    return _fly_rows(X, flies, stims, lock['odor_map']['megamat17'], set(), exclude=()), info


def remy_megamat(metric: str = 'cosine') -> Tuple[pd.Series, dict]:
    """Her megamat somata's mean matrix (`fly_mean` of `megamat_somata`) and what was read."""
    rows, info = megamat_somata(metric)
    return fly_mean(rows), info


def yang_control(metric: str = 'cosine') -> pd.Series:
    """Yang's control flies' mean matrix (`H.yang_mean_rdm('TNT3c', ...)`, 'cosine' or 'correlation' -> 'pearson')."""
    return H.yang_mean_rdm('TNT3c', 'cosine' if metric == 'cosine' else 'pearson')


def level_check(metric: str = 'cosine') -> dict:
    """Her megamat somata against Yang's control flies on the 66 pairs of the 12 odors both have, the same metric: the
    two means, their difference (Remy minus Yang), flagged above LEVEL_TOL; ρ on odor-odor distances between them."""
    rm, info = remy_megamat(metric)
    y = yang_control(info['metric'])
    keep = [p for p in rm.index if p[0] in H.MEGAMAT12 and p[1] in H.MEGAMAT12]
    a, b = H.align_pairs(rm.loc[keep], H.canonical_pairs(y))
    diff = float(a.mean() - b.mean())
    return dict(metric=info['metric'], cells=info['cells'], n_pairs=int(len(a)), remy_mean=float(a.mean()),
                yang_mean=float(b.mean()), difference=diff, flagged=bool(abs(diff) > LEVEL_TOL), tol=LEVEL_TOL,
                rho=H.spearman(a, b))


# ---------------------------------------------------------------------------
# the ORN input

def panel_input(panel: str) -> pd.DataFrame:
    """The model's ORN input for `panel` ('validation2': `validation2_orn_deltas`, 22 odors; 'megamat':
    `megamat_orn_deltas`, 17 odors) on the 36 glomeruli of the lab's input (`H.apl16_orn_deltas`, its row order), its
    (panel, odor) columns as loaded. The columns the lab's input takes from it are identical to the lab's."""
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas
    if panel not in PANELS:
        raise ValueError(f'{panel=}')
    df = validation2_orn_deltas() if panel == 'validation2' else megamat_orn_deltas()
    lab = H.apl16_orn_deltas()
    out = df.loc[list(lab.index)]
    if out.isna().any().any():
        raise RuntimeError(f'{panel}: NaN in the input')
    shared = [c for c in out.columns if c in lab.columns]
    if not shared or not out[shared].equals(lab[shared]):
        raise RuntimeError(f"{panel}: the columns shared with the lab's input differ from it")
    return out


def shared_odors(panel: str) -> List[str]:
    """The odors of `panel`'s input that the fit's panel (APL16) has, as the lab's input has them."""
    lab = H.apl16_orn_deltas()
    return [H.odor_abbrev(c) for c in panel_input(panel).columns if c in lab.columns]


def input_pairs(panel: str, measure: str = 'pearson') -> pd.Series:
    """The ORN input alone as odor-odor distances: 1 - Pearson across the 36 glomeruli (the program's reading;
    'cosine' beside), canonical pairs."""
    return H.canonical_pairs(H.pair_distances(panel_input(panel), measure=measure))


# ---------------------------------------------------------------------------
# part B, the runs

def _frame_sha(df: pd.DataFrame) -> str:
    from experiments.E77_input_screen import frame_sha256
    return frame_sha256(df)


def _abbrev(o: Tuple[pd.DataFrame, pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    return tuple(H.abbrev_columns(x) for x in o)


def run_identity(label: str, o: Dict[str, tuple], panel: str) -> dict:
    """The as-fitted run against the saved fit on the odors `panel` shares with the fit's panel: the number of spike
    counts and responses that differ, APL on and off (0 expected)."""
    from experiments import E83_cross_wiring_fit as E83
    saved = E83.saved_fit(label)['saved']
    odors = shared_odors(panel)
    res = dict(wiring=label, panel=panel, odors=odors)
    for c in ('on', 'off'):
        rs, ss = _abbrev(saved[c])
        r, s = o[c]
        same_index = bool(s.index.equals(ss.index) and r.index.equals(rs.index))
        res[c] = dict(index_equal=same_index,
                      counts_differ=int((s[odors].to_numpy(dtype=float) != ss[odors].to_numpy(dtype=float)).sum()) if same_index else -1,
                      responses_differ=int((r[odors].to_numpy(dtype=float) != rs[odors].to_numpy(dtype=float)).sum()) if same_index else -1)
    res['reproduces'] = bool(all(res[c]['counts_differ'] == 0 and res[c]['responses_differ'] == 0 for c in ('on', 'off')))
    return res


def run_fitted(label: str, panel: str) -> Dict[str, object]:
    """One wiring at its fit on `panel`'s input (E83's `build_simulator`, breadth held at the lab input's, so the
    APL->KC weights are the fit's; nothing re-tuned), APL on and acutely off (the APL->KC and APL->bouton scales 0):
    {'on', 'off': (responses, spike counts), 'params', 'params_off', 'n_sims'}."""
    from experiments import E83_cross_wiring_fit as E83
    P = E83.wiring_parts(label)
    p = dict(P['p'])
    p_off = dict(p, wAPLKC_scale=0.0, wAPLPN_scale=0.0)
    sim, info = E83.build_simulator(label, orn=panel_input(panel), breadth=P['breadth'], params=p)
    on, off = _abbrev(sim(p)), _abbrev(sim(p_off))
    del sim, info
    return dict(on=on, off=off, params=p, params_off=p_off, n_sims=3)


def held_targets(label: str, panel: str) -> Tuple[int, int, float, float]:
    """(target_off, target_on, fraction_off, fraction_on) for the re-held reading on `panel`: the wiring's own APL-off
    fraction on its fit panel and the control flies' response rate (registry ctrl_frac), as counts on the panel's
    KCs x odors (rounded)."""
    from experiments import E83_cross_wiring_fit as E83
    P = E83.saved_fit(label)
    n_kc = P['saved']['on'][0].shape[0]
    n = n_kc * panel_input(panel).shape[1]
    f_off = P['target_off'] / P['size']
    f_on = float(H.data_target('ctrl_frac'))
    return int(round(f_off * n)), int(round(f_on * n)), float(f_off), f_on


def run_held(label: str, panel: str) -> dict:
    """One wiring on `panel`'s input with the response rates re-held (E77's `hold_rates`: the somatic threshold
    bisected with APL off to the wiring's own APL-off fraction, then the tied APL scales to 6.8 % with APL on; breadth
    held at the lab input's): hold_rates' dict, 'on' / 'off' with abbreviated columns."""
    from experiments import E77_input_screen as S
    from experiments import E83_cross_wiring_fit as E83
    P = E83.wiring_parts(label)
    t_off, t_on, _, _ = held_targets(label, panel)
    orn = panel_input(panel)
    res = S.hold_rates(lambda: E83.build_simulator(label, orn=orn, breadth=P['breadth'], params=dict(P['p'])),
                       target_off=t_off, target_on=t_on, params=dict(P['p']), label=f'{label} {panel}')
    res['on'], res['off'] = _abbrev(res['on']), _abbrev(res['off'])
    return res


def run_tag(panel: str, reading: str, wiring: str) -> str:
    return f'{panel}_{reading}_{wiring}'


def runs_dir() -> Path:
    return local_dir('B')


def save_run(tag: str, o: Dict[str, tuple]) -> None:
    d = runs_dir()
    for c, cc in (('on', 'on'), ('off', 'acute_off')):
        o[c][0].to_parquet(d / f'{tag}_{cc}_responses.parquet')
        o[c][1].to_parquet(d / f'{tag}_{cc}_spike_counts.parquet')


def load_run(tag: str) -> Optional[Dict[str, tuple]]:
    """A saved run ({'on', 'off'}: (responses, spike counts)) from the local outputs, else the model-output archive."""
    for d in (runs_dir(), H.shared_dir() / 'data' / 'model_outputs' / EXPERIMENT / 'B'):
        out = {}
        for c, cc in (('on', 'on'), ('off', 'acute_off')):
            ps = [d / f'{tag}_{cc}_{k}.parquet' for k in ('responses', 'spike_counts')]
            if not all(p.exists() for p in ps):
                break
            out[c] = tuple(pd.read_parquet(p) for p in ps)
        if len(out) == 2:
            return out
    return None


def _append_csv(path: Path, row: dict) -> None:
    df = pd.DataFrame([row])
    df.to_csv(path, mode='a', header=not path.exists(), index=False)


def runs_table() -> pd.DataFrame:
    p = out_dir() / 'runs_B.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _commits_clean() -> dict:
    from experiments import E89_apl_breadth_weighting as E89
    return E89._commits_clean()


def part_b(panels: Sequence[str] = PANELS, wirings: Optional[Sequence[str]] = None,
           readings: Sequence[str] = READINGS) -> None:
    """Runs and saves every wiring, panel and reading not yet in runs_B.csv (as-fitted runs first, each checked by
    `run_identity`, which must pass), with a runs_B.csv row (absolute parameters, the factors and targets, the response
    counts and rates, the identity, the process's commits, host, threads, the input's sha256)."""
    import gc
    import os
    from experiments import E83_cross_wiring_fit as E83
    require_lock()
    commits = _commits_clean()
    wirings = tuple(wirings or E83.ALL_WIRINGS)
    done = set(runs_table().get('tag', pd.Series(dtype=str)).astype(str))
    ident_path = out_dir() / 'identity_B.json'
    idents = _read_json(ident_path) or {}
    for reading in readings:
        for panel in panels:
            orn = panel_input(panel)
            for w in wirings:
                tag = run_tag(panel, reading, w)
                if tag in done and load_run(tag) is not None:
                    continue
                t0 = time.time()
                if reading == 'fitted':
                    res = run_fitted(w, panel)
                    ident = run_identity(w, res, panel)
                    idents[tag] = ident
                    _write_json(ident_path, idents)
                    if not ident['reproduces']:
                        raise RuntimeError(f'{tag}: the as-fitted run does not reproduce the saved fit: {ident}')
                    extra = dict(thr_factor=1.0, apl_factor=1.0, target_on=np.nan, target_off=np.nan, at_bound_thr=False,
                                 at_bound_apl=False, identity_reproduces=True)
                else:
                    res = run_held(w, panel)
                    extra = dict(thr_factor=res['thr_factor'], apl_factor=res['apl_factor'], target_on=res['target_on'],
                                 target_off=res['target_off'], at_bound_thr=res['at_bound_thr'],
                                 at_bound_apl=res['at_bound_apl'], identity_reproduces=np.nan)
                save_run(tag, res)
                P = E83.saved_fit(w)
                n = res['on'][0].size
                row = dict(tag=tag, panel=panel, reading=reading, wiring=w, fit_folder=P['folder'], fit_tag=P['tag'],
                           n_kcs=int(res['on'][0].shape[0]), n_odors=int(res['on'][0].shape[1]),
                           n_on=int((res['on'][0].to_numpy(dtype=float) > 0.5).sum()),
                           n_off=int((res['off'][0].to_numpy(dtype=float) > 0.5).sum()),
                           **{k: float(v) for k, v in res['params'].items()},
                           **{f'off_{k}': float(v) for k, v in res['params_off'].items()},
                           n_sims=int(res.get('n_sims', 0)), seconds=time.time() - t0, host=_host(),
                           time=time.strftime('%Y-%m-%d %H:%M'), threads=os.environ.get('OMP_NUM_THREADS'),
                           input_sha256=_frame_sha(orn), breadth='held at the lab input\'s (E83 wiring_parts)', **extra,
                           **{k: v for k, v in commits.items() if not isinstance(v, (dict, list))})
                row['rate_on'], row['rate_off'] = row['n_on'] / n, row['n_off'] / n
                _append_csv(out_dir() / 'runs_B.csv', row)
                print(f"{tag}: {row['seconds']:.0f} s, rate on {row['rate_on']:.4f}, off {row['rate_off']:.4f}", flush=True)
                del res
                gc.collect()


# ---------------------------------------------------------------------------
# the data pull

DATA_COPY = ('fly_counts_stim_rdm_concat.csv', 'stims_per_fly_long.csv', 'Odor-Panels.tex',
             'validation2/kc_soma/stim_to_drop_by_qcthresh.csv', 'figure-05/resampled/column_groups.yaml') + DISTANCE_FILES


def data_dir() -> Path:
    return H.shared_dir() / 'data' / 'remy_validation2_kc'


def copy_data() -> Path:
    """Copies the files read (DATA_COPY; not the per-trial file) to <shared>/data/remy_validation2_kc/ under their
    paths in Remy's folder, each checked against the lock's sha256, with a README.md of source paths, sizes and sha256;
    returns the folder. Existing copies with the same sha256 are kept."""
    import shutil
    lock = require_lock()
    locked = {f['path']: f for kind in ('manifests', 'distance_files') for f in lock['files'][kind]}
    dst = data_dir()
    dst.mkdir(parents=True, exist_ok=True)
    rows = []
    for rel in DATA_COPY:
        src = remy_dir() / rel
        out = dst / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        sha = _sha256(src)
        if sha != locked[rel]['sha256']:
            raise RuntimeError(f"{rel}: sha256 differs from the lock's")
        if not out.exists() or _sha256(out) != sha:
            shutil.copyfile(src, out)
        if _sha256(out) != sha:
            raise RuntimeError(f'{rel}: copy differs')
        rows.append(dict(path=rel, size=int(src.stat().st_size), sha256=sha))
    lines = [
        "# Remy's validation2 KC somata (and the megamat somata), copied for E107",
        '',
        f"Copied by experiments/E107_validation2_panel.py (`copy_data`) from {_host()} on {time.strftime('%Y-%m-%d %H:%M')}, "
        f"from `{remy_dir()}` (read only; nothing written there). Released for use by Tom on 2026-11-17 16:13 "
        "(chat.md, 16:13 hal: \"all that data may be used\").",
        '',
        "Each file keeps its path in Remy's folder. The per-trial file validation2/kc_soma/xrds_trial_rdm_concat.nc "
        '(41 MB) is not copied (not read). The lock of E107 (experiments/E107_part0_lock.json in model_experiment_queue) '
        'records the same sha256.',
        '',
        "What they are (Remy's manifest fly_counts_stim_rdm_concat.csv): the manuscript's per-fly validation2 KC soma "
        'correlation matrices (xrda_..._rmap_nls_cluster__correlation__0.3.nc, 20 flies, after her quality control); '
        'the same before quality control for three cell selections and two metrics (xrds_stim_rdm_concat.nc, '
        'variable Fc_zscore; stim_to_drop_by_qcthresh.csv is the drop list); her mean table of 22 odors '
        '(figure-05/resampled/tidy_abbrev_dists.tsv: KC, ORN, chemical, MetaCyc, LOTUS, VCF and older model distances; '
        'column_groups.yaml names the groups); her megamat KC somata per fly (18 flies; megamat17/kc_soma/'
        'xrds_stim_rdm_concat_mm.nc); the stimulus table (stims_per_fly_long.csv) and her odor panels (Odor-Panels.tex).',
        '',
        '| file | bytes | sha256 |',
        '|---|---:|---|',
    ] + [f"| {r['path']} | {r['size']} | {r['sha256']} |" for r in rows]
    (dst / 'README.md').write_text('\n'.join(lines) + '\n')
    return dst


# ---------------------------------------------------------------------------
# part A, no simulations

def sets() -> Dict[str, List[str]]:
    """{'136': its 17 odors, '231': the 22, '153o': the 18 odors outside the lock's rule (pathway_lca's complete set,
    descriptive only)}; the 136 is the carrying set when identified."""
    s = require_sets()
    out = {k: v['odors'] for k, v in s['sets'].items() if v.get('identified')}
    o18 = [v for v in s.get('outside_rule', {}).values() if len(v) == 18]
    if o18:
        out['153o'] = o18[0]
    return out


def carrying_set() -> str:
    return '136' if '136' in sets() else '231'


def set_pairs(key: str) -> List[Tuple[str, str]]:
    return H.pair_list(sets()[key])


def subset_of(pairs: Sequence[Tuple[str, str]], which: str) -> List[Tuple[str, str]]:
    """'all', 'pb' (pairs with paa or B-myr), 'four' (pairs with any of FOUR), 'nopb' (without paa and B-myr)."""
    if which == 'all':
        return list(pairs)
    if which == 'pb':
        return [p for p in pairs if p[0] in VERDICT_PAIR or p[1] in VERDICT_PAIR]
    if which == 'four':
        return [p for p in pairs if p[0] in FOUR or p[1] in FOUR]
    if which == 'nopb':
        return [p for p in pairs if p[0] not in VERDICT_PAIR and p[1] not in VERDICT_PAIR]
    raise ValueError(which)


def rho(model: pd.Series, data: pd.Series, pairs: Optional[Sequence[Tuple[str, str]]] = None) -> Tuple[float, int]:
    """ρ on odor-odor distances of `model` against `data` over `pairs` (E106's `rho_pairs`), and the number of pairs
    both hold."""
    from experiments import E106_rho_gap_parts as E106
    m = E106._canon(model)
    idx = E106.select_pairs(m.index, pairs)
    a, b = H.align_pairs(m.loc[idx], E106._canon(data))
    return E106.rho_pairs(model, data, pairs=pairs), int(len(a))


_CACHE: Dict[object, object] = {}


def remy_v2(kind: str = 'table_cos') -> pd.Series:
    """Remy's validation2 matrix: 'table_cos' / 'table_corr' (her mean table's kc-cosine- / kc-correlation-
    rmap_nls_cluster), 'flymean_cos' / 'flymean_corr' (the mean over flies of the lock's per-fly reading, pairs held by
    MIN_FLIES or more), 'table_orn_corr' / 'table_orn_cos' (her ORN columns, 5 flies: cross-check only)."""
    if kind not in _CACHE:
        col = {'table_cos': f'kc-cosine-{CELLS}', 'table_corr': f'kc-correlation-{CELLS}',
               'table_orn_corr': 'orn-correlation-01234', 'table_orn_cos': 'orn-cosine-01234'}
        if kind in col:
            _CACHE[kind] = mean_table(col[kind])
        elif kind in ('flymean_cos', 'flymean_corr'):
            _CACHE[kind] = fly_mean(per_fly('cosine' if kind.endswith('cos') else 'correlation'))
        else:
            raise ValueError(kind)
    return _CACHE[kind]


def remy_mm(metric: str = 'cosine') -> pd.Series:
    key = ('mm', metric)
    if key not in _CACHE:
        _CACHE[key] = remy_megamat(metric)
    return _CACHE[key][0]


def megamat_pairs(which: str = '136') -> List[Tuple[str, str]]:
    """The megamat panel's pairs: '136' (its 17 odors) or '66' (the 12 of Yang's panel)."""
    lock = require_lock()
    odors = lock['megamat17']['odors']
    return H.pair_list(odors if which == '136' else [o for o in odors if o in H.MEGAMAT12])


def identities() -> dict:
    """The readers' identities: the lock's per-fly correlation reading against the manuscript's per-fly file (equal),
    her mean table against the mean over flies (correlation, the lock's drop list: equal; cosine with the lock's drop
    list and with her cosine drop list), max |difference| and the pairs that differ."""
    out = {}
    a, b = per_fly('correlation'), manuscript_per_fly()
    a, b = a.reindex(index=b.index, columns=b.columns), b
    out['per_fly_correlation_vs_manuscript'] = dict(
        max_abs_diff=float(np.nanmax(np.abs(a.to_numpy() - b.to_numpy()))), nan_pattern_equal=bool((a.isna() == b.isna()).all().all()),
        n_flies=int(len(b)), n_pairs=int(b.shape[1]))
    for name, rows, col in (('correlation_lock_qc', per_fly('correlation'), 'table_corr'),
                            ('cosine_lock_qc', per_fly('cosine'), 'table_cos'),
                            ('cosine_cosine_qc', per_fly('cosine', qc='cosine'), 'table_cos')):
        fm = rows.mean()
        t = remy_v2(col)
        d = (fm.reindex(t.index) - t).dropna()
        out[f'table_vs_flymean_{name}'] = dict(max_abs_diff=float(d.abs().max()), mean_diff=float(d.mean()),
                                               n_pairs_over_1e9=int((d.abs() > 1e-9).sum()),
                                               pairs_over_0p01=[f'{p[0]} / {p[1]}' for p in d.index[d.abs() > 0.01]])
    return out


def _pct_rank(s: pd.Series, pair: Tuple[str, str]) -> float:
    """The pair's rank among the Series' pairs as a fraction (1 = the largest distance)."""
    return float(s.rank(pct=True).loc[pkey(*pair)])


def a_rhos() -> pd.DataFrame:
    """ρ on odor-odor distances of the no-simulation readers against Remy's matrices, per set and subset of pairs:
    the ORN input alone (1 - Pearson; cosine), her own ORN columns (cross-check), on validation2; on megamat, the
    megamat ORN input against her megamat somata (136, 66) and Yang's control flies against them (66)."""
    rows = []
    readers = {'orn': input_pairs('validation2', 'pearson'), 'orn_cosine': input_pairs('validation2', 'cosine'),
               'remy_orn_corr': remy_v2('table_orn_corr'), 'remy_orn_cos': remy_v2('table_orn_cos')}
    for rk, r in readers.items():
        for dk in ('table_cos', 'table_corr', 'flymean_cos', 'flymean_corr'):
            for sk in sets():
                for sub in ('all', 'pb', 'four', 'nopb'):
                    v, n = rho(r, remy_v2(dk), subset_of(set_pairs(sk), sub))
                    rows.append(dict(panel='validation2', reader=rk, remy=dk, set=sk, subset=sub, rho=v, n_pairs=n))
    mm_readers = {'orn': input_pairs('megamat', 'pearson'), 'orn_cosine': input_pairs('megamat', 'cosine'),
                  'yang_control': yang_control('cosine')}
    for rk, r in mm_readers.items():
        for metric in ('cosine', 'correlation'):
            for which in ('136', '66'):
                rr = yang_control(metric) if rk == 'yang_control' else r
                v, n = rho(rr, remy_mm(metric), megamat_pairs(which))
                rows.append(dict(panel='megamat', reader=rk, remy=f'flymean_{metric[:3] if metric == "cosine" else "corr"}',
                                 set=which, subset='all', rho=v, n_pairs=n))
    return pd.DataFrame(rows)


def a_verdict_pair() -> pd.DataFrame:
    """The distance of paa from B-myr in every measure: Remy's (mean table, mean over flies with the number of flies,
    the per-fly median and quartiles), with its rank among her pairs; the ORN input; Yang's control flies and APL output
    blocked (TNT); the family in Yang's panel (E106's population of the saved fits, APL on)."""
    from experiments import E106_rho_gap_parts as E106
    p = pkey(*VERDICT_PAIR)
    rows = []
    for k in ('table_cos', 'table_corr', 'flymean_cos', 'flymean_corr'):
        s = remy_v2(k)
        rows.append(dict(source='Remy, validation2 KC somata', reading=k, distance=float(s.loc[p]),
                         pct_rank_136=_pct_rank(s.loc[[q for q in set_pairs(carrying_set()) if q in s.index]], p),
                         pct_rank_231=_pct_rank(s, p)))
    for metric in ('cosine', 'correlation'):
        x = per_fly(metric)[p].dropna()
        rows.append(dict(source='Remy, validation2 KC somata, per fly', reading=metric, distance=float(x.median()),
                         q1=float(x.quantile(0.25)), q3=float(x.quantile(0.75)), n_flies=int(len(x))))
    for measure in ('pearson', 'cosine'):
        s = input_pairs('validation2', measure)
        rows.append(dict(source='ORN input (validation2, 36 glomeruli)', reading=measure, distance=float(s.loc[p]),
                         pct_rank_136=_pct_rank(s.loc[set_pairs(carrying_set())], p), pct_rank_231=_pct_rank(s, p)))
    for g, lab in (('TNT3c', 'control flies (Yang)'), ('TNT', 'APL output blocked (TNT) (Yang)')):
        s = H.canonical_pairs(H.yang_mean_rdm(g, 'cosine'))
        rows.append(dict(source=lab, reading='cosine', distance=float(s.loc[p]), pct_rank_120=_pct_rank(s, p)))
    pop = E106.population(E106.family_members('family'))
    s = H.canonical_pairs(pop['on'])
    rows.append(dict(source="the family in Yang's panel (saved fits, APL on)", reading='cosine on spike counts',
                     distance=float(s.loc[p]), pct_rank_120=_pct_rank(s, p)))
    return pd.DataFrame(rows)


def a_by_concentration() -> pd.DataFrame:
    """euc's and -aPine's pairs by concentration (the per-fly matrices, `by_conc`; each concentration's flies, pairs
    held by MIN_FLIES or more), pooled beside: per odor, concentration and metric, the flies, the mean distance to paa
    and to B-myr, the mean over its pairs, and ρ on odor-odor distances over its pairs with the ORN input alone (both
    measures) and against her pooled matrix."""
    rows = []
    lock = require_lock()
    conc = pd.DataFrame(lock['validation2']['concentrations'])
    inp = {m: input_pairs('validation2', m) for m in ('pearson', 'cosine')}
    for metric in ('cosine', 'correlation'):
        pooled = remy_v2('table_cos' if metric == 'cosine' else 'table_corr')
        for odor in ('euc', '-aPine'):
            bc = per_fly(metric, by_conc=(odor,))
            for c in sorted(POOLED[odor]) + [None]:
                name = odor if c is None else f'{odor} @ {c:g}'
                if c is None:
                    s = pooled
                    n_flies = int(conc[(conc.odor == odor)].n_in_rdm.sum())
                else:
                    cols = [p for p in bc.columns if name in p]
                    sub = bc[cols]
                    n = sub.notna().sum()
                    s = sub.mean()[n >= MIN_FLIES]
                    n_flies = int(conc[(conc.odor == odor) & np.isclose(conc.conc, c)].n_in_rdm.sum())
                    # back to the pooled name, so the pairs meet the input's (only `odor` is split)
                    s = pair_series({tuple(odor if x == name else x for x in p): v for p, v in s.items()})
                mine = [p for p in s.index if odor in p]
                d_paa = s.get(pkey(odor, 'paa'), np.nan) if odor != 'paa' else np.nan
                d_bm = s.get(pkey(odor, 'B-myr'), np.nan)
                r = dict(odor=odor, conc='pooled' if c is None else f'{c:g}', metric=metric, n_flies=n_flies,
                         n_pairs=len(mine), to_paa=float(d_paa), to_B_myr=float(d_bm),
                         mean_over_pairs=float(s.loc[mine].mean()) if mine else np.nan,
                         input_conc=float(conc[conc.odor == odor].input_conc.iloc[0]))
                for m, ip in inp.items():
                    r[f'rho_orn_{m}'] = rho(ip, s.loc[mine], mine)[0] if len(mine) > 2 else np.nan
                r['rho_vs_pooled'] = rho(s.loc[mine], pooled, mine)[0] if (c is not None and len(mine) > 2) else np.nan
                r['note'] = (f'7 of 8 flies with euc and 4 of 9 with -aPine had the odor at 10^-2, against an ORN input '
                             f'at 10^-1.75 and 10^-1.5')
                rows.append(r)
    return pd.DataFrame(rows)


def a_ceilings() -> pd.DataFrame:
    """The ceilings from halves of the flies (`split_half_ceiling`) of Remy's validation2 matrix (cosine and
    correlation, each set, all pairs and the pairs of paa and B-myr and of the four odors) and of her megamat somata (136
    and 66)."""
    rows = []
    for metric in ('cosine', 'correlation'):
        pf = per_fly(metric)
        for sk in sets():
            for sub in ('all', 'pb', 'four'):
                r = split_half_ceiling(pf, subset_of(set_pairs(sk), sub))
                rows.append(dict(panel='validation2', metric=metric, set=sk, subset=sub, **r))
        mm, info = megamat_somata(metric)
        for which in ('136', '66'):
            r = split_half_ceiling(mm, megamat_pairs(which))
            rows.append(dict(panel='megamat', metric=info['metric'], cells=info['cells'], set=which, subset='all', **r))
    return pd.DataFrame(rows)


def part_a() -> dict:
    """Part A: the identities, the ρ's of the readers, the distance of paa from B-myr, euc's and -aPine's pairs by
    concentration, the ceilings, the level check (a_*.csv, a_summary.json)."""
    out = out_dir()
    ident = identities()
    _write_json(out / 'a_identities.json', ident)
    tabs = dict(a_rhos=a_rhos(), a_verdict_pair=a_verdict_pair(), a_by_concentration=a_by_concentration(),
                a_ceilings=a_ceilings())
    for k, df in tabs.items():
        df.to_csv(out / f'{k}.csv', index=False)
    lv = {m: level_check(m) for m in ('cosine', 'correlation')}
    mm_info = {m: _CACHE[('mm', m)][1] for m in ('cosine', 'correlation') if ('mm', m) in _CACHE}
    summ = dict(identities=ident, level_check=lv, megamat_read=mm_info, sets=sets(), carrying_set=carrying_set(),
                registry=dict(rho_ctrl_remy_soma=H.data_target('rho_ctrl_remy_soma'),
                              rho_tnt_remy_soma=H.data_target('rho_tnt_remy_soma')))
    _write_json(out / 'a_summary.json', summ)
    provenance('A', level_check=lv)
    return summ


# ---------------------------------------------------------------------------
# part C, the runs read against Remy's matrices

def family_runs(panel: str, reading: str, wirings: Optional[Sequence[str]] = None) -> Dict[str, Dict[str, tuple]]:
    """{wiring: the saved run} of `panel` and `reading` (every wiring of the family by default; raises if one is
    missing)."""
    from experiments import E83_cross_wiring_fit as E83
    out = {}
    for w in (wirings or E83.ALL_WIRINGS):
        o = load_run(run_tag(panel, reading, w))
        if o is None:
            raise RuntimeError(f'{run_tag(panel, reading, w)}: no saved run')
        out[w] = o
    return out


def family_pairs(panel: str, reading: str) -> Dict[str, object]:
    """The family's pair distances (cosine on spike counts) on `panel`: 'pop' / 'pop_off' (the mean over the eight
    wirings, APL on / acutely off; E106's `mean_pairs`) and 'wiring' / 'wiring_off' ({wiring: pairs})."""
    from experiments import E106_rho_gap_parts as E106
    key = ('fam', panel, reading)
    if key not in _CACHE:
        runs = family_runs(panel, reading)
        per = {w: {c: H.canonical_pairs(H.pair_distances(o[c][1], 'cosine')) for c in ('on', 'off')} for w, o in runs.items()}
        pop = E106.mean_pairs(list(per.values()))
        _CACHE[key] = dict(pop=pop['on'], pop_off=pop['off'], wiring={w: d['on'] for w, d in per.items()},
                           wiring_off={w: d['off'] for w, d in per.items()})
    return _CACHE[key]


def quantities(remy: str = 'table', orn: str = 'pearson', reading: str = 'fitted', set_key: Optional[str] = None,
               mm_which: str = '136', mm_metric: str = 'cosine', rho_metric: str = 'cosine',
               member: str = 'pop') -> dict:
    """The outcome rules' quantities for one reading: the distance of paa from B-myr (cosine and correlation, Remy's
    `remy` matrix: 'table' or 'flymean'), and for the ORN input alone ('orn') and the family ('fam'; `member` 'pop' or a
    wiring) ρ against Remy's validation2 matrix (`rho_metric`) over the pairs of paa and B-myr ('pb') and the set
    ('all'), and against her megamat somata (`mm_metric`) over `mm_which` ('megamat')."""
    set_key = set_key or carrying_set()
    p = pkey(*VERDICT_PAIR)
    r_cos, r_corr = remy_v2(f'{remy}_cos'), remy_v2(f'{remy}_corr')
    data = r_cos if rho_metric == 'cosine' else r_corr
    sp = set_pairs(set_key)
    pb = subset_of(sp if all(o in sets()[set_key] for o in VERDICT_PAIR) else set_pairs('231'), 'pb')
    mm = remy_mm(mm_metric)
    mmp = megamat_pairs(mm_which)
    f_v2, f_mm = family_pairs('validation2', reading), family_pairs('megamat', reading)
    fam_v2 = f_v2['pop'] if member == 'pop' else f_v2['wiring'][member]
    fam_mm = f_mm['pop'] if member == 'pop' else f_mm['wiring'][member]
    o_v2, o_mm = input_pairs('validation2', orn), input_pairs('megamat', orn)
    return dict(
        d_cos=float(r_cos.loc[p]), d_corr=float(r_corr.loc[p]),
        orn=dict(pb=rho(o_v2, data, pb)[0], all=rho(o_v2, data, sp)[0], megamat=rho(o_mm, mm, mmp)[0]),
        fam=dict(pb=rho(fam_v2, data, pb)[0], all=rho(fam_v2, data, sp)[0], megamat=rho(fam_mm, mm, mmp)[0],
                 d_on=float(H.canonical_pairs(fam_v2).loc[p])),
        reading=dict(remy=remy, orn=orn, family=reading, set=set_key, megamat_pairs=mm_which, megamat_metric=mm_metric,
                     rho_metric=rho_metric, member=member, n_pb=len(pb), n_all=len(sp)))


PARTNER_KW = dict(held=dict(reading='held'), other_set=dict(set_key='231'), megamat66=dict(mm_which='66'),
                  fly_mean=dict(remy='flymean'), correlation=dict(rho_metric='correlation', mm_metric='correlation'),
                  orn_cosine=dict(orn='cosine'))


def verdict_core(level_flag: bool) -> dict:
    """The carrying class (`classify` on `quantities()`), each partner's (PARTNER_KW; 'other_set' is the 136 when the
    231 carries), each wiring's (descriptive), and whether the class holds (the lock's rule)."""
    from experiments import E83_cross_wiring_fit as E83
    carry_q = quantities()
    carry = classify(carry_q['d_cos'], carry_q['d_corr'], carry_q['orn'], carry_q['fam'], level_flag)
    partners = {}
    for name, kw in PARTNER_KW.items():
        kw = dict(kw)
        if name == 'other_set' and carrying_set() == '231':
            if '136' not in sets():
                continue
            kw['set_key'] = '136'
        q = quantities(**kw)
        partners[name] = dict(quantities=q, **classify(q['d_cos'], q['d_corr'], q['orn'], q['fam'], level_flag))
    per_wiring = {}
    for w in E83.ALL_WIRINGS:
        q = quantities(member=w)
        per_wiring[w] = dict(quantities=q, **classify(q['d_cos'], q['d_corr'], q['orn'], q['fam'], level_flag))
    differ = [k for k, v in partners.items() if v['class'] != carry['class']]
    level_ok = (not level_flag) or carry['class'] == carry['order_class']
    held = (not differ) and level_ok
    return dict(carrying=dict(quantities=carry_q, **carry), partners=partners, per_wiring=per_wiring,
                partners_differing=differ, level_ok=bool(level_ok), held=bool(held),
                outcome=carry['class'] if held else 'reading-dependent',
                n_wirings_same_class=int(sum(v['class'] == carry['class'] for v in per_wiring.values())))


def c_rhos() -> pd.DataFrame:
    """ρ on odor-odor distances of the family (population and each wiring), both readings, against Remy's validation2
    matrices per set and subset, and against her megamat somata (136, 66); the ORN input alone beside."""
    from experiments import E83_cross_wiring_fit as E83
    rows = []
    for reading in READINGS:
        fv, fm = family_pairs('validation2', reading), family_pairs('megamat', reading)
        members = [('pop', fv['pop'], fm['pop'])] + [(w, fv['wiring'][w], fm['wiring'][w]) for w in E83.ALL_WIRINGS]
        for mem, v2, mm in members:
            for dk in ('table_cos', 'table_corr', 'flymean_cos', 'flymean_corr'):
                for sk in sets():
                    for sub in ('all', 'pb', 'four', 'nopb'):
                        v, n = rho(v2, remy_v2(dk), subset_of(set_pairs(sk), sub))
                        rows.append(dict(panel='validation2', reader='family', reading=reading, member=mem, remy=dk,
                                         set=sk, subset=sub, rho=v, n_pairs=n))
            for metric in ('cosine', 'correlation'):
                for which in ('136', '66'):
                    v, n = rho(mm, remy_mm(metric), megamat_pairs(which))
                    rows.append(dict(panel='megamat', reader='family', reading=reading, member=mem,
                                     remy='flymean_cos' if metric == 'cosine' else 'flymean_corr', set=which, subset='all',
                                     rho=v, n_pairs=n))
    a = a_rhos()
    a = a[a.reader.isin(['orn', 'orn_cosine', 'yang_control'])].assign(reading='-', member='-')
    return pd.concat([pd.DataFrame(rows), a], ignore_index=True)


def c_shares() -> pd.DataFrame:
    """E61's shares by odor (E106's `gap_shares_by_odor`) of the ORN input alone and the family (both readings)
    against Remy's matrix (her mean table, cosine) over the carrying set, and against her megamat somata over the 136."""
    from experiments import E106_rho_gap_parts as E106
    rows = []
    sk = carrying_set()
    readers = [('orn', '-', input_pairs('validation2', 'pearson'), input_pairs('megamat', 'pearson'))]
    for reading in READINGS:
        readers.append(('family', reading, family_pairs('validation2', reading)['pop'], family_pairs('megamat', reading)['pop']))
    for name, reading, v2, mm in readers:
        for panel, model, data, pairs in (('validation2', v2, remy_v2('table_cos'), set_pairs(sk)),
                                          ('megamat', mm, remy_mm('cosine'), megamat_pairs('136'))):
            g = E106.gap_shares_by_odor(model, data, pairs=pairs).reset_index()
            g.insert(0, 'reader', name)
            g.insert(1, 'reading', reading)
            g.insert(2, 'panel', panel)
            rows.append(g)
    return pd.concat(rows, ignore_index=True)


def c_pb_pairs() -> pd.DataFrame:
    """The pairs of paa and B-myr in the carrying set: Remy's distances (cosine, correlation, flies), the ORN input
    alone, the family (both readings), and each one's rank residual against Remy's cosine matrix over the set (E106's
    `rank_residuals`)."""
    from experiments import E106_rho_gap_parts as E106
    lock = require_lock()
    nf = {(r['odor_a'], r['odor_b']): r['n_flies'] for r in lock['validation2']['flies_by_pair'] if r['conc_a'] == 'all'}
    sk = carrying_set()
    sp = set_pairs(sk)
    data = remy_v2('table_cos')
    readers = {'orn': input_pairs('validation2', 'pearson'),
               'family_fitted': family_pairs('validation2', 'fitted')['pop'],
               'family_held': family_pairs('validation2', 'held')['pop']}
    res = {k: E106.rank_residuals(v.loc[[p for p in sp if p in v.index]], data) for k, v in readers.items()}
    rows = []
    for p in subset_of(sp, 'pb'):
        r = dict(odor_a=p[0], odor_b=p[1], n_flies=nf.get(p, np.nan), remy_cos=float(data.loc[p]),
                 remy_corr=float(remy_v2('table_corr').loc[p]), remy_rank=float(res['orn'].loc[p, 'rank_data']))
        for k, v in readers.items():
            r[k] = float(v.loc[p])
            r[f'{k}_rank_resid'] = float(res[k].loc[p, 'rank_resid'])
        rows.append(r)
    return pd.DataFrame(rows)


def realism_rows(o: Dict[str, tuple]) -> dict:
    """A run's response rates (two or more spikes) with APL on and acutely off, their ratio, the shares of KCs by the
    number of odors they respond to (OPK_EDGES' classes, over the panel's odors) and the largest spike counts."""
    from pyolfsysm import kc_stats
    r_on, s_on = o['on']
    r_off, s_off = o['off']
    on = float((r_on.to_numpy(dtype=float) > 0.5).mean())
    off = float((r_off.to_numpy(dtype=float) > 0.5).mean())
    out = dict(rate_on=on, rate_off=off, ratio=off / on if on > 0 else np.nan,
               max_count_on=float(s_on.to_numpy(dtype=float).max()), max_count_off=float(s_off.to_numpy(dtype=float).max()),
               n_odors=int(r_on.shape[1]))
    labels = kc_stats.opk_class_labels()
    for c, r in (('on', r_on), ('off', r_off)):
        for lab, v in zip(labels, kc_stats.opk_class_shares(r.to_numpy(dtype=float) > 0.5)):
            out[f'opk_{c}_{lab}'] = float(v)
    return out


def c_realism() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(per run: `realism_rows` for each panel, reading and wiring, and the mean over wirings; per odor: the response
    fraction with APL on and off, the mean over wirings) with no band applied; the family's own fit panel and the KC
    data's values from the registry beside (control flies and APL output blocked (TNT); Remy's files hold distances
    only)."""
    from experiments import E83_cross_wiring_fit as E83
    rows, odor_rows = [], []
    for panel in PANELS:
        for reading in READINGS:
            runs = family_runs(panel, reading)
            per = []
            for w, o in runs.items():
                r = realism_rows(o)
                per.append(r)
                rows.append(dict(panel=panel, reading=reading, member=w, **r))
                for c in ('on', 'off'):
                    fr = (o[c][0].astype(float) > 0.5).mean()
                    for od, v in fr.items():
                        odor_rows.append(dict(panel=panel, reading=reading, member=w, condition=c, odor=od, fraction=float(v)))
            rows.append(dict(panel=panel, reading=reading, member='mean of 8', **pd.DataFrame(per).mean().to_dict()))
    for w in E83.ALL_WIRINGS:
        P = E83.saved_fit(w)
        rows.append(dict(panel="APL16 (the fit's panel)", reading='fitted', member=w, **realism_rows(P['saved'])))
    data = dict(panel='KC data (Yang, registry)', reading='-', member='control flies / APL output blocked (TNT)',
                rate_on=H.data_target('ctrl_frac'), rate_off=H.data_target('tnt_frac'),
                ratio=H.data_target('rate_ratio_fraction'), n_odors=16,
                **{f'opk_on_{k}': H.data_target(f'ctrl_odors_per_kc_{k}') for k in ('0', '1', '2_3', '4_7', '8plus')},
                **{f'opk_off_{k}': H.data_target(f'tnt_odors_per_kc_{k}') for k in ('0', '1', '2_3', '4_7', '8plus')})
    rows.append(data)
    od = pd.DataFrame(odor_rows)
    od_mean = od.groupby(['panel', 'reading', 'condition', 'odor'], sort=False).fraction.mean().reset_index().assign(member='mean of 8')
    return pd.DataFrame(rows), pd.concat([od, od_mean], ignore_index=True)


def distance_band(d: float) -> str:
    """'close' (CLOSE or less), 'apart' (APART or more) or 'middle'."""
    if d <= CLOSE:
        return 'close'
    if d >= APART:
        return 'apart'
    return 'middle'


def classify(d_cos: float, d_corr: float, orn: Dict[str, float], fam: Dict[str, float], level_flag: bool = False) -> dict:
    """The outcome class by the lock's rules (`lock_definitions()['outcome_rules']`): the band of the distance of paa
    from B-myr (cosine, 'middle' if the correlation distance's band differs), and from the ρ's of the ORN input alone
    and of the family ('pb', 'all', 'megamat') whether each orders the pairs within ORDER_TOL of its megamat level.
    Returns the class ('a' to 'd' or 'none'), the band, the class on the pairs' order alone and each condition."""
    b_cos, b_corr = distance_band(d_cos), distance_band(d_corr)
    band = b_cos if b_cos == b_corr else 'middle'
    orn_orders = bool(orn['pb'] >= orn['megamat'] - ORDER_TOL)
    fam_pb = bool(fam['pb'] >= fam['megamat'] - ORDER_TOL)
    fam_all = bool(fam['all'] >= fam['megamat'] - ORDER_TOL)
    if fam_pb and fam_all:
        order_class = 'd'
    elif orn_orders and not fam_pb:
        order_class = 'b'
    elif not orn_orders and not fam_pb:
        order_class = 'c'
    else:
        order_class = 'none'
    cls = 'a' if band == 'close' else order_class
    return {'class': cls, 'band': band, 'band_cos': b_cos, 'band_corr': b_corr, 'order_class': order_class,
            'on_order_alone': band == 'middle', 'orn_orders': orn_orders, 'fam_orders_pb': fam_pb,
            'fam_orders_all': fam_all, 'level_flag': bool(level_flag),
            'margins': dict(orn_pb=float(orn['pb'] - (orn['megamat'] - ORDER_TOL)),
                            fam_pb=float(fam['pb'] - (fam['megamat'] - ORDER_TOL)),
                            fam_all=float(fam['all'] - (fam['megamat'] - ORDER_TOL)))}


def verdict() -> dict:
    """verdict.json: the level check, the outcome class with every partner and wiring, the numbers it rests on, the
    level-adjusted distance (Remy's minus the level difference, beside), the ceilings of the carrying reading."""
    lv = level_check('cosine')
    core = verdict_core(lv['flagged'])
    ce = pd.read_csv(out_dir() / 'a_ceilings.csv') if (out_dir() / 'a_ceilings.csv').exists() else a_ceilings()
    sk = carrying_set()
    ceil = {f'{r.panel}_{r.metric}_{r.set}_{r.subset}': float(r.ceiling) for r in ce.itertuples()}
    q = core['carrying']['quantities']
    out = dict(outcome=core['outcome'], carrying_class=core['carrying']['class'], held=core['held'],
               partners_differing=core['partners_differing'], level_check=lv,
               level_adjusted_distance=float(q['d_cos'] - lv['difference']),
               carrying_set=sk, ceilings=ceil, core=core,
               rules=lock_definitions()['outcome_rules'], time=time.strftime('%Y-%m-%d %H:%M'), host=_host())
    _write_json(out_dir() / 'verdict.json', out)
    return out


def _f(x, spec: str = '.3f') -> str:
    try:
        return '-' if x is None or (isinstance(x, float) and np.isnan(x)) else format(float(x), spec)
    except (TypeError, ValueError):
        return str(x)


def _r(tab: pd.DataFrame, **sel) -> float:
    m = np.ones(len(tab), dtype=bool)
    for k, v in sel.items():
        m &= (tab[k].astype(str) == str(v)).to_numpy()
    x = tab[m]
    return float(x.rho.iloc[0]) if len(x) else np.nan


def comparison_table() -> str:
    """comparison_table.md from part A's and C's saved tables: ρ on odor-odor distances of the ORN input alone, the
    family as fitted and re-held, and Remy's flies' ceilings side by side; the distance of paa from B-myr; the class
    with its partners; euc's and -aPine's pairs by concentration; the realism numbers; the identities."""
    out = out_dir()
    ct = pd.read_csv(out / 'c_rhos.csv')
    ce = pd.read_csv(out / 'a_ceilings.csv')
    vp = pd.read_csv(out / 'a_verdict_pair.csv')
    bc = pd.read_csv(out / 'a_by_concentration.csv')
    rl = pd.read_csv(out / 'c_realism.csv')
    v = _read_json(out / 'verdict.json')
    ident = _read_json(out / 'a_identities.json')
    sk = v['carrying_set']
    other = '231' if sk == '136' else '136'

    def ceil(panel, metric, s, sub='all'):
        x = ce[(ce.panel == panel) & (ce.metric == metric) & (ce.set.astype(str) == s) & (ce.subset == sub)]
        return float(x.ceiling.iloc[0]) if len(x) else np.nan

    L = ['# E107: the four validation2 odors in the panel their ORN input was recorded in', '',
         "Generated by experiments/E107_validation2_panel.py (`comparison_table`) from results/E107_validation2_panel/. "
         "ρ is ρ on odor-odor distances (Spearman over pairs) against Remy's validation2 KC somata (her mean table, "
         "cosine distance of the manuscript's cells, rmap_nls_cluster, unless said); the megamat level is against her "
         f"megamat somata (mean over her 18 flies, cosine, {v['level_check']['cells']}). The family is the mean matrix "
         "over its eight wirings, APL on, cosine on spike counts. Remy's flies are wild type, so only APL on is read.", '',
         f"## The outcome: {v['outcome']}", '',
         f"Carrying class ({v['carrying_class']}); held: {v['held']}; partners giving another class: "
         f"{', '.join(v['partners_differing']) or 'none'}; wirings giving the carrying class: "
         f"{v['core']['n_wirings_same_class']} of 8.", '',
         '| reading | class | band | ORN orders paa / B-myr pairs | family orders them | family orders the set |',
         '|---|---|---|---|---|---|']
    c = v['core']['carrying']
    L.append(f"| carrying | {c['class']} | {c['band']} | {c['orn_orders']} | {c['fam_orders_pb']} | {c['fam_orders_all']} |")
    for k, p in v['core']['partners'].items():
        L.append(f"| {k} | {p['class']} | {p['band']} | {p['orn_orders']} | {p['fam_orders_pb']} | {p['fam_orders_all']} |")
    lv = v['level_check']
    L += ['', f"Level check: her megamat somata average {_f(lv['remy_mean'])} against Yang's control flies' "
          f"{_f(lv['yang_mean'])} on the {lv['n_pairs']} shared pairs (cosine; difference {_f(lv['difference'], '+.3f')}, "
          f"{'flagged' if lv['flagged'] else 'within'} the {LEVEL_TOL} tolerance). The distance of paa from B-myr shifted "
          f"by that difference: {_f(v['level_adjusted_distance'])}.", '']
    # rho table
    L += [f'## ρ on odor-odor distances, side by side ({sk} pairs carrying)', '',
          '| statistic | ORN input alone | the family, as fitted | the family, re-held | Remy\'s flies (ceiling) |',
          '|---|---:|---:|---:|---:|']
    rows = [(f'validation2, the {sk} pairs', 'validation2', sk, 'all'), (f'validation2, the {other} pairs', 'validation2', other, 'all'),
            (f'validation2, the pairs of paa and B-myr ({sk})', 'validation2', sk, 'pb'),
            (f'validation2, the pairs of the four odors ({sk})', 'validation2', sk, 'four'),
            (f'validation2, the {sk} without paa and B-myr', 'validation2', sk, 'nopb'),
            ('validation2, the 153 outside the rule (pathway_lca\'s 18 odors)', 'validation2', '153o', 'all')]
    for lab, panel, s, sub in rows:
        o = _r(ct, panel=panel, reader='orn', remy='table_cos', set=s, subset=sub)
        f1 = _r(ct, panel=panel, reader='family', reading='fitted', member='pop', remy='table_cos', set=s, subset=sub)
        f2 = _r(ct, panel=panel, reader='family', reading='held', member='pop', remy='table_cos', set=s, subset=sub)
        L.append(f'| {lab} | {_f(o)} | {_f(f1)} | {_f(f2)} | {_f(ceil(panel, "cosine", s, sub if sub != "nopb" else "all"))} |')
    for which in ('136', '66'):
        o = _r(ct, panel='megamat', reader='orn', remy='flymean_cos', set=which)
        f1 = _r(ct, panel='megamat', reader='family', reading='fitted', member='pop', remy='flymean_cos', set=which)
        f2 = _r(ct, panel='megamat', reader='family', reading='held', member='pop', remy='flymean_cos', set=which)
        L.append(f'| megamat level, the {which} pairs | {_f(o)} | {_f(f1)} | {_f(f2)} | {_f(ceil("megamat", "cosine", which))} |')
    y = _r(ct, panel='megamat', reader='yang_control', remy='flymean_cos', set='66')
    L += ['', f"Yang's control flies against her megamat somata on the 66 pairs: ρ={_f(y)} on odor-odor distances "
          f"(cosine both; the registry's rho_ctrl_remy_soma, Remy's correlation distance, "
          f"{_f(H.data_target('rho_ctrl_remy_soma'))}).", '']
    pw = pd.DataFrame([dict(wiring=w, pb=d['quantities']['fam']['pb'], all=d['quantities']['fam']['all'],
                            megamat=d['quantities']['fam']['megamat'], d_on=d['quantities']['fam']['d_on'], cls=d['class'])
                       for w, d in v['core']['per_wiring'].items()])
    L += ['Per wiring (as fitted; ρ against her cosine matrix):', '',
          '| wiring | ρ, pairs of paa and B-myr | ρ, the set | megamat level | paa / B-myr distance | class |',
          '|---|---:|---:|---:|---:|---|'] + [f'| {r.wiring} | {_f(r.pb)} | {_f(r.all)} | {_f(r.megamat)} | {_f(r.d_on)} | {r.cls} |'
                                               for r in pw.itertuples()]
    # the distance of paa from B-myr
    L += ['', '## The distance of paa from B-myr', '', '| source | reading | distance | rank among the pairs | flies |',
          '|---|---|---:|---:|---:|']
    for r in vp.itertuples():
        rk = r.pct_rank_136 if pd.notna(getattr(r, 'pct_rank_136', np.nan)) else getattr(r, 'pct_rank_120', np.nan)
        L.append(f"| {r.source} | {r.reading} | {_f(r.distance)} | {_f(rk, '.2f')} | {_f(getattr(r, 'n_flies', np.nan), '.0f')} |")
    for reading in READINGS:
        q = quantities(reading=reading)
        L.append(f"| the family on the validation2 input, {reading} | cosine on spike counts | {_f(q['fam']['d_on'])} | - | - |")
    L += ['', "Rank among the pairs: the pair's percentile among the carrying set's pairs (1 = the largest distance; Yang's "
          'and the family in her panel: among the 120 pairs).', '']
    # by concentration
    L += ['## euc and -aPine by concentration', '',
          '7 of her 8 flies with euc and 4 of her 9 with -aPine had the odor at 10^-2, against an ORN input recorded at '
          '10^-1.75 and 10^-1.5. The verdict does not rest on them.', '',
          '| odor | concentration | metric | flies | to paa | to B-myr | ρ with the ORN input over its pairs |',
          '|---|---|---|---:|---:|---:|---:|']
    for r in bc.itertuples():
        L.append(f'| {r.odor} | {r.conc} | {r.metric} | {r.n_flies} | {_f(r.to_paa)} | {_f(r.to_B_myr)} | {_f(r.rho_orn_pearson)} |')
    # realism
    L += ['', '## Realism (no band applied)', '',
          '| panel | reading | response rate, APL on | APL off | ratio | KCs to no odor, APL on | 8 or more, APL on | largest count, APL on |',
          '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in rl[rl.member.isin(['mean of 8', 'control flies / APL output blocked (TNT)'])].itertuples():
        L.append(f"| {r.panel} | {r.reading} | {_f(r.rate_on, '.4f')} | {_f(r.rate_off, '.4f')} | {_f(r.ratio, '.2f')} | "
                 f"{_f(getattr(r, 'opk_on_0'), '.3f')} | {_f(getattr(r, 'opk_on_8plus'), '.3f')} | {_f(r.max_count_on, '.0f')} |")
    fit = rl[rl.panel.str.startswith('APL16')]
    L.append(f"| APL16 (the fit's panel), mean of 8 | fitted | {_f(fit.rate_on.mean(), '.4f')} | {_f(fit.rate_off.mean(), '.4f')} | "
             f"{_f(fit.ratio.mean(), '.2f')} | {_f(fit.opk_on_0.mean(), '.3f')} | {_f(fit.opk_on_8plus.mean(), '.3f')} | "
             f"{_f(fit.max_count_on.mean(), '.0f')} |")
    L += ['', "The odors-per-KC classes count over each panel's odors (22, 17, 16), so they are not comparable across "
          "panels as a distribution; the KC data's are Yang's (16 odors). Remy's files hold distances only.", '']
    # identities
    L += ['## Identities and readings', '',
          f"- The lock's per-fly correlation reading equals the manuscript's per-fly file: max |difference| "
          f"{_f(ident['per_fly_correlation_vs_manuscript']['max_abs_diff'], '.2g')}, the same missing pairs.",
          f"- Her mean table's correlation column is the mean over flies of that reading (max |difference| "
          f"{_f(ident['table_vs_flymean_correlation_lock_qc']['max_abs_diff'], '.2g')}). Its cosine column is the mean over "
          f"flies with her cosine drop list (max |difference| {_f(ident['table_vs_flymean_cosine_cosine_qc']['max_abs_diff'], '.2g')}), "
          f"not the lock's (the manuscript's correlation list; max |difference| "
          f"{_f(ident['table_vs_flymean_cosine_lock_qc']['max_abs_diff'], '.3f')}, "
          f"{ident['table_vs_flymean_cosine_lock_qc']['n_pairs_over_1e9']} pairs, those of bbenz and Z4-7al): "
          "the partner 'fly_mean' reads the lock's.",
          f"- The 136 is the MetaCyc columns' 17 odors (the lock's rule); the 153 is not identified by the rule "
          "(pathway_lca alone gives 18 odors, the 17 and sab: shown as '153 outside the rule').",
          '- Every as-fitted run reproduces the saved fit on the odors its panel shares with the fit panel (0 differences; '
          'identity_B.json).']
    text = '\n'.join(L) + '\n'
    (out / 'comparison_table.md').write_text(text)
    return text


def figure() -> Path:
    """E107_validation2.png: Remy's validation2 distances (cosine) against the ORN input alone and the family as
    fitted over the carrying set, the pairs of paa and B-myr marked; the ρ's beside their megamat levels."""
    plt = H._plt()
    sk = carrying_set()
    sp = set_pairs(sk)
    data = remy_v2('table_cos')
    fam = family_pairs('validation2', 'fitted')['pop']
    orn = input_pairs('validation2', 'pearson')
    pb = set(subset_of(sp, 'pb'))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (lab, m, xl) in zip(axs[:2], (('ORN input alone', orn, 'ORNs: Pearson correlation distance'),
                                          ('the family, as fitted, APL on', fam, 'model: cosine distance'))):
        x = np.array([m.loc[p] for p in sp])
        y = np.array([data.loc[p] for p in sp])
        is_pb = np.array([p in pb for p in sp])
        ax.scatter(x[~is_pb], y[~is_pb], s=12, c='0.6', label='other pairs')
        ax.scatter(x[is_pb], y[is_pb], s=18, c='tab:red', label='pairs of paa or B-myr')
        i = sp.index(pkey(*VERDICT_PAIR))
        ax.scatter([x[i]], [y[i]], s=70, facecolors='none', edgecolors='k', label='paa / B-myr')
        r_all, r_pb = rho(m, data, sp)[0], rho(m, data, list(pb))[0]
        ax.set_title(f'{lab}\n(ρ on odor-odor distances {r_all:.2f}; paa / B-myr pairs {r_pb:.2f})', fontsize=10)
        ax.set_xlabel(xl)
        ax.set_ylabel("KCs, Remy's wild-type flies: cosine distance")
        ax.legend(fontsize=8, loc='lower right')
    ax = axs[2]
    q = quantities()
    labels = ['ORN input\nalone', 'family,\nas fitted']
    vals_pb = [q['orn']['pb'], q['fam']['pb']]
    vals_all = [q['orn']['all'], q['fam']['all']]
    vals_mm = [q['orn']['megamat'], q['fam']['megamat']]
    xs = np.arange(2)
    ax.bar(xs - 0.25, vals_all, 0.25, label=f'validation2, the {sk} pairs')
    ax.bar(xs, vals_pb, 0.25, label='validation2, pairs of paa and B-myr')
    ax.bar(xs + 0.25, vals_mm, 0.25, label='megamat level (136 pairs)')
    for xx, mm in zip(xs, vals_mm):
        ax.plot([xx - 0.4, xx + 0.4], [mm - ORDER_TOL] * 2, 'k--', lw=1)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1)
    ax.set_ylabel("ρ on odor-odor distances with Remy's flies")
    ax.set_title(f'the pairs of paa and B-myr against the megamat level\n(dashed: the level minus {ORDER_TOL})', fontsize=10)
    ax.legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    p = out_dir() / 'E107_validation2.png'
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def part_c() -> dict:
    """Part C: c_rhos.csv, c_shares.csv, c_pb_pairs.csv, c_realism.csv, c_rates_by_odor.csv, verdict.json,
    comparison_table.md, the figure, provenance_C.json (with the run parameters as absolute values)."""
    out = out_dir()
    if not (out / 'a_ceilings.csv').exists():
        part_a()
    c_rhos().to_csv(out / 'c_rhos.csv', index=False)
    c_shares().to_csv(out / 'c_shares.csv', index=False)
    c_pb_pairs().to_csv(out / 'c_pb_pairs.csv', index=False)
    rl, ro = c_realism()
    rl.to_csv(out / 'c_realism.csv', index=False)
    ro.to_csv(out / 'c_rates_by_odor.csv', index=False)
    v = verdict()
    comparison_table()
    figure()
    rb = runs_table()
    provenance('C', model=H.model_provenance(orn_deltas=panel_input('validation2'), orn_deltas_megamat=panel_input('megamat')),
               outcome=v['outcome'], runs=rb.to_dict(orient='records'),
               input_sha256={p: _frame_sha(panel_input(p)) for p in PANELS},
               part_files=sorted(x.name for x in out.glob('provenance_*.json')))
    (out / 'provenance.json').write_text((out / 'provenance_C.json').read_text())
    return v


def write_variants() -> pd.DataFrame:
    """variants.csv (the archive's params.csv source): one row per part B run with its absolute parameters, the input,
    the factors, the response rates and the largest counts. No model_pair_distances_long.csv is written: the runs are on
    other panels than Yang's, which kcq/apl_effect.py scores (E98's choice)."""
    rb = runs_table()
    rows = []
    for _, r in rb.iterrows():
        o = load_run(r.tag)
        re_ = realism_rows(o)
        rows.append(dict(tag=r.tag, variant=r.tag, part='B', base='bouton', k=2.0, breadth_k=2.0, wiring=r.wiring,
                         config=f'{r.panel}_{r.reading}', orn_input=f'{r.panel}_orn_deltas on the 36 glomeruli of the lab input',
                         kind=(f"C2's family on {r.wiring} with the {r.panel} ORN input ({r.n_odors} odors); breadth held at "
                               f"the lab input's; " + ('the fit\'s parameters, nothing re-tuned' if r.reading == 'fitted' else
                                                      "the somatic threshold and the tied APL scale re-held at the wiring's "
                                                      'APL-off fraction and 6.8 % (E77 hold_rates)')),
                         ref_experiment=r.fit_folder, ref_fit=r.fit_tag, fixed_thr=r.fixed_thr, wAPLKC_scale=r.wAPLKC_scale,
                         wKCAPL_scale=r.wKCAPL_scale, wAPLPN_scale=r.wAPLPN_scale, wPNAPL_scale=r.wPNAPL_scale,
                         wKCKC_scale=r.wKCKC_scale, kc_kc_inh=r.wKCKC_scale, thr_factor=r.thr_factor, apl_factor=r.apl_factor,
                         rate_on=re_['rate_on'], rate_off=re_['rate_off'], rate_acute_off=re_['rate_off'],
                         max_count_on=re_['max_count_on'], max_count_off=re_['max_count_off'], n_kcs=r.n_kcs,
                         stimuli=r.n_odors, subdir='B'))
    df = pd.DataFrame(rows)
    df.to_csv(out_dir() / 'variants.csv', index=False)
    provenance('variants', n_variants=len(df))
    return df


def provenance(part: str, **extra) -> None:
    """provenance_<part>.json: script, host, time, `H.code_commits()`, the lock's and the sets' sha256, Remy's files."""
    lock = require_lock()
    s = _read_json(SETS_PATH) or {}
    _write_json(out_dir() / f'provenance_{part}.json', dict(
        script='experiments/E107_validation2_panel.py', part=part, host=_host(), time=time.strftime('%Y-%m-%d %H:%M'),
        **H.code_commits(), lock_sha256=lock['sha256'], lock_file_sha256=_sha256(LOCK_PATH),
        sets_file_sha256=_sha256(SETS_PATH) if SETS_PATH.exists() else None, sets=s.get('sets'),
        remy_files=lock['files'], **extra))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--part', required=True, choices=('0', 'sets', 'data', 'A', 'B', 'C', 'variants'))
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--panels', default=','.join(PANELS))
    ap.add_argument('--readings', default=','.join(READINGS))
    ap.add_argument('--wirings', default='', help='comma-separated wiring labels (default all eight)')
    args = ap.parse_args()
    if args.part == 'B':
        part_b(panels=[x for x in args.panels.split(',') if x], wirings=[x for x in args.wirings.split(',') if x] or None,
               readings=[x for x in args.readings.split(',') if x])
        return
    if args.part == '0':
        lock = write_lock(force=args.force)
        print(f"-> {LOCK_PATH} (sha256 {lock['sha256'][:12]}); commit it before any other part", flush=True)
        return
    if args.part == 'sets':
        s = identify_sets()
        for k, v in s['sets'].items():
            print(k, v, flush=True)
        print(f'-> {SETS_PATH}; commit it before part A', flush=True)
        return
    if args.part == 'data':
        print(f'-> {copy_data()}', flush=True)
        return
    if args.part == 'C':
        v = part_c()
        print(f"outcome {v['outcome']} (carrying {v['carrying_class']}; partners differing {v['partners_differing']})", flush=True)
        print(f'-> {out_dir()}', flush=True)
        return
    if args.part == 'variants':
        print(f'{len(write_variants())} variants -> {out_dir() / "variants.csv"}', flush=True)
        return
    if args.part == 'A':
        s = part_a()
        print(json.dumps(dict(level_check=s['level_check'], carrying_set=s['carrying_set']), indent=1), flush=True)
        print(f'-> {out_dir()}', flush=True)
        return
    raise NotImplementedError(args.part)


if __name__ == '__main__':
    main()
