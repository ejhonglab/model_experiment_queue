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


def pooled_name(name: str, conc: float, by_conc: bool = False) -> str:
    """A stimulus's odor name in a matrix: the name, or for the pooled odors read by concentration '<name> @ <conc>'."""
    return f'{name} @ {conc:g}' if by_conc and name in POOLED else name


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


def require_lock(check_git: bool = True, check_files: bool = True) -> dict:
    """The committed lock, after checking that its sha256 is its contents', its definitions are the module's
    (`lock_definitions`), it is committed and unchanged in git (`check_git`) and Remy's files still have the locked
    sha256 (`check_files`); raises PermissionError otherwise (the reading code calls this before opening any
    distance)."""
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
              by_conc: bool = False, exclude: Sequence[str] = EXCLUDED_STIMULI) -> pd.DataFrame:
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


def per_fly(measure: str = 'cosine', cells: str = CELLS, by_conc: bool = False, qc: str = 'lock') -> pd.DataFrame:
    """Remy's per-fly validation2 matrices (xrds_stim_rdm_concat.nc, Fc_zscore, at `cells` and `measure`, with the
    lock's drop list, `qc`) as flies x pairs (NaN where a fly lacks a pair); with `by_conc`, euc, -aPine and ger keep
    their concentration in the odor name."""
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
    raise NotImplementedError(args.part)


if __name__ == '__main__':
    main()
