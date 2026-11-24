#!/usr/bin/env python3
"""E134: the recorded populations of every stage and panel, in one interface, with the lock.

A population is the table of trial-mean responses of all units of one stage (glomeruli, bouton regions, claw regions,
KCs) to the odors of one panel, one table per fly (CLAUDE.md, "The main goal for new experiments"). This module loads
every population E134's entry names (item 2 of experiment_priorities.md, design steps 0 and 1), in one layout, and
writes and checks the lock: the files with their sha256, the flies, the odor names, the exclusions, each stage's
response criterion and units, the breadth bins scaled to the panel's size, and for the ORN stage that the recordings
are of axon terminals.

The populations (`POPULATIONS`; keys (stage, panel, cohort)) and where each comes from (`SOURCES`):
- megamat (17 odors at 10^-3): ORNs (olfdb dataset 1, 9 flies), PN dendrites (GH146, al_analysis's
  data/sent_to_remy/2023-10-29, 7 flies), PN boutons (5 flies) and KC claws (4 flies) from data/remy_per_roi/
  (E119 part 0; the bouton fly 2022-10-21 fly02 left out), KC somata of Remy's four 2022 flies (her suite2p trials and
  respvec files, read from her processed folders on matrix and copied, `--copy`), and the one fly recorded at claws and
  somata on the same trials (cohort 'same_fly', data/remy_per_roi/same_fly_trials.npz);
- validation2 (22 odors): ORNs (olfdb dataset 4, 6 flies), KC somata (Remy's 20 flies of the manuscript, her
  per-acquisition respvec files in analysis_outputs/data_by_acq_2/kc_soma_nls, copied);
- the kiwi and control mixtures, two panels of 8 stimuli each ('kmix', 'cmix'; every KC fly had one of the two, every
  ORN fly both): ORNs (olfdb dataset 5, 9 flies), KC somata (olfdb dataset 8: 7 kmix flies, 4 cmix flies);
- the APL16 panel: KC somata of the control flies (cohort 'control', 11 recordings) and of APL output blocked (TNT)
  ('tnt', 8 recordings), data/yang_apl16_kc_noise_units/.

Units (`Population.units`) and criteria (`CRITERIA`), one per stage:
- 'baseline-noise SD' where the recordings' frames are on hand: the window mean minus the pre-odor mean, over the SD
  of the unit's frames within the pre-odor windows (each window's mean removed, pooled over trials). Boutons and claws
  2.0 (E119), KC somata 1.9684 (the realism bands' criterion, notes/response_criterion_decision.md).
- 'trial-noise SD' for the glomeruli (ORNs, PN dendrites), whose processed tables hold one dF/F per trial and no
  frames: the trial mean over a robust SD of the unit's single trials about their odor means (per fly and glomerulus,
  1.4826 times the median absolute residual over the panel's trials, each residual scaled by sqrt(n / (n - 1)) for n
  trials of the odor). Criterion 2.0. The lock's `noise_bridge` gives, on the stages that have both, the ratio of this
  unit to the baseline-noise unit.
- 'session z (Fc_zscore)' for the mixtures' KC somata (olfdb dataset 8 holds trial means of Remy's whole-session
  z-score and no frames or trials): criterion `MIXTURE_Z_CRITERION`, the z at which Remy's four megamat flies have the
  response rate their baseline-noise criterion gives (`z_bridge`, computed by `--copy` and recorded in the lock).

Every population is built by `--build` into one long parquet file, <shared>/data/stage_populations/populations/
<stage>__<panel>__<cohort>.parquet (columns fly, unit, odor, repeat, value; repeat -1 for a trial mean where the source
has no trials), which `load_population` reads. The lock (results/E134_population_baseline/lock.json and its committed
copy experiments/E134_lock.json) holds each such file's sha256 and every source's; `load_population` refuses a file
whose sha256 differs from the lock's.

Usage (from the repository root):
    nice -n 10 ./venv/bin/python -u data_prep/stage_population_data.py --copy    # copy and reduce what the share lacks
    ./venv/bin/python -u data_prep/stage_population_data.py --build              # the population files and checks
    ./venv/bin/python -u data_prep/stage_population_data.py --lock               # the lock (then commit it)
    ./venv/bin/python -u data_prep/stage_population_data.py --check              # every population against the lock
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import socket
import subprocess
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

STAGES = ('orn', 'pn_dendrite', 'pn_bouton', 'kc_claw', 'kc_soma')
PANELS = ('megamat', 'validation2', 'kmix', 'cmix', 'apl16')
MIXTURE_PANELS = ('kmix', 'cmix')       # "the kiwi / control mixtures" of the entry: two panels
# what one unit of each stage is; only glomeruli are named (the fourth view, "the matrix itself", reads them)
UNIT_KIND = {'orn': 'glomerulus', 'pn_dendrite': 'glomerulus', 'pn_bouton': 'bouton_region',
             'kc_claw': 'claw_region', 'kc_soma': 'kc'}
NAMED_UNITS = {s: UNIT_KIND[s] == 'glomerulus' for s in STAGES}
RECORDING_SITE = {'orn': 'axon terminals (pebbled > GCaMP6f, antennal lobe glomeruli)',
                  'pn_dendrite': 'dendrites (GH146 > GCaMP6f, antennal lobe glomeruli)',
                  'pn_bouton': 'boutons (calyx, per suite2p ROI)', 'kc_claw': 'claws (calyx, per suite2p ROI)',
                  'kc_soma': 'somata (nuclear GCaMP6s, per cell)'}
# (stage, panel, cohort): cohort 'wild_type' (Remy's and the lab's flies), 'same_fly' (the one megamat fly recorded at
# claws and somata on the same trials), 'control' / 'tnt' (Yang's APL16 flies)
POPULATIONS: Tuple[Tuple[str, str, str], ...] = (
    ('orn', 'megamat', 'wild_type'),
    ('pn_dendrite', 'megamat', 'wild_type'),
    ('pn_bouton', 'megamat', 'wild_type'),
    ('kc_claw', 'megamat', 'wild_type'),
    ('kc_soma', 'megamat', 'wild_type'),
    ('kc_claw', 'megamat', 'same_fly'),
    ('kc_soma', 'megamat', 'same_fly'),
    ('orn', 'validation2', 'wild_type'),
    ('kc_soma', 'validation2', 'wild_type'),
    ('orn', 'kmix', 'wild_type'),
    ('kc_soma', 'kmix', 'wild_type'),
    ('orn', 'cmix', 'wild_type'),
    ('kc_soma', 'cmix', 'wild_type'),
    ('kc_soma', 'apl16', 'control'),
    ('kc_soma', 'apl16', 'tnt'),
)
# the bouton fly that every analysis leaves out (Remy's dashboards do; data/remy_per_roi/README.md)
EXCLUDED_FLIES = ({'stage': 'pn_bouton', 'panel': 'megamat', 'fly': '2022-10-21__fly02__megamat1_001__dsub3',
                   'reason': "Remy's dashboards and bouton analysis leave it out (median cellprob 0.03); "
                             'data/remy_per_roi/README.md'},)

# the panels' odors (abbreviations; the mixtures' stimuli as olfdb names them, with their log10 concentration)
MEGAMAT_ODORS = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms', 'B-cit',
                 'Lin', 'aa', 'benz', 'va')
VALIDATION2_ODORS = ('paa', '-aPine', 'euc', 'B-myr', '+pul', '-bCar', '1-3ol', '1o3one', '1p3one', '2-mba', '2-mib',
                     'EtOct', 'PEA', 'Z4-7al', 'bbenz', 'geos', 'ger', 'guai', 'long', 'mchav', 'menth', 'sab')
# the mixture panels' stimuli as the model's input names them (`natmix_orn_deltas`' columns through `H.odor_abbrev`):
# the components (kmix: EtOH 10^-2, IAol 10^-3.6, IaA 10^-3.7, ea 10^-4.2, eb 10^-3.5; cmix: ms 10^-3, va 10^-3, fur
# 10^-4, 2h 10^-5, 1o3ol 10^-3) and the mixture at three dilutions (kmix0 = olfdb's 'kmix @ 0', kmix-1, kmix-2)
KMIX_ODORS = ('EtOH', 'IAol', 'IaA', 'ea', 'eb', 'kmix0', 'kmix-1', 'kmix-2')
CMIX_ODORS = ('ms', 'va', 'fur', '2h', '1o3ol', 'cmix0', 'cmix-1', 'cmix-2')
APL16_ODORS = ('2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol', '1-5ol', '1-6ol', 'ms', 'paa', '-aPine',
               'euc', 'B-myr')
PANEL_ODORS = {'megamat': MEGAMAT_ODORS, 'validation2': VALIDATION2_ODORS, 'kmix': KMIX_ODORS, 'cmix': CMIX_ODORS,
               'apl16': APL16_ODORS}

# one response criterion per stage, in the stage's units; the mixtures' somata (session z) take MIXTURE_Z_CRITERION
CRITERIA = {'orn': 2.0, 'pn_dendrite': 2.0, 'pn_bouton': 2.0, 'kc_claw': 2.0, 'kc_soma': 1.9684}
UNITS_BASELINE = 'baseline-noise SD'
UNITS_TRIAL = 'trial-noise SD'
UNITS_SESSION_Z = 'session z (Fc_zscore)'
MAD_SCALE = 1.4826                       # the MAD of a normal variable times this is its SD
# response and pre-odor windows (s after odor onset; closed)
KC_WINDOW_REMY = {'response': (2.0, 8.0), 'pre': (-10.0, -0.5)}     # her validation2 respvec (data_by_acq_2)
KC_WINDOW_YANG = {'response': (2.0, 8.0), 'pre': (-5.0, 0.0)}       # data/yang_apl16_kc_noise_units/README.md

def _shared() -> Path:
    from kcq import harness as H
    try:
        return H.shared_dir()
    except FileNotFoundError:                                   # no share mounted (tests on made-up tables)
        return Path('/local/tensor/Tom/claude_kc_remapping')


SHARED = _shared()                                      # H.shared_dir(): /local/tensor/... or /mnt/tensor/...
HERE = Path(__file__).resolve().parents[1]
RESULTS = SHARED / 'results' / 'E134_population_baseline'
LOCK_PATH = RESULTS / 'lock.json'                       # the lock, in the shared directory ...
LOCK_COPY = HERE / 'experiments' / 'E134_lock.json'     # ... and its committed copy (no model is read before it)
DATA = SHARED / 'data' / 'stage_populations'            # what `--copy` copied and reduced, and `--build`'s files
SOURCES_DIR = DATA / 'sources'
POP_DIR = DATA / 'populations'
OLFDB = 'http://hal:8080'
OLFDB_DATASETS = {'orn_megamat': 1, 'orn_validation2': 4, 'orn_mixtures': 5, 'kc_mixtures': 8}
REMY_TENSOR = SHARED.parents[1] / 'Remy' / 'odor_space_collab' / 'analysis_outputs'
REMY_V2_DIR = REMY_TENSOR / 'data_by_acq_2' / 'kc_soma_nls'
# Remy's processed folders on matrix (her suite2p outputs and respvec files of the 2022 megamat flies; the same share
# E119 part 0b read the same fly's suite2p output from, at viki's user's instruction of 2026-11-20), read only
REMY_MATRIX = Path('/local/matrix/Remy-Data/projects/odor_space_collab/processed_data')
REMY_MEGAMAT_ACQS = ('2022-10-10/1/megamat0', '2022-10-10/2/megamat0', '2022-10-11/1/megamat0',
                     '2022-11-10/1/megamat0__dsub03')
V2_CELL_MASK = 'rmap_nls_cluster'        # the manuscript's cells (E107's CELLS)
V2_QC = dict(cell_mask_coord='rmap_nls_cluster', metric='correlation', qc_thresh=0.3)   # E107's QC_LIST
V2_EXCLUDED_STIMULI = ('MethOH @ -2.0', 'pfo @ 0.0', 'geos @ -5.0')                     # E107's
V2_NAMES = {'PAA': 'paa', '1-prop': '1-3ol'}            # E107's odor map; every other name unchanged
ORN_V2_NAMES = {'1-3ol': '1-3ol'}


# ---------------------------------------------------------------------------
# the population

@dataclass
class Population:
    """One stage's recorded population on one panel: per fly, a units x odors table of trial means.

    Attributes:
        stage, panel, cohort: the key (one of `POPULATIONS`).
        odors: the panel's odors (abbreviations; the mixtures' stimuli with their concentration), in panel order;
            every table's columns are these, NaN where a fly was not given the odor (or Remy's quality control dropped
            it).
        flies: the flies (recordings) kept, in order; excluded flies are not loaded.
        tables: fly -> DataFrame, index the units (glomerulus names where `named_units`, else the recording's ROI or
            cell numbers as strings), columns `odors`, values the trial-mean response in `units`.
        criterion: the response criterion in `units`: a unit responds to an odor if its trial mean is above it.
        units: `UNITS_BASELINE`, `UNITS_TRIAL` or `UNITS_SESSION_Z` (the module docstring).
        graded: False when only binary responses are on hand (no population is binary at present).
        trials: fly -> DataFrame (index the units, columns a MultiIndex (odor, repeat)) of single-trial responses in
            the same units, or None when not loaded or not on hand (`load_population(trials=True)`).
        has_trials: whether single trials are on hand (False for the mixtures' KC somata, olfdb's trial means only:
            their measurement noise cannot be read from them).
        unit_kind: `UNIT_KIND[stage]`.
        named_units: True for glomeruli.
        recording_site: where the indicator was imaged (`RECORDING_SITE`; the ORNs are their axon terminals).
        window: the response window and the pre-odor window in s after odor onset, as (start, end) pairs.
        exclusions: what was left out (flies, units, odors, trials), one dict per item with a reason.
        files: the population file and its sources, as recorded in the lock (path, bytes, sha256).
        notes: anything a reader must know (concentration substitutions, how a fly's odors were pooled).
    """
    stage: str
    panel: str
    cohort: str
    odors: Tuple[str, ...]
    flies: Tuple[str, ...]
    tables: Dict[str, pd.DataFrame]
    criterion: float
    units: str
    graded: bool = True
    trials: Optional[Dict[str, pd.DataFrame]] = None
    has_trials: bool = True
    unit_kind: str = ''
    named_units: bool = False
    recording_site: str = ''
    window: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    exclusions: List[dict] = field(default_factory=list)
    files: List[dict] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def key(self) -> Tuple[str, str, str]:
        return (self.stage, self.panel, self.cohort)

    def responses(self, fly: str) -> pd.DataFrame:
        """The fly's units x odors table of responses (1.0 above the criterion, 0.0 at or below, NaN where the
        odor was not given)."""
        t = self.tables[fly]
        return (t > self.criterion).astype(float).where(t.notna())

    def breadth_bins(self, fly: Optional[str] = None) -> Tuple[Tuple[int, int], ...]:
        """`breadth_bins` of the panel's size, or of the number of odors `fly` was given (validation2's flies were
        each given part of the panel)."""
        return breadth_bins(len(self.odors) if fly is None else len(self.fly_odors(fly)))

    def n_units(self) -> pd.Series:
        """Units per fly."""
        return pd.Series({f: len(t) for f, t in self.tables.items()}, name='n_units')

    def fly_odors(self, fly: str) -> Tuple[str, ...]:
        """The odors the fly was given (columns not all NaN)."""
        t = self.tables[fly]
        return tuple(o for o in self.odors if t[o].notna().any())


def population_keys() -> Tuple[Tuple[str, str, str], ...]:
    """Every (stage, panel, cohort) this module loads (`POPULATIONS`)."""
    return POPULATIONS


def panel_odors(panel: str) -> Tuple[str, ...]:
    """The panel's odors in panel order (`PANEL_ODORS`): megamat 17; validation2 22; the mixtures' 8 stimuli each,
    as olfdb names them ('kmix @ -2', 'IaA @ -3.7', ...); APL16 16 (the solvent left out)."""
    return PANEL_ODORS[panel]


def response_criterion(stage: str, panel: str, cohort: str = 'wild_type') -> float:
    """The stage's response criterion in its units: `CRITERIA[stage]`, except the mixtures' KC somata (session z),
    whose criterion is the z bridge recorded by `--copy` (`mixture_z_criterion`)."""
    if stage == 'kc_soma' and panel in MIXTURE_PANELS:
        return mixture_z_criterion()
    return CRITERIA[stage]


def _round_half_up(x: float) -> int:
    return int(np.floor(x + 0.5))


def breadth_bins(n_odors: int) -> Tuple[Tuple[int, int], ...]:
    """The bins of the number of odors a unit responds to, scaled to the panel's size: 0, 1, 2-3, 4-7, 8-16 on 16
    odors (the realism band's bins), with the upper three edges at 2/16, 4/16 and 8/16 of the panel (rounded half up)
    and kept distinct (each at least one above the last, the first at least 2). Returns closed (lo, hi) pairs covering
    0 to n_odors: 17 odors 0, 1, 2-3, 4-8, 9-17; 22 odors 0, 1-2, 3-5, 6-10, 11-22; 8 odors 0, 1, 2, 3, 4-8."""
    n = int(n_odors)
    if n < 4:
        raise ValueError(f'breadth bins need at least 4 odors, got {n}')
    e1 = max(2, _round_half_up(2 * n / 16))
    e2 = max(e1 + 1, _round_half_up(4 * n / 16))
    e3 = max(e2 + 1, _round_half_up(8 * n / 16))
    e3 = min(e3, n)
    return ((0, 0), (1, e1 - 1), (e1, e2 - 1), (e2, e3 - 1), (e3, n))


def breadth_labels(n_odors: int) -> Tuple[str, ...]:
    """Labels of `breadth_bins(n_odors)` ('0', '1', '2-3', ...)."""
    return tuple(str(lo) if lo == hi else f'{lo}-{hi}' for lo, hi in breadth_bins(n_odors))


# ---------------------------------------------------------------------------
# small helpers

def _host() -> str:
    return socket.gethostname().split('.')[0]


def _now() -> str:
    return time.strftime('%Y-%m-%d %H:%M')


def sha256_file(path: Path, block: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(block), b''):
            h.update(b)
    return h.hexdigest()


def rel_path(path: Path) -> str:
    """A path relative to the shared directory ('<shared>/...'), the repository ('<queue>/...') or as given, so that
    a lock written on one host reads on another (atlas mounts the share at /mnt/tensor, hal and viki at /local/tensor)."""
    path = Path(path)
    for root, tag in ((SHARED, '<shared>'), (HERE, '<queue>')):
        try:
            return f'{tag}/{path.resolve().relative_to(root.resolve()).as_posix()}'
        except ValueError:
            continue
    return str(path)


def file_entry(path: Path, sha: bool = True) -> dict:
    path = Path(path)
    st = path.stat()
    return dict(path=rel_path(path), bytes=int(st.st_size), mtime=time.strftime('%Y-%m-%d %H:%M:%S',
                time.localtime(st.st_mtime)), sha256=sha256_file(path) if sha else None)


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def _read_json(path: Path) -> Optional[dict]:
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def pop_name(key: Tuple[str, str, str]) -> str:
    return '__'.join(key)


def pop_path(key: Tuple[str, str, str], pop_dir: Optional[Path] = None) -> Path:
    return (POP_DIR if pop_dir is None else Path(pop_dir)) / f'{pop_name(key)}.parquet'


def robust_trial_sd(trials: np.ndarray) -> np.ndarray:
    """The trial-noise SD of each unit: trials is (repeats x units x odors), NaN where a trial or odor is missing.
    Each trial's residual from its odor's mean over the unit's trials, times sqrt(n / (n - 1)) for the n trials of that
    odor (so that it has the trial's own SD), then 1.4826 times the median absolute residual over the unit's odors and
    trials (NaN for a unit with no odor of two or more trials)."""
    X = np.asarray(trials, dtype=float)
    n = np.sum(np.isfinite(X), axis=0)                          # units x odors
    import warnings
    with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        m = np.nanmean(X, axis=0)
        f = np.sqrt(n / (n - 1.0))
        r = (X - m[None]) * f[None]
    r[:, n < 2] = np.nan
    a = np.abs(r).transpose(1, 0, 2).reshape(X.shape[1], -1)  # units x (repeats * odors)
    with np.errstate(all='ignore'):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return MAD_SCALE * np.nanmedian(a, axis=1)


def long_from_trials(fly: str, units: Sequence[str], odors: Sequence[str], X: np.ndarray) -> pd.DataFrame:
    """Long rows (fly, unit, odor, repeat, value) of a (repeats x units x odors) array, NaN entries dropped."""
    R, U, O = X.shape
    df = pd.DataFrame(dict(fly=fly, unit=np.tile(np.repeat(np.asarray(units, dtype=str), O), R),
                           odor=np.tile(np.asarray(odors, dtype=object), R * U),
                           repeat=np.repeat(np.arange(R), U * O), value=X.reshape(-1)))
    return df[np.isfinite(df.value)].reset_index(drop=True)


def trials_from_long(df: pd.DataFrame, odors: Sequence[str]) -> Tuple[np.ndarray, List[str]]:
    """(repeats x units x odors) array and the unit labels of one fly's long rows (repeat >= 0)."""
    units = list(pd.unique(df.unit))
    reps = sorted(df.repeat.unique())
    X = np.full((len(reps), len(units), len(odors)), np.nan)
    ui = {u: i for i, u in enumerate(units)}
    oi = {o: i for i, o in enumerate(odors)}
    ri = {r: i for i, r in enumerate(reps)}
    sub = df[df.odor.isin(list(odors))]
    X[sub.repeat.map(ri).to_numpy(), sub.unit.map(ui).to_numpy(), sub.odor.map(oi).to_numpy()] = sub.value.to_numpy()
    return X, units


# ---------------------------------------------------------------------------
# sources: olfdb

def olfdb_download(dataset: int, url: str = OLFDB) -> pd.DataFrame:
    """olfdb's long download of one dataset (`/api/download?ids=<id>&format=long`)."""
    import io
    import urllib.request
    with urllib.request.urlopen(f'{url}/api/download?ids={dataset}&format=long', timeout=300) as r:
        return pd.read_csv(io.BytesIO(r.read()))


def olfdb_datasets(url: str = OLFDB) -> List[dict]:
    import urllib.request
    with urllib.request.urlopen(f'{url}/api/datasets', timeout=60) as r:
        return json.loads(r.read().decode())['datasets']


def olfdb_path(dataset: int) -> Path:
    return SOURCES_DIR / 'olfdb' / f'ds{dataset}_long.csv'


def _fly_label(date, fly_num) -> str:
    return f'{pd.Timestamp(date).date()}/{int(fly_num)}'


def orn_long(dataset: int, panel: str) -> pd.DataFrame:
    """One olfdb glomerulus dataset's trials of `panel` (odor2 'solvent' only, the panel's odors), in trial-noise
    units (`robust_trial_sd` per fly and glomerulus over the panel's trials): long rows fly, unit, odor, repeat, value,
    and the noise SD per row's unit in `noise_sd` (dF/F)."""
    d = pd.read_csv(olfdb_path(dataset))
    olf_panel = {'megamat': 'megamat', 'validation2': 'validation2', 'kmix': 'kmix', 'cmix': 'cmix'}[panel]
    d = d[d.panel == olf_panel]
    if 'odor2' in d.columns:
        d = d[d.odor2 == 'solvent']
    if panel in MIXTURE_PANELS:
        d = d.assign(odor=[_olfdb_mix_name(o) for o in d.odor1])
    else:
        d = d.assign(odor=[V2_NAMES.get(o.split(' @ ')[0], o.split(' @ ')[0]) for o in d.odor1])
    odors = PANEL_ODORS[panel]
    d = d[d.odor.isin(odors)]
    rows = []
    for (date, fly), sub in d.groupby(['date', 'fly_num'], sort=True):
        lab = _fly_label(date, fly)
        sub = sub.assign(fly=lab, unit=sub.roi.astype(str))
        X, units = trials_from_long(sub[['fly', 'unit', 'odor', 'repeat', 'value']], odors)
        sd = robust_trial_sd(X)
        keep = np.isfinite(sd) & (sd > 0)
        out = long_from_trials(lab, np.asarray(units)[keep], odors, X[:, keep] / sd[keep][None, :, None])
        out['noise_sd'] = out.unit.map(dict(zip(np.asarray(units)[keep], sd[keep])))
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def _olfdb_mix_name(o: str) -> str:
    """olfdb's mixture stimuli as the model's input names them (KMIX_ODORS / CMIX_ODORS): a component by its
    abbreviation ('IaA @ -3.7' -> 'IaA'), a mixture with its dilution ('kmix @ -2' -> 'kmix-2', 'kmix @ 0' -> 'kmix0')."""
    name, conc = o.split(' @ ')
    if name in ('kmix', 'cmix'):
        return f'{name}{float(conc):g}'
    return name


def mixture_kc_long(panel: str) -> pd.DataFrame:
    """olfdb dataset 8 (Remy's kiwi / control mixture KCs, trial means of Fc_zscore) of `panel`: long rows fly, unit,
    odor, repeat (-1: trial means), value (session z)."""
    d = pd.read_csv(olfdb_path(OLFDB_DATASETS['kc_mixtures']))
    d = d[d.panel == panel]
    d = d.assign(fly=[_fly_label(a, b) for a, b in zip(d.date, d.fly_num)], unit=d.roi.astype(str),
                 odor=[_olfdb_mix_name(o) for o in d.odor1], repeat=-1)
    d = d[d.odor.isin(PANEL_ODORS[panel])]
    return d[['fly', 'unit', 'odor', 'repeat', 'value']].reset_index(drop=True)


# ---------------------------------------------------------------------------
# sources: the PN dendrites (GH146) and the pebbled ORNs of al_analysis's sent_to_remy tables

def gh146_dir() -> Path:
    from kcq import harness as H
    return H.al_analysis_data_dir() / 'sent_to_remy' / '2023-10-29'


def gh146_copy(kind: str) -> Path:
    return SOURCES_DIR / 'gh146' / f'{kind}_ij_certain-roi_stats.p'


def al_trials_long(kind: str) -> pd.DataFrame:
    """The megamat trials (odor2 'solvent') of the sent_to_remy table `kind` ('GH146' / 'pebbled'), as E90's
    `al_trials` reads them (the copy in sources/gh146/): long rows fly, unit (glomerulus), odor, repeat, value (dF/F)."""
    df = pd.read_pickle(gh146_copy(kind))
    m = df.loc['megamat']
    m = m[m.index.get_level_values('odor2') == 'solvent']
    m = m.droplevel(['is_pair', 'odor2'])
    s = m.stack(level=list(range(m.columns.nlevels)), future_stack=True) if hasattr(m, 'stack') else None
    s = s.dropna()
    s.name = 'value'
    s = s.reset_index()
    s['odor'] = [o.split(' @ ')[0] for o in s['odor1']]
    s['fly'] = [_fly_label(a, b) for a, b in zip(s['date'], s['fly_num'])]
    s['unit'] = s['roi'].astype(str)
    s = s[s.odor.isin(MEGAMAT_ODORS)]
    return s[['fly', 'unit', 'odor', 'repeat', 'value']].reset_index(drop=True)


def pn_dendrite_long() -> pd.DataFrame:
    """The GH146 PN dendrites in trial-noise units (as `orn_long`)."""
    d = al_trials_long('GH146')
    rows = []
    for lab, sub in d.groupby('fly', sort=True):
        X, units = trials_from_long(sub, MEGAMAT_ODORS)
        sd = robust_trial_sd(X)
        keep = np.isfinite(sd) & (sd > 0)
        out = long_from_trials(lab, np.asarray(units)[keep], MEGAMAT_ODORS, X[:, keep] / sd[keep][None, :, None])
        out['noise_sd'] = out.unit.map(dict(zip(np.asarray(units)[keep], sd[keep])))
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# sources: Remy's per-ROI boutons and claws (E119 part 0) and the same fly

def roi_dir() -> Path:
    return SHARED / 'data' / 'remy_per_roi'


def per_roi_long(stage: str) -> pd.DataFrame:
    """The boutons ('pn_bouton') or claws ('kc_claw') in baseline-noise units, as E119's `fly_trials` reads them
    (window 'e90', 0.5 to 3 s; F for the boutons, Fc for the claws), the excluded fly out: long rows fly, unit (ROI),
    odor, repeat, value."""
    from experiments import E119_stages_per_roi as E119
    st = {'pn_bouton': 'bouton', 'kc_claw': 'claw'}[stage]
    rows = []
    for acq, d in E119.fly_trials(st, odors=MEGAMAT_ODORS).items():
        if any(acq == x['fly'] or x['fly'].startswith(acq) for x in EXCLUDED_FLIES):
            continue
        rows.append(long_from_trials(acq, d['rois'].astype(str), MEGAMAT_ODORS, d['noise']))
    return pd.concat(rows, ignore_index=True)


def same_fly_long(stage: str) -> pd.DataFrame:
    """The same fly's claws or somata (data/remy_per_roi/same_fly_trials.npz, baseline-noise units)."""
    z = np.load(roi_dir() / 'same_fly_trials.npz', allow_pickle=True)
    odors = [str(o) for o in z['odors']]
    key = {'kc_claw': 'claw', 'kc_soma': 'soma'}[stage]
    X = z[f'{key}_noise']
    out = long_from_trials('2023-04-11/3/megamat1_redo', z[f'{key}_cells'].astype(str), odors, X)
    return out[out.odor.isin(MEGAMAT_ODORS)].reset_index(drop=True)


# ---------------------------------------------------------------------------
# sources: Yang's APL16 KC somata

def apl16_dir() -> Path:
    return SHARED / 'data' / 'yang_apl16_kc_noise_units'


def apl16_long(cohort: str) -> pd.DataFrame:
    """Yang's APL16 KC somata in baseline-noise units (per_trial_long.parquet; the solvent out): long rows fly
    ('<date>/<fly_num>/<hemibrain>', a recording), unit (olfdb's roi), odor, repeat, value."""
    g = {'control': 'TNT3c', 'tnt': 'TNT'}[cohort]
    p = pd.read_parquet(apl16_dir() / 'per_trial_long.parquet',
                        columns=['genotype', 'date', 'fly_num', 'hemibrain', 'roi', 'odor', 'repeat', 'value'])
    p = p[p.genotype == g]
    odor = p.odor.astype(str).str.split(' @ ').str[0]
    p = p.assign(odor=odor)
    p = p[p.odor.isin(APL16_ODORS)]
    fly = p.date.astype(str) + '/' + p.fly_num.astype(int).astype(str) + '/' + p.hemibrain.astype(str)
    return pd.DataFrame(dict(fly=fly.to_numpy(), unit=p.roi.astype(int).astype(str).to_numpy(), odor=p.odor.to_numpy(),
                             repeat=p.repeat.astype(int).to_numpy(), value=p.value.astype(float).to_numpy()))


# ---------------------------------------------------------------------------
# sources: Remy's validation2 KC somata (data_by_acq_2, her respvec files, copied)

def v2_copy_dir() -> Path:
    return SOURCES_DIR / 'remy_validation2_kc_soma'


def v2_acqs(root: Optional[Path] = None) -> List[str]:
    """The 20 acquisitions of Remy's aggregate validation2_ijk_all_2 (df_acqs.csv), as folder names: from the copy
    (sources/remy_validation2_kc_soma/aggregate/df_acqs.csv) when `root` is None and it exists, else from `root`
    (default her folder, REMY_V2_DIR)."""
    copy = v2_copy_dir() / 'aggregate' / 'df_acqs.csv'
    if root is None and copy.exists():
        f = copy
    else:
        root = REMY_V2_DIR if root is None else Path(root)
        f = root / 'aggregate' / 'df_acqs.csv'
        if not f.exists():
            f = root / 'aggregate' / 'validation2_ijk_all_2' / 'df_acqs.csv'
    a = pd.read_csv(f, sep='\t')
    return [f'{d}__fly{int(n):02d}__{t}' for d, n, t in zip(a.date_imaged, a.fly_num, a.thorimage_name)]


V2_FILES = ('respvec/xrds_peak_amp.nc', 'respvec/xrds_baseline_std.nc', 'stim_rdm/xrds_stim_rdm.nc')


def v2_mask(coord: np.ndarray) -> np.ndarray:
    """Remy's cell mask from a cell coordinate (the rule her stim_rdm files reproduce, `check_v2_rdms`): a cell is in
    the mask if the coordinate is finite and nonzero."""
    c = np.asarray(coord, dtype=float)
    return np.isfinite(c) & (c != 0)


def _v2_stims(stims: Sequence[str]) -> Tuple[List[Optional[str]], List[str]]:
    """Remy's stimulus names -> the panel's odor names (None for an excluded stimulus); also her names as given."""
    out = []
    for s in stims:
        s = str(s)
        name, conc = s.split(' @ ')
        if f'{name} @ {float(conc):.1f}' in V2_EXCLUDED_STIMULI or s in V2_EXCLUDED_STIMULI:
            out.append(None)
            continue
        out.append(V2_NAMES.get(name, name))
    return out, [str(s) for s in stims]


def v2_qc_drops() -> pd.DataFrame:
    """Remy's quality-control drop list at V2_QC (the share's copy, E107's): fly ('date/fly_num'), stim."""
    q = pd.read_csv(SHARED / 'data' / 'remy_validation2_kc' / 'validation2' / 'kc_soma' /
                    'stim_to_drop_by_qcthresh.csv', sep='\t')
    q = q[(q.cell_mask_coord == V2_QC['cell_mask_coord']) & (q.metric == V2_QC['metric']) &
          np.isclose(q.qc_thresh.astype(float), float(V2_QC['qc_thresh']))]
    return pd.DataFrame(dict(fly=[f'{d}/{int(n)}' for d, n in zip(q.date_imaged, q.fly_num)],
                             stim=q.stim_to_drop.astype(str).to_numpy())).drop_duplicates()


def v2_soma_long(copy_dir: Optional[Path] = None) -> Tuple[pd.DataFrame, List[dict]]:
    """Remy's validation2 KC somata in baseline-noise units, from her respvec files (peak_amp of Fc, the 2 to 8 s mean
    minus the -10 to -0.5 s mean; baseline_std of Fc, the SD of each trial's -10 to -0.5 s frames): per cell, the
    response over sqrt(mean over trials of baseline_std^2). Cells of her manuscript's mask (V2_CELL_MASK), stimuli
    mapped to the panel's names (E107's map, pooled concentrations under one name), the excluded stimuli and her
    quality-control drops (V2_QC) out. Returns the long rows (fly 'date/fly_num') and the exclusions."""
    import xarray as xr
    d = v2_copy_dir() if copy_dir is None else Path(copy_dir)
    drops = v2_qc_drops()
    rows, excl = [], []
    for acq in v2_acqs():
        p = xr.open_dataset(d / acq / 'respvec' / 'xrds_peak_amp.nc')
        b = xr.open_dataset(d / acq / 'respvec' / 'xrds_baseline_std.nc')
        date, fly = str(p.attrs['acq.date_imaged']), int(p.attrs['acq.fly_num'])
        lab = f'{date}/{fly}'
        mask = v2_mask(p[V2_CELL_MASK].values)
        amp = p['Fc'].values[:, mask]
        sd = np.sqrt(np.nanmean(b['Fc'].values[:, mask] ** 2, axis=0))
        cells = p['cells'].values[mask]
        names, stims = _v2_stims(p['stim'].values)
        dropped = set(drops[drops.fly == lab].stim)
        occ = p['stim_occ'].values
        ok_cell = np.isfinite(sd) & (sd > 0)
        if (~ok_cell).any():
            excl.append(dict(fly=lab, what='cells', n=int((~ok_cell).sum()), reason='zero or undefined baseline SD'))
        for k, (nm, st) in enumerate(zip(names, stims)):
            if nm is None:
                continue
            if st in dropped:
                continue
            rows.append(pd.DataFrame(dict(fly=lab, unit=cells[ok_cell].astype(str), odor=nm, repeat=int(occ[k]),
                                          value=amp[k, ok_cell] / sd[ok_cell])))
        for st in sorted(dropped & set(stims)):
            excl.append(dict(fly=lab, what='stimulus', stim=st, reason=f"Remy's quality control ({V2_QC})"))
        p.close()
        b.close()
    return pd.concat(rows, ignore_index=True), excl


def check_v2_rdms(copy_dir: Optional[Path] = None, metric: str = 'cosine', tol: float = 1e-6) -> dict:
    """Our per-fly distances over the mask reproduce Remy's stim_rdm files (xrds_stim_rdm.nc, cell_mask_coord
    V2_CELL_MASK, `metric`) from the respvec trial means of Fc_zscore: the largest absolute difference per fly."""
    import xarray as xr
    from scipy.spatial.distance import cdist
    d = v2_copy_dir() if copy_dir is None else Path(copy_dir)
    out = {}
    for acq in v2_acqs():
        p = xr.open_dataset(d / acq / 'respvec' / 'xrds_peak_amp.nc')
        r = xr.open_dataset(d / acq / 'stim_rdm' / 'xrds_stim_rdm.nc')
        mask = v2_mask(p[V2_CELL_MASK].values)
        X = p['Fc_zscore'].to_pandas().iloc[:, mask]
        X.index = [str(s) for s in p['stim'].values]
        M = X.groupby(level=0).mean()
        k = np.flatnonzero((r['cell_mask_coord'].values == V2_CELL_MASK) & (r['metric'].values == metric))
        if len(k) != 1:
            raise ValueError(f'{acq}: {len(k)} rdms at {V2_CELL_MASK}, {metric}')
        R = r['Fc_zscore'].isel(rdms=int(k[0]))
        stims_a = [str(s) for s in R['stim_row'].values]
        stims_b = [str(s) for s in R['stim_col'].values]
        A = M.reindex(stims_a).to_numpy()
        B = M.reindex(stims_b).to_numpy()
        ours = cdist(A, B, metric='cosine' if metric == 'cosine' else 'correlation')
        diff = float(np.nanmax(np.abs(ours - R.values)))
        out[acq] = dict(max_abs_diff=diff, n_rois_remy=int(r['n_rois'].values[k[0]]), n_rois_ours=int(mask.sum()))
        p.close()
        r.close()
    return dict(per_fly=out, passed=bool(all(v['max_abs_diff'] <= tol and v['n_rois_remy'] == v['n_rois_ours']
                                             for v in out.values())), metric=metric, tol=tol)


# ---------------------------------------------------------------------------
# sources: Remy's four 2022 megamat flies (matrix; copied and reduced by `--copy`)

def megamat_copy_dir() -> Path:
    return SOURCES_DIR / 'remy_megamat_kc_soma'


def megamat_acq_label(acq: str) -> str:
    date, fly, thor = acq.split('/')
    return f'{date}__fly{int(fly):02d}__{thor}'


def megamat_src(acq: str) -> Dict[str, Path]:
    c = REMY_MATRIX / acq / 'source_extraction_s2p' / 'suite2p' / 'combined'
    return dict(respvec=c / 'xrds_suite2p_respvec_mean_peak.nc', trials=c / 'iscell_soma' / 'xrds_suite2p_trials.nc')


def reduce_megamat_trials(trials_path: Path, respvec_path: Path,
                          window: Dict[str, Tuple[float, float]] = KC_WINDOW_REMY) -> Tuple[pd.DataFrame, dict]:
    """One 2022 megamat fly's KC somata from her suite2p trials file (Fc frames, -10 to 44.75 s at 0.25 s): per cell
    of the trials file whose cluster is one of her good clusters (respvec's `good_xid`, the cells
    `load_remy_megamat_kc_responses` keeps), per trial the response window mean minus the pre-odor window mean (Fc), and
    per cell the SD of its pre-odor frames (each trial's SD, the root mean square over trials, as her validation2
    baseline_std is pooled in `v2_soma_long`). Returns long rows (cell, trial, stim, dfc, z_respvec) and per-cell noise
    in `sd` (one value per row's cell), with a summary."""
    import xarray as xr
    rv = xr.open_dataset(respvec_path)
    good = set(int(x) for x in np.atleast_1d(rv.attrs['good_xid']))
    keep_rv = np.isin(rv['xid0'].values, list(good)) & (rv['iscell'].values == 1)
    good_cells = set(int(c) for c in rv['cells'].values[keep_rv])
    zrv = rv['Fc_zscore'].to_pandas()                           # trials x cells
    tr = xr.open_dataset(trials_path)
    t = tr['time'].values
    cells = tr['cells'].values
    sel = np.array([int(c) in good_cells for c in cells])
    Fc = tr['Fc'].values                                        # trials x cells x time
    Fc = Fc[:, sel]
    resp = (t >= window['response'][0]) & (t <= window['response'][1])
    pre = (t >= window['pre'][0]) & (t <= window['pre'][1])
    dfc = Fc[:, :, resp].mean(axis=2) - Fc[:, :, pre].mean(axis=2)
    sd = np.sqrt(np.mean(np.std(Fc[:, :, pre], axis=2) ** 2, axis=0))
    stims = [str(s) for s in tr['stim'].values]
    cs = cells[sel]
    z = zrv.reindex(columns=cs).to_numpy()
    rows = pd.DataFrame(dict(cell=np.tile(cs, len(stims)), trial=np.repeat(np.arange(len(stims)), len(cs)),
                             stim=np.repeat(stims, len(cs)), dfc=dfc.reshape(-1),
                             sd=np.tile(sd, len(stims)), z_respvec=z.reshape(-1)))
    summ = dict(n_cells_trials_file=int(len(cells)), n_good_respvec=int(len(good_cells)), n_kept=int(sel.sum()),
                n_frames_response=int(resp.sum()), n_frames_pre=int(pre.sum()), window=window,
                frame_dt=float(np.median(np.diff(t))))
    rv.close()
    tr.close()
    return rows, summ


def megamat_soma_long(copy_dir: Optional[Path] = None) -> pd.DataFrame:
    """Remy's four 2022 megamat flies' KC somata in baseline-noise units (sources/remy_megamat_kc_soma/
    frames_reduced.parquet): long rows fly ('date/fly_num'), unit (cell), odor, repeat, value."""
    d = megamat_copy_dir() if copy_dir is None else Path(copy_dir)
    r = pd.read_parquet(d / 'frames_reduced.parquet')
    r = r[np.isfinite(r.sd) & (r.sd > 0)]
    odor = r.stim.str.split(' @ ').str[0]
    rep = r.groupby(['fly', 'cell', 'stim']).cumcount()
    out = pd.DataFrame(dict(fly=r.fly.to_numpy(), unit=r.cell.astype(str).to_numpy(), odor=odor.to_numpy(),
                            repeat=rep.to_numpy(), value=(r.dfc / r.sd).to_numpy()))
    return out[out.odor.isin(MEGAMAT_ODORS)].reset_index(drop=True)


def z_bridge(copy_dir: Optional[Path] = None, criterion: float = CRITERIA['kc_soma']) -> dict:
    """The session-z criterion equivalent to the KC somata's baseline-noise criterion, on Remy's four megamat flies:
    the z at which the fraction of (cell, odor) trial means of her respvec Fc_zscore above it equals the fraction of
    baseline-noise trial means above `criterion`, pooled over the flies' cells and odors (each fly weighted by its
    cells). The criterion of the mixtures' KC somata, which olfdb holds in session z only."""
    d = megamat_copy_dir() if copy_dir is None else Path(copy_dir)
    r = pd.read_parquet(d / 'frames_reduced.parquet')
    r = r[np.isfinite(r.sd) & (r.sd > 0) & r.stim.str.split(' @ ').str[0].isin(MEGAMAT_ODORS)]
    m = r.assign(v=r.dfc / r.sd).groupby(['fly', 'cell', 'stim'])[['v', 'z_respvec']].mean()
    rate = float((m.v > criterion).mean())
    z = float(np.quantile(m.z_respvec.dropna(), 1 - rate))
    return dict(rate_noise=rate, z=z, criterion=criterion, n=int(len(m)),
                per_fly_rate=m.groupby(level='fly').v.apply(lambda x: float((x > criterion).mean())).to_dict())


@lru_cache(maxsize=None)
def mixture_z_criterion() -> float:
    """The z bridge recorded in sources/remy_megamat_kc_soma/z_bridge.json (`z_bridge`), rounded to 4 decimals."""
    b = _read_json(megamat_copy_dir() / 'z_bridge.json')
    if b is None:
        raise FileNotFoundError('sources/remy_megamat_kc_soma/z_bridge.json: run --copy first')
    return round(float(b['z']), 4)


# ---------------------------------------------------------------------------
# --copy: what the share lacks

def copy_sources(sha: bool = True) -> dict:
    """Copy to <shared>/data/stage_populations/sources/ what E134 reads and the share lacks, and reduce the 2022 megamat
    flies' trial frames: olfdb datasets 1, 4, 5 and 8 (long downloads, with olfdb's metadata), the GH146 and pebbled
    tables of al_analysis's sent_to_remy/2023-10-29, Remy's validation2 respvec and stim_rdm files (data_by_acq_2, 20
    acquisitions), the 2022 megamat flies' respvec files (matrix) and frames_reduced.parquet from their trials files
    (which stay where they are; their sha256 are recorded); then z_bridge.json. Writes files.csv (every source with
    bytes and sha256) and README.md."""
    SOURCES_DIR.mkdir(parents=True, exist_ok=True)
    rec = []
    # olfdb
    od = SOURCES_DIR / 'olfdb'
    od.mkdir(exist_ok=True)
    meta = {x['id']: x for x in olfdb_datasets()}
    for ds in sorted(set(OLFDB_DATASETS.values())):
        df = olfdb_download(ds)
        df.to_csv(olfdb_path(ds), index=False)
        rec.append(dict(kind='olfdb', source=f'{OLFDB}/api/download?ids={ds}&format=long', copy=str(olfdb_path(ds)),
                        olfdb_data_file=meta[ds]['data_file'], olfdb_data_sha256=meta[ds]['data_sha256']))
    _write_json(od / 'datasets.json', {ds: meta[ds] for ds in sorted(set(OLFDB_DATASETS.values()))})
    # GH146 / pebbled
    gd = SOURCES_DIR / 'gh146'
    gd.mkdir(exist_ok=True)
    for kind in ('GH146', 'pebbled'):
        for ext in ('p', 'csv'):
            src = gh146_dir() / f'{kind}_ij_certain-roi_stats.{ext}'
            shutil.copy2(src, gd / src.name)
            rec.append(dict(kind='al_analysis', source=str(src), copy=str(gd / src.name)))
    # validation2 respvec (tensor)
    vd = v2_copy_dir()
    for acq in v2_acqs(REMY_V2_DIR):
        for rel in V2_FILES:
            src = REMY_V2_DIR / acq / rel
            dst = vd / acq / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            rec.append(dict(kind='remy_validation2', source=str(src), copy=str(dst)))
    agg = REMY_V2_DIR / 'aggregate' / 'validation2_ijk_all_2' / 'df_acqs.csv'
    (vd / 'aggregate').mkdir(parents=True, exist_ok=True)
    shutil.copy2(agg, vd / 'aggregate' / 'df_acqs.csv')
    rec.append(dict(kind='remy_validation2', source=str(agg), copy=str(vd / 'aggregate' / 'df_acqs.csv')))
    # the 2022 megamat flies (matrix)
    md = megamat_copy_dir()
    md.mkdir(parents=True, exist_ok=True)
    red, summ = [], {}
    for acq in REMY_MEGAMAT_ACQS:
        s = megamat_src(acq)
        lab = megamat_acq_label(acq)
        dst = md / lab / 'xrds_suite2p_respvec_mean_peak.nc'
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s['respvec'], dst)
        rec.append(dict(kind='remy_megamat', source=str(s['respvec']), copy=str(dst)))
        rows, sm = reduce_megamat_trials(s['trials'], s['respvec'])
        date, fly, _ = acq.split('/')
        rows.insert(0, 'fly', f'{date}/{int(fly)}')
        red.append(rows)
        summ[acq] = sm
        rec.append(dict(kind='remy_megamat_frames', source=str(s['trials']), copy=None))
    pd.concat(red, ignore_index=True).to_parquet(md / 'frames_reduced.parquet', index=False)
    _write_json(md / 'frames_reduced_summary.json', summ)
    _write_json(md / 'z_bridge.json', z_bridge(md))
    mixture_z_criterion.cache_clear()
    # sha256
    for r in rec:
        r['source_entry'] = file_entry(Path(r['source']), sha=sha) if not r['source'].startswith('http') else None
        r['copy_entry'] = file_entry(Path(r['copy']), sha=sha) if r.get('copy') else None
    flat = [dict(kind=r['kind'], source=r['source'], copy=r.get('copy'),
                 source_bytes=(r['source_entry'] or {}).get('bytes'), source_sha256=(r['source_entry'] or {}).get('sha256'),
                 copy_bytes=(r['copy_entry'] or {}).get('bytes'), copy_sha256=(r['copy_entry'] or {}).get('sha256'),
                 olfdb_data_file=r.get('olfdb_data_file'), olfdb_data_sha256=r.get('olfdb_data_sha256')) for r in rec]
    pd.DataFrame(flat).to_csv(SOURCES_DIR / 'files.csv', index=False)
    write_sources_readme()
    return dict(n_files=len(rec), megamat=summ)


def write_sources_readme() -> None:
    (SOURCES_DIR / 'README.md').write_text(f"""# E134's sources: what the share lacked, copied and reduced

Made by `data_prep/stage_population_data.py --copy` (model_experiment_queue) on {_host()}, {_now()}. files.csv lists
every source with its path, bytes and sha256, and its copy here. Read only from the sources; nothing written there.

- `olfdb/ds<id>_long.csv`: olfdb's long downloads (http://hal:8080/api/download?ids=<id>&format=long) of datasets 1
  (megamat ORNs, 9 flies), 4 (validation2 ORNs, 6 flies), 5 (kiwi / control mixture ORNs, 9 flies) and 8 (Remy's kiwi /
  control mixture KCs, 11 flies, trial means of Fc_zscore); datasets.json has olfdb's metadata (its data_sha256 is of
  the uploaded file, not of the download).
- `gh146/`: al_analysis's data/sent_to_remy/2023-10-29 tables, GH146 PN dendrites (7 flies) and pebbled ORNs, per trial
  dF/F (E90's `al_table` / `al_trials` read them).
- `remy_validation2_kc_soma/<acq>/`: Remy's validation2 KC somata, the 20 acquisitions of her aggregate
  validation2_ijk_all_2 (aggregate/df_acqs.csv), from /local/tensor/Remy/odor_space_collab/analysis_outputs/
  data_by_acq_2/kc_soma_nls/: respvec/xrds_peak_amp.nc (per trial and cell, the 2 to 8 s mean minus the -10 to -0.5 s
  mean), respvec/xrds_baseline_std.nc (the SD of each trial's -10 to -0.5 s frames) and stim_rdm/xrds_stim_rdm.nc
  (her per-fly distances, for the check). The 660 MB trial files stay where they are.
- `remy_megamat_kc_soma/`: Remy's four 2022 megamat flies (2022-10-10 flies 1 and 2, 2022-10-11 fly 1, 2022-11-10 fly
  1), from her processed folders on matrix (/local/matrix/Remy-Data/projects/odor_space_collab/processed_data/<date>/
  <fly>/<acq>/source_extraction_s2p/suite2p/combined/): `<acq>/xrds_suite2p_respvec_mean_peak.nc` (the file
  pyolfsysm's `load_remy_megamat_kc_responses` reads, which the olfsysm checkout lacks), and `frames_reduced.parquet`,
  reduced from her iscell_soma/xrds_suite2p_trials.nc (0.85 to 1.55 GB each; not copied, their sha256 in files.csv):
  per good cell (her `good_xid` clusters) and trial, the 2 to 8 s mean minus the -10 to -0.5 s mean of Fc (`dfc`), the
  cell's pre-odor SD (`sd`, the root mean square over trials of each trial's SD) and her respvec Fc_zscore
  (`z_respvec`). `z_bridge.json`: the session-z criterion of the mixtures' KC somata (the module's `z_bridge`).
""")


# ---------------------------------------------------------------------------
# --build: the population files

def population_spec(key: Tuple[str, str, str]) -> dict:
    """How a population is read: units, criterion, window, notes and source paths (for the lock)."""
    stage, panel, cohort = key
    notes = []
    if panel == 'apl16':
        notes.append('paa, euc and -aPine: the model input for these is ORN data within 1 log10 unit of the flies\' '
                     'concentrations (paa 10^-3 in the flies, 10^-4 in the input; -aPine and euc 10^-2 against 10^-1.5 '
                     'and 10^-1.75)')
    if panel == 'validation2' and stage == 'kc_soma':
        notes.append("euc, -aPine and ger: each fly's one concentration pooled under one name (E107's POOLED; "
                     "euc -2 / -1.75, -aPine -2 / -1.5, ger -2 / -3); geos at 10^-4 only (10^-5 left out); "
                     'every fly was given part of the panel')
    if panel == 'validation2' and stage == 'orn':
        notes.append('the ORN recordings at -aPine 10^-1.5, euc 10^-1.75; the KC flies mostly at 10^-2 (within 1 log10 '
                     'unit)')
    if panel in MIXTURE_PANELS and stage == 'kc_soma':
        notes.append('session z (Fc_zscore) trial means from olfdb dataset 8; no frames or single trials on hand, so not '
                     'in baseline-noise units; the criterion is the z bridge on the 2022 megamat flies')
    if panel in MIXTURE_PANELS and stage == 'orn':
        notes.append('the mixture ORN flies (2024) are not the KC flies (2022); measured mixtures used as measured')
    if stage in ('orn', 'pn_dendrite'):
        notes.append('dF/F per trial from processed tables (no frames on hand): trial-noise units')
    if cohort == 'same_fly':
        notes.append('one fly (2023-04-11 fly 3, megamat1_redo), claws and somata on the same trials (E119 part 0b)')
    if key == ('kc_soma', 'megamat', 'wild_type'):
        notes.append("Remy's four 2022 flies with the whole panel; her 14 older flies (2018-2019, partial panels, "
                     'another pipeline) are not used')
    units = UNITS_TRIAL if stage in ('orn', 'pn_dendrite') else (
        UNITS_SESSION_Z if (stage == 'kc_soma' and panel in MIXTURE_PANELS) else UNITS_BASELINE)
    if stage in ('pn_bouton', 'kc_claw') and cohort == 'wild_type':
        window = {'response': (0.5, 3.0), 'pre': (-5.0, -0.2)}
    elif cohort == 'same_fly':
        z = np.load(roi_dir() / 'same_fly_trials.npz', allow_pickle=True)
        window = {'response': tuple(float(x) for x in z['t_window']), 'pre': tuple(float(x) for x in z['pre'])}
    elif panel == 'apl16':
        window = KC_WINDOW_YANG
    elif stage == 'kc_soma' and panel in ('megamat', 'validation2'):
        window = KC_WINDOW_REMY
    else:
        window = {}
    return dict(units=units, criterion=response_criterion(stage, panel, cohort), window=window, notes=notes,
                recording_site=RECORDING_SITE[stage])


def source_paths(key: Tuple[str, str, str]) -> List[Path]:
    """The files `build_population` reads for `key`."""
    stage, panel, cohort = key
    if stage == 'orn':
        ds = {'megamat': 1, 'validation2': 4, 'kmix': 5, 'cmix': 5}[panel]
        return [olfdb_path(ds)]
    if stage == 'pn_dendrite':
        return [gh146_copy('GH146')]
    if cohort == 'same_fly':
        return [roi_dir() / 'same_fly_trials.npz']
    if stage in ('pn_bouton', 'kc_claw'):
        return [roi_dir() / 'window_means.parquet', roi_dir() / 'noise_sd.parquet', roi_dir() / 'lock.json']
    if panel == 'apl16':
        return [apl16_dir() / 'per_trial_long.parquet']
    if panel in MIXTURE_PANELS:
        return [olfdb_path(8)]
    if panel == 'validation2':
        d = v2_copy_dir()
        return [d / 'aggregate' / 'df_acqs.csv'] + [d / a / r for a in v2_acqs() for r in V2_FILES[:2]] + [
            SHARED / 'data' / 'remy_validation2_kc' / 'validation2' / 'kc_soma' / 'stim_to_drop_by_qcthresh.csv']
    if panel == 'megamat':
        return [megamat_copy_dir() / 'frames_reduced.parquet']
    raise KeyError(key)


def build_long(key: Tuple[str, str, str]) -> Tuple[pd.DataFrame, List[dict]]:
    """The population's long rows (fly, unit, odor, repeat, value) and its exclusions."""
    stage, panel, cohort = key
    excl = [dict(x) for x in EXCLUDED_FLIES if x['stage'] == stage and x['panel'] == panel and cohort == 'wild_type']
    if stage == 'orn':
        ds = {'megamat': 1, 'validation2': 4, 'kmix': 5, 'cmix': 5}[panel]
        d = orn_long(ds, panel)
        return d.drop(columns='noise_sd'), excl
    if stage == 'pn_dendrite':
        return pn_dendrite_long().drop(columns='noise_sd'), excl
    if cohort == 'same_fly':
        return same_fly_long(stage), excl
    if stage in ('pn_bouton', 'kc_claw'):
        return per_roi_long(stage), excl
    if panel == 'apl16':
        return apl16_long(cohort), excl
    if panel in MIXTURE_PANELS:
        return mixture_kc_long(panel), excl
    if panel == 'validation2':
        d, e = v2_soma_long()
        return d, excl + e
    if panel == 'megamat':
        excl.append(dict(what='flies', reason="Remy's 14 older megamat flies (2018-2019): partial panels and another "
                         'pipeline (peak_amp of 2024-11-12); not read'))
        return megamat_soma_long(), excl
    raise KeyError(key)


def build_populations(keys: Optional[Sequence[Tuple[str, str, str]]] = None, pop_dir: Optional[Path] = None) -> dict:
    """Write every population's long parquet (`pop_path`) and a summary (populations/summary.json)."""
    keys = POPULATIONS if keys is None else keys
    pd_ = POP_DIR if pop_dir is None else Path(pop_dir)
    pd_.mkdir(parents=True, exist_ok=True)
    summ = _read_json(pd_ / 'summary.json') or {}
    for key in keys:
        t0 = time.time()
        d, excl = build_long(key)
        d = d.astype({'fly': str, 'unit': str, 'odor': str, 'repeat': int, 'value': float})
        d.to_parquet(pop_path(key, pd_), index=False)
        summ[pop_name(key)] = dict(rows=int(len(d)), flies=sorted(d.fly.unique()),
                                   units_per_fly=d.groupby('fly').unit.nunique().to_dict(),
                                   odors_per_fly=d.groupby('fly').odor.nunique().to_dict(), exclusions=excl,
                                   seconds=round(time.time() - t0, 1))
        print(f'{pop_name(key)}: {len(d)} rows, {d.fly.nunique()} flies ({time.time() - t0:.0f} s)', flush=True)
    _write_json(pd_ / 'summary.json', summ)
    if pop_dir is None:
        write_populations_readme()
    return summ


def write_populations_readme() -> None:
    """<shared>/data/stage_populations/README.md."""
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / 'README.md').write_text(f"""# E134's recorded populations, every stage and panel

Made by `data_prep/stage_population_data.py` (model_experiment_queue) on {_host()}, {_now()}; read them with its
`load_population(stage, panel, cohort)`, which checks each file's sha256 against the lock
(results/E134_population_baseline/lock.json, committed as experiments/E134_lock.json).

- `sources/`: what the share lacked, copied and reduced by `--copy` (README.md and files.csv there).
- `populations/<stage>__<panel>__<cohort>.parquet`: one per population (15), long rows fly, unit, odor, repeat, value
  (repeat -1: a trial mean, where the source has no single trials). Units: baseline-noise SDs where the recordings'
  frames exist (boutons, claws, KC somata of megamat, the same fly, validation2 and APL16), trial-noise SDs for the
  glomeruli (ORNs, PN dendrites; dF/F per trial only), Remy's session z for the mixtures' KC somata (olfdb dataset 8).
  The module's docstring defines each.
- `populations/summary.json`: rows, flies, units and odors per fly, exclusions; `populations/checks.json`: the checks
  against the existing loaders (the pebbled table, E119's per-ROI reading, the bands' APL16 table, pyolfsysm's
  `load_remy_megamat_kc_responses`, Remy's per-fly validation2 matrices) and the noise bridge.

Excluded everywhere: the bouton fly 2022-10-21 fly02. Panels: megamat, validation2, kmix and cmix (the kiwi and the
control mixtures; each KC fly had one of the two), apl16.
""")


# ---------------------------------------------------------------------------
# loading

_SHA_CACHE: Dict[str, str] = {}


def _cached_sha(path: Path) -> str:
    p = str(path)
    if p not in _SHA_CACHE:
        _SHA_CACHE[p] = sha256_file(path)
    return _SHA_CACHE[p]


def population_from_long(key: Tuple[str, str, str], d: pd.DataFrame, *, trials: bool = False,
                         spec: Optional[dict] = None, exclusions: Optional[List[dict]] = None,
                         files: Optional[List[dict]] = None) -> Population:
    """A `Population` from long rows (fly, unit, odor, repeat, value): per fly, units x odors trial means (columns the
    panel's odors), and the single trials if `trials` and the rows have them (repeat >= 0)."""
    stage, panel, cohort = key
    odors = PANEL_ODORS[panel]
    spec = population_spec(key) if spec is None else spec
    tables, trs = {}, {}
    flies = sorted(d.fly.unique())
    for fly in flies:
        sub = d[d.fly == fly]
        units = list(pd.unique(sub.unit))
        m = sub.groupby(['unit', 'odor'], sort=False).value.mean().unstack('odor')
        tables[fly] = m.reindex(index=units, columns=list(odors))
        tables[fly].columns.name = 'odor'
        tables[fly].index.name = 'unit'
        if trials and (sub.repeat >= 0).any():
            t = sub[sub.repeat >= 0].pivot_table(index='unit', columns=['odor', 'repeat'], values='value',
                                                  aggfunc='first')
            trs[fly] = t.reindex(index=units)
    return Population(stage=stage, panel=panel, cohort=cohort, odors=tuple(odors), flies=tuple(flies), tables=tables,
                      criterion=float(spec['criterion']), units=spec['units'], graded=True,
                      trials=(trs if trials and trs else None), has_trials=bool((d.repeat >= 0).any()),
                      unit_kind=UNIT_KIND[stage],
                      named_units=NAMED_UNITS[stage], recording_site=spec['recording_site'],
                      window={k: tuple(v) for k, v in spec['window'].items()}, exclusions=list(exclusions or []),
                      files=list(files or []), notes=list(spec['notes']))


def load_population(stage: str, panel: str, cohort: str = 'wild_type', *, trials: bool = False,
                    check_lock: bool = True, pop_dir: Optional[Path] = None) -> Population:
    """One recorded population (a key of `POPULATIONS`), excluded flies left out.

    Args:
        stage, panel, cohort: the key ('mixtures' is not a panel: use 'kmix' or 'cmix').
        trials: also load the single-trial tables (`Population.trials`) where they exist.
        check_lock: once the lock exists, refuse (ValueError) if the population file differs from its sha256 in the
            lock; the criterion, units and window are then the lock's.
        pop_dir: where the population files are (default <shared>/data/stage_populations/populations/).
    """
    key = (stage, panel, cohort)
    if key not in POPULATIONS:
        raise KeyError(f'{key} is not one of POPULATIONS')
    path = pop_path(key, pop_dir)
    if not path.exists():
        raise FileNotFoundError(f'{path}: run --build first')
    lock = read_lock() if check_lock else None
    spec, excl, files = None, None, None
    if lock is not None and pop_name(key) in lock.get('populations', {}):
        L = lock['populations'][pop_name(key)]
        sha = _cached_sha(path)
        if sha != L['file']['sha256']:
            raise ValueError(f"{path.name}: sha256 {sha[:12]} differs from the lock's {L['file']['sha256'][:12]}")
        spec = dict(units=L['units'], criterion=L['criterion'], window=L['window'], notes=L['notes'],
                    recording_site=L['recording_site'])
        excl = L['exclusions']
        files = [L['file']] + L['sources']
    d = pd.read_parquet(path)
    return population_from_long(key, d, trials=trials, spec=spec, exclusions=excl, files=files)


def load_populations(keys: Optional[Sequence[Tuple[str, str, str]]] = None, **kwargs) -> Dict[Tuple[str, str, str],
                                                                                                Population]:
    """`load_population` for every key (all of `POPULATIONS` by default); a population whose file is missing is left
    out, and listed by `missing_populations`."""
    out = {}
    for key in (POPULATIONS if keys is None else keys):
        try:
            out[key] = load_population(*key, **kwargs)
        except FileNotFoundError:
            continue
    return out


def missing_populations(pop_dir: Optional[Path] = None) -> Dict[Tuple[str, str, str], str]:
    """Keys whose population file or sources are not on hand, each with what would be needed (outcome (f))."""
    out = {}
    for key in POPULATIONS:
        if not pop_path(key, pop_dir).exists():
            miss = [str(p) for p in source_paths(key) if not Path(p).exists()]
            out[key] = ('sources missing: ' + ', '.join(miss) + ' (run --copy)') if miss else 'run --build'
    return out


# ---------------------------------------------------------------------------
# checks against the existing loaders

def check_orn_megamat_vs_pebbled(tol: float = 1e-6) -> dict:
    """olfdb dataset 1's megamat trials equal al_analysis's pebbled table (E90's ORN reading) on the flies, glomeruli,
    odors and repeats both hold."""
    a = pd.read_csv(olfdb_path(1))
    a = a[a.panel == 'megamat']
    a = a.assign(fly=[_fly_label(x, y) for x, y in zip(a.date, a.fly_num)], odor=a.odor1.str.split(' @ ').str[0],
                 unit=a.roi.astype(str))
    b = al_trials_long('pebbled')
    m = a.merge(b, on=['fly', 'unit', 'odor', 'repeat'], suffixes=('_olfdb', '_al'))
    diff = float(np.max(np.abs(m.value_olfdb - m.value_al))) if len(m) else float('nan')
    return dict(n_matched=int(len(m)), n_olfdb=int(len(a)), n_al=int(len(b)), max_abs_diff=diff,
                passed=bool(len(m) > 0 and diff <= tol))


def check_per_roi_vs_e119() -> dict:
    """The boutons' and claws' trial means equal E119's (`fly_trials`) noise reading, by construction; checked on the
    built files."""
    from experiments import E119_stages_per_roi as E119
    out = {}
    for stage, st in (('pn_bouton', 'bouton'), ('kc_claw', 'claw')):
        p = load_population(stage, 'megamat', check_lock=False)
        ft = E119.fly_trials(st, odors=MEGAMAT_ODORS)
        dmax = 0.0
        for fly in p.flies:
            ref = np.nanmean(ft[fly]['noise'], axis=0)
            ours = p.tables[fly].reindex(index=ft[fly]['rois'].astype(str)).to_numpy()
            dmax = max(dmax, float(np.nanmax(np.abs(ours - ref))))
        out[stage] = dmax
    return dict(max_abs_diff=out, passed=bool(all(v < 1e-9 for v in out.values())))


def check_apl16_vs_trialmean(tol: float = 1e-6) -> dict:
    """The APL16 populations' trial means equal trialmean_long.parquet (the realism bands' table)."""
    t = pd.read_parquet(apl16_dir() / 'trialmean_long.parquet')
    out = {}
    for cohort, g in (('control', 'TNT3c'), ('tnt', 'TNT')):
        p = load_population('kc_soma', 'apl16', cohort, check_lock=False)
        s = t[(t.genotype == g) & t.odor.isin(APL16_ODORS)]
        fly = s.date.astype(str) + '/' + s.fly_num.astype(int).astype(str) + '/' + s.hemibrain.astype(str)
        s = s.assign(fly=fly, unit=s.roi.astype(int).astype(str))
        dmax = 0.0
        for f, sub in s.groupby('fly'):
            ref = sub.pivot(index='unit', columns='odor', values='value')
            ours = p.tables[f].reindex(index=ref.index, columns=ref.columns)
            dmax = max(dmax, float(np.nanmax(np.abs(ours.to_numpy() - ref.to_numpy()))))
        out[cohort] = dict(max_abs_diff=dmax, n_flies=len(p.flies))
    return dict(per_cohort=out, passed=bool(all(v['max_abs_diff'] <= tol for v in out.values())))


def check_megamat_vs_loader(tol: float = 1e-9) -> dict:
    """The 2022 megamat flies' cells are those pyolfsysm's `load_remy_megamat_kc_responses` keeps (iscell, cluster in
    `good_xid`) that the trials file holds, and her respvec trial means of Fc_zscore over those cells are what the
    loader returns (the loader is re-run on a folder holding the copies, the layout it expects)."""
    import tempfile
    import pyolfsysm.mb_model as mb
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / 'from_remy'
        for acq in REMY_MEGAMAT_ACQS:
            lab = megamat_acq_label(acq)
            dst = root / 'megamat17' / 'per_fly' / lab / 'respvec'
            dst.mkdir(parents=True)
            (dst / 'xrds_suite2p_respvec_mean_peak.nc').symlink_to(megamat_copy_dir() / lab /
                                                                    'xrds_suite2p_respvec_mean_peak.nc')
        (root / '2024-11-12').mkdir(parents=True)
        old = mb.remy_data_dir
        try:
            mb.remy_data_dir = root
            ref = mb.load_remy_megamat_kc_responses()
        except Exception as e:                                  # the loader's other checks, reported
            return dict(passed=False, error=repr(e))
        finally:
            mb.remy_data_dir = old
    r = pd.read_parquet(megamat_copy_dir() / 'frames_reduced.parquet')
    ours = r.groupby(['fly', 'cell', 'stim']).z_respvec.mean().unstack('stim')
    ours.columns = [c.split(' @ ')[0] for c in ours.columns]
    lv = list(ref.index.names)
    ref = ref.copy()
    ref.columns = [str(c).split(' @ ')[0] for c in ref.columns]
    dmax, n, n_cells_differ = 0.0, 0, 0
    for fly, sub in ours.groupby(level='fly'):
        # the loader numbers a fly's cells 0..n-1 in its respvec order (the good cells by suite2p number)
        sub = sub.sort_index(level='cell')
        r = ref[ref.index.get_level_values('datefly') == fly]
        if len(r) != len(sub):
            n_cells_differ += 1
            continue
        A = r.reindex(columns=sub.columns).to_numpy(dtype=float)
        dmax = max(dmax, float(np.max(np.abs(A - sub.to_numpy(dtype=float)))))
        n += len(sub)
    return dict(index_levels=lv, n_ref_cells=int(len(ref)), n_ours=int(len(ours)), n_matched=n,
                n_flies_cell_count_differs=n_cells_differ, max_abs_diff=dmax,
                passed=bool(n == len(ours) == len(ref) and dmax <= tol))


def noise_bridge() -> dict:
    """On every population with single trials in baseline-noise units, the median over units of the trial-noise SD
    (`robust_trial_sd`, in baseline-noise units): how the glomeruli's trial-noise unit compares with the baseline-noise
    unit of the stages that have frames (1 if the trial-to-trial noise were the baseline noise)."""
    out = {}
    for key in POPULATIONS:
        spec = population_spec(key)
        if spec['units'] != UNITS_BASELINE or not pop_path(key).exists():
            continue
        p = load_population(*key, trials=True, check_lock=False)
        if p.trials is None:
            continue
        meds = []
        for fly in p.flies:
            t = p.trials[fly]
            reps = sorted(set(t.columns.get_level_values('repeat')))
            X = np.stack([t.xs(r, axis=1, level='repeat').reindex(columns=list(p.odors)).to_numpy() for r in reps])
            meds.append(float(np.nanmedian(robust_trial_sd(X))))
        out[pop_name(key)] = dict(median_over_flies=float(np.median(meds)), per_fly=meds)
    return out


def run_checks() -> dict:
    """Every check; written to populations/checks.json."""
    checks = {}
    for name, fn in (('orn_megamat_vs_pebbled', check_orn_megamat_vs_pebbled), ('per_roi_vs_e119', check_per_roi_vs_e119),
                     ('apl16_vs_trialmean', check_apl16_vs_trialmean), ('megamat_vs_loader', check_megamat_vs_loader),
                     ('v2_rdms_cosine', lambda: check_v2_rdms(metric='cosine')),
                     ('v2_rdms_correlation', lambda: check_v2_rdms(metric='correlation'))):
        try:
            checks[name] = fn()
        except Exception as e:
            checks[name] = dict(passed=False, error=repr(e))
        print(name, checks[name].get('passed'), flush=True)
    try:
        checks['noise_bridge'] = noise_bridge()
    except Exception as e:
        checks['noise_bridge'] = dict(error=repr(e))
    checks['passed'] = bool(all(v.get('passed', True) for k, v in checks.items() if isinstance(v, dict)
                                and k != 'noise_bridge'))
    _write_json(POP_DIR / 'checks.json', checks)
    return checks


# ---------------------------------------------------------------------------
# the lock

def _commits() -> dict:
    try:
        from kcq import harness as H
        return H.code_commits()
    except Exception as e:                                      # provenance only
        return dict(error=repr(e))


def response_rates(pop: Population) -> dict:
    """The population's overall response rate (the fraction of units and odors above the criterion, over the odors
    each fly was given), per fly and as the mean over flies (each fly weighted equally), and per odor (mean over flies)."""
    per_fly, per_odor = {}, []
    for f in pop.flies:
        r = pop.responses(f)
        per_fly[f] = float(np.nanmean(r.to_numpy(dtype=float)))
        per_odor.append(r.mean(axis=0))
    po = pd.concat(per_odor, axis=1).mean(axis=1)
    return dict(mean_over_flies=float(np.mean(list(per_fly.values()))), per_fly=per_fly,
                per_odor={o: (None if not np.isfinite(v) else float(v)) for o, v in po.items()})


def remy_criterion_rate() -> Optional[dict]:
    """Remy's four 2022 megamat flies' overall response rate by her own criterion (pyolfsysm's
    `load_remy_megamat_kc_binary_responses`, rates within each fly then averaged, as `remy_megamat_kc_targets`): the
    operating point the main goal names (0.0915)."""
    try:
        import pyolfsysm.mb_model as mb
        b = mb.load_remy_megamat_kc_binary_responses()
        rates = mb.mean_response_rates(b, fly_level=mb.remy_fly_id)
        per_fly = b.astype(float).groupby(level=mb.remy_fly_id, sort=False).mean().mean(axis=1)
        return dict(mean=float(rates.mean()), per_fly={str(k): float(v) for k, v in per_fly.items()},
                    source="pyolfsysm load_remy_megamat_kc_binary_responses (Remy's df_stim_responders_from_peak.csv)")
    except Exception as e:                                      # reported, not fatal
        return dict(error=repr(e))


def model_readings() -> Optional[dict]:
    """atlas's `kcq.stage_readings.reading_definitions()` (each stage's default reading of the model and its
    alternatives), recorded in the lock so that no reading changes after it."""
    try:
        from kcq import stage_readings as R
        return R.reading_definitions()
    except Exception as e:                                      # reported, not fatal
        return dict(error=repr(e))


def lock_contents(sha_sources: bool = True) -> dict:
    """The lock as a dict: per population the file (path, bytes, sha256), its sources (likewise), flies, odors (the
    panel's, and per fly), units per fly, exclusions, criterion, units, window, breadth bins (of the panel's size, and
    the sizes of partial flies), recording site and notes; the copy sources (files.csv), the checks and the commits."""
    summ = _read_json(POP_DIR / 'summary.json') or {}
    checks = _read_json(POP_DIR / 'checks.json')
    pops = {}
    for key in POPULATIONS:
        path = pop_path(key)
        if not path.exists():
            continue
        spec = population_spec(key)
        s = summ.get(pop_name(key), {})
        odors = PANEL_ODORS[key[1]]
        sizes = sorted(set(s.get('odors_per_fly', {}).values()) | {len(odors)})
        pop = load_population(*key, check_lock=False)
        pops[pop_name(key)] = dict(
            key=list(key), file=file_entry(path), sources=[file_entry(p, sha=sha_sources) for p in source_paths(key)],
            has_trials=pop.has_trials, response_rate=response_rates(pop),
            flies=s.get('flies'), n_flies=len(s.get('flies') or []), odors=list(odors),
            odors_per_fly=s.get('odors_per_fly'), units_per_fly=s.get('units_per_fly'),
            exclusions=s.get('exclusions', []), criterion=spec['criterion'], units=spec['units'],
            window={k: list(v) for k, v in spec['window'].items()}, unit_kind=UNIT_KIND[key[0]],
            named_units=NAMED_UNITS[key[0]], recording_site=spec['recording_site'], notes=spec['notes'],
            breadth_bins={str(n): [list(b) for b in breadth_bins(n)] for n in sizes if n >= 4})
    files_csv = SOURCES_DIR / 'files.csv'
    return dict(
        experiment='E134_population_baseline', step='0, the lock', host=_host(), time=_now(), commits=_commits(),
        populations=pops,
        stages=dict(criteria=CRITERIA, mixture_z_criterion=mixture_z_criterion() if
                    (megamat_copy_dir() / 'z_bridge.json').exists() else None,
                    units={UNITS_BASELINE: 'window mean minus pre-odor mean over the SD of the pre-odor frames',
                           UNITS_TRIAL: 'trial mean over the robust SD of single trials about their odor means '
                                        '(1.4826 MAD, residuals scaled by sqrt(n/(n-1)))',
                           UNITS_SESSION_Z: "Remy's whole-session z of Fc"},
                    orn_recordings='the ORN recordings are of the ORNs\' axon terminals in the antennal lobe '
                                   '(pebbled > GCaMP6f), not of their somata or dendrites',
                    recording_site=RECORDING_SITE),
        exclusions=dict(flies=list(EXCLUDED_FLIES),
                        rules=['the bouton fly 2022-10-21 fly02 in every analysis',
                               "Remy's validation2 quality-control drops at " + json.dumps(V2_QC),
                               'validation2 stimuli ' + ', '.join(V2_EXCLUDED_STIMULI),
                               'the solvent (pfo) in every panel',
                               "the mixture ORN panels' binary mixtures (odor2 not 'solvent') and pfo",
                               'units with a zero or undefined noise SD']),
        sources_files_csv=file_entry(files_csv) if files_csv.exists() else None,
        sources_note='files.csv lists every copied and reduced source (path, bytes, sha256), the large trial files on '
                     'tensor and matrix included; paths in this lock are relative to the shared directory (<shared>) '
                     'or the repository (<queue>) where they are under one',
        checks=checks,
        readings=model_readings(),
        remy_megamat_criterion_rate=remy_criterion_rate(),
        note='Written by data_prep/stage_population_data.py --lock before any model is read on these populations; '
             'load_population refuses a population file whose sha256 differs from this.')


def write_lock(path: Path = LOCK_PATH, copy: Path = LOCK_COPY) -> dict:
    """Write `lock_contents()` to the shared directory's lock and its committed copy in experiments/."""
    lock = lock_contents()
    _write_json(path, lock)
    _write_json(copy, lock)
    read_lock.cache_clear()
    return lock


@lru_cache(maxsize=None)
def read_lock(path: Path = LOCK_COPY) -> Optional[dict]:
    """The committed lock (experiments/E134_lock.json), or None before it exists."""
    return _read_json(path)


def _committed(path: Path) -> bool:
    repo = Path(path).resolve().parent
    r = subprocess.run(['git', '-C', str(repo), 'ls-files', '--error-unmatch', str(path)], capture_output=True)
    if r.returncode != 0:
        return False
    r = subprocess.run(['git', '-C', str(repo), 'diff', '--quiet', 'HEAD', '--', str(path)], capture_output=True)
    return r.returncode == 0


def require_lock() -> dict:
    """The committed lock; raises RuntimeError if it is missing or uncommitted (no model is read before it)."""
    lock = read_lock()
    if lock is None:
        raise RuntimeError(f'{LOCK_COPY}: no lock; run --lock and commit it first')
    if not _committed(LOCK_COPY):
        raise RuntimeError(f'{LOCK_COPY} is not committed (or has uncommitted changes)')
    return lock


def check_against_lock() -> dict:
    """Every population loads against the lock (its file's sha256)."""
    out = {}
    for key in POPULATIONS:
        try:
            p = load_population(*key)
            out[pop_name(key)] = dict(ok=True, n_flies=len(p.flies))
        except Exception as e:
            out[pop_name(key)] = dict(ok=False, error=repr(e))
    return out


def main(argv: Optional[Sequence[str]] = None) -> None:
    """--copy, --build, --lock, --check (see the module docstring)."""
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--copy', action='store_true')
    ap.add_argument('--no-sha', action='store_true', help='--copy without sha256 of the sources')
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--only', nargs='*', help='--build these populations only (stage__panel__cohort)')
    ap.add_argument('--checks', action='store_true', help='the checks against the existing loaders')
    ap.add_argument('--lock', action='store_true')
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args(argv)
    if a.copy:
        print(json.dumps(copy_sources(sha=not a.no_sha), indent=1, default=str)[:4000])
    if a.build:
        keys = None if not a.only else [tuple(x.split('__')) for x in a.only]
        build_populations(keys)
    if a.checks or a.build:
        print(json.dumps(run_checks(), indent=1, default=str)[:6000])
    if a.lock:
        L = write_lock()
        print(f"lock: {len(L['populations'])} populations -> {LOCK_PATH} and {LOCK_COPY}")
    if a.check:
        print(json.dumps(check_against_lock(), indent=1))


if __name__ == '__main__':
    main()
