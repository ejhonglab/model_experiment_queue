#!/usr/bin/env python3
"""Shared harness for the ORN->KC remapping / APL experiments (E1-E10).

Loads the data targets (Yang's APL16 KC data: TNT = APL output blocked, TNT3c = control;
per recording, from the shared directory), the model input (16 APL16 odors from the
megamat and validation2 ORN spike-delta estimates), and the reference pair distances the
Nov 3 report used, and computes the report's statistics for any pair of odor x odor
distance matrices (APL on vs off, model vs data), with bootstrap-over-recordings error
bars.

Conventions (see experiment_priorities.md in the shared directory):
- The data matrices are 1 - cosine similarity over all KCs of a recording (trial-mean
  responses), averaged over recordings. `pair_distances` computes the same for model
  outputs (and 1 - Pearson, for comparison).
- Pairs are indexed by (odor_a, odor_b) tuples in `APL16_ODORS` order (`PAIRS16`).
- Model ORN input for paa, euc and -aPine is at the lab's validation2 concentrations, not
  Yang's (`APPROX_CONC_ODORS`, `APPROX_CONC_NOTE`): say so wherever results involving
  those odors are reported.

Since 2026-11-12 the generally useful functions are in hong2p (hong2p.rdm, hong2p.stats,
hong2p.imaging, hong2p.viz, hong2p.util) and pyolfsysm (pyolfsysm.readout, kc_stats, effects);
their names here are re-exports, or thin wrappers with this harness's defaults (APL16 odor
order, abbreviated odor labels, Yang's 1 z criterion, `yang_apl_effect()` as the data's effect).
test/test_frozen_step3.py checks they give the outputs the code here gave before.

Run `python kcq/harness.py` to reproduce the report's data-vs-data statistics
(and print bootstrap CIs over recordings) as a self-test.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.isotonic import IsotonicRegression

# Every model run from these experiments must use a clean olfsysm build of its checkout's
# C++ sources (CLAUDE.md, "Model provenance"): importing the harness makes `fit_mb_model`
# refuse anything else, regardless of AL_ANALYSIS_ALLOW_DIRTY_OLFSYSM. `model_provenance()`
# is what the scripts put in their provenance.json. The model's Python code is pyolfsysm
# (in the olfsysm repository; al_analysis.mb_model until 2026-11-12, which still forwards
# there).
from pyolfsysm import mb_model as _mb_model
_mb_model.STRICT_CLEAN_OLFSYSM = True

# The generally useful functions of this harness moved to hong2p (the odor-pair matrices, their statistics,
# the imaging time course and noise, the colormap, the checkouts' git state) and pyolfsysm (the model read
# as imaging sees it, KC population statistics, the APL effect as the effect of a manipulation), 2026-11-12:
# the names below stay, as re-exports or as thin wrappers with this harness's defaults (the APL16 odor
# order, odor labels abbreviated, Yang's response criterion of 1 z, the KC data's APL effect).
from hong2p import imaging as _imaging
from hong2p import rdm as _rdm
from hong2p import stats as _stats
from hong2p.util import git_commits as _git_commits
from hong2p.viz import oklab_diverging_cmap as _oklab_diverging_cmap
from pyolfsysm import effects as _effects
from pyolfsysm import kc_stats as _kc_stats
from pyolfsysm import readout as _readout


# this repository's root (model_experiment_queue; this file is kcq/harness.py); until 2026-11-11 the repository was
# the folder scripts/claude_pair_remapping_experiments of al_analysis
HERE = Path(__file__).resolve().parents[1]

# the report's odor order (= Yang's original matrix files)
APL16_ODORS: List[str] = ['2h', 'IaA', 'pa', '2-but', 'eb', 'ep', '6al', 't2h', '1-8ol',
    '1-5ol', '1-6ol', 'ms', 'paa', '-aPine', 'euc', 'B-myr'
]
MEGAMAT12: List[str] = APL16_ODORS[:12]
VALIDATION2_4: List[str] = APL16_ODORS[12:]
SOLVENT = 'pfo'
# index level names of fit_mb_model outputs
KC_ID_LEVEL = 'kc_id'
KC_TYPE_LEVEL = 'kc_type'

# chemical groups used for the report's block analysis (null_models.py)
ODOR_GROUPS: Dict[str, str] = {'2h': 'fruity', 'IaA': 'fruity', 'pa': 'fruity',
    '2-but': 'small_ester', 'eb': 'small_ester', 'ep': 'small_ester', '6al': 'ald',
    't2h': 'ald', '1-8ol': 'alc', '1-5ol': 'alc', '1-6ol': 'alc', 'ms': 'ms', 'paa': 'paa',
    '-aPine': 'terp', 'euc': 'terp', 'B-myr': 'terp',
}
FRUITY_SMALL_BLOCK = 'fruity-small_ester'

# log10 v/v concentrations: Yang's APL16 panel vs the ORN data available as model input
YANG_CONC: Dict[str, float] = {**{o: -3.0 for o in MEGAMAT12}, 'paa': -3.0, '-aPine': -2.0,
    'euc': -2.0, 'B-myr': -2.0
}
ORN_INPUT_CONC: Dict[str, float] = {**{o: -3.0 for o in MEGAMAT12}, 'paa': -4.0,
    '-aPine': -1.5, 'euc': -1.75, 'B-myr': -2.0
}
APPROX_CONC_ODORS: Dict[str, Tuple[float, float]] = {o: (YANG_CONC[o], ORN_INPUT_CONC[o])
    for o in APL16_ODORS if YANG_CONC[o] != ORN_INPUT_CONC[o]
}
APPROX_CONC_NOTE: str = ('Model ORN input for ' + ', '.join(
    f'{o} ({y:g} in data vs {m:g} in ORN input)' for o, (y, m) in APPROX_CONC_ODORS.items())
    + ' is within 1 log10 unit of, but not at, the concentration Yang used (allowed by '
    "Tom's rule); results involving these odors are approximate."
)

PAIRS16: List[Tuple[str, str]] = [(a, b) for i, a in enumerate(APL16_ODORS)
    for b in APL16_ODORS[i + 1:]
]
MEGAMAT_PAIRS66: List[Tuple[str, str]] = [p for p in PAIRS16
    if p[0] in MEGAMAT12 and p[1] in MEGAMAT12
]

# fixed colors for conditions (same hues as Yang's earlier Claude analysis)
CONDITION_COLORS: Dict[str, str] = {'control': '#2a78d6', 'TNT': '#c8553d',
    'APL on': '#2a78d6', 'APL off': '#c8553d', 'model': '#555555', 'data': '#999999',
}


# ---------------------------------------------------------------------------
# paths
def olfsysm_dir() -> Path:
    """Returns the olfsysm checkout the installed olfsysm was built from.

    ~/src/olfsysm on Tom's hosts, ~/src/kc_queue/olfsysm on lab members' machines.
    """
    return _mb_model.olfsysm_checkout_path()


def pyolfsysm_dir() -> Path:
    """Returns the checkout pyolfsysm, the model's Python code, is imported from.

    The olfsysm repository (the same checkout as `olfsysm_dir()`, both installed editable in
    this repository's venv): ~/src/olfsysm on Tom's hosts, ~/src/kc_queue/olfsysm on lab members'
    machines. Until 2026-11-12 the model's Python code was in al_analysis (`al_analysis_dir()`).
    """
    return Path(_mb_model.__file__).resolve().parents[1]


def pyolfsysm_data_dir() -> Path:
    """Returns the model's data directory (pyolfsysm/data in the olfsysm checkout).

    Yang's and Remy's KC data, Pratyush's connectome tables, the ORN spike deltas and the rest
    of the data the model reads, in al_analysis's data directory until 2026-11-12.
    """
    from pyolfsysm.util import data_root
    return Path(data_root)


def al_analysis_dir() -> Path:
    """Returns the al_analysis checkout the installed al_analysis package comes from.

    ~/src/al_analysis on Tom's hosts, ~/src/kc_queue/al_analysis on lab members' machines (both
    installed editable in this repository's venv, beside it). Found without importing it.
    """
    spec = importlib.util.find_spec('al_analysis')
    if spec is None or not spec.submodule_search_locations:
        raise ModuleNotFoundError('al_analysis is not installed in this venv', name='al_analysis')
    # al_analysis is a namespace package (no __init__.py): its directories instead of a file
    return next(Path(p) for p in spec.submodule_search_locations if Path(p).is_dir()
        ).resolve().parent


def al_analysis_data_dir() -> Path:
    """Returns al_analysis's data directory (al_analysis/data in its checkout): the AL data
    (the consensus ORN means, hemibrain_neuprint, ...); the model's data is in
    `pyolfsysm_data_dir()`."""
    return al_analysis_dir() / 'al_analysis' / 'data'


def script_path(file: Union[str, Path]) -> str:
    """The path of a file of this repository relative to its root (`script_path(__file__)` in
    experiments/E10_joint_fit.py gives 'experiments/E10_joint_fit.py'), for provenance and for git
    commands run from the root (cwd=HERE)."""
    return Path(file).resolve().relative_to(HERE).as_posix()


def hong2p_dir() -> Path:
    """Returns the hong2p checkout the installed hong2p is imported from (~/src/hong2p on Tom's hosts,
    ~/src/kc_queue/hong2p on lab members' machines, installed editable in this repository's venv)."""
    import hong2p
    return Path(hong2p.__file__).resolve().parents[1]


def code_commits() -> Dict[str, object]:
    """Returns the commits and dirty states of this repository, pyolfsysm's, al_analysis's and hong2p's.

    Keys: 'queue_commit' and 'queue_dirty' (this repository, model_experiment_queue: the
    experiment scripts and this harness), 'pyolfsysm_commit' and 'pyolfsysm_dirty' (the model's
    Python code, the olfsysm checkout `pyolfsysm_dir()`), 'al_analysis_commit' and
    'al_analysis_dirty' (al_analysis, whose AL data some scripts read), 'hong2p_commit' and
    'hong2p_dirty' (hong2p, `hong2p_dir()`, where the odor-pair, statistics and imaging functions
    of this harness are since 2026-11-12). Dirty means uncommitted changes to tracked files
    (`hong2p.util.git_commits`). Scripts put these in their provenance.json. Until 2026-11-12 the
    model's Python code was al_analysis's (no 'pyolfsysm_commit') and there was no 'hong2p_commit',
    and before the move to a repository of its own (2026-11-11), 'al_analysis_commit' was the commit
    of the scripts too.
    """
    return _git_commits({'queue': HERE, 'pyolfsysm': pyolfsysm_dir(), 'al_analysis': al_analysis_dir(),
        'hong2p': hong2p_dir()})


def model_provenance(**inputs) -> Dict[str, object]:
    """`mb_model.model_provenance` (olfsysm, pyolfsysm and hong2p) plus this repository and
    al_analysis.

    Args:
        **inputs: passed on (inputs to record by content hash).

    Returns its dict with 'queue' added: `mb_model.repo_info` of this repository (git hash,
    dirty flag, diff of a dirty checkout), and 'al_analysis' (the same, of the al_analysis
    checkout, whose AL data some scripts read) if pyolfsysm did not already record it.
    """
    prov = _mb_model.model_provenance(**inputs)
    prov['queue'] = _mb_model.repo_info(HERE)
    if 'al_analysis' not in prov:
        prov['al_analysis'] = _mb_model.repo_info(al_analysis_dir())
    return prov


def shared_dir() -> Path:
    """The shared directory on tensor (Tom/claude_kc_remapping).

    Set KC_REMAPPING_SHARED to override. Raises FileNotFoundError if none found.
    """
    cands = [os.environ.get('KC_REMAPPING_SHARED')] + [
        f'/{root}/tensor/Tom/claude_kc_remapping' for root in ('local', 'mnt', 'media')
    ]
    for c in cands:
        if c and Path(c).is_dir():
            return Path(c)
    raise FileNotFoundError('shared directory not found; set KC_REMAPPING_SHARED')


def apl_pair_dir() -> Path:
    """<shared>/apl_pair_analysis/: the original report's analysis and the matrices it used
    (kept in the shared directory, not in this repository)."""
    return shared_dir() / 'apl_pair_analysis'


def results_dir(experiment: str) -> Path:
    """<shared>/results/<experiment>/ (compact outputs), created if needed."""
    d = shared_dir() / 'results' / experiment
    d.mkdir(parents=True, exist_ok=True)
    return d


def outputs_dir(experiment: str) -> Path:
    """./outputs/<experiment>/ (large local outputs, gitignored), created if needed."""
    d = HERE / 'outputs' / experiment
    d.mkdir(parents=True, exist_ok=True)
    return d


def data_targets_path() -> Path:
    """The data-target registry, <shared>/results/data_targets/data_targets.json (written by
    data_targets.py)."""
    return shared_dir() / 'results' / 'data_targets' / 'data_targets.json'


_DATA_TARGETS: Dict[str, dict] = {}


def data_target(name: str) -> float:
    """Value of one entry of the data-target registry (e.g. 'ctrl_ceiling', 'block').

    Figure scripts and loss functions read KC / ORN data values with this rather than from
    literals. The registry is read once per process; run kcq/data_targets.py to rebuild it.

    Raises:
        KeyError: `name` is not in the registry.
    """
    if not _DATA_TARGETS:
        _DATA_TARGETS.update(json.loads(data_targets_path().read_text()))
    if name not in _DATA_TARGETS:
        raise KeyError(f'{name!r} is not in the data-target registry ({data_targets_path()})')
    return float(_DATA_TARGETS[name]['value'])


def data_target_ci(name: str) -> Tuple[float, float]:
    """The 90 % interval (lo, hi) of a data-target registry entry that has one (`data_target`).

    Raises:
        KeyError: `name` is not in the registry, or has no interval.
    """
    data_target(name)
    ci = _DATA_TARGETS[name].get('ci')
    if ci is None:
        raise KeyError(f'{name!r} has no interval in the data-target registry')
    return float(ci[0]), float(ci[1])


# ---------------------------------------------------------------------------
# odors and pairs
def odor_abbrev(odor) -> str:
    """'1-5ol @ -3' -> '1-5ol' (last element if a tuple, e.g. (panel, odor))."""
    if isinstance(odor, tuple):
        odor = odor[-1]
    return str(odor).split(' @ ')[0]


def abbrev_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Returns copy with odor columns renamed to abbreviations (must stay unique)."""
    out = df.copy()
    out.columns = pd.Index([odor_abbrev(c) for c in df.columns], name='odor')
    assert not out.columns.duplicated().any(), out.columns[out.columns.duplicated()]
    return out


def odor_order(odors: Sequence[str]) -> List[str]:
    """Orders odors as in APL16_ODORS, then others alphabetically (`hong2p.rdm.order_items`)."""
    return _rdm.order_items(odors, APL16_ODORS)


def pair_list(odors: Sequence[str]) -> List[Tuple[str, str]]:
    """Unordered pairs of odors, each (a, b) with a before b in `odor_order` (`hong2p.rdm.pair_list`)."""
    return _rdm.pair_list(odors, APL16_ODORS)


def pairs_from_square(sq: pd.DataFrame) -> pd.Series:
    """Upper triangle of a symmetric odor x odor DataFrame, as a pair Series (its rows taken in the
    order of its columns; `hong2p.rdm.pairs_from_square`)."""
    sq = abbrev_columns(sq)
    sq.index = sq.columns
    return _rdm.pairs_from_square(sq, APL16_ODORS)


def square_from_pairs(pairs: pd.Series) -> pd.DataFrame:
    """Symmetric odor x odor DataFrame (zero diagonal) from a pair Series (`hong2p.rdm.square_from_pairs`)."""
    return _rdm.square_from_pairs(pairs, APL16_ODORS)


def canonical_pairs(pairs: pd.Series) -> pd.Series:
    """Reindexes a pair Series (keys (a, b) in any order) to `pair_list` order (`hong2p.rdm.canonical_pairs`)."""
    return _rdm.canonical_pairs(pairs, APL16_ODORS)


# subsets two canonical pair Series to their common pairs (same order), dropping NaN
align_pairs = _rdm.align_pairs


# pairs with both odors in `odors`
subset_pairs = _rdm.subset_pairs


# ---------------------------------------------------------------------------
# model input
def apl16_orn_deltas() -> pd.DataFrame:
    """ORN spike-delta input for the 16 APL16 odors: glomeruli x (panel, odor).

    The 12 megamat odors from `megamat_orn_deltas` and paa, -aPine, euc, B-myr from
    `validation2_orn_deltas`, restricted to the 36 glomeruli both have (megamat alone
    has VA4 and VL1 too; validation2 alone has DL4 and VM4). Concentrations are the lab's
    (see APPROX_CONC_ODORS for the three that differ from Yang's).
    """
    from pyolfsysm.mb_model import megamat_orn_deltas, validation2_orn_deltas

    m = megamat_orn_deltas()
    v = validation2_orn_deltas()
    gl = sorted(set(m.index) & set(v.index))

    def pick(df: pd.DataFrame, panel: str, odors: Sequence[str]) -> pd.DataFrame:
        cols = [c for c in df.columns if c[0] == panel and odor_abbrev(c) in odors]
        cols = sorted(cols, key=lambda c: APL16_ODORS.index(odor_abbrev(c)))
        return df.loc[gl, cols]

    orn = pd.concat([pick(m, 'megamat', MEGAMAT12), pick(v, 'validation2', VALIDATION2_4)],
        axis='columns'
    )
    assert [odor_abbrev(c) for c in orn.columns] == APL16_ODORS
    for c in orn.columns:
        assert float(c[1].split(' @ ')[1]) == ORN_INPUT_CONC[odor_abbrev(c)], c
    assert orn.shape[0] == 36 and not orn.isna().any().any()
    return orn


def orn_input_drive(orn_deltas: pd.DataFrame) -> pd.Series:
    """Per-odor mean spike delta across input glomeruli (a measure of total ORN drive)."""
    return abbrev_columns(orn_deltas).mean(axis='index')


def orn_input_pairs(orn_deltas: pd.DataFrame) -> pd.Series:
    """1 - Pearson between odors across input glomeruli."""
    return pair_distances(orn_deltas, measure='pearson')


# ---------------------------------------------------------------------------
# data: Yang's APL16 KC data (from the shared directory), the report's tables
def yang_data_dir() -> Path:
    return shared_dir() / 'data' / 'yang_apl16_kc_olfdb'


_MEASURE_COL = {'cosine': 'one_minus_cosine', 'pearson': 'one_minus_pearson'}


def _recording_rdm_long_csv(binarize: Optional[float]) -> Path:
    """Path of the per-recording long RDM table (computed and cached if `binarize`)."""
    comp = yang_data_dir() / 'computed'
    if binarize is None:
        return comp / 'rdm_per_recording_long.csv'
    path = comp / f'rdm_per_recording_binarized_gt{binarize:g}_long.csv'
    if path.exists():
        return path
    tm = yang_trialmeans()
    rows = []
    keys = ['genotype', 'date', 'fly_num', 'hemibrain']
    for rec, sub in tm.groupby(keys, sort=True):
        mat = sub.pivot(index='roi', columns='odor', values='value')
        mat = (mat > binarize).astype(float)
        for odor_set, odors in (('no_pfo', [o for o in mat.columns if o != SOLVENT]),
                ('with_pfo', list(mat.columns))):
            cos = pair_distances(mat[odors], 'cosine')
            pear = pair_distances(mat[odors], 'pearson')
            d = pd.DataFrame({'one_minus_pearson': pear, 'one_minus_cosine': cos})
            d.index.names = ['odor_a', 'odor_b']
            d = d.reset_index()
            d.insert(0, 'odor_set', odor_set)
            for k, v in reversed(list(zip(keys, rec))):
                d.insert(0, k, v)
            rows.append(d)
    pd.concat(rows, ignore_index=True).to_csv(path, index=False)
    return path


def yang_recording_rdms(genotype: str, measure: str = 'cosine', with_pfo: bool = False,
    binarize: Optional[float] = None) -> pd.DataFrame:
    """Per-recording pair distances: rows (date, fly_num, hemibrain), columns pairs.

    Args:
        genotype: 'TNT3c' (control) or 'TNT' (APL output blocked)
        measure: 'cosine' (1 - cosine similarity; what the report's matrices are) or
            'pearson' (1 - Pearson r), over all KCs of the recording
        with_pfo: include the solvent control as a 17th odor
        binarize: if passed, distances are computed from binarized trial-mean responses
            (response > this threshold; Yang uses 1), the fairer target for binarized
            model `responses`. Cached under computed/.
    """
    df = pd.read_csv(_recording_rdm_long_csv(binarize), dtype={'date': str, 'hemibrain': str})
    df = df[(df.genotype == genotype) &
        (df.odor_set == ('with_pfo' if with_pfo else 'no_pfo'))
    ]
    col = _MEASURE_COL[measure]
    wide = df.pivot_table(index=['date', 'fly_num', 'hemibrain'],
        columns=['odor_a', 'odor_b'], values=col
    )
    # to canonical pair order
    rows = {rec: canonical_pairs(wide.loc[rec]) for rec in wide.index}
    out = pd.DataFrame(rows).T
    out.index = wide.index
    return out


def yang_mean_rdm(genotype: str, measure: str = 'cosine', how: str = 'recording',
    with_pfo: bool = False, binarize: Optional[float] = None) -> pd.Series:
    """Mean pair distances over recordings ('recording') or flies ('fly').

    With how='recording', measure='cosine', binarize=None, this equals the report's
    matrices. See `yang_recording_rdms` for `binarize`.
    """
    recs = yang_recording_rdms(genotype, measure=measure, with_pfo=with_pfo,
        binarize=binarize
    )
    if how == 'recording':
        return recs.mean()
    elif how == 'fly':
        return recs.groupby(level=['date', 'fly_num']).mean().mean()
    raise ValueError(how)


def yang_response_fractions(genotype: str, thresh: float = 1.0) -> pd.DataFrame:
    """Per-recording fraction of KCs with trial-mean response > `thresh` to each odor.

    Rows recordings (date, fly_num, hemibrain), columns odors. Yang's values are z-scored
    fluorescence; she has used thresholds of 0.75 and 1 (see notes/yang-summary.md).
    """
    if thresh == 1.0:
        df = pd.read_csv(yang_data_dir() / 'computed' /
            'response_fraction_per_recording.csv', dtype={'date': str, 'hemibrain': str}
        )
        df = df[df.genotype == genotype]
        wide = df.pivot_table(index=['date', 'fly_num', 'hemibrain'], columns='odor',
            values='frac_responding'
        )
    else:
        tm = yang_trialmeans()
        tm = tm[tm.genotype == genotype].copy()
        tm['responding'] = (tm.value > thresh).astype(float)
        wide = tm.pivot_table(index=['date', 'fly_num', 'hemibrain'], columns='odor',
            values='responding', aggfunc='mean'
        )
    return wide[odor_order(wide.columns)]


def yang_trialmeans() -> pd.DataFrame:
    """Trial-mean response of each KC to each odor, long (all recordings, both genotypes)."""
    return pd.read_parquet(yang_data_dir() / 'computed' / 'trialmean_responses_long.parquet')


def report_matrices() -> Dict[str, pd.Series]:
    """The two matrices the report analyzed, as pair Series: {'TNT3c': ..., 'TNT': ...}."""
    return {g: pairs_from_square(pd.read_csv(apl_pair_dir() / f'{g}_dist_apl16.csv',
        index_col=0)) for g in ('TNT3c', 'TNT')
    }


def orn16_consensus_pairs() -> pd.Series:
    """The 16-odor ORN odor-odor distances: 1 - Pearson r across the 36 glomeruli of the
    model's ORN input, rebuilt from olfdb (datasets 1 and 4) by orn16_distances_from_olfdb.py
    (shared directory, data/model_outputs/inputs/orn16_corrdist_olfdb_36glom.csv). Within-panel
    pairs are the mean over flies of per-fly 1 - r (al_util.mean_of_fly_corrs); cross-panel pairs
    (no fly imaged both panels) use the fly-mean vectors.

    Since 2026-11-06; before, this returned the report's matrix
    (`orn16_consensus_pairs_zero_filled()`), which counted glomeruli a panel did not image as
    zero responses (unresolved_issues.md).
    """
    return pairs_from_square(pd.read_csv(shared_dir() / 'data' / 'model_outputs' / 'inputs' /
        'orn16_corrdist_olfdb_36glom.csv', index_col=0
    ))


def orn16_consensus_pairs_zero_filled() -> pd.Series:
    """The original report's 16-odor ORN matrix (1 - r of consensus glomerular mean dF/F over
    the union of both panels' glomeruli, unimaged glomeruli as 0). Kept for comparison only."""
    return pairs_from_square(pd.read_csv(apl_pair_dir() / 'orn16_consensus_corrdist.csv',
        index_col=0
    ))


def consensus_csv() -> Path:
    """The consensus glomerular mean dF/F of both panels (al_analysis's
    data/sent_to_vlad/2024-06-21/consensus_mean.csv)."""
    return al_analysis_data_dir() / 'sent_to_vlad' / '2024-06-21' / 'consensus_mean.csv'


def consensus_orn_drive() -> pd.Series:
    """Per-odor mean consensus dF/F across glomeruli (the report's 'ORN drive')."""
    cons = pd.read_csv(consensus_csv())
    cons['abbrev'] = cons['odor1'].str.split(' @ ').str[0].replace({'PAA': 'paa'})
    panel = {**{o: 'megamat' for o in MEGAMAT12}, **{o: 'validation2' for o in VALIDATION2_4}}
    gl = [c for c in cons.columns if c not in ('panel', 'odor1', 'abbrev')]
    drive = {}
    for o in APL16_ODORS:
        r = cons[(cons.abbrev == o) & (cons.panel == panel[o])]
        assert len(r) == 1, (o, len(r))
        drive[o] = r.iloc[0][gl].astype(float).fillna(0).mean()
    return pd.Series(drive, name='orn_drive')


def reference_pairs() -> pd.DataFrame:
    """The report's 66-pair table of megamat12 distances in other spaces (shared12_pairs.csv).

    Columns: chem, orn, pn_dend, pn_bout, kc_claw, kc_soma (Remy's wild-type megamat
    KC somata), kc_soma_new, model_uniform, model_hemibrain, model_hemibrain_wd20,
    vcf_raw, vcf_hyp (hyperbolic VCF embedding, mean of 10 seeds), kc_ctrl, kc_tnt.
    """
    df = pd.read_csv(apl_pair_dir() / 'shared12_pairs.csv', index_col=[0, 1])
    df.index.names = ['odor_a', 'odor_b']
    keyed = {p: df.loc[p] for p in df.index}
    canon = canonical_pairs(pd.Series({p: 0.0 for p in keyed}))
    out = pd.DataFrame([keyed[p] if p in keyed else keyed[(p[1], p[0])]
        for p in canon.index], index=canon.index
    )
    return out.astype(float)


# ---------------------------------------------------------------------------
# model outputs -> pair distances
def pair_distances(mat: pd.DataFrame, measure: str = 'cosine', drop_silent: bool = False
    ) -> pd.Series:
    """Odor x odor distances from a (units x odors) matrix, across units (`hong2p.rdm.pair_distances`).

    Args:
        mat: rows units (KCs, glomeruli, ...), columns odors (any strings `odor_abbrev`
            can parse). E.g. `spike_counts` or `responses` from `fit_mb_model`.
        measure: 'cosine' -> 1 - cosine similarity (the data recipe); 'pearson' -> 1 - r
        drop_silent: drop units that are 0 for every odor first (irrelevant for cosine;
            matters for Pearson)

    Returns pair Series (NaN where an odor has no response at all).
    """
    return _rdm.pair_distances(abbrev_columns(mat), measure, drop_silent, APL16_ODORS)


# name -> (source, measure, drop_silent) (pyolfsysm.readout's)
MODEL_RDM_VARIANTS: Dict[str, Tuple[str, str, bool]] = _readout.MODEL_RDM_VARIANTS
PRIMARY_MODEL_RDM = _readout.PRIMARY_MODEL_RDM


def model_rdms(responses: pd.DataFrame, spike_counts: pd.DataFrame) -> Dict[str, pd.Series]:
    """All `MODEL_RDM_VARIANTS` pair Series for one model output (`pyolfsysm.readout.model_rdms`)."""
    return _readout.model_rdms(abbrev_columns(responses), abbrev_columns(spike_counts), APL16_ODORS)


def response_fractions(responses: pd.DataFrame) -> pd.Series:
    """Fraction of model KCs responding to each odor (columns -> abbreviations;
    `pyolfsysm.readout.response_fractions`)."""
    return _readout.response_fractions(abbrev_columns(responses))


# ---------------------------------------------------------------------------
# plotting
# olfdb's (hong-data's) default blue-white-red colormap (in hong2p.viz since 2026-11-12)
oklab_diverging_cmap = _oklab_diverging_cmap


# ---------------------------------------------------------------------------
# statistics
# residual of y around its OLS fit on x, plus fit summary (slope, intercept, r2, resid_sd)
linear_resid = _stats.linear_resid


# residual of y around the best increasing (isotonic) fit on x, plus r2 and resid_sd
isotonic_resid = _stats.isotonic_resid


spearman = _stats.spearman


pearson = _stats.pearson


def block_label(a: str, b: str) -> str:
    """The pair's block under `ODOR_GROUPS` (`hong2p.rdm.block_label`)."""
    return _rdm.block_label(a, b, ODOR_GROUPS)


def block_means(pairs: pd.Series) -> pd.Series:
    """Mean of a pair Series within each between-group block ('within' = same group;
    `hong2p.rdm.block_means` with `ODOR_GROUPS`)."""
    return _rdm.block_means(pairs, ODOR_GROUPS)


def per_odor_means(pairs: pd.Series) -> pd.Series:
    """Mean over the pairs containing each odor (`hong2p.rdm.per_odor_means`)."""
    return _rdm.per_odor_means(pairs, APL16_ODORS)


def compare_on_off(on: pd.Series, off: pd.Series, *, orn: Optional[pd.Series] = None,
    drive: Optional[pd.Series] = None, gap_ceiling: float = 0.99
    ) -> Tuple[Dict[str, float], Dict[str, pd.DataFrame]]:
    """The report's statistics relating an APL-on (control) to an APL-off (TNT) matrix.

    Args:
        on, off: pair distance Series (1 - similarity) with APL intact / blocked
        orn: ORN pair distances (e.g. `orn16_consensus_pairs()`), for the
            ORN-agreement statistics
        drive: per-odor ORN drive (e.g. `consensus_orn_drive()`), for the per-odor
            APL effect vs drive correlation

    Returns:
        stats: flat dict. Keys include mean_on / mean_off, spearman_on_off,
            lin_r2 / iso_r2 (variance of `on` explained by a linear / any monotone
            function of `off`), lin_resid_sd / iso_resid_sd, block_fruity_small_iso
            (mean isotonic residual of the fruity-ester x small-ester block; the data
            give about +0.055), with `orn`: spearman_orn_on / _off, gap_closed
            (fraction of the gap from spearman_orn_on to `gap_ceiling` closed by APL
            block), apl_share_resid_var (1 - var(resid_off) / var(resid_on) around the
            linear ORN prediction), spearman_resid_orn_on_off, and with `drive`:
            spearman_per_odor_iso_vs_drive.
        tables: 'pairs' (per-pair on, off, delta, lin_resid, iso_resid, [orn,
            resid_orn_on, resid_orn_off, apl_eff]), 'blocks', 'per_odor'.
    """
    on, off = align_pairs(on, off)
    out: Dict[str, float] = dict(n_pairs=len(on), mean_on=float(on.mean()),
        mean_off=float(off.mean()), mean_sim_on=float(1 - on.mean()),
        mean_sim_off=float(1 - off.mean()), spearman_on_off=spearman(on, off),
        pearson_on_off=pearson(on, off),
    )
    lin_res, lin = linear_resid(on, off)
    iso_res, iso = isotonic_resid(on, off)
    out.update(lin_r2=lin['r2'], lin_resid_sd=lin['resid_sd'], lin_slope=lin['slope'],
        lin_intercept=lin['intercept'], iso_r2=iso['r2'], iso_resid_sd=iso['resid_sd'],
    )
    pairs = pd.DataFrame({'on': on, 'off': off, 'delta': on - off, 'lin_resid': lin_res,
        'iso_resid': iso_res
    })
    labels = pd.Series([block_label(a, b) for a, b in on.index], index=on.index)
    blocks = pd.DataFrame({'iso_resid': block_means(iso_res),
        'lin_resid': block_means(lin_res), 'on': block_means(on), 'off': block_means(off),
        'n_pairs': labels.value_counts(),
    })
    out['block_fruity_small_iso'] = float(blocks.iso_resid.get(FRUITY_SMALL_BLOCK, np.nan))
    out['block_fruity_small_lin'] = float(blocks.lin_resid.get(FRUITY_SMALL_BLOCK, np.nan))
    between = [b for b in blocks.index if b not in ('within', FRUITY_SMALL_BLOCK)]
    out['block_other_between_iso'] = float(blocks.loc[between, 'iso_resid'].mean())
    per_odor = pd.DataFrame({'on': per_odor_means(on), 'off': per_odor_means(off),
        'delta': per_odor_means(on - off), 'iso_resid': per_odor_means(iso_res),
        'lin_resid': per_odor_means(lin_res),
    })
    if orn is not None:
        orn_a = orn.loc[orn.index.intersection(on.index)]
        pairs['orn'] = orn_a
        out['spearman_orn_on'] = spearman(orn_a, on)
        out['spearman_orn_off'] = spearman(orn_a, off)
        out['pearson_orn_on'] = pearson(orn_a, on)
        out['pearson_orn_off'] = pearson(orn_a, off)
        r_on, f_on = linear_resid(on, orn_a)
        r_off, f_off = linear_resid(off, orn_a)
        pairs['resid_orn_on'] = r_on
        pairs['resid_orn_off'] = r_off
        pairs['apl_eff'] = lin_res
        out['resid_orn_sd_on'] = f_on['resid_sd']
        out['resid_orn_sd_off'] = f_off['resid_sd']
        out['apl_share_resid_var'] = float(1 - r_off.var() / r_on.var())
        out['gap_closed'] = float((out['spearman_orn_off'] - out['spearman_orn_on']) /
            (gap_ceiling - out['spearman_orn_on'])
        )
        out['spearman_resid_orn_on_off'] = spearman(r_on, r_off)
        out['spearman_apl_eff_vs_resid_orn_on'] = spearman(lin_res, r_on)
    if drive is not None:
        d = drive.reindex(per_odor.index)
        per_odor['drive'] = d
        ok = d.notna()
        out['spearman_per_odor_iso_vs_drive'] = float(stats.spearmanr(per_odor.iso_resid[ok],
            d[ok])[0]
        )
        out['spearman_per_odor_delta_vs_drive'] = float(stats.spearmanr(per_odor.delta[ok],
            d[ok])[0]
        )
    return out, {'pairs': pairs, 'blocks': blocks, 'per_odor': per_odor}


REFERENCE_COLS: List[str] = ['orn', 'pn_dend', 'pn_bout', 'kc_claw', 'kc_soma', 'kc_soma_new',
    'model_uniform', 'model_hemibrain_wd20', 'chem', 'vcf_raw', 'vcf_hyp'
]


def compare_to_references(pairs: pd.Series, ref: Optional[pd.DataFrame] = None
    ) -> Dict[str, float]:
    """Spearman of a pair Series with each megamat12 reference space (66 pairs), and the
    residual-vs-VCF statistics of the report.

    Keys: spearman_<ref> for `REFERENCE_COLS`; resid_orn_sd; spearman_resid_vs_vcf_hyp,
    spearman_resid_vs_vcf_raw, spearman_resid_vs_chem; spearman_resid_vs_remap_wt
    (deviations from the ORN prediction vs those of Remy's wild-type KC somata),
    spearman_resid_vs_remap_ctrl / _tnt (vs Yang's control / TNT deviations); mean_66.
    """
    if ref is None:
        ref = reference_pairs()
    p = pairs.loc[pairs.index.intersection(ref.index)]
    ref = ref.loc[p.index]
    out = {'n_pairs': len(p), 'mean_66': float(p.mean())}
    for c in REFERENCE_COLS:
        out[f'spearman_{c}'] = spearman(p, ref[c])
    resid, fit = linear_resid(p, ref['orn'])
    out['resid_orn_sd'] = fit['resid_sd']
    for c in ('vcf_hyp', 'vcf_raw', 'chem'):
        out[f'spearman_resid_vs_{c}'] = spearman(resid, ref[c])
    for name, col in (('wt', 'kc_soma'), ('ctrl', 'kc_ctrl'), ('tnt', 'kc_tnt')):
        r_ref, _ = linear_resid(ref[col], ref['orn'])
        out[f'spearman_resid_vs_remap_{name}'] = spearman(resid, r_ref)
    return out


def compare_model_to_data(model: pd.Series, data: pd.Series, *,
    orn: Optional[pd.Series] = None) -> Dict[str, float]:
    """How well a model pair-distance Series matches a data one (same condition).

    Keys: n_pairs, spearman, pearson, mean_model, mean_data, rmse, and with `orn`:
    spearman_dev (Spearman between model and data deviations from the linear ORN
    prediction) and pearson_dev, dev_sd_model, dev_sd_data.
    """
    m, d = align_pairs(model, data)
    out = dict(n_pairs=len(m), spearman=spearman(m, d), pearson=pearson(m, d),
        mean_model=float(m.mean()), mean_data=float(d.mean()),
        rmse=float(np.sqrt(np.mean((m - d)**2))),
    )
    if orn is not None:
        o = orn.loc[orn.index.intersection(m.index)]
        rm, fm = linear_resid(m, o)
        rd, fd = linear_resid(d, o)
        out.update(spearman_dev=spearman(rm, rd), pearson_dev=pearson(rm, rd),
            dev_sd_model=fm['resid_sd'], dev_sd_data=fd['resid_sd'],
        )
    return out


def apl_effect(on: pd.Series, off: pd.Series) -> pd.Series:
    """Per-pair APL effect on odor-odor distances: the APL-on distance minus the APL-off
    distance of each pair (canonical pair order; pairs missing from either side dropped;
    `pyolfsysm.effects.pair_effect`)."""
    return _effects.pair_effect(on, off, APL16_ODORS, name='apl_effect')


_YANG_APL_EFFECT: Dict[tuple, pd.Series] = {}


YANG_EFFECT_VALUES = ('graded', 'thresholded', 'noise_corrected', 'nonneg', 'floor')


def yang_recording_rdms_values(genotype: str, values: str = 'graded') -> pd.DataFrame:
    """Per-recording cosine pair distances (rows recordings, columns pairs) of the KC data
    transformed as `values` says: 'graded' (as `yang_recording_rdms`), 'thresholded' (trial
    means <= 1 set to 0 before the distance, as a model KC that does not respond contributes 0),
    'noise_corrected' (each recording's cosine corrected for trial-to-trial noise,
    `kcq.noise_floor.analyse`: the distance a noise-free recording would show),
    'nonneg' (negative trial means set to 0) or 'floor' (trial means at or below the noise-defined
    floor `kc_noise_floor()` set to 0: twice the sd of a trial mean's noise, from the control
    flies' trial-to-trial sd of the non-responding entries, the same floor for both genotypes)."""
    if values == 'graded':
        return yang_recording_rdms(genotype, 'cosine')
    if values in ('thresholded', 'nonneg', 'floor'):
        thresh = {'thresholded': KC_RESPONSE_Z, 'nonneg': 0.0}.get(values)
        if thresh is None:
            thresh = kc_noise_floor()
        return recording_rdms_from_matrices({k: transform_matrix(m, values='thresholded', thresh=thresh)
            for k, m in yang_recording_matrices(genotype).items()})
    if values == 'noise_corrected':
        # a module of this package that imports the harness, hence imported here
        from kcq.noise_floor import analyse, load_trials, recording_arrays
        rows = {}
        for rec, by_odor in recording_arrays(load_trials(genotype)).items():
            _, pr = analyse(by_odor, 'graded')
            rows[rec] = canonical_pairs(pr['corrected'])
        return _stack_recordings(rows)
    raise ValueError(f'unknown values {values!r} (one of {YANG_EFFECT_VALUES})')


def yang_apl_effect(measure: str = 'cosine', binarize: Optional[float] = None,
    values: str = 'graded') -> pd.Series:
    """The KC data's per-pair APL effect: the control flies' (TNT3c) distance minus the
    distance with APL output blocked (TNT), each the mean over recordings (`yang_mean_rdm`).
    `binarize` (e.g. 1, Yang's responder threshold on z-scored trial means) computes it from
    binarized responses instead, the like-for-like target of a binarized model readout
    (`readout(..., 'binary')`). `values` 'thresholded', 'noise_corrected', 'nonneg' or 'floor'
    (cosine only, without `binarize`) computes it from the data transformed as
    `yang_recording_rdms_values` says.
    Cached per (measure, binarize, values)."""
    key = (measure, binarize, values)
    if key not in _YANG_APL_EFFECT:
        if values == 'graded':
            _YANG_APL_EFFECT[key] = apl_effect(yang_mean_rdm('TNT3c', measure, binarize=binarize),
                yang_mean_rdm('TNT', measure, binarize=binarize))
        else:
            if measure != 'cosine' or binarize is not None:
                raise ValueError('values other than graded are computed for cosine distances without binarize')
            _YANG_APL_EFFECT[key] = apl_effect(yang_recording_rdms_values('TNT3c', values).mean(),
                yang_recording_rdms_values('TNT', values).mean())
    return _YANG_APL_EFFECT[key]


def apl_effect_stats(on: pd.Series, off: pd.Series, data: Optional[pd.Series] = None
    ) -> Dict[str, float]:
    """How well a model's per-pair APL effect (APL-on minus APL-off distance) matches the
    KC data's, over the pairs both have (the 120 pairs of the APL16 panel;
    `pyolfsysm.effects.effect_stats`).

    Args:
        on, off: the model's pair distances with APL on and with APL off (acute), cosine
            distance on spike counts as everywhere else.
        data: the data's per-pair APL effect; default `yang_apl_effect()` (control flies
            minus APL output blocked (TNT), cosine, mean over recordings).

    Keys (`apl_effect_` prefix): `rmse`, the APL-effect loss (root mean squared difference
    between the model's and the data's per-pair effect; lower is better); `rmse_centered`, the
    same with each side's mean removed (the pattern part: rmse^2 = mean_error^2 +
    rmse_centered^2); `mean` (the model's mean effect), `mean_error` (model minus data),
    `spearman` and `pearson` (the per-pair pattern), `n_pairs`.
    """
    data = yang_apl_effect() if data is None else data
    return _effects.effect_stats(on, off, data, APL16_ODORS, prefix='apl_effect_')


# the additive design over odors of a set of odor pairs (one row per pair, 1 for its two odors)
pair_design = _rdm.pair_design


# each odor's share of a per-pair quantity under value(a, b) = s_a + s_b, and the variance explained
odor_shares = _rdm.odor_shares


def apl_effect_pair_parts(on: pd.Series, off: pd.Series, data: Optional[pd.Series] = None) -> pd.DataFrame:
    """Per pair, the parts of the APL-effect loss (`apl_effect_decomposition`): columns `model` and
    `data` (each side's per-pair effect), `error` (model minus data, each with its mean removed: the
    pattern part's error), `per_odor` (its projection onto the additive space over odors) and
    `pair_specific` (the residual; `error = per_odor + pair_specific`). Index the pairs both sides
    have; `data` defaults to `yang_apl_effect()` (`pyolfsysm.effects.effect_pair_parts`)."""
    data = yang_apl_effect() if data is None else data
    return _effects.effect_pair_parts(on, off, data, APL16_ODORS)


CURVE_EPS = 1e-6


def power_curve(s_off: np.ndarray, s_on: np.ndarray) -> Tuple[float, float]:
    """Least-squares fit of `s_on = a * s_off ** b` over pairs (similarities, one minus cosine distance; values below
    CURVE_EPS enter as CURVE_EPS), a and b positive (E77). Returns (a, b)."""
    from scipy.optimize import curve_fit
    x = np.maximum(np.asarray(s_off, dtype=float), CURVE_EPS)
    y = np.asarray(s_on, dtype=float)
    (a, b), _ = curve_fit(lambda x, a, b: a * x ** b, x, y, p0=(1.0, 1.0), bounds=([1e-9, 1e-9], [np.inf, np.inf]), maxfev=20000)
    return float(a), float(b)


def apply_power_curve(s_off, curve) -> np.ndarray:
    """`a * s_off ** b` for `curve` (a, b), or per pair for an (n, 2) array of curves."""
    c = np.asarray(curve, dtype=float)
    x = np.maximum(np.asarray(s_off, dtype=float), CURVE_EPS)
    if c.ndim == 1:
        return c[0] * x ** c[1]
    return c[:, 0] * x ** c[:, 1]


def apl_effect_inherited_split(m_on: pd.Series, m_off: pd.Series, d_on: pd.Series, d_off: pd.Series,
                               data_curve=None, model_curve=None) -> pd.DataFrame:
    """A model's per-pair APL-effect error split into what it inherits from its APL-off distances, its on-against-off
    curve and a pair-specific part (E77). With similarities s = 1 - distance and monotone curves s_on ~ f(s_off),
    f(x) = a x^b fitted to the data (f_d) and to the model (f_m), and g_d(x) = x - f_d(x):

        error = (m_on - m_off) - (d_on - d_off)
              = [g_d(s_off_m) - g_d(s_off_d)]          inherited: the data's own curve applied to both APL-off sets
              + [f_d(s_off_m) - f_m(s_off_m)]          curve: the model's curve against the data's
              + [r_d - r_m]                             pair-specific: each side's residual about its own curve,

    with r = s_on - f(s_off), exactly (to rounding). The `_c` columns are the same terms with their means over pairs
    removed (the pattern part's error is the sum of those).

    Args:
        m_on, m_off, d_on, d_off: pair distances (any pair order; aligned in canonical order).
        data_curve, model_curve: (a, b), or one (a, b) per pair as an (n, 2) array in the canonical pair order (e.g. a
            data curve fitted leaving each pair's odors out); fitted in sample when None.

    Returns:
        DataFrame by pair: error, inherited, curve, pair_specific and their `_c` versions; `attrs` holds the two curves
        when fitted here.
    """
    mo, mf = align_pairs(canonical_pairs(m_on), canonical_pairs(m_off))
    do, df_ = align_pairs(canonical_pairs(d_on), canonical_pairs(d_off))
    mo, do = align_pairs(mo, do)
    mf, df_ = mf.reindex(mo.index), df_.reindex(mo.index)
    som, sfm, sod, sfd = (1 - x.to_numpy(dtype=float) for x in (mo, mf, do, df_))
    attrs = {}
    if data_curve is None:
        data_curve = power_curve(sfd, sod)
        attrs['data_curve'] = data_curve
    if model_curve is None:
        model_curve = power_curve(sfm, som)
        attrs['model_curve'] = model_curve
    fd_m, fd_d = apply_power_curve(sfm, data_curve), apply_power_curve(sfd, data_curve)
    fm_m = apply_power_curve(sfm, model_curve)
    inherited = (sfm - fd_m) - (sfd - fd_d)
    curve = fd_m - fm_m
    pair_specific = (sod - fd_d) - (som - fm_m)
    error = (mo.to_numpy(dtype=float) - mf.to_numpy(dtype=float)) - (do.to_numpy(dtype=float) - df_.to_numpy(dtype=float))
    out = pd.DataFrame({'error': error, 'inherited': inherited, 'curve': curve, 'pair_specific': pair_specific}, index=mo.index)
    for c in list(out.columns):
        out[f'{c}_c'] = out[c] - out[c].mean()
    out.attrs.update(attrs)
    return out


def inherited_shares(split: pd.DataFrame, mask=None) -> Dict[str, float]:
    """Each term's share of the squared pattern error (its projection on the centered error, over `mask`'s pairs if
    given; the three add to 1) and the pattern part's rms (E77)."""
    e = split['error_c'].to_numpy(dtype=float)
    m = np.ones(len(e), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    den = float(e[m] @ e[m])
    return {**{f'share_{k}': float(split[f'{k}_c'].to_numpy(dtype=float)[m] @ e[m]) / den for k in ('inherited', 'curve', 'pair_specific')},
            'pattern_rms': float(np.sqrt(np.mean(e[m] ** 2)))}


def apl_effect_decomposition(on: pd.Series, off: pd.Series, data: Optional[pd.Series] = None
    ) -> Tuple[Dict[str, float], pd.DataFrame]:
    """The APL-effect loss (`apl_effect_stats`) split into three parts whose squares add:
    `rmse**2 = mean_error**2 + per_odor**2 + pair_specific**2` (`pyolfsysm.effects.effect_decomposition`).
    The pattern part's error (model minus data, each with its mean removed) is projected onto the
    additive space over odors (a pair's value the sum of its two odors' shares, `odor_shares`): the
    projection's root mean square is the per-odor part, the residual's the pair-specific part.

    Args:
        on, off: the model's pair distances with APL on and APL off (any readout, e.g.
            `readout_pair_distances`).
        data: the data's per-pair effect; default `yang_apl_effect()`.

    Returns (dict with the `apl_effect_stats` keys plus `apl_effect_per_odor`,
    `apl_effect_pair_specific`, `apl_effect_share_pearson` (the model's per-odor shares against the
    data's), `apl_effect_model_additive_r2` and `apl_effect_data_additive_r2` (variance over pairs
    each side's additive model explains); DataFrame indexed by odor with the model's and the data's
    shares and their difference).
    """
    data = yang_apl_effect() if data is None else data
    return _effects.effect_decomposition(on, off, data, APL16_ODORS, prefix='apl_effect_')


def apl_effect_spread_shape(on: pd.Series, off: pd.Series, data: Optional[pd.Series] = None
    ) -> Dict[str, float]:
    """The pattern part of the APL-effect loss (`apl_effect_stats`' `apl_effect_rmse_centered`) split
    exactly into an amplitude part and a shape part, `rmse_centered**2 = amplitude**2 + shape**2`
    (`pyolfsysm.effects.effect_spread_shape`; Yang's request of 2026-11-13, from
    analyses/apl_effect_metrics.py): `apl_effect_amplitude` = |sd_model - sd_data| (the model's
    per-pair effects too spread out or too compressed over the pairs), `apl_effect_shape` = sqrt(2
    sd_model sd_data (1 - r)) (a different pattern at matched spread; r Pearson's), and
    `apl_effect_sd_model`, `apl_effect_sd_data`, `apl_effect_sd_ratio`. `data` defaults to
    `yang_apl_effect()`."""
    data = yang_apl_effect() if data is None else data
    return _effects.effect_spread_shape(on, off, data, APL16_ODORS, prefix='apl_effect_')


# ---------------------------------------------------------------------------
# readouts of the model's spike counts (E65): the model's side of the per-pair APL effect is a
# cosine distance on raw spike counts, the data's on z-scored GCaMP fluorescence, which
# compresses large responses; these map a (KCs x odors) spike-count matrix onto a compressed
# signal before the distances are taken
READOUT_KINDS = _readout.READOUT_KINDS


# a readout of a model's spike counts (rows KCs, columns odors), elementwise (pyolfsysm.readout)
readout = _readout.readout


def readout_pair_distances(counts: pd.DataFrame, kind: str = 'raw', param=None, floor: float = 0.0
    ) -> pd.Series:
    """Cosine odor-odor distances across KCs after `readout` (`floor` as there): the model's
    side of the per-pair APL effect under that readout ('raw' gives the `spike_counts_cosine`
    distances; `pyolfsysm.readout.readout_pair_distances`)."""
    return _readout.readout_pair_distances(abbrev_columns(counts), kind, param, floor, APL16_ODORS)


def apl_effect_stats_readout(on_counts: pd.DataFrame, off_counts: pd.DataFrame, kind: str = 'raw',
    param=None, data: Optional[pd.Series] = None) -> Dict[str, float]:
    """`apl_effect_stats` of a model whose APL-on and acute APL-off spike counts are read out
    the same way (`readout`, one `kind` and `param` for both conditions, so the APL effect stays
    a prediction). `data` defaults to the KC data's effect, or for kind 'binary' to the effect
    of the binarized data (`yang_apl_effect(binarize=1.0)`). The keys of `apl_effect_stats`,
    plus `mean_dist_on` and `mean_dist_off` (the mean pair distance in each condition)."""
    if data is None:
        data = yang_apl_effect(binarize=1.0) if kind == 'binary' else yang_apl_effect()
    on = readout_pair_distances(on_counts, kind, param)
    off = readout_pair_distances(off_counts, kind, param)
    out = apl_effect_stats(on, off, data)
    out.update(mean_dist_on=float(on.mean()), mean_dist_off=float(off.mean()))
    return out


# how concentrated each odor's response vector is across units (pyolfsysm.kc_stats)
amplitude_stats = _kc_stats.amplitude_stats


responder_cv = _kc_stats.responder_cv


readout_stat = _readout.readout_stat


# the parameter of a one-parameter readout at which the APL-on readout_stat equals a target
calibrate_readout = _readout.calibrate_readout


rank_match_map = _readout.rank_match_map


# ---------------------------------------------------------------------------
# the flies' recording protocol (E82): the flies' response is the mean calcium signal 2.0 to 8.1 s
# after the odor's onset minus the 5 s before it (a 2 s pulse, nuclear GCaMP6s), so a model spike
# at time s counts by the area inside that window of the indicator kernel it starts,
# (1 - exp(-(t - s) / rise)) exp(-(t - s) / decay) for t > s.

# each spike's weight under the flies' protocol, and a model's counts read through it (E82; pyolfsysm.readout)
kernel_window_weights = _readout.kernel_window_weights


kernel_weighted_counts = _readout.kernel_weighted_counts


# ---------------------------------------------------------------------------
# the flies' traces (E84): per-trial excerpts of a KC's signal around each odor onset, the response
# in any window after the onset (Yang's is the mean over volumes round(2 v) .. round(8 v) - 1 after the
# onset volume minus the mean over the round(5 v) volumes before it, v the volume rate), and the
# indicator kernel (1 - exp(-u / rise)) exp(-u / decay) fitted to mean traces under a drive (an onset
# burst, a steady drive over the pulse, or a histogram of firing).

window_volumes = _imaging.window_volumes


trace_excerpts = _imaging.trace_excerpts


windowed_responses = _imaging.windowed_responses


indicator_kernel = _imaging.indicator_kernel


indicator_kernel_integral = _imaging.indicator_kernel_integral


DRIVES = _imaging.DRIVES


_drive_bins = _imaging._drive_bins


kernel_signal = _imaging.kernel_signal


fit_indicator_kernel = _imaging.fit_indicator_kernel


# ---------------------------------------------------------------------------
# the low end of the readout (E68): measurement noise below the data's response criterion, the
# saturating readout's one-spike floor (`readout(..., floor=b)`), the conversion of readout units
# to the data's z units, and the split of a cosine into responder / sub-threshold terms. The
# data's criterion is a trial mean (of 3 trials) above 1 z; its measurement noise is the
# trial-to-trial sd of the non-responding entries (`kc_trial_noise_sd`).
KC_RESPONSE_Z = 1.0
LOW_END_BAND = (0.5, 1.0)


def expected_pair_distances_noisy(mat: pd.DataFrame, noise_sd: float, n_trials: int = 1) -> pd.Series:
    """Cosine distances of a (units x odors) matrix, in expectation over independent zero-mean
    measurement noise of sd `noise_sd` per entry and trial, after averaging `n_trials` trials
    (`pyolfsysm.readout.expected_pair_distances_noisy`). `noise_sd` 0 gives
    `pair_distances(mat, 'cosine')`.
    """
    return _readout.expected_pair_distances_noisy(abbrev_columns(mat), noise_sd, n_trials, APL16_ODORS)


# pseudo-recordings of a model's readout through measurement noise (pyolfsysm.readout)
sample_pseudo_recordings = _readout.sample_pseudo_recordings


def noise_corrected_pair_distances(trials: np.ndarray, odors: Sequence[str],
    min_reliability: Optional[float] = None) -> pd.Series:
    """The noise-corrected cosine distances of one recording from its trials (n_trials, units,
    odors), as `kcq.noise_floor.analyse(..., 'graded')` computes them for the data
    (`hong2p.imaging.noise_corrected_pair_distances`): each odor's signal power is the squared length
    of its trial mean minus the trial-to-trial variance summed over units divided by the number of
    trials (floored at 1e-12), and the corrected distance is 1 - (dot product of the trial means) /
    sqrt(signal power a x signal power b). Returns a canonical pair Series.

    `min_reliability`: pairs with an odor whose signal power is less than this share of its trial
    mean's squared length are NaN (the correction divides by a near-zero signal there; a model odor
    with almost no responders). The KC data's shares are all above 0.5, so 0.1 leaves them intact.
    """
    return _imaging.noise_corrected_pair_distances(trials, odors, min_reliability, APL16_ODORS)


def low_end_stats(trial_signal: np.ndarray, noise_sd: float, *, band: Tuple[float, float] = LOW_END_BAND,
    thresh: float = KC_RESPONSE_Z, weights: Optional[np.ndarray] = None, n_units: Optional[np.ndarray] = None
    ) -> Dict[str, float]:
    """The control-fly statistics the low end is calibrated on, in expectation over Gaussian
    measurement noise (closed form, so a calibration is smooth and needs no sampling): the
    response rate (share of entries whose trial mean exceeds `thresh`), the share of entries
    whose trial mean lies in `band` (lo, hi], and the responders-only participation ratio
    ((sum x^2)^2 / (n sum x^4) over the entries above `thresh`, others 0; the ratio of the
    expectations of the two sums, per group), each the mean over groups (odors, or recordings x
    odors).

    Args:
        trial_signal: (groups, entries) noise-free trial means in z (an entry's mean over trials of
            its readout times the conversion c).
        noise_sd: sd of the noise on a trial mean (the per-trial sd / sqrt(n_trials)); 0 for none.
        weights: (groups, entries) multiplicities of the entries (identical entries grouped);
            default 1.
        n_units: (groups,) number of units per group (the participation ratio's n and the shares'
            denominator); default the summed weights.

    Returns dict with 'rate', 'band_share', 'participation_ratio_resp'.
    """
    mu = np.atleast_2d(np.asarray(trial_signal, dtype=float))
    w = np.ones_like(mu) if weights is None else np.atleast_2d(np.asarray(weights, dtype=float))
    n = w.sum(axis=1) if n_units is None else np.asarray(n_units, dtype=float)
    lo, hi = band
    if noise_sd > 0:
        from scipy.special import ndtr
        sd = float(noise_sd)
        a = (thresh - mu) / sd
        tail = ndtr(-a)
        p_resp = tail
        p_band = ndtr((hi - mu) / sd) - ndtr((lo - mu) / sd)
        # partial moments of the standard normal above a: E[Z^j 1{Z > a}], j = 0..4
        phi = np.exp(-0.5 * a ** 2) / np.sqrt(2 * np.pi)
        m = [tail, phi, a * phi + tail, (a ** 2 + 2) * phi, (a ** 3 + 3 * a) * phi + 3 * tail]
        # E[(mu + sd Z)^k 1{Z > a}] for k = 2, 4, binomially from the partial moments
        e2 = mu ** 2 * m[0] + 2 * mu * sd * m[1] + sd ** 2 * m[2]
        e4 = (mu ** 4 * m[0] + 4 * mu ** 3 * sd * m[1] + 6 * mu ** 2 * sd ** 2 * m[2] + 4 * mu * sd ** 3 * m[3]
            + sd ** 4 * m[4])
    else:
        resp = mu > thresh
        p_resp = resp.astype(float)
        p_band = ((mu > lo) & (mu <= hi)).astype(float)
        e2, e4 = np.where(resp, mu ** 2, 0.0), np.where(resp, mu ** 4, 0.0)
    s2, s4 = (w * e2).sum(axis=1), (w * e4).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        pr = np.where(s4 > 0, s2 ** 2 / (n * s4), np.nan)
    return dict(rate=float(np.mean((w * p_resp).sum(axis=1) / n)), band_share=float(np.mean((w * p_band).sum(axis=1) / n)),
        participation_ratio_resp=float(np.nanmean(pr)))


def _unique_groups(mu: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(groups x entries) values -> (values, weights, n_units), identical values of a group merged
    (zero weights pad the groups to one width)."""
    groups = [np.unique(row, return_counts=True) for row in mu]
    width = max(len(c) for _, c in groups)
    vals, w = np.zeros((len(groups), width)), np.zeros((len(groups), width))
    for g, (u, c) in enumerate(groups):
        vals[g, :len(c)], w[g, :len(c)] = u, c
    return vals, w, np.array([mu.shape[1]] * len(groups), dtype=float)


def _solve_c(grouped, target_rate: float, sd_mean: float, thresh: float = KC_RESPONSE_Z,
    bounds: Tuple[float, float] = (1e-3, 1e4)) -> float:
    """The conversion c at which the expected response rate of c x `grouped` (values, weights,
    n_units) through noise of sd `sd_mean` on a trial mean equals `target_rate`."""
    from scipy.optimize import brentq
    vals, w, n = grouped
    f = lambda lc: low_end_stats(np.exp(lc) * vals, sd_mean, thresh=thresh, weights=w, n_units=n)['rate'] - target_rate
    lo, hi = np.log(bounds[0]), np.log(bounds[1])
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return float(np.exp(lo) if abs(flo) < abs(fhi) else np.exp(hi))
    return float(np.exp(brentq(f, lo, hi, xtol=1e-12, rtol=1e-12)))


def z_scale(signal_on, target_rate: float, noise_sd: float, *, n_trials: int = 3, thresh: float = KC_RESPONSE_Z,
    bounds: Tuple[float, float] = (1e-3, 1e4)) -> float:
    """The conversion c of readout units to the data's z units: the factor at which the APL-on
    response rate under the data's criterion (trial mean of `n_trials` trials above `thresh` z,
    through independent Gaussian measurement noise of sd `noise_sd` per trial; `low_end_stats`)
    equals `target_rate` (e.g. `data_target('ctrl_frac')`, the control flies' 6.8 %).

    Args:
        signal_on: the APL-on readout (units x odors DataFrame or array), or per-trial readouts
            (n_recordings, n_trials, units, odors), whose trial means are taken first.

    Root of the (monotone) expected rate in log c (Brent); returns the nearer bound when the target
    is not bracketed.
    """
    sig = signal_on.to_numpy(dtype=float) if isinstance(signal_on, pd.DataFrame) else np.asarray(signal_on, dtype=float)
    if sig.ndim == 4:
        n_trials = sig.shape[1]
        sig = sig.mean(axis=1)                              # recordings x units x odors
        mu = np.concatenate([sig[r].T for r in range(sig.shape[0])])
    else:
        mu = sig.T                                          # odors x units
    return _solve_c(_unique_groups(mu), target_rate, float(noise_sd) / np.sqrt(n_trials), thresh, bounds)


def calibrate_low_end(on_counts, *, noise_sd: float, target_rate: float, target_band: float,
    target_pr: float, n_trials: int = 3, fit_floor: bool = True, K0: float = 10.0, floor0: float = 0.0,
    K_bounds: Tuple[float, float] = (1e-3, 1e3), floor_bounds: Tuple[float, float] = (-1.0, 50.0),
    c_bounds: Tuple[float, float] = (1e-2, 1e3), tol: float = 1e-5, max_iter: int = 50,
    band: Tuple[float, float] = LOW_END_BAND, thresh: float = KC_RESPONSE_Z) -> Dict[str, object]:
    """The saturating readout with a low end, calibrated jointly on the APL-on counts alone, so the
    APL-off amplitudes and the APL effect stay predictions: the half-saturation K, the one-spike
    floor b (`readout(..., 'saturating', K, floor=b)`) and the conversion c of readout units to z,
    through measurement noise of sd `noise_sd` per trial (`low_end_stats`, closed form), such that
    the APL-on response rate under the data's criterion equals `target_rate` (sets c), the share of
    entries with a trial mean in `band` equals `target_band` (sets b) and the responders-only
    participation ratio equals `target_pr` (sets K). Iterated to a fixed point: each round solves b
    (c re-solved at every b) and then K (c re-solved at every K) until K and b change by less than
    `tol` (relative for K, absolute for b). The band share need not be monotone in b (with c
    re-solved), so b is the root nearest to `floor0` among the sign changes on a grid over
    `floor_bounds`; K is the root of the (monotone) participation ratio in log K.

    Args:
        on_counts: APL-on spike counts, (units x odors) DataFrame or array, or per-trial counts
            (n_recordings, n_trials, units, odors) (jittered trials; `n_trials` is then taken from
            the array).
        fit_floor: False holds b at `floor0` (K and c still calibrated jointly).
        c_bounds: the range of c; floors at which the target rate needs c outside it are
            excluded, and when no floor inside the range reaches the band target, the one
            nearest to it is returned with `floor_at_bound` True.

    Returns dict: K, floor, c, iterations, converged, history (list of per-round dicts), the
    achieved rate / band_share / participation_ratio_resp, and K_at_bound / floor_at_bound.
    """
    arr = on_counts.to_numpy(dtype=float) if isinstance(on_counts, pd.DataFrame) else np.asarray(on_counts, dtype=float)
    if arr.ndim == 2:
        trials = arr.T[:, :, None]                          # groups (odors) x units x trials
    elif arr.ndim == 4:
        n_trials = arr.shape[1]
        trials = np.concatenate([np.moveaxis(arr[r], 0, -1).transpose(1, 0, 2) for r in range(arr.shape[0])])
    else:
        raise ValueError(f'on_counts must be (units, odors) or (recordings, trials, units, odors), not {arr.shape}')
    # identical count tuples read out identically: group them per odor (and recording)
    groups = []
    for g in range(trials.shape[0]):
        u, cnt = np.unique(trials[g], axis=0, return_counts=True)
        groups.append((u, cnt))
    n_units = np.array([trials.shape[1]] * trials.shape[0], dtype=float)
    width = max(len(c) for _, c in groups)
    counts = np.zeros((len(groups), width, trials.shape[2]))
    weights = np.zeros((len(groups), width))
    for g, (u, cnt) in enumerate(groups):
        counts[g, :len(cnt)] = u
        weights[g, :len(cnt)] = cnt
    sd_mean = float(noise_sd) / np.sqrt(n_trials)

    def signal(K, b):
        x = counts
        with np.errstate(invalid='ignore', divide='ignore'):
            r = np.where(x > 0, (x + b) / (x + b + K), 0.0)
        return r.mean(axis=2)

    def solve_c(K, b):
        return _solve_c((signal(K, b), weights, n_units), target_rate, sd_mean, thresh, c_bounds)

    def c_ok(c):
        return c_bounds[0] * 1.0001 < c < c_bounds[1] * 0.9999

    def band_err(bb):
        # nan where the rate cannot be reached with c inside its bounds (e.g. b near -1, where a
        # model tuned below the target rate reaches it only by inflating c without limit)
        c = solve_c(K, bb)
        return stats_at(K, bb, c)['band_share'] - target_band if c_ok(c) else np.nan

    def stats_at(K, b, c):
        return low_end_stats(c * signal(K, b), sd_mean, band=band, thresh=thresh, weights=weights, n_units=n_units)

    def root(f, lo, hi):
        flo, fhi = f(lo), f(hi)
        if not np.isfinite(flo) or not np.isfinite(fhi) or flo * fhi > 0:
            return (lo if abs(flo) < abs(fhi) else hi), True
        from scipy.optimize import brentq
        return float(brentq(f, lo, hi, xtol=1e-10, rtol=1e-12)), False

    def root_nearest(f, grid, x0):
        # the band share need not be monotone in b once c is re-solved at every b (a model tuned
        # below the target rate reaches it only through its one-spike entries), so the root nearest
        # to the starting floor is taken among the sign changes on a grid
        vals = np.array([f(x) for x in grid])
        ok = np.isfinite(vals)
        brackets = [i for i in range(len(grid) - 1) if ok[i] and ok[i + 1] and vals[i] * vals[i + 1] <= 0]
        if not brackets:
            if not ok.any():
                # the target rate is out of reach at every floor (c at its bounds): keep the start
                return float(x0), True
            i = int(np.nanargmin(np.abs(np.where(ok, vals, np.nan))))
            return float(grid[i]), True
        from scipy.optimize import brentq
        roots = [float(brentq(f, grid[i], grid[i + 1], xtol=1e-10, rtol=1e-12)) if vals[i] != 0 else float(grid[i])
            for i in brackets]
        return min(roots, key=lambda r: abs(r - x0)), False

    b_grid = np.unique(np.concatenate([np.linspace(floor_bounds[0], 0.0, 21),
        np.geomspace(0.05, floor_bounds[1], 40)]))
    K, b = float(K0), float(floor0)
    history = []
    K_at_bound = floor_at_bound = False
    converged = False
    for it in range(1, max_iter + 1):
        K_prev, b_prev = K, b
        if fit_floor:
            b, floor_at_bound = root_nearest(band_err, b_grid, float(floor0))
        fk = lambda lk: stats_at(np.exp(lk), b, solve_c(np.exp(lk), b))['participation_ratio_resp'] - target_pr
        lk, K_at_bound = root(fk, np.log(K_bounds[0]), np.log(K_bounds[1]))
        K = float(np.exp(lk))
        c = solve_c(K, b)
        st = stats_at(K, b, c)
        history.append(dict(iteration=it, K=K, floor=b, c=c, **st))
        if abs(np.log(K) - np.log(K_prev)) < tol and abs(b - b_prev) < tol and it > 1 or (not fit_floor and it > 1):
            converged = True
            break
    return dict(K=K, floor=b, c=c, iterations=it, converged=converged, history=history, K_at_bound=K_at_bound,
        floor_at_bound=floor_at_bound, c_at_bound=not c_ok(c), noise_sd=float(noise_sd),
        n_trials=n_trials, **st)


def cosine_terms(mat: pd.DataFrame, thresh: float = KC_RESPONSE_Z) -> pd.DataFrame:
    """Each odor pair's cosine similarity across units split into three terms by unit: both
    entries responders (above `thresh`), one a responder and the other sub-threshold, both
    sub-threshold; each term is that part of the dot product divided by the product of the two
    odors' lengths, so the three add up to the cosine. Columns 'both', 'one', 'neither', 'cosine';
    index canonical pairs (`pyolfsysm.kc_stats.cosine_terms`)."""
    return _kc_stats.cosine_terms(abbrev_columns(mat), thresh, APL16_ODORS)


_KC_NOISE_SD: Dict[str, float] = {}


def kc_trial_noise_sd(genotype: str = 'TNT3c') -> float:
    """The KC data's measurement noise: the trial-to-trial sd of the non-responding entries
    (trial mean <= 1 z), the root of the mean over those entries of their variance across the
    three trials (ddof 1), pooled over a genotype's recordings (Yang's per-trial values,
    `kcq.noise_floor.load_trials`; `hong2p.imaging.nonresponder_trial_sd`). Cached."""
    if genotype not in _KC_NOISE_SD:
        from kcq.noise_floor import load_trials, recording_arrays
        recs = (np.stack([by_odor[o] for o in APL16_ODORS if o in by_odor], axis=-1)   # trials x KCs x odors
            for by_odor in recording_arrays(load_trials(genotype)).values())
        _KC_NOISE_SD[genotype] = _imaging.nonresponder_trial_sd(recs, KC_RESPONSE_Z)
    return _KC_NOISE_SD[genotype]


def kc_noise_floor(n_sd: float = 2.0) -> float:
    """The noise-defined floor of the 'floor' reading of the KC data: `n_sd` times the sd of a
    trial mean's noise, kc_trial_noise_sd('TNT3c') / sqrt(3) (the control flies' measurement
    noise, used for both genotypes and for the models, as every low-end parameter is set on the
    control flies)."""
    return float(n_sd * kc_trial_noise_sd('TNT3c') / np.sqrt(3))


# ---------------------------------------------------------------------------
# the recordings' per-KC normalization on the model's side (E73): Yang's values appear to be each
# KC's fluorescence divided by its sd over the whole recording (a z-score over the session), so a
# KC's response is relative to its own session; these read a model's spike counts the same way
ZCELL_MAPS = _readout.ZCELL_MAPS
ZCELL_READINGS = _readout.ZCELL_READINGS
OPK_EDGES = _kc_stats.OPK_EDGES


_zcell_amplitude = _readout._zcell_amplitude


# a model's KC responses read as each KC z-scored over its session (E73; pyolfsysm.readout)
session_zscore = _readout.session_zscore


_tail_moments = _readout.tail_moments


def _zcell_odors(z: Dict[str, object], n_odors: int) -> List[str]:
    """The odor labels of a `session_zscore` reading: its columns, abbreviated, or the APL16 odors for an
    unlabeled array (which must then have 16 columns)."""
    odors = z.get('odors')
    if odors is None:
        if n_odors != len(APL16_ODORS):
            raise ValueError('a reading of an unlabeled array needs 16 columns (the APL16 odors)')
        return list(APL16_ODORS)
    return [odor_abbrev(o) for o in odors]


def zcell_pair_distances(z: Dict[str, object], reading: str = 'graded', thresh: float = KC_RESPONSE_Z) -> pd.Series:
    """Cosine odor-odor distances of a `session_zscore` reading, in expectation over its measurement noise
    (closed form, the ratio of the expectations as `expected_pair_distances_noisy`: independent noise adds
    to each odor's squared length and nothing to the products of different odors), mean over its
    recordings: 'graded' (all values, negative ones included: the data's cosine), 'thresholded' (trial
    means <= `thresh` set to 0) or 'binary' (trial mean > `thresh`). Canonical pair Series
    (`pyolfsysm.readout.zcell_pair_distances`)."""
    return _readout.zcell_pair_distances(z, reading, thresh, odors=_zcell_odors(z, z['center'].shape[2]),
        order=APL16_ODORS)


poisson_binomial_bins = _kc_stats.poisson_binomial_bins


def zcell_stats(z: Dict[str, object], thresh: float = KC_RESPONSE_Z) -> Dict[str, object]:
    """Statistics of a `session_zscore` reading under the data's criterion (trial mean > `thresh`), in
    expectation over the measurement noise, pooled over its recordings (`pyolfsysm.readout.zcell_stats`, whose
    docstring lists the keys; `rate_per_odor` is indexed by the abbreviated odors, or the APL16 odors for an
    unlabeled array): the checks against the KC data's (notes/yang_response_normalization/)."""
    return _readout.zcell_stats(z, thresh, odors=_zcell_odors(z, z['center'].shape[2]))


_zcell_rate = _readout._zcell_rate


# the parameters of session_zscore set on the APL-on counts alone (E73; pyolfsysm.readout, criterion 1 z)
calibrate_session_zscore = _readout.calibrate_session_zscore


# ---------------------------------------------------------------------------
# the per-pair APL effect split by the number of odors each KC responds to (E74): a pair's cosine
# similarity is a sum over units, so each condition's similarity splits into the contributions of the
# units in each odors-per-KC class (classes set per condition), and the APL effect (APL-on distance minus
# APL-off distance = APL-off similarity minus APL-on similarity) into per-class effects
OPK_CLASSES = _kc_stats.OPK_CLASSES
CLASS_TERMS = _kc_stats.CLASS_TERMS


opk_class_labels = _kc_stats.opk_class_labels


def _abbrev_matrix(mat, odors: Optional[Sequence[str]] = None):
    """(a DataFrame with its odor columns abbreviated, None) or (an array, its abbreviated odors: `odors`,
    default `APL16_ODORS`), for pyolfsysm.kc_stats."""
    if isinstance(mat, pd.DataFrame):
        return mat.set_axis([odor_abbrev(c) for c in mat.columns], axis=1), None
    return mat, [odor_abbrev(o) for o in (APL16_ODORS if odors is None else odors)]


def _ordered_matrix(mat, odors: Optional[Sequence[str]] = None) -> Tuple[np.ndarray, List[str], List[int]]:
    """(units x odors float array with its columns in `odor_order`, that order, the column positions taken)
    of a DataFrame (odor labels abbreviated) or an array whose columns are `odors` (default `APL16_ODORS`;
    `pyolfsysm.kc_stats.ordered_matrix`)."""
    m, o = _abbrev_matrix(mat, odors)
    return _kc_stats.ordered_matrix(m, o, APL16_ODORS)


class_cosine_terms_array = _kc_stats.class_cosine_terms_array


opk_class_shares = _kc_stats.opk_class_shares


_class_columns = _kc_stats.class_columns


def _abbrev_resp(resp):
    return abbrev_columns(resp) if isinstance(resp, pd.DataFrame) else resp


def _responders(X: np.ndarray, resp, idx: List[int], order: List[str], thresh: float) -> np.ndarray:
    """The responders of an `_ordered_matrix` X: X > thresh, or `resp` in X's column order (a DataFrame by its
    odor labels, an array by the positions `idx` of the matrix it belongs to)."""
    return _kc_stats.responders(X, _abbrev_resp(resp), idx, order, thresh)


def class_cosine_terms(mat, resp=None, *, thresh: float = KC_RESPONSE_Z, edges=OPK_EDGES,
    odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Each odor pair's cosine similarity across units split by the units' odors-per-KC class and, within a
    class, into `cosine_terms`' three terms: both entries responders ('both'), one a responder and the other
    sub-threshold ('one'), both sub-threshold ('neither'), with their sum ('total') (`pyolfsysm.kc_stats.
    class_cosine_terms`). Each part is that part of the dot product over the product of the two odors'
    lengths over all units, so the classes' totals add up to the cosine (class 'all', computed over all units).

    Args:
        mat: (units x odors) values: a DataFrame (odor labels, abbreviated) or an array whose columns are
            `odors` (default `APL16_ODORS`).
        resp: (units x odors) responders, in `mat`'s column order (a bool array or DataFrame); default
            `mat > thresh` (e.g. a model read out as n / (n + K) whose responders are its two-spike entries
            passes resp=counts >= 2).
        thresh: the response criterion when `resp` is None.
        edges: the classes, (lo, hi) inclusive numbers of odors responded to (hi None: no upper limit).

    Returns DataFrame indexed by canonical pairs, columns (class, term) (`opk_class_labels(edges)` + 'all'
    x `CLASS_TERMS`); NaN for pairs with an odor of zero length.
    """
    m, o = _abbrev_matrix(mat, odors)
    return _kc_stats.class_cosine_terms(m, _abbrev_resp(resp), thresh=thresh, edges=edges, odors=o,
        order=APL16_ODORS)


def apl_effect_by_class(on, off, *, resp_on=None, resp_off=None, thresh_on: float = KC_RESPONSE_Z,
    thresh_off: float = KC_RESPONSE_Z, edges=OPK_EDGES, odors: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """The per-pair APL effect split by odors-per-KC class (E74; `pyolfsysm.effects.effect_by_class`): each
    condition's `class_cosine_terms` (classes set by that condition's responses), averaged over its
    recordings, and the effect, APL-off similarity minus APL-on similarity (= APL-on distance minus APL-off
    distance, as `apl_effect`), per class and term. Class 'all' term 'total' of the effect is `apl_effect` of
    the two conditions' cosine distances (mean over recordings).

    Args:
        on, off: each condition's (units x odors) values, or a list of them (recordings; averaged with
            equal weight, NaN pairs of a recording skipped).
        resp_on, resp_off: each condition's responders (one per matrix, as `class_cosine_terms`' `resp`), or
            None for values above `thresh_on` / `thresh_off`.
        edges, odors: as `class_cosine_terms`.

    Returns DataFrame indexed by canonical pairs, columns (block, class, term), block 'on', 'off' (the
    similarity contributions) and 'effect' (off minus on); attrs['class_shares'] {'on', 'off'}: the share of
    units in each class (mean over recordings).
    """
    ab_odors = [odor_abbrev(o) for o in (APL16_ODORS if odors is None else odors)]

    def prep(mats, resps):
        if isinstance(mats, (list, tuple)):
            return ([_abbrev_matrix(m)[0] for m in mats],
                None if resps is None else [_abbrev_resp(r) for r in resps])
        return _abbrev_matrix(mats)[0], _abbrev_resp(resps)

    on, resp_on = prep(on, resp_on)
    off, resp_off = prep(off, resp_off)
    return _effects.effect_by_class(on, off, resp_intact=resp_on, resp_manipulated=resp_off,
        thresh_intact=thresh_on, thresh_manipulated=thresh_off, edges=edges, odors=ab_odors,
        order=APL16_ODORS)


def _unit_cluster(unit: str) -> Optional[Tuple[str, ...]]:
    if unit == 'recording':
        return None
    if unit == 'fly':
        return ('date', 'fly_num')
    raise ValueError(unit)


def resample_recording_means(rec_rdms: Dict[str, pd.DataFrame], rng: np.random.Generator,
    *, unit: str = 'recording') -> Dict[str, pd.Series]:
    """One bootstrap replicate of {genotype: mean pair Series} (`hong2p.stats.resample_group_means`).

    Args:
        rec_rdms: {genotype: per-recording pair DataFrame} (rows (date, fly_num,
            hemibrain), as `yang_recording_rdms`)
        unit: 'recording' draws recordings with replacement within each genotype;
            'fly' draws flies (date, fly_num) with replacement, keeping all hemispheres
            of each drawn fly (the two hemispheres of a fly are not independent)

    Recordings are then averaged with equal weight, as in the report's recipe.
    """
    return _stats.resample_group_means(rec_rdms, rng, cluster=_unit_cluster(unit))


def bootstrap_recordings(fn: Callable[[Dict[str, pd.Series]], Dict[str, float]],
    rec_rdms: Dict[str, pd.DataFrame], *, n_boot: int = 1000, seed: int = 0,
    ci: Tuple[float, float] = (5, 95), unit: str = 'recording') -> pd.DataFrame:
    """Bootstraps statistics over recordings (or flies; `hong2p.stats.bootstrap_groups`).

    Args:
        fn: maps {genotype: mean pair Series} -> {stat: value}
        rec_rdms: {genotype: per-recording pair DataFrame} (`yang_recording_rdms`)
        n_boot: resamples (drawn with replacement within each genotype)
        unit: 'recording' or 'fly' (see `resample_recording_means`)

    Returns DataFrame indexed by stat: observed, boot_mean, boot_sd, lo, hi.
    """
    return _stats.bootstrap_groups(fn, rec_rdms, n_boot=n_boot, seed=seed, ci=ci, cluster=_unit_cluster(unit))


# reliability of a mean pair-distance matrix, from random split halves of recordings (hong2p.stats)
split_half_reliability = _stats.split_half_reliability


# ---------------------------------------------------------------------------
# per-recording KC x odor matrices and RDM pipeline variants (E7; atlas)
RECORDING_KEYS: List[str] = ['date', 'fly_num', 'hemibrain']


def yang_recording_matrices(genotype: str, with_pfo: bool = False
    ) -> Dict[Tuple[str, int, str], pd.DataFrame]:
    """{(date, fly_num, hemibrain): KCs x odors trial-mean matrix} for one genotype.

    Odor columns are in `APL16_ODORS` order (plus `SOLVENT` last if `with_pfo`). Values
    are Yang's z-scored fluorescence trial means (negatives included).
    """
    tm = yang_trialmeans()
    tm = tm[tm.genotype == genotype]
    odors = list(APL16_ODORS) + ([SOLVENT] if with_pfo else [])
    mats = {}
    for rec, sub in tm.groupby(RECORDING_KEYS, sort=True):
        mat = sub.pivot(index='roi', columns='odor', values='value')[odors]
        assert not mat.isna().any().any(), rec
        mats[(str(rec[0]), int(rec[1]), str(rec[2]))] = mat
    return mats


def transform_matrix(mat: pd.DataFrame, *, values: str = 'graded', cells: str = 'all',
    thresh: float = 1.0) -> pd.DataFrame:
    """One RDM pipeline variant applied to a KCs x odors matrix (`hong2p.rdm.transform_matrix`).

    Args:
        values: 'graded' (as is), 'binary' (response > thresh -> 1.0 else 0.0), or
            'thresholded' (graded values kept where > thresh, else 0)
        cells: 'all' KCs, or 'responsive' (KCs with a response > thresh to at least one
            odor column of `mat`)
        thresh: response threshold (Yang's data are z-scores; she has used 0.75 and 1)
    """
    return _rdm.transform_matrix(mat, values=values, cells=cells, thresh=thresh)


def _stack_recordings(rows: Dict[tuple, pd.Series]) -> pd.DataFrame:
    df = pd.DataFrame(rows).T
    df.index = pd.MultiIndex.from_tuples(df.index, names=RECORDING_KEYS)
    return df


def recording_rdms_from_matrices(mats: Dict[tuple, pd.DataFrame], *,
    measure: str = 'cosine', drop_silent: bool = False, **transform_kws) -> pd.DataFrame:
    """Per-recording pair distances (rows recordings, columns canonical pairs) for one
    pipeline variant (`transform_kws` -> `transform_matrix`; `measure`, `drop_silent` ->
    `pair_distances`). Defaults reproduce `yang_recording_rdms(genotype)`.
    """
    return _stack_recordings({rec: pair_distances(transform_matrix(m, **transform_kws),
        measure=measure, drop_silent=drop_silent) for rec, m in mats.items()
    })


def recording_response_fractions(mats: Dict[tuple, pd.DataFrame], thresh: float = 1.0
    ) -> pd.DataFrame:
    """Rows recordings, columns odors: fraction of KCs with response > thresh."""
    return _stack_recordings({rec: (m > thresh).mean(axis='index') for rec, m in mats.items()})


def match_threshold(mats: Dict[tuple, pd.DataFrame], target_frac: float, *,
    lo: float = -1.0, hi: float = 8.0, tol: float = 1e-4) -> float:
    """Threshold at which the recording-mean response fraction (mean over odors) equals
    `target_frac`, by bisection. Pass a single-recording dict to match one recording.
    """
    def frac(th: float) -> float:
        return float(recording_response_fractions(mats, th).mean(axis='columns').mean())

    assert frac(lo) > target_frac > frac(hi), (frac(lo), target_frac, frac(hi))
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if frac(mid) > target_frac:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def density_null_similarity(mat: pd.DataFrame) -> pd.Series:
    """Expected cosine similarity between odor columns of a KCs x odors matrix if KC
    identities were shuffled independently within each odor (per-odor response
    distributions, hence response counts, kept; `hong2p.rdm.density_null_similarity`).
    Returned as a canonical pair Series of similarities (not distances).
    """
    return _rdm.density_null_similarity(abbrev_columns(mat), APL16_ODORS)


def recording_density_null_rdms(mats: Dict[tuple, pd.DataFrame], **transform_kws
    ) -> pd.DataFrame:
    """Per-recording expected 1 - cosine under the within-odor KC-identity shuffle null."""
    return _stack_recordings({rec: 1 - density_null_similarity(
        transform_matrix(m, **transform_kws)) for rec, m in mats.items()
    })


def shuffled_pair_distances(mat: pd.DataFrame, rng: np.random.Generator, *,
    measure: str = 'cosine', **transform_kws) -> pd.Series:
    """Pair distances after independently permuting KC identities within each odor
    column (one draw from the density null; `density_null_similarity` is its mean;
    `hong2p.rdm.shuffled_pair_distances`, `transform_kws` as `transform_matrix`'s)."""
    return _rdm.shuffled_pair_distances(abbrev_columns(mat), rng, measure=measure, order=APL16_ODORS,
        **dict(dict(thresh=1.0), **transform_kws))


# ---------------------------------------------------------------------------
# the three repeats of the recordings (E81, the flies' variability across the three repeats; atlas):
# single-trial distances by kind of trial pair, the same-repeat gap, noise-free distances from products
# of different repeats, and the parts of a response's variability across repeats. A recording (or a
# model's pseudo-recording) is an array (n_repeats, units, odors) of single-trial values; in the model,
# whose trials are exchangeable, the trial index stands in for the repeat.
REPEAT_KINDS: Tuple[str, ...] = _imaging.REPEAT_PAIRS
TRIAL_VALUES: Tuple[str, ...] = _imaging.TRIAL_VALUES
YANG_DATASETS: Dict[str, int] = {'TNT': 19, 'TNT3c': 20}
_YANG_TRIALS: Dict[tuple, Dict[tuple, Tuple[np.ndarray, List[str]]]] = {}


def yang_recording_trials(genotype: str, units: str = 'z', odors: Optional[Sequence[str]] = None
    ) -> Dict[Tuple[str, int, str], np.ndarray]:
    """One genotype's single-trial KC responses, per recording.

    Args:
        genotype: 'TNT3c' (control flies) or 'TNT' (APL output blocked (TNT)).
        units: 'z' (Yang's whole-recording z, olfdb's values: data/yang_apl16_kc_olfdb/<genotype>/
            dataset_<id>_long.csv) or 'noise' (baseline-noise units: data/yang_apl16_kc_noise_units/
            per_trial_long.parquet, the realism bands' units).
        odors: odor abbreviations, in the order of the last axis (default `APL16_ODORS`; 'pfo' is the
            solvent).

    Returns {(date, fly_num, hemibrain): array (n_repeats, KCs, odors)}, repeats in the order of olfdb's
    `repeat` (taken as the repeat's number in the session), KCs sorted by roi. Cached per (genotype,
    units); the arrays returned are copies.
    """
    odors = list(APL16_ODORS if odors is None else odors)
    key = (genotype, units)
    if key not in _YANG_TRIALS:
        cols = RECORDING_KEYS + ['roi', 'odor', 'repeat', 'value']
        if units == 'z':
            d = pd.read_csv(yang_data_dir() / genotype / f'dataset_{YANG_DATASETS[genotype]}_long.csv',
                usecols=cols)
        elif units == 'noise':
            d = pd.read_parquet(shared_dir() / 'data' / 'yang_apl16_kc_noise_units' / 'per_trial_long.parquet',
                columns=['genotype'] + cols)
            d = d[d['genotype'] == genotype].drop(columns='genotype')
        else:
            raise ValueError(f"units must be 'z' or 'noise', not {units!r}")
        d['abbrev'] = d['odor'].map(odor_abbrev)
        out = {}
        for rec, sub in d.groupby(RECORDING_KEYS, sort=True):
            w = sub.set_index(['repeat', 'roi', 'abbrev'])['value'].unstack('abbrev').sort_index()
            reps = sorted(w.index.get_level_values('repeat').unique())
            rois = sorted(w.index.get_level_values('roi').unique())
            w = w.reindex(pd.MultiIndex.from_product([reps, rois], names=['repeat', 'roi']))
            if w.isna().any().any():
                raise ValueError(f'{genotype} {rec}: missing trials')
            arr = w.to_numpy(dtype=float).reshape(len(reps), len(rois), w.shape[1])
            out[(str(rec[0]), int(rec[1]), str(rec[2]))] = (arr, list(w.columns))
        _YANG_TRIALS[key] = out
    res = {}
    for rec, (arr, cols) in _YANG_TRIALS[key].items():
        missing = [o for o in odors if o not in cols]
        if missing:
            raise KeyError(f'{genotype} {rec}: no odors {missing}')
        res[rec] = arr[:, :, [cols.index(o) for o in odors]].copy()
    return res


def trial_values(trials: np.ndarray, values: str = 'graded', thresh: float = KC_RESPONSE_Z) -> np.ndarray:
    """Single trials read as `values`: 'graded' (unchanged), 'thresholded' (values at or below `thresh`
    set to 0, trial by trial, as a non-responding model KC contributes 0) or 'binary' (1 above `thresh`,
    else 0). Works on any array shape (`hong2p.imaging.trial_values`)."""
    return _imaging.trial_values(trials, values, thresh)


# cosine similarities between every pair of single-trial response vectors of one recording
trial_cosines = _imaging.trial_cosines


_repeat_masks = _imaging.repeat_pair_masks


_kind_means = _imaging.mean_over_trial_pairs


def distances_by_kind(trials: np.ndarray, odors: Sequence[str], *, values: str = 'graded',
    thresh: float = KC_RESPONSE_Z) -> pd.DataFrame:
    """One recording's single-trial cosine distances (1 - cosine), averaged by kind of trial pair.

    Args:
        trials: (n_repeats, units, odors), odors in the order of `odors`.
        values, thresh: how the single trials are read (`trial_values`).

    Returns a DataFrame indexed (odor_a, odor_b): the canonical two-odor pairs (in `odors` order), then
    one same-odor row per odor (odor_a == odor_b); columns `REPEAT_KINDS`: 'same_repeat' (the two odors
    in the same repeat; NaN for a same-odor row), 'different_repeat' (any two different repeats),
    'adjacent' (repeats one apart) and 'first_third' (repeats two apart; NaN with fewer than three
    repeats). Each entry is the mean over the trial pairs of that kind (NaN cosines skipped).
    """
    odors = list(odors)
    X = trial_values(trials, values, thresh)
    C = trial_cosines(X)
    masks = _repeat_masks(X.shape[0])
    D = {k: 1 - _kind_means(C, m) for k, m in masks.items()}
    pl = pair_list(odors)
    idx = {o: i for i, o in enumerate(odors)}
    rows = [(a, b) for a, b in pl] + [(o, o) for o in odors]
    out = pd.DataFrame({k: [D[k][idx[a], idx[b]] for a, b in rows] for k in REPEAT_KINDS},
        index=pd.MultiIndex.from_tuples(rows, names=['odor_a', 'odor_b']))
    out.loc[[(o, o) for o in odors], 'same_repeat'] = np.nan
    return out


def same_repeat_gap(by_kind: pd.DataFrame, pairs: Optional[Sequence[Tuple[str, str]]] = None) -> float:
    """The same-repeat gap of one recording: the mean over two-odor pairs (default all of `by_kind`'s;
    `pairs` restricts them, e.g. to the six pairs with near-identical ORN input; either order) of the
    single-trial distance between different repeats minus the distance within the same repeat
    (`distances_by_kind`). Positive when two odors of one repeat are more alike than two odors of
    different repeats, i.e. when part of the trial-to-trial variability is shared by the odors of a
    repeat."""
    a = by_kind.index.get_level_values(0)
    b = by_kind.index.get_level_values(1)
    keep = np.asarray(a != b)
    if pairs is not None:
        want = {tuple(sorted(p)) for p in pairs}
        keep &= np.array([tuple(sorted(x)) in want for x in by_kind.index])
        if not keep.any():
            raise KeyError(f'none of {sorted(want)} in by_kind')
    g = by_kind.loc[keep, 'different_repeat'] - by_kind.loc[keep, 'same_repeat']
    return float(g.mean())


_noise_free_matrix = _imaging.noise_free_matrix


def noise_free_pair_distances(trials: np.ndarray, odors: Sequence[str], *, repeats: str = 'different',
    values: str = 'graded', thresh: float = KC_RESPONSE_Z) -> pd.Series:
    """Noise-free cosine distances of one recording from products of different repeats
    (`hong2p.imaging.noise_free_pair_distances`).

    For each pair (a, b): 1 - C_ab / sqrt(C_aa C_bb), C_ab the mean single-trial cosine between odor a in
    repeat i and odor b in repeat j over the repeat pairs i != j selected by `repeats`, C_aa the same for
    odor a with itself. Noise independent between repeats adds nothing, in expectation, to a product of
    two different repeats, whether or not it is shared by the odors of a repeat; a drift common to an
    odor's repeats acts on C_ab and C_aa alike and cancels in the ratio. Not clipped.

    Args:
        trials: (n_repeats, units, odors).
        repeats: 'different' (every i != j), 'adjacent' (|i - j| == 1) or 'first_third' (|i - j| == 2).
        values, thresh: how the single trials are read (`trial_values`).

    Returns a canonical pair Series (index (odor_a, odor_b)).
    """
    return _imaging.noise_free_pair_distances(trials, odors, repeats=repeats, values=values, thresh=thresh,
        order=APL16_ODORS)


def repeat_pair_table(trials: np.ndarray, odors: Sequence[str], *, values: str = 'graded',
    thresh: float = KC_RESPONSE_Z) -> pd.DataFrame:
    """Every pair distance of one recording that E81 compares, one row per canonical two-odor pair.

    Columns: 'trial_mean' (1 - cosine of the two odors' trial means, the main target's distance; the
    trial mean of the trials read as `values`), 'corrected' (the noise correction in use,
    `noise_corrected_pair_distances`: noise assumed independent between odors), 'noise_free',
    'noise_free_adjacent', 'noise_free_first_third' (`noise_free_pair_distances` with each `repeats`),
    and the single-trial distances by kind (`REPEAT_KINDS`, from `distances_by_kind`).
    """
    odors = list(odors)
    X = trial_values(trials, values, thresh)
    C = trial_cosines(X)
    masks = _repeat_masks(X.shape[0])
    pl = pair_list(odors)
    idx = {o: i for i, o in enumerate(odors)}
    ia = [idx[a] for a, _ in pl]
    ib = [idx[b] for _, b in pl]
    xbar = X.mean(axis=0)
    nrm = np.linalg.norm(xbar, axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        S = (xbar.T @ xbar) / np.outer(nrm, nrm)
    cols = {'trial_mean': 1 - S[ia, ib],
        'corrected': noise_corrected_pair_distances(X, odors).to_numpy() if X.shape[0] > 1 else np.full(len(pl), np.nan)}
    for name, kind in (('noise_free', 'different_repeat'), ('noise_free_adjacent', 'adjacent'),
            ('noise_free_first_third', 'first_third')):
        cols[name] = _noise_free_matrix(C, masks[kind])[ia, ib]
    for kind in REPEAT_KINDS:
        cols[kind] = (1 - _kind_means(C, masks[kind]))[ia, ib]
    return pd.DataFrame(cols, index=pd.MultiIndex.from_tuples(pl, names=['odor_a', 'odor_b']))


def presentation_pairs(trials: np.ndarray, odors: Sequence[str], order: pd.DataFrame, *, values: str = 'graded',
    thresh: float = KC_RESPONSE_Z, volume_rate: Optional[float] = None) -> pd.DataFrame:
    """One recording's single-trial distance between every two presentations of two different odors, with the
    time between them (E81 part A6, the stimulus order).

    Args:
        trials: (n_repeats, units, odors), odors in the order of `odors` (`yang_recording_trials`).
        odors: the odors of the last axis. Presentations of other stimuli (the solvent) form no pairs but keep
            their places in the order.
        order: the recording's presentations, one row each: 'presentation' (index in time order), 'block',
            'odor' (abbreviation), 'repeat' (the index of the presentation's trial on the first axis of
            `trials`) and, with `volume_rate`, 'onset_volume'.
        values, thresh: how the single trials are read (`trial_values`).
        volume_rate: volumes per second, for 'separation_s'.

    Returns a DataFrame, one row per pair of presentations of two different odors of `odors`: 'odor_a',
    'odor_b' (the canonical pair, as `pair_list` orders it), 'presentation_a', 'presentation_b', 'block_a',
    'block_b' (of odor_a's and odor_b's presentation), 'position_a', 'position_b' (in the block, when `order` has
    'position'), 'same_block', 'separation' (presentations apart, solvent presentations counted), 'separation_s'
    (seconds between the onsets; NaN without `volume_rate`) and 'distance' (1 - cosine of the two single trials;
    NaN where either is all zero).
    """
    odors = list(odors)
    idx = {o: i for i, o in enumerate(odors)}
    o = order[order['odor'].isin(odors)].sort_values('presentation')
    if o.duplicated(['odor', 'repeat']).any():
        raise ValueError('order: an (odor, repeat) given to two presentations')
    if (o['repeat'] < 0).any() or (o['repeat'] >= trials.shape[0]).any():
        raise ValueError('order: a repeat outside the trials')
    X = trial_values(trials, values, thresh)
    C = trial_cosines(X)
    rank = {p: i for i, p in enumerate(odor_order(odors))}
    pres = o['presentation'].to_numpy()
    blk = o['block'].to_numpy()
    rep = o['repeat'].to_numpy()
    oi = np.array([idx[x] for x in o['odor']])
    onset = o['onset_volume'].to_numpy(dtype=float) if volume_rate is not None else np.full(len(o), np.nan)
    i, j = np.triu_indices(len(o), k=1)
    keep = oi[i] != oi[j]
    i, j = i[keep], j[keep]
    swap = np.array([rank[odors[oi[a]]] > rank[odors[oi[b]]] for a, b in zip(i, j)], dtype=bool)
    ia, ib = np.where(swap, j, i), np.where(swap, i, j)
    with np.errstate(invalid='ignore'):
        dist = 1 - C[rep[ia], oi[ia], rep[ib], oi[ib]]
    out = {'odor_a': [odors[k] for k in oi[ia]], 'odor_b': [odors[k] for k in oi[ib]],
        'presentation_a': pres[ia], 'presentation_b': pres[ib], 'block_a': blk[ia], 'block_b': blk[ib]}
    if 'position' in o.columns:
        pos = o['position'].to_numpy()
        out.update(position_a=pos[ia], position_b=pos[ib])
    out.update(same_block=blk[ia] == blk[ib], separation=np.abs(pres[ib] - pres[ia]),
        separation_s=np.abs(onset[ib] - onset[ia]) / (volume_rate if volume_rate else np.nan), distance=dist)
    return pd.DataFrame(out)


def repeat_response_stats(trials: np.ndarray, odors: Sequence[str], *, thresh: float = KC_RESPONSE_Z
    ) -> pd.DataFrame:
    """Per repeat and odor of one recording: 'rate' (fraction of units whose single-trial value is above
    `thresh`), 'amplitude' (the mean single-trial value of those units; NaN if none), 'n_resp', 'n_units'.
    Indexed (repeat (0-based), odor)."""
    X = np.asarray(trials, dtype=float)
    R, U, O = X.shape
    resp = X > thresh
    n_resp = resp.sum(axis=1)                                   # (R, O)
    with np.errstate(invalid='ignore', divide='ignore'):
        amp = np.where(resp, X, 0.0).sum(axis=1) / n_resp
    odors = list(odors)
    idx = pd.MultiIndex.from_product([range(R), odors], names=['repeat', 'odor'])
    return pd.DataFrame({'rate': (n_resp / U).ravel(), 'amplitude': amp.ravel(), 'n_resp': n_resp.ravel(),
        'n_units': np.full(R * O, U)}, index=idx)


def _pair_cov_sums(h: np.ndarray, A: np.ndarray, B: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """Sum over units, repeats and ordered odor pairs o != p with (unit, o) in A and (unit, p) in B
    (B = A if None) of h[r, k, o] h[r, k, p], and the number of such (unit, o, p) triples.
    A, B: (units, odors) bool; when B is given, A must be a subset of B."""
    hA = np.where(A[None], h, 0.0)
    if B is None:
        s = hA.sum(axis=2)
        q = (hA ** 2).sum(axis=2)
        nA = A.sum(axis=1)
        return float((s ** 2 - q).sum()), float((nA * (nA - 1)).sum())
    hB = np.where(B[None], h, 0.0)
    total = float((hA.sum(axis=2) * hB.sum(axis=2) - (hA ** 2).sum(axis=2)).sum())
    return total, float((A.sum(axis=1) * (B.sum(axis=1) - 1)).sum())


def _repeat_gains(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The gain of each presentation and the deviations of one recording (`repeat_variance_parts`).

    The dot products over units between two different repeats of one odor, P_rs = x_r . x_s (r != s), are
    g_r g_s |m|^2 in expectation whatever noise is independent between repeats, so each repeat's gain follows
    from them: a_r = sqrt(g_r^2 |m|^2) from log P_rs = log a_r + log a_s by least squares (with three repeats
    exactly a_r^2 = P_rs P_rt / P_st), and the relative gain is a_r over the mean of the a's.

    Returns (d, gdev, m, Q): d (R, U, O) the deviations from each entry's mean over repeats; gdev (R, O) each
    presentation's relative gain minus 1 (0 for every presentation of an odor with a product not positive, and
    with fewer than three repeats, where gain and noise cannot be told apart); m (U, O) the trial means; Q
    (U, O) the mean over i != j of x_i x_j per entry, an estimate of the entry's squared noise-free mean.
    """
    R, U, O = X.shape
    m = X.mean(axis=0)
    d = X - m[None]
    gdev = np.zeros((R, O))
    if R >= 3:
        P = np.einsum('iko,jko->ijo', X, X)                      # (R, R, O)
        iu, ju = np.triu_indices(R, 1)
        A = np.zeros((len(iu), R))
        A[np.arange(len(iu)), iu] = 1.0
        A[np.arange(len(iu)), ju] = 1.0
        Pu = P[iu, ju]                                           # (n_pairs, O)
        ok = (Pu > 0).all(axis=0)
        if ok.any():
            la = np.linalg.lstsq(A, np.log(Pu[:, ok]), rcond=None)[0]     # (R, n_ok)
            a = np.exp(la)
            gdev[:, ok] = a / a.mean(axis=0) - 1.0
    S = X.sum(axis=0)
    Q = (S ** 2 - (X ** 2).sum(axis=0)) / (R * (R - 1))
    return d, gdev, m, Q


def repeat_variance_parts(trials: np.ndarray, odors: Sequence[str], *, thresh: float = KC_RESPONSE_Z,
    n_amp_bins: int = 4) -> pd.DataFrame:
    """The parts of one recording's variability across repeats.

    Each entry's deviation from its own mean over repeats, d[r, k, o] = x[r, k, o] - mean_r x[., k, o],
    is split into three parts, with x the units' (KCs') single-trial values:
    - 'gain': common to all units of a presentation (repeat r of odor o), proportional to the units'
      noise-free response to that odor (a cosine does not see it; it also takes up a run-down common to the
      whole repeat). The presentation's relative gain comes from products of different repeats
      (`_repeat_gains`), so the noise in the trial mean does not attenuate it as a regression on the trial
      mean would; the part's variance in an entry is the variance over repeats of its odor's relative gains
      times the entry's squared noise-free mean (estimated from products of different repeats);
    - 'shared': common to one unit's odors within a repeat: the covariance between two different odors'
      residual deviations (the gain part, relative gain deviation times trial mean, removed) of the same unit
      in the same repeat;
    - 'rest': the remainder, var_total - var_gain - var_shared (independent between the odors of a repeat).
    Variances are per trial (sums over repeats of squared deviations, or of products, divided by
    n_repeats - 1, which is unbiased for a part independent between repeats), averaged over the entries
    (or, for the shared part, the unit x odor-pair products) of a class.

    Args:
        trials: (n_repeats, units, odors).
        thresh: an entry (unit, odor) is a responder when its trial mean is above `thresh`, silent
            otherwise; the shared part of 'all', 'responders' and 'silent' is taken over pairs of odors
            whose two entries are both in the class.
        n_amp_bins: the responders' entries are also split into this many quantile bins of their trial
            mean ('responders_amp1' the lowest); a bin's shared part is the covariance of its entries
            with the same unit's other responder entries, to show whether each part scales with the
            amplitude.

    Returns a DataFrame indexed by class ('all', 'responders', 'silent', 'responders_amp1', ...): n_entries,
    n_pairs (unit x ordered odor pairs behind the shared part), mean_amplitude (mean trial mean),
    var_total, var_gain, var_shared, var_rest, share_gain, share_shared, share_rest (of var_total),
    sd_trial (sqrt(var_total)) and corr_shared (var_shared over var_total - var_gain: the correlation of
    two odors' deviations within a repeat once the gain part is removed).
    """
    X = np.asarray(trials, dtype=float)
    R = X.shape[0]
    if R < 2:
        raise ValueError('need at least two repeats')
    d, gdev, m, Q = _repeat_gains(X)
    h = d - gdev[:, None, :] * m[None]
    v_tot = (d ** 2).sum(axis=0) / (R - 1)                      # per entry
    v_gain = ((gdev ** 2).sum(axis=0) / (R - 1))[None, :] * Q
    resp = m > thresh
    classes = {'all': (np.ones_like(resp), None), 'responders': (resp, None), 'silent': (~resp, None)}
    if n_amp_bins and resp.sum() >= n_amp_bins:
        edges = np.quantile(m[resp], np.linspace(0, 1, n_amp_bins + 1))
        b = np.clip(np.searchsorted(edges, m, side='right') - 1, 0, n_amp_bins - 1)
        for i in range(n_amp_bins):
            classes[f'responders_amp{i + 1}'] = (resp & (b == i), resp)
    rows = {}
    for name, (A, B) in classes.items():
        n = int(A.sum())
        if n == 0:
            rows[name] = dict(n_entries=0, n_pairs=0)
            continue
        s, npairs = _pair_cov_sums(h, A, B)
        vt, vg = float(v_tot[A].mean()), float(v_gain[A].mean())
        vs = s / (R - 1) / npairs if npairs > 0 else np.nan
        rows[name] = dict(n_entries=n, n_pairs=int(npairs), mean_amplitude=float(m[A].mean()), var_total=vt,
            var_gain=vg, var_shared=vs, var_rest=vt - vg - vs, share_gain=vg / vt, share_shared=vs / vt,
            share_rest=(vt - vg - vs) / vt, sd_trial=float(np.sqrt(vt)), corr_shared=vs / (vt - vg))
    return pd.DataFrame(rows).T.rename_axis('class')


def coresponder_deviation_corr(trials: np.ndarray, odors: Sequence[str], *, thresh: float = KC_RESPONSE_Z,
    min_units: int = 5, remove_gain: bool = False) -> pd.Series:
    """Per two-odor pair of one recording: the correlation, over the units responding to both odors (trial
    mean above `thresh`) and the repeats, between the two odors' deviations from each unit's mean over
    repeats (the shared state among co-responders; with `remove_gain`, after the gain part of
    `repeat_variance_parts`). NaN for pairs with fewer than `min_units` such units.
    Returns a canonical pair Series."""
    X = np.asarray(trials, dtype=float)
    d, gdev, m, _ = _repeat_gains(X)
    if remove_gain:
        d = d - gdev[:, None, :] * m[None]
    resp = m > thresh
    odors = list(odors)
    pl = pair_list(odors)
    idx = {o: i for i, o in enumerate(odors)}
    vals = []
    for a, b in pl:
        ia, ib = idx[a], idx[b]
        k = resp[:, ia] & resp[:, ib]
        if k.sum() < min_units:
            vals.append(np.nan)
            continue
        da, db = d[:, k, ia].ravel(), d[:, k, ib].ravel()
        den = np.sqrt((da ** 2).sum() * (db ** 2).sum())
        vals.append(float((da * db).sum() / den) if den > 0 else np.nan)
    return pd.Series(vals, index=pd.MultiIndex.from_tuples(pl, names=['odor_a', 'odor_b']))


# ---------------------------------------------------------------------------
# plots (matplotlib; sequential single-hue maps for distances, diverging for residuals)
def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def plot_square(ax, pairs: pd.Series, title: str = '', vmin: float = 0.0,
    vmax: float = 1.0, cmap: str = 'Blues', cbar_label: str = '1 - similarity') -> None:
    """Heatmap of a pair Series as an odor x odor matrix."""
    sq = square_from_pairs(pairs)
    im = ax.imshow(sq.to_numpy(), vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(len(sq)))
    ax.set_xticklabels(sq.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(sq)))
    ax.set_yticklabels(sq.index, fontsize=7)
    ax.set_title(title, fontsize=9)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(cbar_label, fontsize=7)
    cb.ax.tick_params(labelsize=7)


def plot_on_vs_off(ax, on: pd.Series, off: pd.Series, *, label: str = 'model',
    color: str = CONDITION_COLORS['model'], data_on: Optional[pd.Series] = None,
    data_off: Optional[pd.Series] = None, highlight_block: str = FRUITY_SMALL_BLOCK
    ) -> None:
    """Scatter of APL-on vs APL-off pair distances with the isotonic fit; block pairs marked.

    Optionally overlays the data (control vs TNT) in a second color.
    """
    on, off = align_pairs(on, off)
    iso = IsotonicRegression(increasing=True, out_of_bounds='clip').fit(off, on)
    xs = np.linspace(off.min(), off.max(), 100)
    blk = np.array([block_label(a, b) == highlight_block for a, b in on.index])
    ax.scatter(off[~blk], on[~blk], s=12, color=color, alpha=0.7, label=label,
        edgecolors='none'
    )
    ax.scatter(off[blk], on[blk], s=26, facecolors='none', edgecolors=color, linewidths=1.2,
        label=f'{label}: {highlight_block} pairs'
    )
    ax.plot(xs, iso.predict(xs), color=color, lw=1.5, alpha=0.9)
    if data_on is not None and data_off is not None:
        d_on, d_off = align_pairs(data_on, data_off)
        dc = CONDITION_COLORS['data']
        ax.scatter(d_off, d_on, s=12, color=dc, alpha=0.7, label='data (TNT3c vs TNT)',
            edgecolors='none'
        )
        iso_d = IsotonicRegression(increasing=True, out_of_bounds='clip').fit(d_off, d_on)
        xd = np.linspace(d_off.min(), d_off.max(), 100)
        ax.plot(xd, iso_d.predict(xd), color=dc, lw=1.5, alpha=0.9)
    lim = [0, 1.05]
    ax.plot(lim, lim, color='#bbbbbb', lw=0.8, ls='--')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('1 - similarity, APL off (TNT)')
    ax.set_ylabel('1 - similarity, APL on (control)')
    ax.legend(fontsize=7, frameon=False)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def plot_per_odor_vs_drive(ax, per_odor: pd.DataFrame, *, col: str = 'iso_resid',
    color: str = CONDITION_COLORS['model'], label: str = 'model') -> None:
    """Per-odor APL effect (mean isotonic residual of its pairs) vs ORN drive, labeled."""
    d = per_odor.dropna(subset=['drive', col])
    ax.scatter(d['drive'], d[col], s=18, color=color, label=label, edgecolors='none')
    for o, row in d.iterrows():
        ax.annotate(o, (row['drive'], row[col]), fontsize=6, xytext=(2, 2),
            textcoords='offset points', color='#444444'
        )
    ax.axhline(0, color='#bbbbbb', lw=0.8, ls='--')
    ax.set_xlabel('ORN drive (mean over glomeruli)')
    ax.set_ylabel(f'per-odor APL effect ({col})')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def plot_model_vs_data(ax, model: pd.Series, data: pd.Series, *, label: str = '',
    color: str = CONDITION_COLORS['model']) -> None:
    """Scatter of model vs data pair distances, with unity line."""
    m, d = align_pairs(model, data)
    ax.scatter(d, m, s=12, color=color, alpha=0.75, edgecolors='none', label=label)
    lim = [0, 1.05]
    ax.plot(lim, lim, color='#bbbbbb', lw=0.8, ls='--')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('data: 1 - similarity')
    ax.set_ylabel('model: 1 - similarity')
    if label:
        ax.legend(fontsize=7, frameon=False)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------------------
def selftest(n_boot: int = 500) -> None:
    """Reproduces the report's data-vs-data numbers and prints recording-bootstrap CIs."""
    pd.set_option('display.width', 160)
    rep = report_matrices()
    recs = {g: yang_recording_rdms(g, 'cosine') for g in ('TNT3c', 'TNT')}
    for g in recs:
        diff = (recs[g].mean() - rep[g]).abs().max()
        print(f'{g}: {len(recs[g])} recordings; recording-mean cosine vs report matrix: '
            f'max |diff| = {diff:.2e}'
        )
    orn = orn16_consensus_pairs()
    drive = consensus_orn_drive()
    st, tabs = compare_on_off(rep['TNT3c'], rep['TNT'], orn=orn, drive=drive)
    print('\ncontrol (TNT3c) vs TNT, 120 pairs (report: spearman 0.869, lin R2 0.848, '
        'iso R2 0.94, block +0.055, spearman ORN 0.655 -> 0.748, per-odor vs drive 0.50):'
    )
    for k, v in st.items():
        print(f'  {k:>34s}: {v:.4f}' if isinstance(v, float) else f'  {k:>34s}: {v}')
    print('\nblocks (mean isotonic residual):')
    print(tabs['blocks'].sort_values('iso_resid', ascending=False).round(3).to_string())
    print('\nreference comparisons (66 megamat pairs):')
    ref = reference_pairs()
    rows = {g: compare_to_references(rep[g], ref) for g in ('TNT3c', 'TNT')}
    rows['WT (Remy)'] = compare_to_references(ref['kc_soma'], ref)
    print(pd.DataFrame(rows).round(3).to_string())

    def fn(means: Dict[str, pd.Series]) -> Dict[str, float]:
        s, _ = compare_on_off(means['TNT3c'], means['TNT'], orn=orn, drive=drive)
        keys = ['mean_on', 'mean_off', 'spearman_on_off', 'lin_r2', 'iso_r2',
            'block_fruity_small_iso', 'block_other_between_iso', 'spearman_orn_on',
            'spearman_orn_off', 'gap_closed', 'apl_share_resid_var',
            'spearman_per_odor_iso_vs_drive',
        ]
        out = {k: s[k] for k in keys}
        for g in ('TNT3c', 'TNT'):
            r = compare_to_references(means[g], ref)
            out[f'{g}_spearman_vcf_hyp'] = r['spearman_vcf_hyp']
            out[f'{g}_spearman_resid_vs_vcf_hyp'] = r['spearman_resid_vs_vcf_hyp']
            out[f'{g}_spearman_resid_vs_vcf_raw'] = r['spearman_resid_vs_vcf_raw']
            out[f'{g}_spearman_kc_claw'] = r['spearman_kc_claw']
        return out

    print(f'\nbootstrap over recordings ({n_boot} draws; 90% CI):')
    print(bootstrap_recordings(fn, recs, n_boot=n_boot).round(3).to_string())
    print('\nsplit-half reliability of the mean matrices (Spearman; ceiling for a model):')
    for g in ('TNT3c', 'TNT'):
        for bz in (None, 1.0):
            r = split_half_reliability(yang_recording_rdms(g, 'cosine', binarize=bz))
            print(f"  {g} {'binarized' if bz else 'raw':9s}: split-half {r['split_half']:.3f} "
                f"[{r['split_half_lo']:.3f}, {r['split_half_hi']:.3f}], full-sample "
                f"{r['spearman_brown_full']:.3f}, model ceiling {r['ceiling']:.3f}")
    print('\nbinarized data (trial-mean > 1), recording-mean 1 - cosine: '
        + ', '.join(f'{g} {yang_mean_rdm(g, binarize=1.0).mean():.3f}' for g in ('TNT3c', 'TNT'))
        + ' (raw: ' + ', '.join(f'{g} {rep[g].mean():.3f}' for g in ('TNT3c', 'TNT')) + ')'
    )
    print('\n' + APPROX_CONC_NOTE)


if __name__ == '__main__':
    selftest()
